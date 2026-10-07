program usrneko
  use LightKrylov
  use LightKrylov, only: wp => dp
  use LightKrylov_AbstractVectors, only: init_like_basis
  use LightKrylov_Constants
  use LightKrylov_Logger
  use field_math, only: field_copy
  use field_output, only: field_output_t
  use fld_file, only: fld_file_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use json_utils_ext, only: json_read_file
  use neko, only: neko_finalize, neko_init
  use case, only: case_t
  use fixed_point_system, only: fixed_point_jacobian_t, fixed_point_system_t
  use linop_fixed_baseflow, only: linop_fixed_baseflow_t
  use neko_top, only: neko_top_register_types
  use neko_vector, only: state_vector_t
  use num_types, only: neko_rp => rp, neko_sp => sp
  use perturb_case, only: perturb_case_t
  use stability_case_utils, only: stability_case_free, stability_case_init
  use utils, only: neko_error
  implicit none

  character(len=*), parameter :: this_module = 'Cylinder Re40 fixed point'
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case
  type(fixed_point_system_t) :: non_linear
  type(state_vector_t), allocatable :: bf
  type(linop_fixed_baseflow_t), allocatable :: A
  type(state_vector_t), allocatable :: X(:), X_writer
  type(state_vector_t) :: vector_mold
  complex(kind=wp), allocatable :: lambda(:)
  real(kind=wp), allocatable :: residuals(:)
  real(kind=wp) :: newton_tol, eigs_tol, tau
  integer :: argc, nev, info, i
  logical :: compute_eigs, write_outer_iterations
  character(len=256) :: parameter_file
  character(len=:), allocatable :: outer_output_prefix
  type(json_file) :: parameters

  call neko_init()
  call neko_top_register_types()
  call logger_setup()

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  parameters = json_read_file(trim(parameter_file))
  call json_get_or_default(parameters, &
       'case.lightkrylov.newton.relative_tolerance', newton_tol, 1.0e-6_wp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.eigensolver.enabled', compute_eigs, .false.)
  call json_get_or_default(parameters, &
       'case.lightkrylov.eigensolver.number_of_eigenvalues', nev, 2)
  call json_get_or_default(parameters, &
       'case.lightkrylov.eigensolver.tolerance', eigs_tol, rtol_dp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.newton.write_outer_iterations', &
       write_outer_iterations, .false.)
  call json_get_or_default(parameters, &
       'case.lightkrylov.newton.outer_iteration_output_prefix', &
       outer_output_prefix, 'outer_baseflow_newton_Re40')

  if (newton_tol <= 0.0_wp) call neko_error('Newton tolerance must be positive')
  if (nev < 1) call neko_error('Number of eigenvalues must be positive')
  if (eigs_tol <= 0.0_wp) call neko_error('Eigensolver tolerance must be positive')

  call stability_case_init(base_case, parameters)
  lns_case%if_adjoint = .false.
  call lns_case%init(base_case)
  tau = real(base_case%time%end_time, kind=wp)

  call non_linear%init(base_case)
  if (write_outer_iterations) then
     call non_linear%enable_outer_iteration_output(trim(outer_output_prefix))
  end if
  vector_mold%coef => base_case%fluid%c_Xh
  vector_mold%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)

  allocate(fixed_point_jacobian_t :: non_linear%jacobian)
  select type (f => non_linear%jacobian)
  type is (fixed_point_jacobian_t)
     call f%init(base_case, lns_case)
  end select

  allocate(bf)
  call bf%init_like(vector_mold)
  call field_copy(bf%u, base_case%fluid%u)
  call field_copy(bf%v, base_case%fluid%v)
  call field_copy(bf%w, base_case%fluid%w)

  call newton(non_linear, bf, gmres_rdp, info, &
       scheduler=dynamic_tol_dp, rtol=newton_tol)
  if (info < 0) call check_info(info, 'newton', module=this_module, &
       procedure='main')

  call field_copy(base_case%fluid%u, bf%u)
  call field_copy(base_case%fluid%v, bf%v)
  call field_copy(base_case%fluid%w, bf%w)
  call write_baseflow_fields('baseflow_newton_Re40', bf, 0)

  if (compute_eigs) then
     allocate(A)
     call A%init(base_case, lns_case)
     call A%set_baseflow(bf)

     allocate(X(nev))
     call init_like_basis(X, bf)
     call initialize_krylov_subspace(X)

     allocate(X_writer)
     call X_writer%init_like(bf)

     call eigs(A, X, lambda, residuals, info, tolerance=eigs_tol)
     call check_info(info, 'eigs', module=this_module, procedure='main')
     lambda = log(lambda) / tau

     do i = 1, nev
        call X_writer%axpby(1.0_wp, X(i), 0.0_wp)
        call X_writer%write(i - 1)
        call write_state_fields('mode_newton_Re40', X_writer, i - 1)
     end do
  end if

  if (allocated(A)) then
     call A%free()
     deallocate(A)
  end if
  call non_linear%free()
  if (allocated(X_writer)) then
     call X_writer%free()
     deallocate(X_writer)
  end if
  if (allocated(bf)) then
     call bf%free()
     deallocate(bf)
  end if
  if (allocated(X)) then
     do i = 1, size(X)
        call X(i)%free()
     end do
     deallocate(X)
  end if
  if (allocated(lambda)) deallocate(lambda)
  if (allocated(residuals)) deallocate(residuals)
  call lns_case%free()
  call stability_case_free(base_case)
  call neko_finalize()
contains

  subroutine write_state_fields(name, vec, idx)
    character(len=*), intent(in) :: name
    type(state_vector_t), intent(inout) :: vec
    integer, intent(in) :: idx
    type(field_output_t) :: output

    call execute_command_line('mkdir -p viz')
    call output%init(trim(name), 3, precision=neko_sp, path='viz/', format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .true.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, vec%u)
    call output%fields%assign_to_field(2, vec%v)
    call output%fields%assign_to_field(3, vec%w)
    call output%set_counter(idx - 1)
    call output%sample(real(idx, kind=neko_rp))
    call output%free()
  end subroutine write_state_fields

  subroutine write_baseflow_fields(name, vec, idx)
    character(len=*), intent(in) :: name
    type(state_vector_t), intent(inout) :: vec
    integer, intent(in) :: idx
    type(field_output_t) :: output

    call execute_command_line('mkdir -p viz')
    call output%init(trim(name), 4, precision=neko_sp, path='viz/', &
         format='nek5000')
    select type (ft => output%file_%file_type)
    type is (fld_file_t)
       ft%skip_pressure = .false.
       ft%skip_velocity = .false.
       ft%skip_temperature = .true.
    end select
    call output%fields%assign_to_field(1, base_case%fluid%p)
    call output%fields%assign_to_field(2, vec%u)
    call output%fields%assign_to_field(3, vec%v)
    call output%fields%assign_to_field(4, vec%w)
    call output%set_counter(idx - 1)
    call output%sample(real(idx, kind=neko_rp))
    call output%free()
  end subroutine write_baseflow_fields

end program usrneko
