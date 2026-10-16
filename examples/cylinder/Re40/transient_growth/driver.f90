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
  use linop_fixed_baseflow, only: linop_fixed_baseflow_t
  use neko_top, only: neko_top_register_types
  use neko_vector, only: state_vector_t
  use num_types, only: neko_rp => rp, neko_sp => sp
  use perturb_case, only: perturb_case_t
  use stability_case_utils, only: stability_case_free, stability_case_init
  use utils, only: neko_error
  implicit none

  character(len=*), parameter :: this_module = 'Cylinder transient growth'
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case
  type(perturb_case_t), target :: adjoint_case
  type(linop_fixed_baseflow_t) :: A
  type(state_vector_t), allocatable :: U(:), V(:), baseflow
  type(state_vector_t) :: vector_mold
  real(kind=wp), allocatable :: S(:), residuals(:)
  real(kind=wp) :: svds_tol
  integer :: argc, nsv, kdim, info, i
  logical :: write_intermediate
  character(len=256) :: parameter_file
  type(json_file) :: parameters

  call neko_init()
  call neko_top_register_types()
  call logger_setup()

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  parameters = json_read_file(trim(parameter_file))
  call json_get_or_default(parameters, &
       'case.lightkrylov.transient_growth.number_of_singular_values', nsv, 4)
  call json_get_or_default(parameters, &
       'case.lightkrylov.transient_growth.krylov_dimension', kdim, 4 * nsv)
  call json_get_or_default(parameters, &
       'case.lightkrylov.transient_growth.tolerance', svds_tol, rtol_dp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.transient_growth.write_intermediate', &
       write_intermediate, .true.)

  if (nsv < 1) call neko_error('Number of singular values must be positive')
  if (kdim < nsv) call neko_error('Krylov dimension must be at least nsv')
  if (svds_tol <= 0.0_wp) call neko_error('SVD tolerance must be positive')

  call stability_case_init(base_case, parameters)

  lns_case%if_adjoint = .false.
  call lns_case%init(base_case)
  adjoint_case%if_adjoint = .true.
  call adjoint_case%init(base_case)

  vector_mold%coef => base_case%fluid%c_Xh
  vector_mold%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)

  allocate(baseflow)
  call baseflow%init_like(vector_mold)
  call field_copy(baseflow%u, base_case%fluid%u)
  call field_copy(baseflow%v, base_case%fluid%v)
  call field_copy(baseflow%w, base_case%fluid%w)
  call write_state_fields('baseflow_transient_growth', baseflow, 0)

  call A%init(base_case, lns_case, adjoint_case)
  call A%set_baseflow(baseflow)

  allocate(U(nsv), V(nsv))
  call init_like_basis(U, vector_mold)
  call init_like_basis(V, vector_mold)
  call initialize_krylov_subspace(U)
  call initialize_krylov_subspace(V)

  call svds(A, U, S, V, residuals, info, kdim=kdim, tolerance=svds_tol, &
       write_intermediate=write_intermediate)

  call write_singular_spectrum(S, residuals)

  do i = 1, nsv
     call write_state_fields('prt', V(i), i - 1)
     call write_state_fields('rsp', U(i), i - 1)
  end do

  call A%free()
  do i = 1, nsv
     call U(i)%free()
     call V(i)%free()
  end do
  deallocate(U, V)
  call baseflow%free()
  deallocate(baseflow)
  if (allocated(S)) deallocate(S)
  if (allocated(residuals)) deallocate(residuals)
  call lns_case%free()
  call adjoint_case%free()
  call stability_case_free(base_case)
  call neko_finalize()

contains

  subroutine write_singular_spectrum(sigma, res)
    real(kind=wp), intent(in) :: sigma(:), res(:)
    integer :: j

    open(unit=1234, file='singular_spectrum.dat', status='replace', &
         action='write')
    write(1234,'(A16,1x,A16,1x,A16)') '#', 'SIGMA', 'RESIDUAL'
    do j = 1, size(sigma)
       write(1234,'(I16,1x,ES16.8,1x,ES16.8)') j, sigma(j), res(j)
    end do
    close(1234)
  end subroutine write_singular_spectrum

  subroutine write_state_fields(name, vec, idx)
    character(len=*), intent(in) :: name
    type(state_vector_t), intent(inout) :: vec
    integer, intent(in) :: idx
    type(field_output_t) :: output

    call execute_command_line('mkdir -p viz')
    call output%init(trim(name), 3, precision=neko_sp, path='viz/', &
         format='nek5000')
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

end program usrneko
