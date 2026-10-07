program usrneko

!  use stdlib_io_npy, only : save_npy
  use LightKrylov
  use LightKrylov, only: wp => dp
  use LightKrylov_AbstractVectors, only: init_like_basis
  use LightKrylov_Logger
  use LightKrylov_Constants
  use neko, only: neko_init, neko_finalize
  use neko_vector, only: state_vector_t
  use fixed_point_system, only: fixed_point_jacobian_t, fixed_point_system_t
  use linop_fixed_baseflow, only: linop_fixed_baseflow_t
  use case, only: case_t
  use perturb_case, only: perturb_case_t
  use stability_case_utils, only: stability_case_free, stability_case_init
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file
  use neko_top, only: neko_top_register_types
  use field_math, only: field_copy
  implicit none

  character(len=128), parameter :: this_module = 'Example cylinder'

  !------------------------------------------------
  !-----                 NEWTON               -----
  !------------------------------------------------

  !> simulation
  type(case_t), target :: base_case
  type(perturb_case_t), target :: lns_case
  !> nonlinear propagator.
  type(fixed_point_system_t) :: non_linear
  !> State vectors
  type(state_vector_t), allocatable :: bf

  !---------------------------------------------------
  !-----     KRYLOV-BASED EIGENDECOMPOSITION     -----
  !---------------------------------------------------
  !> Exponential propagator.
  type(linop_fixed_baseflow_t), allocatable :: A
  !> Sampling time.
  real(kind=wp) :: tau
  !> Number of eigenvalues we wish to converge.
  integer :: nev
  !> Eigensolver convergence tolerance.
  real(kind=wp) :: eigs_tol
  !> Krylov subspace.
  type(state_vector_t), allocatable :: X(:)
  !> Eigenvalues.
  complex(kind=wp), allocatable :: lambda(:)
  !> Residual.
  real(kind=wp), allocatable    :: residuals(:)
  !> Information flag.
  integer          :: info
  !> writer
  type(state_vector_t), allocatable :: X_writer
  !> Mold used to initialize Neko-backed LightKrylov vectors.
  type(state_vector_t) :: vector_mold
  !> Miscellaneous.
  integer :: i
  logical :: compute_eigs
  ! JSON related arguments
  integer :: argc
  character(len=256) :: parameter_file
  type(json_file) :: parameters

  real(kind=wp) :: newton_tol

  !=============================================================================

  !----------------------------------
  !-----     INITIALIZATION     -----
  !----------------------------------

  ! Initialize the Neko environment before any LightKrylov/MPI setup.
  call neko_init()
  call neko_top_register_types()
  call logger_setup()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))
  call json_get_or_default(parameters, &
       'case.lightkrylov.newton.relative_tolerance', newton_tol, 1.0e-3_wp)
  call json_get_or_default(parameters, &
       'case.lightkrylov.eigensolver.enabled', compute_eigs, .false.)
  call json_get_or_default(parameters, &
       'case.lightkrylov.eigensolver.number_of_eigenvalues', nev, 2)
  call json_get_or_default(parameters, &
       'case.lightkrylov.eigensolver.tolerance', eigs_tol, rtol_dp)

  if (newton_tol <= 0.0_wp) then
     call neko_error('LightKrylov Newton relative_tolerance must be positive')
  end if
  if (nev < 1) then
     call neko_error('LightKrylov number_of_eigenvalues must be positive')
  end if
  if (eigs_tol <= 0.0_wp) then
     call neko_error('LightKrylov eigensolver tolerance must be positive')
  end if

  !> Initialize propagators.
  ! -------------------------------------------------------------------------- !
  ! do a linear simulation
  call stability_case_init(base_case, parameters)
  lns_case%if_adjoint = .false.
  call lns_case%init(base_case)
  tau = real(base_case%time%end_time, kind=wp)
  call non_linear%init(base_case)
  vector_mold%coef => base_case%fluid%c_Xh
  vector_mold%if_2d = (base_case%fluid%c_Xh%msh%gdim .eq. 2)
  allocate(fixed_point_jacobian_t :: non_linear%jacobian)
  select type (f => non_linear%jacobian)
  type is (fixed_point_jacobian_t)
     call f%init(base_case, lns_case)
  end select

  !> initial guess is baseflow loaded
  allocate(bf)
  call bf%init_like(vector_mold)
  call field_copy(bf%u, base_case%fluid%u)
  call field_copy(bf%v, base_case%fluid%v)
  call field_copy(bf%w, base_case%fluid%w)

  call newton(non_linear, bf, gmres_rdp, info, &
       scheduler=dynamic_tol_dp, rtol=newton_tol)
  if (info < 0) call check_info(info, 'newton', module=this_module, &
       procedure='main')

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
     end do
  end if

  !> Clean up
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

end program usrneko
