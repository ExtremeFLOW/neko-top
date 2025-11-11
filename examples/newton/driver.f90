program usrneko

!  use stdlib_io_npy, only : save_npy
  use LightKrylov
  use LightKrylov, only: wp => dp, rtol => rtol_dp
  use LightKrylov_Logger
  use LightKrylov_Constants
  use global_coef, only: global_coef_t, global_coef_getter
  use neko_vector, only: state_vector_t
  use neko_system, only: non_linear_propagator_t
  use neko_jacobian, only: jacobian_t
  use neko_linop, only: linear_propagator_t
  use simulation_m, only: simulation_t
  use LightKrylov_IterativeSolvers, only: write_results_cdp
  use json_module, only: json_file
  use json_utils, only: json_get
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
  type(simulation_t), target :: simulation
  !> nonlinear propagator.
  type(non_linear_propagator_t) :: non_linear
  !> a jacobian.
  type(jacobian_t) :: jacobian
  !> A way to access coef globally
  type(global_coef_t), target :: my_global_coef_getter
  !> State vectors
   type(state_vector_t), allocatable :: bf, dx, residual

  !---------------------------------------------------
  !-----     KRYLOV-BASED EIGENDECOMPOSITION     -----
  !---------------------------------------------------
  !> Exponential propagator.
  type(linear_propagator_t), allocatable :: A
  !> Sampling time.
  real(kind=wp) :: tau
  !> Number of eigenvalues we wish to converge.
  integer, parameter :: nev = 15
  !> Size of Krylov subspace.
  integer, parameter :: kdim = 126
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
  !> Miscellaneous.
  integer :: i, j
  ! JSON related arguments
  integer :: argc
  character(len=256) :: parameter_file
  type(json_file) :: parameters, design_parameters
  ! MPI parameters
  integer :: ierr

  !=============================================================================

  !----------------------------------
  !-----     INITIALIZATION     -----
  !----------------------------------

  !> Set up logging
  call logger_setup()

  ! Initialize the MPI environment

  call MPI_Init(ierr)
  call neko_top_register_types()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))

  !> Initialize propagators.
  ! -------------------------------------------------------------------------- !
  ! do a linear simulation
  simulation%adjoint_case%if_adjoint = .false.
  call simulation%init(parameters)
  call non_linear%init(simulation)
  non_linear%jacobian = jacobian
    select type (f => non_linear%jacobian)
    type is (jacobian_t)
       call f%init(simulation)
    end select

!> Extract the global coef from neko
  my_global_coef_getter%global_coef = non_linear%simulation%neko_case%fluid%c_Xh
  global_coef_getter => my_global_coef_getter

  !> initial guess is baseflow loaded
  allocate(bf); call bf%init()
  call field_copy(bf%u, non_linear%simulation%neko_case%fluid%u)
  call field_copy(bf%v, non_linear%simulation%neko_case%fluid%v)
  call field_copy(bf%w, non_linear%simulation%neko_case%fluid%w)
  call field_copy(bf%p, non_linear%simulation%neko_case%fluid%p)

  ! call newton(non_linear, bf, gmres_rdp, info)
  ! call newton(non_linear, bf, gmres_rdp, info, scheduler=dynamic_tol_dp)

  ! Now we have a steady baseflow, let's compute the spectra
  allocate(A)
  call A%init(simulation)

  !> Initialize Krylov subspace.
  allocate(X(nev)); call initialize_krylov_subspace(X)

  !> initialize writer
  allocate(X_writer); call X_writer%init()

  !> Call to LightKrylov.
  call eigs(A, X, lambda, residuals, info)

  call check_info(info, 'eigs', module=this_module, procedure='main')
  !> Transform eigenspectrum from unit-disk to standard complex plane.
  lambda = log(lambda) / tau

  !--------------------------------
  !-----     SAVE TO DISK     -----
  !--------------------------------

  !> Save eigenvectors to disk.
  do i = 1, nev
      call X_writer%axpby(1.0_wp, X(i), 0.0_wp)
      call X_writer%write(i)
  enddo

  !> Clean up
  call A%free()
  do i = 1, nev
    call X(i)%free()
  enddo

end program usrneko
