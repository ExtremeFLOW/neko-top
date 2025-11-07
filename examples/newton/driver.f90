program usrneko

!  use stdlib_io_npy, only : save_npy
  use LightKrylov
  use LightKrylov, only: wp => dp, rtol => rtol_dp
  use LightKrylov_Logger
  use LightKrylov_Constants
  use newton, only: linear_propagator, non_linear_propagator, state_vector
  use global_coef, only: global_coef_t, global_coef_getter
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
  !-----     LINEAR OPERATOR INVESTIGATED     -----
  !------------------------------------------------

  !> nonlinear propagator.
  type(non_linear_propagator) :: non_linear
  !> nonlinear propagator.
  type(linear_propagator) :: linear
  !> A way to access coef globally
  type(global_coef_t), target :: my_global_coef_getter
  !> State vectors
   type(state_vector), allocatable :: bf, dx, residual

  integer :: info
  !> writer
  type(state_vector), allocatable :: X_writer
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
  call non_linear%init(parameters)
  non_linear%jacobian = linear
    select type (f => non_linear%jacobian)
    type is (linear_propagator)
       call f%init(non_linear%simulation)
    end select

!> Extract the global coef from neko
  my_global_coef_getter%global_coef = non_linear%simulation%neko_case%fluid%c_Xh
  global_coef_getter => my_global_coef_getter

  !> initialize writer
  allocate(X_writer); call X_writer%init()
  allocate(bf); call bf%init()
  call field_copy(bf%u, non_linear%simulation%neko_case%fluid%u)
  call field_copy(bf%v, non_linear%simulation%neko_case%fluid%v)
  call field_copy(bf%w, non_linear%simulation%neko_case%fluid%w)
  call field_copy(bf%p, non_linear%simulation%neko_case%fluid%p)

  ! call newton(non_linear, bf, gmres_rdp, info)
  call newton(non_linear, bf, gmres_rdp, info, scheduler=dynamic_tol_dp)


end program usrneko
