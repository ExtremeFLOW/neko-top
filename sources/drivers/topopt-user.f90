program topopt_user
  use simulation_m, only: simulation_t
  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory
  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file
  use user, only: user_setup
  use design, only: design_t, design_factory
  use neko_top, only: neko_top_register_types

  use mpi_f08, only: MPI_Init

  implicit none

  ! JSON related arguments
  integer :: argc
  type(json_file) :: parameters
  character(len=256) :: parameter_file

  ! MPI parameters
  integer :: ierr

  !> The simulation we are working with
  type(simulation_t) :: sim
  !> The problem type
  type(problem_t) :: prob
  !> The design type
  class(design_t), allocatable :: des
  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt

  ! -------------------------------------------------------------------------- !
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

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  ! initialize the user additions for the forward (through the neko interface)
  call user_setup(sim%neko_case%usr)

  ! initialize the simulation
  call sim%init(parameters)

  ! initialize the design
  call design_factory(des, parameters, sim)

  ! initialize the problem
  call prob%init(parameters, des, sim)

  ! initialize the optimizer
  call optimizer_factory(opt, parameters, prob, des, sim)

  ! -------------------------------------------------------------------------- !
  ! Execute the optimization

  call opt%run(prob, des, sim)

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call opt%free()
  if (allocated(opt)) deallocate(opt)

  call prob%free()
  call des%free()
  call sim%free()

end program topopt_user
