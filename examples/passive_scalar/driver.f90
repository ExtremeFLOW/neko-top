program usrneko
  use simulation_m, only: simulation_t
  use problem, only: problem_t
  use optimizer, only : optimizer_t, optimizer_factory
  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file
  use user, only: user_setup
  use design, only: design_t, design_factory

  use mpi_f08, only: MPI_Init

  implicit none

  ! JSON related arguments
  integer :: argc
  type(json_file) :: parameters
  character(len=256) :: parameter_file

  ! MPI parameters
  integer :: ierr

  !> The simulation we are working with
  type(simulation_t) :: simulation
  !> The problem type
  type(problem_t) :: problem
  !> The design type
  class(design_t), allocatable :: design
  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: optimizer

  ! -------------------------------------------------------------------------- !
  ! Initialize the MPI environment
  call MPI_Init(ierr)

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
  call user_setup(simulation%neko_case%usr)

  ! initialize the simulation
  call simulation%init(parameters)

  ! initialize the design
  call design_factory(design, parameters, simulation)

  ! initialize the problem
  call problem%init(parameters, design, simulation)

  ! initialize the optimizer
  call optimizer_factory(optimizer, parameters, problem, design, simulation)

  ! -------------------------------------------------------------------------- !
  ! Execute the optimization

  call optimizer%run(problem, design, simulation)

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call optimizer%free()
  if (allocated(optimizer)) deallocate(optimizer)

  call problem%free()
  call design%free()
  call simulation%free()

end program usrneko
