program usrneko
  use simulation, only: simulation_t
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

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
  !> The design type
  class(design_t), allocatable :: design
  !> The problem type
  type(problem_t) :: problem
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

  call simulation%init(parameters)

  call design_factory(design, parameters, simulation)

  call problem%init(parameters, design, simulation)
  call optimizer_factory(optimizer, parameters, problem, design, simulation)

  ! -------------------------------------------------------------------------- !
  ! Execute the optimization

  call optimizer%run(problem, design, simulation)

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call optimizer%free()
  call problem%free()
  call design%free()
  call simulation%free()

  if (allocated(design)) deallocate(design)
  if (allocated(optimizer)) deallocate(optimizer)

end program usrneko
