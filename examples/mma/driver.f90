program usrneko
  use neko, only: neko_init, neko_finalize
  use design, only: design_t, design_factory
  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory
  use example_problem, only: mma_obj, mma_con

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  implicit none

  ! JSON related arguments
  integer :: argc
  type(json_file) :: parameters
  character(len=256) :: parameter_file

  !> The design type
  class(design_t), allocatable :: design
  !> The problem type
  type(problem_t) :: problem
  type(mma_obj) :: obj
  type(mma_con) :: con_1, con_2
  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: optimizer

  ! -------------------------------------------------------------------------- !
  ! Initialize the MPI environment

  call neko_init()

  ! -------------------------------------------------------------------------- !
  ! Read the parameters file as the first terminal argument

  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))

  ! -------------------------------------------------------------------------- !
  ! Initialization of the components

  call design_factory(design, parameters)
  call problem%init(parameters, design)

  ! -------------------------------------------------------------------------- !
  ! Construct the problem
  !
  ! This subroutine calculates function values and gradients
  ! for "toy problem 3":
  !
  !   minimize \f$\sum_(j = 1,..,n) x_j/n \f$
  ! subject to \f$\sum_(j = 1,..,n) (x_j - X_{j,GLL})^2 = 0 \f$

  call obj%init_from_components("Objective", design)
  call con_1%init_from_components("Positive", design, 1)
  call con_2%init_from_components("Negative", design, -1)

  call problem%add_objective(obj)
  call problem%add_constraint(con_1)
  call problem%add_constraint(con_2)

  ! -------------------------------------------------------------------------- !
  ! Execute the optimization
  call optimizer_factory(optimizer, parameters, problem, design)

  call optimizer%run(problem, design)

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call optimizer%free()
  call problem%free()
  call design%free()

  if (allocated(design)) deallocate(design)
  if (allocated(optimizer)) deallocate(optimizer)

  call neko_finalize()
end program usrneko
