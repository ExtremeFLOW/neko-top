program usrneko
  use mma_optimizer, only : mma_optimizer_t
  use steady_state_problem, only: steady_state_problem_t
  use topopt_design, only: topopt_design_t
  use num_types, only : rp
  use mask_ops, only: mask_exterior_const
  use json_module, only: json_file
  use mpi_f08, only: MPI_COMM_WORLD, MPI_Init, MPI_Comm_rank
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  integer :: argc
  type(json_file) :: parameters
  character(len=256) :: parameter_file
  real(kind=rp) :: tolerance

  integer :: mpi_rank, ierr

  !> The problem type
  type(steady_state_problem_t) :: problem
  !> The design type
  type(topopt_design_t) :: design
  !> The optimizer (in this case mma)
  type(mma_optimizer_t) :: optimizer

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)

  ! Read the parameters file as the first terminal argument
  argc = command_argument_count()
  if (argc .lt. 1) call neko_error('Missing parameter file')
  call get_command_argument(1, parameter_file)

  ! Read the parameters file
  parameters = json_read_file(trim(parameter_file))

  ! init the problem (base)
  call problem%init()

  ! init the problem, with the design
  call problem%init_design(design)

  ! init the optimizer
  call optimizer%init(problem, design)

  tolerance = 1.0e-3_rp
  max_iter = 100
  call optimizer%run(problem, tolerance, max_iter)

  call problem%free()
  call design%free()
  call optimizer%free()

end program usrneko
