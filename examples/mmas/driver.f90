program usrneko

  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  use mpi_f08, only: MPI_Init


  use example_problem, only: mma_obj, mma_con

  use simplefield_design, only: simplefield_design_t
  use neko, only: neko_init, neko_finalize, neko_solve
  use case, only: case_t
  use field, only: field_t
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t

  implicit none

  ! JSON related arguments
  integer :: argc
  type(json_file) :: parameters
  character(len=256) :: parameter_file

  ! MPI parameters
  integer :: ierr

  !> neko case and field types
  type(case_t) :: neko_case
  type(field_t) :: neko_field

  ! !> The design
  type(simplefield_design_t) :: design

  !> The problem type
  type(problem_t) :: problem
  type(mma_obj), allocatable :: obj
  type(mma_con), allocatable :: con_1, con_2

  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: optimizer

  integer :: nloc
  type(vector_t) :: xcoord, ycoord, zcoord,  finaldesign, initdesign

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
  call neko_init(neko_case)
  call neko_field%init(neko_case%msh, neko_case%fluid%Xh, "neko_field")
  nloc = neko_field%dof%size()
  call xcoord%init(nloc)
  call ycoord%init(nloc)
  call zcoord%init(nloc)
  xcoord%x = reshape(neko_field%dof%x, [nloc])
  ycoord%x = reshape(neko_field%dof%y, [nloc])
  zcoord%x = reshape(neko_field%dof%z, [nloc])
  call design%init_from_components(nloc, xcoord, ycoord, zcoord, neko_field)

  print *, "nloc=", nloc, "number of design variables=", design%size_global()
  print *, "max(xcoord%x)=", maxval(xcoord%x), "min(xcoord%x)=", minval(xcoord%x), &
       "max(ycoord%x)=", maxval(ycoord%x), "min(ycoord%x)=", minval(ycoord%x), &
       "max(zcoord%x)=", maxval(zcoord%x), "min(zcoord%x)=", minval(zcoord%x)
  ! initialize the design
  call initdesign%init(design%size())
  initdesign%x = 1.0_rp
  call design%update_design(initdesign)

   ! -------------------------------------------------------------------------- !
  ! Construct the problem
  !
  ! This subroutine calculates function values and gradients
  ! for "toy problem 3":
  !
  !   minimize \f$\sum_(j = 1,..,n) x_j/n \f$
  ! subject to \f$\sum_(j = 1,..,n) (x_j - X_{j,GLL})^2 = 0 \f$

  allocate(obj)
  allocate(con_1)
  allocate(con_2)

  call obj%init_from_components("Objective", design)
  call con_1%init_from_components("Positive", design, 1)
  call con_2%init_from_components("Negative", design, -1)
  
  ! update obj and cons and sensitivities for the init design
  call obj%update_value(design)
  call obj%update_sensitivity(design)
  call con_1%update_value(design)
  call con_1%update_sensitivity(design)
  call con_2%update_value(design)
  call con_2%update_sensitivity(design)
  
  print *, "objective value for the initial design=", obj%value

  ! initialize the problem
  call problem%init(parameters, design)
  
  call problem%add_objective(obj)
  call problem%add_constraint(con_1)
  call problem%add_constraint(con_2)

  ! -------------------------------------------------------------------------- !
  ! Execute the optimization
  call optimizer_factory(optimizer, parameters, problem, design)


  call optimizer%run(problem, design)

  call finaldesign%init(design%size_global())
  finaldesign = design%get_values()
  print *, "min(design%values)=", minval(finaldesign%x), &
       "max(design%values)=", maxval(finaldesign%x)
  ! -------------------------------------------------------------------------- !
  ! Clean up the components


  call neko_finalize(neko_case)
  call optimizer%free()
  call problem%free()
  call design%free()

  if (allocated(optimizer)) deallocate(optimizer)

end program usrneko
