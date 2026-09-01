program usrneko

  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  use mpi_f08, only: MPI_Init, MPI_Wtime, MPI_COMM_WORLD


  use example_problem_mma, only: mma_obj
  use objective, only: objective_t

  use simplefield_design, only: simplefield_design_t
  use neko, only: neko_init, neko_finalize, neko_solve
  use nekotop_logger, only: nekotop_log
  use case, only: case_t
  use field, only: field_t
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t

  use comm, only: pe_rank
  use device, only: device_memcpy, HOST_TO_DEVICE
  use neko_config, only: NEKO_BCKND_DEVICE

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
  type(simplefield_design_t) :: des

  !> The problem type
  type(problem_t) :: prob
  class(objective_t), allocatable :: obj

  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt

  integer :: nloc
  type(vector_t) :: xcoord, ycoord, zcoord, initdesign
  real(kind=rp) :: t_start, t_end

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
  call nekotop_log%init(env_prefix = "NEKO_TOP")
  call neko_field%init(neko_case%msh, neko_case%fluid%Xh, "neko_field")
  nloc = neko_field%dof%size()
  call xcoord%init(nloc)
  call ycoord%init(nloc)
  call zcoord%init(nloc)
  xcoord%x = reshape(neko_field%dof%x, [nloc])
  ycoord%x = reshape(neko_field%dof%y, [nloc])
  zcoord%x = reshape(neko_field%dof%z, [nloc])

  if (NEKO_BCKND_DEVICE .eq. 1) then
     call device_memcpy(xcoord%x, xcoord%x_d, nloc, &
          HOST_TO_DEVICE, sync = .true.)
     call device_memcpy(ycoord%x, ycoord%x_d, nloc, &
          HOST_TO_DEVICE, sync = .true.)
     call device_memcpy(zcoord%x, zcoord%x_d, nloc, &
          HOST_TO_DEVICE, sync = .true.)
  end if
  call des%init_from_components(nloc, xcoord, ycoord, zcoord, neko_field)

  if (pe_rank == 0) then
     print *, "nloc=", nloc, "number of design variables=", des%size_global()
     print *, "max(xcoord%x)=", maxval(xcoord%x), "min(xcoord%x)=", &
          minval(xcoord%x), "max(ycoord%x)=", maxval(ycoord%x), &
          "min(ycoord%x)=", minval(ycoord%x), "max(zcoord%x)=", &
          maxval(zcoord%x), "min(zcoord%x)=", minval(zcoord%x)
  end if

  ! initialize the design
  call initdesign%init(des%size())
  initdesign = 1.0_rp
  call des%update_design(initdesign)

  ! -------------------------------------------------------------------------- !
  ! Construct the problem
  !
  ! This subroutine calculates function values and gradients
  ! for the unconstrained problem:
  !   minimize \f$\sum_(j = 1,..,n) (x_j - X_{j,GLL})^2/nglobal \f$

  ! allocate as subtype mma_obj
  allocate(mma_obj :: obj)
  call obj%init_json(parameters, des)

  ! update obj and sensitivities for the init design
  call obj%update_value(des)
  call obj%update_sensitivity(des)

  if (pe_rank == 0) then
     print *, "objective value for the initial design=", obj%value
  end if

  ! initialize the problem
  call prob%init(parameters, des)
  call prob%add_objective(obj)

  ! -------------------------------------------------------------------------- !
  ! Execute the optimization
  call optimizer_factory(opt, parameters, prob, des)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_start = MPI_Wtime()

  call opt%run(prob, des)

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_end = MPI_Wtime()

  if (pe_rank == 0) then
     print *, "opt%run execution time:", t_end - t_start, "seconds"
  end if

  ! -------------------------------------------------------------------------- !
  ! Clean up the components

  call opt%free()
  call prob%free()
  call des%free()

  if (allocated(opt)) deallocate(opt)

  call nekotop_log%free()
  call neko_finalize(neko_case)
end program usrneko
