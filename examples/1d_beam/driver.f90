program usrneko

  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  use mpi_f08, only: MPI_Init,  MPI_Wtime, MPI_COMM_WORLD, mpi_in_place, &
       mpi_allreduce, mpi_exscan


  use example_problem, only: mma_obj, beamweight_obj, mma_con

  use design_3dto1d , only: design_3dto1d_t
  use neko, only: neko_init, neko_finalize, neko_solve
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
  type(design_3dto1d_t) :: des

  !> The problem type
  type(problem_t) :: prob
  type(mma_obj), allocatable :: deflection
  type(beamweight_obj), allocatable :: beamweight
  type(mma_con), allocatable :: con_1 !, con_2

  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt

  integer :: nloc
  type(vector_t) :: initdesign
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
  call neko_field%init(neko_case%msh, neko_case%fluid%Xh, "neko_field")
  nloc = neko_field%dof%size()

  call des%init_from_components(nloc)

  if (pe_rank == 0) then
     print *, "Global number of design variables=", des%size_global()
  end if
   
  ! initialize the design
  call initdesign%init(des%size())
  initdesign = 0.5_rp
  call des%update_design(initdesign)

  ! -------------------------------------------------------------------------- !
  ! Construct the problem

  allocate(deflection)
  allocate(beamweight)
  allocate(con_1)


  call deflection%init_from_components("tip_deflection", des, 1.0_rp)
  call beamweight%init_from_components("beamweight", des, 1.0_rp)
  
!   call con_1%init_from_components("Constraint", des, 1)
  
  ! update obj and cons and sensitivities for the init design
  call deflection%update_value(des)
  call beamweight%update_value(des)

!   call obj%update_sensitivity(des)
!   call con_1%update_value(des)
!   call con_1%update_sensitivity(des)
   
!   if (pe_rank == 0) then
!      print *, "objective value for the initial design=", obj%value, &
!         "Positive=", con_1%value
!   end if

!   ! initialize the problem
!   call prob%init(parameters, des)
  
!   call prob%add_objective(obj)
!   call prob%add_constraint(con_1)

!   ! -------------------------------------------------------------------------- !
!   ! Execute the optimization
!   call optimizer_factory(opt, parameters, prob, des)

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
!   t_start = MPI_Wtime()
  
!   call opt%run(prob, des)

!   call MPI_Barrier(MPI_COMM_WORLD, ierr)
!   t_end = MPI_Wtime()



!   if (pe_rank == 0) then
!      print *, "opt%run execution time:", t_end - t_start, "seconds"
!   end if

  ! -------------------------------------------------------------------------- !
  ! Clean up the components


!   call neko_finalize(neko_case)
!   call opt%free()
!   call prob%free()
!   call des%free()

!   if (allocated(opt)) deallocate(opt)

end program usrneko