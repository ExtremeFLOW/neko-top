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
  use case, only: case_t
  use field, only: field_t
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t

  use comm, only: pe_rank
  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use neko_config, only: NEKO_BCKND_DEVICE
  
  use device_mma_math, only: device_solve_linear_system, device_gpu_solve_system
  use iso_c_binding

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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  type(vector_t) :: bvec, bvec_original, rhs
  type(matrix_t) :: Amat
  integer :: info, systemsize, i, j
  real(kind=rp) :: time_custom, time_cusolver, error_custom, error_cusolver
  integer, parameter :: sizes(8) = [10, 50, 100, 500, 1000, 2000, 5000, 10000]
 
  call MPI_Init(ierr)

do i = 1, size(sizes)
   systemsize = sizes(i)
   call bvec%init(systemsize)
   call bvec_original%init(systemsize)  ! Store original RHS
   call Amat%init(systemsize, systemsize)
   
   ! Initialize matrix
   do j = 1, systemsize
      Amat%x(j,j) = 0.5_rp 
      bvec%x(j) = 1.0_rp
      bvec_original%x(j) = 1.0_rp  ! Store original RHS
   end do
   
   call device_memcpy(Amat%x, Amat%x_d, systemsize*systemsize, HOST_TO_DEVICE, sync=.true.)
   call device_memcpy(bvec%x, bvec%x_d, systemsize, HOST_TO_DEVICE, sync=.true.)
   
   ! Time custom solver
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_start = MPI_Wtime()
   call device_solve_linear_system(Amat%x_d, bvec%x_d, systemsize, info)
   call device_memcpy(bvec%x, bvec%x_d, systemsize, DEVICE_TO_HOST, sync=.true.)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_end = MPI_Wtime()
   time_custom = t_end - t_start
   
   ! Correct error calculation
   rhs%x = matmul(Amat%x, bvec%x)  ! A*x
   rhs%x = bvec_original%x - rhs%x  ! b - A*x
   error_custom = maxval(abs(rhs%x))
   
   print *, "Size:", systemsize, "custom info=", info, "error:", error_custom
   
   ! Reset for cuSOLVER
   do j = 1, systemsize
      Amat%x(j,j) = 0.5_rp 
      bvec%x(j) = 1.0_rp
   end do
   
   call device_memcpy(Amat%x, Amat%x_d, systemsize*systemsize, HOST_TO_DEVICE, sync=.true.)
   call device_memcpy(bvec%x, bvec%x_d, systemsize, HOST_TO_DEVICE, sync=.true.)
   
   ! Time cuSOLVER
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_start = MPI_Wtime()
   call device_gpu_solve_system(Amat%x_d, bvec%x_d, systemsize, info)
   call device_memcpy(bvec%x, bvec%x_d, systemsize, DEVICE_TO_HOST, sync=.true.)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   t_end = MPI_Wtime()
   time_cusolver = t_end - t_start
   
   ! Correct error calculation
   rhs%x = matmul(Amat%x, bvec%x)  ! A*x
   rhs%x = bvec_original%x - rhs%x  ! b - A*x
   error_cusolver = maxval(abs(rhs%x))
   
   print *, "Size:", systemsize, "cuSOLVER info=", info, "error:", error_cusolver
   
   if (pe_rank == 0) then
      print *, "Size:", systemsize, "Custom:", time_custom, "cuSOLVER:", time_cusolver, "Ratio:", time_cusolver/time_custom
      print *, "---------------------------------------------------------------"
   end if
end do
call neko_error('Done comparing cuSOLVER and my custom solver!')


  ! -------------------------------------------------------------------------- !
  ! Initialize the MPI environment


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

  call neko_finalize(neko_case)
  call opt%free()
  call prob%free()
  call des%free()

  if (allocated(opt)) deallocate(opt)

end program usrneko
