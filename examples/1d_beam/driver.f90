program usrneko

  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  use mpi_f08, only: MPI_Init,  MPI_Wtime, MPI_COMM_WORLD, mpi_in_place, &
       mpi_allreduce, mpi_exscan


  use example_problem, only: mma_obj, beamweight_obj, stress_con

  use design_3dto1d , only: design_3dto1d_t
  use neko, only: neko_init, neko_finalize, neko_solve
  use case, only: case_t
  use field, only: field_t
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t

  use comm, only: pe_rank, neko_comm
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
  type(stress_con), allocatable :: stress_constraints(:)

  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt

  integer :: nloc, i
  type(vector_t) :: initdesign
  real(kind=rp) :: t_start, t_end

  !> Stress constraints
  character(len=20) :: index_str
  integer :: stress_global_indices(9)
  real(rp) :: stress_sigma_max(9)

  !> For getting objectives and constraints values though getters in problem_t
  type(vector_t) :: all_objectives, constraint_value
  real(rp) :: objective_value

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

  ! add constraints on global indices
  stress_global_indices = [1, 2, 3, 3, 4, 5, des%size_global()-1 , des%size_global(), 150000]
  stress_sigma_max = 250e6_rp  ! Same max stress for all
  
  allocate(stress_constraints(size(stress_global_indices)))

  call deflection%init_from_components("tip_deflection", des, 1.0_rp)
  call beamweight%init_from_components("beamweight", des, 1.0_rp)
  
  ! Add each constraint to the problem
  do i = 1, size(stress_constraints)
     write(index_str, '(I0)') stress_global_indices(i)  ! I0 format for minimal length
     !   write(index_str, '(I0)') i  ! I0 format for minimal length
     call stress_constraints(i)%init_stress_con( &
          "stress_con_"//trim(index_str), &
          des, stress_global_indices(i), stress_sigma_max(i))

     call prob%add_constraint(stress_constraints(i))
  end do

  ! Add objectives to the problem
  call prob%add_objective(deflection)
  call prob%add_objective(beamweight)


  call MPI_Barrier(neko_comm, ierr)
  if (pe_rank == 0) print *, "Constraints distribution complete!"
  if (pe_rank == 0) print *, "number of problem objectives=", prob%get_n_objectives(), "constraints=", prob%get_n_constraints()
  
  ! update obj and cons and sensitivities for the init design
   !   call deflection%update_value(des)
   !   call beamweight%update_value(des)
   !   call deflection%update_sensitivity(des)
   !   call beamweight%update_sensitivity(des)
  ! -------------------------------------------------------------------------- !
  ! Perform finite difference validation
  !   if (pe_rank == 0) then
  !      print *, "Performing finite difference validation..."
  !   endif
  !   call finite_difference_validation(des, 1, 1.0e-6_rp)

  ! Update values and sensitivities
  !   do i = 1, size(stress_constraints)
  !      call stress_constraints(i)%update_value(des)
  !      call stress_constraints(i)%update_sensitivity(des)
  !   end do
   call prob%update_objectives(des)
   call prob%update_objective_sensitivities(des)
   call prob%update_constraints(des)
   call prob%update_constraint_sensitivities(des)

  !   call con_1%update_value(des)
  !   call con_1%update_sensitivity(des)
  call prob%get_all_objective_values(all_objectives)
  call prob%get_constraint_values(constraint_value)
  call prob%get_objective_value(objective_value)
  if (pe_rank == 0) then
     print *, "nobject=", prob%get_n_objectives(), "nconstraint=", prob%get_n_constraints(), &
        "total objective=", objective_value, &
        "all_objectives%x=", all_objectives%x, "constraint_value", constraint_value%x
  end if

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



! ========================================================================== !
! Finite difference validation subroutines
! ========================================================================== !

subroutine finite_difference_validation(des, k_test, delta)
  use mpi_f08, only: MPI_Allreduce, MPI_SUM, MPI_DOUBLE_PRECISION
  use comm, only: pe_rank, pe_size, neko_comm
  use example_problem, only: mma_obj
  use design_3dto1d, only: design_3dto1d_t
  use num_types, only: rp
  use vector, only: vector_t
  
  type(design_3dto1d_t), intent(in) :: des
  integer, intent(in) :: k_test
  real(rp), intent(in) :: delta
  
  type(vector_t) :: designvec
  type(design_3dto1d_t) :: pert_design
  type(mma_obj) :: obj
  real(rp) :: f_original, f_perturbed, fd_derivative, analytical_derivative
  real(rp) :: error, rel_error
  integer :: n, i, ierr
  real(rp), allocatable :: sensitivities(:)
  
  ! Initialize objective
  call obj%init_from_components("test_obj", des, 1.0_rp)
  
  ! Get original value and sensitivities
  call obj%update_value(des)
  call obj%update_sensitivity(des)
  f_original = obj%value
  
  n = des%size()
  allocate(sensitivities(n))
  sensitivities = obj%sensitivity%x
  
  ! Create perturbed design
  call pert_design%init_from_components(n)
  call designvec%init(n)
  designvec = des%get_values()
  if (k_test >= 1 .and. k_test <= n) then
    designvec%x(k_test) = designvec%x(k_test) + delta
  endif
  call pert_design%update_design(designvec)


  ! Compute perturbed value
  call obj%update_value(pert_design)
  f_perturbed = obj%value
  
  ! Finite difference derivative
  fd_derivative = (f_perturbed - f_original) / delta
  analytical_derivative = sensitivities(k_test)
  
  ! Calculate error
  error = abs(fd_derivative - analytical_derivative)
  if (abs(analytical_derivative) > 1e-12) then
    rel_error = error / abs(analytical_derivative)
  else
    rel_error = error
  endif
  
  ! Output results
  if (pe_rank == 0) then
    print *, "=============================================="
    print *, "FINITE DIFFERENCE VALIDATION"
    print *, "=============================================="
    print *, "Test variable index:      ", k_test
    print *, "Perturbation size (delta):", delta
    print *, "Original function value:  ", f_original
    print *, "Perturbed function value: ", f_perturbed
    print *, "Finite difference deriv:  ", fd_derivative
    print *, "Analytical derivative:    ", analytical_derivative
    print *, "Absolute error:           ", error
    print *, "Relative error:           ", rel_error
    print *, "=============================================="
  endif
  
  deallocate(sensitivities)
  call obj%free()
end subroutine finite_difference_validation