program usrneko

  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  use mpi_f08, only: MPI_Init, MPI_Wtime, MPI_COMM_WORLD
  use constraint, only: constraint_t
  use objective, only: objective_t


  use example_problem_1d_beam, only: deflection_con, beamweight_obj, stress_con

  use design_3dto1d , only: design_3dto1d_t
  use neko, only: neko_init, neko_finalize, neko_solve
  use nekotop_logger, only: nekotop_log
  use case, only: case_t
  use field, only: field_t
  use num_types, only: rp
  use vector, only: vector_t
  use matrix, only: matrix_t

  use comm, only: pe_rank, neko_comm
  use device, only: device_memcpy

  implicit none

  ! ========================================================================== !
  ! Set up distributed stress constraints

  ! number of elements with stress constraints
  integer :: num_constraints = 10

  ! number of beam sections to distribute the constraint
  integer :: num_constraint_partitions = 10
  ! ========================================================================== !

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
  class(constraint_t), allocatable :: deflection
  class(objective_t), allocatable :: beamweight
  class(constraint_t), allocatable :: tmp_constraint

  !> The optimizer (in this case mma)
  class(optimizer_t), allocatable :: opt

  integer :: nloc, i
  type(vector_t) :: initdesign
  real(kind=rp) :: t_start, t_end

  !> Stress constraints
  character(len=20) :: index_str
  integer, allocatable :: stress_global_indices(:)
  real(rp), allocatable :: stress_sigma_max(:)

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
  call nekotop_log%init(env_prefix = "NEKO_TOP")
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

  ! initialize the problem
  call prob%init(parameters, des)

  allocate(beamweight_obj::beamweight)
  allocate(deflection_con::deflection)

  select type(beamweight)
  type is (beamweight_obj)
     call beamweight%beamweight_obj_init( 1.0_rp, des)
  class default
     call neko_error("beamweight is not beamweight_obj!")
  end select

  select type(deflection)
  type is (deflection_con)
     call deflection%deflection_con_init(des)
  class default
     call neko_error("deflection is not deflection_con!")
  end select
  call prob%add_constraint(deflection)

  if (num_constraint_partitions .gt. num_constraints) then
     num_constraint_partitions = num_constraints
  end if
  if (num_constraints .gt. 0) then
     allocate(stress_global_indices(num_constraints))
     allocate(stress_sigma_max(num_constraints))
     ! Add constraints on global indices
     call fill_constraint_indices(stress_global_indices, num_constraints, &
          num_constraint_partitions, des%size_global())

     stress_sigma_max = 250e6_rp ! Same max stress for all

     ! Add each constraint to the problem
     do i = 1, size(stress_global_indices)
        allocate(stress_con::tmp_constraint)
        write(index_str, '(I0)') i

        select type(c => tmp_constraint)
        type is (stress_con)
           call c%init_stress_con("stress_con_" // trim(index_str), des, &
                stress_global_indices(i), stress_sigma_max(i))
        class default
           call neko_error("tmp_constraint is not stress_con!")
        end select

        call prob%add_constraint(tmp_constraint)

        if (allocated(tmp_constraint)) then
           deallocate(tmp_constraint)
        end if
     end do
  end if

  ! Add objectives to the problem
  call prob%add_objective(beamweight)

  call MPI_Barrier(neko_comm, ierr)
  ! -------------------------------------------------------------------------- !
  ! Perform finite difference validation
  ! update obj and cons and sensitivities for the init design
  !   call deflection%update_value(des)
  !   call beamweight%update_value(des)
  !   call deflection%update_sensitivity(des)
  !   call beamweight%update_sensitivity(des)
  !   if (pe_rank == 0) then
  !      print *, "Performing finite difference validation..."
  !   endif
  !   call finite_difference_validation(des, 1, 1.0e-6_rp)

  call prob%update_objectives(des)
  call prob%update_objective_sensitivities(des)
  call prob%update_constraints(des)
  call prob%update_constraint_sensitivities(des)

  call all_objectives%init(prob%get_n_objectives())
  call constraint_value%init(prob%get_n_constraints())

  call prob%get_all_objective_values(all_objectives)
  call prob%get_constraint_values(constraint_value)
  call prob%get_objective_value(objective_value)

  if (pe_rank == 0) then
     print *, "nobject=", prob%get_n_objectives(), "nconstraint=", &
          prob%get_n_constraints(), "total objective=", objective_value, &
          "all_objectives%x=", all_objectives%x, "constraint_value", &
          constraint_value%x
  end if

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
  call nekotop_log%free()
  call neko_finalize(neko_case)
  call opt%free()
  call prob%free()
  call des%free()

  if (allocated(opt)) deallocate(opt)
  if (allocated(beamweight)) deallocate(beamweight)
  if (allocated(deflection)) deallocate(deflection)
  call all_objectives%free()
  call constraint_value%free()
  call initdesign%free()
end program usrneko


! ========================================================================== !
! Finite difference validation subroutines
! ========================================================================== !

subroutine finite_difference_validation(des, k_test, delta)
  use comm, only: pe_rank
  use example_problem_1d_beam, only: deflection_con
  use design_3dto1d, only: design_3dto1d_t
  use num_types, only: rp
  use vector, only: vector_t

  type(design_3dto1d_t), intent(in) :: des
  integer, intent(in) :: k_test
  real(rp), intent(in) :: delta

  type(vector_t) :: designvec
  type(design_3dto1d_t) :: pert_design
  type(deflection_con) :: obj
  real(rp) :: f_original, f_perturbed, fd_derivative, analytical_derivative
  real(rp) :: error, rel_error
  integer :: n
  real(rp), allocatable :: sensitivities(:)

  ! Initialize objective
  call obj%deflection_con_init(des)

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
  call des%get_values(designvec)
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

subroutine fill_constraint_indices(stress_global_indices, num_constraints, &
     num_constraint_partitions, design_size)
  implicit none

  integer, intent(in) :: num_constraints, num_constraint_partitions, design_size
  integer, intent(out) :: stress_global_indices(num_constraints)

  integer :: base_size, remainder_size
  integer :: constraints_per_partition, remainder_constraints
  integer :: i, j, idx, partition_start, partition_end
  integer :: partition_constraints, partition_size
  integer :: cum_start

  ! --- partitioning of constraints ---
  constraints_per_partition = num_constraints / num_constraint_partitions
  remainder_constraints = mod(num_constraints, num_constraint_partitions)

  ! --- partitioning of design variables ---
  base_size = design_size / num_constraint_partitions
  remainder_size = mod(design_size, num_constraint_partitions)

  idx = 1
  cum_start = 1

  do i = 1, num_constraint_partitions
     ! how many constraints in this partition
     if (i <= remainder_constraints) then
        partition_constraints = constraints_per_partition + 1
     else
        partition_constraints = constraints_per_partition
     endif

     ! how many design variables in this partition
     if (i <= remainder_size) then
        partition_size = base_size + 1
     else
        partition_size = base_size
     endif

     partition_start = cum_start
     partition_end = partition_start + partition_size - 1

     ! pick consecutive indices in this partition
     do j = 1, partition_constraints
        if (idx > num_constraints) exit
        stress_global_indices(idx) = min(partition_start + (j-1), partition_end)
        idx = idx + 1
     end do

     cum_start = partition_end + 1
  end do

end subroutine fill_constraint_indices


