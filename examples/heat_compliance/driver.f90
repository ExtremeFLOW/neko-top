program usrneko

  use problem, only: problem_t
  use optimizer, only: optimizer_t, optimizer_factory
  use objective, only: objective_t
  use design, only: design_t

  use json_module, only: json_file
  use utils, only: neko_error
  use json_utils_ext, only: json_read_file

  use mpi_f08, only: MPI_Init, MPI_Wtime, MPI_COMM_WORLD
  use constraint, only: constraint_t


  use heat_compliance, only: heat_compliance_t, thermal_conductivity_design_t, thermal_volume_constraint_t

  use design_3dto1d , only: design_3dto1d_t
  use neko, only: neko_init, neko_finalize, neko_solve
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
  class(design_t), allocatable :: des

  !> The problem type
  type(problem_t) :: prob
  class(objective_t), allocatable :: heat_obj
  class(constraint_t), allocatable :: volume_constraint
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
  call neko_field%init(neko_case%msh, neko_case%fluid%Xh, "neko_field")
  nloc = neko_field%dof%size()

  allocate(thermal_conductivity_design_t :: des)
  select type(des)
  type is (thermal_conductivity_design_t)
     call des%init(parameters, neko_case%fluid%c_Xh)
  class default
     call neko_error("??!")
  end select
  

  ! -------------------------------------------------------------------------- !
  ! Construct the problem

  ! initialize the problem
  call prob%init(parameters, des)

  allocate(heat_compliance_t :: heat_obj)

  select type(heat_obj)
  type is (heat_compliance_t)
     call heat_obj%init_from_attributes(des, neko_case%fluid%c_Xh, parameters)
  class default
     call neko_error("??!")
  end select

!   ! -------------------------------------------------------------------------- !
!   ! Perform finite difference validation
!   ! update obj and cons and sensitivities for the init design
!     if (pe_rank == 0) then
!        print *, "Performing finite difference validation..."
!     endif
!     select type(des)
!   type is (thermal_conductivity_design_t)
!     call finite_difference_validation(des, 110, 1.0e-6_rp, neko_case%fluid%c_Xh, parameters)
!   class default
!      call neko_error("??!")
!   end select

!    call neko_error("done")
!   ! -------------------------------------------------------------------------- !
  ! Add objectives to the problem
  call prob%add_objective(heat_obj)

  allocate(thermal_volume_constraint_t :: volume_constraint)

  select type(volume_constraint)
  type is (thermal_volume_constraint_t)
     call volume_constraint%init_from_components(des, neko_case%fluid%c_Xh, &
     name   = 'Volume constraint', &
     is_max = .true., &
     limit  = 0.15_rp)
  class default
     call neko_error("??!")
  end select
  ! Add constraint to the problem
  call prob%add_constraint(volume_constraint)

  call MPI_Barrier(neko_comm, ierr)


  !call prob%update_objectives(des)
  !call prob%update_objective_sensitivities(des)
  !call prob%update_constraints(des)
  !call prob%update_constraint_sensitivities(des)


  !call prob%get_all_objective_values(all_objectives)
  !call prob%get_constraint_values(constraint_value)
  !call prob%get_objective_value(objective_value)

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

  ! ========================================================================== !
! Finite difference validation subroutines
! ========================================================================== !

subroutine finite_difference_validation(des, k_test, delta, coef, parameters)
  use comm, only: pe_rank
  use heat_compliance, only: heat_compliance_t, thermal_conductivity_design_t, thermal_volume_constraint_t
  use num_types, only: rp
  use vector, only: vector_t
  use coefs, only: coef_t
  use json_module, only: json_file
  use design, only: design_t

  type(thermal_conductivity_design_t), intent(in) :: des
  integer, intent(in) :: k_test
  real(rp), intent(in) :: delta
  type(coef_t), intent(in) :: coef
  type(json_file), intent(inout) :: parameters

  type(vector_t) :: designvec
  type(thermal_conductivity_design_t) :: pert_design
  type(heat_compliance_t) :: obj
  real(rp) :: f_original, f_perturbed, fd_derivative, analytical_derivative
  real(rp) :: error, rel_error
  integer :: n
  real(rp), allocatable :: sensitivities(:)

  ! Initialize objective
  call obj%init_from_attributes(des, coef, parameters)

  ! Get original value and sensitivities
  call obj%update_value(des)
  call obj%update_sensitivity(des)
  f_original = obj%value

  n = des%size()
  allocate(sensitivities(n))
  sensitivities = obj%sensitivity%x

  ! Create perturbed design
  call pert_design%init(parameters, coef)
  call designvec%init(n)
  call des%get_values(designvec)
  if (pe_rank == 0) then !only have one rank perturb
  if (k_test >= 1 .and. k_test <= n) then
     designvec%x(k_test) = designvec%x(k_test) + delta
  endif
  end if
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


