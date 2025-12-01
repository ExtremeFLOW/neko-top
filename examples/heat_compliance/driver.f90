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
  class(objective_t), allocatable :: heat_compliance
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

  allocate(heat_compliance_t :: heat_compliance)

  select type(heat_compliance)
  type is (heat_compliance_t)
     call heat_compliance%init_from_attributes(des, neko_case%fluid%c_Xh, parameters)
  class default
     call neko_error("??!")
  end select
  ! Add objectives to the problem
  call prob%add_objective(heat_compliance)

  allocate(thermal_volume_constraint_t :: volume_constraint)

  select type(volume_constraint)
  type is (thermal_volume_constraint_t)
     call volume_constraint%init_from_components(des, neko_case%fluid%c_Xh, &
     name   = 'Volume constraint', &
     is_max = .true., &
     limit  = 0.4_rp)
  class default
     call neko_error("??!")
  end select
  ! Add constraint to the problem
  call prob%add_constraint(volume_constraint)

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
  !call prob%update_constraints(des)
  !call prob%update_constraint_sensitivities(des)


  call prob%get_all_objective_values(all_objectives)
  !call prob%get_constraint_values(constraint_value)
  call prob%get_objective_value(objective_value)

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
