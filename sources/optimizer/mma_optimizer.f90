module mma_optimizer
  use optimizer, only: optimizer_t
  use mma, only: mma_t
  use problem, only: problem_t
  use num_types, only: rp
  use utils, only: neko_error
  use json_utils, only: json_get, json_get_or_default
  use simulation_m, only: simulation_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use constraint, only: constraint_t
  use dummy_constraint, only: dummy_constraint_t

  ! External modules
  use json_module, only: json_file
  use vector, only: vector_t
  use matrix, only: matrix_t
  use comm, only: pe_rank
  use neko_config, only: NEKO_BCKND_DEVICE
  use scratch_registry, only: neko_scratch_registry
  use profiler, only: profiler_start_region, profiler_end_region
  use logger, only: neko_log
  use csv_file, only: csv_file_t
  use vector_math, only: vector_cmult
  use matrix_math, only: matrix_cmult
  use device, only: device_memcpy, HOST_TO_DEVICE
  implicit none
  private

  public :: mma_optimizer_t

  ! Concrete type for MMA optimizer
  type, extends(optimizer_t) :: mma_optimizer_t

     type(mma_t) :: mma

     !> Scaling constraint_value%x and constraint_sensitivities%x.
     !! (if auto_scale then constraint_value%x=scale else
     !! constraint_value%x=scale*constraint_value%x)
     !! When auto_scale is true, we use an adaptable scale for
     !! constraint_value%x and constraint_sensitivities%x
     !! in every iteration (variable scale factors)
     real(kind=rp) :: scale = 1.0_rp
     real(kind=rp) :: scaling_factor = 1.0_rp
     logical :: auto_scale = .false.

     ! Set to flase to remove logging for optimal performance
     logical :: unconstrained_problem = .false.

     !> A file writer to document the convergence history
     logical :: enable_output = .true.
     type(csv_file_t) :: csv_log
   contains

     ! Override the deferred methods
     generic :: init => init_from_json, init_from_components
     procedure, pass(this) :: init_from_json => mma_optimizer_init_from_json
     procedure, pass(this) :: init_from_components => &
          mma_optimizer_init_from_components

     procedure, pass(this) :: run => mma_optimizer_run
     procedure, pass(this) :: validate => mma_optimizer_validate
     procedure, pass(this) :: write => mma_optimizer_write
     procedure, pass(this) :: free => mma_optimizer_free

  end type mma_optimizer_t

contains

  !> Initialize the MMA optimizer from JSON file
  subroutine mma_optimizer_init_from_json(this, parameters, problem, design, &
       simulation)
    class(mma_optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), optional, intent(in) :: simulation

    ! Variables for settings
    type(json_file) :: solver_parameters
    logical :: enable_output
    integer :: max_iterations
    real(kind=rp) :: tolerance

    ! Read the solver properties from the JSON file
    call json_get(parameters, 'optimization.solver', solver_parameters)
    call json_get_or_default(solver_parameters, 'max_iterations', &
         max_iterations, 100)
    call json_get_or_default(solver_parameters, 'tolerance', &
         tolerance, 1.0e-3_rp)
    call json_get_or_default(solver_parameters, 'enable_output', &
         enable_output, .true.)

    call this%init_from_components(problem, design, max_iterations, tolerance, &
         enable_output, solver_parameters, simulation)

  end subroutine mma_optimizer_init_from_json

  !> Initialize the MMA optimizer from JSON file
  subroutine mma_optimizer_init_from_components(this, problem, design, &
       max_iterations, tolerance, enable_output, solver_parameters, simulation)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(in) :: design
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance
    logical, intent(in) :: enable_output
    type(json_file), intent(inout), optional :: solver_parameters
    type(simulation_t), intent(in), optional :: simulation

    ! Local variables
    type(vector_t), pointer :: x
    integer :: ind

    ! Local variables
    class(constraint_t), allocatable :: dummy_con

    call neko_log%section('Optimizer Initialization')

    ! Check if the problem is unconstrained
    this%unconstrained_problem = (problem%get_n_constraints() .eq. 0)
    if (this%unconstrained_problem) then
       call neko_log%message('Unconstrained problem detected. ' // &
            'Adding a dummy constraint to enable MMA optimization.')

       allocate(dummy_constraint_t::dummy_con)
       select type (con => dummy_con)
       type is (dummy_constraint_t)
          call con%init_from_attributes(design)
       end select

       call problem%add_constraint(dummy_con)
    end if

    ! Initialize mma_t, handling the dummy_constraint added for unconstrained
    ! problems in mma_optimizer_run()
    call neko_scratch_registry%request(x, ind, design%size(), .false.)

    call design%get_values(x)
    call this%mma%init(x, design%size(), problem%get_n_constraints(), &
         solver_parameters, this%scale, this%auto_scale)

    call neko_scratch_registry%relinquish_vector(ind)

    !set the enable_output flag
    this%enable_output = enable_output
    this%scaling_factor = this%scale

    ! Initialize the logger
    if (this%enable_output) then
       call this%csv_log%init('optimization_data.csv')
    end if

    call this%init_base(max_iterations, tolerance)

    call neko_log%end_section()

  end subroutine mma_optimizer_init_from_components

  !> Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, problem, design, simulation)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation

    type(vector_t), pointer :: x

    integer :: iter, n

    real(kind=rp) :: objective_value
    type(vector_t), pointer :: constraint_value
    type(vector_t), pointer :: objective_sensitivities
    type(matrix_t), pointer :: constraint_sensitivities
    integer :: ind(4)

    n = design%size()

    ! Initialize the vectors
    call neko_scratch_registry%request(x, ind(1), &
         n, .false.)
    call neko_scratch_registry%request(constraint_value, ind(2), &
         problem%get_n_constraints(), .false.)
    call neko_scratch_registry%request(objective_sensitivities, ind(3), &
         n, .false.)
    call neko_scratch_registry%request(constraint_sensitivities, ind(4), &
         problem%get_n_constraints(), n, .false.)

    !>initializing the scaling factor
    if (pe_rank .eq. 0) then
       print *, 'max_iterations for the optimization loop = ', &
            this%max_iterations
    end if
    call profiler_start_region('Optimizer iteration')

    call problem%compute(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_forward(0)
    end if
    call problem%compute_sensitivity(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_adjoint(0)
    end if

    call problem%get_objective_value(objective_value)
    call problem%get_constraint_values(constraint_value)

    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(objective_sensitivities)
    class default
       call problem%get_objective_sensitivities(objective_sensitivities)
    end select

    call problem%get_constraint_sensitivities(constraint_sensitivities)

    call profiler_end_region('Optimizer iteration')

    call this%write(0, problem)
    call design%write(0)

    do iter = 1, this%max_iterations
       if (this%mma%get_residumax() .lt. this%tolerance) exit

       call profiler_start_region('Optimizer iteration')

       ! Scaling
       if (this%auto_scale) then
          this%scaling_factor = abs(this%scale / constraint_value%x(1))
       end if

       call vector_cmult(constraint_value, this%scaling_factor)
       call matrix_cmult(constraint_sensitivities, this%scaling_factor)

       ! Use scaled sensitivities to update the design variable
       call design%get_values(x)
       call this%mma%update(iter, x, objective_sensitivities, &
            constraint_value, constraint_sensitivities)

       call design%update_design(x)

       call problem%compute(design, simulation)
       if (present(simulation) .and. this%enable_output) then
          call simulation%write_forward(iter)
       end if
       call problem%compute_sensitivity(design, simulation)
       if (present(simulation) .and. this%enable_output) then
          call simulation%write_adjoint(iter)
       end if

       call problem%get_constraint_values(constraint_value)

       select type (des => design)
       type is (brinkman_design_t)
          call des%get_sensitivity(objective_sensitivities)
       class default
          call problem%get_objective_sensitivities(objective_sensitivities)
       end select

       call problem%get_constraint_sensitivities(constraint_sensitivities)

       call this%mma%KKT(x, objective_sensitivities, &
            constraint_value, constraint_sensitivities)

       call profiler_end_region('Optimizer iteration')

       ! Log the progress and outputs
       call this%write(iter, problem)
       call design%write(iter)
    end do

    call this%validate(problem, design)

    ! Final state after optimization
    if (pe_rank .eq. 0) then
       print *, 'MMA Optimization completed after', iter-1, 'iterations.'
    end if

    ! Free local resources
    call neko_scratch_registry%relinquish(ind)

  end subroutine mma_optimizer_run

  !> Validate the solution for the MMA optimizer
  subroutine mma_optimizer_validate(this, problem, design)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design

    type(vector_t), pointer :: constraint_values
    integer :: ind

    call neko_scratch_registry%request( &
         constraint_values, ind, problem%get_n_constraints(), .false.)

    call problem%get_constraint_values(constraint_values)
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(constraint_values%x, constraint_values%x_d, &
            constraint_values%size(), HOST_TO_DEVICE, .true.)
    end if

    if (any(constraint_values%x .gt. 0.0_rp)) then
       call neko_error('MMA optimizer validation failed: ' // &
            'Constraints are not satisfied.')
    end if

    ! Free local resources
    call neko_scratch_registry%relinquish_vector(ind)

  end subroutine mma_optimizer_validate

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%free_base()
    call this%mma%free()
  end subroutine mma_optimizer_free

  subroutine mma_optimizer_write(this, iter, problem)
    class(mma_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(in) :: problem

    type(vector_t), pointer :: log_data
    type(vector_t), pointer :: all_objectives
    type(vector_t), pointer :: constraint_value
    real(kind=rp) :: objective_value
    character(len=1024) :: header

    integer :: log_size, ind(3), n, m, i_tmp1, i_tmp2

    if (.not. this%enable_output) return
    call profiler_start_region('Optimizer logging')

    n = problem%get_n_objectives()
    m = problem%get_n_constraints()
    if (this%unconstrained_problem) then
       log_size = 5 + n
    else
       log_size = 5 + n + m
    endif

    call neko_scratch_registry%request(log_data, ind(1), log_size, .false.)
    call neko_scratch_registry%request(all_objectives, ind(2), n, .false.)
    call neko_scratch_registry%request(constraint_value, ind(3), m, .false.)

    if (iter .eq. 0) then
       header = 'iter, ' // &
            trim(problem%get_log_header()) // &
            ', KKTmax, KKTnorm2, scaling factor, ' // &
            trim(this%mma%get_backend_and_subsolver())

       call this%csv_log%set_header(trim(header))
    end if

    ! Prepare data for logging
    call problem%get_objective_value(objective_value)
    call problem%get_all_objective_values(all_objectives)
    call problem%get_constraint_values(constraint_value)

    ! Assemble the log data
    log_data%x(1) = real(iter, kind=rp)

    ! total objective
    log_data%x(2) = objective_value

    ! individual objectives
    i_tmp1 = 3
    i_tmp2 = i_tmp1 + n - 1
    log_data%x(i_tmp1 : i_tmp2) = all_objectives%x

    ! constraints
    if (.not. this%unconstrained_problem) then
       i_tmp1 = i_tmp2 + 1
       i_tmp2 = i_tmp1 + m - 1
       log_data%x(i_tmp1 : i_tmp2) = constraint_value%x
    end if

    ! convergence stuff
    if (iter .eq. 0) then
       log_data%x(i_tmp2 + 1) = 0.0_rp
       log_data%x(i_tmp2 + 2) = 0.0_rp
    else
       log_data%x(i_tmp2 + 1) = this%mma%get_residumax()
       log_data%x(i_tmp2 + 2) = this%mma%get_residunorm()
    end if
    log_data%x(i_tmp2 + 3) = this%scaling_factor

    call this%csv_log%write(log_data)

    ! Free local resources
    call neko_scratch_registry%relinquish_vector(ind)

    call profiler_end_region('Optimizer logging')
  end subroutine mma_optimizer_write

end module mma_optimizer

