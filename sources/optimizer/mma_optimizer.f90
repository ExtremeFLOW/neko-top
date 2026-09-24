module mma_optimizer
  use optimizer, only: optimizer_t
  use problem, only : problem_t
  use mma, only: mma_t
  use problem, only: problem_t
  use num_types, only: rp
  use utils, only: neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use simulation_m, only: simulation_t
  use design, only: design_t
  use field, only: field_t
  use field_registry, only: neko_field_registry

  use vector, only: vector_t
  use matrix, only: matrix_t

  !only to print nglobal when running in parallel
  use comm, only: neko_comm, pe_rank
  use mpi_f08, only: MPI_INTEGER, mpi_sum, MPI_Allreduce

  use neko_config, only: NEKO_BCKND_DEVICE
  ! Inclusions from external dependencies and standard libraries
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit

  use math, only: copy, cmult
  use device_math, only: device_copy
  use field_math, only: field_rzero
  use neko_ext, only: reset
  use mask_ops, only: mask_exterior_const
  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use json_utils_ext, only: json_get_subdict

  use device_math, only: device_copy
  implicit none
  private
  public :: mma_optimizer_t

  ! Concrete type for MMA optimizer
  type, extends(optimizer_t) :: mma_optimizer_t

     type(mma_t) :: mma

     !> Scaling constraint_value%x and constraint_sensitivities%x.
     !! Note that the values are not updated but they are scaled when passed
     !! to the optimizer.
     !! (if auto_scale then constraint_value%x=scale else constraint_value%x=scale*constraint_value%x)
     !! When auto_scale is true, we use an adaptable scale for
     !! constraint_value%x and constraint_sensitivities%x in every iteration (variable scale factors)
     real(kind=rp) :: scale
     logical :: auto_scale

   contains

     ! Override the deferred methods
     generic :: init => init_from_json, init_from_components
     procedure, pass(this) :: init_from_json => mma_optimizer_init_from_json
     procedure, pass(this) :: init_from_components => &
          mma_optimizer_init_from_components

     procedure :: run => mma_optimizer_run
     procedure :: free => mma_optimizer_free

  end type mma_optimizer_t

contains

  !> Initialize the MMA optimizer from JSON file
  subroutine mma_optimizer_init_from_json(this, parameters, problem, design, &
       simulation, max_iterations, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), intent(in) :: simulation
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance

    character(len=1024) :: optimization_header
    character(len=1024) :: problem_header

    type(vector_t) :: x
    type(json_file) :: solver_parameters

    ! Initialize the logger
    call this%logger%init('optimization_data.csv')

    ! Write the header
    problem_header = problem%get_log_header()
    optimization_header = 'iter, ' // trim(problem_header) // &
         ', KKTmax, KKTnorm2, scaling factor'
    call this%logger%set_header(trim(optimization_header))

    x = design%get_design()

    if (pe_rank .eq. 0) then
       print *, "Initializing mma_optimizer with steady_state_problem_t."
    end if

    call json_get_subdict(parameters, "optimization.solver", solver_parameters)
    call this%mma%init(x%x, design%size(), problem%get_n_constraints(), &
         solver_parameters, this%scale, this%auto_scale)

    call this%init_from_components(problem, design, simulation, &
         max_iterations, tolerance)

    call x%free()
  end subroutine mma_optimizer_init_from_json

  !> Initialize the MMA optimizer from JSON file
  subroutine mma_optimizer_init_from_components(this, problem, design, &
       simulation, max_iterations, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), intent(in) :: simulation
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance

    call this%init_base(max_iterations, tolerance)

  end subroutine mma_optimizer_init_from_components

  !> Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, problem, design, simulation)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), intent(inout) :: simulation

    type(vector_t) :: x

    integer :: iter, ierr, nglobal, n
    real(kind=rp) :: scaling_factor

    real(kind=rp) :: objective_value
    type(vector_t) :: all_objectives
    type(vector_t) :: constraint_value
    type(vector_t) :: objective_sensitivities
    type(matrix_t) :: constraint_sensitivities

    type(vector_t) :: log_data


    n = design%size()
    call MPI_Allreduce(n, nglobal, 1, MPI_INTEGER, mpi_sum, neko_comm, ierr)

    !>initializing the scaling factor
    scaling_factor = 1.0_rp
    if (pe_rank .eq. 0) then
       print *, "max_iterations for the optimization loop = ", &
            this%max_iterations
    end if

    call simulation%run_forward()
    call problem%compute(design)

    call simulation%run_backward()
    call problem%compute_sensitivity(design)

    call problem%get_objective_value(objective_value)
    call problem%get_constraint_values(constraint_value)
    call problem%get_objective_sensitivities(objective_sensitivities)
    call problem%get_constraint_sensitivities(constraint_sensitivities)
    call problem%get_all_objective_values(all_objectives)

    ! Stamp the initial condition
    call mma_logger_assemble_data(log_data, 0, objective_value, &
         all_objectives, constraint_value, 0.0_rp, 0.0_rp, scaling_factor, &
         problem%get_n_objectives(), problem%get_n_constraints())
    call this%logger%write(log_data)

    call simulation%write(0)
    call design%write(0)

    do iter = 1, this%max_iterations
       if (this%mma%get_residumax() .lt. this%tolerance) exit

       ! Scaling
       if (this%auto_scale .eqv. .true.) then
          scaling_factor = abs(this%scale/constraint_value%x(1))
       else
          scaling_factor = abs(this%scale)
       end if

       x = design%get_design()

       ! Use scaled sensitivities to update the design variable
       call this%mma%update(iter, x%x, objective_sensitivities, &
            scaling_factor * constraint_value, &
            scaling_factor * constraint_sensitivities)

       call design%update_design(x)

       call simulation%run_forward()
       call problem%compute(design)

       call simulation%run_backward()
       call problem%compute_sensitivity(design)

       call problem%get_objective_value(objective_value)
       call problem%get_constraint_values(constraint_value)
       call problem%get_objective_sensitivities(objective_sensitivities)
       call problem%get_constraint_sensitivities(constraint_sensitivities)
       call problem%get_all_objective_values(all_objectives)

       call this%mma%KKT(x%x, objective_sensitivities, &
            constraint_value, constraint_sensitivities)

       ! Stamp the i^th iteration
       call mma_logger_assemble_data(log_data, iter, objective_value, &
            all_objectives, constraint_value, this%mma%get_residumax(), &
            this%mma%get_residunorm(), scaling_factor, &
            problem%get_n_objectives(), problem%get_n_constraints())
       call this%logger%write(log_data)

       call simulation%write(iter)
       call design%write(iter)
       call simulation%reset()
    end do

    ! Final state after optimization
    if (pe_rank .eq. 0) then
       print *, "MMA Optimization completed after", iter-1, "iterations."
    end if

    call constraint_value%free()
    call objective_sensitivities%free()
    call constraint_sensitivities%free()

  end subroutine mma_optimizer_run

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%mma%free()
  end subroutine mma_optimizer_free

  ! package up the log data
  subroutine mma_logger_assemble_data(log_data, iter, objective_value, &
       all_objectives, constraint_value, residumax, residunorm, &
       scaling_factor, n, m)
    type(vector_t), intent(out) :: log_data
    integer, intent(in) :: iter
    real(kind=rp), intent(in) ::objective_value
    type(vector_t), intent(in) :: all_objectives
    type(vector_t), intent(in) :: constraint_value
    real(kind=rp), intent(in) :: residumax, residunorm, scaling_factor
    integer, intent(in) :: n, m
    integer :: i_tmp1, i_tmp2

    ! initialize the logger data
    ! iter | tot F | F_1 | .. |F_n | C_1 | ... | C_n | KKT | KKT2 | scale |
    call log_data%init(5 + n + m)

    ! iteration
    log_data%x(1) = real(iter, kind=rp)

    ! total objective
    log_data%x(2) = objective_value

    ! individual objectives
    i_tmp1 = 3
    i_tmp2 = i_tmp1 + n - 1
    log_data%x(i_tmp1 : i_tmp2) = all_objectives%x

    ! constraints
    i_tmp1 = i_tmp2 + 1
    i_tmp2 = i_tmp1 + m - 1
    log_data%x(i_tmp1 : i_tmp2) = constraint_value%x

    ! convergence stuff
    log_data%x(i_tmp2 + 1) = residumax
    log_data%x(i_tmp2 + 2) = residunorm
    log_data%x(i_tmp2 + 3) = scaling_factor

  end subroutine mma_logger_assemble_data
end module mma_optimizer

