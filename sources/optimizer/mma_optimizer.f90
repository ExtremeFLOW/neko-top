!> @file mma_optimizer.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.

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
  use math, only: abscmp
  use profiler, only: profiler_start_region, profiler_end_region
  use logger, only: neko_log
  use vector_math, only: vector_cmult
  use matrix_math, only: matrix_cmult
  use device, only: device_memcpy, DEVICE_TO_HOST
  use scratch_registry, only: neko_scratch_registry
  use comm, only: pe_rank, NEKO_COMM
  use mpi_f08, only: MPI_Barrier

  implicit none
  private

  public :: mma_optimizer_t

  ! Concrete type for MMA optimizer
  type, extends(optimizer_t) :: mma_optimizer_t

     type(mma_t), private :: mma

     !> Scaling constraint_value%x and constraint_sensitivities%x.
     !! (if auto_scale then constraint_value%x=scale else
     !! constraint_value%x=scale*constraint_value%x)
     !! When auto_scale is true, we use an adaptable scale for
     !! constraint_value%x and constraint_sensitivities%x
     !! in every iteration (variable scale factors)
     real(kind=rp), private :: scale = 1.0_rp
     real(kind=rp), private :: scaling_factor = 1.0_rp
     logical, private :: auto_scale = .false.
     real(kind=rp) :: tolerance = 0.0_rp

     ! Set to flags to remove logging for optimal performance
     logical, private :: unconstrained_problem = .false.

     !> A file writer to document the convergence history
     logical, private :: enable_output = .true.
   contains

     ! Override the deferred methods
     generic :: init => init_from_json, init_from_components
     procedure, pass(this) :: init_from_json => mma_optimizer_init_from_json
     procedure, pass(this) :: init_from_components => &
          mma_optimizer_init_from_components

     procedure, pass(this) :: initialize => mma_optimizer_initialize
     procedure, pass(this) :: step => mma_optimizer_step
     procedure, pass(this) :: validate => mma_optimizer_validate
     procedure, pass(this) :: write => mma_optimizer_write
     procedure, pass(this) :: free => mma_optimizer_free

     procedure, pass(this) :: save_checkpoint_components => &
          mma_optimizer_save_checkpoint_components
     procedure, pass(this) :: load_checkpoint_components => &
          mma_optimizer_load_checkpoint_components

  end type mma_optimizer_t

contains

  ! -------------------------------------------------------------------------- !
  ! Allocator and deallocator methods for the MMA optimizer

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
    call this%read_base_settings(solver_parameters)

    call this%init_from_components(problem, design, max_iterations, tolerance, &
         enable_output, solver_parameters, simulation)

  end subroutine mma_optimizer_init_from_json

  !> Initialize the MMA optimizer from JSON file
  subroutine mma_optimizer_init_from_components(this, problem, design, &
       max_iterations, tolerance, enable_output, &
       solver_parameters, simulation)
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
    character(len=32) :: extra_headers(3)
    class(constraint_t), allocatable :: dummy_con

    call neko_log%section('Optimizer Initialization')

    ! Check if the problem is unconstrained
    this%unconstrained_problem = problem%get_n_constraints() .eq. 0
    if (this%unconstrained_problem) then
       call neko_log%message('Unconstrained problem detected. ' // &
            'Switching to explicit closed-form MMA subsolver and KKT.')

       allocate(dummy_constraint_t::dummy_con)
       select type (con => dummy_con)
       type is (dummy_constraint_t)
          call con%init_from_attributes(design)
       end select

       call problem%add_constraint(dummy_con)
       if (allocated(dummy_con)) deallocate(dummy_con)
    end if

    ! Initialize mma_t, handling the dummy_constraint added for unconstrained
    ! problems in mma_optimizer_run()
    call neko_scratch_registry%request(x, ind, design%size(), .false.)

    call design%get_values(x)
    call this%mma%init(x, design%size(), problem%get_n_constraints(), &
         solver_parameters, this%scale, this%auto_scale, &
         this%unconstrained_problem)

    call neko_scratch_registry%relinquish(ind)

    !set the enable_output flag
    this%enable_output = enable_output
    this%scaling_factor = this%scale
    this%tolerance = tolerance

    ! Initialize the logger
    if (this%enable_output) then
       extra_headers(1) = 'KKTmax'
       extra_headers(2) = 'KKTnorm2'
       extra_headers(3) = 'scaling factor'
       call this%init_log(problem, extra_headers = extra_headers, &
            include_constraints = .not. this%unconstrained_problem, &
            filename = 'optimization_data.csv')
    end if

    call this%init_base('MMA', max_iterations)

    call neko_log%end_section()

  end subroutine mma_optimizer_init_from_components

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%free_base()
    call this%mma%free()
  end subroutine mma_optimizer_free

  ! -------------------------------------------------------------------------- !
  ! Implementation of the deferred methods for the MMA optimizer

  !> Prepare the MMA optimizer before starting the optimization loop
  subroutine mma_optimizer_initialize(this, problem, design, simulation)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation

    type(vector_t), pointer :: x
    type(vector_t), pointer :: constraint_value
    type(vector_t), pointer :: objective_sensitivities
    type(matrix_t), pointer :: constraint_sensitivities
    integer :: n_design, n_constraint, indices(4)

    n_design = design%size()
    n_constraint = problem%get_n_constraints()

    ! Grab some local pointers
    call neko_scratch_registry%request(x, indices(1), n_design, .false.)
    call neko_scratch_registry%request(constraint_value, indices(2), &
         n_constraint, .false.)
    call neko_scratch_registry%request(objective_sensitivities, indices(3), &
         n_design, .false.)
    call neko_scratch_registry%request(constraint_sensitivities, indices(4), &
         n_constraint, n_design, .false.)

    ! Evaluate the problem based on the updated design
    call problem%compute(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_forward(0)
    end if
    call problem%compute_sensitivity(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_adjoint(0)
    end if

    ! Retrieve the updated objective and constraint values and sensitivities
    call design%get_values(x)
    call problem%get_constraint_values(constraint_value)
    call problem%get_constraint_sensitivities(constraint_sensitivities)

    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(objective_sensitivities)
       ! Convert gradient to directional derivative
       call des%project_sensitivity(objective_sensitivities)
       call des%project_sensitivity(constraint_sensitivities)
    class default
       call problem%get_objective_sensitivities(objective_sensitivities)
    end select

    ! Check the KKT conditions and check for convergence
    call this%mma%KKT(x, objective_sensitivities, &
         constraint_value, constraint_sensitivities)

    call neko_scratch_registry%relinquish(indices)
  end subroutine mma_optimizer_initialize

  !> Function for computing a step in the optimization loop
  function mma_optimizer_step(this, iter, problem, design, simulation) &
       result(converged)
    class(mma_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation

    type(vector_t), pointer :: x
    type(vector_t), pointer :: constraint_value
    type(vector_t), pointer :: objective_sensitivities
    type(matrix_t), pointer :: constraint_sensitivities
    integer :: n_design, n_constraint, indices(4)

    logical :: converged

    n_design = design%size()
    n_constraint = problem%get_n_constraints()

    ! Grab some local pointers
    call neko_scratch_registry%request(x, indices(1), n_design, .false.)
    call neko_scratch_registry%request(constraint_value, indices(2), &
         n_constraint, .false.)
    call neko_scratch_registry%request(objective_sensitivities, indices(3), &
         n_design, .false.)
    call neko_scratch_registry%request(constraint_sensitivities, indices(4), &
         n_constraint, n_design, .false.)

    !  Retrieve the current objective and constraint values and sensitivities
    call design%get_values(x)
    call problem%get_constraint_values(constraint_value)
    call problem%get_constraint_sensitivities(constraint_sensitivities)

    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(objective_sensitivities)
       ! Convert gradient to directional derivative
       call des%project_sensitivity(objective_sensitivities)
       call des%project_sensitivity(constraint_sensitivities)
    class default
       call problem%get_objective_sensitivities(objective_sensitivities)
    end select

    ! Execute the scaling
    if (this%auto_scale) then
       call constraint_value%copy_from(DEVICE_TO_HOST, sync = .true.)
       this%scaling_factor = abs(this%scale / constraint_value%x(1))
    end if

    if (.not. abscmp(this%scaling_factor, 1.0_rp)) then
       call vector_cmult(constraint_value, this%scaling_factor)
       call matrix_cmult(constraint_sensitivities, this%scaling_factor)
    end if

    ! Update the design variable
    call this%mma%update(iter, x, objective_sensitivities, &
         constraint_value, constraint_sensitivities)
    call design%update_design(x)

    ! Evaluate the problem based on the updated design
    call problem%compute(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_forward(iter)
    end if
    call problem%compute_sensitivity(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_adjoint(iter)
    end if

    ! Retrieve the updated objective and constraint values and sensitivities
    call problem%get_constraint_values(constraint_value)
    call problem%get_constraint_sensitivities(constraint_sensitivities)

    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(objective_sensitivities)
       ! Convert gradient to directional derivative
       call des%project_sensitivity(objective_sensitivities)
       call des%project_sensitivity(constraint_sensitivities)
    class default
       call problem%get_objective_sensitivities(objective_sensitivities)
    end select

    ! Check the KKT conditions and check for convergence
    call this%mma%KKT(x, objective_sensitivities, &
         constraint_value, constraint_sensitivities)
    converged = this%mma%get_residumax() .lt. this%tolerance

    ! Free local resources
    call neko_scratch_registry%relinquish(indices)

  end function mma_optimizer_step

  !> Validate the solution for the MMA optimizer
  subroutine mma_optimizer_validate(this, problem, design)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design

    type(vector_t), pointer :: constraint_values
    integer :: ind

    call neko_scratch_registry%request(constraint_values, ind, &
         problem%get_n_constraints(), .false.)

    call problem%get_constraint_values(constraint_values)
    call constraint_values%copy_from(DEVICE_TO_HOST, sync = .true.)

    if (any(constraint_values%x .gt. 0.0_rp)) then
       call neko_error('MMA optimizer validation failed: ' // &
            'Constraints are not satisfied.')
    end if

    ! Free local resources
    call neko_scratch_registry%relinquish(ind)

  end subroutine mma_optimizer_validate

  ! -------------------------------------------------------------------------- !
  ! Logging and IO methods for the MMA optimizer

  !> Write the progress of the MMA optimizer to the log file
  !! This subroutine logs the current iteration, objective values,
  !! constraint values, and convergence metrics to a CSV file.
  !! @param this The MMA optimizer object.
  !! @param iter The current iteration number.
  !! @param problem The problem object.
  subroutine mma_optimizer_write(this, iter, problem)
    class(mma_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(inout) :: problem
    real(kind=rp) :: extras(3)

    if (.not. this%enable_output) return
    call profiler_start_region('Optimizer logging')

    if (iter .eq. 0) then
       extras(1) = 0.0_rp
       extras(2) = 0.0_rp
    else
       extras(1) = this%mma%get_residumax()
       extras(2) = this%mma%get_residunorm()
    end if
    extras(3) = this%scaling_factor

    call this%write_log(iter, problem, extras)

    call profiler_end_region('Optimizer logging')
  end subroutine mma_optimizer_write

  ! -------------------------------------------------------------------------- !
  ! Checkpointing methods for the MMA optimizer

  !> Save the MMA optimizer-specific checkpoint data
  subroutine mma_optimizer_save_checkpoint_components(this, filename, overwrite)
    class(mma_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite

    call this%mma%save_checkpoint(filename, overwrite)
  end subroutine mma_optimizer_save_checkpoint_components

  !> Restore the MMA optimizer-specific checkpoint data
  subroutine mma_optimizer_load_checkpoint_components(this, filename)
    class(mma_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: filename

    call this%mma%load_checkpoint(filename)
  end subroutine mma_optimizer_load_checkpoint_components

end module mma_optimizer
