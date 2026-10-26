!> @file steepest_descent_optimizer.f90
!! @copyright
!! Copyright (c) 2026, The Neko-TOP Authors
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

!> Steepest descent optimizer with move limits.
!! @details
!! Implements the steepest descent method with move limits (SDM) of
!! Appendix D in da Silva, Aage, Beck & Sigmund, "Large scale
!! three-dimensional manufacturing tolerant stress-constrained topology
!! optimization", IJNME (2021), Eqs. D22-D24. Each iteration
!! 1. takes the steepest-descent direction S = -G and zeroes the components
!!    that point out of the box at an active bound (D22),
!! 2. normalises it by its global maximum absolute value,
!! 3. adapts the per-variable move limits from the sign of the last two
!!    design changes, from the third iteration on (D23),
!! 4. takes a move-limited, box-constrained step of length `step_length`
!!    (D24).
!!
!! The method handles box bounds only; problems with constraints are
!! rejected in favour of the "mma" optimizer. There is no tolerance-based
!! convergence criterion: the step reports convergence only when the
!! bound-projected gradient vanishes, and otherwise the run is stopped by the
!! base criteria (`max_iterations`, `stop_design_change`, `max_runtime`).
module steepest_descent_optimizer
  use optimizer, only: optimizer_t
  use problem, only: problem_t
  use design, only: design_t
  use brinkman_design, only: brinkman_design_t
  use simulation_m, only: simulation_t
  use num_types, only: rp
  use utils, only: neko_error
  use json_utils, only: json_get, json_get_or_default
  use logger, only: neko_log
  use profiler, only: profiler_start_region, profiler_end_region
  use scratch_registry, only: neko_scratch_registry
  use steepest_descent_math, only: sd_direction, sd_update_move_limits, &
       sd_update_design

  ! External modules
  use json_module, only: json_file
  use vector, only: vector_t
  use vector_math, only: vector_copy, vector_cfill, vector_sub3, &
       vector_absval, vector_glmax, vector_glmin, vector_glsc2, &
       vector_glsubnorm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite

  implicit none
  private

  public :: steepest_descent_optimizer_t

  !> Steepest descent optimizer with move limits.
  type, extends(optimizer_t) :: steepest_descent_optimizer_t

     !> Design variables one iteration back.
     type(vector_t), private :: xold1
     !> Design variables two iterations back.
     type(vector_t), private :: xold2
     !> Per-variable move limits.
     type(vector_t), private :: delta

     !> Lower bound of the design variables.
     real(kind=rp), private :: xmin = 0.0_rp
     !> Upper bound of the design variables.
     real(kind=rp), private :: xmax = 1.0_rp
     !> Step length applied to the normalised direction.
     real(kind=rp), private :: step_length = 1.0_rp
     !> Initial move limit.
     real(kind=rp), private :: move_limit_initial = 0.1_rp
     !> Smallest allowed move limit.
     real(kind=rp), private :: move_limit_min = 0.001_rp
     !> Largest allowed move limit.
     real(kind=rp), private :: move_limit_max = 0.1_rp
     !> Factor applied to a move limit when its variable oscillates.
     real(kind=rp), private :: move_limit_decrease = 0.7_rp
     !> Factor applied to a move limit when its variable moves monotonically.
     real(kind=rp), private :: move_limit_increase = 1.1_rp
     !> Global maximum absolute value of the last projected gradient.
     real(kind=rp), private :: last_max_gradient = 0.0_rp

     !> Write the convergence history and the simulation output.
     logical, private :: enable_output = .true.

   contains

     ! Override the deferred methods
     procedure, pass(this) :: init_from_json => sd_optimizer_init_from_json
     procedure, pass(this) :: free => sd_optimizer_free
     procedure, pass(this) :: initialize => sd_optimizer_initialize
     procedure, pass(this) :: step => sd_optimizer_step
     procedure, pass(this) :: validate => sd_optimizer_validate
     procedure, pass(this) :: write => sd_optimizer_write

     procedure, pass(this) :: save_checkpoint_components => &
          sd_optimizer_save_checkpoint_components
     procedure, pass(this) :: load_checkpoint_components => &
          sd_optimizer_load_checkpoint_components

  end type steepest_descent_optimizer_t

  ! -------------------------------------------------------------------------- !
  ! IO routines for HDF5 checkpoints

  interface
     ! Save the steepest descent state to an HDF5 checkpoint
     module subroutine sd_optimizer_save_checkpoint_hdf5(object, filename, &
          overwrite)
       class(steepest_descent_optimizer_t), intent(inout) :: object
       character(len=*), intent(in) :: filename
       logical, intent(in), optional :: overwrite
     end subroutine sd_optimizer_save_checkpoint_hdf5

     ! Load the steepest descent state from an HDF5 checkpoint
     module subroutine sd_optimizer_load_checkpoint_hdf5(object, filename)
       class(steepest_descent_optimizer_t), intent(inout) :: object
       character(len=*), intent(in) :: filename
     end subroutine sd_optimizer_load_checkpoint_hdf5
  end interface

contains

  ! -------------------------------------------------------------------------- !
  ! Constructor and destructor

  !> Initialise the steepest descent optimizer from the case file.
  !! @details Reads `optimization.solver` and its `steepest_descent` object.
  !! The base settings are read after `free`, which resets them, and before
  !! `init_base`, which keeps them.
  !! @param this The steepest descent optimizer.
  !! @param parameters The JSON case file.
  !! @param problem The problem to optimise, which must have no constraints.
  !! @param design The design whose values start the optimisation.
  !! @param simulation The simulation, unused.
  subroutine sd_optimizer_init_from_json(this, parameters, problem, design, &
       simulation)
    class(steepest_descent_optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), optional, intent(in) :: simulation

    type(json_file) :: solver_parameters
    integer :: max_iterations, n

    call this%free()

    call neko_log%section('Optimizer Initialization')

    call json_get(parameters, 'optimization.solver', solver_parameters)
    call this%read_base_settings(solver_parameters)

    if (problem%get_n_constraints() .gt. 0) then
       call neko_error('steepest_descent: the problem has constraints, ' // &
            'which this optimizer does not support; use the "mma" ' // &
            'optimizer instead')
    end if

    ! Read and validate the optimizer settings
    call json_get_or_default(solver_parameters, 'max_iterations', &
         max_iterations, 100)
    call json_get_or_default(solver_parameters, 'enable_output', &
         this%enable_output, .true.)

    call json_get_or_default(solver_parameters, 'steepest_descent.xmin', &
         this%xmin, 0.0_rp)
    call json_get_or_default(solver_parameters, 'steepest_descent.xmax', &
         this%xmax, 1.0_rp)
    call json_get_or_default(solver_parameters, &
         'steepest_descent.step_length', this%step_length, 1.0_rp)
    call json_get_or_default(solver_parameters, &
         'steepest_descent.move_limit.initial', this%move_limit_initial, &
         0.1_rp)
    call json_get_or_default(solver_parameters, &
         'steepest_descent.move_limit.min', this%move_limit_min, 0.001_rp)
    call json_get_or_default(solver_parameters, &
         'steepest_descent.move_limit.max', this%move_limit_max, 0.1_rp)
    call json_get_or_default(solver_parameters, &
         'steepest_descent.move_limit.decrease', this%move_limit_decrease, &
         0.7_rp)
    call json_get_or_default(solver_parameters, &
         'steepest_descent.move_limit.increase', this%move_limit_increase, &
         1.1_rp)

    if (this%xmin .ge. this%xmax) then
       call neko_error('steepest_descent: "steepest_descent.xmin" must ' // &
            'be smaller than "steepest_descent.xmax"')
    end if
    if (this%step_length .le. 0.0_rp) then
       call neko_error('steepest_descent: ' // &
            '"steepest_descent.step_length" must be positive')
    end if
    if (this%move_limit_min .le. 0.0_rp) then
       call neko_error('steepest_descent: ' // &
            '"steepest_descent.move_limit.min" must be positive')
    end if
    if (this%move_limit_initial .lt. this%move_limit_min) then
       call neko_error('steepest_descent: ' // &
            '"steepest_descent.move_limit.initial" must not be smaller ' // &
            'than "steepest_descent.move_limit.min"')
    end if
    if (this%move_limit_max .lt. this%move_limit_initial) then
       call neko_error('steepest_descent: ' // &
            '"steepest_descent.move_limit.max" must not be smaller ' // &
            'than "steepest_descent.move_limit.initial"')
    end if
    if (this%move_limit_decrease .le. 0.0_rp .or. &
         this%move_limit_decrease .gt. 1.0_rp) then
       call neko_error('steepest_descent: ' // &
            '"steepest_descent.move_limit.decrease" must lie in (0, 1]')
    end if
    if (this%move_limit_increase .lt. 1.0_rp) then
       call neko_error('steepest_descent: ' // &
            '"steepest_descent.move_limit.increase" must be at least 1')
    end if

    ! Start the design history from the initial design
    n = design%size()
    call this%xold1%init(n)
    call this%xold2%init(n)
    call this%delta%init(n)
    call design%get_values(this%xold1)
    call design%get_values(this%xold2)
    call vector_cfill(this%delta, this%move_limit_initial)

    if (this%enable_output) then
       call this%init_log(problem, extra_headers = ['max_abs_gradient'], &
            include_constraints = .false.)
    end if

    call this%init_base('steepest_descent', max_iterations)

    call solver_parameters%destroy()
    call neko_log%end_section()

  end subroutine sd_optimizer_init_from_json

  !> Free the steepest descent optimizer.
  !! @param this The steepest descent optimizer.
  subroutine sd_optimizer_free(this)
    class(steepest_descent_optimizer_t), intent(inout) :: this

    call this%xold1%free()
    call this%xold2%free()
    call this%delta%free()
    this%last_max_gradient = 0.0_rp

    call this%free_base()

  end subroutine sd_optimizer_free

  ! -------------------------------------------------------------------------- !
  ! Implementation of the deferred methods

  !> Evaluate the problem and its sensitivity at the starting design.
  !! @param this The steepest descent optimizer.
  !! @param problem The problem to optimise.
  !! @param design The design.
  !! @param simulation The simulation, if any.
  subroutine sd_optimizer_initialize(this, problem, design, simulation)
    class(steepest_descent_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation

    call problem%compute(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_forward(0)
    end if
    call problem%compute_sensitivity(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_adjoint(0)
    end if

  end subroutine sd_optimizer_initialize

  !> Take one steepest descent step and re-evaluate the problem.
  !! @details The sensitivity of the current design must be up to date on
  !! entry; it is on exit, for the updated design.
  !! @param this The steepest descent optimizer.
  !! @param iter The current iteration, starting from 1.
  !! @param problem The problem to optimise.
  !! @param design The design, updated in place.
  !! @param simulation The simulation, if any.
  !! @return `.true.` only when the bound-projected gradient vanishes.
  function sd_optimizer_step(this, iter, problem, design, simulation) &
       result(converged)
    class(steepest_descent_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation
    logical :: converged

    type(vector_t), pointer :: x, x_new, s
    integer :: n, indices(3)
    real(kind=rp) :: smax

    converged = .false.
    n = design%size()

    call neko_scratch_registry%request(x, indices(1), n, .false.)
    call neko_scratch_registry%request(x_new, indices(2), n, .false.)
    call neko_scratch_registry%request(s, indices(3), n, .false.)

    ! Retrieve the current design and the gradient of the objective
    call design%get_values(x)
    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(s)
       ! Convert gradient to directional derivative
       call des%project_sensitivity(s)
    class default
       call problem%get_objective_sensitivities(s)
    end select

    ! Bound-projected steepest-descent direction
    call sd_direction(s, x, this%xmin, this%xmax)

    ! The global max and min drop NaN entries, while the sum of squares
    ! propagates them
    if (ieee_is_nan(vector_glsc2(s, s))) then
       call neko_scratch_registry%relinquish(indices)
       call neko_error('steepest_descent: the projected gradient is NaN')
    end if

    ! Maximum norm of the direction
    smax = max(vector_glmax(s), -vector_glmin(s))
    if (.not. ieee_is_finite(smax)) then
       call neko_scratch_registry%relinquish(indices)
       call neko_error('steepest_descent: the projected gradient is ' // &
            'not finite')
    end if

    ! A vanishing projected gradient is a KKT point of the box-constrained
    ! problem, so there is no step to take
    if (smax .le. tiny(1.0_rp)) then
       this%max_design_change = 0.0_rp
       this%norm2_design_change = 0.0_rp
       this%last_max_gradient = 0.0_rp
       converged = .true.
       nullify(x, x_new, s)
       call neko_scratch_registry%relinquish(indices)
       return
    end if
    this%last_max_gradient = smax

    ! Adapt the move limits once two previous designs exist
    if (iter .ge. 3) then
       call sd_update_move_limits(this%delta, x, this%xold1, this%xold2, &
            this%move_limit_decrease, this%move_limit_increase, &
            this%move_limit_min, this%move_limit_max)
    end if

    ! Update the design variables, which also applies any design mask
    call sd_update_design(x_new, x, s, this%delta, this%step_length, smax, &
         this%xmin, this%xmax)
    call design%update_design(x_new)

    ! Evaluate the problem based on the updated design
    call problem%compute(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_forward(iter)
    end if
    call problem%compute_sensitivity(design, simulation)
    if (present(simulation) .and. this%enable_output) then
       call simulation%write_adjoint(iter)
    end if

    ! Compute L2 norm of design change
    this%norm2_design_change = vector_glsubnorm(x_new, x)

    ! Compute maximum absolute change, reusing s as workspace
    call vector_sub3(s, x_new, x)
    call vector_absval(s)
    this%max_design_change = vector_glmax(s)

    ! Shift the design history, now that the move limits have used it
    call vector_copy(this%xold2, this%xold1)
    call vector_copy(this%xold1, x)

    ! Free local resources
    nullify(x, x_new, s)
    call neko_scratch_registry%relinquish(indices)

  end function sd_optimizer_step

  !> Validate the final design.
  !! @details Nothing to check: the update keeps the design within its box
  !! bounds by construction, and problems with constraints are rejected at
  !! initialisation.
  !! @param this The steepest descent optimizer.
  !! @param problem The problem.
  !! @param design The design.
  subroutine sd_optimizer_validate(this, problem, design)
    class(steepest_descent_optimizer_t), intent(inout) :: this
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design

  end subroutine sd_optimizer_validate

  ! -------------------------------------------------------------------------- !
  ! Logging and IO methods

  !> Write the progress of the optimizer to the log file.
  !! @param this The steepest descent optimizer.
  !! @param iter The current iteration.
  !! @param problem The problem.
  subroutine sd_optimizer_write(this, iter, problem)
    class(steepest_descent_optimizer_t), intent(inout) :: this
    integer, intent(in) :: iter
    class(problem_t), intent(inout) :: problem

    if (.not. this%enable_output) return
    call profiler_start_region('Optimizer logging')

    call this%write_log(iter, problem, [this%last_max_gradient])

    call profiler_end_region('Optimizer logging')
  end subroutine sd_optimizer_write

  !> Save the steepest descent state to a checkpoint.
  !! @param this The steepest descent optimizer.
  !! @param filename The checkpoint file.
  !! @param overwrite Replace an existing steepest descent group.
  subroutine sd_optimizer_save_checkpoint_components(this, filename, &
       overwrite)
    class(steepest_descent_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite

    call sd_optimizer_save_checkpoint_hdf5(this, filename, overwrite)
  end subroutine sd_optimizer_save_checkpoint_components

  !> Restore the steepest descent state from a checkpoint.
  !! @param this The steepest descent optimizer.
  !! @param filename The checkpoint file.
  subroutine sd_optimizer_load_checkpoint_components(this, filename)
    class(steepest_descent_optimizer_t), intent(inout) :: this
    character(len=*), intent(in) :: filename

    call sd_optimizer_load_checkpoint_hdf5(this, filename)
  end subroutine sd_optimizer_load_checkpoint_components

  ! ========================================================================== !
  ! Dummy implementations for module procedures

#if !HAVE_HDF5
  module subroutine sd_optimizer_save_checkpoint_hdf5(object, filename, &
       overwrite)
    class(steepest_descent_optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: overwrite
    call neko_error('steepest_descent: HDF5 support not enabled rebuild ' // &
         'with HAVE_HDF5')
  end subroutine sd_optimizer_save_checkpoint_hdf5

  module subroutine sd_optimizer_load_checkpoint_hdf5(object, filename)
    class(steepest_descent_optimizer_t), intent(inout) :: object
    character(len=*), intent(in) :: filename
    call neko_error('steepest_descent: HDF5 support not enabled rebuild ' // &
         'with HAVE_HDF5')
  end subroutine sd_optimizer_load_checkpoint_hdf5
#endif

end module steepest_descent_optimizer
