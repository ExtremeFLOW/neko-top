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
  use neko_config, only: NEKO_BCKND_DEVICE
  use constraint, only: constraint_t
  use dummy_constraint, only: dummy_constraint_t

  ! External modules
  use json_module, only: json_file
  use vector, only: vector_t
  use matrix, only: matrix_t
  use math, only: abscmp, copy, glsum
  use profiler, only: profiler_start_region, profiler_end_region
  use logger, only: neko_log
  use csv_file, only: csv_file_t
  use vector_math, only: vector_cmult, vector_col2, vector_invcol2
  use matrix_math, only: matrix_cmult
  use device, only: DEVICE_TO_HOST, HOST_TO_DEVICE
  use scratch_registry, only: neko_scratch_registry

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
     type(csv_file_t), private :: csv_log
     !> MMA variable transform weights (z = w * x).
     type(vector_t), private :: mma_weights
     !> Inverse squared weights (1 / w^2), used for sensitivity transform.
     type(vector_t), private :: mma_inv_weight_sq
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
     procedure, pass(this), private :: init_variable_weights => &
          mma_optimizer_init_variable_weights
     procedure, pass(this), private :: transform_variables => &
          mma_optimizer_transform_variables
     procedure, pass(this), private :: transform_sensitivities => &
          mma_optimizer_transform_sensitivities

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
    character(len=1024) :: header
    class(constraint_t), allocatable :: dummy_con

    call neko_log%section('Optimizer Initialization')

    ! Check if the problem is unconstrained
    this%unconstrained_problem = problem%get_n_constraints() .eq. 0
    if (this%unconstrained_problem) then
       call neko_log%message('Unconstrained problem detected. ' // &
            'Adding a dummy constraint to enable MMA optimization.')

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

    call this%init_variable_weights(design, simulation)

    call design%get_values(x)
    call this%transform_variables(x, .true.)
    call this%mma%init(x, design%size(), problem%get_n_constraints(), &
         solver_parameters, this%scale, this%auto_scale)
    call this%mma%scale_variable_bounds(this%mma_weights)

    call neko_scratch_registry%relinquish(ind)

    !set the enable_output flag
    this%enable_output = enable_output
    this%scaling_factor = this%scale
    this%tolerance = tolerance

    ! Initialize the logger
    if (this%enable_output) then
       call this%csv_log%init('optimization_data.csv')
       header = 'iter, ' // trim(problem%get_log_header()) // &
            ', KKTmax, KKTnorm2, scaling factor, ' // &
            trim(this%mma%get_backend_and_subsolver())

       call this%csv_log%set_header(trim(header))
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
    call this%mma_weights%free()
    call this%mma_inv_weight_sq%free()
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
    call problem%compute_sensitivity(design, simulation)

    ! Retrieve the updated objective and constraint values and sensitivities
    call design%get_values(x)
    call this%transform_variables(x, .true.)
    call problem%get_constraint_values(constraint_value)

    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(objective_sensitivities)
    class default
       call problem%get_objective_sensitivities(objective_sensitivities)
    end select

    call problem%get_constraint_sensitivities(constraint_sensitivities)
    call this%transform_sensitivities(objective_sensitivities, &
         constraint_sensitivities)

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
    call this%transform_variables(x, .true.)
    call problem%get_constraint_values(constraint_value)

    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(objective_sensitivities)
    class default
       call problem%get_objective_sensitivities(objective_sensitivities)
    end select

    call problem%get_constraint_sensitivities(constraint_sensitivities)
    call this%transform_sensitivities(objective_sensitivities, &
         constraint_sensitivities)

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
    call this%transform_variables(x, .false.)
    call design%update_design(x)
    call this%transform_variables(x, .true.)

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

    select type (des => design)
    type is (brinkman_design_t)
       call des%get_sensitivity(objective_sensitivities)
    class default
       call problem%get_objective_sensitivities(objective_sensitivities)
    end select

    call problem%get_constraint_sensitivities(constraint_sensitivities)
    call this%transform_sensitivities(objective_sensitivities, &
         constraint_sensitivities)

    ! Check the KKT conditions and check for convergence
    call this%mma%KKT(x, objective_sensitivities, &
         constraint_value, constraint_sensitivities)
    converged = this%mma%get_residumax() .lt. this%tolerance

    ! Free local resources
    call neko_scratch_registry%relinquish(indices)

  end function mma_optimizer_step

  !> Initialize MMA variable transform weights and inverse-squared weights.
  !! Uses identity weights by default and, when available, sets
  !! `weights = sqrt(B) / avg(sqrt(B))` from the simulation mass-matrix
  !! coefficients.
  !! @param[inout] this The MMA optimizer.
  !! @param[in] design The design object used to determine vector size.
  !! @param[in] simulation Optional simulation with coefficient data.
  subroutine mma_optimizer_init_variable_weights(this, design, simulation)
    class(mma_optimizer_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(simulation_t), optional, intent(in) :: simulation
    integer :: n, i, n_coef
    real(kind=rp) :: local_sum, global_sum, weight_avg, global_count
    real(kind=rp) :: local_count(1)

    n = design%size()
    call this%mma_weights%init(n)
    call this%mma_inv_weight_sq%init(n)
    this%mma_weights = 1.0_rp
    this%mma_inv_weight_sq = 1.0_rp

    if (present(simulation)) then
      n_coef = simulation%fluid%c_Xh%dof%size()
      if (n_coef .eq. n) then
         call copy(this%mma_weights%x, simulation%fluid%c_Xh%B, n)
         do i = 1, n
            if (this%mma_weights%x(i) .le. 0.0_rp) then
               call neko_error('mma_optimizer: non-positive mass-matrix entry ' // &
                    'encountered when building variable weights')
            end if
            this%mma_weights%x(i) = sqrt(this%mma_weights%x(i))
         end do

         global_sum = glsum(this%mma_weights%x, n)
         local_count(1) = real(n, rp)
         global_count = glsum(local_count, 1)

         if (global_count .le. 0.0_rp) then
            call neko_error('mma_optimizer: invalid global design size when ' // &
                 'normalizing variable weights')
         end if

         weight_avg = global_sum / global_count
         if (weight_avg .le. 0.0_rp) then
            call neko_error('mma_optimizer: non-positive average weight ' // &
                 'encountered when normalizing variable weights')
         end if

         !------- DELETE -------------
         ! weight_avg = 1.0_rp
         !----------------------------
         this%mma_weights%x = this%mma_weights%x / weight_avg
         this%mma_inv_weight_sq%x = 1.0_rp / &
              (this%mma_weights%x * this%mma_weights%x * this%mma_weights%x)
         !------ DELETE -------------
         ! this%mma_weights%x = 1.0_rp
         !---------------------------
      else
         call neko_log%message('mma_optimizer: design size and coefficient ' // &
              'size differ; using identity variable weights.')
      end if
    else
      call neko_log%message('mma_optimizer: simulation not present; using ' // &
           'identity variable weights.')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call this%mma_weights%copy_from(HOST_TO_DEVICE, .true.)
       call this%mma_inv_weight_sq%copy_from(HOST_TO_DEVICE, .true.)
    end if
  end subroutine mma_optimizer_init_variable_weights

  !> Transform design variables between internal and MMA spaces.
  !! Applies `x = w*x` when `multiply=.true.` and `x = x/w` otherwise.
  !! @param[inout] this The MMA optimizer.
  !! @param[inout] x The variable vector to transform.
  !! @param[in] multiply Whether to multiply (`.true.`) or divide (`.false.`).
  subroutine mma_optimizer_transform_variables(this, x, multiply)
    class(mma_optimizer_t), intent(inout) :: this
    type(vector_t), intent(inout) :: x
    logical, intent(in) :: multiply

    if (multiply) then
       call vector_col2(x, this%mma_weights)
    else
       call vector_invcol2(x, this%mma_weights)
    end if
  end subroutine mma_optimizer_transform_variables

  !> Transform objective and constraint sensitivities for MMA-space variables.
  !! Applies elementwise weighting with `1/w^2` (equivalent to `B^-1` for
  !! `w=sqrt(B)`).
  !! @param[inout] this The MMA optimizer.
  !! @param[inout] objective_sensitivities Objective sensitivity vector.
  !! @param[inout] constraint_sensitivities Constraint sensitivity matrix.
  subroutine mma_optimizer_transform_sensitivities(this, objective_sensitivities, &
       constraint_sensitivities)
    class(mma_optimizer_t), intent(inout) :: this
    type(vector_t), intent(inout) :: objective_sensitivities
    type(matrix_t), intent(inout) :: constraint_sensitivities
    integer :: j, n

    call vector_col2(objective_sensitivities, this%mma_inv_weight_sq)

    n = constraint_sensitivities%get_ncols()
    if (n .ne. this%mma_inv_weight_sq%size()) then
       call neko_error('mma_optimizer: constraint sensitivity width does not ' // &
            'match weight vector size')
    end if

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call constraint_sensitivities%copy_from(DEVICE_TO_HOST, .true.)
    end if

    do j = 1, n
       constraint_sensitivities%x(:, j) = constraint_sensitivities%x(:, j) * &
            this%mma_inv_weight_sq%x(j)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call constraint_sensitivities%copy_from(HOST_TO_DEVICE, .true.)
    end if
  end subroutine mma_optimizer_transform_sensitivities

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
       !call neko_error('MMA optimizer validation failed: ' // &
       !     'Constraints are not satisfied.')
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
    class(problem_t), intent(in) :: problem

    type(vector_t), pointer :: log_data
    type(vector_t), pointer :: all_objectives
    type(vector_t), pointer :: constraint_value
    real(kind=rp) :: objective_value

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

    ! Prepare data for logging
    call problem%get_objective_value(objective_value)
    call problem%get_all_objective_values(all_objectives)
    call problem%get_constraint_values(constraint_value)

    call all_objectives%copy_from(DEVICE_TO_HOST, sync = .true.)
    call constraint_value%copy_from(DEVICE_TO_HOST, sync = .true.)

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
    call neko_scratch_registry%relinquish(ind)

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
