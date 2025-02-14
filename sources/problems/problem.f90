!> @file problem.f90
!! @copyright (c) 2024-2025, The Neko-TOP Authors
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

!> Module for handling the optimization problem.
module problem
  use num_types, only: rp
  use fld_file_output, only: fld_file_output_t
  use design, only: design_t
  use objective, only: objective_t, objective_wrapper_t, objective_factory
  use constraint, only: constraint_t, constraint_wrapper_t, constraint_factory
  use vector, only: vector_t
  use matrix, only: matrix_t
  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  use neko_config, only: NEKO_BCKND_DEVICE
  use json_module, only: json_file
  use json_utils, only: json_extract_item, json_get
  use simulation, only: simulation_t
  use logger, only: neko_log
  use device_math, only: device_copy

  implicit none
  private

  !> The abstract problem type.
  !!
  !! This module defines the `problem_t` type which is the main interface for
  !! the optimization problem. The problem is defined by a set of objectives and
  !! constraints that are evaluated based on the design variables. The problem
  !! also handles the output of the problem and the simulation.
  type, abstract, public :: problem_t
     private

     !> The number of design variables.
     integer :: n_design
     !> Number of objectives in the problem.
     integer :: n_objectives
     !> Number of constraints in the problem.
     integer :: n_constraints

     !> The objective of the problem.
     class(objective_wrapper_t), allocatable, dimension(:) :: objective_list
     !> The constraints of the problem.
     class(constraint_wrapper_t), allocatable, dimension(:) :: constraint_list

     !> An output sampler for the problem. This should probably be an output
     !! controller at some point instead.
     type(fld_file_output_t), public :: output

   contains

     ! ----------------------------------------------------------------------- !
     ! Interfaces

     !> Constructor for physics of the problem.
     procedure(problem_init), pass(this), public, deferred :: init
     !> Destructor.
     procedure(problem_free), pass(this), deferred, public :: free

     !> Evaluate the optimization problem.
     !! This is the main function that evaluates the problem. It should be
     !! implemented in the derived classes.
     procedure(problem_compute), pass(this), public, deferred :: compute

     !> Evaluate the sensitivity of the optimization problem.
     !! This is the main function that evaluates the problem sensitivity to the
     !! design. It should be implemented in the derived classes.
     procedure(problem_compute_sensitivity), pass(this), public, deferred :: &
          compute_sensitivity

     ! ----------------------------------------------------------------------- !
     ! Base class methods

     !> Constructor for the base class.
     procedure, pass(this) :: init_base => problem_init_base
     !> Destructor for the base class.
     procedure, pass(this) :: free_base => problem_free_base

     !> Read objective json-file.
     procedure, pass(this), public :: read_objectives => problem_read_objectives
     !> Read constraint json-file.
     procedure, pass(this), public :: read_constraints => &
          problem_read_constraints

     ! ----------------------------------------------------------------------- !
     ! Actual methods

     !> Sample the problem.
     procedure, pass(this), public :: write => problem_write

     !> Add an objective to the list.
     procedure, pass(this), public :: add_objective => problem_add_objective
     !> Add a constraint to the list.
     procedure, pass(this), public :: add_constraint => problem_add_constraint

     ! ----------------------------------------------------------------------- !
     ! Internal Updater methods

     !> Update the objective function.
     procedure, pass(this) :: update_objectives => &
          problem_update_objectives
     !> Update the volume constraint.
     procedure, pass(this) :: update_constraints => &
          problem_update_constraints
     !> Update the objective sensitivities.
     procedure, pass(this) :: update_objective_sensitivities => &
          problem_update_objective_sensitivities
     !> Update the constraint sensitivities.
     procedure, pass(this) :: update_constraint_sensitivities => &
          problem_update_constraint_sensitivities

     ! ----------------------------------------------------------------------- !
     ! Public Getters

     !> Return the objective.
     procedure, pass(this), public :: get_objective_value => &
          problem_get_objective_value
     !> Return all components of the objective.
     procedure, pass(this), public :: get_all_objective_values => &
          problem_get_all_objective_values
     !> Return the constraints.
     procedure, pass(this), public :: get_constraint_values => &
          problem_get_constraint_values
     !> Return the sensitivity of the objective.
     procedure, pass(this), public :: get_objective_sensitivities => &
          problem_get_objective_sensitivities
     !> Return the sensitivity of the constraints.
     procedure, pass(this), public :: get_constraint_sensitivities => &
          problem_get_constraint_sensitivities

     !> Return the number of objectives.
     procedure, pass(this) :: get_n_objectives => problem_get_num_objectives
     !> Return the number of constraints.
     procedure, pass(this) :: get_n_constraints => problem_get_num_constraints

     !> Return the logfile header
     procedure, pass(this) :: get_log_header => problem_get_log_header

  end type problem_t

  ! -------------------------------------------------------------------------- !
  ! Interfaces for the derived types

  abstract interface
     !> Constructor for physics of the problem.
     !! This is the main constructor for a problem. This should be defined in
     !! the derived types to initialize the problem. This is based on the
     !! abstract design type, We suggest that a switch statement is used to
     !! initialize the problem based on the design type.
     subroutine problem_init(this, parameters, design, simulation)
       import problem_t, json_file, simulation_t, design_t
       class(problem_t), intent(inout) :: this
       type(json_file), intent(inout) :: parameters
       class(design_t), intent(in) :: design
       type(simulation_t), intent(inout) :: simulation
     end subroutine problem_init

     !> Compute the problem.
     subroutine problem_compute(this, design)
       import problem_t, design_t
       class(problem_t), intent(inout) :: this
       class(design_t), intent(inout) :: design
     end subroutine problem_compute

     !> Compute the problem.
     subroutine problem_compute_sensitivity(this, design)
       import problem_t, design_t
       class(problem_t), intent(inout) :: this
       class(design_t), intent(inout) :: design
     end subroutine problem_compute_sensitivity

     !> Destructor
     subroutine problem_free(this)
       import problem_t
       class(problem_t), intent(inout) :: this
     end subroutine problem_free
  end interface

contains

  ! ========================================================================== !
  ! Base class methods

  !> Constructor for the base class
  subroutine problem_init_base(this, n_design)
    class(problem_t), intent(inout) :: this
    integer, intent(in) :: n_design

    this%n_design = n_design
    this%n_objectives = 0
    this%n_constraints = 0

  end subroutine problem_init_base

  !> Destructor for the base class
  subroutine problem_free_base(this)
    class(problem_t), intent(inout) :: this
    integer :: i

    ! Free the objective list
    if (allocated(this%objective_list)) then
       do i = 1, size(this%objective_list)
          call this%objective_list(i)%free()
       end do
       deallocate(this%objective_list)
    end if

    ! Free the constraint list
    if (allocated(this%constraint_list)) then
       do i = 1, size(this%constraint_list)
          call this%constraint_list(i)%free()
       end do
       deallocate(this%constraint_list)
    end if
  end subroutine problem_free_base

  !> Sample the fields/design.
  subroutine problem_write(this, idx)
    class(problem_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output%sample(real(idx, kind=rp))
  end subroutine problem_write

  ! ========================================================================== !
  ! Handling constraints and objectives

  !> Read the objective from a parameters file.
  subroutine problem_read_objectives(this, parameters, simulation, design)
    class(problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation
    class(design_t), intent(in) :: design
    class(objective_t), allocatable :: objective

    ! A single objective term as its own json_file.
    character(len=:), allocatable :: path, type
    type(json_file) :: objective_json
    integer :: n_objectives, i

    call neko_log%section("Reading objectives")

    ! Get the number of objectives.
    path = "optimization.objectives"
    call parameters%info(path, n_children = n_objectives)

    ! Grab a single parameters entry and create a constraint from it.
    do i = 1, n_objectives
       call json_extract_item(parameters, path, i, objective_json)
       call json_get(objective_json, "type", type)
       call neko_log%message(type)

       call objective_factory(objective, objective_json, design, simulation)
       call this%add_objective(objective)

    end do

    call neko_log%end_section()

  end subroutine problem_read_objectives

  !> Read the constraint from a parameters file.
  subroutine problem_read_constraints(this, parameters, simulation, design)
    class(problem_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(inout) :: simulation
    class(design_t), intent(in) :: design
    class(constraint_t), allocatable :: constraint

    ! A single constraint term as its own json_file.
    character(len=:), allocatable :: path, type
    type(json_file) :: constraint_json
    integer :: n_constraints, i

    call neko_log%section("Reading constraints")

    ! Get the number of constraints.
    path = "optimization.constraints"
    call parameters%info(path, n_children = n_constraints)

    ! Grab a single parameters entry and create a constraint from it.
    do i = 1, n_constraints
       call json_extract_item(parameters, path, i, constraint_json)
       call json_get(constraint_json, "type", type)
       call neko_log%message(type)

       call constraint_factory(constraint, constraint_json, design, simulation)
       call this%add_constraint(constraint)

    end do

    call neko_log%end_section()

  end subroutine problem_read_constraints

  !> Add an objective to the list.
  subroutine problem_add_objective(this, objective)
    class(problem_t), intent(inout) :: this
    class(objective_t), allocatable, intent(inout) :: objective
    class(objective_wrapper_t), allocatable, dimension(:) :: temp_list
    integer :: i, n

    n = 0
    if (allocated(this%objective_list)) then
       n = size(this%objective_list)
       call move_alloc(this%objective_list, temp_list)
       allocate(this%objective_list(n + 1))
       if (allocated(temp_list)) then
          do i = 1, n
             call move_alloc(temp_list(i)%objective, &
                  this%objective_list(i)%objective)
          end do
       end if
    else
       allocate(this%objective_list(1))
    end if

    call move_alloc(objective, this%objective_list(n + 1)%objective)
    this%n_objectives = n + 1
  end subroutine problem_add_objective

  !> Add an objective to the list.
  subroutine problem_add_constraint(this, constraint)
    class(problem_t), intent(inout) :: this
    class(constraint_t), allocatable, intent(inout) :: constraint
    class(constraint_wrapper_t), allocatable, dimension(:) :: temp_list
    integer :: i, n

    n = 0
    if (allocated(this%constraint_list)) then
       n = size(this%constraint_list)
       call move_alloc(this%constraint_list, temp_list)
       allocate(this%constraint_list(n + 1))
       if (allocated(temp_list)) then
          do i = 1, n
             call move_alloc(temp_list(i)%constraint, &
                  this%constraint_list(i)%constraint)
          end do
       end if
    else
       allocate(this%constraint_list(1))
    end if

    call move_alloc(constraint, this%constraint_list(n + 1)%constraint)
    this%n_constraints = n + 1
  end subroutine problem_add_constraint

  ! ========================================================================== !
  ! Update the objectives and constraints

  !> Update the objectives.
  !!
  !! This function should be called after the design has been updated.
  !! It will update the value of all the objectives.
  !! @param[in] design The design to update the objectives with.
  subroutine problem_update_objectives(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    integer :: i

    do i = 1, this%n_objectives
       call this%objective_list(i)%objective%update_value(design)
    end do
  end subroutine problem_update_objectives

  !> Update the constraints.
  !!
  !! This function should be called after the design has been updated.
  !! It will update the value of all the constraints.
  !! @param[in] design The design to update the constraints with.
  subroutine problem_update_constraints(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    integer :: i

    do i = 1, this%n_constraints
       call this%constraint_list(i)%constraint%update_value(design)
    end do
  end subroutine problem_update_constraints

  !> Update the sensitivity of the objectives.
  !!
  !! This function should be called after the design has been updated.
  !! It will update the sensitivity of all the objectives.
  !! @param[in] design The design to update the objectives with.
  subroutine problem_update_objective_sensitivities(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    integer :: i

    do i = 1, this%n_objectives
       call this%objective_list(i)%objective%update_sensitivity(design)
    end do
  end subroutine problem_update_objective_sensitivities

  !> Update the sensitivity of the constraints.
  !!
  !! This function should be called after the design has been updated.
  !! It will update the sensitivity of all the constraints.
  !! @param[in] design The design to update the constraints with.
  subroutine problem_update_constraint_sensitivities(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    integer :: i

    do i = 1, this%n_constraints
       call this%constraint_list(i)%constraint%update_sensitivity(design)
    end do
  end subroutine problem_update_constraint_sensitivities

  ! ========================================================================== !
  ! Problem part getters

  !> Construct and get the objective.
  !!
  !! This function constructs the objective value from the individual
  !! objectives and their weights.
  !! @param[out] objective_value The weighted sum of all objective values.
  subroutine problem_get_objective_value(this, objective_value)
    class(problem_t), intent(inout) :: this
    real(kind=rp), intent(out) :: objective_value
    integer :: i

    objective_value = 0.0_rp
    do i = 1, this%n_objectives
       objective_value = objective_value + &
            this%objective_list(i)%objective%weight * &
            this%objective_list(i)%objective%value
    end do

  end subroutine problem_get_objective_value

  !> Construct and get the objective.
  !!
  !! This function returns all the indivual objectives comprising the
  !! objective function
  !! @param[out] all_objective_values A vector containing all objectives
  subroutine problem_get_all_objective_values(this, all_objective_values)
    class(problem_t), intent(inout) :: this
    type(vector_t), intent(out) :: all_objective_values
    integer :: i

    call all_objective_values%init(this%n_objectives)

    do i = 1, this%n_objectives
       all_objective_values%x(i) = this%objective_list(i)%objective%value
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(all_objective_values%x, all_objective_values%x_d, &
            this%n_objectives, HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine problem_get_all_objective_values

  !> Construct and get the constraints.
  !!
  !! This function constructs the constraint values from the individual
  !! constraints.
  !! @param[out] constraint_value The vector of all constraint values.
  subroutine problem_get_constraint_values(this, constraint_value)
    class(problem_t), intent(inout) :: this
    type(vector_t), intent(out) :: constraint_value
    integer :: i

    call constraint_value%init(this%n_constraints)

    do i = 1, this%n_constraints
       constraint_value%x(i) = this%constraint_list(i)%constraint%value
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(constraint_value%x, constraint_value%x_d, &
            this%n_constraints, HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine problem_get_constraint_values

  !> Construct and get the sensitivity of the objective.
  !!
  !! This function constructs the sensitivity of the objective value from the
  !! individual objectives and their weights.
  !! @param[out] sensitivity The weighted sum of all objective sensitivities.
  subroutine problem_get_objective_sensitivities(this, sensitivity)
    class(problem_t), intent(inout) :: this
    type(vector_t), intent(out) :: sensitivity
    integer :: i

    call sensitivity%init(this%n_design)

    do i = 1, this%n_objectives
       sensitivity = sensitivity + this%objective_list(i)%objective%sensitivity
    end do

  end subroutine problem_get_objective_sensitivities

  !> Construct and get the sensitivity of the constraints.
  !!
  !! This function constructs the sensitivity of the constraint values from the
  !! individual constraints.
  !! @param[out] sensitivity The matrix of all constraint sensitivities.
  subroutine problem_get_constraint_sensitivities(this, sensitivity)
    class(problem_t), intent(inout) :: this
    type(matrix_t), intent(out) :: sensitivity
    integer :: i, j

    call sensitivity%init(this%n_design, this%n_constraints)

    do i = 1, this%n_constraints
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy( &
               this%constraint_list(i)%constraint%sensitivity%x, &
               this%constraint_list(i)%constraint%sensitivity%x_d, &
               this%n_design, DEVICE_TO_HOST, sync = .true.)
       end if
       do j = 1, this%n_design
          sensitivity%x(j, i) = &
               this%constraint_list(i)%constraint%sensitivity%x(j)
       end do
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(sensitivity%x, sensitivity%x_d, &
            this%n_design * this%n_constraints, HOST_TO_DEVICE, sync = .true.)
    end if

  end subroutine problem_get_constraint_sensitivities

  ! ========================================================================== !
  ! Simple getters

  !> Return the number of objectives.
  pure function problem_get_num_objectives(this) result(n)
    class(problem_t), intent(in) :: this
    integer :: n

    n = this%n_objectives
  end function problem_get_num_objectives

  !> Return the number of constraints.
  pure function problem_get_num_constraints(this) result(n)
    class(problem_t), intent(in) :: this
    integer :: n

    n = this%n_constraints
  end function problem_get_num_constraints

  !> Return the header for the problem.
  function problem_get_log_header(this) result(buff)
    class(problem_t), intent(in) :: this
    character(len=1024) :: buff
    character(len=50) :: mini_buff
    integer :: i

    ! When it comes to multi-objective optimization
    ! (handled in the way that we do) we want to know the value of each
    ! objective individually, not just the combined effect.
    !
    ! my vision is:
    !
    !      | Total F | F_1 | F_2 | ... | F_n | C_1 | C_2 | ... | C_m |
    !
    ! And then if we also want things like thie iteration or KKT they can be
    ! appended to the begining or end of this by the optimizer.
    !
    ! iter | Total F | F_1 | F_2 | ... | F_n | C_1 | C_2 | ... | C_m | KKT
    buff = "Total objective function"
    do i = 1, this%get_n_objectives()
       mini_buff = ""
       write(mini_buff, '(", ", A)') this%objective_list(i)%objective%name
       buff = trim(buff)//trim(mini_buff)
    end do

    do i = 1, this%get_n_constraints()
       mini_buff = ""
       write(mini_buff, '(", ", A)') &
            this%constraint_list(i)%constraint%name
       buff = trim(buff)//trim(mini_buff)
    end do

  end function problem_get_log_header
end module problem
