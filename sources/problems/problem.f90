!> @file problem.f90
!! @Copyright (c) 2023, The Neko Authors
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

!>
module problem
  use num_types, only: rp
  use fld_file_output, only: fld_file_output_t
  use design, only: design_t
  use utils, only: neko_error

  use vector, only: vector_t

  use objective, only: objective_t
  use constraint, only: constraint_t

  implicit none
  private

  !> implements the problem type.
  ! Currently very abstract, could include unsteady problems etc.
  ! Also, dependingo on the type of optimizer used, we may require
  ! different functionality.
  ! Right now, all that is required in base class is to init and
  ! evaluate the problem.
  type, abstract, public :: problem_t
     private

     !> The number of design variables
     integer :: n_design
     !> Number of objectives in the problem
     integer :: n_objectives
     !> Number of constraints in the problem
     integer :: n_constraints

     !> TODO
     ! we need a `objective_list` which is allocatable and contains a factory
     ! to fill itself up with from the JSON
     ! for now, I'm hardcoding these two
     class(objective_t), public, allocatable :: objective_function
     class(constraint_t), public, allocatable :: volume_constraint

     !> An output sampler for the problem. This should probably be an output
     !! controller at some point intead.
     type(fld_file_output_t), public :: output

   contains

     ! ----------------------------------------------------------------------- !
     ! Interfaces

     !> Constructor for physics of the problem
     procedure(problem_init), pass(this), public, deferred :: init
     !> Additional constructor specific to a design
     procedure(problem_init_design), pass(this), public, deferred :: init_design
     !> Destructor.
     procedure(problem_free), pass(this), deferred, public :: free

     !> Evaluate the optimization problem.
     !! This is the main function that evaluates the problem. It should be
     !! implemented in the derived classes.
     procedure(problem_compute), pass(this), public, deferred :: compute

     ! ----------------------------------------------------------------------- !
     ! Base class methods

     !> Constructor for the base class
     procedure, pass(this) :: init_base => problem_init_base
     !> Destructor for the base class
     procedure, pass(this) :: free_base => problem_free_base

     ! ----------------------------------------------------------------------- !
     ! Actual methods

     !> Sample the problem
     procedure, pass(this), public :: write => problem_write

     !> Update the objective function
     procedure, pass(this) :: update_objectives => &
          problem_update_objectives
     !> Update the volume constraint
     procedure, pass(this) :: update_constraints => &
          problem_update_constraints
     !> Update the objective sensitivities
     procedure, pass(this) :: update_objective_sensitivities => &
          problem_update_objective_sensitivities
     !> Update the constraint sensitivities
     procedure, pass(this) :: update_constraint_sensitivities => &
          problem_update_constraint_sensitivities

     !> Evaluate the objective.
     procedure, pass(this), public :: get_objective_value => &
          problem_get_objective_value
     !> Evaluate the constraints.
     procedure, pass(this), public :: get_constraint_values => &
          problem_get_constraint_values
     !  !> Evaluate the sensitivity of the objective.
     !  procedure, pass(this), public :: get_objective_sensitivities => &
     !       problem_get_objective_sensitivities
     !  !> Evaluate the sensitivity of the constraints.
     !  procedure, pass(this), public :: get_constraint_sensitivities => &
     !       problem_get_constraint_sensitivities


     !> Return the number of design variables.
     procedure, pass(this) :: get_n_design => problem_get_num_design_variables
     !> Return the number of objectives.
     procedure, pass(this) :: get_n_objectives => problem_get_num_objectives
     !> Return the number of constraints.
     procedure, pass(this) :: get_n_constraints => problem_get_num_constraints
  end type problem_t

  ! -------------------------------------------------------------------------- !
  ! Interfaces for the derived types

  abstract interface
     !> Constructor for physics of the problem.
     !! This is the main constructor for a problem. This should be defined in
     !! the derived types to initialize the problem. This is based on the
     !! abstract design type, We suggest that a swticth statement is used to
     !! initialize the problem based on the design type.
     subroutine problem_init(this)
       import problem_t
       class(problem_t), intent(inout) :: this
     end subroutine problem_init

     !> Additional constructor based on a design
     subroutine problem_init_design(this, design)
       import problem_t, design_t
       class(problem_t), intent(inout) :: this
       ! class(design_variable_t), intent(in) :: design
       ! we also only have the `topopt_design_t` but this should take the more
       ! abstract `design_variable_t` and initialize differently according to
       ! the type entering here.
       class(design_t), target, intent(inout) :: design

       ! This is confusing to me..
       ! The `problem` and the `design` seem very coupled in my mind.
       ! I want to argue it's coupled one way, since the problem depends on the
       ! design representation.

       ! In principle we could have our design represented with
       ! - splines
       ! - levelset
       ! - etc
       !
       ! BUT, for density based topology optimization, because we get all our mesh
       ! info etc from neko, our design representation is based on the fluid.
       ! (of course this isn't 100% true, it's just the dofmap. We could define
       ! our design on a different set of basis functions too... but I guess that
       ! is rather far out of scope now...)
       !
       ! So it's sort of coupled both ways.. :/
       !
       ! Tim you may need to untagle this, for now I dont see an option other than
       ! - initialising the fluid first.
       !
       ! - The initializing the design
       !
       ! - Then coming here and intializing the impact of the design on the fluid
       !
     end subroutine problem_init_design


     !> Compute the problem
     subroutine problem_compute(this, design)
       import problem_t
       import design_t

       class(problem_t), intent(inout) :: this
       class(design_t), intent(inout) :: design

     end subroutine problem_compute

     !> Destructor
     subroutine problem_free(this)
       import problem_t
       class(problem_t), intent(inout) :: this
     end subroutine problem_free
  end interface

contains

  ! -------------------------------------------------------------------------- !
  ! Base class methods

  !> Constructor for the base class
  subroutine problem_init_base(this, n_design, n_objectives, n_constraints)
    class(problem_t), intent(inout) :: this
    integer, intent(in) :: n_design, n_objectives, n_constraints

    this%n_design = n_design
    this%n_objectives = n_objectives
    this%n_constraints = n_constraints
  end subroutine problem_init_base

  subroutine problem_free_base(this)
    class(problem_t), intent(inout) :: this
  end subroutine problem_free_base

  !> Sample the fields/design.
  subroutine problem_write(this, idx)
    class(problem_t), intent(inout) :: this
    integer, intent(in) :: idx

    call this%output%sample(real(idx, kind=rp))
  end subroutine problem_write


  ! -------------------------------------------------------------------------- !
  ! Updater methods

  subroutine problem_update_objectives(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design

    call this%objective_function%update_value(design)
  end subroutine problem_update_objectives

  subroutine problem_update_constraints(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design

    call this%volume_constraint%update_value(design)
  end subroutine problem_update_constraints

  subroutine problem_update_objective_sensitivities(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design

    call this%objective_function%update_sensitivity(design)
  end subroutine problem_update_objective_sensitivities

  subroutine problem_update_constraint_sensitivities(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(inout) :: design

    call this%volume_constraint%update_sensitivity(design)
  end subroutine problem_update_constraint_sensitivities

  ! -------------------------------------------------------------------------- !
  ! Getter methods

  subroutine problem_get_objective_value(this, value)
    class(problem_t), intent(inout) :: this
    real(kind=rp), intent(out) :: value

    value = this%objective_function%value
  end subroutine problem_get_objective_value

  subroutine problem_get_constraint_values(this, values)
    class(problem_t), intent(inout) :: this
    real(kind=rp), intent(out) :: values

    values = this%volume_constraint%value

  end subroutine problem_get_constraint_values

  ! -------------------------------------------------------------------------- !
  ! Simple getters

  !> Return the number of design variables.
  pure function problem_get_num_design_variables(this) result(n)
    class(problem_t), intent(in) :: this
    integer :: n

    n = this%n_design
  end function problem_get_num_design_variables

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

end module problem
