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
  use design, only: design_t
  use simulation, only: simulation_t

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

     !> The simulation.
     type(simulation_t), public :: simulation

   contains

     ! ----------------------------------------------------------------------- !
     ! Interfaces

     !> Constructor for physics of the problem.
     procedure(problem_init), pass(this), public, deferred :: init
     !> Additional constructor specific to a design.
     procedure(problem_init_design), pass(this), public, deferred :: init_design
     !> Destructor.
     procedure(problem_free), pass(this), deferred, public :: free

     !> Evaluate the optimization problem.
     procedure(problem_compute), pass(this), deferred :: compute
     !> Sample the problem
     procedure(problem_sample), pass(this), deferred :: sample

  end type problem_t

  !> Constructor for physics of the problem
  abstract interface
     subroutine problem_init(this)
       import problem_t
       class(problem_t), target, intent(inout) :: this
     end subroutine problem_init
  end interface

  !> Additional constructor based on a design
  abstract interface
     subroutine problem_init_design(this, design)
       import problem_t, design_t
       class(problem_t), intent(inout) :: this
       ! class(design_variable_t), intent(in) :: design
       ! we also only have the `design_t` but this should take the more
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
  end interface

  !> Destructor
  abstract interface
     subroutine problem_free(this)
       import problem_t
       class(problem_t), intent(inout) :: this
       ! TODO
     end subroutine problem_free
  end interface

  !> Compute
  abstract interface
     subroutine problem_compute(this)
       import problem_t
       class(problem_t), intent(inout) :: this

     end subroutine problem_compute
  end interface

  !> Sample
  abstract interface
     subroutine problem_sample(this, t)
       import problem_t, rp
       class(problem_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t

     end subroutine problem_sample
  end interface
end module problem
