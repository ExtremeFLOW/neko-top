! Copyright (c) 2024, The Neko Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

!> Defines the abstract the `base_functional_t` type.
module base_functional
  use design, only: design_t
  use json_module, only: json_file
  use num_types, only: rp
  use point_zone, only: point_zone_t
  use simulation, only: simulation_t
  use vector, only: vector_t
  implicit none
  private

  !> The base functional type
  !!
  !! This is the base class for objectives and constraints alike.
  !! A base functional should be able to evaluate itself and its sensitivity
  !! with respect to the design variables.
  !!
  !! The base functional is also responsible for managing the adjoint forcing
  !! terms that are required for the adjoint problem, any source terms
  !! simulation components that are required to evaluate the base functional.
  !! All of which should be prepared in the `init` method.
  type, abstract, public :: base_functional_t

     !> Value of the base_functional
     real(kind=rp) :: value
     !> Sensitivity field
     type(vector_t) :: sensitivity
     !> Name of constraint/objective in the logfile
     character(len=25) :: log_name
     !> containing a mask
     logical :: has_mask
     !> A mask for where the objective function is evaluated
     class(point_zone_t), pointer :: mask => null()

   contains

     ! ----------------------------------------------------------------------- !
     ! Derived class interfaces

     !> Constructor
     procedure(functional_init), pass(this), deferred :: init_json
     !> Destructor
     procedure(functional_free), pass(this), deferred :: free

     !> Update the value of the function
     procedure(functional_update_value), pass(this), deferred :: update_value
     !> Update the sensitivity of the function
     procedure(functional_update_sensitivity), pass(this), deferred :: &
          update_sensitivity
  end type base_functional_t

  ! -------------------------------------------------------------------------- !
  ! Interface specifications for the derived types, these are the constructors
  ! for the different types of objective functions.

  abstract interface

     !> Initialize the objective function
     subroutine functional_init(this, json, design, simulation)
       import base_functional_t, design_t, simulation_t, json_file
       class(base_functional_t), intent(inout) :: this
       type(json_file), intent(inout) :: json
       class(design_t), intent(in) :: design
       type(simulation_t), target, intent(inout) :: simulation
     end subroutine functional_init

     !> Destructor
     subroutine functional_free(this)
       import base_functional_t
       class(base_functional_t), intent(inout) :: this
     end subroutine functional_free

     !> Compute the objective function
     subroutine functional_update_value(this, design)
       import base_functional_t, design_t
       class(base_functional_t), intent(inout) :: this
       class(design_t), intent(in) :: design
     end subroutine functional_update_value

     !> Compute the sensitivity
     subroutine functional_update_sensitivity(this, design)
       import base_functional_t, design_t
       class(base_functional_t), intent(inout) :: this
       class(design_t), intent(in) :: design
     end subroutine functional_update_sensitivity

  end interface

end module base_functional
