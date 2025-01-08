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

!> Implements the `base_functional_t` type.
module base_functional
  use num_types, only: rp, dp
  use field, only: field_t
  use simulation, only: simulation_t
  use design, only: design_t
  use topopt_design, only: topopt_design_t
  use point_zone_registry, only: neko_point_zone_registry
  use point_zone, only: point_zone_t
  use utils, only: neko_error
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
  !! simulation components that are required to evaluate the base functional. All of
  !! which should be prepared in the `init` method.
  type, abstract, public :: base_functional_t

     !> Value of the base_functional
     real(kind=rp) :: value
     !> Sensitivity field
     type(field_t) :: sensitivity

     !> containing a mask?
     logical :: if_mask
     !> A mask for where the objective function is evaluated
     class(point_zone_t), pointer :: mask

   contains
     !> Initializers of the base class
     procedure, pass(this) :: init_base => functional_init_base
     procedure, pass(this) :: free_base => functional_free_base

     !> Constructor
     procedure(functional_init), pass(this), deferred :: init
     !> Destructor
     procedure(functional_free), pass(this), deferred :: free


     !> Compute the objective function
     procedure(functional_compute), pass(this), deferred :: compute
     !> Compute the sensitivity
     procedure(functional_sensitivity), pass(this), deferred :: &
          compute_sensitivity
  end type base_functional_t

  ! -------------------------------------------------------------------------- !
  ! Interfaces for the derived types, these are the constructors for the
  ! different types of objective functions.
  abstract interface

     !> Initialize the objective function
     subroutine functional_init(this, design, simulation)
       import base_functional_t, design_t, simulation_t
       class(base_functional_t), intent(inout) :: this
       type(simulation_t), target, intent(inout) :: simulation
       class(design_t), intent(in) :: design
     end subroutine functional_init

     !> Destructor
     subroutine functional_free(this)
       import base_functional_t
       class(base_functional_t), intent(inout) :: this
     end subroutine functional_free

     !> Compute the objective function
     subroutine functional_compute(this, design)
       import base_functional_t, design_t
       class(base_functional_t), intent(inout) :: this
       class(design_t), intent(in) :: design
     end subroutine functional_compute

     !> Compute the sensitivity
     subroutine functional_sensitivity(this, design)
       import base_functional_t, design_t
       class(base_functional_t), intent(inout) :: this
       class(design_t), intent(in) :: design
     end subroutine functional_sensitivity

  end interface

contains

  ! -------------------------------------------------------------------------- !
  ! Implementations for the base class

  !> Initialize the base class
  !! @param design The design variable
  !! @param[optional] if_mask Whether the base_functional is masked
  !! @param[optional] mask_name The name design the mask
  subroutine functional_init_base(this, design, if_mask, mask_name)
    class(base_functional_t), target, intent(inout) :: this
    class(design_t), intent(in) :: design
    logical, intent(in), optional :: if_mask
    character(len=*), intent(in), optional :: mask_name

    this%value = 0.0_rp
    select type(design)
      type is (topopt_design_t)
       call this%sensitivity%init(design%design_indicator%dof)
      class default
       call neko_error('Unknown design type')
    end select

    if (present(if_mask)) this%if_mask = if_mask
    if (.not. present(if_mask)) this%if_mask = .false.

    if (this%if_mask) then
       if (.not. present(mask_name)) then
          call neko_error('Mask name not provided')
       end if

       this%mask => neko_point_zone_registry%get_point_zone(mask_name)
    end if

  end subroutine functional_init_base

  !> Free the base class
  subroutine functional_free_base(this)
    class(base_functional_t), target, intent(inout) :: this

    if (associated(this%mask)) nullify(this%mask)

    call this%sensitivity%free()

  end subroutine functional_free_base

end module base_functional

