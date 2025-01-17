! Copyright (c) 2025, The Neko Authors
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

!> Implements the `objective_t` type.
module objective
  use base_functional, only: base_functional_t
  use simulation, only: simulation_t
  use design, only: design_t
  use num_types, only: rp
  use point_zone_registry, only: neko_point_zone_registry
  implicit none
  private

  public :: objective_t, objective_factory, objective_wrapper_t

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
  type, abstract, extends(base_functional_t) :: objective_t
     real(kind=rp) :: weight = 1.0_rp

   contains

     !> Initialize the base class
     procedure, pass(this) :: init_base => objective_init_base
     !> Free the base class
     procedure, pass(this) :: free_base => objective_free_base

  end type objective_t

  !> Wrapper for objectives for use in lists.
  type :: objective_wrapper_t
     class(objective_t), allocatable :: objective
   contains
     procedure, pass(this) :: free => objective_wrapper_free
  end type objective_wrapper_t

  ! -------------------------------------------------------------------------- !
  ! Explicit interfaces

  !> Factory function interface
  interface
     module subroutine objective_factory(object, type, design, simulation)
       class(objective_t), allocatable, intent(inout) :: object
       character(len=*), intent(in) :: type
       class(design_t), intent(in) :: design
       type(simulation_t), target, intent(inout) :: simulation
     end subroutine objective_factory
  end interface

contains

  ! -------------------------------------------------------------------------- !
  ! Implementations for the base class

  !> Initialize the objective base class.
  !! @param design_size The number of design variables.
  !! @param weight The weight of the objective function.
  !! @param[optional] mask_name The name design the mask.
  subroutine objective_init_base(this, design_size, weight, mask_name)
    class(objective_t), intent(inout) :: this
    integer, intent(in) :: design_size
    real(kind=rp), intent(in) :: weight
    character(len=*), intent(in), optional :: mask_name

    this%value = 0.0_rp
    call this%sensitivity%init(design_size)

    this%weight = weight

    this%has_mask = .false.
    if (trim(mask_name) .ne. "") then
       this%has_mask = .true.
       this%mask => neko_point_zone_registry%get_point_zone(mask_name)
    end if

  end subroutine objective_init_base

  !> Free the base class
  subroutine objective_free_base(this)
    class(objective_t), target, intent(inout) :: this

    this%value = 0.0_rp
    call this%sensitivity%free()
    this%weight = 1.0_rp
    this%has_mask = .false.
    if (associated(this%mask)) nullify(this%mask)

  end subroutine objective_free_base

  ! -------------------------------------------------------------------------- !
  ! Implementations for the wrapper

  !> Free the objective wrapper.
  subroutine objective_wrapper_free(this)
    class(objective_wrapper_t), intent(inout) :: this
    if (allocated(this%objective)) then
       call this%objective%free()
       deallocate(this%objective)
    end if
  end subroutine objective_wrapper_free

end module objective

