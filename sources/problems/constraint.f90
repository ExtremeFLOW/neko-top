!> @file constraint.f90
!! @copyright (c) 2025, The Neko Authors
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

!> Implements the `constraint_t` type.
module constraint
  use base_functional, only: base_functional_t
  use simulation_m, only: simulation_t
  use design, only: design_t
  use num_types, only: rp
  use point_zone_registry, only: neko_point_zone_registry
  use json_module, only: json_file
  implicit none
  private

  public :: constraint_t, constraint_factory, constraint_wrapper_t

  !> The abstract constraint type.
  !!
  !! This is the base class for constraints, which is a type of base functional.
  type, abstract, extends(base_functional_t) :: constraint_t

   contains

     !> Initializer for the base class
     procedure, pass(this) :: init_base => constraint_init_base
     !> Destructor of the base class
     procedure, pass(this) :: free_base => constraint_free_base

  end type constraint_t

  !> Wrapper for constraints for use in lists.
  type :: constraint_wrapper_t
     class(constraint_t), allocatable :: constraint
   contains
     procedure, pass(this) :: free => constraint_wrapper_free
  end type constraint_wrapper_t

  ! -------------------------------------------------------------------------- !
  ! Explicit interfaces

  !> Factory function interface
  interface
     module subroutine constraint_factory(object, json, design, simulation)
       class(constraint_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: json
       class(design_t), intent(in) :: design
       type(simulation_t), target, intent(inout) :: simulation
     end subroutine constraint_factory
  end interface

contains

  ! -------------------------------------------------------------------------- !
  ! Implementations for the base class

  !> Initialize the constraint base class.
  !! @param this The constraint object.
  !! @param name The name of the constraint.
  !! @param design_size The number of design variables.
  !! @param mask_name The name design the mask. [optional]
  subroutine constraint_init_base(this, name, design_size, mask_name)
    class(constraint_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: design_size
    character(len=*), intent(in), optional :: mask_name

    this%name = name
    this%value = 0.0_rp
    call this%sensitivity%init(design_size)

    this%has_mask = .false.
    if (trim(mask_name) .ne. "") then
       this%has_mask = .true.
       this%mask => neko_point_zone_registry%get_point_zone(mask_name)
    end if

  end subroutine constraint_init_base

  !> Free the base class
  subroutine constraint_free_base(this)
    class(constraint_t), target, intent(inout) :: this

    this%value = 0.0_rp
    call this%sensitivity%free()
    this%has_mask = .false.
    if (associated(this%mask)) nullify(this%mask)

  end subroutine constraint_free_base

  ! -------------------------------------------------------------------------- !
  ! Implementations for the wrapper

  !> Free the constraint wrapper.
  subroutine constraint_wrapper_free(this)
    class(constraint_wrapper_t), intent(inout) :: this
    if (allocated(this%constraint)) then
       call this%constraint%free()
       deallocate(this%constraint)
    end if
  end subroutine constraint_wrapper_free

end module constraint

