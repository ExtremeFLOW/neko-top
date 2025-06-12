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

module example_problem
  use num_types, only: rp

  use objective, only: objective_t
  use constraint, only: constraint_t

  use design, only: design_t
  use math, only: glsum
  use json_module, only: json_file
  use vector, only: vector_t

  use device, only: device_memcpy, DEVICE_TO_HOST
  use neko_config, only: NEKO_BCKND_DEVICE

  implicit none
  private

  type, public, extends(objective_t) :: mma_obj

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json => mma_obj_init_from_json
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_components => &
          mma_obj_init_from_components
     !> Destructor.
     procedure, public, pass(this) :: free => mma_obj_free

     !> Computes the value of the objective function.
     procedure, public, pass(this) :: update_value => mma_obj_update_value

     !> Computes the sensitivity with respect to the coefficient $\chi$.
     procedure, public, pass(this) :: update_sensitivity => &
          mma_obj_update_sensitivity

  end type mma_obj

  type, public, extends(constraint_t) :: mma_con

     real(kind=rp) :: sign = 1.0_rp

   contains
     !> The common constructor using a JSON object.
     procedure, public, pass(this) :: init_json => mma_con_init_from_json
     !> The actual constructor.
     procedure, public, pass(this) :: init_from_components => &
          mma_con_init_from_components
     !> Destructor.
     procedure, public, pass(this) :: free => mma_con_free

     !> Computes the value of the constraint function.
     procedure, public, pass(this) :: update_value => mma_con_update_value

     !> Computes the sensitivity with respect to the coefficient $\chi$.
     procedure, public, pass(this) :: update_sensitivity => &
          mma_con_update_sensitivity
  end type mma_con

contains

  ! ========================================================================== !
  ! Methods for the Objective Function

  subroutine mma_obj_init_from_json(this, json, design)
    class(mma_obj), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    character(len=256), parameter :: name = 'mma_obj'

    call this%init_from_components(name, design)

  end subroutine mma_obj_init_from_json

  subroutine mma_obj_init_from_components(this, name, design)
    class(mma_obj), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(design_t), intent(in) :: design

    call this%init_base(name, design%size(), 1.0_rp)

  end subroutine mma_obj_init_from_components

  subroutine mma_obj_free(this)
    class(mma_obj), intent(inout) :: this
    call this%free_base()
  end subroutine mma_obj_free

  subroutine mma_obj_update_value(this, design)
    class(mma_obj), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(vector_t) :: design_values

    design_values = design%get_values()

    if (NEKO_BCKND_DEVICE .eq. 1) then
        call device_memcpy(design_values%x, &
             design_values%x_d, design%size(), &
             DEVICE_TO_HOST, sync = .false.)
    end if
    this%value = glsum(design_values%x, design%size()) &
         / real(design%size_global(), kind=rp)

  end subroutine mma_obj_update_value

  subroutine mma_obj_update_sensitivity(this, design)
    class(mma_obj), intent(inout) :: this
    class(design_t), intent(in) :: design

    this%sensitivity = 1.0_rp / real(design%size_global(), kind=rp)

  end subroutine mma_obj_update_sensitivity

  ! ========================================================================== !
  ! Methods for the Constraint Function

  subroutine mma_con_init_from_json(this, json, design)
    class(mma_con), intent(inout) :: this
    type(json_file), intent(inout) :: json
    class(design_t), intent(in) :: design
    character(len=256), parameter :: name = 'mma_con'
    integer :: sign

    call this%init_from_components(name, design, sign)

  end subroutine mma_con_init_from_json

  subroutine mma_con_init_from_components(this, name, design, sign)
    class(mma_con), intent(inout) :: this
    character(len=*), intent(in) :: name
    class(design_t), intent(in) :: design
    integer, intent(in), optional :: sign
    integer :: sign_ = 1

    call this%init_base(name, design%size())

    if (present(sign)) sign_ = sign
    if (sign_ .lt. 0) this%sign = -1.0_rp
    if (sign_ .ge. 0) this%sign = 1.0_rp

  end subroutine mma_con_init_from_components

  subroutine mma_con_free(this)
    class(mma_con), intent(inout) :: this
    call this%free_base()
  end subroutine mma_con_free

  subroutine mma_con_update_value(this, design)
    class(mma_con), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(vector_t) :: difference

    difference = design%get_values() - design%get_x()
    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(difference%x, &
            difference%x_d, design%size(), &
            DEVICE_TO_HOST, sync = .false.)
    end if

    difference%x = difference%x * difference%x

    ! Note that we divide by n to make a fair comparison when ploting the 
    ! scaling plots with respect to n (KKT tol is effected by this)
    this%value = this%sign * glsum(difference%x, design%size()) / &
         real(design%size_global(), kind=rp)

  end subroutine mma_con_update_value

  subroutine mma_con_update_sensitivity(this, design)
    class(mma_con), intent(inout) :: this
    class(design_t), intent(in) :: design
    type(vector_t) :: difference

    difference = design%get_values() - design%get_x()

    this%sensitivity = this%sign * 2.0_rp * difference / &
         real(design%size_global(), kind=rp)

  end subroutine mma_con_update_sensitivity
end module example_problem
