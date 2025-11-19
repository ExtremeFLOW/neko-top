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

module example_problem_mma
  use num_types, only: rp

  use objective, only: objective_t
  use constraint, only: constraint_t

  use design, only: design_t
  use json_module, only: json_file
  use vector, only: vector_t

  use device, only: device_memcpy, DEVICE_TO_HOST
  use neko_config, only: NEKO_BCKND_DEVICE
  use vector_math, only: vector_sub2, vector_col2, vector_glsum, vector_cmult

  use vector_scratch_registry, only: neko_vector_scratch_registry

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

    type(vector_t), pointer :: difference
    type(vector_t), pointer :: x_coordinate
    integer :: ind(2)

    call neko_vector_scratch_registry%request_vector(design%size(), &
         x_coordinate, ind(1))
    call neko_vector_scratch_registry%request_vector(design%size(), &
         difference, ind(2))

    call design%get_x(x_coordinate)
    call design%get_values(difference)
    call vector_sub2(difference, x_coordinate, design%size())

    call vector_col2(difference, difference, design%size())

    this%value = vector_glsum(difference, design%size()) / &
         real(design%size_global(), kind=rp)

    call neko_vector_scratch_registry%relinquish_vector(ind)
  end subroutine mma_obj_update_value

  subroutine mma_obj_update_sensitivity(this, design)
    class(mma_obj), intent(inout) :: this
    class(design_t), intent(in) :: design

    type(vector_t), pointer :: difference
    type(vector_t), pointer :: x_coordinate
    integer :: ind(2)

    call neko_vector_scratch_registry%request_vector(design%size(), &
         x_coordinate, ind(1))
    call neko_vector_scratch_registry%request_vector(design%size(), &
         difference, ind(2))

    call design%get_x(x_coordinate)
    call design%get_values(difference)
    call vector_sub2(difference, x_coordinate, design%size())

    call vector_cmult(difference , 2.0_rp / &
         real(design%size_global(), kind=rp), design%size())
    this%sensitivity = difference

    call neko_vector_scratch_registry%relinquish_vector(ind)
  end subroutine mma_obj_update_sensitivity

end module example_problem_mma
