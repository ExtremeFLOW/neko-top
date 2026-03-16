!> @file neko_objective.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
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

!> Implements the `neko_objective_t` type.
module neko_objective
  use num_types, only: rp
  use objective, only: objective_t
  use design, only: design_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use time_state, only: time_state_t
  use vector_math, only: vector_copy, vector_add2s1
  implicit none
  private

  public :: neko_objective_t

  !> Intermediate base type for objectives requiring a neko simulation.
  type, public, abstract, extends(objective_t) :: neko_objective_t
     real(kind=rp) :: start_time = 0.0_rp
     real(kind=rp) :: end_time = huge(0.0_rp)
     type(time_state_t), pointer :: time => null()
   contains
     procedure, pass(this) :: init_time_window_json => &
          neko_objective_init_time_window_json
     procedure, pass(this) :: bind_time => neko_objective_bind_time
     procedure, pass(this) :: accumulate_value => &
          neko_objective_accumulate_value
     procedure, pass(this) :: accumulate_sensitivity => &
          neko_objective_accumulate_sensitivity
     procedure, private, pass(this) :: is_active => neko_objective_is_active
  end type neko_objective_t

contains

  !> Read the active time window for this objective from JSON.
  subroutine neko_objective_init_time_window_json(this, json)
    class(neko_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    call json_get_or_default(json, "start_time", this%start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", this%end_time, huge(0.0_rp))
  end subroutine neko_objective_init_time_window_json

  !> Bind the objective to the primal simulation time state.
  subroutine neko_objective_bind_time(this, time)
    class(neko_objective_t), intent(inout) :: this
    type(time_state_t), target, intent(inout) :: time

    this%time => time
  end subroutine neko_objective_bind_time

  !> Accumulate the value only while the objective time window is active.
  subroutine neko_objective_accumulate_value(this, design, dt)
    class(neko_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    real(kind=rp), intent(in) :: dt

    if (.not. this%is_active()) return

    this%value_old = this%value
    call this%update_value(design)
    this%value = this%value_old + this%value * dt
  end subroutine neko_objective_accumulate_value

  !> Accumulate the sensitivity while the objective time window is active.
  subroutine neko_objective_accumulate_sensitivity(this, design, dt)
    class(neko_objective_t), intent(inout) :: this
    class(design_t), intent(in) :: design
    real(kind=rp), intent(in) :: dt

    if (.not. this%is_active()) return

    call vector_copy(this%sensitivity_old, this%sensitivity)
    call this%update_sensitivity(design)
    call vector_add2s1(this%sensitivity, this%sensitivity_old, dt)
  end subroutine neko_objective_accumulate_sensitivity

  !> Return true when the current primal time is inside the objective window.
  logical function neko_objective_is_active(this)
    class(neko_objective_t), intent(in) :: this

    if (.not. associated(this%time)) then
       neko_objective_is_active = .true.
       return
    end if

    neko_objective_is_active = this%time%t .ge. this%start_time .and. &
         this%time%t .le. this%end_time
  end function neko_objective_is_active

end module neko_objective
