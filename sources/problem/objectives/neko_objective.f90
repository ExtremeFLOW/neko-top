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
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  implicit none
  private

  public :: neko_objective_t

  !> Intermediate base type for objectives requiring a neko simulation.
  type, public, abstract, extends(objective_t) :: neko_objective_t
     real(kind=rp) :: start_time = 0.0_rp
     real(kind=rp) :: end_time = huge(0.0_rp)
   contains
     procedure, pass(this) :: init_time_window_json => &
          neko_objective_init_time_window_json
  end type neko_objective_t

contains

  !> Read the active time window for this objective from JSON.
  subroutine neko_objective_init_time_window_json(this, json)
    class(neko_objective_t), intent(inout) :: this
    type(json_file), intent(inout) :: json

    call json_get_or_default(json, "start_time", this%start_time, 0.0_rp)
    call json_get_or_default(json, "end_time", this%end_time, huge(0.0_rp))
  end subroutine neko_objective_init_time_window_json

end module neko_objective
