!> @file state_recover.f90
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
!
module state_recover
  use case, only: case_t
  use json_file_module, only: json_file
  use time_state, only: time_state_t
  implicit none
  private

  type, abstract, public :: state_recover_t
     private
     !> number of time steps in the forward simulation
     integer :: n_timesteps = 0
   contains
     procedure(state_recover_init), pass(this), public, deferred :: init
     procedure(state_recover_free), pass(this), public, deferred :: free
     procedure(state_recover_reset), pass(this), public, deferred :: reset
     procedure(state_recover_save), pass(this), public, deferred :: save
     procedure(state_recover_restore), pass(this), public, deferred :: restore
     procedure, pass(this), public :: get_n_timesteps => &
          state_recover_get_n_timesteps
     procedure, pass(this), public :: set_n_timesteps => &
          state_recover_set_n_timesteps
  end type state_recover_t

  abstract interface
     subroutine state_recover_init(this, neko_case, params)
       import state_recover_t, case_t, json_file
       class(state_recover_t), intent(inout) :: this
       class(case_t), target, intent(inout) :: neko_case
       type(json_file), target, intent(inout) :: params
     end subroutine state_recover_init

     subroutine state_recover_free(this)
       import state_recover_t
       class(state_recover_t), intent(inout) :: this
     end subroutine state_recover_free

     subroutine state_recover_reset(this)
       import state_recover_t
       class(state_recover_t), intent(inout) :: this
     end subroutine state_recover_reset

     subroutine state_recover_save(this, neko_case, time)
       import state_recover_t, case_t, time_state_t
       class(state_recover_t), intent(inout) :: this
       class(case_t), intent(inout) :: neko_case
       type(time_state_t), intent(in) :: time
     end subroutine state_recover_save

     subroutine state_recover_restore(this, neko_case, time)
       import state_recover_t, case_t, time_state_t
       class(state_recover_t), intent(inout) :: this
       class(case_t), target, intent(inout) :: neko_case
       type(time_state_t), intent(in) :: time
     end subroutine state_recover_restore
  end interface

contains

  pure function state_recover_get_n_timesteps(this) result(n)
    class(state_recover_t), intent(in) :: this
    integer :: n

    n = this%n_timesteps
  end function state_recover_get_n_timesteps

  subroutine state_recover_set_n_timesteps(this, n)
    class(state_recover_t), intent(inout) :: this
    integer, intent(in) :: n

    this%n_timesteps = n
  end subroutine state_recover_set_n_timesteps


end module state_recover
