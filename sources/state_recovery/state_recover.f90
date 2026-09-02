!> @file state_recover.f90
!! @copyright
!! Copyright (c) 2025-2026, The Neko-TOP Authors
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
!> @brief Abstract interface for state recovery strategies.
module state_recover
  use case, only: case_t
  use json_file_module, only: json_file
  implicit none
  private

  !> Abstract base type for state recovery implementations.
  type, abstract, public :: state_recover_t
     private
     !> Case associated with this recovery strategy at initialization.
     class(case_t), pointer, public :: neko_case => null()
     !> Number of states recorded by this recovery strategy.
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
     !> Initialize state recovery from JSON parameters.
     !! @param[inout] this State recovery instance.
     !! @param[inout] neko_case Case data structure.
     !! @param[inout] params JSON parameters.
     subroutine state_recover_init(this, neko_case, params)
       import state_recover_t, case_t, json_file
       class(state_recover_t), intent(inout) :: this
       class(case_t), target, intent(inout) :: neko_case
       type(json_file), target, intent(inout) :: params
     end subroutine state_recover_init

     !> Free state recovery resources.
     !! @param[inout] this State recovery instance.
     subroutine state_recover_free(this)
       import state_recover_t
       class(state_recover_t), intent(inout) :: this
     end subroutine state_recover_free

     !> Reset state recovery bookkeeping.
     !! @param[inout] this State recovery instance.
     subroutine state_recover_reset(this)
       import state_recover_t
       class(state_recover_t), intent(inout) :: this
     end subroutine state_recover_reset

     !> Save forward state for recovery.
     !! @param[inout] this State recovery instance.
     subroutine state_recover_save(this)
       import state_recover_t
       class(state_recover_t), intent(inout) :: this
     end subroutine state_recover_save

     !> Restore forward state for adjoint.
     !! @param[inout] this State recovery instance.
     !! @param[in] tstep Timestep to restore.
     subroutine state_recover_restore(this, tstep)
       import state_recover_t
       class(state_recover_t), intent(inout) :: this
       integer, intent(in) :: tstep
     end subroutine state_recover_restore
  end interface

contains

  !> Get the number of states recorded by this recovery strategy.
  !! @param[in] this State recovery instance.
  pure function state_recover_get_n_timesteps(this) result(n)
    class(state_recover_t), intent(in) :: this
    integer :: n

    n = this%n_timesteps
  end function state_recover_get_n_timesteps

  !> Set the number of states recorded by this recovery strategy.
  !! @param[inout] this State recovery instance.
  !! @param[in] n Number of recorded states.
  subroutine state_recover_set_n_timesteps(this, n)
    class(state_recover_t), intent(inout) :: this
    integer, intent(in) :: n

    this%n_timesteps = n
  end subroutine state_recover_set_n_timesteps

end module state_recover
