!> @file state_recover_factory.f90
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
!> @brief Factory for state recovery implementations.
module state_recover_factory
  use case, only: case_t
  use json_file_module, only: json_file
  use json_utils, only: json_get_or_default
  use state_recover, only: state_recover_t
  use simulation_checkpoint, only: simulation_checkpoint_t
#if HAVE_ADIOS2
  use simulation_POD_state_recover, only: POD_state_recover_t
#endif
  use utils, only: neko_error
  implicit none
  private

  public :: state_recover_create

contains

  !> Create and initialize a state recovery object from JSON parameters.
  !! @param[inout] recover Allocatable state recovery instance.
  !! @param[inout] neko_case Case data structure.
  !! @param[inout] params JSON parameters for state recovery.
  subroutine state_recover_create(recover, neko_case, params)
    class(state_recover_t), allocatable, intent(inout) :: recover
    class(case_t), target, intent(inout) :: neko_case
    type(json_file), intent(inout) :: params
    character(len=:), allocatable :: recover_type

    call json_get_or_default(params, "type", recover_type, "checkpoint")

    if (allocated(recover)) then
       call recover%free()
       deallocate(recover)
    end if

    select case (trim(recover_type))
    case ("checkpoint", "simulation_checkpoint")
       allocate(simulation_checkpoint_t :: recover)
    case ("pod")
#if HAVE_ADIOS2
       allocate(POD_state_recover_t :: recover)
#else
       call neko_error("POD state recovery requires ADIOS2. Rebuild with " // &
            "ADIOS2 enabled.")
#endif
    case default
       call neko_error("Unknown state recover type: " // trim(recover_type))
    end select

    call recover%init(neko_case, params)
  end subroutine state_recover_create

end module state_recover_factory
