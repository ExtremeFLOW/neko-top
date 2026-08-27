!> @file checkpoint_linear.f90
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

!> Implementation for the Linear Checkpointing algorithm.
!! In this case, we save the state of the simulation every `n_saves_memory`
!! time steps. When restoring to a given time step, we load the nearest
!! checkpoint and then we fill our cache with the following
!! `n_saves_memory` time steps. Finally, we copy the required time step from
!! our cache.
!!
!! This algorithm is the simplest one and do a minimum of re-computation. But
!! requires large amounts of memory and disk space.
submodule (simulation_checkpoint) checkpoint_linear
  use num_types, only: dp
  use simulation, only: simulation_step, simulation_restart
  use time_step_controller, only: time_step_controller_t
  use profiler, only: profiler_start_region, profiler_end_region

contains

  !> Save the current state of the simulation in a linear fashion.
  !! We save every `n_saves_memory` time steps to disc and we always save
  !! any timestep leading up to the `first_valid_timestep` time steps to disc.
  module subroutine checkpoint_save_linear(this, neko_case)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    integer :: index, tstep, counter, n_total
    real(kind=rp) :: time

    time = neko_case%time%t
    tstep = neko_case%time%tstep

    ! We save to disc only every n_saves_memory time steps
    index = modulo(tstep, this%n_saves_memory)
    if (index .eq. 0 .or. tstep .le. this%first_valid_timestep) then
       this%loaded_checkpoint = tstep

       counter = determine_counter(tstep, this%n_saves_memory, &
            this%first_valid_timestep)

       call this%chkp_output%set_counter(counter)
       call profiler_start_region("Checkpoint write to disk")
       call this%chkp_output%sample(time)
       call profiler_end_region("Checkpoint write to disk")
       this%n_saves_disc = this%n_saves_disc + 1
    end if

    ! Only save to RAM from the last disc checkpoint to the end of the forward
    ! simulation. With fixed timesteps, the total count and the last disc-save
    ! timestep are known from the time object.
    ! Note: the plus 0.5 is to round up to the next integer, as the division can
    ! be slightly smaller than the actual number of steps due to floating point
    ! errors.
    n_total = int(((neko_case%time%end_time - neko_case%time%start_time) &
         / neko_case%time%dt) + 0.5_rp)
    if (tstep .ge. n_total - modulo(n_total, this%n_saves_memory)) then
       call this%save_data(index + 1)
    end if

  end subroutine checkpoint_save_linear

  !> Restore the forward simulation state in a linear fashion.
  !! If the requested time step is not in memory, we load the nearest
  !! checkpoint from disc and then we step forward in time to fill our cache.
  !! Finally, we copy the requested time step from our cache.
  module subroutine checkpoint_restore_linear(this, neko_case, tstep)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: tstep
    type(time_step_controller_t) :: dt_controller
    real(kind=dp) :: loop_start
    integer :: k, previous_save, next_save, local_idx, counter

    loop_start = MPI_WTIME()

    ! Determine the nearest save states on both sides
    previous_save = tstep - modulo(tstep, this%n_saves_memory)
    next_save = previous_save + this%n_saves_memory

    ! Before the first valid state we always load from disc and we do not step
    ! forward in time.
    if (tstep .lt. this%first_valid_timestep) then
       previous_save = tstep
       next_save = previous_save + 1
    else if (previous_save .lt. this%first_valid_timestep) then
       previous_save = this%first_valid_timestep
    end if

    ! Load a new batch of checkpoints if needed
    if (this%loaded_checkpoint .ne. previous_save) then

       ! Restart the simulation form the checkpoint file
       counter = determine_counter(previous_save, this%n_saves_memory, &
            this%first_valid_timestep)
       call this%chkp_output%set_counter(counter)
       call profiler_start_region("Checkpoint read from disk")
       call this%chkp_output%file_%read(neko_case%chkp)
       call profiler_end_region("Checkpoint read from disk")
       call simulation_restart(neko_case, neko_case%chkp)

       ! Initialize the time step controller and set the time step
       call dt_controller%init(neko_case%params)
       neko_case%time%tstep = previous_save
       this%loaded_checkpoint = neko_case%time%tstep

       call profiler_start_region("Checkpoint recompute")
       ! Step through the simulation and store field states in memory
       do k = previous_save, min(next_save - 1, this%get_n_timesteps())

          ! Do not run simulation step on the first iteration
          if (k .ne. previous_save) then
             if (neko_case%time%t .ge. neko_case%time%end_time) exit
             call simulation_step(neko_case, dt_controller, loop_start)
          end if

          ! Save the restored state in memory
          local_idx = modulo(k, this%n_saves_memory) + 1
          call this%save_data(local_idx)
       end do
       call profiler_end_region("Checkpoint recompute")
    end if

    ! Restore the required time step from memory
    local_idx = modulo(tstep, this%n_saves_memory) + 1
    call this%load_data(local_idx)
  end subroutine checkpoint_restore_linear

  pure function determine_counter(tstep, n_memory, first) result(counter)
    integer, intent(in) :: tstep, n_memory, first
    integer :: counter

    if (tstep .le. first) then
       counter = tstep
    else
       counter = first + tstep / n_memory
    end if

  end function determine_counter

end submodule checkpoint_linear
