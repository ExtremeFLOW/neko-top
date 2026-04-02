!> @file checkpoint_ram.f90
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

!> Implementation for the RAM Checkpointing algorithm.
!! In this case, we save every time step in memory and never write to disk.
submodule (simulation_checkpoint) checkpoint_ram

contains

  !> Save the current state of the simulation in RAM only.
  !! @param[inout] this Checkpointing implementation.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] time Current time state.
  module subroutine checkpoint_save_ram(this, neko_case, time)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    integer :: tstep

    tstep = time%tstep
    if (tstep .lt. 1) then
       call neko_error("Requested timestep is out of range for RAM checkpoint")
    end if

    call checkpoint_ensure_ram_capacity(this, tstep)
    call this%save_data(tstep)
  end subroutine checkpoint_save_ram

  !> Restore the forward simulation state from RAM only.
  !! @param[inout] this Checkpointing implementation.
  !! @param[inout] neko_case Case data structure.
  !! @param[in] time Target time state.
  module subroutine checkpoint_restore_ram(this, neko_case, time)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    integer :: tstep

    tstep = time%tstep
    if (tstep .lt. 1 .or. tstep .gt. this%get_n_timesteps()) then
       call neko_error("Requested timestep is out of range for RAM checkpoint")
    end if

    if (.not. allocated(this%state_storage)) then
       call neko_error("Requested RAM checkpoint is not available")
    else if (tstep .gt. size(this%state_storage, 1)) then
       call neko_error("Requested RAM checkpoint is not available")
    end if

    call this%load_data(tstep)
  end subroutine checkpoint_restore_ram

  !> Grow the RAM checkpoint storage to hold @a required_steps snapshots.
  subroutine checkpoint_ensure_ram_capacity(this, required_steps)
    class(simulation_checkpoint_t), intent(inout) :: this
    integer, intent(in) :: required_steps
    type(host_array), allocatable :: storage_new(:, :)
    integer :: current_capacity, new_capacity
    integer :: i, j, n_fields

    if (allocated(this%state_storage)) then
      current_capacity = size(this%state_storage, 1)
      n_fields = size(this%state_storage, 2)
    else
      current_capacity = 0
      n_fields = this%state_list%size()
    end if

    if (required_steps .le. current_capacity) return

    new_capacity = max(required_steps, max(1, max(this%n_saves_memory, &
         2 * current_capacity)))
    allocate(storage_new(new_capacity, n_fields))

    if (allocated(this%state_storage)) then
       do i = 1, current_capacity
          do j = 1, n_fields
             if (this%state_storage(i, j)%is_allocated()) then
                call storage_new(i, j)%init(this%state_storage(i, j)%size)
                call copy(storage_new(i, j)%data, &
                     this%state_storage(i, j)%data, &
                     this%state_storage(i, j)%size)
                call this%state_storage(i, j)%free()
             end if
          end do
       end do
       deallocate(this%state_storage)
    end if

    call move_alloc(storage_new, this%state_storage)
    this%n_saves_memory = new_capacity
  end subroutine checkpoint_ensure_ram_capacity

end submodule checkpoint_ram
