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
  use field, only: field_t
  use field_math, only: field_copy
  use scalar_scheme, only: scalar_scheme_t

contains

  !> Save the current state of the simulation in RAM only.
  module subroutine checkpoint_save_ram(this, neko_case, time)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    integer :: tstep
    integer :: i_scalars, idx
    type(field_t), pointer :: u, v, w, p, s

    tstep = time%tstep
    if (tstep .lt. 1) then
       call neko_error("Requested timestep is out of range for RAM checkpoint")
    end if

    call checkpoint_ensure_ram_capacity(this, neko_case, tstep)

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    s => null()

    call field_copy(this%p_list(tstep), p)
    call field_copy(this%u_list(tstep), u)
    call field_copy(this%v_list(tstep), v)
    call field_copy(this%w_list(tstep), w)

    do i_scalars = 1, this%n_scalars
       idx = (tstep - 1) * this%n_scalars + i_scalars
       s => neko_case%scalars%scalar_fields(i_scalars)%s
       call field_copy(this%s_list(idx), s)
    end do
  end subroutine checkpoint_save_ram

  !> Restore the forward simulation state from RAM only.
  module subroutine checkpoint_restore_ram(this, neko_case, time)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    type(time_state_t), intent(in) :: time
    integer :: tstep
    integer :: i_scalars, idx
    type(field_t), pointer :: u, v, w, p, s

    tstep = time%tstep
    if (tstep .lt. 1 .or. tstep .gt. this%get_n_timesteps()) then
       call neko_error("Requested timestep is out of range for RAM checkpoint")
    end if
    if (.not. allocated(this%p_list)) then
       call neko_error("Requested RAM checkpoint is not available")
    else if (tstep .gt. size(this%p_list)) then
       call neko_error("Requested RAM checkpoint is not available")
    end if

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    s => null()

    call field_copy(p, this%p_list(tstep))
    call field_copy(u, this%u_list(tstep))
    call field_copy(v, this%v_list(tstep))
    call field_copy(w, this%w_list(tstep))
    do i_scalars = 1, this%n_scalars
       idx = (tstep - 1) * this%n_scalars + i_scalars
       s => neko_case%scalars%scalar_fields(i_scalars)%s
       call field_copy(s, this%s_list(idx))
    end do
  end subroutine checkpoint_restore_ram

  subroutine checkpoint_ensure_ram_capacity(this, neko_case, required_steps)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case
    integer, intent(in) :: required_steps
    integer :: current_capacity, new_capacity
    integer :: i, j, idx
    character(len=80) :: str
    type(field_t), allocatable :: p_new(:), u_new(:), v_new(:), w_new(:)
    type(field_t), allocatable :: s_new(:)

    if (allocated(this%p_list)) then
       current_capacity = size(this%p_list)
    else
       current_capacity = 0
    end if

    if (required_steps .le. current_capacity) return

    if (current_capacity .gt. 0) then
       new_capacity = max(required_steps, current_capacity + &
            max(1, this%n_saves_memory))
    else
       new_capacity = max(required_steps, max(1, this%n_saves_memory))
    end if

    allocate(p_new(new_capacity))
    allocate(u_new(new_capacity))
    allocate(v_new(new_capacity))
    allocate(w_new(new_capacity))
    if (this%n_scalars .gt. 0) then
       allocate(s_new(new_capacity * this%n_scalars))
    end if

    do i = 1, new_capacity
       write(str, '(A,I0)') "p_chkp_", i
       call p_new(i)%init(neko_case%fluid%p%dof, str)
       write(str, '(A,I0)') "u_chkp_", i
       call u_new(i)%init(neko_case%fluid%u%dof, str)
       write(str, '(A,I0)') "v_chkp_", i
       call v_new(i)%init(neko_case%fluid%v%dof, str)
       write(str, '(A,I0)') "w_chkp_", i
       call w_new(i)%init(neko_case%fluid%w%dof, str)
       if (this%n_scalars .gt. 0) then
          do j = 1, this%n_scalars
             write(str, '(A,I0,A,I0)') "s_chkp_", i, "_", j
             call s_new((i - 1) * this%n_scalars + j)%init( &
                  neko_case%scalars%scalar_fields(j)%s%dof, str)
          end do
       end if
    end do

    if (current_capacity .gt. 0) then
       do i = 1, current_capacity
          call field_copy(p_new(i), this%p_list(i))
          call field_copy(u_new(i), this%u_list(i))
          call field_copy(v_new(i), this%v_list(i))
          call field_copy(w_new(i), this%w_list(i))
          if (this%n_scalars .gt. 0) then
             do j = 1, this%n_scalars
                idx = (i - 1) * this%n_scalars + j
                call field_copy(s_new(idx), this%s_list(idx))
             end do
          end if
       end do
    end if

    if (allocated(this%p_list)) then
       do i = 1, size(this%p_list)
          call this%p_list(i)%free()
          call this%u_list(i)%free()
          call this%v_list(i)%free()
          call this%w_list(i)%free()
       end do
    end if
    if (allocated(this%s_list)) then
       do i = 1, size(this%s_list)
          call this%s_list(i)%free()
       end do
    end if

    if (allocated(this%p_list)) deallocate(this%p_list)
    if (allocated(this%u_list)) deallocate(this%u_list)
    if (allocated(this%v_list)) deallocate(this%v_list)
    if (allocated(this%w_list)) deallocate(this%w_list)
    if (allocated(this%s_list)) deallocate(this%s_list)

    call move_alloc(p_new, this%p_list)
    call move_alloc(u_new, this%u_list)
    call move_alloc(v_new, this%v_list)
    call move_alloc(w_new, this%w_list)
    if (allocated(s_new)) then
       call move_alloc(s_new, this%s_list)
    end if
  end subroutine checkpoint_ensure_ram_capacity

end submodule checkpoint_ram
