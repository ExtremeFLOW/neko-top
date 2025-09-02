! Copyright (c) 2025, The Neko Authors
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
!
submodule (simulation_checkpoint) checkpoint_linear
  use simulation, only: simulation_step, simulation_restart
  use file, only: file_t, file_free
  use time_step_controller, only: time_step_controller_t

contains

  !> Save the current state of the simulation in a linear fashion
  module subroutine checkpoint_save_linear(this, neko_case)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), intent(inout) :: neko_case

    ! Sample the checkpoint if needed
    if (modulo(neko_case%time%tstep, this%n_saves_memory) .eq. 0 .or. &
         neko_case%time%tstep .le. this%first_valid_timestep) then

       call this%chkp_output%set_counter(neko_case%time%tstep)
       call this%chkp_output%sample(neko_case%time%t)
       this%n_saves_disc = this%n_saves_disc + 1
    end if
  end subroutine checkpoint_save_linear

  !> Restore the forward simulation state in a linear fashion
  module subroutine checkpoint_restore_linear(this, neko_case, tstep)
    class(simulation_checkpoint_t), intent(inout) :: this
    class(case_t), target, intent(inout) :: neko_case
    integer, intent(in) :: tstep
    type(time_step_controller_t) :: dt_controller
    type(file_t) :: chkp_file
    character(len=256) :: chkp_file_name
    real(kind=dp) :: loop_start
    integer :: j, k, previous_save, next_save
    integer :: i_scalars, n_scalars
    type(field_t), pointer :: u, v, w, p, s

    loop_start = MPI_WTIME()

    u => neko_case%fluid%u
    v => neko_case%fluid%v
    w => neko_case%fluid%w
    p => neko_case%fluid%p
    s => null()

    ! Determine the nearest save states
    previous_save = tstep - modulo(tstep, this%n_saves_memory)
    next_save = previous_save + this%n_saves_memory

    if (previous_save .le. this%first_valid_timestep) then
       previous_save = min(tstep, this%first_valid_timestep)
    end if

    if (previous_save .lt. this%first_valid_timestep) then
       next_save = previous_save + 1
    end if

    if (this%loaded_checkpoint .ne. previous_save) then
       write(chkp_file_name, '(A,I5.5,A)') trim(this%filename), &
            previous_save, '.chkp'

       call chkp_file%init(chkp_file_name)
       call chkp_file%read(neko_case%chkp)
       call file_free(chkp_file)

       call dt_controller%init(neko_case%params)
       call simulation_restart(neko_case, neko_case%chkp)
       neko_case%time%tstep = previous_save
       this%loaded_checkpoint = neko_case%time%tstep

       do k = previous_save, min(next_save - 1, this%n_timesteps)

          if (k .ne. previous_save) then
             if (neko_case%time%t .ge. neko_case%time%end_time) exit
             call simulation_step(neko_case, dt_controller, loop_start)
          end if

          j = modulo(k, this%n_saves_memory) + 1
          call field_copy(this%p_list(j), p)
          call field_copy(this%u_list(j), u)
          call field_copy(this%v_list(j), v)
          call field_copy(this%w_list(j), w)
          if (this%have_scalar) then
             n_scalars = size(neko_case%scalars%scalar_fields)
             do i_scalars = 1, n_scalars
                s => neko_case%scalars%scalar_fields(i_scalars)%s
                call field_copy(this%s_list((j - 1) * n_scalars + i_scalars), s)
             end do
          end if
       end do

    end if

    j = modulo(tstep, this%n_saves_memory) + 1

    call field_copy(p, this%p_list(j))
    call field_copy(u, this%u_list(j))
    call field_copy(v, this%v_list(j))
    call field_copy(w, this%w_list(j))
    if (this%have_scalar) then
       n_scalars = size(neko_case%scalars%scalar_fields)
       do i_scalars = 1, n_scalars
          s => neko_case%scalars%scalar_fields(i_scalars)%s
          call field_copy(s, this%s_list((j - 1) * n_scalars + i_scalars))
       end do
    end if

  end subroutine checkpoint_restore_linear

end submodule checkpoint_linear
