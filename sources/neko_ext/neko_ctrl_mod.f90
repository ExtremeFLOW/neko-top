!> @file neko_ctrl_mod.f90
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
!!   * Redistributions in binary form must reproduce the above copyright
!!     notice, this list of conditions and the following disclaimer in the
!!     documentation and/or other materials provided with the distribution.
!!
!!   * Neither the name of the authors nor the names of its contributors may
!!     be used to endorse or promote products derived from this software
!!     without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
!> @brief Control stream helpers for POD in-situ coordination.
module neko_ctrl_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use comm, only: neko_comm
  use utils, only: neko_error
  use mpi_f08, only: MPI_Bcast, MPI_Comm_rank, MPI_Comm_size, &
       MPI_INTEGER4, MPI_COMM_WORLD, MPI_Send, MPI_Recv, &
       MPI_REAL8, MPI_STATUS_IGNORE
  implicit none
  private

  ! Control protocol:
  ! - mode describes the coarse solver stage, while phase describes the
  !   synchronization point inside that stage.
  ! - MODE_IDLE + PHASE_INIT is a one-off startup tick used to prove that the
  !   control channel is alive. The Python driver may ignore it.
  ! - MODE_FORWARD + PHASE_FWD_RUNNING means Neko is in the forward solve and
  !   is currently streaming snapshots to Python.
  ! - MODE_FORWARD + PHASE_FWD_DONE means the forward window is complete. Neko
  !   is now blocked at the forward-to-adjoint boundary while Python finalizes
  !   the POD basis, writes outputs, and prepares the adjoint data.
  ! - MODE_ADJOINT + PHASE_ADJ_RUNNING is Python's reply that POD modes and
  !   time coefficients are ready, so Neko may enter adjoint
  !   restore/reconstruction.
  ! - MODE_ADJOINT + PHASE_ADJ_DONE means the adjoint window is complete and
  !   the Python side should discard the old POD state and wait for a new
  !   forward window.
  ! - MODE_STOP is a terminal shutdown message. The current implementation
  !   pairs it with PHASE_ADJ_DONE as a placeholder phase, but the phase is not
  !   semantically important once MODE_STOP is seen.
  integer(int32), public, parameter :: MODE_IDLE = 0_int32
  integer(int32), public, parameter :: MODE_FORWARD = 1_int32
  integer(int32), public, parameter :: MODE_ADJOINT = 2_int32
  integer(int32), public, parameter :: MODE_STOP = 9_int32

  integer(int32), public, parameter :: PHASE_INIT = 0_int32
  integer(int32), public, parameter :: PHASE_FWD_RUNNING = 10_int32
  integer(int32), public, parameter :: PHASE_FWD_DONE = 11_int32
  integer(int32), public, parameter :: PHASE_ADJ_RUNNING = 20_int32
  integer(int32), public, parameter :: PHASE_ADJ_DONE = 21_int32

  integer, parameter :: CTRL_TAG_STATE_INT = 4101
  integer, parameter :: CTRL_TAG_STATE_REAL = 4102
  integer, parameter :: CTRL_TAG_CMD = 4103

  !> Lightweight control channel for coordinating POD streaming phases.
  type, public :: ctrl_stream_t
     logical :: inited = .false.
     logical :: debug = .false.
     integer :: peer_root = -1
   contains
     procedure, pass(this) :: init => ctrl_stream_init
     procedure, pass(this) :: free => ctrl_stream_free
     procedure, pass(this) :: send => ctrl_stream_send
     procedure, pass(this) :: recieve => ctrl_stream_recieve
  end type ctrl_stream_t

contains

  !> Debug print for control stream (rank-tagged).
  subroutine ctrl_dbg_print(this, msg)
    class(ctrl_stream_t), intent(in) :: this
    character(len=*), intent(in) :: msg
    integer :: r, s, ierr

    if (.not. this%debug) return
    call MPI_Comm_rank(neko_comm, r, ierr)
    call MPI_Comm_size(neko_comm, s, ierr)
    write(*, '(A,I0,A,I0,A,A)') '[neko_ctrl r=', r, '/', s, '] ', trim(msg)
  end subroutine ctrl_dbg_print

  !> Convert mode enum to human-readable name.
  function mode_name(m) result(nm)
    integer(int32), intent(in) :: m
    character(len=16) :: nm

    select case (m)
    case (MODE_IDLE)
       nm = 'IDLE'
    case (MODE_FORWARD)
       nm = 'FORWARD'
    case (MODE_ADJOINT)
       nm = 'ADJOINT'
    case (MODE_STOP)
       nm = 'STOP'
    case default
       nm = 'UNKNOWN'
    end select
  end function mode_name

  !> Convert phase enum to human-readable name.
  function phase_name(p) result(nm)
    integer(int32), intent(in) :: p
    character(len=16) :: nm

    select case (p)
    case (PHASE_INIT)
       nm = 'INIT'
    case (PHASE_FWD_RUNNING)
       nm = 'FWD_RUNNING'
    case (PHASE_FWD_DONE)
       nm = 'FWD_DONE'
    case (PHASE_ADJ_RUNNING)
       nm = 'ADJ_RUNNING'
    case (PHASE_ADJ_DONE)
       nm = 'ADJ_DONE'
    case default
       nm = 'UNKNOWN'
    end select
  end function phase_name

  !> Read the peer root rank from the environment.
  subroutine ctrl_stream_init_peer_root(this)
    class(ctrl_stream_t), intent(inout) :: this
    character(len=32) :: env_val
    integer :: env_len
    integer :: ios

    call get_environment_variable("NEKO_CTRL_PEER_ROOT", env_val, env_len)
    if (env_len <= 0) then
       call neko_error('NEKO_CTRL_PEER_ROOT must be set for POD MPI control.')
    end if

    read(env_val(1:env_len), *, iostat=ios) this%peer_root
    if (ios /= 0 .or. this%peer_root < 0) then
       call neko_error('Invalid NEKO_CTRL_PEER_ROOT for POD MPI control.')
    end if
  end subroutine ctrl_stream_init_peer_root

  !> Initialize the MPI control stream.
  subroutine ctrl_stream_init(this)
    class(ctrl_stream_t), intent(inout) :: this

    if (this%inited) return
    call ctrl_stream_init_peer_root(this)
    call ctrl_dbg_print(this, 'ctrl_init: MPI control ready')
    this%inited = .true.
  end subroutine ctrl_stream_init

  !> Finalize the MPI control stream.
  subroutine ctrl_stream_free(this)
    class(ctrl_stream_t), intent(inout) :: this

    if (.not. this%inited) return
    call ctrl_dbg_print(this, 'ctrl_finalize: MPI control done')
    this%inited = .false.
    this%peer_root = -1
  end subroutine ctrl_stream_free

  !> Send current mode/phase/step/time over the control stream.
  subroutine ctrl_stream_send(this, mode, phase, step, time)
    class(ctrl_stream_t), intent(inout) :: this
    integer(int32), intent(in) :: mode, phase, step
    real(real64), intent(in) :: time
    integer(int32) :: state_i(3)
    character(len=128) :: msg
    integer :: ierr
    integer :: rank

    if (.not. this%inited) return

    call MPI_Comm_rank(neko_comm, rank, ierr)
    if (rank /= 0) return

    state_i(1) = int(mode, int32)
    state_i(2) = int(phase, int32)
    state_i(3) = int(step, int32)

    write(msg, '(A,A,A,A,A,I0,A,ES12.4)') 'ctrl_send: mode=', &
         trim(mode_name(mode)), ' phase=', trim(phase_name(phase)), &
         ' step=', int(step), ' t=', real(time, kind=real64)
    call ctrl_dbg_print(this, msg)
    call MPI_Send(state_i, size(state_i), MPI_INTEGER4, this%peer_root, &
         CTRL_TAG_STATE_INT, MPI_COMM_WORLD, ierr)
    call MPI_Send(time, 1, MPI_REAL8, this%peer_root, &
         CTRL_TAG_STATE_REAL, MPI_COMM_WORLD, ierr)
  end subroutine ctrl_stream_send

  !> Recieve a control command and broadcast it to all ranks.
  subroutine ctrl_stream_recieve(this, mode_cmd, phase_cmd)
    class(ctrl_stream_t), intent(inout) :: this
    integer(int32), intent(inout) :: mode_cmd, phase_cmd
    integer :: ierr, rank
    integer(int32) :: mode_i, phase_i
    integer(int32) :: cmd_i(2)
    character(len=128) :: msg

    if (.not. this%inited) return

    call MPI_Comm_rank(neko_comm, rank, ierr)

    write(msg, '(A,A,A,A)') 'ctrl_recieve: enter with defaults mode=', &
         trim(mode_name(mode_cmd)), ' phase=', trim(phase_name(phase_cmd))
    call ctrl_dbg_print(this, msg)

    if (rank == 0) then
       call ctrl_dbg_print(this, 'ctrl_recieve: rank0 waiting on MPI cmd')
       call MPI_Recv(cmd_i, size(cmd_i), MPI_INTEGER4, this%peer_root, &
            CTRL_TAG_CMD, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
       mode_i = cmd_i(1)
       phase_i = cmd_i(2)
    else
       call ctrl_dbg_print(this, &
            'ctrl_recieve: non-root waiting for Bcast from rank0')
       mode_i = 0_int32
       phase_i = 0_int32
    end if

    call ctrl_dbg_print(this, 'ctrl_recieve: MPI_Bcast(mode)')
    call MPI_Bcast(mode_i, 1, MPI_INTEGER4, 0, neko_comm, ierr)
    call ctrl_dbg_print(this, 'ctrl_recieve: MPI_Bcast(phase)')
    call MPI_Bcast(phase_i, 1, MPI_INTEGER4, 0, neko_comm, ierr)

    mode_cmd = mode_i
    phase_cmd = phase_i

    write(msg, '(A,A,A,A)') 'ctrl_recieve: exit with mode=', &
         trim(mode_name(mode_cmd)), ' phase=', trim(phase_name(phase_cmd))
    call ctrl_dbg_print(this, msg)
  end subroutine ctrl_stream_recieve

end module neko_ctrl_mod
