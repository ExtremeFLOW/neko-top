module neko_ctrl_mod
  use, intrinsic :: iso_c_binding, only: c_int, c_double
  use comm, only: neko_comm
  use mpi_f08, only: MPI_Bcast, MPI_Comm_rank, MPI_Comm_size, MPI_INTEGER
  implicit none
  private

  ! -------------------------
  ! Public enums / constants
  ! -------------------------
  integer(c_int), public, parameter :: MODE_IDLE    = 0_c_int
  integer(c_int), public, parameter :: MODE_FORWARD = 1_c_int
  integer(c_int), public, parameter :: MODE_ADJOINT = 2_c_int
  integer(c_int), public, parameter :: MODE_STOP    = 9_c_int

  integer(c_int), public, parameter :: PHASE_INIT        = 0_c_int
  integer(c_int), public, parameter :: PHASE_FWD_RUNNING = 10_c_int
  integer(c_int), public, parameter :: PHASE_FWD_DONE    = 11_c_int
  integer(c_int), public, parameter :: PHASE_ADJ_RUNNING = 20_c_int
  integer(c_int), public, parameter :: PHASE_ADJ_DONE    = 21_c_int

  type, public :: ctrl_stream_t
     logical :: inited = .false.
     logical :: debug = .false.
     integer(c_int) :: comm_int = 0_c_int
   contains
     procedure, pass(this) :: init => ctrl_stream_init
     procedure, pass(this) :: free => ctrl_stream_free
     procedure, pass(this) :: put => ctrl_stream_put
     procedure, pass(this) :: wait_cmd => ctrl_stream_wait_cmd
  end type ctrl_stream_t

  interface
    subroutine adios2_ctrl_initialize_(comm_int) bind(C, name="adios2_ctrl_initialize_")
      import :: c_int
      integer(c_int), intent(in) :: comm_int
    end subroutine

    subroutine adios2_ctrl_finalize_() bind(C, name="adios2_ctrl_finalize_")
    end subroutine

    subroutine adios2_ctrl_put_state_(mode, phase, step, time) bind(C, name="adios2_ctrl_put_state_")
      import :: c_int, c_double
      integer(c_int), intent(in) :: mode, phase, step
      real(c_double), intent(in) :: time
    end subroutine

    subroutine adios2_ctrl_wait_cmd_(mode_cmd, phase_cmd) bind(C, name="adios2_ctrl_wait_cmd_")
      import :: c_int
      integer(c_int), intent(inout) :: mode_cmd, phase_cmd
    end subroutine
  end interface

contains

  subroutine ctrl_dbg_print(this, msg)
    class(ctrl_stream_t), intent(in) :: this
    character(len=*), intent(in) :: msg
    integer :: r, s, ierr
    if (.not. this%debug) return
    call MPI_Comm_rank(neko_comm, r, ierr)
    call MPI_Comm_size(neko_comm, s, ierr)
    write(*,'(A,I0,A,I0,A,A)') '[neko_ctrl r=', r, '/', s, '] ', trim(msg)
  end subroutine ctrl_dbg_print

  function mode_name(m) result(nm)
    integer(c_int), intent(in) :: m
    character(len=16) :: nm
    select case (m)
    case (MODE_IDLE);    nm = 'IDLE'
    case (MODE_FORWARD); nm = 'FORWARD'
    case (MODE_ADJOINT); nm = 'ADJOINT'
    case (MODE_STOP);    nm = 'STOP'
    case default;        nm = 'UNKNOWN'
    end select
  end function mode_name

  function phase_name(p) result(nm)
    integer(c_int), intent(in) :: p
    character(len=16) :: nm
    select case (p)
    case (PHASE_INIT);        nm = 'INIT'
    case (PHASE_FWD_RUNNING); nm = 'FWD_RUNNING'
    case (PHASE_FWD_DONE);    nm = 'FWD_DONE'
    case (PHASE_ADJ_RUNNING); nm = 'ADJ_RUNNING'
    case (PHASE_ADJ_DONE);    nm = 'ADJ_DONE'
    case default;             nm = 'UNKNOWN'
    end select
  end function phase_name

  subroutine ctrl_stream_init(this, comm_int)
    class(ctrl_stream_t), intent(inout) :: this
    integer(c_int), intent(in) :: comm_int
    if (this%inited) return
    this%comm_int = comm_int
    call ctrl_dbg_print(this, 'ctrl_init: calling adios2_ctrl_initialize_')
    call adios2_ctrl_initialize_(comm_int)
    call ctrl_dbg_print(this, 'ctrl_init: returned from adios2_ctrl_initialize_')
    this%inited = .true.
  end subroutine ctrl_stream_init

  subroutine ctrl_stream_free(this)
    class(ctrl_stream_t), intent(inout) :: this
    if (.not. this%inited) return
    call ctrl_dbg_print(this, 'ctrl_finalize: calling adios2_ctrl_finalize_')
    call adios2_ctrl_finalize_()
    call ctrl_dbg_print(this, 'ctrl_finalize: returned from adios2_ctrl_finalize_')
    this%inited = .false.
  end subroutine ctrl_stream_free

  subroutine ctrl_stream_put(this, mode, phase, step, time)
    class(ctrl_stream_t), intent(inout) :: this
    integer(c_int), intent(in) :: mode, phase, step
    real(c_double), intent(in) :: time
    character(len=128) :: msg
    if (.not. this%inited) return
    write(msg,'(A,A,A,A,A,I0,A,ES12.4)') 'ctrl_put: mode=', trim(mode_name(mode)), &
         ' phase=', trim(phase_name(phase)), ' step=', int(step), ' t=', real(time, kind=c_double)
    call ctrl_dbg_print(this, msg)
    call adios2_ctrl_put_state_(mode, phase, step, time)
    call ctrl_dbg_print(this, 'ctrl_put: returned from adios2_ctrl_put_state_')
  end subroutine ctrl_stream_put

  subroutine ctrl_stream_wait_cmd(this, mode_cmd, phase_cmd)
    class(ctrl_stream_t), intent(inout) :: this
    integer(c_int), intent(inout) :: mode_cmd, phase_cmd
    integer :: ierr, rank
    integer :: mode_i, phase_i
    character(len=128) :: msg
    if (.not. this%inited) return

    call MPI_Comm_rank(neko_comm, rank, ierr)

    write(msg,'(A,A,A,A)') 'ctrl_wait_cmd: enter with defaults mode=', trim(mode_name(mode_cmd)), &
         ' phase=', trim(phase_name(phase_cmd))
    call ctrl_dbg_print(this, msg)

    if (rank == 0) then
       call ctrl_dbg_print(this, 'ctrl_wait_cmd: rank0 calling C++ adios2_ctrl_wait_cmd_ (BLOCKING)')
       call adios2_ctrl_wait_cmd_(mode_cmd, phase_cmd)
       call ctrl_dbg_print(this, 'ctrl_wait_cmd: rank0 returned from C++ wait_cmd')
       mode_i  = int(mode_cmd)
       phase_i = int(phase_cmd)
    else
       call ctrl_dbg_print(this, 'ctrl_wait_cmd: non-root waiting for Bcast from rank0')
       mode_i  = 0
       phase_i = 0
    end if

    call ctrl_dbg_print(this, 'ctrl_wait_cmd: MPI_Bcast(mode)')
    call MPI_Bcast(mode_i,  1, MPI_INTEGER, 0, neko_comm, ierr)
    call ctrl_dbg_print(this, 'ctrl_wait_cmd: MPI_Bcast(phase)')
    call MPI_Bcast(phase_i, 1, MPI_INTEGER, 0, neko_comm, ierr)

    mode_cmd  = int(mode_i,  c_int)
    phase_cmd = int(phase_i, c_int)

    write(msg,'(A,A,A,A)') 'ctrl_wait_cmd: exit with mode=', trim(mode_name(mode_cmd)), &
         ' phase=', trim(phase_name(phase_cmd))
    call ctrl_dbg_print(this, msg)
  end subroutine ctrl_stream_wait_cmd

end module neko_ctrl_mod
