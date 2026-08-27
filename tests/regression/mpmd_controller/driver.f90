program mpmd_controller_driver
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use neko, only: neko_finalize, neko_init
  use neko_ctrl_mod, only: MODE_ADJOINT, MODE_FORWARD, PHASE_ADJ_RUNNING, &
       PHASE_FWD_DONE, ctrl_stream_t
  implicit none

  type(ctrl_stream_t) :: controller
  integer(int32) :: mode_cmd
  integer(int32) :: phase_cmd

  call neko_init()
  call controller%init()
  call controller%send(MODE_FORWARD, PHASE_FWD_DONE, 7_int32, 1.25_real64)

  mode_cmd = 0_int32
  phase_cmd = 0_int32
  call controller%recieve(mode_cmd, phase_cmd)

  if (mode_cmd /= MODE_ADJOINT .or. phase_cmd /= PHASE_ADJ_RUNNING) then
     error stop 'Unexpected command from Python MPMD controller peer.'
  end if

  call controller%free()
  call neko_finalize()
end program mpmd_controller_driver
