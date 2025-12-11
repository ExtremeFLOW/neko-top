program usrneko
  use neko
  use user       ! provided by insitu_turb_pipe.f90
  implicit none

  type(case_t), target :: C

  call user_setup(C%user)
  call neko_init(C)
  call neko_user_access%init(C)
  call neko_solve(C)
  call neko_finalize(C)
end program usrneko
