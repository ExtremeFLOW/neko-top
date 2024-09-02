program usrneko
  use neko, only: neko_solve, neko_init, neko_finalize
  use case, only: case_t
  use user, only: user_setup
  type(case_t), target :: C

  call user_setup(C%usr)
  call neko_init(C)
  call neko_solve(C)
  call neko_finalize(C)


end program usrneko
