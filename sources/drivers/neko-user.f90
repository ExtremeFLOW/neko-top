program usrneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use user, only: user_setup
  use neko_ext, only: neko_top_register_types

  type(case_t), target :: C

  call neko_top_register_types()
  call user_setup(C%usr)
  call neko_init(C)
  call neko_solve(C)
  call neko_finalize(C)


end program usrneko
