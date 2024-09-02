program usrneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use user, only: user_setup
  ! use topopt, only: topopt_init, topopt_finalize, topopt_t

  type(case_t) :: C

  call user_setup(C%usr)
  call neko_init(C)

  call neko_solve(C)

  call neko_finalize(C)


end program usrneko
