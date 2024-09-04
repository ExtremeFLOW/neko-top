program usrneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use user, only: user_setup
  use adjoint_mod, only: adjoint_obj
  use solve_adjoint_mod, only: solve_adjoint
  ! use topopt, only: topopt_init, topopt_finalize, topopt_t

  type(case_t) :: C
  type(adjoint_obj) :: adj

  call user_setup(C%usr)
  call neko_init(C)
  call adj%init(C)

  call neko_solve(C)
  call solve_adjoint(adj)

  call neko_finalize(C)


end program usrneko
