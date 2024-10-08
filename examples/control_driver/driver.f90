program usrneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use user, only: user_setup
  use adjoint_case, only: adjoint_case_t, adjoint_init, adjoint_free
  use simulation_adjoint, only: solve_adjoint
  ! use topopt, only: topopt_init, topopt_finalize, topopt_t

  type(case_t) :: C
  type(adjoint_case_t) :: adj

  call user_setup(C%usr)
  call neko_init(C)

  call adjoint_init(adj, C)

  call neko_solve(C)
  call solve_adjoint(adj)

  call adjoint_free(adj)
  call neko_finalize(C)
end program usrneko
