program plainneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t

  type(case_t) :: C

  call neko_init(C)
  call neko_solve(C)
  call neko_finalize(C)

end program plainneko
