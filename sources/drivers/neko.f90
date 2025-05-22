program plainneko
  use neko, only: neko_init, neko_solve, neko_finalize
  use case, only: case_t
  use neko_ext, only: neko_top_register_types

  type(case_t) :: C

  call neko_top_register_types()
  call neko_init(C)
  call neko_solve(C)
  call neko_finalize(C)

end program plainneko
