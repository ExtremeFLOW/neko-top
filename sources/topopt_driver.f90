! ============================================================================ !
! Neko Top program
! ============================================================================ !

program nekotop
  use neko_top
  use case, only: case_t

  implicit none

  type(case_t) :: neko_case

  call neko_top_init(neko_case)
  call neko_top_solve(neko_case)
  call neko_top_finalize(neko_case)

end program nekotop
