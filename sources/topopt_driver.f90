! ============================================================================ !
! Neko Top program
! ============================================================================ !

program nekotop
   use neko_top, only: case_t, topology_t
   use neko_top, only: neko_top_init, neko_top_solve, neko_top_finalize

   implicit none

   type(case_t) :: neko_case
   type(topology_t) :: topology

   call neko_top_init(neko_case, topology)
   call neko_top_solve(neko_case, topology)
   call neko_top_finalize(neko_case, topology)

end program nekotop
