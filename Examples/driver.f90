program plainneko
   use neko
   type(case_t) :: C

   call neko_init(C)
   call neko_solve(C)
   call neko_finalize(C)


end program plainneko
