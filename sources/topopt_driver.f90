program topneko
   use neko
   use user
   use neko_top
   type(case_t) :: C

   call user_setup(C%usr)
   call neko_init(C)

   call top_setup(C)

   call neko_solve(C)
   call neko_finalize(C)


end program topneko

