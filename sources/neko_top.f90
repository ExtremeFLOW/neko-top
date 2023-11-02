module neko_top
   use neko
   use json_module
   implicit none

contains

   subroutine top_setup(C)
      type(case_t), target, intent(inout) :: C
      logical :: value

      call json_get(C%params, 'case.topopt.test_mode', value)

      if (value) then
         print *, "topopt test mode"
      else
         print *, "topopt normal mode"
      end if

   end subroutine top_setup

end module neko_top

