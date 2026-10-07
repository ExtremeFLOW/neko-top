module stability_case_utils
   use case, only: case_t
   use json_file_module, only: json_file
   use profiler, only: profiler_start, profiler_stop
   use runtime_stats, only: neko_rt_stats
   use simcomp_executor, only: neko_simcomps
   use user_access_singleton, only: neko_user_access

   implicit none
   private

   public :: stability_case_init, stability_case_free

contains

   subroutine stability_case_init(base_case, parameters)
      type(case_t), intent(inout) :: base_case
      type(json_file), intent(inout) :: parameters

      call base_case%init(parameters)
      call neko_user_access%init(base_case)
      call neko_rt_stats%init(parameters)
      call neko_simcomps%init(base_case)
      call profiler_start
   end subroutine stability_case_init

   subroutine stability_case_free(base_case)
      type(case_t), intent(inout) :: base_case

      call profiler_stop
      call neko_user_access%free()
      call base_case%free()
      call neko_simcomps%free()
   end subroutine stability_case_free

end module stability_case_utils
