module system_base
   use LightKrylov, only: abstract_system_rdp
   use case, only: case_t

   implicit none
   private

   type, abstract, extends(abstract_system_rdp), public :: stability_system_t
      type(case_t), pointer, public :: base_case => null()
      logical :: if_2d = .false.
      character(len=:), allocatable :: base_output_base_fname
   contains
      private
      procedure, pass(self), public :: init => stability_system_init
      procedure, pass(self), public :: free => stability_system_free
   end type stability_system_t

   public :: stability_system_init, stability_system_free

contains

   subroutine stability_system_init(self, base_case)
      class(stability_system_t), intent(inout) :: self
      type(case_t), target, intent(inout) :: base_case

      self%base_case => base_case
      self%base_output_base_fname = trim(base_case%f_out%file_%get_base_fname())
      self%if_2d = (self%base_case%fluid%c_Xh%msh%gdim .eq. 2)
   end subroutine stability_system_init

   subroutine stability_system_free(self)
      class(stability_system_t), intent(inout) :: self

      nullify(self%base_case)
      self%if_2d = .false.
      if (allocated(self%base_output_base_fname)) then
         deallocate(self%base_output_base_fname)
      end if
   end subroutine stability_system_free

end module system_base
