module topopt
  use case, only: case_t

  implicit none
  private

  type, public :: topopt_t
     !> The case
     type(case_t), pointer :: neko_case
   contains
     !> Constructor
     procedure, pass(this) :: init => topopt_init
     !> Destructor
     procedure, pass(this) :: free => topopt_free
  end type topopt_t

contains

  !> Constructor
  subroutine topopt_init(this, C)
    class(topopt_t), intent(inout) :: this
    type(case_t), pointer :: C
    this%neko_case => C
  end subroutine topopt_init

  !> Destructor
  subroutine topopt_free(this)
    class(topopt_t), intent(inout) :: this
    this%neko_case => null()
  end subroutine topopt_free

end module topopt
