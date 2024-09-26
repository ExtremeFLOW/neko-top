module problem
  use case, only: case_t

  implicit none
  private

  type, abstract :: problem_t
     !> In the future I would imagine the problem would own the fluid and adjoint and passive scalars etc
     ! depending on what we're doing. For now, I'm keeping them outside and hardcoding a specific example
     
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
