module optimizer
  !-----------------------------------------------------------!
  ! An abstract type optimizer is defined to solve the design !
  ! optimization problem using a specific type of optimizer   !
  ! algorithm, e.g. MMA.                                      !
  !-----------------------------------------------------------!

  use problem, only: problem_t
  use topopt_design, only: topopt_design_t
  use num_types, only : rp
  implicit none
  private

  !> Abstract optimizer class.
  type, abstract, public :: optimizer_t

   contains
     !> Initialize the optimizer, associate it with a specific problem
     procedure(optimizer_init), pass(this), deferred :: init
     !> Run the optimization loop
     procedure(optimizer_run), pass(this), deferred :: run
     !> Free resources.
     procedure(optimizer_free), pass(this), deferred :: free
  end type optimizer_t

  !> Interface for optimizer initialization
  abstract interface
     subroutine optimizer_init(this, problem, design)
       import optimizer_t, problem_t, topopt_design_t
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(inout) :: problem
       type(topopt_design_t), intent(inout) :: design
     end subroutine optimizer_init
  end interface

  !> Interface for running the optimization loop
  abstract interface
     subroutine optimizer_run(this, problem, design, tolerance)
       import optimizer_t, problem_t, rp, topopt_design_t
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(inout) :: problem
       type(topopt_design_t), intent(inout) :: design
       real(kind=rp), intent(in) :: tolerance
     end subroutine optimizer_run
  end interface

  !> Interface for freeing resources
  abstract interface
     subroutine optimizer_free(this)
       import optimizer_t
       class(optimizer_t), intent(inout) :: this
     end subroutine optimizer_free
  end interface

end module optimizer
