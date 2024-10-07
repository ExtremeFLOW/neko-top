module problem
!use design_variable, only: design_variable_t
  use topopt_design, only: topopt_design_t
  use num_types, only : rp, sp, dp

  implicit none
  private

  !> implements the problem type.
  ! Currently very abstract, could include unsteady problems etc.
  ! Also, dependingo on the type of optimizer used, we may require
  ! different functionality.
  ! Right now, all that is required in base class is to init and
  ! evaluate the problem.
  type, abstract, public :: problem_t
     !> pointer to the design
     ! TODO
     ! should be abstract!
     ! class(design_variable_t), pointer :: design
     type(topopt_design_t), pointer :: design

   contains
     !> Constructor for physics of the problem
     procedure(problem_init_base), pass(this), deferred :: init_base
     !> Additional constructor specific to a design
     procedure(problem_init_design), pass(this), deferred :: init_design
     !> Destructor.
     procedure(problem_free), pass(this), deferred :: free
     !> Evaluate the optimization problem.
     procedure(problem_compute), pass(this), deferred :: compute
     !> Sample the problem
     procedure(problem_sample), pass(this), deferred :: sample

  end type problem_t

  !> Constructor for physics of the problem
  abstract interface
     subroutine problem_init_base(this)
       import problem_t
       class(problem_t), intent(inout) :: this


     end subroutine problem_init_base
  end interface

  !> Additional constructor based on a design
  abstract interface
     subroutine problem_init_design(this, design)
       import problem_t, topopt_design_t
       class(problem_t), intent(inout) :: this
       ! class(design_variable_t), intent(in) :: design
       ! we also only have the `topopt_design_t` but this should take the more
       ! abstract `design_variable_t` and initialize differently according to
       ! the type entering here.
       type(topopt_design_t),target, intent(inout) :: design

       ! This is confusing to me..
       ! The `problem` and the `design` seem very coupled in my mind.
       ! I want to argue it's coupled one way, since the problem depends on the
       ! design representation.

       ! In principle we could have our design represented with
       ! - splines
       ! - levelset
       ! - etc
       !
       ! BUT, for density based topology optimization, because we get all our mesh
       ! info etc from neko, our design representation is based on the fluid.
       ! (of course this isn't 100% true, it's just the dofmap. We could define
       ! our design on a different set of basis functions too... but I guess that
       ! is rather far out of scope now...)
       !
       ! So it's sort of coupled both ways.. :/
       !
       ! Tim you may need to untagle this, for now I dont see an option other than
       ! - initialising the fluid first.
       !
       ! - The initializing the design
       !
       ! - Then coming here and intializing the impact of the design on the fluid
       !
     end subroutine problem_init_design
  end interface

  !> Destructor
  abstract interface
     subroutine problem_free(this)
       import problem_t
       class(problem_t), intent(inout) :: this
       ! TODO
     end subroutine problem_free
  end interface

  !> Compute
  abstract interface
     subroutine problem_compute(this)
       import problem_t
       class(problem_t), intent(inout) :: this

     end subroutine problem_compute
  end interface

  !> Sample
  abstract interface
     subroutine problem_sample(this, t)
       import problem_t, rp
       class(problem_t), intent(inout) :: this
       real(kind=rp), intent(in) :: t

     end subroutine problem_sample
  end interface
end module problem
