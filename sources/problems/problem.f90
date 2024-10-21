module problem
  use num_types, only: rp
  use topopt_design, only: topopt_design_t

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
     !  type(topopt_design_t), pointer :: design

     !> a sampler
     ! Fuck the internal samplers, they don't make sense in this context,
     ! we should have our own.
     ! - design (\rho)
     ! - mapped (\chi)
     ! - forward (u,v,w,p)
     ! - adjoint (u,v,w,p)
     ! - sensitivity to coefficients (dF/d\chi and dC/d\chi)
     ! (maybe this is redundant... but I want it for debugging)
     ! - sensitivity (dF/d\rho and dC/d\rho)
     type(fld_file_output_t) :: output

   contains

     ! ----------------------------------------------------------------------- !
     ! Interfaces

     !> Constructor for physics of the problem
     procedure(problem_init_base), pass(this), deferred :: init_base
     !> Additional constructor specific to a design
     procedure(problem_init_design), pass(this), deferred :: init_design
     !> Destructor.
     procedure(problem_free), pass(this), deferred, public :: free

     ! ----------------------------------------------------------------------- !
     ! Actual methods

     !> Sample the problem
     procedure, pass(this), public :: write => problem_write

     !> Evaluate the optimization problem.
     procedure, pass(this), public, non_overridable :: compute => &
          problem_compute

     ! ----------------------------------------------------------------------- !
     ! Dummy methods for unimplemented design types

     !> Dummy pointer, ensuring safe exit if design type not defined
     procedure, pass(this), public :: compute_topopt => dummy_compute_topopt

  end type problem_t

  abstract interface
     !> Constructor for physics of the problem
     subroutine problem_init_base(this)
       import problem_t
       class(problem_t), intent(inout) :: this
     end subroutine problem_init_base

     !> Additional constructor based on a design
     subroutine problem_init_design(this, design)
       import problem_t, topopt_design_t
       class(problem_t), intent(inout) :: this
       ! class(design_variable_t), intent(in) :: design
       ! we also only have the `topopt_design_t` but this should take the more
       ! abstract `design_variable_t` and initialize differently according to
       ! the type entering here.
       type(topopt_design_t), target, intent(inout) :: design

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

     !> Destructor
     subroutine problem_free(this)
       import problem_t
       class(problem_t), intent(inout) :: this
       ! TODO
     end subroutine problem_free

     !> Compute
     subroutine problem_compute_topopt(this)
       import problem_t
       import topopt_design_t
       class(problem_t), intent(inout) :: this
     end subroutine problem_compute_topopt
  end interface

contains

  !> Compute
  subroutine problem_compute(this, design)
    class(problem_t), intent(inout) :: this
    class(design_t), intent(in) :: design

    select type(design)
      type is(topopt_design_t)
       call this%compute_topopt(design)
      class default
       call neko_error("Design type not supported.")
    end select

  end subroutine problem_compute

  !> Sample the fields/design.
  subroutine problem_write(this, t)
    class(problem_t), intent(inout) :: this
    real(kind=rp), intent(in) :: t

    call this%output%sample(t)
  end subroutine problem_write

  !> Dummy compute function
  subroutine dummy_compute_topopt(this, design)
    class(problem_t), intent(inout) :: this
    type(topopt_design_t), intent(in) :: design
    call neko_error("problem_compute_topopt not implemented.")
  end subroutine dummy_compute_topopt
end module problem
