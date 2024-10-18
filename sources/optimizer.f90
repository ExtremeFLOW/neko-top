module optimizer
    !-----------------------------------------------------------!
    ! An abstract type optimizer is defined to solve the design !
    ! optimization problem using a specific type of optimizer   !
    ! algorithm, e.g. MMA.                                      !
    !-----------------------------------------------------------!

    use problem, only: problem_t
    use num_types, only : rp

    implicit none
    private

    !> Abstract optimizer class.
    type, abstract, public :: optimizer_t
        !> Pointer to the problem that the optimizer will work on
        class(problem_t), pointer :: prob

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
        subroutine optimizer_init(this, prob)
            import optimizer_t, problem_t
            class(optimizer_t), intent(inout) :: this
            class(problem_t), target, intent(inout) :: prob
        end subroutine optimizer_init
    end interface

    !> Interface for running the optimization loop
    abstract interface
        subroutine optimizer_run(this,tolerance,max_iter)
            import optimizer_t, rp
            class(optimizer_t), intent(inout) :: this
            real(kind=rp), intent(in) :: tolerance
            integer, intent(in) :: max_iter
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
