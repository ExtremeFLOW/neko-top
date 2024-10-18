module mma_optimizer
  use optimizer, only: optimizer_t
  use steady_state_problem, only : steady_state_problem_t
  use mma, only: mma_t
  use problem, only: problem_t
  use num_types, only : rp
  use utils, only : neko_error

  implicit none
  private
  public :: mma_optimizer_t

  ! Concrete type for MMA optimizer
  type, extends(optimizer_t) :: mma_optimizer_t

      type(mma_t) :: mma

      !Pointer to the specific problem_type(steady_state_problem_t in this case)
      type(steady_state_problem_t), pointer :: steady_state_prob => null()
      ! type(unsteady_problem_t), pointer :: unsteady_prob => null()

  contains
      ! Override the deferred methods
      procedure :: init => mma_optimizer_init
      procedure :: run => mma_optimizer_run
      procedure :: free => mma_optimizer_free
  end type mma_optimizer_t

contains

  !> Initialize the MMA optimizer, associate it with a specific problem
  subroutine mma_optimizer_init(this, prob)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), target, intent(inout) :: prob
    real(kind=rp), allocatable :: xmax(:), xmin(:)


    ! Associate the problem with the optimizer
    this%prob => prob

    !setting xmax, xmin using neko_scratch_registry
    allocate(xmax(prob%design%design_indicator%size()))
    allocate(xmin(prob%design%design_indicator%size()))
    xmin=0_rp
    xmax=1_rp

    ! Initialize MMA solver
    ! Check the type of the problem using select type
    select type(prob)
    type is (steady_state_problem_t)
      ! Now we know prob is of type steady_state_problem_t, assign the pointer
      this%steady_state_prob => prob
      print *, "Initializing mma_optimizer with steady_state_problem_t."
      ! mma_init_json(this, x, n, m, a0, a, c, d, xmin, xmax, json)
      call this%mma%init_json( prob%design%design_indicator%x, &
        prob%design%design_indicator%size(), 1, 1.0_rp, [0.0_rp], [100.0_rp], &
        [0.0_rp], xmin, xmax, prob%C%params)


    class default
      !Unknown problem
      call neko_error('Unknown problem type in the mma_optimizer_init')
    end select
  end subroutine mma_optimizer_init

  ! Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, tolerance, max_iter)
    class(mma_optimizer_t), intent(inout) :: this
    integer, intent(in) :: max_iter
    real(kind=rp), intent(in) :: tolerance
    integer :: iter
    


    !check if there is a drived type for prob
    !if not, then prob is of abstract type problem_t and therefore we get errors
    !later on we will add other types of problem here as well:
    ! if (.not. (associated(this%steady_state_prob) .and. &
    !           (associated(this%unsteady_prob) .and.
    !     ......  ))
    if (.not. associated(this%steady_state_prob)) then
      call neko_error('steady_state_prob not initialized in mma_optimizer_run')
    endif

    call this%prob%compute()
    ! print *, "initial objective function value=" , &
    !   this%steady_state_prob%volume_constraint%objective_function_value
    ! print *, "size(this%prob%design%design_indicator%x)=", &
    !   size(this%steady_state_prob%design%design_indicator%x)
    ! print *, "size(this%prob%volume_constraint%sensitivity_to_coefficient%x)=",&
    !   size(this%steady_state_prob%volume_constraint%sensitivity_to_coefficient%x)
          
    associate(x => this%steady_state_prob%design%design_indicator%x, &
      f0val => this%steady_state_prob%objective_function%dissipation, &
      fval => this%steady_state_prob%volume_constraint%objective_function_value, &
      df0dx => this%steady_state_prob%design%sensitivity%x, &
      dfdx => this%steady_state_prob%volume_constraint%sensitivity_to_coefficient%x)
    do iter = 1, 2 !max_iter
      if (this%mma%get_residumax() .lt. tolerance) exit

      ! mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
      call this%mma%mma_update_cpu( iter, x, df0dx, &
        reshape([fval],[this%mma%get_m()]) , dfdx)
      call this%prob%compute()
      
      call this%mma%KKT(x,df0dx,reshape([fval],[this%mma%get_m()]),dfdx)
      print *, 'iter=', iter,&
        '-------,f0val= ', f0val, ',   fval= ', fval, &
        ',  KKTmax=', this%mma%get_residumax(), ', KKTnorm2=', this%mma%get_residunorm()

    end do
    end associate

    ! mma_update_cpu(this, iter, x, df0dx, fval, dfdx)
    !call this%mma%update(this%prob%design%design_indicator%x)

    ! Get the initial objective value
    ! call this%prob%compute()
    ! print *, "initial objective function value=" , &
    !     this%prob%volume_constraint%objective_function_value

    ! print *, "size(this%prob%design%design_indicator%x)=", &
    !     size(this%prob%design%design_indicator%x)

    ! print *, "size(this%prob%volume_constraint%sensitivity_to_coefficient%x)=",&
    !     size(this%prob%volume_constraint%sensitivity_to_coefficient%x)

    ! ! Optimization loop
    ! do iter = 1, max_iter
    !     ! Call the update method (MMA-specific update)
    !     call this%update()

    !     ! Recompute the objective after the design update
    !     call this%prob%compute()
    !     obj_current = this%prob%objective_value

    !     ! Check for convergence
    !     if (abs(obj_current - obj_prev) < tolerance) exit

    !     obj_prev = obj_current
    ! end do

    ! Final state after optimization
    print*, 'MMA Optimization completed after', iter, 'iterations.'

  end subroutine mma_optimizer_run

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%mma%free()
  end subroutine mma_optimizer_free

end module mma_optimizer

