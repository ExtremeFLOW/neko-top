module mma_optimizer
  use optimizer, only: optimizer_t
  use steady_state_problem, only : steady_state_problem_t
  use mma, only: mma_t, mma_factory
  use problem, only: problem_t
  use topopt_design, only: topopt_design_t
  use num_types, only : rp
  use utils, only : neko_error

  use vector, only: vector_t
  use matrix, only: matrix_t

  !only to print nglobal when running in parallel
  use comm, only: neko_comm
  use mpi_f08, only: MPI_INTEGER, mpi_sum, MPI_Allreduce


  use neko_config, only: NEKO_BCKND_DEVICE
  ! Inclusions from external dependencies and standard libraries
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit

  use math, only: copy, cmult
  use field_math, only: field_rzero
  use neko_ext, only: reset
  use mask_ops, only: mask_exterior_const

  use comm, only : pe_rank


  implicit none
  private
  public :: mma_optimizer_t

  ! Concrete type for MMA optimizer
  type, extends(optimizer_t) :: mma_optimizer_t

      ! type(mma_t) :: mma
      class(mma_t), allocatable :: mma
      !> Scaling fval and dfdx.
      !! Note that the values are not updated but they are scaled when passed
      !! to the optimizer.
      !! (if auto_scale then fval=scale else fval=scale*fval)
      !! When auto_scale is true, we use an adaptable scale for
      !! fval and dfdx in every iteration (variable scale factors)
      real(kind=rp) :: scale
      logical :: auto_scale
  contains
      ! Override the deferred methods
      procedure :: init => mma_optimizer_init
      procedure :: run => mma_optimizer_run
      procedure :: free => mma_optimizer_free

     procedure, pass(this) :: run_ss => mma_optimizer_run_steady_state_prob
  end type mma_optimizer_t

contains

  !> Initialize the MMA optimizer, associate it with a specific problem
  subroutine mma_optimizer_init(this, problem, design)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    type(topopt_design_t), intent(inout) :: design

    ! Initialize MMA solver
    ! Check the type of the problem using select type
    select type (problem)
    type is (steady_state_problem_t)

       print *, "Initializing mma_optimizer with steady_state_problem_t."

       call this%mma%init_json(this%mma, design%design_indicator%x, &
            problem%get_n_design(), problem%get_n_constraints(), &
            problem%simulation%neko_case%params, this%scale, &
            this%auto_scale)

       print *, "scale = ", this%scale
       print *, "auto_scale = ", this%auto_scale
    class default

       call neko_error('Unknown problem type in the mma_optimizer_init')
    end select

  end subroutine mma_optimizer_init

  ! Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, problem, design, tolerance, max_iter)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    type(topopt_design_t), intent(inout) :: design
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: max_iter

    ! Check the type of the problem using select type
    select type (problem)
      type is (steady_state_problem_t)
       call this%run_ss(problem, design, tolerance, max_iter)

    class default
       call neko_error('Unknown problem type in the mma_optimizer_run')
    end select
  end subroutine mma_optimizer_run

  subroutine mma_optimizer_run_steady_state_prob(this, problem, design, &
       tolerance, max_iter)
    class(mma_optimizer_t), intent(inout) :: this
    class(steady_state_problem_t), intent(inout) :: problem
    type(topopt_design_t), intent(inout) :: design
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: max_iter
    integer :: iter, ierr, nglobal
    real(kind=rp) :: scalingfactor

    real(kind=rp) :: objective_value
    type(vector_t) :: constraint_value
    type(vector_t) :: objective_sensitivities
    type(matrix_t) :: constraint_sensitivities

    call MPI_Allreduce(this%mma%get_n(), nglobal, 1, &
         MPI_INTEGER, mpi_sum, neko_comm, ierr)

    !>initializing the scaling factor
    scalingfactor = 1.0_rp
    print *, "max_iter for the optimization loop = ", max_iter

    call problem%compute(design)
    call problem%compute_sensitivity(design)

    call problem%get_objective_value(objective_value)
    call problem%get_constraint_values(constraint_value)
    call problem%get_objective_sensitivities(objective_sensitivities)
    call problem%get_constraint_sensitivities(constraint_sensitivities)

    associate(x => design%design_indicator%x)

      do iter = 1, max_iter
         if (this%mma%get_residumax() .lt. tolerance) exit

         !Scaling
         if (this%auto_scale .eqv. .true.) then
            scalingfactor = abs(this%scale/constraint_value%x(1))
         else
            scalingfactor = abs(this%scale)
         end if

         call this%mma%update(iter, x, objective_sensitivities, &
              scalingfactor*constraint_value, &
              scalingfactor*constraint_sensitivities)

         call problem%compute(design)
         call problem%compute_sensitivity(design)

         call problem%get_objective_value(objective_value)
         call problem%get_constraint_values(constraint_value)
         call problem%get_objective_sensitivities(objective_sensitivities)
         call problem%get_constraint_sensitivities(constraint_sensitivities)

         call this%mma%KKT(x, objective_sensitivities, &
              constraint_value, constraint_sensitivities)

         call problem%write(iter)

         call design%map_forward()
         call reset(problem%simulation%neko_case)
         ! TODO
         ! reset for the adjoint
         call field_rzero(problem%simulation%adjoint_case%scheme%u_adj)
         call field_rzero(problem%simulation%adjoint_case%scheme%v_adj)
         call field_rzero(problem%simulation%adjoint_case%scheme%w_adj)
         problem%simulation%neko_case%fluid%freeze = .false.
      end do
    end associate


    ! Final state after optimization
    if (pe_rank .eq. 0) then
       print*, "MMA Optimization completed after", iter-1, "iterations."
    end if

    call constraint_value%free()
    call objective_sensitivities%free()
    call constraint_sensitivities%free()

  end subroutine mma_optimizer_run_steady_state_prob

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%mma%free()
  end subroutine mma_optimizer_free

end module mma_optimizer

