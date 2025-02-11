module mma_optimizer
  use optimizer, only: optimizer_t
  use steady_state_problem, only : steady_state_problem_t
  use mma, only: mma_t
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

     type(mma_t) :: mma
     type(topopt_design_t), pointer :: design

     !> Scaling constraint_value%x and constraint_sensitivities%x.
     !! Note that the values are not updated but they are scaled when passed
     !! to the optimizer.
     !! (if auto_scale then constraint_value%x=scale else constraint_value%x=scale*constraint_value%x)
     !! When auto_scale is true, we use an adaptable scale for
     !! constraint_value%x and constraint_sensitivities%x in every iteration (variable scale factors)
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
    type(topopt_design_t), target, intent(in) :: design

    this%design => design

    ! Initialize MMA solver
    ! Check the type of the problem using select type
    select type (problem)
    type is (steady_state_problem_t)

       print *, "Initializing mma_optimizer with steady_state_problem_t."

       call this%mma%init_json(design%design_indicator%x, &
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
  subroutine mma_optimizer_run(this, problem, tolerance, max_iter)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: max_iter

    ! Check the type of the problem using select type
    select type (problem)
    type is (steady_state_problem_t)
       call this%run_ss(problem, tolerance, max_iter)

    class default
       call neko_error('Unknown problem type in the mma_optimizer_run')
    end select
  end subroutine mma_optimizer_run

  subroutine mma_optimizer_run_steady_state_prob(this, problem, tolerance, &
       max_iter)
    class(mma_optimizer_t), intent(inout) :: this
    class(steady_state_problem_t), intent(inout) :: problem
    real(kind=rp), intent(in) :: tolerance
    integer, intent(in) :: max_iter
    integer :: iter, ierr, nglobal
    real(kind=rp) :: scalingfactor

    real(kind=rp) :: objective_value
    type(vector_t) :: constraint_value
    type(vector_t) :: objective_sensitivities
    type(matrix_t) :: constraint_sensitivities

    ! These would all be owned by the logger_t
    ! ------------------------------------------------------------------------
    !> format for writing real values in the log file
    character(len=10) :: real_format 
    ! logging
    character(len=1024) :: problem_log ! I assume 1024 is enough?
    character(len=1024) :: optimization_log ! I assume 1024 is enough?
    character(len=100) :: log_format 
    integer :: log_unit
    ! ------------------------------------------------------------------------


    ! call MPI_Comm_rank(neko_comm, rank, ierr)
    call MPI_Allreduce(this%mma%get_n(), nglobal, 1, &
         MPI_INTEGER, mpi_sum, neko_comm, ierr)

    !>initializing the scaling factor
    scalingfactor = 1.0_rp
    print *, "max_iter for the optimization loop = ", max_iter

    call problem%compute(this%design)
    call problem%compute_sensitivity(this%design)

    call problem%get_objective_value(objective_value)
    call problem%get_constraint_values(constraint_value)
    call problem%get_objective_sensitivities(objective_sensitivities)
    call problem%get_constraint_sensitivities(constraint_sensitivities)

    ! Should I write a logger_t such that this becomes a logger%init()
    ! ------------------------------------------------------------------------
    !Writing the optimization data in a separate file
    log_unit = 1368
    real_format = 'ES25.17' 
    open(log_unit, file = "optimization_data.txt", status = "replace")
    log_format = '(I3, A, ", ", '//real_format//', ", ", '//real_format//&
       ', ", ", '//real_format//')'
    ! Write n, m, and tolerance in the first line of optimization_data.txt
    ! do we really want this?
    if (pe_rank .eq. 0) then
       write(log_unit, '("n =", I10, ", m =", I10, ", tolerance =", ES25.17)') &
            nglobal, this%mma%get_m(), tolerance
    end if

    ! Write the header for the remaining data
    problem_log = problem%get_log_header()
    optimization_log = 'iter, '//trim(problem_log)//&
       ', KKTmax, KKTnorm2, scalaing factor'
    if (pe_rank .eq. 0) then
       write(log_unit, '(A)') trim(optimization_log)
    end if
    ! ------------------------------------------------------------------------

    associate(x => this%design%design_indicator%x)


      ! Should I write a logger_t such that this becomes a logger%write()
      ! ----------------------------------------------------------------------
      ! Write the data row-by-row
      problem_log = problem%get_log_state(real_format)
      write(optimization_log, log_format) 0, trim(problem_log), &
         this%mma%get_residumax(), this%mma%get_residunorm(), scalingfactor
      if (pe_rank .eq. 0) then
         write(1368, '(A)') trim(optimization_log)
      end if
      ! ----------------------------------------------------------------------

      call problem%write(0)

      do iter = 1, max_iter
         if (this%mma%get_residumax() .lt. tolerance) exit

         !Scaling
         if (this%auto_scale .eqv. .true.) then
            scalingfactor = abs(this%scale/constraint_value%x(1))
         else
            scalingfactor = abs(this%scale)
         end if

         if (NEKO_BCKND_DEVICE .eq. 0) then
            call this%mma%mma_update_cpu( iter, x, objective_sensitivities%x, &
                 constraint_value%x * scalingfactor, &
                 constraint_sensitivities%x*scalingfactor)
         else
            write(stderr, *) "Device not supported in mma_optimizer.f90."
            error stop
         end if

         call this%design%map_forward()
         call problem%compute(this%design)
         call problem%compute_sensitivity(this%design)

         call problem%get_objective_value(objective_value)
         call problem%get_constraint_values(constraint_value)
         call problem%get_objective_sensitivities(objective_sensitivities)
         call problem%get_constraint_sensitivities(constraint_sensitivities)

         call this%mma%KKT(x, objective_sensitivities%x, &
              reshape([constraint_value%x], [this%mma%get_m()]), &
              constraint_sensitivities%x)

         ! Write the data row-by-row
         ! should I implement a logger_t, such that this is a logger%write()
         ! -------------------------------------------------------------------
         problem_log = problem%get_log_state(real_format)
         write(optimization_log, log_format) iter, trim(problem_log), &
            this%mma%get_residumax(), this%mma%get_residunorm(), scalingfactor
         if (pe_rank .eq. 0) then
            write(log_unit, '(A)') trim(optimization_log)
         end if
         ! Flush the buffer to write the data during the run
         flush(log_unit)
         ! -------------------------------------------------------------------

         call problem%write(iter)

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
    print*, "MMA Optimization completed after", iter-1, "iterations."

    close(log_unit)
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

