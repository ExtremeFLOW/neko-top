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


  implicit none
  private
  public :: mma_optimizer_t

  ! Concrete type for MMA optimizer
  type, extends(optimizer_t) :: mma_optimizer_t

     type(mma_t) :: mma

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
    type(topopt_design_t), intent(in) :: design

    ! Initialize MMA solver
    ! Check the type of the problem using select type
    select type (problem)
      type is (steady_state_problem_t)

       print *, "Initializing mma_optimizer with steady_state_problem_t."

       call this%mma%init_json( design%design_indicator%x, &
            design%design_indicator%size(), problem%simulation%neko_case%params, this%scale, &
            this%auto_scale)
       print *, "scale = ", this%scale
       print *, "auto_scale = ", this%auto_scale
      class default

       call neko_error('Unknown problem type in the mma_optimizer_init')
    end select

  end subroutine mma_optimizer_init

  ! Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, problem, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    real(kind=rp), intent(in) :: tolerance

    ! Check the type of the problem using select type
    select type (problem)
      type is (steady_state_problem_t)
       call this%run_ss(problem, tolerance)

      class default
       call neko_error('Unknown problem type in the mma_optimizer_run')
    end select
  end subroutine mma_optimizer_run

  subroutine mma_optimizer_run_steady_state_prob(this, problem, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    class(steady_state_problem_t), intent(inout) :: problem
    real(kind=rp), intent(in) :: tolerance
    integer :: max_iter
    integer :: iter, rank, ierr, nglobal
    real(kind=rp) :: scalingfactor

    real(kind=rp) :: objective_value
    type(vector_t) :: constraint_value
    type(vector_t) :: objective_sensitivities
    type(matrix_t) :: constraint_sensitivities

    max_iter = this%mma%get_max_iter()
    ! call MPI_Comm_rank(neko_comm, rank, ierr)
    call MPI_Allreduce(this%mma%get_n(), nglobal, 1, &
         MPI_INTEGER, mpi_sum, neko_comm, ierr)

    !>initializing the scaling factor
    scalingfactor = 1.0_rp
    print *, "max_iter for the optimization loop = ", max_iter

    call problem%compute()
    print *, "initial objective function value = " , &
         problem%volume_constraint%objective_function_value
    print *, "size(problem%design%design_indicator%x) = ", &
         size(problem%design%design_indicator%x)
    print *, "size(&
         &problem%volume_constraint%sensitivity_to_coefficient%x) = ",&
         size(&
         problem%volume_constraint%sensitivity_to_coefficient%x)


    !Writing the optimization data in a separate file
    open(1368, file = "optimization_data.txt", status = "replace")

    associate(x => problem%design%design_indicator%x, &
         f0val => &
         problem%objective_function%objective_function_value, &
         fval => &
         problem%volume_constraint%objective_function_value, &
         df0dx => &
         problem%design%sensitivity%x, &
         dfdx => &
         problem%volume_constraint%sensitivity_to_coefficient%x)

      ! Write n, m, and tolerance in the first line of optimization_data.txt
      write(1368, '("n =", I10, ", m =", I10, ", tolerance =", ES25.17)') &
           nglobal, this%mma%get_m(), tolerance

      ! Write the header for the remaining data
      write(1368, '(A)') "iter, f0val, fval(1), KKTmax, KKTnorm2, scalingfactor"

      ! Write the data row-by-row
      write(1368, '(I3, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
           & ES25.17, ",", ES25.17)') 0, f0val, fval, this%mma%get_residumax(), &
           this%mma%get_residunorm(), scalingfactor

      do iter = 1, max_iter
         if (this%mma%get_residumax() .lt. tolerance) exit
         !Scaling
         if (this%auto_scale .eqv. .true.) then
            scalingfactor = abs(this%scale/fval)
         else
            scalingfactor = abs(this%scale)
         end if

         if (NEKO_BCKND_DEVICE .eq. 0) then
            call this%mma%mma_update_cpu( iter, x, df0dx, &
                 reshape([fval*scalingfactor],[this%mma%get_m()]) , dfdx*scalingfactor)
         else
            write(stderr, *) "Device not supported in mma_optimizer.f90."
            error stop
         end if

         call problem%compute()
         call problem%compute_sensitivity()
         if (problem%design%if_mask) then
            call mask_exterior_const(&
                 problem%volume_constraint%sensitivity_to_coefficient, &
                 problem%design%optimization_domain, 0.0_rp)
         end if

         call this%mma%KKT(x, df0dx, reshape([fval], [this%mma%get_m()]), dfdx)

         print *, 'iter =', iter,&
              '-------, f0val = ', f0val, ',   fval = ', fval, &
              ',  KKTmax =', this%mma%get_residumax(), ', KKTnorm2 =',&
              this%mma%get_residunorm()

         write(1368, '(I3, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
              & ES25.17, ",", ES25.17)') iter, f0val, fval, &
              this%mma%get_residumax(), this%mma%get_residunorm(), scalingfactor
         ! Flush the buffer to write the data during the run
         flush(1368)

         call problem%sample(real(iter, rp))

         call problem%design%map_forward()
         call reset(problem%simulation%neko_case)
         ! TODO
         ! reset for the adjoint
         call field_rzero(problem%simulation%adjoint_case%scheme%u_adj)
         call field_rzero(problem%simulation%adjoint_case%scheme%v_adj)
         call field_rzero(problem%simulation%adjoint_case%scheme%w_adj)
         problem%simulation%neko_case%fluid%freeze = .false.
      end do
    end associate


    close(1368)

    ! Final state after optimization
    print*, "MMA Optimization completed after", iter-1, "iterations."
  end subroutine mma_optimizer_run_steady_state_prob

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%mma%free()
  end subroutine mma_optimizer_free

end module mma_optimizer

