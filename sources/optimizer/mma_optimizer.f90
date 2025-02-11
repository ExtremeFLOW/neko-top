module mma_optimizer
  use optimizer, only: optimizer_t
  use steady_state_problem, only : steady_state_problem_t
  use mma, only: mma_t
  use problem, only: problem_t
  use topopt_design, only: topopt_design_t
  use num_types, only : rp
  use utils, only : neko_error
  use json_module, only: json_file
  use simulation, only: simulation_t
  use design, only: design_t

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

     !> Maximum number of iterations
     integer :: max_iterations
     !> Tolerance for the optimization loop
     real(kind=rp) :: tolerance

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
  subroutine mma_optimizer_init(this, parameters, simulation, problem, design)
    class(mma_optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    type(simulation_t), intent(in) :: simulation
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design

    type(topopt_design_t), pointer :: top_des
    this%tolerance = 1.0e-3_rp

    top_des => null()
    select type (design)
    type is (topopt_design_t)
       top_des => design
    class default
       call neko_error('Unknown design type in the mma_optimizer_init')
    end select

    ! Initialize MMA solver
    ! Check the type of the problem using select type
    select type (problem)
    type is (steady_state_problem_t)

       print *, "Initializing mma_optimizer with steady_state_problem_t."

       call this%mma%init_json(top_des%design_indicator%x, &
            design%size(), problem%get_n_constraints(), &
            parameters, this%scale, &
            this%auto_scale)

       print *, "scale = ", this%scale
       print *, "auto_scale = ", this%auto_scale
    class default

       call neko_error('Unknown problem type in the mma_optimizer_init')
    end select

  end subroutine mma_optimizer_init

  ! Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, simulation, problem, design)
    class(mma_optimizer_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design

    ! Check the type of the problem using select type
    select type (problem)
    type is (steady_state_problem_t)

       select type (design)
       type is (topopt_design_t)
          call this%run_ss(simulation, problem, design)
       class default
          call neko_error('Unknown design type in the mma_optimizer_run')
       end select

    class default
       call neko_error('Unknown problem type in the mma_optimizer_run')
    end select
  end subroutine mma_optimizer_run

  subroutine mma_optimizer_run_steady_state_prob(this, simulation, problem, design)
    class(mma_optimizer_t), intent(inout) :: this
    type(simulation_t), intent(inout) :: simulation
    class(steady_state_problem_t), intent(inout) :: problem
    class(topopt_design_t), intent(inout) :: design
    integer :: max_iter
    integer :: iter, ierr, nglobal
    real(kind=rp) :: scalingfactor

    real(kind=rp) :: objective_value
    type(vector_t) :: constraint_value
    type(vector_t) :: objective_sensitivities
    type(matrix_t) :: constraint_sensitivities

    ! call MPI_Comm_rank(neko_comm, rank, ierr)
    call MPI_Allreduce(this%mma%get_n(), nglobal, 1, &
         MPI_INTEGER, mpi_sum, neko_comm, ierr)

    !>initializing the scaling factor
    scalingfactor = 1.0_rp
    max_iter = 5
    print *, "max_iter for the optimization loop = ", max_iter

    call simulation%run_forward()
    call problem%compute(design)

    call simulation%run_backward()
    call problem%compute_sensitivity(design)

    call problem%get_objective_value(objective_value)
    call problem%get_constraint_values(constraint_value)
    call problem%get_objective_sensitivities(objective_sensitivities)
    call problem%get_constraint_sensitivities(constraint_sensitivities)

    !Writing the optimization data in a separate file
    open(1368, file = "optimization_data.txt", status = "replace")

    associate(x => design%design_indicator%x)

      ! Write n, m, and tolerance in the first line of optimization_data.txt
      write(1368, '("n =", I10, ", m =", I10, ", tolerance =", ES25.17)') &
           nglobal, this%mma%get_m(), this%tolerance

      ! Write the header for the remaining data
      write(1368, '(A)') "iter, objective_value, constraint_value, KKTmax, &
           &KKTnorm2, scalingfactor"

      ! Write the data row-by-row
      write(1368, '(I3, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
           & ES25.17, ",", ES25.17)') &
           0, objective_value, constraint_value%x, &
           this%mma%get_residumax(), this%mma%get_residunorm(), scalingfactor

      call problem%write(0)

      do iter = 1, max_iter
         if (this%mma%get_residumax() .lt. this%tolerance) exit

         !Scaling
         if (this%auto_scale .eqv. .true.) then
            scalingfactor = abs(this%scale/constraint_value%x(1))
         else
            scalingfactor = abs(this%scale)
         end if

         if (NEKO_BCKND_DEVICE .eq. 0) then
            call this%mma%mma_update_cpu(iter, x, objective_sensitivities%x, &
                 constraint_value%x * scalingfactor, &
                 constraint_sensitivities%x*scalingfactor)
            call design%map_forward()
         else
            write(stderr, *) "Device not supported in mma_optimizer.f90."
            error stop
         end if

         call simulation%run_forward()
         call problem%compute(design)

         call simulation%run_backward()
         call problem%compute_sensitivity(design)

         call problem%get_objective_value(objective_value)
         call problem%get_constraint_values(constraint_value)
         call problem%get_objective_sensitivities(objective_sensitivities)
         call problem%get_constraint_sensitivities(constraint_sensitivities)

         call this%mma%KKT(x, objective_sensitivities%x, &
              reshape([constraint_value%x], [this%mma%get_m()]), &
              constraint_sensitivities%x)

         write(1368, '(I3, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
              & ES25.17, ",", ES25.17)') iter, objective_value, constraint_value%x, &
              this%mma%get_residumax(), this%mma%get_residunorm(), scalingfactor

         ! Flush the buffer to write the data during the run
         flush(1368)

         call problem%write(iter)

         call simulation%reset()
      end do
    end associate


    close(1368)

    ! Final state after optimization
    print*, "MMA Optimization completed after", iter-1, "iterations."


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

