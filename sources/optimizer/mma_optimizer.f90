module mma_optimizer
  use optimizer, only: optimizer_t
  use steady_state_problem, only : steady_state_problem_t
  use mma, only: mma_t
  use problem, only: problem_t
  use topopt_design, only: topopt_design_t
  use num_types, only : rp
  use utils, only : neko_error
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
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
     generic :: init => init_from_json, init_from_components
     procedure, pass(this) :: init_from_json => mma_optimizer_init_from_json
     procedure, pass(this) :: init_from_components => &
          mma_optimizer_init_from_components

     procedure :: run => mma_optimizer_run
     procedure :: free => mma_optimizer_free

  end type mma_optimizer_t

contains

  !> Initialize the MMA optimizer from JSON file
  subroutine mma_optimizer_init_from_json(this, parameters, problem, design, &
       simulation, max_iterations, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), intent(in) :: simulation
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance
    type(vector_t) :: x

    call x%init(design%size())

    select type (d => design)
    type is (topopt_design_t)
       call copy(x%x, d%design_indicator%x, d%size())
    class default
       call neko_error('Unknown design type for MMA Optimizer')
    end select

    print *, "Initializing mma_optimizer with steady_state_problem_t."

    call this%mma%init_json(x%x, &
         design%size(), problem%get_n_constraints(), &
         parameters, this%scale, &
         this%auto_scale)

    call this%init_from_components(problem, design, simulation, &
         max_iterations, tolerance)

    call x%free()
  end subroutine mma_optimizer_init_from_json

  !> Initialize the MMA optimizer from JSON file
  subroutine mma_optimizer_init_from_components(this, problem, design, &
       simulation, max_iterations, tolerance)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), intent(in) :: simulation
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance

    call this%init_base(max_iterations, tolerance)

  end subroutine mma_optimizer_init_from_components

  !> Define the optimization loop for MMA
  subroutine mma_optimizer_run(this, problem, design, simulation)
    class(mma_optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), intent(inout) :: simulation

    type(vector_t) :: x

    integer :: iter, ierr, nglobal, n
    real(kind=rp) :: scaling_factor

    real(kind=rp) :: objective_value
    type(vector_t) :: constraint_value
    type(vector_t) :: objective_sensitivities
    type(matrix_t) :: constraint_sensitivities

    n = design%size()
    call MPI_Allreduce(n, nglobal, 1, MPI_INTEGER, mpi_sum, neko_comm, ierr)

    call x%init(n)

    select type (d => design)
    type is (topopt_design_t)
       call copy(x%x, d%design_indicator%x, n)
    class default
       call neko_error('Unknown design type for MMA Optimizer')
    end select

    !>initializing the scaling factor
    scaling_factor = 1.0_rp
    print *, "max_iterations for the optimization loop = ", this%max_iterations

    call simulation%run_forward()
    call problem%compute(design)

    call simulation%run_backward()
    call problem%compute_sensitivity(design)

    call problem%get_objective_value(objective_value)
    call problem%get_constraint_values(constraint_value)
    call problem%get_objective_sensitivities(objective_sensitivities)
    call problem%get_constraint_sensitivities(constraint_sensitivities)

    call problem%write(0)

    !Writing the optimization data in a separate file
    open(1368, file = "optimization_data.txt", status = "replace")

    ! Write n, m, and tolerance in the first line of optimization_data.txt
    write(1368, '("n =", I10, ", m =", I10, ", tolerance =", ES25.17)') &
         nglobal, problem%get_n_constraints(), this%tolerance

    ! Write the header for the remaining data
    write(1368, '(A4, ",", A25, ",", A25, ",", A25, ",", A25, ",", A25)') &
         "iter", "objective_value", "constraint_value", "KKTmax", &
         "KKTnorm2", "scalingfactor"

    ! Write the data row-by-row
    write(1368, '(I4, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
         & ES25.17, ",", ES25.17)') &
         0, objective_value, constraint_value%x, &
         this%mma%get_residumax(), this%mma%get_residunorm(), scaling_factor

    do iter = 1, this%max_iterations
       if (this%mma%get_residumax() .lt. this%tolerance) exit

       !Scaling
       if (this%auto_scale .eqv. .true.) then
          scaling_factor = abs(this%scale/constraint_value%x(1))
       else
          scaling_factor = abs(this%scale)
       end if

       ! Scale the constraint value and sensitivities
       constraint_value = scaling_factor * constraint_value
       constraint_sensitivities = scaling_factor * constraint_sensitivities

       call this%mma%update(iter, x, objective_sensitivities, &
            constraint_value, constraint_sensitivities)

       select type (d => design)
       type is (topopt_design_t)

          call copy(d%design_indicator%x, x%x, n)
          call d%map_forward()
          call copy(x%x, d%design_indicator%x, n)

       class default
          call neko_error('Unknown design type for MMA Optimizer')
       end select

       call simulation%run_forward()
       call problem%compute(design)

       call simulation%run_backward()
       call problem%compute_sensitivity(design)

       call problem%get_objective_value(objective_value)
       call problem%get_constraint_values(constraint_value)
       call problem%get_objective_sensitivities(objective_sensitivities)
       call problem%get_constraint_sensitivities(constraint_sensitivities)

       call this%mma%KKT(x%x, objective_sensitivities%x, &
            constraint_value%x, constraint_sensitivities%x)

       write(1368, '(I4, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
            & ES25.17, ",", ES25.17)') iter, objective_value, constraint_value%x, &
            this%mma%get_residumax(), this%mma%get_residunorm(), scaling_factor

       ! Flush the buffer to write the data during the run
       flush(1368)

       call problem%write(iter)
       call simulation%reset()
    end do


    close(1368)

    ! Final state after optimization
    print*, "MMA Optimization completed after", iter-1, "iterations."

    call constraint_value%free()
    call objective_sensitivities%free()
    call constraint_sensitivities%free()

  end subroutine mma_optimizer_run

  ! Free resources associated with the MMA optimizer
  subroutine mma_optimizer_free(this)
    class(mma_optimizer_t), intent(inout) :: this

    ! Free MMA-specific data
    call this%mma%free()
  end subroutine mma_optimizer_free

end module mma_optimizer

