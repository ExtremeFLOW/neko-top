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
       simulation)
    class(mma_optimizer_t), intent(inout) :: this
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(in) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), intent(in) :: simulation

    integer :: max_iterations
    real(kind=rp) :: tolerance

    real(kind=rp), dimension(:), allocatable :: x

    call json_get_or_default(parameters, "optimizer.max_iterations", &
         max_iterations, 5)
    call json_get_or_default(parameters, "optimizer.tolerance", &
         tolerance, 1.0e-3_rp)

    if (allocated(x)) deallocate(x)
    allocate(x(design%size()))
    select type (d => design)
    type is (topopt_design_t)
       call copy(x, d%design_indicator%x, d%size())
    class default
       call neko_error('Unknown design type for MMA Optimizer')
    end select

    print *, "Initializing mma_optimizer with steady_state_problem_t."

    call this%mma%init_json(x, &
         design%size(), problem%get_n_constraints(), &
         parameters, this%scale, &
         this%auto_scale)

    print *, "scale = ", this%scale
    print *, "auto_scale = ", this%auto_scale

    call this%init_from_components(problem, design, simulation, &
         max_iterations, tolerance)

    if (allocated(x)) deallocate(x)
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

    real(kind=rp), dimension(:), allocatable :: x

    integer :: iter, ierr, nglobal, n
    real(kind=rp) :: scalingfactor

    real(kind=rp) :: objective_value
    type(vector_t) :: constraint_value
    type(vector_t) :: objective_sensitivities
    type(matrix_t) :: constraint_sensitivities

    if (allocated(x)) deallocate(x)

    n = design%size()
    call MPI_Allreduce(n, nglobal, 1, MPI_INTEGER, mpi_sum, neko_comm, ierr)

    if (n .ge. 0) then
       allocate(x(n))
    end if

    select type (d => design)
    type is (topopt_design_t)
       call copy(x, d%design_indicator%x, n)
    class default
       call neko_error('Unknown design type for MMA Optimizer')
    end select

    !>initializing the scaling factor
    scalingfactor = 1.0_rp
    print *, "this%max_iterations for the optimization loop = ", this%max_iterations

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
         nglobal, this%mma%get_m(), this%tolerance

    ! Write the header for the remaining data
    write(1368, '(A4, ",", A25, ",", A25, ",", A25, ",", A25, ",", A25)') &
         "iter", "objective_value", "constraint_value", "KKTmax", &
         "KKTnorm2", "scalingfactor"

    ! Write the data row-by-row
    write(1368, '(I4, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
         & ES25.17, ",", ES25.17)') &
         0, objective_value, constraint_value%x, &
         this%mma%get_residumax(), this%mma%get_residunorm(), scalingfactor

    do iter = 1, this%max_iterations
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
       else
          write(stderr, *) "Device not supported in mma_optimizer.f90."
          error stop
       end if

       select type (d => design)
       type is (topopt_design_t)

          call copy(d%design_indicator%x, x, n)
          call d%map_forward()
          call copy(x, d%design_indicator%x, n)

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

       call this%mma%KKT(x, objective_sensitivities%x, &
            reshape([constraint_value%x], [this%mma%get_m()]), &
            constraint_sensitivities%x)

       write(1368, '(I4, ",", ES25.17, ",", ES25.17, ",", ES25.17, ",", &
            & ES25.17, ",", ES25.17)') iter, objective_value, constraint_value%x, &
            this%mma%get_residumax(), this%mma%get_residunorm(), scalingfactor

       ! Flush the buffer to write the data during the run
       flush(1368)

       call problem%write(iter)
       call simulation%reset()
    end do


    close(1368)

    ! Final state after optimization
    print*, "MMA Optimization completed after", iter-1, "iterations."

    if (allocated(x)) deallocate(x)
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

