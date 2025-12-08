module optimizer
  !-----------------------------------------------------------!
  ! An abstract type optimizer is defined to solve the design !
  ! optimization problem using a specific type of optimizer   !
  ! algorithm, e.g. MMA.                                      !
  !-----------------------------------------------------------!

  use json_module, only: json_file
  use simulation_m, only: simulation_t
  use problem, only: problem_t
  use design, only: design_t
  use num_types, only: rp
  use logger, only: neko_log
  use profiler, only: profiler_start_region, profiler_end_region

  implicit none
  private

  !> Abstract optimizer class.
  type, abstract :: optimizer_t
     private

     !> The maximum number of iterations
     integer, public :: max_iterations
     !> The tolerance for the optimization loop
     real(kind=rp), public :: tolerance

   contains
     !> Initialize the optimizer, associate it with a specific problem
     procedure(optimizer_init_from_json), pass(this), public, deferred :: &
          init_from_json
     !> Run the optimization loop
     procedure, pass(this), public :: run => optimizer_run
     !> Perform a single optimization step
     procedure(optimizer_step), pass(this), public, deferred :: step
     !> Free resources.
     procedure(optimizer_free), pass(this), public, deferred :: free
     !> Free base resources.
     procedure, pass(this) :: free_base => optimizer_free_base

     !> Validate the solution
     procedure(optimizer_validate), pass(this), public, deferred :: validate
     !> Write the progress of the optimizer to the log file
     procedure(optimizer_write), pass(this), public, deferred :: write

     !> The base initializer
     procedure, pass(this) :: init_base => optimizer_init_base


  end type optimizer_t

  ! -------------------------------------------------------------------------- !
  ! Interface for the optimizer module.

  abstract interface
     !> Interface for optimizer initialization
     subroutine optimizer_init_from_json(this, parameters, problem, design, &
          simulation)
       import optimizer_t, json_file, simulation_t, problem_t, design_t, rp
       class(optimizer_t), intent(inout) :: this
       type(json_file), intent(inout) :: parameters
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(in) :: design
       type(simulation_t), optional, intent(in) :: simulation
     end subroutine optimizer_init_from_json

     !> Interface for running an optimization step
     logical function optimizer_step(this, iter, problem, design, simulation)
       import optimizer_t, simulation_t, problem_t, design_t
       integer, intent(in) :: iter
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(inout) :: design
       type(simulation_t), optional, intent(inout) :: simulation
     end function optimizer_step

     !> Interface for writing the optimizer progress
     subroutine optimizer_write(this, iter, problem)
       import optimizer_t, simulation_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       integer, intent(in) :: iter
       class(problem_t), intent(in) :: problem
     end subroutine optimizer_write

     !> Interface for freeing resources
     subroutine optimizer_free(this)
       import optimizer_t
       class(optimizer_t), intent(inout) :: this
     end subroutine optimizer_free

     !> Interface for validating the solution
     subroutine optimizer_validate(this, problem, design)
       import optimizer_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       class(problem_t), intent(in) :: problem
       class(design_t), intent(in) :: design
     end subroutine optimizer_validate
  end interface

  ! -------------------------------------------------------------------------- !
  ! Interfaces for the factory functions

  !> Factory function for the optimizer
  !! @param object The optimizer object to be created.
  !! @param parameters The JSON file containing the optimizer parameters.
  !! @param problem The problem object.
  !! @param design The design object.
  !! @param simulation The simulation object.
  interface optimizer_factory
     module subroutine optimizer_factory(object, parameters, problem, design, &
          simulation)
       class(optimizer_t), allocatable, intent(inout) :: object
       type(json_file), intent(inout) :: parameters
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(in) :: design
       type(simulation_t), optional, intent(in) :: simulation
     end subroutine optimizer_factory
  end interface optimizer_factory

  public :: optimizer_t, optimizer_factory

contains

  !> Base initializer for the optimizer
  !! @param this The optimizer object.
  !! @param max_iterations The maximum number of iterations.
  !! @param tolerance The tolerance for the optimization loop.
  subroutine optimizer_init_base(this, max_iterations, tolerance)
    class(optimizer_t), intent(inout) :: this
    integer, intent(in) :: max_iterations
    real(kind=rp), intent(in) :: tolerance

    this%max_iterations = max_iterations
    this%tolerance = tolerance

  end subroutine optimizer_init_base

  !> Define the optimization loop
  subroutine optimizer_run(this, problem, design, simulation)
    class(optimizer_t), intent(inout) :: this
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(inout) :: design
    type(simulation_t), optional, intent(inout) :: simulation
    logical :: converged
    character(len=256) :: msg
    integer :: iter

    call neko_log%section('Optimization Loop')
    converged = .false.
    do iter = 1, this%max_iterations

       call profiler_start_region('Optimizer iteration')
       converged = this%step(iter, problem, design, simulation)
       call profiler_end_region('Optimizer iteration')

       ! Log the progress and outputs
       call this%write(iter, problem)
       call design%write(iter)

       if (converged) exit
    end do
    call neko_log%end_section()

    call this%validate(problem, design)

    if (.not. converged) then
       write(msg, '(A,I0,A)') 'Optimizer did not converge in ', &
            this%max_iterations, ' iterations.'
       call neko_log%warning(trim(msg))
    else
       write(msg, '(A,I0,A)') 'Optimizer converged after ', iter, &
            ' iterations.'
       call neko_log%message(trim(msg))
    end if

  end subroutine optimizer_run

  !> Base free routine for the optimizer
  !! @param this The optimizer object.
  subroutine optimizer_free_base(this)
    class(optimizer_t), intent(inout) :: this
  end subroutine optimizer_free_base

end module optimizer
