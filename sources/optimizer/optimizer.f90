module optimizer
  !-----------------------------------------------------------!
  ! An abstract type optimizer is defined to solve the design !
  ! optimization problem using a specific type of optimizer   !
  ! algorithm, e.g. MMA.                                      !
  !-----------------------------------------------------------!

  use json_module, only: json_file
  use simulation, only: simulation_t
  use problem, only: problem_t
  use design, only: design_t
  use num_types, only : rp
  implicit none
  private

  !> Abstract optimizer class.
  type, abstract, public :: optimizer_t
     private


   contains
     !> Initialize the optimizer, associate it with a specific problem
     procedure(optimizer_init), pass(this), public, deferred :: init
     !> Run the optimization loop
     procedure(optimizer_run), pass(this), public, deferred :: run
     !> Free resources.
     procedure(optimizer_free), pass(this), public, deferred :: free
  end type optimizer_t

  ! -------------------------------------------------------------------------- !
  ! Interface for the optimizer module.

  abstract interface
     !> Interface for optimizer initialization
     subroutine optimizer_init(this, parameters, simulation, problem, design)
       import optimizer_t, json_file, simulation_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       type(json_file), intent(inout) :: parameters
       type(simulation_t), intent(in) :: simulation
       class(problem_t), intent(in) :: problem
       class(design_t), intent(in) :: design
     end subroutine optimizer_init

     !> Interface for running the optimization loop
     subroutine optimizer_run(this, simulation, problem, design)
       import optimizer_t, simulation_t, problem_t, design_t
       class(optimizer_t), intent(inout) :: this
       type(simulation_t), intent(inout) :: simulation
       class(problem_t), intent(inout) :: problem
       class(design_t), intent(inout) :: design
     end subroutine optimizer_run

     !> Interface for freeing resources
     subroutine optimizer_free(this)
       import optimizer_t
       class(optimizer_t), intent(inout) :: this
     end subroutine optimizer_free
  end interface

end module optimizer
