

submodule(optimizer) optimizer_factory_mod
  use json_utils, only: json_get
  use utils, only: neko_type_error
  use mma_optimizer, only: mma_optimizer_t


  implicit none

  !> Known function types
  character(len=25), parameter :: KNOWN_TYPES(1) = [ character(len=25) :: &
       "mma"]

contains

  ! -------------------------------------------------------------------------- !
  ! Interfaces for the factory functions

  !> @brief Factory function for the optimizer.
  !! @details This function creates an optimizer object based on the type
  !! specified in the JSON file under "optimizer.type".
  !!
  !! @param object The optimizer object to be created.
  !! @param parameters The JSON file containing the optimizer parameters.
  !! @param problem The problem object.
  !! @param design The design object.
  !! @param simulation The simulation object.
  module subroutine optimizer_factory(object, parameters, problem, design, &
       simulation)
    class(optimizer_t), allocatable, intent(inout) :: object
    type(json_file), intent(inout) :: parameters
    class(problem_t), intent(inout) :: problem
    class(design_t), intent(in) :: design
    type(simulation_t), optional, intent(in) :: simulation

    character(len=:), allocatable :: type

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    ! Get the type of the optimizer
    call json_get(parameters, "optimization.solver.type", type)

    ! Select the optimizer type
    select case (trim(type))
    case ("mma")
       allocate(mma_optimizer_t::object)
    case default
       call neko_type_error("Optimizer", type, KNOWN_TYPES)
    end select

    call object%init_from_json(parameters, problem, design, simulation)

  end subroutine optimizer_factory


end submodule optimizer_factory_mod
