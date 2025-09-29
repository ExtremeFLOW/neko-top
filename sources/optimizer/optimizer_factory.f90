

submodule(optimizer) optimizer_factory_mod
  use json_utils, only: json_get, json_get_or_default
  use utils, only: neko_type_error
  use dummy_constraint, only: dummy_constraint_t
  use mma_optimizer, only: mma_optimizer_t
  use constraint, only: constraint_t

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
    class(simulation_t), optional, intent(in) :: simulation

    character(len=:), allocatable :: type
    integer :: max_iterations
    real(kind=rp) :: tolerance
    logical :: performance
    class(constraint_t), allocatable :: dummy_con

    if (allocated(object)) then
       call object%free()
       deallocate(object)
    end if

    ! Get the type of the optimizer
    call json_get(parameters, "optimization.solver.type", type)
    call json_get_or_default(parameters, "optimization.solver.max_iterations", &
         max_iterations, 100)
    call json_get_or_default(parameters, "optimization.solver.tolerance", &
         tolerance, 1.0e-3_rp)
    call json_get_or_default(parameters, "optimization.solver.performance", &
         performance, .false.)

    ! Select the optimizer type
    select case (trim(type))
    case ("mma")
       allocate(mma_optimizer_t::object)

    case default
       call neko_type_error("Optimizer", type, KNOWN_TYPES)
    end select

    !Check if we are solving an unconstrained problem and add a dummy contraint
    if (problem%get_n_constraints() .eq. 0) then
       allocate(dummy_constraint_t::dummy_con)
       call dummy_con%init(parameters, design)
       call problem%add_constraint(dummy_con)
    end if

    if (present(simulation)) then
       call object%init_from_json(parameters, problem, design, &
            max_iterations, tolerance, performance, simulation)
    else
       call object%init_from_json(parameters, problem, design, &
            max_iterations, tolerance, performance)
    end if

  end subroutine optimizer_factory


end submodule optimizer_factory_mod
