! User module for the user defined simulation component
module user
  use user_intf, only: user_t
  use json_module, only: json_file
  use mma_simcomp, only: mma_comp_t
  use simulation_component, only: simulation_component_t, &
       simulation_component_allocate, register_simulation_component
  use simcomp_executor, only: neko_simcomps
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    procedure(simulation_component_allocate), pointer :: mma

    ! Assign the pointers
    mma => mma_comp_allocate

    ! Register the simulation components
    call register_simulation_component('mma', mma)

  end subroutine user_setup

  ! ========================================================================== !
  ! Definitions of the simulation component allocators

  !> Allocator for the steady simulation component.
  subroutine mma_comp_allocate(obj)
    class(simulation_component_t), allocatable, intent(inout) :: obj
    allocate(mma_comp_t::obj)
  end subroutine mma_comp_allocate

end module user
