! User module for the user defined simulation component
module user
  use user_intf, only: user_t, simulation_component_user_settings
  use json_module, only: json_file
  use simcomp_example, only: adjoint_t
  use simcomp_executor, only: neko_simcomps
  implicit none

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%init_user_simcomp => user_simcomp
  end subroutine user_setup

  subroutine user_simcomp(params)
    type(json_file), intent(inout) :: params
    type(adjoint_t), allocatable :: my_simcomp
    type(json_file) :: simcomp_settings

    ! Allocate a simulation component
    allocate(my_simcomp)
    simcomp_settings = simulation_component_user_settings("adjoint", params)

    call neko_simcomps%add_user_simcomp(my_simcomp, simcomp_settings)

  end subroutine user_simcomp

end module user
