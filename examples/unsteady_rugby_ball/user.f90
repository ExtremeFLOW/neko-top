! User module for the user defined simulation component
module user
  use user_intf, only: user_t
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use field_dirichlet, only : field_dirichlet_t
  use time_state, only : time_state_t
  use field_registry, only : neko_field_registry
  use neko_config, only: NEKO_BCKND_DEVICE
  implicit none
  !> Case parameters
  real(kind=rp), parameter :: T_fin = 0.05_rp
  real(kind=rp), parameter :: U_max = 1.0_rp

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%dirichlet_conditions => dirichlet_conditions
  end subroutine user_setup

  !> user-defined boundary condition
  subroutine dirichlet_conditions(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    integer :: i, idx
    logical :: is_fluid
    type(field_t), pointer :: u, v, w

    is_fluid = (fields%items(1)%ptr%name .eq. 'u')

    if (is_fluid) then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          ! Inflow velocity profile is a ramp up
          if (time%t .lt. T_fin) then
               u%x(idx, 1, 1, 1) = U_max * time%t / T_fin
          else
               u%x(idx, 1, 1, 1) = U_max
          end if
          v%x(idx, 1, 1, 1) = 0.0_rp
          w%x(idx, 1, 1, 1) = 0.0_rp
       end do
    end if
  end subroutine dirichlet_conditions

end module user
