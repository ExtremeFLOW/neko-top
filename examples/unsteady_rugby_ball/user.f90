! User module for the user defined simulation component
module user
  use user_intf, only: user_t
  use json_module, only: json_file
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use field_dirichlet, only : field_dirichlet_t
  use time_state, only : time_state_t
  use field_registry, only : neko_field_registry
  use math, only : rzero, copy, chsign
  use device_math, only: device_copy, device_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use scratch_registry, only : neko_scratch_registry
  implicit none
  !> Case parameters
  real(kind=rp), parameter :: T_fin = 1.0_rp
  real(kind=rp), parameter :: U_max = 1.0_rp

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%dirichlet_conditions => user_bc
  end subroutine user_setup

  !> user-defined boundary condition
  subroutine user_bc(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    integer :: i, idx
    logical :: is_fluid

    is_fluid = (fields%items(1)%ptr%name .eq. 'u')
    print *, "poo"

    if (is_fluid) then
     associate(u => field_bc_list%items(1)%ptr, &
               v => field_bc_list%items(2)%ptr, &
               w => field_bc_list%items(3)%ptr)

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          ! Inflow velocity profile is a ramp up
          print *, "yof"
          if (time%t .lt. T_fin) then
          print *, "1"
               u%x(idx, 1, 1, 1) = U_max * time%t / T_fin
          else
          print *, "2"
               u%x(idx, 1, 1, 1) = U_max
          end if
          u%x(idx, 1, 1, 1) = 1.0_rp
          v%x(idx, 1, 1, 1) = 0.0_rp
          w%x(idx, 1, 1, 1) = 0.0_rp
       end do
      end associate
    end if
  end subroutine user_bc

end module user
