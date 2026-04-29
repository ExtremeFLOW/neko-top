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
  use registry, only : neko_registry
  use math, only : rzero, copy, chsign
  use device_math, only: device_copy, device_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: curl
  use scratch_registry, only : neko_scratch_registry
  use device, only: HOST_TO_DEVICE, device_memcpy
  implicit none
  !> Case parameters
  ! To define the initial boundary conditions we don't wish to introduce a
  ! discontinuity,
  !
  !    DONT                    DO (but smoother)
  !        _______               ______
  !       |                     /
  !       |                    /
  ! -------             -------
  !
  ! hence,
  !> we use a logistic function defined by
  !! $ \frac{L}{1 + exp{-k(z-z_0)}}$
  real(kind=rp), parameter :: L = 1.0_rp
  real(kind=rp), parameter :: k = 20.0_rp
  real(kind=rp), parameter :: z_0 = 1.0_rp

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%dirichlet_conditions => user_bc
    user%initial_conditions => initial_conditions
  end subroutine user_setup

  elemental function inflow_velocity(y, z) result(u_in)
    real(kind=rp), intent(in) :: y, z
    real(kind=rp) :: u_in

    u_in = -0.5_rp * (y - 1.0_rp)**2 - 0.5_rp * (z - 1.0_rp)**2 + 1.0_rp
  end function inflow_velocity

  elemental function scalar_profile(z) result(s_in)
    real(kind=rp), intent(in) :: z
    real(kind=rp) :: s_in

    s_in = L / (1.0_rp + exp(-k*(z - z_0)))
  end function scalar_profile

  !> user-defined boundary condition
  subroutine user_bc(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: y, z
    integer :: i, idx
    logical :: is_fluid

    is_fluid = (fields%items(1)%ptr%name .eq. 'u')

    if (is_fluid) then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          y = u%dof%y(idx, 1, 1, 1)
          z = u%dof%z(idx, 1, 1, 1)

          ! Inflow velocity profile is a paraboloid
          u%x(idx, 1, 1, 1) = inflow_velocity(y, z)
          v%x(idx, 1, 1, 1) = 0._rp
          w%x(idx, 1, 1, 1) = 0._rp
       end do

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(u%x, u%x_d, u%size(), HOST_TO_DEVICE, sync=.false.)
          call device_memcpy(v%x, v%x_d, v%size(), HOST_TO_DEVICE, sync=.false.)
          call device_memcpy(w%x, w%x_d, w%size(), HOST_TO_DEVICE, sync=.false.)
       end if

    else
       s => fields%get("s")

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          z = s%dof%z(idx, 1, 1, 1)
          ! Inflow scalar profile is a sigmoid separating the two species
          s%x(idx, 1, 1, 1) = scalar_profile(z)
       end do
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(s%x, s%x_d, s%size(), HOST_TO_DEVICE, sync=.false.)
       end if
    end if
  end subroutine user_bc

  !> user-defined initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: u, v, w, s
    integer :: i

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       do i = 1, u%dof%size()
          u%x(i,1,1,1) = inflow_velocity(u%dof%y(i,1,1,1), u%dof%z(i,1,1,1))
          v%x(i,1,1,1) = 0.0_rp
          w%x(i,1,1,1) = 0.0_rp
       end do

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(u%x, u%x_d, u%size(), HOST_TO_DEVICE, sync=.false.)
          call device_memcpy(v%x, v%x_d, v%size(), HOST_TO_DEVICE, sync=.false.)
          call device_memcpy(w%x, w%x_d, w%size(), HOST_TO_DEVICE, sync=.false.)
       end if

       return
    end if

    ! Initial scalar profile is a sigmoid separating the two species
    s => fields%get("s")
    do i = 1, s%dof%size()
       s%x(i,1,1,1) = scalar_profile(s%dof%z(i,1,1,1))
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(s%x, s%x_d, s%size(), HOST_TO_DEVICE, sync=.false.)
    end if

  end subroutine initial_conditions

end module user

