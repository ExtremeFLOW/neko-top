! User module for the user defined simulation component
module user
  use user_intf, only: user_t
  use json_module, only: json_file
  use json_utils, only: json_get_or_default
  use steady_simcomp, only: steady_simcomp_t
  use simcomp_executor, only: neko_simcomps
  use num_types, only : rp
  use field, only : field_t
  use field_list, only : field_list_t
  use field_math, only : field_cfill, field_sub2, field_col2, field_glmax, &
       field_copy, field_cpwmax2
  use field_dirichlet, only : field_dirichlet_t
  use time_state, only : time_state_t
  use registry, only : neko_registry
  use math, only : rzero, copy, chsign, NEKO_EPS
  use device_math, only: device_copy, device_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use operators, only: curl
  use scratch_registry, only : neko_scratch_registry
  use device, only: HOST_TO_DEVICE, device_memcpy
  use utils, only: neko_error
  use logger, only: neko_log
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
  real(kind=rp), parameter :: k = 40.0_rp
  real(kind=rp), parameter :: z_0 = 1.0_rp
  logical :: zero_ramp = .false.

contains

  ! Register user-defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%startup => user_read_params
    user%initialize => user_intialize
    user%dirichlet_conditions => user_bc
    user%initial_conditions => scalar_ic
  end subroutine user_setup

  !> Read parameters from the json file
  subroutine user_read_params(json)
    type(json_file), intent(inout) :: json

    call json_get_or_default(json, "optimization.zero_fields", zero_ramp, .false.)

  end subroutine user_read_params

  !> user-defined initialization function
  subroutine user_intialize(time)
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w
    type(field_t), pointer :: design, tmp
    integer :: tmp_idx

    if (.not. zero_ramp) return

    ! Retrieve pointers to the fields.
    u => neko_registry%get_field("u")
    v => neko_registry%get_field("v")
    w => neko_registry%get_field("w")

    design => neko_registry%get_field("RAMP_mapping")

    call neko_scratch_registry%request(tmp, tmp_idx, clear=.true.)

    ! We need to zero out all velocities when design is 1 (solid)
    call field_cfill(tmp, 1.0_rp)
    call field_sub2(tmp, design)
    call field_cpwmax2(tmp, 10.0_rp * NEKO_EPS)

    call field_col2(u, tmp)

    call neko_scratch_registry%relinquish(tmp_idx)

  end subroutine user_intialize

  !> user-defined boundary condition
  subroutine user_bc(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: u, v, w, s
    real(kind=rp) :: x, y, z
    integer :: i, idx
    logical :: is_fluid

    is_fluid = (fields%items(1)%ptr%name .eq. 'u')

    if (is_fluid) then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          x = u%dof%x(idx, 1, 1, 1)
          y = u%dof%y(idx, 1, 1, 1)
          z = u%dof%z(idx, 1, 1, 1)

          ! Inflow velocity profile is a paraboloid
          u%x(idx, 1, 1, 1) = -0.5_rp * (y - 1.0_rp)**2 - &
               0.5_rp * (z - 1.0_rp)**2 + 1.0_rp
          v%x(idx, 1, 1, 1) = 0._rp
          w%x(idx, 1, 1, 1) = 0._rp
       end do

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(u%x, u%x_d, u%size(), HOST_TO_DEVICE, sync=.false.)
          call device_memcpy(v%x, v%x_d, v%size(), HOST_TO_DEVICE, sync=.false.)
          call device_memcpy(w%x, w%x_d, w%size(), HOST_TO_DEVICE, sync=.false.)
       end if

    else
       s => fields%get("Scalar")

       do i = 1, bc%msk(0)
          idx = bc%msk(i)
          z = s%dof%z(idx, 1, 1, 1)
          ! Inflow scalar profile is a sigmoid separating the two species
          s%x(idx, 1, 1, 1) = logistic(z)
       end do
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(s%x, s%x_d, s%size(), HOST_TO_DEVICE, sync=.false.)
       end if
    end if
  end subroutine user_bc

  !> user-defined initial condition
  subroutine scalar_ic(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields
    type(field_t), pointer :: s, u, v, w
    real(kind=rp) :: x, y, z
    integer :: i

    if(trim(scheme_name) .eq. 'Scalar') then

       ! Initial scalar profile is a sigmoid separating the two species
       s => fields%get("Scalar")
       s%x = logistic(s%dof%z)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(s%x, s%x_d, s%size(), HOST_TO_DEVICE, sync=.true.)
       end if

    else if (trim(scheme_name) .eq. 'Fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       do i = 1, u%size()
          x = u%dof%x(i, 1, 1, 1)
          y = u%dof%y(i, 1, 1, 1)
          z = u%dof%z(i, 1, 1, 1)

          ! Inflow velocity profile is a paraboloid
          u%x(i, 1, 1, 1) = -0.5_rp * (y - 1.0_rp)**2 - &
               0.5_rp * (z - 1.0_rp)**2 + 1.0_rp
       end do

       call field_cfill(v, 0.0_rp)
       call field_cfill(w, 0.0_rp)

       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_memcpy(w%x, w%x_d, w%size(), HOST_TO_DEVICE, sync=.true.)
       end if
    end if

  end subroutine scalar_ic

  elemental function logistic(x) result(res)
    real(kind=rp), intent(in) :: x
    real(kind=rp) :: res
    res = L / (1.0_rp + exp(-k*(x - z_0)))
  end function logistic

end module user
