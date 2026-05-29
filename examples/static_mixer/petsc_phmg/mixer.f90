module user
  use neko
  use user_initial_conditions, only: scalar_z_split_ic
  use field_math, only: field_rzero
  implicit none

  real(kind=rp) :: Re, element_size

contains

  ! Register user defined functions (see user_intf.f90)
  subroutine user_setup(user)
    type(user_t), intent(inout) :: user
    user%initial_conditions => initial_conditions
    user%dirichlet_conditions => dirichlet_update
    user%startup => startup
  end subroutine user_setup

  subroutine startup(params)
    type(json_file), intent(inout) :: params
    call json_get(params, "case.fluid.Re", Re)
    call json_get(params, "case.element_size", element_size)
  end subroutine startup

  !> User initial condition
  subroutine initial_conditions(scheme_name, fields)
    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: fields

    type(field_t), pointer :: u, v, w, s, brinkman_indicator

    brinkman_indicator => neko_registry%get_field("brinkman_indicator")

    ! See scalar.name in the case file, makes sure that we only
    ! run this for the scalar field.
    if (scheme_name .eq. 'temperature') then
       s => fields%items(1)%ptr
       call scalar_z_split_ic(s, 0.5_rp, 0.0_rp, 1.0_rp)
    end if

    if (scheme_name .eq. 'fluid') then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       call field_cfill(u, 1.0_rp)
       call field_subcol3(u, u, brinkman_indicator)
       call field_rzero(v)
       call field_rzero(w)
    end if
  end subroutine initial_conditions

  subroutine dirichlet_update(fields, bc, time)
    type(field_list_t), intent(inout) :: fields
    type(field_dirichlet_t), intent(in) :: bc
    type(time_state_t), intent(in) :: time

    type(field_t), pointer :: u, v, w, s
    integer :: i
    real(kind=rp) :: band_size
    real(kind=rp) :: y, z, val_y0, val_z0, val_y1, val_z1
    logical :: is_velocity_bc, is_scalar_bc

    if (time%tstep .ne. 1) return

    is_velocity_bc = fields%items(1)%ptr%name == "u"
    is_scalar_bc = fields%items(1)%ptr%name == "temperature"

    ! ------------------------------------------------------------------------ !
    ! Assign a boundary condition for the velocity field here if needed.

    if (is_velocity_bc) then
       u => fields%get("u")
       v => fields%get("v")
       w => fields%get("w")

       band_size = 1.0_rp / sqrt(Re)
       band_size = max(band_size, element_size)

       do i = 1, bc%msk(0)
          y = bc%dof%y(bc%msk(i), 1, 1, 1)
          z = bc%dof%z(bc%msk(i), 1, 1, 1)

          val_y0 = smooth_step(y, 0.0_rp, band_size)
          val_z0 = smooth_step(z, 0.0_rp, band_size)
          val_y1 = smooth_step(y, 1.0_rp, 1.0_rp - band_size)
          val_z1 = smooth_step(z, 1.0_rp, 1.0_rp - band_size)

          u%x(bc%msk(i), 1, 1, 1) = val_y0 * val_z0 * val_y1 * val_z1
       end do

       call u%copy_from(HOST_TO_DEVICE, .true.)
       call field_rzero(v)
       call field_rzero(w)
    end if

    ! ------------------------------------------------------------------------ !
    ! Assign the boundary condition for the scalar field here.

    if (is_scalar_bc) then
       s => fields%get("temperature")

       do i = 1, bc%msk(0)
          z = bc%dof%z(bc%msk(i), 1, 1, 1)
          s%x(bc%msk(i), 1, 1, 1) = smooth_step(z, 0.45_rp, 0.55_rp)
       end do

       call s%copy_from(HOST_TO_DEVICE, .true.)
    end if

  end subroutine dirichlet_update


  !> @brief Apply a smooth step function to a scalar.
  elemental function smooth_step(x, edge0, edge1) result(res)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: edge0, edge1
    real(kind=rp) :: res, t

    t = clamp((x - edge0) / (edge1 - edge0), 0.0_rp, 1.0_rp)

    res = t**3 * (t * (6.0_rp * t - 15.0_rp) + 10.0_rp)

  end function smooth_step

  !> @brief Clamp a value between two limits.
  elemental function clamp(x, lowerlimit, upperlimit) result(res)
    real(kind=rp), intent(in) :: x
    real(kind=rp), intent(in) :: lowerlimit, upperlimit
    real(kind=rp) :: res

    res = max(lowerlimit, min(upperlimit, x))
  end function clamp

end module user
