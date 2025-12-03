!> @file user.f90
!> @brief User defined user region
!>
!> This file is part of Neko-TOP.

!> @brief User defined user region
module user
  use neko

  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr

  implicit none

  real(kind=rp), dimension(:), allocatable :: resistance
  type(c_ptr) :: resistance_d = c_null_ptr
  logical :: is_initialized = .false.

  real(kind=rp) :: target_temperature = 0.0_rp
  real(kind=rp) :: ramp_time_end = 0.0_rp

  !> Variables to store the Rayleigh and Peclet numbers
  real(kind=rp) :: Re = 0.0_rp
  real(kind=rp) :: Pe = 0.0_rp

contains

  !> Register user defined functions (see nekos user_intf.f90)
  subroutine user_setup(usr)
    type(user_t), intent(inout) :: usr
    usr%source_term => heat_source
    usr%startup => user_startup
  end subroutine user_setup

  !> Initialize the user module
  subroutine user_startup(params)
    type(json_file), intent(inout) :: params

    call json_get(params, "case.fluid.Re", Re)
    call json_get(params, "case.scalar.Pe", Pe)
    call json_get(params, "case.scalar.target_temperature", target_temperature)
    call json_get(params, "case.time.end_time", ramp_time_end)
    ramp_time_end = ramp_time_end * 0.01_rp

  end subroutine user_startup

  !> Read the material properties from the JSON file
  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params

    mu = 1.0_rp/Re
    lambda = 1.0_rp/Pe
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

  !> Heat source
  subroutine heat_source(scheme_name, rhs, time)
    use math_ext, only: sub3_mask
    use device_math_ext, only: device_sub3_mask

    character(len=*), intent(in) :: scheme_name
    type(field_list_t), intent(inout) :: rhs
    type(time_state_t), intent(in) :: time
    type(field_t), pointer :: s, rhs_s

    integer :: i
    class(point_zone_t), pointer :: cyl, ball

    real(kind=rp) :: current_temperature

    if (scheme_name .ne. 'scalar') return

    cyl => neko_point_zone_registry%get_point_zone("cylinder")
    ball => neko_point_zone_registry%get_point_zone("ball")
    rhs_s => rhs%get('s')

    if (.not. allocated(resistance)) then

       allocate (resistance(rhs_s%dof%size()))
       call rzero(resistance, rhs_s%dof%size())

       ! Set the resistance to 1.0 in the cylinder
       do i = 1, cyl%size
          resistance(cyl%mask%get(i)) = target_temperature
       end do

       ! Set the resistance to 1.0 in the ball
       do i = 1, ball%size
          resistance(ball%mask%get(i)) = target_temperature
       end do

       ! Copy to device
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(resistance, resistance_d, rhs_s%dof%size())
          call device_memcpy(resistance, resistance_d, rhs_s%dof%size(), &
               host_to_device, .true.)
       end if
    end if

    s => neko_registry%get_field('s')

    current_temperature = min(1.0_rp, time%t / ramp_time_end) * target_temperature

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill_mask(s%x_d, current_temperature, rhs_s%dof%size(), &
            cyl%mask%get_d(), cyl%size)
       call device_cfill_mask(s%x_d, current_temperature, rhs_s%dof%size(), &
            ball%mask%get_d(), ball%size)


       !  call device_cfill(rhs_s%s_d, 0.0_rp, rhs_s%dof%size())
       !  call device_sub3_mask(rhs_s%s_d, resistance_d, s%x_d, rhs_s%dof%size(), &
       !       cyl%mask_d, cyl%size)
       !  call device_sub3_mask(rhs_s%s_d, resistance_d, s%x_d, rhs_s%dof%size(), &
       !       ball%mask_d, ball%size)
    else

       call cfill_mask(s%x, current_temperature, rhs_s%dof%size(), cyl%mask%get(), &
            cyl%size)
       call cfill_mask(s%x, current_temperature, rhs_s%dof%size(), ball%mask%get(), &
            ball%size)

       !  call cfill(rhs_s%s, 0.0_rp, rhs_s%dof%size())
       !  call sub3_mask(rhs_s%s, resistance, s%x, rhs_s%dof%size(), cyl%mask, &
       !       cyl%size)
       !  call sub3_mask(rhs_s%s, resistance, s%x, rhs_s%dof%size(), ball%mask, &
       !       ball%size)
    end if

  end subroutine heat_source
end module user
