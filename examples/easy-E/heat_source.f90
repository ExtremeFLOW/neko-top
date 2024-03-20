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
  real(kind=rp) :: Pe, target_temperature

  real(kind=rp) :: ramp_time_end = 0.0_rp

contains

  !> Register user defined functions (see nekos user_intf.f90)
  subroutine user_setup(usr)
    type(user_t), intent(inout) :: usr
    usr%scalar_user_f_vector => heat_source
    usr%material_properties => set_material_properties
  end subroutine user_setup

  !> Read the material properties from the JSON file
  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params

    real(kind=rp) :: Re

    call json_get(params, "case.fluid.Re", Re)
    call json_get(params, "case.scalar.Pe", Pe)
    call json_get(params, "case.scalar.target_temperature", target_temperature)
    call json_get(params, "case.end_time", ramp_time_end)
    ramp_time_end = ramp_time_end * 0.01_rp

    mu = 1.0_rp/Re
    lambda = 1.0_rp/Pe
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

  !> Heat source
  subroutine heat_source(rhs, t)
    use math_ext, only: sub3_mask, cfill_mask
    use device_math_ext, only: device_sub3_mask, device_cfill_mask

    class(scalar_user_source_term_t), intent(inout) :: rhs
    real(kind=rp), intent(in) :: t
    type(field_t), pointer :: s

    integer :: i
    class(point_zone_t), pointer :: cyl, ball

    real(kind=rp) :: current_temperature

    cyl => neko_point_zone_registry%get_point_zone("cylinder")
    ball => neko_point_zone_registry%get_point_zone("ball")

    if (.not. allocated(resistance)) then

       allocate (resistance(rhs%dm%size()))
       call rzero(resistance, rhs%dm%size())

       ! Set the resistance to 1.0 in the cylinder
       do i = 1, cyl%size
          resistance(cyl%mask(i)) = target_temperature
       end do

       ! Set the resistance to 1.0 in the ball
       do i = 1, ball%size
          resistance(ball%mask(i)) = target_temperature
       end do

       ! Copy to device
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(resistance, resistance_d, rhs%dm%size())
          call device_memcpy(resistance, resistance_d, rhs%dm%size(), &
                                                                    host_to_device, .true.)
       end if
    end if

    s => neko_field_registry%get_field('s')

    current_temperature = min(1.0_rp, t / ramp_time_end) * target_temperature

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill_mask(s%x_d, current_temperature, rhs%dm%size(), cyl%mask_d, cyl%size)
       call device_cfill_mask(s%x_d, current_temperature, rhs%dm%size(), ball%mask_d, ball%size)


       !  call device_cfill(rhs%s_d, 0.0_rp, rhs%dm%size())
       !  call device_sub3_mask(rhs%s_d, resistance_d, s%x_d, rhs%dm%size(), cyl%mask_d, cyl%size)
       !  call device_sub3_mask(rhs%s_d, resistance_d, s%x_d, rhs%dm%size(), ball%mask_d, ball%size)
    else

       call cfill_mask(s%x, current_temperature, rhs%dm%size(), cyl%mask, cyl%size)
       call cfill_mask(s%x, current_temperature, rhs%dm%size(), ball%mask, ball%size)

       !  call cfill(rhs%s, 0.0_rp, rhs%dm%size())
       !  call sub3_mask(rhs%s, resistance, s%x, rhs%dm%size(), cyl%mask, cyl%size)
       !  call sub3_mask(rhs%s, resistance, s%x, rhs%dm%size(), ball%mask, ball%size)
    end if

  end subroutine heat_source
end module user
