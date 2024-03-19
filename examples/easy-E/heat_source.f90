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
  real(kind=rp) :: Pe

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

    mu = 1.0_rp/Re
    lambda = 1.0_rp/Pe
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

  !> Set the initial condition for the scalar field
  !! @details This function will initialize the scalar field with a two part
  !! uniform value. Above z=0 the scalar field will be 0.0 and below z=0 the
  !! scalar field will be 1.0.
  !!
  !! @param[inout] s Scalar field
  !! @param[inout] params JSON file
  subroutine initial_condition(s, params)
    type(field_t), intent(inout) :: s
    type(json_file), intent(inout) :: params

    real(kind=rp) :: z_value
    integer :: i

    do i = 1, s%dof%size()
       z_value = s%dof%z(i, 1, 1, 1)

       if (z_value > 0.0_rp) then
          s%x(i, 1, 1, 1) = 0.0_rp
       else
          s%x(i, 1, 1, 1) = 1.0_rp
       end if

    end do

  end subroutine initial_condition

  !> Heat source
  subroutine heat_source(rhs, t)
    use math_ext, only: sub3_mask
    use device_math_ext, only: device_sub3_mask

    class(scalar_user_source_term_t), intent(inout) :: rhs
    real(kind=rp), intent(in) :: t
    type(field_t), pointer :: u, v, w, s

    integer :: i
    real(kind=rp) :: target_temp
    class(point_zone_t), pointer :: cyl, ball

    target_temp = 100.0_rp
    cyl => neko_point_zone_registry%get_point_zone("cylinder")
    ball => neko_point_zone_registry%get_point_zone("ball")

    if (.not. allocated(resistance)) then

       allocate (resistance(rhs%dm%size()))
       call rzero(resistance, rhs%dm%size())

       ! Set the resistance to 1.0 in the cylinder
       do i = 1, cyl%size
          resistance(cyl%mask(i)) = target_temp
       end do

       ! Set the resistance to 1.0 in the ball
       do i = 1, ball%size
          resistance(ball%mask(i)) = target_temp
       end do

       ! Copy to device
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(resistance, resistance_d, rhs%dm%size())
          call device_memcpy(resistance, resistance_d, rhs%dm%size(), &
                                                                    host_to_device, .true.)
       end if
    end if

    s => neko_field_registry%get_field('s')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_cfill(rhs%s_d, 0.0_rp, rhs%dm%size())
       call device_sub3_mask(rhs%s_d, resistance_d, s%x_d, rhs%dm%size(), cyl%mask_d, cyl%size)
       call device_sub3_mask(rhs%s_d, resistance_d, s%x_d, rhs%dm%size(), ball%mask_d, ball%size)
    else
       call cfill(rhs%s, 0.0_rp, rhs%dm%size())
       call sub3_mask(rhs%s, resistance, s%x, rhs%dm%size(), cyl%mask, cyl%size)
       call sub3_mask(rhs%s, resistance, s%x, rhs%dm%size(), ball%mask, ball%size)
    end if

  end subroutine heat_source
end module user
