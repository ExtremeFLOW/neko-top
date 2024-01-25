!> @file user.f90
!> @brief User defined user region
!>
!> This file is part of Neko-TOP.

!> @brief User defined user region
module user
  use neko

  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr

  implicit none

  real(kind=rp) :: perm = 1.0_rp
  real(kind=rp), allocatable :: resistance(:)
  type(c_ptr) :: resistance_d = c_null_ptr
  logical :: is_initialized = .false.

contains

  !> Register user defined functions (see nekos user_intf.f90)
  subroutine user_setup(u)
    type(user_t), intent(inout) :: u
    u%fluid_user_f_vector => fluid_permeability
    u%material_properties => set_material_properties
    u%scalar_user_ic => initial_condition
  end subroutine user_setup

  subroutine set_material_properties(t, tstep, rho, mu, cp, lambda, params)
    real(kind=rp), intent(in) :: t
    integer, intent(in) :: tstep
    real(kind=rp), intent(inout) :: rho, mu, cp, lambda
    type(json_file), intent(inout) :: params

    real(kind=rp) :: Re, Pe

    call json_get(params, "case.fluid.Re", Re)
    call json_get(params, "case.fluid.perm", perm)
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

  !> Forcing
  subroutine fluid_permeability(f, t)
    class(fluid_user_source_term_t), intent(inout) :: f
    real(kind=rp), intent(in) :: t
    type(field_t), pointer :: u, v, w, s

    integer :: i
    class(point_zone_t), pointer :: my_point_zone

    if (.not. is_initialized) then
       is_initialized = .true.

       my_point_zone => neko_point_zone_registry%get_point_zone("lowperm")

       allocate (resistance(f%dm%size()))
       call rzero(resistance, f%dm%size())

       do i = 1, my_point_zone%size
          resistance(my_point_zone%mask(i)) = 1.0_rp - 1.0_rp/perm
       end do

       ! Copy to device
       if (NEKO_BCKND_DEVICE .eq. 1) then
          call device_map(resistance, resistance_d, f%dm%size())
          call device_memcpy(resistance, resistance_d, f%dm%size(),&
                                                                  host_to_device, .true.)
       end if
    end if

    u => neko_field_registry%get_field('u')
    v => neko_field_registry%get_field('v')
    w => neko_field_registry%get_field('w')

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col3(f%u_d, u%x_d, resistance_d, f%dm%size())
       call device_col3(f%v_d, v%x_d, resistance_d, f%dm%size())
       call device_col3(f%w_d, w%x_d, resistance_d, f%dm%size())
    else
       call col3(f%u, u%x, resistance, f%dm%size())
       call col3(f%v, v%x, resistance, f%dm%size())
       call col3(f%w, w%x, resistance, f%dm%size())
    end if
  end subroutine fluid_permeability
end module user
