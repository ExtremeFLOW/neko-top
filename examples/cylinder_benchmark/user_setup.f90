!> @file user.f90
!> @brief User defined user region
!>
!> This file is part of Neko-TOP.

!> @brief User defined user region
module user
  use neko
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr

  implicit none

contains

  !> Register user defined functions (see nekos user_intf.f90)
  subroutine user_setup(usr)
    type(user_t), intent(inout) :: usr
    usr%fluid_user_ic => cylinder_ic
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
    mu = 1.0_rp/Re
    lambda = 1.0_rp
    rho = 1.0_rp
    cp = 1.0_rp
  end subroutine set_material_properties

  !> Set the initial condition
  !! @details Set the initial condition for the velocity and pressure fields.
  !! The conditions is set to 1 everywhere and 0 in the cylinder.
  !! Additionally we wish to add some noise to the initial condition.
  !! @param u velocity field
  !! @param v velocity field
  !! @param w velocity field
  !! @param p pressure field
  !! @param params JSON file
  subroutine cylinder_ic(u, v, w, p, params)
    type(field_t), intent(inout) :: u, v, w, p
    type(json_file), intent(inout) :: params

    class(point_zone_t), pointer :: cylinder
    real(kind=rp) :: noise, noise_scale
    integer :: i

    ! Set a uniform flow field.
    call cfill(u%x, 1.0_rp, u%dof%size())
    call cfill(v%x, 0.0_rp, v%dof%size())
    call cfill(w%x, 0.0_rp, w%dof%size())

    ! Apply a random perturbation to the initial condition.
    noise = 0.0_rp
    noise_scale = 1e-2_rp
    do i = 1, u%dof%size()
       call random_number(noise)
       u%x(i, 1, 1, 1) = u%x(i, 1, 1, 1) + noise_scale * noise
       call random_number(noise)
       v%x(i, 1, 1, 1) = v%x(i, 1, 1, 1) + noise_scale * noise
       call random_number(noise)
       w%x(i, 1, 1, 1) = w%x(i, 1, 1, 1) + noise_scale * noise
    end do

    ! Set flow to zero in the cylinder.
    if (neko_point_zone_registry%point_zone_exists("cylinder")) then
       cylinder => neko_point_zone_registry%get_point_zone("cylinder")

       call cfill_mask(u%x, 0.0_rp, u%dof%size(), cylinder%mask, cylinder%size)
    end if

  end subroutine cylinder_ic
end module user
