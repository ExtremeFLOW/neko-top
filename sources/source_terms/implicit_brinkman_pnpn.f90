!> @file implicit_brinkman_pnpn.f90
!! Helpers for the experimental implicit Brinkman Pn/Pn projection.
module implicit_brinkman_pnpn
  use num_types, only: rp
  use coefs, only: coef_t
  use field, only: field_t
  use pnpn_residual, only: pnpn_pressure_coef_hook, &
       pnpn_pressure_rhs_hook, pnpn_pressure_bc_rhs_hook, &
       pnpn_velocity_res_hook
  use adjoint_pnpn_residual, only: adjoint_pnpn_pressure_coef_hook, &
       adjoint_pnpn_pressure_rhs_hook, adjoint_pnpn_projection_hook
  use neko_config, only: NEKO_BCKND_DEVICE
  use device, only: device_memcpy, HOST_TO_DEVICE, DEVICE_TO_HOST
  implicit none
  private

  public :: setup_implicit_brinkman_pnpn, clear_implicit_brinkman_pnpn

  type(field_t), pointer :: brinkman_chi => null()
  type(field_t), pointer :: brinkman_u_sens => null()
  type(field_t), pointer :: brinkman_v_sens => null()
  type(field_t), pointer :: brinkman_w_sens => null()

contains

  subroutine setup_implicit_brinkman_pnpn(chi, u_sens, v_sens, w_sens)
    type(field_t), target, intent(in) :: chi
    type(field_t), target, intent(inout), optional :: u_sens
    type(field_t), target, intent(inout), optional :: v_sens
    type(field_t), target, intent(inout), optional :: w_sens

    brinkman_chi => chi
    if (present(u_sens)) brinkman_u_sens => u_sens
    if (present(v_sens)) brinkman_v_sens => v_sens
    if (present(w_sens)) brinkman_w_sens => w_sens
    pnpn_pressure_coef_hook => implicit_brinkman_pressure_coef
    pnpn_pressure_rhs_hook => implicit_brinkman_scale_vector
    pnpn_pressure_bc_rhs_hook => implicit_brinkman_scale_vector
    pnpn_velocity_res_hook => implicit_brinkman_velocity_residual
    adjoint_pnpn_pressure_coef_hook => implicit_brinkman_pressure_coef
    adjoint_pnpn_pressure_rhs_hook => implicit_brinkman_scale_vector
    adjoint_pnpn_projection_hook => implicit_brinkman_scale_vector
  end subroutine setup_implicit_brinkman_pnpn

  subroutine clear_implicit_brinkman_pnpn()
    nullify(brinkman_chi)
    nullify(brinkman_u_sens)
    nullify(brinkman_v_sens)
    nullify(brinkman_w_sens)
    nullify(pnpn_pressure_coef_hook)
    nullify(pnpn_pressure_rhs_hook)
    nullify(pnpn_pressure_bc_rhs_hook)
    nullify(pnpn_velocity_res_hook)
    nullify(adjoint_pnpn_pressure_coef_hook)
    nullify(adjoint_pnpn_pressure_rhs_hook)
    nullify(adjoint_pnpn_projection_hook)
  end subroutine clear_implicit_brinkman_pnpn

  subroutine implicit_brinkman_pressure_coef(c_Xh, bd, dt, mu, rho, n)
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    real(kind=rp) :: rho_val
    integer :: i

    if (.not. associated(brinkman_chi)) return

    rho_val = rho%x(1,1,1,1)

    do i = 1, n
       c_Xh%h1(i,1,1,1) = (1.0_rp / rho_val) * &
            bd / (bd + dt * brinkman_chi%x(i,1,1,1))
       c_Xh%h2(i,1,1,1) = 0.0_rp
    end do
    c_Xh%ifh2 = .false.

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(c_Xh%h1, c_Xh%h1_d, n, HOST_TO_DEVICE, &
            sync = .false.)
       call device_memcpy(c_Xh%h2, c_Xh%h2_d, n, HOST_TO_DEVICE, &
            sync = .true.)
    end if
  end subroutine implicit_brinkman_pressure_coef

  subroutine implicit_brinkman_scale_vector(x, y, z, c_Xh, bd, dt, mu, rho, n)
    type(field_t), intent(inout) :: x
    type(field_t), intent(inout) :: y
    type(field_t), intent(inout) :: z
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    real(kind=rp) :: b
    integer :: i

    if (.not. associated(brinkman_chi)) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(x%x, x%x_d, n, DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(y%x, y%x_d, n, DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(z%x, z%x_d, n, DEVICE_TO_HOST, sync = .true.)
    end if

    do i = 1, n
       b = bd / (bd + dt * brinkman_chi%x(i,1,1,1))
       x%x(i,1,1,1) = b * x%x(i,1,1,1)
       y%x(i,1,1,1) = b * y%x(i,1,1,1)
       z%x(i,1,1,1) = b * z%x(i,1,1,1)
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(x%x, x%x_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(y%x, y%x_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(z%x, z%x_d, n, HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine implicit_brinkman_scale_vector

  subroutine implicit_brinkman_velocity_residual(u_res, v_res, w_res, gradp_x, &
       gradp_y, gradp_z, f_x, f_y, f_z, c_Xh, bd, dt, mu, rho, n)
    type(field_t), intent(inout) :: u_res
    type(field_t), intent(inout) :: v_res
    type(field_t), intent(inout) :: w_res
    type(field_t), intent(inout) :: gradp_x
    type(field_t), intent(inout) :: gradp_y
    type(field_t), intent(inout) :: gradp_z
    type(field_t), intent(in) :: f_x
    type(field_t), intent(in) :: f_y
    type(field_t), intent(in) :: f_z
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    real(kind=rp) :: b, u_projected, v_projected, w_projected
    integer :: i

    if (.not. associated(brinkman_chi)) return

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_memcpy(u_res%x, u_res%x_d, n, DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(v_res%x, v_res%x_d, n, DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(w_res%x, w_res%x_d, n, DEVICE_TO_HOST, sync = .false.)
       call device_memcpy(gradp_x%x, gradp_x%x_d, n, DEVICE_TO_HOST, &
            sync = .false.)
       call device_memcpy(gradp_y%x, gradp_y%x_d, n, DEVICE_TO_HOST, &
            sync = .false.)
       call device_memcpy(gradp_z%x, gradp_z%x_d, n, DEVICE_TO_HOST, &
            sync = .true.)
    end if

    do i = 1, n
       b = bd / (bd + dt * brinkman_chi%x(i,1,1,1))
       u_projected = b * (f_x%x(i,1,1,1) - gradp_x%x(i,1,1,1))
       v_projected = b * (f_y%x(i,1,1,1) - gradp_y%x(i,1,1,1))
       w_projected = b * (f_z%x(i,1,1,1) - gradp_z%x(i,1,1,1))
       if (associated(brinkman_u_sens)) then
          brinkman_u_sens%x(i,1,1,1) = u_projected / bd
          brinkman_v_sens%x(i,1,1,1) = v_projected / bd
          brinkman_w_sens%x(i,1,1,1) = w_projected / bd
       end if
       u_res%x(i,1,1,1) = -u_res%x(i,1,1,1) + u_projected
       v_res%x(i,1,1,1) = -v_res%x(i,1,1,1) + v_projected
       w_res%x(i,1,1,1) = -w_res%x(i,1,1,1) + w_projected
    end do

    if (NEKO_BCKND_DEVICE .eq. 1) then
       if (associated(brinkman_u_sens)) then
          call device_memcpy(brinkman_u_sens%x, brinkman_u_sens%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(brinkman_v_sens%x, brinkman_v_sens%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
          call device_memcpy(brinkman_w_sens%x, brinkman_w_sens%x_d, n, &
               HOST_TO_DEVICE, sync = .false.)
       end if
       call device_memcpy(u_res%x, u_res%x_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(v_res%x, v_res%x_d, n, HOST_TO_DEVICE, sync = .false.)
       call device_memcpy(w_res%x, w_res%x_d, n, HOST_TO_DEVICE, sync = .true.)
    end if
  end subroutine implicit_brinkman_velocity_residual

end module implicit_brinkman_pnpn
