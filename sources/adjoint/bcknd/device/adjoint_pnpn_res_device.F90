!> Residuals in the Pn-Pn formulation (device backend)
module adjoint_pnpn_res_device
  use gather_scatter, only : gs_t, GS_OP_ADD
  use operators, only : cdtp
  use field, only : field_t
  use ax_product, only : ax_t
  use coefs, only : coef_t
  use facet_normal, only : facet_normal_t
  use adjoint_pnpn_residual, only : adjoint_pnpn_prs_res_t, adjoint_pnpn_vel_res_t
  use scratch_registry, only: neko_scratch_registry
  use mesh, only : mesh_t
  use num_types, only : rp, c_rp
  use space, only : space_t
  use device_math, only : device_cfill, device_col2
  use utils, only : neko_error
  use device, only : device_event_sync
  use, intrinsic :: iso_c_binding, only : c_ptr, c_int
  implicit none
  private

  type, public, extends(adjoint_pnpn_prs_res_t) :: adjoint_pnpn_prs_res_device_t
   contains
     procedure, nopass :: compute => adjoint_pnpn_prs_res_device_compute
  end type adjoint_pnpn_prs_res_device_t

  type, public, extends(adjoint_pnpn_vel_res_t) :: adjoint_pnpn_vel_res_device_t
   contains
     procedure, nopass :: compute => adjoint_pnpn_vel_res_device_compute
  end type adjoint_pnpn_vel_res_device_t

contains

#ifdef HAVE_HIP
  interface
     subroutine adjoint_prs_res_part1_hip(ta1_d, ta2_d, ta3_d, &
          f_u_d, f_v_d, f_w_d, B_d, inv_rho, n) &
          bind(c, name = 'adjoint_prs_res_part1_hip')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d
       real(c_rp) :: inv_rho
       integer(c_int) :: n
     end subroutine adjoint_prs_res_part1_hip
  end interface

  interface
     subroutine adjoint_prs_res_part2_hip(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name = 'adjoint_prs_res_part2_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine adjoint_prs_res_part2_hip
  end interface

  interface
     subroutine adjoint_vel_res_update_hip(u_res_d, v_res_d, w_res_d, &
          f_u_d, f_v_d, f_w_d, n) &
          bind(c, name = 'adjoint_vel_res_update_hip')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine adjoint_vel_res_update_hip
  end interface
#elif HAVE_CUDA
  interface
     subroutine adjoint_prs_res_part1_cuda(ta1_d, ta2_d, ta3_d, &
          f_u_d, f_v_d, f_w_d, B_d, inv_rho, n) &
          bind(c, name = 'adjoint_prs_res_part1_cuda')
       use, intrinsic :: iso_c_binding
       import c_rp
       implicit none
       type(c_ptr), value :: ta1_d, ta2_d, ta3_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       type(c_ptr), value :: B_d
       real(c_rp) :: inv_rho
       integer(c_int) :: n
     end subroutine adjoint_prs_res_part1_cuda
  end interface

  interface
     subroutine adjoint_prs_res_part2_cuda(p_res_d, wa1_d, wa2_d, wa3_d, n) &
          bind(c, name = 'adjoint_prs_res_part2_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: p_res_d, wa1_d, wa2_d, wa3_d
       integer(c_int) :: n
     end subroutine adjoint_prs_res_part2_cuda
  end interface

  interface
     subroutine adjoint_vel_res_update_cuda(u_res_d, v_res_d, w_res_d, &
          f_u_d, f_v_d, f_w_d, n) &
          bind(c, name = 'adjoint_vel_res_update_cuda')
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: u_res_d, v_res_d, w_res_d
       type(c_ptr), value :: f_u_d, f_v_d, f_w_d
       integer(c_int) :: n
     end subroutine adjoint_vel_res_update_cuda
  end interface
#endif

  subroutine adjoint_pnpn_prs_res_device_compute(p, p_res, u, v, w, u_e, v_e, w_e, f_x, &
       f_y, f_z, c_Xh, gs_Xh, bc_prs_surface, bc_sym_surface, Ax, bd, dt, mu, &
       rho, event)
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(in) :: u_e, v_e, w_e
    type(field_t), intent(inout) :: p_res
    type(field_t), intent(in) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(gs_t), intent(inout) :: gs_Xh
    type(facet_normal_t), intent(in) :: bc_prs_surface
    type(facet_normal_t), intent(in) :: bc_sym_surface
    class(ax_t), intent(inout) :: Ax
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    type(c_ptr), intent(inout) :: event
    real(kind=rp) :: rho_val, inv_rho
    integer :: n
    type(field_t), pointer :: ta1, ta2, ta3, wa1, wa2, wa3
    integer :: temp_indices(6)

    call neko_scratch_registry%request_field(ta1, temp_indices(1), .false.)
    call neko_scratch_registry%request_field(ta2, temp_indices(2), .false.)
    call neko_scratch_registry%request_field(ta3, temp_indices(3), .false.)
    call neko_scratch_registry%request_field(wa1, temp_indices(4), .false.)
    call neko_scratch_registry%request_field(wa2, temp_indices(5), .false.)
    call neko_scratch_registry%request_field(wa3, temp_indices(6), .false.)

    n = u%dof%size()

    ! We assume the material properties are constant
    rho_val = rho%x(1,1,1,1)
    inv_rho = 1.0_rp / rho_val

    c_Xh%ifh2 = .false.
    call device_cfill(c_Xh%h1_d, inv_rho, n)
    call device_cfill(c_Xh%h2_d, 0.0_rp, n)

#ifdef HAVE_HIP
    call adjoint_prs_res_part1_hip(ta1%x_d, ta2%x_d, ta3%x_d, &
         f_x%x_d, f_y%x_d, f_z%x_d, c_Xh%B_d, inv_rho, n)
#elif HAVE_CUDA
    call adjoint_prs_res_part1_cuda(ta1%x_d, ta2%x_d, ta3%x_d, &
         f_x%x_d, f_y%x_d, f_z%x_d, c_Xh%B_d, inv_rho, n)
#elif HAVE_OPENCL
    call neko_error('adjoint_prs_res_part1_opencl not implemented')
#else
    call neko_error('No device backend configured')
#endif

    call gs_Xh%op(ta1, GS_OP_ADD, event)
    call device_event_sync(event)
    call gs_Xh%op(ta2, GS_OP_ADD, event)
    call device_event_sync(event)
    call gs_Xh%op(ta3, GS_OP_ADD, event)
    call device_event_sync(event)

    call device_col2(ta1%x_d, c_Xh%Binv_d, n)
    call device_col2(ta2%x_d, c_Xh%Binv_d, n)
    call device_col2(ta3%x_d, c_Xh%Binv_d, n)

    call cdtp(wa1%x, ta1%x, c_Xh%drdx, c_Xh%dsdx, c_Xh%dtdx, c_Xh)
    call cdtp(wa2%x, ta2%x, c_Xh%drdy, c_Xh%dsdy, c_Xh%dtdy, c_Xh)
    call cdtp(wa3%x, ta3%x, c_Xh%drdz, c_Xh%dsdz, c_Xh%dtdz, c_Xh)

    call Ax%compute(p_res%x, p%x, c_Xh, p%msh, p%Xh)

#ifdef HAVE_HIP
    call adjoint_prs_res_part2_hip(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n)
#elif HAVE_CUDA
    call adjoint_prs_res_part2_cuda(p_res%x_d, wa1%x_d, wa2%x_d, wa3%x_d, n)
#elif HAVE_OPENCL
    call neko_error('adjoint_prs_res_part2_opencl not implemented')
#else
    call neko_error('No device backend configured')
#endif

    call neko_scratch_registry%relinquish_field(temp_indices)
  end subroutine adjoint_pnpn_prs_res_device_compute

  subroutine adjoint_pnpn_vel_res_device_compute(Ax, u, v, w, u_res, v_res, w_res, &
       p, f_x, f_y, f_z, c_Xh, msh, Xh, mu, rho, bd, dt, n)
    class(ax_t), intent(in) :: Ax
    type(mesh_t), intent(inout) :: msh
    type(space_t), intent(inout) :: Xh
    type(field_t), intent(inout) :: p, u, v, w
    type(field_t), intent(inout) :: u_res, v_res, w_res
    type(field_t), intent(in) :: f_x, f_y, f_z
    type(coef_t), intent(inout) :: c_Xh
    type(field_t), intent(in) :: mu
    type(field_t), intent(in) :: rho
    real(kind=rp), intent(in) :: bd
    real(kind=rp), intent(in) :: dt
    integer, intent(in) :: n
    real(kind=rp) :: mu_val, rho_val

    ! We assume the material properties are constant
    mu_val = mu%x(1,1,1,1)
    rho_val = rho%x(1,1,1,1)

    call device_cfill(c_Xh%h1_d, mu_val, n)
    call device_cfill(c_Xh%h2_d, rho_val * (bd / dt), n)
    c_Xh%ifh2 = .true.

    call Ax%compute_vector(u_res%x, v_res%x, w_res%x, u%x, v%x, w%x, c_Xh, msh, Xh)

#ifdef HAVE_HIP
    call adjoint_vel_res_update_hip(u_res%x_d, v_res%x_d, w_res%x_d, &
         f_x%x_d, f_y%x_d, f_z%x_d, n)
#elif HAVE_CUDA
    call adjoint_vel_res_update_cuda(u_res%x_d, v_res%x_d, w_res%x_d, &
         f_x%x_d, f_y%x_d, f_z%x_d, n)
#elif HAVE_OPENCL
    call neko_error('adjoint_vel_res_update_opencl not implemented')
#else
    call neko_error('No device backend configured')
#endif
  end subroutine adjoint_pnpn_vel_res_device_compute

end module adjoint_pnpn_res_device
