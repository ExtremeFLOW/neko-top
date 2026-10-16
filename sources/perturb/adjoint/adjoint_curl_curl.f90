!> @file adjoint_curl_curl.f90
!! Exact adjoint of the curl-curl term in the Pn/Pn pressure residual.
module adjoint_curl_curl
  use num_types, only: rp
  use field, only: field_t
  use coefs, only: coef_t
  use gather_scatter, only: gs_t
  use gs_ops, only: GS_OP_ADD
  use operators, only: opgrad, cdtp
  use field_math, only: field_sub2, field_cmult
  use neko_config, only: NEKO_BCKND_DEVICE
  use math, only: col2
  use device_math, only: device_col2
  implicit none
  private

  public :: adjoint_curl_curl_load

contains

  !> Compute L = (mu/rho) Dc^T M^T Dc^T M^T G p_adj.
  !! The load is element-local, not gather-scattered, because it is added to
  !! the adjoint forcing before the velocity residual gather-scatter.
  subroutine adjoint_curl_curl_load(load_x, load_y, load_z, p_adj, coef, gs, &
       mu, rho, work1, work2, work3, work4)
    type(field_t), intent(inout) :: load_x, load_y, load_z
    type(field_t), intent(in) :: p_adj
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    real(kind=rp), intent(in) :: mu, rho
    type(field_t), intent(inout) :: work1, work2, work3, work4

    call opgrad(load_x%x, load_y%x, load_z%x, p_adj%x, coef)
    call adjoint_curl_curl_gs_binv(load_x, load_y, load_z, coef, gs)

    call adjoint_curl_transpose(work1, work2, work3, load_x, load_y, load_z, &
         work4, coef)
    call adjoint_curl_curl_gs_binv(work1, work2, work3, coef, gs)

    call adjoint_curl_transpose(load_x, load_y, load_z, work1, work2, work3, &
         work4, coef)
    call field_cmult(load_x, mu / rho)
    call field_cmult(load_y, mu / rho)
    call field_cmult(load_z, mu / rho)

  end subroutine adjoint_curl_curl_load

  subroutine adjoint_curl_transpose(w1, w2, w3, y1, y2, y3, work, coef)
    type(field_t), intent(inout) :: w1, w2, w3
    type(field_t), intent(inout) :: y1, y2, y3
    type(field_t), intent(inout) :: work
    type(coef_t), intent(in) :: coef

    call cdtp(w1%x, y2%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call cdtp(work%x, y3%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call field_sub2(w1, work)

    call cdtp(w2%x, y3%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call cdtp(work%x, y1%x, coef%drdz, coef%dsdz, coef%dtdz, coef)
    call field_sub2(w2, work)

    call cdtp(w3%x, y1%x, coef%drdy, coef%dsdy, coef%dtdy, coef)
    call cdtp(work%x, y2%x, coef%drdx, coef%dsdx, coef%dtdx, coef)
    call field_sub2(w3, work)

  end subroutine adjoint_curl_transpose

  subroutine adjoint_curl_curl_gs_binv(y1, y2, y3, coef, gs)
    type(field_t), intent(inout) :: y1, y2, y3
    type(coef_t), intent(in) :: coef
    type(gs_t), intent(inout) :: gs
    integer :: n

    n = y1%size()

    call gs%op(y1, GS_OP_ADD)
    call gs%op(y2, GS_OP_ADD)
    call gs%op(y3, GS_OP_ADD)

    if (NEKO_BCKND_DEVICE .eq. 1) then
       call device_col2(y1%x_d, coef%Binv_d, n)
       call device_col2(y2%x_d, coef%Binv_d, n)
       call device_col2(y3%x_d, coef%Binv_d, n)
    else
       call col2(y1%x, coef%Binv, n)
       call col2(y2%x, coef%Binv, n)
       call col2(y3%x, coef%Binv, n)
    end if

  end subroutine adjoint_curl_curl_gs_binv

end module adjoint_curl_curl
