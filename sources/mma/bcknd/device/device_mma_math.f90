!> @file device_mma_math.f90
!! @copyright
!! Copyright (c) 2025, The Neko-TOP Authors
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!!   * Redistributions of source code must retain the above copyright
!!     notice, this list of conditions and the following disclaimer.
!!
!!   * Redistributions in binary form must reproduce the above
!!     copyright notice, this list of conditions and the following
!!     disclaimer in the documentation and/or other materials provided
!!     with the distribution.
!!
!!   * Neither the name of the authors nor the names of its
!!     contributors may be used to endorse or promote products derived
!!     from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
!! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
!! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
!! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
!! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!
module device_mma_math
  use, intrinsic :: iso_c_binding, only: c_ptr, c_int
  use num_types, only: rp, c_rp
  use utils, only: neko_error
  use comm, only: NEKO_COMM, pe_size, MPI_REAL_PRECISION
  use mpi_f08, only: MPI_SUM, MPI_IN_PLACE, MPI_Allreduce
  use cuda_mma_math, only: cuda_mma_max, cuda_max2, cuda_rex, cuda_lcsc2, &
       cuda_relambda, cuda_sub2cons2, cuda_maxval, cuda_norm, cuda_delx, &
       cuda_add2inv2, cuda_GG, cuda_diagx, cuda_bb, cuda_updatebb, cuda_AA, &
       cuda_updateAA, cuda_dx, cuda_dy, cuda_deta, cuda_dxsi, cuda_maxval2, &
       cuda_maxval3, cuda_kkt_rex, mma_gensub1_cuda, mma_gensub2_cuda, &
       mma_gensub3_cuda, mma_gensub4_cuda, mattrans_v_mul_cuda, &
       mma_dipsolvesub1_cuda, mma_Ljjxinv_cuda, cuda_Hess, delta_1dbeam_cuda, &
       cuSOLVER_wrapper, mma_prepare_hessian_cuda, mma_prepare_aa_matrix_cuda, &
       cuda_custom_solver, mma_update_hessian_z_cuda, mma_unconstrained_kkt_cuda
  use hip_mma_math, only: hip_mma_max, hip_max2, hip_rex, hip_lcsc2, &
       hip_relambda, hip_sub2cons2, hip_maxval, hip_norm, hip_delx, &
       hip_add2inv2, hip_GG, hip_diagx, hip_bb, hip_updatebb, hip_AA, &
       hip_updateAA, hip_dx, hip_dy, hip_deta, hip_dxsi, hip_maxval2, &
       hip_maxval3, hip_kkt_rex, mma_gensub1_hip, mma_gensub2_hip, &
       mma_gensub3_hip, mma_gensub4_hip, mattrans_v_mul_hip, &
       mma_dipsolvesub1_hip, mma_Ljjxinv_hip, hip_Hess, delta_1dbeam_hip, &
       hip_custom_solver, mma_prepare_hessian_hip, &
       mma_prepare_aa_matrix_hip, hipSOLVER_wrapper, mma_update_hessian_z_hip!, &
      ! mma_unconstrained_kkt_hip

  implicit none
  private


  public :: device_mma_gensub1, device_mma_gensub2, device_mma_gensub3, &
       device_mma_gensub4, device_mma_max, device_max2, device_rex, &
       device_lcsc2, device_relambda, device_sub2cons2, device_maxval, &
       device_norm, device_delx, device_add2inv2, device_GG, device_diagx, &
       device_bb, device_updatebb, device_AA, device_updateAA, device_dx, &
       device_dy, device_deta, device_dxsi, device_maxval2, device_maxval3, &
       device_kkt_rex, device_mattrans_v_mul, device_mma_dipsolvesub1, &
       device_mma_Ljjxinv, device_Hess, device_delta_1dbeam, &
       device_solve_linear_system, device_prepare_hessian, &
       device_prepare_aa_matrix, device_update_hessian_z, device_unconstrained_kkt

contains
  !> Compute the stationarity residual projected onto the box bounds on device
  subroutine device_unconstrained_kkt(rex_d, x_d, xmin_d, xmax_d, df0dx_d, eps, n)
    type(c_ptr), intent(in) :: rex_d, x_d, xmin_d, xmax_d, df0dx_d
    real(c_rp), intent(in) :: eps
    integer, value :: n
#if HAVE_HIP
    !call mma_unconstrained_kkt_hip(rex_d, x_d, xmin_d, xmax_d, df0dx_d, eps, n)
#elif HAVE_CUDA
    call mma_unconstrained_kkt_cuda(rex_d, x_d, xmin_d, xmax_d, df0dx_d, eps, n)
#elif HAVE_OPENCL
    call neko_error('Unconstrained KKT not implemented for OpenCL')
#else
    call neko_error('No device backend configured for Unconstrained KKT')
#endif
  end subroutine device_unconstrained_kkt

  !> Update Hessian for dual solver with z-term contribution: Hess -= a * a^T
  subroutine device_update_hessian_z(Hess_d, a_d, m)
    use iso_c_binding
    type(c_ptr), intent(in) :: Hess_d
    type(c_ptr), intent(in) :: a_d
    integer, value :: m
#if HAVE_HIP
    call mma_update_hessian_z_hip(Hess_d, a_d, m)
#elif HAVE_CUDA
    call mma_update_hessian_z_cuda(Hess_d, a_d, m)
#elif HAVE_OPENCL
    call neko_error('Z-term Hessian update not implemented for OpenCL')
#else
    call neko_error('No device backend configured for Z-term Hessian update')
#endif
  end subroutine device_update_hessian_z


  !> Prepare AA matrix for dual-primal solver on device
  subroutine device_prepare_aa_matrix(AA_d, s_d, lambda_d, d_d, mu_d, y_d, &
       a_d, zeta, z, m)
    type(c_ptr) :: AA_d, s_d, lambda_d, d_d, mu_d, y_d, a_d
    real(c_rp) :: zeta, z
    integer, value :: m
#if HAVE_HIP
    call mma_prepare_aa_matrix_hip(AA_d, s_d, lambda_d, d_d, mu_d, y_d, a_d, &
         zeta, z, m)
#elif HAVE_CUDA
    call mma_prepare_aa_matrix_cuda(AA_d, s_d, lambda_d, d_d, mu_d, y_d, a_d, &
         zeta, z, m)
#elif HAVE_OPENCL
    call neko_error('AA matrix preparation not implemented for OpenCL')
#else
    call neko_error('No device backend configured for AA matrix preparation')
#endif
  end subroutine device_prepare_aa_matrix

  !> Solve linear system Ax = b on device
  subroutine device_prepare_hessian(Hess_d, y_d, mu_d, lambda_d, m)
    type(c_ptr) :: Hess_d, y_d, mu_d, lambda_d
    integer, value :: m
#if HAVE_HIP
    call mma_prepare_hessian_hip(Hess_d, y_d, mu_d, lambda_d, m)
#elif HAVE_CUDA
    call mma_prepare_hessian_cuda(Hess_d, y_d, mu_d, lambda_d, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_prepare_hessian

  subroutine device_solve_linear_system(A_d, b_d, n, info)
    type(c_ptr) :: A_d, b_d
    integer(c_int), value :: n
    integer(c_int) :: info
#if HAVE_HIP
    call hipSOLVER_wrapper(A_d, b_d, n, info)
    ! call hip_custom_solver(A_d, b_d, n, info)
#elif HAVE_CUDA
    call cuSOLVER_wrapper(A_d, b_d, n, info)
    ! call cuda_custom_solver(A_d, b_d, n, info)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_solve_linear_system


  !> A device support to do the following calculation for 1D beam elements:
  !!   Delta(k) = ((L_total - Le*(offset+k-1))**3 - &
  !!              (L_total - Le*(offset+k))**3) / 3.0_rp
  !! Where k ranges from 1 to n
  subroutine device_delta_1dbeam(Delta_d, L_total, Le, offset, n)
    type(c_ptr) :: Delta_d
    real(c_rp) :: L_total, Le
    integer(c_int), value :: offset, n
#if HAVE_HIP
    call delta_1dbeam_hip(Delta_d, L_total, Le, offset, n)
#elif HAVE_CUDA
    call delta_1dbeam_cuda(Delta_d, L_total, Le, offset, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_delta_1dbeam


  subroutine device_Hess(Hess_d, hijx_d, Ljjxinv_d, n, m)
    type(c_ptr):: Hess_d, hijx_d, Ljjxinv_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_Hess(Hess_d, hijx_d, Ljjxinv_d, n, m)
#elif HAVE_CUDA
    call cuda_Hess(Hess_d, hijx_d, Ljjxinv_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_Hess


  !> A device support to do the following calculation needed for the
  !! dualsubsolve for MMA:
  !!   Ljjxinv= - 1 / ( (2*pjlambda/(upp - x)**3) + &
  !!                                  (2*qjlambda/(x - low)**3))
  !!
  !! And then remove the sensitivity for the active primal constraints
  !! Ljjxinv = merge(0.0_rp, Ljjxinv, x .eq. alpha)
  !! Ljjxinv = merge(0.0_rp, Ljjxinv, x .eq. beta)
  subroutine device_mma_Ljjxinv(Ljjxinv_d,pjlambda_d, qjlambda_d, x_d, &
       low_d, upp_d, alpha_d, beta_d, n)
    type(c_ptr) :: Ljjxinv_d, pjlambda_d, qjlambda_d, x_d, &
         low_d, upp_d, alpha_d, beta_d
    integer(c_int) :: n
#if HAVE_HIP
    call mma_Ljjxinv_hip(Ljjxinv_d, pjlambda_d, qjlambda_d, x_d, &
         low_d, upp_d, alpha_d, beta_d, n)
#elif HAVE_CUDA
    call mma_Ljjxinv_cuda(Ljjxinv_d, pjlambda_d, qjlambda_d, x_d, &
         low_d, upp_d, alpha_d, beta_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mma_Ljjxinv


  !> A device support to do the following calculation needed for the
  !! dualsubsolve for MMA:
  !!   x = (sqrt(pjlambda) * low + sqrt(qjlambda) * upp) /  &
  !!                                  (sqrt(pjlambda) + sqrt(qjlambda))
  subroutine device_mma_dipsolvesub1(x_d, pjlambda_d, qjlambda_d, &
       low_d, upp_d, alpha_d, beta_d, n)
    type(c_ptr) :: x_d, pjlambda_d, qjlambda_d, &
         low_d, upp_d, alpha_d, beta_d
    integer(c_int) :: n
#if HAVE_HIP
    call mma_dipsolvesub1_hip(x_d, pjlambda_d, qjlambda_d, &
         low_d, upp_d, alpha_d, beta_d, n)
#elif HAVE_CUDA
    call mma_dipsolvesub1_cuda(x_d, pjlambda_d, qjlambda_d, &
         low_d, upp_d, alpha_d, beta_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mma_dipsolvesub1


  !> A device support to do the following matrix multiplication
  !!               output = matmul(transpose(pij), lambda)
  !! where matrix pij is mxn, vector lambda is of size m and the output
  !! vector is of size n
  subroutine device_mattrans_v_mul(output_d, pij_d, lambda_d, m, n)
    type(c_ptr) :: output_d, pij_d, lambda_d
    integer :: m, n
#if HAVE_HIP
    call mattrans_v_mul_hip(output_d, pij_d, lambda_d, m, n)
#elif HAVE_CUDA
    call mattrans_v_mul_cuda(output_d, pij_d, lambda_d, m, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mattrans_v_mul



  subroutine device_mma_gensub1(low_d, upp_d, x_d, xmin_d, xmax_d, asyinit, n)
    type(c_ptr) :: low_d, upp_d, x_d, xmin_d, xmax_d
    real(c_rp) :: asyinit
    integer :: n
#if HAVE_HIP
    call mma_gensub1_hip(low_d, upp_d, x_d, xmin_d, xmax_d, asyinit, n)
#elif HAVE_CUDA
    call mma_gensub1_cuda(low_d, upp_d, x_d, xmin_d, xmax_d, asyinit, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mma_gensub1


  subroutine device_mma_gensub2(low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d, &
       asydecr, asyincr, n)
    type(c_ptr) :: low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d
    real(c_rp) :: asydecr, asyincr
    integer :: n
#if HAVE_HIP
    call mma_gensub2_hip(low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d, &
         asydecr, asyincr, n)
#elif HAVE_CUDA
    call mma_gensub2_cuda(low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d, &
         asydecr, asyincr, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mma_gensub2

  subroutine device_mma_gensub3(x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, &
       max_d, alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, n, m)
    type(c_ptr) :: x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, max_d, &
         alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d
    integer :: n, m
#if HAVE_HIP
    call mma_gensub3_hip(x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, max_d, &
         alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, n, m)
#elif HAVE_CUDA
    call mma_gensub3_cuda(x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, max_d, &
         alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured2')
#else
    call neko_error('no device backend configured3')
#endif
  end subroutine device_mma_gensub3

  subroutine device_mma_gensub4(x_d, low_d, upp_d, pij_d, qij_d, n, m, bi_d)
    type(c_ptr) :: x_d, low_d, upp_d, pij_d, qij_d, bi_d
    integer :: n, m
#if HAVE_HIP
    call mma_gensub4_hip(x_d, low_d, upp_d, pij_d, qij_d, n, m, bi_d)
#elif HAVE_CUDA
    call mma_gensub4_cuda(x_d, low_d, upp_d, pij_d, qij_d, n, m, bi_d)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mma_gensub4

  subroutine device_mma_max(xsi_d, x_d, alpha_d, n)
    type(c_ptr) :: xsi_d, x_d, alpha_d
    integer :: n
#if HAVE_HIP
    call hip_mma_max(xsi_d, x_d, alpha_d, n)
#elif HAVE_CUDA
    call cuda_mma_max(xsi_d, x_d, alpha_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_mma_max

  subroutine device_max2(a_d, b, c_d, d, n)
    type(c_ptr) :: a_d, c_d
    real(c_rp) :: b, d
    integer :: n
#if HAVE_HIP
    call hip_max2(a_d, b, c_d, d, n)
#elif HAVE_CUDA
    call cuda_max2(a_d, b, c_d, d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_max2



  subroutine device_rex(rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, &
       q0j_d, lambda_d, xsi_d, eta_d, n, m)
    type(c_ptr) :: rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, q0j_d, &
         lambda_d, xsi_d, eta_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_rex(rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, q0j_d, &
         lambda_d, xsi_d, eta_d, n, m)
#elif HAVE_CUDA
    call cuda_rex(rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, q0j_d, &
         lambda_d, xsi_d, eta_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_rex

  function device_lcsc2(a_d, b_d, n) result(res)
    type(c_ptr) :: a_d, b_d
    integer(c_int) :: n
    real(kind=rp) :: res
    ! Default value in case of no valid backend (resolve compiler warning)
    res = 0.0_rp
#if HAVE_HIP
    res = hip_lcsc2(a_d, b_d, n)
#elif HAVE_CUDA
    res = cuda_lcsc2(a_d, b_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end function device_lcsc2

  subroutine device_relambda(relambda_d, x_d, upp_d, low_d, pij_d, qij_d, n, m)
    type(c_ptr) :: relambda_d, x_d, upp_d, low_d, pij_d, qij_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_relambda(relambda_d, x_d, upp_d, low_d, pij_d, qij_d, n, m)
#elif HAVE_CUDA
    call cuda_relambda(relambda_d, x_d, upp_d, low_d, pij_d, qij_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_relambda

  subroutine device_sub2cons2(rexsi_d, xsi_d, x_d, alpha_d, epsi, n)
    type(c_ptr):: rexsi_d, xsi_d, x_d, alpha_d
    real(kind=rp) :: epsi
    integer(c_int) :: n
#if HAVE_HIP
    call hip_sub2cons2(rexsi_d, xsi_d, x_d, alpha_d, epsi, n)
#elif HAVE_CUDA
    call cuda_sub2cons2(rexsi_d, xsi_d, x_d, alpha_d, epsi, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_sub2cons2

  function device_maxval(rex_d, n) result(res)
    type(c_ptr):: rex_d
    real(kind=rp) :: res
    integer(c_int) :: n
    ! Default value in case of no valid backend (resolve compiler warning)
    res = 0.0_rp
#if HAVE_HIP
    res = hip_maxval(rex_d, n)
#elif HAVE_CUDA
    res = cuda_maxval(rex_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end function device_maxval

  function device_norm(rex_d, n) result(res)
    type(c_ptr):: rex_d
    real(kind=rp) :: res
    integer(c_int) :: n
    ! Default value in case of no valid backend (resolve compiler warning)
    res = 0.0_rp
#if HAVE_HIP
    res = hip_norm(rex_d, n)
#elif HAVE_CUDA
    res = cuda_norm(rex_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end function device_norm

  subroutine device_delx(delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, &
       q0j_d, alpha_d, beta_d, lambda_d, epsi, n, m)
    type(c_ptr):: delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, q0j_d, &
         alpha_d, beta_d, lambda_d
    real(kind=rp) :: epsi
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_delx(delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, q0j_d, &
         alpha_d, beta_d, lambda_d, epsi, n, m)
#elif HAVE_CUDA
    call cuda_delx(delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, q0j_d, &
         alpha_d, beta_d, lambda_d, epsi, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_delx

  subroutine device_add2inv2(a_d, b_d, c, n)
    type(c_ptr):: a_d, b_d
    real(kind=rp) :: c
    integer(c_int) :: n
#if HAVE_HIP
    call hip_add2inv2(a_d, b_d, c, n)
#elif HAVE_CUDA
    call cuda_add2inv2(a_d, b_d, c, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_add2inv2


  subroutine device_GG(GG_d, x_d, low_d, upp_d, pij_d, qij_d, n, m)
    type(c_ptr):: GG_d, x_d, low_d, upp_d, pij_d, qij_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_GG(GG_d, x_d, low_d, upp_d, pij_d, qij_d, n, m)
#elif HAVE_CUDA
    call cuda_GG(GG_d, x_d, low_d, upp_d, pij_d, qij_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_GG

  subroutine device_diagx(diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, &
       pij_d, qij_d, alpha_d, beta_d, eta_d, lambda_d, n, m)
    type(c_ptr):: diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, pij_d, &
         qij_d, alpha_d, &
         beta_d, eta_d, lambda_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_diagx(diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, pij_d, &
         qij_d, alpha_d, beta_d, eta_d, lambda_d, n, m)
#elif HAVE_CUDA
    call cuda_diagx(diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, pij_d, &
         qij_d, alpha_d, beta_d, eta_d, lambda_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_diagx


  subroutine device_bb(bb_d, GG_d, delx_d, diagx_d, n, m)
    type(c_ptr):: bb_d, GG_d, delx_d, diagx_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_bb(bb_d, GG_d, delx_d, diagx_d, n, m)
#elif HAVE_CUDA
    call cuda_bb(bb_d, GG_d, delx_d, diagx_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_bb

  subroutine device_updatebb(bb_d, dellambda_d, dely_d, d_d, mu_d, y_d, delz, m)
    type(c_ptr):: bb_d, dellambda_d, dely_d, d_d, mu_d, y_d
    integer(c_int) :: m
    real(c_rp) :: delz
#if HAVE_HIP
    call hip_updatebb(bb_d, dellambda_d, dely_d, d_d, mu_d, y_d, delz, m)
#elif HAVE_CUDA
    call cuda_updatebb(bb_d, dellambda_d, dely_d, d_d, mu_d, y_d, delz, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_updatebb

  subroutine device_AA(AA_d, GG_d, diagx_d, n, m)
    type(c_ptr):: AA_d, GG_d, diagx_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_AA(AA_d, GG_d, diagx_d, n, m)
#elif HAVE_CUDA
    call cuda_AA(AA_d, GG_d, diagx_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_AA

  subroutine device_updateAA(AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, &
       y_d, a_d, zeta, z, m)
    type(c_ptr):: AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, y_d, a_d
    integer(c_int) :: m
    real(c_rp) :: zeta, z
#if HAVE_HIP
    call hip_updateAA(AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, y_d, &
         a_d, zeta, z, m)
#elif HAVE_CUDA
    call cuda_updateAA(AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, y_d, &
         a_d, zeta, z, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_updateAA

  subroutine device_dx(dx_d, delx_d, diagx_d, GG_d, dlambda_d, n, m)
    type(c_ptr):: dx_d, delx_d, diagx_d, GG_d, dlambda_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_dx(dx_d, delx_d, diagx_d, GG_d, dlambda_d, n, m)
#elif HAVE_CUDA
    call cuda_dx(dx_d, delx_d, diagx_d, GG_d, dlambda_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_dx

  subroutine device_dy(dy_d, dely_d, dlambda_d, d_d, mu_d, y_d, n)
    type(c_ptr):: dy_d, dely_d, dlambda_d, d_d, mu_d, y_d
    integer(c_int) :: n
#if HAVE_HIP
    call hip_dy(dy_d, dely_d, dlambda_d, d_d, mu_d, y_d, n)
#elif HAVE_CUDA
    call cuda_dy(dy_d, dely_d, dlambda_d, d_d, mu_d, y_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_dy

  subroutine device_dxsi(dxsi_d, xsi_d, dx_d, x_d, alpha_d, epsi, n)
    type(c_ptr):: dxsi_d, xsi_d, dx_d, x_d, alpha_d
    integer(c_int) :: n
    real(c_rp) :: epsi
#if HAVE_HIP
    call hip_dxsi(dxsi_d, xsi_d, dx_d, x_d, alpha_d, epsi, n)
#elif HAVE_CUDA
    call cuda_dxsi(dxsi_d, xsi_d, dx_d, x_d, alpha_d, epsi, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_dxsi

  subroutine device_deta(deta_d, eta_d, dx_d, x_d, beta_d, epsi, n)
    type(c_ptr):: deta_d, eta_d, dx_d, x_d, beta_d
    integer(c_int) :: n
    real(c_rp) :: epsi
#if HAVE_HIP
    call hip_deta(deta_d, eta_d, dx_d, x_d, beta_d, epsi, n)
#elif HAVE_CUDA
    call cuda_deta(deta_d, eta_d, dx_d, x_d, beta_d, epsi, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_deta

  function device_maxval2(dxx_d, xx_d, cons, n) result(res)
    type(c_ptr):: dxx_d, xx_d
    integer :: n
    real(kind=rp), intent(in) :: cons
    real(kind=rp) :: res
    ! Default value in case of no valid backend (resolve compiler warning)
    res = 0.0_rp
#if HAVE_HIP
    res = hip_maxval2(dxx_d, xx_d, cons, n)
#elif HAVE_CUDA
    res = cuda_maxval2(dxx_d, xx_d, cons, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end function device_maxval2


  function device_maxval3(dx_d, x_d, alpha_d, cons, n) result(res)
    type(c_ptr):: dx_d, x_d, alpha_d
    real(kind=rp), intent(in) :: cons
    real(kind=rp) :: res
    integer(c_int) :: n
    ! Default value in case of no valid backend (resolve compiler warning)
    res = 0.0_rp
#if HAVE_HIP
    res = hip_maxval3(dx_d, x_d, alpha_d, cons, n)
#elif HAVE_CUDA
    res = cuda_maxval3(dx_d, x_d, alpha_d, cons, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end function device_maxval3



  subroutine device_kkt_rex(rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d, &
       n, m)
    type(c_ptr):: rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d
    integer(c_int) :: n, m
#if HAVE_HIP
    call hip_kkt_rex(rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d, n, m)
#elif HAVE_CUDA
    call cuda_kkt_rex(rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_kkt_rex

! #endif

end module device_mma_math
