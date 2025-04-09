! Copyright (c) 2025, The Neko-TOP Authors
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions
! are met:
!
!   * Redistributions of source code must retain the above copyright
!     notice, this list of conditions and the following disclaimer.
!
!   * Redistributions in binary form must reproduce the above
!     copyright notice, this list of conditions and the following
!     disclaimer in the documentation and/or other materials provided
!     with the distribution.
!
!   * Neither the name of the authors nor the names of its
!     contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
! FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
! COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
! INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
! ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.
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
       mma_gensub3_cuda, mma_gensub4_cuda
  use hip_mma_math, only: hip_mma_max, hip_max2, hip_rex, hip_lcsc2, &
       hip_relambda, hip_sub2cons2, hip_maxval, hip_norm, hip_delx, &
       hip_add2inv2, hip_GG, hip_diagx, hip_bb, hip_updatebb, hip_AA, &
       hip_updateAA, hip_dx, hip_dy, hip_deta, hip_dxsi, hip_maxval2, &
       hip_maxval3, hip_kkt_rex, hip_mma_gensub1, hip_mma_gensub2, &
       hip_mma_gensub3, hip_mma_gensub4

  implicit none
  private

! #ifdef HAVE_CUDA
! #ifdef HAVE_HIP

  public :: device_mma_gensub1, device_mma_gensub2, device_mma_gensub3, &
       device_mma_gensub4, device_mma_max, device_max2, device_rex, &
       device_lcsc2, device_relambda, device_sub2cons2, device_maxval, &
       device_norm, device_delx, device_add2inv2, device_GG, device_diagx, &
       device_bb, device_updatebb, device_AA, device_updateAA, device_dx, &
       device_dy, device_deta, device_dxsi, device_maxval2, device_maxval3, &
       device_kkt_rex

contains

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

  subroutine device_mma_gensub2(low_d, upp_d, x_d, xold1_d, xold2_d, xmin_d, &
       xmax_d, asydecr, asyincr, n)
    type(c_ptr) :: low_d, upp_d, x_d, xold1_d, xold2_d, xmin_d, xmax_d
    real(c_rp) :: asydecr, asyincr
    integer :: n
#if HAVE_HIP
    call mma_gensub2_hip(low_d, upp_d, x_d, xold1_d, xold2_d, xmin_d, xmax_d, &
         asydecr, asyincr, n)
#elif HAVE_CUDA
    call mma_gensub2_cuda(low_d, upp_d, x_d, xold1_d, xold2_d, xmin_d, xmax_d, &
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

  subroutine device_rex(x_d, low_d, upp_d, lambda_d, n, m)
    type(c_ptr) :: x_d, low_d, upp_d, lambda_d
    integer :: n, m
#if HAVE_HIP
    call hip_rex(x_d, low_d, upp_d, lambda_d, n, m)
#elif HAVE_CUDA
    call cuda_rex(x_d, low_d, upp_d, lambda_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_rex

  subroutine device_lcsc2(x_d, low_d, upp_d, n)
    type(c_ptr) :: x_d, low_d, upp_d
    integer :: n
#if HAVE_HIP
    call hip_lcsc2(x_d, low_d, upp_d, n)
#elif HAVE_CUDA
    call cuda_lcsc2(x_d, low_d, upp_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_lcsc2

  subroutine device_relambda(x_d, low_d, upp_d, lambda_d, n)
    type(c_ptr) :: x_d, low_d, upp_d, lambda_d
    integer :: n
#if HAVE_HIP
    call hip_relambda(x_d, low_d, upp_d, lambda_d, n)
#elif HAVE_CUDA
    call cuda_relambda(x_d, low_d, upp_d, lambda_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_relambda

  subroutine device_sub2cons2(x_d, low_d, upp_d, n)
    type(c_ptr) :: x_d, low_d, upp_d
    integer :: n
#if HAVE_HIP
    call hip_sub2cons2(x_d, low_d, upp_d, n)
#elif HAVE_CUDA
    call cuda_sub2cons2(x_d, low_d, upp_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_sub2cons2

  subroutine device_maxval(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_maxval(x_d, n)
#elif HAVE_CUDA
    call cuda_maxval(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_maxval

  subroutine device_norm(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_norm(x_d, n)
#elif HAVE_CUDA
    call cuda_norm(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_norm

  subroutine device_delx(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_delx(x_d, n)
#elif HAVE_CUDA
    call cuda_delx(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_delx

  subroutine device_add2inv2(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_add2inv2(x_d, n)
#elif HAVE_CUDA
    call cuda_add2inv2(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_add2inv2

  subroutine device_GG(x_d, y_d, n)
    type(c_ptr) :: x_d, y_d
    integer :: n
#if HAVE_HIP
    call hip_GG(x_d, y_d, n)
#elif HAVE_CUDA
    call cuda_GG(x_d, y_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_GG

  subroutine device_diagx(x_d, y_d, n)
    type(c_ptr) :: x_d, y_d
    integer :: n
#if HAVE_HIP
    call hip_diagx(x_d, y_d, n)
#elif HAVE_CUDA
    call cuda_diagx(x_d, y_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_diagx

  subroutine device_bb(x_d, y_d, n)
    type(c_ptr) :: x_d, y_d
    integer :: n
#if HAVE_HIP
    call hip_bb(x_d, y_d, n)
#elif HAVE_CUDA
    call cuda_bb(x_d, y_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_bb

  subroutine device_updatebb(x_d, y_d, n)
    type(c_ptr) :: x_d, y_d
    integer :: n
#if HAVE_HIP
    call hip_updatebb(x_d, y_d, n)
#elif HAVE_CUDA
    call cuda_updatebb(x_d, y_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_updatebb

  subroutine device_AA(x_d, y_d, n)
    type(c_ptr) :: x_d, y_d
    integer :: n
#if HAVE_HIP
    call hip_AA(x_d, y_d, n)
#elif HAVE_CUDA
    call cuda_AA(x_d, y_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_AA

  subroutine device_updateAA(x_d, y_d, n)
    type(c_ptr) :: x_d, y_d
    integer :: n
#if HAVE_HIP
    call hip_updateAA(x_d, y_d, n)
#elif HAVE_CUDA
    call cuda_updateAA(x_d, y_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_updateAA

  subroutine device_dx(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_dx(x_d, n)
#elif HAVE_CUDA
    call cuda_dx(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_dx

  subroutine device_dy(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_dy(x_d, n)
#elif HAVE_CUDA
    call cuda_dy(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_dy

  subroutine device_deta(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_deta(x_d, n)
#elif HAVE_CUDA
    call cuda_deta(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_deta

  subroutine device_dxsi(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_dxsi(x_d, n)
#elif HAVE_CUDA
    call cuda_dxsi(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_dxsi

  subroutine device_maxval2(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_maxval2(x_d, n)
#elif HAVE_CUDA
    call cuda_maxval2(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_maxval2

  subroutine device_maxval3(x_d, n)
    type(c_ptr) :: x_d
    integer :: n
#if HAVE_HIP
    call hip_maxval3(x_d, n)
#elif HAVE_CUDA
    call cuda_maxval3(x_d, n)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_maxval3

  subroutine device_kkt_rex(x_d, low_d, upp_d, lambda_d, n, m)
    type(c_ptr) :: x_d, low_d, upp_d, lambda_d
    integer :: n, m
#if HAVE_HIP
    call hip_kkt_rex(x_d, low_d, upp_d, lambda_d, n, m)
#elif HAVE_CUDA
    call cuda_kkt_rex(x_d, low_d, upp_d, lambda_d, n, m)
#elif HAVE_OPENCL
    call neko_error('no device backend configured')
#else
    call neko_error('no device backend configured')
#endif
  end subroutine device_kkt_rex

end module device_mma_math
