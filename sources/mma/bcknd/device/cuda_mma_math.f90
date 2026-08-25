!> @file cuda_mma_math.f90
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
module cuda_mma_math
  use num_types, only: rp, c_rp
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr

  implicit none
  public

  interface
     subroutine mma_unconstrained_kkt_cuda(rex, x, xmin, xmax, df0dx, eps, n) &
          bind(c, name="mma_unconstrained_kkt_cuda")
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: rex, x, xmin, xmax, df0dx
       real(c_rp) :: eps
       integer(c_int) :: n
     end subroutine mma_unconstrained_kkt_cuda


     subroutine mma_update_hessian_z_cuda(Hess_d, a_d, m) &
          bind(C, name="mma_update_hessian_z_cuda")
       use iso_c_binding
       type(c_ptr), value :: Hess_d
       type(c_ptr), value :: a_d
       integer(c_int) :: m
     end subroutine mma_update_hessian_z_cuda


     subroutine mma_prepare_aa_matrix_cuda(AA, s, lambda, d, mu, y, a, zeta, &
          z, m) bind(c, name="mma_prepare_aa_matrix_cuda")
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: AA, s, lambda, d, mu, y, a
       real(c_rp) :: zeta, z
       integer(c_int) :: m
     end subroutine mma_prepare_aa_matrix_cuda

     subroutine mma_prepare_hessian_cuda(Hess, y, mu, lambda, m) &
          bind(c, name="mma_prepare_hessian_cuda")
       import c_ptr, c_int
       type(c_ptr), value :: Hess, y, mu, lambda
       integer(c_int) :: m
     end subroutine mma_prepare_hessian_cuda

     subroutine cuda_custom_solver(A_d, b_d, n, info) &
          bind(c, name = 'cuda_custom_solver')
       import c_int, c_ptr
       type(c_ptr), value :: A_d, b_d
       integer(c_int), value :: n
       integer(c_int) :: info
     end subroutine cuda_custom_solver

     subroutine cuSOLVER_wrapper(A_d, b_d, n, info) &
          bind(c, name = 'cuSOLVER_wrapper')
       import c_int, c_ptr
       type(c_ptr), value :: A_d, b_d
       integer(c_int), value :: n
       integer(c_int) :: info
     end subroutine cuSOLVER_wrapper

     subroutine delta_1dbeam_cuda(Delta_d, L_total, Le, offset, n) &
          bind(c, name = 'delta_1dbeam_cuda')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: Delta_d
       real(c_rp) :: L_total, Le
       integer(c_int) :: offset, n
     end subroutine delta_1dbeam_cuda

     subroutine cuda_Hess(Hess_d, hijx_d, Ljjxinv_d, n, m) bind(c, name = 'cuda_Hess')
       import c_int, c_ptr
       type(c_ptr), value :: Hess_d, hijx_d, Ljjxinv_d
       integer(c_int) :: n, m
     end subroutine cuda_Hess

     subroutine mma_Ljjxinv_cuda(Ljjxinv_d,pjlambda_d, qjlambda_d, x_d, &
          low_d, upp_d, alpha_d, beta_d, n) bind(c, name = 'mma_Ljjxinv_cuda')
       import c_int, c_ptr
       type(c_ptr), value :: Ljjxinv_d, x_d, pjlambda_d, qjlambda_d, low_d, &
            upp_d, alpha_d, beta_d
       integer(c_int) :: n
     end subroutine mma_Ljjxinv_cuda

     subroutine mma_dipsolvesub1_cuda(x_d, pjlambda_d, qjlambda_d, low_d, &
          upp_d, alpha_d, beta_d, n) bind(c, name = 'mma_dipsolvesub1_cuda')
       import c_int, c_ptr
       type(c_ptr), value :: x_d, pjlambda_d, qjlambda_d, low_d, &
            upp_d, alpha_d, beta_d
       integer(c_int) :: n
     end subroutine mma_dipsolvesub1_cuda

     subroutine mattrans_v_mul_cuda(output_d, pij_d, lambda_d, m, n) &
          bind(c, name = 'mattrans_v_mul_cuda')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: output_d, pij_d, lambda_d
       integer(c_int) :: m, n
     end subroutine mattrans_v_mul_cuda

     subroutine mma_gensub1_cuda(low_d, upp_d, x_d, xmin_d, xmax_d, asyinit, n)&
          bind(c, name = 'mma_gensub1_cuda')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: low_d, upp_d, x_d, xmin_d, xmax_d
       real(c_rp) :: asyinit
       integer(c_int) :: n
     end subroutine mma_gensub1_cuda

     subroutine mma_gensub2_cuda(low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d, &
          asydecr, asyincr, n) bind(c, name = 'mma_gensub2_cuda')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d
       real(c_rp) :: asydecr, asyincr
       integer(c_int) :: n
     end subroutine mma_gensub2_cuda

     subroutine mma_gensub3_cuda(x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, &
          max_d, alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, n, m) &
          bind(c, name = 'mma_gensub3_cuda')
       import c_int, c_ptr
       type(c_ptr), value :: x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, max_d, &
            alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine mma_gensub3_cuda

     subroutine mma_gensub4_cuda(x_d, low_d, upp_d, pij_d, qij_d, n, m, bi_d) &
          bind(c, name = 'mma_gensub4_cuda')
       import c_int, c_ptr
       type(c_ptr), value :: x_d, low_d, upp_d, pij_d, qij_d, bi_d
       integer(c_int) :: n, m
     end subroutine mma_gensub4_cuda

     subroutine cuda_mma_max(xsi_d, x_d, alpha_d, n) &
          bind(c, name = 'cuda_mma_max')
       import c_int, c_ptr
       type(c_ptr), value :: xsi_d, x_d, alpha_d
       integer(c_int) :: n
     end subroutine cuda_mma_max

     subroutine cuda_rex(rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, q0j_d, &
          lambda_d, xsi_d, eta_d, n, m) bind(c, name = 'cuda_rex')
       import c_int, c_ptr
       type(c_ptr), value :: rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, &
            q0j_d, lambda_d, xsi_d, eta_d
       integer(c_int) :: n, m
     end subroutine cuda_rex

     subroutine cuda_relambda(relambda_d, x_d, upp_d, low_d, pij_d, qij_d, n, &
          m) bind(c, name = 'cuda_relambda')
       import c_int, c_ptr
       type(c_ptr), value :: relambda_d, x_d, upp_d, low_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine cuda_relambda

     subroutine cuda_sub2cons2(rexsi_d, xsi_d, x_d, alpha_d, epsi, n) &
          bind(c, name = 'cuda_sub2cons2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: rexsi_d, xsi_d, x_d, alpha_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine cuda_sub2cons2

     real(c_rp) function cuda_maxval(rex_d, n) bind(c, name = 'cuda_maxval')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: rex_d
       integer(c_int) :: n
     end function cuda_maxval

     real(c_rp) function cuda_norm(rex_d, n) bind(c, name = 'cuda_norm')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: rex_d
       integer(c_int) :: n
     end function cuda_norm

     subroutine cuda_delx(delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, &
          q0j_d, alpha_d, beta_d, lambda_d, epsi, n, m) &
          bind(c, name = 'cuda_delx')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, &
            q0j_d, alpha_d, beta_d, lambda_d
       real(c_rp) :: epsi
       integer(c_int) :: n, m
     end subroutine cuda_delx



     subroutine cuda_GG(GG_d, x_d, low_d, upp_d, pij_d, qij_d, n, m) &
          bind(c, name = 'cuda_GG')
       import c_int, c_ptr
       type(c_ptr), value :: GG_d, x_d, low_d, upp_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine cuda_GG

     subroutine cuda_diagx(diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, &
          pij_d, qij_d, alpha_d, beta_d, eta_d, lambda_d, n, m) &
          bind(c, name = 'cuda_diagx')
       import c_int, c_ptr
       type(c_ptr), value :: diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, &
            pij_d, qij_d, alpha_d, beta_d, eta_d, lambda_d
       integer(c_int) :: n, m
     end subroutine cuda_diagx

     subroutine cuda_bb(bb_d, GG_d, delx_d, diagx_d, n, m) &
          bind(c, name = 'cuda_bb')
       import c_int, c_ptr
       type(c_ptr), value :: bb_d, GG_d, delx_d, diagx_d
       integer(c_int) :: n, m
     end subroutine cuda_bb

     subroutine cuda_AA(AA_d, GG_d, diagx_d, n, m) bind(c, name = 'cuda_AA')
       import c_int, c_ptr
       type(c_ptr), value :: AA_d, GG_d, diagx_d
       integer(c_int) :: n, m
     end subroutine cuda_AA

     subroutine cuda_dx(dx_d, delx_d, diagx_d, GG_d, dlambda_d, n, m) &
          bind(c, name = 'cuda_dx')
       import c_int, c_ptr
       type(c_ptr), value :: dx_d, delx_d, diagx_d, GG_d, dlambda_d
       integer(c_int) :: n, m
     end subroutine cuda_dx

     subroutine cuda_dxsi(dxsi_d, xsi_d, dx_d, x_d, alpha_d, epsi, n) &
          bind(c, name = 'cuda_dxsi')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: dxsi_d, xsi_d, dx_d, x_d, alpha_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine cuda_dxsi

     subroutine cuda_deta(deta_d, eta_d, dx_d, x_d, beta_d, epsi, n) &
          bind(c, name = 'cuda_deta')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: deta_d, eta_d, dx_d, x_d, beta_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine cuda_deta

     real(c_rp) function cuda_maxval2(dxx_d, xx_d, cons, n) &
          bind(c, name = 'cuda_maxval2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: dxx_d, xx_d
       real(c_rp) :: cons
       integer(c_int) :: n
     end function cuda_maxval2

     real(c_rp) function cuda_maxval3(dx_d, x_d, alpha_d, cons, n) &
          bind(c, name = 'cuda_maxval3')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: dx_d, x_d, alpha_d
       real(c_rp) :: cons
       integer(c_int) :: n
     end function cuda_maxval3

     subroutine cuda_kkt_rex(rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d, &
          n, m) bind(c, name = 'cuda_kkt_rex')
       import c_int, c_ptr
       type(c_ptr), value :: rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d
       integer(c_int) :: n, m
     end subroutine cuda_kkt_rex


     subroutine cuda_maxcons(a_d, b, c, d_d, n) bind(c, name = 'cuda_maxcons')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, d_d
       real(c_rp) :: b, c
       integer(c_int) :: n
     end subroutine cuda_maxcons


     real(c_rp) function cuda_lcsc2(a_d, b_d, n) bind(c, name = 'cuda_lcsc2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end function cuda_lcsc2

     subroutine cuda_mpisum(a_d, n) bind(c, name = 'cuda_mpisum')
       import c_int, c_ptr
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine cuda_mpisum

     subroutine cuda_add2inv2(a_d, b_d, c, n) bind(c, name = 'cuda_add2inv2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
       real(c_rp) :: c
     end subroutine cuda_add2inv2

     subroutine cuda_max2(a_d, b, c_d, d, n) bind(c, name = 'cuda_max2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, c_d
       integer(c_int) :: n
       real(c_rp) :: b, d
     end subroutine cuda_max2

     subroutine cuda_updatebb(bb_d, dellambda_d, dely_d, d_d, mu_d, y_d, delz, &
          m) bind(c, name = 'cuda_updatebb')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: bb_d, dellambda_d, dely_d, d_d, mu_d, y_d
       integer(c_int) :: m
       real(c_rp) :: delz
     end subroutine cuda_updatebb

     subroutine cuda_updateAA(AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, &
          y_d, a_d, zeta, z, m) bind(c, name = 'cuda_updateAA')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, &
            y_d, a_d
       integer(c_int) :: m
       real(c_rp) :: zeta, z
     end subroutine cuda_updateAA

     subroutine cuda_dy(dy_d, dely_d, dlambda_d, d_d, mu_d, y_d, n) &
          bind(c, name = 'cuda_dy')
       import c_int, c_ptr
       type(c_ptr), value :: dy_d, dely_d, dlambda_d, d_d, mu_d, y_d
       integer(c_int) :: n
     end subroutine cuda_dy

  end interface

end module cuda_mma_math
