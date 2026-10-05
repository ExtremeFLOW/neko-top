!> @file hip_mma_math.f90
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
module hip_mma_math
  use num_types, only: rp, c_rp
  use, intrinsic :: iso_c_binding, only: c_int, c_ptr

  implicit none
  public

  interface
     subroutine mma_unconstrained_kkt_hip(rex, x, xmin, xmax, df0dx, eps, n) &
          bind(c, name="mma_unconstrained_kkt_hip")
       import c_rp, c_ptr, c_int
       type(c_ptr), value :: rex, x, xmin, xmax, df0dx
       real(c_rp) :: eps
       integer(c_int) :: n
     end subroutine mma_unconstrained_kkt_hip

     subroutine mma_update_hessian_z_hip(Hess_d, a_d, m) &
          bind(C, name="mma_update_hessian_z_hip")
       use iso_c_binding
       type(c_ptr), value :: Hess_d
       type(c_ptr), value :: a_d
       integer(c_int), value :: m
     end subroutine mma_update_hessian_z_hip

     subroutine hipSOLVER_wrapper(A_d, b_d, n, info) &
          bind(c, name = 'hipSOLVER_wrapper')
       import c_int, c_ptr
       type(c_ptr), value :: A_d, b_d
       integer(c_int), value :: n
       integer(c_int) :: info
     end subroutine hipSOLVER_wrapper

     subroutine hip_custom_solver(A_d, b_d, n, info) &
          bind(c, name = 'hip_custom_solver')
       import c_int, c_ptr
       type(c_ptr), value :: A_d, b_d
       integer(c_int), value :: n
       integer(c_int) :: info
     end subroutine hip_custom_solver

     subroutine mma_prepare_hessian_hip(Hess_d, y_d, mu_d, lambda_d, m) &
          bind(c, name = 'mma_prepare_hessian_hip')
       import c_int, c_ptr
       type(c_ptr), value :: Hess_d, y_d, mu_d, lambda_d
       integer(c_int), value :: m
     end subroutine mma_prepare_hessian_hip

     subroutine mma_prepare_aa_matrix_hip(AA_d, s_d, lambda_d, d_d, mu_d, y_d, &
          a_d, zeta, z, m) bind(c, name = 'mma_prepare_aa_matrix_hip')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: AA_d, s_d, lambda_d, d_d, mu_d, y_d, a_d
       real(c_rp), value :: zeta, z
       integer(c_int), value :: m
     end subroutine mma_prepare_aa_matrix_hip

     subroutine delta_1dbeam_hip(Delta_d, L_total, Le, offset, n) &
          bind(c, name = 'delta_1dbeam_hip')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: Delta_d
       real(c_rp) :: L_total, Le
       integer(c_int) :: offset, n
     end subroutine delta_1dbeam_hip

     subroutine hip_Hess(Hess_d, hijx_d, Ljjxinv_d, n, m) bind(c, name = 'hip_Hess')
       import c_int, c_ptr
       type(c_ptr), value :: Hess_d, hijx_d, Ljjxinv_d
       integer(c_int) :: n, m
     end subroutine hip_Hess

     subroutine mma_Ljjxinv_hip(Ljjxinv_d,pjlambda_d, qjlambda_d, x_d, &
          low_d, upp_d, alpha_d, beta_d, n) bind(c, name = 'mma_Ljjxinv_hip')
       import c_int, c_ptr
       type(c_ptr), value :: Ljjxinv_d, x_d, pjlambda_d, qjlambda_d, low_d, &
            upp_d, alpha_d, beta_d
       integer(c_int) :: n
     end subroutine mma_Ljjxinv_hip

     subroutine mma_dipsolvesub1_hip(x_d, pjlambda_d, qjlambda_d, low_d, &
          upp_d, alpha_d, beta_d, n) bind(c, name = 'mma_dipsolvesub1_hip')
       import c_int, c_ptr
       type(c_ptr), value :: x_d, pjlambda_d, qjlambda_d, low_d, &
            upp_d, alpha_d, beta_d
       integer(c_int) :: n
     end subroutine mma_dipsolvesub1_hip

     subroutine mattrans_v_mul_hip(output_d, pij_d, lambda_d, m, n) &
          bind(c, name = 'mattrans_v_mul_hip')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: output_d, pij_d, lambda_d
       integer(c_int) :: m, n
     end subroutine mattrans_v_mul_hip
     subroutine mma_gensub1_hip(low_d, upp_d, x_d, xmin_d, xmax_d, asyinit, n)&
          bind(c, name = 'mma_gensub1_hip')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: low_d, upp_d, x_d, xmin_d, xmax_d
       real(c_rp) :: asyinit
       integer(c_int) :: n
     end subroutine mma_gensub1_hip

     subroutine mma_gensub2_hip(low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d, &
          asydecr, asyincr, n) bind(c, name = 'mma_gensub2_hip')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: low_d, upp_d, x_d, xold1_d, xold2_d, xdiff_d
       real(c_rp) :: asydecr, asyincr
       integer(c_int) :: n
     end subroutine mma_gensub2_hip

     subroutine mma_gensub3_hip(x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, &
          max_d, alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, n, m) &
          bind(c, name = 'mma_gensub3_hip')
       import c_int, c_ptr
       type(c_ptr), value :: x_d, df0dx_d, dfdx_d, low_d, upp_d, min_d, max_d, &
            alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine mma_gensub3_hip

     subroutine mma_gensub4_hip(x_d, low_d, upp_d, pij_d, qij_d, n, m, bi_d) &
          bind(c, name = 'mma_gensub4_hip')
       import c_int, c_ptr
       type(c_ptr), value :: x_d, low_d, upp_d, pij_d, qij_d, bi_d
       integer(c_int) :: n, m
     end subroutine mma_gensub4_hip

     subroutine hip_mma_max(xsi_d, x_d, alpha_d, n) &
          bind(c, name = 'hip_mma_max')
       import c_int, c_ptr
       type(c_ptr), value :: xsi_d, x_d, alpha_d
       integer(c_int) :: n
     end subroutine hip_mma_max

     subroutine hip_rex(rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, q0j_d, &
          lambda_d, xsi_d, eta_d, n, m) bind(c, name = 'hip_rex')
       import c_int, c_ptr
       type(c_ptr), value :: rex_d, x_d, low_d, upp_d, pij_d, p0j_d, qij_d, &
            q0j_d, lambda_d, xsi_d, eta_d
       integer(c_int) :: n, m
     end subroutine hip_rex

     subroutine hip_relambda(relambda_d, x_d, upp_d, low_d, pij_d, qij_d, n, &
          m) bind(c, name = 'hip_relambda')
       import c_int, c_ptr
       type(c_ptr), value :: relambda_d, x_d, upp_d, low_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine hip_relambda

     subroutine hip_sub2cons2(rexsi_d, xsi_d, x_d, alpha_d, epsi, n) &
          bind(c, name = 'hip_sub2cons2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: rexsi_d, xsi_d, x_d, alpha_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine hip_sub2cons2

     real(c_rp) function hip_maxval(rex_d, n) bind(c, name = 'hip_maxval')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: rex_d
       integer(c_int) :: n
     end function hip_maxval

     real(c_rp) function hip_norm(rex_d, n) bind(c, name = 'hip_norm')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: rex_d
       integer(c_int) :: n
     end function hip_norm

     subroutine hip_delx(delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, &
          q0j_d, alpha_d, beta_d, lambda_d, epsi, n, m) &
          bind(c, name = 'hip_delx')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: delx_d, x_d, low_d, upp_d, pij_d, qij_d, p0j_d, &
            q0j_d, alpha_d, beta_d, lambda_d
       real(c_rp) :: epsi
       integer(c_int) :: n, m
     end subroutine hip_delx



     subroutine hip_GG(GG_d, x_d, low_d, upp_d, pij_d, qij_d, n, m) &
          bind(c, name = 'hip_GG')
       import c_int, c_ptr
       type(c_ptr), value :: GG_d, x_d, low_d, upp_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine hip_GG

     subroutine hip_diagx(diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, &
          pij_d, qij_d, alpha_d, beta_d, eta_d, lambda_d, n, m) &
          bind(c, name = 'hip_diagx')
       import c_int, c_ptr
       type(c_ptr), value :: diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d, &
            pij_d, qij_d, alpha_d, beta_d, eta_d, lambda_d
       integer(c_int) :: n, m
     end subroutine hip_diagx

     subroutine hip_bb(bb_d, GG_d, delx_d, diagx_d, n, m) &
          bind(c, name = 'hip_bb')
       import c_int, c_ptr
       type(c_ptr), value :: bb_d, GG_d, delx_d, diagx_d
       integer(c_int) :: n, m
     end subroutine hip_bb

     subroutine hip_AA(AA_d, GG_d, diagx_d, n, m) bind(c, name = 'hip_AA')
       import c_int, c_ptr
       type(c_ptr), value :: AA_d, GG_d, diagx_d
       integer(c_int) :: n, m
     end subroutine hip_AA

     subroutine hip_dx(dx_d, delx_d, diagx_d, GG_d, dlambda_d, n, m) &
          bind(c, name = 'hip_dx')
       import c_int, c_ptr
       type(c_ptr), value :: dx_d, delx_d, diagx_d, GG_d, dlambda_d
       integer(c_int) :: n, m
     end subroutine hip_dx

     subroutine hip_dxsi(dxsi_d, xsi_d, dx_d, x_d, alpha_d, epsi, n) &
          bind(c, name = 'hip_dxsi')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: dxsi_d, xsi_d, dx_d, x_d, alpha_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine hip_dxsi

     subroutine hip_deta(deta_d, eta_d, dx_d, x_d, beta_d, epsi, n) &
          bind(c, name = 'hip_deta')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: deta_d, eta_d, dx_d, x_d, beta_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine hip_deta

     real(c_rp) function hip_maxval2(dxx_d, xx_d, cons, n) &
          bind(c, name = 'hip_maxval2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: dxx_d, xx_d
       real(c_rp) :: cons
       integer(c_int) :: n
     end function hip_maxval2

     real(c_rp) function hip_maxval3(dx_d, x_d, alpha_d, cons, n) &
          bind(c, name = 'hip_maxval3')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: dx_d, x_d, alpha_d
       real(c_rp) :: cons
       integer(c_int) :: n
     end function hip_maxval3

     subroutine hip_kkt_rex(rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d, &
          n, m) bind(c, name = 'hip_kkt_rex')
       import c_int, c_ptr
       type(c_ptr), value :: rex_d, df0dx_d, dfdx_d, xsi_d, eta_d, lambda_d
       integer(c_int) :: n, m
     end subroutine hip_kkt_rex


     subroutine hip_maxcons(a_d, b, c, d_d, n) bind(c, name = 'hip_maxcons')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, d_d
       real(c_rp) :: b, c
       integer(c_int) :: n
     end subroutine hip_maxcons


     real(c_rp) function hip_lcsc2(a_d, b_d, n) bind(c, name = 'hip_lcsc2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
     end function hip_lcsc2

     subroutine hip_mpisum(a_d, n) bind(c, name = 'hip_mpisum')
       import c_int, c_ptr
       type(c_ptr), value :: a_d
       integer(c_int) :: n
     end subroutine hip_mpisum

     subroutine hip_add2inv2(a_d, b_d, c, n) bind(c, name = 'hip_add2inv2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, b_d
       integer(c_int) :: n
       real(c_rp) :: c
     end subroutine hip_add2inv2

     subroutine hip_max2(a_d, b, c_d, d, n) bind(c, name = 'hip_max2')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: a_d, c_d
       integer(c_int) :: n
       real(c_rp) :: b, d
     end subroutine hip_max2

     subroutine hip_updatebb(bb_d, dellambda_d, dely_d, d_d, mu_d, y_d, delz, &
          m) bind(c, name = 'hip_updatebb')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: bb_d, dellambda_d, dely_d, d_d, mu_d, y_d
       integer(c_int) :: m
       real(c_rp) :: delz
     end subroutine hip_updatebb

     subroutine hip_updateAA(AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, &
          y_d, a_d, zeta, z, m) bind(c, name = 'hip_updateAA')
       import c_rp, c_int, c_ptr
       type(c_ptr), value :: AA_d, globaltmp_mm_d, s_d, lambda_d, d_d, mu_d, &
            y_d, a_d
       integer(c_int) :: m
       real(c_rp) :: zeta, z
     end subroutine hip_updateAA

     subroutine hip_dy(dy_d, dely_d, dlambda_d, d_d, mu_d, y_d, n) &
          bind(c, name = 'hip_dy')
       import c_int, c_ptr
       type(c_ptr), value :: dy_d, dely_d, dlambda_d, d_d, mu_d, y_d
       integer(c_int) :: n
     end subroutine hip_dy

  end interface

end module hip_mma_math
