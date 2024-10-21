module cuda_mma_math
  use num_types, only: rp, c_rp
  implicit none
  public
  
  interface
     subroutine mma_gensub1_gpu(low_d, upp_d,x_d, xmin_d, xmax_d, asyinit, n) &
          bind(c, name = 'mma_gensub1_gpu')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: low_d, upp_d,x_d, xmin_d, xmax_d
       real(c_rp) :: asyinit
       integer(c_int) :: n
     end subroutine mma_gensub1_gpu

     subroutine mma_gensub2_gpu(low_d, upp_d, x_d, xold1_d, xold2_d,xmin_d, xmax_d, asydecr, asyincr, n) &
          bind(c, name = 'mma_gensub2_gpu')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: low_d, upp_d, x_d, xold1_d, xold2_d,xmin_d, xmax_d
       real(c_rp) :: asydecr, asyincr
       integer(c_int) :: n
     end subroutine mma_gensub2_gpu

     subroutine mma_gensub3_gpu(x_d, df0dx_d, dfdx_d,low_d, upp_d, min_d,max_d,alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, n, m) &
          bind(c, name = 'mma_gensub3_gpu')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: x_d, df0dx_d, dfdx_d,low_d, upp_d, min_d, max_d,alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine mma_gensub3_gpu

     subroutine mma_gensub4_gpu(x_d, low_d, upp_d, pij_d, qij_d, n, m, bi_d) &
          bind(c, name = 'mma_gensub4_gpu')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: x_d, low_d, upp_d, pij_d, qij_d, bi_d
       integer(c_int) :: n, m
     end subroutine mma_gensub4_gpu

     subroutine cuda_mma_max(xsi_d,x_d,alpha_d,n) &
          bind(c, name = 'cuda_mma_max')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: xsi_d,x_d,alpha_d
       integer(c_int) :: n
     end subroutine cuda_mma_max

     subroutine cuda_rex(rex_d,  x_d,  low_d, upp_d,  pij_d, p0j_d,qij_d, q0j_d, lambda_d, xsi_d, eta_d, n, m)  &
          bind(c, name = 'cuda_rex')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: rex_d,  x_d,  low_d, upp_d,  pij_d, p0j_d,qij_d, q0j_d, lambda_d, xsi_d, eta_d
       integer(c_int) :: n,m
     end subroutine cuda_rex

     subroutine cuda_relambda(relambda_d, x_d,  upp_d, low_d, pij_d, qij_d,  n, m) &
          bind(c, name = 'cuda_relambda')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: relambda_d, x_d,  upp_d, low_d, pij_d, qij_d
       integer(c_int) :: n,m
     end subroutine cuda_relambda

     subroutine cuda_sub2cons2(rexsi_d,xsi_d,x_d,alpha_d,epsi,n) &
          bind(c, name = 'cuda_sub2cons2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: rexsi_d,xsi_d,x_d,alpha_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine cuda_sub2cons2

     subroutine cuda_maxval(rex_d,n) &
          bind(c, name = 'cuda_maxval')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: rex_d
       integer(c_int) :: n
     end subroutine cuda_maxval

     subroutine cuda_norm(rex_d,n) &
          bind(c, name = 'cuda_norm')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: rex_d
       integer(c_int) :: n
     end subroutine cuda_norm

     subroutine cuda_delx(delx_d, x_d, low_d, upp_d,  pij_d,  qij_d,  p0j_d, q0j_d, alpha_d,  beta_d, lambda_d, epsi, n, m) &
          bind(c, name = 'cuda_delx')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: delx_d, x_d, low_d, upp_d,  pij_d,  qij_d,  p0j_d, q0j_d, alpha_d,  beta_d, lambda_d
       real(c_rp) :: epsi
       integer(c_int) :: n, m
     end subroutine cuda_delx

     subroutine cuda_dellambda( dellambda_d, x_d, low_d, upp_d, pij_d, qij_d, n) &
          bind(c, name = 'cuda_dellambda')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: dellambda_d, x_d, low_d, upp_d, pij_d, qij_d
       integer(c_int) :: n
     end subroutine cuda_dellambda

     subroutine cuda_GG(GG_d,  x_d,  low_d,  upp_d, pij_d, qij_d, n, m)&
          bind(c, name = 'cuda_GG')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: GG_d,  x_d,  low_d,  upp_d, pij_d, qij_d
       integer(c_int) :: n, m
     end subroutine cuda_GG

     subroutine cuda_diagx(diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d,  pij_d, qij_d, &
       alpha_d, beta_d,  eta_d, lambda_d, n, m) bind(c, name = 'cuda_diagx')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: diagx_d, x_d, xsi_d, low_d, upp_d, p0j_d, q0j_d,  pij_d, qij_d,  alpha_d, beta_d,  eta_d, lambda_d
       integer(c_int) :: n, m
     end subroutine cuda_diagx

     subroutine cuda_bb(bb_d, GG_d, delx_d,diagx_d,n,m) &
          bind(c, name = 'cuda_bb')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: bb_d, GG_d, delx_d,diagx_d
       integer(c_int) :: n, m
     end subroutine cuda_bb

     subroutine cuda_AA(AA_d, GG_d,  diagx_d, n, m)  &
          bind(c, name = 'cuda_AA')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: AA_d, GG_d,  diagx_d
       integer(c_int) :: n, m
     end subroutine cuda_AA

     subroutine cuda_dx(dx_d, delx_d, diagx_d, GG_d, dlambda_d, n, m) &
          bind(c, name = 'cuda_dx')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: dx_d, delx_d, diagx_d, GG_d, dlambda_d
       integer(c_int) :: n,m
     end subroutine cuda_dx

     subroutine cuda_dxsi(dxsi_d, xsi_d, dx_d,x_d,alpha_d, epsi, n)  &
          bind(c, name = 'cuda_dxsi')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp
       type(c_ptr), value :: dxsi_d, xsi_d, dx_d,x_d,alpha_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine cuda_dxsi

     subroutine cuda_deta(deta_d, eta_d, dx_d,  x_d, beta_d, epsi,n) &
          bind(c, name = 'cuda_deta')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp       
       type(c_ptr), value :: deta_d, eta_d, dx_d,  x_d, beta_d
       real(c_rp) :: epsi
       integer(c_int) :: n
     end subroutine cuda_deta

     subroutine cuda_maxval2(dxx_d, xx_d, cons, n) &
          bind(c, name = 'cuda_maxval2')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp       
       type(c_ptr), value :: dxx_d, xx_d
       real(c_rp) :: cons
       integer(c_int) :: n
     end subroutine cuda_maxval2

     subroutine cuda_maxval3(dx_d, x_d, alpha_d,cons, n) &
          bind(c, name = 'cuda_maxval3')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       import c_rp       
       type(c_ptr), value :: dx_d, x_d, alpha_d
       real(c_rp) :: cons
       integer(c_int) :: n
     end subroutine cuda_maxval3

     subroutine cuda_kkt_rex(rex_d,  df0dx_d,  dfdx_d, xsi_d, eta_d, lambda_d, n, m) &
          bind(c, name = 'cuda_kkt_rex')
       use, intrinsic :: iso_c_binding, only: c_int, c_ptr
       type(c_ptr), value :: rex_d,  df0dx_d,  dfdx_d, xsi_d, eta_d, lambda_d
       integer(c_int) :: n,m
     end subroutine cuda_kkt_rex


     subroutine cuda_maxcons(a_d,b,c, d_d, n) &
          bind(c, name = 'cuda_maxcons')
          use, intrinsic :: iso_c_binding, only: c_int, c_ptr
          import c_rp       
          type(c_ptr), value :: a_d,d_d
          real(c_rp) :: b,c
          integer(c_int) :: n
     end subroutine cuda_maxcons


     subroutine cuda_lcsum(a_d, n) &
          bind(c, name = 'cuda_lcsum')
          use, intrinsic :: iso_c_binding, only: c_int, c_ptr
          type(c_ptr), value :: a_d
          integer(c_int) :: n
     end subroutine cuda_lcsum

     subroutine  cuda_lcsc2(a_d, b_d, n)  &
          bind(c, name = 'cuda_lcsc2')
          use, intrinsic :: iso_c_binding, only: c_int, c_ptr
          type(c_ptr), value :: a_d, b_d
          integer(c_int) :: n
     end subroutine cuda_lcsc2

     subroutine cuda_mpisum(a_d, n)  &
          bind(c, name = 'cuda_mpisum')
          use, intrinsic :: iso_c_binding, only: c_int, c_ptr
          type(c_ptr), value :: a_d
          integer(c_int) :: n
     end subroutine cuda_mpisum

     subroutine cuda_add2inv2(a_d, b_d, c, n)  &
          bind(c, name = 'cuda_add2inv2')
          use, intrinsic :: iso_c_binding, only: c_int, c_ptr
          import c_rp       
          type(c_ptr), value :: a_d, b_d
          integer(c_int) :: n
          real(c_rp) :: c
     end subroutine cuda_add2inv2


    end interface
end module cuda_mma_math