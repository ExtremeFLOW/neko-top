module cuda_mma
  use num_types, only: rp, c_rp
  implicit none
  public

  interface
     subroutine mma_gensub_gpu( x_d,  xold1_d, xold2_d, df0dx_d, dfdx_d, xlow_d, xupp_d, xmin_d, xmax_d, alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, bi_d,&
	 asyinit, asydecr, asyincr, n, m, iter) &
          bind(c, name = 'mma_gensub_gpu')
		use, intrinsic :: iso_c_binding, only: c_int, c_ptr
		type(c_ptr), value :: x_d,  xold1_d, xold2_d, df0dx_d, dfdx_d, xlow_d, xupp_d, xmin_d, xmax_d, alpha_d, beta_d, p0j_d, q0j_d, pij_d, qij_d, bi_d
		real(c_rp) :: asyinit, asydecr, asyincr
		integer(c_int) :: n, m, iter
     end subroutine mma_gensub_gpu
  end interface
end module cuda_mma