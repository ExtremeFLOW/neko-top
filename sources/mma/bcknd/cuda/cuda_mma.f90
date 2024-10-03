module cuda_mma
  use num_types, only: rp, c_rp
  implicit none
  public

  interface
     subroutine mma_gensub_gpu(void* x, void*  xold1, void* xold2, void* df0dx, void* dfdx, void* xlow, void* xupp, &
	 void *xmin, void* xmax, void* alpha, void* beta, void* p0j, void* q0j, void* pij, void* qij, void *bi,
	 real *asyinit, real* asydecr, real* asyincr, int* n, int* m, int* iter) &
          bind(c, name = 'mma_gensub_gpu')

     end subroutine mma_gensub_gpu
  end interface
end module cuda_mma