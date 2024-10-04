#include <device/device_config.h>
#include <device/cuda/check.h>
#include <mma_gpu_kernel.h>
#include <stdio.h>
#include <stdlib.h>
void mma_gensub_gpu(void* x, void* xold1, void* xold2, void* df0dx, void* dfdx, void* xlow, void* xupp, void* xmin, void* xmax,
  void* alpha, void* beta, void* p0j, void* q0j, void* pij, void* qij, void* bi,
  real* asyinit, real* asydecr, real* asyincr, int* n, int* m, int* iter) {
  int num1 = *n;
  int num2 = *m;

  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks((num1 + 1024 - 1) / 1024, 1, 1);
  real* temp;
  real* temp_sum;
  cudaMalloc(&temp, num1 * num2 * sizeof(real));
  cudaMalloc(&temp_sum, num1 * sizeof(real));


  BoundCalculation_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>> ((real*)x, (real*)xold1, (real*)xold2, 
   (real*)df0dx, (real*)dfdx, (real*)xlow, (real*)xupp, (real*)xmin, (real*)xmax,
   (real*)alpha, (real*)beta, (real*)p0j, (real*)q0j, (real*)pij, (real*)qij, temp,
   *asyinit, *asydecr, *asyincr, *n, *m, *iter);




  for (int i = 0; i < num2; i++) {
    int nb = (num1 + 2048 - 1) / 2048;
    mmasum_kernel <real> <<<nb, 1024 >>> (temp, temp_sum, num1, i);

    mmareduce_kernel<real> <<<1, 1024, 0 >>> (temp_sum, nb);

    cudaMemcpy((real*)bi + i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);
 }
 cudaFree(temp);
 cudaFree(temp_sum);
    ////global communication ?
    /*
    globaltmp_m = 0.0_rp
    call MPI_Allreduce(this%bi%x, globaltmp_m, this%m, &
         mpi_real_precision, mpi_sum, neko_comm, ierr)
    this%bi%x = globaltmp_m - fval
    */

 }

 void cuda_max(void * x,void * alpha,void * beta, void * xsi,void * eta,void * mu,void * c, int * n) {

  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

  max_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>((real *) x, (real *) alpha, (real *) beta,
     (real *) xsi, (real *) eta, (real *) mu,
     (real *) c, *n);
  CUDA_CHECK(cudaGetLastError());

}

void cuda_rex(void* rex, void* x, void* xlow, void* xupp, void* pij, void* p0j,
	void* qij, void* q0j, void* lambda, void* xsi, void* eta, int* n, int* m)
{
	const dim3 nthrds(1024, 1, 1);
	const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

	RexCalculation_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)rex, (real*)x, (real*)xlow, (real*)xupp, (real*)pij, (real*)p0j,
		(real*)qij, (real*)q0j, (real*)lambda, (real*)xsi, (real*)eta, *n, *m);
	CUDA_CHECK(cudaGetLastError());

}

