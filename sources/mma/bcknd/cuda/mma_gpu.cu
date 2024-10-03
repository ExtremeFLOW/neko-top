#include <device/device_config.h>
#include <device/cuda/check.h>
#include <mma_gpu_kernel.h>
#include <stdio.h>
#include <stdlib.h>
void mma_gensub_gpu(void* x, void* xold1, void* xold2, void* df0dx, void* dfdx, void* xlow, void* xupp, void* xmin, void* xmax, void* alpha, void* beta, void* p0j, void* q0j, void* pij, void* qij, void* bi,
	real* asyinit, real* asydecr, real* asyincr, int* n, int* m, int* iter) {
	int num1 = *n;
	int num2 = *m;
	const dim3 nthrds(1024, 1, 1);
	const dim3 nblcks((num1 + 1024 - 1) / 1024, 1, 1);
	real* temp;
	real* temp_sum;
	cudaMalloc(&temp, num1 * num2 * sizeof(real));
	cudaMalloc(&temp_sum, num1 * sizeof(real));

	for (int i = 0; i < num2; i++) {
		int nb = (num1 + 2048 - 1) / 2048;
		mmasum_kernel <real><< <nb, 1024 >> > (temp, temp_sum, num1, i);

		reduce_kernel<real> << <1, 1024, 0 >> > (temp_sum, nb);

		cudaMemcpy((real*)bi + num2, temp_sum + num2, sizeof(real), cudaMemcpyDeviceToDevice);
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