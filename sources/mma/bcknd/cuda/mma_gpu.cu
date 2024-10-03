#include <device/device_config.h>
#include <device/cuda/check.h>
#include <mma_gpu_kernel.h>
#include <stdio.h>
#include <stdlib.h>
void mma_gensub_gpu(void* x, void*  xold1, void* xold2, void* df0dx, void* dfdx, void* xlow, void* xupp, void *xmin, void* xmax, void* alpha, void* beta, void* p0j, void* q0j, void* pij, void* qij, void *bi,
    real *asyinit, real* asydecr, real* asyincr, int* n, int* m, int* iter) {
    int num1 = *n;
    int num2 = *m;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((num1 + 1024 - 1) / 1024, 1, 1);
    real* temp;
    real* temp_sum;
    cudaMalloc(&temp, num1 * num2 * sizeof(real));
    cudaMalloc(&temp_sum, num1* sizeof(real));

    BoundCalculation_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real *)x, (real*)xold1, (real*)xold2, (real*)df0dx, (real*)dfdx, (real*)xlow, (real*)xupp, (real*)xmin,
        (real*)xmax, (real*)alpha, (real*)beta, (real*)p0j, (real*)q0j, (real*)pij, (real*)qij, temp, *asyinit, *asydecr, *asyincr, *n, *m, *iter);
 
    for (int i = 0; i < num2; i++) {
        sum_reductions << <num1 / 2048 + 1, nthrds >> > (temp, temp_sum, num1, i);

        int k = 1;
        while (num1 / pow(2048, k) > 1) {
            sum_reductions_final << <num1 / pow(2048, k) + 1, nthrds >> > (temp_sum, num1 / pow(2048, k - 1) + 1);
            k = k + 1;
        }
        cudaMemcpy((real*)bi + num2, temp_sum + num2, sizeof(real), cudaMemcpyDeviceToDevice);
    }
    cudaFree(temp);
    cudaFree(temp_sum);
}