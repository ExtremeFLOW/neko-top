#include <device/device_config.h>
#include <device/cuda/check.h>
#include <mma_gpu_kernel.h>
#include <stdio.h>
#include <stdlib.h>
#include <math/bcknd/device/device_mpi_reduce.h>
#include <math/bcknd/device/device_mpi_op.h>

void mma_gensub4_gpu(void* x, void* low, void* upp, void* pij, void* qij, int* n, int* m, void* bi) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    real* temp;
    real* temp_sum;
    real* bi_d = (real*)bi;
    cudaMalloc(&temp, (*m) * (*n) * sizeof(real));
    cudaMalloc(&temp_sum, (*n) * sizeof(real));
    mma_sub4_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)x, (real*)low,
        (real*)upp, (real*)pij, (real*)qij, temp, *n, *m);
    for (int i = 0; i < (*m); i++) {
        int nb = ((*n) + 2048 - 1) / 2048;
        mmasum_kernel <real> << <nb, 1024 >> > (temp, temp_sum, (*n), i);
        mmareduce_kernel<real> << <1, 1024, 0 >> > (temp_sum, nb);
        cudaMemcpy(bi_d + i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);

    }
    cudaFree(temp);
    cudaFree(temp_sum);

    //CUDA_CHECK(cudaGetLastError());
}



void mma_gensub3_gpu(void* x, void* df0dx, void* dfdx, void* low, void* upp, void* xmin, void* xmax, void* alpha, void* beta,
    void* p0j, void* q0j, void* pij, void* qij, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    mma_sub3_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)x, (real*)df0dx, (real*)dfdx, (real*)low,
        (real*)upp, (real*)xmin, (real*)xmax, (real*)alpha, (real*)beta, (real*)p0j, (real*)q0j, (real*)pij, (real*)qij, *n, *m);
    //CUDA_CHECK(cudaGetLastError());
}



void mma_gensub2_gpu(void* low, void* upp, void* x, void* xold1, void* xold2, void* xmin, void* xmax,
    real* asydecr, real* asyincr, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    mma_sub2_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)low, (real*)upp, (real*)x, (real*)xold1,
        (real*)xold2, (real*)xmin, (real*)xmax, *asydecr, *asyincr, *n);
    //CUDA_CHECK(cudaGetLastError());
}


void mma_gensub1_gpu(void* low, void* upp, void* x, void* xmin, void* xmax, real* asyinit, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    mma_sub1_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)low, (real*)upp, 
        (real*)x, (real*)xmin, (real*)xmax, *asyinit, *n);
    //CUDA_CHECK(cudaGetLastError());
}





void mma_gensub_gpu(void* x, void* xold1, void* xold2, void* df0dx, void* dfdx, void* xlow, void* xupp, void* xmin, void* xmax,
    void* alpha, void* beta, void* p0j, void* q0j, void* pij, void* qij, void* bi,
    real* asyinit, real* asydecr, real* asyincr, int* n, int* m, int* iter) {
    int num1 = *n;
    int num2 = *m;
    real* bi_d = (real*)bi;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((num1 + 1024 - 1) / 1024, 1, 1);
    real* temp;
    real* temp_sum;
    cudaMalloc(&temp, num1 * num2 * sizeof(real));
    cudaMalloc(&temp_sum, num1 * sizeof(real));


    BoundCalculation_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)x, (real*)xold1, (real*)xold2,
        (real*)df0dx, (real*)dfdx, (real*)xlow, (real*)xupp, (real*)xmin, (real*)xmax,
        (real*)alpha, (real*)beta, (real*)p0j, (real*)q0j, (real*)pij, (real*)qij, temp,
        *asyinit, *asydecr, *asyincr, *n, *m, *iter);


    for (int i = 0; i < num2; i++) {
        int nb = (num1 + 2048 - 1) / 2048;
        mmasum_kernel <real> << <nb, 1024 >> > (temp, temp_sum, num1, i);

        mmareduce_kernel<real> << <1, 1024, 0 >> > (temp_sum, nb);
        cudaMemcpy(bi_d+i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);

    }
    cudaFree(temp);
    cudaFree(temp_sum);
}

void cuda_max(void * x,void * alpha,void * beta, void * xsi,void * eta,void * mu,void * c, int * n) {

	const dim3 nthrds(1024, 1, 1);
	const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

	max_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>((real *) x, (real *) alpha, (real *) beta,
		(real *) xsi, (real *) eta, (real *) mu,
		(real *) c, *n);
	//CUDA_CHECK(cudaGetLastError());
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



void cuda_rey(void* rey, void* c, void* d, void* y, void* lambda, void* mu, int* n)
{
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    rey_calculation_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)rey, (real*)c, (real*)d, (real*)y,
        (real*)lambda, (real*)mu, * n);
    //CUDA_CHECK(cudaGetLastError());
}



void cuda_relambda(void* relambda,  void* x,  void* xupp,  void* xlow,
   void* pij,  void* qij,  int* n,  int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    real* temp;
    real* temp_sum;
    real* relambda_d = (real*)relambda;
    cudaMalloc(&temp, (*n) * (*m) * sizeof(real));
    cudaMalloc(&temp_sum, (*n) * sizeof(real));
    relambda_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > (temp, (real*)x, (real*)xupp, (real*)xlow,
        (real*)pij, (real*)qij, *n, *m);
    for (int i = 0; i < (*m); i++) {
        int nb = ((*n) + 2048 - 1) / 2048;
        mmasum_kernel <real> << <nb, 1024 >> > (temp, temp_sum, (*n), i);

        mmareduce_kernel<real> << <1, 1024, 0 >> > (temp_sum, nb);

        cudaMemcpy(relambda_d + i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);

    }
    cudaFree(temp);
    cudaFree(temp_sum);
    //CUDA_CHECK(cudaGetLastError());
    
}



/////a=b*c-d
void cuda_sub2cons(void * a,void * b,void * c, real *d, int * n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    sub2cons_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>((real *) a, (real *) b, (real *) c,
      *d, *n);
    //CUDA_CHECK(cudaGetLastError());
}

/////a=b*c-d
void cuda_sub2cons2(void * a,void * b,void * c, void * d,real *e, int * n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    sub2cons_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>((real *) a, (real *) b, (real *) c, (real*) d
      *e, *n);
    //CUDA_CHECK(cudaGetLastError());
}



real cuda_maxval(void* a, int* n) {
    real* temp;
    real* temp_cpu=new real[1];
    cudaMalloc(&temp, (*n) * sizeof(real));
    int nb = ((*n) + 2048 - 1) / 2048;
    maxval_kernel <real> << <nb, 1024 >> > ((real*)a, temp, (*n));
    max_reduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
    //CUDA_CHECK(cudaGetLastError());
    cudaMemcpy(temp_cpu, temp, sizeof(real), cudaMemcpyDeviceToHost);
    cudaFree(temp);
    return temp_cpu[0];
}


real cuda_norm(void* a, int* n) {
 real* temp;
 real* temp_cpu=new real[1];
 cudaMalloc(&temp, (*n) * sizeof(real));
 int nb = ((*n) + 2048 - 1) / 2048;
 norm_kernel <real> << <nb, 1024 >> > ((real*)a, temp, (*n));
 mmareduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
    //CUDA_CHECK(cudaGetLastError());
 cudaMemcpy(temp_cpu, temp, sizeof(real), cudaMemcpyDeviceToHost);
 cudaFree(temp);

 return temp_cpu[0];
}



void cuda_delx(void* delx, void* x, void* xlow, void* xupp, void* pij, void* qij, void * p0j, void *q0j, void *alpha, void * beta, 
    void *lambda, real* epsi, int* n, int *m)
{
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    delx_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)delx, (real*)x, (real*)xlow, (real*)xupp,
        (real*)pij, (real*)qij, (real*)p0j, (real*)q0j,(real*)alpha, (real*)beta, (real*) lambda, *epsi, * n,*m);
    //CUDA_CHECK(cudaGetLastError());
}


void cuda_dely(void* dely, void* c, void* d, void* y, void* lambda, real* epsi, int* n)
{
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dely_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)c, (real*)d, (real*)y, (real*)lambda,*epsi, * n);
    //CUDA_CHECK(cudaGetLastError());
}


void cuda_dellambda(void* dellambda, void*x,  void*xlow,void*  xupp,
   void* pij, void* qij, int * n) {
    real *dellambda_d=(real *) dellambda;
    real* temp;
    real* temp_cpu=new real[1];
    cudaMalloc(&temp, (*n) * sizeof(real));
    int nb = ((*n) + 2048 - 1) / 2048;
    dellambda_kernel <real> << <nb, 1024 >> > (temp, (real*)x,  (real*)xlow,(real*) xupp,
        (real*)pij, (real*) qij, *n);
    mmareduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
    cudaFree(temp);
    //CUDA_CHECK(cudaGetLastError());
    cudaMemcpy(dellambda_d, temp, sizeof(real), cudaMemcpyDeviceToDeivce);
}

void cuda_GG(void* GG, void* x, void* xlow, void* xupp,
    void* pij, void* qij, int* n, int* m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    GG_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> >((real*)GG, (real*)x,  (real*)xlow,(real*) xupp,
        (real*)pij, (real*) qij, *n,*m);
    //CUDA_CHECK(cudaGetLastError());
}

void cuda_diagx(void* diagx, void* x, void* xsi,void* xlow, void* xupp,
    void* p0j, void* q0j, void* pij, void* qij, void* alpha, void* beta, void* eta, void* lambda, int *n, int *m) {

   const dim3 nthrds(1024, 1, 1);
   const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
   diagx_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> >((real*)diagx, (real*)x,  (real*)xsi,(real*)xlow,
    (real*) xupp,(real*)p0j, (real*) q0j, (real*)pij, (real*) qij, (real*)alpha, (real*) beta, (real*)eta, (real*) lambda, *n,*m);
    //CUDA_CHECK(cudaGetLastError());
}




void cuda_bb(void* bb, void* GG, void* delx,void* diagx, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    real* temp;
    real* temp_sum;
    real* bb_d = (real*)bb;
    cudaMalloc(&temp, (*n) * (*m) * sizeof(real));
    cudaMalloc(&temp_sum, (*n) * sizeof(real));
    bb_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> >((real*)temp, (real*) GG, (real*) delx, 
       (real*)diagx, *n, *m);
    for (int i = 0; i < (*m); i++) {
        int nb = ((*n) + 2048 - 1) / 2048;
        mmasum_kernel <real> << <nb, 1024 >> > (temp, temp_sum, (*n), i);

        mmareduce_kernel<real> << <1, 1024, 0 >> > (temp_sum, nb);

        cudaMemcpy(bb_d + i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);

    }
    cudaFree(temp);
    cudaFree(temp_sum);
    //CUDA_CHECK(cudaGetLastError());
}

void cuda_AA(void* AA, void* GG, void* diagx, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    real* temp;
    real* temp_sum;
    real* AA_d = (real*)AA;
    cudaMalloc(&temp, (*n) * (*m) * sizeof(real));
    cudaMalloc(&temp_sum, (*n) * sizeof(real));
    AA_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> >((real*)temp, (real*) GG, (real*) diagx, *n, *m);
    for (int i = 0; i < (*m)*(*m); i++) {
        int nb = ((*n) + 2048 - 1) / 2048;
        mmasum_kernel <real> << <nb, 1024 >> > (temp, temp_sum, (*n), i);

        mmareduce_kernel<real> << <1, 1024, 0 >> > (temp_sum, nb);

        cudaMemcpy(AA_d + i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);

    }
    cudaFree(temp);
    cudaFree(temp_sum);
    //CUDA_CHECK(cudaGetLastError());
}



void device_to_host_copy(void *a, real *b, int *n){
    cudaMemcpy(AA_d + i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);
}
