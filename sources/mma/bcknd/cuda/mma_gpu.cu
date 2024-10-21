#include "../../../../external/neko/src/device/cuda/check.h"
#include "mma_gpu_kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include "../../../../external/neko/src/math/bcknd/device/device_mpi_reduce.h"
#include "../../../../external/neko/src/math/bcknd/device/device_mpi_op.h"
#include "../../../../external/neko/src/device/device_config.h"
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

void cuda_mma_max(void* xsi, void* x, void* alpha, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    mma_max2_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>((real *) xsi, (real *) x, (real *) alpha, *n);
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

void cuda_sub2cons2(void* a, void* b, void* c, void* d, real* e, int* n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    sub2cons2_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)a, (real*)b, (real*)c, (real*)d, *e, *n);
    //CUDA_CHECK(cudaGetLastError());
}


//////////////max abs values of input
real cuda_maxval(void* a, int* n) {
    real* temp;
    real* temp_cpu = new real[1];
    cudaMalloc(&temp, (*n) * sizeof(real));
    int nb = ((*n) + 2048 - 1) / 2048;
    maxval_kernel <real> << <nb, 1024 >> > ((real*)a, temp, (*n));
    max_reduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
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
    cudaMemcpy(dellambda_d, temp, sizeof(real), cudaMemcpyDeviceToDevice);
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
        if(i<*m)
            cudaMemcpy(AA_d + i, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);
        else
            cudaMemcpy(AA_d + i+1, temp_sum, sizeof(real), cudaMemcpyDeviceToDevice);

    }
    cudaFree(temp);
    cudaFree(temp_sum);
    //CUDA_CHECK(cudaGetLastError());
}


void cuda_dx(void* dx,void* delx, void* diagx, void* GG, void* dlambda, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dx_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> >((real*)dx, (real*)delx,  (real*)diagx,(real*) GG,
        (real*)dlambda, *n,*m);
    //CUDA_CHECK(cudaGetLastError());
}

void cuda_dxsi(void* dxsi, void* xsi, void* dx,void* x,
    void* alpha, real*epsi, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dxsi_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> >((real*)dxsi, (real*)xsi,  (real*)dx,(real*) x,
        (real*)alpha, *epsi,*n);
    //CUDA_CHECK(cudaGetLastError());
}

void cuda_deta(void* deta, void* eta, void* dx, void* x,
    void* beta, real* epsi, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    deta_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)deta, (real*)eta, (real*)dx, (real*)x,
        (real*)beta, *epsi, *n);
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



void cuda_rex(void* rex, void* x, void* xlow, void* xupp, void* pij, void* p0j,
	void* qij, void* q0j, void* lambda, void* xsi, void* eta, int* n, int* m)
{
	const dim3 nthrds(1024, 1, 1);
	const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

	RexCalculation_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)rex, (real*)x, (real*)xlow, (real*)xupp, (real*)pij, (real*)p0j,
		(real*)qij, (real*)q0j, (real*)lambda, (real*)xsi, (real*)eta, *n, *m);
	//CUDA_CHECK(cudaGetLastError());

}



void cuda_rey(void* rey, void* c, void* d, void* y, void* lambda, void* mu, int* n)
{
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    rey_calculation_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)rey, (real*)c, (real*)d, (real*)y,
        (real*)lambda, (real*)mu, * n);
    //CUDA_CHECK(cudaGetLastError());
}





/////a_d=b_d*c_d-d
void cuda_sub2cons(void * a,void * b,void * c, real *d, int * n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);

    sub2cons_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>((real *) a, (real *) b, (real *) c,
      *d, *n);
    //CUDA_CHECK(cudaGetLastError());
}




/////sum(a^2)
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






void cuda_dely(void* dely, void* c, void* d, void* y, void* lambda, real* epsi, int* n)
{
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dely_kernel<real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)dely,(real*)c, (real*)d, (real*)y, (real*)lambda,*epsi, * n);
    //CUDA_CHECK(cudaGetLastError());
}















real cuda_maxval2(void* a, void* b, real* cons, int* n) {
    real* temp;
    real* temp_cpu = new real[1];
    cudaMalloc(&temp, (*n) * sizeof(real));
    int nb = ((*n) + 2048 - 1) / 2048;
    maxval2_kernel <real> << <nb, 1024 >> > ((real*)a, (real*)b, temp, (*cons), (*n));
    max_reduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
    //CUDA_CHECK(cudaGetLastError());
    cudaMemcpy(temp_cpu, temp, sizeof(real), cudaMemcpyDeviceToHost);
    cudaFree(temp);
    return temp_cpu[0];
}

real cuda_maxval3(void* a, void* b, void* c, real* cons, int* n) {
    real* temp;
    real* temp_cpu = new real[1];
    cudaMalloc(&temp, (*n) * sizeof(real));
    int nb = ((*n) + 2048 - 1) / 2048;
    maxval3_kernel <real> << <nb, 1024 >> > ((real*)a, (real*)b, (real*)c, temp, (*cons), (*n));
    max_reduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
    //CUDA_CHECK(cudaGetLastError());
    cudaMemcpy(temp_cpu, temp, sizeof(real), cudaMemcpyDeviceToHost);
    cudaFree(temp);
    return temp_cpu[0];
}


void cuda_kkt_rex(void* rex, void* df0dx, void* dfdx, void* xsi,
    void* eta, void* lambda, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    kkt_rex_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)rex, (real*)df0dx, (real*)dfdx, (real*)xsi,
        (real*)eta, (real*)lambda, *n, *m);
    //CUDA_CHECK(cudaGetLastError());
}


//////a_d=max(b,c*d_d)
void cuda_maxcons(void* a, real* b, real* c, void* d, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    maxcons_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)a, *b, *c, (real*)d, *n);
    //CUDA_CHECK(cudaGetLastError());
}


real cuda_lcsum(void *a, int *n) {
   real* temp;
   real* temp_cpu = new real[1];
   cudaMalloc(&temp, (*n) * sizeof(real));
   int nb = ((*n) + 2048 - 1) / 2048;
   glsum_kernel <real> << <nb, 1024 >> > ((real*)a,  temp, (*n));
   mmareduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
    //CUDA_CHECK(cudaGetLastError());
   cudaMemcpy(temp_cpu, temp, sizeof(real), cudaMemcpyDeviceToHost);
   cudaFree(temp);
   return temp_cpu[0];
}

real cuda_lcsc2(void *a, void*b, int *n) {
   real* temp;
   real* temp_cpu = new real[1];
   cudaMalloc(&temp, (*n) * sizeof(real));
   int nb = ((*n) + 2048 - 1) / 2048;
   glsc2_kernel <real> << <nb, 1024 >> > ((real*)a, (real*)b, temp, (*n));
   mmareduce_kernel<real> << <1, 1024, 0 >> > (temp, nb);
    //CUDA_CHECK(cudaGetLastError());
   cudaMemcpy(temp_cpu, temp, sizeof(real), cudaMemcpyDeviceToHost);
   cudaFree(temp);
   return temp_cpu[0];
}


void cuda_mpisum(void *a, int *n) {
    real* temp=(real*)a;
 #ifdef HAVE_DEVICE_MPI
    cudaStreamSynchronize(stream);
    device_mpi_allreduce_inplace(temp, *n, sizeof(real), DEVICE_MPI_SUM);
 #endif
}
  
 void cuda_add2inv2(void* a, void *b, real* c, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    add2inv2_kernel <real> << <nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >> > ((real*)a, (real*) b, *c, *n);
    //CUDA_CHECK(cudaGetLastError());
}



