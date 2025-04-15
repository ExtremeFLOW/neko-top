#include "device/cuda/check.h"
#include "mma_kernel.h"
#include <stdio.h>
#include <stdlib.h>


extern "C" {
#include "math/bcknd/device/device_mpi_reduce.h"
#include "math/bcknd/device/device_mpi_op.h"
#include "device/device_config.h"


  int mma_red_s = 0;
  real * mma_bufred = NULL;
  real * mma_bufred_d = NULL;
 //////
 
  void mma_dipsolvesub1_cuda(void* x, void* pjlambda, void* qjlambda,
     void* low, void* upp, void* alpha, void* beta, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    mma_dipsolvesub1_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
         ((real*)x, (real*)pjlambda, (real*)qjlambda, (real*)low, (real*)upp,
         (real*)alpha, (real*)beta, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  void mattrans_v_mul_cuda(void* output, void* pij, void* lambda,
       int* m, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    mattrans_v_mul_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
         ((real*)output, (real*)pij, (real*)lambda, *m, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  void mma_gensub4_cuda(void* x, void* low, void* upp, void* pij, void* qij, 
       int* n, int* m, void* bi) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1) / 1024;
    const cudaStream_t stream = (cudaStream_t)glb_cmd_queue;

    if (nb > mma_red_s) {
       mma_red_s = nb;
       if (mma_bufred != NULL) {
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred, 
            nb * sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, 
            nb * sizeof(real)));
    }

    real* temp;
    real* bi_d = (real*)bi;
    cudaMalloc(&temp, (*m) * (*n) * sizeof(real));

    mma_sub4_kernel<real><<<nblcks, nthrds, 0, stream>>>(
         (real*)x, (real*)low, (real*)upp, (real*)pij, (real*)qij, 
         temp, *n, *m);

    for (int i = 0; i < (*m); i++) {
       mmasum_kernel<real><<<nblcks, nthrds, 0, stream>>>(
            temp, mma_bufred_d, (*n), (*m), i);
       CUDA_CHECK(cudaGetLastError());

       mmareduce_kernel<real><<<1, 1024, 0, stream>>>(
            mma_bufred_d, nb);
       CUDA_CHECK(cudaGetLastError());

       CUDA_CHECK(cudaMemcpyAsync(
            bi_d + i, mma_bufred_d, sizeof(real), 
            cudaMemcpyDeviceToDevice, stream));
       
       cudaStreamSynchronize(stream);
    }

    cudaFree(temp);
  }
 
 //////
   void mma_gensub3_cuda(void* x, void* df0dx, void* dfdx, void* low, 
       void* upp, void* xmin, void* xmax, void* alpha, void* beta, 
       void* p0j, void* q0j, void* pij, void* qij, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    mma_sub3_kernel<real><<<nblcks, nthrds, 0, 
         (cudaStream_t)glb_cmd_queue>>>(
         (real*)x, (real*)df0dx, (real*)dfdx, (real*)low, 
         (real*)upp, (real*)xmin, (real*)xmax, (real*)alpha, 
         (real*)beta, (real*)p0j, (real*)q0j, (real*)pij, 
         (real*)qij, *n, *m);

    CUDA_CHECK(cudaGetLastError());
  }

  void mma_gensub2_cuda(void* low, void* upp, void* x, void* xold1, 
       void* xold2, void* xmin, void* xmax, real* asydecr, 
       real* asyincr, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    mma_sub2_kernel<real><<<nblcks, nthrds, 0, 
         (cudaStream_t)glb_cmd_queue>>>(
         (real*)low, (real*)upp, (real*)x, (real*)xold1, 
         (real*)xold2, (real*)xmin, (real*)xmax, *asydecr, 
         *asyincr, *n);

    CUDA_CHECK(cudaGetLastError());
  }

 
 
  void mma_gensub1_cuda(void* low, void* upp, void* x, void* xmin, void* xmax,
       real* asyinit, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    mma_sub1_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
         ((real*)low, (real*)upp, (real*)x, (real*)xmin, (real*)xmax,
         *asyinit, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  void cuda_mma_max(void* xsi, void* x, void* alpha, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    mma_max2_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
         ((real*)xsi, (real*)x, (real*)alpha, *n);
    CUDA_CHECK(cudaGetLastError());
  }
 
  void cuda_relambda(void* relambda,  void* x,  void* xupp,  void* xlow,
       void* pij,  void* qij,  int* n,  int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;   

    if ( nb > mma_red_s){
       mma_red_s = nb;
       if (mma_bufred != NULL) {
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));        
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred,nb*sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb*sizeof(real)));
    }
    real* temp;
    cudaMalloc(&temp, (*n) * (*m) * sizeof(real));
    relambda_kernel<real> <<<nblcks, nthrds, 0, stream >>> (temp, (real*)x, 
         (real*)xupp, (real*)xlow, (real*)pij, (real*)qij, *n, *m);
    for (int i = 0; i < (*m); i++) {
       mmasum_kernel <real> <<<nblcks, nthrds, 0, stream >>>
            (temp, mma_bufred_d, (*n),(*m), i);
       CUDA_CHECK(cudaGetLastError());
       mmareduce_kernel<real> <<<1, 1024, 0, stream  >>>
            (mma_bufred_d, nb);
       CUDA_CHECK(cudaGetLastError());
       mma_copy_kernel<<<1, 1, 0, stream>>>((real*)relambda, mma_bufred_d, 1, i);
       CUDA_CHECK(cudaGetLastError());
       cudaStreamSynchronize(stream);
    }
    cudaFree(temp);
  }
 
  void cuda_sub2cons2(void* a, void* b, void* c, void* d, real* e, int* n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  
    sub2cons2_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
         (real*)a, (real*)b, (real*)c, (real*)d, *e, *n);
    CUDA_CHECK(cudaGetLastError());
  }
 
 //////////////max abs values of input
 real cuda_maxval(void* a, int* n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1) / 1024;
    const cudaStream_t stream = (cudaStream_t)glb_cmd_queue;
  
    if (nb > mma_red_s) {
       mma_red_s = nb;
       if (mma_bufred != NULL) {
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred, nb * sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb * sizeof(real)));
    }
  
    maxval_kernel<real><<<nblcks, nthrds, 0, stream>>>(
         (real*)a, mma_bufred_d, (*n));
    CUDA_CHECK(cudaGetLastError());
  
    max_reduce_kernel<real><<<1, 1024, 0, stream>>>(
         mma_bufred_d, nb);
    CUDA_CHECK(cudaGetLastError());
  
    CUDA_CHECK(cudaMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
         cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
  
    return mma_bufred[0];
  }
 

  void cuda_delx(void* delx, void* x, void* xlow, void* xupp, void* pij,
       void* qij, void* p0j, void* q0j, void* alpha, void* beta, void* lambda, 
       real* epsi, int* n, int* m) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    delx_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
         (real*)delx, (real*)x, (real*)xlow, (real*)xupp, (real*)pij, 
         (real*)qij, (real*)p0j, (real*)q0j, (real*)alpha, (real*)beta,
         (real*)lambda, *epsi, *n, *m);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_GG(void* GG, void* x, void* xlow, void* xupp,
       void* pij, void* qij, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    GG_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)GG, (real*)x,  (real*)xlow,(real*) xupp,
         (real*)pij, (real*) qij, *n,*m);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_diagx(void* diagx, void* x, void* xsi,void* xlow, void* xupp,
       void* p0j, void* q0j, void* pij, void* qij, void* alpha, void* beta,
       void* eta, void* lambda, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    diagx_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)diagx, (real*)x,  (real*)xsi,(real*)xlow,
         (real*) xupp,(real*)p0j, (real*) q0j, (real*)pij, (real*) qij,
         (real*)alpha, (real*) beta, (real*)eta, (real*) lambda, *n,*m);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_bb(void* bb, void* GG, void* delx,void* diagx, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;      
    cudaStreamSynchronize(stream);
    if(nb > mma_red_s){
       mma_red_s = nb;
       if(mma_bufred != NULL){
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred,nb*sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb*sizeof(real)));
    }
    for (int i = 0; i < (*m); i++) {
       mmasumbb_kernel <real> <<<nblcks, nthrds, 0, stream>>> 
            ((real*)GG,(real*)delx,(real*)diagx, mma_bufred_d, (*n),(*m), i);
       CUDA_CHECK(cudaGetLastError());
       mmareduce_kernel<real> <<<1, 1024, 0, stream>>> (mma_bufred_d, nb);
       CUDA_CHECK(cudaGetLastError());
       mma_copy_kernel<<<1, 1, 0, stream>>>((real*)bb, mma_bufred_d, 1, i);
       CUDA_CHECK(cudaGetLastError());
       cudaStreamSynchronize(stream);
    }
  }


  void cuda_AA(void* AA, void* GG, void* diagx, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;   
    cudaStreamSynchronize(stream);
    if(nb > mma_red_s){
       mma_red_s = nb;
       if(mma_bufred != NULL){
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred,nb*sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb*sizeof(real)));
    }
    for (int i = 0; i < (*m); i++){
       for (int j=0; j<(*m);j++){
          mmasumAA_kernel <real> <<<nblcks, nthrds, 0, stream>>>
               ((real*)GG,(real*)diagx, mma_bufred_d, (*n),(*m), i, j);
          CUDA_CHECK(cudaGetLastError());
          mmareduce_kernel<real> <<<1, 1024, 0, stream>>> (mma_bufred_d, nb);
          CUDA_CHECK(cudaGetLastError());
          mma_copy_kernel<<<1, 1, 0, stream>>>((real*)AA, mma_bufred_d, 1,
               i+j*(*m+1));
          CUDA_CHECK(cudaGetLastError());
          cudaStreamSynchronize(stream);
        }
    }
  }
 
 
  void cuda_dx(void* dx,void* delx, void* diagx, void* GG, void* dlambda,
       int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dx_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)dx, (real*)delx,  (real*)diagx,(real*) GG, (real*)dlambda,
         *n,*m);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_dxsi(void* dxsi, void* xsi, void* dx,void* x,
    void* alpha, real*epsi, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dxsi_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)dxsi, (real*)xsi,  (real*)dx,(real*) x,
         (real*)alpha, *epsi,*n);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_deta(void* deta, void* eta, void* dx, void* x,
       void* beta, real* epsi, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    deta_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)deta, (real*)eta, (real*)dx, (real*)x,
         (real*)beta, *epsi, *n);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_rex(void* rex, void* x, void* xlow, void* xupp, void* pij,
       void* p0j, void* qij, void* q0j, void* lambda, void* xsi, void* eta, 
       int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    RexCalculation_kernel<real> <<<nblcks, nthrds, 0,
          (cudaStream_t)glb_cmd_queue >>> ((real*)rex, (real*)x, (real*)xlow,
          (real*)xupp, (real*)pij, (real*)p0j, (real*)qij, (real*)q0j,
          (real*)lambda, (real*)xsi, (real*)eta, *n, *m);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_rey(void* rey, void* c, void* d, void* y, void* lambda, void* mu,
       int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    rey_calculation_kernel<real> <<<nblcks, nthrds, 0,
         (cudaStream_t)glb_cmd_queue >>> ((real*)rey, (real*)c,
         (real*)d, (real*)y, (real*)lambda, (real*)mu, * n);
    CUDA_CHECK(cudaGetLastError());
  }


  /////a_d=b_d*c_d-d
  void cuda_sub2cons(void * a,void * b,void * c, real *d, int * n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    sub2cons_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>
         ((real *) a, (real *) b, (real *) c, *d, *n);
    CUDA_CHECK(cudaGetLastError());
  }


  /////sum(a^2)
  real cuda_norm(void* a, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
    if(nb > mma_red_s){
       mma_red_s = nb;
       if(mma_bufred != NULL){
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred,nb*sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb*sizeof(real)));
    }
    norm_kernel <real> <<<nblcks, nthrds, 0, stream >>> 
         ((real*)a, mma_bufred_d, (*n));
    CUDA_CHECK(cudaGetLastError());
    mmareduce_kernel<real> <<<1, 1024, 0, stream >>> (mma_bufred_d, nb);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
         cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
    return mma_bufred[0];
  }


  void cuda_dely(void* dely, void* c, void* d, void* y, void* lambda,
       real* epsi, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dely_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>> 
         ((real*)dely,(real*)c, (real*)d, (real*)y, (real*)lambda,*epsi, * n);
    CUDA_CHECK(cudaGetLastError());
  }


  real cuda_maxval2(void* a, void* b, real* cons, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;
    if(nb > mma_red_s){
       mma_red_s = nb;
       if(mma_bufred != NULL) {
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred,nb*sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb*sizeof(real)));
    }
    maxval2_kernel <real> <<<nblcks, nthrds, 0,stream >>> 
         ((real*)a, (real*)b, mma_bufred_d, *cons, *n);
    CUDA_CHECK(cudaGetLastError());
    max_reduce_kernel<real> <<<1, 1024, 0,stream >>> (mma_bufred_d, nb);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
         cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
    return mma_bufred[0];
  }


  real cuda_maxval3(void* a, void* b, void* c, real* cons, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;      
    if(nb > mma_red_s){
       mma_red_s = nb;
       if(mma_bufred != NULL) {
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred,nb*sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb*sizeof(real)));
    }
    maxval3_kernel <real> <<<nblcks, nthrds, 0,stream>>> 
         ((real*)a, (real*)b, (real*)c, mma_bufred_d, *cons, *n);
    max_reduce_kernel<real> <<<1, 1024, 0,stream >>> (mma_bufred_d, nb);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
         cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
    return mma_bufred[0];
  }
 
 
  void cuda_kkt_rex(void* rex, void* df0dx, void* dfdx, void* xsi,
       void* eta, void* lambda, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    kkt_rex_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)rex, (real*)df0dx, (real*)dfdx, (real*)xsi,
         (real*)eta, (real*)lambda, *n, *m);
    CUDA_CHECK(cudaGetLastError());
  }

 
  //////a_d=max(b,c*d_d)
  void cuda_maxcons(void* a, real* b, real* c, void* d, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    maxcons_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)a, *b, *c, (real*)d, *n);
    CUDA_CHECK(cudaGetLastError());
  }


  real cuda_lcsc2(void *a, void*b, int *n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const cudaStream_t stream = (cudaStream_t) glb_cmd_queue;  
    if ( nb > mma_red_s){
       mma_red_s = nb;
       if (mma_bufred != NULL) {
          CUDA_CHECK(cudaFreeHost(mma_bufred));
          CUDA_CHECK(cudaFree(mma_bufred_d));        
       }
       CUDA_CHECK(cudaMallocHost(&mma_bufred,nb*sizeof(real)));
       CUDA_CHECK(cudaMalloc(&mma_bufred_d, nb*sizeof(real)));
    }
    glsc2_kernel <real> <<<nblcks, nthrds, 0, stream>>> 
         ((real*)a, (real*)b, mma_bufred_d, (*n));
    CUDA_CHECK(cudaGetLastError());
    mmareduce_kernel<real> <<<1, 1024, 0, stream>>> (mma_bufred_d, nb);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
         cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);
    return mma_bufred[0];
  }


  void cuda_mpisum(void *a, int *n) {
#ifdef HAVE_DEVICE_MPI
    real* temp=(real*)a;
    cudaStreamSynchronize(stream);
    device_mpi_allreduce_inplace(temp, *n, sizeof(real), DEVICE_MPI_SUM);
#endif
  }


  void cuda_add2inv2(void* a, void *b, real* c, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    add2inv2_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)a, (real*) b, *c, *n);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_max2(void* a, real* b, void* c, real* d, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    max2_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*)a, *b, (real*)c,*d, *n);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_updatebb(void* bb, void* dellambda, void* dely,void* d,
       void* mu, void* y, real* delz, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m+1) + 1024 - 1) / 1024, 1, 1);
    updatebb_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*) bb, (real*) dellambda, (real*) dely,(real*) d,
         (real*) mu, (real*) y, *delz, *m);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_updateAA(void* AA, void* globaltmp_mm, void* s, void* lambda,
       void* d, void*mu,void* y,void* a, real* zeta, real* z, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m+1) + 1024 - 1) / 1024, 1, 1);
    updateAA_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*) AA,(real*) globaltmp_mm, (real*) s, (real*) lambda,(real*) d,
         (real*)mu,(real*) y, (real*)a, *zeta, *z, *m);
    CUDA_CHECK(cudaGetLastError());
  }


  void cuda_dy(void* dy, void* dely, void* dlambda,void* d, void* mu, 
       void* y,  int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    dy_kernel <real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue >>>
         ((real*) dy, (real*) dely, (real*) dlambda, (real*) d,
         (real*) mu,(real*) y, *n);
    CUDA_CHECK(cudaGetLastError());
  }
 
}/* extern "C" */
 