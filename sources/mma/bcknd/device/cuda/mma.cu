/**
 * @file mma.cu
 * @copyright
 * Copyright (c) 2025, The Neko-TOP Authors
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 * 
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 * 
 *   * Neither the name of the authors nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

// System includes
#include <stdio.h>
#include <stdlib.h>

// Device includes
#include <cuda_runtime.h>
#include <cusolverDn.h>

// Neko includes
#include <neko/device/device_config.h>
#include <neko/device/cuda/check.h>
#include <neko/math/bcknd/device/device_mpi_reduce.h>
#include <neko/math/bcknd/device/device_mpi_op.h>

// Local includes
#include "mma_kernel.h"


extern "C" {

  int mma_red_s = 0;
  real* mma_bufred = NULL;
  real* mma_bufred_d = NULL;

  void mma_unconstrained_kkt_cuda(void* rex, void* x, void* xmin, void* xmax,
                                  void* df0dx, real* eps, int* n) {
    const int N = *n;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((N + 1024 - 1) / 1024, 1, 1);

    mma_unconstrained_kkt_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
        (real*)rex, (const real*)x, (const real*)xmin, (const real*)xmax,
        (const real*)df0dx, *eps, N);

    CUDA_CHECK(cudaGetLastError());
  }

  void mma_update_hessian_z_cuda(void* Hess, void* a, int* m) {
    const int M = *m;
    // This function is called ONLY if dot(lambda,a) > 0

    const int total = M * M;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((total + 1024 - 1) / 1024, 1, 1);

    mma_update_hessian_z_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
        (real*)Hess, (real*)a, M);

    CUDA_CHECK(cudaGetLastError());
  }

  void mma_prepare_aa_matrix_cuda(void* AA, void* s, void* lambda,
                               void* d, void* mu, void* y,
                               void* a, real* zeta, real* z, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m) + 1024 - 1) / 1024, 1, 1);

    // Launch kernel to prepare AA matrix
    mma_prepare_aa_matrix_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
        (real*)AA, (real*)s, (real*)lambda, (real*)d,
        (real*)mu, (real*)y, (real*)a,
        *zeta, *z, *m);

    CUDA_CHECK(cudaGetLastError());
  }

  void mma_prepare_hessian_cuda(void* Hess, void* y,
                             void* mu, void* lambda, int* m) {
    const int M = *m;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((M + 1024 - 1) / 1024, 1, 1);
    const cudaStream_t stream = (cudaStream_t)glb_cmd_queue;

    // Update diagonal elements
    mma_update_hessian_diagonal_kernel<real><<<nblcks, nthrds, 0, stream>>>(
        (real*)Hess, (real*)y, (real*)mu, (real*)lambda, M);
    CUDA_CHECK(cudaGetLastError());

    // Synchronize to ensure diagonal updates are complete
    CUDA_CHECK(cudaStreamSynchronize(stream));

    // Choose kernel based on problem size
    if (M <= 1024) {
        // Single-block version (fast for small m)
        const dim3 stab_nblcks(1, 1, 1);
        mma_stabilize_hessian_single_kernel<real><<<stab_nblcks, nthrds, 0, stream>>>(
            (real*)Hess, M);
        CUDA_CHECK(cudaGetLastError());
    } else {
        // Multi-block version (for large m)
        // Compute trace on host (simple and reliable)
        real* h_Hess = (real*)malloc(M * sizeof(real));

        // Extract diagonal elements
        for (int i = 0; i < M; i++) {
            CUDA_CHECK(cudaMemcpyAsync(&h_Hess[i],
                                      (real*)Hess + i * M + i,
                                      sizeof(real),
                                      cudaMemcpyDeviceToHost, stream));
        }
        CUDA_CHECK(cudaStreamSynchronize(stream));

        // Compute trace and LM factor
        real trace = 0.0;
        for (int i = 0; i < M; i++) {
            trace += h_Hess[i];
        }
        real lm_factor = fmax(-1.0e-4 * trace / M, 1.0e-7);

        // Apply stabilization in parallel
        mma_stabilize_hessian_multi_kernel<real><<<nblcks, nthrds, 0, stream>>>(
            (real*)Hess, lm_factor, M);
        CUDA_CHECK(cudaGetLastError());

        free(h_Hess);
    }
  }

 void cuda_custom_solver(void* A, void* b, int n, int* info) {
    const cudaStream_t stream = (cudaStream_t)glb_cmd_queue;

    if (n <= 0) {
        *info = -1; // Use CPU fallback
        return;
    }
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(1, 1, 1);

    mma_small_lu_kernel<real><<<nblcks, nthrds, 0, stream>>>(
        (real*)A, (real*)b, n);

    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess) {
        *info = 0; // GPU solver succeeded
    } else {
        *info = -1; // GPU failed
    }
 }

  void cuSOLVER_wrapper(void* A, void* b, int n, int* jj) {
    cusolverDnHandle_t handle;
    cusolverStatus_t status;
    cusolverDnCreate(&handle);

    int lwork;
    double* workspace;
    int* ipiv;
    int* info;  // Device pointer for cuSOLVER info
    int host_info = 0;  // Host variable to store the info

    // Workspace query
    status = cusolverDnDgetrf_bufferSize(handle, n, n, (double*)A, n, &lwork);
    if (status != CUSOLVER_STATUS_SUCCESS) {
        exit(EXIT_FAILURE);
    }
    cudaMalloc(&workspace, lwork * sizeof(double));
    cudaMalloc(&ipiv, n * sizeof(int));
    cudaMalloc(&info, sizeof(int));

    // LU factorization and solve
    cusolverDnDgetrf(handle, n, n, (double*)A, n, workspace, ipiv, info);

    // Copy info from device to host to check if factorization succeeded
    cudaMemcpy(&host_info, info, sizeof(int), cudaMemcpyDeviceToHost);

    if (host_info == 0) {
        // Only solve if factorization was successful
        cusolverDnDgetrs(handle, CUBLAS_OP_N, n, 1, (double*)A, n, ipiv, (double*)b, n, info);
        // Copy the final info value
        cudaMemcpy(&host_info, info, sizeof(int), cudaMemcpyDeviceToHost);
    }


    // Return the actual info value through jj
    *jj = host_info;

    // Cleanup
    cudaFree(workspace);
    cudaFree(ipiv);
    cudaFree(info);
    cusolverDnDestroy(handle);
  }

  void custom_solve_linear_system(void* A, void* b, int n, int* info) {
    const cudaStream_t stream = (cudaStream_t)glb_cmd_queue;

    if (n <= 0) {
        *info = -1; // Use CPU fallback
        return;
    }
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(1, 1, 1);

    mma_small_lu_kernel<real><<<nblcks, nthrds, 0, stream>>>(
        (real*)A, (real*)b, n);

    cudaError_t err = cudaGetLastError();
    if (err == cudaSuccess) {
        *info = 0; // GPU solver succeeded
    } else {
        *info = -1; // GPU failed, use CPU fallback
    }
  }

  void delta_1dbeam_cuda(void* Delta, real* L_total, real* Le,
                       int* offset, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    delta_1dbeam_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
        ((real*)Delta, *L_total, *Le, *offset, *n);
    CUDA_CHECK(cudaGetLastError());
  }

  void cuda_Hess(void* Hess, void* hijx, void* Ljjxinv, int *n, int *m) {
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
           mmasumHess_kernel <real> <<<nblcks, nthrds, 0, stream>>>
                ((real*)hijx,(real*)Ljjxinv, mma_bufred_d, (*n),(*m), i, j);
           CUDA_CHECK(cudaGetLastError());
           mmareduce_kernel<real> <<<1, 1024, 0, stream>>> (mma_bufred_d, nb);
           CUDA_CHECK(cudaGetLastError());
           mma_copy_kernel<<<1, 1, 0, stream>>>((real*)Hess, mma_bufred_d, 1,
                i+j*(*m));
           CUDA_CHECK(cudaGetLastError());
           cudaStreamSynchronize(stream);
         }
     }
  }

  void mma_Ljjxinv_cuda(void* Ljjxinv, void* pjlambda, void* qjlambda, void* x,
     void* low, void* upp, void* alpha, void* beta, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    mma_Ljjxinv_kernel<real> <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
         ((real*)Ljjxinv, (real*)pjlambda, (real*)qjlambda, (real*)x, (real*)low,
         (real*)upp, (real*)alpha, (real*)beta, *n);
    CUDA_CHECK(cudaGetLastError());
  }

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

  void mma_gensub4_cuda(const void* x, const void* low, const void* upp,
                      const void* pij, const void* qij,
                      const int* n, const int* m, void* bi) {

    const int N = *n;
    const int M = *m;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((N + 1023) / 1024, 1, 1);
    const int nb = (N + 1023) / 1024;
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

    real* temp;
    real* bi_d = static_cast<real*>(bi);
    CUDA_CHECK(cudaMalloc(&temp, M * N * sizeof(real)));

    mma_sub4_kernel<real><<<nblcks, nthrds, 0, stream>>>(
        static_cast<const real*>(x),
        static_cast<const real*>(low),
        static_cast<const real*>(upp),
        static_cast<const real*>(pij),
        static_cast<const real*>(qij),
        temp, N, M);

    for (int i = 0; i < M; ++i) {
        mmasum_kernel<real><<<nblcks, nthrds, 0, stream>>>(
            temp, mma_bufred_d, N, M, i);
        CUDA_CHECK(cudaGetLastError());

        mmareduce_kernel<real><<<1, 1024, 0, stream>>>(mma_bufred_d, nb);
        CUDA_CHECK(cudaGetLastError());

        CUDA_CHECK(cudaMemcpyAsync(
            bi_d + i, mma_bufred_d, sizeof(real),
            cudaMemcpyDeviceToDevice, stream));

        CUDA_CHECK(cudaStreamSynchronize(stream));
    }

    CUDA_CHECK(cudaFree(temp));
  }

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
       void* xold2, void* xdiff, real* asydecr, real* asyincr, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    mma_sub2_kernel<real><<<nblcks, nthrds, 0,
         (cudaStream_t)glb_cmd_queue>>>(
         (real*)low, (real*)upp, (real*)x, (real*)xold1,
         (real*)xold2, (real*)xdiff, *asydecr, *asyincr, *n);

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

 //max abs values of input
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


  void cuda_sub2cons(void * a,void * b,void * c, real *d, int * n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    sub2cons_kernel<real><<<nblcks, nthrds, 0,(cudaStream_t) glb_cmd_queue>>>
         ((real *) a, (real *) b, (real *) c, *d, *n);
    CUDA_CHECK(cudaGetLastError());
  }

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
