/**
 * @file mma.hip
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
#include <hip/hip_runtime.h>
#include <hipsolver/hipsolver.h>

// Neko includes
#include <neko/device/device_config.h>
#include <neko/device/hip/check.h>
#include <neko/math/bcknd/device/device_mpi_op.h>

// Local includes
#include "mma_kernel.h"

extern "C" {

int mma_red_s = 0;
real * mma_bufred = NULL;
real * mma_bufred_d = NULL;

void mma_unconstrained_kkt_hip(void* rex, void* x, void* xmin, void* xmax, 
                               void* df0dx, real* eps, int* n) {
    const int N = *n;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((N + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipLaunchKernelGGL(
        mma_unconstrained_kkt_kernel<real>, nblcks, nthrds, 0, stream,
        (real*)rex, (const real*)x, (const real*)xmin, (const real*)xmax, 
        (const real*)df0dx, *eps, N);

    HIP_CHECK(hipGetLastError());
}

void mma_update_hessian_z_hip(void* Hess, void* a, int m) {
    const int M = m;
    const int total = M * M;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((total + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;  // same as CUDA stream

    hipLaunchKernelGGL(
        mma_update_hessian_z_kernel<real>, nblcks, nthrds, 0, stream,
             (real*)Hess, (real*)a, M);
    HIP_CHECK(hipGetLastError());
}

void hipSOLVER_wrapper(void* A, void* b, int n, int* jj) {
    hipsolverDnHandle_t handle;
    hipsolverStatus_t status;
    hipsolverDnCreate(&handle);

    int lwork;
    double *workspace;
    int *ipiv;
    int *info;  // Device pointer for hipSOLVER info
    int host_info = 0;  // Host variable to store the info

    // Workspace query
    status = hipsolverDnDgetrf_bufferSize(handle, n, n, (double*)A, n, &lwork);
    hipMalloc(&workspace, lwork * sizeof(double));
    hipMalloc(&ipiv, n * sizeof(int));
    hipMalloc(&info, sizeof(int));

    // LU factorization and solve
    hipsolverDnDgetrf(handle, n, n, (double*)A, n, workspace, ipiv, info);

    // Copy info from device to host to check if factorization succeeded
    hipMemcpy(&host_info, info, sizeof(int), hipMemcpyDeviceToHost);

    if (host_info == 0) {
        // Only solve if factorization was successful
        hipsolverDnDgetrs(handle, HIPSOLVER_OP_N, n, 1, (double*)A, n, ipiv, (double*)b, n, info);
        // Copy the final info value
        hipMemcpy(&host_info, info, sizeof(int), hipMemcpyDeviceToHost);
    }

    // Return the actual info value through jj
    *jj = host_info;

    // Cleanup
    hipFree(workspace);
    hipFree(ipiv);
    hipFree(info);
    hipsolverDnDestroy(handle);
}

void mma_prepare_aa_matrix_hip(void* AA, void* s, void* lambda,
                                         void* d, void* mu, void* y,
                                         void* a, real zeta, real z, int m) {
    const int M = m;
    const int matrix_size = M + 1;
    const dim3 nthrds(256, 1, 1);
    const dim3 nblcks((M + 256 - 1) / 256, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    // Launch kernel to prepare AA matrix entirely on device
    hipLaunchKernelGGL(mma_prepare_aa_matrix_kernel<real>,
                      nblcks, nthrds, 0, stream,
                      (real*)AA, (real*)s, (real*)lambda, (real*)d,
                      (real*)mu, (real*)y, (real*)a, zeta, z, M);

    HIP_CHECK(hipGetLastError());
}

void mma_prepare_hessian_hip(void* Hess, void* y,
                                       void* mu, void* lambda, int m) {
    const int M = m;
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks((M + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    // Update diagonal elements
    hipLaunchKernelGGL(mma_update_hessian_diagonal_kernel<real>,
                      nblcks, nthrds, 0, stream,
                      (real*)Hess, (real*)y, (real*)mu, (real*)lambda, M);
    HIP_CHECK(hipGetLastError());

    // Synchronize to ensure diagonal updates are complete
    HIP_CHECK(hipStreamSynchronize(stream));

    // Choose kernel based on problem size
    if (M <= 1024) {
        // Single-block version (fast for small m)
        const dim3 stab_nblcks(1, 1, 1);
        hipLaunchKernelGGL(mma_stabilize_hessian_single_kernel<real>,
                          stab_nblcks, nthrds, 0, stream,
                          (real*)Hess, M);
        HIP_CHECK(hipGetLastError());
    } else {
        // Multi-block version (for large m)
        // Compute trace on host (simple and reliable)
        real* h_Hess = (real*)malloc(M * sizeof(real));

        // Extract diagonal elements
        for (int i = 0; i < M; i++) {
            HIP_CHECK(hipMemcpyAsync(&h_Hess[i],
                                   (real*)Hess + i * M + i,
                                   sizeof(real),
                                   hipMemcpyDeviceToHost, stream));
        }
        HIP_CHECK(hipStreamSynchronize(stream));

        // Compute trace and LM factor
        real trace = 0.0;
        for (int i = 0; i < M; i++) {
            trace += h_Hess[i];
        }
        real lm_factor = fmax(-1.0e-4 * trace / M, 1.0e-7);

        // Apply stabilization in parallel
        hipLaunchKernelGGL(mma_stabilize_hessian_multi_kernel<real>,
                          nblcks, nthrds, 0, stream,
                          (real*)Hess, lm_factor, M);
        HIP_CHECK(hipGetLastError());

        free(h_Hess);
    }
}

// Custom linear solver using kernel
extern "C" void hip_custom_solver(void* A, void* b, int n, int* info) {
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    if (n <= 0) {
        *info = -1; // Use CPU fallback
        return;
    }
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(1, 1, 1);

    hipLaunchKernelGGL(mma_small_lu_kernel<real>, nblcks, nthrds, 0, stream,
        (real*)A, (real*)b, n);

    hipError_t err = hipGetLastError();
    if (err == hipSuccess) {
        *info = 0; // GPU solver succeeded
    } else {
        *info = -1; // GPU failed
    }
}

void delta_1dbeam_hip(void* Delta, real* L_total, real* Le,
                       int* offset, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    hipLaunchKernelGGL(delta_1dbeam_kernel<real>,
                      nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue,
                      (real*)Delta, *L_total, *Le, *offset, *n);
    HIP_CHECK(hipGetLastError());
}

void hip_Hess(void* Hess, void* hijx, void* Ljjxinv, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1) / 1024;
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;
    hipStreamSynchronize(stream);

    if (nb > mma_red_s) {
        mma_red_s = nb;
        if (mma_bufred != NULL) {
            HIP_CHECK(hipHostFree(mma_bufred));
            HIP_CHECK(hipFree(mma_bufred_d));
        }
        HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
        HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
    }

    for (int i = 0; i < (*m); i++) {
        for (int j = 0; j < (*m); j++) {
            hipLaunchKernelGGL(mmasumHess_kernel<real>, nblcks, nthrds, 0, stream,
                (real*)hijx, (real*)Ljjxinv, mma_bufred_d, (*n), (*m), i, j);
            HIP_CHECK(hipGetLastError());

            hipLaunchKernelGGL(mmareduce_kernel<real>, dim3(1), dim3(1024), 0, stream,
                mma_bufred_d, nb);
            HIP_CHECK(hipGetLastError());

            hipLaunchKernelGGL(mma_copy_kernel, dim3(1), dim3(1), 0, stream,
                (real*)Hess, mma_bufred_d, 1, i + j * (*m));
            HIP_CHECK(hipGetLastError());

            hipStreamSynchronize(stream);
        }
    }
}

void mma_Ljjxinv_hip(void* Ljjxinv, void* pjlambda, void* qjlambda, void* x,
    void* low, void* upp, void* alpha, void* beta, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    hipLaunchKernelGGL(mma_Ljjxinv_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue, (real*)Ljjxinv, (real*)pjlambda, (real*)qjlambda,
        (real*)x, (real*)low, (real*)upp, (real*)alpha, (real*)beta, *n);
    HIP_CHECK(hipGetLastError());
}

void mma_dipsolvesub1_hip(void* x, void* pjlambda, void* qjlambda, void* low,
    void* upp, void* alpha, void* beta, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    hipLaunchKernelGGL(mma_dipsolvesub1_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue, (real*)x, (real*)pjlambda, (real*)qjlambda,
        (real*)low, (real*)upp, (real*)alpha, (real*)beta, *n);
    HIP_CHECK(hipGetLastError());
}

void mattrans_v_mul_hip(void* output, void* pij, void* lambda, int* m, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    hipLaunchKernelGGL(mattrans_v_mul_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue, (real*)output, (real*)pij, (real*)lambda, *m, *n);
    HIP_CHECK(hipGetLastError());
}

void mma_gensub4_hip(const void* x, const void* low, const void* upp,
                     const void* pij, const void* qij,
                     const int* n, const int* m, void* bi) {

  const int N = *n;
  const int M = *m;

  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks((N + 1023) / 1024, 1, 1);
  const int nb = (N + 1023) / 1024;
  const hipStream_t stream = (hipStream_t)glb_cmd_queue;

  if (nb > mma_red_s) {
    mma_red_s = nb;

    if (mma_bufred != nullptr) {
      HIP_CHECK(hipFreeHost(mma_bufred));
      HIP_CHECK(hipFree(mma_bufred_d));
    }

    HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
    HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
  }

  real* temp;
  real* bi_d = static_cast<real*>(bi);
  HIP_CHECK(hipMalloc(&temp, M * N * sizeof(real)));

  hipLaunchKernelGGL(mma_sub4_kernel<real>, nblcks, nthrds, 0, stream,
                     static_cast<const real*>(x),
                     static_cast<const real*>(low),
                     static_cast<const real*>(upp),
                     static_cast<const real*>(pij),
                     static_cast<const real*>(qij),
                     temp, N, M);

  for (int i = 0; i < M; ++i) {
    hipLaunchKernelGGL(mmasum_kernel<real>, nblcks, nthrds, 0, stream,
                       temp, mma_bufred_d, N, M, i);
    HIP_CHECK(hipGetLastError());

    hipLaunchKernelGGL(mmareduce_kernel<real>, dim3(1), dim3(1024), 0, stream,
                       mma_bufred_d, nb);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpyAsync(
        bi_d + i, mma_bufred_d, sizeof(real),
        hipMemcpyDeviceToDevice, stream));

    HIP_CHECK(hipStreamSynchronize(stream));
  }

  HIP_CHECK(hipFree(temp));
}


void mma_gensub3_hip(void* x, void* df0dx, void* dfdx, void* low,
                     void* upp, void* xmin, void* xmax, void* alpha,
                     void* beta, void* p0j, void* q0j, void* pij,
                     void* qij, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(mma_sub3_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)x, (real*)df0dx, (real*)dfdx, (real*)low,
        (real*)upp, (real*)xmin, (real*)xmax, (real*)alpha,
        (real*)beta, (real*)p0j, (real*)q0j, (real*)pij,
        (real*)qij, *n, *m);

    HIP_CHECK(hipGetLastError());
}

void mma_gensub2_hip(void* low, void* upp, void* x, void* xold1,
                     void* xold2, void* xdiff, real* asydecr,
                     real* asyincr, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(mma_sub2_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)low, (real*)upp, (real*)x, (real*)xold1,
        (real*)xold2, (real*)xdiff, *asydecr, *asyincr, *n);

    HIP_CHECK(hipGetLastError());
}

void mma_gensub1_hip(void* low, void* upp, void* x, void* xmin, void* xmax,
                     real* asyinit, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(mma_sub1_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)low, (real*)upp, (real*)x, (real*)xmin, (real*)xmax,
        *asyinit, *n);

    HIP_CHECK(hipGetLastError());
}

void hip_mma_max(void* xsi, void* x, void* alpha, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(mma_max2_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)xsi, (real*)x, (real*)alpha, *n);

    HIP_CHECK(hipGetLastError());
}

void hip_relambda(void* relambda, void* x, void* xupp, void* xlow,
                  void* pij, void* qij, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = nblcks.x;
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    if (nb > mma_red_s) {
        mma_red_s = nb;
        if (mma_bufred != NULL) {
            HIP_CHECK(hipHostFree(mma_bufred));
            HIP_CHECK(hipFree(mma_bufred_d));
        }
        HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
        HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
    }

    real* temp;
    hipMalloc(&temp, (*n) * (*m) * sizeof(real));

    hipLaunchKernelGGL(relambda_kernel<real>, nblcks, nthrds, 0, stream,
        temp, (real*)x, (real*)xupp, (real*)xlow,
        (real*)pij, (real*)qij, *n, *m);

    for (int i = 0; i < (*m); i++) {
        hipLaunchKernelGGL(mmasum_kernel<real>, nblcks, nthrds, 0, stream,
            temp, mma_bufred_d, (*n), (*m), i);
        HIP_CHECK(hipGetLastError());

        hipLaunchKernelGGL(mmareduce_kernel<real>, dim3(1), dim3(1024), 0,
            stream, mma_bufred_d, nb);
        HIP_CHECK(hipGetLastError());

        hipLaunchKernelGGL(mma_copy_kernel, dim3(1), dim3(1), 0, stream,
            (real*)relambda, mma_bufred_d, 1, i);
        HIP_CHECK(hipGetLastError());

        hipStreamSynchronize(stream);
    }

    hipFree(temp);
}

void hip_sub2cons2(void* a, void* b, void* c, void* d, real* e, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(sub2cons2_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)a, (real*)b, (real*)c, (real*)d, *e, *n);

    HIP_CHECK(hipGetLastError());
}

real hip_maxval(void* a, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const int nb = nblcks.x;
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    if (nb > mma_red_s) {
        mma_red_s = nb;
        if (mma_bufred != NULL) {
            HIP_CHECK(hipHostFree(mma_bufred));
            HIP_CHECK(hipFree(mma_bufred_d));
        }
        HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
        HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
    }

    hipLaunchKernelGGL(maxval_kernel<real>, nblcks, nthrds, 0, stream,
        (real*)a, mma_bufred_d, (*n));
    HIP_CHECK(hipGetLastError());

    hipLaunchKernelGGL(max_reduce_kernel<real>, dim3(1), dim3(1024), 0, stream,
        mma_bufred_d, nb);
    HIP_CHECK(hipGetLastError());

    HIP_CHECK(hipMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
                             hipMemcpyDeviceToHost, stream));
    hipStreamSynchronize(stream);

    return mma_bufred[0];
}


void hip_delx(void* delx, void* x, void* xlow, void* xupp, void* pij,
       void* qij, void* p0j, void* q0j, void* alpha, void* beta, void* lambda,
       real* epsi, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(delx_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)delx, (real*)x, (real*)xlow, (real*)xupp, (real*)pij,
        (real*)qij, (real*)p0j, (real*)q0j, (real*)alpha, (real*)beta,
        (real*)lambda, *epsi, *n, *m);
    HIP_CHECK(hipGetLastError());
}

void hip_GG(void* GG, void* x, void* xlow, void* xupp,
       void* pij, void* qij, int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(GG_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)GG, (real*)x, (real*)xlow, (real*)xupp, (real*)pij,
        (real*)qij, *n, *m);
    HIP_CHECK(hipGetLastError());
}

void hip_diagx(void* diagx, void* x, void* xsi, void* xlow, void* xupp,
       void* p0j, void* q0j, void* pij, void* qij, void* alpha, void* beta,
       void* eta, void* lambda, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(diagx_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)diagx, (real*)x, (real*)xsi, (real*)xlow, (real*)xupp,
        (real*)p0j, (real*)q0j, (real*)pij, (real*)qij, (real*)alpha,
        (real*)beta, (real*)eta, (real*)lambda, *n, *m);
    HIP_CHECK(hipGetLastError());
}

void hip_bb(void* bb, void* GG, void* delx, void* diagx, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipStreamSynchronize(stream);

    if (nb > mma_red_s) {
        mma_red_s = nb;
        if (mma_bufred != NULL) {
            HIP_CHECK(hipHostFree(mma_bufred));
            HIP_CHECK(hipFree(mma_bufred_d));
        }
        HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
        HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
    }

    for (int i = 0; i < (*m); i++) {
        hipLaunchKernelGGL(mmasumbb_kernel<real>, nblcks, nthrds, 0, stream,
            (real*)GG, (real*)delx, (real*)diagx, mma_bufred_d, *n, *m, i);
        HIP_CHECK(hipGetLastError());

        hipLaunchKernelGGL(mmareduce_kernel<real>, 1, 1024, 0, stream,
            mma_bufred_d, nb);
        HIP_CHECK(hipGetLastError());

        hipLaunchKernelGGL(mma_copy_kernel, 1, 1, 0, stream, (real*)bb,
            mma_bufred_d, 1, i);
        HIP_CHECK(hipGetLastError());

        hipStreamSynchronize(stream);
    }
}

void hip_AA(void* AA, void* GG, void* diagx, int *n, int *m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n)+1024 - 1)/ 1024, 1, 1);
    const int nb = ((*n) + 1024 - 1)/ 1024;
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipStreamSynchronize(stream);

    if (nb > mma_red_s) {
        mma_red_s = nb;
        if (mma_bufred != NULL) {
            HIP_CHECK(hipHostFree(mma_bufred));
            HIP_CHECK(hipFree(mma_bufred_d));
        }
        HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
        HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
    }

    for (int i = 0; i < (*m); i++) {
        for (int j = 0; j < (*m); j++) {
            hipLaunchKernelGGL(mmasumAA_kernel<real>, nblcks, nthrds, 0, stream,
                (real*)GG, (real*)diagx, mma_bufred_d, *n, *m, i, j);
            HIP_CHECK(hipGetLastError());

            hipLaunchKernelGGL(mmareduce_kernel<real>, 1, 1024, 0, stream,
                mma_bufred_d, nb);
            HIP_CHECK(hipGetLastError());

            hipLaunchKernelGGL(mma_copy_kernel, 1, 1, 0, stream,
                (real*)AA, mma_bufred_d, 1, i + j * (*m + 1));
            HIP_CHECK(hipGetLastError());

            hipStreamSynchronize(stream);
        }
    }
}

void hip_dx(void* dx, void* delx, void* diagx, void* GG, void* dlambda,
       int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(dx_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)dx, (real*)delx, (real*)diagx, (real*)GG, (real*)dlambda, *n, *m);
    HIP_CHECK(hipGetLastError());
}

void hip_dxsi(void* dxsi, void* xsi, void* dx, void* x,
    void* alpha, real* epsi, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(dxsi_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)dxsi, (real*)xsi, (real*)dx, (real*)x, (real*)alpha, *epsi, *n);
    HIP_CHECK(hipGetLastError());
}

void hip_deta(void* deta, void* eta, void* dx, void* x,
       void* beta, real* epsi, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(deta_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)deta, (real*)eta, (real*)dx, (real*)x, (real*)beta, *epsi, *n);
    HIP_CHECK(hipGetLastError());
}

void hip_rex(void* rex, void* x, void* xlow, void* xupp, void* pij,
       void* p0j, void* qij, void* q0j, void* lambda, void* xsi, void* eta,
       int* n, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(RexCalculation_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)rex, (real*)x, (real*)xlow, (real*)xupp, (real*)pij, (real*)p0j,
        (real*)qij, (real*)q0j, (real*)lambda, (real*)xsi, (real*)eta, *n, *m);
    HIP_CHECK(hipGetLastError());
}

void hip_rey(void* rey, void* c, void* d, void* y, void* lambda, void* mu, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(rey_calculation_kernel<real>, nblcks, nthrds, 0,
        (hipStream_t)glb_cmd_queue,
        (real*)rey, (real*)c, (real*)d, (real*)y, (real*)lambda, (real*)mu, *n);
    HIP_CHECK(hipGetLastError());
}


 // a_d = b_d * c_d - d
void hip_sub2cons(void *a, void *b, void *c, real *d, int *n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  hipLaunchKernelGGL(sub2cons_kernel<real>, nblcks, nthrds, 0,
       (hipStream_t)glb_cmd_queue,
       (real *)a, (real *)b, (real *)c, *d, *n);
  HIP_CHECK(hipGetLastError());
}


// sum(a^2)
real hip_norm(void* a, int* n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const int nb = ((*n) + 1024 - 1) / 1024;
  const hipStream_t stream = (hipStream_t)glb_cmd_queue;

  if (nb > mma_red_s) {
    mma_red_s = nb;
    if (mma_bufred != NULL) {
      HIP_CHECK(hipFreeHost(mma_bufred));
      HIP_CHECK(hipFree(mma_bufred_d));
    }
    HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
    HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
  }

  hipLaunchKernelGGL(norm_kernel<real>, nblcks, nthrds, 0, stream,
       (real*)a, mma_bufred_d, (*n));
  HIP_CHECK(hipGetLastError());

  hipLaunchKernelGGL(mmareduce_kernel<real>, dim3(1), dim3(1024), 0, stream,
       mma_bufred_d, nb);
  HIP_CHECK(hipGetLastError());

  HIP_CHECK(hipMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
       hipMemcpyDeviceToHost, stream));

  hipStreamSynchronize(stream);

  return mma_bufred[0];
}


void hip_dely(void* dely, void* c, void* d, void* y, void* lambda,
       real* epsi, int* n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  hipLaunchKernelGGL(dely_kernel<real>, nblcks, nthrds, 0,
       (hipStream_t)glb_cmd_queue,
       (real*)dely, (real*)c, (real*)d, (real*)y, (real*)lambda, *epsi, *n);
  HIP_CHECK(hipGetLastError());
}


real hip_maxval2(void* a, void* b, real* cons, int* n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const int nb = ((*n) + 1024 - 1) / 1024;
  const hipStream_t stream = (hipStream_t)glb_cmd_queue;

  if (nb > mma_red_s) {
    mma_red_s = nb;
    if (mma_bufred != NULL) {
      HIP_CHECK(hipFreeHost(mma_bufred));
      HIP_CHECK(hipFree(mma_bufred_d));
    }
    HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
    HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
  }

  hipLaunchKernelGGL(maxval2_kernel<real>, nblcks, nthrds, 0, stream,
       (real*)a, (real*)b, mma_bufred_d, *cons, *n);
  HIP_CHECK(hipGetLastError());

  hipLaunchKernelGGL(max_reduce_kernel<real>, dim3(1), dim3(1024), 0, stream,
       mma_bufred_d, nb);
  HIP_CHECK(hipGetLastError());

  HIP_CHECK(hipMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
       hipMemcpyDeviceToHost, stream));

  hipStreamSynchronize(stream);

  return mma_bufred[0];
}


real hip_maxval3(void* a, void* b, void* c, real* cons, int* n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const int nb = ((*n) + 1024 - 1) / 1024;
  const hipStream_t stream = (hipStream_t)glb_cmd_queue;

  if (nb > mma_red_s) {
    mma_red_s = nb;
    if (mma_bufred != NULL) {
      HIP_CHECK(hipFreeHost(mma_bufred));
      HIP_CHECK(hipFree(mma_bufred_d));
    }
    HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
    HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
  }

  hipLaunchKernelGGL(maxval3_kernel<real>, nblcks, nthrds, 0, stream,
       (real*)a, (real*)b, (real*)c, mma_bufred_d, *cons, *n);
  hipLaunchKernelGGL(max_reduce_kernel<real>, dim3(1), dim3(1024), 0, stream,
       mma_bufred_d, nb);
  HIP_CHECK(hipGetLastError());

  HIP_CHECK(hipMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
       hipMemcpyDeviceToHost, stream));

  hipStreamSynchronize(stream);

  return mma_bufred[0];
}


void hip_kkt_rex(void* rex, void* df0dx, void* dfdx, void* xsi,
       void* eta, void* lambda, int* n, int* m) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  hipLaunchKernelGGL(kkt_rex_kernel<real>, nblcks, nthrds, 0,
       (hipStream_t)glb_cmd_queue,
       (real*)rex, (real*)df0dx, (real*)dfdx, (real*)xsi,
       (real*)eta, (real*)lambda, *n, *m);
  HIP_CHECK(hipGetLastError());
}


// a_d = max(b, c * d_d)
void hip_maxcons(void* a, real* b, real* c, void* d, int* n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  hipLaunchKernelGGL(maxcons_kernel<real>, nblcks, nthrds, 0,
       (hipStream_t)glb_cmd_queue,
       (real*)a, *b, *c, (real*)d, *n);
  HIP_CHECK(hipGetLastError());
}


real hip_lcsc2(void *a, void*b, int *n) {
  const dim3 nthrds(1024, 1, 1);
  const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
  const int nb = ((*n) + 1024 - 1) / 1024;
  const hipStream_t stream = (hipStream_t)glb_cmd_queue;

  if (nb > mma_red_s) {
    mma_red_s = nb;
    if (mma_bufred != NULL) {
      HIP_CHECK(hipFreeHost(mma_bufred));
      HIP_CHECK(hipFree(mma_bufred_d));
    }
    HIP_CHECK(hipHostMalloc(&mma_bufred, nb * sizeof(real)));
    HIP_CHECK(hipMalloc(&mma_bufred_d, nb * sizeof(real)));
  }

  hipLaunchKernelGGL(glsc2_kernel<real>, nblcks, nthrds, 0, stream,
       (real*)a, (real*)b, mma_bufred_d, (*n));
  HIP_CHECK(hipGetLastError());

  hipLaunchKernelGGL(mmareduce_kernel<real>, dim3(1), dim3(1024), 0, stream,
       mma_bufred_d, nb);
  HIP_CHECK(hipGetLastError());

  HIP_CHECK(hipMemcpyAsync(mma_bufred, mma_bufred_d, sizeof(real),
       hipMemcpyDeviceToHost, stream));

  hipStreamSynchronize(stream);

  return mma_bufred[0];
}


void hip_mpisum(void *a, int *n) {
#ifdef HAVE_DEVICE_MPI
  real* temp = (real*)a;
  hipStreamSynchronize(stream);
  device_mpi_allreduce_inplace(temp, *n, sizeof(real), DEVICE_MPI_SUM);
#endif
}


void hip_add2inv2(void* a, void* b, real* c, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipLaunchKernelGGL(add2inv2_kernel<real>, nblcks, nthrds, 0, stream,
                       (real*)a, (real*)b, *c, *n);
    HIP_CHECK(hipGetLastError());
}

void hip_max2(void* a, real* b, void* c, real* d, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipLaunchKernelGGL(max2_kernel<real>, nblcks, nthrds, 0, stream,
                       (real*)a, *b, (real*)c, *d, *n);
    HIP_CHECK(hipGetLastError());
}

void hip_updatebb(void* bb, void* dellambda, void* dely, void* d,
                  void* mu, void* y, real* delz, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m + 1) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipLaunchKernelGGL(updatebb_kernel<real>, nblcks, nthrds, 0, stream,
                       (real*)bb, (real*)dellambda, (real*)dely, (real*)d,
                       (real*)mu, (real*)y, *delz, *m);
    HIP_CHECK(hipGetLastError());
}

void hip_updateAA(void* AA, void* globaltmp_mm, void* s, void* lambda,
                  void* d, void* mu, void* y, void* a,
                  real* zeta, real* z, int* m) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*m + 1) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipLaunchKernelGGL(updateAA_kernel<real>, nblcks, nthrds, 0, stream,
                       (real*)AA, (real*)globaltmp_mm, (real*)s,
                       (real*)lambda, (real*)d, (real*)mu,
                       (real*)y, (real*)a, *zeta, *z, *m);
    HIP_CHECK(hipGetLastError());
}

void hip_dy(void* dy, void* dely, void* dlambda, void* d,
            void* mu, void* y, int* n) {
    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t)glb_cmd_queue;

    hipLaunchKernelGGL(dy_kernel<real>, nblcks, nthrds, 0, stream,
                       (real*)dy, (real*)dely, (real*)dlambda, (real*)d,
                       (real*)mu, (real*)y, *n);
    HIP_CHECK(hipGetLastError());
}
}
