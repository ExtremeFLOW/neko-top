/**
 * @file steepest_descent.hip
 * @brief Fortran-callable HIP wrappers of the steepest-descent kernels.
 * @details Scalars are passed by reference and device pointers by value.
 * Every wrapper returns early for an empty range, so that a rank without
 * design variables never launches a zero-block grid.
 * @copyright
 * Copyright (c) 2026, The Neko-TOP Authors
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

// Neko includes
#include <neko/device/device_config.h>
#include <neko/device/hip/check.h>

// Local includes
#include "steepest_descent_kernel.h"

extern "C" {

/**
 * Fortran wrapper for the bound-projected steepest-descent direction (D22).
 * @param s_d Device pointer to the gradient on entry, and to the projected
 * direction on exit.
 * @param x_d Device pointer to the design variables.
 * @param xmin Lower bound of the design variables.
 * @param xmax Upper bound of the design variables.
 * @param n Number of design variables.
 */
void hip_sd_direction(void* s_d, void* x_d, real* xmin, real* xmax,
    int* n) {
    if (*n < 1) return;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sd_direction_kernel<real>),
        nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue,
        (real*)s_d, (real*)x_d, *xmin, *xmax, *n);
    HIP_CHECK(hipGetLastError());
}

/**
 * Fortran wrapper for the move-limit adaptation (D23).
 * @param delta_d Device pointer to the move limits, updated in place.
 * @param x_d Device pointer to the design variables at the current
 * iteration.
 * @param x1_d Device pointer to the design variables one iteration back.
 * @param x2_d Device pointer to the design variables two iterations back.
 * @param dec Factor applied on oscillation.
 * @param inc Factor applied on monotone progress.
 * @param dmin Smallest allowed move limit.
 * @param dmax Largest allowed move limit.
 * @param n Number of design variables.
 */
void hip_sd_update_move_limits(void* delta_d, void* x_d, void* x1_d,
    void* x2_d, real* dec, real* inc, real* dmin, real* dmax, int* n) {
    if (*n < 1) return;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sd_update_move_limits_kernel<real>),
        nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue,
        (real*)delta_d, (real*)x_d, (real*)x1_d, (real*)x2_d,
        *dec, *inc, *dmin, *dmax, *n);
    HIP_CHECK(hipGetLastError());
}

/**
 * Fortran wrapper for the move-limited, box-constrained update (D24).
 * @param x_new_d Device pointer to the updated design variables.
 * @param x_d Device pointer to the current design variables.
 * @param s_d Device pointer to the bound-projected direction.
 * @param delta_d Device pointer to the move limits.
 * @param psi Step length.
 * @param smax Global maximum of the absolute value of the direction.
 * @param xmin Lower bound of the design variables.
 * @param xmax Upper bound of the design variables.
 * @param n Number of design variables.
 */
void hip_sd_update_design(void* x_new_d, void* x_d, void* s_d,
    void* delta_d, real* psi, real* smax, real* xmin, real* xmax, int* n) {
    if (*n < 1) return;

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(sd_update_design_kernel<real>),
        nblcks, nthrds, 0, (hipStream_t)glb_cmd_queue,
        (real*)x_new_d, (real*)x_d, (real*)s_d, (real*)delta_d,
        *psi, *smax, *xmin, *xmax, *n);
    HIP_CHECK(hipGetLastError());
}

}
