/**
 * @file RAMP_mapping.cu
 * @copyright
 * Copyright (c) 2024-2025, The Neko-TOP Authors
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

// Neko includes
#include <neko/device/cuda/check.h>
#include <neko/device/device_config.h>

// Local includes
#include "RAMP_mapping_kernel.h"

extern "C" {

/** Fortran wrapper for RAMP (convex down) mapping
 */
void cuda_convex_down_RAMP_mapping_apply(real* f_min, real* f_max, real* q,
    void* X_out_d, void* X_in_d, int* n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    convex_down_RAMP_mapping_apply_kernel<real>
        <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
        (*f_min, *f_max, *q, (real*)X_out_d, (real*)X_in_d, *n);
    CUDA_CHECK(cudaGetLastError());
}

/** Fortran wrapper for RAMP (convex down) chain rule
 */
void cuda_convex_down_RAMP_mapping_apply_backward(real* f_min, real* f_max,
    real* q, void* sens_out_d, void* sens_in_d, void* X_in_d, int* n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);

    convex_down_RAMP_mapping_apply_backward_kernel<real>
        <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>
        (*f_min, *f_max, *q, (real*)sens_out_d,(real*)sens_in_d,
        (real*)X_in_d, *n);
    CUDA_CHECK(cudaGetLastError());
}
}
