/**
 * @file Borrvall_Peterson_mapping_kernel.h
 * @copyright
 * Copyright (c) 2024-2026, The Neko-TOP Authors
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

#ifndef __NEKO_HIP_BORRVALL_PETERSON_MAPPING_KERNELS__
#define __NEKO_HIP_BORRVALL_PETERSON_MAPPING_KERNELS__

/**
 * Device kernel for Borrvall & Peterson mapping
 */
template <typename T>
__global__ void Borrvall_Peterson_mapping_apply_kernel(
    const T f_min, const T f_max, const T q, T* __restrict__ X_out_d,
    T* __restrict__ X_in_d, const int n) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < n; i += str) {
        X_out_d[i] = f_min
            + (f_max - f_min) * X_in_d[i] * (1.0 + q) / (X_in_d[i] + q);
    }
}

/**
 * Device kernel for Borrvall & Peterson chain rule
 */
template <typename T>
__global__ void Borrvall_Peterson_mapping_apply_backward_kernel(
    const T f_min, const T f_max, const T q, T* __restrict__ sens_out_d,
    T* __restrict__ sens_in_d, T* __restrict__ X_in_d, const int n) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < n; i += str) {
        sens_out_d[i] = (f_max - f_min) * q
            * (q + 1.0) / ( (X_in_d[i] + q) * (X_in_d[i] + q)) * sens_in_d[i];
    }
}
#endif // __NEKO_HIP_BORRVALL_PETERSON_MAPPING_KERNELS__
