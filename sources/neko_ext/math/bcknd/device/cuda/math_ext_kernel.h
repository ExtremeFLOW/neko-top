/*
 Copyright (c) 2021-2023, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef __NEKO_CUDA_MATH_EXT_KERNELS__
#define __NEKO_CUDA_MATH_EXT_KERNELS__

/**
 * Device kernel for cadd_mask
 */
template <typename T>
__global__ void cadd_mask_kernel(
    T* __restrict__ a, const T c, const int size, int* __restrict__ mask,
    const int mask_size) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < mask_size; i += str) { a[mask[i]] = a[mask[i]] + c; }
}

/**
 * Device kernel for invcol1_mask
 */
template <typename T>
__global__ void invcol1_mask_kernel(
    T* __restrict__ a, const int size, int* __restrict__ mask,
    const int mask_size) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < mask_size; i += str) {
        a[mask[i]] = 1.0 / a[mask[i]];
    }
}

/**
 * Device kernel for col2_mask
 */
template <typename T>
__global__ void col2_mask_kernel(
    T* __restrict__ a, T* __restrict__ b, const int size,
    int* __restrict__ mask, const int mask_size) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < mask_size; i += str) {
        a[mask[i]] = a[mask[i]] * b[mask[i]];
    }
}

/**
 * Device kernel for col3_mask
 */
template <typename T>
__global__ void col3_mask_kernel(
    T* __restrict__ a, T* __restrict__ b, T* __restrict__ c, const int size,
    int* __restrict__ mask, const int mask_size) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < mask_size; i += str) {
        a[mask[i]] = b[mask[i]] * c[mask[i]];
    }
}

#endif // __NEKO_CUDA_MATH_EXT_KERNELS__