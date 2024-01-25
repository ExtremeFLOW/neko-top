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

#include "math_ext_kernel.h"
#include <device/cuda/check.h>
#include <device/device_config.h>
#include <stdio.h>
#include <stdlib.h>

extern "C" {

#include <math/bcknd/device/device_mpi_op.h>
#include <math/bcknd/device/device_mpi_reduce.h>

/** Fortran wrapper for cadd_mask
 * Add a scalar to vector \f$ a_i = a_i + s, for i \in mask \f$
 */
void cuda_cadd_mask(void* a, real* c, int* size, int* mask, int* mask_size) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*mask_size) + 1024 - 1) / 1024, 1, 1);

    cadd_mask_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
        (real*)a, *c, *size, mask, *mask_size);
    CUDA_CHECK(cudaGetLastError());
}
/** Fortran wrapper for invcol1_mask
 * Invert elements of vector \f$ a_i = 1.0 / a_i, for i \in mask \f$
 */
void cuda_invcol1_mask(void* a, int* size, int* mask, int* mask_size) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*mask_size) + 1024 - 1) / 1024, 1, 1);

    invcol1_mask_kernel<real>
        <<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
            (real*)a, *size, mask, *mask_size);
    CUDA_CHECK(cudaGetLastError());
}
/** Fortran wrapper for col2_mask
 * Invert elements of vector \f$ a_i = b_i * c_i, for i \in mask \f$
 */
void cuda_col2_mask(void* a, void* b, int* size, int* mask, int* mask_size) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*mask_size) + 1024 - 1) / 1024, 1, 1);

    col2_mask_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
        (real*)a, (real*)b, *size, mask, *mask_size);
    CUDA_CHECK(cudaGetLastError());
}
/** Fortran wrapper for col3_mask
 * Invert elements of vector \f$ a_i = b_i * c_i, for i \in mask \f$
 */
void cuda_col3_mask(
    void* a, void* b, void* c, int* size, int* mask, int* mask_size) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*mask_size) + 1024 - 1) / 1024, 1, 1);

    col3_mask_kernel<real><<<nblcks, nthrds, 0, (cudaStream_t)glb_cmd_queue>>>(
        (real*)a, (real*)b, (real*)c, *size, mask, *mask_size);
    CUDA_CHECK(cudaGetLastError());
}
}
