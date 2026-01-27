/*
 Copyright (c) 2025, The Neko-TOP Authors
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

#ifndef __ADJOINT_PRS_RES_KERNEL__
#define __ADJOINT_PRS_RES_KERNEL__

template< typename T >
__global__ void adjoint_prs_res_part1_kernel(T * __restrict__ ta1,
                                             T * __restrict__ ta2,
                                             T * __restrict__ ta3,
                                             const T * __restrict__ f_u,
                                             const T * __restrict__ f_v,
                                             const T * __restrict__ f_w,
                                             const T * __restrict__ B,
                                             const T inv_rho,
                                             const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    const T scaled = inv_rho * B[i];
    ta1[i] = f_u[i] * scaled;
    ta2[i] = f_v[i] * scaled;
    ta3[i] = f_w[i] * scaled;
  }
}


template< typename T >
__global__ void adjoint_prs_res_part2_kernel(T * __restrict__ p_res,
                                             const T * __restrict__ wa1,
                                             const T * __restrict__ wa2,
                                             const T * __restrict__ wa3,
                                             const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  for (int i = idx; i < n; i += str) {
    p_res[i] = (-p_res[i]) + (wa1[i] + wa2[i] + wa3[i]);
  }
}

#endif
