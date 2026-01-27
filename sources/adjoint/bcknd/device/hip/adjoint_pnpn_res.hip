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

#include <device/device_config.h>
#include <device/hip/check.h>

#include "adjoint_prs_res_kernel.h"
#include "adjoint_vel_res_update_kernel.h"

extern "C" {

  void adjoint_prs_res_part1_hip(void *ta1, void *ta2, void *ta3,
                                 void *f_u, void *f_v, void *f_w,
                                 void *B, real *inv_rho, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hipLaunchKernelGGL(adjoint_prs_res_part1_kernel<real>, nblcks, nthrds, 0, stream,
                       (real *) ta1, (real *) ta2, (real *) ta3,
                       (real *) f_u, (real *) f_v, (real *) f_w,
                       (real *) B, *inv_rho, *n);
    HIP_CHECK(hipGetLastError());
  }

  void adjoint_prs_res_part2_hip(void *p_res, void *wa1, void *wa2,
                                 void *wa3, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hipLaunchKernelGGL(adjoint_prs_res_part2_kernel<real>, nblcks, nthrds, 0, stream,
                       (real *) p_res, (real *) wa1, (real *) wa2,
                       (real *) wa3, *n);
    HIP_CHECK(hipGetLastError());
  }

  void adjoint_vel_res_update_hip(void *u_res, void *v_res, void *w_res,
                                  void *f_u, void *f_v, void *f_w, int *n) {

    const dim3 nthrds(1024, 1, 1);
    const dim3 nblcks(((*n) + 1024 - 1) / 1024, 1, 1);
    const hipStream_t stream = (hipStream_t) glb_cmd_queue;

    hipLaunchKernelGGL(adjoint_vel_res_update_kernel<real>, nblcks, nthrds, 0, stream,
                       (real *) u_res, (real *) v_res, (real *) w_res,
                       (real *) f_u, (real *) f_v, (real *) f_w, *n);
    HIP_CHECK(hipGetLastError());
  }

}
