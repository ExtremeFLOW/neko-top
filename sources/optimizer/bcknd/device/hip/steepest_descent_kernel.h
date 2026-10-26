/**
 * @file steepest_descent_kernel.h
 * @brief HIP kernels of the steepest-descent optimizer.
 * @details Elementwise kernels of the steepest descent method with move
 * limits, Appendix D of da Silva, Aage, Beck & Sigmund (2021), Eqs. D22-D24.
 * They implement the same formulas as `steepest_descent_cpu.f90`.
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

#ifndef __NEKO_TOP_HIP_STEEPEST_DESCENT_KERNELS__
#define __NEKO_TOP_HIP_STEEPEST_DESCENT_KERNELS__

/**
 * Device kernel for the bound-projected steepest-descent direction (D22).
 * On entry s holds the gradient G, on exit S = -G with every component that
 * points out of the box at an active bound set to zero.
 * @param s Gradient on entry, projected direction on exit.
 * @param x Design variables.
 * @param xmin Lower bound of the design variables.
 * @param xmax Upper bound of the design variables.
 * @param n Number of design variables.
 */
template <typename T>
__global__ void sd_direction_kernel(
    T* __restrict__ s, const T* __restrict__ x, const T xmin, const T xmax,
    const int n) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < n; i += str) {
        const T xi = x[i];
        T si = -s[i];
        if ((xi >= xmax && si > T(0.0)) || (xi <= xmin && si < T(0.0))) {
            si = T(0.0);
        }
        s[i] = si;
    }
}

/**
 * Device kernel for the move-limit adaptation (D23).
 * With d = (x - x1) * (x1 - x2), delta is scaled by dec where d < 0, by inc
 * where d > 0, left unchanged where d = 0, and then clamped to [dmin, dmax].
 * @param delta Move limits, updated in place.
 * @param x Design variables at the current iteration.
 * @param x1 Design variables one iteration back.
 * @param x2 Design variables two iterations back.
 * @param dec Factor applied on oscillation.
 * @param inc Factor applied on monotone progress.
 * @param dmin Smallest allowed move limit.
 * @param dmax Largest allowed move limit.
 * @param n Number of design variables.
 */
template <typename T>
__global__ void sd_update_move_limits_kernel(
    T* __restrict__ delta, const T* __restrict__ x, const T* __restrict__ x1,
    const T* __restrict__ x2, const T dec, const T inc, const T dmin,
    const T dmax, const int n) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < n; i += str) {
        const T d = (x[i] - x1[i]) * (x1[i] - x2[i]);
        T di = delta[i];
        if (d < T(0.0)) {
            di = dec * di;
        } else if (d > T(0.0)) {
            di = inc * di;
        }
        delta[i] = min(dmax, max(dmin, di));
    }
}

/**
 * Device kernel for the move-limited, box-constrained update (D24).
 * x_new = min(xmax, max(xmin, min(x + delta, max(x - delta,
 * x + psi * (s / smax))))).
 * @param x_new Updated design variables.
 * @param x Current design variables.
 * @param s Bound-projected steepest-descent direction.
 * @param delta Move limits.
 * @param psi Step length.
 * @param smax Global maximum of the absolute value of s.
 * @param xmin Lower bound of the design variables.
 * @param xmax Upper bound of the design variables.
 * @param n Number of design variables.
 */
template <typename T>
__global__ void sd_update_design_kernel(
    T* __restrict__ x_new, const T* __restrict__ x, const T* __restrict__ s,
    const T* __restrict__ delta, const T psi, const T smax, const T xmin,
    const T xmax, const int n) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    for (int i = idx; i < n; i += str) {
        const T xi = x[i];
        const T di = delta[i];
        const T trial = xi + psi * (s[i] / smax);
        x_new[i] = min(xmax, max(xmin, min(xi + di, max(xi - di, trial))));
    }
}

#endif // __NEKO_TOP_HIP_STEEPEST_DESCENT_KERNELS__
