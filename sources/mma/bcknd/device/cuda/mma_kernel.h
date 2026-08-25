/**
 * @file mma_kernel.h
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

#ifndef MMA_CUDA_KERNEL_H
#define MMA_CUDA_KERNEL_H

template <typename T>
__global__ void mma_unconstrained_kkt_kernel(
    T* __restrict__ rex,
    const T* __restrict__ x,
    const T* __restrict__ xmin,
    const T* __restrict__ xmax,
    const T* __restrict__ df0dx,
    const T eps,
    const int n)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j >= n) return;

    if (x[j] - xmin[j] <= eps) {
        // At lower bound: residual is non-zero only if gradient is negative
        rex[j] = (df0dx[j] < (T)0.0) ? df0dx[j] : (T)0.0;
    } 
    else if (xmax[j] - x[j] <= eps) {
        // At upper bound: residual is non-zero only if gradient is positive
        rex[j] = (df0dx[j] > (T)0.0) ? df0dx[j] : (T)0.0;
    } 
    else {
        // Strictly internal
        rex[j] = df0dx[j];
    }
}
	
template <typename T>
__global__ void mma_update_hessian_z_kernel(
    T* __restrict__ Hess,
    const T* __restrict__ a,
    const int m)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total = m * m;

    if (tid >= total) return;

    int i = tid % m;
    int j = tid / m;

    Hess[tid] -= a[i] * a[j];
}

template<typename T>
__global__ void mma_prepare_aa_matrix_kernel(T* __restrict__ AA,
    const T* __restrict__ s, const T* __restrict__ lambda,
    const T* __restrict__ d, const T* __restrict__ mu,
    const T* __restrict__ y, const T* __restrict__ a,
    const T zeta, const T z, const int m) {
  const int tj = blockIdx.x * blockDim.x + threadIdx.x;
  const int matrix_size = m + 1;

  if (tj >= m) return;
  AA[tj * matrix_size + tj] += s[tj] / lambda[tj] +
       (T)1.0 / (d[tj] + mu[tj] / y[tj]);
  AA[tj * matrix_size + m] = a[tj]; // column m+1
  AA[m * matrix_size + tj] = a[tj];// row m+1

  // Only first thread updates the bottom-right  corner element.
  if (tj == 0)
    AA[m * matrix_size + m] = -zeta / z;
}


//Update Hessian diagonal elements (y contributions for dip subsolve)
template<typename T>
__global__ void mma_update_hessian_diagonal_kernel(T* __restrict__ Hess,
     const T* __restrict__ y, const T* __restrict__ mu, 
     const T* __restrict__ lambda, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= m) return;

  T diag = Hess[tj * m + tj];
  // Contribution from y terms (inactive constraints)
  if (y[tj] > (T)0.0) {
     diag -= (T)1.0;
  }

  // Contribution from -mu/lambda (eq 10)
  diag -= mu[tj] / lambda[tj];
  Hess[tj * m + tj] = diag;
}

// Levenberg-Marquardt algorithm (heuristically)
// Single-block version for m <= 1024
template<typename T>
__global__ void mma_stabilize_hessian_single_kernel(T* __restrict__ Hess, const int m) {
    const int tid = threadIdx.x;

    // Single thread computes trace and LM factor
    if (tid == 0) {
        T trace = (T)0.0;
        for (int j = 0; j < m; j++) {
            trace += Hess[j * m + j];
        }
        T lm_factor = max((T)(-1.0e-4) * trace / m, (T)1.0e-7);

        // Apply to all diagonal elements
        for (int j = 0; j < m; j++) {
            Hess[j * m + j] -= lm_factor;
        }
    }
}

// Levenberg-Marquardt algorithm (heuristically)
// Multi-block version for m > 1024
template<typename T>
__global__ void mma_stabilize_hessian_multi_kernel(T* __restrict__ Hess, const T lm_factor, const int m) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < m) {
        Hess[i * m + i] -= lm_factor;
    }
}

// Small linear solver for MMA (n <= 50)
template<typename T>
__global__ void mma_small_lu_kernel(T* __restrict__ A, T* __restrict__ b, const int n) {
    const int tid = threadIdx.x;

    // Handle 1x1 case
    if (n == 1) {
        if (tid == 0 && abs(A[0]) > (T)1e-12) b[0] /= A[0];
        return;
    }

    // LU decomposition with partial pivoting
    for (int k = 0; k < n; k++) {
        // Pivoting - single thread
        if (tid == 0) {
            int max_row = k;
            T max_val = abs(A[k * n + k]);
            for (int i = k + 1; i < n; i++) {
                T val = abs(A[i * n + k]);
                if (val > max_val) {
                    max_val = val;
                    max_row = i;
                }
            }
            if (max_val > (T)1e-12 && max_row != k) {
                // Swap rows
                for (int j = k; j < n; j++) {
                    T temp = A[k * n + j];
                    A[k * n + j] = A[max_row * n + j];
                    A[max_row * n + j] = temp;
                }
                // Swap rhs
                T temp_b = b[k];
                b[k] = b[max_row];
                b[max_row] = temp_b;
            }
        }
        __syncthreads();

        // Parallel elimination
        T diag = A[k * n + k];
        if (abs(diag) > (T)1e-12) {
            for (int i = tid + k + 1; i < n; i += blockDim.x) {
                T factor = A[i * n + k] / diag;
                A[i * n + k] = factor;
                for (int j = k + 1; j < n; j++) {
                    A[i * n + j] -= factor * A[k * n + j];
                }
            }
        }
        __syncthreads();
    }

    // Parallel forward substitution
    for (int i = tid; i < n; i += blockDim.x) {
        T sum = b[i];
        for (int j = 0; j < i; j++) {
            sum -= A[i * n + j] * b[j];
        }
        b[i] = sum;
    }
    __syncthreads();

    // Parallel backward substitution
    for (int i = n - 1 - tid; i >= 0; i -= blockDim.x) {
        if (i >= 0) {
            T sum = b[i];
            for (int j = i + 1; j < n; j++) {
                sum -= A[i * n + j] * b[j];
            }
            if (abs(A[i * n + i]) > (T)1e-12) {
                b[i] = sum / A[i * n + i];
            }
        }
    }
    __syncthreads();
}

template <typename T>
__global__ void delta_1dbeam_kernel(T* __restrict__ Delta,
     const T L_total, const T Le, const int offset, const int n) {
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  if (k >= n) return;

  // Convert to 1-based indexing for the calculation
  int idx = k + 1;

  // Calculate first term: (L_total - Le*(offset+idx-1))^3
  T x1 = L_total - Le * static_cast<T>(offset + idx - 1);
  T term1 = x1 * x1 * x1;

  // Calculate second term: (L_total - Le*(offset+idx))^3
  T x2 = L_total - Le * static_cast<T>(offset + idx);
  T term2 = x2 * x2 * x2;

  // Final result
  Delta[k] = (term1 - term2) / static_cast<T>(3.0);
}

template <typename T>
__global__ void mma_Ljjxinv_kernel(T* __restrict__ Ljjxinv,
     const T* __restrict__ pjlambda, const T* __restrict__ qjlambda,
     const T* __restrict__ x, const T* __restrict__ low, const T* __restrict__ upp,
     const T* __restrict__ alpha, const T* __restrict__ beta,
     const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  // Load inputs into registers
  const T xt    = x[tj];
  const T low_j = low[tj];
  const T upp_j = upp[tj];
  const T pj    = pjlambda[tj];
  const T qj    = qjlambda[tj];

  // Precompute reused differences
  const T diff_u = upp_j - xt;
  const T diff_l = xt - low_j;

  // Cube once (avoiding pow for speed)
  const T diff_u3 = diff_u * diff_u * diff_u;
  const T diff_l3 = diff_l * diff_l * diff_l;

  // Compute inverse value safely
  T denom = 2.0 * pj / diff_u3 + 2.0 * qj / diff_l3;
  T val = -1.0 / denom;

  // Mask out active primal constraints
  bool active = (fabs(xt - alpha[tj]) <= T(1e-16)) ||
              (fabs(xt - beta[tj])  <= T(1e-16));

  Ljjxinv[tj] = active ? T(0.0) : val;
}

template <typename T>
__global__ void mma_dipsolvesub1_kernel(T* __restrict__ x,
     const T* __restrict__ pjlambda, const T* __restrict__ qjlambda,
     const T* __restrict__ low, const T* __restrict__ upp,
     const T* __restrict__ alpha, const T* __restrict__ beta,
     const int n) {

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  // Load inputs into registers
  const T pj_sqrt = sqrt(pjlambda[tj]);
  const T qj_sqrt = sqrt(qjlambda[tj]);
  const T denom   = pj_sqrt + qj_sqrt;

  const T low_j   = low[tj];
  const T upp_j   = upp[tj];
  const T alpha_j = alpha[tj];
  const T beta_j  = beta[tj];

  // Weighted average
  const T val = (pj_sqrt * low_j + qj_sqrt * upp_j) / denom;

  // Clamp x between alpha and beta using branchless min/max
  x[tj] = fmin(fmax(val, alpha_j), beta_j);
}

template <typename T>
__global__ void mattrans_v_mul_kernel(T* __restrict__ output,
     const T* __restrict__ pij, const T* __restrict__ lambda,
     const int m, const int n) {

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  T acc = 0.0;

  // Use register accumulation
  for (int i = 0; i < m; ++i) {
    acc += pij[i + tj * m] * lambda[i];
  }

  output[tj] = acc;
}

template <typename T>
__global__ void mma_sub1_kernel(
    T* __restrict__ xlow,
    T* __restrict__ xupp,
    const T* __restrict__ x,
    const T* __restrict__ xmin,
    const T* __restrict__ xmax,
    const T asyinit,
    const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  // Load values once into registers (global memory is slow)
  const T xt    = x[tj];
  const T xmin_j = xmin[tj];
  const T xmax_j = xmax[tj];

  // Reuse xgap calculation
  const T xgap = xmax_j - xmin_j;
  const T offset = asyinit * xgap;

  xlow[tj] = xt - offset;
  xupp[tj] = xt + offset;
}



template< typename T >
__global__ void mma_sub2_kernel(T* __restrict__ low, T* __restrict__ upp,
     const T* __restrict__ x, const T* __restrict__ xold1,
     const T* __restrict__ xold2, const T* __restrict__ xdiff,
     const T asydecr, const T asyincr, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  // Load data into registers for faster accessing compare to global memory
  // when accessing repeatedly)
  const T xval     = x[tj];
  const T xold1val = xold1[tj];
  const T xold2val = xold2[tj];
  const T lowval   = low[tj];
  const T uppval   = upp[tj];
  const T xdiffval = xdiff[tj];

  // Compute the product
  const T prod = (xval - xold1val) * (xold1val - xold2val);

  // Compute asy_factor without branching
  T asy_factor = (prod < T(0)) ? asydecr :
                 (prod > T(0)) ? asyincr : T(1);

  // Update low and upp using fma (fused multiply-add) for numerical stability
  T new_low = fma(-asy_factor, (xold1val - lowval), xval);
  T new_upp = fma(asy_factor,  (uppval - xold1val), xval);

  // Apply bounds
  new_low = max(new_low, xval - T(10.0) * xdiffval);
  new_low = min(new_low, xval - T(0.01) * xdiffval);

  new_upp = min(new_upp, xval + T(10.0) * xdiffval);
  new_upp = max(new_upp, xval + T(0.01) * xdiffval);

  // Write results back
  low[tj] = new_low;
  upp[tj] = new_upp;
}

template <typename T>
__global__ void mma_sub3_kernel( const T* __restrict__ x,
    const T* __restrict__ df0dx, const T* __restrict__ dfdx,
    T* __restrict__ low, T* __restrict__ upp, const T* __restrict__ xmin,
    const T* __restrict__ xmax, T* __restrict__ alpha, T* __restrict__ beta,
    T* __restrict__ p0j, T* __restrict__ q0j, T* __restrict__ pij,
    T* __restrict__ qij, const int n, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  // Load into registers once
  const T xt    = x[tj];
  const T xmin_j = xmin[tj];
  const T xmax_j = xmax[tj];
  const T low_j = low[tj];
  const T upp_j = upp[tj];
  const T df0 = df0dx[tj];
  const T xgap = xmax_j - xmin_j;

  // Clamp helpers
  const T half_xgap = 0.5 * xgap;
  const T tenth_low_diff = T(0.1) * (xt - low_j);
  const T tenth_upp_diff = T(0.1) * (upp_j - xt);

  // Compute alpha and beta with fused max/min and fewer calls
  T alpha_val = max(max(xmin_j, low_j + tenth_low_diff), xt - half_xgap);
  T beta_val = min(min(xmax_j, upp_j - tenth_upp_diff), xt + half_xgap);

  alpha[tj] = alpha_val;
  beta[tj] = beta_val;

  // Avoid multiple pow calls, compute once
  const T upp_minus_x = upp_j - xt;
  const T x_minus_low = xt - low_j;
  const T upp_minus_x_sq = upp_minus_x * upp_minus_x;
  const T x_minus_low_sq = x_minus_low * x_minus_low;

  // Small epsilon for numerical stability
  const T eps = T(1e-5);
  const T inv_xgap = T(1.0) / max(eps, xgap);

  // Compute terms reused multiple times
  const T max_df0_pos = max(df0, T(0));
  const T max_df0_neg = max(-df0, T(0));

  p0j[tj] = upp_minus_x_sq * (T(1.001) * max_df0_pos +
       T(0.001) * max_df0_neg + eps * inv_xgap);
  q0j[tj] = x_minus_low_sq * (T(0.001) * max_df0_pos +
       T(1.001) * max_df0_neg + eps * inv_xgap);

  // Loop over m for pij and qij
  for (int i = 0; i < m; ++i) {
    int idx = i + tj * m;

    T dfdx_val = dfdx[idx];
    T max_pos = max(dfdx_val, T(0));
    T max_neg = max(-dfdx_val, T(0));

    pij[idx] = upp_minus_x_sq * (T(1.001) * max_pos +
         T(0.001) * max_neg + eps * inv_xgap);
    qij[idx] = x_minus_low_sq * (T(0.001) * max_pos +
         T(1.001) * max_neg + eps * inv_xgap);
  }
}

template <typename T>
__global__ void mma_sub4_kernel(
    const T* __restrict__ x,
    const T* __restrict__ low,
    const T* __restrict__ upp,
    const T* __restrict__ pij,
    const T* __restrict__ qij,
    T* __restrict__ temp,
    const int n,
    const int m) {

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  // Register caching
  const T xt     = x[tj];
  const T low_j  = low[tj];
  const T upp_j  = upp[tj];

  const T denom_upp = upp_j - xt;
  const T denom_low = xt - low_j;

  const T eps = T(1e-12);  // Precision-dependent epsilon
  const T inv_denom_upp = T(1) / max(denom_upp, eps);
  const T inv_denom_low = T(1) / max(denom_low, eps);

  const int base_idx = tj * m;

  for (int i = 0; i < m; ++i) {
    int idx = base_idx + i;
    temp[idx] = pij[idx] * inv_denom_upp + qij[idx] * inv_denom_low;
  }
}


template <typename T>
__global__ void mma_max2_kernel(
    T* __restrict__ xsi,
    const T* __restrict__ x,
    const T* __restrict__ alpha,
    const int n) {

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  const T eps = T(1e-12);
  T denom = max(x[tj] - alpha[tj], eps);
  xsi[tj] = max(T(1), T(1) / denom);
}




template <typename T>
__global__ void relambda_kernel(
    T* __restrict__ temp,
    const T* __restrict__ x,
    const T* __restrict__ xupp,
    const T* __restrict__ xlow,
    const T* __restrict__ pij,
    const T* __restrict__ qij,
    const int n,
    const int m) {

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj >= n) return;

  const T xt = x[tj];
  const T xup = xupp[tj];
  const T xlo = xlow[tj];

  // Prevent divide-by-zero using small epsilon
  const T eps = T(1e-12);
  const T inv_denom_upp = T(1) / max(xup - xt, eps);
  const T inv_denom_low = T(1) / max(xt - xlo, eps);

  int base_idx = tj * m;
  for (int i = 0; i < m; ++i) {
    int idx = base_idx + i;
    temp[idx] = pij[idx] * inv_denom_upp + qij[idx] * inv_denom_low;
  }
}


template <typename T>
__global__ void sub2cons2_kernel(
    T* __restrict__ a,
    const T* __restrict__ b,
    const T* __restrict__ c,
    const T* __restrict__ d,
    const T e,
    const int n) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    const T bt = b[idx];
    const T ct = c[idx];
    const T dt = d[idx];

    a[idx] = bt * (ct - dt) - e;
}


template< typename T>
__inline__ __device__ T max_reduce_warp(T val) {
  val = max(val, __shfl_down_sync(0xffffffff, val, 16));
  val = max(val, __shfl_down_sync(0xffffffff, val, 8));
  val = max(val, __shfl_down_sync(0xffffffff, val, 4));
  val = max(val, __shfl_down_sync(0xffffffff, val, 2));
  val = max(val, __shfl_down_sync(0xffffffff, val, 1));
  return val;
}



template <typename T>
__global__ void maxval_kernel(const T* __restrict__ a, T* temp, const int n) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;

    const unsigned int lane = threadIdx.x % warpSize;
    const unsigned int warp_id = threadIdx.x / warpSize;

    // Shared memory to store warp-level maxima
    __shared__ T shared[32];  // Assumes max 1024 threads per block

    T local_max = T(0);
    for (int i = idx; i < n; i += stride) {
        local_max = max(local_max, abs(a[i]));
    }

    // Warp-level reduction
    local_max = max_reduce_warp<T>(local_max);

    // Write warp-level results to shared memory
    if (lane == 0) {
        shared[warp_id] = local_max;
    }
    __syncthreads();

    // Let the first warp reduce shared values
    local_max = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : T(0);

    if (warp_id == 0) {
        local_max = max_reduce_warp<T>(local_max);
    }

    if (threadIdx.x == 0) {
        temp[blockIdx.x] = local_max;
    }
}


template <typename T>
__global__ void max_reduce_kernel(T* __restrict__ bufred, const int n) {

    T maxval = T(0);
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;

    // Grid-stride loop to cover all elements
    for (int i = idx; i < n; i += stride) {
        maxval = max(maxval, bufred[i]);
    }

    __shared__ T shared[32];  // One slot per warp (max 1024 threads/block)

    unsigned int lane = threadIdx.x % warpSize;
    unsigned int wid = threadIdx.x / warpSize;

    // Warp-level max reduction
    maxval = max_reduce_warp<T>(maxval);

    // Store each warp's max value in shared memory
    if (lane == 0) {
        shared[wid] = maxval;
    }
    __syncthreads();

    // Let the first warp reduce the warp-level maxima
    maxval = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : T(0);
    if (wid == 0) {
        maxval = max_reduce_warp<T>(maxval);
    }

    // Write the block's maximum value back to bufred[blockIdx.x]
    if (threadIdx.x == 0) {
        bufred[blockIdx.x] = maxval;
    }
}

template <typename T>
__global__ void delx_kernel(
    T* __restrict__ delx,
    const T* __restrict__ x,
    const T* __restrict__ xlow,
    const T* __restrict__ xupp,
    const T* __restrict__ pij,
    const T* __restrict__ qij,
    const T* __restrict__ p0j,
    const T* __restrict__ q0j,
    const T* __restrict__ alpha,
    const T* __restrict__ beta,
    const T* __restrict__ lambda,
    const T epsi,
    const int n,
    const int m)
{
    int tj = blockIdx.x * blockDim.x + threadIdx.x;
    if (tj < n) {
        T xt = x[tj];
        T xlow_j = xlow[tj];
        T xupp_j = xupp[tj];
        T alpha_j = alpha[tj];
        T beta_j = beta[tj];

        // Precompute denominators squared for better performance
        T denom_low = xt - xlow_j;
        T denom_upp = xupp_j - xt;
        T denom_alpha = xt - alpha_j;
        T denom_beta = beta_j - xt;

        const T small_eps = T(1e-12);
        denom_low = (abs(denom_low) < small_eps) ? small_eps : denom_low;
        denom_upp = (abs(denom_upp) < small_eps) ? small_eps : denom_upp;
        denom_alpha = (abs(denom_alpha) < small_eps) ? small_eps : denom_alpha;
        denom_beta = (abs(denom_beta) < small_eps) ? small_eps : denom_beta;

        T sum = T(0);
        for (int i = 0; i < m; i++) {
            T lambda_i = lambda[i];
            sum += pij[i + tj * m] * lambda_i / (denom_upp * denom_upp)
                 - qij[i + tj * m] * lambda_i / (denom_low * denom_low);
        }
        sum += p0j[tj] / (denom_upp * denom_upp)
             - q0j[tj] / (denom_low * denom_low)
             - epsi / denom_alpha
             + epsi / denom_beta;

        delx[tj] = sum;
    }
}

template <typename T>
__global__ void GG_kernel(
    T* __restrict__ GG,
    const T* __restrict__ x,
    const T* __restrict__ xlow,
    const T* __restrict__ xupp,
    const T* __restrict__ pij,
    const T* __restrict__ qij,
    const int n,
    const int m) {

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    T xt = x[tj];
    T xlow_j = xlow[tj];
    T xupp_j = xupp[tj];

    // Distances from bounds
    T diff_upper = xupp_j - xt;
    T diff_lower = xt - xlow_j;

    // Squared distances
    T diff_upper2 = diff_upper * diff_upper;
    T diff_lower2 = diff_lower * diff_lower;

    for (int i = 0; i < m; i++) {
      int idx = i + tj * m;
      GG[idx] = pij[idx] / diff_upper2 - qij[idx] / diff_lower2;
    }
  }
}


template <typename T>
__global__ void diagx_kernel(T* __restrict__ diagx, const T* __restrict__ x,
     const T* __restrict__ xsi, const T* __restrict__ xlow,
     const T* __restrict__ xupp, const T* __restrict__ p0j,
   const T* __restrict__ q0j, const T* __restrict__ pij,
   const T* __restrict__ qij, const T* alpha, const T*  beta,
   const T*  eta, const T* lambda, const int n, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    T sum = 0;
    T sum1 = 0;
    for (int i = 0; i < m; i++) {
      sum = sum + pij[tj *m+ i] * lambda[i];
      sum1 = sum1 + qij[tj*m + i] * lambda[i];
    }
    diagx[tj] = (p0j[tj] + sum) / pow(xupp[tj] - x[tj], 3) +
       (q0j[tj] + sum1) / pow(x[tj] - xlow[tj], 3);
    diagx[tj] = 2.0 * diagx[tj] + xsi[tj] / (x[tj] - alpha[tj]) +
       eta[tj] / (beta[tj] - x[tj]);
  }
}


template <typename T>
__inline__ __device__ T reduce_warp(T val) {
  volatile T v = val;
  v += __shfl_down_sync(0xffffffff, v, 16);
  v += __shfl_down_sync(0xffffffff, v, 8);
  v += __shfl_down_sync(0xffffffff, v, 4);
  v += __shfl_down_sync(0xffffffff, v, 2);
  v += __shfl_down_sync(0xffffffff, v, 1);
  return v;
}


template <typename T>
__global__ void mmareduce_kernel(T* __restrict__ bufred, const int n) {

  T sum = 0;
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;
  for (int i = idx; i < n; i += str)
  {
    sum += bufred[i];
  }

  __shared__ T shared[32];
  unsigned int lane = threadIdx.x % warpSize;
  unsigned int wid = threadIdx.x / warpSize;

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    bufred[blockIdx.x] = sum;
}


template< typename T >
__global__ void mmasum_kernel(const T*  __restrict__ a, T*  __restrict__ buf_h,
     const int n, const int m, const int k) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T shared[32];
  T sum = 0;
  for (int i = idx; i < n; i += str)
  {
    sum += a[m * i + k ];
  }

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum;

}
template< typename T >
__global__ void mmasumbb_kernel(const T*  __restrict__ GG,
     const T*  __restrict__ delx, const T*  __restrict__ diagx,
     T*  __restrict__ buf_h, const int n, const int m, const int k) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T shared[32];
  T sum = 0;
  for (int i = idx; i < n; i += str)
  {
    sum += GG[ k + i * m] * delx[i] / diagx[i];
  }

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum;

}

template <typename T>
__global__ void mmasumHess_kernel(
    const T* __restrict__ hijx,
    const T* __restrict__ Ljjxinv,
    T* __restrict__ buf_h,
    const int n,
    const int m,
    const int k0,
    const int k1)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int stride = blockDim.x * gridDim.x;

    // Warp lane and warp id within the block
    const unsigned int lane = threadIdx.x % warpSize;
    const unsigned int wid = threadIdx.x / warpSize;

    __shared__ T shared[32];

    T sum = T(0);

    // Grid-stride loop for global reduction
    for (int i = idx; i < n; i += stride) {
        // hijx is indexed as hijx[offset + row * m]
        T val0 = hijx[k0 + i * m];
        T val1 = hijx[k1 + i * m];
        sum += val0 * Ljjxinv[i] * val1;
    }

    // Warp-level reduction using your reduce_warp implementation
    sum = reduce_warp<T>(sum);

    // Write each warp's partial sum to shared memory
    if (lane == 0) {
        shared[wid] = sum;
    }
    __syncthreads();

    // First warp reduces the warp sums in shared memory
    sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : T(0);
    if (wid == 0) {
        sum = reduce_warp<T>(sum);
    }

    // Write final result from first thread of block
    if (threadIdx.x == 0) {
        buf_h[blockIdx.x] = sum;
    }
}


template< typename T >
__global__ void mmasumAA_kernel(const T*  __restrict__ GG,
     const T*  __restrict__ diagx, T*  __restrict__ buf_h, const int n,
   const int m, const int k0, const int k1) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T shared[32];
  T sum = 0;
  for (int i = idx; i < n; i += str)
  {
    sum += GG[ k0 + i * m] /diagx[i]  * GG[ k1 + i * m];
  }

  sum = reduce_warp<T>(sum);
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum;

}


template <typename T>
__global__ void mma_copy_kernel(T* __restrict__ a, const T* __restrict__ b,
     const int n, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
    if(tj<n)
      a[tj+m]=b[tj];
}



template <typename T>
__global__ void AA_kernel(T* __restrict__ temp, const T* __restrict__ GG,
     const T* __restrict__ diagx, const int n, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    for (int i0 = 0; i0 < m; i0++) {
      for (int i1 = 0; i1 < m; i1++) {
        temp[tj + i0 * (n + 1) + i1 * (m + 1) * (n + 1)] = GG[i0 * n + tj] *
         (1.0 / diagx[tj]) * GG[i1 * n + tj];
      }
    }
  }
}


template <typename T>
__global__ void dx_kernel(T* __restrict__ dx, const T* __restrict__ delx,
     const T* __restrict__ diagx, const T* __restrict__ GG,
     const T* __restrict__ dlambda, const int n, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    dx[tj] = -delx[tj]/diagx[tj];
    for(int i=0;i<m;i++){
      dx[tj] =dx[tj] - GG[tj*m+i]*dlambda[i]/diagx[tj];
    }
  }
}



template <typename T>
__global__ void dxsi_kernel(T* __restrict__ dxsi, const T* __restrict__ xsi,
     const T* __restrict__ dx, const T* __restrict__ x,
     const T* __restrict__ alpha, const T epsi, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    dxsi[tj]= -xsi[tj] + (epsi-dx[tj]*xsi[tj])/(x[tj] - alpha[tj]);
  }
}


template <typename T>
__global__ void deta_kernel(T* __restrict__ deta, const T* __restrict__ eta,
     const T* __restrict__ dx, const T* __restrict__ x,
     const T* __restrict__ beta, const T epsi, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    deta[tj] = -eta[tj] + (epsi + dx[tj] * eta[tj]) / (beta[tj] - x[tj]);
  }
}


template <typename T>
__global__ void RexCalculation_kernel(
    T* __restrict__ rex,
    const T* __restrict__ x,
    const T* __restrict__ xlow,
    const T* __restrict__ xupp,
    const T* __restrict__ pij,
    const T* __restrict__ p0j,
    const T* __restrict__ qij,
    const T* __restrict__ q0j,
    const T* __restrict__ lambda,
    const T* __restrict__ xsi,
    const T* __restrict__ eta,
    const int n,
    const int m)
{
    int tj = blockIdx.x * blockDim.x + threadIdx.x;
    if (tj < n) {
        T upp_diff = xupp[tj] - x[tj];
        T low_diff = x[tj] - xlow[tj];
        T upp_diff_sq = upp_diff * upp_diff;
        T low_diff_sq = low_diff * low_diff;

        T sum = 0.0;
        for (int i = 0; i < m; ++i) {
            sum += pij[i + tj * m] * lambda[i] / upp_diff_sq
                 - qij[i + tj * m] * lambda[i] / low_diff_sq;
        }

        rex[tj] = sum
                 + p0j[tj] / upp_diff_sq
                 - q0j[tj] / low_diff_sq
                 - xsi[tj] + eta[tj];
    }
}


template <typename T>
__global__ void rey_calculation_kernel(T* __restrict__ rey,
     const T* __restrict__ c, const T* __restrict__ d, const T* __restrict__ y,
     const T* __restrict__ lambda, const T* __restrict__ mu, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    rey[tj] = c[tj] + d[tj] * y[tj] - lambda[tj] - mu[tj];
  }
}


template <typename T>
__global__ void norm_kernel(const T* __restrict__ a,
      T* __restrict__ buf_h, const int n) {
  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T shared[32];
  T sum = T(0);

  for (int i = idx; i < n; i += str) {
    T val = a[i];
    sum += val * val;  // faster than pow(val, 2)
  }

  sum = reduce_warp<T>(sum);  // your warp-level sum reduction function
  if (lane == 0)
    shared[wid] = sum;
  __syncthreads();

  sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : T(0);
  if (wid == 0)
    sum = reduce_warp<T>(sum);

  if (threadIdx.x == 0)
    buf_h[blockIdx.x] = sum;
}


template <typename T>
__global__ void sub2cons_kernel(T* __restrict__ a, const T* __restrict__ b,
     const T* __restrict__ c,
  const T d, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    a[tj] = b[tj]*c[tj]-d;
  }
}


template <typename T>
__global__ void dely_kernel(T* __restrict__ dely, const T* __restrict__ c,
     const T* __restrict__ d, const T* __restrict__ y,
     const T* __restrict__ lambda, const T epsi, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    dely[tj] = c[tj] + d[tj]*y[tj] - lambda[tj] - epsi/y[tj];
  }
}



template< typename T >
__global__ void maxval2_kernel(const T* __restrict__ a, const T* __restrict__ b,
     T* __restrict__ temp, const T cons, const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T shared[32];
  T maxval = cons * a[0] / b[0];
  for (int i = idx; i < n; i += str)
  {
    maxval = max(maxval, cons * a[i] / b[i]);
  }

  maxval = max_reduce_warp<T>(maxval);
  if (lane == 0)
    shared[wid] = maxval;
  __syncthreads();

  maxval = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0.0;
  if (wid == 0)
    maxval = max_reduce_warp<T>(maxval);

  if (threadIdx.x == 0)
    temp[blockIdx.x] = maxval;

}


template< typename T >
__global__ void maxval3_kernel(const T* __restrict__ a, const T* __restrict__ b,
     const T* __restrict__ c, T* __restrict__ temp, const T cons, const int n) {

  const int idx = blockIdx.x * blockDim.x + threadIdx.x;
  const int str = blockDim.x * gridDim.x;

  const unsigned int lane = threadIdx.x % warpSize;
  const unsigned int wid = threadIdx.x / warpSize;

  __shared__ T shared[32];
  T maxval = cons * a[0] / b[0];
  for (int i = idx; i < n; i += str)
  {
    maxval = max(maxval, cons * a[i] / (b[i] - c[i]));
  }

  maxval = max_reduce_warp<T>(maxval);
  if (lane == 0)
    shared[wid] = maxval;
  __syncthreads();

  maxval = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
  if (wid == 0)
    maxval = max_reduce_warp<T>(maxval);

  if (threadIdx.x == 0)
    temp[blockIdx.x] = maxval;

}


template <typename T>
__global__ void kkt_rex_kernel(T* __restrict__ rex, const T* __restrict__ df0dx,
     const T* __restrict__ dfdx, const T* __restrict__ xsi,
     const T* __restrict__ eta, const T* __restrict__ lambda, const int n,
   const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    rex[tj] = 0.0;
    for (int i = 0; i < m; i++) {
      rex[tj] = rex[tj] + dfdx[i + tj*m] * lambda[i];
    }
    rex[tj] += df0dx[tj] - xsi[tj] + eta[tj];
  }
}


template <typename T>
__global__ void maxcons_kernel(T* __restrict__ a, const T b,
     const T c, const T* __restrict__ d, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    a[tj] = max(b, c * d[tj]);
  }
}

///////////////////////////////////
//////defined
 template< typename T >
 __global__ void glsum_kernel(const T * a, T * buf_h, const int n) {
   const int idx = blockIdx.x * blockDim.x + threadIdx.x;
   const int str = blockDim.x * gridDim.x;

   const unsigned int lane = threadIdx.x % warpSize;
   const unsigned int wid = threadIdx.x / warpSize;

   __shared__ T shared[32];
   T sum = 0;
   for (int i = idx; i<n ; i += str)
   {
     sum += a[i];
   }

   sum = reduce_warp<T>(sum);
   if (lane == 0)
     shared[wid] = sum;
   __syncthreads();

   sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
   if (wid == 0)
     sum = reduce_warp<T>(sum);

   if (threadIdx.x == 0)
     buf_h[blockIdx.x] = sum;

 }
  template< typename T >
__global__ void glsc2_kernel(const T * a,
                              const T * b,
                              T * buf_h,
                              const int n) {

   const int idx = blockIdx.x * blockDim.x + threadIdx.x;
   const int str = blockDim.x * gridDim.x;

   const unsigned int lane = threadIdx.x % warpSize;
   const unsigned int wid = threadIdx.x / warpSize;

   __shared__ T shared[32];
   T sum = 0.0;
   for (int i = idx; i < n; i+= str) {
     sum += a[i] * b[i];
   }

   sum = reduce_warp<T>(sum);
   if (lane == 0)
     shared[wid] = sum;
   __syncthreads();

   sum = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
   if (wid == 0)
     sum = reduce_warp<T>(sum);

   if (threadIdx.x == 0)
     buf_h[blockIdx.x] = sum;

 }





template <typename T>
__global__ void add2inv2_kernel(T* __restrict__ a, const T* __restrict__ b,
     const T c, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    a[tj] = a[tj]+c/b[tj];
  }
}

template <typename T>
__global__ void max2_kernel(T* __restrict__ a, const T b,
     const T* __restrict__ c, const T d, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if (tj < n) {
    a[tj]=max(b, d*c[tj]);
  }
}

template <typename T>
__global__ void updatebb_kernel(T* __restrict__ bb,
     const T* __restrict__ dellambda, const T* __restrict__ dely,
   const T* __restrict__ d, const T* __restrict__ mu,
   const T* __restrict__ y, const T delz, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if(tj<m)
    bb[tj]=dellambda[tj] + dely[tj]/(d[tj] + mu[tj]/y[tj]) - bb[tj];
  else if(tj<m+1)
    bb[tj]=delz;
}



template <typename T>
__global__ void updateAA_kernel(T* __restrict__ AA,
     const T* __restrict__ globaltmp_mm, const T* __restrict__ s,
   const T* __restrict__ lambda, const T* __restrict__ d,
     const T* __restrict__ mu, const T* __restrict__ y, const T* __restrict__ a,
   const T zeta, const T z, const int m) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if(tj<m)
    {
      AA[tj+tj*(m+1)]=globaltmp_mm[tj+tj*m] + (s[tj] / lambda[tj] +
         1.0/ (d[tj] + mu[tj] / y[tj]));
      AA[tj+m*(m+1)]=a[tj];
      AA[m+tj*(m+1)]=a[tj];
    }
  else if(tj<m+1)
    AA[tj+tj*(m+1)]= -zeta/z;
}

template <typename T>
__global__ void dy_kernel(T* __restrict__ dy, const T* __restrict__ dely,
     const T* __restrict__ dlambda, const T* __restrict__ d,
     const T* __restrict__ mu, const T* __restrict__ y, const int n) {
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  if(tj<n)
    dy[tj] = (-dely[tj]+dlambda[tj])/(d[tj] + mu[tj]/y[tj]);
}

#endif

