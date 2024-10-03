#ifndef MMA_GPU_KERNEL_H   
#define MMA_GPU_KERNEL_H   
#endif


template< typename T >
__global__ void BoundCalculation_kernel(const T* __restrict__ x, const T* __restrict__ xold1, const T* __restrict__ xold2, const T* __restrict__ df0dx, const T* __restrict__ dfdx,
    T* __restrict__ xlow, T* __restrict__ xupp, const T* __restrict__ xmin, const T* __restrict__ xmax,
    T* __restrict__ alpha, T* __restrict__ beta, T* __restrict__ p0j, T* __restrict__ q0j, T* __restrict__ pij, T* __restrict__ qij, T* __restrict__ temp,
    const T asyinit, const T asydecr, const T asyincr,
    const int n, const int m, const int iter) {
    int tj = blockIdx.x * blockDim.x + threadIdx.x;
    if (tj < n) {
        T xgap = xmax[tj] - xmin[tj];
        if (iter < 3) {
            xlow[tj] = x[tj] - asyinit * xgap;
            xupp[tj] = x[tj] + asyinit * xgap;
        }

        else {
            T xdiff = (x[tj] - xold1[tj]) * (xold1[tj] - xold2[tj]);
            if (xdiff < 0)
            {
                xlow[tj] = x[tj] - asydecr * (xold1[tj] - xlow[tj]);
                xupp[tj] = x[tj] + asydecr * (xupp[tj] - xold1[tj]);
            }
            else if (xdiff > 0)
            {
                xlow[tj] = x[tj] - asyincr * (xold1[tj] - xlow[tj]);
                xupp[tj] = x[tj] + asyincr * (xupp[tj] - xold1[tj]);
            }
            else {
                xlow[tj] = x[tj] - (xold1[tj] - xlow[tj]);
                xupp[tj] = x[tj] + (xupp[tj] - xold1[tj]);
            }
            xlow[tj] = max(xlow[tj], x[tj] - 10 * xgap);
            xlow[tj] = min(xlow[tj], x[tj] - 0.01 * xgap);
            xupp[tj] = min(xupp[tj], x[tj] + 10 * xgap);
            xupp[tj] = max(xupp[tj], x[tj] - 0.01 * xgap);

        }
        alpha[tj] = max(max(xmin[tj], xlow[tj] + 0.1 * (x[tj] - xlow[tj])), x[tj] - 0.5 * xgap);
        beta[tj] = min(min(xmax[tj], xupp[tj] - 0.1 * (xupp[tj] - x[tj])), x[tj] + 0.5 * xgap);

        //Calculate p0j, q0j, pij, qij
       ///where j = 1, 2, ..., n and i = 1, 2, ..., m(eq(2.3) - eq(2.5))
        p0j[tj] = pow(xupp[tj] - x[tj], 2) * (1.001 * max(df0dx[tj], 0.0) + 0.001 * max(-df0dx[tj], 0.0) + 0.00001 / max(0.00001, xgap));

        q0j[tj] = pow(x[tj] - xlow[tj], 2) * (0.001 * max(df0dx[tj], 0.0) + 1.001 * max(-df0dx[tj], 0.0) + 0.00001 / max(0.00001, xgap));
        for (int i = 0; i < m; i++) {
            pij[tj + i * n] = pow(xupp[tj] - x[tj], 2) * (1.001 * max(dfdx[tj + i * n], 0.0) + 0.001 * max(-dfdx[tj + i * n], 0.0) + 0.00001 / max(0.00001, xgap));
            qij[tj + i * n] = pow(x[tj] - xlow[tj], 2) * (0.001 * max(dfdx[tj + i * n], 0.0) + 1.001 * max(-dfdx[tj + i * n], 0.0) + 0.00001 / max(0.00001, xgap));
            temp[tj + i * n] = pij[tj + i * n] / (xupp[tj] - x[tj]) + qij[tj + i * n] / (x[tj] - xlow[tj]);
        }
    }
}

template <typename T>
__device__ void warpReduce(volatile T* sdata, int tid) {
    sdata[tid] += sdata[tid + 32];
    sdata[tid] += sdata[tid + 16];
    sdata[tid] += sdata[tid + 8];
    sdata[tid] += sdata[tid + 4];
    sdata[tid] += sdata[tid + 2];
    sdata[tid] += sdata[tid + 1];
}
template <typename T>
__global__ void sum_reductions(const T* __restrict__ input, T* __restrict__ output, int size, int n) {
    __shared__ T sdata[1024];
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
    int gridSize = 1024 * 2 * gridDim.x;
    sdata[tid] = 0;
    while (i < size) {
        sdata[tid] += input[i + size * n] + ((i + blockDim.x < size) ? input[i + size * n + blockDim.x] : 0);
        i += gridSize;
    }
    __syncthreads();
    if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads();
    if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads();
    if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads();
    if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads();
    if (tid < 32) warpReduce(sdata, tid);
    if (tid == 0) output[blockIdx.x] = sdata[0];
}



template <typename T>
__global__ void sum_reductions_final(T* __restrict__ input, int size) {
    __shared__ T sdata[1024];
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
    int gridSize = 1024 * 2 * gridDim.x;
    sdata[tid] = 0;

    while (i < size) {
        sdata[tid] = input[i] + input[i + blockDim.x];
        i += gridSize;
    }
    __syncthreads();
    if (tid < 512) { sdata[tid] += sdata[tid + 512]; } __syncthreads();
    if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads();
    if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads();
    if (tid < 64) { sdata[tid] += sdata[tid + 64]; } __syncthreads();
    if (tid < 32) warpReduce(sdata, tid);
    if (tid == 0) input[blockIdx.x] = sdata[0];
}



template< typename T>
__inline__ __device__ T reduce_warp(T val) {
	val += __shfl_down_sync(0xffffffff, val, 16);
	val += __shfl_down_sync(0xffffffff, val, 8);
	val += __shfl_down_sync(0xffffffff, val, 4);
	val += __shfl_down_sync(0xffffffff, val, 2);
	val += __shfl_down_sync(0xffffffff, val, 1);
	return val;
}

template <typename T>
__global__ void mmareduce_kernel(T* bufred, const int n) {

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
__global__ void mmasum_kernel(const T* a, T* buf_h, const int n, const int k) {

    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int str = blockDim.x * gridDim.x;

    const unsigned int lane = threadIdx.x % warpSize;
    const unsigned int wid = threadIdx.x / warpSize;

    __shared__ T shared[32];
    T sum = 0;
    for (int i = idx; i < n; i += str)
    {
        sum += a[i + k * n];
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
__global__ void mma_max_kernel(const T* __restrict__ x, T* __restrict__ alpha, T* __restrict__ beta,
    T* __restrict__ xsi, T* __restrict__ eta, T* __restrict__ mu, T* __restrict__ c, const int n) {
    int tj = blockIdx.x * blockDim.x + threadIdx.x;
    if (tj < n) {
        xsi[tj] = max(1.0, 1.0 / (x[tj] - alpha[tj]));
        eta[tj] = max(1.0, 1.0 / (beta[tj] - x[tj]));
        mu[tj] = max(1.0, 0.5 * c[tj]);
    }
}

template <typename T>
__global__ void RexCalculation_kernel(T* __restrict__ rex, const T* __restrict__ x, const T* __restrict__ xlow, const T* __restrict__ xupp, const T* __restrict__ pij, const T* __restrict__ p0j,
    const T* __restrict__ qij, const T* __restrict__ q0j, const T* __restrict__ lambda, const T* __restrict__ xsi, const T* __restrict__ eta, const int n, const int m) {
    int tj = blockIdx.x * blockDim.x + threadIdx.x;
    if (tj < n) {
        rex[tj] = 0.0;
        for (int i = 0; i < m; i++) {
            rex[tj] = rex[tj] + pij[tj + i * n] * lambda[i] / pow(xupp[tj] - x[tj], 2) - qij[tj + i * n] * lambda[i] / pow(x[tj] - xlow[tj], 2);
        }
        rex[tj] = rex[tj] + p0j[tj] / pow(xupp[tj] - x[tj], 2) - q0j[tj] / pow(x[tj] - xlow[tj], 2) - xsi[tj] + eta[tj];
    }
}

template <typename T>
__global__ void rey_calculation_kernel(T* __restrict__ rey, const T* __restrict__ c, const T* __restrict__ d, const T* __restrict__ y,
	const T* __restrict__ lambda, const T* __restrict__ mu, const int n) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		rey[tj] = c[tj] + d[tj] * y[tj] - lambda[tj] - mu[tj];
	}
}

template <typename T>
__global__ void dot_product_kernel(T* __restrict__ temp, const T* __restrict__ lambda, const T* __restrict__ mu,
	const T* __restrict__ a, const int n) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		temp[tj] = lambda[tj] * a[tj];
	}
}

template <typename T>
__global__ void relambda_kernel(T* __restrict__ temp, const T* __restrict__ x, const T* __restrict__ xupp, const T* __restrict__ xlow,
	const T* __restrict__ pij, const T* __restrict__ qij, const int n, const int m) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		for (int i = 0; i < m; i++) {
			temp[tj + i * n] = pij[tj + i * n] / (xupp[tj] - x[tj]) + qij[tj + i * n] / (x[tj] - xlow[tj]);
		}
	}
}


template <typename T>
__global__ void residu_kernel(T* __restrict__ rexsi,  T* __restrict__ reeta,  T* __restrict__ remu,  T* __restrict__ res,
	const T* __restrict__ x, const T* __restrict__ y, const T* __restrict__ xsi, const T* __restrict__ alpha, const T* __restrict__ beta,
	const T* __restrict__ eta, const T* __restrict__ mu, const T* __restrict__ lambda, const T* __restrict__ s,
	const T epsi, const int n) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		rexsi[tj] = xsi[tj] * (x[tj] - alpha[tj]) - epsi;
		reeta[tj] = eta[tj] * (beta[tj] - x[tj]) - epsi;
		remu[tj] = mu[tj] * y[tj] - epsi;
		res[tj] = lambda[tj] * s[tj] - epsi;
	}
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



template< typename T >
__global__ void maxval_kernel(const T* a, T *temp, const int n) {

	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const int str = blockDim.x * gridDim.x;

	const unsigned int lane = threadIdx.x % warpSize;
	const unsigned int wid = threadIdx.x / warpSize;

	__shared__ T shared[32];
	T maxval = 0;
	for (int i = idx; i < n; i += str)
	{
		maxval = max(maxval,a[i]);
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
__global__ void max_reduce_kernel(T* bufred, const int n) {

	T maxval = 0;
	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const int str = blockDim.x * gridDim.x;
	for (int i = idx; i < n; i += str)
	{
		maxval =max(maxval, bufred[i]);
	}

	__shared__ T shared[32];
	unsigned int lane = threadIdx.x % warpSize;
	unsigned int wid = threadIdx.x / warpSize;

	maxval = max_reduce_warp<T>(maxval);
	if (lane == 0)
		shared[wid] = maxval;
	__syncthreads();

	maxval = (threadIdx.x < blockDim.x / warpSize) ? shared[lane] : 0;
	if (wid == 0)
		maxval = max_reduce_warp<T>(maxval);

	if (threadIdx.x == 0)
		bufred[blockIdx.x] = maxval;
}



template< typename T >
__global__ void norm_kernel(const T* a, T* buf_h, const int n) {

	const int idx = blockIdx.x * blockDim.x + threadIdx.x;
	const int str = blockDim.x * gridDim.x;

	const unsigned int lane = threadIdx.x % warpSize;
	const unsigned int wid = threadIdx.x / warpSize;

	__shared__ T shared[32];
	T sum = 0;
	for (int i = idx; i < n; i += str)
	{
		sum += pow(a[i], 2);
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
__global__ void delx_kernel(T* __restrict__ delx, const T* __restrict__ x, const T* __restrict__ xlow, const T* __restrict__ xupp,
	const T* __restrict__ pij, const T* __restrict__ qij, const T* __restrict__ p0j, const T* __restrict__ q0j,
	const T* __restrict__ alpha, const T* __restrict__ beta, const T* __restrict__ lambda, const T epsi, const int n, const int m) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		for (int i = 0; i < m; i++) {
			delx[tj] = delx[tj] + pij[tj + i * n] * lambda[i] / pow(xupp[tj] - x[tj], 2) -
				qij[tj + i * n] * lambda[i] / pow(x[tj] - xlow[tj], 2);
		}
		delx[tj] = delx[tj] + p0j[tj] / pow(xupp[tj] - x[tj], 2) - q0j[tj] / pow(x[tj] - xlow[tj], 2) - epsi / (x[tj] - alpha[tj])
			+ epsi / (beta[tj] - x[tj]);
	}
}

template <typename T>
__global__ void dellambda_kernel(T* __restrict__ temp, const T* __restrict__ x, const T* __restrict__ xlow, const T* __restrict__ xupp,
	const T* __restrict__ pij, const T* __restrict__ qij, const int n, const int m) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		temp[tj] = pij[tj] / (xupp[tj] - x[tj]) + qij[tj] / (x[tj] - xlow[tj]);
	}
}

template <typename T>
__global__ void GG_kernel(T* __restrict__ GG, const T* __restrict__ x, const T* __restrict__ xlow, const T* __restrict__ xupp,
	const T* __restrict__ pij, const T* __restrict__ qij, const int n, const int m) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		for (int ggdumiter = 0; ggdumiter < m; ggdumiter++) {
			GG[ggdumiter * n + tj] = pij[ggdumiter * n + tj] / pow(xupp[tj] - x[tj], 2) -
				qij[ggdumiter * n + tj] / pow(x[tj] - xlow[tj], 2);
		}
	}
}

template <typename T>
__global__ void diagx_kernel(T* __restrict__ diagx, const T* __restrict__ x, const T* __restrict__ xsi,const T* __restrict__ xlow, const T* __restrict__ xupp,
	const T* __restrict__ p0j, const T* __restrict__ q0j, const T* __restrict__ pij, const T* __restrict__ qij, 
	const T* alpha, const T*  beta, const T*  eta, const T* lambda, const int n, const int m) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		T sum = 0;
		T sum1 = 0;
		for (int i = 0; i < m; i++) {
			sum = sum + pij[tj + i * n] * lambda[i];
			sum1 = sum1 + qij[tj + i * n] * lambda[i];
		}
		diagx[tj] = (p0j[tj] + sum) / pow(xupp[tj] - x[tj], 3) + (q0j[tj] + sum1) / pow(x[tj] - xlow[tj], 3);
		diagx[tj] = 2 * diagx[tj] + xsi[tj] / (x[tj] - alpha[tj]) + eta[tj] / (beta[tj] - x[tj]);
	}
}


template <typename T>
__global__ void bb_kernel(T* __restrict__ temp, const T* __restrict__ GG, const T* __restrict__ delx, const T* __restrict__ diagx,
	const int n) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		temp[tj] = GG[tj] * (delx[tj] / diagx[tj]);
		temp[tj + n] = GG[tj + n] * (delx[tj] / diagx[tj]);
	}
}

template <typename T>
__global__ void AA_kernel(T* __restrict__ temp, const T* __restrict__ GG, const T* __restrict__ diagx,
	const int n, const int m) {
	int tj = blockIdx.x * blockDim.x + threadIdx.x;
	if (tj < n) {
		for (int i0 = 0; i0 < m; i0++) {
			for (int i1 = 0; i1 < m; i1++) {
				temp[tj + i0 * n + i1 * m * n] = GG[i0 * n + tj] * (1.0 / diagx[tj]) * GG[i1 * n + tj];
			}
		}
	}
}
