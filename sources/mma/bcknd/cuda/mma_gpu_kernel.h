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
__global__ void reduce_kernel(T* bufred, const int n) {

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
