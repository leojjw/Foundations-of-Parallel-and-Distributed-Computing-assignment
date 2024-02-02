#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include "SpMV.hpp"
#include "VectorOperator.hpp"

const int threads_per_block = 256;

__global__
void gpu_SpMV_kernel(double* d_A, double* d_V, double* d_result, unsigned int M, unsigned int Xnode) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i == 0 || i >= M) return;
  d_result[i] = d_A[i] * d_V[i] + d_A[M + i] * d_V[i - 1 + Xnode] +
                d_A[2 * M + i] * d_V[i + Xnode] + d_A[3 * M + i] * d_V[i + 1 + Xnode] +
                d_A[4 * M + i] * d_V[i + Xnode + Xnode];
}

void gpu_SpMV(double** A, double* V, double* result, unsigned int M, unsigned int Xnode) {
  int size = M * sizeof(double);
  double *d_A, *d_V, *d_result;

  cudaMalloc((void **) &d_A, 5 * size);
  cudaMalloc((void **) &d_V, size + 2 * Xnode * sizeof(double));
  cudaMalloc((void **) &d_result, size);
  for (int i = 0; i < 5; i++){
    cudaMemcpy(d_A + M * i, A[i], size, cudaMemcpyHostToDevice);
  }
  cudaMemcpy(d_V, V, size + 2 * Xnode * sizeof(double), cudaMemcpyHostToDevice);

  dim3 grid_dim(ceil((double)M/threads_per_block), 1, 1);
  dim3 block_dim(threads_per_block, 1, 1);
  gpu_SpMV_kernel<<<grid_dim, block_dim>>>(d_A, d_V, d_result, M, Xnode);

  cudaMemcpy(result, d_result, size, cudaMemcpyDeviceToHost);
  cudaFree(d_A); cudaFree(d_V); cudaFree(d_result);
}

__global__
void gpu_VectorDotVector_kernel(double* d_a, double* d_b, double* d_sum, unsigned int L) {
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if (i >= L) return;
  d_sum[i] = d_a[i] * d_b[i];
}

double gpu_VectorDotVector(double* a, double* b, unsigned int L) {
  int size = L * sizeof(double);
  double *d_a, *d_b, *d_sum, *sum;
  sum = new double[L];

  cudaMalloc((void **) &d_a, size);
  cudaMalloc((void **) &d_b, size);
  cudaMalloc((void **) &d_sum, size);
  cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);

  dim3 grid_dim(ceil((double)L/threads_per_block), 1, 1);
  dim3 block_dim(threads_per_block, 1, 1);
  gpu_VectorDotVector_kernel<<<grid_dim, block_dim>>>(d_a, d_b, d_sum, L);
  cudaDeviceSynchronize();

  cudaMemcpy(sum, d_sum, size, cudaMemcpyDeviceToHost);

  double result = 0;
  for (int i = 0; i < L; i++) {
    result += sum[i];
  }
  cudaFree(d_a); cudaFree(d_b); cudaFree(d_sum);
  delete[] sum;

  return result;
}

double gpu_VectorNorm(double* a, unsigned int L) {
  return gpu_VectorDotVector(a, a, L);
}

__device__
void lock(int* mutex) {
  while (atomicCAS(mutex, 0, 1) != 0);
}

__device__
void unlock(int* mutex) {
  atomicExch(mutex, 0);
}

__global__
void gpu_VectorDotVector_optimized_kernel(double* d_a, double* d_b, double* d_sum, unsigned int L, int* mutex) {
  __shared__ double thread_sum[threads_per_block];
  
  int idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx >= L) return;
  thread_sum[threadIdx.x] = d_a[idx] * d_b[idx];
  __syncthreads();

  int i = threads_per_block / 2;
  while(i != 0) {
    if (threadIdx.x < i && idx + i < L) {
      thread_sum[threadIdx.x] += thread_sum[threadIdx.x + i];
    }
    __syncthreads();
    i = i / 2;
  }

  if (threadIdx.x == 0) {
    lock(mutex);
    *d_sum += thread_sum[0];
    unlock(mutex);
  }
}

double gpu_VectorDotVector_optimized(double* a, double* b, unsigned int L) {
  int size = L * sizeof(double);
  double sum = 0;
  double *d_a, *d_b, *d_sum;

  cudaMalloc((void **) &d_a, size);
  cudaMalloc((void **) &d_b, size);
  cudaMalloc((void **) &d_sum, sizeof(double));
  cudaMemcpy(d_a, a, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b, size, cudaMemcpyHostToDevice);
  cudaMemcpy(d_sum, &sum, sizeof(double), cudaMemcpyHostToDevice);
  
  int *mutex = NULL;
  cudaMalloc((void **)&mutex, sizeof(int));
  cudaMemset(mutex, 0, sizeof(int));

  dim3 grid_dim(ceil((double)L/threads_per_block), 1, 1);
  dim3 block_dim(threads_per_block, 1, 1);
  gpu_VectorDotVector_optimized_kernel<<<grid_dim, block_dim>>>(d_a, d_b, d_sum, L, mutex);
  cudaDeviceSynchronize();
  cudaMemcpy(&sum, d_sum, sizeof(double), cudaMemcpyDeviceToHost);
  cudaFree(d_a); cudaFree(d_b); cudaFree(d_sum); cudaFree(mutex);

  return sum;
}

double gpu_VectorNorm_optimized(double* a, unsigned int L) {
  return gpu_VectorDotVector(a, a, L);
}