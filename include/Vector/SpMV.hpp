#ifndef _SPMV_
#define _SPMV_

void gpu_SpMV<double>(double** A, double* V, double* result, unsigned int M, unsigned int Xnode);

template <class T>
void SpMV(T** A, T* V, T* result, unsigned int M, unsigned int Xnode) {
    for (int i = 1; i < M; i++) {
        result[i] = A[0][i] * V[i] + A[1][i] * V[i - 1 + Xnode] +
                    A[2][i] * V[i + Xnode] + A[3][i] * V[i + 1 + Xnode] +
                    A[4][i] * V[i + Xnode + Xnode];
    }
}
#endif