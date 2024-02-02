#ifndef _GEMM_
#define _GEMM_
template <class T>
void Gemm(T** A, T* V, T* result, unsigned int M, unsigned int N) {
    for (int i = 0; i < M; i++) {
        T sum = 0;
        for (int j = 0; j < N; j++) {
            sum += A[i][j] * V[j];
        }
        result[i] = sum;
    }
}
#endif