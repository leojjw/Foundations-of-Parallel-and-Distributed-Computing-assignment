
#include <ArrayAllocate/ArrayAllocate.hpp>
#include <Vector/SpMV.hpp>
#include <Vector/VectorOperator.hpp>
#include <cmath>
#include <iostream>
#include <omp.h>
void CG(double** A,
        double* x,
        double* b,
        unsigned int Xnode,
        unsigned int Ynode) {
    std::cout << "Begin to calculate Ax=b using CG method" << std::endl;

    unsigned int LEN = Xnode * Ynode;

    /*declar paramters*/
    double* r = Allocate1D<double>(LEN);
    double* r_tmp = Allocate1D<double>(LEN);
    double* tmp = Allocate1D<double>(LEN);
    double* p_all = Allocate1D<double>(LEN + 2 * Xnode);
    double* p = p_all + Xnode;

    for (int i = 0; i < LEN; i++) {
        r[i] = x[i] = 0;
    }

    /*initializition*/
    /*tmp=A*x0*/
    gpu_SpMV<double>(A, x - Xnode, tmp, LEN, Xnode);

    int nthreads, tlen, j;

    #pragma omp parallel default(shared)
    {
        nthreads = omp_get_num_threads();
        tlen = LEN / nthreads;

        /*initializition*/
        #pragma omp for
            for (j = 0; j < LEN ; j += tlen) {
                /*r0=b-tmp*/
                VectorOperation([](int i, auto a, auto b) { return a[i] - b[i]; }, 
                                (j + tlen) < LEN ? tlen : (LEN - j), r + j, b + j, tmp + j);
                /*p0=r0*/
                VectorOperation([](int i, auto a) { return a[i]; },
                                (j + tlen) < LEN ? tlen : (LEN - j), p + j, r + j);
            }

    }

    int iteration = 0;
    while (true) {
        iteration++;

        auto tmp_v = gpu_VectorNorm_optimized(r, LEN);

        gpu_SpMV(A, p - Xnode, tmp, LEN, Xnode);

        auto alpha = tmp_v / gpu_VectorDotVector_optimized(p, tmp, LEN);

        #pragma omp parallel for default(shared)
            for (j = 0; j < LEN ; j += tlen) {
                /* x_k+1=x_k+alpha*p */
                VectorOperation(
                    [alpha](int i, auto a, auto b) { return a[i] + alpha * b[i]; }, (j + tlen) < LEN ? tlen : (LEN - j),
                    x + j, x + j, p + j);

                /* r_tmp=r-alpha*tmp */
                VectorOperation(
                    [alpha](int i, auto a, auto b) { return a[i] - alpha * b[i]; }, (j + tlen) < LEN ? tlen : (LEN - j),
                    r_tmp + j, r + j, tmp + j);
            }

        auto beta = gpu_VectorNorm_optimized(r_tmp, LEN) / gpu_VectorNorm_optimized(r, LEN);

        #pragma omp parallel for default(shared)
            for (j = 0; j < LEN ; j += tlen) {
                /* p=r_tmp+beta*p */
                VectorOperation(
                    [beta](int i, auto a, auto b) { return a[i] + beta * b[i]; }, (j + tlen) < LEN ? tlen : (LEN - j),
                    p + j, r_tmp + j, p + j);

                /* r=r_tmp */
                VectorOperation([](int i, auto a) { return a[i]; }, (j + tlen) < LEN ? tlen : (LEN - j), r + j, r_tmp + j);
            }

        /*the loop exit condition*/
        if (iteration % 10 == 0)
            std::cout << "Iterator " << iteration << ": \t the absolute err is "
                      << std::sqrt(gpu_VectorNorm_optimized(r, LEN)) << std::endl;

        if (std::sqrt(gpu_VectorNorm_optimized(r, LEN)) < 1.e-3) {
            std::cout << "The CG process has been Convergent, and CG method is "
                         "completed!"
                      << std::endl;
            break;
        }
    }

    Delte1D(r);
    Delte1D(r_tmp);
    Delte1D(tmp);
    Delte1D(p_all);
}