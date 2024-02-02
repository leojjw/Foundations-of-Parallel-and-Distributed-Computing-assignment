
#include <ArrayAllocate/ArrayAllocate.hpp>
#include <Vector/Gemm.hpp>
#include <Vector/VectorOperator.hpp>
#include <cmath>
#include <iostream>
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
    double* p = Allocate1D<double>(LEN);
    for (int i = 0; i < LEN; i++) {
        r[i] = x[i] = 0;
    }

    /*initializition*/
    /*tmp=A*x0*/
    Gemm<double>(A, x, tmp, LEN, LEN);
    /*r0=b-tmp*/
    VectorOperation<double, MINUS>(b, tmp, r, LEN);
    /*p0=r0*/
    VectorOperation<double, ASIG>(r, nullptr, p, LEN);

    int iteration = 0;
    while (true) {
        /*the loop exit condition*/
        if (iteration % 10 == 0)
            std::cout << "Iterator " << iteration << ": \t the absolute err is "
                      << std::sqrt(VectorNorm(r, LEN)) << std::endl;

        if (std::sqrt(VectorNorm(r, LEN)) < 1.e-6) {
            std::cout << "The CG process has been Convergent, and CG method is "
                         "completed!"
                      << std::endl;
            break;
        }

        iteration++;

        auto tmp_v = VectorNorm(r, LEN);
        Gemm(A, p, tmp, LEN, LEN);
        auto alpha = tmp_v / VectorDotVector(p, tmp, LEN);

        /*x_k+1=x_k+alpha*p*/
        VectorPlusKVector(x, p, x, alpha, LEN);

        VectorPlusKVector(r, tmp, r_tmp, -alpha, LEN);

        auto beta = VectorNorm(r_tmp, LEN) / VectorNorm(r, LEN);

        VectorPlusKVector(r_tmp, p, p, beta, LEN);

        VectorOperation<double, ASIG>(r_tmp, nullptr, r, LEN);
    }

    Delte1D(r);
    Delte1D(r_tmp);
    Delte1D(tmp);
    Delte1D(p);
}