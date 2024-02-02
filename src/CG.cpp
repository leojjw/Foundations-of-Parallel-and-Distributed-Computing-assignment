#include <omp.h> 
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

    /*declare paramters*/
    double* r = Allocate1D<double>(LEN);
    double* r_tmp = Allocate1D<double>(LEN);
    double* tmp = Allocate1D<double>(LEN);
    double* p = Allocate1D<double>(LEN);

    double t0, t1;
    int nthreads, tlen, i;
    
    t0 = omp_get_wtime();

    #pragma omp parallel default(shared)
    {
        nthreads = omp_get_num_threads();
        tlen = LEN / nthreads;

        #pragma omp for
            for (i = 0; i < LEN; i++) {
                r[i] = x[i] = 0;
            }

        /*initializition*/
        #pragma omp for
            for (i = 0; i < LEN ; i += tlen) {
                /*tmp=A*x0*/
                Gemm<double>(A + i, x, tmp + i, (i + tlen) < LEN ? tlen : (LEN - i), LEN);
                /*r0=b-tmp*/
                VectorOperation<double, MINUS>(b + i, tmp + i, r + i, (i + tlen) < LEN ? tlen : (LEN - i));
                /*p0=r0*/
                VectorOperation<double, ASIG>(r + i, nullptr, p + i, (i + tlen) < LEN ? tlen : (LEN - i));
            }

    }

    int iteration = 0;

    while (true) {
        /*the loop exit condition*/
        if (iteration % 10 == 0) {
            double error = 0.0;
            #pragma omp parallel for reduction(+, error)
                for (i = 0; i < LEN; i += tlen) {
                    error += VectorNorm(r + i, (i + tlen) < LEN ? tlen : (LEN - i));
                }
            std::cout << "Iterator " << iteration << ": \t the absolute err is "
                      << std::sqrt(error) << std::endl;
        }

        if (std::sqrt(VectorNorm(r, LEN)) < 1.e-6) {
            double error = 0.0;
            #pragma omp parallel for reduction(+, error)
                for (i = 0; i < LEN; i += tlen) {
                    error += VectorNorm(r + i, (i + tlen) < LEN ? tlen : (LEN - i));
                }
            std::cout << "Iterator " << iteration << ": \t the absolute err is "
                      << std::sqrt(error) << std::endl;
            std::cout << "The CG process has been Convergent, and CG method is "
                         "completed!"
                      << std::endl;
            break;
        }

        iteration++;
        double tmp_v, alpha, beta_tmp1, beta_tmp2;
        #pragma omp parallel default(shared)
        {
            tmp_v = 0.0, alpha = 0.0, beta_tmp1 = 0.0, beta_tmp2 = 0.0;
            #pragma omp for reduction(+:tmp_v, alpha)
                for (i = 0; i < LEN; i += tlen) {
                    tmp_v += VectorNorm(r + i, (i + tlen) < LEN ? tlen : (LEN - i));
                    Gemm(A + i, p, tmp + i, (i + tlen) < LEN ? tlen : (LEN - i), LEN);
                    alpha += VectorDotVector(p + i, tmp + i, (i + tlen) < LEN ? tlen : (LEN - i));
                }
            alpha = tmp_v / alpha;

            #pragma omp for reduction(+:beta_tmp1, beta_tmp2)
                for (i = 0; i < LEN ; i += tlen) {
                    VectorPlusKVector(x + i, p + i, x + i, alpha, (i + tlen) < LEN ? tlen : (LEN - i));
                    VectorPlusKVector(r + i, tmp + i, r_tmp + i, -alpha, (i + tlen) < LEN ? tlen : (LEN - i));
                    beta_tmp1 += VectorNorm(r_tmp + i, (i + tlen) < LEN ? tlen : (LEN - i));
                    beta_tmp2 += VectorNorm(r + i, (i + tlen) < LEN ? tlen : (LEN - i));
                }

            auto beta = beta_tmp1 / beta_tmp2;
            
            #pragma omp for
                for (i = 0; i < LEN ; i += tlen) {
                    VectorPlusKVector(r_tmp + i, p + i, p + i, beta, (i + tlen) < LEN ? tlen : (LEN - i));
                    VectorOperation<double, ASIG>(r_tmp + i, nullptr, r + i, (i + tlen) < LEN ? tlen : (LEN - i));
                }
        }
    }

    t1 = omp_get_wtime();
    std::cout << "Wall clock time = " << t1 - t0 << std::endl;

    Delte1D(r);
    Delte1D(r_tmp);
    Delte1D(tmp);
    Delte1D(p);
}
