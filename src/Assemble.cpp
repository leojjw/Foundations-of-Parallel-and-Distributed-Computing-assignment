#include <Setting.hpp>
#include <iostream>
void PrintArray2D(double** A, int M, int N);
void assemble() {
    /*set b*/
    for (int i = 0; i < Xnode; i++) {
        for (int j = 0; j < Ynode; j++) {
            int index = i * Ynode + j;
            if (i == 0 || i == Xnode - 1 || j == 0 || j == Ynode - 1)
                b[index] = 0;
            else b[index] = -h * h * Q;
        }
    }

    /*Set A*/
    for (int i = 0; i < Xnode; i++) {
        for (int j = 0; j < Ynode; j++) {
            int index = i * Ynode + j;
            if (i == 0 || i == Xnode - 1 || j == 0 || j == Ynode - 1)
                A[index][index] = 1;
            else {
                A[index][index - 1] = 1;
                A[index][index + 1] = 1;
                A[index][index] = -4;
                A[index][index - Xnode] = 1;
                A[index][index + Xnode] = 1;
            }
        }
    }
    std::cout << "Matrix A and vector b have been assembeled" << std::endl;
}

void PrintArray2D(double** A, int M, int N) {
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++)
            std::cout << A[i][j] << '\t';
        std::cout << std::endl;
    }
}