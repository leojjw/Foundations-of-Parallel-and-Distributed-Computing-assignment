#include <Setting.hpp>
void initialization();
void finalization();
void assemble();
void CG(double** A,
        double* x,
        double* b,
        unsigned int Xnode,
        unsigned int Ynode);
void Output();

int main() {
    initialization();
    assemble();

    CG(A, x, b, Xnode, Ynode);

    Output();

    finalization();
    return 0;
}