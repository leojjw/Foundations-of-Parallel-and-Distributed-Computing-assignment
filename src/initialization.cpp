#include <ArrayAllocate/ArrayAllocate.hpp>
#include <Setting.hpp>

void initialization() {
    /*determine the size of matrix*/
    unsigned int SIZE = Xnode * Ynode;

    /*allcate matrix A*/
    A = Allocate2D<double>(SIZE, SIZE);

    /*allocate vector X*/
    x = Allocate1D<double>(SIZE);

    /*allocate vector b*/
    b = Allocate1D<double>(SIZE);
}