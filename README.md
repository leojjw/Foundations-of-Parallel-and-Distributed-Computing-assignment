**// load mpi**

module add mpich

**// for c**

mpicc mpi_hello.c -o hello

**// for cpp**

mpicxx mpi_hello.cpp -o hello

**// run**

mpiexec -n 4 ./hello
