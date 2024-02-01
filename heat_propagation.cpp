/* 
  Foundations of Parallel and Distributed Computing , Fall 2023.
  Instructor: Prof. Chao Yang @ Peking University.
  This is a naive implement of heat propagation equation.
  Version 0.1 (date: 10/30/2023)
*/
# include "mpi.h"
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>
# include <time.h>
# include <cstring>

#define ROOT 0

const int M = 6000;
const int N = 8000;
const double epsilon = 0.01;

void init_data(double (&t) [M][N]) {
    for (int i = 1; i < M - 1; i++) {
        t[i][0] = 100.0;
        t[i][N-1] = 100.0;
    }
    for (int j = 0; j < N; j++) {
        t[M-1][j] = 100.0;
        t[0][j] = 0.0;
    }
}

void init_mean(double t[M][N], double& mean) {
    for (int i = 1; i < M - 1; i++ ) {
        mean = mean + t[i][0] + t[i][N-1];
    }
    for (int j = 0; j < N; j++ ) {
        mean = mean + t[M-1][j] + t[0][j];
    }
    mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
    printf ("  MEAN = %f\n", mean );
}

void init_interior(double (&t) [M][N], double mean) {
    for (int i = 1; i < M - 1; i++ ) {
        for (int j = 1; j < N - 1; j++ ) {
            t[i][j] = mean;
        }
    }
}

int main (int argc, char *argv[]) {
    // param
    double odt[M][N];
    double t[M][N];
    double mean = 0.0;
    double t0, t1;
    int iterations = 0;
    int iterations_print = 1;
    double diff = epsilon, my_diff_reduce;
    int size, rank;
    int size_sqrt;
    
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // init   
    init_data(t);
    init_mean(t, mean);
    init_interior(t, mean);

    // solve    
    t0 = MPI_Wtime();

    size_sqrt = sqrt(size);

    for (int i = 0; i < M; i++ ) {
        for (int j = 0; j < N; j++ ) {
            odt[i][j] = t[i][j];
        }
    }

    while(epsilon <= diff) {

        memset(t, 0, sizeof(t));
        
        for (int i = M / size_sqrt * (rank / size_sqrt); i < M / size_sqrt * (rank / size_sqrt + 1); i++ ) {
            for (int j = N / size_sqrt * (rank % size_sqrt); j < N / size_sqrt * (rank % size_sqrt + 1); j++ ) {
                if (i == M-1 || (i != 0 && (j == 0 || j == N-1))){
                    t[i][j] = 100.0;
                    continue;
                }
		        if (i == 0) continue;
                t[i][j] = ( odt[i-1][j] + odt[i+1][j] + odt[i][j-1] + odt[i][j+1] ) / 4.0;
            }
        }

        double my_diff = 0.0;
        diff = 0.0;

        for (int i = M / size_sqrt * (rank / size_sqrt); i < M / size_sqrt * (rank / size_sqrt + 1); i++ ) {
            for (int j = N / size_sqrt * (rank % size_sqrt); j < N / size_sqrt * (rank % size_sqrt + 1); j++ ) {
                if (i == 0 || i == M-1 || j == 0 || j == N-1) continue;
                if ( my_diff < fabs ( t[i][j] - odt[i][j] ) ) {
                    my_diff = fabs ( t[i][j] - odt[i][j] );
                }
            }
        }

    	MPI_Allreduce(t, odt, M*N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&my_diff, &my_diff_reduce, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	    if (diff < my_diff_reduce) {
	        diff = my_diff_reduce;
	    }

        // print
        if (rank == ROOT){
            iterations++;
            if (iterations == iterations_print ) {
                printf ( "  %8d  %f\n", iterations, diff);
                iterations_print = 2 * iterations_print;
            }
        }
    }
    t1 = MPI_Wtime();
    if (rank == ROOT){
        printf ( "\n" );
        printf ( "  %8d  %f\n", iterations, diff );
        printf ( "\n" );
        printf ( "  Error tolerance achieved.\n" );
        printf ( "  Wallclock time = %f\n", t1-t0 );
    }

    MPI_Finalize();
    
    return 0;
}
