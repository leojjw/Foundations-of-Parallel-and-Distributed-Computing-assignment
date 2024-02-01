/* 
  Foundations of Parallel and Distributed Computing , Fall 2023.
  Instructor: Prof. Chao Yang @ Peking University.
  This is a naive implement of heat propagation equation.
  Version 0.1 (date: 10/30/2023)
*/
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <omp.h>
# include <time.h>

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
    int iterations = 0;
    int iterations_print = 1;
    double diff = epsilon;
    clock_t time0, time1;
    double duration;
    
    // init    
    init_data(t);
    init_mean(t, mean);
    init_interior(t, mean);
    
    // solve    
    time0 = clock();

    while(epsilon <= diff) {

        for (int i = 0; i < M; i++ ) {
            for (int j = 0; j < N; j++ ) {
                odt[i][j] = t[i][j];
            }
        }
        for (int i = 1; i < M - 1; i++ ) {
            for (int j = 1; j < N - 1; j++ ) {
                t[i][j] = ( odt[i-1][j] + odt[i+1][j] + odt[i][j-1] + odt[i][j+1] ) / 4.0;
            }
        }

        double my_diff = 0.0;
        diff = 0.0;

        for (int i = 1; i < M - 1; i++ ) {
            for (int j = 1; j < N - 1; j++ ) {
                if ( my_diff < fabs ( t[i][j] - odt[i][j] ) ) {
                    my_diff = fabs ( t[i][j] - odt[i][j] );
                }
            }
        }

        if ( diff < my_diff ) {
            diff = my_diff;
        }

        // print
        iterations++;
        if ( iterations == iterations_print ) {
            printf ( "  %8d  %f\n", iterations, diff);
            iterations_print = 2 * iterations_print;
        }
    }
    duration = (double)(clock()- time0) / CLOCKS_PER_SEC;
    printf ( "\n" );
    printf ( "  %8d  %f\n", iterations, diff );
    printf ( "\n" );
    printf ( "  Error tolerance achieved.\n" );
    printf ( "  Wallclock time = %f\n", duration );
    return 0;
}