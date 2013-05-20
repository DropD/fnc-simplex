/*! \file main.cpp
 * Simplex Application Main
 * solves a linear program in standard form
 * i.e.
 * max c' x
 * A x = b
 * x >= 0
 */
#include "../include/simplex.hpp"
#include <stdlib.h>
#include <stdio.h>

int main (int argc, char const* argv[])
{
    double *c, *A, *b, *x;
    // test problem.
    int n = 6, m = 3;
    c = (double*)malloc(n*sizeof(double));
    A = (double*)malloc(n*m*sizeof(double));
    b = (double*)malloc(m*sizeof(double));
    x = (double*)malloc(n*sizeof(double));

    //c = { -3., -9., 0., 0., 0. };
    c[0] = 2.;
    c[1] = 1.;
    c[2] = 2.;
    c[3] = 0.;
    c[4] = 0.;
    c[5] = 0.;
    //A = {  3.,  1., 1., 0., 0.,
    //       2.,  3., 0., 1., 0.,
    //      -2.,  3., 0., 0., 1. };
    A[0] = 1.;  A[6] = -1.; A[12] = -2.;
    A[1] = 5.;  A[7] = -2.; A[13] = -4.; 
    A[2] = 1.;  A[8] = -1.; A[14] = -1.; 
    A[3] = 1.;  A[9] =  0.; A[15] =  0.; 
    A[4] = 0.; A[10] =  1.; A[16] =  0.; 
    A[5] = 0.; A[11] =  0.; A[17] =  1.; 
    //b = { 15., 18., 6.};
    b[0] = 100.;
    b[1] = 50.;
    b[2] = 80.;

    fnc_simplex(x, c, A, b, n, m);

    printf("x = ");
    for(int i = 0; i < n; ++i){
        printf("%lf ", x[i]);
    }
    printf("\n");
    
    return 0;
}
