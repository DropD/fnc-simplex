#include <iostream>

#include "../include/simplex.hpp"

void fnc_mvmm(double * a, double * b, double * C, double * D, double * e, int n, int m, int * I);
void fnc_mvm(double * a, double * B, double * C, int n, int m, int s);

void fnc_simplex(double * x, double * c, double * A, double * b, int n, int m) {
    double eps = 1e-9;

    double *Ai, *cb;
    Ai = new double[m*m];
    cb = new double[m];

    int *B;
    B = new int[m];

    // initial basis
    int i, j;
    for(i = 0; i < m; ++i) {
        B[i] = i + n - m;
    }

    for(i = 0; i < m; ++i) {
        for(j = 0; j < m; ++j) {
            Ai[i*m + j] = (i == j ? 1 : 0);
        }
    }

    for(i = 0; i < n-m; ++i) {
        x[i] = 0.;
    }
    for(i = 0; i < m; ++i) {
        x[i + n - m] = b[i];
    }

    int go_on = 100;

    // temp for $\bar{c}$
    double *p, *u;
    p  = new double[m];
    u  = new double[m];
    double d, xbu, min, ut;
    int s, count, t, tmp;

    while(go_on) {
        // calculate $\bar{c}' = c' - c_B' A_B^-1 A$
        fnc_mvmm(cb, c, Ai, A, p, n, m, B);

        // optimality condition
        // s corresponds j in ZF
        s = -1;
        for(i = 0; i < n; ++i) {
            if(cb[i] < -eps) {
                s = i;
                break;
            }
        }
        if(s == -1) {
            printf("done.\n");
            return;
        }
        
        // calculate $u = A_B^-1 * A[:,s]$
        fnc_mvm(u, Ai, A, n, m, s);

        // unbounded ?
        count = 0;
        for(i = 0; i < m; ++i) {
            count += u[i] > eps ? 1 : 0;
        }
        if(count == 0) {
            printf("unbounded\n");
            return;
        }

        // find t (B[t] corresponds to l in Zf) and d
        min = -1;
        for(i = 0; i < m; ++i) {
            if(u[i] > eps) {
                xbu = x[B[i]] / u[i];
                if(xbu < min || min == -1) {
                    t = i;
                    min = xbu;
                }
            }
        }
        d = min;

        // update x
        for(i = 0; i < m; ++i) {
            x[B[i]] -= d * u[i];
        }
        x[s] = d;

        // debug output
#ifdef VERBOSE
        //printf("B = [");
        //for(i = 0; i < m; ++i) {
        //    printf(" %d", B[i]);
        //}
        //printf(" ] ");
        std::cout << "B = ";
        fnc_printvec(B, m);
        std::cout << " "
                  << "Bt = " << B[t] << " "
                  << "s = " << s << " ";

        std::cout << "u = ";
        fnc_printvec(u, m);
        std::cout << " ";

        std::cout << "x = ";
        fnc_printvec(x, n);
        std::cout << " ";

        std::cout << "cb = ";
        fnc_printvec(cb, n);
        std::cout << std::endl;

#endif

        // update B
        B[t] = s;
        //// sort B
        //for(i = 0; i < n - m - 1; ++i) {
        //    if(B[i+1] < B[i]) {
        //        tmp = B[i+1];
        //        B[i+1] = B[i];
        //        B[i] = tmp;
        //    }
        //}

        // update $A_B^{-1}$
        ut = 1. / u[t];
        for(i = 0; i < m; ++i) {
            for(j = 0; j < m; ++j) {
                if(i != t) {
                    Ai[i*m + j] -= u[i] * ut * Ai[t*m + j];
                }
            }
        }
        for(i = 0; i < m; ++i) {
            Ai[t*m + i] *= ut;
        }

#ifdef VERBOSE
        for(i = 0; i < m; ++i) {
            std::cout << "Ai[" << i << "] = ";
            fnc_printvec(Ai + i*m, m);
            std::cout << std::endl;
        }
#endif

        --go_on;
    }
    /*
    */
}

// a' = b' - b_I' * C * D
void fnc_mvmm(double * a, double * b, double * C, double * D, double * e, int n, int m, int * I) {
    int i, j, k;
    for(i = 0; i < m; ++i) {
        double sum = 0;
        for(j = 0; j < m; ++j) {
            sum += b[I[j]] * C[j * m + i];
        }
        e[i] = sum;
    }

    for(i = 0; i < n; ++i) {
        double sum = 0;
        for(j = 0; j < m; ++j) {
            sum += e[j] * D[j * n + i];
        }
        a[i] = b[i] - sum;
    }
}

// a = B C[:,j]
void fnc_mvm(double * a, double * B, double * C, int n, int m, int s) {
    int i, j;
    for(i = 0; i < m; ++i) {
        double sum = 0;
        for(j = 0; j < m; ++j) {
            sum += B[i * m + j] * C[j * n + s];
        }
        a[i] = sum;
    }
}