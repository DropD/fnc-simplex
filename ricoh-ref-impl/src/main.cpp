/*! \file main.cpp
 * Simplex Application Main
 * solves a linear program in standard form
 * i.e.
 * max c' x
 * A x = b
 * x >= 0
 */
#include "../include/simplex.hpp"
#include "../include/time_fnc_simplex.hpp"
#include <iostream>
#include <iomanip>
#include <exception>

int main (int argc, char const* argv[])
{
    double *c, *A, *b, *x;
    double r, freq = 2.2e9;

    int n, m;
    std::string lp_file;

    if(argc != 2) {
        std::cout << "usage: lp <inputfile>" << std::endl;
        return 0;
    }
    else {
        lp_file = argv[1];
    }
    
    try {
        std::cout << "reading input file ... ";
        fnc_read_lp(x, c, A, b, n, m, lp_file);
        std::cout << "done." << std::endl;

        std::cout << "optimizing ... ";
        r = time_fnc_simplex(x, c, A, b, n, m);

        std::cout << "x = ";
        fnc_printvec(x, n);
        std::cout << std::endl;
    }
    catch(std::exception & e) {
        std::cout << e.what() << std::endl;
    }

    double obj = 0;
    for(int i = 0; i < n-m; ++i) {
        obj += x[i] * c[i];
    }
    std::cout << "Objective value: " << obj << std::endl;

    std::cout << std::scientific 
              << "computed in " << r << " cycles" 
              << " or " << r / freq << " s" 
              << std::endl;
    
    return 0;
}
