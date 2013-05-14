/*! \file simplex.hpp
 * Simplex Algorithm - straightforward implementation
 */
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

#pragma once

void fnc_simplex(double * x, double * c, double * A, double * b, int n, int m);
void fnc_read_lp(double *& x, double *& c, double *& A, double *& b, int & n, int & m, std::string filename);

template<class T>
void fnc_printvec(T * vec, int n) {
    std::cout << "[ ";
    for(int i = 0; i < n; ++i) {
        std::cout << std::setprecision(3) << std::fixed << vec[i] << " ";
    }
    std::cout << "]";
}

template<class T>
std::string fnc_strvec(T * vec, int n) {
    std::stringstream output;
    output << "[ ";
    for(int i = 0; i < n; ++i) {
        output << std::setprecision(3) << std::fixed << vec[i] << " ";
    }
    output << "]";
    return output.str();
}
