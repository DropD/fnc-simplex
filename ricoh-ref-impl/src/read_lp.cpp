#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <exception>
#include <stdlib.h>

#include "../include/simplex.hpp"

class NonMinException : public std::exception {
    virtual const char* what() const throw() {
        return "Only minimization implemented";
    }
} nonMinException;

std::vector<double> split_vars(const std::string str) {

    std::vector<double> v;
    std::istringstream ss(str);
    double d;
    while(ss >> d) { v.push_back(d); }
    return v;

}

void fnc_read_lp(double *& x, double *& c, double *& A, double *& b, int & n, int & m, std::string filename) {
    std::ifstream fp(filename.c_str());

    if( !fp.is_open() ) throw("failed to open file " + filename);


    std::vector<double> costs;
    std::vector< std::vector<double> > constraints;


    std::string line;
    while( !fp.eof() && fp.good() ) {

        getline(fp, line);
        if(line == "Minimize") {

            getline(fp, line);                     // get objective function
            costs = split_vars(line);              // split_vars is size n

        } else if(line == "Subject To") {

            while( !fp.eof() && fp.good() ) {
                getline(fp, line);                 // get the constraints
                if(line == "") break;
                constraints.push_back(split_vars(line));   // split_vars is size n+1
            }

        } else if (line == "Maximize") {
            throw nonMinException;
        }

    }


    m = constraints.size();
    n = costs.size() + m;


    // set cost vector
    // and allocate and initialize x
    c = new double[n];
    x = new double[n];
    for(int i = 0; i < n - m; ++i) {
        c[i] = costs[i];
        x[i] = 0.;
    }
    for(int i = n - m; i < n; ++i) {
        c[i] = 0.;
        x[i] = 0.;
    }

    // set constraint matrix and b vector
    A = new double[n * m];
    b = new double[m];
    for(int i = 0; i < m; ++i)  {                    // set the rest
        // constraint coeffs
        for(int j = 0; j < n - m; ++j) {
            A[i * n + j] = constraints[i][j];
        }
        // slack identity matrix
        for(int j = n - m; j < n; ++j) {
            A[i * n + j] = (i == j - n + m) ? 1 : 0;
        }
        b[i] = constraints[i][n-m];        // b vector
    }
}
