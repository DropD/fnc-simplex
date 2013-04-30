/*******************************************************************************
 * 
 * Simplex for standard feasible maximisation LPs
 * 
 * 
 * tableau format:
 * 
 * x1 y1 z1 s1 t1  0 | B1
 * x2 y2 z2 s2 t2  0 | B2
 * c1 c2 c3  0  0  1 |  0
 * 
 ******************************************************************************/

/*
Assumptions:
  We have a feasible maximisation problem
  'Maximize' comes before 'Subject To'
  No newlines between 'Maximize'/'Subject To' and data, or within data
  All constraints are <=
  Variables determined by position (0 for empty)
  Variables are >= 0
*/



#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <unistd.h>

#include "DebugPrinter.hpp"

template <typename T>
class Simplex {

    public:

    std::vector< std::vector<T> > data;
    std::vector< int > active;
    int n, m, width;

    void load(std::string fname) {

        std::ifstream fp(fname.c_str());
        if( !fp.is_open() ) throw("failed to open file " + fname);


        std::vector<T> costs;
        std::vector< std::vector<T> > constraints;

        std::string line;
        while( !fp.eof() && fp.good() ) {

            getline(fp, line);
            if(line == "Maximize") {

                getline(fp, line);                     // get objective function
                costs = split_vars(line);              // split_vars is size n

            } else if(line == "Subject To") {

                while( !fp.eof() && fp.good() ) {
                    getline(fp, line);                 // get the constraints
                    if(line == "") break;
                    constraints.push_back(split_vars(line));   // split_vars is size n+1
                }

            }

        }

        n = costs.size();
        m = constraints.size();
        width = n+1 + m+1;
        data = std::vector< std::vector<T> >
            (m+1, std::vector<T>( width ));
        active = std::vector<int>(m);

        data[m][n+m] = 1.;                               // set last row
        for(int i = 0; i < n; ++i)
            data[m][i] = -costs[i];

        for(int i = 0; i < m; ++i)  {                    // set the rest
            for(int j = 0; j < n; ++j) {
                data[i][j] = constraints[i][j];
            }
            data[i][width-1] = constraints[i][n];        // b vector
            data[i][n+i] = 1;                            // diagonal slack
            active[i] = n+i;
        }


    }


    std::vector<double> split_vars(const std::string str) {

        std::vector<double> v;
        std::istringstream ss(str);
        double d;
        while(ss >> d) { v.push_back(d); }
        return v;

    }

    void print() {

        std::cout << "Tableau for " << n << " variables and "
                  << m << " constraints:" << std::endl;
        for(int i = 0; i < data.size(); ++i) {
            for(int j = 0; j < width; ++j) {
                std::cout << "\t" << data[i][j];
            }
            std::cout << std::endl;
        }

    }

    void solve() {

        while(true) {

            int col = pivot_col();
            if(col >= width) break;  // no negative -> finished
            int row = pivot_row(col);
            if(row >= m) throw("invalid pivot row");

            T pivot = data[row][col];
            for(int i = 0; i < m+1; ++i) {  // clear column
                if(i == row)
                    continue;
                T fac = data[i][col]/pivot;   // TODO: check for instability
                for(int j = 0; j < width; ++j)
                    data[i][j] -= fac*data[row][j];
            }
            active[row] = col;

            #ifdef VERBOSE
                std::cout << "Pivoting element (" << row << ", " << col << ")" << std::endl;
                print();
                usleep(0.2*1000*1000);
            #endif // VERBOSE
        }

    }


    int pivot_col() {
        T min = 0;              // min must be negative
        int idx = width;
        for(int i = 0; i < width; ++i)      // TODO: move into comparer struct
            if(data[m][i] < min) {
                min = data[m][i];
                idx = i;
            }
        return idx;
    }

    int pivot_row(int x) {
        T min = 0;              // min can be arbitrary
        bool any = false;
        int idx = m;
        for(int i = 0; i < m; ++i) {
            if(data[i][x] <= 0)
                continue;
            T ratio = data[i][width-1]/data[i][x];
            if(ratio <= min || any == false) {
                min = ratio;
                idx = i;
                any = true;
            }
        }
        return idx;
    }

    std::vector<T> solutions() {
        std::vector<T> sol(2*n+1);
        sol[0] = data[m][width-1] / data[m][2*n];
        for(int i = 0; i < m; ++i)
            sol[active[i]+1] = data[i][width-1] / data[i][active[i]];
        return sol;
            
    }

};
