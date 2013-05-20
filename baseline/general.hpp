/*******************************************************************************
 * 
 * Simplex for standard feasible maximisation LPs
 * 
 * 
 * tab format:
 * 
 * x1 y1 z1 s1 t1  0 | B1
 * x2 y2 z2 s2 t2  0 | B2
 * c1 c2 c3  0  0  1 |  0
 * 
 ******************************************************************************/

/*
Assumptions:
  We have a maximisation problem
  'Maximize' comes before 'Subject To'
  No newlines between 'Maximize'/'Subject To' and the data, or within the data
  All constraint inequalities are <= and rhs is positive
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
#include "rdtsc.h"

//#include "DebugPrinter.hpp"

template <typename T>
class Simplex {

    private:

    std::vector< int > active;
    std::vector< int > nonstandard;

    std::vector< std::vector<T> > tab;
    int n, m, width;
    public:


    void load(std::string fname) {

        std::ifstream fp(fname.c_str());
        if( !fp.is_open() ) {
            std::cerr << "failed to open file " + fname << std::endl;
            throw;
        }


        std::vector<T> costs;
        std::vector< std::vector<T> > constraints;

        std::string line;
        while( !fp.eof() && fp.good() ) {

            getline(fp, line);
            if(line == "Maximize") {

                getline(fp, line);                     // get objective function
                costs = split_vars(line);              // split_vars is size n

            //~ } else if(line == "Minimize") {
//~ 
                //~ getline(fp, line);
                //~ costs = split_vars(line);
                //~ for(int i = 0; i < costs.size(); ++i)
                    //~ costs[i] *= -1;

            } else if(line == "Subject To") {

                int nr = 0;
                while( !fp.eof() && fp.good() ) {
                    getline(fp, line);                 // get the constraints
                    if(line == "") break;
                    constraints.push_back(split_vars(line, nr));   // split_vars is size n+1
                    ++nr;
                }

            }

        }

        n = costs.size();
        m = constraints.size();
        width = n+1 + m+1;
        tab = std::vector< std::vector<T> >
            (m+1, std::vector<T>( width ));
        active = std::vector<int>(m);

        tab[m][n+m] = 1.;                               // set last row
        for(int i = 0; i < n; ++i)
            tab[m][i] = -costs[i];

        for(int i = 0; i < m; ++i)  {                    // set the rest
            for(int j = 0; j < n; ++j) {
                tab[i][j] = constraints[i][j];
            }
            tab[i][width-1] = constraints[i][n];        // b vector
            tab[i][n+i] = 1;                            // diagonal slack
            active[i] = n+i;
        }
        for(int i = 0; i < nonstandard.size(); ++i)
            tab[nonstandard[i]][n+nonstandard[i]] *= -1;


    }


    std::vector<double> split_vars(const std::string str, int nr = -1) {

        std::vector<double> v;
        std::istringstream ss(str);
        T d;
        do {
            for(T d; ss >> d; )
                v.push_back(d);
            if(ss.fail()) {
                ss.clear();
                std::string s;
                ss >> s;
                if(s == "<=") {}
                else if(s == ">=" && nr >= 0) {
                    nonstandard.push_back(nr);
                }
                else if(s.size()) {
                    std::cerr << "Bad input constraints" << std::endl;
                    throw;
                }
            }
        } while(!ss.eof());
        return v;

    }

    void print() {

        std::cout << "Tableau for " << n << " variables and "
                  << m << " constraints:" << std::endl;
        for(int i = 0; i < tab.size(); ++i) {
            for(int j = 0; j < width; ++j) {
                std::cout << "\t" << tab[i][j];
            }
            std::cout << std::endl;
        }

    }

    void solve() {

        for(int i = 0; i < nonstandard.size(); ++i) {

            int col = pivot_col_nonstandard(nonstandard[i]);
            if(col >= width) break;  // no negative -> finished
            int row = pivot_row(col);
            if(row >= m) {
                std::cerr << "Invalid pivot row (problem infeasible)" << std::endl;
                throw;
            }
            basis_exchange(row, col);

            #ifdef VERBOSE
                std::cout << "Pivoting NONSTANDARD element (" << row << ", " << col << ")" << std::endl;
                print();
                std::cout << "active:";
                for(int i = 0; i < active.size(); ++i)
                    std::cout << " " << active[i];
                std::cout << std::endl;
                usleep(0.2*1000*1000);
            #endif // VERBOSE

        }

        while(true) {

            int col = pivot_col();
            if(col >= width) break;  // no negative -> finished
            int row = pivot_row(col);
            if(row >= m) {
                std::cerr << "Invalid pivot row (problem infeasible)" << std::endl;
                throw;
            }

            basis_exchange(row, col);

            #ifdef VERBOSE
                std::cout << "Pivoting element (" << row << ", " << col << ")" << std::endl;
                print();
                std::cout << "active:";
                for(int i = 0; i < active.size(); ++i)
                    std::cout << " " << active[i];
                std::cout << std::endl;
                usleep(0.2*1000*1000);
            #endif // VERBOSE
        }

    }


    int pivot_col() {
        T min = 0;              // min must be negative
        int idx = width;
        for(int i = 0; i < width-1; ++i)      // TODO: move into comparer struct
            if(tab[m][i] < min) {
                min = tab[m][i];
                idx = i;
            }
        return idx;
    }

    int pivot_col_nonstandard(int row) { // split from pivot_col since different optimisations might apply
        T max = 0;
        int idx = width;
        for(int i = 0; i < width-1; ++i)
            if(tab[row][i] > max) {
                max = tab[row][i];
                idx = i;
            }
        return idx;
    }

    int pivot_row(int x) {
        T min = 0;              // min can be arbitrary
        bool any = false;
        int idx = m;
        for(int i = 0; i < m; ++i) {
            if(tab[i][x] <= 0)
                continue;
            T ratio = tab[i][width-1]/tab[i][x];
            if(ratio <= min || any == false) {
                min = ratio;
                idx = i;
                any = true;
            }
        }
        return idx;
    }

    inline void basis_exchange(int row, int col) {
        T pivot = tab[row][col];
        for(int i = 0; i < m+1; ++i) {
            if(i == row)
                continue;
            T fac = tab[i][col]/pivot;   // TODO: check for instability
            for(int j = 0; j < width; ++j)
                tab[i][j] -= fac*tab[row][j];
        }
        active[row] = col;
    }

    std::vector<T> solutions() {
        std::vector<T> sol(2*n+1);
        sol[0] = tab[m][width-1] / tab[m][2*n];
        for(int i = 0; i < m; ++i)
            sol[active[i]+1] = tab[i][width-1] / tab[i][active[i]];
        return sol;
    }

    int problem_size() { return n; }

};
