/*******************************************************************************
 * 
 * Simplex base class
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
  We have a maximisation problem
  'Maximize' comes before 'Subject To'
  No newlines between 'Maximize'/'Subject To' and the respective data
  'Subject To' is closed with either 'End' or EOF
  ** All constraint inequalities are <= and rhs is positive
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
#include <cstring>

#include "../misc/DebugPrinter.hpp"

typedef unsigned int uint;

template <typename T>
class SimplexBase {

    protected:

    std::vector< int > active, backup_active;
    std::vector< int > nonstandard;

    std::vector< std::vector<T> > tab, backup_tab;
    T * tabp, * backup_tabp;
    int n, m, width;
    int iter = 0;

    inline std::vector<double> split_vars(const std::string str, int nr = -1) {

        std::vector<double> v;
        std::istringstream ss(str);
        do {
            for(T d; ss >> d; )
                v.push_back(d);
            if(ss.fail()) {
                ss.clear();
                std::string s;
                ss >> s;
                if(s == "<=") {}
                else if(s == ">=" && nr >= 0) {
                    nonstandard.push_back(nr);    // index nonstd constraints
                }
                else if(s.size()) {
                    std::cerr << "Bad input constraints" << std::endl;
                    throw;
                }
            }
        } while(!ss.eof());
        return v;

    }



    public:

    unsigned int PERFC_MEM, PERFC_ADDMUL, PERFC_DIV;

    virtual void load(std::string fname) {

        std::ifstream fp(fname.c_str());
        if( !fp.is_open() ) {
            std::cerr << "failed to open file " + fname << std::endl;
            throw;
        }

        iter = -1;
        PERFC_MEM = 0;
        PERFC_ADDMUL = 0;
        PERFC_DIV = 0;

        std::vector<T> costs;
        std::vector< std::vector<T> > constraints;
        nonstandard.clear();

        int nr=0;
        std::string line;
        bool constraint_section = false;
        while( !fp.eof() && fp.good() ) {              // read line by line

            getline(fp, line);
            if(line == "Maximize") {
                getline(fp, line);                     // get objective function
                costs = split_vars(line);              // split_vars is size n
            } else if(line == "Minimize") {
                std::cerr << "Minimisation not supported" << std::endl;
                throw;
                //~ getline(fp, line);
                //~ costs = split_vars(line);
                //~ for(int i = 0; i < costs.size(); ++i)
                    //~ costs[i] *= -1;
            } else if(line == "Subject To") {
                nr = 0;
                constraint_section = true;
            } else if(line == "End") {
                break;
            } else if(constraint_section == true && line != "") {
                constraints.push_back(split_vars(line, nr));   // split_vars is size n+1
                ++nr;   // track current constraint
            }

        }

        n = costs.size();
        m = constraints.size();
        width = n+1 + m+1;
        tab = std::vector< std::vector<T> >
            (m+1, std::vector<T>( width ));
        active = std::vector<int>(m);

        tab[m][n+m] = 1.;                               // set last row (obj)
        for(int i = 0; i < n; ++i)
            tab[m][i] = -costs[i];

        for(int i = 0; i < m; ++i)  {                    // set the rest
            for(int j = 0; j < n; ++j) {
                tab[i][j] = constraints[i][j];
            }
            tab[i][width-1] = constraints[i][n];        // b vector
            tab[i][n+i] = 1;                            // diagonal slack
            active[i] = n+i;                            // set active
        }
        for(uint i = 0; i < nonstandard.size(); ++i)
            tab[nonstandard[i]][n+nonstandard[i]] *= -1;

        backup_tableau();

    }

    virtual bool restore_tableau() {
        if(backup_tab.size() == 0)
            return false;
        tab = backup_tab;
        active = backup_active;
        iter = -1;
        return true;
    }


    virtual void backup_tableau() {
        backup_tab = tab;
        backup_active = active;
    }

    void load_array(std::string fname) {

        std::ifstream fp(fname.c_str());
        if( !fp.is_open() ) {
            std::cerr << "failed to open file " + fname << std::endl;
            throw;
        }

        iter = -1;
        PERFC_MEM = 0;
        PERFC_ADDMUL = 0;
        PERFC_DIV = 0;

        std::vector<T> costs;
        std::vector< std::vector<T> > constraints;
        nonstandard.clear();

        int nr=0;
        std::string line;
        bool constraint_section = false;
        while( !fp.eof() && fp.good() ) {              // read line by line

            getline(fp, line);
            if(line == "Maximize") {
                getline(fp, line);                     // get objective function
                costs = this->split_vars(line);              // split_vars is size n
            } else if(line == "Minimize") {
                std::cerr << "Minimisation not supported" << std::endl;
                throw;
                //~ getline(fp, line);
                //~ costs = split_vars(line);
                //~ for(int i = 0; i < costs.size(); ++i)
                    //~ costs[i] *= -1;
            } else if(line == "Subject To") {
                nr = 0;
                constraint_section = true;
            } else if(line == "End") {
                break;
            } else if(constraint_section == true && line != "") {
                constraints.push_back(this->split_vars(line, nr));   // split_vars is size n+1
                ++nr;   // track current constraint
            }

        }

        n = costs.size();
        m = constraints.size();
        width = n+1 + m + 1;                  // #variables + cost col + #constraints + (<= col)
        tabp = (T*)malloc( (m+1)*width * sizeof(T) );
        for(int i = 0; i < (m+1)*width; ++i) tabp[i] = 0;     // zero out the table
        active = std::vector<int>(m);

        tabp[m*width+n+m] = 1.;                               // set last row (obj)
        for(int i = 0; i < n; ++i)
            tabp[m*width+i] = -costs[i];

        for(int i = 0; i < m; ++i)  {                    // set the rest
            for(int j = 0; j < n; ++j) {
                tabp[i*width+j] = constraints[i][j];
            }
            tabp[i*width+width-1] = constraints[i][n];        // b vector
            tabp[i*width+n+i] = 1;                            // diagonal slack
            active[i] = n+i;                            // set active
        }
        for(uint i = 0; i < nonstandard.size(); ++i)
            tabp[nonstandard[i]*width+n+nonstandard[i]] *= -1;

        backup_tableau_array();

    }

    bool restore_tableau_array() {
        if(backup_tabp == NULL)
            return false;
        memcpy(tabp, backup_tabp, sizeof(T)*width*(m+1));
        active = backup_active;
        iter = -1;
        return true;
    }

    void backup_tableau_array() {
        memcpy(backup_tabp, tabp, sizeof(T)*width*(m+1));
        backup_active = active;
    }

    virtual void print() {

        dout.precision(1);
        std::cout << "Tableau for " << n << " variables and "
                  << m << " constraints:" << std::endl;
        for(uint i = 0; i < tab.size(); ++i) {
            for(int j = 0; j < width; ++j) {
                dout << "\t" << tab[i][j];
            }
            std::cout << std::endl;
        }

    }

    void print_array() {

        dout.precision(1);
        std::cout << "Tableau for " << n << " variables and "
                  << m << " constraints:" << std::endl;
        for(int i = 0; i < m+1; ++i) {
            for(int j = 0; j < width; ++j) {
                dout << "\t" << tabp[i*width+j];
            }
            std::cout << std::endl;
        }
    }


    /*
     * Return vector of solutions: 
     */
    virtual std::vector<T> solutions() {
        std::vector<T> sol(2*n+1);
        sol[0] = tab[m][width-1] / tab[m][2*n];
        for(int i = 0; i < m; ++i)
            sol[active[i]+1] = tab[i][width-1] / tab[i][active[i]];
        return sol;
    }

    std::vector<T> solutions_array() {
        std::vector<T> sol(2*n+1);
        sol[0] = tabp[m*width+width-1] / tabp[m*width+2*n];
        for(int i = 0; i < m; ++i)
            sol[active[i]+1] = tabp[i*width+width-1] / tabp[i*width+active[i]];
        return sol;
    }

    /*
     * Return kilobytes of memory used
     */
    unsigned int memusage() {
        //~ return width * m * sizeof(T) / 1000;   // only this during solve
                                                   // but we swamp the cache
                                                   // with necessary loads

        return ( width * m
               + active.size()
             )  * sizeof(T) / 1000;
    }

    int get_iter() { return iter; }
    int get_tabn() { return (m+1)*width; }

    virtual void solve() = 0;
    virtual std::string get_identifier() = 0;


};
