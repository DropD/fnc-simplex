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
  We have a maximisation problem
  'Maximize' comes before 'Subject To'
  No newlines between 'Maximize'/'Subject To' and the data, or within the data
  All constraint inequalities are <= and rhs is positive
  Variables determined by position (0 for empty)
  Variables are >= 0
*/



#pragma once

#include "base_general.hpp"

template <typename T>
class Simplex : public SimplexBase<T> {

    using SimplexBase<T>::m;
    using SimplexBase<T>::width;
    using SimplexBase<T>::tab;
    using SimplexBase<T>::nonstandard;
    using SimplexBase<T>::active;

    public:

    using SimplexBase<T>::PERFC_MEM;
    using SimplexBase<T>::PERFC_ADDMUL;
    using SimplexBase<T>::PERFC_DIV;

    std::string get_identifier() { return "baseline"; }

    void solve() {

        for(int i = 0; i < nonstandard.size(); ++i) {       // can be ignored in cost measure, as we solve highly standard problems in lpbench

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
                this->print();
                std::cout << "active:";
                for(int i = 0; i < active.size(); ++i)
                    std::cout << " " << active[i];
                std::cout << std::endl;
                usleep(0.2*1000*1000);
            #endif // VERBOSE

        }

        while(true) {

            int col = pivot_col();        // width unit-stride memory accesses and comparisons
            if(col >= width) break;       // no negative -> finished
            int row = pivot_row(col);     // m width-stride memory accesses and comparisons, 1-m flop divisons in fairly bad stride
            if(row >= m) {
                std::cerr << "Invalid pivot row (problem infeasible)" << std::endl;
                throw;
            }

            basis_exchange(row, col);

            #ifdef VERBOSE
                std::cout << "Pivoting element (" << row << ", " << col << ")" << std::endl;
                this->print();
                std::cout << "active:";
                for(int i = 0; i < active.size(); ++i)
                    std::cout << " " << active[i];
                std::cout << std::endl;
                usleep(0.2*1000*1000);
            #endif // VERBOSE
        }

    }


    protected:

    inline int pivot_col() {
        T min = 0;              // min must be negative
        int idx = width;
        for(int i = 0; i < width-1; ++i) {
            ++PERFC_MEM;
            if(tab[m][i] < min) {
                min = tab[m][i];
                idx = i;
            }
        }
        //~ if(tab[m][idx] > -0.0001) { dout.precision(13); dout_HERE dout,tab[m][idx]; return width; }
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

    inline int pivot_row(int col) {
        T min = 0;              // min can be arbitrary
        bool any = false;
        int idx = m;
        for(int i = 0; i < m; ++i) {
            ++PERFC_MEM;
            if(tab[i][col] <= 0)
                continue;
            T ratio = tab[i][width-1]/tab[i][col];
            ++PERFC_DIV;
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
        for(int i = 0; i < m+1; ++i) {            // m-1 iterations
            if(i == row)
                continue;
            ++PERFC_DIV;
            T fac = tab[i][col]/pivot;
            for(int j = 0; j < width; ++j) {      // 2*width addmul in unit stride
                PERFC_ADDMUL += 2;
                tab[i][j] -= fac*tab[row][j];
            }
        }
        active[row] = col;                        // 1 random memory access
    }

};
