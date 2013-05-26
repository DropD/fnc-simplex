/*
Assumptions:
    See base class.
*/



#pragma once

#include "Simplex.hpp"

template <typename T>
class SimplexSSA : public SimplexBase<T> {

    using SimplexBase<T>::m;
    using SimplexBase<T>::n;
    using SimplexBase<T>::width;
    using SimplexBase<T>::tabp;
    using SimplexBase<T>::nonstandard;
    using SimplexBase<T>::active;


    public:

    using SimplexBase<T>::PERFC_MEM;
    using SimplexBase<T>::PERFC_ADDMUL;
    using SimplexBase<T>::PERFC_DIV;

    std::string get_identifier() { return "ssa"; }

    void solve() {

        for(uint i = 0; i < nonstandard.size(); ++i) {       // can be ignored in cost measure, as we solve highly standard problems in lpbench
            std::cerr << "Reduction to standard LPs." << std::endl;
            throw;
        }

        while(true) {

            int col = pivot_col();        // width unit-stride memory accesses and comparisons
            if(col >= width) break;       // no negative -> finished
            int row = pivot_row(col);     // m width-stride memory accesses and comparisons, 1..m flop divisons in fairly bad stride
            if(row >= m) {
                std::cerr << "Invalid pivot row (problem infeasible)" << std::endl;
                throw;
            }

            basis_exchange(row, col);     // m*(width+1) mem access, m*(2*width+1) float ops, unit stride

            #ifdef VERBOSE
                std::cout << "Pivoting element (" << row << ", " << col << ")" << std::endl;
                this->print();
                std::cout << "active:";
                for(uint i = 0; i < active.size(); ++i)
                    std::cout << " " << active[i];
                std::cout << std::endl;
            #endif // VERBOSE
        }

    }


    protected:

    inline int pivot_col() {
        T min = 0;              // min must be negative
        int idx = width;
        for(int i = 0; i < width-1; ++i) {
            ++PERFC_MEM;
            if(tabp[m*width+i] < min) {
                min = tabp[m*width+i];
                idx = i;
            }
        }
        if(tabp[m*width+idx] > -1e-9) return width;   // prevent annihilation
        return idx;
    }

    inline int pivot_row(int col) {
        T min = 0;              // min can be arbitrary
        bool any = false;
        int idx = m;
        for(int i = 0; i < m; ++i) {
            ++PERFC_MEM;
            if(tabp[i*width+col] <= 0)
                continue;
            T ratio = tabp[i*width+width-1]/tabp[i*width+col];
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
        ++PERFC_MEM;
        T pivot = tabp[row*width+col];
        for(int i = 0; i < m+1; ++i) {
            if(i == row)
                continue;
            ++PERFC_DIV; ++PERFC_MEM;
            T fac = tabp[i*width+col]/pivot;
            PERFC_ADDMUL += 2*width; PERFC_MEM += width;
            for(int j = 0; j < width-(width%4); j += 4) {

                T l1 = tabp[i*width+j];
                T r1 = tabp[row*width+j];
                T l2 = tabp[i*width+j+1];
                T r2 = tabp[row*width+j+1];
                T l3 = tabp[i*width+j+2];
                T r3 = tabp[row*width+j+2];
                T l4 = tabp[i*width+j+3];
                T r4 = tabp[row*width+j+3];

                T p1 = l1 - fac*r1;
                T p2 = l2 - fac*r2;
                T p3 = l3 - fac*r3;
                T p4 = l4 - fac*r4;

                tabp[i*width+j] = p1;
                tabp[i*width+j+1] = p2;
                tabp[i*width+j+2] = p3;
                tabp[i*width+j+3] = p4;

            }

            for(int j = width-(width%4); j < width; ++j) {
                tabp[i*width+j] -= fac*tabp[row*width+j];
            }



        }
        active[row] = col;
    }

    void load(std::string fname) {
        this->load_array(fname);
    }

    void print() {
        this->print_array();
    }

    std::vector<T> solutions() {
        return this->solutions_array();
    }


};
