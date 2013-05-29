/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "Simplex.hpp"

template <typename T>
class SimplexBlock_swap : public SimplexBase<T> {

    using SimplexBase<T>::m;
    using SimplexBase<T>::n;
    using SimplexBase<T>::width;
    using SimplexBase<T>::tabp;
    using SimplexBase<T>::nonstandard;
    using SimplexBase<T>::active;
    using SimplexBase<T>::iter;


    public:

    using SimplexBase<T>::PERFC_MEM;
    using SimplexBase<T>::PERFC_ADDMUL;
    using SimplexBase<T>::PERFC_DIV;

    std::string get_identifier() { return "block_swap"; }

    void solve() {

        for(uint i = 0; i < nonstandard.size(); ++i) {       // can be ignored in cost measure, as we solve highly standard problems in lpbench
            std::cerr << "Reduction to standard LPs." << std::endl;
            throw;
        }

        while(true) {

            ++iter;

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
        // swap pivot row and last
        PERFC_MEM += 2 * width;
        for(int i = 0; i < width; ++i) {
            T tmp = tabp[row*width+i];
            tabp[row*width+i] = tabp[m*width+i];
            tabp[m*width+i] = tmp;
        }
        ++PERFC_DIV;
        T ipiv = 1. / pivot;
        for(int i = 0; i < m; i += 2) {     // peel off m+1, as we expect an even amount of dofs
            PERFC_MEM+=2; PERFC_ADDMUL+=2;
            T fac1 = tabp[i*width+col] * ipiv;
            T fac2 = tabp[(i+1)*width+col] * ipiv;

            for(int j = 0; j < width-(width%4); j += 4) {

                //~ PERFC_MEM += 4;  // not from memory, as width*sizeof(T) ~ 1-24 kB should easily fit into L2/L3
                T r1 = tabp[m*width+j];
                T r2 = tabp[m*width+j+1];
                T r3 = tabp[m*width+j+2];
                T r4 = tabp[m*width+j+3];

                PERFC_MEM += 4;
                T la1 = tabp[i*width+j];
                T la2 = tabp[i*width+j+1];
                T la3 = tabp[i*width+j+2];
                T la4 = tabp[i*width+j+3];

                PERFC_ADDMUL += 8;
                T pa1 = la1 - fac1*r1;
                T pa2 = la2 - fac1*r2;
                T pa3 = la3 - fac1*r3;
                T pa4 = la4 - fac1*r4;

                tabp[i*width+j]   = pa1;
                tabp[i*width+j+1] = pa2;
                tabp[i*width+j+2] = pa3;
                tabp[i*width+j+3] = pa4;

                PERFC_MEM += 4;
                T lb1 = tabp[(i+1)*width+j];
                T lb2 = tabp[(i+1)*width+j+1];
                T lb3 = tabp[(i+1)*width+j+2];
                T lb4 = tabp[(i+1)*width+j+3];

                PERFC_ADDMUL += 8;
                T pb1 = lb1 - fac2*r1;
                T pb2 = lb2 - fac2*r2;
                T pb3 = lb3 - fac2*r3;
                T pb4 = lb4 - fac2*r4;

                tabp[(i+1)*width+j]   = pb1;
                tabp[(i+1)*width+j+1] = pb2;
                tabp[(i+1)*width+j+2] = pb3;
                tabp[(i+1)*width+j+3] = pb4;
            }

            for(int j = width-(width%4); j < width; ++j) {
                PERFC_ADDMUL+=4; PERFC_MEM+=2;
                tabp[i*width+j] -= fac1*tabp[m*width+j];
                tabp[(i+1)*width+j] -= fac2*tabp[m*width+j];
            }

        }
        // swap back
        PERFC_MEM += 2 * width;
        for(int i = 0; i < width; ++i) {
            T tmp = tabp[row*width+i];
            tabp[row*width+i] = tabp[m*width+i];
            tabp[m*width+i] = tmp;
        }

        active[row] = col;
    }

    void load(std::string fname) {
        this->load_array(fname);
    }

    virtual bool restore_tableau() {
        return this->restore_tableau_array();
    }

    void print() {
        this->print_array();
    }

    std::vector<T> solutions() {
        return this->solutions_array();
    }


};
