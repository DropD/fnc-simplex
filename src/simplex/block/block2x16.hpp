/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "../Simplex.hpp"

template <typename T>
class Simplex_block2x16 : public SimplexBase<T> {

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

    std::string get_identifier() { return "block2x16"; }

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
        
        ++PERFC_DIV;
        T ipiv = 1. / pivot;

        for(int i = 0; i < m-(m%2); i += 2) {
            PERFC_MEM+=2; PERFC_ADDMUL+=2;
            T fac0 = tabp[(i+0)*width+col] * ipiv;
            T fac1 = tabp[(i+1)*width+col] * ipiv;

            PERFC_ADDMUL += 2*2 * width;
            PERFC_MEM += 2*width;

            for(int j = 0; j < width-(width%16); j += 16) {
                T r0 = tabp[row*width+j+0];
                T r1 = tabp[row*width+j+1];
                T r2 = tabp[row*width+j+2];
                T r3 = tabp[row*width+j+3];
                T r4 = tabp[row*width+j+4];
                T r5 = tabp[row*width+j+5];
                T r6 = tabp[row*width+j+6];
                T r7 = tabp[row*width+j+7];
                T r8 = tabp[row*width+j+8];
                T r9 = tabp[row*width+j+9];
                T r10 = tabp[row*width+j+10];
                T r11 = tabp[row*width+j+11];
                T r12 = tabp[row*width+j+12];
                T r13 = tabp[row*width+j+13];
                T r14 = tabp[row*width+j+14];
                T r15 = tabp[row*width+j+15];

                //---------- i + 0 ----------
                if(i+0 != row) {
                    T l_0_0 = tabp[(i+0)*width+j+0];
                    T l_0_1 = tabp[(i+0)*width+j+1];
                    T l_0_2 = tabp[(i+0)*width+j+2];
                    T l_0_3 = tabp[(i+0)*width+j+3];
                    T l_0_4 = tabp[(i+0)*width+j+4];
                    T l_0_5 = tabp[(i+0)*width+j+5];
                    T l_0_6 = tabp[(i+0)*width+j+6];
                    T l_0_7 = tabp[(i+0)*width+j+7];
                    T l_0_8 = tabp[(i+0)*width+j+8];
                    T l_0_9 = tabp[(i+0)*width+j+9];
                    T l_0_10 = tabp[(i+0)*width+j+10];
                    T l_0_11 = tabp[(i+0)*width+j+11];
                    T l_0_12 = tabp[(i+0)*width+j+12];
                    T l_0_13 = tabp[(i+0)*width+j+13];
                    T l_0_14 = tabp[(i+0)*width+j+14];
                    T l_0_15 = tabp[(i+0)*width+j+15];

                    T p_0_0 = l_0_0 - fac0*r0;
                    T p_0_1 = l_0_1 - fac0*r1;
                    T p_0_2 = l_0_2 - fac0*r2;
                    T p_0_3 = l_0_3 - fac0*r3;
                    T p_0_4 = l_0_4 - fac0*r4;
                    T p_0_5 = l_0_5 - fac0*r5;
                    T p_0_6 = l_0_6 - fac0*r6;
                    T p_0_7 = l_0_7 - fac0*r7;
                    T p_0_8 = l_0_8 - fac0*r8;
                    T p_0_9 = l_0_9 - fac0*r9;
                    T p_0_10 = l_0_10 - fac0*r10;
                    T p_0_11 = l_0_11 - fac0*r11;
                    T p_0_12 = l_0_12 - fac0*r12;
                    T p_0_13 = l_0_13 - fac0*r13;
                    T p_0_14 = l_0_14 - fac0*r14;
                    T p_0_15 = l_0_15 - fac0*r15;

                    tabp[(i+0)*width+j+0] = p_0_0;
                    tabp[(i+0)*width+j+1] = p_0_1;
                    tabp[(i+0)*width+j+2] = p_0_2;
                    tabp[(i+0)*width+j+3] = p_0_3;
                    tabp[(i+0)*width+j+4] = p_0_4;
                    tabp[(i+0)*width+j+5] = p_0_5;
                    tabp[(i+0)*width+j+6] = p_0_6;
                    tabp[(i+0)*width+j+7] = p_0_7;
                    tabp[(i+0)*width+j+8] = p_0_8;
                    tabp[(i+0)*width+j+9] = p_0_9;
                    tabp[(i+0)*width+j+10] = p_0_10;
                    tabp[(i+0)*width+j+11] = p_0_11;
                    tabp[(i+0)*width+j+12] = p_0_12;
                    tabp[(i+0)*width+j+13] = p_0_13;
                    tabp[(i+0)*width+j+14] = p_0_14;
                    tabp[(i+0)*width+j+15] = p_0_15;
                }

                //---------- i + 1 ----------
                if(i+1 != row) {
                    T l_1_0 = tabp[(i+1)*width+j+0];
                    T l_1_1 = tabp[(i+1)*width+j+1];
                    T l_1_2 = tabp[(i+1)*width+j+2];
                    T l_1_3 = tabp[(i+1)*width+j+3];
                    T l_1_4 = tabp[(i+1)*width+j+4];
                    T l_1_5 = tabp[(i+1)*width+j+5];
                    T l_1_6 = tabp[(i+1)*width+j+6];
                    T l_1_7 = tabp[(i+1)*width+j+7];
                    T l_1_8 = tabp[(i+1)*width+j+8];
                    T l_1_9 = tabp[(i+1)*width+j+9];
                    T l_1_10 = tabp[(i+1)*width+j+10];
                    T l_1_11 = tabp[(i+1)*width+j+11];
                    T l_1_12 = tabp[(i+1)*width+j+12];
                    T l_1_13 = tabp[(i+1)*width+j+13];
                    T l_1_14 = tabp[(i+1)*width+j+14];
                    T l_1_15 = tabp[(i+1)*width+j+15];

                    T p_1_0 = l_1_0 - fac1*r0;
                    T p_1_1 = l_1_1 - fac1*r1;
                    T p_1_2 = l_1_2 - fac1*r2;
                    T p_1_3 = l_1_3 - fac1*r3;
                    T p_1_4 = l_1_4 - fac1*r4;
                    T p_1_5 = l_1_5 - fac1*r5;
                    T p_1_6 = l_1_6 - fac1*r6;
                    T p_1_7 = l_1_7 - fac1*r7;
                    T p_1_8 = l_1_8 - fac1*r8;
                    T p_1_9 = l_1_9 - fac1*r9;
                    T p_1_10 = l_1_10 - fac1*r10;
                    T p_1_11 = l_1_11 - fac1*r11;
                    T p_1_12 = l_1_12 - fac1*r12;
                    T p_1_13 = l_1_13 - fac1*r13;
                    T p_1_14 = l_1_14 - fac1*r14;
                    T p_1_15 = l_1_15 - fac1*r15;

                    tabp[(i+1)*width+j+0] = p_1_0;
                    tabp[(i+1)*width+j+1] = p_1_1;
                    tabp[(i+1)*width+j+2] = p_1_2;
                    tabp[(i+1)*width+j+3] = p_1_3;
                    tabp[(i+1)*width+j+4] = p_1_4;
                    tabp[(i+1)*width+j+5] = p_1_5;
                    tabp[(i+1)*width+j+6] = p_1_6;
                    tabp[(i+1)*width+j+7] = p_1_7;
                    tabp[(i+1)*width+j+8] = p_1_8;
                    tabp[(i+1)*width+j+9] = p_1_9;
                    tabp[(i+1)*width+j+10] = p_1_10;
                    tabp[(i+1)*width+j+11] = p_1_11;
                    tabp[(i+1)*width+j+12] = p_1_12;
                    tabp[(i+1)*width+j+13] = p_1_13;
                    tabp[(i+1)*width+j+14] = p_1_14;
                    tabp[(i+1)*width+j+15] = p_1_15;
                }
            }

            for(int j = width-(width%16); j < width; ++j) {
                T r1 = tabp[row*width+j];

                if(i+0 != row) {
                    tabp[(i+0)*width+j] -= fac0*r1;
                }
                if(i+1 != row) {
                    tabp[(i+1)*width+j] -= fac1*r1;
                }
            }
        }

        for(int i = m-(m%2); i < m+1; ++i) {
            T fac = tabp[i*width+col] * ipiv;
            for(int j = 0; j < width; ++j) {
                PERFC_ADDMUL += 2; ++PERFC_MEM;
                tabp[i*width+j] -= fac*tabp[row*width+j];
            }
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
