/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "../Simplex.hpp"

template <typename T>
class Simplex_block8x4_swap : public SimplexBase<T> {

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

    std::string get_identifier() { return "block8x4_swap"; }

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

        for(int i = 0; i < m-(m%8); i += 8) {
            PERFC_MEM+=8; PERFC_ADDMUL+=8;
            T fac0 = tabp[(i+0)*width+col] * ipiv;
            T fac1 = tabp[(i+1)*width+col] * ipiv;
            T fac2 = tabp[(i+2)*width+col] * ipiv;
            T fac3 = tabp[(i+3)*width+col] * ipiv;
            T fac4 = tabp[(i+4)*width+col] * ipiv;
            T fac5 = tabp[(i+5)*width+col] * ipiv;
            T fac6 = tabp[(i+6)*width+col] * ipiv;
            T fac7 = tabp[(i+7)*width+col] * ipiv;
            
            for(int j = 0; j < width-(width%4); j += 4) {
                T r0 = tabp[m*width+j+0];
                T r1 = tabp[m*width+j+1];
                T r2 = tabp[m*width+j+2];
                T r3 = tabp[m*width+j+3];

                //---------- i + 0 ----------
                PERFC_MEM += 4;
                T l_0_0 = tabp[(i+0)*width+j+0];
                T l_0_1 = tabp[(i+0)*width+j+1];
                T l_0_2 = tabp[(i+0)*width+j+2];
                T l_0_3 = tabp[(i+0)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_0_0 = l_0_0 - fac0*r0;
                T p_0_1 = l_0_1 - fac0*r1;
                T p_0_2 = l_0_2 - fac0*r2;
                T p_0_3 = l_0_3 - fac0*r3;

                tabp[(i+0)*width+j+0] = p_0_0;
                tabp[(i+0)*width+j+1] = p_0_1;
                tabp[(i+0)*width+j+2] = p_0_2;
                tabp[(i+0)*width+j+3] = p_0_3;

                //---------- i + 1 ----------
                PERFC_MEM += 4;
                T l_1_0 = tabp[(i+1)*width+j+0];
                T l_1_1 = tabp[(i+1)*width+j+1];
                T l_1_2 = tabp[(i+1)*width+j+2];
                T l_1_3 = tabp[(i+1)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_1_0 = l_1_0 - fac1*r0;
                T p_1_1 = l_1_1 - fac1*r1;
                T p_1_2 = l_1_2 - fac1*r2;
                T p_1_3 = l_1_3 - fac1*r3;

                tabp[(i+1)*width+j+0] = p_1_0;
                tabp[(i+1)*width+j+1] = p_1_1;
                tabp[(i+1)*width+j+2] = p_1_2;
                tabp[(i+1)*width+j+3] = p_1_3;

                //---------- i + 2 ----------
                PERFC_MEM += 4;
                T l_2_0 = tabp[(i+2)*width+j+0];
                T l_2_1 = tabp[(i+2)*width+j+1];
                T l_2_2 = tabp[(i+2)*width+j+2];
                T l_2_3 = tabp[(i+2)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_2_0 = l_2_0 - fac2*r0;
                T p_2_1 = l_2_1 - fac2*r1;
                T p_2_2 = l_2_2 - fac2*r2;
                T p_2_3 = l_2_3 - fac2*r3;

                tabp[(i+2)*width+j+0] = p_2_0;
                tabp[(i+2)*width+j+1] = p_2_1;
                tabp[(i+2)*width+j+2] = p_2_2;
                tabp[(i+2)*width+j+3] = p_2_3;

                //---------- i + 3 ----------
                PERFC_MEM += 4;
                T l_3_0 = tabp[(i+3)*width+j+0];
                T l_3_1 = tabp[(i+3)*width+j+1];
                T l_3_2 = tabp[(i+3)*width+j+2];
                T l_3_3 = tabp[(i+3)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_3_0 = l_3_0 - fac3*r0;
                T p_3_1 = l_3_1 - fac3*r1;
                T p_3_2 = l_3_2 - fac3*r2;
                T p_3_3 = l_3_3 - fac3*r3;

                tabp[(i+3)*width+j+0] = p_3_0;
                tabp[(i+3)*width+j+1] = p_3_1;
                tabp[(i+3)*width+j+2] = p_3_2;
                tabp[(i+3)*width+j+3] = p_3_3;

                //---------- i + 4 ----------
                PERFC_MEM += 4;
                T l_4_0 = tabp[(i+4)*width+j+0];
                T l_4_1 = tabp[(i+4)*width+j+1];
                T l_4_2 = tabp[(i+4)*width+j+2];
                T l_4_3 = tabp[(i+4)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_4_0 = l_4_0 - fac4*r0;
                T p_4_1 = l_4_1 - fac4*r1;
                T p_4_2 = l_4_2 - fac4*r2;
                T p_4_3 = l_4_3 - fac4*r3;

                tabp[(i+4)*width+j+0] = p_4_0;
                tabp[(i+4)*width+j+1] = p_4_1;
                tabp[(i+4)*width+j+2] = p_4_2;
                tabp[(i+4)*width+j+3] = p_4_3;

                //---------- i + 5 ----------
                PERFC_MEM += 4;
                T l_5_0 = tabp[(i+5)*width+j+0];
                T l_5_1 = tabp[(i+5)*width+j+1];
                T l_5_2 = tabp[(i+5)*width+j+2];
                T l_5_3 = tabp[(i+5)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_5_0 = l_5_0 - fac5*r0;
                T p_5_1 = l_5_1 - fac5*r1;
                T p_5_2 = l_5_2 - fac5*r2;
                T p_5_3 = l_5_3 - fac5*r3;

                tabp[(i+5)*width+j+0] = p_5_0;
                tabp[(i+5)*width+j+1] = p_5_1;
                tabp[(i+5)*width+j+2] = p_5_2;
                tabp[(i+5)*width+j+3] = p_5_3;

                //---------- i + 6 ----------
                PERFC_MEM += 4;
                T l_6_0 = tabp[(i+6)*width+j+0];
                T l_6_1 = tabp[(i+6)*width+j+1];
                T l_6_2 = tabp[(i+6)*width+j+2];
                T l_6_3 = tabp[(i+6)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_6_0 = l_6_0 - fac6*r0;
                T p_6_1 = l_6_1 - fac6*r1;
                T p_6_2 = l_6_2 - fac6*r2;
                T p_6_3 = l_6_3 - fac6*r3;

                tabp[(i+6)*width+j+0] = p_6_0;
                tabp[(i+6)*width+j+1] = p_6_1;
                tabp[(i+6)*width+j+2] = p_6_2;
                tabp[(i+6)*width+j+3] = p_6_3;

                //---------- i + 7 ----------
                PERFC_MEM += 4;
                T l_7_0 = tabp[(i+7)*width+j+0];
                T l_7_1 = tabp[(i+7)*width+j+1];
                T l_7_2 = tabp[(i+7)*width+j+2];
                T l_7_3 = tabp[(i+7)*width+j+3];

                PERFC_ADDMUL += 2*4;
                T p_7_0 = l_7_0 - fac7*r0;
                T p_7_1 = l_7_1 - fac7*r1;
                T p_7_2 = l_7_2 - fac7*r2;
                T p_7_3 = l_7_3 - fac7*r3;

                tabp[(i+7)*width+j+0] = p_7_0;
                tabp[(i+7)*width+j+1] = p_7_1;
                tabp[(i+7)*width+j+2] = p_7_2;
                tabp[(i+7)*width+j+3] = p_7_3;
            }

            for(int j = width-(width%4); j < width; ++j) {
                PERFC_MEM += 1;
                T r1 = tabp[m*width+j];

                PERFC_ADDMUL += 2*8;
                tabp[(i+0)*width+j] -= fac0*r1;
                tabp[(i+1)*width+j] -= fac1*r1;
                tabp[(i+2)*width+j] -= fac2*r1;
                tabp[(i+3)*width+j] -= fac3*r1;
                tabp[(i+4)*width+j] -= fac4*r1;
                tabp[(i+5)*width+j] -= fac5*r1;
                tabp[(i+6)*width+j] -= fac6*r1;
                tabp[(i+7)*width+j] -= fac7*r1;
            }
        }

        for(int i = m-(m%8); i < m; ++i) {
            T fac = tabp[i*width+col] * ipiv;
            for(int j = 0; j < width; ++j) {
                PERFC_ADDMUL += 2; ++PERFC_MEM;
                tabp[i*width+j] -= fac*tabp[m*width+j];
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
