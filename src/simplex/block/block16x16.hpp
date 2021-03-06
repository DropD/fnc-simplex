/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "../Simplex.hpp"

template <typename T>
class Simplex_block16x16 : public SimplexBase<T> {

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

    std::string get_identifier() { return "block16x16"; }

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

        for(int i = 0; i < m-(m%16); i += 16) {
            PERFC_MEM+=16; PERFC_ADDMUL+=16;
            T fac0 = tabp[(i+0)*width+col] * ipiv;
            T fac1 = tabp[(i+1)*width+col] * ipiv;
            T fac2 = tabp[(i+2)*width+col] * ipiv;
            T fac3 = tabp[(i+3)*width+col] * ipiv;
            T fac4 = tabp[(i+4)*width+col] * ipiv;
            T fac5 = tabp[(i+5)*width+col] * ipiv;
            T fac6 = tabp[(i+6)*width+col] * ipiv;
            T fac7 = tabp[(i+7)*width+col] * ipiv;
            T fac8 = tabp[(i+8)*width+col] * ipiv;
            T fac9 = tabp[(i+9)*width+col] * ipiv;
            T fac10 = tabp[(i+10)*width+col] * ipiv;
            T fac11 = tabp[(i+11)*width+col] * ipiv;
            T fac12 = tabp[(i+12)*width+col] * ipiv;
            T fac13 = tabp[(i+13)*width+col] * ipiv;
            T fac14 = tabp[(i+14)*width+col] * ipiv;
            T fac15 = tabp[(i+15)*width+col] * ipiv;

            PERFC_ADDMUL += 2*16 * width;
            PERFC_MEM += 16*width;

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

                //---------- i + 2 ----------
                if(i+2 != row) {
                    T l_2_0 = tabp[(i+2)*width+j+0];
                    T l_2_1 = tabp[(i+2)*width+j+1];
                    T l_2_2 = tabp[(i+2)*width+j+2];
                    T l_2_3 = tabp[(i+2)*width+j+3];
                    T l_2_4 = tabp[(i+2)*width+j+4];
                    T l_2_5 = tabp[(i+2)*width+j+5];
                    T l_2_6 = tabp[(i+2)*width+j+6];
                    T l_2_7 = tabp[(i+2)*width+j+7];
                    T l_2_8 = tabp[(i+2)*width+j+8];
                    T l_2_9 = tabp[(i+2)*width+j+9];
                    T l_2_10 = tabp[(i+2)*width+j+10];
                    T l_2_11 = tabp[(i+2)*width+j+11];
                    T l_2_12 = tabp[(i+2)*width+j+12];
                    T l_2_13 = tabp[(i+2)*width+j+13];
                    T l_2_14 = tabp[(i+2)*width+j+14];
                    T l_2_15 = tabp[(i+2)*width+j+15];

                    T p_2_0 = l_2_0 - fac2*r0;
                    T p_2_1 = l_2_1 - fac2*r1;
                    T p_2_2 = l_2_2 - fac2*r2;
                    T p_2_3 = l_2_3 - fac2*r3;
                    T p_2_4 = l_2_4 - fac2*r4;
                    T p_2_5 = l_2_5 - fac2*r5;
                    T p_2_6 = l_2_6 - fac2*r6;
                    T p_2_7 = l_2_7 - fac2*r7;
                    T p_2_8 = l_2_8 - fac2*r8;
                    T p_2_9 = l_2_9 - fac2*r9;
                    T p_2_10 = l_2_10 - fac2*r10;
                    T p_2_11 = l_2_11 - fac2*r11;
                    T p_2_12 = l_2_12 - fac2*r12;
                    T p_2_13 = l_2_13 - fac2*r13;
                    T p_2_14 = l_2_14 - fac2*r14;
                    T p_2_15 = l_2_15 - fac2*r15;

                    tabp[(i+2)*width+j+0] = p_2_0;
                    tabp[(i+2)*width+j+1] = p_2_1;
                    tabp[(i+2)*width+j+2] = p_2_2;
                    tabp[(i+2)*width+j+3] = p_2_3;
                    tabp[(i+2)*width+j+4] = p_2_4;
                    tabp[(i+2)*width+j+5] = p_2_5;
                    tabp[(i+2)*width+j+6] = p_2_6;
                    tabp[(i+2)*width+j+7] = p_2_7;
                    tabp[(i+2)*width+j+8] = p_2_8;
                    tabp[(i+2)*width+j+9] = p_2_9;
                    tabp[(i+2)*width+j+10] = p_2_10;
                    tabp[(i+2)*width+j+11] = p_2_11;
                    tabp[(i+2)*width+j+12] = p_2_12;
                    tabp[(i+2)*width+j+13] = p_2_13;
                    tabp[(i+2)*width+j+14] = p_2_14;
                    tabp[(i+2)*width+j+15] = p_2_15;
                }

                //---------- i + 3 ----------
                if(i+3 != row) {
                    T l_3_0 = tabp[(i+3)*width+j+0];
                    T l_3_1 = tabp[(i+3)*width+j+1];
                    T l_3_2 = tabp[(i+3)*width+j+2];
                    T l_3_3 = tabp[(i+3)*width+j+3];
                    T l_3_4 = tabp[(i+3)*width+j+4];
                    T l_3_5 = tabp[(i+3)*width+j+5];
                    T l_3_6 = tabp[(i+3)*width+j+6];
                    T l_3_7 = tabp[(i+3)*width+j+7];
                    T l_3_8 = tabp[(i+3)*width+j+8];
                    T l_3_9 = tabp[(i+3)*width+j+9];
                    T l_3_10 = tabp[(i+3)*width+j+10];
                    T l_3_11 = tabp[(i+3)*width+j+11];
                    T l_3_12 = tabp[(i+3)*width+j+12];
                    T l_3_13 = tabp[(i+3)*width+j+13];
                    T l_3_14 = tabp[(i+3)*width+j+14];
                    T l_3_15 = tabp[(i+3)*width+j+15];

                    T p_3_0 = l_3_0 - fac3*r0;
                    T p_3_1 = l_3_1 - fac3*r1;
                    T p_3_2 = l_3_2 - fac3*r2;
                    T p_3_3 = l_3_3 - fac3*r3;
                    T p_3_4 = l_3_4 - fac3*r4;
                    T p_3_5 = l_3_5 - fac3*r5;
                    T p_3_6 = l_3_6 - fac3*r6;
                    T p_3_7 = l_3_7 - fac3*r7;
                    T p_3_8 = l_3_8 - fac3*r8;
                    T p_3_9 = l_3_9 - fac3*r9;
                    T p_3_10 = l_3_10 - fac3*r10;
                    T p_3_11 = l_3_11 - fac3*r11;
                    T p_3_12 = l_3_12 - fac3*r12;
                    T p_3_13 = l_3_13 - fac3*r13;
                    T p_3_14 = l_3_14 - fac3*r14;
                    T p_3_15 = l_3_15 - fac3*r15;

                    tabp[(i+3)*width+j+0] = p_3_0;
                    tabp[(i+3)*width+j+1] = p_3_1;
                    tabp[(i+3)*width+j+2] = p_3_2;
                    tabp[(i+3)*width+j+3] = p_3_3;
                    tabp[(i+3)*width+j+4] = p_3_4;
                    tabp[(i+3)*width+j+5] = p_3_5;
                    tabp[(i+3)*width+j+6] = p_3_6;
                    tabp[(i+3)*width+j+7] = p_3_7;
                    tabp[(i+3)*width+j+8] = p_3_8;
                    tabp[(i+3)*width+j+9] = p_3_9;
                    tabp[(i+3)*width+j+10] = p_3_10;
                    tabp[(i+3)*width+j+11] = p_3_11;
                    tabp[(i+3)*width+j+12] = p_3_12;
                    tabp[(i+3)*width+j+13] = p_3_13;
                    tabp[(i+3)*width+j+14] = p_3_14;
                    tabp[(i+3)*width+j+15] = p_3_15;
                }

                //---------- i + 4 ----------
                if(i+4 != row) {
                    T l_4_0 = tabp[(i+4)*width+j+0];
                    T l_4_1 = tabp[(i+4)*width+j+1];
                    T l_4_2 = tabp[(i+4)*width+j+2];
                    T l_4_3 = tabp[(i+4)*width+j+3];
                    T l_4_4 = tabp[(i+4)*width+j+4];
                    T l_4_5 = tabp[(i+4)*width+j+5];
                    T l_4_6 = tabp[(i+4)*width+j+6];
                    T l_4_7 = tabp[(i+4)*width+j+7];
                    T l_4_8 = tabp[(i+4)*width+j+8];
                    T l_4_9 = tabp[(i+4)*width+j+9];
                    T l_4_10 = tabp[(i+4)*width+j+10];
                    T l_4_11 = tabp[(i+4)*width+j+11];
                    T l_4_12 = tabp[(i+4)*width+j+12];
                    T l_4_13 = tabp[(i+4)*width+j+13];
                    T l_4_14 = tabp[(i+4)*width+j+14];
                    T l_4_15 = tabp[(i+4)*width+j+15];

                    T p_4_0 = l_4_0 - fac4*r0;
                    T p_4_1 = l_4_1 - fac4*r1;
                    T p_4_2 = l_4_2 - fac4*r2;
                    T p_4_3 = l_4_3 - fac4*r3;
                    T p_4_4 = l_4_4 - fac4*r4;
                    T p_4_5 = l_4_5 - fac4*r5;
                    T p_4_6 = l_4_6 - fac4*r6;
                    T p_4_7 = l_4_7 - fac4*r7;
                    T p_4_8 = l_4_8 - fac4*r8;
                    T p_4_9 = l_4_9 - fac4*r9;
                    T p_4_10 = l_4_10 - fac4*r10;
                    T p_4_11 = l_4_11 - fac4*r11;
                    T p_4_12 = l_4_12 - fac4*r12;
                    T p_4_13 = l_4_13 - fac4*r13;
                    T p_4_14 = l_4_14 - fac4*r14;
                    T p_4_15 = l_4_15 - fac4*r15;

                    tabp[(i+4)*width+j+0] = p_4_0;
                    tabp[(i+4)*width+j+1] = p_4_1;
                    tabp[(i+4)*width+j+2] = p_4_2;
                    tabp[(i+4)*width+j+3] = p_4_3;
                    tabp[(i+4)*width+j+4] = p_4_4;
                    tabp[(i+4)*width+j+5] = p_4_5;
                    tabp[(i+4)*width+j+6] = p_4_6;
                    tabp[(i+4)*width+j+7] = p_4_7;
                    tabp[(i+4)*width+j+8] = p_4_8;
                    tabp[(i+4)*width+j+9] = p_4_9;
                    tabp[(i+4)*width+j+10] = p_4_10;
                    tabp[(i+4)*width+j+11] = p_4_11;
                    tabp[(i+4)*width+j+12] = p_4_12;
                    tabp[(i+4)*width+j+13] = p_4_13;
                    tabp[(i+4)*width+j+14] = p_4_14;
                    tabp[(i+4)*width+j+15] = p_4_15;
                }

                //---------- i + 5 ----------
                if(i+5 != row) {
                    T l_5_0 = tabp[(i+5)*width+j+0];
                    T l_5_1 = tabp[(i+5)*width+j+1];
                    T l_5_2 = tabp[(i+5)*width+j+2];
                    T l_5_3 = tabp[(i+5)*width+j+3];
                    T l_5_4 = tabp[(i+5)*width+j+4];
                    T l_5_5 = tabp[(i+5)*width+j+5];
                    T l_5_6 = tabp[(i+5)*width+j+6];
                    T l_5_7 = tabp[(i+5)*width+j+7];
                    T l_5_8 = tabp[(i+5)*width+j+8];
                    T l_5_9 = tabp[(i+5)*width+j+9];
                    T l_5_10 = tabp[(i+5)*width+j+10];
                    T l_5_11 = tabp[(i+5)*width+j+11];
                    T l_5_12 = tabp[(i+5)*width+j+12];
                    T l_5_13 = tabp[(i+5)*width+j+13];
                    T l_5_14 = tabp[(i+5)*width+j+14];
                    T l_5_15 = tabp[(i+5)*width+j+15];

                    T p_5_0 = l_5_0 - fac5*r0;
                    T p_5_1 = l_5_1 - fac5*r1;
                    T p_5_2 = l_5_2 - fac5*r2;
                    T p_5_3 = l_5_3 - fac5*r3;
                    T p_5_4 = l_5_4 - fac5*r4;
                    T p_5_5 = l_5_5 - fac5*r5;
                    T p_5_6 = l_5_6 - fac5*r6;
                    T p_5_7 = l_5_7 - fac5*r7;
                    T p_5_8 = l_5_8 - fac5*r8;
                    T p_5_9 = l_5_9 - fac5*r9;
                    T p_5_10 = l_5_10 - fac5*r10;
                    T p_5_11 = l_5_11 - fac5*r11;
                    T p_5_12 = l_5_12 - fac5*r12;
                    T p_5_13 = l_5_13 - fac5*r13;
                    T p_5_14 = l_5_14 - fac5*r14;
                    T p_5_15 = l_5_15 - fac5*r15;

                    tabp[(i+5)*width+j+0] = p_5_0;
                    tabp[(i+5)*width+j+1] = p_5_1;
                    tabp[(i+5)*width+j+2] = p_5_2;
                    tabp[(i+5)*width+j+3] = p_5_3;
                    tabp[(i+5)*width+j+4] = p_5_4;
                    tabp[(i+5)*width+j+5] = p_5_5;
                    tabp[(i+5)*width+j+6] = p_5_6;
                    tabp[(i+5)*width+j+7] = p_5_7;
                    tabp[(i+5)*width+j+8] = p_5_8;
                    tabp[(i+5)*width+j+9] = p_5_9;
                    tabp[(i+5)*width+j+10] = p_5_10;
                    tabp[(i+5)*width+j+11] = p_5_11;
                    tabp[(i+5)*width+j+12] = p_5_12;
                    tabp[(i+5)*width+j+13] = p_5_13;
                    tabp[(i+5)*width+j+14] = p_5_14;
                    tabp[(i+5)*width+j+15] = p_5_15;
                }

                //---------- i + 6 ----------
                if(i+6 != row) {
                    T l_6_0 = tabp[(i+6)*width+j+0];
                    T l_6_1 = tabp[(i+6)*width+j+1];
                    T l_6_2 = tabp[(i+6)*width+j+2];
                    T l_6_3 = tabp[(i+6)*width+j+3];
                    T l_6_4 = tabp[(i+6)*width+j+4];
                    T l_6_5 = tabp[(i+6)*width+j+5];
                    T l_6_6 = tabp[(i+6)*width+j+6];
                    T l_6_7 = tabp[(i+6)*width+j+7];
                    T l_6_8 = tabp[(i+6)*width+j+8];
                    T l_6_9 = tabp[(i+6)*width+j+9];
                    T l_6_10 = tabp[(i+6)*width+j+10];
                    T l_6_11 = tabp[(i+6)*width+j+11];
                    T l_6_12 = tabp[(i+6)*width+j+12];
                    T l_6_13 = tabp[(i+6)*width+j+13];
                    T l_6_14 = tabp[(i+6)*width+j+14];
                    T l_6_15 = tabp[(i+6)*width+j+15];

                    T p_6_0 = l_6_0 - fac6*r0;
                    T p_6_1 = l_6_1 - fac6*r1;
                    T p_6_2 = l_6_2 - fac6*r2;
                    T p_6_3 = l_6_3 - fac6*r3;
                    T p_6_4 = l_6_4 - fac6*r4;
                    T p_6_5 = l_6_5 - fac6*r5;
                    T p_6_6 = l_6_6 - fac6*r6;
                    T p_6_7 = l_6_7 - fac6*r7;
                    T p_6_8 = l_6_8 - fac6*r8;
                    T p_6_9 = l_6_9 - fac6*r9;
                    T p_6_10 = l_6_10 - fac6*r10;
                    T p_6_11 = l_6_11 - fac6*r11;
                    T p_6_12 = l_6_12 - fac6*r12;
                    T p_6_13 = l_6_13 - fac6*r13;
                    T p_6_14 = l_6_14 - fac6*r14;
                    T p_6_15 = l_6_15 - fac6*r15;

                    tabp[(i+6)*width+j+0] = p_6_0;
                    tabp[(i+6)*width+j+1] = p_6_1;
                    tabp[(i+6)*width+j+2] = p_6_2;
                    tabp[(i+6)*width+j+3] = p_6_3;
                    tabp[(i+6)*width+j+4] = p_6_4;
                    tabp[(i+6)*width+j+5] = p_6_5;
                    tabp[(i+6)*width+j+6] = p_6_6;
                    tabp[(i+6)*width+j+7] = p_6_7;
                    tabp[(i+6)*width+j+8] = p_6_8;
                    tabp[(i+6)*width+j+9] = p_6_9;
                    tabp[(i+6)*width+j+10] = p_6_10;
                    tabp[(i+6)*width+j+11] = p_6_11;
                    tabp[(i+6)*width+j+12] = p_6_12;
                    tabp[(i+6)*width+j+13] = p_6_13;
                    tabp[(i+6)*width+j+14] = p_6_14;
                    tabp[(i+6)*width+j+15] = p_6_15;
                }

                //---------- i + 7 ----------
                if(i+7 != row) {
                    T l_7_0 = tabp[(i+7)*width+j+0];
                    T l_7_1 = tabp[(i+7)*width+j+1];
                    T l_7_2 = tabp[(i+7)*width+j+2];
                    T l_7_3 = tabp[(i+7)*width+j+3];
                    T l_7_4 = tabp[(i+7)*width+j+4];
                    T l_7_5 = tabp[(i+7)*width+j+5];
                    T l_7_6 = tabp[(i+7)*width+j+6];
                    T l_7_7 = tabp[(i+7)*width+j+7];
                    T l_7_8 = tabp[(i+7)*width+j+8];
                    T l_7_9 = tabp[(i+7)*width+j+9];
                    T l_7_10 = tabp[(i+7)*width+j+10];
                    T l_7_11 = tabp[(i+7)*width+j+11];
                    T l_7_12 = tabp[(i+7)*width+j+12];
                    T l_7_13 = tabp[(i+7)*width+j+13];
                    T l_7_14 = tabp[(i+7)*width+j+14];
                    T l_7_15 = tabp[(i+7)*width+j+15];

                    T p_7_0 = l_7_0 - fac7*r0;
                    T p_7_1 = l_7_1 - fac7*r1;
                    T p_7_2 = l_7_2 - fac7*r2;
                    T p_7_3 = l_7_3 - fac7*r3;
                    T p_7_4 = l_7_4 - fac7*r4;
                    T p_7_5 = l_7_5 - fac7*r5;
                    T p_7_6 = l_7_6 - fac7*r6;
                    T p_7_7 = l_7_7 - fac7*r7;
                    T p_7_8 = l_7_8 - fac7*r8;
                    T p_7_9 = l_7_9 - fac7*r9;
                    T p_7_10 = l_7_10 - fac7*r10;
                    T p_7_11 = l_7_11 - fac7*r11;
                    T p_7_12 = l_7_12 - fac7*r12;
                    T p_7_13 = l_7_13 - fac7*r13;
                    T p_7_14 = l_7_14 - fac7*r14;
                    T p_7_15 = l_7_15 - fac7*r15;

                    tabp[(i+7)*width+j+0] = p_7_0;
                    tabp[(i+7)*width+j+1] = p_7_1;
                    tabp[(i+7)*width+j+2] = p_7_2;
                    tabp[(i+7)*width+j+3] = p_7_3;
                    tabp[(i+7)*width+j+4] = p_7_4;
                    tabp[(i+7)*width+j+5] = p_7_5;
                    tabp[(i+7)*width+j+6] = p_7_6;
                    tabp[(i+7)*width+j+7] = p_7_7;
                    tabp[(i+7)*width+j+8] = p_7_8;
                    tabp[(i+7)*width+j+9] = p_7_9;
                    tabp[(i+7)*width+j+10] = p_7_10;
                    tabp[(i+7)*width+j+11] = p_7_11;
                    tabp[(i+7)*width+j+12] = p_7_12;
                    tabp[(i+7)*width+j+13] = p_7_13;
                    tabp[(i+7)*width+j+14] = p_7_14;
                    tabp[(i+7)*width+j+15] = p_7_15;
                }

                //---------- i + 8 ----------
                if(i+8 != row) {
                    T l_8_0 = tabp[(i+8)*width+j+0];
                    T l_8_1 = tabp[(i+8)*width+j+1];
                    T l_8_2 = tabp[(i+8)*width+j+2];
                    T l_8_3 = tabp[(i+8)*width+j+3];
                    T l_8_4 = tabp[(i+8)*width+j+4];
                    T l_8_5 = tabp[(i+8)*width+j+5];
                    T l_8_6 = tabp[(i+8)*width+j+6];
                    T l_8_7 = tabp[(i+8)*width+j+7];
                    T l_8_8 = tabp[(i+8)*width+j+8];
                    T l_8_9 = tabp[(i+8)*width+j+9];
                    T l_8_10 = tabp[(i+8)*width+j+10];
                    T l_8_11 = tabp[(i+8)*width+j+11];
                    T l_8_12 = tabp[(i+8)*width+j+12];
                    T l_8_13 = tabp[(i+8)*width+j+13];
                    T l_8_14 = tabp[(i+8)*width+j+14];
                    T l_8_15 = tabp[(i+8)*width+j+15];

                    T p_8_0 = l_8_0 - fac8*r0;
                    T p_8_1 = l_8_1 - fac8*r1;
                    T p_8_2 = l_8_2 - fac8*r2;
                    T p_8_3 = l_8_3 - fac8*r3;
                    T p_8_4 = l_8_4 - fac8*r4;
                    T p_8_5 = l_8_5 - fac8*r5;
                    T p_8_6 = l_8_6 - fac8*r6;
                    T p_8_7 = l_8_7 - fac8*r7;
                    T p_8_8 = l_8_8 - fac8*r8;
                    T p_8_9 = l_8_9 - fac8*r9;
                    T p_8_10 = l_8_10 - fac8*r10;
                    T p_8_11 = l_8_11 - fac8*r11;
                    T p_8_12 = l_8_12 - fac8*r12;
                    T p_8_13 = l_8_13 - fac8*r13;
                    T p_8_14 = l_8_14 - fac8*r14;
                    T p_8_15 = l_8_15 - fac8*r15;

                    tabp[(i+8)*width+j+0] = p_8_0;
                    tabp[(i+8)*width+j+1] = p_8_1;
                    tabp[(i+8)*width+j+2] = p_8_2;
                    tabp[(i+8)*width+j+3] = p_8_3;
                    tabp[(i+8)*width+j+4] = p_8_4;
                    tabp[(i+8)*width+j+5] = p_8_5;
                    tabp[(i+8)*width+j+6] = p_8_6;
                    tabp[(i+8)*width+j+7] = p_8_7;
                    tabp[(i+8)*width+j+8] = p_8_8;
                    tabp[(i+8)*width+j+9] = p_8_9;
                    tabp[(i+8)*width+j+10] = p_8_10;
                    tabp[(i+8)*width+j+11] = p_8_11;
                    tabp[(i+8)*width+j+12] = p_8_12;
                    tabp[(i+8)*width+j+13] = p_8_13;
                    tabp[(i+8)*width+j+14] = p_8_14;
                    tabp[(i+8)*width+j+15] = p_8_15;
                }

                //---------- i + 9 ----------
                if(i+9 != row) {
                    T l_9_0 = tabp[(i+9)*width+j+0];
                    T l_9_1 = tabp[(i+9)*width+j+1];
                    T l_9_2 = tabp[(i+9)*width+j+2];
                    T l_9_3 = tabp[(i+9)*width+j+3];
                    T l_9_4 = tabp[(i+9)*width+j+4];
                    T l_9_5 = tabp[(i+9)*width+j+5];
                    T l_9_6 = tabp[(i+9)*width+j+6];
                    T l_9_7 = tabp[(i+9)*width+j+7];
                    T l_9_8 = tabp[(i+9)*width+j+8];
                    T l_9_9 = tabp[(i+9)*width+j+9];
                    T l_9_10 = tabp[(i+9)*width+j+10];
                    T l_9_11 = tabp[(i+9)*width+j+11];
                    T l_9_12 = tabp[(i+9)*width+j+12];
                    T l_9_13 = tabp[(i+9)*width+j+13];
                    T l_9_14 = tabp[(i+9)*width+j+14];
                    T l_9_15 = tabp[(i+9)*width+j+15];

                    T p_9_0 = l_9_0 - fac9*r0;
                    T p_9_1 = l_9_1 - fac9*r1;
                    T p_9_2 = l_9_2 - fac9*r2;
                    T p_9_3 = l_9_3 - fac9*r3;
                    T p_9_4 = l_9_4 - fac9*r4;
                    T p_9_5 = l_9_5 - fac9*r5;
                    T p_9_6 = l_9_6 - fac9*r6;
                    T p_9_7 = l_9_7 - fac9*r7;
                    T p_9_8 = l_9_8 - fac9*r8;
                    T p_9_9 = l_9_9 - fac9*r9;
                    T p_9_10 = l_9_10 - fac9*r10;
                    T p_9_11 = l_9_11 - fac9*r11;
                    T p_9_12 = l_9_12 - fac9*r12;
                    T p_9_13 = l_9_13 - fac9*r13;
                    T p_9_14 = l_9_14 - fac9*r14;
                    T p_9_15 = l_9_15 - fac9*r15;

                    tabp[(i+9)*width+j+0] = p_9_0;
                    tabp[(i+9)*width+j+1] = p_9_1;
                    tabp[(i+9)*width+j+2] = p_9_2;
                    tabp[(i+9)*width+j+3] = p_9_3;
                    tabp[(i+9)*width+j+4] = p_9_4;
                    tabp[(i+9)*width+j+5] = p_9_5;
                    tabp[(i+9)*width+j+6] = p_9_6;
                    tabp[(i+9)*width+j+7] = p_9_7;
                    tabp[(i+9)*width+j+8] = p_9_8;
                    tabp[(i+9)*width+j+9] = p_9_9;
                    tabp[(i+9)*width+j+10] = p_9_10;
                    tabp[(i+9)*width+j+11] = p_9_11;
                    tabp[(i+9)*width+j+12] = p_9_12;
                    tabp[(i+9)*width+j+13] = p_9_13;
                    tabp[(i+9)*width+j+14] = p_9_14;
                    tabp[(i+9)*width+j+15] = p_9_15;
                }

                //---------- i + 10 ----------
                if(i+10 != row) {
                    T l_10_0 = tabp[(i+10)*width+j+0];
                    T l_10_1 = tabp[(i+10)*width+j+1];
                    T l_10_2 = tabp[(i+10)*width+j+2];
                    T l_10_3 = tabp[(i+10)*width+j+3];
                    T l_10_4 = tabp[(i+10)*width+j+4];
                    T l_10_5 = tabp[(i+10)*width+j+5];
                    T l_10_6 = tabp[(i+10)*width+j+6];
                    T l_10_7 = tabp[(i+10)*width+j+7];
                    T l_10_8 = tabp[(i+10)*width+j+8];
                    T l_10_9 = tabp[(i+10)*width+j+9];
                    T l_10_10 = tabp[(i+10)*width+j+10];
                    T l_10_11 = tabp[(i+10)*width+j+11];
                    T l_10_12 = tabp[(i+10)*width+j+12];
                    T l_10_13 = tabp[(i+10)*width+j+13];
                    T l_10_14 = tabp[(i+10)*width+j+14];
                    T l_10_15 = tabp[(i+10)*width+j+15];

                    T p_10_0 = l_10_0 - fac10*r0;
                    T p_10_1 = l_10_1 - fac10*r1;
                    T p_10_2 = l_10_2 - fac10*r2;
                    T p_10_3 = l_10_3 - fac10*r3;
                    T p_10_4 = l_10_4 - fac10*r4;
                    T p_10_5 = l_10_5 - fac10*r5;
                    T p_10_6 = l_10_6 - fac10*r6;
                    T p_10_7 = l_10_7 - fac10*r7;
                    T p_10_8 = l_10_8 - fac10*r8;
                    T p_10_9 = l_10_9 - fac10*r9;
                    T p_10_10 = l_10_10 - fac10*r10;
                    T p_10_11 = l_10_11 - fac10*r11;
                    T p_10_12 = l_10_12 - fac10*r12;
                    T p_10_13 = l_10_13 - fac10*r13;
                    T p_10_14 = l_10_14 - fac10*r14;
                    T p_10_15 = l_10_15 - fac10*r15;

                    tabp[(i+10)*width+j+0] = p_10_0;
                    tabp[(i+10)*width+j+1] = p_10_1;
                    tabp[(i+10)*width+j+2] = p_10_2;
                    tabp[(i+10)*width+j+3] = p_10_3;
                    tabp[(i+10)*width+j+4] = p_10_4;
                    tabp[(i+10)*width+j+5] = p_10_5;
                    tabp[(i+10)*width+j+6] = p_10_6;
                    tabp[(i+10)*width+j+7] = p_10_7;
                    tabp[(i+10)*width+j+8] = p_10_8;
                    tabp[(i+10)*width+j+9] = p_10_9;
                    tabp[(i+10)*width+j+10] = p_10_10;
                    tabp[(i+10)*width+j+11] = p_10_11;
                    tabp[(i+10)*width+j+12] = p_10_12;
                    tabp[(i+10)*width+j+13] = p_10_13;
                    tabp[(i+10)*width+j+14] = p_10_14;
                    tabp[(i+10)*width+j+15] = p_10_15;
                }

                //---------- i + 11 ----------
                if(i+11 != row) {
                    T l_11_0 = tabp[(i+11)*width+j+0];
                    T l_11_1 = tabp[(i+11)*width+j+1];
                    T l_11_2 = tabp[(i+11)*width+j+2];
                    T l_11_3 = tabp[(i+11)*width+j+3];
                    T l_11_4 = tabp[(i+11)*width+j+4];
                    T l_11_5 = tabp[(i+11)*width+j+5];
                    T l_11_6 = tabp[(i+11)*width+j+6];
                    T l_11_7 = tabp[(i+11)*width+j+7];
                    T l_11_8 = tabp[(i+11)*width+j+8];
                    T l_11_9 = tabp[(i+11)*width+j+9];
                    T l_11_10 = tabp[(i+11)*width+j+10];
                    T l_11_11 = tabp[(i+11)*width+j+11];
                    T l_11_12 = tabp[(i+11)*width+j+12];
                    T l_11_13 = tabp[(i+11)*width+j+13];
                    T l_11_14 = tabp[(i+11)*width+j+14];
                    T l_11_15 = tabp[(i+11)*width+j+15];

                    T p_11_0 = l_11_0 - fac11*r0;
                    T p_11_1 = l_11_1 - fac11*r1;
                    T p_11_2 = l_11_2 - fac11*r2;
                    T p_11_3 = l_11_3 - fac11*r3;
                    T p_11_4 = l_11_4 - fac11*r4;
                    T p_11_5 = l_11_5 - fac11*r5;
                    T p_11_6 = l_11_6 - fac11*r6;
                    T p_11_7 = l_11_7 - fac11*r7;
                    T p_11_8 = l_11_8 - fac11*r8;
                    T p_11_9 = l_11_9 - fac11*r9;
                    T p_11_10 = l_11_10 - fac11*r10;
                    T p_11_11 = l_11_11 - fac11*r11;
                    T p_11_12 = l_11_12 - fac11*r12;
                    T p_11_13 = l_11_13 - fac11*r13;
                    T p_11_14 = l_11_14 - fac11*r14;
                    T p_11_15 = l_11_15 - fac11*r15;

                    tabp[(i+11)*width+j+0] = p_11_0;
                    tabp[(i+11)*width+j+1] = p_11_1;
                    tabp[(i+11)*width+j+2] = p_11_2;
                    tabp[(i+11)*width+j+3] = p_11_3;
                    tabp[(i+11)*width+j+4] = p_11_4;
                    tabp[(i+11)*width+j+5] = p_11_5;
                    tabp[(i+11)*width+j+6] = p_11_6;
                    tabp[(i+11)*width+j+7] = p_11_7;
                    tabp[(i+11)*width+j+8] = p_11_8;
                    tabp[(i+11)*width+j+9] = p_11_9;
                    tabp[(i+11)*width+j+10] = p_11_10;
                    tabp[(i+11)*width+j+11] = p_11_11;
                    tabp[(i+11)*width+j+12] = p_11_12;
                    tabp[(i+11)*width+j+13] = p_11_13;
                    tabp[(i+11)*width+j+14] = p_11_14;
                    tabp[(i+11)*width+j+15] = p_11_15;
                }

                //---------- i + 12 ----------
                if(i+12 != row) {
                    T l_12_0 = tabp[(i+12)*width+j+0];
                    T l_12_1 = tabp[(i+12)*width+j+1];
                    T l_12_2 = tabp[(i+12)*width+j+2];
                    T l_12_3 = tabp[(i+12)*width+j+3];
                    T l_12_4 = tabp[(i+12)*width+j+4];
                    T l_12_5 = tabp[(i+12)*width+j+5];
                    T l_12_6 = tabp[(i+12)*width+j+6];
                    T l_12_7 = tabp[(i+12)*width+j+7];
                    T l_12_8 = tabp[(i+12)*width+j+8];
                    T l_12_9 = tabp[(i+12)*width+j+9];
                    T l_12_10 = tabp[(i+12)*width+j+10];
                    T l_12_11 = tabp[(i+12)*width+j+11];
                    T l_12_12 = tabp[(i+12)*width+j+12];
                    T l_12_13 = tabp[(i+12)*width+j+13];
                    T l_12_14 = tabp[(i+12)*width+j+14];
                    T l_12_15 = tabp[(i+12)*width+j+15];

                    T p_12_0 = l_12_0 - fac12*r0;
                    T p_12_1 = l_12_1 - fac12*r1;
                    T p_12_2 = l_12_2 - fac12*r2;
                    T p_12_3 = l_12_3 - fac12*r3;
                    T p_12_4 = l_12_4 - fac12*r4;
                    T p_12_5 = l_12_5 - fac12*r5;
                    T p_12_6 = l_12_6 - fac12*r6;
                    T p_12_7 = l_12_7 - fac12*r7;
                    T p_12_8 = l_12_8 - fac12*r8;
                    T p_12_9 = l_12_9 - fac12*r9;
                    T p_12_10 = l_12_10 - fac12*r10;
                    T p_12_11 = l_12_11 - fac12*r11;
                    T p_12_12 = l_12_12 - fac12*r12;
                    T p_12_13 = l_12_13 - fac12*r13;
                    T p_12_14 = l_12_14 - fac12*r14;
                    T p_12_15 = l_12_15 - fac12*r15;

                    tabp[(i+12)*width+j+0] = p_12_0;
                    tabp[(i+12)*width+j+1] = p_12_1;
                    tabp[(i+12)*width+j+2] = p_12_2;
                    tabp[(i+12)*width+j+3] = p_12_3;
                    tabp[(i+12)*width+j+4] = p_12_4;
                    tabp[(i+12)*width+j+5] = p_12_5;
                    tabp[(i+12)*width+j+6] = p_12_6;
                    tabp[(i+12)*width+j+7] = p_12_7;
                    tabp[(i+12)*width+j+8] = p_12_8;
                    tabp[(i+12)*width+j+9] = p_12_9;
                    tabp[(i+12)*width+j+10] = p_12_10;
                    tabp[(i+12)*width+j+11] = p_12_11;
                    tabp[(i+12)*width+j+12] = p_12_12;
                    tabp[(i+12)*width+j+13] = p_12_13;
                    tabp[(i+12)*width+j+14] = p_12_14;
                    tabp[(i+12)*width+j+15] = p_12_15;
                }

                //---------- i + 13 ----------
                if(i+13 != row) {
                    T l_13_0 = tabp[(i+13)*width+j+0];
                    T l_13_1 = tabp[(i+13)*width+j+1];
                    T l_13_2 = tabp[(i+13)*width+j+2];
                    T l_13_3 = tabp[(i+13)*width+j+3];
                    T l_13_4 = tabp[(i+13)*width+j+4];
                    T l_13_5 = tabp[(i+13)*width+j+5];
                    T l_13_6 = tabp[(i+13)*width+j+6];
                    T l_13_7 = tabp[(i+13)*width+j+7];
                    T l_13_8 = tabp[(i+13)*width+j+8];
                    T l_13_9 = tabp[(i+13)*width+j+9];
                    T l_13_10 = tabp[(i+13)*width+j+10];
                    T l_13_11 = tabp[(i+13)*width+j+11];
                    T l_13_12 = tabp[(i+13)*width+j+12];
                    T l_13_13 = tabp[(i+13)*width+j+13];
                    T l_13_14 = tabp[(i+13)*width+j+14];
                    T l_13_15 = tabp[(i+13)*width+j+15];

                    T p_13_0 = l_13_0 - fac13*r0;
                    T p_13_1 = l_13_1 - fac13*r1;
                    T p_13_2 = l_13_2 - fac13*r2;
                    T p_13_3 = l_13_3 - fac13*r3;
                    T p_13_4 = l_13_4 - fac13*r4;
                    T p_13_5 = l_13_5 - fac13*r5;
                    T p_13_6 = l_13_6 - fac13*r6;
                    T p_13_7 = l_13_7 - fac13*r7;
                    T p_13_8 = l_13_8 - fac13*r8;
                    T p_13_9 = l_13_9 - fac13*r9;
                    T p_13_10 = l_13_10 - fac13*r10;
                    T p_13_11 = l_13_11 - fac13*r11;
                    T p_13_12 = l_13_12 - fac13*r12;
                    T p_13_13 = l_13_13 - fac13*r13;
                    T p_13_14 = l_13_14 - fac13*r14;
                    T p_13_15 = l_13_15 - fac13*r15;

                    tabp[(i+13)*width+j+0] = p_13_0;
                    tabp[(i+13)*width+j+1] = p_13_1;
                    tabp[(i+13)*width+j+2] = p_13_2;
                    tabp[(i+13)*width+j+3] = p_13_3;
                    tabp[(i+13)*width+j+4] = p_13_4;
                    tabp[(i+13)*width+j+5] = p_13_5;
                    tabp[(i+13)*width+j+6] = p_13_6;
                    tabp[(i+13)*width+j+7] = p_13_7;
                    tabp[(i+13)*width+j+8] = p_13_8;
                    tabp[(i+13)*width+j+9] = p_13_9;
                    tabp[(i+13)*width+j+10] = p_13_10;
                    tabp[(i+13)*width+j+11] = p_13_11;
                    tabp[(i+13)*width+j+12] = p_13_12;
                    tabp[(i+13)*width+j+13] = p_13_13;
                    tabp[(i+13)*width+j+14] = p_13_14;
                    tabp[(i+13)*width+j+15] = p_13_15;
                }

                //---------- i + 14 ----------
                if(i+14 != row) {
                    T l_14_0 = tabp[(i+14)*width+j+0];
                    T l_14_1 = tabp[(i+14)*width+j+1];
                    T l_14_2 = tabp[(i+14)*width+j+2];
                    T l_14_3 = tabp[(i+14)*width+j+3];
                    T l_14_4 = tabp[(i+14)*width+j+4];
                    T l_14_5 = tabp[(i+14)*width+j+5];
                    T l_14_6 = tabp[(i+14)*width+j+6];
                    T l_14_7 = tabp[(i+14)*width+j+7];
                    T l_14_8 = tabp[(i+14)*width+j+8];
                    T l_14_9 = tabp[(i+14)*width+j+9];
                    T l_14_10 = tabp[(i+14)*width+j+10];
                    T l_14_11 = tabp[(i+14)*width+j+11];
                    T l_14_12 = tabp[(i+14)*width+j+12];
                    T l_14_13 = tabp[(i+14)*width+j+13];
                    T l_14_14 = tabp[(i+14)*width+j+14];
                    T l_14_15 = tabp[(i+14)*width+j+15];

                    T p_14_0 = l_14_0 - fac14*r0;
                    T p_14_1 = l_14_1 - fac14*r1;
                    T p_14_2 = l_14_2 - fac14*r2;
                    T p_14_3 = l_14_3 - fac14*r3;
                    T p_14_4 = l_14_4 - fac14*r4;
                    T p_14_5 = l_14_5 - fac14*r5;
                    T p_14_6 = l_14_6 - fac14*r6;
                    T p_14_7 = l_14_7 - fac14*r7;
                    T p_14_8 = l_14_8 - fac14*r8;
                    T p_14_9 = l_14_9 - fac14*r9;
                    T p_14_10 = l_14_10 - fac14*r10;
                    T p_14_11 = l_14_11 - fac14*r11;
                    T p_14_12 = l_14_12 - fac14*r12;
                    T p_14_13 = l_14_13 - fac14*r13;
                    T p_14_14 = l_14_14 - fac14*r14;
                    T p_14_15 = l_14_15 - fac14*r15;

                    tabp[(i+14)*width+j+0] = p_14_0;
                    tabp[(i+14)*width+j+1] = p_14_1;
                    tabp[(i+14)*width+j+2] = p_14_2;
                    tabp[(i+14)*width+j+3] = p_14_3;
                    tabp[(i+14)*width+j+4] = p_14_4;
                    tabp[(i+14)*width+j+5] = p_14_5;
                    tabp[(i+14)*width+j+6] = p_14_6;
                    tabp[(i+14)*width+j+7] = p_14_7;
                    tabp[(i+14)*width+j+8] = p_14_8;
                    tabp[(i+14)*width+j+9] = p_14_9;
                    tabp[(i+14)*width+j+10] = p_14_10;
                    tabp[(i+14)*width+j+11] = p_14_11;
                    tabp[(i+14)*width+j+12] = p_14_12;
                    tabp[(i+14)*width+j+13] = p_14_13;
                    tabp[(i+14)*width+j+14] = p_14_14;
                    tabp[(i+14)*width+j+15] = p_14_15;
                }

                //---------- i + 15 ----------
                if(i+15 != row) {
                    T l_15_0 = tabp[(i+15)*width+j+0];
                    T l_15_1 = tabp[(i+15)*width+j+1];
                    T l_15_2 = tabp[(i+15)*width+j+2];
                    T l_15_3 = tabp[(i+15)*width+j+3];
                    T l_15_4 = tabp[(i+15)*width+j+4];
                    T l_15_5 = tabp[(i+15)*width+j+5];
                    T l_15_6 = tabp[(i+15)*width+j+6];
                    T l_15_7 = tabp[(i+15)*width+j+7];
                    T l_15_8 = tabp[(i+15)*width+j+8];
                    T l_15_9 = tabp[(i+15)*width+j+9];
                    T l_15_10 = tabp[(i+15)*width+j+10];
                    T l_15_11 = tabp[(i+15)*width+j+11];
                    T l_15_12 = tabp[(i+15)*width+j+12];
                    T l_15_13 = tabp[(i+15)*width+j+13];
                    T l_15_14 = tabp[(i+15)*width+j+14];
                    T l_15_15 = tabp[(i+15)*width+j+15];

                    T p_15_0 = l_15_0 - fac15*r0;
                    T p_15_1 = l_15_1 - fac15*r1;
                    T p_15_2 = l_15_2 - fac15*r2;
                    T p_15_3 = l_15_3 - fac15*r3;
                    T p_15_4 = l_15_4 - fac15*r4;
                    T p_15_5 = l_15_5 - fac15*r5;
                    T p_15_6 = l_15_6 - fac15*r6;
                    T p_15_7 = l_15_7 - fac15*r7;
                    T p_15_8 = l_15_8 - fac15*r8;
                    T p_15_9 = l_15_9 - fac15*r9;
                    T p_15_10 = l_15_10 - fac15*r10;
                    T p_15_11 = l_15_11 - fac15*r11;
                    T p_15_12 = l_15_12 - fac15*r12;
                    T p_15_13 = l_15_13 - fac15*r13;
                    T p_15_14 = l_15_14 - fac15*r14;
                    T p_15_15 = l_15_15 - fac15*r15;

                    tabp[(i+15)*width+j+0] = p_15_0;
                    tabp[(i+15)*width+j+1] = p_15_1;
                    tabp[(i+15)*width+j+2] = p_15_2;
                    tabp[(i+15)*width+j+3] = p_15_3;
                    tabp[(i+15)*width+j+4] = p_15_4;
                    tabp[(i+15)*width+j+5] = p_15_5;
                    tabp[(i+15)*width+j+6] = p_15_6;
                    tabp[(i+15)*width+j+7] = p_15_7;
                    tabp[(i+15)*width+j+8] = p_15_8;
                    tabp[(i+15)*width+j+9] = p_15_9;
                    tabp[(i+15)*width+j+10] = p_15_10;
                    tabp[(i+15)*width+j+11] = p_15_11;
                    tabp[(i+15)*width+j+12] = p_15_12;
                    tabp[(i+15)*width+j+13] = p_15_13;
                    tabp[(i+15)*width+j+14] = p_15_14;
                    tabp[(i+15)*width+j+15] = p_15_15;
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
                if(i+2 != row) {
                    tabp[(i+2)*width+j] -= fac2*r1;
                }
                if(i+3 != row) {
                    tabp[(i+3)*width+j] -= fac3*r1;
                }
                if(i+4 != row) {
                    tabp[(i+4)*width+j] -= fac4*r1;
                }
                if(i+5 != row) {
                    tabp[(i+5)*width+j] -= fac5*r1;
                }
                if(i+6 != row) {
                    tabp[(i+6)*width+j] -= fac6*r1;
                }
                if(i+7 != row) {
                    tabp[(i+7)*width+j] -= fac7*r1;
                }
                if(i+8 != row) {
                    tabp[(i+8)*width+j] -= fac8*r1;
                }
                if(i+9 != row) {
                    tabp[(i+9)*width+j] -= fac9*r1;
                }
                if(i+10 != row) {
                    tabp[(i+10)*width+j] -= fac10*r1;
                }
                if(i+11 != row) {
                    tabp[(i+11)*width+j] -= fac11*r1;
                }
                if(i+12 != row) {
                    tabp[(i+12)*width+j] -= fac12*r1;
                }
                if(i+13 != row) {
                    tabp[(i+13)*width+j] -= fac13*r1;
                }
                if(i+14 != row) {
                    tabp[(i+14)*width+j] -= fac14*r1;
                }
                if(i+15 != row) {
                    tabp[(i+15)*width+j] -= fac15*r1;
                }
            }
        }

        for(int i = m-(m%16); i < m+1; ++i) {
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
