/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "Simplex.hpp"

template <typename T>
class Simplex_block16x4_swap_avx : public SimplexBase<T> {

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

    std::string get_identifier() { return "block16x4_swap_avx"; }

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

        for(int i = 0; i < m-(m%16); i += 16) {
            PERFC_MEM+=2*16; PERFC_ADDMUL+=16;
            T fac0 = tabp[(i+0)*width+col] * ipiv;
            __m256d f0 = _mm256_set1_pd(fac0);
            T fac1 = tabp[(i+1)*width+col] * ipiv;
            __m256d f1 = _mm256_set1_pd(fac1);
            T fac2 = tabp[(i+2)*width+col] * ipiv;
            __m256d f2 = _mm256_set1_pd(fac2);
            T fac3 = tabp[(i+3)*width+col] * ipiv;
            __m256d f3 = _mm256_set1_pd(fac3);
            T fac4 = tabp[(i+4)*width+col] * ipiv;
            __m256d f4 = _mm256_set1_pd(fac4);
            T fac5 = tabp[(i+5)*width+col] * ipiv;
            __m256d f5 = _mm256_set1_pd(fac5);
            T fac6 = tabp[(i+6)*width+col] * ipiv;
            __m256d f6 = _mm256_set1_pd(fac6);
            T fac7 = tabp[(i+7)*width+col] * ipiv;
            __m256d f7 = _mm256_set1_pd(fac7);
            T fac8 = tabp[(i+8)*width+col] * ipiv;
            __m256d f8 = _mm256_set1_pd(fac8);
            T fac9 = tabp[(i+9)*width+col] * ipiv;
            __m256d f9 = _mm256_set1_pd(fac9);
            T fac10 = tabp[(i+10)*width+col] * ipiv;
            __m256d f10 = _mm256_set1_pd(fac10);
            T fac11 = tabp[(i+11)*width+col] * ipiv;
            __m256d f11 = _mm256_set1_pd(fac11);
            T fac12 = tabp[(i+12)*width+col] * ipiv;
            __m256d f12 = _mm256_set1_pd(fac12);
            T fac13 = tabp[(i+13)*width+col] * ipiv;
            __m256d f13 = _mm256_set1_pd(fac13);
            T fac14 = tabp[(i+14)*width+col] * ipiv;
            __m256d f14 = _mm256_set1_pd(fac14);
            T fac15 = tabp[(i+15)*width+col] * ipiv;
            __m256d f15 = _mm256_set1_pd(fac15);
            
            for(int j = 0; j < width-(width%4); j += 4) {
                __m256d r0 = _mm256_load_pd(tabp+m*width+j+0);

                //---------- i + 0 ----------
                PERFC_MEM += 4;
                __m256d l_0_0 = _mm256_load_pd(tabp+(i+0)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_0_0 = _mm256_mul_pd(f0, r0);
                __m256d q_0_0 = _mm256_sub_pd(l_0_0, p_0_0);

                _mm256_store_pd(tabp+(i+0)*width+j+0, q_0_0);

                //---------- i + 1 ----------
                PERFC_MEM += 4;
                __m256d l_1_0 = _mm256_load_pd(tabp+(i+1)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_1_0 = _mm256_mul_pd(f1, r0);
                __m256d q_1_0 = _mm256_sub_pd(l_1_0, p_1_0);

                _mm256_store_pd(tabp+(i+1)*width+j+0, q_1_0);

                //---------- i + 2 ----------
                PERFC_MEM += 4;
                __m256d l_2_0 = _mm256_load_pd(tabp+(i+2)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_2_0 = _mm256_mul_pd(f2, r0);
                __m256d q_2_0 = _mm256_sub_pd(l_2_0, p_2_0);

                _mm256_store_pd(tabp+(i+2)*width+j+0, q_2_0);

                //---------- i + 3 ----------
                PERFC_MEM += 4;
                __m256d l_3_0 = _mm256_load_pd(tabp+(i+3)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_3_0 = _mm256_mul_pd(f3, r0);
                __m256d q_3_0 = _mm256_sub_pd(l_3_0, p_3_0);

                _mm256_store_pd(tabp+(i+3)*width+j+0, q_3_0);

                //---------- i + 4 ----------
                PERFC_MEM += 4;
                __m256d l_4_0 = _mm256_load_pd(tabp+(i+4)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_4_0 = _mm256_mul_pd(f4, r0);
                __m256d q_4_0 = _mm256_sub_pd(l_4_0, p_4_0);

                _mm256_store_pd(tabp+(i+4)*width+j+0, q_4_0);

                //---------- i + 5 ----------
                PERFC_MEM += 4;
                __m256d l_5_0 = _mm256_load_pd(tabp+(i+5)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_5_0 = _mm256_mul_pd(f5, r0);
                __m256d q_5_0 = _mm256_sub_pd(l_5_0, p_5_0);

                _mm256_store_pd(tabp+(i+5)*width+j+0, q_5_0);

                //---------- i + 6 ----------
                PERFC_MEM += 4;
                __m256d l_6_0 = _mm256_load_pd(tabp+(i+6)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_6_0 = _mm256_mul_pd(f6, r0);
                __m256d q_6_0 = _mm256_sub_pd(l_6_0, p_6_0);

                _mm256_store_pd(tabp+(i+6)*width+j+0, q_6_0);

                //---------- i + 7 ----------
                PERFC_MEM += 4;
                __m256d l_7_0 = _mm256_load_pd(tabp+(i+7)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_7_0 = _mm256_mul_pd(f7, r0);
                __m256d q_7_0 = _mm256_sub_pd(l_7_0, p_7_0);

                _mm256_store_pd(tabp+(i+7)*width+j+0, q_7_0);

                //---------- i + 8 ----------
                PERFC_MEM += 4;
                __m256d l_8_0 = _mm256_load_pd(tabp+(i+8)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_8_0 = _mm256_mul_pd(f8, r0);
                __m256d q_8_0 = _mm256_sub_pd(l_8_0, p_8_0);

                _mm256_store_pd(tabp+(i+8)*width+j+0, q_8_0);

                //---------- i + 9 ----------
                PERFC_MEM += 4;
                __m256d l_9_0 = _mm256_load_pd(tabp+(i+9)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_9_0 = _mm256_mul_pd(f9, r0);
                __m256d q_9_0 = _mm256_sub_pd(l_9_0, p_9_0);

                _mm256_store_pd(tabp+(i+9)*width+j+0, q_9_0);

                //---------- i + 10 ----------
                PERFC_MEM += 4;
                __m256d l_10_0 = _mm256_load_pd(tabp+(i+10)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_10_0 = _mm256_mul_pd(f10, r0);
                __m256d q_10_0 = _mm256_sub_pd(l_10_0, p_10_0);

                _mm256_store_pd(tabp+(i+10)*width+j+0, q_10_0);

                //---------- i + 11 ----------
                PERFC_MEM += 4;
                __m256d l_11_0 = _mm256_load_pd(tabp+(i+11)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_11_0 = _mm256_mul_pd(f11, r0);
                __m256d q_11_0 = _mm256_sub_pd(l_11_0, p_11_0);

                _mm256_store_pd(tabp+(i+11)*width+j+0, q_11_0);

                //---------- i + 12 ----------
                PERFC_MEM += 4;
                __m256d l_12_0 = _mm256_load_pd(tabp+(i+12)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_12_0 = _mm256_mul_pd(f12, r0);
                __m256d q_12_0 = _mm256_sub_pd(l_12_0, p_12_0);

                _mm256_store_pd(tabp+(i+12)*width+j+0, q_12_0);

                //---------- i + 13 ----------
                PERFC_MEM += 4;
                __m256d l_13_0 = _mm256_load_pd(tabp+(i+13)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_13_0 = _mm256_mul_pd(f13, r0);
                __m256d q_13_0 = _mm256_sub_pd(l_13_0, p_13_0);

                _mm256_store_pd(tabp+(i+13)*width+j+0, q_13_0);

                //---------- i + 14 ----------
                PERFC_MEM += 4;
                __m256d l_14_0 = _mm256_load_pd(tabp+(i+14)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_14_0 = _mm256_mul_pd(f14, r0);
                __m256d q_14_0 = _mm256_sub_pd(l_14_0, p_14_0);

                _mm256_store_pd(tabp+(i+14)*width+j+0, q_14_0);

                //---------- i + 15 ----------
                PERFC_MEM += 4;
                __m256d l_15_0 = _mm256_load_pd(tabp+(i+15)*width+j+0);

                PERFC_ADDMUL += 2*4;
                __m256d p_15_0 = _mm256_mul_pd(f15, r0);
                __m256d q_15_0 = _mm256_sub_pd(l_15_0, p_15_0);

                _mm256_store_pd(tabp+(i+15)*width+j+0, q_15_0);
            }

            for(int j = width-(width%4); j < width; ++j) {
                PERFC_MEM += 1;
                T r1 = tabp[m*width+j];

                PERFC_ADDMUL += 2*16;
                tabp[(i+0)*width+j] -= fac0*r1;
                tabp[(i+1)*width+j] -= fac1*r1;
                tabp[(i+2)*width+j] -= fac2*r1;
                tabp[(i+3)*width+j] -= fac3*r1;
                tabp[(i+4)*width+j] -= fac4*r1;
                tabp[(i+5)*width+j] -= fac5*r1;
                tabp[(i+6)*width+j] -= fac6*r1;
                tabp[(i+7)*width+j] -= fac7*r1;
                tabp[(i+8)*width+j] -= fac8*r1;
                tabp[(i+9)*width+j] -= fac9*r1;
                tabp[(i+10)*width+j] -= fac10*r1;
                tabp[(i+11)*width+j] -= fac11*r1;
                tabp[(i+12)*width+j] -= fac12*r1;
                tabp[(i+13)*width+j] -= fac13*r1;
                tabp[(i+14)*width+j] -= fac14*r1;
                tabp[(i+15)*width+j] -= fac15*r1;
            }
        }

        for(int i = m-(m%16); i < m; ++i) {
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
