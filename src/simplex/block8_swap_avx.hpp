/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "Simplex.hpp"

template <typename T>
class Simplex_block8_swap_avx : public SimplexBase<T> {

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

    std::string get_identifier() { return "block8_swap_avx"; }

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
        for(int i = 0; i < m-(m%8); i += 8) {     // peel off m+1, as we temporarily store the pivot row there
            PERFC_MEM+=8; PERFC_ADDMUL+=8;
            T fac1 = tabp[i*width+col] * ipiv;
            T fac2 = tabp[(i+1)*width+col] * ipiv;
            T fac3 = tabp[(i+2)*width+col] * ipiv;
            T fac4 = tabp[(i+3)*width+col] * ipiv;
            T fac5 = tabp[(i+4)*width+col] * ipiv;
            T fac6 = tabp[(i+5)*width+col] * ipiv;
            T fac7 = tabp[(i+6)*width+col] * ipiv;
            T fac8 = tabp[(i+7)*width+col] * ipiv;
            PERFC_MEM+=8;
            __m256d f1 = _mm256_set1_pd(fac1);
            __m256d f2 = _mm256_set1_pd(fac2);
            __m256d f3 = _mm256_set1_pd(fac3);
            __m256d f4 = _mm256_set1_pd(fac4);
            __m256d f5 = _mm256_set1_pd(fac5);
            __m256d f6 = _mm256_set1_pd(fac6);
            __m256d f7 = _mm256_set1_pd(fac7);
            __m256d f8 = _mm256_set1_pd(fac8);

            PERFC_ADDMUL += 16*width; PERFC_MEM += 8*width;

            int peel = (long long)tabp & 0x1f; /* tabp % 32 */
            if(peel != 0) {
                peel = (32 - peel)/sizeof(T);
                for (int j = 0; j < peel; j++) {
                    tabp[i*width+j] -= fac1*tabp[m*width+j];
                    tabp[(i+1)*width+j] -= fac2*tabp[m*width+j];
                }
            }

            int aligned_end = width - (width%4) - peel;

            for(int j = peel; j < aligned_end; j += 4) {

                __m256d r = _mm256_load_pd(tabp+m*width+j);

                // row 1
                __m256d l1 = _mm256_load_pd(tabp+i*width+j);
                __m256d t1a = _mm256_mul_pd(f1, r);
                __m256d t1b = _mm256_sub_pd(l1, t1a);
                _mm256_store_pd(tabp+i*width+j, t1b);

                // row 2
                __m256d l2 = _mm256_load_pd(tabp+(i+1)*width+j);
                __m256d t2a = _mm256_mul_pd(f2, r);
                __m256d t2b = _mm256_sub_pd(l2, t2a);
                _mm256_store_pd(tabp+(i+1)*width+j, t2b);

                // row 3
                __m256d l3 = _mm256_load_pd(tabp+(i+2)*width+j);
                __m256d t3a = _mm256_mul_pd(f3, r);
                __m256d t3b = _mm256_sub_pd(l3, t3a);
                _mm256_store_pd(tabp+(i+2)*width+j, t3b);

                // row 4
                __m256d l4 = _mm256_load_pd(tabp+(i+3)*width+j);
                __m256d t4a = _mm256_mul_pd(f4, r);
                __m256d t4b = _mm256_sub_pd(l4, t4a);
                _mm256_store_pd(tabp+(i+3)*width+j, t4b);

                // row 5
                __m256d l5 = _mm256_load_pd(tabp+(i+4)*width+j);
                __m256d t5a = _mm256_mul_pd(f5, r);
                __m256d t5b = _mm256_sub_pd(l5, t5a);
                _mm256_store_pd(tabp+(i+4)*width+j, t5b);

                // row 6
                __m256d l6 = _mm256_load_pd(tabp+(i+5)*width+j);
                __m256d t6a = _mm256_mul_pd(f6, r);
                __m256d t6b = _mm256_sub_pd(l6, t6a);
                _mm256_store_pd(tabp+(i+5)*width+j, t6b);

                // row 7
                __m256d l7 = _mm256_load_pd(tabp+(i+6)*width+j);
                __m256d t7a = _mm256_mul_pd(f7, r);
                __m256d t7b = _mm256_sub_pd(l7, t7a);
                _mm256_store_pd(tabp+(i+6)*width+j, t7b);

                // row 8
                __m256d l8 = _mm256_load_pd(tabp+(i+7)*width+j);
                __m256d t8a = _mm256_mul_pd(f8, r);
                __m256d t8b = _mm256_sub_pd(l8, t8a);
                _mm256_store_pd(tabp+(i+7)*width+j, t8b);

            }

            for(int j = aligned_end; j < width; ++j) {
                tabp[i*width+j] -= fac1*tabp[m*width+j];
                tabp[(i+1)*width+j] -= fac2*tabp[m*width+j];
            }

        }

        for(int i = m-(m%8); i < m; ++i) {  // FIXME: we could gain a bit by SSA+AVX on this tail loop
            T fac = tabp[i*width+col] * ipiv;
            for(int j = 0; j < width; ++j) {
                PERFC_ADDMUL += 2;
                ++PERFC_MEM;
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
