/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "Simplex.hpp"

template <typename T>
class SimplexBlockAVX : public SimplexBase<T> {

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

    std::string get_identifier() { return "block_avx"; }

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
        for(int i = 0; i < m; i += 2) {     // peel off m+1, as we temporarily store the pivot row there
            PERFC_MEM+=2; PERFC_ADDMUL+=2;
            T fac1 = tabp[i*width+col] * ipiv;
            T fac2 = tabp[(i+1)*width+col] * ipiv;
            PERFC_MEM+=2;
            __m256d f1 = _mm256_set1_pd(fac1);
            __m256d f2 = _mm256_set1_pd(fac2);

            PERFC_ADDMUL += 4*width; PERFC_MEM += 2*width;

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

                __m256d r = _mm256_loadu_pd(tabp+m*width+j);

                __m256d la = _mm256_loadu_pd(tabp+i*width+j);

                __m256d pa1 = _mm256_mul_pd(f1, r);
                __m256d pa2 = _mm256_sub_pd(la, pa1);

                _mm256_storeu_pd(tabp+i*width+j, pa2); // segfaults

                __m256d lb = _mm256_loadu_pd(tabp+(i+1)*width+j);

                __m256d pb1 = _mm256_mul_pd(f2, r);
                __m256d pb2 = _mm256_sub_pd(lb, pb1);

                _mm256_storeu_pd(tabp+(i+1)*width+j, pb2);
            }

            for(int j = aligned_end; j < width; ++j) {
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
