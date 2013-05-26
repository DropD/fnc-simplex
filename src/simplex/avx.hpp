/*
Assumptions:
    See base class.
*/



#pragma once

//~ #include <immintrin.h>  // AVX, but we're aligining to 128 bits
#include <x86intrin.h>  // pulls depending on march

#include "Simplex.hpp"

template <typename T>
class SimplexAVX : public SimplexBase<T> {

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

    std::string get_identifier() { return "avx"; }

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
        for(int i = 0; i < m+1; ++i) {
            if(i == row)
                continue;
            ++PERFC_DIV; ++PERFC_MEM;
            T fac = tabp[i*width+col]/pivot;
            PERFC_ADDMUL += 2*width; PERFC_MEM += width;
            __m256d l, r, f;
            f = _mm256_set1_pd(fac);

            int peel = (long long)tabp & 0x1f; /* tabp % 32 */
            if (peel != 0) {
                peel = (32 - peel)/sizeof(T);
                for (int j = 0; j < peel; j++)
                    tabp[i*width+j] -= fac*tabp[row*width+j];
            }

            int aligned_end = width - (width%4) - peel;

            for(int j = peel; j < aligned_end; j += 4) {

                l = _mm256_loadu_pd(tabp+i*width+j);
                r = _mm256_loadu_pd(tabp+row*width+j);

                r = _mm256_mul_pd(r, f);
                l = _mm256_sub_pd(l, r);

                _mm256_storeu_pd(tabp + i*width+j, l);

            }

            for(int j = aligned_end; j < width; ++j) {
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
