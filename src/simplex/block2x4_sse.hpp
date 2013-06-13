/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "Simplex.hpp"

template <typename T>
class Simplex_block2x4_sse : public SimplexBase<T> {

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

    std::string get_identifier() { return "block2x4_sse"; }

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
        for(int i = 0; i < m; i += 2) {     // peel off m+1, as we expect an even amount of dofs
            PERFC_MEM+=2; PERFC_ADDMUL+=2;
            T fac1 = tabp[i*width+col] * ipiv;
            T fac2 = tabp[(i+1)*width+col] * ipiv;
            PERFC_MEM+=2;
            __m128d f1 = _mm_set1_pd(fac1);
            __m128d f2 = _mm_set1_pd(fac2);

            for(int j = 0; j < width-(width%4); j += 4) {

                //~ PERFC_MEM += 4;  // not from memory, as width*sizeof(T) ~ 1-24 kB should easily fit into L2/L3
                __m128d r1 = _mm_load_pd(tabp+row*width+j);
                __m128d r2 = _mm_load_pd(tabp+row*width+j+2);

                if(i != row) {
                    PERFC_MEM += 4;
                    __m128d l1 = _mm_load_pd(tabp+i*width+j);
                    __m128d l2 = _mm_load_pd(tabp+i*width+j+2);

                    PERFC_ADDMUL += 8;
                    __m128d p1 = _mm_mul_pd(f1, r1);
                    __m128d p2 = _mm_mul_pd(f1, r2);
                    __m128d p3 = _mm_sub_pd(l1, p1);
                    __m128d p4 = _mm_sub_pd(l2, p2);

                    _mm_store_pd(tabp+i*width+j, p3);
                    _mm_store_pd(tabp+i*width+j+2, p4);
                }
                if(i+1 != row) {
                    PERFC_MEM += 4;
                    __m128d l1 = _mm_load_pd(tabp+(i+1)*width+j);
                    __m128d l2 = _mm_load_pd(tabp+(i+1)*width+j+2);

                    PERFC_ADDMUL += 8;
                    __m128d p1 = _mm_mul_pd(f2, r1);
                    __m128d p2 = _mm_mul_pd(f2, r2);
                    __m128d p3 = _mm_sub_pd(l1, p1);
                    __m128d p4 = _mm_sub_pd(l2, p2);

                    _mm_store_pd(tabp+(i+1)*width+j, p3);
                    _mm_store_pd(tabp+(i+1)*width+j+2, p4);
                }
            }

            for(int j = width-(width%4); j < width; ++j) {
                if(i != row) {
                    PERFC_ADDMUL+=2; PERFC_MEM+=2;
                    tabp[i*width+j] -= fac1*tabp[row*width+j];
                }
                if(i+1 != row) {
                    PERFC_ADDMUL+=2; PERFC_MEM+=2;
                    tabp[(i+1)*width+j] -= fac2*tabp[row*width+j];
                }
            }

        }
        active[row] = col;
        ++PERFC_ADDMUL; ++PERFC_MEM;
        T fac = tabp[m*width+col]*ipiv;
        for(int j = 0; j < width; ++j) {
            PERFC_ADDMUL += 2; ++PERFC_MEM;
            tabp[m*width+j] -= fac*tabp[row*width+j];
        }


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
