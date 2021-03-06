/*
Assumptions:
    See base class.
    The number of degrees of freedom must be even.
*/



#pragma once

#include "../Simplex.hpp"

template <typename T>
class Simplex_block8x4_avx : public SimplexBase<T> {

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

    std::string get_identifier() { return "block8x4_avx"; }

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

        for(int i = 0; i < m-(m%8); i += 8) {
            PERFC_MEM+=2*8; PERFC_ADDMUL+=8;
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

            PERFC_ADDMUL += 2*8 * width;
            PERFC_MEM += 8*width;

            int peel = (long long)tabp & 0x1f; /* tabp % 32 */
            if(peel != 0) {
                peel = (32 - peel)/sizeof(T);
                for (int j = 0; j < peel; j++) {
                    tabp[(i+0)*width+j] -= fac0*tabp[m*width+j];
                    tabp[(i+1)*width+j] -= fac1*tabp[m*width+j];
                    tabp[(i+2)*width+j] -= fac2*tabp[m*width+j];
                    tabp[(i+3)*width+j] -= fac3*tabp[m*width+j];
                    tabp[(i+4)*width+j] -= fac4*tabp[m*width+j];
                    tabp[(i+5)*width+j] -= fac5*tabp[m*width+j];
                    tabp[(i+6)*width+j] -= fac6*tabp[m*width+j];
                    tabp[(i+7)*width+j] -= fac7*tabp[m*width+j];
                }
            }

            int aligned_end = width - (width%4) - peel;

            for(int j = peel; j < aligned_end; j += 4) {
                __m256d r0 = _mm256_load_pd(tabp+row*width+j+0);

                //---------- i + 0 ----------
                if(i+0 != row) {
		            __m256d l_0_0 = _mm256_load_pd(tabp+(i+0)*width+j+0);

		            __m256d p_0_0 = _mm256_mul_pd(f0, r0);
		            __m256d q_0_0 = _mm256_sub_pd(l_0_0, p_0_0);

		            _mm256_store_pd(tabp+(i+0)*width+j+0, q_0_0);
				}

                //---------- i + 1 ----------
                if(i+1 != row) {
		            __m256d l_1_0 = _mm256_load_pd(tabp+(i+1)*width+j+0);

		            __m256d p_1_0 = _mm256_mul_pd(f1, r0);
		            __m256d q_1_0 = _mm256_sub_pd(l_1_0, p_1_0);

		            _mm256_store_pd(tabp+(i+1)*width+j+0, q_1_0);
				}

                //---------- i + 2 ----------
                if(i+2 != row) {
		            __m256d l_2_0 = _mm256_load_pd(tabp+(i+2)*width+j+0);

		            __m256d p_2_0 = _mm256_mul_pd(f2, r0);
		            __m256d q_2_0 = _mm256_sub_pd(l_2_0, p_2_0);

		            _mm256_store_pd(tabp+(i+2)*width+j+0, q_2_0);
				}

                //---------- i + 3 ----------
                if(i+3 != row) {
		            __m256d l_3_0 = _mm256_load_pd(tabp+(i+3)*width+j+0);

		            __m256d p_3_0 = _mm256_mul_pd(f3, r0);
		            __m256d q_3_0 = _mm256_sub_pd(l_3_0, p_3_0);

		            _mm256_store_pd(tabp+(i+3)*width+j+0, q_3_0);
				}

                //---------- i + 4 ----------
                if(i+4 != row) {
		            __m256d l_4_0 = _mm256_load_pd(tabp+(i+4)*width+j+0);

		            __m256d p_4_0 = _mm256_mul_pd(f4, r0);
		            __m256d q_4_0 = _mm256_sub_pd(l_4_0, p_4_0);

		            _mm256_store_pd(tabp+(i+4)*width+j+0, q_4_0);
				}

                //---------- i + 5 ----------
                if(i+5 != row) {
		            __m256d l_5_0 = _mm256_load_pd(tabp+(i+5)*width+j+0);

		            __m256d p_5_0 = _mm256_mul_pd(f5, r0);
		            __m256d q_5_0 = _mm256_sub_pd(l_5_0, p_5_0);

		            _mm256_store_pd(tabp+(i+5)*width+j+0, q_5_0);
				}

                //---------- i + 6 ----------
                if(i+6 != row) {
		            __m256d l_6_0 = _mm256_load_pd(tabp+(i+6)*width+j+0);

		            __m256d p_6_0 = _mm256_mul_pd(f6, r0);
		            __m256d q_6_0 = _mm256_sub_pd(l_6_0, p_6_0);

		            _mm256_store_pd(tabp+(i+6)*width+j+0, q_6_0);
				}

                //---------- i + 7 ----------
                if(i+7 != row) {
		            __m256d l_7_0 = _mm256_load_pd(tabp+(i+7)*width+j+0);

		            __m256d p_7_0 = _mm256_mul_pd(f7, r0);
		            __m256d q_7_0 = _mm256_sub_pd(l_7_0, p_7_0);

		            _mm256_store_pd(tabp+(i+7)*width+j+0, q_7_0);
				}
            }

            for(int j = aligned_end; j < width; ++j) {
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
            }
        }

        for(int i = m-(m%8); i < m+1; ++i) {
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
