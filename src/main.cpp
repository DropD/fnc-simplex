/*******************************************************************************
 * 
 * Running Simplex on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

#define NO_GLPK
#define NO_GUROBI
#define NO_SOPLEX

//~ #define RDTSC_CYCLES_REQUIRED 0                 // cold
//~ #define RDTSC_CYCLES_REQUIRED 1E3
#define RDTSC_CYCLES_REQUIRED 1E7               // warm enough
//~ #define RDTSC_CYCLES_REQUIRED 1E8
//~ #define RDTSC_CYCLES_REQUIRED 1E9               // warm
#include "misc/rdtsc_testing.hpp"


const bool INFO = false;
//~ #define VERBOSE         // in-algorithm info

#include "simplex/baseline.hpp"
#include "simplex/array.hpp"
#include "simplex/block.hpp"
#include "simplex/block_swap.hpp"
#include "simplex/block_avx.hpp"
#include "simplex/block_swap_avx.hpp"

#include "simplex/ssa.hpp"
#include "simplex/sse.hpp"
#include "simplex/avx.hpp"
#include "simplex/block2x4_sse.hpp"
#include "simplex/block2x4_swap_nta.hpp"
#include "simplex/block2x4_swap_nta_pf.hpp"
#include "simplex/block8_swap_avx.hpp"
#include "simplex/nta.hpp"



using namespace std;
typedef double s_type;
typedef unsigned int uint;

#include "run_simplex.hpp"

#ifndef NO_GLPK
#include "run_glpk.hpp"
#endif // NO_GLPK

#ifndef NO_GUROBI
#include "run_gurobi.hpp"
#endif // NO_GUROBI

#ifndef NO_SOPLEX
#include "run_soplex.hpp"
#endif // NO_SOPLEX


#define SIMPLEX_IMPL(id)                                                       \
    {                                                                          \
        Simplex_ ## id   <s_type> s ## id;                                     \
        run(&s ## id , fname);                                                 \
    }


int main(int argc, char ** argv) {

    string fname;

    if(argc != 2) {
        cout << "Requires input file" << endl;
        return 0;
    }
    else
        fname = argv[1];

    remove("rdtsc");

    Simplex_baseline<s_type> s1;  // this one shouldn't go out of scope
    run(&s1, fname);

    //~ SIMPLEX_IMPL(array)

    //~ SIMPLEX_IMPL(block1x1_swap)
    //~ SIMPLEX_IMPL(block1x2_swap)
    //~ SIMPLEX_IMPL(block1x4_swap)
    //~ SIMPLEX_IMPL(block1x8_swap)
    //~ SIMPLEX_IMPL(block1x16_swap)
    //~ SIMPLEX_IMPL(block2x1_swap)
    //~ SIMPLEX_IMPL(block2x2_swap)
    SIMPLEX_IMPL(block2x4_swap)
    //~ SIMPLEX_IMPL(block2x8_swap)
    //~ SIMPLEX_IMPL(block2x16_swap)
    //~ SIMPLEX_IMPL(block4x1_swap)
    //~ SIMPLEX_IMPL(block4x2_swap)
    //~ SIMPLEX_IMPL(block4x4_swap)
    //~ SIMPLEX_IMPL(block4x8_swap)
    //~ SIMPLEX_IMPL(block4x16_swap)
    //~ SIMPLEX_IMPL(block8x1_swap)
    //~ SIMPLEX_IMPL(block8x2_swap)
    //~ SIMPLEX_IMPL(block8x4_swap)
    //~ SIMPLEX_IMPL(block8x8_swap)
    //~ SIMPLEX_IMPL(block8x16_swap)
    //~ SIMPLEX_IMPL(block16x1_swap)
    //~ SIMPLEX_IMPL(block16x2_swap)
    //~ SIMPLEX_IMPL(block16x4_swap)
    //~ SIMPLEX_IMPL(block16x8_swap)
    //~ SIMPLEX_IMPL(block16x16_swap)

    //~ SIMPLEX_IMPL(block1x4_swap_avx)
    //~ SIMPLEX_IMPL(block1x8_swap_avx)
    //~ SIMPLEX_IMPL(block1x16_swap_avx)
    //~ SIMPLEX_IMPL(block2x4_swap_avx)
    //~ SIMPLEX_IMPL(block2x8_swap_avx)
    //~ SIMPLEX_IMPL(block2x16_swap_avx)
    //~ SIMPLEX_IMPL(block4x4_swap_avx)
    //~ SIMPLEX_IMPL(block4x8_swap_avx)
    //~ SIMPLEX_IMPL(block4x16_swap_avx)
    //~ SIMPLEX_IMPL(block8x4_swap_avx)
    //~ SIMPLEX_IMPL(block8x8_swap_avx)
    //~ SIMPLEX_IMPL(block8x16_swap_avx)
    //~ SIMPLEX_IMPL(block16x4_swap_avx)
    //~ SIMPLEX_IMPL(block16x8_swap_avx)
    //~ SIMPLEX_IMPL(block16x16_swap_avx)

    //~ SIMPLEX_IMPL(block2x4_sse)
    SIMPLEX_IMPL(block2x4_swap_nta)
    SIMPLEX_IMPL(block2x4_swap_nta_pf)
    //~ SIMPLEX_IMPL(ssa)
    //~ SIMPLEX_IMPL(sse)
    //~ SIMPLEX_IMPL(avx)
    //~ SIMPLEX_IMPL(nta)

    //replace file extension
    string lname = fname.substr(0, fname.length()-3);
    lname.append("lp");

//~ #ifndef NO_GLPK
    //~ run_glpk(lname, &s1);
//~ #endif // NO_GLPK
//~ #ifndef NO_GUROBI
    //~ run_gurobi(lname, &s1);
//~ #endif // NO_GUROBI
//~ #ifndef NO_SOPLEX
    //~ run_soplex(lname, &s1);
//~ #endif // NO_SOPLEX

    return 0;

}
