/*******************************************************************************
 * 
 * Running Simplex on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

//~ #define NO_GLPK
//~ #define NO_GUROBI
//~ #define NO_SOPLEX

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
#include "simplex/block2x4.hpp"
//** start blockswap
#include "simplex/block1x1_swap.hpp"
#include "simplex/block1x2_swap.hpp"
#include "simplex/block1x4_swap.hpp"
#include "simplex/block1x8_swap.hpp"
#include "simplex/block1x16_swap.hpp"
#include "simplex/block2x1_swap.hpp"
#include "simplex/block2x2_swap.hpp"
#include "simplex/block2x4_swap.hpp"
#include "simplex/block2x8_swap.hpp"
#include "simplex/block2x16_swap.hpp"
#include "simplex/block4x1_swap.hpp"
#include "simplex/block4x2_swap.hpp"
#include "simplex/block4x4_swap.hpp"
#include "simplex/block4x8_swap.hpp"
#include "simplex/block4x16_swap.hpp"
#include "simplex/block8x1_swap.hpp"
#include "simplex/block8x2_swap.hpp"
#include "simplex/block8x4_swap.hpp"
#include "simplex/block8x8_swap.hpp"
#include "simplex/block8x16_swap.hpp"
#include "simplex/block16x1_swap.hpp"
#include "simplex/block16x2_swap.hpp"
#include "simplex/block16x4_swap.hpp"
#include "simplex/block16x8_swap.hpp"
#include "simplex/block16x16_swap.hpp"
//** end blockswap
//** start blockswap_avx
#include "simplex/block1x4_swap_avx.hpp"
#include "simplex/block1x8_swap_avx.hpp"
#include "simplex/block1x16_swap_avx.hpp"
#include "simplex/block2x4_swap_avx.hpp"
#include "simplex/block2x8_swap_avx.hpp"
#include "simplex/block2x16_swap_avx.hpp"
#include "simplex/block4x4_swap_avx.hpp"
#include "simplex/block4x8_swap_avx.hpp"
#include "simplex/block4x16_swap_avx.hpp"
#include "simplex/block8x4_swap_avx.hpp"
#include "simplex/block8x8_swap_avx.hpp"
#include "simplex/block8x16_swap_avx.hpp"
#include "simplex/block16x4_swap_avx.hpp"
#include "simplex/block16x8_swap_avx.hpp"
#include "simplex/block16x16_swap_avx.hpp"
//** end blockswap_avx
#include "simplex/ssa.hpp"
#include "simplex/sse.hpp"
#include "simplex/avx.hpp"
#include "simplex/block2x4_sse.hpp"
#include "simplex/block2x4_swap_avx.hpp"
#include "simplex/block4x4_swap_avx.hpp"
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
    {
    //~ Simplex_array<s_type> s;
    //~ run(&s, fname);
    }{
    Simplex_block2x4<s_type> s;
    run(&s, fname);
    }{
    //** start blockswap
    Simplex_block1x1_swap<s_type> s; 
    run(&s, fname);
    }{
    Simplex_block1x2_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block1x4_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block1x8_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block1x16_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x1_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x2_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x4_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x8_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x16_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x1_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x2_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x4_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x8_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x16_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x1_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x2_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x4_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x8_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x16_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x1_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x2_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x4_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x8_swap<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x16_swap<s_type> s;
    run(&s, fname);
    }{ 
    //** end blockswap
    //** start blockswap_avx_avx
    Simplex_block1x4_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block1x8_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block1x16_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x4_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x8_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block2x16_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x4_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x8_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block4x16_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x4_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x8_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block8x16_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x4_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x8_swap_avx<s_type> s;
    run(&s, fname);
    }{
    Simplex_block16x16_swap_avx<s_type> s;
    run(&s, fname);
    }{ 
    //** end blockswap_avx_avx
    //~ Simplex_block2x4_sse<s_type> s;
    //~ run(&s, fname);
    }{
    //~ Simplex_block2x4_sse<s_type> s;
    //~ run(&s, fname);
    }{
    //~ Simplex_block2x4_swap_avx<s_type> s;
    //~ run(&s, fname);
    }{
    //~ Simplex_block4x4_swap_avx<s_type> s;
    //~ run(&s, fname);
    }{
    //~ Simplex_block8_swap_avx<s_type> s;
    //~ run(&s, fname);
    }{
    //~ Simplex_ssa<s_type> s
    //~ run(&s fname);
    }{
    //~ Simplex_sse<s_type> s
    //~ run(&s fname);
    }{
    //~ Simplex_avx<s_type> s
    //~ run(&s fname);
    }{
    //~ Simplex_nta<s_type> s
    //~ run(&s fname);
    }

    //replace file extension
    string lname = fname.substr(0, fname.length()-3);
    lname.append("lp");

#ifndef NO_GLPK
    run_glpk(lname, &s1);
#endif // NO_GLPK
#ifndef NO_GUROBI
    run_gurobi(lname, &s1);
#endif // NO_GUROBI
#ifndef NO_SOPLEX
    run_soplex(lname, &s1);
#endif // NO_SOPLEX

    return 0;

}
