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
#include "simplex/block2.hpp"
#include "simplex/block2_swap.hpp"
#include "simplex/ssa.hpp"
#include "simplex/sse.hpp"
#include "simplex/avx.hpp"
#include "simplex/block2_sse.hpp"
#include "simplex/block2_swap_avx.hpp"
#include "simplex/block4_swap_avx.hpp"
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

    Simplex_baseline<s_type> s1;
    run(&s1, fname);
    //~ Simplex_array<s_type> s2;
    //~ run(&s2, fname);
    //~ Simplex_block2<s_type> s3;
    //~ run(&s3, fname);
    //~ Simplex_block2_swap<s_type> s10;
    //~ run(&s10, fname);
    //~ Simplex_block2_sse<s_type> s9;
    //~ run(&s9, fname);
    //~ Simplex_block2_swap_avx<s_type> s11;
    //~ run(&s11, fname);
    //~ Simplex_block4_swap_avx<s_type> s12;
    //~ run(&s12, fname);
    //~ Simplex_block8_swap_avx<s_type> s13;
    //~ run(&s13, fname);
    //~ Simplex_ssa<s_type> s4;
    //~ run(&s4, fname);
    //~ Simplex_sse<s_type> s5;
    //~ run(&s5, fname);
    //~ Simplex_avx<s_type> s6;
    //~ run(&s6, fname);
    //~ Simplex_nta<s_type> s7;
    //~ run(&s7, fname);

    //replace file extension
    string lname = fname.substr(0, fname.length()-3);
    lname.append("lp");

//~ #ifndef NO_GLPK
    //~ run_glpk(lname, &s1);
//~ #endif // NO_GLPK
#ifndef NO_GUROBI
    run_gurobi(lname, &s1);
#endif // NO_GUROBI
#ifndef NO_SOPLEX
    run_soplex(lname, &s1);
#endif // NO_SOPLEX

    return 0;

}
