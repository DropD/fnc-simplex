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
//~ #define RDTSC_CYCLES_REQUIRED 1E6
//~ #define RDTSC_CYCLES_REQUIRED 1E7               // warm enough
#define RDTSC_CYCLES_REQUIRED 1E8
//~ #define RDTSC_CYCLES_REQUIRED 1E9               // warm
#include "misc/rdtsc_testing.hpp"


const bool INFO = false;
//~ #define VERBOSE         // in-algorithm info

#include "simplex/baseline.hpp"
#include "simplex/array.hpp"
#include "simplex/block2.hpp"
#include "simplex/block_swap.hpp"
#include "simplex/ssa.hpp"
#include "simplex/sse.hpp"
#include "simplex/avx.hpp"
#include "simplex/block_avx.hpp"
#include "simplex/block_sse.hpp"
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

    Simplex<s_type> s1;
    run(&s1, fname);
    SimplexArray<s_type> s2;
    run(&s2, fname);
    SimplexBlock2<s_type> s3;
    run(&s3, fname);
    SimplexBlock_swap<s_type> s10;
    run(&s10, fname);
    //~ SimplexBlockSSE<s_type> s9;
    //~ run(&s9, fname);
    SimplexBlockAVX<s_type> s11;
    run(&s11, fname);
    SimplexSSA<s_type> s4;
    run(&s4, fname);
    //~ SimplexSSE<s_type> s5;
    //~ run(&s5, fname);
    SimplexAVX<s_type> s6;
    run(&s6, fname);
    SimplexNTA<s_type> s7;
    run(&s7, fname);

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
