/*******************************************************************************
 * 
 * Running Simplex on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

//~ #define RDTSC_CYCLES_REQUIRED 0                 // cold
#define RDTSC_CYCLES_REQUIRED 1E6
//~ #define RDTSC_CYCLES_REQUIRED 1E7               // warm enough
//~ #define RDTSC_CYCLES_REQUIRED 1E9               // warm
#include "misc/rdtsc_testing.hpp"


const bool INFO = false;
//~ #define VERBOSE         // in-algorithm info

#include "simplex/baseline.hpp"
#include "simplex/array.hpp"
#include "simplex/block_2.hpp"
#include "simplex/ssa.hpp"
#include "simplex/sse.hpp"
#include "simplex/avx.hpp"
#include "simplex/nta.hpp"
#include <glpk.h>


using namespace std;
typedef double s_type;
typedef unsigned int uint;

void run(SimplexBase<s_type> * s, string fname) {

    string id = s->get_identifier();
    cout << "\033[0;31m" << "Running " << id  << "\033[0m" << std::endl;

    if(INFO) {
        s->load(fname);
        s->print();
    }

    int n = rdtsc_warmup(s, fname);
    double cycles = rdtsc_measure(n, s, fname);    // solve

    vector<double> sol = s->solutions();
    cout << "Optimal value: " << sol[0] << endl;
    if(INFO) {
        cout << "Variables:";
        for(uint i = 1; i < sol.size(); ++i)
            cout << "  " << sol[i];
        cout << endl;
    }

    double fpc = (s->PERFC_ADDMUL + s->PERFC_DIV) / cycles;
    double ci = (s->PERFC_ADDMUL+s->PERFC_DIV)/8./s->PERFC_MEM;

    cout << "Iterations: " << s->get_iter() << endl;
    cout << "Memory used: " << s->memusage() << " kB" << endl;
    cout << "RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
    cout << "Memory accesses: " << s->PERFC_MEM
         << " (theory: " << s->get_iter() * s->get_tabn() << ")" << endl;
    cout << "Float add/mul: " << s->PERFC_ADDMUL << endl;
    cout << "Float div: " << s->PERFC_DIV << endl;
    cout << "FLOP/C: " << fpc << endl;
    cout << "Op Intensity: " << ci << endl;

    ofstream fp("rdtsc", fstream::app);
    if(fp.is_open()) {
        fp << id << ','
           << cycles << ','
           << fpc << ','
           << ci << endl;
        fp.close();
    } else
        cout << "Error: unable to write to file rdtsc!" << endl;

}

void run_glpk(string fname, SimplexBase<s_type> * s) {
    cout << "\033[0;31m" << "Running glpk" << "\033[0m" << std::endl;

    glp_term_out(GLP_OFF);
    glp_prob *lp;
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    parm.meth = GLP_PRIMAL;
    lp = glp_create_prob();
    int n = rdtsc_warmup(lp, &parm, fname);
    double cycles = rdtsc_measure(n, lp, &parm, fname);
    double obj = glp_get_obj_val(lp);
    cout << "Optimal value: " << obj << endl;
    glp_delete_prob(lp);

    double fpc = (s->PERFC_ADDMUL + s->PERFC_DIV) / cycles;
    double ci = (s->PERFC_ADDMUL+s->PERFC_DIV)/8./s->PERFC_MEM;

    cout << "Memory used: " << s->memusage() << " kB" << endl;
    cout << "RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
    cout << "Memory accesses: " << s->PERFC_MEM
         << " (theory: " << s->get_iter() * s->get_tabn() << ")" << endl;
    cout << "Float add/mul: " << s->PERFC_ADDMUL << endl;
    cout << "Float div: " << s->PERFC_DIV << endl;
    cout << "FLOP/C: " << fpc << endl;
    cout << "Op Intensity: " << ci << endl;

    ofstream fp("rdtsc", fstream::app);
    if(fp.is_open()) {
        fp << "glpk" << ','
           << cycles << ','
           << fpc << ','
           << ci << endl;
        fp.close();
    } else
        cout << "Error: unable to write to file rdtsc!" << endl;
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

    Simplex<s_type> s1;
    run(&s1, fname);
    SimplexArray<s_type> s2;
    run(&s2, fname);
    SimplexBlock_2<s_type> s3;
    run(&s3, fname);
    SimplexSSA<s_type> s4;
    run(&s4, fname);
    SimplexSSE<s_type> s5;
    run(&s5, fname);
    SimplexAVX<s_type> s6;
    run(&s6, fname);
    SimplexNTA<s_type> s7;
    run(&s7, fname);

    //replace file extension
    string lname = fname.substr(0, fname.length()-3);
    lname.append("lp");

    run_glpk(lname, &s1);

    return 0;

}
