/*******************************************************************************
 * 
 * Running Simplex on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

#define RDTSC_CYCLES_REQUIRED 1E1           // cold
//~ #define RDTSC_CYCLES_REQUIRED 1E6
//~ #define RDTSC_CYCLES_REQUIRED 1E8           // warm
#include "misc/rdtsc_testing.hpp"


const bool INFO = false;
//~ #define VERBOSE         // in-algorithm info

#include "simplex_baseline.hpp"
#include "simplex_array.hpp"



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

    cout << "Memory used: " << s->memusage() << " kB" << endl;
    cout << "RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
    cout << "Memory accesses: " << s->PERFC_MEM << endl;
    cout << "Float add/mul: " << s->PERFC_ADDMUL << endl;
    cout << "Float div: " << s->PERFC_DIV << endl;
    cout << "FLOP/C: " << fpc << endl;
    cout << "Comp Intensity: " << ci << endl;

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


int main(int argc, char ** argv) {

    string fname;

    if(argc != 2) {
        cout << "Requires input file" << endl;
        return 0;
    }
    else
        fname = argv[1];

    remove("rdtsc");

    Simplex<s_type> s;
    run(&s, fname);
    SimplexArray<s_type> s2;
    run(&s2, fname);

    return 0;

}
