/*******************************************************************************
 * 
 * Running Simplex on specified input file
 * 
 ******************************************************************************/

const bool INFO = false;
//~ #define VERBOSE         // in-algorithm info

#include "simplex_baseline.hpp"
typedef Simplex<double> s_type;


//~ #define RDTSC_CYCLES_REQUIRED 1E1           // cold
#define RDTSC_CYCLES_REQUIRED 1E6
//~ #define RDTSC_CYCLES_REQUIRED 1E8           // warm
#include "misc/rdtsc_testing.hpp"


#include <iostream>
#include <string>


using namespace std;


int main(int argc, char ** argv) {

    string fname;

    if(argc != 2) {
        cout << "Requires input file" << endl;
        return 0;
    }
    else
        fname = argv[1];


    s_type s;
    cout << "Running " << s.get_identifier() << endl;

    if(INFO) {
        s.load(fname);
        s.print();
    }

    int n = rdtsc_warmup(&s, fname);
    double cycles = rdtsc_measure(n, &s, fname);    // solve

    vector<double> sol = s.solutions();
    cout << "Optimal value: " << sol[0] << endl;
    if(INFO) {
        cout << "Variables:";
        for(int i = 1; i < sol.size(); ++i)
            cout << "  " << sol[i];
        cout << endl;
    }

    cout << "Memory used: " << s.memusage() << " kB" << endl;
    cout << "RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
    cout << "Memory accesses: " << s.PERFC_MEM << endl;
    cout << "Float add/mul: " << s.PERFC_ADDMUL << endl;
    cout << "Float div: " << s.PERFC_DIV << endl;
    double fpc = (s.PERFC_ADDMUL + s.PERFC_DIV) / cycles;
    cout << "FLOP/C: " << fpc << endl;

    ofstream fp("rdtsc");
    if(fp.is_open()) {
        fp << cycles << endl << fpc << endl;
        fp.close();
    } else
        cout << "Error: unable to write to file rdtsc!" << endl;

    return 0;

}
