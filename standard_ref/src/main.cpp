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
#define RDTSC_CYCLES_REQUIRED 1E8           // warm
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
    double cycles = rdtsc_measure(222, &s, fname);    // solve

    if(INFO) {
        vector<double> sol = s.solutions();
        cout << "Optimal value: " << sol[0] << endl;
        cout << "Variables:";
        for(int i = 1; i < sol.size(); ++i)
            cout << "  " << sol[i];
        cout << endl;
    }

    cout << "RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
    cout << "Memory accesses: " << s.PERFC_MEM << endl;
    cout << "Float add/mul: " << s.PERFC_ADDMUL << endl;
    cout << "Float div: " << s.PERFC_DIV << endl;
    cout << "FLOP/C: " << (s.PERFC_ADDMUL + s.PERFC_DIV) / cycles << endl;

    return 0;

}
