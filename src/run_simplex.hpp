/*******************************************************************************
 * 
 * Running Gurobi on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

#include "simplex/baseline.hpp"

void run(SimplexBase<s_type> * s, string fname) {

    string id = s->get_identifier();
    cout << "\033[0;31m" << "Running " << id  << "\033[0m" << std::endl;

    if(INFO) {
        s->load(fname);
        s->print();
    }

    int n = rdtsc_warmup(s, fname);
    std::pair<double, double> res = rdtsc_measure(n, s, fname);    // solve
    double cycles = res.first;
    double walltime = res.second;

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
    cout << "Wall time: " << walltime << endl;
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
           << ci << ','
           << walltime
           << endl;
        fp.close();
    } else
        cout << "Error: unable to write to file rdtsc!" << endl;

}
