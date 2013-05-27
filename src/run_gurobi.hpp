/*******************************************************************************
 * 
 * Running Gurobi on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

#include "simplex/baseline.hpp"
#include "../../gurobi/include/gurobi_c++.h"

using namespace std;
typedef double s_type;
typedef unsigned int uint;

void run_gurobi(string fname, SimplexBase<s_type> * s) {

    cout << "\033[0;31m" << "Running " << "Gurobi" << "\033[0m" << std::endl;

    if(INFO) {
        s->load(fname);
        s->print();
    }
    try {

        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_Threads, 1);
        env.set(GRB_IntParam_Presolve, 0);
        //~ env.set(GRB_IntParam_LogToConsole, 0);

        //~ GRBModel model = GRBModel(env, fname);
        //~ model.optimize();
        int n = rdtsc_warmup(env, fname);
        //~ double cycles = rdtsc_measure(n, env, fname);    // solve
        std::pair<double, double> res = rdtsc_measure(n, env, fname);    // solve
        double cycles = res.first;
        double walltime = res.second;

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
            fp << "gurobi" << ','
               << cycles << ','
               << fpc << ','
               << ci << ','
               << walltime
               << endl;
            fp.close();
        } else
            cout << "Error: unable to write to file rdtsc!" << endl;

    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Error during optimization" << endl;
    }

}
