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

        int n = rdtsc_warmup(env, fname);
        std::pair<double, double> res = rdtsc_measure(n, env, fname);    // solve
        double cycles = res.first;
        double walltime = res.second;

        double fpc = (s->PERFC_ADDMUL + s->PERFC_DIV) / cycles;
        double ci = (s->PERFC_ADDMUL+s->PERFC_DIV)/8./s->PERFC_MEM;

        cout << "Wall time: " << walltime << endl;
        cout << "(pseudo) RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
        cout << "(pseudo) FLOP/C: " << fpc << endl;
        cout << "(pseudo) Op Intensity: " << ci << endl;

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
