/*******************************************************************************
 * 
 * Running Gurobi on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

#include "simplex/baseline.hpp"
#include <glpk.h>

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
    std::pair<double, double> res = rdtsc_measure(n, lp, &parm, fname);
    double cycles = res.first;
    double walltime = res.second;
    double obj = glp_get_obj_val(lp);
    cout << "Optimal value: " << obj << endl;
    glp_delete_prob(lp);

    double fpc = (s->PERFC_ADDMUL + s->PERFC_DIV) / cycles;
    double ci = (s->PERFC_ADDMUL+s->PERFC_DIV)/8./s->PERFC_MEM;

    cout << "Wall time: " << walltime << endl;
    cout << "(pseudo) RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
    cout << "(pseudo) FLOP/C: " << fpc << endl;
    cout << "(pseudo) Op Intensity: " << ci << endl;

    ofstream fp("rdtsc", fstream::app);
    if(fp.is_open()) {
        fp << "glpk" << ','
           << cycles << ','
           << fpc << ','
           << ci << ','
           << walltime
           << endl;
        fp.close();
    } else
        cout << "Error: unable to write to file rdtsc!" << endl;
}
