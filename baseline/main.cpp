//#define VERBOSE
#include "general.hpp"

#include <iostream>
#include <string>
#include <glpk.h>
#include <fstream>
#include <vector>

using namespace std;

vector<double> time_impl(string fname) {
    Simplex<double> s;
    s.load(fname);

    tsc_counter start, end;
    CPUID(); RDTSC(start);
    s.solve();
    RDTSC(end); CPUID();
    double cycles = (double)(COUNTER_DIFF(end, start));

    vector<double> sol = s.solutions();

    vector<double> res;
    res.push_back(cycles);
    res.push_back(sol[0]);

    return res;
}

vector<double> time_glpk(string lname) {
    glp_term_out(GLP_OFF);
    glp_prob *lp;
    glp_smcp parm;
    glp_init_smcp(&parm);
    parm.msg_lev = GLP_MSG_OFF;
    lp = glp_create_prob();
    glp_read_lp(lp, NULL, lname.c_str());

    tsc_counter start, end;
    CPUID(); RDTSC(start);
    glp_simplex(lp, NULL);
    RDTSC(end); CPUID();
    double cycles = (double)(COUNTER_DIFF(end, start));
    double obj = glp_get_obj_val(lp);
    glp_delete_prob(lp);

    vector<double> res;
    res.push_back(cycles);
    res.push_back(obj);

    return res;
}

int main(int argc, char ** argv) {
    string fname, lname;

    if(argc != 3) {
        cout << "Requires dlp and lp input files" << endl;
        return 0;
    }
    else {
        fname = argv[1];
        lname = argv[2];
    }

    vector<double> impl_res, glpk_res;
    impl_res = time_impl(fname);
    glpk_res = time_glpk(lname);

    //double tol = 1e-3;

    double delta = impl_res[1] - glpk_res[1];
    //if(delta < tol && delta > -tol) {
    //    cout << impl_res[0] << " " << glpk_res[0] << endl;
    //}
    //else {
    //    cout << "objective value disagreement" << endl;
    //    return -1;
    //}
    cout << impl_res[0] << " " << glpk_res[0] << " " << delta << endl;

    return 0;

}
