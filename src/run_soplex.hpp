/*******************************************************************************
 * 
 * Running SoPlex on specified input file
 * 
 ******************************************************************************/

#include <iostream>
#include <string>

#include "simplex/baseline.hpp"
#include <soplex.h>

void run_soplex(string fname, SimplexBase<s_type> * s) {

    cout << "\033[0;31m" << "Running " << "soplex" << "\033[0m" << std::endl;

    soplex::SoPlex soplex(soplex::SPxSolver::LEAVE, soplex::SPxSolver::COLUMN);
    soplex.setUtype             ( soplex::SLUFactor::FOREST_TOMLIN );
    soplex.setFeastol           ( DEFAULT_BND_VIOL );
    soplex.setOpttol            ( DEFAULT_BND_VIOL );
    soplex.setIrthreshold       ( DEFAULT_BND_VIOL * 1e-6 );
    soplex.setTerminationTime   ( 10.0 );
    soplex.setTerminationIter   ( -1 );


    std::pair<double, double> res = std::make_pair(0, 0);
    int n = 0;
    try {
        n = rdtsc_warmup(soplex, fname);
        //~ n = 1;
        res = rdtsc_measure(n, soplex, fname);    // solve
    } catch(...) {
        
    }
    double cycles = res.first;
    double walltime = res.second;

    cout << "Optimal value: " << soplex.objValue() << endl;

    double fpc = (cycles==0) ? 0 : (s->PERFC_ADDMUL + s->PERFC_DIV) / cycles;
    double ci = (s->PERFC_ADDMUL+s->PERFC_DIV)/8./s->PERFC_MEM;

    cout << "Wall time: " << walltime << endl;
    cout << "(pseudo) RDTSC cycles: " << cycles << " (avg over " << n << " runs)" << endl;
    cout << "(pseudo) FLOP/C: " << fpc << endl;
    cout << "(pseudo) Op Intensity: " << ci << endl;

    ofstream fp("rdtsc", fstream::app);
    if(fp.is_open()) {
        fp << "soplex" << ','
           << cycles << ','
           << fpc << ','
           << ci << ','
           << walltime
           << endl;
        fp.close();
    } else
        cout << "Error: unable to write to file rdtsc!" << endl;

}
