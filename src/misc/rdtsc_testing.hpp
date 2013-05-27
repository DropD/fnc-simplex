/*******************************************************************************
 * 
 * Measurement functions for RDTSC
 * 
 * Simplified version for testing the Simplex algorithm.
 * 
 ******************************************************************************/

#ifndef RDTSC_TESTING_H
#define RDTSC_TESTING_H

#include "rdtsc.h"  /* actual timing macros */
#include "Timer.hpp"

#include <iostream>

#ifndef RDTSC_CYCLES_REQUIRED
#define RDTSC_CYCLES_REQUIRED 1
#endif

/*
* Cache warm-up for at least RDTSC_CYCLES_REQUIRED cycles,
* recording the required number of executions
*/

#include "../simplex/Simplex.hpp"
template <typename T>
int rdtsc_warmup(SimplexBase<T> * s, std::string fname) {
  if(RDTSC_CYCLES_REQUIRED <= 1) return 1;
  double cycles = 0;
  tsc_counter start, end;
  int i;
  int num_runs = 1;
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  while(1) {
    for(i = 0; i < num_runs; i++) {
      if( s->restore_tableau() == false )
        s->load(fname);
      CPUID(); RDTSC(start);
      s->solve();
      RDTSC(end); CPUID();
      cycles += ((double)COUNTER_DIFF(end, start));
    }
    if(cycles >= RDTSC_CYCLES_REQUIRED) break;
    num_runs *= 2;
  }
  return num_runs;
}

template <typename T>
std::pair<double, double> rdtsc_measure(int num_runs, SimplexBase<T> * s, std::string fname) {
  double cycles = 0;
  double walltime = 0;
  Timer tim;
  tsc_counter start, end;
  int i;

  for(i = 0; i < num_runs; i++) {
    if( s->restore_tableau() == false )
      s->load(fname);
    tim.reset();
    CPUID(); RDTSC(start);
    s->solve();
    RDTSC(end); CPUID();
    walltime += tim.check();
    cycles += ((double)COUNTER_DIFF(end, start));
  }
  cycles = cycles / ((double) num_runs);
  walltime = walltime / ((double) num_runs);
  return std::make_pair(cycles, walltime);
}


#ifndef NO_GLPK
#include <glpk.h>
int rdtsc_warmup(glp_prob * lp, glp_smcp * parm, std::string fname) {
  if(RDTSC_CYCLES_REQUIRED <= 1) return 1;
  double cycles = 0;
  tsc_counter start, end;
  int i;
  int num_runs = 1;
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  while(1) {
    for(i = 0; i < num_runs; i++) {
      glp_read_lp(lp, NULL, fname.c_str());
      CPUID(); RDTSC(start);
      glp_simplex(lp, parm);
      RDTSC(end); CPUID();
      cycles += ((double)COUNTER_DIFF(end, start));
    }
    if(cycles >= RDTSC_CYCLES_REQUIRED) break;
    num_runs *= 2;
  }
  return num_runs;
}
std::pair<double, double> rdtsc_measure(int num_runs, glp_prob * lp, glp_smcp * parm, std::string fname) {
  double cycles = 0;
  double walltime = 0;
  Timer tim;
  tsc_counter start, end;
  int i;

  for(i = 0; i < num_runs; i++) {
    glp_read_lp(lp, NULL, fname.c_str());
    tim.reset();
    CPUID(); RDTSC(start);
    glp_simplex(lp, parm);
    RDTSC(end); CPUID();
    walltime += tim.check();
    cycles += ((double)COUNTER_DIFF(end, start));
  }
  cycles = cycles / ((double) num_runs);
  walltime = walltime / ((double) num_runs);
  return std::make_pair(cycles, walltime);
}
#endif // GLPK

#ifndef NO_GUROBI
#include <gurobi_c++.h>
int rdtsc_warmup(GRBEnv & env, std::string fname) {
  if(RDTSC_CYCLES_REQUIRED <= 1) return 1;
  double cycles = 0;
  tsc_counter start, end;
  int i;
  int num_runs = 1;
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  while(1) {
    for(i = 0; i < num_runs; i++) {
      GRBModel model = GRBModel(env, fname);
      CPUID(); RDTSC(start);
      model.optimize();
      RDTSC(end); CPUID();
      cycles += ((double)COUNTER_DIFF(end, start));
    }
    if(cycles >= RDTSC_CYCLES_REQUIRED) break;
    num_runs *= 2;
  }
  return num_runs;
}
std::pair<double, double> rdtsc_measure(int num_runs, GRBEnv & env, std::string fname) {
  double cycles = 0;
  double walltime = 0;
  Timer tim;
  tsc_counter start, end;
  int i;

  for(i = 0; i < num_runs; i++) {
    GRBModel model = GRBModel(env, fname);
    tim.reset();
    CPUID(); RDTSC(start);
    model.optimize();
    RDTSC(end); CPUID();
    walltime += tim.check();
    cycles += ((double)COUNTER_DIFF(end, start));
  }
  cycles = cycles / ((double) num_runs);
  walltime = walltime / ((double) num_runs);
  return std::make_pair(cycles, walltime);
}
#endif // NO_GUROBI


#ifndef NO_SOPLEX
#include <soplex.h>

int rdtsc_warmup(soplex::SoPlex& s, std::string fname) {
  if(RDTSC_CYCLES_REQUIRED <= 1) return 1;
  double cycles = 0;
  tsc_counter start, end;
  int i;
  int num_runs = 1;
  soplex::NameSet rownames, colnames;
  soplex::SPxSolver::Status stat;
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  CPUID(); RDTSC(start); RDTSC(end);
  while(1) {
    for(i = 0; i < num_runs; i++) {
      s.readFile(fname.c_str(), &rownames, &colnames, 0);
      CPUID(); RDTSC(start);
      stat = s.solve();
      RDTSC(end); CPUID();
      cycles += ((double)COUNTER_DIFF(end, start));
    }
    if(cycles >= RDTSC_CYCLES_REQUIRED) break;
    num_runs *= 2;
  }
  if( stat != soplex::SPxSolver::OPTIMAL )
    std::cout << "Warning: soplex has no optimal solution!" << std::endl;
  return num_runs;
}
std::pair<double, double> rdtsc_measure(int num_runs, soplex::SoPlex& s, std::string fname) {
  double cycles = 0;
  double walltime = 0;
  Timer tim;
  tsc_counter start, end;
  int i;
  soplex::NameSet rownames, colnames;
  soplex::SPxSolver::Status stat = soplex::SPxSolver::OPTIMAL;

  for(i = 0; i < num_runs; i++) {
    s.readFile(fname.c_str(), &rownames, &colnames, 0);
    tim.reset();
    CPUID(); RDTSC(start);
    stat = s.solve();
    RDTSC(end); CPUID();
    walltime += tim.check();
    cycles += ((double)COUNTER_DIFF(end, start));
  }
  cycles = cycles / ((double) num_runs);
  walltime = walltime / ((double) num_runs);
  if( stat != soplex::SPxSolver::OPTIMAL )
    std::cout << "Warning: soplex has no optimal solution!" << std::endl;
  return std::make_pair(cycles, walltime);
}
#endif // NO_SOPLEX


#endif /* RDTSC_TESTING_H */
