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

#include "../simplex/Simplex.hpp"
#include <glpk.h>
#include <gurobi_c++.h>

#include <iostream>

#ifndef RDTSC_CYCLES_REQUIRED
#define RDTSC_CYCLES_REQUIRED 1E9
#endif


/*
* Cache warm-up for at least RDTSC_CYCLES_REQUIRED cycles,
* recording the required number of executions
*/
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





template <typename T>
//~ double rdtsc_measure(int num_runs, SimplexBase<T> * s, std::string fname) {
std::pair<double, double> rdtsc_measure(int num_runs, SimplexBase<T> * s, std::string fname) {
  double cycles = 0;
  double walltime = 0;
  Timer tim;
  tsc_counter start, end;
  int i;

  for(i = 0; i < num_runs; i++) {
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
  //~ return cycles;
  return std::make_pair(cycles, walltime);
}

//~ double rdtsc_measure(int num_runs, glp_prob * lp, glp_smcp * parm, std::string fname) {
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
  //~ return cycles;
  return std::make_pair(cycles, walltime);
}

//~ double rdtsc_measure(int num_runs, GRBEnv & env, std::string fname) {
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
  //~ return cycles;
  return std::make_pair(cycles, walltime);
}

#endif /* RDTSC_TESTING_H */
