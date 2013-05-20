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
#include "../base_general.hpp"

#ifndef RDTSC_CYCLES_REQUIRED
#define RDTSC_CYCLES_REQUIRED 1E9
#endif


/*
* Cache warm-up for at least RDTSC_CYCLES_REQUIRED cycles,
* recording the required number of executions
*/
template <typename T>
int rdtsc_warmup(SimplexBase<T> * s, std::string fname) {
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

template <typename T>
double rdtsc_measure(int num_runs, SimplexBase<T> * s, std::string fname) {
  double cycles = 0;
  tsc_counter start, end;
  int i;

  for(i = 0; i < num_runs; i++) {
    s->load(fname);
    CPUID(); RDTSC(start);
    s->solve();
    RDTSC(end); CPUID();
    cycles += ((double)COUNTER_DIFF(end, start));
  }
  cycles = cycles / ((double) num_runs);
  return cycles;
}

#endif /* RDTSC_TESTING_H */
