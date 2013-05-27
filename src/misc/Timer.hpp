/*******************************************************************************
 * 
 * Simple timer
 * 
 * The reset routine saves the current time for comparison, the check routine
 * returns the elapsed time. Both can be used arbitrarily often.
 * 
 ******************************************************************************/

#include <sys/time.h>

class Timer {

  private:

  struct timeval a, b;
  struct timezone tz;

  public:

  Timer() {
    reset();
  }

  void reset() {
    gettimeofday(&a,  &tz);
  }

  double check() {
    gettimeofday(&b,  &tz);
    return (b.tv_sec - a.tv_sec) + 1e-6*(b.tv_usec - a.tv_usec);
  }
};
