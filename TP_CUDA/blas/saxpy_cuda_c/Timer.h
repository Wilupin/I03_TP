#ifndef _TIMER_H_
#define _TIMER_H_

#include <time.h>
#include <sys/time.h> // for gettimeofday and struct timeval

typedef struct timeval timeval_t;

/**
 * \brief a simple Timer class.
 */
class Timer
{
public:
  /** default constructor, timing starts rightaway */
  Timer();
  
  Timer(double t);
  Timer(Timer const& aTimer);
  ~Timer();
  
  /** start time measure */
  void start();
  
  /** stop time measure and add result to total_time */
  void stop();
  
  /** return elapsed time in seconds (as stored in total_time) */
  double elapsed() const;
  
protected:
  timeval_t start_time;
  
  /** store total accumulated timings */
  double    total_time;
  
}; // class Timer

#endif // _TIMER_H_
