/*
 * Timer.cxx
 *
 *  Created on: 5 janv. 2016
 *      Author: tajchman
 */

#include "timer.hpp"

#ifdef __cplusplus
extern "C" {
#endif
  
#if defined(_WIN32)
  
#include <windows.h>
  
  struct sTime {
	  LARGE_INTEGER frequency;        // ticks per second
	  LARGE_INTEGER t;
	};
  
  void initTime(struct sTime *t)
  {
    QueryPerformanceFrequency(&(t->frequency));
    QueryPerformanceCounter(&(t->t));
  }
  
  double elapsedTime(struct sTime *t0)
  {
    LARGE_INTEGER t1; 
    QueryPerformanceCounter(&t1);
    return ((double) (t1.QuadPart - t0->t.QuadPart)) / t0->frequency.QuadPart;
  }
  
#elif defined(__unix)
  
#include <sys/time.h>
#define sTime timeval
  
  
  void initTime(struct sTime *t){
    gettimeofday(t, nullptr);
  }

  double elapsedTime(struct sTime *t0)
  {
    struct sTime t1;
    double dt;
    gettimeofday(&t1, nullptr);
    
    dt = t1.tv_sec - t0->tv_sec + 1e-6*(t1.tv_usec - t0->tv_usec);
    return dt;
  }
  
#endif
  
#ifdef __cplusplus
}
#endif

#include <cstdlib>

Timer::Timer() : m_running(false), m_elapsed(0.0) {
  m_t0 = new sTime;
}

Timer::~Timer() {
  delete (sTime *) m_t0;
}

void Timer::reset()
{
  m_elapsed = 0.0;
  m_running = false;
}

void Timer::start()
{
  if (m_running == true) return;

  initTime((sTime *) m_t0);
  m_running = true;
}

void Timer::stop()
{
  if (m_running == false) return;

  m_elapsed += elapsedTime((sTime *) m_t0);
  m_running = false;
}

double Timer::elapsed() const
{
  return  m_elapsed;

}

#include <iomanip>

void displayTimers(std::ostream & f, const std::vector<Timer> & T)
{
  double T_total;
  size_t i;
  
  for (i=0, T_total=0; i<T.size(); i++)
    T_total += T[i].elapsed();
      
  f << std::endl << std::endl;
  for (i=0; i<T.size(); i++)
    f  << std::setw(5) << i << " " << std::setw(9) << T[i].name()
       << std::setw(11) << std::setprecision(5) << std::fixed
       << T[i].elapsed() << " s"
       << " (" << std::setw(6) << std::setprecision(2) << std::fixed
       << T[i].elapsed()/T_total * 100.0 << " %)"<< std::endl;
  f << "Total "
       << std::setw(20) << std::setprecision(5) << std::fixed
       << T_total << " s"
       << " (" << std::setprecision(2) << std::fixed
       << 100.0 << " %)"<< std::endl << std::endl;
}
