/*
 * Timer.hxx
 *
 *  Created on: 5 janv. 2016
 *      Author: tajchman
 */

#ifndef TIMER_HXX_
#define TIMER_HXX_

#include <string>

class Timer {
public:
  Timer();
  ~Timer();

  const std::string & name() const  { return m_name; }
  void name(const std::string & name) { m_name = name; }
  
  void reset();
  void start();
  void stop();
  double elapsed() const;

private:
  std::string m_name;
  void *m_t0;
  bool m_running;
  double m_elapsed;
};

#include <vector>
#include <iostream>

void displayTimers(std::ostream & f, const std::vector<Timer> & T);

#endif /* TIMER_HXX_ */
