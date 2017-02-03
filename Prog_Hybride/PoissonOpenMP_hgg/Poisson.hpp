#ifndef Poisson_HPP
#define Poisson_HPP

#include <fstream>
#include <vector>
#include "timer.hpp"
#include "values.hpp"
#include "Poisson_Parameters.hpp"

class Poisson {
  
public:
  
  Poisson();
  virtual ~Poisson();

  void init(Poisson_Parameters & P);

  bool solve(unsigned int nsteps = 1);

  double variation();
  double present();
  
  void terminate();
  
  void setInput(const Values &);
  const Values & getOutput();
  
  size_t getDomainSize(int dim) const;
  
  void save(const char *fName);

  Timer & timer(int i) { return T[i]; }
  std::vector<double> timers();

protected:

  std::string codeName;

  size_t kStep;
  Values m_u, m_v;
  int m_n[3], m_dx[3], m_di[3], m_imin[3], m_imax[3];
  double m_duv;
  double m_t;
  
  std::vector<Timer> T;
  bool debug;
  std::ofstream fDebug;

  Poisson_Parameters *m_P; 

};

#endif
