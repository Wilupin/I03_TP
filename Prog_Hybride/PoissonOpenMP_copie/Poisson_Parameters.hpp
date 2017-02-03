#ifndef __POISSON_PARAMETERS__
#define __POISSON_PARAMETERS__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <memory>
#include "GetPot"

class Poisson_Parameters {
public:

  Poisson_Parameters(int *, char ***);
  ~Poisson_Parameters();
  std::ostream & out() { return *m_out; }
  void info();

  const int *n() const { return m_n; }
  const double *dx() const { return m_dx; }
  int n(int i) const { return m_n[i]; }
  double dx(int i) const { return m_dx[i]; }

  int imin(int i) const { return m_imin[i]; }
  int imax(int i) const { return m_imax[i]; }
  int di(int i) const { return m_di[i]; }
  
  double dt() const { return m_dt; }
  int it() const { return m_it; }
  int nthreads() const { return m_nthreads; };

  int freq() const { return m_freq; }
  std::string resultPath() { return m_path; }
  
  template <typename T>
  T getKey(const char *name, const T& deft) {
    T k = m_prm(name, deft);
    return k;
  }

  std::string getKey(const char *name, const char * deft) {
    return m_prm(name, deft);
  }

  template <typename T>
  void setKey(const char *name, const T& value)
  {
    m_prm.set(name, value);
  }

  bool help();
  bool debug() { return m_debug; }

private:

  std::ostream * m_out;

  int m_n[3];
  double m_dx[3];
  int m_imin[3], m_imax[3], m_di[3];
  
  double m_dt;

  int m_it;
  int m_freq;
  int m_nthreads;

  GetPot m_prm;

  bool m_debug;
  std::string m_command, m_path;
  bool m_help;

};

std::ostream & operator<<(std::ostream &f, const Poisson_Parameters & p);


#endif
