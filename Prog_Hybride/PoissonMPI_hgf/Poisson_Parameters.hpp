#ifndef __POISSON_PARAMETERS__
#define __POISSON_PARAMETERS__

#include <mpi.h>
#include <omp.h>

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

  MPI_Comm & comm() { return m_comm; }

  int size() const { return m_size; }
  int rank() const { return m_rank; }
  int neighbor(int idim, int j) const { return m_neigh[idim][j]; }
  
  const int *n() const { return m_n; }
  const double *dx() const { return m_dx; }
  const double *xmin() const { return m_xmin; }

  int p(int i) const { return m_p[i]; }
  int n(int i) const { return m_n[i]; }
  int nmax(int i) const { return m_nmax[i]; }
  int p0(int i) const { return m_p0[i]; }
  double dx(int i) const { return m_dx[i]; }

  int imin(int i) const { return m_imin[i]; }
  int imax(int i) const { return m_imax[i]; }
  int di(int i) const { return m_di[i]; }
  
  double dt() const { return m_dt; }
  int it() const { return m_it; }

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
  int m_rank, m_size, m_p[3];

  MPI_Comm m_comm;

  int m_neigh[3][2];
  int m_p0[3], m_n[3], m_nmax[3];
  double m_dx[3], m_xmin[3];
  int m_imin[3], m_imax[3], m_di[3];
  
  double m_dt;

  int m_it;
  int m_freq;

  GetPot m_prm;

  bool m_debug;
  std::string m_command, m_path;
  bool m_help;

};

std::ostream & operator<<(std::ostream &f, const Poisson_Parameters & p);


#endif
