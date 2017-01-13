#ifndef __VALUES__
#define __VALUES__

#include <vector>

class Values {

public:

  Values();
  Values(const Values & v);
  void operator=(const Values & other);
  
  void init(const int *n, const double *dx,
	    double (*f)(double, double, double) = 0L);
  void resize(int *n);

  double & operator() (int i,int j,int k) {
    return m_u[i + m_n[0]*j + m_n[0]*m_n[1]*k];
  }
  double operator() (int i,int j,int k) const {
    return m_u[i + m_n[0]*j + m_n[0]*m_n[1]*k];
  }

  
  void swap(Values & other);
  int size(int i) const { return m_n[i]; }
  double dx(int i) const { return m_dx[i]; }
  
private:
  
  std::vector<double> m_u;
  int m_n[3];
  double m_dx[3];
  
};
				   

#endif
