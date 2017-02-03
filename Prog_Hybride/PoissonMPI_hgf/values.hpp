#ifndef __VALUES__
#define __VALUES__

#include <vector>
#include "Poisson_Parameters.hpp"

class Values {

  public:

    Values();
    Values(const Values & v);
    void operator=(const Values & other);

    void init(const int *n, const double *dx, const double * xmin,
        double (*f)(double, double, double) = 0L);
    void resize(int *n);

    double & operator()(int i, int j, int k) {
      return m_u[i + m_n[0] * j + m_n[0] * m_n[1] * k];
    }
    double operator()(int i, int j, int k) const {
      return m_u[i + m_n[0] * j + m_n[0] * m_n[1] * k];
    }

    void swap(Values & other);
    int size(int i) const {
      return m_n[i];
    }
    double dx(int i) const {
      return m_dx[i];
    }

    void synchronize(Poisson_Parameters & P);

  private:

    std::vector<double> m_u;
    int m_n[3];
    double m_dx[3];
    double m_xmin[3];

    double & coef(int i, int j, int k) {
      return m_u[i + m_n[0] * j + m_n[0] * m_n[1] * k];
    }
};


#endif
