#include "values.hpp"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>

Values::Values() 
{
  m_dx[0] =1;  m_dx[1] = 1;  m_dx[2] = 1;
  m_n[0] = 0;  m_n[1] = 0;  m_n[2] = 0;
}

Values::Values(const Values &other) 
{
  m_dx[0] = other.m_dx[0];  m_dx[1] = other.m_dx[1];  m_dx[2] = other.m_dx[2];
  m_n[0] = other.m_n[0];  m_n[1] = other.m_n[1];  m_n[2] = other.m_n[2];
  m_u = other.m_u;
}

void Values::operator=(const Values & other)
{
  m_dx[0] = other.m_dx[0];  m_dx[1] = other.m_dx[1];  m_dx[2] = other.m_dx[2];
  m_n[0] = other.m_n[0];  m_n[1] = other.m_n[1];  m_n[2] = other.m_n[2];
  m_u = other.m_u;
}

void Values::init(const int *pn, const double *pdx,
		  double (*f)(double, double, double))
{
  int i, nn = 1;
  for (i=0; i<3; i++)
    nn *= (m_n[i] = pn[i]);
  m_dx[0] = pdx[0]; m_dx[1] = pdx[1]; m_dx[2] = pdx[2];

  if (f) {
    m_u.resize(nn);

    int j, k;

    double xmin =  0,
      ymin = 0,
      zmin = 0;
    
    for (i=0; i<m_n[0]; i++)
      for (j=0; j<m_n[1]; j++)
        for (k=0; k<m_n[2]; k++)
          operator()(i,j,k) = f(xmin + i*m_dx[0],
                                ymin + j*m_dx[1],
                                zmin + k*m_dx[2]);
  }
  else
    m_u.assign(nn, 0.0);
}

void Values::resize(int *n)
{
  int i, nn = 1;
  
  for (i=0; i<3; i++)
    nn *= (m_n[i] = n[i]);

  m_u.resize(nn);
}

void Values::swap(Values & other)
{
  m_u.swap(other.m_u);
  int i, tempd;
  double tempx;
  for (i=0; i<3; i++) {
    tempd = m_n[i];
    m_n[i] = other.m_n[i];
    other.m_n[i] = tempd;
    tempx = m_dx[i];
    m_dx[i] = other.m_dx[i];
    other.m_dx[i] = tempx;
  }
}


