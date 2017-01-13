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

void Values::init(const int *pn, const double *pdx, const double * pxmin,
		  double (*f)(double, double, double))
{
  int i, nn = 1;
  for (i=0; i<3; i++)
    nn *= (m_n[i] = pn[i]);
  m_dx[0] = pdx[0]; m_dx[1] = pdx[1]; m_dx[2] = pdx[2];
  m_xmin[0] = pxmin[0];  m_xmin[1] = pxmin[1];  m_xmin[2] = pxmin[2];

  if (f) {
    m_u.resize(nn);

    int j, k;

    double xmin = m_xmin[0],
      ymin =  m_xmin[1],
      zmin =  m_xmin[2];

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
    tempx = m_xmin[i];
    m_xmin[i] = other.m_xmin[i];
    other.m_xmin[i] = tempx;
  }
}

#include <mpi.h>

void Values::synchronize(Poisson_Parameters & P) {

  int i, j, k, l, n = m_n[0], m = m_n[1], p = m_n[2];

  MPI_Status status;

  if (P.neighbor(0, 0) != MPI_PROC_NULL) {
    std::vector<double> buf_in(m * p);
    std::vector<double> buf_out(m * p);
    for (l = 0, k = 0; k < p; k++)
      for (j = 0; j < m; j++)
        buf_in[l++] = coef(1, j, k);

    MPI_Sendrecv(buf_in.data(), m * p, MPI_DOUBLE, P.neighbor(0, 0), 0,
        buf_out.data(), m * p, MPI_DOUBLE, P.neighbor(0, 0), 0, P.comm(),
        &status);

    for (l = 0, k = 0; k < p; k++)
      for (j = 0; j < m; j++) {
        coef(0, j, k) = buf_out[l++];
      }
  }

  if (P.neighbor(0, 1) != MPI_PROC_NULL) {
    std::vector<double> buf_in(m * p);
    std::vector<double> buf_out(m * p);
    for (l = 0, k = 0; k < p; k++)
      for (j = 0; j < m; j++) {
        buf_in[l++] = coef(n - 2, j, k);
      }

    MPI_Sendrecv(buf_in.data(), m * p, MPI_DOUBLE, P.neighbor(0, 1), 0,
        buf_out.data(), m * p, MPI_DOUBLE, P.neighbor(0, 1), 0, P.comm(),
        &status);

    for (l = 0, k = 0; k < p; k++)
      for (j = 0; j < m; j++)
        coef(n - 1, j, k) = buf_out[l++];
  }

  if (P.neighbor(1, 0) != MPI_PROC_NULL) {
    std::vector<double> buf_in(n * p);
    std::vector<double> buf_out(n * p);
    for (l = 0, k = 0; k < p; k++)
      for (i = 0; i < n; i++)
        buf_in[l++] = coef(i, 1, k);

    MPI_Sendrecv(buf_in.data(), n * p, MPI_DOUBLE, P.neighbor(1, 0), 0,
        buf_out.data(), n * p, MPI_DOUBLE, P.neighbor(1, 0), 0, P.comm(),
        &status);

    for (l = 0, k = 0; k < p; k++)
      for (i = 0; i < n; i++)
        coef(i, 0, k) = buf_out[l++];
  }

  if (P.neighbor(1, 1) != MPI_PROC_NULL) {
    std::vector<double> buf_in(n * p);
    std::vector<double> buf_out(n * p);
    for (l = 0, k = 0; k < p; k++)
      for (i = 0; i < n; i++)
        buf_in[l++] = coef(i, m - 2, k);

    MPI_Sendrecv(buf_in.data(), n * p, MPI_DOUBLE, P.neighbor(1, 1), 0,
        buf_out.data(), n * p, MPI_DOUBLE, P.neighbor(1, 1), 0, P.comm(),
        &status);

    for (l = 0, k = 0; k < p; k++)
      for (i = 0; i < n; i++)
        coef(i, m - 1, k) = buf_out[l++];
  }

  if (P.neighbor(2, 0) != MPI_PROC_NULL) {
    std::vector<double> buf_in(n * m);
    std::vector<double> buf_out(n * m);
    for (l = 0, j = 0; j < m; j++)
      for (i = 0; i < n; i++)
        buf_in[l++] = coef(i, j, 1);

    MPI_Sendrecv(buf_in.data(), n * m, MPI_DOUBLE, P.neighbor(2, 0), 0,
        buf_out.data(), n * m, MPI_DOUBLE, P.neighbor(2, 0), 0, P.comm(),
        &status);

    for (l = 0, j = 0; j < m; j++)
      for (i = 0; i < n; i++)
        coef(i, j, 0) = buf_out[l++];
  }

  if (P.neighbor(2, 1) != MPI_PROC_NULL) {
    std::vector<double> buf_in(n * m);
    std::vector<double> buf_out(n * m);
    for (l = 0, j = 0; j < m; j++)
      for (i = 0; i < n; i++)
        buf_in[l++] = coef(i, j, p - 2);

    MPI_Sendrecv(buf_in.data(), n * m, MPI_DOUBLE, P.neighbor(2, 1), 0,
        buf_out.data(), n * m, MPI_DOUBLE, P.neighbor(2, 1), 0, P.comm(),
        &status);

    for (l = 0, j = 0; j < m; j++)
      for (i = 0; i < n; i++)
        coef(i, j, p - 1) = buf_out[l++];
  }
}


