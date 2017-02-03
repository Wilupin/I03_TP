#include "Poisson.hpp"
#include "Poisson_Parameters.hpp"

#include <mpi.h>
#include <sstream>

#include "version.hpp"

Poisson::Poisson() : codeName(version), debug(false)
{
  kStep = 1;
  m_t = 0.0;
  int i;
  for (i=0; i<3; i++) {
    m_n[i] = 0;
    m_dx[i] = 0.0;
    m_di[i] = 0;
  }

  T.resize(3);
  T[0].name("init");
  T[1].name("solve");
  T[2].name("comm");

  m_duv = 0.0;
  m_P = NULL;
}

Poisson::~Poisson()
{
  if (m_P->rank() == 0)
    std::cerr << "delete " << codeName << std::endl;
}

double Poisson::present()
{
  return m_t;
}

void Poisson::init(Poisson_Parameters & P)
{
  int i;
  if (P.rank() == 0)
    std::cerr << "init " << codeName << std::endl;
  debug = P.debug();

  kStep = 1;
  m_t = 0;

  for (i=0; i<3; i++) {
    m_n[i] = P.n(i);
    m_dx[i] = P.dx(i);
    m_di[i] = (m_n[i] < 2) ? 0 : 1;
  }

  m_u.resize(m_n);
  m_v.resize(m_n);

  m_P = &P;
}

size_t Poisson::getDomainSize(int dim) const
{
  size_t d;
  switch (dim) {
    case 0:
      d = m_n[0];
      break;
    case 1:
      d = m_n[1];
      break;
    case 2:
      d = m_n[2];
      break;
    default:
      d = 1;
  }
  return d;
}


bool Poisson::solve(unsigned int nSteps)
{
  T[1].start();

  if (debug) {
    std::ostringstream s;
    s << "output_solve_" << kStep;
    fDebug.open(s.str().c_str());
  }

  double du_max, du_max_loc, du;
  double dx2 = m_dx[0]*m_dx[0] + m_dx[1]*m_dx[1] + m_dx[2]*m_dx[2];
  double dt = 0.5*(dx2 + 1e-12);
  double lambda = 0.25*dt/(dx2 + 1e-12);

  int i, j, k;
  size_t iStep;
  int   di = m_di[0],     dj = m_di[1],     dk = m_di[2];
  int imin = 1, imax = m_n[0] - 1;
  int jmin = 1, jmax = m_n[1] - 1;
  int kmin = 1, kmax = m_n[2] - 1;

  for (iStep=0; iStep < nSteps; iStep++) {

    du_max = 0.0;
    du_max_loc = 0.0;

    for (i = imin; i < imax; i++)
      for (j = jmin; j < jmax; j++)
        for (k = kmin; k < kmax; k++) {
          du = 6 * m_u(i, j, k)
              - m_u(i + di, j, k) - m_u(i - di, j, k)
              - m_u(i, j + dj, k) - m_u(i, j - dj, k)
              - m_u(i, j, k + dk) - m_u(i, j, k - dk);
          du *= lambda;
          m_v(i, j, k) = m_u(i, j, k) - du;
          du_max_loc += du > 0 ? du : -du;
        }

    T[1].stop();
    T[2].start();


    // Reduction de la variable du_max sur le noeud 0
    MPI_Reduce(&du_max_loc, &du_max, 1, MPI_DOUBLE, MPI_SUM, 0, m_P->comm());


    // Syncronisation de la solution
    m_v.synchronize(*m_P);

    T[2].stop();
    T[1].start();

    m_u.swap(m_v);
    m_t += dt;
    kStep++;

    //std::cout << " Mon rang est : " << m_P->rank() << std::endl;

    if (m_P->rank() == 0)
      //std::cout << "Ici" << std::endl;
      std::cerr << "solve    " << codeName
        << " iteration " << kStep
              << " variation " << du_max
              << "    \r";

  }

  m_duv = du_max;

  if (debug) {
    fDebug.close();
  }
  T[1].stop();
  return true;
}

double Poisson::variation()
{
  return m_duv;
}

void Poisson::terminate() {
  if (m_P->rank() == 0)
    std::cerr << "terminate " << codeName << std::endl;
}

const Values & Poisson::getOutput()
{
  return m_u;
}

void Poisson::setInput(const Values & u)
{
  m_u = u;
  m_v = u;
}

void Poisson::save(const char * /*fName*/)
{
}

std::vector<double> Poisson::timers()
{
  displayTimers(std::cerr, T);
  size_t i, n = T.size();
  std::vector<double> ts(n);
  for (i=0; i<n; i++)
    ts[i] = T[i].elapsed();
  return ts;
}


