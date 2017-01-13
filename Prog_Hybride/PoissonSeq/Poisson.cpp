#include "Poisson.hpp"
#include "Poisson_Parameters.hpp"

#include <fstream>
#include <sstream>

#include "version.hpp"

Poisson::Poisson() : codeName(version), debug(false)
{
  int i;
  for (i=0; i<3; i++) {
    m_n[i] = 0;
    m_dx[i] = 0.0;
    m_di[i] = 0;
    m_imin[i] = 1;
    m_imax[i] = 0;
  }
  
  T.resize(2);
  T[0].name("init");
  T[1].name("solve");
}

Poisson::~Poisson()
{
  std::cerr << "delete " << codeName << std::endl;
}

double Poisson::present()
{
  return m_t;
}

void Poisson::init(Poisson_Parameters & P)
{
  int i;
  std::cerr << "init " << codeName << std::endl;
  debug = P.debug();

  kStep = 1;  
  m_t = 0;
  
  for (i=0; i<3; i++) {
    m_n[i] = P.n(i);
    m_dx[i] = m_n[i]>1 ? 1.0/(m_n[i]-1) : 0.0;
    m_di[i] = 1;
    m_imin[i] = 1;
    m_imax[i] = m_n[i]-1;
    if (m_n[i] < 2) {
      m_imin[i]=0; m_imax[i] = 1; m_di[i] = 0;
    }
  }

  m_u.resize(m_n);
  m_v.resize(m_n);
}

size_t Poisson::getDomainSize(int dim) const
{
  size_t d;
  switch (dim) {
  case 0:
    d = m_n[0];
    break;
  case 1:
    d = m_n[1] ;
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

  double du_max = 0.0, du;
  double dx2 = m_dx[0]*m_dx[0] + m_dx[1]*m_dx[1] + m_dx[2]*m_dx[2];
  double dt = 0.5*(dx2 + 1e-12);
  double lambda = 0.25*dt/(dx2 + 1e-12);

  int i, j, k;
  size_t iStep;
  int   di = m_di[0],     dj = m_di[1],     dk = m_di[2];
  int imin = m_imin[0], jmin = m_imin[1], kmin = m_imin[2];
  int imax = m_imax[0], jmax = m_imax[1], kmax = m_imax[2];

  for (iStep=0; iStep < nSteps; iStep++) {

    du_max = 0.0;
    
    for (i=imin; i<imax; i++)
      for (j=jmin; j<jmax; j++)
	for (k=kmin; k<kmax; k++) {
	  du = 6*m_u(i,j,k)
	    - m_u(i+di,j,k) - m_u(i-di,j,k)
	    - m_u(i,j+dj,k) - m_u(i,j-dj,k)
	    - m_u(i,j,k+dk) - m_u(i,j,k-dk);
	  du *= lambda;
	  m_v(i,j,k) = m_u(i,j,k) - du;
	  du_max += du > 0 ? du : -du;
	}

    m_u.swap(m_v);
    m_t += dt;
    kStep++;
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


