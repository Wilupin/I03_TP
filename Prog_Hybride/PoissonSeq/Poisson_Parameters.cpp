#include "osutils.hpp"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <limits>
#include "Poisson_Parameters.hpp"

Poisson_Parameters::Poisson_Parameters(int *argc, char *** argv)
  : m_prm(*argc, *argv)
{
  int i;
  
  m_command = (*argv)[0];  
  m_help = m_prm.search(2, "-h", "--help");
  m_debug = m_prm("debug", 0) != 0;

  m_n[0] = m_prm("n", 100);
  m_n[1] = m_prm("m", 100);
  m_n[2] = m_prm("p", 100);
  m_it = m_prm("it", 10);
  m_freq = m_prm("out", 0);

  if (m_help) return;

  m_dt = -1.0;
   
  if (m_dt <= 0.0)
    m_dt = 0.5/(m_n[0]*m_n[0]
		+ m_n[1]*m_n[1]
		+ m_n[2]*m_n[2]);

  for (i=0; i<3; i++) {
    m_dx[i] = m_n[i]>1 ? 1.0/(m_n[i]-1) : 0.0;
    m_di[i] = 1;
    m_imin[i] = 1;
    m_imax[i] = m_n[i]-1;
    if (m_n[i] < 2) {
      m_imin[i]=0; m_imax[i] = 1; m_di[i] = 0;
    }
  }
  
  if (m_freq > 0) {
    char buffer[256];
    stime(buffer, 256);
    
    std::ostringstream pth;
    pth << "results"
	<< "_n_" << m_n[0] << "x" << m_n[1] << "x" << m_n[2]
	<< "_" << buffer << "/";
    m_path = pth.str();
    
    mkdir_p(m_path.c_str());
    
    std::ostringstream s;
    s << m_path << "/out.txt";
    m_out = new std::ofstream(s.str().c_str());
  }
  else
    m_out = NULL;
}



bool Poisson_Parameters::help()
{
  if (m_help) {
    std::cerr << "Usage : ./PoissonSeq <list of options>\n\n";
    std::cerr << "Options:\n\n"
              << "-h|--help : display this message\n"
              << "n=<int>       : number of internal points in the X direction (default: " << m_n[0] << ")\n"
              << "m=<int>       : number of internal points in the Y direction (default: " << m_n[1] << ")\n"
              << "p=<int>       : number of internal points in the Z direction (default: " << m_n[2] << ")\n"
              << "dt=<real>     : time step size (default : value to assure stable computations)\n"
              << "it=<int>      : number of time steps (default : " << m_it << ")\n"
              << "out=<int>     : number of time steps between saving the solution on files\n"
              << "                (default : 0 or no output)\n\n";
  }
  return m_help;
}

Poisson_Parameters::~Poisson_Parameters()
{
  if (!m_help) {
    delete m_out;
  }
}

std::ostream & operator<<(std::ostream &f, const Poisson_Parameters & p)
{
  f << "Domain :   "
    << "[" << 0 << "," << p.n(0) - 1  << "] x "
    << "[" << 0 << "," << p.n(1) - 1  << "] x "
    << "[" << 0 << "," << p.n(2) - 1  << "]"
    << "\n\n";

  f << "Dt :         " << p.dt() << "\n";
  f << "Iterations : " << p.it() << "\n";

  return f;
}
