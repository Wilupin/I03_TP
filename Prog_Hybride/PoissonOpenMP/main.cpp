#include "Poisson_Parameters.hpp"
#include "Poisson.hpp"
#include "plot.hpp"

double f(double x, double y, double z)
{
  x -= 0.5;
  y -= 0.5;
  z -= 0.5;
  if (x*x+y*y+z*z < 0.1)
    return 1.0;
  else
    return 0.0;
}

int main(int argc, char *argv[])
{
  Poisson_Parameters Prm(&argc, &argv);
  if (Prm.help()) return 0;
  std::cout << Prm << std::endl;
  
  int freq = Prm.freq();
  int itMax = Prm.it();

  int nsteps = freq > 0 ? itMax/freq : 1;
  int ksteps = freq > 0 ? freq : itMax;
  
  Poisson C;
  Plot P;
  Values u_0;

  C.timer(0).start();
  u_0.init(Prm.n(), Prm.dx(), f);
  C.init(Prm);
  C.setInput(u_0);
  C.timer(0).stop();
  
  if (freq > 0) {
    int n = C.getDomainSize(0);
    int m = C.getDomainSize(1);
    int p = C.getDomainSize(2);
    P.init(n, m, p, Prm);
    P.process(C.getOutput());
  }

  int i;
  for (i=0; i<nsteps; i++) {
    C.solve(ksteps);
    if (freq > 0) {
      P.process(C.getOutput());
    }
  }
  
  C.timers();
  C.terminate();
  if (freq > 0)
    P.terminate();
  
  return 0;
}
