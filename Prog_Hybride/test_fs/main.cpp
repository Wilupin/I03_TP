#include <iostream>
#include "timer.hpp"
#include <math.h>
#include <omp.h>

int main(int argc, char **argv)
{

  const size_t n = 4;
  const size_t m = 40000000L;
  const size_t offset = argc > 1 ? atoi(argv[1]) : 1;

  size_t i, j;
  double **A;
  A = new double*[n];

#pragma omp parallel for private(i)
  for (i=0; i<n; i++)
    A[i] = new double[m];

  double s[n*offset];
  for (i=0; i<n; i++) s[i*offset] = 0.0;


  Timer t1, t2;

  t1.start();

#pragma omp parallel for private(i,j)
  for (i=0; i<n; i++)
    for (j=0; j<m; j++)
      A[i][j] = exp(-(i+j-1.0)) * sin(j*2.0);
  t1.stop();

  std::cout << "initialisation : " << t1.elapsed() << " s" << std::endl;

  t2.start();
#pragma omp parallel for private(i,j)
  for (i=0; i<n; i++) {
    double d;
    s[i*offset] = 0;
    for (j=0; j<m; j++) {
       d = A[i][j] * cos(3*j);
       s[i * offset] += d;
    }
  }
  t2.stop();

  std::cout << "calcul         : " << t2.elapsed() << " s" 
	          << std::endl << std::endl;
  return 0;
}

