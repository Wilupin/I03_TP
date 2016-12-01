/*
 * config #1
 * nvcc --generate-code arch=compute_13,code=sm_13 --ptxas-options -v limitingFactor.cu -o limitingFactor_v1
 *
 * config #2
 * nvcc --generate-code arch=compute_13,code=sm_13 --ptxas-options -v --use_fast_math limitingFactor.cu -o limitingFactor_v2
 *
 * Note that fastmath option triggers the use of the special unit hardware; less IEEE754 compliant;
 * denormal numbers are flushed to zero; division is approximated; square root is approximated
 *
 */
#include <stdlib.h>
#include <stdio.h>

typedef float real_t;
//typedef double real_t;

__global__ void base (real_t* a, real_t* b) 
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  a[i] = sin ( b[i] ) + cos ( b[i] );
} // end base

__global__  void memory (real_t* a, real_t* b)
{

  int i = blockIdx.x * blockDim.x + threadIdx.x;
  a[i] = b[i];
} // end memory

  
__global__ void math (real_t* a, real_t b, int flag )
{

  real_t v;
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  v = sin ( b ) + cos ( b );
  if ( v * flag == 1) a [ i ] = v;
} // end math

__global__ void initArray(real_t* data, real_t value)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  data[i] = value;
}

int main(int argc, char* argv[])
{

  const int n = 8*1024*1024 , blockSize = 256;
  real_t* a;
  real_t* a_d;
  real_t* b_d;

  // memory allocation
  a = (real_t*) malloc(n*sizeof(real_t));
  cudaMalloc((void**)&a_d, n*sizeof(real_t));
  cudaMalloc((void**)&b_d, n*sizeof(real_t));

  // memory initialization
  initArray<<< n / blockSize , blockSize >>> ( b_d , 1.0 );
  base<<< n / blockSize , blockSize >>> ( a_d , b_d );
  memory   <<< n / blockSize , blockSize >>> ( a_d , b_d );
  math     <<< n / blockSize , blockSize >>> ( a_d , 1.0 , 0);

  // copy back
  cudaMemcpy(a,a_d,n*sizeof(real_t),cudaMemcpyDeviceToHost);
 
  printf("%f\n", a[0]);

  return 0;
}
