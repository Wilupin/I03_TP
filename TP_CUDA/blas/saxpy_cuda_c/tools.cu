// collection of useful routines

#include "tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <cuda_runtime.h>
#include <cublas.h>


/////////////////////////////////////
// initialises CUDA and directs all
// computations to the given
// CUDA device
/////////////////////////////////////
void initCuda(const int selectedDevice)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0)
  {
    fprintf(stderr, "Sorry, no CUDA device fount");
    exit(1);
  }
  if (selectedDevice >= deviceCount)
  {
    fprintf(stderr, "Choose device ID between 0 and %d\n", deviceCount-1);
    exit(2);
  }
  cudaSetDevice(selectedDevice);
  checkErrors("initCuda");

  cublasInit();
}



/////////////////////////////////////
// error checking 
/////////////////////////////////////
void checkErrors(char *label)
{
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }
}




/////////////////////////////////////
// CPU based timings 
/////////////////////////////////////
struct timeval timer_start;


// clears timer and restarts clock ticking
void startTimer(void)
{
  gettimeofday(&timer_start, NULL);
}

// stops timer and returns elapsed time in secs
double stopTimer(void)
{
  struct timeval end;
  double elapsed;
  gettimeofday(&end, NULL);
  elapsed = (end.tv_sec-timer_start.tv_sec) + (end.tv_usec - timer_start.tv_usec)/1000000.0;
  return elapsed;
}
