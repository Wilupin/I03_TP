/* 
 * nvcc -O3 -arch=sm_13 --ptxas-options -v -o helloworld_copy helloworld_copy.cu
 */

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

// load our GPU specififc timer
#include "../CudaTimer.h"


/*
 * handle CUDA-related error messages
 */
static void HandleError( const char* kernelName,
                         const char *file,
                         int line ) {
  cudaError_t err = cudaPeekAtLastError();
  
  if (err != cudaSuccess) {
    printf("Kernel %s FAILED in %s at line %d with error message:\n%s\n", 
	   kernelName,
	   file, 
	   line, 
	   cudaGetErrorString(cudaGetLastError()));
    exit( EXIT_FAILURE );
  }
}
#define HANDLE_ERROR( kernelName ) (HandleError( kernelName, __FILE__, __LINE__ ))


/**
 * a simple CUDA kernel to copy a into b.
 *
 * Using a 1-to-1 mapping between threads index and memory addresses.
 *
 * \param[in]  a input integer pointer
 * \param[in]  b output integer pointer
 */
__global__ void copy( int *a, int *b, int n ) {

  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i<n)
    b[i] = a[i];
}

/**
 * a simple CUDA kernel to copy a into b.
 *
 * Using a fixed size grid of block of threads.
 *
 * \param[in]  a input integer pointer
 * \param[in]  b output integer pointer
 */
__global__ void copy2( int *a, int *b, int n ) {

  int i = threadIdx.x + blockIdx.x*blockDim.x;

  while (i<n) {
    b[i] = a[i];

    i += blockDim.x*gridDim.x;
  }

} // end copy2

/**
 * a simple CUDA kernel to copy a into b.
 *
 * Using a fixed size grid of block of threads + loop unrolling technique.
 *
 * \param[in]  a input integer pointer
 * \param[in]  b output integer pointer
 */
__global__ void copy3( int *a, int *b, int n ) {

  const int nUnroll = 4;

  for (int i = threadIdx.x + blockIdx.x*blockDim.x * nUnroll;
       i < n;
       i += nUnroll*blockDim.x*gridDim.x) {

    for(int j=0; j<nUnroll; j++) {
      int index = i+j*blockDim.x;
      if (index<n)
	b[index] = a[index];
    }

  }

} // end copy3

/*
 * use a random number generator to pick a GPU (0, 1, 2 or 3)
 */
int get_GPU_device_number() {

  // generate a double value in [0, 3[
  double randomNum = 4.0 * rand() / (RAND_MAX+1.0);
  
  // return integer
  return (int) randomNum;

} // get GPU device number


/*
 * main
 */
int main( int argc, char* argv[] ) {

  // array size
  int N;

  // host variables
  int *a, *b;

  // device variables
  int *dev_a, *dev_b;

  if (argc>1) {
    N = atoi(argv[1]);
  } else {
    // assign a default value
    N = 1000000;
  }

  // initialize random number generator (used in get_GPU_device_number
  srand(time(NULL));
  
  // choose a random GPU out of GPU available
  int deviceNumber = get_GPU_device_number();
  cudaSetDevice( deviceNumber );
  printf("We are using GPU number : %d\n", deviceNumber);


  // get device properties of the chosen device
  cudaDeviceProp prop ;
  cudaGetDeviceProperties (&prop , deviceNumber ); 

  printf("nombre de multi-processor %d sur le GPU %d\n",prop.multiProcessorCount,deviceNumber);

  

  // CPU memory allocation / initialization
  a = (int *) malloc(N*sizeof(int));
  b = (int *) malloc(N*sizeof(int));
  for (int i=0; i<N; i++) {
    a[i]=i;
  }

  // GPU device memory allocation / initialization
  cudaMalloc( (void**)&dev_a, N*sizeof(int) );
  cudaMalloc( (void**)&dev_b, N*sizeof(int) );
  cudaMemcpy( dev_a, a, N*sizeof(int),
	      cudaMemcpyHostToDevice );
  cudaMemset( (void*)dev_b, 0, N*sizeof(int) );

  // print maximum achievable peak bandwidth
  // memoryClockRate is given in kHz, so divide by 10^6 to have giga transfert per sec
  double peakBandwidth = 2.0*prop.memoryClockRate/1000000*prop.memoryBusWidth/8;
  printf("Peak bandwidth : %6.2f GB/s\n", peakBandwidth);


  // perform computation on GPU

  /*
   * test copy on GPU
   */
  if (1) {
    CudaTimer copyTimer;

    int nbThreads = 256;
    dim3 blockSize(nbThreads,1,1);
    dim3 gridSize(N/nbThreads+1,1,1);
    copyTimer.start();
    copy<<<gridSize,blockSize>>>( dev_a, dev_b, N );
    copyTimer.stop();
    HANDLE_ERROR( "copy" );

    printf("copy  kernel time   : %f\n",copyTimer.elapsed_in_second());
    double bandWidth = (double) 2.0*N*sizeof(int)/copyTimer.elapsed_in_second()/(1024*1024*1024);
    printf("Effective bandwidth : %6.2f GBytes/s (%5.2f %% of peak bandwidth)\n", bandWidth ,100*bandWidth/peakBandwidth);

   
    // check results are OK
    // get back computation result into host CPU memory
    cudaMemcpy( b, dev_b, N*sizeof(int),
		cudaMemcpyDeviceToHost);

    // output result on screen
    int passed=1;
    for (int i=0; i<N; i++) {
      if (b[i] != i) {
	passed = 0;
	printf("wrong value : %d %d\n",i,b[i]);
      }
    }
    if (passed) {
      printf("test succeeded !\n");
    } else {
      printf("test failed !\n");
    }
  } // end test copy

  /*
   * test copy2 on GPU
   */
  if (1) {
    CudaTimer copy2Timer;

    int nbThreads = 256;
    dim3 blockSize(nbThreads,1,1);
    dim3 gridSize(100*prop.multiProcessorCount,1,1);
    copy2Timer.start();
    copy2<<<gridSize,blockSize>>>( dev_a, dev_b, N );
    copy2Timer.stop();
    HANDLE_ERROR( "copy2" );
    
    printf("copy2 kernel time   : %f\n",copy2Timer.elapsed_in_second());
    double bandWidth = (double) 2.0*N*sizeof(int)/copy2Timer.elapsed_in_second()/(1024*1024*1024);
    printf("Effective bandwidth : %6.2f GBytes/s (%5.2f %% of peak bandwidth)\n", bandWidth ,100*bandWidth/peakBandwidth);

    // check results are OK
    // get back computation result into host CPU memory
    cudaMemcpy( b, dev_b, N*sizeof(int),
		cudaMemcpyDeviceToHost);

    // output result on screen
    int passed=1;
    for (int i=0; i<N; i++) {
      if (b[i] != i) {
	passed = 0;
	printf("wrong value : %d %d\n",i,b[i]);
      }
    }
    if (passed) {
      printf("test succeeded !\n");
    } else {
      printf("test failed !\n");
    }
  } // end test copy2

  /*
   * use cudaMemcpy for reference
   */
  if (1) {
    CudaTimer copy3Timer;

    copy3Timer.start();
    cudaMemcpy(dev_b, dev_a, N*sizeof(int), cudaMemcpyDeviceToDevice);
    copy3Timer.stop();

    printf("copy3 kernel time   : %f\n",copy3Timer.elapsed_in_second());
    double bandWidth = (double) 2.0*N*sizeof(int)/copy3Timer.elapsed_in_second()/(1024*1024*1024);
    printf("Effective bandwidth : %6.2f GBytes/s (%5.2f %% of peak bandwidth)\n", bandWidth ,100*bandWidth/peakBandwidth);

  }


  // de-allocate CPU host memory
  free(b);
  free(a);

  // de-allocate GPU device memory
  cudaFree( dev_b );
  cudaFree( dev_a );

  return 0;
}
