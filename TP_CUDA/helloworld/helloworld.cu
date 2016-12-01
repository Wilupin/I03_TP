/* 
 * nvcc -arch=sm_13 --ptxas-options -v -o helloworld helloworld.cu
 */

#include <stdio.h>
#include <cuda_runtime.h>

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
 * a simple CUDA kernel
 *
 * \param[in]  a input integer
 * \param[in]  b input integer
 * \param[out] c pointer-to-integer for result
 */
__global__ void add( int *a, int *b, int *c, int n ) {

  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if (i<n)
    c[i] = a[i] + b[i];
}

/*
 * main
 */
int main( void ) {
  // array size
  int N = 30000000;

  // host variables
  int *a, *b, *c;

  // device variables
  int *dev_a, *dev_b, *dev_c;
  
  // CPU memory allocation / initialization
  a = (int *) malloc(N*sizeof(int));
  b = (int *) malloc(N*sizeof(int));
  c = (int *) malloc(N*sizeof(int));
  for (int i=0; i<N; i++) {
    a[i]=i;
    b[i]=N-i;
  }

  // GPU device memory allocation / initialization
  cudaMalloc( (void**)&dev_a, N*sizeof(int) );
  cudaMalloc( (void**)&dev_b, N*sizeof(int) );
  cudaMalloc( (void**)&dev_c, N*sizeof(int) );
  cudaMemcpy( dev_a, a, N*sizeof(int),
	      cudaMemcpyHostToDevice );
  cudaMemcpy( dev_b, b, N*sizeof(int),
	      cudaMemcpyHostToDevice );

  // perform computation on GPU
  int nbThreads = 32;
  dim3 blockSize(nbThreads,1,1);
  dim3 gridSize(N/nbThreads+1,1,1);
  printf("Using %d blocks of %d threads\n",gridSize.x,blockSize.x);
  add<<<gridSize,blockSize>>>( dev_a, dev_b, dev_c, N );
  HANDLE_ERROR( "add" );


  // get back computation result into host CPU memory
  cudaMemcpy( c, dev_c, N*sizeof(int),
	      cudaMemcpyDeviceToHost);

  // output result on screen
  int passed=1;
  for (int i=0; i<N; i++) {
    if (c[i] != N) {
      passed = 0;
      printf("wrong value : %d %d\n",i,c[i]);
    }
  }
  if (passed) {
    printf("test succeeded !\n");
  } else {
    printf("test failed !\n");
  }

  // de-allocate CPU host memory
  free(c);
  free(b);
  free(a);

  // de-allocate GPU device memory
  cudaFree( dev_c );
  cudaFree( dev_b );
  cudaFree( dev_a );

  return 0;
}
