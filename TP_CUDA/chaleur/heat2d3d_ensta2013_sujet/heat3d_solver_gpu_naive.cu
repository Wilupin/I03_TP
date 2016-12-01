/**
 * \file heat3d_solver_gpu_naive.cu
 * \brief Solve 3D heat equation (finite difference method). GPU version (naive).
 *
 * We solve the 3D Heat equation \f$\partial_t \phi = \alpha \left[
 * \partial^2_x \phi + \partial^2_y \phi + \partial^2_z \ phi \right] \f$, \f$ 0 \leq x
 * \leq L_x \f$, \f$ 0 \leq y \leq L_y \f$, \f$ 0 \leq t\f$.\\
 *
 * Method : Finite Difference, FTCS scheme
 *
 * GPU Features: use only global memory
 *
 * boundary condition : Dirichlet
 *
 * GPU version : naive
 *
 * \date 27-dec-2009.
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <sys/time.h> // gettimeofday
#include <assert.h>

// includes, project
#include <helper_functions.h>
#include "CudaTimer.h"
#include "Timer.h"

// parameters and real_t typedef
#include "param.h"

// for output results
#include "output.h"

// GPU solver
#include "heat3d_kernel_gpu_naive.cu"

// CPU solver
#include "heat_kernel_cpu.h"

// initial conditions
#include "misc.h"

// cuda helper
#include "cuda_helper.cu"

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest( int argc, char** argv);


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main(int argc, char** argv) 
{
  runTest(argc, argv);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void
runTest(int argc, char** argv) 
{
  int devID;
  cudaDeviceProp deviceProps;
  
  devID = findCudaDevice(argc, (const char **)argv);
  
  // get number of SMs on this GPU
  checkCudaErrors(cudaGetDeviceProperties(&deviceProps, devID));
  printf("CUDA device [%s] has %d Multi-Processors\n", deviceProps.name, deviceProps.multiProcessorCount);
  
  /*
   * read and print parameters
   */
  // default parameter file
  std::string paramFile("heatEqSolver.par");

  // if argv[1] exists use it as a parameter file
  if (argc>1) {
    printf("trying to read parameters from file %s ...\n",argv[1]);
    paramFile = std::string(argv[1]);
  }

  // read parameter file
  readParamFile(paramFile);

  if (NZ<=1) {
    printf("NZ should be larger than 1 in the 3D version\n");
    cudaThreadExit();
  }

  // print parameters on screen
  printParameters("HEAT 3D - GPU (NAIVE)");

  CudaTimer gpuTimer;

  unsigned int mem_size = sizeof(real_t)*NX*NY*NZ;

  // allocate host memory
  real_t* data1 = (real_t*) malloc( mem_size);
  real_t* data2 = (real_t*) malloc( mem_size);
  
  ///////////////////////////////////////////////////
  // compute GPU solution to 3D heat equation
  ///////////////////////////////////////////////////
  
  // inital condition
  initCondition3D (data1);
  
  // allocate device memory
  real_t* d_data1;
  real_t* d_data2;
  
  // naive kernel memory allocation (using cudaMalloc)
  checkCudaErrors( cudaMalloc( (void**) &d_data1, mem_size));
  checkCudaErrors( cudaMalloc( (void**) &d_data2, mem_size));

  // copy host memory to device
  checkCudaErrors( cudaMemcpy( d_data1, data1, mem_size,
			     cudaMemcpyHostToDevice) );
  checkCudaErrors( cudaMemcpy( d_data2, data1, mem_size,
			     cudaMemcpyHostToDevice) );


   
  // setup execution parameters for cuda kernel
  // grid dimension for naive kernel
  unsigned int threadsPerBlockX=32;
  unsigned int threadsPerBlockY=10;
  dim3  threads(threadsPerBlockX,threadsPerBlockY);
  dim3  grid( (NX+threads.x-1)/threads.x, (NY+threads.y-1)/threads.y );
    
  printf("grid  size : %u %u %u\n",grid.x,grid.y,grid.z);
  printf("block size : %u %u %u\n",threads.x,threads.y,threads.z);

  // start timer
  gpuTimer.start();

  // time loop executing naive kernel
  unsigned int iTime=0;
  int iOutput=-1;
  for (real_t t=0.0f; t<TMAX; t+=(2*DT), iTime+=2) {

    if (useOrder2) { // use the 2nd order accurate scheme
      
      heat3d_ftcs_naive_order2_kernel<<< grid, threads >>>( d_data1, d_data2, 
							    NX, NY, NZ,
							    o2.R, o2.R3);
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
      
      heat3d_ftcs_naive_order2_kernel<<< grid, threads >>>( d_data2, d_data1, 
							    NX, NY, NZ,
							    o2.R, o2.R3);   
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");

    } else if (useOrder2b) { 
      
      heat3d_ftcs_naive_order2b_kernel<<< grid, threads >>>( d_data1, d_data2, 
							     NX, NY, NZ,
							     o2.R, o2.R3b);
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
      
      heat3d_ftcs_naive_order2b_kernel<<< grid, threads >>>( d_data2, d_data1, 
							    NX, NY, NZ,
							    o2.R, o2.R3b);   
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");

    } else {

      heat3d_ftcs_naive_order4_kernel<<< grid, threads >>>( d_data1, d_data2, 
							    NX, NY, NZ,
							    o4.S, o4.S3);
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
      
      heat3d_ftcs_naive_order4_kernel<<< grid, threads >>>( d_data2, d_data1, 
							    NX, NY, NZ,
							    o4.S, o4.S3);   
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");

    }
    
    /* save output (just for cross-checking, do not save when
       measuring computing time */
    if (ENABLE_GPU_SAVE) {
      
      if (iTime%T_OUTPUT == 0) {
	iOutput++;
	checkCudaErrors( cudaMemcpy( data1, d_data1, NX*NY*NZ*sizeof( real_t),
				   cudaMemcpyDeviceToHost) );      
      }
      // VTK output
      if (SAVE_VTK and iTime%T_OUTPUT == 0)
	save_vtk(data1, "heat3d_gpu_naive_",iOutput);

      // HDF5 output
      if (SAVE_HDF5 and iTime%T_OUTPUT == 0)
	save_hdf5(data1, "heat3d_gpu_naive_",iOutput);

    }

  } // end for loop
  
  // stop timer
  gpuTimer.stop();

  real_t gpu_time = gpuTimer.elapsed();
   printf( "GPU Processing time: %f (s)\n", gpu_time);
  
  // copy result from device to host
  real_t *resGPU = (real_t*) malloc( mem_size);
  checkCudaErrors( cudaMemcpy( resGPU, d_data1, mem_size,
			     cudaMemcpyDeviceToHost) );
    
  if (SAVE_HDF5)
    write_xdmf_wrapper("heat3d_gpu_naive",N_ITER,T_OUTPUT);

  ////////////////////////////////////////////////////////
  // compute reference (CPU) solution to 3D heat equation
  // for performance comparison
  ////////////////////////////////////////////////////////
  initCondition3D (data1);
  initCondition3D (data2);
  
  Timer cpuTimer;
  cpuTimer.start();
  
  // time loop
  iTime=0;

  for (real_t t=0.0f; t<TMAX; t+=(2*DT), iTime+=2) {
    
    // compute next 2 time steps
    if (useOrder2) {
      heat3d_ftcs_cpu_order2( data1, data2);
      heat3d_ftcs_cpu_order2( data2, data1);
    } /*else if (useOrder2b) {
      heat3d_ftcs_cpu_order2b( data1, data2);
      heat3d_ftcs_cpu_order2b( data2, data1);
      } */ else {
      heat3d_ftcs_cpu_order4( data1, data2);
      heat3d_ftcs_cpu_order4( data2, data1);
    }
  }
  
  // stop timer
  cpuTimer.stop();
  real_t cpu_time = cpuTimer.elapsed();
  
  printf( "CPU Processing time: %g (s)\n", cpu_time);
  printf( "Speedup GPU/CPU : %f\n",cpu_time/gpu_time);
   
  printf("...comparing the results\n");
  double sum = 0, delta = 0;
  for(unsigned i = 0; i < NX*NY*NZ; i++){
    delta += (resGPU[i] - data1[i]) * (resGPU[i] - data1[i]);
    sum   += data1[i] * data1[i];
  }
  double L2norm = sqrt(delta / sum);
  printf("iteration %d relative L2 norm: %E\n", iTime, L2norm);
 
  // cleanup memory
  free(data1);
  free(data2);
  free(resGPU);
  
  checkCudaErrors(cudaFree(d_data1));
  checkCudaErrors(cudaFree(d_data2));
  
  cudaDeviceSynchronize();
  cudaDeviceReset();

  exit(0);
}
