/**
 * \file heat2d_solver_gpu_naive.cu
 * \brief Solve 2D heat equation (finite difference method). GPU version (naive).
 *
 * We solve the 2D Heat equation \f$\partial_t \phi = \alpha \left[
 * \partial^2_x \phi + \partial^2_y \phi \right] \f$, \f$ 0 \leq x
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
 * \date 17-dec-2009.
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <sys/time.h> // gettimeofday

// includes, project
#include <helper_functions.h>
#include "CudaTimer.h"
#include "Timer.h"

// parameters + real_t typedef
#include "param.h"

// for output results
#include "output.h"

// GPU solver
#include "heat2d_kernel_gpu_naive.cu"

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

  exit(0);
}

////////////////////////////////////////////////////////////////////////////////
//! Run solver on GPU
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

  // print parameters on screen
  printParameters("HEAT 2D - GPU (NAIVE)");


  CudaTimer gpuTimer;

  unsigned int mem_size = sizeof(real_t)*NX*NY;

  // allocate host memory
  real_t* data1 = (real_t*) malloc( mem_size);
  real_t* data2 = (real_t*) malloc( mem_size);
  
  ///////////////////////////////////////////////////
  // compute GPU solution to 2D heat equation
  ///////////////////////////////////////////////////
  
  // inital condition
  initCondition2D (data1);
  
  // allocate device memory
  real_t* d_data1;
  real_t* d_data2;
  
  // naive kernel memory allocation (using cudaMalloc)
  /* 
   * TODO
   */

  // copy host memory to device
  /*
   * TODO
   */
    
   
  // setup execution parameters for cuda kernel
  // grid dimension for naive kernel
  unsigned int threadsPerBlockX=16;
  unsigned int threadsPerBlockY=16;
  dim3  threads(/* TODO */);
  dim3  grid(/* TODO */);
    
  printf("grid  size : %u %u\n",grid.x,grid.y);
  printf("block size : %u %u\n",threads.x,threads.y);

  // start timer
  gpuTimer.start();

  // time loop executing naive kernel
  int iTime=0;
  int iOutput=-1;
  for (real_t t=0.0f; t<TMAX; t+=(2*DT), iTime+=2) {
    
    if (useOrder2) { // use the 2nd order accurate scheme
      
      heat2d_ftcs_naive_order2_kernel<<< grid, threads >>>( d_data1, d_data2, 
							    NX, NY,
							    o2.R, o2.R2);
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
      
      heat2d_ftcs_naive_order2_kernel<<< grid, threads >>>( d_data2, d_data1, 
							    NX, NY,
							    o2.R, o2.R2);   
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
      
    } else if (useOrder2b) { // use the 2nd order accurate scheme
      
      heat2d_ftcs_naive_order2b_kernel<<< grid, threads >>>( d_data1, d_data2, 
							     NX, NY,
							     o2.R, o2.R2b);
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
      
      heat2d_ftcs_naive_order2b_kernel<<< grid, threads >>>( d_data2, d_data1, 
							     NX, NY,
							     o2.R, o2.R2b);   
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");

    } else { // use the 4th order accurate scheme
      
      heat2d_ftcs_naive_order4_kernel<<< grid, threads >>>( d_data1, d_data2, 
							    NX, NY,
							    o4.S, o4.S2);
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
      
      heat2d_ftcs_naive_order4_kernel<<< grid, threads >>>( d_data2, d_data1, 
							     NX, NY,
							     o4.S, o4.S2);   
      // check if kernel execution generated an error
      getLastCudaError("Kernel execution failed");
    }

    /* save output (just for cross-checking, do not save when
       measuring computing time */
    if (ENABLE_GPU_SAVE) {

      if (iTime%T_OUTPUT == 0) {
	iOutput++;
	checkCudaErrors( cudaMemcpy( /* TODO */ ) );      
      }
      // PGM output
      if (SAVE_PGM and iTime%T_OUTPUT == 0)
	save_pgm(data1, "heat2d_gpu_naive_",iOutput,NX,NY);
      
      // MathGL save (3D view)
      if (SAVE_MGL and iTime%T_OUTPUT == 0)
	save_mgl(data1, "heat2d_gpu_naive_",iOutput,NX,NY);

      // VTK output
      if (SAVE_VTK and iTime%T_OUTPUT == 0)
	save_vtk(data1, "heat2d_gpu_naive_",iOutput);

      // HDF5 output
      if (SAVE_HDF5 and iTime%T_OUTPUT == 0)
	save_hdf5(data1, "heat2d_gpu_naive_",iOutput);

    }

  } // end for loop 
  
  // stop timer
  gpuTimer.stop();

  real_t gpu_time = gpuTimer.elapsed();
  printf( "GPU Processing time: %f (s)\n", gpu_time);
  
  // copy result from device to host
  real_t *resGPU = (real_t*) malloc( mem_size);
  checkCudaErrors( cudaMemcpy( /* TODO */ ) );
    
  if (SAVE_HDF5)
      write_xdmf_wrapper("heat2d_gpu_naive",N_ITER,T_OUTPUT);
 
  ////////////////////////////////////////////////////////
  // compute reference (CPU) solution to 2D heat equation
  // for performance comparison
  ////////////////////////////////////////////////////////
  initCondition2D (data1);
  initCondition2D (data2);
  
  Timer cpuTimer;
  cpuTimer.start();
  
  // time loop
  iTime=0;

  for (real_t t=0.0f; t<TMAX; t+=(2*DT), iTime+=2) {
    
    // compute next 2 time steps
    if (useOrder2) {
      heat2d_ftcs_cpu_order2( data1, data2);
      heat2d_ftcs_cpu_order2( data2, data1);
    } else if (useOrder2b) {
      heat2d_ftcs_cpu_order2b( data1, data2);
      heat2d_ftcs_cpu_order2b( data2, data1);
    } else {
      heat2d_ftcs_cpu_order4( data1, data2);
      heat2d_ftcs_cpu_order4( data2, data1);
    }
  }

  // stop timer
  cpuTimer.stop();
  real_t cpu_time = cpuTimer.elapsed();
  
  printf( "CPU Processing time: %g (s)\n", cpu_time);
  printf( "Speedup GPU/CPU : %f\n",cpu_time/gpu_time);

  printf("...comparing the results\n");
  double sum = 0, delta = 0;
  for(unsigned i = 0; i < NX*NY; i++){
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
