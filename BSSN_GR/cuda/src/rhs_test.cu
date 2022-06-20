/**
 * @brief CPU-GPU rhs comparison. 
*/
#include <iostream>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include "device.h"
#include "point.h"
#include "derivs_cu.cuh"
#include "derivs.h"
#include "block_gpu.h"
#include "profiler.h"
#include "grDef.h"
#include "grUtils.h"
#include "rhs.h"
#include "test_utils.h"
#include "bssn_kernels.cuh"


int main(int argc, char** argv)
{

  const unsigned int nblocks   = atoi(argv[1]);
  const unsigned int nReps     = atoi(argv[2]);
  const unsigned int bflag     = atoi(argv[3]);

  const unsigned int blk_sz_1d = 13;
  const unsigned int pencils   = 13;
  const Point grid_min(-0.5,-0.5,-0.5);
  const unsigned int sz[3]={blk_sz_1d,blk_sz_1d,blk_sz_1d};
  const Point dx(1.0e-2/sz[0],1.0e-2/sz[1],1.0e-2/sz[2]);
  const unsigned int dof     = bssn::BSSN_NUM_VARS;

  int devicesCount;
  cudaGetDeviceCount(&devicesCount);
  printf("number of cuda devices: %d\n",devicesCount);

  cudaSetDevice(0);
  cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
  
  std::vector<BlockGPU3D> blk;
  blk.resize(nblocks);

  blk[0].m_offset=0;
  for(unsigned int i=0; i < blk.size(); i++)
  {
    blk[i].m_sz[0]  = sz[0];
    blk[i].m_sz[1]  = sz[1];
    blk[i].m_sz[2]  = sz[2];
    blk[i].m_aligned_sz[0] = sz[0];
    blk[i].m_aligned_sz[1] = sz[1];
    blk[i].m_aligned_sz[2] = sz[2];
    blk[i].m_dx[0]  = dx.x()/(i%3 + 1);
    blk[i].m_dx[1]  = dx.y()/(i%3 + 1);
    blk[i].m_dx[2]  = dx.z()/(i%3 + 1);

    blk[i].m_ptMin[0] = grid_min.x();  //+ 0.12* (i%10);
    blk[i].m_ptMin[1] = grid_min.y();  //+ 0.12* (i%10);
    blk[i].m_ptMin[2] = grid_min.z();  //+ 0.12* (i%10);
    blk[i].m_bflag  = bflag;


    if( i > 0 )
      blk[i].m_offset = blk[i-1].m_offset + blk[i-1].m_aligned_sz[0] * blk[i-1].m_aligned_sz[1] * blk[i-1].m_aligned_sz[2];

  }
  const double ko_sigma = 4e-1;
  const unsigned int szp_dof    = (blk.back().m_offset +  blk.back().m_aligned_sz[0] * blk.back().m_aligned_sz[1] * blk.back().m_aligned_sz[2]);
  const unsigned int blk_vec_sz = (blk.back().m_offset +  blk.back().m_aligned_sz[0] * blk.back().m_aligned_sz[1] * blk.back().m_aligned_sz[2]) * dof;
  double * u          = create_block<double>(blk.data(),blk.size(),dof);
  double * Fu_cpu     = new double[blk_vec_sz];
  double * Fu_gpu     = new double[blk_vec_sz];

  double * dptr_u  = GPUDevice::device_malloc<double>(blk_vec_sz);
  double * dptr_Fu = GPUDevice::device_malloc<double>(blk_vec_sz);

  BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(nblocks);
  GPUDevice::host_to_device(blk.data(),dptr_blk,nblocks);
  GPUDevice::host_to_device(u,dptr_u,blk_vec_sz);

  profiler_t t_cpu;
  profiler_t t_gpu;

  const unsigned int thread_qty = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
  const unsigned int num_bssn_derivs=210;
  const unsigned int blk_batch_size =1;
  std::cout<<"OMP_NUM_THREADS="<<thread_qty<<std::endl;
  std::cout<<"nblocks="<<nblocks<<std::endl;
  double * deriv_workspace = (double*)malloc(sizeof(double) * thread_qty * num_bssn_derivs * blk_sz_1d * blk_sz_1d * blk_sz_1d);
  bssnrhs_cpu<3,pencils,blk_sz_1d>(Fu_cpu,u,blk.data(),blk.size(),deriv_workspace, szp_dof, ko_sigma);

  
  #pragma nounroll
  for(unsigned int i=0; i < nReps; i++){
    t_cpu.start();
    bssnrhs_cpu<3,pencils,blk_sz_1d>(Fu_cpu,u,blk.data(),blk.size(),deriv_workspace, szp_dof, ko_sigma);
    t_cpu.stop();
    //printf("iter = %d time = %f snap= %f\n",i, (double)t_cpu.seconds, (double)t_cpu.snap);
    // Fu_cpu[0]+=0;
    // u[0]+=1;
  }
  
  free(deriv_workspace);

  // Print device and precision
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("\nDevice Name: %s\n", prop.name);
  printf("Compute Capability: %d.%d\n\n", prop.major, prop.minor);

  
  const unsigned int nstreams             = 1;
  BSSN_EVAR_DERIVS* deriv_evars_ws        = new BSSN_EVAR_DERIVS[nstreams];
  BSSN_EVAR_DERIVS* dptr_deriv_evars_ws   = GPUDevice::device_malloc<BSSN_EVAR_DERIVS>(nstreams);
  const unsigned int BLK_SZ               = blk_sz_1d * blk_sz_1d * blk_sz_1d;
  const unsigned int sz_p_stream          = num_bssn_derivs * blk_batch_size * BLK_SZ;
  DEVICE_REAL* deriv_base_ws              = GPUDevice::device_malloc<DEVICE_REAL>(nstreams * sz_p_stream);
  
  for(unsigned int i=0; i < nstreams; i++)
  { 
    const unsigned int BATCHED_BLOCKS_SZ = blk_batch_size;
    DEVICE_REAL* deriv_base = deriv_base_ws  + i * sz_p_stream;
    BSSN_EVAR_DERIVS* deriv_evars = deriv_evars_ws + i;
    #include "../scripts/bssnrhs_memalloc.cuh"
  }
  GPUDevice::host_to_device(deriv_evars_ws,dptr_deriv_evars_ws,nstreams);
  // device::bssnrhs_gpu<3,pencils,blk_sz_1d, blk_batch_size, blk_sz_1d>(dptr_Fu,dptr_u,dptr_blk, blk.size(),dptr_deriv_evars, BLK_SZ, ko_sigma);
  // cudaDeviceSynchronize();
  // cudaFuncSetCacheConfig(&device::bssnrhs<3,pencils,blk_sz_1d, blk_sz_1d>, cudaFuncCachePreferL1);
  // dim3 gb = dim3(blk_batch_size, 1, 1);
  // dim3 tb = dim3(blk_sz_1d, blk_sz_1d, 1);
  // device::bssnrhs<3,pencils,blk_sz_1d, blk_sz_1d><<<gb,tb>>>(dptr_Fu,dptr_u,dptr_blk, blk.size(),dptr_deriv_evars, szp_dof, ko_sigma);
  // cudaDeviceSynchronize();

  // t_gpu.start();
  // for(unsigned int i=0; i < nReps; i++){
  //   device::bssnrhs<3,pencils,blk_sz_1d, blk_sz_1d><<<gb,tb>>>(dptr_Fu,dptr_u,dptr_blk, blk.size(),dptr_deriv_evars, szp_dof, ko_sigma);
  //   cudaDeviceSynchronize();
  // }
  // t_gpu.stop();
  int maxbytes = 13*13*13 *8 * 3; // 96 KB
  cudaFuncSetAttribute(device::eval_rhs3<3,blk_sz_1d, blk_sz_1d>, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
  //cudaFuncSetAttribute(device::eval_rhs4<3,blk_sz_1d, blk_sz_1d>, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
  //cudaFuncSetCacheConfig(&device::eval_rhs3<3,blk_sz_1d, blk_sz_1d>, cudaFuncCachePreferEqual);
  device::bssnrhs3<3, blk_sz_1d, blk_sz_1d, blk_batch_size,nstreams> (dptr_Fu,dptr_u,dptr_blk, blk.size(), szp_dof, ko_sigma);
  cudaDeviceSynchronize();
  GPUDevice::check_last_error();

  
  for(unsigned int i=0; i < nReps; i++){
    t_gpu.start();
    device::bssnrhs3<3, blk_sz_1d, blk_sz_1d, blk_batch_size, nstreams> (dptr_Fu,dptr_u,dptr_blk, blk.size(), szp_dof, ko_sigma);
    cudaDeviceSynchronize();
    t_gpu.stop();
  }

  // device::bssnrhs4<3, blk_sz_1d, blk_sz_1d, blk_batch_size,nstreams> (dptr_Fu,dptr_u,dptr_blk, blk.size(), dptr_deriv_evars_ws , szp_dof, ko_sigma);
  // cudaDeviceSynchronize();
  // GPUDevice::check_last_error();

  
  // for(unsigned int i=0; i < nReps; i++){
  //   t_gpu.start();
  //   device::bssnrhs4<3, blk_sz_1d, blk_sz_1d, blk_batch_size, nstreams> (dptr_Fu,dptr_u,dptr_blk, blk.size(),dptr_deriv_evars_ws, szp_dof, ko_sigma);
  //   cudaDeviceSynchronize();
  //   t_gpu.stop();
  // }


  // device::bssnrhs2<3, blk_sz_1d, blk_sz_1d, blk_batch_size,nstreams> (dptr_Fu,dptr_u,dptr_blk, blk.size(), dptr_deriv_evars_ws, szp_dof, ko_sigma);
  // cudaDeviceSynchronize();
  // GPUDevice::check_last_error();

  
  // for(unsigned int i=0; i < nReps; i++){
  //   t_gpu.start();
  //   device::bssnrhs2<3, blk_sz_1d, blk_sz_1d, blk_batch_size, nstreams> (dptr_Fu,dptr_u,dptr_blk, blk.size(), dptr_deriv_evars_ws, szp_dof, ko_sigma);
  //   cudaDeviceSynchronize();
  //   t_gpu.stop();
  // }

  


  
  // dim3 gb = dim3(1, blk_sz_1d, blk_batch_size);
  // dim3 tb = dim3(blk_sz_1d, blk_sz_1d, 1);
  // {
  //   const unsigned int nblocks = blk.size();
  //   void* kernelArgs[] = {&dptr_Fu, &dptr_u, &dptr_blk, (void * )&nblocks, &dptr_deriv_evars, (void *)&BLK_SZ,  (void * )&ko_sigma/* add kernel args */ };
  //   //device::bssnrhs1<3,pencils,blk_sz_1d, blk_batch_size, blk_sz_1d><<<gb,tb>>>(dptr_Fu,dptr_u,dptr_blk, blk.size(),dptr_deriv_evars, BLK_SZ, ko_sigma);
  //   cudaLaunchCooperativeKernel((void*)device::bssnrhs1<3,pencils,blk_sz_1d, blk_batch_size, blk_sz_1d>, gb, tb, kernelArgs);
  //   GPUDevice::check_last_error();
  //   cudaDeviceSynchronize();
  // }
  
  // t_gpu.start();
  // for(unsigned int i=0; i < nReps; i++)
  // {
  //   const unsigned int nblocks = blk.size();
  //   void* kernelArgs[] = {&dptr_Fu, &dptr_u, &dptr_blk, (void * )&nblocks, &dptr_deriv_evars, (void *)&BLK_SZ,  (void * )&ko_sigma/* add kernel args */ };
  //   //device::bssnrhs1<3,pencils,blk_sz_1d, blk_batch_size, blk_sz_1d><<<gb,tb>>>(dptr_Fu,dptr_u,dptr_blk, blk.size(),dptr_deriv_evars, BLK_SZ, ko_sigma);
  //   cudaLaunchCooperativeKernel((void*)device::bssnrhs1<3,pencils,blk_sz_1d, blk_batch_size, blk_sz_1d>, gb, tb, kernelArgs);
  //   GPUDevice::check_last_error();
  //   cudaDeviceSynchronize();
  // }
  // t_gpu.stop();

  GPUDevice::device_to_host(Fu_gpu, dptr_Fu, blk_vec_sz);

  printf("cpu_time\tgpu_time\tspeedup\n");
  printf("%.4f\t%.4f\t%.2f\n",(double)t_cpu.seconds,(double)t_gpu.seconds,(double)t_cpu.seconds/(double)t_gpu.seconds);

  GPUDevice::device_free(deriv_base_ws);
  GPUDevice::device_free(dptr_deriv_evars_ws);
  delete [] deriv_evars_ws;
  GPUDevice::device_free(dptr_u);
  GPUDevice::device_free(dptr_Fu);
  
  const unsigned int PW=3;
  double l_inf_error[dof];
  for(unsigned int v=0; v < dof; v++)
    l_inf_error[v]=0.0;

  for(unsigned int blk_id =0; blk_id < blk.size(); blk_id++)
  {
    const unsigned int offset  = blk[blk_id].m_offset; 
    const unsigned int a_nx      = blk[blk_id].m_sz[0];
    const unsigned int a_ny      = blk[blk_id].m_sz[1];
    const unsigned int a_nz      = blk[blk_id].m_sz[2];

    const unsigned int nx      = blk[blk_id].m_aligned_sz[0];
    const unsigned int ny      = blk[blk_id].m_aligned_sz[1];
    const unsigned int nz      = blk[blk_id].m_aligned_sz[2];

    const unsigned int blkSz   = nx * ny * nz; 

    for(unsigned int k=PW; k < a_nz-PW; k++)
    for(unsigned int j=PW; j < a_ny-PW; j++)
    for(unsigned int i=PW; i < a_nz-PW; i++)
    for(unsigned int v=0; v < dof; v++)
    {
      const unsigned int pp     = v * szp_dof  + offset + k * nx * ny + j * nx + i;
      const double error_metric = fabs((Fu_cpu[pp]-Fu_gpu[pp])/Fu_cpu[pp]);
      if (l_inf_error[v] < error_metric )
      {
        l_inf_error[v]=error_metric;
        // if(v== bssn::VAR::U_SYMGT0)
        //   printf("idx: %d\n ", k*nx*ny + j*nx + i);
      }
        
    }

  }

  for(unsigned int v=0; v < dof; v++)
    std::cout<<"l_inf var "<<bssn::BSSN_VAR_NAMES[v]<<" = "<<l_inf_error[v]<<std::endl;
  
  
  delete [] u;
  delete [] Fu_cpu;
  delete [] Fu_gpu;

  return 0;
}