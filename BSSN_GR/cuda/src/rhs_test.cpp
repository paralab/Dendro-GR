/**
 * @file rhs_test.cpp
 * @brief Host only rhs
 * @version 0.1
 * @date 2022-01-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <iostream>
#include "test_utils.h"
#include "point.h"
#include "block_gpu.h"
#include "parameters.h"

int main(int argc, char** argv)
{

  const unsigned int nblocks   = atoi(argv[1]);
  const unsigned int nReps     = atoi(argv[2]);
  const unsigned int bflag     = atoi(argv[3]);
  const unsigned int blk_sz_1d = 32;
  const unsigned int pencils   = 8;
  const Point grid_min(-0.5,-0.5,-0.5);
  const unsigned int sz[3]={blk_sz_1d,blk_sz_1d,blk_sz_1d};
  const Point dx(1.0/sz[0],1.0/sz[1],1.0/sz[2]);
  const unsigned int dof     = bssn::BSSN_NUM_VARS;
  

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
    blk[i].m_dx[0]  = dx.x() + 1e-4 * (i%3);
    blk[i].m_dx[1]  = dx.y() + 1e-4 * (i%3);
    blk[i].m_dx[2]  = dx.z() + 1e-4 * (i%3);

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
  
  profiler_t t_cpu;
  profiler_t t_gpu;
  unsigned int thread_qty = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
  unsigned int num_bssn_derivs=138;
  std::cout<<"OMP_NUM_THREADS="<<thread_qty<<std::endl;
  std::cout<<"blk sz (1d) : "<<blk_sz_1d<<std::endl;
  std::cout<<"num_iterations : "<<nReps<<std::endl;
  std::cout<<"pencils (gpu only) : "<<pencils<<std::endl;

  double * deriv_workspace = (double*)malloc(sizeof(double) * thread_qty * num_bssn_derivs * blk_sz_1d * blk_sz_1d * blk_sz_1d);
  
  
  bssnrhs_cpu<3,pencils,blk_sz_1d>(Fu_cpu,u,blk.data(),blk.size(),deriv_workspace, szp_dof, ko_sigma);
  t_cpu.start();
  for(unsigned int i=0; i < nReps; i++)
    bssnrhs_cpu<3,pencils,blk_sz_1d>(Fu_cpu,u,blk.data(),blk.size(),deriv_workspace, szp_dof, ko_sigma);
  t_cpu.stop();

  free(deriv_workspace);

  printf("cpu_time\tgpu_time\tspeedup\n");
  printf("%.4f\t%.4f\t%.2f\n",(double)t_cpu.seconds,(double)t_gpu.seconds,(double)t_cpu.seconds/(double)t_gpu.seconds);
  delete [] u;
  delete [] Fu_cpu;
  return 0;

}