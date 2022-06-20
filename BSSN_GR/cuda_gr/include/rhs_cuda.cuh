// * Created by milinda on 8/10/18
/**
 * @author Milinda Fernando.
 * School of Computing, University of Utah.
 * @brief Computation of the rhs using cuda
 * Created by milinda on 8/10/18
 * */

#ifndef SFCSORTBENCH_CUDARHS_H
#define SFCSORTBENCH_CUDARHS_H

#include "cuda_runtime.h"
#include <device_launch_parameters.h>
#include "block_cu.h"
#include "bssn_rhs_deriv_mem_cuda.h"
#include "cudaUtils.h"
#include "params_cu.h"
#include "rhs_bssn.cuh"
#include "profile_gpu.h"
#include "cudaUtils.cuh"
#include "gpuTest.cuh"
#include "parameters.h"
#include "grDef.h"
#include "rhs.h"



namespace cuda
{


    /***
     * @brief performs kernel pre-launch tasks and launch the bssnrhs kernel
     *
     **/
     void computeRHS(double **unzipVarsRHS, const double **uZipVars,const ot::Block* dendroBlockList,unsigned int numBlocks,const cuda::BSSNComputeParams* bssnPars,dim3 blockDim,const Point & pt_min, const Point & pt_max,unsigned int numStreams,unsigned int device=0);




}

#endif //SFCSORTBENCH_CUDARHS_H
