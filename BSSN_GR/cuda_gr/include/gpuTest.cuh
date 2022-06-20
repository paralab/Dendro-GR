//
// Created by milinda on 9/18/18.
//
/**
 * @brief constains some gpu test functions for computations.
 * @author Milinda Fernando, School of Computing, University of Utah
 *
 * */

#ifndef DENDRO_5_0_GPUTEST_H
#define DENDRO_5_0_GPUTEST_H

#include<iostream>
#include"cuda_runtime.h"
#include<device_launch_parameters.h>
#include "block_cu.h"
#include "params_cu.h"
#include "bssn_rhs_deriv_mem_cuda.h"
#include "cudaUtils.cuh"
#include "derivs.cuh"
#include "derivs_bssn.cuh"
#include "cudaUtils.h"

namespace cuda
{
    __global__ void __compute1D_derivs(double** __unzipOutVar, const double**__unzipInVar, MemoryDerivs* __derivWorkspace, const cuda::_Block* __dendroBlkList, const BSSNComputeParams* __bssnPar, const cudaDeviceProp*__deviceProperties);


}



#endif //DENDRO_5_0_GPUTEST_H
