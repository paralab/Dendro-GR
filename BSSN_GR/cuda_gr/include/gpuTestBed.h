//
// Created by milinda on 9/19/18.
//

/**
 * @brief This is a simple test bed to figure out the best configurations for the GPU computation. Bench mark GPU for different operations.
 * @author Milinda Fernando.
 * School of Computing, University of Utah.
 * */
#ifndef DENDRO_5_0_GPUTESTBED_H
#define DENDRO_5_0_GPUTESTBED_H

#include <iostream>
#include <vector>
#include "block.h"
#include "grUtils.h"
#include "rhs.h"
#include "parameters.h"
#include "profile_params.h"
#include "mpi.h"
#include <chrono>
#include "random"


typedef std::chrono::high_resolution_clock Time;
typedef std::chrono::milliseconds ms;
typedef std::chrono::duration<float> fsec;


#ifdef BSSN_ENABLE_CUDA
    #include "rhs_cuda.cuh"
    #include "params_cu.h"
    #include "profile_gpu.h"
    #include "gpuTest.cuh"
#endif






#endif //DENDRO_5_0_GPUTESTBED_H
