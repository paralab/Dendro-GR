//
// Created by milinda on 8/9/18.
//

#include "params_cu.h"

namespace cuda {
/**stores the device properties*/
cudaDeviceProp* __CUDA_DEVICE_PROPERTIES = NULL;

/**pointer to the block list*/
_Block* __DENDRO_BLOCK_LIST              = NULL;

/** number of blocks*/
unsigned int* __DENDRO_NUM_BLOCKS        = NULL;

unsigned int* __GPU_BLOCK_MAP            = NULL;

unsigned int* __NUM_GPU_BLOCKS           = NULL;

/** number of evol vars */
unsigned int* __BSSN_NUM_VARS            = NULL;

/** number of constraint vars */
unsigned int* __BSSN_CONSTRAINT_NUM_VARS = NULL;

/** x% of block shared memory utilised for the bssn computations*/
double* __GPU_BLOCK_SHARED_MEM_UTIL      = NULL;

/**max block size*/
unsigned int* __DENDRO_BLK_MAX_SZ        = NULL;

MemoryDerivs* __BSSN_DERIV_WORKSPACE     = NULL;

/**unzip input */
double** __UNZIP_INPUT                   = NULL;

/**unzip output*/
double** __UNZIP_OUTPUT                  = NULL;

/**bssn compute parameters needed bssn equations*/
BSSNComputeParams* __BSSN_COMPUTE_PARMS  = NULL;

}  // end of namespace cuda
