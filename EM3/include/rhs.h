#ifndef RHS_H
#define RHS_H

#include <cmath>
#include <iostream>
#include <time.h>
#include "derivs.h"
#include "parameters.h"
#include "profile_params.h"
#include "mathUtils.h"

#ifdef EM3_ENABLE_CUDA
#include "rhs_cuda.cuh"
    #include "params_cu.h"
    #include "profile_gpu.h"
#endif

/**@brief computes complete RHS iteratiing over all the blocks.
 * @param[out] unzipVarsRHS: unzipped variables computed RHS
 * @param[in]  unzipVars: unzipped variables. 
 * @param[in]  blkList: block list. 
 * @param[in]  numBlocks: number of blocks. 
 */
void em3rhs(double **uzipVarsRHS, const double **uZipVars, 
             const ot::Block* blkList, unsigned int numBlocks);

void em3rhs(double **uzipVarsRHS, const double **uZipVars,
             const unsigned int &offset,
             const double *ptmin, const double *ptmax, const unsigned int *sz,
             const unsigned int &bflag);

void em3_bcs(double *f_rhs, const double *f,
              const double *dxf, const double *dyf, const double *dzf,
              const double *pmin, const double *pmax,
              const double f_falloff, const double f_asymptotic,
              const unsigned int *sz, const unsigned int &bflag);

#endif
