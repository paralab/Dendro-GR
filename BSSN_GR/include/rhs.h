#ifndef RHS_H
#define RHS_H

#include <time.h>

#include <cmath>
#include <iostream>

#include "block.h"
#include "derivs.h"
#include "grDef.h"
#include "mathUtils.h"
#include "parameters.h"
#include "profile_params.h"

#ifdef BSSN_ENABLE_CUDA
#include "params_cu.h"
#include "profile_gpu.h"
#include "rhs_cuda.cuh"
#endif

/**@brief computes complete RHS iteratiing over all the blocks.
 * @param[out] unzipVarsRHS: unzipped variables computed RHS
 * @param[in]  unzipVars: unzipped variables.
 * @param[in]  blkList: block list.
 * @param[in]  numBlocks: number of blocks.
 */
void bssnRHS(double **uzipVarsRHS, const double **uZipVars,
             const ot::Block *blkList, unsigned int numBlocks,
             const double curr_time, const double **uZipConstVars);


void bssnrhs(double **uzipVarsRHS, const double **uZipVars,
             const unsigned int &offset, const double *ptmin,
             const double *ptmax, const unsigned int *sz,
             const unsigned int &bflag, const double t,
             const double **uZipConstVars);

// void bssnrhs_sep(double **uzipVarsRHS, const double **uZipVars,
//              const unsigned int &offset,
//              const double *ptmin, const double *ptmax, const unsigned int
//              *sz, const unsigned int &bflag);

void bssn_bcs(double *f_rhs, const double *f, const double *dxf,
              const double *dyf, const double *dzf, const double *pmin,
              const double *pmax, const double f_falloff,
              const double f_asymptotic, const unsigned int *sz,
              const unsigned int &bflag);

void freeze_bcs(double *f_rhs, const unsigned int *sz,
                const unsigned int &bflag);

void fake_initial_data(double x, double y, double z, double *u);

void max_spacetime_speeds(double *const lambda1max, double *const lambda2max,
                          double *const lambda3max, const double *const alpha,
                          const double *const beta1, const double *const beta2,
                          const double *const beta3, const double *const gtd11,
                          const double *const gtd12, const double *const gtd13,
                          const double *const gtd22, const double *const gtd23,
                          const double *const gtd33, const double *const chi,
                          const unsigned int *sz);

void call_HAD_rhs();

#endif
