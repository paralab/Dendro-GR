//
// Created by milinda on 7/25/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief This file contains all the parameters related to EM2 simulation.
 */
//

#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include <string.h>

#include <iostream>

#include "dendro.h"
#include "memory_pool.h"

namespace em2 {

extern mem::memory_pool<DendroScalar> EM2_MEM_POOL;

/** @brief These variable indexes are based on the variables defined in rkEM2.h
 */
enum VAR {
    U_E0    = 0,
    U_E1    = 1,
    U_E2    = 2,
    U_A0    = 3,
    U_A1    = 4,
    U_A2    = 5,
    U_PSI   = 6,
    U_GAMMA = 7
};

enum VAR_CONSTRAINT { C_DIVE = 0, C_DIVA = 1 };

/**@brief variable names. */
static const char* EM2_VAR_NAMES[] = {"U_E0", "U_E1", "U_E2",  "U_A0",
                                      "U_A1", "U_A2", "U_PSI", "U_GAMMA"};

static const char* EM2_CONSTRAINT_VAR_NAMES[] = {"C_DIVE", "C_DIVA"};

/**
 * @brief Refinement mode types.
 * WAMR : Wavelet based refinement.
 * EH : black hole event horizon based refinement.
 * EH_WAMR: both even horizon as well as WAMR based refinement.
 */
// this is not what is in NLSigma or BSSN_GR
// enum RefinementMode{WAMR=0, EH, EH_WAMR};
enum RefinementMode { WAMR = 0, FR, WAMR_FR };
// the following is what is in NLSigma:
// enum RefineMode {WAMR=0,FR,WAMR_FR};

/**@brief element order*/
extern unsigned int EM2_ELE_ORDER;

/**@brief assume padding begin*/
extern unsigned int EM2_PADDING_WIDTH;

/**@brief number of variables*/
static const unsigned int EM2_NUM_VARS            = 8;

/**@brief number of constraint variables*/
static const unsigned int EM2_CONSTRAINT_NUM_VARS = 2;

/***@brief number of RK45 stages*/
static const unsigned int EM2_RK45_STAGES         = 6;

/***@brief number of RK4 stages*/
static const unsigned int EM2_RK4_STAGES          = 4;

/**@brief number of rk4 stages*/
static const unsigned int EM2_RK3_STAGES          = 3;

/**@brief: parameter used for adaptive time step update. */
static const double EM2_SAFETY_FAC                = 0.8;

/**@brief number of internal variables*/
static const unsigned int EM2_NUM_VARS_INTENL =
    (EM2_RK45_STAGES + 1) * EM2_NUM_VARS;

/**@brief min domain add these to the parameter file.*/
extern double EM2_COMPD_MIN[3];
/**@brief min domain @todo add these to the parameter file. */
extern double EM2_COMPD_MAX[3];

/**@brief CFL stability number (specifies dt=EM2_CFL_FACTOR*dx)*/
extern double EM2_CFL_FACTOR;

/**@brief min coords of the OCTREE */
extern double EM2_OCTREE_MIN[3];
/**@brief max coords of the OCTREE */
extern double EM2_OCTREE_MAX[3];

/**@brief solution output frequency*/
extern unsigned int EM2_IO_OUTPUT_FREQ;

/**@brief timestep norms out put freq.*/
extern unsigned int EM2_TIME_STEP_OUTPUT_FREQ;

/**@brief remesh test frequency*/
extern unsigned int EM2_REMESH_TEST_FREQ;

/**@brief checkpoint store frequency*/
extern unsigned int EM2_CHECKPT_FREQ;

/**@brief restore the solver from check point if set to 1. */
extern unsigned int EM2_RESTORE_SOLVER;

/**@brief use the block adaptivity and disable the AMR*/
extern unsigned int EM2_ENABLE_BLOCK_ADAPTIVITY;

/**@brief file prefix for VTU*/
extern std::string EM2_VTU_FILE_PREFIX;

/**@brief file prefix for write check point*/
extern std::string EM2_CHKPT_FILE_PREFIX;

/**@brief file prefix to write profile info.*/
extern std::string EM2_PROFILE_FILE_PREFIX;

/**@brief number of refine variables*/
extern unsigned int EM2_NUM_REFINE_VARS;

/**@brief indices of refine var ids*/
extern unsigned int EM2_REFINE_VARIABLE_INDICES[EM2_NUM_VARS];

/**@brief number of evolution variables written to vtu files*/
extern unsigned int EM2_NUM_EVOL_VARS_VTU_OUTPUT;

/**@brief number of constrint variables written to vtu files*/
extern unsigned int EM2_NUM_CONST_VARS_VTU_OUTPUT;

/**@brief evolution variable IDs written to vtu files*/
extern unsigned int EM2_VTU_OUTPUT_EVOL_INDICES[EM2_NUM_VARS];

/**@brief constraint variable IDs written to vtu files*/
extern unsigned int EM2_VTU_OUTPUT_CONST_INDICES[EM2_CONSTRAINT_NUM_VARS];

/**@brief solution output gap (instead of freq. we can use to output the
 * solution if currentTime > lastIOOutputTime + EM2_IO_OUTPUT_GAP)*/
extern double EM2_IO_OUTPUT_GAP;

/**@brief prefered grain sz to use when selecting active npes*/
extern unsigned int EM2_DENDRO_GRAIN_SZ;

/**@brief AMR coarsening factor (we coarsen if
 * tol<EM2_DENDRO_AMR_FAC*EM2_WAVELET_TOL)*/
extern double EM2_DENDRO_AMR_FAC;

/**@brief wavelet tolerance value. */
extern double EM2_WAVELET_TOL;
/**@brief load-imbalance tolerance value. */
extern double EM2_LOAD_IMB_TOL;
/**@brief: Splitter fix value*/
extern unsigned int EM2_SPLIT_FIX;

/**@brief: async. communication at a time. (upper bound shoud be EM2_NUM_VARS)
 */
extern unsigned int EM2_ASYNC_COMM_K;

/**@brief simulation begin time. */
extern double EM2_RK45_TIME_BEGIN;
/**@brief simulation end time*/
extern double EM2_RK45_TIME_END;
/**@brief rk time step size. */
extern double EM2_RK45_TIME_STEP_SIZE;

/** desired tolerance value for the rk45 method (adaptive time stepping. )*/
extern double EM2_RK45_DESIRED_TOL;

/**@brief rk method type*/
// this one is not in nlsm
extern unsigned int EM2_RK_TYPE;

/**@brief BBH initial data type */
extern unsigned int EM2_ID_TYPE;

/**@brief physical coordinates for grid, x_min */
extern double EM2_GRID_MIN_X;

/**@brief physical coordinates for grid, x_max */
extern double EM2_GRID_MAX_X;

/**@brief physical coordinates for grid, y_min */
extern double EM2_GRID_MIN_Y;

/**@brief physical coordinates for grid, y_max */
extern double EM2_GRID_MAX_Y;

/**@brief physical coordinates for grid, z_min */
extern double EM2_GRID_MIN_Z;

/**@brief physical coordinates for grid, z_max */
extern double EM2_GRID_MAX_Z;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double EM2_BLK_MIN_X;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double EM2_BLK_MIN_Y;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double EM2_BLK_MIN_Z;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double EM2_BLK_MAX_X;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double EM2_BLK_MAX_Y;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double EM2_BLK_MAX_Z;

/**@brief: dimension of the grid*/
extern unsigned int EM2_DIM;

/**@brief: max refinement level*/
extern unsigned int EM2_MAXDEPTH;

/**@brief: Kreiss-Oliger dissipation */
extern double KO_DISS_SIGMA;

/**@brief: Initial data Gaussian amplitude */
extern double EM2_ID_AMP1;

/**@brief: Initial data Gaussian amplitude */
extern double EM2_ID_LAMBDA1;

/**@brief: Initial data Gaussian amplitude */
extern double EM2_ID_AMP2;

/**@brief: Initial data Gaussian width */
// extern double EM2_ID_DELTA1;

/**@brief: Initial data Gaussian width */
// extern double EM2_ID_DELTA2;

/**@brief: Initial data Gaussian x offset */
// extern double EM2_ID_XC1;
// extern double EM2_ID_YC1;
// extern double EM2_ID_ZC1;

/**@brief: Initial data Gaussian x offset */
// extern double EM2_ID_XC2;
// extern double EM2_ID_YC2;
// extern double EM2_ID_ZC2;

/**@brief: Initial data Gaussian elliptic x factor */
// extern double EM2_ID_EPSX1;

/**@brief: Initial data Gaussian elliptic y factor */
// extern double EM2_ID_EPSY1;

/**@brief: Initial data Gaussian elliptic z factor */
// extern double EM2_ID_EPSZ1;

/**@brief: Initial data Gaussian elliptic x factor */
// extern double EM2_ID_EPSX2;

/**@brief: Initial data Gaussian elliptic y factor */
// extern double EM2_ID_EPSY2;

/**@brief: Initial data Gaussian elliptic z factor */
// extern double EM2_ID_EPSZ2;

/**@brief: Initial data Gaussian R */
// extern double EM2_ID_R1;

/**@brief: Initial data Gaussian R */
// extern double EM2_ID_R2;

/**@brief: Initial data Gaussian nu */
// extern double EM2_ID_NU1;

/**@brief: Initial data Gaussian nu */
// extern double EM2_ID_NU2;

/**@brief: Initial data Gaussian Omega */
// extern double EM2_ID_OMEGA;

/**@brief: waves force refinement threshold for */
extern double EM2_CHI_REFINE_VAL;

/**@brief: waves force coarsen threshold for */
extern double EM2_CHI_COARSEN_VAL;

///**@brief: waves specify the refinement mode.  */
// extern RefinementMode EM2_REFINEMENT_MODE;

/**@brief: dissipation type */
// this one is not in nlsm
extern unsigned int DISSIPATION_TYPE;

/**@brief: Kreiss-Oliger dissipation */
extern double KO_DISS_SIGMA;

/**@brief: EM2_USE_WAVELET_TOL_FUNCTION */
extern unsigned int EM2_USE_WAVELET_TOL_FUNCTION;

/**@brief: EM2_WAVELET_TOL_FUNCTION_R0 */
extern double EM2_WAVELET_TOL_FUNCTION_R0;

/**@brief: EM2_WAVELET_TOL_FUNCTION_R0 */
extern double EM2_WAVELET_TOL_FUNCTION_R1;

/**@brief: EM2_WAVELET_TOL_MAX */
extern double EM2_WAVELET_TOL_MAX;

/**@brief: @david can you please add some comeents for these parameters. */
extern unsigned int EM2_DISSIPATION_NC;

/**@brief: @david can you please add some comeents for these parameters. */
extern unsigned int EM2_DISSIPATION_S;

/**@brief: if true it will use finite differnce like grid transfer*/
extern bool EM2_USE_FD_GRID_TRANSFER;

/**@brief: tolerance for refinement based on EH */
extern double EM2_EH_REFINE_VAL;

/**@brief: tolerance for coarsen based on EH */
extern double EM2_EH_COARSEN_VAL;

/**@brief: refinement mode for the application*/
extern RefinementMode EM2_REFINEMENT_MODE;

///**@brief: if true output only the z slice*/
extern bool EM2_VTU_Z_SLICE_ONLY;

}  // namespace em2

#endif  // SFCSORTBENCH_PARAMETERS_H
