#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include <string.h>

#include <iostream>

namespace maxwell {

static const unsigned int MAXWELL_ELE_ORDER   = 4;

static const unsigned int MAXWELL_RK45_STAGES = 6;

static const unsigned int MAXWELL_RK4_STAGES  = 4;

static const unsigned int MAXWELL_RK3_STAGES  = 3;

static const double MAXWELL_SAFETY_FAC        = 0.8;

extern unsigned int MAXWELL_TIME_STEP_OUTPUT_FREQ;
extern unsigned int MAXWELL_SPLIT_FIX;
extern double MAXWELL_COMPD_MIN[3];

extern double MAXWELL_COMPD_MAX[3];

extern double MAXWELL_OCTREE_MIN[3];

extern double MAXWELL_OCTREE_MAX[3];

static const unsigned int MAXWELL_NUM_VARS = 8;

extern unsigned int MAXWELL_RESTORE_SOLVER;

extern unsigned int MAXWELL_IO_OUTPUT_FREQ;

extern unsigned int MAXWELL_REMESH_TEST_FREQ;

extern unsigned int MAXWELL_CHECKPT_FREQ;

extern double MAXWELL_IO_OUTPUT_GAP;

extern std::string MAXWELL_VTU_FILE_PREFIX;

extern std::string MAXWELL_CHKPT_FILE_PREFIX;

extern std::string MAXWELL_PROFILE_FILE_PREFIX;

extern unsigned int MAXWELL_NUM_EVOL_VARS_VTU_OUTPUT;

extern unsigned int MAXWELL_VTU_OUTPUT_EVOL_INDICES[MAXWELL_NUM_VARS];

extern unsigned int MAXWELL_DENDRO_GRAIN_SZ;

extern unsigned int MAXWELL_ASYNC_COMM_K;

extern double MAXWELL_DENDRO_AMR_FAC;

extern double MAXWELL_LOAD_IMB_TOL;

extern unsigned int MAXWELL_DIM;

extern unsigned int MAXWELL_MAXDEPTH;

extern double MAXWELL_WAVELET_TOL;

extern unsigned int MAXWELL_NUM_REFINE_VARS;

extern unsigned int MAXWELL_REFINE_VARIABLE_INDICES[MAXWELL_NUM_VARS];

extern double MAXWELL_RK45_TIME_BEGIN;

extern double MAXWELL_RK45_TIME_END;

extern double MAXWELL_RK45_TIME_STEP_SIZE;

extern double MAXWELL_RK45_DESIRED_TOL;

extern double KO_DISS_SIGMA;

extern unsigned int MAXWELL_ENABLE_BLOCK_ADAPTIVITY;

extern double MAXWELL_BLK_MIN_X;

extern double MAXWELL_BLK_MIN_Y;

extern double MAXWELL_BLK_MIN_Z;

extern double MAXWELL_BLK_MAX_X;

extern double MAXWELL_BLK_MAX_Y;

extern double MAXWELL_BLK_MAX_Z;

extern double MAXWELL_GRID_MIN_X;

extern double MAXWELL_GRID_MAX_X;

extern double MAXWELL_GRID_MIN_Y;

extern double MAXWELL_GRID_MAX_Y;

extern double MAXWELL_GRID_MIN_Z;

extern double MAXWELL_GRID_MAX_Z;

static const double MAXWELL_CFL_FACTOR = 0.1;

extern unsigned int MAXWELL_ID_TYPE;

static const unsigned int MAXWELL_NUM_VARS_INTENL =
    (MAXWELL_RK45_STAGES + 1) * MAXWELL_NUM_VARS;

}  // namespace maxwell

#endif  // SFCSORTBENCH_PARAMETERS_H
