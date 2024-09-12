#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include <string.h>

#include <iostream>

namespace fluid {

static const unsigned int FLUID_ELE_ORDER   = 4;

static const unsigned int FLUID_RK45_STAGES = 6;

static const unsigned int FLUID_RK4_STAGES  = 4;

static const unsigned int FLUID_RK3_STAGES  = 3;

static const double FLUID_SAFETY_FAC        = 0.8;

extern unsigned int FLUID_TIME_STEP_OUTPUT_FREQ;
extern unsigned int FLUID_SPLIT_FIX;
extern double FLUID_COMPD_MIN[3];

extern double FLUID_COMPD_MAX[3];

extern double FLUID_OCTREE_MIN[3];

extern double FLUID_OCTREE_MAX[3];

extern double FLUID_VACUUM_RESET;
extern double FLUID_VACUUM_D_RESET;
extern double FLUID_VACUUM_TAU_RESET;
static const unsigned int FLUID_NUM_VARS = 5;

extern unsigned int FLUID_RESTORE_SOLVER;

extern unsigned int FLUID_IO_OUTPUT_FREQ;

extern unsigned int FLUID_REMESH_TEST_FREQ;

extern unsigned int FLUID_CHECKPT_FREQ;

extern double FLUID_IO_OUTPUT_GAP;

extern std::string FLUID_VTU_FILE_PREFIX;

extern std::string FLUID_CHKPT_FILE_PREFIX;

extern std::string FLUID_PROFILE_FILE_PREFIX;

extern unsigned int FLUID_NUM_EVOL_VARS_VTU_OUTPUT;

extern unsigned int FLUID_VTU_OUTPUT_EVOL_INDICES[FLUID_NUM_VARS];

extern unsigned int FLUID_DENDRO_GRAIN_SZ;

extern unsigned int FLUID_ASYNC_COMM_K;

extern double FLUID_DENDRO_AMR_FAC;

extern double FLUID_LOAD_IMB_TOL;

extern unsigned int FLUID_DIM;

extern unsigned int FLUID_MAXDEPTH;

extern double FLUID_WAVELET_TOL;

extern unsigned int FLUID_NUM_REFINE_VARS;

extern unsigned int FLUID_REFINE_VARIABLE_INDICES[FLUID_NUM_VARS];

extern double FLUID_RK45_TIME_BEGIN;

extern double FLUID_RK45_TIME_END;

extern double FLUID_RK45_TIME_STEP_SIZE;

extern double FLUID_RK45_DESIRED_TOL;

extern double KO_DISS_SIGMA;

extern unsigned int FLUID_ENABLE_BLOCK_ADAPTIVITY;

extern double FLUID_BLK_MIN_X;

extern double FLUID_BLK_MIN_Y;

extern double FLUID_BLK_MIN_Z;

extern double FLUID_BLK_MAX_X;

extern double FLUID_BLK_MAX_Y;

extern double FLUID_BLK_MAX_Z;

extern double FLUID_GRID_MIN_X;

extern double FLUID_GRID_MAX_X;

extern double FLUID_GRID_MIN_Y;

extern double FLUID_GRID_MAX_Y;

extern double FLUID_GRID_MIN_Z;

extern double FLUID_GRID_MAX_Z;

static const double FLUID_CFL_FACTOR = 0.1;

extern unsigned int FLUID_ID_TYPE;

extern unsigned int FLUID_COORDS;

extern unsigned int FLUID_RECON_METHOD;

extern double FLUID_GAMMA;

extern double FLUID_VACUUM;

extern double FLUID_VACUUM_D;

extern double FLUID_VACUUM_TAU;

extern unsigned int FLUID_EOS_SOLVER;

extern unsigned int FLUID_USE_WAVE_SPEEDS;

extern unsigned int FLUID_CONTOPRIMWARN;

extern unsigned int FLUID_AVERAGING_FREQUENCY;

extern double FLUID_UNITS_CLU;

extern double FLUID_UNITS_CMU;

extern double FLUID_UNITS_CKB;

static const unsigned int FLUID_NUM_VARS_INTENL =
    (FLUID_RK45_STAGES + 1) * FLUID_NUM_VARS;

}  // namespace fluid

#endif  // SFCSORTBENCH_PARAMETERS_H