#ifndef DENDRO_5_FLUID_PARAMETERS_H
#define DENDRO_5_FLUID_PARAMETERS_H

#include <string.h>

#include <iostream>

namespace fluid {
/**@brief:element order*/
static const unsigned int FLUID_ELE_ORDER   = 6;

/**@brief:number of stages for rk45*/
static const unsigned int FLUID_RK45_STAGES = 6;

/**@brief: number of stages for RK4*/
static const unsigned int FLUID_RK4_STAGES  = 4;

/**@brief: number of stages for RK3*/
static const unsigned int FLUID_RK3_STAGES  = 3;

/**@brief: parameter for adaptive time step update (rk45)*/
static const double FLUID_SAFETY_FAC        = 0.8;

/**@brief timestep norms out put freq.*/
extern unsigned int FLUID_TIME_STEP_OUTPUT_FREQ;

/**@brief: dendro splitter fix value */
extern unsigned int FLUID_SPLIT_FIX;

/**@brief min of computational domain*/
extern double FLUID_COMPD_MIN[3];

/**@brief max of computational domain*/
extern double FLUID_COMPD_MAX[3];

/**@brief: min coords of the octree*/
extern double FLUID_OCTREE_MIN[3];

/**@brief: max coords of the octree*/
extern double FLUID_OCTREE_MAX[3];

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_VACUUM_RESET;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_VACUUM_D_RESET;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_VACUUM_TAU_RESET;

/**@brief: number of primitive fluid variables*/
static const unsigned int FLUID_NUM_PRIM_VARS = 8;

/**@brief: number of evolv variables*/
static const unsigned int FLUID_NUM_EVOL_VARS = 5;

/**@brief: number of constraint varibales */
// static const unsigned int FLUID_NUM_CONS_VARS=1;

/**@brief: 1 indicates try to restore from the closest check point. */
extern unsigned int FLUID_RESTORE_SOLVER;

/**@brief: VTU output frequncy*/
extern unsigned int FLUID_IO_OUTPUT_FREQ;

/**@brief: remesh test frequency */
extern unsigned int FLUID_REMESH_TEST_FREQ;

/**@brief: checkpoint write frequency */
extern unsigned int FLUID_CHECKPT_FREQ;

/**@brief: vtu file write based on the rk time (frequency) */
extern double FLUID_IO_OUTPUT_GAP;

/**@brief: vtu file prefix*/
extern std::string FLUID_VTU_FILE_PREFIX;

/**@brief: checkpoint file prefix*/
extern std::string FLUID_CHKPT_FILE_PREFIX;

/**@brief: profile file prefix*/
extern std::string FLUID_PROFILE_FILE_PREFIX;

/**@brief: number of evolve vars to output*/
extern unsigned int FLUID_NUM_EVOL_VARS_VTU_OUTPUT;

/**@brief: evolve var ids to vtu output*/
extern unsigned int FLUID_VTU_OUTPUT_EVOL_INDICES[FLUID_NUM_EVOL_VARS];

/**@brief: primitive variable ids*/
extern unsigned int FLUID_VTU_OUTPUT_PRIM_INDICES[FLUID_NUM_PRIM_VARS];

/**@brief: constraint variable ids*/
// extern unsigned int FLUID_VTU_OUTPUT_CONST_INDICES[FLUID_NUM_CONS_VARS];

/**@brief: dendro grain sz (number of octants per core)*/
extern unsigned int FLUID_DENDRO_GRAIN_SZ;

/**@brief: number of variables for asynchronous communication*/
extern unsigned int FLUID_ASYNC_COMM_K;

/**@brief: AMR factor used for octree coarsening <1*/
extern double FLUID_DENDRO_AMR_FAC;

/**@brief:load imblance tolerance (for flexible partitioning)*/
extern double FLUID_LOAD_IMB_TOL;

/**@brief:dimension of the problem*/
extern unsigned int FLUID_DIM;

/**@brief: maximum depth (refinement levels for the octree)*/
extern unsigned int FLUID_MAXDEPTH;

/**@brief: Wavelet tolerance for WAMR*/
extern double FLUID_WAVELET_TOL;

/**@brief: number of variables to consider for the refinemenet. */
extern unsigned int FLUID_NUM_REFINE_VARS;

/**@brief: variable IDs considered for the refinement. */
extern unsigned int FLUID_REFINE_VARIABLE_INDICES[FLUID_NUM_PRIM_VARS];

/**@brief:rk time begin. */
extern double FLUID_RK_TIME_BEGIN;

/**@brief: rk time end*/
extern double FLUID_RK_TIME_END;

/**@brief: initial time step size fo the rk45 method. */
extern double FLUID_RK45_TIME_STEP_SIZE;

/**@brief: tolerance value for adaptive time step size. */
extern double FLUID_RK45_DESIRED_TOL;

/**@brief: ko dissipation value*/
extern double KO_DISS_SIGMA;

/**@brief: if 1 enable the block adaptivity*/
extern unsigned int FLUID_ENABLE_BLOCK_ADAPTIVITY;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double FLUID_BLK_MIN_X;

/**@brief physical coordinates for the blk adaptive y_min*/
extern double FLUID_BLK_MIN_Y;

/**@brief physical coordinates for the blk adaptive z_min*/
extern double FLUID_BLK_MIN_Z;

/**@brief physical coordinates for the blk adaptive x_max*/
extern double FLUID_BLK_MAX_X;

/**@brief physical coordinates for the blk adaptive y_max*/
extern double FLUID_BLK_MAX_Y;

/**@brief physical coordinates for the blk adaptive z_max*/
extern double FLUID_BLK_MAX_Z;

/**@brief physical coordinates for grid, x_min (blk adptivity) */
extern double FLUID_GRID_MIN_X;

/**@brief physical coordinates for grid, x_max (blk adptivity) */
extern double FLUID_GRID_MAX_X;

/**@brief physical coordinates for grid, y_min (blk adptivity) */
extern double FLUID_GRID_MIN_Y;

/**@brief physical coordinates for grid, y_max (blk adptivity) */
extern double FLUID_GRID_MAX_Y;

/**@brief physical coordinates for grid, z_min (blk adptivity) */
extern double FLUID_GRID_MIN_Z;

/**@brief physical coordinates for grid, z_max (blk adptivity) */
extern double FLUID_GRID_MAX_Z;

/**@brief:CFL factor for time stepping*/
static const double FLUID_CFL_FACTOR = 0.1;

/**@brief: fluid initial condition type*/
extern unsigned int FLUID_ID_TYPE;

/**@brief: @jacob can you please add a discription for the parameter*/
extern unsigned int FLUID_COORDS;

/**@brief: @jacob can you please add a discription for the parameter*/
extern unsigned int FLUID_RECON_METHOD;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_GAMMA;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_VACUUM;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_VACUUM_D;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_VACUUM_TAU;

/**@brief: which axis the cylindrical Gaussian is aligned to.*/
extern unsigned int FLUID_CYL_GAUSSIAN_AXIS;

/**@brief: @jacob can you please add a discription for the parameter*/
extern unsigned int FLUID_EOS_SOLVER;

/**@brief: @jacob can you please add a discription for the parameter*/
extern unsigned int FLUID_USE_WAVE_SPEEDS;

/**@brief: @jacob can you please add a discription for the parameter*/
extern unsigned int FLUID_CONTOPRIMWARN;

/**@brief: @jacob can you please add a discription for the parameter*/
extern unsigned int FLUID_AVERAGING_FREQUENCY;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_UNITS_CLU;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_UNITS_CMU;

/**@brief: @jacob can you please add a discription for the parameter*/
extern double FLUID_UNITS_CKB;

/**@brief: Write to 2d VTU slice instead of 3d VTU block*/
extern bool FLUID_VTU_Z_SLICE_ONLY;

// Duffell-MacFadyen initial data parameters
/**@brief: Reference density at r0.*/
extern double FLUID_DUFF_MAC_RHO0;

/**@brief: Density profile inside point mass region.*/
extern double FLUID_DUFF_MAC_K0;

/**@brief: Density profile of ambient medium.*/
extern double FLUID_DUFF_MAC_K;

/**@brief: Radius determining constant profile.*/
extern double FLUID_DUFF_MAC_RMIN;

/**@brief: Radius determining end of point mass region.*/
extern double FLUID_DUFF_MAC_R0;

/**@brief: Radius determining explosion region, i.e., region of high pressure.*/
extern double FLUID_DUFF_MAC_REXP;

/**@brief: Radius determining energy/ratio between pressure and density.*/
extern double FLUID_DUFF_MAC_E0;

/**@brief:*/
static const unsigned int FLUID_NUM_VARS_INTENL =
    (FLUID_RK45_STAGES + 1) * FLUID_NUM_PRIM_VARS;

}  // namespace fluid

#endif  // DENDRO_5_FLUID_PARAMETERS_H
