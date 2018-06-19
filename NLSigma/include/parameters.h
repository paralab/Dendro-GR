//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains all the parameters related to NLSM simulation.
*/
//

#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include <string.h>
#include <iostream>




namespace nlsm
{

    /**@brief element order*/
    static const unsigned int NLSM_ELE_ORDER=4;

    /**@brief number of variables*/
    static const unsigned int NLSM_NUM_VARS=2;

    /***@brief number of RK45 stages*/
    static const unsigned int NLSM_RK45_STAGES=6;

    /***@brief number of RK4 stages*/
    static const unsigned int NLSM_RK4_STAGES=4;

    /**@brief number of rk4 stages*/
    static const unsigned int NLSM_RK3_STAGES=3;

    /**@brief CFL stability number number (specifies how dt=NLSM_CFL_FACTOR*dx)*/
    static const double NLSM_CFL_FACTOR=0.1;

    /**@brief: parameter used for adaptive time step update. */
    static const double NLSM_SAFETY_FAC=0.8;

    /**@brief number of internal variables*/
    static const unsigned int NLSM_NUM_VARS_INTENL=(NLSM_RK45_STAGES+1)*NLSM_NUM_VARS;

    /**@brief min bh domain add these to the parameter file.*/
    extern double NLSM_COMPD_MIN[3];
    /**@brief min bh domain @todo add these to the parameter file. */
    extern double NLSM_COMPD_MAX[3];

    /**@brief min coords of the OCTREE */
    extern double NLSM_OCTREE_MIN[3];
    /**@brief max coords of the OCTREE */
    extern double NLSM_OCTREE_MAX[3];

    /**@brief solution output frequency*/
    extern unsigned int NLSM_IO_OUTPUT_FREQ;

    /**@brief timestep norms out put freq.*/
    extern unsigned int NLSM_TIME_STEP_OUTPUT_FREQ;

    /**@brief remesh test frequency*/
    extern unsigned int NLSM_REMESH_TEST_FREQ;

    /**@brief checkpoint store frequency*/
    extern unsigned int NLSM_CHECKPT_FREQ;

    /**@brief restore the solver from check point if set to 1. */
    extern unsigned int NLSM_RESTORE_SOLVER;

    /**@brief use the block adaptivity and disable the AMR*/
    extern unsigned int NLSM_ENABLE_BLOCK_ADAPTIVITY;

    /**@brief file prefix for VTU*/
    extern std::string NLSM_VTU_FILE_PREFIX;

    /**@brief file prefix for write check point*/
    extern std::string NLSM_CHKPT_FILE_PREFIX;

    /**@brief file prefix to write profile info.*/
    extern std::string NLSM_PROFILE_FILE_PREFIX;

    /**@brief number of refine variables*/
    extern unsigned int NLSM_NUM_REFINE_VARS;

    /**@brief indices of refine var ids*/
    extern unsigned int NLSM_REFINE_VARIABLE_INDICES[NLSM_NUM_VARS];

    /**@brief number of evolution variables written to vtu files*/
    extern unsigned int NLSM_NUM_EVOL_VARS_VTU_OUTPUT;

    /**@brief evolution variable IDs written to vtu files*/
    extern unsigned int NLSM_VTU_OUTPUT_EVOL_INDICES[NLSM_NUM_VARS];

    /**@brief solution output gap (instead of freq. we can use to output the solution if currentTime > lastIOOutputTime + NLSM_IO_OUTPUT_GAP)*/
    extern  double NLSM_IO_OUTPUT_GAP;

    /**@brief prefered grain sz to use when selecting active npes*/
    extern unsigned int NLSM_DENDRO_GRAIN_SZ;

    /**@brief AMR coarsening factor (we coarsen if tol<NLSM_DENDRO_AMR_FAC*NLSM_WAVELET_TOL)*/
    extern double NLSM_DENDRO_AMR_FAC;

    /**@brief wavelet tolerance value. */
    extern  double NLSM_WAVELET_TOL;
    /**@brief load-imbalance tolerance value. */
    extern  double NLSM_LOAD_IMB_TOL;
    /**@brief: Splitter fix value*/
    extern unsigned int NLSM_SPLIT_FIX;

    /**@brief: async. communication at a time. (upper bound shoud be NLSM_NUM_VARS) */
    extern unsigned int NLSM_ASYNC_COMM_K;


    /**@brief simulation begin time. */
    extern double NLSM_RK45_TIME_BEGIN;
    /**@brief simulation end time*/
    extern double NLSM_RK45_TIME_END;
    /**@brief rk time step size. */
    extern double NLSM_RK45_TIME_STEP_SIZE;

    /** desired tolerance value for the rk45 method (adaptive time stepping. )*/
    extern double NLSM_RK45_DESIRED_TOL;

    /**@brief BBH initial data type */
    extern unsigned int NLSM_ID_TYPE;

    /**@brief physical coordinates for grid, x_min */
    extern double NLSM_GRID_MIN_X;

    /**@brief physical coordinates for grid, x_max */
    extern double NLSM_GRID_MAX_X;

    /**@brief physical coordinates for grid, y_min */
    extern double NLSM_GRID_MIN_Y;

    /**@brief physical coordinates for grid, y_max */
    extern double NLSM_GRID_MAX_Y;

    /**@brief physical coordinates for grid, z_min */
    extern double NLSM_GRID_MIN_Z;

    /**@brief physical coordinates for grid, z_max */
    extern double NLSM_GRID_MAX_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double NLSM_BLK_MIN_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double NLSM_BLK_MIN_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double NLSM_BLK_MIN_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double NLSM_BLK_MAX_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double NLSM_BLK_MAX_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double NLSM_BLK_MAX_Z;

    /**@brief: dimension of the grid*/
    extern unsigned int NLSM_DIM;

    /**@brief: max refinement level*/
    extern unsigned int NLSM_MAXDEPTH;

    /**@brief: Kreiss-Oliger dissipation */
    extern double KO_DISS_SIGMA;

    /**@brief: Kreiss-Oliger dissipation */
    extern double KO_DISS_SIGMA;



    /**@brief: Initial data Gaussian amplitude */
    extern double NLSM_ID_AMP1;

    /**@brief: Initial data Gaussian amplitude */
    extern double NLSM_ID_AMP2;

    /**@brief: Initial data Gaussian width */
    extern double NLSM_ID_DELTA1;

    /**@brief: Initial data Gaussian width */
    extern double NLSM_ID_DELTA2;

    /**@brief: Initial data Gaussian x offset */
    extern double NLSM_ID_XC1;
    extern double NLSM_ID_YC1;
    extern double NLSM_ID_ZC1;

    /**@brief: Initial data Gaussian x offset */
    extern double NLSM_ID_XC2;
    extern double NLSM_ID_YC2;
    extern double NLSM_ID_ZC2;

    /**@brief: Initial data Gaussian elliptic x factor */
    extern double NLSM_ID_EPSX1;

    /**@brief: Initial data Gaussian elliptic y factor */
    extern double NLSM_ID_EPSY1;

    /**@brief: Initial data Gaussian elliptic x factor */
    extern double NLSM_ID_EPSX2;

    /**@brief: Initial data Gaussian elliptic y factor */
    extern double NLSM_ID_EPSY2;

    /**@brief: Initial data Gaussian R */
    extern double NLSM_ID_R1;

    /**@brief: Initial data Gaussian R */
    extern double NLSM_ID_R2;

    /**@brief: Initial data Gaussian nu */
    extern double NLSM_ID_NU1;

    /**@brief: Initial data Gaussian nu */
    extern double NLSM_ID_NU2;

    /**@brief: Initial data Gaussian Omega */
    extern double NLSM_ID_OMEGA;

}

#endif //SFCSORTBENCH_PARAMETERS_H
