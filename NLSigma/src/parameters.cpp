//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace nlsm
{
    
    unsigned int NLSM_ELE_ORDER =4;
    unsigned int NLSM_PADDING_WIDTH = NLSM_ELE_ORDER>>1u;
    unsigned int NLSM_IO_OUTPUT_FREQ=10;
    unsigned int NLSM_TIME_STEP_OUTPUT_FREQ=10;
    unsigned int NLSM_REMESH_TEST_FREQ=10;
    double NLSM_IO_OUTPUT_GAP=1.0;

    double NLSM_WAVELET_TOL=0.0001;

    double NLSM_LOAD_IMB_TOL=0.1;
    unsigned int NLSM_SPLIT_FIX=256;
    unsigned int NLSM_ASYNC_COMM_K=4;
    double NLSM_RK45_TIME_BEGIN=0;
    double NLSM_RK45_TIME_END=10;

    double NLSM_RK45_DESIRED_TOL=1e-6;


    unsigned int NLSM_CHECKPT_FREQ=10;
    unsigned int NLSM_RESTORE_SOLVER=0;
    unsigned int NLSM_ENABLE_BLOCK_ADAPTIVITY=0;

    std::string NLSM_VTU_FILE_PREFIX="nlsm_gr";
    std::string NLSM_CHKPT_FILE_PREFIX="nlsm_cp";
    std::string NLSM_PROFILE_FILE_PREFIX="nlsm_prof";


    unsigned int NLSM_DIM=3;
    unsigned int NLSM_MAXDEPTH=8;
    unsigned int NLSM_MINDEPTH=3;

    unsigned int NLSM_ID_TYPE=0;
    double NLSM_CFL_FACTOR=0.1;

    double NLSM_GRID_MIN_X=-50.0;
    double NLSM_GRID_MAX_X=50.0;
    double NLSM_GRID_MIN_Y=-50.0;
    double NLSM_GRID_MAX_Y=50.0;
    double NLSM_GRID_MIN_Z=-50.0;
    double NLSM_GRID_MAX_Z=50.0;

    double NLSM_BLK_MIN_X=-6.0;
    double NLSM_BLK_MIN_Y=-6.0;
    double NLSM_BLK_MIN_Z=-6.0;

    double NLSM_BLK_MAX_X=6.0;
    double NLSM_BLK_MAX_Y=6.0;
    double NLSM_BLK_MAX_Z=6.0;

    double NLSM_COMPD_MIN[3]={NLSM_GRID_MIN_X,NLSM_GRID_MIN_Y,NLSM_GRID_MIN_Z};
    double NLSM_COMPD_MAX[3]={NLSM_GRID_MAX_X,NLSM_GRID_MAX_Y,NLSM_GRID_MAX_Z};

    double NLSM_OCTREE_MIN[3]={0.0,0.0,0.0};
    double NLSM_OCTREE_MAX[3]={(double)(1u<<NLSM_MAXDEPTH),(double)(1u<<NLSM_MAXDEPTH),(double)(1u<<NLSM_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double NLSM_RK45_TIME_STEP_SIZE=NLSM_CFL_FACTOR*(NLSM_COMPD_MAX[0]-NLSM_COMPD_MIN[0])*(1.0/(double)(1u<<NLSM_MAXDEPTH));

    double KO_DISS_SIGMA=0.01;


    unsigned int NLSM_DENDRO_GRAIN_SZ=1000;

    double NLSM_DENDRO_AMR_FAC=0.1;

    unsigned int NLSM_NUM_REFINE_VARS=2;
    unsigned int NLSM_REFINE_VARIABLE_INDICES[NLSM_NUM_VARS]={0,1};

    unsigned int NLSM_NUM_EVOL_VARS_VTU_OUTPUT=2;
    unsigned int NLSM_VTU_OUTPUT_EVOL_INDICES[NLSM_NUM_VARS]={0,1};

    double NLSM_ID_AMP1 = 0.5;
    double NLSM_ID_AMP2 = 0.5;
    double NLSM_ID_DELTA1 = 1.0;
    double NLSM_ID_DELTA2 = 1.0;
    double NLSM_ID_XC1 = 0.0;
    double NLSM_ID_YC1 = 0.0;
    double NLSM_ID_ZC1 = 0.0;
    double NLSM_ID_XC2 = 0.0;
    double NLSM_ID_YC2 = 0.0;
    double NLSM_ID_ZC2 = 0.0;
    double NLSM_ID_EPSX1 = 1.0;
    double NLSM_ID_EPSY1 = 1.0;
    double NLSM_ID_EPSZ1 = 1.0;
    double NLSM_ID_EPSX2 = 1.0;
    double NLSM_ID_EPSY2 = 1.0;
    double NLSM_ID_EPSZ2 = 1.0;
    double NLSM_ID_R1 = 0.0;
    double NLSM_ID_R2 = 0.0;
    double NLSM_ID_NU1 = 0.0;
    double NLSM_ID_NU2 = 0.0;
    double NLSM_ID_OMEGA = 0.0;

    double NLSM_WAVE_SPEED_X = 1.0;
    double NLSM_WAVE_SPEED_Y = 0.0;
    double NLSM_WAVE_SPEED_Z = 0.0;

    double NLSM_CHI_REFINE_VAL= std::max(nlsm::NLSM_ID_AMP1,nlsm::NLSM_ID_AMP2);
    double NLSM_CHI_COARSEN_VAL=0.1;

    RefineMode NLSM_REFINE_MODE=RefineMode::WAMR;

    
}
