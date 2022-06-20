//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace em3
{
    mem::memory_pool<DendroScalar> EM3_MEM_POOL = mem::memory_pool<DendroScalar>(0,16);

    unsigned int EM3_ELE_ORDER =6;
    unsigned int EM3_PADDING_WIDTH=EM3_ELE_ORDER>>1u;
    unsigned int EM3_IO_OUTPUT_FREQ=10;
    unsigned int EM3_TIME_STEP_OUTPUT_FREQ=10;
    unsigned int EM3_REMESH_TEST_FREQ=10;
    double EM3_IO_OUTPUT_GAP=1.0;

    double EM3_WAVELET_TOL=0.0001;

    double EM3_LOAD_IMB_TOL=0.1;
    unsigned int EM3_SPLIT_FIX=256;
    unsigned int EM3_ASYNC_COMM_K=1;
    double EM3_RK45_TIME_BEGIN=0;
    double EM3_RK45_TIME_END=10;

    double EM3_RK45_DESIRED_TOL=1e-6;


    unsigned int EM3_CHECKPT_FREQ=10;
    unsigned int EM3_RESTORE_SOLVER=0;
    unsigned int EM3_ENABLE_BLOCK_ADAPTIVITY=0;

    std::string EM3_VTU_FILE_PREFIX="em3";
    std::string EM3_CHKPT_FILE_PREFIX="em3_cp";
    std::string EM3_PROFILE_FILE_PREFIX="em3_prof";


    unsigned int EM3_DIM=3;
    unsigned int EM3_MAXDEPTH=8;

    unsigned int EM3_ID_TYPE=0;
    double EM3_CFL_FACTOR=0.1;

    double EM3_GRID_MIN_X=-50.0;
    double EM3_GRID_MAX_X=50.0;
    double EM3_GRID_MIN_Y=-50.0;
    double EM3_GRID_MAX_Y=50.0;
    double EM3_GRID_MIN_Z=-50.0;
    double EM3_GRID_MAX_Z=50.0;

    double EM3_BLK_MIN_X=-6.0;
    double EM3_BLK_MIN_Y=-6.0;
    double EM3_BLK_MIN_Z=-6.0;

    double EM3_BLK_MAX_X=6.0;
    double EM3_BLK_MAX_Y=6.0;
    double EM3_BLK_MAX_Z=6.0;

    double EM3_COMPD_MIN[3]={EM3_GRID_MIN_X,EM3_GRID_MIN_Y,EM3_GRID_MIN_Z};
    double EM3_COMPD_MAX[3]={EM3_GRID_MAX_X,EM3_GRID_MAX_Y,EM3_GRID_MAX_Z};

    double EM3_OCTREE_MIN[3]={0.0,0.0,0.0};
    double EM3_OCTREE_MAX[3]={(double)(1u<<EM3_MAXDEPTH),(double)(1u<<EM3_MAXDEPTH),(double)(1u<<EM3_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double EM3_RK45_TIME_STEP_SIZE=EM3_CFL_FACTOR*(EM3_COMPD_MAX[0]-EM3_COMPD_MIN[0])*(1.0/(double)(1u<<EM3_MAXDEPTH));

    double KO_DISS_SIGMA=0.01;
    unsigned int DISSIPATION_TYPE=0; 

    unsigned int EM3_DENDRO_GRAIN_SZ=1000;

    double EM3_DENDRO_AMR_FAC=0.1;

    unsigned int EM3_NUM_REFINE_VARS=6;
    unsigned int EM3_REFINE_VARIABLE_INDICES[EM3_NUM_VARS]={0,1,2,3,4,5};

    unsigned int EM3_NUM_EVOL_VARS_VTU_OUTPUT=6;
    unsigned int EM3_NUM_CONST_VARS_VTU_OUTPUT=2;
    unsigned int EM3_VTU_OUTPUT_EVOL_INDICES[EM3_NUM_VARS]={0,1,2,3,4,5};
    unsigned int EM3_VTU_OUTPUT_CONST_INDICES[EM3_CONSTRAINT_NUM_VARS]={0,1};

    double EM3_ID_AMP1 = 0.5;
    double EM3_ID_AMP2 = 0.5;
    double EM3_ID_LAMBDA1 = 1.0;
    //double EM3_ID_DELTA1 = 1.0;
    //double EM3_ID_DELTA2 = 1.0;
    //double EM3_ID_XC1 = 0.0;
    //double EM3_ID_YC1 = 0.0;
    //double EM3_ID_ZC1 = 0.0;
    //double EM3_ID_XC2 = 0.0;
    //double EM3_ID_YC2 = 0.0;
    //double EM3_ID_ZC2 = 0.0;
    //double EM3_ID_EPSX1 = 1.0;
    //double EM3_ID_EPSY1 = 1.0;
    //double EM3_ID_EPSZ1 = 1.0;
    //double EM3_ID_EPSX2 = 1.0;
    //double EM3_ID_EPSY2 = 1.0;
    //double EM3_ID_EPSZ2 = 1.0;
    //double EM3_ID_R1 = 0.0;
    //double EM3_ID_R2 = 0.0;
    //double EM3_ID_NU1 = 0.0;
    //double EM3_ID_NU2 = 0.0;
    //double EM3_ID_OMEGA = 0.0;

    double EM3_CHI_REFINE_VAL= std::max(em3::EM3_ID_AMP1,em3::EM3_ID_AMP2);
    double EM3_CHI_COARSEN_VAL=0.1;

    RefinementMode EM3_REFINEMENT_MODE=RefinementMode::WAMR;

    bool EM3_VTU_Z_SLICE_ONLY = true; 

}
