//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace em2
{
    
    mem::memory_pool<DendroScalar> EM2_MEM_POOL = mem::memory_pool<DendroScalar>(0,16);

    unsigned int EM2_ELE_ORDER =6;
    unsigned int EM2_PADDING_WIDTH=EM2_ELE_ORDER>>1u;
    unsigned int EM2_IO_OUTPUT_FREQ=1;
    unsigned int EM2_TIME_STEP_OUTPUT_FREQ=1;
    unsigned int EM2_REMESH_TEST_FREQ=10;
    double EM2_IO_OUTPUT_GAP=1.0;

    double EM2_WAVELET_TOL=0.0001;

    double EM2_LOAD_IMB_TOL=0.1;
    unsigned int EM2_SPLIT_FIX=256;
    unsigned int EM2_ASYNC_COMM_K=1;
    double EM2_RK45_TIME_BEGIN=0;
    double EM2_RK45_TIME_END=10;

    double EM2_RK45_DESIRED_TOL=1e-6;


    unsigned int EM2_CHECKPT_FREQ=10;
    unsigned int EM2_RESTORE_SOLVER=0;
    unsigned int EM2_ENABLE_BLOCK_ADAPTIVITY=0;

    std::string EM2_VTU_FILE_PREFIX="em2";
    std::string EM2_CHKPT_FILE_PREFIX="em2_cp";
    std::string EM2_PROFILE_FILE_PREFIX="em2_prof";


    unsigned int EM2_DIM=3;
    unsigned int EM2_MAXDEPTH=8;

    unsigned int EM2_ID_TYPE=0;
    double EM2_CFL_FACTOR=0.1;

    double EM2_GRID_MIN_X=-50.0;
    double EM2_GRID_MAX_X=50.0;
    double EM2_GRID_MIN_Y=-50.0;
    double EM2_GRID_MAX_Y=50.0;
    double EM2_GRID_MIN_Z=-50.0;
    double EM2_GRID_MAX_Z=50.0;

    double EM2_BLK_MIN_X=-6.0;
    double EM2_BLK_MIN_Y=-6.0;
    double EM2_BLK_MIN_Z=-6.0;

    double EM2_BLK_MAX_X=6.0;
    double EM2_BLK_MAX_Y=6.0;
    double EM2_BLK_MAX_Z=6.0;

    double EM2_COMPD_MIN[3]={EM2_GRID_MIN_X,EM2_GRID_MIN_Y,EM2_GRID_MIN_Z};
    double EM2_COMPD_MAX[3]={EM2_GRID_MAX_X,EM2_GRID_MAX_Y,EM2_GRID_MAX_Z};

    double EM2_OCTREE_MIN[3]={0.0,0.0,0.0};
    double EM2_OCTREE_MAX[3]={(double)(1u<<EM2_MAXDEPTH),(double)(1u<<EM2_MAXDEPTH),(double)(1u<<EM2_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double EM2_RK45_TIME_STEP_SIZE=EM2_CFL_FACTOR*(EM2_COMPD_MAX[0]-EM2_COMPD_MIN[0])*(1.0/(double)(1u<<EM2_MAXDEPTH));

    double KO_DISS_SIGMA=0.01;
    unsigned int DISSIPATION_TYPE=0; 

    unsigned int EM2_DENDRO_GRAIN_SZ=1000;

    double EM2_DENDRO_AMR_FAC=0.1;

    unsigned int EM2_NUM_REFINE_VARS=8;
    unsigned int EM2_REFINE_VARIABLE_INDICES[EM2_NUM_VARS]={0,1,2,3,4,5,6,7};

    unsigned int EM2_NUM_EVOL_VARS_VTU_OUTPUT=8;
    unsigned int EM2_NUM_CONST_VARS_VTU_OUTPUT=8;
    unsigned int EM2_VTU_OUTPUT_EVOL_INDICES[EM2_NUM_VARS]={0,1,2,3,4,5,6,7};
    unsigned int EM2_VTU_OUTPUT_CONST_INDICES[EM2_CONSTRAINT_NUM_VARS]={0,1};

    double EM2_ID_AMP1 = 0.5;
    double EM2_ID_AMP2 = 0.5;
    double EM2_ID_LAMBDA1 = 1.0;
    //double EM2_ID_DELTA1 = 1.0;
    //double EM2_ID_DELTA2 = 1.0;
    //double EM2_ID_XC1 = 0.0;
    //double EM2_ID_YC1 = 0.0;
    //double EM2_ID_ZC1 = 0.0;
    //double EM2_ID_XC2 = 0.0;
    //double EM2_ID_YC2 = 0.0;
    //double EM2_ID_ZC2 = 0.0;
    //double EM2_ID_EPSX1 = 1.0;
    //double EM2_ID_EPSY1 = 1.0;
    //double EM2_ID_EPSZ1 = 1.0;
    //double EM2_ID_EPSX2 = 1.0;
    //double EM2_ID_EPSY2 = 1.0;
    //double EM2_ID_EPSZ2 = 1.0;
    //double EM2_ID_R1 = 0.0;
    //double EM2_ID_R2 = 0.0;
    //double EM2_ID_NU1 = 0.0;
    //double EM2_ID_NU2 = 0.0;
    //double EM2_ID_OMEGA = 0.0;

    double EM2_CHI_REFINE_VAL= std::max(em2::EM2_ID_AMP1,em2::EM2_ID_AMP2);
    double EM2_CHI_COARSEN_VAL=0.1;

    RefinementMode EM2_REFINEMENT_MODE=RefinementMode::WAMR;

    bool EM2_VTU_Z_SLICE_ONLY = true; 

}
