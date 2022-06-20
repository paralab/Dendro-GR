//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace em1
{
    mem::memory_pool<DendroScalar> EM1_MEM_POOL = mem::memory_pool<DendroScalar>(0,16);
    
    unsigned int EM1_ELE_ORDER =6;
    unsigned int EM1_PADDING_WIDTH = EM1_ELE_ORDER>>1u;
    unsigned int EM1_IO_OUTPUT_FREQ=10;
    unsigned int EM1_TIME_STEP_OUTPUT_FREQ=10;
    unsigned int EM1_REMESH_TEST_FREQ=10;
    double EM1_IO_OUTPUT_GAP=1.0;

    double EM1_WAVELET_TOL=0.0001;

    double EM1_LOAD_IMB_TOL=0.1;
    unsigned int EM1_SPLIT_FIX=256;
    unsigned int EM1_ASYNC_COMM_K=1;
    double EM1_RK45_TIME_BEGIN=0;
    double EM1_RK45_TIME_END=10;

    double EM1_RK45_DESIRED_TOL=1e-6;


    unsigned int EM1_CHECKPT_FREQ=10;
    unsigned int EM1_RESTORE_SOLVER=0;
    unsigned int EM1_ENABLE_BLOCK_ADAPTIVITY=0;

    std::string EM1_VTU_FILE_PREFIX="em1";
    std::string EM1_CHKPT_FILE_PREFIX="em1_cp";
    std::string EM1_PROFILE_FILE_PREFIX="em1_prof";


    unsigned int EM1_DIM=3;
    unsigned int EM1_MAXDEPTH=8;

    unsigned int EM1_ID_TYPE=0;
    double EM1_CFL_FACTOR=0.1;

    double EM1_GRID_MIN_X=-50.0;
    double EM1_GRID_MAX_X=50.0;
    double EM1_GRID_MIN_Y=-50.0;
    double EM1_GRID_MAX_Y=50.0;
    double EM1_GRID_MIN_Z=-50.0;
    double EM1_GRID_MAX_Z=50.0;

    double EM1_BLK_MIN_X=-6.0;
    double EM1_BLK_MIN_Y=-6.0;
    double EM1_BLK_MIN_Z=-6.0;

    double EM1_BLK_MAX_X=6.0;
    double EM1_BLK_MAX_Y=6.0;
    double EM1_BLK_MAX_Z=6.0;

    double EM1_COMPD_MIN[3]={EM1_GRID_MIN_X,EM1_GRID_MIN_Y,EM1_GRID_MIN_Z};
    double EM1_COMPD_MAX[3]={EM1_GRID_MAX_X,EM1_GRID_MAX_Y,EM1_GRID_MAX_Z};

    double EM1_OCTREE_MIN[3]={0.0,0.0,0.0};
    double EM1_OCTREE_MAX[3]={(double)(1u<<EM1_MAXDEPTH),(double)(1u<<EM1_MAXDEPTH),(double)(1u<<EM1_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double EM1_RK45_TIME_STEP_SIZE=EM1_CFL_FACTOR*(EM1_COMPD_MAX[0]-EM1_COMPD_MIN[0])*(1.0/(double)(1u<<EM1_MAXDEPTH));

    double KO_DISS_SIGMA=0.01;
    unsigned int DISSIPATION_TYPE=0; 

    unsigned int EM1_DENDRO_GRAIN_SZ=1000;

    double EM1_DENDRO_AMR_FAC=0.1;

    unsigned int EM1_NUM_REFINE_VARS=6;
    unsigned int EM1_REFINE_VARIABLE_INDICES[EM1_NUM_VARS]={0,1,2,3,4,5,6};

    unsigned int EM1_NUM_EVOL_VARS_VTU_OUTPUT=7;
    unsigned int EM1_NUM_CONST_VARS_VTU_OUTPUT=1;
    unsigned int EM1_VTU_OUTPUT_EVOL_INDICES[EM1_NUM_VARS]={0,1,2,3,4,5,6};
    unsigned int EM1_VTU_OUTPUT_CONST_INDICES[EM1_CONSTRAINT_NUM_VARS]={0};

    double EM1_ID_AMP1 = 0.5;
    double EM1_ID_AMP2 = 0.5;
    double EM1_ID_LAMBDA1 = 1.0;
    //double EM1_ID_DELTA1 = 1.0;
    //double EM1_ID_DELTA2 = 1.0;
    //double EM1_ID_XC1 = 0.0;
    //double EM1_ID_YC1 = 0.0;
    //double EM1_ID_ZC1 = 0.0;
    //double EM1_ID_XC2 = 0.0;
    //double EM1_ID_YC2 = 0.0;
    //double EM1_ID_ZC2 = 0.0;
    //double EM1_ID_EPSX1 = 1.0;
    //double EM1_ID_EPSY1 = 1.0;
    //double EM1_ID_EPSZ1 = 1.0;
    //double EM1_ID_EPSX2 = 1.0;
    //double EM1_ID_EPSY2 = 1.0;
    //double EM1_ID_EPSZ2 = 1.0;
    //double EM1_ID_R1 = 0.0;
    //double EM1_ID_R2 = 0.0;
    //double EM1_ID_NU1 = 0.0;
    //double EM1_ID_NU2 = 0.0;
    //double EM1_ID_OMEGA = 0.0;

    double EM1_CHI_REFINE_VAL= std::max(em1::EM1_ID_AMP1,em1::EM1_ID_AMP2);
    double EM1_CHI_COARSEN_VAL=0.1;

    RefinementMode EM1_REFINEMENT_MODE=RefinementMode::WAMR;

    bool EM1_VTU_Z_SLICE_ONLY = true; 

}
