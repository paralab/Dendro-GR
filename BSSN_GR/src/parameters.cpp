//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace bssn
{

    mem::memory_pool<double> BSSN_MEM_POOL = mem::memory_pool<double>(0,16);
    unsigned int BSSN_ELE_ORDER =6;
    unsigned int BSSN_PADDING_WIDTH=BSSN_ELE_ORDER>>1u;
    
    unsigned int BSSN_IO_OUTPUT_FREQ=10;
    unsigned int BSSN_TIME_STEP_OUTPUT_FREQ=10;

    unsigned int BSSN_GW_EXTRACT_FREQ=std::max(1u,BSSN_IO_OUTPUT_FREQ>>1u);
    
    unsigned int BSSN_REMESH_TEST_FREQ=10;
    unsigned int BSSN_REMESH_TEST_FREQ_AFTER_MERGER=10;
    unsigned int BSSN_GW_EXTRACT_FREQ_AFTER_MERGER=10;
    double BSSN_IO_OUTPUT_GAP=1.0;

    unsigned int BSSN_USE_WAVELET_TOL_FUNCTION=0;
    double BSSN_WAVELET_TOL=0.0001;
    double BSSN_GW_REFINE_WTOL=1e-4;
    double BSSN_WAVELET_TOL_MAX=0.001;
    double BSSN_WAVELET_TOL_FUNCTION_R0=10.0;
    double BSSN_WAVELET_TOL_FUNCTION_R1=50.0;

    double BSSN_CFL_FACTOR=0.1;

    double BSSN_LOAD_IMB_TOL=0.1;
    unsigned int BSSN_SPLIT_FIX=2;
    unsigned int BSSN_ASYNC_COMM_K=4;
    double BSSN_RK_TIME_BEGIN=0;
    double BSSN_RK_TIME_END=10;

    double BSSN_RK45_DESIRED_TOL=1e-6;

    unsigned int BSSN_RK_TYPE;


    unsigned int BSSN_CHECKPT_FREQ=10;
    unsigned int BSSN_RESTORE_SOLVER=0;
    unsigned int BSSN_ENABLE_BLOCK_ADAPTIVITY=0;

    double BSSN_ETA_R0=1.31;
    double BSSN_ETA_POWER[]={2.0,2.0};

    std::string BSSN_VTU_FILE_PREFIX="bssn_gr";
    std::string BSSN_CHKPT_FILE_PREFIX="bssn_cp";
    std::string BSSN_PROFILE_FILE_PREFIX="bssn_prof";

    double BSSN_BH1_AMR_R=2.0;
    double BSSN_BH2_AMR_R=2.0;

    double BSSN_BH1_CONSTRAINT_R=5.0;
    double BSSN_BH2_CONSTRAINT_R=5.0;

    double BSSN_BH1_MASS;
    double BSSN_BH2_MASS;


    unsigned int BSSN_BH1_MAX_LEV;
    unsigned int BSSN_BH2_MAX_LEV;

    unsigned int BSSN_INIT_GRID_ITER=10;

    BH BH1;
    BH BH2;
    Point BSSN_BH_LOC[2];
    unsigned int BSSN_DIM=3;
    unsigned int BSSN_MAXDEPTH=8;
    unsigned int BSSN_MINDEPTH=3;


    unsigned int BSSN_ID_TYPE=0;

    double BSSN_GRID_MIN_X=-50.0;
    double BSSN_GRID_MAX_X=50.0;
    double BSSN_GRID_MIN_Y=-50.0;
    double BSSN_GRID_MAX_Y=50.0;
    double BSSN_GRID_MIN_Z=-50.0;
    double BSSN_GRID_MAX_Z=50.0;

    double BSSN_BLK_MIN_X=-6.0;
    double BSSN_BLK_MIN_Y=-6.0;
    double BSSN_BLK_MIN_Z=-6.0;

    double BSSN_BLK_MAX_X=6.0;
    double BSSN_BLK_MAX_Y=6.0;
    double BSSN_BLK_MAX_Z=6.0;

    double BSSN_COMPD_MIN[3]={BSSN_GRID_MIN_X,BSSN_GRID_MIN_Y,BSSN_GRID_MIN_Z};
    double BSSN_COMPD_MAX[3]={BSSN_GRID_MAX_X,BSSN_GRID_MAX_Y,BSSN_GRID_MAX_Z};

    double BSSN_OCTREE_MIN[3]={0.0,0.0,0.0};
    double BSSN_OCTREE_MAX[3]={(double)(1u<<BSSN_MAXDEPTH),(double)(1u<<BSSN_MAXDEPTH),(double)(1u<<BSSN_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double BSSN_RK45_TIME_STEP_SIZE=BSSN_CFL_FACTOR*(BSSN_COMPD_MAX[0]-BSSN_COMPD_MIN[0])*(1.0/(double)(1u<<BSSN_MAXDEPTH));

    unsigned int BSSN_LAMBDA[4]={1, 1, 1, 1};
    double BSSN_LAMBDA_F[2]={1.0, 0.0};
    double BSSN_TRK0=0.0;
    double ETA_CONST=2.0;
    double ETA_R0=50.0;
    double ETA_DAMPING=1.0;
    double ETA_DAMPING_EXP=1.0;
    double CHI_FLOOR=0.1;
    double KO_DISS_SIGMA=0.01;

    unsigned int RIT_ETA_FUNCTION = 1;
    double RIT_ETA_OUTER = 0.25;
    double RIT_ETA_CENTRAL = 2.0;
    double RIT_ETA_WIDTH = 40.0;

    unsigned int DISSIPATION_TYPE=0;

    unsigned int BSSN_DENDRO_GRAIN_SZ=1000;

    double BSSN_DENDRO_AMR_FAC=0.1;

    unsigned int BSSN_NUM_REFINE_VARS=1;
    unsigned int BSSN_REFINE_VARIABLE_INDICES[BSSN_NUM_VARS]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

    unsigned int BSSN_NUM_EVOL_VARS_VTU_OUTPUT=1;
    unsigned int BSSN_NUM_CONST_VARS_VTU_OUTPUT=1;
    unsigned int BSSN_VTU_OUTPUT_EVOL_INDICES[BSSN_NUM_VARS]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
    unsigned int BSSN_VTU_OUTPUT_CONST_INDICES[BSSN_CONSTRAINT_NUM_VARS]={0,1,2,3,4,5};

    unsigned int BSSN_XI[3]={0 , 0 , 0};

    unsigned int BSSN_DISSIPATION_NC=0;

    unsigned int BSSN_DISSIPATION_S=10;

    double BSSN_EH_REFINE_VAL = 0.3 ;
    double BSSN_EH_COARSEN_VAL = 0.4;

    // by default use WAMR refinement. 
    RefinementMode BSSN_REFINEMENT_MODE = RefinementMode::WAMR;

    bool BSSN_VTU_Z_SLICE_ONLY = true;

    unsigned int BSSN_LTS_TS_OFFSET=4;

    bool BSSN_MERGED_CHKPT_WRITTEN=false;

    double BSSN_CURRENT_RK_COORD_TIME  = 0;
    unsigned int  BSSN_CURRENT_RK_STEP = 0;

    /***@brief: derivs workspace*/
    double* BSSN_DERIV_WORKSPACE=nullptr;


}

namespace TPID {
  double target_M_plus=1.0;
  double target_M_minus=1.0;
  double par_m_plus=1.0;
  double par_m_minus=1.0;
  double par_b=4.0;
  double par_P_plus[3]={0.0,0.0,0.0};
  double par_P_minus[3]={0.0,0.0,0.0};
  double par_S_plus[3]={0.0,0.0,0.0};
  double par_S_minus[3]={0.0,0.0,0.0};
  double center_offset[3]={0.0,0.0,0.00014142135623730951};
  double initial_lapse_psi_exponent=-2;
  int npoints_A=30;
  int npoints_B=30;
  int npoints_phi=16;
  int give_bare_mass=0;
  int initial_lapse=2;
  int solve_momentum_constraint=1;
  int grid_setup_method=1;
  int verbose = 1;
  double adm_tol=1.0e-10;
  double Newton_tol=1.0e-10;
  std::string FILE_PREFIX="tpid";
}


namespace BHLOC
{
    unsigned int EXTRACTION_VAR_ID=bssn::VAR::U_ALPHA;
    double EXTRACTION_TOL=0.3;
}


namespace GW
{
    
    unsigned int BSSN_GW_NUM_RADAII;
    
    unsigned int BSSN_GW_NUM_LMODES;
    
    double BSSN_GW_RADAII[BSSN_GW_MAX_RADAII];
    
    unsigned int BSSN_GW_SPIN=2;
    
    unsigned int BSSN_GW_L_MODES[BSSN_GW_MAX_LMODES];
}

namespace AEH
{
    unsigned int AEH_LMAX        = 6;
    unsigned int AEH_Q_THETA     = 32;
    unsigned int AEH_Q_PHI       = 32 ;
    unsigned int AEH_MAXITER     = 50;
    double AEH_ATOL              = 1e-8;
    double AEH_RTOL              = 1e-8;
    unsigned int AEH_SOLVER_FREQ = 0;
    
}
