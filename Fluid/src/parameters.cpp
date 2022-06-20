#include "parameters.h"

namespace fluid
{

	unsigned int FLUID_RESTORE_SOLVER = 0;

	unsigned int FLUID_IO_OUTPUT_FREQ = 10;

	unsigned int FLUID_REMESH_TEST_FREQ = 10;

	unsigned int FLUID_CHECKPT_FREQ = 10000;

	double FLUID_IO_OUTPUT_GAP = 1;

	std::string FLUID_VTU_FILE_PREFIX = "fluid_gr";

	std::string FLUID_CHKPT_FILE_PREFIX = "fluid_cp";

	std::string FLUID_PROFILE_FILE_PREFIX = "fluid_prof";

	unsigned int FLUID_NUM_EVOL_VARS_VTU_OUTPUT = 5;

	unsigned int FLUID_VTU_OUTPUT_EVOL_INDICES[FLUID_NUM_EVOL_VARS] = {0,1,2,3,4};
    
    unsigned int FLUID_VTU_OUTPUT_PRIM_INDICES[FLUID_NUM_PRIM_VARS]={0,1,2,3,4,5,6,7};
    
    //unsigned int FLUID_VTU_OUTPUT_CONST_INDICES[FLUID_NUM_CONS_VARS]={0};

	unsigned int FLUID_DENDRO_GRAIN_SZ = 100;

	unsigned int FLUID_ASYNC_COMM_K = 1;

	double FLUID_DENDRO_AMR_FAC = 1;

	double FLUID_LOAD_IMB_TOL = 0.1;

	unsigned int FLUID_DIM = 3;

	unsigned int FLUID_MAXDEPTH = 8;

	double FLUID_WAVELET_TOL = 1.00E-04;

	unsigned int FLUID_NUM_REFINE_VARS = 2;

	unsigned int FLUID_REFINE_VARIABLE_INDICES[FLUID_NUM_PRIM_VARS] = {0,4,5,6,7};

	double FLUID_RK_TIME_BEGIN = 0;

	double FLUID_RK_TIME_END = 30;

	double FLUID_RK45_TIME_STEP_SIZE = 0.01;

	double FLUID_RK45_DESIRED_TOL = 1.00E-03;

	double KO_DISS_SIGMA = 1.00E-01;

	unsigned int FLUID_ENABLE_BLOCK_ADAPTIVITY = 0;

	double FLUID_BLK_MIN_X = -6;

	double FLUID_BLK_MIN_Y = -6;

	double FLUID_BLK_MIN_Z = -6;

	double FLUID_BLK_MAX_X = 6;

	double FLUID_BLK_MAX_Y = 6;

	double FLUID_BLK_MAX_Z = 6;

	double FLUID_GRID_MIN_X = -200;

	double FLUID_GRID_MAX_X = 200;

	double FLUID_GRID_MIN_Y = -200;

	double FLUID_GRID_MAX_Y = 200;

	double FLUID_GRID_MIN_Z = -200;

	double FLUID_GRID_MAX_Z = 200;

	unsigned int FLUID_ID_TYPE = 1;

	unsigned int FLUID_COORDS = 0;

	unsigned int FLUID_RECON_METHOD = 2;

	double FLUID_GAMMA = 1.6667;

	double FLUID_VACUUM = 1.00E-10;

	double FLUID_VACUUM_D = 1.00E-10;

	double FLUID_VACUUM_TAU = 1.00E-10;

  unsigned int FLUID_CYL_GAUSSIAN_AXIS = 2;

	unsigned int FLUID_EOS_SOLVER = 1;

	unsigned int FLUID_USE_WAVE_SPEEDS = 1;

	unsigned int FLUID_CONTOPRIMWARN = 1;

	unsigned int FLUID_AVERAGING_FREQUENCY = 1;

	double FLUID_UNITS_CLU = 1;

	double FLUID_UNITS_CMU = 1;

	double FLUID_UNITS_CKB = 1;

	double FLUID_COMPD_MIN[3]={FLUID_GRID_MIN_X,FLUID_GRID_MIN_Y,FLUID_GRID_MIN_Z};
	double FLUID_COMPD_MAX[3]={FLUID_GRID_MAX_X,FLUID_GRID_MAX_Y,FLUID_GRID_MAX_Z};

	double FLUID_OCTREE_MIN[3]={0.0,0.0,0.0};
	double FLUID_OCTREE_MAX[3]={(double)(1u<<FLUID_MAXDEPTH),(double)(1u<<FLUID_MAXDEPTH),(double)(1u<<FLUID_MAXDEPTH)};

	unsigned int FLUID_TIME_STEP_OUTPUT_FREQ=10;

	unsigned int FLUID_SPLIT_FIX=2;

	double FLUID_VACUUM_RESET=FLUID_VACUUM;
	double FLUID_VACUUM_D_RESET=FLUID_VACUUM_D;
	double FLUID_VACUUM_TAU_RESET=FLUID_VACUUM_TAU;

  bool FLUID_VTU_Z_SLICE_ONLY = false;

  double FLUID_DUFF_MAC_RHO0 = 0.1;
  double FLUID_DUFF_MAC_K0 = 4.0;
  double FLUID_DUFF_MAC_K = 2.0;
  double FLUID_DUFF_MAC_RMIN = 0.001;
  double FLUID_DUFF_MAC_R0 = 0.1;
  double FLUID_DUFF_MAC_REXP = 0.003;
  double FLUID_DUFF_MAC_E0 = 6;
}
