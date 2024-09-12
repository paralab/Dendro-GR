#include "parameters.h"

namespace maxwell {

unsigned int MAXWELL_RESTORE_SOLVER                            = 0;

unsigned int MAXWELL_IO_OUTPUT_FREQ                            = 10;

unsigned int MAXWELL_REMESH_TEST_FREQ                          = 10;

unsigned int MAXWELL_CHECKPT_FREQ                              = 10000;

double MAXWELL_IO_OUTPUT_GAP                                   = 1;

std::string MAXWELL_VTU_FILE_PREFIX                            = "maxwell_gr";

std::string MAXWELL_CHKPT_FILE_PREFIX                          = "maxwell_cp";

std::string MAXWELL_PROFILE_FILE_PREFIX                        = "maxwell_prof";

unsigned int MAXWELL_NUM_EVOL_VARS_VTU_OUTPUT                  = 2;

unsigned int MAXWELL_VTU_OUTPUT_EVOL_INDICES[MAXWELL_NUM_VARS] = {0, 1};

unsigned int MAXWELL_DENDRO_GRAIN_SZ                           = 100;

unsigned int MAXWELL_ASYNC_COMM_K                              = 2;

double MAXWELL_DENDRO_AMR_FAC                                  = 1.0;

double MAXWELL_LOAD_IMB_TOL                                    = 0.1;

unsigned int MAXWELL_DIM                                       = 3;

unsigned int MAXWELL_MAXDEPTH                                  = 8;

double MAXWELL_WAVELET_TOL                                     = 1e-4;

unsigned int MAXWELL_NUM_REFINE_VARS                           = 2;

unsigned int MAXWELL_REFINE_VARIABLE_INDICES[MAXWELL_NUM_VARS] = {0, 1};

double MAXWELL_RK45_TIME_BEGIN                                 = 0;

double MAXWELL_RK45_TIME_END                                   = 30;

double MAXWELL_RK45_TIME_STEP_SIZE                             = 0.01;

double MAXWELL_RK45_DESIRED_TOL                                = 1e-3;

double KO_DISS_SIGMA                                           = 1e-1;

unsigned int MAXWELL_ENABLE_BLOCK_ADAPTIVITY                   = 0;

double MAXWELL_BLK_MIN_X                                       = -6.0;

double MAXWELL_BLK_MIN_Y                                       = -6.0;

double MAXWELL_BLK_MIN_Z                                       = -6.0;

double MAXWELL_BLK_MAX_X                                       = 6.0;

double MAXWELL_BLK_MAX_Y                                       = 6.0;

double MAXWELL_BLK_MAX_Z                                       = 6.0;

double MAXWELL_GRID_MIN_X                                      = -200.0;

double MAXWELL_GRID_MAX_X                                      = 200.0;

double MAXWELL_GRID_MIN_Y                                      = -200.0;

double MAXWELL_GRID_MAX_Y                                      = 200.0;

double MAXWELL_GRID_MIN_Z                                      = -200.0;

double MAXWELL_GRID_MAX_Z                                      = 200.0;

unsigned int MAXWELL_ID_TYPE                                   = 1;

double MAXWELL_COMPD_MIN[3]  = {MAXWELL_GRID_MIN_X, MAXWELL_GRID_MIN_Y,
                                MAXWELL_GRID_MIN_Z};
double MAXWELL_COMPD_MAX[3]  = {MAXWELL_GRID_MAX_X, MAXWELL_GRID_MAX_Y,
                                MAXWELL_GRID_MAX_Z};

double MAXWELL_OCTREE_MIN[3] = {0.0, 0.0, 0.0};
double MAXWELL_OCTREE_MAX[3] = {(double)(1u << MAXWELL_MAXDEPTH),
                                (double)(1u << MAXWELL_MAXDEPTH),
                                (double)(1u << MAXWELL_MAXDEPTH)};

unsigned int MAXWELL_TIME_STEP_OUTPUT_FREQ = 10;

unsigned int MAXWELL_SPLIT_FIX             = 2;

}  // namespace maxwell