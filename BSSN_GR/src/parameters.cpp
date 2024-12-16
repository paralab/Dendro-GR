//
// Created by milinda on 8/23/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief
 */
//

#include "parameters.h"

namespace bssn {

mem::memory_pool<double> BSSN_MEM_POOL  = mem::memory_pool<double>(0, 16);
unsigned int BSSN_ELE_ORDER             = 6;
unsigned int BSSN_PADDING_WIDTH         = BSSN_ELE_ORDER >> 1u;

unsigned int BSSN_IO_OUTPUT_FREQ        = 10;
unsigned int BSSN_TIME_STEP_OUTPUT_FREQ = 10;

unsigned int BSSN_GW_EXTRACT_FREQ  = std::max(1u, BSSN_IO_OUTPUT_FREQ >> 1u);

unsigned int BSSN_REMESH_TEST_FREQ = 10;

unsigned int BSSN_REMESH_TEST_FREQ_AFTER_MERGER        = 10;
unsigned int BSSN_GW_EXTRACT_FREQ_AFTER_MERGER         = 10;
double BSSN_IO_OUTPUT_GAP                              = 1.0;

unsigned int BSSN_USE_WAVELET_TOL_FUNCTION             = 0;
double BSSN_WAVELET_TOL                                = 0.0001;
double BSSN_GW_REFINE_WTOL                             = 1e-4;
double BSSN_WAVELET_TOL_MAX                            = 0.001;
double BSSN_WAVELET_TOL_FUNCTION_R0                    = 10.0;
double BSSN_WAVELET_TOL_FUNCTION_R1                    = 50.0;

double BSSN_CFL_FACTOR                                 = 0.1;

bool BSSN_KO_SIGMA_SCALE_BY_CONFORMAL                  = false;
bool BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY = false;
double BSSN_EPSILON_CAKO_GAUGE                         = 0.99;
double BSSN_EPSILON_CAKO_OTHER                         = 0.3;
bool BSSN_CAKO_ENABLED                                 = false;
double BSSN_CAHD_C                                     = 0.0;

double BSSN_LOAD_IMB_TOL                               = 0.1;
unsigned int BSSN_SPLIT_FIX                            = 2;
unsigned int BSSN_ASYNC_COMM_K                         = 4;
double BSSN_RK_TIME_BEGIN                              = 0;
double BSSN_RK_TIME_END                                = 10;

double BSSN_RK45_DESIRED_TOL                           = 1e-6;

unsigned int BSSN_RK_TYPE;

unsigned int BSSN_CHECKPT_FREQ            = 10;
unsigned int BSSN_RESTORE_SOLVER          = 0;
unsigned int BSSN_ENABLE_BLOCK_ADAPTIVITY = 0;

double BSSN_ETA_R0                        = 1.31;
double BSSN_ETA_POWER[]                   = {2.0, 2.0};

std::string BSSN_VTU_FILE_PREFIX          = "bssn_gr";
std::string BSSN_CHKPT_FILE_PREFIX        = "bssn_cp";
std::string BSSN_PROFILE_FILE_PREFIX      = "bssn_prof";

double BSSN_BH1_AMR_R                     = 2.0;
double BSSN_BH2_AMR_R                     = 2.0;
double BSSN_AMR_R_RATIO =
    2.5;  // ratio for the near to far portion, was originally 2.5

double BSSN_BH1_CONSTRAINT_R = 5.0;
double BSSN_BH2_CONSTRAINT_R = 5.0;

double BSSN_BH1_MASS;
double BSSN_BH2_MASS;

unsigned int BSSN_BH1_MAX_LEV;
unsigned int BSSN_BH2_MAX_LEV;

unsigned int BSSN_INIT_GRID_ITER = 10;

BH BH1;
BH BH2;
Point BSSN_BH_LOC[2];
unsigned int BSSN_DIM      = 3;
unsigned int BSSN_MAXDEPTH = 8;
unsigned int BSSN_MINDEPTH = 3;

unsigned int BSSN_ID_TYPE  = 0;

double BSSN_GRID_MIN_X     = -50.0;
double BSSN_GRID_MAX_X     = 50.0;
double BSSN_GRID_MIN_Y     = -50.0;
double BSSN_GRID_MAX_Y     = 50.0;
double BSSN_GRID_MIN_Z     = -50.0;
double BSSN_GRID_MAX_Z     = 50.0;

double BSSN_BLK_MIN_X      = -6.0;
double BSSN_BLK_MIN_Y      = -6.0;
double BSSN_BLK_MIN_Z      = -6.0;

double BSSN_BLK_MAX_X      = 6.0;
double BSSN_BLK_MAX_Y      = 6.0;
double BSSN_BLK_MAX_Z      = 6.0;

double BSSN_COMPD_MIN[3]  = {BSSN_GRID_MIN_X, BSSN_GRID_MIN_Y, BSSN_GRID_MIN_Z};
double BSSN_COMPD_MAX[3]  = {BSSN_GRID_MAX_X, BSSN_GRID_MAX_Y, BSSN_GRID_MAX_Z};

double BSSN_OCTREE_MIN[3] = {0.0, 0.0, 0.0};
double BSSN_OCTREE_MAX[3] = {(double)(1u << BSSN_MAXDEPTH),
                             (double)(1u << BSSN_MAXDEPTH),
                             (double)(1u << BSSN_MAXDEPTH)};

//@note assumes the computational domain is a cube as well.
double BSSN_RK45_TIME_STEP_SIZE = BSSN_CFL_FACTOR *
                                  (BSSN_COMPD_MAX[0] - BSSN_COMPD_MIN[0]) *
                                  (1.0 / (double)(1u << BSSN_MAXDEPTH));

// calculate the minimum dx
double BSSN_CURRENT_MIN_DX = (BSSN_COMPD_MAX[0] - BSSN_COMPD_MIN[0]) *
                             (1.0 / (double)(1u << BSSN_MAXDEPTH));

unsigned int BSSN_LAMBDA[4]                              = {1, 1, 1, 1};
double BSSN_LAMBDA_F[2]                                  = {1.0, 0.0};
double BSSN_TRK0                                         = 0.0;
double ETA_CONST                                         = 2.0;
double ETA_R0                                            = 50.0;
double ETA_DAMPING                                       = 1.0;
double ETA_DAMPING_EXP                                   = 1.0;
double CHI_FLOOR                                         = 0.1;
double KO_DISS_SIGMA                                     = 0.01;

unsigned int RIT_ETA_FUNCTION                            = 1;
double RIT_ETA_OUTER                                     = 0.25;
double RIT_ETA_CENTRAL                                   = 2.0;
double RIT_ETA_WIDTH                                     = 40.0;

unsigned int DISSIPATION_TYPE                            = 0;

unsigned int BSSN_DENDRO_GRAIN_SZ                        = 1000;

double BSSN_DENDRO_AMR_FAC                               = 0.1;

unsigned int BSSN_NUM_REFINE_VARS                        = 1;
unsigned int BSSN_REFINE_VARIABLE_INDICES[BSSN_NUM_VARS] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

unsigned int BSSN_NUM_EVOL_VARS_VTU_OUTPUT               = 1;
unsigned int BSSN_NUM_CONST_VARS_VTU_OUTPUT              = 1;
unsigned int BSSN_VTU_OUTPUT_EVOL_INDICES[BSSN_NUM_VARS] = {
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
unsigned int BSSN_VTU_OUTPUT_CONST_INDICES[BSSN_CONSTRAINT_NUM_VARS] = {
    0, 1, 2, 3, 4, 5};

unsigned int BSSN_XI[3]                         = {0, 0, 0};

unsigned int BSSN_DISSIPATION_NC                = 0;

unsigned int BSSN_DISSIPATION_S                 = 10;

double BSSN_EH_REFINE_VAL                       = 0.3;
double BSSN_EH_COARSEN_VAL                      = 0.4;

// by default use WAMR refinement.
RefinementMode BSSN_REFINEMENT_MODE             = RefinementMode::WAMR;

bool BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE = false;

bool BSSN_VTU_Z_SLICE_ONLY                      = true;

unsigned int BSSN_LTS_TS_OFFSET                 = 4;

bool BSSN_MERGED_CHKPT_WRITTEN                  = false;

double BSSN_CURRENT_RK_COORD_TIME               = 0;
unsigned int BSSN_CURRENT_RK_STEP               = 0;

unsigned int BSSN_NYQUIST_M                     = 0;

bool BSSN_SCALE_VTU_AND_GW_EXTRACTION           = false;

unsigned int BSSN_GW_EXTRACT_FREQ_TRUE          = 0;

unsigned int BSSN_IO_OUTPUT_FREQ_TRUE           = 0;

double BSSN_SSL_SIGMA                           = 20.0;
double BSSN_SSL_H                               = 0.6;

/***@brief: derivs workspace*/
double* BSSN_DERIV_WORKSPACE                    = nullptr;

void readParamTOMLFile(const char* fName, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    auto parFile                = toml::parse(fName);

    bssn::BSSN_IO_OUTPUT_FREQ   = parFile["BSSN_IO_OUTPUT_FREQ"].as_integer();
    bssn::BSSN_REMESH_TEST_FREQ = parFile["BSSN_REMESH_TEST_FREQ"].as_integer();
    bssn::BSSN_CHECKPT_FREQ     = parFile["BSSN_CHECKPT_FREQ"].as_integer();
    bssn::BSSN_IO_OUTPUT_GAP    = parFile["BSSN_IO_OUTPUT_GAP"].as_integer();
    bssn::BSSN_VTU_FILE_PREFIX  = parFile["BSSN_VTU_FILE_PREFIX"].as_string();
    bssn::BSSN_CHKPT_FILE_PREFIX =
        parFile["BSSN_CHKPT_FILE_PREFIX"].as_string();

    bssn::BSSN_PROFILE_FILE_PREFIX =
        parFile["BSSN_PROFILE_FILE_PREFIX"].as_string();
    bssn::BSSN_RESTORE_SOLVER = parFile["BSSN_RESTORE_SOLVER"].as_integer();
    bssn::BSSN_ID_TYPE        = parFile["BSSN_ID_TYPE"].as_integer();

    bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY =
        parFile["BSSN_ENABLE_BLOCK_ADAPTIVITY"].as_integer();

    bssn::BSSN_BLK_MIN_X       = parFile["BSSN_BLK_MIN_X"].as_floating();
    bssn::BSSN_BLK_MIN_Y       = parFile["BSSN_BLK_MIN_Y"].as_floating();
    bssn::BSSN_BLK_MIN_Z       = parFile["BSSN_BLK_MIN_Z"].as_floating();
    bssn::BSSN_BLK_MAX_X       = parFile["BSSN_BLK_MAX_X"].as_floating();
    bssn::BSSN_BLK_MAX_Y       = parFile["BSSN_BLK_MAX_Y"].as_floating();
    bssn::BSSN_BLK_MAX_Z       = parFile["BSSN_BLK_MAX_Z"].as_floating();

    bssn::BSSN_DENDRO_GRAIN_SZ = parFile["BSSN_DENDRO_GRAIN_SZ"].as_integer();
    bssn::BSSN_ASYNC_COMM_K    = parFile["BSSN_ASYNC_COMM_K"].as_integer();
    bssn::BSSN_DENDRO_AMR_FAC  = parFile["BSSN_DENDRO_AMR_FAC"].as_floating();
    bssn::BSSN_LOAD_IMB_TOL    = parFile["BSSN_LOAD_IMB_TOL"].as_floating();
    bssn::BSSN_RK_TIME_BEGIN   = parFile["BSSN_RK_TIME_BEGIN"].as_floating();
    bssn::BSSN_RK_TIME_END     = parFile["BSSN_RK_TIME_END"].as_floating();
    bssn::BSSN_RK_TYPE         = parFile["BSSN_RK_TYPE"].as_integer();
    bssn::BSSN_RK45_TIME_STEP_SIZE =
        parFile["BSSN_RK45_TIME_STEP_SIZE"].as_floating();
    bssn::BSSN_RK45_DESIRED_TOL =
        parFile["BSSN_RK45_DESIRED_TOL"].as_floating();
    bssn::BSSN_DIM        = parFile["BSSN_DIM"].as_integer();
    bssn::BSSN_MAXDEPTH   = parFile["BSSN_MAXDEPTH"].as_integer();

    bssn::BH1             = BH(parFile["BSSN_BH1"]["MASS"].as_floating(),
                               parFile["BSSN_BH1"]["X"].as_floating(),
                               parFile["BSSN_BH1"]["Y"].as_floating(),
                               parFile["BSSN_BH1"]["Z"].as_floating(),
                               parFile["BSSN_BH1"]["V_X"].as_floating(),
                               parFile["BSSN_BH1"]["V_Y"].as_floating(),
                               parFile["BSSN_BH1"]["V_Z"].as_floating(),
                               parFile["BSSN_BH1"]["SPIN"].as_floating(),
                               parFile["BSSN_BH1"]["SPIN_THETA"].as_floating(),
                               parFile["BSSN_BH1"]["SPIN_PHI"].as_floating());
    bssn::BH2             = BH(parFile["BSSN_BH2"]["MASS"].as_floating(),
                               parFile["BSSN_BH2"]["X"].as_floating(),
                               parFile["BSSN_BH2"]["Y"].as_floating(),
                               parFile["BSSN_BH2"]["Z"].as_floating(),
                               parFile["BSSN_BH2"]["V_X"].as_floating(),
                               parFile["BSSN_BH2"]["V_Y"].as_floating(),
                               parFile["BSSN_BH2"]["V_Z"].as_floating(),
                               parFile["BSSN_BH2"]["SPIN"].as_floating(),
                               parFile["BSSN_BH2"]["SPIN_THETA"].as_floating(),
                               parFile["BSSN_BH2"]["SPIN_PHI"].as_floating());

    bssn::BSSN_GRID_MIN_X = parFile["BSSN_GRID_MIN_X"].as_floating();
    bssn::BSSN_GRID_MAX_X = parFile["BSSN_GRID_MAX_X"].as_floating();
    bssn::BSSN_GRID_MIN_Y = parFile["BSSN_GRID_MIN_Y"].as_floating();
    bssn::BSSN_GRID_MAX_Y = parFile["BSSN_GRID_MAX_Y"].as_floating();
    bssn::BSSN_GRID_MIN_Z = parFile["BSSN_GRID_MIN_Z"].as_floating();
    bssn::BSSN_GRID_MAX_Z = parFile["BSSN_GRID_MAX_Z"].as_floating();

    bssn::ETA_CONST       = parFile["ETA_CONST"].as_floating();
    bssn::ETA_R0          = parFile["ETA_R0"].as_floating();
    bssn::ETA_DAMPING     = parFile["ETA_DAMPING"].as_floating();
    bssn::ETA_DAMPING_EXP = parFile["ETA_DAMPING_EXP"].as_floating();

    if (parFile.contains("RIT_ETA_FUNCTION")) {
        bssn::RIT_ETA_FUNCTION = parFile["RIT_ETA_FUNCTION"].as_integer();
    }
    if (parFile.contains("RIT_ETA_OUTER")) {
        bssn::RIT_ETA_OUTER = parFile["RIT_ETA_OUTER"].as_floating();
    }
    if (parFile.contains("RIT_ETA_CENTRAL")) {
        bssn::RIT_ETA_CENTRAL = parFile["RIT_ETA_CENTRAL"].as_floating();
    }
    if (parFile.contains("RIT_ETA_WIDTH")) {
        bssn::RIT_ETA_WIDTH = parFile["RIT_ETA_WIDTH"].as_floating();
    }

    if (parFile.contains("BSSN_AMR_R_RATIO")) {
        bssn::BSSN_AMR_R_RATIO = parFile["BSSN_AMR_R_RATIO"].as_floating();
    }

    if (parFile.contains("BSSN_KO_SIGMA_SCALE_BY_CONFORMAL")) {
        bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL =
            parFile["BSSN_KO_SIGMA_SCALE_BY_CONFORMAL"].as_boolean();
    }

    if (parFile.contains("BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY")) {
        bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY =
            parFile["BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY"]
                .as_boolean();
    }

    if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY) {
        bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL = false;
    }
    if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL) {
        bssn::BSSN_CAKO_ENABLED = true;
    }

    if (parFile.contains("BSSN_EPSILON_CAKO_GAUGE")) {
        bssn::BSSN_EPSILON_CAKO_GAUGE =
            parFile["BSSN_EPSILON_CAKO_GAUGE"].as_floating();
    }
    if (parFile.contains("BSSN_EPSILON_CAKO_OTHER")) {
        bssn::BSSN_EPSILON_CAKO_OTHER =
            parFile["BSSN_EPSILON_CAKO_OTHER"].as_floating();
    }

    if (parFile.contains("BSSN_CAHD_C")) {
        bssn::BSSN_CAHD_C = parFile["BSSN_CAHD_C"].as_floating();
    }

    bssn::BSSN_LAMBDA[0]   = parFile["BSSN_LAMBDA"][0].as_integer();
    bssn::BSSN_LAMBDA[1]   = parFile["BSSN_LAMBDA"][1].as_integer();
    bssn::BSSN_LAMBDA[2]   = parFile["BSSN_LAMBDA"][2].as_integer();
    bssn::BSSN_LAMBDA[3]   = parFile["BSSN_LAMBDA"][3].as_integer();
    bssn::BSSN_LAMBDA_F[0] = parFile["BSSN_LAMBDA_F"][0].as_floating();
    bssn::BSSN_LAMBDA_F[1] = parFile["BSSN_LAMBDA_F"][1].as_floating();

    bssn::BSSN_XI[0]       = (unsigned int)parFile["BSSN_XI"][0].as_integer();
    bssn::BSSN_XI[1]       = (unsigned int)parFile["BSSN_XI"][1].as_integer();
    bssn::BSSN_XI[2]       = (unsigned int)parFile["BSSN_XI"][2].as_integer();

    if (parFile.contains("BSSN_ELE_ORDER"))
        bssn::BSSN_ELE_ORDER = parFile["BSSN_ELE_ORDER"].as_integer();
    bssn::CHI_FLOOR = parFile["CHI_FLOOR"].as_floating();
    bssn::BSSN_TRK0 = parFile["BSSN_TRK0"].as_floating();

    if (parFile.contains("DISSIPATION_TYPE"))
        bssn::DISSIPATION_TYPE = parFile["DISSIPATION_TYPE"].as_integer();

    bssn::KO_DISS_SIGMA     = parFile["KO_DISS_SIGMA"].as_floating();

    bssn::BSSN_ETA_R0       = parFile["BSSN_ETA_R0"].as_floating();
    bssn::BSSN_ETA_POWER[0] = parFile["BSSN_ETA_POWER"][0].as_floating();
    bssn::BSSN_ETA_POWER[1] = parFile["BSSN_ETA_POWER"][1].as_floating();

    bssn::BSSN_USE_WAVELET_TOL_FUNCTION =
        parFile["BSSN_USE_WAVELET_TOL_FUNCTION"].as_integer();
    bssn::BSSN_WAVELET_TOL     = parFile["BSSN_WAVELET_TOL"].as_floating();
    bssn::BSSN_WAVELET_TOL_MAX = parFile["BSSN_WAVELET_TOL_MAX"].as_floating();
    bssn::BSSN_WAVELET_TOL_FUNCTION_R0 =
        parFile["BSSN_WAVELET_TOL_FUNCTION_R0"].as_floating();
    bssn::BSSN_WAVELET_TOL_FUNCTION_R1 =
        parFile["BSSN_WAVELET_TOL_FUNCTION_R1"].as_floating();

    bssn::BSSN_NUM_REFINE_VARS = parFile["BSSN_NUM_REFINE_VARS"].as_integer();
    for (unsigned int i = 0; i < bssn::BSSN_NUM_REFINE_VARS; i++)
        bssn::BSSN_REFINE_VARIABLE_INDICES[i] =
            parFile["BSSN_REFINE_VARIABLE_INDICES"][i].as_integer();

    bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT =
        parFile["BSSN_NUM_EVOL_VARS_VTU_OUTPUT"].as_integer();
    bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT =
        parFile["BSSN_NUM_CONST_VARS_VTU_OUTPUT"].as_integer();

    for (unsigned int i = 0; i < bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT; i++)
        bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[i] =
            parFile["BSSN_VTU_OUTPUT_EVOL_INDICES"][i].as_integer();

    for (unsigned int i = 0; i < bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT; i++)
        bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i] =
            parFile["BSSN_VTU_OUTPUT_CONST_INDICES"][i].as_integer();

    if (parFile.contains("BSSN_CFL_FACTOR"))
        bssn::BSSN_CFL_FACTOR = parFile["BSSN_CFL_FACTOR"].as_floating();

    if (parFile.contains("BSSN_VTU_Z_SLICE_ONLY"))
        bssn::BSSN_VTU_Z_SLICE_ONLY =
            parFile["BSSN_VTU_Z_SLICE_ONLY"].as_boolean();

    if (parFile.contains("BSSN_GW_EXTRACT_FREQ"))
        bssn::BSSN_GW_EXTRACT_FREQ =
            parFile["BSSN_GW_EXTRACT_FREQ"].as_integer();
    else
        bssn::BSSN_GW_EXTRACT_FREQ =
            std::max(1u, bssn::BSSN_IO_OUTPUT_FREQ >> 1u);

    if (parFile.contains("BSSN_TIME_STEP_OUTPUT_FREQ")) {
        bssn::BSSN_TIME_STEP_OUTPUT_FREQ =
            parFile["BSSN_TIME_STEP_OUTPUT_FREQ"].as_integer();
    } else {
        bssn::BSSN_TIME_STEP_OUTPUT_FREQ = bssn::BSSN_GW_EXTRACT_FREQ;
    }

    if (parFile.contains("BSSN_BH1_AMR_R"))
        bssn::BSSN_BH1_AMR_R = parFile["BSSN_BH1_AMR_R"].as_floating();

    if (parFile.contains("BSSN_BH2_AMR_R"))
        bssn::BSSN_BH2_AMR_R = parFile["BSSN_BH2_AMR_R"].as_floating();

    if (parFile.contains("BSSN_AMR_R_RATIO")) {
        bssn::BSSN_AMR_R_RATIO = parFile["BSSN_AMR_R_RATIO"].as_floating();
    }

    if (parFile.contains("BSSN_BH1_MAX_LEV"))
        bssn::BSSN_BH1_MAX_LEV = parFile["BSSN_BH1_MAX_LEV"].as_integer();
    else
        bssn::BSSN_BH1_MAX_LEV = bssn::BSSN_MAXDEPTH;

    if (parFile.contains("BSSN_BH2_MAX_LEV"))
        bssn::BSSN_BH2_MAX_LEV = parFile["BSSN_BH2_MAX_LEV"].as_integer();
    else
        bssn::BSSN_BH2_MAX_LEV = bssn::BSSN_MAXDEPTH;

    if (parFile.contains("BSSN_INIT_GRID_ITER"))
        bssn::BSSN_INIT_GRID_ITER = parFile["BSSN_INIT_GRID_ITER"].as_integer();

    if (parFile.contains("BSSN_GW_REFINE_WTOL"))
        bssn::BSSN_GW_REFINE_WTOL =
            parFile["BSSN_GW_REFINE_WTOL"].as_floating();

    if (parFile.contains("BSSN_MINDEPTH"))
        bssn::BSSN_MINDEPTH = parFile["BSSN_MINDEPTH"].as_integer();

    if (parFile.contains("BSSN_BH1_CONSTRAINT_R"))
        bssn::BSSN_BH1_CONSTRAINT_R =
            parFile["BSSN_BH1_CONSTRAINT_R"].as_floating();

    if (parFile.contains("BSSN_BH2_CONSTRAINT_R"))
        bssn::BSSN_BH2_CONSTRAINT_R =
            parFile["BSSN_BH2_CONSTRAINT_R"].as_floating();

    if (parFile.contains("BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE"))
        bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE =
            parFile["BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE"].as_boolean();

    if (parFile.contains("BSSN_NYQUIST_M")) {
        bssn::BSSN_NYQUIST_M = parFile["BSSN_NYQUIST_M"].as_integer();
    }

    if (parFile.contains("BSSN_SCALE_VTU_AND_GW_EXTRACTION")) {
        bssn::BSSN_SCALE_VTU_AND_GW_EXTRACTION =
            parFile["BSSN_SCALE_VTU_AND_GW_EXTRACTION"].as_boolean();
    }
    bssn::BSSN_IO_OUTPUT_FREQ_TRUE  = bssn::BSSN_IO_OUTPUT_FREQ;
    bssn::BSSN_GW_EXTRACT_FREQ_TRUE = bssn::BSSN_GW_EXTRACT_FREQ;

    if (parFile.contains("BSSN_SSL_SIGMA")) {
        bssn::BSSN_SSL_SIGMA = parFile["BSSN_SSL_SIGMA"].as_floating();
    }

    if (parFile.contains("BSSN_SSL_H")) {
        bssn::BSSN_SSL_H = parFile["BSSN_SSL_H"].as_floating();
    }

    /* Parameters for TPID */
    TPID::target_M_plus  = parFile["TPID_TARGET_M_PLUS"].as_floating();
    TPID::target_M_minus = parFile["TPID_TARGET_M_MINUS"].as_floating();
    TPID::par_m_plus     = TPID::target_M_plus;
    TPID::par_m_minus    = TPID::target_M_minus;
    TPID::par_b          = parFile["TPID_PAR_B"].as_floating();

    TPID::par_P_plus[0] =
        bssn::BH1.getVx();  // parFile["TPID_PAR_P_PLUS"]["X"];
    TPID::par_P_plus[1] =
        bssn::BH1.getVy();  // parFile["TPID_PAR_P_PLUS"]["Y"];
    TPID::par_P_plus[2] =
        bssn::BH1.getVz();  // parFile["TPID_PAR_P_PLUS"]["Z"];

    TPID::par_P_minus[0] =
        bssn::BH2.getVx();  // parFile["TPID_PAR_P_MINUS"]["X"];
    TPID::par_P_minus[1] =
        bssn::BH2.getVy();  // parFile["TPID_PAR_P_MINUS"]["Y"];
    TPID::par_P_minus[2] =
        bssn::BH2.getVz();  // parFile["TPID_PAR_P_MINUS"]["Z"];

    TPID::par_S_plus[0] =
        bssn::BH1.getBHSpin() * sin(bssn::BH1.getBHSpinTheta()) *
        cos(bssn::BH1.getBHSpinPhi());  // parFile["TPID_PAR_S_PLUS"]["X"];
    TPID::par_S_plus[1] =
        bssn::BH1.getBHSpin() * sin(bssn::BH1.getBHSpinTheta()) *
        sin(bssn::BH1.getBHSpinPhi());  // parFile["TPID_PAR_S_PLUS"]["Y"];
    TPID::par_S_plus[2] =
        bssn::BH1.getBHSpin() *
        cos(bssn::BH1.getBHSpinTheta());  // parFile["TPID_PAR_S_PLUS"]["Z"];

    TPID::par_S_minus[0] =
        bssn::BH2.getBHSpin() * sin(bssn::BH2.getBHSpinTheta()) *
        cos(bssn::BH2.getBHSpinPhi());  // parFile["TPID_PAR_S_MINUS"]["X"];
    TPID::par_S_minus[1] =
        bssn::BH2.getBHSpin() * sin(bssn::BH2.getBHSpinTheta()) *
        sin(bssn::BH2.getBHSpinPhi());  // parFile["TPID_PAR_S_MINUS"]["Y"];
    TPID::par_S_minus[2] =
        bssn::BH2.getBHSpin() *
        cos(bssn::BH2.getBHSpinTheta());  // parFile["TPID_PAR_S_MINUS"]["Z"];

    TPID::center_offset[0] = parFile["TPID_CENTER_OFFSET"]["X"].as_floating();
    TPID::center_offset[1] = parFile["TPID_CENTER_OFFSET"]["Y"].as_floating();
    TPID::center_offset[2] = parFile["TPID_CENTER_OFFSET"]["Z"].as_floating();

    TPID::initial_lapse_psi_exponent =
        parFile["TPID_INITIAL_LAPSE_PSI_EXPONENT"].as_floating();
    TPID::npoints_A      = parFile["TPID_NPOINTS_A"].as_integer();
    TPID::npoints_B      = parFile["TPID_NPOINTS_B"].as_integer();
    TPID::npoints_phi    = parFile["TPID_NPOINTS_PHI"].as_integer();

    TPID::give_bare_mass = parFile["TPID_GIVE_BARE_MASS"].as_integer();
    TPID::initial_lapse  = parFile["INITIAL_LAPSE"].as_integer();
    TPID::solve_momentum_constraint =
        parFile["TPID_SOLVE_MOMENTUM_CONSTRAINT"].as_integer();
    TPID::grid_setup_method = parFile["TPID_GRID_SETUP_METHOD"].as_integer();
    TPID::verbose           = parFile["TPID_VERBOSE"].as_integer();
    TPID::adm_tol           = parFile["TPID_ADM_TOL"].as_floating();
    TPID::Newton_tol        = parFile["TPID_NEWTON_TOL"].as_floating();

    if (parFile.contains("TPID_FILEPREFIX"))
        TPID::FILE_PREFIX = parFile["TPID_FILEPREFIX"].as_string();

    if (parFile.contains("TPID_REPLACE_LAPSE_WITH_SQRT_CHI"))
        TPID::replace_lapse_with_sqrt_chi =
            parFile["TPID_REPLACE_LAPSE_WITH_SQRT_CHI"].as_boolean();

    if (parFile.contains("EXTRACTION_VAR_ID"))
        BHLOC::EXTRACTION_VAR_ID = parFile["EXTRACTION_VAR_ID"].as_integer();

    if (parFile.contains("EXTRACTION_TOL"))
        BHLOC::EXTRACTION_TOL = parFile["EXTRACTION_TOL"].as_floating();

    GW::BSSN_GW_NUM_RADAII = parFile["BSSN_GW_NUM_RADAII"].as_integer();
    GW::BSSN_GW_NUM_LMODES = parFile["BSSN_GW_NUM_LMODES"].as_integer();

    for (unsigned int i = 0; i < GW::BSSN_GW_NUM_RADAII; i++)
        GW::BSSN_GW_RADAII[i] = parFile["BSSN_GW_RADAII"][i].as_floating();

    for (unsigned int i = 0; i < GW::BSSN_GW_NUM_LMODES; i++)
        GW::BSSN_GW_L_MODES[i] = parFile["BSSN_GW_L_MODES"][i].as_integer();

    if (parFile.contains("BSSN_EH_COARSEN_VAL"))
        bssn::BSSN_EH_COARSEN_VAL =
            parFile["BSSN_EH_COARSEN_VAL"].as_floating();

    if (parFile.contains("BSSN_EH_REFINE_VAL"))
        bssn::BSSN_EH_REFINE_VAL = parFile["BSSN_EH_REFINE_VAL"].as_floating();

    if (parFile.contains("BSSN_REFINEMENT_MODE"))
        bssn::BSSN_REFINEMENT_MODE = static_cast<bssn::RefinementMode>(
            parFile["BSSN_REFINEMENT_MODE"].as_integer());

    BSSN_OCTREE_MAX[0] = (double)(1u << bssn::BSSN_MAXDEPTH);
    BSSN_OCTREE_MAX[1] = (double)(1u << bssn::BSSN_MAXDEPTH);
    BSSN_OCTREE_MAX[2] = (double)(1u << bssn::BSSN_MAXDEPTH);

    BSSN_COMPD_MIN[0]  = bssn::BSSN_GRID_MIN_X;
    BSSN_COMPD_MIN[1]  = bssn::BSSN_GRID_MIN_Y;
    BSSN_COMPD_MIN[2]  = bssn::BSSN_GRID_MIN_Z;

    BSSN_COMPD_MAX[0]  = bssn::BSSN_GRID_MAX_X;
    BSSN_COMPD_MAX[1]  = bssn::BSSN_GRID_MAX_Y;
    BSSN_COMPD_MAX[2]  = bssn::BSSN_GRID_MAX_Z;

    if (BSSN_NUM_REFINE_VARS > BSSN_NUM_VARS) {
        std::cout << "Error[parameter file]: Number of refine variables should "
                     "be less than number of BSSN_NUM_VARS"
                  << std::endl;
        exit(0);
    }
    if (BSSN_NUM_EVOL_VARS_VTU_OUTPUT > BSSN_NUM_VARS) {
        std::cout << "Error[parameter file]: Number of evolution VTU variables "
                     "should be less than number of BSSN_NUM_VARS"
                  << std::endl;
        exit(0);
    }
    if (BSSN_NUM_CONST_VARS_VTU_OUTPUT > BSSN_CONSTRAINT_NUM_VARS) {
        std::cout
            << "Error[parameter file]: Number of constraint VTU variables "
               "should be less than number of BSSN_CONSTRAINT_NUM_VARS"
            << std::endl;
        exit(0);
    }

    BSSN_PADDING_WIDTH = BSSN_ELE_ORDER >> 1u;
    bssn::BSSN_BH_LOC[0] =
        Point(BH1.getBHCoordX(), BH1.getBHCoordY(), BH1.getBHCoordZ());
    bssn::BSSN_BH_LOC[1] =
        Point(BH2.getBHCoordX(), BH2.getBHCoordY(), BH2.getBHCoordZ());
    bssn::BSSN_BH1_MASS = BH1.getBHMass();
    bssn::BSSN_BH2_MASS = BH2.getBHMass();

    // AH parameters
    if (parFile.contains("AEH_LMAX"))
        AEH::AEH_LMAX = parFile["AEH_LMAX"].as_integer();

    if (parFile.contains("AEH_Q_THETA"))
        AEH::AEH_Q_THETA = parFile["AEH_Q_THETA"].as_integer();

    if (parFile.contains("AEH_Q_PHI"))
        AEH::AEH_Q_PHI = parFile["AEH_Q_PHI"].as_integer();

    if (parFile.contains("AEH_MAXITER"))
        AEH::AEH_MAXITER = parFile["AEH_MAXITER"].as_integer();

    if (parFile.contains("AEH_ATOL"))
        AEH::AEH_ATOL = parFile["AEH_ATOL"].as_floating();

    if (parFile.contains("AEH_RTOL"))
        AEH::AEH_RTOL = parFile["AEH_RTOL"].as_floating();

    if (parFile.contains("AEH_SOLVER_FREQ"))
        AEH::AEH_SOLVER_FREQ = parFile["AEH_SOLVER_FREQ"].as_integer();

    if (parFile.contains("AEH_ALPHA"))
        AEH::AEH_ALPHA = parFile["AEH_ALPHA"].as_floating();

    if (parFile.contains("AEH_BETA"))
        AEH::AEH_BETA = parFile["AEH_BETA"].as_floating();

    MPI_Barrier(comm);
}

}  // namespace bssn

namespace TPID {
double target_M_plus              = 1.0;
double target_M_minus             = 1.0;
double par_m_plus                 = 1.0;
double par_m_minus                = 1.0;
double par_b                      = 4.0;
double par_P_plus[3]              = {0.0, 0.0, 0.0};
double par_P_minus[3]             = {0.0, 0.0, 0.0};
double par_S_plus[3]              = {0.0, 0.0, 0.0};
double par_S_minus[3]             = {0.0, 0.0, 0.0};
double center_offset[3]           = {0.0, 0.0, 0.00014142135623730951};
double initial_lapse_psi_exponent = -2;
int npoints_A                     = 30;
int npoints_B                     = 30;
int npoints_phi                   = 16;
int give_bare_mass                = 0;
int initial_lapse                 = 2;
int solve_momentum_constraint     = 1;
int grid_setup_method             = 1;
int verbose                       = 1;
double adm_tol                    = 1.0e-10;
double Newton_tol                 = 1.0e-10;
std::string FILE_PREFIX           = "tpid";
bool replace_lapse_with_sqrt_chi  = false;
}  // namespace TPID

namespace BHLOC {
unsigned int EXTRACTION_VAR_ID = bssn::VAR::U_ALPHA;
double EXTRACTION_TOL          = 0.3;
}  // namespace BHLOC

namespace GW {

unsigned int BSSN_GW_NUM_RADAII;

unsigned int BSSN_GW_NUM_LMODES;

double BSSN_GW_RADAII[BSSN_GW_MAX_RADAII];

unsigned int BSSN_GW_SPIN = 2;

unsigned int BSSN_GW_L_MODES[BSSN_GW_MAX_LMODES];
}  // namespace GW

namespace AEH {
unsigned int AEH_LMAX        = 6;
unsigned int AEH_Q_THETA     = 32;
unsigned int AEH_Q_PHI       = 32;
unsigned int AEH_MAXITER     = 50;
double AEH_ATOL              = 1e-8;
double AEH_RTOL              = 1e-8;
unsigned int AEH_SOLVER_FREQ = 0;

double AEH_ALPHA             = 1.0;
double AEH_BETA              = 0.1;

}  // namespace AEH
