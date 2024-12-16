//
// Created by milinda on 7/26/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Contains utility functions for BSSN simulation.
 */
//

#include "grUtils.h"

#include <tuple>

#include "base.h"
#include "git_version_and_date.h"
#include "parameters.h"

namespace bssn {

void printGitInformation(int rank, std::vector<std::string> arg_s) {
    if (!rank) {
        std::cout << YLW << "  COMPILED ON  -  " << compile_info::compileDate
                  << NRM << std::endl;
        std::cout << YLW << "  LATEST GIT HASH - " << compile_info::currGitHash
                  << compile_info::dirtyStatus << NRM << std::endl;
    }

    for (size_t ii = 1; ii < arg_s.size(); ++ii) {
        if (arg_s[ii] == "--compile-info") {
            if (!rank)
                std::cout << "Compile info only flag found, exiting..." << NRM
                          << std::endl
                          << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }
}

void readParamFile(const char* fName, MPI_Comm comm) {
    std::string fNameStr(fName);
    std::string tomlSuffix = ".toml";

    if (fNameStr.size() >= tomlSuffix.size() &&
        fNameStr.compare(fNameStr.size() - tomlSuffix.size(), tomlSuffix.size(),
                         tomlSuffix) == 0) {
        // we found a toml file!
        readParamTOMLFile(fName, comm);

    } else {
        // fall back to JSON file reading
        readParamJSONFile(fName, comm);
    }
}

void readParamJSONFile(const char* fName, MPI_Comm comm) {
    json parFile;
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    unsigned int vtu_len;
    unsigned int chp_len;
    unsigned int prf_len;
    unsigned int tpf_len;

    std::ifstream infile(fName);
    if (!infile) {
        std::cout << fName << " parameter file open failed by " << rank << "/"
                  << npes << std::endl;
    }
    infile >> parFile;

    bssn::BSSN_IO_OUTPUT_FREQ   = parFile["BSSN_IO_OUTPUT_FREQ"];
    bssn::BSSN_REMESH_TEST_FREQ = parFile["BSSN_REMESH_TEST_FREQ"];
    bssn::BSSN_CHECKPT_FREQ     = parFile["BSSN_CHECKPT_FREQ"];
    bssn::BSSN_IO_OUTPUT_GAP    = parFile["BSSN_IO_OUTPUT_GAP"];
    bssn::BSSN_VTU_FILE_PREFIX =
        parFile["BSSN_VTU_FILE_PREFIX"].get<std::string>();
    bssn::BSSN_CHKPT_FILE_PREFIX =
        parFile["BSSN_CHKPT_FILE_PREFIX"].get<std::string>();

    bssn::BSSN_PROFILE_FILE_PREFIX =
        parFile["BSSN_PROFILE_FILE_PREFIX"].get<std::string>();
    bssn::BSSN_RESTORE_SOLVER = parFile["BSSN_RESTORE_SOLVER"];
    bssn::BSSN_ID_TYPE        = parFile["BSSN_ID_TYPE"];

    bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY =
        parFile["BSSN_ENABLE_BLOCK_ADAPTIVITY"];
    bssn::BSSN_BLK_MIN_X           = parFile["BSSN_BLK_MIN_X"];
    bssn::BSSN_BLK_MIN_Y           = parFile["BSSN_BLK_MIN_Y"];
    bssn::BSSN_BLK_MIN_Z           = parFile["BSSN_BLK_MIN_Z"];
    bssn::BSSN_BLK_MAX_X           = parFile["BSSN_BLK_MAX_X"];
    bssn::BSSN_BLK_MAX_Y           = parFile["BSSN_BLK_MAX_Y"];
    bssn::BSSN_BLK_MAX_Z           = parFile["BSSN_BLK_MAX_Z"];

    bssn::BSSN_DENDRO_GRAIN_SZ     = parFile["BSSN_DENDRO_GRAIN_SZ"];
    bssn::BSSN_ASYNC_COMM_K        = parFile["BSSN_ASYNC_COMM_K"];
    bssn::BSSN_DENDRO_AMR_FAC      = parFile["BSSN_DENDRO_AMR_FAC"];
    bssn::BSSN_LOAD_IMB_TOL        = parFile["BSSN_LOAD_IMB_TOL"];
    bssn::BSSN_RK_TIME_BEGIN       = parFile["BSSN_RK_TIME_BEGIN"];
    bssn::BSSN_RK_TIME_END         = parFile["BSSN_RK_TIME_END"];
    bssn::BSSN_RK_TYPE             = parFile["BSSN_RK_TYPE"];
    bssn::BSSN_RK45_TIME_STEP_SIZE = parFile["BSSN_RK45_TIME_STEP_SIZE"];
    bssn::BSSN_RK45_DESIRED_TOL    = parFile["BSSN_RK45_DESIRED_TOL"];
    bssn::BSSN_DIM                 = parFile["BSSN_DIM"];
    bssn::BSSN_MAXDEPTH            = parFile["BSSN_MAXDEPTH"];

    bssn::BH1                      = BH(
        (double)parFile["BSSN_BH1"]["MASS"], (double)parFile["BSSN_BH1"]["X"],
        (double)parFile["BSSN_BH1"]["Y"], (double)parFile["BSSN_BH1"]["Z"],
        (double)parFile["BSSN_BH1"]["V_X"], (double)parFile["BSSN_BH1"]["V_Y"],
        (double)parFile["BSSN_BH1"]["V_Z"], (double)parFile["BSSN_BH1"]["SPIN"],
        (double)parFile["BSSN_BH1"]["SPIN_THETA"],
        (double)parFile["BSSN_BH1"]["SPIN_PHI"]);
    bssn::BH2 = BH(
        (double)parFile["BSSN_BH2"]["MASS"], (double)parFile["BSSN_BH2"]["X"],
        (double)parFile["BSSN_BH2"]["Y"], (double)parFile["BSSN_BH2"]["Z"],
        (double)parFile["BSSN_BH2"]["V_X"], (double)parFile["BSSN_BH2"]["V_Y"],
        (double)parFile["BSSN_BH2"]["V_Z"], (double)parFile["BSSN_BH2"]["SPIN"],
        (double)parFile["BSSN_BH2"]["SPIN_THETA"],
        (double)parFile["BSSN_BH2"]["SPIN_PHI"]);

    bssn::BSSN_GRID_MIN_X = parFile["BSSN_GRID_MIN_X"];
    bssn::BSSN_GRID_MAX_X = parFile["BSSN_GRID_MAX_X"];
    bssn::BSSN_GRID_MIN_Y = parFile["BSSN_GRID_MIN_Y"];
    bssn::BSSN_GRID_MAX_Y = parFile["BSSN_GRID_MAX_Y"];
    bssn::BSSN_GRID_MIN_Z = parFile["BSSN_GRID_MIN_Z"];
    bssn::BSSN_GRID_MAX_Z = parFile["BSSN_GRID_MAX_Z"];

    bssn::ETA_CONST       = parFile["ETA_CONST"];
    bssn::ETA_R0          = parFile["ETA_R0"];
    bssn::ETA_DAMPING     = parFile["ETA_DAMPING"];
    bssn::ETA_DAMPING_EXP = parFile["ETA_DAMPING_EXP"];

    if (parFile.find("RIT_ETA_FUNCTION") != parFile.end()) {
        bssn::RIT_ETA_FUNCTION = parFile["RIT_ETA_FUNCTION"];
    }
    if (parFile.find("RIT_ETA_OUTER") != parFile.end()) {
        bssn::RIT_ETA_OUTER = parFile["RIT_ETA_OUTER"];
    }
    if (parFile.find("RIT_ETA_CENTRAL") != parFile.end()) {
        bssn::RIT_ETA_CENTRAL = parFile["RIT_ETA_CENTRAL"];
    }
    if (parFile.find("RIT_ETA_WIDTH") != parFile.end()) {
        bssn::RIT_ETA_WIDTH = parFile["RIT_ETA_WIDTH"];
    }

    if (parFile.find("BSSN_AMR_R_RATIO") != parFile.end()) {
        bssn::BSSN_AMR_R_RATIO = parFile["BSSN_AMR_R_RATIO"];
    }

    bssn::BSSN_LAMBDA[0] =
        (unsigned int)parFile["BSSN_LAMBDA"]["BSSN_LAMBDA_1"];
    bssn::BSSN_LAMBDA[1] =
        (unsigned int)parFile["BSSN_LAMBDA"]["BSSN_LAMBDA_2"];
    bssn::BSSN_LAMBDA[2] =
        (unsigned int)parFile["BSSN_LAMBDA"]["BSSN_LAMBDA_3"];
    bssn::BSSN_LAMBDA[3] =
        (unsigned int)parFile["BSSN_LAMBDA"]["BSSN_LAMBDA_4"];
    bssn::BSSN_LAMBDA_F[0] = parFile["BSSN_LAMBDA_F"]["BSSN_LAMBDA_F0"];
    bssn::BSSN_LAMBDA_F[1] = parFile["BSSN_LAMBDA_F"]["BSSN_LAMBDA_F1"];

    bssn::BSSN_XI[0]       = (unsigned int)parFile["BSSN_XI"]["BSSN_XI_0"];
    bssn::BSSN_XI[1]       = (unsigned int)parFile["BSSN_XI"]["BSSN_XI_1"];
    bssn::BSSN_XI[2]       = (unsigned int)parFile["BSSN_XI"]["BSSN_XI_2"];

    if (parFile.find("BSSN_ELE_ORDER") != parFile.end())
        bssn::BSSN_ELE_ORDER = parFile["BSSN_ELE_ORDER"];

    bssn::CHI_FLOOR = parFile["CHI_FLOOR"];
    bssn::BSSN_TRK0 = parFile["BSSN_TRK0"];
    if (parFile.find("DISSIPATION_TYPE") != parFile.end()) {
        bssn::DISSIPATION_TYPE = parFile["DISSIPATION_TYPE"];
    }
    bssn::KO_DISS_SIGMA = parFile["KO_DISS_SIGMA"];

    if (parFile.find("BSSN_KO_SIGMA_SCALE_BY_CONFORMAL") != parFile.end()) {
        bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL =
            parFile["BSSN_KO_SIGMA_SCALE_BY_CONFORMAL"];
    }

    if (parFile.find("BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY") !=
        parFile.end()) {
        bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY =
            parFile["BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY"];
    }

    if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY) {
        bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL = false;
    }
    if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL) {
        bssn::BSSN_CAKO_ENABLED = true;
    }

    if (parFile.find("BSSN_EPSILON_CAKO_GAUGE") != parFile.end()) {
        bssn::BSSN_EPSILON_CAKO_GAUGE = parFile["BSSN_EPSILON_CAKO_GAUGE"];
    }
    if (parFile.find("BSSN_EPSILON_CAKO_OTHER") != parFile.end()) {
        bssn::BSSN_EPSILON_CAKO_OTHER = parFile["BSSN_EPSILON_CAKO_OTHER"];
    }

    if (parFile.find("BSSN_CAHD_C") != parFile.end()) {
        bssn::BSSN_CAHD_C = parFile["BSSN_CAHD_C"];
    }

    // Parameters for eta_damping function
    bssn::BSSN_ETA_R0       = parFile["BSSN_ETA_R0"];
    bssn::BSSN_ETA_POWER[0] = parFile["BSSN_ETA_POWER"]["BSSN_ETA_POWER_1"];
    bssn::BSSN_ETA_POWER[1] = parFile["BSSN_ETA_POWER"]["BSSN_ETA_POWER_2"];

    bssn::BSSN_USE_WAVELET_TOL_FUNCTION =
        parFile["BSSN_USE_WAVELET_TOL_FUNCTION"];
    bssn::BSSN_WAVELET_TOL     = parFile["BSSN_WAVELET_TOL"];
    bssn::BSSN_WAVELET_TOL_MAX = parFile["BSSN_WAVELET_TOL_MAX"];
    bssn::BSSN_WAVELET_TOL_FUNCTION_R0 =
        parFile["BSSN_WAVELET_TOL_FUNCTION_R0"];
    bssn::BSSN_WAVELET_TOL_FUNCTION_R1 =
        parFile["BSSN_WAVELET_TOL_FUNCTION_R1"];

    bssn::BSSN_NUM_REFINE_VARS = parFile["BSSN_NUM_REFINE_VARS"];
    for (unsigned int i = 0; i < bssn::BSSN_NUM_REFINE_VARS; i++)
        bssn::BSSN_REFINE_VARIABLE_INDICES[i] =
            parFile["BSSN_REFINE_VARIABLE_INDICES"][i];

    bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT =
        parFile["BSSN_NUM_EVOL_VARS_VTU_OUTPUT"];
    bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT =
        parFile["BSSN_NUM_CONST_VARS_VTU_OUTPUT"];

    for (unsigned int i = 0; i < bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT; i++)
        bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[i] =
            parFile["BSSN_VTU_OUTPUT_EVOL_INDICES"][i];

    for (unsigned int i = 0; i < bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT; i++)
        bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i] =
            parFile["BSSN_VTU_OUTPUT_CONST_INDICES"][i];

    if (parFile.find("BSSN_CFL_FACTOR") != parFile.end()) {
        bssn::BSSN_CFL_FACTOR = parFile["BSSN_CFL_FACTOR"];
    }

    if (parFile.find("BSSN_VTU_Z_SLICE_ONLY") != parFile.end())
        bssn::BSSN_VTU_Z_SLICE_ONLY = parFile["BSSN_VTU_Z_SLICE_ONLY"];

    if (parFile.find("BSSN_GW_EXTRACT_FREQ") != parFile.end()) {
        bssn::BSSN_GW_EXTRACT_FREQ = parFile["BSSN_GW_EXTRACT_FREQ"];
    } else {
        bssn::BSSN_GW_EXTRACT_FREQ =
            std::max(1u, bssn::BSSN_IO_OUTPUT_FREQ >> 1u);
    }

    if (parFile.find("BSSN_TIME_STEP_OUTPUT_FREQ") != parFile.end()) {
        bssn::BSSN_TIME_STEP_OUTPUT_FREQ =
            parFile["BSSN_TIME_STEP_OUTPUT_FREQ"];
    } else {
        bssn::BSSN_TIME_STEP_OUTPUT_FREQ = bssn::BSSN_GW_EXTRACT_FREQ;
    }

    if (parFile.find("BSSN_BH1_AMR_R") != parFile.end())
        bssn::BSSN_BH1_AMR_R = parFile["BSSN_BH1_AMR_R"];

    if (parFile.find("BSSN_BH2_AMR_R") != parFile.end())
        bssn::BSSN_BH2_AMR_R = parFile["BSSN_BH2_AMR_R"];

    if (parFile.find("BSSN_AMR_R_RATIO") != parFile.end())
        bssn::BSSN_AMR_R_RATIO = parFile["BSSN_AMR_R_RATIO"];

    if (parFile.find("BSSN_BH1_MAX_LEV") != parFile.end())
        bssn::BSSN_BH1_MAX_LEV = parFile["BSSN_BH1_MAX_LEV"];
    else
        bssn::BSSN_BH1_MAX_LEV = bssn::BSSN_MAXDEPTH;

    if (parFile.find("BSSN_BH2_MAX_LEV") != parFile.end())
        bssn::BSSN_BH2_MAX_LEV = parFile["BSSN_BH2_MAX_LEV"];
    else
        bssn::BSSN_BH2_MAX_LEV = bssn::BSSN_MAXDEPTH;

    if (parFile.find("BSSN_INIT_GRID_ITER") != parFile.end())
        bssn::BSSN_INIT_GRID_ITER = parFile["BSSN_INIT_GRID_ITER"];

    if (parFile.find("BSSN_GW_REFINE_WTOL") != parFile.end())
        bssn::BSSN_GW_REFINE_WTOL = parFile["BSSN_GW_REFINE_WTOL"];

    if (parFile.find("BSSN_MINDEPTH") != parFile.end())
        bssn::BSSN_MINDEPTH = parFile["BSSN_MINDEPTH"];

    if (parFile.find("BSSN_BH1_CONSTRAINT_R") != parFile.end())
        bssn::BSSN_BH1_CONSTRAINT_R = parFile["BSSN_BH1_CONSTRAINT_R"];

    if (parFile.find("BSSN_BH2_CONSTRAINT_R") != parFile.end())
        bssn::BSSN_BH2_CONSTRAINT_R = parFile["BSSN_BH2_CONSTRAINT_R"];

    if (parFile.find("BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE") !=
        parFile.end())
        bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE =
            parFile["BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE"];

    /* Parameters for TPID */
    TPID::target_M_plus  = parFile["TPID_TARGET_M_PLUS"];
    TPID::target_M_minus = parFile["TPID_TARGET_M_MINUS"];
    TPID::par_m_plus     = TPID::target_M_plus;
    TPID::par_m_minus    = TPID::target_M_minus;
    TPID::par_b          = parFile["TPID_PAR_B"];

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

    TPID::center_offset[0] = parFile["TPID_CENTER_OFFSET"]["X"];
    TPID::center_offset[1] = parFile["TPID_CENTER_OFFSET"]["Y"];
    TPID::center_offset[2] = parFile["TPID_CENTER_OFFSET"]["Z"];

    TPID::initial_lapse_psi_exponent =
        parFile["TPID_INITIAL_LAPSE_PSI_EXPONENT"];
    TPID::npoints_A                 = parFile["TPID_NPOINTS_A"];
    TPID::npoints_B                 = parFile["TPID_NPOINTS_B"];
    TPID::npoints_phi               = parFile["TPID_NPOINTS_PHI"];

    TPID::give_bare_mass            = parFile["TPID_GIVE_BARE_MASS"];
    TPID::initial_lapse             = parFile["INITIAL_LAPSE"];
    TPID::solve_momentum_constraint = parFile["TPID_SOLVE_MOMENTUM_CONSTRAINT"];
    TPID::grid_setup_method         = parFile["TPID_GRID_SETUP_METHOD"];
    TPID::verbose                   = parFile["TPID_VERBOSE"];
    TPID::adm_tol                   = parFile["TPID_ADM_TOL"];
    TPID::Newton_tol                = parFile["TPID_NEWTON_TOL"];

    if (parFile.find("TPID_FILEPREFIX") != parFile.end())
        TPID::FILE_PREFIX = parFile["TPID_FILEPREFIX"].get<std::string>();

    if (parFile.find("TPID_REPLACE_LAPSE_WITH_SQRT_CHI") != parFile.end())
        TPID::replace_lapse_with_sqrt_chi =
            parFile["TPID_REPLACE_LAPSE_WITH_SQRT_CHI"];

    if (parFile.find("EXTRACTION_VAR_ID") != parFile.end())
        BHLOC::EXTRACTION_VAR_ID = parFile["EXTRACTION_VAR_ID"];

    if (parFile.find("EXTRACTION_TOL") != parFile.end())
        BHLOC::EXTRACTION_TOL = parFile["EXTRACTION_TOL"];

    if (parFile.find("BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE") !=
        parFile.end())
        bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE =
            parFile["BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE"];

    GW::BSSN_GW_NUM_RADAII = parFile["BSSN_GW_NUM_RADAII"];
    GW::BSSN_GW_NUM_LMODES = parFile["BSSN_GW_NUM_LMODES"];

    for (unsigned int i = 0; i < GW::BSSN_GW_NUM_RADAII; i++)
        GW::BSSN_GW_RADAII[i] = parFile["BSSN_GW_RADAII"][i];

    for (unsigned int i = 0; i < GW::BSSN_GW_NUM_LMODES; i++)
        GW::BSSN_GW_L_MODES[i] = parFile["BSSN_GW_L_MODES"][i];

    if (parFile.find("BSSN_EH_COARSEN_VAL") != parFile.end())
        bssn::BSSN_EH_COARSEN_VAL = parFile["BSSN_EH_COARSEN_VAL"];

    if (parFile.find("BSSN_EH_REFINE_VAL") != parFile.end())
        bssn::BSSN_EH_REFINE_VAL = parFile["BSSN_EH_REFINE_VAL"];

    if (parFile.find("BSSN_REFINEMENT_MODE") != parFile.end())
        bssn::BSSN_REFINEMENT_MODE =
            static_cast<bssn::RefinementMode>(parFile["BSSN_REFINEMENT_MODE"]);

    if (parFile.find("BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE") !=
        parFile.end()) {
        bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE =
            parFile["BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE"];
    }

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
        std::cout << "Error[parameter file]: Number of refine variables "
                     "should be less than number of BSSN_NUM_VARS"
                  << std::endl;
        exit(0);
    }
    if (BSSN_NUM_EVOL_VARS_VTU_OUTPUT > BSSN_NUM_VARS) {
        std::cout << "Error[parameter file]: Number of evolution VTU "
                     "variables should be less than number of BSSN_NUM_VARS"
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
    if (parFile.find("AEH_LMAX") != parFile.end())
        AEH::AEH_LMAX = parFile["AEH_LMAX"];

    if (parFile.find("AEH_Q_THETA") != parFile.end())
        AEH::AEH_Q_THETA = parFile["AEH_Q_THETA"];

    if (parFile.find("AEH_Q_PHI") != parFile.end())
        AEH::AEH_Q_PHI = parFile["AEH_Q_PHI"];

    if (parFile.find("AEH_MAXITER") != parFile.end())
        AEH::AEH_MAXITER = parFile["AEH_MAXITER"];

    if (parFile.find("AEH_ATOL") != parFile.end())
        AEH::AEH_ATOL = parFile["AEH_ATOL"];

    if (parFile.find("AEH_RTOL") != parFile.end())
        AEH::AEH_RTOL = parFile["AEH_RTOL"];

    if (parFile.find("AEH_SOLVER_FREQ") != parFile.end())
        AEH::AEH_SOLVER_FREQ = parFile["AEH_SOLVER_FREQ"];

    if (parFile.find("AEH_ALPHA") != parFile.end())
        AEH::AEH_ALPHA = parFile["AEH_ALPHA"];

    if (parFile.find("AEH_BETA") != parFile.end())
        AEH::AEH_BETA = parFile["AEH_BETA"];

    MPI_Barrier(comm);
}

void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (rank == root) {
        sout << "parameters read: " << std::endl;
        sout << YLW << "\tnpes :" << npes << NRM << std::endl;
        sout << YLW << "\tBSSN_DIM :" << bssn::BSSN_DIM << NRM << std::endl;
        sout << YLW << "\tBSSN_ELE_ORDER :" << bssn::BSSN_ELE_ORDER << NRM
             << std::endl;
        sout << YLW << "\tBSSN_PADDING_WIDTH :" << bssn::BSSN_PADDING_WIDTH
             << NRM << std::endl;
        sout << YLW << "\tBSSN_CFL_FACTOR :" << bssn::BSSN_CFL_FACTOR << NRM
             << std::endl;
        sout << YLW << "\tBSSN_IO_OUTPUT_FREQ :" << bssn::BSSN_IO_OUTPUT_FREQ
             << NRM << std::endl;
        sout << YLW << "\tBSSN_GW_EXTRACT_FREQ :" << bssn::BSSN_GW_EXTRACT_FREQ
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_REMESH_TEST_FREQ :" << bssn::BSSN_REMESH_TEST_FREQ
             << NRM << std::endl;
        sout << YLW << "\tBSSN_CHECKPT_FREQ :" << bssn::BSSN_CHECKPT_FREQ << NRM
             << std::endl;
        sout << YLW << "\tBSSN_RESTORE_SOLVER :" << bssn::BSSN_RESTORE_SOLVER
             << NRM << std::endl;
        sout << YLW << "\tBSSN_ENABLE_BLOCK_ADAPTIVITY :"
             << bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY << NRM << std::endl;
        sout << YLW << "\tBSSN_VTU_FILE_PREFIX :" << bssn::BSSN_VTU_FILE_PREFIX
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_CHKPT_FILE_PREFIX :" << bssn::BSSN_CHKPT_FILE_PREFIX
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_PROFILE_FILE_PREFIX :" << bssn::BSSN_PROFILE_FILE_PREFIX
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_VTU_Z_SLICE_ONLY :" << bssn::BSSN_VTU_Z_SLICE_ONLY
             << NRM << std::endl;
        sout << YLW << "\tBSSN_IO_OUTPUT_GAP :" << bssn::BSSN_IO_OUTPUT_GAP
             << NRM << std::endl;
        sout << YLW << "\tBSSN_DENDRO_GRAIN_SZ :" << bssn::BSSN_DENDRO_GRAIN_SZ
             << NRM << std::endl;
        sout << YLW << "\tBSSN_ASYNC_COMM_K :" << bssn::BSSN_ASYNC_COMM_K << NRM
             << std::endl;
        sout << YLW << "\tBSSN_DENDRO_AMR_FAC :" << bssn::BSSN_DENDRO_AMR_FAC
             << NRM << std::endl;
        sout << YLW << "\tBSSN_USE_WAVELET_TOL_FUNCTION :"
             << bssn::BSSN_USE_WAVELET_TOL_FUNCTION << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL :" << bssn::BSSN_WAVELET_TOL << NRM
             << std::endl;
        sout << YLW << "\tBSSN_GW_REFINE_WTOL:" << bssn::BSSN_GW_REFINE_WTOL
             << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL_MAX:" << bssn::BSSN_WAVELET_TOL_MAX
             << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL_FUNCTION_R0: "
             << bssn::BSSN_WAVELET_TOL_FUNCTION_R0 << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL_FUNCTION_R1: "
             << bssn::BSSN_WAVELET_TOL_FUNCTION_R1 << NRM << std::endl;
        sout << YLW << "\tBSSN_LOAD_IMB_TOL :" << bssn::BSSN_LOAD_IMB_TOL << NRM
             << std::endl;
        sout << YLW << "\tBSSN_RK_TIME_BEGIN :" << bssn::BSSN_RK_TIME_BEGIN
             << NRM << std::endl;
        sout << YLW << "\tBSSN_RK_TIME_END :" << bssn::BSSN_RK_TIME_END << NRM
             << std::endl;
        sout << YLW << "\tBSSN_RK_TYPE :" << bssn::BSSN_RK_TYPE << NRM
             << std::endl;
        sout << YLW
             << "\tBSSN_RK45_TIME_STEP_SIZE :" << bssn::BSSN_RK45_TIME_STEP_SIZE
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_RK45_DESIRED_TOL :" << bssn::BSSN_RK45_DESIRED_TOL
             << NRM << std::endl;
        sout << YLW << "\tBSSN_COMPD_MIN : ( " << bssn::BSSN_COMPD_MIN[0]
             << " ," << bssn::BSSN_COMPD_MIN[1] << ","
             << bssn::BSSN_COMPD_MIN[2] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_COMPD_MAX : ( " << bssn::BSSN_COMPD_MAX[0]
             << " ," << bssn::BSSN_COMPD_MAX[1] << ","
             << bssn::BSSN_COMPD_MAX[2] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_BLK_MIN : ( " << bssn::BSSN_BLK_MIN_X << " ,"
             << bssn::BSSN_BLK_MIN_Y << "," << bssn::BSSN_BLK_MIN_Z << " )"
             << NRM << std::endl;
        sout << YLW << "\tBSSN_BLK_MAX : ( " << bssn::BSSN_BLK_MAX_X << " ,"
             << bssn::BSSN_BLK_MAX_Y << "," << bssn::BSSN_BLK_MAX_Z << " )"
             << NRM << std::endl;
        sout << YLW << "\tBSSN_OCTREE_MIN : ( " << bssn::BSSN_OCTREE_MIN[0]
             << " ," << bssn::BSSN_OCTREE_MIN[1] << ","
             << bssn::BSSN_OCTREE_MIN[2] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_OCTREE_MAX : ( " << bssn::BSSN_OCTREE_MAX[0]
             << " ," << bssn::BSSN_OCTREE_MAX[1] << ","
             << bssn::BSSN_OCTREE_MAX[2] << " )" << NRM << std::endl;
        sout << YLW << "\tETA_CONST :" << bssn::ETA_CONST << NRM << std::endl;
        sout << YLW << "\tETA_R0 :" << bssn::ETA_R0 << NRM << std::endl;
        sout << YLW << "\tETA_DAMPING :" << bssn::ETA_DAMPING << NRM
             << std::endl;
        sout << YLW << "\tETA_DAMPING_EXP :" << bssn::ETA_DAMPING_EXP << NRM
             << std::endl;
        sout << YLW << "\tBSSN_ETA_R0 :" << bssn::BSSN_ETA_R0 << NRM
             << std::endl;
        sout << YLW << "\tBSSN_ETA_POWER : (" << bssn::BSSN_ETA_POWER[0] << " ,"
             << bssn::BSSN_ETA_POWER[1] << " )" << NRM << std::endl;
        sout << YLW << "\tRIT_ETA_FUNCTION: " << bssn::RIT_ETA_FUNCTION << NRM
             << std::endl;
        sout << YLW << "\tRIT_ETA_CENTRAL: " << bssn::RIT_ETA_CENTRAL << NRM
             << std::endl;
        sout << YLW << "\tRIT_ETA_WIDTH: " << bssn::RIT_ETA_WIDTH << NRM
             << std::endl;
        sout << YLW << "\tRIT_ETA_OUTER: " << bssn::RIT_ETA_OUTER << NRM
             << std::endl;
        sout << YLW << "\tBSSN_LAMBDA : (" << bssn::BSSN_LAMBDA[0] << " ,"
             << bssn::BSSN_LAMBDA[1] << "," << bssn::BSSN_LAMBDA[2]
             << bssn::BSSN_LAMBDA[3] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_LAMBDA_F : (" << bssn::BSSN_LAMBDA_F[0] << " ,"
             << bssn::BSSN_LAMBDA_F[1] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_XI : (" << bssn::BSSN_XI[0] << " ,"
             << bssn::BSSN_XI[1] << " ," << bssn::BSSN_XI[2] << " )" << NRM
             << std::endl;
        sout << YLW << "\tCHI_FLOOR :" << bssn::CHI_FLOOR << NRM << std::endl;
        sout << YLW << "\tBSSN_TRK0 :" << bssn::BSSN_TRK0 << NRM << std::endl;
        sout << YLW << "\tDISSIPATION_TYPE :" << bssn::DISSIPATION_TYPE << NRM
             << std::endl;
        sout << YLW << "\tKO_DISS_SIGMA :" << bssn::KO_DISS_SIGMA << NRM
             << std::endl;

        sout << YLW << "\tBH1 MASS :" << bssn::BH1.getBHMass() << NRM
             << std::endl;
        sout << YLW << "\tBH1 POSITION (x,y,z) : (" << bssn::BH1.getBHCoordX()
             << ", " << bssn::BH1.getBHCoordY() << ", "
             << bssn::BH1.getBHCoordZ() << " )" << NRM << std::endl;
        sout << YLW << "\tBH1 VELOCITY (x,y,z) : (" << bssn::BH1.getVx() << ", "
             << bssn::BH1.getVy() << ", " << bssn::BH1.getVz() << " )" << NRM
             << std::endl;
        sout << YLW << "\tBH1 SPIN (||,theta,phi): ( " << bssn::BH1.getBHSpin()
             << ", " << bssn::BH1.getBHSpinTheta() << ", "
             << bssn::BH1.getBHSpinPhi() << " )" << NRM << std::endl;

        sout << YLW << "\tBH2 MASS :" << bssn::BH2.getBHMass() << NRM
             << std::endl;
        sout << YLW << "\tBH2 POSITION (x,y,z) : (" << bssn::BH2.getBHCoordX()
             << ", " << bssn::BH2.getBHCoordY() << ", "
             << bssn::BH2.getBHCoordZ() << " )" << NRM << std::endl;
        sout << YLW << "\tBH2 VELOCITY (x,y,z) : (" << bssn::BH2.getVx() << ", "
             << bssn::BH2.getVy() << ", " << bssn::BH2.getVz() << " )" << NRM
             << std::endl;
        sout << YLW << "\tBH2 SPIN (||,theta,phi): ( " << bssn::BH2.getBHSpin()
             << ", " << bssn::BH2.getBHSpinTheta() << ", "
             << bssn::BH2.getBHSpinPhi() << " )" << NRM << std::endl;

        sout << YLW << "\tBSSN_DIM :" << bssn::BSSN_DIM << NRM << std::endl;
        sout << YLW << "\tBSSN_MAXDEPTH :" << bssn::BSSN_MAXDEPTH << NRM
             << std::endl;
        sout << YLW << "\tBSSN_MINDEPTH :" << bssn::BSSN_MINDEPTH << NRM
             << std::endl;

        sout << YLW << "\tBSSN_NUM_REFINE_VARS :" << bssn::BSSN_NUM_REFINE_VARS
             << NRM << std::endl;
        sout << YLW << "\tBSSN_REFINE_VARIABLE_INDICES :[";
        for (unsigned int i = 0; i < bssn::BSSN_NUM_REFINE_VARS - 1; i++)
            sout << bssn::BSSN_REFINE_VARIABLE_INDICES[i] << ", ";
        sout << bssn::BSSN_REFINE_VARIABLE_INDICES[bssn::BSSN_NUM_REFINE_VARS -
                                                   1]
             << "]" << NRM << std::endl;

        sout << YLW << "\tBSSN_REFINEMENT_MODE :" << bssn::BSSN_REFINEMENT_MODE
             << NRM << std::endl;
        sout << YLW << "\tBSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE :"
             << bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE << NRM
             << std::endl;

        sout << YLW << "\tBSSN_BH1_AMR_R: " << bssn::BSSN_BH1_AMR_R << NRM
             << std::endl;
        sout << YLW << "\tBSSN_BH2_AMR_R: " << bssn::BSSN_BH2_AMR_R << NRM
             << std::endl;
        sout << YLW << "\tBSSN_AMR_R_RATIO: " << bssn::BSSN_AMR_R_RATIO << NRM
             << std::endl;

        sout << YLW << "\tBSSN_BH1_CONSTRAINT_R:" << bssn::BSSN_BH1_CONSTRAINT_R
             << NRM << std::endl;
        sout << YLW << "\tBSSN_BH2_CONSTRAINT_R:" << bssn::BSSN_BH2_CONSTRAINT_R
             << NRM << std::endl;

        sout << YLW << "\tBSSN_BH1_MAX_LEV:" << bssn::BSSN_BH1_MAX_LEV << NRM
             << std::endl;
        sout << YLW << "\tBSSN_BH2_MAX_LEV:" << bssn::BSSN_BH2_MAX_LEV << NRM
             << std::endl;
        sout << YLW << "\tBSSN_INIT_GRID_ITER:" << bssn::BSSN_INIT_GRID_ITER
             << NRM << std::endl;

#ifdef BSSN_REFINE_BASE_EH
        sout << YLW << "\tBSSN_EH_REFINE_VAL  : " << bssn::BSSN_EH_REFINE_VAL
             << NRM << std::endl;
        sout << YLW << "\tBSSN_EH_COARSEN_VAL : " << bssn::BSSN_EH_COARSEN_VAL
             << NRM << std::endl;
#endif

        sout << YLW << "\tBSSN_NUM_EVOL_VARS_VTU_OUTPUT :"
             << bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT << NRM << std::endl;
        sout << YLW << "\tBSSN_VTU_OUTPUT_EVOL_INDICES :[";
        for (unsigned int i = 0; i < bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT - 1;
             i++)
            sout << bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[i] << ", ";
        sout << bssn::BSSN_VTU_OUTPUT_EVOL_INDICES
                    [bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT - 1]
             << "]" << NRM << std::endl;

        sout << YLW << "\tBSSN_NUM_CONST_VARS_VTU_OUTPUT :"
             << bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT << NRM << std::endl;
        sout << YLW << "\tBSSN_VTU_OUTPUT_CONST_INDICES :[";
        for (unsigned int i = 0; i < bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT - 1;
             i++)
            sout << bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i] << ", ";
        sout << bssn::BSSN_VTU_OUTPUT_CONST_INDICES
                    [bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT - 1]
             << "]" << NRM << std::endl;

        sout << YLW << "\tTPID_TARGET_M_PLUS :" << TPID::target_M_plus << NRM
             << std::endl;
        sout << YLW << "\tTPID_TARGET_M_MINUS :" << TPID::target_M_minus << NRM
             << std::endl;
        sout << YLW << "\tTPID_PAR_B :" << TPID::par_b << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_P_PLUS : ( " << TPID::par_P_plus[0] << ", "
             << TPID::par_P_plus[1] << ", " << TPID::par_P_plus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_P_MINUS : ( " << TPID::par_P_minus[0] << ", "
             << TPID::par_P_minus[1] << ", " << TPID::par_P_minus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_S_PLUS : ( " << TPID::par_S_plus[0] << ", "
             << TPID::par_S_plus[1] << ", " << TPID::par_S_plus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_S_MINUS : ( " << TPID::par_S_minus[0] << ", "
             << TPID::par_S_minus[1] << ", " << TPID::par_S_minus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_CENTER_OFFSET : ( " << TPID::center_offset[0]
             << ", " << TPID::center_offset[1] << ", " << TPID::center_offset[2]
             << " )" << NRM << std::endl;

        sout << YLW << "\tTPID_INITIAL_LAPSE_PSI_EXPONENT :"
             << TPID::initial_lapse_psi_exponent << NRM << std::endl;
        sout << YLW << "\tTPID_NPOINTS_A :" << TPID::npoints_A << NRM
             << std::endl;
        sout << YLW << "\tTPID_NPOINTS_B :" << TPID::npoints_B << NRM
             << std::endl;
        sout << YLW << "\tTPID_NPOINTS_PHI :" << TPID::npoints_phi << NRM
             << std::endl;
        sout << YLW << "\tTPID_GIVE_BARE_MASS :" << TPID::give_bare_mass << NRM
             << std::endl;
        sout << YLW << "\tINITIAL_LAPSE :" << TPID::initial_lapse << NRM
             << std::endl;
        sout << YLW << "\tTPID_SOLVE_MOMENTUM_CONSTRAINT :"
             << TPID::solve_momentum_constraint << NRM << std::endl;
        sout << YLW << "\tTPID_GRID_SETUP_METHOD :" << TPID::grid_setup_method
             << NRM << std::endl;
        sout << YLW << "\tTPID_VERBOSE :" << TPID::verbose << NRM << std::endl;
        sout << YLW << "\tTPID_ADM_TOL :" << TPID::adm_tol << NRM << std::endl;
        sout << YLW << "\tTPID_NEWTON_TOL :" << TPID::Newton_tol << NRM
             << std::endl;

        sout << YLW << "\tEXTRACTION_VAR_ID :" << BHLOC::EXTRACTION_VAR_ID
             << NRM << std::endl;
        sout << YLW << "\tEXTRACTION_TOL :" << BHLOC::EXTRACTION_TOL << NRM
             << std::endl;

        sout << YLW << "\tBSSN_GW_NUM_RADAII: " << GW::BSSN_GW_NUM_RADAII << NRM
             << std::endl;
        sout << YLW << "\tBSSN_GW_NUM_LMODES: " << GW::BSSN_GW_NUM_LMODES << NRM
             << std::endl;

        sout << YLW << "\tBSSN_GW_RADAII: {";
        for (unsigned int i = 0; i < GW::BSSN_GW_NUM_RADAII; i++)
            sout << " ," << GW::BSSN_GW_RADAII[i];
        sout << "}" << NRM << std::endl;

        sout << YLW << "\tBSSN_GW_L_MODES: {";
        for (unsigned int i = 0; i < GW::BSSN_GW_NUM_LMODES; i++)
            sout << " ," << GW::BSSN_GW_L_MODES[i];
        sout << "}" << NRM << std::endl;

        sout << YLW << "\tAEH_SOLVER_FREQ: " << AEH::AEH_SOLVER_FREQ
             << std::endl;
        sout << YLW << "\tAEH_LMAX: " << AEH::AEH_LMAX << std::endl;
        sout << YLW << "\tAEH_Q_THETA: " << AEH::AEH_Q_THETA << std::endl;
        sout << YLW << "\tAEH_Q_PHI: " << AEH::AEH_Q_PHI << std::endl;
        sout << YLW << "\tAEH_MAXITER: " << AEH::AEH_MAXITER << std::endl;
        sout << YLW << "\tAEH_ATOL: " << AEH::AEH_ATOL << std::endl;
        sout << YLW << "\tAEH_RTOL: " << AEH::AEH_RTOL << NRM << std::endl;

        sout << YLW << "\tAEH_ALPHA: " << AEH::AEH_ALPHA << std::endl;
        sout << YLW << "\tAEH_BETA: " << AEH::AEH_BETA << NRM << std::endl;

        sout << YLW << "\tBSSN_KO_SIGMA_SCALE_BY_CONFORMAL: "
             << (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL ? "true" : "false")
             << NRM << std::endl;
        if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL) {
            sout << YLW << "\t\tBSSN_PSILON_CAKO_GAUGE: "
                 << bssn::BSSN_EPSILON_CAKO_GAUGE << NRM << std::endl;
            sout << YLW << "\t\tBSSN_PSILON_CAKO_OTHER: "
                 << bssn::BSSN_EPSILON_CAKO_OTHER << NRM << std::endl;
        }
        sout << YLW << "\tBSSN_NYQUIST_M: " << bssn::BSSN_NYQUIST_M << NRM
             << std::endl;

        sout << YLW << "\tBSSN_SSL_H: " << bssn::BSSN_SSL_H << NRM << std::endl;
        sout << YLW << "\tBSSN_SSL_SIGMA: " << bssn::BSSN_SSL_SIGMA << NRM
             << std::endl;
    }
}

double CalTolHelper(const double t, const double r, const double rad[],
                    const double eps[], const double toffset) {
    const double R0     = rad[0];
    const double R1     = rad[1];
    const double RGW    = rad[2];
    const double tol    = eps[0];
    const double tolGW  = eps[1];
    const double tolMax = eps[2];
    const double WRR    = std::min(tolGW, tolMax);
    if (r < R0) {
        return tol;
    }
    if (t > (R1 + toffset)) {
        const double RR         = std::min(t - toffset, RGW + 10.0);
        const double WTolExpFac = (RR - R0) / log10(WRR / tol);
        return std::min(tolMax, tol * pow(10.0, ((r - R0) / WTolExpFac)));
    } else {
        const double WTolExpFac = (R1 - R0) / log10(WRR / tol);
        return std::min(tolMax, tol * pow(10.0, ((r - R0) / WTolExpFac)));
    }
}

void initialDataFunctionWrapper(const double xx_grid, const double yy_grid,
                                const double zz_grid, double* var) {
    // code to convert grid functions to non-grid, if the function doesn't
    // already support it const double xx = GRIDX_TO_X(xx_grid); const double yy
    // = GRIDY_TO_Y(yy_grid); const double zz = GRIDZ_TO_Z(zz_grid);

    switch (bssn::BSSN_ID_TYPE) {
        case 0:
            // NOTE: this is the TwoPunctures code! For **pure**
            // initialization done when building up the grid before pure
            // population, this calls bssn::punctureData which is simply the
            // two puncture intiial data from HAD code. In bssnCtx.cpp's
            // init_grid the bssn::BSSN_ID_TYPE switch will call the true
            // TwoPunctures code which is the proper TPID initial data
            bssn::punctureData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 1:
            // ID 1 is the puncture data function on its own
            bssn::punctureData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 2:
            // ID 2 is the KerrSchild Data
            bssn::KerrSchildData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 3:
            // ID 3 is a noise data
            bssn::noiseData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 4:
            // ID 4 is a fake initial data
            bssn::fake_initial_data(xx_grid, yy_grid, zz_grid, var);

            break;

        case 5:
        // minkowski initial data is flat space!
        bssn:
            minkowskiInitialData(xx_grid, yy_grid, zz_grid, var);

            break;
        case 6:
            // minkowski initial data is flat space!
            bssn::kerrData(xx_grid, yy_grid, zz_grid, var);

            break;
            // MORE CAN BE ADDED HERE

        default:
            int rank;
            int npes;
            MPI_Comm comm = MPI_COMM_WORLD;
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &npes);
            if (!rank) {
                std::cerr << RED << "ERROR::: Invalid initial data ID: "
                          << bssn::BSSN_ID_TYPE << NRM << std::endl;
            }

            MPI_Abort(comm, 0);

            break;
    }
}

void punctureDataPhysicalCoord(const double xx, const double yy,
                               const double zz, double* var) {
    /* Define the Levi-Cevita pseudo-tensor and Kroneckar delta */
    double epijk[3][3][3];
    int i, j, k;
    for (k = 0; k < 3; k++) {
        for (j = 0; j < 3; j++) {
            for (i = 0; i < 3; i++) {
                epijk[k][j][i] = 0.0;
            }
        }
    }
    epijk[0][1][2] = 1.0;
    epijk[1][2][0] = 1.0;
    epijk[2][0][1] = 1.0;
    epijk[0][2][1] = -1.0;
    epijk[2][1][0] = -1.0;
    epijk[1][0][2] = -1.0;

    double deltaij[3][3];
    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) {
            deltaij[j][i] = 0.0;
        }
    }

    deltaij[0][0] = 1.0;
    deltaij[1][1] = 1.0;
    deltaij[2][2] = 1.0;

    double x1, y1, z1, rv1;
    double x2, y2, z2, rv2;
    double vn1[3], vn2[3];

    double vpsibl;
    double v_u_corr, amp_capj, amp_capr, l_r, u0_j, u2_j, mu_j, p2_mu_j, v_u_j1;
    double v1, v2, v3, v4, vt1, vt2;

    int i1, i2, i3, i4;
    double amp_capp, u0_p, u2_p, mu_p, p2_mu_p;
    double v_u_p1, v_u_c1, v_u_j2, v_u_p2;
    double v_u_c2, vpsibl_u, vpsibl_u2;

    // bh 1
    double mass1 = BH1.getBHMass();
    double bh1x  = BH1.getBHCoordX();
    double bh1y  = BH1.getBHCoordY();
    double bh1z  = BH1.getBHCoordZ();

    double vp1[3];
    vp1[0]          = BH1.getVx();
    vp1[1]          = BH1.getVy();
    vp1[2]          = BH1.getVz();

    double vp1tot   = sqrt(vp1[0] * vp1[0] + vp1[1] * vp1[1] + vp1[2] * vp1[2]);
    double spin1    = BH1.getBHSpin();
    double spin1_th = BH1.getBHSpinTheta();
    double spin1_phi = BH1.getBHSpinPhi();
    double vs1[3];

    vs1[0]       = spin1 * sin(spin1_th) * cos(spin1_phi);
    vs1[1]       = spin1 * sin(spin1_th) * sin(spin1_phi);
    vs1[2]       = spin1 * cos(spin1_th);

    // bh 2
    double mass2 = BH2.getBHMass();
    double bh2x  = BH2.getBHCoordX();
    double bh2y  = BH2.getBHCoordY();
    double bh2z  = BH2.getBHCoordZ();

    double vp2[3];
    vp2[0]          = BH2.getVx();
    vp2[1]          = BH2.getVy();
    vp2[2]          = BH2.getVz();

    double vp2tot   = sqrt(vp2[0] * vp2[0] + vp2[1] * vp2[1] + vp2[2] * vp2[2]);
    double spin2    = BH2.getBHSpin();
    double spin2_th = BH2.getBHSpinTheta();
    double spin2_phi = BH2.getBHSpinPhi();

    double vs2[3];
    vs2[0]   = spin2 * sin(spin2_th) * cos(spin2_phi);
    vs2[1]   = spin2 * sin(spin2_th) * sin(spin2_phi);
    vs2[2]   = spin2 * cos(spin2_th);

    // coordinates with respect to center of bh1
    x1       = xx - bh1x;
    y1       = yy - bh1y;
    z1       = zz - bh1z;

    // locating as a radial form
    rv1      = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    vn1[0]   = x1 / rv1;
    vn1[1]   = y1 / rv1;
    vn1[2]   = z1 / rv1;

    // same as BH2
    x2       = xx - bh2x;
    y2       = yy - bh2y;
    z2       = zz - bh2z;

    rv2      = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    vn2[0]   = x2 / rv2;
    vn2[1]   = y2 / rv2;
    vn2[2]   = z2 / rv2;

    // Initial data is related with the paper: http://arxiv.org/abs/0711.1165
    // Brill-Lindquist conformal factor
    vpsibl   = 1.0 + mass1 / (2.0 * rv1);
    vpsibl   = vpsibl + mass2 / (2.0 * rv2);

    v_u_corr = 0.0;
    // bh 1

    // For spinning puncture
    if (fabs(spin1) > 1.e-6) {
        amp_capj = 4.0 * spin1 / (mass1 * mass1);
        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_j =
            (l_r + l_r * l_r + l_r * l_r * l_r - 4.0 * l_r * l_r * l_r * l_r +
             2.0 * l_r * l_r * l_r * l_r * l_r) /
            40.0;
        u2_j    = -pow(l_r, 5) / 20.0;
        mu_j    = vn1[0] * vs1[0];
        mu_j    = mu_j + vn1[1] * vs1[1];
        mu_j    = (mu_j + vn1[2] * vs1[2]) / fabs(spin1);
        p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;
        v_u_j1 =
            amp_capj * amp_capj * (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);
        v_u_corr = v_u_corr + v_u_j1;
    }
    // For boosting puncture
    if (vp1tot > 1.e-6) {
        amp_capp = 2.0 * vp1tot / mass1;
        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_p     = l_r - 2.0 * l_r * l_r + 2.0 * pow(l_r, 3);
        u0_p     = (u0_p - pow(l_r, 4) + 0.20 * pow(l_r, 5)) * (5.0 / 32.0);
        u2_p     = 15.0 * l_r + 132.0 * l_r * l_r + 53.0 * pow(l_r, 3);
        u2_p     = u2_p + 96.0 * pow(l_r, 4) + 82.0 * pow(l_r, 5);
        u2_p = u2_p + (84.0 / amp_capr) * (pow(l_r, 5) + log(l_r) / amp_capr);
        u2_p = (u2_p) / (80.0 * amp_capr);
        mu_p = vn1[0] * vp1[0] / vp1tot;
        mu_p = mu_p + vn1[1] * vp1[1] / vp1tot;
        mu_p = mu_p + vn1[2] * vp1[2] / vp1tot;
        p2_mu_p  = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;
        v_u_p1   = pow(amp_capp, 2) * (u0_p + u2_p * p2_mu_p);
        v_u_corr = v_u_corr + v_u_p1;
    }
    // For spinning boosted pucture
    if (vp1tot > 1.e-6 && fabs(spin1) > 1.e-6) {
        v1       = (vp1[1] * vs1[2] - vp1[2] * vs1[1]) * vn1[0];
        v1       = v1 + (vp1[2] * vs1[0] - vp1[0] * vs1[2]) * vn1[1];
        v1       = v1 + (vp1[0] * vs1[1] - vp1[1] * vs1[0]) * vn1[2];
        v1       = v1 * (16.0 / pow(mass1, 4)) * rv1;

        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);

        v2       = 1.0 + 5.0 * amp_capr + 10.0 * pow(amp_capr, 2);

        v_u_c1   = (v1 * v2 * pow(l_r, 5)) / 80.0;
        v_u_corr = v_u_corr + v_u_c1;
    }
    // bh 2 same puncture as bh 1
    if (fabs(spin2) > 1.e-6) {
        amp_capj = 4.0 * spin2 / (mass2 * mass2);
        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_j =
            (l_r + l_r * l_r + l_r * l_r * l_r - 4.0 * l_r * l_r * l_r * l_r +
             2.0 * l_r * l_r * l_r * l_r * l_r) /
            40.0;
        u2_j    = -pow(l_r, 5) / 20.0;
        mu_j    = vn2[0] * vs2[0];
        mu_j    = mu_j + vn2[1] * vs2[1];
        mu_j    = (mu_j + vn2[2] * vs2[2]) / fabs(spin2);
        p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;
        v_u_j2 =
            amp_capj * amp_capj * (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);
        v_u_corr = v_u_corr + v_u_j2;
    }

    if (vp2tot > 1.e-6) {
        amp_capp = 2.0 * vp2tot / mass2;
        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);
        u0_p     = l_r - 2.0 * l_r * l_r + 2.0 * pow(l_r, 3);
        u0_p     = (u0_p - pow(l_r, 4) + 0.20 * pow(l_r, 5)) * (5.0 / 32.0);
        u2_p     = 15.0 * l_r + 132.0 * l_r * l_r + 53.0 * pow(l_r, 3);
        u2_p     = u2_p + 96.0 * pow(l_r, 4) + 82.0 * pow(l_r, 5);
        u2_p = u2_p + (84.0 / amp_capr) * (pow(l_r, 5) + log(l_r) / amp_capr);
        u2_p = (u2_p) / (80.0 * amp_capr);
        mu_p = vn2[0] * vp2[0] / vp2tot;
        mu_p = mu_p + vn2[1] * vp2[1] / vp2tot;
        mu_p = mu_p + vn2[2] * vp2[2] / vp2tot;
        p2_mu_p  = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;
        v_u_p2   = pow(amp_capp, 2) * (u0_p + u2_p * p2_mu_p);
        v_u_corr = v_u_corr + v_u_p2;
    }

    if (vp2tot > 1.e-6 && fabs(spin2) > 1.e-6) {
        v1       = (vp2[1] * vs2[2] - vp2[2] * vs2[1]) * vn2[0];
        v1       = v1 + (vp2[2] * vs2[0] - vp2[0] * vs2[2]) * vn2[1];
        v1       = v1 + (vp2[0] * vs2[1] - vp2[1] * vs2[0]) * vn2[2];
        v1       = v1 * (16.0 / pow(mass2, 4)) * rv2;

        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);

        v2       = 1.0 + 5.0 * amp_capr + 10.0 * pow(amp_capr, 2);

        v_u_c2   = (v1 * v2 * pow(l_r, 5)) / 80.0;
        v_u_corr = v_u_corr + v_u_c2;
    }

    // vpsibl_u will be used for the conformal factor,
    vpsibl_u          = vpsibl + v_u_corr;
    // vpsibl_u2 is for the Aij terms...
    // ! since the corrections are first order...
    // ! adding half of the correction seems to give the best results...
    // ! update - do a fit for spin = 0.6...
    vpsibl_u2         = vpsibl + v_u_corr;

    var[VAR::U_ALPHA] = 1.0 / (vpsibl_u * vpsibl_u);
    // std::cout<<"Alpha: "<<u[U_ALPHA]<<" vpsibl_u: "<< vpsibl_u<<std::endl;
    var[VAR::U_ALPHA] = std::max(var[VAR::U_ALPHA], CHI_FLOOR);

    v2                = 1.0 / pow(vpsibl_u, 4);
    var[VAR::U_CHI]   = v2;

    if (var[VAR::U_CHI] < CHI_FLOOR) var[VAR::U_CHI] = CHI_FLOOR;

    var[VAR::U_K]      = 0.0;

    var[VAR::U_BETA0]  = 0.0;
    var[VAR::U_BETA1]  = 0.0;
    var[VAR::U_BETA2]  = 0.0;

    var[VAR::U_GT0]    = 0.0;
    var[VAR::U_GT1]    = 0.0;
    var[VAR::U_GT2]    = 0.0;

    var[VAR::U_B0]     = 0.0;
    var[VAR::U_B1]     = 0.0;
    var[VAR::U_B2]     = 0.0;

    var[VAR::U_SYMGT0] = 1.0;  // XX
    var[VAR::U_SYMGT1] = 0.0;  // XY
    var[VAR::U_SYMGT2] = 0.0;  // XZ
    var[VAR::U_SYMGT3] = 1.0;  // YY
    var[VAR::U_SYMGT4] = 0.0;  // YZ
    var[VAR::U_SYMGT5] = 1.0;  // ZZ

    for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++) {
            // first BH
            v2 = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                for (i4 = 0; i4 < 3; i4++) {
                    vt1 = epijk[i1][i3][i4] * vs1[i3] * vn1[i4] * vn1[i2];
                    vt2 = epijk[i2][i3][i4] * vs1[i3] * vn1[i4] * vn1[i1];
                    v2  = v2 + vt1 + vt2;
                }
            }

            v3  = vp1[i1] * vn1[i2] + vp1[i2] * vn1[i1];
            vt1 = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                vt1 = vt1 + vp1[i3] * vn1[i3];
            }
            vt1 = vt1 * (vn1[i1] * vn1[i2] - deltaij[i1][i2]);
            v3  = v3 + vt1;

            v1  = 3.0 / (pow(vpsibl_u2, 6) * pow(rv1, 3));
            v4  = v1 * (v2 + (rv1 / 2.0) * v3);

            // second BH
            v2  = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                for (i4 = 0; i4 < 3; i4++) {
                    vt1 = epijk[i1][i3][i4] * vs2[i3] * vn2[i4] * vn2[i2];
                    vt2 = epijk[i2][i3][i4] * vs2[i3] * vn2[i4] * vn2[i1];
                    v2  = v2 + vt1 + vt2;
                }
            }

            v3  = vp2[i1] * vn2[i2] + vp2[i2] * vn2[i1];
            vt1 = 0.0;
            for (i3 = 0; i3 < 3; i3++) {
                vt1 = vt1 + vp2[i3] * vn2[i3];
            }
            vt1 = vt1 * (vn2[i1] * vn2[i2] - deltaij[i1][i2]);
            v3  = v3 + vt1;

            v1  = 3.0 / (pow(vpsibl_u2, 6) * pow(rv2, 3));
            v4  = v4 + v1 * (v2 + (rv2 / 2.0) * v3);

            if (i1 == 0 && i2 == 0) {
                var[VAR::U_SYMAT0] = v4;  // XX
            } else if (i1 == 0 && i2 == 1) {
                var[VAR::U_SYMAT1] = v4;  // XY
            } else if (i1 == 0 && i2 == 2) {
                var[VAR::U_SYMAT2] = v4;  // XZ
            } else if (i1 == 1 && i2 == 1) {
                var[VAR::U_SYMAT3] = v4;  // YY
            } else if (i1 == 1 && i2 == 2) {
                var[VAR::U_SYMAT4] = v4;  // YZ
            } else if (i1 == 2 && i2 == 2) {
                var[VAR::U_SYMAT5] = v4;  // ZZ
            }
        }
    }
}

void punctureData(const double xx1, const double yy1, const double zz1,
                  double* var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    punctureDataPhysicalCoord(xx, yy, zz, var);
}
void kerrData(const double xx1, const double yy1, const double zz1,
              double* var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    // parameters for the BH (mass, location, spin parameter)
    double M        = BH1.getBHMass();
    double bh1x     = BH1.getBHCoordX();
    double bh1y     = BH1.getBHCoordY();
    double bh1z     = BH1.getBHCoordZ();
    double spin1    = BH1.getBHSpin();

    // coordinates relative to the center of the BH
    double x        = xx - bh1x;
    double y        = yy - bh1y;
    double z        = zz - bh1z;

    // locating as a radial form
    double r        = sqrt(x * x + y * y + z * z);

    // HL : Angular momentum parameter will be added as param file after
    // testing
    double a        = spin1;

    double gtd[3][3], Atd[3][3];
    double alpha, Gamt[3];
    double Chi, TrK, Betau[3];

#include "Kerr.cpp"
#include "kerr_vars.cpp"
}

void KerrSchildData(const double xx1, const double yy1, const double zz1,
                    double* var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    // parameters for the BH (mass, location, spin parameter)
    double M        = BH1.getBHMass();
    double bh1x     = BH1.getBHCoordX();
    double bh1y     = BH1.getBHCoordY();
    double bh1z     = BH1.getBHCoordZ();
    double spin1    = BH1.getBHSpin();

    // coordinates relative to the center of the BH
    double x        = xx - bh1x;
    double y        = yy - bh1y;
    double z        = zz - bh1z;

    // locating as a radial form
    double r        = sqrt(x * x + y * y + z * z);

    // HL : Angular momentum parameter will be added as param file after testing
    double a        = spin1;

    double gtd[3][3], Atd[3][3];
    double alpha, Gamt[3];
    double Chi, TrK, Betau[3];

#include "ks_vars.cpp"
#include "ksinit.cpp"

    var[VAR::U_ALPHA]  = alpha;
    var[VAR::U_CHI]    = Chi;
    var[VAR::U_K]      = TrK;

    var[VAR::U_BETA0]  = Betau[0];
    var[VAR::U_BETA1]  = Betau[1];
    var[VAR::U_BETA2]  = Betau[2];

    var[VAR::U_GT0]    = Gamt[0];
    var[VAR::U_GT1]    = Gamt[1];
    var[VAR::U_GT2]    = Gamt[2];

    var[VAR::U_B0]     = 0.0;
    var[VAR::U_B1]     = 0.0;
    var[VAR::U_B2]     = 0.0;

    var[VAR::U_SYMGT0] = gtd[0][0];
    var[VAR::U_SYMGT1] = gtd[0][1];
    var[VAR::U_SYMGT2] = gtd[0][2];
    var[VAR::U_SYMGT3] = gtd[1][1];
    var[VAR::U_SYMGT4] = gtd[1][2];
    var[VAR::U_SYMGT5] = gtd[2][2];

    var[VAR::U_SYMAT0] = Atd[0][0];
    var[VAR::U_SYMAT1] = Atd[0][1];
    var[VAR::U_SYMAT2] = Atd[0][2];
    var[VAR::U_SYMAT3] = Atd[1][1];
    var[VAR::U_SYMAT4] = Atd[1][2];
    var[VAR::U_SYMAT5] = Atd[2][2];

    // std::cout<<"KS init data: (x,y,z) = ( "<<x<<", "<<y<<", "<<z<<"), alpha =
    // "<<alpha<<std::endl;

#if 0
            //BSSN vars for Kerr-Schild
            var[VAR::U_ALPHA] = sqrt(rv1/(2.0*M+rv1));
            var[VAR::U_CHI] = 1.0/pow(1.0+2.0*M/rv1, 1.0/3.0);
            var[VAR::U_K] = 2.0*M*sqrt(rv1/(2.0*M+rv1))*(rv1+3.0*M)/(rv1*rv1*(2.0*M+rv1));

            var[VAR::U_BETA0] = 2.0*M*x1/(rv1*(2.0*M+rv1));
            var[VAR::U_BETA1] = 2.0*M*y1/(rv1*(2.0*M+rv1));
            var[VAR::U_BETA2] = 2.0*M*z1/(rv1*(2.0*M+rv1));

            var[VAR::U_GT0] = pow(2,8.0/3.0)*M*x1*(M+rv1*rv1+(3*M*M*rv1+x1*x1+y1*y1)/5.0)/(pow(5.0,1.0/3.0)*rv1*pow(M/5.0+rv1,2.0)*pow(M/(5.0*rv1)+1.0,2.0/3.0));
            var[VAR::U_GT1] = pow(2,8.0/3.0)*M*y1*(M+rv1*rv1+(3*M*M*rv1+x1*x1+y1*y1)/5.0)/(pow(5.0,1.0/3.0)*rv1*pow(M/5.0+rv1,2.0)*pow(M/(5.0*rv1)+1.0,2.0/3.0));
            var[VAR::U_GT2] = pow(2,8.0/3.0)*M*z1*(M+rv1*rv1+(3*M*M*rv1+x1*x1+y1*y1)/5.0)/(pow(5.0,1.0/3.0)*rv1*pow(M/5.0+rv1,2.0)*pow(M/(5.0*rv1)+1.0,2.0/3.0));

            var[VAR::U_B0] = M*x1/(500.0*rv1*(M/5.0+rv1));
            var[VAR::U_B1] = M*y1/(500.0*rv1*(M/5.0+rv1));
            var[VAR::U_B2] = M*z1/(500.0*rv1*(M/5.0+rv1));

            var[VAR::U_SYMGT0] = (M*x1*x1/2.0+rv1*rv1)/pow(10,5.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0); //XX
            var[VAR::U_SYMGT1] = M*x1*y1/(50.0*pow(10.0,2.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0)); //XY
            var[VAR::U_SYMGT2] = M*x1*z1/(50.0*pow(10.0,2.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0)); //XZ
            var[VAR::U_SYMGT3] = (M*y1*y1/2.0+rv1*rv1)/pow(10,5.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0); //YY
            var[VAR::U_SYMGT4] = M*y1*z1/(50.0*pow(10.0,2.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0)); //YZ
            var[VAR::U_SYMGT5] = (M*z1*z1/2.0+rv1*rv1)/pow(10,5.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0); //ZZ

            var[VAR::U_SYMAT0] = (sqrt(rv1/(M/5.0+rv1))*(rv1-(1.0/5.0+M/(10.0*rv1))*x1*x1)-M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0+rv1)*(1.0/10.0+M*x1*x1/(50.0*rv1*rv1))/(150.0*sqrt(10.0)*(M/5.0+rv1)*rv1))/pow(0.1+M/(50.0*rv1),1.0/3.0); //XX
            var[VAR::U_SYMAT1] = -(((1.0/5.0+M/(10.0*rv1))*sqrt(rv1/(M/5.0+rv1))*x1*y1)/(50.0*sqrt(10.0)*M*rv1)+(M*M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0*rv1)*x1*y1)/(750.0*sqrt(10.0)*rv1*rv1*(M/5.0*rv1)))/pow(0.1+M/(50.0*rv1),1.0/3.0); //XY
            var[VAR::U_SYMAT2] = -(((1.0/5.0+M/(10.0*rv1))*sqrt(rv1/(M/5.0+rv1))*x1*y1)/(50.0*sqrt(10.0)*M*rv1)+(M*M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0*rv1)*x1*z1)/(750.0*sqrt(10.0)*rv1*rv1*(M/5.0*rv1)))/pow(0.1+M/(50.0*rv1),1.0/3.0); //XZ
            var[VAR::U_SYMAT3] = (sqrt(rv1/(M/5.0+rv1))*(rv1-(1.0/5.0+M/(10.0*rv1))*y1*y1)-M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0+rv1)*(1.0/10.0+M*x1*x1/(50.0*rv1*rv1))/(150.0*sqrt(10.0)*(M/5.0+rv1)*rv1))/pow(0.1+M/(50.0*rv1),1.0/3.0); //YY
            var[VAR::U_SYMAT4] = -(((1.0/5.0+M/(10.0*rv1))*sqrt(rv1/(M/5.0+rv1))*x1*y1)/(50.0*sqrt(10.0)*M*rv1)+(M*M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0*rv1)*y1*z1)/(750.0*sqrt(10.0)*rv1*rv1*(M/5.0*rv1)))/pow(0.1+M/(50.0*rv1),1.0/3.0); //YZ
            var[VAR::U_SYMAT5] = (sqrt(rv1/(M/5.0+rv1))*(rv1-(1.0/5.0+M/(10.0*rv1))*z1*z1)-M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0+rv1)*(1.0/10.0+M*x1*x1/(50.0*rv1*rv1))/(150.0*sqrt(10.0)*(M/5.0+rv1)*rv1))/pow(0.1+M/(50.0*rv1),1.0/3.0); //ZZ
#endif
}

void noiseData(const double xx1, const double yy1, const double zz1,
               double* var) {
    // const double xx=GRIDX_TO_X(xx1);
    // const double yy=GRIDY_TO_Y(yy1);
    // const double zz=GRIDZ_TO_Z(zz1);

    // call random number generator between -1 and 1)
    double random_variable[30];
    int i;
    for (i = 0; i < 30; i++) {
        random_variable[i] = 2.0 * rand() / ((double)RAND_MAX) - 1.0;
    }

    // set a (uniform) amplitude for the noise
    double noise_amp   = 1.0e-8;

    var[VAR::U_ALPHA]  = 1.0 + noise_amp * random_variable[0];

    var[VAR::U_CHI]    = 1.0 + noise_amp * random_variable[1];

    var[VAR::U_K]      = noise_amp * random_variable[2];

    var[VAR::U_GT0]    = noise_amp * random_variable[3];
    var[VAR::U_GT1]    = noise_amp * random_variable[4];
    var[VAR::U_GT2]    = noise_amp * random_variable[5];

    var[VAR::U_BETA0]  = noise_amp * random_variable[6];
    var[VAR::U_BETA1]  = noise_amp * random_variable[7];
    var[VAR::U_BETA2]  = noise_amp * random_variable[8];

    var[VAR::U_B0]     = noise_amp * random_variable[9];
    var[VAR::U_B1]     = noise_amp * random_variable[10];
    var[VAR::U_B2]     = noise_amp * random_variable[11];

    var[VAR::U_SYMGT0] = 1.0 + noise_amp * random_variable[12];  // XX
    var[VAR::U_SYMGT1] = noise_amp * random_variable[13];        // XY
    var[VAR::U_SYMGT2] = noise_amp * random_variable[14];        // XZ
    var[VAR::U_SYMGT3] = 1.0 + noise_amp * random_variable[15];  // YY
    var[VAR::U_SYMGT4] = noise_amp * random_variable[16];        // YZ
    var[VAR::U_SYMGT5] = 1.0 + noise_amp * random_variable[17];  // ZZ

    var[VAR::U_SYMAT0] = noise_amp * random_variable[18];  // XX
    var[VAR::U_SYMAT1] = noise_amp * random_variable[19];  // XY
    var[VAR::U_SYMAT2] = noise_amp * random_variable[20];  // XZ
    var[VAR::U_SYMAT3] = noise_amp * random_variable[21];  // YY
    var[VAR::U_SYMAT4] = noise_amp * random_variable[22];  // YZ
    var[VAR::U_SYMAT5] = noise_amp * random_variable[23];  // ZZ
}

void fake_initial_data(double xx1, double yy1, double zz1, double* u) {
    /* const double x=GRIDX_TO_X(xx1);
     const double y=GRIDY_TO_Y(yy1);
     const double z=GRIDZ_TO_Z(zz1);


     const double pi = acos(-1.0);
     const double f1 = 31.0/17.0;
     const double f2 = 37.0/11.0;

     u[VAR::U_ALPHA] = 1.0 - 0.25*sin(f1*x);
     //u[F_ALPHA][pp] = 1.0;
     u[VAR::U_BETA0] = 4.0/17.0*sin(x)*cos(z);
     u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
     u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

     u[VAR::U_B0] = 31.0*x*cos(f1*z+y);
     u[VAR::U_B1] = 7.0*y*sin(f1*x+y) + 3.0*cos(z);
     u[VAR::U_B2] = 5.0*z*cos(f1*x+y) + 7.0*sin(z+y+x) + 1.0;

     u[VAR::U_GT0] = 5.0*cos(x)/(10.0*sin(x+z)+26.0-1.0*cos(x*z)*cos(x));
     u[VAR::U_GT1] = -5.0*sin(y)/(25.0+10.0*cos(y+z)+cos(y)*cos(y*z));
     u[VAR::U_GT2] = -5.0*sin(z)/(25.0+10.0*cos(y+x)+cos(y*x)*cos(z));

     u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));
     //u[F_CHI][pp] = 2.0;

     u[VAR::U_SYMGT0] = 1.00+0.2*sin(x+z)*cos(y);
     u[VAR::U_SYMGT3] = 1.00+0.2*cos(y)*cos(z+ x);
     u[VAR::U_SYMGT5] = 1.00 / ( u[VAR::U_SYMGT0] + u[VAR::U_SYMGT3]);
     u[VAR::U_SYMGT1] = 0.7*cos(x*x + y*y);
     u[VAR::U_SYMGT2] = 0.3*sin(z)*cos(x);
     u[VAR::U_SYMGT4] = -0.5*sin(x*x)*cos(y)*cos(z);

     u[VAR::U_K] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
                   +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
                   +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
                    *exp(-4.0*cos(x)*sin(y))*cos(z);

     u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x)
     -0.3333333333*exp(4.0*cos(x)*sin(y))
     *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y))
     /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))
     /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)
     +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
     u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
     u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
     u[VAR::U_SYMAT3] =
     exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
     u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
     u[VAR::U_SYMAT5] =
     exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
     */

    const double x  = GRIDX_TO_X(xx1);
    const double y  = GRIDY_TO_Y(yy1);
    const double z  = GRIDZ_TO_Z(zz1);

    const double pi = acos(-1.0);
    const double f1 = 31.0 / 17.0;
    const double f2 = 37.0 / 11.0;

    u[VAR::U_ALPHA] = 1.0 - 0.25 * sin(f1 * x);
    // u[F_ALPHA][pp] = 1.0;
    u[VAR::U_BETA0] = 4.0 / 17.0 * sin(x) * cos(z);
    u[VAR::U_BETA1] = pi / 5.0 * cos(y) * sin(z + x);
    u[VAR::U_BETA2] = 4.0 / 17.0 * sin(f2 * x) * sin(z);

    u[VAR::U_B0]    = 31.0 * x * cos(f1 * z + y);
    u[VAR::U_B1]    = 7.0 * y * sin(f1 * x + y) + 3.0 * cos(z);
    u[VAR::U_B2]    = 5.0 * z * cos(f1 * x + y) + 7.0 * sin(z + y + x) + 1.0;

    u[VAR::U_GT0] =
        5.0 * cos(x) / (10.0 * sin(x + z) + 26.0 - 1.0 * cos(x * z) * cos(x));
    u[VAR::U_GT1] =
        -5.0 * sin(y) / (25.0 + 10.0 * cos(y + z) + cos(y) * cos(y * z));
    u[VAR::U_GT2] =
        -5.0 * sin(z) / (25.0 + 10.0 * cos(y + x) + cos(y * x) * cos(z));

    u[VAR::U_CHI]    = 1.0 + exp(-4.0 * cos(x) * sin(y));
    // u[F_CHI][pp] = 2.0;

    u[VAR::U_SYMGT0] = 1.00 + 0.2 * sin(x + z) * cos(y);
    u[VAR::U_SYMGT3] = 1.00 + 0.2 * cos(y) * cos(z + x);
    u[VAR::U_SYMGT5] = 1.00 / (u[VAR::U_SYMGT0] + u[VAR::U_SYMGT3]);
    u[VAR::U_SYMGT1] = 0.07 * (2.0 + cos(x * x + y * y));
    u[VAR::U_SYMGT2] = 0.1 * (3.0 + sin(z) * cos(x));
    u[VAR::U_SYMGT4] = 0.15 * (1.751 - sin(x * x) * cos(y) * cos(z));

    u[VAR::U_K] = 5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
                  5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
                  0.4 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                      exp(-4.0 * cos(x) * sin(y)) * cos(z);
    u[VAR::U_K] *= 0.01234;

    u[VAR::U_SYMAT0] =
        exp(-4.0 * cos(x) * sin(y)) *
        (cos(x) -
         0.3333333333 * exp(4.0 * cos(x) * sin(y)) * (1.0 + 0.2 * sin(x)) *
             (5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
              5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
              0.04 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                  exp(-4.0 * cos(x) * sin(y)) * cos(z)));
    u[VAR::U_SYMAT1] = 1.0 + x * z / (0.1 + x * x + y * y + z * z);
    u[VAR::U_SYMAT2] =
        1.3 - x * y / (3.0 + x * x + 2.0 * y * y + z * z) * (x * x + z * z);
    u[VAR::U_SYMAT3] =
        exp(-4.0 * cos(x) * sin(y)) *
        (cos(y) -
         0.33333333330 * exp(4 * cos(x) * sin(y)) * (1 + 0.2 * cos(y)) *
             (5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
              5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
              0.04 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                  exp(-4.0 * cos(x) * sin(y)) * cos(z)));
    u[VAR::U_SYMAT4] = -1.0 + y * z / (1.0 + 3.0 * x * x + y * y + z * z);
    u[VAR::U_SYMAT5] =
        exp(-4.0 * cos(x) * sin(y)) *
        (cos(z) -
         0.3333333333 * exp(4 * cos(x) * sin(y)) / (1 + 0.2 * sin(x)) /
             (1 + 0.2 * cos(y)) *
             (5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
              5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
              0.04 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                  exp(-4.0 * cos(x) * sin(y)) * cos(z)));

    /* Enforce BSSN constraints */
    double gtd[3][3], Atd[3][3];

    gtd[0][0]              = u[VAR::U_SYMGT0];
    gtd[0][1]              = u[VAR::U_SYMGT1];
    gtd[0][2]              = u[VAR::U_SYMGT2];
    gtd[1][0]              = gtd[0][1];
    gtd[1][1]              = u[VAR::U_SYMGT3];
    gtd[1][2]              = u[VAR::U_SYMGT4];
    gtd[2][0]              = gtd[0][2];
    gtd[2][1]              = gtd[1][2];
    gtd[2][2]              = u[VAR::U_SYMGT5];

    Atd[0][0]              = u[VAR::U_SYMAT0];
    Atd[0][1]              = u[VAR::U_SYMAT1];
    Atd[0][2]              = u[VAR::U_SYMAT2];
    Atd[1][0]              = Atd[0][1];
    Atd[1][1]              = u[VAR::U_SYMAT3];
    Atd[1][2]              = u[VAR::U_SYMAT4];
    Atd[2][0]              = Atd[0][2];
    Atd[2][1]              = Atd[1][2];
    Atd[2][2]              = u[VAR::U_SYMAT5];

    const double one_third = 1.0 / 3.0;
    double det_gtd =
        gtd[0][0] * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]) -
        gtd[0][1] * gtd[0][1] * gtd[2][2] +
        2.0 * gtd[0][1] * gtd[0][2] * gtd[1][2] -
        gtd[0][2] * gtd[0][2] * gtd[1][1];

    if (det_gtd < 0.0) {
        /* FIXME What to do here? The metric is not physical. Do we reset the
         * metric to be flat? */
        gtd[0][0] = 1.0;
        gtd[0][1] = 0.0;
        gtd[0][2] = 0.0;
        gtd[1][0] = 0.0;
        gtd[1][1] = 1.0;
        gtd[1][2] = 0.0;
        gtd[2][0] = 0.0;
        gtd[2][1] = 0.0;
        gtd[2][2] = 1.0;
        det_gtd   = 1.0;
    }
    double det_gtd_to_neg_third = 1.0 / pow(det_gtd, one_third);

    for (unsigned int j = 0; j < 3; j++) {
        for (unsigned int i = 0; i < 3; i++) {
            gtd[i][j] *= det_gtd_to_neg_third;
        }
    }

    det_gtd = gtd[0][0] * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]) -
              gtd[0][1] * gtd[0][1] * gtd[2][2] +
              2.0 * gtd[0][1] * gtd[0][2] * gtd[1][2] -
              gtd[0][2] * gtd[0][2] * gtd[1][1];

    double detgt_m1 = det_gtd - 1.0;

    if (fabs(detgt_m1) > 1.0e-6) {
        std::cout.precision(14);
        std::cout << "enforce_bssn_constraint: det(gtd) != 1. det="
                  << std::fixed << det_gtd << std::endl;
        std::cout << "      gtd(1,1)=" << gtd[0][0] << std::endl;
        std::cout << "      gtd(1,2)=" << gtd[0][1] << std::endl;
        std::cout << "      gtd(1,3)=" << gtd[0][2] << std::endl;
        std::cout << "      gtd(2,2)=" << gtd[1][1] << std::endl;
        std::cout << "      gtd(2,3)=" << gtd[1][2] << std::endl;
        std::cout << "      gtd(3,3)=" << gtd[2][2] << std::endl;
    }

    double gtu[3][3];
    double idet_gtd = 1.0 / det_gtd;
    gtu[0][0] = idet_gtd * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]);
    gtu[0][1] = idet_gtd * (-gtd[0][1] * gtd[2][2] + gtd[0][2] * gtd[1][2]);
    gtu[0][2] = idet_gtd * (gtd[0][1] * gtd[1][2] - gtd[0][2] * gtd[1][1]);
    gtu[1][0] = gtu[0][1];
    gtu[1][1] = idet_gtd * (gtd[0][0] * gtd[2][2] - gtd[0][2] * gtd[0][2]);
    gtu[1][2] = idet_gtd * (-gtd[0][0] * gtd[1][2] + gtd[0][1] * gtd[0][2]);
    gtu[2][0] = gtu[0][2];
    gtu[2][1] = gtu[1][2];
    gtu[2][2] = idet_gtd * (gtd[0][0] * gtd[1][1] - gtd[0][1] * gtd[0][1]);

    /* Require Atd to be traceless. */
    double one_third_trace_Atd =
        one_third *
        (Atd[0][0] * gtu[0][0] + Atd[1][1] * gtu[1][1] + Atd[2][2] * gtu[2][2] +
         2.0 * (Atd[0][1] * gtu[0][1] + Atd[0][2] * gtu[0][2] +
                Atd[1][2] * gtu[1][2]));

    Atd[0][0] -= one_third_trace_Atd * gtd[0][0];
    Atd[0][1] -= one_third_trace_Atd * gtd[0][1];
    Atd[0][2] -= one_third_trace_Atd * gtd[0][2];
    Atd[1][1] -= one_third_trace_Atd * gtd[1][1];
    Atd[1][2] -= one_third_trace_Atd * gtd[1][2];
    Atd[2][2] -= one_third_trace_Atd * gtd[2][2];

    double tr_A = Atd[0][0] * gtu[0][0] + Atd[1][1] * gtu[1][1] +
                  Atd[2][2] * gtu[2][2] +
                  2.0 * (Atd[0][1] * gtu[0][1] + Atd[0][2] * gtu[0][2] +
                         Atd[1][2] * gtu[1][2]);

    if (fabs(tr_A) > 1.0e-6) {
        std::cout << "enforce_bssn_constraint: tr_A != 0. tr_A=" << tr_A
                  << std::endl;
        std::cout << "      Atd(1,1)=" << Atd[0][0] << std::endl;
        std::cout << "      Atd(1,2)=" << Atd[0][1] << std::endl;
        std::cout << "      Atd(1,3)=" << Atd[0][2] << std::endl;
        std::cout << "      Atd(2,2)=" << Atd[1][1] << std::endl;
        std::cout << "      Atd(2,3)=" << Atd[1][2] << std::endl;
        std::cout << "      Atd(3,3)=" << Atd[2][2] << std::endl;
    }

    u[VAR::U_SYMAT0] = Atd[0][0];
    u[VAR::U_SYMAT1] = Atd[0][1];
    u[VAR::U_SYMAT2] = Atd[0][2];
    u[VAR::U_SYMAT3] = Atd[1][1];
    u[VAR::U_SYMAT4] = Atd[1][2];
    u[VAR::U_SYMAT5] = Atd[2][2];

    u[VAR::U_SYMGT0] = gtd[0][0];
    u[VAR::U_SYMGT1] = gtd[0][1];
    u[VAR::U_SYMGT2] = gtd[0][2];
    u[VAR::U_SYMGT3] = gtd[1][1];
    u[VAR::U_SYMGT4] = gtd[1][2];
    u[VAR::U_SYMGT5] = gtd[2][2];
}

void minkowskiInitialData(const double xx1, const double yy1, const double zz1,
                          double* var) {
    // Flat space initialization!
    var[VAR::U_ALPHA]  = 1;  // lapse
    var[VAR::U_CHI]    = 1;  // chi
    var[VAR::U_K]      = 0;  // trace K
    var[VAR::U_GT0]    = 0;  // Gt0
    var[VAR::U_GT1]    = 0;  // Gt1
    var[VAR::U_GT2]    = 0;  // Gt2
    var[VAR::U_BETA0]  = 0;  // shift 0
    var[VAR::U_BETA1]  = 0;  // shift 1
    var[VAR::U_BETA2]  = 0;  // shift 2
    var[VAR::U_B0]     = 0;  // gaugeB0
    var[VAR::U_B1]     = 0;  // gaugeB1
    var[VAR::U_B2]     = 0;  // gaugeB2
    var[VAR::U_SYMGT0] = 1;  // gt11
    var[VAR::U_SYMGT1] = 0;  // gt12
    var[VAR::U_SYMGT2] = 0;  // gt13
    var[VAR::U_SYMGT3] = 1;  // gt22
    var[VAR::U_SYMGT4] = 0;  // gt23
    var[VAR::U_SYMGT5] = 1;  // gt33
    var[VAR::U_SYMAT0] = 0;  // At11
    var[VAR::U_SYMAT1] = 0;  // At12
    var[VAR::U_SYMAT2] = 0;  // At13
    var[VAR::U_SYMAT3] = 0;  // At22
    var[VAR::U_SYMAT4] = 0;  // At23
    var[VAR::U_SYMAT5] = 0;  // At33
}

void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,
                         const Point& pt_min, const Point& pt_max,
                         const unsigned int regLev, const unsigned int maxDepth,
                         MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    double pt_g_min[3];
    double pt_g_max[3];

    pt_g_min[0] = X_TO_GRIDX(pt_min.x());
    pt_g_min[1] = Y_TO_GRIDY(pt_min.y());
    pt_g_min[2] = Z_TO_GRIDZ(pt_min.z());

    pt_g_max[0] = X_TO_GRIDX(pt_max.x());
    pt_g_max[1] = Y_TO_GRIDY(pt_max.y());
    pt_g_max[2] = Z_TO_GRIDZ(pt_max.z());

    assert(pt_g_min[0] >= 0 && pt_g_min[0] <= (1u << maxDepth));
    assert(pt_g_min[1] >= 0 && pt_g_min[1] <= (1u << maxDepth));
    assert(pt_g_min[2] >= 0 && pt_g_min[2] <= (1u << maxDepth));

    assert(pt_g_max[0] >= 0 && pt_g_max[0] <= (1u << maxDepth));
    assert(pt_g_max[1] >= 0 && pt_g_max[1] <= (1u << maxDepth));
    assert(pt_g_max[2] >= 0 && pt_g_max[2] <= (1u << maxDepth));

    unsigned int xRange_b, xRange_e;
    unsigned int yRange_b = pt_g_min[1], yRange_e = pt_g_max[1];
    unsigned int zRange_b = pt_g_min[2], zRange_e = pt_g_max[2];

    xRange_b =
        pt_g_min[0];  //(rank*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];
    xRange_e =
        pt_g_max[1];  //((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

    unsigned int stepSz = 1u << (maxDepth - regLev);

    /* std::cout<<" x min: "<<xRange_b<<" x_max: "<<xRange_e<<std::endl;
     std::cout<<" y min: "<<yRange_b<<" y_max: "<<yRange_e<<std::endl;
     std::cout<<" z min: "<<zRange_b<<" z_max: "<<zRange_e<<std::endl;*/

    for (unsigned int x = xRange_b; x < xRange_e; x += stepSz)
        for (unsigned int y = yRange_b; y < yRange_e; y += stepSz)
            for (unsigned int z = zRange_b; z < zRange_e; z += stepSz) {
                if (x >= (1u << maxDepth)) x = x - 1;
                if (y >= (1u << maxDepth)) y = y - 1;
                if (z >= (1u << maxDepth)) z = z - 1;

                tmpNodes.push_back(
                    ot::TreeNode(x, y, z, regLev, m_uiDim, maxDepth));
            }

    return;
}

double computeWTol(double x, double y, double z, double tolMin) {
    double origin[3];
    origin[0] = (double)(1u << bssn::BSSN_MAXDEPTH - 1);
    origin[1] = (double)(1u << bssn::BSSN_MAXDEPTH - 1);
    origin[2] = (double)(1u << bssn::BSSN_MAXDEPTH - 1);

    double r =
        sqrt(GRIDX_TO_X(x) * GRIDX_TO_X(x) + GRIDY_TO_Y(y) * GRIDY_TO_Y(y) +
             GRIDZ_TO_Z(z) * GRIDZ_TO_Z(z));

    const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
    const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
    const double R1     = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;

    if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 1) {
        return std::min(
            tolMax, std::max(tolMin, (tolMax - tolMin) / (R1 - R0) * (r - R0) +
                                         tolMin));
    } else {
        return tolMin;
    }
}

double computeWTolDCoords(double x, double y, double z, double* hx) {
    // set up a few useful values for computing the wavelet tolerances
    // element order: how many points we have in each grid
    const unsigned int eleOrder = bssn::BSSN_ELE_ORDER;
    // current simulation time
    const double T_CURRENT      = bssn::BSSN_CURRENT_RK_COORD_TIME;
    // radius (from the center of the grid)
    const double r              = sqrt(x * x + y * y + z * z);
    // distance between the BHs
    const double dbh = (bssn::BSSN_BH_LOC[0] - bssn::BSSN_BH_LOC[1]).abs();
    // set up grid point for relative distances to each BH
    Point grid_p(x, y, z);
    // distance from BH0
    const double dbh0 = (grid_p - bssn::BSSN_BH_LOC[0]).abs();
    // distance from BH1
    const double dbh1 = (grid_p - bssn::BSSN_BH_LOC[1]).abs();

    if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 1) {
        const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
        const double tolMin = bssn::BSSN_WAVELET_TOL;

        const double R0     = bssn::BSSN_BH1_AMR_R;
        const double R1     = bssn::BSSN_BH2_AMR_R;

        // R_Max is defined based on the initial separation.
        const double R_MAX =
            (bssn::BH1.getBHCoord() - bssn::BH2.getBHCoord()).abs() + R0 + R1;

#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        if (dbh < 0.1) {
            if ((dbh0 > R_MAX) && (dbh1 > R_MAX)) {
                if (r < (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10))
                    return BSSN_GW_REFINE_WTOL;
                else
                    return tolMax;
            } else {
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            const double xx = x + i * hx[0];
                            const double yy = y + j * hx[1];
                            const double zz = z + k * hx[2];

                            const Point grid_pp(xx, yy, zz);

                            const double dd0 =
                                (grid_pp - bssn::BSSN_BH_LOC[0]).abs();
                            const double dd1 =
                                (grid_pp - bssn::BSSN_BH_LOC[1]).abs();

                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;

                            if (dd0 < R0 || dd1 < R1) {
                                // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                                // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                                // "<<hx[0]<<std::endl;
                                return tolMin;
                            }
                        }

                if (r < (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10))
                    return BSSN_GW_REFINE_WTOL;
                else
                    return tolMax;
            }

        } else {
            if ((dbh0 > R_MAX) &&
                (dbh1 >
                 R_MAX))  // no need to check individual points in the element
                return tolMax;
            else {
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            const double xx = x + i * hx[0];
                            const double yy = y + j * hx[1];
                            const double zz = z + k * hx[2];

                            const Point grid_pp(xx, yy, zz);

                            const double dd0 =
                                (grid_pp - bssn::BSSN_BH_LOC[0]).abs();
                            const double dd1 =
                                (grid_pp - bssn::BSSN_BH_LOC[1]).abs();

                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;

                            if (dd0 < R0 || dd1 < R1) {
                                // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                                // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                                // "<<hx[0]<<std::endl;
                                return tolMin;
                            }
                        }

                //@milinda 21/11/2020 - smooth transition of the wtol.
                if (dbh0 < dbh1)
                    return std::min(
                        tolMax,
                        std::max(tolMin, (tolMax - tolMin) / (R_MAX - R0) *
                                                 (dbh0 - R0) +
                                             tolMin));
                else
                    return std::min(
                        tolMax,
                        std::max(tolMin, (tolMax - tolMin) / (R_MAX - R1) *
                                                 (dbh1 - R1) +
                                             tolMin));
            }
        }

#else
        if ((dbh0 > R_MAX) &&
            (dbh1 >
             R_MAX))  // no need to check individual points in the element
            return tolMax;
        else {
            for (unsigned int k = 0; k < (eleOrder + 1); k++)
                for (unsigned int j = 0; j < (eleOrder + 1); j++)
                    for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                        const double xx = x + i * hx[0];
                        const double yy = y + j * hx[1];
                        const double zz = z + k * hx[2];

                        const Point grid_pp(xx, yy, zz);

                        const double dd0 =
                            (grid_pp - bssn::BSSN_BH_LOC[0]).abs();
                        const double dd1 =
                            (grid_pp - bssn::BSSN_BH_LOC[1]).abs();

                        // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<" dd0:
                        // "<<dd0<<" dd1: "<<dd1<<" hx: "<<hx[0]<<std::endl;

                        if (dd0 < R0 || dd1 < R1) {
                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;
                            return tolMin;
                        }
                    }

            //@milinda 21/11/2020 - smooth transition of the wtol.
            if (dbh0 < dbh1)
                return std::min(
                    tolMax, std::max(tolMin, (tolMax - tolMin) / (R_MAX - R0) *
                                                     (dbh0 - R0) +
                                                 tolMin));
            else
                return std::min(
                    tolMax, std::max(tolMin, (tolMax - tolMin) / (R_MAX - R1) *
                                                     (dbh1 - R1) +
                                                 tolMin));
        }
#endif

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 2) {
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        if (dbh < 0.1) {
            const double R0 =
                (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10);
            const double R1 =
                (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 20);

            const double tolMin =
                std::min(bssn::BSSN_WAVELET_TOL, bssn::BSSN_GW_REFINE_WTOL);
            const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;

            return std::min(
                tolMax,
                std::max(tolMin,
                         (tolMax - tolMin) / (R1 - R0) * (r - R0) + tolMin));

        } else {
            const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double R1     = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
            const double tolMin = bssn::BSSN_WAVELET_TOL;
            return std::min(
                tolMax,
                std::max(tolMin,
                         (tolMax - tolMin) / (R1 - R0) * (r - R0) + tolMin));
        }
#else
        const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
        const double R1     = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
        const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
        const double tolMin = bssn::BSSN_WAVELET_TOL;
        return std::min(
            tolMax, std::max(tolMin, (tolMax - tolMin) / (R1 - R0) * (r - R0) +
                                         tolMin));
#endif

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 3) {
        const double GW_R_SAFETY_FAC = 10.0;
        const double TIME_OFFSET_FAC = 5.0;

        if (T_CURRENT > bssn::BSSN_WAVELET_TOL_FUNCTION_R1 + TIME_OFFSET_FAC) {
            const double RR =
                std::min(T_CURRENT - TIME_OFFSET_FAC,
                         GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10.0);
            const double GW_TOL = bssn::BSSN_GW_REFINE_WTOL;
            const double W_RR   = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double WTOL_EXP = 10.0;
            // const double R1 = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double WTOL_EXP_FAC =
                (RR - R0) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);
            if (r < R0)
                return bssn::BSSN_WAVELET_TOL;
            else
                return std::min(bssn::BSSN_WAVELET_TOL_MAX,
                                ((std::pow(WTOL_EXP, (r - R0) / WTOL_EXP_FAC)) *
                                 bssn::BSSN_WAVELET_TOL));

        } else {
            const double R0       = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double R1       = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double GW_TOL   = bssn::BSSN_GW_REFINE_WTOL;
            const double WTOL_EXP = 10.0;
            const double W_RR = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double WTOL_EXP_FAC =
                (R1 - R0) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);

            if (r < R0)
                return bssn::BSSN_WAVELET_TOL;
            else
                return std::min(bssn::BSSN_WAVELET_TOL_MAX,
                                ((std::pow(WTOL_EXP, (r - R0) / WTOL_EXP_FAC)) *
                                 bssn::BSSN_WAVELET_TOL));
        }

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 4) {
        const double GW_R_SAFETY_FAC = 10.0;
        const double TIME_OFFSET_FAC = 20.0;

        if (T_CURRENT > bssn::BSSN_WAVELET_TOL_FUNCTION_R1 + TIME_OFFSET_FAC) {
            const double RR =
                std::min(T_CURRENT - TIME_OFFSET_FAC,
                         GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10.0);
            const double GW_TOL = bssn::BSSN_GW_REFINE_WTOL;
            const double W_RR   = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double WTOL_EXP = 10.0;
            // const double R1 = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double WTOL_EXP_FAC =
                (RR - R0) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);
            if (r < R0)
                return bssn::BSSN_WAVELET_TOL;
            else
                return std::min(bssn::BSSN_WAVELET_TOL_MAX,
                                ((std::pow(WTOL_EXP, (r - R0) / WTOL_EXP_FAC)) *
                                 bssn::BSSN_WAVELET_TOL));

        } else {
            const double R01      = 3.0 * bssn::BSSN_BH1_MASS;
            const double R02      = 3.0 * bssn::BSSN_BH2_MASS;
            const double R11      = 4.0 * R01;
            const double R12      = 4.0 * R02;
            const double GW_TOL   = bssn::BSSN_GW_REFINE_WTOL;
            const double WTOL_EXP = 10.0;
            const double W_RR = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double W1 =
                (R11 - R01) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);
            const double W2 =
                (R12 - R02) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);

            if (dbh0 < R01 || dbh1 < R02)
                return bssn::BSSN_WAVELET_TOL;
            else {
                double minbheps =
                    std::min(((std::pow(WTOL_EXP, (dbh0 - R01) / W1)) *
                              bssn::BSSN_WAVELET_TOL),
                             ((std::pow(WTOL_EXP, (dbh1 - R02) / W2)) *
                              bssn::BSSN_WAVELET_TOL));
                return std::min(bssn::BSSN_WAVELET_TOL_MAX, minbheps);
            }
        }

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 5) {
        Point grid_p(x, y, z);
        const double d1      = (grid_p - bssn::BSSN_BH_LOC[0]).abs();
        const double d2      = (grid_p - bssn::BSSN_BH_LOC[1]).abs();
        const double m1      = bssn::BSSN_BH1_MASS;
        const double m2      = bssn::BSSN_BH2_MASS;
        const double toffset = 16.0;

        const double eps[3]  = {bssn::BSSN_WAVELET_TOL,
                                bssn::BSSN_GW_REFINE_WTOL,
                                bssn::BSSN_WAVELET_TOL_MAX};
        double rad[3];
        rad[0]    = 3.0 * m1;
        rad[1]    = 4.0 * rad[0];
        rad[3]    = GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1];

        double e1 = CalTolHelper(T_CURRENT, d1, rad, eps, toffset);

        rad[0]    = 3.0 * m2;
        rad[1]    = 4.0 * rad[0];
        rad[3]    = GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1];
        double e2 = CalTolHelper(T_CURRENT, d2, rad, eps, toffset);

        return std::min(e1, e2);

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 6) {
        // WKB Aug 2024
        // use different sensitivities for regions of spacetime which
        // are causally connected to the BHs as cf regions which are
        // spacelike, causally disconnected from the BHs.

        ////////////////////////////////////////////////////////////////
        // set up constants used in this function

        // (max) orbital radius; use strictest refinement here
        const double R_orbit     = 8;
        // outer radius of simulation
        const double R_max       = 400;

        // expected lapse wave tail length (M) + backreflections
        const double L           = 120;
        // calculate the time after which a given radius's relationship
        // with the grid center is both time-like & clean of lapse noise
        const double t_lim       = std::max(r, (r + L) / std::sqrt(2));

        // wavelet tolerance in acausal (or dirty) regions.
        const double eps_disable = .001;
        // time to fade from eps_disable to eps_goal
        const double t_fade      = 100;

        ////////////////////////////////////////////////////////////////
        // set up goal resolution to hit in causal clean regions
        // linearly interpolate log tolerances vs log radii

        double eps_goal;
        if (r <= R_orbit) {
            eps_goal = bssn::BSSN_WAVELET_TOL;
        } else {  // log falloff
            // power we're raising the next expression to, scaling out radius
            const double pwr =
                std::log(r / R_orbit) / std::log(R_max / R_orbit);
            // goal wavelet tolerance at end times
            eps_goal =
                bssn::BSSN_WAVELET_TOL *
                std::pow(bssn::BSSN_WAVELET_TOL_MAX / bssn::BSSN_WAVELET_TOL,
                         pwr);
        }

        ////////////////////////////////////////////////////////////////
        // return time-delayed & smoothed wavelet tolerance

        if (T_CURRENT < t_lim) {  // in spacelike or dirty region
            // return max permissible wavelet tolerance
            // effectively disabling / kneecapping WAMR
            return eps_disable;
        } else if (T_CURRENT > t_lim + t_fade) {  // in clean timelike region
            // return standard wavelet tolerance
            return eps_goal;
        } else {  // in transition region
            // return linear transition between log tolerance values
            // slope of transition region
            const double slope = std::log10(eps_goal / eps_disable) / t_fade;
            const double lg_eps =
                std::log10(eps_disable) + slope * (T_CURRENT - t_lim);
            return std::pow(10.0, lg_eps);
        }
    } else {
        // return global wavelet tolerance, irrespective of position
        return bssn::BSSN_WAVELET_TOL;
    }
}

void writeBLockToBinary(const double** unzipVarsRHS, unsigned int offset,
                        const double* pmin, const double* pmax, double* bxMin,
                        double* bxMax, const unsigned int* sz,
                        unsigned int blkSz, double dxFactor,
                        const char* fprefix) {
    const unsigned int nx       = sz[0];
    const unsigned int ny       = sz[1];
    const unsigned int nz       = sz[2];

    const unsigned int ib       = 3;
    const unsigned int jb       = 3;
    const unsigned int kb       = 3;

    const unsigned int ie       = nx - 3;
    const unsigned int je       = ny - 3;
    const unsigned int ke       = nz - 3;

    const unsigned int blkInlSz = (nx - 3) * (ny - 3) * (nz - 3);

    double hx                   = (pmax[0] - pmin[0]) / (nx - 1);
    double hy                   = (pmax[1] - pmin[1]) / (ny - 1);
    double hz                   = (pmax[2] - pmin[2]) / (nz - 1);

    const double dx = (bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                      (1.0 / (double)(1u << bssn::BSSN_MAXDEPTH));
    unsigned int level = bssn::BSSN_MAXDEPTH - ((unsigned int)(hx / dx) - 1);

    MPI_Comm comm      = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    // std::cout<<"ranl: "<<rank<<"npes: "<<npes<<std::endl;

    // std::cout<<"nx: "<<nx<<" level: "<<level<<" hx: "<<hx<<" dx:
    // "<<dx<<std::endl;

    if ((hx > (dxFactor * dx)) ||
        (pmin[0] < bxMin[0] || pmin[1] < bxMin[1] || pmin[2] < bxMin[2]) ||
        (pmax[0] > bxMax[0] || pmax[1] > bxMax[1] || pmax[2] > bxMax[2]))
        return;

    double* blkInternal = new double[blkInlSz];
    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++) {
        char fName[256];
        sprintf(fName, "%s_%s_n_%d_r_%d_p_%d.bin", fprefix,
                bssn::BSSN_VAR_NAMES[var], nx, rank, npes);
        FILE* outfile = fopen(fName, "w");
        if (outfile == NULL) {
            std::cout << fName << " file open failed " << std::endl;
        }

        for (unsigned int k = kb; k < ke; k++)
            for (unsigned int j = jb; j < je; j++)
                for (unsigned int i = ib; i < ie; i++)
                    blkInternal[k * (ny - 3) * (nx - 3) + j * (nx - 3) + i] =
                        unzipVarsRHS[var]
                                    [offset + k * (ny * nx) + j * (ny) + i];

        fwrite(blkInternal, sizeof(double), blkInlSz,
               outfile);  // write out the number of elements.
        fclose(outfile);
    }

    delete[] blkInternal;
}

unsigned int getOctantWeight(const ot::TreeNode* pNode) {
    return (1u << (3 * pNode->getLevel())) * 1;
}

void computeBHLocations(const ot::Mesh* pMesh, const Point* in, Point* out,
                        double** zipVars, double dt) {
    MPI_Comm commActive = pMesh->getMPICommunicator();

    Point grid_limits[2];
    Point domain_limits[2];

    grid_limits[0]   = Point(bssn::BSSN_OCTREE_MIN[0], bssn::BSSN_OCTREE_MIN[1],
                             bssn::BSSN_OCTREE_MIN[2]);
    grid_limits[1]   = Point(bssn::BSSN_OCTREE_MAX[0], bssn::BSSN_OCTREE_MAX[1],
                             bssn::BSSN_OCTREE_MAX[2]);

    domain_limits[0] = Point(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                             bssn::BSSN_COMPD_MIN[2]);
    domain_limits[1] = Point(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                             bssn::BSSN_COMPD_MAX[2]);

    double beta0[2];
    double beta1[2];
    double beta2[2];

    std::vector<unsigned int> validIndex_beta0;
    std::vector<unsigned int> validIndex_beta1;
    std::vector<unsigned int> validIndex_beta2;

    double beta3vec[6] = {0, 0, 0, 0, 0, 0};
    double bh_pts[6]   = {0, 0, 0, 0, 0, 0};

    bh_pts[0]          = in[0].x();
    bh_pts[1]          = in[0].y();
    bh_pts[2]          = in[0].z();

    bh_pts[3]          = in[1].x();
    bh_pts[4]          = in[1].y();
    bh_pts[5]          = in[1].z();

    if (pMesh->isActive()) {
        unsigned int activeRank = pMesh->getMPIRank();

        ot::da::interpolateToCoords(pMesh, zipVars[VAR::U_BETA0], bh_pts, 6,
                                    grid_limits, domain_limits, beta0,
                                    validIndex_beta0);
        ot::da::interpolateToCoords(pMesh, zipVars[VAR::U_BETA1], bh_pts, 6,
                                    grid_limits, domain_limits, beta1,
                                    validIndex_beta1);
        ot::da::interpolateToCoords(pMesh, zipVars[VAR::U_BETA2], bh_pts, 6,
                                    grid_limits, domain_limits, beta2,
                                    validIndex_beta2);

        assert(validIndex_beta0.size() == validIndex_beta1.size());
        assert(validIndex_beta1.size() == validIndex_beta2.size());
    }

    unsigned int red_ranks[2]   = {0, 0};
    unsigned int red_ranks_g[2] = {0, 0};
    // global bcast
    for (unsigned int ind = 0; ind < validIndex_beta0.size(); ind++) {
        assert(validIndex_beta0[ind] == validIndex_beta1[ind]);
        assert(validIndex_beta0[ind] == validIndex_beta2[ind]);
        const unsigned int gRank                = pMesh->getMPIRankGlobal();

        beta3vec[validIndex_beta0[ind] * 3 + 0] = beta0[validIndex_beta0[ind]];
        beta3vec[validIndex_beta1[ind] * 3 + 1] = beta1[validIndex_beta1[ind]];
        beta3vec[validIndex_beta2[ind] * 3 + 2] = beta2[validIndex_beta2[ind]];

        // std::cout<<"rank: "<<gRank<<"beta["<<(validIndex_beta0[ind]*3)<<"]: (
        // "<<beta3vec[validIndex_beta0[ind]*3 + 0]<<",
        // "<<beta3vec[validIndex_beta0[ind]*3 + 1]<<",
        // "<<beta3vec[validIndex_beta0[ind]*3 + 2]<<")"<<std::endl;
        red_ranks[validIndex_beta2[ind]]        = gRank;
    }

    par::Mpi_Allreduce(red_ranks, red_ranks_g, 2, MPI_MAX,
                       pMesh->getMPIGlobalCommunicator());
    MPI_Bcast(&beta3vec[0], 3, MPI_DOUBLE, red_ranks_g[0],
              pMesh->getMPIGlobalCommunicator());
    MPI_Bcast(&beta3vec[3], 3, MPI_DOUBLE, red_ranks_g[1],
              pMesh->getMPIGlobalCommunicator());

    // if(!pMesh->getMPIRankGlobal())
    //     std::cout<<"beta bh0: ( "<<beta3vec[0]<<", "<<beta3vec[1]<<",
    //     "<<beta3vec[2]<<") :  beta 1 ( "<<beta3vec[3]<<", "<<beta3vec[4]<<",
    //     "<<beta3vec[5]<<") "<<std::endl;

    double x[2], y[2], z[2];
    for (unsigned int bh = 0; bh < 2; bh++) {
        x[bh]   = in[bh].x() - beta3vec[bh * 3 + 0] * dt;
        y[bh]   = in[bh].y() - beta3vec[bh * 3 + 1] * dt;
        z[bh]   = in[bh].z() - beta3vec[bh * 3 + 2] * dt;

        out[bh] = Point(x[bh], y[bh], z[bh]);
    }

    return;
}

ot::Mesh* weakScalingReMesh(ot::Mesh* pMesh, unsigned int target_npes) {
    int rank, npes;
    MPI_Comm comm = pMesh->getMPIGlobalCommunicator();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (target_npes > npes) {
        if (!rank)
            RAISE_ERROR("target npes "
                        << target_npes
                        << " is larger than global npes:" << npes);

        MPI_Abort(comm, 0);
    }

    const double R_RES_FAC         = 10;
    const double dr1               = bssn::BSSN_BH1_AMR_R / R_RES_FAC;
    const double dr2               = bssn::BSSN_BH2_AMR_R / R_RES_FAC;

    const DendroIntL ELE_SZ_REQ_GL = bssn::BSSN_DENDRO_GRAIN_SZ * target_npes;
    const unsigned int MAX_ITER    = 30;

    unsigned int iter_count        = 0;
    bool is_converged              = true;

    // Sequential split of each level, such that elemental grain size is
    // satiesfied.
    ot::Mesh* current_mesh         = NULL;

    do {
        if (current_mesh == NULL) current_mesh = pMesh;

        unsigned int LMIN, LMAX;
        current_mesh->computeMinMaxLevel(LMIN, LMAX);

        DendroIntL localSz  = current_mesh->getNumLocalMeshElements();
        DendroIntL globalSz = 0;

        par::Mpi_Allreduce(&localSz, &globalSz, 1, MPI_SUM, comm);

        DendroIntL req_l_splits =
            std::max(1ll, (ELE_SZ_REQ_GL - globalSz) / (14 * npes));

        std::vector<unsigned int> ref_flags;
        if (current_mesh->isActive()) {
            ref_flags.resize(current_mesh->getNumLocalMeshElements(),
                             OCT_NO_CHANGE);
            const unsigned int active_rank = current_mesh->getMPIRank();
            const unsigned int active_npes = current_mesh->getMPICommSize();
            const MPI_Comm active_comm     = current_mesh->getMPICommunicator();
            // const DendroIntL lsplit_b =
            // (active_rank)*req_l_splits/active_npes; const DendroIntL lsplit_e
            // = (active_rank+1)*req_l_splits/active_npes; const DendroIntL
            // num_split_rank = lsplit_e-lsplit_b;

            const ot::TreeNode* pNodes = current_mesh->getAllElements().data();

            // DendroIntL lcount =0;
            // for (unsigned int ele = current_mesh->getElementLocalBegin(); ele
            // < current_mesh->getElementLocalEnd(); ele++)
            // {
            //     if(pNodes[ele].getLevel() == LMAX-1)
            //         lcount++;
            // }

            // std::vector<DendroIntL> num_lev_counts;
            // std::vector<DendroIntL> num_lev_offsets;

            // num_lev_counts.resize(active_npes,0);
            // num_lev_offsets.resize(active_npes,0);

            // par::Mpi_Allgather(&lcount,num_lev_counts.data(),1,active_comm);
            // omp_par::scan(num_lev_counts.data(),num_lev_offsets.data(),active_npes);

            // auto it_upper =
            // std::upper_bound(num_lev_offsets.begin(),num_lev_offsets.end(),req_l_splits);
            // unsigned int valid_npes=active_npes;

            // if(it_upper != num_lev_offsets.end())
            //     valid_npes = std::distance(num_lev_offsets.begin(),it_upper)
            //     +1;

            // if (active_rank < valid_npes)
            // {
            //     unsigned int ele_offset =
            //     current_mesh->getElementLocalBegin(); lcount =0; for
            //     (unsigned int ele = current_mesh->getElementLocalBegin(); ele
            //     < current_mesh->getElementLocalEnd(); ele++)
            //     {
            //         if( lcount < req_l_splits  && pNodes[ele].getLevel() ==
            //         LMAX-1)
            //         {
            //             ref_flags[ele-ele_offset] = OCT_SPLIT;
            //             lcount++;
            //         }

            //     }

            // }
            unsigned int ele_offset    = current_mesh->getElementLocalBegin();
            DendroIntL lcount          = 0;
            for (unsigned int ele = current_mesh->getElementLocalBegin();
                 ele < current_mesh->getElementLocalEnd(); ele++) {
                if (lcount < req_l_splits &&
                    pNodes[ele].getLevel() == LMAX - 1) {
                    ref_flags[ele - ele_offset] = OCT_SPLIT;
                    lcount++;
                }
            }
        }

        bool is_refine   = current_mesh->setMeshRefinementFlags(ref_flags);
        bool is_refine_g = false;
        MPI_Allreduce(&is_refine, &is_refine_g, 1, MPI_CXX_BOOL, MPI_LOR, comm);
        if (is_refine_g) {
            ot::Mesh* new_mesh =
                current_mesh->ReMesh(bssn::BSSN_DENDRO_GRAIN_SZ);
            if (current_mesh == pMesh)
                current_mesh = new_mesh;
            else {
                std::swap(current_mesh, new_mesh);
                delete new_mesh;
            }

            localSz = current_mesh->getNumLocalMeshElements();
            par::Mpi_Allreduce(&localSz, &globalSz, 1, MPI_SUM, comm);

            if (!rank)
                std::cout << "weak scaling remesh iter : " << iter_count
                          << " num global elements : " << globalSz << std::endl;

            is_converged =
                (globalSz >=
                 ELE_SZ_REQ_GL);  //|| (fabs(globalSz/(double)npes -
                                  // bssn::BSSN_DENDRO_GRAIN_SZ)/(double)bssn::BSSN_DENDRO_GRAIN_SZ
                                  //< 0.1);//(globalSz >= ELE_SZ_REQ_GL);
            iter_count++;

        } else {
            // no remesh triggered hence, is_converged true;
            is_converged = true;
        }

    } while (!is_converged && iter_count < MAX_ITER);

    return current_mesh;
}

void allocate_bssn_deriv_workspace(const ot::Mesh* pMesh, unsigned int s_fac) {
    deallocate_bssn_deriv_workspace();

    if (!pMesh->isActive()) return;

    // gets the largest block size.
    const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
    unsigned int max_blk_sz               = 0;
    for (unsigned int i = 0; i < blkList.size(); i++) {
        unsigned int blk_sz = blkList[i].getAllocationSzX() *
                              blkList[i].getAllocationSzY() *
                              blkList[i].getAllocationSzZ();
        if (blk_sz > max_blk_sz) max_blk_sz = blk_sz;
    }

    if (bssn::BSSN_DERIV_WORKSPACE != nullptr) {
        delete[] bssn::BSSN_DERIV_WORKSPACE;
        bssn::BSSN_DERIV_WORKSPACE = nullptr;
    }

    bssn::BSSN_DERIV_WORKSPACE =
        new double[s_fac * max_blk_sz * bssn::BSSN_NUM_DERIVS];
}

void deallocate_bssn_deriv_workspace() {
    if (bssn::BSSN_DERIV_WORKSPACE != nullptr) {
        delete[] bssn::BSSN_DERIV_WORKSPACE;
        bssn::BSSN_DERIV_WORKSPACE = nullptr;
    }
}

std::tuple<std::string, std::string, std::string> encode_bh_locs(
    const std::vector<std::pair<Point, Point>>& bh_history,
    const std::vector<double>& bh_times) {
    std::vector<unsigned char> bh1_bytes;
    std::vector<unsigned char> bh2_bytes;

    for (const auto& pair : bh_history) {
        double coords1[3] = {pair.first.x(), pair.first.y(), pair.first.z()};
        double coords2[3] = {pair.second.x(), pair.second.y(), pair.second.z()};

        bh1_bytes.insert(
            bh1_bytes.end(), reinterpret_cast<unsigned char*>(coords1),
            reinterpret_cast<unsigned char*>(coords1) + sizeof(coords1));

        bh2_bytes.insert(
            bh2_bytes.end(), reinterpret_cast<unsigned char*>(coords2),
            reinterpret_cast<unsigned char*>(coords2) + sizeof(coords2));
    }

    // with bytes available, now we can encode with base91_encoding
    std::string bh1_str  = base<91>::encode(std::string(
        reinterpret_cast<const char*>(bh1_bytes.data()), bh1_bytes.size()));

    std::string bh2_str  = base<91>::encode(std::string(
        reinterpret_cast<const char*>(bh2_bytes.data()), bh2_bytes.size()));

    std::string time_str = base<91>::encode(
        std::string(reinterpret_cast<const char*>(bh_times.data()),
                    bh_times.size() * sizeof(double)));

    return std::make_tuple(bh1_str, bh2_str, time_str);
}

std::tuple<std::vector<std::pair<Point, Point>>, std::vector<double>>
decode_bh_locs(const std::string& bh1_str, const std::string& bh2_str,
               const std::string& time_str) {
    const size_t double_size = sizeof(double);
    // decode the strings
    std::string bh1_bytes    = base<91>::decode(bh1_str);
    std::string bh2_bytes    = base<91>::decode(bh2_str);
    std::string time_bytes   = base<91>::decode(time_str);

    const size_t num_entries = time_bytes.size() / double_size;

    // with the bytes back in place, we need to do a reinterpret cast for time
    std::vector<double> time_vector;
    std::vector<std::pair<Point, Point>> bh_locs;
    for (size_t i = 0; i < num_entries; ++i) {
        double value;
        std::memcpy(&value, time_bytes.data() + i * double_size, double_size);
        time_vector.push_back(value);

        double bh_temp[3];

        // create the point for b1
        std::memcpy(&bh_temp, bh1_bytes.data() + i * double_size * 3,
                    double_size * 3);
        Point bh1Pt = Point(bh_temp[0], bh_temp[1], bh_temp[2]);

        // then do it for bh2
        std::memcpy(&bh_temp, bh2_bytes.data() + i * double_size * 3,
                    double_size * 3);
        Point bh2Pt = Point(bh_temp[0], bh_temp[1], bh_temp[2]);

        bh_locs.push_back(std::make_pair(bh1Pt, bh2Pt));
    }

    return std::make_tuple(bh_locs, time_vector);
}

}  // end of namespace bssn

namespace bssn {

namespace timer {
void initFlops() {
    total_runtime.start();
    t_f2o.start();
    t_cons.start();
    t_bal.start();
    t_mesh.start();
    t_rkSolve.start();
    t_ghostEx_sync.start();
    t_unzip_sync.start();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_sync_edge.start();
    dendro::timer::t_unzip_sync_vtex.start();
    dendro::timer::t_unzip_p2c.start();
    dendro::timer::t_unzip_sync_nodalval.start();
    dendro::timer::t_unzip_sync_cpy.start();
    dendro::timer::t_unzip_sync_f_c1.start();
    dendro::timer::t_unzip_sync_f_c2.start();
    dendro::timer::t_unzip_sync_f_c3.start();

    t_unzip_async.start();
    dendro::timer::t_unzip_async_comm.start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_async_external.start();
    dendro::timer::t_unzip_async_comm.start();
    t_deriv.start();
    t_rhs.start();

    t_rhs_a.start();
    t_rhs_b.start();
    t_rhs_gt.start();
    t_rhs_chi.start();
    t_rhs_At.start();
    t_rhs_K.start();
    t_rhs_Gt.start();
    t_rhs_B.start();

    t_bdyc.start();

    t_zip.start();
    t_rkStep.start();
    t_isReMesh.start();
    t_gridTransfer.start();
    t_ioVtu.start();
    t_ioCheckPoint.start();
}

void resetSnapshot() {
    total_runtime.snapreset();
    t_f2o.snapreset();
    t_cons.snapreset();
    t_bal.snapreset();
    t_mesh.snapreset();
    t_rkSolve.snapreset();
    t_ghostEx_sync.snapreset();
    t_unzip_sync.snapreset();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].snapreset();

    dendro::timer::t_unzip_sync_internal.snapreset();
    dendro::timer::t_unzip_sync_edge.snapreset();
    dendro::timer::t_unzip_sync_vtex.snapreset();
    dendro::timer::t_unzip_p2c.snapreset();
    dendro::timer::t_unzip_sync_nodalval.snapreset();
    dendro::timer::t_unzip_sync_cpy.snapreset();

    dendro::timer::t_unzip_sync_f_c1.snapreset();
    dendro::timer::t_unzip_sync_f_c2.snapreset();
    dendro::timer::t_unzip_sync_f_c3.snapreset();

    t_unzip_async.snapreset();
    dendro::timer::t_unzip_async_internal.snapreset();
    dendro::timer::t_unzip_async_external.snapreset();
    dendro::timer::t_unzip_async_comm.snapreset();

    t_deriv.snapreset();
    t_rhs.snapreset();

    t_rhs_a.snapreset();
    t_rhs_b.snapreset();
    t_rhs_gt.snapreset();
    t_rhs_chi.snapreset();
    t_rhs_At.snapreset();
    t_rhs_K.snapreset();
    t_rhs_Gt.snapreset();
    t_rhs_B.snapreset();

    t_bdyc.snapreset();

    t_zip.snapreset();
    t_rkStep.snapreset();
    t_isReMesh.snapreset();
    t_gridTransfer.snapreset();
    t_ioVtu.snapreset();
    t_ioCheckPoint.snapreset();
}

void profileInfo(const char* filePrefix, const ot::Mesh* pMesh) {
    int activeRank, activeNpes, globalRank, globalNpes;

    MPI_Comm commActive;
    MPI_Comm commGlobal;

    if (pMesh->isActive()) {
        commActive = pMesh->getMPICommunicator();
        activeRank = pMesh->getMPIRank();
        activeNpes = pMesh->getMPICommSize();
    }

    globalRank = pMesh->getMPIRankGlobal();
    globalNpes = pMesh->getMPICommSizeGlobal();
    commGlobal = pMesh->getMPIGlobalCommunicator();

    double t_stat;
    double t_stat_g[3];

    const char separator = ' ';
    const int nameWidth  = 30;
    const int numWidth   = 10;

    char fName[256];
    std::ofstream outfile;

    DendroIntL localSz, globalSz;

    if (!activeRank) {
        sprintf(fName, "%s_final.prof", filePrefix);
        outfile.open(fName);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        outfile << "active npes : " << activeNpes << std::endl;
        outfile << "global npes : " << globalNpes << std::endl;
        outfile << "partition tol : " << bssn::BSSN_LOAD_IMB_TOL << std::endl;
        outfile << "wavelet tol : " << bssn::BSSN_WAVELET_TOL << std::endl;
        outfile << "maxdepth : " << bssn::BSSN_MAXDEPTH << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    localSz           = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "Elements : " << globalSz << std::endl;

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(zip) : " << globalSz << std::endl;

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(unzip) : " << globalSz << std::endl;

    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(s)" << std::endl;

    t_stat = total_runtime.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "+runtime(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_f2o.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++f2o";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_cons.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++construction";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkSolve.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++rkSolve";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bal.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --2:1 balance";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_mesh.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --mesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkStep.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --rkstep";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ghostEx_sync.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --ghostExchge.";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_sync.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_async.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_async";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
    t_stat = dendro::timer::t_unzip_async_internal.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_internal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_async_external.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_external";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_async_comm.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_comm (comm) ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;
#endif

    t_stat = t_deriv.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --deriv ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bdyc.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --boundary con ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_zip.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --zip";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioVtu.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --vtu";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioCheckPoint.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --checkpoint";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank) outfile.close();
}

void profileInfoIntermediate(const char* filePrefix, const ot::Mesh* pMesh,
                             const unsigned int currentStep) {
    int activeRank, activeNpes, globalRank, globalNpes;

    MPI_Comm commActive;
    MPI_Comm commGlobal;

    if (pMesh->isActive()) {
        commActive = pMesh->getMPICommunicator();
        activeRank = pMesh->getMPIRank();
        activeNpes = pMesh->getMPICommSize();
    }

    globalRank = pMesh->getMPIRankGlobal();
    globalNpes = pMesh->getMPICommSizeGlobal();
    commGlobal = pMesh->getMPIGlobalCommunicator();

    double t_stat;
    double t_stat_g[3];

    const char separator = ' ';
    const int nameWidth  = 30;
    const int numWidth   = 10;

    char fName[256];
    std::ofstream outfile;

    DendroIntL localSz, globalSz;

    DendroIntL ghostElements;
    DendroIntL localElements;

    DendroIntL ghostNodes;
    DendroIntL localNodes;

    DendroIntL totalSendNode;
    DendroIntL totalRecvNode;

    DendroIntL numCalls;

#ifdef BSSN_PROFILE_HUMAN_READABLE
    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        outfile << "active npes : " << activeNpes << std::endl;
        outfile << "global npes : " << globalNpes << std::endl;
        outfile << "current step : " << currentStep << std::endl;
        outfile << "partition tol : " << bssn::BSSN_LOAD_IMB_TOL << std::endl;
        outfile << "wavelet tol : " << bssn::BSSN_WAVELET_TOL << std::endl;
        outfile << "maxdepth : " << bssn::BSSN_MAXDEPTH << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    localSz           = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "Elements : " << globalSz << std::endl;

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(zip) : " << globalSz << std::endl;

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(unzip) : " << globalSz << std::endl;

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    ghostNodes    = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes    = pMesh->getNumLocalMeshNodes();

    if (!rank)
        outfile << "========================= MESH "
                   "==========================================================="
                   "============ "
                << std::endl;

    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(#)" << std::endl;

    t_stat = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = ghostNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "send Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "recv Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank)
        outfile << "========================= RUNTIME "
                   "==========================================================="
                   "======== "
                << std::endl;
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(s)" << std::endl;

    /* t_stat=total_runtime.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"+runtime(s)"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_f2o.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++f2o"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_cons.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++construction"; if(!rank)outfile << std::left
    << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_rkSolve.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++rkSolve"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;*/

    t_stat = t_bal.snap;
    // numCalls=t_bal.num_calls;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++2:1 balance";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++mesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++rkstep";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++ghostExchge.";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_sync";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_async.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_async";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_comm_wait (comm) ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_nodalval.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_nodalVal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c1.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c1";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c2.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c2";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c3.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c3";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_cpy.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_cpy";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_internal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[0].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_left";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[1].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_right";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[2].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_down";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[3].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_up";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[4].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_back";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[5].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_front";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_edge.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_edge";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_vtex.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_vtex";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_p2c.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_p2c";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;
#endif

    /*
    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
    t_stat=t_unzip_async_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_internal"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

    t_stat=t_unzip_async_external.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_external"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_comm (comm) "; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    #endif
    */
    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++isReMesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++gridTransfer";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++deriv ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++compute_rhs ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_a.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_a ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_b.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_b ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_gt.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_gt ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_chi.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_chi ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_At.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_At ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_K.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_K ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_Gt.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_Gt ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_B.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_B ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bdyc.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++boundary con ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_zip.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++zip";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioVtu.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++vtu";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioCheckPoint.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++checkpoint";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank) outfile.close();
#else

    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        // writes the header
        if (currentStep == 0)
            outfile << "step\t act_npes\t glb_npes\t part_tol\t wave_tol\t "
                       "maxdepth\t numOcts\t dof_zip\t dof_unzip\t"
                    << "element_ghost_min\t element_ghost_mean\t "
                       "element_ghost_max\t"
                    << "element_local_min\t element_local_mean\t "
                       "element_local_max\t"
                    << "nodes_local_min\t nodes_local_mean\t nodes_local|max\t"
                    << "send_nodes_min\t send_nodes_mean\t send_nodes_max\t"
                    << "recv_nodes_min\t recv_nodes_mean\t recv_nodes_max\t"
                    << "bal_min\t bal_mean\t bal_max\t"
                    << "mesh_min\t mesh_mean\t mesh_max\t"
                    << "rkstep_min\t rkstep_mean\t rkstep_max\t"
                    << "ghostEx_min\t ghostEx_mean\t ghostEx_max\t"
                    << "unzip_sync_min\t unzip_sync_mean\t unzip_sync_max\t"
                    << "unzip_async_min\t unzip_async_mean\t unzip_async_max\t"
                    << "unzip_async_wait_min\t unzip_async_wait_mean\t "
                       "unzip_async_wait_max\t"
                    << "isRemesh_min\t isRemesh_mean\t isRemesh_max\t"
                    << "GT_min\t GT_mean\t GT_max\t"
                    << "deriv_min\t deriv_mean\t deriv_max\t"
                    << "rhs_min\t rhs_mean\t rhs_max\t" << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    if (!rank) outfile << currentStep << "\t ";
    if (!rank) outfile << activeNpes << "\t ";
    if (!rank) outfile << globalNpes << "\t ";
    if (!rank) outfile << bssn::BSSN_LOAD_IMB_TOL << "\t ";
    if (!rank) outfile << bssn::BSSN_WAVELET_TOL << "\t ";
    if (!rank) outfile << bssn::BSSN_MAXDEPTH << "\t ";

    localSz = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    t_stat        = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    ghostNodes = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes = pMesh->getNumLocalMeshNodes();

    /*t_stat=ghostNodes;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t
    ";*/

    t_stat     = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_bal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_async.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    if (!rank) outfile << std::endl;
    if (!rank) outfile.close();
#endif
}

}  // namespace timer

}  // namespace bssn

namespace GW {
void psi4ShpereDump(const ot::Mesh* mesh, DendroScalar** cVar,
                    unsigned int timestep, double time, double dtheta,
                    double dphi) {
    unsigned int rankGlobal   = mesh->getMPIRankGlobal();
    unsigned int npesGlobal   = mesh->getMPICommSizeGlobal();
    MPI_Comm commGlobal       = mesh->getMPIGlobalCommunicator();

    const unsigned int nTheta = (M_PI) / dtheta;
    const unsigned int nPhi   = (2 * M_PI) / dphi;
    const unsigned int numPts = nTheta * nPhi;

    unsigned int totalModes   = 0;
    for (unsigned int l = 0; l < BSSN_GW_NUM_LMODES; l++)
        totalModes += 2 * BSSN_GW_L_MODES[l] + 1;

    const unsigned int TOTAL_MODES = totalModes;

    DendroComplex* swsh_coeff =
        new DendroComplex[BSSN_GW_NUM_RADAII * TOTAL_MODES];
    DendroComplex* swsh_coeff_g =
        new DendroComplex[BSSN_GW_NUM_RADAII * TOTAL_MODES];

    std::vector<unsigned int> lmCounts;
    std::vector<unsigned int> lmOffset;

    lmCounts.resize(BSSN_GW_NUM_LMODES);
    lmOffset.resize(BSSN_GW_NUM_LMODES);

    for (unsigned int l = 0; l < BSSN_GW_NUM_LMODES; l++)
        lmCounts[l] = 2 * BSSN_GW_L_MODES[l] + 1;

    lmOffset[0] = 0;
    omp_par::scan(&(*(lmCounts.begin())), &(*(lmOffset.begin())),
                  BSSN_GW_NUM_LMODES);

    if (mesh->isActive()) {
        const unsigned int rankActive = mesh->getMPIRank();
        const unsigned int npesActive = mesh->getMPICommSize();

        std::vector<double> coords;
        coords.reserve(3 * numPts);

        std::vector<double> psi4_real;
        psi4_real.resize(numPts);

        std::vector<double> psi4_imag;
        psi4_imag.resize(numPts);

        Point grid_limits[2];
        Point domain_limits[2];

        grid_limits[0] =
            Point(bssn::BSSN_OCTREE_MIN[0], bssn::BSSN_OCTREE_MIN[1],
                  bssn::BSSN_OCTREE_MIN[2]);
        grid_limits[1] =
            Point(bssn::BSSN_OCTREE_MAX[0], bssn::BSSN_OCTREE_MAX[1],
                  bssn::BSSN_OCTREE_MAX[2]);

        domain_limits[0] =
            Point(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                  bssn::BSSN_COMPD_MIN[2]);
        domain_limits[1] =
            Point(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                  bssn::BSSN_COMPD_MAX[2]);

        std::vector<unsigned int> validIndex;

        for (unsigned int k = 0; k < BSSN_GW_NUM_RADAII; k++) {
            for (unsigned int i = 0; i < nTheta; i++)
                for (unsigned int j = 0; j < nPhi; j++) {
                    double x =
                        BSSN_GW_RADAII[k] * sin(j * dtheta) * cos(i * dphi);
                    double y =
                        BSSN_GW_RADAII[k] * sin(j * dtheta) * sin(i * dphi);
                    double z = BSSN_GW_RADAII[k] * cos(j * dtheta);

                    coords.push_back(x);
                    coords.push_back(y);
                    coords.push_back(z);
                }

            validIndex.clear();
            ot::da::interpolateToCoords(
                mesh, cVar[bssn::VAR_CONSTRAINT::C_PSI4_REAL],
                &(*(coords.begin())), coords.size(), grid_limits, domain_limits,
                &(*(psi4_real.begin())), validIndex);

            validIndex.clear();
            ot::da::interpolateToCoords(
                mesh, cVar[bssn::VAR_CONSTRAINT::C_PSI4_IMG],
                &(*(coords.begin())), coords.size(), grid_limits, domain_limits,
                &(*(psi4_imag.begin())), validIndex);
        }
    }
}

}  // namespace GW
