//
// Created by milinda on 7/25/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief This file contains all the parameters related to BSSN simulation.
 */
//

#pragma once
#include <string.h>

#include <iostream>
#include <toml.hpp>

#include "bh.h"
#include "dendro.h"
#include "grDef.h"
#include "memory_pool.h"

namespace bssn {

extern mem::memory_pool<double> BSSN_MEM_POOL;
/**@brief element order*/
extern unsigned int BSSN_ELE_ORDER;

/**@brief padding width for the blocks*/
extern unsigned int BSSN_PADDING_WIDTH;

/**@brief number of variables*/
static const unsigned int BSSN_NUM_VARS            = 24;

/**@brief number of constraints variables*/
static const unsigned int BSSN_CONSTRAINT_NUM_VARS = 6;

/***@brief number of RK45 stages*/
static const unsigned int BSSN_RK45_STAGES         = 6;

/***@brief number of RK4 stages*/
static const unsigned int BSSN_RK4_STAGES          = 4;

/**@brief number of rk4 stages*/
static const unsigned int BSSN_RK3_STAGES          = 3;

/**@brief: parameter used for adaptive time step update. */
static const double BSSN_SAFETY_FAC                = 0.8;

/**@brief number of internal variables*/
static const unsigned int BSSN_NUM_VARS_INTENL =
    (BSSN_RK45_STAGES + 1) * BSSN_NUM_VARS;

/**@brief CFL stability number number (specifies how dt=BSSN_CFL_FACTOR*dx)*/
extern double BSSN_CFL_FACTOR;

/** @brief the minimum dx spacing of the grid, this needs to be stored and
 * updated */
extern double BSSN_CURRENT_MIN_DX;

/** @brief The ability to scale KO Diss by the conformal factor
 *
 * This comes from the followign paper: https://arxiv.org/pdf/2404.01137.pdf
 */
extern bool BSSN_KO_SIGMA_SCALE_BY_CONFORMAL;

/** @brief The ability to scale KO Diss by the conformal factor
 *
 * This comes from the followign paper: https://arxiv.org/pdf/2404.01137.pdf
 */
extern bool BSSN_KO_SIGMA_SCALE_BY_CONFORMAL_POST_MERGER_ONLY;

/** @brief The amount to scale each of the gauge variables by when conformal KO
 * scaling is on */
extern double BSSN_EPSILON_CAKO_GAUGE;

/** @brief The amount to scale each of the other variables by when conformal KO
 * scaling is on */
extern double BSSN_EPSILON_CAKO_OTHER;

/** @brief The boolean that controls whether or not the CAKO should be enabled
 *
 * This is only an internal variable. The BSSN_KO_SIGMA_SCALE_BY_CONFORMAL
 * controls if it should turn on
 */
extern bool BSSN_CAKO_ENABLED;

/** @brief the parameter that controls the CAHD strength */
extern double BSSN_CAHD_C;

/**@brief min bh domain add these to the parameter file.*/
extern double BSSN_COMPD_MIN[3];
/**@brief min bh domain @todo add these to the parameter file. */
extern double BSSN_COMPD_MAX[3];

/**@brief min coords of the OCTREE */
extern double BSSN_OCTREE_MIN[3];
/**@brief max coords of the OCTREE */
extern double BSSN_OCTREE_MAX[3];

/**@brief solution output frequency*/
extern unsigned int BSSN_IO_OUTPUT_FREQ;

/**@brief Gravitational wave data extraction frequency*/
extern unsigned int BSSN_GW_EXTRACT_FREQ;

/**@brief Gravitational wave data extraction frequency (after merger)*/
extern unsigned int BSSN_GW_EXTRACT_FREQ_AFTER_MERGER;

/**@brief timestep norms out put freq.*/
extern unsigned int BSSN_TIME_STEP_OUTPUT_FREQ;

/**@brief remesh test frequency*/
extern unsigned int BSSN_REMESH_TEST_FREQ;

/**@brief remesh test freq. after the merger*/
extern unsigned int BSSN_REMESH_TEST_FREQ_AFTER_MERGER;

/**@brief checkpoint store frequency*/
extern unsigned int BSSN_CHECKPT_FREQ;

/**@brief restore the solver from check point if set to 1. */
extern unsigned int BSSN_RESTORE_SOLVER;

/**@brief use the block adaptivity and disable the AMR*/
extern unsigned int BSSN_ENABLE_BLOCK_ADAPTIVITY;

/**@brief file prefix for VTU*/
extern std::string BSSN_VTU_FILE_PREFIX;

/**@brief file prefix for write check point*/
extern std::string BSSN_CHKPT_FILE_PREFIX;

/**@brief file prefix to write profile info.*/
extern std::string BSSN_PROFILE_FILE_PREFIX;

/**@brief number of refine variables*/
extern unsigned int BSSN_NUM_REFINE_VARS;

/**@brief indices of refine var ids*/
extern unsigned int BSSN_REFINE_VARIABLE_INDICES[BSSN_NUM_VARS];

/**@brief number of evolution variables written to vtu files*/
extern unsigned int BSSN_NUM_EVOL_VARS_VTU_OUTPUT;

/**@brief number of constrint variables written to vtu files*/
extern unsigned int BSSN_NUM_CONST_VARS_VTU_OUTPUT;

/**@brief evolution variable IDs written to vtu files*/
extern unsigned int BSSN_VTU_OUTPUT_EVOL_INDICES[BSSN_NUM_VARS];

/**@brief constraint variable IDs written to vtu files*/
extern unsigned int BSSN_VTU_OUTPUT_CONST_INDICES[BSSN_CONSTRAINT_NUM_VARS];

/**@brief solution output gap (instead of freq. we can use to output the
 * solution if currentTime > lastIOOutputTime + BSSN_IO_OUTPUT_GAP)*/
extern double BSSN_IO_OUTPUT_GAP;

/**@brief prefered grain sz to use when selecting active npes*/
extern unsigned int BSSN_DENDRO_GRAIN_SZ;

/**@brief AMR coarsening factor (we coarsen if
 * tol<BSSN_DENDRO_AMR_FAC*BSSN_WAVELET_TOL)*/
extern double BSSN_DENDRO_AMR_FAC;

/**@brief: AMR radius for the BH location based refinement. (BH1)*/
extern double BSSN_BH1_AMR_R;

/**@brief: AMR radius for the BH location based refinement. (BH2) */
extern double BSSN_BH2_AMR_R;

/** @brief: AMR ratio for the far factor, this was originally just 2.5, which is
 * why it's the default **/
extern double BSSN_AMR_R_RATIO;

/**@brief: BH Mass. (BH1)*/
extern double BSSN_BH1_MASS;

/**@brief: BH Mass. (BH2)*/
extern double BSSN_BH2_MASS;

/**@brief: skip grid point distance from bh < BSSN_BH1_CONSTRAINT_R when
 * computing the constraint norms for BH1*/
extern double BSSN_BH1_CONSTRAINT_R;

/**@brief: skip grid point distance from bh < BSSN_BH2_CONSTRAINT_R when
 * computing the constraint norms  for BH2*/
extern double BSSN_BH2_CONSTRAINT_R;

/**@brief: Maximum refinement level for BH1 */
extern unsigned int BSSN_BH1_MAX_LEV;

/**@brief: Maximum refinement level for BH2 */
extern unsigned int BSSN_BH2_MAX_LEV;

/**@brief: grid iterations till the grid converges*/
extern unsigned int BSSN_INIT_GRID_ITER;

/**@brief wavelet tolerance value. */
extern double BSSN_WAVELET_TOL;

/**@brief: wabelet tolernace for GW extration(refinement after the merger)*/
extern double BSSN_GW_REFINE_WTOL;

/**@brief load-imbalance tolerance value. */
extern double BSSN_LOAD_IMB_TOL;
/**@brief: Splitter fix value*/
extern unsigned int BSSN_SPLIT_FIX;

/**@brief: async. communication at a time. (upper bound shoud be BSSN_NUM_VARS)
 */
extern unsigned int BSSN_ASYNC_COMM_K;

/**@brief simulation begin time. */
extern double BSSN_RK_TIME_BEGIN;
/**@brief simulation end time*/
extern double BSSN_RK_TIME_END;
/**@brief rk time step size. */
extern double BSSN_RK45_TIME_STEP_SIZE;

/**@brief rk method type*/
extern unsigned int BSSN_RK_TYPE;

/** desired tolerance value for the rk45 method (adaptive time stepping. )*/
extern double BSSN_RK45_DESIRED_TOL;

/**@brief Black hole 1 */
extern BH BH1;
/**@brief Black hole 2 */
extern BH BH2;

/**@brief: evolved locations of the 0-BH1, 1 BH2*/
extern Point BSSN_BH_LOC[2];

/**@brief BBH initial data type */
extern unsigned int BSSN_ID_TYPE;

/**@brief physical coordinates for grid, x_min */
extern double BSSN_GRID_MIN_X;

/**@brief physical coordinates for grid, x_max */
extern double BSSN_GRID_MAX_X;

/**@brief physical coordinates for grid, y_min */
extern double BSSN_GRID_MIN_Y;

/**@brief physical coordinates for grid, y_max */
extern double BSSN_GRID_MAX_Y;

/**@brief physical coordinates for grid, z_min */
extern double BSSN_GRID_MIN_Z;

/**@brief physical coordinates for grid, z_max */
extern double BSSN_GRID_MAX_Z;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double BSSN_BLK_MIN_X;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double BSSN_BLK_MIN_Y;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double BSSN_BLK_MIN_Z;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double BSSN_BLK_MAX_X;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double BSSN_BLK_MAX_Y;

/**@brief physical coordinates for the blk adaptive x_min*/
extern double BSSN_BLK_MAX_Z;

/**@brief: dimension of the grid*/
extern unsigned int BSSN_DIM;

/**@brief: minimum refinement level default=2*/
extern unsigned int BSSN_MINDEPTH;

/**@brief: max refinement level*/
extern unsigned int BSSN_MAXDEPTH;

/**@brief: lambda values for evolution */
extern unsigned int BSSN_LAMBDA[4];

/**@brief: lambda values for evolution */
extern double BSSN_LAMBDA_F[2];

/**@brief : parameters for eta_damping function */
extern double BSSN_ETA_R0;
extern double BSSN_ETA_POWER[2];

/**@brief : parameters for RIT eta_damping function */
extern unsigned int RIT_ETA_FUNCTION;
extern double RIT_ETA_OUTER;
extern double RIT_ETA_CENTRAL;
extern double RIT_ETA_WIDTH;

/**@brief: lambda values for evolution */
extern double BSSN_TRK0;

/**@brief: base value for eta in evolution */
extern double ETA_CONST;

/**@brief: eta_R0, radius where eta is damped for evolution */
extern double ETA_R0;

/**@brief: eta damping for evolution */
extern double ETA_DAMPING;

/**@brief: eta damping exponent for evolution */
extern double ETA_DAMPING_EXP;

/**@brief: chi floor value */
extern double CHI_FLOOR;

/**@brief: dissipation type */
extern unsigned int DISSIPATION_TYPE;

/**@brief: Kreiss-Oliger dissipation */
extern double KO_DISS_SIGMA;

/**@brief: BSSN_USE_WAVELET_TOL_FUNCTION */
extern unsigned int BSSN_USE_WAVELET_TOL_FUNCTION;

/**@brief: BSSN_WAVELET_TOL_FUNCTION_R0 */
extern double BSSN_WAVELET_TOL_FUNCTION_R0;

/**@brief: BSSN_WAVELET_TOL_FUNCTION_R0 */
extern double BSSN_WAVELET_TOL_FUNCTION_R1;

/**@brief: BSSN_WAVELET_TOL_MAX */
extern double BSSN_WAVELET_TOL_MAX;

/**@brief: eta function parameters*/
extern double BSSN_ETA_R0;

/**@brief: eta function parameters (powers)*/
extern double BSSN_ETA_POWER[2];

/**@brief: xi parameters to change between different gauge conditions*/
extern unsigned int BSSN_XI[3];

/**@brief: @david can you please add some comeents for these parameters. */
extern unsigned int BSSN_DISSIPATION_NC;

/**@brief: @david can you please add some comeents for these parameters. */
extern unsigned int BSSN_DISSIPATION_S;

/**@brief: tolerance for refinement based on EH */
extern double BSSN_EH_REFINE_VAL;

/**@brief: tolerance for coarsen based on EH */
extern double BSSN_EH_COARSEN_VAL;

/**@brief: refinement mode for the application*/
extern RefinementMode BSSN_REFINEMENT_MODE;

/**@brief: option to enable if the set refinement mode should be used for
 * initial grid converge */
extern bool BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE;

/**@brief: if true output only the z slice*/
extern bool BSSN_VTU_Z_SLICE_ONLY;

/**@brief TS off set for LTS in BSSN*/
extern unsigned int BSSN_LTS_TS_OFFSET;

/**@brief to track if a merged checkpoint file is writte.*/
extern bool BSSN_MERGED_CHKPT_WRITTEN;

/**@brief BSSN RK coordinate time*/
extern double BSSN_CURRENT_RK_COORD_TIME;

/**@brief BSSN RK step*/
extern unsigned int BSSN_CURRENT_RK_STEP;

/***@brief: derivs workspace*/
extern double* BSSN_DERIV_WORKSPACE;

/** @brief: Nyquist goal for nyquist-based refinement **/
extern unsigned int BSSN_NYQUIST_M;

extern bool BSSN_SCALE_VTU_AND_GW_EXTRACTION;

extern unsigned int BSSN_GW_EXTRACT_FREQ_TRUE;

extern unsigned int BSSN_IO_OUTPUT_FREQ_TRUE;

extern double BSSN_SSL_SIGMA;
extern double BSSN_SSL_H;

// note ko derivs are not included
#ifdef BSSN_USE_ADVECTIVE_DERIVS
const unsigned int BSSN_NUM_DERIVS = 138 + 74;
#else
const unsigned int BSSN_NUM_DERIVS = 138;
#endif

void readParamTOMLFile(const char* fName, MPI_Comm comm);

}  // namespace bssn

namespace TPID {
static const double TP_epsilon              = 1.0e-6;
static const int swap_xz                    = 0;
static const int use_sources                = 0;
static const int rescale_sources            = 0;
static const int use_external_initial_guess = 0;
static const int do_residuum_debug_output   = 1;
static const int do_initial_debug_output    = 1;
static const int multiply_old_lapse         = 0;
static const double TP_Tiny                 = 1.0e-15;
static const double TP_Extend_Radius        = 0.0;
static const int Newton_maxit               = 5;

extern double target_M_plus;
extern double target_M_minus;
extern double par_m_plus;
extern double par_m_minus;
extern double par_b;
extern double par_P_plus[3];
extern double par_P_minus[3];
extern double par_S_plus[3];
extern double par_S_minus[3];
extern double center_offset[3];
extern double initial_lapse_psi_exponent;
extern int npoints_A;
extern int npoints_B;
extern int npoints_phi;
extern int give_bare_mass;
extern int initial_lapse;
extern int solve_momentum_constraint;
extern int grid_setup_method;
extern int verbose;
extern double adm_tol;
extern double Newton_tol;
extern std::string FILE_PREFIX;
extern bool replace_lapse_with_sqrt_chi;
}  // namespace TPID

/**@brief parameters related to BH location extraction*/
namespace BHLOC {
// Note: Current implementation of the BH location extraction is simple
// clustering property based on the initial location of the BH Later we can add
// the code to compute the BH locations based on the event horizon extraction.

/**@brief variable ID used for BH location extraction*/
extern unsigned int EXTRACTION_VAR_ID;

/**@brief tolerance for BH extraction*/
extern double EXTRACTION_TOL;

}  // namespace BHLOC

/**@brief parameters related to GW extraction*/
namespace GW {

/**@brief max allowed radii values*/
static const unsigned int BSSN_GW_MAX_RADAII       = 20;

/**@brief max allowed lmode values*/
static const unsigned int BSSN_GW_MAX_LMODES       = 8;

static const unsigned int BSSN_GW_LEBEDEV_PREC     = 25;

/**@brief: GW output precision in .dat files*/
static const unsigned int BSSN_GW_OUTPUT_PRECISION = 10;

/**@brief: number of different radius values for psi4 poly fit*/
extern unsigned int BSSN_GW_NUM_RADAII;

/**@brief: number of L mode values*/
extern unsigned int BSSN_GW_NUM_LMODES;

/**@brief: vallues of extraction radaii*/
extern double BSSN_GW_RADAII[BSSN_GW_MAX_RADAII];

/**@brief: value of the spin*/
extern unsigned int BSSN_GW_SPIN;

/**@brief values for l modes in SWSH*/
extern unsigned int BSSN_GW_L_MODES[BSSN_GW_MAX_LMODES];

}  // namespace GW

namespace AEH {
/**@brief lmax used for AH surface parameterization*/
extern unsigned int AEH_LMAX;

/**@brief quadrature points in the theta direction*/
extern unsigned int AEH_Q_THETA;

/**@brief quadrature points in the phi direction*/
extern unsigned int AEH_Q_PHI;

/**@brief number of max. iterations for AH solver*/
extern unsigned int AEH_MAXITER;

/**@brief absolute tolerance for AH convergence*/
extern double AEH_ATOL;

/**@brief relative tolerance for AH convergence*/
extern double AEH_RTOL;

extern unsigned int AEH_SOLVER_FREQ;

/***@brief AEH update A = (alpha / (lmax * (lmax + 1)))  + beta  and B=
 * beta/alpha and dlambda = A / (1 + B * l * (l + 1) ) */
extern double AEH_ALPHA;

/***@brief AEH update A = (alpha / (lmax * (lmax + 1)))  + beta  and B=
 * beta/alpha and dlambda = A / (1 + B * l * (l + 1) ) */
extern double AEH_BETA;

}  // namespace AEH
