/**
 * @file bssnCtx.h
 * @author Milinda Fernando
 * @brief Application context class for solving the Einstein equations in BSSNKO
 * formulation.
 * @version 0.1
 * @date 2019-12-20
 *
 * @copyright Copyright (c) 2019, University of Utah.
 *
 */

#pragma once
#include <filesystem>
#include <iostream>

#include "TwoPunctures.h"
#include "aeh.h"
#include "bssn_constraints.h"
#include "checkPoint.h"
#include "ctx.h"
#include "dataUtils.h"
#include "grDef.h"
#include "grUtils.h"
#include "gwExtract.h"
#include "mathMeshUtils.h"
#include "oct2vtk.h"
#include "parUtils.h"
#include "parameters.h"
#include "physcon.h"
#include "rhs.h"

namespace bssn {

/**@brief smoothing modes avail for LTS recomended for LTS time stepping. */
enum LTS_SMOOTH_MODE { KO = 0, WEIGHT_FUNC };
enum VL {
    CPU_EV = 0,
    CPU_CV,
    CPU_EV_UZ_IN,
    CPU_EV_UZ_OUT,
    CPU_CV_UZ_IN,
    GPU_EV,
    GPU_EV_UZ_IN,
    GPU_EV_UZ_OUT,
    END
};
typedef ot::DVector<DendroScalar, unsigned int> DVec;

class BSSNCtx : public ts::Ctx<BSSNCtx, DendroScalar, unsigned int> {
   protected:
    /**@brief: evolution var (zip)*/
    DVec m_var[VL::END];

    Point m_uiBHLoc[2];

   private:
    // TODO: move these into the main Ctx object and have the remesh logic work
    // there

    DendroIntL m_uiGlobalMeshElements;
    DendroIntL m_uiGlobalGridPoints;
    bool m_uiWroteGridInfoHeader = false;

   public:
    /**@brief: default constructor*/
    BSSNCtx(ot::Mesh* pMesh);

    /**@brief: default deconstructor*/
    ~BSSNCtx();

    /**@brief get bh locations*/
    const Point& get_bh0_loc() const { return m_uiBHLoc[0]; }
    const Point& get_bh1_loc() const { return m_uiBHLoc[1]; }

    /**
     * @brief sets time adaptive offset
     * @param tadapoffst
     */
    void set_time_adap_offset(unsigned int tadapoffst) {
        BSSN_LTS_TS_OFFSET = tadapoffst;
    }

    /**@brief : returns the time adaptive offset value*/
    unsigned int get_time_adap_offset() { return BSSN_LTS_TS_OFFSET; }

    /**@brief: initial solution and grid convergence calls init_grid()*/
    int initialize();

    /**@brief: initialize the grid, solution. */
    int init_grid();

    /**
     * @brief computes the BSSN rhs
     *
     * @param in : zipped input
     * @param out : zipped output
     * @param sz  : number of variables.
     * @param time : current time.
     * @return int : status. (0) on success.
     */
    int rhs(DVec* in, DVec* out, unsigned int sz, DendroScalar time);

    /**
     * @brief computes AEH expansion function on the spacelike hypersurface.
     *
     * @param in : zipped input
     * @param out : zipped output
     * @param sz  : number of variables.
     * @param time : current time.
     * @return int : status. (0) on success.
     */
    int aeh_expansion(const Point& origin, aeh::AEH_VARS* m_aeh_vars,
                      DVec& aeh_f, DVec& aeh_h, const DendroScalar* const rlim);

    /**
     * @brief block wise RHS.
     *
     * @param in : input vector (unzip version)
     * @param out : output vector (unzip version)
     * @param blkIDs : local block ids where the rhs is computed.
     * @param sz : size of the block ids
     * @param blk_time : block time  corresponding to the block ids.
     * @return int
     */
    int rhs_blkwise(DVec in, DVec out, const unsigned int* const blkIDs,
                    unsigned int numIds, DendroScalar* blk_time);

    int rhs_blk(const DendroScalar* in, DendroScalar* out, unsigned int dof,
                unsigned int local_blk_id, DendroScalar blk_time);

    int pre_stage_blk(DendroScalar* in, unsigned int dof,
                      unsigned int local_blk_id, DendroScalar blk_time);

    int post_stage_blk(DendroScalar* in, unsigned int dof,
                       unsigned int local_blk_id, DendroScalar blk_time);

    int pre_timestep_blk(DendroScalar* in, unsigned int dof,
                         unsigned int local_blk_id, DendroScalar blk_time);

    int post_timestep_blk(DendroScalar* in, unsigned int dof,
                          unsigned int local_blk_id, DendroScalar blk_time);

    /**@brief: function execute before each stage
     * @param sIn: stage var in.
     */
    inline int pre_stage(DVec& sIn) { return 0; }

    /**@brief: function execute after each stage
     * @param sIn: stage var in.
     */
    int post_stage(DVec& sIn);

    /**@brief: function execute before each step*/
    inline int pre_timestep(DVec& sIn) { return 0; }

    /**@brief: function execute after each step*/
    int post_timestep(DVec& sIn);

    /**@brief: function execute after each step*/
    bool is_remesh();

    /**@brief: write to vtu. */
    int write_vtu();

    /**@brief: writes checkpoint*/
    int write_checkpt();

    /**@brief: restore from check point*/
    int restore_checkpt();

    /**@brief: should be called for free up the contex memory. */
    int finalize();

    void compute_constraint_variables();

    /**@brief: pack and returns the evolution variables to one DVector*/
    DVec& get_evolution_vars();

    /**@brief: pack and returns the constraint variables to one DVector*/
    DVec& get_constraint_vars();

    /**@brief: pack and returns the primitive variables to one DVector*/
    DVec& get_primitive_vars();

    /**@brief: prints any messages to the terminal output. */
    int terminal_output();

    /**
     * @brief: writes information about the simulation to a .dat file
     *
     * Some of the information saved here:
     *   - Step
     *   - Wall Time from last step
     *   - Simulation time ( t / M)
     *   - Active communicator size
     *   - Number of mesh elements
     *   - Current dt
     */
    void write_grid_summary_data();

    /**@brief: returns the async communication batch size. */
    unsigned int get_async_batch_sz() { return bssn::BSSN_ASYNC_COMM_K; }

    /**@brief: returns the number of variables considered when performing
     * refinement*/
    unsigned int get_num_refine_vars() { return BSSN_NUM_REFINE_VARS; }

    /**@brief: return the pointer for containing evolution refinement variable
     * ids*/
    const unsigned int* get_refine_var_ids() {
        return BSSN_REFINE_VARIABLE_INDICES;
    }

    /**@brief return the wavelet tolerance function / value*/
    std::function<double(double, double, double, double* hx)>
    get_wtol_function() {
        double wtol = BSSN_WAVELET_TOL;
        std::function<double(double, double, double, double*)> waveletTolFunc =
            [](double x, double y, double z, double* hx) {
                return bssn::computeWTolDCoords(x, y, z, hx);
            };
        return waveletTolFunc;
    }

    /**@brief computes the LTS TS offset based on the eta damping.*/
    unsigned int compute_lts_ts_offset();

    /**@brief : blk time step factor. */
    static unsigned int getBlkTimestepFac(unsigned int blev, unsigned int lmin,
                                          unsigned int lmax);

    /**@biref: evolve bh locations. */
    void evolve_bh_loc(DVec sIn, double dt);

    /**
     * @brief LTS smooth mode.
     *
     * @param sIn : time synced evolution vector.
     * @param mode : smoothing mode.
     */
    void lts_smooth(DVec sIn, LTS_SMOOTH_MODE mode);

    /**@brief: return true if the BH are merged. */
    bool is_bh_merged(double tol) const {
        return ((bssn::BSSN_BH_LOC[0] - bssn::BSSN_BH_LOC[1]).abs() < tol);
    };

    int grid_transfer(const ot::Mesh* m_new);

    void calculate_full_grid_size();
};

}  // end of namespace bssn
