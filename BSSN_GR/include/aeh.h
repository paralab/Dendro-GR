/**
 * @brief main apparent event horizon solvers code
 * @version 1.0
 * @date 2023-06-23
 * @author Milinda Fernando (milinda@oden.utexas.edu)
 * @copyright Copyright (c) 2019
 */

#pragma once

#include <fstream>

#include "daUtils.h"
#include "derivs.h"
#include "dvec.h"
#include "grDef.h"
#include "mesh.h"
#include "mpi.h"
#include "parameters.h"
#include "sph.h"

namespace aeh {
enum AEHSolverType { NKO = 0, FAST_FLOW };
enum AEHErrorType { SUCCESS, MAX_ITERATIONS_REACHED, FAIL };
typedef std::pair<int, int> LM_MODE;
typedef ot::DVector<DendroScalar, unsigned int> DVec;

typedef struct {
    DVec grad_gt;
    DVec grad_chi;

    DVec gt;
    DVec At;
    DVec chi;
    DVec K;

} AEH_VARS;

typedef struct {
    DendroScalar area;
    DendroScalar Jx;
    DendroScalar Jy;
    DendroScalar Jz;
    DendroScalar rmin;
    DendroScalar rmean;
    DendroScalar rmax;
    DendroScalar expansion_L2;

} AEH_QoI;

template <typename Ctx, typename T>
class SpectralAEHSolver {
   private:
    /**@brief pointer to the GR context class*/
    Ctx* m_ctx          = nullptr;

    /***@breif maximum l allowed.*/
    unsigned int m_lmax = 0;

    /**@breif currently using lmax value */
    unsigned int m_ll   = 0;

    /**@brief spherical modes used for the AEH solver*/
    std::vector<LM_MODE> m_sph_modes;

    /**@brief number of theta (polar angle) points used for quadrature*/
    unsigned int m_num_theta;

    /**@brief number of phi (azimuthal angle) points used for quadrature*/
    unsigned int m_num_phi;

    /**@brief gls workspace for Gauss quadrature*/
    gsl_integration_fixed_workspace* m_quad_theta  = nullptr;

    /**@brief gls workspace for Gauss quadrature without sin(theta) correction*/
    gsl_integration_fixed_workspace* m_quad_theta1 = nullptr;

    /**@brief gls workspace for Gauss quadrature*/
    gsl_integration_fixed_workspace* m_quad_phi    = nullptr;

    /**@brief AEH setup vars to evaluate the expansion*/
    AEH_VARS* m_aeh_vars;

    std::vector<T> m_sin_qtheta;
    std::vector<T> m_cos_qtheta;

    std::vector<T> m_sin_qphi;
    std::vector<T> m_cos_qphi;

    std::vector<DendroIntL> m_active_ele_local_id;
    std::vector<DendroIntL> m_active_blk_local_id;
    DendroScalar m_rlim[2];

    Point m_origin                      = Point();

    AEHSolverType m_solver_type         = AEHSolverType::FAST_FLOW;

    /**@brief number of quasi measure variables.
     * 0 - det(induced metric)
     * 1 - J_x (angular mom. x )
     * 2 - J_y (angular mom. y )
     * 3 - J_z (angular mom. z )
     */
    unsigned int NUM_QUASI_MEASURE_VARS = 4;

   public:
    /**
     * @brief Construct a new Spectral AEH solver object
     *
     * @param ctx   : Application Ctx, i.e., GR spatial variables.
     * @param l_max : Spectral expansion using l=0 to l=l_max, and corresponding
     * m-modes
     * @param azimuthal_symmetry : Use azimuthal symmetry or not. (m=0) modes
     */
    SpectralAEHSolver(Point origin, Ctx* ctx, unsigned int l_max,
                      unsigned int q_theta, unsigned int q_phi,
                      const DendroScalar* const rlim,
                      bool azimuthal_symmetry = false, bool verbose = true);

    /**@brief deconstructor for AEH solver object*/
    ~SpectralAEHSolver();

    /**@breif set the current active lmodes */
    void set_lmodes(unsigned int ll) { m_ll = ll; }

    /**@brief get the current l modes */
    unsigned int get_lmodes() { return m_ll; }

    /**@brief set aeh solver type*/
    void set_solver_type(AEHSolverType type) { m_solver_type = type; }

    /**@brief get aeh solver type*/
    AEHSolverType get_solver_type() { return m_solver_type; }

    const std::vector<LM_MODE>& get_lm_modes() const { return m_sph_modes; }

    unsigned int get_num_lm_modes() const { return m_sph_modes.size(); }

    /**
     * @brief Level set function computation for the AEH surface.
     * @param m_uiMesh
     * @param aeh_f_ptr
     * @param h_qs
     * @return int
     */
    int eval_aeh_level_set(ot::Mesh* m_uiMesh, T* const aeh_f_ptr,
                           const T* const h_qs);

    /**
     * @brief computes the induced metric on the AEH surface,
     * @param m_uiMesh
     * @param ctx
     * @param aeh_f_ptr
     * @param h_qs
     * @return int
     */
    int eval_quasi_measurements_integrands(Ctx* ctx, ot::Mesh* m_uiMesh,
                                           T* const u_ptr, const T* const h_qs);

    /**
     * @brief main AEH solver update
     * @param ctx
     * @param r_init
     * @param max_iter
     * @return AEHErrorType
     */
    AEHErrorType solve(Ctx* ctx, const T* const h_init, T* h_qs,
                       unsigned int max_iter, DendroScalar rel_eps,
                       DendroScalar abs_eps, double alpha, double beta,
                       unsigned int verbose = 0);

    /**
     * @brief evaluate the expansion norm
     * @param origin
     * @param ctx
     * @param h_qs
     * @param aeh_f : workspace var for level set evaluation
     * @param aeh_h : workspace var for expansion evaluation
     * @param interp_coords : workspace for interpolation coordinates
     * @param norm_type
     * @return DendroScalar
     */
    DendroScalar eval_expansion_norm(Ctx* ctx, const T* const h_qs, DVec& aeh_f,
                                     DVec& aeh_h, std::vector<T>& interp_coords,
                                     unsigned int norm_type);

    /**
     * @brief
     *
     * @param origin
     * @param ctx
     * @param h_qs
     * @param qoi
     * @param ws_v1
     */
    void aeh_quasi_measurements(Ctx* ctx, const T* const h_qs, AEH_QoI& qoi);

    /**
     * @brief write the horizon to json
     * @param origin coordinate of the AEH parameterization
     * @param ctx gr ctx pointer
     * @param h_qs AEH expansion coefficients
     * @param fname file name
     * @param mode file write mode
     */
    int aeh_to_json(Ctx* ctx, T* h_qs, char* fname,
                    std::ios_base::openmode mode);

   private:
    DendroScalar rhs_00(Ctx* ctx, DendroScalar a, T* h_qs, DVec& aeh_f,
                        DVec& aeh_h, std::vector<T>& interp_coords) {
        ot::Mesh* m_uiMesh = ctx->get_mesh();
        if (!(m_uiMesh->isActive())) return -1;

        const Point& origin = m_origin;
        T* aeh_f_ptr        = aeh_f.get_vec_ptr();
        T* aeh_h_ptr        = aeh_h.get_vec_ptr();

        h_qs[0]             = a;
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

        this->eval_aeh_level_set(m_uiMesh, aeh_f_ptr, h_qs);
        ctx->aeh_expansion(origin, m_aeh_vars, aeh_f, aeh_h, m_rlim);

        m_uiMesh->readFromGhostBegin(aeh_h.get_vec_ptr(), 1);

        unsigned int num_angular_pts = m_num_theta * m_num_phi;
        const double* const m_qtheta = m_quad_theta->x;
        const double* const m_qphi   = m_quad_phi->x;

        for (unsigned int qt = 0; qt < m_num_theta; qt++)
            for (unsigned int qp = 0; qp < m_num_phi; qp++) {
                const unsigned int q_idx = qt * m_num_phi + qp;

                T r_val                  = (T)0;
                for (unsigned int lm_idx = 0; lm_idx < m_sph_modes.size();
                     lm_idx++) {
                    int l = m_sph_modes[lm_idx].first;
                    if (l > m_ll) continue;

                    r_val += real_spherical_harmonic(m_sph_modes[lm_idx].first,
                                                     m_sph_modes[lm_idx].second,
                                                     m_qtheta[qt], m_qphi[qp]) *
                             h_qs[lm_idx];
                }

                interp_coords[3 * q_idx + 0] =
                    origin.x() + (r_val)*m_sin_qtheta[qt] * m_cos_qphi[qp];
                interp_coords[3 * q_idx + 1] =
                    origin.y() + (r_val)*m_sin_qtheta[qt] * m_sin_qphi[qp];
                interp_coords[3 * q_idx + 2] =
                    origin.z() + (r_val)*m_cos_qtheta[qt];
            }

        m_uiMesh->readFromGhostEnd(aeh_h.get_vec_ptr(), 1);

        std::vector<unsigned int> valid_idx;
        valid_idx.clear();

        std::vector<T> aeh_h_inp;
        aeh_h_inp.resize(num_angular_pts);

        ot::da::interpolateToCoords(m_uiMesh, aeh_h_ptr, interp_coords.data(),
                                    interp_coords.size(), grid_limits,
                                    domain_limits, aeh_h_inp.data(), valid_idx);
        // printf("rhs_{00}\n");
        // printArray_1D(interp_coords.data(), interp_coords.size());
        // printArray_1D(aeh_h_inp.data(), aeh_h_inp.size());

        DendroScalar result   = 0;
        DendroScalar result_g = 0;
        for (unsigned int idx = 0; idx < valid_idx.size(); idx++) {
            const unsigned int qp = valid_idx[idx] % m_num_phi;
            const unsigned int qt = (valid_idx[idx] - qp) / m_num_phi;

            DendroScalar s =
                aeh_h_inp[valid_idx[idx]] * aeh_h_inp[valid_idx[idx]];
            // for(unsigned int lm_idx=0; lm_idx < m_sph_modes.size(); lm_idx++)
            // {
            //     int l  = m_sph_modes[lm_idx].first;
            //     int m  = m_sph_modes[lm_idx].second;

            //     s += -(l)*(l+1) * real_spherical_harmonic(l, m, m_qtheta[qt],
            //     m_qphi[qp]) * h_qs[lm_idx];

            // }
            result += (s)*m_quad_theta->weights[qt] * m_quad_phi->weights[qp];
        }

        par::Mpi_Allreduce(&result, &result_g, 1, MPI_SUM,
                           m_uiMesh->getMPICommunicator());
        // par::Mpi_Allreduce(&result, &result_g, 1, MPI_MAX,
        // m_uiMesh->getMPICommunicator());
        return result_g;
    }

    DendroScalar solve_00_bisection(DendroScalar a, DendroScalar b, Ctx* ctx,
                                    T* h_qs, DVec& aeh_f, DVec& aeh_h,
                                    std::vector<T>& interp_coords,
                                    DendroScalar eps) {
        ot::Mesh* m_uiMesh = ctx->get_mesh();
        if (!(m_uiMesh->isActive())) return -1.0;

        const Point& origin = m_origin;
        int rank            = m_uiMesh->getMPIRank();
        int npes            = m_uiMesh->getMPICommSize();

        DendroScalar h_00   = h_qs[0];

        DendroScalar f_a =
            rhs_00(origin, ctx, a, h_qs, aeh_f, aeh_h, interp_coords);
        DendroScalar f_b =
            rhs_00(origin, ctx, b, h_qs, aeh_f, aeh_h, interp_coords);

        if (f_a * f_b >= 0) {
            if (!rank)
                printf(
                    "lm=(0,0) initial interval (a,b) = (%.8E, %.8E) = (f(a), "
                    "f(b)) = (%.8E, %.8E)\n",
                    a, b, f_a, f_b);

            return -1.0;
        }

        DendroScalar c = a;
        while ((b - a) >= eps) {
            // Find middle point
            c = (a + b) / 2;
            DendroScalar f_c =
                rhs_00(origin, ctx, c, h_qs, aeh_f, aeh_h, interp_coords);

            if (!rank)
                printf("lm=(0,0)  h_00 = %.8E ,  abs(H(h)) = %.8E\n", c,
                       abs(f_c));

            // Check if middle point is root
            if (abs(f_c) < eps) {
                return c;
            }
            DendroScalar f_a =
                rhs_00(origin, ctx, a, h_qs, aeh_f, aeh_h, interp_coords);
            // Decide the side to repeat the steps
            if (f_c * f_a < 0)
                b = c;
            else
                a = c;
        }

        return c;
    }

    DendroScalar solve_00_newton(Ctx* ctx, T* h_qs, DVec& aeh_f, DVec& aeh_h,
                                 std::vector<T>& interp_coords,
                                 DendroScalar eps) {
        ot::Mesh* m_uiMesh = ctx->get_mesh();
        if (!(m_uiMesh->isActive())) return -1.0;

        const Point& origin         = m_origin;
        int rank                    = m_uiMesh->getMPIRank();
        int npes                    = m_uiMesh->getMPICommSize();

        DendroScalar h_00           = h_qs[0];
        DendroScalar x              = h_00;

        DendroScalar relative_error = 1.0;
        unsigned int iter           = 0;
        unsigned int max_iter       = 100;
        DendroScalar grad_f         = 0.0;
        DendroScalar fa = rhs_00(ctx, x, h_qs, aeh_f, aeh_h, interp_coords);
        do {
            DendroScalar dh =
                sqrt(std::numeric_limits<DendroScalar>::epsilon()) * x;
            DendroScalar fp =
                rhs_00(ctx, x + dh, h_qs, aeh_f, aeh_h, interp_coords);
            DendroScalar fm =
                rhs_00(ctx, x - dh, h_qs, aeh_f, aeh_h, interp_coords);

            grad_f             = 0.5 * (fp - fm) / dh;
            DendroScalar alpha = 1e0;

            DendroScalar b, fb;

            if (!rank)
                printf(
                    "  Newton iter=%03d\t h_00=%.8E\t abs(f)=%.8E\t "
                    "relative_error=%.8E |grad_f|=%.8E\n",
                    iter, x, abs(fa), relative_error, abs(grad_f));

            if ((abs(grad_f) < eps) || (abs(fa) < eps)) break;

            do {
                b = x - alpha * (fa / grad_f);
                if (b < m_rlim[0] || b > m_rlim[1])
                    fb = 1e10;
                else
                    fb = rhs_00(ctx, abs(b), h_qs, aeh_f, aeh_h, interp_coords);

                alpha = alpha * 0.1;

            } while ((alpha > 1e-16) && ((abs(fb) > abs(fa))));

            relative_error = abs(b - x) / abs(x);

            iter += 1;
            x  = b;
            fa = fb;

        } while ((iter < max_iter) && (relative_error > eps) &&
                 (abs(grad_f) > eps));

        return x;
    }
};

template <typename Ctx, typename T>
SpectralAEHSolver<Ctx, T>::SpectralAEHSolver(
    Point origin, Ctx* ctx, unsigned int l_max, unsigned int q_theta,
    unsigned int q_phi, const DendroScalar* const rlim, bool azimuthal_symmetry,
    bool verbose) {
    m_ctx              = ctx;
    ot::Mesh* m_uiMesh = ctx->get_mesh();

    if (!m_uiMesh->isActive()) return;

    const int rank             = m_uiMesh->getMPIRank();
    m_lmax                     = l_max;
    m_origin                   = origin;

    m_rlim[0]                  = rlim[0];
    m_rlim[1]                  = rlim[1];

    const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin()));
    for (unsigned int elem = m_uiMesh->getElementLocalBegin();
         elem < m_uiMesh->getElementLocalEnd(); elem++) {
        const DendroScalar lx =
            (Rx / RgX) *
            (double)(1u << (m_uiMaxDepth - pNodes[elem].getLevel()));
        const DendroScalar ly =
            (Ry / RgY) *
            (double)(1u << (m_uiMaxDepth - pNodes[elem].getLevel()));
        const DendroScalar lz =
            (Rz / RgZ) *
            (double)(1u << (m_uiMaxDepth - pNodes[elem].getLevel()));

        const DendroScalar ll   = sqrt(lx * lx + ly * ly + lz * lz);

        const DendroScalar x    = pNodes[elem].getX();
        const DendroScalar y    = pNodes[elem].getY();
        const DendroScalar z    = pNodes[elem].getZ();

        const DendroScalar xx   = GRIDX_TO_X(x) - origin.x();
        const DendroScalar yy   = GRIDY_TO_Y(y) - origin.y();
        const DendroScalar zz   = GRIDZ_TO_Z(z) - origin.z();

        const DendroScalar p_rr = sqrt(xx * xx + yy * yy + zz * zz);

        if ((p_rr < (m_rlim[1] + ll)) && (p_rr > (m_rlim[0] - ll)))
            m_active_ele_local_id.push_back(elem);
    }

    // printf("rank = %04d active elements =%d\n", rank,
    // m_active_ele_local_id.size());

    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW        = bssn::BSSN_PADDING_WIDTH;

    const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        const unsigned int offset = blkList[blk].getOffset();
        const unsigned int sz[3]  = {blkList[blk].getAllocationSzX(),
                                     blkList[blk].getAllocationSzY(),
                                     blkList[blk].getAllocationSzZ()};
        const unsigned int bflag  = blkList[blk].getBlkNodeFlag();

        const DendroScalar dx     = blkList[blk].computeDx(pt_min, pt_max);
        const DendroScalar dy     = blkList[blk].computeDy(pt_min, pt_max);
        const DendroScalar dz     = blkList[blk].computeDz(pt_min, pt_max);

        DendroScalar ptmin[3], ptmax[3];

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

        const DendroScalar xx   = ptmin[0] - origin.x();
        const DendroScalar yy   = ptmin[1] - origin.y();
        const DendroScalar zz   = ptmin[2] - origin.z();

        const DendroScalar lx   = (ptmax[0] - ptmin[0]);
        const DendroScalar ly   = (ptmax[1] - ptmin[1]);
        const DendroScalar lz   = (ptmax[2] - ptmin[2]);

        const DendroScalar ll   = sqrt(lx * lx + ly * ly + lz * lz);

        const DendroScalar p_rr = sqrt(xx * xx + yy * yy + zz * zz);

        if ((p_rr < (m_rlim[1] + ll)) && (p_rr > (m_rlim[0] - ll)))
            m_active_blk_local_id.push_back(blk);
    }

    if (azimuthal_symmetry) {
        for (unsigned int l = 0; l <= l_max; l++)
            m_sph_modes.push_back(LM_MODE(l, 0));
    } else {
        for (int l = 0; l <= l_max; l++)
            for (int m = -l; m <= l; m++) m_sph_modes.push_back(LM_MODE(l, m));
    }

    // need to initialize the quadrature grid
    m_num_theta   = q_theta;
    m_num_phi     = q_phi;

    // need to initialize the q_theta based on the Gauss-Legendre quadrature
    m_quad_theta  = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
                                                m_num_theta, -1, 1, 0.0, 0.0);
    m_quad_phi    = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
                                                m_num_phi, 0, 2 * M_PI, 0.0, 0.0);
    m_quad_theta1 = gsl_integration_fixed_alloc(gsl_integration_fixed_legendre,
                                                m_num_phi, 0, M_PI, 0.0, 0.0);

    for (unsigned int i = 0; i < m_num_theta; i++)
        m_quad_theta->x[i] = acos(m_quad_theta->x[i]);

    m_sin_qtheta.resize(m_num_theta);
    m_cos_qtheta.resize(m_num_theta);

    m_sin_qphi.resize(m_num_phi);
    m_cos_qphi.resize(m_num_phi);

    const double* const m_qtheta = m_quad_theta->x;
    const double* const m_qphi   = m_quad_phi->x;

    for (unsigned int qt = 0; qt < m_num_theta; qt++) {
        m_sin_qtheta[qt] = sin(m_qtheta[qt]);
        m_cos_qtheta[qt] = cos(m_qtheta[qt]);
    }

    for (unsigned int qp = 0; qp < m_num_phi; qp++) {
        m_sin_qphi[qp] = sin(m_qphi[qp]);
        m_cos_qphi[qp] = cos(m_qphi[qp]);
    }

    // printArray_1D(m_quad_theta->x, m_num_theta);
    // printArray_1D(m_quad_theta->weights, m_num_theta);

    // printArray_1D(m_quad_phi->x, m_num_phi);
    // printArray_1D(m_quad_phi->weights, m_num_phi);

    m_aeh_vars = new AEH_VARS();

    m_aeh_vars->grad_gt.create_vector(m_uiMesh,
                                      ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                      ot::DVEC_LOC::HOST, 18, true);
    m_aeh_vars->grad_chi.create_vector(m_uiMesh,
                                       ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                       ot::DVEC_LOC::HOST, 3, true);

    m_aeh_vars->gt.create_vector(m_uiMesh,
                                 ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                 ot::DVEC_LOC::HOST, 6, true);
    m_aeh_vars->At.create_vector(m_uiMesh,
                                 ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                 ot::DVEC_LOC::HOST, 6, true);
    m_aeh_vars->chi.create_vector(m_uiMesh,
                                  ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                  ot::DVEC_LOC::HOST, 1, true);
    m_aeh_vars->K.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                ot::DVEC_LOC::HOST, 1, true);

    DVec& gr_vars             = ctx->get_evolution_vars();
    DendroScalar* gr_vars_ptr = gr_vars.get_vec_ptr();

    const unsigned int cg_sz  = m_uiMesh->getDegOfFreedom();
    const unsigned int uz_sz  = m_uiMesh->getDegOfFreedomUnZip();

    DendroScalar* gt          = &gr_vars_ptr[bssn::VAR::U_SYMGT0 * cg_sz];
    DendroScalar* At          = &gr_vars_ptr[bssn::VAR::U_SYMAT0 * cg_sz];
    DendroScalar* chi         = &gr_vars_ptr[bssn::VAR::U_CHI * cg_sz];
    DendroScalar* K           = &gr_vars_ptr[bssn::VAR::U_K * cg_sz];

    DendroScalar* gt_uz       = m_aeh_vars->gt.get_vec_ptr();
    DendroScalar* At_uz       = m_aeh_vars->At.get_vec_ptr();
    DendroScalar* chi_uz      = m_aeh_vars->chi.get_vec_ptr();
    DendroScalar* K_uz        = m_aeh_vars->K.get_vec_ptr();

    m_uiMesh->readFromGhostBegin<DendroScalar>(gt, 6);
    m_uiMesh->readFromGhostBegin<DendroScalar>(At, 6);
    m_uiMesh->readFromGhostBegin<DendroScalar>(chi, 1);
    m_uiMesh->readFromGhostBegin<DendroScalar>(K, 1);

    m_uiMesh->readFromGhostEnd<DendroScalar>(gt, 6);
    m_uiMesh->unzip(gt, gt_uz, 6);

    m_uiMesh->readFromGhostEnd<DendroScalar>(At, 6);
    m_uiMesh->unzip(At, At_uz, 6);

    m_uiMesh->readFromGhostEnd<DendroScalar>(chi, 1);
    m_uiMesh->readFromGhostEnd<DendroScalar>(K, 1);

    m_uiMesh->unzip(chi, chi_uz, 1);
    m_uiMesh->unzip(K, K_uz, 1);

    // const Point
    // pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
    // const Point
    // pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);
    // const unsigned int PW=bssn::BSSN_PADDING_WIDTH;

    // const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
    // const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        const unsigned int offset = blkList[blk].getOffset();
        const unsigned int sz[3]  = {blkList[blk].getAllocationSzX(),
                                     blkList[blk].getAllocationSzY(),
                                     blkList[blk].getAllocationSzZ()};
        const unsigned int bflag  = blkList[blk].getBlkNodeFlag();

        const DendroScalar dx     = blkList[blk].computeDx(pt_min, pt_max);
        const DendroScalar dy     = blkList[blk].computeDy(pt_min, pt_max);
        const DendroScalar dz     = blkList[blk].computeDz(pt_min, pt_max);

        DendroScalar ptmin[3], ptmax[3];

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

        const unsigned int n      = sz[0] * sz[1] * sz[2];
        const unsigned int BLK_SZ = n;

        const unsigned int nx     = sz[0];
        const unsigned int ny     = sz[1];
        const unsigned int nz     = sz[2];

        const double hx           = (ptmax[0] - ptmin[0]) / (nx - 1);
        const double hy           = (ptmax[1] - ptmin[1]) / (ny - 1);
        const double hz           = (ptmax[2] - ptmin[2]) / (nz - 1);

        double* grad_0_chi =
            &m_aeh_vars->grad_chi.get_vec_ptr()[0 * uz_sz + offset];
        double* grad_1_chi =
            &m_aeh_vars->grad_chi.get_vec_ptr()[1 * uz_sz + offset];
        double* grad_2_chi =
            &m_aeh_vars->grad_chi.get_vec_ptr()[2 * uz_sz + offset];

        double* grad_0_gt0 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[0 * 3 * uz_sz + 0 * uz_sz + offset];
        double* grad_1_gt0 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[0 * 3 * uz_sz + 1 * uz_sz + offset];
        double* grad_2_gt0 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[0 * 3 * uz_sz + 2 * uz_sz + offset];

        double* grad_0_gt1 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[1 * 3 * uz_sz + 0 * uz_sz + offset];
        ;
        double* grad_1_gt1 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[1 * 3 * uz_sz + 1 * uz_sz + offset];
        ;
        double* grad_2_gt1 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[1 * 3 * uz_sz + 2 * uz_sz + offset];
        ;

        double* grad_0_gt2 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[2 * 3 * uz_sz + 0 * uz_sz + offset];
        double* grad_1_gt2 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[2 * 3 * uz_sz + 1 * uz_sz + offset];
        double* grad_2_gt2 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[2 * 3 * uz_sz + 2 * uz_sz + offset];

        double* grad_0_gt3 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[3 * 3 * uz_sz + 0 * uz_sz + offset];
        double* grad_1_gt3 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[3 * 3 * uz_sz + 1 * uz_sz + offset];
        double* grad_2_gt3 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[3 * 3 * uz_sz + 2 * uz_sz + offset];

        double* grad_0_gt4 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[4 * 3 * uz_sz + 0 * uz_sz + offset];
        double* grad_1_gt4 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[4 * 3 * uz_sz + 1 * uz_sz + offset];
        double* grad_2_gt4 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[4 * 3 * uz_sz + 2 * uz_sz + offset];

        double* grad_0_gt5 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[5 * 3 * uz_sz + 0 * uz_sz + offset];
        double* grad_1_gt5 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[5 * 3 * uz_sz + 1 * uz_sz + offset];
        double* grad_2_gt5 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[5 * 3 * uz_sz + 2 * uz_sz + offset];

        const double* const gt0 = &gt_uz[0 * uz_sz + offset];
        const double* const gt1 = &gt_uz[1 * uz_sz + offset];
        const double* const gt2 = &gt_uz[2 * uz_sz + offset];
        const double* const gt3 = &gt_uz[3 * uz_sz + offset];
        const double* const gt4 = &gt_uz[4 * uz_sz + offset];
        const double* const gt5 = &gt_uz[5 * uz_sz + offset];

        const double* const chi = &chi_uz[0 * uz_sz + offset];

        deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
        deriv_y(grad_1_gt0, gt0, hy, sz, bflag);
        deriv_z(grad_2_gt0, gt0, hz, sz, bflag);

        deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
        deriv_y(grad_1_gt1, gt1, hy, sz, bflag);
        deriv_z(grad_2_gt1, gt1, hz, sz, bflag);

        deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
        deriv_y(grad_1_gt2, gt2, hy, sz, bflag);
        deriv_z(grad_2_gt2, gt2, hz, sz, bflag);

        deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
        deriv_y(grad_1_gt3, gt3, hy, sz, bflag);
        deriv_z(grad_2_gt3, gt3, hz, sz, bflag);

        deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
        deriv_y(grad_1_gt4, gt4, hy, sz, bflag);
        deriv_z(grad_2_gt4, gt4, hz, sz, bflag);

        deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
        deriv_y(grad_1_gt5, gt5, hy, sz, bflag);
        deriv_z(grad_2_gt5, gt5, hz, sz, bflag);

        deriv_x(grad_0_chi, chi, hx, sz, bflag);
        deriv_y(grad_1_chi, chi, hy, sz, bflag);
        deriv_z(grad_2_chi, chi, hz, sz, bflag);
    }

    if (verbose && !rank) {
        std::cout << "========================================================="
                     "============================================"
                  << std::endl;
        std::cout << "                      AEH Solver Initialized             "
                     "                                            "
                  << std::endl;
        std::cout << "========================================================="
                     "============================================"
                  << std::endl;
        std::cout << "l max = " << m_sph_modes[m_sph_modes.size() - 1].first;
        // std::cout<<" lm : ";
        // for (unsigned int lm_idx=0; lm_idx<m_sph_modes.size(); lm_idx++)
        //     std::cout<<"( "<<m_sph_modes[lm_idx].first<<",
        //     "<<m_sph_modes[lm_idx].second<<" )  ";
        std::cout << std::endl;
        std::cout << "quad theta : " << m_num_theta << std::endl;
        std::cout << "quad phi   : " << m_num_phi << std::endl;
        std::cout << "========================================================="
                     "============================================"
                  << std::endl;
    }
}

template <typename Ctx, typename T>
SpectralAEHSolver<Ctx, T>::~SpectralAEHSolver() {
    ot::Mesh* m_uiMesh = m_ctx->get_mesh();
    if (!m_uiMesh->isActive()) return;

    m_sph_modes.clear();
    gsl_integration_fixed_free(m_quad_theta);
    gsl_integration_fixed_free(m_quad_theta1);
    gsl_integration_fixed_free(m_quad_phi);

    m_aeh_vars->gt.destroy_vector();
    m_aeh_vars->At.destroy_vector();
    m_aeh_vars->chi.destroy_vector();
    m_aeh_vars->K.destroy_vector();

    m_aeh_vars->grad_gt.destroy_vector();
    m_aeh_vars->grad_chi.destroy_vector();

    delete m_aeh_vars;

    m_sin_qtheta.clear();
    m_cos_qtheta.clear();

    m_sin_qphi.clear();
    m_cos_qphi.clear();

    return;
}

template <typename Ctx, typename T>
int SpectralAEHSolver<Ctx, T>::eval_aeh_level_set(ot::Mesh* m_uiMesh,
                                                  T* const aeh_f_ptr,
                                                  const T* const h_qs) {
    if (!(m_uiMesh->isActive())) return 0;

    const Point& origin         = m_origin;
    const ot::TreeNode* pNodes  = &(*(m_uiMesh->getAllElements().begin()));
    const unsigned int eleOrder = m_uiMesh->getElementOrder();
    const unsigned int* e2n_cg  = &(*(m_uiMesh->getE2NMapping().begin()));
    const unsigned int* e2n_dg  = &(*(m_uiMesh->getE2NMapping_DG().begin()));
    const unsigned int nPe      = m_uiMesh->getNumNodesPerElement();
    const unsigned int nodeLocalBegin = m_uiMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd   = m_uiMesh->getNodeLocalEnd();

    const DendroScalar rmin           = m_rlim[0];
    const DendroScalar rmax           = m_rlim[1];

    // for (unsigned int elem = m_uiMesh->getElementLocalBegin(); elem <
    // m_uiMesh->getElementLocalEnd(); elem++)
    for (unsigned int idx = 0; idx < m_active_ele_local_id.size(); idx++) {
        const unsigned int elem = m_active_ele_local_id[idx];
        for (unsigned int k = 0; k < (eleOrder + 1); k++)
            for (unsigned int j = 0; j < (eleOrder + 1); j++)
                for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                    const unsigned int idx =
                        elem * nPe + k * (eleOrder + 1) * (eleOrder + 1) +
                        j * (eleOrder + 1) + i;
                    const unsigned int nodeLookUp_CG = e2n_cg[idx];

                    if (nodeLookUp_CG >= nodeLocalBegin &&
                        nodeLookUp_CG < nodeLocalEnd) {
                        const unsigned int nodeLookUp_DG = e2n_dg[idx];
                        unsigned int ownerID, ii_x, jj_y, kk_z;

                        m_uiMesh->dg2eijk(nodeLookUp_DG, ownerID, ii_x, jj_y,
                                          kk_z);
                        const DendroScalar len =
                            (double)(1u << (m_uiMaxDepth -
                                            pNodes[ownerID].getLevel()));

                        const DendroScalar x =
                            pNodes[ownerID].getX() + ii_x * (len / (eleOrder));
                        const DendroScalar y =
                            pNodes[ownerID].getY() + jj_y * (len / (eleOrder));
                        const DendroScalar z =
                            pNodes[ownerID].getZ() + kk_z * (len / (eleOrder));

                        const DendroScalar xx = GRIDX_TO_X(x) - origin.x();
                        const DendroScalar yy = GRIDY_TO_Y(y) - origin.y();
                        const DendroScalar zz = GRIDZ_TO_Z(z) - origin.z();

                        const DendroScalar p_rr =
                            sqrt(xx * xx + yy * yy + zz * zz);

                        DendroScalar h_tp = (T)0;
                        if (p_rr > rmin && p_rr < rmax) {
                            const DendroScalar p_pt = acos(zz / p_rr);
                            const DendroScalar p_ph =
                                std::fmod(atan2(yy, xx), 2 * M_PI);
                            for (unsigned int lm_idx = 0;
                                 lm_idx < m_sph_modes.size(); lm_idx++) {
                                int l = m_sph_modes[lm_idx].first;
                                if (l > m_ll) continue;

                                h_tp +=
                                    h_qs[lm_idx] *
                                    real_spherical_harmonic(
                                        m_sph_modes[lm_idx].first,
                                        m_sph_modes[lm_idx].second, p_pt, p_ph);
                            }
                        }

                        aeh_f_ptr[nodeLookUp_CG] = p_rr - h_tp;
                    }
                }
    }

    return 0;
}

template <typename Ctx, typename T>
int SpectralAEHSolver<Ctx, T>::eval_quasi_measurements_integrands(
    Ctx* ctx, ot::Mesh* m_uiMesh, T* const u_ptr, const T* const h_qs) {
    if (!(m_uiMesh->isActive())) return 0;

    const Point& origin         = m_origin;
    const ot::TreeNode* pNodes  = &(*(m_uiMesh->getAllElements().begin()));
    const unsigned int eleOrder = m_uiMesh->getElementOrder();
    const unsigned int* e2n_cg  = &(*(m_uiMesh->getE2NMapping().begin()));
    const unsigned int* e2n_dg  = &(*(m_uiMesh->getE2NMapping_DG().begin()));
    const unsigned int nPe      = m_uiMesh->getNumNodesPerElement();
    const unsigned int nodeLocalBegin = m_uiMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd   = m_uiMesh->getNodeLocalEnd();

    const DendroScalar rmin           = m_rlim[0];
    const DendroScalar rmax           = m_rlim[1];

    DVec& gr_vars                     = ctx->get_evolution_vars();
    DendroScalar* gr_vars_ptr         = gr_vars.get_vec_ptr();
    const unsigned int cg_sz          = m_uiMesh->getDegOfFreedom();

    const DendroScalar* const gt0 = &gr_vars_ptr[bssn::VAR::U_SYMGT0 * cg_sz];
    const DendroScalar* const gt1 = &gr_vars_ptr[bssn::VAR::U_SYMGT1 * cg_sz];
    const DendroScalar* const gt2 = &gr_vars_ptr[bssn::VAR::U_SYMGT2 * cg_sz];
    const DendroScalar* const gt3 = &gr_vars_ptr[bssn::VAR::U_SYMGT3 * cg_sz];
    const DendroScalar* const gt4 = &gr_vars_ptr[bssn::VAR::U_SYMGT4 * cg_sz];
    const DendroScalar* const gt5 = &gr_vars_ptr[bssn::VAR::U_SYMGT5 * cg_sz];

    const DendroScalar* const At0 = &gr_vars_ptr[bssn::VAR::U_SYMAT0 * cg_sz];
    const DendroScalar* const At1 = &gr_vars_ptr[bssn::VAR::U_SYMAT1 * cg_sz];
    const DendroScalar* const At2 = &gr_vars_ptr[bssn::VAR::U_SYMAT2 * cg_sz];
    const DendroScalar* const At3 = &gr_vars_ptr[bssn::VAR::U_SYMAT3 * cg_sz];
    const DendroScalar* const At4 = &gr_vars_ptr[bssn::VAR::U_SYMAT4 * cg_sz];
    const DendroScalar* const At5 = &gr_vars_ptr[bssn::VAR::U_SYMAT5 * cg_sz];

    const DendroScalar* const K   = &gr_vars_ptr[bssn::VAR::U_K * cg_sz];
    const DendroScalar* const chi = &gr_vars_ptr[bssn::VAR::U_CHI * cg_sz];

    T* const sqrt_det_m_ab        = &u_ptr[0 * cg_sz];
    T* const J_x                  = &u_ptr[1 * cg_sz];
    T* const J_y                  = &u_ptr[2 * cg_sz];
    T* const J_z                  = &u_ptr[3 * cg_sz];

    // for (unsigned int elem = m_uiMesh->getElementLocalBegin(); elem <
    // m_uiMesh->getElementLocalEnd(); elem++) {
    for (unsigned int idx = 0; idx < m_active_ele_local_id.size(); idx++) {
        const unsigned int elem = m_active_ele_local_id[idx];
        for (unsigned int k = 0; k < (eleOrder + 1); k++)
            for (unsigned int j = 0; j < (eleOrder + 1); j++)
                for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                    const unsigned int idx =
                        elem * nPe + k * (eleOrder + 1) * (eleOrder + 1) +
                        j * (eleOrder + 1) + i;
                    const unsigned int nodeLookUp_CG = e2n_cg[idx];

                    if (nodeLookUp_CG >= nodeLocalBegin &&
                        nodeLookUp_CG < nodeLocalEnd) {
                        const unsigned int nodeLookUp_DG = e2n_dg[idx];
                        unsigned int ownerID, ii_x, jj_y, kk_z;

                        m_uiMesh->dg2eijk(nodeLookUp_DG, ownerID, ii_x, jj_y,
                                          kk_z);
                        const DendroScalar len =
                            (double)(1u << (m_uiMaxDepth -
                                            pNodes[ownerID].getLevel()));

                        const DendroScalar x =
                            pNodes[ownerID].getX() + ii_x * (len / (eleOrder));
                        const DendroScalar y =
                            pNodes[ownerID].getY() + jj_y * (len / (eleOrder));
                        const DendroScalar z =
                            pNodes[ownerID].getZ() + kk_z * (len / (eleOrder));

                        const DendroScalar xx = GRIDX_TO_X(x) - origin.x();
                        const DendroScalar yy = GRIDY_TO_Y(y) - origin.y();
                        const DendroScalar zz = GRIDZ_TO_Z(z) - origin.z();

                        const DendroScalar p_rr =
                            sqrt(xx * xx + yy * yy + zz * zz);

                        if (p_rr > rmin && p_rr < rmax) {
                            const DendroScalar p_pt = acos(zz / p_rr);
                            const DendroScalar p_ph =
                                std::fmod(atan2(yy, xx), 2 * M_PI);

                            DendroScalar h_tp      = (T)0;
                            DendroScalar d_theta_h = (T)0;
                            DendroScalar d_phi_h   = (T)0;

                            // computing the spectral derivatives for h(theta,
                            // phi)
                            for (unsigned int lm_idx = 0;
                                 lm_idx < m_sph_modes.size(); lm_idx++) {
                                int l = m_sph_modes[lm_idx].first;
                                if (l > m_ll) continue;

                                // Need to know the gradient of f wrt the 3
                                // cartesian components df/dx, df/dy, df/dz. df=
                                // r-h(theta, phi) Nees to get r = sqrt(x^2 +
                                // y^2 + z^2). How do we differentiate h? dh/dxi
                                // = dh/dtheta dtheta/dxi + dh/dphi dphi/dxi h
                                // as a function of theta and phi and the
                                // folllowing 2 lines are the derivatives

                                h_tp +=
                                    h_qs[lm_idx] *
                                    real_spherical_harmonic(
                                        m_sph_modes[lm_idx].first,
                                        m_sph_modes[lm_idx].second, p_pt, p_ph);
                                d_theta_h += h_qs[lm_idx] *
                                             real_spherical_harmonic_d_theta(
                                                 m_sph_modes[lm_idx].first,
                                                 m_sph_modes[lm_idx].second,
                                                 p_pt, p_ph, 1);
                                d_phi_h += h_qs[lm_idx] *
                                           real_spherical_harmonic_d_phi(
                                               m_sph_modes[lm_idx].first,
                                               m_sph_modes[lm_idx].second, p_pt,
                                               p_ph, 1);
                            }

                            // Compute r
                            DendroScalar r = sqrt(xx * xx + yy * yy + zz * zz);

                            // Partial derivatives of theta
                            DendroScalar dtheta_dx =
                                xx * zz /
                                (r * r * r * sqrt(1 - zz * zz / (r * r)));
                            DendroScalar dtheta_dy =
                                yy * zz /
                                (r * r * r * sqrt(1 - zz * zz / (r * r)));
                            DendroScalar dtheta_dz =
                                (1 / (r * sqrt(1 - zz * zz / (r * r)))) -
                                (zz * zz /
                                 (r * r * r * sqrt(1 - zz * zz / (r * r))));

                            // Partial derivatives of phi
                            DendroScalar dphi_dx = -yy / (xx * xx + yy * yy);
                            DendroScalar dphi_dy = xx / (xx * xx + yy * yy);
                            DendroScalar dphi_dz = 0;

                            // Gradients of h
                            DendroScalar dh_dx =
                                d_theta_h * dtheta_dx + d_phi_h * dphi_dx;
                            DendroScalar dh_dy =
                                d_theta_h * dtheta_dy + d_phi_h * dphi_dy;
                            DendroScalar dh_dz =
                                d_theta_h * dtheta_dz + d_phi_h * dphi_dz;

                            // Gradients of f
                            DendroScalar grad_0_F = x / r - dh_dx;
                            DendroScalar grad_1_F = y / r - dh_dy;
                            DendroScalar grad_2_F = z / r - dh_dz;

                            // embedding coordinate maps
                            // const DendroScalar A0 = h_tp * sin(p_pt) *
                            // cos(p_ph); const DendroScalar A1 = h_tp *
                            // sin(p_pt) * sin(p_ph); const DendroScalar A2 =
                            // h_tp * cos(p_pt);

                            const DendroScalar grad_0_A0 =
                                (h_tp * cos(p_pt) + sin(p_pt) * d_theta_h) *
                                cos(p_ph);
                            const DendroScalar grad_0_A1 =
                                (h_tp * cos(p_pt) + sin(p_pt) * d_theta_h) *
                                sin(p_ph);
                            const DendroScalar grad_0_A2 =
                                (-h_tp * sin(p_pt) + cos(p_pt) * d_theta_h);

                            const DendroScalar grad_1_A0 =
                                (-h_tp * sin(p_ph) + cos(p_ph) * d_phi_h) *
                                sin(p_pt);
                            const DendroScalar grad_1_A1 =
                                (h_tp * cos(p_ph) + sin(p_ph) * d_phi_h) *
                                sin(p_pt);
                            const DendroScalar grad_1_A2 = d_phi_h * cos(p_pt);

                            const unsigned int pp        = nodeLookUp_CG;
#include "../src/det_metric_aeh.cpp"

                        } else {
                            sqrt_det_m_ab[nodeLookUp_CG] = (T)0.0;
                            J_x[nodeLookUp_CG]           = (T)0.0;
                            J_y[nodeLookUp_CG]           = (T)0.0;
                            J_z[nodeLookUp_CG]           = (T)0.0;
                        }
                    }
                }
    }

    return 0;
}

template <typename Ctx, typename T>
DendroScalar SpectralAEHSolver<Ctx, T>::eval_expansion_norm(
    Ctx* ctx, const T* const h_qs, DVec& aeh_f, DVec& aeh_h,
    std::vector<T>& interp_coords, unsigned int norm_type) {
    ot::Mesh* m_uiMesh = ctx->get_mesh();
    if (!(m_uiMesh->isActive())) return 0.0;

    int rank            = m_uiMesh->getMPIRank();
    // int npes = m_uiMesh->getMPICommSize();

    const Point& origin = m_origin;
    T* aeh_f_ptr        = aeh_f.get_vec_ptr();
    T* aeh_h_ptr        = aeh_h.get_vec_ptr();

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

    this->eval_aeh_level_set(m_uiMesh, aeh_f_ptr, h_qs);
    ctx->aeh_expansion(origin, m_aeh_vars, aeh_f, aeh_h, m_rlim);

    m_uiMesh->readFromGhostBegin(aeh_h.get_vec_ptr(), 1);

    unsigned int num_angular_pts = m_num_theta * m_num_phi;
    const double* const m_qtheta = m_quad_theta->x;
    const double* const m_qphi   = m_quad_phi->x;

    for (unsigned int qt = 0; qt < m_num_theta; qt++)
        for (unsigned int qp = 0; qp < m_num_phi; qp++) {
            const unsigned int q_idx = qt * m_num_phi + qp;

            T r_val                  = (T)0;
            for (unsigned int lm_idx = 0; lm_idx < m_sph_modes.size();
                 lm_idx++) {
                int l = m_sph_modes[lm_idx].first;
                if (l > m_ll) continue;

                r_val += real_spherical_harmonic(m_sph_modes[lm_idx].first,
                                                 m_sph_modes[lm_idx].second,
                                                 m_qtheta[qt], m_qphi[qp]) *
                         h_qs[lm_idx];
            }

            interp_coords[3 * q_idx + 0] =
                origin.x() + (r_val)*m_sin_qtheta[qt] * m_cos_qphi[qp];
            interp_coords[3 * q_idx + 1] =
                origin.y() + (r_val)*m_sin_qtheta[qt] * m_sin_qphi[qp];
            interp_coords[3 * q_idx + 2] =
                origin.z() + (r_val)*m_cos_qtheta[qt];
        }

    m_uiMesh->readFromGhostEnd(aeh_h.get_vec_ptr(), 1);

    std::vector<unsigned int> valid_idx;
    valid_idx.clear();

    std::vector<T> aeh_h_inp;
    aeh_h_inp.resize(num_angular_pts);
    ot::da::interpolateToCoords(m_uiMesh, aeh_h_ptr, interp_coords.data(),
                                interp_coords.size(), grid_limits,
                                domain_limits, aeh_h_inp.data(), valid_idx);
    DendroScalar result       = 0;
    DendroScalar result_g     = 0;

    DendroScalar min_max[2]   = {1e30, 0};
    DendroScalar min_max_g[2] = {0, 0};

    if (norm_type == 0) {
        for (unsigned int idx = 0; idx < valid_idx.size(); idx++) {
            DendroScalar s = abs(aeh_h_inp[valid_idx[idx]]);
            if (s > result) result = s;

            if (s > min_max[1]) min_max[1] = s;

            if (min_max[0] > s) min_max[0] = s;
        }

        par::Mpi_Allreduce(&result, &result_g, 1, MPI_MAX,
                           m_uiMesh->getMPICommunicator());
        // par::Mpi_Allreduce(&min_max[0], &min_max_g[0], 1, MPI_MIN,
        // m_uiMesh->getMPICommunicator()); par::Mpi_Allreduce(&min_max[1],
        // &min_max_g[1], 1, MPI_MAX, m_uiMesh->getMPICommunicator()); if(!rank)
        //     printf("expansion abs min =%.8E max = %.8E\n",min_max_g[0],
        //     min_max_g[1]);

    } else if (norm_type == 1) {
        for (unsigned int idx = 0; idx < valid_idx.size(); idx++) {
            const unsigned int qp = valid_idx[idx] % m_num_phi;
            const unsigned int qt = (valid_idx[idx] - qp) / m_num_phi;

            DendroScalar s        = abs(aeh_h_inp[valid_idx[idx]]);
            result += (s)*m_quad_theta->weights[qt] * m_quad_phi->weights[qp];
        }
        par::Mpi_Allreduce(&result, &result_g, 1, MPI_SUM,
                           m_uiMesh->getMPICommunicator());

    } else if (norm_type == 2) {
        for (unsigned int idx = 0; idx < valid_idx.size(); idx++) {
            const unsigned int qp = valid_idx[idx] % m_num_phi;
            const unsigned int qt = (valid_idx[idx] - qp) / m_num_phi;

            DendroScalar s        = aeh_h_inp[valid_idx[idx]];
            result +=
                (s * s) * m_quad_theta->weights[qt] * m_quad_phi->weights[qp];
        }
        par::Mpi_Allreduce(&result, &result_g, 1, MPI_SUM,
                           m_uiMesh->getMPICommunicator());
    }
    return result_g;
}

template <typename Ctx, typename T>
AEHErrorType SpectralAEHSolver<Ctx, T>::solve(Ctx* ctx, const T* const h_init,
                                              T* h_qs, unsigned int max_iter,
                                              double rel_eps, double abs_eps,
                                              double alpha, double beta,
                                              unsigned int verbose) {
    ot::Mesh* m_uiMesh = ctx->get_mesh();
    if (!m_uiMesh->isActive()) return AEHErrorType::SUCCESS;

    const Point& origin       = m_origin;
    const int rank            = m_uiMesh->getMPIRank();
    const unsigned int num_lm = m_sph_modes.size();
    const unsigned int cg_sz  = m_uiMesh->getDegOfFreedom();

    DendroScalar* h_qs0       = new DendroScalar[num_lm];
    DendroScalar* h_qs1       = new DendroScalar[num_lm];

    std::memcpy(h_qs0, h_init, sizeof(T) * num_lm);
    std::memcpy(h_qs1, h_init, sizeof(T) * num_lm);

    // DVec aeh_v;
    // aeh_v.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
    // ot::DVEC_LOC::HOST, 2, true);

    DVec aeh_f;
    DVec aeh_h;

    // DendroScalar * aeh_v_ptr_0 = aeh_v.get_vec_ptr();
    // DendroScalar * aeh_v_ptr_1 = aeh_v.get_vec_ptr() + cg_sz;

    // aeh_f.set_vec_ptr(aeh_v_ptr_0 , m_uiMesh,
    // ot::DVEC_TYPE::OCT_SHARED_NODES, ot::DVEC_LOC::HOST, 1, true);
    // aeh_h.set_vec_ptr(aeh_v_ptr_1 , m_uiMesh,
    // ot::DVEC_TYPE::OCT_SHARED_NODES, ot::DVEC_LOC::HOST, 1, true);

    aeh_f.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                        ot::DVEC_LOC::HOST, 1, true);
    aeh_h.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                        ot::DVEC_LOC::HOST, 1, true);

    T* aeh_h_ptr                 = aeh_h.get_vec_ptr();
    T* aeh_f_ptr                 = aeh_f.get_vec_ptr();

    unsigned int num_angular_pts = m_num_theta * m_num_phi;

    std::vector<T> interp_coords;
    interp_coords.resize(3 * num_angular_pts);

    std::vector<T> aeh_r;
    aeh_r.resize(num_angular_pts);

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

    unsigned int iter            = 0;
    const double* const m_qtheta = m_quad_theta->x;
    const double* const m_qphi   = m_quad_phi->x;

    std::vector<T> aeh_h_inp;
    aeh_h_inp.resize(num_angular_pts);
    double relative_error = 0;
    double absolute_error = 0;

    this->eval_aeh_level_set(m_uiMesh, aeh_f_ptr, h_qs0);
    ctx->aeh_expansion(origin, m_aeh_vars, aeh_f, aeh_h, m_rlim);
    const double abs_r0 =
        eval_expansion_norm(ctx, h_qs0, aeh_f, aeh_h, interp_coords, 2);

    do {
        // printArray_1D(interp_coords.data(), interp_coords.size());
        this->eval_aeh_level_set(m_uiMesh, aeh_f_ptr, h_qs0);
        ctx->aeh_expansion(origin, m_aeh_vars, aeh_f, aeh_h, m_rlim);

        m_uiMesh->readFromGhostBegin(aeh_h.get_vec_ptr(), 1);

        for (unsigned int qt = 0; qt < m_num_theta; qt++)
            for (unsigned int qp = 0; qp < m_num_phi; qp++) {
                const unsigned int q_idx = qt * m_num_phi + qp;
                T r_val                  = (T)0;
                for (unsigned int lm_idx = 0; lm_idx < m_sph_modes.size();
                     lm_idx++) {
                    int l = m_sph_modes[lm_idx].first;
                    if (l > m_ll) continue;

                    r_val += real_spherical_harmonic(m_sph_modes[lm_idx].first,
                                                     m_sph_modes[lm_idx].second,
                                                     m_qtheta[qt], m_qphi[qp]) *
                             h_qs0[lm_idx];
                }

                aeh_r[q_idx] = r_val;

                interp_coords[3 * q_idx + 0] =
                    origin.x() + (r_val)*m_sin_qtheta[qt] * m_cos_qphi[qp];
                interp_coords[3 * q_idx + 1] =
                    origin.y() + (r_val)*m_sin_qtheta[qt] * m_sin_qphi[qp];
                interp_coords[3 * q_idx + 2] =
                    origin.z() + (r_val)*m_cos_qtheta[qt];
            }

        m_uiMesh->readFromGhostEnd(aeh_h.get_vec_ptr(), 1);

        std::vector<unsigned int> valid_idx;
        valid_idx.clear();

        ot::da::interpolateToCoords(m_uiMesh, aeh_h_ptr, interp_coords.data(),
                                    interp_coords.size(), grid_limits,
                                    domain_limits, aeh_h_inp.data(), valid_idx);
        if (m_solver_type == AEHSolverType::NKO) {
            const int l_max =
                std::max(10, 4 * m_sph_modes[m_sph_modes.size() - 1].first);
            const DendroScalar dlambda =
                (1.0 / (DendroScalar)(l_max * (l_max + 1)));
            for (unsigned int lm_idx = 1; lm_idx < m_sph_modes.size();
                 lm_idx++) {
                DendroScalar tmp  = 0.0;
                DendroScalar a_qs = 0.0;

                int l             = m_sph_modes[lm_idx].first;
                int m             = m_sph_modes[lm_idx].second;

                if (l > m_ll) continue;

                for (unsigned int idx = 0; idx < valid_idx.size(); idx++) {
                    const unsigned int qp = valid_idx[idx] % m_num_phi;
                    const unsigned int qt = (valid_idx[idx] - qp) / m_num_phi;
                    tmp += aeh_h_inp[valid_idx[idx]] *
                           real_spherical_harmonic(l, m, m_qtheta[qt],
                                                   m_qphi[qp]) *
                           m_quad_theta->weights[qt] * m_quad_phi->weights[qp];
                }

                par::Mpi_Allreduce(&tmp, &a_qs, 1, MPI_SUM,
                                   m_uiMesh->getMPICommunicator());
                // if(!rank)
                //     printf("lm=(%d, %d) integral =%.8E\n",l,m,a_qs);
                h_qs1[lm_idx] = h_qs0[lm_idx] - dlambda * a_qs;
            }

            h_qs1[0] = h_qs0[0];
            h_qs1[0] = this->solve_00_newton(ctx, h_qs1, aeh_f, aeh_h,
                                             interp_coords, 1e-6);
            if (h_qs1[0] == -1.0) return AEHErrorType::FAIL;

        } else if (m_solver_type == AEHSolverType::FAST_FLOW) {
            const int l_max =
                2 * m_sph_modes[m_sph_modes.size() - 1]
                        .first;  // std::max(10, 4 *
                                 // m_sph_modes[m_sph_modes.size()-1].first);
            const DendroScalar A =
                (alpha / (DendroScalar)(l_max * (l_max + 1))) + beta;
            const DendroScalar B        = beta / alpha;
            const DendroScalar dlambda0 = A / (1 + B * l_max * (l_max + 1));
            DendroScalar dlambda        = dlambda0;

            for (unsigned int lm_idx = 0; lm_idx < m_sph_modes.size();
                 lm_idx++) {
                DendroScalar tmp  = 0.0;
                DendroScalar a_qs = 0.0;

                const int l       = m_sph_modes[lm_idx].first;
                const int m       = m_sph_modes[lm_idx].second;

                if (l > m_ll) continue;

                for (unsigned int idx = 0; idx < valid_idx.size(); idx++) {
                    const unsigned int qp = valid_idx[idx] % m_num_phi;
                    const unsigned int qt = (valid_idx[idx] - qp) / m_num_phi;
                    tmp += aeh_h_inp[valid_idx[idx]] *
                           real_spherical_harmonic(l, m, m_qtheta[qt],
                                                   m_qphi[qp]) *
                           m_quad_theta->weights[qt] * m_quad_phi->weights[qp];
                }

                par::Mpi_Allreduce(&tmp, &a_qs, 1, MPI_SUM,
                                   m_uiMesh->getMPICommunicator());
                // if(!rank)
                //     printf("lm=(%d, %d) integral =%.8E\n",l,m,a_qs);
                // if (lm_idx==0)
                //     dlambda = std::min(0.1 * h_qs0[0]/a_qs, dlambda0);

                h_qs1[lm_idx] = h_qs0[lm_idx] - dlambda * a_qs;
            }

        } else
            return AEHErrorType::FAIL;

        absolute_error =
            eval_expansion_norm(ctx, h_qs1, aeh_f, aeh_h, interp_coords, 2);
        relative_error = absolute_error / abs_r0;
        const DendroScalar converged =
            normL2(h_qs1, h_qs0, m_sph_modes.size()) /
            normL2(h_qs0, m_sph_modes.size());

        if (!rank && verbose == 1) {
            if ((iter + 1) % 10 == 0)
                printf(
                    "l = %02d / lmax = %02d AEH iter = %04d / %04d relative "
                    "error = %.8E  absolute error = %.8E ||hqs_1 - "
                    "hqs_0||/||hqs_0|| = %.8E\n",
                    m_ll, m_lmax, iter + 1, max_iter, relative_error,
                    absolute_error, converged);
            // for(unsigned int lm_idx=0; lm_idx < m_sph_modes.size(); lm_idx++)
            //     printf("h_{%d, %d} = %.8E ,",m_sph_modes[lm_idx].first,
            //     m_sph_modes[lm_idx].second, h_qs1[lm_idx]);
            // printf("\n");
        }

        std::swap(h_qs0, h_qs1);
        iter += 1;

        // {
        //     DendroScalar* pData[16];

        //     this->eval_aeh_level_set(m_uiMesh, aeh_f_ptr, h_qs0);
        //     ctx->aeh_expansion(origin, m_aeh_vars, aeh_f, aeh_h, m_rlim);

        //     m_uiMesh->readFromGhostBegin(aeh_h.get_vec_ptr(), 1);
        //     m_uiMesh->readFromGhostEnd(aeh_h.get_vec_ptr(), 1);

        //     pData[0] = aeh_f.get_vec_ptr();
        //     pData[1] = aeh_h.get_vec_ptr();

        //     DendroScalar* gr_evar  = ctx->get_evolution_vars().get_vec_ptr();
        //     pData[2] = &gr_evar[bssn::VAR::U_SYMGT0 * cg_sz];
        //     pData[3] = &gr_evar[bssn::VAR::U_SYMGT1 * cg_sz];
        //     pData[4] = &gr_evar[bssn::VAR::U_SYMGT2 * cg_sz];
        //     pData[5] = &gr_evar[bssn::VAR::U_SYMGT3 * cg_sz];
        //     pData[6] = &gr_evar[bssn::VAR::U_SYMGT4 * cg_sz];
        //     pData[7] = &gr_evar[bssn::VAR::U_SYMGT5 * cg_sz];

        //     pData[8]  = &gr_evar[bssn::VAR::U_SYMAT0 * cg_sz];
        //     pData[9]  = &gr_evar[bssn::VAR::U_SYMAT1 * cg_sz];
        //     pData[10] = &gr_evar[bssn::VAR::U_SYMAT2 * cg_sz];
        //     pData[11] = &gr_evar[bssn::VAR::U_SYMAT3 * cg_sz];
        //     pData[12] = &gr_evar[bssn::VAR::U_SYMAT4 * cg_sz];
        //     pData[13] = &gr_evar[bssn::VAR::U_SYMAT5 * cg_sz];

        //     pData[14] = &gr_evar[bssn::VAR::U_CHI * cg_sz];
        //     pData[15] = &gr_evar[bssn::VAR::U_K * cg_sz];

        //     const char * pNames[]={"F", "H", "gt0", "gt1", "gt2", "gt3",
        //     "gt4", "gt5", "At0", "At1", "At2", "At3", "At4", "At5", "At5",
        //     "chi", "K" }; char fname[256]; sprintf(fname,
        //     "%s_iter_%03d","vtu/aeh", iter); io::vtk::mesh2vtuFine(m_uiMesh,
        //     fname, 0, nullptr, nullptr, 16, pNames, (const double**)pData, 0,
        //     nullptr, nullptr, false);

        // }

        if (converged < rel_eps) break;

    } while ((iter < max_iter) &&
             (relative_error > rel_eps && absolute_error > abs_eps));

    std::memcpy(h_qs, h_qs0, sizeof(T) * num_lm);

    delete[] h_qs1;
    delete[] h_qs0;

    aeh_f.destroy_vector();
    aeh_h.destroy_vector();

    if (iter == max_iter &&
        (relative_error > rel_eps && absolute_error > abs_eps)) {
        if (m_ll == m_lmax) std::memcpy(h_qs, h_init, sizeof(T) * num_lm);

        if (!rank && verbose == 1) {
            printf(
                "AEH solver did not convege for abs tol = %.8E / %.8E rel tol "
                "= %.8E / %.8E with max iterations = %04d\n",
                absolute_error, abs_eps, relative_error, rel_eps, iter);
        }
        return AEHErrorType::FAIL;
    } else {
        if (!rank && verbose == 1)
            printf(
                "ll = %02d AEH solver converged with iter = %04d abs = %.8E "
                "rel = %.8E\n",
                m_ll, iter, absolute_error, relative_error);

        return AEHErrorType::SUCCESS;
    }
}

template <typename Ctx, typename T>
void SpectralAEHSolver<Ctx, T>::aeh_quasi_measurements(Ctx* ctx,
                                                       const T* const h_qs,
                                                       AEH_QoI& qoi) {
    ot::Mesh* m_uiMesh = ctx->get_mesh();
    if (!(m_uiMesh->isActive())) return;

    const Point& origin      = m_origin;
    int rank                 = m_uiMesh->getMPIRank();
    int npes                 = m_uiMesh->getMPICommSize();
    const unsigned int cg_sz = m_uiMesh->getDegOfFreedom();

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

    unsigned int num_angular_pts = m_num_theta * m_num_phi;
    const double* const m_qtheta = m_quad_theta1->x;
    const double* const m_qphi   = m_quad_phi->x;

    std::vector<T> interp_coords;
    interp_coords.resize(3 * num_angular_pts);

    DVec aeh_f;
    DVec aeh_h;

    DVec aeh_qi;

    aeh_f.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                        ot::DVEC_LOC::HOST, 1, true);
    aeh_h.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                        ot::DVEC_LOC::HOST, 1, true);
    aeh_qi.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                         ot::DVEC_LOC::HOST, NUM_QUASI_MEASURE_VARS, true);

    qoi.expansion_L2 =
        this->eval_expansion_norm(ctx, h_qs, aeh_f, aeh_h, interp_coords, 2);

    T* aeh_qi_ptr = aeh_qi.get_vec_ptr();
    this->eval_quasi_measurements_integrands(ctx, m_uiMesh, aeh_qi_ptr, h_qs);
    m_uiMesh->readFromGhostBegin(aeh_qi_ptr, NUM_QUASI_MEASURE_VARS);

    DendroScalar r_min_mean_max[3] = {1e6, 0, 0};
    DendroScalar f1                = 1.0 / (DendroScalar)num_angular_pts;

    for (unsigned int qt = 0; qt < m_num_theta; qt++)
        for (unsigned int qp = 0; qp < m_num_phi; qp++) {
            const unsigned int q_idx = qt * m_num_phi + qp;

            T r_val                  = (T)0;
            for (unsigned int lm_idx = 0; lm_idx < m_sph_modes.size();
                 lm_idx++) {
                int l = m_sph_modes[lm_idx].first;
                if (l > m_ll) continue;
                r_val += real_spherical_harmonic(m_sph_modes[lm_idx].first,
                                                 m_sph_modes[lm_idx].second,
                                                 m_qtheta[qt], m_qphi[qp]) *
                         h_qs[lm_idx];
            }

            if (r_min_mean_max[0] > r_val) r_min_mean_max[0] = r_val;

            if (r_min_mean_max[2] < r_val) r_min_mean_max[2] = r_val;

            r_min_mean_max[1] += r_val * f1;

            interp_coords[3 * q_idx + 0] =
                origin.x() +
                (r_val)*sin(m_qtheta[qt]) *
                    cos(m_qphi[qp]);  //* m_sin_qtheta[qt] * m_cos_qphi[qp];
            interp_coords[3 * q_idx + 1] =
                origin.y() +
                (r_val)*sin(m_qtheta[qt]) *
                    sin(m_qphi[qp]);  //* m_sin_qtheta[qt] * m_sin_qphi[qp];
            interp_coords[3 * q_idx + 2] =
                origin.z() + (r_val)*cos(m_qtheta[qt]);  //* m_cos_qtheta[qt];
        }

    m_uiMesh->readFromGhostEnd(aeh_qi_ptr, NUM_QUASI_MEASURE_VARS);

    std::vector<T> ws_v1_inp;
    ws_v1_inp.resize(num_angular_pts);
    DendroScalar result_g[NUM_QUASI_MEASURE_VARS];

    for (unsigned int var_idx = 0; var_idx < NUM_QUASI_MEASURE_VARS;
         var_idx++) {
        std::vector<unsigned int> valid_idx;
        valid_idx.clear();
        ot::da::interpolateToCoords(m_uiMesh, aeh_qi_ptr + var_idx * cg_sz,
                                    interp_coords.data(), interp_coords.size(),
                                    grid_limits, domain_limits,
                                    ws_v1_inp.data(), valid_idx);

        DendroScalar result = 0;

        for (unsigned int idx = 0; idx < valid_idx.size(); idx++) {
            const unsigned int qp = valid_idx[idx] % m_num_phi;
            const unsigned int qt = (valid_idx[idx] - qp) / m_num_phi;

            result += ws_v1_inp[valid_idx[idx]] * m_quad_theta1->weights[qt] *
                      m_quad_phi->weights[qp];
        }

        par::Mpi_Allreduce(&result, &result_g[var_idx], 1, MPI_SUM,
                           m_uiMesh->getMPICommunicator());
    }

    aeh_f.destroy_vector();
    aeh_h.destroy_vector();
    aeh_qi.destroy_vector();

    qoi.area  = result_g[0];

    qoi.Jx    = result_g[1];
    qoi.Jy    = result_g[2];
    qoi.Jz    = result_g[3];

    qoi.rmin  = r_min_mean_max[0];
    qoi.rmean = r_min_mean_max[1];
    qoi.rmax  = r_min_mean_max[2];
    return;
}

template <typename Ctx, typename T>
int SpectralAEHSolver<Ctx, T>::aeh_to_json(Ctx* ctx, T* h_qs, char* fname,
                                           std::ios_base::openmode mode) {
    ot::Mesh* m_uiMesh = ctx->get_mesh();
    if (!m_uiMesh->isActive()) return 0;

    const Point& origin = m_origin;
    int rank            = m_uiMesh->getMPIRank();
    int npes            = m_uiMesh->getMPICommSize();

    AEH_QoI ah_qoi;
    aeh_quasi_measurements(ctx, h_qs, ah_qoi);

    if (!rank) {
        try {
            std::ofstream file_aeh;
            file_aeh.open(fname, std::ofstream::app);
            file_aeh.precision(12);
            file_aeh << std::scientific;
            file_aeh << "{";
            file_aeh << "\"time\":" << bssn::BSSN_CURRENT_RK_COORD_TIME << ", ";
            file_aeh << "\"origin\":[" << origin.x() << "," << origin.y() << ","
                     << origin.z() << "], ";

            file_aeh << "\"lm\":[[0,0]";
            for (unsigned int lm_idx = 1; lm_idx < m_sph_modes.size(); lm_idx++)
                file_aeh << ",[" << m_sph_modes[lm_idx].first << ","
                         << m_sph_modes[lm_idx].second << "]";
            file_aeh << "], ";

            file_aeh << "\"coeff\":[" << h_qs[0];
            for (unsigned int lm_idx = 1; lm_idx < m_sph_modes.size(); lm_idx++)
                file_aeh << "," << h_qs[lm_idx];
            file_aeh << "], ";

            file_aeh << "\"expansion_L2\":" << ah_qoi.expansion_L2 << ", ";
            file_aeh << "\"area\":" << ah_qoi.area << ", ";
            file_aeh << "\"Jx\":" << ah_qoi.Jx << ", ";
            file_aeh << "\"Jy\":" << ah_qoi.Jy << ", ";
            file_aeh << "\"Jz\":" << ah_qoi.Jz << ", ";
            file_aeh << "\"rmin\":" << ah_qoi.rmin << ", ";
            file_aeh << "\"rmean\":" << ah_qoi.rmean << ", ";
            file_aeh << "\"rmax\":" << ah_qoi.rmax << "";

            // file_aeh<<"\"quad\":["<<m_qtheta[0]<<", "<<m_qphi[0];
            // for(unsigned int qt=0; qt < m_num_theta; qt++)
            //     for(unsigned int qp=0; qp < m_num_phi; qp++)
            //     {
            //         if(qt==0 && qp==0)
            //             continue;

            //         file_aeh<<","<<m_qtheta[qt]<<", "<<m_qphi[qp];

            //     }
            // file_aeh<<"],";

            // file_aeh<<"\"expansion\":["<<aeh_h_inp[0];
            // for(unsigned int qt=0; qt < m_num_theta; qt++)
            //     for(unsigned int qp=0; qp < m_num_phi; qp++)
            //     {
            //         if(qt==0 && qp==0)
            //             continue;

            //         const unsigned int q_idx = qt * m_num_phi + qp;
            //         file_aeh<<","<<aeh_h_inp[q_idx];

            //     }
            // file_aeh<<"]";
            file_aeh << "}" << std::endl;
            file_aeh.close();

        } catch (const std::exception& e) {
            std::cout << "Error occurred while writing AEH data" << std::endl;
            return 0;
        }
    }

    return 0;
}

}  // namespace aeh
