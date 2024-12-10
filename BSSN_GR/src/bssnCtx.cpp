/**
 * @file bssnCtx.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief BSSN contex file.
 * @version 0.1
 * @date 2019-12-20
 *
 * School of Computing, University of Utah.
 * @copyright Copyright (c) 2019
 *
 */

#include "bssnCtx.h"

#include <mpi.h>
#include <sys/types.h>

#include <cstdint>

#include "grUtils.h"
#include "parUtils.h"
#include "parameters.h"

namespace bssn {
BSSNCtx::BSSNCtx(ot::Mesh* pMesh) : Ctx() {
    m_uiMesh = pMesh;

    m_var[VL::CPU_EV].create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST, BSSN_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_IN].create_vector(
        m_uiMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_OUT].create_vector(
        m_uiMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_NUM_VARS, true);

    m_var[VL::CPU_CV].create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST,
                                    BSSN_CONSTRAINT_NUM_VARS, true);
    m_var[VL::CPU_CV_UZ_IN].create_vector(
        m_uiMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_CONSTRAINT_NUM_VARS, true);

    m_uiTinfo._m_uiStep = 0;
    m_uiTinfo._m_uiT    = 0;
    m_uiTinfo._m_uiTb   = bssn::BSSN_RK_TIME_BEGIN;
    m_uiTinfo._m_uiTe   = bssn::BSSN_RK_TIME_END;
    m_uiTinfo._m_uiTh   = bssn::BSSN_RK45_TIME_STEP_SIZE;

    m_uiElementOrder    = bssn::BSSN_ELE_ORDER;

    m_uiMinPt           = Point(bssn::BSSN_GRID_MIN_X, bssn::BSSN_GRID_MIN_Y,
                                bssn::BSSN_GRID_MIN_Z);
    m_uiMaxPt           = Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,
                                bssn::BSSN_GRID_MAX_Z);

    deallocate_bssn_deriv_workspace();
    allocate_bssn_deriv_workspace(m_uiMesh, 1);

    ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS,
                                      BSSN_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS,
                                    BSSN_ASYNC_COMM_K);

    return;
}

BSSNCtx::~BSSNCtx() {
    for (unsigned int i = 0; i < VL::END; i++) m_var[i].destroy_vector();

    deallocate_bssn_deriv_workspace();
    ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS,
                                      BSSN_ASYNC_COMM_K);
}

int BSSNCtx::rhs(DVec* in, DVec* out, unsigned int sz, DendroScalar time) {
    // all the variables should be packed together.
    // assert(sz==1);
    // DendroScalar * sVar[BSSN_NUM_VARS];
    // in->to_2d(sVar);

    this->unzip(*in, m_var[VL::CPU_EV_UZ_IN], bssn::BSSN_ASYNC_COMM_K);

    // AT THIS POINT THE CONSTRAINT VARIABLES SHOULD BE CALCULATED.
    // Just need to set up the proper variable representation.
    DVec& m_cvar_unz = m_var[VL::CPU_CV_UZ_IN];
    DendroScalar* consUnzipVar[bssn::BSSN_CONSTRAINT_NUM_VARS];
    m_cvar_unz.to_2d(consUnzipVar);

#ifdef __PROFILE_CTX__
    this->m_uiCtxpt[ts::CTXPROFILE::RHS].start();
#endif

    DendroScalar* unzipIn[BSSN_NUM_VARS];
    DendroScalar* unzipOut[BSSN_NUM_VARS];

    m_var[CPU_EV_UZ_IN].to_2d(unzipIn);
    m_var[CPU_EV_UZ_OUT].to_2d(unzipOut);

    const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

    bssnRHS(unzipOut, (const DendroScalar**)unzipIn, blkList, numBlocks, time,
            (const DendroScalar**)consUnzipVar);

#ifdef __PROFILE_CTX__
    this->m_uiCtxpt[ts::CTXPROFILE::RHS].stop();
#endif

    this->zip(m_var[CPU_EV_UZ_OUT], *out);

    return 0;
}

#if 0
int BSSNCtx::rhs_blkwise(DVec in, DVec out, const unsigned int* const blkIDs,
                         unsigned int numIds, DendroScalar* blk_time) {
    DendroScalar* unzipIn[BSSN_NUM_VARS];
    DendroScalar* unzipOut[BSSN_NUM_VARS];

    in.to_2d(unzipIn);
    out.to_2d(unzipOut);

    const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;

    for (unsigned int i = 0; i < numIds; i++) {
        const unsigned int blk = blkIDs[i];
        assert(blk < numBlocks);

        offset = blkList[blk].getOffset();
        sz[0] = blkList[blk].getAllocationSzX();
        sz[1] = blkList[blk].getAllocationSzY();
        sz[2] = blkList[blk].getAllocationSzZ();

        bflag = blkList[blk].getBlkNodeFlag();

        dx = blkList[blk].computeDx(pt_min, pt_max);
        dy = blkList[blk].computeDy(pt_min, pt_max);
        dz = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

#ifdef BSSN_RHS_STAGED_COMP
        bssnrhs_sep(unzipOut, (const double**)unzipIn, offset, ptmin, ptmax, sz,
                    bflag);
#else
        bssnrhs(unzipOut, (const double**)unzipIn, offset, ptmin, ptmax, sz,
                bflag);
#endif
    }

    return 0;
}

int BSSNCtx::rhs_blk(const DendroScalar* in, DendroScalar* out,
                     unsigned int dof, unsigned int local_blk_id,
                     DendroScalar blk_time) {
#ifdef __PROFILE_CTX__
    this->m_uiCtxpt[ts::CTXPROFILE::RHS_BLK].start();
#endif

    DendroScalar* unzipIn[dof];
    DendroScalar* unzipOut[dof];

    const unsigned int blk = local_blk_id;

    const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;

    sz[0] = blkList[blk].getAllocationSzX();
    sz[1] = blkList[blk].getAllocationSzY();
    sz[2] = blkList[blk].getAllocationSzZ();

    const unsigned int NN = sz[0] * sz[1] * sz[2];

    for (unsigned int v = 0; v < dof; v++) {
        unzipIn[v] = (DendroScalar*)(in + v * NN);
        unzipOut[v] = (DendroScalar*)(out + v * NN);
    }

    bflag = blkList[blk].getBlkNodeFlag();
    const unsigned int pw = blkList[blk].get1DPadWidth();

    dx = blkList[blk].computeDx(pt_min, pt_max);
    dy = blkList[blk].computeDy(pt_min, pt_max);
    dz = blkList[blk].computeDz(pt_min, pt_max);

    ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
    ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
    ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

    ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
    ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
    ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

#ifdef BSSN_RHS_STAGED_COMP
    bssnrhs_sep(unzipOut, (const DendroScalar**)unzipIn, 0, ptmin, ptmax, sz,
                bflag);
#else
    bssnrhs(unzipOut, (const DendroScalar**)unzipIn, 0, ptmin, ptmax, sz,
            bflag);
#endif

#ifdef __PROFILE_CTX__
    this->m_uiCtxpt[ts::CTXPROFILE::RHS_BLK].stop();
#endif

    return 0;
}

int BSSNCtx::pre_stage_blk(DendroScalar* in, unsigned int dof,
                           unsigned int local_blk_id, DendroScalar blk_time) {
    DendroScalar** unzipIn = new DendroScalar*[dof];
    const unsigned int blk = local_blk_id;

    const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);

    sz[0] = blkList[blk].getAllocationSzX();
    sz[1] = blkList[blk].getAllocationSzY();
    sz[2] = blkList[blk].getAllocationSzZ();

    const unsigned int NN = sz[0] * sz[1] * sz[2];

    for (unsigned int v = 0; v < dof; v++) {
        unzipIn[v] = (DendroScalar*)(in + v * NN);
    }

    bflag = blkList[blk].getBlkNodeFlag();
    const unsigned int pw = blkList[blk].get1DPadWidth();

    if (!bflag) {
        for (unsigned int node = 0; node < NN; node++)
            enforce_bssn_constraints(unzipIn, node);

    } else {
        // note that we can apply enforce bssn constraints in the right padd, at
        // the left boundary block, currently we only apply internal parts of
        // the boundary blocks.
        for (unsigned int k = pw; k < sz[2] - pw; k++)
            for (unsigned int j = pw; j < sz[1] - pw; j++)
                for (unsigned int i = pw; i < sz[0] - pw; i++) {
                    const unsigned nid = k * sz[1] * sz[0] + j * sz[0] + i;
                    enforce_bssn_constraints(unzipIn, nid);
                }
    }

    delete[] unzipIn;
    return 0;
}

int BSSNCtx::post_stage_blk(DendroScalar* in, unsigned int dof,
                            unsigned int local_blk_id, DendroScalar blk_time) {
    return 0;
}

int BSSNCtx::pre_timestep_blk(DendroScalar* in, unsigned int dof,
                              unsigned int local_blk_id,
                              DendroScalar blk_time) {
    return 0;
}

int BSSNCtx::post_timestep_blk(DendroScalar* in, unsigned int dof,
                               unsigned int local_blk_id,
                               DendroScalar blk_time) {
    DendroScalar** unzipIn = new DendroScalar*[dof];
    const unsigned int blk = local_blk_id;

    const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);

    sz[0] = blkList[blk].getAllocationSzX();
    sz[1] = blkList[blk].getAllocationSzY();
    sz[2] = blkList[blk].getAllocationSzZ();

    const unsigned int NN = sz[0] * sz[1] * sz[2];

    for (unsigned int v = 0; v < dof; v++) {
        unzipIn[v] = (DendroScalar*)(in + v * NN);
    }

    bflag = blkList[blk].getBlkNodeFlag();
    const unsigned int pw = blkList[blk].get1DPadWidth();

    if (!bflag) {
        for (unsigned int node = 0; node < NN; node++)
            enforce_bssn_constraints(unzipIn, node);
        // for(unsigned int k=pw; k < sz[2]-pw; k++)
        // for(unsigned int j=pw; j < sz[1]-pw; j++)
        // for(unsigned int i=pw; i < sz[0]-pw; i++)
        // {
        //     const unsigned nid = k*sz[1]*sz[0] + j*sz[0] + i;
        //     enforce_bssn_constraints(unzipIn,nid);
        // }

    } else {
        // note that we can apply enforce bssn constraints in the right padd, at
        // the left boundary block, currently we only apply internal parts of
        // the boundary blocks.
        for (unsigned int k = pw; k < sz[2] - pw; k++)
            for (unsigned int j = pw; j < sz[1] - pw; j++)
                for (unsigned int i = pw; i < sz[0] - pw; i++) {
                    const unsigned nid = k * sz[1] * sz[0] + j * sz[0] + i;
                    enforce_bssn_constraints(unzipIn, nid);
                }
    }

    delete[] unzipIn;
    return 0;
}
#endif

int BSSNCtx::initialize() {
    if (bssn::BSSN_RESTORE_SOLVER) {
        if (!this->restore_checkpt()) {
            return 0;
        }
        if (!m_uiMesh->getMPIRankGlobal()) {
            std::cout << "[BSSNCtx : initialization]  Continuing normal "
                         "initialization due to failed checkpoint restoration."
                      << std::endl;
        }
    }

    this->init_grid();

    bool isRefine = false;
    DendroIntL oldElements, oldElements_g;
    DendroIntL newElements, newElements_g;

    DendroIntL oldGridPoints, oldGridPoints_g;
    DendroIntL newGridPoints, newGridPoints_g;

    unsigned int iterCount         = 1;
    const unsigned int max_iter    = bssn::BSSN_INIT_GRID_ITER;
    const unsigned int rank_global = m_uiMesh->getMPIRankGlobal();
    MPI_Comm gcomm                 = m_uiMesh->getMPIGlobalCommunicator();

    DendroScalar* unzipVar[bssn::BSSN_NUM_VARS];
    unsigned int refineVarIds[bssn::BSSN_NUM_REFINE_VARS];

    for (unsigned int vIndex = 0; vIndex < bssn::BSSN_NUM_REFINE_VARS; vIndex++)
        refineVarIds[vIndex] = bssn::BSSN_REFINE_VARIABLE_INDICES[vIndex];

    double wTol = bssn::BSSN_WAVELET_TOL;
    std::function<double(double, double, double, double* hx)> waveletTolFunc =
        [](double x, double y, double z, double* hx) {
            return bssn::computeWTolDCoords(x, y, z, hx);
        };

    DVec& m_evar     = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];

    // calculate the initial size of the grid
    this->calculate_full_grid_size();

    do {
        // hold on to the "old numbers" for final check
        oldElements_g   = m_uiGlobalMeshElements;
        oldGridPoints_g = m_uiGlobalGridPoints;

        this->unzip(m_evar, m_evar_unz, bssn::BSSN_ASYNC_COMM_K);
        m_evar_unz.to_2d(unzipVar);
        // isRefine=this->is_remesh();
        // NOTE: use a parameter to determine if initial grid covergence should
        // be based on WAMR or the input method
        if (bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE) {
            // is remesh will use whatever the parameter was set to, not the
            // default WAMR though... WAMR may actually be selected! But this
            // way, any of the other initial refinement modes can be called
            // instead during the initial convergence!
            isRefine = this->is_remesh();
        } else {
            // enforce WMAR refinement based refinement initially.
            isRefine =
                bssn::isReMeshWAMR(m_uiMesh, (const double**)unzipVar,
                                   refineVarIds, bssn::BSSN_NUM_REFINE_VARS,
                                   waveletTolFunc, bssn::BSSN_DENDRO_AMR_FAC);
        }

        if (isRefine) {
            ot::Mesh* newMesh =
                this->remesh(bssn::BSSN_DENDRO_GRAIN_SZ,
                             bssn::BSSN_LOAD_IMB_TOL, bssn::BSSN_SPLIT_FIX);

            newElements   = newMesh->getNumLocalMeshElements();
            newGridPoints = newMesh->getNumLocalMeshNodes();

            par::Mpi_Allreduce(&newElements, &newElements_g, 1, MPI_SUM, gcomm);
            par::Mpi_Allreduce(&newGridPoints, &newGridPoints_g, 1, MPI_SUM,
                               gcomm);

            if (!rank_global) {
                std::cout << "[bssnCtx] iter : " << iterCount
                          << " (Remesh triggered) ->  old mesh : "
                          << m_uiGlobalMeshElements
                          << " new mesh : " << newElements_g << std::endl;
                std::cout << "[bssnCtx] iter : " << iterCount
                          << " (Remesh triggered) ->  old mesh (zip nodes) : "
                          << m_uiGlobalGridPoints
                          << " new mesh (zip nodes) : " << newGridPoints_g
                          << std::endl;
            }

            this->grid_transfer(newMesh);

            std::swap(m_uiMesh, newMesh);
            delete newMesh;

            // then update the size of the grid, no need to recompute
            m_uiGlobalMeshElements = newElements_g;
            m_uiGlobalGridPoints   = newGridPoints_g;

#ifdef __CUDACC__
            device::MeshGPU*& dptr_mesh = this->get_meshgpu_device_ptr();
            device::MeshGPU* mesh_gpu   = this->get_meshgpu_host_handle();

            mesh_gpu->dealloc_mesh_on_device(dptr_mesh);
            dptr_mesh = mesh_gpu->alloc_mesh_on_device(m_uiMesh);
#endif
        }

        iterCount += 1;

    } while (isRefine &&
             (newElements_g != oldElements_g ||
              newGridPoints_g != oldGridPoints_g) &&
             (iterCount < max_iter));

    this->init_grid();

    // // realloc bssn deriv space
    deallocate_bssn_deriv_workspace();
    allocate_bssn_deriv_workspace(m_uiMesh, 1);

    unsigned int lmin, lmax;
    m_uiMesh->computeMinMaxLevel(lmin, lmax);

    // calculate the minimum dx
    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));

    bssn::BSSN_RK45_TIME_STEP_SIZE =
        bssn::BSSN_CFL_FACTOR *
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    m_uiTinfo._m_uiTh = bssn::BSSN_RK45_TIME_STEP_SIZE;

    if (bssn::BSSN_SCALE_VTU_AND_GW_EXTRACTION) {
        // REMEMBER: the true max depth of the array is two minus m_uiMaxDepth
        bssn::BSSN_IO_OUTPUT_FREQ_TRUE =
            bssn::BSSN_IO_OUTPUT_FREQ >> (m_uiMaxDepth - 2 - lmax);
        bssn::BSSN_GW_EXTRACT_FREQ_TRUE =
            bssn::BSSN_GW_EXTRACT_FREQ >> (m_uiMaxDepth - 2 - lmax);
    }

    if (!m_uiMesh->getMPIRankGlobal()) {
        const DendroScalar dx_finest =
            ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
             ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
             ((double)(1u << (m_uiMaxDepth))));
        const DendroScalar dt_finest = bssn::BSSN_CFL_FACTOR * dx_finest;

        std::cout << "================= Grid Info (After init grid "
                     "converge):==============================================="
                     "========"
                  << std::endl;
        std::cout << "lmin: " << lmin << " lmax:" << lmax << std::endl;
        std::cout << "dx: " << dx_finest << std::endl;
        std::cout << "dt: " << dt_finest << std::endl;
        std::cout << "VTU IO Output Freq: " << bssn::BSSN_IO_OUTPUT_FREQ_TRUE
                  << std::endl;
        std::cout << "GW IO Output Freq: " << bssn::BSSN_GW_EXTRACT_FREQ_TRUE
                  << std::endl;
        std::cout << "========================================================="
                     "======================================================"
                  << std::endl;
    }

    return 0;
}

int BSSNCtx::init_grid() {
    DVec& m_evar                = m_var[VL::CPU_EV];
    DVec& m_dptr_evar           = m_var[VL::GPU_EV];

    const ot::TreeNode* pNodes  = &(*(m_uiMesh->getAllElements().begin()));
    const unsigned int eleOrder = m_uiMesh->getElementOrder();
    const unsigned int* e2n_cg  = &(*(m_uiMesh->getE2NMapping().begin()));
    const unsigned int* e2n_dg  = &(*(m_uiMesh->getE2NMapping_DG().begin()));
    const unsigned int nPe      = m_uiMesh->getNumNodesPerElement();
    const unsigned int nodeLocalBegin = m_uiMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd   = m_uiMesh->getNodeLocalEnd();

    DendroScalar* zipIn[bssn::BSSN_NUM_VARS];
    m_evar.to_2d(zipIn);

    DendroScalar var1[bssn::BSSN_NUM_VARS];

    DendroScalar mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
    // set the TP communicator.
    if (bssn::BSSN_ID_TYPE == 0) {
        TP_MPI_COMM = m_uiMesh->getMPIGlobalCommunicator();
        TwoPunctures((double)0, (double)0, (double)0, var1, &mp, &mm, &mp_adm,
                     &mm_adm, &E, &J1, &J2, &J3);
    }

    for (unsigned int elem = m_uiMesh->getElementLocalBegin();
         elem < m_uiMesh->getElementLocalEnd(); elem++) {
        DendroScalar var[bssn::BSSN_NUM_VARS];
        for (unsigned int k = 0; k < (eleOrder + 1); k++)
            for (unsigned int j = 0; j < (eleOrder + 1); j++)
                for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                    const unsigned int nodeLookUp_CG =
                        e2n_cg[elem * nPe +
                               k * (eleOrder + 1) * (eleOrder + 1) +
                               j * (eleOrder + 1) + i];
                    if (nodeLookUp_CG >= nodeLocalBegin &&
                        nodeLookUp_CG < nodeLocalEnd) {
                        const unsigned int nodeLookUp_DG =
                            e2n_dg[elem * nPe +
                                   k * (eleOrder + 1) * (eleOrder + 1) +
                                   j * (eleOrder + 1) + i];
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

                        if (bssn::BSSN_ID_TYPE == 0) {
                            const DendroScalar xx = GRIDX_TO_X(x);
                            const DendroScalar yy = GRIDY_TO_Y(y);
                            const DendroScalar zz = GRIDZ_TO_Z(z);

                            TwoPunctures((double)xx, (double)yy, (double)zz,
                                         var, &mp, &mm, &mp_adm, &mm_adm, &E,
                                         &J1, &J2, &J3);
                        } else {
                            // all other values are handled in the initial data
                            // wrapper including an error message
                            initialDataFunctionWrapper((double)x, (double)y,
                                                       (double)z, var);
                        }
                        for (unsigned int v = 0; v < bssn::BSSN_NUM_VARS; v++)
                            zipIn[v][nodeLookUp_CG] = var[v];
                    }
                }
    }

    for (unsigned int node = m_uiMesh->getNodeLocalBegin();
         node < m_uiMesh->getNodeLocalEnd(); node++)
        enforce_bssn_constraints(zipIn, node);

#ifdef BSSN_EXTRACT_BH_LOCATIONS
    m_uiBHLoc[0] = Point(bssn::BH1.getBHCoordX(), bssn::BH1.getBHCoordY(),
                         bssn::BH1.getBHCoordZ());
    m_uiBHLoc[1] = Point(bssn::BH2.getBHCoordX(), bssn::BH2.getBHCoordY(),
                         bssn::BH2.getBHCoordZ());
#endif

    return 0;
}

int BSSNCtx::finalize() { return 0; }

void BSSNCtx::compute_constraint_variables() {
    if (m_bConstraintsComputed) return;

    DVec& m_evar     = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];
    DVec& m_cvar     = m_var[VL::CPU_CV];
    DVec& m_cvar_unz = m_var[VL::CPU_CV_UZ_IN];

    this->unzip(m_evar, m_evar_unz, BSSN_ASYNC_COMM_K);

    DendroScalar* consUnzipVar[bssn::BSSN_CONSTRAINT_NUM_VARS];
    DendroScalar* consVar[bssn::BSSN_CONSTRAINT_NUM_VARS];

    DendroScalar* evolUnzipVar[bssn::BSSN_NUM_VARS];
    DendroScalar* evolVar[bssn::BSSN_NUM_VARS];

    m_evar_unz.to_2d(evolUnzipVar);
    m_cvar_unz.to_2d(consUnzipVar);

    m_evar.to_2d(evolVar);
    m_cvar.to_2d(consVar);

#if BSSN_COMPUTE_CONSTRAINTS

    if (!(m_uiMesh->getMPIRankGlobal())) {
        std::cout << BLU << "[BSSN] - Now computing constraints" << NRM
                  << std::endl;
    }

    const std::vector<ot::Block> blkList = m_uiMesh->getLocalBlockList();

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;

    for (unsigned int blk = 0; blk < blkList.size(); blk++) {
        offset   = blkList[blk].getOffset();
        sz[0]    = blkList[blk].getAllocationSzX();
        sz[1]    = blkList[blk].getAllocationSzY();
        sz[2]    = blkList[blk].getAllocationSzZ();

        bflag    = blkList[blk].getBlkNodeFlag();

        dx       = blkList[blk].computeDx(pt_min, pt_max);
        dy       = blkList[blk].computeDy(pt_min, pt_max);
        dz       = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

        physical_constraints(consUnzipVar, (const DendroScalar**)evolUnzipVar,
                             offset, ptmin, ptmax, sz, bflag);
    }

    /*double consVecMin[bssn::BSSN_CONSTRAINT_NUM_VARS];
    double consVecMax[bssn::BSSN_CONSTRAINT_NUM_VARS];*/
    double constraintMaskedL2[bssn::BSSN_CONSTRAINT_NUM_VARS];
    this->zip(m_cvar_unz, m_cvar);
    m_uiMesh->readFromGhostBegin(m_cvar.get_vec_ptr(), m_cvar.get_dof());
    m_uiMesh->readFromGhostEnd(m_cvar.get_vec_ptr(), m_cvar.get_dof());

    if (!(m_uiMesh->getMPIRankGlobal())) {
        std::cout << BLU << "[BSSN] - Finished computing constraints!" << NRM
                  << std::endl;
    }

#endif

    m_bConstraintsComputed = true;
}

int BSSNCtx::extract_constraints() {
    // make sure constraints are computed
    this->compute_constraint_variables();

    // variables should be ghost-exchanged and unzipped thanks to
    // compute_constraint_variables

    DVec& m_evar     = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];
    DVec& m_cvar     = m_var[VL::CPU_CV];
    DVec& m_cvar_unz = m_var[VL::CPU_CV_UZ_IN];

    DendroScalar* evolVar[bssn::BSSN_NUM_VARS];
    DendroScalar* consVar[bssn::BSSN_CONSTRAINT_NUM_VARS];

    m_evar.to_2d(evolVar);
    m_cvar.to_2d(consVar);

    bssn::extractConstraints(m_uiMesh, (const DendroScalar**)consVar,
                             evolVar[BHLOC::EXTRACTION_VAR_ID],
                             BHLOC::EXTRACTION_TOL, m_uiTinfo._m_uiStep,
                             m_uiTinfo._m_uiT);
    return 0;
}

int BSSNCtx::extract_gravitational_waves() {
#ifndef BSSN_KERR_SCHILD_TEST
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
    // make sure constraints are compuated
    this->compute_constraint_variables();

    DVec& m_cvar = m_var[VL::CPU_CV];
    DendroScalar* consVar[bssn::BSSN_CONSTRAINT_NUM_VARS];
    m_cvar.to_2d(consVar);

    GW::extractFarFieldPsi4(m_uiMesh, (const DendroScalar**)consVar,
                            m_uiTinfo._m_uiStep, m_uiTinfo._m_uiT);
#endif
#endif
    return 0;
}

int BSSNCtx::write_vtu() {
    if (!m_uiMesh->isActive()) return 0;

    DVec& m_evar     = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];
    DVec& m_cvar     = m_var[VL::CPU_CV];
    DVec& m_cvar_unz = m_var[VL::CPU_CV_UZ_IN];

    this->unzip(m_evar, m_evar_unz, BSSN_ASYNC_COMM_K);

    DendroScalar* consUnzipVar[bssn::BSSN_CONSTRAINT_NUM_VARS];
    DendroScalar* consVar[bssn::BSSN_CONSTRAINT_NUM_VARS];

    DendroScalar* evolUnzipVar[bssn::BSSN_NUM_VARS];
    DendroScalar* evolVar[bssn::BSSN_NUM_VARS];

    m_evar_unz.to_2d(evolUnzipVar);
    m_cvar_unz.to_2d(consUnzipVar);

    m_evar.to_2d(evolVar);
    m_cvar.to_2d(consVar);

    // make sure the constraint variables are computed, they should be, but at
    // least it'll early exit if they are this time step
    this->compute_constraint_variables();

#ifdef BSSN_ENABLE_VTU_OUTPUT
    if (!(m_uiMesh->getMPIRankGlobal())) {
        std::cout << GRN << "=== Now Writing VTU Output Files! ===" << NRM
                  << std::endl;
    }

    std::vector<std::string> pDataNames;
    const unsigned int numConstVars = bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT;
    const unsigned int numEvolVars  = bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT;

    double* pData[(numConstVars + numEvolVars)];

    for (unsigned int i = 0; i < numEvolVars; i++) {
        pDataNames.push_back(
            std::string(bssn::BSSN_VAR_NAMES[BSSN_VTU_OUTPUT_EVOL_INDICES[i]]));
        pData[i] = evolVar[BSSN_VTU_OUTPUT_EVOL_INDICES[i]];
    }

    for (unsigned int i = 0; i < numConstVars; i++) {
        pDataNames.push_back(std::string(
            bssn::BSSN_CONSTRAINT_VAR_NAMES[BSSN_VTU_OUTPUT_CONST_INDICES[i]]));
        pData[numEvolVars + i] = consVar[BSSN_VTU_OUTPUT_CONST_INDICES[i]];
    }

    std::vector<char*> pDataNames_char;
    pDataNames_char.reserve(pDataNames.size());

    for (unsigned int i = 0; i < pDataNames.size(); i++)
        pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

    const char* fDataNames[] = {"Time", "Cycle"};
    const double fData[]     = {m_uiTinfo._m_uiT, (double)m_uiTinfo._m_uiStep};

    char fPrefix[256];
    sprintf(fPrefix, "%s_%d", bssn::BSSN_VTU_FILE_PREFIX.c_str(),
            m_uiTinfo._m_uiStep);

    if (bssn::BSSN_VTU_Z_SLICE_ONLY) {
        unsigned int s_val[3]  = {1u << (m_uiMaxDepth - 1),
                                  1u << (m_uiMaxDepth - 1),
                                  1u << (m_uiMaxDepth - 1)};
        unsigned int s_norm[3] = {0, 0, 1};
        io::vtk::mesh2vtu_slice(m_uiMesh, s_val, s_norm, fPrefix, 2, fDataNames,
                                fData, (numEvolVars + numConstVars),
                                (const char**)&pDataNames_char[0],
                                (const double**)pData);
    } else
        io::vtk::mesh2vtuFine(m_uiMesh, fPrefix, 2, fDataNames, fData,
                              (numEvolVars + numConstVars),
                              (const char**)&pDataNames_char[0],
                              (const double**)pData);

    if (!(m_uiMesh->getMPIRankGlobal())) {
        std::cout << GRN << "=== Finished Writing the VTU Files! ===" << NRM
                  << std::endl;
    }

#endif

    return 0;
}

int BSSNCtx::write_bh_coords() {
#ifdef BSSN_EXTRACT_BH_LOCATIONS
    // ensure that the black hole locations have been updated, it does early
    // exit if it's happened
    this->evolve_bh_loc();
    bssn::writeBHCoordinates((const ot::Mesh*)m_uiMesh, (const Point*)m_uiBHLoc,
                             2, m_uiTinfo._m_uiStep, m_uiTinfo._m_uiT);
#endif
    return 0;
}

int BSSNCtx::write_checkpt() {
    if (!m_uiMesh->isActive()) return 0;

    // every other checkpoint index should be 0 or 1, this allows "alternate"
    // file writing
    unsigned int cpIndex = (m_uiTinfo._m_uiStep / bssn::BSSN_CHECKPT_FREQ) % 2;

    const bool is_merged =
        ((bssn::BSSN_BH_LOC[0] - bssn::BSSN_BH_LOC[1]).abs() < 0.1);
    if (is_merged && !bssn::BSSN_MERGED_CHKPT_WRITTEN) {
        cpIndex                         = 3;
        bssn::BSSN_MERGED_CHKPT_WRITTEN = true;
    }

    unsigned int rank = m_uiMesh->getMPIRank();
    unsigned int npes = m_uiMesh->getMPICommSize();

    DendroScalar* eVar[BSSN_NUM_VARS];
    DVec& m_evar = m_var[VL::CPU_EV];
    m_evar.to_2d(eVar);

    char fName[256];
    const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin() +
                                     m_uiMesh->getElementLocalBegin()));
    sprintf(fName, "%s_%d_octree_%d.oct", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),
            cpIndex, rank);
    io::checkpoint::writeOctToFile(fName, pNodes,
                                   m_uiMesh->getNumLocalMeshElements());

    unsigned int numVars  = bssn::BSSN_NUM_VARS;
    const char** varNames = bssn::BSSN_VAR_NAMES;

    sprintf(fName, "%s_%d_%d.var", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),
            cpIndex, rank);
    io::checkpoint::writeVecToFile(fName, m_uiMesh, (const double**)eVar,
                                   bssn::BSSN_NUM_VARS);

    if (!rank)
        std::cout << "   ...Finished writing the octree checkpoints!"
                  << std::endl;

    if (!rank) {
        sprintf(fName, "%s_%d_step.cp", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),
                cpIndex);
        std::cout << "[BSSNCtx] \t writing checkpoint file : " << fName
                  << std::endl;
        std::ofstream outfile(fName);
        if (!outfile) {
            std::cout << fName << " file open failed " << std::endl;
            return 0;
        }

        json checkPoint;
        checkPoint["DENDRO_TS_TIME_BEGIN"]         = m_uiTinfo._m_uiTb;
        checkPoint["DENDRO_TS_TIME_END"]           = m_uiTinfo._m_uiTe;
        checkPoint["DENDRO_TS_ELEMENT_ORDER"]      = m_uiElementOrder;

        checkPoint["DENDRO_TS_TIME_CURRENT"]       = m_uiTinfo._m_uiT;
        checkPoint["DENDRO_TS_STEP_CURRENT"]       = m_uiTinfo._m_uiStep;
        checkPoint["DENDRO_TS_TIME_STEP_SIZE"]     = m_uiTinfo._m_uiTh;
        checkPoint["DENDRO_TS_LAST_IO_TIME"]       = m_uiTinfo._m_uiT;

        checkPoint["DENDRO_TS_WAVELET_TOLERANCE"]  = bssn::BSSN_WAVELET_TOL;
        checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"] = bssn::BSSN_LOAD_IMB_TOL;
        checkPoint["DENDRO_TS_NUM_VARS"] =
            numVars;  // number of variables to restore.
        checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"] =
            m_uiMesh->getMPICommSize();  // (note that rank 0 is always active).

        checkPoint["DENDRO_BH1_X"]              = m_uiBHLoc[0].x();
        checkPoint["DENDRO_BH1_Y"]              = m_uiBHLoc[0].y();
        checkPoint["DENDRO_BH1_Z"]              = m_uiBHLoc[0].z();

        checkPoint["DENDRO_BH2_X"]              = m_uiBHLoc[1].x();
        checkPoint["DENDRO_BH2_Y"]              = m_uiBHLoc[1].y();
        checkPoint["DENDRO_BH2_Z"]              = m_uiBHLoc[1].z();

        checkPoint["DENDRO_BSSN_BH_MERGE"]      = m_bIsBHMerged;
        checkPoint["DENDRO_BSSN_BH_MERGE_TIME"] = m_dMergeTime;
        checkPoint["DENDRO_BSSN_BH_MERGE_STEP"] = m_uiMergeStep;

        // then also the BH time history
        checkPoint["DENDRO_BSSN_BH_LOC_TIMES"]  = m_uiBHTimeHistory;
        // then the X, Y, Z points for the BSSN BH locations

        std::tuple<std::string, std::string, std::string> temp =
            encode_bh_locs(m_uiBHLocHistory, m_uiBHTimeHistory);
        auto [bh1_str, bh2_str, time_str]   = temp;

        checkPoint["DENDRO_BSSN_BH_LOC_T"]  = time_str;
        checkPoint["DENDRO_BSSN_BH_LOC_B1"] = bh1_str;
        checkPoint["DENDRO_BSSN_BH_LOC_B2"] = bh2_str;

        outfile << std::setw(4) << checkPoint << std::endl;
        outfile.close();

        std::cout << "   ...Finished writing the plain-text checkpoint!"
                  << std::endl;
        std::cout << "[BSSNCtx] \t finished writing checkpoint! " << std::endl;
    }

    return 0;
}

int BSSNCtx::restore_checkpt() {
    unsigned int numVars = 0;
    std::vector<ot::TreeNode> octree;
    json checkPoint;

    int rank;
    int npes;
    MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    unsigned int activeCommSz;

    char fName[256];

    // 0 -> everything is fine
    // 1 -> file failed to open
    // 2 -> file doesn't exist
    unsigned int restoreStatus = 0;
    unsigned int restoreStatusGlobal =
        0;  // 0 indicates successfully restorable.

    ot::Mesh* newMesh;
    unsigned int restoreStep[2];
    restoreStep[0]                = 0;
    restoreStep[1]                = 0;

    unsigned int restoreFileIndex = 0;

    for (unsigned int cpIndex = 0; cpIndex < 2; cpIndex++) {
        restoreStatus = 0;

        if (!rank) {
            sprintf(fName, "%s_%d_step.cp",
                    bssn::BSSN_CHKPT_FILE_PREFIX.c_str(), cpIndex);

            std::cout << "    Checking to see if " << fName
                      << " is a valid checkpoint file..." << std::endl;

            if (!std::filesystem::exists(fName)) {
                // std::cout << YLW << "WARNING: " << NRM << "Checkpoint
                // filename " << fName << " does not exist!" << std::endl;
                restoreStatus = 2;
                continue;
            }

            std::ifstream infile(fName);
            if (!infile) {
                std::cout << fName << " file open failed " << std::endl;
                restoreStatus = 1;
            }

            if (restoreStatus == 0) {
                infile >> checkPoint;
                m_uiTinfo._m_uiTb   = checkPoint["DENDRO_TS_TIME_BEGIN"];
                m_uiTinfo._m_uiTe   = checkPoint["DENDRO_TS_TIME_END"];
                m_uiTinfo._m_uiT    = checkPoint["DENDRO_TS_TIME_CURRENT"];
                m_uiTinfo._m_uiStep = checkPoint["DENDRO_TS_STEP_CURRENT"];
                m_uiTinfo._m_uiTh   = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                m_uiElementOrder    = checkPoint["DENDRO_TS_ELEMENT_ORDER"];

                bssn::BSSN_WAVELET_TOL =
                    checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                bssn::BSSN_LOAD_IMB_TOL =
                    checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];

                numVars      = checkPoint["DENDRO_TS_NUM_VARS"];
                activeCommSz = checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];

                m_uiBHLoc[0] = Point((double)checkPoint["DENDRO_BH1_X"],
                                     (double)checkPoint["DENDRO_BH1_Y"],
                                     (double)checkPoint["DENDRO_BH1_Z"]);
                m_uiBHLoc[1] = Point((double)checkPoint["DENDRO_BH2_X"],
                                     (double)checkPoint["DENDRO_BH2_Y"],
                                     (double)checkPoint["DENDRO_BH2_Z"]);

                // if this key is in, then all three keys should be
                if (checkPoint.find("DENDRO_BSSN_BH_MERGE") !=
                    checkPoint.end()) {
                    // restore bh merge and merge time information
                    m_bIsBHMerged = checkPoint["DENDRO_BSSN_BH_MERGE"];
                    m_dMergeTime  = checkPoint["DENDRO_BSSN_BH_MERGE_TIME"];
                    m_uiMergeStep = checkPoint["DENDRO_BSSN_BH_MERGE_STEP"];
                }

                // OLD: DEPRECIATED
                if (checkPoint.find("DENDRO_BSSN_BH_LOC_TIMES") !=
                    checkPoint.end()) {
                    // restore bh location data

                    // clear bhLocHistory vector to start fresh
                    m_uiBHLocHistory.clear();

                    for (const auto& pair_json :
                         checkPoint["DENDRO_BSSN_BH_LOC_HISTORY"]) {
                        Point bh1pt =
                            Point(pair_json["bh1"]["x"].get<double>(),
                                  pair_json["bh1"]["y"].get<double>(),
                                  pair_json["bh1"]["z"].get<double>());
                        Point bh2pt =
                            Point(pair_json["bh2"]["x"].get<double>(),
                                  pair_json["bh2"]["y"].get<double>(),
                                  pair_json["bh2"]["z"].get<double>());

                        m_uiBHLocHistory.emplace_back(bh1pt, bh2pt);
                    }

                    // then restore the times the bh's were output
                    m_uiBHTimeHistory = checkPoint["DENDRO_BSSN_BH_LOC_TIMES"]
                                            .get<std::vector<double>>();
                }

                if (checkPoint.find("DENDRO_BSSN_BH_LOC_T") !=
                    checkPoint.end()) {
                    auto bh_decoded = decode_bh_locs(
                        checkPoint["DENDRO_BSSN_BH_LOC_B1"].get<std::string>(),
                        checkPoint["DENDRO_BSSN_BH_LOC_B2"].get<std::string>(),
                        checkPoint["DENDRO_BSSN_BH_LOC_T"].get<std::string>());

                    m_uiBHLocHistory.clear();

                    m_uiBHLocHistory  = std::get<0>(bh_decoded);
                    m_uiBHTimeHistory = std::get<1>(bh_decoded);
                }

                restoreStep[cpIndex] = m_uiTinfo._m_uiStep;
            }
        }
    }

    if (!rank) {
        if (restoreStep[0] < restoreStep[1])
            restoreFileIndex = 1;
        else
            restoreFileIndex = 0;
    }

    par::Mpi_Bcast(&restoreFileIndex, 1, 0, comm);

    restoreStatus = 0;
    octree.clear();
    if (!rank)
        std::cout << "[BSSNCtx] :  Trying to restore from checkpoint index : "
                  << restoreFileIndex << std::endl;

    // now we do the "true" restore, where every process knows which one is the
    // right one
    sprintf(fName, "%s_%d_step.cp", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),
            restoreFileIndex);

    // first check to see if the file even exists
    if (!std::filesystem::exists(fName)) {
        std::cout << YLW << "WARNING: " << NRM << "Checkpoint filename "
                  << fName << " does not exist!" << std::endl;
        restoreStatus = 2;
    }

    if (restoreStatus == 0) {
        std::ifstream infile(fName);
        if (!infile) {
            std::cout << fName << " file open failed " << std::endl;
            restoreStatus = 1;
        }

        if (restoreStatus == 0) {
            infile >> checkPoint;
            m_uiTinfo._m_uiTb      = checkPoint["DENDRO_TS_TIME_BEGIN"];
            m_uiTinfo._m_uiTe      = checkPoint["DENDRO_TS_TIME_END"];
            m_uiTinfo._m_uiT       = checkPoint["DENDRO_TS_TIME_CURRENT"];
            m_uiTinfo._m_uiStep    = checkPoint["DENDRO_TS_STEP_CURRENT"];
            m_uiTinfo._m_uiTh      = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
            m_uiElementOrder       = checkPoint["DENDRO_TS_ELEMENT_ORDER"];

            bssn::BSSN_WAVELET_TOL = checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
            bssn::BSSN_LOAD_IMB_TOL =
                checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];

            numVars      = checkPoint["DENDRO_TS_NUM_VARS"];
            activeCommSz = checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];

            m_uiBHLoc[0] = Point((double)checkPoint["DENDRO_BH1_X"],
                                 (double)checkPoint["DENDRO_BH1_Y"],
                                 (double)checkPoint["DENDRO_BH1_Z"]);
            m_uiBHLoc[1] = Point((double)checkPoint["DENDRO_BH2_X"],
                                 (double)checkPoint["DENDRO_BH2_Y"],
                                 (double)checkPoint["DENDRO_BH2_Z"]);

            // if this key is in, then all three keys should be
            if (checkPoint.find("DENDRO_BSSN_BH_MERGE") != checkPoint.end()) {
                // restore bh merge and merge time information
                m_bIsBHMerged = checkPoint["DENDRO_BSSN_BH_MERGE"];
                m_dMergeTime  = checkPoint["DENDRO_BSSN_BH_MERGE_TIME"];
                m_uiMergeStep = checkPoint["DENDRO_BSSN_BH_MERGE_STEP"];
            }

            if (checkPoint.find("DENDRO_BSSN_BH_LOC_TIMES") !=
                checkPoint.end()) {
                // restore bh location data

                // make sure BHLocHistory is completely empty upon reading!
                m_uiBHLocHistory.clear();

                for (const auto& pair_json :
                     checkPoint["DENDRO_BSSN_BH_LOC_HISTORY"]) {
                    Point bh1pt = Point(pair_json["bh1"]["x"].get<double>(),
                                        pair_json["bh1"]["y"].get<double>(),
                                        pair_json["bh1"]["z"].get<double>());
                    Point bh2pt = Point(pair_json["bh2"]["x"].get<double>(),
                                        pair_json["bh2"]["y"].get<double>(),
                                        pair_json["bh2"]["z"].get<double>());

                    m_uiBHLocHistory.emplace_back(bh1pt, bh2pt);
                }

                // then restore the times the bh's were output
                m_uiBHTimeHistory = checkPoint["DENDRO_BSSN_BH_LOC_TIMES"]
                                        .get<std::vector<double>>();
            }

            if (checkPoint.find("DENDRO_BSSN_BH_LOC_T") != checkPoint.end()) {
                auto bh_decoded = decode_bh_locs(
                    checkPoint["DENDRO_BSSN_BH_LOC_B1"].get<std::string>(),
                    checkPoint["DENDRO_BSSN_BH_LOC_B2"].get<std::string>(),
                    checkPoint["DENDRO_BSSN_BH_LOC_T"].get<std::string>());

                m_uiBHLocHistory.clear();

                m_uiBHLocHistory  = std::get<0>(bh_decoded);
                m_uiBHTimeHistory = std::get<1>(bh_decoded);
            }

            restoreStep[restoreFileIndex] = m_uiTinfo._m_uiStep;
        }
    }

    par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX, comm);
    if (restoreStatusGlobal == 1) {
        if (!rank)
            std::cout
                << "[BSSNCtx] : Restore step failed, restore file corrupted. "
                << std::endl;
        MPI_Abort(comm, 0);
    } else if (restoreStatusGlobal == 2) {
        if (!rank) {
            std::cout << "[BSSNCtx] : " << YLW << "WARNING:" << NRM
                      << " No checkpoint files could be found. The program "
                      << GRN << "will continue" << NRM
                      << " as BSSN_RESTORE_SOLVER were false!" << std::endl;
        }
        return 2;
    }

    bssn::BSSN_BH_LOC[0] = m_uiBHLoc[0];
    bssn::BSSN_BH_LOC[1] = m_uiBHLoc[1];

    if (activeCommSz > npes) {
        if (!rank)
            std::cout
                << " [BSSNCtx] : checkpoint file written from  a larger "
                   "communicator than the current global comm. (i.e. "
                   "communicator shrinking not allowed in the restore step. )"
                << std::endl;

        MPI_Abort(comm, 0);
    }

    bool isActive = (rank < activeCommSz);

    MPI_Comm newComm;
    par::splitComm2way(isActive, &newComm, comm);

    if (isActive) {
        int activeRank;
        int activeNpes;

        MPI_Comm_rank(newComm, &activeRank);
        MPI_Comm_size(newComm, &activeNpes);
        assert(activeNpes == activeCommSz);

        sprintf(fName, "%s_%d_octree_%d.oct",
                bssn::BSSN_CHKPT_FILE_PREFIX.c_str(), restoreFileIndex,
                activeRank);
        restoreStatus = io::checkpoint::readOctFromFile(fName, octree);
        assert(par::test::isUniqueAndSorted(octree, newComm));
    }

    par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX, comm);
    if (restoreStatusGlobal == 1) {
        if (!rank)
            std::cout << "[BSSNCtx]: octree (*.oct) restore file is corrupted "
                      << std::endl;
        MPI_Abort(comm, 0);
    }

    newMesh = new ot::Mesh(octree, 1, m_uiElementOrder, activeCommSz, comm);
    newMesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X, bssn::BSSN_GRID_MIN_Y,
                                   bssn::BSSN_GRID_MIN_Z),
                             Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,
                                   bssn::BSSN_GRID_MAX_Z));
    // no need to transfer data only to resize the contex variables.
    // this->grid_transfer(newMesh);
    for (unsigned int i = 0; i < VL::END; i++) m_var[i].destroy_vector();

    m_var[VL::CPU_EV].create_vector(newMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST, BSSN_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_IN].create_vector(
        newMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_OUT].create_vector(
        newMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_NUM_VARS, true);

    m_var[VL::CPU_CV].create_vector(newMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST,
                                    BSSN_CONSTRAINT_NUM_VARS, true);
    m_var[VL::CPU_CV_UZ_IN].create_vector(
        newMesh, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_CONSTRAINT_NUM_VARS, true);

    ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS,
                                      BSSN_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(newMesh, m_mpi_ctx, BSSN_NUM_VARS,
                                    BSSN_ASYNC_COMM_K);

    // only reads the evolution variables.
    if (isActive) {
        int activeRank;
        int activeNpes;

        DendroScalar* inVec[BSSN_NUM_VARS];
        DVec& m_evar = m_var[VL::CPU_EV];
        m_evar.to_2d(inVec);

        MPI_Comm_rank(newComm, &activeRank);
        MPI_Comm_size(newComm, &activeNpes);
        assert(activeNpes == activeCommSz);

        sprintf(fName, "%s_%d_%d.var", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),
                restoreFileIndex, activeRank);
        restoreStatus = io::checkpoint::readVecFromFile(fName, newMesh, inVec,
                                                        bssn::BSSN_NUM_VARS);
    }

    MPI_Comm_free(&newComm);
    par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX, comm);
    if (restoreStatusGlobal == 1) {
        if (!rank)
            std::cout << "[BSSNCtx]: varible (*.var) restore file currupted "
                      << std::endl;
        MPI_Abort(comm, 0);
    }

    std::swap(m_uiMesh, newMesh);
    delete newMesh;

    // realloc bssn deriv space
    deallocate_bssn_deriv_workspace();
    allocate_bssn_deriv_workspace(m_uiMesh, 1);

    unsigned int localSz    = m_uiMesh->getNumLocalMeshElements();
    unsigned int totalElems = 0;
    par::Mpi_Allreduce(&localSz, &totalElems, 1, MPI_SUM, comm);

    // NOTE: this chunk is only needed to restore the minimum dx size that's
    // needed in some computations. Since the initialization ends as soon as the
    // restore is complete, it's important to recalculate this since it won't be
    // called again until a remesh
    unsigned int lmin, lmax;
    m_uiMesh->computeMinMaxLevel(lmin, lmax);
    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));

    if (!rank) {
        std::cout << GRN << "---------------------------------------"
                  << std::endl;
        std::cout << "CHECKPOINT RESTORE SUCCESSFUL: " << NRM << std::endl;
        std::cout << "   checkpoint at step : " << m_uiTinfo._m_uiStep
                  << " | active Comm. sz: " << activeCommSz
                  << " | restored mesh size: " << totalElems << std::endl
                  << std::endl;
        std::cout << " restored mesh min dx: " << bssn::BSSN_CURRENT_MIN_DX
                  << std::endl;
        std::cout << GRN << "---------------------------------------" << NRM
                  << std::endl;
    }

    m_uiIsETSSynced = false;
    return 0;
}

int BSSNCtx::post_stage(DVec& sIn) { return 0; }

int BSSNCtx::post_timestep(DVec& sIn) {
    DendroScalar* evar[BSSN_NUM_VARS];
    sIn.to_2d(evar);
    for (unsigned int node = m_uiMesh->getNodeLocalBegin();
         node < m_uiMesh->getNodeLocalEnd(); node++)
        enforce_bssn_constraints(evar, node);

    return 0;
}

bool BSSNCtx::is_remesh() {
    bool isRefine = false;
    if (bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY) return false;

    MPI_Comm comm    = m_uiMesh->getMPIGlobalCommunicator();

    DVec& m_evar     = m_var[VL::CPU_EV];
    DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];

    this->unzip(m_evar, m_evar_unz, bssn::BSSN_ASYNC_COMM_K);

    DendroScalar* unzipVar[BSSN_NUM_VARS];
    m_evar_unz.to_2d(unzipVar);

    unsigned int refineVarIds[bssn::BSSN_NUM_REFINE_VARS];
    for (unsigned int vIndex = 0; vIndex < bssn::BSSN_NUM_REFINE_VARS; vIndex++)
        refineVarIds[vIndex] = bssn::BSSN_REFINE_VARIABLE_INDICES[vIndex];

    double wTol = bssn::BSSN_WAVELET_TOL;
    std::function<double(double, double, double, double* hx)> waveletTolFunc =
        [](double x, double y, double z, double* hx) {
            return bssn::computeWTolDCoords(x, y, z, hx);
        };

    if (bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::WAMR) {
        isRefine =
            bssn::isReMeshWAMR(m_uiMesh, (const double**)unzipVar, refineVarIds,
                               bssn::BSSN_NUM_REFINE_VARS, waveletTolFunc,
                               bssn::BSSN_DENDRO_AMR_FAC);

    } else if (bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::EH) {
        isRefine = bssn::isRemeshEH(
            m_uiMesh, (const double**)unzipVar, bssn::VAR::U_ALPHA,
            bssn::BSSN_EH_REFINE_VAL, bssn::BSSN_EH_COARSEN_VAL, true);

    } else if (bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::EH_WAMR) {
        const bool isR1 =
            bssn::isReMeshWAMR(m_uiMesh, (const double**)unzipVar, refineVarIds,
                               bssn::BSSN_NUM_REFINE_VARS, waveletTolFunc,
                               bssn::BSSN_DENDRO_AMR_FAC);
        const bool isR2 = bssn::isRemeshEH(
            m_uiMesh, (const double**)unzipVar, bssn::VAR::U_ALPHA,
            bssn::BSSN_EH_REFINE_VAL, bssn::BSSN_EH_COARSEN_VAL, false);

        isRefine = (isR1 || isR2);
    } else if (bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::BH_LOC) {
        isRefine = bssn::isRemeshBH(m_uiMesh, m_uiBHLoc);
    } else if (bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::BH_WAMR) {
        // BHLB refinement for baseline
        const bool isR1 = bssn::isRemeshBH(m_uiMesh, m_uiBHLoc);
        // WAMR for additional refinement
        const bool isR2 =
            bssn::addRemeshWAMR(m_uiMesh, (const double**)unzipVar,
                                refineVarIds, bssn::BSSN_NUM_REFINE_VARS,
                                waveletTolFunc, bssn::BSSN_DENDRO_AMR_FAC);

        isRefine = (isR1 || isR2);
    }

    return isRefine;
}

DVec& BSSNCtx::get_evolution_vars() { return m_var[CPU_EV]; }

DVec& BSSNCtx::get_constraint_vars() { return m_var[CPU_CV]; }

int BSSNCtx::terminal_output() {
    if (m_uiMesh->isActive()) {
        DendroScalar min = 0, max = 0;
        DVec& m_evar = m_var[VL::CPU_EV];
        min = vecMin(m_uiMesh, m_evar.get_vec_ptr(), ot::VEC_TYPE::CG_NODAL,
                     true);
        max = vecMax(m_uiMesh, m_evar.get_vec_ptr(), ot::VEC_TYPE::CG_NODAL,
                     true);

        if (!(m_uiMesh->getMPIRankGlobal())) {
            std::cout << "[BSSNCtx]:  "
                      << bssn::BSSN_VAR_NAMES[bssn::VAR::U_ALPHA]
                      << " (min,max) : \t ( " << min << ", " << max << " ) "
                      << std::endl;
            if (std::isnan(min) || std::isnan(max)) {
                std::cout << "[Error]: NAN detected " << std::endl;
                MPI_Abort(m_uiMesh->getMPICommunicator(), 0);
            }
        }
    }

    return 0;
}

void BSSNCtx::write_grid_summary_data() {
    if (m_uiMesh->isActive()) {
        if (!m_uiMesh->getMPIRankGlobal()) {
            std::string fname =
                bssn::BSSN_PROFILE_FILE_PREFIX + "_GridInfo.dat";
            try {
                std::ofstream file_grid_data;
                file_grid_data.open(fname, std::ofstream::app);
                file_grid_data.precision(12);
                file_grid_data << std::scientific;

                if (!m_uiWroteGridInfoHeader) {
                    file_grid_data << "timeStep,simTime,commSize,wTime,"
                                      "meshSize,totalGridPoints,stepSize\n";

                    m_uiWroteGridInfoHeader = true;
                }

                file_grid_data << bssn::BSSN_CURRENT_RK_STEP << ",";
                file_grid_data << bssn::BSSN_CURRENT_RK_COORD_TIME << ",";
                file_grid_data << m_uiMesh->getMPICommSize() << ",";
                file_grid_data << MPI_Wtime() << ",";
                file_grid_data << m_uiGlobalMeshElements << ",";
                file_grid_data << m_uiGlobalGridPoints << ",";
                file_grid_data << bssn::BSSN_RK45_TIME_STEP_SIZE << "\n";
                file_grid_data.close();
            } catch (const std::exception& e) {
                std::cout << "Error occured while writing grid summary data!"
                          << std::endl;
                return;
            }
            std::cout << "[BSSNCtx]: Finished writing grid summary data to "
                      << fname << std::endl;
        }
    }

    return;
}

void BSSNCtx::calculate_full_grid_size() {
    if (m_uiMesh->isActive()) {
        // number of mesh elements
        DendroIntL mesh_elements = m_uiMesh->getNumLocalMeshElements();

        DendroIntL grid_points   = m_uiMesh->getNumLocalMeshNodes();

        // perform an all reduce on the mesh
        par::Mpi_Allreduce(&mesh_elements, &m_uiGlobalMeshElements, 1, MPI_SUM,
                           m_uiMesh->getMPICommunicator());

        par::Mpi_Allreduce(&grid_points, &m_uiGlobalGridPoints, 1, MPI_SUM,
                           m_uiMesh->getMPICommunicator());
    }
}

int BSSNCtx::grid_transfer(const ot::Mesh* m_new) {
#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].start();
#endif
    DVec& m_evar = m_var[VL::CPU_EV];
    DVec::grid_transfer(m_uiMesh, m_new, m_evar);
    // printf("igt ended\n");

    m_var[VL::CPU_CV].destroy_vector();
    m_var[VL::CPU_CV_UZ_IN].destroy_vector();

    m_var[VL::CPU_EV_UZ_IN].destroy_vector();
    m_var[VL::CPU_EV_UZ_OUT].destroy_vector();

    m_var[VL::CPU_CV].create_vector(m_new, ot::DVEC_TYPE::OCT_SHARED_NODES,
                                    ot::DVEC_LOC::HOST,
                                    BSSN_CONSTRAINT_NUM_VARS, true);
    m_var[VL::CPU_CV_UZ_IN].create_vector(
        m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_CONSTRAINT_NUM_VARS, true);

    m_var[VL::CPU_EV_UZ_IN].create_vector(
        m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_NUM_VARS, true);
    m_var[VL::CPU_EV_UZ_OUT].create_vector(
        m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING, ot::DVEC_LOC::HOST,
        BSSN_NUM_VARS, true);

    ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS,
                                      BSSN_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(m_new, m_mpi_ctx, BSSN_NUM_VARS,
                                    BSSN_ASYNC_COMM_K);

    m_uiIsETSSynced = false;

#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].stop();
#endif
    return 0;
}

unsigned int BSSNCtx::compute_lts_ts_offset() {
    const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
    const ot::TreeNode* pNodes   = m_uiMesh->getAllElements().data();

    const unsigned int ldiff     = 0;

    unsigned int lmin, lmax;
    m_uiMesh->computeMinMaxLevel(lmin, lmax);
    BSSN_LTS_TS_OFFSET = 0;

    const double dt_min =
        bssn::BSSN_CFL_FACTOR *
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    const double dt_eta_fac = dt_min * ETA_CONST;

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        const unsigned int blev =
            pNodes[blkList[blk].getLocalElementBegin()].getLevel();
        ot::TreeNode blkNode = blkList[blk].getBlockNode();
        const unsigned int szby2 =
            (1u << (m_uiMaxDepth - blkNode.getLevel() - 1));

        Point pt_oct(blkNode.minX() + szby2, blkNode.minY() + szby2,
                     blkNode.minZ() + szby2);
        Point pt_domain;
        m_uiMesh->octCoordToDomainCoord(pt_oct, pt_domain);
        const double r_coord = pt_domain.abs();

        double eta;
        if (bssn::RIT_ETA_FUNCTION == 0) {
            // HAD eta function
            eta = ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
            }
        } else {
            // RIT eta function
            double w   = r_coord / bssn::RIT_ETA_WIDTH;
            double arg = -w * w * w * w;
            eta = (bssn::RIT_ETA_CENTRAL - bssn::RIT_ETA_OUTER) * exp(arg) +
                  bssn::RIT_ETA_OUTER;
        }

        const double dt_eta      = dt_eta_fac * (1 / eta);
        const double dt_cfl      = (1u << (lmax - blev)) * dt_min;

        const double dt_feasible = std::min(dt_eta, dt_cfl);

        if (dt_feasible > dt_min) {
            unsigned int lts_offset =
                lmax - blev - std::floor(std::log2(dt_feasible / dt_min));
            if (BSSN_LTS_TS_OFFSET < lts_offset)
                BSSN_LTS_TS_OFFSET = lts_offset;
        }
    }

    unsigned int lts_offset_max = 0;
    par::Mpi_Allreduce(&BSSN_LTS_TS_OFFSET, &lts_offset_max, 1, MPI_MAX,
                       m_uiMesh->getMPIGlobalCommunicator());
    BSSN_LTS_TS_OFFSET = lts_offset_max;

    if (m_uiMesh->isActive() && (!(m_uiMesh->getMPIRankGlobal())))
        std::cout << "LTS offset : " << BSSN_LTS_TS_OFFSET << std::endl;

    return BSSN_LTS_TS_OFFSET;
}

unsigned int BSSNCtx::getBlkTimestepFac(unsigned int blev, unsigned int lmin,
                                        unsigned int lmax) {
    const unsigned int ldiff = BSSN_LTS_TS_OFFSET;
    if ((lmax - blev) <= ldiff)
        return 1;
    else {
        return 1u << (lmax - blev - ldiff);
    }
}

void BSSNCtx::evolve_bh_loc() {
#ifdef BSSN_EXTRACT_BH_LOCATIONS
    // early exit if we're at time step zero!
    if (m_bBHEvolved) return;

    // if time step is zero, don't evolve, just store and return.
    if (this->m_uiTinfo._m_uiStep == 0) {
        // make sure initial locations are stored, we need the "last" time
        // history for dt calculation
        this->store_bh_loc_history();
        m_bBHEvolved = true;
        return;
    }

    this->compute_constraint_variables();

    // compute how long it's been since the last time we calculatd it, thanks to
    // storing history!
    const double dt = m_uiTinfo._m_uiT - m_uiBHTimeHistory.back();
    DVec sIn        = this->get_evolution_vars();

    // m_uiMesh->readFromGhostBegin(sIn.GetVecArray()+ VAR::U_BETA0 *
    // m_uiMesh->getDegOfFreedom(),3);
    // m_uiMesh->readFromGhostEnd(sIn.GetVecArray()  + VAR::U_BETA0 *
    // m_uiMesh->getDegOfFreedom(),3);
    Point bhLoc[2];
    DendroScalar* evar[bssn::BSSN_NUM_VARS];
    sIn.to_2d(evar);
    bssn::computeBHLocations((const ot::Mesh*)m_uiMesh, m_uiBHLoc, bhLoc, evar,
                             dt);
    // if(!m_uiMesh->getMPIRankGlobal())
    // {
    //     std::cout<<"bh0 "<<bhLoc[0]<<std::endl;
    //     std::cout<<"bh1 "<<bhLoc[1]<<std::endl;

    // }
    m_uiBHLoc[0]         = bhLoc[0];
    m_uiBHLoc[1]         = bhLoc[1];

    bssn::BSSN_BH_LOC[0] = m_uiBHLoc[0];
    bssn::BSSN_BH_LOC[1] = m_uiBHLoc[1];

    // append the bssn location points to our history
    this->store_bh_loc_history();

// old bh location extractor.
#if 0
    DendroScalar *evar[bssn::BSSN_NUM_VARS];
    Point bhLoc[2];
    sIn.to_2d(evar);
    bssn::extractBHCoords((const ot::Mesh *)m_uiMesh,(const DendroScalar*)evar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,(const Point *) m_uiBHLoc,2,(Point*)bhLoc);
    
    m_uiBHLoc[0] = bhLoc[0];
    m_uiBHLoc[1] = bhLoc[1];
#endif
#endif

    // make sure we know that it's been evolved this time step if we're wanting
    // to evolve it
    m_bBHEvolved = true;

    return;
}

void BSSNCtx::store_bh_loc_history() {
    // simple call that stores bh loc history based on bssn::BSSN_BH_LOC
    m_uiBHLocHistory.push_back(
        std::make_pair(bssn::BSSN_BH_LOC[0], bssn::BSSN_BH_LOC[1]));
    m_uiBHTimeHistory.push_back(m_uiTinfo._m_uiT);
}

int BSSNCtx::aeh_expansion(const Point& origin, aeh::AEH_VARS* m_aeh_vars,
                           DVec& aeh_f, DVec& aeh_h,
                           const DendroScalar* const rlim) {
    const unsigned int cg_sz = m_uiMesh->getDegOfFreedom();
    const unsigned int uz_sz = m_uiMesh->getDegOfFreedomUnZip();

    // temporary unzip storage (using evar unzip vectors)
    DendroScalar* unzip_0    = m_var[VL::CPU_EV_UZ_IN].get_vec_ptr();
    DendroScalar* unzip_1    = m_var[VL::CPU_EV_UZ_OUT].get_vec_ptr();

    DendroScalar* F          = aeh_f.get_vec_ptr();
    DendroScalar* F_uz       = &unzip_1[0 * uz_sz];
    DendroScalar* H_uz       = &unzip_1[1 * uz_sz];

    m_uiMesh->readFromGhostBegin<DendroScalar>(F, 1);
    m_uiMesh->readFromGhostEnd<DendroScalar>(F, 1);
    m_uiMesh->unzip(F, F_uz, 1);

    const DendroScalar* const gt_uz  = m_aeh_vars->gt.get_vec_ptr();
    const DendroScalar* const At_uz  = m_aeh_vars->At.get_vec_ptr();
    const DendroScalar* const chi_uz = m_aeh_vars->chi.get_vec_ptr();
    const DendroScalar* const K_uz   = m_aeh_vars->K.get_vec_ptr();

    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW        = bssn::BSSN_PADDING_WIDTH;

    const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

    const DendroScalar r_min     = 1e-2;

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        DendroScalar ptmin[3], ptmax[3];
        unsigned int sz[3];

        const unsigned int offset = blkList[blk].getOffset();
        sz[0]                     = blkList[blk].getAllocationSzX();
        sz[1]                     = blkList[blk].getAllocationSzY();
        sz[2]                     = blkList[blk].getAllocationSzZ();

        const unsigned int bflag  = blkList[blk].getBlkNodeFlag();

        const DendroScalar dx     = blkList[blk].computeDx(pt_min, pt_max);
        const DendroScalar dy     = blkList[blk].computeDy(pt_min, pt_max);
        const DendroScalar dz     = blkList[blk].computeDz(pt_min, pt_max);

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

        if (!((p_rr < (rlim[1] + ll)) && (p_rr > (rlim[0] - ll)))) continue;

        // if(sqrt(ptmax[0] * ptmax[0] + ptmax[1]*ptmax[1] + ptmax[2] *
        // ptmax[2]) > rlim[1])
        //     continue;

        const unsigned int n           = sz[0] * sz[1] * sz[2];
        const unsigned int BLK_SZ      = n;

        const unsigned int nx          = sz[0];
        const unsigned int ny          = sz[1];
        const unsigned int nz          = sz[2];

        const double hx                = (ptmax[0] - ptmin[0]) / (nx - 1);
        const double hy                = (ptmax[1] - ptmin[1]) / (ny - 1);
        const double hz                = (ptmax[2] - ptmin[2]) / (nz - 1);

        DendroScalar* const deriv_base = bssn::BSSN_DERIV_WORKSPACE;

        double* grad_0_F               = deriv_base + 0 * BLK_SZ;
        double* grad_1_F               = deriv_base + 1 * BLK_SZ;
        double* grad_2_F               = deriv_base + 2 * BLK_SZ;

        double* grad2_0_0_F            = deriv_base + 3 * BLK_SZ;
        double* grad2_0_1_F            = deriv_base + 4 * BLK_SZ;
        double* grad2_0_2_F            = deriv_base + 5 * BLK_SZ;

        double* grad2_1_1_F            = deriv_base + 6 * BLK_SZ;
        double* grad2_1_2_F            = deriv_base + 7 * BLK_SZ;
        double* grad2_2_2_F            = deriv_base + 8 * BLK_SZ;

        const double* const F          = &F_uz[offset];
        double* const H                = &H_uz[offset];

        deriv_x(grad_0_F, F, hx, sz, bflag);
        deriv_xx(grad2_0_0_F, F, hx, sz, bflag);

        deriv_y(grad_1_F, F, hy, sz, bflag);
        deriv_yy(grad2_1_1_F, F, hy, sz, bflag);

        deriv_z(grad_2_F, F, hz, sz, bflag);
        deriv_zz(grad2_2_2_F, F, hz, sz, bflag);

        deriv_y(grad2_0_1_F, grad_0_F, hy, sz, bflag);
        deriv_z(grad2_0_2_F, grad_0_F, hz, sz, bflag);
        deriv_z(grad2_1_2_F, grad_1_F, hz, sz, bflag);

        const double* const grad_0_chi =
            &m_aeh_vars->grad_chi.get_vec_ptr()[0 * uz_sz + offset];
        const double* const grad_1_chi =
            &m_aeh_vars->grad_chi.get_vec_ptr()[1 * uz_sz + offset];
        const double* const grad_2_chi =
            &m_aeh_vars->grad_chi.get_vec_ptr()[2 * uz_sz + offset];

        const double* const grad_0_gt0 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[0 * 3 * uz_sz + 0 * uz_sz + offset];
        const double* const grad_1_gt0 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[0 * 3 * uz_sz + 1 * uz_sz + offset];
        const double* const grad_2_gt0 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[0 * 3 * uz_sz + 2 * uz_sz + offset];

        const double* const grad_0_gt1 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[1 * 3 * uz_sz + 0 * uz_sz + offset];
        ;
        const double* const grad_1_gt1 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[1 * 3 * uz_sz + 1 * uz_sz + offset];
        ;
        const double* const grad_2_gt1 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[1 * 3 * uz_sz + 2 * uz_sz + offset];
        ;

        const double* const grad_0_gt2 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[2 * 3 * uz_sz + 0 * uz_sz + offset];
        const double* const grad_1_gt2 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[2 * 3 * uz_sz + 1 * uz_sz + offset];
        const double* const grad_2_gt2 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[2 * 3 * uz_sz + 2 * uz_sz + offset];

        const double* const grad_0_gt3 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[3 * 3 * uz_sz + 0 * uz_sz + offset];
        const double* const grad_1_gt3 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[3 * 3 * uz_sz + 1 * uz_sz + offset];
        const double* const grad_2_gt3 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[3 * 3 * uz_sz + 2 * uz_sz + offset];

        const double* const grad_0_gt4 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[4 * 3 * uz_sz + 0 * uz_sz + offset];
        const double* const grad_1_gt4 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[4 * 3 * uz_sz + 1 * uz_sz + offset];
        const double* const grad_2_gt4 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[4 * 3 * uz_sz + 2 * uz_sz + offset];

        const double* const grad_0_gt5 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[5 * 3 * uz_sz + 0 * uz_sz + offset];
        const double* const grad_1_gt5 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[5 * 3 * uz_sz + 1 * uz_sz + offset];
        const double* const grad_2_gt5 =
            &m_aeh_vars->grad_gt
                 .get_vec_ptr()[5 * 3 * uz_sz + 2 * uz_sz + offset];

        const double* const gt0 = &gt_uz[0 * uz_sz + offset];
        const double* const gt1 = &gt_uz[1 * uz_sz + offset];
        const double* const gt2 = &gt_uz[2 * uz_sz + offset];
        const double* const gt3 = &gt_uz[3 * uz_sz + offset];
        const double* const gt4 = &gt_uz[4 * uz_sz + offset];
        const double* const gt5 = &gt_uz[5 * uz_sz + offset];

        const double* const At0 = &At_uz[0 * uz_sz + offset];
        const double* const At1 = &At_uz[1 * uz_sz + offset];
        const double* const At2 = &At_uz[2 * uz_sz + offset];
        const double* const At3 = &At_uz[3 * uz_sz + offset];
        const double* const At4 = &At_uz[4 * uz_sz + offset];
        const double* const At5 = &At_uz[5 * uz_sz + offset];

        const double* const K   = &K_uz[offset];
        const double* const chi = &chi_uz[offset];

        for (unsigned int k = PW; k < nz - PW; k++) {
            for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
                for (unsigned int i = PW; i < nx - PW; i++) {
                    const unsigned int pp = i + nx * (j + ny * k);
                    const double xx       = ptmin[0] + i * hx - origin.x();
                    const double yy       = ptmin[1] + j * hy - origin.y();
                    const double zz       = ptmin[2] + k * hz - origin.z();

                    const double rr       = sqrt(xx * xx + yy * yy + zz * zz);
                    {
#include "expansion_aeh.cpp"
                        // H[pp] *= sqrt(grad_0_F[pp] * grad_0_F[pp] +
                        // grad_1_F[pp] * grad_1_F[pp] + grad_2_F[pp] *
                        // grad_2_F[pp]);
                    }
                }
            }
        }
    }

    m_uiMesh->zip(H_uz, aeh_h.get_vec_ptr());
    return 0;
}

#if 0
void BSSNCtx::lts_smooth(DVec sIn, LTS_SMOOTH_MODE mode) {
    if (mode == LTS_SMOOTH_MODE::KO) {
        this->unzip(sIn, m_uiEUnzip[0], bssn::BSSN_ASYNC_COMM_K);

        DendroScalar** uZipVars;
        DendroScalar** unzipVarsRHS;

        m_uiEUnzip[0].Get2DArray(uZipVars, false);
        m_uiEUnzip[1].Get2DArray(unzipVarsRHS, false);

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

        unsigned int sz[3];
        const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                           bssn::BSSN_COMPD_MIN[2]);
        const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                           bssn::BSSN_COMPD_MAX[2]);
        double ptmin[3], ptmax[3];
        const unsigned int PW = bssn::BSSN_PADDING_WIDTH;

        mem::memory_pool<double>* __mem_pool = &BSSN_MEM_POOL;

        for (unsigned int blk = 0; blk < numBlocks; blk++) {
            const unsigned int offset = blkList[blk].getOffset();
            sz[0] = blkList[blk].getAllocationSzX();
            sz[1] = blkList[blk].getAllocationSzY();
            sz[2] = blkList[blk].getAllocationSzZ();

            const unsigned int bflag = blkList[blk].getBlkNodeFlag();

            const double dx = blkList[blk].computeDx(pt_min, pt_max);
            const double dy = blkList[blk].computeDy(pt_min, pt_max);
            const double dz = blkList[blk].computeDz(pt_min, pt_max);

            ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
            ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
            ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

            ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
            ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
            ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

            const double* alpha = &uZipVars[VAR::U_ALPHA][offset];
            const double* chi = &uZipVars[VAR::U_CHI][offset];
            const double* K = &uZipVars[VAR::U_K][offset];
            const double* gt0 = &uZipVars[VAR::U_SYMGT0][offset];
            const double* gt1 = &uZipVars[VAR::U_SYMGT1][offset];
            const double* gt2 = &uZipVars[VAR::U_SYMGT2][offset];
            const double* gt3 = &uZipVars[VAR::U_SYMGT3][offset];
            const double* gt4 = &uZipVars[VAR::U_SYMGT4][offset];
            const double* gt5 = &uZipVars[VAR::U_SYMGT5][offset];
            const double* beta0 = &uZipVars[VAR::U_BETA0][offset];
            const double* beta1 = &uZipVars[VAR::U_BETA1][offset];
            const double* beta2 = &uZipVars[VAR::U_BETA2][offset];
            const double* At0 = &uZipVars[VAR::U_SYMAT0][offset];
            const double* At1 = &uZipVars[VAR::U_SYMAT1][offset];
            const double* At2 = &uZipVars[VAR::U_SYMAT2][offset];
            const double* At3 = &uZipVars[VAR::U_SYMAT3][offset];
            const double* At4 = &uZipVars[VAR::U_SYMAT4][offset];
            const double* At5 = &uZipVars[VAR::U_SYMAT5][offset];
            const double* Gt0 = &uZipVars[VAR::U_GT0][offset];
            const double* Gt1 = &uZipVars[VAR::U_GT1][offset];
            const double* Gt2 = &uZipVars[VAR::U_GT2][offset];
            const double* B0 = &uZipVars[VAR::U_B0][offset];
            const double* B1 = &uZipVars[VAR::U_B1][offset];
            const double* B2 = &uZipVars[VAR::U_B2][offset];

            double* a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
            double* chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
            double* K_rhs = &unzipVarsRHS[VAR::U_K][offset];
            double* gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
            double* gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
            double* gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
            double* gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
            double* gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
            double* gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
            double* b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
            double* b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
            double* b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
            double* At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
            double* At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
            double* At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
            double* At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
            double* At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
            double* At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
            double* Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
            double* Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
            double* Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
            double* B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
            double* B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
            double* B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

            const unsigned int nx = sz[0];
            const unsigned int ny = sz[1];
            const unsigned int nz = sz[2];

            double hx = (ptmax[0] - ptmin[0]) / (nx - 1);
            double hy = (ptmax[1] - ptmin[1]) / (ny - 1);
            double hz = (ptmax[2] - ptmin[2]) / (nz - 1);

// clang-format off
//             const unsigned int n = sz[0] * sz[1] * sz[2];
// #include "bssnrhs_lts_smooth_ko.h"

//             const double sigma = KO_DISS_SIGMA;

//             for (unsigned int k = 3; k < nz - 3; k++) {
//                 for (unsigned int j = 3; j < ny - 3; j++) {
//                     for (unsigned int i = 3; i < nx - 3; i++) {
//                         const unsigned int pp = i + nx * (j + ny * k);

//                         a_rhs[pp] = alpha[pp] + sigma * (grad_0_alpha[pp] +
//                                                          grad_1_alpha[pp] +
//                                                          grad_2_alpha[pp]);
//                         b_rhs0[pp] = beta0[pp] + sigma * (grad_0_beta0[pp] +
//                                                           grad_1_beta0[pp] +
//                                                           grad_2_beta0[pp]);
//                         b_rhs1[pp] = beta1[pp] + sigma * (grad_0_beta1[pp] +
//                                                           grad_1_beta1[pp] +
//                                                           grad_2_beta1[pp]);
//                         b_rhs2[pp] = beta2[pp] + sigma * (grad_0_beta2[pp] +
//                                                           grad_1_beta2[pp] +
//                                                           grad_2_beta2[pp]);

//                         gt_rhs00[pp] =
//                             gt0[pp] + sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] +
//                                                grad_2_gt0[pp]);
//                         gt_rhs01[pp] =
//                             gt1[pp] + sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] +
//                                                grad_2_gt1[pp]);
//                         gt_rhs02[pp] =
//                             gt2[pp] + sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] +
//                                                grad_2_gt2[pp]);
//                         gt_rhs11[pp] =
//                             gt3[pp] + sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] +
//                                                grad_2_gt3[pp]);
//                         gt_rhs12[pp] =
//                             gt4[pp] + sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] +
//                                                grad_2_gt4[pp]);
//                         gt_rhs22[pp] =
//                             gt5[pp] + sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] +
//                                                grad_2_gt5[pp]);

//                         chi_rhs[pp] =
//                             chi[pp] + sigma * (grad_0_chi[pp] + grad_1_chi[pp] +
//                                                grad_2_chi[pp]);

//                         At_rhs00[pp] =
//                             At0[pp] + sigma * (grad_0_At0[pp] + grad_1_At0[pp] +
//                                                grad_2_At0[pp]);
//                         At_rhs01[pp] =
//                             At1[pp] + sigma * (grad_0_At1[pp] + grad_1_At1[pp] +
//                                                grad_2_At1[pp]);
//                         At_rhs02[pp] =
//                             At2[pp] + sigma * (grad_0_At2[pp] + grad_1_At2[pp] +
//                                                grad_2_At2[pp]);
//                         At_rhs11[pp] =
//                             At3[pp] + sigma * (grad_0_At3[pp] + grad_1_At3[pp] +
//                                                grad_2_At3[pp]);
//                         At_rhs12[pp] =
//                             At4[pp] + sigma * (grad_0_At4[pp] + grad_1_At4[pp] +
//                                                grad_2_At4[pp]);
//                         At_rhs22[pp] =
//                             At5[pp] + sigma * (grad_0_At5[pp] + grad_1_At5[pp] +
//                                                grad_2_At5[pp]);

//                         K_rhs[pp] =
//                             K[pp] + sigma * (grad_0_K[pp] + grad_1_K[pp] +
//                                              grad_2_K[pp]);

//                         Gt_rhs0[pp] =
//                             Gt0[pp] + sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] +
//                                                grad_2_Gt0[pp]);
//                         Gt_rhs1[pp] =
//                             Gt1[pp] + sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] +
//                                                grad_2_Gt1[pp]);
//                         Gt_rhs2[pp] =
//                             Gt2[pp] + sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] +
//                                                grad_2_Gt2[pp]);

//                         B_rhs0[pp] =
//                             B0[pp] + sigma * (grad_0_B0[pp] + grad_1_B0[pp] +
//                                               grad_2_B0[pp]);
//                         B_rhs1[pp] =
//                             B1[pp] + sigma * (grad_0_B1[pp] + grad_1_B1[pp] +
//                                               grad_2_B1[pp]);
//                         B_rhs2[pp] =
//                             B2[pp] + sigma * (grad_0_B2[pp] + grad_1_B2[pp] +
//                                               grad_2_B2[pp]);
//                     }
//                 }
//             }
//         }
// clang-format on

        this->zip(m_uiEUnzip[1], sIn, bssn::BSSN_ASYNC_COMM_K);
        delete[] uZipVars;
        delete[] unzipVarsRHS;
    }
}
#endif

}  // end of namespace bssn.
