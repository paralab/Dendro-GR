/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief
 */
//

#include "rkfluid.h"

namespace ode {
namespace solver {

RK_FLUID::RK_FLUID(ot::Mesh* pMesh, double pTBegin, double pTEnd, double pTh,
                   RKType rkType)
    : RK(pMesh, pTBegin, pTEnd, pTh) {
    m_uiRKType = rkType;

    if (m_uiRKType == RKType::RK3)
        m_uiNumRKStages = fluid::FLUID_RK3_STAGES;
    else if (m_uiRKType == RKType::RK4)
        m_uiNumRKStages = fluid::FLUID_RK4_STAGES;
    else if (m_uiRKType == RK45)
        m_uiNumRKStages = fluid::FLUID_RK45_STAGES;
    else {
        if (!(pMesh->getMPIRankGlobal()))
            std::cout << "[RK Solver Error]: undefined rk solver type"
                      << std::endl;
        MPI_Abort(m_uiComm, 0);
    }

    // allocate memory for the variables.
    m_uiEvolVar         = new DendroScalar*[fluid::FLUID_NUM_EVOL_VARS];
    m_uiEvolPrevVar     = new DendroScalar*[fluid::FLUID_NUM_EVOL_VARS];
    m_uiEvolVarIm       = new DendroScalar*[fluid::FLUID_NUM_EVOL_VARS];
    m_uiEvolUnzipVar    = new DendroScalar*[fluid::FLUID_NUM_EVOL_VARS];
    m_uiEvolUnzipVarRHS = new DendroScalar*[fluid::FLUID_NUM_EVOL_VARS];

    for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS; index++) {
        m_uiEvolVar[index]     = m_uiMesh->createVector<DendroScalar>();
        m_uiEvolPrevVar[index] = m_uiMesh->createVector<DendroScalar>();
        m_uiEvolVarIm[index]   = m_uiMesh->createVector<DendroScalar>();

        m_uiEvolUnzipVar[index] =
            m_uiMesh->createUnZippedVector<DendroScalar>();
        m_uiEvolUnzipVarRHS[index] =
            m_uiMesh->createUnZippedVector<DendroScalar>();
    }

    m_uiPrimVar      = new DendroScalar*[fluid::FLUID_NUM_PRIM_VARS];
    m_uiPrimPrevVar  = new DendroScalar*[fluid::FLUID_NUM_PRIM_VARS];
    m_uiPrimVarIm    = new DendroScalar*[fluid::FLUID_NUM_PRIM_VARS];
    m_uiPrimUnzipVar = new DendroScalar*[fluid::FLUID_NUM_PRIM_VARS];

    for (unsigned int index = 0; index < fluid::FLUID_NUM_PRIM_VARS; index++) {
        m_uiPrimVar[index]     = m_uiMesh->createVector<DendroScalar>();
        m_uiPrimPrevVar[index] = m_uiMesh->createVector<DendroScalar>();
        m_uiPrimVarIm[index]   = m_uiMesh->createVector<DendroScalar>();

        m_uiPrimUnzipVar[index] =
            m_uiMesh->createUnZippedVector<DendroScalar>();
    }

    m_uiEvolStage = new DendroScalar**[m_uiNumRKStages];

    for (unsigned int stage = 0; stage < m_uiNumRKStages; stage++) {
        m_uiEvolStage[stage] = new DendroScalar*[fluid::FLUID_NUM_EVOL_VARS];

        for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
             index++)
            m_uiEvolStage[stage][index] =
                m_uiMesh->createVector<DendroScalar>();
    }

    // mpi communication
    m_uiSendNodeBuf = new double*[fluid::FLUID_ASYNC_COMM_K];
    m_uiRecvNodeBuf = new double*[fluid::FLUID_ASYNC_COMM_K];

    m_uiSendReqs    = new MPI_Request*[fluid::FLUID_ASYNC_COMM_K];
    m_uiRecvReqs    = new MPI_Request*[fluid::FLUID_ASYNC_COMM_K];
    m_uiSendSts     = new MPI_Status*[fluid::FLUID_ASYNC_COMM_K];
    m_uiRecvSts     = new MPI_Status*[fluid::FLUID_ASYNC_COMM_K];

    for (unsigned int index = 0; index < fluid::FLUID_ASYNC_COMM_K; index++) {
        m_uiSendNodeBuf[index] = NULL;
        m_uiRecvNodeBuf[index] = NULL;

        m_uiSendReqs[index]    = NULL;
        m_uiRecvReqs[index]    = NULL;
        m_uiSendSts[index]     = NULL;
        m_uiRecvSts[index]     = NULL;
    }

    if (m_uiMesh->isActive()) {
        // allocate mpi comm. reqs and status
        for (unsigned int index = 0; index < fluid::FLUID_ASYNC_COMM_K;
             index++) {
            if (m_uiMesh->getGhostExcgTotalSendNodeCount() != 0)
                m_uiSendNodeBuf[index] =
                    new DendroScalar[m_uiMesh
                                         ->getGhostExcgTotalSendNodeCount()];
            if (m_uiMesh->getGhostExcgTotalRecvNodeCount() != 0)
                m_uiRecvNodeBuf[index] =
                    new DendroScalar[m_uiMesh
                                         ->getGhostExcgTotalRecvNodeCount()];

            if (m_uiMesh->getSendProcListSize() != 0) {
                m_uiSendReqs[index] =
                    new MPI_Request[m_uiMesh->getSendProcListSize()];
                m_uiSendSts[index] =
                    new MPI_Status[m_uiMesh->getSendProcListSize()];
            }

            if (m_uiMesh->getRecvProcListSize() != 0) {
                m_uiRecvReqs[index] =
                    new MPI_Request[m_uiMesh->getRecvProcListSize()];
                m_uiRecvSts[index] =
                    new MPI_Status[m_uiMesh->getRecvProcListSize()];
            }
        }
    }
}

RK_FLUID::~RK_FLUID() {
    for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS; index++) {
        delete[] m_uiEvolPrevVar[index];
        delete[] m_uiEvolVar[index];
        delete[] m_uiEvolVarIm[index];
        delete[] m_uiEvolUnzipVarRHS[index];
        delete[] m_uiEvolUnzipVar[index];
    }

    delete[] m_uiEvolPrevVar;
    delete[] m_uiEvolVar;
    delete[] m_uiEvolVarIm;
    delete[] m_uiEvolUnzipVar;
    delete[] m_uiEvolUnzipVarRHS;

    for (unsigned int index = 0; index < fluid::FLUID_NUM_PRIM_VARS; index++) {
        delete[] m_uiPrimPrevVar[index];
        delete[] m_uiPrimVarIm[index];
        delete[] m_uiPrimVar[index];

        delete[] m_uiPrimUnzipVar[index];
    }

    delete[] m_uiPrimPrevVar;
    delete[] m_uiPrimVarIm;
    delete[] m_uiPrimVar;
    delete[] m_uiPrimUnzipVar;

    for (unsigned int stage = 0; stage < m_uiNumRKStages; stage++) {
        for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
             index++)
            delete[] m_uiEvolStage[stage][index];
    }

    for (unsigned int stage = 0; stage < m_uiNumRKStages; stage++) {
        delete[] m_uiEvolStage[stage];
    }

    delete[] m_uiEvolStage;

    // mpi communication
    for (unsigned int index = 0; index < fluid::FLUID_ASYNC_COMM_K; index++) {
        delete[] m_uiSendNodeBuf[index];
        delete[] m_uiRecvNodeBuf[index];

        delete[] m_uiSendReqs[index];
        delete[] m_uiRecvReqs[index];

        delete[] m_uiSendSts[index];
        delete[] m_uiRecvSts[index];
    }

    delete[] m_uiSendNodeBuf;
    delete[] m_uiRecvNodeBuf;

    delete[] m_uiSendReqs;
    delete[] m_uiSendSts;
    delete[] m_uiRecvReqs;
    delete[] m_uiRecvSts;
}

void RK_FLUID::applyInitialConditions(DendroScalar** zipIn,
                                      DendroScalar** zipPrimIn) {
    unsigned int nodeLookUp_CG;
    unsigned int nodeLookUp_DG;
    double x, y, z, len;
    const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin()));
    unsigned int ownerID, ii_x, jj_y, kk_z;
    unsigned int eleOrder      = m_uiMesh->getElementOrder();
    const unsigned int* e2n_cg = &(*(m_uiMesh->getE2NMapping().begin()));
    const unsigned int* e2n_dg = &(*(m_uiMesh->getE2NMapping_DG().begin()));
    const unsigned int nPe     = m_uiMesh->getNumNodesPerElement();
    const unsigned int nodeLocalBegin = m_uiMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd   = m_uiMesh->getNodeLocalEnd();

    double* var                       = new double[fluid::FLUID_NUM_EVOL_VARS];
    double* pvar                      = new double[fluid::FLUID_NUM_PRIM_VARS];

    double mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
    double gu[3][3], gd[3][3], detg, sdetg;

    for (unsigned int elem = m_uiMesh->getElementLocalBegin();
         elem < m_uiMesh->getElementLocalEnd(); elem++) {
        for (unsigned int k = 0; k < (eleOrder + 1); k++)
            for (unsigned int j = 0; j < (eleOrder + 1); j++)
                for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                    nodeLookUp_CG = e2n_cg[elem * nPe +
                                           k * (eleOrder + 1) * (eleOrder + 1) +
                                           j * (eleOrder + 1) + i];
                    if (nodeLookUp_CG >= nodeLocalBegin &&
                        nodeLookUp_CG < nodeLocalEnd) {
                        nodeLookUp_DG =
                            e2n_dg[elem * nPe +
                                   k * (eleOrder + 1) * (eleOrder + 1) +
                                   j * (eleOrder + 1) + i];
                        m_uiMesh->dg2eijk(nodeLookUp_DG, ownerID, ii_x, jj_y,
                                          kk_z);
                        len = (double)(1u << (m_uiMaxDepth -
                                              pNodes[ownerID].getLevel()));
                        x   = pNodes[ownerID].getX() +
                            ii_x * (len / (double)(eleOrder));
                        y = pNodes[ownerID].getY() +
                            jj_y * (len / (double)(eleOrder));
                        z = pNodes[ownerID].getZ() +
                            kk_z * (len / (double)(eleOrder));
                        assert(nodeLookUp_CG < m_uiMesh->getDegOfFreedom());
                        fluid::initData((DendroScalar)x, (DendroScalar)y,
                                        (DendroScalar)z, pvar);
                        // Calculate the metric and apply prim to con.
                        const double pos[] = {GRIDX_TO_X((DendroScalar)x),
                                              GRIDY_TO_Y((DendroScalar)y),
                                              GRIDZ_TO_Z((DendroScalar)z)};
                        metric_vars(gu, &detg, gd, pos);
                        sdetg = sqrt(detg);
                        fluid_math::prim_to_con_pt(pvar, var, gd, sdetg);

                        for (unsigned int v = 0; v < fluid::FLUID_NUM_EVOL_VARS;
                             v++) {
                            zipIn[v][nodeLookUp_CG] = var[v];
                        }

                        for (unsigned int v = 0; v < fluid::FLUID_NUM_PRIM_VARS;
                             v++) {
                            zipPrimIn[v][nodeLookUp_CG] = pvar[v];
                        }
                    }
                }
    }

    delete[] pvar;
    delete[] var;
}

void RK_FLUID::initialGridConverge() {
    applyInitialConditions(m_uiEvolPrevVar, m_uiPrimPrevVar);

    bool isRefine = false;
    unsigned int oldElements, oldElements_g;
    unsigned int newElements, newElements_g;

    // refine based on all the variables
    const unsigned int refineNumVars = fluid::FLUID_NUM_REFINE_VARS;
    unsigned int refineVarIds[refineNumVars];
    for (unsigned int vIndex = 0; vIndex < refineNumVars; vIndex++)
        refineVarIds[vIndex] = fluid::FLUID_REFINE_VARIABLE_INDICES[vIndex];

    double wTol = fluid::FLUID_WAVELET_TOL;
    std::function<double(double, double, double)> waveletTolFunc =
        [wTol](double x, double y, double z) {
            return fluid::computeWTol(x, y, z, wTol);
        };
    unsigned int iterCount = 1;
    do {
        unzipVars_async(m_uiEvolPrevVar, m_uiEvolUnzipVar,
                        fluid::FLUID_NUM_EVOL_VARS);
        unzipPrimVars(m_uiPrimPrevVar, m_uiPrimUnzipVar);
        // unzipVars_async(m_uiPrimPrevVar,m_uiPrimUnzipVar,fluid::FLUID_NUM_PRIM_VARS);

        // enforce boundary conditions
        fluid::fluid_bcs_cons(m_uiMesh, m_uiEvolUnzipVar);
        fluid::fluid_bcs_prim(m_uiMesh, m_uiPrimUnzipVar);

        std::vector<unsigned int> refine_flags;
        if (fluid::FLUID_ENABLE_BLOCK_ADAPTIVITY)
            isRefine = false;
        else {
            // isRefine=m_uiMesh->isReMeshUnzip((const double
            // **)m_uiEvolUnzipVar,refineVarIds,refineNumVars,waveletTolFunc,fluid::FLUID_DENDRO_AMR_FAC);
            isRefine = fluid::computeFluidRemeshFlags(
                m_uiMesh, refine_flags, (const double**)m_uiPrimUnzipVar,
                refineVarIds, refineNumVars, waveletTolFunc,
                0.0 * fluid::FLUID_DENDRO_AMR_FAC);
            m_uiMesh->setMeshRefinementFlags(refine_flags);
        }

        if (isRefine) {
            ot::Mesh* newMesh = m_uiMesh->ReMesh(fluid::FLUID_DENDRO_GRAIN_SZ,
                                                 fluid::FLUID_LOAD_IMB_TOL,
                                                 fluid::FLUID_SPLIT_FIX);

            oldElements       = m_uiMesh->getNumLocalMeshElements();
            newElements       = newMesh->getNumLocalMeshElements();

            par::Mpi_Allreduce(&oldElements, &oldElements_g, 1, MPI_SUM,
                               m_uiMesh->getMPIGlobalCommunicator());
            par::Mpi_Allreduce(&newElements, &newElements_g, 1, MPI_SUM,
                               m_uiMesh->getMPIGlobalCommunicator());

            if (!(m_uiMesh->getMPIRankGlobal()))
                std::cout << "initial grid iteration : " << iterCount
                          << " old mesh: " << oldElements_g
                          << " new mesh: " << newElements_g << std::endl;

            // performs the inter-grid transfer
            intergridTransferVars(m_uiEvolPrevVar, newMesh,
                                  fluid::FLUID_NUM_EVOL_VARS);
            intergridTransferPrimitiveVars(m_uiPrimPrevVar, newMesh);
            fluid_math::prim_to_con(newMesh, m_uiPrimPrevVar, m_uiEvolPrevVar);

            for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                 index++) {
                delete[] m_uiEvolVar[index];
                delete[] m_uiEvolVarIm[index];
                delete[] m_uiEvolUnzipVar[index];
                delete[] m_uiEvolUnzipVarRHS[index];

                m_uiEvolVar[index]         = NULL;
                m_uiEvolVarIm[index]       = NULL;
                m_uiEvolUnzipVar[index]    = NULL;
                m_uiEvolUnzipVarRHS[index] = NULL;

                m_uiEvolVar[index]   = newMesh->createVector<DendroScalar>();
                m_uiEvolVarIm[index] = newMesh->createVector<DendroScalar>();

                m_uiEvolUnzipVar[index] =
                    newMesh->createUnZippedVector<DendroScalar>();
                m_uiEvolUnzipVarRHS[index] =
                    newMesh->createUnZippedVector<DendroScalar>();
            }

            for (unsigned int index = 0; index < fluid::FLUID_NUM_PRIM_VARS;
                 index++) {
                delete[] m_uiPrimVar[index];
                delete[] m_uiPrimVarIm[index];
                delete[] m_uiPrimUnzipVar[index];

                m_uiPrimVar[index]      = NULL;
                m_uiPrimVarIm[index]    = NULL;
                m_uiPrimUnzipVar[index] = NULL;

                m_uiPrimVar[index]      = newMesh->createVector<DendroScalar>();
                m_uiPrimVarIm[index]    = newMesh->createVector<DendroScalar>();

                m_uiPrimUnzipVar[index] =
                    newMesh->createUnZippedVector<DendroScalar>();
            }

            for (unsigned int stage = 0; stage < m_uiNumRKStages; stage++) {
                for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                     index++) {
                    delete[] m_uiEvolStage[stage][index];
                    m_uiEvolStage[stage][index] = NULL;
                    m_uiEvolStage[stage][index] =
                        newMesh->createVector<DendroScalar>();
                }
            }

            std::swap(newMesh, m_uiMesh);
            delete newMesh;

            reallocateMPIResources();

            applyInitialConditions(m_uiEvolPrevVar, m_uiPrimPrevVar);

            if (m_uiMesh->isActive()) {
                double l_min = vecMin(m_uiEvolPrevVar[fluid::PVAR::V_RHO] +
                                          m_uiMesh->getNodeLocalBegin(),
                                      (m_uiMesh->getNumLocalMeshNodes()),
                                      m_uiMesh->getMPICommunicator());
                double l_max = vecMax(m_uiEvolPrevVar[fluid::PVAR::V_RHO] +
                                          m_uiMesh->getNodeLocalBegin(),
                                      (m_uiMesh->getNumLocalMeshNodes()),
                                      m_uiMesh->getMPICommunicator());
                if (!(m_uiMesh->getMPIRank())) {
                    std::cout << "transfer completed:    ||PVAR::V_RHO|| (min, "
                                 "max) : ("
                              << l_min << ", " << l_max << " ) " << std::endl;
                }
            }

            iterCount += 1;
        }

    } while (isRefine && (newElements_g != oldElements_g));
}

void RK_FLUID::reallocateMPIResources() {
    for (unsigned int index = 0; index < fluid::FLUID_ASYNC_COMM_K; index++) {
        delete[] m_uiSendNodeBuf[index];
        delete[] m_uiRecvNodeBuf[index];

        delete[] m_uiSendReqs[index];
        delete[] m_uiRecvReqs[index];

        delete[] m_uiSendSts[index];
        delete[] m_uiRecvSts[index];
    }

    for (unsigned int index = 0; index < fluid::FLUID_ASYNC_COMM_K; index++) {
        m_uiSendNodeBuf[index] = NULL;
        m_uiRecvNodeBuf[index] = NULL;

        m_uiSendReqs[index]    = NULL;
        m_uiRecvReqs[index]    = NULL;
        m_uiSendSts[index]     = NULL;
        m_uiRecvSts[index]     = NULL;
    }

    if (m_uiMesh->isActive()) {
        // allocate mpi comm. reqs and status
        for (unsigned int index = 0; index < fluid::FLUID_ASYNC_COMM_K;
             index++) {
            if (m_uiMesh->getGhostExcgTotalSendNodeCount() != 0)
                m_uiSendNodeBuf[index] =
                    new DendroScalar[m_uiMesh
                                         ->getGhostExcgTotalSendNodeCount()];
            if (m_uiMesh->getGhostExcgTotalRecvNodeCount() != 0)
                m_uiRecvNodeBuf[index] =
                    new DendroScalar[m_uiMesh
                                         ->getGhostExcgTotalRecvNodeCount()];

            if (m_uiMesh->getSendProcListSize() != 0) {
                m_uiSendReqs[index] =
                    new MPI_Request[m_uiMesh->getSendProcListSize()];
                m_uiSendSts[index] =
                    new MPI_Status[m_uiMesh->getSendProcListSize()];
            }

            if (m_uiMesh->getRecvProcListSize() != 0) {
                m_uiRecvReqs[index] =
                    new MPI_Request[m_uiMesh->getRecvProcListSize()];
                m_uiRecvSts[index] =
                    new MPI_Status[m_uiMesh->getRecvProcListSize()];
            }
        }
    }
}

void RK_FLUID::writeToVTU(DendroScalar** evolZipVarIn, DendroScalar** primVars,
                          DendroScalar** constrZipVarIn,
                          unsigned int numEvolVars, unsigned int numPrimVars,
                          unsigned int numConstVars,
                          const unsigned int* evolVarIndices,
                          const unsigned int* primVarIndices,
                          const unsigned int* constVarIndices, bool zslice) {
    fluid::timer::t_ioVtu.start();

    std::vector<std::string> pDataNames;
    double* pData[(numEvolVars + numPrimVars + numConstVars)];

    for (unsigned int i = 0; i < numEvolVars; i++) {
        pDataNames.push_back(
            std::string(fluid::FLUID_EVOL_VAR_NAMES[evolVarIndices[i]]));
        pData[i] = evolZipVarIn[evolVarIndices[i]];
    }

    for (unsigned int i = 0; i < numPrimVars; i++) {
        pDataNames.push_back(
            std::string(fluid::FLUID_PRIM_VAR_NAMES[primVarIndices[i]]));
        pData[numEvolVars + i] = primVars[primVarIndices[i]];
    }

    for (unsigned int i = 0; i < numConstVars; i++) {
        pDataNames.push_back(
            std::string(fluid::FLUID_CONS_VAR_NAMES[constVarIndices[i]]));
        pData[numEvolVars + numPrimVars + i] =
            constrZipVarIn[constVarIndices[i]];
    }

    std::vector<char*> pDataNames_char;
    pDataNames_char.reserve(pDataNames.size());

    for (unsigned int i = 0; i < pDataNames.size(); i++)
        pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

    const char* fDataNames[] = {"Time", "Cycle"};
    const double fData[]     = {m_uiCurrentTime, (double)m_uiCurrentStep};

    char fPrefix[256];
    sprintf(fPrefix, "%s_%d", fluid::FLUID_VTU_FILE_PREFIX.c_str(),
            m_uiCurrentStep);
    if (!zslice) {
        io::vtk::mesh2vtuFine(m_uiMesh, fPrefix, 2, fDataNames, fData,
                              (numEvolVars + numPrimVars + numConstVars),
                              (const char**)&pDataNames_char[0],
                              (const double**)pData);
    } else {
        unsigned int s_val[3]  = {1u << (m_uiMaxDepth - 1),
                                  1u << (m_uiMaxDepth - 1),
                                  1u << (m_uiMaxDepth - 1)};
        unsigned int s_norm[3] = {0, 0, 1};

        io::vtk::mesh2vtu_slice(
            m_uiMesh, s_val, s_norm, fPrefix, 2, fDataNames, fData,
            (numEvolVars + numPrimVars + numConstVars),
            (const char**)&pDataNames_char[0], (const double**)pData);
    }

    fluid::timer::t_ioVtu.stop();
}

void RK_FLUID::performGhostExchangeVars(DendroScalar** zipIn,
                                        unsigned int numVars) {
    fluid::timer::t_ghostEx_sync.start();

    for (unsigned int v = 0; v < numVars; v++)
        m_uiMesh->performGhostExchange(zipIn[v]);

    fluid::timer::t_ghostEx_sync.stop();
}

void RK_FLUID::intergridTransferVars(DendroScalar**& zipIn,
                                     const ot::Mesh* pnewMesh,
                                     unsigned int numVars) {
    fluid::timer::t_gridTransfer.start();

    for (unsigned int v = 0; v < numVars; v++)
        m_uiMesh->interGridTransfer(zipIn[v], pnewMesh);

    fluid::timer::t_gridTransfer.stop();
}

void RK_FLUID::intergridTransferPrimitiveVars(DendroScalar**& zipPrimIn,
                                              ot::Mesh* pnewMesh) {
    fluid::timer::t_gridTransfer.start();

    m_uiMesh->interGridTransfer(zipPrimIn[fluid::PVAR::V_RHO], pnewMesh);
    m_uiMesh->interGridTransfer(zipPrimIn[fluid::PVAR::V_P], pnewMesh);
    m_uiMesh->interGridTransfer(zipPrimIn[fluid::PVAR::V_U1], pnewMesh);
    m_uiMesh->interGridTransfer(zipPrimIn[fluid::PVAR::V_U2], pnewMesh);
    m_uiMesh->interGridTransfer(zipPrimIn[fluid::PVAR::V_U3], pnewMesh);

    // release mem for previous 3 velocity vector
    delete[] zipPrimIn[fluid::PVAR::V_VX];
    delete[] zipPrimIn[fluid::PVAR::V_VY];
    delete[] zipPrimIn[fluid::PVAR::V_VZ];

    zipPrimIn[fluid::PVAR::V_VX]      = pnewMesh->createVector<DendroScalar>();
    zipPrimIn[fluid::PVAR::V_VY]      = pnewMesh->createVector<DendroScalar>();
    zipPrimIn[fluid::PVAR::V_VZ]      = pnewMesh->createVector<DendroScalar>();

    const unsigned int nodeLocalBegin = pnewMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd   = pnewMesh->getNodeLocalEnd();

    double v4[3];
    double v3[3];

    double gd[3][3], gu[3][3], detg, sdetg;
    double pos[3];
    const unsigned int* cg_to_dg = &(*(pnewMesh->getCG2DGMap().begin()));
    unsigned int lookUp;
    unsigned int ownerID, ii_x, jj_y, kk_z, sz;
    const ot::TreeNode* allElements = &(*(pnewMesh->getAllElements().begin()));
    const unsigned int elementOrder = pnewMesh->getElementOrder();
    double alpha                    = 1.0;
    double Beta[3]                  = {0.0, 0.0, 0.0};

    pnewMesh->performGhostExchange(zipPrimIn[fluid::PVAR::V_U1]);
    pnewMesh->performGhostExchange(zipPrimIn[fluid::PVAR::V_U2]);
    pnewMesh->performGhostExchange(zipPrimIn[fluid::PVAR::V_U3]);

    for (unsigned int node = nodeLocalBegin; node < nodeLocalEnd; node++) {
        lookUp = cg_to_dg[node];
        pnewMesh->dg2eijk(lookUp, ownerID, ii_x, jj_y, kk_z);
        sz = (1u << (m_uiMaxDepth - allElements[ownerID].getLevel()));

        pos[0] =
            allElements[ownerID].minX() + ii_x * (sz / ((double)elementOrder));
        pos[1] =
            allElements[ownerID].minY() + jj_y * (sz / ((double)elementOrder));
        pos[2] =
            allElements[ownerID].minZ() + kk_z * (sz / ((double)elementOrder));

        metric_vars(gu, &detg, gd, pos);
        sdetg = sqrt(detg);

        v4[0] = zipPrimIn[fluid::PVAR::V_U1][node];
        v4[1] = zipPrimIn[fluid::PVAR::V_U2][node];
        v4[2] = zipPrimIn[fluid::PVAR::V_U3][node];

        fluid_math::Vvec4_to_Vvec3(v4, v3, gd, alpha, Beta);

        zipPrimIn[fluid::PVAR::V_VX][node] = v3[0];
        zipPrimIn[fluid::PVAR::V_VY][node] = v3[1];
        zipPrimIn[fluid::PVAR::V_VZ][node] = v3[2];
    }

    fluid::timer::t_gridTransfer.stop();
}

void RK_FLUID::unzipVars(DendroScalar** zipIn, DendroScalar** uzipOut,
                         unsigned int numVars) {
    fluid::timer::t_unzip_sync.start();

    for (unsigned int index = 0; index < numVars; index++)
        m_uiMesh->unzip(zipIn[index], uzipOut[index]);

    fluid::timer::t_unzip_sync.stop();
}

void RK_FLUID::unzipVars_async(DendroScalar** zipIn, DendroScalar** uzipOut,
                               unsigned int numVars) {
    fluid::timer::t_unzip_async.start();

    // for(unsigned int var=0;var<numVars;var+=fluid::FLUID_ASYNC_COMM_K){

    //     for(unsigned int i=0;i<fluid::FLUID_ASYNC_COMM_K;i++)
    //     {
    //         if((var + i)< numVars)
    //             m_uiMesh->ghostExchangeStart(zipIn[var+i],m_uiSendNodeBuf[i],m_uiRecvNodeBuf[i],m_uiSendReqs[i],m_uiRecvReqs[i]);
    //     }

    //     for(unsigned int i=0;i<fluid::FLUID_ASYNC_COMM_K;i++)
    //     {
    //         if((var + i)< numVars)
    //         {
    //             m_uiMesh->ghostExchangeRecvSync(zipIn[var + i],
    //             m_uiRecvNodeBuf[i],m_uiRecvReqs[i], m_uiRecvSts[i]);
    //             m_uiMesh->unzip(zipIn[var+i],uzipOut[var+i]);
    //         }

    //     }

    //     for(unsigned int i=0;i<fluid::FLUID_ASYNC_COMM_K;i++)
    //     {
    //         if((var + i)< numVars)
    //             m_uiMesh->ghostExchangeSendSync(m_uiSendReqs[i],
    //             m_uiSendSts[i]);
    //     }

    // }

    for (unsigned int v = 0; v < numVars; v++) {
        m_uiMesh->performGhostExchange(zipIn[v]);
        m_uiMesh->unzip(zipIn[v], uzipOut[v], 1);
    }

    /*const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();
    unsigned int sz[3],offset,bflag;
    for(unsigned int blk=0; blk<blkList.size(); blk++)
    {

        offset=blkList[blk].getOffset();
        sz[0]=blkList[blk].getAllocationSzX();
        sz[1]=blkList[blk].getAllocationSzY();
        sz[2]=blkList[blk].getAllocationSzZ();

        bflag=blkList[blk].getBlkNodeFlag();

        const unsigned int nx=sz[0];
        const unsigned int ny=sz[1];
        const unsigned int nz=sz[2];

        for(unsigned int var =0;var<numVars;var++)
        {
            for(unsigned int k=0;k<nz;k++)
                for(unsigned int j=0;j<ny;j++)
                    uzipOut[var][offset+ IDX(0,j,k)]=0.0;

            for(unsigned int k=0;k<nz;k++)
                for(unsigned int j=0;j<ny;j++)
                    uzipOut[var][offset+ IDX(nx-1,j,k)]=0.0;


            for(unsigned int k=0;k<nz;k++)
                for(unsigned int i=0;i<nx;i++)
                    uzipOut[var][offset+ IDX(i,0,k)]=0.0;

            for(unsigned int k=0;k<nz;k++)
                for(unsigned int i=0;i<nx;i++)
                    uzipOut[var][offset+ IDX(i,ny-1,k)]=0.0;

            for(unsigned int j=0;j<ny;j++)
                for(unsigned int i=0;i<nx;i++)
                    uzipOut[var][offset+ IDX(i,j,0)]=0.0;

            for(unsigned int j=0;j<ny;j++)
                for(unsigned int i=0;i<nx;i++)
                    uzipOut[var][offset+ IDX(i,j,nz-1)]=0.0;



        }





    }*/

    fluid::timer::t_unzip_async.stop();
}

void RK_FLUID::unzipPrimVars(DendroScalar** zipIn, DendroScalar** unZipOut) {
    double* zipPrim[]   = {zipIn[fluid::PVAR::V_RHO], zipIn[fluid::PVAR::V_U1],
                           zipIn[fluid::PVAR::V_U2], zipIn[fluid::PVAR::V_U3],
                           zipIn[fluid::PVAR::V_P]};
    double* unzipPrim[] = {
        unZipOut[fluid::PVAR::V_RHO], unZipOut[fluid::PVAR::V_U1],
        unZipOut[fluid::PVAR::V_U2], unZipOut[fluid::PVAR::V_U3],
        unZipOut[fluid::PVAR::V_P]};

    // unzipVars_async(zipPrim,unzipPrim,5);
    unzipVars_async(zipIn, unZipOut, 8);

    // compute  unzip 3 vector from 4 vector.
    unsigned int offset, bflag;
    unsigned int sz[3];
    DendroScalar dx, dy, dz;
    DendroScalar ptmin[3];
    DendroScalar ptmax[3];

    const Point pt_min(fluid::FLUID_COMPD_MIN[0], fluid::FLUID_COMPD_MIN[1],
                       fluid::FLUID_COMPD_MIN[2]);
    const Point pt_max(fluid::FLUID_COMPD_MAX[0], fluid::FLUID_COMPD_MAX[1],
                       fluid::FLUID_COMPD_MAX[2]);

    const std::vector<ot::Block>& blkList = m_uiMesh->getLocalBlockList();
    double gd[3][3], gu[3][3], detg, sdetg, pos[3];

    double v4[3];
    double v3[3];
    double alpha   = 1.0;
    double Beta[3] = {0.0, 0.0, 0.0};

    for (unsigned int blk = 0; blk < blkList.size(); blk++) {
        offset   = blkList[blk].getOffset();
        sz[0]    = blkList[blk].getAllocationSzX();
        sz[1]    = blkList[blk].getAllocationSzY();
        sz[2]    = blkList[blk].getAllocationSzZ();

        bflag    = blkList[blk].getBlkNodeFlag();

        dx       = blkList[blk].computeDx(pt_min, pt_max);
        dy       = blkList[blk].computeDy(pt_min, pt_max);
        dz       = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - 3 * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - 3 * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - 3 * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + 3 * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + 3 * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + 3 * dz;

        for (unsigned int k = 0; k < sz[2]; k++) {
            pos[2] = ptmin[2] + k * dz;
            for (unsigned int j = 0; j < sz[1]; j++) {
                pos[1] = ptmin[1] + j * dy;
                for (unsigned int i = 0; i < sz[0]; i++) {
                    pos[0] = ptmin[0] + i * dx;
                    metric_vars(gu, &detg, gd, pos);
                    const unsigned int pp = i + sz[0] * (j + k * sz[1]);

                    v4[0] = unZipOut[fluid::PVAR::V_U1][offset + pp];
                    v4[1] = unZipOut[fluid::PVAR::V_U2][offset + pp];
                    v4[2] = unZipOut[fluid::PVAR::V_U3][offset + pp];

                    fluid_math::Vvec4_to_Vvec3(v4, v3, gd, alpha, Beta);

                    unZipOut[fluid::PVAR::V_VX][offset + pp] = v3[0];
                    unZipOut[fluid::PVAR::V_VY][offset + pp] = v3[1];
                    unZipOut[fluid::PVAR::V_VZ][offset + pp] = v3[2];

                    // Floor the density and pressure.
                    if (unZipOut[fluid::PVAR::V_RHO][offset + pp] <
                        fluid::FLUID_VACUUM_RESET) {
                        unZipOut[fluid::PVAR::V_RHO][offset + pp] =
                            fluid::FLUID_VACUUM_RESET;
                        unZipOut[fluid::PVAR::V_VX][offset + pp] = 0.0;
                        unZipOut[fluid::PVAR::V_VY][offset + pp] = 0.0;
                        unZipOut[fluid::PVAR::V_VZ][offset + pp] = 0.0;
                    }
                    unZipOut[fluid::PVAR::V_P][offset + pp] =
                        fmax(unZipOut[fluid::PVAR::V_P][offset + pp],
                             fluid::FLUID_VACUUM_TAU_RESET);
                }
            }
        }
    }
}

void RK_FLUID::zipVars(DendroScalar** uzipIn, DendroScalar** zipOut,
                       unsigned int numVars) {
    fluid::timer::t_zip.start();

    for (unsigned int index = 0; index < numVars; index++)
        m_uiMesh->zip(uzipIn[index], zipOut[index]);

    fluid::timer::t_zip.stop();
}

void RK_FLUID::applyBoundaryConditions() {}

void RK_FLUID::performSingleIteration() {
    if (m_uiMesh->isActive()) {
        double current_t     = m_uiCurrentTime;
        double current_t_adv = current_t;

        // unzip both primitive and evolve variables.
        unzipVars_async(m_uiEvolPrevVar, m_uiEvolUnzipVar,
                        fluid::FLUID_NUM_EVOL_VARS);
        unzipPrimVars(m_uiPrimPrevVar, m_uiPrimUnzipVar);
        // FIXME: Do we need to call prim_to_con here?
        // unzipVars_async(m_uiPrimPrevVar,m_uiPrimUnzipVar,fluid::FLUID_NUM_PRIM_VARS);

        // enforce boundary conditions
        fluid::fluid_bcs_cons(m_uiMesh, m_uiEvolUnzipVar);
        fluid::fluid_bcs_prim(
            m_uiMesh, m_uiPrimUnzipVar);  // FIXME: Possibly unnecessary?

        int rank                              = m_uiMesh->getMPIRank();

        const unsigned int nodeLocalBegin     = m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd       = m_uiMesh->getNodeLocalEnd();

        const std::vector<ot::Block>& blkList = m_uiMesh->getLocalBlockList();
        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx, dy, dz;
        const Point pt_min(fluid::FLUID_COMPD_MIN[0], fluid::FLUID_COMPD_MIN[1],
                           fluid::FLUID_COMPD_MIN[2]);
        const Point pt_max(fluid::FLUID_COMPD_MAX[0], fluid::FLUID_COMPD_MAX[1],
                           fluid::FLUID_COMPD_MAX[2]);
        std::set<unsigned int> badBounds;

        // use prev prim var as a current prim var guess.
        for (unsigned int var = 0; var < fluid::FLUID_NUM_PRIM_VARS; var++)
            for (unsigned int node = nodeLocalBegin; node < nodeLocalEnd;
                 node++)
                m_uiPrimVar[var][node] = m_uiPrimPrevVar[var][node];

        if (m_uiRKType == RKType::RK3) {
#ifdef FLUID_RK_DEBUG_NAN
            for (unsigned int var = 0; var < fluid::FLUID_NUM_EVOL_VARS; var++)
                ot::test::isUnzipNaN(m_uiMesh, m_uiEvolUnzipVar[var]);

            for (unsigned int var = 0; var < fluid::FLUID_NUM_PRIM_VARS; var++)
                ot::test::isUnzipNaN(m_uiMesh, m_uiPrimUnzipVar[var]);
#endif

            // compute fluid rhs (evol vars)
            fluid::fluid_rhs(m_uiMesh, current_t_adv, m_uiEvolUnzipVarRHS,
                             (const DendroScalar**)m_uiEvolUnzipVar,
                             (const DendroScalar**)m_uiPrimUnzipVar);
            zipVars(m_uiEvolUnzipVarRHS, m_uiEvolStage[0],
                    fluid::FLUID_NUM_EVOL_VARS);

#ifdef FLUID_RK_DEBUG_NAN
            for (unsigned int var = 0; var < fluid::FLUID_NUM_EVOL_VARS; var++)
                ot::test::isZipNAN(m_uiMesh, m_uiEvolStage[0][var]);
#endif

            for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                 index++) {
                for (unsigned int node = nodeLocalBegin; node < nodeLocalEnd;
                     node++) {
                    m_uiEvolStage[0][index][node] =
                        m_uiEvolPrevVar[index][node] +
                        m_uiT_h * m_uiEvolStage[0][index][node];
                }
            }

            // FIXME: Application of boundary conditions should come before
            // con_to_prim and grid sync.
            // fluid::fluid_bcs_cons(m_uiMesh,m_uiEvolStage);

            fluid_math::con_to_prim(m_uiMesh, m_uiEvolStage[0], m_uiPrimVar,
                                    badBounds);

            // unzip both primitive and evolve variables.
            unzipVars_async(m_uiEvolStage[0], m_uiEvolUnzipVar,
                            fluid::FLUID_NUM_EVOL_VARS);
            unzipPrimVars(m_uiPrimVar, m_uiPrimUnzipVar);

            // enforce boundary conditions
            fluid::fluid_bcs_cons(m_uiMesh, m_uiEvolUnzipVar);
            fluid::fluid_bcs_prim(m_uiMesh, m_uiPrimUnzipVar);

#ifdef FLUID_RK_DEBUG_NAN
            for (unsigned int var = 0; var < fluid::FLUID_NUM_EVOL_VARS; var++)
                ot::test::isUnzipNaN(m_uiMesh, m_uiEvolUnzipVar[var]);

            for (unsigned int var = 0; var < fluid::FLUID_NUM_PRIM_VARS; var++)
                ot::test::isUnzipNaN(m_uiMesh, m_uiPrimUnzipVar[var]);
#endif

            current_t_adv += m_uiT_h;

            fluid::fluid_rhs(m_uiMesh, current_t_adv, m_uiEvolUnzipVarRHS,
                             (const DendroScalar**)m_uiEvolUnzipVar,
                             (const DendroScalar**)m_uiPrimUnzipVar);
            zipVars(m_uiEvolUnzipVarRHS, m_uiEvolStage[1],
                    fluid::FLUID_NUM_EVOL_VARS);

            for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                 index++) {
                for (unsigned int node = nodeLocalBegin; node < nodeLocalEnd;
                     node++) {
                    m_uiEvolStage[1][index][node] =
                        (0.75) * m_uiEvolPrevVar[index][node] +
                        0.25 * m_uiEvolStage[0][index][node] +
                        m_uiT_h * 0.25 * m_uiEvolStage[1][index][node];
                }
            }

            // FIXME: Application of boundary conditions should come before
            // con_to_prim and grid sync.

            fluid_math::con_to_prim(m_uiMesh, m_uiEvolStage[1], m_uiPrimVar,
                                    badBounds);
            // unzip both primitive and evolve variables.
            unzipVars_async(m_uiEvolStage[1], m_uiEvolUnzipVar,
                            fluid::FLUID_NUM_EVOL_VARS);
            unzipPrimVars(m_uiPrimVar, m_uiPrimUnzipVar);

            // enforce boundary conditions
            fluid::fluid_bcs_cons(m_uiMesh, m_uiEvolUnzipVar);
            fluid::fluid_bcs_prim(m_uiMesh, m_uiPrimUnzipVar);

#ifdef FLUID_RK_DEBUG_NAN
            for (unsigned int var = 0; var < fluid::FLUID_NUM_EVOL_VARS; var++)
                ot::test::isUnzipNaN(m_uiMesh, m_uiEvolUnzipVar[var]);

            for (unsigned int var = 0; var < fluid::FLUID_NUM_PRIM_VARS; var++)
                ot::test::isUnzipNaN(m_uiMesh, m_uiPrimUnzipVar[var]);
#endif

            current_t_adv += m_uiT_h;
            fluid::fluid_rhs(m_uiMesh, current_t_adv, m_uiEvolUnzipVarRHS,
                             (const DendroScalar**)m_uiEvolUnzipVar,
                             (const DendroScalar**)m_uiPrimUnzipVar);
            zipVars(m_uiEvolUnzipVarRHS, m_uiEvolVar,
                    fluid::FLUID_NUM_EVOL_VARS);

            for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                 index++) {
                for (unsigned int node = nodeLocalBegin; node < nodeLocalEnd;
                     node++)
                    m_uiEvolVar[index][node] =
                        (1.0 / 3.0) * m_uiEvolPrevVar[index][node] +
                        (2.0 / 3.0) * m_uiEvolStage[1][index][node] +
                        m_uiT_h * (2.0 / 3.0) * m_uiEvolVar[index][node];
            }

            fluid_math::con_to_prim(m_uiMesh, m_uiEvolVar, m_uiPrimVar,
                                    badBounds);

        } else if (m_uiRKType == RKType::RK4) {
            for (unsigned int stage = 0; stage < (fluid::FLUID_RK4_STAGES - 1);
                 stage++) {
                // compute fluid rhs (evol vars)
                fluid::fluid_rhs(m_uiMesh, current_t_adv, m_uiEvolUnzipVarRHS,
                                 (const DendroScalar**)m_uiEvolUnzipVar,
                                 (const DendroScalar**)m_uiPrimUnzipVar);
                zipVars(m_uiEvolUnzipVarRHS, m_uiEvolStage[stage],
                        fluid::FLUID_NUM_EVOL_VARS);

                for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                     index++) {
                    for (unsigned int node = nodeLocalBegin;
                         node < nodeLocalEnd; node++) {
                        m_uiEvolVarIm[index][node] =
                            m_uiEvolPrevVar[index][node];
                        m_uiEvolVarIm[index][node] +=
                            (RK4_U[stage + 1] * m_uiT_h *
                             m_uiEvolStage[stage][index][node]);
                        /*if(index==fluid::VAR::U_TAU){
                          if(m_uiEvolVarIm[fluid::VAR::U_TAU][node] < 0.0){
                            double D = m_uiEvolVarIm[fluid::VAR::U_D][node];
                            double Sx = m_uiEvolVarIm[fluid::VAR::U_SX][node];
                            double Sy = m_uiEvolVarIm[fluid::VAR::U_SY][node];
                            double Sz = m_uiEvolVarIm[fluid::VAR::U_SZ][node];
                            double tau = m_uiEvolVarIm[fluid::VAR::U_TAU][node];
                            double Ssq = Sx*Sx + Sy*Sy + Sz*Sz;
                            double a = tau + D;

                            std::cout << "Tau is negative after RK update!\n";
                            std::cout << "  D = " << D << "\n";
                            std::cout << "  tau = " << tau << "\n";
                            std::cout << "  rho = " <<
                        m_uiPrimVar[fluid::PVAR::V_RHO][node] << "\n"; std::cout
                        << "  P   = " << m_uiPrimVar[fluid::PVAR::V_P][node] <<
                        "\n"; std::cout << "  Ssq = " << Ssq << "\n"; std::cout
                        << "  tau + D " << a << "\n"; std::cout << "  Ssq - (tau
                        + D)^2 " << Ssq - a*a << "\n"; std::cout << "  Ssq -
                        (tau + D)^2 - D^2 " << Ssq - a*a - D*D << "\n";
                            std::cout << "  RK Stage: " << stage << "\n";
                          }
                        }*/
                    }
                }

                // fluid_math::con_to_prim(m_uiMesh,m_uiEvolStage[stage],m_uiPrimVar);
                fluid_math::con_to_prim(m_uiMesh, m_uiEvolVarIm, m_uiPrimVar,
                                        badBounds);

                // unzip both primitive and evolve variables.
                unzipVars_async(m_uiEvolVarIm, m_uiEvolUnzipVar,
                                fluid::FLUID_NUM_EVOL_VARS);
                unzipPrimVars(m_uiPrimVar, m_uiPrimUnzipVar);

                // enforce boundary conditions
                fluid::fluid_bcs_cons(m_uiMesh, m_uiEvolUnzipVar);
                fluid::fluid_bcs_prim(m_uiMesh, m_uiPrimUnzipVar);

                current_t_adv = current_t + RK4_T[stage + 1] * m_uiT_h;
            }

            // last stage
            // compute fluid rhs (evol vars)
            fluid::fluid_rhs(m_uiMesh, current_t_adv, m_uiEvolUnzipVarRHS,
                             (const DendroScalar**)m_uiEvolUnzipVar,
                             (const DendroScalar**)m_uiPrimUnzipVar);
            zipVars(m_uiEvolUnzipVarRHS,
                    m_uiEvolStage[(fluid::FLUID_RK4_STAGES - 1)],
                    fluid::FLUID_NUM_EVOL_VARS);

            for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                 index++) {
                for (unsigned int node = nodeLocalBegin; node < nodeLocalEnd;
                     node++) {
                    m_uiEvolVar[index][node] = m_uiEvolPrevVar[index][node];
                    for (unsigned int s = 0; s < (fluid::FLUID_RK4_STAGES); s++)
                        m_uiEvolVar[index][node] +=
                            (RK4_C[s] * m_uiT_h *
                             m_uiEvolStage[s][index][node]);
                    if (m_uiEvolVarIm[fluid::VAR::U_TAU][node] < 0.0) {
                        double D   = m_uiEvolVarIm[fluid::VAR::U_D][node];
                        double Sx  = m_uiEvolVarIm[fluid::VAR::U_SX][node];
                        double Sy  = m_uiEvolVarIm[fluid::VAR::U_SY][node];
                        double Sz  = m_uiEvolVarIm[fluid::VAR::U_SZ][node];
                        double tau = m_uiEvolVarIm[fluid::VAR::U_TAU][node];
                        double Ssq = Sx * Sx + Sy * Sy + Sz * Sz;
                        double a   = tau + D;

                        std::cout << "Tau is negative after RK update!\n";
                        std::cout << "  D = " << D << "\n";
                        std::cout << "  tau = " << tau << "\n";
                        std::cout << "  rho = "
                                  << m_uiPrimVar[fluid::PVAR::V_RHO][node]
                                  << "\n";
                        std::cout
                            << "  P   = " << m_uiPrimVar[fluid::PVAR::V_P][node]
                            << "\n";
                        std::cout << "  Ssq = " << Ssq << "\n";
                        std::cout << "  tau + D " << a << "\n";
                        std::cout << "  Ssq - (tau + D)^2 " << Ssq - a * a
                                  << "\n";
                        std::cout << "  Ssq - (tau + D)^2 - D^2 "
                                  << Ssq - a * a - D * D << "\n";
                        std::cout << "  RK Stage: Final \n";
                    }
                }
            }

            fluid_math::con_to_prim(m_uiMesh, m_uiEvolVar, m_uiPrimVar,
                                    badBounds);

        } else if (m_uiRKType ==
                   RKType::RK45) {  // rk45 solver
                                    //@jacob may be you can write this, if not
                                    //tell me I will write this - Milinda.
        }
    }

    m_uiMesh->waitAll();

    m_uiCurrentStep++;
    m_uiCurrentTime += m_uiT_h;
}

void RK_FLUID::rkSolve() {
    if (m_uiCurrentStep == 0) {
        // applyInitialConditions(m_uiEvolPrevVar,m_uiPrimPrevVar);
        initialGridConverge();
    }

    bool isRefine = true;
    unsigned int oldElements, oldElements_g;
    unsigned int newElements, newElements_g;

    // refine based on all the variables
    const unsigned int refineNumVars = fluid::FLUID_NUM_REFINE_VARS;
    unsigned int refineVarIds[refineNumVars];
    for (unsigned int vIndex = 0; vIndex < refineNumVars; vIndex++)
        refineVarIds[vIndex] = fluid::FLUID_REFINE_VARIABLE_INDICES[vIndex];

    double wTol = fluid::FLUID_WAVELET_TOL;
    std::function<double(double, double, double)> waveletTolFunc =
        [wTol](double x, double y, double z) {
            return fluid::computeWTol(x, y, z, wTol);
        };

    double l_min, l_max;
    for (double t = m_uiCurrentTime; t < m_uiTimeEnd; t = t + m_uiT_h) {
        // checkpoint the previous solution value before going to the next step.
        fluid::timer::t_ioCheckPoint.start();
        if ((m_uiMesh->isActive()) &&
            (m_uiCurrentStep % fluid::FLUID_CHECKPT_FREQ) == 0)
            storeCheckPoint(fluid::FLUID_CHKPT_FILE_PREFIX.c_str());
        fluid::timer::t_ioCheckPoint.stop();
        // write sol to vtu.
        if ((m_uiMesh->isActive()) &&
            (m_uiCurrentStep % fluid::FLUID_TIME_STEP_OUTPUT_FREQ) == 0) {
            if (!(m_uiMesh->getMPIRank()))
                std::cout << "executing step: " << m_uiCurrentStep
                          << " dt: " << m_uiT_h
                          << " rk_time : " << m_uiCurrentTime << std::endl;

            for (unsigned int var = 0; var < fluid::FLUID_NUM_EVOL_VARS;
                 var++) {
                l_min =
                    vecMin(m_uiEvolPrevVar[var] + m_uiMesh->getNodeLocalBegin(),
                           (m_uiMesh->getNumLocalMeshNodes()),
                           m_uiMesh->getMPICommunicator());
                l_max =
                    vecMax(m_uiEvolPrevVar[var] + m_uiMesh->getNodeLocalBegin(),
                           (m_uiMesh->getNumLocalMeshNodes()),
                           m_uiMesh->getMPICommunicator());
                if (!(m_uiMesh->getMPIRank())) {
                    std::cout << "\tconserve variable:    ||"
                              << fluid::FLUID_EVOL_VAR_NAMES[var]
                              << "|| (min, max) : (" << l_min << ", " << l_max
                              << " ) " << std::endl;
                }
            }

            if (!(m_uiMesh->getMPIRank())) std::cout << std::endl;

            for (unsigned int var = 0; var < fluid::FLUID_NUM_PRIM_VARS;
                 var++) {
                l_min =
                    vecMin(m_uiPrimPrevVar[var] + m_uiMesh->getNodeLocalBegin(),
                           (m_uiMesh->getNumLocalMeshNodes()),
                           m_uiMesh->getMPICommunicator());
                l_max =
                    vecMax(m_uiPrimPrevVar[var] + m_uiMesh->getNodeLocalBegin(),
                           (m_uiMesh->getNumLocalMeshNodes()),
                           m_uiMesh->getMPICommunicator());
                if (!(m_uiMesh->getMPIRank())) {
                    std::cout << "\tprimitive variable:    ||"
                              << fluid::FLUID_PRIM_VAR_NAMES[var]
                              << "|| (min, max) : (" << l_min << ", " << l_max
                              << " ) " << std::endl;
                }
            }

            fluid::timer::profileInfoIntermediate(
                fluid::FLUID_PROFILE_FILE_PREFIX.c_str(), m_uiMesh,
                m_uiCurrentStep);
        }

        if ((m_uiCurrentStep % fluid::FLUID_TIME_STEP_OUTPUT_FREQ) == 0)
            fluid::timer::resetSnapshot();

        if ((m_uiCurrentStep % fluid::FLUID_REMESH_TEST_FREQ) == 0) {
            unzipVars_async(m_uiEvolPrevVar, m_uiEvolUnzipVar,
                            fluid::FLUID_NUM_EVOL_VARS);
            unzipPrimVars(m_uiPrimPrevVar, m_uiPrimUnzipVar);
            // unzipVars_async(m_uiPrimPrevVar,m_uiPrimUnzipVar,fluid::FLUID_NUM_PRIM_VARS);

            // enforce boundary conditions
            fluid::fluid_bcs_cons(m_uiMesh, m_uiEvolUnzipVar);
            fluid::fluid_bcs_prim(
                m_uiMesh, m_uiPrimUnzipVar);  // FIXME: Possible unnecessary?

#ifdef DEBUG_RK_SOLVER
            if (m_uiMesh->isActive()) {
                if (!m_uiMesh->getMPIRank())
                    std::cout << " isRemesh Unzip : " << std::endl;
                for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                     index++)
                    ot::test::isUnzipNaN(m_uiMesh, m_uiUnzipVar[index]);
            }
#endif

            fluid::timer::t_isReMesh.start();
            std::vector<unsigned int> refine_flags;
            if (fluid::FLUID_ENABLE_BLOCK_ADAPTIVITY)
                isRefine = false;
            else {
                // isRefine=m_uiMesh->isReMeshUnzip((const double
                // **)m_uiEvolUnzipVar,refineVarIds,refineNumVars,waveletTolFunc,fluid::FLUID_DENDRO_AMR_FAC);
                isRefine = fluid::computeFluidRemeshFlags(
                    m_uiMesh, refine_flags, (const double**)m_uiPrimUnzipVar,
                    refineVarIds, refineNumVars, waveletTolFunc,
                    fluid::FLUID_DENDRO_AMR_FAC);
                m_uiMesh->setMeshRefinementFlags(refine_flags);
            }
            fluid::timer::t_isReMesh.stop();

            if (isRefine) {
#ifdef DEBUG_IS_REMESH
                unsigned int rank   = m_uiMesh->getMPIRankGlobal();
                MPI_Comm globalComm = m_uiMesh->getMPIGlobalCommunicator();
                std::vector<ot::TreeNode> unChanged;
                std::vector<ot::TreeNode> refined;
                std::vector<ot::TreeNode> coarsened;
                std::vector<ot::TreeNode> localBlocks;

                const ot::Block* blkList =
                    &(*(m_uiMesh->getLocalBlockList().begin()));
                for (unsigned int ele = 0;
                     ele < m_uiMesh->getLocalBlockList().size(); ele++) {
                    localBlocks.push_back(blkList[ele].getBlockNode());
                }

                const ot::TreeNode* pNodes =
                    &(*(m_uiMesh->getAllElements().begin()));
                for (unsigned int ele = m_uiMesh->getElementLocalBegin();
                     ele < m_uiMesh->getElementLocalEnd(); ele++) {
                    if ((pNodes[ele].getFlag() >> NUM_LEVEL_BITS) ==
                        OCT_NO_CHANGE) {
                        unChanged.push_back(pNodes[ele]);
                    } else if ((pNodes[ele].getFlag() >> NUM_LEVEL_BITS) ==
                               OCT_SPLIT) {
                        refined.push_back(pNodes[ele]);
                    } else {
                        assert((pNodes[ele].getFlag() >> NUM_LEVEL_BITS) ==
                               OCT_COARSE);
                        coarsened.push_back(pNodes[ele]);
                    }
                }

                char fN1[256];
                char fN2[256];
                char fN3[256];
                char fN4[256];

                sprintf(fN1, "unchanged_%d", m_uiCurrentStep);
                sprintf(fN2, "refined_%d", m_uiCurrentStep);
                sprintf(fN3, "coarsend_%d", m_uiCurrentStep);
                sprintf(fN4, "blocks_%d", m_uiCurrentStep);

                DendroIntL localSz = unChanged.size();
                DendroIntL globalSz;
                par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, globalComm);
                if (!rank)
                    std::cout << " total unchanged: " << globalSz << std::endl;

                localSz = refined.size();
                par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, globalComm);
                if (!rank)
                    std::cout << " total refined: " << globalSz << std::endl;

                localSz = coarsened.size();
                par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, globalComm);
                if (!rank)
                    std::cout << " total coarsend: " << globalSz << std::endl;

                io::vtk::oct2vtu(&(*(unChanged.begin())), unChanged.size(), fN1,
                                 globalComm);
                io::vtk::oct2vtu(&(*(refined.begin())), refined.size(), fN2,
                                 globalComm);
                io::vtk::oct2vtu(&(*(coarsened.begin())), coarsened.size(), fN3,
                                 globalComm);
                io::vtk::oct2vtu(&(*(localBlocks.begin())), localBlocks.size(),
                                 fN4, globalComm);

#endif

                fluid::timer::t_mesh.start();
                ot::Mesh* newMesh = m_uiMesh->ReMesh(
                    fluid::FLUID_DENDRO_GRAIN_SZ, fluid::FLUID_LOAD_IMB_TOL,
                    fluid::FLUID_SPLIT_FIX);
                fluid::timer::t_mesh.stop();

                oldElements = m_uiMesh->getNumLocalMeshElements();
                newElements = newMesh->getNumLocalMeshElements();

                par::Mpi_Reduce(&oldElements, &oldElements_g, 1, MPI_SUM, 0,
                                m_uiMesh->getMPIGlobalCommunicator());
                par::Mpi_Reduce(&newElements, &newElements_g, 1, MPI_SUM, 0,
                                newMesh->getMPIGlobalCommunicator());

                if (!(m_uiMesh->getMPIRankGlobal()))
                    std::cout << "step : " << m_uiCurrentStep
                              << " time : " << m_uiCurrentTime
                              << " old mesh: " << oldElements_g
                              << " new mesh: " << newElements_g << std::endl;

                // performs the inter-grid transfer
                intergridTransferVars(m_uiEvolPrevVar, newMesh,
                                      fluid::FLUID_NUM_EVOL_VARS);
                intergridTransferPrimitiveVars(m_uiPrimPrevVar, newMesh);
                fluid_math::prim_to_con(newMesh, m_uiPrimPrevVar,
                                        m_uiEvolPrevVar);

                for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                     index++) {
                    delete[] m_uiEvolVar[index];
                    delete[] m_uiEvolVarIm[index];
                    delete[] m_uiEvolUnzipVar[index];
                    delete[] m_uiEvolUnzipVarRHS[index];

                    m_uiEvolVar[index]         = NULL;
                    m_uiEvolVarIm[index]       = NULL;
                    m_uiEvolUnzipVar[index]    = NULL;
                    m_uiEvolUnzipVarRHS[index] = NULL;

                    m_uiEvolVar[index] = newMesh->createVector<DendroScalar>();
                    m_uiEvolVarIm[index] =
                        newMesh->createVector<DendroScalar>();

                    m_uiEvolUnzipVar[index] =
                        newMesh->createUnZippedVector<DendroScalar>();
                    m_uiEvolUnzipVarRHS[index] =
                        newMesh->createUnZippedVector<DendroScalar>();
                }

                for (unsigned int index = 0; index < fluid::FLUID_NUM_PRIM_VARS;
                     index++) {
                    delete[] m_uiPrimVar[index];
                    delete[] m_uiPrimVarIm[index];
                    delete[] m_uiPrimUnzipVar[index];

                    m_uiPrimVar[index]      = NULL;
                    m_uiPrimVarIm[index]    = NULL;
                    m_uiPrimUnzipVar[index] = NULL;

                    m_uiPrimVar[index] = newMesh->createVector<DendroScalar>();
                    m_uiPrimVarIm[index] =
                        newMesh->createVector<DendroScalar>();

                    m_uiPrimUnzipVar[index] =
                        newMesh->createUnZippedVector<DendroScalar>();
                }

                for (unsigned int stage = 0; stage < m_uiNumRKStages; stage++) {
                    for (unsigned int index = 0;
                         index < fluid::FLUID_NUM_EVOL_VARS; index++) {
                        delete[] m_uiEvolStage[stage][index];
                        m_uiEvolStage[stage][index] = NULL;
                        m_uiEvolStage[stage][index] =
                            newMesh->createVector<DendroScalar>();
                    }
                }

                std::swap(newMesh, m_uiMesh);
                delete newMesh;

                reallocateMPIResources();

                if (m_uiMesh->isActive()) {
                    for (unsigned int var = 0; var < fluid::FLUID_NUM_EVOL_VARS;
                         var++) {
                        l_min = vecMin(m_uiEvolPrevVar[var] +
                                           m_uiMesh->getNodeLocalBegin(),
                                       (m_uiMesh->getNumLocalMeshNodes()),
                                       m_uiMesh->getMPICommunicator());
                        l_max = vecMax(m_uiEvolPrevVar[var] +
                                           m_uiMesh->getNodeLocalBegin(),
                                       (m_uiMesh->getNumLocalMeshNodes()),
                                       m_uiMesh->getMPICommunicator());
                        if (!(m_uiMesh->getMPIRank())) {
                            std::cout << "transfer completed:    ||"
                                      << fluid::FLUID_EVOL_VAR_NAMES[var]
                                      << "|| (min, max) : (" << l_min << ", "
                                      << l_max << " ) " << std::endl;
                        }
                    }

                    for (unsigned int var = 0; var < fluid::FLUID_NUM_PRIM_VARS;
                         var++) {
                        l_min = vecMin(m_uiPrimPrevVar[var] +
                                           m_uiMesh->getNodeLocalBegin(),
                                       (m_uiMesh->getNumLocalMeshNodes()),
                                       m_uiMesh->getMPICommunicator());
                        l_max = vecMax(m_uiPrimPrevVar[var] +
                                           m_uiMesh->getNodeLocalBegin(),
                                       (m_uiMesh->getNumLocalMeshNodes()),
                                       m_uiMesh->getMPICommunicator());
                        if (!(m_uiMesh->getMPIRank())) {
                            std::cout << "transfer completed:    ||"
                                      << fluid::FLUID_PRIM_VAR_NAMES[var]
                                      << "|| (min, max) : (" << l_min << ", "
                                      << l_max << " ) " << std::endl;
                        }
                    }
                }
            }
        }

        if ((m_uiCurrentStep % fluid::FLUID_IO_OUTPUT_FREQ) == 0) {
            for (unsigned int var = 0; var < fluid::FLUID_NUM_EVOL_VARS; var++)
                m_uiMesh->performGhostExchange(m_uiEvolPrevVar[var]);

            for (unsigned int var = 0; var < fluid::FLUID_NUM_PRIM_VARS; var++)
                m_uiMesh->performGhostExchange(m_uiPrimPrevVar[var]);

            const unsigned int evarIndices[] = {0, 1, 2, 3, 4, 5};
            const unsigned int pvarIndices[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

            writeToVTU(m_uiEvolPrevVar, m_uiPrimPrevVar, NULL,
                       fluid::FLUID_NUM_EVOL_VARS, fluid::FLUID_NUM_PRIM_VARS,
                       0, evarIndices, pvarIndices, NULL,
                       fluid::FLUID_VTU_Z_SLICE_ONLY);
        }

        fluid::timer::t_rkStep.start();
        performSingleIteration();
        fluid::timer::t_rkStep.stop();

        std::swap(m_uiEvolVar, m_uiEvolPrevVar);
        std::swap(m_uiPrimVar, m_uiPrimPrevVar);
    }
}

void RK_FLUID::storeCheckPoint(const char* fNamePrefix) {
    if (m_uiMesh->isActive()) {
        unsigned int cpIndex;
        (m_uiCurrentStep % (2 * fluid::FLUID_CHECKPT_FREQ) == 0)
            ? cpIndex = 0
            : cpIndex = 1;  // to support alternate file writing.
        unsigned int rank = m_uiMesh->getMPIRank();
        unsigned int npes = m_uiMesh->getMPICommSize();

        char fName[256];
        const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin() +
                                         m_uiMesh->getElementLocalBegin()));
        sprintf(fName, "%s_octree_%d_%d.oct", fNamePrefix, cpIndex, rank);
        io::checkpoint::writeOctToFile(fName, pNodes,
                                       m_uiMesh->getNumLocalMeshElements());

        unsigned int numEvolVars  = fluid::FLUID_NUM_EVOL_VARS;
        const char** varEvolNames = fluid::FLUID_EVOL_VAR_NAMES;

        for (unsigned int i = 0; i < numEvolVars; i++) {
            sprintf(fName, "%s_%s_%d_%d.var", fNamePrefix, varEvolNames[i],
                    cpIndex, rank);
            io::checkpoint::writeVecToFile(fName, m_uiMesh, m_uiEvolPrevVar[i]);
        }

        unsigned int numPrimVars  = fluid::FLUID_NUM_PRIM_VARS;
        const char** varPrimNames = fluid::FLUID_PRIM_VAR_NAMES;

        for (unsigned int i = 0; i < numPrimVars; i++) {
            sprintf(fName, "%s_%s_%d_%d.var", fNamePrefix, varPrimNames[i],
                    cpIndex, rank);
            io::checkpoint::writeVecToFile(fName, m_uiMesh, m_uiPrimPrevVar[i]);
        }

        sprintf(fName, "%s_%d_%d.var", fNamePrefix, cpIndex, rank);
        io::checkpoint::writeVecToFile(fName, m_uiMesh,
                                       (const double**)m_uiPrimPrevVar,
                                       fluid::FLUID_NUM_EVOL_VARS);

        if (!rank) {
            sprintf(fName, "%s_step_%d.cp", fNamePrefix, cpIndex);
            // std::cout<<"writing : "<<fName<<std::endl;
            std::ofstream outfile(fName);
            if (!outfile) {
                std::cout << fName << " file open failed " << std::endl;
                return;
            }

            json checkPoint;
            checkPoint["DENDRO_RK_TIME_BEGIN"]     = m_uiTimeBegin;
            checkPoint["DENDRO_RK_TIME_END"]       = m_uiTimeEnd;
            checkPoint["DENDRO_RK_ELEMENT_ORDER"]  = m_uiOrder;

            checkPoint["DENDRO_RK_TIME_CURRENT"]   = m_uiCurrentTime;
            checkPoint["DENDRO_RK_STEP_CURRENT"]   = m_uiCurrentStep;
            checkPoint["DENDRO_RK_TIME_STEP_SIZE"] = m_uiT_h;
            checkPoint["DENDRO_RK_LAST_IO_TIME"]   = m_uiCurrentTime;

            checkPoint["DENDRO_RK_WAVELET_TOLERANCE"] =
                fluid::FLUID_WAVELET_TOL;
            checkPoint["DENDRO_RK_LOAD_IMB_TOLERANCE"] =
                fluid::FLUID_LOAD_IMB_TOL;
            checkPoint["DENDRO_RK_NUM_EVOL_VARS"] = numEvolVars;
            checkPoint["DENDRO_RK_NUM_PRIM_VARS"] = numPrimVars;
            checkPoint["DENDRO_RK45_ACTIVE_COMM_SZ"] =
                m_uiMesh
                    ->getMPICommSize();  // (note that rank 0 is always active).

            outfile << std::setw(4) << checkPoint << std::endl;
            outfile.close();
        }
    }
}

void RK_FLUID::restoreCheckPoint(const char* fNamePrefix, MPI_Comm comm) {
    unsigned int numPrimVars = 0;
    unsigned int numEvolVars = 0;
    std::vector<ot::TreeNode> octree;
    json checkPoint;

    int rank;
    int npes;
    m_uiComm = comm;
    MPI_Comm_rank(m_uiComm, &rank);
    MPI_Comm_size(m_uiComm, &npes);

    unsigned int activeCommSz;

    char fName[256];
    unsigned int restoreStatus = 0;
    unsigned int restoreStatusGlobal =
        0;  // 0 indicates successfully restorable.

    ot::Mesh* newMesh;

    for (unsigned int cpIndex = 0; cpIndex < 2; cpIndex++) {
        restoreStatus = 0;
        octree.clear();
        if (!rank)
            std::cout << " Trying to restore from checkpoint index : "
                      << cpIndex << std::endl;

        if (!rank) {
            sprintf(fName, "%s_step_%d.cp", fNamePrefix, cpIndex);
            std::ifstream infile(fName);
            if (!infile) {
                std::cout << fName << " file open failed " << std::endl;
                restoreStatus = 1;
            }

            if (restoreStatus == 0) {
                infile >> checkPoint;
                m_uiTimeBegin   = checkPoint["DENDRO_RK_TIME_BEGIN"];
                m_uiTimeEnd     = checkPoint["DENDRO_RK_TIME_END"];

                m_uiOrder       = checkPoint["DENDRO_RK_ELEMENT_ORDER"];
                m_uiCurrentTime = checkPoint["DENDRO_RK_TIME_CURRENT"];
                m_uiCurrentStep = checkPoint["DENDRO_RK_STEP_CURRENT"];
                m_uiT_h         = checkPoint["DENDRO_RK_TIME_STEP_SIZE"];

                fluid::FLUID_WAVELET_TOL =
                    checkPoint["DENDRO_RK_WAVELET_TOLERANCE"];
                fluid::FLUID_LOAD_IMB_TOL =
                    checkPoint["DENDRO_RK_LOAD_IMB_TOLERANCE"];
                numEvolVars  = checkPoint["DENDRO_RK_NUM_EVOL_VARS"];
                numPrimVars  = checkPoint["DENDRO_RK_NUM_PRIM_VARS"];
                activeCommSz = checkPoint["DENDRO_RK_ACTIVE_COMM_SZ"];
            }
        }

        par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX,
                           m_uiComm);
        if (restoreStatusGlobal == 1) continue;

        par::Mpi_Bcast(&m_uiTimeBegin, 1, 0, comm);
        par::Mpi_Bcast(&m_uiTimeEnd, 1, 0, comm);

        par::Mpi_Bcast(&m_uiCurrentTime, 1, 0, comm);
        par::Mpi_Bcast(&m_uiCurrentStep, 1, 0, comm);

        par::Mpi_Bcast(&m_uiT_h, 1, 0, comm);

        par::Mpi_Bcast(&fluid::FLUID_WAVELET_TOL, 1, 0, comm);
        par::Mpi_Bcast(&fluid::FLUID_LOAD_IMB_TOL, 1, 0, comm);

        par::Mpi_Bcast(&numEvolVars, 1, 0, comm);
        par::Mpi_Bcast(&numPrimVars, 1, 0, comm);
        par::Mpi_Bcast(&m_uiOrder, 1, 0, comm);
        par::Mpi_Bcast(&m_uiT_h, 1, 0, comm);

        par::Mpi_Bcast(&activeCommSz, 1, 0, comm);

        if (activeCommSz > npes) {
            if (!rank)
                std::cout << " [Error] : checkpoint file written from  a "
                             "larger communicator than the current global "
                             "comm. (i.e. communicator shrinking not allowed "
                             "in the restore step. )"
                          << std::endl;
            exit(0);
        }

        bool isActive = (rank < activeCommSz);

        MPI_Comm newComm;
        par::splitComm2way(isActive, &newComm, m_uiComm);

        if (isActive) {
            int activeRank;
            int activeNpes;

            MPI_Comm_rank(newComm, &activeRank);
            MPI_Comm_size(newComm, &activeNpes);
            assert(activeNpes == activeCommSz);

            sprintf(fName, "%s_octree_%d_%d.oct", fNamePrefix, cpIndex,
                    activeRank);
            restoreStatus = io::checkpoint::readOctFromFile(fName, octree);
            assert(par::test::isUniqueAndSorted(octree, newComm));
        }

        par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX,
                           m_uiComm);
        if (restoreStatusGlobal == 1) {
            if (!rank)
                std::cout << "[Error]: octree (*.oct) restore file currupted "
                          << std::endl;
            continue;
        }

        newMesh = new ot::Mesh(octree, 1, m_uiOrder, activeCommSz, m_uiComm);

        for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
             index++) {
            delete[] m_uiEvolVar[index];
            delete[] m_uiEvolVarIm[index];
            delete[] m_uiEvolUnzipVar[index];
            delete[] m_uiEvolUnzipVarRHS[index];

            m_uiEvolVar[index]         = NULL;
            m_uiEvolVarIm[index]       = NULL;
            m_uiEvolUnzipVar[index]    = NULL;
            m_uiEvolUnzipVarRHS[index] = NULL;

            m_uiEvolVar[index]         = newMesh->createVector<DendroScalar>();
            m_uiEvolVarIm[index]       = newMesh->createVector<DendroScalar>();

            m_uiEvolUnzipVar[index] =
                newMesh->createUnZippedVector<DendroScalar>();
            m_uiEvolUnzipVarRHS[index] =
                newMesh->createUnZippedVector<DendroScalar>();
        }

        for (unsigned int index = 0; index < fluid::FLUID_NUM_PRIM_VARS;
             index++) {
            delete[] m_uiPrimVar[index];
            delete[] m_uiPrimVarIm[index];
            delete[] m_uiPrimUnzipVar[index];

            m_uiPrimVar[index]      = NULL;
            m_uiPrimVarIm[index]    = NULL;
            m_uiPrimUnzipVar[index] = NULL;

            m_uiPrimVar[index]      = newMesh->createVector<DendroScalar>();
            m_uiPrimVarIm[index]    = newMesh->createVector<DendroScalar>();

            m_uiPrimUnzipVar[index] =
                newMesh->createUnZippedVector<DendroScalar>();
        }

        for (unsigned int stage = 0; stage < m_uiNumRKStages; stage++) {
            for (unsigned int index = 0; index < fluid::FLUID_NUM_EVOL_VARS;
                 index++) {
                delete[] m_uiEvolStage[stage][index];
                m_uiEvolStage[stage][index] = NULL;
                m_uiEvolStage[stage][index] =
                    newMesh->createVector<DendroScalar>();
            }
        }

        std::swap(newMesh, m_uiMesh);
        delete newMesh;

        reallocateMPIResources();

        const char** varEvolNames = fluid::FLUID_EVOL_VAR_NAMES;
        const char** varPrimNames = fluid::FLUID_PRIM_VAR_NAMES;

        if (isActive) {
            int activeRank;
            int activeNpes;

            MPI_Comm_rank(newComm, &activeRank);
            MPI_Comm_size(newComm, &activeNpes);
            assert(activeNpes == activeCommSz);

            /*for(unsigned int i=0;i<numVars;i++)
            {
                sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,activeRank);
                restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,m_uiPrevVar[i]);
                if(restoreStatus==1) break;

            }*/
            sprintf(fName, "%s_%d_%d.var", fNamePrefix, cpIndex, activeRank);
            restoreStatus = io::checkpoint::readVecFromFile(
                fName, newMesh, m_uiEvolPrevVar, fluid::FLUID_NUM_EVOL_VARS);
            restoreStatus = io::checkpoint::readVecFromFile(
                fName, newMesh, m_uiPrimPrevVar, fluid::FLUID_NUM_PRIM_VARS);
        }

        MPI_Comm_free(&newComm);
        par::Mpi_Allreduce(&restoreStatus, &restoreStatusGlobal, 1, MPI_MAX,
                           m_uiComm);
        if (restoreStatusGlobal == 1) {
            if (!rank)
                std::cout << "[Error]: varible (*.var) restore file currupted "
                          << std::endl;
            continue;
        }

        std::swap(m_uiMesh, newMesh);
        delete newMesh;
        reallocateMPIResources();
        if (restoreStatusGlobal == 0) break;
    }

    if (restoreStatusGlobal == 1) {
        std::cout << "rank: " << rank << "[Error]: rk solver restore error "
                  << std::endl;
        exit(0);
    }

    unsigned int localSz = m_uiMesh->getNumLocalMeshElements();
    unsigned int totalElems;

    par::Mpi_Allreduce(&localSz, &totalElems, 1, MPI_SUM, m_uiComm);

    if (!rank)
        std::cout << " checkpoint at step : " << m_uiCurrentStep
                  << "active Comm. sz: " << activeCommSz
                  << " restore successful: " << " restored mesh size: "
                  << totalElems << std::endl;
    return;
}

}  // end of namespace solver
}  // end of namespace ode
