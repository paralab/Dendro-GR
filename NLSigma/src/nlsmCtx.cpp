#include "nlsmCtx.h"

namespace nlsm {
NLSMCtx::NLSMCtx(ot::Mesh* pMesh) : Ctx() {
    m_uiMesh = pMesh;
    // variable allocation for evolution variables
    m_evar   = DVec();
    m_evar.create_vector(m_uiMesh, ot::DVEC_TYPE::OCT_SHARED_NODES,
                         ot::DVEC_LOC::HOST, NLSM_NUM_VARS, true);

    m_evar_unzip[0] = DVec();
    m_evar_unzip[1] = DVec();

    m_evar_unzip[0].create_vector(m_uiMesh,
                                  ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                  ot::DVEC_LOC::HOST, NLSM_NUM_VARS, true);
    m_evar_unzip[1].create_vector(m_uiMesh,
                                  ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                  ot::DVEC_LOC::HOST, NLSM_NUM_VARS, true);

    m_uiTinfo._m_uiStep = 0;
    m_uiTinfo._m_uiT    = 0;
    m_uiTinfo._m_uiTb   = NLSM_RK45_TIME_BEGIN;
    m_uiTinfo._m_uiTe   = NLSM_RK45_TIME_END;
    m_uiTinfo._m_uiTh   = NLSM_RK45_TIME_STEP_SIZE;

    m_uiElementOrder    = NLSM_ELE_ORDER;

    m_uiMinPt = Point(NLSM_GRID_MIN_X, NLSM_GRID_MIN_Y, NLSM_GRID_MIN_Z);
    m_uiMaxPt = Point(NLSM_GRID_MAX_X, NLSM_GRID_MAX_Y, NLSM_GRID_MAX_Z);

    ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, NLSM_NUM_VARS,
                                      NLSM_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, NLSM_NUM_VARS,
                                    NLSM_ASYNC_COMM_K);

    return;
}

NLSMCtx::~NLSMCtx() {
    m_evar.destroy_vector();
    m_evar_unzip[0].destroy_vector();
    m_evar_unzip[1].destroy_vector();

    ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, NLSM_NUM_VARS,
                                      NLSM_ASYNC_COMM_K);
}

int NLSMCtx::initialize() {
    if (NLSM_RESTORE_SOLVER) {
        this->restore_checkpt();
        return 0;
    }

    unsigned int nodeLookUp_CG;
    unsigned int nodeLookUp_DG;
    DendroScalar x, y, z, len;
    const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin()));
    unsigned int ownerID, ii_x, jj_y, kk_z;
    unsigned int eleOrder      = m_uiMesh->getElementOrder();
    const unsigned int* e2n_cg = &(*(m_uiMesh->getE2NMapping().begin()));
    const unsigned int* e2n_dg = &(*(m_uiMesh->getE2NMapping_DG().begin()));
    const unsigned int nPe     = m_uiMesh->getNumNodesPerElement();
    const unsigned int nodeLocalBegin = m_uiMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd   = m_uiMesh->getNodeLocalEnd();

    DendroScalar var[NLSM_NUM_VARS];
    DendroScalar* zipIn[NLSM_NUM_VARS];

    m_evar.to_2d(zipIn);

    DendroScalar mp, mm, mp_adm, mm_adm, E, J1, J2, J3;

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
                        x = pNodes[ownerID].getX() + ii_x * (len / (eleOrder));
                        y = pNodes[ownerID].getY() + jj_y * (len / (eleOrder));
                        z = pNodes[ownerID].getZ() + kk_z * (len / (eleOrder));

                        initData(x, y, z, var);

                        for (unsigned int v = 0; v < NLSM_NUM_VARS; v++)
                            zipIn[v][nodeLookUp_CG] = var[v];
                    }
                }
    }

    return 0;
}

int NLSMCtx::finalize() { return 0; }

int NLSMCtx::rhs(DVec* in, DVec* out, unsigned int sz, DendroScalar time) {
    // all the variables should be packed together.
    assert(sz == 1);

    DendroScalar* sVar[NLSM_NUM_VARS];
    in[0].to_2d(sVar);

    this->unzip(in[0], m_evar_unzip[0], NLSM_ASYNC_COMM_K);

    DendroScalar* unzipIn[NLSM_NUM_VARS];
    DendroScalar* unzipOut[NLSM_NUM_VARS];

    m_evar_unzip[0].to_2d(unzipIn);
    m_evar_unzip[1].to_2d(unzipOut);

    const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int lsz[3];
    unsigned int bflag;
    double dx, dy, dz;

    const Point pt_min(nlsm::NLSM_COMPD_MIN[0], nlsm::NLSM_COMPD_MIN[1],
                       nlsm::NLSM_COMPD_MIN[2]);
    const Point pt_max(nlsm::NLSM_COMPD_MAX[0], nlsm::NLSM_COMPD_MAX[1],
                       nlsm::NLSM_COMPD_MAX[2]);
    const unsigned int PW = nlsm::NLSM_PADDING_WIDTH;

    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        offset   = blkList[blk].getOffset();
        lsz[0]   = blkList[blk].getAllocationSzX();
        lsz[1]   = blkList[blk].getAllocationSzY();
        lsz[2]   = blkList[blk].getAllocationSzZ();

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

#ifdef __PROFILE_CTX__
        m_uiCtxpt[ts::CTXPROFILE::RHS].start();
#endif

        nlsmRhs(unzipOut, (const DendroScalar**)unzipIn, offset, ptmin, ptmax,
                lsz, bflag);

        // std::cout<<":::\n";
        // for(unsigned int n=0; n < lsz[0]*lsz[1]*lsz[2]; n++)
        //     std::cout<<"unzipOut: "<<n<<" val: "<<unzipOut[0][n]<<std::endl;

#ifdef __PROFILE_CTX__
        m_uiCtxpt[ts::CTXPROFILE::RHS].stop();
#endif
    }

    this->zip(m_evar_unzip[1], out[0]);

    return 0;
}

int NLSMCtx::rhs_blk(const DendroScalar* in, DendroScalar* out,
                     unsigned int dof, unsigned int local_blk_id,
                     DendroScalar blk_time) {
#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::RHS_BLK].start();
#endif

    // all the variables should be packed together.
    const ot::Block* blkList     = m_uiMesh->getLocalBlockList().data();
    const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
    assert(local_blk_id < numBlocks);

    const unsigned int blk  = local_blk_id;
    DendroScalar** unzipIn  = new DendroScalar*[dof];
    DendroScalar** unzipOut = new DendroScalar*[dof];

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int lsz[3];
    unsigned int bflag;
    double dx, dy, dz;

    const Point pt_min(nlsm::NLSM_COMPD_MIN[0], nlsm::NLSM_COMPD_MIN[1],
                       nlsm::NLSM_COMPD_MIN[2]);
    const Point pt_max(nlsm::NLSM_COMPD_MAX[0], nlsm::NLSM_COMPD_MAX[1],
                       nlsm::NLSM_COMPD_MAX[2]);
    const unsigned int PW = nlsm::NLSM_PADDING_WIDTH;

    offset                = blkList[blk].getOffset();
    lsz[0]                = blkList[blk].getAllocationSzX();
    lsz[1]                = blkList[blk].getAllocationSzY();
    lsz[2]                = blkList[blk].getAllocationSzZ();

    const unsigned int NN = lsz[0] * lsz[1] * lsz[2];

    for (unsigned int v = 0; v < dof; v++) {
        unzipIn[v]  = (DendroScalar*)(in + v * NN);
        unzipOut[v] = (DendroScalar*)(out + v * NN);
    }

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

    // note that the offset zero is important since it is the block vector.
    nlsmRhs(unzipOut, (const DendroScalar**)unzipIn, 0, ptmin, ptmax, lsz,
            bflag);

    // for(unsigned int v =0; v < dof; v++)
    // {
    //     unsigned int nid=0;
    //     for(unsigned int k=3; k < lsz[2]-3; k++)
    //     for(unsigned int j=3; j < lsz[1]-3; j++)
    //     for(unsigned int i=3; i < lsz[0]-3; i++,nid++)
    //     {
    //         std::cout<<" blk::v: "<<v<<" n: "<<nid<< " val:
    //         "<<unzipOut[v][k*lsz[1]*lsz[0] +j*lsz[0] +i ]<<std::endl;;
    //     }
    // }

    delete[] unzipIn;
    delete[] unzipOut;

#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::RHS_BLK].stop();
#endif

    return 0;
}

int NLSMCtx::write_vtu() {
    DendroScalar* evolUnzipVar = NULL;
    DendroScalar* consUnzipVar = NULL;
    DendroScalar* consVar      = NULL;
    DendroScalar* evolVar[NLSM_NUM_VARS];

    m_uiMesh->readFromGhostBegin(m_evar.get_vec_ptr(), m_evar.get_dof());
    m_uiMesh->readFromGhostEnd(m_evar.get_vec_ptr(), m_evar.get_dof());

    m_evar.to_2d(evolVar);

    std::vector<std::string> pDataNames;
    const unsigned int numConstVars = 0;
    const unsigned int numEvolVars  = NLSM_NUM_EVOL_VARS_VTU_OUTPUT;

    double* pData[(numConstVars + numEvolVars)];

    for (unsigned int i = 0; i < numEvolVars; i++) {
        pDataNames.push_back(
            std::string(NLSM_VAR_NAMES[NLSM_VTU_OUTPUT_EVOL_INDICES[i]]));
        pData[i] = evolVar[NLSM_VTU_OUTPUT_EVOL_INDICES[i]];
    }

    std::vector<char*> pDataNames_char;
    pDataNames_char.reserve(pDataNames.size());

    for (unsigned int i = 0; i < pDataNames.size(); i++)
        pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

    const char* fDataNames[] = {"Time", "Cycle"};
    const double fData[]     = {m_uiTinfo._m_uiT, (double)m_uiTinfo._m_uiStep};

    char fPrefix[256];
    sprintf(fPrefix, "%s_%d", NLSM_VTU_FILE_PREFIX.c_str(),
            m_uiTinfo._m_uiStep);

    io::vtk::mesh2vtuFine(
        m_uiMesh, fPrefix, 2, fDataNames, fData, (numEvolVars + numConstVars),
        (const char**)&pDataNames_char[0], (const double**)pData);

    return 0;
}

int NLSMCtx::write_checkpt() {
    if (m_uiMesh->isActive()) {
        unsigned int cpIndex;
        (m_uiTinfo._m_uiStep % (2 * NLSM_CHECKPT_FREQ) == 0)
            ? cpIndex = 0
            : cpIndex = 1;  // to support alternate file writing.
        unsigned int rank = m_uiMesh->getMPIRank();
        unsigned int npes = m_uiMesh->getMPICommSize();

        char fName[256];
        const ot::TreeNode* pNodes = &(*(m_uiMesh->getAllElements().begin() +
                                         m_uiMesh->getElementLocalBegin()));
        sprintf(fName, "%s_octree_%d_%d.oct", NLSM_CHKPT_FILE_PREFIX.c_str(),
                cpIndex, rank);
        io::checkpoint::writeOctToFile(fName, pNodes,
                                       m_uiMesh->getNumLocalMeshElements());

        unsigned int numVars   = nlsm::NLSM_NUM_VARS;
        const char** varNames  = nlsm::NLSM_VAR_NAMES;

        /*for(unsigned int i=0;i<numVars;i++)
        {
            sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,rank);
            io::checkpoint::writeVecToFile(fName,m_uiMesh,m_uiPrevVar[i]);
        }*/

        const unsigned int dof = m_evar.get_dof();
        DendroScalar* eVar[dof];
        m_evar.to_2d(eVar);

        sprintf(fName, "%s_%d_%d.var", NLSM_CHKPT_FILE_PREFIX.c_str(), cpIndex,
                rank);
        io::checkpoint::writeVecToFile(fName, m_uiMesh, (const double**)eVar,
                                       nlsm::NLSM_NUM_VARS);

        if (!rank) {
            sprintf(fName, "%s_step_%d.cp", NLSM_CHKPT_FILE_PREFIX.c_str(),
                    cpIndex);
            // std::cout<<"writing : "<<fName<<std::endl;
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

            checkPoint["DENDRO_TS_WAVELET_TOLERANCE"]  = NLSM_WAVELET_TOL;
            checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"] = NLSM_LOAD_IMB_TOL;
            checkPoint["DENDRO_TS_NUM_VARS"] =
                numVars;  // number of variables to restore.
            checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"] =
                m_uiMesh
                    ->getMPICommSize();  // (note that rank 0 is always active).

            outfile << std::setw(4) << checkPoint << std::endl;
            outfile.close();
        }
    }

    return 0;
}

int NLSMCtx::restore_checkpt() {
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
            sprintf(fName, "%s_step_%d.cp", NLSM_CHKPT_FILE_PREFIX.c_str(),
                    cpIndex);
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

                NLSM_WAVELET_TOL    = checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                NLSM_LOAD_IMB_TOL = checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];

                numVars           = checkPoint["DENDRO_TS_NUM_VARS"];
                activeCommSz      = checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];

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

    if (!rank) {
        sprintf(fName, "%s_step_%d.cp", NLSM_CHKPT_FILE_PREFIX.c_str(),
                restoreFileIndex);
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

            NLSM_WAVELET_TOL    = checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
            NLSM_LOAD_IMB_TOL   = checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];

            numVars             = checkPoint["DENDRO_TS_NUM_VARS"];
            activeCommSz        = checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];

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
    }

    MPI_Bcast(&m_uiTinfo, sizeof(ts::TSInfo), MPI_BYTE, 0, comm);
    par::Mpi_Bcast(&NLSM_WAVELET_TOL, 1, 0, comm);
    par::Mpi_Bcast(&NLSM_LOAD_IMB_TOL, 1, 0, comm);

    par::Mpi_Bcast(&numVars, 1, 0, comm);
    par::Mpi_Bcast(&m_uiElementOrder, 1, 0, comm);
    par::Mpi_Bcast(&activeCommSz, 1, 0, comm);

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

        sprintf(fName, "%s_octree_%d_%d.oct", NLSM_CHKPT_FILE_PREFIX.c_str(),
                restoreFileIndex, activeRank);
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
    // no need to transfer data only to resize the contex variables.
    this->grid_transfer(newMesh);

    // only reads the evolution variables.
    if (isActive) {
        int activeRank;
        int activeNpes;

        DendroScalar* inVec[NLSM_NUM_VARS];
        m_evar.to_2d(inVec);

        MPI_Comm_rank(newComm, &activeRank);
        MPI_Comm_size(newComm, &activeNpes);
        assert(activeNpes == activeCommSz);

        sprintf(fName, "%s_%d_%d.var", NLSM_CHKPT_FILE_PREFIX.c_str(),
                restoreFileIndex, activeRank);
        restoreStatus = io::checkpoint::readVecFromFile(fName, newMesh, inVec,
                                                        NLSM_NUM_VARS);
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

    unsigned int localSz    = m_uiMesh->getNumLocalMeshElements();
    unsigned int totalElems = 0;
    par::Mpi_Allreduce(&localSz, &totalElems, 1, MPI_SUM, comm);

    if (!rank)
        std::cout << " checkpoint at step : " << m_uiTinfo._m_uiStep
                  << "active Comm. sz: " << activeCommSz
                  << " restore successful: " << " restored mesh size: "
                  << totalElems << std::endl;

    m_uiIsETSSynced = false;
    return 0;
}

int NLSMCtx::pre_timestep(DVec sIn) { return 0; }

int NLSMCtx::pre_stage(DVec sIn) { return 0; }

int NLSMCtx::post_stage(DVec sIn) { return 0; }

int NLSMCtx::post_timestep(DVec sIn) { return 0; }

bool NLSMCtx::is_remesh() {
#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::IS_REMESH].start();
#endif

    bool isRefine = false;

    if (NLSM_ENABLE_BLOCK_ADAPTIVITY) return isRefine;

    MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();
    this->unzip(m_evar, m_evar_unzip[0], NLSM_ASYNC_COMM_K);

    DendroScalar* unzipVar[NLSM_NUM_VARS];
    m_evar_unzip[0].to_2d(unzipVar);

    unsigned int refineVarIds[NLSM_NUM_REFINE_VARS];
    for (unsigned int vIndex = 0; vIndex < NLSM_NUM_REFINE_VARS; vIndex++)
        refineVarIds[vIndex] = NLSM_REFINE_VARIABLE_INDICES[vIndex];

    double wTol = NLSM_WAVELET_TOL;
    std::function<double(double, double, double, double*)> waveletTolFunc =
        [wTol](double x, double y, double z, double* hx) {
            return computeWTol(x, y, z, hx);
        };

    isRefine = m_uiMesh->isReMeshUnzip((const double**)unzipVar, refineVarIds,
                                       NLSM_NUM_REFINE_VARS, waveletTolFunc,
                                       NLSM_DENDRO_AMR_FAC);

#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::IS_REMESH].stop();
#endif

    return isRefine;
}

DVec& NLSMCtx::get_evolution_vars() { return m_evar; }

DVec& NLSMCtx::get_constraint_vars() { return m_cvar; }

DVec& NLSMCtx::get_primitive_vars() { return m_pvar; }

int NLSMCtx::terminal_output() {
    if (m_uiMesh->isActive()) {
        for (unsigned int v = 0; v < m_evar.get_dof(); v++) {
            DendroScalar min = 0, max = 0;
            min = vecMin(m_uiMesh, m_evar.get_vec_ptr(), ot::VEC_TYPE::CG_NODAL,
                         true);
            max = vecMax(m_uiMesh, m_evar.get_vec_ptr(), ot::VEC_TYPE::CG_NODAL,
                         true);
            if (!(m_uiMesh->getMPIRank()))
                std::cout << "[NLSMCtx]:  " << NLSM_VAR_NAMES[v]
                          << " (min,max) : \t ( " << min << ", " << max << " ) "
                          << std::endl;
        }

#ifdef NLSM_COMPARE_WITH_ANALYTICAL_SOL
        double* evolVar[2];
        m_uiEVar.Get2DVec(evolVar);
        double* chiAnalytical = m_uiMesh->createVector<double>();
        double* diffVec       = m_uiMesh->createVector<double>();

        std::function<void(double, double, double, double, double*)> u_x_t =
            [](double x, double y, double z, double t, double* var) {
                nlsm::analyticalSol(x, y, z, t, var);
            };

        // initialize diff begin.
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

        double var[2];

        double mp, mm, mp_adm, mm_adm, E, J1, J2, J3;

        for (unsigned int elem = m_uiMesh->getElementLocalBegin();
             elem < m_uiMesh->getElementLocalEnd(); elem++) {
            for (unsigned int k = 0; k < (eleOrder + 1); k++)
                for (unsigned int j = 0; j < (eleOrder + 1); j++)
                    for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                        nodeLookUp_CG =
                            e2n_cg[elem * nPe +
                                   k * (eleOrder + 1) * (eleOrder + 1) +
                                   j * (eleOrder + 1) + i];
                        if (nodeLookUp_CG >= nodeLocalBegin &&
                            nodeLookUp_CG < nodeLocalEnd) {
                            nodeLookUp_DG =
                                e2n_dg[elem * nPe +
                                       k * (eleOrder + 1) * (eleOrder + 1) +
                                       j * (eleOrder + 1) + i];
                            m_uiMesh->dg2eijk(nodeLookUp_DG, ownerID, ii_x,
                                              jj_y, kk_z);
                            len = (double)(1u << (m_uiMaxDepth -
                                                  pNodes[ownerID].getLevel()));
                            x   = pNodes[ownerID].getX() +
                                ii_x * (len / (eleOrder));
                            y = pNodes[ownerID].getY() +
                                jj_y * (len / (eleOrder));
                            z = pNodes[ownerID].getZ() +
                                kk_z * (len / (eleOrder));

                            u_x_t((double)x, (double)y, (double)z,
                                  m_uiTinfo._m_uiT, var);
                            diffVec[nodeLookUp_CG] =
                                var[nlsm::VAR::U_CHI] -
                                evolVar[nlsm::VAR::U_CHI][nodeLookUp_CG];
                            chiAnalytical[nodeLookUp_CG] =
                                var[nlsm::VAR::U_CHI];
                        }
                    }
        }

        m_uiMesh->performGhostExchange(diffVec);
        m_uiMesh->performGhostExchange(chiAnalytical);

        double l_rs          = rsNormLp<double>(m_uiMesh, diffVec, 2);
        double l_min         = vecMin(diffVec + m_uiMesh->getNodeLocalBegin(),
                                      (m_uiMesh->getNumLocalMeshNodes()),
                                      m_uiMesh->getMPICommunicator());
        double l_max         = vecMax(diffVec + m_uiMesh->getNodeLocalBegin(),
                                      (m_uiMesh->getNumLocalMeshNodes()),
                                      m_uiMesh->getMPICommunicator());
        double l2_norm       = normL2(diffVec + m_uiMesh->getNodeLocalBegin(),
                                      (m_uiMesh->getNumLocalMeshNodes()),
                                      m_uiMesh->getMPICommunicator());
        DendroIntL local_dof = m_uiMesh->getNumLocalMeshNodes();
        DendroIntL total_dof = 0;
        par::Mpi_Reduce(&local_dof, &total_dof, 1, MPI_SUM, 0,
                        m_uiMesh->getMPICommunicator());

        if (!m_uiMesh->getMPIRank()) {
            // std::cout << "executing step: " << m_uiCurrentStep << " dt: " <<
            // m_uiT_h << " rk_time : "<< m_uiCurrentTime << std::endl;
            l2_norm =
                sqrt((l2_norm * l2_norm) / (double)(total_dof * total_dof));
            std::cout << YLW << "\t ||VAR::DIFF|| (min, max,l2,l_2rs) : ("
                      << l_min << ", " << l_max << ", " << l2_norm << ", "
                      << l_rs << " ) " << NRM << std::endl;

            std::ofstream fileGW;
            char fName[256];
            sprintf(fName, "%s_error.dat",
                    nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
            fileGW.open(fName, std::ofstream::app);
            // writes the header
            if (m_uiTinfo._m_uiStep == 0)
                fileGW << "TimeStep\t" << " time\t" << " min\t" << " max\t"
                       << " l2\t cgNodes\t l2_rs" << std::endl;

            fileGW << m_uiTinfo._m_uiStep << "\t" << m_uiTinfo._m_uiT << "\t"
                   << l_min << "\t" << l_max << "\t" << l2_norm << "\t"
                   << total_dof << "\t " << l_rs << std::endl;
            fileGW.close();
        }
        // initialize diff end
        delete[] chiAnalytical;
        delete[] diffVec;
#endif
    }

    return 0;
}

unsigned int NLSMCtx::getBlkTimestepFac(unsigned int blev, unsigned int lmin,
                                        unsigned int lmax) {
    const unsigned int ldiff = 0;
    if ((lmax - blev) <= ldiff)
        return 1;
    else {
        return 1u << (lmax - blev - ldiff);
    }
}

int NLSMCtx::grid_transfer(const ot::Mesh* m_new) {
#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].start();
#endif

    DVec::grid_transfer(m_uiMesh, m_new, m_evar);

    // resize other unzip variables.
    m_evar_unzip[0].destroy_vector();
    m_evar_unzip[1].destroy_vector();

    m_evar_unzip[0].create_vector(m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                  ot::DVEC_LOC::HOST, NLSM_NUM_VARS, true);
    m_evar_unzip[1].create_vector(m_new, ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,
                                  ot::DVEC_LOC::HOST, NLSM_NUM_VARS, true);

    m_uiIsETSSynced = false;

    ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, NLSM_NUM_VARS,
                                      NLSM_ASYNC_COMM_K);
    ot::alloc_mpi_ctx<DendroScalar>(m_new, m_mpi_ctx, NLSM_NUM_VARS,
                                    NLSM_ASYNC_COMM_K);

#ifdef __PROFILE_CTX__
    m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].stop();
#endif
    return 0;
}

}  // namespace nlsm