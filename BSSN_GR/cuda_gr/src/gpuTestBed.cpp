//
// Created by milinda on 9/19/18.
//

#include "gpuTestBed.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0]
                  << " low_level high_level numBlocks(for each level) dim "
                     "threadX=2,threadY=2,threadZ=2 numstreams useStaged(0 "
                     "unstage 1 staged) dist (0-uniform 1-gauss)"
                  << std::endl;
        return 0;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    m_uiMaxDepth                 = 8;

    const unsigned int lLev      = atoi(argv[1]);
    const unsigned int hLev      = atoi(argv[2]);

    const unsigned int numBlocks = atoi(argv[3]);

    const unsigned int dim       = atoi(argv[4]);

    unsigned int threadX         = 2;
    unsigned int threadY         = 2;
    unsigned int threadZ         = 2;
    unsigned int numStreams      = 1;
    unsigned int useStaged       = 1;
    unsigned int dist            = 1;

    if (argc > 5) threadX = atoi(argv[5]);
    if (argc > 6) threadY = atoi(argv[6]);
    if (argc > 7) threadZ = atoi(argv[7]);

    bool useAsync = false;

    if (argc > 8) {
        numStreams = atoi(argv[8]);
        useAsync   = true;
    }

    if (argc > 9) useStaged = atoi(argv[9]);

    if (argc > 10) dist = atoi(argv[10]);

    if (!rank) {
        printf(
            "=================================================================="
            "=================================================\n");
        printf(
            "low : %d high %d numBlocks %d dim %d tx %d ty %d tz %d numstreams "
            "%d staged %d  distribution %d \n",
            lLev, hLev, numBlocks, dim, threadX, threadY, threadZ, numStreams,
            useStaged, dist);
        printf(
            "=================================================================="
            "=================================================\n");
    }

    std::vector<ot::Block> blkList;
    ot::Block tmpBlk;
    ot::TreeNode rootNode(0, 0, 0, 0, m_uiDim, m_uiMaxDepth);
    unsigned long unzipSz = 0;
    if (dist == 0) {
        for (unsigned int lev = lLev; lev < hLev; lev++) {
            for (unsigned int i = 0; i < numBlocks; i++) {
                tmpBlk.setBlkNodeFlag(0);

                if (dim == 3) {
                    tmpBlk = ot::Block(rootNode, 0, lev, 0, 0, 4);

                } else if (dim == 1) {
                    tmpBlk.setAllocationSzX((1u << (lev)));
                    tmpBlk.setAllocationSzY(1);
                    tmpBlk.setAllocationSzZ(1);
                    tmpBlk.setSiz1D((1u << (lev)));
                }

                tmpBlk.setOffset(unzipSz);
                blkList.push_back(tmpBlk);
                unzipSz += tmpBlk.getAlignedBlockSz();
            }
        }
    } else if (dist == 1) {
        const double mean = 0.5 * (hLev - lLev);
        const double sd   = (hLev - mean);
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(mean, sd);

        unsigned int *blkLevs = new unsigned int[numBlocks];
        unsigned int count    = 0;
        while (count < numBlocks) {
            int level = (int)distribution(generator);

            if (level < lLev) level = lLev;

            if (level > hLev) level = (hLev - 1);

            if ((level >= lLev) && (level < hLev)) {
                tmpBlk.setBlkNodeFlag(0);
                if (dim == 3) {
                    tmpBlk = ot::Block(rootNode, 0, level, 0, 0, 4);

                } else if (dim == 1) {
                    tmpBlk.setAllocationSzX((1u << (level)));
                    tmpBlk.setAllocationSzY(1);
                    tmpBlk.setAllocationSzZ(1);
                    tmpBlk.setSiz1D((1u << (level)));
                }

                tmpBlk.setOffset(unzipSz);
                blkList.push_back(tmpBlk);
                unzipSz += tmpBlk.getAlignedBlockSz();
                count++;
            }
        }
    }

    const unsigned int UNZIP_DOF = unzipSz;

    // variable input
    double **varUnzipIn          = new double *[bssn::BSSN_NUM_VARS];
    double **varUnzipOutGPU      = new double *[bssn::BSSN_NUM_VARS];
    double **varUnzipOutCPU      = new double *[bssn::BSSN_NUM_VARS];

    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++) {
        varUnzipIn[var]     = new double[UNZIP_DOF];
        varUnzipOutGPU[var] = new double[UNZIP_DOF];
        varUnzipOutCPU[var] = new double[UNZIP_DOF];
    }

    bssn::BSSN_COMPD_MIN[0]  = -1e-3;
    bssn::BSSN_COMPD_MIN[1]  = -1e-3;
    bssn::BSSN_COMPD_MIN[2]  = -1e-3;

    bssn::BSSN_COMPD_MAX[0]  = 1e-3;
    bssn::BSSN_COMPD_MAX[1]  = 1e-3;
    bssn::BSSN_COMPD_MAX[2]  = 1e-3;

    bssn::BSSN_OCTREE_MIN[0] = 0;
    bssn::BSSN_OCTREE_MIN[1] = 0;
    bssn::BSSN_OCTREE_MIN[2] = 0;

    bssn::BSSN_OCTREE_MAX[0] = 1u << m_uiMaxDepth;
    bssn::BSSN_OCTREE_MAX[1] = 1u << m_uiMaxDepth;
    bssn::BSSN_OCTREE_MAX[2] = 1u << m_uiMaxDepth;

    bssn::KO_DISS_SIGMA      = 0;
    bssn::BH1 = bssn::BH(0.48, -0.01, 0, 1e-6, 0.1, 0, 0, 0, 0, 0);
    bssn::BH2 = bssn::BH(0.52, 0.01, 0, 1e-6, -0.1, 0, 0, 0, 0, 0);

    Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                 bssn::BSSN_COMPD_MIN[2]);
    Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                 bssn::BSSN_COMPD_MAX[2]);

#pragma omp parallel for default(none)                          \
    shared(blkList, pt_min, pt_max, varUnzipIn, varUnzipOutCPU, \
               bssn::BSSN_COMPD_MAX, bssn::BSSN_COMPD_MIN,      \
               bssn::BSSN_OCTREE_MAX, bssn::BSSN_OCTREE_MIN) schedule(dynamic)
    for (unsigned int blk = 0; blk < blkList.size(); blk++) {
        const unsigned int offset = blkList[blk].getOffset();
        unsigned int sz[3];
        double hx[3];
        double ptmin[3];
        double ptmax[3];

        double x, y, z;
        double varVal[bssn::BSSN_NUM_VARS];

        sz[0]                    = blkList[blk].getAllocationSzX();
        sz[1]                    = blkList[blk].getAllocationSzY();
        sz[2]                    = blkList[blk].getAllocationSzZ();

        const unsigned int bflag = blkList[blk].getBlkNodeFlag();

        hx[0]                    = blkList[blk].computeDx(pt_min, pt_max);
        hx[1]                    = blkList[blk].computeDy(pt_min, pt_max);
        hx[2]                    = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - 3 * hx[0];
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - 3 * hx[1];
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - 3 * hx[2];

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + 3 * hx[0];
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + 3 * hx[1];
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + 3 * hx[2];

        /*printf("grid x min to x : %f grid x max to x: %f
           \n",GRIDX_TO_X(blkList[blk].getBlockNode().minX()),GRIDX_TO_X(blkList[blk].getBlockNode().maxX()));
            printf(" x min to grid x : %f x max to grid x: %f
           \n",X_TO_GRIDX(GRIDX_TO_X(blkList[blk].getBlockNode().minX())),X_TO_GRIDX(GRIDX_TO_X(blkList[blk].getBlockNode().maxX())));
            printf("dx: %f \n",dx);*/

        for (unsigned int k = 0; k < sz[2]; k++) {
            z = ptmin[2] + k * hx[2];
            for (unsigned int j = 0; j < sz[1]; j++) {
                y = ptmin[1] + j * hx[1];
                for (unsigned int i = 0; i < sz[0]; i++) {
                    x = ptmin[0] + i * hx[0];

                    // bssn::punctureData(x,y,z,varVal);
                    // printf("x %f y %f z
                    // %f\n",X_TO_GRIDX(x),Y_TO_GRIDY(y),Z_TO_GRIDZ(z));
                    bssn::fake_initial_data(X_TO_GRIDX(x), Y_TO_GRIDY(y),
                                            Z_TO_GRIDZ(z), varVal);

                    for (unsigned int v = 0; v < bssn::BSSN_NUM_VARS; v++)
                        varUnzipIn[v][offset + k * sz[0] * sz[1] + j * sz[0] +
                                      i] = varVal[v];
                }
            }
        }
    }

    if (useStaged == 1) {
        bssn::timer::initFlops();

        auto t1 = Time::now();

#pragma omp parallel for default(none)                          \
    shared(blkList, pt_min, pt_max, varUnzipIn, varUnzipOutCPU, \
               bssn::BSSN_COMPD_MAX, bssn::BSSN_COMPD_MIN,      \
               bssn::BSSN_OCTREE_MAX, bssn::BSSN_OCTREE_MIN) schedule(dynamic)
        for (unsigned int blk = 0; blk < blkList.size(); blk++) {
            const unsigned int offset = blkList[blk].getOffset();
            unsigned int sz[3];
            double hx[3];
            double ptmin[3];
            double ptmax[3];

            sz[0]                    = blkList[blk].getAllocationSzX();
            sz[1]                    = blkList[blk].getAllocationSzY();
            sz[2]                    = blkList[blk].getAllocationSzZ();

            const unsigned int bflag = blkList[blk].getBlkNodeFlag();

            hx[0]                    = blkList[blk].computeDx(pt_min, pt_max);
            hx[1]                    = blkList[blk].computeDy(pt_min, pt_max);
            hx[2]                    = blkList[blk].computeDz(pt_min, pt_max);

            ptmin[0] =
                GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - 3 * hx[0];
            ptmin[1] =
                GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - 3 * hx[1];
            ptmin[2] =
                GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - 3 * hx[2];

            ptmax[0] =
                GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + 3 * hx[0];
            ptmax[1] =
                GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + 3 * hx[1];
            ptmax[2] =
                GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + 3 * hx[2];

            bssnrhs(varUnzipOutCPU, (const double **)varUnzipIn, offset, ptmin,
                    ptmax, sz, bflag);
        }

        auto t2 = Time::now();
        fsec fs = t2 - t1;

        bssn::timer::total_runtime.stop();
        std::cout << "CPU compute staged time : " << fs.count() << std::endl;
        std::cout << "Derivative  time : " << bssn::timer::t_deriv.seconds
                  << std::endl;
        std::cout << "RHS         time : " << bssn::timer::t_rhs.seconds
                  << std::endl;

    } else if (useStaged == 0) {
        bssn::timer::initFlops();

        auto t1 = Time::now();

#pragma omp parallel for default(none)                          \
    shared(blkList, pt_min, pt_max, varUnzipIn, varUnzipOutCPU, \
               bssn::BSSN_COMPD_MAX, bssn::BSSN_COMPD_MIN,      \
               bssn::BSSN_OCTREE_MAX, bssn::BSSN_OCTREE_MIN) schedule(dynamic)
        for (unsigned int blk = 0; blk < blkList.size(); blk++) {
            const unsigned int offset = blkList[blk].getOffset();
            unsigned int sz[3];
            double hx[3];
            double ptmin[3];
            double ptmax[3];

            sz[0]                    = blkList[blk].getAllocationSzX();
            sz[1]                    = blkList[blk].getAllocationSzY();
            sz[2]                    = blkList[blk].getAllocationSzZ();

            const unsigned int bflag = blkList[blk].getBlkNodeFlag();

            hx[0]                    = blkList[blk].computeDx(pt_min, pt_max);
            hx[1]                    = blkList[blk].computeDy(pt_min, pt_max);
            hx[2]                    = blkList[blk].computeDz(pt_min, pt_max);

            ptmin[0] =
                GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - 3 * hx[0];
            ptmin[1] =
                GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - 3 * hx[1];
            ptmin[2] =
                GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - 3 * hx[2];

            ptmax[0] =
                GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + 3 * hx[0];
            ptmax[1] =
                GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + 3 * hx[1];
            ptmax[2] =
                GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + 3 * hx[2];

            bssnrhs(varUnzipOutCPU, (const double **)varUnzipIn, offset, ptmin,
                    ptmax, sz, bflag);
        }
        auto t2 = Time::now();
        fsec fs = t2 - t1;

        std::cout << "CPU compute unstaged time : " << fs.count() << std::endl;
        std::cout << "Derivative  time : " << bssn::timer::t_deriv.seconds
                  << std::endl;
        std::cout << "RHS         time : " << bssn::timer::t_rhs.seconds
                  << std::endl;
    }

#ifdef BSSN_ENABLE_CUDA
    cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
    CUDA_CHECK_ERROR();
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
    CUDA_CHECK_ERROR();
    cuda::BSSNComputeParams bssnParams;
    bssnParams.BSSN_LAMBDA[0]    = bssn::BSSN_LAMBDA[0];
    bssnParams.BSSN_LAMBDA[1]    = bssn::BSSN_LAMBDA[1];
    bssnParams.BSSN_LAMBDA[2]    = bssn::BSSN_LAMBDA[2];
    bssnParams.BSSN_LAMBDA[3]    = bssn::BSSN_LAMBDA[3];

    bssnParams.BSSN_LAMBDA_F[0]  = bssn::BSSN_LAMBDA_F[0];
    bssnParams.BSSN_LAMBDA_F[1]  = bssn::BSSN_LAMBDA_F[1];

    bssnParams.BSSN_ETA_POWER[0] = bssn::BSSN_ETA_POWER[0];
    bssnParams.BSSN_ETA_POWER[1] = bssn::BSSN_ETA_POWER[1];

    bssnParams.ETA_R0            = bssn::ETA_R0;
    // printf("R0 %f \n",bssnParams.BSSN_ETA_R0);
    bssnParams.ETA_CONST         = bssn::ETA_CONST;
    bssnParams.ETA_DAMPING       = bssn::ETA_DAMPING;
    bssnParams.ETA_DAMPING_EXP   = bssn::ETA_DAMPING_EXP;
    bssnParams.KO_DISS_SIGMA     = bssn::KO_DISS_SIGMA;

    dim3 threadBlock(threadX, threadY, threadZ);

    cuda::profile::initialize();
    cuda::computeRHS(varUnzipOutGPU, (const double **)varUnzipIn,
                     &(*(blkList.begin())), blkList.size(),
                     (const cuda::BSSNComputeParams *)&bssnParams, threadBlock,
                     pt_min, pt_max, numStreams, rank);

    cuda::profile::printOutput(blkList);

    std::cout << YLW << " ================================" << NRM << std::endl;
    std::cout << YLW << "     GPU end             " << NRM << std::endl;
    std::cout << YLW << " ================================" << NRM << std::endl;

#endif

    double l_inf;
    for (unsigned int var = 0; var < 24 /*bssn::BSSN_NUM_VARS*/; var++) {
        l_inf = 0;
        for (unsigned int blk = 0; blk < blkList.size(); blk++) {
            const unsigned int offset = blkList[blk].getOffset();
            unsigned int sz[3];

            sz[0] = blkList[blk].getAllocationSzX();
            sz[1] = blkList[blk].getAllocationSzY();
            sz[2] = blkList[blk].getAllocationSzZ();

            for (unsigned int k = 3; k < sz[2] - 3; k++)
                for (unsigned int j = 3; j < sz[1] - 3; j++)
                    for (unsigned int i = 3; i < sz[0] - 3; i++)
                        if (l_inf <
                            fabs(
                                varUnzipOutCPU[var][offset + k * sz[0] * sz[1] +
                                                    j * sz[0] + i] -
                                varUnzipOutGPU
                                    [var]
                                    [offset + k * sz[0] * sz[1] + j * sz[0] +
                                     i]) /*varUnzipOutGPU[var][offset+k*sz[0]*sz[1]+j*sz[0]+i]!=1*/)
                            l_inf = fabs(
                                varUnzipOutCPU[var][offset + k * sz[0] * sz[1] +
                                                    j * sz[0] + i] -
                                varUnzipOutGPU[var][offset + k * sz[0] * sz[1] +
                                                    j * sz[0] + i]);
            // std::cout<<"blk: "<<blk<<" offset : "<<offset<<"i,j,k:
            // ("<<i<<","<<j<<", "<<k<<")"<<" cpu:
            // "<<varUnzipOutCPU[var][offset+k*sz[0]*sz[1]+j*sz[0]+i]<<" gpu :
            // "<<varUnzipOutGPU[var][offset+k*sz[0]*sz[1]+j*sz[0]+i]<<" diff :
            // "<<fabs(varUnzipOutCPU[var][offset+k*sz[0]*sz[1]+j*sz[0]+i]-varUnzipOutGPU[var][offset+k*sz[0]*sz[1]+j*sz[0]+i])<<std::endl;
        }

        std::cout << "comparison for var: " << var << bssn::BSSN_VAR_NAMES[var]
                  << " l_inf : " << l_inf << std::endl;
    }

    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++) {
        delete[] varUnzipIn[var];
        delete[] varUnzipOutGPU[var];
        delete[] varUnzipOutCPU[var];
    }

    delete[] varUnzipIn;
    delete[] varUnzipOutCPU;
    delete[] varUnzipOutGPU;

    MPI_Finalize();

    return 0;
}
