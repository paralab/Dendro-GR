//
// Created by milinda on 10/25/18.
//

#include "rhsTestBed.h"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " low_level high_level numBlocks"
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

    if (!rank) {
        printf(
            "=================================================================="
            "=================================================\n");
        printf("low : %d high %d numBlocks %d \n", lLev, hLev, numBlocks);
        printf(
            "=================================================================="
            "=================================================\n");
    }

    std::vector<ot::Block> blkList;
    ot::Block tmpBlk;
    ot::TreeNode rootNode(0, 0, 0, 0, m_uiDim, m_uiMaxDepth);

    const double mean = 0.5 * (hLev - lLev);
    const double sd   = (hLev - mean);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mean, sd);

    unsigned int* blkLevs = new unsigned int[numBlocks];
    unsigned int count    = 0;
    unsigned int unzipSz  = 0;
    while (count < numBlocks) {
        int level = (int)distribution(generator);

        if (level < lLev) level = lLev;

        if (level > hLev) level = (hLev - 1);

        if ((level >= lLev) && (level < hLev)) {
            tmpBlk.setBlkNodeFlag(0);
            tmpBlk = ot::Block(rootNode, 0, level, 0, 0, 4);

            tmpBlk.setOffset(unzipSz);
            blkList.push_back(tmpBlk);
            unzipSz += tmpBlk.getAlignedBlockSz();
            count++;
        }
    }

    const unsigned int UNZIP_DOF = unzipSz;

    // variable input
    double** varUnzipIn          = new double*[bssn::BSSN_NUM_VARS];
    double** varUnzipOutCPU0     = new double*[bssn::BSSN_NUM_VARS];  // staged
    double** varUnzipOutCPU1 = new double*[bssn::BSSN_NUM_VARS];  // unstaged

    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++) {
        varUnzipIn[var]      = new double[UNZIP_DOF];
        varUnzipOutCPU0[var] = new double[UNZIP_DOF];
        varUnzipOutCPU1[var] = new double[UNZIP_DOF];
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

    unsigned int PW = bssn::BSSN_PADDING_WIDTH;

#pragma omp parallel for default(none)                                \
    shared(blkList, pt_min, pt_max, varUnzipIn, bssn::BSSN_COMPD_MAX, \
               bssn::BSSN_COMPD_MIN, bssn::BSSN_OCTREE_MAX,           \
               bssn::BSSN_OCTREE_MIN, PW) schedule(dynamic)
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

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * hx[0];
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * hx[1];
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * hx[2];

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * hx[0];
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * hx[1];
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * hx[2];

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

    Time::time_point t1, t2;
    fsec fs;

    t1 = Time::now();
#pragma omp parallel for default(none) shared(                                 \
        blkList, pt_min, pt_max, varUnzipIn, varUnzipOutCPU0,                  \
            bssn::BSSN_COMPD_MAX, bssn::BSSN_COMPD_MIN, bssn::BSSN_OCTREE_MAX, \
            bssn::BSSN_OCTREE_MIN, PW) schedule(dynamic)
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

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * hx[0];
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * hx[1];
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * hx[2];

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * hx[0];
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * hx[1];
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * hx[2];

        bssnrhs_sep(varUnzipOutCPU0, (const double**)varUnzipIn, offset, ptmin,
                    ptmax, sz, bflag);
    }

    t2 = Time::now();
    fs = t2 - t1;

    bssn::timer::total_runtime.stop();
    std::cout << "CPU compute staged time : " << fs.count() << std::endl;
    std::cout << "Derivative  time : " << bssn::timer::t_deriv.seconds
              << std::endl;
    std::cout << "RHS         time : " << bssn::timer::t_rhs.seconds
              << std::endl;

    t1 = Time::now();

#pragma omp parallel for default(none) shared(                                 \
        blkList, pt_min, pt_max, varUnzipIn, varUnzipOutCPU1,                  \
            bssn::BSSN_COMPD_MAX, bssn::BSSN_COMPD_MIN, bssn::BSSN_OCTREE_MAX, \
            bssn::BSSN_OCTREE_MIN, PW) schedule(dynamic)
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

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * hx[0];
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * hx[1];
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * hx[2];

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * hx[0];
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * hx[1];
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * hx[2];

        bssnrhs(varUnzipOutCPU1, (const double**)varUnzipIn, offset, ptmin,
                ptmax, sz, bflag);
    }
    t2 = Time::now();
    fs = t2 - t1;

    std::cout << "CPU compute unstaged time : " << fs.count() << std::endl;
    std::cout << "Derivative  time : " << bssn::timer::t_deriv.seconds
              << std::endl;
    std::cout << "RHS         time : " << bssn::timer::t_rhs.seconds
              << std::endl;

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
                            fabs(varUnzipOutCPU0[var]
                                                [offset + k * sz[0] * sz[1] +
                                                 j * sz[0] + i] -
                                 varUnzipOutCPU1[var]
                                                [offset + k * sz[0] * sz[1] +
                                                 j * sz[0] + i]))
                            l_inf =
                                fabs(varUnzipOutCPU0[var][offset +
                                                          k * sz[0] * sz[1] +
                                                          j * sz[0] + i] -
                                     varUnzipOutCPU1[var][offset +
                                                          k * sz[0] * sz[1] +
                                                          j * sz[0] + i]);
        }

        std::cout << "comparison for var: " << var << bssn::BSSN_VAR_NAMES[var]
                  << " l_inf : " << l_inf << std::endl;
    }

    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++) {
        delete[] varUnzipIn[var];
        delete[] varUnzipOutCPU0[var];
        delete[] varUnzipOutCPU1[var];
    }

    delete[] varUnzipIn;
    delete[] varUnzipOutCPU1;
    delete[] varUnzipOutCPU0;

    MPI_Finalize();

    return 0;
}
