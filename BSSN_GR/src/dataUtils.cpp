//
// Created by milinda on 1/16/19.
//

#include "dataUtils.h"

#include "parameters.h"

namespace bssn {

void extractBHCoords(const ot::Mesh* pMesh, const DendroScalar* var,
                     double tolerance, const Point* ptIn, unsigned int numPt,
                     Point* ptOut) {
    if ((pMesh->isActive())) {
        MPI_Comm commActive     = pMesh->getMPICommunicator();
        unsigned int rankActive = pMesh->getMPIRank();
        unsigned int npesActive = pMesh->getMPICommSize();

        double v_min = vecMin((DendroScalar*)(var + pMesh->getNodeLocalBegin()),
                              pMesh->getNumLocalMeshNodes(), commActive);
        par::Mpi_Bcast(&v_min, 1, 0, commActive);

        assert(numPt == 2);

        const double extraction_tol     = tolerance;  // 10*v_min;

        const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
        const unsigned int* e2n         = &(*(pMesh->getE2NMapping().begin()));
        const unsigned int* e2n_dg = &(*(pMesh->getE2NMapping_DG().begin()));
        const unsigned int* cgToDg = &(*(pMesh->getCG2DGMap().begin()));

        unsigned int lookup        = 0;
        unsigned int ownerID, ii_x, jj_y, kk_z;

        ot::TreeNode tmpOct;
        const unsigned int eleOrder = pMesh->getElementOrder();
        double hx, x, y, z;

        std::vector<Point> ptList;
        for (unsigned int node = pMesh->getNodeLocalBegin();
             node < pMesh->getNodeLocalEnd(); node++) {
            if (var[node] < extraction_tol) {
                lookup = cgToDg[node];
                pMesh->dg2eijk(lookup, ownerID, ii_x, jj_y, kk_z);
                tmpOct = allElements[ownerID];
                hx     = (tmpOct.maxX() - tmpOct.minX()) / ((double)eleOrder);
                x      = tmpOct.minX() + ii_x * (hx);
                y      = tmpOct.minY() + jj_y * (hx);
                z      = tmpOct.minZ() + kk_z * (hx);

                ptList.push_back(
                    Point(GRIDX_TO_X(x), GRIDY_TO_Y(y), GRIDZ_TO_Z(z)));

                // std::cout<<" x: "<<GRIDX_TO_X(x)<<" y: "<<GRIDY_TO_Y(y)<<" z:
                // "<<GRIDZ_TO_Z(z)<<std::endl;
            }
        }

        std::vector<Point>* ptCluster = new std::vector<Point>[numPt];

        double min0, min1;
        unsigned int cID;
        for (unsigned int pt = 0; pt < ptList.size(); pt++) {
            cID  = 0;
            min0 = (ptIn[0] - ptList[pt]).abs();
            min1 = (ptIn[1] - ptList[pt]).abs();

            if (fabs(min0 - min1) <
                1e-6) {  // implies the point is closer to the both clusters
                ptCluster[0].push_back(ptList[pt]);
                ptCluster[1].push_back(ptList[pt]);
            } else if (min0 < min1) {
                ptCluster[0].push_back(ptList[pt]);
            } else {
                ptCluster[1].push_back(ptList[pt]);
            }
        }

        ptList.clear();
        Point* ptMean          = new Point[numPt];
        DendroIntL* ptCounts   = new DendroIntL[numPt];
        DendroIntL* ptCounts_g = new DendroIntL[numPt];

        for (unsigned int c = 0; c < numPt; c++) {
            ptMean[c] = Point(0, 0, 0);
            ptOut[c]  = Point(0, 0, 0);
        }

        for (unsigned int c = 0; c < numPt; c++) {
            ptCounts[c] = ptCluster[c].size();
            for (unsigned int pt = 0; pt < ptCluster[c].size(); pt++)
                ptMean[c] += ptCluster[c][pt];
        }

        par::Mpi_Allreduce(ptCounts, ptCounts_g, numPt, MPI_SUM, commActive);
        par::Mpi_Allreduce(ptMean, ptOut, numPt,
                           par::Mpi_datatype<Point>::_SUM(), commActive);

        for (unsigned int c = 0; c < numPt; c++)
            ptOut[c] /= (double)ptCounts_g[c];

        // if(pMesh->getMPIRank()==0)std::cout<<"bh1 in : "<<ptIn[0].x()<<",
        // "<<ptIn[0].y()<<", "<<ptIn[0].z()<<std::endl;
        // if(pMesh->getMPIRank()==0)std::cout<<"bh2 in : "<<ptIn[1].x()<<",
        // "<<ptIn[1].y()<<", "<<ptIn[1].z()<<std::endl;

        // if(pMesh->getMPIRank()==0)std::cout<<"bh1 out: "<<ptOut[0].x()<<",
        // "<<ptOut[0].y()<<", "<<ptOut[0].z()<<std::endl;
        // if(pMesh->getMPIRank()==0)std::cout<<"bh2 out: "<<ptOut[1].x()<<",
        // "<<ptOut[1].y()<<", "<<ptOut[1].z()<<std::endl;

        delete[] ptCounts;
        delete[] ptCounts_g;
        delete[] ptMean;
        delete[] ptCluster;
    }
}

void writeBHCoordinates(const ot::Mesh* pMesh, const Point* ptLocs,
                        unsigned int numPt, unsigned int timestep,
                        double time) {
    unsigned int rankGlobal = pMesh->getMPIRankGlobal();
    if (!rankGlobal) {
        std::ofstream fileGW;
        char fName[256];
        sprintf(fName, "%s_BHLocations.dat",
                bssn::BSSN_PROFILE_FILE_PREFIX.c_str());
        fileGW.open(fName, std::ofstream::app);

        // writes the header
        if (timestep == 0)
            fileGW << "TimeStep\t" << " time\t" << " bh1_x\t" << " bh1_y\t"
                   << " bh1_z\t" << " bh2_x\t" << " bh2_y\t" << " bh2_z\t"
                   << std::endl;

        fileGW << timestep << "\t" << time << "\t" << ptLocs[0].x() << "\t"
               << ptLocs[0].y() << "\t" << ptLocs[0].z() << "\t"
               << ptLocs[1].x() << "\t" << ptLocs[1].y() << "\t"
               << ptLocs[1].z() << std::endl;
        fileGW.close();
        return;
    }
}

bool isRemeshBH(ot::Mesh* pMesh, const Point* bhLoc) {
    const double r_near[2]       = {bssn::BSSN_BH1_AMR_R, bssn::BSSN_BH2_AMR_R};
    // const double r_far[2]  =  {2.5 * r_near[0], 2.5 * r_near[1] };
    const double r_far[2]        = {bssn::BSSN_AMR_R_RATIO * r_near[0],
                                    bssn::BSSN_AMR_R_RATIO * r_near[1]};
    // set up level offsets near black holes
    // (we don't actually refine to BSSN_BH?_MAX_LEV)
    // level offset immediately about the BHs
    const unsigned int LVL_OFF = MAXDEAPTH_LEVEL_DIFF + 1;

    // consider BHs merged if punctures are less than this value
    const double BH_MERGED_SEP_TOL      = 0.1;

    const unsigned int eleLocalBegin    = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd      = pMesh->getElementLocalEnd();
    bool isOctChange                    = false;
    bool isOctChange_g                  = false;
    Point d1, d2, temp;
    // distance btw the black holes
    const double dBH = (bhLoc[0] - bhLoc[1]).abs();
    // lower of two max depths
    const unsigned int refLevMin =
        std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);

    std::vector<unsigned int> refine_flags;

    // NOTE: these are the previous refine flags,
    // they're usually "no change" on remesh
    std::vector<unsigned int> prev_refine_flags =
        pMesh->getAllRefinementFlags();

    if (pMesh->isActive()) {
        // if(!pMesh->getMPIRank())
        //     std::cout<<"bh distance: "<<dBH<<std::endl;

        const ot::TreeNode* pNodes = pMesh->getAllElements().data();
        refine_flags.resize(pMesh->getNumLocalMeshElements(), OCT_NO_CHANGE);

        // refine pass: iterate over all elements
        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++) {
            // calculate which region of the grid we're in
            const unsigned int ln = 1u
                                    << (m_uiMaxDepth - pNodes[ele].getLevel());

            // Set obnoxiously large values for minimum radii
            // (to be overwritten) so as we go through box edges
            // we find minimum distances to BHs and to grid center
            double r1_min = 1000000;
            double r2_min = 1000000;
            double r_min  = 1000000; // distance to grid center

            // measure minimum distances between both BHs
            for (unsigned int kk = 0; kk < 2; kk++)
                for (unsigned int jj = 0; jj < 2; jj++)
                    for (unsigned int ii = 0; ii < 2; ii++) {
                        const double x      = pNodes[ele].minX() + ii * ln;
                        const double y      = pNodes[ele].minY() + jj * ln;
                        const double z      = pNodes[ele].minZ() + kk * ln;
                        const Point oct_mid = Point(x, y, z);
                        pMesh->octCoordToDomainCoord(oct_mid, temp);

                        // vectors pointing toward each BH
                        d1     = temp - bhLoc[0];
                        d2     = temp - bhLoc[1];
                        // update minimum distance from each BH
                        r1_min = std::min(r1_min, d1.abs());
                        r2_min = std::min(r2_min, d2.abs());
                        // update minimum distance from grid center
                        r_min  = std::min(r_min, temp.abs());
                    }

            ////////////////////////////////////////////////////////////
            // wkb 5 Sept 2024: make this into nice functions

            // set default of coarsening; use BHLB before any others!
            refine_flags[ele - eleLocalBegin] = OCT_COARSE;

            // function which ensures we're at least at a given level
            auto setLevelFloor                = [&](int l_min) {
                int currentLevel = pNodes[ele].getLevel();

                if (currentLevel < l_min) {
                    // if below desired refinement level, split
                    refine_flags[ele - eleLocalBegin] = OCT_SPLIT;
                } else if (currentLevel == l_min &&
                           refine_flags[ele - eleLocalBegin] == OCT_COARSE) {
                    // if at desired level, prevent coarsening
                    refine_flags[ele - eleLocalBegin] = OCT_NO_CHANGE;
                }
                // else: sufficiently refined
            };

            ////////////////////////////////////////////////////////////
            // wkb 29 Aug 2024: add refinement based on radius
            // if radius r <= r*, keep level l >= l*
            // 9/9/24: change to depend on current BH separation dist.

            // set up orbital scale
            const double m1      = bssn::BSSN_BH1_MASS;  // mass of BH1
            const double m2      = bssn::BSSN_BH2_MASS;  // mass of BH2
            // calculate relative distance to each BH
            const double f1      = m2 / (m1 + m2);
            const double f2      = m1 / (m1 + m2);
            const double f       = std::max(f1, f2);
            const double R_orbit = f * dBH + 8;  // M; resolve scale
            const int l_orbit    = 9;  // desired refinement level within
            if (r_min <= R_orbit) {
                // set up orbital radius scale
                setLevelFloor(l_orbit);
            }

            ////////////////////////////////////////////////////////////
            // wkb 2 Dec 2024: Onion refinement about the BHs 
            auto onionLevel = [R_orbit, l_orbit](double radius, double r_AMR, int maxLevel, double ratio = 2.0) -> int {
                /*
                // TEMPORARY: disable all but innermost level
                if (radius < r_AMR) {
                    return maxLevel;
                } else {
                    return 0;
                }
                */
                if (radius > R_orbit) {
                    // don't enforce onion outside orbital radius
                    return 0;
                } else {
                    // Start with the AMR radius and maximum level
                    double currentRadius = r_AMR;
                    int currentLevel = maxLevel;
                    // Loop until we reach or go below the orbit level
                    while (currentLevel > l_orbit) {
                        // If input radius w/i current radius,
                        // return current level requirement
                        if (radius <= currentRadius) {
                            return currentLevel;
                        }
                        // Otherwise increase the radius limit
                        currentRadius *= ratio;
                        // and decrement the refinement level
                        currentLevel--;
                    }
                    // Shouldn't be here.
                    return -1;
                }
            };

            // Onion refinement immediately about the black holes
            if (dBH > BH_MERGED_SEP_TOL) {
                // if not merged yet, handle BHs separately to set up onion 
                double ratio;
                int shift; 
                if (bssn::BSSN_CURRENT_RK_COORD_TIME < -10) {
                    // boost initial refinement, bolstering onion
                    ratio = 2.0;
                    shift = 1;
                } else { // relax onion later
                    // ratio = 1.618033988749; // golden ratio
                    ratio = bssn::BSSN_AMR_R_RATIO; // softcode
                    shift = 0;
                }
                const int l_goal_0 = onionLevel(r1_min,r_near[0],bssn::BSSN_BH1_MAX_LEV - LVL_OFF + shift,ratio);
                const int l_goal_1 = onionLevel(r2_min,r_near[1],bssn::BSSN_BH2_MAX_LEV - LVL_OFF + shift,ratio);
                setLevelFloor(l_goal_0);
                setLevelFloor(l_goal_1);
            } else {
                // if merged, handle BHs together
                // calculate minimum distance to either BH
                const double rBH_min = std::min(r1_min,r2_min);
                // calculate outer radius to which we should refine
                // ensure it captures both BHs and is >= than before
                const double rBH_lim = std::max(std::max(r_near[0],r_near[1]),1.0 * (m1 + m2)); 
                // calculate level floor due to onion structure
                const int l_post = 14; // hardcoding innermost level 
                const int l_goal = onionLevel(rBH_min,rBH_lim,l_post - LVL_OFF);
                // const int l_goal = onionLevel(rBH_min,rBH_lim,refLevMin - LVL_OFF);
                // set level floor
                setLevelFloor(l_goal);
            }
            
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
            ////////////////////////////////////////////////////////////
            // @wkb 4 Sept 2024:
            // add refinement based on expected gravitational wavelength
            const double eta     = f1 * f2;  // symmetric mass ratio
            // estimate ending GW frequency from quasi-normal modes 
            const double lam_qnm = 4 * M_PI / (.37009 + .6475 * eta);
            // clang-format off
            #if 1
            // q=1 parameters
            constexpr double A    = 18.8;
            // More nuanced fit: A tau^(3/8) * (1 + B tau^(-2/8))
            // constexpr double B    = 0.0; 
            constexpr double tau0 = 445.0;
            // constexpr double mn = 10.5; // old
            constexpr double mn = 14.49;
            // constexpr double t_end = 642.0; // time (2,2) hits R=50
            // constexpr double t_end = 608.14; // ~150 past merger
            constexpr double t_end = 580.14; // ~122 past merger
            #else
            // q=4 parameters
            constexpr double A    = 17;
            // More nuanced fit: A tau^(3/8) * (1 + B tau^(-2/8))
            // constexpr double B    = .75;
            constexpr double tau0 = 490.0;
            // constexpr double mn = 22.90; // old
            constexpr double mn = 18.10;
            // constexpr double t_end = 660.0;
            constexpr double t_end = 670.4; // ~122 past merger
            #endif
            // decide whether to use QNM or orbital freq for ending
            const double lam_min = mn; // use orbital freq @ ds = .1
            // clang-format on

            auto get_ell = [A, tau0, lam_min](double t_ret,
                                              unsigned int m = 8) -> int {
                // Get refinement level necessary from wavelength
                double lambda = (t_ret < tau0)
                                    ? A * std::pow(tau0 - t_ret, 3.0 / 8.0)
                                    : lam_min;
                lambda        = std::max(lambda, lam_min);
                return static_cast<int>(
                    std::ceil(std::log2(800.0 * m / (3.0 * lambda))));
            };

            // calculate retarded time
            const double t_ret = bssn::BSSN_CURRENT_RK_COORD_TIME - r_min;
            // const double R_GW  = 50.0; // inner GW radius; TODO: should use outer probably
            const double R_GW  = 100.0; // outer GW radius
            int ell_star;

            /*
            // min refinement level required from GWs
            // use this code here to have different criteria beyond
            // the GW extraction radii to cut down costs
            if (r_min <= R_GW) {
              // if w/i region of capturing GWs, resolve highly
              ell_star = get_ell(t_ret, bssn::BSSN_NYQUIST_M);
            } else {
              // otherwise, only prevent major backreflections
              ell_star = get_ell(t_ret, bssn::BSSN_NYQUIST_M - 2);
            }
            */

            // goal spherical harmonic order m to refine to
            // DFVK NOTE: put the variable is now in the global parameter
            const bool using_nyquist = bssn::BSSN_NYQUIST_M > 0;
            const double speed = std::sqrt(2);
            const bool past_of_end = std::abs(r_min - R_GW) < speed * (t_end - bssn::BSSN_CURRENT_RK_COORD_TIME);
            if (using_nyquist && past_of_end) {
                // same ell_star everywhere
                ell_star = get_ell(t_ret, bssn::BSSN_NYQUIST_M);
                setLevelFloor(ell_star);
            }
#endif

            ////////////////////////////////////////////////////////////
            // @wkb 12 June 2024: Don't allow min depth violation
            setLevelFloor(bssn::BSSN_MINDEPTH);
        }

        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    bool isOctChanged_g;
    MPI_Allreduce(&isOctChange, &isOctChanged_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChanged_g;
}

bool isRemeshEH(ot::Mesh* pMesh, const double** unzipVec, unsigned int vIndex,
                double refine_th, double coarsen_th, bool isOverwrite) {
    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
    bool isOctChange                 = false;
    bool isOctChange_g               = false;
    const unsigned int eOrder        = pMesh->getElementOrder();

    std::vector<unsigned int> refine_flags;

    // NOTE: this is what the flags are set to, note that
    const std::vector<unsigned int> prev_refine_flags =
        pMesh->getAllRefinementFlags();

    if (pMesh->isActive()) {
        ot::TreeNode* pNodes =
            (ot::TreeNode*)&(*(pMesh->getAllElements().begin()));

        if (isOverwrite)
            for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++)
                pNodes[ele].setFlag(((OCT_NO_CHANGE << NUM_LEVEL_BITS) |
                                     pNodes[ele].getLevel()));

        const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
        unsigned int sz[3];
        unsigned int ei[3];

        // refine test
        for (unsigned int b = 0; b < blkList.size(); b++) {
            const ot::TreeNode blkNode = blkList[b].getBlockNode();

            sz[0]                      = blkList[b].getAllocationSzX();
            sz[1]                      = blkList[b].getAllocationSzY();
            sz[2]                      = blkList[b].getAllocationSzZ();

            const unsigned int bflag   = blkList[b].getBlkNodeFlag();
            const unsigned int offset  = blkList[b].getOffset();

            const unsigned int regLev  = blkList[b].getRegularGridLev();
            const unsigned int eleIndexMax =
                (1u << (regLev - blkNode.getLevel())) - 1;
            const unsigned int eleIndexMin = 0;

            for (unsigned int ele = blkList[b].getLocalElementBegin();
                 ele < blkList[b].getLocalElementEnd(); ele++) {
                ei[0] = (pNodes[ele].getX() - blkNode.getX()) >>
                        (m_uiMaxDepth - regLev);
                ei[1] = (pNodes[ele].getY() - blkNode.getY()) >>
                        (m_uiMaxDepth - regLev);
                ei[2] = (pNodes[ele].getZ() - blkNode.getZ()) >>
                        (m_uiMaxDepth - regLev);

                if ((bflag & (1u << OCT_DIR_LEFT)) && ei[0] == eleIndexMin)
                    continue;
                if ((bflag & (1u << OCT_DIR_DOWN)) && ei[1] == eleIndexMin)
                    continue;
                if ((bflag & (1u << OCT_DIR_BACK)) && ei[2] == eleIndexMin)
                    continue;

                if ((bflag & (1u << OCT_DIR_RIGHT)) && ei[0] == eleIndexMax)
                    continue;
                if ((bflag & (1u << OCT_DIR_UP)) && ei[1] == eleIndexMax)
                    continue;
                if ((bflag & (1u << OCT_DIR_FRONT)) && ei[2] == eleIndexMax)
                    continue;

                // refine test.
                for (unsigned int k = 3; k < eOrder + 1 + 3; k++)
                    for (unsigned int j = 3; j < eOrder + 1 + 3; j++)
                        for (unsigned int i = 3; i < eOrder + 1 + 3; i++) {
                            if (unzipVec[vIndex]
                                        [offset +
                                         (ei[2] * eOrder + k) * sz[0] * sz[1] +
                                         (ei[1] * eOrder + j) * sz[0] +
                                         (ei[0] * eOrder + i)] < refine_th) {
                                if ((pNodes[ele].getLevel() +
                                     MAXDEAPTH_LEVEL_DIFF + 1) < m_uiMaxDepth)
                                    pNodes[ele].setFlag(
                                        ((OCT_SPLIT << NUM_LEVEL_BITS) |
                                         pNodes[ele].getLevel()));
                            }
                        }
            }
        }

        // coarsen test.
        for (unsigned int b = 0; b < blkList.size(); b++) {
            const ot::TreeNode blkNode = blkList[b].getBlockNode();

            sz[0]                      = blkList[b].getAllocationSzX();
            sz[1]                      = blkList[b].getAllocationSzY();
            sz[2]                      = blkList[b].getAllocationSzZ();

            const unsigned int bflag   = blkList[b].getBlkNodeFlag();
            const unsigned int offset  = blkList[b].getOffset();

            const unsigned int regLev  = blkList[b].getRegularGridLev();
            const unsigned int eleIndexMax =
                (1u << (regLev - blkNode.getLevel())) - 1;
            const unsigned int eleIndexMin = 0;

            if ((eleIndexMax == 0) || (bflag != 0))
                continue;  // this implies the blocks with only 1 child and
                           // boundary blocks.

            for (unsigned int ele = blkList[b].getLocalElementBegin();
                 ele < blkList[b].getLocalElementEnd(); ele++) {
                assert(pNodes[ele].getParent() ==
                       pNodes[ele + NUM_CHILDREN - 1].getParent());
                bool isCoarsen = true;

                for (unsigned int child = 0; child < NUM_CHILDREN; child++) {
                    if ((pNodes[ele + child].getFlag() >> NUM_LEVEL_BITS) ==
                        OCT_SPLIT) {
                        isCoarsen = false;
                        break;
                    }
                }

                if (isCoarsen && pNodes[ele].getLevel() > 1) {
                    bool coarse = true;
                    for (unsigned int child = 0; child < NUM_CHILDREN;
                         child++) {
                        ei[0] = (pNodes[ele + child].getX() - blkNode.getX()) >>
                                (m_uiMaxDepth - regLev);
                        ei[1] = (pNodes[ele + child].getY() - blkNode.getY()) >>
                                (m_uiMaxDepth - regLev);
                        ei[2] = (pNodes[ele + child].getZ() - blkNode.getZ()) >>
                                (m_uiMaxDepth - regLev);

                        for (unsigned int k = 3; k < eOrder + 1 + 3; k++)
                            for (unsigned int j = 3; j < eOrder + 1 + 3; j++)
                                for (unsigned int i = 3; i < eOrder + +3; i++) {
                                    if (!((refine_th <
                                           unzipVec[vIndex]
                                                   [offset +
                                                    (ei[2] * eOrder + k) *
                                                        sz[0] * sz[1] +
                                                    (ei[1] * eOrder + j) *
                                                        sz[0] +
                                                    (ei[0] * eOrder + i)]) &&
                                          (unzipVec[vIndex]
                                                   [offset +
                                                    (ei[2] * eOrder + k) *
                                                        sz[0] * sz[1] +
                                                    (ei[1] * eOrder + j) *
                                                        sz[0] +
                                                    (ei[0] * eOrder + i)] <=
                                           coarsen_th)))
                                        coarse = false;
                                }
                    }

                    if (coarse)
                        for (unsigned int child = 0; child < NUM_CHILDREN;
                             child++)
                            pNodes[ele + child].setFlag(
                                ((OCT_COARSE << NUM_LEVEL_BITS) |
                                 pNodes[ele].getLevel()));
                }

                ele = ele + NUM_CHILDREN - 1;
            }
        }

        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++)
            if ((pNodes[ele].getFlag() >> NUM_LEVEL_BITS) ==
                OCT_SPLIT)  // trigger remesh only when some refinement occurs
                            // (laid back remesh :)  )
            {
                isOctChange = true;
                break;
            }
    }

    bool isOctChanged_g;
    MPI_Allreduce(&isOctChange, &isOctChanged_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    // if(!m_uiGlobalRank) std::cout<<"is oct changed:
    // "<<isOctChanged_g<<std::endl;
    return isOctChanged_g;
}

bool isReMeshWAMR(
    ot::Mesh* pMesh, const double** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double, double*)> wavelet_tol,
    double amr_coarse_fac) {
    // if(!(pMesh->isReMeshUnzip((const double
    // **)unzippedVec,varIds,numVars,wavelet_tol,bssn::BSSN_DENDRO_AMR_FAC)))
    //     return false;

    // if(bssn::BSSN_CURRENT_RK_COORD_TIME > 0 &&
    // bssn::BSSN_CURRENT_RK_COORD_TIME < 80)
    //     return bssn::isReMeshBHRadial(pMesh);

    std::vector<unsigned int> refine_flags;
    const double r_near[2] = {bssn::BSSN_BH1_AMR_R, bssn::BSSN_BH2_AMR_R};

    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
    bool isOctChange                 = false;
    bool isOctChange_g               = false;
    Point d1, d2, temp;

    const unsigned int eOrder = pMesh->getElementOrder();
    const double dBH          = (BSSN_BH_LOC[0] - BSSN_BH_LOC[1]).abs();
    const unsigned int refLevMin =
        std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);

    // BH considered merged if the distance between punctures are less than the
    // specified value.
    const double BH_MERGED_SEP_TOL = 0.1;

    if (pMesh->isActive()) {
        if (!pMesh->getMPIRank()) printf("BH coord sep: %.8E \n", dBH);
        // std::cout<<"BH coord sep: "<<dBH<<std::endl;

        const RefElement* refEl    = pMesh->getReferenceElement();
        wavelet::WaveletEl* wrefEl = new wavelet::WaveletEl((RefElement*)refEl);

        refine_flags.resize(pMesh->getNumLocalMeshElements(), OCT_NO_CHANGE);
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();

        std::vector<double> wtol_vals;
        wtol_vals.resize(BSSN_NUM_VARS, 0);

        const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
        const unsigned int eOrder             = pMesh->getElementOrder();

        const unsigned int nx                 = (2 * eOrder + 1);
        const unsigned int ny                 = (2 * eOrder + 1);
        const unsigned int nz                 = (2 * eOrder + 1);

        const unsigned int sz_per_dof         = nx * ny * nz;
        const unsigned int isz[]              = {nx, ny, nz};
        std::vector<double> eVecTmp;
        eVecTmp.resize(sz_per_dof);

        std::vector<double> wCout;
        wCout.resize(sz_per_dof);

        for (unsigned int blk = 0; blk < blkList.size(); blk++) {
            const unsigned int pw    = blkList[blk].get1DPadWidth();
            const unsigned int bflag = blkList[blk].getBlkNodeFlag();
            assert(pw == (eOrder >> 1u));

            for (unsigned int ele = blkList[blk].getLocalElementBegin();
                 ele < blkList[blk].getLocalElementEnd(); ele++) {
                const bool isBdyOct = pMesh->isBoundaryOctant(ele);
                const double oct_dx =
                    (1u << (m_uiMaxDepth - pNodes[ele].getLevel())) /
                    (double(eOrder));

                Point oct_pt1 = Point(pNodes[ele].minX(), pNodes[ele].minY(),
                                      pNodes[ele].minZ());
                Point oct_pt2 = Point(pNodes[ele].minX() + oct_dx,
                                      pNodes[ele].minY() + oct_dx,
                                      pNodes[ele].minZ() + oct_dx);
                Point domain_pt1, domain_pt2, dx_domain;
                pMesh->octCoordToDomainCoord(oct_pt1, domain_pt1);
                pMesh->octCoordToDomainCoord(oct_pt2, domain_pt2);
                dx_domain    = domain_pt2 - domain_pt1;
                double hx[3] = {dx_domain.x(), dx_domain.y(), dx_domain.z()};
                const double tol_ele = wavelet_tol(
                    domain_pt1.x(), domain_pt1.y(), domain_pt1.z(), hx);

                // initialize all the wavelet errors to zero initially.
                for (unsigned int v = 0; v < BSSN_NUM_VARS; v++)
                    wtol_vals[v] = 0;

                // calculate L2 norm of wavelet coefficients
                // for each variable; if any exceed limit,
                // break out early
                for (unsigned int v = 0; v < numVars; v++) {
                    const unsigned int vid = varIds[v];
                    pMesh->getUnzipElementalNodalValues(
                        unzippedVec[vid], blk, ele, eVecTmp.data(), true);

                    // computes the wavelets.
                    wrefEl->compute_wavelets_3D((double*)(eVecTmp.data()), isz,
                                                wCout, isBdyOct);
                    // calculate the L-infinity norm of values on the grid
                    // double Linf = normLInfty(eVecTmp.data(), eVecTmp.size());
                    double Linf    = 1.0;  // remove this for relWAMR
                    // ensure not exploding with the relative values (div by 0)
                    Linf           = std::max(Linf, 1e-2);
                    // renormalize waveleth values as relative to L-infty norm
                    wtol_vals[vid] = (normL2(wCout.data(), wCout.size())) /
                                     sqrt(wCout.size()) / Linf;

                    // early bail if the computed tolerance value is large.
                    if (wtol_vals[vid] > tol_ele) break;
                }
                // compute maximum wavelet coefficient of the compiled set
                const double l_max = vecMax(wtol_vals.data(), wtol_vals.size());

                // use max wavelet coefficient to decide whether
                // to coaresen / refine / no change
                if (l_max > tol_ele) {
                    // if under-refined, then refine
                    refine_flags[(ele - eleLocalBegin)] = OCT_SPLIT;
                } else if (l_max < amr_coarse_fac * tol_ele) {
                    // if over-refined, then coarsen
                    refine_flags[(ele - eleLocalBegin)] = OCT_COARSE;
                } else {
                    // Goldilox zone - no changes needed
                    refine_flags[(ele - eleLocalBegin)] = OCT_NO_CHANGE;
                }
            }
        }

        delete wrefEl;
        // end of WAMR core calculation.

        ////////////////////////////////////////////////////////////////
        // Below code enforces a certain level of refinement at the BHs,
        // overriding what's currently set by the wavelets.
        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++) {
            // refine_flags[ele-eleLocalBegin] =
            // (pNodes[ele].getFlag()>>NUM_LEVEL_BITS); std::cout<<"ref flag:
            // "<<(pNodes[ele].getFlag()>>NUM_LEVEL_BITS)<<std::endl;
            // if(refine_flags[ele-eleLocalBegin]==OCT_SPLIT)
            pMesh->octCoordToDomainCoord(
                Point((double)pNodes[ele].minX(), (double)pNodes[ele].minY(),
                      (double)pNodes[ele].minZ()),
                temp);
            d1 = temp - BSSN_BH_LOC[0];
            d2 = temp - BSSN_BH_LOC[1];

            //@milinda: 11/21/2020 : Don't allow to violate the min depth
            if (pNodes[ele].getLevel() < bssn::BSSN_MINDEPTH) {
                refine_flags[ele - eleLocalBegin] = OCT_SPLIT;
            } else if (pNodes[ele].getLevel() == bssn::BSSN_MINDEPTH &&
                       refine_flags[ele - eleLocalBegin] == OCT_COARSE) {
                refine_flags[ele - eleLocalBegin] = OCT_NO_CHANGE;
            }

            // don't overide things away from puntures let wavelets handle that.
            if (d1.abs() > 10 && d2.abs() > 10)
                continue;
            else {
                const unsigned int ln =
                    1u << (m_uiMaxDepth - pNodes[ele].getLevel());
                const double hx = ln / (double)(eOrder);
                for (unsigned int k = 0; k < (eOrder + 1); k++)
                    for (unsigned int j = 0; j < (eOrder + 1); j++)
                        for (unsigned int i = 0; i < (eOrder + 1); i++) {
                            const double x      = pNodes[ele].minX() + k * hx;
                            const double y      = pNodes[ele].minY() + j * hx;
                            const double z      = pNodes[ele].minZ() + i * hx;
                            const Point oct_mid = Point(x, y, z);

                            pMesh->octCoordToDomainCoord(oct_mid, temp);

                            d1                     = temp - BSSN_BH_LOC[0];
                            d2                     = temp - BSSN_BH_LOC[1];

                            // std::cout<<"d1: "<<d1 <<
                            // "BHLOC_0:"<<BSSN_BH_LOC[0]<<std::endl;
                            // std::cout<<"d2: "<<d2<<std::endl;

                            const double rd1       = d1.abs();
                            const double rd2       = d2.abs();

                            const bool isNearTobh1 = (rd1 <= r_near[0]);
                            const bool isNearTobh2 = (rd2 <= r_near[1]);

                            const bool isMidNearTobh1 =
                                (rd1 > r_near[0] && rd1 <= 10.0 * r_near[0]);
                            const bool isMidNearTobh2 =
                                (rd2 > r_near[1] && rd1 <= 10.0 * r_near[1]);

                            const bool isFarTobh1 = (rd1 > 2.0 * r_near[0]);
                            const bool isFarTobh2 = (rd2 > 2.0 * r_near[1]);

                            if (dBH < BH_MERGED_SEP_TOL) {
                                if (isNearTobh1 || isNearTobh2) {
                                    // std::cout<<"d1:
                                    // "<<d1.abs()<<"BHLOC_0:"<<BSSN_BH_LOC[0]<<std::endl;
                                    // std::cout<<"d2:
                                    // "<<d2.abs()<<"BHLOC_1:"<<BSSN_BH_LOC[1]<<std::endl;

                                    if ((pNodes[ele].getLevel() +
                                         MAXDEAPTH_LEVEL_DIFF + 1) < refLevMin)
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_SPLIT;
                                    else if ((pNodes[ele].getLevel() +
                                              MAXDEAPTH_LEVEL_DIFF + 1) >
                                             refLevMin)
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                    else
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_NO_CHANGE;

                                    // if( ( pNodes[ele].getLevel() +
                                    // MAXDEAPTH_LEVEL_DIFF +1)== refLevMin )
                                    //     refine_flags[ele-eleLocalBegin] =
                                    //     OCT_NO_CHANGE;
                                    // else if(( pNodes[ele].getLevel() +
                                    // MAXDEAPTH_LEVEL_DIFF +1)> refLevMin)
                                    //     refine_flags[ele-eleLocalBegin] =
                                    //     OCT_COARSE;
                                }

                            } else {
                                if (bssn::BSSN_BH1_MAX_LEV == refLevMin) {
                                    if (isNearTobh1) {
                                        // std::cout<<"d1:
                                        // "<<d1.abs()<<"BHLOC_0:"<<BSSN_BH_LOC[0]<<"
                                        // rnear: "<<r_near[0]<<std::endl;
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;

                                    } else {
                                        if (refine_flags[ele - eleLocalBegin] ==
                                                OCT_SPLIT &&
                                            (pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) ==
                                                bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                    }

                                    // changes in bh 1 will get overidden by lev
                                    // 2
                                    if (isNearTobh2) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    }

                                } else {
                                    assert(bssn::BSSN_BH2_MAX_LEV == refLevMin);
                                    if (isNearTobh2) {
                                        // std::cout<<"d1:
                                        // "<<d1.abs()<<"BHLOC_0:"<<BSSN_BH_LOC[0]<<"
                                        // rnear: "<<r_near[0]<<std::endl;
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;

                                    } else {
                                        if (refine_flags[ele - eleLocalBegin] ==
                                                OCT_SPLIT &&
                                            (pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) ==
                                                bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                    }

                                    // changes in bh 2 will get overidden by lev
                                    // 1 which is the higher level than bh2.
                                    if (isNearTobh1) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    }
                                }
                            }
                        }
            }
        }

        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    // communicate refinement between cores.
    MPI_Allreduce(&isOctChange, &isOctChange_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChange_g;
}

// Don't use this : 1) This is expensive
// 2). It has a coarsening bug (doesn't coarsen properly)
bool isReMeshBHRadial(ot::Mesh* pMesh) {
    std::vector<unsigned int> refine_flags;
    const double r_near[2] = {bssn::BSSN_BH1_AMR_R, bssn::BSSN_BH2_AMR_R};

    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
    bool isOctChange                 = false;
    bool isOctChange_g               = false;
    Point d1, d2, temp;

    const unsigned int eOrder = pMesh->getElementOrder();
    const double dBH          = (BSSN_BH_LOC[0] - BSSN_BH_LOC[1]).abs();
    const unsigned int refLevMin =
        std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);

    // BH considered merged if the distance between punctures are less than the
    // specified value.
    const double BH_MERGED_SEP_TOL        = 0.1;
    const unsigned int NUM_REFINE_SPHERES = 5;

    std::vector<double> bh1_amr_r;
    std::vector<double> bh2_amr_r;

    bh1_amr_r.push_back(0.0);
    bh2_amr_r.push_back(0.0);

    for (unsigned int i = 0; i < NUM_REFINE_SPHERES; i++) {
        const unsigned int rfac = (1u << i);
        bh1_amr_r.push_back(rfac * bssn::BSSN_BH1_AMR_R);
        bh2_amr_r.push_back(rfac * bssn::BSSN_BH2_AMR_R);
    }

    if (pMesh->isActive()) {
        refine_flags.resize(pMesh->getNumLocalMeshElements(), OCT_NO_CHANGE);
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();

        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++) {
            // refine_flags[ele-eleLocalBegin] =
            // (pNodes[ele].getFlag()>>NUM_LEVEL_BITS); std::cout<<"ref flag:
            // "<<(pNodes[ele].getFlag()>>NUM_LEVEL_BITS)<<std::endl;
            // if(refine_flags[ele-eleLocalBegin]==OCT_SPLIT)
            pMesh->octCoordToDomainCoord(
                Point((double)pNodes[ele].minX(), (double)pNodes[ele].minY(),
                      (double)pNodes[ele].minZ()),
                temp);
            d1 = temp - BSSN_BH_LOC[0];
            d2 = temp - BSSN_BH_LOC[1];

            //@milinda: 11/21/2020 : Don't allow to violate the min depth
            if (pNodes[ele].getLevel() < bssn::BSSN_MINDEPTH) {
                refine_flags[ele - eleLocalBegin] = OCT_SPLIT;
            } else if (pNodes[ele].getLevel() == bssn::BSSN_MINDEPTH &&
                       refine_flags[ele - eleLocalBegin] == OCT_COARSE) {
                refine_flags[ele - eleLocalBegin] = OCT_NO_CHANGE;
            }

            const unsigned int ln = 1u
                                    << (m_uiMaxDepth - pNodes[ele].getLevel());
            const double hx = ln / (double)(eOrder);

            for (unsigned int k = 0; k < (eOrder + 1); k++)
                for (unsigned int j = 0; j < (eOrder + 1); j++)
                    for (unsigned int i = 0; i < (eOrder + 1); i++) {
                        const double x      = pNodes[ele].minX() + i * hx;
                        const double y      = pNodes[ele].minY() + j * hx;
                        const double z      = pNodes[ele].minZ() + k * hx;
                        const Point oct_mid = Point(x, y, z);

                        pMesh->octCoordToDomainCoord(oct_mid, temp);

                        d1               = temp - BSSN_BH_LOC[0];
                        d2               = temp - BSSN_BH_LOC[1];

                        // std::cout<<"d1: "<<d1 <<
                        // "BHLOC_0:"<<BSSN_BH_LOC[0]<<std::endl;
                        // std::cout<<"d2: "<<d2<<std::endl;

                        const double rd1 = d1.abs();
                        const double rd2 = d2.abs();

                        if (dBH < BH_MERGED_SEP_TOL) {
                            for (unsigned int rs = 1;
                                 rs < NUM_REFINE_SPHERES + 1; rs++) {
                                if ((rd1 > bh1_amr_r[rs - 1]) &&
                                    (rd1 <= bh1_amr_r[rs])) {
                                    if ((pNodes[ele].getLevel() +
                                         MAXDEAPTH_LEVEL_DIFF + 1) <
                                        std::max(refLevMin - (rs - 1),
                                                 bssn::BSSN_MINDEPTH))
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_SPLIT;
                                    else if ((pNodes[ele].getLevel() +
                                              MAXDEAPTH_LEVEL_DIFF + 1) >
                                             std::max(refLevMin - (rs - 1),
                                                      bssn::BSSN_MINDEPTH))
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                    else
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_NO_CHANGE;
                                }
                                if (rd1 >
                                    bh1_amr_r.back())  // note : It is important
                                                       // to keep this inside
                                                       // the for loop to ensure
                                                       // proper refinement.
                                    refine_flags[ele - eleLocalBegin] =
                                        OCT_COARSE;
                            }

                        } else {
                            if (bssn::BSSN_BH1_MAX_LEV == refLevMin) {
                                for (unsigned int rs = 1;
                                     rs < NUM_REFINE_SPHERES + 1; rs++) {
                                    if ((rd1 > bh1_amr_r[rs - 1]) &&
                                        (rd1 <= bh1_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;

                                    } else if (rd1 > bh1_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;

                                    // overide by smaller bh - BH2 is smaller
                                    // from the depth parameter.
                                    if ((rd2 > bh2_amr_r[rs - 1]) &&
                                        (rd2 <= bh2_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    } else if (rd2 > bh2_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                }

                            } else {
                                assert(bssn::BSSN_BH2_MAX_LEV == refLevMin);
                                for (unsigned int rs = 1;
                                     rs < NUM_REFINE_SPHERES + 1; rs++) {
                                    if ((rd2 > bh2_amr_r[rs - 1]) &&
                                        (rd2 <= bh2_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    } else if (rd2 > bh2_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;

                                    // overide by smaller bh - BH1 is smaller
                                    // from the depth parameter.
                                    if ((rd1 > bh1_amr_r[rs - 1]) &&
                                        (rd1 <= bh1_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    } else if (rd1 > bh1_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                }
                            }
                        }
                    }
        }
        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    MPI_Allreduce(&isOctChange, &isOctChange_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChange_g;
}

bool addRemeshWAMR(
    ot::Mesh* pMesh, const double** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double, double*)> wavelet_tol,
    double amr_coarse_fac, bool relative_WAMR) {
    // WKB Aug 2024: add new refinement helper to add on refinement
    // in regions where wavelets say so, but don't use wavelets
    // directly for de-refinement. This adds up to this equalling
    // isReMeshWAMR, but without calling 'COARSEN' (nor trying to do
    // its own BHLB method

    // pull in current refinement flags
    std::vector<unsigned int> refine_flags = pMesh->getAllRefinementFlags();

    const double r_near[2] = {bssn::BSSN_BH1_AMR_R, bssn::BSSN_BH2_AMR_R};

    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
    bool isOctChange                 = false;
    bool isOctChange_g               = false;
    Point d1, d2, temp;

    const unsigned int eOrder = pMesh->getElementOrder();
    const double dBH          = (BSSN_BH_LOC[0] - BSSN_BH_LOC[1]).abs();
    const unsigned int refLevMin =
        std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);

    // BH considered merged if the distance between punctures are less than the
    // specified value.
    const double BH_MERGED_SEP_TOL = 0.1;

    if (pMesh->isActive()) {
        if (!pMesh->getMPIRank()) printf("BH coord sep: %.8E \n", dBH);
        // std::cout<<"BH coord sep: "<<dBH<<std::endl;

        const RefElement* refEl    = pMesh->getReferenceElement();
        wavelet::WaveletEl* wrefEl = new wavelet::WaveletEl((RefElement*)refEl);

        // pull current refinement flags
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();

        // set up vector to collect wavelet tolerance values
        std::vector<double> wtol_vals(BSSN_NUM_VARS, 0.0);

        const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
        const unsigned int eOrder             = pMesh->getElementOrder();

        const unsigned int nx                 = (2 * eOrder + 1);
        const unsigned int ny                 = (2 * eOrder + 1);
        const unsigned int nz                 = (2 * eOrder + 1);

        const unsigned int sz_per_dof         = nx * ny * nz;
        const unsigned int isz[]              = {nx, ny, nz};
        std::vector<double> eVecTmp;
        eVecTmp.resize(sz_per_dof);

        std::vector<double> wCout;
        wCout.resize(sz_per_dof);

        for (unsigned int blk = 0; blk < blkList.size(); blk++) {
            const unsigned int pw    = blkList[blk].get1DPadWidth();
            const unsigned int bflag = blkList[blk].getBlkNodeFlag();
            assert(pw == (eOrder >> 1u));

            for (unsigned int ele = blkList[blk].getLocalElementBegin();
                 ele < blkList[blk].getLocalElementEnd(); ele++) {
                const bool isBdyOct = pMesh->isBoundaryOctant(ele);
                const double oct_dx =
                    (1u << (m_uiMaxDepth - pNodes[ele].getLevel())) /
                    (double(eOrder));

                Point oct_pt1 = Point(pNodes[ele].minX(), pNodes[ele].minY(),
                                      pNodes[ele].minZ());
                Point oct_pt2 = Point(pNodes[ele].minX() + oct_dx,
                                      pNodes[ele].minY() + oct_dx,
                                      pNodes[ele].minZ() + oct_dx);
                Point domain_pt1, domain_pt2, dx_domain;
                pMesh->octCoordToDomainCoord(oct_pt1, domain_pt1);
                pMesh->octCoordToDomainCoord(oct_pt2, domain_pt2);
                dx_domain    = domain_pt2 - domain_pt1;
                double hx[3] = {dx_domain.x(), dx_domain.y(), dx_domain.z()};
                double tol_ele = wavelet_tol(
                    domain_pt1.x(), domain_pt1.y(), domain_pt1.z(), hx);
                
                ////////////////////////////////////////////////////////
                // wkb 17 Oct: amplify wavelet tolerance if w/i BH radii
                // calculate distances to each black hole
                const double r_BH1 = (domain_pt1 - BSSN_BH_LOC[0]).abs() / bssn::BSSN_BH1_AMR_R;
                const double r_BH2 = (domain_pt1 - BSSN_BH_LOC[1]).abs() / bssn::BSSN_BH2_AMR_R;
                // overwrite if we're w/i the AMR radius of either BH
                if (std::min(r_BH1,r_BH2) <= 1) {
                    // large amplification of tolerance should 
                    // effectively disable WAMR w/i the BHs.
                    tol_ele *= 1e12; 
                }
                
                // initialize all the wavelet errors to zero initially.
                for (unsigned int v = 0; v < BSSN_NUM_VARS; v++)
                    wtol_vals[v] = 0;

                // calculate L2 norm of wavelet tolerances for each
                // variable; if any exceed limit, break out early
                for (unsigned int v = 0; v < numVars; v++) {
                    const unsigned int vid = varIds[v];
                    pMesh->getUnzipElementalNodalValues(
                        unzippedVec[vid], blk, ele, eVecTmp.data(), true);

                    // computes the wavelets.
                    wrefEl->compute_wavelets_3D((double*)(eVecTmp.data()), isz,
                                                wCout, isBdyOct);

                    // toggle using relative wavelets or not
                    double Linf;  // initialize normalization
                    if (relative_WAMR) {
                        // calculate the L-infinity norm of values on the grid
                        Linf = normLInfty(eVecTmp.data(), eVecTmp.size());
                        // ensure not exploding with the relative values (div by
                        // 0)
                        Linf = std::max(Linf, 1e-2);
                    } else {
                        Linf = 1.0;
                    }

                    // renormalize wavelet tolerance values
                    wtol_vals[vid] = (normL2(wCout.data(), wCout.size())) /
                                     sqrt(wCout.size()) / Linf;

                    // early bail if the computed tolerance value is large.
                    if (wtol_vals[vid] > tol_ele) break;
                }
                const double l_max = vecMax(wtol_vals.data(), wtol_vals.size());

                // if we exceed the tolerance, then add refinement,
                // overwriting whatever the current flag is
                if (l_max > tol_ele) {
                    refine_flags[(ele - eleLocalBegin)] = OCT_SPLIT;
                } else if (refine_flags[(ele - eleLocalBegin)] == OCT_COARSE &&
                           l_max >= amr_coarse_fac * tol_ele) {
                    // if it wanted to coarsen, but the current level of
                    // refinement is actually necessary, then we're in
                    // the Goldilocks zone: Don't change!
                    refine_flags[(ele - eleLocalBegin)] = OCT_NO_CHANGE;
                }
            }
        }

        delete wrefEl;

        // check whether any changes have been made to the grid here
        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    // check whether changes have been made *anywhere* on the grid
    MPI_Allreduce(&isOctChange, &isOctChange_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChange_g;
}

}  // end of namespace bssn
