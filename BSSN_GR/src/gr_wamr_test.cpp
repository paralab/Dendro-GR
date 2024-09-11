/**
 * @file gr_wamr_test.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Explore the Wavelet AMR errors and quantifications on that for BSSN
 * system.
 * @version 0.1
 * @date 2021-03-22
 * @copyright Copyright (c) 2021
 *
 */

#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "gr.h"
#include "grUtils.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"
#include "rkBSSN.h"

int main(int argc, char** argv) {
    if (argc < 2)
        std::cout << "Usage: " << argv[0] << " paramFile" << std::endl;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    // Print out CMAKE options
    if (!rank) {
#ifdef BSSN_COMPUTE_CONSTRAINTS
        std::cout << GRN << "  Compiled with BSSN_COMPUTE_CONSTRAINTS" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without BSSN_COMPUTE_CONSTRAINTS" << NRM
                  << std::endl;
#endif
#ifdef BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT
        std::cout << GRN << "  Compiled with BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"
                  << NRM << std::endl;
#endif
#ifdef BSSN_ENABLE_VTU_OUTPUT
        std::cout << GRN << "  Compiled with BSSN_ENABLE_VTU_OUTPUT" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without BSSN_ENABLE_VTU_OUTPUT" << NRM
                  << std::endl;
#endif
#ifdef BSSN_ETA_FUNCTION
        std::cout << GRN << "  Compiled with  BSSN_ETA_FUNCTION" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_ETA_FUNCTION" << NRM
                  << std::endl;
#endif
#ifdef BSSN_EXTRACT_BH_LOCATIONS
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_BH_LOCATIONS" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_EXTRACT_BH_LOCATIONS"
                  << NRM << std::endl;
#endif
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#endif
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        std::cout << GRN << "  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#else
        std::cout << RED
                  << "  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"
                  << NRM << std::endl;
#endif
#ifdef BSSN_GAUGE_ROCHESTER
        std::cout << GRN << "  Compiled with  BSSN_GAUGE_ROCHESTER" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_GAUGE_ROCHESTER" << NRM
                  << std::endl;
#endif
#ifdef BSSN_KERR_SCHILD_TEST
        std::cout << GRN << "  Compiled with  BSSN_KERR_SCHILD_TEST" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_KERR_SCHILD_TEST" << NRM
                  << std::endl;
#endif

#ifdef BSSN_REFINE_BASE_EH
        std::cout << GRN << "  Compiled with  BSSN_REFINE_BASE_EH" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  BSSN_REFINE_BASE_EH" << NRM
                  << std::endl;
#endif

#ifdef USE_FD_INTERP_FOR_UNZIP
        std::cout << GRN << "  Compiled with  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#else
        std::cout << RED << "  Compiled without  USE_FD_INTERP_FOR_UNZIP" << NRM
                  << std::endl;
#endif

#ifdef BSSN_USE_4TH_ORDER_DERIVS
        std::cout
            << GRN
            << "  Using 4th order FD interpolations (2nd order on boundary)"
            << NRM << std::endl;
#endif

#ifdef BSSN_USE_6TH_ORDER_DERIVS
        std::cout
            << GRN
            << "  Using 6th order FD interpolations (2nd order on boundary)"
            << NRM << std::endl;
#endif

#ifdef BSSN_USE_8TH_ORDER_DERIVS
        std::cout
            << GRN
            << "  Using 8th order FD interpolations (2nd order on boundary)"
            << NRM << std::endl;
#endif
    }

    // 1 . read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    bssn::readParamFile(argv[1], comm);

    int root = std::min(1, npes - 1);
    bssn::dumpParamFile(std::cout, root, comm);

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth = bssn::BSSN_MAXDEPTH;

    if (bssn::BSSN_NUM_VARS % bssn::BSSN_ASYNC_COMM_K != 0) {
        if (!rank)
            std::cout << "[overlap communication error]: total BSSN_NUM_VARS: "
                      << bssn::BSSN_NUM_VARS
                      << " is not divisable by BSSN_ASYNC_COMM_K: "
                      << bssn::BSSN_ASYNC_COMM_K << std::endl;
        MPI_Abort(comm, 0);
    }

    const unsigned int interpVars = bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) varIndex[i] = i;

    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::punctureData(x, y, z, var);
        };
    std::function<void(double, double, double, double*)> f_init_phy =
        [](double x, double y, double z, double* var) {
            bssn::punctureDataPhysicalCoord(x, y, z, var);
        };
    std::function<double(double, double, double)> f_init_alpha =
        [](double x, double y, double z) {
            double var[24];
            bssn::punctureData(x, y, z, var);
            return var[0];
        };

    if (!rank)
        std::cout << YLW << " Generating block AMR mesh " << NRM << std::endl;
    const unsigned int f2olmin =
        std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);
    unsigned int lmin, lmax;

    // const Point
    // pt_min(bssn::BSSN_BLK_MIN_X,bssn::BSSN_BLK_MIN_Y,bssn::BSSN_BLK_MIN_Z);
    // const Point
    // pt_max(bssn::BSSN_BLK_MAX_X,bssn::BSSN_BLK_MAX_Y,bssn::BSSN_BLK_MAX_Z);

    // bssn::blockAdaptiveOctree(tmpNodes,pt_min,pt_max, f2olmin-2,
    // m_uiMaxDepth,comm); ot::Mesh * mesh_blk =
    // ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
    // mesh_blk->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z),
    // Point(bssn::BSSN_GRID_MAX_X,
    // bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
    // mesh_blk->computeMinMaxLevel(lmin,lmax);

    // if(!rank)
    // {
    //     std::cout<<"================= Grid Info (Before init grid
    //     converge):======================================================="<<std::endl;
    //     std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
    //     std::cout<<"dx:
    //     "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double)
    //     bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
    //     std::cout<<"dt:
    //     "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double)
    //     bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
    //     std::cout<<"==============================================================================================================="<<std::endl;
    // }

    // DendroScalar * blk_v1  = mesh_blk ->
    // createCGVector<DendroScalar>(f_init,bssn::BSSN_NUM_VARS);
    // mesh_blk->readFromGhostBegin(blk_v1,bssn::BSSN_NUM_VARS);
    // mesh_blk->readFromGhostEnd(blk_v1,bssn::BSSN_NUM_VARS);

    // {

    //     char fPrefix[256];
    //     sprintf(fPrefix,"%s_blk",bssn::BSSN_VTU_FILE_PREFIX.c_str());

    //     std::vector<std::string> pDataNames;
    //     double *pData[bssn::BSSN_NUM_VARS];

    //     for(unsigned int i=0; i < bssn::BSSN_NUM_VARS; i++)
    //     {
    //         pDataNames.push_back(std::string(bssn::BSSN_VAR_NAMES[i]));
    //         pData[i] = bssn_blk  +  i * mesh_blk->getDegOfFreedom();
    //     }

    //     std::vector<char*> pDataNames_char;
    //     pDataNames_char.reserve(pDataNames.size());

    //     for(unsigned int  i = 0; i < pDataNames.size(); i++)
    //         pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

    //     io::vtk::mesh2vtuFine(mesh_blk,fPrefix,0,NULL,NULL,bssn::BSSN_NUM_VARS,(const
    //     char **)&pDataNames_char[0],(const double **)pData);

    // }

    // mesh_blk->destroyVector(bssn_blk);
    // delete mesh_blk;

    if (!rank)
        std::cout << YLW << "Using function2Octree. Wavelet AMR enabled " << NRM
                  << std::endl;
    tmpNodes.clear();
    function2Octree(f_init, bssn::BSSN_NUM_VARS, varIndex, interpVars, tmpNodes,
                    (f2olmin - 2), bssn::BSSN_WAVELET_TOL, bssn::BSSN_ELE_ORDER,
                    comm);

    ot::Mesh* mesh_wamr =
        ot::createMesh(tmpNodes.data(), tmpNodes.size(), bssn::BSSN_ELE_ORDER,
                       comm, 1, ot::SM_TYPE::FDM, bssn::BSSN_DENDRO_GRAIN_SZ,
                       bssn::BSSN_LOAD_IMB_TOL, bssn::BSSN_SPLIT_FIX);
    mesh_wamr->setDomainBounds(
        Point(bssn::BSSN_GRID_MIN_X, bssn::BSSN_GRID_MIN_Y,
              bssn::BSSN_GRID_MIN_Z),
        Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,
              bssn::BSSN_GRID_MAX_Z));
    mesh_wamr->computeMinMaxLevel(lmin, lmax);
    if (!rank) {
        std::cout << "================= Grid Info (Before init grid "
                     "converge):==============================================="
                     "========"
                  << std::endl;
        std::cout << "lmin: " << lmin << " lmax:" << lmax << std::endl;
        std::cout << "dx: "
                  << ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                      ((1u << (m_uiMaxDepth - lmax)) /
                       ((double)bssn::BSSN_ELE_ORDER)) /
                      ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "dt: "
                  << bssn::BSSN_CFL_FACTOR *
                         ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                          ((1u << (m_uiMaxDepth - lmax)) /
                           ((double)bssn::BSSN_ELE_ORDER)) /
                          ((double)(1u << (m_uiMaxDepth))))
                  << std::endl;
        std::cout << "========================================================="
                     "======================================================"
                  << std::endl;
    }

    DendroScalar* wamr_v1 = mesh_wamr->createCGVector<DendroScalar>(
        f_init_phy, bssn::BSSN_NUM_VARS);

    mesh_wamr->readFromGhostBegin(wamr_v1, bssn::BSSN_NUM_VARS);
    mesh_wamr->readFromGhostEnd(wamr_v1, bssn::BSSN_NUM_VARS);

    ot::Mesh* mesh_wamr_refined = NULL;
    // refine mesh_wamr lmax-1 octant to lmax
    {
        bool isRefine_g, isRefine = false;
        ot::Mesh* mesh = mesh_wamr;
        std::vector<unsigned int> refine_flag;
        refine_flag.reserve(mesh->getNumLocalMeshElements());
        const ot::TreeNode* pNodes = mesh->getAllElements().data();
        for (unsigned int ele = mesh->getElementLocalBegin();
             ele < mesh->getElementLocalEnd(); ele++) {
            if (pNodes[ele].getLevel() == (lmax - 1)) {
                refine_flag.push_back(OCT_SPLIT);
                isRefine = true;
            } else
                refine_flag.push_back(OCT_NO_CHANGE);
        }

        MPI_Allreduce(&isRefine, &isRefine_g, 1, MPI_C_BOOL, MPI_LOR, comm);

        if (isRefine_g) {
            mesh->setMeshRefinementFlags(refine_flag);
            mesh_wamr_refined  = mesh->ReMesh();

            DendroIntL localSz = mesh->getNumLocalMeshElements();
            DendroIntL gSz_new, gSz_old;

            par::Mpi_Reduce(&localSz, &gSz_old, 1, MPI_SUM, 0, comm);
            localSz = mesh_wamr_refined->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz, &gSz_new, 1, MPI_SUM, 0, comm);

            if (!rank)
                std::cout << "old mesh size: " << gSz_old
                          << " new mesh size: " << gSz_new << std::endl;
        }
    }

    if (mesh_wamr_refined != NULL) {
        DendroScalar* wamr_v2 = mesh_wamr_refined->createCGVector<DendroScalar>(
            f_init_phy, bssn::BSSN_NUM_VARS);
        DendroScalar* wamr_v3 = mesh_wamr_refined->createCGVector<DendroScalar>(
            0, bssn::BSSN_NUM_VARS);
        mesh_wamr->interGridTransfer(wamr_v1, wamr_v3, mesh_wamr_refined,
                                     ot::INJECTION, bssn::BSSN_NUM_VARS);

        mesh_wamr_refined->readFromGhostBegin(wamr_v2, bssn::BSSN_NUM_VARS);
        mesh_wamr_refined->readFromGhostEnd(wamr_v2, bssn::BSSN_NUM_VARS);

        mesh_wamr_refined->readFromGhostBegin(wamr_v3, bssn::BSSN_NUM_VARS);
        mesh_wamr_refined->readFromGhostEnd(wamr_v3, bssn::BSSN_NUM_VARS);

        {
            char fPrefix[256];
            sprintf(fPrefix, "%s_wamr_refined",
                    bssn::BSSN_VTU_FILE_PREFIX.c_str());

            std::vector<std::string> pDataNames;
            double* pData[2 * bssn::BSSN_NUM_VARS];

            for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) {
                pDataNames.push_back(std::string(bssn::BSSN_VAR_NAMES[i]));
                pData[i] = wamr_v2 + i * mesh_wamr_refined->getDegOfFreedom();
            }

            for (unsigned int i = 0; i < bssn::BSSN_NUM_VARS; i++) {
                pDataNames.push_back("IGT_" +
                                     std::string(bssn::BSSN_VAR_NAMES[i]));
                pData[bssn::BSSN_NUM_VARS + i] =
                    wamr_v3 + i * mesh_wamr_refined->getDegOfFreedom();
            }

            std::vector<char*> pDataNames_char;
            pDataNames_char.reserve(pDataNames.size());

            for (unsigned int i = 0; i < pDataNames.size(); i++)
                pDataNames_char.push_back(
                    const_cast<char*>(pDataNames[i].c_str()));

            io::vtk::mesh2vtuFine(mesh_wamr_refined, fPrefix, 0, NULL, NULL,
                                  2 * bssn::BSSN_NUM_VARS,
                                  (const char**)&pDataNames_char[0],
                                  (const double**)pData);
        }

        mesh_wamr_refined->destroyVector(wamr_v2);
        mesh_wamr_refined->destroyVector(wamr_v3);
    }

    // {

    //     char fPrefix[256];
    //     sprintf(fPrefix,"%s_wamr",bssn::BSSN_VTU_FILE_PREFIX.c_str());

    //     std::vector<std::string> pDataNames;
    //     double *pData[bssn::BSSN_NUM_VARS];

    //     for(unsigned int i=0; i < bssn::BSSN_NUM_VARS; i++)
    //     {
    //         pDataNames.push_back(std::string(bssn::BSSN_VAR_NAMES[i]));
    //         pData[i] = bssn_wamr  +  i * mesh_wamr->getDegOfFreedom();
    //     }

    //     std::vector<char*> pDataNames_char;
    //     pDataNames_char.reserve(pDataNames.size());

    //     for(unsigned int  i = 0; i < pDataNames.size(); i++)
    //         pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

    //     io::vtk::mesh2vtuFine(mesh_wamr,fPrefix,0,NULL,NULL,bssn::BSSN_NUM_VARS,(const
    //     char **)&pDataNames_char[0],(const double **)pData);

    // }

    mesh_wamr->destroyVector(wamr_v1);

    delete mesh_wamr;
    delete mesh_wamr_refined;

    MPI_Finalize();
    return 0;
}
