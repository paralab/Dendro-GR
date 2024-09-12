//
// Created by milinda on 7/25/17.
/**
 *@author Jacob Fields
 *Department of Physics and Astronomy, Brigham Young University
 *@brief Header file for the relativistic fluid simulation.
 */
//

#include "fluid.h"

#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "fluidUtils.h"
#include "mesh.h"
#include "mpi.h"
#include "octUtils.h"
#include "rkfluid.h"

int main(int argc, char** argv) {
    if (argc < 2)
        std::cout << "Usage: " << argv[0] << " paramFile" << std::endl;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    fluid::timer::initFlops();

    fluid::timer::total_runtime.start();

    // 1 . read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    fluid::readParamFile(argv[1], comm);

    if (rank == 1 || npes == 1) {
        std::cout << "parameters read: " << std::endl;
        std::cout << YLW << "\tfluid::FLUID_IO_OUTPUT_FREQ: "
                  << fluid::FLUID_IO_OUTPUT_FREQ << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_REMESH_TEST_FREQ: "
                  << fluid::FLUID_REMESH_TEST_FREQ << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_CHECKPT_FREQ: "
                  << fluid::FLUID_CHECKPT_FREQ << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_IO_OUTPUT_GAP: "
                  << fluid::FLUID_IO_OUTPUT_GAP << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_VTU_FILE_PREFIX: "
                  << fluid::FLUID_VTU_FILE_PREFIX << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_CHKPT_FILE_PREFIX: "
                  << fluid::FLUID_CHKPT_FILE_PREFIX << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_PROFILE_FILE_PREFIX: "
                  << fluid::FLUID_PROFILE_FILE_PREFIX << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_RESTORE_SOLVER: "
                  << fluid::FLUID_RESTORE_SOLVER << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_ENABLE_BLOCK_ADAPTIVITY: "
                  << fluid::FLUID_ENABLE_BLOCK_ADAPTIVITY << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_BLK_MIN_X: " << fluid::FLUID_BLK_MIN_X
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_BLK_MIN_Y: " << fluid::FLUID_BLK_MIN_Y
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_BLK_MIN_Z: " << fluid::FLUID_BLK_MIN_Z
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_BLK_MAX_X: " << fluid::FLUID_BLK_MAX_X
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_BLK_MAX_Y: " << fluid::FLUID_BLK_MAX_Y
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_BLK_MAX_Z: " << fluid::FLUID_BLK_MAX_Z
                  << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_DENDRO_GRAIN_SZ: "
                  << fluid::FLUID_DENDRO_GRAIN_SZ << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_ASYNC_COMM_K: "
                  << fluid::FLUID_ASYNC_COMM_K << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_DENDRO_AMR_FAC: "
                  << fluid::FLUID_DENDRO_AMR_FAC << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_LOAD_IMB_TOL: "
                  << fluid::FLUID_LOAD_IMB_TOL << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_RK_TIME_BEGIN: "
                  << fluid::FLUID_RK_TIME_BEGIN << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_RK_TIME_END: " << fluid::FLUID_RK_TIME_END
                  << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_RK45_TIME_STEP_SIZE: "
                  << fluid::FLUID_RK45_TIME_STEP_SIZE << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_RK45_DESIRED_TOL: "
                  << fluid::FLUID_RK45_DESIRED_TOL << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_DIM: " << fluid::FLUID_DIM << NRM
                  << std::endl;
        std::cout << YLW << "\tfluid::FLUID_MAXDEPTH: " << fluid::FLUID_MAXDEPTH
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_GRID_MIN_X: " << fluid::FLUID_GRID_MIN_X
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_GRID_MAX_X: " << fluid::FLUID_GRID_MAX_X
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_GRID_MIN_Y: " << fluid::FLUID_GRID_MIN_Y
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_GRID_MAX_Y: " << fluid::FLUID_GRID_MAX_Y
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_GRID_MIN_Z: " << fluid::FLUID_GRID_MIN_Z
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_GRID_MAX_Z: " << fluid::FLUID_GRID_MAX_Z
                  << NRM << std::endl;
        std::cout << YLW << "\tfluid::KO_DISS_SIGMA: " << fluid::KO_DISS_SIGMA
                  << NRM << std::endl;
        std::cout << YLW << "\tfluid::FLUID_ID_TYPE: " << fluid::FLUID_ID_TYPE
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tfluid::FLUID_WAVELET_TOL: " << fluid::FLUID_WAVELET_TOL
                  << NRM << std::endl;
    }

    _InitializeHcurve(fluid::FLUID_DIM);
    m_uiMaxDepth = fluid::FLUID_MAXDEPTH;

    /*if(fluid::FLUID_NUM_PRIM_VARS%fluid::FLUID_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total
    FLUID_NUM_PRIM_VARS: "<<fluid::FLUID_NUM_PRIM_VARS<<" is not divisable by
    FLUID_ASYNC_COMM_K: "<<fluid::FLUID_ASYNC_COMM_K<<std::endl; exit(0);
    }*/

    fluid::FLUID_RK45_TIME_STEP_SIZE =
        fluid::FLUID_CFL_FACTOR *
        (fluid::FLUID_COMPD_MAX[0] - fluid::FLUID_COMPD_MIN[0]) *
        (1.0 / (double)(1u << fluid::FLUID_MAXDEPTH));
    if (rank == 0) {
        std::cout << "CFL factor: " << fluid::FLUID_CFL_FACTOR << "\n";
    }

    // 2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            fluid::initData(x, y, z, var);
        };
    // std::function<void(double,double,double,double*)> f_init=[](double
    // x,double y,double z,double*var){fluid::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars = fluid::FLUID_NUM_PRIM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < fluid::FLUID_NUM_PRIM_VARS; i++)
        varIndex[i] = i;

    /*varIndex[0]=fluid::VAR::U_ALPHA;
    varIndex[1]=fluid::VAR::U_CHI;*/
    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    fluid::timer::t_f2o.start();

    if (fluid::FLUID_ENABLE_BLOCK_ADAPTIVITY) {
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled "
                      << NRM << std::endl;
        const Point pt_min(fluid::FLUID_BLK_MIN_X, fluid::FLUID_BLK_MIN_Y,
                           fluid::FLUID_BLK_MIN_Z);
        const Point pt_max(fluid::FLUID_BLK_MAX_X, fluid::FLUID_BLK_MAX_Y,
                           fluid::FLUID_BLK_MAX_Z);

        fluid::blockAdaptiveOctree(tmpNodes, pt_min, pt_max, m_uiMaxDepth - 2,
                                   m_uiMaxDepth, comm);
    } else {
        if (!rank)
            std::cout << YLW << "Using function2Octree. AMR enabled " << NRM
                      << std::endl;
        function2Octree(f_init, fluid::FLUID_NUM_PRIM_VARS, varIndex,
                        interpVars, tmpNodes, m_uiMaxDepth,
                        fluid::FLUID_WAVELET_TOL, fluid::FLUID_ELE_ORDER, comm);
    }

    fluid::timer::t_f2o.stop();

    t_stat = fluid::timer::t_f2o.seconds;
    par::Mpi_Reduce(&t_stat, t_stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(&t_stat, t_stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(&t_stat, t_stat_g + 2, 1, MPI_MAX, 0, comm);
    t_stat_g[1] = t_stat_g[1] / (double)npes;

    localSz     = tmpNodes.size();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);

    if (!rank)
        std::cout << GRN << " function to octree max (s): " << t_stat_g[2]
                  << NRM << std::endl;
    if (!rank)
        std::cout << GRN << " function to octree # octants : " << globalSz
                  << NRM << std::endl;

    par::Mpi_Bcast(&globalSz, 1, 0, comm);
    const unsigned int grainSz =
        fluid::FLUID_DENDRO_GRAIN_SZ;  // DENDRO_DEFAULT_GRAIN_SZ;

    bool isActive;
    MPI_Comm commActive;
    const int p_npes_prev =
        binOp::getPrevHighestPowerOfTwo((globalSz / grainSz));
    const int p_npes_next =
        binOp::getNextHighestPowerOfTwo((globalSz / grainSz));

    int p_npes = globalSz / grainSz;
    (std::abs(p_npes_prev - p_npes) <= std::abs(p_npes_next - p_npes))
        ? p_npes = p_npes_prev
        : p_npes = p_npes_next;

    if (p_npes > npes) p_npes = npes;
    // quick fix to enforce the npes>=2 for any given grain size.
    if (p_npes <= 1 && npes > 1) p_npes = 2;

    if (p_npes == npes) {
        MPI_Comm_dup(comm, &commActive);
        isActive = true;

    } else {
        // isActive=(rank*grainSz<globalSz);
        isActive = isRankSelected(npes, rank, p_npes);
        par::splitComm2way(isActive, &commActive, comm);
    }

    shrinkOrExpandOctree(tmpNodes, fluid::FLUID_LOAD_IMB_TOL,
                         DENDRO_DEFAULT_SF_K, isActive, commActive, comm);

    if (!isActive)
        if (tmpNodes.size() != 0)
            std::cout << " rank_g: " << rank << " isActive: " << isActive
                      << " f2O octants: " << tmpNodes.size() << std::endl;

    std::vector<ot::TreeNode> balOct;
    localSz = 0;
    if (isActive) {
        int rank_active, npes_active;

        MPI_Comm_size(commActive, &npes_active);
        MPI_Comm_rank(commActive, &rank_active);

        if (!rank_active)
            std::cout << "[MPI_COMM_SWITCH]: " << npes_active << std::endl;

        ot::TreeNode root(fluid::FLUID_DIM, fluid::FLUID_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        fluid::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes, tmpVec, tmpVec, tmpVec,
                                   fluid::FLUID_LOAD_IMB_TOL, m_uiMaxDepth,
                                   root, ROOT_ROTATION, 1, TS_REMOVE_DUPLICATES,
                                   fluid::FLUID_SPLIT_FIX, commActive);
        std::swap(tmpNodes, tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes, tmpVec, tmpVec, tmpVec,
                                   fluid::FLUID_LOAD_IMB_TOL, m_uiMaxDepth,
                                   root, ROOT_ROTATION, 1, TS_CONSTRUCT_OCTREE,
                                   fluid::FLUID_SPLIT_FIX, commActive);
        std::swap(tmpNodes, tmpVec);
        tmpVec.clear();

        fluid::timer::t_cons.stop();
        t_stat = fluid::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat, t_stat_g, 1, MPI_MIN, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 1, 1, MPI_SUM, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 2, 1, MPI_MAX, 0, commActive);
        t_stat_g[1] = t_stat_g[1] / (double)rank_active;

        localSz     = tmpNodes.size();
        par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, commActive);

        if (!rank_active)
            std::cout << GRN << "remove duplicates + octree construction (s): "
                      << t_stat_g[2] << NRM << std::endl;
        if (!rank_active)
            std::cout << GRN << " # const. octants: " << globalSz << NRM
                      << std::endl;

        fluid::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes, balOct, balOct, balOct,
                                   fluid::FLUID_LOAD_IMB_TOL, m_uiMaxDepth,
                                   root, ROOT_ROTATION, 1, TS_BALANCE_OCTREE,
                                   fluid::FLUID_SPLIT_FIX, commActive);
        tmpNodes.clear();

        fluid::timer::t_bal.stop();

        t_stat = fluid::timer::t_bal.seconds;
        par::Mpi_Reduce(&t_stat, t_stat_g, 1, MPI_MIN, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 1, 1, MPI_SUM, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 2, 1, MPI_MAX, 0, commActive);
        t_stat_g[1] = t_stat_g[1] / (double)rank_active;

        if (!rank_active)
            std::cout << GRN << " 2:1 balancing max (s): " << t_stat_g[2] << NRM
                      << std::endl;
        localSz = balOct.size();
    }
    MPI_Comm_free(&commActive);

    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz, &globalSz, 1, MPI_SUM, comm);
    if (!rank)
        std::cout << GRN << " balanced # octants : " << globalSz << NRM
                  << std::endl;

    fluid::timer::t_mesh.start();

    ot::Mesh* mesh =
        new ot::Mesh(balOct, 1, fluid::FLUID_ELE_ORDER, comm, true,
                     ot::SM_TYPE::FDM, fluid::FLUID_DENDRO_GRAIN_SZ,
                     fluid::FLUID_LOAD_IMB_TOL, fluid::FLUID_SPLIT_FIX);

    fluid::timer::t_mesh.stop();

    t_stat = fluid::timer::t_mesh.seconds;
    par::Mpi_Reduce(&t_stat, t_stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(&t_stat, t_stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(&t_stat, t_stat_g + 2, 1, MPI_MAX, 0, comm);
    t_stat_g[1] = t_stat_g[1] / (double)npes;

    mesh->setDomainBounds(
        Point(fluid::FLUID_GRID_MIN_X, fluid::FLUID_GRID_MIN_Y,
              fluid::FLUID_GRID_MIN_Z),
        Point(fluid::FLUID_GRID_MAX_X, fluid::FLUID_GRID_MAX_Y,
              fluid::FLUID_GRID_MAX_Z));
    localSz = mesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank)
        std::cout << GRN << " # of CG nodes (vertices) : " << globalSz << NRM
                  << std::endl;
    if (!rank) {
        std::cout << GRN << "Mesh generation time (max): " << t_stat_g[2] << NRM
                  << std::endl;
        std::cout << "\t" << GRN << " e2e (min,mean,max): " << "( "
                  << t_e2e_g[0] << "\t" << t_e2e_g[1] << "\t" << t_e2e_g[2]
                  << " )" << NRM << std::endl;
        std::cout << "\t" << GRN << " e2n (min,mean,max): " << "( "
                  << t_e2n_g[0] << "\t" << t_e2n_g[1] << "\t" << t_e2n_g[2]
                  << " )" << NRM << std::endl;
        std::cout << "\t" << GRN << " sm (min,mean,max): " << "( " << t_sm_g[0]
                  << "\t" << t_sm_g[1] << "\t" << t_sm_g[2] << " )" << NRM
                  << std::endl;
        std::cout << "\t" << GRN << " blk (min,mean,max): " << "( "
                  << t_blk_g[0] << "\t" << t_blk_g[1] << "\t" << t_blk_g[2]
                  << " )" << NRM << std::endl;
    }

    ode::solver::RK_FLUID rk_fluid(
        mesh, fluid::FLUID_RK_TIME_BEGIN, fluid::FLUID_RK_TIME_END,
        fluid::FLUID_RK45_TIME_STEP_SIZE, RKType::RK4);

    if (fluid::FLUID_RESTORE_SOLVER == 1)
        rk_fluid.restoreCheckPoint(fluid::FLUID_CHKPT_FILE_PREFIX.c_str(),
                                   comm);

    fluid::timer::t_rkSolve.start();
    rk_fluid.rkSolve();
    fluid::timer::t_rkSolve.stop();

    fluid::timer::total_runtime.stop();
    rk_fluid.freeMesh();
    // fluid::timer::profileInfo(fluid::FLUID_PROFILE_FILE_PREFIX.c_str(),mesh);

    MPI_Finalize();

    return 0;
}
