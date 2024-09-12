//
// Created by milinda on 7/25/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Header file for the GR simulation.
 */
//

#include "em1.h"

#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "em1Utils.h"
#include "mesh.h"
#include "mpi.h"
#include "octUtils.h"
#include "rk4em1.h"

int main(int argc, char** argv) {
    if (argc < 2)
        std::cout << "Usage: " << argv[0] << " paramFile" << std::endl;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    em1::timer::initFlops();

    em1::timer::total_runtime.start();

    // 1 . read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    em1::readParamFile(argv[1], comm);

    if (rank == 1 || npes == 1) {
        std::cout << "parameters read: " << std::endl;

        std::cout << YLW << "\tnpes :" << npes << NRM << std::endl;
        std::cout << YLW << "\tEM1_ELE_ORDER :" << em1::EM1_ELE_ORDER << NRM
                  << std::endl;
        std::cout << YLW << "\tEM1_PADDING_WIDTH :" << em1::EM1_PADDING_WIDTH
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_DIM :" << em1::EM1_DIM << NRM << std::endl;
        std::cout << YLW << "\tEM1_IO_OUTPUT_FREQ :" << em1::EM1_IO_OUTPUT_FREQ
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM1_REMESH_TEST_FREQ :" << em1::EM1_REMESH_TEST_FREQ
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_CHECKPT_FREQ :" << em1::EM1_CHECKPT_FREQ
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_RESTORE_SOLVER :" << em1::EM1_RESTORE_SOLVER
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_ENABLE_BLOCK_ADAPTIVITY :"
                  << em1::EM1_ENABLE_BLOCK_ADAPTIVITY << NRM << std::endl;
        std::cout << YLW
                  << "\tEM1_VTU_FILE_PREFIX :" << em1::EM1_VTU_FILE_PREFIX
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM1_CHKPT_FILE_PREFIX :" << em1::EM1_CHKPT_FILE_PREFIX
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_PROFILE_FILE_PREFIX :"
                  << em1::EM1_PROFILE_FILE_PREFIX << NRM << std::endl;
        std::cout << YLW
                  << "\tEM1_VTU_Z_SLICE_ONLY :" << em1::EM1_VTU_Z_SLICE_ONLY
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_IO_OUTPUT_GAP :" << em1::EM1_IO_OUTPUT_GAP
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM1_DENDRO_GRAIN_SZ :" << em1::EM1_DENDRO_GRAIN_SZ
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_ASYNC_COMM_K :" << em1::EM1_ASYNC_COMM_K
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_DENDRO_AMR_FAC :" << em1::EM1_DENDRO_AMR_FAC
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_CFL_FACTOR:" << em1::EM1_CFL_FACTOR << NRM
                  << std::endl;
        std::cout << YLW << "\tEM1_WAVELET_TOL :" << em1::EM1_WAVELET_TOL << NRM
                  << std::endl;
        std::cout << YLW << "\tEM1_LOAD_IMB_TOL :" << em1::EM1_LOAD_IMB_TOL
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM1_RK45_TIME_BEGIN :" << em1::EM1_RK45_TIME_BEGIN
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_RK45_TIME_END :" << em1::EM1_RK45_TIME_END
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_RK45_TIME_STEP_SIZE :"
                  << em1::EM1_RK45_TIME_STEP_SIZE << NRM << std::endl;
        std::cout << YLW
                  << "\tEM1_RK45_DESIRED_TOL :" << em1::EM1_RK45_DESIRED_TOL
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_COMPD_MIN : ( :" << em1::EM1_COMPD_MIN[0]
                  << " ," << em1::EM1_COMPD_MIN[1] << ","
                  << em1::EM1_COMPD_MIN[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tEM1_COMPD_MAX : ( :" << em1::EM1_COMPD_MAX[0]
                  << " ," << em1::EM1_COMPD_MAX[1] << ","
                  << em1::EM1_COMPD_MAX[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tEM1_BLK_MIN : ( :" << em1::EM1_BLK_MIN_X << " ,"
                  << em1::EM1_BLK_MIN_Y << "," << em1::EM1_BLK_MIN_Z << " )"
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_BLK_MAX : ( :" << em1::EM1_BLK_MAX_X << " ,"
                  << em1::EM1_BLK_MAX_Y << "," << em1::EM1_BLK_MAX_Z << " )"
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_OCTREE_MIN : ( :" << em1::EM1_OCTREE_MIN[0]
                  << " ," << em1::EM1_OCTREE_MIN[1] << ","
                  << em1::EM1_OCTREE_MIN[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tEM1_OCTREE_MAX : ( :" << em1::EM1_OCTREE_MAX[0]
                  << " ," << em1::EM1_OCTREE_MAX[1] << ","
                  << em1::EM1_OCTREE_MAX[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tKO_DISS_SIGMA :" << em1::KO_DISS_SIGMA << NRM
                  << std::endl;
        std::cout << YLW << "\tEM1_ID_TYPE:" << em1::EM1_ID_TYPE << NRM
                  << std::endl;
        std::cout << YLW << "\tEM1_ID_AMP1:" << em1::EM1_ID_AMP1 << NRM
                  << std::endl;
        std::cout << YLW << "\tEM1_ID_LAMBDA1:" << em1::EM1_ID_LAMBDA1 << NRM
                  << std::endl;
        std::cout << YLW << "\tEM1_ID_AMP2:" << em1::EM1_ID_AMP2 << NRM
                  << std::endl;
        // std::cout<<YLW<<"\tEM1_ID_DELTA1:"<<em1::EM1_ID_DELTA1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_DELTA2:"<<em1::EM1_ID_DELTA2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_XC1:"<<em1::EM1_ID_XC1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_YC1:"<<em1::EM1_ID_YC1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_ZC1:"<<em1::EM1_ID_ZC1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_XC2:"<<em1::EM1_ID_XC2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_YC2:"<<em1::EM1_ID_YC2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_ZC2:"<<em1::EM1_ID_ZC2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_EPSX1:"<<em1::EM1_ID_EPSX1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_EPSY1:"<<em1::EM1_ID_EPSY1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_EPSZ1:"<<em1::EM1_ID_EPSY1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_EPSX2:"<<em1::EM1_ID_EPSX2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_EPSY2:"<<em1::EM1_ID_EPSY2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_EPSZ2:"<<em1::EM1_ID_EPSY2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_R1:"<<em1::EM1_ID_R1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_R2:"<<em1::EM1_ID_R2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_NU1:"<<em1::EM1_ID_NU1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_NU2:"<<em1::EM1_ID_NU2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM1_ID_OMEGA:"<<em1::EM1_ID_OMEGA<<NRM<<std::endl;

        // std::cout<<YLW<<"\tEM1_DIM :"<<em1::EM1_DIM<<NRM<<std::endl;
        std::cout << YLW << "\tEM1_MAXDEPTH :" << em1::EM1_MAXDEPTH << NRM
                  << std::endl;

        std::cout << YLW
                  << "\tEM1_NUM_REFINE_VARS :" << em1::EM1_NUM_REFINE_VARS
                  << NRM << std::endl;
        std::cout << YLW << "\tEM1_REFINE_VARIABLE_INDICES :[";
        for (unsigned int i = 0; i < em1::EM1_NUM_REFINE_VARS - 1; i++)
            std::cout << em1::EM1_REFINE_VARIABLE_INDICES[i] << ", ";
        std::cout
            << em1::EM1_REFINE_VARIABLE_INDICES[em1::EM1_NUM_REFINE_VARS - 1]
            << "]" << NRM << std::endl;

        std::cout << YLW
                  << "\tEM1_REFINEMENT_MODE :" << em1::EM1_REFINEMENT_MODE
                  << NRM << std::endl;

        std::cout << YLW << "\tEM1_NUM_EVOL_VARS_VTU_OUTPUT :"
                  << em1::EM1_NUM_EVOL_VARS_VTU_OUTPUT << NRM << std::endl;
        std::cout << YLW << "\tEM1_VTU_OUTPUT_EVOL_INDICES :[";
        for (unsigned int i = 0; i < em1::EM1_NUM_EVOL_VARS_VTU_OUTPUT - 1; i++)
            std::cout << em1::EM1_VTU_OUTPUT_EVOL_INDICES[i] << ", ";
        std::cout << em1::EM1_VTU_OUTPUT_EVOL_INDICES
                         [em1::EM1_NUM_EVOL_VARS_VTU_OUTPUT - 1]
                  << "]" << NRM << std::endl;

#ifdef EM1_USE_4TH_ORDER_DERIVS
        std::cout << "Using 4th order FD stencils. " << std::endl;
#endif

#ifdef EM1_USE_6TH_ORDER_DERIVS
        std::cout << "Using 6th order FD stencils. " << std::endl;
#endif
    }

    _InitializeHcurve(em1::EM1_DIM);
    m_uiMaxDepth = em1::EM1_MAXDEPTH;

    if (em1::EM1_NUM_VARS % em1::EM1_ASYNC_COMM_K != 0) {
        if (!rank)
            std::cout << "[overlap communication error]: total EM1_NUM_VARS: "
                      << em1::EM1_NUM_VARS
                      << " is not divisable by EM1_ASYNC_COMM_K: "
                      << em1::EM1_ASYNC_COMM_K << std::endl;
        exit(0);
    }

    // 2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            em1::initData(x, y, z, var);
        };
    // std::function<void(double,double,double,double,double*)> u_x_t=[](double
    // x,double y,double z,double
    // t,double*var){em1::analyticalSol(x,y,z,t,var);};

    const unsigned int interpVars = em1::EM1_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < em1::EM1_NUM_VARS; i++) varIndex[i] = i;

    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    em1::timer::t_f2o.start();

    if (em1::EM1_ENABLE_BLOCK_ADAPTIVITY) {
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled "
                      << NRM << std::endl;
        const Point pt_min(em1::EM1_BLK_MIN_X, em1::EM1_BLK_MIN_Y,
                           em1::EM1_BLK_MIN_Z);
        const Point pt_max(em1::EM1_BLK_MAX_X, em1::EM1_BLK_MAX_Y,
                           em1::EM1_BLK_MAX_Z);

        em1::blockAdaptiveOctree(
            tmpNodes, pt_min, pt_max,
            m_uiMaxDepth - (binOp::fastLog2(em1::EM1_ELE_ORDER)), m_uiMaxDepth,
            comm);
    } else {
        if (!rank)
            std::cout << YLW << "Using function2Octree. AMR enabled " << NRM
                      << std::endl;
        function2Octree(f_init, em1::EM1_NUM_VARS,
                        em1::EM1_REFINE_VARIABLE_INDICES,
                        em1::EM1_NUM_REFINE_VARS, tmpNodes, m_uiMaxDepth,
                        em1::EM1_WAVELET_TOL, em1::EM1_ELE_ORDER, comm);
        // std::cout<<"f2o else end"<<std::endl;
    }

    em1::timer::t_f2o.stop();

    t_stat = em1::timer::t_f2o.seconds;
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
        em1::EM1_DENDRO_GRAIN_SZ;  // DENDRO_DEFAULT_GRAIN_SZ;

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

    shrinkOrExpandOctree(tmpNodes, em1::EM1_LOAD_IMB_TOL, DENDRO_DEFAULT_SF_K,
                         isActive, commActive, comm);

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

        ot::TreeNode root(em1::EM1_DIM, em1::EM1_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        em1::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes, tmpVec, tmpVec, tmpVec,
                                   em1::EM1_LOAD_IMB_TOL, m_uiMaxDepth, root,
                                   ROOT_ROTATION, 1, TS_REMOVE_DUPLICATES,
                                   em1::EM1_SPLIT_FIX, commActive);
        std::swap(tmpNodes, tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes, tmpVec, tmpVec, tmpVec,
                                   em1::EM1_LOAD_IMB_TOL, m_uiMaxDepth, root,
                                   ROOT_ROTATION, 1, TS_CONSTRUCT_OCTREE,
                                   em1::EM1_SPLIT_FIX, commActive);
        std::swap(tmpNodes, tmpVec);
        tmpVec.clear();

        em1::timer::t_cons.stop();
        t_stat = em1::timer::t_cons.seconds;

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

        em1::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes, balOct, balOct, balOct,
                                   em1::EM1_LOAD_IMB_TOL, m_uiMaxDepth, root,
                                   ROOT_ROTATION, 1, TS_BALANCE_OCTREE,
                                   em1::EM1_SPLIT_FIX, commActive);
        tmpNodes.clear();

        em1::timer::t_bal.stop();

        t_stat = em1::timer::t_bal.seconds;
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

    em1::timer::t_mesh.start();

    ot::Mesh* mesh = new ot::Mesh(balOct, 1, em1::EM1_ELE_ORDER, comm, true,
                                  ot::SM_TYPE::FDM, em1::EM1_DENDRO_GRAIN_SZ,
                                  em1::EM1_LOAD_IMB_TOL, em1::EM1_SPLIT_FIX);
    mesh->setDomainBounds(
        Point(em1::EM1_GRID_MIN_X, em1::EM1_GRID_MIN_Y, em1::EM1_GRID_MIN_Z),
        Point(em1::EM1_GRID_MAX_X, em1::EM1_GRID_MAX_Y, em1::EM1_GRID_MAX_Z));
    em1::timer::t_mesh.stop();

    t_stat = em1::timer::t_mesh.seconds;
    par::Mpi_Reduce(&t_stat, t_stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(&t_stat, t_stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(&t_stat, t_stat_g + 2, 1, MPI_MAX, 0, comm);
    t_stat_g[1] = t_stat_g[1] / (double)npes;

    localSz     = mesh->getNumLocalMeshNodes();
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

    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);
    em1::EM1_RK45_TIME_STEP_SIZE =
        em1::EM1_CFL_FACTOR *
        ((em1::EM1_COMPD_MAX[0] - em1::EM1_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)em1::EM1_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    par::Mpi_Bcast(&em1::EM1_RK45_TIME_STEP_SIZE, 1, 0, comm);
    // std::cout<<" lmin: "<<lmin<<" lmax: "<<lmax<<std::endl;

    ode::solver::RK4_EM1 rk_em1(mesh, em1::EM1_RK45_TIME_BEGIN,
                                em1::EM1_RK45_TIME_END,
                                em1::EM1_RK45_TIME_STEP_SIZE);
    // ode::solver::RK3_EM1
    // rk_em1(mesh,em1::EM1_RK45_TIME_BEGIN,em1::EM1_RK45_TIME_END,em1::EM1_RK45_TIME_STEP_SIZE);

    if (em1::EM1_RESTORE_SOLVER == 1)
        rk_em1.restoreCheckPoint(em1::EM1_CHKPT_FILE_PREFIX.c_str(), comm);

    em1::timer::t_rkSolve.start();
    rk_em1.rkSolve();
    em1::timer::t_rkSolve.stop();

    em1::timer::total_runtime.stop();
    rk_em1.freeMesh();
    // em1::timer::profileInfo(em1::EM1_PROFILE_FILE_PREFIX.c_str(),mesh);
    // delete mesh;
    MPI_Finalize();

    return 0;
}
