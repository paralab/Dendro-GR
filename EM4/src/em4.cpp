//
// Created by milinda on 7/25/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Header file for the GR simulation.
 */
//

#include "em4.h"

#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "em4Utils.h"
#include "mesh.h"
#include "mpi.h"
#include "octUtils.h"
#include "rk4em4.h"

int main(int argc, char** argv) {
    // feenableexcept( FE_DIVBYZERO );
    // feenableexcept( FE_DIVBYZERO | FE_INVALID );
    // feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

    if (argc < 2)
        std::cout << "Usage: " << argv[0] << " paramFile" << std::endl;

    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    em4::timer::initFlops();

    em4::timer::total_runtime.start();

    // 1 . read the parameter file.
    if (!rank) std::cout << " reading parameter file :" << argv[1] << std::endl;
    em4::readParamFile(argv[1], comm);

    if (rank == 1 || npes == 1) {
        std::cout << "parameters read: " << std::endl;

        std::cout << YLW << "\tnpes :" << npes << NRM << std::endl;
        std::cout << YLW << "\tEM4_ELE_ORDER :" << em4::EM4_ELE_ORDER << NRM
                  << std::endl;
        std::cout << YLW << "\tEM4_PADDING_WIDTH :" << em4::EM4_PADDING_WIDTH
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_DIM :" << em4::EM4_DIM << NRM << std::endl;
        std::cout << YLW << "\tEM4_IO_OUTPUT_FREQ :" << em4::EM4_IO_OUTPUT_FREQ
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM4_REMESH_TEST_FREQ :" << em4::EM4_REMESH_TEST_FREQ
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_CHECKPT_FREQ :" << em4::EM4_CHECKPT_FREQ
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_RESTORE_SOLVER :" << em4::EM4_RESTORE_SOLVER
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_ENABLE_BLOCK_ADAPTIVITY :"
                  << em4::EM4_ENABLE_BLOCK_ADAPTIVITY << NRM << std::endl;
        std::cout << YLW
                  << "\tEM4_VTU_FILE_PREFIX :" << em4::EM4_VTU_FILE_PREFIX
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM4_CHKPT_FILE_PREFIX :" << em4::EM4_CHKPT_FILE_PREFIX
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_PROFILE_FILE_PREFIX :"
                  << em4::EM4_PROFILE_FILE_PREFIX << NRM << std::endl;
        std::cout << YLW
                  << "\tEM4_VTU_Z_SLICE_ONLY :" << em4::EM4_VTU_Z_SLICE_ONLY
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_IO_OUTPUT_GAP :" << em4::EM4_IO_OUTPUT_GAP
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM4_DENDRO_GRAIN_SZ :" << em4::EM4_DENDRO_GRAIN_SZ
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_ASYNC_COMM_K :" << em4::EM4_ASYNC_COMM_K
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_DENDRO_AMR_FAC :" << em4::EM4_DENDRO_AMR_FAC
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_CFL_FACTOR:" << em4::EM4_CFL_FACTOR << NRM
                  << std::endl;
        std::cout << YLW << "\tEM4_WAVELET_TOL :" << em4::EM4_WAVELET_TOL << NRM
                  << std::endl;
        std::cout << YLW << "\tEM4_LOAD_IMB_TOL :" << em4::EM4_LOAD_IMB_TOL
                  << NRM << std::endl;
        std::cout << YLW
                  << "\tEM4_RK45_TIME_BEGIN :" << em4::EM4_RK45_TIME_BEGIN
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_RK45_TIME_END :" << em4::EM4_RK45_TIME_END
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_RK45_TIME_STEP_SIZE :"
                  << em4::EM4_RK45_TIME_STEP_SIZE << NRM << std::endl;
        std::cout << YLW
                  << "\tEM4_RK45_DESIRED_TOL :" << em4::EM4_RK45_DESIRED_TOL
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_COMPD_MIN : ( :" << em4::EM4_COMPD_MIN[0]
                  << " ," << em4::EM4_COMPD_MIN[1] << ","
                  << em4::EM4_COMPD_MIN[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tEM4_COMPD_MAX : ( :" << em4::EM4_COMPD_MAX[0]
                  << " ," << em4::EM4_COMPD_MAX[1] << ","
                  << em4::EM4_COMPD_MAX[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tEM4_BLK_MIN : ( :" << em4::EM4_BLK_MIN_X << " ,"
                  << em4::EM4_BLK_MIN_Y << "," << em4::EM4_BLK_MIN_Z << " )"
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_BLK_MAX : ( :" << em4::EM4_BLK_MAX_X << " ,"
                  << em4::EM4_BLK_MAX_Y << "," << em4::EM4_BLK_MAX_Z << " )"
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_OCTREE_MIN : ( :" << em4::EM4_OCTREE_MIN[0]
                  << " ," << em4::EM4_OCTREE_MIN[1] << ","
                  << em4::EM4_OCTREE_MIN[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tEM4_OCTREE_MAX : ( :" << em4::EM4_OCTREE_MAX[0]
                  << " ," << em4::EM4_OCTREE_MAX[1] << ","
                  << em4::EM4_OCTREE_MAX[2] << " )" << NRM << std::endl;
        std::cout << YLW << "\tKO_DISS_SIGMA :" << em4::KO_DISS_SIGMA << NRM
                  << std::endl;
        std::cout << YLW << "\tEM4_ID_TYPE:" << em4::EM4_ID_TYPE << NRM
                  << std::endl;
        std::cout << YLW << "\tEM4_ID_AMP1:" << em4::EM4_ID_AMP1 << NRM
                  << std::endl;
        std::cout << YLW << "\tEM4_ID_LAMBDA1:" << em4::EM4_ID_LAMBDA1 << NRM
                  << std::endl;
        std::cout << YLW << "\tEM4_ID_AMP2:" << em4::EM4_ID_AMP2 << NRM
                  << std::endl;
        // std::cout<<YLW<<"\tEM4_ID_DELTA1:"<<em4::EM4_ID_DELTA1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_DELTA2:"<<em4::EM4_ID_DELTA2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_XC1:"<<em4::EM4_ID_XC1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_YC1:"<<em4::EM4_ID_YC1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_ZC1:"<<em4::EM4_ID_ZC1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_XC2:"<<em4::EM4_ID_XC2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_YC2:"<<em4::EM4_ID_YC2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_ZC2:"<<em4::EM4_ID_ZC2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_EPSX1:"<<em4::EM4_ID_EPSX1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_EPSY1:"<<em4::EM4_ID_EPSY1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_EPSZ1:"<<em4::EM4_ID_EPSY1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_EPSX2:"<<em4::EM4_ID_EPSX2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_EPSY2:"<<em4::EM4_ID_EPSY2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_EPSZ2:"<<em4::EM4_ID_EPSY2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_R1:"<<em4::EM4_ID_R1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_R2:"<<em4::EM4_ID_R2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_NU1:"<<em4::EM4_ID_NU1<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_NU2:"<<em4::EM4_ID_NU2<<NRM<<std::endl;
        // std::cout<<YLW<<"\tEM4_ID_OMEGA:"<<em4::EM4_ID_OMEGA<<NRM<<std::endl;

        // std::cout<<YLW<<"\tEM4_DIM :"<<em4::EM4_DIM<<NRM<<std::endl;
        std::cout << YLW << "\tEM4_MAXDEPTH :" << em4::EM4_MAXDEPTH << NRM
                  << std::endl;

        std::cout << YLW
                  << "\tEM4_NUM_REFINE_VARS :" << em4::EM4_NUM_REFINE_VARS
                  << NRM << std::endl;
        std::cout << YLW << "\tEM4_REFINE_VARIABLE_INDICES :[";
        for (unsigned int i = 0; i < em4::EM4_NUM_REFINE_VARS - 1; i++)
            std::cout << em4::EM4_REFINE_VARIABLE_INDICES[i] << ", ";
        std::cout
            << em4::EM4_REFINE_VARIABLE_INDICES[em4::EM4_NUM_REFINE_VARS - 1]
            << "]" << NRM << std::endl;

        std::cout << YLW
                  << "\tEM4_REFINEMENT_MODE :" << em4::EM4_REFINEMENT_MODE
                  << NRM << std::endl;

        std::cout << YLW << "\tEM4_NUM_EVOL_VARS_VTU_OUTPUT :"
                  << em4::EM4_NUM_EVOL_VARS_VTU_OUTPUT << NRM << std::endl;
        std::cout << YLW << "\tEM4_VTU_OUTPUT_EVOL_INDICES :[";
        for (unsigned int i = 0; i < em4::EM4_NUM_EVOL_VARS_VTU_OUTPUT - 1; i++)
            std::cout << em4::EM4_VTU_OUTPUT_EVOL_INDICES[i] << ", ";
        std::cout << em4::EM4_VTU_OUTPUT_EVOL_INDICES
                         [em4::EM4_NUM_EVOL_VARS_VTU_OUTPUT - 1]
                  << "]" << NRM << std::endl;

#ifdef EM4_USE_4TH_ORDER_DERIVS
        std::cout << "Using 4th order FD stencils. " << std::endl;
#endif

#ifdef EM4_USE_6TH_ORDER_DERIVS
        std::cout << "Using 6th order FD stencils. " << std::endl;
#endif

#ifdef EM4_USE_8TH_ORDER_DERIVS
        std::cout << "Using 8th order FD stencils. " << std::endl;
#endif
    }

    _InitializeHcurve(em4::EM4_DIM);
    m_uiMaxDepth = em4::EM4_MAXDEPTH;

    if (em4::EM4_NUM_VARS % em4::EM4_ASYNC_COMM_K != 0) {
        if (!rank)
            std::cout << "[overlap communication error]: total EM4_NUM_VARS: "
                      << em4::EM4_NUM_VARS
                      << " is not divisable by EM4_ASYNC_COMM_K: "
                      << em4::EM4_ASYNC_COMM_K << std::endl;
        exit(0);
    }

    // 2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            em4::initData(x, y, z, var);
        };
    // std::function<void(double,double,double,double,double*)> u_x_t=[](double
    // x,double y,double z,double
    // t,double*var){em4::analyticalSol(x,y,z,t,var);};

    const unsigned int interpVars = em4::EM4_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < em4::EM4_NUM_VARS; i++) varIndex[i] = i;

    DendroIntL localSz, globalSz;
    double t_stat;
    double t_stat_g[3];

    em4::timer::t_f2o.start();

    if (em4::EM4_ENABLE_BLOCK_ADAPTIVITY) {
        if (!rank)
            std::cout << YLW << "Using block adaptive mesh. AMR disabled "
                      << NRM << std::endl;
        const Point pt_min(em4::EM4_BLK_MIN_X, em4::EM4_BLK_MIN_Y,
                           em4::EM4_BLK_MIN_Z);
        const Point pt_max(em4::EM4_BLK_MAX_X, em4::EM4_BLK_MAX_Y,
                           em4::EM4_BLK_MAX_Z);

        em4::blockAdaptiveOctree(
            tmpNodes, pt_min, pt_max,
            m_uiMaxDepth - (binOp::fastLog2(em4::EM4_ELE_ORDER)), m_uiMaxDepth,
            comm);
    } else {
        if (!rank)
            std::cout << YLW << "Using function2Octree. AMR enabled " << NRM
                      << std::endl;
        function2Octree(f_init, em4::EM4_NUM_VARS,
                        em4::EM4_REFINE_VARIABLE_INDICES,
                        em4::EM4_NUM_REFINE_VARS, tmpNodes, m_uiMaxDepth,
                        em4::EM4_WAVELET_TOL, em4::EM4_ELE_ORDER, comm);
        // std::cout<<"f2o else end"<<std::endl;
    }

    em4::timer::t_f2o.stop();

    t_stat = em4::timer::t_f2o.seconds;
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
        em4::EM4_DENDRO_GRAIN_SZ;  // DENDRO_DEFAULT_GRAIN_SZ;

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

    shrinkOrExpandOctree(tmpNodes, em4::EM4_LOAD_IMB_TOL, DENDRO_DEFAULT_SF_K,
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

        ot::TreeNode root(em4::EM4_DIM, em4::EM4_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        em4::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes, tmpVec, tmpVec, tmpVec,
                                   em4::EM4_LOAD_IMB_TOL, m_uiMaxDepth, root,
                                   ROOT_ROTATION, 1, TS_REMOVE_DUPLICATES,
                                   em4::EM4_SPLIT_FIX, commActive);
        std::swap(tmpNodes, tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes, tmpVec, tmpVec, tmpVec,
                                   em4::EM4_LOAD_IMB_TOL, m_uiMaxDepth, root,
                                   ROOT_ROTATION, 1, TS_CONSTRUCT_OCTREE,
                                   em4::EM4_SPLIT_FIX, commActive);
        std::swap(tmpNodes, tmpVec);
        tmpVec.clear();

        em4::timer::t_cons.stop();
        t_stat = em4::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat, t_stat_g, 1, MPI_MIN, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 1, 1, MPI_SUM, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 2, 1, MPI_MAX, 0, commActive);
        t_stat_g[1] = t_stat_g[1] / (double)npes_active;

        localSz     = tmpNodes.size();
        par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, commActive);

        if (!rank_active)
            std::cout << GRN << "remove duplicates + octree construction (s): "
                      << t_stat_g[2] << NRM << std::endl;
        if (!rank_active)
            std::cout << GRN << " # const. octants: " << globalSz << NRM
                      << std::endl;

        em4::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes, balOct, balOct, balOct,
                                   em4::EM4_LOAD_IMB_TOL, m_uiMaxDepth, root,
                                   ROOT_ROTATION, 1, TS_BALANCE_OCTREE,
                                   em4::EM4_SPLIT_FIX, commActive);
        tmpNodes.clear();

        em4::timer::t_bal.stop();

        t_stat = em4::timer::t_bal.seconds;
        par::Mpi_Reduce(&t_stat, t_stat_g, 1, MPI_MIN, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 1, 1, MPI_SUM, 0, commActive);
        par::Mpi_Reduce(&t_stat, t_stat_g + 2, 1, MPI_MAX, 0, commActive);
        t_stat_g[1] = t_stat_g[1] / (double)npes_active;

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

    em4::timer::t_mesh.start();

    ot::Mesh* mesh = new ot::Mesh(balOct, 1, em4::EM4_ELE_ORDER, comm, true,
                                  ot::SM_TYPE::FDM, em4::EM4_DENDRO_GRAIN_SZ,
                                  em4::EM4_LOAD_IMB_TOL, em4::EM4_SPLIT_FIX);
    mesh->setDomainBounds(
        Point(em4::EM4_GRID_MIN_X, em4::EM4_GRID_MIN_Y, em4::EM4_GRID_MIN_Z),
        Point(em4::EM4_GRID_MAX_X, em4::EM4_GRID_MAX_Y, em4::EM4_GRID_MAX_Z));
    em4::timer::t_mesh.stop();

    t_stat = em4::timer::t_mesh.seconds;
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
    em4::EM4_RK45_TIME_STEP_SIZE =
        em4::EM4_CFL_FACTOR *
        ((em4::EM4_COMPD_MAX[0] - em4::EM4_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)em4::EM4_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    par::Mpi_Bcast(&em4::EM4_RK45_TIME_STEP_SIZE, 1, 0, comm);
    // std::cout<<" lmin: "<<lmin<<" lmax: "<<lmax<<std::endl;

    ode::solver::RK4_EM4 rk_em4(mesh, em4::EM4_RK45_TIME_BEGIN,
                                em4::EM4_RK45_TIME_END,
                                em4::EM4_RK45_TIME_STEP_SIZE);
    // ode::solver::RK3_EM4
    // rk_em4(mesh,em4::EM4_RK45_TIME_BEGIN,em4::EM4_RK45_TIME_END,em4::EM4_RK45_TIME_STEP_SIZE);

    if (em4::EM4_RESTORE_SOLVER == 1)
        rk_em4.restoreCheckPoint(em4::EM4_CHKPT_FILE_PREFIX.c_str(), comm);

    em4::timer::t_rkSolve.start();
    rk_em4.rkSolve();
    em4::timer::t_rkSolve.stop();

    em4::timer::total_runtime.stop();
    rk_em4.freeMesh();
    // em4::timer::profileInfo(em4::EM4_PROFILE_FILE_PREFIX.c_str(),mesh);
    // delete mesh;
    MPI_Finalize();

    return 0;
}
