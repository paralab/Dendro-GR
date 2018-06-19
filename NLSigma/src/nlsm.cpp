//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Header file for the GR simulation.
*/
//

#include "nlsm.h"
#include "nlsmUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rk4nlsm.h"
#include "octUtils.h"

int main (int argc, char** argv)
{


    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    nlsm::timer::initFlops();

    nlsm::timer::total_runtime.start();


    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    nlsm::readParamFile(argv[1],comm);



    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;

        std::cout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_DIM :"<<nlsm::NLSM_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_IO_OUTPUT_FREQ :"<<nlsm::NLSM_IO_OUTPUT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_REMESH_TEST_FREQ :"<<nlsm::NLSM_REMESH_TEST_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_CHECKPT_FREQ :"<<nlsm::NLSM_CHECKPT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_RESTORE_SOLVER :"<<nlsm::NLSM_RESTORE_SOLVER<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ENABLE_BLOCK_ADAPTIVITY :"<<nlsm::NLSM_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_VTU_FILE_PREFIX :"<<nlsm::NLSM_VTU_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_CHKPT_FILE_PREFIX :"<<nlsm::NLSM_CHKPT_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_PROFILE_FILE_PREFIX :"<<nlsm::NLSM_PROFILE_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_IO_OUTPUT_GAP :"<<nlsm::NLSM_IO_OUTPUT_GAP<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_DENDRO_GRAIN_SZ :"<<nlsm::NLSM_DENDRO_GRAIN_SZ<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ASYNC_COMM_K :"<<nlsm::NLSM_ASYNC_COMM_K<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_DENDRO_AMR_FAC :"<<nlsm::NLSM_DENDRO_AMR_FAC<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_WAVELET_TOL :"<<nlsm::NLSM_WAVELET_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_LOAD_IMB_TOL :"<<nlsm::NLSM_LOAD_IMB_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_RK45_TIME_BEGIN :"<<nlsm::NLSM_RK45_TIME_BEGIN<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_RK45_TIME_END :"<<nlsm::NLSM_RK45_TIME_END<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_RK45_TIME_STEP_SIZE :"<<nlsm::NLSM_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_RK45_DESIRED_TOL :"<<nlsm::NLSM_RK45_DESIRED_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_COMPD_MIN : ( :"<<nlsm::NLSM_COMPD_MIN[0]<<" ,"<<nlsm::NLSM_COMPD_MIN[1]<<","<<nlsm::NLSM_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_COMPD_MAX : ( :"<<nlsm::NLSM_COMPD_MAX[0]<<" ,"<<nlsm::NLSM_COMPD_MAX[1]<<","<<nlsm::NLSM_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_BLK_MIN : ( :"<<nlsm::NLSM_BLK_MIN_X<<" ,"<<nlsm::NLSM_BLK_MIN_Y<<","<<nlsm::NLSM_BLK_MIN_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_BLK_MAX : ( :"<<nlsm::NLSM_BLK_MAX_X<<" ,"<<nlsm::NLSM_BLK_MAX_Y<<","<<nlsm::NLSM_BLK_MAX_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_OCTREE_MIN : ( :"<<nlsm::NLSM_OCTREE_MIN[0]<<" ,"<<nlsm::NLSM_OCTREE_MIN[1]<<","<<nlsm::NLSM_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_OCTREE_MAX : ( :"<<nlsm::NLSM_OCTREE_MAX[0]<<" ,"<<nlsm::NLSM_OCTREE_MAX[1]<<","<<nlsm::NLSM_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<nlsm::KO_DISS_SIGMA<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_TYPE:"<<nlsm::NLSM_ID_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_AMP1:"<<nlsm::NLSM_ID_AMP1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_AMP2:"<<nlsm::NLSM_ID_AMP2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_DELTA1:"<<nlsm::NLSM_ID_DELTA1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_DELTA2:"<<nlsm::NLSM_ID_DELTA2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_XC1:"<<nlsm::NLSM_ID_XC1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_YC1:"<<nlsm::NLSM_ID_YC1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_ZC1:"<<nlsm::NLSM_ID_ZC1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_XC2:"<<nlsm::NLSM_ID_XC2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_YC2:"<<nlsm::NLSM_ID_YC2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_ZC2:"<<nlsm::NLSM_ID_ZC2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_EPSX1:"<<nlsm::NLSM_ID_EPSX1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_EPSY1:"<<nlsm::NLSM_ID_EPSY1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_EPSX2:"<<nlsm::NLSM_ID_EPSX2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_EPSY2:"<<nlsm::NLSM_ID_EPSY2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_R1:"<<nlsm::NLSM_ID_R1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_R2:"<<nlsm::NLSM_ID_R2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_NU1:"<<nlsm::NLSM_ID_NU1<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_NU2:"<<nlsm::NLSM_ID_NU2<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_ID_OMEGA:"<<nlsm::NLSM_ID_OMEGA<<NRM<<std::endl;

        std::cout<<YLW<<"\tNLSM_DIM :"<<nlsm::NLSM_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_MAXDEPTH :"<<nlsm::NLSM_MAXDEPTH<<NRM<<std::endl;

        std::cout<<YLW<<"\tNLSM_NUM_REFINE_VARS :"<<nlsm::NLSM_NUM_REFINE_VARS<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_REFINE_VARIABLE_INDICES :[";
        for(unsigned int i=0;i<nlsm::NLSM_NUM_REFINE_VARS-1;i++)
            std::cout<<nlsm::NLSM_REFINE_VARIABLE_INDICES[i]<<", ";
        std::cout<<nlsm::NLSM_REFINE_VARIABLE_INDICES[nlsm::NLSM_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tNLSM_NUM_EVOL_VARS_VTU_OUTPUT :"<<nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tNLSM_VTU_OUTPUT_EVOL_INDICES :[";
        for(unsigned int i=0;i<nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
            std::cout<<nlsm::NLSM_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
        std::cout<<nlsm::NLSM_VTU_OUTPUT_EVOL_INDICES[nlsm::NLSM_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;



    }

    _InitializeHcurve(nlsm::NLSM_DIM);
    m_uiMaxDepth=nlsm::NLSM_MAXDEPTH;

    if(nlsm::NLSM_NUM_VARS%nlsm::NLSM_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total NLSM_NUM_VARS: "<<nlsm::NLSM_NUM_VARS<<" is not divisable by NLSM_ASYNC_COMM_K: "<<nlsm::NLSM_ASYNC_COMM_K<<std::endl;
        exit(0);
    }

    nlsm::NLSM_RK45_TIME_STEP_SIZE=nlsm::NLSM_CFL_FACTOR*(nlsm::NLSM_COMPD_MAX[0]-nlsm::NLSM_COMPD_MIN[0])*(1.0/(double)(1u<<nlsm::NLSM_MAXDEPTH));

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){nlsm::initData(x,y,z,var);};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){nlsm::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=nlsm::NLSM_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<nlsm::NLSM_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=nlsm::VAR::U_ALPHA;
    varIndex[1]=nlsm::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    nlsm::timer::t_f2o.start();

    if(nlsm::NLSM_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(nlsm::NLSM_BLK_MIN_X,nlsm::NLSM_BLK_MIN_Y,nlsm::NLSM_BLK_MIN_Z);
        const Point pt_max(nlsm::NLSM_BLK_MAX_X,nlsm::NLSM_BLK_MAX_Y,nlsm::NLSM_BLK_MAX_Z);

        nlsm::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,nlsm::NLSM_NUM_VARS,varIndex,interpVars,tmpNodes,m_uiMaxDepth,nlsm::NLSM_WAVELET_TOL,nlsm::NLSM_ELE_ORDER,comm);

    }

    nlsm::timer::t_f2o.stop();

    t_stat=nlsm::timer::t_f2o.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=nlsm::NLSM_DENDRO_GRAIN_SZ;//DENDRO_DEFAULT_GRAIN_SZ;

    bool isActive;
    MPI_Comm commActive;
    const int p_npes_prev=binOp::getPrevHighestPowerOfTwo((globalSz/grainSz));
    const int p_npes_next=binOp::getNextHighestPowerOfTwo((globalSz/grainSz));

    int p_npes=globalSz/grainSz;
    (std::abs(p_npes_prev-p_npes)<=std::abs(p_npes_next-p_npes)) ? p_npes=p_npes_prev : p_npes=p_npes_next;

    if(p_npes>npes) p_npes=npes;
    // quick fix to enforce the npes>=2 for any given grain size.
    if(p_npes<=1 && npes>1) p_npes=2;

    if(p_npes==npes)
    {
        MPI_Comm_dup(comm,&commActive);
        isActive=true;

    }else
    {
        //isActive=(rank*grainSz<globalSz);
        isActive=isRankSelected(npes,rank,p_npes);
        par::splitComm2way(isActive,&commActive,comm);

    }

    shrinkOrExpandOctree(tmpNodes,nlsm::NLSM_LOAD_IMB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

    if(!isActive)
        if(tmpNodes.size()!=0)
            std::cout<<" rank_g: "<<rank<<" isActive: "<<isActive<<" f2O octants: "<<tmpNodes.size()<<std::endl;




    std::vector<ot::TreeNode> balOct;
    localSz=0;
    if(isActive)
    {

        int rank_active,npes_active;

        MPI_Comm_size(commActive,&npes_active);
        MPI_Comm_rank(commActive,&rank_active);

        if(!rank_active) std::cout<<"[MPI_COMM_SWITCH]: "<<npes_active<<std::endl;

        ot::TreeNode root(nlsm::NLSM_DIM,nlsm::NLSM_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        nlsm::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,nlsm::NLSM_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,nlsm::NLSM_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,nlsm::NLSM_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,nlsm::NLSM_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        nlsm::timer::t_cons.stop();
        t_stat=nlsm::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


        nlsm::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,nlsm::NLSM_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,nlsm::NLSM_SPLIT_FIX,commActive);
        tmpNodes.clear();

        nlsm::timer::t_bal.stop();

        t_stat=nlsm::timer::t_bal.seconds;
        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        if(!rank_active) std::cout<<GRN<<" 2:1 balancing max (s): "<<t_stat_g[2]<<NRM<<std::endl;
        localSz=balOct.size();


    }
    MPI_Comm_free(&commActive);

    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;

    nlsm::timer::t_mesh.start();

    ot::Mesh * mesh=new ot::Mesh(balOct,1,nlsm::NLSM_ELE_ORDER,comm,nlsm::NLSM_DENDRO_GRAIN_SZ,nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);

    nlsm::timer::t_mesh.stop();

    t_stat=nlsm::timer::t_mesh.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=mesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
    if(!rank) std::cout<<GRN<<" # of CG nodes (vertices) : "<<globalSz<<NRM<<std::endl;
    if(!rank)
    {
        std::cout<< GRN<<"Mesh generation time (max): "<<t_stat_g[2]<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;
    }


    ode::solver::RK4_NLSM rk_nlsm(mesh,nlsm::NLSM_RK45_TIME_BEGIN,nlsm::NLSM_RK45_TIME_END,nlsm::NLSM_RK45_TIME_STEP_SIZE);
    //ode::solver::RK3_NLSM rk_nlsm(mesh,nlsm::NLSM_RK45_TIME_BEGIN,nlsm::NLSM_RK45_TIME_END,nlsm::NLSM_RK45_TIME_STEP_SIZE);

    if(nlsm::NLSM_RESTORE_SOLVER==1)
        rk_nlsm.restoreCheckPoint(nlsm::NLSM_CHKPT_FILE_PREFIX.c_str(),comm);

    nlsm::timer::t_rkSolve.start();
    rk_nlsm.rkSolve();
    nlsm::timer::t_rkSolve.stop();

    nlsm::timer::total_runtime.stop();
    rk_nlsm.freeMesh();
    //nlsm::timer::profileInfo(nlsm::NLSM_PROFILE_FILE_PREFIX.c_str(),mesh);
    //delete mesh;
    MPI_Finalize();

    return 0;
}
