//
// Created by milinda on 7/25/17.
/**
*@author David Maughan
*Department of Physics, Utah State University
*@brief Header file for the Maxwell simulation.
*/
//

#include "maxwell.h"
#include "maxwellUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rk4maxwell.h"
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

    maxwell::timer::initFlops();

    maxwell::timer::total_runtime.start();


    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    maxwell::readParamFile(argv[1],comm);



    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;
// include file goes here to print parameters to screen.

    }

    _InitializeHcurve(maxwell::MAXWELL_DIM);
    m_uiMaxDepth=maxwell::MAXWELL_MAXDEPTH;

    if(maxwell::MAXWELL_NUM_VARS%maxwell::MAXWELL_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total MAXWELL_NUM_VARS: "<<maxwell::MAXWELL_NUM_VARS<<" is not divisable by MAXWELL_ASYNC_COMM_K: "<<maxwell::MAXWELL_ASYNC_COMM_K<<std::endl;
        exit(0);
    }

    maxwell::MAXWELL_RK45_TIME_STEP_SIZE=maxwell::MAXWELL_CFL_FACTOR*(maxwell::MAXWELL_COMPD_MAX[0]-maxwell::MAXWELL_COMPD_MIN[0])*(1.0/(double)(1u<<maxwell::MAXWELL_MAXDEPTH));

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){maxwell::initData(x,y,z,var);};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){maxwell::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=maxwell::MAXWELL_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<maxwell::MAXWELL_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=maxwell::VAR::U_ALPHA;
    varIndex[1]=maxwell::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    maxwell::timer::t_f2o.start();

    if(maxwell::MAXWELL_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(maxwell::MAXWELL_BLK_MIN_X,maxwell::MAXWELL_BLK_MIN_Y,maxwell::MAXWELL_BLK_MIN_Z);
        const Point pt_max(maxwell::MAXWELL_BLK_MAX_X,maxwell::MAXWELL_BLK_MAX_Y,maxwell::MAXWELL_BLK_MAX_Z);

        maxwell::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,maxwell::MAXWELL_NUM_VARS,varIndex,interpVars,tmpNodes,m_uiMaxDepth,maxwell::MAXWELL_WAVELET_TOL,maxwell::MAXWELL_ELE_ORDER,comm);

    }

    maxwell::timer::t_f2o.stop();

    t_stat=maxwell::timer::t_f2o.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=maxwell::MAXWELL_DENDRO_GRAIN_SZ;//DENDRO_DEFAULT_GRAIN_SZ;

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

    shrinkOrExpandOctree(tmpNodes,maxwell::MAXWELL_LOAD_IMB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

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

        ot::TreeNode root(maxwell::MAXWELL_DIM,maxwell::MAXWELL_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        maxwell::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,maxwell::MAXWELL_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,maxwell::MAXWELL_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,maxwell::MAXWELL_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,maxwell::MAXWELL_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        maxwell::timer::t_cons.stop();
        t_stat=maxwell::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


        maxwell::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,maxwell::MAXWELL_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,maxwell::MAXWELL_SPLIT_FIX,commActive);
        tmpNodes.clear();

        maxwell::timer::t_bal.stop();

        t_stat=maxwell::timer::t_bal.seconds;
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

    maxwell::timer::t_mesh.start();

    ot::Mesh * mesh=new ot::Mesh(balOct,1,maxwell::MAXWELL_ELE_ORDER,comm,maxwell::MAXWELL_DENDRO_GRAIN_SZ,maxwell::MAXWELL_LOAD_IMB_TOL,maxwell::MAXWELL_SPLIT_FIX);

    maxwell::timer::t_mesh.stop();

    t_stat=maxwell::timer::t_mesh.seconds;
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


    ode::solver::RK4_MAXWELL rk_maxwell(mesh,maxwell::MAXWELL_RK45_TIME_BEGIN,maxwell::MAXWELL_RK45_TIME_END,maxwell::MAXWELL_RK45_TIME_STEP_SIZE);
    //ode::solver::RK3_MAXWELL rk_maxwell(mesh,maxwell::MAXWELL_RK45_TIME_BEGIN,maxwell::MAXWELL_RK45_TIME_END,maxwell::MAXWELL_RK45_TIME_STEP_SIZE);

    if(maxwell::MAXWELL_RESTORE_SOLVER==1)
        rk_maxwell.restoreCheckPoint(maxwell::MAXWELL_CHKPT_FILE_PREFIX.c_str(),comm);

    maxwell::timer::t_rkSolve.start();
    rk_maxwell.rkSolve();
    maxwell::timer::t_rkSolve.stop();

    maxwell::timer::total_runtime.stop();
    rk_maxwell.freeMesh();
    //maxwell::timer::profileInfo(maxwell::MAXWELL_PROFILE_FILE_PREFIX.c_str(),mesh);
    //delete mesh;
    MPI_Finalize();

    return 0;
}
