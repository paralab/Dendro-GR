//
// Created by milinda on 11/30/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains a basic laplace or poisson equation and solve that using FEM
 *
 * \Delta u = f
 * with zero Diriclet boundary conditions.
 *
*/
//


#include "TreeNode.h"
#include "mpi.h"
#include "genPts_par.h"
#include "sfcSort.h"
#include "mesh.h"
#include "dendro.h"
#include "dendroIO.h"
#include "octUtils.h"
#include "functional"
#include "fdCoefficient.h"
#include "stencil.h"
#include "rkTransport.h"
#include "refel.h"
#include "operators.h"
#include "cg.h"


int main (int argc, char** argv) {


    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (argc < 4) {
        if (!rank)
            std::cout << "Usage: " << argv[0]
                      << " maxDepth wavelet_tol partition_tol eleOrder "
                      << std::endl;
        return 0;
    }

    m_uiMaxDepth = atoi(argv[1]);
    double wavelet_tol = atof(argv[2]);
    double partition_tol = atof(argv[3]);
    unsigned int eOrder = atoi(argv[4]);


    double tBegin = 0, tEnd = 10, th = 0.01;

    if (!rank) {
        std::cout << YLW << "maxDepth: " << m_uiMaxDepth << NRM << std::endl;
        std::cout << YLW << "wavelet_tol: " << wavelet_tol << NRM << std::endl;
        std::cout << YLW << "partition_tol: " << partition_tol << NRM << std::endl;
        std::cout << YLW << "eleOrder: " << eOrder << NRM << std::endl;

    }

    _InitializeHcurve(m_uiDim);

    RefElement refEl(m_uiDim,eOrder);

    // function that we need to interpolate.
    fem::domain::grid_min=Point(0, 0, 0);
    fem::domain::grid_max=Point((1u << m_uiMaxDepth), (1u << m_uiMaxDepth), (1u << m_uiMaxDepth));

    fem::domain::domain_min=Point(-0.5,-0.5,-0.5);
    fem::domain::domain_max=Point(0.5,0.5,0.5);

    fem::domain::Rg_x=(fem::domain::grid_max.x()-fem::domain::grid_min.x());
    fem::domain::Rg_y=(fem::domain::grid_max.y()-fem::domain::grid_min.y());
    fem::domain::Rg_z=(fem::domain::grid_max.z()-fem::domain::grid_min.z());

    fem::domain::Rd_x=(fem::domain::domain_max.x()-fem::domain::domain_min.x());
    fem::domain::Rd_y=(fem::domain::domain_max.y()-fem::domain::domain_min.y());
    fem::domain::Rd_z=(fem::domain::domain_max.z()-fem::domain::domain_min.z());

    const Point d_min=fem::domain::domain_min;
    const Point d_max=fem::domain::domain_max;

    const Point g_min=fem::domain::grid_min;
    const Point g_max=fem::domain::grid_max;

    const double Rg_x=fem::domain::Rg_x;
    const double Rg_y=fem::domain::Rg_y;
    const double Rg_z=fem::domain::Rg_z;

    const double Rd_x=fem::domain::Rd_x;
    const double Rd_y=fem::domain::Rd_y;
    const double Rd_z=fem::domain::Rd_z;


    std::function<double(double,double,double)> f_rhs =[d_min,d_max,g_min,g_max,Rg_x,Rg_y,Rg_z,Rd_x,Rd_y,Rd_z](const double x,const double y,const double z){
        return (-12*M_PI*M_PI*sin(2*M_PI*(((x-g_min.x())/(Rg_x))*(Rd_x)+d_min.x()))*sin(2*M_PI*(((y-g_min.y())/(Rg_y))*(Rd_y)+d_min.y()))*sin(2*M_PI*(((z-g_min.z())/(Rg_z))*(Rd_z)+d_min.z())));
    };

    std::vector<ot::TreeNode> tmpNodes;
    DendroIntL localSz,globalSz;
    double t_stat_g[3];
    double t_stat;

    double t1=MPI_Wtime();
    function2Octree(f_rhs, tmpNodes, m_uiMaxDepth, wavelet_tol, eOrder, comm);
    double t2=MPI_Wtime();

    t_stat=t2-t1;
    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=DENDRO_DEFAULT_GRAIN_SZ;

    bool isActive;
    MPI_Comm commActive;

    if((globalSz/grainSz)>=npes)
    {
        MPI_Comm_dup(comm,&commActive);
        isActive=true;

    }else
    {
        isActive=(rank*grainSz<globalSz);
        par::splitComm2way(isActive,&commActive,comm);

    }

    shrinkOrExpandOctree(tmpNodes,partition_tol,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

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

        ot::TreeNode root(m_uiDim,m_uiMaxDepth);
        std::vector<ot::TreeNode> tmpVec;
        t1=MPI_Wtime();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,partition_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,DENDRO_DEFAULT_SF_K,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,partition_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        t2=MPI_Wtime();
        t_stat=t2-t1;

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


        t1=MPI_Wtime();

        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,partition_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
        tmpNodes.clear();

        t2=MPI_Wtime();

        t_stat=t2-t1;
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


    t1=MPI_Wtime();
    ot::Mesh * mesh=new ot::Mesh(balOct,1,eOrder,comm,true,ot::SM_TYPE::FEM_CG,grainSz,partition_tol);
    t2=MPI_Wtime();

    t_stat=t2-t1;
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

    const unsigned int dof=mesh->getDegOfFreedom();
    const unsigned int nodeLocalBegin=mesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd=mesh->getNodeLocalEnd();
    const unsigned int local_dof=nodeLocalEnd-nodeLocalBegin;
    int maxIter=1000;
    double cg_tol=1e-5;


    double * rhs=mesh->createVector<double>(f_rhs);
    double * Mrhs=mesh->createVector<double>(0);
    double * u=mesh->createVector<double>(0);

    mesh->performGhostExchange(rhs);
    mesh->performGhostExchange(Mrhs);
    mesh->performGhostExchange(u);

    MPI_Comm activeComm;
    int activeRank;
    int activeNpes;

/*    for(unsigned int i=nodeLocalBegin;i<nodeLocalEnd;i++)
        if(isnanf(rhs[i]))std::cout<<"rank: "<<rank<<" node: "<<i<<" value: "<<rhs[i]<<std::endl;*/

    if(mesh->isActive())
    {
        activeComm=mesh->getMPICommunicator();
        activeRank=mesh->getMPIRank();
        activeNpes=mesh->getMPICommSize();

        const unsigned int nPe=mesh->getNumNodesPerElement();
        // allocate the work space vars for the matvec
        fem::operators::poisson::allocateWorkSpace(nPe);

        double vec_min,vec_max;
        fem::operators::poisson::computeRHS(mesh,rhs,Mrhs);

        vec_min=vecMin(rhs+nodeLocalBegin,local_dof,activeComm);
        vec_max=vecMax(rhs+nodeLocalBegin,local_dof,activeComm);
        if(!activeRank) std::cout<<"(f_min , f_max): ( "<<vec_min<<" , "<<vec_max<<" ) "<<std::endl;

        vec_min=vecMin(Mrhs+nodeLocalBegin,local_dof,activeComm);
        vec_max=vecMax(Mrhs+nodeLocalBegin,local_dof,activeComm);
        if(!activeRank) std::cout<<"(Mf_min , Mf_max): ( "<<vec_min<<" , "<<vec_max<<" ) "<<std::endl;

        linalg::CG(mesh,u,Mrhs,maxIter,cg_tol);

        fem::operators::poisson::deallocateWorkSpace();

        char fPrefix[256];
        sprintf(fPrefix,"%s","possoin");
        const char * varNames[]={"U","rhs","Mrhs"};
        const double * var[]={u,rhs,Mrhs};
        io::vtk::mesh2vtuFine(mesh,fPrefix,0,NULL,NULL,3,varNames,var);



    }






    delete [] rhs;
    delete [] u;

    //delete [] mesh;

    MPI_Finalize();


}
