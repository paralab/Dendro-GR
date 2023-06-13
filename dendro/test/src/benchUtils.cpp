/**
 * @file benchUtils.cpp
 * @brief run/test simple mesh utility funtions
 * @version 0.1
 * @date 2022-01-21
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include "dendro.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include "meshUtils.h"
#include "oct2vtk.h"
#include "profiler.h"
#include "meshTestUtils.h"



int main(int argc, char** argv)
{   
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(argc<4)
    {
        if(!rank) std::cout<<"Usage: "<<argv[0]<<" maxDepth wavelet_tol partition_tol eleOrder"<<std::endl;
        MPI_Abort(comm,0);
    }

    m_uiMaxDepth=atoi(argv[1]);
    double wavelet_tol=atof(argv[2]);
    double partition_tol=atof(argv[3]);
    unsigned int eOrder=atoi(argv[4]);
    const unsigned int iter=16;

    if(!rank)
    {
        std::cout<<YLW<<"maxDepth: "<<m_uiMaxDepth<<NRM<<std::endl;
        std::cout<<YLW<<"wavelet_tol: "<<wavelet_tol<<NRM<<std::endl;
        std::cout<<YLW<<"partition_tol: "<<partition_tol<<NRM<<std::endl;
        std::cout<<YLW<<"eleOrder: "<<eOrder<<NRM<<std::endl;

    }

    _InitializeHcurve(m_uiDim);

    // function that we need to interpolate.
    const double d_min=-10;
    const double d_max=10;
    double dMin[]={d_min,d_min,d_min};
    double dMax[]={d_max,d_max,d_max};

    Point pt_min(d_min,d_min,d_min);
    Point pt_max(d_max,d_max,d_max);
    //@note that based on how the functions are defined (f(x), dxf(x), etc) the compuatational domain is equivalent to the grid domain.
    std::function<void(double,double,double,double*)> func =[d_min,d_max](double x,double y,double z, double * var){

        double ca[]={-2,0,0};
        double cb[]={2,0,0};
        double rra = ((x-ca[0]) * (x-ca[0]) + (y-ca[1])* (y-ca[1]) + (z-ca[2])* (z-ca[2]));
        double rrb = ((x-cb[0]) * (x-cb[0]) + (y-cb[1])* (y-cb[1]) + (z-cb[2])* (z-cb[2]));
        var[0] = exp(-rra) + exp(-rrb); 
        //(sin(2*M_PI*xx)*sin(2*M_PI*yy)*sin(2*M_PI*zz));
        return;
         
    };
    
    std::function<double(double,double,double)> dx_func=[d_min,d_max](const double x,const double y,const double z){ 

        double xx = (x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double yy = (y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double zz = (z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;

        // if((xx < -0.5 || xx > 0.5) || ( yy < -0.5 || yy > 0.5) || (zz < -0.5 || zz > 0.5) )
        //     return 0.0;
    
        return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(cos(2*M_PI*xx)*sin(2*M_PI*yy)*sin(2*M_PI*zz));

    };
    
    std::function<double(double,double,double)> dy_func=[d_min,d_max](const double x,const double y,const double z){ 

        double xx = (x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double yy = (y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double zz = (z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;

        // if((xx < -0.5 || xx > 0.5) || ( yy < -0.5 || yy > 0.5) || (zz < -0.5 || zz > 0.5) )
        //     return 0.0;
        
        return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*xx)*cos(2*M_PI*yy)*sin(2*M_PI*zz));

    
    
    };

    std::function<double(double,double,double)> dz_func=[d_min,d_max](const double x,const double y,const double z){ 

        double xx = (x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double yy = (y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double zz = (z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;

        // if((xx < -0.5 || xx > 0.5) || ( yy < -0.5 || yy > 0.5) || (zz < -0.5 || zz > 0.5) )
        //     return 0.0;
        
        return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*xx)*sin(2*M_PI*yy)*cos(2*M_PI*zz));

    };


    std::vector<ot::TreeNode> tmpNodes;
    std::function<double(double,double,double)> fr = [func,d_min,d_max](double x,double y,double z){
        double xx = (x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double yy = (y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double zz = (z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double var; 
        func(xx,yy,zz,&var); 
        return var;
    };
    function2Octree(fr,tmpNodes,m_uiMaxDepth,wavelet_tol,eOrder,comm);

    ot::Mesh * mesh= ot::createMesh(tmpNodes.data(),tmpNodes.size(),eOrder,comm,1,ot::SM_TYPE::FDM,DENDRO_DEFAULT_GRAIN_SZ,partition_tol,DENDRO_DEFAULT_SF_K);
    mesh->setDomainBounds(pt_min, pt_max);
    //io::vtk::mesh2vtuFine(mesh,"begin",0,NULL,NULL,0,NULL,NULL,0,NULL,NULL,false);
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);
    if(!rank)
        printf("created mesh with lev_min = %d and lev_max=%d\n",lmin,lmax);
    
    double * u = mesh->createCGVector<double>(func,1);
    double * v_dg    = mesh->createDGVector(func,1);
    
    double * u_unzip = mesh->createUnZippedVector<double>(1);
    double * v_unzip = mesh->createUnZippedVector<double>(1);
    
    profiler_t t_unzip_cg;
    profiler_t t_unzip_dg;

    mesh->readFromGhostBegin(u,1);
    mesh->readFromGhostEnd(u,1);

    mesh->readFromGhostBeginEleDGVec(v_dg,1);
    mesh->readFromGhostEndEleDGVec(v_dg,1);

    mesh->unzip_scatter(u,u_unzip,1);
    bool valid = ot::test::isUnzipValid(mesh,u_unzip,fr,1e-3);
    if(!rank)
        std::cout<<"unzip valid scatter: "<<valid<<std::endl;

    mesh->unzipDG_scatter(v_dg,v_unzip,1);
    valid = ot::test::isUnzipValid(mesh,v_unzip,fr,1e-3);
    if(!rank)
        std::cout<<"unzip valid dg scatter: "<<valid<<std::endl;

    t_unzip_cg.start();
    for(unsigned int i=0; i < iter;i++)
        mesh->unzip_scatter(u,u_unzip,1);
    t_unzip_cg.stop();

    t_unzip_dg.start();
    for(unsigned int i=0; i < iter;i++)
        mesh->unzipDG_scatter(v_dg,v_unzip,1);
    t_unzip_dg.stop();

    

    double t_local;
    double t_stat[3];
    
    t_local=t_unzip_cg.seconds/(double)iter;
    par::computeOverallStats(&t_local,t_stat,comm,"cg unzip_scatter");
    
    t_local=t_unzip_dg.seconds/(double)iter;
    par::computeOverallStats(&t_local,t_stat,comm,"dg unzip");
    
    t_unzip_cg.clear();
    t_unzip_dg.clear();

    mesh->unzip(u,u_unzip,1);
    valid = ot::test::isUnzipValid(mesh,u_unzip,fr,1e-3);
    if(!rank)
        std::cout<<"unzip valid gather: "<<valid<<std::endl;


    mesh->unzipDG(v_dg,v_unzip,1);
    valid = ot::test::isUnzipValid(mesh,v_unzip,fr,1e-3);
    if(!rank)
        std::cout<<"unzip valid dg gather: "<<valid<<std::endl;

    
    t_unzip_cg.start();
    for(unsigned int i=0; i < iter;i++)
        mesh->unzip(u,u_unzip,1);
    t_unzip_cg.stop();

    t_unzip_dg.start();
    for(unsigned int i=0; i < iter;i++)
        mesh->unzipDG(v_dg,v_unzip,1);
    t_unzip_dg.stop();
    
    t_local=t_unzip_cg.seconds/(double)iter;;
    par::computeOverallStats(&t_local,t_stat,comm,"cg unzip_gather");
    
    t_local=t_unzip_dg.seconds/(double)iter;
    par::computeOverallStats(&t_local,t_stat,comm,"dg unzip_gather");
    

    mesh->destroyVector(u);
    mesh->destroyVector(v_dg);

    mesh->destroyVector(u_unzip);
    mesh->destroyVector(v_unzip);

    MPI_Finalize();

    return 0;
}