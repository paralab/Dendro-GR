//
// Created by milinda on 5/24/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
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



int main (int argc, char** argv) {


    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(argc<4)
    {
        if(!rank) std::cout<<"Usage: "<<argv[0]<<" maxDepth wavelet_tol partition_tol eleOrder fPrefix=[sol_step] begin=[0] end=[10] th=[0.01] restore_step[=0] restore_rk[=0]"<<std::endl;
        return 0;
    }

    m_uiMaxDepth=atoi(argv[1]);
    double wavelet_tol=atof(argv[2]);
    double partition_tol=atof(argv[3]);
    unsigned int eOrder=atoi(argv[4]);


    double tBegin=0,tEnd=10,th=0.01;

    char * fPrefix=NULL;
    bool restore_rk=false;
    unsigned int restore_step=0;

    if(argc>5) fPrefix=argv[5];
    if(argc>6) tBegin=atof(argv[6]);
    if(argc>7) tEnd=atof(argv[7]);
    if(argc>8) th=atof(argv[8]);
    if(argc>9) restore_step=atoi(argv[9]);
    if(argc>10) restore_rk=(bool)(atoi(argv[10]));


    if(!rank)
    {
        std::cout<<YLW<<"maxDepth: "<<m_uiMaxDepth<<NRM<<std::endl;
        std::cout<<YLW<<"wavelet_tol: "<<wavelet_tol<<NRM<<std::endl;
        std::cout<<YLW<<"partition_tol: "<<partition_tol<<NRM<<std::endl;
        std::cout<<YLW<<"eleOrder: "<<eOrder<<NRM<<std::endl;

    }

    _InitializeHcurve(m_uiDim);

    // function that we need to interpolate.
    const Point grid_min(0,0,0);
    const Point grid_max((1u<<m_uiMaxDepth),(1u<<m_uiMaxDepth),(1u<<m_uiMaxDepth));

    const Point domain_min=adv_param::domain_min;
    const Point domain_max=adv_param::domain_max;

    std::function<double(double,double,double)> g =[domain_min,domain_max,grid_min,grid_max](const double x,const double y,const double z){
        return (sin((((x-grid_min.x())/(grid_max.x()-grid_min.x()))*(domain_max.x()-domain_min.x())+domain_min.x()))*sin((((y-grid_min.y())/(grid_max.y()-grid_min.y()))*(domain_max.y()-domain_min.y())+domain_min.y()))*sin((((z-grid_min.z())/(grid_max.z()-grid_min.z()))*(domain_max.z()-domain_min.z())+domain_min.z())));
    };
    std::function<double(double,double,double,double)> f =[domain_min,domain_max](const double x,const double y,const double z,const double t){ return 0;};


    //std::function<double(double,double,double)> dx_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(cos(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
    //std::function<double(double,double,double)> dy_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
    //std::function<double(double,double,double)> dz_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};

    /*std::function<double(double,double,double)> func =[](const double x,const double  y,const double z){ return (x*x+y*y+z*z);};
    std::function<double(double,double,double)> dx_func =[](const double x,const double  y,const double z){ return (2*x);};
    std::function<double(double,double,double)> dy_func =[](const double x,const double  y,const double z){ return (2*y);};
    std::function<double(double,double,double)> dz_func =[](const double x,const double  y,const double z){ return (2*z);};*/
    //std::function<double(double,double,double)> func =[](const double x,const double  y,const double z){ return (x+y+z);};
    std::vector<ot::TreeNode> tmpNodes;

    function2Octree(g,tmpNodes,m_uiMaxDepth,wavelet_tol,eOrder,comm);
    std::vector<ot::TreeNode> balOctree;
    ot::TreeNode rootNode;
    unsigned int sf_k=2;

    SFC::parSort::SFC_treeSort(tmpNodes,balOctree,balOctree,balOctree,partition_tol,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,2,comm);
    std::swap(tmpNodes,balOctree);
    balOctree.clear();

    SFC::parSort::SFC_treeSort(tmpNodes,balOctree,balOctree,balOctree,partition_tol,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,2,comm);
    std::swap(tmpNodes,balOctree);
    balOctree.clear();

    SFC::parSort::SFC_treeSort(tmpNodes,balOctree,balOctree,balOctree,partition_tol,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_BALANCE_OCTREE,2,comm);
    tmpNodes.clear();

    DendroIntL balOctSize=balOctree.size();
    DendroIntL balOctSize_g;

    par::Mpi_Reduce(&(balOctSize),&balOctSize_g,1,MPI_SUM,0,comm);
    if(!rank)
        std::cout<<"Total number of bal octs: "<<balOctSize_g<<std::endl;

    ot::Mesh* mesh=new ot::Mesh(balOctree,1,eOrder,comm);

    DendroIntL numVertices=mesh->getNumLocalMeshNodes();
    DendroIntL numVertices_g;

    par::Mpi_Reduce(&numVertices,&numVertices_g,1,MPI_SUM,0,comm);

    if(!rank)
        std::cout<<"Total number of vertices: "<<numVertices_g<<std::endl;


    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) std::cout<<"mesh generation ended "<<std::endl;


    ode::solver::RK45Transport rk45Transport(mesh,tBegin,tEnd,th);
    double b[3]={1.0,1.0,1.0};
    rk45Transport.setParameters(b,g,f,fPrefix,wavelet_tol);

    if(restore_rk)
        rk45Transport.restoreRK45Solver(fPrefix,restore_step,comm);

    rk45Transport.rkSolve();

    delete mesh;

    MPI_Finalize();

    return 0;

}
