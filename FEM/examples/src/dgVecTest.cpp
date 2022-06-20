/**
 * @file dgVecTest.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief : Simple code to perform some examples in DG functionalities, in Dendro. 
 * @version 0.1
 * @date 2020-06-28
 * 
 * @copyright Copyright (c) 2020
 * 
 */

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
#include "oct2vtk.h"
#include "rkTransportUtils.h"
#include "meshTestUtils.h"
#include "rawIO.h"
#include "oda.h"
#include "waveletRefEl.h"
#include "waveletAMR.h"



int main (int argc, char** argv) {

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

    if(!rank)
    {
        std::cout<<YLW<<"maxDepth: "<<m_uiMaxDepth<<NRM<<std::endl;
        std::cout<<YLW<<"wavelet_tol: "<<wavelet_tol<<NRM<<std::endl;
        std::cout<<YLW<<"partition_tol: "<<partition_tol<<NRM<<std::endl;
        std::cout<<YLW<<"eleOrder: "<<eOrder<<NRM<<std::endl;

    }

    _InitializeHcurve(m_uiDim);

    // function that we need to interpolate.
    const double d_min=-0.5;
    const double d_max=0.5;
    double dMin[]={d_min,d_min,d_min};
    double dMax[]={d_max,d_max,d_max};

    Point pt_min(d_min,d_min,d_min);
    Point pt_max(d_max,d_max,d_max);
    //@note that based on how the functions are defined (f(x), dxf(x), etc) the compuatational domain is equivalent to the grid domain.
    std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){

        double xx = (x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double yy = (y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double zz = (z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;

        // if((xx < -0.5 || xx > 0.5) || ( yy < -0.5 || yy > 0.5) || (zz < -0.5 || zz > 0.5) )
        //     return 0.0;
        
        return (sin(2*M_PI*xx)*sin(2*M_PI*yy)*sin(2*M_PI*zz));
    };

    std::vector<ot::TreeNode> tmpNodes;
    function2Octree(func,tmpNodes,m_uiMaxDepth,wavelet_tol,eOrder,comm);

    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),eOrder,comm,1,ot::SM_TYPE::FDM);
    const unsigned int dof=2;
    double* dg_vec = mesh->createDGVector(0.0,dof);
    unsigned int dg_sz = mesh->getDegOfFreedomDG();
    const ot::TreeNode* pNodes = mesh->getAllElements().data();
    const unsigned int nPe = mesh->getNumNodesPerElement();

    for(unsigned int v=0; v < dof; v++)
    {
        for(unsigned int ele=mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele++)
        {
            for(unsigned int k=0; k < eOrder+1 ; k++)
             for(unsigned int j=0; j < eOrder+1; j++)
              for(unsigned int i=0; i < eOrder+1; i++)
              {
                  double x = pNodes[ele].minX() + i*((1u<<m_uiMaxDepth - pNodes[ele].getLevel())/eOrder);
                  double y = pNodes[ele].minY() + j*((1u<<m_uiMaxDepth - pNodes[ele].getLevel())/eOrder);
                  double z = pNodes[ele].minZ() + k*((1u<<m_uiMaxDepth - pNodes[ele].getLevel())/eOrder);
                
                  (dg_vec + v*dg_sz)[ele*nPe + k*(eOrder+1)*(eOrder+1) + j*(eOrder+1) + i ] = func(x,y,z);  

              }
           
        }
    }

    

    std::vector<unsigned int> refine_flags;
    refine_flags.resize(mesh->getNumLocalMeshElements(),0);
    
    for(unsigned int ele = mesh->getElementLocalBegin() ; ele < mesh->getElementLocalEnd(); ele++)
    {
        if( (ele-mesh->getElementLocalBegin()) % 5==0)
            refine_flags[ele-mesh->getElementLocalBegin()]=OCT_SPLIT;
    }

    mesh->setMeshRefinementFlags(refine_flags);
    ot::Mesh* newMesh = mesh->ReMesh();

    double* dg_vecNew = newMesh->createDGVector(0.0,dof);

    mesh->readFromGhostBeginEleDGVec(dg_vec,dof);
    mesh->readFromGhostEndEleDGVec(dg_vec,dof);

    mesh->interGridTransfer_DG(dg_vec,dg_vecNew,newMesh,dof);
    
    const char   * vNames [2]  = {"dgV0","dgV1"};
    const double * pData  [2]  = {dg_vecNew , dg_vecNew + newMesh->getDegOfFreedomDG()};

    newMesh->readFromGhostBeginEleDGVec(dg_vecNew,dof);
    newMesh->readFromGhostEndEleDGVec(dg_vecNew,dof);
    
    io::vtk::mesh2vtuFine(newMesh,"dgVec",0,NULL,NULL,dof,vNames,pData,0,NULL,NULL,true);

    newMesh->destroyVector(dg_vecNew);
    mesh->destroyVector(dg_vec);
    delete newMesh;
    delete mesh;

    MPI_Finalize();


}