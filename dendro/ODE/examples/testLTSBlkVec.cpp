/**
 * @file testLTSBlkVec.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief File to test LTS block vector routines, unzip form DG vector zipDG and zipCG when all the blocks are synced. 
 * @version 0.1
 * @date 2021-03-19
 * @copyright Copyright (c) 2021
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
#include "enutsUtils.tcc"




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

    std::function<void(double,double,double,double*)> func_1 =[](const double xx,const double yy,const double zz,double * var)
    {
        var[0] = sin(2*M_PI*xx)*sin(2*M_PI*yy)*sin(2*M_PI*zz);
        var[1] = 1 * sin(2*M_PI*xx)*sin(2*M_PI*yy)*sin(2*M_PI*zz); //1/(1+ sin(2*M_PI*xx)*sin(2*M_PI*yy)*sin(2*M_PI*zz));
        return;
    };

    std::vector<ot::TreeNode> tmpNodes;
    function2Octree(func,tmpNodes,m_uiMaxDepth,wavelet_tol,eOrder,comm);

    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),eOrder,comm,1,ot::SM_TYPE::FDM);
    mesh->setDomainBounds(pt_min,pt_max);

    const unsigned int DOF = 2;
    DendroScalar* funcValCG = mesh->createCGVector<DendroScalar>(func_1,DOF);
    DendroScalar* funcValDG = mesh->createDGVector<DendroScalar>(func_1,DOF);

    DendroScalar* funcValCGUnzip = mesh->createUnZippedVector<DendroScalar>(DOF);
    DendroScalar* funcValDGUnzip = mesh->createUnZippedVector<DendroScalar>(DOF);

    DendroScalar* fVec1_CG = mesh->createCGVector<DendroScalar>(0,DOF);
    DendroScalar* fVec1_DG = mesh->createDGVector<DendroScalar>(0,DOF);

    mesh->readFromGhostBeginEleDGVec(funcValDG,DOF);
    mesh->readFromGhostEndEleDGVec(funcValDG  ,DOF);

    mesh->readFromGhostBegin(funcValCG,DOF);
    mesh->readFromGhostEnd(funcValCG,DOF);

    // unzip from the CG vector. 
    mesh->unzip(funcValCG,funcValCGUnzip,DOF);
    
    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        UnZip Test Check Begin  (unzip from a DG vec)    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    { 
        // unzip by DG vector test
        const unsigned int dof=DOF;
        
        const unsigned int uzSz = mesh->getDegOfFreedomUnZip();
        const unsigned int nPe = mesh->getNumNodesPerElement();
        const unsigned int dgSz = mesh->getAllElements().size()*nPe;
        
        mesh->unzipDG(funcValDG,funcValDGUnzip,dof);
        
        ot::test::isUnzipValid(mesh,            funcValDGUnzip, func,1e-3);
        ot::test::isUnzipValid(mesh, funcValDGUnzip + 1* uzSz , func,1e-3);
        

    }

    

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        UnZip Test Check End  (unzip from a DG vec)    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }


    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Blk async vec  (begin)    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    {

        if(mesh->isActive())
        {
            unsigned int active_rank = mesh->getMPIRank();
            MPI_Comm aComm  = mesh->getMPICommunicator();

            std::vector<ts::BlockTimeStep<DendroScalar>> m_uiBVec;

            const std::vector<ot::Block>& blkList = mesh->getLocalBlockList();
            m_uiBVec.resize(blkList.size());

            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                const unsigned int sz[3] = { blkList[blk].getAllocationSzX(), blkList[blk].getAllocationSzY(), blkList[blk].getAllocationSzZ()};
                m_uiBVec[blk].alloc_vec(2, blk, sz, DOF);
            }

            const unsigned int INDEX_PT_SYNC=1;
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                m_uiBVec[blk]._time = 0;
                m_uiBVec[blk]._rks  = 0;
                m_uiBVec[blk]._vec[0].copyFromUnzip(mesh,funcValDGUnzip, true, DOF);
                m_uiBVec[blk]._vec[0].mark_synced();

                // second vector copy without the padding region. 
                m_uiBVec[blk]._vec[INDEX_PT_SYNC].copyFromUnzip(mesh,funcValDGUnzip, false, DOF);
                m_uiBVec[blk]._vec[INDEX_PT_SYNC].mark_unsynced();
            
            }

            // sync m_uiBVec[blk]._vec[1] using partial block sync
            for(unsigned int blk =0; blk < blkList.size(); blk++)
                ts::sync_blk_padding(mesh,funcValDG,m_uiBVec,blk,INDEX_PT_SYNC,DOF);
            
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                const unsigned int bflag = blkList[blk].getBlkNodeFlag();

                if(!bflag)
                {
                    const unsigned int lx     =  blkList[blk].getAllocationSzX();
                    const unsigned int ly     =  blkList[blk].getAllocationSzY();
                    const unsigned int lz     =  blkList[blk].getAllocationSzZ();
                    const unsigned int offset =  blkList[blk].getOffset();

                    const DendroScalar* bdata = m_uiBVec[blk]._vec[INDEX_PT_SYNC].data() ;
                    for(unsigned int node = 0; node < (lx*ly*lz); node++)
                    {
                        for (unsigned int v=0; v < DOF; v++)
                        {   
                            const DendroScalar v1 = bdata[(v*lx*ly*lz) + node];
                            const DendroScalar v2 = funcValDGUnzip[ v * mesh->getDegOfFreedomUnZip() + offset +node];
                            const DendroScalar diff_abs = fabs(v1-v2);
                            if(diff_abs > 1e-6)
                            {
                                std::cout<<RED<<"Failed: [Error in parital block sync]: rank: "<<active_rank<<" blk: "<<blk<<" v1: "<<v1<<" v2: "<<v2<<NRM<<std::endl;
                            }
                        }
                    }
                }

            }

            // zip it back to cg vec and dg vec
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                m_uiBVec[blk]._vec[0].zip(mesh,fVec1_CG,DOF);
                m_uiBVec[blk]._vec[0].zipDG(mesh,fVec1_DG,DOF);
            }
                
            mesh->readFromGhostBegin(fVec1_CG,DOF);
            mesh->readFromGhostEnd(fVec1_CG,DOF);

            mesh->readFromGhostBeginEleDGVec(fVec1_DG,DOF);
            mesh->readFromGhostEndEleDGVec(fVec1_DG,DOF);

            const unsigned int dg_sz = mesh->getDegOfFreedomDG();
            const unsigned int cg_sz = mesh->getDegOfFreedom();
            const unsigned int nPe   = mesh->getNumNodesPerElement();
            DendroScalar linf[DOF];

            for(unsigned int v=0; v < DOF; v++)
                linf[v] = normLInfty(fVec1_DG + v*dg_sz + mesh->getElementLocalBegin()*nPe , funcValDG + v*dg_sz + mesh->getElementLocalBegin()*nPe, mesh->getNumLocalMeshElements()*nPe , aComm);

            if(!active_rank)
            {
                for(unsigned int v=0; v < DOF; v++)
                    std::cout<<" var dg : "<<v<<" linf : "<<linf[v]<<std::endl;
            }

            for(unsigned int v=0; v < DOF; v++)
                linf[v] = normLInfty(fVec1_CG + v*cg_sz + mesh->getNodeLocalBegin(), funcValCG + v*cg_sz + mesh->getNodeLocalBegin(), mesh->getNumLocalMeshNodes(), aComm);

            if(!active_rank)
            {
                for(unsigned int v=0; v < DOF; v++)
                    std::cout<<" var cg : "<<v<<" linf : "<<linf[v]<<std::endl;
            }


            
        }

        

    }




    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Blk async vec  (end)    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        
    }





    mesh->destroyVector(funcValDG);
    mesh->destroyVector(funcValCG);
    mesh->destroyVector(funcValCGUnzip);
    mesh->destroyVector(funcValDGUnzip);
    mesh->destroyVector(fVec1_CG);
    mesh->destroyVector(fVec1_DG);


    delete mesh;
    MPI_Finalize();
    return 0;
    
}