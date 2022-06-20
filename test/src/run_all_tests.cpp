//
// Created by milinda on 3/31/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains collection of test cases, to check the all mesh, interp, ghost sync functionality. 
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


    //std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){ return (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
    /*std::function<double(double,double,double)> func =[](const double x,const double  y,const double z){ return (x*x+y*y+z*z);};
    std::function<double(double,double,double)> dx_func =[](const double x,const double  y,const double z){ return (2*x);};
    std::function<double(double,double,double)> dy_func =[](const double x,const double  y,const double z){ return (2*y);};
    std::function<double(double,double,double)> dz_func =[](const double x,const double  y,const double z){ return (2*z);};*/
    //std::function<double(double,double,double)> func =[](const double x,const double  y,const double z){ return (x+y+z);};

    std::vector<ot::TreeNode> tmpNodes;
    function2Octree(func,tmpNodes,m_uiMaxDepth,wavelet_tol,eOrder,comm);

    DendroIntL localSz,globalSz;
    localSz=tmpNodes.size();
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" # f20. octants: "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=100;//DENDRO_DEFAULT_GRAIN_SZ;

    bool isActive;
    MPI_Comm commActive;
    const int p_npes_prev=binOp::getPrevHighestPowerOfTwo(std::max((globalSz/grainSz),1LL));
    const int p_npes_next=binOp::getNextHighestPowerOfTwo(std::max((globalSz/grainSz),1LL));

    int p_npes=std::max((globalSz/grainSz),1LL);
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

    shrinkOrExpandOctree(tmpNodes,DENDRO_DEFAULT_LB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

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


        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,DENDRO_DEFAULT_LB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,DENDRO_DEFAULT_SF_K,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,DENDRO_DEFAULT_LB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();


        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;




        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,DENDRO_DEFAULT_LB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
        tmpNodes.clear();

        localSz=balOct.size();


    }
    MPI_Comm_free(&commActive);

    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;

    ot::Mesh* mesh=new ot::Mesh(balOct,1,eOrder,comm);
    mesh->setDomainBounds(pt_min,pt_max);

    DendroIntL numVertices=mesh->getNumLocalMeshNodes();
    DendroIntL numVertices_g;

    par::Mpi_Reduce(&numVertices,&numVertices_g,1,MPI_SUM,0,comm);

    if(!rank)
        std::cout<<"Total number of vertices: "<<numVertices_g<<std::endl;


    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) std::cout<<"mesh generation ended "<<std::endl;
    //io::vtk::mesh2vtu(&mesh, "balOct",0,NULL,NULL,0,NULL,NULL);

    std::vector<double> funcVal;
    std::vector<double> funcValUnZip;
    std::vector<double> dx_funcVal;
    std::vector<double> dx_funcVal1;

    std::vector<double> dy_funcVal;
    std::vector<double> dy_funcVal1;

    std::vector<double> dz_funcVal;
    std::vector<double> dz_funcVal1;

    


    // refine some elements to make it adaptive, 
    // this is explicit refinement , for the given function above it sometimes generate regular grid by wavelets, this below code is to 
    // enforce artificial refinement, so that, the unzip and other routines can be tested on adaptive grid. 
    bool isRefine1 = false;
    unsigned int refineCount=0;
    do{

        mesh->createVector(funcVal,func);
        mesh->createUnZippedVector(funcValUnZip);
        mesh->performGhostExchange(funcVal);
        mesh->unzip(funcVal.data(),funcValUnZip.data(),1);
        
        std::vector<unsigned int> refine_flag;
        const double * refvar [] = {funcValUnZip.data()};
        const unsigned int varids [] ={0};
        std::function<double(double,double,double)> waveletTolFunc =[wavelet_tol](double x,double y, double z) {
            return wavelet_tol;
        };
        bool includebdy=false;
        isRefine1 = wavelet::compute_wavelet_remesh_flags(mesh,refine_flag,refvar,varids,1,waveletTolFunc,0.01,includebdy);
        //mesh->isReMeshUnzip((const double **)vars ,(const unsigned int *)vIDs,1,w_tol_coarsen);
        // if(!rank)
        //     std::cout<<"is refine triggered for coarsening tolerance : "<<isRefine<<std::endl;

        // const ot::TreeNode * pNodes = mesh->getAllElements().data();
        // refine_flag.reserve(mesh->getNumLocalMeshElements());
        // for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele++)
        // {
        //     // if( ((ele + NUM_CHILDREN-1) < mesh->getElementLocalEnd()) && (pNodes[ele].getParent() == pNodes[ele + NUM_CHILDREN-1].getParent()) )
        //     // {
        //     //     for(unsigned int c=0; c < NUM_CHILDREN; c++)
        //     //         refine_flag.push_back(OCT_COARSE);
                
        //     //     ele += (NUM_CHILDREN-1);
        //     // }else
        //     if( (ele % 3) == 0 )
        //         refine_flag.push_back(OCT_SPLIT);
        //     else
        //         refine_flag.push_back(OCT_NO_CHANGE);
            
        // }

        
        if(isRefine1)
        {
            mesh->setMeshRefinementFlags(refine_flag);
            ot::Mesh* newMesh = mesh->ReMesh();
            DendroIntL localSz = mesh->getNumLocalMeshElements();
            DendroIntL gSz_new, gSz_old;

            par::Mpi_Reduce(&localSz,&gSz_old,1,MPI_SUM,0,comm);
            localSz = newMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&gSz_new,1,MPI_SUM,0,comm);
            if(!rank)
                std::cout<<"old mesh size: "<<gSz_old<<" new mesh size: "<<gSz_new<<std::endl;

            std::swap(newMesh,mesh);
            delete newMesh;

            refineCount++;

        }

        funcVal.clear();
        funcValUnZip.clear();
        

    }while( refineCount < 2 && isRefine1);


    

    mesh->createVector(funcVal,func);
    mesh->createUnZippedVector(funcValUnZip);
    mesh->performGhostExchange(funcVal);

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal Value Check Begin     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    //assert(ot::test::isElementalNodalValuesValid(&mesh,&(*(funcVal.begin())),func,1e-3));
    ot::test::isElementalNodalValuesValid(mesh,&(*(funcVal.begin())),func,1e-3);
    

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal Value Check end     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }




    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal to Hanging Computation Check Begin     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    //assert(ot::test::isElementalNodalValuesValid(&mesh,&(*(funcVal.begin())),func,1e-3));
    ot::test::isElementalContributionValid(mesh,func,func,1e-3);

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal to Hanging Computation Check End     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }



    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        UnZip Test Check Begin     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    
    {
        mesh->unzip(&(*(funcVal.begin())),&(*(funcValUnZip.begin())));
        ot::test::isUnzipValid(mesh,&(*(funcValUnZip.begin())),func,1e-3);
    }
    

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        UnZip Test Check End     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }
    

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        UnZip Test Check Begin  (unzip from a DG vec)    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }


    { 
        // unzip by DG vector test
        const unsigned int dof=2;
        double * funcValDG = mesh->createDGVector(0.0,dof);
        double * funcValDGUnzip = mesh->createUnZippedVector<double>(dof);
        const unsigned int nPe = mesh->getNumNodesPerElement();
        const unsigned int uzSz = mesh->getDegOfFreedomUnZip();
        const unsigned int dgSz = mesh->getAllElements().size()*nPe;
        std::vector<double> nvec;
        nvec.resize(mesh->getNumNodesPerElement());
        for(unsigned int ele=mesh->getElementLocalBegin(); ele< mesh->getElementLocalEnd(); ele++)
        {
            mesh->getElementNodalValues(funcVal.data(),nvec.data(),ele,false);
            
            for(unsigned int v=0; v < dof; v++)
            {
                for(unsigned int node=0; node < nPe ; node ++)
                    funcValDG[v*dgSz + ele*nPe + node] = nvec[node];
            }
            
        }

        mesh->readFromGhostBeginEleDGVec(funcValDG,dof);
        mesh->readFromGhostEndEleDGVec(funcValDG,dof);

        mesh->unzipDG(funcValDG,funcValDGUnzip,dof);
        // //std::cout<<"unzip ended "<<std::endl;
        // //assert(ot::test::isUnzipValid(&mesh,&(*(funcValUnZip.begin())),func,1e-3));
        ot::test::isUnzipValid(mesh,            funcValDGUnzip, func,1e-3);
        ot::test::isUnzipValid(mesh, funcValDGUnzip + 1* uzSz , func,1e-3);
        mesh->destroyVector(funcValDG);
        mesh->destroyVector(funcValDGUnzip);

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
        std::cout<<"        Wavelet Check Begin     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }



    {
        wavelet::WaveletEl wEle((RefElement*)mesh->getReferenceElement());
        std::vector<double> blkIn;
        std::vector<double> wCout;

        const std::vector<ot::Block>& blkList = mesh->getLocalBlockList();
        for(unsigned int blk=0; blk <blkList.size(); blk++)
        {   

            const unsigned int pw = blkList[blk].get1DPadWidth();
            
            if((eOrder>>1u) != pw)
            {
                std::cout<<" padding width should be half the eleOrder for generic wavelet computations. "<<std::endl;
                MPI_Abort(comm,0);
            }
            
            const unsigned int nx = (2*eOrder+1);
            const unsigned int ny = (2*eOrder+1);
            const unsigned int nz = (2*eOrder+1);

            //std::cout<<"nx "<<nx<<std::endl;
            
            blkIn.resize(nx*ny*nz);
            const unsigned int isz[] = {nx,ny,nz};

            const unsigned int bflag = blkList[blk].getBlkNodeFlag();
            
            if(bflag!=0) 
                continue;

            for(unsigned int ele =blkList[blk].getLocalElementBegin(); ele < blkList[blk].getLocalElementEnd(); ele++)
            {
                mesh->getUnzipElementalNodalValues(funcValUnZip.data(),blk,ele,blkIn.data(),true);
                wEle.compute_wavelets_3D(blkIn.data(),isz,wCout);
            
                // double l_max = normL2(wCout.data(),wCout.size());

                // if(l_max > 1e-2)
                // {
                //     std::cout<<"rank: "<<mesh->getMPIRank()<<" blk : "<<blk<<" ele ID: "<<ele<<" wavelets(max) : "<<l_max<<std::endl;
                //     for(unsigned int k=0; k < wCout.size(); k++)
                //      std::cout<<"wCout["<<k<<"]: "<<wCout[k]<<std::endl;
                // }
                
            }

        }
    
    }
    
    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Wavelet Check End     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        interpolation to sphere check begin    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }
    
    //ot::test::isInterpToSphereValid(mesh,&(*(funcVal.begin())),func,0.1,dMin,dMax,1e-3);
    //ot::test::isSphereInterpValid(mesh,&(*(funcVal.begin())),func,0.1,1e-2,pt_min,pt_max);

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        interpolation to sphere check end     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    // test elemental ghost exchange
    
    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        elemental ghost sync test begin    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    {
        unsigned int * cellVec = mesh->createElementVector(0u);
        const ot::TreeNode* const pNodes = mesh->getAllElements().data();

        for (unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele++)
            cellVec[ele] = pNodes[ele].getLevel();

        mesh->readFromGhostBeginElementVec(cellVec,1);
        mesh->readFromGhostEndElementVec(cellVec,1);

        for (unsigned int ele = mesh->getElementPreGhostBegin(); ele < mesh->getElementPreGhostEnd(); ele++)
        {
            if( !(cellVec[ele] ==0 || cellVec[ele] == pNodes[ele].getLevel()))
                std::cout<<"rank: "<<mesh->getMPIRank()<<" elemental ghost sync error in pre ghost "<<std::endl;
        }

        for (unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele++)
        {
            if( !(cellVec[ele] == pNodes[ele].getLevel()))
                std::cout<<"rank: "<<mesh->getMPIRank()<<" elemental ghost sync currupts local values "<<std::endl;
        }

        for (unsigned int ele = mesh->getElementPostGhostBegin(); ele < mesh->getElementPostGhostEnd(); ele++)
        {
            if( !(cellVec[ele] ==0 || cellVec[ele] == pNodes[ele].getLevel()))
                std::cout<<"rank: "<<mesh->getMPIRank()<<" elemental ghost sync error in post ghost cell val "<<cellVec[ele]<<" nodel lev: "<<pNodes[ele].getLevel()<<std::endl;
        }

        mesh->destroyVector(cellVec);
    }


    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        elemental ghost sync test end    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    
    
    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        DG ghost sync test begin    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }


    {
        const ot::TreeNode* const pNodes = mesh->getAllElements().data();
        double * dg_vec = mesh->createDGVector<double>(0,2);
        unsigned int dg_sz = mesh->getDegOfFreedomDG();

        const unsigned int nPe = mesh->getNumNodesPerElement();
        for(unsigned int v=0; v < 2; v++)
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

        mesh->readFromGhostBeginEleDGVec(dg_vec,2);
        mesh->readFromGhostEndEleDGVec(dg_vec,2);

        
        if(npes >1 && mesh->isActive())
        {
            const std::vector<unsigned int> recvEleOffset = mesh->getElementRecvOffsets();
            const std::vector<unsigned int> recvEleCounts = mesh->getElementRecvCounts();
            const std::vector<unsigned int> recvESM       = mesh->getRecvElementSM(); 

            for(unsigned int v=0; v < 2; v++)
            {

                for(unsigned int p=0; p < mesh->getMPICommSize(); p++)
                {
                    for(unsigned int ee=recvEleOffset[p]; ee < (recvEleOffset[p]+ recvEleCounts[p]); ee++)
                    {   
                        
                        const unsigned int ele = recvESM[ee];
                        for(unsigned int k=0; k < eOrder+1 ; k++)
                        for(unsigned int j=0; j < eOrder+1; j++)
                        for(unsigned int i=0; i < eOrder+1; i++)
                        {
                            double x = pNodes[ele].minX() + i*((1u<<m_uiMaxDepth - pNodes[ele].getLevel())/eOrder);
                            double y = pNodes[ele].minY() + j*((1u<<m_uiMaxDepth - pNodes[ele].getLevel())/eOrder);
                            double z = pNodes[ele].minZ() + k*((1u<<m_uiMaxDepth - pNodes[ele].getLevel())/eOrder);

                            double diff = std::fabs((dg_vec + v*dg_sz)[ele*nPe + k*(eOrder+1)*(eOrder+1) + j*(eOrder+1) + i ] - func(x,y,z));
                            
                            if(diff > 1e-3)
                                std::cout<<"rank: "<<rank<<" ele: "<<ele<<" DG vec:"<<(dg_vec + v*dg_sz)[ele*nPe + k*(eOrder+1)*(eOrder+1) + j*(eOrder+1) + i ]<<" func: "<<func(x,y,z)<<" diff : "<<diff<<std::endl;
                        
                        }
                    }
                }
            }

        }
        

        mesh->destroyVector(dg_vec);

    }


    

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        DG ghost sync test end    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }




    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        FD stencil test begin    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    {

        mesh->createVector(dx_funcVal1,dx_func);
        mesh->createVector(dy_funcVal1,dy_func);
        mesh->createVector(dz_funcVal1,dz_func);

        mesh->createVector(dx_funcVal);
        mesh->createVector(dy_funcVal);
        mesh->createVector(dz_funcVal);


        Stencil<double,5,2> D1_Order4StencilCentered_x(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_X);
        Stencil<double,5,4> D1_Order4StencilBackward_x(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_X);
        Stencil<double,5,0> D1_Order4StencilForward_x(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_X);

        Stencil<double,5,1> D1_Order4StencilUpWind_x(fd::D1_ORDER_4_UPWIND,5,StencilDirection::STENCIL_DIR_X);
        Stencil<double,5,3> D1_Order4StencilDownWind_x(fd::D1_ORDER_4_DOWNWIND,5,StencilDirection::STENCIL_DIR_X);


        Stencil<double,5,2> D1_Order4StencilCentered_y(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_Y);
        Stencil<double,5,4> D1_Order4StencilBackward_y(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_Y);
        Stencil<double,5,0> D1_Order4StencilForward_y(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_Y);

        Stencil<double,5,1> D1_Order4StencilUpWind_y(fd::D1_ORDER_4_UPWIND,5,StencilDirection::STENCIL_DIR_Y);
        Stencil<double,5,3> D1_Order4StencilDownWind_y(fd::D1_ORDER_4_DOWNWIND,5,StencilDirection::STENCIL_DIR_Y);

        Stencil<double,5,2> D1_Order4StencilCentered_z(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_Z);
        Stencil<double,5,4> D1_Order4StencilBackward_z(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_Z);
        Stencil<double,5,0> D1_Order4StencilForward_z(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_Z);

        Stencil<double,5,1> D1_Order4StencilUpWind_z(fd::D1_ORDER_4_UPWIND,5,StencilDirection::STENCIL_DIR_Z);
        Stencil<double,5,3> D1_Order4StencilDownWind_z(fd::D1_ORDER_4_DOWNWIND,5,StencilDirection::STENCIL_DIR_Z);



        mesh->performGhostExchange(funcVal);

        #if 0
            double* vars[3];
            vars[0]=&(*(funcVal.begin()));
            vars[1]=&(*(dy_funcVal1.begin()));
            vars[2]=&(*(dz_funcVal1.begin()));

            std::cout<<"write begin: "<<std::endl;
            io::varToRawData((const ot::Mesh*)mesh,(const double **)vars,3,NULL,"dendro_raw");
            std::cout<<"write end: "<<std::endl;
        #endif
        
        /*std::vector<double> unZipIn;
        std::vector<double> unZipOut;
        //@milinda : If you are using grad functions defined in rkTransport utills, make sure to use the correct domain paramerters when computing the h.
        mesh.createUnZippedVector(unZipIn,0.0);
        mesh.createUnZippedVector(unZipOut,0.0);

        mesh.unzip(&(*(funcVal.begin())),&(*(unZipIn.begin())));

        grad(&mesh,0,&(*(unZipIn.begin())),&(*(unZipOut.begin())));
        mesh.zip(&(*(unZipOut.begin())),&(*(dx_funcVal.begin())));

        grad(&mesh,1,&(*(unZipIn.begin())),&(*(unZipOut.begin())));
        mesh.zip(&(*(unZipOut.begin())),&(*(dy_funcVal.begin())));

        grad(&mesh,2,&(*(unZipIn.begin())),&(*(unZipOut.begin())));
        mesh.zip(&(*(unZipOut.begin())),&(*(dz_funcVal.begin())));

        //std::cout<<"zip vec: "<<funcVal.size()<<" unZipVec: "<<unZipIn.size()<<std::endl;

        unZipOut.clear();
        unZipIn.clear();*/

        const std::vector<ot::Block>& blockList=mesh->getLocalBlockList();
        std::vector<ot::TreeNode> localBlocks;

        // for(unsigned int e=0;e<blockList.size();e++)
        //     localBlocks.push_back(blockList[e].getBlockNode());
        //treeNodesTovtk(localBlocks,rank,"blocks");

        mesh->applyStencil(funcVal, dx_funcVal, D1_Order4StencilCentered_x, D1_Order4StencilBackward_x, D1_Order4StencilForward_x);
        mesh->applyStencil(funcVal, dy_funcVal, D1_Order4StencilCentered_y, D1_Order4StencilBackward_y, D1_Order4StencilForward_y);
        mesh->applyStencil(funcVal, dz_funcVal, D1_Order4StencilCentered_z, D1_Order4StencilBackward_z, D1_Order4StencilForward_z);

        // Need to do the gost exchange before wirting vtu files.
        mesh->performGhostExchange(funcVal);
        mesh->performGhostExchange(dx_funcVal1);
        mesh->performGhostExchange(dx_funcVal);

        unsigned int numCVars=2;
        const char * cvarNames [] = {"rank","level"};

        double* cVec = mesh->createElementVector(0.0,2);

        if(mesh->isActive())
        {
            for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele ++)
                cVec[ele] = mesh->getMPIRank();

            for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele ++)
                cVec[mesh->getAllElements().size() + ele] = mesh->getAllElements()[ele].getLevel();
        }
        


        unsigned int numPVars=3;
        const char * pVarNames [] ={"F(x)","dxF(x)","SxF(x)"};
        double * pVarVal []={&(*(funcVal.begin())),&(*(dx_funcVal1.begin())),&(*(dx_funcVal.begin()))};
        double * cVarVal []={cVec, cVec+mesh->getAllElements().size() };
        
        const char *fVarNames []={"Time","Cycle"};
        double fVarVal []={0.1,1};

        unsigned int s_val[] ={1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1)};
        unsigned int s_normal[]={0,0,1};
        mesh->setDomainBounds(pt_min,pt_max);
        io::vtk::mesh2vtu_slice(mesh, s_val,s_normal, "f_val_slice",2,(const char **)fVarNames,(const double *)fVarVal,numPVars,(const char** )pVarNames,(const double **)pVarVal);
        io::vtk::mesh2vtuFine(mesh, "f_val",2,(const char **)fVarNames,(const double *)fVarVal,numPVars,(const char** )pVarNames,(const double **)pVarVal,numCVars,cvarNames,(const double**)cVarVal,false);


        mesh->destroyVector(cVec);

        if(mesh->isActive())
        {

            unsigned int nLocalBegin=mesh->getNodeLocalBegin();
            unsigned int nLocalEnd=mesh->getNodeLocalEnd();

            double l2_x=normL2(&(*(dx_funcVal1.begin()+nLocalBegin)),&(*(dx_funcVal.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),mesh->getMPICommunicator());
            double l_inf_x= normLInfty(&(*(dx_funcVal1.begin() + nLocalBegin)), &(*(dx_funcVal.begin() + nLocalBegin)),
                                    (nLocalEnd - nLocalBegin), mesh->getMPICommunicator());

            if(!mesh->getMPIRank()) std::cout<<"l2 norm (x derivative): "<<l2_x<<std::endl;
            if(!mesh->getMPIRank()) std::cout<<"l_inf norm (x derivative): "<<l_inf_x<<std::endl;


            double l2_y=normL2(&(*(dy_funcVal1.begin()+nLocalBegin)),&(*(dy_funcVal.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),mesh->getMPICommunicator());
            double l_inf_y= normLInfty(&(*(dy_funcVal1.begin() + nLocalBegin)), &(*(dy_funcVal.begin() + nLocalBegin)),(nLocalEnd - nLocalBegin), mesh->getMPICommunicator());

            if(!mesh->getMPIRank()) std::cout<<"l2 norm (y derivative): "<<l2_y<<std::endl;
            if(!mesh->getMPIRank()) std::cout<<"l_inf norm (y derivative): "<<l_inf_y<<std::endl;


            double l2_z=normL2(&(*(dz_funcVal1.begin()+nLocalBegin)),&(*(dz_funcVal.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),mesh->getMPICommunicator());
            double l_inf_z= normLInfty(&(*(dz_funcVal1.begin() + nLocalBegin)), &(*(dz_funcVal.begin() + nLocalBegin)),(nLocalEnd - nLocalBegin), mesh->getMPICommunicator());

            if(!mesh->getMPIRank()) std::cout<<"l2 norm (z derivative): "<<l2_z<<std::endl;
            if(!mesh->getMPIRank()) std::cout<<"l_inf norm (z derivative): "<<l_inf_z<<std::endl;

        }

    }

    
    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        FD stencil test end    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }


    

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"       Intergrid transfer test begin    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }


    {
        const ot::TreeNode* const pNodes = mesh->getAllElements().data();
        const double* vars[] = {funcValUnZip.data()};
        const unsigned int vIDs[]={0};
        std::function<double(double,double,double)>w_tol_coarsen=[wavelet_tol](double x,double y, double z) {return wavelet_tol*1e8;};
        std::function<double(double,double,double)>w_tol_refine=[wavelet_tol](double x,double y, double z) {return wavelet_tol;};
        
        // Note: change this to false to preserve the grid generated by the wavelets. 
        bool isRefine = true;
        //mesh->isReMeshUnzip((const double **)vars ,(const unsigned int *)vIDs,1,w_tol_coarsen);
        // if(!rank)
        //     std::cout<<"is refine triggered for coarsening tolerance : "<<isRefine<<std::endl;
        
        std::vector<unsigned int> refine_flag;
        refine_flag.reserve(mesh->getNumLocalMeshElements());
        for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele++)
        {
            if( ((ele + NUM_CHILDREN-1) < mesh->getElementLocalEnd()) && (pNodes[ele].getParent() == pNodes[ele + NUM_CHILDREN-1].getParent()) )
            {
                for(unsigned int c=0; c < NUM_CHILDREN; c++)
                    refine_flag.push_back(OCT_COARSE);
                
                ele += (NUM_CHILDREN-1);
            }else if( (ele % 10) == 0 )
                refine_flag.push_back(OCT_SPLIT);
            else
                refine_flag.push_back(OCT_NO_CHANGE);
            
        }


        // dof 2 function value.
        unsigned int dof=2;
        double * fv2 = mesh->createCGVector<double>(0.0,dof);

        for(unsigned int v=0; v < dof; v++)
            std::memcpy(fv2 + v*mesh->getDegOfFreedom(), funcVal.data(), sizeof(double)*mesh->getDegOfFreedom());

        if(isRefine)
        {

            mesh->setMeshRefinementFlags(refine_flag);
            ot::Mesh* newMesh = mesh->ReMesh();

            DendroIntL localSz = mesh->getNumLocalMeshElements();
            DendroIntL gSz_new, gSz_old;

            par::Mpi_Reduce(&localSz,&gSz_old,1,MPI_SUM,0,comm);
            localSz = newMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&gSz_new,1,MPI_SUM,0,comm);

            if(!rank)
                std::cout<<"old mesh size: "<<gSz_old<<" new mesh size: "<<gSz_new<<std::endl;

            
            funcValUnZip.clear();
            newMesh->createUnZippedVector(funcValUnZip);

            mesh->interGridTransfer(funcVal,newMesh,ot::INTERGRID_TRANSFER_MODE::INJECTION);
            
            // igt for dof=2
            mesh->interGridTransfer(fv2,newMesh,ot::INTERGRID_TRANSFER_MODE::INJECTION,dof);
            
            std::swap(newMesh,mesh);
            delete newMesh;

        }

        // copy the 1st dof to the funcVal to test the multiple dof is working. 
        std::memcpy(funcVal.data(), fv2 + 1*(mesh->getDegOfFreedom()) , sizeof(double)*mesh->getDegOfFreedom());


        mesh->performGhostExchange(funcVal);
        mesh->unzip(funcVal.data(),funcValUnZip.data(),1);

        delete [] fv2;

        if(isRefine)
        {


            if(!rank)
            {
                std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
                std::cout<<"        UnZip Test Check (after remesh) Begin     "<<std::endl;
                std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

            }

                ot::test::isUnzipValid(mesh,&(*(funcValUnZip.begin())),func,1e-3);

            if(!rank)
            {
                std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
                std::cout<<"        UnZip Test Check (after remesh) End     "<<std::endl;
                std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

            }


        }

        unsigned int numPVars1=1;
        const char * pVarNames1 [] ={"F(x)"};
        double * pVarVal1 []={&(*(funcVal.begin()))};

        //unsigned int s_val[] ={1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1)};
        //unsigned int s_normal[]={0,0,1};
        mesh->setDomainBounds(pt_min,pt_max);
        //io::vtk::mesh2vtu_slice(mesh, s_val,s_normal, "f_val_slice",2,(const char **)fVarNames,(const double *)fVarVal,numPVars,(const char** )pVarNames,(const double **)pVarVal);
        io::vtk::mesh2vtuFine(mesh, "f_val_rmesh",0,NULL,NULL,numPVars1,(const char** )pVarNames1,(const double **)pVarVal1);

    }

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"       Intergrid transfer test end    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }
    
    delete mesh;
    if(!rank) std::cout<<" all tests ended "<<std::endl;
    MPI_Finalize();

    return 0;

}


