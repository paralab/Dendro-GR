/**
 * @file meshUtils.cpp
 * @author Milinda Fernando
 * @brief Contains utility function related to easily generate a mesh.
 * @version 0.1
 * @date 2020-01-01
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */

#include "meshUtils.h"

namespace ot 
{

    Mesh* createMesh(const ot::TreeNode* oct, unsigned int num,unsigned int eleOrder, MPI_Comm comm , unsigned int verbose,  ot::SM_TYPE sm_type, unsigned int grain_sz, double ld_tol, unsigned int sf_k, unsigned int (*getWeight)(const ot::TreeNode *))
    {
        

        int rank, npes;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &npes);

        std::vector<ot::TreeNode> tmpNodes;
        
        if(num > 0)
        {
            tmpNodes.resize(num);
            std::memcpy(tmpNodes.data(),oct,sizeof(ot::TreeNode)*num);

        }
        
        DendroIntL localSz=tmpNodes.size();
        DendroIntL globalSz = 0; 
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);


        par::Mpi_Bcast(&globalSz,1,0,comm);
        
        if(!rank) 
            std::cout<<" number of oct: createMesh: "<<globalSz<<std::endl;

        unsigned int grainSz=grain_sz;
        
        if(grain_sz >=globalSz)
            grainSz = globalSz/2;
        
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

        shrinkOrExpandOctree(tmpNodes,ld_tol,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

        if(!isActive)
            if(tmpNodes.size()!=0)
            {
                std::cout<<"Error: "<<__LINE__<<" rank_g: "<<rank<<" isActive: "<<isActive<<" f2O octants: "<<tmpNodes.size()<<std::endl;
                MPI_Abort(comm,0);
            }
        
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
            
            SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,ld_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,sf_k,commActive);
            std::swap(tmpNodes,tmpVec);
            tmpVec.clear();

            SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,ld_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,sf_k,commActive);
            std::swap(tmpNodes,tmpVec);
            tmpVec.clear();

            
            localSz=tmpNodes.size();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

            if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


        
            SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,ld_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,sf_k,commActive);
            tmpNodes.clear();

            localSz=balOct.size();

        }

        MPI_Comm_free(&commActive);

        // all reduce act as barrier to sync all procs.
        par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
        if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;


        ot::Mesh * mesh=new ot::Mesh(balOct,1,eleOrder,comm,true,sm_type,grain_sz,ld_tol,sf_k,getWeight);
        localSz=mesh->getNumLocalMeshNodes();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
        if(!rank) std::cout<<GRN<<" # of CG nodes (vertices) : "<<globalSz<<NRM<<std::endl;

        return mesh;
    }

    Mesh* createWAMRMesh(std::function<void(double,double,double,double*)> func, double wtol, unsigned int numVars,unsigned int eleOrder, MPI_Comm comm, unsigned int verbose,ot::SM_TYPE sm_type , unsigned int * refIds, unsigned int sz, unsigned int grain_sz , double ld_tol , unsigned int sf_k )
    {
        std::vector<TreeNode> tmpNodes;
        std::vector<unsigned int> refVars;
        refVars.resize(numVars);

        double t_stat;
        double t_stat_g[3];
        
        for(unsigned int i=0; i< refVars.size(); i++)
            refVars[i] = i; 
        
        if(refIds == NULL || sz==0)
            function2Octree(func,numVars,refVars.data(),refVars.size(),tmpNodes,m_uiMaxDepth,wtol,eleOrder,comm);
        else
            function2Octree(func,numVars,refIds,sz,tmpNodes,m_uiMaxDepth,wtol,eleOrder,comm);
        
        return createMesh(tmpNodes.data(),tmpNodes.size(),eleOrder,comm,verbose,sm_type,grain_sz,ld_tol,sf_k);    

        
    }

    void meshWAMRConvergence(ot::Mesh*& pMesh, std::function<void(double,double,double,double*)> func,  double wtol, unsigned int numVars, unsigned int eleOrder,unsigned int * refIds, unsigned int sz, unsigned int maxiter)
    {
        typedef ot::DVector<DendroScalar, unsigned int> DVec;
        MPI_Comm gcomm = pMesh->getMPIGlobalCommunicator();
        int rank, npes;

        MPI_Comm_size(gcomm,&npes);
        MPI_Comm_rank(gcomm,&rank);

        std::vector<unsigned int> refVars;
        unsigned int * rid=NULL;
        unsigned int   nr =0;

        if(refIds==NULL || sz==0)
        {
            refVars.resize(numVars);
        
            for(unsigned int i=0; i< refVars.size(); i++)
                refVars[i] = i;

            rid = refVars.data();
            nr = refVars.size();

        }else
        {
            rid = refIds;
            nr  = sz;
        }

        bool isRefine1=false;
        unsigned int iter=0;
        
        std::function<double(double,double,double)> waveletTolFunc =[wtol](double x,double y, double z) {
            return wtol;
        };
        
        do
        {
            DVec vecCG;
            DVec vecUz;
            vecCG.create_vector(pMesh,DVEC_TYPE::OCT_SHARED_NODES,DVEC_LOC::HOST,numVars,true);
            vecCG.create_vector(pMesh,DVEC_TYPE::OCT_LOCAL_WITH_PADDING,DVEC_LOC::HOST,numVars,true);
            
            DendroScalar* zipIn[numVars];
            DendroScalar* unzipIn[numVars];

            vecCG.to_2d(zipIn);
            vecUz.to_2d(unzipIn);

            if(pMesh->isActive())
            {
                
                unsigned int nodeLookUp_CG;
                unsigned int nodeLookUp_DG;
                double x,y,z,len;
                const ot::TreeNode * pNodes=&(*(pMesh->getAllElements().begin()));
                unsigned int ownerID,ii_x,jj_y,kk_z;
                unsigned int eleOrder=pMesh->getElementOrder();
                const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
                const unsigned int * e2n_dg = &(*(pMesh->getE2NMapping_DG().begin()));
                const unsigned int nPe = pMesh->getNumNodesPerElement();
                const unsigned int nodeLocalBegin = pMesh->getNodeLocalBegin();
                const unsigned int nodeLocalEnd   = pMesh->getNodeLocalEnd();


                double var[numVars];

                for(unsigned int elem=pMesh->getElementLocalBegin();elem<pMesh->getElementLocalEnd();elem++)
                {

                    for(unsigned int k=0;k<(eleOrder+1);k++)
                        for(unsigned int j=0;j<(eleOrder+1);j++ )
                            for(unsigned int i=0;i<(eleOrder+1);i++)
                            {
                                nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                                {
                                    nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                    pMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                                    len= (double) (1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                                    x=pNodes[ownerID].getX()+ ii_x*(len/(double)(eleOrder));
                                    y=pNodes[ownerID].getY()+ jj_y*(len/(double)(eleOrder));
                                    z=pNodes[ownerID].getZ()+ kk_z*(len/(double)(eleOrder));
                                    
                                    func((double)x,(double)y,(double)z,var);
                                    for(unsigned int v=0;v< numVars ; v++)
                                        zipIn[v][nodeLookUp_CG]=var[v];

                                }
                            }
                    
                }


                pMesh->readFromGhostBegin(vecCG.get_vec_ptr(),numVars);
                pMesh->readFromGhostEnd(vecCG.get_vec_ptr(),numVars);

                pMesh->unzip(vecCG.get_vec_ptr(), vecUz.get_vec_ptr(), numVars);

            }

            
            std::vector<unsigned int> refine_flags;
            refine_flags.resize(pMesh->getNumLocalMeshElements());
            isRefine1 = wavelet::compute_wavelet_remesh_flags(pMesh, refine_flags, (const DendroScalar**) unzipIn , rid, nr ,waveletTolFunc,0.1,true);

            if(isRefine1)
            {
                
                
                pMesh->setMeshRefinementFlags(refine_flags);
                ot::Mesh* newMesh = pMesh->ReMesh();
                
                DendroIntL localSz = pMesh->getNumLocalMeshElements();
                DendroIntL gSz_new, gSz_old;

                par::Mpi_Reduce(&localSz,&gSz_old,1,MPI_SUM,0,gcomm);
                localSz = newMesh->getNumLocalMeshElements();
                par::Mpi_Reduce(&localSz,&gSz_new,1,MPI_SUM,0,gcomm);
                if(!rank)
                    std::cout<<"[WAMR Grid Convergence]  old mesh size: "<<gSz_old<<" new mesh size: "<<gSz_new<<std::endl;

                std::swap(newMesh,pMesh);
                delete newMesh;

                iter++; 
            }

            vecCG.destroy_vector();
            vecUz.destroy_vector();

        }while(isRefine1 && iter < maxiter);
        
        return ;

    }

    void computeBlockUnzipGhostNodes(const ot::Mesh* pMesh, unsigned int blk, std::vector<unsigned int>& gid)
    {

        if(pMesh->isActive())
        {
            const ot::TreeNode * pNodes= &(*(pMesh->getAllElements().begin()));

            const unsigned int nodeLocalBegin = pMesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd = pMesh->getNodeLocalEnd();

            const unsigned int* e2n_cg = &(*(pMesh->getE2NMapping().begin()));
            const unsigned int* e2e = &(*(pMesh->getE2EMapping().begin()));

            const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
            const unsigned int nPe = pMesh->getNumNodesPerElement();

            unsigned int lookup,node_cg;
            unsigned int child[NUM_CHILDREN];

            if(blk> blkList.size()) 
                return;

            const unsigned int pWidth = blkList[blk].get1DPadWidth();
            const ot::TreeNode blkNode = blkList[blk].getBlockNode();
            const unsigned int regLevel = blkList[blk].getRegularGridLev();

            for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
            {
                const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                const unsigned int emin = 0;
                const unsigned int emax = (1u<<(regLevel-blkNode.getLevel()))-1;

                if(((ei == emin) || (ej == emin) || (ek== emin) || (ei==emax) || (ej==emax) || (ek==emax)))
                {
                    // this is block boundary element. 
                    for(unsigned int node =0; node < nPe; node++)
                    {
                        node_cg = e2n_cg[elem*nPe + node];
                        if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                            gid.push_back(node_cg); 
                    }


                }

                if(pWidth > 0)
                {   
                    // we need to look for the boundary neigbours only when the padding width is > 0 . 
                    if(ei==emin)
                    { 
                        // OCT_DIR_LEFT
                        lookup = e2e[elem*NUM_FACES + OCT_DIR_LEFT];
                        if(lookup!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookup].getLevel() > regLevel )
                            {  
                                // neighbour octant is smaller. 
                                assert(pNodes[lookup].getLevel()==(regLevel+1));
                                //child.resize(NUM_CHILDREN,LOOK_UP_TABLE_DEFAULT);
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                child[1]=lookup;
                                child[3]=e2e[child[1]*NUM_FACES+OCT_DIR_UP];
                                assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                                child[5]=e2e[child[1]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                                child[7]=e2e[child[3]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                child[0]=LOOK_UP_TABLE_DEFAULT; 
                                child[2]=LOOK_UP_TABLE_DEFAULT; 
                                child[4]=LOOK_UP_TABLE_DEFAULT; 
                                child[6]=LOOK_UP_TABLE_DEFAULT; 

                                for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                {
                                    if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                    {
                                        for(unsigned int node =0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[child[c]*nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                            gid.push_back(node_cg);
                                        }
                                
                                    }
                                }
                        
                            }else
                            {   
                                // neighbour octant is same lev or coarser
                                assert(pNodes[lookup].getLevel() <= regLevel );
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }


                            }

                        }

                            
                    }

                    if(ei==emax)
                    { // OCT_DIR_RIGHT

                        lookup = e2e[elem*NUM_FACES + OCT_DIR_RIGHT];
                        if(lookup!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookup].getLevel() > regLevel )
                            { // neighbour octant is smaller. 
                                
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                child[0]=lookup;
                                child[2]=e2e[child[0]*NUM_FACES+OCT_DIR_UP];
                                assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                                child[4]=e2e[child[0]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                                child[6]=e2e[child[2]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[6]!=LOOK_UP_TABLE_DEFAULT);

                                child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_RIGHT];

                                for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                {
                                    if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                    {
                                        for(unsigned int node =0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[child[c]*nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                            gid.push_back(node_cg);
                                        }
                                
                                    }
                                }
                                

                            }else
                            {   
                                // neighbour octant is same lev or coarser
                                assert(pNodes[lookup].getLevel() <= regLevel );
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                            }

                        }
                    }

                    if(ej==emin)
                    {   // OCT_DIR_DOWN
                        lookup = e2e[elem*NUM_FACES + OCT_DIR_DOWN];
                        if(lookup!=LOOK_UP_TABLE_DEFAULT)
                        {

                            if(pNodes[lookup].getLevel() > regLevel )
                            { // neighbour octant is smaller. 

                                child[2]=lookup;
                                child[3]=e2e[child[2]*NUM_FACES+OCT_DIR_RIGHT];
                                assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                                child[6]=e2e[child[2]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                                child[7]=e2e[child[3]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_DOWN];

                                for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                {
                                    if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                    {
                                        for(unsigned int node =0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[child[c]*nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                gid.push_back(node_cg);
                                        }
                                
                                    }
                                }

                            }else
                            { // neighbour octant is same lev or coarser
                                
                                assert(pNodes[lookup].getLevel() <= regLevel );
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                            }

                        }
                    }

                    if(ej==emax)
                    {   // OCT_DIR_UP
                        lookup = e2e[elem*NUM_FACES + OCT_DIR_UP];
                        if(lookup!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookup].getLevel() > regLevel )
                            { // neighbour octant is smaller. 

                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                child[0]=lookup;
                                child[1]=e2e[child[0]*NUM_FACES+OCT_DIR_RIGHT];
                                assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                                child[4]=e2e[child[0]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                                child[5]=e2e[child[1]*NUM_FACES+OCT_DIR_FRONT];
                                assert(child[5]!=LOOK_UP_TABLE_DEFAULT);

                                child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                                child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];

                                for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                {
                                    if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                    {
                                        for(unsigned int node =0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[child[c]*nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                gid.push_back(node_cg);
                                        }
                                
                                        }
                                }

                            }else
                            { // neighbour octant is same lev or coarser

                                assert(pNodes[lookup].getLevel() <= regLevel );
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                            }
                        }
                    }


                    if(ek==emin)
                    {   // OCT_DIR_BACK
                        lookup = e2e[elem*NUM_FACES + OCT_DIR_BACK];
                        if(lookup!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookup].getLevel() > regLevel )
                            { // neighbour octant is smaller. 

                                child[4]=lookup;
                                child[5]=e2e[child[4]*NUM_FACES+OCT_DIR_RIGHT];
                                assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                                child[6]=e2e[child[4]*NUM_FACES+OCT_DIR_UP];
                                assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                                child[7]=e2e[child[5]*NUM_FACES+OCT_DIR_UP];
                                assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_BACK];
                                child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_BACK];
                                child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_BACK];
                                child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_BACK];

                                for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                {
                                    if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                    {
                                        for(unsigned int node =0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[child[c]*nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                gid.push_back(node_cg);
                                        }
                                
                                        }
                                }

                            }else
                            { // neighbour octant is same lev or coarser

                                assert(pNodes[lookup].getLevel() <= regLevel );
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                            }


                        }
                    }

                    if(ek==emax)
                    {  // OCT_DIR_FRONT
                        lookup = e2e[elem*NUM_FACES + OCT_DIR_FRONT];
                        if(lookup!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookup].getLevel() > regLevel )
                            { // neighbour octant is smaller. 

                                child[0]=lookup;
                                child[1]=e2e[child[0]*NUM_FACES+OCT_DIR_RIGHT];
                                assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                                child[2]=e2e[child[0]*NUM_FACES+OCT_DIR_UP];
                                assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                                child[3]=e2e[child[1]*NUM_FACES+OCT_DIR_UP];
                                assert(child[3]!=LOOK_UP_TABLE_DEFAULT);

                                child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];

                                for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                {
                                    if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                    {
                                        for(unsigned int node =0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[child[c]*nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                gid.push_back(node_cg);
                                        }
                                
                                    }
                                }

                            }else
                            { // neighbour octant is same lev or coarser

                                assert(pNodes[lookup].getLevel() <= regLevel );
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                            }
                        }
                    }
                
                }




            }

            // now look for edge neighbors and vertex neighbors of the block, this is only needed when the padding width is >0
            if(pWidth>0)
            {
                const std::vector<unsigned int> blk2Edge_map = blkList[blk].getBlk2DiagMap_vec();
                const std::vector<unsigned int> blk2Vert_map = blkList[blk].getBlk2VertexMap_vec();
                const unsigned int blk_ele_1D = blkList[blk].getElemSz1D();

                for(unsigned int edir =0; edir < NUM_EDGES; edir++)
                {

                    for(unsigned int k=0; k< blk_ele_1D; k++)
                    {

                        if(blk2Edge_map[edir*(2*blk_ele_1D) +  2*k]!= LOOK_UP_TABLE_DEFAULT)
                        {
                            if(blk2Edge_map[edir*(2*blk_ele_1D) + 2*k+0] == blk2Edge_map[edir*(2*blk_ele_1D) + 2*k+1])
                            {
                                lookup = blk2Edge_map[edir*(2*blk_ele_1D) + 2*k+0];
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                            }else
                            {
                                lookup = blk2Edge_map[edir*(2*blk_ele_1D) + 2*k+0];
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                                lookup = blk2Edge_map[edir*(2*blk_ele_1D) + 2*k+1];
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        gid.push_back(node_cg);
                                }

                            }
                        }
                        
                    }

                }

                
                for(unsigned int k=0; k < blk2Vert_map.size(); k++)
                {
                    lookup = blk2Vert_map[k];
                    if(lookup!=LOOK_UP_TABLE_DEFAULT)
                    {
                        for(unsigned int node = 0; node < nPe; node++)
                        {
                            node_cg = e2n_cg[lookup * nPe + node];
                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                gid.push_back(node_cg);
                        }
                    }
                }



            }

            
        }
    }

    void computeBlockUnzipDepElements(const ot::Mesh* pMesh, unsigned int blk, std::vector<unsigned int>& eid)
    {
        pMesh->blkUnzipElementIDs(blk,eid);        
        return;
    }

    unsigned int computeTLevel(const ot::Mesh* const pMesh, unsigned int blk)
    {
        if(!(pMesh->isActive()))
            return 0;
        
        const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
        const unsigned int * const e2e = pMesh->getE2EMapping().data();
        const ot::TreeNode* const pNodes = pMesh->getAllElements().data();
        const ot::Block block = blkList[blk];

        unsigned int b_lev = pNodes[block.getLocalElementBegin()].getLevel();
        unsigned int eijk[3];

        const std::vector<unsigned int> blkEdgeMap = block.getBlk2DiagMap_vec();
        const std::vector<unsigned int> blkVertMap = block.getBlk2VertexMap_vec();

        for(unsigned int i=0; i < blkEdgeMap.size(); i++)
        {
            if( blkEdgeMap[i]!=LOOK_UP_TABLE_DEFAULT && (b_lev < pNodes[blkEdgeMap[i]].getLevel()))
            {
                b_lev = pNodes[blkEdgeMap[i]].getLevel();
                return b_lev; // return since 2:1 balance octree there for the smallest neighbour octant. 
            }
                
        }

        for(unsigned int i=0; i < blkVertMap.size(); i++)
        {
            if( blkVertMap[i]!=LOOK_UP_TABLE_DEFAULT && (b_lev < pNodes[blkVertMap[i]].getLevel()))
            {
                b_lev = pNodes[blkVertMap[i]].getLevel();
                return b_lev; // return since 2:1 balance octree there for the smallest neighbour octant. 
            }
                
        }

        for(unsigned int ele = block.getLocalElementBegin(); ele < block.getLocalElementEnd(); ele ++)
        {
            if(block.isBlockInternalEle(pNodes[ele]))
                continue;

            for(unsigned int dir =0; dir< NUM_FACES; dir++)
            {
                const unsigned int lookUp = e2e[ ele * NUM_FACES + dir ];
                if( lookUp !=LOOK_UP_TABLE_DEFAULT  && (b_lev < pNodes[lookUp].getLevel()) )
                {
                    b_lev = pNodes[lookUp].getLevel();
                    return b_lev; // return since 2:1 balance octree there for the smallest neighbour octant.  
                }
                    
            }
            
        }

        return b_lev;
        
    }

    int slice_mesh(const ot::Mesh* pMesh, unsigned int s_val[3], unsigned int s_normal[3], std::vector<unsigned int>& sids)
    {
        if(!(pMesh->isActive()))
            return 0;

        const unsigned int nLocalElements = pMesh->getNumLocalMeshElements();
        sids.reserve(nLocalElements);

        const ot::TreeNode* const pNodes = pMesh->getAllElements().data();
        int pvec[3];

        for(unsigned int e=pMesh->getElementLocalBegin(); e < pMesh->getElementLocalEnd(); e++ )
        {
            pvec[0] = pNodes[e].minX() - s_val[0];
            pvec[1] = pNodes[e].minY() - s_val[1];
            pvec[2] = pNodes[e].minZ() - s_val[2];

            const int dot_p = pvec[0]*s_normal[0] + pvec[1]*s_normal[1] + pvec[2]*s_normal[2];
            if(dot_p == 0)
                sids.push_back(e);
        }

        return 0;


    }

    Mesh* createSplitMesh(unsigned int eleOrder, unsigned int lmin, unsigned int lmax,MPI_Comm comm)
    {
        std::vector<ot::TreeNode> tmpOct;
        createRegularOctree(tmpOct,lmin,m_uiDim,m_uiMaxDepth,comm);

        std::cout<<"tmpOct: "<<tmpOct.size()<<std::endl;

        const unsigned int mid_x = (1u<<(m_uiMaxDepth-1));
        std::vector<ot::TreeNode> rfTmpOct;
        
        for(unsigned int l = lmin; l < lmax ; l++)
        {
            
            for(unsigned int ele =0; ele < tmpOct.size(); ele++)
            {
                if(tmpOct[ele].minX() < mid_x /*&& tmpOct[ele].minY() < mid_x && tmpOct[ele].minZ() < mid_x*/ )
                    tmpOct[ele].addChildren(rfTmpOct);
                else
                    rfTmpOct.push_back(tmpOct[ele]);
               
            }

            std::swap(rfTmpOct,tmpOct);
            rfTmpOct.clear();

            std::cout<<"tmpOct: at lev "<<l<<" size : "<<tmpOct.size()<<std::endl;

        }


        return createMesh(tmpOct.data(), tmpOct.size(),eleOrder,comm);
        

    }


}// end of namespace ot. 