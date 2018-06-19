//
// Created by milinda on 9/2/16.
//

/**
 * @author: Milinda Shayamal Fernando
 * School of Computing , University of Utah
 *
 * @breif Contains the functions to generate the mesh data structure from the 2:1 balanced linear octree.
 *
 * Assumptions:
 * 1). Assumes that octree is balanced and sorted.
 * 2). Assumes that there is no duplicate nodes.
 *
 *
 * Communicator swithch assumptions.
 *
 ** global rank=0 is always active.
 ** communicator split should be always contigous.
 *
 *
 *
 * */


#include "mesh.h"
double t_e2e; // e2e map generation time
double t_e2n; // e2n map generation time
double t_sm; // sm map generation time
double t_blk; // perform blk setup time.

double t_e2e_g[3];
double t_e2n_g[3];
double t_sm_g[3];
double t_blk_g[3];


//#define DEBUG_E2N_MAPPING_SM
//#define DEBUG_MESH_GENERATION

namespace ot {



    Mesh::Mesh(std::vector<ot::TreeNode> &in, unsigned int k_s, unsigned int pOrder,unsigned int activeNpes,MPI_Comm comm,unsigned int grainSz,double ld_tol,unsigned int sf_k)
    {
        m_uiCommGlobal=comm;
        MPI_Comm_rank(m_uiCommGlobal,&m_uiGlobalRank);
        MPI_Comm_size(m_uiCommGlobal,&m_uiGlobalNpes);

        DendroIntL localSz=in.size();
        DendroIntL globalSz;

        par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,m_uiCommGlobal);

        m_uiIsActive=m_uiGlobalRank<activeNpes;

        par::splitComm2way(m_uiIsActive,&m_uiCommActive,m_uiCommGlobal);


        m_uiMeshDomain_min=0;
        m_uiMeshDomain_max=(1u<<(m_uiMaxDepth));

        m_uiStensilSz=k_s;
        m_uiElementOrder=pOrder;
        m_uiEL_i=0;

        if(m_uiDim==2)
            m_uiNpE=(m_uiElementOrder+1)*(m_uiElementOrder+1);
        else if(m_uiDim==3)
            m_uiNpE=(m_uiElementOrder+1)*(m_uiElementOrder+1)*(m_uiElementOrder+1);

        m_uiNumDirections=(1u << m_uiDim) - 2 * (m_uiDim - 2);

        m_uiRefEl=RefElement(m_uiDim,m_uiElementOrder);

        if(!m_uiIsActive)
        { // set internal data structures for inactive mesh instances.

            assert(in.size()==0);
            if(in.size()!=0)
            {
                std::cout<<"[COMM Shrink/Expansion error]: global_rank: "<<m_uiGlobalRank<<" is active: "<<m_uiIsActive<<" balOct sz: "<<in.size()<<std::endl;
                exit(0);
            }


            m_uiElementPreGhostBegin=0;
            m_uiElementPreGhostEnd=0;
            m_uiElementLocalBegin=0;
            m_uiElementLocalEnd=0;
            m_uiElementPostGhostBegin=0;
            m_uiElementPostGhostEnd=0;

            m_uiNumPreGhostElements=m_uiElementPreGhostEnd-m_uiElementPreGhostBegin;
            m_uiNumLocalElements=m_uiElementLocalEnd=m_uiElementLocalBegin;
            m_uiNumPostGhostElements=m_uiElementPostGhostEnd-m_uiElementPostGhostBegin;

            m_uiNodePreGhostBegin=0;
            m_uiNodePreGhostEnd=0;
            m_uiNodeLocalBegin=0;
            m_uiNodeLocalEnd=0;
            m_uiNodePostGhostBegin=0;
            m_uiNodePostGhostEnd=0;

            m_uiNumActualNodes=0;
            m_uiUnZippedVecSz=0;

            m_uiAllElements.clear();
            m_uiLocalSplitterElements.clear();
            m_uiAllLocalNode.clear();
            m_uiE2NMapping_CG.clear();
            m_uiE2NMapping_DG.clear();
            m_uiE2EMapping.clear();
            m_uiLocalBlockList.clear();






        }else {

            MPI_Comm_rank(m_uiCommActive, &m_uiActiveRank);
            MPI_Comm_size(m_uiCommActive, &m_uiActiveNpes);

            if (!m_uiActiveRank)
                std::cout << " [MPI_COMM_SWITCH]: Selected comm.size: " << m_uiActiveNpes << std::endl;

            if (in.size() <= 1) {
                std::cout << "rank: " << m_uiActiveRank << " input octree of size " << in.size()
                          << " is too small for the current comm.  " << std::endl;
                exit(0);
            }



            /*m_uiElementOrder=pOrder;
            m_uiRefEl=RefElement(1,m_uiElementOrder);*/

            double t_e2e_begin = MPI_Wtime();
            if (m_uiActiveNpes > 1)buildE2EMap(in, m_uiCommActive);
            else buildE2EMap(in);
            double t_e2e_end = MPI_Wtime();
            t_e2e = t_e2e_end - t_e2e_begin;

            double t_e2n_begin = MPI_Wtime();
            buildE2NMap();
            double t_e2n_end = MPI_Wtime();
            t_e2n = t_e2n_end - t_e2n_begin;

            double t_sm_begin = MPI_Wtime();
            //if(m_uiActiveNpes>1)computeNodeScatterMaps(m_uiCommActive);
            //if(m_uiActiveNpes>1)computeNodalScatterMap(m_uiCommActive);
            //if(m_uiActiveNpes>1)computeNodalScatterMap1(m_uiCommActive);
            //if(m_uiActiveNpes>1)computeNodalScatterMap2(m_uiCommActive);
            if (m_uiActiveNpes > 1)computeNodalScatterMap4(m_uiCommActive);

            double t_sm_end = MPI_Wtime();
            t_sm = t_sm_end - t_sm_begin;

            double t_blk_begin = MPI_Wtime();
            performBlocksSetup();
            double t_blk_end = MPI_Wtime();
            t_blk = t_blk_end - t_blk_begin;


            par::Mpi_Reduce(&t_e2e, t_e2e_g, 1, MPI_MIN, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_e2e, t_e2e_g + 1, 1, MPI_SUM, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_e2e, t_e2e_g + 2, 1, MPI_MAX, 0, m_uiCommActive);
            t_e2e_g[1] /= (double) m_uiActiveNpes;

            par::Mpi_Reduce(&t_e2n, t_e2n_g, 1, MPI_MIN, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_e2n, t_e2n_g + 1, 1, MPI_SUM, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_e2n, t_e2n_g + 2, 1, MPI_MAX, 0, m_uiCommActive);
            t_e2n_g[1] /= (double) m_uiActiveNpes;

            par::Mpi_Reduce(&t_sm, t_sm_g, 1, MPI_MIN, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_sm, t_sm_g + 1, 1, MPI_SUM, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_sm, t_sm_g + 2, 1, MPI_MAX, 0, m_uiCommActive);
            t_sm_g[1] /= (double) m_uiActiveNpes;

            par::Mpi_Reduce(&t_blk, t_blk_g, 1, MPI_MIN, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_blk, t_blk_g + 1, 1, MPI_SUM, 0, m_uiCommActive);
            par::Mpi_Reduce(&t_blk, t_blk_g + 2, 1, MPI_MAX, 0, m_uiCommActive);
            t_blk_g[1] /= (double) m_uiActiveNpes;

            for(unsigned int p=0;p<m_uiActiveNpes;p++)
            {
                if(m_uiSendNodeCount[p]!=0)
                    m_uiSendProcList.push_back(p);

                if(m_uiRecvNodeCount[p]!=0)
                    m_uiRecvProcList.push_back(p);

            }

            m_uiSendBufferNodes.resize(m_uiSendNodeOffset[m_uiActiveNpes-1]+m_uiSendNodeCount[m_uiActiveNpes-1]);
            m_uiRecvBufferNodes.resize(m_uiRecvNodeOffset[m_uiActiveNpes-1]+m_uiRecvNodeCount[m_uiActiveNpes-1]);

            // release comm counter memory
            if (m_uiActiveNpes > 1) {


                delete[] m_uiSplitterNodes;
                delete[] m_uiSendKeyCount;
                delete[] m_uiSendKeyOffset;
                delete[] m_uiSendOctCountRound1;
                delete[] m_uiSendOctOffsetRound1;
                delete[] m_uiSendOctCountRound2;
                delete[] m_uiSendOctOffsetRound2;


                delete[] m_uiSendKeyDiagCount;
                delete[] m_uiRecvKeyDiagCount;
                delete[] m_uiSendKeyDiagOffset;
                delete[] m_uiRecvKeyDiagOffset;

                delete[] m_uiSendOctCountRound1Diag;
                delete[] m_uiRecvOctCountRound1Diag;
                delete[] m_uiSendOctOffsetRound1Diag;
                delete[] m_uiRecvOctOffsetRound1Diag;


                delete[] m_uiRecvKeyCount;
                delete[] m_uiRecvKeyOffset;
                delete[] m_uiRecvOctCountRound1;
                delete[] m_uiRecvOctOffsetRound1;
                delete[] m_uiRecvOctCountRound2;
                delete[] m_uiRecvOctOffsetRound2;

                // note: @milinda I have moved sencNode count and offsets to std::vector<unsigned int > because we need them to perform the ghost exchange.

            }

        }
    }


    Mesh::Mesh(std::vector<ot::TreeNode> &in, unsigned int k_s, unsigned int pOrder,MPI_Comm comm,unsigned int grainSz,double ld_tol,unsigned int sf_k)
    {

        m_uiCommGlobal=comm;
        MPI_Comm_rank(m_uiCommGlobal,&m_uiGlobalRank);
        MPI_Comm_size(m_uiCommGlobal,&m_uiGlobalNpes);

        DendroIntL localSz=in.size();
        DendroIntL globalSz;

        par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,m_uiCommGlobal);
        int p_npes_prev=binOp::getPrevHighestPowerOfTwo((globalSz/grainSz));
        int p_npes_next=binOp::getNextHighestPowerOfTwo((globalSz/grainSz));
        
        int p_npes=globalSz/grainSz;
        //if(!m_uiGlobalRank) std::cout<<"p_npes_prev: "<<p_npes_prev<<" p_npes_next: "<<p_npes_next<<" p_npes: "<<p_npes<<" diff1: "<<std::abs(p_npes_prev-p_npes)<<" diff2: "<<std::abs(p_npes_next-p_npes)<<std::endl;
        (std::abs(p_npes_prev-p_npes)<=std::abs(p_npes_next-p_npes)) ? p_npes=p_npes_prev : p_npes=p_npes_next;
        
        if(p_npes>m_uiGlobalNpes) p_npes=m_uiGlobalNpes;
        // quick fix to enforce the npes>=2 for any given grain size.
        if(p_npes<=1 && m_uiGlobalNpes>1) p_npes=2;
        if(p_npes==m_uiGlobalNpes)
        {
            //m_uiCommActive=m_uiCommGlobal; // note : use MPI_Comm_dup which is more safe than the assignment operator. (and MPI_Comm_free is always possible for m_uiCommActive)
            MPI_Comm_dup(m_uiCommGlobal,&m_uiCommActive);
            m_uiIsActive=true;
        }else
        {
            assert(p_npes<m_uiGlobalNpes);
            //m_uiIsActive=(m_uiGlobalRank<(globalSz/grainSz));
            //m_uiIsActive=(m_uiGlobalRank<p_npes);
            m_uiIsActive=isRankSelected(m_uiGlobalNpes,m_uiGlobalRank,p_npes);
            par::splitComm2way(m_uiIsActive,&m_uiCommActive,m_uiCommGlobal);
        }

        shrinkOrExpandOctree(in,ld_tol,sf_k,m_uiIsActive,m_uiCommActive,m_uiCommGlobal);

        m_uiMeshDomain_min=0;
        m_uiMeshDomain_max=(1u<<(m_uiMaxDepth));

        m_uiStensilSz=k_s;
        m_uiElementOrder=pOrder;
        m_uiEL_i=0;

        if(m_uiDim==2)
            m_uiNpE=(m_uiElementOrder+1)*(m_uiElementOrder+1);
        else if(m_uiDim==3)
            m_uiNpE=(m_uiElementOrder+1)*(m_uiElementOrder+1)*(m_uiElementOrder+1);

        m_uiNumDirections=(1u << m_uiDim) - 2 * (m_uiDim - 2);

        m_uiRefEl=RefElement(m_uiDim,m_uiElementOrder);

        if(!m_uiIsActive)
        { // set internal data structures for inactive mesh instances.

            assert(in.size()==0);
            if(in.size()!=0)
            {
                std::cout<<"[COMM Shrink/Expansion error]: global_rank: "<<m_uiGlobalRank<<" is active: "<<m_uiIsActive<<" balOct sz: "<<in.size()<<std::endl;
                exit(0);
            }


            m_uiElementPreGhostBegin=0;
            m_uiElementPreGhostEnd=0;
            m_uiElementLocalBegin=0;
            m_uiElementLocalEnd=0;
            m_uiElementPostGhostBegin=0;
            m_uiElementPostGhostEnd=0;

            m_uiNumPreGhostElements=m_uiElementPreGhostEnd-m_uiElementPreGhostBegin;
            m_uiNumLocalElements=m_uiElementLocalEnd=m_uiElementLocalBegin;
            m_uiNumPostGhostElements=m_uiElementPostGhostEnd-m_uiElementPostGhostBegin;

            m_uiNodePreGhostBegin=0;
            m_uiNodePreGhostEnd=0;
            m_uiNodeLocalBegin=0;
            m_uiNodeLocalEnd=0;
            m_uiNodePostGhostBegin=0;
            m_uiNodePostGhostEnd=0;

            m_uiNumActualNodes=0;
            m_uiUnZippedVecSz=0;

            m_uiAllElements.clear();
            m_uiLocalSplitterElements.clear();
            m_uiAllLocalNode.clear();
            m_uiE2NMapping_CG.clear();
            m_uiE2NMapping_DG.clear();
            m_uiE2EMapping.clear();
            m_uiLocalBlockList.clear();






        }else
        {

            MPI_Comm_rank(m_uiCommActive, &m_uiActiveRank);
            MPI_Comm_size(m_uiCommActive, &m_uiActiveNpes);

            if(!m_uiActiveRank)
                std::cout<<" [MPI_COMM_SWITCH]: Selected comm.size: "<<m_uiActiveNpes<<std::endl;

            if(in.size()<=1)
            {
                std::cout<<"rank: "<<m_uiActiveRank<<" input octree of size "<<in.size()<<" is too small for the current comm.  "<<std::endl;
                exit(0);
            }



            /*m_uiElementOrder=pOrder;
            m_uiRefEl=RefElement(1,m_uiElementOrder);*/

            double t_e2e_begin=MPI_Wtime();
            if(m_uiActiveNpes>1 )buildE2EMap(in,m_uiCommActive);
            else buildE2EMap(in);
            double t_e2e_end=MPI_Wtime();
            t_e2e=t_e2e_end-t_e2e_begin;

            double t_e2n_begin=MPI_Wtime();
            buildE2NMap();
            double t_e2n_end=MPI_Wtime();
            t_e2n=t_e2n_end-t_e2n_begin;

            double t_sm_begin=MPI_Wtime();
            //if(m_uiActiveNpes>1)computeNodeScatterMaps(m_uiCommActive);
            //if(m_uiActiveNpes>1)computeNodalScatterMap(m_uiCommActive);
            //if(m_uiActiveNpes>1)computeNodalScatterMap1(m_uiCommActive);
            //if(m_uiActiveNpes>1)computeNodalScatterMap2(m_uiCommActive);
            if(m_uiActiveNpes>1)computeNodalScatterMap4(m_uiCommActive);

            double t_sm_end=MPI_Wtime();
            t_sm=t_sm_end-t_sm_begin;

            double t_blk_begin=MPI_Wtime();
            performBlocksSetup();
            double t_blk_end=MPI_Wtime();
            t_blk=t_blk_end-t_blk_begin;


            par::Mpi_Reduce(&t_e2e,t_e2e_g,1,MPI_MIN,0,m_uiCommActive);
            par::Mpi_Reduce(&t_e2e,t_e2e_g+1,1,MPI_SUM,0,m_uiCommActive);
            par::Mpi_Reduce(&t_e2e,t_e2e_g+2,1,MPI_MAX,0,m_uiCommActive);
            t_e2e_g[1]/=(double)m_uiActiveNpes;

            par::Mpi_Reduce(&t_e2n,t_e2n_g,1,MPI_MIN,0,m_uiCommActive);
            par::Mpi_Reduce(&t_e2n,t_e2n_g+1,1,MPI_SUM,0,m_uiCommActive);
            par::Mpi_Reduce(&t_e2n,t_e2n_g+2,1,MPI_MAX,0,m_uiCommActive);
            t_e2n_g[1]/=(double)m_uiActiveNpes;

            par::Mpi_Reduce(&t_sm,t_sm_g,1,MPI_MIN,0,m_uiCommActive);
            par::Mpi_Reduce(&t_sm,t_sm_g+1,1,MPI_SUM,0,m_uiCommActive);
            par::Mpi_Reduce(&t_sm,t_sm_g+2,1,MPI_MAX,0,m_uiCommActive);
            t_sm_g[1]/=(double)m_uiActiveNpes;

            par::Mpi_Reduce(&t_blk,t_blk_g,1,MPI_MIN,0,m_uiCommActive);
            par::Mpi_Reduce(&t_blk,t_blk_g+1,1,MPI_SUM,0,m_uiCommActive);
            par::Mpi_Reduce(&t_blk,t_blk_g+2,1,MPI_MAX,0,m_uiCommActive);
            t_blk_g[1]/=(double)m_uiActiveNpes;


            
            
            if(m_uiActiveNpes>1)
            {
                for(unsigned int p=0;p<m_uiActiveNpes;p++)
                {
                    if(m_uiSendNodeCount[p]!=0)
                        m_uiSendProcList.push_back(p);
                    
                    if(m_uiRecvNodeCount[p]!=0)
                        m_uiRecvProcList.push_back(p);
                }
                
                m_uiSendBufferNodes.resize(m_uiSendNodeOffset[m_uiActiveNpes-1]+m_uiSendNodeCount[m_uiActiveNpes-1]);
                m_uiRecvBufferNodes.resize(m_uiRecvNodeOffset[m_uiActiveNpes-1]+m_uiRecvNodeCount[m_uiActiveNpes-1]);
            }
            

            // release comm counter memory
            if(m_uiActiveNpes>1) {


                delete[] m_uiSplitterNodes;
                delete[] m_uiSendKeyCount;
                delete[] m_uiSendKeyOffset;
                delete[] m_uiSendOctCountRound1;
                delete[] m_uiSendOctOffsetRound1;
                delete[] m_uiSendOctCountRound2;
                delete[] m_uiSendOctOffsetRound2;




                delete [] m_uiSendKeyDiagCount;
                delete [] m_uiRecvKeyDiagCount;
                delete [] m_uiSendKeyDiagOffset;
                delete [] m_uiRecvKeyDiagOffset;

                delete [] m_uiSendOctCountRound1Diag;
                delete [] m_uiRecvOctCountRound1Diag;
                delete [] m_uiSendOctOffsetRound1Diag;
                delete [] m_uiRecvOctOffsetRound1Diag;




                delete[] m_uiRecvKeyCount;
                delete[] m_uiRecvKeyOffset;
                delete[] m_uiRecvOctCountRound1;
                delete[] m_uiRecvOctOffsetRound1;
                delete[] m_uiRecvOctCountRound2;
                delete[] m_uiRecvOctOffsetRound2;

                // note: @milinda I have moved sencNode count and offsets to std::vector<unsigned int > because we need them to perform the ghost exchange.

            }





        }


    }

    Mesh::~Mesh() {


        m_uiPreGhostOctants.clear();
        m_uiPostGhostOctants.clear();

        m_uiEmbeddedOctree.clear();
        m_uiGhostOctants.clear();
        m_uiAllElements.clear();

        m_uiE2NMapping_CG.clear();
        m_uiE2NMapping_DG.clear();
        m_uiE2EMapping.clear();

        m_uiSendBufferElement.clear();
        m_uiScatterMapElementRound1.clear();
        m_uiSendBufferNodes.clear();
        m_uiRecvBufferNodes.clear();
        m_uiScatterMapActualNodeSend.clear();
        m_uiScatterMapActualNodeRecv.clear();

        m_uiSendNodeOffset.clear();
        m_uiSendNodeCount.clear();
        m_uiRecvNodeOffset.clear();
        m_uiRecvNodeCount.clear();
        m_uiLocalSplitterElements.clear();

        m_uiSendProcList.clear();
        m_uiRecvProcList.clear();

        MPI_Comm_free(&m_uiCommActive);


    }


    void Mesh::generateSearchKeys()
    {
        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        std::vector<SearchKey> skeys;
        std::vector<SearchKey>::iterator hint;
        ot::TreeNode *inPtr = (&(*(m_uiEmbeddedOctree.begin())));
        unsigned int domain_max = 1u<<(m_uiMaxDepth);
        SearchKey skey;
        register unsigned int mySz;
        register unsigned int myX;
        register unsigned int myY;
        register unsigned int myZ;
        register unsigned int myLev;
        const unsigned int K=1;


        for (int i = 0; i < m_uiEmbeddedOctree.size(); i++) {
            myLev = m_uiEmbeddedOctree[i].getLevel();

            mySz = (1u << (m_uiMaxDepth - myLev));
            myX = inPtr[i].getX();
            myY = inPtr[i].getY();
            myZ = inPtr[i].getZ();



            /* Below orientation is used when generating keys.
            *
            * [up]
            * Y
            * |     Z [front]
            * |    /
            * |   /
            * |  /
            * | /
            * -------------> X [right]
                             */
            // Key generation along X axis.
            if ((myX + K * mySz) < domain_max) {
                hint = skeys.emplace(skeys.end(),SearchKey((myX + K * mySz), myY, myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT);

            }
            if (myX >0) {
                hint = skeys.emplace(skeys.end(),SearchKey((myX - 1), myY, myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT);
            }

            // Key generation along Y axis.
            if ((myY + K * mySz) < domain_max) {
                hint = skeys.emplace(skeys.end(),SearchKey(myX, (myY + K * mySz), myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP);

            }
            if (myY >0) {
                hint = skeys.emplace(skeys.end(),SearchKey(myX, (myY - 1), myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN);
            }

            if (m_uiDim == 3) {

                if ((myZ + K * mySz) < domain_max) {
                    hint = skeys.emplace(skeys.end(),SearchKey(myX, myY, (myZ + K * mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(i);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_FRONT);
                }


                if (myZ >0) {
                    hint = skeys.emplace(skeys.end(),SearchKey(myX, myY, (myZ - 1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(i);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_BACK);
                }

            }

        }


        if(m_uiActiveNpes>1) {
            for (unsigned int i = 0; i < 2*m_uiActiveNpes; i++) {
                skeys.emplace(skeys.end(),SearchKey(m_uiLocalSplitterElements[i]));
            }
        }

       SearchKey rootSkey(m_uiDim,m_uiMaxDepth);
       std::vector<SearchKey> tmpSKeys;
       SFC::seqSort::SFC_treeSort(&(*(skeys.begin())),skeys.size(),tmpSKeys,tmpSKeys,tmpSKeys,m_uiMaxDepth,m_uiMaxDepth,rootSkey,ROOT_ROTATION,1,TS_SORT_ONLY);
       assert(seq::test::isSorted(skeys));

       Key tmpKey;
       unsigned int skip=0;
       for(unsigned int e=0;e<(skeys.size());e++)
       {
           tmpKey=Key(skeys[e].getX(),skeys[e].getY(),skeys[e].getZ(),skeys[e].getLevel(),m_uiDim,m_uiMaxDepth);
           if(skeys[e].getOwner()>=0){
               tmpKey.addOwner(skeys[e].getOwner());
               tmpKey.addStencilIndexAndDirection(K-1,skeys[e].getStencilIndexDirectionList());
           }

           skip=1;
           while(((e+skip)<skeys.size()) && (skeys[e]==skeys[e+skip]))
           {
               if(skeys[e+skip].getOwner()>=0){
                   tmpKey.addOwner(skeys[e+skip].getOwner());
                   tmpKey.addStencilIndexAndDirection(K-1,skeys[e+skip].getStencilIndexDirectionList());
               }
               skip++;
           }

           m_uiKeys.push_back(tmpKey);
           e+=(skip-1);

       }


       skeys.clear();
       //if(!m_uiActiveRank) std::cout<<"key gen 0 ended "<<std::endl;

    }


    void Mesh::generateGhostElementSearchKeys()
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        std::vector<SearchKey> skeys;
        std::vector<SearchKey>::iterator hint;
        ot::TreeNode *inPtr = (&(*(m_uiAllElements.begin())));
        unsigned int domain_max = 1u<<(m_uiMaxDepth);
        SearchKey skey;
        register unsigned int mySz;
        register unsigned int myX;
        register unsigned int myY;
        register unsigned int myZ;
        register unsigned int myLev;
        const unsigned int K=1;

        for (unsigned int i = m_uiElementPreGhostBegin; i < m_uiElementPreGhostEnd; i++) {

            myLev = inPtr[i].getLevel();


            mySz= (1u << (m_uiMaxDepth - myLev));
            myX = inPtr[i].getX();
            myY = inPtr[i].getY();
            myZ = inPtr[i].getZ();
            //domain_max=(1u<<my)
            // We can skip the morton index -0  because that key is mapped to the current element. So No need to search for that.
            // Note: we do not need to perform any boundary checks when generating the keys because, we are skipping all the level one octatns.

            //for (unsigned int K = 1; K <= m_uiStensilSz; K++) {


            /** *
             * Below orientation is used when generating keys.
             *
             * [up]
             * Y
             * |     Z [front]
             * |    /
             * |   /
             * |  /
             * | /
             * -------------> X [right]
             */

            // Key generation along X axis.
            if ((myX + K * mySz) < domain_max) {
                hint = skeys.emplace(skeys.end(),SearchKey((myX + K * mySz), myY, myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT);

            }
            if (myX >0) {
                hint = skeys.emplace(skeys.end(),SearchKey((myX - 1), myY, myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT);

            }
            // Key generation along Y axis.
            if ((myY + K * mySz) < domain_max) {
                hint = skeys.emplace(skeys.end(),SearchKey(myX, (myY + K * mySz), myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP);

            }
            if (myY >0) {
                hint = skeys.emplace(skeys.end(),SearchKey(myX, (myY - 1), myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN);

            }

            if (m_uiDim == 3) {

                if ((myZ + K * mySz) < domain_max) {
                    hint = skeys.emplace(skeys.end(),SearchKey(myX, myY, (myZ + K * mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(i);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_FRONT);

                }
                if (myZ >0) {
                    hint = skeys.emplace(skeys.end(),SearchKey(myX, myY, (myZ - 1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(i);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_BACK);

                }

            }


        }//end for i



        for (unsigned int i = m_uiElementPostGhostBegin; i < m_uiElementPostGhostEnd; i++) {

            myLev = inPtr[i].getLevel();

            mySz = (1u << (m_uiMaxDepth - myLev));
            myX = inPtr[i].getX();
            myY = inPtr[i].getY();
            myZ = inPtr[i].getZ();
            //domain_max=(1u<<my)
            // We can skip the morton index -0  because that key is mapped to the current element. So No need to search for that.
            // Note: we do not need to perform any boundary checks when generating the keys because, we are skipping all the level one octatns.

            /** *
             * Below orientation is used when generating keys.
             *
             * [up]
             * Y
             * |     Z [front]
             * |    /
             * |   /
             * |  /
             * | /
             * -------------> X [right]
             */

            // Key generation along X axis.
            if ((myX + K * mySz) < domain_max) {
                hint = skeys.emplace(skeys.end(),SearchKey((myX + K * mySz), myY, myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT);


            }
            if (myX >0) {
                hint = skeys.emplace(skeys.end(),SearchKey((myX - 1), myY, myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT);


            }
            // Key generation along Y axis.
            if ((myY + K * mySz) < domain_max) {
                hint = skeys.emplace(skeys.end(),SearchKey(myX, (myY + K * mySz), myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP);

            }
            if (myY >0) {
                hint = skeys.emplace(skeys.end(),SearchKey(myX, (myY - 1), myZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(i);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN);


            }

            if (m_uiDim == 3) {

                if ((myZ + K * mySz) < domain_max) {
                    hint = skeys.emplace(skeys.end(),SearchKey(myX, myY, (myZ + K * mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(i);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_FRONT);

                }
                if (myZ >0) {
                    hint = skeys.emplace(skeys.end(),SearchKey(myX, myY, (myZ - 1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(i);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_BACK);

                }

            }

        }


        SearchKey rootSkey(m_uiDim,m_uiMaxDepth);
        std::vector<SearchKey> tmpSKeys;
        SFC::seqSort::SFC_treeSort(&(*(skeys.begin())),skeys.size(),tmpSKeys,tmpSKeys,tmpSKeys,m_uiMaxDepth,m_uiMaxDepth,rootSkey,ROOT_ROTATION,1,TS_SORT_ONLY);
        assert(seq::test::isSorted(skeys));


        //std::cout<<"rank: "<<m_uiActiveRank<<" skeys: "<<skeys.size()<<" pre local post: "<<m_uiElementPreGhostEnd<<" "<<m_uiElementLocalEnd<<" "<<m_uiElementPostGhostEnd<<std::endl;

        Key tmpKey;
        unsigned int skip=0;
        for(unsigned int e=0;e<(skeys.size());e++)
        {
            tmpKey=Key(skeys[e].getX(),skeys[e].getY(),skeys[e].getZ(),skeys[e].getLevel(),m_uiDim,m_uiMaxDepth);
            if(skeys[e].getOwner()>=0){
                tmpKey.addOwner(skeys[e].getOwner());
                tmpKey.addStencilIndexAndDirection(K-1,skeys[e].getStencilIndexDirectionList());
            }

            skip=1;
            while(((e+skip)<skeys.size()) && (skeys[e]==skeys[e+skip]))
            {
                if(skeys[e+skip].getOwner()>=0){
                    tmpKey.addOwner(skeys[e+skip].getOwner());
                    tmpKey.addStencilIndexAndDirection(K-1,skeys[e+skip].getStencilIndexDirectionList());
                }
                skip++;
            }

            m_uiGhostKeys.push_back(tmpKey);
            e+=(skip-1);

        }



        skeys.clear();
        //if(!m_uiActiveRank) std::cout<<"key gen 1 ended "<<std::endl;

    }


    void Mesh::generateBdyElementDiagonalSearchKeys()
    {


        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        std::vector<SearchKey> skeys;
        std::vector<SearchKey>::iterator hint;
        ot::TreeNode *inPtr = (&(*(m_uiAllElements.begin())));
        unsigned int domain_max = 1u<<(m_uiMaxDepth);
        SearchKey skey;
        unsigned int elementLookUp;
        register unsigned int mySz;
        register unsigned int myX;
        register unsigned int myY;
        register unsigned int myZ;
        register unsigned int myLev;
        const unsigned int K=1;


        // note : this is just to find the boundary of the local octree. boundary octants should be face 2 distance from the ghsot elemnts.
        std::vector<unsigned int > bdyID_L1;
        std::vector<unsigned int > bdyID_L2;

        for (unsigned int i = m_uiElementLocalBegin; i < m_uiElementLocalEnd; i++)
        {
            for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
            {
                elementLookUp=m_uiE2EMapping[i*m_uiNumDirections+dir];
                if((elementLookUp!=LOOK_UP_TABLE_DEFAULT) && ((elementLookUp<m_uiElementLocalBegin) ||  (elementLookUp>=m_uiElementLocalEnd))) {
                    bdyID_L1.push_back(i);
                    break;
                }
            }


        }

        for(unsigned int i=0;i<bdyID_L1.size();i++)
        {
            for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
            {

                elementLookUp=m_uiE2EMapping[bdyID_L1[i]*m_uiNumDirections+dir];
                if((elementLookUp>=m_uiElementLocalBegin) && (elementLookUp<m_uiElementLocalEnd)) {
                    bdyID_L2.push_back(elementLookUp);
                }
            }

        }

        // merger level 1 & level 2 bdy octants.
        bdyID_L1.insert(bdyID_L1.end(),bdyID_L2.begin(),bdyID_L2.end());

        // remove duplicates.
        std::sort(bdyID_L1.begin(),bdyID_L1.end());
        bdyID_L1.erase(std::unique(bdyID_L1.begin(),bdyID_L1.end()),bdyID_L1.end());


        const unsigned int * bdyID=&(*(bdyID_L1.begin()));

        for(unsigned int e=0;e<bdyID_L1.size();e++)
        {
            myLev = inPtr[bdyID[e]].getLevel();

            mySz= (1u << (m_uiMaxDepth - myLev));
            myX = inPtr[bdyID[e]].getX();
            myY = inPtr[bdyID[e]].getY();
            myZ = inPtr[bdyID[e]].getZ();

            // Edge keys

            if(myX>0 && myY>0)
            {


                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY-1),(myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_DOWN);

                if((myZ+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY-1),(myZ+((K*mySz)>>1u)), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_DOWN);

                    if((myZ+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY-1),(myZ+((K*mySz))), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_DOWN);
                    }

                }


            }



            if(myX>0 && (myY+K*mySz)<domain_max)
            {


                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+K*mySz),(myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_UP);

                if((myZ+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+K*mySz),(myZ+((K*mySz)>>1u)), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_UP);

                    if((myZ+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+K*mySz),(myZ+((K*mySz))), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_UP);
                    }

                }


            }



            if(myX>0 && myZ>0)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_BACK);

                if((myY+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+((K*mySz)>>1u)),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_BACK);

                    if((myY+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+((K*mySz))),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_BACK);
                    }

                }


            }



            if(myX>0 && (myZ+K*mySz)<domain_max)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_FRONT);


                if((myY+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+((K*mySz)>>1u)),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_FRONT);

                    if((myY+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+((K*mySz))),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_FRONT);
                    }

                }


            }



            if((myX+K*mySz) < domain_max && myY>0)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY-1),(myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_DOWN);

                if((myZ+((K*mySz)>>1u))<domain_max)
                {

                    hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY-1),(myZ+((K*mySz)>>1u)), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_DOWN);

                    if((myZ+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY-1),(myZ+((K*mySz))), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_DOWN);

                    }

                }



            }



            if((myX+K*mySz)<domain_max && (myY+K*mySz)<domain_max)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+K*mySz),(myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_UP);

                if((myZ+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+K*mySz),(myZ+((K*mySz)>>1u)), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_UP);

                    if((myZ+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+K*mySz),(myZ+((K*mySz))), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_UP);

                    }

                }

            }



            if((myX+K*mySz)<domain_max && myZ>0)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_BACK);

                if((myY+((K*mySz)>>1u))<domain_max)
                {

                    hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+((K*mySz)>>1u)),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_BACK);


                    if((myY+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+((K*mySz))),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_BACK);
                    }

                }


            }


            if((myX+K*mySz)<domain_max && (myZ+K*mySz)<domain_max)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_FRONT);

                if((myY+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+((K*mySz)>>1u)),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_FRONT);

                    if((myY+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+((K*mySz))),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_FRONT);
                    }

                }


            }


            if(myY>0 && myZ>0)
            {


                hint=skeys.emplace(skeys.end(),SearchKey((myX),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN_BACK);

                if((myX+((K*mySz)>>1u))<domain_max)
                {

                    hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz)>>1u)),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN_BACK);

                    if((myX+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz))),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN_BACK);
                    }

                }



            }


            if(myY > 0 && (myZ+K*mySz)<domain_max)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX),(myY-1),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN_FRONT);

                if((myX+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz)>>1u)),(myY-1),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN_FRONT);

                    if((myX+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz))),(myY-1),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_DOWN_FRONT);
                    }
                }


            }


            if((myY+K*mySz)<domain_max && myZ>0)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX),(myY+K*mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP_BACK);

                if((myX+((K*mySz)>>1u))<domain_max)
                {
                    hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz)>>1u)),(myY+K*mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP_BACK);

                    if((myX+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz))),(myY+K*mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP_BACK);

                    }

                }

            }



            if((myY+K*mySz)<domain_max && (myZ+K*mySz)<domain_max)
            {

                hint=skeys.emplace(skeys.end(),SearchKey((myX),(myY+K*mySz),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP_FRONT);

                if((myX+((K*mySz)>>1u))<domain_max)
                {

                    hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz)>>1u)),(myY+K*mySz),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                    hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                    hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP_FRONT);

                    if((myX+((K*mySz)))<domain_max)
                    {
                        hint=skeys.emplace(skeys.end(),SearchKey((myX+((K*mySz))),(myY+K*mySz),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                        hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                        hint->addStencilIndexAndDirection(K - 1, OCT_DIR_UP_FRONT);
                    }
                }


            }


            // Vertex Keys.
            if( (myX>0) && (myY>0) && (myZ>0) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_DOWN_BACK);

            }

            if( ((myX+K*mySz)<domain_max) && (myY>0) && (myZ>0) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_DOWN_BACK);

            }

            if( (myX>0) && ((myY+K*mySz)<domain_max) && (myZ>0) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+K*mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_UP_BACK);

            }

            if( ((myX+K*mySz)<domain_max) && ((myY+K*mySz)<domain_max) && (myZ>0) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+K*mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_UP_BACK);

            }


            if( (myX>0) && (myY>0) && ((myZ+K*mySz)<domain_max) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY-1),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_DOWN_FRONT);

            }

            if( ((myX+K*mySz)<domain_max) && (myY>0) && ((myZ+K*mySz)<domain_max) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY-1),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_DOWN_FRONT);

            }

            if( (myX>0) && ((myY+K*mySz)<domain_max) && ((myZ+K*mySz)<domain_max) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX-1),(myY+K*mySz),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_LEFT_UP_FRONT);

            }

            if( ((myX+K*mySz)<domain_max) && ((myY+K*mySz)<domain_max) && ((myZ+K*mySz)<domain_max) )
            {
                hint=skeys.emplace(skeys.end(),SearchKey((myX+K*mySz),(myY+K*mySz),(myZ+K*mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
                hint->addOwner(bdyID[e]-m_uiElementLocalBegin);
                hint->addStencilIndexAndDirection(K - 1, OCT_DIR_RIGHT_UP_FRONT);

            }






        }


       if(m_uiActiveNpes>1) {
           for (unsigned int i = 0; i < 2*m_uiActiveNpes; i++) {
               //tmpKey = Key(m_uiSplitterElementsWGhost[i].getX(), m_uiSplitterElementsWGhost[i].getY(), m_uiSplitterElementsWGhost[i].getZ(), m_uiMaxDepth,m_uiDim, m_uiMaxDepth);
               //tmpKey=Key(m_uiLocalSplitterElements[i]);
               skeys.emplace(skeys.end(),SearchKey(m_uiLocalSplitterElements[i]));
           }
       }

        SearchKey rootSkey(m_uiDim,m_uiMaxDepth);
        std::vector<SearchKey> tmpSKeys;
        SFC::seqSort::SFC_treeSort(&(*(skeys.begin())),skeys.size(),tmpSKeys,tmpSKeys,tmpSKeys,m_uiMaxDepth,m_uiMaxDepth,rootSkey,ROOT_ROTATION,1,TS_SORT_ONLY);
        assert(seq::test::isSorted(skeys));

        Key tmpKey;
        unsigned int skip=0;
        for(unsigned int e=0;e<(skeys.size());e++)
        {
            tmpKey=Key(skeys[e].getX(),skeys[e].getY(),skeys[e].getZ(),skeys[e].getLevel(),m_uiDim,m_uiMaxDepth);
            if(skeys[e].getOwner()>=0){
                tmpKey.addOwner(skeys[e].getOwner());
                tmpKey.addStencilIndexAndDirection(K-1,skeys[e].getStencilIndexDirectionList());
            }

            skip=1;
            while(((e+skip)<skeys.size()) && (skeys[e]==skeys[e+skip]))
            {
                if(skeys[e+skip].getOwner()>=0){
                    tmpKey.addOwner(skeys[e+skip].getOwner());
                    tmpKey.addStencilIndexAndDirection(K-1,skeys[e+skip].getStencilIndexDirectionList());
                }
                skip++;
            }

            m_uiKeysDiag.push_back(tmpKey);
            e+=(skip-1);

        }


       skeys.clear();
       //if(!m_uiActiveRank) std::cout<<"key gen 2 ended "<<std::endl;

    }


    void Mesh::buildE2EMap(std::vector<ot::TreeNode> &in,MPI_Comm comm)
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        int rank, npes;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &npes);


#ifdef DEBUG_MESH_GENERATION
        //treeNodesTovtk(in,rank,"balOct");
#endif
        std::swap(m_uiEmbeddedOctree, in);
        in.clear();

        // Below Key vector and ot::TreeNode vector is being use for sorting Keys and treeNodes repeatedly. So make sure you clear them after using them in SFC_TreeSort.
        std::vector<Key> tmpKeys;
        std::vector<ot::TreeNode> tmpNodes;

        Key rootKey(0, 0, 0, (OCT_KEY_NONE | 0), m_uiDim, m_uiMaxDepth);
        ot::TreeNode rootNode(0,0,0,0,m_uiDim,m_uiMaxDepth);

        assert(m_uiEmbeddedOctree.size()>1); // m_uiEmbedded octree cannot be empty.  (Remove this assertion once we handle this case. )
        assert(par::test::isUniqueAndSorted(m_uiEmbeddedOctree, comm));


#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(m_uiEmbeddedOctree,rank,"m_uiEmbeddedOctree");
#endif



        // AllGather of the the max local octant. This will be used to split the keys among the processors.
        m_uiLocalSplitterElements.resize(2*npes);
        ot::TreeNode localMinMaxElement[2];
        localMinMaxElement[0]=m_uiEmbeddedOctree.front();
        localMinMaxElement[1]=m_uiEmbeddedOctree.back();
        par::Mpi_Allgather(localMinMaxElement, &(*(m_uiLocalSplitterElements.begin())), 2, comm);


#ifdef DEBUG_MESH_GENERATION
        if(!rank)
        {
            std::vector<ot::TreeNode> splitter_vec;
            for(int i=0;i<npes;i++)
            {
                splitter_vec.push_back(m_uiLocalSplitterElements[i]);
            }
            treeNodesTovtk(splitter_vec,rank,"splitters");
        }
#endif


        // 2- generate the keys. (These are the keys for the local elements \tau_{loc})
       generateSearchKeys();

#ifdef DEBUG_MESH_GENERATION
        if(!rank) std::cout<<"Key generation completed "<<std::endl;
#endif


#ifdef DEBUG_MESH_GENERATION
       /* unsigned int localKeySz=keys_vec.size();
        unsigned int globalKeySz=0;
        par::Mpi_Reduce(&localKeySz,&globalKeySz,1,MPI_SUM,0,comm);
        if(!rank) std::cout<<" Total number of keys generated: "<<globalKeySz<<std::endl;*/
#endif



        //4- Compute Face Neighbors (By sending owners of the key to the correct proc. )===================================================================================================================================================

        SFC::seqSort::SFC_treeSort(&(*(m_uiKeys.begin())), m_uiKeys.size(), tmpKeys, tmpKeys, tmpKeys, m_uiMaxDepth,m_uiMaxDepth, rootKey, 0, 1, TS_SORT_ONLY);
        tmpKeys.clear();
        assert(seq::test::isUniqueAndSorted(m_uiKeys));


#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(m_uiKeys,rank,"m_uiKeys",false);
#endif


        m_uiSendKeyCount = new unsigned int[npes];
        m_uiRecvKeyCount = new unsigned int[npes];

        m_uiSendKeyOffset = new unsigned int[npes];
        m_uiRecvKeyOffset = new unsigned int[npes];


        m_uiSendOctCountRound1 = new unsigned int[npes];
        m_uiRecvOctCountRound1 = new unsigned int[npes];

        m_uiSendOctOffsetRound1 = new unsigned int[npes];
        m_uiRecvOctOffsetRound1 = new unsigned int[npes];

        Key * m_uiKeysPtr=&(*(m_uiKeys.begin()));


        std::vector<Key> splitterElements;
        splitterElements.resize(2*npes);
        for(unsigned int p=0;p<2*npes;p++)
           splitterElements[p]=Key(m_uiLocalSplitterElements[p]);//Key(m_uiLocalSplitterElements[p].getX(),m_uiLocalSplitterElements[p].getY(),m_uiLocalSplitterElements[p].getZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth);

        assert(seq::test::isUniqueAndSorted(splitterElements));


        // search element splitters in the keys, to determine who owns the keys.
        SFC::seqSearch::SFC_treeSearch(&(*(splitterElements.begin())),&(*(m_uiKeys.begin())),0,splitterElements.size(),0,m_uiKeys.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


        // send owners of the keys to the key owner. (R1-a ghost exchange.)
        m_uiGhostElementIDsToBeSent.clear();
        unsigned int sBegin=0;
        unsigned int sEnd;

        for(unsigned int p=0;p<npes;p++)
        {

            m_uiSendKeyCount[p]=0;
            assert((splitterElements[2*p].getFlag() & OCT_FOUND));
            assert((splitterElements[2*p+1].getFlag() & OCT_FOUND));;
            assert(m_uiKeys[splitterElements[2*p].getSearchResult()]==splitterElements[2*p]);
            assert(m_uiKeys[splitterElements[2*p+1].getSearchResult()]==splitterElements[2*p+1]);
            sBegin=splitterElements[2*p].getSearchResult();
            //sEnd=splitterElements[2*p+1].getSearchResult()+1;
            (p<(npes-1))? sEnd=splitterElements[2*p+2].getSearchResult()+1: sEnd=m_uiKeys.size();
            if(p!=m_uiActiveRank)
            for(unsigned int k=sBegin;k<sEnd;k++)
            {
                for (unsigned int w = 0; w < m_uiKeysPtr[k].getOwnerList()->size(); w++)
                    m_uiGhostElementIDsToBeSent.push_back((*(m_uiKeysPtr[k].getOwnerList()))[w]);

                 m_uiSendKeyCount[p]+=m_uiKeysPtr[k].getOwnerList()->size();

            }


        }


        par::Mpi_Alltoall(m_uiSendKeyCount,m_uiRecvKeyCount,1,comm);

        m_uiSendKeyOffset[0]=0;
        m_uiRecvKeyOffset[0]=0;
        omp_par::scan(m_uiSendKeyCount,m_uiSendKeyOffset,npes);
        omp_par::scan(m_uiRecvKeyCount,m_uiRecvKeyOffset,npes);

        //std::cout<<rank <<" ghostEleSz: "<<m_uiGhostElementIDsToBeSent.size()<<" computed Size: "<<(m_uiSendKeyOffset[npes-1]+m_uiSendKeyCount[npes-1])<<std::endl;
        assert(m_uiGhostElementIDsToBeSent.size()==(m_uiSendKeyOffset[npes-1]+m_uiSendKeyCount[npes-1]));

        // we need to collect both face edge & vertex neighbors as level 1 ghost layer.
        // here we collect the face neighbours later we will collect the edge and vertex neighbors.
        std::vector<unsigned int > * scatterMapSend_R1=new std::vector<unsigned int >[npes];

        std::set<unsigned  int > tmpSendElementIds;
        for(unsigned int p=0;p<npes;p++)
        {
            scatterMapSend_R1[p]=std::vector<unsigned int >();
            tmpSendElementIds.insert((m_uiGhostElementIDsToBeSent.begin()+m_uiSendKeyOffset[p]),(m_uiGhostElementIDsToBeSent.begin()+m_uiSendKeyCount[p]+m_uiSendKeyOffset[p]));
            m_uiSendOctCountRound1[p]=tmpSendElementIds.size();
            scatterMapSend_R1[p].insert(scatterMapSend_R1[p].end(),tmpSendElementIds.begin(),tmpSendElementIds.end());
            tmpSendElementIds.clear();
        }




#ifdef DEBUG_MESH_GENERATION
       /* MPI_Barrier(comm);
        if(!rank) std::cout<<"Ghost Element Size: "<<m_uiScatterMapElementRound1.size()<<" without removing duplicates: "<<m_uiGhostElementIDsToBeSent.size()<<std::endl;*/
#endif
        m_uiGhostElementIDsToBeSent.clear();


#ifdef DEBUG_MESH_GENERATION
         if(!rank)
         for(unsigned int k=0;k<npes;k++)
         {
             std::cout<<"Processor "<<k<<" sendCnt: "<<m_uiSendKeyCount[k]<<" recvCnt: "<<m_uiRecvKeyCount[k]<<std::endl;
         }
#endif

        par::Mpi_Alltoall(m_uiSendOctCountRound1,m_uiRecvOctCountRound1,1,comm);

        m_uiSendOctOffsetRound1[0]=0;
        m_uiRecvOctOffsetRound1[0]=0;

        omp_par::scan(m_uiSendOctCountRound1,m_uiSendOctOffsetRound1,npes);
        omp_par::scan(m_uiRecvOctCountRound1,m_uiRecvOctOffsetRound1,npes);


        for(unsigned int p=0;p<npes;p++)
            for(unsigned int k=0;k<scatterMapSend_R1[p].size();k++)
                m_uiSendBufferElement.push_back(m_uiEmbeddedOctree[scatterMapSend_R1[p][k]]);


        m_uiGhostOctants.resize(m_uiRecvOctOffsetRound1[npes - 1] + m_uiRecvOctCountRound1[npes - 1]);

        par::Mpi_Alltoallv(&(*(m_uiSendBufferElement.begin())), (int *) m_uiSendOctCountRound1, (int *) m_uiSendOctOffsetRound1,
                           &(*(m_uiGhostOctants.begin())), (int *) m_uiRecvOctCountRound1, (int *) m_uiRecvOctOffsetRound1, comm);




        std::swap(m_uiAllElements, m_uiEmbeddedOctree);
        m_uiEmbeddedOctree.clear();
        m_uiAllElements.insert(m_uiAllElements.end(), m_uiGhostOctants.begin(), m_uiGhostOctants.end());

        SFC::seqSort::SFC_treeSort(&(*(m_uiAllElements.begin())), m_uiAllElements.size(), tmpNodes, tmpNodes, tmpNodes, m_uiMaxDepth,m_uiMaxDepth, rootNode, 0, 1, TS_REMOVE_DUPLICATES);
        std::swap(m_uiAllElements, tmpNodes);
        tmpNodes.clear();

        assert(seq::test::isUniqueAndSorted(m_uiAllElements));

        m_uiGhostOctants.clear();

        m_uiElementPreGhostBegin = 0;
        m_uiElementPreGhostEnd = 0;
        m_uiElementLocalBegin = 0;
        m_uiElementLocalEnd = 0;
        m_uiElementPostGhostBegin = 0;
        m_uiElementPostGhostEnd = 0;

        for (unsigned int k = 0; k < m_uiAllElements.size(); k++) {
            if (m_uiAllElements[k] == localMinMaxElement[0]) {
                m_uiElementLocalBegin = k;
                break;
            }
        }

        for (unsigned int k = (m_uiAllElements.size() - 1); k > 0; k--) {
            if (m_uiAllElements[k] == localMinMaxElement[1])
                m_uiElementLocalEnd = k + 1;
        }

        m_uiElementPreGhostEnd = m_uiElementLocalBegin;
        m_uiElementPostGhostBegin = m_uiElementLocalEnd;
        m_uiElementPostGhostEnd = m_uiAllElements.size();

        m_uiNumLocalElements=(m_uiElementLocalEnd-m_uiElementLocalBegin);
        m_uiNumPreGhostElements=(m_uiElementPreGhostEnd-m_uiElementPreGhostBegin);
        m_uiNumPostGhostElements=(m_uiElementPostGhostEnd-m_uiElementPostGhostBegin);

        m_uiNumTotalElements=m_uiNumPreGhostElements+m_uiNumLocalElements+m_uiNumPostGhostElements;

        const unsigned int FACE_EXC_ELE_BEGIN=m_uiElementLocalBegin;
        const unsigned int FACE_EXC_ELE_END=m_uiElementLocalEnd;



        // E2E Mapping for the Round 1 Ghost Exchange. Later we will perform another round of ghost exchange, (for the ghost elements that are hanging ) and build the correct complete E2E mapping.

        SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeys.begin())), &(*(m_uiAllElements.begin())), 0, m_uiKeys.size(), 0,m_uiAllElements.size(), m_uiMaxDepth, m_uiMaxDepth, 0);


        std::vector<unsigned int> *ownerList;
        std::vector<unsigned int> *stencilIndexDirection;
        unsigned int result;
        unsigned int direction;

        m_uiKeysPtr=&(*(m_uiKeys.begin()));
        m_uiE2EMapping.resize(m_uiAllElements.size()*m_uiNumDirections,LOOK_UP_TABLE_DEFAULT);

        /**Neighbour list is made from
        *  first x axis. (left to right)
        *  second y axis. (down to up)
        *  third z axis. (back to front)**/

        for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();
            if(!(OCT_FOUND & m_uiKeysPtr[k].getFlag())) continue;// Note that some keys might not be found due to absence of diagonal keys at this stage. (but all the keys should be found after diagonal neighbour exchange. )
            stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
            result = m_uiKeysPtr[k].getSearchResult();
            if(ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiKeys[k]) || m_uiAllElements[result]==m_uiKeys[k]); // To check the result found in the treeSearch is correct or not.

            for (unsigned int w = 0; w < ownerList->size(); w++) {
                direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                m_uiE2EMapping[(((*ownerList)[w]) + m_uiElementLocalBegin)*m_uiNumDirections+direction] = result;
            }

        }


        m_uiGhostKeys.clear();
        generateGhostElementSearchKeys();


        SFC::seqSearch::SFC_treeSearch(&(*(m_uiGhostKeys.begin())),&(*(m_uiAllElements.begin())),0,m_uiGhostKeys.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


        m_uiKeysPtr=&(*(m_uiGhostKeys.begin()));
        // Note : Since the ghost elements are not complete it is not required to find all the keys in the ghost.
        for (unsigned int k = 0; k < m_uiGhostKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();

            if((OCT_FOUND & m_uiKeysPtr[k].getFlag())) {

                stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
                result = m_uiKeysPtr[k].getSearchResult();
                if (ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiGhostKeys[k]) || m_uiAllElements[result] ==  m_uiGhostKeys[k]); // To check the result found in the treeSearch is correct or not.

                for (unsigned int w = 0; w < ownerList->size(); w++) {
                    direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                    m_uiE2EMapping[(((*ownerList)[w])) * m_uiNumDirections + direction] = result;
                }

            }

        }


        //===========================================================================================================================================================================================================================================

        //5- Compute missing face keys ==============================================================================================================================================================================================================
        // computing the owner ship of the ghost elements.
        std::vector<unsigned int> elementOwner;
        computeElementOwnerRanks(elementOwner);

        unsigned int lookUp;
        unsigned int gOwner;
        std::vector<unsigned int > * missedSendID=new std::vector<unsigned int >[npes]; // store the missed id after sending the ghost keys to it self.
        // for local elements.
        for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
        {
            for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
            {
                lookUp=m_uiE2EMapping[ele*m_uiNumDirections+dir];
                if(lookUp!=LOOK_UP_TABLE_DEFAULT && ((lookUp<m_uiElementLocalBegin) || (lookUp>=m_uiElementLocalEnd)))
                {
                    gOwner=elementOwner[lookUp];
                    assert(gOwner!=rank);
                    assert(gOwner<npes && gOwner>=0);
                    if(!std::binary_search(scatterMapSend_R1[gOwner].begin(),scatterMapSend_R1[gOwner].end(),(ele-m_uiElementLocalBegin)))
                        missedSendID[gOwner].push_back(ele-m_uiElementLocalBegin);

                }
            }
        }

        // for pre ghost
        for(unsigned int ele=m_uiElementPreGhostBegin;ele<m_uiElementPreGhostEnd;ele++)
        {
            for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
            {
                lookUp=m_uiE2EMapping[ele*m_uiNumDirections+dir];
                if(lookUp>=m_uiElementLocalBegin && lookUp<m_uiElementLocalEnd)
                {
                    gOwner=elementOwner[ele];
                    assert(gOwner!=rank);
                    assert(gOwner<npes && gOwner>=0);
                    if(!std::binary_search(scatterMapSend_R1[gOwner].begin(),scatterMapSend_R1[gOwner].end(),(lookUp-m_uiElementLocalBegin)))
                        missedSendID[gOwner].push_back(lookUp-m_uiElementLocalBegin);
                }

            }
        }

        // for post ghost
        for(unsigned int ele=m_uiElementPostGhostBegin;ele<m_uiElementPostGhostEnd;ele++)
        {
            for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
            {
                lookUp=m_uiE2EMapping[ele*m_uiNumDirections+dir];
                if(lookUp>=m_uiElementLocalBegin && lookUp<m_uiElementLocalEnd)
                {
                    gOwner=elementOwner[ele];
                    assert(gOwner!=rank);
                    assert(gOwner<npes && gOwner>=0);
                    if(!std::binary_search(scatterMapSend_R1[gOwner].begin(),scatterMapSend_R1[gOwner].end(),(lookUp-m_uiElementLocalBegin)))
                        missedSendID[gOwner].push_back(lookUp-m_uiElementLocalBegin);
                }

            }
        }

        m_uiSendBufferElement.clear();
        for(unsigned int p=0;p<npes;p++)
        {
            if(missedSendID[p].size()!=0)
            {
                std::sort(missedSendID[p].begin(),missedSendID[p].end());
                missedSendID[p].erase(std::unique(missedSendID[p].begin(),missedSendID[p].end()),missedSendID[p].end());
            }

            m_uiSendOctCountRound1[p]=missedSendID[p].size();
            for(unsigned int i=0;i<missedSendID[p].size();i++)
                m_uiSendBufferElement.push_back(m_uiAllElements[missedSendID[p][i]+m_uiElementLocalBegin]);

            scatterMapSend_R1[p].insert(scatterMapSend_R1[p].end(),missedSendID[p].begin(),missedSendID[p].end());
            std::sort(scatterMapSend_R1[p].begin(),scatterMapSend_R1[p].end());
            missedSendID[p].clear();
        }


        elementOwner.clear();

        par::Mpi_Alltoall(m_uiSendOctCountRound1,m_uiRecvOctCountRound1,1,comm);

        m_uiSendOctOffsetRound1[0]=0;
        m_uiRecvOctOffsetRound1[0]=0;

        omp_par::scan(m_uiSendOctCountRound1,m_uiSendOctOffsetRound1,npes);
        omp_par::scan(m_uiRecvOctCountRound1,m_uiRecvOctOffsetRound1,npes);


        m_uiGhostOctants.clear();
        m_uiGhostOctants.resize(m_uiRecvOctOffsetRound1[npes - 1] + m_uiRecvOctCountRound1[npes - 1]);

        // exchange missed keys.
        par::Mpi_Alltoallv(&(*(m_uiSendBufferElement.begin())), (int *) m_uiSendOctCountRound1, (int *) m_uiSendOctOffsetRound1,
                           &(*(m_uiGhostOctants.begin())), (int *) m_uiRecvOctCountRound1, (int *) m_uiRecvOctOffsetRound1, comm);


        m_uiAllElements.insert(m_uiAllElements.end(), m_uiGhostOctants.begin(), m_uiGhostOctants.end());

        SFC::seqSort::SFC_treeSort(&(*(m_uiAllElements.begin())), m_uiAllElements.size(), tmpNodes, tmpNodes, tmpNodes, m_uiMaxDepth,m_uiMaxDepth, rootNode, 0, 1, TS_REMOVE_DUPLICATES);
        std::swap(m_uiAllElements, tmpNodes);
        tmpNodes.clear();

        assert(seq::test::isUniqueAndSorted(m_uiAllElements));

        m_uiGhostOctants.clear();

        m_uiElementPreGhostBegin = 0;
        m_uiElementPreGhostEnd = 0;
        m_uiElementLocalBegin = 0;
        m_uiElementLocalEnd = 0;
        m_uiElementPostGhostBegin = 0;
        m_uiElementPostGhostEnd = 0;

        for (unsigned int k = 0; k < m_uiAllElements.size(); k++) {
            if (m_uiAllElements[k] == localMinMaxElement[0]) {
                m_uiElementLocalBegin = k;
                break;
            }
        }

        for (unsigned int k = (m_uiAllElements.size() - 1); k > 0; k--) {
            if (m_uiAllElements[k] == localMinMaxElement[1])
                m_uiElementLocalEnd = k + 1;
        }

        m_uiElementPreGhostEnd = m_uiElementLocalBegin;
        m_uiElementPostGhostBegin = m_uiElementLocalEnd;
        m_uiElementPostGhostEnd = m_uiAllElements.size();

        m_uiNumLocalElements=(m_uiElementLocalEnd-m_uiElementLocalBegin);
        m_uiNumPreGhostElements=(m_uiElementPreGhostEnd-m_uiElementPreGhostBegin);
        m_uiNumPostGhostElements=(m_uiElementPostGhostEnd-m_uiElementPostGhostBegin);

        m_uiNumTotalElements=m_uiNumPreGhostElements+m_uiNumLocalElements+m_uiNumPostGhostElements;

        const unsigned int FACE_MISNG_EXC_ELE_BEGIN=m_uiElementLocalBegin;
        const unsigned int FACE_MISNG_EXC_ELE_END=m_uiElementLocalEnd;


#ifdef DEBUG_MESH_GENERATION

        for(unsigned int ele=m_uiElementPreGhostBegin;ele<m_uiElementPreGhostEnd;ele++)
            m_uiGhostOctants.push_back(m_uiAllElements[ele]);

        for(unsigned int ele=m_uiElementPostGhostBegin;ele<m_uiElementPostGhostEnd;ele++)
            m_uiGhostOctants.push_back(m_uiAllElements[ele]);

        treeNodesTovtk(m_uiGhostOctants,rank,"ghostR1");
        m_uiGhostOctants.clear();

#endif


        // E2E Mapping for the Round face-1 & face-2 face ghost exchange.

        SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeys.begin())), &(*(m_uiAllElements.begin())), 0, m_uiKeys.size(), 0,m_uiAllElements.size(), m_uiMaxDepth, m_uiMaxDepth, 0);

        m_uiKeysPtr=&(*(m_uiKeys.begin()));
        m_uiE2EMapping.resize(m_uiAllElements.size()*m_uiNumDirections,LOOK_UP_TABLE_DEFAULT);

        /**Neighbour list is made from
        *  first x axis. (left to right)
        *  second y axis. (down to up)
        *  third z axis. (back to front)**/

        for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();
            if(ownerList->size() && (!((OCT_FOUND & m_uiKeysPtr[k].getFlag())))) {
                std::cout<<"rank: "<<m_uiActiveRank<<"[E2E Error]: Local  face key missing after R1 face-1 & face-2 ghost exchange: "<<m_uiKeysPtr[k]<<std::endl;
                exit(0); // no point in continuing if this fails E2N fails for sure :)
                //std::cout<<"rank: "<<m_uiActiveRank<<" key: "<<m_uiKeysPtr[k]<<" missing "<<" status: "<<(splitterElements[4]<=m_uiKeysPtr[k] && m_uiKeysPtr[k]<=splitterElements[5])<<" k: "<<k<<" sbegin: "<<splitterElements[44].getSearchResult()<<" sEnd: "<<splitterElements[45].getSearchResult()<<std::endl;
                //missingKeys.push_back(m_uiKeysPtr[k]);
                //continue;
            }

            stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
            result = m_uiKeysPtr[k].getSearchResult();
            if(ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiKeys[k]) || m_uiAllElements[result]==m_uiKeys[k]); // To check the result found in the treeSearch is correct or not.

            for (unsigned int w = 0; w < ownerList->size(); w++) {
                direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                m_uiE2EMapping[(((*ownerList)[w]) + m_uiElementLocalBegin)*m_uiNumDirections+direction] = result;
            }

        }

        m_uiGhostKeys.clear();
        generateGhostElementSearchKeys();

        SFC::seqSearch::SFC_treeSearch(&(*(m_uiGhostKeys.begin())),&(*(m_uiAllElements.begin())),0,m_uiGhostKeys.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


        m_uiKeysPtr=&(*(m_uiGhostKeys.begin()));
        // Note : Since the ghost elements are not complete it is not required to find all the keys in the ghost.
        for (unsigned int k = 0; k < m_uiGhostKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();

            if((OCT_FOUND & m_uiKeysPtr[k].getFlag())) {
                stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
                result = m_uiKeysPtr[k].getSearchResult();
                if (ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiGhostKeys[k]) || m_uiAllElements[result] ==  m_uiGhostKeys[k]); // To check the result found in the treeSearch is correct or not.

                for (unsigned int w = 0; w < ownerList->size(); w++) {
                    direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                    m_uiE2EMapping[(((*ownerList)[w])) * m_uiNumDirections + direction] = result;
                    // Note : Following is done to enforce that the local elements that points to ghost elements has the inverse mapping.
                }

            }

        }


    //===========================================================================================================================================================================================================================================


    // 6- Perform exchange for edge and vertex neighbors. =======================================================================================================================================================================================


        m_uiKeysDiag.clear();
        generateBdyElementDiagonalSearchKeys();

#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(m_uiKeysDiag,rank,"m_uiKeyDiag");
#endif



        tmpKeys.clear();
        SFC::seqSort::SFC_treeSort(&(*(m_uiKeysDiag.begin())),m_uiKeysDiag.size(),tmpKeys,tmpKeys,tmpKeys,m_uiMaxDepth,m_uiMaxDepth,rootKey,ROOT_ROTATION,1,TS_SORT_ONLY);
        tmpKeys.clear();

        for(unsigned int i=0;i<2*npes;i++)
            splitterElements[i]=ot::Key(m_uiLocalSplitterElements[i]);

        assert(seq::test::isUniqueAndSorted(splitterElements));

        SFC::seqSearch::SFC_treeSearch(&(*(splitterElements.begin())),&(*(m_uiKeysDiag.begin())),0,splitterElements.size(),0,m_uiKeysDiag.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


        m_uiKeysPtr=&(*(m_uiKeysDiag.begin()));
        m_uiGhostElementIDsToBeSent.clear();
        m_uiSendBufferElement.clear();

        m_uiSendKeyDiagCount = new unsigned int[npes];
        m_uiRecvKeyDiagCount = new unsigned int[npes];

        m_uiSendKeyDiagOffset = new unsigned int[npes];
        m_uiRecvKeyDiagOffset = new unsigned int[npes];


        m_uiSendOctCountRound1Diag = new unsigned int[npes];
        m_uiRecvOctCountRound1Diag = new unsigned int[npes];

        m_uiSendOctOffsetRound1Diag = new unsigned int[npes];
        m_uiRecvOctOffsetRound1Diag = new unsigned int[npes];


        for(unsigned int p=0;p<npes;p++)
            m_uiSendKeyDiagCount[p]=0;


        m_uiSendBufferElement.clear();

        for(unsigned int p=0;p<npes;p++)
        {
            assert((splitterElements[2*p].getFlag() & OCT_FOUND));
            assert((splitterElements[2*p+1].getFlag() & OCT_FOUND));;
            assert(m_uiKeysDiag[splitterElements[2*p].getSearchResult()]==splitterElements[2*p]);
            assert(m_uiKeysDiag[splitterElements[2*p+1].getSearchResult()]==splitterElements[2*p+1]);

            sBegin=splitterElements[2*p].getSearchResult();
            //sEnd=splitterElements[2*p+1].getSearchResult()+1;
            (p<(m_uiActiveNpes-1)) ? sEnd=splitterElements[2*p+2].getSearchResult()+1: sEnd=m_uiKeysDiag.size();
            if(p!=m_uiActiveRank)
            for(unsigned int i=sBegin;i<sEnd;i++)
            {
               if(m_uiKeysDiag[i].getOwnerList()->size())
               {
                   m_uiSendBufferElement.push_back(m_uiKeysDiag[i]);
                   m_uiSendKeyDiagCount[p]++;
               }


            }

        }


        par::Mpi_Alltoall(m_uiSendKeyDiagCount,m_uiRecvKeyDiagCount,1,comm);

        m_uiSendKeyDiagOffset[0]=0;
        m_uiRecvKeyDiagOffset[0]=0;

        omp_par::scan(m_uiSendKeyDiagCount,m_uiSendKeyDiagOffset,npes);
        omp_par::scan(m_uiRecvKeyDiagCount,m_uiRecvKeyDiagOffset,npes);


        std::vector<ot::TreeNode> recvDiagKeyOct;
        recvDiagKeyOct.resize(m_uiRecvKeyDiagOffset[npes-1]+m_uiRecvKeyDiagCount[npes-1]);

        par::Mpi_Alltoallv(&(*(m_uiSendBufferElement.begin())), (int *) m_uiSendKeyDiagCount, (int *) m_uiSendKeyDiagOffset,
                           &(*(recvDiagKeyOct.begin())), (int *) m_uiRecvKeyDiagCount, (int *) m_uiRecvKeyDiagOffset, comm);


#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(recvDiagKeyOct,rank,"recvDiagKeyOct");
#endif

        std::vector<ot::Key> recvDiagKey_keys;
        std::vector<ot::Key>::iterator itKey;
        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=m_uiRecvKeyDiagOffset[p];e<(m_uiRecvKeyDiagOffset[p]+m_uiRecvKeyDiagCount[p]);e++)
            {
                itKey=recvDiagKey_keys.emplace(recvDiagKey_keys.end(),ot::Key(recvDiagKeyOct[e]));
                itKey->addOwner(p);
            }
        }


        SFC::seqSearch::SFC_treeSearch(&(*(recvDiagKey_keys.begin())),&(*(m_uiAllElements.begin())),0,recvDiagKey_keys.size(),m_uiElementLocalBegin,m_uiElementLocalEnd,m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);

        for(unsigned int p=0;p<npes;p++)
        {
            m_uiSendOctCountRound1Diag[p]=0;
            missedSendID[p].clear();
        }

        // note that not all the recv keys needs to be found since recvkey send is overlapped among the processors.
        for(unsigned int e=0;e<recvDiagKey_keys.size();e++)
        {
            if(!(recvDiagKey_keys[e].getFlag() & OCT_FOUND))
                continue;

            ownerList=recvDiagKey_keys[e].getOwnerList();
            result=recvDiagKey_keys[e].getSearchResult();
            for (unsigned int w = 0; w < ownerList->size(); w++)
            {
                missedSendID[(*ownerList)[w]].push_back(result-m_uiElementLocalBegin);
            }

        }

        m_uiSendBufferElement.clear();
        for(unsigned int p=0;p<npes;p++)
        {
            std::sort(missedSendID[p].begin(),missedSendID[p].end());
            missedSendID[p].erase(std::unique(missedSendID[p].begin(),missedSendID[p].end()),missedSendID[p].end());
            m_uiSendOctCountRound1Diag[p]=missedSendID[p].size();

            for(unsigned int e=0;e<missedSendID[p].size();e++)
            {
                m_uiSendBufferElement.push_back(m_uiAllElements[missedSendID[p][e]+m_uiElementLocalBegin]);
            }

            scatterMapSend_R1[p].insert(scatterMapSend_R1[p].end(),missedSendID[p].begin(),missedSendID[p].end());
            missedSendID[p].clear();


        }

        delete [] missedSendID;

        par::Mpi_Alltoall(m_uiSendOctCountRound1Diag,m_uiRecvOctCountRound1Diag,1,comm);

        m_uiSendOctOffsetRound1Diag[0]=0;
        m_uiRecvOctOffsetRound1Diag[0]=0;

        omp_par::scan(m_uiSendOctCountRound1Diag,m_uiSendOctOffsetRound1Diag,npes);
        omp_par::scan(m_uiRecvOctCountRound1Diag,m_uiRecvOctOffsetRound1Diag,npes);

        m_uiGhostOctants.clear();

        m_uiGhostOctants.resize(m_uiRecvOctOffsetRound1Diag[npes - 1] + m_uiRecvOctCountRound1Diag[npes - 1]);
        par::Mpi_Alltoallv(&(*(m_uiSendBufferElement.begin())), (int *) m_uiSendOctCountRound1Diag, (int *) m_uiSendOctOffsetRound1Diag,
                           &(*(m_uiGhostOctants.begin())), (int *) m_uiRecvOctCountRound1Diag, (int *) m_uiRecvOctOffsetRound1Diag, comm);




#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(m_uiGhostOctants,rank,"ghostR1Diag");
        std::vector<ot::TreeNode> missingKeys;
#endif

        m_uiAllElements.insert(m_uiAllElements.end(),m_uiGhostOctants.begin(),m_uiGhostOctants.end());
        m_uiGhostOctants.clear();


        SFC::seqSort::SFC_treeSort(&(*(m_uiAllElements.begin())), m_uiAllElements.size(), tmpNodes, tmpNodes, tmpNodes, m_uiMaxDepth,m_uiMaxDepth, rootNode, 0, 1, TS_REMOVE_DUPLICATES);
        std::swap(m_uiAllElements, tmpNodes);
        tmpNodes.clear();

        assert(seq::test::isUniqueAndSorted(m_uiAllElements));

        m_uiElementPreGhostBegin = 0;
        m_uiElementPreGhostEnd = 0;
        m_uiElementLocalBegin = 0;
        m_uiElementLocalEnd = 0;
        m_uiElementPostGhostBegin = 0;
        m_uiElementPostGhostEnd = 0;

        for (unsigned int k = 0; k < m_uiAllElements.size(); k++) {
            if (m_uiAllElements[k] == localMinMaxElement[0]) {
                m_uiElementLocalBegin = k;
                break;
            }
        }

        for (unsigned int k = (m_uiAllElements.size() - 1); k > 0; k--) {
            if (m_uiAllElements[k] == localMinMaxElement[1])
                m_uiElementLocalEnd = k + 1;
        }

        m_uiElementPreGhostEnd = m_uiElementLocalBegin;
        m_uiElementPostGhostBegin = m_uiElementLocalEnd;
        m_uiElementPostGhostEnd = m_uiAllElements.size();

        m_uiNumLocalElements=(m_uiElementLocalEnd-m_uiElementLocalBegin);
        m_uiNumPreGhostElements=(m_uiElementPreGhostEnd-m_uiElementPreGhostBegin);
        m_uiNumPostGhostElements=(m_uiElementPostGhostEnd-m_uiElementPostGhostBegin);

        m_uiNumTotalElements=m_uiNumPreGhostElements+m_uiNumLocalElements+m_uiNumPostGhostElements;


        const unsigned int EDGE_VERTEX_EXC_ELE_BEGIN=m_uiElementLocalBegin;
        const unsigned int EDGE_VERTEX_EXC_ELE_END=m_uiElementLocalEnd;


        // E2E Mapping for the Round 1 Ghost Exchange. Later we will perform another round of ghost exchange, (for the ghost elements that are hanging ) and build the correct complete E2E mapping.

        SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeys.begin())), &(*(m_uiAllElements.begin())), 0, m_uiKeys.size(), 0,m_uiAllElements.size(), m_uiMaxDepth, m_uiMaxDepth, 0);


        m_uiKeysPtr=&(*(m_uiKeys.begin()));
        m_uiE2EMapping.resize(m_uiAllElements.size()*m_uiNumDirections,LOOK_UP_TABLE_DEFAULT);

        /**Neighbour list is made from
        *  first x axis. (left to right)
        *  second y axis. (down to up)
        *  third z axis. (back to front)**/

        for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();

            if(ownerList->size() && (!((OCT_FOUND & m_uiKeysPtr[k].getFlag())))) {
                std::cout<<"rank: "<<m_uiActiveRank<<"[E2E Error]: Local  face key missing after R1 (face edge vertex) ghost exchange: "<<m_uiKeysPtr[k]<<std::endl;
                exit(0); // no point in continuing if this fails E2N fails for sure :)
                //std::cout<<"rank: "<<m_uiActiveRank<<" key: "<<m_uiKeysPtr[k]<<" missing "<<" status: "<<(splitterElements[4]<=m_uiKeysPtr[k] && m_uiKeysPtr[k]<=splitterElements[5])<<" k: "<<k<<" sbegin: "<<splitterElements[44].getSearchResult()<<" sEnd: "<<splitterElements[45].getSearchResult()<<std::endl;
                //missingKeys.push_back(m_uiKeysPtr[k]);
                //continue;
            }

            if(ownerList->size()) assert((OCT_FOUND & m_uiKeysPtr[k].getFlag()));// Note that all the keys should be found locally due to the fact that we have exchanged ghost elements.
            stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
            result = m_uiKeysPtr[k].getSearchResult();
            if(ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiKeys[k]) || m_uiAllElements[result]==m_uiKeys[k]); // To check the result found in the treeSearch is correct or not.
            for (unsigned int w = 0; w < ownerList->size(); w++) {
                direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                m_uiE2EMapping[(((*ownerList)[w]) + m_uiElementLocalBegin)*m_uiNumDirections+direction] = result;

            }

        }

        m_uiGhostKeys.clear();
        generateGhostElementSearchKeys();

        SFC::seqSearch::SFC_treeSearch(&(*(m_uiGhostKeys.begin())),&(*(m_uiAllElements.begin())),0,m_uiGhostKeys.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


        m_uiKeysPtr=&(*(m_uiGhostKeys.begin()));
        // Note : Since the ghost elements are not complete it is not required to find all the keys in the ghost.
        for (unsigned int k = 0; k < m_uiGhostKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();

            if((OCT_FOUND & m_uiKeysPtr[k].getFlag())) {
                stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
                result = m_uiKeysPtr[k].getSearchResult();
                if (ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiGhostKeys[k]) || m_uiAllElements[result] ==  m_uiGhostKeys[k]); // To check the result found in the treeSearch is correct or not.

                for (unsigned int w = 0; w < ownerList->size(); w++) {
                    direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                    m_uiE2EMapping[(((*ownerList)[w])) * m_uiNumDirections + direction] = result;
                    // Note : Following is done to enforce that the local elements that points to ghost elements has the inverse mapping.
                }

            }

        }



        SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeysDiag.begin())),&(*(m_uiAllElements.begin())),0,m_uiKeysDiag.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);

        m_uiKeysPtr=&(*(m_uiKeysDiag.begin()));
        for (unsigned int k = 0; k < m_uiKeysDiag.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();

            if(ownerList->size() && (!((OCT_FOUND & m_uiKeysPtr[k].getFlag())))) {
                std::cout<<"rank: "<<m_uiActiveRank<<"[E2E Error]: Local edge or vertex key missing after R1 (face edge vertex) ghost exchange: "<<m_uiKeysPtr[k]<<std::endl;
                //missingKeys.push_back(m_uiKeysPtr[k]);
                //continue;
                exit(0); // no point in continuing if this fails E2N fails for sure :)
            }

        }



        // merge the face neighbors with edge and vertex neighbors,
        // Note that m_uiScatterMapElementRound1 should only contain the local element IDs.
        // (pure local ID (ID's before even merging with ghost. )) these IDs indicate that my local element is an ghost element to some other processor. Since we have all the neighours for the boundary nodes we can
        // determine whether the each local element is hanging or not. If is it hanging it will be participated in the layer 2 ghost level exchange.

        m_uiScatterMapElementRound1.clear();
        for(unsigned int p=0;p<npes;p++)
        {
            std::sort(scatterMapSend_R1[p].begin(),scatterMapSend_R1[p].end());
            scatterMapSend_R1[p].erase(std::unique(scatterMapSend_R1[p].begin(),scatterMapSend_R1[p].end()),scatterMapSend_R1[p].end());
            m_uiScatterMapElementRound1.insert(m_uiScatterMapElementRound1.end(),scatterMapSend_R1[p].begin(),scatterMapSend_R1[p].end());
            m_uiSendOctCountRound1[p]=scatterMapSend_R1[p].size();
            scatterMapSend_R1[p].clear();
        }

        m_uiSendOctOffsetRound1[0]=0;
        omp_par::scan(m_uiSendOctCountRound1,m_uiSendOctOffsetRound1,npes);

        delete [] scatterMapSend_R1;



        // push the diagonal ghost layer 1 keys. this includes the face
        std::vector<ot::TreeNode> gKeys_R1;
        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPreGhostEnd;e++)
            gKeys_R1.push_back(m_uiAllElements[e]);

        for(unsigned int e=m_uiElementPostGhostBegin;e<m_uiElementPostGhostEnd;e++)
            gKeys_R1.push_back(m_uiAllElements[e]);



        //8 R2 ghost exchange =====================================================================================================================================================================================================

        m_uiSendOctCountRound2 = new unsigned int[npes];
        m_uiRecvOctCountRound2 = new unsigned int[npes];

        m_uiSendOctOffsetRound2 = new unsigned int[npes];
        m_uiRecvOctOffsetRound2 = new unsigned int[npes];

        tmpSendElementIds.clear();
        m_uiSendBufferElement.clear();
        unsigned int elementLookup;
        std::set<unsigned int > * tmpSendEleIdR2=new std::set<unsigned int>[npes];
        std::pair<std::set<unsigned int>::iterator,bool> setHintUint;


        for(unsigned int p=0;p<npes;p++)
        {
            m_uiSendOctCountRound2[p]=0;
            for(unsigned int ele=m_uiSendOctOffsetRound1[p];ele<(m_uiSendOctOffsetRound1[p]+m_uiSendOctCountRound1[p]);ele++)
            {

                for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
                {
                    elementLookup=m_uiE2EMapping[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)*m_uiNumDirections+dir];

                    if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                        setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);
                    }

                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_LEFT,OCT_DIR_DOWN,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_LEFT,OCT_DIR_UP,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_LEFT,OCT_DIR_BACK,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_LEFT,OCT_DIR_FRONT,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }


                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_RIGHT,OCT_DIR_DOWN,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_RIGHT,OCT_DIR_UP,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }


                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_RIGHT,OCT_DIR_BACK,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }


                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_RIGHT,OCT_DIR_FRONT,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }


                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_DOWN,OCT_DIR_BACK,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_DOWN,OCT_DIR_FRONT,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);
                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_UP,OCT_DIR_BACK,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }

                OCT_DIR_DIAGONAL_E2E((m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin),OCT_DIR_UP,OCT_DIR_FRONT,elementLookup);
                if((elementLookup!=LOOK_UP_TABLE_DEFAULT) && (m_uiAllElements[elementLookup].getLevel()<=m_uiAllElements[(m_uiScatterMapElementRound1[ele]+m_uiElementLocalBegin)].getLevel())) {
                    setHintUint = tmpSendEleIdR2[p].emplace(elementLookup);

                }

            }
        }

        m_uiSendBufferElement.clear();
        std::vector<unsigned int> common_data;
        for(unsigned int pi=0;pi<npes;pi++)
        {
            tmpSendElementIds.clear();
            tmpSendElementIds.insert(tmpSendEleIdR2[pi].begin(),tmpSendEleIdR2[pi].end());


             for(auto it=tmpSendElementIds.begin();it!=tmpSendElementIds.end();++it)
             {
                m_uiSendBufferElement.push_back(m_uiAllElements[*it]);
             }

            m_uiSendOctCountRound2[pi]=tmpSendElementIds.size();



        }

        delete[] tmpSendEleIdR2;


            par::Mpi_Alltoall(m_uiSendOctCountRound2, m_uiRecvOctCountRound2, 1, comm);
            m_uiSendOctOffsetRound2[0] = 0;
            m_uiRecvOctOffsetRound2[0] = 0;

            omp_par::scan(m_uiSendOctCountRound2,m_uiSendOctOffsetRound2,npes);
            omp_par::scan(m_uiRecvOctCountRound2,m_uiRecvOctOffsetRound2,npes);


            m_uiGhostOctants.clear();
            m_uiGhostOctants.resize(m_uiRecvOctOffsetRound2[npes-1]+m_uiRecvOctCountRound2[npes-1]);
            assert(m_uiSendBufferElement.size()==(m_uiSendOctOffsetRound2[npes-1]+m_uiSendOctCountRound2[npes-1]));
            par::Mpi_Alltoallv(&(*(m_uiSendBufferElement.begin())),(int *)m_uiSendOctCountRound2,(int *)m_uiSendOctOffsetRound2,&(*(m_uiGhostOctants.begin())),(int *)m_uiRecvOctCountRound2,(int *)m_uiRecvOctOffsetRound2,comm);



#ifdef DEBUG_MESH_GENERATION
            treeNodesTovtk(m_uiGhostOctants,rank,"ghostR2");
#endif


            m_uiAllElements.insert(m_uiAllElements.end(),m_uiGhostOctants.begin(),m_uiGhostOctants.end());


            tmpNodes.clear();
            SFC::seqSort::SFC_treeSort(&(*(m_uiAllElements.begin())), m_uiAllElements.size(), tmpNodes, tmpNodes, tmpNodes, m_uiMaxDepth,m_uiMaxDepth, rootNode, 0, 1, TS_REMOVE_DUPLICATES);
            std::swap(m_uiAllElements, tmpNodes);
            tmpNodes.clear();

            assert(seq::test::isUniqueAndSorted(m_uiAllElements));
            m_uiGhostOctants.clear();

            m_uiElementPreGhostBegin = 0;
            m_uiElementPreGhostEnd = 0;
            m_uiElementLocalBegin = 0;
            m_uiElementLocalEnd = 0;
            m_uiElementPostGhostBegin = 0;
            m_uiElementPostGhostEnd = 0;

            for (unsigned int k = 0; k < m_uiAllElements.size(); k++) {
                if (m_uiAllElements[k] == localMinMaxElement[0]) {
                    m_uiElementLocalBegin = k;
                    break;
                }
            }

            for (unsigned int k = (m_uiAllElements.size() - 1); k > 0; k--) {
                if (m_uiAllElements[k] == localMinMaxElement[1])
                    m_uiElementLocalEnd = k + 1;
            }

            m_uiElementPreGhostEnd = m_uiElementLocalBegin;
            m_uiElementPostGhostBegin = m_uiElementLocalEnd;
            m_uiElementPostGhostEnd = m_uiAllElements.size();

            m_uiNumLocalElements=(m_uiElementLocalEnd-m_uiElementLocalBegin);
            m_uiNumPreGhostElements=(m_uiElementPreGhostEnd-m_uiElementPreGhostBegin);
            m_uiNumPostGhostElements=(m_uiElementPostGhostEnd-m_uiElementPostGhostBegin);

            m_uiNumTotalElements=m_uiNumPreGhostElements+m_uiNumLocalElements+m_uiNumPostGhostElements;


            // E2E Mapping for the Round 2 Ghost Exchange. (Final E2E Mapping )

            SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeys.begin())), &(*(m_uiAllElements.begin())), 0, m_uiKeys.size(), 0,m_uiAllElements.size(), m_uiMaxDepth, m_uiMaxDepth, 0);


            m_uiKeysPtr=&(*(m_uiKeys.begin()));
            m_uiE2EMapping.clear();
            m_uiE2EMapping.resize(m_uiAllElements.size()*m_uiNumDirections,LOOK_UP_TABLE_DEFAULT);

            // * Neighbour list is made from
            // *  first x axis. (left to right)
            // *  second y axis. (down to up)
            // *  third z axis. (back to front)

            for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
                ownerList = m_uiKeysPtr[k].getOwnerList();

                if(ownerList->size()) assert((OCT_FOUND & m_uiKeysPtr[k].getFlag()));// Note that all the keys should be found locally due to the fact that we have exchanged ghost elements.
                stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
                result = m_uiKeysPtr[k].getSearchResult();
                if(ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiKeys[k]) || m_uiAllElements[result]==m_uiKeys[k]); // To check the result found in the treeSearch is correct or not.
                for (unsigned int w = 0; w < ownerList->size(); w++) {
                    direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                    m_uiE2EMapping[(((*ownerList)[w]) + m_uiElementLocalBegin)*m_uiNumDirections+direction] = result;
                    // Note : Following is done to enforce that the local elements that points to ghost elements has the inverse mapping.
                    if(result<m_uiElementLocalBegin || result>=m_uiElementLocalEnd)
                        m_uiE2EMapping[result*m_uiNumDirections+ (1u^direction)]=(((*ownerList)[w]) + m_uiElementLocalBegin); // Note This depends on the OCT_DIR numbering.


                }

            }

        //}



        assert(seq::test::checkE2EMapping(m_uiE2EMapping,m_uiAllElements,m_uiElementLocalBegin,m_uiElementLocalEnd,1,m_uiNumDirections));



        m_uiGhostKeys.clear();
        generateGhostElementSearchKeys();

        SFC::seqSearch::SFC_treeSearch(&(*(m_uiGhostKeys.begin())),&(*(m_uiAllElements.begin())),0,m_uiGhostKeys.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


        m_uiKeysPtr=&(*(m_uiGhostKeys.begin()));
        // Note : Since the ghost elements are not complete it is not required to find all the keys in the ghost.
        for (unsigned int k = 0; k < m_uiGhostKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();
            if((OCT_FOUND & m_uiKeysPtr[k].getFlag())) {
                stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
                result = m_uiKeysPtr[k].getSearchResult();
                if (ownerList->size()) assert(m_uiAllElements[result].isAncestor(m_uiGhostKeys[k]) || m_uiAllElements[result] ==  m_uiGhostKeys[k]); // To check the result found in the treeSearch is correct or not.

                for (unsigned int w = 0; w < ownerList->size(); w++) {
                    direction = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                    m_uiE2EMapping[(((*ownerList)[w])) * m_uiNumDirections + direction] = result;
                    // Note : Following is done to enforce that the local elements that points to ghost elements has the inverse mapping.
                }

            }

        }




//===========================================================================================================================================================================================================================================

#ifdef DEBUG_MESH_GENERATION
        /* if(!m_uiActiveRank)
        {
            for(unsigned int ele=0;ele<m_uiE2HangingFaces.size();ele++)
            {
                printf("ele: %d value: %d \n",ele,m_uiE2HangingFaces[ele]);
            }
        }*/
        treeNodesTovtk(m_uiAllElements,rank,"m_uiAllElements");
        //std::cout<<"rank: "<<rank<<"pre begin: "<<m_uiElementPreGhostBegin<<" pre end: "<<m_uiElementPreGhostEnd<<" local begin: "<<m_uiElementLocalBegin<<" local end: "<<m_uiElementLocalEnd<<" post begin: "<<m_uiElementPostGhostBegin<<" post end: "<<m_uiElementPostGhostEnd<<std::endl;
#endif

        // Note: Note that m_uiAlllNodes need not to be sorted and contains duplicates globally.


#ifdef DEBUG_MESH_GENERATION
        std::vector<ot::TreeNode> missedKeys;
            for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
                if (!(m_uiAllElements[result].isAncestor(m_uiKeys[k]) || m_uiAllElements[result]==m_uiKeys[k])) {
                    /*std::cout << "rank: " << rank << " key : " << m_uiKeys[k] << " not found!" << " Search index count: "
                              << searchResults.size() << std::endl;*/
                    missedKeys.push_back(m_uiKeys[k]);
                }
            }

        treeNodesTovtk(missedKeys,rank,"missedKeys");
#endif



#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(notFoundKeys,rank,"notFoundKeys");
        treeNodesTovtk(m_uiGhostKeys,rank,"m_uiGhostKeys");
        treeNodesTovtk(foundResult,rank,"foundResult");
#endif


        // to identify the true ghost (level 1) elements.
        // Defn: True level 1 ghost element defined as any ghost element who shares a face (now edge and vertex) with local element.
        m_uiGhostElementRound1Index.clear();
        unsigned int r1_count=0;
        m_uiGhostElementRound1Index.resize(gKeys_R1.size());
        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPreGhostEnd;e++)
        {
            if(m_uiAllElements[e]==gKeys_R1[r1_count])
            {
                m_uiGhostElementRound1Index[r1_count]=e;
                r1_count++;
            }

            if(r1_count==gKeys_R1.size()) break;
        }

        if(r1_count<gKeys_R1.size())
        {

            for(unsigned int e=m_uiElementPostGhostBegin;e<m_uiElementPostGhostEnd;e++)
            {
                if(m_uiAllElements[e]==gKeys_R1[r1_count])
                {
                    m_uiGhostElementRound1Index[r1_count]=e;
                    r1_count++;
                }

                if(r1_count==gKeys_R1.size()) break;
            }


        }

        //std::cout<<" rank: "<<m_uiActiveRank<<" m_uiGR1 Size: "<<m_uiGhostElementRound1Index.size()<<std::endl;


        // BELOW CODE IS RISKY. YOU SHOULD NOT ADD R2 GHOST AS R1 GHOST WHEN GENERATING SM. THIS IS NOT CORRECT. IF WE HAVE DONE THE R1 GHOST EXCHANGE EXCHANGE CORRECTLY WE SHOULD NOT NEED THE BELOW CODE.
        // KEEPING THIS IF WE NEEDED A QUICK FIX. AGAIN DO NOT ENABLE THE BELOW CODE. !!!!
        /*
        //m_uiGhostElementRound1Index.clear();
        SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeys.begin())),&(*(m_uiAllElements.begin())),0,m_uiKeys.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);
        SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeysDiag.begin())),&(*(m_uiAllElements.begin())),0,m_uiKeysDiag.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


        m_uiKeysPtr=&(*(m_uiKeys.begin()));
        for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();

            if(ownerList->size() && (!((OCT_FOUND & m_uiKeysPtr[k].getFlag())))) {
                std::cout<<"rank: "<<m_uiActiveRank<<"[E2E Error]: Local key missing after R1 & R1 Diag & R2  ghost exchange: "<<m_uiKeysPtr[k]<<std::endl;
                exit(0); // no point in continuing if this fails E2N fails for sure :)
            }

            if(ownerList->size())
            {
                result=m_uiKeysPtr[k].getSearchResult();
                if((result<m_uiElementLocalBegin) || (result>=m_uiElementLocalEnd))
                    m_uiGhostElementRound1Index.push_back(result);
            }


        }

        m_uiKeysPtr=&(*(m_uiKeysDiag.begin()));
        for (unsigned int k = 0; k < m_uiKeysDiag.size(); k++) {
            ownerList = m_uiKeysPtr[k].getOwnerList();

            if(ownerList->size() && (!((OCT_FOUND & m_uiKeysPtr[k].getFlag())))) {
                std::cout<<"rank: "<<m_uiActiveRank<<"[E2E Error]: Local Diag and corner key missing after R1 & R1 Diag & R2  ghost exchange: "<<m_uiKeysPtr[k]<<std::endl;
                exit(0); // no point in continuing if this fails E2N fails for sure :)
            }

            if(ownerList->size())
            {
                result=m_uiKeysPtr[k].getSearchResult();
                if((result<m_uiElementLocalBegin) || (result>=m_uiElementLocalEnd))
                    m_uiGhostElementRound1Index.push_back(result);
            }

        }
        std::sort(m_uiGhostElementRound1Index.begin(),m_uiGhostElementRound1Index.end());
        m_uiGhostElementRound1Index.erase(std::unique(m_uiGhostElementRound1Index.begin(),m_uiGhostElementRound1Index.end()),m_uiGhostElementRound1Index.end());*/


#ifdef DEBUG_MESH_GENERATION
        std::vector<ot::TreeNode> r1GhostElements;
        std::vector<ot::TreeNode> localElements;
        std::vector<ot::TreeNode> ghostElements;

        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPreGhostEnd;e++)
            ghostElements.push_back(m_uiAllElements[e]);

        for(unsigned int e=m_uiElementPostGhostBegin;e<m_uiElementPostGhostEnd;e++)
            ghostElements.push_back(m_uiAllElements[e]);

        for(unsigned int e=m_uiElementLocalBegin;e<m_uiElementLocalEnd;e++)
            localElements.push_back(m_uiAllElements[e]);

        for(unsigned int e=0;e<m_uiGhostElementRound1Index.size();e++)
            r1GhostElements.push_back(m_uiAllElements[m_uiGhostElementRound1Index[e]]);

        treeNodesTovtk(r1GhostElements,rank,"r1GhostElements");
        treeNodesTovtk(localElements,rank,"localElements");
        treeNodesTovtk(ghostElements,rank,"ghostElements");
#endif
        //std::cout<<" rank: "<<m_uiActiveRank<<" m_uiGR2 Size: "<<m_uiGhostElementRound1Index.size()<<std::endl;

        if (!rank)  std::cout << "E2E mapping Ended" << std::endl;


    }

    void Mesh::buildE2EMap(std::vector<ot::TreeNode> &in)
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

       /* addBoundaryNodesType1(in, m_uiPositiveBoundaryOctants, m_uiDim, m_uiMaxDepth);
        m_uiMaxDepth = m_uiMaxDepth + 1;*/

        std::swap(m_uiEmbeddedOctree, in);
       /*  m_uiEmbeddedOctree.insert(m_uiEmbeddedOctree.end(), m_uiPositiveBoundaryOctants.begin(),m_uiPositiveBoundaryOctants.end());*/

        std::vector<ot::TreeNode> tmpNodes;
        ot::TreeNode rootNode(0, 0, 0, 0, m_uiDim, m_uiMaxDepth);
       // clear the in and positive boundary octants;
       // m_uiPositiveBoundaryOctants.clear();
        in.clear();

        //std::cout<<"rank: "<<m_uiActiveRank<<" par::sort begin: "<<m_uiMaxDepth<<" "<<m_uiDim<<" embedded octreeSize: "<<m_uiEmbeddedOctree.size()<<std::endl;
        SFC::seqSort::SFC_treeSort(&(*(m_uiEmbeddedOctree.begin())), m_uiEmbeddedOctree.size(), tmpNodes, tmpNodes, tmpNodes,m_uiMaxDepth, m_uiMaxDepth, rootNode, 0, 1, TS_REMOVE_DUPLICATES);

        std::swap(tmpNodes, m_uiEmbeddedOctree);
        tmpNodes.clear();

        generateSearchKeys(); // generates keys for sequential case.

        std::swap(m_uiAllElements,m_uiEmbeddedOctree);
        m_uiEmbeddedOctree.clear();

        //update the loop counters.
        m_uiElementPreGhostBegin=0;
        m_uiElementPreGhostEnd=0;
        m_uiElementLocalBegin=0;
        m_uiElementLocalEnd=m_uiAllElements.size();
        m_uiElementPostGhostBegin=m_uiElementLocalEnd;
        m_uiElementPostGhostEnd=m_uiElementLocalEnd;

        m_uiNumLocalElements=(m_uiElementLocalEnd-m_uiElementLocalBegin);
        m_uiNumPreGhostElements=(m_uiElementPreGhostEnd-m_uiElementPreGhostBegin);
        m_uiNumPostGhostElements=(m_uiElementPostGhostEnd-m_uiElementPostGhostBegin);

        m_uiNumTotalElements=m_uiNumPreGhostElements+m_uiNumLocalElements+m_uiNumPostGhostElements;

        m_uiLocalSplitterElements.resize(2);
        m_uiLocalSplitterElements[0]=m_uiAllElements.front();
        m_uiLocalSplitterElements[1]=m_uiAllElements.back();



        //1b - allocate  & initialize E2E mapping.
        m_uiE2EMapping.resize(m_uiAllElements.size()*m_uiNumDirections,LOOK_UP_TABLE_DEFAULT);
        SFC::seqSearch::SFC_treeSearch(&(*(m_uiKeys.begin())), &(*(m_uiAllElements.begin())), 0, m_uiKeys.size(), 0, m_uiAllElements.size(), m_uiMaxDepth, m_uiMaxDepth, 0);

        std::vector<unsigned int> *ownerList;
        std::vector<unsigned int> *stencilIndexDirection;
        unsigned int result;
        unsigned int dir;
        //unsigned int stencilIndex;
        Key *m_uiKeysPtr=&(*(m_uiKeys.begin()));

        m_uiE2EMapping.resize(m_uiAllElements.size()*m_uiNumDirections,LOOK_UP_TABLE_DEFAULT);

        /**Neighbour list is made from
        *  first x axis. (left to right)
        *  second y axis. (down to up)
        *  third z axis. (back to front)**/

        for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
            if(!(OCT_FOUND & m_uiKeysPtr[k].getFlag())) std::cout<<"key: "<<m_uiKeysPtr[k]<<" not found: "<<std::endl;
            assert((OCT_FOUND & m_uiKeysPtr[k].getFlag()));// Note that all the keys should be found locally due to the fact that we have exchanged ghost elements.
            ownerList = m_uiKeysPtr[k].getOwnerList();
            stencilIndexDirection = m_uiKeysPtr[k].getStencilIndexDirectionList();
            result = m_uiKeysPtr[k].getSearchResult();
            assert(m_uiAllElements[result].isAncestor(m_uiKeys[k]) || m_uiAllElements[result]==m_uiKeys[k]); // To check the result found in the treeSearch is correct or not.
            for (unsigned int w = 0; w < ownerList->size(); w++) {
                dir = ((*stencilIndexDirection)[w]) & KEY_DIR_OFFSET;
                //stencilIndex = (((*stencilIndexDirection)[w]) & (KS_MAX << 3u)) >> 3u;
                //std::cout<<"dir: "<<dir<<" stencil Index: "<<stencilIndex<<"owner index: "<<(*ownerList)[w]<<std::endl;
                m_uiE2EMapping[(((*ownerList)[w]) + m_uiElementLocalBegin)*m_uiNumDirections+dir] = result;
                // Note : Following is done to enforce that the local elements that points to ghost elements has the inverse mapping.
                if(result<m_uiElementLocalBegin || result>=m_uiElementLocalEnd) {
                    m_uiE2EMapping[result * m_uiNumDirections + (1u ^ dir)] = (((*ownerList)[w]) +
                                                                               m_uiElementLocalBegin); // Note This depends on the OCT_DIR numbering.
                    assert(false);// for sequential case this cannot be true.
                }


            }

        }

        assert(seq::test::checkE2EMapping(m_uiE2EMapping,m_uiAllElements,m_uiElementPreGhostBegin,m_uiElementPostGhostEnd,1,m_uiNumDirections));
#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(m_uiAllElements,m_uiActiveRank,"m_uiAllElements");
        //std::cout<<"rank: "<<rank<<"pre begin: "<<m_uiElementPreGhostBegin<<" pre end: "<<m_uiElementPreGhostEnd<<" local begin: "<<m_uiElementLocalBegin<<" local end: "<<m_uiElementLocalEnd<<" post begin: "<<m_uiElementPostGhostBegin<<" post end: "<<m_uiElementPostGhostEnd<<std::endl;
#endif

        // Note: Note that m_uiAlllNodes need not to be sorted and contains duplicates globally.


#ifdef DEBUG_MESH_GENERATION
      /*  std::vector<ot::TreeNode> missedKeys;
            for (unsigned int k = 0; k < m_uiKeys.size(); k++) {
                if (!(OCT_FOUND & m_uiKeys[k].getFlag())) {
                    std::cout << "rank: " << m_uiActiveRank << " key : " << m_uiKeys[k] << " not found!" << " Search index count: "
                              << searchResults.size() << std::endl;
                    missedKeys.push_back(m_uiKeys[k]);
                }
            }

        treeNodesTovtk(missedKeys,rank,"missedKeys");*/
#endif

        std::cout << "Seq: E2E mapping Ended" << std::endl;

    }

    void Mesh::buildE2NMap() {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

#ifdef DEBUG_E2N_MAPPING
        std::vector<ot::TreeNode> invalidatedPreGhost;
        std::vector<ot::TreeNode> invalidatedPostGhost;
        //if(!m_uiActiveRank)  std::cout<<"E2E  rank : "<<m_uiActiveRank<<std::endl;
        //if(!m_uiActiveRank)
            for(unsigned int e=0;e<m_uiAllElements.size();e++)
            {
                //if(m_uiAllElements[e].getLevel()!=1) {
                //std::cout << "Element : "<<e<<" " << m_uiAllElements[e] << " : Node List :";
                for (unsigned int k = 0; k < m_uiNumDirections; k++) {

                  //  std::cout << " " << m_uiE2EMapping[e * m_uiNumDirections + k];
                   if(m_uiE2EMapping[e * m_uiNumDirections + k]!=LOOK_UP_TABLE_DEFAULT) assert(m_uiE2EMapping[e * m_uiNumDirections + k]<m_uiAllElements.size());

                }

                //std::cout << std::endl;
                //}

            }
#endif

        // update the E2E mapping with fake elements.
        unsigned int lookUp = 0;
        unsigned int lev1 = 0;
        unsigned int lev2 = 0;

        unsigned int child;
        unsigned int parent;

#ifdef DEBUG_E2N_MAPPING
        for(unsigned int ge=m_uiElementPreGhostBegin;ge<m_uiElementPreGhostEnd;ge++)
        {
            for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
              if(m_uiE2EMapping[ge*m_uiNumDirections+dir]!=LOOK_UP_TABLE_DEFAULT)
                  assert((m_uiE2EMapping[ge*m_uiNumDirections+dir]>=m_uiElementLocalBegin) && (m_uiE2EMapping[ge*m_uiNumDirections+dir]<m_uiElementLocalEnd) );
        }

        for(unsigned int ge=m_uiElementPostGhostBegin;ge<m_uiElementPostGhostEnd;ge++)
        {
            for(unsigned int dir=0;dir<m_uiNumDirections;dir++)
                if(m_uiE2EMapping[ge*m_uiNumDirections+dir]!=LOOK_UP_TABLE_DEFAULT)
                    assert((m_uiE2EMapping[ge*m_uiNumDirections+dir]>=m_uiElementLocalBegin) && (m_uiE2EMapping[ge*m_uiNumDirections+dir]<m_uiElementLocalEnd) );
        }

#endif
        assert(m_uiNumTotalElements == m_uiAllElements.size());
        assert((m_uiElementPostGhostEnd - m_uiElementPreGhostBegin) > 0);
        assert(m_uiNumTotalElements == ((m_uiElementPostGhostEnd - m_uiElementPreGhostBegin)));

        m_uiE2NMapping_CG.resize(m_uiNumTotalElements * m_uiNpE);
        m_uiE2NMapping_DG.resize(m_uiNumTotalElements * m_uiNpE);

        // initialize the DG mapping. // this order is mandotory.
        for (unsigned int e = 0; e < (m_uiNumTotalElements); e++)
            for (unsigned int k = 0; k < (m_uiElementOrder + 1); k++) //z coordinate
                for (unsigned int j = 0; j < (m_uiElementOrder + 1); j++) // y coordinate
                    for (unsigned int i = 0; i < (m_uiElementOrder + 1); i++) // x coordinate
                        m_uiE2NMapping_CG[e * m_uiNpE + k * (m_uiElementOrder + 1) * (m_uiElementOrder + 1) +
                                          j * (m_uiElementOrder + 1) + i] =
                                e * m_uiNpE + k * (m_uiElementOrder + 1) * (m_uiElementOrder + 1) +
                                j * (m_uiElementOrder + 1) + i;


#ifdef DEBUG_E2N_MAPPING
        MPI_Barrier(MPI_COMM_WORLD);
        if (!m_uiActiveRank) std::cout << "Invalid nodes removed " << std::endl;
        //treeNodesTovtk(invalidatedPreGhost,m_uiActiveRank,"invalidPre");
        //treeNodesTovtk(invalidatedPostGhost,m_uiActiveRank,"invalidPost");
#endif

        // 2. Removing the duplicate nodes from the mapping.
        unsigned int ownerIndexChild;
        unsigned int ownerIndexParent;

        unsigned int child_i,child_j,child_k;
        unsigned int parent_i,parent_j,parent_k;


#ifdef DEBUG_E2N_MAPPING
        std::vector<ot::TreeNode> cusE2ECheck;

        if(!m_uiActiveRank && m_uiActiveNpes>1){

            unsigned int eleID=2;
            //cusE2ECheck.push_back(m_uiAllElements[eleID]);
            for(unsigned int dir=0;dir<(m_uiNumDirections);dir++){
                if(m_uiE2EMapping[eleID*m_uiNumDirections+dir]!=LOOK_UP_TABLE_DEFAULT)
                    cusE2ECheck.push_back(m_uiAllElements[m_uiE2EMapping[eleID*m_uiNumDirections+dir]]);
            }

            treeNodesTovtk(cusE2ECheck,m_uiActiveRank,"cusE2ECheck");
        }

        unsigned int eleVtk=2;
        unsigned int eleVtkAt=1;

        unsigned  int lookUp1,lookUp2,lookUp3,lookUp4;
        unsigned int edgeOwner;
        unsigned int cornerNodeOwner;

        unsigned int cornerNodeOwnerIndex;
        unsigned int cornerNodeChildIndex;
#endif

        bool parentChildLevEqual=false;
        std::vector<unsigned int> faceChildIndex; // indices that are being updated for a face
        std::vector<unsigned int> faceOwnerIndex; // indices that are being used for updates. for a face

        std::vector<unsigned int> edgeChildIndex; // indices that are being updated for an edge
        std::vector<unsigned int> edgeOwnerIndex; // indices that are being used for updates. for an edge



        for (unsigned int e = m_uiElementPreGhostBegin; e < (m_uiElementPostGhostEnd); e++) {


            // 1. All local nodes for a given element is automatically mapped for a given element by construction of the DG indexing.

            //2. Map internal nodes on each face,  with the corresponding face.
            parentChildLevEqual=false;
            lookUp = m_uiE2EMapping[e * m_uiNumDirections + OCT_DIR_LEFT];
            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
            {
                lev1 = m_uiAllElements[e].getLevel();
                lev2 = m_uiAllElements[lookUp].getLevel();
                if (lev1 == lev2) {
                    parentChildLevEqual=true;
                    lev1 = e;
                    lev2 = lookUp;
                    assert(e != lookUp);
                }

                child = ((lev1 > lev2) ? e : lookUp);
                parent = ((lev1 < lev2) ? e : lookUp);

               if(child==e){

                   //if(m_uiActiveRank==0 && e==324)std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" LEFT ele: "<<e<<std::endl;
                   assert(parent==lookUp);

                   faceNodesIndex(child,OCT_DIR_LEFT,faceChildIndex,true);
                   faceNodesIndex(parent,OCT_DIR_RIGHT,faceOwnerIndex,true);

                   assert(faceChildIndex.size()==faceOwnerIndex.size());

                   for(unsigned int index=0;index<faceChildIndex.size();index++)
                       m_uiE2NMapping_CG[faceChildIndex[index]]=m_uiE2NMapping_CG[faceOwnerIndex[index]];


                   OCT_DIR_LEFT_INTERNAL_EDGE_MAP(child,parent,parentChildLevEqual,edgeChildIndex,edgeOwnerIndex); // maps all the edges in the left face


               }


            }


            parentChildLevEqual=false;
            lookUp = m_uiE2EMapping[e * m_uiNumDirections + OCT_DIR_RIGHT];
            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
            {
                lev1 = m_uiAllElements[e].getLevel();
                lev2 = m_uiAllElements[lookUp].getLevel();
                if (lev1 == lev2) {
                    lev1 = e;
                    lev2 = lookUp;
                    parentChildLevEqual=true;
                    assert(e != lookUp);
                }

                child = ((lev1 > lev2) ? e : lookUp);
                parent = ((lev1 < lev2) ? e : lookUp);

               if(child==e)
                {
                    //if(m_uiActiveRank==0 && e==324)std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" RIGHT ele: "<<e<<std::endl;
                    assert(parent==lookUp);
                    faceNodesIndex(child,OCT_DIR_RIGHT,faceChildIndex,true);
                    faceNodesIndex(parent,OCT_DIR_LEFT,faceOwnerIndex,true);

                    assert(faceChildIndex.size()==faceOwnerIndex.size());

                    for(unsigned int index=0;index<faceChildIndex.size();index++)
                        m_uiE2NMapping_CG[faceChildIndex[index]]=m_uiE2NMapping_CG[faceOwnerIndex[index]];


                    OCT_DIR_RIGHT_INTERNAL_EDGE_MAP(child,parent,parentChildLevEqual,edgeChildIndex,edgeOwnerIndex); // maps all the edges in the right face


                }


            }



            parentChildLevEqual=false;
            lookUp = m_uiE2EMapping[e * m_uiNumDirections + OCT_DIR_DOWN];
            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
            {
                lev1 = m_uiAllElements[e].getLevel();
                lev2 = m_uiAllElements[lookUp].getLevel();
                if (lev1 == lev2) {
                    lev1 = e;
                    lev2 = lookUp;
                    parentChildLevEqual=true;
                    assert(e != lookUp);
                }

                child = ((lev1 > lev2) ? e : lookUp);
                parent = ((lev1 < lev2) ? e : lookUp);

                if(child==e)
                {
                    //if(m_uiActiveRank==0 && e==324)std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" DOWN ele: "<<e<<std::endl;
                    assert(parent==lookUp);
                    faceNodesIndex(child,OCT_DIR_DOWN,faceChildIndex,true);
                    faceNodesIndex(parent,OCT_DIR_UP,faceOwnerIndex,true);

                    assert(faceChildIndex.size()==faceOwnerIndex.size());

                    for(unsigned int index=0;index<faceChildIndex.size();index++)
                        m_uiE2NMapping_CG[faceChildIndex[index]]=m_uiE2NMapping_CG[faceOwnerIndex[index]];



                    OCT_DIR_DOWN_INTERNAL_EDGE_MAP(child,parent,parentChildLevEqual,edgeChildIndex,edgeOwnerIndex); // maps all the edges in the DOWN face


                }


            }



            parentChildLevEqual=false;
            lookUp = m_uiE2EMapping[e * m_uiNumDirections + OCT_DIR_UP];
            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
            {
                lev1 = m_uiAllElements[e].getLevel();
                lev2 = m_uiAllElements[lookUp].getLevel();
                if (lev1 == lev2) {
                    lev1 = e;
                    lev2 = lookUp;
                    parentChildLevEqual=true;
                    assert(e != lookUp);
                }

                child = ((lev1 > lev2) ? e : lookUp);
                parent = ((lev1 < lev2) ? e : lookUp);

                if(child==e)
                {
                    //if(m_uiActiveRank==0 && e==324)std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" UP ele: "<<e<<std::endl;
                    assert(parent==lookUp);
                    faceNodesIndex(child,OCT_DIR_UP,faceChildIndex,true);
                    faceNodesIndex(parent,OCT_DIR_DOWN,faceOwnerIndex,true);

                    assert(faceChildIndex.size()==faceOwnerIndex.size());

                    for(unsigned int index=0;index<faceChildIndex.size();index++)
                        m_uiE2NMapping_CG[faceChildIndex[index]]=m_uiE2NMapping_CG[faceOwnerIndex[index]];


                    OCT_DIR_UP_INTERNAL_EDGE_MAP(child,parent,parentChildLevEqual,edgeChildIndex,edgeOwnerIndex); // maps all the edges in the UP face


                }



            }


            parentChildLevEqual=false;
            lookUp = m_uiE2EMapping[e * m_uiNumDirections + OCT_DIR_BACK];
            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
            {
                lev1 = m_uiAllElements[e].getLevel();
                lev2 = m_uiAllElements[lookUp].getLevel();
                if (lev1 == lev2) {
                    lev1 = e;
                    lev2 = lookUp;
                    parentChildLevEqual=true;
                    assert(e != lookUp);
                }

                child = ((lev1 > lev2) ? e : lookUp);
                parent = ((lev1 < lev2) ? e : lookUp);

                if(child==e)
                {
                    //if(m_uiActiveRank==0 && e==324)std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" BACK ele: "<<e<<std::endl;
                    assert(parent==lookUp);
                    faceNodesIndex(child,OCT_DIR_BACK,faceChildIndex,true);
                    faceNodesIndex(parent,OCT_DIR_FRONT,faceOwnerIndex,true);

                    assert(faceChildIndex.size()==faceOwnerIndex.size());

                    for(unsigned int index=0;index<faceChildIndex.size();index++)
                        m_uiE2NMapping_CG[faceChildIndex[index]]=m_uiE2NMapping_CG[faceOwnerIndex[index]];

                    OCT_DIR_BACK_INTERNAL_EDGE_MAP(child,parent,parentChildLevEqual,edgeChildIndex,edgeOwnerIndex); // maps all the edges in the back face


                }


            }



            parentChildLevEqual=false;
            lookUp = m_uiE2EMapping[e * m_uiNumDirections + OCT_DIR_FRONT];
            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
            {
                lev1 = m_uiAllElements[e].getLevel();
                lev2 = m_uiAllElements[lookUp].getLevel();
                if (lev1 == lev2) {
                    lev1 = e;
                    lev2 = lookUp;
                    parentChildLevEqual=true;
                    assert(e != lookUp);
                }

                child = ((lev1 > lev2) ? e : lookUp);
                parent = ((lev1 < lev2) ? e : lookUp);

                if(child==e)
                {

                    //if(m_uiActiveRank==0 && e==324)std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<"FRONT ele: "<<e<<std::endl;

                    assert(parent==lookUp);
                    faceNodesIndex(child,OCT_DIR_FRONT,faceChildIndex,true);
                    faceNodesIndex(parent,OCT_DIR_BACK,faceOwnerIndex,true);

                    assert(faceChildIndex.size()==faceOwnerIndex.size());

                    for(unsigned int index=0;index<faceChildIndex.size();index++)
                        m_uiE2NMapping_CG[faceChildIndex[index]]=m_uiE2NMapping_CG[faceOwnerIndex[index]];

                    OCT_DIR_FRONT_INTERNAL_EDGE_MAP(child,parent,parentChildLevEqual,edgeChildIndex,edgeOwnerIndex); // maps all the edges in the front face


                }


            }



            CORNER_NODE_MAP(e);


            /*if(m_uiActiveRank==0 && e==284)
            {

                std::vector<ot::TreeNode> cusEleCheck;
                unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
                unsigned int x,y,z,sz;
                cusEleCheck.push_back(m_uiAllElements[e]);
                for(unsigned int node=0;node<m_uiNpE;node++)
                {

                    dg2eijk(m_uiE2NMapping_CG[e*m_uiNpE+node],ownerID,ii_x,jj_y,kk_z);

                    x=m_uiAllElements[ownerID].getX();
                    y=m_uiAllElements[ownerID].getY();
                    z=m_uiAllElements[ownerID].getZ();
                    sz=1u<<(m_uiMaxDepth-m_uiAllElements[ownerID].getLevel());
                    cusEleCheck.push_back(m_uiAllElements[ownerID]);

                    cusEleCheck.push_back(ot::TreeNode((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth,m_uiDim, m_uiMaxDepth));

                }

                treeNodesTovtk(cusEleCheck,e,"cusE2N_1");

            }*/


        }



       /* for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPostGhostEnd;e++)
        {
            if(m_uiActiveRank==0 && e==284)
            {

                    std::vector<ot::TreeNode> cusEleCheck;
                    unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
                    unsigned int x,y,z,sz;
                    cusEleCheck.push_back(m_uiAllElements[e]);
                    for(unsigned int node=0;node<m_uiNpE;node++)
                    {

                        dg2eijk(m_uiE2NMapping_CG[e*m_uiNpE+node],ownerID,ii_x,jj_y,kk_z);

                        x=m_uiAllElements[ownerID].getX();
                        y=m_uiAllElements[ownerID].getY();
                        z=m_uiAllElements[ownerID].getZ();
                        sz=1u<<(m_uiMaxDepth-m_uiAllElements[ownerID].getLevel());
                        cusEleCheck.push_back(m_uiAllElements[ownerID]);

                        cusEleCheck.push_back(ot::TreeNode((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth,m_uiDim, m_uiMaxDepth));

                    }

                    treeNodesTovtk(cusEleCheck,e,"cusE2N_2");

            }

        }*/



#ifdef DEBUG_E2N_MAPPING
         MPI_Barrier(MPI_COMM_WORLD);
         unsigned int eleIndex;
         if(m_uiActiveRank==0)  std::cout<<"E2N  rank : "<<m_uiActiveRank<<std::endl;
         if(m_uiActiveRank==0)
         for(unsigned int e=0;e<m_uiAllElements.size();e++)
         {
             //if(m_uiAllElements[e].getLevel()!=1) {
                 std::cout << "Element : "<<e<<" " << m_uiAllElements[e] << " : Node List :";
                 for (unsigned int k = 0; k < m_uiNpE; k++) {

                     std::cout << " " << m_uiE2NMapping_CG[e * m_uiNpE + k];
          }

                 std::cout << std::endl;
             //}

         }

#endif




        //assert(seq::test::checkE2NMapping(m_uiE2EMapping, m_uiE2NMapping_CG,m_uiAllElements,m_uiNumDirections,m_uiElementOrder));
        std::vector<unsigned int > E2N_DG_Sorted;
        std::vector<unsigned int> dg2dg_p; // dg to dg prime
        E2N_DG_Sorted.resize(m_uiE2NMapping_CG.size());
        E2N_DG_Sorted.assign(m_uiE2NMapping_CG.begin(),m_uiE2NMapping_CG.end());

        m_uiDG2CG.resize(m_uiAllElements.size()*m_uiNpE,LOOK_UP_TABLE_DEFAULT);
        dg2dg_p.resize(m_uiAllElements.size()*m_uiNpE,LOOK_UP_TABLE_DEFAULT);

// 3. Update DG indexing with CG indexing.
        std::sort( E2N_DG_Sorted.begin(), E2N_DG_Sorted.end() );
        E2N_DG_Sorted.erase( std::unique( E2N_DG_Sorted.begin(), E2N_DG_Sorted.end() ), E2N_DG_Sorted.end() );

        unsigned int owner1,ii_x1,jj_y1,kk_z1;
        unsigned int owner2,ii_x2,jj_y2,kk_z2;
        unsigned int old_val;
        unsigned int new_val;
        unsigned int nsz;

        SearchKey tmpSKey;
        Key tmpKey;
        SearchKey rootSKey(m_uiDim,m_uiMaxDepth);
        std::vector<SearchKey> tmpSKeys;
        std::vector<Key> cgNodes;
        std::vector<SearchKey> skeys_cg;
        std::vector<SearchKey>::iterator hintSKey;
        unsigned int skip=1;
        unsigned int i_cg,i_dg;
        std::vector<unsigned int> * ownerList_ptr;

        for(unsigned int index=0;index<E2N_DG_Sorted.size();index++)
        {
            dg2eijk(E2N_DG_Sorted[index],owner1,ii_x1,jj_y1,kk_z1);
            assert(owner1<m_uiAllElements.size());
            nsz=1u<<(m_uiMaxDepth-m_uiAllElements[owner1].getLevel());
            assert(nsz%m_uiElementOrder==0);
            hintSKey=skeys_cg.emplace(skeys_cg.end(),SearchKey((m_uiAllElements[owner1].getX())+(ii_x1*nsz/m_uiElementOrder),(m_uiAllElements[owner1].getY())+(jj_y1*nsz/m_uiElementOrder),(m_uiAllElements[owner1].getZ())+(kk_z1*nsz/m_uiElementOrder),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
            hintSKey->addOwner(E2N_DG_Sorted[index]);


        }


        SFC::seqSort::SFC_treeSort(&(*(skeys_cg.begin())),skeys_cg.size(),tmpSKeys,tmpSKeys,tmpSKeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        for(unsigned int e=0;e<skeys_cg.size();e++)
        {
            skip=1;
            tmpSKey=skeys_cg[e];

            tmpKey=Key(skeys_cg[e].getX(),skeys_cg[e].getY(),skeys_cg[e].getZ(),skeys_cg[e].getLevel(),m_uiDim,m_uiMaxDepth+1);
            tmpKey.addOwner(skeys_cg[e].getOwner());

            while(((e+skip)<skeys_cg.size()) && (skeys_cg[e]==skeys_cg[e+skip]))
            {
                tmpKey.addOwner(skeys_cg[e+skip].getOwner());
                dg2eijk(tmpSKey.getOwner(),owner1,ii_x1,jj_y1,kk_z1);
                dg2eijk(skeys_cg[e+skip].getOwner(),owner2,ii_x2,jj_y2,kk_z2);

                lev1=m_uiAllElements[owner1].getLevel();
                lev2=m_uiAllElements[owner2].getLevel();

                if(lev1==lev2) {
                    lev1 = owner1;
                    lev2 = owner2;
                }
                assert(lev1!=lev2);

                if(lev1>lev2)
                    tmpSKey.addOwner(skeys_cg[e+skip].getOwner());



                skip++;

            }
            m_uiCG2DG.push_back(tmpSKey.getOwner());
            tmpKey.setSearchResult(tmpSKey.getOwner());
            cgNodes.push_back(tmpKey);
            e+=(skip-1);

        }

        std::sort(m_uiCG2DG.begin(),m_uiCG2DG.end());


        for(unsigned int i=0;i<m_uiCG2DG.size();i++)
            m_uiDG2CG[m_uiCG2DG[i]]=i;

        for(unsigned int i=0;i<cgNodes.size();i++)
        {
            ownerList_ptr=cgNodes[i].getOwnerList();
            if(ownerList_ptr->size()>1)
            {
                i_cg=(std::lower_bound(m_uiCG2DG.begin(), m_uiCG2DG.end(), cgNodes[i].getSearchResult()) -m_uiCG2DG.begin());
                assert(i_cg < m_uiCG2DG.size());
                for(unsigned int w=0;w<ownerList_ptr->size();w++)
                {
                    m_uiDG2CG[(*(ownerList_ptr))[w]]=i_cg;
                    dg2dg_p[(*(ownerList_ptr))[w]]=cgNodes[i].getSearchResult();
                }
            }

        }

        for(unsigned int i=0;i<m_uiE2NMapping_CG.size();i++)
         if(dg2dg_p[m_uiE2NMapping_CG[i]]!=LOOK_UP_TABLE_DEFAULT) m_uiE2NMapping_CG[i]=dg2dg_p[m_uiE2NMapping_CG[i]];






#ifdef DEBUG_E2N_MAPPING
        //MPI_Barrier(MPI_COMM_WORLD);
        if(!m_uiActiveRank) std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<"Number of actual nodes: "<<(E2N_DG_Sorted.size())<<std::endl;
#endif

        m_uiNodePreGhostBegin=UINT_MAX;
        m_uiNodeLocalBegin=UINT_MAX;
        m_uiNodePostGhostBegin=UINT_MAX;

        unsigned int preOwner=UINT_MAX;
        unsigned int localOwner=UINT_MAX;
        unsigned int postOwner=UINT_MAX;




        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPostGhostEnd;e++)
        {
            unsigned int tmpIndex;
            for(unsigned int k=0;k<m_uiNpE;k++)
            {

                tmpIndex = (m_uiE2NMapping_CG[e * m_uiNpE + k]/m_uiNpE);
                assert(tmpIndex==(((m_uiE2NMapping_CG[e * m_uiNpE + k]) / (m_uiElementOrder + 1)) /
                                                                              (m_uiElementOrder + 1)) / (m_uiElementOrder + 1));
                if ((tmpIndex >= m_uiElementPreGhostBegin) && (tmpIndex < m_uiElementPreGhostEnd) && /*(preOwner>=(m_uiE2NMapping_CG[e * m_uiNpE + k])/m_uiNpE) &&*/ (m_uiNodePreGhostBegin>m_uiE2NMapping_CG[e * m_uiNpE + k])){
                    //preOwner = m_uiE2NMapping_CG[e * m_uiNpE + k]/m_uiNpE;
                    m_uiNodePreGhostBegin = m_uiE2NMapping_CG[e * m_uiNpE + k];
                }

                if ((tmpIndex >= m_uiElementLocalBegin) && (tmpIndex < m_uiElementLocalEnd) && /*(localOwner >=(m_uiE2NMapping_CG[e * m_uiNpE + k])/m_uiNpE) &&*/ (m_uiNodeLocalBegin > m_uiE2NMapping_CG[e * m_uiNpE + k])) {
                    //localOwner = m_uiE2NMapping_CG[e * m_uiNpE + k]/m_uiNpE;
                    m_uiNodeLocalBegin = m_uiE2NMapping_CG[e * m_uiNpE + k];
                }

                if ((tmpIndex >= m_uiElementPostGhostBegin) && (tmpIndex < m_uiElementPostGhostEnd) && /*(postOwner >=(m_uiE2NMapping_CG[e * m_uiNpE + k])/m_uiNpE) &&*/ (m_uiNodePostGhostBegin>m_uiE2NMapping_CG[e * m_uiNpE + k])) {
                    //postOwner = m_uiE2NMapping_CG[e * m_uiNpE + k]/m_uiNpE;
                    m_uiNodePostGhostBegin = m_uiE2NMapping_CG[e * m_uiNpE + k];
                }

            }

        }


        assert(m_uiNodeLocalBegin!=UINT_MAX); // local node begin should be found.
        assert(m_uiDG2CG[m_uiNodeLocalBegin]!=LOOK_UP_TABLE_DEFAULT);
        m_uiNodeLocalBegin=m_uiDG2CG[m_uiNodeLocalBegin];//(std::lower_bound(E2N_DG_Sorted.begin(),E2N_DG_Sorted.end(),m_uiNodeLocalBegin)-E2N_DG_Sorted.begin());
        if(m_uiNodePreGhostBegin==UINT_MAX) {
            m_uiNodePreGhostBegin=0;
            m_uiNodePreGhostEnd=0;
            assert(m_uiNodeLocalBegin==0);
        }
        else{
            assert(m_uiDG2CG[m_uiNodePreGhostBegin]!=LOOK_UP_TABLE_DEFAULT);
            m_uiNodePreGhostBegin=m_uiDG2CG[m_uiNodePreGhostBegin];//(std::lower_bound(E2N_DG_Sorted.begin(),E2N_DG_Sorted.end(),m_uiNodePreGhostBegin)-E2N_DG_Sorted.begin());
            m_uiNodePreGhostEnd=m_uiNodeLocalBegin;
        }

        if(m_uiNodePostGhostBegin==UINT_MAX) {
            m_uiNodeLocalEnd=m_uiCG2DG.size(); //E2N_DG_Sorted.size();
            m_uiNodePostGhostBegin=m_uiNodeLocalEnd;
            m_uiNodePostGhostEnd=m_uiNodeLocalEnd;
        }
        else
        {
            assert(m_uiDG2CG[m_uiNodePostGhostBegin]!=LOOK_UP_TABLE_DEFAULT);
            m_uiNodePostGhostBegin=m_uiDG2CG[m_uiNodePostGhostBegin];//(std::lower_bound(E2N_DG_Sorted.begin(),E2N_DG_Sorted.end(),m_uiNodePostGhostBegin)-E2N_DG_Sorted.begin());
            m_uiNodeLocalEnd=m_uiNodePostGhostBegin;
            m_uiNodePostGhostEnd=m_uiCG2DG.size();//E2N_DG_Sorted.size();

        }

        m_uiNumActualNodes=cgNodes.size();

        m_uiE2NMapping_DG.assign(m_uiE2NMapping_CG.begin(),m_uiE2NMapping_CG.end());

        for(unsigned int i=0;i<m_uiE2NMapping_CG.size();i++)
        {
            assert(m_uiDG2CG[m_uiE2NMapping_CG[i]]!=LOOK_UP_TABLE_DEFAULT);
            m_uiE2NMapping_CG[i]=m_uiDG2CG[m_uiE2NMapping_CG[i]];
        }


        dg2dg_p.clear();
        E2N_DG_Sorted.clear();


#ifdef DEBUG_E2N_MAPPING
        MPI_Barrier(MPI_COMM_WORLD);
        if(m_uiActiveRank) std::cout<<" DG to CG index updated "<<std::endl;
#endif



    /*    for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPostGhostEnd;e++)
        {
            if(m_uiActiveRank==2 && e==28)
            {

                std::vector<ot::TreeNode> cusEleCheck;
                unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
                unsigned int x,y,z,sz;
                cusEleCheck.push_back(m_uiAllElements[e]);
                for(unsigned int node=0;node<m_uiNpE;node++)
                {

                    dg2eijk(m_uiE2NMapping_DG[e*m_uiNpE+node],ownerID,ii_x,jj_y,kk_z);

                    x=m_uiAllElements[ownerID].getX();
                    y=m_uiAllElements[ownerID].getY();
                    z=m_uiAllElements[ownerID].getZ();
                    sz=1u<<(m_uiMaxDepth-m_uiAllElements[ownerID].getLevel());
                    cusEleCheck.push_back(m_uiAllElements[ownerID]);

                    cusEleCheck.push_back(ot::TreeNode((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth,m_uiDim, m_uiMaxDepth));

                }

                treeNodesTovtk(cusEleCheck,e,"cusE2N_3");

            }

        }*/




          /*unsigned int eleIndex;
          if(!m_uiActiveRank)  std::cout<<"E2N  rank : "<<m_uiActiveRank<<std::endl;
          if(!m_uiActiveRank)
              for(unsigned int e=0;e<m_uiAllElements.size();e++)
              {
                  if(m_uiAllElements[e]==ot::TreeNode(24, 12, 40, 4,m_uiDim,m_uiMaxDepth)) {
                  std::cout << "rank: "<<m_uiActiveRank<<" Element : "<<e<<" " << m_uiAllElements[e] << " : Node List :";
                  for (unsigned int k = 0; k < m_uiNpE; k++) {

                      std::cout << " " << m_uiE2NMapping_DG[e * m_uiNpE + k];

                  }

                  std::cout << std::endl;
                  }

              }*/
//--------------------------------------------------------PRINT THE E2N MAP------------------------------------------------------------------------------------
        /*for(unsigned int w=0;w<m_uiE2NMapping_CG.size();w++)
            std::cout<<"w: "<<w<<" -> : "<<m_uiE2NMapping_CG[w]<<std::endl;*/

//--------------------------------------------------------PRINT THE E2N MAP------------------------------------------------------------------------------------


#ifdef DEBUG_E2N_MAPPING
        //MPI_Barrier(MPI_COMM_WORLD);
        if(!m_uiActiveRank) std::cout<<"[NODE] rank:  "<<m_uiActiveRank<<" pre ( "<<m_uiNodePreGhostBegin<<", "<<m_uiNodePreGhostEnd<<") local ( "<<m_uiNodeLocalBegin<<", "<<m_uiNodeLocalEnd<<")"<<" post ("<<m_uiNodePostGhostBegin<<" , "<<m_uiNodePostGhostEnd<<")"<<std::endl;
        if(!m_uiActiveRank) std::cout<<"[ELEMENT] rank:  "<<m_uiActiveRank<<" pre ( "<<m_uiElementPreGhostBegin<<", "<<m_uiElementPreGhostEnd<<") local ( "<<m_uiElementLocalBegin<<", "<<m_uiElementLocalEnd<<")"<<" post ("<<m_uiElementPostGhostBegin<<" , "<<m_uiElementPostGhostEnd<<")"<<std::endl;

        if(m_uiActiveRank) std::cout<<"[NODE] rank:  "<<m_uiActiveRank<<" pre ( "<<m_uiNodePreGhostBegin<<", "<<m_uiNodePreGhostEnd<<") local ( "<<m_uiNodeLocalBegin<<", "<<m_uiNodeLocalEnd<<")"<<" post ("<<m_uiNodePostGhostBegin<<" , "<<m_uiNodePostGhostEnd<<")"<<std::endl;
        if(m_uiActiveRank) std::cout<<"[ELEMENT] rank:  "<<m_uiActiveRank<<" pre ( "<<m_uiElementPreGhostBegin<<", "<<m_uiElementPreGhostEnd<<") local ( "<<m_uiElementLocalBegin<<", "<<m_uiElementLocalEnd<<")"<<" post ("<<m_uiElementPostGhostBegin<<" , "<<m_uiElementPostGhostEnd<<")"<<std::endl;

#endif

//---------------------------------------print out the E2N mapping of fake elements. (This is done for only the fake elements. ) ---------------------------------------------------------------------------------


       /* MPI_Barrier(MPI_COMM_WORLD);
        if(!m_uiActiveRank){
            std::cout<<"rank: "<<m_uiActiveRank<<"fake element e2n mapping. "<<std::endl;
            std::cout<<"number of Fake Elements : "<<fakeElements_vec.size()<<std::endl;
            std::cout<<"number of FakeElement Nodes: "<<m_uiNumFakeNodes<<std::endl;
            std::cout<<"[Fake ELEMENT] rank:  "<<m_uiActiveRank<<" pre ( "<<m_uiFElementPreGhostBegin<<", "<<m_uiFElementPreGhostEnd<<") local ( "<<m_uiFElementLocalBegin<<", "<<m_uiFElementLocalEnd<<")"<<" post ("<<m_uiFElementPostGhostBegin<<" , "<<m_uiFElementPostGhostEnd<<")"<<std::endl;
          for(unsigned int e=0;e<fakeElements_vec.size();e++)
          {
                  std::cout << "Element : "<<e<<" " << fakeElements_vec[e] << " : Node List :";
                  for (unsigned int k = 0; k < m_uiNpE; k++) {

                      std::cout << " " << fakeElement2Node_CG[e * m_uiNpE + k];

                  }

                  std::cout << std::endl;


          }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if(m_uiActiveRank==1){
            std::cout<<"rank: "<<m_uiActiveRank<<"fake element e2n mapping. "<<std::endl;
            std::cout<<"number of Fake Elements : "<<fakeElements_vec.size()<<std::endl;
            std::cout<<"number of FakeElement Nodes: "<<m_uiNumFakeNodes<<std::endl;
            std::cout<<"[Fake ELEMENT] rank:  "<<m_uiActiveRank<<" pre ( "<<m_uiFElementPreGhostBegin<<", "<<m_uiFElementPreGhostEnd<<") local ( "<<m_uiFElementLocalBegin<<", "<<m_uiFElementLocalEnd<<")"<<" post ("<<m_uiFElementPostGhostBegin<<" , "<<m_uiFElementPostGhostEnd<<")"<<std::endl;
            for(unsigned int e=0;e<fakeElements_vec.size();e++)
            {
                std::cout << "Element : "<<e<<" " << fakeElements_vec[e] << " : Node List :";
                for (unsigned int k = 0; k < m_uiNpE; k++) {

                    std::cout << " " << fakeElement2Node_CG[e * m_uiNpE + k];

                }

                std::cout << std::endl;


            }
        }
        MPI_Barrier(MPI_COMM_WORLD);*/

//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        // 8. Change the size of the original E2N mapping and copy fake element to node mapping at the end of the actual element to node mapping.
        //m_uiE2NMapping_CG.resize(m_uiNumTotalElements*m_uiNpE+(m_uiFElementPostGhostEnd-m_uiFElementPreGhostBegin)*m_uiNpE);
        //memcpy(&(*(m_uiE2NMapping_CG.begin()+(m_uiNumTotalElements*m_uiNpE))),&(*(fakeElement2Node_CG.begin())),sizeof(unsigned int )*(m_uiFElementPostGhostEnd-m_uiFElementPreGhostBegin)*m_uiNpE);

//---------------------------------------print out the final e2n mapping of all, actual and fake element to node mapping.--------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        if(!m_uiActiveRank) std::cout<<"E2N Mapping ended"<<std::endl;


    }


    void Mesh::computeNodeScatterMaps(MPI_Comm comm)
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);


        ot::TreeNode minMaxLocalNode[2];

        ot::TreeNode rootNode(0,0,0,0,m_uiDim,m_uiMaxDepth);
        std::vector<ot::TreeNode>tmpNode;


        unsigned int x,y,z,sz; // x y z and size of an octant.
        unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
        std::set<unsigned int > nodeIndexVisited; // To keep track of the nodes already include in the sendNode maps.
        std::pair<std::set<unsigned int >::iterator,bool> setHintUint;

        unsigned int nodeIndex;
        unsigned int nodeIndex_DG;
        unsigned int elementLookUp;

#ifdef DEBUG_E2N_MAPPING_SM
        std::vector<ot::TreeNode> cusEleCheck;
         /*  for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
        {

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++) {

                        if(m_uiActiveRank==0)
                        {
                            dg2eijk(m_uiE2NMapping_DG[ele * m_uiNpE +
                                                      k * (m_uiElementOrder + 1) * (m_uiElementOrder + 1) +
                                                      j * (m_uiElementOrder + 1) + i], ownerID, ii_x, jj_y, kk_z);
                            x = m_uiAllElements[ownerID].getX();
                            y = m_uiAllElements[ownerID].getY();
                            z = m_uiAllElements[ownerID].getZ();
                            sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());

                            cusEleCheck.push_back(ot::TreeNode((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth,m_uiDim, m_uiMaxDepth));

                        }


                    }

            if(m_uiActiveRank==0){
                cusEleCheck.push_back(m_uiAllElements[ele]);
                treeNodesTovtk(cusEleCheck,ele,"cusEleCheck");
                cusEleCheck.clear();
            }


        }*/
#endif


        // 1. compute the local nodes.
        ///@todo We need not to generate all the nodes (which is expensive). we just need to find min and the max nodes.

        for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
        {

            for(unsigned int i=0;i<m_uiElementOrder+1;i++)
            for(unsigned int j=0;j<m_uiElementOrder+1;j++)
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
            {
                if(((i>1) && (i<(m_uiElementOrder-1))) || ((j>1) && (j<(m_uiElementOrder-1))) || ((k>1) && (k<(m_uiElementOrder-1)))) continue;
                nodeIndex=m_uiE2NMapping_CG[ele*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                if(nodeIndex>=m_uiNodeLocalBegin && nodeIndex<m_uiNodeLocalEnd){
                    setHintUint=nodeIndexVisited.emplace(nodeIndex);
                    if(setHintUint.second)
                    {   nodeIndex_DG=m_uiE2NMapping_DG[ele*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        dg2eijk(nodeIndex_DG,ownerID,ii_x,jj_y,kk_z);
                        x = m_uiAllElements[ownerID].getX();
                        y = m_uiAllElements[ownerID].getY();
                        z = m_uiAllElements[ownerID].getZ();
                        sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
                        assert(sz%m_uiElementOrder==0);
                        m_uiAllLocalNode.push_back(ot::TreeNode((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                    }
                }
            }
        }



        // 2. Find the local max of the local nodes to compute local splitters.
        ///@todo: Can use optimal Local node computation to find the min max.
        //NOTE: We need to use m_uiMaxDepth+1 to sort because this generates the all the possible nodes. Hence it contains 1u<<m_uiMaxDepth for (x y z) values.
        //SFC::seqSort::SFC_treeSort(&(*(m_uiAllLocalNode.begin())),m_uiAllLocalNode.size(),tmpNode,tmpNode,tmpNode,m_uiMaxDepth+1,m_uiMaxDepth+1,rootNode,ROOT_ROTATION,1,TS_SORT_ONLY);
        //SFC::seqSort::SFC_treeSortLocalOptimal(&(*(m_uiAllLocalNode.begin())),m_uiAllLocalNode.size(),m_uiMaxDepth+1,m_uiMaxDepth+1,rootNode,ROOT_ROTATION,true,minMaxLocalNode[0]);
        //SFC::seqSort::SFC_treeSortLocalOptimal(&(*(m_uiAllLocalNode.begin())),m_uiAllLocalNode.size(),m_uiMaxDepth+1,m_uiMaxDepth+1,rootNode,ROOT_ROTATION,false,minMaxLocalNode[1]);
        SFC::seqSort::SFC_treeSortLocalOptimal(&(*(m_uiAllLocalNode.begin())),m_uiAllLocalNode.size(),m_uiMaxDepth+1,m_uiMaxDepth+1,rootNode,ROOT_ROTATION,minMaxLocalNode[0],minMaxLocalNode[1]);
        //minMaxLocalNode[0] = m_uiAllLocalNode.front();
        //minMaxLocalNode[1] = m_uiAllLocalNode.back();

#ifdef DEBUG_E2N_MAPPING_SM
        treeNodesTovtk(m_uiAllLocalNode,m_uiActiveRank,"m_uiAllLocalNode");
#endif

        // 3. Gather the splitter max
        m_uiSplitterNodes=new ot::TreeNode[2*npes]; // both min and max
        std::vector<unsigned int > minMaxIDs;
        minMaxIDs.resize(2*npes);

        par::Mpi_Allgather(minMaxLocalNode,m_uiSplitterNodes,2,comm);
#ifdef DEBUG_E2N_MAPPING_SM
        std::vector<ot::TreeNode> splitterNodes;
        std::vector<ot::TreeNode> splitterElements;

        for(unsigned int p=0;p<npes;p++)
        {
            splitterNodes.push_back(m_uiSplitterNodes[p]);
            splitterElements.push_back(m_uiLocalSplitterElements[p]);
        }

        if(!rank) treeNodesTovtk(splitterNodes,rank,"splitterNodes");
        if(!rank) treeNodesTovtk(splitterElements,rank,"splitterElements");

        assert(seq::test::isUniqueAndSorted(splitterElements));
#endif

        m_uiScatterMapActualNodeSend.clear();
        m_uiSendNodeCount.resize(npes);//=new unsigned int [npes];
        m_uiRecvNodeCount.resize(npes);//=new unsigned int [npes];
        m_uiSendNodeOffset.resize(npes);//=new unsigned int [npes];
        m_uiRecvNodeOffset.resize(npes);//=new unsigned int [npes];

        std::set<unsigned int >* scatterMapNodeSet=new std::set<unsigned int > [npes]; // To keep track of the nodes to send to each processor.
        std::vector<ot::Key>* sendNodeOctants=new std::vector<ot::Key>[npes];
        std::vector<SearchKey> allocatedGhostNodes;

        //assert((m_uiNumPreGhostElements+m_uiNumPostGhostElements)==(m_uiRecvOctOffsetRound1[npes-1]+m_uiRecvOctCountRound1[npes-1]));


#ifdef DEBUG_E2N_MAPPING_SM
        std::vector<ot::TreeNode> customElements;
        std::vector<ot::TreeNode> sendNodes[npes];

        /*if(!m_uiActiveRank) customElements.push_back(m_uiAllElements[206]);
        if(!m_uiActiveRank) customElements.push_back(m_uiAllElements[57]);
        if(!m_uiActiveRank) treeNodesTovtk(customElements,m_uiActiveRank,"customElement");*/

        std::vector<ot::TreeNode> allocatedGNodes;

#endif

        // 3a. Compute send nodes based on the send elements that processor p has sent to other processors in the round 1 ghost communication.
        for(unsigned int p=0;p<npes;p++)
        {
            m_uiSendNodeCount[p]=0;
            scatterMapNodeSet[p]=std::set<unsigned int >();
            sendNodeOctants[p]=std::vector<ot::Key>();
        }


        std::vector<ot::TreeNode> neighbourElement;
        std::vector<ot::SearchKey> * gEleChained=new std::vector<ot::SearchKey>[npes];
        std::vector<ot::SearchKey>::iterator hintSK;
        std::vector<Key> ghostElementChained;
        std::vector<Key> tmpKeys;
        Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth);
        SearchKey rootSKey(0,0,0,0,m_uiDim,m_uiMaxDepth);
        unsigned int nodeFlag=0;
        unsigned int tmpEleLev=0;
        ot::TreeNode tmpElement;
        ot::SearchKey tmpSKey;
        ot::Key tmpKey;

        unsigned int* sendChainedGCount=new unsigned int [npes];
        unsigned int* recvChainedGCount=new unsigned int [npes];
        unsigned int* sendChainedGOffset=new unsigned int [npes];
        unsigned int* recvChainedGOffset=new unsigned int [npes];
        unsigned int* recvChainedKeyCount=new unsigned int [npes];
        unsigned int* recvChainedKeyOffset=new unsigned int [npes];

        std::vector<unsigned int > sendChainedGBuffer;
        std::set<unsigned  int>* sendChainedGhostIDSet=new std::set<unsigned int>[npes];


        // 4. Allocation of the node corresponding to pre & post level 1 ghost elements.

        nodeIndexVisited.clear();


        for(unsigned int ele=0;ele<m_uiGhostElementRound1Index.size();ele++)
        {
            for(unsigned int node=0;node<m_uiNpE;node++)
            {
                nodeIndex=m_uiE2NMapping_CG[m_uiGhostElementRound1Index[ele]*m_uiNpE+node];
                if(!(nodeIndex>=m_uiNodeLocalBegin && nodeIndex<m_uiNodeLocalEnd))
                {
                    assert(nodeIndex!=LOOK_UP_TABLE_DEFAULT);
                    setHintUint=nodeIndexVisited.emplace(nodeIndex);
                    if(setHintUint.second)
                    {
                        nodeIndex_DG=m_uiE2NMapping_DG[m_uiGhostElementRound1Index[ele]*m_uiNpE+node];
                        dg2eijk(nodeIndex_DG, ownerID, ii_x, jj_y, kk_z);

                        x = m_uiAllElements[ownerID].getX();
                        y = m_uiAllElements[ownerID].getY();
                        z = m_uiAllElements[ownerID].getZ();
                        sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
                        assert(sz%m_uiElementOrder==0);
                        tmpSKey=SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1);
                        tmpSKey.addOwner((*(setHintUint.first)));

                        allocatedGhostNodes.push_back(tmpSKey);
                        assert(ownerID<m_uiAllElements.size());

                        tmpElement=m_uiAllElements[ownerID];
                        assert(tmpElement.getLevel()>=m_uiAllElements[ownerID].getLevel());
                        nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);

                        if(nodeFlag==OCT_DIR_INTERNAL)assert(ownerID==m_uiGhostElementRound1Index[ele]);

                        tmpEleLev=tmpElement.getLevel();
                        tmpElement.setFlag((tmpEleLev) | (1u<<(nodeFlag+CHAINED_GHOST_OFFSET)));
                        assert(tmpElement.getFlag()>>(CHAINED_GHOST_OFFSET) & (1u<<nodeFlag));

                        hintSK=gEleChained[m_uiActiveRank].emplace(gEleChained[m_uiActiveRank].end(),SearchKey(tmpElement));
                        hintSK->addOwner(m_uiActiveRank);

                    }
                }

            }
        }


#ifdef DEBUG_E2N_MAPPING_SM
        treeNodesTovtk(allocatedGhostNodes,rank,"allocatedGhostNodes");
#endif


        ghostElementChained.clear();
        std::vector<SearchKey> tmpSearchKeyVec;
        unsigned int skip=1;
        for(unsigned int p=0;p<npes;p++)
        {

            SFC::seqSort::SFC_treeSort(&(*(gEleChained[p].begin())),gEleChained[p].size(),tmpSearchKeyVec,tmpSearchKeyVec,tmpSearchKeyVec,m_uiMaxDepth,m_uiMaxDepth,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);
            for(unsigned int e=0;e<(gEleChained[p].size());e++)
            {
                tmpElement=gEleChained[p][e];
                skip=1;
                while(((e+skip)<gEleChained[p].size()) && (gEleChained[p][e]==gEleChained[p][e+skip]))
                {
                    tmpElement.setFlag((tmpElement.getFlag()) | (gEleChained[p][e+skip].getFlag()));
                    skip++;
                }
                e+=(skip-1);


                tmpKey=ot::Key(tmpElement);
                tmpKey.addOwner(p);
                tmpKey.setFlag(tmpElement.getFlag());
                ghostElementChained.push_back(tmpKey);


            }


        }


#ifdef DEBUG_E2N_MAPPING_SM
        std::cout<<"m_uiActiveRank "<<m_uiActiveRank<<" gElementChained.size(): "<<ghostElementChained.size()<<std::endl;
        treeNodesTovtk(ghostElementChained,m_uiActiveRank,"ghostElementChained");
#endif

        std::vector<unsigned int> recvGhostChainedFlag;
        std::vector<Key> recvGhostChained;
        std::vector<unsigned int > recvChainedGBuffer;
        std::vector<Key> recvGhostChainedSearchKeys;
        std::vector<unsigned int> * ownerList;
        unsigned int result;

        unsigned int tmpCornerIndex;
        std::vector<unsigned int> tmpCornerIndexVec;

        std::vector<Key> missingNodes;
        std::vector<SearchKey> missingNodesSkey;


        while(ghostElementChained.size())
        {

            sendChainedGBuffer.clear();
            recvGhostChainedFlag.clear();
            recvGhostChained.clear();
            recvChainedGBuffer.clear();
            recvGhostChainedSearchKeys.clear();
            missingNodes.clear();
            missingNodesSkey.clear();

            for(unsigned int p=0;p<npes;p++)
            {
                gEleChained[p].clear();
                sendChainedGhostIDSet[p].clear();
                sendChainedGCount[p]=0;
                recvChainedKeyCount[p]=0;
            }


            unsigned int myX,myY,myZ,mySz;

            for(unsigned int e=0;e<ghostElementChained.size();e++)
            {

                myX=ghostElementChained[e].getX();
                myY=ghostElementChained[e].getY();
                myZ=ghostElementChained[e].getZ();
                mySz=1u<<(m_uiMaxDepth-ghostElementChained[e].getLevel());

                //assert(mySz%m_uiElementOrder ==0);

                nodeFlag=ghostElementChained[e].getFlag();
                nodeFlag=nodeFlag>>(CHAINED_GHOST_OFFSET);

                //assert(mySz%m_uiElementOrder==0);
                // Corner nodes.
                if((nodeFlag & (1u<<(OCT_DIR_LEFT_DOWN_BACK))))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+0*mySz),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                if((nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_BACK)))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+0*mySz),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                if((nodeFlag & (1u<<OCT_DIR_LEFT_UP_BACK)))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+1*mySz),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                if((nodeFlag & (1u<<OCT_DIR_RIGHT_UP_BACK)))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+1*mySz),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                if((nodeFlag & (1u<<OCT_DIR_LEFT_DOWN_FRONT)))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+0*mySz),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                if((nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_FRONT)))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+0*mySz),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                if((nodeFlag & (1u<<OCT_DIR_LEFT_UP_FRONT)))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+1*mySz),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                if((nodeFlag & (1u<<OCT_DIR_RIGHT_UP_FRONT)))
                {
                    hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+1*mySz),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    hintSK->addOwner(e);
                }

                //if(!m_uiActiveRank) std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" corner missing nodes: "<<missingNodesSet.size()<<std::endl;

                if(m_uiElementOrder>1)
                {

                    // internal Nodes;
                    if(nodeFlag & (1u<<OCT_DIR_INTERNAL))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+(mySz/2)),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                        /*if(m_uiActiveRank==0 && setHintKey.second)
                        {
                            std::cout<<" internal node added for element e: "<<e<<" val: "<<ghostElementChained[e]<<std::endl;
                            std::vector<ot::TreeNode> flagedElements;
                            flagedElements.push_back(ghostElementChained[e]);
                            flagedElements.push_back(*setHintKey.first);

                            treeNodesTovtk(flagedElements,e,"flagedKeys");
                        }*/

                    }

                    // internal edges.  can happen if only order is >1
                    if((nodeFlag & (1u<<OCT_DIR_LEFT_DOWN)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+0*mySz),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_LEFT_UP)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+1*mySz),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_LEFT_BACK)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+(mySz/2)),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_LEFT_FRONT)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+(mySz/2)),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+0*mySz),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_RIGHT_UP)))
                    {
                        //if(!m_uiActiveRank) std::cout<<" m_uiActiveRank: "<<m_uiActiveRank<<" RIGHT_UP MISSING NODE ADDED "<<std::endl;
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+1*mySz),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_RIGHT_BACK)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+(mySz/2)),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_RIGHT_FRONT)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+(mySz/2)),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_DOWN_BACK)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+0*mySz),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    }  if((nodeFlag & (1u<<OCT_DIR_DOWN_FRONT)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+0*mySz),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_UP_BACK)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+1*mySz),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_UP_FRONT)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+1*mySz),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    }


                    // internal faces.

                    if((nodeFlag & (1u<<OCT_DIR_LEFT)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+0*mySz),(myY+(mySz/2)),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    } if((nodeFlag & (1u<<OCT_DIR_RIGHT)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+1*mySz),(myY+(mySz/2)),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    }if((nodeFlag & (1u<<OCT_DIR_DOWN)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+0*mySz),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    }if((nodeFlag & (1u<<OCT_DIR_UP)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+1*mySz),(myZ+(mySz/2)),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);


                    }if((nodeFlag & (1u<<OCT_DIR_BACK)))
                    {

                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+(mySz/2)),(myZ+0*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);

                    }if((nodeFlag & (1u<<OCT_DIR_FRONT)))
                    {
                        hintSK= missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey((myX+(mySz/2)),(myY+(mySz/2)),(myZ+1*mySz),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        hintSK->addOwner(e);
                    }

                    //if(!m_uiActiveRank) std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" internal missing nodes: "<<missingNodesSet.size()<<" missing Node: "<<*setHintKey.first<<std::endl;


                }




            }


#ifdef DEBUG_E2N_MAPPING_SM
            missingNodes.clear();
            missingNodes.insert(missingNodes.end(),missingNodesSet.begin(),missingNodesSet.end());
            std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" missing node size: "<<missingNodes.size()<<std::endl;
            treeNodesTovtk(missingNodes,m_uiActiveRank,"missingNodes");
            treeNodesTovtk(nonLocalNodes,rank,"nonLocal");
            missingNodes.clear();
#endif



            for(unsigned int p=0;p<2*npes;p++) {
                missingNodesSkey.emplace(missingNodesSkey.end(),SearchKey(m_uiSplitterNodes[p]));
            }

            missingNodes.clear();
            SFC::seqSort::SFC_treeSort(&(*(missingNodesSkey.begin())),missingNodesSkey.size(),tmpSearchKeyVec,tmpSearchKeyVec,tmpSearchKeyVec,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

            for(unsigned int e=0;e<(missingNodesSkey.size());e++)
            {
                skip=1;
                tmpKey=Key(missingNodesSkey[e].getX(),missingNodesSkey[e].getY(),missingNodesSkey[e].getZ(),missingNodesSkey[e].getLevel(),m_uiDim,m_uiMaxDepth+1);
                if(missingNodesSkey[e].getOwner()>=0)tmpKey.addOwner(missingNodesSkey[e].getOwner());
                while(((e+skip)<missingNodesSkey.size()) && (missingNodesSkey[e]==missingNodesSkey[e+skip]))
                {
                    if(missingNodesSkey[e+skip].getOwner()>=0)tmpKey.addOwner(missingNodesSkey[e+skip].getOwner());
                    skip++;
                }
                missingNodes.push_back(tmpKey);
                e+=(skip-1);
            }



            /*missingNodes.insert(missingNodes.end(),missingNodesSet.begin(),missingNodesSet.end());

            //NOTE: We need to use m_uiMaxDepth+1 to sort because this generates the all the possible nodes. Hence it contains 1u<<m_uiMaxDepth for (x y z) values.
            SFC::seqSort::SFC_treeSort(&(*(missingNodes.begin())),missingNodes.size(),tmpKeys,tmpKeys,tmpKeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootKey,ROOT_ROTATION,1,TS_SORT_ONLY);*/

            std::vector<ot::Key> splitterNode_keys;
            splitterNode_keys.resize(2*npes);
            for(unsigned int p=0;p<2*npes;p++)
                splitterNode_keys[p]=ot::Key(m_uiSplitterNodes[p]);


            m_uiMaxDepth++;
            searchKeys(splitterNode_keys,missingNodes);
            m_uiMaxDepth--;


            for(unsigned int p=0;p<2*npes;p++)
            {
                assert(splitterNode_keys[p].getFlag() & OCT_FOUND);
                minMaxIDs[p]=splitterNode_keys[p].getSearchResult();
                /*minMaxIDs[p]=(std::find(missingNodes.begin(),missingNodes.end(),m_uiSplitterNodes[p])-missingNodes.begin());
                if(minMaxIDs[p]!=splitterNode_keys[p].getSearchResult())
                    std::cout<<" minMax ID: "<<minMaxIDs[p]<<" sfcTSearch: "<<splitterNode_keys[p].getSearchResult()<<std::endl;*/
                assert(minMaxIDs[p]<missingNodes.size());
            }




            Key * missingNodesPtr=&(*(missingNodes.begin()));
            unsigned int sBegin,sEnd;
            for(unsigned int p=0;p<npes;p++)
            {
                if(p==m_uiActiveRank) continue;
                sBegin=minMaxIDs[2*p];
                sEnd=minMaxIDs[2*p+1]+1;
                for(unsigned int e=sBegin;e<sEnd;e++)
                {
                    for (unsigned int w = 0; w < missingNodesPtr[e].getOwnerList()->size(); w++)
                    {
                        setHintUint=sendChainedGhostIDSet[p].emplace((*(missingNodesPtr[e].getOwnerList()))[w]);
                    }
                }

            }

            // std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" splitterNodeCount: "<<splitterNodeCount<<std::endl;
            for(unsigned int p=0;p<npes;p++)
            {
                for(auto it=sendChainedGhostIDSet[p].begin();it!=sendChainedGhostIDSet[p].end();++it)
                {

                    sendChainedGBuffer.push_back(ghostElementChained[*it].getX());
                    sendChainedGBuffer.push_back(ghostElementChained[*it].getY());
                    sendChainedGBuffer.push_back(ghostElementChained[*it].getZ());
                    sendChainedGBuffer.push_back(ghostElementChained[*it].getFlag());
                    sendChainedGBuffer.push_back(ghostElementChained[*it].getOwnerList()->size());
                    sendChainedGCount[p]+=5+ghostElementChained[*it].getOwnerList()->size();
                    for (unsigned int w = 0; w < ghostElementChained[*it].getOwnerList()->size(); w++)
                        sendChainedGBuffer.push_back((*(ghostElementChained[*it].getOwnerList()))[w]);

                }

            }

            par::Mpi_Alltoall(sendChainedGCount,recvChainedGCount,1,comm);


            sendChainedGOffset[0]=0;
            recvChainedGOffset[0]=0;

            omp_par::scan(sendChainedGCount,sendChainedGOffset,npes);
            omp_par::scan(recvChainedGCount,recvChainedGOffset,npes);

            //std::cout<<" m_uiActiveRank: "<<m_uiActiveRank<<" sendBuf ID size: "<<sendChainedGBuffer.size()<<" sendCount Size: "<<(sendChainedGOffset[npes-1]+sendChainedGCount[npes-1])<<std::endl;
            assert(sendChainedGBuffer.size()==(sendChainedGOffset[npes-1]+sendChainedGCount[npes-1]));


            recvChainedGBuffer.resize(recvChainedGOffset[npes-1]+recvChainedGCount[npes-1]);
            par::Mpi_Alltoallv(&(*(sendChainedGBuffer.begin())),(int *) sendChainedGCount,(int *) sendChainedGOffset,&(*(recvChainedGBuffer.begin())),(int *) recvChainedGCount,(int *) recvChainedGOffset,comm);



            unsigned int recvGindex=0;
            unsigned int pCount=0;
            unsigned int recvKeyCount=0;

            for(unsigned int p=0;p<npes;p++)
            {   recvKeyCount=0;
                while(recvGindex<(recvChainedGOffset[p]+recvChainedGCount[p]))
                {
                    tmpKey=Key(recvChainedGBuffer[recvGindex],recvChainedGBuffer[recvGindex+1],recvChainedGBuffer[recvGindex+2],(recvChainedGBuffer[recvGindex+3] & ot::TreeNode::MAX_LEVEL),m_uiDim,m_uiMaxDepth);
                    recvGhostChainedFlag.push_back(recvChainedGBuffer[recvGindex+3]);
                    tmpKey.getOwnerList()->resize(recvChainedGBuffer[recvGindex+4]);
                    tmpKey.getOwnerList()->assign((recvChainedGBuffer.begin()+(recvGindex+5)),(recvChainedGBuffer.begin()+(recvGindex+5)+recvChainedGBuffer[recvGindex+4]));
                    recvGindex=recvGindex+5+recvChainedGBuffer[recvGindex+4];
                    recvGhostChained.push_back(tmpKey);
                    recvKeyCount++;

                }

                recvChainedKeyCount[p]=recvKeyCount;

            }

            recvChainedKeyOffset[0]=0;
            omp_par::scan(recvChainedKeyCount,recvChainedKeyOffset,npes);

#ifdef DEBUG_E2N_MAPPING_SM
            treeNodesTovtk(recvGhostChained,m_uiActiveRank,"recvGChained");
#endif

            recvGhostChainedSearchKeys.resize(recvGhostChained.size());
            for(unsigned int ele=0;ele<recvGhostChained.size();ele++)
            {
                recvGhostChainedSearchKeys[ele]=Key(recvGhostChained[ele].getX(),recvGhostChained[ele].getY(),recvGhostChained[ele].getZ(),recvGhostChained[ele].getLevel(),m_uiDim,m_uiMaxDepth);
                recvGhostChainedSearchKeys[ele].addOwner(ele);
                recvGhostChained[ele].setSearchResult(LOOK_UP_TABLE_DEFAULT);

            }


            /*std::cout<<"m_uiActiveRank "<<m_uiActiveRank<<" recvKeySize: "<<recvGhostChained.size()<<std::endl;*/
            /* if(m_uiActiveRank==12)
             for(unsigned int p=0;p<npes;p++)
             {   //std::cout<<"recvKeyGCount p : "<<p<<" : "<< (recvChainedGOffset[p]+recvChainedGCount[p])<<std::endl;
                 std::cout<<"recvKeyCount p: "<<p<<" : "<<recvChainedKeyCount[p]<<std::endl;
             }*/

            SFC::seqSearch::SFC_treeSearch(&(*(recvGhostChainedSearchKeys.begin())),&(*(m_uiAllElements.begin())),0,recvGhostChainedSearchKeys.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);
            for(unsigned int ele=0;ele<recvGhostChained.size();ele++)
            {
                assert(recvGhostChained[ele].getOwnerList()->size()==1);
                if(recvGhostChainedSearchKeys[ele].getFlag() & OCT_FOUND)
                {
                    recvGhostChained[(*(recvGhostChainedSearchKeys[ele].getOwnerList()))[0]].setSearchResult(recvGhostChainedSearchKeys[ele].getSearchResult());
                }

            }


            ghostElementChained.clear();
            std::vector<unsigned int > edgeInternalIndex;
            std::vector<unsigned int > faceInternalIndex;
            std::vector<unsigned int > elementInternalIndex;



            for(unsigned int ele=0;ele<recvGhostChained.size();ele++)
            {

                if(recvGhostChained[ele].getSearchResult()!=LOOK_UP_TABLE_DEFAULT)
                { // implies that current chainedGhost is one of my local node.

                    tmpCornerIndexVec.clear();
                    ownerList=recvGhostChained[ele].getOwnerList();
                    result=recvGhostChained[ele].getSearchResult();
                    nodeFlag=recvGhostChainedFlag[ele];
                    nodeFlag=nodeFlag>>(CHAINED_GHOST_OFFSET);

                    /*if(*//*m_uiActiveRank==12*//*(m_uiActiveRank==12 && ele==615) || ( m_uiActiveRank==12 && (recvGhostChained[ele]==recvGhostChained[615])))
                {

                    std::cout<<"ele: "<<ele<<" Key: "<<recvGhostChained[ele]<<" found: "<<m_uiAllElements[result]<<" nodeFlag: "<<recvGhostChainedFlag[ele]<<std::endl;
                    for(unsigned int node=0;node<m_uiNpE;node++)
                    {
                        nodeIndex=m_uiE2NMapping_CG[result*m_uiNpE+node];
                        if(nodeIndex>=m_uiNodeLocalBegin && nodeIndex<m_uiNodeLocalEnd)
                        {
                            std::cout<<"E2N: "<<ele<<" found: "<<result<<" node: "<<"node: "<<node;
                        }
                        std::cout<<std::endl;


                    }
                    std::vector<ot::TreeNode> cusElement;
                    cusElement.push_back(recvGhostChained[ele]);
                    treeNodesTovtk(cusElement,ele,"cusElement");
                    cusElement.clear();


                }
*/

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN_BACK))
                    {
                        cornerNodeIndex(result,0,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_BACK))
                    {
                        cornerNodeIndex(result,1,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_UP_BACK))
                    {
                        cornerNodeIndex(result,2,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP_BACK))
                    {
                        cornerNodeIndex(result,3,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN_FRONT))
                    {
                        cornerNodeIndex(result,4,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_FRONT))
                    {
                        cornerNodeIndex(result,5,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_UP_FRONT))
                    {
                        cornerNodeIndex(result,6,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP_FRONT))
                    {
                        cornerNodeIndex(result,7,tmpCornerIndex);
                        tmpCornerIndexVec.push_back(tmpCornerIndex);
                    }


                    if(m_uiElementOrder>1)
                    {



                        if(nodeFlag & (1u<<OCT_DIR_INTERNAL))
                        {
                            elementNodeIndex(result,elementInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),elementInternalIndex.begin(),elementInternalIndex.end());
                        }


                        if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN))
                        {

                            edgeNodeIndex(result,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }

                        if(nodeFlag & (1u<<OCT_DIR_LEFT_UP))
                        {
                            edgeNodeIndex(result,OCT_DIR_LEFT,OCT_DIR_UP,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }

                        if(nodeFlag & (1u<<OCT_DIR_LEFT_BACK))
                        {
                            edgeNodeIndex(result,OCT_DIR_LEFT,OCT_DIR_BACK,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }

                        if(nodeFlag & (1u<<OCT_DIR_LEFT_FRONT))
                        {
                            edgeNodeIndex(result,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }


                        if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN))
                        {
                            edgeNodeIndex(result,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }

                        if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP))
                        {
                            edgeNodeIndex(result,OCT_DIR_RIGHT,OCT_DIR_UP,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }

                        if(nodeFlag & (1u<<OCT_DIR_RIGHT_BACK))
                        {
                            edgeNodeIndex(result,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }

                        if(nodeFlag & (1u<<OCT_DIR_RIGHT_FRONT))
                        {
                            edgeNodeIndex(result,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());

                        }


                        if(nodeFlag & (1u<<OCT_DIR_DOWN_BACK))
                        {
                            edgeNodeIndex(result,OCT_DIR_DOWN,OCT_DIR_BACK,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());
                        }


                        if(nodeFlag & (1u<<OCT_DIR_DOWN_FRONT))
                        {
                            edgeNodeIndex(result,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());
                        }


                        if(nodeFlag & (1u<<OCT_DIR_UP_BACK))
                        {
                            edgeNodeIndex(result,OCT_DIR_UP,OCT_DIR_BACK,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());
                        }


                        if(nodeFlag & (1u<<OCT_DIR_UP_FRONT))
                        {
                            edgeNodeIndex(result,OCT_DIR_UP,OCT_DIR_FRONT,edgeInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),edgeInternalIndex.begin(),edgeInternalIndex.end());
                        }


                        if(nodeFlag & (1u<<OCT_DIR_LEFT))
                        {
                            faceNodesIndex(result,OCT_DIR_LEFT,faceInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),faceInternalIndex.begin(),faceInternalIndex.end());
                        }


                        if(nodeFlag & (1u<<OCT_DIR_RIGHT))
                        {
                            faceNodesIndex(result,OCT_DIR_RIGHT,faceInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),faceInternalIndex.begin(),faceInternalIndex.end());
                        }
                        if(nodeFlag & (1u<<OCT_DIR_DOWN))
                        {
                            faceNodesIndex(result,OCT_DIR_DOWN,faceInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),faceInternalIndex.begin(),faceInternalIndex.end());
                        }

                        if(nodeFlag & (1u<<OCT_DIR_UP))
                        {
                            faceNodesIndex(result,OCT_DIR_UP,faceInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),faceInternalIndex.begin(),faceInternalIndex.end());
                        }

                        if(nodeFlag & (1u<<OCT_DIR_BACK))
                        {
                            faceNodesIndex(result,OCT_DIR_BACK,faceInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),faceInternalIndex.begin(),faceInternalIndex.end());
                        }

                        if(nodeFlag & (1u<<OCT_DIR_FRONT))
                        {
                            faceNodesIndex(result,OCT_DIR_FRONT,faceInternalIndex,true);
                            tmpCornerIndexVec.insert(tmpCornerIndexVec.end(),faceInternalIndex.begin(),faceInternalIndex.end());
                        }


                    }





                    for(unsigned int cornerNodeIndex=0;cornerNodeIndex<tmpCornerIndexVec.size();cornerNodeIndex++)
                    {

                        nodeIndex=m_uiE2NMapping_CG[tmpCornerIndexVec[cornerNodeIndex]];
                        nodeIndex_DG=m_uiE2NMapping_DG[tmpCornerIndexVec[cornerNodeIndex]];
                        dg2eijk(nodeIndex_DG,ownerID,ii_x,jj_y,kk_z);
                        //if(!m_uiActiveRank && getDIROfANode(ii_x,jj_y,kk_z)==OCT_DIR_LEFT_DOWN) std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" OwnerID: "<<ownerID<<" ii_x: "<<ii_x<<" jj_y: "<<jj_y<<" kk_z: "<<kk_z<<std::endl;
                        nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);

                        for(unsigned int w=0;w<ownerList->size();w++)
                        {
                            if((nodeIndex>=m_uiNodeLocalBegin) && (nodeIndex<m_uiNodeLocalEnd))
                            {

                                assert((*ownerList)[w]!=m_uiActiveRank);
                                setHintUint=scatterMapNodeSet[(*ownerList)[w]].emplace(nodeIndex);
                                if(setHintUint.second)
                                {
                                    dg2eijk(nodeIndex_DG,ownerID,ii_x,jj_y,kk_z);
                                    x=m_uiAllElements[ownerID].getX();
                                    y=m_uiAllElements[ownerID].getY();
                                    z=m_uiAllElements[ownerID].getZ();
                                    sz=1u<<(m_uiMaxDepth-m_uiAllElements[ownerID].getLevel());
                                    assert(sz%m_uiElementOrder==0);
                                    tmpKey=ot::Key((x+ii_x*sz/m_uiElementOrder),(y+jj_y*sz/m_uiElementOrder),(z+kk_z*sz/m_uiElementOrder),m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1);
                                    tmpKey.addOwner((*setHintUint.first));
                                    assert((*setHintUint.first)==nodeIndex);
                                    sendNodeOctants[(*ownerList)[w]].push_back(tmpKey);
#ifdef DEBUG_E2N_MAPPING_SM
                                    x=m_uiAllElements[ownerID].getX();
                                    y=m_uiAllElements[ownerID].getY();
                                    z=m_uiAllElements[ownerID].getZ();
                                    sz=1u<<(m_uiMaxDepth-m_uiAllElements[ownerID].getLevel());
                                    sendNodes[(*ownerList)[w]].push_back(ot::TreeNode((x+ii_x*sz/m_uiElementOrder),(y+jj_y*sz/m_uiElementOrder),(z+kk_z*sz/m_uiElementOrder),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
                                    //if((*ownerList)[w]==2) std::cout<<" m_uiActiveRank: "<<m_uiActiveRank<<" SendNode: "<<sendNodes[(*ownerList)[w]].back()<<" to rank 2"<<std::endl;
#endif
                                    m_uiSendNodeCount[(*ownerList)[w]]++;
                                    //std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" R2 SendNodes Executed. "<<std::endl;

                                }

                            }else
                            {


                                assert(!(nodeIndex>=m_uiNodeLocalBegin && nodeIndex<m_uiNodeLocalEnd));
                                assert(ownerID<m_uiAllElements.size());
                                tmpElement=m_uiAllElements[ownerID];
                                assert(tmpElement.getLevel()>=m_uiAllElements[ownerID].getLevel());
                                nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);
                                //if(!m_uiActiveRank && nodeFlag==(OCT_DIR_LEFT_DOWN)) std::cout<<RED<<" nodeFlag: "<<nodeFlag<<NRM<<std::endl;
                                //assert(nodeFlag!=OCT_DIR_INTERNAL);
                                tmpEleLev=tmpElement.getLevel();
                                tmpElement.setFlag((tmpEleLev) | (1u<<(nodeFlag+CHAINED_GHOST_OFFSET)));
                                assert(tmpElement.getFlag()>>(CHAINED_GHOST_OFFSET) & (1u<<nodeFlag));

                                if((ownerID>=m_uiElementLocalBegin && ownerID<m_uiElementLocalEnd))
                                {

                                    hintSK=gEleChained[(*ownerList)[w]].emplace(gEleChained[(*ownerList)[w]].end(),SearchKey(tmpElement));
                                    hintSK->addOwner((*ownerList)[w]);


                                } else
                                {
                                    assert(ownerID!=LOOK_UP_TABLE_DEFAULT);


                                }



                            }


                        }

                    }

                }

            }


            unsigned int skip=1;
            for(unsigned int p=0;p<npes;p++)
            {


                SFC::seqSort::SFC_treeSort(&(*(gEleChained[p].begin())),gEleChained[p].size(),tmpSearchKeyVec,tmpSearchKeyVec,tmpSearchKeyVec,m_uiMaxDepth,m_uiMaxDepth,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);
                for(unsigned int e=0;e<(gEleChained[p].size());e++)
                {
                    tmpElement=gEleChained[p][e];
                    skip=1;
                    while(((e+skip)<gEleChained[p].size()) && (gEleChained[p][e]==gEleChained[p][e+skip]))
                    {
                        tmpElement.setFlag((tmpElement.getFlag())|(gEleChained[p][e+skip].getFlag()));
                        skip++;
                    }
                    e+=(skip-1);

                    tmpKey=ot::Key(tmpElement);
                    tmpKey.addOwner(p);
                    tmpKey.setFlag(tmpElement.getFlag());
                    ghostElementChained.push_back(tmpKey);

                }



            }


        }







#ifdef DEBUG_E2N_MAPPING_SM


        for(unsigned int p=0;p<npes;p++)
        {
            tmpNode.clear();
            SFC::seqSort::SFC_treeSort(&(*(sendNodes[p].begin())),sendNodes[p].size(),tmpNode,tmpNode,tmpNode,m_uiMaxDepth,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_SORT_ONLY);
            assert(seq::test::isUniqueAndSorted(sendNodes[p]));
        }


        unsigned int * sendNodeCount=new unsigned int [npes];
        unsigned int * recvNodeCount=new unsigned int [npes];
        unsigned int * sendNodeOffset=new unsigned int [npes];
        unsigned int * recvNodeOffset=new unsigned int [npes];

        std::vector<ot::TreeNode> sendNodeBuffer;
        std::vector<ot::TreeNode> recvNodeBuffer;


        for(unsigned int p=0;p<npes;p++) {
            sendNodeCount[p] = scatterMapNodeSet[p].size();
            sendNodeBuffer.insert(sendNodeBuffer.end(),sendNodes[p].begin(),sendNodes[p].end());
        }

        par::Mpi_Alltoall(sendNodeCount,recvNodeCount,1,comm);

        sendNodeOffset[0]=0;
        recvNodeOffset[0]=0;

        omp_par::scan(sendNodeCount,sendNodeOffset,npes);
        omp_par::scan(recvNodeCount,recvNodeOffset,npes);

        recvNodeBuffer.resize(recvNodeCount[npes-1]+recvNodeOffset[npes-1]);
        par::Mpi_Alltoallv(&(*(sendNodeBuffer.begin())),(int *) sendNodeCount,(int *) sendNodeOffset, &(*(recvNodeBuffer.begin())),(int *) recvNodeCount,(int *) recvNodeOffset,comm);

        tmpNode.clear();
        SFC::seqSort::SFC_treeSort(&(*(allocatedGNodes.begin())),allocatedGNodes.size(),tmpNode,tmpNode,tmpNode,m_uiMaxDepth,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_SORT_ONLY);
        assert(seq::test::isUniqueAndSorted(allocatedGNodes));


       /* if(m_uiActiveRank==2)
        {
            for(unsigned int ele=0;ele<recvNodeBuffer.size();ele++)
            {
                std::cout<<"recvNode: "<<recvNodeBuffer[ele]<<std::endl;
            }
        }*/

        tmpNode.clear();
        SFC::seqSort::SFC_treeSort(&(*(recvNodeBuffer.begin())),recvNodeBuffer.size(),tmpNode,tmpNode,tmpNode,m_uiMaxDepth,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_SORT_ONLY);
        assert(seq::test::isUniqueAndSorted(recvNodeBuffer));


        if(recvNodeBuffer.size()!=allocatedGNodes.size())
        {
            std::vector<Key> recvNodeKeys;
            recvNodeKeys.resize(recvNodeBuffer.size());
            std::vector<ot::TreeNode> missmatchedNodes;


            for(unsigned int e=0;e<recvNodeBuffer.size();e++)
            {
               unsigned int findIndex=std::find(allocatedGNodes.begin(),allocatedGNodes.end(),recvNodeBuffer[e])-allocatedGNodes.begin();
               if(findIndex >=allocatedGNodes.size())
                   missmatchedNodes.push_back(recvNodeBuffer[e]);
            }

            treeNodesTovtk(allocatedGNodes,rank,"allGNode");
            treeNodesTovtk(recvNodeBuffer,rank,"recvNodeBuffer");


            treeNodesTovtk(missmatchedNodes,rank,"missmatchedNodes");

        }







        delete [] sendNodeCount;
        delete [] recvNodeCount;
        delete [] sendNodeOffset;
        delete [] recvNodeOffset;




        for(unsigned int p=0;p<npes;p++) {
            char filename[256];
            sprintf(filename,"sendNode_R2%d",p);
            treeNodesTovtk(sendNodes[p],rank, filename);
            sendNodes[p].clear();
        }
#endif



        delete [] sendChainedGCount;
        delete [] recvChainedGCount;
        delete [] sendChainedGOffset;
        delete [] recvChainedGOffset;
        delete [] recvChainedKeyCount;
        delete [] recvChainedKeyOffset;
        delete [] gEleChained;





#ifdef DEBUG_E2N_MAPPING
        treeNodesTovtk(neighbourElement,m_uiActiveRank,"neighBourElement");
        std::cout<<"======== m_uiActiveRank: "<<m_uiActiveRank<<" : "<<allocatedGNodes.size()<<std::endl;
        treeNodesTovtk(ownerElement1,rank,"owner1Elements");
        treeNodesTovtk(ownerElement2,rank,"owner2Elements");
        treeNodesTovtk(nonLocalNodes,rank,"nonLocal");
        treeNodesTovtk(LocalNodes,rank,"Local");
#endif

        //prepare send recv scatter maps.
        m_uiScatterMapActualNodeSend.clear();
        m_uiScatterMapActualNodeRecv.clear();


        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),(&(*(m_uiRecvNodeOffset.begin()))),npes);

        std::vector<ot::TreeNode> sendNodeBuffer;
        std::vector<ot::TreeNode> recvNodeBuffer;

        recvNodeBuffer.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);
        //sendNodeBuffer.resize(m_uiSendNodeOffset[npes-1]+m_uiSendNodeCount[npes-1]);
        for(unsigned int p=0;p<npes;p++) {

            /*m_uiScatterMapActualNodeSend.insert(m_uiScatterMapActualNodeSend.end(), scatterMapNodeSet[p].begin(),
                                                scatterMapNodeSet[p].end());*/
            for(unsigned int k=0;k<sendNodeOctants[p].size();k++) {
                assert((sendNodeOctants[p][k].getOwnerList()->size()==1));
                m_uiScatterMapActualNodeSend.push_back((sendNodeOctants[p][k].getOwnerList()->front()));
                sendNodeBuffer.push_back(sendNodeOctants[p][k]);
            }
            assert(m_uiSendNodeCount[p]==scatterMapNodeSet[p].size());
            assert(scatterMapNodeSet[p].size()==sendNodeOctants[p].size());
            scatterMapNodeSet[p].clear();
            sendNodeOctants[p].clear();

        }

        delete [] scatterMapNodeSet;
        delete [] sendChainedGhostIDSet;
        delete [] sendNodeOctants;

        par::Mpi_Alltoallv(&(*(sendNodeBuffer.begin())),(int *) (&(*(m_uiSendNodeCount.begin()))), (int *) (&(*(m_uiSendNodeOffset.begin()))),(&(*(recvNodeBuffer.begin()))),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *) (&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);

        std::vector<ot::Key> recvNodeKeys;
        recvNodeKeys.resize(recvNodeBuffer.size());

        for(unsigned int k=0;k<recvNodeBuffer.size();k++)
        {
            tmpKey=ot::Key(recvNodeBuffer[k]);
            tmpKey.addOwner(k);
            recvNodeKeys[k]=tmpKey;
        }

#ifdef DEBUG_MESH_GENERATION
        treeNodesTovtk(allocatedGhostNodes,m_uiActiveRank,"allocatedNodes");
        treeNodesTovtk(recvNodeBuffer,m_uiActiveRank,"recvNodes");
#endif
        if((m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1])!=allocatedGhostNodes.size())
        {
            std::cout<<RED<<"[Error]"<<" m_uiActiveRank: "<<m_uiActiveRank<<" Total Ghost Elements allocated: "<<allocatedGhostNodes.size()<<" Number of elements will get recieved: "<<(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1])<<NRM<<std::endl;
            assert(false);
        }


       /* std::sort(allocatedGhostNodes.begin(),allocatedGhostNodes.end(),OctreeComp<ot::Key>());
        std::sort(recvNodeKeys.begin(),recvNodeKeys.end(),OctreeComp<ot::Key>());*/


        SFC::seqSort::SFC_treeSort(&(*(allocatedGhostNodes.begin())),allocatedGhostNodes.size(),tmpSearchKeyVec,tmpSearchKeyVec,tmpSearchKeyVec,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);
        SFC::seqSort::SFC_treeSort(&(*(recvNodeKeys.begin())),recvNodeKeys.size(),tmpKeys,tmpKeys,tmpKeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootKey,ROOT_ROTATION,1,TS_SORT_ONLY);



        m_uiScatterMapActualNodeRecv.clear();
        m_uiScatterMapActualNodeRecv.resize((m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]));
        assert(allocatedGhostNodes.size()==recvNodeKeys.size());
        for(unsigned int k=0;k<allocatedGhostNodes.size();k++)
        {
            if(allocatedGhostNodes[k]!=recvNodeKeys[k]) {
                std::cout << RED << "[ERROR]: " << "m_uiActiveRank : " << m_uiActiveRank << " allocated [" << k << "]: "<< allocatedGhostNodes[k] << " recieved[" << k << "]: " << recvNodeKeys[k] << std::endl;
                assert(false);
            }

            //m_uiScatterMapActualNodeRecv[(*(recvNodeKeys[k].getOwnerList()))[0]]=(*(allocatedGhostNodes[k].getOwnerList()))[0];
            //assert(allocatedGhostNodes[k].getOwnerList()->size()==1);
            m_uiScatterMapActualNodeRecv[(*(recvNodeKeys[k].getOwnerList()))[0]]=allocatedGhostNodes[k].getOwner();

        }





#ifdef DEBUG_E2N_MAPPING
        MPI_Barrier(comm);
        if(!rank)
            for(unsigned int p=0;p<npes;p++)
                std::cout<<"rank: "<<rank<<" recv nodes from : ["<<p<<"]: begin:  "<<m_uiRecvNodeOffset[p]<<" end: "<<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p])<<std::endl;


        MPI_Barrier(comm);
        if(rank==1)
            for(unsigned int p=0;p<npes;p++)
                std::cout<<"rank: "<<rank<<" recv nodes from : ["<<p<<"]: begin:  "<<m_uiRecvNodeOffset[p]<<" end: "<<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p])<<std::endl;


        MPI_Barrier(comm);
#endif






    }


    void Mesh::computeNodalScatterMap(MPI_Comm comm)
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(npes<=1) return; // nothing to do in the sequential case. (No scatter map required.)

        unsigned int x,y,z,sz; // x y z and size of an octant.
        unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
        unsigned int nodeIndex;

        std::vector<SearchKey> allocatedNodes;
        std::vector<SearchKey> localNodes;
        std::vector<SearchKey>::iterator it;
        std::vector<SearchKey> tmpSkeys;
        std::vector<Key> tmpKeys;

        SearchKey rootSKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);
        Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);

        ot::TreeNode minMaxLocalNode[2];
        m_uiSplitterNodes=new ot::TreeNode[2*npes]; // both min and max (nodal)
        std::vector<unsigned int > minMaxIDs;
        minMaxIDs.resize(2*npes);


        for(unsigned int e=m_uiNodeLocalBegin;e<m_uiNodeLocalEnd;e++)
        {
            dg2eijk(m_uiCG2DG[e],ownerID,ii_x,jj_y,kk_z);
            x = m_uiAllElements[ownerID].getX();
            y = m_uiAllElements[ownerID].getY();
            z = m_uiAllElements[ownerID].getZ();
            sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
            assert(sz%m_uiElementOrder==0);
            it=localNodes.emplace(localNodes.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
            it->addOwner(e);
        }


        SFC::seqSort::SFC_treeSort(&(*(localNodes.begin())),localNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(localNodes));
        m_uiMaxDepth--;

        minMaxLocalNode[0]=localNodes.front();
        minMaxLocalNode[1]=localNodes.back();

        par::Mpi_Allgather(minMaxLocalNode,m_uiSplitterNodes,2,comm);


        std::vector<bool > g1Visited;
        g1Visited.resize(m_uiCG2DG.size(),false);

        for(unsigned int ele=0;ele<m_uiGhostElementRound1Index.size();ele++) {
            for (unsigned int node = 0; node < m_uiNpE; node++) {
                nodeIndex = m_uiE2NMapping_CG[m_uiGhostElementRound1Index[ele] * m_uiNpE + node];
                if ((!(nodeIndex >= m_uiNodeLocalBegin && nodeIndex < m_uiNodeLocalEnd)) && (!g1Visited[nodeIndex])) {
                    assert(nodeIndex<g1Visited.size());
                    dg2eijk(m_uiCG2DG[nodeIndex],ownerID,ii_x,jj_y,kk_z);
                    x = m_uiAllElements[ownerID].getX();
                    y = m_uiAllElements[ownerID].getY();
                    z = m_uiAllElements[ownerID].getZ();
                    sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
                    it=allocatedNodes.emplace(allocatedNodes.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                    it->addOwner(nodeIndex);
                    g1Visited[nodeIndex]= true;
                }
            }
        }

        // number of allocated nodes without splitters. (By construction this should be unique)
        unsigned int allocatedNodeSz=allocatedNodes.size();

        for(unsigned int p=0;p<2*npes;p++)
            allocatedNodes.emplace(allocatedNodes.end(),m_uiSplitterNodes[p]);

        SFC::seqSort::SFC_treeSort(&(*(allocatedNodes.begin())),allocatedNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


        tmpSkeys.clear();
        SearchKey tmpSkey;
        unsigned int skip;
        for(unsigned int e=0;e<(allocatedNodes.size());e++)
        {
            skip=1;
            tmpSkey=allocatedNodes[e];
            while(((e+skip)<allocatedNodes.size()) && (allocatedNodes[e]==allocatedNodes[e+skip]))
            {
                if(allocatedNodes[e+skip].getOwner()>=0)tmpSkey.addOwner(allocatedNodes[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(allocatedNodes,tmpSkeys);
        tmpSkeys.clear();

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(allocatedNodes));
        m_uiMaxDepth--;


        std::vector<ot::Key> splitterNode_keys;
        splitterNode_keys.resize(2*npes);
        for(unsigned int p=0;p<2*npes;p++)
            splitterNode_keys[p]=ot::Key(m_uiSplitterNodes[p]);


        m_uiMaxDepth++;
        searchKeys(splitterNode_keys,allocatedNodes);
        m_uiMaxDepth--;

        for(unsigned int p=0;p<2*npes;p++)
        {
            assert(splitterNode_keys[p].getFlag() & OCT_FOUND);
            minMaxIDs[p]=splitterNode_keys[p].getSearchResult();
            assert(minMaxIDs[p]<allocatedNodes.size());
        }


        m_uiScatterMapActualNodeSend.clear();
        m_uiScatterMapActualNodeRecv.clear();

        m_uiSendNodeCount.resize(npes);
        m_uiRecvNodeCount.resize(npes);
        m_uiSendNodeOffset.resize(npes);
        m_uiRecvNodeOffset.resize(npes);

        for(unsigned int p=0;p<npes;p++)
            m_uiSendNodeCount[p]=0;


        std::vector<ot::TreeNode> sendNodes;
        std::vector<ot::TreeNode> recvNodes;


        for(unsigned int p=0;p<npes;p++)
        {
            if(p==rank) continue;
            for(unsigned int e=minMaxIDs[2*p];e<(minMaxIDs[2*p+1]+1);e++)
            {
                if(allocatedNodes[e].getOwner()>=0){
                    sendNodes.push_back(allocatedNodes[e]);
                    m_uiSendNodeCount[p]++;
                }
            }
        }

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        assert(sendNodes.size()==(m_uiSendNodeOffset[npes-1]+m_uiSendNodeCount[npes-1]));
        recvNodes.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);

        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNodes.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);

        std::vector<Key>recvNodekeys;
        std::vector<Key>::iterator itKey;
        unsigned int sResult;


       /* m_uiMaxDepth++;
        treeNodesTovtk(allocatedNodes,rank,"allocatedNodes");
        treeNodesTovtk(recvNodes,rank,"recvNodes");
        treeNodesTovtk(localNodes,rank,"localNodes");
        m_uiMaxDepth--;*/

        for(unsigned int p=0;p<npes;p++)
            m_uiSendNodeCount[p]=0;

        sendNodes.clear();

        for(unsigned int p=0;p<npes;p++) {
            recvNodekeys.clear();

            for (unsigned int e = m_uiRecvNodeOffset[p]; e < (m_uiRecvNodeOffset[p] + m_uiRecvNodeCount[p]); e++) {
                itKey = recvNodekeys.emplace(recvNodekeys.end(), Key(recvNodes[e]));
                itKey->addOwner(p);
            }

            SFC::seqSearch::SFC_treeSearch(&(*(recvNodekeys.begin())), &(*(localNodes.begin())), 0, recvNodekeys.size(),0, localNodes.size(), m_uiMaxDepth + 1, m_uiMaxDepth + 1, ROOT_ROTATION);

            for (unsigned int e = 0; e < (recvNodekeys.size()); e++) {

            //NOTE: recvNodes can contain duplicates but recvNodeKeys cannot contain duplicates since we traverse by p.
            if ((recvNodekeys[e].getFlag() & OCT_FOUND)) {
                    sResult = recvNodekeys[e].getSearchResult();
                    assert(sResult >= 0 && sResult < localNodes.size());
                    m_uiScatterMapActualNodeSend.push_back(localNodes[sResult].getOwner());
                    sendNodes.push_back(localNodes[sResult]);
                    m_uiSendNodeCount[p]++;
                }
            }

        }

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);


        if(allocatedNodeSz!=(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]))
            std::cout<<"rank: "<<rank<<"[SM Error]: allocated nodes: "<<allocatedNodeSz<<" received nodes: "<<(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1])<<std::endl;


        recvNodes.clear();
        recvNodes.resize((m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]));

        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNodes.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);

        recvNodekeys.clear();
        std::vector<SearchKey> recvNodeSKeys;
        recvNodeSKeys.resize(recvNodes.size());
        for(unsigned int e=0;e<recvNodes.size();e++)
        {
            recvNodeSKeys[e]=SearchKey(recvNodes[e]);
            recvNodeSKeys[e].addOwner(e);
        }


        SFC::seqSort::SFC_treeSort(&(*(recvNodeSKeys.begin())),recvNodeSKeys.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(recvNodeSKeys));
        m_uiMaxDepth--;

        m_uiScatterMapActualNodeRecv.resize(recvNodes.size());
        unsigned int alCount=0;
        for(int e=0;e<recvNodeSKeys.size();e++)
        {
            if(allocatedNodes[alCount].getOwner()<0)
            {
                e--;
                alCount++;
                continue;
            }

            if(allocatedNodes[alCount]!=recvNodeSKeys[e]) {
               std::cout << "rank: " << rank << " allocated[" << alCount << "]: " << allocatedNodes[alCount]<< " received[" << e << "]: " << recvNodeSKeys[e] << std::endl;
               exit(0);
            }

            m_uiScatterMapActualNodeRecv[recvNodeSKeys[e].getOwner()]=allocatedNodes[alCount].getOwner();
            alCount++;

        }

        m_uiCG2DG.clear();
        m_uiDG2CG.clear();
        localNodes.clear();
        allocatedNodes.clear();


    }


    void Mesh::computeNodalScatterMap1(MPI_Comm comm)
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(npes<=1) return; // nothing to do in the sequential case. (No scatter map required.)

        unsigned int x,y,z,sz; // x y z and size of an octant.
        unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
        unsigned int nodeIndex;

        std::vector<SearchKey> allocatedNodes;
        std::vector<SearchKey> localNodes;
        std::vector<SearchKey>::iterator it;
        std::vector<SearchKey> tmpSkeys;
        std::vector<Key> tmpKeys;

        SearchKey rootSKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);
        Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);

        ot::TreeNode minMaxLocalNode[2];
        m_uiSplitterNodes=new ot::TreeNode[2*npes]; // both min and max (nodal)
        std::vector<unsigned int > minMaxIDs;
        minMaxIDs.resize(2*npes);
        unsigned int domain_max=1u<<(m_uiMaxDepth);

        for(unsigned int e=m_uiNodeLocalBegin;e<m_uiNodeLocalEnd;e++)
        {
            dg2eijk(m_uiCG2DG[e],ownerID,ii_x,jj_y,kk_z);
            x = m_uiAllElements[ownerID].getX();
            y = m_uiAllElements[ownerID].getY();
            z = m_uiAllElements[ownerID].getZ();
            sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
            assert(sz%m_uiElementOrder==0);
            it=localNodes.emplace(localNodes.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
            it->addOwner(e);
        }


        SFC::seqSort::SFC_treeSort(&(*(localNodes.begin())),localNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(localNodes));
        m_uiMaxDepth--;

        minMaxLocalNode[0]=localNodes.front();
        minMaxLocalNode[1]=localNodes.back();



        par::Mpi_Allgather(minMaxLocalNode,m_uiSplitterNodes,2,comm);


        std::vector<bool > g1Visited;
        g1Visited.resize(m_uiCG2DG.size(),false);

        for(unsigned int ele=0;ele<m_uiGhostElementRound1Index.size();ele++) {
            for (unsigned int node = 0; node < m_uiNpE; node++) {
                nodeIndex = m_uiE2NMapping_CG[m_uiGhostElementRound1Index[ele] * m_uiNpE + node];
                if ((!(nodeIndex >= m_uiNodeLocalBegin && nodeIndex < m_uiNodeLocalEnd)) && (!g1Visited[nodeIndex])) {
                    assert(nodeIndex<g1Visited.size());
                    dg2eijk(m_uiCG2DG[nodeIndex],ownerID,ii_x,jj_y,kk_z);
                    x = m_uiAllElements[ownerID].getX();
                    y = m_uiAllElements[ownerID].getY();
                    z = m_uiAllElements[ownerID].getZ();
                    sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
                    it=allocatedNodes.emplace(allocatedNodes.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                    it->addOwner(nodeIndex);
                    g1Visited[nodeIndex]= true;
                }
            }
        }

        // number of allocated nodes without splitters. (By construction this should be unique)
        unsigned int allocatedNodeSz=allocatedNodes.size();

        for(unsigned int p=0;p<2*npes;p++)
            allocatedNodes.emplace(allocatedNodes.end(),m_uiSplitterNodes[p]);

        SFC::seqSort::SFC_treeSort(&(*(allocatedNodes.begin())),allocatedNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


        tmpSkeys.clear();
        SearchKey tmpSkey;
        unsigned int skip;
        for(unsigned int e=0;e<(allocatedNodes.size());e++)
        {
            skip=1;
            tmpSkey=allocatedNodes[e];
            while(((e+skip)<allocatedNodes.size()) && (allocatedNodes[e]==allocatedNodes[e+skip]))
            {
                if(allocatedNodes[e+skip].getOwner()>=0)tmpSkey.addOwner(allocatedNodes[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(allocatedNodes,tmpSkeys);
        tmpSkeys.clear();

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(allocatedNodes));
        m_uiMaxDepth--;


        std::vector<ot::Key> splitterNode_keys;
        splitterNode_keys.resize(2*npes);
        for(unsigned int p=0;p<2*npes;p++)
            splitterNode_keys[p]=ot::Key(m_uiSplitterNodes[p]);


        m_uiMaxDepth++;
        searchKeys(splitterNode_keys,allocatedNodes);
        m_uiMaxDepth--;

        for(unsigned int p=0;p<2*npes;p++)
        {
            assert(splitterNode_keys[p].getFlag() & OCT_FOUND);
            minMaxIDs[p]=splitterNode_keys[p].getSearchResult();
            assert(minMaxIDs[p]<allocatedNodes.size());
        }


        m_uiScatterMapActualNodeSend.clear();
        m_uiScatterMapActualNodeRecv.clear();

        m_uiSendNodeCount.resize(npes);
        m_uiRecvNodeCount.resize(npes);
        m_uiSendNodeOffset.resize(npes);
        m_uiRecvNodeOffset.resize(npes);

        for(unsigned int p=0;p<npes;p++)
            m_uiSendNodeCount[p]=0;


        std::vector<ot::TreeNode> sendNodes;
        std::vector<ot::TreeNode> recvNodes;

        std::vector<ot::TreeNode> sendElements;
        std::vector<ot::TreeNode> recvElements;
        std::vector<ot::TreeNode>::iterator itTN;
        std::vector<ot::TreeNode> tmpOcts;
        ot::TreeNode tmpOct;
        ot::TreeNode rootOct(m_uiDim,m_uiMaxDepth);
        unsigned int nodeFlag;


        for(unsigned int p=0;p<npes;p++)
        {
            if(p==rank) continue;
            sendElements.clear();
            for(unsigned int e=minMaxIDs[2*p];e<(minMaxIDs[2*p+1]+1);e++)
            {
                if(allocatedNodes[e].getOwner()>=0){
                    dg2eijk(m_uiCG2DG[allocatedNodes[e].getOwner()],ownerID,ii_x,jj_y,kk_z);
                    nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);
                    itTN=sendElements.emplace(sendElements.end(),m_uiAllElements[ownerID]);
                    itTN->setFlag((itTN->getLevel())|(1u<<(nodeFlag+CHAINED_GHOST_OFFSET)));
                }
            }


            SFC::seqSort::SFC_treeSort(&(*(sendElements.begin())),sendElements.size(),tmpOcts,tmpOcts,tmpOcts,m_uiMaxDepth,m_uiMaxDepth,rootOct,ROOT_ROTATION,1,TS_SORT_ONLY);
            for(unsigned int e=0;e<(sendElements.size());e++)
            {
                itTN=sendNodes.emplace(sendNodes.end(),sendElements[e]);
                skip=1;
                while(((e+skip)<sendElements.size()) && (sendElements[e]==sendElements[e+skip]))
                {
                    itTN->setFlag((itTN->getFlag()) | (sendElements[e+skip].getFlag()));
                    skip++;
                }
                m_uiSendNodeCount[p]++;

                e+=(skip-1);

            }

        }


        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        assert(sendNodes.size()==(m_uiSendNodeOffset[npes-1]+m_uiSendNodeCount[npes-1]));
        recvElements.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);
        double t1,t2,t_stat,t_stat_g;
        DendroIntL localSz;
        DendroIntL stat_sz[3];

        localSz=sendNodes.size();

        t1=MPI_Wtime();
        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvElements.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);

        par::Mpi_Reduce(&localSz,&stat_sz[0],1,MPI_MIN,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[1],1,MPI_SUM,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[2],1,MPI_MAX,0,comm);
        stat_sz[1]=stat_sz[1]/npes;

        if(!rank) std::cout<<"a2a_ 1 time max: "<<t_stat_g<<std::endl;
        if(!rank) std::cout<<"a2a_ 1 sz  min mean max : "<<stat_sz[0]<<", "<<stat_sz[1]<<", "<<stat_sz[2]<<std::endl;

        recvNodes.clear();
        unsigned int hSz;
        unsigned int* recvNodeCount=new unsigned int [npes];
        for(unsigned int p=0;p<npes;p++)
            recvNodeCount[p]=0;


        t1=MPI_Wtime();
        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=m_uiRecvNodeOffset[p];e<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p]);e++)
            {
                nodeFlag=recvElements[e].getFlag();
                nodeFlag=nodeFlag>>(CHAINED_GHOST_OFFSET);
                sz=1u<<(m_uiMaxDepth-recvElements[e].getLevel());
                assert((sz%m_uiElementOrder)==0);
                hSz=sz/m_uiElementOrder;
                x=recvElements[e].getX();
                y=recvElements[e].getY();
                z=recvElements[e].getZ();

                if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_LEFT_UP_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+sz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+sz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }


                if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_LEFT_UP_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+sz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+sz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                // face internal.
                if(m_uiElementOrder>1)
                {
                    if(nodeFlag & (1u<<OCT_DIR_LEFT))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int j=1;j<m_uiElementOrder;j++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+j*hSz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);

                    }


                    if(nodeFlag & (1u<<OCT_DIR_RIGHT))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int j=1;j<m_uiElementOrder;j++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+j*hSz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);


                    }
                    if(nodeFlag & (1u<<OCT_DIR_DOWN))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_UP))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+sz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_BACK))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+j*hSz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_FRONT))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+j*hSz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN))
                    {

                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_UP))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+sz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_BACK))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+j*hSz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_FRONT))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+j*hSz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }


                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+sz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_BACK))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+j*hSz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_FRONT))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+j*hSz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);


                    }


                    if(nodeFlag & (1u<<OCT_DIR_DOWN_BACK))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_DOWN_FRONT))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_UP_BACK))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+sz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_UP_FRONT))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+sz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_INTERNAL))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int j=1;j<m_uiElementOrder;j++)
                                for(unsigned int i=1;i<m_uiElementOrder;i++)
                                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+j*hSz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));

                        recvNodeCount[p]+=((m_uiElementOrder-1)*(m_uiElementOrder-1)*(m_uiElementOrder-1));
                    }



                }
            }
        }

        for(unsigned int p=0;p<npes;p++)
            m_uiRecvNodeCount[p]=recvNodeCount[p];

        delete [] recvNodeCount;
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" unzip time: "<<t_stat_g<<std::endl;


        m_uiRecvNodeOffset[0]=0;
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);


        std::vector<Key>recvNodekeys;
        std::vector<Key>::iterator itKey;
        unsigned int sResult;

        /* m_uiMaxDepth++;
         treeNodesTovtk(allocatedNodes,rank,"allocatedNodes");
         treeNodesTovtk(recvNodes,rank,"recvNodes");
         treeNodesTovtk(localNodes,rank,"localNodes");
         m_uiMaxDepth--;*/

        for(unsigned int p=0;p<npes;p++)
            m_uiSendNodeCount[p]=0;

        sendNodes.clear();

        for(unsigned int p=0;p<npes;p++) {
            for (unsigned int e = m_uiRecvNodeOffset[p]; e < (m_uiRecvNodeOffset[p] + m_uiRecvNodeCount[p]); e++) {
                itKey = recvNodekeys.emplace(recvNodekeys.end(), Key(recvNodes[e]));
                itKey->addOwner(p);
            }
        }

        t1=MPI_Wtime();
        SFC::seqSearch::SFC_treeSearch(&(*(recvNodekeys.begin())), &(*(localNodes.begin())), 0, recvNodekeys.size(),0, localNodes.size(), m_uiMaxDepth + 1, m_uiMaxDepth + 1, ROOT_ROTATION);
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" search time: "<<t_stat_g<<std::endl;


        t1=MPI_Wtime();
        std::vector<unsigned int > * sendResultID=new std::vector<unsigned int >[npes];
        for (unsigned int e = 0; e < (recvNodekeys.size()); e++) {

            //NOTE: recvNodes can contain duplicates but recvNodeKeys cannot contain duplicates since we traverse by p.
            if ((recvNodekeys[e].getFlag() & OCT_FOUND)) {
                sResult = recvNodekeys[e].getSearchResult();
                assert(sResult >= 0 && sResult < localNodes.size());
                sendResultID[recvNodekeys[e].getOwnerList()->front()].push_back(sResult);
            }
        }

        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=0;e<sendResultID[p].size();e++)
            {
                m_uiScatterMapActualNodeSend.push_back(localNodes[sendResultID[p][e]].getOwner());
                sendNodes.push_back(localNodes[sendResultID[p][e]]);
            }
            m_uiSendNodeCount[p]=sendResultID[p].size();
            sendResultID[p].clear();
        }

        delete [] sendResultID;
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" send sm time: "<<t_stat_g<<std::endl;

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);


        if(allocatedNodeSz!=(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]))
            std::cout<<"rank: "<<rank<<"[SM Error]: allocated nodes: "<<allocatedNodeSz<<" received nodes: "<<(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1])<<std::endl;


        recvNodes.clear();
        recvNodes.resize((m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]));
        t1=MPI_Wtime();
        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNodes.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);
        t2=MPI_Wtime();
        t_stat=t2-t1;

        localSz=sendNodes.size();
        par::Mpi_Reduce(&localSz,&stat_sz[0],1,MPI_MIN,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[1],1,MPI_SUM,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[2],1,MPI_MAX,0,comm);
        stat_sz[1]=stat_sz[1]/npes;

        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<"a2a_ 2 time max: "<<t_stat_g<<std::endl;
        if(!rank) std::cout<<"a2a_ 2 sz  min mean max : "<<stat_sz[0]<<", "<<stat_sz[1]<<", "<<stat_sz[2]<<std::endl;

        recvNodekeys.clear();
        std::vector<SearchKey> recvNodeSKeys;
        recvNodeSKeys.resize(recvNodes.size());
        for(unsigned int e=0;e<recvNodes.size();e++)
        {
            recvNodeSKeys[e]=SearchKey(recvNodes[e]);
            recvNodeSKeys[e].addOwner(e);
        }

        t1=MPI_Wtime();
        SFC::seqSort::SFC_treeSort(&(*(recvNodeSKeys.begin())),recvNodeSKeys.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" sort time: "<<t_stat_g<<std::endl;


        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(recvNodeSKeys));
        m_uiMaxDepth--;

        t1=MPI_Wtime();
        m_uiScatterMapActualNodeRecv.resize(recvNodes.size());
        unsigned int alCount=0;
        for(int e=0;e<recvNodeSKeys.size();e++)
        {
            if(allocatedNodes[alCount].getOwner()<0)
            {
                e--;
                alCount++;
                continue;
            }

            if(allocatedNodes[alCount]!=recvNodeSKeys[e]) {
                std::cout << "rank: " << rank << " allocated[" << alCount << "]: " << allocatedNodes[alCount]<< " received[" << e << "]: " << recvNodeSKeys[e] << std::endl;
                exit(0);
            }

            m_uiScatterMapActualNodeRecv[recvNodeSKeys[e].getOwner()]=allocatedNodes[alCount].getOwner();
            alCount++;

        }

        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" sm recv update time: "<<t_stat_g<<std::endl;

        m_uiCG2DG.clear();
        m_uiDG2CG.clear();
        localNodes.clear();
        allocatedNodes.clear();

    }


    void Mesh::computeNodalScatterMap2(MPI_Comm comm)
    {
        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(npes<=1) return; // nothing to do in the sequential case. (No scatter map required.)

        unsigned int x,y,z,sz; // x y z and size of an octant.
        unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
        unsigned int nodeIndex;
        unsigned int lookUp;

        std::vector<SearchKey> localNodes;
        std::vector<SearchKey>::iterator it;
        std::vector<SearchKey> tmpSkeys;
        std::vector<Key> tmpKeys;

        SearchKey rootSKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);
        Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);
        const unsigned int domain_max=1u<<(m_uiMaxDepth);

        // 1. generate the local & sort the local nodes. (this should be unique and sorted)
        for(unsigned int e=m_uiNodeLocalBegin;e<m_uiNodeLocalEnd;e++)
        {
            dg2eijk(m_uiCG2DG[e],ownerID,ii_x,jj_y,kk_z);
            x = m_uiAllElements[ownerID].getX();
            y = m_uiAllElements[ownerID].getY();
            z = m_uiAllElements[ownerID].getZ();
            sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
            assert(sz%m_uiElementOrder==0);
            it=localNodes.emplace(localNodes.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
            it->addOwner(e);
        }


        SFC::seqSort::SFC_treeSort(&(*(localNodes.begin())),localNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(localNodes));
        m_uiMaxDepth--;

        //2. compute the local splitters. We compute 3 splitters. min max(in m_uiMaxDepth domain), & max (in m_uiMaxDepth +1)
        ot::TreeNode localNodeSpliters[3];
        assert((localNodes[0].getX()<domain_max) && (localNodes[0].getY()<domain_max) && (localNodes[0].getZ()<domain_max));
        localNodeSpliters[0]=localNodes.front();
        localNodeSpliters[2]=localNodes.back();

        for(int i=(localNodes.size()-1);i>=0;i--)
        {
            if((localNodes[i].getX()<domain_max) && (localNodes[i].getY()<domain_max) && (localNodes[i].getZ()<domain_max))
            {
                localNodeSpliters[1]=localNodes[i];
                break;
            }

        }

        //3. gather all the local splitters.
        m_uiSplitterNodes=new ot::TreeNode[3*npes];
        par::Mpi_Allgather(localNodeSpliters,m_uiSplitterNodes,3,comm);

        //4. compute the ownership (which processor it belongs to) all the ghost elements.
        std::vector<unsigned int> elementOwner;
        elementOwner.resize(m_uiAllElements.size(),rank);

        std::vector<ot::SearchKey> ghostElements;
        std::vector<ot::SearchKey>::iterator itSKey;
        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPreGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int e=m_uiElementPostGhostBegin;e<m_uiElementPostGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int p=0;p<npes;p++)
            ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiLocalSplitterElements[2*p+1]));

        SFC::seqSort::SFC_treeSort(&(*(ghostElements.begin())),ghostElements.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth,m_uiMaxDepth,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        tmpSkeys.clear();
        SearchKey tmpSkey;
        unsigned int skip;
        for(unsigned int e=0;e<(ghostElements.size());e++)
        {
            skip=1;
            tmpSkey=ghostElements[e];
            while(((e+skip)<ghostElements.size()) && (ghostElements[e]==ghostElements[e+skip]))
            {
                if(ghostElements[e+skip].getOwner()>=0)tmpSkey.addOwner(ghostElements[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(ghostElements,tmpSkeys);
        tmpSkeys.clear();


        unsigned int gCount=0;
        for(unsigned int p=0;p<npes;p++)
        {

            while(gCount<ghostElements.size() && (ghostElements[gCount]!=m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;

                gCount++;
            }

            if(gCount<ghostElements.size() && (ghostElements[gCount]==m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;
                gCount++;
            }


        }


        std::vector<SearchKey>* allocated_p=new std::vector<SearchKey>[npes];

        std::vector<unsigned int > lookUps;
        lookUps.resize(NUM_CHILDREN);
        std::set<unsigned int> ownerPoc;

        unsigned int nodeFlag;
        std::vector<unsigned int > internalIndex;
        std::vector<unsigned int > faceIndex;   // face internal
        std::vector<unsigned int > edgeIndex;   // edge internal
        std::vector<unsigned int > vertexIndex; // vertex internal

        unsigned int edge_dir[]={OCT_DIR_LEFT_BACK,OCT_DIR_LEFT_FRONT,OCT_DIR_LEFT_DOWN,OCT_DIR_LEFT_FRONT,OCT_DIR_RIGHT_BACK,OCT_DIR_RIGHT_FRONT,OCT_DIR_RIGHT_DOWN,OCT_DIR_RIGHT_FRONT,OCT_DIR_UP_BACK,OCT_DIR_UP_FRONT,OCT_DIR_DOWN_BACK,OCT_DIR_DOWN_FRONT};
        unsigned int face_dir[]={OCT_DIR_LEFT,OCT_DIR_RIGHT,OCT_DIR_DOWN,OCT_DIR_UP,OCT_DIR_BACK,OCT_DIR_FRONT};
        unsigned int vertex_dir[]={OCT_DIR_LEFT_DOWN_BACK,OCT_DIR_RIGHT_DOWN_BACK,OCT_DIR_LEFT_UP_BACK,OCT_DIR_RIGHT_UP_BACK, OCT_DIR_LEFT_DOWN_FRONT,OCT_DIR_RIGHT_DOWN_FRONT,OCT_DIR_LEFT_UP_FRONT,OCT_DIR_RIGHT_UP_FRONT};

        std::vector<bool > g1Visited;
        g1Visited.resize(m_uiCG2DG.size(),false);

        std::vector<SearchKey> allocated1; // allocated nodes where the owner of the nodes are undecided.
        std::vector<SearchKey> allocatedNodes; // actual allocated nodes.

        for(unsigned int ele=0;ele<m_uiGhostElementRound1Index.size();ele++) {
              for (unsigned int node = 0; node < m_uiNpE; node++) {
                nodeIndex = m_uiE2NMapping_CG[m_uiGhostElementRound1Index[ele] * m_uiNpE + node];
                if ((!(nodeIndex >= m_uiNodeLocalBegin && nodeIndex < m_uiNodeLocalEnd)) && (!g1Visited[nodeIndex])) {
                    assert(nodeIndex<g1Visited.size());
                    dg2eijk(m_uiCG2DG[nodeIndex],ownerID,ii_x,jj_y,kk_z);
                    nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);
                    x = m_uiAllElements[ownerID].getX();
                    y = m_uiAllElements[ownerID].getY();
                    z = m_uiAllElements[ownerID].getZ();
                    sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());

                    if(nodeFlag == OCT_DIR_INTERNAL)
                    { // for internal nodes we can directly determine the ownership.
                        it=allocated_p[elementOwner[ownerID]].emplace(allocated_p[elementOwner[ownerID]].end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                        it->addOwner(nodeIndex);
                    }else
                    { // for other nodes we use the modified splitter approach.
                        it=allocated1.emplace(allocated1.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                        it->addOwner(nodeIndex);
                    }
                    g1Visited[nodeIndex]= true;
                }
            }
        }

        allocatedNodes.clear();
        allocatedNodes.resize(allocated1.size());
        allocatedNodes.assign(allocated1.begin(),allocated1.end());

        for(unsigned int p=0;p<npes;p++)
            allocatedNodes.insert(allocatedNodes.end(),allocated_p[p].begin(),allocated_p[p].end());

        unsigned  int allocatedNodeSz=allocatedNodes.size();


        for(unsigned int p=0;p<npes;p++)
        {
            allocated1.emplace(allocated1.end(),SearchKey(m_uiSplitterNodes[3*p]));
            allocated1.emplace(allocated1.end(),SearchKey(m_uiSplitterNodes[3*p+1]));
            allocated1.emplace(allocated1.end(),SearchKey(m_uiSplitterNodes[3*p+2]));
        }


        SFC::seqSort::SFC_treeSort(&(*(allocated1.begin())),allocated1.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


        tmpSkeys.clear();
        for(unsigned int e=0;e<(allocated1.size());e++)
        {
            skip=1;
            tmpSkey=allocated1[e];
            while(((e+skip)<allocated1.size()) && (allocated1[e]==allocated1[e+skip]))
            {
                if(allocated1[e+skip].getOwner()>=0)tmpSkey.addOwner(allocated1[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(allocated1,tmpSkeys);
        tmpSkeys.clear();

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(allocated1));
        m_uiMaxDepth--;

        std::vector<ot::Key> splitterNode_keys;
        splitterNode_keys.resize(3*npes);
        std::vector<unsigned int > nodeSplitterID;
        nodeSplitterID.resize(3*npes);

        for(unsigned int p=0;p<3*npes;p++)
            splitterNode_keys[p]=ot::Key(m_uiSplitterNodes[p]);


        m_uiMaxDepth++;
        searchKeys(splitterNode_keys,allocated1);
        m_uiMaxDepth--;

        for(unsigned int p=0;p<3*npes;p++)
        {
            assert(splitterNode_keys[p].getFlag() & OCT_FOUND);
            nodeSplitterID[p]=splitterNode_keys[p].getSearchResult();
            assert(nodeSplitterID[p]<allocated1.size());
        }


        for(unsigned int p=0;p<npes;p++)
        {
            if(p==rank) continue;
            for(unsigned int e=nodeSplitterID[3*p];e<(nodeSplitterID[3*p+1]+1);e++)
            {
                if(allocated1[e].getOwner()>=0 ){
                    assert((allocated1[e].getX()<domain_max) && (allocated1[e].getY()<domain_max) && (allocated1[e].getZ()<domain_max));
                    allocated_p[p].push_back(allocated1[e]);
                }
            }


            for(unsigned int e=nodeSplitterID[3*p+1]+1;e<(nodeSplitterID[3*p+2]+1);e++)
            {
                if(allocated1[e].getOwner()>=0 && (!((allocated1[e].getX()<domain_max) && (allocated1[e].getY()<domain_max) && (allocated1[e].getZ()<domain_max))) ){
                    allocated_p[p].push_back(allocated1[e]);
                }
            }

        }

        allocated1.clear();


        m_uiScatterMapActualNodeSend.clear();
        m_uiScatterMapActualNodeRecv.clear();

        m_uiSendNodeCount.resize(npes);
        m_uiRecvNodeCount.resize(npes);
        m_uiSendNodeOffset.resize(npes);
        m_uiRecvNodeOffset.resize(npes);

        std::vector<ot::TreeNode> sendNodes;
        std::vector<ot::TreeNode> recvNodes;

        double t1,t2,t_stat,t_stat_g;
        DendroIntL localSz;
        DendroIntL stat_sz[3];

        for(unsigned int p=0;p<npes;p++)
        {
            m_uiSendNodeCount[p]=allocated_p[p].size();
            sendNodes.insert(sendNodes.end(),allocated_p[p].begin(),allocated_p[p].end());
            allocated_p[p].clear();
        }

        delete [] allocated_p;


        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        assert(sendNodes.size()==(m_uiSendNodeOffset[npes-1]+m_uiSendNodeCount[npes-1]));
        recvNodes.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);

        localSz=sendNodes.size();
        t1=MPI_Wtime();
        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNodes.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);

        par::Mpi_Reduce(&localSz,&stat_sz[0],1,MPI_MIN,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[1],1,MPI_SUM,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[2],1,MPI_MAX,0,comm);
        stat_sz[1]=stat_sz[1]/npes;

        if(!rank) std::cout<<"a2a_ 1 time max: "<<t_stat_g<<std::endl;
        if(!rank) std::cout<<"a2a_ 1 sz  min mean max : "<<stat_sz[0]<<", "<<stat_sz[1]<<", "<<stat_sz[2]<<std::endl;


        std::vector<Key>recvNodekeys;
        std::vector<Key>::iterator itKey;
        unsigned int sResult;


        /* m_uiMaxDepth++;
         treeNodesTovtk(allocatedNodes,rank,"allocatedNodes");
         treeNodesTovtk(recvNodes,rank,"recvNodes");
         treeNodesTovtk(localNodes,rank,"localNodes");
         m_uiMaxDepth--;*/

        for(unsigned int p=0;p<npes;p++)
            m_uiSendNodeCount[p]=0;

        sendNodes.clear();

        for(unsigned int p=0;p<npes;p++) {
            recvNodekeys.clear();

            for (unsigned int e = m_uiRecvNodeOffset[p]; e < (m_uiRecvNodeOffset[p] + m_uiRecvNodeCount[p]); e++) {
                itKey = recvNodekeys.emplace(recvNodekeys.end(), Key(recvNodes[e]));
                itKey->addOwner(p);
            }

            SFC::seqSearch::SFC_treeSearch(&(*(recvNodekeys.begin())), &(*(localNodes.begin())), 0, recvNodekeys.size(),0, localNodes.size(), m_uiMaxDepth + 1, m_uiMaxDepth + 1, ROOT_ROTATION);

            for (unsigned int e = 0; e < (recvNodekeys.size()); e++) {

                //NOTE: recvNodes can contain duplicates but recvNodeKeys cannot contain duplicates since we traverse by p.
                if ((recvNodekeys[e].getFlag() & OCT_FOUND)) {
                    sResult = recvNodekeys[e].getSearchResult();
                    assert(sResult >= 0 && sResult < localNodes.size());
                    m_uiScatterMapActualNodeSend.push_back(localNodes[sResult].getOwner());
                    sendNodes.push_back(localNodes[sResult]);
                    m_uiSendNodeCount[p]++;
                }/*else
                {
                    std::cout<<" key: recv form "<<p<<" to rank: "<<rank<<" key: "<<recvNodekeys[e]<<" not found "<<std::endl;
                }*/
            }

        }

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);


        if(allocatedNodeSz!=(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]))
            std::cout<<"rank: "<<rank<<"[SM Error]: allocated nodes: "<<allocatedNodeSz<<" received nodes: "<<(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1])<<std::endl;


        recvNodes.clear();
        recvNodes.resize((m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]));

        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNodes.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);

        recvNodekeys.clear();
        std::vector<SearchKey> recvNodeSKeys;
        recvNodeSKeys.resize(recvNodes.size());
        for(unsigned int e=0;e<recvNodes.size();e++)
        {
            recvNodeSKeys[e]=SearchKey(recvNodes[e]);
            recvNodeSKeys[e].addOwner(e);
        }


        SFC::seqSort::SFC_treeSort(&(*(recvNodeSKeys.begin())),recvNodeSKeys.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);
        SFC::seqSort::SFC_treeSort(&(*(allocatedNodes.begin())),allocatedNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(recvNodeSKeys));
        m_uiMaxDepth--;

        m_uiScatterMapActualNodeRecv.resize(recvNodes.size());
        unsigned int alCount=0;
        for(int e=0;e<recvNodeSKeys.size();e++)
        {
            if(allocatedNodes[alCount].getOwner()<0)
            {
                e--;
                alCount++;
                continue;
            }

            if(allocatedNodes[alCount]!=recvNodeSKeys[e]) {
                std::cout << "rank: " << rank << " allocated[" << alCount << "]: " << allocatedNodes[alCount]<< " received[" << e << "]: " << recvNodeSKeys[e] << std::endl;
                exit(0);
            }

            m_uiScatterMapActualNodeRecv[recvNodeSKeys[e].getOwner()]=allocatedNodes[alCount].getOwner();
            alCount++;

        }

        m_uiCG2DG.clear();
        m_uiDG2CG.clear();
        localNodes.clear();
        allocatedNodes.clear();


    }

    void Mesh::computeNodalScatterMap3(MPI_Comm comm)
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(npes<=1) return; // nothing to do in the sequential case. (No scatter map required.)

        unsigned int x,y,z,sz; // x y z and size of an octant.
        unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
        unsigned int ownerID1;
        unsigned int nodeIndex;
        unsigned int lookUp;


        double t1,t2,t_stat,t_stat_g;
        DendroIntL localSz;
        DendroIntL stat_sz[3];


        std::vector<SearchKey> localNodes;
        std::vector<SearchKey>::iterator it;
        std::vector<SearchKey> tmpSkeys;
        std::vector<Key> tmpKeys;

        SearchKey rootSKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);
        Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth+1);
        const unsigned int domain_max=1u<<(m_uiMaxDepth);

        t1=MPI_Wtime();
        // 1. generate the local & sort the local nodes. (this should be unique and sorted)
        for(unsigned int e=m_uiNodeLocalBegin;e<m_uiNodeLocalEnd;e++)
        {
            dg2eijk(m_uiCG2DG[e],ownerID,ii_x,jj_y,kk_z);
            x = m_uiAllElements[ownerID].getX();
            y = m_uiAllElements[ownerID].getY();
            z = m_uiAllElements[ownerID].getZ();
            sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
            assert(sz%m_uiElementOrder==0);
            it=localNodes.emplace(localNodes.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
            it->addOwner(e);
        }


        SFC::seqSort::SFC_treeSort(&(*(localNodes.begin())),localNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" local node generation + sort time (max) (s): "<<t_stat_g<<std::endl;


        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(localNodes));
        m_uiMaxDepth--;

        //2. compute the local splitters. We compute 3 splitters. min max(in m_uiMaxDepth domain), & max (in m_uiMaxDepth +1)
        ot::TreeNode localNodeSpliters[3];
        assert((localNodes[0].getX()<domain_max) && (localNodes[0].getY()<domain_max) && (localNodes[0].getZ()<domain_max));
        localNodeSpliters[0]=localNodes.front();
        localNodeSpliters[2]=localNodes.back();

        for(int i=(localNodes.size()-1);i>=0;i--)
        {
            if((localNodes[i].getX()<domain_max) && (localNodes[i].getY()<domain_max) && (localNodes[i].getZ()<domain_max))
            {
                localNodeSpliters[1]=localNodes[i];
                break;
            }

        }

        //3. gather all the local splitters.
        m_uiSplitterNodes=new ot::TreeNode[3*npes];
        par::Mpi_Allgather(localNodeSpliters,m_uiSplitterNodes,3,comm);


        t1=MPI_Wtime();
        //4. compute the ownership (which processor it belongs to) all the ghost elements.
        std::vector<unsigned int> elementOwner;
        elementOwner.resize(m_uiAllElements.size(),rank);

        std::vector<ot::SearchKey> ghostElements;
        std::vector<ot::SearchKey>::iterator itSKey;
        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPreGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int e=m_uiElementPostGhostBegin;e<m_uiElementPostGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int p=0;p<npes;p++)
            ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiLocalSplitterElements[2*p+1]));

        SFC::seqSort::SFC_treeSort(&(*(ghostElements.begin())),ghostElements.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth,m_uiMaxDepth,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        tmpSkeys.clear();
        SearchKey tmpSkey;
        unsigned int skip;
        for(unsigned int e=0;e<(ghostElements.size());e++)
        {
            skip=1;
            tmpSkey=ghostElements[e];
            while(((e+skip)<ghostElements.size()) && (ghostElements[e]==ghostElements[e+skip]))
            {
                if(ghostElements[e+skip].getOwner()>=0)tmpSkey.addOwner(ghostElements[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(ghostElements,tmpSkeys);
        tmpSkeys.clear();


        unsigned int gCount=0;
        for(unsigned int p=0;p<npes;p++)
        {

            while(gCount<ghostElements.size() && (ghostElements[gCount]!=m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;

                gCount++;
            }

            if(gCount<ghostElements.size() && (ghostElements[gCount]==m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;
                gCount++;
            }


        }
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" ghost element ownership build time (max) (s): "<<t_stat_g<<std::endl;


        std::vector<SearchKey>* allocated_p=new std::vector<SearchKey>[npes];

        std::vector<unsigned int > lookUps;
        lookUps.resize(NUM_CHILDREN);
        std::set<unsigned int> ownerPoc;

        unsigned int nodeFlag;
        std::vector<unsigned int > internalIndex;
        std::vector<unsigned int > faceIndex;   // face internal
        std::vector<unsigned int > edgeIndex;   // edge internal
        std::vector<unsigned int > vertexIndex; // vertex internal

        /*unsigned int edge_dir[]={OCT_DIR_LEFT_BACK,OCT_DIR_LEFT_FRONT,OCT_DIR_LEFT_DOWN,OCT_DIR_LEFT_FRONT,OCT_DIR_RIGHT_BACK,OCT_DIR_RIGHT_FRONT,OCT_DIR_RIGHT_DOWN,OCT_DIR_RIGHT_FRONT,OCT_DIR_UP_BACK,OCT_DIR_UP_FRONT,OCT_DIR_DOWN_BACK,OCT_DIR_DOWN_FRONT};
        unsigned int face_dir[]={OCT_DIR_LEFT,OCT_DIR_RIGHT,OCT_DIR_DOWN,OCT_DIR_UP,OCT_DIR_BACK,OCT_DIR_FRONT};
        unsigned int vertex_dir[]={OCT_DIR_LEFT_DOWN_BACK,OCT_DIR_RIGHT_DOWN_BACK,OCT_DIR_LEFT_UP_BACK,OCT_DIR_RIGHT_UP_BACK, OCT_DIR_LEFT_DOWN_FRONT,OCT_DIR_RIGHT_DOWN_FRONT,OCT_DIR_LEFT_UP_FRONT,OCT_DIR_RIGHT_UP_FRONT};*/

        std::vector<bool > g1Visited;
        g1Visited.resize(m_uiCG2DG.size(),false);

        std::vector<SearchKey> allocated1; // allocated nodes where the owner of the nodes are undecided.
        std::vector<SearchKey> allocatedNodes; // actual allocated nodes.
        t1=MPI_Wtime();
        for(unsigned int ele=0;ele<m_uiGhostElementRound1Index.size();ele++) {
            for (unsigned int node = 0; node < m_uiNpE; node++) {
                nodeIndex = m_uiE2NMapping_CG[m_uiGhostElementRound1Index[ele] * m_uiNpE + node];
                if ((!(nodeIndex >= m_uiNodeLocalBegin && nodeIndex < m_uiNodeLocalEnd)) && (!g1Visited[nodeIndex])) {
                    assert(nodeIndex<g1Visited.size());
                    dg2eijk(m_uiCG2DG[nodeIndex],ownerID,ii_x,jj_y,kk_z);
                    nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);
                    x = m_uiAllElements[ownerID].getX();
                    y = m_uiAllElements[ownerID].getY();
                    z = m_uiAllElements[ownerID].getZ();
                    sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());

                    if(nodeFlag == OCT_DIR_INTERNAL)
                    { // for internal nodes we can directly determine the ownership.
                        it=allocated_p[elementOwner[ownerID]].emplace(allocated_p[elementOwner[ownerID]].end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                        it->addOwner(nodeIndex);
                    }else
                    { // for other nodes we use the modified splitter approach.
                        it=allocated1.emplace(allocated1.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                        it->addOwner(nodeIndex);
                    }
                    g1Visited[nodeIndex]= true;
                }
            }
        }

        allocatedNodes.clear();
        allocatedNodes.resize(allocated1.size());
        allocatedNodes.assign(allocated1.begin(),allocated1.end());

        for(unsigned int p=0;p<npes;p++)
            allocatedNodes.insert(allocatedNodes.end(),allocated_p[p].begin(),allocated_p[p].end());


        SFC::seqSort::SFC_treeSort(&(*(allocatedNodes.begin())),allocatedNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


        unsigned  int allocatedNodeSz=allocatedNodes.size();
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" allocated node generation + sort time (max) (s): "<<t_stat_g<<std::endl;

        for(unsigned int p=0;p<npes;p++)
        {
            allocated1.emplace(allocated1.end(),SearchKey(m_uiSplitterNodes[3*p]));
            allocated1.emplace(allocated1.end(),SearchKey(m_uiSplitterNodes[3*p+1]));
            allocated1.emplace(allocated1.end(),SearchKey(m_uiSplitterNodes[3*p+2]));
        }

        t1=MPI_Wtime();

        SFC::seqSort::SFC_treeSort(&(*(allocated1.begin())),allocated1.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


        tmpSkeys.clear();
        for(unsigned int e=0;e<(allocated1.size());e++)
        {
            skip=1;
            tmpSkey=allocated1[e];
            while(((e+skip)<allocated1.size()) && (allocated1[e]==allocated1[e+skip]))
            {
                if(allocated1[e+skip].getOwner()>=0)tmpSkey.addOwner(allocated1[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(allocated1,tmpSkeys);
        tmpSkeys.clear();

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(allocated1));
        m_uiMaxDepth--;

        std::vector<ot::Key> splitterNode_keys;
        splitterNode_keys.resize(3*npes);
        std::vector<unsigned int > nodeSplitterID;
        nodeSplitterID.resize(3*npes);

        for(unsigned int p=0;p<3*npes;p++)
            splitterNode_keys[p]=ot::Key(m_uiSplitterNodes[p]);


        m_uiMaxDepth++;
        searchKeys(splitterNode_keys,allocated1);
        m_uiMaxDepth--;

        for(unsigned int p=0;p<3*npes;p++)
        {
            assert(splitterNode_keys[p].getFlag() & OCT_FOUND);
            nodeSplitterID[p]=splitterNode_keys[p].getSearchResult();
            assert(nodeSplitterID[p]<allocated1.size());
        }


        for(unsigned int p=0;p<npes;p++)
        {
            if(p==rank) continue;
            for(unsigned int e=nodeSplitterID[3*p];e<(nodeSplitterID[3*p+1]+1);e++)
            {
                if(allocated1[e].getOwner()>=0 ){
                    assert((allocated1[e].getX()<domain_max) && (allocated1[e].getY()<domain_max) && (allocated1[e].getZ()<domain_max));
                    allocated_p[p].push_back(allocated1[e]);
                }
            }


            for(unsigned int e=nodeSplitterID[3*p+1]+1;e<(nodeSplitterID[3*p+2]+1);e++)
            {
                if(allocated1[e].getOwner()>=0 && (!((allocated1[e].getX()<domain_max) && (allocated1[e].getY()<domain_max) && (allocated1[e].getZ()<domain_max))) ){
                    allocated_p[p].push_back(allocated1[e]);
                }
            }

        }

        allocated1.clear();
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" allocated1 ownership computation (max) (s): "<<t_stat_g<<std::endl;


        m_uiScatterMapActualNodeSend.clear();
        m_uiScatterMapActualNodeRecv.clear();

        m_uiSendNodeCount.resize(npes);
        m_uiRecvNodeCount.resize(npes);
        m_uiSendNodeOffset.resize(npes);
        m_uiRecvNodeOffset.resize(npes);

        std::vector<ot::TreeNode> sendNodes;
        std::vector<ot::TreeNode> recvNodes;

        std::vector<ot::TreeNode> sendElements;
        std::vector<ot::TreeNode> recvElements;
        std::vector<ot::TreeNode> tmpOcts;
        ot::TreeNode rootOct(m_uiDim,m_uiMaxDepth);

        std::vector<ot::TreeNode>::iterator itTN;

        assert(allocated_p[rank].size()==0);

        t1=MPI_Wtime();
        for(unsigned int p=0;p<npes;p++)
        {
            m_uiSendNodeCount[p]=0;
            if(p==rank) continue;

            sendElements.clear();
            for(unsigned int e=0;e<allocated_p[p].size();e++)
            {
                assert(allocated_p[p][e].getOwner()>=0);
                dg2eijk(m_uiCG2DG[allocated_p[p][e].getOwner()],ownerID,ii_x,jj_y,kk_z);
                nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);
                itTN=sendElements.emplace(sendElements.end(),m_uiAllElements[ownerID]);
                itTN->setFlag((itTN->getLevel())|(1u<<(nodeFlag+CHAINED_GHOST_OFFSET)));

            }

            allocated_p[p].clear();

            SFC::seqSort::SFC_treeSort(&(*(sendElements.begin())),sendElements.size(),tmpOcts,tmpOcts,tmpOcts,m_uiMaxDepth,m_uiMaxDepth,rootOct,ROOT_ROTATION,1,TS_SORT_ONLY);
            for(unsigned int e=0;e<(sendElements.size());e++)
            {
                itTN=sendNodes.emplace(sendNodes.end(),sendElements[e]);
                skip=1;
                while(((e+skip)<sendElements.size()) && (sendElements[e]==sendElements[e+skip]))
                {
                    itTN->setFlag((itTN->getFlag()) | (sendElements[e+skip].getFlag()));
                    skip++;
                }
                m_uiSendNodeCount[p]++;
                e+=(skip-1);
            }

        }

        delete [] allocated_p;

        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" allocated node zip time (max) (s): "<<t_stat_g<<std::endl;



        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        assert(sendNodes.size()==(m_uiSendNodeOffset[npes-1]+m_uiSendNodeCount[npes-1]));
        recvElements.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);

        localSz=sendNodes.size();
        t1=MPI_Wtime();
        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvElements.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);

        par::Mpi_Reduce(&localSz,&stat_sz[0],1,MPI_MIN,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[1],1,MPI_SUM,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[2],1,MPI_MAX,0,comm);
        stat_sz[1]=stat_sz[1]/npes;

        if(!rank) std::cout<<"a2a_ 1 time max: "<<t_stat_g<<std::endl;
        if(!rank) std::cout<<"a2a_ 1 sz  min mean max : "<<stat_sz[0]<<", "<<stat_sz[1]<<", "<<stat_sz[2]<<std::endl;


        recvNodes.clear();
        unsigned int hSz;
        unsigned int* recvNodeCount=new unsigned int [npes];
        for(unsigned int p=0;p<npes;p++)
            recvNodeCount[p]=0;

       t1=MPI_Wtime();
        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=m_uiRecvNodeOffset[p];e<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p]);e++)
            {
                nodeFlag=recvElements[e].getFlag();
                nodeFlag=nodeFlag>>(CHAINED_GHOST_OFFSET);
                sz=1u<<(m_uiMaxDepth-recvElements[e].getLevel());
                assert((sz%m_uiElementOrder)==0);
                hSz=sz/m_uiElementOrder;
                x=recvElements[e].getX();
                y=recvElements[e].getY();
                z=recvElements[e].getZ();

                if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_LEFT_UP_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+sz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP_BACK))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+sz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }


                if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_LEFT_UP_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+sz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP_FRONT))
                {
                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+sz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                    recvNodeCount[p]++;
                }

                // face internal.
                if(m_uiElementOrder>1)
                {
                    if(nodeFlag & (1u<<OCT_DIR_LEFT))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int j=1;j<m_uiElementOrder;j++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+j*hSz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);

                    }


                    if(nodeFlag & (1u<<OCT_DIR_RIGHT))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int j=1;j<m_uiElementOrder;j++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+j*hSz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);


                    }
                    if(nodeFlag & (1u<<OCT_DIR_DOWN))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_UP))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+sz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_BACK))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+j*hSz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_FRONT))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            for(unsigned int i=1;i<m_uiElementOrder;i++)
                                recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+j*hSz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1)*(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_LEFT_DOWN))
                    {

                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_UP))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+sz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_BACK))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+j*hSz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_LEFT_FRONT))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x,y+j*hSz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }


                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_DOWN))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_UP))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+sz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_BACK))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+j*hSz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);

                    }

                    if(nodeFlag & (1u<<OCT_DIR_RIGHT_FRONT))
                    {
                        for(unsigned int j=1;j<m_uiElementOrder;j++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+sz,y+j*hSz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);


                    }


                    if(nodeFlag & (1u<<OCT_DIR_DOWN_BACK))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_DOWN_FRONT))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_UP_BACK))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+sz,z,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_UP_FRONT))
                    {
                        for(unsigned int i=1;i<m_uiElementOrder;i++)
                            recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+sz,z+sz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));
                        recvNodeCount[p]+=(m_uiElementOrder-1);
                    }


                    if(nodeFlag & (1u<<OCT_DIR_INTERNAL))
                    {
                        for(unsigned int k=1;k<m_uiElementOrder;k++)
                            for(unsigned int j=1;j<m_uiElementOrder;j++)
                                for(unsigned int i=1;i<m_uiElementOrder;i++)
                                    recvNodes.emplace(recvNodes.end(),ot::TreeNode(x+i*hSz,y+j*hSz,z+k*hSz,m_uiMaxDepth+1,m_uiDim,m_uiMaxDepth+1));

                        recvNodeCount[p]+=((m_uiElementOrder-1)*(m_uiElementOrder-1)*(m_uiElementOrder-1));
                    }



                }
            }
        }

        for(unsigned int p=0;p<npes;p++)
            m_uiRecvNodeCount[p]=recvNodeCount[p];

        delete [] recvNodeCount;

        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" allocated node unzip time (max) (s): "<<t_stat_g<<std::endl;



        m_uiRecvNodeOffset[0]=0;
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        std::vector<Key>recvNodekeys;
        std::vector<Key>::iterator itKey;
        unsigned int sResult;

        for(unsigned int p=0;p<npes;p++)
            m_uiSendNodeCount[p]=0;

        sendNodes.clear();


        for(unsigned int p=0;p<npes;p++) {
            for (unsigned int e = m_uiRecvNodeOffset[p]; e < (m_uiRecvNodeOffset[p] + m_uiRecvNodeCount[p]); e++) {
                itKey = recvNodekeys.emplace(recvNodekeys.end(), Key(recvNodes[e]));
                itKey->addOwner(p);
            }
        }


        t1=MPI_Wtime();
        SFC::seqSearch::SFC_treeSearch(&(*(recvNodekeys.begin())), &(*(localNodes.begin())), 0, recvNodekeys.size(),0, localNodes.size(), m_uiMaxDepth + 1, m_uiMaxDepth + 1, ROOT_ROTATION);
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" search time: "<<t_stat_g<<std::endl;

        t1=MPI_Wtime();
        std::vector<unsigned int > * sendResultID=new std::vector<unsigned int >[npes];

        for (unsigned int e = 0; e < (recvNodekeys.size()); e++) {

            //NOTE: recvNodes can contain duplicates but recvNodeKeys cannot contain duplicates since we traverse by p.
            if ((recvNodekeys[e].getFlag() & OCT_FOUND)) {
                sResult = recvNodekeys[e].getSearchResult();
                assert(sResult >= 0 && sResult < localNodes.size());
                sendResultID[recvNodekeys[e].getOwnerList()->front()].push_back(sResult);
            }
        }

        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=0;e<sendResultID[p].size();e++)
            {
                m_uiScatterMapActualNodeSend.push_back(localNodes[sendResultID[p][e]].getOwner());
                sendNodes.push_back(localNodes[sendResultID[p][e]]);
            }
            m_uiSendNodeCount[p]=sendResultID[p].size();
            sendResultID[p].clear();
        }

        delete [] sendResultID;
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" send sm time: "<<t_stat_g<<std::endl;

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);


        if(allocatedNodeSz!=(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]))
            std::cout<<"rank: "<<rank<<"[SM Error]: allocated nodes: "<<allocatedNodeSz<<" received nodes: "<<(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1])<<std::endl;


        recvNodes.clear();
        recvNodes.resize((m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]));

        localSz=sendNodes.size();

        t1=MPI_Wtime();
        par::Mpi_Alltoallv(&(*(sendNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNodes.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);

        par::Mpi_Reduce(&localSz,&stat_sz[0],1,MPI_MIN,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[1],1,MPI_SUM,0,comm);
        par::Mpi_Reduce(&localSz,&stat_sz[2],1,MPI_MAX,0,comm);
        stat_sz[1]=stat_sz[1]/npes;

        if(!rank) std::cout<<"a2a_ 2 (nodal) time max: "<<t_stat_g<<std::endl;
        if(!rank) std::cout<<"a2a_ 2 (nodal) sz  min mean max : "<<stat_sz[0]<<", "<<stat_sz[1]<<", "<<stat_sz[2]<<std::endl;

        recvNodekeys.clear();
        std::vector<SearchKey> recvNodeSKeys;
        recvNodeSKeys.resize(recvNodes.size());
        for(unsigned int e=0;e<recvNodes.size();e++)
        {
            recvNodeSKeys[e]=SearchKey(recvNodes[e]);
            recvNodeSKeys[e].addOwner(e);
        }

        t1=MPI_Wtime();
        SFC::seqSort::SFC_treeSort(&(*(recvNodeSKeys.begin())),recvNodeSKeys.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);
        t2=MPI_Wtime();

        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<"recvnodes sort time  max (s): "<<t_stat_g<<std::endl;

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(recvNodeSKeys));
        m_uiMaxDepth--;

        m_uiScatterMapActualNodeRecv.resize(recvNodes.size());
        t1=MPI_Wtime();
        unsigned int alCount=0;
        for(int e=0;e<recvNodeSKeys.size();e++)
        {
            if(allocatedNodes[alCount].getOwner()<0)
            {
                e--;
                alCount++;
                continue;
            }

            if(allocatedNodes[alCount]!=recvNodeSKeys[e]) {
                std::cout << "rank: " << rank << " allocated[" << alCount << "]: " << allocatedNodes[alCount]<< " received[" << e << "]: " << recvNodeSKeys[e] << std::endl;
                exit(0);
            }

            m_uiScatterMapActualNodeRecv[recvNodeSKeys[e].getOwner()]=allocatedNodes[alCount].getOwner();
            alCount++;

        }

        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<"recv sm time  max (s): "<<t_stat_g<<std::endl;


        m_uiCG2DG.clear();
        m_uiDG2CG.clear();
        localNodes.clear();
        allocatedNodes.clear();


    }


    void Mesh::computeNodalScatterMap4(MPI_Comm comm)
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        if(npes<=1) return; // nothing to do in the sequential case. (No scatter map required.)

#ifdef PROFILE_SM
        double t1,t2,t_stat,t_stat_g;
        DendroIntL localSz;
        DendroIntL stat_sz[3];
#endif

        unsigned int x,y,z,sz; // x y z and size of an octant.
        unsigned int ownerID,ii_x,jj_y,kk_z; // DG index to ownerID and ijk decomposition variable.
        unsigned int lookUp;
        unsigned int nodeIndex;
        unsigned int nodeFlag;

        std::vector<ot::TreeNode> sendTNElements;
        std::vector<ot::TreeNode> recvTNElements;
        std::vector<ot::Node> sendNNodal;
        std::vector<ot::Node> recvNNodal;

        ot::TreeNode rootTN(m_uiDim,m_uiMaxDepth);
        ot::Node rootNN(m_uiDim,m_uiMaxDepth);
        ot::SearchKey rootSKey(m_uiDim,m_uiMaxDepth);

        std::vector<ot::TreeNode>::iterator itTN;
        std::vector<ot::Node>::iterator itNN;
        std::vector<ot::SearchKey>::iterator itSK;
        std::vector<ot::Key>::iterator itKK;

        std::vector<SearchKey> tmpSkeys;
        std::vector<ot::Node> tmpNN;

        ot::Node tmpNode;

        m_uiSendNodeCount.resize(npes);
        m_uiRecvNodeCount.resize(npes);
        m_uiSendNodeOffset.resize(npes);
        m_uiRecvNodeOffset.resize(npes);



        //1. compute the ownership (which processor it belongs to) all the ghost elements.
#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        std::vector<unsigned int> elementOwner;
        elementOwner.resize(m_uiAllElements.size(),rank);

        std::vector<ot::SearchKey> ghostElements;
        std::vector<ot::SearchKey>::iterator itSKey;
        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPreGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int e=m_uiElementPostGhostBegin;e<m_uiElementPostGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int p=0;p<npes;p++)
            ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiLocalSplitterElements[2*p+1]));

        SFC::seqSort::SFC_treeSort(&(*(ghostElements.begin())),ghostElements.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth,m_uiMaxDepth,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        tmpSkeys.clear();
        SearchKey tmpSkey;
        unsigned int skip;
        for(unsigned int e=0;e<(ghostElements.size());e++)
        {
            skip=1;
            tmpSkey=ghostElements[e];
            while(((e+skip)<ghostElements.size()) && (ghostElements[e]==ghostElements[e+skip]))
            {
                if(ghostElements[e+skip].getOwner()>=0)tmpSkey.addOwner(ghostElements[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(ghostElements,tmpSkeys);
        tmpSkeys.clear();


        unsigned int gCount=0;
        for(unsigned int p=0;p<npes;p++)
        {

            while(gCount<ghostElements.size() && (ghostElements[gCount]!=m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;

                gCount++;
            }

            if(gCount<ghostElements.size() && (ghostElements[gCount]==m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;
                gCount++;
            }


        }

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" ghost element ownership build time (max) (s): "<<t_stat_g<<std::endl;
#endif

        std::vector<unsigned int>* sendIDs=new std::vector<unsigned int >[npes];


        for(unsigned int e=0;e<m_uiGhostElementRound1Index.size();e++)
            sendIDs[elementOwner[m_uiGhostElementRound1Index[e]]].push_back(m_uiGhostElementRound1Index[e]);

        sendTNElements.clear();
        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=0;e<sendIDs[p].size();e++)
                sendTNElements.push_back(m_uiAllElements[sendIDs[p][e]]);

            m_uiSendNodeCount[p]=sendIDs[p].size();
            sendIDs[p].clear();
        }

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        recvTNElements.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);

        par::Mpi_Alltoallv(&(*(sendTNElements.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvTNElements.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);

        //2. generate recvTNElement keys and send the local nodes to the owner.
        std::vector<ot::Key> recvTNElem_keys;
        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=m_uiRecvNodeOffset[p];e<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p]);e++)
            {
                itKK=recvTNElem_keys.emplace(recvTNElem_keys.end(),ot::Key(recvTNElements[e]));
                itKK->addOwner(p);
            }
        }

        SFC::seqSearch::SFC_treeSearch(&(*(recvTNElem_keys.begin())), &(*(m_uiAllElements.begin())), 0, recvTNElem_keys.size(),m_uiElementLocalBegin, m_uiElementLocalEnd, m_uiMaxDepth, m_uiMaxDepth, ROOT_ROTATION);
        // local nodes is done.



        //for bdy nodes.
        //3. create local bdy nodes.
        std::vector<ot::Node> localBdy;
        std::vector<ot::SearchKey> localBdy1; // original local Bdy1.

        for(unsigned int e=m_uiNodeLocalBegin;e<m_uiNodeLocalEnd;e++)
        {
            dg2eijk(m_uiCG2DG[e],ownerID,ii_x,jj_y,kk_z);
            if(getDIROfANode(ii_x,jj_y,kk_z)!=OCT_DIR_INTERNAL)
            {
                x = m_uiAllElements[ownerID].getX();
                y = m_uiAllElements[ownerID].getY();
                z = m_uiAllElements[ownerID].getZ();
                sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
                assert(sz%m_uiElementOrder==0);

                itSK=localBdy1.emplace(localBdy1.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                itSK->addOwner(e);
            }

        }

        localBdy.resize(localBdy1.size());
        for(unsigned int e=0;e<localBdy1.size();e++)
        {
            localBdy[e]=localBdy1[e];
            localBdy[e].setOwner(rank);
        }



        SFC::seqSort::SFC_treeSort(&(*(localBdy1.begin())),localBdy1.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        //3.a par sort of local bdy nodes.
#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        SFC::parSort::SFC_treeSort(localBdy,tmpNN,tmpNN,tmpNN,0.1,m_uiMaxDepth+1,rootNN,ROOT_ROTATION,1,TS_SORT_ONLY,2,comm);

#ifdef DEBUG_SM
        m_uiMaxDepth++;
        treeNodesTovtk(localBdy,rank,"localBdy");
        m_uiMaxDepth--;
#endif



#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" par::sort local_bdy nodes (max) (s): "<<t_stat_g<<std::endl;
#endif
        m_uiMaxDepth++;
        assert(par::test::isUniqueAndSorted(localBdy,comm));
        m_uiMaxDepth--;

        ot::TreeNode minMax[2];
        minMax[0]=localBdy.front();
        minMax[1]=localBdy.back();

        m_uiSplitterNodes=new ot::TreeNode[2*npes];
        par::Mpi_Allgather(minMax,m_uiSplitterNodes,2,comm);



        std::vector<bool > g1Visited;
        g1Visited.resize(m_uiCG2DG.size(),false);

        //4. Generated allocated nodes with boundary allocation.

        std::vector<ot::Node> allocatedBdy; // allocated nodes where the owner of the nodes are undecided.
        std::vector<SearchKey> allocatedNodes; // actual allocated nodes.
#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        for(unsigned int ele=0;ele<m_uiGhostElementRound1Index.size();ele++) {
            for (unsigned int node = 0; node < m_uiNpE; node++) {
                nodeIndex = m_uiE2NMapping_CG[m_uiGhostElementRound1Index[ele] * m_uiNpE + node];
                if ((!(nodeIndex >= m_uiNodeLocalBegin && nodeIndex < m_uiNodeLocalEnd)) && (!g1Visited[nodeIndex])) {
                    assert(nodeIndex<g1Visited.size());
                    dg2eijk(m_uiCG2DG[nodeIndex],ownerID,ii_x,jj_y,kk_z);
                    nodeFlag=getDIROfANode(ii_x,jj_y,kk_z);
                    x = m_uiAllElements[ownerID].getX();
                    y = m_uiAllElements[ownerID].getY();
                    z = m_uiAllElements[ownerID].getZ();
                    sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
                    itSK=allocatedNodes.emplace(allocatedNodes.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                    itSK->addOwner(nodeIndex);

                    if(nodeFlag != OCT_DIR_INTERNAL)
                    {   // for internal nodes we can directly determine the ownership.
                        itNN=allocatedBdy.emplace(allocatedBdy.end(),SearchKey((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                        itNN->setOwner(rank);
                    }
                    g1Visited[nodeIndex]= true;
                }
            }
        }

        const unsigned int totAllocated=allocatedNodes.size();
        SFC::seqSort::SFC_treeSort(&(*(allocatedNodes.begin())),allocatedNodes.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

#ifdef DEBUG_SM
        m_uiMaxDepth++;
        treeNodesTovtk(allocatedBdy,rank,"allocatedBdy");
        m_uiMaxDepth--;
#endif
        for(unsigned int p=0;p<2*npes;p++)
            allocatedBdy.emplace(allocatedBdy.end(),m_uiSplitterNodes[p]);


        SFC::seqSort::SFC_treeSort(&(*(allocatedBdy.begin())),allocatedBdy.size(),tmpNN,tmpNN,tmpNN,m_uiMaxDepth+1,m_uiMaxDepth+1,rootNN,ROOT_ROTATION,1,TS_SORT_ONLY);

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" allocated node generation + sort time (max) (s): "<<t_stat_g<<std::endl;
#endif
        // remove duplicates for allocated bdy octants.
        tmpNN.clear();
        for(unsigned int e=0;e<(allocatedBdy.size());e++)
        {
            skip=1;
            tmpNode=allocatedBdy[e];
            while(((e+skip)<allocatedBdy.size()) && (allocatedBdy[e]==allocatedBdy[e+skip]))
            {
                if(allocatedBdy[e+skip].getOwner()>=0)tmpNode.setOwner(allocatedBdy[e+skip].getOwner());
                skip++;
            }
            tmpNN.push_back(tmpNode);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(allocatedBdy,tmpNN);
        tmpNN.clear();

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(allocatedBdy));
        m_uiMaxDepth--;


        std::vector<ot::Key> splitterNode_keys;
        splitterNode_keys.resize(2*npes);
        std::vector<unsigned int > nodeSplitterID;
        nodeSplitterID.resize(2*npes);

        for(unsigned int p=0;p<2*npes;p++)
            splitterNode_keys[p]=ot::Key(m_uiSplitterNodes[p]);

        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(splitterNode_keys));
        m_uiMaxDepth--;


        m_uiMaxDepth++;
        searchKeys(splitterNode_keys,allocatedBdy);
        m_uiMaxDepth--;

        for(unsigned int p=0;p<2*npes;p++)
        {
            assert(splitterNode_keys[p].getFlag() & OCT_FOUND);
            assert(allocatedBdy[splitterNode_keys[p].getSearchResult()]==splitterNode_keys[p]);
            nodeSplitterID[p]=splitterNode_keys[p].getSearchResult();
            assert(nodeSplitterID[p]<allocatedBdy.size());
        }

        // 5. send bdy allocated nodes to the correct processor (this is according to the par::sort of localBdy).
        sendNNodal.clear();
        unsigned int sBegin;
        unsigned int sEnd;
        for(unsigned int p=0;p<npes;p++)
        {
            sBegin=nodeSplitterID[2*p];
            sEnd=nodeSplitterID[2*p+1]+1;
            m_uiSendNodeCount[p]=0;
            for(unsigned int e=sBegin;e<sEnd;e++)
            {
                if(allocatedBdy[e].getOwner()>=0)
                {
                    sendNNodal.push_back(allocatedBdy[e]);
                    m_uiSendNodeCount[p]++;
                }

            }

        }


        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        recvNNodal.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);
#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        par::Mpi_Alltoallv(&(*(sendNNodal.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNNodal.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" R1 all2all (max) (s): "<<t_stat_g<<std::endl;
#endif

        std::vector<ot::Key> recvNodal_keys;
        recvNodal_keys.resize(recvNNodal.size());
        for(unsigned int e=0;e<recvNNodal.size();e++)
        {
            recvNodal_keys[e]=ot::Key(recvNNodal[e]);
            recvNodal_keys[e].addOwner(e);
        }

#ifdef DEBUG_SM
        m_uiMaxDepth++;
        treeNodesTovtk(recvNodal_keys,rank,"recvNNodal");
        m_uiMaxDepth--;
#endif


#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        SFC::seqSearch::SFC_treeSearch(&(*(recvNodal_keys.begin())), &(*(localBdy.begin())), 0, recvNodal_keys.size(),0, localBdy.size(), m_uiMaxDepth + 1, m_uiMaxDepth + 1, ROOT_ROTATION);

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" R1 seq::search max (s): "<<t_stat_g<<std::endl;
#endif

        sendNNodal.clear();
        unsigned int result;

        for(unsigned int p=0;p<npes;p++)
            sendIDs[p].clear();

        // all the recv nodes should be found in localBdy
#ifdef DEBUG_SM
        std::vector<ot::TreeNode> missingKeys;
#endif

        for(unsigned int e=0;e<recvNodal_keys.size();e++)
        {
            if(!(recvNodal_keys[e].getFlag() & OCT_FOUND))
            {

#ifdef DEBUG_SM
                unsigned int  found=std::find(localBdy.begin(),localBdy.end(),recvNodal_keys[e])-localBdy.begin();
                missingKeys.push_back(recvNodal_keys[e]);
                m_uiMaxDepth++;
                std::cout<<"rank: "<<rank<<" recvNodalKey : "<<recvNodal_keys[e]<<" not found: status: "<<(recvNodal_keys[e]>=m_uiSplitterNodes[32] && recvNodal_keys[e]<=m_uiSplitterNodes[33])<<"found: "<<found<<" of "<<localBdy.size()<<std::endl;
                m_uiMaxDepth--;
                continue;
#endif
                std::cout<<"rank: "<<m_uiActiveRank<<"SM4 Error: Allocated key<"<<recvNodal_keys[e]<<" not found. "<<std::endl;
                exit(0);


            }

            assert(recvNodal_keys[e].getFlag() & OCT_FOUND);
            result=recvNodal_keys[e].getSearchResult();
            ownerID=recvNodal_keys[e].getOwnerList()->front();
            sendIDs[localBdy[result].getOwner()].push_back(ownerID);
        }

#ifdef DEBUG_SM
        m_uiMaxDepth++;
        if(missingKeys.size()) treeNodesTovtk(missingKeys,rank,"missingKeys");
        m_uiMaxDepth--;
#endif

        sendNNodal.clear();
        // now send the recv nodes to actual correct proceesor (based on localBdy1)
        for(unsigned int p=0;p<npes;p++)
        {
            for(unsigned int e=0;e<sendIDs[p].size();e++)
                sendNNodal.push_back(recvNNodal[sendIDs[p][e]]);

            m_uiSendNodeCount[p]=sendIDs[p].size();
            sendIDs[p].clear();

        }

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        recvNNodal.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);
#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        par::Mpi_Alltoallv(&(*(sendNNodal.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvNNodal.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" R2 all2all max (s): "<<t_stat_g<<std::endl;
#endif
        recvNodal_keys.resize(recvNNodal.size());
        for(unsigned int e=0;e<recvNNodal.size();e++)
        {
            recvNodal_keys[e]=ot::Key(recvNNodal[e]);
            recvNodal_keys[e].addOwner(e);
        }

#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        SFC::seqSearch::SFC_treeSearch(&(*(recvNodal_keys.begin())), &(*(localBdy1.begin())), 0, recvNodal_keys.size(),0, localBdy1.size(), m_uiMaxDepth + 1, m_uiMaxDepth + 1, ROOT_ROTATION);

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" R2 seq::search max (s): "<<t_stat_g<<std::endl;
#endif

        sendNNodal.clear();


        m_uiScatterMapActualNodeSend.clear();
        m_uiScatterMapActualNodeRecv.clear();
        sendTNElements.clear();

        for(unsigned int p=0;p<npes;p++)
            sendIDs[p].clear();

#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        // put internal nodes to the scattermap.
        for(unsigned int e=0;e<recvTNElem_keys.size();e++)
        {
            assert(recvTNElem_keys[e].getFlag() & OCT_FOUND);
            result=recvTNElem_keys[e].getSearchResult();
            ownerID=recvTNElem_keys[e].getOwnerList()->front();
            assert(ownerID!=rank);
            assert(result>=m_uiElementLocalBegin && result<=m_uiElementLocalEnd);

            for(unsigned int k=1;k<m_uiElementOrder;k++)
               for(unsigned int j=1;j<m_uiElementOrder;j++)
                  for(unsigned int i=1;i<m_uiElementOrder;i++)
                   sendIDs[ownerID].push_back(m_uiE2NMapping_CG[result*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]);

        }



        // all the recv nodal should be found in localBdy1
        for(unsigned int e=0;e<recvNodal_keys.size();e++)
        {
            assert(recvNodal_keys[e].getFlag() & OCT_FOUND);
            result=recvNodal_keys[e].getSearchResult();
            ownerID=recvNNodal[recvNodal_keys[e].getOwnerList()->front()].getOwner();
            assert(ownerID!=rank);
            sendIDs[ownerID].push_back(localBdy1[result].getOwner());

        }


        for(unsigned int p=0;p<npes;p++)
        {
            //Note: This should be unique.
            std::sort(sendIDs[p].begin(),sendIDs[p].end());
            /*if(!rank) std::cout<<" rank : "<<rank<<" p: "<<p<<" bf rmd : "<<sendIDs[p].size()<<std::endl;
            sendIDs[p].erase(std::unique(sendIDs[p].begin(),sendIDs[p].end()),sendIDs[p].end());
            if(!rank) std::cout<<" rank : "<<rank<<" p: "<<p<<" af rmd : "<<sendIDs[p].size()<<std::endl;*/
            for(unsigned int e=0;e<sendIDs[p].size();e++)
            {
                //if(!rank) std::cout<<" e: "<<e<<" val: "<<sendIDs[p][e]<<std::endl;
                dg2eijk(m_uiCG2DG[sendIDs[p][e]],ownerID,ii_x,jj_y,kk_z);
                x = m_uiAllElements[ownerID].getX();
                y = m_uiAllElements[ownerID].getY();
                z = m_uiAllElements[ownerID].getZ();
                sz = 1u << (m_uiMaxDepth - m_uiAllElements[ownerID].getLevel());
                sendTNElements.emplace(sendTNElements.end(),ot::TreeNode((x + ii_x * sz/m_uiElementOrder), (y + jj_y * sz/m_uiElementOrder), (z + kk_z * sz/m_uiElementOrder), m_uiMaxDepth+1,m_uiDim, m_uiMaxDepth+1));
                m_uiScatterMapActualNodeSend.push_back(sendIDs[p][e]);

            }

            m_uiSendNodeCount[p]=sendIDs[p].size();
            sendIDs[p].clear();

        }

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" send sm setup max (s): "<<t_stat_g<<std::endl;
#endif

        par::Mpi_Alltoall(&(*(m_uiSendNodeCount.begin())),&(*(m_uiRecvNodeCount.begin())),1,comm);

        m_uiSendNodeOffset[0]=0;
        m_uiRecvNodeOffset[0]=0;

        omp_par::scan(&(*(m_uiSendNodeCount.begin())),&(*(m_uiSendNodeOffset.begin())),npes);
        omp_par::scan(&(*(m_uiRecvNodeCount.begin())),&(*(m_uiRecvNodeOffset.begin())),npes);

        recvTNElements.resize(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1]);

#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        par::Mpi_Alltoallv(&(*(sendTNElements.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(recvTNElements.begin())),(int *)(&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),comm);

#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" final all2all max (s): "<<t_stat_g<<std::endl;
#endif
        if(totAllocated!=recvTNElements.size())
        {
            std::cout<<"rank: "<<rank<<"[SM Error]: allocated nodes: "<<totAllocated<<" received nodes: "<<(m_uiRecvNodeOffset[npes-1]+m_uiRecvNodeCount[npes-1])<<std::endl;
            exit(0);
        }

        std::vector<SearchKey> recvNodeSKeys;
        recvNodeSKeys.resize(recvTNElements.size());
        for(unsigned int e=0;e<recvTNElements.size();e++)
        {
            recvNodeSKeys[e]=SearchKey(recvTNElements[e]);
            recvNodeSKeys[e].addOwner(e);
        }

#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif
        SFC::seqSort::SFC_treeSort(&(*(recvNodeSKeys.begin())),recvNodeSKeys.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth+1,m_uiMaxDepth+1,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);


#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<" recv node seq::sort max (s): "<<t_stat_g<<std::endl;
#endif
        m_uiMaxDepth++;
        assert(seq::test::isUniqueAndSorted(recvNodeSKeys));
        m_uiMaxDepth--;

#ifdef PROFILE_SM
        t1=MPI_Wtime();
#endif

        m_uiScatterMapActualNodeRecv.resize(recvTNElements.size());
        unsigned int alCount=0;
        for(int e=0;e<recvNodeSKeys.size();e++)
        {
            if(allocatedNodes[alCount].getOwner()<0)
            {
                e--;
                alCount++;
                continue;
            }

            if(allocatedNodes[alCount]!=recvNodeSKeys[e]) {
                std::cout << "rank: " << rank << " allocated[" << alCount << "]: " << allocatedNodes[alCount]<< " received[" << e << "]: " << recvNodeSKeys[e] << std::endl;
                exit(0);
            }

            m_uiScatterMapActualNodeRecv[recvNodeSKeys[e].getOwner()]=allocatedNodes[alCount].getOwner();
            alCount++;

        }


#ifdef PROFILE_SM
        t2=MPI_Wtime();
        t_stat=t2-t1;
        par::Mpi_Reduce(&t_stat,&t_stat_g,1,MPI_MAX,0,comm);
        if(!rank) std::cout<<"recv sm time  max (s): "<<t_stat_g<<std::endl;
#endif

        delete [] sendIDs;




    }


    void Mesh::performBlocksSetup()
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        // assumes that E2E and E2N mapping is done and m_uiAllElements should be sorted otherwise this will chnage the order of elements in m_uiAllElements.
        assert(seq::test::isUniqueAndSorted(m_uiAllElements));
        octree2BlockDecomposition(m_uiAllElements,m_uiLocalBlockList,m_uiMaxDepth,m_uiDmin,m_uiDmax,m_uiElementLocalBegin,m_uiElementLocalEnd,m_uiElementOrder);
        assert(ot::test::isBlockListValid(m_uiAllElements,m_uiLocalBlockList,m_uiDmin,m_uiDmax,m_uiElementLocalBegin,m_uiElementLocalEnd));

        std::vector<DendroIntL> blkSz;
        std::vector<DendroIntL> blkSzOffset;

        blkSz.resize(m_uiLocalBlockList.size());
        blkSzOffset.resize(m_uiLocalBlockList.size());

        for(unsigned int k=0;k<m_uiLocalBlockList.size();k++)
            blkSz[k]=m_uiLocalBlockList[k].get1DArraySize()*m_uiLocalBlockList[k].get1DArraySize()*m_uiLocalBlockList[k].get1DArraySize();

        blkSzOffset[0]=0;
        omp_par::scan(&(*(blkSz.begin())),&(*(blkSzOffset.begin())),m_uiLocalBlockList.size());

        for(unsigned int k=0;k<m_uiLocalBlockList.size();k++)
            m_uiLocalBlockList[k].setOffset(blkSzOffset[k]);

        m_uiUnZippedVecSz=blkSzOffset[m_uiLocalBlockList.size()-1]+blkSz.back();

        const unsigned int dmin=0;
        const unsigned int dmax=1u<<(m_uiMaxDepth);
        ot::TreeNode blkNode;

        std::vector<ot::Key> blkKeys;
        std::vector<ot::SearchKey> blkSkeys;

        unsigned int sz;
        unsigned int regLev;
        unsigned int blkElem_1D;

        std::vector<unsigned int >* ownerList;
        std::vector<unsigned int >* directionList;
        unsigned int result;


        for(unsigned int e=0;e<m_uiLocalBlockList.size();e++)
        {


            blkNode=m_uiLocalBlockList[e].getBlockNode();

            if(blkNode.minX()==dmin)
            {
                blkNode.setFlag(((blkNode.getFlag())|((1u<<(OCT_DIR_LEFT+NUM_LEVEL_BITS))|blkNode.getLevel())));
                assert((blkNode.getFlag()>>NUM_LEVEL_BITS)&(1u<<OCT_DIR_LEFT));
            }

            if(blkNode.minY()==dmin)
            {
                blkNode.setFlag(((blkNode.getFlag())|((1u<<(OCT_DIR_DOWN+NUM_LEVEL_BITS))|blkNode.getLevel())));
                assert((blkNode.getFlag()>>NUM_LEVEL_BITS)&(1u<<OCT_DIR_DOWN));
            }

            if(blkNode.minZ()==dmin)
            {
                blkNode.setFlag(((blkNode.getFlag())|((1u<<(OCT_DIR_BACK+NUM_LEVEL_BITS))|blkNode.getLevel())));
                assert((blkNode.getFlag()>>NUM_LEVEL_BITS)&(1u<<OCT_DIR_BACK));
            }


            if(blkNode.maxX()==dmax)
            {
                blkNode.setFlag(((blkNode.getFlag())|((1u<<(OCT_DIR_RIGHT+NUM_LEVEL_BITS))|blkNode.getLevel())));
                assert((blkNode.getFlag()>>NUM_LEVEL_BITS)&(1u<<OCT_DIR_RIGHT));
            }

            if(blkNode.maxY()==dmax)
            {
                blkNode.setFlag(((blkNode.getFlag())|((1u<<(OCT_DIR_UP+NUM_LEVEL_BITS))|blkNode.getLevel())));
                assert((blkNode.getFlag()>>NUM_LEVEL_BITS)&(1u<<OCT_DIR_UP));
            }

            if(blkNode.maxZ()==dmax)
            {
                blkNode.setFlag(((blkNode.getFlag())|((1u<<(OCT_DIR_FRONT+NUM_LEVEL_BITS))|blkNode.getLevel())));
                assert((blkNode.getFlag()>>NUM_LEVEL_BITS)&(1u<<OCT_DIR_FRONT));
            }

            assert(blkNode.getLevel()==m_uiLocalBlockList[e].getBlockNode().getLevel());
            m_uiLocalBlockList[e].setBlkNodeFlag(blkNode.getFlag());

            regLev=m_uiLocalBlockList[e].getRegularGridLev();
            sz=1u<<(m_uiMaxDepth-regLev);

            m_uiLocalBlockList[e].initializeBlkDiagMap(LOOK_UP_TABLE_DEFAULT);
            m_uiLocalBlockList[e].initializeBlkVertexMap(LOOK_UP_TABLE_DEFAULT);

            blkSkeys.clear();
            blkKeys.clear();
            generateBlkEdgeSKeys(m_uiLocalBlockList[e],blkSkeys);
            generateBlkVertexSKeys(m_uiLocalBlockList[e],blkSkeys);
            mergeKeys(blkSkeys,blkKeys);
            blkSkeys.clear();
            SFC::seqSearch::SFC_treeSearch(&(*(blkKeys.begin())),&(*(m_uiAllElements.begin())),0,blkKeys.size(),0,m_uiAllElements.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);


            for(unsigned int i=0;i<blkKeys.size();i++)
            {
                assert(blkKeys[i].getFlag()& OCT_FOUND);
                if(!(blkKeys[i].getFlag()& OCT_FOUND)) {std::cout<<RED<<"block diagonal key not found"<<NRM<<std::endl;}
                ownerList=blkKeys[i].getOwnerList();
                directionList=blkKeys[i].getStencilIndexDirectionList();
                result=blkKeys[i].getSearchResult();


                assert(ownerList->size()==directionList->size());
                for(unsigned int w=0;w<ownerList->size();w++)
                {
                    if((*directionList)[w]<VERTEX_OFFSET)
                    {
                        assert((*directionList)[w]>=EDGE_OFFSET);
                        m_uiLocalBlockList[e].setBlk2DiagMap((*ownerList)[w],((*directionList)[w]-EDGE_OFFSET),result);
                    }
                    else
                    {// this is an vertex neighbour.
                        assert((*directionList)[w]>=VERTEX_OFFSET);
                        m_uiLocalBlockList[e].setBlk2VertexMap(((*directionList)[w]-VERTEX_OFFSET),result);
                    }

                }


            }



        }





        double * zippVec=createVector<double>(NAN);
        for(unsigned int node=m_uiNodeLocalBegin;node<m_uiNodeLocalEnd;node++)
            zippVec[node]=0;

        double * unzipVec=createUnZippedVector<double>(NAN);

        // need to do this without performing the ghost exchange.
        unzip(zippVec,unzipVec);

        unsigned int offset;
        unsigned int lx,ly,lz;
        unsigned int pW;
        bool haveanyGhost=false;
        unsigned int ib,ie,jb,je,kb,ke;
        unsigned int bflag;

        for(unsigned int e=0;e<m_uiLocalBlockList.size();e++)
        {

            haveanyGhost=false;

            offset=m_uiLocalBlockList[e].getOffset();
            lx=m_uiLocalBlockList[e].getAllocationSzX();
            ly=m_uiLocalBlockList[e].getAllocationSzY();
            lz=m_uiLocalBlockList[e].getAllocationSzZ();

            pW=m_uiLocalBlockList[e].get1DPadWidth();
            bflag=m_uiLocalBlockList[e].getBlkNodeFlag();

            /*for(unsigned int k=3;k<(nz-3);k++)
                for(unsigned int j=3;j<(ny-3);j++)
                    for(unsigned int i=3;i<(nx-3);i++)
                        std::cout<<m_uiActiveRank<<": unzip value: "<<unzipVec[offset+k*(nz*ny)+j*(ny)+i]<<std::endl;*/

            //internal
            ib=pW;
            ie=lx-pW;

            jb=pW;
            je=ly-pW;

            kb=pW;
            ke=lz-pW;

            for(unsigned int k=kb; k < ke; k++)
                for(unsigned int j=jb;j < je; j++)
                    for(unsigned int i=ib;i < ie; i++)
                        if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                            haveanyGhost=true;


            if(!(bflag & (1u<<OCT_DIR_LEFT)))
            {
                ib=0;
                ie=pW;
                jb=pW;
                je=ly-pW;
                kb=pW;
                ke=lz-pW;

                for(unsigned int k=kb; k < ke; k++)
                    for(unsigned int j=jb;j < je; j++)
                        for(unsigned int i=ib;i < ie; i++)
                            if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                haveanyGhost=true;


                if(!(bflag & (1u<<OCT_DIR_DOWN)))
                {
                    ib=1;
                    ie=pW;
                    jb=1;
                    je=pW;
                    kb=pW;
                    ke=lz-pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }

                if(!(bflag & (1u<<OCT_DIR_UP)))
                {
                    ib=1;
                    ie=pW;
                    jb=ly-pW;
                    je=ly-1;
                    kb=pW;
                    ke=lz-pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }


                if(!(bflag & (1u<<OCT_DIR_BACK)))
                {
                    ib=1;
                    ie=pW;
                    jb=pW;
                    je=ly-pW;
                    kb=1;
                    ke=pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }


                if(!(bflag & (1u<<OCT_DIR_FRONT)))
                {
                    ib=1;
                    ie=pW;
                    jb=pW;
                    je=ly-pW;

                    kb=lz-pW;
                    ke=lz-1;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }




            }


            if(!(bflag & (1u<<OCT_DIR_RIGHT)))
            {


                ib=lx-pW;
                ie=lx;
                jb=pW;
                je=ly-pW;
                kb=pW;
                ke=lz-pW;

                for(unsigned int k=kb; k < ke; k++)
                    for(unsigned int j=jb;j < je; j++)
                        for(unsigned int i=ib;i < ie; i++)
                            if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                haveanyGhost=true;



                if(!(bflag & (1u<<OCT_DIR_DOWN)))
                {
                    ib=lx-pW;
                    ie=lx-1;
                    jb=1;
                    je=pW;
                    kb=pW;
                    ke=lz-pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }

                if(!(bflag & (1u<<OCT_DIR_UP)))
                {
                    ib=lx-pW;
                    ie=lx-1;
                    jb=ly-pW;
                    je=ly-1;
                    kb=pW;
                    ke=lz-pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }


                if(!(bflag & (1u<<OCT_DIR_BACK)))
                {
                    ib=lx-pW;
                    ie=lx-1;
                    jb=pW;
                    je=ly-pW;
                    kb=1;
                    ke=pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }


                if(!(bflag & (1u<<OCT_DIR_FRONT)))
                {
                    ib=lx-pW;
                    ie=lx-1;
                    jb=pW;
                    je=ly-pW;

                    kb=lz-pW;
                    ke=lz-1;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }


            }



            if(!(bflag & (1u<<OCT_DIR_DOWN)))
            {


                ib=pW;
                ie=lx-pW;
                jb=0;
                je=pW;
                kb=pW;
                ke=lz-pW;

                for(unsigned int k=kb; k < ke; k++)
                    for(unsigned int j=jb;j < je; j++)
                        for(unsigned int i=ib;i < ie; i++)
                            if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                haveanyGhost=true;



                if(!(bflag & (1u<<OCT_DIR_BACK)))
                {
                    ib=pW;
                    ie=lx-pW;
                    jb=1;
                    je=pW;
                    kb=1;
                    ke=pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }


                if(!(bflag & (1u<<OCT_DIR_FRONT)))
                {
                    ib=pW;
                    ie=lx-pW;
                    jb=1;
                    je=pW;

                    kb=lz-pW;
                    ke=lz-1;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }



            }


            if(!(bflag & (1u<<OCT_DIR_UP)))
            {


                ib=pW;
                ie=lx-pW;
                jb=ly-pW;
                je=ly;
                kb=pW;
                ke=lz-pW;

                for(unsigned int k=kb; k < ke; k++)
                    for(unsigned int j=jb;j < je; j++)
                        for(unsigned int i=ib;i < ie; i++)
                            if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                haveanyGhost=true;



                if(!(bflag & (1u<<OCT_DIR_BACK)))
                {
                    ib=pW;
                    ie=lx-pW;
                    jb=ly-pW;
                    je=ly-1;
                    kb=1;
                    ke=pW;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }


                if(!(bflag & (1u<<OCT_DIR_FRONT)))
                {
                    ib=pW;
                    ie=lx-pW;
                    jb=ly-pW;
                    je=ly-1;

                    kb=lz-pW;
                    ke=lz-1;

                    for(unsigned int k=kb; k < ke; k++)
                        for(unsigned int j=jb;j < je; j++)
                            for(unsigned int i=ib;i < ie; i++)
                                if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                    haveanyGhost=true;

                }



            }


            if(!(bflag & (1u<<OCT_DIR_BACK)))
            {


                ib=pW;
                ie=lx-pW;
                jb=pW;
                je=ly-pW;
                kb=0;
                ke=pW;

                for(unsigned int k=kb; k < ke; k++)
                    for(unsigned int j=jb;j < je; j++)
                        for(unsigned int i=ib;i < ie; i++)
                            if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                haveanyGhost=true;

            }


            if(!(bflag & (1u<<OCT_DIR_FRONT)))
            {


                ib=pW;
                ie=lx-pW;
                jb=pW;
                je=ly-pW;
                kb=lz-pW;
                ke=lz;

                for(unsigned int k=kb; k < ke; k++)
                    for(unsigned int j=jb;j < je; j++)
                        for(unsigned int i=ib;i < ie; i++)
                            if(std::isnan(unzipVec[offset+k*(ly*lx)+j*(lx)+i]))
                                haveanyGhost=true;

            }


           // note: add the vertex corner paddings when you have done them.

           m_uiLocalBlockList[e].setIsInternal(!haveanyGhost);

           /*if(!haveanyGhost)
              std::cout<<"rank: "<<m_uiActiveRank<<" internal block: "<<m_uiLocalBlockList[e].getBlockNode()<<std::endl;*/


        }



        delete [] zippVec;
        delete [] unzipVec;


    }


    bool Mesh::isEdgeHanging(unsigned int elementId,unsigned int edgeId,unsigned int &cnum) const
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return false;

        unsigned int nodeLookUp_DG;
        bool isHanging=false;
        assert(elementId>=0 && elementId<m_uiAllElements.size());

        switch (edgeId)
        {

            case OCT_DIR_LEFT_DOWN:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(0)*(m_uiElementOrder+1)+(0)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minZ()==m_uiAllElements[elementId].minZ())
                    cnum=0;
                else
                {
                   assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxZ()==m_uiAllElements[elementId].maxZ());
                   cnum=1;
                }
                break;

            case OCT_DIR_LEFT_UP:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(m_uiElementOrder)*(m_uiElementOrder+1)+(0)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minZ()==m_uiAllElements[elementId].minZ())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxZ()==m_uiAllElements[elementId].maxZ());
                    cnum=1;
                }
                break;

            case OCT_DIR_LEFT_BACK:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(0)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+(0)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minY()==m_uiAllElements[elementId].minY())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxY()==m_uiAllElements[elementId].maxY());
                    cnum=1;
                }
                break;

            case OCT_DIR_LEFT_FRONT:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(m_uiElementOrder)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+(0)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minY()==m_uiAllElements[elementId].minY())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxY()==m_uiAllElements[elementId].maxY());
                    cnum=1;
                }
                break;
                break;

            case OCT_DIR_RIGHT_DOWN:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(0)*(m_uiElementOrder+1)+(m_uiElementOrder)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minZ()==m_uiAllElements[elementId].minZ())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxZ()==m_uiAllElements[elementId].maxZ());
                    cnum=1;
                }
                break;

            case OCT_DIR_RIGHT_UP:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(m_uiElementOrder)*(m_uiElementOrder+1)+(m_uiElementOrder)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minZ()==m_uiAllElements[elementId].minZ())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxZ()==m_uiAllElements[elementId].maxZ());
                    cnum=1;
                }
                break;

            case OCT_DIR_RIGHT_BACK:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(0)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+(m_uiElementOrder)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minY()==m_uiAllElements[elementId].minY())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxY()==m_uiAllElements[elementId].maxY());
                    cnum=1;
                }
                break;

            case OCT_DIR_RIGHT_FRONT:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(m_uiElementOrder)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+(m_uiElementOrder)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minY()==m_uiAllElements[elementId].minY())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxY()==m_uiAllElements[elementId].maxY());
                    cnum=1;
                }
                break;

            case OCT_DIR_DOWN_BACK:
                 nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(0)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(0)*(m_uiElementOrder+1)+(m_uiElementOrder>>1u)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minX()==m_uiAllElements[elementId].minX())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxX()==m_uiAllElements[elementId].maxX());
                    cnum=1;
                }
                 break;

            case OCT_DIR_DOWN_FRONT:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(m_uiElementOrder)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(0)*(m_uiElementOrder+1)+(m_uiElementOrder>>1u)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minX()==m_uiAllElements[elementId].minX())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxX()==m_uiAllElements[elementId].maxX());
                    cnum=1;
                }
                break;


            case OCT_DIR_UP_BACK:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(0)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(m_uiElementOrder)*(m_uiElementOrder+1)+(m_uiElementOrder>>1u)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minX()==m_uiAllElements[elementId].minX())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxX()==m_uiAllElements[elementId].maxX());
                    cnum=1;
                }
                break;
            case OCT_DIR_UP_FRONT:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(m_uiElementOrder)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(m_uiElementOrder)*(m_uiElementOrder+1)+(m_uiElementOrder>>1u)];
                if(m_uiAllElements[nodeLookUp_DG/m_uiNpE].minX()==m_uiAllElements[elementId].minX())
                    cnum=0;
                else
                {
                    assert(m_uiAllElements[nodeLookUp_DG/m_uiNpE].maxX()==m_uiAllElements[elementId].maxX());
                    cnum=1;
                }
                break;


        }

        isHanging=m_uiAllElements[nodeLookUp_DG/m_uiNpE].getLevel()<m_uiAllElements[elementId].getLevel();
        return isHanging;


    }


    bool Mesh::isFaceHanging(unsigned int elementId, unsigned int faceId,unsigned int &cnum) const
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return false;

        unsigned int nodeLookUp_DG;
        bool isHanging=false;
        unsigned int ownerID;
        unsigned int mid_bit;
        assert(elementId>=0 && elementId<m_uiAllElements.size());
        unsigned int pSz,cSz;

        switch (faceId)
        {

            case OCT_DIR_LEFT:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+(0)];
                ownerID=nodeLookUp_DG/m_uiNpE;
                mid_bit=m_uiMaxDepth - m_uiAllElements[ownerID].getLevel()-1;
                cnum=( (((((m_uiAllElements[elementId].getZ())-(m_uiAllElements[ownerID].getZ())) >> mid_bit) & 1u) << 1u) | ((((m_uiAllElements[elementId].getY())-(m_uiAllElements[ownerID].getY())) >>mid_bit) & 1u));
                break;

            case OCT_DIR_RIGHT:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+(m_uiElementOrder)];
                ownerID=nodeLookUp_DG/m_uiNpE;
                mid_bit=m_uiMaxDepth - m_uiAllElements[ownerID].getLevel()-1;
                cnum=( (((((m_uiAllElements[elementId].getZ())-(m_uiAllElements[ownerID].getZ())) >> mid_bit) & 1u) << 1u) | ((((m_uiAllElements[elementId].getY())-(m_uiAllElements[ownerID].getY())) >>mid_bit) & 1u));

                break;

            case OCT_DIR_DOWN:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(0)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))];
                ownerID=nodeLookUp_DG/m_uiNpE;
                mid_bit=m_uiMaxDepth - m_uiAllElements[ownerID].getLevel()-1;
                cnum=( (((((m_uiAllElements[elementId].getZ())-(m_uiAllElements[ownerID].getZ())) >> mid_bit) & 1u) << 1u) | ((((m_uiAllElements[elementId].getX())-(m_uiAllElements[ownerID].getX())) >>mid_bit) & 1u));
                break;

            case OCT_DIR_UP:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(m_uiElementOrder)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))];
                ownerID=nodeLookUp_DG/m_uiNpE;
                mid_bit=m_uiMaxDepth - m_uiAllElements[ownerID].getLevel()-1;
                cnum=( (((((m_uiAllElements[elementId].getZ())-(m_uiAllElements[ownerID].getZ())) >> mid_bit) & 1u) << 1u) | ((((m_uiAllElements[elementId].getX())-(m_uiAllElements[ownerID].getX())) >>mid_bit) & 1u));
                break;

            case OCT_DIR_BACK:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(0)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))];
                ownerID=nodeLookUp_DG/m_uiNpE;
                mid_bit=m_uiMaxDepth - m_uiAllElements[ownerID].getLevel()-1;
                cnum=( (((((m_uiAllElements[elementId].getY())-(m_uiAllElements[ownerID].getY())) >> mid_bit) & 1u) << 1u) | ((((m_uiAllElements[elementId].getX())-(m_uiAllElements[ownerID].getX())) >>mid_bit) & 1u));
                break;

            case OCT_DIR_FRONT:
                nodeLookUp_DG=m_uiE2NMapping_DG[elementId*m_uiNpE+(m_uiElementOrder)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))*(m_uiElementOrder+1)+((m_uiElementOrder>>1u))];
                ownerID=nodeLookUp_DG/m_uiNpE;
                mid_bit=m_uiMaxDepth - m_uiAllElements[ownerID].getLevel()-1;
                cnum=( (((((m_uiAllElements[elementId].getY())-(m_uiAllElements[ownerID].getY())) >> mid_bit) & 1u) << 1u) | ((((m_uiAllElements[elementId].getX())-(m_uiAllElements[ownerID].getX())) >>mid_bit) & 1u));
                break;

        }

        isHanging=m_uiAllElements[nodeLookUp_DG/m_uiNpE].getLevel()<m_uiAllElements[elementId].getLevel();
        return isHanging;

    }

    bool Mesh::isNodeHanging(unsigned int eleID,unsigned int ix,unsigned int jy,unsigned int kz)const
    {
        // should not be called if the mesh is not active
        if(!m_uiIsActive) return false;

        return m_uiAllElements[(m_uiE2NMapping_DG[eleID*m_uiNpE+kz*(m_uiElementOrder+1)*(m_uiElementOrder+1)+jy*(m_uiElementOrder+1)+ix]/m_uiNpE)].getLevel()<m_uiAllElements[eleID].getLevel();
    }


    ot::Mesh* Mesh::ReMesh(unsigned int grainSz,double ld_tol,unsigned int sfK)
    {

        std::vector<ot::TreeNode> balOct1; //new balanced octree.

        if(m_uiIsActive)
        {

            // 1. build the unbalanced octree from the remesh flags.
            std::vector<ot::TreeNode> unBalancedOctree;
            unsigned int remeshFlag;
            unsigned int sz;
            for(unsigned int ele=m_uiElementLocalBegin;ele<(m_uiElementLocalEnd);ele++)
            {
                remeshFlag=(m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS);
                assert(m_uiAllElements[ele].getLevel()!=0);
                if(remeshFlag==OCT_SPLIT)
                {
                    m_uiAllElements[ele].addChildren(unBalancedOctree);

                }else if(remeshFlag==OCT_COARSE)
                {
                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                    unBalancedOctree.push_back(m_uiAllElements[ele].getParent());
                    ele=ele+NUM_CHILDREN-1;

                }else
                {
                    assert(remeshFlag==OCT_NO_CHANGE);
                    unBalancedOctree.push_back(m_uiAllElements[ele]);
                }

            }


            std::vector<ot::TreeNode> unBalOctSplitters;
            unBalOctSplitters.resize(m_uiActiveNpes);

            par::Mpi_Allgather(&(unBalancedOctree.back()),&(*(unBalOctSplitters.begin())),1,m_uiCommActive);



            //std::cout<<"rank: "<<m_uiActiveRank<<" unbalnced oct size: "<<unBalancedOctree.size()<<std::endl;
            //treeNodesTovtk(unBalancedOctree,m_uiActiveRank,"unbalncedOctree");
            // unbalancedOctree should be complete. Due to 1 level coarsen and refinement.

            // 2. balanced the octree.

            ot::TreeNode rootNode(m_uiDim,m_uiMaxDepth);
            SFC::parSort::SFC_treeSort(unBalancedOctree,balOct1,balOct1,balOct1,ld_tol,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_BALANCE_OCTREE,sfK,m_uiCommActive);
            assert(par::test::isUniqueAndSorted(balOct1,m_uiCommActive));
            unBalancedOctree.clear();


            if(m_uiActiveNpes==1)
            {
                // sequential case to synchronize flags.

                unsigned int count2=0; // element iterator for the mesh2.
                const ot::TreeNode* pNodes2=&(*(balOct1.begin()));
                ot::TreeNode *pNodes1=&(*(m_uiAllElements.begin()));

                assert(m_uiAllElements[m_uiElementLocalBegin]==balOct1.front() || m_uiAllElements[m_uiElementLocalBegin].getParent()==balOct1.front() || m_uiAllElements[m_uiElementLocalBegin]==balOct1.front().getParent());
                assert(m_uiAllElements[m_uiElementLocalEnd-1]==balOct1.back() || m_uiAllElements[m_uiElementLocalEnd-1].getParent()==balOct1.back() || m_uiAllElements[m_uiElementLocalEnd-1]==balOct1.back().getParent());

                /* while( (count2<balOct1.size()) && (pNodes1[m_uiElementLocalBegin]!=pNodes2[count2]) && (pNodes1[m_uiElementLocalBegin].getParent()!=pNodes2[count2]) && (pNodes1[m_uiElementLocalBegin]!=pNodes2[count2].getParent())) count2++;
                 assert(count2<balOct1.size());*/

                for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                {
                    /* if(!(pNodes1[ele]==pNodes2[count2] || pNodes1[ele].getParent()==pNodes2[count2] || pNodes1[ele]==pNodes2[count2].getParent()))
                     {
                         std::cout<<"rank: "<<m_uiActiveRank<<" ele: "<<ele<<"pNodes1: "<<pNodes1[ele]<<" pNodes2: "<<pNodes2[count2]<<std::endl;
                     }*/
                    assert(count2<balOct1.size());
                    assert(pNodes1[ele]==pNodes2[count2] || pNodes1[ele].getParent()==pNodes2[count2] || pNodes1[ele]==pNodes2[count2].getParent() );

                    if(pNodes1[ele]==pNodes2[count2].getParent())
                    { // old elements have splitted.
                        assert(pNodes2[count2].getParent()==pNodes2[count2+NUM_CHILDREN-1].getParent());
                        pNodes1[ele].setFlag(((OCT_SPLIT<<NUM_LEVEL_BITS)| pNodes1[ele].getLevel()));
                        count2=count2+NUM_CHILDREN;

                    }else if(pNodes1[ele].getParent()==pNodes2[count2])
                    { // old elements have coarsen
                        assert(pNodes1[ele].getParent()==pNodes1[ele+NUM_CHILDREN-1].getParent());
                        pNodes1[ele].setFlag(((OCT_COARSE<<NUM_LEVEL_BITS)| pNodes1[ele].getLevel()));
                        ele=ele+NUM_CHILDREN-1;
                        count2++;

                    }else
                    {
                        assert(pNodes1[ele]==pNodes2[count2]);
                        pNodes1[ele].setFlag(((OCT_NO_CHANGE<<NUM_LEVEL_BITS) | pNodes1[ele].getLevel()));
                        count2++;
                    }

                }

                assert(par::test::isUniqueAndSorted(balOct1,m_uiCommActive));
                //std::cout<<"rank: "<<m_uiActiveRank<<" balanced oct size: "<<balOct1.size()<<std::endl;


            }else
            {

                assert(par::test::isUniqueAndSorted(balOct1,m_uiCommActive));


                int * sendOctCount= new int[m_uiActiveNpes];
                int * recvOctCount= new int[m_uiActiveNpes];
                int * sendOctOffset = new int [m_uiActiveNpes];
                int * recvOctOffset = new int [m_uiActiveNpes];



                unsigned int sBegin=0,sEnd=0,sResult=0;
                for(unsigned int p=0;p<m_uiActiveNpes;p++)
                    sendOctCount[p]=0;

                unsigned int pCount=0;
                std::vector<unsigned int > searchResultIndex;
                searchResultIndex.resize(m_uiActiveNpes,LOOK_UP_TABLE_DEFAULT);

                if(balOct1.size())
                {
                    while((pCount<m_uiActiveNpes) && (!unBalOctSplitters[pCount].isAncestor(balOct1.front())) &&  (balOct1.front()>unBalOctSplitters[pCount]) ) pCount++;

                    for(unsigned int e=0;((e<balOct1.size()) && (pCount<m_uiActiveNpes));e++)
                    {
                        if(balOct1[e]==unBalOctSplitters[pCount])
                        {
                            searchResultIndex[pCount]=e+1;
                            pCount++;

                        }else if(unBalOctSplitters[pCount].isAncestor(balOct1[e]))
                        {
                            while(((e+1)<balOct1.size()) && (unBalOctSplitters[pCount].isAncestor(balOct1[e+1]))) e++;
                            searchResultIndex[pCount]=e+1;
                            pCount++;

                        }else if(balOct1[e]==unBalOctSplitters[pCount].getParent())
                        {
                            searchResultIndex[pCount]=e+1;
                            pCount++;

                        }else if((e==(balOct1.size()-1)))
                        {
                            searchResultIndex[pCount]=balOct1.size();
                        }

                    }

                    for(unsigned int p=0;p<m_uiActiveNpes;p++)
                    {

                        if(searchResultIndex[p]!=LOOK_UP_TABLE_DEFAULT)
                        {
                            sBegin=sEnd;
                            sEnd=searchResultIndex[p];
                            assert(sBegin<=sEnd);
                            sendOctCount[p]=sEnd-sBegin;
                        }

                    }

                }


                par::Mpi_Alltoall(sendOctCount,recvOctCount,1,m_uiCommActive);

                sendOctOffset[0]=0;
                recvOctOffset[0]=0;

                omp_par::scan(sendOctCount,sendOctOffset,m_uiActiveNpes);
                omp_par::scan(recvOctCount,recvOctOffset,m_uiActiveNpes);

                //std::cout<<" rank: "<<m_uiActiveRank<<" balOct1 size: "<<balOct1.size()<<" sendCount: "<<(sendOctOffset[m_uiActiveNpes-1]+sendOctCount[m_uiActiveNpes-1])<<std::endl;
                assert(balOct1.size()==(sendOctOffset[m_uiActiveNpes-1]+sendOctCount[m_uiActiveNpes-1]));

                std::vector<ot::TreeNode> recvOctBuffer;
                recvOctBuffer.resize(recvOctOffset[m_uiActiveNpes-1]+recvOctCount[m_uiActiveNpes-1]);

                par::Mpi_Alltoallv(&(*(balOct1.begin())),sendOctCount,sendOctOffset,&(*(recvOctBuffer.begin())),recvOctCount,recvOctOffset,m_uiCommActive);
                std::swap(balOct1,recvOctBuffer);
                recvOctBuffer.clear();
                assert(par::test::isUniqueAndSorted(balOct1,m_uiCommActive));

                delete [] sendOctCount;
                delete [] recvOctCount;
                delete [] sendOctOffset;
                delete [] recvOctOffset;

                // 3. synchronize the wavelet remesh flags with 2:1 balanced octree.
                if(balOct1.size())
                {
                    unsigned int count2=0; // element iterator for the mesh2.
                    const ot::TreeNode* pNodes2=&(*(balOct1.begin()));
                    ot::TreeNode *pNodes1=&(*(m_uiAllElements.begin()));

                    assert(m_uiAllElements[m_uiElementLocalBegin]==balOct1.front() || m_uiAllElements[m_uiElementLocalBegin].getParent()==balOct1.front() || m_uiAllElements[m_uiElementLocalBegin]==balOct1.front().getParent());
                    assert(m_uiAllElements[m_uiElementLocalEnd-1]==balOct1.back() || m_uiAllElements[m_uiElementLocalEnd-1].getParent()==balOct1.back() || m_uiAllElements[m_uiElementLocalEnd-1]==balOct1.back().getParent());
                    if(!(m_uiAllElements[m_uiElementLocalBegin]==balOct1.front() || m_uiAllElements[m_uiElementLocalBegin].getParent()==balOct1.front() || m_uiAllElements[m_uiElementLocalBegin]==balOct1.front().getParent()))
                    {
                        std::cout<<"[Remesh Error]: rank: "<<m_uiActiveRank<<" M1 & M2 front alignment failed "<<std::endl;
                        exit(0);
                    }

                    if(!(m_uiAllElements[m_uiElementLocalEnd-1]==balOct1.back() || m_uiAllElements[m_uiElementLocalEnd-1].getParent()==balOct1.back() || m_uiAllElements[m_uiElementLocalEnd-1]==balOct1.back().getParent()))
                    {
                        std::cout<<"[Remesh Error]: rank: "<<m_uiActiveRank<<" M1 & M2 back alignment failed "<<std::endl;
                        exit(0);
                    }

                    /* while( (count2<balOct1.size()) && (pNodes1[m_uiElementLocalBegin]!=pNodes2[count2]) && (pNodes1[m_uiElementLocalBegin].getParent()!=pNodes2[count2]) && (pNodes1[m_uiElementLocalBegin]!=pNodes2[count2].getParent())) count2++;
                     assert(count2<balOct1.size());*/

                    for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                    {
                        if(!(pNodes1[ele]==pNodes2[count2] || pNodes1[ele].getParent()==pNodes2[count2] || pNodes1[ele]==pNodes2[count2].getParent()))
                        {
                            std::cout<<"[Remesh Error]: rank: "<<m_uiActiveRank<<" ele: "<<ele<<"pNodes1: "<<pNodes1[ele]<<" count2: "<<count2<<" pNodes2: "<<pNodes2[count2]<<std::endl;
                            exit(0);
                        }

                        assert(count2<balOct1.size());
                        assert(pNodes1[ele]==pNodes2[count2] || pNodes1[ele].getParent()==pNodes2[count2] || pNodes1[ele]==pNodes2[count2].getParent() );

                        if(pNodes1[ele]==pNodes2[count2].getParent())
                        { // old elements have splitted.
                            assert(pNodes2[count2].getParent()==pNodes2[count2+NUM_CHILDREN-1].getParent());
                            pNodes1[ele].setFlag(((OCT_SPLIT<<NUM_LEVEL_BITS)| pNodes1[ele].getLevel()));
                            count2=count2+NUM_CHILDREN;

                        }else if(pNodes1[ele].getParent()==pNodes2[count2])
                        { // old elements have coarsen
                            assert(pNodes1[ele].getParent()==pNodes1[ele+NUM_CHILDREN-1].getParent());
                            pNodes1[ele].setFlag(((OCT_COARSE<<NUM_LEVEL_BITS)| pNodes1[ele].getLevel()));
                            ele=ele+NUM_CHILDREN-1;
                            count2++;

                        }else
                        {
                            assert(pNodes1[ele]==pNodes2[count2]);
                            pNodes1[ele].setFlag(((OCT_NO_CHANGE<<NUM_LEVEL_BITS) | pNodes1[ele].getLevel()));
                            count2++;
                        }

                    }

                 }

                std::vector<ot::TreeNode> balOct2;
                SFC::parSort::SFC_treeSort(balOct1,balOct2,balOct2,balOct2,ld_tol,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,sfK,m_uiCommActive);
                std::swap(balOct1,balOct2);
                balOct2.clear();

                // repartition balOct1 to ensure that it is not partitioned across, children of the same parent.
                enforceSiblingsAreNotPartitioned(balOct1,m_uiCommActive);
                assert(par::test::isUniqueAndSorted(balOct1,m_uiCommActive));


            }





        }

        ot::Mesh * pMesh = new ot::Mesh(balOct1,1,m_uiElementOrder,m_uiCommGlobal,grainSz,ld_tol,sfK);
        return pMesh;


    }


    void ot::Mesh::getElementalFaceNeighbors(const unsigned int eID,const unsigned int dir, unsigned int* lookup)const
    {
        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        lookup[0]=eID;
        lookup[1]=LOOK_UP_TABLE_DEFAULT;

        unsigned int lk=m_uiE2EMapping[eID*m_uiNumDirections+dir];
        if(lk!=LOOK_UP_TABLE_DEFAULT && m_uiAllElements[lk].getLevel()<=m_uiAllElements[eID].getLevel())
            lookup[1]=lk;

        return;


    }

    void ot::Mesh::getElementalEdgeNeighbors(const unsigned int eID,const unsigned int dir, unsigned int* lookup) const
    {

        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        lookup[0]=eID;
        lookup[1]=LOOK_UP_TABLE_DEFAULT;
        lookup[2]=LOOK_UP_TABLE_DEFAULT;
        lookup[3]=LOOK_UP_TABLE_DEFAULT;
        unsigned level=m_uiAllElements[eID].getLevel();

        unsigned int dir1,dir2;
        unsigned int lk=LOOK_UP_TABLE_DEFAULT;

        if(dir==OCT_DIR_LEFT_DOWN)
        {
            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_DOWN;

        }else if(dir==OCT_DIR_LEFT_UP)
        {
            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_UP;

        }else if(dir==OCT_DIR_LEFT_FRONT)
        {
            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_FRONT;

        }else if(dir==OCT_DIR_LEFT_BACK)
        {
            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_RIGHT_DOWN)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_DOWN;

        }else if(dir==OCT_DIR_RIGHT_UP)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_UP;

        }else if(dir==OCT_DIR_RIGHT_BACK)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_RIGHT_FRONT)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_FRONT;

        }else if(dir==OCT_DIR_UP_BACK)
        {
            dir1=OCT_DIR_UP;
            dir2=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_UP_FRONT)
        {
            dir1=OCT_DIR_UP;
            dir2=OCT_DIR_FRONT;

        }else if(dir==OCT_DIR_DOWN_BACK)
        {
            dir1=OCT_DIR_DOWN;
            dir2=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_DOWN_FRONT)
        {
            dir1=OCT_DIR_DOWN;
            dir2=OCT_DIR_FRONT;

        }


        lookup[1]=m_uiE2EMapping[lookup[0]*m_uiNumDirections+dir1];

        for(unsigned int i=0;i<2;i++)
        {
            if(lookup[i]!=LOOK_UP_TABLE_DEFAULT)
                lookup[i+2]=m_uiE2EMapping[lookup[i]*m_uiNumDirections+dir2];
        }


        for(unsigned int i=1;i<4;i++)
        {
            if(lookup[i]!=LOOK_UP_TABLE_DEFAULT && (m_uiAllElements[lookup[i]].getLevel()>level))
                lookup[i]=LOOK_UP_TABLE_DEFAULT;
        }

        return;


    }

    void ot::Mesh::getElementalVertexNeighbors(const unsigned int eID,const unsigned int dir, unsigned int* lookup)const
    {
        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        lookup[0]=eID;
        for(unsigned int i=1;i<NUM_CHILDREN;i++)
            lookup[i]=LOOK_UP_TABLE_DEFAULT;

        unsigned int level=m_uiAllElements[eID].getLevel();
        unsigned int dir1,dir2,dir3;

        if(dir==OCT_DIR_LEFT_DOWN_BACK)
        {
            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_DOWN;
            dir3=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_RIGHT_DOWN_BACK)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_DOWN;
            dir3=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_LEFT_UP_BACK)
        {

            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_UP;
            dir3=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_RIGHT_UP_BACK)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_UP;
            dir3=OCT_DIR_BACK;

        }else if(dir==OCT_DIR_LEFT_DOWN_FRONT)
        {
            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_DOWN;
            dir3=OCT_DIR_FRONT;

        }else if(dir==OCT_DIR_RIGHT_DOWN_FRONT)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_DOWN;
            dir3=OCT_DIR_FRONT;

        }else if(dir==OCT_DIR_LEFT_UP_FRONT)
        {
            dir1=OCT_DIR_LEFT;
            dir2=OCT_DIR_UP;
            dir3=OCT_DIR_FRONT;

        }else if(dir==OCT_DIR_RIGHT_UP_FRONT)
        {
            dir1=OCT_DIR_RIGHT;
            dir2=OCT_DIR_UP;
            dir3=OCT_DIR_FRONT;

        }

        lookup[1]=m_uiE2EMapping[lookup[0]*m_uiNumDirections+dir1];

        for(unsigned int i=0;i<2;i++)
        {
            if(lookup[i]!=LOOK_UP_TABLE_DEFAULT)
                lookup[i+2]=m_uiE2EMapping[lookup[i]*m_uiNumDirections+dir2];
        }

        for(unsigned int i=0;i<4;i++)
        {
            if(lookup[i]!=LOOK_UP_TABLE_DEFAULT)
                lookup[i+4]=m_uiE2EMapping[lookup[i]*m_uiNumDirections+dir3];
        }


        for(unsigned int i=1;i<NUM_CHILDREN;i++)
        {
            if(lookup[i]!=LOOK_UP_TABLE_DEFAULT && (m_uiAllElements[lookup[i]].getLevel()>level))
                lookup[i]=LOOK_UP_TABLE_DEFAULT;
        }




        return;

    }


    void ot::Mesh::computeElementOwnerRanks(std::vector<unsigned int > & elementOwner)
    {
        // should not be called if the mesh is not active
        if(!m_uiIsActive) return;

        elementOwner.resize(m_uiAllElements.size());
        for(unsigned int e=0;e<m_uiAllElements.size();e++)
            elementOwner[e]=m_uiActiveRank;

        std::vector<ot::SearchKey> ghostElements;
        std::vector<ot::SearchKey>::iterator itSKey;
        for(unsigned int e=m_uiElementPreGhostBegin;e<m_uiElementPreGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int e=m_uiElementPostGhostBegin;e<m_uiElementPostGhostEnd;e++)
        {
            itSKey=ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiAllElements[e]));
            itSKey->addOwner(e);
        }


        for(unsigned int p=0;p<m_uiActiveNpes;p++)
            ghostElements.emplace(ghostElements.end(),ot::SearchKey(m_uiLocalSplitterElements[2*p+1]));

        std::vector<ot::SearchKey> tmpSkeys;
        ot::SearchKey rootSKey(m_uiDim,m_uiMaxDepth);

        SFC::seqSort::SFC_treeSort(&(*(ghostElements.begin())),ghostElements.size(),tmpSkeys,tmpSkeys,tmpSkeys,m_uiMaxDepth,m_uiMaxDepth,rootSKey,ROOT_ROTATION,1,TS_SORT_ONLY);

        tmpSkeys.clear();
        SearchKey tmpSkey;
        unsigned int skip;
        for(unsigned int e=0;e<(ghostElements.size());e++)
        {
            skip=1;
            tmpSkey=ghostElements[e];
            while(((e+skip)<ghostElements.size()) && (ghostElements[e]==ghostElements[e+skip]))
            {
                if(ghostElements[e+skip].getOwner()>=0)tmpSkey.addOwner(ghostElements[e+skip].getOwner());
                skip++;
            }
            tmpSkeys.push_back(tmpSkey);
            assert(skip<=2);
            e+=(skip-1);
        }

        std::swap(ghostElements,tmpSkeys);
        tmpSkeys.clear();


        unsigned int gCount=0;
        for(unsigned int p=0;p<m_uiActiveNpes;p++)
        {

            while(gCount<ghostElements.size() && (ghostElements[gCount]!=m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;

                gCount++;
            }

            if(gCount<ghostElements.size() && (ghostElements[gCount]==m_uiLocalSplitterElements[2*p+1]))
            {
                if(ghostElements[gCount].getOwner()>=0)
                    elementOwner[ghostElements[gCount].getOwner()]=p;
                gCount++;
            }


        }
    }


}
