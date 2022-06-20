/**
 * @file blkAsync.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Fully asynchronous storage to ot::Block data structure for non-uniform time stepping. 
 * @version 0.1
 * @date 2020-03-17
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */

#pragma once
#include "ts.h"
#include "mesh.h"
#include "enutsOp.h"

namespace ts
{
    enum BLK_ASYNC_VEC_MODE{BLK_CG=0, BLK_DG, BLK_UNZIP};


    template<typename T>
    class BlockAsyncVector
    {

        protected: 

            /**@brief: block ID (local)*/
            unsigned int m_uiBlkID = LOOK_UP_TABLE_DEFAULT;

            /**@brief : data vector (rhs in)*/
            T* m_uiVec = NULL;

            unsigned int m_uiSz[3];

            /**@brief: is block synced and ready to be evolved. */
            bool m_uiIsSynced = false;

            /**@brief: type of the vector. */
            BLK_ASYNC_VEC_MODE m_uiMode;

            /**@brief: number of dof. */
            unsigned int m_uiDof;

            
            
        public:
            
            /**@brief: default constructor.*/
            BlockAsyncVector(){};

            /**@brief: default destructor*/
            ~BlockAsyncVector()
            {
                if(m_uiVec!=NULL);
                    delete [] m_uiVec;

                m_uiVec==NULL;
                
            };
            
            /**
             * @brief allocate memory for a block vector
             * @param blk :block id
             * @param sz  : size of the block x,y,z
             * @param synced : true if the block vector is synced 
             * @param mode : BLK_ASYNC_VEC_MODE 
             * @param dof : number of dof
             */
            void createVec(unsigned int blk, const unsigned int* sz, bool synced=false, BLK_ASYNC_VEC_MODE mode = BLK_ASYNC_VEC_MODE::BLK_UNZIP, unsigned int dof=1)
            {
                m_uiBlkID = blk ;
                m_uiSz[0] = sz[0]; m_uiSz[1] = sz[1]; m_uiSz[2] = sz[2]; 
                const unsigned int NN = (m_uiSz[0]*m_uiSz[1]*m_uiSz[2]);
                m_uiIsSynced = synced;
                m_uiMode = mode;
                m_uiDof = dof;

                if(NN==0)
                {
                    m_uiVec=NULL;
                    return;
                }

                m_uiVec = new T [m_uiDof * NN ];

                return;

            }

            /**
             * @brief : destroy a block vector. 
             */
            void destroyVec()
            {
                m_uiBlkID = LOOK_UP_TABLE_DEFAULT;
                m_uiIsSynced = false;
        
                delete [] m_uiVec;
                m_uiVec=NULL;

            }

            /**
             * @brief copy the block buffer from an unzip vector. 
             * @param pMesh : mesh data strucuture. 
             * @param vUnzip : unzip vector. 
             * @param isWithPadding : if true the copy will happen with padding. otherwise padding region is not coppied. 
             * @param dof : number of dofs
             */
            void copyFromUnzip(const ot::Mesh*pMesh, const T* const vUnzip , bool isWithPadding, unsigned int dof=1)
            {
                if(!(pMesh->isActive()))
                    return;
                
                const unsigned int blk = m_uiBlkID;
                const ot::Block* blkList  = pMesh->getLocalBlockList().data();
                
                if(m_uiMode == BLK_ASYNC_VEC_MODE::BLK_UNZIP)
                {
                    const unsigned int unzipSz = pMesh->getDegOfFreedomUnZip();
                    const unsigned int offset  = blkList[blk].getOffset();
                    const unsigned int nx      = blkList[blk].getAllocationSzX();
                    const unsigned int ny      = blkList[blk].getAllocationSzY();
                    const unsigned int nz      = blkList[blk].getAllocationSzZ();
                    const unsigned int pw      = blkList[blk].get1DPadWidth();

                    assert(m_uiSz[0]==nx && m_uiSz[1] ==ny && m_uiSz[2]==nz);

                    if(isWithPadding)
                    {
                        for(unsigned int v =0; v < dof; v++ )
                        {
                            const T* v1 = vUnzip + v * unzipSz;
                            std::memcpy(m_uiVec + v*nx*ny*nz ,&v1[offset],sizeof(T)*nx*ny*nz);
                        }


                    }else
                    {

                        for(unsigned int v =0; v < dof; v++ )
                        {
                            const T* v1 = vUnzip + v * unzipSz;
                            for(unsigned int k= pw; k < nz-pw; k++)
                            for(unsigned int j= pw; j < ny-pw; j++)
                            for(unsigned int i= pw; i < nx-pw; i++)
                                m_uiVec[(v*nx*ny*nz) + k*ny*nx + j*nx + i] = v1[offset + k * ny*nx +  j * nx + i];
                        }
                        
                    }
                    


                }

                return ;
            }

            /**
             * @brief Copy block vector form the vecDG. 
             * 
             * @param pMesh : pointer to the mesh 
             * @param vecDG : pointer to the DG vector.   
             * @param dof   : number of DOF  
             */
            void copyFromVecDG(const ot::Mesh*pMesh, const T* const vecDG,unsigned int dof=1)
            {
                if(!(pMesh->isActive()))
                    return;
                
                const unsigned int blk = m_uiBlkID;
                const ot::Block* blkList  = pMesh->getLocalBlockList().data();
                const ot::TreeNode* pNodes = pMesh->getAllElements().data();

                
                const unsigned int unzipSz   = pMesh->getDegOfFreedomUnZip();
                const unsigned int offset    = blkList[blk].getOffset();
                const unsigned int nx        = blkList[blk].getAllocationSzX();
                const unsigned int ny        = blkList[blk].getAllocationSzY();
                const unsigned int nz        = blkList[blk].getAllocationSzZ();
                const unsigned int pw        = blkList[blk].get1DPadWidth();
                const unsigned int regLev    = blkList[blk].getRegularGridLev();
                const unsigned int paddWidth = blkList[blk].get1DPadWidth();

                const unsigned int lx   =   blkList[blk].getAllocationSzX();
                const unsigned int ly   =   blkList[blk].getAllocationSzY();
                const unsigned int lz   =   blkList[blk].getAllocationSzZ();
                const ot::TreeNode blkNode = blkList[blk].getBlockNode();

                const unsigned int eOrder = pMesh->getElementOrder();
                const unsigned int nPe  = pMesh->getNumNodesPerElement(); 

                const unsigned int dgSz = pMesh->getDegOfFreedomDG();
                for(unsigned int v=0; v < dof; v++)
                {
                    for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
                    {
                        const unsigned int ei  =  (pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                        const unsigned int ej  =  (pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                        const unsigned int ek  =  (pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                        for(unsigned int k=0; k < (eOrder+1); k++)
                         for(unsigned int j=0; j < (eOrder+1); j++)
                          for(unsigned int i=0; i < (eOrder+1); i++)
                            m_uiVec[ (v*lx*ly*lz)  + (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)] = vecDG[ (v*dgSz + elem*nPe)  +  k*(eOrder+1)*(eOrder+1) +  j*(eOrder+1) + i] ;
                    }
                
                }

                
            }

            /**
             * @brief copy a given other vector to the this vector. 
             * @param other : BlockAsyncVec that need to be coppied. 
             */
            void vecCopy(const BlockAsyncVector<T>& other)
            {
                m_uiBlkID = other.m_uiBlkID;
                m_uiIsSynced = other.isSynced();
                m_uiMode = other.m_uiMode;
                m_uiDof = other.m_uiDof;

                m_uiSz[0] = other.m_uiSz[0];
                m_uiSz[1] = other.m_uiSz[1];
                m_uiSz[2] = other.m_uiSz[2];

                const unsigned int NN = m_uiSz[0]*m_uiSz[1]*m_uiSz[2];
                
                if(m_uiMode == BLK_ASYNC_VEC_MODE::BLK_UNZIP)
                    std::memcpy(m_uiVec,other.m_uiVec,sizeof(T)*NN*m_uiDof);

                return;
            }


            /**
             * @brief perform zip operation for the block
             * @note: !!! only copy the owned elemental indices. 
             * @param pMesh : octree Mesh object.
             * @param zipVec : zipped vector
             * @param dof : number of dof. 
             */
            void zip(ot::Mesh*pMesh, T* zipVec,  unsigned int dof=1) const
            {

                if(!(pMesh->isActive()))
                    return;

                const ot::Block* blkList    = pMesh->getLocalBlockList().data();
                const unsigned int * e2n    = pMesh->getE2NMapping().data();
                const unsigned int * e2n_dg = pMesh->getE2NMapping_DG().data(); 
                const unsigned int lx = m_uiSz[0];
                const unsigned int ly = m_uiSz[1];
                const unsigned int lz = m_uiSz[2];
                
                const unsigned int paddWidth=blkList[m_uiBlkID].get1DPadWidth();

                const unsigned int eOrder = pMesh->getElementOrder();
                const unsigned int nPe = pMesh->getNumNodesPerElement();

                const unsigned int vsz_cg = pMesh->getDegOfFreedom();
                const ot::TreeNode* pNodes = pMesh->getAllElements().data();
                

                
                const ot::TreeNode blkNode =  blkList[m_uiBlkID].getBlockNode();
                const unsigned int regLev  =  blkList[m_uiBlkID].getRegularGridLev();

                unsigned int oEid, oXi, oYj, oZk;
                
                if(m_uiMode == BLK_ASYNC_VEC_MODE::BLK_UNZIP)
                {
                    for(unsigned int v=0; v < dof; v++)
                    {
                        
                        for(unsigned int elem = blkList[m_uiBlkID].getLocalElementBegin(); elem < blkList[m_uiBlkID].getLocalElementEnd(); elem++)
                        {
                            const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                            const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                            const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                            for(unsigned int k=0; k < (eOrder+1); k++)
                            for(unsigned int j=0; j < (eOrder+1); j++)
                            for(unsigned int i=0; i < (eOrder+1); i++)
                            {

                                if((e2n_dg[elem*nPe + k*(eOrder+1)*(eOrder+1)+j*(eOrder+1)+i ] / nPe)==elem)
                                    zipVec[ (v*vsz_cg) + e2n[elem*nPe + k*(eOrder+1)*(eOrder+1) + j* (eOrder+1)+ i]] = m_uiVec[ (v*lx*ly*lz) + (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)];

                                // Note 10/28/2020 : The below zip is wrong you copy the nodal values in the wrong place for non hanging case, 
                                // node index mapping is not that simple.
                                // const unsigned int node_dg = e2n_dg[elem*nPe + k*(eOrder+1)*(eOrder+1)+j*(eOrder+1)+i ];
                                // pMesh->dg2eijk(node_dg,oEid,oXi,oYj,oZk);
                                // const bool isHanging = pMesh->isNodeHanging(elem,i,j,k);
                                // if(!isHanging)
                                //     zipVec[ (v*vsz_cg) + e2n[elem*nPe + k*(eOrder+1)*(eOrder+1) + j* (eOrder+1)+ i]] = m_uiVec[ (v*lx*ly*lz) + (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)];
                                // else
                                // {
                                //     const unsigned int cnum=pNodes[(elem)].getMortonIndex();
                                //     const unsigned int iix = eOrder * (int) (cnum & 1u)  +  i;
                                //     const unsigned int jjy = eOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                //     const unsigned int kkz = eOrder * (int) ((cnum & 4u)>>2u)  +  k;

                                //     if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                //         zipVec[(v*vsz_cg) + e2n[elem*nPe + (kkz>>1u)*(eOrder+1)*(eOrder+1) + (jjy>>1u) * (eOrder+1)+ (iix>>1u)]] = m_uiVec[ (v*lx*ly*lz)  +   (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)];

                                // }
                                
                            }
            
                        }
                    
                    }
                
                }

            }

            /**
             * @brief copy block vector to a DG vector. 
             * @param pMesh 
             * @param dgVec 
             * @param dof 
             */
            void zipDG(ot::Mesh*pMesh, T* dgVec,  unsigned int dof=1) const
            {

                if(!(pMesh->isActive()))
                    return;

                const ot::Block* blkList     =   pMesh->getLocalBlockList().data();
                const unsigned int regLev    =   blkList[m_uiBlkID].getRegularGridLev();
                const ot::TreeNode* pNodes   =   pMesh->getAllElements().data();
                const ot::TreeNode blkNode   =   blkList[m_uiBlkID].getBlockNode();
                const unsigned int eOrder    =   pMesh->getElementOrder();
                const unsigned int paddWidth =   blkList[m_uiBlkID].get1DPadWidth();

                const unsigned int nPe       =  (eOrder+1)*(eOrder+1)*(eOrder+1);

                const unsigned int lx   =   blkList[m_uiBlkID].getAllocationSzX();
                const unsigned int ly   =   blkList[m_uiBlkID].getAllocationSzY();
                const unsigned int lz   =   blkList[m_uiBlkID].getAllocationSzZ();

                const unsigned int dgSz = pMesh->getAllElements().size() * pMesh->getNumNodesPerElement();
        

                for(unsigned int v=0; v < dof; v++)
                {
                    for(unsigned int elem = blkList[m_uiBlkID].getLocalElementBegin(); elem < blkList[m_uiBlkID].getLocalElementEnd(); elem++)
                    {
                        const unsigned int ei  =  (pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                        const unsigned int ej  =  (pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                        const unsigned int ek  =  (pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                        for(unsigned int k=0; k < (eOrder+1); k++)
                         for(unsigned int j=0; j < (eOrder+1); j++)
                          for(unsigned int i=0; i < (eOrder+1); i++)
                            dgVec[ (v*dgSz + elem*nPe)  +  k*(eOrder+1)*(eOrder+1) +  j*(eOrder+1) + i] = m_uiVec[ (v*lx*ly*lz)  + (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)];
                    }
                
                }

                return;

            }

           

            /**@brief: returns the m_uiIsSynced status*/
            inline bool isSynced() const { return m_uiIsSynced; }
            
            /**
             * @brief compute the rhs vector based on the given application ctx. 
             * @param appCtx : const pointer to the application contex. 
             */
            template<typename ACtx>
            void computeVec(const ACtx * appCtx, const BlockAsyncVector<T>& in, T time)
            {
                const ot::Mesh* pMesh = appCtx->get_mesh();
                if(!(pMesh->isActive()))
                    return;
                
                m_uiMode = in.m_uiMode;
                m_uiDof = in.m_uiDof;
                m_uiIsSynced = false;
                ((ACtx*)appCtx)->rhs_blk(in.m_uiVec,m_uiVec,m_uiDof,m_uiBlkID,time);
                return;
            }

            /**@brief: returns the m_uiVec pointer. */
            inline const T* data() const { return m_uiVec ; }

            /**@brief: return the block id*/
            inline unsigned int getBlkID() const {return m_uiBlkID;}
            
            /**@brief: return the DOF*/
            inline unsigned int getDOF() const {return m_uiDof;}

            /**@brief: mark the block vector synced. */
            inline void mark_synced() { m_uiIsSynced =true;}

            /**@brief: mark the block vector unsynced. */
            inline void mark_unsynced() { m_uiIsSynced =false;}

            /**@brief: returns the size of the block vector. */
            inline unsigned int getSz() const { return (m_uiSz[0]*m_uiSz[1]*m_uiSz[2]); }

            /**@brief: dump out the vector. */
            void dump_vec(std::ostream & sout)
            {
                const unsigned int nn = m_uiSz[0]*m_uiSz[1]*m_uiSz[2];

                for(unsigned int v =0; v < m_uiDof; v++)
                {
                    sout<<"var : "<<v<<std::endl;
                    for(unsigned int n=0; n < nn; n++)
                        sout<<"m_uiVec["<<n<<"] : "<<m_uiVec[ v*nn + n ]<<std::endl;
                }
            }

             /**
             * @brief Get the Unzip Elemental values from the block async vector. 
             * @param pMesh : Pointer to the mesh class. 
             * @param ele : element ID.
             * @param eleVec : element vector (unzip). 
             */
            
            void getUnzipElementalValues(const ot::Mesh* pMesh, unsigned int ele ,T* eleVec)
            {

                if(pMesh->isActive())
                {

                    const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
                    const ot::TreeNode* pNodes = pMesh->getAllElements().data();

                    if( !(ele >= blkList[m_uiBlkID].getLocalElementBegin()  && ele < blkList[m_uiBlkID].getLocalElementEnd())  )
                    {
                        std::cout<<"Error : "<<__LINE__<<" block async. get unzip elemental values called on the invalid BVec"<<std::endl;
                        return;
                    }
                        

                    ot::TreeNode blkNode = blkList[m_uiBlkID].getBlockNode();

                    const unsigned int lx        = blkList[m_uiBlkID].getAllocationSzX();
                    const unsigned int ly        = blkList[m_uiBlkID].getAllocationSzY();
                    const unsigned int lz        = blkList[m_uiBlkID].getAllocationSzZ();
                    const unsigned int offset    = blkList[m_uiBlkID].getOffset();
                    const unsigned int paddWidth = blkList[m_uiBlkID].get1DPadWidth();

                    const unsigned int regLev    = blkList[m_uiBlkID].getRegularGridLev();

                    const unsigned int ei=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    const unsigned int ej=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    const unsigned int ek=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);
                    const unsigned int eleIDMax = blkList[m_uiBlkID].getElemSz1D();

                    const unsigned int eOrder = pMesh->getElementOrder();

                    const unsigned int ib = ei*eOrder;
                    const unsigned int ie = ei*eOrder + (eOrder+1) + 2*paddWidth;

                    const unsigned int jb = ej*eOrder;
                    const unsigned int je = ej*eOrder + (eOrder+1) + 2*paddWidth;

                    const unsigned int kb = ek*eOrder;
                    const unsigned int ke = ek*eOrder + (eOrder+1) + 2*paddWidth;

                    const unsigned int en[3] = {(eOrder+1) + 2*paddWidth , (eOrder+1) + 2*paddWidth, (eOrder+1) + 2*paddWidth };


                    for(unsigned int v=0; v < m_uiDof; v++)
                    {
                        T* uzipVec = m_uiVec + v*m_uiSz[0]*m_uiSz[1]*m_uiSz[2];
                        T* out     = eleVec  + v*en[0]*en[1]*en[2];


                        for(unsigned int k=kb; k< ke; k++)
                        for(unsigned int j=jb; j< je; j++)
                        for(unsigned int i=ib; i< ie; i++)
                            out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[k*ly*lx + j*lx + i];

            

                        // copy the unzip element last point to the padding region. 
                        // if(pNodes[ele].minX()==0)
                        // {
                        //     assert(ei==0);

                        //     for(unsigned int k=kb; k< ke; k++)
                        //     for(unsigned int j=jb; j< je; j++)
                        //     for(unsigned int i=ib; i< paddWidth; i++)
                        //         out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[k*ly*lx + j*lx + paddWidth];

                        // }


                        // if(pNodes[ele].minY()==0)
                        // {
                        //     assert(ej==0);

                        //     for(unsigned int k=kb; k< ke; k++)
                        //     for(unsigned int j=jb; j< paddWidth; j++)
                        //     for(unsigned int i=ib; i< ie; i++)
                        //         out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[k*ly*lx + paddWidth*lx + i];

                        // }

                        // if(pNodes[ele].minZ()==0)
                        // {
                        //     assert(ek==0);

                        //     for(unsigned int k=kb; k< paddWidth; k++)
                        //     for(unsigned int j=jb; j< je; j++)
                        //     for(unsigned int i=ib; i< ie; i++)
                        //         out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[paddWidth*ly*lx + j*lx + i];

                        // }


                        // if(pNodes[ele].maxX()==(1u<<m_uiMaxDepth))
                        // {
                        //     assert(ei==(eleIDMax-1));

                        //     for(unsigned int k=kb; k< ke; k++)
                        //     for(unsigned int j=jb; j< je; j++)
                        //     for(unsigned int i=(ie-paddWidth); i< ie; i++)
                        //         out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[k*ly*lx + j*lx + (ie-paddWidth-1)];

                        // }

                        // if(pNodes[ele].maxY()==(1u<<m_uiMaxDepth))
                        // {
                        //     assert(ej==(eleIDMax-1));

                        //     for(unsigned int k=kb; k< ke; k++)
                        //     for(unsigned int j=(je-paddWidth); j< je; j++)
                        //     for(unsigned int i=ib; i< ie; i++)
                        //         out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[k*ly*lx + (je-paddWidth-1)*lx + i];

                        // }

                        // if(pNodes[ele].maxZ()==(1u<<m_uiMaxDepth))
                        // {
                        //     assert(ek==(eleIDMax-1));

                        //     for(unsigned int k=(ke-paddWidth); k< ke; k++)
                        //     for(unsigned int j=jb; j< je; j++)
                        //     for(unsigned int i=ib; i< ie; i++)
                        //         out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[(ke-paddWidth-1)*ly*lx + j*lx + i];

                        // }



                    }

                    


                }


                return;
            
            }
            

            void getElementalValues(const ot::Mesh* pMesh, unsigned int ele, T* eleVec)
            {
                
                
                if(pMesh->isActive())
                {

                    const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
                    const ot::TreeNode* pNodes = pMesh->getAllElements().data();

                    if( !(ele >= blkList[m_uiBlkID].getLocalElementBegin()  && ele < blkList[m_uiBlkID].getLocalElementEnd())  )
                    {
                        std::cout<<"Error : "<<__LINE__<<" block async. get unzip elemental values called on the invalid BVec"<<std::endl;
                        return;
                    }
                        

                    ot::TreeNode blkNode = blkList[m_uiBlkID].getBlockNode();

                    const unsigned int lx        = blkList[m_uiBlkID].getAllocationSzX();
                    const unsigned int ly        = blkList[m_uiBlkID].getAllocationSzY();
                    const unsigned int lz        = blkList[m_uiBlkID].getAllocationSzZ();
                    const unsigned int offset    = blkList[m_uiBlkID].getOffset();
                    const unsigned int paddWidth = blkList[m_uiBlkID].get1DPadWidth();

                    const unsigned int regLev    = blkList[m_uiBlkID].getRegularGridLev();

                    const unsigned int ei=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    const unsigned int ej=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    const unsigned int ek=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);
                    const unsigned int eleIDMax = blkList[m_uiBlkID].getElemSz1D();

                    const unsigned int eOrder = pMesh->getElementOrder();

                    const unsigned int ib = ei*eOrder + paddWidth  ;
                    const unsigned int ie = ei*eOrder + (eOrder+1) ;

                    const unsigned int jb = ej*eOrder + paddWidth  ;
                    const unsigned int je = ej*eOrder + (eOrder+1) ;

                    const unsigned int kb = ek*eOrder + paddWidth  ;
                    const unsigned int ke = ek*eOrder + (eOrder+1) ;

                    const unsigned int en[3] = {(eOrder+1) , (eOrder+1) , (eOrder+1) };


                    for(unsigned int v=0; v < m_uiDof; v++)
                    {
                        T* uzipVec = m_uiVec + v*m_uiSz[0]*m_uiSz[1]*m_uiSz[2];
                        T* out     = eleVec  + v*en[0]*en[1]*en[2];


                        for(unsigned int k=kb; k< ke; k++)
                        for(unsigned int j=jb; j< je; j++)
                        for(unsigned int i=ib; i< ie; i++)
                            out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[k*ly*lx + j*lx + i];

                    }

                }


                return;
            }


    };





    



}// end of namespace ts.


