/**
 * @file enutsUtils.tcc
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief LTS utility functions 
 * @version 0.1
 * @date 2021-03-21
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "enuts.h"

namespace ts
{
    template<typename T>
    void sync_blk_padding(const ot::Mesh* pMesh, const T* const dgWVec,std::vector<ts::BlockTimeStep<T>>& bVec,  unsigned int blk, unsigned int bVecIndex,unsigned int dof)
    {
        if((!(pMesh->isActive())))
            return;
        
        const ot::TreeNode* pNodes   =   pMesh->getAllElements().data();
        const ot::Block* blkList     =   pMesh->getLocalBlockList().data();

        const unsigned int regLevel  =   blkList[blk].getRegularGridLev();
        const ot::TreeNode blkNode   =   blkList[blk].getBlockNode();
        const unsigned int PW        =   blkList[blk].get1DPadWidth();


        const unsigned int eOrder    =   pMesh->getElementOrder();
        const unsigned int nPe       =   pMesh->getNumNodesPerElement();
        
        MPI_Comm comm = pMesh->getMPICommunicator();
        
        const unsigned int lx     =  blkList[blk].getAllocationSzX();
        const unsigned int ly     =  blkList[blk].getAllocationSzY();
        const unsigned int lz     =  blkList[blk].getAllocationSzZ();
        const unsigned int offset =  blkList[blk].getOffset(); 

        const unsigned int dgSz  =  pMesh->getAllElements().size() * pMesh->getNumNodesPerElement();
        const unsigned int cgSz  =  pMesh->getDegOfFreedom();
        const unsigned int unSz  =  pMesh->getDegOfFreedomUnZip();   

        const unsigned int* e2n  =  pMesh->getE2NMapping().data();
        const unsigned int* e2e  =  pMesh->getE2EMapping().data();


        const unsigned int bLev  =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
        
        T* m_uiVec = (T*) bVec[blk]._vec[bVecIndex].data();


        // now need to copy to the block unzip/ block asyncVector 
        const double  hx = (1u<<(m_uiMaxDepth-bLev))/(double)eOrder;
        
        const double xmin = blkNode.minX() - PW*hx; const double xmax = blkNode.maxX() + PW*hx;
        const double ymin = blkNode.minY() - PW*hx; const double ymax = blkNode.maxY() + PW*hx;
        const double zmin = blkNode.minZ() - PW*hx; const double zmax = blkNode.maxZ() + PW*hx;

        std::vector<ot::TreeNode> childOct;
        childOct.reserve(NUM_CHILDREN);

        std::vector<T> p2cI;
        p2cI.resize(nPe);

        const double d_compar_tol=1e-6;
        std::vector<unsigned int> eid;
        computeBlockUnzipDepElements(pMesh, blk, eid);

        for(unsigned int m=0; m < eid.size(); m++)
        {
            const unsigned int ele = eid[m];
            
            // no interpolation needed just copy. 
            if(pNodes[ele].getLevel()==bLev)
            {
                const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                const double invhh = 1.0/hh;
                
                for(unsigned int k=0; k < eOrder+1; k++)
                {
                    double zz  = pNodes[ele].minZ() + k*hh;
                    
                    if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                    if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                    if(zz < zmin || zz > zmax) 
                        continue;
                    const unsigned int kkz = std::round((zz-zmin)*invhh);
                    assert( std::fabs(zz-zmin-kkz*hh) < 1e-6);
                    assert(kkz >= 0 && kkz < lz);

                    for(unsigned int j=0; j < eOrder+1; j++)
                    {   
                        double yy  = pNodes[ele].minY() + j*hh;

                        if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                        if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                        if(yy < ymin || yy > ymax) 
                            continue;
                        const unsigned int jjy = std::round((yy-ymin)*invhh);
                        //std::cout<<"yy: "<<yy<<" (ymin + hh*jjy): "<<(ymin + hh*jjy)<<std::endl;
                        assert( std::fabs(yy-ymin-jjy*hh) < 1e-6);
                        assert(jjy>=0 && jjy<ly);

                        for(unsigned int i=0; i < eOrder+1; i++)
                        {
                            double xx = pNodes[ele].minX() + i*hh;

                            if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                            if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                            
                            if(xx < xmin || xx > xmax) 
                                continue;
                            const unsigned int iix = std::round((xx-xmin)*invhh);
                            assert( std::fabs(xx-xmin-iix*hh) < 1e-6);
                            assert(iix>=0 && iix<lx);

                            //std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<kkz<<" , "<<jjy<<", "<<iix<<")"<<" of : "<<lx<<std::endl;

                            // if(blkNode.isAncestor(pNodes[ele]))
                            // {
                            //     const unsigned int ei=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ej=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ek=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                            //     std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" ei: "<<ei<<" ej: "<<ej<<" ek: "<<ek << " xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;

                            // }

                            for(unsigned int v=0; v < dof; v++)
                            {
                                // double v1 = m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix];
                                // double v2 = dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                // if(rk_s==0 && fabs(v2-v1)>1e-3)
                                // {
                                //     std::cout<<"var : "<<v<<" equal level copy: "<<blk<<" old : "<<v1<<" new : "<<v2<<std::endl;    
                                // }
                                m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                            }
                                

                        }
                    
                    }
                
                }


            }
            else if(pNodes[ele].getLevel() > bLev)
            {
                assert((bLev+1) == pNodes[ele].getLevel());
                const unsigned int cnum = pNodes[ele].getMortonIndex();
                ot::TreeNode tmpParent = pNodes[ele].getParent();
                
                const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                const double invhh = 1.0/(2*hh);

                assert(eOrder>1);
                const unsigned int cb =(eOrder%2==0) ? 0 : 1;

                for(unsigned int k=cb; k < eOrder+1; k+=2)
                {
                    double zz  = (pNodes[ele].minZ() + k*hh);
                    if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                    if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                    if(zz < zmin || zz > zmax) 
                        continue;
                    const unsigned int kkz = std::round((zz-zmin)*invhh);
                    assert(kkz >= 0 && kkz < lz);

                    for(unsigned int j=cb; j < eOrder+1; j+=2)
                    {   
                        double yy  = pNodes[ele].minY() + j*hh;
                        
                        if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                        if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                        if(yy < ymin || yy > ymax) 
                            continue;

                        const unsigned int jjy = std::round((yy-ymin)*invhh);
                        assert(jjy>=0 && jjy<ly);

                        for(unsigned int i=cb; i < eOrder+1; i+=2)
                        {
                            double xx = pNodes[ele].minX() + i*hh;

                            if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                            if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                            
                            if(xx < xmin || xx > xmax) 
                                continue;
                            const unsigned int iix = std::round((xx-xmin)*invhh);
                            assert(iix>=0 && iix<lx);


                            // if(blkNode.isAncestor(pNodes[ele]))
                            // {
                            //     const unsigned int ei=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ej=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ek=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                            //     std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" ei: "<<ei<<" ej: "<<ej<<" ek: "<<ek << " xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;

                            // }

                            //std::cout<<"blk: "<<blk<<" blk copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;
                            for(unsigned int v=0; v < dof; v++)
                            {
                                // double v1 = m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix];
                                // double v2 = dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                // if(rk_s==0 && fabs(v2-v1)>1e-3)
                                // {
                                //     std::cout<<"var : "<<v<<" copy from finer octant: "<<blk<<" old : "<<v1<<" new : "<<v2<<std::endl;    
                                // }
                                m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix] = dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                            }
                                
                            

                        }
                    
                    }
                
                }


              
            }
            else
            {
                assert((bLev) == (pNodes[ele].getLevel()+1));
                childOct.clear();
                pNodes[ele].addChildren(childOct); // note this is the ordering of SFC (depends on Hilbert or Morton. )

                for(unsigned int child = 0; child < NUM_CHILDREN; child++)
                {

                    if( (childOct[child].maxX() < xmin || childOct[child].minX() >=xmax)  || (childOct[child].maxY() < ymin || childOct[child].minY() >=ymax) || (childOct[child].maxZ() < zmin || childOct[child].minZ() >=zmax) )
                        continue;


                    //std::cout<<"blk: "<<blk<<" blkNode: "<<blkNode<<" child: "<<child<<" child node "<<childOct[child]<<" parent : "<<pNodes[ele]<<std::endl;
                    const double hh = (1u<<(m_uiMaxDepth - childOct[child].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/hh;

                    for(unsigned int v=0; v < dof; v++)
                    {
                        const unsigned int cnum = childOct[child].getMortonIndex();
                        pMesh->parent2ChildInterpolation(&dgWVec[v*dgSz + ele*nPe],p2cI.data(),cnum,m_uiDim);

                        for(unsigned int k=0; k < eOrder+1; k++)
                        {
                            double zz  = childOct[child].minZ() + k*hh;
                            
                            if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                            if(fabs(zz-zmax)<d_compar_tol) zz=zmax;
                            
                            if(zz < zmin || zz > zmax) 
                                continue;
                            const unsigned int kkz = std::round((zz-zmin)*invhh);
                            assert(kkz >= 0 && kkz < lz);

                            for(unsigned int j=0; j < eOrder+1; j++)
                            {   
                                double yy  = childOct[child].minY() + j*hh;
                                
                                if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                                if(fabs(yy-ymax)<d_compar_tol) yy=ymax;
                                
                                if(yy < ymin || yy > ymax) 
                                    continue;

                                const unsigned int jjy = std::round((yy-ymin)*invhh);
                                assert(jjy>=0 && jjy<ly);

                                for(unsigned int i=0; i < eOrder+1; i++)
                                {
                                    double xx = childOct[child].minX() + i*hh;
                                    
                                    if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                    if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                    
                                    if(xx < xmin || xx > xmax) 
                                        continue;
                                    const unsigned int iix = std::round((xx-xmin)*invhh);
                                    assert(iix>=0 && iix<lx);
                                    //std::cout<<"blk: "<<blk<<" blk copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<" child: "<<child<<std::endl;
                                    // double v1 = m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix];
                                    // double v2 = p2cI[k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                    // if(rk_s==0 && fabs(v2-v1)>1e-3)
                                    // {
                                    //     std::cout<<"var : "<<v<<" copy from coaser octant(p2c): "<<blk<<" old : "<<v1<<" new : "<<v2<<std::endl;    
                                    // }
                                    m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix] =  p2cI[k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];

                                }
                            
                            }
                        
                        }
                        
                    }

                }

            }



        }
        
        for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
        {
            const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
            const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
            const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

            const unsigned int emin = 0;
            const unsigned int emax = (1u<<(regLevel-blkNode.getLevel()))-1;

            //assert(etVec[elem] == bTime );

            // #pragma unroll
            // for(unsigned int v=0; v < dof; v++)
            //     std::memcpy(dgWVec + v*dgSz + elem * nPe , dgStages[rk_s-1] + v*dgSz + elem * nPe, sizeof(T)*nPe );
            
            for(unsigned int v=0; v < dof; v++)
            for(unsigned int k=0;k<(eOrder+1);k++)
            for(unsigned int j=0;j<(eOrder+1);j++)
            for(unsigned int i=0;i<(eOrder+1);i++)
            {
                // double v1 = m_uiVec[v*lx*ly*lz + (ek*eOrder+k+PW)*(ly*lx)+(ej*eOrder+j+PW)*(lx)+(ei*eOrder+i+PW)];
                // double v2 = dgEVar[ v*dgSz + elem * nPe + k*(eOrder+1)*(eOrder+1) + j*(eOrder+1) + i];

                // if(fabs(v1-v2)>1e-3)
                //     std::cout<<"i: "<<i<< "j: "<<j<<" k: "<<k<<" [internal]: v1(old): "<<v1<<" v2:(new): "<<v2<<std::endl;

                m_uiVec[v*lx*ly*lz + (ek*eOrder+k+PW)*(ly*lx)+(ej*eOrder+j+PW)*(lx)+(ei*eOrder+i+PW)]=dgWVec[ v*dgSz + elem * nPe + k*(eOrder+1)*(eOrder+1) + j*(eOrder+1) + i];
            }

        }


        
    }



}