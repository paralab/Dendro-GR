//
// Created by milinda on 7/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains gradient computation functions for the advection equation
*/
//

#ifndef SFCSORTBENCH_RKTRANSPORTUTILS_H
#define SFCSORTBENCH_RKTRANSPORTUTILS_H

#include "mesh.h"
#include "block.h"
#include <iostream>
#include <vector>
#include "fdCoefficient.h"

namespace adv_param
{
    static const Point domain_min(-M_PI,-M_PI,-M_PI);
    static const Point domain_max(M_PI,M_PI,M_PI);
}


/**
 * @param [in] blk: target block to compute the gradient.
 * @param [in] blkID: index of the block.
 * @param [in] uZipIn: unzipped version of the input vector.
 * @param [out] uZipOut: unzigned version of the output vector.
 *
 * */
template <typename T>
void grad(const ot::Mesh* pMesh,unsigned int dir,const T* uZipIn,T* uZipOut);


template <typename T>
void grad(const ot::Mesh* pMesh,unsigned int dir,const T* uZipIn,T* uZipOut)
{


    const std::vector<ot::Block> blkList=pMesh->getLocalBlockList();
    ot::TreeNode blkNode;
    double h;
    unsigned int blkNpe_1D;

    unsigned int grid_min=0;
    unsigned int grid_max=(1u<<m_uiMaxDepth);
    unsigned int paddWidth;
    const unsigned int stencilSz=5;
    unsigned int lx,ly,lz,offset;


    if(dir==0)
    { // compute the derivative in w.r.t x direction



        for(unsigned int blk=0;blk<blkList.size();blk++)
        {

            blkNode=blkList[blk].getBlockNode();
            paddWidth=blkList[blk].get1DPadWidth();
            blkNpe_1D=blkList[blk].get1DArraySize();
            h=1.0/(blkList[blk].computeDx(adv_param::domain_min,adv_param::domain_max));

            lx=blkList[blk].getAllocationSzX();
            ly=blkList[blk].getAllocationSzY();
            lz=blkList[blk].getAllocationSzZ();

            offset=blkList[blk].getOffset();


            assert(blkNpe_1D>paddWidth);

            if(blkNode.minX()==grid_min)
            {  //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;

                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    {
                        uZipOut[offset+k*(ly*lx)+j*(lx)+paddWidth]=0;
                        for(unsigned int index=0;index<stencilSz;index++)
                            uZipOut[offset+k*(ly*lx)+j*(lx)+paddWidth]+=fd::D1_ORDER_4_FORWARD[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+paddWidth+index]*h;
                    }


                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=(paddWidth+1);i<(2*paddWidth);i++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_UPWIND[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+i+index-1]*h;
                        }



                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=2*paddWidth;i<(blkNpe_1D-2*paddWidth);i++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+i+index-2]*h;
                        }



                if(blkNode.maxX()==grid_max)
                {
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth-1);i++)
                            {
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<stencilSz;index++)
                                    uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_DOWNWIND[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+i+index-3]*h;
                            }



                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+(blkNpe_1D-paddWidth-1)]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+(blkNpe_1D-paddWidth-1)]+=fd::D1_ORDER_4_BACKWARD[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+(blkNpe_1D-paddWidth-1)+index-4]*h;
                        }



                }else
                {
                    assert(blkNode.maxX()<grid_max);
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth);i++)
                            {
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<stencilSz;index++)
                                    uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+i+index-2]*h;
                            }


                }

            }else if(blkNode.maxX()==grid_max)
            {

                assert(blkNode.minX()>grid_min);
                assert((blkNpe_1D-2*paddWidth));
                //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-2*paddWidth);i++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+i+index-2]*h;

                        }


                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth-1);i++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_DOWNWIND[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+i+index-3]*h;
                        }



                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    {
                        uZipOut[offset+k*(ly*lx)+j*(lx)+(blkNpe_1D-paddWidth-1)]=0;
                        for(unsigned int index=0;index<stencilSz;index++)
                            uZipOut[offset+k*(ly*lx)+j*(lx)+(blkNpe_1D-paddWidth-1)]+=fd::D1_ORDER_4_BACKWARD[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+(blkNpe_1D-paddWidth-1)+index-4]*h;
                    }



            }else
            {
                assert(blkNode.minX()>grid_min && blkNode.maxX()<grid_max);

                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+j*(lx)+(i+index-2)]*h;
                        }


            }


        }
    }else if(dir==1)
    { // compute the derivative in w.r.t y direction

        for(unsigned int blk=0;blk<blkList.size();blk++)
        {

            blkNode=blkList[blk].getBlockNode();
            paddWidth=blkList[blk].get1DPadWidth();
            blkNpe_1D=blkList[blk].get1DArraySize();
            h=1.0/(blkList[blk].computeDy(adv_param::domain_min,adv_param::domain_max));

            lx=blkList[blk].getAllocationSzX();
            ly=blkList[blk].getAllocationSzY();
            lz=blkList[blk].getAllocationSzZ();

            offset=blkList[blk].getOffset();


            assert(blkNpe_1D>paddWidth);

            if(blkNode.minY()==grid_min)
            {  //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;

                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                    {
                        uZipOut[offset+k*(ly*lx)+(paddWidth)*(lx)+i]=0;
                        for(unsigned int index=0;index<stencilSz;index++)
                            uZipOut[offset+k*(ly*lx)+(paddWidth)*(lx)+i]+=fd::D1_ORDER_4_FORWARD[index]*uZipIn[offset+k*(ly*lx)+(paddWidth+index)*(lx)+i]*h;
                    }


                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int j=(paddWidth+1);j<(2*paddWidth);j++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_UPWIND[index]*uZipIn[offset+k*(ly*lx)+(j+index-1)*(lx)+i]*h;
                        }



                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int j=2*paddWidth;j<(blkNpe_1D-2*paddWidth);j++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+(j+index-2)*(lx)+i]*h;
                        }



                if(blkNode.maxY()==grid_max)
                {
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth-1);j++)
                            {
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<stencilSz;index++)
                                    uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_DOWNWIND[index]*uZipIn[offset+k*(ly*lx)+(j+index-3)*(lx)+i]*h;
                            }



                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        {
                            uZipOut[offset+k*(ly*lx)+(blkNpe_1D-paddWidth-1)*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+(blkNpe_1D-paddWidth-1)*(lx)+i]+=fd::D1_ORDER_4_BACKWARD[index]*uZipIn[offset+k*(ly*lx)+((blkNpe_1D-paddWidth-1)+index-4)*(lx)+i]*h;
                        }



                }else
                {
                    assert(blkNode.maxY()<grid_max);
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth);j++)
                            {
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<stencilSz;index++)
                                    uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+(j+index-2)*(lx)+i]*h;
                            }


                }

            }else if(blkNode.maxY()==grid_max)
            {

                assert(blkNode.minY()>grid_min);
                assert((blkNpe_1D-2*paddWidth));
                //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-2*paddWidth);j++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+(j+index-2)*(lx)+i]*h;

                        }


                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth-1);j++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_DOWNWIND[index]*uZipIn[offset+k*(ly*lx)+(j+index-3)*(lx)+i]*h;
                        }



                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                    {
                        uZipOut[offset+k*(ly*lx)+(blkNpe_1D-paddWidth-1)*(lx)+i]=0;
                        for(unsigned int index=0;index<stencilSz;index++)
                            uZipOut[offset+k*(ly*lx)+(blkNpe_1D-paddWidth-1)*(lx)+i]+=fd::D1_ORDER_4_BACKWARD[index]*uZipIn[offset+k*(ly*lx)+((blkNpe_1D-paddWidth-1)+index-4)*(lx)+i]*h;
                    }



            }else
            {
                assert(blkNode.minY()>grid_min && blkNode.maxY()<grid_max);

                for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+k*(ly*lx)+(j+index-2)*(lx)+i]*h;
                        }


            }


        }

    }else if(dir==2)
    { // compute the derivative in w.r.t z direction

        for(unsigned int blk=0;blk<blkList.size();blk++)
        {

            blkNode=blkList[blk].getBlockNode();
            paddWidth=blkList[blk].get1DPadWidth();
            blkNpe_1D=blkList[blk].get1DArraySize();
            h=1.0/(blkList[blk].computeDz(adv_param::domain_min,adv_param::domain_max));

            lx=blkList[blk].getAllocationSzX();
            ly=blkList[blk].getAllocationSzY();
            lz=blkList[blk].getAllocationSzZ();

            offset=blkList[blk].getOffset();


            assert(blkNpe_1D>paddWidth);

            if(blkNode.minZ()==grid_min)
            {  //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;

                for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                    {
                        uZipOut[offset+(paddWidth)*(ly*lx)+(j)*(lx)+i]=0;
                        for(unsigned int index=0;index<stencilSz;index++)
                            uZipOut[offset+(paddWidth)*(ly*lx)+(j)*(lx)+i]+=fd::D1_ORDER_4_FORWARD[index]*uZipIn[offset+(paddWidth+index)*(ly*lx)+j*(lx)+i]*h;
                    }


                for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int k=(paddWidth+1);k<(2*paddWidth);k++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_UPWIND[index]*uZipIn[offset+(k+index-1)*(ly*lx)+j*(lx)+i]*h;
                        }



                for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int k=2*paddWidth;k<(blkNpe_1D-2*paddWidth);k++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+(k+index-2)*(ly*lx)+j*(lx)+i]*h;
                        }



                if(blkNode.maxZ()==grid_max)
                {
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth-1);k++)
                            {
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<stencilSz;index++)
                                    uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_DOWNWIND[index]*uZipIn[offset+(k+index-3)*(ly*lx)+j*(lx)+i]*h;
                            }



                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        {
                            uZipOut[offset+(blkNpe_1D-paddWidth-1)*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+(blkNpe_1D-paddWidth-1)*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_BACKWARD[index]*uZipIn[offset+((blkNpe_1D-paddWidth-1)+index-4)*(ly*lx)+j*(lx)+i]*h;
                        }



                }else
                {
                    assert(blkNode.maxY()<grid_max);
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth);k++)
                            {
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<stencilSz;index++)
                                    uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+(k+index-2)*(ly*lx)+j*(lx)+i]*h;
                            }


                }

            }else if(blkNode.maxZ()==grid_max)
            {

                assert(blkNode.minZ()>grid_min);
                assert((blkNpe_1D-2*paddWidth));
                //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-2*paddWidth);k++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+(k+index-2)*(ly*lx)+j*(lx)+i]*h;

                        }


                for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth-1);k++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_DOWNWIND[index]*uZipIn[offset+(k+index-3)*(ly*lx)+j*(lx)+i]*h;
                        }



                for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                    {
                        uZipOut[offset+(blkNpe_1D-paddWidth-1)*(ly*lx)+j*(lx)+i]=0;
                        for(unsigned int index=0;index<stencilSz;index++)
                            uZipOut[offset+(blkNpe_1D-paddWidth-1)*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_BACKWARD[index]*uZipIn[offset+((blkNpe_1D-paddWidth-1)+index-4)*(ly*lx)+j*(lx)+i]*h;
                    }



            }else
            {
                assert(blkNode.minZ()>grid_min && blkNode.maxZ()<grid_max);

                for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                    for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        {
                            uZipOut[offset+k*(ly*lx)+j*(lx)+i]=0;
                            for(unsigned int index=0;index<stencilSz;index++)
                                uZipOut[offset+k*(ly*lx)+j*(lx)+i]+=fd::D1_ORDER_4_CENTERED[index]*uZipIn[offset+(k+index-2)*(ly*lx)+j*(lx)+i]*h;
                        }


            }


        }



    }else
    {
        std::cout<<" unknown stencil direction "<<std::endl;
    }



}



#endif //SFCSORTBENCH_RKTRANSPORTUTILS_H
