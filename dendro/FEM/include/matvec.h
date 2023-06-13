//
// Created by milinda on 11/20/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains the functionality to perform SFC- based traversal for performing FEM matvec for
 * mesh free and well structured memory access to the array.
 *
 * This file also contains the matvec based on the mesh just to
 * compare against to the sfc-based mesh free matvec.
 *
*/
//

#ifndef SFCSORTBENCH_FEMMATVEC_H
#define SFCSORTBENCH_FEMMATVEC_H

#define NUM_VERTEX_BUCKETS 27
#define NUM_EDGE_BUCKETS 54
#define NUM_FACE_BUCKETS 36
#define NUM_INTERNAL_BUCKETS 8
#define NUM_TOTAL_BUCKETS 125

#define IDX(i,j,k) i+nx*(j+k*ny)
#define IDX2D(i,j) j*ny+i

#include "binUtils.h"
#include "dendro.h"
#include "refel.h"
#include "assert.h"
#include "TreeNode.h"
#include "mesh.h"
#include "operators.h"

#define PT_EXTERNAL 0 // if the point is external to a given element
#define PT_INTERNAL 1 // if the point is internal but not a DG node.
#define PT_INTERNAL_DG_NODE 2 // if the point is internal and it is a DG node.

#define SFC_MVEC_TABLE_DEFAULT UCHAR_MAX
#define SFC_MVEC_TABLE_Rz 9
#define SFC_MVEC_TABLE_Ry 3


namespace fem
{
    namespace seq
    {
        /**@brief generated sfc matvec bucking table, do not change !!! */
        //m_uiDim  3
        static const unsigned char SFC_MATVEC_BUCKET_TABLE[m_uiDim*m_uiDim*m_uiDim][NUM_CHILDREN]={
                {0,255,255,255,255,255,255,255},
                {0,1,255,255,255,255,255,255},
                {1,255,255,255,255,255,255,255},
                {0,2,255,255,255,255,255,255},
                {0,1,2,3,255,255,255,255},
                {1,3,255,255,255,255,255,255},
                {2,255,255,255,255,255,255,255},
                {2,3,255,255,255,255,255,255},
                {3,255,255,255,255,255,255,255},
                {0,4,255,255,255,255,255,255},
                {0,1,4,5,255,255,255,255},
                {1,5,255,255,255,255,255,255},
                {0,2,4,6,255,255,255,255},
                {0,1,2,3,4,5,6,7},
                {1,3,5,7,255,255,255,255},
                {2,6,255,255,255,255,255,255},
                {2,3,6,7,255,255,255,255},
                {3,7,255,255,255,255,255,255},
                {4,255,255,255,255,255,255,255},
                {4,5,255,255,255,255,255,255},
                {5,255,255,255,255,255,255,255},
                {4,6,255,255,255,255,255,255},
                {4,5,6,7,255,255,255,255},
                {5,7,255,255,255,255,255,255},
                {6,255,255,255,255,255,255,255},
                {6,7,255,255,255,255,255,255},
                {7,255,255,255,255,255,255,255}
        };

    }
}


namespace fem
{

    // Note: if possible remove this sequenction namespace once the FEM mesh free matvec implemented for parallel case as well.
    namespace seq
    {

        /***
         * @brief : computes i,j,k for a given nodal point relative to an given element.
         *
         * */
        template<typename T>
        inline unsigned int  computeIJK(const T& pt,const ot::TreeNode elem,const unsigned int dx,const unsigned int logDx,unsigned int *ijk);


        template<typename T>
        inline bool isInternal(const T& pt,const ot::TreeNode elem,const unsigned int dx,const unsigned int logDx);

        /**@brief computes the FEM matrix vector multiplication in the SFC induced traversal ordering.
         * @param [in] pt_coords : spatial coordinates corresponding to vector u
         * @param [in] in : input vector for MATVEC
         * @param [in] dof: number of points (degrees of freedom).
         * @param [out] out: result in vector after matrix-vector product.
         * */
        template <typename T,typename da>
        inline void matvec(const T *pt_coords,const da*  in,unsigned int pBegin, unsigned int pEnd,const unsigned int lev,const unsigned int maxDepth,const unsigned int order,const ot::TreeNode& pOctant,const RefElement * refEl,da *  out,unsigned int mVecType);




    }





}

template<typename T>
inline unsigned int  fem::seq::computeIJK(const T& pt,const ot::TreeNode elem,const unsigned int dx,const unsigned int logDx,unsigned int *ijk)
{
    /*ijk[0]=UINT_MAX;
    ijk[1]=UINT_MAX;
    ijk[2]=UINT_MAX;*/

    //dendro::timer::sfcmatvec::t_computeIJK.start();
    /*const unsigned int upperX=elem.maxX();
    const unsigned int upperY=elem.maxY();
    const unsigned int upperZ=elem.maxZ();

    const unsigned int lowerX=elem.minX();
    const unsigned int lowerY=elem.minY();
    const unsigned int lowerZ=elem.minZ();

    const unsigned int x=pt.xint();
    const unsigned int y=pt.yint();
    const unsigned int z=pt.zint();*/

    /*int p = (x | ~upper) & ((x ^ upper) | (~upper + x));
    int q = (lower | ~x) & ((lower ^ x) | (~x + lower));
    return 1 & ((p & q) >> 31);*/

    //if((((pt.xint() | ~elem.maxX()) & ((pt.xint() ^ elem.maxX()) | (~elem.maxX() + pt.xint())))>>31u) && (((elem.minX() | ~pt.xint()) & ((elem.minX() ^ pt.xint()) | (~pt.xint() + elem.minX())))>>31u)  && (((pt.yint() | ~elem.maxY()) & ((pt.yint() ^ elem.maxY()) | (~elem.maxY() + pt.yint())))>>31u) && (((elem.minY() | ~pt.yint()) & ((elem.minY() ^ pt.yint()) | (~pt.yint() + elem.minY())))>>31u) && (((pt.zint() | ~elem.maxZ()) & ((pt.zint() ^ elem.maxZ()) | (~elem.maxZ() + pt.zint()))) && ((elem.minZ() | ~pt.zint()) & ((elem.minZ() ^ pt.zint()) | (~pt.zint() + elem.minZ())))>>31u) )
    if(elem.minX()<=pt.xint() && elem.minY()<=pt.yint() && elem.minZ()<=pt.zint() && pt.xint()<=elem.maxX() && pt.yint()<=elem.maxY() && pt.zint() <=elem.maxZ())
    {
       // check if diff %dx ==0
       if(!((pt.xint()-elem.minX())&((1u<<logDx)-1)) && !((pt.yint()-elem.minY())&((1u<<logDx)-1)) && !((pt.zint()-elem.minZ())&((1u<<(logDx))-1)) )
       {
           ijk[0]=(pt.xint()-elem.minX())>>logDx;
           ijk[1]=(pt.yint()-elem.minY())>>logDx;
           ijk[2]=(pt.zint()-elem.minZ())>>logDx;

           //dendro::timer::sfcmatvec::t_computeIJK.stop();

           return PT_INTERNAL_DG_NODE;
       }else {
               //dendro::timer::sfcmatvec::t_computeIJK.stop();
               return PT_INTERNAL;
       }


    }else {
        //dendro::timer::sfcmatvec::t_computeIJK.stop();
        return PT_EXTERNAL;
    }


}


template<typename T>
inline bool fem::seq::isInternal(const T& pt,const ot::TreeNode elem,const unsigned int dx,const unsigned int logDx)
{

    return (elem.minX()<=pt.xint() && elem.minY()<=pt.yint() && elem.minZ()<=pt.zint() && pt.xint()<=elem.maxX() && pt.yint()<=elem.maxY() && pt.zint() <=elem.maxZ()) ?  true : false;
}



template <typename T,typename da>
inline void fem::seq::matvec(const T *pt_coords,const da* in,unsigned int pBegin, unsigned int pEnd,const unsigned int lev,const unsigned int maxDepth,const unsigned int order,const ot::TreeNode& pOctant,const RefElement * refEl,da *  out,unsigned int mVecType){


    const unsigned int nPe=(order+1)*(order+1)*(order+1);
    const unsigned int pMaxDepthBit=maxDepth-lev-1;
    assert(lev==pOctant.getLevel());
    assert(((1u<<(maxDepth-lev)) % order)==0);

    const unsigned int sz=(1u<<(maxDepth-pOctant.getLevel()));
    const unsigned int szb2=sz>>1u;

    const unsigned int dx=(1u<<(maxDepth-lev))/order;
    const unsigned int dxb2=(1u<<(maxDepth-lev-1))/order;
    const unsigned int logDx=binOp::fastLog2(dx);
    const unsigned int logDxb2=binOp::fastLog2(dxb2);

    const unsigned int nx=order+1;
    const unsigned int ny=order+1;
    const unsigned int nz=order+1;

    register unsigned int ijk[3];
    unsigned int cnum;
    const unsigned int pMin[]={pOctant.minX(),pOctant.minY(),pOctant.minZ()};
    const unsigned int pMax[]={pOctant.maxX(),pOctant.maxY(),pOctant.maxZ()};
    const unsigned int pMid[]={(pMax[0]+pMin[0])>>1u,(pMax[1]+pMin[1])>>1u,(pMax[2]+pMin[2])>>1u};

    unsigned int bCounts[NUM_CHILDREN];
    for(unsigned int child=0;child<NUM_CHILDREN;child++)
        bCounts[child]=0;



    for(unsigned int p=pBegin;p<pEnd;p++)
          out[p]=(da)0; // initialize the out vector to zero.



    ot::TreeNode childElem[NUM_CHILDREN];


    register unsigned int pt_status;
    childElem[0]=ot::TreeNode(pOctant.minX(),pOctant.minY(),pOctant.minZ(),pOctant.getLevel()+1,m_uiDim,maxDepth);
    childElem[1]=ot::TreeNode(pOctant.minX()+szb2,pOctant.minY(),pOctant.minZ(),pOctant.getLevel()+1,m_uiDim,maxDepth);
    childElem[2]=ot::TreeNode(pOctant.minX(),pOctant.minY()+szb2,pOctant.minZ(),pOctant.getLevel()+1,m_uiDim,maxDepth);
    childElem[3]=ot::TreeNode(pOctant.minX()+szb2,pOctant.minY()+szb2,pOctant.minZ(),pOctant.getLevel()+1,m_uiDim,maxDepth);

    childElem[4]=ot::TreeNode(pOctant.minX(),pOctant.minY(),pOctant.minZ()+szb2,pOctant.getLevel()+1,m_uiDim,maxDepth);
    childElem[5]=ot::TreeNode(pOctant.minX()+szb2,pOctant.minY(),pOctant.minZ()+szb2,pOctant.getLevel()+1,m_uiDim,maxDepth);
    childElem[6]=ot::TreeNode(pOctant.minX(),pOctant.minY()+szb2,pOctant.minZ()+szb2,pOctant.getLevel()+1,m_uiDim,maxDepth);
    childElem[7]=ot::TreeNode(pOctant.minX()+szb2,pOctant.minY()+szb2,pOctant.minZ()+szb2,pOctant.getLevel()+1,m_uiDim,maxDepth);

    // keep track of the interpolations needed.


#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p1_count.start();
#endif

    std::vector<unsigned int> duplicateCoordIndex;
    //pass 1
    for(unsigned int p=pBegin;p<pEnd;p++)
    {

       ijk[2]=0;
       ijk[1]=0;
       ijk[0]=0;

       if(pt_coords[p].zint()>pMid[2]) ijk[2]=2;
       if(pt_coords[p].yint()>pMid[1]) ijk[1]=2;
       if(pt_coords[p].xint()>pMid[0]) ijk[0]=2;

       if(pt_coords[p].zint()<pMid[2]) ijk[2]=0;
       if(pt_coords[p].yint()<pMid[1]) ijk[1]=0;
       if(pt_coords[p].xint()<pMid[0]) ijk[0]=0;

       if(pt_coords[p].zint()==pMid[2]) ijk[2]=1;
       if(pt_coords[p].yint()==pMid[1]) ijk[1]=1;
       if(pt_coords[p].xint()==pMid[0]) ijk[0]=1;

       if(ijk[0]!=1 && ijk[1]!=1 && ijk[2]!=1)
       {

           cnum=(((ijk[2]>>1u)<<2u) | ((ijk[1]>>1u)<<1u) | ((ijk[0]>>1u)));
           bCounts[cnum]++;

       }else{
           duplicateCoordIndex.push_back(p);
       }

    }

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p1_count.stop();
#endif

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p2_count.start();
#endif

    //pass 2
    for(unsigned int index=0;index<duplicateCoordIndex.size();index++)
    {
        ijk[2]=0;
        ijk[1]=0;
        ijk[0]=0;

        if(pt_coords[duplicateCoordIndex[index]].zint()>pMid[2]) ijk[2]=2;
        if(pt_coords[duplicateCoordIndex[index]].yint()>pMid[1]) ijk[1]=2;
        if(pt_coords[duplicateCoordIndex[index]].xint()>pMid[0]) ijk[0]=2;

        if(pt_coords[duplicateCoordIndex[index]].zint()<pMid[2]) ijk[2]=0;
        if(pt_coords[duplicateCoordIndex[index]].yint()<pMid[1]) ijk[1]=0;
        if(pt_coords[duplicateCoordIndex[index]].xint()<pMid[0]) ijk[0]=0;

        if(pt_coords[duplicateCoordIndex[index]].zint()==pMid[2]) ijk[2]=1;
        if(pt_coords[duplicateCoordIndex[index]].yint()==pMid[1]) ijk[1]=1;
        if(pt_coords[duplicateCoordIndex[index]].xint()==pMid[0]) ijk[0]=1;

        for(unsigned int child=0;child<NUM_CHILDREN;child++)
        {

          cnum=SFC_MATVEC_BUCKET_TABLE[(ijk[2]*SFC_MVEC_TABLE_Rz+ijk[1]*SFC_MVEC_TABLE_Ry+ijk[0])][child];
          if(cnum==UCHAR_MAX)
              break;

          //std::cout<<"pt : "<<pt_coords[p].xint()<<" "<<pt_coords[p].yint()<<" "<<pt_coords[p].zint()<<" ijk: "<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<" child: "<<child <<" cnum "<<cnum<<" pMid: "<<pMid[2]<<std::endl;
          //  std::cout<<"pt : "<<pt_coords[p].xint()<<" "<<pt_coords[p].yint()<<" "<<pt_coords[p].zint()<<" ijk: "<<ijk[0]<<" "<<ijk[1]<<" "<<ijk[2]<<" cnum "<<cnum<<" pMid: "<<pMid[2]<<"table 0 : "<<(int)SFC_MATVEC_BUCKET_TABLE[(ijk[2]*SFC_MVEC_TABLE_Rz+ijk[1]*SFC_MVEC_TABLE_Ry+ijk[0])][0]<<" table 1: "<<(int)SFC_MATVEC_BUCKET_TABLE[(ijk[2]*SFC_MVEC_TABLE_Rz+ijk[1]*SFC_MVEC_TABLE_Ry+ijk[0])][1]<<std::endl;
          bCounts[cnum]++;
        }

    }


#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p2_count.stop();
#endif




#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_malloc.start();
#endif

    da** u_internal=new da*[NUM_CHILDREN]; // in vector internal
    da** v_internal=new da*[NUM_CHILDREN]; // in vector internal
    T** pt_internal=new T*[NUM_CHILDREN];


    for(unsigned int child=0;child<NUM_CHILDREN;child++)
    {
        u_internal[child]=NULL;
        v_internal[child]=NULL;
        pt_internal[child]=NULL;

        if(bCounts[child]!=0)
        {
            u_internal[child]=new da[binOp::getNextHighestPowerOfTwo(bCounts[child])];
            v_internal[child]=new da[binOp::getNextHighestPowerOfTwo(bCounts[child])];
            pt_internal[child]=new T[binOp::getNextHighestPowerOfTwo(bCounts[child])];
        }

    }

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_malloc.stop();
#endif

    DendroIntL counts[NUM_CHILDREN];
    for(unsigned int child=0;child<NUM_CHILDREN;child++)
        counts[child]=0;

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p1_cpy.start();
#endif

    //pass 1
    for(unsigned int p=pBegin;p<pEnd;p++)
    {

        ijk[2]=0;
        ijk[1]=0;
        ijk[0]=0;

        if(pt_coords[p].zint()>pMid[2]) ijk[2]=2;
        if(pt_coords[p].yint()>pMid[1]) ijk[1]=2;
        if(pt_coords[p].xint()>pMid[0]) ijk[0]=2;

        if(pt_coords[p].zint()<pMid[2]) ijk[2]=0;
        if(pt_coords[p].yint()<pMid[1]) ijk[1]=0;
        if(pt_coords[p].xint()<pMid[0]) ijk[0]=0;

        if(pt_coords[p].zint()==pMid[2]) ijk[2]=1;
        if(pt_coords[p].yint()==pMid[1]) ijk[1]=1;
        if(pt_coords[p].xint()==pMid[0]) ijk[0]=1;


        if(ijk[0]!=1 && ijk[1]!=1 && ijk[2]!=1)
        {
            cnum=(((ijk[2]>>1u)<<2u) | ((ijk[1]>>1u)<<1u) | ((ijk[0]>>1u)));
            u_internal[cnum][counts[cnum]]=in[p];
            pt_internal[cnum][counts[cnum]]=pt_coords[p];
            counts[cnum]++;
        }


    }
#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p1_cpy.stop();
#endif

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p2_cpy.start();
#endif

    // pass 2;
    for(unsigned int index=0;index<duplicateCoordIndex.size();index++) {

        ijk[2] = 0;
        ijk[1] = 0;
        ijk[0] = 0;

        if (pt_coords[duplicateCoordIndex[index]].zint() > pMid[2]) ijk[2] = 2;
        if (pt_coords[duplicateCoordIndex[index]].yint() > pMid[1]) ijk[1] = 2;
        if (pt_coords[duplicateCoordIndex[index]].xint() > pMid[0]) ijk[0] = 2;

        if (pt_coords[duplicateCoordIndex[index]].zint() < pMid[2]) ijk[2] = 0;
        if (pt_coords[duplicateCoordIndex[index]].yint() < pMid[1]) ijk[1] = 0;
        if (pt_coords[duplicateCoordIndex[index]].xint() < pMid[0]) ijk[0] = 0;

        if (pt_coords[duplicateCoordIndex[index]].zint() == pMid[2]) ijk[2] = 1;
        if (pt_coords[duplicateCoordIndex[index]].yint() == pMid[1]) ijk[1] = 1;
        if (pt_coords[duplicateCoordIndex[index]].xint() == pMid[0]) ijk[0] = 1;

        for(unsigned int child=0;child<NUM_CHILDREN;child++)
        {
            cnum=SFC_MATVEC_BUCKET_TABLE[(ijk[2]*SFC_MVEC_TABLE_Rz+ijk[1]*SFC_MVEC_TABLE_Ry+ijk[0])][child];
            if(cnum==UCHAR_MAX)
                break;

            u_internal[cnum][counts[cnum]]=in[duplicateCoordIndex[index]];
            pt_internal[cnum][counts[cnum]]=pt_coords[duplicateCoordIndex[index]];
            counts[cnum]++;
        }

    }

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p2_cpy.stop();
#endif


#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_malloc.start();
#endif
    da * parentVal=new da[nPe];

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_malloc.stop();
#endif

    unsigned int pNodeCount=0;
    bool computeParent=false;

    for(unsigned int child=0;child<NUM_CHILDREN;child++)
    {
        if(bCounts[child]<nPe)
        {
            computeParent=true;
            break;
        }


    }


#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_parent_bucket.start();
#endif

    if(computeParent)
    {



        // find parent nodal index and corresponding values
        for(unsigned int p=pBegin;p<pEnd;p++)
        {
            if(fem::seq::computeIJK(pt_coords[p],pOctant,dx,logDx,ijk)==PT_INTERNAL_DG_NODE)
            {

                parentVal[IDX(ijk[0],ijk[1],ijk[2])]=in[p];
                pNodeCount++;

            }

        }
    }

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_parent_bucket.stop();
#endif



    for(unsigned int child=0;child<NUM_CHILDREN;child++)
    {


        counts[child]=0;


        if(bCounts[child]>nPe)
        { // this is a non-leaf node. hence recurse.
            fem::seq::matvec(pt_internal[child],u_internal[child],0,bCounts[child],lev+1,maxDepth,order,childElem[child],refEl,v_internal[child],mVecType);
        }else{
            // only reached by leaf nodes.



            da* interp3D_out=NULL;

            da* interp2D_in=NULL;
            da* interp2D_out=NULL;

            da* interp1D_in=NULL;
            da* interp1D_out=NULL;

            da* u_unzip=new da[nPe];
            da* v_unzip=new da[nPe];

            bool* unzipStatus=new bool[nPe];
            bool interpReq[OCT_DIR_TOTAL];
            bool interpCpy[OCT_DIR_TOTAL];

            for(unsigned int i=0;i<nPe;i++)
                unzipStatus[i]=false;

            for(unsigned int i=0;i<OCT_DIR_TOTAL;i++)
            {
                interpReq[i]=false;
                interpCpy[i]=false;
            }


#ifdef MATVEC_PROFILE
            dendro::timer::sfcmatvec::t_pts_bucket.start();
#endif

            for(unsigned int p=0;p<bCounts[child];p++)
            {
                v_internal[child][p]=(da)0;
                pt_status=fem::seq::computeIJK(pt_internal[child][p],childElem[child],dxb2,logDxb2,ijk);
                if(pt_status==PT_INTERNAL_DG_NODE)
                {
                    u_unzip[IDX(ijk[0],ijk[1],ijk[2])]=u_internal[child][p];
                    unzipStatus[IDX(ijk[0],ijk[1],ijk[2])]=true;
                }


            }

#ifdef MATVEC_PROFILE
            dendro::timer::sfcmatvec::t_pts_bucket.stop();
#endif



            if(bCounts[child]!=nPe) {
                // we need to figure out the missing nodes and perform interpolations.

#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_p2cInterp.start();
#endif

                interp3D_out=new da[nPe];
                interp2D_in=new da[(order+1)*(order+1)];
                interp2D_out=new da[(order+1)*(order+1)];

                interp1D_in=new da[(order+1)];
                interp1D_out=new da[(order+1)];





                const unsigned int intepCheckIndex=1;
                //check for missing nodes
                //OCT_DIR_LEFT_DOWN
                if(!unzipStatus[IDX(0,0,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_LEFT_DOWN]=true;
                }


                //OCT_DIR_LEFT_UP
                if(!unzipStatus[IDX(0,order,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_LEFT_UP]=true;
                }


                //OCT_DIR_LEFT_BACK
                if(!unzipStatus[IDX(0,intepCheckIndex,0)])
                {
                    interpReq[OCT_DIR_LEFT_BACK]=true;
                }


                //OCT_DIR_LEFT_FRONT
                if(!unzipStatus[IDX(0,intepCheckIndex,order)])
                {
                    interpReq[OCT_DIR_LEFT_FRONT]=true;
                }


                //OCT_DIR_RIGHT_DOWN
                if(!unzipStatus[IDX(order,0,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_RIGHT_DOWN]=true;
                }


                //OCT_DIR_RIGHT_UP
                if(!unzipStatus[IDX(order,order,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_RIGHT_UP]=true;
                }


                //OCT_DIR_RIGHT_BACK
                if(!unzipStatus[IDX(order,intepCheckIndex,0)])
                {
                    interpReq[OCT_DIR_RIGHT_BACK]=true;
                }


                //OCT_DIR_RIGHT_FRONT
                if(!unzipStatus[IDX(order,intepCheckIndex,order)])
                {
                    interpReq[OCT_DIR_RIGHT_FRONT]=true;
                }


                //OCT_DIR_DOWN_BACK
                if(!unzipStatus[IDX(intepCheckIndex,0,0)])
                {
                    interpReq[OCT_DIR_DOWN_BACK]=true;
                }


                //OCT_DIR_UP_BACK
                if(!unzipStatus[IDX(intepCheckIndex,order,0)])
                {
                    interpReq[OCT_DIR_UP_BACK]=true;
                }


                //OCT_DIR_DOWN_FRONT
                if(!unzipStatus[IDX(intepCheckIndex,0,order)])
                {
                    interpReq[OCT_DIR_DOWN_FRONT]=true;
                }


                //OCT_DIR_UP_FRONT
                if(!unzipStatus[IDX(intepCheckIndex,order,order)])
                {
                    interpReq[OCT_DIR_UP_FRONT]=true;
                }

                // check for faces.

                //OCT_DIR_LEFT
                if(!unzipStatus[IDX(0,intepCheckIndex,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_LEFT]=true;

                }


                //OCT_DIR_RIGT
                if(!unzipStatus[IDX(order,intepCheckIndex,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_RIGHT]=true;


                }


                //OCT_DIR_DOWN
                if(!unzipStatus[IDX(intepCheckIndex,0,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_DOWN]=true;

                }


                // OCT_DIR_UP
                if(!unzipStatus[IDX(intepCheckIndex,order,intepCheckIndex)])
                {
                    interpReq[OCT_DIR_UP]=true;

                }


                //OCT_DIR_BACK
                if(!unzipStatus[IDX(intepCheckIndex,intepCheckIndex,0)])
                {
                    interpReq[OCT_DIR_BACK]=true;

                }


                //OCT_DIR_FRONT
                if(!unzipStatus[IDX(intepCheckIndex,intepCheckIndex,order)])
                {
                    interpReq[OCT_DIR_FRONT]=true;
                }




                assert(pNodeCount==nPe);
                if(pNodeCount!=nPe) std::cout<<"[Error]"<<__func__<<" pNodeCount: "<<pNodeCount<<std::endl;


                //face copy.
                if(interpReq[OCT_DIR_LEFT])
                {

                    assert(binOp::getBit(child,0)==0);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,1));
                    //std::cout<<"cnumL: "<<cnum<<" binOp::getBit(child,2) "<<binOp::getBit(child,2)<<" binOp::getBit(child,1)"<<binOp::getBit(child,1)<<std::endl;
                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int j=0;j<(order+1);j++)
                            interp2D_in[IDX2D(j,k)]=parentVal[IDX(0,j,k)];

                    refEl->I2D_Parent2Child(interp2D_in,interp2D_out,cnum);


                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int j=0;j<(order+1);j++)
                            u_unzip[IDX(0,j,k)]=interp2D_out[IDX2D(j,k)];


                    interpCpy[OCT_DIR_LEFT_DOWN]=true;
                    interpCpy[OCT_DIR_LEFT_UP]=true;
                    interpCpy[OCT_DIR_LEFT_BACK]=true;
                    interpCpy[OCT_DIR_LEFT_FRONT]=true;

                }


                if(interpReq[OCT_DIR_RIGHT])
                {

                    assert(binOp::getBit(child,0)==1);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,1));

                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int j=0;j<(order+1);j++)
                            interp2D_in[IDX2D(j,k)]=parentVal[IDX(order,j,k)];

                    refEl->I2D_Parent2Child(interp2D_in,interp2D_out,cnum);

                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int j=0;j<(order+1);j++)
                            u_unzip[IDX(order,j,k)]=interp2D_out[IDX2D(j,k)];

                    interpCpy[OCT_DIR_RIGHT_DOWN]=true;
                    interpCpy[OCT_DIR_RIGHT_UP]=true;
                    interpCpy[OCT_DIR_RIGHT_BACK]=true;
                    interpCpy[OCT_DIR_RIGHT_FRONT]=true;

                }


                if(interpReq[OCT_DIR_DOWN])
                {

                    assert(binOp::getBit(child,1)==0);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,0));


                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,k)]=parentVal[IDX(i,0,k)];

                    refEl->I2D_Parent2Child(interp2D_in,interp2D_out,cnum);

                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int i=0;i<(order+1);i++)
                            u_unzip[IDX(i,0,k)]=interp2D_out[IDX2D(i,k)];

                    interpCpy[OCT_DIR_RIGHT_DOWN]=true;
                    interpCpy[OCT_DIR_LEFT_DOWN]=true;
                    interpCpy[OCT_DIR_DOWN_BACK]=true;
                    interpCpy[OCT_DIR_DOWN_FRONT]=true;

                }


                if(interpReq[OCT_DIR_UP])
                {

                    assert(binOp::getBit(child,1)==1);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,0));


                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,k)]=parentVal[IDX(i,order,k)];

                    refEl->I2D_Parent2Child(interp2D_in,interp2D_out,cnum);

                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int i=0;i<(order+1);i++)
                            u_unzip[IDX(i,order,k)]=interp2D_out[IDX2D(i,k)];


                    interpCpy[OCT_DIR_RIGHT_UP]=true;
                    interpCpy[OCT_DIR_LEFT_UP]=true;
                    interpCpy[OCT_DIR_UP_BACK]=true;
                    interpCpy[OCT_DIR_UP_FRONT]=true;

                }


                if(interpReq[OCT_DIR_BACK])
                {

                    assert(binOp::getBit(child,2)==0);
                    cnum=(binOp::getBit(child,1)<<1)|(binOp::getBit(child,0));

                    for(unsigned int j=0;j<(order+1);j++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,j)]=parentVal[IDX(i,j,0)];

                    refEl->I2D_Parent2Child(interp2D_in,interp2D_out,cnum);


                    for(unsigned int j=0;j<(order+1);j++)
                        for(unsigned int i=0;i<(order+1);i++)
                            u_unzip[IDX(i,j,0)]=interp2D_out[IDX2D(i,j)];

                    interpCpy[OCT_DIR_RIGHT_BACK]=true;
                    interpCpy[OCT_DIR_LEFT_BACK]=true;
                    interpCpy[OCT_DIR_UP_BACK]=true;
                    interpCpy[OCT_DIR_DOWN_BACK]=true;

                }


                if(interpReq[OCT_DIR_FRONT])
                {

                    assert(binOp::getBit(child,2)==1);
                    cnum=(binOp::getBit(child,1)<<1)|(binOp::getBit(child,0));

                    for(unsigned int j=0;j<(order+1);j++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,j)]=parentVal[IDX(i,j,order)];

                    refEl->I2D_Parent2Child(interp2D_in,interp2D_out,cnum);


                    for(unsigned int j=0;j<(order+1);j++)
                        for(unsigned int i=0;i<(order+1);i++)
                            u_unzip[IDX(i,j,order)]=interp2D_out[IDX2D(i,j)];


                    interpCpy[OCT_DIR_RIGHT_FRONT]=true;
                    interpCpy[OCT_DIR_LEFT_FRONT]=true;
                    interpCpy[OCT_DIR_UP_FRONT]=true;
                    interpCpy[OCT_DIR_DOWN_FRONT]=true;

                }

                //edges copy.
                if( interpReq[OCT_DIR_LEFT_DOWN] && (!interpCpy[OCT_DIR_LEFT_DOWN]))
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,1)==0);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=parentVal[IDX(0,0,k)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=0;k<(order+1);k++)
                        u_unzip[IDX(0,0,k)]=interp1D_out[k];

                }

                if( interpReq[OCT_DIR_LEFT_UP] && (!interpCpy[OCT_DIR_LEFT_UP]))
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,1)==1);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=parentVal[IDX(0,order,k)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=0;k<(order+1);k++)
                        u_unzip[IDX(0,order,k)]=interp1D_out[k];

                }

                if(interpReq[OCT_DIR_LEFT_BACK] && (!interpCpy[OCT_DIR_LEFT_BACK]))
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=parentVal[IDX(0,j,0)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int j=0;j<(order+1);j++)
                        u_unzip[IDX(0,j,0)]=interp1D_out[j];


                }

                if(interpReq[OCT_DIR_LEFT_FRONT] && (!interpCpy[OCT_DIR_LEFT_FRONT]))
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=parentVal[IDX(0,j,order)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);


                    for(unsigned int j=0;j<(order+1);j++)
                        u_unzip[IDX(0,j,order)]=interp1D_out[j];


                }



                if(interpReq[OCT_DIR_RIGHT_DOWN] && (!interpCpy[OCT_DIR_RIGHT_DOWN]))
                {

                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,1)==0);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=parentVal[IDX(order,0,k)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=0;k<(order+1);k++)
                        u_unzip[IDX(order,0,k)]=interp1D_out[k];

                }

                if(interpReq[OCT_DIR_RIGHT_UP] && (!interpCpy[OCT_DIR_RIGHT_UP]))
                {

                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,1)==1);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=parentVal[IDX(order,order,k)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=0;k<(order+1);k++)
                        u_unzip[IDX(order,order,k)]=interp1D_out[k];

                }

                if(interpReq[OCT_DIR_RIGHT_BACK] && (!interpCpy[OCT_DIR_RIGHT_BACK]))
                {
                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=parentVal[IDX(order,j,0)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);


                    for(unsigned int j=0;j<(order+1);j++)
                        u_unzip[IDX(order,j,0)]=interp1D_out[j];

                }

                if(interpReq[OCT_DIR_RIGHT_FRONT] && (!interpCpy[OCT_DIR_RIGHT_FRONT]))
                {

                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=parentVal[IDX(order,j,order)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int j=0;j<(order+1);j++)
                        u_unzip[IDX(order,j,order)]=interp1D_out[j];

                }


                if(interpReq[OCT_DIR_DOWN_BACK] && (!interpCpy[OCT_DIR_DOWN_BACK]))
                {

                    assert(binOp::getBit(child,1)==0 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=parentVal[IDX(i,0,0)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int i=0;i<(order+1);i++)
                        u_unzip[IDX(i,0,0)]=interp1D_out[i];


                }

                if(interpReq[OCT_DIR_DOWN_FRONT] && (!interpCpy[OCT_DIR_DOWN_FRONT]))
                {

                    assert(binOp::getBit(child,1)==0 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=parentVal[IDX(i,0,order)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int i=0;i<(order+1);i++)
                        u_unzip[IDX(i,0,order)]=interp1D_out[i];

                }

                if(interpReq[OCT_DIR_UP_BACK] && (!interpCpy[OCT_DIR_UP_BACK]))
                {

                    assert(binOp::getBit(child,1)==1 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=parentVal[IDX(i,order,0)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int i=0;i<(order+1);i++)
                        u_unzip[IDX(i,order,0)]=interp1D_out[i];

                }

                if(interpReq[OCT_DIR_UP_FRONT] && (!interpCpy[OCT_DIR_UP_FRONT]))
                {

                    assert(binOp::getBit(child,1)==1 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=parentVal[IDX(i,order,order)];

                    refEl->I1D_Parent2Child(interp1D_in,interp1D_out,cnum);

                    for(unsigned int i=0;i<(order+1);i++)
                        u_unzip[IDX(i,order,order)]=interp1D_out[i];

                }

#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_p2cInterp.stop();
#endif

            }


          if(mVecType==0) // mass matvec
                fem::operators::poisson::elementMassMatvec(u_unzip,childElem[child],refEl,v_unzip);
            else if(mVecType==1) // stifness matvec
                fem::operators::poisson::elementStiffnessMatvec(u_unzip,childElem[child],refEl,v_unzip);



            bool * accumStatus = new bool [nPe];
            bool * c2pInterp=new bool[nPe];

            for(unsigned int node=0;node<nPe;node++)
            {
                accumStatus[node]=false;
                c2pInterp[node]=false;
            }


            if(bCounts[child]!=nPe)
            {

#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_c2pInterp.start();
#endif

                unsigned int cnum;

                for(unsigned int node=0;node<nPe;node++)
                    parentVal[node]=(da)0;

                bool vertexAcc[NUM_CHILDREN];
                for(unsigned int i=0;i<NUM_CHILDREN;i++)
                    vertexAcc[i]=false;

                //face copy.
                if(interpReq[OCT_DIR_LEFT])
                {

                    assert(binOp::getBit(child,0)==0);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,1));
                    //std::cout<<"cnumL: "<<cnum<<" binOp::getBit(child,2) "<<binOp::getBit(child,2)<<" binOp::getBit(child,1)"<<binOp::getBit(child,1)<<std::endl;
                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int j=0;j<(order+1);j++)
                            interp2D_in[IDX2D(j,k)]=v_unzip[IDX(0,j,k)];

                    refEl->I2D_Child2Parent(interp2D_in,interp2D_out,cnum);


                    for(unsigned int k=1;k<(order);k++)
                        for(unsigned int j=1;j<(order);j++)
                        {
                            parentVal[IDX(0,j,k)]=interp2D_out[IDX2D(j,k)];
                            c2pInterp[IDX(0,j,k)]=true;
                        }


                }


                if(interpReq[OCT_DIR_RIGHT])
                {

                    assert(binOp::getBit(child,0)==1);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,1));

                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int j=0;j<(order+1);j++)
                            interp2D_in[IDX2D(j,k)]=v_unzip[IDX(order,j,k)];

                    refEl->I2D_Child2Parent(interp2D_in,interp2D_out,cnum);

                    for(unsigned int k=1;k<(order);k++)
                        for(unsigned int j=1;j<(order);j++)
                        {
                            parentVal[IDX(order,j,k)]=interp2D_out[IDX2D(j,k)];
                            c2pInterp[IDX(order,j,k)]=true;
                        }


                }


                if(interpReq[OCT_DIR_DOWN])
                {

                    assert(binOp::getBit(child,1)==0);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,0));


                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,k)]=v_unzip[IDX(i,0,k)];

                    refEl->I2D_Child2Parent(interp2D_in,interp2D_out,cnum);

                    for(unsigned int k=1;k<(order);k++)
                        for(unsigned int i=1;i<(order);i++)
                        {
                            parentVal[IDX(i,0,k)]=interp2D_out[IDX2D(i,k)];
                            c2pInterp[IDX(i,0,k)]=true;
                        }



                }


                if(interpReq[OCT_DIR_UP])
                {

                    assert(binOp::getBit(child,1)==1);
                    cnum=(binOp::getBit(child,2)<<1)|(binOp::getBit(child,0));


                    for(unsigned int k=0;k<(order+1);k++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,k)]=v_unzip[IDX(i,order,k)];

                    refEl->I2D_Child2Parent(interp2D_in,interp2D_out,cnum);

                    for(unsigned int k=1;k<(order);k++)
                        for(unsigned int i=1;i<(order);i++)
                        {
                            parentVal[IDX(i,order,k)]=interp2D_out[IDX2D(i,k)];
                            c2pInterp[IDX(i,order,k)]=true;
                        }



                }


                if(interpReq[OCT_DIR_BACK])
                {

                    assert(binOp::getBit(child,2)==0);
                    cnum=(binOp::getBit(child,1)<<1)|(binOp::getBit(child,0));

                    for(unsigned int j=0;j<(order+1);j++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,j)]=v_unzip[IDX(i,j,0)];

                    refEl->I2D_Child2Parent(interp2D_in,interp2D_out,cnum);


                    for(unsigned int j=1;j<(order);j++)
                        for(unsigned int i=1;i<(order);i++)
                        {
                            parentVal[IDX(i,j,0)]=interp2D_out[IDX2D(i,j)];
                            c2pInterp[IDX(i,j,0)]=true;
                        }




                }


                if(interpReq[OCT_DIR_FRONT])
                {

                    assert(binOp::getBit(child,2)==1);
                    cnum=(binOp::getBit(child,1)<<1)|(binOp::getBit(child,0));

                    for(unsigned int j=0;j<(order+1);j++)
                        for(unsigned int i=0;i<(order+1);i++)
                            interp2D_in[IDX2D(i,j)]=v_unzip[IDX(i,j,order)];

                    refEl->I2D_Child2Parent(interp2D_in,interp2D_out,cnum);


                    for(unsigned int j=1;j<(order);j++)
                        for(unsigned int i=1;i<(order);i++)
                        {
                            parentVal[IDX(i,j,order)]=interp2D_out[IDX2D(i,j)];
                            c2pInterp[IDX(i,j,order)]=true;
                        }




                }

                //edges copy.
                if( interpReq[OCT_DIR_LEFT_DOWN])
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,1)==0);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=v_unzip[IDX(0,0,k)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=1;k<(order);k++)
                    {
                        parentVal[IDX(0,0,k)]=interp1D_out[k];
                        c2pInterp[IDX(0,0,k)]=true;
                    }


                    if(!vertexAcc[0])
                    {
                        parentVal[IDX(0,0,0)]=interp1D_out[0];
                        c2pInterp[IDX(0,0,0)]=true;
                    }


                    if(!vertexAcc[4])
                    {
                        parentVal[IDX(0,0,order)]=interp1D_out[order];
                        c2pInterp[IDX(0,0,order)]=true;
                    }


                    vertexAcc[0]=true;
                    vertexAcc[4]=true;

                }

                if( interpReq[OCT_DIR_LEFT_UP])
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,1)==1);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=v_unzip[IDX(0,order,k)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=1;k<(order);k++)
                    {
                        parentVal[IDX(0,order,k)]=interp1D_out[k];
                        c2pInterp[IDX(0,order,k)]=true;
                    }


                    if(!vertexAcc[2])
                    {
                        parentVal[IDX(0,order,0)]=interp1D_out[0];
                        c2pInterp[IDX(0,order,0)]=true;
                    }


                    if(!vertexAcc[6])
                    {
                        parentVal[IDX(0,order,order)]=interp1D_out[order];
                        c2pInterp[IDX(0,order,order)]=true;
                    }


                    vertexAcc[2]=true;
                    vertexAcc[6]=true;

                }

                if(interpReq[OCT_DIR_LEFT_BACK])
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=v_unzip[IDX(0,j,0)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int j=1;j<(order);j++)
                    {
                        parentVal[IDX(0,j,0)]=interp1D_out[j];
                        c2pInterp[IDX(0,j,0)]=true;
                    }



                    if(!vertexAcc[0])
                    {
                        parentVal[IDX(0,0,0)]=interp1D_out[0];
                        c2pInterp[IDX(0,0,0)]=true;
                    }


                    if(!vertexAcc[2])
                    {
                        parentVal[IDX(0,order,0)]=interp1D_out[order];
                        c2pInterp[IDX(0,order,0)]=true;
                    }



                    vertexAcc[0]=true;
                    vertexAcc[2]=true;

                }

                if(interpReq[OCT_DIR_LEFT_FRONT])
                {

                    assert(binOp::getBit(child,0)==0 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=v_unzip[IDX(0,j,order)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);


                    for(unsigned int j=1;j<(order);j++)
                    {
                        parentVal[IDX(0,j,order)]=interp1D_out[j];
                        c2pInterp[IDX(0,j,order)]=true;
                    }


                    if(!vertexAcc[4])
                    {
                        parentVal[IDX(0,0,order)]=interp1D_out[0];
                        c2pInterp[IDX(0,0,order)]=true;
                    }


                    if(!vertexAcc[6])
                    {
                        parentVal[IDX(0,order,order)]=interp1D_out[order];
                        c2pInterp[IDX(0,order,order)]=true;
                    }


                    vertexAcc[4]=true;
                    vertexAcc[6]=true;

                }



                if(interpReq[OCT_DIR_RIGHT_DOWN])
                {

                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,1)==0);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=v_unzip[IDX(order,0,k)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=1;k<(order);k++)
                    {
                        parentVal[IDX(order,0,k)]=interp1D_out[k];
                        c2pInterp[IDX(order,0,k)]=true;
                    }



                    if(!vertexAcc[1])
                    {
                        parentVal[IDX(order,0,0)]=interp1D_out[0];
                        c2pInterp[IDX(order,0,0)]=true;
                    }


                    if(!vertexAcc[5])
                    {
                        parentVal[IDX(order,0,order)]=interp1D_out[order];
                        c2pInterp[IDX(order,0,order)]=true;
                    }


                    vertexAcc[1]=true;
                    vertexAcc[5]=true;

                }

                if(interpReq[OCT_DIR_RIGHT_UP])
                {

                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,1)==1);
                    cnum=binOp::getBit(child,2);

                    for(unsigned int k=0;k<(order+1);k++)
                        interp1D_in[k]=v_unzip[IDX(order,order,k)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int k=1;k<order;k++)
                    {
                        parentVal[IDX(order,order,k)]=interp1D_out[k];
                        c2pInterp[IDX(order,order,k)]=true;
                    }



                    if(!vertexAcc[3])
                    {
                        parentVal[IDX(order,order,0)]=interp1D_out[0];
                        c2pInterp[IDX(order,order,0)]=true;
                    }


                    if(!vertexAcc[7])
                    {
                        parentVal[IDX(order,order,order)]=interp1D_out[order];
                        c2pInterp[IDX(order,order,order)]=true;
                    }



                    vertexAcc[3]=true;
                    vertexAcc[7]=true;

                }

                if(interpReq[OCT_DIR_RIGHT_BACK])
                {
                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=v_unzip[IDX(order,j,0)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);


                    for(unsigned int j=1;j<(order);j++)
                    {
                        parentVal[IDX(order,j,0)]=interp1D_out[j];
                        c2pInterp[IDX(order,j,0)]=true;
                    }


                    if(!vertexAcc[1])
                    {
                        parentVal[IDX(order,0,0)]=interp1D_out[0];
                        c2pInterp[IDX(order,0,0)]=true;
                    }


                    if(!vertexAcc[3])
                    {
                        parentVal[IDX(order,order,0)]=interp1D_out[order];
                        c2pInterp[IDX(order,order,0)]=true;
                    }



                    vertexAcc[1]=true;
                    vertexAcc[3]=true;


                }

                if(interpReq[OCT_DIR_RIGHT_FRONT])
                {

                    assert(binOp::getBit(child,0)==1 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,1);

                    for(unsigned int j=0;j<(order+1);j++)
                        interp1D_in[j]=v_unzip[IDX(order,j,order)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int j=1;j<(order);j++)
                    {
                        parentVal[IDX(order,j,order)]=interp1D_out[j];
                        c2pInterp[IDX(order,j,order)]=true;
                    }


                    if(!vertexAcc[5])
                    {
                        parentVal[IDX(order,0,order)]=interp1D_out[0];
                        c2pInterp[IDX(order,0,order)]=true;
                    }


                    if(!vertexAcc[7])
                    {
                        parentVal[IDX(order,order,order)]=interp1D_out[order];
                        c2pInterp[IDX(order,order,order)]=true;
                    }



                    vertexAcc[5]=true;
                    vertexAcc[7]=true;

                }


                if(interpReq[OCT_DIR_DOWN_BACK])
                {

                    assert(binOp::getBit(child,1)==0 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=v_unzip[IDX(i,0,0)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int i=1;i<(order);i++)
                    {
                        parentVal[IDX(i,0,0)]=interp1D_out[i];
                        c2pInterp[IDX(i,0,0)]=true;
                    }


                    if(!vertexAcc[0])
                    {
                        parentVal[IDX(0,0,0)]=interp1D_out[0];
                        c2pInterp[IDX(0,0,0)]=true;
                    }


                    if(!vertexAcc[1])
                    {
                        parentVal[IDX(order,0,0)]=interp1D_out[order];
                        c2pInterp[IDX(order,0,0)]=true;
                    }



                    vertexAcc[0]=true;
                    vertexAcc[1]=true;

                }

                if(interpReq[OCT_DIR_DOWN_FRONT])
                {

                    assert(binOp::getBit(child,1)==0 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=v_unzip[IDX(i,0,order)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int i=1;i<(order);i++)
                    {
                        parentVal[IDX(i,0,order)]=interp1D_out[i];
                        c2pInterp[IDX(i,0,order)]=true;
                    }


                    if(!vertexAcc[4])
                    {
                        parentVal[IDX(0,0,order)]=interp1D_out[0];
                        c2pInterp[IDX(0,0,order)]=true;
                    }


                    if(!vertexAcc[5])
                    {
                        parentVal[IDX(order,0,order)]=interp1D_out[order];
                        c2pInterp[IDX(order,0,order)]=true;
                    }



                    vertexAcc[4]=true;
                    vertexAcc[5]=true;


                }

                if(interpReq[OCT_DIR_UP_BACK])
                {

                    assert(binOp::getBit(child,1)==1 && binOp::getBit(child,2)==0);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=v_unzip[IDX(i,order,0)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);

                    for(unsigned int i=1;i<(order);i++)
                    {
                        parentVal[IDX(i,order,0)]=interp1D_out[i];
                        c2pInterp[IDX(i,order,0)]=true;
                    }



                    if(!vertexAcc[2])
                    {
                        parentVal[IDX(0,order,0)]=interp1D_out[0];
                        c2pInterp[IDX(0,order,0)]=true;
                    }



                    if(!vertexAcc[3])
                    {
                        parentVal[IDX(order,order,0)]=interp1D_out[order];
                        c2pInterp[IDX(order,order,0)]=true;
                    }



                    vertexAcc[2]=true;
                    vertexAcc[3]=true;


                }

                if(interpReq[OCT_DIR_UP_FRONT])
                {

                    assert(binOp::getBit(child,1)==1 && binOp::getBit(child,2)==1);
                    cnum=binOp::getBit(child,0);

                    for(unsigned int i=0;i<(order+1);i++)
                        interp1D_in[i]=v_unzip[IDX(i,order,order)];

                    refEl->I1D_Child2Parent(interp1D_in,interp1D_out,cnum);


                    for(unsigned int i=1;i<(order);i++)
                    {
                        parentVal[IDX(i,order,order)]=interp1D_out[i];
                        c2pInterp[IDX(i,order,order)]=true;
                    }




                    if(!vertexAcc[6])
                    {
                        parentVal[IDX(0,order,order)]=interp1D_out[0];
                        c2pInterp[IDX(0,order,order)]=true;
                    }


                    if(!vertexAcc[7])
                    {
                        parentVal[IDX(order,order,order)]=interp1D_out[order];
                        c2pInterp[IDX(order,order,order)]=true;
                    }



                    vertexAcc[6]=true;
                    vertexAcc[7]=true;


                }



                /*for(unsigned int node=0;node<nPe;node++)
                    out[parentIndex[node]]+=parentVal[node];*/


#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_c2pInterp.stop();
#endif


#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_accum.start();
#endif
                for(unsigned int p=pBegin;p<pEnd;p++) {
                    if (fem::seq::computeIJK(pt_coords[p], pOctant, dx, logDx, ijk) == PT_INTERNAL_DG_NODE)
                    {
                        out[p]+=parentVal[IDX(ijk[0],ijk[1],ijk[2])];
                        if(c2pInterp[IDX(ijk[0],ijk[1],ijk[2])])accumStatus[IDX(ijk[0],ijk[1],ijk[2])]=true;
                    }

                }

#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_accum.stop();
#endif



            }


#ifdef MATVEC_PROFILE
            dendro::timer::sfcmatvec::t_accum.start();
#endif

            for(unsigned int p=0;p<bCounts[child];p++)
            {
                pt_status=fem::seq::computeIJK(pt_internal[child][p],childElem[child],dxb2,logDxb2,ijk);
                if(pt_status==PT_INTERNAL_DG_NODE && !accumStatus[IDX(ijk[0],ijk[1],ijk[2])])
                    v_internal[child][p]=v_unzip[IDX(ijk[0],ijk[1],ijk[2])];

            }


#ifdef MATVEC_PROFILE
            dendro::timer::sfcmatvec::t_accum.stop();
#endif


#ifdef MATVEC_PROFILE
            dendro::timer::sfcmatvec::t_malloc.start();
#endif

            delete [] accumStatus;
            delete [] c2pInterp;

            delete [] unzipStatus;
            delete [] interp3D_out;

            delete [] interp2D_in;
            delete [] interp2D_out;

            delete [] interp1D_in;
            delete [] interp1D_out;

            delete [] u_unzip;
            delete [] v_unzip;

#ifdef MATVEC_PROFILE
            dendro::timer::sfcmatvec::t_malloc.stop();
#endif




        }

/*
#ifdef MATVEC_PROFILE
        dendro::timer::sfcmatvec::t_accum.start();
#endif

        for(unsigned int p=pBegin;p<pEnd;p++)
        {
            pt_status=fem::seq::computeIJK(pt_coords[p],childElem[child],dxb2,logDxb2,ijk);
            if(pt_status!=PT_EXTERNAL)
            {
                out[p]+=v_internal[child][counts[child]];
                counts[child]++;
            }

        }


#ifdef MATVEC_PROFILE
        dendro::timer::sfcmatvec::t_accum.stop();
#endif







        delete [] u_internal[child];
        delete [] v_internal[child];
        delete [] pt_internal[child];
*/



    }


#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p1_accum.start();
#endif
        // pass 1
        for(unsigned int p=pBegin;p<pEnd;p++)
        {

            ijk[2]=0;
            ijk[1]=0;
            ijk[0]=0;

            if(pt_coords[p].zint()>pMid[2]) ijk[2]=2;
            if(pt_coords[p].yint()>pMid[1]) ijk[1]=2;
            if(pt_coords[p].xint()>pMid[0]) ijk[0]=2;

            if(pt_coords[p].zint()<pMid[2]) ijk[2]=0;
            if(pt_coords[p].yint()<pMid[1]) ijk[1]=0;
            if(pt_coords[p].xint()<pMid[0]) ijk[0]=0;

            if(pt_coords[p].zint()==pMid[2]) ijk[2]=1;
            if(pt_coords[p].yint()==pMid[1]) ijk[1]=1;
            if(pt_coords[p].xint()==pMid[0]) ijk[0]=1;


            if(ijk[0]!=1 && ijk[1]!=1 && ijk[2]!=1) {

                cnum = (((ijk[2] >> 1u) << 2u) | ((ijk[1] >> 1u) << 1u) | ((ijk[0] >> 1u)));
                out[p]+=v_internal[cnum][counts[cnum]];
                counts[cnum]++;

            }

        }
#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p1_accum.stop();
#endif

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p2_accum.start();
#endif
    // pass 2;
    for(unsigned int index=0;index<duplicateCoordIndex.size();index++) {

        ijk[2] = 0;
        ijk[1] = 0;
        ijk[0] = 0;

        if (pt_coords[duplicateCoordIndex[index]].zint() > pMid[2]) ijk[2] = 2;
        if (pt_coords[duplicateCoordIndex[index]].yint() > pMid[1]) ijk[1] = 2;
        if (pt_coords[duplicateCoordIndex[index]].xint() > pMid[0]) ijk[0] = 2;

        if (pt_coords[duplicateCoordIndex[index]].zint() < pMid[2]) ijk[2] = 0;
        if (pt_coords[duplicateCoordIndex[index]].yint() < pMid[1]) ijk[1] = 0;
        if (pt_coords[duplicateCoordIndex[index]].xint() < pMid[0]) ijk[0] = 0;

        if (pt_coords[duplicateCoordIndex[index]].zint() == pMid[2]) ijk[2] = 1;
        if (pt_coords[duplicateCoordIndex[index]].yint() == pMid[1]) ijk[1] = 1;
        if (pt_coords[duplicateCoordIndex[index]].xint() == pMid[0]) ijk[0] = 1;

        for(unsigned int child=0;child<NUM_CHILDREN;child++)
        {
            cnum=SFC_MATVEC_BUCKET_TABLE[(ijk[2]*SFC_MVEC_TABLE_Rz+ijk[1]*SFC_MVEC_TABLE_Ry+ijk[0])][child];
            if(cnum==UCHAR_MAX)
                break;

            out[duplicateCoordIndex[index]]+=v_internal[cnum][counts[cnum]];
            counts[cnum]++;


        }

    }

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_pts_p2_accum.stop();
#endif



#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_malloc.start();
#endif


    for(unsigned int child=0;child<NUM_CHILDREN;child++)
    {
        delete [] u_internal[child];
        delete [] v_internal[child];
        delete [] pt_internal[child];
    }


    delete [] parentVal;
    delete [] u_internal;
    delete [] v_internal;
    delete [] pt_internal;

#ifdef MATVEC_PROFILE
    dendro::timer::sfcmatvec::t_malloc.stop();
#endif



}


#endif //SFCSORTBENCH_FEMMATVEC_H
