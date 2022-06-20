//
// Created by milinda on 9/6/16.
//

#ifndef SFCSORTBENCH_OCTUTILS_H
#define SFCSORTBENCH_OCTUTILS_H

#include "key.h"
#include "sfcSearch.h"
#include "sfcSort.h"
#include "dendro.h"

/**
  @brief A collection of simple functions for manipulating octrees.
Examples: Regular Refinements, Linearizing an octree, I/O,
Nearest Common Ancestor, adding positive boundaries, marking hanging nodes
@author Rahul S. Sampath, rahul.sampath@gmail.com
@author Hari Sundar, hsundar@gmail.com
@author Milinda Fernando ,milinda@cs.utah.edu

 @remarks Most of the functions used for the mesh generation. Most of the implementations are based on the previous implementation of dendro version 4.0

*/


#include "TreeNode.h"
#include <vector>
#include <assert.h>
#include "mpi.h"
#include "parUtils.h"
#include <functional>
#include "refel.h"
#include "mathUtils.h"
#include "block.h"
#include "sfcSort.h"
#include "dendro.h"
#include "skey.h"

#define OCT2BLK_DECOMP_BLK_FILL_RATIO 0.5 // gurantees how fraction of the block covered by regular octants.
#define OCT2BLK_DECOMP_LEV_GAP 0


/**
 * @author Rahul Sampath
 * @author Hari Sundar
 *
 * @breif add the positive boundary octants for the given octree input.
 * */
void addBoundaryNodesType1(std::vector<ot::TreeNode> &in,
                           std::vector<ot::TreeNode>& bdy,
                           unsigned int dim, unsigned int maxDepth);




int refineOctree(const std::vector<ot::TreeNode> & inp,
                 std::vector<ot::TreeNode> &out);

int refineAndPartitionOctree(const std::vector<ot::TreeNode> & inp,
                             std::vector<ot::TreeNode> &out, MPI_Comm comm);


int createRegularOctree(std::vector<ot::TreeNode>& out, unsigned int lev,unsigned int dim, unsigned int maxDepth, MPI_Comm comm);



/**
   * @author  Hari Sundar
   * @author  Milinda Fernando
   * @brief   Generates an octree based on a function provided by the user based on the Wavelet method to decide on adaptivity.
   * @param   fx        the function that taxes $x,y,z$ coordinates and returns the the value at that point
   * @param   maxDepth  The maximum depth that the octree should be refined to.
   * @param   tol       user specified tolerance for the wavelet representation of the function.
   * @param   elementOrder order of the element when defining the wavelet representation of the function.
   * @param   comm      The MPI communicator to be use for parallel computation.
   *
   * Generates an octree based on a function provided by the user. The function is expected to return the
   * signed distance to the surface that needs to be meshed. The coordinates are expected to be in [0,1]^3.
   */

int function2Octree(std::function<double(double,double,double)> fx, std::vector<ot::TreeNode> & nodes,unsigned int maxDepth, const double & tol ,unsigned int elementOrder,MPI_Comm comm );



/**
   * @author  Hari Sundar
   * @author  Milinda Fernando
   * @brief   Generates an octree based on a function provided by the user based on the Wavelet method to decide on adaptivity.
   * @param[in] fx:        the function that taxes $x,y,z$ coordinates and returns the the value at that point
   * @param[in] numVars: number of total variables computed from fx function.
   * @param[in] varIndex: variable index location that can be used to ascess double * pointer in fx, to determine the refinement of octree.
   * @param[in] numInterpVars: Number of variables considered during the refinement
   * @param[in] maxDepth:  The maximum depth that the octree should be refined to.
   * @param[in] tol       user specified tolerance for the wavelet representation of the function.
   * @param[in] elementOrder order of the element when defining the wavelet representation of the function.
   * @param[in] comm      The MPI communicator to be use for parallel computation.
   *
   * Generates an octree based on a function provided by the user. The function is expected to return the
   * signed distance to the surface that needs to be meshed. The coordinates are expected to be in [0,1]^3.
   */

int function2Octree(std::function<void(double,double,double,double*)> fx,const unsigned int numVars,const unsigned int* varIndex,const unsigned int numInterpVars, std::vector<ot::TreeNode> & nodes,unsigned int maxDepth, const double & tol ,unsigned int elementOrder,MPI_Comm comm );

/**
 * @author Milinda Fernando
 * @brief ensures that the children of the same parent not partioned arcross different partitions.
 * @param [in] in: input octree.
 * @param [in] comm: MPI communicator.
 * @param [out] in : octree with enforced partition constraint.
 *
 * */
void enforceSiblingsAreNotPartitioned(std::vector<ot::TreeNode> & in,MPI_Comm comm);


/**
 * @author Milinda Fernando
 * @brief Performs octree to sequence of finite regular blocks decomposition.
 *
 * @param [in] pNodes 2:1 balance sorted octree. Note that pNodes should be a sorted octree.
 * @param [in] maxDepth maximum depth of the adaptive octree
 * @param [in] coarsetLev: coarset block level allowed
 * @param [out] d_min min depth of the given octree
 * @param [out] d_max max depth of the given octree
 * @param [out] blockList finite sequence of blocks.
 * */

void octree2BlockDecomposition(std::vector<ot::TreeNode>& pNodes, std::vector<ot::Block>& blockList,unsigned int maxDepth,unsigned int & d_min, unsigned int & d_max, DendroIntL localBegin, DendroIntL localEnd,unsigned int eleOrder,unsigned int coarsetLev=0, unsigned int *tag=NULL, unsigned int tsz=0);



/**
 * @author Milinda Fernando
 * @brief Output the blocks, nodes and ghost layer for each block.
 * @param[in] blkList list of blocks.
 * @param[in] pNodes octree which used to compute the block list.
 * @param[in] fNamePrefix prefix name of vtk file
 * @param[in] comm MPI communicator.
 *
 * */

void blockListToVtk(std::vector<ot::Block>& blkList, const std::vector<ot::TreeNode>& pNodes,char* fNamePrefix, MPI_Comm comm);

/**
 * @author Milinda Fernando
 * @brief Computes blucket splitters for a sorted element array.
 * @param [in] pNodes: const ptr to sorted element list
 * @param [in] lev: level of the parent bucket
 * @param [in] maxDepth: max depth of the octree
 * @param [in] rot_id: rotation id of the bucket
 * @param [in] begin: begining of the bucket.
 * @param [in] end : end of the bucket
 * @param [in] splitters: Splitter counts.
 *
 * */
template <typename  T>
void computeSFCBucketSplitters(const T *pNodes, int lev, unsigned int maxDepth,unsigned char rot_id,DendroIntL &begin, DendroIntL &end, DendroIntL *splitters);

/**
 *
 * @param[in] elem: octant that need the edge search keys
 * @param[out] generate the search keys and push them to sKeys.
 *  *
 * */
template<typename T>
void genEdgeSearchKeys(const T& elem,std::vector<ot::SearchKey>& sKeys);


/**
 * @param [in] sKeys: Generated SearchKeys
 * @param [out] keys; Merged keys.
 *
 * */
void mergeKeys(std::vector<ot::SearchKey>& sKeys,std::vector<ot::Key>& keys);



/**@brief : Generates the edge keys for each block. */
void generateBlkEdgeSKeys(const ot::Block & blk, std::vector<ot::SearchKey>& sKeys);


/**@brief : Generates the vertex keys for each block. */
void generateBlkVertexSKeys(const ot::Block & blk, std::vector<ot::SearchKey>& sKeys);



template<typename  T>
bool linearSearch(const T * pNodes, const T& key,unsigned int n,unsigned int sWidth,unsigned int &result);



/**
 * @author Milinda Fernando
 * @brief Performs repartition of the pNodes array based on the provided splitters.
 * @note assumption 1: Assumes that pNodes is sorted and unique, and pNodes contains the splitters. (because of that this might not be a complete octree.)
 * @note splitters are unique and numSplitters equalt to the comm size.
 *
 * @param [in, out] pNodes: Nodes to be repartiioned.
 * @param [in] numNodes: number of nodes
 * @param [in] splitters : splitter elements to determine the partions
 * @param [in] numSplitters: number of splitters (should be equal to the comm size. )
 * @param [in] comm: MPI communicator.
 *
 * */
template <typename  T>
void partitionBasedOnSplitters(std::vector<T>& pNodes, const T* splitters, unsigned int numSplitters,MPI_Comm comm);



/**
 * @brief shrink or expand the input octree while and enforce sorting afterwards.
 * @param[in] in: input octree
 * @param[in] comm1: Number of partitions needed. (or after this method input octree is going to end up in comm1)
 * @para[in] comm2: current communicator.
 * */
template <typename T>
void shrinkOrExpandOctree(std::vector<T> & in,const double ld_tol,const unsigned int sf_k,bool isActive,MPI_Comm activeComm, MPI_Comm globalComm,unsigned int (*getWeight)(const ot::TreeNode *)=NULL);


/**
 * @brief: Performs blocks histogram and compute unzip times for each block.
 *
 * */
template<typename Blk>
void printBlockStats(const Blk* blkList, unsigned int n,unsigned int maxDepth,MPI_Comm comm);


/**
 *@brief: computes the global ranks selected for the local comm. 
 *@param[in] size_global : global comm. size
 *@param[in] rank_global : global rank
 *@param[in] size_local : local size
 *@param[in] rank_i : ith rank of local comm.  
 *@return returns the global rank of ith rank in local comm.  
 */
unsigned int rankSelectRule(unsigned int size_global,unsigned int rank_global, unsigned int size_local,unsigned int rank_i);


/**
 *@brief returns true if the rank_global is selected for size_local comm. based on the rankSelectRule
 *@param[in] size_global : global comm. size
 *@param[in] rank_global : global rank
 *@param[in] size_local : local size
 * @return returns true if the rank_global is selected for size_local comm. based on the rankSelectRule
  */
inline bool isRankSelected(unsigned int size_global,unsigned int rank_global, unsigned int size_local)
{
    bool isSelected=false;
    for(unsigned int p=0;p<size_local;p++)
        if(rank_global==rankSelectRule(size_global,rank_global,size_local,p))
        {
           isSelected=true;
           break;
        }
    return isSelected;
    
}


template <typename  T>
void computeSFCBucketSplitters(const T *pNodes, int lev, unsigned int maxDepth,unsigned char rot_id,DendroIntL &begin, DendroIntL &end, DendroIntL *splitters)
{
    if ((lev >= maxDepth) || (begin == end)) {
        // Special Case when the considering level exceeds the max depth.

        for (int ii = 0; ii < NUM_CHILDREN; ii++) {
            int index = (rotations[2 * NUM_CHILDREN * rot_id + ii] - '0');
            int nextIndex = 0;
            if (ii == (NUM_CHILDREN-1))
                nextIndex = ii + 1;
            else
                nextIndex = (rotations[2 * NUM_CHILDREN * rot_id + ii + 1] - '0');

            if (ii == 0) {
                splitters[index] = begin;
                splitters[nextIndex] = end;
                continue;
            }
            splitters[nextIndex] = splitters[index];
        }
        //std::cout<<"End return "<<"maxDepth "<<maxDepth<<" Lev: "<<lev<< " Begin "<<begin <<" End "<<end<<std::endl;
        return;

    }

    register unsigned int cnum;
    register unsigned int cnum_prev=0;
    DendroIntL num_elements=0;
    unsigned int rotation=0;
    DendroIntL count[(NUM_CHILDREN+2)]={};
    //unsigned int pMaxDepth=(lev);
    //pMaxDepth--;
    unsigned int mid_bit = maxDepth - lev - 1;
    count[0]=begin;
    for (DendroIntL i=begin; i<end; ++i) {

        /*cnum = (lev < pNodes[i].getLevel())? 1 +(((((pNodes[i].getZ() & (1u << mid_bit)) >> mid_bit) << 2u) |
                                                  (((pNodes[i].getY() & (1u << mid_bit)) >> mid_bit) << 1u) |
                                                  ((pNodes[i].getX() & (1u << mid_bit)) >>
                                                   mid_bit))):0;*/

        cnum = (lev < pNodes[i].getLevel())? 1 +( (((pNodes[i].getZ() >> mid_bit) & 1u) << 2u) | (((pNodes[i].getY() >> mid_bit) & 1u) << 1u) | ((pNodes[i].getX() >>mid_bit) & 1u)):0;
        count[cnum+1]++;


    }

    DendroIntL loc[NUM_CHILDREN+1];
    T unsorted[NUM_CHILDREN+1];
    unsigned int live = 0;

    //if(count[1]>0) std::cout<<"For rank: "<<rank<<" count [1]:  "<<count[1]<<std::endl;

    for (unsigned int ii = 0; ii < NUM_CHILDREN; ii++) {
        int index = (rotations[2 * NUM_CHILDREN * rot_id + ii] - '0');
        int nextIndex = 0;
        if (ii == (NUM_CHILDREN-1))
            nextIndex = ii + 1;
        else
            nextIndex = (rotations[2 * NUM_CHILDREN * rot_id + ii + 1] - '0');

        if (ii == 0) {
            splitters[index] = begin;
            splitters[nextIndex] = splitters[index]+count[1]+ count[(index+2)]; // number of elements which needs to come before the others due to level constraint.

        }else {
            splitters[nextIndex] = splitters[index] + count[(index + 2)];
        }
        //   if(count[1]>0 & !rank) std::cout<<" Spliter B:"<<index <<" "<<splitters[index]<<" Splitters E "<<nextIndex<<" "<<splitters[nextIndex]<<std::endl;

    }
}



template <typename  T>
void partitionBasedOnSplitters(std::vector<T>& pNodes, const T* splitters, unsigned int numSplitters,MPI_Comm comm)
{

    int rank, npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);


    std::vector<ot::Key> splitterKeys;
    splitterKeys.resize(numSplitters);

    assert(npes==numSplitters);

    for(unsigned int p=0;p<npes;p++) {
        splitterKeys[p] = ot::Key(splitters[p]);
        pNodes.push_back(splitters[p]);
    }


    std::vector<T> tmpVec;
    T rootNode(m_uiDim,m_uiMaxDepth);



    SFC::seqSort::SFC_treeSort(&(*(pNodes.begin())),pNodes.size(),tmpVec,tmpVec,tmpVec,m_uiMaxDepth,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES);
    std::swap(pNodes,tmpVec);
    tmpVec.clear();

    assert(seq::test::isUniqueAndSorted(pNodes));


    SFC::seqSearch::SFC_treeSearch(&(*(splitterKeys.begin())),&(*(pNodes.begin())),0,numSplitters,0,pNodes.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);

   /* if(!rank)
    {
        for(unsigned int p=0;p<npes;p++)
        {
            std::cout<<" p: "<<p<<" splitterKey: "<<splitterKeys[p]<<" searchResult: "<<splitterKeys[p].getSearchResult()<<std::endl;
        }
    }*/



    int * sendCount=new int[npes];
    int * recvCount=new int[npes];
    int * sendOffset=new int [npes];
    int * recvOffset=new int [npes];

    unsigned int sResult;
    unsigned int sResultPrev=0;
    std::vector<T> sendBuffer;


    assert((splitterKeys[0].getFlag() & OCT_FOUND));
    sendBuffer.resize(sendBuffer.size()+(splitterKeys[0].getSearchResult()));

    for(unsigned int p=0;p<npes;p++) sendCount[p]=0;

    for(unsigned int ele=0;ele<splitterKeys[0].getSearchResult();ele++)
    {
        sendBuffer[ele]=pNodes[ele];
        sendCount[0]++;
    }

    for(unsigned int p=1;p<npes;p++)
    {
        assert( (splitterKeys[p].getFlag() & OCT_FOUND));
        sResultPrev=splitterKeys[p-1].getSearchResult();
        sResult=splitterKeys[p].getSearchResult();

        if((sResultPrev+1)<pNodes.size() && ((sResultPrev+1)<sResult))
        {
            sendBuffer.resize(sendBuffer.size()+(sResult-sResultPrev-1));
            for(unsigned int ele=(sResultPrev+1);ele<(sResult);ele++)
            {
              sendBuffer[sendCount[p-1]+sendCount[p]]=pNodes[ele];
              sendCount[p]++;
            }
        }
    }


    par::Mpi_Alltoall(sendCount,recvCount,1,comm);

    sendOffset[0]=0;
    recvOffset[0]=0;

    omp_par::scan(sendCount,sendOffset,npes);
    omp_par::scan(recvCount,recvOffset,npes);

    pNodes.clear();
    pNodes.resize(recvCount[npes-1]+recvOffset[npes-1]);

    //if(!rank) std::cout<<"rank: "<<rank<<" pNodes size: "<<pNodes.size()<<std::endl;

    assert(sendBuffer.size()==(sendCount[npes-1]+sendOffset[npes-1]));
    par::Mpi_Alltoallv(&(*(sendBuffer.begin())),sendCount,sendOffset,&(*(pNodes.begin())),recvCount,recvOffset,comm);


    for(unsigned int ele=0;ele<pNodes.size();ele++)
    {
      if(pNodes[ele].getLevel()==m_uiMaxDepth)  std::cout<<"rank: "<<rank<<" ele: "<<ele<<" pNodes: "<<pNodes[ele]<<std::endl;

    }



    assert(par::test::isUniqueAndSorted(pNodes,comm));



    delete [] sendCount;
    delete [] recvCount;
    delete [] sendOffset;
    delete [] recvOffset;





}



template<typename T>
void genEdgeSearchKeys(const T& elem,unsigned int blkId,std::vector<ot::SearchKey>& sKeys)
{
    const unsigned int domain_max = 1u<<(m_uiMaxDepth);

    const unsigned int myX=elem.getX();
    const unsigned int myY=elem.getY();
    const unsigned int myZ=elem.getZ();
    const unsigned int mySz=1u<<(m_uiMaxDepth-elem.getLevel());
    std::vector<ot::SearchKey>::iterator hint;
    if(myX>0 && myY>0)
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY-1),(myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_DOWN);
    }

    if(myX>0 && (myY+mySz)<domain_max) {

        hint = sKeys.emplace(sKeys.end(),ot::SearchKey((myX - 1), (myY + mySz), (myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_UP);
    }

    if(myX>0 && myZ>0) {

        hint = sKeys.emplace(sKeys.end(), ot::SearchKey((myX - 1), (myY), (myZ - 1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_BACK);
    }

    if(myX>0 && (myZ+mySz)<domain_max)
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_FRONT);

    }


    if((myX+mySz) < domain_max && myY>0)
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY-1),(myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_DOWN);
    }

    if((myX+mySz)<domain_max && (myY+mySz)<domain_max)
    {

        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY+mySz),(myZ), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_UP);

    }


    if((myX+mySz)<domain_max && myZ>0) {

        hint = sKeys.emplace(sKeys.end(), ot::SearchKey((myX + mySz), (myY), (myZ - 1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_BACK);
    }

    if((myX+mySz)<domain_max && (myZ+mySz)<domain_max)
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_FRONT);

    }

    if(myY>0 && myZ>0)
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_DOWN_BACK);
    }

    if(myY > 0 && (myZ+mySz)<domain_max)
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX),(myY-1),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_DOWN_FRONT);

    }


    if((myY+mySz)<domain_max && myZ>0)
    {

        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX),(myY+mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_UP_BACK);

    }

    if((myY+mySz)<domain_max && (myZ+mySz)<domain_max) {

        hint = sKeys.emplace(sKeys.end(), ot::SearchKey((myX), (myY + mySz), (myZ + mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(blkId);
        hint->addStencilIndexAndDirection(OCT_DIR_UP_FRONT);
    }



}



template<typename  T>
bool linearSearch(const T * pNodes, const T& key,unsigned int n,unsigned int sWidth,unsigned int &result)
{
    unsigned int sBegin=(int)std::max(0,(int)(n-sWidth));
    unsigned int sEnd=(int)std::min(n,n+sWidth);

    for(unsigned int e=sBegin;e<sEnd;e++)
        if(pNodes[e]==key)
        {
            result=e;
            return true;
        }
    result=LOOK_UP_TABLE_DEFAULT;
    return false;

}



template <typename T>
void shrinkOrExpandOctree(std::vector<T> & in,const double ld_tol,const unsigned int sf_k,bool isActive,MPI_Comm activeComm, MPI_Comm globalComm,unsigned int (*getWeight)(const ot::TreeNode *))
{

    int rank_g,npes_g;
    MPI_Comm_rank(globalComm,&rank_g);
    MPI_Comm_size(globalComm,&npes_g);

    if(!rank_g)
    {
        if(!isActive)
        {
            std::cout<<"[Shrink/Expand Error]: active communicator does not include global rank=0. "<<std::endl;
            exit(0);
        }
    }

    int activeCommSz=0;
    if(!rank_g)
        MPI_Comm_size(activeComm,&activeCommSz);

    par::Mpi_Bcast(&activeCommSz,1,0,globalComm);

    assert(activeCommSz<=npes_g);
    if(activeCommSz>npes_g)
    {
        std::cout<<"[Shrink/Expand Error]: active communicator size is larger than the global comm. "<<std::endl;
        exit(0);
    }

    int * sendCount=new int [npes_g];
    int * recvCount=new int [npes_g];
    int * sendOffset=new int [npes_g];
    int * recvOffset=new int [npes_g];


    for(unsigned int i=0;i<npes_g;i++)
        sendCount[i]=0;

    unsigned int localSz=in.size();

    for(unsigned int i=0;i<activeCommSz;i++)
        sendCount[rankSelectRule(npes_g,rank_g,activeCommSz,i)]=(((i+1)*localSz)/activeCommSz) - ((i*localSz)/activeCommSz);

    par::Mpi_Alltoall(sendCount,recvCount,1,globalComm);

    sendOffset[0]=0;
    recvOffset[0]=0;

    omp_par::scan(sendCount,sendOffset,npes_g);
    omp_par::scan(recvCount,recvOffset,npes_g);

    std::vector<T> recvBuf;
    recvBuf.resize(recvOffset[npes_g-1]+recvCount[npes_g-1]);
    
    par::Mpi_Alltoallv_sparse(&(*(in.begin())),sendCount,sendOffset,&(*(recvBuf.begin())),recvCount,recvOffset,globalComm);
    std::swap(in,recvBuf);
    recvBuf.clear();

    delete [] sendCount;
    delete [] recvCount;
    delete [] sendOffset;
    delete [] recvOffset;

    // now std::vector<in> should be in comm1.
    //std::cout<<rank<<" in: "<<in.size()<<std::endl;
    if(isActive)
    {
        T rootTN(m_uiDim,m_uiMaxDepth);
        // @note: should not need a remove duplicates. but we ge duplicates if we did not do it. Just find out why .
        SFC::parSort::SFC_treeSort(in,recvBuf,recvBuf,recvBuf,ld_tol,m_uiMaxDepth,rootTN,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,sf_k,activeComm);
        std::swap(in,recvBuf);
        recvBuf.clear();
        if(getWeight!=NULL ) 
            par::partitionW(in,getWeight,activeComm);
        assert(par::test::isUniqueAndSorted(in,activeComm));
    }



}

template<typename Blk>
void printBlockStats(const Blk* blkList, unsigned int n,unsigned int maxDepth,MPI_Comm comm)
{
    int rank,npes;

    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    unsigned int blk_counts[maxDepth];
    unsigned int blk_counts_g[maxDepth];

    for(unsigned int i=0;i<maxDepth;i++)
        blk_counts[i]=0;

    unsigned int ele1D,index;
    for(unsigned int blk=0;blk<n;blk++)
    {
       ele1D=blkList[blk].get1DArraySize();
       index=(blkList[blk].get1DArraySize()-2*blkList[blk].get1DPadWidth()-1)/blkList[blk].getElementOrder();
       blk_counts[index]++;
    }

    par::Mpi_Reduce(blk_counts,blk_counts_g,maxDepth,MPI_SUM,0,comm);

    if(!rank)
    {
        for(unsigned int k=0;k<maxDepth;k++)
        {
            std::cout<<"blk_lev["<<k<<"]: ";
            for(unsigned int w=0;w<blk_counts_g[k];w++)
                std::cout<<"*";
            std::cout<<std::endl;
        }
    }

    return ;

}

/**@brief compute octree statistics
 * @param[in] in: input octree
 * @param[in] n: size of in array
 * @param[out] octsByLevLocal : array size of m_uiMaxDepth where octsByLev[i] denotes the number of octants for level i for each rank. 
 * @param[out] octsByLevGlobal: reduction of octsByLevLocal with comm, and operator +
 * @param[out] regOcts: precentage of regular octants in the entrie octree. 
*/
template<typename T>
void computeOctreeStats(const T* in, unsigned int n, unsigned int * octsByLevLocal,unsigned int * octsByLevGlobal, double& regOcts,MPI_Comm comm)
{

    unsigned int octsScanByLev[m_uiMaxDepth];

    for(unsigned int i=0;i<m_uiMaxDepth;i++)
    {
        octsByLevLocal[i]=0;
        octsByLevGlobal[i]=0;
    }

    for(unsigned int i=0;i<n;i++)
        octsByLevLocal[in[i].getLevel()]++;

    par::Mpi_Allreduce(octsByLevLocal,octsByLevGlobal,m_uiMaxDepth,MPI_SUM,comm);

    octsScanByLev[0]=0;
    omp_par::scan(octsByLevGlobal,octsScanByLev,m_uiMaxDepth);

    DendroIntL totalOcts=octsScanByLev[m_uiMaxDepth-1]+octsByLevGlobal[m_uiMaxDepth-1];

    regOcts=totalOcts/(double)(1u<<(3*(m_uiMaxDepth-2)));
    return;

}


#endif //SFCSORTBENCH_OCTUTILS_H
