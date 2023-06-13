//
// Created by milinda on 6/29/16.

/**
 * @author: Milinda Fernando
 * School of Computing, University of Utah
 * @breif Contains SFC based sorting functionality (finalized) for Morton and Hilbert Curve.
 * @detail This is combined with the
 * + remove duplicates
 * + octree construction
 * + octree balancing
 *
 * */


#ifndef SFCSORTBENCH_SFCSORT_H
#define SFCSORTBENCH_SFCSORT_H


#include <iostream>
#include <vector>
#include <assert.h>
#include "hcurvedata.h"
#include "TreeNode.h"
#include "treenode2vtk.h"
#include "dendro.h"
#include "testUtils.h"
#include "ompUtils.h"

#include <mpi.h>
#include <chrono>
#include "dtypes.h"
#include "parUtils.h"
#include <set>
#include <unordered_set>
#include "dollar.hpp"
#include <stdio.h>

#ifdef PROFILE_TREE_SORT
#include <chrono>
// for timer
std::vector<double> stats_previous; // stats for the 1st pass of the treeSort.
std::vector<double> stats; // contains the last run of the SFC::parSort::SFC_TreeSort.

/*
 * For option 1, 2, 4:
 *  inputSz (min mean max)
 *  splitter_fix_all (min mean max)
 *  splitter_time (min mean max)
 *  all2all1_time (min mean max)
 *  all2all2_time (min mean max)
 *  localSort_time (min mean max)
 *  remove_duplicates_seq (min mean max)
 *  remove duplicates_par (min mean max)
 *  auxBalOct (min mean max)
 *  2ndPass (min mean max)
 *  Overall total (min mean mex)
 *
 * */


unsigned int inputSz[3];

double splitter_fix_all=0;
double splitter_time=0;
double all2all1_time=0;
double all2all2_time=0;

double localSort_time=0;
double remove_duplicates_seq=0;
double remove_duplicates_par=0;

double auxBalOCt_time=0;

double construction_time=0;
double balancing_time=0;
double total_rd=0;

auto t1=std::chrono::high_resolution_clock::now();
auto t2=std::chrono::high_resolution_clock::now();
auto t3=std::chrono::high_resolution_clock::now();
auto t4=std::chrono::high_resolution_clock::now();
auto t5_sf_staged=std::chrono::high_resolution_clock::now();

double stat_property[3]={};

#endif


typedef unsigned __int128 uint128_t;

#ifdef DIM_2
    #define NUM_CHILDREN 4
    #define ROTATION_OFFSET 8
    #define m_uiDim 2
#else
#define NUM_CHILDREN 8
#define ROTATION_OFFSET 16
#define m_uiDim 3
#endif

#define MILLISECOND_CONVERSION 1e3
#define ROOT_ROTATION 0


/**
 * @detail
 * OPTIONS For TreeSort.
 *
 * NOTE: Consider these options as bits, in the little endian ordering,
 *
 * (2^0 bit location) TS_REMOVE_DUPLICATES : if selected ensures that the output of the algorithm is sorted and unique.
 * (2^1 bit location) TS_CONSTRUCT_OCTREE  : if selected ensures that the output of the algorithm is sorted and completed octree
 * (2^2 bit location) TS_BALANCE_OCTREE    : if selected ensures that the output of the algorithm is sorted completed and balanced.
 *
 *
 * */

#define TS_SORT_ONLY 0
#define TS_REMOVE_DUPLICATES 1
#define TS_CONSTRUCT_OCTREE 2
#define TS_BALANCE_OCTREE 4 // which ensures that balance octree cannote be called without construct octree function.


template <typename T>
struct BucketInfo
{
    unsigned char rot_id;
    unsigned char lev;
    DendroIntL begin;
    DendroIntL end;

    BucketInfo()
    {
        rot_id=0;
        lev=0;

    }

    BucketInfo(unsigned char p_rot_id,unsigned char p_lev,DendroIntL p_begin, DendroIntL p_end)
    {
        rot_id=p_rot_id;
        lev=p_lev;
        begin=p_begin;
        end=p_end;

    }

};


template<typename T>
struct OctreeComp
{

    inline uint64_t splitBy3(unsigned int a) const{
        uint64_t x=a;
        x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
        x = (x | (x << 8)) & 0x00FF00FF00FF00FF;
        x = (x | (x << 4)) & 0x0F0F0F0F0F0F0F0F;
        x = (x | (x << 2)) & 0x3333333333333333;
        x = (x | (x << 1)) & 0x5555555555555555;
        return x;
    }


     inline bool operator()(const T & first, const T &other) const
    {//$

/*        return (((first.getLevel()!=other.getLevel()) && (first.getLevel()>other.getLevel()))
               || ((first.getLevel()==other.getLevel()) && ((first.getX() != other.getX())) && (first.getX() < other.getX()))
               || ((first.getLevel()==other.getLevel()) &&(first.getX() == other.getX()) && (first.getY() != other.getY()) && (first.getY() < other.getY()))
               || ((first.getLevel()==other.getLevel()) && (first.getX() == other.getX()) && (first.getY() == other.getY())  && (first.getZ() < other.getZ())));*/


         /*return (((first.getLevel()==other.getLevel()) || (first.getLevel()<other.getLevel()))
                 && ((first.getLevel()!=other.getLevel()) || ((first.getX() == other.getX())) || (first.getX() > other.getX()))
                 && ((first.getLevel()!=other.getLevel()) || (first.getX() != other.getX()) || (first.getY() == other.getY()) || (first.getY() > other.getY()))
                 && ((first.getLevel()!=other.getLevel()) || (first.getX() != other.getX()) || (first.getY() != other.getY())  || (first.getZ() > other.getZ())));*/


        uint128_t a1=0;
        uint128_t a2=0;
        a1=(((uint128_t)first.getLevel())<<96u) | (((uint128_t)first.getX())<<64u) | (((uint128_t)first.getY())<<32u) |(((uint128_t)first.getZ()));
        a2=(((uint128_t)other.getLevel())<<96u) | (((uint128_t)other.getX())<<64u) | (((uint128_t)other.getY())<<32u) |(((uint128_t)other.getZ()));
        return a1<a2;

/*
        if(first.getLevel()!=other.getLevel())
        {
          if(state!=(first.getLevel()>other.getLevel())) std::cout<<" lev first"<<first<<" other: "<<other<<std::endl;
          return (first.getLevel()>other.getLevel());
        }else
        {
          if(first.getX()!=other.getX())
          {
            if(state!=(first.getX()<other.getX())) std::cout<<" x first"<<first<<" other: "<<other<<std::endl;
            return first.getX()<other.getX();
          }else
          {
             if(first.getY()!=other.getY())
             {
               if(state!=(first.getY()<other.getY())) std::cout<<" y first"<<first<<" other: "<<other<<std::endl;
               return first.getY()<other.getY();
             }else
             {
               if(state!=(first.getZ()<other.getZ())) std::cout<<" z first"<<first<<" other: "<<other<<std::endl;
               return first.getZ()<other.getZ();
             }
          }
        }*/

    }

};
template<typename T>
struct OctreeHash{

    inline uint128_t operator() (const T &first) const {

/*        uint64_t answer1 = 0;
        answer1 |= splitBy3(node.getX()) | splitBy3(node.getY()) << 1 | splitBy3(node.getZ()) << 2;
        return answer1;*/
        /*size_t h1 = std::hash<unsigned int>{}(node.getX());
        size_t h2 = std::hash<unsigned int>{}(node.getY());
        size_t h3 = std::hash<unsigned int>{}(node.getZ());
        size_t h4 = std::hash<unsigned int>{}(node.getLevel());
        return (h4^(h3^((h2 ^ (h1 << 16))<<16))<<16); // or use boost::hash_combine*/
        uint128_t a1=0;
        a1=(((uint128_t)first.getLevel())<<96u) | (((uint128_t)first.getX())<<64u) | (((uint128_t)first.getY())<<32u) |(((uint128_t)first.getZ()));
        return a1;
        //return (((uint64_t)h4<<48u)|((uint64_t)h3<<32u)|((uint64_t)h2<<16u)|((uint64_t)h1));
    }
};




namespace SFC {

    namespace seqSort {


        /*
         *
         * Contains the sequential utils and treeSort function implementations.
         *
         */


        //========================================================= Function declaration begin.=========================================================================================

        /**
         * @author Milinda Fernando
         * @breif Bottom up construction of the auxilary octants which will be needed in the balancing stage.
         * * Assumes that ,
         * 1) input is sorted and complete.
         * */
        template<typename T>
        inline void SFC_bottomUpBalOctantCreation(std::vector<T> & pNodes);


        /**
         * @author Milinda Fernando
         * @breif Sequential version of the tree sort algorithm.
         * @detail
            * OPTIONS For TreeSort.
            *
             * NOTE: Consider these options as bits, in the little endian ordering,
             *
             * (2^0 bit location) TS_REMOVE_DUPLICATES : if selected ensures that the output of the algorithm is sorted and unique.
             * (2^1 bit location) TS_CONSTRUCT_OCTREE  : if selected ensures that the output of the algorithm is sorted and completed octree
             * (2^2 bit location) TS_BALANCE_OCTREE    : if selected ensures that the output of the algorithm is sorted completed and balanced.
         * */
        template<typename T>
        void SFC_treeSort(T* pNodes , DendroIntL n ,std::vector<T>& pOutSorted,std::vector<T>& pOutConstruct,std::vector<T>& pOutBalanced, unsigned int pMaxDepthBit,unsigned int pMaxDepth, T& parent, unsigned int rot_id,unsigned int k, unsigned int options);


        /**
       * @author Milinda Fernando
       * @breif Sequential Bucketing function which will be needed in adaptive load balancing and parallel tree sort implmentation.
       * */

        template<typename T>
        inline void SFC_bucketing(T *pNodes, int lev, unsigned int maxDepth,unsigned char rot_id,DendroIntL &begin, DendroIntL &end, DendroIntL *splitters);

        /**
         *@brief Finds the local optimal node, in pNodes array.
         * @param[in] pNodes : pNodes in consideration
         * @param[in] n: number of elements in the pNodes array
         * @param[in] pMaxDepthBit initial value should be max depth of the tree. this is automatically decrease as we traverse deper in the tree.
         * @param[in] pMaxDepth maximum depth of the tree
         * @param[in] parent: Element which is an ancesotor of all the pNodes elements.
         * @param[in] rot_id: rotaion id of the parent
         * @param[in] minimum returns minimum if true. Else return the local maximum
         * @param[out] optimal: local optimal of the array.
         * */

        template<typename T>
        void SFC_treeSortLocalOptimal(T* pNodes , DendroIntL n , unsigned int pMaxDepthBit,unsigned int pMaxDepth, T& parent, unsigned int rot_id,bool minimum,T& optimal);


        /**
         *@brief Finds the local optimal node, in pNodes array.
         * @param[in] pNodes : pNodes in consideration
         * @param[in] n: number of elements in the pNodes array
         * @param[in] pMaxDepthBit initial value should be max depth of the tree. this is automatically decrease as we traverse deper in the tree.
         * @param[in] pMaxDepth maximum depth of the tree
         * @param[in] parent: Element which is an ancesotor of all the pNodes elements.
         * @param[in] rot_id: rotaion id of the parent
         * @param[out] min: local min of the array.
         * @param[out] max: local max of the array.
         * */

        template<typename T>
        void SFC_treeSortLocalOptimal(T* pNodes , DendroIntL n , unsigned int pMaxDepthBit,unsigned int pMaxDepth, T& parent, unsigned int rot_id,T&min,T&max);


        //========================================================= Function declaration end.==========================================================================================




        template<typename T>
        inline void SFC_bottomUpBalOctantCreation(std::vector<T> & pNodes)
        {//$

            if(pNodes.empty()){ return; }

            T tmpParent;
            /*std::vector<T> tmpSorted;*/
            T root(m_uiDim,m_uiMaxDepth);
            std::set<T,OctreeComp<T> > auxOct(pNodes.begin(),pNodes.end());
            std::set<T,OctreeComp<T> > parentAux;
            unsigned int neighbourCount=0;

            unsigned int newCount=0;
            unsigned int preSz=0;
            unsigned int curr=0;
            unsigned int d_max=(1u<<m_uiMaxDepth);

#ifdef DIM_2
            neighbourCount=8;
#else
            neighbourCount=26;
#endif

            std::pair<typename std::set<T,OctreeComp<T>>::iterator,bool> hint[neighbourCount];
            std::pair<typename std::set<T,OctreeComp<T>>::iterator,bool> hintParent;

           /* std::unordered_set<T,OctreeHash<T>> auxOct(pNodes.begin(),pNodes.end());
            std::pair<typename std::unordered_set<T,OctreeHash<T>>::iterator,bool> hint;*/


            for(unsigned int w=0;w<pNodes.size();++w){
                hintParent=parentAux.emplace(pNodes[w].getParent());
                if(hintParent.second) {
                    tmpParent=*(hintParent.first);

#ifdef DIM_2

                auxOct.emplace(tmpParent.getLeft());
                auxOct.emplace(tmpParent.getRight());
                auxOct.emplace(tmpParent.getFront());
                auxOct.emplace(tmpParent.getBack());
                auxOct.emplace(tmpParent.getLeftBack());
                auxOct.emplace(tmpParent.getRightBack());
                auxOct.emplace(tmpParent.getLeftFront());
                auxOct.emplace(tmpParent.getRightFront());

#else


                    preSz=auxOct.size();
                    hint[0] = auxOct.emplace(tmpParent.getLeft());
                    hint[1] = auxOct.emplace(tmpParent.getLeftBack());
                    hint[2] = auxOct.emplace(tmpParent.getLeftFront());
                    hint[3] = auxOct.emplace(tmpParent.getRight());
                    hint[4] = auxOct.emplace(tmpParent.getRightBack());
                    hint[5] = auxOct.emplace(tmpParent.getRightFront());
                    hint[6] = auxOct.emplace(tmpParent.getBack());
                    hint[7] = auxOct.emplace(tmpParent.getFront());
                    hint[8] = auxOct.emplace(tmpParent.getBottom());
                    hint[9] = auxOct.emplace(tmpParent.getBottomLeft());
                    hint[10] = auxOct.emplace(tmpParent.getBottomLeftBack());
                    hint[11] = auxOct.emplace(tmpParent.getBottomLeftFront());
                    hint[12] = auxOct.emplace(tmpParent.getBottomRight());
                    hint[13] = auxOct.emplace(tmpParent.getBottomRightBack());
                    hint[14] = auxOct.emplace(tmpParent.getBottomRightFront());
                    hint[15] = auxOct.emplace(tmpParent.getBottomBack());
                    hint[16] = auxOct.emplace(tmpParent.getBottomFront());
                    hint[17] = auxOct.emplace(tmpParent.getTop());
                    hint[18] = auxOct.emplace(tmpParent.getTopLeft());
                    hint[19] = auxOct.emplace(tmpParent.getTopLeftBack());
                    hint[20] = auxOct.emplace(tmpParent.getTopLeftFront());
                    hint[21] = auxOct.emplace(tmpParent.getTopRight());
                    hint[22] = auxOct.emplace(tmpParent.getTopRightBack());
                    hint[23] = auxOct.emplace(tmpParent.getTopRightFront());
                    hint[24] = auxOct.emplace(tmpParent.getTopBack());
                    hint[25] = auxOct.emplace(tmpParent.getTopFront());
#endif

                    if(auxOct.size()>preSz)
                    {
                        curr=0;
                        pNodes.resize(pNodes.size()+auxOct.size()-preSz);
                        for(unsigned int kk=0;kk<neighbourCount;kk++)
                        {
                            if(hint[kk].second) {
                                pNodes[preSz + curr] = *(hint[kk].first);
                                //pNodes.push_back(*(hint[kk].first));
                                curr++;
                            }

                        }
                    }
                    //std::cout<<"aux Size: "<<auxOct.size()<<" pNodes size: "<<pNodes.size()<<std::endl;
                }


            }



        }


        template<typename T>
        void SFC_treeSort(T* pNodes , DendroIntL n ,std::vector<T>& pOutSorted,std::vector<T>& pOutConstruct,std::vector<T>& pOutBalanced, unsigned int pMaxDepthBit,unsigned int pMaxDepth, T& parent, unsigned int rot_id,unsigned int k, unsigned int options)
        {

            if(n==0) return;
            register unsigned int cnum;
            register unsigned int cnum_prev=0;
            //register unsigned int n=0;
            unsigned int rotation=0;
            DendroIntL count[(NUM_CHILDREN+2)]={};
            unsigned int lev=pMaxDepth-pMaxDepthBit;
            pMaxDepthBit--;
            unsigned int x,y,z;
            T temp;
            count[0]=0;

            for (DendroIntL i=0; i< n; ++i) {

              cnum = (lev < pNodes[i].getLevel())? 1 +( (((pNodes[i].getZ() >> pMaxDepthBit) & 1u) << 2u) | (((pNodes[i].getY() >> pMaxDepthBit) & 1u) << 1u) | ((pNodes[i].getX() >>pMaxDepthBit) & 1u)):0;
              count[cnum+1]++;

            }
            //std::cout<<"initial Count finished"<<std::endl;
            DendroIntL loc[NUM_CHILDREN+1];
            T unsorted[NUM_CHILDREN+1];
            unsigned int live = 0;


            for (unsigned int i=0; i<(NUM_CHILDREN+1); ++i) {
                if(i==0)
                {
                    loc[0]=count[0];
                    count[1]+=count[0];
                    unsorted[live] = pNodes[loc[0]];
                    if (loc[0] < count[1]) {live++; /*std::cout<<i<<" Live: "<<live<<std::endl;*/}
                }else
                {
                    cnum=(rotations[ROTATION_OFFSET * rot_id+ i-1] - '0');
                    (i>1) ? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id+i-2] - '0')+2): cnum_prev=1;
                    loc[cnum+1]=count[cnum_prev];
                    count[cnum+2] += count[cnum_prev];
                    //std::cout<<" loc[cnum+1]: "<<loc[cnum+1]<<std::endl;
                    (loc[cnum+1]==n) ? unsorted[live] = pNodes[loc[cnum+1]-1]: unsorted[live] = pNodes[loc[cnum+1]];
                    if (loc[cnum+1] < count[cnum+2]) {live++; /*std::cout<<i<<" Live after increment: "<<live<<std::endl;*/}
                }

            }
            //std::cout<<"second pass finished "<<std::endl;

            if(live>0)
            {

                live--;

                for (DendroIntL i=0; i < n ; ++i) {


                    cnum = (lev < unsorted[live].getLevel()) ? ((((unsorted[live].getZ() >> pMaxDepthBit) & 1u) << 2u) | (((unsorted[live].getY() >> pMaxDepthBit) & 1u) << 1u) | ((unsorted[live].getX() >> pMaxDepthBit) & 1u))+ 1: 0 ;

                    pNodes[loc[(cnum )]++] = unsorted[live];
                    (loc[cnum]==n) ? unsorted[live] = pNodes[loc[cnum]-1] : unsorted[live] = pNodes[loc[cnum]];
                    if ((loc[cnum] == count[cnum + 1])) {
                     /* if(i<(n-1)) assert(live>0);
                        if(live==0) assert(i==(n-1));*/
                          live--;

                    }


                }
            }

            if((options==TS_SORT_ONLY) || (options==TS_REMOVE_DUPLICATES))
            {
                if (pMaxDepthBit > 0) {

                    DendroIntL numElements=0;
                    for (unsigned int i=1; i<(NUM_CHILDREN+1); i++) {
                        cnum=(rotations[ROTATION_OFFSET*rot_id+i-1]-'0');
                        (i>1)? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id+i-2] - '0')+2) : cnum_prev=1;
                        numElements = count[cnum+2] - count[cnum_prev];
                        if (numElements > k) {
                            rotation=HILBERT_TABLE[NUM_CHILDREN * rot_id + cnum];
                            if(numElements)SFC_treeSort(pNodes+count[cnum_prev],numElements,pOutSorted,pOutConstruct,pOutBalanced,(pMaxDepthBit),pMaxDepth,temp,rotation,k,options);
                        }

                    }
                }

            }else
            {
                if (pMaxDepthBit > MAXDEAPTH_LEVEL_DIFF) {

                    DendroIntL numElements=0;
                    for (unsigned int i=1; i<(NUM_CHILDREN+1); i++) {
                        cnum=(rotations[ROTATION_OFFSET*rot_id+i-1]-'0');
                        (i>1)? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id+i-2] - '0')+2) : cnum_prev=1;
                        numElements = count[cnum+2] - count[cnum_prev];
                        if((options & TS_CONSTRUCT_OCTREE) | (options & TS_BALANCE_OCTREE))
                        {
                            x=parent.getX() +(((int)((bool)(cnum & 1u)))<<(pMaxDepthBit));
                            y=parent.getY() +(((int)((bool)(cnum & 2u)))<<(pMaxDepthBit));
                            z=parent.getZ() +(((int)((bool)(cnum & 4u)))<<(pMaxDepthBit));
                            temp=T(x,y,z,(lev+1),parent.getDim(),pMaxDepth);

                        }
                        if (numElements > k) {
                            rotation=HILBERT_TABLE[NUM_CHILDREN * rot_id + cnum];
                            SFC_treeSort(pNodes+count[cnum_prev],numElements,pOutSorted,pOutConstruct,pOutBalanced,(pMaxDepthBit),pMaxDepth,temp,rotation,k,options);
                        }else if((options & TS_CONSTRUCT_OCTREE) | (options & TS_BALANCE_OCTREE))
                        {
                            if(options & TS_CONSTRUCT_OCTREE) {
                                pOutConstruct.push_back(temp);
                            }

                            if (options & TS_BALANCE_OCTREE) {
                                //generate all the neighbours and add them to unordered map. use pOutBalaced to push balanced octree.
                                pOutBalanced.push_back(temp);


                            }

                        }

                    }
                }

            }



            if(lev==0)
            {

                // !!!! Note: Please note that all the code here executed only once. In the final stage of the recursion.

                if((options & TS_REMOVE_DUPLICATES)) {

#ifdef PROFILE_TREE_SORT
                    t1=std::chrono::high_resolution_clock::now();//MPI_Wtime();
#endif

                    // Note: This is executed only once. In the final stage of the recursion.
                    // Do the remove duplicates here.
                    if (n >= 1) {
                        std::vector<T> tmp(n);
                        T *tmpPtr = (&(*(tmp.begin())));

                        tmpPtr[0] = pNodes[0];

                        unsigned int tmpSize = 1;
                        unsigned int vecTsz = static_cast<unsigned int>(n);

                        for (unsigned int i = 1; i < n; i++) {
                            if ( /*(!tmpPtr[tmpSize-1].isAncestor(pNodes[i])) &*/  (tmpPtr[tmpSize - 1] != pNodes[i])) { // It is efficient to do this rather than marking all the elements in sorting. (Which will cause a performance degradation. )
                                tmpPtr[tmpSize] = pNodes[i];
                                tmpSize++;
                            }
                        }//end for


                        // Remove ancestor loop for removing local ancestors.
                        // Assumes that we have removed all the duplicates after the first iteration.

                        tmp.resize(tmpSize);
                        std::vector<T> tmp_rmvAncestors(tmp.size());
                        tmpPtr = (&(*(tmp_rmvAncestors.begin())));
                        tmpPtr[0]=tmp[0];
                        tmpSize=0;

                        for(unsigned int i=1;i<tmp.size();i++)
                        {
                            if(tmpPtr[tmpSize].isAncestor(tmp[i]))
                                tmpPtr[tmpSize]=tmp[i];
                            else {
                                tmpPtr[tmpSize+1]=tmp[i];
                                tmpSize++;
                            }

                        }
                        tmp_rmvAncestors.resize(tmpSize+1);
                        std::swap(pOutSorted, tmp_rmvAncestors);

                        tmp_rmvAncestors.clear();
                        tmp.clear();
                    }

#ifdef PROFILE_TREE_SORT
                    remove_duplicates_seq=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count();

#endif

                }

                if(options & TS_BALANCE_OCTREE)
                {
                    // Bottom up balancing octant creation.
                    //std::cout<<"balOCt: before Aux octants: "<<pOutBalanced.size()<<std::endl;
#ifdef PROFILE_TREE_SORT
                    t1=std::chrono::high_resolution_clock::now();//MPI_Wtime();
#endif
                    /*int rank;
                    MPI_Comm_rank(MPI_COMM_WORLD,&rank);*/

                    SFC::seqSort::SFC_bottomUpBalOctantCreation(pOutBalanced);

#ifdef PROFILE_TREE_SORT
                    auxBalOCt_time=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t1).count();

#endif
                    std::vector<T> tmpSorted;
                    std::vector<T> tmpConstruct;
                    std::vector<T> tmpBalanced;
                    T root =T(0,0,0,0,m_uiDim,pMaxDepth);
                    //std::cout<<"bal with aux octants: "<<pOutBalanced.size()<<std::endl;
                    SFC::seqSort::SFC_treeSort(&(*(pOutBalanced.begin())),pOutBalanced.size(),tmpSorted,tmpConstruct,tmpBalanced,pMaxDepth,pMaxDepth,root,0,k,2);
                    std::swap(tmpConstruct,pOutBalanced);
                    tmpConstruct.clear();
                }


            }

        } // end of function SFC_treeSort


        template<typename T>
        inline void SFC_bucketing(T *pNodes, int lev, unsigned int maxDepth,unsigned char rot_id,DendroIntL &begin, DendroIntL &end, DendroIntL *splitters)
        {


           /* int rank,npes;
            MPI_Comm_rank(MPI_COMM_WORLD,&rank);
            MPI_Comm_size(MPI_COMM_WORLD,&npes);

            MPI_Barrier(MPI_COMM_WORLD);*/

            //std::cout<<rank<<" Bucketing : begin: "<<begin<<" end: "<<end<<" lev: "<<lev<<" rotaion: "<<(int)rot_id<<std::endl;


            // Special case that needs to be handled when performing the load balancing.
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


            for (unsigned int  i=0; i<(NUM_CHILDREN+1); ++i) {
                if(i==0)
                {
                    loc[0]=count[0];
                    count[1]+=count[0];
                    unsorted[live] = pNodes[loc[0]];
                    if (loc[0] < count[1]) {live++; /*std::cout<<" level buck Live ++ : "<<live<<std::endl;*/}
                }else {

                    cnum = (rotations[ROTATION_OFFSET * rot_id + i-1] - '0');
                    (i > 1) ? cnum_prev = ((rotations[ROTATION_OFFSET * rot_id + i - 2] - '0') + 2) : cnum_prev = 1;

                    loc[cnum+1] = count[cnum_prev];
                    count[cnum + 2] += count[cnum_prev];
                    if (loc[cnum+1] < count[cnum + 2]) { unsorted[live] = pNodes[loc[(cnum+1)]]; live++; /*std::cout<<i<<" Live ++:  "<<live<<std::endl;*/}
                }
            }

            live--;
            //if(live <0) std::cout<<"live count overflow. " <<std::endl; // Note: This could be a potential bug.
            for (DendroIntL i=begin; i<end; ++i) {

                //cnum = (lev < unsorted[live].getLevel()) ? (((((unsorted[live].getZ() & (1u << mid_bit)) >> mid_bit) << 2u) | (((unsorted[live].getY() & (1u << mid_bit)) >> mid_bit) << 1u) | ((unsorted[live].getX() & (1u << mid_bit)) >> mid_bit)))+ 1: 0 ;
                cnum=(lev < unsorted[live].getLevel()) ? ((((unsorted[live].getZ() >> mid_bit) & 1u) << 2u) | (((unsorted[live].getY() >> mid_bit) & 1u) << 1u) | ((unsorted[live].getX() >> mid_bit) & 1u))+ 1: 0 ;
                if(loc[cnum]<count[cnum+1]) {

                    pNodes[loc[cnum]++] = unsorted[live];
                    if ((loc[cnum] == count[cnum + 1])) {
                        live--;
                        continue;
                    }
                    unsorted[live] = pNodes[loc[cnum]];

                }
                /*if(live<0)
                {
                    std::cout<<"begin: "<<begin<<" end: "<<end<<" i: "<<i<<" loc[0]: "<<loc[0]<<" live : "<<live<<std::endl;
                }*/

            }



        }// end of function SFC_bucketing.


        template<typename T>
        void SFC_treeSortLocalOptimal(T* pNodes , DendroIntL n , unsigned int pMaxDepthBit,unsigned int pMaxDepth, T& parent, unsigned int rot_id,bool minimum,T& optimal)
        {

            if(n==0 && minimum)
            {
                optimal=T(0,0,0,0,m_uiDim,m_uiMaxDepth);
                return ;
            }else if(n==0 && (!minimum))
            {
                optimal=T(((1u<<m_uiMaxDepth)-1),((1u<<m_uiMaxDepth)-1),((1u<<m_uiMaxDepth)-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth);
                return ;
            }

            unsigned int lev=pMaxDepth-pMaxDepthBit;
            DendroIntL nNodeBegin=0;
            DendroIntL nNodeEnd=n;

            unsigned int hindex = 0;
            unsigned int hindexN = 0;
            unsigned int index = 0;
            DendroIntL splitterNodes[(NUM_CHILDREN + 1)];
            unsigned int bucketIndex=0;
            if(!minimum) bucketIndex=(NUM_CHILDREN-1);

            if(((nNodeEnd-nNodeBegin)==1))
            {
                optimal=pNodes[nNodeBegin];
                return ;
            }

            if( (pMaxDepthBit) && ((nNodeEnd-nNodeBegin)>1)) {

                SFC::seqSort::SFC_bucketing(pNodes,lev,pMaxDepth,rot_id,nNodeBegin,nNodeEnd,splitterNodes);
                hindex = (rotations[2 * NUM_CHILDREN * rot_id + bucketIndex] - '0');
                index=HILBERT_TABLE[NUM_CHILDREN * rot_id + bucketIndex];
                (bucketIndex == (NUM_CHILDREN - 1)) ? hindexN = bucketIndex + 1: hindexN = (rotations[2 * NUM_CHILDREN * rot_id + bucketIndex + 1] - '0');

                assert(splitterNodes[hindex] <= splitterNodes[hindexN]);
                if(splitterNodes[hindex]!=splitterNodes[hindexN]) SFC_treeSortLocalOptimal(pNodes+splitterNodes[hindex],(splitterNodes[hindexN]-splitterNodes[hindex]),(pMaxDepthBit-1),pMaxDepth,parent,index,minimum,optimal);
                else
                {
                    //std::cout<<"Calling Seq: sort at Lev: "<<lev<<"size : "<<n<<std::endl;
                    std::vector<T>tmpVec;
                    SFC::seqSort::SFC_treeSort(pNodes,n,tmpVec,tmpVec,tmpVec,pMaxDepth,pMaxDepth,parent,ROOT_ROTATION,1,TS_SORT_ONLY);
                    (minimum) ? optimal=pNodes[0]: optimal=pNodes[n-1];
                    return ;
                }



            }



        }

        template<typename T>
        void SFC_treeSortLocalOptimal(T* pNodes , DendroIntL n , unsigned int pMaxDepthBit,unsigned int pMaxDepth, T& parent, unsigned int rot_id,T&min,T&max)
        {
            unsigned int lev=pMaxDepth-pMaxDepthBit;
            DendroIntL splitterNodes[(NUM_CHILDREN + 1)];
            unsigned int hindex = 0;
            unsigned int hindexN = 0;
            unsigned int index = 0;
            DendroIntL nNodeBegin=0;
            DendroIntL nNodeEnd=n;

            SFC::seqSort::SFC_bucketing(pNodes,lev,pMaxDepth,rot_id,nNodeBegin,nNodeEnd,splitterNodes);

            unsigned int min_b=0,max_b=NUM_CHILDREN-1;
            for(unsigned int i=0;i<NUM_CHILDREN;i++)
            {
                hindex = (rotations[2 * NUM_CHILDREN * rot_id + min_b] - '0');
                (min_b == (NUM_CHILDREN - 1)) ? hindexN = min_b + 1: hindexN = (rotations[2 * NUM_CHILDREN * rot_id + min_b + 1] - '0');
                if(splitterNodes[hindex]!=splitterNodes[hindexN])break;
                min_b++;


            }

            for(unsigned int i=0;i<(NUM_CHILDREN);i++)
            {
                hindex = (rotations[2 * NUM_CHILDREN * rot_id + max_b] - '0');
                (max_b == (NUM_CHILDREN - 1)) ? hindexN = max_b + 1: hindexN = (rotations[2 * NUM_CHILDREN * rot_id + max_b + 1] - '0');
                if(splitterNodes[hindex]!=splitterNodes[hindexN])break;
                max_b--;
            }


            hindex = (rotations[2 * NUM_CHILDREN * rot_id + min_b] - '0');
            (min_b == (NUM_CHILDREN - 1)) ? hindexN = min_b + 1: hindexN = (rotations[2 * NUM_CHILDREN * rot_id + min_b + 1] - '0');
            index=HILBERT_TABLE[NUM_CHILDREN * rot_id + min_b];
            DendroIntL numElem=splitterNodes[hindexN]-splitterNodes[hindex];
            SFC::seqSort::SFC_treeSortLocalOptimal(&pNodes[splitterNodes[hindex]],numElem,pMaxDepthBit-1,pMaxDepth,parent,index,true,min);
            //SFC::seqSort::SFC_treeSortLocalOptimal(pNodes,n,pMaxDepthBit,pMaxDepth,parent,rot_id,true,min);

            hindex = (rotations[2 * NUM_CHILDREN * rot_id + max_b] - '0');
            (max_b == (NUM_CHILDREN - 1)) ? hindexN = max_b + 1: hindexN = (rotations[2 * NUM_CHILDREN * rot_id + max_b + 1] - '0');
            index=HILBERT_TABLE[NUM_CHILDREN * rot_id + max_b];
            numElem=splitterNodes[hindexN]-splitterNodes[hindex];
            SFC::seqSort::SFC_treeSortLocalOptimal(&pNodes[splitterNodes[hindex]],numElem,pMaxDepthBit-1,pMaxDepth,parent,index,false,max);
            //SFC::seqSort::SFC_treeSortLocalOptimal(pNodes,n,pMaxDepthBit,pMaxDepth,parent,rot_id,false,max);







        }


        //========================================================= Function definition end =========================================================================================


    }
}

namespace SFC
{

    namespace parSort
    {

        //========================================================= Function declaration begin.=========================================================================================
        /**
             * @author Milinda Fernando
             * @breif Staged version of the splitter selection. This is used when the number of mpi tasks are high.
             * */
        template<typename T>
        inline void SFC_SplitterFix(std::vector<T>& pNodes,unsigned int pMaxDepth,double loadFlexibility,unsigned int sf_k,MPI_Comm comm,MPI_Comm * newComm);


        /**
             * @author Milinda Fernando
             * @breif Staged version of the splitter selection. This is used when the number of mpi tasks are high.
             * @detail
            * OPTIONS For TreeSort.
            *
             * NOTE: Consider these options as bits, in the little endian ordering,
             *
             * (2^0 bit location) TS_REMOVE_DUPLICATES : if selected ensures that the output of the algorithm is sorted and unique.
             * (2^1 bit location) TS_CONSTRUCT_OCTREE  : if selected ensures that the output of the algorithm is sorted and completed octree
             * (2^2 bit location) TS_BALANCE_OCTREE    : if selected ensures that the output of the algorithm is sorted completed and balanced.
             * */
        template <typename T>
        void SFC_treeSort(std::vector<T> &pNodes, std::vector<T>& pOutSorted,std::vector<T>& pOutConstruct,std::vector<T>& pOutBalanced , double loadFlexibility,unsigned int pMaxDepth, T& parent, unsigned int rot_id,unsigned int k, unsigned int options, unsigned int sf_k,MPI_Comm pcomm);


       /* template <typename T>
        void SFC_PartitionW(std::vector<T>&pNodes,double loadFlexibility, unsigned int maxDepth,MPI_Comm comm);*/


        //========================================================= Function declaration end.=========================================================================================



        //========================================================= Function definition begin.=========================================================================================

        template<typename T>
        inline void SFC_SplitterFix(std::vector<T>& pNodes,unsigned int pMaxDepth,double loadFlexibility,unsigned int sf_k,MPI_Comm comm,MPI_Comm * newComm) {

            #ifdef SPLITTER_SELECTION_FIX

            int rank, npes;
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &npes);



            //#pragma message ("Splitter selection FIX ON")

            if (npes > NUM_NPES_THRESHOLD) {

                unsigned int npes_sqrt = 1u << (binOp::fastLog2(npes) / 2);

                const unsigned int a=sf_k;
                const unsigned int b=npes/a;
                const unsigned int p_mod_a=npes%a;

                // if(!rank)
                //     std::cout<<" a : "<<a <<" b : "<<b<<" p_mod_a "<<p_mod_a<<std::endl;


                unsigned int dim = 3;
                #ifdef DIM_2
                    dim=2;
                #else
                    dim = 3;
                #endif

                unsigned int firstSplitLevel = std::ceil(binOp::fastLog2(a) / (double) dim);
                unsigned int totalNumBuckets = 1u << (dim * firstSplitLevel);


                DendroIntL localSz = pNodes.size();
                DendroIntL globalSz = 0;
                MPI_Allreduce(&localSz, &globalSz, 1, MPI_LONG_LONG, MPI_SUM, comm);

                // Number of initial buckets. This should be larger than npes.

                // maintain the splitters and buckets for splitting for further splitting.
                std::vector<DendroIntL> bucketCounts;
                std::vector<BucketInfo<T>> bucketInfo;  // Stores the buckets info of the buckets where the initial buckets was splitted.
                std::vector<DendroIntL> bucketSplitter;


                std::vector<BucketInfo<T>> nodeStack; // rotation id stack
                BucketInfo<T> root(0, 0, 0, pNodes.size());
                nodeStack.push_back(root);
                BucketInfo<T> tmp = root;
                unsigned int levSplitCount = 0;

                // Used repetitively  in rotation computations.
                unsigned int hindex = 0;
                unsigned int hindexN = 0;

                unsigned int index = 0;
                //bool *updateState = new bool[pNodes.size()];
                unsigned int numLeafBuckets = 0;

                unsigned int begin_loc = 0;


                DendroIntL spliterstemp[(NUM_CHILDREN + 1)];
                while (numLeafBuckets < totalNumBuckets) {

                    tmp = nodeStack[0];
                    nodeStack.erase(nodeStack.begin());

                    SFC::seqSort::SFC_bucketing(&(*(pNodes.begin())), tmp.lev, pMaxDepth, tmp.rot_id, tmp.begin,
                                                tmp.end, spliterstemp);

                    for (int i = 0; i < NUM_CHILDREN; i++) {
                        hindex = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i] - '0');
                        if (i == (NUM_CHILDREN - 1))
                            hindexN = i + 1;
                        else
                            hindexN = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i + 1] - '0');
                        assert(spliterstemp[hindex] <= spliterstemp[hindexN]);
                        index = HILBERT_TABLE[NUM_CHILDREN * tmp.rot_id + hindex];

                        BucketInfo<T> child(index, (tmp.lev + 1), spliterstemp[hindex], spliterstemp[hindexN]);
                        nodeStack.push_back(child);

                        if (tmp.lev == (firstSplitLevel - 1)) {
                            BucketInfo<T> bucket(index, (tmp.lev + 1), spliterstemp[hindex], spliterstemp[hindexN]);
                            bucketCounts.push_back((spliterstemp[hindexN] - spliterstemp[hindex]));
                            bucketSplitter.push_back(spliterstemp[hindex]);
                            bucketInfo.push_back(bucket);
                            numLeafBuckets++;

                        }

                    }


                }

                /*MPI_Barrier(MPI_COMM_WORLD);
                if(!rank)std::cout<<" intial buckets created: "<<std::endl;*/

                std::vector<DendroIntL> bucketCounts_g(bucketCounts.size());
                std::vector<DendroIntL> bucketCounts_gScan(bucketCounts.size());

                par::Mpi_Allreduce<DendroIntL>(&(*(bucketCounts.begin())), &(*(bucketCounts_g.begin())),
                                               bucketCounts.size(), MPI_SUM, comm);

                #ifdef DEBUG_TREE_SORT
                    assert(totalNumBuckets);
                #endif
                bucketCounts_gScan[0] = bucketCounts_g[0];
                for (int k = 1; k < bucketCounts_g.size(); k++) {
                    bucketCounts_gScan[k] = bucketCounts_gScan[k - 1] + bucketCounts_g[k];
                }


                #ifdef DEBUG_TREE_SORT
                    if(!rank) {
                        for (int i = 0; i < totalNumBuckets; i++) {
                            //std::cout << "Bucket count G : " << i << " : " << bucketCounts_g[i] << std::endl;
                            std::cout << "Bucket initial count scan G : " << i << " : " << bucketCounts_gScan[i] << std::endl;
                        }
                    }
                #endif
                DendroIntL *localSplitterTmp = new DendroIntL[a];
                DendroIntL idealLoadBalance = 0;
                begin_loc = 0;

                std::vector<unsigned int> splitBucketIndex;
                begin_loc = 0;
                for (int i = 0; i < a - 1; i++) {
                    idealLoadBalance += ((i + 1) * globalSz / a - i * globalSz / a);
                    DendroIntL toleranceLoadBalance = ((i + 1) * globalSz / a - i * globalSz / a) * loadFlexibility;

                    unsigned int loc = (
                            std::lower_bound(bucketCounts_gScan.begin(), bucketCounts_gScan.end(), idealLoadBalance) -
                            bucketCounts_gScan.begin());
                    //std::cout<<rank<<" Searching: "<<idealLoadBalance<<"found: "<<loc<<std::endl;

                    if (abs(bucketCounts_gScan[loc] - idealLoadBalance) > toleranceLoadBalance) {

                        if (splitBucketIndex.empty() || splitBucketIndex.back() != loc)
                            splitBucketIndex.push_back(loc);
                        /*if(!rank)
                          std::cout<<"Bucket index :  "<<loc << " Needs a split "<<std::endl;*/
                    } else {
                        if ((loc + 1) < bucketSplitter.size())
                            localSplitterTmp[i] = bucketSplitter[loc + 1];
                        else
                            localSplitterTmp[i] = bucketSplitter[loc];
                    }

                    /* if(loc+1<bucketCounts_gScan.size())
                         begin_loc=loc+1;
                     else
                         begin_loc=loc;*/

                }
                localSplitterTmp[a - 1] = pNodes.size();


                #ifdef DEBUG_TREE_SORT
                    for(int i=0;i<splitBucketIndex.size()-1;i++)
                    {
                        assert(pNodes[bucketSplitter[splitBucketIndex[i]]]<pNodes[bucketSplitter[splitBucketIndex[i+1]]]);
                    }
                #endif


                std::vector<DendroIntL> newBucketCounts;
                std::vector<DendroIntL> newBucketCounts_g;
                std::vector<BucketInfo<T>> newBucketInfo;
                std::vector<DendroIntL> newBucketSplitters;


                std::vector<DendroIntL> bucketCount_gMerge;
                std::vector<DendroIntL> bucketSplitterMerge;
                std::vector<BucketInfo<T>> bucketInfoMerge;

                DendroIntL splitterTemp[(NUM_CHILDREN + 1)];
                while (!splitBucketIndex.empty()) {


                    newBucketCounts.clear();
                    newBucketCounts_g.clear();
                    newBucketInfo.clear();
                    newBucketSplitters.clear();

                    bucketCount_gMerge.clear();
                    bucketSplitterMerge.clear();
                    bucketInfoMerge.clear();

                    if (bucketInfo[splitBucketIndex[0]].lev < pMaxDepth) {

                        BucketInfo<T> tmp;
                        // unsigned int numSplitBuckets = NUM_CHILDREN * splitBucketIndex.size();

                        #ifdef DEBUG_TREE_SORT
                        if (!rank)
                            for (int i = 0; i < splitBucketIndex.size(); i++)
                                std::cout << "Splitter Bucket Index: " << i << "  : " << splitBucketIndex[i] << std::endl;
                        #endif


                        //unsigned int maxDepthBuckets = 0;

                        #ifdef DEBUG_TREE_SORT
                        if(!rank) {
                            for (int i = 0; i <bucketSplitter.size(); i++) {
                                std::cout<<" Bucket Splitter : "<<i<<" : "<<bucketSplitter[i]<<std::endl;
                            }
                        }
                        #endif
                        for (int k = 0; k < splitBucketIndex.size(); k++) {
                            #ifdef DEBUG_TREE_SORT
                                if(!rank)
                                std::cout<<"Splitting Bucket index "<<splitBucketIndex[k]<<std::endl;
                            #endif

                            tmp = bucketInfo[splitBucketIndex[k]];
                            SFC::seqSort::SFC_bucketing(&(*(pNodes.begin())), tmp.lev, pMaxDepth, tmp.rot_id, tmp.begin,
                                                        tmp.end, splitterTemp);


                            for (int i = 0; i < NUM_CHILDREN; i++) {
                                hindex = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i] - '0');
                                if (i == (NUM_CHILDREN - 1))
                                    hindexN = i + 1;
                                else
                                    hindexN = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i + 1] - '0');

                                //newBucketCounts[NUM_CHILDREN * k + i] = (splitterTemp[hindexN] - splitterTemp[hindex]);
                                newBucketCounts.push_back((splitterTemp[hindexN] - splitterTemp[hindex]));

                                index = HILBERT_TABLE[NUM_CHILDREN * tmp.rot_id + hindex];
                                BucketInfo<T> bucket(index, (tmp.lev + 1), splitterTemp[hindex], splitterTemp[hindexN]);
                                //                          newBucketInfo[NUM_CHILDREN * k + i] = bucket;
                                //                          newBucketSplitters[NUM_CHILDREN * k + i] = splitterTemp[hindex];
                                newBucketInfo.push_back(bucket);
                                newBucketSplitters.push_back(splitterTemp[hindex]);
                                #ifdef DEBUG_TREE_SORT
                                    assert(pNodes[splitterTemp[hindex]]<=pNodes[splitterTemp[hindexN]]);
                                #endif
                            }

                        }

                        newBucketCounts_g.resize(newBucketCounts.size());
                        par::Mpi_Allreduce(&(*(newBucketCounts.begin())), &(*(newBucketCounts_g.begin())),
                                           newBucketCounts.size(), MPI_SUM, comm);
                        //MPI_Allreduce(&newBucketCounts[0], &newBucketCounts_g[0], newBucketCounts.size(), MPI_LONG_LONG, MPI_SUM, comm);


                        #ifdef DEBUG_TREE_SORT
                            for(int i=0;i<newBucketSplitters.size()-1;i++)
                            {
                                assert(pNodes[newBucketSplitters[i]]<=pNodes[newBucketSplitters[i+1]]);
                            }
                        #endif


                        #ifdef DEBUG_TREE_SORT
                        for (int k = 0; k < splitBucketIndex.size(); k++) {
                            unsigned int sum = 0;
                            for (int i = 0; i < NUM_CHILDREN; i++)
                                sum += newBucketCounts_g[NUM_CHILDREN * k + i];
                                if (!rank) if (bucketCounts_g[splitBucketIndex[k]] != sum) {
                                assert(false);
                            }

                        }
                        #endif

                        //int count = 0;
                        for (int k = 0; k < splitBucketIndex.size(); k++) {

                            unsigned int bucketBegin = 0;

                            (k == 0) ? bucketBegin = 0 : bucketBegin = splitBucketIndex[k - 1] + 1;


                            for (int i = bucketBegin; i < splitBucketIndex[k]; i++) {

                                bucketCount_gMerge.push_back(bucketCounts_g[i]);
                                bucketInfoMerge.push_back(bucketInfo[i]);
                                bucketSplitterMerge.push_back(bucketSplitter[i]);

                            }

                            tmp = bucketInfo[splitBucketIndex[k]];
                            for (int i = 0; i < NUM_CHILDREN; i++) {
                                bucketCount_gMerge.push_back(newBucketCounts_g[NUM_CHILDREN * k + i]);
                                bucketInfoMerge.push_back(newBucketInfo[NUM_CHILDREN * k + i]);
                                bucketSplitterMerge.push_back(newBucketSplitters[NUM_CHILDREN * k + i]);

                            }


                            if (k == (splitBucketIndex.size() - 1)) {
                                for (int i = splitBucketIndex[k] + 1; i < bucketCounts_g.size(); i++) {
                                    bucketCount_gMerge.push_back(bucketCounts_g[i]);
                                    bucketInfoMerge.push_back(bucketInfo[i]);
                                    bucketSplitterMerge.push_back(bucketSplitter[i]);
                                }


                            }


                        }

                        std::swap(bucketCounts_g, bucketCount_gMerge);
                        std::swap(bucketInfo, bucketInfoMerge);
                        std::swap(bucketSplitter, bucketSplitterMerge);


                        bucketCounts_gScan.resize(bucketCounts_g.size());

                        bucketCounts_gScan[0] = bucketCounts_g[0];
                        for (int k = 1; k < bucketCounts_g.size(); k++) {
                            bucketCounts_gScan[k] = bucketCounts_gScan[k - 1] + bucketCounts_g[k];
                        }
                        #ifdef DEBUG_TREE_SORT
                            assert(bucketCounts_gScan.back()==globalSz);
                        #endif
                        splitBucketIndex.clear();
                        #ifdef DEBUG_TREE_SORT
                            for(int i=0;i<bucketSplitter.size()-2;i++)
                            {
                                std::cout<<"Bucket Splitter : "<<bucketSplitter[i]<<std::endl;
                                assert( bucketSplitter[i+1]!=(pNodes.size()) && pNodes[bucketSplitter[i]]<=pNodes[bucketSplitter[i+1]]);
                            }
                        #endif

                        idealLoadBalance = 0;
                        //begin_loc=0;
                        for (unsigned int i = 0; i < a - 1; i++) {
                            idealLoadBalance += ((i + 1) * globalSz / a - i * globalSz / a);
                            DendroIntL toleranceLoadBalance = (((i + 1) * globalSz / a - i * globalSz / a) *
                                                               loadFlexibility);
                            unsigned int loc = (std::lower_bound(bucketCounts_gScan.begin(), bucketCounts_gScan.end(),
                                                                 idealLoadBalance) -
                                                bucketCounts_gScan.begin());

                            if (abs(bucketCounts_gScan[loc] - idealLoadBalance) > toleranceLoadBalance) {
                                if (splitBucketIndex.empty() || splitBucketIndex.back() != loc)
                                    splitBucketIndex.push_back(loc);


                            } else {
                                if ((loc + 1) < bucketSplitter.size())
                                    localSplitterTmp[i] = bucketSplitter[loc + 1];
                                else
                                    localSplitterTmp[i] = bucketSplitter[loc];

                            }

                            /* if(loc+1<bucketCounts_gScan.size())
                                 begin_loc=loc+1;
                             else
                                 begin_loc=loc;*/

                        }
                        localSplitterTmp[a - 1] = pNodes.size();

                    } else {
                        //begin_loc=0;
                        idealLoadBalance = 0;
                        for (unsigned int i = 0; i < a - 1; i++) {

                            idealLoadBalance += ((i + 1) * globalSz / a - i * globalSz / a);
                            //DendroIntL toleranceLoadBalance = ((i + 1) * globalSz / npes - i * globalSz / npes) * loadFlexibility;
                            unsigned int loc = (
                                    std::lower_bound(bucketCounts_gScan.begin(), bucketCounts_gScan.end(),
                                                     idealLoadBalance) -
                                    bucketCounts_gScan.begin());


                            if ((loc + 1) < bucketSplitter.size())
                                localSplitterTmp[i] = bucketSplitter[loc + 1];
                            else
                                localSplitterTmp[i] = bucketSplitter[loc];

                            /* if(loc+1<bucketCounts_gScan.size())
                                 begin_loc=loc+1;
                             else
                                 begin_loc=loc;*/

                        }
                        localSplitterTmp[a - 1] = pNodes.size();


                        break;
                    }

                }

                bucketCount_gMerge.clear();
                bucketCounts.clear();
                bucketCounts_gScan.clear();
                bucketInfoMerge.clear();
                bucketInfo.clear();
                bucketCounts_g.clear();
                bucketSplitter.clear();
                bucketSplitterMerge.clear();
                newBucketCounts.clear();
                newBucketCounts_g.clear();
                newBucketInfo.clear();
                newBucketSplitters.clear();


                #ifdef DEBUG_TREE_SORT
                    if(!rank) {
                        for (int i = 0; i <a; i++) {
                            //std::cout << "Bucket count G : " << i << " : " << bucketCounts_g[i] << std::endl;
                            std::cout << "Local Splitter Tmp : " << i << " : " << localSplitterTmp[i] << std::endl;
                        }
                    }
                #endif


                std::vector<unsigned int> blockCounts;
                blockCounts.resize(a,b);

                for(unsigned int k=0;k<p_mod_a;k++)
                    blockCounts[k]=(b+1);




                std::vector<unsigned int > blockOffset;
                blockOffset.resize(a,0);

                blockOffset[0]=0;
                omp_par::scan(&(*(blockCounts.begin())),&(*(blockOffset.begin())),a);



                // compute the block ids.
                unsigned int blk_id;
                for(unsigned int k=0;k<a;k++)
                {
                    if( (rank>=blockOffset[k]) && rank<(blockOffset[k]+blockCounts[k]))
                    {
                        (blockOffset[k]!=0) ? blk_id=rank%blockOffset[k] : blk_id=rank;
                    }
                }

               /* if(!rank)
                    for(unsigned int k=0;k<a;k++)
                        std::cout<<" bck count: "<<k<<" val: "<<blockCounts[k]<<" offsets: "<<blockOffset[k]<<std::endl;*/


                //std::cout<<"rank: "<<rank<<" block_id: "<<blk_id<<std::endl;


                unsigned int *sendCnt = new unsigned int[a];
                unsigned int owner_blk=0;
                (rank<(b+1)*p_mod_a) ? owner_blk=rank/(b+1) : owner_blk =(((rank-(b+1)*p_mod_a)/b) + p_mod_a);

                //MPI_Barrier(comm);
		        //std::cout<<" rank: "<<rank<<" owner_blk: "<<owner_blk<<std::endl;

                sendCnt[0] = localSplitterTmp[0];
                for (unsigned int i = 1; i < a; i++)
                {
                    assert(localSplitterTmp[i]>=localSplitterTmp[i-1]);
                    sendCnt[i]= (localSplitterTmp[i] - localSplitterTmp[i - 1]);
                }


                // computed a-splitters.

                std::vector<unsigned int> send_count;
                std::vector<unsigned int> send_offset;
                std::vector<unsigned int> recv_count;
                std::vector<unsigned int> recv_offset;

                send_count.resize(npes,0);
                send_offset.resize(npes,0);
                recv_count.resize(npes,0);
                recv_offset.resize(npes,0);

                /*if(rank==4)
                    for(unsigned int p=0;p<a;p++)
                        std::cout<<" sendCnt["<<p<<" ]: "<<sendCnt[p]<<" spliter: ["<<p<<" ] : value: "<<localSplitterTmp[p]<<std::endl;*/

                for(unsigned int k=0;k<a;k++)
                {
                    if(owner_blk==k) continue;
                    if( (rank<(b+1)*p_mod_a) && blockCounts[k]==(b+1))
                    {
                        assert((blockOffset[k]+blk_id)<npes);
                        if(k>0)send_offset[(blockOffset[k]+blk_id)]=localSplitterTmp[k-1];
                        send_count[(blockOffset[k]+blk_id)]=sendCnt[k];
                    }else
                    {
                        assert((blockOffset[k]+(blk_id%b))<npes);
                        if(k>0)send_offset[(blockOffset[k]+(blk_id%b))]=localSplitterTmp[k-1];
                        send_count[(blockOffset[k]+(blk_id%b))]=sendCnt[k];
                    }

                }


                par::Mpi_Alltoall(&(*(send_count.begin())),&(*(recv_count.begin())),1,comm);

                recv_offset[0]=0;
                omp_par::scan(&(*(recv_count.begin())),&(*(recv_offset.begin())),npes);

                std::vector<T> recvBuf;
                recvBuf.resize(recv_offset[npes-1]+recv_count[npes-1]);


                par::Mpi_Alltoallv_sparse(&(*(pNodes.begin())),(int *) &(*(send_count.begin())),(int *) &(*(send_offset.begin())),&(*(recvBuf.begin())),(int *) &(*(recv_count.begin())),(int *) &(*(recv_offset.begin())),comm);
                (owner_blk>0) ? recvBuf.insert(recvBuf.end(),pNodes.begin()+localSplitterTmp[owner_blk-1],pNodes.begin()+localSplitterTmp[owner_blk-1]+sendCnt[owner_blk]) : recvBuf.insert(recvBuf.end(),pNodes.begin(),pNodes.begin()+sendCnt[owner_blk]);


                delete [] sendCnt;

                std::swap(pNodes,recvBuf);
                recvBuf.clear();



                #ifdef DEBUG_TREE_SORT
                if(!rank) {
                    for (int i = 0; i <a; i++) {
                        //std::cout << "Bucket count G : " << i << " : " << bucketCounts_g[i] << std::endl;
                        std::cout << "Send Cnt : " << i << " : " << sendCnt[i] << std::endl;
                    }
                }
                #endif



                #ifdef PROFILE_TREE_SORT
                    t1=std::chrono::high_resolution_clock::now();//MPI_Wtime();
                #endif



                unsigned int col=owner_blk;
                MPI_Comm_split(comm,col,rank,newComm);
		        // MPI_Barrier(comm);
                // if(!rank) std::cout<<" comm split works "<<std::endl;




            }

            #endif



            return;

        }// end of function.



        template <typename T>
        void SFC_treeSort(std::vector<T> &pNodes, std::vector<T>& pOutSorted,std::vector<T>& pOutConstruct,std::vector<T>& pOutBalanced , double loadFlexibility,unsigned int pMaxDepth, T& parent, unsigned int rot_id,unsigned int k, unsigned int options, unsigned int sf_k,MPI_Comm pcomm)
        {

            int rank, npes;
            MPI_Comm_rank(pcomm, &rank);
            MPI_Comm_size(pcomm, &npes);

            MPI_Comm comm;
            MPI_Comm_dup(pcomm,&comm);

            #ifdef PROFILE_TREE_SORT
                stats.clear();

                splitter_fix_all=0;
                splitter_time=0;
                all2all1_time=0;
                all2all2_time=0;
                localSort_time=0;
                remove_duplicates_seq=0;
                remove_duplicates_par=0;
                auxBalOCt_time=0;
                construction_time=0;
                balancing_time=0;
                total_rd=0;

                //MPI_Barrier(pcomm);

                t4=std::chrono::high_resolution_clock::now();//MPI_Wtime();
                t2=std::chrono::high_resolution_clock::now();//MPI_Wtime();

            #endif


            //SFC_SplitterFix(pNodes,pMaxDepth,loadFlexibility,pcomm,&comm);


            MPI_Comm SF_comm=pcomm;
            unsigned int SF_Stages=0;

            #ifdef PROFILE_TREE_SORT
                double * sf_full;
                double * sf_all2all;
                double * sf_splitters;
            #endif


            if(npes==1)
            {
                MPI_Comm_free(&comm);
                //call the sequential case
                SFC::seqSort::SFC_treeSort(&(*(pNodes.begin())),pNodes.size(),pOutSorted,pOutConstruct,pOutBalanced,pMaxDepth,pMaxDepth,parent,rot_id,k,options);
                return ;

            }



            #ifdef SPLITTER_SELECTION_FIX
                if(npes > NUM_NPES_THRESHOLD)
                {

                    SF_Stages= std::ceil((binOp::fastLog2(npes)/(double)binOp::fastLog2(sf_k))) - 1;
                    //if(!rank) std::cout<<" sf_stages: "<<SF_Stages<<std::endl;

                    #ifdef PROFILE_TREE_SORT
                        sf_full=new double[SF_Stages];
                        sf_all2all=new double[SF_Stages];
                        sf_splitters=new double[SF_Stages];
                    #endif

                    MPI_Comm * sf_comms=new MPI_Comm[SF_Stages];
                    unsigned int sf_i_break;
                    for(int i=0;i<SF_Stages;i++)
                    {
	               	    int sf_npes;
			            MPI_Comm_size(SF_comm,&sf_npes);	
                        //if(!rank) std::cout<<"sf_stage: "<<i<<"of "<<SF_Stages<<" sf_comm size "<<sf_npes<<std::endl;
	
                        if(sf_npes > NUM_NPES_THRESHOLD)
                        {
                            #ifdef PROFILE_TREE_SORT
                                t5_sf_staged=std::chrono::high_resolution_clock::now();
                            #endif

                            SFC_SplitterFix(pNodes,pMaxDepth,loadFlexibility,sf_k,SF_comm,&sf_comms[i]);
                                    
                            #ifdef PROFILE_TREE_SORT
                                sf_full[i]=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t5_sf_staged).count();
                                sf_all2all[i]=all2all1_time;
                                sf_splitters[i]=sf_full[i]-sf_all2all[i];
                            #endif

                            sf_i_break=i;                                    
                            SF_comm=sf_comms[i];
                            // MPI_Barrier(comm);
                            // int r1,n1;
                            // MPI_Comm_rank(SF_comm,&r1);
                            // MPI_Comm_size(SF_comm,&n1);
                            // if(!r1) std::cout<<" SF_comm size: "<<n1<<std::endl;
                        
                        }
			
                    }
                    
                    MPI_Comm_free(&comm);
                    MPI_Comm_dup(SF_comm,&comm);
                    for(int i=0;i<(sf_i_break+1);i++)
                    {
                        MPI_Comm_free(&sf_comms[i]);
                    }
                    
                    delete [] sf_comms;
                    
                }
            #endif

            #ifdef PROFILE_TREE_SORT
                splitter_fix_all=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t2).count();
            #endif

            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &npes);

            unsigned int dim=3;
            #ifdef DIM_2
                dim=2;
            #else
                dim=3;
            #endif

            #ifdef PROFILE_TREE_SORT
                //MPI_Barrier(pcomm);
                t2=std::chrono::high_resolution_clock::now();//MPI_Wtime();
            #endif
            
            unsigned int firstSplitLevel = std::ceil(binOp::fastLog2(npes)/(double)(dim));
            unsigned int totalNumBuckets =1u << (dim * firstSplitLevel);
            DendroIntL localSz=pNodes.size();
            DendroIntL globalSz=0;
            MPI_Allreduce(&localSz,&globalSz,1,MPI_LONG_LONG,MPI_SUM,comm);
            //if(!rank) std::cout<<"First Split Level : "<<firstSplitLevel<<" Total number of buckets: "<<totalNumBuckets <<std::endl;
            //if(!rank) std::cout<<"NUM_CHILDREN: "<<NUM_CHILDREN<<std::endl;
            // Number of initial buckets. This should be larger than npes.

            // maintain the splitters and buckets for splitting for further splitting.
            std::vector<DendroIntL> bucketCounts;
            std::vector<BucketInfo<T>> bucketInfo;  // Stores the buckets info of the buckets where the initial buckets was splitted.
            std::vector<DendroIntL > bucketSplitter;

            std::vector<BucketInfo<T>> nodeStack; // rotation id stack
            BucketInfo<T> root(0, 0, 0, pNodes.size());
            nodeStack.push_back(root);
            BucketInfo<T> tmp = root;
            unsigned int levSplitCount = 0;

            // Used repetitively  in rotation computations.
            unsigned int hindex = 0;
            unsigned int hindexN = 0;

            unsigned int index = 0;
            //bool *updateState = new bool[pNodes.size()];
            unsigned int numLeafBuckets =0;

            unsigned int begin_loc=0;
            /*
           * We need to split the array into buckets until it is we number of buckets is slightly higher than the number of processors.
           */


            // 1. ===================Initial Splitting Start===============================
            DendroIntL spliterstemp[(NUM_CHILDREN+1)];
            while(numLeafBuckets<totalNumBuckets) {

                tmp = nodeStack[0];
                nodeStack.erase(nodeStack.begin());


                SFC::seqSort::SFC_bucketing(&(*(pNodes.begin())),tmp.lev,pMaxDepth,tmp.rot_id,tmp.begin,tmp.end,spliterstemp);


                for (int i = 0; i < NUM_CHILDREN; i++) {
                    hindex = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i] - '0');
                    if (i == (NUM_CHILDREN-1))
                        hindexN = i + 1;
                    else
                        hindexN = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i + 1] - '0');
                    assert(spliterstemp[hindex] <= spliterstemp[hindexN]);
                    index = HILBERT_TABLE[NUM_CHILDREN * tmp.rot_id + hindex];

                    BucketInfo<T> child(index, (tmp.lev + 1), spliterstemp[hindex], spliterstemp[hindexN]);
                    nodeStack.push_back(child);


                    if(tmp.lev==(firstSplitLevel-1))
                    {
                        BucketInfo<T> bucket(index, (tmp.lev + 1), spliterstemp[hindex], spliterstemp[hindexN]);
                        bucketCounts.push_back((spliterstemp[hindexN] - spliterstemp[hindex]));
                        bucketSplitter.push_back(spliterstemp[hindex]);
                        bucketInfo.push_back(bucket);
                        numLeafBuckets++;

                    }

                }


            }



            #ifdef DEBUG_TREE_SORT
                for(int i=0;i<bucketSplitter.size()-2;i++)
                {
                    std::cout<<"Bucket Splitter : "<<bucketSplitter[i]<<std::endl;
                    //assert(pNodes[bucketSplitter[i]]<=pNodes[bucketSplitter[i+1]]);
                }
            #endif


            // #ifdef DEBUG_TREE_SORT
            //std::cout<<rank<<" Initial Splitter Calculation ended "<<numLeafBuckets<<std::endl;
            //#endif

            //1=================== Initial Splitting END=========================================================================


            // (2) =================== Pick npes splitters form the bucket splitters.
            std::vector<DendroIntL >bucketCounts_g(bucketCounts.size());
            std::vector<DendroIntL >bucketCounts_gScan(bucketCounts.size());


            par::Mpi_Allreduce<DendroIntL>(&(*(bucketCounts.begin())),&(*(bucketCounts_g.begin())),bucketCounts.size(),MPI_SUM,comm);

            /* if(!rank)
                 std::cout<<"All Reduction Time for : "<<npes<<"  : "<<allReduceTime<<std::endl;*/

            //MPI_Allreduce(&bucketCounts[0], &bucketCounts_g[0], bucketCounts.size(), MPI_LONG_LONG, MPI_SUM, comm);
            //std::cout<<"All to all ended. "<<rank<<std::endl;

            #ifdef DEBUG_TREE_SORT
                assert(totalNumBuckets);
            #endif
            bucketCounts_gScan[0]=bucketCounts_g[0];
            for(int k=1;k<bucketCounts_g.size();k++){
                bucketCounts_gScan[k]=bucketCounts_gScan[k-1]+bucketCounts_g[k];
            }

            #ifdef DEBUG_TREE_SORT
                if(!rank) {
                    for (int i = 0; i < totalNumBuckets; i++) {
                        //std::cout << "Bucket count G : " << i << " : " << bucketCounts_g[i] << std::endl;
                        std::cout << "Bucket initial count scan G : " << i << " : " << bucketCounts_gScan[i] << std::endl;
                    }
                }
            #endif

            #ifdef DEBUG_TREE_SORT
                assert(bucketCounts_gScan.back()==globalSz);
            #endif
            DendroIntL* localSplitter=new DendroIntL[npes];
            std::vector<unsigned int> splitBucketIndex;
            DendroIntL idealLoadBalance=0;
            //begin_loc=0;
            for(int i=0;i<npes-1;i++) {
                idealLoadBalance+=((i+1)*globalSz/npes -i*globalSz/npes);
                DendroIntL toleranceLoadBalance = ((i+1)*globalSz/npes -i*globalSz/npes) * loadFlexibility;

                unsigned int  loc=(std::lower_bound(bucketCounts_gScan.begin(), bucketCounts_gScan.end(), idealLoadBalance) - bucketCounts_gScan.begin());
                //std::cout<<rank<<" Searching: "<<idealLoadBalance<<"found: "<<loc<<std::endl;

                if(abs(bucketCounts_gScan[loc]-idealLoadBalance) > toleranceLoadBalance)
                {

                    if(splitBucketIndex.empty()  || splitBucketIndex.back()!=loc)
                        splitBucketIndex.push_back(loc);
                    /*if(!rank)
                      std::cout<<"Bucket index :  "<<loc << " Needs a split "<<std::endl;*/
                }else
                {
                    if ((loc + 1) < bucketSplitter.size())
                        localSplitter[i] = bucketSplitter[loc + 1];
                    else
                        localSplitter[i] = bucketSplitter[loc];
                }

                /* if(loc+1<bucketCounts_gScan.size())
                     begin_loc=loc+1;
                 else
                     begin_loc=loc;*/

            }

            localSplitter[npes-1]=pNodes.size();


            #ifdef DEBUG_TREE_SORT
                for(int i=0;i<splitBucketIndex.size()-1;i++)
                {
                    assert(pNodes[bucketSplitter[splitBucketIndex[i]]]<pNodes[bucketSplitter[splitBucketIndex[i+1]]]);
                }
            #endif
            std::vector<DendroIntL> newBucketCounts;
            std::vector<DendroIntL> newBucketCounts_g;
            std::vector<BucketInfo<T>> newBucketInfo;
            std::vector<DendroIntL> newBucketSplitters;


            std::vector<DendroIntL> bucketCount_gMerge;
            std::vector<DendroIntL> bucketSplitterMerge;
            std::vector<BucketInfo<T>> bucketInfoMerge;

            DendroIntL splitterTemp[(NUM_CHILDREN+1)];
            while(!splitBucketIndex.empty()) {


                newBucketCounts.clear();
                newBucketCounts_g.clear();
                newBucketInfo.clear();
                newBucketSplitters.clear();

                bucketCount_gMerge.clear();
                bucketSplitterMerge.clear();
                bucketInfoMerge.clear();

                if (bucketInfo[splitBucketIndex[0]].lev < pMaxDepth) {

                    BucketInfo<T> tmp;
                    // unsigned int numSplitBuckets = NUM_CHILDREN * splitBucketIndex.size();

                    #ifdef DEBUG_TREE_SORT
                        if (!rank)
                        for (int i = 0; i < splitBucketIndex.size(); i++)
                            std::cout << "Splitter Bucket Index: " << i << "  : " << splitBucketIndex[i] << std::endl;
                    #endif

                    //unsigned int maxDepthBuckets = 0;

                    #ifdef DEBUG_TREE_SORT
                        if(!rank) {
                                for (int i = 0; i <bucketSplitter.size(); i++) {
                                    std::cout<<" Bucket Splitter : "<<i<<" : "<<bucketSplitter[i]<<std::endl;
                                }
                        }
                    #endif
                    for (int k = 0; k < splitBucketIndex.size(); k++) {

                        tmp = bucketInfo[splitBucketIndex[k]];
                        #ifdef DEBUG_TREE_SORT
                            if(!rank)
                            std::cout<<"Splitting Bucket index "<<splitBucketIndex[k]<<"begin: "<<tmp.begin <<" end: "<<tmp.end<<" rot_id: "<< (int)tmp.rot_id<<std::endl;
                        #endif

                        SFC::seqSort::SFC_bucketing(&(*(pNodes.begin())),tmp.lev,pMaxDepth,tmp.rot_id,tmp.begin,tmp.end,splitterTemp);



                        for (int i = 0; i < NUM_CHILDREN; i++) {
                            hindex = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i] - '0');
                            if (i == (NUM_CHILDREN-1))
                                hindexN = i + 1;
                            else
                                hindexN = (rotations[2 * NUM_CHILDREN * tmp.rot_id + i + 1] - '0');

                            //newBucketCounts[NUM_CHILDREN * k + i] = (splitterTemp[hindexN] - splitterTemp[hindex]);
                            newBucketCounts.push_back((splitterTemp[hindexN] - splitterTemp[hindex]));

                            index = HILBERT_TABLE[NUM_CHILDREN * tmp.rot_id + hindex];
                            BucketInfo<T> bucket(index, (tmp.lev + 1), splitterTemp[hindex], splitterTemp[hindexN]);
                            //newBucketInfo[NUM_CHILDREN * k + i] = bucket;
                            //newBucketSplitters[NUM_CHILDREN * k + i] = splitterTemp[hindex];
                            newBucketInfo.push_back(bucket);
                            newBucketSplitters.push_back(splitterTemp[hindex]);
                            #ifdef DEBUG_TREE_SORT
                                assert(pNodes[splitterTemp[hindex]]<=pNodes[splitterTemp[hindexN]]);
                            #endif
                        }

                    }

                    newBucketCounts_g.resize(newBucketCounts.size());
                    par::Mpi_Allreduce(&(*(newBucketCounts.begin())),&(*(newBucketCounts_g.begin())),newBucketCounts.size(),MPI_SUM,comm);



                    #ifdef DEBUG_TREE_SORT
                    for(int i=0;i<newBucketSplitters.size()-1;i++)
                       {
                           assert(pNodes[newBucketSplitters[i]]<=pNodes[newBucketSplitters[i+1]]);
                       }
                    #endif



                    #ifdef DEBUG_TREE_SORT

                        for (int k = 0; k < splitBucketIndex.size(); k++) {
                            unsigned int sum = 0;
                            for (int i = 0; i < NUM_CHILDREN; i++)
                                sum += newBucketCounts_g[NUM_CHILDREN * k + i];
                            if (!rank) if (bucketCounts_g[splitBucketIndex[k]] != sum) {
                            assert(false);
                            }

                        }
                    #endif

                    //int count = 0;
                    for (int k = 0; k < splitBucketIndex.size(); k++) {

                        unsigned int bucketBegin = 0;

                        (k == 0) ? bucketBegin = 0 : bucketBegin = splitBucketIndex[k - 1] + 1;


                        for (int i = bucketBegin; i < splitBucketIndex[k]; i++) {

                            bucketCount_gMerge.push_back(bucketCounts_g[i]);
                            bucketInfoMerge.push_back(bucketInfo[i]);
                            bucketSplitterMerge.push_back(bucketSplitter[i]);

                        }

                        tmp = bucketInfo[splitBucketIndex[k]];
                        for (int i = 0; i < NUM_CHILDREN; i++) {
                            bucketCount_gMerge.push_back(newBucketCounts_g[NUM_CHILDREN * k + i]);
                            bucketInfoMerge.push_back(newBucketInfo[NUM_CHILDREN * k + i]);
                            bucketSplitterMerge.push_back(newBucketSplitters[NUM_CHILDREN * k + i]);

                        }


                        if (k == (splitBucketIndex.size() - 1)) {
                            for (int i = splitBucketIndex[k] + 1; i < bucketCounts_g.size(); i++) {
                                bucketCount_gMerge.push_back(bucketCounts_g[i]);
                                bucketInfoMerge.push_back(bucketInfo[i]);
                                bucketSplitterMerge.push_back(bucketSplitter[i]);
                            }


                        }


                    }

                    std::swap(bucketCounts_g,bucketCount_gMerge);
                    std::swap(bucketInfo,bucketInfoMerge);
                    std::swap(bucketSplitter,bucketSplitterMerge);


                    bucketCounts_gScan.resize(bucketCounts_g.size());

                    bucketCounts_gScan[0] = bucketCounts_g[0];
                    for (int k = 1; k < bucketCounts_g.size(); k++) {
                        bucketCounts_gScan[k] = bucketCounts_gScan[k - 1] + bucketCounts_g[k];
                    }
                    
                    #ifdef DEBUG_TREE_SORT
                        assert(bucketCounts_gScan.back()==globalSz);
                    #endif
                    splitBucketIndex.clear();
                    #ifdef DEBUG_TREE_SORT
                        for(int i=0;i<bucketSplitter.size()-2;i++)
                        {
                            std::cout<<"Bucket Splitter : "<<bucketSplitter[i]<<std::endl;
                            assert( bucketSplitter[i+1]!=(pNodes.size()) && pNodes[bucketSplitter[i]]<=pNodes[bucketSplitter[i+1]]);
                        }
                    #endif

                    idealLoadBalance = 0;
                    //begin_loc=0;
                    for (unsigned int i = 0; i < npes-1; i++) {
                        idealLoadBalance += ((i + 1) * globalSz / npes - i * globalSz / npes);
                        DendroIntL toleranceLoadBalance = (((i + 1) * globalSz / npes - i * globalSz / npes) * loadFlexibility);
                        unsigned int loc = (std::lower_bound(bucketCounts_gScan.begin(), bucketCounts_gScan.end(), idealLoadBalance) -
                                            bucketCounts_gScan.begin());

                        if (abs(bucketCounts_gScan[loc] - idealLoadBalance) > toleranceLoadBalance) {
                            if (splitBucketIndex.empty() || splitBucketIndex.back() != loc)
                                splitBucketIndex.push_back(loc);


                        } else {
                            if ((loc + 1) < bucketSplitter.size())
                                localSplitter[i] = bucketSplitter[loc + 1];
                            else
                                localSplitter[i] = bucketSplitter[loc];

                        }

                        /* if(loc+1<bucketCounts_gScan.size())
                             begin_loc=loc+1;
                         else
                             begin_loc=loc;*/

                    }
                    localSplitter[npes-1]=pNodes.size();

                } else {
                    //begin_loc=0;
                    idealLoadBalance = 0;
                    for (unsigned int i = 0; i < npes-1; i++) {

                        idealLoadBalance += ((i + 1) * globalSz / npes - i * globalSz / npes);
                        //DendroIntL toleranceLoadBalance = ((i + 1) * globalSz / npes - i * globalSz / npes) * loadFlexibility;
                        unsigned int loc = (
                                std::lower_bound(bucketCounts_gScan.begin(), bucketCounts_gScan.end(), idealLoadBalance) -
                                bucketCounts_gScan.begin());


                        if ((loc + 1) < bucketSplitter.size())
                            localSplitter[i] = bucketSplitter[loc + 1];
                        else
                            localSplitter[i] = bucketSplitter[loc];

                        /* if(loc+1<bucketCounts_gScan.size())
                             begin_loc=loc+1;
                         else
                             begin_loc=loc;*/

                    }
                    localSplitter[npes-1]=pNodes.size();



                    break;
                }

            }



            bucketCount_gMerge.clear();
            bucketCounts.clear();
            bucketCounts_gScan.clear();
            bucketInfoMerge.clear();
            bucketInfo.clear();
            bucketCounts_g.clear();
            bucketSplitter.clear();
            bucketSplitterMerge.clear();
            newBucketCounts.clear();
            newBucketCounts_g.clear();
            newBucketInfo.clear();
            newBucketSplitters.clear();


            //#ifdef DEBUG_TREE_SORT
            // if(!rank) std::cout<<"Splitter Calculation ended "<<std::endl;
            //#endif

            #ifdef DEBUG_TREE_SORT
                for(int i=0;i<npes;i++)
                {
                    for(int j=i+1 ;j<npes -1;j++)
                        assert(pNodes[localSplitter[i]]<=pNodes[localSplitter[j]]);
                }
            #endif


            #ifdef DEBUG_TREE_SORT
            if(!rank)
            {
                for(int i=0;i<npes;i++)
                    std::cout<<"Rank "<<rank<<" Local Splitter: "<<i<<": "<<localSplitter[i]<<std::endl;
            }
            #endif

            // 3. All to all communication

            #ifdef PROFILE_TREE_SORT
                splitter_time=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t2).count();
                //MPI_Barrier(pcomm);
                t2=std::chrono::high_resolution_clock::now();//MPI_Wtime();
            #endif


            int * sendCounts = new  int[npes];
            int * recvCounts = new  int[npes];


            sendCounts[0] = localSplitter[0];

            for(unsigned int i=1;i<npes; ++i)
            {
                sendCounts[i] = localSplitter[i] - localSplitter[i-1];
            }



            par::Mpi_Alltoall(sendCounts,recvCounts,1,comm);
            //MPI_Alltoall(sendCounts, 1, MPI_INT,recvCounts,1,MPI_INT,comm);
            //std::cout<<"rank "<<rank<<" MPI_ALL TO ALL END"<<std::endl;

            int * sendDispl =new  int [npes];
            int * recvDispl =new  int [npes];


            sendDispl[0] = 0;
            recvDispl[0] = 0;

            for(int i=1;i<npes;i++)
            {
                sendDispl[i] = sendCounts[i-1] + sendDispl[i - 1];
                recvDispl[i] =recvCounts[i-1] +recvDispl[i-1];
            }


            #ifdef DEBUG_TREE_SORT
            /*if (!rank)*/ std::cout << rank << " : send = " << sendCounts[0] << ", " << sendCounts[1] << std::endl;
            /*if (!rank)*/ std::cout << rank << " : recv = " << recvCounts[0] << ", " << recvCounts[1] << std::endl;

            /* if (!rank) std::cout << rank << " : send offset  = " << sendDispl[0] << ", " << sendDispl[1] << std::endl;
            if (!rank) std::cout << rank << " : recv offset  = " << recvDispl[0] << ", " << recvDispl[1] << std::endl;*/
            #endif



            std::vector<T> pNodesRecv;
            DendroIntL recvTotalCnt=recvDispl[npes-1]+recvCounts[npes-1];
            if(recvTotalCnt) pNodesRecv.resize(recvTotalCnt);

            //par::Mpi_Alltoallv(&pNodes[0],sendCounts,sendDispl,&pNodesRecv[0],recvCounts,recvDispl,comm);
            //MPI_Alltoallv(&pNodes[0],sendCounts,sendDispl,MPI_TREENODE,&pNodesRecv[0],recvCounts,recvDispl,MPI_TREENODE,comm);
            par::Mpi_Alltoallv_sparse(&pNodes[0],sendCounts,sendDispl,&pNodesRecv[0],recvCounts,recvDispl,comm);
            // MPI_Request req;
            // MPI_Ialltoallv(&pNodes[0],sendCounts,sendDispl,par::Mpi_datatype<ot::TreeNode>::value(),&pNodesRecv[0],recvCounts,recvDispl,par::Mpi_datatype<ot::TreeNode>::value(),comm,&req);
            // MPI_Wait(&req,MPI_STATUS_IGNORE);

            #ifdef PROFILE_TREE_SORT
                all2all2_time=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t2).count();
                //MPI_Barrier(pcomm);
            #endif

            #ifdef DEBUG_TREE_SORT
                if(!rank) std::cout<<"All2All Communication Ended "<<std::endl;
            #endif

            delete[](localSplitter);
            localSplitter=NULL;


            pNodes.clear();
            //pNodes=pNodesRecv;
            std::swap(pNodes,pNodesRecv);
            pNodesRecv.clear();

            delete[](sendCounts);
            delete[](sendDispl);
            delete[](recvCounts);
            delete[](recvDispl);


            //std::cout<<"Rank: "<<rank<<"executing local Sort for size: "<<pNodes.size()<<std::endl;
            #ifdef PROFILE_TREE_SORT
                //MPI_Barrier(pcomm);
                t2=std::chrono::high_resolution_clock::now();//MPI_Wtime();
            #endif
            SFC::seqSort::SFC_treeSort(&(*(pNodes.begin())),pNodes.size(),pOutSorted,pOutConstruct,pOutBalanced,pMaxDepth,pMaxDepth,parent,rot_id,k,options);

            #ifdef PROFILE_TREE_SORT
                localSort_time=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t2).count();//MPI_Wtime()-t2;
                //MPI_Barrier(pcomm);
            #endif
	    
            /*assert(seq::test::isUniqueAndSorted(pOutSorted));
            assert(seq::test::isUniqueAndSorted(pOutConstruct));
            assert(seq::test::isUniqueAndSorted(pOutBalanced));*/
	        // freee the comm from the last stage. 
	        MPI_Comm_free(&comm);

            if(options & TS_REMOVE_DUPLICATES) {

                //if(!rank) std::cout<<"Executing  par::RD begin"<<std::endl;

                #ifdef PROFILE_TREE_SORT
                    //MPI_Barrier(pcomm);
                    t2=std::chrono::high_resolution_clock::now();//MPI_Wtime();
                #endif

                int new_rank, new_size;
                MPI_Comm new_comm;
                // very quick and dirty solution -- assert that tmpVec is non-emply at every processor (repetetive calls to splitComm2way exhaust MPI resources)
                par::splitComm2way(pOutSorted.empty(), &new_comm, pcomm);
                //new_comm = pcomm;
                //assert(!pNodes.empty());
                MPI_Comm_rank(new_comm, &new_rank);
                MPI_Comm_size(new_comm, &new_size);


                #ifdef __DEBUG_PAR__
                MPI_Barrier(comm);
                if(!rank) {
                    std::cout<<"RemDup: Stage-4 passed."<<std::endl;
                }
                MPI_Barrier(comm);
                #endif

                //Checking boundaries...
                if (!pOutSorted.empty()) {

                    T begin =pOutSorted[0];
                    T end = pOutSorted[pOutSorted.size() - 1];
                    T endRecv;
                    T beginRecv;

                    //communicate end to the next processor.
                    MPI_Status status;

                    par::Mpi_Sendrecv<T, T>(&end, 1, ((new_rank < (new_size - 1)) ? (new_rank + 1) : 0), 1, &endRecv,
                                            1, ((new_rank > 0) ? (new_rank - 1) : (new_size - 1)), 1, new_comm,
                                            &status);



                    //Remove endRecv if it exists (There can be no more than one copy of this)
                    if (new_rank) {
                        typename std::vector<T>::iterator Iter = std::find(pOutSorted.begin(), pOutSorted.end(), endRecv);
                        if (Iter != pOutSorted.end()) {
                            pOutSorted.erase(Iter);
                        }//end if found



                    }//end if p not 0


                    bool state=true;
		            bool state_global=true;
                    unsigned int count=0;
		            //@milindasf : Possible location for MPI hang, if for a one processor if above is not true, then followign sendrecv will get hanged. 

                    while(count<pOutSorted.size() & state_global ) {

                       begin=pOutSorted[count];
                       end=pOutSorted.back();
                       state=false;
                       par::Mpi_Sendrecv<T, T>(&begin, 1, ((new_rank > 0) ? (new_rank - 1) : (new_size - 1)), 1,
                                                &beginRecv, 1, ((new_rank < (new_size - 1)) ? (new_rank + 1) : 0), 1, new_comm,
                                                &status);

                        /*if(beginRecv.isAncestor(end))
                        {
                           std::cout<<"for rank: "<<new_rank<<" beginRecv: "<<beginRecv<<" is ancestor to : "<<end<<std::endl;
                        }*/

                       while (end.isAncestor(beginRecv)) {
                            state = true;
                            pOutSorted.pop_back();
                            end = pOutSorted.back();
                       }
                       count++;
		               MPI_Allreduce(&state,&state_global,1,MPI_CXX_BOOL,MPI_LOR,new_comm);


                    }

                }//end if not empty
                
                // free the new comm. 
                MPI_Comm_free(&new_comm);

                //if(!rank) std::cout<<"Executing  par::RD end"<<std::endl;
                #ifdef PROFILE_TREE_SORT
                    remove_duplicates_par=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t2).count();//MPI_Wtime()-t2;
                //MPI_Barrier(pcomm);
                #endif

            }

            #ifdef PROFILE_TREE_SORT
                total_rd=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t4).count();


                int rank_g,npes_g;
                MPI_Comm_rank(pcomm,&rank_g);
                MPI_Comm_size(pcomm,&npes_g);

                par::Mpi_Reduce(&splitter_fix_all,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&splitter_fix_all,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&splitter_fix_all,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g) {
                    stat_property[1] = stat_property[1] / (double) npes_g;
                    //std::cout<<"Rank: "<<rank_g<<"splitter_fix_all: "<<stat_property[0]<< ": "<<stat_property[1]<<" : "<<stat_property[2]<<std::endl;
                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);
                }


                par::Mpi_Reduce(&splitter_time,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&splitter_time,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&splitter_time,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);

                }


                par::Mpi_Reduce(&all2all1_time,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&all2all1_time,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&all2all1_time,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);

                }


                par::Mpi_Reduce(&all2all2_time,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&all2all2_time,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&all2all2_time,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);

                }

                par::Mpi_Reduce(&localSort_time,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&localSort_time,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&localSort_time,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);

                }


                par::Mpi_Reduce(&remove_duplicates_seq,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&remove_duplicates_seq,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&remove_duplicates_seq,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);

                }


                par::Mpi_Reduce(&auxBalOCt_time,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&auxBalOCt_time,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&auxBalOCt_time,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);

                }


                par::Mpi_Reduce(&remove_duplicates_par,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&remove_duplicates_par,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&remove_duplicates_par,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);


                }



                par::Mpi_Reduce(&total_rd,&stat_property[0],1,MPI_MIN,0,pcomm);
                par::Mpi_Reduce(&total_rd,&stat_property[1],1,MPI_SUM,0,pcomm);
                par::Mpi_Reduce(&total_rd,&stat_property[2],1,MPI_MAX,0,pcomm);

                if(!rank_g)
                {
                    stat_property[1] = stat_property[1] / (double) npes_g;

                    stats.push_back(stat_property[0]);
                    stats.push_back(stat_property[1]);
                    stats.push_back(stat_property[2]);


                }

                for(int i=0;i<SF_Stages;i++)
                {
                    par::Mpi_Reduce(&sf_full[i],&stat_property[0],1,MPI_MIN,0,pcomm);
                    par::Mpi_Reduce(&sf_full[i],&stat_property[1],1,MPI_SUM,0,pcomm);
                    par::Mpi_Reduce(&sf_full[i],&stat_property[2],1,MPI_MAX,0,pcomm);

                    stat_property[1] = stat_property[1] / (double) npes_g;

                    if(!rank_g)
                    {
                        stats.push_back(stat_property[0]);
                        stats.push_back(stat_property[1]);
                        stats.push_back(stat_property[2]);
                    }


                    par::Mpi_Reduce(&sf_all2all[i],&stat_property[0],1,MPI_MIN,0,pcomm);
                    par::Mpi_Reduce(&sf_all2all[i],&stat_property[1],1,MPI_SUM,0,pcomm);
                    par::Mpi_Reduce(&sf_all2all[i],&stat_property[2],1,MPI_MAX,0,pcomm);

                    stat_property[1] = stat_property[1] / (double) npes_g;

                    if(!rank_g)
                    {
                        stats.push_back(stat_property[0]);
                        stats.push_back(stat_property[1]);
                        stats.push_back(stat_property[2]);
                    }


                    par::Mpi_Reduce(&sf_splitters[i],&stat_property[0],1,MPI_MIN,0,pcomm);
                    par::Mpi_Reduce(&sf_splitters[i],&stat_property[1],1,MPI_SUM,0,pcomm);
                    par::Mpi_Reduce(&sf_splitters[i],&stat_property[2],1,MPI_MAX,0,pcomm);

                    stat_property[1] = stat_property[1] / (double) npes_g;

                    if(!rank_g)
                    {
                        stats.push_back(stat_property[0]);
                        stats.push_back(stat_property[1]);
                        stats.push_back(stat_property[2]);
                    }

                }



            #endif



            if((options & TS_CONSTRUCT_OCTREE) | (options & TS_BALANCE_OCTREE))
            {


                MPI_Comm_rank(pcomm,&rank);
                MPI_Comm_size(pcomm,&npes);

                /*std::cout<<"rank: "<<rank<<" Remove Duplicates"<<std::endl;*/
                #ifdef PROFILE_TREE_SORT
                    stats_previous.clear();
                    stats_previous.insert(stats_previous.end(),stats.begin(),stats.end());
                #endif
                std::vector<T> tmpSorted;
                std::vector<T> tmpConstructed;
                std::vector<T> tmpBalanced;
                T root(0,0,0,0,m_uiDim,pMaxDepth);
                if(options & TS_CONSTRUCT_OCTREE) {


                #ifdef PROFILE_TREE_SORT
                    //MPI_Barrier(pcomm);
                    t3=std::chrono::high_resolution_clock::now();//MPI_Wtime();
                #endif
                    par::partitionW<T>(pOutConstruct,NULL,pcomm);
                    SFC::parSort::SFC_treeSort(pOutConstruct,tmpSorted,tmpConstructed,tmpBalanced,loadFlexibility,pMaxDepth,root,rot_id,k,1,sf_k,pcomm);

                #ifdef PROFILE_TREE_SORT
                    construction_time=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t3).count();
                            //MPI_Wtime()-t3;

                    //MPI_Barrier(pcomm);

                    int rank_g,npes_g;
                    MPI_Comm_rank(pcomm,&rank_g);
                    MPI_Comm_size(pcomm,&npes_g);


                    par::Mpi_Reduce(&construction_time,&stat_property[0],1,MPI_MIN,0,pcomm);
                    par::Mpi_Reduce(&construction_time,&stat_property[1],1,MPI_SUM,0,pcomm);
                    par::Mpi_Reduce(&construction_time,&stat_property[2],1,MPI_MAX,0,pcomm);

                    if(!rank_g)
                    {
                        stat_property[1] = stat_property[1] / (double) npes_g;

                        stats.push_back(stat_property[0]);
                        stats.push_back(stat_property[1]);
                        stats.push_back(stat_property[2]);

                    }
                #endif

                    std::swap(tmpSorted,pOutConstruct);
                    tmpSorted.clear();

                    //if(!rank) std::cout<<"Executing  Construct end "<<std::endl;

                }
                if(options & TS_BALANCE_OCTREE) {
                    //if(!rank) std::cout<<"Executing balancing"<<std::endl;
                    #ifdef PROFILE_TREE_SORT
                        //MPI_Barrier(pcomm);
                        t3=std::chrono::high_resolution_clock::now();//MPI_Wtime();
                    #endif
                    par::partitionW<T>(pOutBalanced,NULL,pcomm);
                    SFC::parSort::SFC_treeSort(pOutBalanced,tmpSorted,tmpConstructed,tmpBalanced,loadFlexibility,pMaxDepth,root,rot_id,k,1,sf_k,pcomm);

                    #ifdef PROFILE_TREE_SORT
                        balancing_time=std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t3).count();//MPI_Wtime()-t3;
                        int rank_g,npes_g;
                        MPI_Comm_rank(pcomm,&rank_g);
                        MPI_Comm_size(pcomm,&npes_g);

                        par::Mpi_Reduce(&balancing_time,&stat_property[0],1,MPI_MIN,0,pcomm);
                        par::Mpi_Reduce(&balancing_time,&stat_property[1],1,MPI_SUM,0,pcomm);
                        par::Mpi_Reduce(&balancing_time,&stat_property[2],1,MPI_MAX,0,pcomm);

                        if(!rank_g)
                        {
                            stat_property[1] = stat_property[1] / (double) npes_g;

                            stats.push_back(stat_property[0]);
                            stats.push_back(stat_property[1]);
                            stats.push_back(stat_property[2]);

                        }

                        //MPI_Barrier(pcomm);

                    #endif
                    std::swap(tmpSorted,pOutBalanced);
                    tmpSorted.clear();


                }

            }



        }





       /* template <typename T>
        void SFC_PartitionW(std::vector<T>&pNodes,double loadFlexibility, unsigned int maxDepth,MPI_Comm comm)
        {
            T root=T(m_uiDim,maxDepth);
            std::vector<T> tmp;
            SFC_treeSort(pNodes,tmp,tmp,tmp,loadFlexibility,maxDepth,root,ROOT_ROTATION,1,0,NUM_NPES_THRESHOLD,comm);
            tmp.clear();
        }*/





        //========================================================= Function definition end.=========================================================================================



    }// end of namespace parSort


}// end of namespace SFC





#endif //SFCSORTBENCH_SFCSORT_H
