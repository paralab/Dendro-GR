//
// Created by milinda on 9/25/16.
//

/**
 * @author: Milinda Fernando
 * School of Computing, University of Utah
 * @breif Contains SFC based searching functionality (finalized) for Morton and Hilbert Curve.
 * this searching functionality is mainly used the mesh generation stage.
 *
 * */


#ifndef SFCSORTBENCH_SFCSEARCH_H
#define SFCSORTBENCH_SFCSEARCH_H

#include "sfcSort.h"
#include "dendro.h"
#include "TreeNode.h"


namespace  SFC
{
    namespace seqSearch
    {
        /** @author Milinda Fernando
         *  @breif Performs sequential tree search oof pKeys in the octree pNodes.
         *  Assumptions: pNodes can be sorted or unsorted, but pNodes must be unique. If some keys were not found they need to be sent to the correst processor using par::treeSort.
         *  Not found keys are marked
         *
         * */
        template <typename TKey, typename TOctant>
        void SFC_treeSearch(TKey* pKeys , TOctant* pNodes,std::vector<unsigned int>& pSearchIndex,DendroIntL nKeyBegin,DendroIntL nKeyEnd,DendroIntL nNodeBegin,DendroIntL nNodeEnd,unsigned int pMaxDepthBit,unsigned int pMaxDepth, unsigned int rot_id);


        template <typename TKey, typename TOctant>
        void SFC_treeSearch(TKey* pKeys , TOctant* pNodes,DendroIntL nKeyBegin,DendroIntL nKeyEnd,DendroIntL nNodeBegin,DendroIntL nNodeEnd,unsigned int pMaxDepthBit,unsigned int pMaxDepth, unsigned int rot_id);




        /** @author Milinda Fernando
         *  @breif Performs sequential tree search oof pKeys in the octree pNodes.
         *  Assumptions: pNodes can be sorted or unsorted, but pNodes must be unique. If some keys were not found they need to be sent to the correst processor using par::treeSort.
         *  Keys which are not found will be flags in the m_uiLevel. (max level is 31 and rest of the bits are used to indicate the flags for mesh generation part. )
         *
         * */
        template <typename TKey, typename TOctant>
        void SFC_treeSearch(TKey* pKeys , TOctant* pNodes,std::vector<unsigned int> & pSearchIndex,DendroIntL nKeyBegin,DendroIntL nKeyEnd,DendroIntL nNodeBegin,DendroIntL nNodeEnd,unsigned int pMaxDepthBit,unsigned int pMaxDepth, unsigned int rot_id)
        {

            unsigned int lev=pMaxDepth-pMaxDepthBit;
            DendroIntL splitterKeys[(NUM_CHILDREN + 1)];
            DendroIntL splitterNodes[(NUM_CHILDREN + 1)];

            unsigned int hindex = 0;
            unsigned int hindexN = 0;
            unsigned int index = 0;

            hindex = (rotations[2 * NUM_CHILDREN * rot_id] - '0');
            hindexN = (rotations[2 * NUM_CHILDREN * rot_id +1] - '0');

            SFC::seqSort::SFC_bucketing(pKeys,lev,pMaxDepth,rot_id,nKeyBegin,nKeyEnd,splitterKeys);
            SFC::seqSort::SFC_bucketing(pNodes,lev,pMaxDepth,rot_id,nNodeBegin,nNodeEnd,splitterNodes);


            if((pMaxDepthBit) && (((nNodeEnd-nNodeBegin)>=1) && (pNodes[splitterNodes[hindex]].getLevel()>lev))  ) {

                    /*for (unsigned int k = splitterKeys[hindex]; k < splitterKeys[hindexN]; k++) {
                        if ((pKeys[k].getLevel() == lev) & (pNodes[splitterNodes[hindex]].getLevel() == lev)) {
                            pKeys[k].setFlag((pKeys[k].getFlag() | OCT_FOUND));
                            assert(pKeys[k].getFlag() & OCT_FOUND);
                            assert(pKeys[k] == pNodes[splitterNodes[hindex]]);
                            pSearchIndex.push_back(splitterNodes[hindex]);
                            std::cout<<"maxDepth Reached"<<std::endl;
                        }
                    }*/

                for (int i = 0; i < NUM_CHILDREN; i++) {
                    hindex = (rotations[2 * NUM_CHILDREN * rot_id + i] - '0');
                    if (i == (NUM_CHILDREN - 1))
                        hindexN = i + 1;
                    else
                        hindexN = (rotations[2 * NUM_CHILDREN * rot_id + i + 1] - '0');

                    assert(splitterKeys[hindex] <= splitterKeys[hindexN]);
                    assert(splitterNodes[hindex] <= splitterNodes[hindexN]);

                    index = HILBERT_TABLE[NUM_CHILDREN * rot_id + hindex];
                    SFC_treeSearch(pKeys,pNodes,pSearchIndex,splitterKeys[hindex],splitterKeys[hindexN],splitterNodes[hindex],splitterNodes[hindexN],(pMaxDepthBit-1),pMaxDepth,index);

                }

            }else
            {
                if(((nNodeEnd-nNodeBegin)==1) )
                {
                    assert(pNodes[nNodeBegin].getLevel()==lev);
                    for(unsigned int k=nKeyBegin;k<nKeyEnd;k++)
                    {
                        assert(((pNodes[nNodeBegin].isAncestor(pKeys[k]))||(pNodes[nNodeBegin]==pKeys[k])));
                        pKeys[k].setFlag((pKeys[k].getFlag() | OCT_FOUND));
                        assert(pKeys[k].getFlag() & OCT_FOUND);
                        pSearchIndex.push_back(nNodeBegin);
                    }

                }
            }

        }



        template <typename TKey, typename TOctant>
        void SFC_treeSearch(TKey* pKeys , TOctant* pNodes,DendroIntL nKeyBegin,DendroIntL nKeyEnd,DendroIntL nNodeBegin,DendroIntL nNodeEnd,unsigned int pMaxDepthBit,unsigned int pMaxDepth, unsigned int rot_id)
        {

            unsigned int lev=pMaxDepth-pMaxDepthBit;
            if(nKeyEnd==nKeyBegin) return;
            //std::cout<<"call: "<<pMaxDepthBit<<": NBegin: "<<nNodeBegin<<" NEnd: "<<nNodeEnd<<" KBegin: "<<nKeyBegin<<" KEnd: "<<nKeyEnd<<" NodeLev: "<<pNodes[nNodeBegin].getLevel()<<" lev: "<<lev<< std::endl;

            //if( pMaxDepthBit && ((nNodeEnd-nNodeBegin)>=1) && (pNodes[nNodeBegin].getLevel()>lev))
            if( pMaxDepthBit && ((nNodeEnd-nNodeBegin)>1) && (pNodes[nNodeBegin].getLevel()>lev))
            {

                unsigned int hindex = 0;
                unsigned int hindexN = 0;
                unsigned int index = 0;
                DendroIntL splitterKeys[(NUM_CHILDREN + 1)];
                DendroIntL splitterNodes[(NUM_CHILDREN + 1)];

                SFC::seqSort::SFC_bucketing(pKeys, lev, pMaxDepth, rot_id, nKeyBegin, nKeyEnd, splitterKeys);
                SFC::seqSort::SFC_bucketing(pNodes, lev, pMaxDepth, rot_id, nNodeBegin, nNodeEnd, splitterNodes);


                for (int i = 0; i < NUM_CHILDREN; i++) {
                    hindex = (rotations[2 * NUM_CHILDREN * rot_id + i] - '0');
                    if (i == (NUM_CHILDREN - 1))
                        hindexN = i + 1;
                    else
                        hindexN = (rotations[2 * NUM_CHILDREN * rot_id + i + 1] - '0');

                    assert(splitterKeys[hindex] <= splitterKeys[hindexN]);
                    assert(splitterNodes[hindex] <= splitterNodes[hindexN]);

                    index = HILBERT_TABLE[NUM_CHILDREN * rot_id + hindex];
                    if(splitterNodes[hindex]!=splitterNodes[hindexN]) SFC_treeSearch(pKeys,pNodes,splitterKeys[hindex],splitterKeys[hindexN],splitterNodes[hindex],splitterNodes[hindexN],(pMaxDepthBit-1),pMaxDepth,index);

                }

            }else
            {
                if(((nNodeEnd-nNodeBegin)==1) )
                {
                    assert(pNodes[nNodeBegin].getLevel()>=lev);
                    for(unsigned int k=nKeyBegin;k<nKeyEnd;k++)
                    {
                        //assert(((pNodes[nNodeBegin].isAncestor(pKeys[k]))||(pNodes[nNodeBegin]==pKeys[k]))); // Note: Since we are sarching for the ghost elements we might not find a given key. Hence this should be disabled
                        if((pNodes[nNodeBegin].isAncestor(pKeys[k])) || (pNodes[nNodeBegin]==pKeys[k])){
                            pKeys[k].setFlag((pKeys[k].getFlag() | OCT_FOUND));
                            assert(pKeys[k].getFlag() &  OCT_FOUND);  // Note: Since we are sarching for the ghost elements we might not find a given key. Hence this should be disabled
                            pKeys[k].setSearchResult(nNodeBegin);
                        }


                    }
                    //std::cout<<"key update: "<<nKeyBegin<<" to : "<<nKeyEnd<<std::endl;

                }else if(((nNodeEnd-nNodeBegin)>1))
                {

                    //std::cout<<"Else call: "<<pMaxDepthBit<<": NBegin: "<<nNodeBegin<<" NEnd: "<<nNodeEnd<<" KBegin: "<<nKeyBegin<<" KEnd: "<<nKeyEnd<<" NodeLev: "<<pNodes[nNodeBegin].getLevel()<<" lev: "<<lev<< std::endl;
                    unsigned int nCnt=nNodeBegin;
                    for(unsigned int k=nKeyBegin;k<nKeyEnd;k++)
                    {
                        while( (nCnt<nNodeEnd) && (!(pNodes[nCnt].isAncestor(pKeys[k]))&& !(pNodes[nCnt]==pKeys[k])))
                            nCnt++;

                        //assert(((pNodes[nNodeBegin].isAncestor(pKeys[k]))||(pNodes[nNodeBegin]==pKeys[k])));  // Note: Since we are sarching for the ghost elements we might not find a given key. Hence this should be disabled
                        if(nCnt<nNodeEnd) {
                            pKeys[k].setFlag((pKeys[k].getFlag() | OCT_FOUND));
                            //assert(pKeys[k].getFlag() & OCT_FOUND);  // Note: Since we are sarching for the ghost elements we might not find a given key. Hence this should be disabled
                            assert(((pNodes[nCnt].isAncestor(pKeys[k]))||(pNodes[nCnt]==pKeys[k])));
                            pKeys[k].setSearchResult(nCnt);
                            nCnt=nNodeBegin;
                        }

                    }


                }
            }





        }





    } // namespace seqSearch


} // namespace SFC end



#endif //SFCSORTBENCH_SFCSEARCH_H
