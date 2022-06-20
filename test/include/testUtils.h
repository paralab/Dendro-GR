
/**
  @file testUtils.h
  @brief  A Set of utilities to test octrees.	
  @author	Rahul S. Sampath, rahul.sampath@gmail.com
  @author Milinda Shayamal Fernando.
  School of COmputing, University of Utah
  milinda@cs.utah.edu
  */ 

#ifndef _TESTUTILS_H_
#define _TESTUTILS_H_

#include <mpi.h>
#include <vector>
#include "TreeNode.h"
#include "treenode2vtk.h"
#include <climits>
#include <math.h>
#include <block.h>
#include <functional>


#ifdef DIM_2
#define NUM_CHILDREN 4
    #define ROTATION_OFFSET 8
#else
#define NUM_CHILDREN 8
#define ROTATION_OFFSET 16
#endif


namespace seq {

  /**
    @namespace test
    @author Rahul Sampath
    @brief A collection of functions for debugging
    */
  namespace test {

      /**
        @fn
        @param nodes[in] The vector of nodes that have to be tested.
        @return true if it is sorted and false otherwise
      **/ 
      template<typename T>
      bool isSorted(const std::vector<T >& nodes);

      template<typename T>
      bool isSorted_all_pairs(const std::vector<T > & nodes);

      template<typename T>
      bool isSorted(T * nodes, unsigned int sz);

      /**@brief checks whether the numeric array constains nan values. */
      template<typename T>
      bool isNAN(const T * in, const unsigned int sz);

      /**@brief checks whether the numeric array constains nan values. for block local region */
      template<typename T>
      bool isBlockNAN(const T * in, const unsigned int* sz);

      /**@brief checks whether the numeric array constains nan values. for block local region */
      template<typename T>
      bool isBlockNAN(const T * in, const unsigned int* sz,unsigned int flag);


      /**
        @fn
        @param nodes[in] The vector of nodes that have to be tested.
        @return true if it is sorted and unique and false otherwise
      **/ 
      template<typename T>
      bool isUniqueAndSorted(const std::vector<T >& nodes);


      template <typename T>
      bool containsAncestor(const std::vector<T > & nodes);


      template <typename T>
      bool  isComplete (const std::vector<T>& nodes);

      /** @author Milinda Fernando
       * @breif Check whether the local E2E mapping is correct.
       * */
      template <typename T>
      bool checkE2EMapping( const std::vector<unsigned int >&  E2EMap,  const std::vector<T>& allNodes, unsigned int localBegin, unsigned int localEnd ,unsigned int k_s,unsigned numDirections);


      /**
       * @author: Milinda Fernando.
       * @brief: Check whether the E2N mapping is correct based on the DG indexing.
       * This assumes that StensilSz is 1
       */
      template <typename T>
      bool checkE2NMapping( const std::vector<unsigned int >&  E2EMap , const std::vector<unsigned int >& E2NMap, const std::vector<T>& allNodes,unsigned int numDirections,unsigned int elementOrder);


  }//end namespace
}//end namespace

namespace par {

  /**
    @namespace test
    @author Rahul Sampath
    @brief A collection of functions for debugging
    */
  namespace test {

      template<typename T>
      bool isSorted(const std::vector<T>& nodes, MPI_Comm comm);

      template<typename T>
      bool isUniqueAndSorted(const std::vector<T >& nodes,MPI_Comm comm) ;


      template <typename T>
      bool containsAncestor(const std::vector<T > &nodes, MPI_Comm comm);

      //@author: Milinda Fernando. Assess the completeness of the octree based on the volumes of the octants.
      template <typename T>
      bool  isComplete (const std::vector<T>& nodes, MPI_Comm comm);


  }//end namespace
}//end namespace

namespace ot {

  //class TreeNode;

  namespace test {
      /**
      @fn
      @param nodes[in] The vector of nodes that have to be tested.
      @return false if any ancestor of any element in the vector is also present in the vector and true otherwise.
      @remark The function uses Binary Search and hence expects the vector of nodes to be sorted.
    */

    template<typename T>
    bool isLinear(const std::vector<T >& nodes) ;

    template<typename T>
    bool isBalanced(unsigned int dim, unsigned int maxDepth, char* failFileName,
        const std::vector<T>& nodes, bool incCorn, unsigned int maxLevDiff) ;

    template<typename T>
    bool isBalancedInternal(unsigned int dim, unsigned int maxDepth,
        char*failFileName,	const std::vector<T> & nodes,
        TreeNode holder, bool incCorn, unsigned int maxLevDiff) ;


    /**
     * @author Milinda Fernanado
     * @brief for a given sorted 2:1 balanced octree and a block list this procedure checks that for each block, number of octants at level l (reg grid level) cross all blocks
     * should equal to the level l octants across the octree.
     * This ensures that blocks have captured all the octants corresponding to regular grid level l .
     *
     * @param[in] pNodes 2:1 balance sorted octree
     * @param[in] blockList set of computed block list.
     * @param[in] d_min min depth of the octree
     * @param[in] d_max max depth of the octree
     *
     * */
    template<typename T,typename B>
    bool isBlockListValid(const std::vector<T>& pNodes, std::vector<T> & blockList,unsigned int d_min, unsigned int d_max,unsigned int nBegin,unsigned int nEnd);



 




  }//end namespace
}//end namespace










#include "testUtils.tcc"

#endif

