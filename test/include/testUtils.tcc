
/**
  @file testUtils.txx
  @author Rahul S. Sampath, rahul.sampath@gmail.com

  @author: Milinda Shayamal Fernando. milinda@cs.utah.edu
  School of Computing
  University of Utah

  */

#include "dtypes.h"
#include "parUtils.h"
#include "seqUtils.h"
#include <iostream>
#include <stdint.h>

namespace seq { 
  namespace test {

    template <typename T>
    bool isSorted(const std::vector<T > & nodes) {
        for (unsigned int i = 1; i < nodes.size(); i++) {
          if ( nodes[i] < nodes[i-1] ) {
            std::cout<<"\n Local Sort Check failed for: "<<nodes[i]<<" and "
              <<nodes[i-1]<<std::endl<<std::endl;
            return false;
          }
        }
        return true;
      }

    template <typename T>
    bool isSorted_all_pairs(const std::vector<T > & nodes) {
        for (unsigned int i = 0; i < nodes.size(); i++) {
          for(unsigned int j=i+1;j<nodes.size();j++)
          {
            if(nodes[i]>nodes[j]) {
              std::cout << "\n Local Sort_pair Check failed for: " << nodes[i] << " and " << nodes[j] << std::endl << std::endl;
              return false;
            }
          }

        }
        return true;
      }


    template <typename T>
    bool isSorted(T* nodes, unsigned int sz) {
        for (unsigned int i = 1; i < sz; i++) {
          if ( nodes[i] < nodes[i-1] ) {
            std::cout<<"\n Local Sort Check failed for: "<<nodes[i]<<" and "
              <<nodes[i-1]<<std::endl<<std::endl;
            std::cout<<"nodes[i] < nodes[i-1]:"<<(nodes[i] < nodes[i-1])<<std::endl;
            std::cout<<"nodes[i] > nodes[i-1]:"<<(nodes[i] > nodes[i-1])<<std::endl;
            return false;
          }
        }
        return true;
      }


    template<typename T>
    bool isNAN(const T * in, const unsigned int sz)
      {
          for(unsigned int i=0;i<sz;i++)
              if(std::isnan(in[i]))
                  return true;

          return false;

      }


    template<typename T>
    bool isBlockNAN(const T * in,const unsigned int* sz)
      {
          for(unsigned int k=3;k<sz[2]-3;k++)
              for(unsigned int j=3;j<sz[1]-3;j++)
                 for(unsigned int i=3;i<sz[0]-3;i++)
                     if(std::isnan(in[sz[0]*sz[1]*k+j*sz[0]+i]))
                         return true;

          return false;


      }


    template<typename T>
    bool isBlockNAN(const T * in, const unsigned int* sz,unsigned int flag)
      {

          if(flag==0)
          { // implies this is an internal block.
              for(unsigned int k=0;k<sz[2];k++)
                  for(unsigned int j=0;j<sz[1];j++)
                      for(unsigned int i=0;i<sz[0];i++)
                          if(std::isnan(in[sz[0]*sz[1]*k+j*sz[0]+i]))
                              return true;

              return false;

          }else
          {
              return isBlockNAN(in,sz);
          }

      }



    template <typename T>
    bool isUniqueAndSorted(const std::vector<T > & nodes) {
        for (unsigned int i = 1; i < nodes.size(); i++) {
          if ( nodes[i] <= nodes[i-1] ) {
            std::cout<<"\n Local Sort+Unique Check failed for: "<<nodes[i]<<" and "
              <<nodes[i-1]<<std::endl;
            std::cout<<"nodes[i] <= nodes[i-1]:"<<(nodes[i] <= nodes[i-1])<<std::endl;
            std::cout<<"nodes[i] >= nodes[i-1]:"<<(nodes[i] >= nodes[i-1])<<std::endl;
            return false;
          }
        }
        return true;
      }


    template <typename T>
    bool containsAncestor(const std::vector<T > & nodes) {
        for (unsigned int i = 1; i < nodes.size(); i++) {
          if (nodes[i-1].isAncestor(nodes[i]) ) {
            std::cout<<"\n Local contains ancestors failed for: "<<nodes[i]<<" and "
            <<nodes[i-1]<<std::endl;
            std::cout<<"nodes[i] <= nodes[i-1]:"<<(nodes[i] <= nodes[i-1])<<std::endl;
            std::cout<<"nodes[i] >= nodes[i-1]:"<<(nodes[i] >= nodes[i-1])<<std::endl;
            return false;
          }
        }
        return true;
      }



    template <typename T>
    bool  isComplete (const std::vector<T>& nodes) {

    //std::cout<<"size of DendroIntL: "<<sizeof(DendroIntL)<<std::endl;
    // Note: Current implementation of isComplete works only with maxDepth <21 due to capacity of DendroIntL.
    DendroIntL maxDepth = 0;
    DendroIntL volMaxDepth=0;
    DendroIntL dim =3;
    DendroIntL c=1;
    if (nodes.size() != 0) {
      maxDepth = nodes[0].getMaxDepth();
      dim = nodes[0].getDim();
      volMaxDepth=(c<<(dim*maxDepth));
      DendroIntL vol = 0;
      DendroIntL vol1 = 0;

      for (DendroIntL i = 0; i < nodes.size(); i++) {
        vol1=(dim*(maxDepth - nodes[i].getLevel()));
        vol1 =c<<vol1;//(1ull << ((DendroIntL)(dim * (maxDepth - nodes[i].getLevel()))));
        if (i < (nodes.size() - 1) && (nodes[i].isAncestor(nodes[i + 1]))) {
          std::cout << i << " contains duplicate nodes. (children)" << std::endl;
          return false;
        }

        vol = vol + vol1;
        //if(vol1<0)
          //std::cout<<"vol: "<<vol<<" vol1: "<<vol1 <<" len:  "<<(maxDepth - nodes[i].getLevel())<<" node : "<<nodes[i]<<std::endl;
      }

      bool state=(vol==volMaxDepth);

      if(!state)
      {
        std::cout<<"Volume based on MaxDepth: "<<volMaxDepth<<std::endl;
        std::cout<<"Computed Volume based on OCtants: "<<vol<<std::endl;
        std::cout<<"vol diff: "<<abs(volMaxDepth-vol)<<std::endl;
      }

      return state;

    } else {
      return true;
    }

  }


    template <typename T>
    bool checkE2EMapping( const std::vector<unsigned int >&  E2EMap,  const std::vector<T>& allNodes, unsigned int localBegin, unsigned int localEnd ,unsigned int k_s,unsigned numDirections) {
          unsigned int domain_max = (1u<<(m_uiMaxDepth-1))+1;
          const ot::TreeNode *inPtr = (&(*(allNodes.begin())));

          ot::TreeNode tmp;
          ot::TreeNode lookUp;

          std::vector<ot::TreeNode> keys;
          std::vector<ot::TreeNode> result;

          int rank;
          unsigned int lookUPIndex;
          MPI_Comm_rank(MPI_COMM_WORLD,&rank);
          //unsigned int numDirections=(1u << m_uiDim) - 2 * (m_uiDim - 2);

          for (unsigned int i = localBegin; i < localEnd; i++) {
              unsigned int myLev = inPtr[i].getLevel();

              if (myLev == 1) continue;

              unsigned int mySz = (1u << (m_uiMaxDepth - myLev));
              unsigned int myX = inPtr[i].getX();
              unsigned int myY = inPtr[i].getY();
              unsigned int myZ = inPtr[i].getZ();



              for (unsigned int k = 1; k <= k_s; k++) {

                  if ((myX + k * mySz) < domain_max) {
                      tmp = ot::TreeNode((myX + k * mySz), myY, myZ, (OCT_KEY_RIGHT | m_uiMaxDepth), m_uiDim, m_uiMaxDepth);
                      lookUPIndex=E2EMap[i*k_s*numDirections+ (OCT_DIR_RIGHT * k_s + k - 1)];
                      //assert(lookUPIndex!=UINT_MAX);
                      if(lookUPIndex!=UINT_MAX) {
                          assert(lookUPIndex<allNodes.size());
                          lookUp = allNodes[lookUPIndex];
                          if (!(tmp.isAncestor(lookUp)) && (tmp!=lookUp) && (!lookUp.isAncestor(tmp))) {

                              std::cout << "rank: " << rank << "owner: "<<inPtr[i]<<" key right: " << tmp << " lookup: " << lookUp<< " lookupIndex: "<<lookUPIndex << std::endl;
                              return false;
                          }
                      }
                  }


                  if (myX >0) {
                      tmp = ot::TreeNode((myX - 1), myY, myZ, (OCT_KEY_LEFT | m_uiMaxDepth), m_uiDim, m_uiMaxDepth);
                      lookUPIndex=E2EMap[i*k_s*numDirections+OCT_DIR_LEFT * k_s + k - 1];
                      //assert(lookUPIndex!=UINT_MAX);
                      if(lookUPIndex!=UINT_MAX) {
                          assert(lookUPIndex<allNodes.size());
                          lookUp = allNodes[lookUPIndex];
                          if (!(tmp.isAncestor(lookUp)) && (tmp!=lookUp) && (!lookUp.isAncestor(tmp))) {
                              std::cout << "rank: " << rank << "owner: "<<inPtr[i]<<" key left: " << tmp << " lookup: " << lookUp<< " lookupIndex: "<<lookUPIndex << std::endl;
                              return false;
                          }
                      }
                  }
                  if ((myY + k * mySz) < domain_max) {
                      tmp = ot::TreeNode(myX, (myY+ k * mySz), myZ, (OCT_KEY_UP | m_uiMaxDepth), m_uiDim, m_uiMaxDepth);
                      lookUPIndex=E2EMap[i*k_s*numDirections+OCT_DIR_UP * k_s + k - 1];
                      //assert(lookUPIndex!=UINT_MAX);
                      if(lookUPIndex!=UINT_MAX) {
                          assert(lookUPIndex<allNodes.size());
                          lookUp = allNodes[lookUPIndex];
                          if (!(tmp.isAncestor(lookUp)) && (tmp!=lookUp) && (!lookUp.isAncestor(tmp))) {
                              std::cout << "rank: " << rank << " key up: " << tmp << " lookup: " << lookUp << std::endl;
                              return false;
                          }
                      }

                  }

                  if (myY >0) {
                      tmp = ot::TreeNode(myX, (myY-1), myZ, (OCT_KEY_DOWN | m_uiMaxDepth), m_uiDim, m_uiMaxDepth);
                      lookUPIndex=E2EMap[i*k_s*numDirections+OCT_DIR_DOWN * k_s + k - 1];
                      //assert(lookUPIndex!=UINT_MAX);
                      if(lookUPIndex!=UINT_MAX) {
                          assert(lookUPIndex<allNodes.size());
                          lookUp = allNodes[lookUPIndex];
                          if (!(tmp.isAncestor(lookUp)) && (tmp!=lookUp) && (!lookUp.isAncestor(tmp))) {
                              if(!rank)std::cout << "rank: " << rank << "element ID: "<<(i-localBegin)<<" element : "<<inPtr[i]<<" key down: " << tmp << " lookup: " << lookUp <<" lookupIndex: "<<lookUPIndex << std::endl;
                              return false;
                          }
                      }
                  }


                if(m_uiDim==3) {
                    if ((myZ + k * mySz) < domain_max) {
                        tmp = ot::TreeNode(myX, myY, (myZ + k * mySz), (OCT_KEY_FRONT | m_uiMaxDepth), m_uiDim, m_uiMaxDepth);
                        lookUPIndex = E2EMap[i*k_s*numDirections+OCT_DIR_FRONT * k_s + k - 1];
                        //assert(lookUPIndex!=UINT_MAX);
                        if (lookUPIndex != UINT_MAX) {
                            assert(lookUPIndex < allNodes.size());
                            lookUp = allNodes[lookUPIndex];
                            if (!(tmp.isAncestor(lookUp)) && (tmp != lookUp) && (!lookUp.isAncestor(tmp))) {
                                std::cout << "rank: " << rank << " key front: " << tmp << " lookup: " << lookUp << std::endl;
                                return false;
                            }
                        }
                    }


                    if (myZ >0) {
                        tmp = ot::TreeNode(myX, myY, (myZ - 1), (OCT_KEY_BACK | m_uiMaxDepth), m_uiDim, m_uiMaxDepth);
                        lookUPIndex = E2EMap[i*k_s*numDirections+OCT_DIR_BACK * k_s + k - 1];
                        //assert(lookUPIndex!=UINT_MAX);
                        if (lookUPIndex != UINT_MAX) {
                            assert(lookUPIndex < allNodes.size());
                            lookUp = allNodes[lookUPIndex];
                            if (!(tmp.isAncestor(lookUp)) && (tmp != lookUp) && (!lookUp.isAncestor(tmp))) {
                                std::cout << "rank: " << rank << " key back: " << tmp << " lookup: " << lookUp << std::endl;
                                return false;
                            }
                        }
                    }
                }

              }

          }

          return true;
      }




    template <typename T>
    bool checkE2NMapping( const std::vector<unsigned int >&  E2EMap , const std::vector<unsigned int >& E2NMap, const std::vector<T>& allNodes,unsigned int numDirections,unsigned int elementOrder)
      {

          unsigned int numNpE;
          if(m_uiDim==2)
              numNpE=(elementOrder+1)*(elementOrder+1);
          else if(m_uiDim==3)
              numNpE=(elementOrder+1)*(elementOrder+1)*(elementOrder+1);

          // note that E2NMap should be an array size of allNodes.size()*numNpE;
          unsigned  int nodeIndex=0;
          unsigned  int elementIndex=0;
          bool neighState;

          // perform iteration over local nodes and their neighbours, to see whether we have any invalid nodes.

          for(unsigned int e=0;e<allNodes.size();e++)
          {
              if(E2NMap[e*numNpE]!=UINT_MAX) { // This is because we invalidate some ghost element nodes who doesn't have any neighbours.
                  for (unsigned int k = 0; k < numNpE; k++) {
                      nodeIndex = E2NMap[e * numNpE + k];
                      assert(nodeIndex!=UINT_MAX); // This is to check whether we have partial invalidations. Partial invalidations are not allowed.
                      /*Note : This is true, because we can grantee that (k/(d+1)+j/(d+1)^2 +i/(d+1)^3 < 1)*/
                      elementIndex = nodeIndex /((elementOrder+1)*(elementOrder+1)*(elementOrder+1)); //(elementOrder + 1) / (elementOrder + 1) / (elementOrder + 1);
                      assert(allNodes[elementIndex].getLevel() <= allNodes[e].getLevel());
                      if (e != elementIndex) {
                          neighState = false;
                          for (unsigned int dir = 0; dir < numDirections; dir++) {
                              if (elementIndex == E2EMap[e * numDirections + dir]) {
                                  neighState = true;
                                  break;
                              }
                              /*  if(e==E2EMap[elementIndex*numDirections+dir])
                                {
                                    neighState=true;
                                    break;
                                }*/
                          }

                          /* if(!neighState)
                           {
                              std::cout<<"e1: "<<e<<" octant: "<<allNodes[e]<<std::endl;
                               for(unsigned int dir=0;dir<numDirections;dir++)
                               {
                                   std::cout<<" "<<E2EMap[e*numDirections+dir];
                               }
                               std::cout<<std::endl;
                               std::cout<<"e2: "<<elementIndex<<" octant: "<<allNodes[elementIndex]<<std::endl;
                               for(unsigned int dir=0;dir<numDirections;dir++)
                               {
                                   std::cout<<" "<<E2EMap[elementIndex*numDirections+dir];
                               }
                               std::cout<<std::endl;

                           }*/

                          // Note: This is due to the fact that, if DG index of a node doesn't belong to that element those elements should have a commom boundary between smaller octant and larger octant.
                          assert(neighState);
                          if(!neighState) return false;
                      }

                  }
              }

          }
          return  true;


      }



}//end namespace
}//end namespace

namespace par {
namespace test {

    template <typename T>
    bool isSorted(const std::vector<T > &nodes, MPI_Comm comm) {
      bool localPassed = seq::test::isSorted<T>(nodes);
      bool allLocalsPassed;

      par::Mpi_Allreduce<bool>(&localPassed, &allLocalsPassed, 1,
                               par::Mpi_datatype<bool>::LAND(), comm);

      if(allLocalsPassed) {
        bool failedParCheck = false;

        MPI_Comm new_comm;
        par::splitComm2way(nodes.empty(), &new_comm, comm);

        if(!nodes.empty()) {
          int rank;
          int npes;
          MPI_Request request;
          MPI_Status status;
          MPI_Comm_rank(new_comm,&rank);
          MPI_Comm_size(new_comm,&npes);

          //Send last to the next proc.
          T end = nodes[nodes.size()-1];
          if(rank < (npes -1)) {
            par::Mpi_Issend<T>(&end, 1, rank+1, 0, new_comm, &request );
          }

          T prev, me;
          me = nodes[0];

          //Recv prev from the prev proc.d
          if(rank) {
            par::Mpi_Recv<T>( &prev, 1, rank-1, 0, new_comm, &status );
            if(prev > me) {
              std::cout<<"rank "<<rank<<" prev: "<<prev<<" me: "<<me<<std::endl;
              failedParCheck = true;
            }
          }

          if(rank < (npes-1)) {
            MPI_Status statusWait;
            MPI_Wait(&request, &statusWait);
          }
        }//end if nodes not empty
        MPI_Comm_free(&new_comm);
        bool anyProcFailed;
        par::Mpi_Allreduce<bool>(&failedParCheck, &anyProcFailed, 1,
                                 par::Mpi_datatype<bool>::LOR(), comm);

        return (!anyProcFailed);
      }

      return allLocalsPassed;
    }


    template <typename T>
    bool isUniqueAndSorted(const std::vector<T > &nodes, MPI_Comm comm) {

      bool localPassed = seq::test::isUniqueAndSorted<T>(nodes);
      bool allLocalsPassed;

      par::Mpi_Allreduce<bool>(&localPassed, &allLocalsPassed, 1,
          par::Mpi_datatype<bool>::LAND(), comm);

      if(allLocalsPassed) {
        bool failedParCheck = false;

        MPI_Comm new_comm;
        par::splitComm2way(nodes.empty(), &new_comm, comm);

        if(!nodes.empty()) {
          int rank;
          int npes;
          MPI_Request request;
          MPI_Status status;
          MPI_Comm_rank(new_comm,&rank);
          MPI_Comm_size(new_comm,&npes);

          //Send last to the next proc.
          T end = nodes[nodes.size()-1];
          if(rank < (npes-1)) {
            par::Mpi_Issend<T>(&end, 1, rank+1, 0, new_comm, &request );
          }

          T prev, me;
          me = nodes[0];

          //Recv prev from the prev proc.
          if(rank) {
            par::Mpi_Recv<T>( &prev, 1, rank-1, 0, new_comm, &status );
            if((prev >= me)) {
              std::cout<<rank<<" from isUnique test : prev: "<<prev<<" me: "<<me<<std::endl;
              failedParCheck = true;
            }

          }

          if(rank < (npes-1)) {
            MPI_Status statusWait;
            MPI_Wait(&request, &statusWait);
          }
        }//end if nodes not empty
        MPI_Comm_free(&new_comm);
        bool anyProcFailed;
        par::Mpi_Allreduce<bool>(&failedParCheck, &anyProcFailed, 1,
            par::Mpi_datatype<bool>::LOR(), comm);

        return (!anyProcFailed);
      }

      return allLocalsPassed;
    }



    template <typename T>
    bool containsAncestor(const std::vector<T > &nodes, MPI_Comm comm) {

      bool localPassed = seq::test::containsAncestor<T>(nodes);
      bool allLocalsPassed;

      par::Mpi_Allreduce<bool>(&localPassed, &allLocalsPassed, 1,
                               par::Mpi_datatype<bool>::LAND(), comm);

      if(allLocalsPassed) {
        bool failedParCheck = false;

        MPI_Comm new_comm;
        par::splitComm2way(nodes.empty(), &new_comm, comm);

        if(!nodes.empty()) {
          int rank;
          int npes;
          MPI_Request request;
          MPI_Status status;
          MPI_Comm_rank(new_comm,&rank);
          MPI_Comm_size(new_comm,&npes);

          //Send last to the next proc.
          T end = nodes[nodes.size()-1];
          if(rank < (npes-1)) {
            par::Mpi_Issend<T>(&end, 1, rank+1, 0, new_comm, &request );
          }

          T prev, me;
          me = nodes[0];

          //Recv prev from the prev proc.
          if(rank) {
            par::Mpi_Recv<T>( &prev, 1, rank-1, 0, new_comm, &status );
            if((prev.isAncestor(me))) {
              std::cout<<rank<<" par::test::containsAncestor failed: prev: "<<prev<<" me: "<<me<<std::endl;
              failedParCheck = true;
            }

          }

          if(rank < (npes-1)) {
            MPI_Status statusWait;
            MPI_Wait(&request, &statusWait);
          }
        }//end if nodes not empty
        MPI_Comm_free(&new_comm);
        bool anyProcFailed;
        par::Mpi_Allreduce<bool>(&failedParCheck, &anyProcFailed, 1,
                                 par::Mpi_datatype<bool>::LOR(), comm);

        return (!anyProcFailed);
      }

      return allLocalsPassed;
    }



    template <typename T>
    bool  isComplete (const std::vector<T>& nodes, MPI_Comm comm) {

        DendroIntL vol=0;

        int rank,npes;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm,&npes);

        int maxDepth=0;
        int maxDepth_g=0;
        if(nodes.size()!=0) {

        maxDepth = nodes[0].getMaxDepth();
        int dim = nodes[0].getDim();

        vol = 0;
        int len = 0;

        for (int i = 0; i < nodes.size(); i++) {
            len = 1u << (maxDepth - nodes[i].getLevel());
            if (i < (nodes.size() - 1) && (nodes[i].isAncestor(nodes[i + 1]))) {
            std::cout<<"Rank:"<<rank<<" contains duplicate nodes."<<std::endl;
            return false;
            }
            vol += len * len * len;
        }
        }else{
        vol=0;
        }


        DendroIntL g_vol=0;
        par::Mpi_Reduce(&vol,&g_vol,1,MPI_SUM,0,comm);
        par::Mpi_Reduce(&maxDepth,&maxDepth_g,1,MPI_MAX,0,comm);

        if(!rank)
        {
            std::cout<<"MaxDepth of the Octree:"<<maxDepth_g<<std::endl;
            assert(maxDepth_g!=0);
            long long int max_len=1u<<maxDepth_g;
            if(g_vol==(max_len*max_len*max_len))
            {
            return true;
            }else
            {
            std::cout<<"Volume of the Complete Octree:"<<(max_len*max_len*max_len)<<std::endl;
            std::cout<<"Computed octree volume:"<<g_vol<<std::endl;
            return false;
            }
        }

    }


  }//end namespace
}//end namespace



namespace ot
{
    namespace test
    {

        template<typename T>
        bool isBalancedInternal(unsigned int dim, unsigned int maxDepth,
                                char*failFileName,	const std::vector<T> & nodes,
                                TreeNode holder, bool incCorn, unsigned int maxLevDiff)
        {


            bool yesBalanced = true;
            std::vector<TreeNode> failedCorners;
            std::vector<TreeNode> failedEdges;
            std::vector<TreeNode> failedFaces;
            TreeNode root (dim, maxDepth);
            unsigned int retIdx;
            double nxz = (double) nodes.size();
            for (unsigned int i=0; i<nodes.size(); i++) {
                /*
                   if (!(i%100)) {
                   printf("%4.2f\r", 100.0*((double)i)/nxz);
                   }
                   */

                // cout<<RED<<"R"<<NRM<<endl;
                {
                    TreeNode it = nodes[i].getRight();

                    if(holder.isAncestor(it)){

                        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                        if(found) {
                            unsigned int retLev = nodes[retIdx].getLevel();
                            unsigned int myLev = nodes[i].getLevel();
                            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                found = false;
                            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                //Found a very big Neighbour
                                found = false;
                            }
                        }
                        if(!found) {
                            yesBalanced = false;
                            std::cout<<nodes[i]<<": (R) ->"<<nodes[retIdx]<<std::endl;
                            failedFaces.push_back(nodes[i]);
                            //break;
                        }
                    }

                }

                // cout<<RED<<"L"<<NRM<<endl;
                {
                    TreeNode it = nodes[i].getLeft();

                    if(holder.isAncestor(it)){

                        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                        if(found) {
                            unsigned int retLev = nodes[retIdx].getLevel();
                            unsigned int myLev = nodes[i].getLevel();
                            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                found = false;
                            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                //Found a very big Neighbour
                                found = false;
                            }
                        }
                        if(!found) {
                            yesBalanced = false;
                            std::cout<<nodes[i]<<": (L) ->"<<nodes[retIdx]<<std::endl;
                            failedFaces.push_back(nodes[i]);
                            //break;
                        }
                    }

                }


                // cout<<RED<<"T"<<NRM<<endl;
                {
                    TreeNode it = nodes[i].getTop();

                    if(holder.isAncestor(it)){

                        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                        if(found) {
                            unsigned int retLev = nodes[retIdx].getLevel();
                            unsigned int myLev = nodes[i].getLevel();
                            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                found = false;
                            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                //Found a very big Neighbour
                                found = false;
                            }
                        }
                        if(!found) {
                            yesBalanced = false;
                            std::cout<<nodes[i]<<": (T) ->"<<nodes[retIdx]<<std::endl;
                            failedFaces.push_back(nodes[i]);
                            //break;
                        }
                    }
                }

                // cout<<RED<<"Bo"<<NRM<<endl;
                {
                    TreeNode it = nodes[i].getBottom();


                    if(holder.isAncestor(it)){

                        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                        if(found) {
                            unsigned int retLev = nodes[retIdx].getLevel();
                            unsigned int myLev = nodes[i].getLevel();
                            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                found = false;
                            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                //Found a very big Neighbour
                                found = false;
                            }
                        }
                        if(!found) {
                            yesBalanced = false;
                            std::cout<<nodes[i]<<": (Bo) ->"<<nodes[retIdx]<<std::endl;
                            failedFaces.push_back(nodes[i]);
                            // break;
                        }
                    }

                }

                // cout<<RED<<"F"<<NRM<<endl;
                {
                    TreeNode it = nodes[i].getFront();


                    if(holder.isAncestor(it)){

                        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                        if(found) {
                            unsigned int retLev = nodes[retIdx].getLevel();
                            unsigned int myLev = nodes[i].getLevel();
                            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                found = false;
                            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                //Found a very big Neighbour
                                found = false;
                            }
                        }
                        if(!found) {
                            yesBalanced = false;
                            std::cout<<nodes[i]<<": (F) ->"<<nodes[retIdx]<<std::endl;
                            failedFaces.push_back(nodes[i]);
                            // break;
                        }
                    }

                }

                // cout<<RED<<"Bk"<<NRM<<endl;
                {
                    TreeNode it = nodes[i].getBack();


                    if(holder.isAncestor(it)){

                        bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                        if(found) {
                            unsigned int retLev = nodes[retIdx].getLevel();
                            unsigned int myLev = nodes[i].getLevel();
                            if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                found = false;
                            }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                //Found a very big Neighbour
                                found = false;
                            }
                        }
                        if(!found) {
                            yesBalanced = false;
                            std::cout<<nodes[i]<<": (Bk) ->"<<nodes[retIdx]<<std::endl;
                            failedFaces.push_back(nodes[i]);
                            //break;
                        }
                    }

                }

                if(dim == 3 || incCorn) {
                    // cout<<RED<<"TR"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopRight();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TR) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"TL"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopLeft();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TL) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                // break;
                            }
                        }

                    }

                    // cout<<RED<<"BoL"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomLeft();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoL) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"BoR"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomRight();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoR) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                // break;
                            }
                        }

                    }

                    // cout<<RED<<"RBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getRightBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (RBk) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"TBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TBk) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                // break;
                            }
                        }

                    }

                    // cout<<RED<<"LF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getLeftFront();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (LF) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"LBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getLeftBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (LBk) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"RF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getRightFront();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (RF) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"TF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopFront();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TF) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"BoBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoBk) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"BoF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomFront();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoF) ->"<<nodes[retIdx]<<std::endl;
                                failedEdges.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                }//end if dim=3 or incCorn

                if(incCorn) {
                    // cout<<RED<<"TRBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopRightBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TRBk) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"TLF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopLeftFront();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TLF) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"BoRBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomRightBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoRBk) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"TLBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopLeftBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TLBk) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                // break;
                            }
                        }

                    }

                    // cout<<RED<<"BoRF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomRightFront();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoRF) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }

                    // cout<<RED<<"TRF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getTopRightFront();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (TRF) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                //break;
                            }
                        }

                    }


                    // cout<<RED<<"BoLBk"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomLeftBack();


                        if(holder.isAncestor(it)){

                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoLBk) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                //break;
                            }
                        }


                    }


                    // cout<<RED<<"BoLF"<<NRM<<endl;
                    {
                        TreeNode it = nodes[i].getBottomLeftFront();


                        if(holder.isAncestor(it)){


                            bool found = seq::maxLowerBound(nodes, it.getDFDMorton(), retIdx, NULL, NULL);
                            if(found) {
                                unsigned int retLev = nodes[retIdx].getLevel();
                                unsigned int myLev = nodes[i].getLevel();
                                if( (it.getAnchor() != nodes[retIdx].getAnchor()) && (!(nodes[retIdx].isAncestor(it))) ) {
                                    found = false;
                                }else if( (retLev < myLev) && (myLev - retLev) > maxLevDiff  ) {
                                    //Found a very big Neighbour
                                    found = false;
                                }
                            }
                            if(!found) {
                                yesBalanced = false;
                                std::cout<<nodes[i]<<": (BoLF) ->"<<nodes[retIdx]<<std::endl;
                                failedCorners.push_back(nodes[i]);
                                // break;
                            }
                        }

                    }

                }//end if incCorn

            }//end for i

/*                seq::makeVectorUnique(failedFaces,false);
                seq::makeVectorUnique(failedEdges,false);
                seq::makeVectorUnique(failedCorners,false);*/

            treeNodesTovtk(failedFaces,0,"failedFaces");
            treeNodesTovtk(failedEdges,0,"failedEdges");
            treeNodesTovtk(failedCorners,0,"failedCorners");


            if(!yesBalanced)
            {
                std::cout << BLU << "===============================================" << NRM << std::endl;
                std::cout << RED " Balance test failed. " NRM << std::endl;
                std::cout << YLW << "Failed Face Cases:"<<failedFaces.size()<<NRM<<std::endl;
                std::cout << YLW << "Failed Edges Cases:"<<failedEdges.size()<<NRM<<std::endl;
                std::cout << YLW << "Failed Corners Cases:"<<failedCorners.size()<<NRM<<std::endl;
                std::cout << BLU << "===============================================" << NRM << std::endl;

            }


            char failCornerFileName[100];
            char failEdgeFileName[100];
            char failFaceFileName[100];
            strcpy(failFaceFileName,failFileName );
            strcpy(failEdgeFileName,failFileName);
            strcpy(failCornerFileName,failFileName);
            strcat(failFaceFileName,"_Faces.ot\0");
            strcat(failEdgeFileName,"_Edges.ot\0");
            strcat(failCornerFileName,"_Corners.ot\0");
            /*writeNodesToFile(failCornerFileName,failedCorners);
            writeNodesToFile(failEdgeFileName,failedEdges);
            writeNodesToFile(failFaceFileName,failedFaces);*/
            return yesBalanced;


        }



        template <typename T>
        bool isBalanced(unsigned int dim, unsigned int maxDepth, char* failFileName,
                        const std::vector<T>& nodes, bool incCorn, unsigned int maxLevDiff)
        {

            TreeNode root (dim, maxDepth);
            return	isBalancedInternal(dim,maxDepth,failFileName,nodes,root, incCorn, maxLevDiff);

        }


        template<typename T,typename B>
        bool isBlockListValid(const std::vector<T>& pNodes, std::vector<B> & blockList,unsigned int d_min, unsigned int d_max,unsigned int nBegin,unsigned int nEnd)
        {


           unsigned int numOctLevOctree=0;
           unsigned int numOctLevBlock=0;

           int rank;
           MPI_Comm_rank(MPI_COMM_WORLD,&rank);
           /*std::vector<ot::TreeNode> blockLevOcts;
           std::vector<ot::TreeNode> octreeLevOcts;*/

           unsigned int numIdealOct=0;
           for(unsigned int lev=d_min; lev<=d_max;lev++)
           {
               numOctLevOctree=0;
               numOctLevBlock=0;

              /* blockLevOcts.clear();
               octreeLevOcts.clear();*/

               for(unsigned int k=nBegin;k<nEnd;k++)
               {
                   if(pNodes[k].getLevel()==lev)
                   {
                       numOctLevOctree++;
                       //octreeLevOcts.push_back(pNodes[k]);
                   }


               }

               for(unsigned int k=0;k<blockList.size();k++)
               {

                  if(blockList[k].getRegularGridLev()==lev) {
                      for (unsigned int j = blockList[k].getLocalElementBegin();
                           j < blockList[k].getLocalElementEnd(); j++) {

                          if(pNodes[j].getLevel()==lev)
                          {
                              numOctLevBlock++;
                              //blockLevOcts.push_back(pNodes[j]);
                          }

                      }
                  }

               }

               //assert(numOctLevBlock==numOctLevOctree);
               if(numOctLevBlock!=numOctLevOctree) {

                  /* treeNodesTovtk(octreeLevOcts,rank,"octreeLevOcts");
                   treeNodesTovtk(blockLevOcts,rank,"blockLevOcts");*/
                   std::cout << "rank: "<<rank<< " num oct mismatch for lev: " << lev << " numOctLevBlock: " << numOctLevBlock << " numOctLevOctree: " << numOctLevOctree << std::endl;

                   return false;
               }

            }
            return true;


        }




    } //end namespace test


} //end namespace ot.




