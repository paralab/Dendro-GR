
/**
@brief A collection of simple functions for manipulating octrees.
Examples: Regular Refinements, Linearizing an octree, I/O, 
Nearest Common Ancestor, adding positive boundaries, marking hanging nodes 
@author Rahul S. Sampath, rahul.sampath@gmail.com
@author Hari Sundar, hsundar@gmail.com
*/

#include "TreeNode.h"
#include "parUtils.h"
#include "seqUtils.h"
#include <cstring>

#ifdef __DEBUG__
#ifndef __DEBUG_OCT__
#define __DEBUG_OCT__
#endif
#endif

#ifdef __DEBUG_OCT__
#ifndef __MEASURE_FLAG_NODES__
#define __MEASURE_FLAG_NODES__
#endif
#endif

namespace ot {

  //inOct1 and inOct2 must be globally sorted and complete
  int mergeOctrees(std::vector<TreeNode>& inOct1, std::vector<TreeNode>& inOct2, std::vector<TreeNode>& outOct, MPI_Comm comm) {
    
    PROF_MERGE_OCTREES_BEGIN
    int rank;
    int npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    assert(!(inOct1.empty()));

    //Get the mins of inOct1
    ot::TreeNode* mins = new ot::TreeNode[npes];
    assert(mins);
    par::Mpi_Allgather<ot::TreeNode>(&(*(inOct1.begin())), mins, 1, comm);

    //Distribute inOct2 to align with inOct1
    int* sendCnts = new int[npes];
    assert(sendCnts);
    for(int i = 0; i < npes; i++) {
      sendCnts[i] = 0;
    }//end for i 

    if(!(inOct2.empty())) {
      //Since inOct1 is complete and sorted, mins[0] is a deepest first
      //decendant of root.
      //Since inOct2 is also complete and sorted, on the first non-empty processor
      //inOct2[0] will also be a deepest first decendant of root.
      //Now mins[0] can be an ancestor, decendant or equal to inOct2[0] (global
      //first element, so only happens on the first processor). So i=0
      //is handled separately. 
      int minsCnt = 0;
      if(inOct2[0] < mins[0]) {
        sendCnts[0]++;
      } else {
        while( (minsCnt < npes) && (mins[minsCnt] <= inOct2[0]) ) {
          minsCnt++;
        }
        sendCnts[minsCnt - 1]++;
      }
      for(int i = 1; i < inOct2.size(); i++) {
        while( (minsCnt < npes) && (mins[minsCnt] <= inOct2[i]) ) {
          minsCnt++;
        }
        sendCnts[minsCnt - 1]++;
      }//end for i
    }

    int* recvCnts = new int[npes];
    assert(recvCnts);

    par::Mpi_Alltoall<int>(sendCnts, recvCnts, 1, comm);

    int* sendDisps = new int[npes];
    assert(sendDisps);

    int* recvDisps = new int[npes];
    assert(recvDisps);

    sendDisps[0] = 0;
    recvDisps[0] = 0;
    for(int i = 1; i < npes; i++) {
      sendDisps[i] = sendDisps[i - 1] + sendCnts[i - 1];
      recvDisps[i] = recvDisps[i - 1] + recvCnts[i - 1];
    }//end for i

    outOct.resize(recvDisps[npes - 1] + recvCnts[npes - 1]);

    par::Mpi_Alltoallv_sparse<ot::TreeNode>(&(*(inOct2.begin())), sendCnts, sendDisps,
        &(*(outOct.begin())), recvCnts, recvDisps, comm);

    delete [] mins;
    delete [] sendCnts;
    delete [] sendDisps;
    delete [] recvCnts;
    delete [] recvDisps;

    //Merge outOct and inOct1 locally
    std::vector<ot::TreeNode> tmpOct;
    int cnt1 = 0;
    int cnt2 = 0;
    while( (cnt1 < inOct1.size()) && (cnt2 < outOct.size()) ) {
      if(inOct1[cnt1] < outOct[cnt2]) {
        tmpOct.push_back(inOct1[cnt1]);
        cnt1++;
      } else {
        tmpOct.push_back(outOct[cnt2]);
        cnt2++;
      }
    }

    while(cnt1 < inOct1.size()) {
      tmpOct.push_back(inOct1[cnt1]);
      cnt1++;
    }

    while(cnt2 < outOct.size()) {
      tmpOct.push_back(outOct[cnt2]);
      cnt2++;
    }

    outOct = tmpOct;
    tmpOct.clear();

    //Linearize outOct. This will remove duplicates and ancestors
    lineariseList(outOct, comm);

    PROF_MERGE_OCTREES_END
  }//end function

  int refineOctree(const std::vector<ot::TreeNode> & inp, std::vector<ot::TreeNode> &out) {
    out.clear();
    for(unsigned int i = 0; i < inp.size(); i++) {
      if(inp[i].getLevel() < inp[i].getMaxDepth()) {
        inp[i].addChildren(out);
      } else {
        out.push_back(inp[i]);
      }
    }
    return 1;
  }//end function 

  int refineAndPartitionOctree(const std::vector<ot::TreeNode> & inp, std::vector<ot::TreeNode> &out, MPI_Comm comm) {
    refineOctree(inp,out);
    par::partitionW<ot::TreeNode>(out, NULL,comm);
    return 1;
  }//end function

  int createRegularOctree(std::vector<ot::TreeNode>& out, unsigned int lev, unsigned int dim, unsigned int maxDepth, MPI_Comm comm) {
    TreeNode root(dim,maxDepth);
    out.clear();
    int rank;
    MPI_Comm_rank(comm,&rank);
    if(!rank) {
      out.push_back(root);
    }
    for(int i = 0; i < lev; i++) {
      std::vector<ot::TreeNode> tmp;
      refineAndPartitionOctree(out,tmp,comm);
      out = tmp;
      tmp.clear();
    }
    return 1;
  }

  //list must be sorted.
  int lineariseList(std::vector<ot::TreeNode> & list, bool skipLast) {
    std::vector<ot::TreeNode> tmp;
    if(!(list.empty())) {
      for(unsigned int i = 0; i < (list.size()-1); i++) {
    #ifdef __DEBUG_OCT__
        assert(areComparable(list[i], list[i+1]));
    #endif
        if( (!(list[i].isAncestor(list[i+1]))) && (list[i] != list[i+1]) ) {
          tmp.push_back(list[i]);
        }
      }
      if(!skipLast) {
        tmp.push_back(list[list.size()-1]);
      }
    }

    list = tmp;
    tmp.clear();
    return 1;
  }//end fn.

  //list must be sorted.
  int lineariseList(std::vector<ot::TreeNode> & list, MPI_Comm comm) {
    int rank,size;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&size);

    if(size == 1) {
      lineariseList(list, false);
      return 1;
    }

    //Remove empty processors...
    int new_rank, new_size; 
    MPI_Comm   new_comm;
    par::splitComm2way(list.empty(), &new_comm,comm);

    MPI_Comm_rank (new_comm, &new_rank);
    MPI_Comm_size (new_comm, &new_size);
    if(!list.empty()) {
      //Send the last octant to the next processor. 
      ot::TreeNode lastOctant = list[list.size()-1];
      ot::TreeNode lastOnPrev;

      MPI_Request recvRequest;
      MPI_Request sendRequest;

      if(new_rank) {
        par::Mpi_Irecv<ot::TreeNode>( &lastOnPrev, 1, new_rank-1, 1, new_comm, &recvRequest);
      }
      if(new_rank < (new_size-1)) {
        par::Mpi_Issend<ot::TreeNode>( &lastOctant, 1, new_rank+1, 1, new_comm,  &sendRequest);
      }

      if(new_rank) {
        std::vector<ot::TreeNode> tmp(list.size()+1);
        for(int i = 0; i < list.size(); i++) {
          tmp[i+1] = list[i];
        }

        MPI_Status statusWait;
        MPI_Wait(&recvRequest, &statusWait);
        tmp[0] = lastOnPrev;

        list = tmp;
        tmp.clear();
      }

      if(new_rank == (new_size-1)) {
        lineariseList(list, false);
      }else {
        lineariseList(list, true);
      }

      if(new_rank < (new_size-1)) {
        MPI_Status statusWait;
        MPI_Wait(&sendRequest, &statusWait);
      }      
    }//not empty procs only

    MPI_Comm_free(&new_comm);

    return 1;
  }//end fn.

  bool lessThanUsingWts ( TreeNode  const & a,  TreeNode  const & b)  {
    if(a == b){
      return (a.getWeight() < b.getWeight());
    }else {
      return (a < b);
    }
  }//end function

  unsigned int getNodeWeight(const TreeNode * t) {
    return t->getWeight();
  }

  bool bPartComparator(TreeNode a, TreeNode b) {
    if ( (a.getWeight()) != (b.getWeight()) ) {
      //pick denser
      return( (a.getWeight()) > (b.getWeight()) );
    } else if (a.getWeight() > 1) {
      //pick finer
      return( (a.getLevel()) > (b.getLevel()) );
    } else {
      //pick coarser
      return( (a.getLevel()) < (b.getLevel()) );
    }
  }//end function

  //If one of first and second is an ancestor of the other, it is returned.  
  TreeNode getNCA(TreeNode first, TreeNode second) {
    #ifdef __DEBUG_OCT__
      assert(areComparable(first,second));
      assert(first != second);
    #endif
    unsigned int fx = first.getX();
    unsigned int sx = second.getX();
    unsigned int fy = first.getY();
    unsigned int sy = second.getY();
    unsigned int fz = first.getZ();
    unsigned int sz = second.getZ();
    unsigned int maxDepth = first.getMaxDepth(); 
    unsigned int dim = first.getDim(); 
    unsigned int maxDiff = (unsigned int)(std::max((std::max((fx^sx),(fy^sy))),(fz^sz)));
    unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
    //Eliminate the last maxDiffBinLen bits.
    unsigned int ncaX = ((fx>>maxDiffBinLen)<<maxDiffBinLen);
    unsigned int ncaY = ((fy>>maxDiffBinLen)<<maxDiffBinLen);
    unsigned int ncaZ = ((fz>>maxDiffBinLen)<<maxDiffBinLen);
    unsigned int ncaLev = (maxDepth - maxDiffBinLen);
    assert(ncaLev<std::min(first.getLevel(),second.getLevel()));

    //if(ncaLev>std::min(first.getLevel(),second.getLevel()))
    //  ncaLev=std::min(first.getLevel(),second.getLevel());

    TreeNode nca(ncaX,ncaY,ncaZ,ncaLev,dim,maxDepth);
    return nca;
  }//end function

  int readPtsFromFile(char* filename, std::vector<double>& pts) {
    FILE* infile;
    size_t res;
    infile = fopen(filename,"rb");
    unsigned int temp;
    res = fread(&temp,sizeof(unsigned int),1,infile);

    double* ptsTemp = NULL;

    if(temp) {
      ptsTemp = new double[3*temp];
      assert(ptsTemp);
      res = fread(ptsTemp, sizeof(double),3*temp,infile);
    }

    fclose(infile);
    pts.resize(3*temp);

    for (int i=0; i < (3*temp); i++) {
      pts[i] = ptsTemp[i];
    }//end for

    if(ptsTemp) {
      delete [] ptsTemp;
      ptsTemp = NULL;
    }

    // std::cout << __func__ << ": size " << temp << ", " << pts.size() << std::endl;
    return 1;
  }//end function



  int readDataPtsFromFile(char* filename, std::vector<double>& pts,
      std::vector<double>& ptVals) {
    // file format:
    // 4 bytes (unsigned int?)  number of points N
    // 3*N*8 bytes coordinates of point (double);  X1 Y1 Z1 X2 Y2 Z2 ....
    // N*8 weights attached to points (double)
    // 
    FILE* infile;
    size_t res;
    unsigned int temp;

    assert(sizeof(double)==8);
    assert(sizeof(unsigned int)==4);

    infile = fopen(filename,"rb");
    res = fread(&temp,sizeof(unsigned int),1,infile);

    pts.resize(3*temp);
    ptVals.resize(temp);

    res = fread(&(pts[0]), sizeof(double),3*temp,infile);
    res = fread(&(ptVals[0]), sizeof(double),temp,infile);

    fclose(infile);

    return 1;
  }//end function

  int writePtsToFile(char* filename, std::vector<double>& pts) {
    FILE* outfile = fopen(filename,"wb");
    unsigned int ptsLen = pts.size();
    double * ptsTemp = NULL;
    if(!pts.empty()) {
      ptsTemp = (&(*(pts.begin())));
    }

    if (ptsLen >0) {
      unsigned int numPts = ptsLen/3;
      fwrite(&numPts,sizeof(unsigned int),1,outfile);
      fwrite(ptsTemp, sizeof(double),ptsLen,outfile);
    }
    fclose(outfile);
    return 1;
  }//end function


  int writeDataPtsToFile(char* filename, std::vector<double>& pts,
      std::vector<double>& data) {
    FILE* outfile = fopen(filename,"wb");
    unsigned int ptsLen = pts.size();
    unsigned int numPts=ptsLen/3;

    assert(sizeof(double)==8);
    assert(sizeof(unsigned int)==4);

    fwrite(&numPts,sizeof(unsigned int),1,outfile);
    fwrite(&(pts[0]), sizeof(double),ptsLen,outfile);
    //fwrite(&(data[0]), sizeof(double),numPts,outfile);
    fclose(outfile);
    return 1;
  }//end function


  int readNodesFromFile_binary (char* filename, std::vector<TreeNode > & nodes) {

     int res;
     FILE* infile = fopen(filename,"r");

     ot::TreeNode tmp;
     do
     {
       fread(&tmp,sizeof(ot::TreeNode),1,infile);
       if(tmp.getDim()!=3)
         std::cout<<"Dim error:"<<tmp<<std::endl;
       nodes.push_back(tmp);

     }while(!feof(infile));





  }

  int readNodesFromFile (char* filename, std::vector<TreeNode > & nodes) {
    int res;
    FILE* infile = fopen(filename,"r");
    unsigned int numNode;
    unsigned int dim, maxDepth;
    res = fscanf(infile,"%u",&dim);
    res = fscanf(infile,"%u",&maxDepth); 
    res = fscanf(infile,"%u",&numNode);
    nodes.resize(numNode) ;   

    for (unsigned int i =0; i< nodes.size(); i++) {
      unsigned int x,y,z,d;
      res = fscanf(infile,"%u",&x);
      res = fscanf(infile,"%u",&y);
      res = fscanf(infile,"%u",&z);
      res = fscanf(infile,"%u",&d); 
      nodes[i] = ot::TreeNode (x,y,z,d, dim, maxDepth);
    }
    fclose(infile);
    return 1;
  }//end function

  int writeNodesToFile_binary(char * filename,const std::vector<TreeNode> & nodes)
  {
    FILE* outfile = fopen(filename,"w");

    for(int i=0;i<nodes.size();i++)
    {
      fwrite(&nodes[i],sizeof(ot::TreeNode),1,outfile);
    }

    fclose(outfile);

    return 0;
  }

  int writeNodesToFile (char* filename, const std::vector<TreeNode> & nodes) {
    FILE* outfile = fopen(filename,"w");
    if (!nodes.empty()) {
      unsigned int dim = nodes[0].getDim();
      unsigned int maxDepth = nodes[0].getMaxDepth();
      fprintf(outfile,"%u %u\n",dim,maxDepth); 
      fprintf(outfile,"%u\n",static_cast<unsigned int>(nodes.size()));
      for (unsigned int i =0; i< nodes.size(); i++) {
        assert(nodes[i].getDim() == dim);
        assert(nodes[i].getMaxDepth() == maxDepth);
        fprintf(outfile,"%u %u %u %u\n",nodes[i].getX(),nodes[i].getY(),nodes[i].getZ(),nodes[i].getLevel());
      }
    }
    fclose(outfile);
    return 1;
  }//end function

  bool areComparable (TreeNode first, TreeNode second) {
    return( ( (first.getDim()) == (second.getDim()) ) && ( (first.getMaxDepth()) == (second.getMaxDepth()) ) );
  }

  int int2str(int n,char*s) {
    int tmpd[20];
    int i=0;
    int j;   
    if (n==0) {
      strcpy(s,"0\0");
    } else {
      while (n>0) {
        tmpd[i]= (n%10);
        n= (int)(n/10);
        i++;
      }
      for (j=i-1;j>=0;j--) {
        s[i-j-1]=int2char(tmpd[j]);
      }
      s[i]='\0';
    }
    return 1;
  }//end function

  char int2char(int d) {
    switch (d) {
      case 0: return '0';
      case 1: return '1';
      case 2: return '2';
      case 3: return '3';
      case 4: return '4';
      case 5: return '5';
      case 6: return '6';
      case 7: return '7';
      case 8: return '8';
      case 9: return '9';
      default: return '\0';
    }
  }//end function

  // This will add boundary nodes and will also embed the octree one level higher
  // to enable the addition of the boundary nodes.
  void addBoundaryNodesType2(std::vector<ot::TreeNode> &in, 
      std::vector<ot::TreeNode>& bdy, 
      unsigned int dim, unsigned int maxDepth) {

    PROF_ADD_BDY_BEGIN
      assert(bdy.empty());

    for (unsigned int i = 0; i < in.size(); i++) {
      // get basic info ...
      unsigned int d   = in[i].getLevel();
      unsigned int x = in[i].getX();
      unsigned int y = in[i].getY();
      unsigned int z = in[i].getZ();

      unsigned char bdyFlags;
      // check if this is a boundary octant or not ...
      if ( in[i].isBoundaryOctant(ot::TreeNode::POSITIVE, &bdyFlags) ) {
        // bdy flags tells us which octants to add ...

        //NOTE: == is important since a&(b+c) will be true if
        //a=b, a=c and a=b+c

        // +x and more ... add additional octant in +x dir
        if ( bdyFlags & ot::TreeNode::X_POS_BDY ) {
          bdy.push_back(ot::TreeNode( (1u << maxDepth), y, z, (d+1), dim, maxDepth+1));
        }

        // +y and more ... add additional octant in +y dir
        if ( bdyFlags & ot::TreeNode::Y_POS_BDY ) {
          bdy.push_back(ot::TreeNode(x, (1u << maxDepth), z,
                (d+1), dim, maxDepth+1));
        }

        // +z and more ... add additional octant in +z dir
        if ( bdyFlags & ot::TreeNode::Z_POS_BDY ) {
          bdy.push_back(ot::TreeNode(x, y, (1u << maxDepth),
                (d+1), dim, maxDepth+1));
        }

        //+x+y and more
        if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY)) == 
            (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY) ) { 
          bdy.push_back(ot::TreeNode((1u << maxDepth),(1u << maxDepth), z,
                (d+1), dim, maxDepth+1));
        }

        //+x+z and more
        if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Z_POS_BDY)) ==
            (ot::TreeNode::X_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {
          bdy.push_back(ot::TreeNode((1u << maxDepth), y, (1u << maxDepth),
                (d+1), dim, maxDepth+1));
        }

        //+y+z and more
        if ( (bdyFlags & (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY)) == 
            (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {
          bdy.push_back(ot::TreeNode(x, (1u << maxDepth), (1u << maxDepth),
                (d+1), dim, maxDepth+1));
        }

        // if global corner ...
        //+x +y and +z only
        if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY +
                ot::TreeNode::Z_POS_BDY)) == (ot::TreeNode::X_POS_BDY +
                ot::TreeNode::Y_POS_BDY +  ot::TreeNode::Z_POS_BDY) ) {
          bdy.push_back(ot::TreeNode((1u << maxDepth), (1u << maxDepth), (1u << maxDepth),
                (d+1), dim, maxDepth+1));
        }
      }//end if boundary

      // Embed the actual octant in one level higher ...
      in[i] = ot::TreeNode(x, y, z, d+1, dim, maxDepth+1);

    }//end for i

    // A Parallel Sort for the bdy nodes follows.  
    //Then in and bdy will be merged.

    PROF_ADD_BDY_END
  }//end function

  void markBoundaryNodesAtAllLevels(std::vector<ot::TreeNode>& finestOctree, unsigned int nlevels,
      std::vector<ot::TreeNode>* coarserOctrees, unsigned int maxDepth) {
    //coarser octrees
    for(int lev = 0; lev < (nlevels - 1); lev++) {
      for(unsigned int i = 0; i < coarserOctrees[lev].size(); i++) {
        unsigned int myX = coarserOctrees[lev][i].getX();
        unsigned int myY = coarserOctrees[lev][i].getY();
        unsigned int myZ = coarserOctrees[lev][i].getZ();
        if( (myX == (1u << (maxDepth - 1))) ||
            (myY == (1u << (maxDepth - 1))) ||
            (myZ == (1u << (maxDepth - 1))) ) {
          coarserOctrees[lev][i].orFlag(ot::TreeNode::BOUNDARY);
        }
      }//end for i
    }//end for lev

    //finest octree
    for(unsigned int i = 0; i < finestOctree.size(); i++) {
      unsigned int myX = finestOctree[i].getX();
      unsigned int myY = finestOctree[i].getY();
      unsigned int myZ = finestOctree[i].getZ();
      if( (myX == (1u << (maxDepth - 1))) ||
          (myY == (1u << (maxDepth - 1))) ||
          (myZ == (1u << (maxDepth - 1))) ) {
        finestOctree[i].orFlag(ot::TreeNode::BOUNDARY);
      }
    }//end for i

  }//end function

  void discardExtraBoundaryOctants(std::vector<ot::TreeNode>& in,
      unsigned int dim, unsigned int maxDepth) {

    std::vector<ot::TreeNode> tmpOctree;
    for(int i = 0; i < in.size(); i++) {
      unsigned int myX = in[i].getX();
      unsigned int myY = in[i].getY();
      unsigned int myZ = in[i].getZ();
      if( (myX <= (1u << (maxDepth - 1))) &&
          (myY <= (1u << (maxDepth - 1))) &&
          (myZ <= (1u << (maxDepth - 1))) ) {
        tmpOctree.push_back(in[i]);
      }
    }//end for i

    if(tmpOctree.size() < in.size()) {
      in = tmpOctree;
    }
  }//end function

  //A simple implementation for now. This just calls flagNodesType3 for all
  //levels. A smarter implementation would use the fact that regular coarse
  //grid nodes remain regular on all finer octrees 
  void markHangingNodesAtAllLevels(std::vector<ot::TreeNode>& finestOctree, unsigned int nlevels, 
      std::vector<ot::TreeNode>* coarserOctrees, MPI_Comm* activeComms,
      unsigned int dim, unsigned int maxDepth) {

    for(int lev = 0; lev < (nlevels - 1); lev++) {
      if( !(coarserOctrees[lev].empty()) ) {
        ot::flagNodesType3(coarserOctrees[lev], activeComms[lev + 1]);
      }
    }//end for lev

    if( !(finestOctree.empty()) ) {
      ot::flagNodesType3(finestOctree, activeComms[0]);
    }

  }//end function

  // This will add boundary nodes and will also embed the octree one level higher
  // to enable the addition of the boundary nodes. The positive boundary nodes
  // are also marked as BOUNDARY.
  void addBoundaryNodesType1(std::vector<ot::TreeNode> &in, 
      std::vector<ot::TreeNode>& bdy, 
      unsigned int dim, unsigned int maxDepth) {
    PROF_ADD_BDY_BEGIN

      assert(bdy.empty());

    for (unsigned int i = 0; i < in.size(); i++) {
      // get basic info ...
      unsigned int d   = in[i].getLevel();
      unsigned int x = in[i].getX();
      unsigned int y = in[i].getY();
      unsigned int z = in[i].getZ();

      unsigned char bdyFlags;
      // check if this is a boundary octant or not ...
      if ( in[i].isBoundaryOctant(ot::TreeNode::POSITIVE, &bdyFlags) ) {
        // bdy flags tells us which octants to add ...

        //NOTE: == is important since a&(b+c) will be true if
        //a=b, a=c and a=b+c

        // +x and more ... add additional octant in +x dir
        if ( bdyFlags & ot::TreeNode::X_POS_BDY ) {
          bdy.push_back(ot::TreeNode( (1u << maxDepth), y, z, (d+1) |
                ot::TreeNode::BOUNDARY, dim, maxDepth+1));
        }

        // +y and more ... add additional octant in +y dir
        if ( bdyFlags & ot::TreeNode::Y_POS_BDY ) {
          bdy.push_back(ot::TreeNode(x, (1u << maxDepth), z,
                (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
        }

        // +z and more ... add additional octant in +z dir
        if ( bdyFlags & ot::TreeNode::Z_POS_BDY ) {
          bdy.push_back(ot::TreeNode(x, y, (1u << maxDepth),
                (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
        }

        //+x+y and more
        if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY)) == 
            (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY) ) { 
          bdy.push_back(ot::TreeNode((1u << maxDepth),(1u << maxDepth), z,
                (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
        }

        //+x+z and more
        if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Z_POS_BDY)) ==
            (ot::TreeNode::X_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {
          bdy.push_back(ot::TreeNode((1u << maxDepth), y, (1u << maxDepth),
                (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
        }

        //+y+z and more
        if ( (bdyFlags & (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY)) == 
            (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {
          bdy.push_back(ot::TreeNode(x, (1u << maxDepth),(1u << maxDepth),
                (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
        }

        // if global corner ...
        //+x +y and +z only
        if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY +
                ot::TreeNode::Z_POS_BDY)) == (ot::TreeNode::X_POS_BDY +
                ot::TreeNode::Y_POS_BDY +  ot::TreeNode::Z_POS_BDY) ) {
          bdy.push_back(ot::TreeNode((1u << maxDepth), (1u << maxDepth), (1u << maxDepth),
                (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
        }
      }//end if boundary

      // Embed the actual octant in one level higher ...
      in[i] = ot::TreeNode(x, y, z, d+1, dim, maxDepth+1);

    }//end for i

    // A Parallel Sort for the bdy nodes follows in the constructor.  
    //Then in and bdy will be merged.

    PROF_ADD_BDY_END
  }//end function

  void flagNodesType1(std::vector<ot::TreeNode> & in, MPI_Comm comm) {

    #ifdef __PROF_WITH_BARRIER__
        MPI_Barrier(comm);
    #endif

    PROF_MARK_HANGING_BEGIN

      PROF_FLN_STAGE1_BEGIN

      int npes;
    int rank;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    std::vector<unsigned int > keysCount(in.size());
    std::vector<ot::TreeNode > keys;

    assert(!in.empty());
    unsigned int maxD = in[0].getMaxDepth();
    unsigned int dim  = in[0].getDim();
    ot::TreeNode* inPtr = (&(*(in.begin())));

    //1. Generate Keys
    for (int i = 0; i < in.size(); i++) {
      keysCount[i] = 0;
      unsigned int myLev = inPtr[i].getLevel();
      if (myLev == 1) {
        continue;
      }
      unsigned int childNum = inPtr[i].getChildNumber();   
      unsigned int mySz = (1u << (maxD - myLev));
      unsigned int myX = inPtr[i].getX();
      unsigned int myY = inPtr[i].getY();
      unsigned int myZ = inPtr[i].getZ();    

      switch (childNum) {
        case 0:
          {
            break;
          }
        case 7:
          {
            break;
          }
        case 1:
          {
            //-y,-z,-yz
            if (myY) {
              TreeNode tmp(myX,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            if (myZ) {
              TreeNode tmp(myX,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            if (myY && myZ) {
              TreeNode tmp(myX,myY-mySz,myZ-mySz,myLev,dim,maxD);         
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            break;
          }
        case 2:
          {
            //-x,-z,-xz
            if (myX) {
              TreeNode tmp(myX-mySz,myY,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            if (myZ) {
              TreeNode tmp(myX,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            if (myX && myZ) {
              TreeNode tmp(myX-mySz,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            break;
          }
        case 3:
          {
            //-z
            if (myZ) {
              TreeNode tmp(myX,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            break;
          }
        case 4:
          {
            //-x,-y,-xy
            if (myX) {
              TreeNode tmp(myX-mySz,myY,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            if (myY) {
              TreeNode tmp(myX,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            if (myX && myY) {
              TreeNode tmp(myX-mySz,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            break;
          }
        case 5:
          {
            //-y
            if (myY) {
              TreeNode tmp(myX,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
              keysCount[i]++;
            }
            break;
          }
        case 6:
          {
            //-x
            if (myX) {
              TreeNode tmp(myX-mySz,myY,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());    
              keysCount[i]++;
            }
            break;      
          }
        default: assert(false);
      }//end switch
    }//end for i	

    PROF_FLN_STAGE1_END
      PROF_FLN_STAGE2_BEGIN

      // allocate memory for the mins array
      std::vector<ot::TreeNode> mins(npes);

    par::Mpi_Allgather<ot::TreeNode>(inPtr, &(*(mins.begin())), 1, comm);

    unsigned int *part = NULL;
    if(!keys.empty()) {
      part = new unsigned int[keys.size()];    
      assert(part);
    }

    for (unsigned int i = 0; i < keys.size(); i++) {
      unsigned int idx;
      //maxLB returns the last index in a sorted array
      //such that a[ind] <= key and  a[index +1] > key
      bool found = seq::maxLowerBound<TreeNode >(mins, keys[i], idx, NULL, NULL);
      if (!found ) {
        part[i] = rank;
      } else {
        part[i] = idx;
      }
    }//end for i
    mins.clear();

    PROF_FLN_STAGE2_END
      PROF_FLN_STAGE3_BEGIN

      int *numKeysSend = new int[npes];
    assert(numKeysSend);

    int *numKeysRecv = new int[npes];
    assert(numKeysRecv);
    for (int i = 0; i < npes; i++) {
      numKeysSend[i] = 0;
    }

    // calculate the number of keys to send ...
    for (unsigned int i = 0; i < keys.size(); i++) {
      assert(part[i] < npes);
      numKeysSend[part[i]]++;
    }

    // Now do an All2All to get inumKeysRecv
    par::Mpi_Alltoall<int>(numKeysSend, numKeysRecv, 1, comm);

    unsigned int totalKeys=0; // total number of local keys ...
    for (int i = 0; i < npes; i++) {
      totalKeys += numKeysRecv[i];
    }

    #ifdef __MEASURE_FLAG_NODES__
    MPI_Barrier(comm);
    int numProcsSend = 0;
    int numProcsRecv = 0;
    for(int i = 0; i < npes; i++) {
      if(numKeysSend[i]) {
        numProcsSend++;
      }
      if(numKeysRecv[i]) {
        numProcsRecv++;
      }
    }//end for i
    int* allNumProcsSend = new int[npes];
    assert(allNumProcsSend);

    int* allNumProcsRecv = new int[npes];
    assert(allNumProcsRecv);

    unsigned int* allKeysSz = new unsigned int[npes];
    assert(allKeysSz);

    unsigned int* allTotalRecv = new unsigned int[npes]; 
    assert(allTotalRecv);

    unsigned int localKeysSize = keys.size(); 
    par::Mpi_Gather<int>(&numProcsSend, allNumProcsSend, 1, 0, comm);
    par::Mpi_Gather<int>(&numProcsRecv, allNumProcsRecv, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&localKeysSize, allKeysSz, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalKeys, allTotalRecv, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < npes; i++) {
        std::cout<<"In Flag Nodes:  allNumProcsSend["<<i<<"] = "<<allNumProcsSend[i]<<std::endl;
        std::cout<<"In Flag Nodes:  allNumProcsRecv["<<i<<"] = "<<allNumProcsRecv[i]<<std::endl;
        std::cout<<"In Flag Nodes:  allKeysSz["<<i<<"] = "<<allKeysSz[i]<<std::endl;
        std::cout<<"In Flag Nodes:  allTotalRecv["<<i<<"] = "<<allTotalRecv[i]<<std::endl;
      }//end for i
    }

    delete [] allNumProcsSend;
    delete [] allNumProcsRecv;
    delete [] allKeysSz;
    delete [] allTotalRecv;
    MPI_Barrier(comm);
    #endif

    // create the send and recv buffers ...
    std::vector<ot::TreeNode> sendK (keys.size());
    std::vector<ot::TreeNode> recvK (totalKeys);
    // the mapping ..
    unsigned int * comm_map = NULL;
    if(!keys.empty()) {
      comm_map = new unsigned int [keys.size()];
      assert(comm_map);
    }

    // Now create sendK
    int *sendOffsets = new int[npes]; 
    assert(sendOffsets);
    sendOffsets[0] = 0;

    int *recvOffsets = new int[npes]; 
    assert(recvOffsets);
    recvOffsets[0] = 0;

    int *numKeysTmp = new int[npes]; 
    assert(numKeysTmp);
    numKeysTmp[0] = 0; 

    // compute offsets ...
    for (int i = 1; i < npes; i++) {
      sendOffsets[i] = sendOffsets[i-1] + numKeysSend[i-1];
      recvOffsets[i] = recvOffsets[i-1] + numKeysRecv[i-1];
      numKeysTmp[i] = 0; 
    }

    for (unsigned int i = 0; i < keys.size(); i++) {
      unsigned int ni = numKeysTmp[part[i]];
      numKeysTmp[part[i]]++;
      // set entry ...
      assert((sendOffsets[part[i]] + ni) < keys.size());
      sendK[sendOffsets[part[i]] + ni] = keys[i];
      // save mapping .. will need it later ...
      comm_map[i] = sendOffsets[part[i]] + ni;
    }
    unsigned int keysSz = keys.size();
    keys.clear();

    if(part) {
      delete [] part;
      part = NULL;
    }

    delete [] numKeysTmp;
    numKeysTmp = NULL;

    ot::TreeNode* sendKptr = NULL;
    ot::TreeNode* recvKptr = NULL;
    if(!sendK.empty()) {
      sendKptr = &(*(sendK.begin()));
    }
    if(!recvK.empty()) {
      recvKptr = &(*(recvK.begin()));
    }

    par::Mpi_Alltoallv_sparse<ot::TreeNode>(sendKptr, numKeysSend, sendOffsets,
        recvKptr, numKeysRecv, recvOffsets, comm);

    sendK.clear();

    PROF_FLN_STAGE3_END
      PROF_FLN_STAGE4_BEGIN

      std::vector<char>  resSend (totalKeys);
    std::vector<char>  resRecv (keysSz);

    for (unsigned int i = 0; i < totalKeys; i++) {
      unsigned int idx;
      resSend[i] = static_cast<char>(
          seq::BinarySearch<ot::TreeNode>(inPtr,
            in.size(),  recvK[i], &idx)) ;   
    }//end for i
    recvK.clear();

    PROF_FLN_STAGE4_END
      PROF_FLN_STAGE5_BEGIN

      char* resSendPtr = NULL;
    char* resRecvPtr = NULL;
    if(!resSend.empty()) {
      resSendPtr = &(*(resSend.begin()));
    }
    if(!resRecv.empty()) {
      resRecvPtr = &(*(resRecv.begin()));
    }

    par::Mpi_Alltoallv_sparse<char>( resSendPtr, numKeysRecv, recvOffsets, 
        resRecvPtr, numKeysSend, sendOffsets, comm);

    resSend.clear();

    delete [] sendOffsets;
    sendOffsets = NULL;

    delete [] recvOffsets;
    recvOffsets = NULL;

    delete [] numKeysSend;
    numKeysSend = NULL;

    delete [] numKeysRecv;
    numKeysRecv = NULL;

    unsigned int st = 0;
    for (unsigned int i = 0; i < in.size(); i++) {
      bool isHanging = false;
      for (unsigned int k = 0; k < keysCount[i]; k++) {
        assert(comm_map[st + k] < resRecv.size());
        if (resRecv[comm_map[st + k]]) {
          isHanging = true;
          break;
        }
      }//end for k
      if (!isHanging) {
        in[i].orFlag(ot::TreeNode::NODE);
      }
      st += keysCount[i];
    }//end for i
    assert(st == keysSz);
    // Clean up ...
    if(comm_map) {
      delete [] comm_map; 
      comm_map = NULL;
    }

    resRecv.clear();

    PROF_FLN_STAGE5_END

      PROF_MARK_HANGING_END
  }//end function


  void flagNodesType2(std::vector<ot::TreeNode> & in, MPI_Comm comm) {

    #ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
    #endif

    PROF_MARK_HANGING_BEGIN

      PROF_FLN_STAGE1_BEGIN

      int npes;
    int rank;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    std::vector<ot::TreeNode > keys;

    assert(!in.empty());
    unsigned int maxD = in[0].getMaxDepth();
    unsigned int dim  = in[0].getDim();

    ot::TreeNode* inPtr = (&(*in.begin()));

    //1. Generate Keys
    //Only true octants need to generate keys. Pseudo-octants don't generate keys.
    //This is because the positive boundary nodes can only be edge-hanging, they
    //can't be face hanging. And the edge hanging node will be identified by the
    //key generated by the true element instead
    for (int i = 0; i < in.size(); i++) {
      unsigned int myLev = inPtr[i].getLevel();
      if (myLev == maxD) {
        continue;
      }
      unsigned int mySz = (1u << (maxD - myLev));
      unsigned int myX = inPtr[i].getX();
      unsigned int myY = inPtr[i].getY();
      unsigned int myZ = inPtr[i].getZ();    

      if( (myX + mySz) <= (1u << (maxD-1)) ) {
        //Add C6 of my +x neighbour (Face Hanging)
        TreeNode tmp6x((myX + mySz), (myY + (mySz>>1)), (myZ + (mySz>>1)), (myLev + 1), dim, maxD);
        keys.push_back(tmp6x);

        //Add C4 of my +x neighbour (Edge Hanging)
        TreeNode tmp4x((myX + mySz), myY, (myZ + (mySz>>1)), (myLev + 1), dim, maxD);
        keys.push_back(tmp4x);

        //Add C2 of my +x neighbour (Edge Hanging)
        TreeNode tmp2x((myX + mySz), (myY + (mySz>>1)), myZ, (myLev + 1), dim, maxD);
        keys.push_back(tmp2x);

        if( (myY + mySz) <= (1u << (maxD-1)) ) {
          //Add C4 of my +xy neighbour (Edge Hanging)
          TreeNode tmp4xy((myX + mySz), (myY + mySz), (myZ + (mySz>>1)), (myLev + 1), dim, maxD);
          keys.push_back(tmp4xy);
        }

        if( (myZ + mySz) <= (1u << (maxD-1)) ) {
          //Add C2 of my +xz neighbour (Edge Hanging)
          TreeNode tmp2xz((myX + mySz), (myY + (mySz>>1)), (myZ + mySz), (myLev + 1), dim, maxD);
          keys.push_back(tmp2xz);
        }
      } 

      if( (myY + mySz) <= (1u << (maxD-1)) ) {
        //Add C5 of my +y neighbour (Face Hanging)
        TreeNode tmp5y((myX + (mySz>>1)), (myY + mySz), (myZ + (mySz>>1)), (myLev + 1), dim, maxD);
        keys.push_back(tmp5y);

        //Add C1 of my +y neighbour (Edge Hanging)
        TreeNode tmp1y((myX + (mySz>>1)), (myY + mySz), myZ, (myLev + 1), dim, maxD);
        keys.push_back(tmp1y);

        //Add C4 of my +y neighbour (Edge Hanging)
        TreeNode tmp4y(myX, (myY + mySz), (myZ + (mySz>>1)), (myLev + 1), dim, maxD);
        keys.push_back(tmp4y);

        if( (myZ + mySz) <= (1u << (maxD-1)) ) {
          //Add C1 of my +yz neighbour (Edge Hanging)
          TreeNode tmp1yz((myX + (mySz>>1)), (myY + mySz), (myZ + mySz), (myLev + 1), dim, maxD);
          keys.push_back(tmp1yz);
        }
      } 

      if( (myZ + mySz) <= (1u << (maxD-1)) ) {
        //Add C3 of my +z neighbour (Face Hanging)
        TreeNode tmp3z((myX + (mySz>>1)), (myY + (mySz>>1)), (myZ + mySz), (myLev + 1), dim, maxD);
        keys.push_back(tmp3z);

        //Add C1 of my +z neighbour (Edge Hanging)
        TreeNode tmp1z((myX + (mySz>>1)), myY, (myZ + mySz), (myLev + 1), dim, maxD);
        keys.push_back(tmp1z);

        //Add C2 of my +z neighbour (Edge Hanging)
        TreeNode tmp2z(myX, (myY + (mySz>>1)), (myZ + mySz), (myLev + 1), dim, maxD);
        keys.push_back(tmp2z);
      }
    }//end for i	

    #ifdef __MEASURE_FLAG_NODES__
    unsigned int forwardKeysCount = 0;
    for (int i = 0; i < in.size(); i++) {
      unsigned int myLev = inPtr[i].getLevel();
      if (myLev == 1) {
        continue;
      }
      unsigned int childNum = inPtr[i].getChildNumber();   
      unsigned int mySz = (1u << (maxD - myLev));
      unsigned int myX = inPtr[i].getX();
      unsigned int myY = inPtr[i].getY();
      unsigned int myZ = inPtr[i].getZ();    

      switch (childNum) {
        case 0:
          {
            break;
          }
        case 7:
          {
            break;
          }
        case 1:
          {
            //-y,-z,-yz
            if (myY) {
              forwardKeysCount++;
            }
            if (myZ) {
              forwardKeysCount++;
            }
            if (myY && myZ) {
              forwardKeysCount++;
            }
            break;
          }
        case 2:
          {
            //-x,-z,-xz
            if (myX) {
              forwardKeysCount++;
            }
            if (myZ) {
              forwardKeysCount++;
            }
            if (myX && myZ) {
              forwardKeysCount++;
            }
            break;
          }
        case 3:
          {
            //-z
            if (myZ) {
              forwardKeysCount++;
            }
            break;
          }
        case 4:
          {
            //-x,-y,-xy
            if (myX) {
              forwardKeysCount++;
            }
            if (myY) {
              forwardKeysCount++;
            }
            if (myX && myY) {
              forwardKeysCount++;
            }
            break;
          }
        case 5:
          {
            //-y
            if (myY) {
              forwardKeysCount++;
            }
            break;
          }
        case 6:
          {
            //-x
            if (myX) {
              forwardKeysCount++;
            }
            break;      
          }
        default: assert(false);
      }//end switch
    }//end for i	
    #endif

    PROF_FLN_STAGE1_END
      PROF_FLN_STAGE2_BEGIN

    #ifdef __MEASURE_FLAG_NODES__
      unsigned int keysSzBefore = keys.size();
    #endif

    //Make keys sorted and unique locally. There could still be duplicates
    //globally and keys need not be globally sorted
    seq::makeVectorUnique<ot::TreeNode>(keys, false);

    #ifdef __MEASURE_FLAG_NODES__
    unsigned int keysSzAfter = keys.size();
    unsigned int* allKeysSzBefore = new unsigned int[npes];
    unsigned int* allKeysSzAfter = new unsigned int[npes];
    unsigned int* allForwardKeysCount = new unsigned int[npes];
    par::Mpi_Gather<unsigned int>(&keysSzBefore, allKeysSzBefore, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&keysSzAfter, allKeysSzAfter, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&forwardKeysCount, allForwardKeysCount, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < npes; i++) {
        std::cout<<"rank = "<<i<<"keys Before = "<<allKeysSzBefore[i]
          <<" After = "<<allKeysSzAfter[i]
          <<" Forward = "<<allForwardKeysCount[i]
          <<std::endl; 
      }//end for i
    }
    delete [] allKeysSzBefore;
    delete [] allKeysSzAfter;
    delete [] allForwardKeysCount;
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE2_END
      PROF_FLN_STAGE3_BEGIN

      // allocate memory for the mins array
      ot::TreeNode* mins = new ot::TreeNode[npes];
    assert(mins);

    par::Mpi_Allgather<ot::TreeNode>(inPtr, mins, 1, comm);
    //Mins will be sorted and unique

    MPI_Request* requests = new MPI_Request[npes<<1];
    assert(requests);
    ot::TreeNode* keysPtr = NULL;
    if(!(keys.empty())) {
      keysPtr = (&(*keys.begin()));
    }

    int* sendOffsets = new int[npes];
    assert(sendOffsets);
    int* sendCounts = new int[npes];
    assert(sendCounts);
    for(int i = 0; i < npes; i++) {
      sendOffsets[i] = 0;
      sendCounts[i] = 0;
    }

    for(int i = 0; i < (npes - 1); i++) {
      for(int j = sendOffsets[i]; j < keys.size(); j++) {
        if(keysPtr[j] >= mins[i+1]) {
          break;
        } else {
          sendCounts[i]++;
        }
      }
      sendOffsets[i + 1] = (sendOffsets[i] + sendCounts[i]);
    }//end for i

    sendCounts[npes - 1] = (keys.size() - sendOffsets[npes - 1]);

    delete [] mins;

    int* recvCounts = new int[npes];
    assert(recvCounts);
    par::Mpi_Alltoall(sendCounts, recvCounts, 1, comm); 

    std::vector<ot::TreeNode>* recvNodes = new std::vector<ot::TreeNode>[npes];
    assert(recvNodes);
    for(int i = 0; i < npes; i++) {
      recvNodes[i].resize(recvCounts[i]);
    }//end for i

    //Post Recvs
    for(int i = 0; i < rank; i++) {
      if(recvCounts[i]) {
        par::Mpi_Irecv<ot::TreeNode>( (&(*((recvNodes[i]).begin()))) , recvCounts[i],
            i, 1, comm, (requests + (i<<1)) );
      }
    }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(recvCounts[i]) {
        par::Mpi_Irecv<ot::TreeNode>( (&(*((recvNodes[i]).begin()))) , recvCounts[i],
            i, 1, comm, (requests + (i<<1)) );
      }
    }//end for i

    //Post Sends
    for(int i = 0; i < rank; i++) {
      if(sendCounts[i]) {
        par::Mpi_Issend<ot::TreeNode>( (keysPtr + sendOffsets[i]), sendCounts[i], 
            i, 1, comm, (requests + ((i<<1) + 1)) );
      }
    }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(sendCounts[i]) {
        par::Mpi_Issend<ot::TreeNode>( (keysPtr + sendOffsets[i]), sendCounts[i], 
            i, 1, comm, (requests + ((i<<1) + 1)) );
      }
    }//end for i

    PROF_FLN_STAGE3_END
      PROF_FLN_STAGE4_BEGIN

      bool* isHanging = new bool[in.size()];
    assert(isHanging);
    for(int i = 0; i < in.size(); i++) {
      isHanging[i] = false; 
    }//end for i

    for(int i = sendOffsets[rank]; i < (sendOffsets[rank] + sendCounts[rank]); i++) {
      unsigned int retIdx;
      bool found =  seq::BinarySearch(inPtr, in.size(), keysPtr[i], &retIdx);
      if(found) {
        //hanging
        isHanging[retIdx] = true;
      }
    }//end for i

    delete [] sendOffsets;

    PROF_FLN_STAGE4_END
      PROF_FLN_STAGE5_BEGIN

      //Wait for the recvs to complete
      for(int i = 0; i < rank; i++) {
        if(recvCounts[i]) {
          MPI_Status status;
          MPI_Wait( (requests + (i<<1)), &status);
        }
      }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(recvCounts[i]) {
        MPI_Status status;
        MPI_Wait( (requests + (i<<1)), &status);
      }
    }//end for i

    PROF_FLN_STAGE5_END
      PROF_FLN_STAGE6_BEGIN

      for(int i = 0; i < rank; i++) {
        for(int j = 0; j < recvCounts[i]; j++) {
          unsigned int retIdx;
          bool found =  seq::BinarySearch(inPtr, in.size(), recvNodes[i][j], &retIdx);
          if(found) {
            //hanging
            isHanging[retIdx] = true;
          }
        }//end for j
      }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      for(int j = 0; j < recvCounts[i]; j++) {
        unsigned int retIdx;
        bool found =  seq::BinarySearch(inPtr, in.size(), recvNodes[i][j], &retIdx);
        if(found) {
          //hanging
          isHanging[retIdx] = true;
        }
      }//end for j
    }//end for i

    for(int i = 0; i < in.size(); i++) {
      if(!isHanging[i]) {
        inPtr[i].orFlag(ot::TreeNode::NODE);
      } 
    }//end for i

    delete [] isHanging;
    delete [] recvNodes;
    delete [] recvCounts;

    PROF_FLN_STAGE6_END
      PROF_FLN_STAGE7_BEGIN

      //Wait for the sends to complete
      for(int i = 0; i < rank; i++) {
        if(sendCounts[i]) {
          MPI_Status status;
          MPI_Wait( (requests + ((i<<1) + 1)), &status);
        }
      }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(sendCounts[i]) {
        MPI_Status status;
        MPI_Wait( (requests + ((i<<1) + 1)), &status);
      }
    }//end for i

    keys.clear();
    delete [] sendCounts;
    delete [] requests;

    PROF_FLN_STAGE7_END

      PROF_MARK_HANGING_END

  }//end function

  void flagNodesType3(std::vector<ot::TreeNode> & in, MPI_Comm comm) {

    #ifdef __PROF_WITH_BARRIER__
    MPI_Barrier(comm);
    #endif

    PROF_MARK_HANGING_BEGIN

    PROF_FLN_STAGE1_BEGIN

    int npes;
    int rank;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    std::vector<ot::TreeNode > keys;

    assert(!in.empty());
    assert(par::test::isUniqueAndSorted(in,MPI_COMM_WORLD));
    unsigned int maxD = in[0].getMaxDepth();
    unsigned int dim  = in[0].getDim();
    ot::TreeNode* inPtr = (&(*(in.begin())));

    //1. Generate Keys
    for (int i = 0; i < in.size(); i++) {
      //keysCount[i] = 0;
      unsigned int myLev = inPtr[i].getLevel();
      if (myLev == 1) {
        continue;
      }

      unsigned int currKeySz = static_cast<unsigned int>(keys.size());

      //unsigned int childNum = inPtr[i].getChildNumber();
      //oda_change_milinda
      unsigned int childNum=inPtr[i].getMortonIndex();

      unsigned int mySz = (1u << (maxD - myLev));
      unsigned int myX = inPtr[i].getX();
      unsigned int myY = inPtr[i].getY();
      unsigned int myZ = inPtr[i].getZ();

     //@hari: I think this is incorrect for Hilbert if we use the getChildNumber function,
      switch (childNum) {
        case 0:
          {
            break;
          }
        case 7:
          {
            break;
          }
        case 1:
          {
            //-y,-z,-yz
            if (myY) {
              TreeNode tmp(myX,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            if (myZ) {
              TreeNode tmp(myX,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            if (myY && myZ) {
              TreeNode tmp(myX,myY-mySz,myZ-mySz,myLev,dim,maxD);         
              keys.push_back(tmp.getParent());
            }
            break;
          }
        case 2:
          {
            //-x,-z,-xz
            if (myX) {
              TreeNode tmp(myX-mySz,myY,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            if (myZ) {
              TreeNode tmp(myX,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            if (myX && myZ) {
              TreeNode tmp(myX-mySz,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            break;
          }
        case 3:
          {
            //-z
            if (myZ) {
              TreeNode tmp(myX,myY,myZ-mySz,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            break;
          }
        case 4:
          {
            //-x,-y,-xy
            if (myX) {
              TreeNode tmp(myX-mySz,myY,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            if (myY) {
              TreeNode tmp(myX,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            if (myX && myY) {
              TreeNode tmp(myX-mySz,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            break;
          }
        case 5:
          {
            //-y
            if (myY) {
              TreeNode tmp(myX,myY-mySz,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());
            }
            break;
          }
        case 6:
          {
            //-x
            if (myX) {
              TreeNode tmp(myX-mySz,myY,myZ,myLev,dim,maxD);
              keys.push_back(tmp.getParent());    
            }
            break;      
          }
        default: assert(false);
      }//end switch
      for(int  j = currKeySz; j < keys.size(); j++) {
        keys[j].setWeight(i);
      }//end for j 
    }//end for i	

    PROF_FLN_STAGE1_END
    //std::sort(keys.begin(),keys.end()); @Hari: Do we need to sort these keys? In my view it is not necessary
    //treeNodesTovtk(keys,rank,"oda_keys");

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 1 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE2_BEGIN

      // allocate memory for the mins array
      std::vector<ot::TreeNode> mins(npes);

    par::Mpi_Allgather<ot::TreeNode>(inPtr, (&(*(mins.begin()))), 1, comm);

    unsigned int *part = NULL;
    if(!keys.empty()) {
      part = new unsigned int[keys.size()];    
      assert(part);
    }

    assert(seq::test::isUniqueAndSorted(mins));

    for (unsigned int i = 0; i < keys.size(); i++) {
      unsigned int idx;
      //maxLB returns the last index in a sorted array
      //such that a[ind] <= key and  a[index +1] > key
      bool found = seq::maxLowerBound<TreeNode >(mins, keys[i], idx, NULL, NULL);
      if (!found ) {
        part[i] = rank;
      } else {
        part[i] = idx;
      }
    }//end for i
    mins.clear();

    PROF_FLN_STAGE2_END 

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 2 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE3_BEGIN

    int *numKeysSend = new int[npes];
    int *numKeysRecv = new int[npes];
    for (int i = 0; i < npes; i++) {
      numKeysSend[i] = 0;
    }

    // calculate the number of keys to send ...
    for (unsigned int i = 0; i < keys.size(); i++) {
      assert(part[i] < npes);
      numKeysSend[part[i]]++;
    }

    // Now do an All2All to get inumKeysRecv
    par::Mpi_Alltoall<int>(numKeysSend, numKeysRecv, 1, comm);

    unsigned int totalKeys=0; // total number of local keys ...
    for (int i = 0; i < npes; i++) {
      totalKeys += numKeysRecv[i];
    }

    #ifdef __MEASURE_FLAG_NODES__
    MPI_Barrier(comm);
    int numProcsSend = 0;
    int numProcsRecv = 0;
    for(int i = 0; i < npes; i++) {
      if(numKeysSend[i]) {
        numProcsSend++;
      }
      if(numKeysRecv[i]) {
        numProcsRecv++;
      }
    }//end for i
    int* allNumProcsSend = new int[npes];
    int* allNumProcsRecv = new int[npes];
    unsigned int* allKeysSz = new unsigned int[npes];
    unsigned int* allTotalRecv = new unsigned int[npes]; 
    unsigned int localKeysSize = keys.size(); 
    par::Mpi_Gather<int>(&numProcsSend, allNumProcsSend, 1, 0, comm);
    par::Mpi_Gather<int>(&numProcsRecv, allNumProcsRecv, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&localKeysSize, allKeysSz, 1, 0, comm);
    par::Mpi_Gather<unsigned int>(&totalKeys, allTotalRecv, 1, 0, comm);
    if(!rank) {
      for(int i = 0; i < npes; i++) {
        std::cout<<"In Flag Nodes:  allNumProcsSend["<<i<<"] = "<<allNumProcsSend[i]<<std::endl;
        std::cout<<"In Flag Nodes:  allNumProcsRecv["<<i<<"] = "<<allNumProcsRecv[i]<<std::endl;
        std::cout<<"In Flag Nodes:  allKeysSz["<<i<<"] = "<<allKeysSz[i]<<std::endl;
        std::cout<<"In Flag Nodes:  allTotalRecv["<<i<<"] = "<<allTotalRecv[i]<<std::endl;
      }//end for i
    }

    delete [] allNumProcsSend;
    delete [] allNumProcsRecv;
    delete [] allKeysSz;
    delete [] allTotalRecv;
    MPI_Barrier(comm);
    #endif

    // create the send and recv buffers ...
    ot::TreeNode* sendK = NULL;
    if(!(keys.empty())) {
      sendK = new ot::TreeNode[keys.size()];
    }

    ot::TreeNode* recvK = NULL;
    if(totalKeys) {
      recvK = new ot::TreeNode[totalKeys];
    }

    // the mapping ..
    unsigned int* comm_map = NULL;
    if(!keys.empty()) {
      comm_map = new unsigned int[keys.size()];
    }

    // Now create sendK
    int *sendOffsets = new int[npes]; sendOffsets[0] = 0;
    int *recvOffsets = new int[npes]; recvOffsets[0] = 0;
    int *numKeysTmp = new int[npes]; numKeysTmp[0] = 0; 

    // compute offsets ...
    for (int i = 1; i < npes; i++) {
      sendOffsets[i] = sendOffsets[i-1] + numKeysSend[i-1];
      recvOffsets[i] = recvOffsets[i-1] + numKeysRecv[i-1];
      numKeysTmp[i] = 0; 
    }

    for (unsigned int i = 0; i < keys.size(); i++) {
      unsigned int ni = numKeysTmp[part[i]];
      numKeysTmp[part[i]]++;
      // set entry ...
      assert((sendOffsets[part[i]] + ni) < keys.size());
      sendK[sendOffsets[part[i]] + ni] = keys[i];
      // save mapping .. will need it later ...
      comm_map[sendOffsets[part[i]] + ni] = keys[i].getWeight();
    }
    unsigned int keysSz = keys.size();
    keys.clear();

    if(part) {
      delete [] part;
      part = NULL;
    }

    delete [] numKeysTmp;
    numKeysTmp = NULL;

    MPI_Request* requests1 = new MPI_Request[npes<<1];

    //Post Recvs
    for(int i = 0; i < rank; i++) {
      if(numKeysRecv[i]) {
        par::Mpi_Irecv<ot::TreeNode>( (recvK + recvOffsets[i]) , numKeysRecv[i],
            i, 1, comm, (requests1 + (i<<1)) );
      }
    }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysRecv[i]) {
        par::Mpi_Irecv<ot::TreeNode>( (recvK + recvOffsets[i]) , numKeysRecv[i],
            i, 1, comm, (requests1 + (i<<1)) );
      }
    }//end for i

    //Post Sends
    for(int i = 0; i < rank; i++) {
      if(numKeysSend[i]) {
        par::Mpi_Issend<ot::TreeNode>( (sendK + sendOffsets[i]), numKeysSend[i], 
            i, 1, comm, (requests1 + ((i<<1) + 1)) );
      }
    }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysSend[i]) {
        par::Mpi_Issend<ot::TreeNode>( (sendK + sendOffsets[i]), numKeysSend[i], 
            i, 1, comm, (requests1 + ((i<<1) + 1)) );
      }
    }//end for i

    PROF_FLN_STAGE3_END

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 3 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE4_BEGIN

      bool* resSend = NULL;
    if(totalKeys) {
      resSend = new bool[totalKeys];
    }

    bool* resRecv = NULL;
    if(keysSz) {
      resRecv = new bool[keysSz];
    }

    assert(par::test::isUniqueAndSorted(in,MPI_COMM_WORLD));

    for (unsigned int i = sendOffsets[rank];
        i < (sendOffsets[rank] + numKeysSend[rank]); i++) {
      unsigned int idx;
      resSend[recvOffsets[rank] + (i - sendOffsets[rank])] =
        seq::BinarySearch<ot::TreeNode>(inPtr, in.size(), sendK[i], &idx);   
    }//end for i

    PROF_FLN_STAGE4_END 

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 4 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE5_BEGIN

      //Wait for the recvs to complete
      for(int i = 0; i < rank; i++) {
        if(numKeysRecv[i]) {
          MPI_Status status;
          MPI_Wait( (requests1 + (i<<1)), &status);
        }
      }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysRecv[i]) {
        MPI_Status status;
        MPI_Wait( (requests1 + (i<<1)), &status);
      }
    }//end for i

    PROF_FLN_STAGE5_END 

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 5 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE6_BEGIN

      for (unsigned int i = 0; i < recvOffsets[rank]; i++) {
        unsigned int idx;
        resSend[i] = seq::BinarySearch<ot::TreeNode>(inPtr, in.size(), recvK[i], &idx);   
      }//end for i

    for (unsigned int i = (recvOffsets[rank] + numKeysRecv[rank]); i < totalKeys; i++) {
      unsigned int idx;
      resSend[i] = seq::BinarySearch<ot::TreeNode>(inPtr, in.size(), recvK[i], &idx);   
    }//end for i

    PROF_FLN_STAGE6_END

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 6 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE7_BEGIN

      MPI_Request* requests2 = new MPI_Request[npes<<1];

    //Post Recvs
    for(int i = 0; i < rank; i++) {
      if(numKeysSend[i]) {
        par::Mpi_Irecv<bool>( (resRecv + sendOffsets[i]) , numKeysSend[i],
            i, 1, comm, (requests2 + (i<<1)) );
      }
    }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysSend[i]) {
        par::Mpi_Irecv<bool>( (resRecv + sendOffsets[i]) , numKeysSend[i],
            i, 1, comm, (requests2 + (i<<1)) );
      }
    }//end for i

    //Post Sends
    for(int i = 0; i < rank; i++) {
      if(numKeysRecv[i]) {
        par::Mpi_Issend<bool>( (resSend + recvOffsets[i]), numKeysRecv[i], 
            i, 1, comm, (requests2 + ((i<<1) + 1)) );
      }
    }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysRecv[i]) {
        par::Mpi_Issend<bool>( (resSend + recvOffsets[i]), numKeysRecv[i], 
            i, 1, comm, (requests2 + ((i<<1) + 1)) );
      }
    }//end for i

    PROF_FLN_STAGE7_END

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 7 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE8_BEGIN
      bool* isHanging = new bool[in.size()];
    for(unsigned int i = 0; i < in.size(); i++) {
      isHanging[i] = false;
    }//end for i

    for(int i = recvOffsets[rank]; i < (recvOffsets[rank] + numKeysRecv[rank]); i++) {
      if(resSend[i]) {
        isHanging[comm_map[sendOffsets[rank] + i - recvOffsets[rank]]] = true;
      }
    }//end for i

    PROF_FLN_STAGE8_END

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 8 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE9_BEGIN

      //wait for recvs
      for(int i = 0; i < rank; i++) {
        if(numKeysSend[i]) {
          MPI_Status status;
          MPI_Wait( (requests2 + (i<<1)), &status);
        }
      }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysSend[i]) {
        MPI_Status status;
        MPI_Wait( (requests2 + (i<<1)), &status);
      }
    }//end for i

    PROF_FLN_STAGE9_END

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 9 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE10_BEGIN

      for(int i = 0; i < sendOffsets[rank]; i++) {
        if(resRecv[i]) {
          isHanging[comm_map[i]] = true;
        }
      }//end for i

    for(int i = (sendOffsets[rank] + numKeysSend[rank]); i < keysSz; i++) {
      if(resRecv[i]) {
        isHanging[comm_map[i]] = true;
      }
    }//end for i

    for(unsigned int i = 0; i < in.size(); i++) {
      if (!isHanging[i]) {
        in[i].orFlag(ot::TreeNode::NODE);
      }
    }//end for i

    PROF_FLN_STAGE10_END

    #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 10 passed."<<std::endl;
    }
    MPI_Barrier(comm);
    #endif

    PROF_FLN_STAGE11_BEGIN

      //Wait for the sends to complete
      for(int i = 0; i < rank; i++) {
        if(numKeysSend[i]) {
          MPI_Status status;
          MPI_Wait( (requests1 + ((i<<1) + 1)), &status);
        }
      }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysSend[i]) {
        MPI_Status status;
        MPI_Wait( (requests1 + ((i<<1) + 1)), &status);
      }
    }//end for i

    for(int i = 0; i < rank; i++) {
      if(numKeysRecv[i]) {
        MPI_Status status;
        MPI_Wait( (requests2 + ((i<<1) + 1)), &status);
      }
    }//end for i

    for(int i = (rank + 1); i < npes; i++) {
      if(numKeysRecv[i]) {
        MPI_Status status;
        MPI_Wait( (requests2 + ((i<<1) + 1)), &status);
      }
    }//end for i

    // Clean up ...
    if(comm_map) {
      delete [] comm_map; 
    }

    if(resRecv) {
      delete [] resRecv;
    }

    if(resSend) {
      delete []resSend;
    }

    if(sendK) {
      delete [] sendK;
    }

    if(recvK) {
      delete [] recvK;
    }

    delete [] sendOffsets;
    delete [] recvOffsets;
    delete [] numKeysSend;
    delete [] numKeysRecv;
    delete [] requests1;
    delete [] requests2;
    delete [] isHanging;

    PROF_FLN_STAGE11_END

  #ifdef __DEBUG_DA_PUBLIC__
      MPI_Barrier(comm);
    if(!rank) {
      std::cout<<"FLN Stage 11 passed."<<std::endl;
    }
    MPI_Barrier(comm);
  #endif

    PROF_MARK_HANGING_END
  }//end function

}//end namespace



