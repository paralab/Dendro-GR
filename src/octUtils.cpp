//
// Created by milinda on 9/6/16.
//
/**
  @brief A collection of simple functions for manipulating octrees.
Examples: Regular Refinements, Linearizing an octree, I/O,
Nearest Common Ancestor, adding positive boundaries, marking hanging nodes
@author Rahul S. Sampath, rahul.sampath@gmail.com
@author Hari Sundar, hsundar@gmail.com
@author Milinda Fernando ,milinda@cs.utah.edu

 @remarks Most of the functions used for the mesh generation. Most of the implementations are based on the previous implementation of dendro version 4.0

*/

#include "octUtils.h"


// This will add boundary nodes and will also embed the octree one level higher
// to enable the addition of the boundary nodes. The positive boundary nodes
// are also marked as BOUNDARY.
void addBoundaryNodesType1(std::vector<ot::TreeNode> &in,
                           std::vector<ot::TreeNode>& bdy,
                           unsigned int dim, unsigned int maxDepth) {

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

}//end function


int refineOctree(const std::vector<ot::TreeNode> & inp,
                 std::vector<ot::TreeNode> &out) {
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


int refineAndPartitionOctree(const std::vector<ot::TreeNode> & inp,
                             std::vector<ot::TreeNode> &out, MPI_Comm comm) {
    refineOctree(inp,out);
    par::partitionW<ot::TreeNode>(out, NULL,comm);
    return 1;
}//end function

int createRegularOctree(std::vector<ot::TreeNode>& out, unsigned int lev,
                        unsigned int dim, unsigned int maxDepth, MPI_Comm comm) {
    ot::TreeNode root(0,0,0,0,dim,maxDepth);
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


int function2Octree(std::function<double(double,double,double)> fx, std::vector<ot::TreeNode> & nodes,unsigned int max_ref_level, const double & tol ,unsigned int elementOrder, MPI_Comm comm )
{


    int size, rank;


    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    nodes.clear();
    std::vector<ot::TreeNode> nodes_new;

    unsigned int depth = 1;
    unsigned int num_intersected=1;
    unsigned int num_intersected_g=1;
    const unsigned int  nodesPerElement=(elementOrder+1)*(elementOrder+1)*(elementOrder+1);

    double h, h1,h2;
    double* dist_parent=new double[nodesPerElement];
    double* dist_child=new double[nodesPerElement];
    double* dist_child_ip=new double[nodesPerElement];
    Point pt;
    Point pt_child;

    const unsigned int maxDepth = m_uiMaxDepth;

    h = 1.0/(1<<(maxDepth));
    unsigned int mySz;
    RefElement refEl(m_uiDim,elementOrder);
    double l2_norm=0;
    bool splitOctant=false;

    if (!rank) {
        // root does the initial refinement
        //std::cout<<"initial ref:"<<std::endl;
        ot::TreeNode root = ot::TreeNode(m_uiDim, maxDepth);
        root.addChildren(nodes);

        while ( (num_intersected > 0 ) && (num_intersected < size/**size*/ ) && (depth < max_ref_level) ) {
            std::cout << "Depth: " << depth << " n = " << nodes.size() << std::endl;
            num_intersected = 0;

            for (auto elem: nodes ){
                splitOctant=false;
                if ( elem.getLevel() != depth ) {
                    nodes_new.push_back(elem);
                    continue;
                }

                mySz=( 1 << (maxDepth - elem.getLevel()));
                h1 = mySz/(double)elementOrder;

                // check and split
                pt = elem.getAnchor();

                for( int k=0;k<(elementOrder+1);k++)
                  for( int j=0;j<(elementOrder+1);j++)
                    for( int i=0;i<(elementOrder+1);i++)
                       dist_parent[k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=fx(pt.x()+i*h1,pt.y()+j*h1,pt.z()+k*h1);


                for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
                {
                    refEl.I3D_Parent2Child(dist_parent,dist_child_ip,cnum);
                    pt_child =Point((pt.x() +(((int)((bool)(cnum & 1u)))<<(m_uiMaxDepth-elem.getLevel()-1))),(pt.y()+(((int)((bool)(cnum & 2u)))<<(m_uiMaxDepth-elem.getLevel()-1))), (pt.z() +(((int)((bool)(cnum & 4u)))<<(m_uiMaxDepth-elem.getLevel()-1))));
                    //std::cout<<"parent: "<<elem<<" child "<<cnum<<" value: "<<pt_child.x()<<" , "<<pt_child.y()<<" , "<<pt_child.z()<<std::endl;
                    h2=h1/2.0;
                    for( int k=0;k<(elementOrder+1);k++)
                        for( int j=0;j<(elementOrder+1);j++)
                            for( int i=0;i<(elementOrder+1);i++)
                                dist_child[k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=fx(pt_child.x()+i*h2,pt_child.y()+j*h2,pt_child.z()+k*h2);

                    l2_norm=normLInfty(dist_child,dist_child_ip,nodesPerElement);//normL2(dist_child,dist_child_ip,nodesPerElement)/normL2(dist_child_ip,nodesPerElement);
                    //std::cout<<"l2: "<<l2_norm<<std::endl;
                    if(l2_norm>tol){
                        splitOctant=true;
                        break;
                    }

                }

               if (!splitOctant) {
                    // if (!skipInternal)
                    nodes_new.push_back(elem);
                }else {
                    // intersection.
                    elem.addChildren(nodes_new);
                    num_intersected++;
                }
            }
            depth++;
            std::swap(nodes, nodes_new);
            nodes_new.clear();
        }
    } // !rank

    // now scatter the elements.
    DendroIntL totalNumOcts = nodes.size(), numOcts;

    par::Mpi_Bcast<DendroIntL>(&totalNumOcts, 1, 0, comm);

    // TODO do proper load balancing.
    numOcts = totalNumOcts/size + (rank < totalNumOcts%size);
    par::scatterValues<ot::TreeNode>(nodes, nodes_new, numOcts, comm);
    std::swap(nodes, nodes_new);
    nodes_new.clear();

    // now refine in parallel.
    par::Mpi_Bcast(&depth, 1, 0, comm);
    num_intersected=1;

    ot::TreeNode root(m_uiDim,m_uiMaxDepth);

    while ( (num_intersected > 0 ) && (depth < max_ref_level) ) {
        if(!rank)std::cout << "Depth: " << depth << " n = " << nodes.size() << std::endl;
        num_intersected = 0;

        for (auto elem: nodes ){
            splitOctant=false;
            if ( elem.getLevel() != depth ) {
                nodes_new.push_back(elem);
                continue;
            }

            mySz=( 1 << (maxDepth - elem.getLevel()));
            h1 = mySz/(double)elementOrder;

            // check and split
            pt = elem.getAnchor();

            for( int k=0;k<(elementOrder+1);k++)
                for( int j=0;j<(elementOrder+1);j++)
                    for( int i=0;i<(elementOrder+1);i++)
                        dist_parent[k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=fx(pt.x()+i*h1,pt.y()+j*h1,pt.z()+k*h1);


            for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
            {
                refEl.I3D_Parent2Child(dist_parent,dist_child_ip,cnum);
                pt_child =Point((pt.x() +(((int)((bool)(cnum & 1u)))<<(m_uiMaxDepth-elem.getLevel()-1))),(pt.y()+(((int)((bool)(cnum & 2u)))<<(m_uiMaxDepth-elem.getLevel()-1))), (pt.z() +(((int)((bool)(cnum & 4u)))<<(m_uiMaxDepth-elem.getLevel()-1))));
                //std::cout<<"parent: "<<elem<<" child "<<cnum<<" value: "<<pt_child.x()<<" , "<<pt_child.y()<<" , "<<pt_child.z()<<std::endl;
                h2=h1/2.0;
                for( int k=0;k<(elementOrder+1);k++)
                    for( int j=0;j<(elementOrder+1);j++)
                        for( int i=0;i<(elementOrder+1);i++)
                            dist_child[k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=fx(pt_child.x()+i*h2,pt_child.y()+j*h2,pt_child.z()+k*h2);

                l2_norm=normLInfty(dist_child,dist_child_ip,nodesPerElement);//normL2(dist_child,dist_child_ip,nodesPerElement)/normL2(dist_child_ip,nodesPerElement);
                //std::cout<<"l2: "<<l2_norm<<std::endl;
                if(l2_norm>tol){
                    splitOctant=true;
                    break;
                }

            }

            if (!splitOctant) {
                // if (!skipInternal)
                nodes_new.push_back(elem);
            }else {
                // intersection.
                elem.addChildren(nodes_new);
                num_intersected++;
            }
        }
        depth++;
        std::swap(nodes, nodes_new);
        nodes_new.clear();

        SFC::parSort::SFC_treeSort(nodes,nodes_new,nodes_new,nodes_new,0.1,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,2,comm);
        std::swap(nodes,nodes_new);
        nodes_new.clear();

        par::Mpi_Allreduce(&num_intersected,&num_intersected_g,1,MPI_MAX,comm);
        num_intersected=num_intersected_g;



    }

    delete[] dist_child;
    delete[] dist_child_ip;
    delete[] dist_parent;
    
    return 0;


}

int function2Octree(std::function<void(double,double,double,double*)> fx,const unsigned int numVars,const unsigned int* varIndex,const unsigned int numInterpVars, std::vector<ot::TreeNode> & nodes,unsigned int max_ref_level, const double & tol ,unsigned int elementOrder,MPI_Comm comm )
{


    int size, rank;


    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    nodes.clear();
    std::vector<ot::TreeNode> nodes_new;

    unsigned int depth = 1;
    unsigned int num_intersected=1;
    unsigned int num_intersected_g=1;
    const unsigned int  nodesPerElement=(elementOrder+1)*(elementOrder+1)*(elementOrder+1);

    double h, h1,h2;
    double* varVal=new double [numVars];
    double* dist_parent=new double[numVars*nodesPerElement];
    double* dist_child=new double[numVars*nodesPerElement];
    double* dist_child_ip=new double[numVars*nodesPerElement];
    Point pt;
    Point pt_child;
    
    const unsigned int maxDepth=m_uiMaxDepth;

    h = 1.0/(1<<(maxDepth));
    unsigned int mySz;
    RefElement refEl(m_uiDim,elementOrder);
    double l2_norm=0;
    bool splitOctant=false;


    if (!rank) {
        // root does the initial refinement
        //std::cout<<"initial ref:"<<std::endl;
        ot::TreeNode root = ot::TreeNode(m_uiDim, maxDepth);
        root.addChildren(nodes);

        while ( (num_intersected > 0 ) && (num_intersected < size/**size*/ ) && (depth < max_ref_level) ) {
            std::cout << "Depth: " << depth << " n = " << nodes.size() << std::endl;
            num_intersected = 0;

            for (auto elem: nodes ){
                splitOctant=false;
                if ( elem.getLevel() != depth ) {
                    nodes_new.push_back(elem);
                    continue;
                }

                mySz=( 1 << (maxDepth - elem.getLevel()));
                h1 = mySz/(double)elementOrder;
                h2=h1/2.0;
                // check and split
                pt = elem.getAnchor();

                for( int k=0;k<(elementOrder+1);k++)
                    for( int j=0;j<(elementOrder+1);j++)
                        for( int i=0;i<(elementOrder+1);i++)
                        {
                            fx(pt.x()+i*h1,pt.y()+j*h1,pt.z()+k*h1,varVal);
                            for(unsigned int var=0;var<numInterpVars;var++)
                            dist_parent[varIndex[var]*nodesPerElement+k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=varVal[varIndex[var]];
                        }



                for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
                {


                    pt_child =Point((pt.x() +(((int)((bool)(cnum & 1u)))<<(m_uiMaxDepth-elem.getLevel()-1))),(pt.y()+(((int)((bool)(cnum & 2u)))<<(m_uiMaxDepth-elem.getLevel()-1))), (pt.z() +(((int)((bool)(cnum & 4u)))<<(m_uiMaxDepth-elem.getLevel()-1))));
                    //std::cout<<"parent: "<<elem<<" child "<<cnum<<" value: "<<pt_child.x()<<" , "<<pt_child.y()<<" , "<<pt_child.z()<<std::endl;
                    for( int k=0;k<(elementOrder+1);k++)
                        for( int j=0;j<(elementOrder+1);j++)
                            for( int i=0;i<(elementOrder+1);i++)
                            {
                                fx(pt_child.x()+i*h2,pt_child.y()+j*h2,pt_child.z()+k*h2,varVal);
                                for(unsigned int var=0;var<numInterpVars;var++)
                                    dist_child[varIndex[var]*nodesPerElement+k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=varVal[varIndex[var]];
                            }


                    for(unsigned int var=0;var<numInterpVars;var++)
                    {
                        refEl.I3D_Parent2Child(dist_parent+varIndex[var]*nodesPerElement,dist_child_ip+varIndex[var]*nodesPerElement,cnum);
                        l2_norm=normLInfty(dist_child+varIndex[var]*nodesPerElement,dist_child_ip+varIndex[var]*nodesPerElement,nodesPerElement);
                        if(l2_norm>tol)
                        {
                            splitOctant=true;
                            break;
                        }

                    }

                    if(splitOctant) break;

                }

                if (!splitOctant) {
                    nodes_new.push_back(elem);
                }else {
                    elem.addChildren(nodes_new);
                    num_intersected++;
                }
            }
            depth++;
            std::swap(nodes, nodes_new);
            nodes_new.clear();
        }
    } // !rank

    // now scatter the elements.
    DendroIntL totalNumOcts = nodes.size(), numOcts;

    par::Mpi_Bcast<DendroIntL>(&totalNumOcts, 1, 0, comm);

    // TODO do proper load balancing.
    numOcts = totalNumOcts/size + (rank < totalNumOcts%size);
    par::scatterValues<ot::TreeNode>(nodes, nodes_new, numOcts, comm);
    std::swap(nodes, nodes_new);
    nodes_new.clear();


    // now refine in parallel.
    par::Mpi_Bcast(&depth, 1, 0, comm);
    num_intersected=1;

    ot::TreeNode root(m_uiDim,m_uiMaxDepth);

    while ( (num_intersected > 0 ) && (depth < max_ref_level) ) {
        if(!rank)std::cout << "Depth: " << depth << " n = " << nodes.size() << std::endl;
        num_intersected = 0;

        for (auto elem: nodes ){
            splitOctant=false;
            if ( elem.getLevel() != depth ) {
                nodes_new.push_back(elem);
                continue;
            }

            mySz=( 1 << (maxDepth - elem.getLevel()));
            h1 = mySz/(double)elementOrder;
            h2=h1/2.0;

            // check and split
            pt = elem.getAnchor();

            for( int k=0;k<(elementOrder+1);k++)
                for( int j=0;j<(elementOrder+1);j++)
                    for( int i=0;i<(elementOrder+1);i++)
                    {
                        fx(pt.x()+i*h1,pt.y()+j*h1,pt.z()+k*h1,varVal);
                        for(unsigned int var=0;var<numInterpVars;var++)
                            dist_parent[varIndex[var]*nodesPerElement+k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=varVal[varIndex[var]];
                    }


            for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
            {
                pt_child =Point((pt.x() +(((int)((bool)(cnum & 1u)))<<(m_uiMaxDepth-elem.getLevel()-1))),(pt.y()+(((int)((bool)(cnum & 2u)))<<(m_uiMaxDepth-elem.getLevel()-1))), (pt.z() +(((int)((bool)(cnum & 4u)))<<(m_uiMaxDepth-elem.getLevel()-1))));
                //std::cout<<"parent: "<<elem<<" child "<<cnum<<" value: "<<pt_child.x()<<" , "<<pt_child.y()<<" , "<<pt_child.z()<<std::endl;
                for( int k=0;k<(elementOrder+1);k++)
                    for( int j=0;j<(elementOrder+1);j++)
                        for( int i=0;i<(elementOrder+1);i++)
                        {
                            fx(pt_child.x()+i*h2,pt_child.y()+j*h2,pt_child.z()+k*h2,varVal);
                            for(unsigned int var=0;var<numInterpVars;var++)
                                dist_child[varIndex[var]*nodesPerElement+k*(elementOrder+1)*(elementOrder+1)+j*(elementOrder+1)+i]=varVal[varIndex[var]];
                        }


                for(unsigned int var=0;var<numInterpVars;var++)
                {
                    refEl.I3D_Parent2Child(dist_parent+varIndex[var]*nodesPerElement,dist_child_ip+varIndex[var]*nodesPerElement,cnum);
                    l2_norm=normLInfty(dist_child+varIndex[var]*nodesPerElement,dist_child_ip+varIndex[var]*nodesPerElement,nodesPerElement);
                    //std::cout<<"rank: "<<rank<<" node: "<<elem<<" l2 norm : "<<l2_norm<<" var: "<<varIndex[var]<<std::endl;
                    if(l2_norm>tol)
                    {
                        splitOctant=true;
                        break;
                    }

                }

                if(splitOctant) break;

            }

            if (!splitOctant) {
                nodes_new.push_back(elem);
            }else {
                elem.addChildren(nodes_new);
                num_intersected++;
            }
        }
        depth++;
        std::swap(nodes, nodes_new);
        nodes_new.clear();
        //par::partitionW<ot::TreeNode>(nodes,NULL,comm);

        SFC::parSort::SFC_treeSort(nodes,nodes_new,nodes_new,nodes_new,0.1,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,2,comm);
        std::swap(nodes,nodes_new);
        nodes_new.clear();

        par::Mpi_Allreduce(&num_intersected,&num_intersected_g,1,MPI_MAX,comm);
        num_intersected=num_intersected_g;


    }

    delete[] dist_child;
    delete[] dist_child_ip;
    delete[] dist_parent;
    delete[] varVal;

    return 0;
    

}

void octree2BlockDecomposition(std::vector<ot::TreeNode>& pNodes, std::vector<ot::Block>& blockList,unsigned int maxDepth,unsigned int & d_min, unsigned int & d_max,DendroIntL localBegin, DendroIntL localEnd,unsigned int eleOrder, unsigned int coarsetLev, unsigned int* tag, unsigned int tsz)
{

    // Note that we assume pnodes to be sorted.
    assert(seq::test::isUniqueAndSorted(pNodes));

    // Note: Commented out code is for debugging purposes.
    #ifdef OCT2BLK_DEBUG
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    treeNodesTovtk(pNodes,rank,"balOct");
    #endif

    unsigned int x,y,z,hindex,hindexN,index;
    unsigned int rot_id=ROOT_ROTATION;

    std::vector<ot::Block> initialBLocks;
    ot::TreeNode rootNode(0,0,0,0,m_uiDim,maxDepth);

    d_min=maxDepth;
    d_max=0;
    // Computes the dmin and dmax of the tree.
    for(unsigned int k=0;k<pNodes.size();k++)
    {
        if(d_min>pNodes[k].getLevel())
            d_min=pNodes[k].getLevel();

        if(d_max < pNodes[k].getLevel())
            d_max=pNodes[k].getLevel();
    }

    DendroIntL nBegin,nEnd;

    ot::Block rootBlock(rootNode,ROOT_ROTATION,d_min,localBegin,localEnd,eleOrder);
    initialBLocks.push_back(rootBlock);

    unsigned int currRegGridLev=d_min;
    ot::TreeNode parent;
    ot::Block tmpBlock;


    DendroIntL splitters[NUM_CHILDREN+1];
    unsigned int childHasRegLev[NUM_CHILDREN]; // 0 if child i does not have any octants at the reg grid lev , and 1 otherwise.
    //unsigned int numChildHasRegLev=0;
    unsigned int  numRegGridOcts=0; // total number of children that has given reg grid levels.
    unsigned int pMaxDepthBit=0;

    //unsigned  int debug_rank=1;

    DendroIntL numIdealRegGridOct=0;
    double blockFillRatio=0.0; // ratio between number of octants in ideal regular grid and actually available.

    DendroUInt_128 blockVolume=0; // 128-bit integer to store the volume of the block. .
    DendroUInt_128 octVolume=0; // 128-bit integer to store oct volume inside a block.

    

    while(!initialBLocks.empty())
    {

        tmpBlock=initialBLocks.back();
        initialBLocks.pop_back();

        parent=tmpBlock.getBlockNode();

        currRegGridLev=tmpBlock.getRegularGridLev();
        rot_id=tmpBlock.getRotationID();
        nBegin=tmpBlock.getLocalElementBegin();
        nEnd=tmpBlock.getLocalElementEnd();
        assert(parent.getLevel()<=currRegGridLev);
        bool octLevelGap=true;
        bool isTagValid=true;
        if(parent.getLevel()==currRegGridLev)
        {
            assert((nEnd-nBegin)==1);
            assert(pNodes[nBegin]==parent);
            blockList.push_back(tmpBlock);
            continue;
        }

        numRegGridOcts=0;
        numIdealRegGridOct=(1u<<(currRegGridLev-parent.getLevel()));
        blockVolume=1u<<((maxDepth-parent.getLevel())*3);
        (m_uiDim==3)? numIdealRegGridOct=numIdealRegGridOct*numIdealRegGridOct*numIdealRegGridOct : numIdealRegGridOct=numIdealRegGridOct*numIdealRegGridOct;
        octVolume=0;

        if(tag != NULL)
        {
            assert(tsz == (localEnd-localBegin));
            for(unsigned int elem=nBegin;elem<nEnd;elem++)
            {
                if(pNodes[elem].getLevel()==currRegGridLev)
                    numRegGridOcts++;
                else if(abs((int)pNodes[elem].getLevel()-(int)currRegGridLev)>OCT2BLK_DECOMP_LEV_GAP){
                    octLevelGap=false;
                    break;
                }

                if(tag[nBegin-localBegin] != tag[elem-localBegin])
                    isTagValid = false;

                octVolume+=1u<<(3*(maxDepth-pNodes[elem].getLevel()));

            }

        }else
        {
            for(unsigned int elem=nBegin;elem<nEnd;elem++)
            {
                if(pNodes[elem].getLevel()==currRegGridLev)
                    numRegGridOcts++;
                else if(abs((int)pNodes[elem].getLevel()-(int)currRegGridLev)>OCT2BLK_DECOMP_LEV_GAP){
                    octLevelGap=false;
                    break;
                }

                octVolume+=1u<<(3*(maxDepth-pNodes[elem].getLevel()));
            }

        }

        //if(rank==debug_rank) std::cout<<"current parent: "<<parent<<" rot id: "<<(int)rot_id<<" curRegGridLev: "<<currRegGridLev<<" begin: "<<nBegin<<" end: "<<nEnd<<" numRegGridOcts: "<<numRegGridOcts<<std::endl;
        blockFillRatio=(double) numRegGridOcts/numIdealRegGridOct;
        if((parent.getLevel()>=coarsetLev) && (isTagValid) && (octLevelGap) && (blockFillRatio>=OCT2BLK_DECOMP_BLK_FILL_RATIO) && (octVolume==blockVolume))
        {
           blockList.push_back(tmpBlock);
           if((currRegGridLev+1)<=d_max) initialBLocks.push_back(ot::Block(parent,rot_id,(currRegGridLev+1),nBegin,nEnd,eleOrder));

        }else
        { // implies that we need to split the tmpBlock.
            assert(parent.getLevel()<maxDepth);
            pMaxDepthBit=maxDepth-parent.getLevel()-1;
            SFC::seqSort::SFC_bucketing(&(*(pNodes.begin())),parent.getLevel(),maxDepth,rot_id,nBegin,nEnd,splitters);

            for (int i = 0; i < NUM_CHILDREN; i++) {
                childHasRegLev[i]=0;
                hindex = (rotations[2 * NUM_CHILDREN * rot_id + i] - '0');
                if (i == (NUM_CHILDREN-1))
                    hindexN = i + 1;
                else
                    hindexN = (rotations[2 * NUM_CHILDREN * rot_id + i + 1] - '0');
                assert(splitters[hindex] <= splitters[hindexN]);

                for(unsigned int elem=splitters[hindex];elem<splitters[hindexN];elem++)
                {
                    if(pNodes[elem].getLevel()==currRegGridLev)
                    {
                        childHasRegLev[i]=1;
                        break;
                    }
                }

            }

            for(unsigned int i=0;i<(NUM_CHILDREN);i++)
            {
                hindex = (rotations[2 * NUM_CHILDREN * rot_id + i] - '0');
                if (i == (NUM_CHILDREN-1))
                    hindexN = i + 1;
                else
                    hindexN = (rotations[2 * NUM_CHILDREN * rot_id + i + 1] - '0');
                assert(splitters[hindex] <= splitters[hindexN]);
                index = HILBERT_TABLE[NUM_CHILDREN * rot_id + hindex];


                x=parent.getX() +(((int)((bool)(hindex & 1u)))<<(pMaxDepthBit));
                y=parent.getY() +(((int)((bool)(hindex & 2u)))<<(pMaxDepthBit));
                z=parent.getZ() +(((int)((bool)(hindex & 4u)))<<(pMaxDepthBit));

                if((childHasRegLev[i]==1))
                {
                    if((parent.getLevel()+1)<=currRegGridLev)
                    {
                        tmpBlock=ot::Block(ot::TreeNode(x,y,z,parent.getLevel()+1,m_uiDim,maxDepth),index,currRegGridLev,splitters[hindex],splitters[hindexN],eleOrder);
                        initialBLocks.push_back(tmpBlock);
                    }
                }else if(((childHasRegLev[i]==0 && (splitters[hindex]!= splitters[hindexN]))))
                {

                    if((currRegGridLev+1)<=d_max && ((parent.getLevel()+1) <=(currRegGridLev+1)))
                    {
                        tmpBlock=ot::Block(ot::TreeNode(x,y,z,parent.getLevel()+1,m_uiDim,maxDepth),index,(currRegGridLev+1),splitters[hindex],splitters[hindexN],eleOrder);
                        //if(rank==debug_rank) std::cout<<"block node pushed (initialBlocks): "<<tmpBlock.getBlockNode()<<std::endl;
                        initialBLocks.push_back(tmpBlock);
                    }

                }

            }

        }
    }

   std::reverse(blockList.begin(),blockList.end());

    #ifdef OCT2BLK_DEBUG
    std::vector<ot::TreeNode> blockNodes;
    blockNodes.resize(blockList.size());
    for(unsigned int k=0;k<blockList.size();k++)
    {
        blockNodes[k]=blockList[k].getBlockNode();
    }

    treeNodesTovtk(blockNodes,rank,"blockNodes");
    #endif

    #ifdef OCT2BLK_DEBUG
    unsigned int numIdealRegOcts=0;
    unsigned int numActualRegOcts=0;
    unsigned int singleOctBlockCount=0;
    for(unsigned int i=0;i<blockList.size();i++)
    {

        numIdealRegOcts=1u<<(blockList[i].getRegularGridLev()-blockList[i].getBlockNode().getLevel());
        numIdealRegOcts=numIdealRegOcts*numIdealRegOcts*numIdealRegOcts;

        numActualRegOcts=0;

        for(unsigned int j=blockList[i].getLocalElementBegin();j<(blockList[i].getLocalElementEnd());j++)
        {
            if(pNodes[j].getLevel()==blockList[i].getRegularGridLev())
                numActualRegOcts++;
        }

       if(numActualRegOcts==1)
            singleOctBlockCount++;

       //std::cout<<"rank: "<<rank<<" block ID: "<<i<<" : "<<blockList[i].getBlockNode()<<" reg lev: "<<blockList[i].getRegularGridLev()<<" ideal reg: "<<numIdealRegOcts<<" actual reg oct: "<<numActualRegOcts<<" ratio: "<<((double)numActualRegOcts/numIdealRegOcts)<<std::endl;

    }
    std::cout<<"rank: "<<rank<<" singleOctBlocks: "<<singleOctBlockCount<<" pNodes size: "<<pNodes.size()<<" ratio: "<<(double(singleOctBlockCount)/pNodes.size())<<std::endl;
    #endif



}


void blockListToVtk(std::vector<ot::Block>& blkList, const std::vector<ot::TreeNode>& pNodes,char* fNamePrefix, MPI_Comm comm)
{

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    std::vector<ot::TreeNode> blockNodeList;
    std::vector<ot::TreeNode> octsEmbedded;

    char fNameBlk[256];
    char fNameEmbed[256];

    for(unsigned int k=0;k<blkList.size();k++)
    {
        octsEmbedded.clear();
        blockNodeList.clear();
        blockNodeList.push_back(blkList[k].getBlockNode());

        for(unsigned int j=blkList[k].getLocalElementBegin();j<blkList[k].getLocalElementEnd();j++)
            octsEmbedded.push_back(pNodes[j]);

        sprintf(fNameEmbed,"%s_embed_%d",fNamePrefix,k);
        sprintf(fNameBlk,"%s_blk_%d",fNamePrefix,k);

        if(rank==0)treeNodesTovtk(blockNodeList,rank,fNameBlk);
        if(rank==0)treeNodesTovtk(octsEmbedded,rank,fNameEmbed);


    }


}


void enforceSiblingsAreNotPartitioned(std::vector<ot::TreeNode> & in,MPI_Comm comm)
{

    // handles the empty octree case.
    MPI_Comm newComm;
    par::splitComm2way(in.empty(),&newComm,comm);

    if(!in.empty())
    {
        comm=newComm;
        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int prev,next;
        (rank<(npes-1)) ? next=rank+1: next=0;
        (rank>0) ? prev=rank-1: prev=(npes-1);

        //if(in.size()<NUM_CHILDREN) return ;

        unsigned int blCount=0;
        int sendCount=0;
        int recvCount=0;

        MPI_Request req;
        MPI_Status status;

        ot::TreeNode maxNode=in.back();
        ot::TreeNode prev_max;
        if(rank<(npes-1)) sendCount=1;
        if(rank>0) recvCount=1;

        par::Mpi_Sendrecv(&maxNode,sendCount,next,0,&prev_max,recvCount,prev,0,comm,&status);

        sendCount=0;
        recvCount=0;

        for(unsigned int elem=0;elem<NUM_CHILDREN;elem++) {
            if ((elem < in.size()) && (in.front().getParent() == in[elem].getParent()))
                sendCount++;
            else
                break;
        }

        if((rank>0) && (sendCount==1) && (in.front().getParent()!=prev_max.getParent())) sendCount=0;
        if(sendCount==NUM_CHILDREN || rank==0) sendCount=0;

        //std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" sendCount: "<<sendCount<<std::endl;
        par::Mpi_Sendrecv(&sendCount,1,prev,1,&recvCount,1,next,1,comm,&status);
        /*
        MPI_Isend(&sendCount,1,MPI_INT,prev,0,m_uiCommActive,&req);
        MPI_Recv(&recvCount,1,MPI_INT,next,0,m_uiCommActive,&status);*/

        //std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" recvCount: "<<recvCount<<std::endl;

        assert(in.size()>=sendCount);
        std::vector<ot::TreeNode> recvBuffer;
        recvBuffer.resize(recvCount);
        par::Mpi_Sendrecv(&(*(in.begin())),sendCount,prev,2,&(*(recvBuffer.begin())),recvCount,next,2,comm,&status);
        //std::cout<<"rank: "<<m_uiActiveRank<<" send recv ended "<<std::endl;

        if(sendCount)
        {
            std::vector<ot::TreeNode> tmpElements;
            std::swap(in,tmpElements);
            in.clear();
            in.resize(tmpElements.size()-sendCount);
            for(unsigned int ele=sendCount;ele<tmpElements.size();ele++)
                in[ele-sendCount]=tmpElements[ele];

            tmpElements.clear();
        }

        for(unsigned int ele=0;ele<recvBuffer.size();ele++)
            in.push_back(recvBuffer[ele]);

        recvBuffer.clear();

        assert(seq::test::isUniqueAndSorted(in));

    }

    MPI_Comm_free(&newComm);
    return ;

}


void mergeKeys(std::vector<ot::SearchKey>& sKeys,std::vector<ot::Key> & keys)
{

    if(sKeys.size()==0) return;

    ot::SearchKey rootSkey(m_uiDim,m_uiMaxDepth);
    std::vector<ot::SearchKey> tmpSKeys;
    SFC::seqSort::SFC_treeSort(&(*(sKeys.begin())),sKeys.size(),tmpSKeys,tmpSKeys,tmpSKeys,m_uiMaxDepth,m_uiMaxDepth,rootSkey,ROOT_ROTATION,1,TS_SORT_ONLY);
    assert(seq::test::isSorted(sKeys));

    ot::Key tmpKey;
    unsigned int skip=0;
    const unsigned int K=1;
    for(unsigned int e=0;e<(sKeys.size());e++)
    {
        tmpKey=ot::Key(sKeys[e].getX(),sKeys[e].getY(),sKeys[e].getZ(),sKeys[e].getLevel(),m_uiDim,m_uiMaxDepth);
        if(sKeys[e].getOwner()>=0){
            tmpKey.addOwner(sKeys[e].getOwner());
            tmpKey.addStencilIndexAndDirection(K-1,sKeys[e].getStencilIndexDirectionList());
        }

        skip=1;
        while(((e+skip)<sKeys.size()) && (sKeys[e]==sKeys[e+skip]))
        {
            if(sKeys[e+skip].getOwner()>=0){
                tmpKey.addOwner(sKeys[e+skip].getOwner());
                tmpKey.addStencilIndexAndDirection(K-1,sKeys[e+skip].getStencilIndexDirectionList());
            }
            skip++;
        }

        keys.push_back(tmpKey);
        e+=(skip-1);

    }


}


void generateBlkEdgeSKeys(const ot::Block & blk, std::vector<ot::SearchKey>& sKeys)
{
    const unsigned int domain_max = 1u<<(m_uiMaxDepth);
    const ot::TreeNode blkNode=blk.getBlockNode();
    const unsigned int regLev=blk.getRegularGridLev();

    const unsigned int blkElem_1D=(1u<<(regLev-blkNode.getLevel()))*2;

    const unsigned int myX=blkNode.getX();
    const unsigned int myY=blkNode.getY();
    const unsigned int myZ=blkNode.getZ();
    const unsigned int mySz=1u<<(m_uiMaxDepth-blkNode.getLevel());
    const unsigned int hsz=1u<<(m_uiMaxDepth-regLev-1); //hx/2
    std::vector<ot::SearchKey>::iterator hint;

    if(myX>0 && myY>0)
    {
        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY-1),(myZ+k*hsz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_LEFT_DOWN);
        }

    }

    if(myX>0 && (myY+mySz)<domain_max) {

        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint = sKeys.emplace(sKeys.end(),ot::SearchKey((myX - 1), (myY + mySz), (myZ+k*hsz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_LEFT_UP);
        }


    }

    if(myX>0 && myZ>0) {

        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint = sKeys.emplace(sKeys.end(), ot::SearchKey((myX - 1), (myY+k*hsz), (myZ - 1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_LEFT_BACK);
        }

    }

    if(myX>0 && (myZ+mySz)<domain_max)
    {
        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY+k*hsz),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_LEFT_FRONT);
        }


    }


    if((myX+mySz) < domain_max && myY>0)
    {
        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY-1),(myZ+k*hsz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_DOWN);
        }

    }

    if((myX+mySz)<domain_max && (myY+mySz)<domain_max)
    {
        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY+mySz),(myZ+k*hsz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_UP);
        }


    }


    if((myX+mySz)<domain_max && myZ>0) {

        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint = sKeys.emplace(sKeys.end(), ot::SearchKey((myX + mySz), (myY+k*hsz), (myZ - 1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_BACK);
        }


    }

    if((myX+mySz)<domain_max && (myZ+mySz)<domain_max)
    {
        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY+k*hsz),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_FRONT);
        }



    }

    if(myY>0 && myZ>0)
    {
        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+k*hsz),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_DOWN_BACK);
        }

    }

    if(myY > 0 && (myZ+mySz)<domain_max)
    {

        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+k*hsz),(myY-1),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_DOWN_FRONT);
        }



    }


    if((myY+mySz)<domain_max && myZ>0)
    {
        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+k*hsz),(myY+mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_UP_BACK);
        }


    }

    if((myY+mySz)<domain_max && (myZ+mySz)<domain_max) {


        for(unsigned int k=0;k<blkElem_1D;k++)
        {
            hint = sKeys.emplace(sKeys.end(), ot::SearchKey((myX+k*hsz), (myY + mySz), (myZ + mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
            hint->addOwner(k);
            hint->addStencilIndexAndDirection(OCT_DIR_UP_FRONT);
        }


    }

}


void generateBlkVertexSKeys(const ot::Block & blk, std::vector<ot::SearchKey>& sKeys)
{
    const unsigned int domain_max = 1u<<(m_uiMaxDepth);
    const ot::TreeNode blkNode=blk.getBlockNode();
    const unsigned int regLev=blk.getRegularGridLev();

    const unsigned int blkElem_1D=(1u<<(regLev-blkNode.getLevel()))*2;

    const unsigned int myX=blkNode.getX();
    const unsigned int myY=blkNode.getY();
    const unsigned int myZ=blkNode.getZ();
    const unsigned int mySz=1u<<(m_uiMaxDepth-blkNode.getLevel());

    std::vector<ot::SearchKey>::iterator hint;


    if((myX>0) && (myY>0) && (myZ>0))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(0);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_DOWN_BACK);
    }

    if(((myX+mySz)<domain_max) && (myY>0) && (myZ>0))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY-1),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(1);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_DOWN_BACK);
    }

    if((myX>0) && ((myY+mySz)<domain_max) && (myZ>0))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY+mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(2);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_UP_BACK);
    }


    if(((myX+mySz)<domain_max) && ((myY+mySz)<domain_max) && (myZ>0))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY+mySz),(myZ-1), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(3);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_UP_BACK);
    }

    if((myX>0) && (myY>0) && ((myZ+mySz)<domain_max))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY-1),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(4);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_DOWN_FRONT);
    }

    if(((myX+mySz)<domain_max) && (myY>0) && ((myZ+mySz)<domain_max))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY-1),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(5);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_DOWN_FRONT);
    }

    if((myX>0) && ((myY+mySz)<domain_max) && ((myZ+mySz)<domain_max))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX-1),(myY+mySz),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(6);
        hint->addStencilIndexAndDirection(OCT_DIR_LEFT_UP_FRONT);
    }


    if(((myX+mySz)<domain_max) && ((myY+mySz)<domain_max) && ((myZ+mySz)<domain_max))
    {
        hint=sKeys.emplace(sKeys.end(),ot::SearchKey((myX+mySz),(myY+mySz),(myZ+mySz), m_uiMaxDepth, m_uiDim, m_uiMaxDepth));
        hint->addOwner(7);
        hint->addStencilIndexAndDirection(OCT_DIR_RIGHT_UP_FRONT);
    }

}


unsigned int rankSelectRule(unsigned int size_global,unsigned int rank_global, unsigned int size_local,unsigned int rank_i)
{
    
    
    if(size_local>size_global){std::cout<<"[Error] : "<<__func__<<" rank: "<<rank_global<<" size_local > size_global "<<std::endl; exit(0);}
    
    // Rule 1 [enable this code to choose consecative ranks (deafult) works with any size_global and size_local]
    #ifdef BSSN_CONSEC_COMM_SELECT
        return rank_i;
    #else
    // Rule 2 [select ranks which is equivalent to complete binary tree fashion (size_global) & (rank_global) needs to be power of two.]
    if((!binOp::isPowerOfTwo(size_global)) || (!binOp::isPowerOfTwo(size_local)))
        return rank_i;
    else
    {
        const unsigned int commMaxLevel=binOp::fastLog2(size_global);
        const unsigned int commLevel=binOp::fastLog2(size_local);
        const unsigned int commK=1u<<(commMaxLevel-commLevel);
        
        return rank_i*commK;
        
        
    }
    #endif   
    
    
}

