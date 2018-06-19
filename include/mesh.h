//
// Created by milinda on 9/2/16.

/**
 * @author: Milinda Shayamal Fernando
 *
 * School of Computing , University of Utah
 *
 * @breif Contains the functions to generate the mesh data structure from the 2:1 balanced linear octree.
 * @details
 * Assumptions:
 * 1). Assumes that octree is balanced and sorted.
 * 2). Assumes that there is no duplicate nodes.
 * 3). Inorder to make the implementation much simpler we use the Morton ordering to number the each vertex of a given particular element. Note that, this doesn't cause any dependency between SFC curve used and the mesh
 * generation part.
 *
 * */


//

#ifndef SFCSORTBENCH_MESH_H
#define SFCSORTBENCH_MESH_H


#include <vector>
#include "mpi.h"
#include "octUtils.h"
#include "treenode2vtk.h"
#include "TreeNode.h"
#include "sfcSort.h"
#include "sfcSearch.h"
#include "testUtils.h"

#include <memory>
#include <fdCoefficient.h>
#include "refel.h"
#include "block.h"
#include "stencil.h"
#include "key.h"
#include "skey.h"
#include "node.h"
#include "dendro.h"

#include "wavelet.h"
#include "dendroProfileParams.h" // only need to profile unzip_asyn for bssn. remove this header file later.



extern double t_e2e; // e2e map generation time
extern double t_e2n; // e2n map generation time
extern double t_sm; // sm map generation time
extern double t_blk;

extern double t_e2e_g[3];
extern double t_e2n_g[3];
extern double t_sm_g[3];
extern double t_blk_g[3];

#include "dendroProfileParams.h"



/**
 * How the oct flags are used in the mesh generation part.
 *
 * bit 0 - 4 : Used to store the level.
 * bit 5: 1- if the key was found in the tree 0- otherwise.
 * bit 6: Key type NONE
 * bit 7: Key type SPLITTER
 * bit 8: Key type UP
 * bit 9: Key type DOWN
 * bit 10: Key type FRONT
 * bit 11: Key type BACK
 * bit 12: Key type LEFT
 * bit 13: key type RIGHT
 *
 * */

enum KeyType{NONE,SPLITTER,UP,DOWN,FRONT,BACK,LEFT,RIGHT};




//#define DEBUG_MESH_GENERATION


#define KS_MAX 31 // maximum number of points in any direction that wavelet DA supports.
#define KEY_DIR_OFFSET 7
#define CHAINED_GHOST_OFFSET 5u

#define OCT_NO_CHANGE 0u
#define OCT_SPLIT 1u
#define OCT_COARSE 2u

namespace ot
{

    namespace WaveletDA
    {
        enum LoopType{ALL,INDEPENDENT,DEPENDENT};
    }


}




namespace ot
{
    /** Structure to order the send recv nodes based on the (e,i,j,k ) ordering. */
    struct NodeTuple
    {
        private:
             /** i index of the node*/
             unsigned int i;
            /** j index of the node*/
            unsigned int j;
            /** k index of the node*/
            unsigned int k;
            /** element id of the node*/
            unsigned int e;

    public:
          NodeTuple(unsigned int xi,unsigned int yj,unsigned int zk,unsigned int owner)
          {
              i=xi;
              j=yj;
              k=zk;
              e=owner;
          }
          inline unsigned int getX() const {return i;}
          inline unsigned int getY() const {return j;}
          inline unsigned int getZ() const {return k;}
          inline unsigned int getLevel()const {return e;}
    };


}




namespace ot {

    /**
     * @breif Contains all the information needed to build neighbourhood information for the balanced octree based on treeSearch.
     * */


 class Mesh {

    private :

        /** Element to Element mapping data. Array size: [m_uiAllNodes.size()*m_uiStensilSz*m_uiNumDirections];  But this is done for m_uiStencilSz=1*/
        std::vector<unsigned int> m_uiE2EMapping;
        /** Element ot Node mapping data for continous Galerkin methods. Array size: [m_uiAllNodes.size()*m_uiNpE];*/
        std::vector<unsigned int> m_uiE2NMapping_CG;
        /** Element to Node mapping with DG indexing after removing duplicates. This is used for debugging. */
        std::vector<unsigned int> m_uiE2NMapping_DG;
        /** cg to dg mapping*/
        std::vector<unsigned int> m_uiCG2DG;
        /** dg to cg mapping*/
        std::vector<unsigned int> m_uiDG2CG;

        /** splitter element for each processor. */
        std::vector<ot::TreeNode> m_uiLocalSplitterElements;                                 // used to spit the keys to the correct nodes.

        /**Splitter Node for each processor*/
        ot::TreeNode *m_uiSplitterNodes;

        // Pre  and Post ghost octants.

        /** pre ghost elements (will be freed after use)*/
        std::vector<ot::TreeNode> m_uiPreGhostOctants;
        /**post ghost elements (will be freed after use)*/
        std::vector<ot::TreeNode> m_uiPostGhostOctants;

        // Keys Related attributes
        /**search keys generated for local elements. */
        std::vector<Key> m_uiKeys;

        /**search keys generated for ghost elements*/
        std::vector<Key>m_uiGhostKeys;

        /**search keys generated for missing diagonal elements in the ghost*/
        std::vector<Key>m_uiKeysDiag;

        /**input to the mesh generation (will be freed after use)*/
        std::vector<ot::TreeNode> m_uiEmbeddedOctree;

        /** Ghost Elements. (will be freed after use)*/
        std::vector<ot::TreeNode> m_uiGhostOctants;

        /** store the indices of round 1 communication.
         *
         * Note: !! You cannot build a m_uiGhostElementRound2Index, because at the round 2 exchange of ghost elements we might
         * need to send a round1 ghost element to another processor.
         *
         * */
        std::vector<unsigned int > m_uiGhostElementRound1Index;


        /** stores all the pre + local + post ELEMENTS. */
        std::vector<ot::TreeNode> m_uiAllElements;
        /** stores the local nodes */
        std::vector<ot::TreeNode> m_uiAllLocalNode;

        /** block list for the current local elements in order to apply the stencil. */
        std::vector<ot::Block> m_uiLocalBlockList;

        /**minimum depth of the octree*/
        unsigned int m_uiDmin;
        /**maximum depth of the octree*/
        unsigned int m_uiDmax;

        /**begin of the pre ghost elements*/
        unsigned int m_uiElementPreGhostBegin = 0; // Not mandatory to store this.
        /** end of pre ghost elements*/
        unsigned int m_uiElementPreGhostEnd;
        /** begin of locat octants (INDEPENDENT)*/
        unsigned int m_uiElementLocalBegin;
        /** end of the local octants (INDEPENDET)*/
        unsigned int m_uiElementLocalEnd;
        /**begin location for the post ghost octants*/
        unsigned int m_uiElementPostGhostBegin;
        /** end location for the post ghost octants*/
        unsigned int m_uiElementPostGhostEnd;
        /**begin of the pre ghost elements of  fake elements*/
        unsigned int m_uiFElementPreGhostBegin;
        /**end of the pre ghost elements of fake elements. */
        unsigned int m_uiFElementPreGhostEnd;
        /**begin of the local fake elements*/
        unsigned int m_uiFElementLocalBegin;
        /**end of the local fake elements*/
        unsigned int m_uiFElementLocalEnd;
        /**begin of the fake post ghost elements*/
        unsigned int m_uiFElementPostGhostBegin;
        /**end of the fake post ghost elements*/
        unsigned int m_uiFElementPostGhostEnd;

        /**denotes the mesh domain min*/
        unsigned int m_uiMeshDomain_min;
        /**denotes the mesh domain max*/
        unsigned int m_uiMeshDomain_max;


        /**Number of local elements. */
        unsigned int m_uiNumLocalElements;
        /**Number of pre ghost elements*/
        unsigned int m_uiNumPreGhostElements;
        /**Number of post ghost elements*/
        unsigned int m_uiNumPostGhostElements;

        /*** Total number of actual elements in the mesh (preG+ local + postG ) */
        unsigned int m_uiNumTotalElements;

        /** Total number of actual nodes. Defn: Actual nodes are the ones that are owned by actual elememts. */
        unsigned int m_uiNumActualNodes;

        /** unzipped vector size. */
        DendroIntL m_uiUnZippedVecSz;

        /** Total number of fake nodes. Defn: these nodes are owned my the fake elements. */
        unsigned int m_uiNumFakeNodes;


        /** begin location  of the pre ghost nodes in CG indexing*/
        unsigned int m_uiNodePreGhostBegin;
        /** end location of the pre ghost nodes in CG indexing*/
        unsigned int m_uiNodePreGhostEnd;
        /** begin location of the local nodes in CG indexing*/
        unsigned int m_uiNodeLocalBegin;
        /** end location of the local nodes in CG indexing*/
        unsigned int m_uiNodeLocalEnd;
        /** begin location of the post ghost nodes in CG indexing*/
        unsigned int m_uiNodePostGhostBegin;
        /** end location of the post ghost nodes in CG indexing*/
        unsigned int m_uiNodePostGhostEnd;

        // Stuff for communication and synchronization ...
        //-------------------------------------------------------------------------

        /**Rank of the current process*/
        int m_uiActiveRank;
        /** size of the active  Comm */
        int m_uiActiveNpes;

        /** MPI active communicator. (this is a subcomm of the m_uiCommGlobal) */
        MPI_Comm m_uiCommActive;

        /** MPI global communicator*/
        MPI_Comm m_uiCommGlobal;

        /**Rank of the current process*/
        int m_uiGlobalRank;
        /** size of the active  Comm */
        int m_uiGlobalNpes;

        /** indicates whether the current mesh is active under the m_uiCommActive*/
        bool m_uiIsActive;

        /**Number of keys should be sent to each processor. (Note that instead of keys we are sending the owners of keys to the required processor)*/
        unsigned int *m_uiSendKeyCount;
        /**Number of keys received by each processor*/
        unsigned int *m_uiRecvKeyCount;
        /** send key offset array*/
        unsigned int *m_uiSendKeyOffset;
        /**receive key offset array*/
        unsigned int *m_uiRecvKeyOffset;


        /**send diagonal key count*/
        unsigned int *m_uiSendKeyDiagCount;
        /**recv diagonal key count*/
        unsigned int *m_uiRecvKeyDiagCount;
        /** send key diagonal offset*/
        unsigned int *m_uiSendKeyDiagOffset;
        /** recv key diagonal offset*/
        unsigned int *m_uiRecvKeyDiagOffset;



        /**SendOct counts related to round 1 of ghost exchange. */
        /**Number of ghost octants(elements) that needed to be sent to each processor */
        unsigned int *m_uiSendOctCountRound1;
        /**Number of ghost elemets recieved from each processor. */
        unsigned int *m_uiRecvOctCountRound1;
        /**Ghost send count offset, used in all2allv*/
        unsigned int *m_uiSendOctOffsetRound1;
        /**Ghost recieve count offset, used in all2allv*/
        unsigned int *m_uiRecvOctOffsetRound1;


        unsigned int *m_uiSendOctCountRound1Diag;
        /**Number of ghost elemets recieved from each processor. */
        unsigned int *m_uiRecvOctCountRound1Diag;
        /**Ghost send count offset, used in all2allv*/
        unsigned int *m_uiSendOctOffsetRound1Diag;
        /**Ghost recieve count offset, used in all2allv*/
        unsigned int *m_uiRecvOctOffsetRound1Diag;



        /**SendOct counts related to round 2 of ghost exchange. */
        /**Number of ghost octants(elements) that needed to be sent to each processor */
        unsigned int *m_uiSendOctCountRound2;
        /**Number of ghost elemets recieved from each processor. */
        unsigned int *m_uiRecvOctCountRound2;
        /**Ghost send count offset, used in all2allv*/
        unsigned int *m_uiSendOctOffsetRound2;
        /**Ghost recieve count offset, used in all2allv*/
        unsigned int *m_uiRecvOctOffsetRound2;




        /**Number of nodes that needed to be sent to each processor*/
        std::vector<unsigned int>  m_uiSendNodeCount;
        /**Number of nodes that recieved from each processor*/
        std::vector<unsigned int> m_uiRecvNodeCount;
        /**Send node count offset*/
        std::vector<unsigned int> m_uiSendNodeOffset;
        /**Recv node count offset*/
        std::vector<unsigned int> m_uiRecvNodeOffset;

        /**Send processor list for ghost exchange*/
        std::vector<unsigned int> m_uiSendProcList;

        /**recv processor list for the ghost exchange*/
        std::vector<unsigned int> m_uiRecvProcList;

        /** local element ids that needs to be send for other processors. */
        std::vector<unsigned int> m_uiGhostElementIDsToBeSent;
        /** local element ids that recieved from other processors after ghost element exchange*/
        std::vector<unsigned int> m_uiGhostElementIDsToBeRecv;

        /**stores the CG indecies of pre ghost nodes which are hanging, */
        std::vector<unsigned int > m_uiPreGhostHangingNodeCGID;

        /**stores the CG indices of post ghost nodes which are hanging.*/
        std::vector<unsigned int > m_uiPostGhostHangingNodeCGID;



        /**Send buffer of ghost elements (octants)*/
        std::vector<ot::TreeNode> m_uiSendBufferElement;

        /**Actual exchange of node values. This is the actual exhange that happens in element looping*/
        std::vector<double> m_uiSendBufferNodes;
        /**recv buffer for the ghost node values receiveing from others. */
        std::vector<double> m_uiRecvBufferNodes;

        /**Scatter map for the elements. Keeps track of which local elements need to be sent to which processor. (This has to be derived with m_uiSendOctCount values)
         * Element ID's that is been exchanged at round 1 of ghost exchange.
         * Property: For a given processor rank p the scatter map elements that is set to processor p is sorted and unique.
         * */
        std::vector<unsigned int> m_uiScatterMapElementRound1;


        /**Scatter map for the actual nodes. Keeps track of which local node needs to be sent to which processor. */
        std::vector<unsigned int > m_uiScatterMapActualNodeSend;

        /** Scatter map for the actual nodes, recieving from other processors. */
        std::vector<unsigned int > m_uiScatterMapActualNodeRecv;


        // variables to manage loop access over elements.
        /**counter for the current element*/
        unsigned int m_uiEL_i;

        /**order of the octant elements. */
        unsigned int m_uiElementOrder;

        /**Number of nodes per element*/
        unsigned int m_uiNpE;

        /**Number of neighbours that need to be searched in each direction. */
        unsigned int m_uiStensilSz;

        /**number of total directions that the search is performed*/
        unsigned int m_uiNumDirections;



        /**Reference element to perform interpolations. */
        RefElement m_uiRefEl;


    private:


        /**
         * @author Milinda Fernando
         * @brief generates search key elements for local elements.
         *
         *
         * */

        void generateSearchKeys();

        /**
         * @author Milinda Fernando
         * @brief generates search keys for ghost elements, inorder to build the E2E mapping between ghost elements.
         * */

        void generateGhostElementSearchKeys();

        /**
         * @brief generates diagonal keys for ghost if it is not already in the current processor.
         * */

        void generateBdyElementDiagonalSearchKeys();



        /**
         * @author Milinda Fernando
         * @brief Builds the E2E mapping, in the mesh
         * @param [in] in: 2:1 balanced octree (assumes that the input is 2:1 balanced unique and sorted)
         * @param [in] k_s : Stencil size, how many elements to search in each direction
         * @param [in] comm: MPI Communicator.
         * */
        void buildE2EMap(std::vector<ot::TreeNode> &in,MPI_Comm comm );


        void computeElementOwnerRanks(std::vector<unsigned int > & elementOwner);




        /**
         * @author Milinda Fernando
         * @brief Builds the E2E mapping, in the mesh (sequential case. No ghost nodes. )
         * @param [in] in: 2:1 balanced octree (assumes that the input is 2:1 balanced unique and sorted)
         * @param [in] k_s : Stencil size, how many elements to search in each direction
         * */
         void buildE2EMap(std::vector<ot::TreeNode> &in);

        /**
         *
         * @author Milinda Fernando
         * @brief Builds the Element to Node (E2N) mapping to enforce the continuity of the solution. (Needed in continous Galerkin methods. )
         * This function assumes that E2E mapping is already built.
         *
         * */

        void buildE2NMap();

        /**
         * @author Milinda Fernando
         * @brief inorder to make the E2N consistent across all processors we have exchanged additional layer of ghost elements. This procedure remove
         * those additional elements( layer 2 ghost) from global octant array and update E2E and E2N accordingly.
         *
         * */
        void shrinkE2EAndE2N();


        /**
         * @author Milinda Fernando
         * @brief Perform the blocks initialization so that we can apply the stencil for the grid as a sequnce of finite number of regular grids.
         * note that this should be called after performing all E2N and E2N mapping.
         * */
        void performBlocksSetup();




        /**
         *
         * @author: Milinda Fernando.
         * @def Fake Elements(Local) are teh fake elements where all it's nodal lies in the nodal local and nodal L1 ghost layer.
         * @brief: Generates local fake elements and update E2E and E2N mapping. (This is not used. )
         * */

        void generateFakeElementsAndUpdateE2EAndE2N();


        /**
         * @author  Milinda Fernando
         * @brief Computes the overlapping nodes between given two elements. Note that idx ,idy, idz should be integer array size of (m_uiElementOrder + 1).
         * @param [in] Parent : denotes the larger(ancestor) element of the two considering elements.
         * @param [in] child : denotes the desendent from the parent element.
         * @param [out] idx: mapping between x nodes indecies of parent element to x node indecies for the child element.
         * @param [out] idy: mapping between y nodes indecies of parent element to y node indecies for the child element.
         * @param [out] idz: mapping between z nodes indecies of parent element to z node indecies for the child element.
         * @example idx[0]= 2 implies that the parent element x index 0 is given from the child x index 2.
         *
         * */
        inline bool computeOveralppingNodes(const ot::TreeNode & parent, const ot::TreeNode & child, int * idx,int * idy, int * idz);

        /**
         *  @breif Computes the ScatterMap and the send node counts for ghost node exchage.
         *  We do the scatter map exchange in two different steps. Hence we need to compute 2 different scatter maps.
         *      1. Scatter map for actual nodes.
         *      2. Scatter map for fake nodes.
         *  @param[in] MPI_Communicator for scatter map node exchange.
         *  @note: depreciated remove this method later.
         *
         */
        void computeNodeScatterMaps(MPI_Comm comm);

        /**
         *  @breif Computes the ScatterMap and the send node counts for ghost node exchage.
         *  @param[in] comm : MPI communicator.
         * */
        void computeNodalScatterMap(MPI_Comm comm);


         /**
               *  @breif Computes the ScatterMap and the send node counts for ghost node exchage. (uses the compression of the all to all data exchange to reduce the total number of nodes get exchanged. )
               *  @param[in] comm : MPI communicator.
               *
        * */
       void computeNodalScatterMap1(MPI_Comm comm);

     /**
           *@breif Computes the ScatterMap and the send node counts for ghost node exchage. (uses the compression of the all to all data exchange to reduce the total number of nodes get exchanged. )
           *@param[in] comm : MPI communicator.
           *
        **/
     void computeNodalScatterMap2(MPI_Comm comm);


     /**
          *@breif Computes the ScatterMap and the send node counts for ghost node exchage. (uses the compression of the all to all data exchange to reduce the total number of nodes get exchanged. )
          *@param[in] comm : MPI communicator. scattermap2 with compression.
          *
       **/
     void computeNodalScatterMap3(MPI_Comm comm);


     /**
          *@breif Computes the ScatterMap and the send node counts for ghost node exchage. (uses the compression of the all to all data exchange to reduce the total number of nodes get exchanged. )
          *@param[in] comm : MPI communicator. scattermap2 with compression.
          *
       **/
     void computeNodalScatterMap4(MPI_Comm comm);


        /**
         * childf1 and childf2 denotes the edge that happens ar the intersection of those two planes.
         * */

        /*
         * NOTE: These private functions does not check for m_uiIsActive
         * sicne they are called in buildE2E and buildE2N mappings.
         *
         * */

        /**
         * @breif This function is to map internal edge nodes in the LEFT face.
         * @param[in] child element in cordieration for mapping.
         * @param[in] parent or the LEFT neighbour of child
         * @param[in] parentChildLevEqual specifies whether the parent child levels are equal (true) if their equal.
         * @param[in] edgeChildIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
         * @param[in] edgeOwnerIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
         * */
        inline void OCT_DIR_LEFT_INTERNAL_EDGE_MAP(unsigned int child,unsigned int parent,bool parentChildLevEqual,std::vector<unsigned int >& edgeChildIndex,std::vector<unsigned int>& edgeOwnerIndex);
        /**
        * @breif This function is to map internal edge nodes in the RIGHT face.
        * @param[in] child element in cordieration for mapping.
        * @param[in] parent or the RIGHT neighbour of child
        * @param[in] parentChildLevEqual specifies whether the parent child levels are equal (true) if their equal.
        * @param[in] edgeChildIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * @param[in] edgeOwnerIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * */
        inline void OCT_DIR_RIGHT_INTERNAL_EDGE_MAP(unsigned int child,unsigned int parent,bool parentChildLevEqual,std::vector<unsigned int >& edgeChildIndex,std::vector<unsigned int>& edgeOwnerIndex);
        /**
        * @breif This function is to map internal edge nodes in the DOWN face.
        * @param[in] child element in cordieration for mapping.
        * @param[in] parent or the DOWN neighbour of child
        * @param[in] parentChildLevEqual specifies whether the parent child levels are equal (true) if their equal.
        * @param[in] edgeChildIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * @param[in] edgeOwnerIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * */
        inline void OCT_DIR_DOWN_INTERNAL_EDGE_MAP(unsigned int child,unsigned int parent,bool parentChildLevEqual,std::vector<unsigned int >& edgeChildIndex,std::vector<unsigned int>& edgeOwnerIndex);
        /**
        * @breif This function is to map internal edge nodes in the UP face.
        * @param[in] child element in cordieration for mapping.
        * @param[in] parent or the UP neighbour of child
        * @param[in] parentChildLevEqual specifies whether the parent child levels are equal (true) if their equal.
        * @param[in] edgeChildIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * @param[in] edgeOwnerIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * */
        inline void OCT_DIR_UP_INTERNAL_EDGE_MAP(unsigned int child,unsigned int parent,bool parentChildLevEqual,std::vector<unsigned int >& edgeChildIndex,std::vector<unsigned int>& edgeOwnerIndex);
        /**
        * @breif This function is to map internal edge nodes in the FRONT face.
        * @param[in] child element in cordieration for mapping.
        * @param[in] parent or the FRONT neighbour of child
        * @param[in] parentChildLevEqual specifies whether the parent child levels are equal (true) if their equal.
        * @param[in] edgeChildIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * @param[in] edgeOwnerIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
        * */
        inline void OCT_DIR_FRONT_INTERNAL_EDGE_MAP(unsigned int child,unsigned int parent,bool parentChildLevEqual,std::vector<unsigned int >& edgeChildIndex,std::vector<unsigned int>& edgeOwnerIndex);
        /**
         * @breif This function is to map internal edge nodes in the BACK face.
         * @param[in] child element in cordieration for mapping.
         * @param[in] parent or the BACK neighbour of child
         * @param[in] parentChildLevEqual specifies whether the parent child levels are equal (true) if their equal.
         * @param[in] edgeChildIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
         * @param[in] edgeOwnerIndex std::vector (pre-allocated) for the size of the internal nodes of the edge.
         * */
        inline void OCT_DIR_BACK_INTERNAL_EDGE_MAP(unsigned int child,unsigned int parent,bool parentChildLevEqual,std::vector<unsigned int >& edgeChildIndex,std::vector<unsigned int>& edgeOwnerIndex);


        /**
         * @brief maps the corner nodes of an elements to a  corresponsing elements.
         * @param[in] child element id in corsideration.
         *  */
        inline void CORNER_NODE_MAP(unsigned int child);


        /**
         * @brief Returns the diagonal element ID if considering two faces are simillar in size of elementID.
         * */
        inline void OCT_DIR_DIAGONAL_E2E(unsigned int elementID, unsigned int face1,unsigned int face2,unsigned int& lookUp) const ;


       /**
             * @brief Search the given set of keys in given nodes and update the search result of the keys. The speciality of this method from SFC_TreeSearch() is that this won't change the ordering of the
             * Keys. Hence this is expensive than the SFC_TreeSearch. So use it only when the ordering of the keys matters.
             *
             * @param [in] pKeys: Keys needed to be searched.
             * @param [in] pNodes: Nodes (represents the space where keys being searched)
             * */

        template <typename pKey, typename pNode>
        void searchKeys(std::vector<pKey>& pKeys,std::vector<pNode>& pNodes);

         /**
             * @brief Returns the direction of the node when i, j, k index of the node given.
             * */
         inline unsigned int getDIROfANode(unsigned int ii_x,unsigned int jj_y,unsigned int kk_z);

         /**
          * @brief decompose the direction to it's ijk values.
          * */

         inline void directionToIJK(unsigned int direction,std::vector<unsigned int > &ii_x, std::vector<unsigned int >& jj_y, std::vector<unsigned int > &kk_z);

         /**
         * @brief: Specific to GR code. Interpolation function written to interpolate the 3rd points for 4th order element.
         * @param [in] upWind: advection up wind stencil coefficients
         * @param [in] element : unzip element ID
         * @param [in] lookUp: face lookup ID
         * @param [in] vecLookUp:  allocated vector to lookup Nodal values.
         * @param [in] cnum: injection element cnum
         * @param [in] parentInterpIn: parent values.
         * @param [in] parentInterpOut: vector to perform parent to child interpolation
         * @param [in] padDir : direction of the padding.
         * @param [in] padWidth: width of the padding.
         * @param [in] zippedVec: zipped vector
         * @param [out] out : child to parent injection filled with the correct value.
         * */
        template<typename T>
        inline void interpUpWind(const double *upWind,const unsigned int element,const unsigned int lookup, T* vecLookUp, const unsigned int cnum, const T* parentInterpIn, T* parentInterpOut,const unsigned int padDir,const unsigned int padWidth,const T* zippedVec,T* out);

         /**
         * @brief: Specific to GR code. Interpolation function written to interpolate the 3rd points for 4th order element.
         * @param [in] upWind: advection up wind stencil coefficients
         * @param [in] element : unzip element ID
         * @param [in] lookUp: face lookup ID
         * @param [in] vecLookUp:  allocated vector to lookup Nodal values.
         * @param [in] cnum: injection element cnum
         * @param [in] parentInterpIn: parent values.
         * @param [in] parentInterpOut: vector to perform parent to child interpolation
         * @param [in] padDir : direction of the padding.
         * @param [in] padWidth: width of the padding.
         * @param [in] zippedVec: zipped vector
         * @param [out] out : child to parent injection filled with the correct value.
         * */
        template<typename T>
        inline void interpDownWind(const double *downWind,const unsigned int element,const unsigned int lookup, T* vecLookUp, const unsigned int cnum, const T* parentInterpIn, T* parentInterpOut,const unsigned int padDir,const unsigned int padWidth,const T* zippedVec,T* out);

        /**@brief: Performs block padding along the diagonal direction.
         * @param[in] blk: block to perform diagonal direction padding.
         * @param[in] zippedVec: zipped vector
         * @param[out] unzippedVec: updated unzip vec.
         * */
        template<typename T>
        void blockDiagonalUnZip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the vertex direction.
         * @param[in] blk: block to perform diagonal direction padding.
         * @param[in] zippedVec: zipped vector
         * @param[out] unzippedVec: updated unzip vec.
         * */
        template<typename T>
        void blockVertexUnZip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);


         /**@brief: Performs block padding along the OCT_DIR_LEFT_DOWN direction.
         * @param[in] blk: block to perform diagonal direction padding.
         * @param[in] zippedVec: zipped vector
         * @param[out] unzippedVec: updated unzip vec.
         * */
        template<typename T>
        void OCT_DIR_LEFT_DOWN_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

        /**@brief: Performs block padding along the OCT_DIR_LEFT_UP direction.
         * @param[in] blk: block to perform diagonal direction padding.
         * @param[in] zippedVec: zipped vector
         * @param[out] unzippedVec: updated unzip vec.
         * */
        template<typename T>
        void OCT_DIR_LEFT_UP_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

        /**@brief: Performs block padding along the OCT_DIR_LEFT_BACK direction.
         * @param[in] blk: block to perform diagonal direction padding.
         * @param[in] zippedVec: zipped vector
         * @param[out] unzippedVec: updated unzip vec.
         * */
        template<typename T>
        void OCT_DIR_LEFT_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_LEFT_FRONT direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
        template<typename T>
        void OCT_DIR_LEFT_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);


         /**@brief: Performs block padding along the OCT_DIR_RIGHT_DOWN direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_RIGHT_DOWN_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_RIGHT_UP direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_RIGHT_UP_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_RIGHT_BACK direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_RIGHT_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_RIGHT_FRONT direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_RIGHT_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_UP_BACK direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_UP_FRONT direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_DOWN_BACK direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_DOWN_FRONT direction.
          * @param[in] blk: block to perform diagonal direction padding.
          * @param[in] zippedVec: zipped vector
          * @param[out] unzippedVec: updated unzip vec.
          * */
         template<typename T>
         void OCT_DIR_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);



         /**@brief: Performs block padding along the OCT_DIR_LEFT_DOWN_BACK direction.
              * @param[in] blk: block to perform diagonal direction padding.
              * @param[in] zippedVec: zipped vector
              * @param[out] unzippedVec: updated unzip vec.
              * */
         template<typename T>
         void OCT_DIR_LEFT_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_RIGHT_DOWN_BACK direction.
              * @param[in] blk: block to perform diagonal direction padding.
              * @param[in] zippedVec: zipped vector
              * @param[out] unzippedVec: updated unzip vec.
              * */
        template<typename T>
        void OCT_DIR_RIGHT_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);


         /**@brief: Performs block padding along the OCT_DIR_LEFT_UP_BACK direction.
                 * @param[in] blk: block to perform diagonal direction padding.
                 * @param[in] zippedVec: zipped vector
                 * @param[out] unzippedVec: updated unzip vec.
                 * */
         template<typename T>
         void OCT_DIR_LEFT_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);


         /**@brief: Performs block padding along the OCT_DIR_RIGHT_UP_BACK direction.
                 * @param[in] blk: block to perform diagonal direction padding.
                 * @param[in] zippedVec: zipped vector
                 * @param[out] unzippedVec: updated unzip vec.
                 * */
         template<typename T>
         void OCT_DIR_RIGHT_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);


         /**@brief: Performs block padding along the OCT_DIR_LEFT_DOWN_FRONT direction.
                  * @param[in] blk: block to perform diagonal direction padding.
                  * @param[in] zippedVec: zipped vector
                  * @param[out] unzippedVec: updated unzip vec.
                  * */
         template<typename T>
         void OCT_DIR_LEFT_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);

         /**@brief: Performs block padding along the OCT_DIR_RIGHT_DOWN_FRONT direction.
              * @param[in] blk: block to perform diagonal direction padding.
              * @param[in] zippedVec: zipped vector
              * @param[out] unzippedVec: updated unzip vec.
              * */
         template<typename T>
         void OCT_DIR_RIGHT_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);


         /**@brief: Performs block padding along the OCT_DIR_LEFT_UP_FRONT direction.
                 * @param[in] blk: block to perform diagonal direction padding.
                 * @param[in] zippedVec: zipped vector
                 * @param[out] unzippedVec: updated unzip vec.
                 * */
         template<typename T>
         void OCT_DIR_LEFT_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);


         /**@brief: Performs block padding along the OCT_DIR_RIGHT_UP_FRONT direction.
                 * @param[in] blk: block to perform diagonal direction padding.
                 * @param[in] zippedVec: zipped vector
                 * @param[out] unzippedVec: updated unzip vec.
                 * */
         template<typename T>
         void OCT_DIR_RIGHT_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec);



    public:


     /**@brief parallel mesh constructor
         * @param[in] in: complete sorted 2:1 balanced octree to generate mesh
         * @param[in] k_s: how many neighbours to check on each direction (used =1)
         * @param[in] pOrder: order of an element.
         * @param[in] commActive: MPI active communicator
         * @param[in] comm: MPI communicator (global)

         * */
     Mesh(std::vector<ot::TreeNode> &in, unsigned int k_s, unsigned int pOrder,unsigned int activeNpes,MPI_Comm comm,unsigned int grainSz=DENDRO_DEFAULT_GRAIN_SZ,double ld_tol=DENDRO_DEFAULT_LB_TOL,unsigned int sf_k=DENDRO_DEFAULT_SF_K);


        /**@brief parallel mesh constructor
         * @param[in] in: complete sorted 2:1 balanced octree to generate mesh
         * @param[in] k_s: how many neighbours to check on each direction (used =1)
         * @param[in] pOrder: order of an element.
         * @param[in] comm: MPI communicator (global)
         * @param[in] grainSz: prefered grain sz. (this parameter is used to perform automatic comm expansion and shrinking)
         * @param[in] ld_tol: load imbalance tolerance for comm expansion and shrinking
         * @param[in] sf_k: splitter fix _k value. (Needed by SFC_partitioinng for large p>=10,000)
         * */
        Mesh(std::vector<ot::TreeNode> &in, unsigned int k_s, unsigned int pOrder,MPI_Comm comm,unsigned int grainSz=DENDRO_DEFAULT_GRAIN_SZ,double ld_tol=DENDRO_DEFAULT_LB_TOL,unsigned int sf_k=DENDRO_DEFAULT_SF_K);

        /**@brief destructor for mesh (releases the allocated variables in the class. )*/
        ~Mesh();

        // Setters and getters.
        /** @breif Returns the number of local elements in the grid. (local to the current considering processor)*/
        inline unsigned int getNumLocalMeshElements() const {return m_uiNumLocalElements;}
        /** @breif Returns the number of pre-ghost elements */
        inline unsigned int getNumPreGhostElements() const {return m_uiNumPreGhostElements;}
        /** @brief Returns the number of post-ghost elements */
        inline unsigned int getNumPostGhostElements() const {return m_uiNumPostGhostElements;}
        // get nodal information.

        /**@brief return the number of nodes local to the mesh*/
        inline unsigned int getNumLocalMeshNodes() const {return (m_uiNodeLocalEnd-m_uiNodeLocalBegin);}
        /**@brief return the number of pre ghost mesh nodes*/
        inline unsigned int getNumPreMeshNodes() const {return (m_uiNodePreGhostEnd-m_uiNodePreGhostBegin);}
        /**@brief return the number of post ghost mesh nodes*/
        inline unsigned int getNumPostMeshNodes()const {return (m_uiNodePostGhostEnd-m_uiNodePostGhostBegin);}
        /**@brief return the begin location of element pre ghost*/
        inline unsigned int getElementPreGhostBegin() const {return m_uiElementPreGhostBegin;}
        /**@brief return the end location of element pre ghost*/
        inline unsigned int getElementPreGhostEnd()const {return m_uiElementPreGhostEnd;}
        /**@brief return the begin location of element local*/
        inline unsigned int getElementLocalBegin() const {return m_uiElementLocalBegin;}
        /**@brief return the end location of element local*/
        inline unsigned int getElementLocalEnd()const {return m_uiElementLocalEnd;}
        /**@brief return the begin location of element post ghost */
        inline unsigned int getElementPostGhostBegin()const {return m_uiElementPostGhostBegin;}
        /**@brief return the end location of element post ghost*/
        inline unsigned int getElementPostGhostEnd()const {return m_uiElementPostGhostEnd;}

        /**@brief return the begin location of pre ghost nodes*/
        inline unsigned int getNodePreGhostBegin()const {return m_uiNodePreGhostBegin;}
        /**@brief return the end location of pre ghost nodes*/
        inline unsigned int getNodePreGhostEnd()const {return m_uiNodePreGhostEnd;}
        /**@brief return the location of local node begin*/
        inline unsigned int getNodeLocalBegin()const {return m_uiNodeLocalBegin;}
        /**@brief return the location of local node end*/
        inline unsigned int getNodeLocalEnd() const {return m_uiNodeLocalEnd;}
        /**@brief return the location of post node begin*/
        inline unsigned int getNodePostGhostBegin()const {return m_uiNodePostGhostBegin;}
        /**@brief return the location of post node end*/
        inline unsigned int getNodePostGhostEnd()const { return m_uiNodePostGhostEnd;}

        /**@brief returns the dof for a partition*/
        inline unsigned int getDegOfFreedom() const {return m_uiNumActualNodes;}

        /**@brief returns the dof for a partition*/
        inline unsigned int getDegOfFreedomUnZip() const {return m_uiUnZippedVecSz;}

        /**@brief returns the pointer to All elements array. */
        inline const std::vector<ot::TreeNode>& getAllElements()const {return m_uiAllElements;}

        /**@brief returns the Level 1 ghost element indices. */
        inline const std::vector<unsigned int >& getLevel1GhostElementIndices(){return m_uiGhostElementRound1Index;}

        /**@brief returns the splitter elements of the mesh local elements. */
        inline const std::vector<ot::TreeNode>& getSplitterElements() const {return m_uiLocalSplitterElements;}

        /**@brief returns all local nodes(vertices)*/
        inline const std::vector<ot::TreeNode>& getAllLocalNodes() const { return m_uiAllLocalNode;}

        /**@brief returns const e2e mapping instance. */
        inline const std::vector<unsigned int> & getE2EMapping()const {return m_uiE2EMapping;}

        /**@breif returns const e2n mapping instance */
        inline const std::vector<unsigned int> & getE2NMapping()const {return m_uiE2NMapping_CG;}

        /**@breif returns const e2n mapping instance (debuging purposes only) */
        inline const std::vector<unsigned int> & getE2NMapping_DG()const {return m_uiE2NMapping_DG;}
        /**@brief returns const list of local blocks (regular grids) for the consdering mesh. */
        inline const std::vector<ot::Block> & getLocalBlockList() const {return m_uiLocalBlockList;}

        /**@brief return the number of directions in the E2E mapping. */
        inline unsigned int getNumDirections()const {return m_uiNumDirections;}

        /**@brief return the number of nodes per element.*/
        inline unsigned int getNumNodesPerElement()const {return m_uiNpE;}

        /**@brief returns the order of an element*/
        inline unsigned int getElementOrder() const { return m_uiElementOrder;}

        /**@brief returns the communicator (acitve) */
        inline MPI_Comm getMPICommunicator() const {return m_uiCommActive;}

        /**@brief returns the global communicator*/
        inline MPI_Comm getMPIGlobalCommunicator() const {return m_uiCommGlobal;}

        /**@brief returns the rank w.r.t. global comm*/
        inline unsigned int getMPIRankGlobal() const {return m_uiGlobalRank;}

        /** @brief returns the comm size w.r.t. global comm*/
        inline unsigned int getMPICommSizeGlobal() const {return m_uiGlobalNpes;}

        /** @brief returns the rank */
        inline unsigned int getMPIRank() const {return m_uiActiveRank;}

        /** @brief returns the comm size: */
        inline unsigned int getMPICommSize() const {return m_uiActiveNpes;}

        /** @brief returns const pointer to reference element */
        inline const RefElement * getReferenceElement() const {return &m_uiRefEl;}

        /**@brief returns the send proc list size*/
        inline const unsigned int getSendProcListSize() const {return m_uiSendProcList.size();}

        /**@brief returns the recv proc list size*/
        inline const unsigned int getRecvProcListSize() const {return m_uiRecvProcList.size();}


        /**
         * @brief: Decompose the DG index to element id and it's i,j,k values.
         * */
        inline void dg2eijk(unsigned int dg_index,unsigned int& e, unsigned int & i,unsigned int & j, unsigned int &k) const
        {
            e=dg_index/m_uiNpE;
            k=0;j=0;i=0;

            //std::cout<<"e: "<<e<<std::endl;
            if(dg_index>e*m_uiNpE) k=(dg_index-e*m_uiNpE)/((m_uiElementOrder+1)*(m_uiElementOrder+1));
            //std::cout<<"k: "<<k<<std::endl;
            if((dg_index + k*((m_uiElementOrder+1)*(m_uiElementOrder+1)))> (e*m_uiNpE)) j=(dg_index-e*m_uiNpE-k*((m_uiElementOrder+1)*(m_uiElementOrder+1)))/(m_uiElementOrder+1);
            //std::cout<<"j: "<<j<<std::endl;
            if((dg_index+ k*((m_uiElementOrder+1)*(m_uiElementOrder+1))+j*(m_uiElementOrder+1)) >(e*m_uiNpE) ) i=(dg_index-e*m_uiNpE-k*((m_uiElementOrder+1)*(m_uiElementOrder+1))-j*(m_uiElementOrder+1));
            //std::cout<<"i: "<<i<<std::endl;

        }

        /**@brief returns true if mesh is active*/
        inline bool isActive() const {return m_uiIsActive;}

        /**@brief waiting for all the mesh instances both active and inactive. This should not be called if not needed. this is a BARRIER. */
        inline void waitAll() const {MPI_Barrier(m_uiCommGlobal);}



            /**
             * @brief returns the face neighours in specidied direction.
             * @param[in] eID:  Element ID.
             * @param[in] dir: direction of the face.
             * @param[out] lookUp: result if the result is a same level or lower level than element ID. (only 2 neighbor)
             **/

        void getElementalFaceNeighbors(const unsigned int eID,const unsigned int dir, unsigned int* lookup)const;
            /**
           * @brief returns the edge neighours in specidied direction.
           * @param[in] eID:  Element ID.
           * @param[in] dir: direction of the edge.
           * @param[out] lookUp: result if the result is a same level or lower level than element ID. (only 4 neighbor)
           * */

        void getElementalEdgeNeighbors(const unsigned int eID,const unsigned int dir, unsigned int* lookup) const ;
            /**
            * @brief returns the vertex neighours in specidied direction.
            * @param[in] eID:  Element ID.
            * @param[in] dir: direction of the vertex.
            * @param[out] lookUp: result if the result is a same level or lower level than element ID. (only 8 neighbor)
            * */
        void getElementalVertexNeighbors(const unsigned int eID,const unsigned int dir, unsigned int* lookup)const ;

         // Wavelet Init functions
         /**
          * @brief Initilizes the wavelet DA loop depending on the WaveletDA flags specified
          * */
         template <ot::WaveletDA::LoopType type>
         void init();

         /**
          * @brief Check whether the next element is available.
          * */

         template <ot::WaveletDA::LoopType type>
         bool nextAvailable();

         /**
          * @brief Increment the counters to access the next element in the mesh.
          * */

         template <ot::WaveletDA::LoopType  type>
         void next();

         /**
          * @brief Return the current element as an octant.
          * */
         inline const ot::TreeNode& currentOctant();

         /**
          * @brief Returns the current element index.
          * */
         inline unsigned  int currentIndex();

         /**
          * @brief Returns the current neighbour list (Element) information. Note that for 3D it is 8 neighbours for each octant and for 2D it is 4.
          * */

         inline void currentElementNeighbourIndexList(unsigned int *neighList);

         /**
          * @brief: Returns the node index list belongs to current element.
          * NodeList size should be m_uiNpE;
          * */

         inline void currentElementNodeList(unsigned int * nodeList);


         /**
          * @brief Returns the node index list belogns to the currentl element in DG indexing.
          * NodeList size should be m_uiNpE;
          *
          * */
         inline void currentElementNodeList_DG(unsigned int * nodeList);



     // functions to access, faces and edges of a given element.

        /**
         * @brief Returns the index of m_uiE2NMapping (CG or DG) for a specified face.
         * @param [in] elementID: element ID of the face, that the face belongs to .
         * @param [in] face : Face that you need the indexing,
         * @param [in] isInternal: return only the internal nodes if true hence the index size would be ((m_uiElementOrder-1)*(m_uiElementOrder-1))
         * @param [out] index: returns the indecies of the requested face.
         * */

        inline void faceNodesIndex(unsigned int elementID, unsigned int face, std::vector<unsigned int>& index,
                                   bool isInternal) const ;

        /**
         * @brief Returns the index of m_uiE2NMapping (CG or DG) for a specified Edge.
         * @param [in] elementID: element ID of the edge, that the face belongs to .
         * @param [in] face1 : one of the face, that an edge belongs to
         * @param [in] face2: second face that an edge belongs to . Hence the edge in cosideration will be intersection of the two faces, face1 and face2.
         * @param [in] isInternal: return only the internal nodes if true hence the index size would be ((m_uiElementOrder-1)*(m_uiElementOrder-1))
         * @param [out] index: returns the indecies of the requested face.
         * */

        inline void edgeNodeIndex(unsigned int elementID,unsigned int face1,unsigned int face2, std::vector<unsigned int >& index, bool isInternal) const ;


        /**
        * @brief Returns the index of m_uiE2NMapping (CG or DG) for a specified Edge.
        * @param [in] elementID: element ID of the edge, that the face belongs to .
        * @param [in] mortonIndex: morton ID of the corner node in the cordinate change in the order of x y z.
        * @param [out] index: returns the indecies of the requested face.
        * */

        inline void cornerNodeIndex(unsigned int elementID,unsigned int mortonIndex, unsigned int &index) const ;

        /**
         * @brief Returns the all the internal node indices of an element.
         * @param [in] elementID: element ID of the edge, that the face belongs to .
         * @param [in] isInternal: return only the internal nodes if true hence the index size would be ((m_uiElementOrder-1)*(m_uiElementOrder-1))
         * @param [out] index: returns the indecies of the requested face.
         *
         *
         * */

        inline void elementNodeIndex(unsigned int elementID,std::vector<unsigned int > & index, bool isInternal) const ;

        /**
         * @brief: Returns true or false (bases on specified edge is hanging or not.)for a given element id , and edge id.
         * @param[in] elementId: element ID of the octant.
         * @param[in] edgeId: edge id
         * */
         bool isEdgeHanging(unsigned int elementId,unsigned int edgeId,unsigned int & cnum) const ;

         /**
          * @brief: Returns true or false (bases on specified edge is hanging or not.)for a given element id , and edge id.
          * @param[in] elementId: element ID of the octant.
          * @param[in] faceId: face id
          * */
         bool isFaceHanging(unsigned int elementId,unsigned int faceId,unsigned int & cnum)const ;

        /**
         * @brief: Returns true if the specified node (e,i,j,k) is hanging.
         * @param[in] eleID: element ID
         * @param[in] ix: i-index of the node.
         * @param[in] jy: j-index of the node.
         * @param[in] kz: k-index of the node.
         * */
        bool isNodeHanging(unsigned int eleID,unsigned int ix,unsigned int jy,unsigned int kz) const ;

         /**
          * @brief: Returns true if the specified node (e,i,j,k) is local.
          * @param[in] eleID: element ID
          * @param[in] ix: i-index of the node.
          * @param[in] jy: j-index of the node.
          * @param[in] kz: k-index of the node.
          * */
        inline bool isNodeLocal(unsigned int eleID,unsigned int ix,unsigned int jy,unsigned int kz) const {
             return ((m_uiE2NMapping_CG[eleID*m_uiNpE+kz*(m_uiElementOrder+1)*(m_uiElementOrder+1)+jy*(m_uiElementOrder+1)+ix]>=m_uiNodeLocalBegin) && (m_uiE2NMapping_CG[eleID*m_uiNpE+kz*(m_uiElementOrder+1)*(m_uiElementOrder+1)+jy*(m_uiElementOrder+1)+ix]<m_uiNodeLocalEnd));
         };

        inline DendroIntL getGhostExcgTotalSendNodeCount() const
        {
            if(m_uiGlobalNpes==1)
                return 0;

            if(m_uiIsActive)
             return (m_uiSendNodeOffset[m_uiActiveNpes-1]+m_uiSendNodeCount[m_uiActiveNpes-1]);
            else
              return 0;
        }

         inline DendroIntL getGhostExcgTotalRecvNodeCount() const
         {
             if(m_uiGlobalNpes==1)
                 return 0;
             
             if(m_uiIsActive)
                 return (m_uiRecvNodeOffset[m_uiActiveNpes-1]+m_uiRecvNodeCount[m_uiActiveNpes-1]);
             else
                 return 0;
         }

        // Methods needed for PDE & ODE and other solvers.
        /**@brief allocate memory for variable array based on the adaptive mesh*/
        template <typename T>
        T* createVector() const;

        /**@brief allocate memory for variable array based on the adaptive mesh.*/
        template <typename T>
        void createVector(std::vector<T> &vec)const ;

        /**@brief allocate memory for variable array based on the adaptive mesh.
          * @param[in] initValue: initialize the vector to the given value.
          * */
        template <typename T>
        T* createVector(const T initValue) const;


         /**@brief allocate memory for variable array based on the adaptive mesh.
         * @param[in] initValue: initialize the vector to the given value.
         * */
        template <typename T>
        T* createVector(std::function<T(T,T,T)> func) const;

         /**@brief allocate memory for variable array based on the adaptive mesh.
          * @param[in] vec: allocate memory for vec
          * @param[in] initValue: initialize the vector to the given value.
          * */
        template <typename T>
        void createVector(std::vector<T> &vec,const T initValue)const ;

        /**
         * @brief create and initialize local elements based on a given function.
         * @param[in] vec: mesh varaible vector
         * @param[in] func: function to initialize with.
         * */
        template  <typename T>
        void createVector(std::vector<T> & vec, std::function<T(T,T,T)> func) const;

        /**
         * @brief: creates (memory allocation) for the unzipped version of the vector.
         * @param[in]: vector to allocate memory for unzipped version
         * @param[out]: vector allocated the memory for unzipped version.
         * */
        template <typename T>
        void createUnZippedVector(std::vector<T> & uvec) const ;

         /**
           * @brief: creates (memory allocation) for the unzipped version of the vector.
           * @param[in]: vector to allocate memory for unzipped version
           * @param[out]: vector allocated the memory for unzipped version.
         * */
        template <typename T>
        T* createUnZippedVector()const ;

         /**@brief allocate memory for variable array based on the adaptive mesh.
          * @param[in] uvec: allocate memory for uvec (unzipped version)
          * @param[in] initValue: initialize the vector to the given value.
          * */
        template <typename T>
        void createUnZippedVector(std::vector<T> & uvec,const T initValue) const ;

        /**@brief allocate memory for variable array based on the adaptive mesh.
          * @param[in] uvec: allocate memory for uvec (unzipped version)
          * @param[in] initValue: initialize the vector to the given value.
          * */
        template <typename T>
        T* createUnZippedVector(const T initValue)const ;

        /**
         *
         * @brief Performs all parent to child interpolations for the m_uiEl_i element in order to apply the stencil.
         * @param[in] parent: function values
         * @param[in] cnum: child number to interpolate.
         * @param[in] dim: dim of the interpolation. (dim=1 for edge interpolation , dim=2 for face interpolation, dim=3 for octant to child interpolation.)
         * @param[out] out: interpolated values.
         *
         * */
        inline void parent2ChildInterpolation(const double * in ,double * out, unsigned int cnum,unsigned int dim=3) const ;

        /**
         * @brief performs the child to parent contribution (only from a single child).
         * @param[in] in: child function values
         * @param[in] cnum: morton ID of the current child.
         * @param[in] dim: dim of the interpolation. (dim=1 for edge interpolation , dim=2 for face interpolation, dim=3 for octant to child interpolation.)
         * @param[out] out: interpolated values. (child to parent contribution)
         *
         * */

        inline void child2ParentInterpolation(const double * in ,double * out, unsigned int cnum,unsigned int dim=3) const ;

        /**
         * @brief Performs child to parent injection.
         * @param [in] children function values. (all the function values of children ordered according to the SFC ordering)
         * @param [in] isHanging  array of boolean variables specifying each node is hanging ot not.
         * @param [out] injected values back to parent.
         *
         * */
        inline void child2ParentInjection(double * in, double * out, bool* isHanging) const;

         /**
           * @author Milinda Fernando
           * @brief Creates the decomposition of adaptive octree variables into blocklist variables that we computed.
           * @param [in] zippedVec adaptive representation of the variable array. (created by createVec function)
           * @param [out] unzippedVec decomposed representation of the adaptive array.
           * @note this routine assumes that for both arrays memory has been allocated. Routine is responsible only to fill up the unzipped entries.
           * */
         template <typename T>
         void unzip(const T* zippedVec, T* unzippedVec);


         /**
           * @brief perform the unzip operation overlap with the communication. GHOST exchage done inside. User doesn't need to perofrm ghost exchange.
           * @param [in] zippedVec adaptive representation of the variable array. (created by createVec function)
           * @param [out] unzippedVec decomposed representation of the adaptive array.
           * @note this routine assumes that for both arrays memory has been allocated (No prior ghost exchange required). Routine is responsible only to fill up the unzipped entries.
           * */
         template <typename T>
         void unzip_async(T* zippedVec, T* unzippedVec,MPI_Request * send_reqs, MPI_Request* recv_reqs, MPI_Status* send_sts, MPI_Status * recv_sts);


         /**@author Milinda Fernando
          * @brief Performs the compression frrom regular block grid varable list to adaptive representation.
          * @param [in] unzippedVec decomposed version of the adaptive array
          * @param [out] compressed version of the unzippedVec.
          * */
         template<typename T>
         void zip(const T* unzippedVec, T* zippedVec);


        /**
         * @brief Apply a given stencil to for provided variable array.
         * @param [in] in : vector that we need to apply the stencil on.
         * @param [in] centered: Stencil that we need to apply on vector in, this is the centered stencil.
         * @param [in] backward: backward version of the centered stencil.(This is used for boundary elements. )
         * @param [in] forward: foward version of the centered stencil. (This is used for boundary elements. )
         * @param [out] out: output vector that after applying the stencil.
         */

        template <typename T,unsigned int length,unsigned int offsetCentered,unsigned int offsetBackward,unsigned int offsetForward>
        void applyStencil(const std::vector<T> &in, std::vector<T> &out,
                          const Stencil<T, length, offsetCentered> &centered,
                          const Stencil<T, length, offsetBackward> &backward,
                          const Stencil<T, length, offsetForward> &forward) ;

        /**
         * @brief Perform the ghost exchange for the vector vec.
         * @param [in] vec: adaptive mesh vector contaiting the values.
         * */
        template<typename T>
        void performGhostExchange(std::vector<T>& vec) ;

         /**
         * @brief Perform the ghost exchange for the vector vec.
         * @param [in] vec: adaptive mesh vector contaiting the values.
         * */
        template<typename T>
        void performGhostExchange(T* vec) ;

         /**
            * @brief Perform the ghost asynchronous send for the vector vec Note: this is a non-blocking asynchronous communication.
            * User is resposible to call the (synchronous call such as MPI_WaitAll) when overlapping the communication and computation.
            * @param [in] vec: adaptive mesh vector contaiting the values.
            * */
         template<typename T>
         void ghostExchangeStart(T* vec,T* sendNodeBuffer,T* recvNodeBuffer, MPI_Request * send_reqs, MPI_Request * recv_reqs) ;

         /**
         * @brief Perform the wait on the recv requests
         * @param [in] vec: adaptive mesh vector containing the values.
         * @param [in] recv_reqs: m_uiRecvProcList.size() recv request
         * @param [in] recv_sts: m_uiRecvProcList.size() recv status
         * */
         template<typename T>
         void ghostExchangeRecvSync(T *vec, T *recvNodeBuffer,MPI_Request *recv_reqs, MPI_Status *recv_sts) ;

         /**
         * @brief Perform the wait on the recv requests
         * @param [in] vec: adaptive mesh vector containing the values.
         * @param [in] send_reqs: m_uiSendProcList.size() send request
         * @param [in] send_sts: m_uiSendProcList.size() send status
         * */
         inline void ghostExchangeSendSync(MPI_Request *send_reqs, MPI_Status *send_sts)
         {
             MPI_Waitall(m_uiSendProcList.size(),send_reqs,send_sts);
         }


        /**
         * @brief write out function values to a vtk file.
         * @param[in] vec: variable vector that needs to be written as a vtk file.
         * @param[in] fprefix: prefix of the output vtk file name.
         * */
        template<typename T>
        void vectorToVTK(const std::vector<T>& vec, char * fprefix,double pTime=0.0, unsigned int nCycle=0) const;



        /**
             * @brief : determine whether any refinement or coarsening need for a specified set of elements. (This uses zipped version of the varibles which needs to satiesfy some constraints. Not every
             * element is eligible for refinemenet or coarsening).
             * @param[in] vec: sequence of varaibles to check
             * @param[in] varIds: variable ids to check. (var ids to index the vec, vec[i] is a T* pointintg to one of the variable in vec. )
             * @param[in] numVars: number of variables to check
             * @param[in] tol: wavelet tolerance
             * Returns true if specified variable violates the specified wavelet toerlance.
             * @note: this method will flag every element in the mesh with OCT_NO_CHANGE, OCT_SPLIT, OCT_COARSE.
             *
             * */
        template <typename T>
        bool isReMesh(const T **vec,const unsigned int * varIds,const unsigned int numVars,double tol, double amr_coarse_fac=DENDRO_AMR_COARSEN_FAC);

       /**
         * @param[in] vec: sequence of varaibles to check
         * @param[in] varIds: variable ids to check. (var ids to index the vec, vec[i] is a T* pointintg to one of the variable in vec. )
         * @param[in] numVars: number of variables to check
         * @param[in] tol: wavelet tolerance
                 * Returns true if specified variable violates the specified wavelet toerlance.
         * @note: this method will flag every element in the mesh with OCT_NO_CHANGE, OCT_SPLIT, OCT_COARSE.
         *
        * */

         template <typename T>
         bool isReMeshUnzip(const T **unzippedVec,const unsigned int * varIds,const unsigned int numVars,std::function<double(double,double,double)>wavelet_tol,double amr_coarse_fac=DENDRO_AMR_COARSEN_FAC, double coarsen_hx=DENDRO_REMESH_UNZIP_SCALE_FAC);

        /**
         * @brief: Remesh the mesh with the new computed elements.
         * @note assumes that refinedIDs and corasenIDs are sorted. (This is automatically done by the isRemesh Fucntion)
         * @param[in] refinedIDs: element IDs need to be refined. (computed by isReMesh function)
         * @param[in] coarsenIDs: element IDs need to be coarsened. (computed by isReMesh function)
         * @param[in] ld_tol: tolerance value used for flexible partitioning
         * @param[in] sfK: spliiter fix parameter (need to specify larger value when run in super large scale)
         * */
        ot::Mesh* ReMesh(unsigned int grainSz=DENDRO_DEFAULT_GRAIN_SZ,double ld_tol=DENDRO_DEFAULT_LB_TOL,unsigned int sfK=DENDRO_DEFAULT_SF_K);

         /**
         * @brief transfer a variable vector form old grid to new grid.
         * @param[in] vec: variable vector needs to be transfered.
         * @param[out] vec: transfered varaible vector
         * @param[in] pMesh: Mesh that we need to transfer the old varaible.
         * */
        template<typename T>
        void interGridTransfer(std::vector<T> & vec,const ot::Mesh* pMesh);

         /**
         * @brief transfer a variable vector form old grid to new grid.
         * @param[in] vec: variable vector needs to be transfered.
         * @param[out] vec: transfered varaible vector
         * @param[in] pMesh: Mesh that we need to transfer the old varaible.
         * */
         template<typename T>
         void interGridTransfer(T*& vec,const ot::Mesh* pMesh);



        /**
         *@brief : Returns the nodal values of a given element for a given variable vector.
         *@param[in] vec: variable vector that we want to get the nodal values.
         *@param[in] elementID: element ID that we need to get the nodal values.
         *@param[out] nodalValues: nodal values of the specified element ID
         *
         * */
        template<typename T>
        void getElementNodalValues(const T* vec,T* nodalValues,unsigned int elementID) const;

        /**
         * @assumption: input is the elemental nodal values.
         * @brief: Computes the contribution of elemental nodal values to the parent elements if it is hanging.
         * Note: internal nodes for the elements cannnot be hagging. Only the face edge nodes are possible for hanging.
         *
         * @param[in] vec: child var vector (nPe)
         * @param[in] elementID: element ID of the current element (or child octant)
         * @param[out] out: add the contributions to the current vector accordingly.
         *
         * Usage: This is needed when performing matrix-free matvec for FEM method.
         *
         * */
        template<typename T>
        void computeElementalContribution(const T* in,T* out,unsigned int elementID) const;



    };


    template<>
    inline void Mesh::init<WaveletDA::LoopType ::ALL>(){
        m_uiEL_i=m_uiElementPreGhostBegin;
    }


    template <>
    inline void Mesh::init<WaveletDA::INDEPENDENT>(){
        m_uiEL_i=m_uiElementLocalBegin;
    }

    template <>
    inline void Mesh::init<WaveletDA::DEPENDENT>() {
        m_uiEL_i=m_uiElementPreGhostBegin;
    }

    template <>
    inline bool Mesh::nextAvailable<WaveletDA::ALL>()
    {
        return (m_uiEL_i<m_uiElementPostGhostEnd);
    }


    template <>
    inline bool Mesh::nextAvailable<WaveletDA::INDEPENDENT>()
    {
        return (m_uiEL_i<m_uiElementLocalEnd);
    }
    template <>

    inline bool Mesh::nextAvailable<WaveletDA::DEPENDENT>()
    {
        return (m_uiEL_i<m_uiElementPreGhostEnd) || ((m_uiEL_i>m_uiElementLocalEnd) && m_uiEL_i<m_uiElementPostGhostEnd);
    }

    template <>
    inline void Mesh::next<WaveletDA::ALL>()
    {
        // interpolation should come here.
        m_uiEL_i++;
    }

    template <>
    inline void Mesh::next<WaveletDA::INDEPENDENT>()
    {
        // interpolation should come here.
        m_uiEL_i++;
    }

    template <>
    inline void Mesh::next<WaveletDA::DEPENDENT>()
    {
        // interpolation should come here.
        m_uiEL_i++;
        if(m_uiEL_i==m_uiElementPreGhostEnd)
            m_uiEL_i=m_uiElementPostGhostBegin;
    }


    inline const ot::TreeNode & Mesh::currentOctant()
    {
        return m_uiAllElements[m_uiEL_i];
    }

    inline unsigned int Mesh::currentIndex() {return m_uiEL_i;}

    inline void Mesh::currentElementNeighbourIndexList(unsigned int *neighList)
    {
        if(!m_uiIsActive) return;

        for(unsigned int k=0;k<m_uiNumDirections;k++)
            neighList[k]=m_uiE2EMapping[m_uiEL_i*m_uiNumDirections+k];
    }


    inline void Mesh::currentElementNodeList(unsigned int * nodeList)
    {
        if(!m_uiIsActive) return;

        for(unsigned int k=0;k<m_uiNpE;k++)
        {
            nodeList[k]=m_uiE2NMapping_CG[m_uiEL_i*m_uiNpE+k];
        }
    }

    inline void Mesh::currentElementNodeList_DG(unsigned int * nodeList)
    {
        if(!m_uiIsActive) return;

        for(unsigned int k=0;k<m_uiNpE;k++)
        {
            nodeList[k]=m_uiE2NMapping_DG[m_uiEL_i*m_uiNpE+k];
        }

    }




    inline void Mesh::parent2ChildInterpolation(const double * in ,double * out,unsigned int cnum,unsigned int dim) const
    {

        dendro::timer::t_unzip_p2c.start();
        if(dim==3)
            m_uiRefEl.I3D_Parent2Child(in,out,cnum);
        else if (dim==2)
            m_uiRefEl.I2D_Parent2Child(in,out,cnum);
        else if (dim==1)
            m_uiRefEl.I1D_Parent2Child(in,out,cnum);
        dendro::timer::t_unzip_p2c.stop();
    }


    inline void Mesh::child2ParentInterpolation(const double * in ,double * out, unsigned int cnum,unsigned int dim) const
    {

        if(dim==3)
            m_uiRefEl.I3D_Child2Parent(in,out,cnum);
        else if (dim==2)
            m_uiRefEl.I2D_Child2Parent(in,out,cnum);
        else if (dim==1)
            m_uiRefEl.I1D_Child2Parent(in,out,cnum);

    }


    inline void Mesh::child2ParentInjection(double * in, double * out, bool * isHanging) const
    {

        assert(m_uiElementOrder%2==0 || m_uiElementOrder==1);
        if(m_uiElementOrder==1)
        { // special case
            for(unsigned int child=0;child<NUM_CHILDREN;child++)
            {
                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                            out[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=in[child*m_uiNpE+k*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
            }


        }else
        {

            for(unsigned int child=0;child<NUM_CHILDREN;child++)
            {

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                        {
                            if((i%2==0) && (j%2==0) && (k%2==0) && isHanging[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i])
                            {
                                out[((((child & 4u)>>2u)*m_uiElementOrder+k)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((child & 2u)>>1u)*m_uiElementOrder+j)>>1)*(m_uiElementOrder+1)+(((((child & (1u)))*m_uiElementOrder+i)>>1))]=in[child*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                //if(!m_uiActiveRank) std::cout<<"rank: "<<m_uiActiveRank<<" child: "<<child<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" inserted to node : "<<((((child & 4u)>>2u)*m_uiElementOrder+k)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((child & 2u)>>1u)*m_uiElementOrder+j)>>1)*(m_uiElementOrder+1)+(((((child & (1u)))*m_uiElementOrder+i)>>1))<<std::endl;
                            }


                        }

            }



        }


    }


    inline bool Mesh::computeOveralppingNodes(const ot::TreeNode & parent, const ot::TreeNode & child, int * idx,int *idy, int * idz)
    {

        unsigned int  Lp =1u<<(m_uiMaxDepth-parent.getLevel());
        unsigned int  Lc =1u<<(m_uiMaxDepth-child.getLevel());
        //intilize the mapping to -1. -1 denotes that mapping is not defined for given k value.

        unsigned int dp,dc;
        dp=(m_uiElementOrder);
        dc=m_uiElementOrder;

        assert(Lp%dp==0);
        assert(Lc%dc==0);


        for (unsigned int k = 0; k < (m_uiElementOrder + 1); k++) {
            idx[k] = -1;
            idy[k] = -1;
            idz[k] = -1;

        }
        bool state=false;
        bool stateX=false;
        bool stateY=false;
        bool stateZ=false;
        if(parent==child) {
            for (unsigned int k = 0; k < (m_uiElementOrder + 1); k++) {
                idx[k] = k;
                idy[k] = k;
                idz[k] = k;

            }
            return true;
        }else if(parent.isAncestor(child))
        {

            /*if((((child.getX()-parent.getX())*m_uiElementOrder)%Lp) || (((child.getY()-parent.getY())*m_uiElementOrder)%Lp) || (((child.getZ()-parent.getZ())*m_uiElementOrder)%Lp)) return false;
            else*/
            {
                unsigned int index[3];
                for (unsigned int k = 0; k < (m_uiElementOrder + 1); k++) {

                    index[0]=(m_uiElementOrder+1);
                    index[1]=(m_uiElementOrder+1);
                    index[2]=(m_uiElementOrder+1);

                    if(!(((child.getX()-parent.getX())*dp*dc + k*Lc*dp)%(Lp*dc)))
                        index[0] = ((child.getX()-parent.getX())*dp*dc + k*Lc*dp)/(Lp*dc);//((child.getX()-parent.getX())*m_uiElementOrder + k*Lc)/Lp;

                    if(!(((child.getY()-parent.getY())*dp*dc + k*Lc*dp)%(Lp*dc)))
                        index[1] = ((child.getY()-parent.getY())*dp*dc + k*Lc*dp)/(Lp*dc);//((child.getY()-parent.getY())*m_uiElementOrder + k*Lc)/Lp;

                    if(!(((child.getZ()-parent.getZ())*dp*dc + k*Lc*dp)%(Lp*dc)))
                        index[2] = ((child.getZ()-parent.getZ())*dp*dc + k*Lc*dp)/(Lp*dc);//((child.getZ()-parent.getZ())*m_uiElementOrder + k*Lc)/Lp;

                    if(!stateX && index[0]<(m_uiElementOrder+1))
                        stateX=true;

                    if(!stateY && index[1]<(m_uiElementOrder+1))
                        stateY=true;

                    if(!stateZ && index[2]<(m_uiElementOrder+1))
                        stateZ=true;


                    if(index[0]<(m_uiElementOrder+1)) {idx[k]=index[0]; assert((parent.getX()+idx[k]*Lp/dp)==(child.getX()+k*Lc/dc));}
                    if(index[1]<(m_uiElementOrder+1)) {idy[k]=index[1]; assert((parent.getY()+idy[k]*Lp/dp)==(child.getY()+k*Lc/dc));}
                    if(index[2]<(m_uiElementOrder+1)) {idz[k]=index[2]; assert((parent.getZ()+idz[k]*Lp/dp)==(child.getZ()+k*Lc/dc));}

                }
                state=stateX & stateY & stateZ;
                return state;

            }

        }else
        {
            return false;
        }




    }







}




#include "mesh.tcc"
#include "meshE2NUtils.tcc"



#endif //SFCSORTBENCH_MESH_H
