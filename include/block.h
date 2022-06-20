//
// Created by milinda on 4/21/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief contains the class for the block. In order to perform a stencil on the adaptive grid, we treat the adaptive grid as a collection of regular blocks.
 *  Let \f$ \tau \f$ be a 2:1 balance complete sorted (according SFC ordering) octree then we can decompose \f$\tau =\{b_i^l\}\f$   sequence of finite number of regular
 *  blocks.
 *
 *  This class contains the block coordinates level of the regular grid embedded, stencil ghost width and number of local points available.
*/

#ifndef SFCSORTBENCH_BLOCK_H
#define SFCSORTBENCH_BLOCK_H

#include "dendro.h"
#include "TreeNode.h"
#include <assert.h>
#include <treenode2vtk.h>


namespace ot
{
   /**
    * @brief Block type 
    * UNSPECIFIED : block flag type is not set
    * UNZIP_INDEPENDENT: unzip operation does not depend on ghost nodes
    * UNZIP_DEPENDENT: unzip operation does depend on at least one ghost node
    * 
    */
   enum BlockType{UNSPECIFIED=0, UNZIP_INDEPENDENT, UNZIP_DEPENDENT};

   class Block
   {



   private:
     /**Coordinates of the block. */
     ot::TreeNode m_uiBlockNode;

     /**rotation id of the block*/
     unsigned int m_uiRotID;

     /** size of the regular grid inside the block. */
     unsigned int m_uiRegGridLev;

      /**regular grid local element begin*/
     DendroIntL m_uiLocalElementBegin;

     /** regular grid local element end. Note that element ids in [localBegin,localEnd] is continous and all those elements are inside the current block. */
     DendroIntL m_uiLocalElementEnd;

     /** offset used for local memory allocation.*/
     DendroIntL m_uiOffset;

     /** array size (1D for the current block. */
     unsigned int m_uiSize1D;

     /**padding width (1D) used for pad the block for neighbour blocks. */
     unsigned int m_uiPaddingWidth;

     /**element order */
     unsigned int m_uiEleOrder;

     /** allocation length on X direction*/
     unsigned int m_uiSzX;

     /** allocation length on Y direction*/
     unsigned int m_uiSzY;

     /** allocation length on Z direction*/
     unsigned int m_uiSzZ;

     /** indecies of the 12 negihbour elems*/
     std::vector<unsigned int> m_uiBLK2DIAG;

     /** indecies of the 8 vertex neighbor elems.*/
     std::vector<unsigned int> m_uiBLKVERTX;

     /** number of elements per block. **/
     unsigned int m_uiBlkElem_1D;

     /** set true after the perform block setpup if the block doesn't depend on the ghost region*/
     bool m_uiIsInternal;
     
     /** block type*/      
     BlockType m_uiBlkType;


   public:
    /**@brief Default constructor*/
    Block();

    /**
     * @brief constructor to initialize and create a block.
     * @param [in] pNode ot::TreeNode for the block.
     * @param [in] rotID rotation ID for the block.
     * @param [in] regLev level of the regular grid embedded by the block.
     * @param [in] regEleBegin Local element begin location for the for the octree embedded by the block.
     * @param [in] regEleEnd Local element end location for the octree embedded by the block .
     * @param [in] eleorder: element order of the mesh.
     * */
    Block(ot::TreeNode pNode, unsigned int rotID ,unsigned int regLev, unsigned int regEleBegin,unsigned int regEleEnd,unsigned int  eleOrder);

    ~Block();

    /**
     * @brief Return the block node
     * */
    inline ot::TreeNode getBlockNode()const {return m_uiBlockNode;}

    /**
     * @brief returns the regular grid lev (m_uiRegGridLev) value.
     * note: In octree2BlockDecomposition m_uiRegGridLev is used to store the rotation id of the block.
     *  */
     inline unsigned int getRegularGridLev()const {return m_uiRegGridLev;}

     /**@brief returns the rotation id of the block*/
     inline unsigned int getRotationID()const {return m_uiRotID;}

     /**
      * @brief returns the local element begin for the block.
      * */
      inline DendroIntL getLocalElementBegin()const {return m_uiLocalElementBegin;}

     /**
      * @brief returns the local element end for the block.
      * */

     inline DendroIntL getLocalElementEnd()const {return m_uiLocalElementEnd;}

     /** @brief returns 1D padding width */
     inline unsigned int get1DPadWidth()const {return m_uiPaddingWidth;}

     /**@brief returns the element order*/
     inline unsigned int getElementOrder() const {return m_uiEleOrder;}

     /**@brief set the block offset*/
     void setOffset(DendroIntL offset);

     inline void setBlk2DiagMap(unsigned int owner,unsigned int dir, unsigned int id){m_uiBLK2DIAG[dir*(2*m_uiBlkElem_1D)+owner]=id;}

     inline void setBlk2VertexMap(unsigned int dir, unsigned int id){m_uiBLKVERTX[dir]=id;}

     inline void setIsInternal(bool isInternal){m_uiIsInternal=isInternal;}

     inline void setBlkType(BlockType btype) { m_uiBlkType = btype; }

     inline BlockType getBlockType() const {return m_uiBlkType;}

     inline bool isInternal(){return m_uiIsInternal;}

     /**@brief set the blkFlag with the correct bdy*/
     inline void setBlkNodeFlag(unsigned int flag){m_uiBlockNode.setFlag(flag);};

     /**@brief set the blkFlag with the correct bdy*/
     inline unsigned int getBlkNodeFlag() const { return (m_uiBlockNode.getFlag()>>NUM_LEVEL_BITS);};

       /** @brief get offset*/
     inline DendroIntL getOffset() const {return m_uiOffset;}

     /** @brief returns the 1D array size*/
     inline unsigned int get1DArraySize() const {return m_uiSize1D ;}

     /**@brief allocation length on X direction*/
     inline unsigned int getAllocationSzX() const {return m_uiSzX;}

     /**@brief allocation length on Y direction*/
     inline unsigned int getAllocationSzY() const {return m_uiSzY;}

     /**@brief allocation length on Z direction*/
     inline unsigned int getAllocationSzZ() const {return m_uiSzZ;}

     /**@brief align the total block size*/
     inline unsigned int getAlignedBlockSz() const
     {
        // unsigned int tmp;
        // ((m_uiSzX & ((1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG)-1))==0)? tmp=m_uiSzX : tmp=((m_uiSzX/(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG))+1)*(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG);
        return m_uiSzX*m_uiSzY*m_uiSzZ;
        // unsigned int ax=binOp::getNextHighestPowerOfTwo(m_uiSzX);
        // unsigned int ay=binOp::getNextHighestPowerOfTwo(m_uiSzY);
        // unsigned int az=binOp::getNextHighestPowerOfTwo(m_uiSzZ);
        // return ax * ay * az;
     }


     inline const unsigned int * getBlk2DiagMap() const {return &(*(m_uiBLK2DIAG.begin()));}

     inline const unsigned int * getBlk2VertexMap() const {return &(*(m_uiBLKVERTX.begin()));}

     inline void setAllocationSzX(unsigned int sz) {m_uiSzX=sz;}
     inline void setAllocationSzY(unsigned int sz) {m_uiSzY=sz;}
     inline void setAllocationSzZ(unsigned int sz) {m_uiSzZ=sz;}
     inline void setSiz1D(unsigned int sz){m_uiSize1D=sz;}

     inline unsigned int getElemSz1D() const { return m_uiBlkElem_1D;}


     inline const std::vector<unsigned int >& getBlk2DiagMap_vec() const {return m_uiBLK2DIAG;}
     inline const std::vector<unsigned int >& getBlk2VertexMap_vec() const {return m_uiBLKVERTX;}

     /**@brief computes and returns the space discretization (grid domain) */
     double computeGridDx () const;

     /**@brief computes and returns the space discretization (grid domain) */
     double computeGridDy () const;

     /**@brief computes and returns the space discretization (grid domain) */
     double computeGridDz () const;

     /**@brief computes and returns the space discretization on x direction (problem domain)*/
     double computeDx(const Point & d_min,const Point & d_max) const ;
     /**@brief computes and returns the space discretization on x direction (problem domain)*/
     double computeDy(const Point & d_min,const Point & d_max) const ;
     /**@brief computes and returns the space discretization on x direction (problem domain)*/
     double computeDz(const Point & d_min,const Point & d_max) const ;

     /*** @brief initialize the block diagonal map. */
     void initializeBlkDiagMap(const unsigned int value);

     /*** @brief initialize the block vertex neighbour map. */
     void initializeBlkVertexMap(const unsigned int value);

     /**@brief compute the eijk for an element inside the block.  */
     void computeEleIJK(ot::TreeNode pNode, unsigned int* eijk) const;

     /**@brief: returns true if the pNode is inside the current block*/
     bool isBlockInternalEle(ot::TreeNode pNode) const ; 




   };

} // end of namespace ot





#endif //SFCSORTBENCH_BLOCK_H
