//
// Created by milinda on 2/8/16.
//

/**
 * @author Milinda Fernando
 * @author Hari Sundar
 * @author Rahul Sampath
 * @breif Refactored version of the TreeNode class with the minimum data required.
 * @remarks m_uiMaxDepth is removed from the arttributre list and listed as a global variable.
 * */

#ifndef SFCSORTBENCH_TREENODE_H
#define SFCSORTBENCH_TREENODE_H

static unsigned int MAX_LEVEL=31;
//#define m_uiMaxDepth 30
extern unsigned int m_uiMaxDepth;


#include "hcurvedata.h"
#include "binUtils.h"
#include "mpi.h"
#include "point.h"
#include <algorithm>
#include "dendro.h"
#include <ostream>


namespace ot {

    /**
      @brief A class to manage octants.
      @author Rahul Sampath
      */
    class TreeNode {
    protected:
        //Level is also used as a flag.
        unsigned int m_uiX, m_uiY, m_uiZ, m_uiLevel;
        //unsigned int m_uiMaxDepth;
        unsigned int m_uiWeight=1;

    public:

        TreeNode (const int dummy, const unsigned int x,const unsigned int y,
                  const unsigned int z,const unsigned int level, const unsigned int dim,
                  const unsigned int maxDepth);

        TreeNode (const unsigned int x,const unsigned int y,
                  const unsigned int z,const unsigned int level, const unsigned int dim,
                  const unsigned int maxDepth);

        TreeNode(const unsigned int dim, const unsigned int maxDepth );

        TreeNode();

    public:

        /**@brief: getweight of the octant (used in the wpart)*/
        unsigned int getWeight() const {return m_uiWeight;}

        /**@brief: set weight of the octant  */
        void setWeight(unsigned int w) {m_uiWeight=w;}

        /**
          @author Milinda Fernando
          @brief return the dimension of the octant
          */
         unsigned int getDim() const;
        /**
          @author Milinda Fernando
          @brief return the maxDepth of the octant based of the global value.
          */
        unsigned int getMaxDepth() const;

        /**
          @author Milinda Fernando
          @brief return the level of the  of the octant
          */
         unsigned int getLevel() const;

        /**
          @author Milinda Fernando
          @brief set the m_uiFlag value. Which is used to store the level and other aditional info.
          */

         inline void setFlag(unsigned int flag);

        /**
          @author Milinda Fernando
          @brief get  the m_uiFlag value. Which is used to store the level and other aditional info.
          */
         inline unsigned int getFlag() const;

        /**
         @author Milinda Fernando
         @brief get integer values of the octree coordinates.
         */
         unsigned int getX() const;
         unsigned int getY() const;
         unsigned int getZ() const;

        /**
          @author Milinda Fernando
          @brief returns true if the current considering node is a root node.
          */
         bool isRoot() const;

        /**
         * @author Milinda Fernando
         * @brief Can be used the increment octant level by 1. Used in SFC treeSearch function.
         * */
        void incrementLevel();

        /**
        @author Rahul Sampath
        @brief Two octants are equal if their respective anchors are equal and their levels are
         equal.
        */
        bool  operator == ( TreeNode const  &other) const;

        /**
          @author Rahul Sampath
          @brief Two octants are equal if their respective anchors are equal and their levels are
          equal.
          */
        bool  operator != (TreeNode const  &other) const;


        /**
         @author Rahul Sampath
         @author Milinda Fernando
         @brief The comparisons are based on the Morton/Hilbert
                ordering of the octants
                @remarks Hilbert based comparision operator is O(nlog(n)). Which is costly than the standard Morton implmentation.
         */

        bool  operator < ( TreeNode const  &other) const;

        /**
          @author Rahul Sampath
          @brief The comparisons are based on the Morton/Hilbert
                 ordering of the octants
          */
        bool  operator > ( TreeNode const  &other) const;

        /**
          @author Rahul Sampath
          @brief The comparisons are based on the Morton/Hilbert
                 ordering of the octants
          */
        bool  operator <= ( TreeNode const  &other) const;

        /**
          @author Rahul Sampath
          @brief The comparisons are based on the Morton/Hilbert
                 ordering of the octants
          */
        bool  operator >= ( TreeNode const  &other) const;

        unsigned int genHkey_Bonsai_sc16() const ;

        TreeNode getNCA(TreeNode const & other ) const ;

        inline bool isAncestor(const TreeNode & other) const;


        /**
         @author Milinda Fernando
         @brief returns true min corner(anchor) and the max corner coordinates of a octant.
         */

        inline unsigned int minX() const {
            return getX();
        } //end fn.

        inline unsigned int minY() const {
            if (m_uiDim < 2) { return 0; }
            return getY();
        } //end fn.

        inline unsigned int minZ() const {
            if (m_uiDim < 3) { return 0; }
            return getZ();
        } //end fn.

        inline unsigned int maxX() const {
            unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
            return (minX() + len);
        } //end fn.

        inline unsigned int maxY() const {
            if (m_uiDim < 2) { return 1; }
            unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
            return (minY() + len);
        } //end fn.

        inline unsigned int maxZ() const {
            if (m_uiDim < 3) { return 1; }
            unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
            return (minZ() + len);
        } //end fn.

        /**
         @author Milinda Fernando
         @brief returns the parent of the current octant
         */
        TreeNode  getParent() const;


        /**
         *@breif : Returns the deepest first child in Morton ordering.
         *
         */
        TreeNode getDFDMorton() const;


        /** Returns the deepest first decent  of the considering octant based on the SFC being used. */
        TreeNode getDFD() const;

        /** Returns the deepest last decent of the considering octant based on the SFC beging used*/
        TreeNode getDLD() const;



        /**
        @author Milinda Fernando
        @author Rahul Sampath
        @brief returns the neighbour octants (right, left top bottom etc)
        */

        TreeNode  getLeft() const;
        TreeNode  getRight() const;
        TreeNode  getTop() const;
        TreeNode  getBottom() const;
        TreeNode  getFront() const;
        TreeNode  getBack() const;
        TreeNode  getTopLeft() const;
        TreeNode  getTopRight() const;
        TreeNode  getBottomLeft() const;
        TreeNode  getBottomRight() const;
        TreeNode  getLeftFront() const;
        TreeNode  getRightFront() const;
        TreeNode  getTopFront() const;
        TreeNode  getBottomFront() const;
        TreeNode  getTopLeftFront() const;
        TreeNode  getTopRightFront() const;
        TreeNode  getBottomLeftFront() const;
        TreeNode  getBottomRightFront() const;
        TreeNode  getLeftBack() const;
        TreeNode  getRightBack() const;
        TreeNode  getTopBack() const;
        TreeNode  getBottomBack() const;
        TreeNode  getTopLeftBack() const;
        TreeNode  getTopRightBack() const;
        TreeNode  getBottomLeftBack() const;
        TreeNode  getBottomRightBack() const;
        std::vector<TreeNode> getAllNeighbours() const;


        int getAnchor(unsigned int &x, unsigned int&y, unsigned int&z) const;
        Point getAnchor() const { return Point(m_uiX, m_uiY, m_uiZ); };

        inline unsigned int getMortonIndex() const;


        /**
         * @author Rahul Sampath
         * @author Milinda Fernando
         * @breif Following are the attributes and functionalities that needed by mesh computation based on the 2:1 balanced octree.
         * Most of the functions taken from the old dendro implementation.
         * */

        //==================================================================================================================

        /**
        @brief The type of boundary
        */
        enum BoundaryType1 { NEGATIVE= 2, POSITIVE= 4};

        /**
          @brief The type of boundary
          */
        enum BoundaryType2 {
            X_NEG_BDY=1, Y_NEG_BDY=2, Z_NEG_BDY=4, NEG_POS_DEMARCATION=8,
            EXTERNAL_BDY=16, X_POS_BDY=32, Y_POS_BDY=64, Z_POS_BDY=128
        };

        /**
          @brief  The type of boundary
          */
        enum BoundaryType3 {
            FACE_BDY=1, EDGE_BDY=2, CORNER_BDY=3
        };

        enum OctantFlagType {
            MAX_LEVEL=31, BOUNDARY=64, NODE=128
        };


        /**
        @author Hari Sundar
        @brief flags is a datastructure which will store which boundaries were
        touched. highest 3 bits are for +Z,+y, and +x axes ... and and
        smallest 3 are for -z,-y and -x axes.
        */
        bool isBoundaryOctant(int type=POSITIVE, unsigned char *flags=NULL) const;

        /**
          @author Hari Sundar
          @brief flags is a datastructure which will store which boundaries were
          touched. highest 3 bits are for +Z,+y, and +x axes ... and and
          smallest 3 are for -z,-y and -x axes.
          */
        bool isBoundaryOctant(const TreeNode &block, int type=POSITIVE, unsigned char *flags=NULL) const;


        //=====================================================Mesh Definitions End ==============================================

        int addChildren(std::vector<ot::TreeNode>& children) const;

        int getChildrenInMortonOrdering(std::vector<ot::TreeNode>& children) const
        {
            unsigned int dim = m_uiDim;
            unsigned int maxDepth = m_uiMaxDepth;
            children.resize((1 << dim));

            //#define MORTON_ORDERING

            if ((m_uiLevel & ot::TreeNode::MAX_LEVEL) == maxDepth) {
                for (int i = 0; i < (1 << dim); i++) {
                    children[i] = *this;
                }
                return 1;
            }
            //The check that lev < maxD is taken care of in the constructor.
            //Order: X first, Y next and Z last

            unsigned int len = (unsigned int)(1u << (maxDepth - ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1)));

            TreeNode   tmpNode0(m_uiX, m_uiY, m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[0] = tmpNode0;

            TreeNode   tmpNode1((m_uiX + len), m_uiY, m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[1] = tmpNode1;

            if (dim >= 2) {
                TreeNode   tmpNode2(m_uiX, (m_uiY + len), m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
                children[2] = tmpNode2;

                TreeNode   tmpNode3((m_uiX + len), (m_uiY + len), m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
                children[3] = tmpNode3;
            }

            if (dim == 3) {
                TreeNode   tmpNode4(m_uiX, m_uiY, (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
                children[4] = tmpNode4;

                TreeNode   tmpNode5((m_uiX + len), m_uiY, (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
                children[5] = tmpNode5;

                TreeNode   tmpNode6(m_uiX, (m_uiY + len), (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
                children[6] = tmpNode6;

                TreeNode   tmpNode7((m_uiX + len), (m_uiY + len), (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
                children[7] = tmpNode7;
            } //end if
            return 1;

        }

    };


    inline unsigned int TreeNode::getX() const {
        return m_uiX;
    }

    inline unsigned int TreeNode::getY() const {
        return m_uiY;
    }

    inline unsigned int TreeNode::getZ() const {
        return m_uiZ;
    }


    inline unsigned int TreeNode::getDim() const {
        return m_uiDim;
    }

    inline unsigned int TreeNode::getMaxDepth() const {
        return m_uiMaxDepth;
    }

    inline unsigned int TreeNode::getLevel() const {
        return (m_uiLevel & MAX_LEVEL);
    }

    inline bool TreeNode::isRoot() const {
        return ((this->getLevel()==0) & (m_uiX==0) & (m_uiY==0) & (m_uiZ==0));
    }

    inline void TreeNode::incrementLevel()
    {
        if(this->getLevel()<m_uiMaxDepth)
            m_uiLevel++;
    }

    inline void  TreeNode::setFlag(unsigned int flag) {  m_uiLevel=flag;  }
    inline unsigned int TreeNode::getFlag() const { return m_uiLevel; }

    inline bool TreeNode::isAncestor(const TreeNode & other) const
    {
        /*unsigned int min1[3], min2[3], max1[3], max2[3];

        min1[0] = this->minX(); min1[1] = this->minY(); min1[2] = this->minZ();
        min2[0] = other.minX(); min2[1] = other.minY(); min2[2] = other.minZ();

        max1[0] = this->maxX(); max1[1] = this->maxY(); max1[2] = this->maxZ();
        max2[0] = other.maxX(); max2[1] = other.maxY(); max2[2] = other.maxZ();

        bool state1=( (this->getLevel() < other.getLevel()) && ( (min2[0] >= min1[0]) && (min2[1] >= min1[1]) && (min2[2] >= min1[2]) && (max2[0] <= max1[0]) && (max2[1] <= max1[1]) && (max2[2] <= max1[2]) ));

        return state1;*/

        return ( (this->getLevel() < other.getLevel()) && ( (other.minX() >= this->minX()) && (other.minY() >= this->minY()) && (other.minZ() >= this->minZ()) && (other.maxX() <= this->maxX()) && (other.maxY() <= this->maxY()) && (other.maxZ() <= this->maxZ()) ));

    }


    std::ostream & operator << (std::ostream & os,TreeNode const & node) ;


    inline bool TreeNode::operator<(TreeNode const &other) const {


#ifdef __DEBUG_TN__
        if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
      std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
      std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
      assert(false);
    }
#endif

        if ((this->m_uiX == other.m_uiX) && (this->m_uiY == other.m_uiY) && (this->m_uiZ == other.m_uiZ)) {
            return ((this->m_uiLevel & MAX_LEVEL) < (other.m_uiLevel & MAX_LEVEL));
        } //end if

#ifdef HILBERT_ORDERING

        // #pragma message "Hilbert NCA"
        // If you need to initilize the Hilbert table and the rotations for 2D you need to define DENDRO_DIM2
        // Default initialization for 3D case.
        // NOTE: To work the Hilbert Ordering You need the Hilbert Table Initialized.
        //std::cout<<"Using H Sort:"<<std::endl;
        unsigned int x1 = m_uiX;
        unsigned int x2 = other.getX();

        unsigned int y1 = m_uiY;
        unsigned int y2 = other.getY();

        unsigned int z1 = m_uiZ;
        unsigned int z2 = other.getZ();

        unsigned len;
        unsigned int maxDepth = m_uiMaxDepth;

        if(this->getLevel()>other.getLevel())
        {
            len=1u<<(this->getMaxDepth()-other.getLevel());
            if(!((x1<x2 || x1>=(x2+len)) || (y1<y2 || y1>=(y2+len)) ||(z1<z2 || z1>=(z2+len))))
                return false;


        }else if(this->getLevel()<other.getLevel())
        {
            len=1u<<(this->getMaxDepth()-this->getLevel());
            if(!((x2<x1 || x2>=(x1+len))||(y2<y1 || y2>=(y1+len))||(z2<z1 || z2>=(z1+len))))
                return true;


        }

        unsigned int maxDiff = (unsigned int)(std::max((std::max((x1^x2),(y1^y2))),(z1^z2)));
        int dim=m_uiDim;




        unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
        //Eliminate the last maxDiffBinLen bits.
        unsigned int ncaX = ((x1>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaY = ((y1>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaZ = ((z1>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaLev = (maxDepth - maxDiffBinLen);

//        if(ncaLev>std::min(this->getLevel(),other.getLevel()))
//        {
//
//            std::cout<<"P1:"<<*(this)<<"\t P2:"<<other<<std::endl;
//            std::cout<<"MaxDiff:"<<maxDiff<<std::endl;
//            std::cout<<"MaxDiffLen:"<<maxDiffBinLen<<std::endl;
//            std::cout<<"NCALEV:"<<ncaLev<<"\t this_lev:"<<this->getLevel()<<"\t other_lev:"<<other.getLevel()<<std::endl;
//
//            std::cout<<std::endl;
//
//        }

        unsigned int index1=0;
        unsigned int index2=0;
        unsigned int num_children=1u<<dim; // This is basically the hilbert table offset
        unsigned int rot_offset=num_children<<1;
        //char index_temp=0;
        int current_rot=0;

        //unsigned int b_x,b_y,b_z;
        //unsigned int a,b,c;
        unsigned int mid_bit=m_uiMaxDepth;

        for(int i=0; i<ncaLev;i++)
        {
            mid_bit=m_uiMaxDepth-i-1;

            //b_x=((ncaX&(1<<mid_bit))>>mid_bit);
            //b_y=((ncaY&(1<<mid_bit))>>mid_bit);
            //b_z=((ncaZ&(1<<mid_bit))>>mid_bit);

            // index1=(b_z<<2) + ((b_x^b_z)<<1) + (b_x^b_y^b_z);
            index1= ((((ncaZ & (1u << mid_bit)) >> mid_bit) << 2u) |(((ncaY & (1u << mid_bit)) >> mid_bit) << 1u) | ((ncaX & (1u << mid_bit)) >> mid_bit));
            //index_temp=rotations[rot_offset*current_rot+num_children+index1]-'0';
            current_rot=HILBERT_TABLE[current_rot*num_children+index1];
        }
        mid_bit--;
        index1= ((((z1 & (1u << mid_bit)) >> mid_bit) << 2u) |(((y1 & (1u << mid_bit)) >> mid_bit) << 1u) | ((x1 & (1u << mid_bit)) >> mid_bit));
        index2= ((((z2 & (1u << mid_bit)) >> mid_bit) << 2u) |(((y2 & (1u << mid_bit)) >> mid_bit) << 1u) | ((x2 & (1u << mid_bit)) >> mid_bit));

        return rotations[rot_offset*current_rot+num_children+index1] < rotations[rot_offset*current_rot+num_children+index2];


#else
        //#ifdef USE_NCA_PROPERTY

	      //  #pragma message "Morton NCA"
	      //  return morton_order_NCA(p1,p2);
        // #else
	        #pragma message "Morton"
          // -- original Morton
          // first compare the x, y, and z to determine which one dominates ...
          //Ancestor is smaller.
          if ((this->m_uiX == other.m_uiX) && (this->m_uiY == other.m_uiY) && (this->m_uiZ == other.m_uiZ)) {
            return ((this->m_uiLevel & MAX_LEVEL) < (other.m_uiLevel & MAX_LEVEL));
          } //end if

          unsigned int x = (m_uiX ^ other.m_uiX);
          unsigned int y = (m_uiY ^ other.m_uiY);
          unsigned int z = (m_uiZ ^ other.m_uiZ);

          //Default pref: z > y > x.
          unsigned int maxC = z;
          unsigned int yOrx = y;
          if (yOrx < x) {if ((x ^ yOrx) >= yOrx) {yOrx = x;}
          }
          if (maxC < yOrx) {if ((maxC ^ yOrx) >= maxC) {maxC = yOrx;}
          }

          if (maxC == z) {return (m_uiZ < other.m_uiZ); } else if (maxC == y) {return (m_uiY < other.m_uiY); } else {return (m_uiX < other.m_uiX); }
        // -- original Morton

    // #endif
#endif

    } //end function

    inline bool TreeNode::operator<=(TreeNode const &other) const {
#ifdef __DEBUG_TN__
        if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
      std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
      std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
      assert(false);
    }
#endif
        return (((*this) < other) || ((*this) == other));
    } //end fn.

    inline bool TreeNode::operator>(TreeNode const &other) const {
#ifdef __DEBUG_TN__
        if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
      std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
      std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
      assert(false);
    }
#endif
        return ((!((*this) < other)) && (!((*this) == other)));
    } //end fn.

    inline bool TreeNode::operator>=(TreeNode const &other) const {
#ifdef __DEBUG_TN__
        if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
      std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
      std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
      assert(false);
    }
#endif
        return (!((*this) < other));
    } //end fn.

    inline bool TreeNode::operator==(TreeNode const &other) const {
#ifdef __DEBUG_TN__
        if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
      std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
      std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
      assert(false);
    }
#endif
        if ((m_uiX == other.m_uiX) && (m_uiY == other.m_uiY) && (m_uiZ == other.m_uiZ) &&
            ((m_uiLevel & MAX_LEVEL) == (other.m_uiLevel & MAX_LEVEL))) {
            return true;
        } else {
            return false;
        }
    } //end fn.

    inline bool TreeNode::operator!=(TreeNode const &other) const {
#ifdef __DEBUG_TN__
        if (((this->m_uiDim) != (other.m_uiDim)) || ((this->m_uiMaxDepth) != (other.m_uiMaxDepth))) {
      std::cout << "Me: " << (*this) << " Other: " << other << std::endl;
      std::cout << "My Dim: " << m_uiDim << " OthDim: " << other.m_uiDim << " My MaxD: " << m_uiMaxDepth << " othMD: " << other.m_uiMaxDepth << std::endl;
      assert(false);
    }
#endif
        return (!((*this) == other));
    } //end fn.


    inline TreeNode TreeNode::getParent() const {
        //For any node at level l, the last (maxD-l) bits are 0.
        //By convention, root's parent is also root.
        unsigned int parX, parY, parZ;
        unsigned int parLev = (((m_uiLevel & MAX_LEVEL) > 0) ? ((m_uiLevel & MAX_LEVEL) - 1) : 0);
        parX = ((m_uiX >> (m_uiMaxDepth - parLev)) << (m_uiMaxDepth - parLev));
        parY = ((m_uiY >> (m_uiMaxDepth - parLev)) << (m_uiMaxDepth - parLev));
        parZ = ((m_uiZ >> (m_uiMaxDepth - parLev)) << (m_uiMaxDepth - parLev));
        return TreeNode(1, parX, parY, parZ, parLev, m_uiDim, m_uiMaxDepth);
    } //end function


    inline TreeNode   TreeNode::getLeft() const {
        //-ve X
        if (minX() > 0) {
            unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
            unsigned int xres = (minX() - len);
            unsigned int yres = minY();
            unsigned int zres = minZ();
            TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
            return res;
        } else {
            TreeNode res(m_uiDim, m_uiMaxDepth);
            return res;
        }
    } //end fn.

    inline TreeNode   TreeNode::getRight() const {
        //+ve X
        if (maxX() < (1u << m_uiMaxDepth)) {
            unsigned int xres = maxX();
            unsigned int yres = minY();
            unsigned int zres = minZ();
            TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
            return res;
        } else {
            TreeNode res(m_uiDim, m_uiMaxDepth);
            return res;
        }
    } //end fn.

    inline TreeNode   TreeNode::getBack() const {
        //+ve Y
        if ((m_uiDim > 1) && (maxY() < (1u << m_uiMaxDepth))) {
            unsigned int xres = minX();
            unsigned int yres = maxY();
            unsigned int zres = minZ();
            TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
            return res;
        } else {
            TreeNode res(m_uiDim, m_uiMaxDepth);
            return res;
        }
    } //end fn.

    inline TreeNode   TreeNode::getFront() const {
        //-ve Y
        if (minY() > 0) {
            unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
            unsigned int xres = minX();
            unsigned int yres = minY() - len;
            unsigned int zres = minZ();
            TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
            return res;
        } else {
            TreeNode res(m_uiDim, m_uiMaxDepth);
            return res;
        }
    } //end fn.

    inline TreeNode   TreeNode::getBottom() const {
        //-ve Z
        if (minZ() > 0) {
            unsigned int len = (1u << (m_uiMaxDepth - getLevel()));
            unsigned int xres = minX();
            unsigned int yres = minY();
            unsigned int zres = minZ() - len;
            TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
            return res;
        } else {
            TreeNode res(m_uiDim, m_uiMaxDepth);
            return res;
        }
    } //end fn.

    inline TreeNode   TreeNode::getTop() const {
        //+ve Z
        if ((m_uiDim == 3) && (maxZ() < (1u << m_uiMaxDepth))) {
            unsigned int xres = minX();
            unsigned int yres = minY();
            unsigned int zres = maxZ();
            TreeNode res(1, xres, yres, zres, getLevel(), m_uiDim, m_uiMaxDepth);
            return res;
        } else {
            TreeNode res(m_uiDim, m_uiMaxDepth);
            return res;
        }
    } //end fn.

    inline TreeNode   TreeNode::getLeftBack() const {
        return (getLeft().getBack());
    } //end fn.

    inline TreeNode    TreeNode::getRightBack() const {
        return (getRight().getBack());
    } //end fn.

    inline TreeNode   TreeNode::getLeftFront() const {
        return (getLeft().getFront());
    } //end fn.

    inline TreeNode   TreeNode::getRightFront() const {
        return (getRight().getFront());
    } //end fn.

    inline TreeNode   TreeNode::getBottomLeft() const {
        return (getBottom().getLeft());
    } //end fn.

    inline TreeNode   TreeNode::getBottomRight() const {
        return (getBottom().getRight());
    } //end fn.

    inline TreeNode   TreeNode::getBottomBack() const {
        return (getBottom().getBack());
    } //end fn.

    inline TreeNode   TreeNode::getBottomFront() const {
        return (getBottom().getFront());
    } //end fn.

    inline TreeNode   TreeNode::getBottomLeftBack() const {
        return (getBottomLeft().getBack());
    } //end fn.

    inline TreeNode   TreeNode::getBottomRightBack() const {
        return (getBottomRight().getBack());
    } //end fn.

    inline TreeNode   TreeNode::getBottomLeftFront() const {
        return (getBottomLeft().getFront());
    } //end fn.

    inline TreeNode   TreeNode::getBottomRightFront() const {
        return (getBottomRight().getFront());
    } //end fn.

    inline TreeNode   TreeNode::getTopLeft() const {
        return (getTop().getLeft());
    } //end fn.

    inline TreeNode    TreeNode::getTopRight() const {
        return (getTop().getRight());
    } //end fn.

    inline TreeNode    TreeNode::getTopBack() const {
        return (getTop().getBack());
    } //end fn.

    inline TreeNode   TreeNode::getTopFront() const {
        return (getTop().getFront());
    } //end fn.

    inline TreeNode   TreeNode::getTopLeftBack() const {
        return (getTopLeft().getBack());
    } //end fn.

    inline TreeNode    TreeNode::getTopRightBack() const {
        return (getTopRight().getBack());
    } //end fn.

    inline TreeNode   TreeNode::getTopLeftFront() const {
        return (getTopLeft().getFront());
    } //end fn.

    inline TreeNode   TreeNode::getTopRightFront() const {
        return (getTopRight().getFront());
    } //end fn.


    inline int TreeNode::getAnchor(unsigned int &x, unsigned int &y, unsigned int &z) const {
        x = m_uiX;
        y = m_uiY;
        z = m_uiZ;
        return 1;
    }



    inline std::vector<TreeNode> TreeNode::getAllNeighbours() const {
        /*
           0 = Left;  1 =  Right;  2 =  Front;  3 = Back;   4 = LeftBack;  5 = RightBack;  6 = LeftFront;  7 = RightFront;  8 = Top;
           9 = TopRight;  10 =  TopBack;  11 =  TopRightBack;  12 =  Bottom;  13 =  BottomBack;  14 =  TopLeft;  15 =  BottomLeft;
           16 =  BottomRight;  17 =  TopFront;  18 =  BottomFront;  19 =  TopLeftFront;  20 =  TopRightFront;  21 =  BottomLeftFront;
           22 =  BottomRightFront;  23 =  TopLeftBack;  24 = BottomLeftBack;  25 = BottomRightBack;
           */
        std::vector<TreeNode> neighList;

        if (m_uiDim == 3) {
            neighList.resize(26);
            neighList[0] = getLeft();
            neighList[1] =  getRight();
            neighList[2] =  getFront();
            neighList[3] =  getBack();
            neighList[4] =  getLeftBack();
            neighList[5] = getRightBack();
            neighList[6] =  getLeftFront();
            neighList[7] =  getRightFront();
            neighList[8] =  getTop();
            neighList[9] = getTopRight();
            neighList[10] =  getTopBack();
            neighList[11] =  getTopRightBack();
            neighList[12] =  getBottom();
            neighList[13] =  getBottomBack();
            neighList[14] =  getTopLeft();
            neighList[15] =  getBottomLeft();
            neighList[16] =  getBottomRight();
            neighList[17] =  getTopFront();
            neighList[18] =  getBottomFront();
            neighList[19] =  getTopLeftFront();
            neighList[20] =  getTopRightFront();
            neighList[21] =   getBottomLeftFront();
            neighList[22] =   getBottomRightFront();
            neighList[23] =  getTopLeftBack();
            neighList[24] = getBottomLeftBack();
            neighList[25] = getBottomRightBack();
        } else if (m_uiDim == 2) {
            neighList.resize(8);
            neighList[0] = getLeft();
            neighList[1] =  getRight();
            neighList[2] =  getFront();
            neighList[3] =  getBack();
            neighList[4] =  getLeftBack();
            neighList[5] = getRightBack();
            neighList[6] =  getLeftFront();
            neighList[7] =  getRightFront();
        } else {
            neighList.resize(2);
            neighList[0] = getLeft();
            neighList[1] =  getRight();
        }
        return neighList;
    }

    inline unsigned int ot::TreeNode::getMortonIndex() const
    {
        // to get the morton child number.
        unsigned int mid_bit = m_uiMaxDepth - this->getLevel();
        return( (((m_uiZ >> mid_bit) & 1u) << 2u) | (((m_uiY >> mid_bit) & 1u) << 1u) | ((m_uiX >>mid_bit) & 1u));

    }

}


namespace par {

    //Forward Declaration
    template <typename T>
    class Mpi_datatype;

/**
@author Rahul Sampath, rahul.sampath@gmail.com
@brief A template specialization of the abstract class "Mpi_datatype" for communicating messages of type "ot::TreeNode".
*/
    template <>
    class Mpi_datatype< ot::TreeNode > {

        static void Node_MAX_LEVEL(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for(int i = 0; i < (*len); i++) {
                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                (static_cast<ot::TreeNode*>(inout))[i] =
                        ( ( (first.getLevel()) > (second.getLevel()) )? first : second );
            }//end for
        }//end function

        static void Node_MAX(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for(int i = 0; i < (*len); i++) {
                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                (static_cast<ot::TreeNode*>(inout))[i] = ( ( first > second )? first : second );
            }//end for
        }//end function

        static void Node_MIN(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for(int i = 0; i < (*len); i++) {
                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                (static_cast<ot::TreeNode*>(inout))[i] = ( ( first < second )? first : second );
            }//end for
        }//end function

       /* static void Node_NCA(void *in, void *inout, int* len, MPI_Datatype * dptr) {
            for(int i = 0; i < (*len); i++) {
                ot::TreeNode first = (static_cast<ot::TreeNode*>(in))[i];
                ot::TreeNode second = (static_cast<ot::TreeNode*>(inout))[i];
                if( first != second ) {
                    (static_cast<ot::TreeNode*>(inout))[i] = ot::getNCA(first, second);
                }//end if
            }//end for
        }//end function*/

    public:
        /**
          @brief User defined MPI_Operation that sets second[i] to first[i] if first[i] is at a greater level than second[i].
          @remark first and second are 2 arrays of type TreeNode.
        **/
        /*static MPI_Op MAX_LEVEL(){
            static bool first = true;
            static MPI_Op maxLev;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MAX_LEVEL ,true ,&maxLev);
            }
            return maxLev;
        }*/

        /**
          @brief User defined MPI_Operation that computes: second[i] = Max(first[i], second[i]),
          @remark first and second are 2 arrays of type TreeNode.
          "MAX" is a macro
        **/
        static MPI_Op _MAX() {
            static bool         first = true;
            static MPI_Op max;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MAX ,true ,&max);
            }
            return max;
        }

        /**
          @brief User defined MPI_Operation that computes: second[i] = Min(first[i], second[i]),
          @remark first and second are 2 arrays of type TreeNode.
          "MIN" is a macro
        **/
        static MPI_Op _MIN() {
            static bool         first = true;
            static MPI_Op min;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_MIN ,true ,&min);
            }
            return min;
        }

        /**
          @brief User defined MPI_Operation that computes: second[i] = NCA(first[i], second[i]),
          @remark first and second are 2 arrays of type TreeNode and
          NCA returns the nearest common ancestor of its 2 arguments.
        **/
        /*static MPI_Op NCA() {
            static bool         first = true;
            static MPI_Op nca;
            if (first) {
                first = false;
                MPI_Op_create(Mpi_datatype<ot::TreeNode>::Node_NCA ,true ,&nca);
            }
            return nca;
        }*/

        /**
         @return The MPI_Datatype corresponding to the datatype "ot::TreeNode".
       */
        static MPI_Datatype value()
        {
            static bool         first = true;
            static MPI_Datatype datatype;

            if (first)
            {
                first = false;
                MPI_Type_contiguous(sizeof(ot::TreeNode), MPI_BYTE, &datatype);
                MPI_Type_commit(&datatype);
            }

            return datatype;
        }

    };

}//end namespace par







#endif //SFCSORTBENCH_TREENODE_H


