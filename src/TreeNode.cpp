//
// Created by milinda on 2/8/16.
//

#include <sfcSort.h>
#include "TreeNode.h"
unsigned int m_uiMaxDepth=30;

namespace ot
{

    TreeNode::TreeNode (const int dummy, const unsigned int x,const unsigned int y,
              const unsigned int z,const unsigned int level, const unsigned int dim,
              const unsigned int maxDepth)
    {


        //m_uiMaxDepth = maxDepth;
        m_uiX = x;
        if (dim > 1) {m_uiY = y; } else {m_uiY = 0; }
        if (dim > 2) {m_uiZ = z; } else {m_uiZ = 0; }

        m_uiLevel = level;
        //m_uiMaxDepth=maxDepth;

    }



    TreeNode::TreeNode(const unsigned int x, const unsigned int y,
                          const unsigned int z, const unsigned int lev, const unsigned int dim, const unsigned int maxDepth) {

        //m_uiMaxDepth = maxDepth;
        m_uiX = x;
        if (dim > 1) {m_uiY = y; } else {m_uiY = 0; }
        if (dim > 2) {m_uiZ = z; } else {m_uiZ = 0; }

        m_uiLevel = lev;
        //m_uiMaxDepth=maxDepth;

    } //end function

    TreeNode::TreeNode()
    {

      //m_uiMaxDepth=31;
      m_uiX=0; m_uiY=0; m_uiZ=0;
      m_uiLevel=0;

    }

    TreeNode:: TreeNode(const unsigned int dim,const unsigned int maxDepth)
    {
        m_uiX=0;
        m_uiY=0;
        m_uiZ=0;
        m_uiLevel=0;
        //m_uiMaxDepth=maxDepth;
    }

    std::ostream& operator<<(std::ostream& os, TreeNode const& other) {
            return (os << other.getX() << " " << other.getY() << " " << other.getZ() << " " << other.getLevel());
    } //end fn.


    TreeNode TreeNode::getNCA(TreeNode const & other) const {
#ifdef __DEBUG_OCT__
        assert(areComparable(first,other));
    assert(first != other);
#endif
        unsigned int fx = this->getX();
        unsigned int sx = other.getX();
        unsigned int fy = this->getY();
        unsigned int sy = other.getY();
        unsigned int fz = this->getZ();
        unsigned int sz = other.getZ();
        unsigned int maxDepth = this->getMaxDepth();
        unsigned int dim = this->getDim();
        unsigned int maxDiff = (unsigned int)(std::max((std::max((fx^sx),(fy^sy))),(fz^sz)));
        unsigned int maxDiffBinLen = binOp::binLength(maxDiff);
        //Eliminate the last maxDiffBinLen bits.
        unsigned int ncaX = ((fx>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaY = ((fy>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaZ = ((fz>>maxDiffBinLen)<<maxDiffBinLen);
        unsigned int ncaLev = (maxDepth - maxDiffBinLen);
        //assert(ncaLev<std::min(first.getLevel(),other.getLevel()));

//    if(ncaLev>std::min(first.getLevel(),other.getLevel()))
//      ncaLev=std::min(first.getLevel(),other.getLevel());

        TreeNode nca(ncaX,ncaY,ncaZ,ncaLev,dim,maxDepth);
        return nca;

    }//end function




    bool TreeNode::isBoundaryOctant(const TreeNode& block, int type, unsigned char *flags) const {
        unsigned char _flags = 0;

        unsigned int _x = block.getX();
        unsigned int _y = block.getY();
        unsigned int _z = block.getZ();
        unsigned int _d = block.getLevel();

        if ((type & NEGATIVE) == NEGATIVE) {
            // test if any of the anchor values matches those of the block ...
            if (m_uiX == _x) _flags |= X_NEG_BDY;
            if (m_uiY == _y) _flags |= Y_NEG_BDY;
            if (m_uiZ == _z) _flags |= Z_NEG_BDY;
        }

        if ((type & POSITIVE) == POSITIVE) {
            unsigned int len  = (unsigned int)(1u << (m_uiMaxDepth - getLevel()));
            unsigned int blen = ((unsigned int)(1u << (m_uiMaxDepth - _d))) - len;

            if (m_uiX == (_x + blen))  _flags |= X_POS_BDY;
            if (m_uiY == (_y + blen))  _flags |= Y_POS_BDY;
            if (m_uiZ == (_z + blen))  _flags |= Z_POS_BDY;
        }

        if (flags) {
            *flags = _flags;
        }
        if (_flags) {
            return true;
        }
        return false;
    } //end function

    bool TreeNode::isBoundaryOctant(int type, unsigned char *flags) const {

        unsigned char _flags = 0;
        if ((type & NEGATIVE) == NEGATIVE) {
            // test if any of the anchor values is zero ...  (sufficient ??? )
            if (!m_uiX) _flags |= X_NEG_BDY;
            if (!m_uiY) _flags |=  Y_NEG_BDY;
            if (!m_uiZ) _flags |=   Z_NEG_BDY;
        }

        if ((type & POSITIVE) == POSITIVE) {
            unsigned int len  = (unsigned int)(1u << (m_uiMaxDepth - getLevel()));
            unsigned int blen = ((unsigned int)(1u << m_uiMaxDepth)) - len;

            if (m_uiX == blen)  _flags |= X_POS_BDY;
            if (m_uiY == blen)  _flags |= Y_POS_BDY;
            if (m_uiZ == blen)  _flags |= Z_POS_BDY;
        }

        if (flags) *flags = _flags;
        if (_flags) return true;

        return false;
    } //end function

    int TreeNode::addChildren(std::vector<ot::TreeNode>& children) const {
        unsigned int dim = m_uiDim;
        unsigned int maxDepth = m_uiMaxDepth;
        unsigned int childrenSz = children.size();
        children.resize(childrenSz + (1 << dim));

        //#define MORTON_ORDERING

        if ((m_uiLevel & ot::TreeNode::MAX_LEVEL) == maxDepth) {
            for (int i = 0; i < (1 << dim); i++) {
                children[childrenSz + i] = *this;
            }
            return 1;
        }
        //The check that lev < maxD is taken care of in the constructor.

        //Order: X first, Y next and Z last

        unsigned int len = (unsigned int)(1u << (maxDepth - ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1)));

        TreeNode   tmpNode0(1, m_uiX, m_uiY, m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
        children[childrenSz + 0] = tmpNode0;

        TreeNode   tmpNode1(1, (m_uiX + len), m_uiY, m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
        children[childrenSz + 1] = tmpNode1;

        if (dim >= 2) {
            TreeNode   tmpNode2(1, m_uiX, (m_uiY + len), m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[childrenSz + 2] = tmpNode2;

            TreeNode   tmpNode3(1, (m_uiX + len), (m_uiY + len), m_uiZ, ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[childrenSz + 3] = tmpNode3;
        }

        if (dim == 3) {
            TreeNode   tmpNode4(1, m_uiX, m_uiY, (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[childrenSz + 4] = tmpNode4;

            TreeNode   tmpNode5(1, (m_uiX + len), m_uiY, (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[childrenSz + 5] = tmpNode5;

            TreeNode   tmpNode6(1, m_uiX, (m_uiY + len), (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[childrenSz + 6] = tmpNode6;

            TreeNode   tmpNode7(1, (m_uiX + len), (m_uiY + len), (m_uiZ + len), ((m_uiLevel & ot::TreeNode::MAX_LEVEL) + 1), m_uiDim, m_uiMaxDepth);
            children[childrenSz + 7] = tmpNode7;
        } //end if

#ifdef HILBERT_ORDERING
#pragma message("===FIX ME===")
        std::sort(children.begin(), children.end());
        /*std::vector<ot::TreeNode> tmp;
        ot::TreeNode root=ot::TreeNode(0,0,0,0,m_uiDim,m_uiMaxDepth);
        SFC::seqSort::SFC_treeSort(&(*(children.begin())),children.size(),tmp,tmp,tmp,m_uiMaxDepth,m_uiMaxDepth,root,0,1,0);*/
#endif
        return 1;
    } //end function


    TreeNode TreeNode::getDFDMorton() const
    {
        TreeNode dfd(m_uiX, m_uiY, m_uiZ, m_uiMaxDepth, m_uiDim, m_uiMaxDepth);
        return dfd;
    }



    TreeNode TreeNode::getDFD() const
    {
        std::vector<ot::TreeNode> keys;
        keys.resize((1u<<m_uiDim));

        keys[0]=(ot::TreeNode(this->minX(),this->minY(),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
        keys[1]=(ot::TreeNode((this->maxX()-1),this->minY(),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
        keys[2]=(ot::TreeNode(this->minX(),(this->maxY()-1),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
        keys[3]=(ot::TreeNode((this->maxX()-1),(this->maxY()-1),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));

        if(m_uiDim==3){

            keys[4]=(ot::TreeNode(this->minX(),this->minY(),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
            keys[5]=(ot::TreeNode((this->maxX()-1),this->minY(),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
            keys[6]=(ot::TreeNode(this->minX(),(this->maxY()-1),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
            keys[7]=(ot::TreeNode((this->maxX()-1),(this->maxY()-1),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));

        }



        std::vector<ot::TreeNode> tmpNodes;
        ot::TreeNode rootNode(0,0,0,0,m_uiDim,m_uiMaxDepth);
        SFC::seqSort::SFC_treeSort(&(*(keys.begin())),keys.size(),tmpNodes,tmpNodes,tmpNodes,m_uiMaxDepth,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_SORT_ONLY);
        assert(seq::test::isUniqueAndSorted(keys));
        return keys.front();

    }

    TreeNode TreeNode::getDLD() const
    {
        std::vector<ot::TreeNode> keys;
        keys.resize((1u<<m_uiDim));

        keys[0]=(ot::TreeNode(this->minX(),this->minY(),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
        keys[1]=(ot::TreeNode((this->maxX()-1),this->minY(),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
        keys[2]=(ot::TreeNode(this->minX(),(this->maxY()-1),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
        keys[3]=(ot::TreeNode((this->maxX()-1),(this->maxY()-1),this->minZ(),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));

        if(m_uiDim==3){

            keys[4]=(ot::TreeNode(this->minX(),this->minY(),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
            keys[5]=(ot::TreeNode((this->maxX()-1),this->minY(),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
            keys[6]=(ot::TreeNode(this->minX(),(this->maxY()-1),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));
            keys[7]=(ot::TreeNode((this->maxX()-1),(this->maxY()-1),(this->maxZ()-1),m_uiMaxDepth,m_uiDim,m_uiMaxDepth));

        }

        std::vector<ot::TreeNode> tmpNodes;
        ot::TreeNode rootNode(0,0,0,0,m_uiDim,m_uiMaxDepth);
        SFC::seqSort::SFC_treeSort(&(*(keys.begin())),keys.size(),tmpNodes,tmpNodes,tmpNodes,m_uiMaxDepth,m_uiMaxDepth,rootNode,ROOT_ROTATION,1,TS_SORT_ONLY);
        assert(seq::test::isUniqueAndSorted(keys));
        return keys.back();
    }





}