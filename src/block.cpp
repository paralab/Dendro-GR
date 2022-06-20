//
// Created by milinda on 4/21/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief class definition of the class Block.
*/
//


#include "block.h"

ot::Block::Block()
{
    m_uiBlockNode=ot::TreeNode(m_uiDim,m_uiMaxDepth);
    m_uiRegGridLev=0;
    m_uiLocalElementBegin=0;
    m_uiLocalElementEnd=0;

    m_uiPaddingWidth=0;

    m_uiEleOrder=0;
    m_uiSize1D=0;

    m_uiSzX=m_uiSize1D;
    m_uiSzY=m_uiSize1D;
    m_uiSzZ=m_uiSize1D;

    m_uiBlkElem_1D=1u<<(m_uiRegGridLev-m_uiBlockNode.getLevel());

    m_uiIsInternal=false;
    m_uiBlkType = BlockType::UNSPECIFIED;


}


ot::Block::Block(ot::TreeNode pNode, unsigned int rotID ,unsigned int regLev, unsigned int regEleBegin,unsigned int regEleEnd,unsigned int eleOrder)
{

    m_uiBlockNode=pNode;
    m_uiRotID=rotID;
    m_uiRegGridLev=regLev;
    m_uiLocalElementBegin=regEleBegin;
    m_uiLocalElementEnd=regEleEnd;

    m_uiPaddingWidth=(eleOrder>>1u);//GHOST_WIDTH set to 1/2 of the element order (currently unzip mainly tested with even element orders);

    m_uiEleOrder=eleOrder;
    m_uiSize1D=m_uiEleOrder*(1u<<(m_uiRegGridLev-m_uiBlockNode.getLevel()))+1+2*m_uiPaddingWidth;

    m_uiSzX=m_uiSize1D;
    m_uiSzY=m_uiSize1D;
    m_uiSzZ=m_uiSize1D;

    m_uiBlkElem_1D=1u<<(m_uiRegGridLev-m_uiBlockNode.getLevel());

    m_uiIsInternal=false;
    m_uiBlkType = BlockType::UNSPECIFIED;

}


ot::Block::~Block()
{
    m_uiBLK2DIAG.clear();
}


double ot::Block::computeGridDx() const
{
    return ((m_uiBlockNode.maxX()-m_uiBlockNode.minX()))/((double)((1u<<(m_uiRegGridLev-m_uiBlockNode.getLevel()))*m_uiEleOrder));
}

double ot::Block::computeGridDy() const
{
    return ((m_uiBlockNode.maxY()-m_uiBlockNode.minY()))/((double)((1u<<(m_uiRegGridLev-m_uiBlockNode.getLevel()))*m_uiEleOrder));
}

double ot::Block::computeGridDz() const
{
    return ((m_uiBlockNode.maxZ()-m_uiBlockNode.minZ()))/((double)((1u<<(m_uiRegGridLev-m_uiBlockNode.getLevel()))*m_uiEleOrder));
}

double ot::Block::computeDx(const Point & d_min,const Point & d_max) const
{
    return ((d_max.x()-d_min.x())/((double)(1u<<m_uiMaxDepth)))*(computeGridDx());
}

double ot::Block::computeDy(const Point & d_min,const Point & d_max) const
{
    return ((d_max.y()-d_min.y())/((double)(1u<<m_uiMaxDepth)))*(computeGridDy());
}

double ot::Block::computeDz(const Point & d_min,const Point & d_max) const
{
    return ((d_max.z()-d_min.z())/((double)(1u<<m_uiMaxDepth)))*(computeGridDz());
}


void ot::Block::setOffset(DendroIntL offset)
{
    m_uiOffset=offset;
}


void ot::Block::initializeBlkDiagMap(const unsigned int value)
{

    m_uiBLK2DIAG.clear();
    m_uiBLK2DIAG.resize(NUM_EDGES*2*m_uiBlkElem_1D,value);


}

void ot::Block::initializeBlkVertexMap(const unsigned int value)
{
    m_uiBLKVERTX.clear();
    m_uiBLKVERTX.resize(NUM_CHILDREN,value);
}


void ot::Block::computeEleIJK(ot::TreeNode pNode, unsigned int* eijk) const 
{
    eijk[0]=(pNode.getX()-m_uiBlockNode.getX())>>(m_uiMaxDepth-m_uiRegGridLev);
    eijk[1]=(pNode.getY()-m_uiBlockNode.getY())>>(m_uiMaxDepth-m_uiRegGridLev);
    eijk[2]=(pNode.getZ()-m_uiBlockNode.getZ())>>(m_uiMaxDepth-m_uiRegGridLev);

}

bool ot::Block::isBlockInternalEle(ot::TreeNode pNode) const 
{
    unsigned int eijk[3];
    this->computeEleIJK(pNode,eijk);
    // note that integer overflow if the element is out of the min bounds. hence give MAX number which is larger than m_uiBlkElem_1D;
    bool s1 = ( (eijk[0] < (m_uiBlkElem_1D-1) ) && (eijk[1] < (m_uiBlkElem_1D-1) ) && (eijk[2] < (m_uiBlkElem_1D-1) ) );
    bool s2 = ( (eijk[0] >  0 ) && (eijk[1] >  0 ) && (eijk[2] >  0 ) );

    return (s1 && s2);
}