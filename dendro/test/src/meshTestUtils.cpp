/**
 * @file meshTestUtils.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief meshTest utility functiuons. 
 * @version 0.1
 * @date 2020-01-16
 * 
 * School of Computing, University of Utah.
 * @copyright Copyright (c) 2020
 * 
 */

#include "meshTestUtils.h"


bool ot::test::isBlkFlagsValid(const ot::Mesh* pMesh)
{
    if(pMesh->isActive())
    {
        const ot::TreeNode * pNodes= &(*(pMesh->getAllElements().begin()));

        const unsigned int nodeLocalBegin = pMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd = pMesh->getNodeLocalEnd();

        const unsigned int* e2n_cg = &(*(pMesh->getE2NMapping().begin()));
        const unsigned int* e2e = &(*(pMesh->getE2EMapping().begin()));

        const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
        const unsigned int nPe = pMesh->getNumNodesPerElement();

        unsigned int lookup,node_cg;
        unsigned int child[NUM_CHILDREN];

        for(unsigned int blk =0; blk < blkList.size(); blk ++)
        {
            std::vector<unsigned int > gid;
            computeBlockUnzipGhostNodes(pMesh,blk,gid);

            if( (gid.empty() && blkList[blk].getBlockType() == ot::BlockType::UNZIP_DEPENDENT ) ||  (!gid.empty() && blkList[blk].getBlockType() == ot::BlockType::UNZIP_INDEPENDENT ))
                return false;
            
        }

        return true;
     
    }else
        return true;

    
}
