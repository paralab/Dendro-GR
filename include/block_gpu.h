/**
 * @file block_gpu.h
 * @author Milinda Fernando
 * @brief Device independent GPU block class for FD/FV computations on adaptive grids. 
 * @version 0.1
 * @date 2021-11-29
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once
#include "dendro.h"
#include "TreeNode.h"
#include "block.h"

/**
 * @brief simple structure to store the TreeNodes. 
 * 
 */
struct TNode3d
{
    public:
        DEVICE_UINT m_x=0;
        DEVICE_UINT m_y=0;
        DEVICE_UINT m_z=0;
        DEVICE_UINT m_level=0;
    
    #ifdef __CUDACC__
        __host__ __device__ 
    #endif
    TNode3d(){}
    TNode3d(const ot::TreeNode& tnode)
    {
        m_x = tnode.minX();
        m_y = tnode.minY();
        m_z = tnode.minZ();
        m_level = tnode.getLevel();
    }
    
};

struct BlockGPU3D
{
        /**@brief: min coordinate of the block*/
        DEVICE_REAL m_ptMin[3]={0,0,0};

        /**@brief: dx parameters dx,dy,dz*/
        DEVICE_REAL m_dx[3]={0,0,0};

        /**@brief: sz parameters */
        DEVICE_UINT m_sz[3]={0,0,0};

        /**@brief: sz parameters with aligned blocks*/
        DEVICE_UINT m_aligned_sz[3]={0,0,0};

        /**@brief boundry flag*/
        DEVICE_UINT m_bflag=0;

        /**@brief offset value for series of blocks*/
        DEVICE_UINT m_offset=0;

        /**@brief capturing tree node corresponding to the block. */
        TNode3d m_blockNode;

        /**@brief number of regular splits, from the block node*/
        DEVICE_UINT m_patch_level =0;

        /**@brief elements bounds captured in the block. */
        DEVICE_UINT m_elem_begin = 0;
        DEVICE_UINT m_elem_end   = 0;
        
        BlockGPU3D(){}
        BlockGPU3D(const ot::Block& block)
        {
                m_offset = block.getOffset();
                
                m_sz[0]  = block.getAllocationSzX();
                m_sz[1]  = block.getAllocationSzY();
                m_sz[2]  = block.getAllocationSzZ();

                m_aligned_sz[0]  = block.getAllocationSzX();
                m_aligned_sz[1]  = block.getAllocationSzY();
                m_aligned_sz[2]  = block.getAllocationSzZ();

                m_bflag  = block.getBlkNodeFlag();

                m_elem_begin = block.getLocalElementBegin();
                m_elem_end   = block.getLocalElementEnd();
                m_patch_level = block.getRegularGridLev();

                m_blockNode = TNode3d(block.getBlockNode());
                
        }

        
};