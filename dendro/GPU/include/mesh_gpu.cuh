/**
 * @file mesh_gpu.cuh
 * @author Milinda Fernando (milinda@oden.utexas.edu)
 * @brief minimal mesh class, for CUDA only, specialization of the mesh_gpu.h class. 
 * @version 0.1
 * @date 2022-01-21
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include "mesh_gpu.h"
#include "device.h"
#include "refel_const.cuh"
#include "asyncExchangeContex.h"

namespace device
{
    #define IDX_DG(i,j,k) ((i) + nx * ( (j) + nx * (k)))

    DEVICE_FUNC
    inline void dg2eijk(DEVICE_UINT p, DEVICE_UINT dg_index, DEVICE_UINT &e, DEVICE_UINT &i, DEVICE_UINT &j, DEVICE_UINT &k)
    {
        const DEVICE_UINT nPe = (p+1)*(p+1)*(p+1);
        e = dg_index / nPe;
        k = 0;
        j = 0;
        i = 0;

        //std::cout<<"e: "<<e<<std::endl;
        if (dg_index > e * nPe)
            k = (dg_index - e * nPe) / ((p + 1) * (p + 1));
        //std::cout<<"k: "<<k<<std::endl;
        if ((dg_index + k * ((p + 1) * (p + 1))) > (e * nPe))
            j = (dg_index - e * nPe - k * ((p + 1) * (p + 1))) / (p + 1);
        //std::cout<<"j: "<<j<<std::endl;
        if ((dg_index + k * ((p + 1) * (p + 1)) + j * (p + 1)) > (e * nPe))
            i = (dg_index - e * nPe - k * ((p + 1) * (p + 1)) - j * (p + 1));
        //std::cout<<"i: "<<i<<std::endl;

        return;
    }

    template<DEVICE_UINT p>
    DEVICE_FUNC 
    inline DEVICE_BOOL is_node_hanging(const MeshGPU* const dptr_mesh, DEVICE_UINT eleID, DEVICE_UINT ix, DEVICE_UINT jy, DEVICE_UINT kz)
    {
        const DEVICE_UINT   nPe      = (p+1) * (p+1) * (p+1);
        const DEVICE_UINT   nx       = p+1;
        const DEVICE_UINT   owner_id = (dptr_mesh->m_e2n_dg[eleID * (nPe) + IDX_DG(ix,jy,kz)]/nPe); 
        return dptr_mesh->m_all_elements[owner_id].m_level <  dptr_mesh->m_all_elements[eleID].m_level;
    }

    template<DEVICE_UINT p>
    DEVICE_FUNC 
    DEVICE_BOOL is_edge_hanging(DEVICE_UINT sx, DEVICE_UINT sy, DEVICE_UINT sz , const MeshGPU* const dptr_mesh, DEVICE_UINT elementId, DEVICE_UINT& cnum)
    {
        assert(p>1);
        const DEVICE_UINT nPe       = (p+1) * (p+1) * (p+1);
        const DEVICE_UINT m_uiNpE   = nPe;
        const DEVICE_UINT nx        = p+1;

        DEVICE_BOOL isHanging;
        DEVICE_BOOL isVertexHanging[2];
        DEVICE_UINT nodeLookUp_DG;
        
        DEVICE_UINT owner[2];
        DEVICE_UINT ii_x[2];
        DEVICE_UINT jj_y[2];
        DEVICE_UINT kk_z[2];
        DEVICE_UINT edge_id;
       
        if(sx==0 && sy==0 && sz==1)
            edge_id=OCT_DIR_LEFT_DOWN;
        else if(sx==0 && sy==2 && sz==1)
            edge_id=OCT_DIR_LEFT_UP;
        else if(sx==0 && sy==1 && sz==0)
            edge_id=OCT_DIR_LEFT_BACK;
        else if(sx==0 && sy==1 && sz==2)
            edge_id=OCT_DIR_LEFT_FRONT;
        else if(sx==2 && sy==0 && sz==1)
            edge_id=OCT_DIR_RIGHT_DOWN;
        else if(sx==2 && sy==2 && sz==1)
            edge_id=OCT_DIR_RIGHT_UP;
        else if(sx==2 && sy==1 && sz==0)
            edge_id=OCT_DIR_RIGHT_BACK;
        else if(sx==2 && sy==1 && sz==2)
            edge_id=OCT_DIR_RIGHT_FRONT;
        else if(sx==1 && sy==0 && sz==0)
            edge_id=OCT_DIR_DOWN_BACK;
        else if(sx==1 && sy==0 && sz==2)
            edge_id=OCT_DIR_DOWN_FRONT;
        else if(sx==1 && sy==2 && sz==0)
            edge_id=OCT_DIR_UP_BACK;
        else if(sx==1 && sy==2 && sz==2)
            edge_id=OCT_DIR_UP_FRONT;
        else
        {
            printf("s = (%d, %d, %d)\n ", sx, sy, sz);
            assert(false);
        }
            

        switch (edge_id)
        {

            case OCT_DIR_LEFT_DOWN:
                    isVertexHanging[0] =  is_node_hanging<p>(dptr_mesh,elementId,0,0,0);
                    isVertexHanging[1] =  is_node_hanging<p>(dptr_mesh,elementId,0,0,p);

                    isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                    if(!isHanging)
                        return false;

                    nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+((p>>1u)) *(p+1)*(p+1) + (0)*(p+1) +(0)];
                    
                    if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_z==dptr_mesh->m_all_elements[elementId].m_z)
                        cnum=0;
                    else
                        cnum=1;
                
                
                break;

            case OCT_DIR_LEFT_UP:
                    isVertexHanging[0] =  is_node_hanging<p> (dptr_mesh,elementId,0,p,0);
                    isVertexHanging[1] =  is_node_hanging<p> (dptr_mesh,elementId,0,p,p);

                    isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                    if(!isHanging)
                        return false;

                    nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+((p>>1u)) *(p+1)*(p+1) + (p) *(p+1) + (0)];
                    if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_z==dptr_mesh->m_all_elements[elementId].m_z)
                        cnum=0;
                    else
                        cnum=1;
                
                break;

            case OCT_DIR_LEFT_BACK:
                    isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,0,0,0);
                    isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,0,p,0);

                    isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                    if(!isHanging)
                        return false;
                    

                    nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(0)*(p+1)*(p+1)+((p>>1u))*(p+1)+(0)];
                    if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_y==dptr_mesh->m_all_elements[elementId].m_y)
                        cnum=0;
                    else
                        cnum=1;
                break;

            case OCT_DIR_LEFT_FRONT:
                    isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,0,0,p);
                    isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,0,p,p);

                    isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                    if(!isHanging)
                        return false;

                    nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(p)*(p+1)*(p+1)+((p>>1u))*(p+1)+(0)];
                    if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_y==dptr_mesh->m_all_elements[elementId].m_y)
                        cnum=0;
                    else
                        cnum=1;

                break;

            case OCT_DIR_RIGHT_DOWN:
                    isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,p,0,0);
                    isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,0,p);

                    isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                    if(!isHanging)
                        return false;

                    
                    nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+((p>>1u))*(p+1)*(p+1)+(0)*(p+1)+(p)];
                    if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_z==dptr_mesh->m_all_elements[elementId].m_z)
                        cnum=0;
                    else
                        cnum=1;

                break;

            case OCT_DIR_RIGHT_UP:
                    isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,p,p,0);
                    isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,p,p);

                    isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                    if(!isHanging)
                        return false;

                    nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+((p>>1u))*(p+1)*(p+1)+(p)*(p+1)+(p)];
                    if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_z==dptr_mesh->m_all_elements[elementId].m_z)
                        cnum=0;
                    else
                        cnum=1;

                
                break;

            case OCT_DIR_RIGHT_BACK:
                    isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,p,0,0);
                    isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,p,0);

                    isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                    if(!isHanging)
                        return false;

                    
                        nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(0)*(p+1)*(p+1)+((p>>1u))*(p+1)+(p)];
                        if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_y==dptr_mesh->m_all_elements[elementId].m_y)
                            cnum=0;
                        else
                            cnum=1;
                        
                break;

            case OCT_DIR_RIGHT_FRONT:
                
                isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,p,0,p);
                isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,p,p);

                isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                if(!isHanging)
                    return false;

                nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(p)*(p+1)*(p+1)+((p>>1u))*(p+1)+(p)];
                if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_y==dptr_mesh->m_all_elements[elementId].m_y)
                    cnum=0;
                else
                    cnum=1;
                
                break;

            case OCT_DIR_DOWN_BACK:
                
                isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,0,0,0);
                isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,0,0);

                isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                if(!isHanging)
                    return false;

                nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(0)*(p+1)*(p+1)+(0)*(p+1)+(p>>1u)];
                if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_x==dptr_mesh->m_all_elements[elementId].m_x)
                    cnum=0;
                else
                    cnum=1;
                
                break;

            case OCT_DIR_DOWN_FRONT:
                
                isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,0,0,p);
                isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,0,p);

                isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                if(!isHanging)
                    return false;
                
                nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(p)*(p+1)*(p+1)+(0)*(p+1)+(p>>1u)];
                if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_x==dptr_mesh->m_all_elements[elementId].m_x)
                    cnum=0;
                else 
                    cnum=1;
                
                
                break;


            case OCT_DIR_UP_BACK:
               
                isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,0,p,0);
                isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,p,0);

                isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                if(!isHanging)
                    return false;
                
                nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(0)*(p+1)*(p+1)+(p)*(p+1)+(p>>1u)];
                if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_x==dptr_mesh->m_all_elements[elementId].m_x)
                    cnum=0;
                else
                    cnum=1;
                
                
                break;

            case OCT_DIR_UP_FRONT:
                
                isVertexHanging[0] =  is_node_hanging<p> ( dptr_mesh,elementId,0,p,p);
                isVertexHanging[1] =  is_node_hanging<p> ( dptr_mesh,elementId,p,p,p);

                isHanging = (isVertexHanging[0] && isVertexHanging[1]);
                if(!isHanging)
                    return false;

                
                nodeLookUp_DG=dptr_mesh->m_e2n_dg[elementId*m_uiNpE+(p)*(p+1)*(p+1)+(p)*(p+1)+(p>>1u)];
                if(dptr_mesh->m_all_elements[nodeLookUp_DG/m_uiNpE].m_x==dptr_mesh->m_all_elements[elementId].m_x)
                    cnum=0;
                else
                    cnum=1;
                
                
                break;


        }
        
        return isHanging;
        
    }
    
    template<DEVICE_UINT p>
    DEVICE_FUNC 
    DEVICE_BOOL is_face_hanging(DEVICE_UINT sx, DEVICE_UINT sy, DEVICE_UINT sz ,const MeshGPU* const dptr_mesh, DEVICE_UINT elementId, DEVICE_UINT& cnum)
    {
        
        DEVICE_BOOL isHanging = false;
        const DEVICE_UINT m_uiMaxDepth = dptr_mesh->m_max_depth;
        DEVICE_UINT ownerID;
        DEVICE_UINT mid_bit;
        cnum=0;

        DEVICE_UINT face_id;
        if(sx==0 && sy==1 && sz==1)
            face_id=OCT_DIR_LEFT;
        else if(sx==2 && sy==1 && sz==1)
            face_id=OCT_DIR_RIGHT;
        else if(sx==1 && sy==0 && sz==1)
            face_id=OCT_DIR_DOWN;
        else if(sx==1 && sy==2 && sz==1)
            face_id=OCT_DIR_UP;
        else if(sx==1 && sy==1 && sz==0)
            face_id=OCT_DIR_BACK;
        else if(sx==1 && sy==1 && sz==2)
            face_id=OCT_DIR_FRONT;
        else
            assert(false);
        

        const DEVICE_UINT lookup=dptr_mesh->m_e2e[elementId * NUM_FACES + face_id];
        if(lookup==LOOK_UP_TABLE_DEFAULT)
            isHanging=false;
        else
            isHanging = dptr_mesh->m_all_elements[lookup].m_level<dptr_mesh->m_all_elements[elementId].m_level;
        
        // Note: if the face is Hanging it is reliable to use the element to element map. 
        if(isHanging)
        {
            switch (face_id)
            {

                case OCT_DIR_LEFT:
                    ownerID=lookup;
                    mid_bit=m_uiMaxDepth - dptr_mesh->m_all_elements[ownerID].m_level-1;
                    cnum=( (((((dptr_mesh->m_all_elements[elementId].m_z)-(dptr_mesh->m_all_elements[ownerID].m_z)) >> mid_bit) & 1u) << 1u) | ((((dptr_mesh->m_all_elements[elementId].m_y)-(dptr_mesh->m_all_elements[ownerID].m_y)) >>mid_bit) & 1u));

                    break;

                case OCT_DIR_RIGHT:
                    ownerID=lookup;
                    mid_bit=m_uiMaxDepth - dptr_mesh->m_all_elements[ownerID].m_level-1;
                    cnum=( (((((dptr_mesh->m_all_elements[elementId].m_z)-(dptr_mesh->m_all_elements[ownerID].m_z)) >> mid_bit) & 1u) << 1u) | ((((dptr_mesh->m_all_elements[elementId].m_y)-(dptr_mesh->m_all_elements[ownerID].m_y)) >>mid_bit) & 1u));

                    break;

                case OCT_DIR_DOWN:
                    ownerID=lookup;
                    mid_bit=m_uiMaxDepth - dptr_mesh->m_all_elements[ownerID].m_level-1;
                    cnum=( (((((dptr_mesh->m_all_elements[elementId].m_z)-(dptr_mesh->m_all_elements[ownerID].m_z)) >> mid_bit) & 1u) << 1u) | ((((dptr_mesh->m_all_elements[elementId].m_x)-(dptr_mesh->m_all_elements[ownerID].m_x)) >>mid_bit) & 1u));
                    break;

                case OCT_DIR_UP:
                    ownerID=lookup;
                    mid_bit=m_uiMaxDepth - dptr_mesh->m_all_elements[ownerID].m_level-1;
                    cnum=( (((((dptr_mesh->m_all_elements[elementId].m_z)-(dptr_mesh->m_all_elements[ownerID].m_z)) >> mid_bit) & 1u) << 1u) | ((((dptr_mesh->m_all_elements[elementId].m_x)-(dptr_mesh->m_all_elements[ownerID].m_x)) >>mid_bit) & 1u));
                    break;

                case OCT_DIR_BACK:
                    ownerID=lookup;
                    mid_bit=m_uiMaxDepth - dptr_mesh->m_all_elements[ownerID].m_level-1;
                    cnum=( (((((dptr_mesh->m_all_elements[elementId].m_y)-(dptr_mesh->m_all_elements[ownerID].m_y)) >> mid_bit) & 1u) << 1u) | ((((dptr_mesh->m_all_elements[elementId].m_x)-(dptr_mesh->m_all_elements[ownerID].m_x)) >>mid_bit) & 1u));
                    break;

                case OCT_DIR_FRONT:
                    ownerID=lookup;
                    mid_bit=m_uiMaxDepth - dptr_mesh->m_all_elements[ownerID].m_level-1;
                    cnum=( (((((dptr_mesh->m_all_elements[elementId].m_y)-(dptr_mesh->m_all_elements[ownerID].m_y)) >> mid_bit) & 1u) << 1u) | ((((dptr_mesh->m_all_elements[elementId].m_x)-(dptr_mesh->m_all_elements[ownerID].m_x)) >>mid_bit) & 1u));
                    break;
            }

        }

        return isHanging;

    }


    DEVICE_FUNC 
    inline DEVICE_BOOL is_hinfo_edge_hanging(DEVICE_UINT sx, DEVICE_UINT sy, DEVICE_UINT sz , const HangingInfo& hinfo_ele, DEVICE_UINT& cnum)
    {
        DEVICE_UINT edge_id;

        if(sx==0 && sy==0 && sz==1)
            edge_id=OCT_DIR_LEFT_DOWN;
        else if(sx==0 && sy==2 && sz==1)
            edge_id=OCT_DIR_LEFT_UP;
        else if(sx==0 && sy==1 && sz==0)
            edge_id=OCT_DIR_LEFT_BACK;
        else if(sx==0 && sy==1 && sz==2)
            edge_id=OCT_DIR_LEFT_FRONT;
        else if(sx==2 && sy==0 && sz==1)
            edge_id=OCT_DIR_RIGHT_DOWN;
        else if(sx==2 && sy==2 && sz==1)
            edge_id=OCT_DIR_RIGHT_UP;
        else if(sx==2 && sy==1 && sz==0)
            edge_id=OCT_DIR_RIGHT_BACK;
        else if(sx==2 && sy==1 && sz==2)
            edge_id=OCT_DIR_RIGHT_FRONT;
        else if(sx==1 && sy==0 && sz==0)
            edge_id=OCT_DIR_DOWN_BACK;
        else if(sx==1 && sy==0 && sz==2)
            edge_id=OCT_DIR_DOWN_FRONT;
        else if(sx==1 && sy==2 && sz==0)
            edge_id=OCT_DIR_UP_BACK;
        else if(sx==1 && sy==2 && sz==2)
            edge_id=OCT_DIR_UP_FRONT;
        else
        {
            printf("s = (%d, %d, %d)\n ", sx, sy, sz);
            assert(false);
        }

        cnum = hinfo_ele.get_ecnum(edge_id);
        return hinfo_ele.is_hanging(edge_id);

    }

    DEVICE_FUNC 
    inline DEVICE_BOOL is_hinfo_face_hanging(DEVICE_UINT sx, DEVICE_UINT sy, DEVICE_UINT sz , const HangingInfo& hinfo_ele, DEVICE_UINT& cnum)
    {
        DEVICE_UINT face_id;
        if(sx==0 && sy==1 && sz==1)
            face_id=OCT_DIR_LEFT;
        else if(sx==2 && sy==1 && sz==1)
            face_id=OCT_DIR_RIGHT;
        else if(sx==1 && sy==0 && sz==1)
            face_id=OCT_DIR_DOWN;
        else if(sx==1 && sy==2 && sz==1)
            face_id=OCT_DIR_UP;
        else if(sx==1 && sy==1 && sz==0)
            face_id=OCT_DIR_BACK;
        else if(sx==1 && sy==1 && sz==2)
            face_id=OCT_DIR_FRONT;
        else
            assert(false);

        cnum = hinfo_ele.get_fcnum(face_id);
        return hinfo_ele.is_hanging(face_id);
    }


    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(512)
    void __unzip_dg__3(const MeshGPU* const __restrict__ dptr_mesh, const T* const __restrict__ dptr_in, T* const __restrict__ dptr_out)
    {
        const DEVICE_UINT ele    = GPUDevice::block_id_x();
        const DEVICE_UINT blk_id = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();

        // implies no contribution for the unzip. 
        if(blk_id < dptr_mesh->m_e2b_unzip_map_count[ele])
        {
            const DEVICE_UINT tx     = GPUDevice::thread_id_x();
            const DEVICE_UINT ty     = GPUDevice::thread_id_y();
            const DEVICE_UINT tz     = GPUDevice::thread_id_z();
            
            const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
            const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
            const DEVICE_UINT dg_sz    = dptr_mesh->m_oct_local_sz;
            const DEVICE_UINT nx       = p+1;
            const DEVICE_UINT nPe      = nx * nx * nx; 

            const DEVICE_UINT blk = dptr_mesh->m_e2b_unzip_map[dptr_mesh->m_e2b_unzip_map_offset[ele] + blk_id];

            const DEVICE_UINT lx     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
            const DEVICE_UINT ly     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
            const DEVICE_UINT lz     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
            const DEVICE_UINT offset =  dptr_mesh->m_blk_list[blk].m_offset;
            
            const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
            const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;
            const DEVICE_UINT bLev         = dptr_mesh->m_all_elements[dptr_mesh->m_blk_list[blk].m_elem_begin].m_level; // or blkNode.m_level + (patch_level -blkNode.m_level);

            const DEVICE_UINT max_depth = dptr_mesh->m_max_depth;
            const TNode3d elem_node = dptr_mesh->m_all_elements[ele]; 
            if(ele>= dptr_mesh->m_blk_list[blk].m_elem_begin && ele < dptr_mesh->m_blk_list[blk].m_elem_end)
            {
                if(!(tx < nx && ty < nx && tz <nx))
                    return;
                
                const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
                const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
                const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);
                
                dptr_out[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)] = dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];
                
                return;
                
            }
            
            const DEVICE_FLOAT32  hx  = (1u<<( max_depth - bLev))/(DEVICE_FLOAT32)p;
            const DEVICE_FLOAT32 xmin = blkNode.m_x - PW*hx; 
            const DEVICE_FLOAT32 ymin = blkNode.m_y - PW*hx; 
            const DEVICE_FLOAT32 zmin = blkNode.m_z - PW*hx; 

            const DEVICE_FLOAT32 xmax = blkNode.m_x + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 ymax = blkNode.m_y + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 zmax = blkNode.m_z + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            
            
            //const DEVICE_FLOAT32 d_compar_tol=1e-6;

            SHARED_MEM T dgVec[p+1][p+1][p+1];
            SHARED_MEM T refEL[2][p+1][p+1];
            SHARED_MEM T v1[p+1][p+1][p+1];
            SHARED_MEM T v2[p+1][p+1][p+1];

            if(elem_node.m_level == bLev)
            {

                if(!(tx < nx && ty < nx && tz <nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/hh;

                const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                if((kkz>=0 && kkz<lz) && (jjy>=0 && jjy<ly) && (iix>=0 && iix<lx))
                    dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];

                return;

            }else if(elem_node.m_level == (bLev + 1))
            {

                if(!(tx < nx && ty < nx && tz <nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/(2*hh);
                const DEVICE_UINT cb =0;//(p%2==0) ? 0 : 1;

                if(tx %2 == cb && ty%2==cb && tz%2==cb)
                {
                    const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                    const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                    const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                    const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                    const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                    const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                    if((kkz>=0 && kkz<lz) && (jjy>=0 && jjy<ly) && (iix>=0 && iix<lx))
                        dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];

                }

                return;

            }else
            {
                assert(bLev == elem_node.m_level+1);
                TNode3d child_node;
                const DEVICE_UINT child_sz = (1u<<(max_depth-elem_node.m_level-1));

                if(tx < nx && ty < nx && tz <nx)
                {
                    dgVec[tz][ty][tx]  = dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];
                    refEL[0][ty][tx]= refel_1d[ 0 * (p+1)*(p+1) + ty * (p+1) + tx];
                    refEL[1][ty][tx]= refel_1d[ 1 * (p+1)*(p+1) + ty * (p+1) + tx];
                }
                
                GPUDevice::sync_threads();

                for(DEVICE_UINT child = 0; child < NUM_CHILDREN; child++)
                {
                    const DEVICE_UINT cbit[3] = {(child & 1u), ((child & 2u)>>1u), ((child & 4u)>>2u)};
                    child_node.m_x = elem_node.m_x + cbit[0] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_y = elem_node.m_y + cbit[1] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_z = elem_node.m_z + cbit[2] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_level = elem_node.m_level+1;

                    if((child_node.m_x + child_sz < xmin || child_node.m_x >=xmax) || 
                    (child_node.m_y + child_sz < ymin || child_node.m_y >=ymax) ||
                    (child_node.m_z + child_sz < zmin || child_node.m_z >=zmax) )
                        continue;
                    
                    const DEVICE_FLOAT32 hh = (1u<<(max_depth - child_node.m_level))/(DEVICE_FLOAT32) p;
                    const DEVICE_FLOAT32 invhh = 1.0/hh;

                    const DEVICE_FLOAT32 zz  = child_node.m_z + tz*hh;
                    const DEVICE_FLOAT32 yy  = child_node.m_y + ty*hh;
                    const DEVICE_FLOAT32 xx  = child_node.m_x + tx*hh;

                    const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                    const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                    const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                    if((tx < nx && ty < nx && tz<nx))
                    {
                        DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[0]][tx][m] * dgVec[tz][ty][m];

                        v1[tz][ty][tx]=val;
                    }
                        
                    GPUDevice::sync_threads();

                    if((tx < nx && ty < nx && tz<nx))
                    {
                        DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[1]][ty][m] * v1[tz][m][tx];
                        v2[tz][ty][tx]=val;
                    }
                        

                    GPUDevice::sync_threads();

                    if((tx < nx && ty < nx && tz<nx))
                    {  
                        DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[2]][tz][m] * v2[m][ty][tx];

                        v1[tz][ty][tx]=val;

                    }
                    
                    GPUDevice::sync_threads();

                    if(!(tx < nx && ty < nx && tz <nx))
                        return;

                    if((kkz>=0 && kkz<lz) && (jjy>=0 && jjy<ly) && (iix>=0 && iix<lx))
                    {
                        // !! the full outer product complexity wise expensive. 
                        // DEVICE_REAL interp_out = 0;
                        // for(DEVICE_UINT mm_z=0; mm_z< (p+1); mm_z++)
                        //     for(DEVICE_UINT mm_y=0; mm_y< (p+1); mm_y++)
                        //         for(DEVICE_UINT mm_x=0; mm_x< (p+1); mm_x++)
                        //             interp_out += refEL[cbit[0]][tx][mm_x] * refEL[cbit[1]][ty][mm_y] * refEL[cbit[2]][tz][mm_z] * dgVec[mm_z][mm_y][mm_x];
                        //printf("interp out : %.8E  val = %.8E\n", interp_out, v1[tz][ty][tx]);
                        dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  v1[tz][ty][tx];
                    }
                        

                }

                return;

            }


        }
        
        return;

    }
    

    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __block_internal_unzip_dg__2(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        // copy only the block internel element,
        const DEVICE_UINT blk    = GPUDevice::block_id_x();
        const DEVICE_UINT ele_id = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();

        const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        const DEVICE_UINT ty     = GPUDevice::thread_id_y();

        const DEVICE_UINT nx       = p+1;
        const DEVICE_UINT nPe      = nx * nx * nx;

        if(!(tx < nx && ty < nx))
            return;

        const DEVICE_UINT ele      = dptr_mesh->m_blk_list[blk].m_elem_begin + ele_id;
        const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
        const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
        const DEVICE_UINT dg_sz    = dptr_mesh->m_oct_local_sz;

        const DEVICE_UINT lx       =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
        const DEVICE_UINT ly       =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
        const DEVICE_UINT lz       =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
        const DEVICE_UINT offset   =  dptr_mesh->m_blk_list[blk].m_offset;
        
        const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
        const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;
        const DEVICE_UINT max_depth    = dptr_mesh->m_max_depth;

        assert(ele>= dptr_mesh->m_blk_list[blk].m_elem_begin && ele < dptr_mesh->m_blk_list[blk].m_elem_end);

        const TNode3d elem_node = dptr_mesh->m_all_elements[ele];
        const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
        const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
        const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);
        
        for(DEVICE_UINT tz=0; tz<nx; tz++)
            dptr_out[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)] = dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];
        
        return;

    }

    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __unzip_dg__2(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        const DEVICE_UINT ele    = GPUDevice::block_id_x();
        const DEVICE_UINT blk_id = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();

        if(blk_id < dptr_mesh->m_e2b_unzip_map_count[ele])
        {

            const DEVICE_UINT tx     = GPUDevice::thread_id_x();
            const DEVICE_UINT ty     = GPUDevice::thread_id_y();

            const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
            const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
            const DEVICE_UINT dg_sz    = dptr_mesh->m_oct_local_sz;
            const DEVICE_UINT nx       = p+1;
            const DEVICE_UINT nPe      = nx * nx * nx;

            const DEVICE_UINT blk = dptr_mesh->m_e2b_unzip_map[dptr_mesh->m_e2b_unzip_map_offset[ele] + blk_id];

            const DEVICE_UINT lx     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
            const DEVICE_UINT ly     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
            const DEVICE_UINT lz     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
            const DEVICE_UINT offset =  dptr_mesh->m_blk_list[blk].m_offset;

            const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
            const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;
            const DEVICE_UINT bLev         = patch_level; //dptr_mesh->m_all_elements[dptr_mesh->m_blk_list[blk].m_elem_begin].m_level; // or blkNode.m_level + (patch_level -blkNode.m_level);
            const DEVICE_UINT max_depth = dptr_mesh->m_max_depth;

            const TNode3d elem_node = dptr_mesh->m_all_elements[ele];

            const DEVICE_FLOAT32  hx  = (1u<<( max_depth - bLev))/(DEVICE_FLOAT32)p;
            const DEVICE_FLOAT32 xmin = blkNode.m_x - PW*hx; 
            const DEVICE_FLOAT32 ymin = blkNode.m_y - PW*hx; 
            const DEVICE_FLOAT32 zmin = blkNode.m_z - PW*hx; 

            const DEVICE_FLOAT32 xmax = blkNode.m_x + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 ymax = blkNode.m_y + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 zmax = blkNode.m_z + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 d_compar_tol=1e-6;
            
            SHARED_MEM T dgVec[p+1][p+1][p+1];
            SHARED_MEM T refEL[2][p+1][p+1];
            SHARED_MEM T v1[p+1][p+1][p+1];
            SHARED_MEM T v2[p+1][p+1][p+1];
            
            if(ele>= dptr_mesh->m_blk_list[blk].m_elem_begin && ele < dptr_mesh->m_blk_list[blk].m_elem_end)
            {
                if(!(tx < nx && ty < nx))
                    return;

                const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
                const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
                const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);
                
                for(DEVICE_UINT tz=0; tz<nx; tz++)
                    dptr_out[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)] = dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];
                
                return;
                
            }else if(elem_node.m_level == bLev)
            {

                if(!(tx < nx && ty < nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/hh;

                const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                {
                    for(DEVICE_UINT tz=0; tz< nx; tz++)
                    {
                        const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                        const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                        
                        if((kkz>=0 && kkz<lz))
                            dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];
                    }

                }

                return;

            }else if(elem_node.m_level == (bLev + 1))
            {
                if(!(tx < nx && ty < nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/(2*hh);
                const DEVICE_UINT cb = 0; //(p%2==0) ? 0 : 1;

                const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                if(tx %2 == cb && ty%2==cb)
                {
                    const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                    const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                    if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                    {
                        for(DEVICE_UINT tz=cb; tz<nx; tz+=2)
                        {
                            const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                            const DEVICE_INT kkz = std::round((zz-zmin)*invhh);

                            if((kkz>=0 && kkz<lz) )
                                dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];

                        }
                    }
                    
                }

                return;

            }else
            {
                assert(bLev == elem_node.m_level+1);
                TNode3d child_node;
                const DEVICE_UINT child_sz = (1u<<(max_depth-elem_node.m_level-1));
                if(tx < nx && ty < nx)
                {
                    for(DEVICE_UINT tz=0; tz<nx; tz++)
                        dgVec[tz][ty][tx]  = dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];

                    refEL[0][ty][tx]= refel_1d[ 0 * (p+1)*(p+1) + ty * (p+1) + tx];
                    refEL[1][ty][tx]= refel_1d[ 1 * (p+1)*(p+1) + ty * (p+1) + tx];
                }
                GPUDevice::sync_threads();

                for(DEVICE_UINT child = 0; child < NUM_CHILDREN; child++)
                {
                    const DEVICE_UINT cbit[3] = {(child & 1u), ((child & 2u)>>1u), ((child & 4u)>>2u)};
                    child_node.m_x = elem_node.m_x + cbit[0] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_y = elem_node.m_y + cbit[1] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_z = elem_node.m_z + cbit[2] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_level = elem_node.m_level+1;

                    if((child_node.m_x + child_sz < xmin || child_node.m_x >=xmax) || 
                    (child_node.m_y + child_sz < ymin || child_node.m_y >=ymax) ||
                    (child_node.m_z + child_sz < zmin || child_node.m_z >=zmax) )
                        continue;
                    
                    // interp on x
                    if((tx < nx && ty < nx))
                    {
                        for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                        {
                            DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[0]][tx][m] * dgVec[tz][ty][m];

                            v1[tz][ty][tx]=val;

                        }
                        
                    }
                            
                    GPUDevice::sync_threads();

                    // interp on y
                    if((tx < nx && ty < nx))
                    {
                        for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                        {
                            DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[1]][ty][m] * v1[tz][m][tx];
                            v2[tz][ty][tx]=val;
                        }
                        
                    }
                    
                    GPUDevice::sync_threads();

                    // interp on z
                    if((tx < nx && ty < nx))
                    {  
                        for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                        {
                            DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[2]][tz][m] * v2[m][ty][tx];

                            v1[tz][ty][tx]=val;

                        }
                    }

                    GPUDevice::sync_threads();

                    if((tx < nx && ty < nx))
                    {
                        const DEVICE_FLOAT32 hh = (1u<<(max_depth - child_node.m_level))/(DEVICE_FLOAT32) p;
                        const DEVICE_FLOAT32 invhh = 1.0/hh;

                        const DEVICE_FLOAT32 yy  = child_node.m_y + ty*hh;
                        const DEVICE_FLOAT32 xx  = child_node.m_x + tx*hh;
                        
                        const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                        const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                        if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                        {
                            
                            for(DEVICE_UINT tz=0; tz<nx; tz++)
                            {
                                const DEVICE_FLOAT32 zz  = child_node.m_z + tz*hh;
                                const DEVICE_INT kkz = std::round((zz-zmin)*invhh);

                                if((kkz>=0 && kkz<lz))
                                    dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  v1[tz][ty][tx];

                            }
                        }

                    }

                }
                
                return;

            }
 
        }

        return;
    }
    
    template<typename T, DEVICE_UINT p>
    DEVICE_FUNC
    inline void __p2c_3d__(const T refel_1d[][p+1][p+1], DEVICE_UINT cnum, const T* const dptr_in, T* const dptr_out, T* v1, T* v2)
    {
        const DEVICE_UINT tx      = GPUDevice::thread_id_x();
        const DEVICE_UINT ty      = GPUDevice::thread_id_y();
        const DEVICE_UINT cbit[3] = {(cnum & 1u), ((cnum & 2u)>>1u), ((cnum & 4u)>>2u)};
        const DEVICE_UINT nx      = p+1;

        // interp on x
        if((tx < nx && ty < nx))
        {
            for(DEVICE_UINT tz=0; tz < (p+1); tz++)
            {
                DEVICE_REAL val=0;
                for(DEVICE_UINT m=0; m< (p+1); m++)
                    val += refel_1d[cbit[0]][tx][m] * dptr_in[IDX_DG(m,ty,tz)]; //[tz][ty][m];

                v1[IDX_DG(tx,ty,tz)]=val;
            }
            
        }
                
        GPUDevice::sync_threads();

        // interp on y
        if((tx < nx && ty < nx))
        {
            for(DEVICE_UINT tz=0; tz < (p+1); tz++)
            {
                DEVICE_REAL val=0;
                for(DEVICE_UINT m=0; m< (p+1); m++)
                    val += refel_1d[cbit[1]][ty][m] * v1[IDX_DG(tx,m,tz)]; //[tz][m][tx];
                v2[IDX_DG(tx,ty,tz)]=val;
            }
            
        }
        
        GPUDevice::sync_threads();

        // interp on z
        if((tx < nx && ty < nx))
        {  
            for(DEVICE_UINT tz=0; tz < (p+1); tz++)
            {
                DEVICE_REAL val=0;
                for(DEVICE_UINT m=0; m< (p+1); m++)
                    val += refel_1d[cbit[2]][tz][m] * v2[IDX_DG(tx,ty,m)]; //[m][ty][tx];

                dptr_out[IDX_DG(tx,ty,tz)]=val;

            }
        }

        GPUDevice::sync_threads();

    }


    template<typename T, DEVICE_UINT p>
    DEVICE_FUNC
    inline void __p2c_2d__(const T refel_1d[][p+1][p+1], DEVICE_UINT cnum, DEVICE_UINT idx_base, DEVICE_UINT inc1, DEVICE_UINT inc2,
                    const T* const dptr_in, T* const dptr_out, T* v1)
    {
        const DEVICE_UINT tx      = GPUDevice::thread_id_x();
        const DEVICE_UINT ty      = GPUDevice::thread_id_y();
        const DEVICE_UINT cbit[3] = {(cnum & 1u), ((cnum & 2u)>>1u), ((cnum & 4u)>>2u)};
        const DEVICE_UINT nx      = p+1; 
        if((tx < nx && ty < nx))
        {
            DEVICE_REAL val=0;
            for(DEVICE_UINT m=0; m< (p+1); m++)
                val +=refel_1d[cbit[0]][tx][m] * dptr_in[idx_base + ty * inc2  +  m * inc1];
            v1[ty * nx  +  tx ]=val;
        }

        GPUDevice::sync_threads();

        if((tx < nx && ty < nx))
        {
            DEVICE_REAL val=0;
            for(DEVICE_UINT m=0; m< (p+1); m++)
                val +=refel_1d[cbit[1]][ty][m] * v1[m * nx  +  tx ];//[0][m][tx];
                        
            dptr_out[idx_base + ty * inc2  +  tx * inc1]=val;
        }

        //GPUDevice::sync_threads();

        return;
    }

    template<typename T, DEVICE_UINT p>
    DEVICE_FUNC
    inline void __p2c_1d__(const T refel_1d[][p+1][p+1], DEVICE_UINT cnum, DEVICE_UINT idx_base, DEVICE_UINT inc1,
                    const T* const dptr_in, T* const dptr_out)
    {
        //const DEVICE_UINT cbit[3]  = {(cnum & 1u), ((cnum & 2u)>>1u), ((cnum & 4u)>>2u)};
        // const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        // const DEVICE_UINT ty     = GPUDevice::thread_id_y();
        const DEVICE_UINT nx     = p+1;
        for(DEVICE_UINT tx =0; tx <nx; tx++)
        //if((tx < nx))
        {
            DEVICE_REAL val=0;
            for(DEVICE_UINT m=0; m< (p+1); m++)
                val +=refel_1d[cnum][tx][m] * dptr_in[idx_base + m * inc1];
            
            dptr_out[idx_base + tx * inc1]=val;
        }

        //GPUDevice::sync_threads();

        return;
    }


    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __unzip_cg__2(const MeshGPU* const dptr_mesh, const T* const dptr_in,T* const dptr_out)
    {
        const DEVICE_UINT ele    = GPUDevice::block_id_x();
        const DEVICE_UINT blk_id = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();

        if(blk_id < dptr_mesh->m_e2b_unzip_map_count[ele])
        {

            const DEVICE_UINT tx     = GPUDevice::thread_id_x();
            const DEVICE_UINT ty     = GPUDevice::thread_id_y();

            const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
            const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
            const DEVICE_UINT cg_sz    = dptr_mesh->m_oct_shared_sz;
            const DEVICE_UINT nx       = p+1;
            const DEVICE_UINT nPe      = nx * nx * nx;

            const DEVICE_UINT blk = dptr_mesh->m_e2b_unzip_map[dptr_mesh->m_e2b_unzip_map_offset[ele] + blk_id];

            const DEVICE_UINT lx     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
            const DEVICE_UINT ly     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
            const DEVICE_UINT lz     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
            const DEVICE_UINT offset =  dptr_mesh->m_blk_list[blk].m_offset;

            const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
            const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;
            const DEVICE_UINT bLev         = patch_level; 
            const DEVICE_UINT max_depth    = dptr_mesh->m_max_depth;
            const DEVICE_UINT* const e2n_cg   = dptr_mesh->m_e2n_cg;

            const TNode3d elem_node = dptr_mesh->m_all_elements[ele];

            const DEVICE_FLOAT32  hx  = (1u<<( max_depth - bLev))/(DEVICE_FLOAT32)p;
            const DEVICE_FLOAT32 xmin = blkNode.m_x - PW*hx; 
            const DEVICE_FLOAT32 ymin = blkNode.m_y - PW*hx; 
            const DEVICE_FLOAT32 zmin = blkNode.m_z - PW*hx; 

            const DEVICE_FLOAT32 xmax = blkNode.m_x + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 ymax = blkNode.m_y + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 zmax = blkNode.m_z + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 d_compar_tol=1e-6;
            
            SHARED_MEM T refEL[2][p+1][p+1];
            SHARED_MEM T dgVec[nPe];
            SHARED_MEM T cgVec[nPe];
            SHARED_MEM T v1[nPe];
            SHARED_MEM T v2[nPe];

            if(tx < nx && ty < nx)
            {
                for(DEVICE_UINT tz=0; tz<nx; tz++)
                    cgVec[IDX_DG(tx,ty,tz)]  = dptr_in[var* cg_sz + e2n_cg[ele * nPe + IDX_DG(tx,ty,tz)]];

                for(DEVICE_UINT tz=0; tz<nx; tz++)
                    dgVec[IDX_DG(tx,ty,tz)]  = cgVec[IDX_DG(tx,ty,tz)]; //dptr_in[var* cg_sz + e2n_cg[ele * nPe + IDX_DG(tx,ty,tz)]];               

                refEL[0][ty][tx]= refel_1d[ 0 * (p+1)*(p+1) + ty * (p+1) + tx];
                refEL[1][ty][tx]= refel_1d[ 1 * (p+1)*(p+1) + ty * (p+1) + tx];
            }
            
            GPUDevice::sync_threads();
            DEVICE_BOOL node_valid[3][3][3];
            for(DEVICE_UINT s3=0; s3 < 3 ; s3+=1)
                for(DEVICE_UINT s2=0; s2 < 3 ; s2+=1)
                    for(DEVICE_UINT s1=0; s1 < 3 ; s1+=1)
                        node_valid[s1][s2][s3]=false;
            
            // face interpolations. 
            for(DEVICE_UINT s1=0; s1 < 3 ; s1+=2)
            {   
                DEVICE_UINT cnum=0;
                DEVICE_BOOL is_hanging = false;
                is_hanging = is_face_hanging<p>(s1, 1, 1, dptr_mesh,ele, cnum);
                if(is_hanging)
                {
                    const DEVICE_UINT inc1     = nx;
                    const DEVICE_UINT inc2     = nx * nx;
                    const DEVICE_UINT idx_base = IDX_DG(((s1>>1u) * p),0,0);
                    
                    __p2c_2d__<DEVICE_REAL,p>(refEL,cnum,idx_base,inc1,inc2,cgVec,dgVec,v1);
                    node_valid[s1][1][1]=true;
                    for(DEVICE_UINT s2=0; s2 < 3 ; s2+=2)
                    {
                        node_valid[s1][s2][1]=true;
                        node_valid[s1][1][s2]=true;
                    }
                        

                }

                is_hanging = is_face_hanging<p> ( 1, s1 , 1, dptr_mesh,ele, cnum);
                if(is_hanging)
                {
                    const DEVICE_UINT inc1     = 1;
                    const DEVICE_UINT inc2     = nx * nx;
                    const DEVICE_UINT idx_base = IDX_DG(0,((s1>>1u) * p),0);
                    
                    __p2c_2d__<DEVICE_REAL,p>(refEL,cnum,idx_base,inc1,inc2,cgVec,dgVec,v1);
                    
                    node_valid[1][s1][1]=true;

                    for(DEVICE_UINT s2=0; s2 < 3 ; s2+=2)
                    {
                        node_valid[1][s1][s2]=true;
                        node_valid[s2][s1][1]=true;
                    }
                }

                is_hanging = is_face_hanging<p> (1, 1 , s1, dptr_mesh,ele, cnum);
                if(is_hanging)
                {
                    const DEVICE_UINT inc1     = 1;
                    const DEVICE_UINT inc2     = nx;
                    const DEVICE_UINT idx_base = IDX_DG(0,0,((s1>>1u) * p));
                    
                    __p2c_2d__<DEVICE_REAL,p>(refEL,cnum,idx_base,inc1,inc2,cgVec,dgVec,v1);

                    node_valid[1][1][s1]=true;
                    for(DEVICE_UINT s2=0; s2 < 3 ; s2+=2)
                    {
                        node_valid[1][s2][s1]=true;
                        node_valid[s2][1][s1]=true;
                    }
                }

            }

            //GPUDevice::sync_threads();

            // edge interpolations. 
            //for(DEVICE_UINT s1=0; s1 < 3 ; s1+=2)
                //for(DEVICE_UINT s2=0; s2 < 3 ; s2+=2)
                if( tx < 2 && ty < 2)
                {   
                    DEVICE_UINT s1 = tx * 2;
                    DEVICE_UINT s2 = ty * 2;
                    DEVICE_UINT cnum       = 0;
                    DEVICE_BOOL is_hanging = false;
                    
                    is_hanging = is_edge_hanging<p>(1,s1,s2, dptr_mesh,ele,cnum);
                    if( (!node_valid[1][s1][s2]) && is_hanging)
                    {
                        const DEVICE_UINT inc1     = 1;
                        const DEVICE_UINT idx_base = IDX_DG(0, ((s1>>1u) * p), ((s2>>1u) * p));
                        
                        __p2c_1d__<DEVICE_REAL,p>(refEL,cnum,idx_base,inc1,cgVec,dgVec);
                    }
                    
                    is_hanging = is_edge_hanging<p>(s1,1,s2, dptr_mesh,ele,cnum);
                    
                    if((!node_valid[s1][1][s2]) && is_hanging)
                    {
                        const DEVICE_UINT inc1     = nx;
                        const DEVICE_UINT idx_base = IDX_DG(((s1>>1u) * p), 0, ((s2>>1u) * p));
                        
                        __p2c_1d__<DEVICE_REAL,p>(refEL,cnum,idx_base,inc1,cgVec,dgVec);
                    }

                    is_hanging = is_edge_hanging<p>(s1,s2,1, dptr_mesh,ele,cnum);
                    
                    if((!node_valid[s1][s2][1]) && is_hanging)
                    {
                        const DEVICE_UINT inc1     = nx*nx;
                        const DEVICE_UINT idx_base = IDX_DG(((s1>>1u) * p),((s2>>1u) * p),0);
                        
                        __p2c_1d__<DEVICE_REAL,p>(refEL,cnum,idx_base,inc1,cgVec,dgVec);
                    }

                }

            GPUDevice::sync_threads();

            if(ele>= dptr_mesh->m_blk_list[blk].m_elem_begin && ele < dptr_mesh->m_blk_list[blk].m_elem_end)
            {
                if(!(tx < nx && ty < nx))
                    return;

                const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
                const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
                const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);
                
                for(DEVICE_UINT tz=0; tz<nx; tz++)
                    dptr_out[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)] = dgVec[IDX_DG(tx,ty,tz)];
                
                return;
                
            }else if(elem_node.m_level == bLev)
            {

                if(!(tx < nx && ty < nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/hh;

                const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                {
                    for(DEVICE_UINT tz=0; tz< nx; tz++)
                    {
                        const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                        const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                        
                        if((kkz>=0 && kkz<lz))
                            dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dgVec[IDX_DG(tx,ty,tz)];
                    }

                }

                return;

            }else if(elem_node.m_level == (bLev + 1))
            {
                if(!(tx < nx && ty < nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/(2*hh);
                const DEVICE_UINT cb = 0; //(p%2==0) ? 0 : 1;

                const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                if(tx %2 == cb && ty%2==cb)
                {
                    const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                    const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                    if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                    {
                        for(DEVICE_UINT tz=cb; tz<nx; tz+=2)
                        {
                            const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                            const DEVICE_INT kkz = std::round((zz-zmin)*invhh);

                            if((kkz>=0 && kkz<lz) )
                                dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dgVec[IDX_DG(tx,ty,tz)];

                        }
                    }
                    
                }

                return;

            }else
            {
                assert(bLev == elem_node.m_level+1);
                TNode3d child_node;
                const DEVICE_UINT child_sz = (1u<<(max_depth-elem_node.m_level-1));
                
                for(DEVICE_UINT child = 0; child < NUM_CHILDREN; child++)
                {
                    const DEVICE_UINT cbit[3] = {(child & 1u), ((child & 2u)>>1u), ((child & 4u)>>2u)};
                    child_node.m_x = elem_node.m_x + cbit[0] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_y = elem_node.m_y + cbit[1] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_z = elem_node.m_z + cbit[2] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_level = elem_node.m_level+1;

                    if((child_node.m_x + child_sz < xmin || child_node.m_x >=xmax) || 
                    (child_node.m_y + child_sz < ymin || child_node.m_y >=ymax) ||
                    (child_node.m_z + child_sz < zmin || child_node.m_z >=zmax) )
                        continue;
                    
                    __p2c_3d__<DEVICE_REAL,p>(refEL,child,dgVec,cgVec,v1,v2);
                    //std::swap(dgVec,cgVec);
                    
                    if((tx < nx && ty < nx))
                    {
                        const DEVICE_FLOAT32 hh = (1u<<(max_depth - child_node.m_level))/(DEVICE_FLOAT32) p;
                        const DEVICE_FLOAT32 invhh = 1.0/hh;

                        const DEVICE_FLOAT32 yy  = child_node.m_y + ty*hh;
                        const DEVICE_FLOAT32 xx  = child_node.m_x + tx*hh;
                        
                        const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                        const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                        if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                        {
                            
                            for(DEVICE_UINT tz=0; tz<nx; tz++)
                            {
                                const DEVICE_FLOAT32 zz  = child_node.m_z + tz*hh;
                                const DEVICE_INT kkz = std::round((zz-zmin)*invhh);

                                if((kkz>=0 && kkz<lz))
                                    dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  cgVec[IDX_DG(tx,ty,tz)];

                            }
                        }

                    }

                }
                
                return;

            }
 
        }

        return;
    }

    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    void __unzip_cg1__2(const MeshGPU* const dptr_mesh, const T* const dptr_in,T* const dptr_out)
    {
        const DEVICE_UINT ele    = GPUDevice::block_id_x();
        //const DEVICE_UINT blk_id = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();
        
        const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        const DEVICE_UINT ty     = GPUDevice::thread_id_y();

        const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
        const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
        const DEVICE_UINT cg_sz    = dptr_mesh->m_oct_shared_sz;
        const DEVICE_UINT nx       = p+1;
        const DEVICE_UINT nPe      = nx * nx * nx;

        SHARED_MEM T refEL[2][p+1][p+1];
        //SHARED_MEM T dgVec[nPe];
        //SHARED_MEM T v1[nPe];
        SHARED_MEM T v2[nPe];
        //v3[0] interp in , v3[1] interp out
        SHARED_MEM T v3[2][nPe];

        // this will increase shared memory leading the low occupancy
        //SHARED_MEM T u_refine[8][nPe];
        //bool is_refine[8]={false, false, false, false, false, false, false, false};

        T* dgVec = v3[0];
        T* v1    = v3[1];
        
        
        const DEVICE_UINT max_depth    = dptr_mesh->m_max_depth;
        const DEVICE_UINT* const e2n_cg   = dptr_mesh->m_e2n_cg;
        const TNode3d elem_node = dptr_mesh->m_all_elements[ele];
        const DEVICE_FLOAT32 d_compar_tol=1e-6;

        HangingInfo hinfo_ele = dptr_mesh->m_hinfo[ele];
        DEVICE_UINT hflag     = hinfo_ele.hflag;

        //if(hflag!=0 || bLev > elem_node.m_level)
        {
            if(tx < nx && ty < nx)
            {
                refEL[0][ty][tx]= refel_1d[ 0 * (p+1)*(p+1) + ty * (p+1) + tx];
                refEL[1][ty][tx]= refel_1d[ 1 * (p+1)*(p+1) + ty * (p+1) + tx];
            }
        }

        GPUDevice::sync_threads();

        if(hflag !=0)
        {
            if(tx < nx && ty < nx)
                for(DEVICE_INT tz=0; tz<nx; tz++)
                    v3[0][IDX_DG(tx,ty,tz)]  = dptr_in[var* cg_sz + e2n_cg[ele * nPe + IDX_DG(tx,ty,tz)]];
            
            GPUDevice::sync_threads();    
            DEVICE_UINT mid_bit = max_depth - elem_node.m_level;
            DEVICE_UINT cnum    = ((((elem_node.m_z >> mid_bit) & 1u) << 2u) | (((elem_node.m_y >> mid_bit) & 1u) << 1u) | ((elem_node.m_x >>mid_bit) & 1u));
            __p2c_3d__<DEVICE_REAL,p>(refEL,cnum,v3[0],v3[1],v1,v2);

            if(tx < nx && ty < nx)
            {
                DEVICE_INT ele_level= elem_node.m_level; 
                DEVICE_INT idx      = IDX_DG(tx,ty,0);
                DEVICE_INT owner_id = (dptr_mesh->m_e2n_dg[ele * (nPe) + idx]/nPe); 
                DEVICE_INT l        = dptr_mesh->m_all_elements[owner_id].m_level < ele_level;
                dgVec[idx]          = v3[l][idx];

                idx                 = IDX_DG(tx,ty,nx-1);
                owner_id = (dptr_mesh->m_e2n_dg[ele * (nPe) + idx]/nPe); 
                l        = dptr_mesh->m_all_elements[owner_id].m_level <  ele_level;
                dgVec[idx]      = v3[l][idx];

                idx = IDX_DG(tx,0,ty);
                owner_id = (dptr_mesh->m_e2n_dg[ele * (nPe) + idx]/nPe); 
                l        = dptr_mesh->m_all_elements[owner_id].m_level <  ele_level;
                dgVec[idx]      = v3[l][idx];

                idx = IDX_DG(tx,nx-1,ty);
                owner_id = (dptr_mesh->m_e2n_dg[ele * (nPe) + idx]/nPe); 
                l        = dptr_mesh->m_all_elements[owner_id].m_level <  ele_level;
                dgVec[idx]      = v3[l][idx];

                idx = IDX_DG(0,tx,ty);
                owner_id = (dptr_mesh->m_e2n_dg[ele * (nPe) + idx]/nPe); 
                l        = dptr_mesh->m_all_elements[owner_id].m_level <  ele_level;
                dgVec[idx]      = v3[l][idx];

                idx = IDX_DG(nx-1,tx,ty);
                owner_id = (dptr_mesh->m_e2n_dg[ele * (nPe) + idx]/nPe); 
                l        = dptr_mesh->m_all_elements[owner_id].m_level <  ele_level;
                dgVec[idx]      = v3[l][idx];
            }
            
            
        }else
        {
            if(tx < nx && ty < nx)
                for(DEVICE_INT tz=0; tz<nx; tz++)
                    dgVec[IDX_DG(tx,ty,tz)]  = dptr_in[var* cg_sz + e2n_cg[ele * nPe + IDX_DG(tx,ty,tz)]];
        }

        GPUDevice::sync_threads();


        //if(blk_id < dptr_mesh->m_e2b_unzip_map_count[ele])
        for(DEVICE_INT blk_id=0; blk_id < dptr_mesh->m_e2b_unzip_map_count[ele]; blk_id++)
        {

            const DEVICE_UINT blk          = dptr_mesh->m_e2b_unzip_map[dptr_mesh->m_e2b_unzip_map_offset[ele] + blk_id];
            const DEVICE_UINT lx           =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
            const DEVICE_UINT ly           =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
            const DEVICE_UINT lz           =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
            const DEVICE_UINT offset       =  dptr_mesh->m_blk_list[blk].m_offset;

            const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
            const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;
            const DEVICE_UINT bLev         = patch_level; 
            
            const DEVICE_FLOAT32  hx  = (1u<<( max_depth - bLev))/(DEVICE_FLOAT32)p;
            const DEVICE_FLOAT32 xmin = blkNode.m_x - PW*hx; 
            const DEVICE_FLOAT32 ymin = blkNode.m_y - PW*hx; 
            const DEVICE_FLOAT32 zmin = blkNode.m_z - PW*hx; 

            const DEVICE_FLOAT32 xmax = blkNode.m_x + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 ymax = blkNode.m_y + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 zmax = blkNode.m_z + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            
            
            
            /*if(ele>= dptr_mesh->m_blk_list[blk].m_elem_begin && ele < dptr_mesh->m_blk_list[blk].m_elem_end)
            {
                if(!(tx < nx && ty < nx))
                    return;

                const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
                const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
                const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);
                
                for(DEVICE_UINT tz=0; tz<nx; tz++)
                    dptr_out[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)] = dgVec[IDX_DG(tx,ty,tz)];
                
                return;
                
            }else*/ 
            if(elem_node.m_level == bLev)
            {

                // if(!(tx < nx && ty < nx))
                //     return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/hh;

                const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                {
                    for(DEVICE_UINT tz=0; tz< nx; tz++)
                    {
                        const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                        const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                        
                        if((kkz>=0 && kkz<lz))
                            dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dgVec[IDX_DG(tx,ty,tz)];
                    }

                }

                //return;

            }else if(elem_node.m_level == (bLev + 1))
            {
                // if(!(tx < nx && ty < nx))
                //     return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/(2*hh);
                const DEVICE_INT cb = 0; //(p%2==0) ? 0 : 1;

                const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;

                if(tx %2 == cb && ty%2==cb)
                {
                    const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                    const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                    if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                    {
                        for(DEVICE_INT tz=cb; tz<nx; tz+=2)
                        {
                            const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                            const DEVICE_INT kkz = std::round((zz-zmin)*invhh);

                            if((kkz>=0 && kkz<lz) )
                                dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dgVec[IDX_DG(tx,ty,tz)];

                        }
                    }
                    
                }

                //return;

            }else
            {
                assert(bLev == elem_node.m_level+1);
                TNode3d child_node;
                const DEVICE_UINT child_sz = (1u<<(max_depth-elem_node.m_level-1));
                
                for(DEVICE_UINT child = 0; child < NUM_CHILDREN; child++)
                {
                    const DEVICE_UINT cbit[3] = {(child & 1u), ((child & 2u)>>1u), ((child & 4u)>>2u)};
                    child_node.m_x = elem_node.m_x + cbit[0] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_y = elem_node.m_y + cbit[1] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_z = elem_node.m_z + cbit[2] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_level = elem_node.m_level+1;

                    if((child_node.m_x + child_sz < xmin || child_node.m_x >=xmax) || 
                    (child_node.m_y + child_sz < ymin || child_node.m_y >=ymax) ||
                    (child_node.m_z + child_sz < zmin || child_node.m_z >=zmax) )
                        continue;
                    
                    // if(!is_refine[child])
                    // {
                    //     __p2c_3d__<DEVICE_REAL,p>(refEL,child,dgVec,u_refine[child],v1,v2);
                    //     is_refine[child]=true;
                    // }
                    __p2c_3d__<DEVICE_REAL,p>(refEL,child,dgVec,v3[1],v1,v2);
                    
                    
                    if((tx < nx && ty < nx))
                    {
                        const DEVICE_FLOAT32 hh = (1u<<(max_depth - child_node.m_level))/(DEVICE_FLOAT32) p;
                        const DEVICE_FLOAT32 invhh = 1.0/hh;

                        const DEVICE_FLOAT32 yy  = child_node.m_y + ty*hh;
                        const DEVICE_FLOAT32 xx  = child_node.m_x + tx*hh;
                        
                        const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                        const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                        if(((jjy>=0 && jjy<ly) && (iix>=0 && iix<lx)))
                        {
                            
                            for(DEVICE_INT tz=0; tz<nx; tz++)
                            {
                                const DEVICE_FLOAT32 zz  = child_node.m_z + tz*hh;
                                const DEVICE_INT kkz = std::round((zz-zmin)*invhh);

                                // if((kkz>=0 && kkz<lz))
                                //     dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  u_refine[child][IDX_DG(tx,ty,tz)];
                                if((kkz>=0 && kkz<lz))
                                    dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  v3[1][IDX_DG(tx,ty,tz)];

                            }
                        }

                    }

                }
                
                //return;

            }
 
        }

        return;
    }
    

    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(8)
    void __unzip_dg__1(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        const DEVICE_UINT ele    = GPUDevice::block_id_x();
        const DEVICE_UINT blk_id = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();

        if(blk_id < dptr_mesh->m_e2b_unzip_map_count[ele])
        {

            const DEVICE_UINT tx     = GPUDevice::thread_id_x();
            //const DEVICE_UINT ty     = GPUDevice::thread_id_y();

            const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
            const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
            const DEVICE_UINT dg_sz    = dptr_mesh->m_oct_local_sz;
            const DEVICE_UINT nx       = p+1;
            const DEVICE_UINT nPe      = nx * nx * nx;

            const DEVICE_UINT blk = dptr_mesh->m_e2b_unzip_map[dptr_mesh->m_e2b_unzip_map_offset[ele] + blk_id];

            const DEVICE_UINT lx     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
            const DEVICE_UINT ly     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
            const DEVICE_UINT lz     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
            const DEVICE_UINT offset =  dptr_mesh->m_blk_list[blk].m_offset;

            const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
            const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;
            const DEVICE_UINT bLev         = patch_level; //dptr_mesh->m_all_elements[dptr_mesh->m_blk_list[blk].m_elem_begin].m_level; // or blkNode.m_level + (patch_level -blkNode.m_level);
            const DEVICE_UINT max_depth = dptr_mesh->m_max_depth;

            const TNode3d elem_node = dptr_mesh->m_all_elements[ele];

            const DEVICE_FLOAT32  hx  = (1u<<( max_depth - bLev))/(DEVICE_FLOAT32)p;
            const DEVICE_FLOAT32 xmin = blkNode.m_x - PW*hx; 
            const DEVICE_FLOAT32 ymin = blkNode.m_y - PW*hx; 
            const DEVICE_FLOAT32 zmin = blkNode.m_z - PW*hx; 

            const DEVICE_FLOAT32 xmax = blkNode.m_x + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 ymax = blkNode.m_y + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 zmax = blkNode.m_z + (1u<<(max_depth-blkNode.m_level)) + PW*hx;
            const DEVICE_FLOAT32 d_compar_tol=1e-6;
            
            SHARED_MEM T dgVec[p+1][p+1][p+1];
            SHARED_MEM T refEL[2][p+1][p+1];
            SHARED_MEM T v1[p+1][p+1][p+1];
            SHARED_MEM T v2[p+1][p+1][p+1];
            
            if(ele>= dptr_mesh->m_blk_list[blk].m_elem_begin && ele < dptr_mesh->m_blk_list[blk].m_elem_end)
            {
                if(!(tx < nx))
                    return;

                const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
                const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
                const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);
                
                for(DEVICE_UINT tz=0; tz<nx; tz++)
                    for(DEVICE_UINT ty=0; ty<nx; ty++)
                        dptr_out[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)] = dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];
                
                return;
                
            }else if(elem_node.m_level == bLev)
            {

                if(!(tx < nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/hh;

                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;
                const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                if((iix>=0 && iix<lx))
                {
                    for(DEVICE_UINT tz=0; tz< nx; tz++)
                    for(DEVICE_UINT ty=0; ty< nx; ty++)
                    {
                        const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                        const DEVICE_INT kkz = std::round((zz-zmin)*invhh);

                        const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;
                        const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                        
                        
                        if((kkz>=0 && kkz<lz) && (jjy>=0 && jjy<ly))
                            dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];
                    }

                }

                return;

            }else if(elem_node.m_level == (bLev + 1))
            {
                if(!(tx < nx))
                    return;

                const DEVICE_FLOAT32 hh = (1u<<(max_depth - elem_node.m_level))/(DEVICE_FLOAT32) p;
                const DEVICE_FLOAT32 invhh = 1.0/(2*hh);
                const DEVICE_UINT cb = 0; //(p%2==0) ? 0 : 1;

                const DEVICE_FLOAT32 xx  = elem_node.m_x + tx*hh;
                const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                
                if(tx %2 == cb)
                {
                    if((iix>=0 && iix<lx))
                    {
                        for(DEVICE_UINT tz=cb; tz<nx; tz+=2)
                        for(DEVICE_UINT ty=cb; ty<nx; ty+=2)
                        {
                            const DEVICE_FLOAT32 zz  = elem_node.m_z + tz*hh;
                            const DEVICE_FLOAT32 yy  = elem_node.m_y + ty*hh;

                            const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                            const DEVICE_INT jjy = std::round((yy-ymin)*invhh);

                            if((jjy>=0 && jjy<ly) && (kkz>=0 && kkz<lz) )
                                dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];

                        }
                    }
                    
                }

                return;

            }else
            {
                assert(bLev == elem_node.m_level+1);
                TNode3d child_node;
                const DEVICE_UINT child_sz = (1u<<(max_depth-elem_node.m_level-1));
                if(tx < nx)
                {
                    for(DEVICE_UINT tz=0; tz<nx; tz++)
                        for(DEVICE_UINT ty=0; ty<nx; ty++)
                            dgVec[tz][ty][tx]  = dptr_in[var* dg_sz + ele * nPe + IDX_DG(tx,ty,tz)];

                    for(DEVICE_UINT ty=0; ty<nx; ty++)
                    {
                        refEL[0][ty][tx]= refel_1d[ 0 * (p+1)*(p+1) + ty * (p+1) + tx];
                        refEL[1][ty][tx]= refel_1d[ 1 * (p+1)*(p+1) + ty * (p+1) + tx];
                    }
                    
                }
                GPUDevice::sync_threads();

                for(DEVICE_UINT child = 0; child < NUM_CHILDREN; child++)
                {
                    const DEVICE_UINT cbit[3] = {(child & 1u), ((child & 2u)>>1u), ((child & 4u)>>2u)};
                    child_node.m_x = elem_node.m_x + cbit[0] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_y = elem_node.m_y + cbit[1] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_z = elem_node.m_z + cbit[2] * (1u<<(max_depth-elem_node.m_level-1));
                    child_node.m_level = elem_node.m_level+1;

                    if((child_node.m_x + child_sz < xmin || child_node.m_x >=xmax) || 
                    (child_node.m_y + child_sz < ymin || child_node.m_y >=ymax) ||
                    (child_node.m_z + child_sz < zmin || child_node.m_z >=zmax) )
                        continue;
                    
                    // interp on x
                    if((tx < nx))
                    {
                        for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                        for(DEVICE_UINT ty=0; ty < (p+1); ty++)
                        {
                            DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[0]][tx][m] * dgVec[tz][ty][m];

                            v1[tz][ty][tx]=val;

                        }
                        
                    }
                            
                    GPUDevice::sync_threads();

                    // interp on y
                    if((tx < nx))
                    {
                        for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                        for(DEVICE_UINT ty=0; ty < (p+1); ty++)
                        {
                            DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[1]][ty][m] * v1[tz][m][tx];
                            v2[tz][ty][tx]=val;
                        }
                        
                    }
                    
                    GPUDevice::sync_threads();

                    // interp on z
                    if((tx < nx))
                    {  
                        for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                        for(DEVICE_UINT ty=0; ty < (p+1); ty++)
                        {
                            DEVICE_REAL val=0;
                            for(DEVICE_UINT m=0; m< (p+1); m++)
                                val += refEL[cbit[2]][tz][m] * v2[m][ty][tx];

                            v1[tz][ty][tx]=val;

                        }
                    }

                    GPUDevice::sync_threads();

                    if((tx < nx))
                    {
                        const DEVICE_FLOAT32 hh = (1u<<(max_depth - child_node.m_level))/(DEVICE_FLOAT32) p;
                        const DEVICE_FLOAT32 invhh = 1.0/hh;

                        
                        const DEVICE_FLOAT32 xx  = child_node.m_x + tx*hh;
                        const DEVICE_INT iix = std::round((xx-xmin)*invhh);

                        if(((iix>=0 && iix<lx)))
                        {
                            
                            for(DEVICE_UINT tz=0; tz<nx; tz++)
                            for(DEVICE_UINT ty=0; ty<nx; ty++)
                            {
                                const DEVICE_FLOAT32 yy  = child_node.m_y + ty*hh;
                                const DEVICE_FLOAT32 zz  = child_node.m_z + tz*hh;
                                
                                
                                const DEVICE_INT jjy = std::round((yy-ymin)*invhh);
                                const DEVICE_INT kkz = std::round((zz-zmin)*invhh);
                                

                                if( (jjy>=0 && jjy<ly) && (kkz>=0 && kkz<lz))
                                    dptr_out[var*unzip_sz + offset + kkz*lx*ly + jjy*lx + iix] =  v1[tz][ty][tx];

                            }
                        }

                    }

                }
                
                return;

            }
 
        }

        return;
    }

    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __zip_dg_enforce_c0__(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        const DEVICE_UINT ele_id = GPUDevice::block_id_x();
        const DEVICE_UINT blk    = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();
        
        // implies no contribution for the unzip. 
        const DEVICE_UINT num_elems_blk = (dptr_mesh->m_blk_list[blk].m_elem_end - dptr_mesh->m_blk_list[blk].m_elem_begin);
        if(num_elems_blk<=ele_id)
            return;

        const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        const DEVICE_UINT ty     = GPUDevice::thread_id_y();
        //const DEVICE_UINT tz   = GPUDevice::thread_id_z();
        const DEVICE_UINT nx     = (p+1);

        if(tx < nx && ty < nx)
        {
            const DEVICE_UINT ele    = dptr_mesh->m_blk_list[blk].m_elem_begin + ele_id;
            const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
            const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
            const DEVICE_UINT dg_sz    = dptr_mesh->m_oct_local_sz;
            const DEVICE_UINT nPe      = nx * nx * nx; 

            for(DEVICE_UINT tz=0; tz < (p+1); tz+=p)
            {
                const DEVICE_UINT dg_index = dptr_mesh->m_e2n_dg[ele  * nPe + IDX_DG(tx,ty,tz)];
                const DEVICE_UINT dg_ele   = ele * nPe + IDX_DG(tx,ty,tz);
                assert(dg_index!=LOOK_UP_TABLE_DEFAULT);
                if(dg_index != dg_ele)
                {
                    const DEVICE_UINT owner = dg_index / nPe;
                    if(dptr_mesh->m_all_elements[owner].m_level == dptr_mesh->m_all_elements[ele].m_level)
                    {
                        DEVICE_UINT ee, k_ee, j_ee, i_ee;
                        dg2eijk(p, dg_index, ee, i_ee, j_ee, k_ee);
                        dptr_out[var * dg_sz + ele*nPe + IDX_DG(tx,ty,tz)] = dptr_out[var * dg_sz + ee*nPe + IDX_DG(i_ee,j_ee,k_ee)];
                    }
                }

            }
        
        }

        return;
            
    }

    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __zip_dg__(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        const DEVICE_UINT ele_id = GPUDevice::block_id_x();
        const DEVICE_UINT blk    = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();
        
        // implies no contribution for the unzip. 
        const DEVICE_UINT num_elems_blk = (dptr_mesh->m_blk_list[blk].m_elem_end - dptr_mesh->m_blk_list[blk].m_elem_begin);
        if(num_elems_blk<=ele_id)
            return;

        const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        const DEVICE_UINT ty     = GPUDevice::thread_id_y();
        //const DEVICE_UINT tz   = GPUDevice::thread_id_z();
        const DEVICE_UINT nx     = (p+1);
        if(tx < nx && ty < nx)
        {
            const DEVICE_UINT ele    = dptr_mesh->m_blk_list[blk].m_elem_begin + ele_id;
            const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
            const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
            const DEVICE_UINT dg_sz    = dptr_mesh->m_oct_local_sz;
            const DEVICE_UINT nPe      = nx * nx * nx; 

            const DEVICE_UINT lx     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
            const DEVICE_UINT ly     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
            const DEVICE_UINT lz     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
            const DEVICE_UINT offset =  dptr_mesh->m_blk_list[blk].m_offset;
            
            const TNode3d elem_node        = dptr_mesh->m_all_elements[ele]; 
            const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
            const DEVICE_UINT max_depth    = dptr_mesh->m_max_depth;
            const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;

            const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
            const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
            const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);

            for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                dptr_out[var* dg_sz + ele*nPe + IDX_DG(tx,ty,tz)] = dptr_in[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)];
                

        }
        return;


    }

    template<typename T, DEVICE_UINT p, DEVICE_UINT PW>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __zip_cg__(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        const DEVICE_UINT ele_id = GPUDevice::block_id_x();
        const DEVICE_UINT blk    = GPUDevice::block_id_y();
        const DEVICE_UINT var    = GPUDevice::block_id_z();
        
        // implies no contribution for the unzip. 
        const DEVICE_UINT num_elems_blk = (dptr_mesh->m_blk_list[blk].m_elem_end - dptr_mesh->m_blk_list[blk].m_elem_begin);
        if(num_elems_blk<=ele_id)
            return;

        const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        const DEVICE_UINT ty     = GPUDevice::thread_id_y();
        //const DEVICE_UINT tz   = GPUDevice::thread_id_z();
        const DEVICE_UINT nx     = (p+1);
        if(tx < nx && ty < nx)
        {
            const DEVICE_UINT ele    = dptr_mesh->m_blk_list[blk].m_elem_begin + ele_id;
            const DEVICE_UINT porder   = dptr_mesh->m_ele_order;
            const DEVICE_UINT unzip_sz = dptr_mesh->m_oct_unzip_sz;
            const DEVICE_UINT cg_sz    = dptr_mesh->m_oct_shared_sz;
            const DEVICE_UINT nPe      = nx * nx * nx; 

            const DEVICE_UINT lx     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[0];
            const DEVICE_UINT ly     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[1];
            const DEVICE_UINT lz     =  dptr_mesh->m_blk_list[blk].m_aligned_sz[2];
            const DEVICE_UINT offset =  dptr_mesh->m_blk_list[blk].m_offset;
            
            const TNode3d elem_node        = dptr_mesh->m_all_elements[ele]; 
            const TNode3d blkNode          = dptr_mesh->m_blk_list[blk].m_blockNode;
            const DEVICE_UINT max_depth    = dptr_mesh->m_max_depth;
            const DEVICE_UINT patch_level  = dptr_mesh->m_blk_list[blk].m_patch_level;

            const DEVICE_UINT ei=(elem_node.m_x-blkNode.m_x)>>(max_depth-patch_level);
            const DEVICE_UINT ej=(elem_node.m_y-blkNode.m_y)>>(max_depth-patch_level);
            const DEVICE_UINT ek=(elem_node.m_z-blkNode.m_z)>>(max_depth-patch_level);

            for(DEVICE_UINT tz=0; tz < (p+1); tz++)
                if((dptr_mesh->m_e2n_dg[ele*nPe + IDX_DG(tx,ty,tz)]/nPe)==ele)
                    dptr_out[var* cg_sz + dptr_mesh->m_e2n_cg[ele*nPe + IDX_DG(tx,ty,tz)]] = dptr_in[var * unzip_sz + offset + (ek*p+tz+PW) * (ly * lx) +(ej * p + ty + PW) * (lx) + (ei*p + tx + PW)];
                

        }
        return;


    }


    template<typename T>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __write_sb_dg__(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        const DEVICE_UINT proc_id   = GPUDevice::block_id_x();
        const DEVICE_UINT ele_id    = GPUDevice::block_id_y();
        const DEVICE_UINT var       = GPUDevice::block_id_z();
        
        const DEVICE_UINT dof       = GPUDevice::grid_dim_z();
        const DEVICE_UINT npes      = GPUDevice::grid_dim_x();

        if(dptr_mesh->m_elem_send_c[proc_id]<=ele_id)
            return;


        const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        const DEVICE_UINT ty     = GPUDevice::thread_id_y();
        //const DEVICE_UINT tz   = GPUDevice::thread_id_z();
        const DEVICE_UINT nx     = dptr_mesh->m_ele_order + 1;
        if(tx < nx && ty < nx)
        {

            const DEVICE_UINT nPe     = nx*nx*nx;
            
            const DEVICE_UINT dg_sz  = dptr_mesh->m_oct_local_sz;
            const DEVICE_UINT lb      = dptr_mesh->m_num_pre_ghost_elements;
            
            const DEVICE_UINT sf     = dptr_mesh->m_elem_send_f[proc_id];
            const DEVICE_UINT sc     = dptr_mesh->m_elem_send_c[proc_id];
            const DEVICE_UINT sm_ele = dptr_mesh->m_elem_send_sm[sf+ele_id] +lb;

            for(DEVICE_UINT tz=0; tz < nx; tz++)
                dptr_out[dof * nPe * sf + var * nPe * sc + ele_id * nPe + IDX_DG(tx,ty,tz) ] = dptr_in[var * dg_sz + sm_ele*nPe + IDX_DG(tx,ty,tz)];

        }

        return;
        
    }

    template<typename T>
    GLOBAL_FUNC 
    __launch_bounds__(64)
    void __read_rb_dg__(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {

        const DEVICE_UINT proc_id   = GPUDevice::block_id_x();
        const DEVICE_UINT ele_id    = GPUDevice::block_id_y();
        const DEVICE_UINT var       = GPUDevice::block_id_z();
        
        const DEVICE_UINT dof       = GPUDevice::grid_dim_z();
        const DEVICE_UINT npes      = GPUDevice::grid_dim_x();

        if(dptr_mesh->m_elem_recv_c[proc_id]<=ele_id)
            return;

        const DEVICE_UINT tx     = GPUDevice::thread_id_x();
        const DEVICE_UINT ty     = GPUDevice::thread_id_y();
        const DEVICE_UINT nx     = dptr_mesh->m_ele_order + 1;

        if(tx < nx && ty < nx)
        {
            const DEVICE_UINT nPe     = nx*nx*nx;

            const DEVICE_UINT dg_sz  = dptr_mesh->m_oct_local_sz;
            const DEVICE_UINT lb     = dptr_mesh->m_num_pre_ghost_elements;
            
            const DEVICE_UINT rf     = dptr_mesh->m_elem_recv_f[proc_id];
            const DEVICE_UINT rc     = dptr_mesh->m_elem_recv_c[proc_id];
            const DEVICE_UINT rm_ele = dptr_mesh->m_elem_recv_sm[rf+ele_id];

            for(DEVICE_UINT tz=0; tz < nx; tz++)
                dptr_out[var * dg_sz + rm_ele * nPe + IDX_DG(tx,ty,tz)] = dptr_in[dof * rf *nPe + var*nPe*rc + ele_id*nPe + IDX_DG(tx,ty,tz) ];

        }
        return;

    }

    template<typename T>
    GLOBAL_FUNC 
    __launch_bounds__(GPU_MAX_THREADS_PER_BLOCK)
    void __write_sb_cg__(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {
        const DEVICE_UINT proc_id   = GPUDevice::block_id_y();
        const DEVICE_UINT var       = GPUDevice::block_id_z();
        
        const DEVICE_UINT dof       = GPUDevice::grid_dim_z();
        const DEVICE_UINT npes      = GPUDevice::grid_dim_x();

        const DEVICE_UINT sf        = dptr_mesh->m_cg_send_f[proc_id];
        const DEVICE_UINT sc        = dptr_mesh->m_cg_send_c[proc_id];
        const DEVICE_UINT cg_sz     = dptr_mesh->m_oct_shared_sz;

        const DEVICE_UINT stride = GPUDevice::block_dim_x();
        const DEVICE_UINT t_id   = GPUDevice::thread_id_x() + GPUDevice::block_id_x() * stride;
        for(unsigned int k=sf + t_id; k < sf + sc; k+=stride)
        {
            const DEVICE_UINT sm_sm       = dptr_mesh->m_cg_send_sm[k];
            dptr_out[dof * sf + var * sc + k-sf]= dptr_in[var * cg_sz + sm_sm];
        }

        return;
        
    }

    template<typename T>
    GLOBAL_FUNC 
    __launch_bounds__(GPU_MAX_THREADS_PER_BLOCK)
    void __read_rb_cg__(const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out)
    {

        const DEVICE_UINT proc_id   = GPUDevice::block_id_y();
        const DEVICE_UINT var       = GPUDevice::block_id_z();
        
        const DEVICE_UINT dof       = GPUDevice::grid_dim_z();
        const DEVICE_UINT npes      = GPUDevice::grid_dim_x();

        const DEVICE_UINT rf        = dptr_mesh->m_cg_recv_f[proc_id];
        const DEVICE_UINT rc        = dptr_mesh->m_cg_recv_c[proc_id];
        const DEVICE_UINT cg_sz     = dptr_mesh->m_oct_shared_sz;

        const DEVICE_UINT stride = GPUDevice::block_dim_x();
        const DEVICE_UINT t_id   = GPUDevice::thread_id_x() + GPUDevice::block_id_x() * stride;

        for(unsigned int k=rf + t_id; k < rf + rc; k+=stride)
        {
            const DEVICE_UINT rm_sm       = dptr_mesh->m_cg_recv_sm[k];
            dptr_out[var * cg_sz + rm_sm] = dptr_in[dof * rf + var * rc + k-rf];
        }

        return;

    }

    template<typename T, vec_type vtype>
    T* MeshGPU::createVector(DEVICE_UINT dof, bool is_local)
    {
        T* vec=nullptr;
        if(vtype == vec_type::device)
            (!is_local) ? vec = GPUDevice::device_malloc<T>(dof* m_oct_shared_sz) : vec =GPUDevice::device_malloc<T>(dof* m_num_local_nodes);
        else if(vtype == vec_type::host)
            (!is_local) ? vec = new T[dof* m_oct_shared_sz] : vec = new T[dof* m_num_local_nodes];
        else if(vtype == vec_type::host_pinned)
            (!is_local) ? vec = GPUDevice::host_malloc<T>(dof* m_oct_shared_sz) : vec =GPUDevice::host_malloc<T>(dof* m_num_local_nodes);
        
        return vec;
    }

    template<typename T, vec_type vtype>
    T* MeshGPU::createUnZippedVector(DEVICE_UINT dof, bool is_local)
    {
        T* vec=nullptr;
        if(vtype == vec_type::device)
            (!is_local) ? vec = GPUDevice::device_malloc<T>(dof* m_oct_unzip_sz) : vec =GPUDevice::device_malloc<T>(dof* m_oct_unzip_sz);
        else if(vtype == vec_type::host)
            (!is_local) ? vec = new T[dof* m_oct_unzip_sz] : vec = new T[dof* m_oct_unzip_sz];
        else if(vtype == vec_type::host_pinned)
            (!is_local) ? vec = GPUDevice::host_malloc<T>(dof* m_oct_unzip_sz) : vec =GPUDevice::host_malloc<T>(dof* m_oct_unzip_sz);
            
        return vec;
    }


    template<typename T, vec_type vtype>
    T* MeshGPU::createDGVector(DEVICE_UINT dof, bool is_local)
    {
        T* vec=nullptr;
        const unsigned int nPe  = (m_ele_order+1) * (m_ele_order+1) * (m_ele_order+1);
        if(vtype == vec_type::device)
            (!is_local) ? vec = GPUDevice::device_malloc<T>(dof* m_oct_local_sz) : vec =GPUDevice::device_malloc<T>(dof * m_num_local_elements * nPe);
        else if(vtype == vec_type::host)
            (!is_local) ? vec = new T[dof* m_oct_local_sz] : vec = new T[dof * m_num_local_elements * nPe ];
        else if(vtype == vec_type::host_pinned)
            (!is_local) ? vec = GPUDevice::host_malloc<T>(dof* m_oct_local_sz) : vec =GPUDevice::host_malloc<T>(dof * m_num_local_elements * nPe);
            
        return vec;
    }


    template<typename T, vec_type vtype>
    T* MeshGPU::createElementVector(DEVICE_UINT dof, bool is_local)
    {
        T* vec=nullptr;
        const unsigned int num_all_elements  = m_num_pre_ghost_elements + m_num_local_elements + m_num_post_ghost_elements;
        if(vtype == vec_type::device)
            (!is_local) ? vec = GPUDevice::device_malloc<T>(dof* num_all_elements ) : vec =GPUDevice::device_malloc<T>(dof* m_num_local_elements);
        else if(vtype == vec_type::host)
            (!is_local) ? vec = new T[dof* num_all_elements] : vec = new T[dof* m_num_local_elements];
        else if(vtype == vec_type::host_pinned)
            (!is_local) ? vec = GPUDevice::host_malloc<T>(dof* num_all_elements) : vec =GPUDevice::host_malloc<T>(dof* m_num_local_elements);

        return vec;
            
    }
    
    template<typename T, vec_type vtype>
    void MeshGPU::destroyVec(T* vec_ptr)
    {
        if(vtype == vec_type::device)
            GPUDevice::device_free<T>(vec_ptr);
        else if(vtype == vec_type::host)
            delete vec_ptr;
        else if(vtype == vec_type::host_pinned)
            GPUDevice::device_free<T>(vec_ptr);

        return;
    }


    template<typename T,typename stream>
    void MeshGPU::unzip_dg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out, unsigned int dof, stream s)
    {
        if(!pMesh->isActive())
            return;

        assert(pMesh->getNumLocalMeshElements() == pMesh->getLocalBlockList().size());
        dim3 gb  = dim3(m_unzip_grid[0],m_unzip_grid[1],dof);
        dim3 gb1 = dim3(m_num_local_elements,1,dof);
        dim3 tb  = dim3(8,8,1);

        //printf("Grid : {%d, %d, %d} blocks. Blocks : {%d, %d, %d} threads.\n", gb.x, gb.y, gb.z, tb.x, tb.y, tb.z);
        //__unzip_dg__3<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh, dptr_in,dptr_out);
        //__block_internal_unzip_dg__2<T,6,3> <<<gb1,tb,0,s>>> (dptr_mesh, dptr_in,dptr_out);
        __unzip_dg__2<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh, dptr_in,dptr_out);
        //__unzip_dg__1<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh, dptr_in,dptr_out);
        return;

    }

    template<typename T,typename stream>
    void MeshGPU::unzip_cg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out, unsigned int dof, stream s)
    {
        if(!pMesh->isActive())
            return;

        assert(pMesh->getNumLocalMeshElements() == pMesh->getLocalBlockList().size());
        //dim3 gb  = dim3(m_unzip_grid[0],m_unzip_grid[1],dof);
        dim3 gb  = dim3(m_unzip_grid[0],1,dof);
        //dim3 gb1 = dim3(m_num_local_elements,1,dof);
        dim3 tb  = dim3(7,7,1);
        //printf("Grid : {%d, %d, %d} blocks. Blocks : {%d, %d, %d} threads.\n", gb.x, gb.y, gb.z, tb.x, tb.y, tb.z);
        //__unzip_cg__2<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh, dptr_in,dptr_out);
        __unzip_cg1__2<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh, dptr_in,dptr_out);

    }

    template<typename T, typename stream>
    void MeshGPU::zip_dg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out,unsigned int dof, stream s)
    {
        dim3 gb  = dim3(m_zip_grid[0],m_zip_grid[1],dof);
        dim3 tb  = dim3(8,8,1);
        //printf("Grid : {%d, %d, %d} blocks. Blocks : {%d, %d, %d} threads.\n", gb.x, gb.y, gb.z, tb.x, tb.y, tb.z);
        __zip_dg__<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh,dptr_in,dptr_out);
        // the __zip_dg_enforce_c0__ does not work. The hanging node remove duplicates is complicated. 
        //__zip_dg_enforce_c0__<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh,dptr_in,dptr_out);
        return;

    }

    template<typename T, typename stream>
    void MeshGPU::zip_cg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out,unsigned int dof, stream s)
    {
        dim3 gb  = dim3(m_zip_grid[0],m_zip_grid[1],dof);
        dim3 tb  = dim3(8,8,1);
        //printf("Grid : {%d, %d, %d} blocks. Blocks : {%d, %d, %d} threads.\n", gb.x, gb.y, gb.z, tb.x, tb.y, tb.z);
        __zip_cg__<T,6,3> <<<gb,tb,0,s>>> (dptr_mesh,dptr_in,dptr_out);
        return;

    }


    template<typename T, typename stream>
    void MeshGPU::read_from_ghost_dg_begin(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device,
                                     const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh,
                                     T* dptr_vec, unsigned int dof, stream s)
    {   
        if(pMesh->getMPICommSizeGlobal()==1 || (!pMesh->isActive()))
            return;

        const std::vector<unsigned int >& se_c  = pMesh->getElementSendCounts();
        const std::vector<unsigned int >& se_f  = pMesh->getElementSendOffsets();

        const std::vector<unsigned int >& re_c  = pMesh->getElementRecvCounts();
        const std::vector<unsigned int >& re_f  = pMesh->getElementRecvOffsets();

        const std::vector<unsigned int>& send_plist = pMesh->getSendEleProcList();
        const std::vector<unsigned int>& recv_plist = pMesh->getRecvEleProcList();
        

        const const int   active_p = pMesh->getMPICommSize();
        const DEVICE_UINT nPe      = pMesh->getNumNodesPerElement();
        const unsigned int sendBSz = (se_f[active_p-1] + se_c[active_p-1]) * nPe;
        const unsigned int recvBSz = (re_f[active_p-1] + re_c[active_p-1]) * nPe;

        MPI_Comm commActive = pMesh->getMPICommunicator();

        
        dim3 gb  = dim3(pMesh->getMPICommSize(),m_max_send_ele,dof);
        dim3 tb  = dim3(8,8,1);

        T* s_buff = (T*)ctx_device.getSendBuffer();
        
        //printf("Grid : {%d, %d, %d} blocks. Blocks : {%d, %d, %d} threads.\n", gb.x, gb.y, gb.z, tb.x, tb.y, tb.z);
        __write_sb_dg__<T> <<<gb,tb,0,s>>>(dptr_mesh,dptr_vec,s_buff);

        #ifdef MPIX_CUDA_AWARE_SUPPORT
            T* const sendB = s_buff;
            if(recvBSz)
            {
                T* recvB=(T*)ctx_device.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p = 0 ;  recv_p <  recv_plist.size(); recv_p++)
                {
                    const unsigned int proc_id = recv_plist[recv_p];
                    par::Mpi_Irecv((recvB+ dof*nPe*re_f[proc_id]), dof * nPe*re_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_recv_req[recv_p]);
                }

            }
            GPUDevice::stream_synchronize(s);
            
            if(sendBSz)
            {
                for(unsigned int send_p = 0; send_p < send_plist.size(); send_p++) 
                {
                    const unsigned int proc_id=send_plist[send_p];
                    par::Mpi_Isend(sendB+dof*nPe*se_f[proc_id],dof*nPe*se_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_send_req[send_p]);
                }
            }
        #else
            T* const sendB=(T*)ctx_host.getSendBuffer();
            GPUDevice::device_to_host_async<DEVICE_REAL,cudaStream_t>(sendB, s_buff, dof*sendBSz,s);
            if(recvBSz)
            {
                T* recvB=(T*)ctx_host.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p = 0 ;  recv_p <  recv_plist.size(); recv_p++)
                {
                    const unsigned int proc_id = recv_plist[recv_p];
                    par::Mpi_Irecv((recvB+ dof*nPe*re_f[proc_id]), dof * nPe*re_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_recv_req[recv_p]);
                }

            }
            GPUDevice::stream_synchronize(s);
            if(sendBSz)
            {
                for(unsigned int send_p = 0; send_p < send_plist.size(); send_p++) 
                {
                    const unsigned int proc_id=send_plist[send_p];
                    par::Mpi_Isend(sendB+dof*nPe*se_f[proc_id],dof*nPe*se_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_send_req[send_p]);
                }
            }
        #endif

       

        m_ctag++;

        return ;
    }

    template<typename T, typename stream>
    void MeshGPU::read_from_ghost_dg_end(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device,
                                         const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh,
                                         T* dptr_vec, unsigned int dof, stream s)
    {
        if(pMesh->getMPICommSizeGlobal()==1 || (!pMesh->isActive()))
            return;

        const std::vector<unsigned int >& se_c  = pMesh->getElementSendCounts();
        const std::vector<unsigned int >& se_f  = pMesh->getElementSendOffsets();

        const std::vector<unsigned int >& re_c  = pMesh->getElementRecvCounts();
        const std::vector<unsigned int >& re_f  = pMesh->getElementRecvOffsets();

        const std::vector<unsigned int>& send_plist = pMesh->getSendEleProcList();
        const std::vector<unsigned int>& recv_plist = pMesh->getRecvEleProcList();

        const const int   active_p = pMesh->getMPICommSize();
        const DEVICE_UINT nPe      = pMesh->getNumNodesPerElement();
        const unsigned int sendBSz = (se_f[active_p-1] + se_c[active_p-1]) * nPe;
        const unsigned int recvBSz = (re_f[active_p-1] + re_c[active_p-1]) * nPe;

        if(recvBSz)
        {
            #ifdef MPIX_CUDA_AWARE_SUPPORT
                T* r_buff = (T*)ctx_device.getRecvBuffer();
            
                MPI_Waitall(send_plist.size(), ctx_host.m_send_req.data(), MPI_STATUSES_IGNORE);
                MPI_Waitall(recv_plist.size(), ctx_host.m_recv_req.data(), MPI_STATUSES_IGNORE);
                dim3 gb  = dim3(pMesh->getMPICommSize(),m_max_recv_ele,dof);
                dim3 tb  = dim3(8,8,1);
                __read_rb_dg__<T> <<<gb,tb,0,s>>>(dptr_mesh,r_buff,dptr_vec);
                GPUDevice::stream_synchronize(s);
            #else
                T* recvB  = (T*)ctx_host.getRecvBuffer();
                T* r_buff = (T*)ctx_device.getRecvBuffer();
                
                MPI_Waitall(send_plist.size(), ctx_host.m_send_req.data(), MPI_STATUSES_IGNORE);
                MPI_Waitall(recv_plist.size(), ctx_host.m_recv_req.data(), MPI_STATUSES_IGNORE);
                
            
                GPUDevice::host_to_device_async(recvB, r_buff, recvBSz * dof, s);
                dim3 gb  = dim3(pMesh->getMPICommSize(),m_max_recv_ele,dof);
                dim3 tb  = dim3(8,8,1);
                __read_rb_dg__<T> <<<gb,tb,0,s>>>(dptr_mesh,r_buff,dptr_vec);
                GPUDevice::stream_synchronize(s);
            #endif

        }

        return;

    }

    template<typename T, typename stream>
    void MeshGPU::read_from_ghost_cg_begin(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device,
                                            const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh,
                                            T* dptr_vec,unsigned int dof, stream s)
    {   
        
        if(pMesh->getMPICommSizeGlobal()==1 || (!pMesh->isActive()))
            return;

        const std::vector<unsigned int >& se_c  = pMesh->getNodalSendCounts();
        const std::vector<unsigned int >& se_f  = pMesh->getNodalSendOffsets();

        const std::vector<unsigned int >& re_c  = pMesh->getNodalRecvCounts();
        const std::vector<unsigned int >& re_f  = pMesh->getNodalRecvOffsets();

        const std::vector<unsigned int>& send_plist = pMesh->getSendProcList();
        const std::vector<unsigned int>& recv_plist = pMesh->getRecvProcList();
        

        const unsigned int active_p = pMesh->getMPICommSize();
        const unsigned int nPe      = pMesh->getNumNodesPerElement();
        const unsigned int sendBSz  = (se_f[active_p-1] + se_c[active_p-1]);
        const unsigned int recvBSz  = (re_f[active_p-1] + re_c[active_p-1]);

        
        dim3 gb  = dim3(m_max_send_cg/GPU_MAX_THREADS_PER_BLOCK + 1, pMesh->getMPICommSize(), dof);
        dim3 tb  = dim3(GPU_MAX_THREADS_PER_BLOCK,1,1);

        T* s_buff = (T*)ctx_device.getSendBuffer();
        //printf("Grid : {%d, %d, %d} blocks. Blocks : {%d, %d, %d} threads.\n", gb.x, gb.y, gb.z, tb.x, tb.y, tb.z);
        __write_sb_cg__<T> <<<gb,tb,0,s>>>(dptr_mesh,dptr_vec,s_buff);
        
        #ifdef MPIX_CUDA_AWARE_SUPPORT
            T* const sendB=s_buff;
            MPI_Comm commActive = pMesh->getMPICommunicator();
            if(recvBSz)
            {
                T* recvB=(T*)ctx_device.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p = 0 ;  recv_p <  recv_plist.size(); recv_p++)
                {
                    const unsigned int proc_id = recv_plist[recv_p];
                    par::Mpi_Irecv((recvB+ dof * re_f[proc_id]), dof * re_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_recv_req[recv_p]);
                }

            }
            GPUDevice::stream_synchronize(s);
            if(sendBSz)
            {
                for(unsigned int send_p = 0; send_p < send_plist.size(); send_p++) 
                {
                    const unsigned int proc_id=send_plist[send_p];
                    par::Mpi_Isend(sendB + dof *se_f[proc_id], dof * se_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_send_req[send_p]);
                }
            }
        #else
            T* const sendB=(T*)ctx_host.getSendBuffer();
            GPUDevice::device_to_host_async<DEVICE_REAL,cudaStream_t>(sendB, s_buff, dof*sendBSz,s);
            MPI_Comm commActive = pMesh->getMPICommunicator();

            if(recvBSz)
            {
                T* recvB=(T*)ctx_host.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p = 0 ;  recv_p <  recv_plist.size(); recv_p++)
                {
                    const unsigned int proc_id = recv_plist[recv_p];
                    par::Mpi_Irecv((recvB+ dof * re_f[proc_id]), dof * re_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_recv_req[recv_p]);
                }

            }
            GPUDevice::stream_synchronize(s);
            if(sendBSz)
            {
                for(unsigned int send_p = 0; send_p < send_plist.size(); send_p++) 
                {
                    const unsigned int proc_id=send_plist[send_p];
                    par::Mpi_Isend(sendB + dof *se_f[proc_id], dof * se_c[proc_id],proc_id,m_ctag,commActive,&ctx_host.m_send_req[send_p]);
                }
            }
        #endif

        m_ctag++;

        return;
    }

    template<typename T, typename stream>
    void MeshGPU::read_from_ghost_cg_end(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device,
                                            const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh,
                                            T* dptr_vec,unsigned int dof, stream s)
    {
        if(pMesh->getMPICommSizeGlobal()==1 || (!pMesh->isActive()))
            return;

        const std::vector<unsigned int >& se_c  = pMesh->getNodalSendCounts();
        const std::vector<unsigned int >& se_f  = pMesh->getNodalSendOffsets();

        const std::vector<unsigned int >& re_c  = pMesh->getNodalRecvCounts();
        const std::vector<unsigned int >& re_f  = pMesh->getNodalRecvOffsets();

        const std::vector<unsigned int>& send_plist = pMesh->getSendProcList();
        const std::vector<unsigned int>& recv_plist = pMesh->getRecvProcList();

        const const int   active_p = pMesh->getMPICommSize();
        const DEVICE_UINT nPe      = pMesh->getNumNodesPerElement();
        const unsigned int sendBSz = (se_f[active_p-1] + se_c[active_p-1]);
        const unsigned int recvBSz = (re_f[active_p-1] + re_c[active_p-1]);

        if(recvBSz)
        {

            #ifdef MPIX_CUDA_AWARE_SUPPORT
                T* r_buff = (T*)ctx_device.getRecvBuffer();
                MPI_Waitall(send_plist.size(), ctx_host.m_send_req.data(), MPI_STATUSES_IGNORE);
                MPI_Waitall(recv_plist.size(), ctx_host.m_recv_req.data(), MPI_STATUSES_IGNORE);
                dim3 gb  = dim3(m_max_recv_cg/GPU_MAX_THREADS_PER_BLOCK + 1, pMesh->getMPICommSize(),dof);
                dim3 tb  = dim3(GPU_MAX_THREADS_PER_BLOCK,1,1);
                __read_rb_cg__<T> <<<gb,tb,0,s>>>(dptr_mesh,r_buff,dptr_vec);
                GPUDevice::stream_synchronize(s);
            #else
                T* recvB  = (T*)ctx_host.getRecvBuffer();
                T* r_buff = (T*)ctx_device.getRecvBuffer();
                
                MPI_Waitall(send_plist.size(), ctx_host.m_send_req.data(), MPI_STATUSES_IGNORE);
                MPI_Waitall(recv_plist.size(), ctx_host.m_recv_req.data(), MPI_STATUSES_IGNORE);
                
                GPUDevice::host_to_device_async(recvB, r_buff, recvBSz * dof, s);
                dim3 gb  = dim3(m_max_recv_cg/GPU_MAX_THREADS_PER_BLOCK + 1, pMesh->getMPICommSize(),dof);
                dim3 tb  = dim3(GPU_MAX_THREADS_PER_BLOCK,1,1);
                __read_rb_cg__<T> <<<gb,tb,0,s>>>(dptr_mesh,r_buff,dptr_vec);
                GPUDevice::stream_synchronize(s);
            #endif

            
            
        }

        return;

    }

    
};