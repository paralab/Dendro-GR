/**
 * @file mesh_gpu.h
 * @author Milinda Fernando (milinda@oden.utexas.edu)
 * @brief minimal mesh class, for GPUs
 * @version 0.1
 * @date 2022-01-21
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include "device.h"
#include "block_gpu.h"
#include "mesh.h"
namespace device
{
    /**@brief where the vectors allocated. */
    enum vec_type{device=0, host, host_pinned};

    /**@brief keep track of the hanging configurations in gpu. */
    struct HangingInfo
    {
        public:
            /**@brief: DIR bit is one is that DIR is hanging.*/
            unsigned int hflag=0u;

            /**@brief: cnum if hanging. */
            unsigned int fcnum=0u;
            unsigned int ecnum=0u;

        /**@get face child number for interpolation if hangging*/
        #ifdef __CUDACC__
        __host__ __device__ 
        #endif
        unsigned int get_fcnum(unsigned int dir) const {return ((fcnum & (3u << (2*dir) )) >> ( (2*dir) )); }
        
       #ifdef __CUDACC__
        __host__ __device__ 
        #endif
        bool is_hanging(unsigned int dir) const { return (hflag & (1u<<dir)); }
        
        /**@get edge child number for interpolation if hangging*/
        #ifdef __CUDACC__
        __host__ __device__ 
        #endif
        unsigned int get_ecnum(unsigned int dir) const {return ((ecnum & ( 1u << (dir-NUM_FACES) ) ) >> ( (dir-NUM_FACES) ) );}

    };

    class MeshGPU
    {
        public:

            /**@brief: number of pre ghost nodes*/
            DEVICE_UINT m_num_pre_ghost_nodes=0;
            
            /**@brief: number of local (owned) nodes*/
            DEVICE_UINT m_num_local_nodes=0;

            /**@brief: number of post ghost nodes*/
            DEVICE_UINT m_num_post_ghost_nodes=0;

            /**@brief: number of pre ghost elements*/
            DEVICE_UINT m_num_pre_ghost_elements=0;

            /**@brief: number of local elements*/
            DEVICE_UINT m_num_local_elements=0;

            /**@brief: number of post elements*/
            DEVICE_UINT m_num_post_ghost_elements=0;

            /**@brief: number DG type nodes, octant local nodes, */
            DEVICE_UINT m_oct_local_sz=0;

            /**@brief: number CG type nodes, octant shared nodes */
            DEVICE_UINT m_oct_shared_sz=0;

            /**@brief: number unzip nodes */
            DEVICE_UINT m_oct_unzip_sz=0;

            /**@brief: element order*/
            DEVICE_UINT m_ele_order=0;

            /**@brief: number of mesh blocks*/
            DEVICE_UINT m_num_blocks=0;

            /**@brief: maxdepth */
            DEVICE_UINT m_max_depth=0;

            DEVICE_UINT m_unzip_grid[3]={1,1,1};
            
            DEVICE_UINT m_zip_grid[3]={1,1,1};

            /**@brief: element to block unzip map*/
            DEVICE_UINT* m_e2b_unzip_map        = nullptr;

            /**@brief: element to block unzip map counts*/
            DEVICE_UINT* m_e2b_unzip_map_count  = nullptr;

            /**@brief: element to block unzip map offset*/
            DEVICE_UINT* m_e2b_unzip_map_offset = nullptr;

            DEVICE_UINT* m_e2e                  = nullptr;
            DEVICE_UINT* m_e2n_cg               = nullptr;
            DEVICE_UINT* m_e2n_dg               = nullptr;

            /**@brief: all elements including ghost*/
            TNode3d*     m_all_elements         = nullptr; 

            /**@brief: block list*/
            BlockGPU3D * m_blk_list             = nullptr;

            DEVICE_UINT * m_elem_send_sm         = nullptr;
            DEVICE_UINT * m_elem_recv_sm         = nullptr;

            DEVICE_UINT * m_elem_send_c          = nullptr;
            DEVICE_UINT * m_elem_send_f          = nullptr;

            DEVICE_UINT * m_elem_recv_c          = nullptr;
            DEVICE_UINT * m_elem_recv_f          = nullptr;


            DEVICE_UINT * m_cg_send_sm           = nullptr;
            DEVICE_UINT * m_cg_recv_sm           = nullptr;

            DEVICE_UINT * m_cg_send_c            = nullptr;
            DEVICE_UINT * m_cg_send_f            = nullptr;

            DEVICE_UINT * m_cg_recv_c            = nullptr;
            DEVICE_UINT * m_cg_recv_f            = nullptr;
            HangingInfo* m_hinfo                 =nullptr;

            DEVICE_UINT m_max_send_ele           = 1;
            DEVICE_UINT m_max_recv_ele           = 1;

            DEVICE_UINT m_max_send_cg            = 1;
            DEVICE_UINT m_max_recv_cg            = 1;

            DEVICE_UINT m_ctag                   = 0;

            

        public:
            /**
             * @brief Construct a new Mesh GPU object
             * 
             */
            MeshGPU();
            /**
             * @brief Destroy the Mesh GPU object
             * 
             */
            ~MeshGPU();
            
            /**
             * @brief allocates GPU mesh data structures and send them to the device. 
             * @param pMesh CPU mesh object
             * @return MeshGPU* 
             */
            MeshGPU* alloc_mesh_on_device(const ot::Mesh* pMesh);
            
            
            /**
             * @brief free GPU mesh data structres allocated on the device
             * @param dptr_mesh 
             */
            void dealloc_mesh_on_device(MeshGPU* dptr_mesh);
            
            template<typename T, vec_type vtype>
            T* createVector(DEVICE_UINT dof=1, bool is_local=false);
            
            template<typename T, vec_type vtype>
            T* createUnZippedVector(DEVICE_UINT dof=1, bool is_local=false);
            
            template<typename T, vec_type vtype>
            T* createDGVector(DEVICE_UINT dof=1, bool is_local=false);

            template<typename T, vec_type vtype>
            T* createElementVector(DEVICE_UINT dof=1, bool is_local=false);
            
            template<typename T, vec_type vtype>
            void destroyVec(T* vec_ptr);

            template<typename T, typename stream>
            void unzip_dg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out, unsigned int dof=1, stream s = (stream)0);

            template<typename T, typename stream>
            void zip_dg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out, unsigned int dof=1, stream s = (stream)0);

            template<typename T, typename stream>
            void unzip_cg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out, unsigned int dof=1, stream s = (stream)0);

            template<typename T, typename stream>
            void zip_cg(const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, const T* const dptr_in, T* const dptr_out, unsigned int dof=1, stream s = (stream)0);

            template<typename T, typename stream>
            void read_from_ghost_dg_begin(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device, const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, T* dptr_vec, unsigned int dof=1, stream s = (stream)0);

            template<typename T, typename stream>
            void read_from_ghost_dg_end(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device, const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, T* dptr_vec, unsigned int dof=1, stream s = (stream)0);


            template<typename T, typename stream>
            void read_from_ghost_cg_begin(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device, const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, T* dptr_vec, unsigned int dof=1, stream s = (stream)0);

            template<typename T, typename stream>
            void read_from_ghost_cg_end(ot::AsyncExchangeContex& ctx_host,ot::AsyncExchangeContex& ctx_device, const ot::Mesh* const pMesh, const MeshGPU* const dptr_mesh, T* dptr_vec, unsigned int dof=1, stream s = (stream)0);

    };


    
}// end of namespace device. 