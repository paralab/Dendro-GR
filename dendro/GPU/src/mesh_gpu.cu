/**
 * @file mesh_gpu.cu
 * @author Milinda Fernando
 * @brief  minimal mesh class for GPUs. 
 * @version 0.1
 * @date 2022-02-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "mesh_gpu.cuh"
namespace device
{
    MeshGPU::MeshGPU()
    {
        m_num_pre_ghost_nodes=0;
        m_num_local_nodes=0;
        m_num_post_ghost_nodes=0;


        m_num_pre_ghost_elements=0;
        m_num_local_elements=0;
        m_num_post_ghost_elements=0;

        m_num_blocks=0;


        m_oct_local_sz  = 0;
        m_oct_shared_sz = 0;
        m_oct_unzip_sz  = 0;
        m_ele_order     = 0;
        m_max_depth     = 0; 
        m_ctag=0;
        m_max_send_ele  = 1;
        m_max_send_ele  = 1;
        
        m_e2b_unzip_map        = nullptr;
        m_e2b_unzip_map_count  = nullptr;
        m_e2b_unzip_map_offset = nullptr;
        m_all_elements         = nullptr;
        m_blk_list             = nullptr;

        m_elem_send_sm         = nullptr;
        m_elem_recv_sm         = nullptr;
        m_elem_send_c          = nullptr;
        m_elem_send_f          = nullptr;
        m_elem_recv_c          = nullptr;
        m_elem_recv_f          = nullptr;

        m_cg_send_sm           = nullptr;
        m_cg_recv_sm           = nullptr;
        m_cg_send_c            = nullptr;
        m_cg_send_f            = nullptr;
        m_cg_recv_c            = nullptr;
        m_cg_recv_f            = nullptr;

        m_e2e                  = nullptr;
        m_e2n_cg               = nullptr;    
        m_e2n_dg               = nullptr;

    }

    MeshGPU::~MeshGPU(){}
    
    MeshGPU* MeshGPU::alloc_mesh_on_device(const ot::Mesh* pMesh)
    {

        MeshGPU* dptr_mesh=nullptr;
        if(!pMesh->isActive())
            return dptr_mesh;

        const std::vector<unsigned int >& m2b   = pMesh->getE2BUnzipMap();
        const std::vector<unsigned int >& m2b_c = pMesh->getE2BUnzipMapCounts();
        const std::vector<unsigned int >& m2b_f = pMesh->getE2BUnzipMapOffsets();
        
        const std::vector<unsigned int >& se_sm = pMesh->getSendElementSM();
        const std::vector<unsigned int >& re_sm = pMesh->getRecvElementSM();

        const std::vector<unsigned int >& se_c  = pMesh->getElementSendCounts();
        const std::vector<unsigned int >& se_f  = pMesh->getElementSendOffsets();

        const std::vector<unsigned int >& re_c  = pMesh->getElementRecvCounts();
        const std::vector<unsigned int >& re_f  = pMesh->getElementRecvOffsets();

        const std::vector<unsigned int >& s_cg_sm = pMesh->getSendNodeSM();
        const std::vector<unsigned int >& r_cg_sm = pMesh->getRecvNodeSM();

        const std::vector<unsigned int >& s_cg_c  = pMesh->getNodalSendCounts();
        const std::vector<unsigned int >& s_cg_f  = pMesh->getNodalSendOffsets();

        const std::vector<unsigned int >& r_cg_c  = pMesh->getNodalRecvCounts();
        const std::vector<unsigned int >& r_cg_f  = pMesh->getNodalRecvOffsets();

        const std::vector<unsigned int> & e2e   = pMesh->getE2EMapping();
        const std::vector<unsigned int> & en_cg = pMesh->getE2NMapping();
        const std::vector<unsigned int> & en_dg = pMesh->getE2NMapping_DG();

        const std::vector<ot::Block>& blk_list  = pMesh->getLocalBlockList();
        const std::vector<ot::TreeNode>& pNodes = pMesh->getAllElements();

        const unsigned int nPe          = pMesh->getNumNodesPerElement();
        const unsigned int ep           = pMesh->getElementOrder(); 

        this->m_num_pre_ghost_elements  = pMesh->getNumPreGhostElements();
        this->m_num_local_elements      = pMesh->getNumLocalMeshElements();
        this->m_num_post_ghost_elements = pMesh->getNumPostGhostElements();

        this->m_num_pre_ghost_nodes     = pMesh->getNumPreMeshNodes();
        this->m_num_local_nodes         = pMesh->getNumLocalMeshNodes();
        this->m_num_post_ghost_nodes    = pMesh->getNumPostMeshNodes();
        this->m_num_blocks              = pMesh->getLocalBlockList().size();
        
        this->m_oct_local_sz            = pMesh->getDegOfFreedomDG();
        this->m_oct_shared_sz           = pMesh->getDegOfFreedom();
        this->m_oct_unzip_sz            = pMesh->getDegOfFreedomUnZip();

        this->m_ele_order               = pMesh->getElementOrder();
        this->m_max_depth               = m_uiMaxDepth;
        this->m_ctag=0;

        this->m_unzip_grid[0]           = pNodes.size();
        this->m_unzip_grid[1]           = 1u;
        if(m2b_c.size()>0)
            this->m_unzip_grid[1] = std::max(1u,*std::max_element(m2b_c.begin(),m2b_c.end()));

        this->m_max_send_ele            = 1u;
        if(se_c.size()>0)
            this->m_max_send_ele = std::max(1u,*std::max_element(se_c.begin(),se_c.end()));

        this->m_max_recv_ele            = 1u;
        if(re_c.size()>0)
            this->m_max_recv_ele = std::max(1u,*std::max_element(re_c.begin(),re_c.end()));

        this->m_max_send_cg             = 1u;
        if(s_cg_c.size()>0)
            this->m_max_send_cg  = std::max(1u,*std::max_element(s_cg_c.begin(),s_cg_c.end()));

        this->m_max_recv_cg            = 1u;
        if(r_cg_c.size()>0)
            this->m_max_recv_cg  = std::max(1u,*std::max_element(r_cg_c.begin(),r_cg_c.end()));

        
        this->m_zip_grid[2]             = 1;
        this->m_zip_grid[1]             = blk_list.size();
        this->m_zip_grid[0]             = 1;
        
        for(unsigned int i=0; i < blk_list.size(); i++)
            if((blk_list[i].getLocalElementEnd()-blk_list[i].getLocalElementBegin()) > m_zip_grid[0] )
                this->m_zip_grid[0] = (blk_list[i].getLocalElementEnd()-blk_list[i].getLocalElementBegin());
        
        this->m_e2b_unzip_map           = GPUDevice::device_malloc<DEVICE_UINT>(m2b.size());
        this->m_e2b_unzip_map_count     = GPUDevice::device_malloc<DEVICE_UINT>(m2b_c.size());
        this->m_e2b_unzip_map_offset    = GPUDevice::device_malloc<DEVICE_UINT>(m2b_f.size());
        this->m_blk_list                = GPUDevice::device_malloc<BlockGPU3D>(blk_list.size());
        this->m_all_elements            = GPUDevice::device_malloc<TNode3d>(pNodes.size());
        this->m_hinfo                   = GPUDevice::device_malloc<HangingInfo>(pNodes.size());
        

        this->m_e2e                     = GPUDevice::device_malloc<DEVICE_UINT>(e2e.size());
        this->m_e2n_cg                  = GPUDevice::device_malloc<DEVICE_UINT>(en_cg.size());
        this->m_e2n_dg                  = GPUDevice::device_malloc<DEVICE_UINT>(en_dg.size());
        

        this->m_elem_send_sm            = GPUDevice::device_malloc<DEVICE_UINT>(se_sm.size());
        this->m_elem_recv_sm            = GPUDevice::device_malloc<DEVICE_UINT>(re_sm.size());
        this->m_elem_send_c             = GPUDevice::device_malloc<DEVICE_UINT>(se_c.size());
        this->m_elem_send_f             = GPUDevice::device_malloc<DEVICE_UINT>(se_f.size());
        this->m_elem_recv_c             = GPUDevice::device_malloc<DEVICE_UINT>(re_c.size());
        this->m_elem_recv_f             = GPUDevice::device_malloc<DEVICE_UINT>(re_f.size());


        this->m_cg_send_sm            = GPUDevice::device_malloc<DEVICE_UINT>(s_cg_sm.size());
        this->m_cg_recv_sm            = GPUDevice::device_malloc<DEVICE_UINT>(r_cg_sm.size());
        this->m_cg_send_c             = GPUDevice::device_malloc<DEVICE_UINT>(s_cg_c.size());
        this->m_cg_send_f             = GPUDevice::device_malloc<DEVICE_UINT>(s_cg_f.size());
        this->m_cg_recv_c             = GPUDevice::device_malloc<DEVICE_UINT>(r_cg_c.size());
        this->m_cg_recv_f             = GPUDevice::device_malloc<DEVICE_UINT>(r_cg_f.size());


        

        dptr_mesh = GPUDevice::device_malloc<MeshGPU>(1);
        GPUDevice::host_to_device(this,dptr_mesh,1);
        
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(m2b.data()   , this->m_e2b_unzip_map,m2b.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(m2b_c.data() , this->m_e2b_unzip_map_count,m2b_c.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(m2b_f.data() , this->m_e2b_unzip_map_offset,m2b_f.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(e2e.data()   , this->m_e2e   , e2e.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(en_cg.data() , this->m_e2n_cg,en_cg.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(en_dg.data() , this->m_e2n_dg,en_dg.size());


        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(se_sm.data(), this->m_elem_send_sm,se_sm.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(re_sm.data(), this->m_elem_recv_sm,re_sm.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(se_c.data() , this->m_elem_send_c,se_c.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(se_f.data() , this->m_elem_send_f,se_f.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(re_c.data() , this->m_elem_recv_c,re_c.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(re_f.data() , this->m_elem_recv_f,re_f.size());


        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(s_cg_sm.data(), this->m_cg_send_sm,s_cg_sm.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(r_cg_sm.data(), this->m_cg_recv_sm,r_cg_sm.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(s_cg_c.data() , this->m_cg_send_c ,s_cg_c.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(s_cg_f.data() , this->m_cg_send_f ,s_cg_f.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(r_cg_c.data() , this->m_cg_recv_c ,r_cg_c.size());
        GPUDevice::host_to_device_async<DEVICE_UINT,cudaStream_t>(r_cg_f.data() , this->m_cg_recv_f ,r_cg_f.size());

        
        BlockGPU3D*  blk_buffer     = GPUDevice::host_malloc<BlockGPU3D>(blk_list.size());
        TNode3d*     tn_buffer      = GPUDevice::host_malloc<TNode3d>(pNodes.size());
        HangingInfo* hanging_info   = GPUDevice::host_malloc<HangingInfo>(pNodes.size());

        // set hanging info. 
        for(unsigned int ele=0; ele < pNodes.size(); ele++)
        {
            for(unsigned int dir=0; dir < NUM_FACES; dir++)
            {
                unsigned int cnum=0;
                if(pMesh->isFaceHanging(ele,dir, cnum))
                {
                    hanging_info[ele].hflag |=(1u<<dir);
                    hanging_info[ele].fcnum |=(cnum<<(2*dir));
                    // if(hanging_info[ele].is_hanging(dir)!=true && hanging_info[ele].get_fcnum(dir)!=cnum)
                    //     printf("dire %d  hflag %d fcnum %d \n",dir, hanging_info[ele].hflag, hanging_info[ele].fcnum);
                }
            }


            for(unsigned int dir=NUM_FACES; dir < NUM_FACES + NUM_EDGES; dir++)
            {
                unsigned int cnum=0;
                if(pMesh->isEdgeHanging(ele,dir,cnum))
                {
                    hanging_info[ele].hflag |=(1u<<(dir));
                    hanging_info[ele].ecnum |=(cnum<<((dir-NUM_FACES)));
                    // if(hanging_info[ele].is_hanging(dir)!=true && hanging_info[ele].get_ecnum(dir)!=cnum)
                    //     printf("dire %d  hflag %d ecnum %d \n",dir, hanging_info[ele].hflag, hanging_info[ele].ecnum);
                }

                
            }
        }

        // for(unsigned int ele=0; ele < pNodes.size(); ele++)
        // {
        //     for(unsigned int dir=NUM_FACES; dir < NUM_FACES + NUM_EDGES; dir++)
        //     {
        //         unsigned int cnum=0;
        //         if(pMesh->isEdgeHanging(ele,dir,cnum))
        //         {
        //             // hanging_info[ele].hflag |=(1u<<(dir));
        //             // hanging_info[ele].ecnum |=(cnum<<((dir-NUM_FACES)));
        //             if(hanging_info[ele].is_hanging(dir)!=true && hanging_info[ele].get_ecnum(dir)!=cnum)
        //                 printf("dire %d  hflag %d ecnum %d \n",dir, hanging_info[ele].hflag, hanging_info[ele].ecnum);
        //         }else
        //         {
        //             if(hanging_info[ele].is_hanging(dir)==true )
        //                 printf("dire %d  hflag %d ecnum %d \n",dir, hanging_info[ele].hflag, hanging_info[ele].ecnum);
        //         }
        //     }
                
        // }



        GPUDevice::host_to_device_async<HangingInfo,cudaStream_t>(hanging_info, this->m_hinfo ,pNodes.size());
        
        const Point pt_min = pMesh->getDomainMinPt();
        const Point pt_max = pMesh->getDomainMaxPt();

        for(unsigned int i=0; i< blk_list.size(); i++)
        {
            blk_buffer[i] = BlockGPU3D(blk_list[i]);
            
            blk_buffer[i].m_dx[0] = blk_list[i].computeDx(pt_min,pt_max);
            blk_buffer[i].m_dx[1] = blk_list[i].computeDy(pt_min,pt_max);
            blk_buffer[i].m_dx[2] = blk_list[i].computeDz(pt_min,pt_max);

            const unsigned int PW = blk_list[i].get1DPadWidth(); 
            
            const Point oct_coord = Point((double)blk_list[i].getBlockNode().minX(),(double)blk_list[i].getBlockNode().minY(),(double)blk_list[i].getBlockNode().minZ());
            
            Point domain_pt;    
            pMesh->octCoordToDomainCoord(oct_coord,domain_pt);

            blk_buffer[i].m_ptMin[0] = domain_pt.x() - PW * blk_buffer[i].m_dx[0];
            blk_buffer[i].m_ptMin[1] = domain_pt.y() - PW * blk_buffer[i].m_dx[1];
            blk_buffer[i].m_ptMin[2] = domain_pt.z() - PW * blk_buffer[i].m_dx[2];
            
        }
        
        GPUDevice::host_to_device_async<BlockGPU3D,cudaStream_t>(blk_buffer,this->m_blk_list,blk_list.size());
            
        for(unsigned int i=0; i< pNodes.size(); i++)
            tn_buffer[i] =  TNode3d(pNodes[i]);

        GPUDevice::host_to_device_async<TNode3d,cudaStream_t>(tn_buffer,this->m_all_elements,pNodes.size());

        const RefElement * refel = pMesh->getReferenceElement();
        cudaMemcpyToSymbolAsync(refel_1d, (const DEVICE_REAL*)refel->getIMTChild0(),sizeof(DEVICE_REAL)*(ep+1)*(ep+1),0);
        cudaMemcpyToSymbolAsync(refel_1d, (const DEVICE_REAL*)refel->getIMTChild1(),sizeof(DEVICE_REAL)*(ep+1)*(ep+1),sizeof(DEVICE_REAL)*(ep+1)*(ep+1));
        

        GPUDevice::device_synchronize();

        GPUDevice::host_free<BlockGPU3D>(blk_buffer);
        GPUDevice::host_free<TNode3d>(tn_buffer);
        
        return dptr_mesh;
        

    }

    void MeshGPU::dealloc_mesh_on_device(MeshGPU* dptr_mesh)
    {
        if(dptr_mesh==nullptr)
            return;

        GPUDevice::device_free<DEVICE_UINT>(this->m_e2b_unzip_map);
        GPUDevice::device_free<DEVICE_UINT>(this->m_e2b_unzip_map_count);
        GPUDevice::device_free<DEVICE_UINT>(this->m_e2b_unzip_map_offset);
        
        GPUDevice::device_free<DEVICE_UINT>(this->m_e2e);
        GPUDevice::device_free<DEVICE_UINT>(this->m_e2n_cg);
        GPUDevice::device_free<DEVICE_UINT>(this->m_e2n_dg);

        GPUDevice::device_free<DEVICE_UINT>(this->m_elem_send_sm);
        GPUDevice::device_free<DEVICE_UINT>(this->m_elem_recv_sm);
        GPUDevice::device_free<DEVICE_UINT>(this->m_elem_send_c);
        GPUDevice::device_free<DEVICE_UINT>(this->m_elem_send_f);
        GPUDevice::device_free<DEVICE_UINT>(this->m_elem_recv_c);
        GPUDevice::device_free<DEVICE_UINT>(this->m_elem_recv_f);

        GPUDevice::device_free<DEVICE_UINT>(this->m_cg_send_sm);
        GPUDevice::device_free<DEVICE_UINT>(this->m_cg_recv_sm);
        GPUDevice::device_free<DEVICE_UINT>(this->m_cg_send_c);
        GPUDevice::device_free<DEVICE_UINT>(this->m_cg_send_f);
        GPUDevice::device_free<DEVICE_UINT>(this->m_cg_recv_c);
        GPUDevice::device_free<DEVICE_UINT>(this->m_cg_recv_f);

        GPUDevice::device_free<BlockGPU3D>(this->m_blk_list);
        GPUDevice::device_free<TNode3d>(this->m_all_elements);
        GPUDevice::device_free<HangingInfo>(this->m_hinfo);
        
        this->m_e2b_unzip_map           = nullptr;
        this->m_e2b_unzip_map_count     = nullptr;
        this->m_e2b_unzip_map_offset    = nullptr;
        
        this->m_elem_send_sm            = nullptr;
        this->m_elem_recv_sm            = nullptr;
        this->m_elem_send_c             = nullptr;
        this->m_elem_send_f             = nullptr;
        this->m_elem_recv_c             = nullptr;
        this->m_elem_recv_f             = nullptr;

        this->m_cg_send_sm              = nullptr;
        this->m_cg_recv_sm              = nullptr;
        this->m_cg_send_c               = nullptr;
        this->m_cg_send_f               = nullptr;
        this->m_cg_recv_c               = nullptr;
        this->m_cg_recv_f               = nullptr;


        this->m_blk_list                = nullptr;
        this->m_all_elements            = nullptr;

        GPUDevice::device_free(dptr_mesh);
        return;
    }

    
}