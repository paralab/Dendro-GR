/**
 * @file device_utils.cuh
 * @brief GPU device util functions. 
 * @version 0.1
 * @date 2022-02-26
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include "device.h"
#include "mesh.h"

namespace device
{

    template<typename T>
    void alloc_mpi_ctx(const ot::Mesh*pMesh, std::vector<ot::AsyncExchangeContex>& ctx_list, int dof, int async_k)
    {

        ctx_list.resize(async_k);

        if(pMesh->getMPICommSizeGlobal() ==1 || !pMesh->isActive())
            return;
        
        
        {
            const std::vector<unsigned int>& nodeSendCount  = pMesh->getNodalSendCounts();
            const std::vector<unsigned int>& nodeSendOffset = pMesh->getNodalSendOffsets();
                
            const std::vector<unsigned int>& e_sf  = pMesh->getElementSendOffsets();
            const std::vector<unsigned int>& e_sc  = pMesh->getElementSendCounts();
            const std::vector<unsigned int>& e_rf  = pMesh->getElementRecvOffsets();
            const std::vector<unsigned int>& e_rc  = pMesh->getElementRecvCounts();

            const std::vector<unsigned int>& nodeRecvCount  = pMesh->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeRecvOffset = pMesh->getNodalRecvOffsets();

            const std::vector<unsigned int>& sendProcList   = pMesh->getSendProcList();
            const std::vector<unsigned int>& recvProcList   = pMesh->getRecvProcList();

            const std::vector<unsigned int>& sendNodeSM     = pMesh->getSendNodeSM();
            const std::vector<unsigned int>& recvNodeSM     = pMesh->getRecvNodeSM();

            const unsigned int activeNpes=pMesh->getMPICommSize();
            const unsigned int nPe       =pMesh->getNumNodesPerElement(); 

            const unsigned int sendBSzCg = (nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1]);
            const unsigned int recvBSzCg = (nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1]);

            const unsigned int sendBSzDg = (e_sf[activeNpes-1] + e_sc[activeNpes-1]) * nPe;
            const unsigned int recvBSzDg = (e_rf[activeNpes-1] + e_rc[activeNpes-1]) * nPe;

            const unsigned int sendBSz = std::max(sendBSzCg,sendBSzDg);
            const unsigned int recvBSz = std::max(recvBSzCg,recvBSzDg);
            
            for(unsigned int i=0 ; i < async_k; i++)
            {
                const unsigned int v_begin  = ((i*dof)/async_k);
                const unsigned int v_end    = (((i+1)*dof)/async_k);
                const unsigned int batch_sz = (v_end-v_begin);

                if(sendBSz)
                    ctx_list[i].getSendBuffer(GPUDevice::device_malloc<T>(batch_sz * sendBSz));

                if(recvBSz) 
                    ctx_list[i].getRecvBuffer(GPUDevice::device_malloc<T>(batch_sz * recvBSz));

                ctx_list[i].m_send_req.resize(pMesh->getMPICommSize(), MPI_Request());
                ctx_list[i].m_recv_req.resize(pMesh->getMPICommSize(), MPI_Request());
            }
        
        }

        return;
        
    }

    template<typename T>
    void dealloc_mpi_ctx(const ot::Mesh*pMesh, std::vector<ot::AsyncExchangeContex>& ctx_list, int dof, int async_k)
    {

        if(pMesh->getMPICommSizeGlobal() ==1 || !pMesh->isActive())
            return;
        
        const std::vector<unsigned int>& nodeSendCount  = pMesh->getNodalSendCounts();
        const std::vector<unsigned int>& nodeSendOffset = pMesh->getNodalSendOffsets();
            
        const std::vector<unsigned int>& e_sf  = pMesh->getElementSendOffsets();
        const std::vector<unsigned int>& e_sc  = pMesh->getElementSendCounts();
        const std::vector<unsigned int>& e_rf  = pMesh->getElementRecvOffsets();
        const std::vector<unsigned int>& e_rc  = pMesh->getElementRecvCounts();

        const std::vector<unsigned int>& nodeRecvCount  = pMesh->getNodalRecvCounts();
        const std::vector<unsigned int>& nodeRecvOffset = pMesh->getNodalRecvOffsets();

        const std::vector<unsigned int>& sendProcList   = pMesh->getSendProcList();
        const std::vector<unsigned int>& recvProcList   = pMesh->getRecvProcList();

        const std::vector<unsigned int>& sendNodeSM     = pMesh->getSendNodeSM();
        const std::vector<unsigned int>& recvNodeSM     = pMesh->getRecvNodeSM();

        const unsigned int activeNpes=pMesh->getMPICommSize();
        const unsigned int nPe       =pMesh->getNumNodesPerElement(); 

        const unsigned int sendBSzCg = (nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1]);
        const unsigned int recvBSzCg = (nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1]);

        const unsigned int sendBSzDg = (e_sf[activeNpes-1] + e_sc[activeNpes-1]) * nPe;
        const unsigned int recvBSzDg = (e_rf[activeNpes-1] + e_rc[activeNpes-1]) * nPe;

        const unsigned int sendBSz = std::max(sendBSzCg,sendBSzDg);
        const unsigned int recvBSz = std::max(recvBSzCg,recvBSzDg);

        for(unsigned int i=0 ; i < ctx_list.size() ; i++)
        {
            const unsigned int v_begin  = ((i*dof)/async_k);
            const unsigned int v_end    = (((i+1)*dof)/async_k);
            const unsigned int batch_sz = (v_end-v_begin);

            if(sendBSz)
                GPUDevice::device_free<T>((T*)ctx_list[i].getSendBuffer());

            if(recvBSz)
                GPUDevice::device_free<T>((T*)ctx_list[i].getRecvBuffer());

            ctx_list[i].m_send_req.clear();
            ctx_list[i].m_recv_req.clear();
        }

        ctx_list.clear();
        return;
    }


}
    

    
