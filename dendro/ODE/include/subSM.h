/**
 * @file sm.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Exchange data by time stepper level. (scatter map by level. )
 * @version 0.1
 * @date 2020-01-16
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include <vector>
#include "mpi.h"
#include <bitset>
#include "TreeNode.h"
#include "mesh.h"
#include "asyncCtxByLev.h"
#include "rawIO.h"

namespace ot
{

    class SubScatterMap
    {
        
        protected:
            
            /**@brief: pointer to mesh */
            const ot::Mesh* m_uiMesh;

            /**@brief: min level of the grid*/
            unsigned int m_uiLMin;

            /**@brief: max level of the grid*/
            unsigned int m_uiLMax;

            /**@brief: pointer to send node tags by level*/
            const std::bitset<ot::TreeNode::OctantFlagType::MAX_LEVEL> * m_uiSendTag;
            
            /**@brief: pointer to recv node tags by level*/
            const std::bitset<ot::TreeNode::OctantFlagType::MAX_LEVEL> * m_uiRecvTag;

            /**@brief : send counts by level*/
            unsigned int ** m_uiSendCount = NULL;

            /**@brief: send offsets by level*/
            unsigned int ** m_uiSendOffset = NULL ;

            /**@brief : recv counts by level*/
            unsigned int ** m_uiRecvCount = NULL;
            
            /**@brief : recv offsets by level*/
            unsigned int ** m_uiRecvOffset = NULL;

            /**@brief: level to send proc list*/
            std::vector<unsigned int>* m_uiL2SendProcList = NULL;

            /**@brief: level to recv proc list*/
            std::vector<unsigned int>* m_uiL2RecvProcList = NULL;

            /**@brief: level to send scatter map*/
            unsigned int ** m_uiL2SendSM = NULL;

            /**@brief: level to recv scatter map*/
            unsigned int ** m_uiL2RecvSM = NULL;

            /**@brief: indicayes whether internal maps are allocated (if set to true)*/
            bool m_uiIsAllocated = false;

            /**@brief: List of async requests by lev*/
            std::vector<ot::AsyncCommByLev> m_uiAsyncCtxList;

            /**@brief: Async communication Tag*/
            unsigned int m_uiAsyncTag=0;


        private:

            inline void print_error(const char* const m) { 
                std::cout<<"[subSM Error] in "<<__FILE__<<" function: "<<__func__<<" at line "<<__LINE__<<" message "<<m<<std::endl;
            }



        protected:
            /**@brief: computes the Scatter maps for all levels. */
            void compute_L2SM();

        public: 
            /**
             * @brief Construct a new Sub Scatter Map object
             * 
             * @param pMesh pointer to the mesh. 
             * @param stag  send node tags
             * @param rtag  recv node tags
             */
            SubScatterMap(const ot::Mesh* pMesh, std::bitset<ot::TreeNode::OctantFlagType::MAX_LEVEL>* stag, std::bitset<ot::TreeNode::OctantFlagType::MAX_LEVEL>* rtag);

            /**@brief: default deallocator for the SM*/
            ~SubScatterMap();

            /**
             * @brief Perform ghost read begin (async) for the givel level 
             * @tparam T type of the vector
             * @param vec : input vector
             * @param l : level of the ghost sync
             * @param dof : number of dof
             */
            template<typename T>
            void readFromGhostBegin(T * vec, unsigned int l ,unsigned int dof=1);

            /**
             * @brief Perform ghost read end (async) for the givel level 
             * @tparam T type of the vector
             * @param vec : input vector
             * @param l : level of the ghost sync
             * @param dof : number of dof
             */
            template<typename T>
            void readFromGhostEnd(T * vec, unsigned int l ,unsigned int dof=1);

            
            
    };


    template<typename T>
    void SubScatterMap::readFromGhostBegin(T* vec, unsigned int l, unsigned int dof)
    {

        const bool isActive = m_uiMesh->isActive();

        if( m_uiMesh -> getMPICommSizeGlobal() == 1 || (!isActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiMesh->isActive())
        {

            MPI_Comm commActive = m_uiMesh->getMPICommunicator();
            
            const unsigned int activeRank = m_uiMesh->getMPIRank();
            const unsigned int activeNpes = m_uiMesh->getMPICommSize();
            const unsigned int total_sz_zip = m_uiMesh->getDegOfFreedom();

            const unsigned int sendBSz  =  m_uiSendOffset[l][activeNpes-1] +  m_uiSendCount[l][activeNpes-1];
            const unsigned int recvBSz  =  m_uiRecvOffset[l][activeNpes-1] +  m_uiRecvCount[l][activeNpes-1];

            ot::AsyncCommByLev ctx(vec,l);
            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz*dof));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p=0;recv_p<m_uiL2RecvProcList[l].size();recv_p++)
                {
                    const unsigned int proc_id = m_uiL2RecvProcList[l][recv_p];
                    
                    MPI_Request* req = new MPI_Request();
                    par::Mpi_Irecv((recvB+dof*m_uiRecvOffset[l][proc_id]),dof*m_uiRecvCount[l][proc_id],proc_id,m_uiAsyncTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*dof*sendBSz);
                sendB=(T*)ctx.getSendBuffer();

                for(unsigned int send_p=0; send_p<m_uiL2SendProcList[l].size(); send_p++) {

                    const unsigned int proc_id = m_uiL2SendProcList[l][send_p];
                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiSendOffset[l][proc_id]; k < (m_uiSendOffset[l][proc_id] + m_uiSendCount[l][proc_id]); k++)
                            sendB[dof*(m_uiSendOffset[l][proc_id]) + (var*m_uiSendCount[l][proc_id])+(k-m_uiSendOffset[l][proc_id])] = (vec + var * total_sz_zip )[m_uiL2SendSM[l][k]];
                    }

                }

                // active send procs
                for(unsigned int send_p=0; send_p<m_uiL2SendProcList[l].size(); send_p++)
                {
                    const unsigned int proc_id=m_uiL2SendProcList[l][send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+dof*m_uiSendOffset[l][proc_id], dof*m_uiSendCount[l][proc_id],proc_id,m_uiAsyncTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }
            }
            m_uiAsyncTag++;
            m_uiAsyncCtxList.push_back(ctx);

        }
        return;

    }


    template<typename T>
    void SubScatterMap::readFromGhostEnd(T* vec, unsigned int l, unsigned int dof)
    {
        const bool isActive = m_uiMesh->isActive();

        if( m_uiMesh -> getMPICommSizeGlobal() == 1 || (!isActive) )
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiMesh->isActive())
        {
            const unsigned int activeRank = m_uiMesh->getMPIRank();   
            const unsigned int activeNpes = m_uiMesh->getMPICommSize();
            const unsigned int total_sz_zip = m_uiMesh->getDegOfFreedom();

            const unsigned int sendBSz  =  m_uiSendOffset[l][activeNpes-1] +  m_uiSendCount[l][activeNpes-1];
            const unsigned int recvBSz  =  m_uiRecvOffset[l][activeNpes-1] +  m_uiRecvCount[l][activeNpes-1];

           
            int ctxIndex=-1;
            for(unsigned int i=0;i<m_uiAsyncCtxList.size();i++)
            {
                if(m_uiAsyncCtxList[i].getBuffer()==vec && m_uiAsyncCtxList[i].getLevel() == l)
                {
                    ctxIndex=i;
                    break;
                }

            }

            if(ctxIndex<0)
            {
                std::cout<<" Rank: "<<m_uiMesh->getMPIRank()<<" [SubSM Error]: "<<__LINE__<<" readFromGhostEnd async context not foound:  "<<std::endl;
                MPI_Abort(m_uiMesh->getMPICommunicator(),0);
            }

            ot::AsyncCommByLev ctx = m_uiAsyncCtxList[ctxIndex];

            MPI_Status status;
            //std::cout<<" ctx.getRequestList().size():  "<<ctx.getRequestList().size()<<std::endl;
            // need to wait for the commns to finish ...
            for (unsigned int i = 0; i < ctx.getRequestList().size(); i++) {
                MPI_Wait(ctx.getRequestList()[i], &status);
            }

            
            if(recvBSz)
            {
                // copy the recv data to the vec
                recvB=(T*)ctx.getRecvBuffer();

                for(unsigned int recv_p=0; recv_p < m_uiL2RecvProcList[l].size(); recv_p++ ){

                    const unsigned int proc_id = m_uiL2RecvProcList[l][recv_p];
                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiRecvOffset[l][proc_id]; k < (m_uiRecvOffset[l][proc_id] + m_uiRecvCount[l][proc_id]); k++)
                        {
                            (vec+var*total_sz_zip)[m_uiL2RecvSM[l][k]]=recvB[dof*(m_uiRecvOffset[l][proc_id]) + (var*m_uiRecvCount[l][proc_id])+(k-m_uiRecvOffset[l][proc_id])];
                        }
                            
                    }
                }
            }

            ctx.deAllocateSendBuffer();
            ctx.deAllocateRecvBuffer();

            for (unsigned int i = 0; i < ctx.getRequestList().size(); i++)
                delete ctx.getRequestList()[i];

            ctx.getRequestList().clear();

            // remove the context ...
            m_uiAsyncCtxList.erase(m_uiAsyncCtxList.begin() + ctxIndex);
            

        }

        return;
    }


} // end of namespace ot







