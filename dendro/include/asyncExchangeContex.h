/**
 * @file asyncExchangeContex.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Context information for asynchronous communication 
 * @version 0.1
 * @date 2018-11-19
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#ifndef DENDRO_5_0_UPDATECTX_H
#define DENDRO_5_0_UPDATECTX_H

#include "mpi.h"
#include <vector>

namespace ot {

    class AsyncExchangeContex {

        protected :
            /** pointer to the variable which perform the ghost exchange */
            void* m_uiBuffer=NULL;

            /** pointer to the send buffer*/
            void* m_uiSendBuf=NULL;

            /** pointer to the send buffer*/
            void* m_uiRecvBuf=NULL;

            std::vector<MPI_Request*>  m_uiRequests;

        public:
            /**@brief batched requests for send */
            std::vector<MPI_Request> m_send_req;
            
            /**@brief batched requests for recv */
            std::vector<MPI_Request> m_recv_req;

            AsyncExchangeContex(){};
            
            /**@brief creates an async ghost exchange contex*/
            AsyncExchangeContex(const void* var)
            {
                m_uiBuffer=(void*)var;
                m_uiSendBuf=nullptr;
                m_uiRecvBuf=nullptr;
                m_uiRequests.clear();
                m_send_req.clear();
                m_recv_req.clear();
            }

            /**@brief : defaut destructor*/
            ~AsyncExchangeContex() {};

            /**@brief allocates send buffer for ghost exchange*/
            inline void allocateSendBuffer(size_t bytes)
            {
                m_uiSendBuf=malloc(bytes);
            }

            /**@brief allocates recv buffer for ghost exchange*/
            inline void allocateRecvBuffer(size_t bytes)
            {
                m_uiRecvBuf=malloc(bytes);
            }

            /**@brief allocates send buffer for ghost exchange*/
            inline void deAllocateSendBuffer()
            {
                free(m_uiSendBuf);
                m_uiSendBuf=NULL;
            }

            /**@brief allocates recv buffer for ghost exchange*/
            inline void deAllocateRecvBuffer()
            {
                free(m_uiRecvBuf);
                m_uiRecvBuf=NULL;
            }

            inline void* getSendBuffer() { return m_uiSendBuf;}
            inline void* getRecvBuffer() { return m_uiRecvBuf;}

            inline void getSendBuffer(void* ptr) { m_uiSendBuf=ptr;}
            inline void getRecvBuffer(void* ptr) { m_uiRecvBuf=ptr;}

            inline const void* getBuffer() {return m_uiBuffer;}

            inline std::vector<MPI_Request*>& getRequestList(){ return m_uiRequests;}

            bool operator== (AsyncExchangeContex other) const{
                return( m_uiBuffer == other.m_uiBuffer );
            }
               

    };

} //end namespace

#endif //DENDRO_5_0_UPDATECTX_H
