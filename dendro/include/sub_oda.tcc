/**
 * @author Milinda Fernando
 * school of computing university of Utah. 
 * @brief Contains template functions for the sub_oda.tcc
 * */

namespace ot
{

    /**@brief init function for DA_FLAG ALL*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::ALL>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementPreGhostBegin;
        m_uiLoopInfo.indexEnd=m_uiElementPostGhostEnd;
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }


    /**@brief init function for DA_FLAG WRITABLE*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::WRITABLE>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementPreGhostBegin;
        m_uiLoopInfo.indexEnd=m_uiElementPostGhostEnd;

        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 

        unsigned int count=m_uiLoopInfo.indexBegin;
        while((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_INDEPENDENT_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_INDEPENDENT_FLAG_BIT)!=1)) count++;
        m_uiLoopInfo.currentIndex=count;

    }


    /**@brief init function for DA_FLAG INDEPENDENT*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::INDEPENDENT>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementLocalBegin;
        m_uiLoopInfo.indexEnd=m_uiElementLocalEnd;

        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 

        unsigned int count=m_uiLoopInfo.indexBegin;
        while((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_INDEPENDENT_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_INDEPENDENT_FLAG_BIT)!=1)) count++;

        m_uiLoopInfo.currentIndex=count;

    }

    /**@brief init function for DA_FLAG W_DEPENDENT*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::W_DEPENDENT>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementPreGhostBegin;
        m_uiLoopInfo.indexEnd=m_uiElementPostGhostEnd;

        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 
        unsigned int count=m_uiLoopInfo.indexBegin;
        
        while((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_W_DEPENDENT_FLAG_BIT)!=1)) count++;

        m_uiLoopInfo.currentIndex=count;
    }

    /**@brief init function for DA_FLAG W_BOUNDARY*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::W_BOUNDARY>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementLocalBegin;
        m_uiLoopInfo.indexEnd=m_uiElementLocalEnd;

        unsigned int count=m_uiLoopInfo.indexBegin;
        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 

        while((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_W_BOUNDARY_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[count]],ODA_W_BOUNDARY_FLAG_BIT)!=1)) count++;

        m_uiLoopInfo.currentIndex=count;
    }

    /**@brief init function for local elements (octants)*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::LOCAL_ELEMENTS>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementLocalBegin;
        m_uiLoopInfo.indexEnd=m_uiElementLocalEnd;
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }

    /**@brief init function for pre-ghost elements (octants)*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::PREGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementPreGhostBegin;
        m_uiLoopInfo.indexEnd=m_uiElementPreGhostEnd;
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }

    /**@brief init function for post-ghost elements (octants)*/
    template<>
    inline void subDA::init<ot::DA_FLAGS::POSTGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.indexBegin=m_uiElementPostGhostBegin;
        m_uiLoopInfo.indexEnd=m_uiElementPostGhostEnd;
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }

    /**@brief end function for DA_FLAG ALL*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::ALL>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG WRITABLE*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::WRITABLE>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG INDEPENDENT*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::INDEPENDENT>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG W_DEPENDENT*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::W_DEPENDENT>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG W_BOUNDARY*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::W_BOUNDARY>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG LOCAL_ELEMENTS*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::LOCAL_ELEMENTS>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG PREGHOST_ELEMENTS*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::PREGHOST_ELEMENTS>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG POSTGHOST_ELEMENTS*/
    template<>
    inline unsigned int subDA::end<ot::DA_FLAGS::POSTGHOST_ELEMENTS>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief next function for ALL*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::ALL>()
    {
        m_uiLoopInfo.currentIndex++;
    }

    /**@brief next function for ALL*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::WRITABLE>()
    {
        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 
        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_INDEPENDENT_FLAG_BIT)!=1) ) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_INDEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

    }

    /**@brief next function for INDEPENDENT*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::INDEPENDENT>()
    {
        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 
        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_INDEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_INDEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

    }

    /**@brief next function for W_DEPENDENT*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::W_DEPENDENT>()
    {
        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 
        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_W_DEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_W_DEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;


    }

    /**@brief next function for W_BOUNDARY*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::W_BOUNDARY>()
    {
        const unsigned int * octFlags = &(*(m_da->getOctFlags().begin())); 
        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_W_BOUNDARY_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(octFlags[m_uip_sub2DA_ElemMap[m_uiLoopInfo.currentIndex]],ODA_W_BOUNDARY_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;


    }

    /**@brief next function for LOCAL_ELEMENTS*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::LOCAL_ELEMENTS>()
    {
        m_uiLoopInfo.currentIndex++;
    }

    /**@brief next function for PREGHOST_ELEMENTS*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::PREGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.currentIndex++;
    }

    /**@brief next function for POSTGHOST_ELEMENTS*/
    template<>
    inline void subDA::next<ot::DA_FLAGS::POSTGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.currentIndex++;
    }


    /**@brief get the current elmenent in the iteration*/
    inline unsigned int subDA::curr()
    {
        return m_uiLoopInfo.currentIndex;
    }

    template <typename T>
    int subDA::createVector(T*& local, bool isElemental, bool isGhosted, unsigned int dof) const
    {

        if(!(m_uiIsActive))
        {
            local=NULL;

        }else {

            if(isElemental)
            {
                if(isGhosted)
                    local=new T[dof*m_uiTotalElementSz];
                else
                    local=new T[dof*m_uiLocalElementSz];

            }else {

                if(isGhosted)
                    local=new T[dof*m_uiTotalNodalSz];
                else
                    local=new T[dof*m_uiLocalNodalSz];
            }
        }

        return 0;


    }

    template<typename T>
    int subDA::createVector(std::vector<T>& local, bool isElemental, bool isGhosted, unsigned int dof) const
    {
        if(!(m_uiIsActive))
        {
            local.clear();

        }else {

            if(isElemental)
            {
                if(isGhosted)
                    local.resize(dof*m_uiTotalElementSz);
                else
                    local.resize(dof*m_uiLocalElementSz);

            }else {

                if(isGhosted)
                    local.resize(dof*m_uiTotalNodalSz);
                else
                    local.resize(dof*m_uiLocalNodalSz);
            }
        }

        return 0;
    }


    template <typename T>
    void subDA::destroyVector(T*& local) const
    {
        delete [] local;
        local=NULL;
    }

    template <typename T>
    void subDA::destroyVector(std::vector<T>& local) const
    {
        local.clear();
    }

    template<typename T>
    void subDA::nodalVecToGhostedNodal(const T* in, T*& out,bool isAllocated,unsigned int dof) const
    {

        if(!(m_uiIsActive))
            return;

        if(!isAllocated)
            createVector<T>(out,false,true,dof);

        for(unsigned int var=0;var<dof;var++)
        {
            std::memcpy((out+var*m_uiTotalNodalSz+m_uiNodeLocalBegin),(in+var*m_uiLocalNodalSz),sizeof(T)*(m_uiLocalNodalSz));
        }


    }

    template<typename T>
    void subDA::ghostedNodalToNodalVec(const T* gVec,T*& local,bool isAllocated,unsigned int dof) const
    {
        if(!(m_uiIsActive))
            return;

        if(!isAllocated)
            createVector(local,false,false,dof);

        for(unsigned int var=0;var<dof;var++)
            std::memcpy((local + var*m_uiLocalNodalSz ),(gVec+(var*m_uiTotalNodalSz)+m_uiNodeLocalBegin),sizeof(T)*(m_uiLocalNodalSz));

    }

    template <typename T>
    void subDA::readFromGhostBegin(T* vec,unsigned int dof)
    {

        if(m_da->getNpesAll()==1)
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiIsActive)
        {
            const unsigned int sendBSz=m_uiSendOffsets[m_uiNpesActive-1] + m_uiSendCounts[m_uiNpesActive-1];
            const unsigned int recvBSz=m_uiRecvOffsets[m_uiNpesActive-1] + m_uiRecvCounts[m_uiNpesActive-1];
            unsigned int proc_id;
            
            AsyncExchangeContex ctx(vec);
            
            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz*dof));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p=0;recv_p<m_uiRecvProcList.size();recv_p++)
                {
                    proc_id=m_uiRecvProcList[recv_p];
                    MPI_Request* req=new MPI_Request();
                    par::Mpi_Irecv((recvB+dof*m_uiRecvOffsets[proc_id]),dof*m_uiRecvCounts[proc_id],proc_id,m_uiCommTag,m_uiCommActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*dof*sendBSz);
                sendB=(T*)ctx.getSendBuffer();

                for(unsigned int send_p=0;send_p<m_uiSendProcList.size();send_p++) {
                    
                    proc_id=m_uiSendProcList[send_p];
                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiSendOffsets[proc_id]; k < (m_uiSendOffsets[proc_id] + m_uiSendCounts[proc_id]); k++)
                        {
                            sendB[dof*(m_uiSendOffsets[proc_id]) + (var*m_uiSendCounts[proc_id])+(k-m_uiSendOffsets[proc_id])] = (vec+var*m_uiTotalNodalSz)[m_uiSendScatterMap[k]];
                        }

                    }



                }

                // active send procs
                for(unsigned int send_p=0;send_p<m_uiSendProcList.size();send_p++)
                {
                    proc_id=m_uiSendProcList[send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+dof*m_uiSendOffsets[proc_id],dof*m_uiSendCounts[proc_id],proc_id,m_uiCommTag,m_uiCommActive,req);
                    ctx.getRequestList().push_back(req);

                }


            }

            m_uiCommTag++;
            m_uiMPIContexts.push_back(ctx);


        }

        return;

    }

    template <typename T>
    void subDA::readFromGhostEnd(T *vec,unsigned int dof)
    {
        if(m_da->getNpesAll()==1)
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiIsActive)
        {
          
            const unsigned int sendBSz=m_uiSendOffsets[m_uiNpesActive-1] + m_uiSendCounts[m_uiNpesActive-1];
            const unsigned int recvBSz=m_uiRecvOffsets[m_uiNpesActive-1] + m_uiRecvCounts[m_uiNpesActive-1];
          
            unsigned int proc_id;

            unsigned int ctxIndex=0;
            for(unsigned int i=0;i<m_uiMPIContexts.size();i++)
            {
                if(m_uiMPIContexts[i].getBuffer()==vec)
                {
                    ctxIndex=i;
                    break;
                }

            }


            MPI_Status status;
            // need to wait for the commns to finish ...
            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++) {
                MPI_Wait(m_uiMPIContexts[ctxIndex].getRequestList()[i], &status);
            }
            // todo: check !!!
            if(recvBSz && m_uiTotalNodalSz)
            {
                // copy the recv data to the vec
                recvB=(T*)m_uiMPIContexts[ctxIndex].getRecvBuffer();

                for(unsigned int recv_p=0;recv_p<m_uiRecvProcList.size();recv_p++){
                    proc_id=m_uiRecvProcList[recv_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiRecvOffsets[proc_id]; k < (m_uiRecvOffsets[proc_id] + m_uiRecvCounts[proc_id]); k++)
                        {
                            (vec+var*m_uiTotalNodalSz)[m_uiRecvScatterMap[k]]=recvB[dof*(m_uiRecvOffsets[proc_id]) + (var*m_uiRecvCounts[proc_id])+(k-m_uiRecvOffsets[proc_id])];
                        }
                    }

                }

            }



            m_uiMPIContexts[ctxIndex].deAllocateSendBuffer();
            m_uiMPIContexts[ctxIndex].deAllocateRecvBuffer();

            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++)
                delete m_uiMPIContexts[ctxIndex].getRequestList()[i];

            m_uiMPIContexts[ctxIndex].getRequestList().clear();

            // remove the context ...
            m_uiMPIContexts.erase(m_uiMPIContexts.begin() + ctxIndex);


        }

        return;


    }


    // template <typename T>
    // void subDA::getElementNodalValues(const T*in, T* eleVecOut,unsigned int eleID,unsigned int dof) const
    // {
    //     //assert((eleID>=m_uiMesh->getElementLocalBegin()) && (eleID<m_uiMesh->getElementLocalEnd()));
    //     const unsigned int nPe=m_da->getNumNodesPerElement();
    //     for(unsigned int var=0;var<dof;var++)
    //         m_da->getElementNodalValues();
    //         //m_uiMesh->getElementNodalValues(&in[var*m_uiTotalNodalSz],&eleVecOut[var*nPe],eleID);

    // }


    // template <typename T>
    // void subDA::eleVecToVecAccumilation(T*out, const T* eleVecIn,unsigned int eleID,unsigned int dof) const
    // {
    //     const unsigned int nPe=m_uiMesh->getNumNodesPerElement();

    //     for(unsigned int var=0;var<dof;var++)
    //         m_uiMesh->computeElementalContribution(&eleVecIn[var*nPe],&out[var*m_uiTotalNodalSz],eleID);

    // }

    // template <typename T>
    // void subDA::setVectorByFunction(T* local,std::function<void(T,T,T,T*)>func,bool isElemental, bool isGhosted, unsigned int dof) const
    // {

    //     ot::TreeNode tmpOct;
    //     const ot::TreeNode* allElements=&(*(m_uiMesh->getAllElements().begin()));
    //     const unsigned int * e2n=&(*(m_uiMesh->getE2NMapping().begin()));
    //     const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
    //     const unsigned int * cgToDg=&(*(m_uiMesh->getCG2DGMap().begin()));

    //     unsigned int offset;
    //     unsigned int arrSz;

    //     unsigned int lookup=0;
    //     unsigned int ownerID, ii_x, jj_y, kk_z;
    //     double hx;
    //     const unsigned int eleOrder=m_uiMesh->getElementOrder();

    //     T x,y,z;
    //     T* val=new T[dof];

    //     if(!isElemental)
    //     {
    //         if(isGhosted) {
    //             offset = 0;
    //             arrSz = m_uiTotalNodalSz;
    //         }else{
    //             offset=m_uiMesh->getNodeLocalBegin();
    //             arrSz=m_uiLocalNodalSz;
    //         }


    //         for(unsigned int node=m_uiMesh->getNodeLocalBegin();node<m_uiMesh->getNodeLocalEnd();node++)
    //         {
    //             lookup=cgToDg[node];
    //             m_uiMesh->dg2eijk(lookup,ownerID,ii_x,jj_y,kk_z);
    //             tmpOct=allElements[ownerID];
    //             hx=(tmpOct.maxX()-tmpOct.minX())/((double) eleOrder);
    //             x=tmpOct.minX() + ii_x*(hx);
    //             y=tmpOct.minY() + jj_y*(hx);
    //             z=tmpOct.minZ() + kk_z*(hx);

    //             func(x,y,z,val);

    //             for(unsigned int var=0;var<dof;var++)
    //                 local[ ( var * arrSz ) + (node - offset) ]=val[var];

    //         }

    //     }else {

    //         if(isGhosted) {
    //             offset = 0;
    //             arrSz = m_uiTotalElementSz;
    //         }else{
    //             offset=m_uiMesh->getElementLocalBegin();
    //             arrSz=m_uiLocalElementSz;
    //         }


    //         for(unsigned int ele=m_uiMesh->getElementLocalBegin();ele<m_uiMesh->getElementLocalEnd();ele++)
    //         {
    //             tmpOct=allElements[ele];
    //             x=tmpOct.minX();
    //             y=tmpOct.minY();
    //             z=tmpOct.minZ();

    //             func(x,y,z,val);
    //             for(unsigned int var=0;var<dof;var++)
    //                 local[ (var*arrSz) + (ele-offset) ] = val[var];

    //         }

    //     }



    //     delete [] val;




    // }


    // template <typename T>
    // void subDA::setVectorByScalar(T* local,const T* value,bool isElemental, bool isGhosted, unsigned int dof) const
    // {

    //     ot::TreeNode tmpOct;
    //     const ot::TreeNode* allElements=&(*(m_uiMesh->getAllElements().begin()));
    //     unsigned int arrSz;

    //     if(!isElemental)
    //     {

    //         if(isGhosted) {
    //             arrSz = m_uiTotalNodalSz;
    //         }else{
    //             arrSz=m_uiLocalNodalSz;
    //         }

    //         for(unsigned int var=0;var<dof;var++)
    //         {
    //             for(unsigned int node=0;node<arrSz;node++)
    //                 local[ (var*arrSz) + node]=value[var];
    //         }


    //     }else{

    //         if(isGhosted) {
    //             arrSz = m_uiTotalElementSz;
    //         }else{
    //             arrSz=m_uiLocalElementSz;
    //         }

    //         for(unsigned int var=0;var<dof;var++)
    //         {
    //             for(unsigned int ele=0;ele<arrSz;ele++)
    //                 local[ (var*arrSz) + ele]=value[var];
    //         }

    //     }


    // }

    template<typename T>
    T* subDA::getVecPointerToDof(T* in ,unsigned int dofInex, bool isElemental,bool isGhosted) const
    {

        if(!(m_uiIsActive))
            return NULL;

        unsigned int arrSz;

        if(!isElemental)
        {
            if(isGhosted) {
                arrSz = m_uiTotalNodalSz;
            }else{
                arrSz=m_uiLocalNodalSz;
            }


        }else{

            if(isGhosted) {
                arrSz = m_uiTotalElementSz;
            }else{
                arrSz=m_uiLocalElementSz;
            }

        }

        return (T*)(&in[dofInex*arrSz]);

    }


    
    template<typename T>
    void subDA::copyVector(T* dest,const T* source,bool isElemental,bool isGhosted) const
    {
        if(!(m_uiIsActive))
            return ;

        unsigned int arrSz;

        if(!isElemental)
        {
            if(isGhosted) {
                arrSz = m_uiTotalNodalSz;
            }else{
                arrSz=m_uiLocalNodalSz;
            }


        }else{

            if(isGhosted) {
                arrSz = m_uiTotalElementSz;
            }else{
                arrSz=m_uiLocalElementSz;
            }

        }


        std::memcpy(dest,source,sizeof(T)*arrSz);

    }

    template<typename T>
    void subDA::copyVectors(T* dest,const T* source,bool isElemental,bool isGhosted,unsigned int dof) const
    {
        if(!(m_uiIsActive))
            return ;

        unsigned int arrSz;

        if(!isElemental)
        {
            if(isGhosted) {
                arrSz = m_uiTotalNodalSz;
            }else{
                arrSz=m_uiLocalNodalSz;
            }


        }else{

            if(isGhosted) {
                arrSz = m_uiTotalElementSz;
            }else{
                arrSz=m_uiLocalElementSz;
            }

        }


        std::memcpy(dest,source,sizeof(T)*arrSz*dof);
    }


    // template<typename T>
    // int subDA::getFaceNeighborValues(unsigned int eleID, const T* in, T* out, T* coords, unsigned int * neighID, unsigned int face, NeighbourLevel & level) const
    // {
    //     if(!m_uiMesh->isActive())
    //         return false;

    //     int state=m_uiMesh->getFaceNeighborValues(eleID,in,out,coords,neighID,face,level);

    //     return state;

    // }


#ifdef BUILD_WITH_PETSC


    template<typename T>
    void subDA::petscSetVectorByFunction(Vec& local,std::function<void(T,T,T,T*)>func,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        PetscScalar * arry=NULL;
        VecGetArray(local,&arry);

        setVectorByFunction(arry,func,isElemental,isGhosted,dof);

        VecRestoreArray(local,&arry);


    }

    template <typename T>
    void subDA::petscSetVectorByScalar(Vec& local,const T* value,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        PetscScalar * arry=NULL;
        VecGetArray(local,&arry);

        setVectorByScalar(arry,value,isElemental,isGhosted,dof);

        VecRestoreArray(local,&arry);

    }

    

#endif




    
} // namespace ot


