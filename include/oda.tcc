/**
 * @author Milinda Fernando
 * @brief Contains template functions for the oda.h
 * */

namespace ot
{

    template<typename T>
    DA::DA(std::function<void(T,T,T,T*)>func,unsigned int dofSz,MPI_Comm comm,unsigned int order,double interp_tol, unsigned int grainSz,double sfc_tol,SM_TYPE smType)
    {


        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        std::vector<ot::TreeNode> tmpNodes;
        DendroIntL localSz,globalSz;
        double t_stat_g[3];
        double t_stat;

        double t1=MPI_Wtime();

        unsigned int numInterpSz=dofSz;
        unsigned int varIndex[numInterpSz];
        for(unsigned int i=0;i<numInterpSz;i++)
            varIndex[i]=i;

        function2Octree(func, dofSz,(unsigned int*)varIndex,numInterpSz,tmpNodes, m_uiMaxDepth, interp_tol, order, comm);
        double t2=MPI_Wtime();

        t_stat=t2-t1;
        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
        t_stat_g[1]=t_stat_g[1]/(double)npes;

        if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

        par::Mpi_Bcast(&globalSz,1,0,comm);


        bool isActive;
        MPI_Comm commActive;

        if((globalSz/grainSz)>=npes)
        {
            MPI_Comm_dup(comm,&commActive);
            isActive=true;

        }else
        {
            isActive=(rank*grainSz<globalSz);
            par::splitComm2way(isActive,&commActive,comm);

        }

        shrinkOrExpandOctree(tmpNodes,sfc_tol,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

        if(!isActive)
            if(tmpNodes.size()!=0)
                std::cout<<" rank_g: "<<rank<<" isActive: "<<isActive<<" f2O octants: "<<tmpNodes.size()<<std::endl;


        std::vector<ot::TreeNode> balOct;
        localSz=0;
        if(isActive)
        {

            int rank_active,npes_active;

            MPI_Comm_size(commActive,&npes_active);
            MPI_Comm_rank(commActive,&rank_active);

            if(!rank_active) std::cout<<"[MPI_COMM_SWITCH]: "<<npes_active<<std::endl;

            ot::TreeNode root(m_uiDim,m_uiMaxDepth);
            std::vector<ot::TreeNode> tmpVec;
            t1=MPI_Wtime();

            SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,sfc_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,DENDRO_DEFAULT_SF_K,commActive);
            std::swap(tmpNodes,tmpVec);
            tmpVec.clear();

            SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,sfc_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
            std::swap(tmpNodes,tmpVec);
            tmpVec.clear();

            t2=MPI_Wtime();
            t_stat=t2-t1;

            par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
            par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
            par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
            t_stat_g[1]=t_stat_g[1]/(double)npes_active;

            localSz=tmpNodes.size();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

            if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
            if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


            t1=MPI_Wtime();

            SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,sfc_tol,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
            tmpNodes.clear();

            t2=MPI_Wtime();

            t_stat=t2-t1;
            par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
            par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
            par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
            t_stat_g[1]=t_stat_g[1]/(double)npes_active;

            if(!rank_active) std::cout<<GRN<<" 2:1 balancing max (s): "<<t_stat_g[2]<<NRM<<std::endl;
            localSz=balOct.size();


        }
        MPI_Comm_free(&commActive);
        // all reduce act as barrier to sync all procs.
        par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
        if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;


        assert(par::test::isUniqueAndSorted(balOct,comm));
        // generates the mesh.
        m_uiMesh=new ot::Mesh(balOct,1,order,comm,false,smType,grainSz,sfc_tol);
        balOct.clear();

        #ifdef DEBUG_ODA
                std::vector<ot::TreeNode> independent;
                std::vector<ot::TreeNode> w_dependent;
                std::vector<ot::TreeNode> w_bdy;
        #endif

        // perform element flags.
        ot::computeODAFlags(m_uiMesh,m_uiOctantFlags);
        // compute the local node to global node mapping.
        ot::computeLocalToGlobalNodalMap(m_uiMesh,m_uiLocalToGlobalNodalMap,m_uiGlobalNodeSz,m_uiNodalOffset);



        #ifdef DEBUG_ODA
                for(init<DA_FLAGS::INDEPENDENT>();curr()<end<DA_FLAGS::INDEPENDENT>();next<DA_FLAGS::INDEPENDENT>())
                    independent.push_back(allElements[curr()]);

                for(init<DA_FLAGS::W_DEPENDENT>();curr()<end<DA_FLAGS::W_DEPENDENT>();next<DA_FLAGS::W_DEPENDENT>())
                    w_dependent.push_back(allElements[curr()]);

                for(init<DA_FLAGS::W_BOUNDARY>();curr()<end<DA_FLAGS::W_BOUNDARY>();next<DA_FLAGS::W_BOUNDARY>())
                    w_bdy.push_back(allElements[curr()]);//std::cout<<"loop ele: "<<curr()<<" : "<<allElements[curr()]<<std::endl;//


                io::vtk::oct2vtu(&(*(independent.begin())),independent.size(),"independent",comm);
                io::vtk::oct2vtu(&(*(w_dependent.begin())),w_dependent.size(),"w_dependent",comm);
                io::vtk::oct2vtu(&(*(w_bdy.begin())),w_bdy.size(),"w_bdy",comm);

        #endif


        m_uiTotalElementSz=m_uiMesh->getAllElements().size();
        m_uiLocalElementSz=m_uiMesh->getNumLocalMeshElements();

        m_uiTotalNodalSz=m_uiMesh->getDegOfFreedom();
        m_uiLocalNodalSz=m_uiMesh->getNumLocalMeshNodes();


        m_uiLoopInfo.currentIndex=0;
        m_uiLoopInfo.indexBegin=0;
        m_uiLoopInfo.indexEnd=0;

        m_uiOctreeLowerBound[0]=0;
        m_uiOctreeLowerBound[1]=0;
        m_uiOctreeLowerBound[2]=0;

        m_uiOctreeUpperBound[0]=1u<<(m_uiMaxDepth);
        m_uiOctreeUpperBound[1]=1u<<(m_uiMaxDepth);
        m_uiOctreeUpperBound[2]=1u<<(m_uiMaxDepth);


        m_uiCommTag=0;
        m_uiMPIContexts.clear();
    }


    /**@brief init function for DA_FLAG ALL*/
    template<>
    inline void DA::init<ot::DA_FLAGS::ALL>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementPreGhostBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementPostGhostEnd();
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }


    /**@brief init function for DA_FLAG WRITABLE*/
    template<>
    inline void DA::init<ot::DA_FLAGS::WRITABLE>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementPreGhostBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementPostGhostEnd();

        unsigned int count=m_uiLoopInfo.indexBegin;
        while((binOp::getBit(m_uiOctantFlags[count],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(m_uiOctantFlags[count],ODA_INDEPENDENT_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(m_uiOctantFlags[count],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(m_uiOctantFlags[count],ODA_INDEPENDENT_FLAG_BIT)!=1)) count++;
        m_uiLoopInfo.currentIndex=count;

    }


    /**@brief init function for DA_FLAG INDEPENDENT*/
    template<>
    inline void DA::init<ot::DA_FLAGS::INDEPENDENT>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementLocalBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementLocalEnd();

        unsigned int count=m_uiLoopInfo.indexBegin;
        while((binOp::getBit(m_uiOctantFlags[count],ODA_INDEPENDENT_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(m_uiOctantFlags[count],ODA_INDEPENDENT_FLAG_BIT)!=1)) count++;

        m_uiLoopInfo.currentIndex=count;

    }

    /**@brief init function for DA_FLAG W_DEPENDENT*/
    template<>
    inline void DA::init<ot::DA_FLAGS::W_DEPENDENT>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementPreGhostBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementPostGhostEnd();

        unsigned int count=m_uiLoopInfo.indexBegin;
        while((binOp::getBit(m_uiOctantFlags[count],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(m_uiOctantFlags[count],ODA_W_DEPENDENT_FLAG_BIT)!=1)) count++;

        m_uiLoopInfo.currentIndex=count;
    }

    /**@brief init function for DA_FLAG W_BOUNDARY*/
    template<>
    inline void DA::init<ot::DA_FLAGS::W_BOUNDARY>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementLocalBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementLocalEnd();

        unsigned int count=m_uiLoopInfo.indexBegin;
        while((binOp::getBit(m_uiOctantFlags[count],ODA_W_BOUNDARY_FLAG_BIT)!=1) && (count<m_uiLoopInfo.indexEnd)) count++;

        if((binOp::getBit(m_uiOctantFlags[count],ODA_W_BOUNDARY_FLAG_BIT)!=1)) count++;

        m_uiLoopInfo.currentIndex=count;
    }

    /**@brief init function for local elements (octants)*/
    template<>
    inline void DA::init<ot::DA_FLAGS::LOCAL_ELEMENTS>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementLocalBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementLocalEnd();
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }

    /**@brief init function for pre-ghost elements (octants)*/
    template<>
    inline void DA::init<ot::DA_FLAGS::PREGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementPreGhostBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementPreGhostEnd();
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }

    /**@brief init function for post-ghost elements (octants)*/
    template<>
    inline void DA::init<ot::DA_FLAGS::POSTGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.indexBegin=m_uiMesh->getElementPostGhostBegin();
        m_uiLoopInfo.indexEnd=m_uiMesh->getElementPostGhostEnd();
        m_uiLoopInfo.currentIndex=m_uiLoopInfo.indexBegin;

    }


    /**@brief end function for DA_FLAG ALL*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::ALL>()
    {
        return m_uiLoopInfo.indexEnd;
    }


    /**@brief end function for DA_FLAG WRITABLE*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::WRITABLE>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG INDEPENDENT*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::INDEPENDENT>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG W_DEPENDENT*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::W_DEPENDENT>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG W_BOUNDARY*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::W_BOUNDARY>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG LOCAL_ELEMENTS*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::LOCAL_ELEMENTS>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG PREGHOST_ELEMENTS*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::PREGHOST_ELEMENTS>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief end function for DA_FLAG POSTGHOST_ELEMENTS*/
    template<>
    inline unsigned int DA::end<ot::DA_FLAGS::POSTGHOST_ELEMENTS>()
    {
        return m_uiLoopInfo.indexEnd;
    }

    /**@brief next function for ALL*/
    template<>
    inline void DA::next<ot::DA_FLAGS::ALL>()
    {
        m_uiLoopInfo.currentIndex++;
    }

    /**@brief next function for ALL*/
    template<>
    inline void DA::next<ot::DA_FLAGS::WRITABLE>()
    {

        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_INDEPENDENT_FLAG_BIT)!=1) ) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_W_DEPENDENT_FLAG_BIT)!=1) && (binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_INDEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

    }

    /**@brief next function for INDEPENDENT*/
    template<>
    inline void DA::next<ot::DA_FLAGS::INDEPENDENT>()
    {

        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_INDEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_INDEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

    }

    /**@brief next function for W_DEPENDENT*/
    template<>
    inline void DA::next<ot::DA_FLAGS::W_DEPENDENT>()
    {

        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_W_DEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_W_DEPENDENT_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;


    }

    /**@brief next function for W_BOUNDARY*/
    template<>
    inline void DA::next<ot::DA_FLAGS::W_BOUNDARY>()
    {

        m_uiLoopInfo.currentIndex++;
        if(m_uiLoopInfo.currentIndex>=m_uiLoopInfo.indexEnd)
            return;

        while(((m_uiLoopInfo.currentIndex)<m_uiLoopInfo.indexEnd) && (binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_W_BOUNDARY_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;

        if((binOp::getBit(m_uiOctantFlags[m_uiLoopInfo.currentIndex],ODA_W_BOUNDARY_FLAG_BIT)!=1)) m_uiLoopInfo.currentIndex++;


    }

    /**@brief next function for LOCAL_ELEMENTS*/
    template<>
    inline void DA::next<ot::DA_FLAGS::LOCAL_ELEMENTS>()
    {
        m_uiLoopInfo.currentIndex++;
    }

    /**@brief next function for PREGHOST_ELEMENTS*/
    template<>
    inline void DA::next<ot::DA_FLAGS::PREGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.currentIndex++;
    }

    /**@brief next function for POSTGHOST_ELEMENTS*/
    template<>
    inline void DA::next<ot::DA_FLAGS::POSTGHOST_ELEMENTS>()
    {
        m_uiLoopInfo.currentIndex++;
    }


    /**@brief get the current elmenent in the iteration*/
    inline unsigned int DA::curr()
    {
        return m_uiLoopInfo.currentIndex;
    }

    template <typename T>
    int DA::createVector(T*& local, bool isElemental, bool isGhosted, unsigned int dof) const
    {

        if(!(m_uiMesh->isActive()))
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
    int DA::createVector(std::vector<T>& local, bool isElemental, bool isGhosted, unsigned int dof) const
    {
        if(!(m_uiMesh->isActive()))
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
    void DA::destroyVector(T*& local) const
    {
        delete [] local;
        local=NULL;
    }

    template <typename T>
    void DA::destroyVector(std::vector<T>& local) const
    {
        local.clear();
    }

    template<typename T>
    void DA::nodalVecToGhostedNodal(const T* in, T*& out,bool isAllocated,unsigned int dof) const
    {

        if(!(m_uiMesh->isActive()))
            return;

        if(!isAllocated)
            createVector<T>(out,false,true,dof);

        for(unsigned int var=0;var<dof;var++)
        {
            std::memcpy((out+var*m_uiTotalNodalSz+m_uiMesh->getNodeLocalBegin()),(in+var*m_uiLocalNodalSz),sizeof(T)*(m_uiLocalNodalSz));
        }


    }

    template<typename T>
    void DA::ghostedNodalToNodalVec(const T* gVec,T*& local,bool isAllocated,unsigned int dof) const
    {
        if(!(m_uiMesh->isActive()))
            return;

        if(!isAllocated)
            createVector(local,false,false,dof);

        for(unsigned int var=0;var<dof;var++)
            std::memcpy((local + var*m_uiLocalNodalSz ),(gVec+(var*m_uiTotalNodalSz)+m_uiMesh->getNodeLocalBegin()),sizeof(T)*(m_uiLocalNodalSz));

    }



    template <typename T>
    void DA::readFromGhostBegin(T* vec,unsigned int dof)
    {

        if(m_uiMesh->getMPICommSizeGlobal()==1)
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiMesh->isActive())
        {
            const std::vector<unsigned int>& nodeSendCount=m_uiMesh->getNodalSendCounts();
            const std::vector<unsigned int>& nodeSendOffset=m_uiMesh->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeRecvCount=m_uiMesh->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeRecvOffset=m_uiMesh->getNodalRecvOffsets();

            const std::vector<unsigned int>& sendProcList=m_uiMesh->getSendProcList();
            const std::vector<unsigned int>& recvProcList=m_uiMesh->getRecvProcList();

            const std::vector<unsigned int>& sendNodeSM=m_uiMesh->getSendNodeSM();
            const std::vector<unsigned int>& recvNodeSM=m_uiMesh->getRecvNodeSM();


            const unsigned int activeNpes=m_uiMesh->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            AsyncExchangeContex ctx(vec);
            MPI_Comm commActive=m_uiMesh->getMPICommunicator();


            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz*dof));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++)
                {
                    proc_id=recvProcList[recv_p];
                    MPI_Request* req=new MPI_Request();
                    par::Mpi_Irecv((recvB+dof*nodeRecvOffset[proc_id]),dof*nodeRecvCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*dof*sendBSz);
                sendB=(T*)ctx.getSendBuffer();

                for(unsigned int send_p=0;send_p<sendProcList.size();send_p++) {
                    proc_id=sendProcList[send_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = nodeSendOffset[proc_id]; k < (nodeSendOffset[proc_id] + nodeSendCount[proc_id]); k++)
                        {
                            sendB[dof*(nodeSendOffset[proc_id]) + (var*nodeSendCount[proc_id])+(k-nodeSendOffset[proc_id])] = (vec+var*m_uiTotalNodalSz)[sendNodeSM[k]];
                        }

                    }



                }

                // active send procs
                for(unsigned int send_p=0;send_p<sendProcList.size();send_p++)
                {
                    proc_id=sendProcList[send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+dof*nodeSendOffset[proc_id],dof*nodeSendCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }


            }

            m_uiCommTag++;
            m_uiMPIContexts.push_back(ctx);


        }

        return;

    }

    template <typename T>
    void DA::readFromGhostEnd(T *vec,unsigned int dof)
    {
        if(m_uiMesh->getMPICommSizeGlobal()==1)
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiMesh->isActive())
        {
            const std::vector<unsigned int>& nodeSendCount=m_uiMesh->getNodalSendCounts();
            const std::vector<unsigned int>& nodeSendOffset=m_uiMesh->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeRecvCount=m_uiMesh->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeRecvOffset=m_uiMesh->getNodalRecvOffsets();

            const std::vector<unsigned int>& sendProcList=m_uiMesh->getSendProcList();
            const std::vector<unsigned int>& recvProcList=m_uiMesh->getRecvProcList();

            const std::vector<unsigned int>& sendNodeSM=m_uiMesh->getSendNodeSM();
            const std::vector<unsigned int>& recvNodeSM=m_uiMesh->getRecvNodeSM();


            const unsigned int activeNpes=m_uiMesh->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
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

            if(recvBSz)
            {
                // copy the recv data to the vec
                recvB=(T*)m_uiMPIContexts[ctxIndex].getRecvBuffer();

                for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++){
                    proc_id=recvProcList[recv_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = nodeRecvOffset[proc_id]; k < (nodeRecvOffset[proc_id] + nodeRecvCount[proc_id]); k++)
                        {
                            (vec+var*m_uiTotalNodalSz)[recvNodeSM[k]]=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
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

    template<typename T>
    void DA::writeToGhostsBegin(T* vec, unsigned int dof) 
    {
        if(m_uiMesh->getMPICommSizeGlobal()==1)
            return;
        
        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiMesh->isActive())
        {
            const std::vector<unsigned int>& nodeRecvCount = m_uiMesh->getNodalSendCounts();
            const std::vector<unsigned int>& nodeRecvOffset = m_uiMesh->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeSendCount = m_uiMesh->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeSendOffset = m_uiMesh->getNodalRecvOffsets();

            const std::vector<unsigned int>& recvProcList=m_uiMesh->getSendProcList();
            const std::vector<unsigned int>& sendProcList=m_uiMesh->getRecvProcList();

            const std::vector<unsigned int>& recvNodeSM = m_uiMesh->getSendNodeSM();
            const std::vector<unsigned int>& sendNodeSM = m_uiMesh->getRecvNodeSM();


            const unsigned int activeNpes=m_uiMesh->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            AsyncExchangeContex ctx(vec);
            MPI_Comm commActive=m_uiMesh->getMPICommunicator();


            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz*dof));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++)
                {
                    proc_id=recvProcList[recv_p];
                    MPI_Request* req=new MPI_Request();
                    par::Mpi_Irecv((recvB+dof*nodeRecvOffset[proc_id]),dof*nodeRecvCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*dof*sendBSz);
                sendB=(T*)ctx.getSendBuffer();

                for(unsigned int send_p=0;send_p<sendProcList.size();send_p++) {
                    proc_id=sendProcList[send_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = nodeSendOffset[proc_id]; k < (nodeSendOffset[proc_id] + nodeSendCount[proc_id]); k++)
                        {
                            sendB[dof*(nodeSendOffset[proc_id]) + (var*nodeSendCount[proc_id])+(k-nodeSendOffset[proc_id])] = (vec+var*m_uiTotalNodalSz)[sendNodeSM[k]];
                        }

                    }



                }

                // active send procs
                for(unsigned int send_p=0;send_p<sendProcList.size();send_p++)
                {
                    proc_id=sendProcList[send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+dof*nodeSendOffset[proc_id],dof*nodeSendCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }


            }

            m_uiCommTag++;
            m_uiMPIContexts.push_back(ctx);


        }

        return;

    }

    template<typename T>
    void DA::writeToGhostsEnd(T* vec, DA_FLAGS::WriteMode mode, unsigned int dof)
    {
        if(m_uiMesh->getMPICommSizeGlobal()==1)
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;
        
        if(m_uiMesh->isActive())
        {
            const std::vector<unsigned int>& nodeRecvCount = m_uiMesh->getNodalSendCounts();
            const std::vector<unsigned int>& nodeRecvOffset = m_uiMesh->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeSendCount = m_uiMesh->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeSendOffset = m_uiMesh->getNodalRecvOffsets();

            const std::vector<unsigned int>& recvProcList=m_uiMesh->getSendProcList();
            const std::vector<unsigned int>& sendProcList=m_uiMesh->getRecvProcList();

            const std::vector<unsigned int>& recvNodeSM = m_uiMesh->getSendNodeSM();
            const std::vector<unsigned int>& sendNodeSM = m_uiMesh->getRecvNodeSM();


            const unsigned int activeNpes=m_uiMesh->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
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

            if(recvBSz)
            {
                // copy the recv data to the vec
                recvB=(T*)m_uiMPIContexts[ctxIndex].getRecvBuffer();

                if(mode==DA_FLAGS::WriteMode::SET_VALUES)
                {
                    for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++){
                    proc_id=recvProcList[recv_p];

                        for(unsigned int var=0;var<dof;var++)
                        {
                            for (unsigned int k = nodeRecvOffset[proc_id]; k < (nodeRecvOffset[proc_id] + nodeRecvCount[proc_id]); k++)
                            {
                                (vec+var*m_uiTotalNodalSz)[recvNodeSM[k]]=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
                            }
                        }

                    }

                }else if(mode ==DA_FLAGS::WriteMode::ADD_VALUES)
                {
                    for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++){
                    proc_id=recvProcList[recv_p];

                        for(unsigned int var=0;var<dof;var++)
                        {
                            for (unsigned int k = nodeRecvOffset[proc_id]; k < (nodeRecvOffset[proc_id] + nodeRecvCount[proc_id]); k++)
                            {
                                (vec+var*m_uiTotalNodalSz)[recvNodeSM[k]]+=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
                            }
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


    template <typename T>
    void DA::getElementNodalValues(const T*in, T* eleVecOut,unsigned int eleID,unsigned int dof) const
    {
        //assert((eleID>=m_uiMesh->getElementLocalBegin()) && (eleID<m_uiMesh->getElementLocalEnd()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();

        for(unsigned int var=0;var<dof;var++)
            m_uiMesh->getElementNodalValues(&in[var*m_uiTotalNodalSz],&eleVecOut[var*nPe],eleID);

    }


    template <typename T>
    void DA::eleVecToVecAccumilation(T*out, const T* eleVecIn,unsigned int eleID,unsigned int dof) const
    {
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();

        for(unsigned int var=0;var<dof;var++)
            m_uiMesh->computeElementalContribution(&eleVecIn[var*nPe],&out[var*m_uiTotalNodalSz],eleID);

    }

    template <typename T>
    void DA::setVectorByFunction(T* local,std::function<void(T,T,T,T*)>func,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        ot::TreeNode tmpOct;
        const ot::TreeNode* allElements=&(*(m_uiMesh->getAllElements().begin()));
        const unsigned int * e2n=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int * cgToDg=&(*(m_uiMesh->getCG2DGMap().begin()));

        unsigned int offset;
        unsigned int arrSz;

        unsigned int lookup=0;
        unsigned int ownerID, ii_x, jj_y, kk_z;
        double hx;
        const unsigned int eleOrder=m_uiMesh->getElementOrder();

        T x,y,z;
        T* val=new T[dof];

        if(!isElemental)
        {
            if(isGhosted) {
                offset = 0;
                arrSz = m_uiTotalNodalSz;
            }else{
                offset=m_uiMesh->getNodeLocalBegin();
                arrSz=m_uiLocalNodalSz;
            }


            for(unsigned int node=m_uiMesh->getNodeLocalBegin();node<m_uiMesh->getNodeLocalEnd();node++)
            {
                lookup=cgToDg[node];
                m_uiMesh->dg2eijk(lookup,ownerID,ii_x,jj_y,kk_z);
                tmpOct=allElements[ownerID];
                hx=(tmpOct.maxX()-tmpOct.minX())/((double) eleOrder);
                x=tmpOct.minX() + ii_x*(hx);
                y=tmpOct.minY() + jj_y*(hx);
                z=tmpOct.minZ() + kk_z*(hx);

                func(x,y,z,val);

                for(unsigned int var=0;var<dof;var++)
                    local[ ( var * arrSz ) + (node - offset) ]=val[var];

            }

        }else {

            if(isGhosted) {
                offset = 0;
                arrSz = m_uiTotalElementSz;
            }else{
                offset=m_uiMesh->getElementLocalBegin();
                arrSz=m_uiLocalElementSz;
            }


            for(unsigned int ele=m_uiMesh->getElementLocalBegin();ele<m_uiMesh->getElementLocalEnd();ele++)
            {
                tmpOct=allElements[ele];
                x=tmpOct.minX();
                y=tmpOct.minY();
                z=tmpOct.minZ();

                func(x,y,z,val);
                for(unsigned int var=0;var<dof;var++)
                    local[ (var*arrSz) + (ele-offset) ] = val[var];

            }

        }



        delete [] val;




    }


    template <typename T>
    void DA::setVectorByScalar(T* local,const T* value,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        ot::TreeNode tmpOct;
        const ot::TreeNode* allElements=&(*(m_uiMesh->getAllElements().begin()));
        unsigned int arrSz;

        if(!isElemental)
        {

            if(isGhosted) {
                arrSz = m_uiTotalNodalSz;
            }else{
                arrSz=m_uiLocalNodalSz;
            }

            for(unsigned int var=0;var<dof;var++)
            {
                for(unsigned int node=0;node<arrSz;node++)
                    local[ (var*arrSz) + node]=value[var];
            }


        }else{

            if(isGhosted) {
                arrSz = m_uiTotalElementSz;
            }else{
                arrSz=m_uiLocalElementSz;
            }

            for(unsigned int var=0;var<dof;var++)
            {
                for(unsigned int ele=0;ele<arrSz;ele++)
                    local[ (var*arrSz) + ele]=value[var];
            }

        }


    }

    template<typename T>
    T* DA::getVecPointerToDof(T* in ,unsigned int dofInex, bool isElemental,bool isGhosted) const
    {

        if(!(m_uiMesh->isActive()))
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


    template <typename T>
    void DA::vecTopvtu(T* local, const char * fPrefix,char** nodalVarNames,bool isElemental,bool isGhosted,unsigned int dof) 
    {

        T* _local=NULL;
        std::vector<std::string> _nodalVarNames;
        bool isDeafaultVarNames=false;
        if(nodalVarNames==NULL)
            isDeafaultVarNames=true;

        if(isGhosted)
            _local=local;
        else{
            createVector(_local,isElemental,true,dof);
            nodalVecToGhostedNodal(local,_local,false,dof);
            std::string vPrefix("m_uiVar_");

            readFromGhostBegin(_local,dof);
            readFromGhostEnd(_local,dof);

            if(nodalVarNames==NULL)
            {
                nodalVarNames=new char*[dof];
                _nodalVarNames.resize(dof);
                for(unsigned int var=0;var<dof;var++) {
                    _nodalVarNames[var] = (vPrefix + std::to_string(var));
                    nodalVarNames[var] = (char *)(_nodalVarNames[var].c_str());
                }

            }

        }



        T** nodalVec=new T*[dof];
        for(unsigned int var=0;var<dof;var++)
            nodalVec[var]=_local+var*m_uiTotalNodalSz;

        io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,0,NULL,NULL,dof,(const char**)nodalVarNames,(const double **)nodalVec);

        // release mem, only if the input vec is ghosted.
        if(!isGhosted)
            destroyVector(_local);

        // delete the pointers to multiple vars.
        delete [] nodalVec;

        // if used the default varNames delete them.

        if(isDeafaultVarNames)
        {
            _nodalVarNames.clear();
            delete [] nodalVarNames;
            nodalVarNames=NULL;
        }
        

    }


    template<typename T>
    void DA::copyVector(T* dest,const T* source,bool isElemental,bool isGhosted) const
    {
        if(!(m_uiMesh->isActive()))
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
    void DA::copyVectors(T* dest,const T* source,bool isElemental,bool isGhosted,unsigned int dof) const
    {
        if(!(m_uiMesh->isActive()))
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


    template<typename T>
    void DA::intergridTransfer(const T* varIn, T*& varOut,const ot::DA* newDA, bool isElemental, bool isGhosted, unsigned int dof,INTERGRID_TRANSFER_MODE mode) 
    {
        T** vIn =new T*[dof];
        const ot::Mesh* newMesh=newDA->getMesh();
        const unsigned int zipSzNew = newMesh->getDegOfFreedom();
        const unsigned int zipSzOld = m_uiMesh->getDegOfFreedom();

        const unsigned int nLocalNodes_old = m_uiMesh->getNumLocalMeshNodes();
        const unsigned int nLocalNodes_new = newMesh->getNumLocalMeshNodes();
        
        for(unsigned int var=0;var<dof;var++)
            vIn[var]=NULL;


        if(!isGhosted)
        {
            for(unsigned int var=0;var<dof;var++)
            {
                // allocation happens inside the function
                this->nodalVecToGhostedNodal(varIn + var*nLocalNodes_old , vIn[var],false , 1);
            }

        }else{

                for(unsigned int var=0;var<dof;var++)
                {
                    this->createVector(vIn[var],false,true,1);
                    std::memcpy(vIn[var],varIn + var*zipSzOld, sizeof(T)*zipSzOld);
                }
        }    

        // ghost exchange initiates
        for(unsigned int var=0;var<dof;var++)
            this->readFromGhostBegin(vIn[var],1);

        for(unsigned int var=0;var<dof;var++)
        {
            // wait till current var ghost exchange completes
            this->readFromGhostEnd(vIn[var],1);
            m_uiMesh->interGridTransfer(vIn[var],newMesh,mode);
        }
           

        if(!isGhosted)
        {

            newDA->createVector(varOut,false,false,dof);
            for(unsigned int var=0;var<dof;var++)
            {
                T* tmpPtr=(varOut+var*nLocalNodes_new);
                newDA->ghostedNodalToNodalVec((const T*)vIn[var],tmpPtr,true,1);
            }
                

        }else{

            newDA->createVector(varOut,false,true,dof);
            for(unsigned int var=0;var<dof;var++)
                std::memcpy(varOut + var*zipSzNew,vIn[var],sizeof(T)*zipSzNew);

        }

        for(unsigned int var=0;var<dof;var++)
            delete [] vIn[var];
        
        delete [] vIn;

    }


    template<typename T>
    int DA::getFaceNeighborValues(unsigned int eleID, const T* in, T* out, T* coords, unsigned int * neighID, unsigned int face, NeighbourLevel & level,unsigned int dof) const
    {
        if(!m_uiMesh->isActive())
            return false;

        int state=0;
        const unsigned int totalNodeSz = m_uiMesh->getDegOfFreedom();
        const unsigned int nPe = m_uiMesh->getNumNodesPerElement();
        const unsigned int outSz = 4*nPe;

        for(unsigned int v=0;v<dof;v++)
            state=m_uiMesh->getFaceNeighborValues(eleID,in + v*totalNodeSz, out + v * outSz,coords,neighID,face,level);
        
        return state;

    }


#ifdef BUILD_WITH_PETSC


    template<typename T>
    void DA::petscSetVectorByFunction(Vec& local,std::function<void(T,T,T,T*)>func,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        PetscScalar * arry=NULL;
        VecGetArray(local,&arry);

        setVectorByFunction(arry,func,isElemental,isGhosted,dof);

        VecRestoreArray(local,&arry);


    }

    template <typename T>
    void DA::petscSetVectorByScalar(Vec& local,const T* value,bool isElemental, bool isGhosted, unsigned int dof) const
    {

        PetscScalar * arry=NULL;
        VecGetArray(local,&arry);

        setVectorByScalar(arry,value,isElemental,isGhosted,dof);

        VecRestoreArray(local,&arry);

    }

    template <typename T>
void DA::petscIntergridTransfer(const Vec & varIn, Vec & varOut, const ot::DA* newDA, bool isElemental, bool isGhosted, unsigned int dof,INTERGRID_TRANSFER_MODE mode) {

        T** vIn =new T*[dof];
        const ot::Mesh* newMesh=newDA->getMesh();
        const unsigned int zipSzNew = newMesh->getDegOfFreedom();
        const unsigned int zipSzOld = m_uiMesh->getDegOfFreedom();

        const unsigned int nLocalNodes_old = m_uiMesh->getNumLocalMeshNodes();
        const unsigned int nLocalNodes_new = newMesh->getNumLocalMeshNodes();

        for(unsigned int var=0;var<dof;var++)
            vIn[var]=NULL;

        const T * inArray;
        if (this->isActive()) {
          VecGetArrayRead(varIn, &inArray);
        }


        if(!isGhosted)
        {
            for(unsigned int var=0;var<dof;var++)
            {
                // allocation happens inside the function
                this->nodalVecToGhostedNodal(inArray + var*nLocalNodes_old , vIn[var],false , 1);
            }

        }else{

            for(unsigned int var=0;var<dof;var++)
            {
                this->createVector(vIn[var],false,true,1);
                std::memcpy(vIn[var],inArray + var*zipSzOld, sizeof(T)*zipSzOld);
            }
        }

        // ghost exchange
        for(unsigned int var=0;var<dof;var++)
            this->readFromGhostBegin(vIn[var],1);

        for(unsigned int var=0;var<dof;var++)
            this->readFromGhostEnd(vIn[var],1);

        for(unsigned int var=0;var<dof;var++)
            m_uiMesh->interGridTransfer(vIn[var],newMesh,mode);

        if(newDA->isActive()) {
          if (!isGhosted) {
            newDA->petscCreateVector(varOut, false, false, dof);
            T *outArray;
            VecGetArray(varOut, &outArray);
            for (unsigned int var = 0; var < dof; var++) {
              T *tmpPtr = (outArray + var * nLocalNodes_new);
              newDA->ghostedNodalToNodalVec((const T *) vIn[var], tmpPtr, true, 1);
            }
            VecRestoreArray(varOut, &outArray);
        } else {
            newDA->petscCreateVector(varOut, false, true, dof);
            T *outArray;
            VecGetArray(varOut, &outArray);
            for (unsigned int var = 0; var < dof; var++)
              std::memcpy(outArray + var * zipSzNew, vIn[var], sizeof(T) * zipSzNew);
            VecRestoreArray(varOut, &outArray);
          }
        }
        
        if (this->isActive()) {
          VecRestoreArrayRead(varIn, &inArray);
        }


      for(unsigned int var=0;var<dof;var++)
            delete [] vIn[var];
        delete [] vIn;
    }


#endif



}



