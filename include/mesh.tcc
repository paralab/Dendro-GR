//
// Created by milinda on 9/22/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains templated functions in the mesh class.
 * (code moved from mesh.h)
*/
//


namespace ot
{
    template <typename T>
    T* Mesh::createVector() const
    {

        if(!m_uiIsActive) return NULL;

        T* vec=NULL;

        try
        {
            vec = new T[m_uiNumActualNodes];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }
        
        return vec;

    }

    template<typename T>
    T* Mesh::createCGVector(T initVal, unsigned int dof) const
    {
        if(!m_uiIsActive) return NULL;

        T* vec =NULL;
        try
        {
            vec=new T[m_uiNumActualNodes*dof];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }

        for(unsigned int i=0; i < m_uiNumActualNodes*dof; i++)
            vec[i] = initVal;
        
        return vec;
    }

    template<typename T>
    T* Mesh::createCGVector(std::function<void(T,T,T,T*)> func, unsigned int dof) const
    {
        if(!m_uiIsActive) return NULL;

        T* vec =NULL;
        try
        {
            vec=new T[m_uiNumActualNodes*dof];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }

        for(unsigned int i=0; i < m_uiNumActualNodes*dof; i++)
            vec[i] = (T)0;
        
        // initialize the vector to the function. 
        T* fvar =new T[dof];
        
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));
        for(unsigned int elem=m_uiElementLocalBegin;elem<m_uiElementLocalEnd;elem++)
        {

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++ )
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        const unsigned int nodeLookUp_CG=m_uiE2NMapping_CG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(nodeLookUp_CG>=m_uiNodeLocalBegin && nodeLookUp_CG<m_uiNodeLocalEnd)
                        {
                            unsigned int ownerID,ii_x,jj_y,kk_z;
                            const unsigned int nodeLookUp_DG=m_uiE2NMapping_DG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            const unsigned int len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                            
                            const double x=pNodes[ownerID].getX()+ ii_x*(len/((double)m_uiElementOrder));
                            const double y=pNodes[ownerID].getY()+ jj_y*(len/((double)m_uiElementOrder));
                            const double z=pNodes[ownerID].getZ()+ kk_z*(len/((double)m_uiElementOrder));
                            
                            Point physical_coord;
                            this->octCoordToDomainCoord(Point(x,y,z),physical_coord);
                            func(physical_coord.x(),physical_coord.y(),physical_coord.z(),fvar);


                            for(unsigned int v=0; v < dof; v++)
                                vec[v*m_uiNumActualNodes + nodeLookUp_CG]=fvar[v];
                        }

                    }

        }


        delete [] fvar;
        return vec;

    }

    template <typename T>
    T* Mesh::createElementVector(T initVal, unsigned int dof) const
    {

        if(!m_uiIsActive) return NULL;

        T* vec=NULL;

        try
        {
            vec=new T[m_uiNumTotalElements*dof];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }
        
        for(unsigned int i=0; i < m_uiNumTotalElements*dof; i++)
            vec[i] = initVal;
        
        return vec;
    }

    template <typename T>
    T* Mesh::createDGVector(T initVal, unsigned int dof) const
    {
        if(!m_uiIsActive) return NULL;
        
        T* vec=NULL;
        try
        {
            vec = new T[m_uiNumTotalElements*m_uiNpE*dof];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }

         for(unsigned int i=0; i < m_uiNumTotalElements*m_uiNpE*dof; i++)
            vec[i] = initVal;

        return vec;
    }

    template <typename T>
    T* Mesh::createDGVector(std::function<void(T,T,T,T*)> func, unsigned int dof) const
    {
        if(!m_uiIsActive) return NULL;
        
        T* vec=NULL;
        try
        {
            vec = new T[m_uiNumTotalElements*m_uiNpE*dof];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }

        for(unsigned int i=0; i < m_uiNumTotalElements*m_uiNpE*dof; i++)
            vec[i] = (T)0;

        T* fvar =new T[dof];
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));
        for(unsigned int elem=m_uiElementLocalBegin;elem<m_uiElementLocalEnd;elem++)
        {

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++ )
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        const unsigned int dg_index = elem * m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j*(m_uiElementOrder+1) + i;
                        const unsigned int len=1u<<(m_uiMaxDepth-pNodes[elem].getLevel());

                        const double x=pNodes[elem].getX()+ i*(len/((double)m_uiElementOrder));
                        const double y=pNodes[elem].getY()+ j*(len/((double)m_uiElementOrder));
                        const double z=pNodes[elem].getZ()+ k*(len/((double)m_uiElementOrder));
                        
                        Point physical_coord;
                        this->octCoordToDomainCoord(Point(x,y,z),physical_coord);
                        func(physical_coord.x(),physical_coord.y(),physical_coord.z(),fvar);


                        for(unsigned int v=0; v < dof; v++)
                            vec[ v*m_uiNumTotalElements*m_uiNpE + dg_index] = fvar[v];

                    }
        }

        delete [] fvar;
        return vec;

    }


    template <typename T>
    T* Mesh::createVector(const T initValue) const
    {
        if(!m_uiIsActive) return NULL;

        T* vec=NULL;

        try
        {
            vec=new T[m_uiNumActualNodes];
        }
        catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }


        for(unsigned int k=0;k<m_uiNumActualNodes;k++)
            vec[k]=initValue;

        return vec;
    }

    template <typename T>
    T* Mesh::createVector(std::function<T(T,T,T)> func) const
    {
        if(!m_uiIsActive) return NULL;

        T* vec = NULL;

        try
        {
            vec=new T[m_uiNumActualNodes];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }


        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        unsigned int len;
        double x,y,z;
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;

        for(unsigned int elem=m_uiElementLocalBegin;elem<m_uiElementLocalEnd;elem++)
        {


            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++ )
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        nodeLookUp_CG=m_uiE2NMapping_CG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(nodeLookUp_CG>=m_uiNodeLocalBegin && nodeLookUp_CG<m_uiNodeLocalEnd)
                        {
                            nodeLookUp_DG=m_uiE2NMapping_DG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                            x=pNodes[ownerID].getX()+ ii_x*(len/((double)m_uiElementOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/((double)m_uiElementOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/((double)m_uiElementOrder));
                            vec[nodeLookUp_CG]=func(x,y,z);
                        }

                    }

        }

        return vec;
    }

    template <typename T>
    void Mesh::createVector(std::vector<T> &vec) const
    {
        if(!m_uiIsActive) { vec.clear(); return; }
        vec.resize(m_uiNumActualNodes);

    }

    template <typename T>
    void Mesh::createVector(std::vector<T> &vec,const T initValue) const
    {
        if(!m_uiIsActive) { vec.clear(); return; }
        vec.resize(m_uiNumActualNodes,initValue);

    }

    template  <typename T>
    void Mesh::createVector(std::vector<T> & vec, std::function<T(T,T,T)> func) const
    {
        if(!m_uiIsActive) { vec.clear(); return; }
        vec.clear();
        vec.resize(m_uiNumActualNodes,0);
        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        unsigned int len;
        double x,y,z;
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;

        for(unsigned int elem=m_uiElementLocalBegin;elem<m_uiElementLocalEnd;elem++)
        {


            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++ )
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        nodeLookUp_CG=m_uiE2NMapping_CG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(nodeLookUp_CG>=m_uiNodeLocalBegin && nodeLookUp_CG<m_uiNodeLocalEnd)
                        {
                            nodeLookUp_DG=m_uiE2NMapping_DG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                            x=pNodes[ownerID].getX()+ ii_x*(len/((double)m_uiElementOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/((double)m_uiElementOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/((double)m_uiElementOrder));
                            vec[nodeLookUp_CG]=func(x,y,z);
                        }

                    }

        }

    }

    template <typename T>
    void Mesh::createUnZippedVector(std::vector<T> & uvec) const
    {
        if(!m_uiIsActive) { uvec.clear(); return; }
        uvec.resize(m_uiUnZippedVecSz);
    }

    template <typename T>
    void Mesh::createUnZippedVector(std::vector<T> & uvec,const T initValue) const
    {
        if(!m_uiIsActive) { uvec.clear(); return; }
        uvec.resize(m_uiUnZippedVecSz,initValue);
    }


    template <typename T>
    T* Mesh::createUnZippedVector(unsigned int dof) const
    {
        if(!m_uiIsActive) return NULL;

        T* uvec=NULL;

        try
        {
            uvec=new T[dof*m_uiUnZippedVecSz];
            
        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }

        return uvec;
    }

    template <typename T>
    T* Mesh::createUnZippedVector(const T initValue,unsigned int dof) const
    {
        if(!m_uiIsActive) return NULL;

        T* uvec=NULL;
        try
        {
            uvec=new T[dof*m_uiUnZippedVecSz];

        }catch (const std::bad_alloc& e) {
            std::cout<<" rank: "<<m_uiActiveRank<<" func: "<<__func__<<" bad allocation error "<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
	    }

        for(unsigned int k=0;k<(dof*m_uiUnZippedVecSz);k++)
            uvec[k]=initValue;

        return uvec;
    }

    template <typename T>
    void Mesh::CG2DGVec(T* cg_vec, T* dg_vec, bool gsynced, unsigned int dof) 
    {
        if(!m_uiIsActive)
            return;

        if(!gsynced)
        {
            this->readFromGhostBegin(cg_vec,dof);
            this->readFromGhostEnd(cg_vec,dof);

        }

        const unsigned int vsz_dg = m_uiNumTotalElements * m_uiNpE; 
        const unsigned int vsz_cg = m_uiNumActualNodes;

        for(unsigned int v=0; v < dof; v++ )
            for(unsigned int ele = m_uiElementLocalBegin; ele < m_uiElementLocalEnd; ele++)
                this->getElementNodalValues(cg_vec + v*vsz_cg , dg_vec + v * vsz_dg + ele*m_uiNpE, ele);
        
        return ;

        
    }
    
    template <typename T>
    void Mesh::DG2CGVec(const T* dg_vec, T* cg_vec, unsigned int dof) const
    {
        if(!m_uiIsActive)
            return;

        const unsigned int vsz_dg = m_uiNumTotalElements * m_uiNpE; 
        const unsigned int vsz_cg = m_uiNumActualNodes;

        bool isHanging;
        unsigned int cnum;
        for(unsigned int v=0; v < dof; v++ )
        {

            for(unsigned int ele = m_uiElementLocalBegin; ele < m_uiElementLocalEnd; ele++)
            {
            
                for(unsigned int k=0; k < (m_uiElementOrder+1); k++)
                for(unsigned int j=0; j < (m_uiElementOrder+1); j++)
                for(unsigned int i=0; i < (m_uiElementOrder+1); i++)
                {
                    isHanging = isNodeHanging(ele,i,j,k);
                    if(!isHanging)
                        cg_vec[v * vsz_cg + m_uiE2NMapping_CG[ele*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j* (m_uiElementOrder+1)+ i]] = dg_vec[ v *vsz_dg + ele*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j* (m_uiElementOrder+1)+ i];
                    else
                    {
                        cnum=m_uiAllElements[(ele)].getMortonIndex();
                        const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                        const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                        const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;

                        if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            cg_vec[v * vsz_cg + m_uiE2NMapping_CG[ele*m_uiNpE + (kkz>>1u)*(m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+ (iix>>1u)]] = dg_vec[ v *vsz_dg + ele*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j* (m_uiElementOrder+1)+ i];

                    }
                        
                }
            
            }

        }

        return ;


    }

    template <typename T>
    void Mesh::DG2CGVec(const T* dg_vec, T*& cg_vec, bool isAllocated, const unsigned int* eleIDs, unsigned int nEle, unsigned int dof) const
    {
         if(!m_uiIsActive)
            return;

        if(!isAllocated)
            cg_vec = this->createCGVector((T)0 , dof);
        
        const unsigned int vsz_dg = m_uiNumTotalElements * m_uiNpE; 
        const unsigned int vsz_cg = m_uiNumActualNodes;

        bool isHanging;
        unsigned int cnum;
        for(unsigned int v=0; v < dof; v++ )
        {

            for(unsigned int i = 0; i < nEle; i++)
            {
                const unsigned int ele = eleIDs[i];
                assert(ele<m_uiAllElements.size());
                
                for(unsigned int k=0; k < (m_uiElementOrder+1); k++)
                for(unsigned int j=0; j < (m_uiElementOrder+1); j++)
                for(unsigned int i=0; i < (m_uiElementOrder+1); i++)
                {
                    isHanging = isNodeHanging(ele,i,j,k);
                    if(!isHanging)
                        cg_vec[v * vsz_cg + m_uiE2NMapping_CG[ele*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j* (m_uiElementOrder+1)+ i]] = dg_vec[ v *vsz_dg + ele*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j* (m_uiElementOrder+1)+ i];
                    else
                    {
                        cnum=m_uiAllElements[(ele)].getMortonIndex();
                        const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                        const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                        const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;

                        if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            cg_vec[v * vsz_cg + m_uiE2NMapping_CG[ele*m_uiNpE + (kkz>>1u)*(m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+ (iix>>1u)]] = dg_vec[ v *vsz_dg + ele*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j* (m_uiElementOrder+1)+ i];

                    }
                        
                }
            
            }

        }

    }

    template<typename T>
    void Mesh::performGhostExchange(std::vector<T> &vec)
    {

        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        for(unsigned int p=0;p<m_uiActiveNpes;p++) {
            for (unsigned int k = m_uiSendNodeOffset[p]; k < (m_uiSendNodeOffset[p] + m_uiSendNodeCount[p]); k++) {
                m_uiSendBufferNodes[k] = (T)vec[m_uiScatterMapActualNodeSend[k]];
            }
        }

        #ifdef ALLTOALL_SPARSE
            par::Mpi_Alltoallv_sparse(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #else
            par::Mpi_Alltoallv(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #endif


        for(unsigned int p=0;p<m_uiActiveNpes;p++)
        {
            for(unsigned int k=m_uiRecvNodeOffset[p];k<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p]);k++)
            {

                //if(fabs(vec[m_uiScatterMapActualNodeRecv[k]]-m_uiRecvBufferNodes[k])>1e-15) std::cout<<"rank: "<<m_uiActiveRank<<" computed: "<<vec[m_uiScatterMapActualNodeRecv[k]]<<" revieved: "<<m_uiRecvBufferNodes[k]<<" recv: from : "<<p<<std::endl;
                vec[m_uiScatterMapActualNodeRecv[k]]=(T)m_uiRecvBufferNodes[k];
            }

        }




    }

    template<typename T>
    void Mesh::performGhostExchange(T*vec)
    {

        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        for(unsigned int p=0;p<m_uiActiveNpes;p++) {
            for (unsigned int k = m_uiSendNodeOffset[p]; k < (m_uiSendNodeOffset[p] + m_uiSendNodeCount[p]); k++) {
                m_uiSendBufferNodes[k] = (T)vec[m_uiScatterMapActualNodeSend[k]];
            }
        }

        #ifdef ALLTOALL_SPARSE
            par::Mpi_Alltoallv_sparse(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #else
            par::Mpi_Alltoallv(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #endif


        for(unsigned int p=0;p<m_uiActiveNpes;p++)
        {
            for(unsigned int k=m_uiRecvNodeOffset[p];k<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p]);k++)
            {

                //if(/*fabs(vec[m_uiScatterMapActualNodeRecv[k]]-m_uiRecvBufferNodes[k])>1e-15*/ isnan(m_uiRecvBufferNodes[k])) std::cout<<"rank: "<<m_uiActiveRank<<" computed: "<<vec[m_uiScatterMapActualNodeRecv[k]]<<" revieved: "<<m_uiRecvBufferNodes[k]<<" recv: from : "<<p<<std::endl;
                vec[m_uiScatterMapActualNodeRecv[k]]=(T)m_uiRecvBufferNodes[k];
            }

        }




    }

    template<typename T>
    void Mesh::ghostExchangeStart(T* vec,T* sendNodeBuffer,T* recvNodeBuffer, MPI_Request * send_reqs, MPI_Request * recv_reqs)
    {

        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        unsigned int proc_id;

        // active recv procs
        for(unsigned int recv_p=0;recv_p<m_uiRecvProcList.size();recv_p++)
        {
            proc_id=m_uiRecvProcList[recv_p];
            recv_reqs[recv_p]=MPI_Request();
            par::Mpi_Irecv((recvNodeBuffer+m_uiRecvNodeOffset[proc_id]),m_uiRecvNodeCount[proc_id],proc_id,0,m_uiCommActive,&recv_reqs[recv_p]);
        }

        for(unsigned int send_p=0;send_p<m_uiSendProcList.size();send_p++)
        {
            proc_id=m_uiSendProcList[send_p];
            for (unsigned int k = m_uiSendNodeOffset[proc_id]; k < (m_uiSendNodeOffset[proc_id] + m_uiSendNodeCount[proc_id]); k++) {
                sendNodeBuffer[k] = (T)vec[m_uiScatterMapActualNodeSend[k]];
            }
        }
        // active send procs
        for(unsigned int send_p=0;send_p<m_uiSendProcList.size();send_p++)
        {
            proc_id=m_uiSendProcList[send_p];
            send_reqs[send_p]=MPI_Request();
            par::Mpi_Isend((sendNodeBuffer+m_uiSendNodeOffset[proc_id]),m_uiSendNodeCount[proc_id],proc_id,0,m_uiCommActive,&send_reqs[send_p]);
        }



    }


    template<typename T>
    void Mesh::ghostExchangeRecvSync(T *vec, T *recvNodeBuffer,MPI_Request *recv_reqs, MPI_Status *recv_sts)
    {
        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        dendro::timer::t_unzip_async_comm.start();
        MPI_Waitall(m_uiRecvProcList.size(),recv_reqs,recv_sts);
        dendro::timer::t_unzip_async_comm.stop();

        unsigned int proc_id=0;
        for(unsigned int recv_p=0;recv_p<m_uiRecvProcList.size();recv_p++)
        {
            proc_id=m_uiRecvProcList[recv_p];
            for(unsigned int k=m_uiRecvNodeOffset[proc_id];k<(m_uiRecvNodeOffset[proc_id]+m_uiRecvNodeCount[proc_id]);k++)
            {
                vec[m_uiScatterMapActualNodeRecv[k]]=(T)recvNodeBuffer[k];
            }
        }

    }

    template<typename T>
    void Mesh::readFromGhostBegin(T* vec, unsigned int dof)
    {
         if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            const std::vector<unsigned int>& nodeSendCount=this->getNodalSendCounts();
            const std::vector<unsigned int>& nodeSendOffset=this->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeRecvCount=this->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeRecvOffset=this->getNodalRecvOffsets();

            const std::vector<unsigned int>& sendProcList=this->getSendProcList();
            const std::vector<unsigned int>& recvProcList=this->getRecvProcList();

            const std::vector<unsigned int>& sendNodeSM=this->getSendNodeSM();
            const std::vector<unsigned int>& recvNodeSM=this->getRecvNodeSM();


            const unsigned int activeNpes=this->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            AsyncExchangeContex ctx(vec);
            MPI_Comm commActive=this->getMPICommunicator();


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
                            sendB[dof*(nodeSendOffset[proc_id]) + (var*nodeSendCount[proc_id])+(k-nodeSendOffset[proc_id])] = (vec+var*m_uiNumActualNodes)[sendNodeSM[k]];
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
    void Mesh::readFromGhostEnd(T* vec, unsigned int dof)
    {
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            const std::vector<unsigned int>& nodeSendCount=this->getNodalSendCounts();
            const std::vector<unsigned int>& nodeSendOffset=this->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeRecvCount=this->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeRecvOffset=this->getNodalRecvOffsets();

            const std::vector<unsigned int>& sendProcList=this->getSendProcList();
            const std::vector<unsigned int>& recvProcList=this->getRecvProcList();

            const std::vector<unsigned int>& sendNodeSM=this->getSendNodeSM();
            const std::vector<unsigned int>& recvNodeSM=this->getRecvNodeSM();


            const unsigned int activeNpes=this->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            int ctxIndex=-1;
            for(unsigned int i=0;i<m_uiMPIContexts.size();i++)
            {
                if(m_uiMPIContexts[i].getBuffer()==vec)
                {
                    ctxIndex=i;
                    break;
                }

            }

            if(ctxIndex==-1)
            {
                std::cout<<"rank: "<<m_uiActiveRank<<" async ctx not found for vec: "<<&vec<<" in async comm end: "<<__LINE__<<std::endl;
                MPI_Abort(m_uiCommActive,0);
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
                            (vec+var*m_uiNumActualNodes)[recvNodeSM[k]]=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
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
    void Mesh::readFromGhostBegin(AsyncExchangeContex& ctx, T* vec, unsigned int dof)
    {
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;
        

        if(this->isActive())
        {
            const std::vector<unsigned int>& nodeSendCount=this->getNodalSendCounts();
            const std::vector<unsigned int>& nodeSendOffset=this->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeRecvCount=this->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeRecvOffset=this->getNodalRecvOffsets();

            const std::vector<unsigned int>& sendProcList=this->getSendProcList();
            const std::vector<unsigned int>& recvProcList=this->getRecvProcList();

            const std::vector<unsigned int>& sendNodeSM=this->getSendNodeSM();
            const std::vector<unsigned int>& recvNodeSM=this->getRecvNodeSM();


            const unsigned int activeNpes=this->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            MPI_Comm commActive=this->getMPICommunicator();


            if(recvBSz)
            {
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++)
                {
                    proc_id=recvProcList[recv_p];
                    par::Mpi_Irecv((recvB+dof*nodeRecvOffset[proc_id]),dof*nodeRecvCount[proc_id],proc_id,m_uiCommTag,commActive,&ctx.m_recv_req[recv_p]);
                }

            }

            if(sendBSz)
            {
                sendB=(T*)ctx.getSendBuffer();
                for(unsigned int send_p=0;send_p<sendProcList.size();send_p++) {
                    proc_id=sendProcList[send_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = nodeSendOffset[proc_id]; k < (nodeSendOffset[proc_id] + nodeSendCount[proc_id]); k++)
                        {
                            sendB[dof*(nodeSendOffset[proc_id]) + (var*nodeSendCount[proc_id])+(k-nodeSendOffset[proc_id])] = (vec+var*m_uiNumActualNodes)[sendNodeSM[k]];
                        }

                    }
                }

                // active send procs
                for(unsigned int send_p=0;send_p<sendProcList.size();send_p++)
                {
                    proc_id=sendProcList[send_p];
                    par::Mpi_Isend(sendB+dof*nodeSendOffset[proc_id],dof*nodeSendCount[proc_id],proc_id,m_uiCommTag,commActive,&ctx.m_send_req[send_p]);

                }
            
            }

            m_uiCommTag++;
            
        }

        return;
    }

    template<typename T>
    void Mesh::readFromGhostEnd(AsyncExchangeContex& ctx, T* vec, unsigned int dof)
    {
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            const std::vector<unsigned int>& nodeSendCount=this->getNodalSendCounts();
            const std::vector<unsigned int>& nodeSendOffset=this->getNodalSendOffsets();

            const std::vector<unsigned int>& nodeRecvCount=this->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeRecvOffset=this->getNodalRecvOffsets();

            const std::vector<unsigned int>& sendProcList=this->getSendProcList();
            const std::vector<unsigned int>& recvProcList=this->getRecvProcList();

            const std::vector<unsigned int>& sendNodeSM=this->getSendNodeSM();
            const std::vector<unsigned int>& recvNodeSM=this->getRecvNodeSM();


            const unsigned int activeNpes=this->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            MPI_Status status;
            // need to wait for the commns to finish ...
            MPI_Waitall(sendProcList.size(),ctx.m_send_req.data(),MPI_STATUSES_IGNORE);
            MPI_Waitall(recvProcList.size(),ctx.m_recv_req.data(),MPI_STATUSES_IGNORE);
            
            if(recvBSz)
            {
                // copy the recv data to the vec
                recvB=(T*)ctx.getRecvBuffer();

                for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++){
                    proc_id=recvProcList[recv_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = nodeRecvOffset[proc_id]; k < (nodeRecvOffset[proc_id] + nodeRecvCount[proc_id]); k++)
                        {
                            (vec+var*m_uiNumActualNodes)[recvNodeSM[k]]=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
                        }
                    }

                }

            }
        }

        return;
    }

    template<typename T>
    void Mesh::readFromGhostBeginElementVec(T* vec, unsigned int dof)
    {
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            const unsigned int activeNpes=m_uiActiveNpes;
            const unsigned int sendBSz = m_uiSendEleOffset[activeNpes-1] + m_uiSendEleCount[activeNpes-1];
            const unsigned int recvBSz = m_uiRecvEleOffset[activeNpes-1] + m_uiRecvEleCount[activeNpes-1];
            unsigned int proc_id;

            AsyncExchangeContex ctx(vec);
            MPI_Comm commActive=this->getMPICommunicator();

            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz*dof));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p = 0 ;  recv_p <  m_uiElementRecvProcList.size(); recv_p++)
                {
                    proc_id=m_uiElementRecvProcList[recv_p];
                    MPI_Request* req=new MPI_Request();
                    par::Mpi_Irecv((recvB+dof*m_uiRecvEleOffset[proc_id]),dof*m_uiRecvEleCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*dof*sendBSz);
                sendB=(T*)ctx.getSendBuffer();

                for(unsigned int send_p = 0; send_p < m_uiElementSendProcList.size(); send_p++) {
                    proc_id=m_uiElementSendProcList[send_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiSendEleOffset[proc_id]; k < (m_uiSendEleOffset[proc_id] + m_uiSendEleCount[proc_id]); k++)
                        {
                            sendB[dof*(m_uiSendEleOffset[proc_id]) + (var*m_uiSendEleCount[proc_id])+(k-m_uiSendEleOffset[proc_id])] = (vec+var*m_uiNumTotalElements)[ m_uiElementLocalBegin +  m_uiScatterMapElementRound1[k]];
                        }

                    }

                }

                
                // active send procs
                for(unsigned int send_p = 0; send_p < m_uiElementSendProcList.size(); send_p++) 
                {
                    proc_id=m_uiElementSendProcList[send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+dof*m_uiSendEleOffset[proc_id],dof*m_uiSendEleCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }


            }


            m_uiCommTag++;
            m_uiMPIContexts.push_back(ctx);
        
        }

        return;
    
    
    }


    template<typename T>
    void Mesh::readFromGhostEndElementVec(T* vec, unsigned int dof)
    {

        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            const unsigned int activeNpes = m_uiActiveNpes ;
            const unsigned int sendBSz = m_uiSendEleOffset[activeNpes-1] + m_uiSendEleCount[activeNpes-1];
            const unsigned int recvBSz = m_uiRecvEleOffset[activeNpes-1] + m_uiRecvEleCount[activeNpes-1];
            unsigned int proc_id;

            
            int ctxIndex=-1;
            for(unsigned int i=0;i<m_uiMPIContexts.size();i++)
            {
                if(m_uiMPIContexts[i].getBuffer()==vec)
                {
                    ctxIndex=i;
                    break;
                }

            }
            
            if(ctxIndex==-1)
            {
                std::cout<<"rank: "<<m_uiActiveRank<<" async ctx not found for vec: "<<&vec<<" in async comm end: "<<__LINE__<<std::endl;
                MPI_Abort(m_uiCommActive,0);
            }

            assert(m_uiMPIContexts[ctxIndex].getBuffer()==vec);

            MPI_Status status;
            // need to wait for the commns to finish ...
            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++) {
                MPI_Wait(m_uiMPIContexts[ctxIndex].getRequestList()[i], &status);
            }

            if(recvBSz)
            {
                // copy the recv data to the vec
                recvB=(T*)m_uiMPIContexts[ctxIndex].getRecvBuffer();

                for(unsigned int recv_p = 0 ; recv_p < m_uiElementRecvProcList.size();recv_p++){

                    proc_id = m_uiElementRecvProcList[recv_p];
                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiRecvEleOffset[proc_id]; k < (m_uiRecvEleOffset[proc_id] + m_uiRecvEleCount[proc_id]); k++)
                         (vec+var*m_uiNumTotalElements)[m_uiGhostElementRound1Index[k]] = recvB[dof*(m_uiRecvEleOffset[proc_id]) + (var*m_uiRecvEleCount[proc_id])+(k-m_uiRecvEleOffset[proc_id])];
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
    void Mesh::readFromGhostBeginEleDGVec(T* vec, unsigned int dof)
    {
        
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            const unsigned int activeNpes=m_uiActiveNpes;
            const unsigned int sendBSz = (m_uiSendEleOffset[activeNpes-1] + m_uiSendEleCount[activeNpes-1])*m_uiNpE;
            const unsigned int recvBSz = (m_uiRecvEleOffset[activeNpes-1] + m_uiRecvEleCount[activeNpes-1])*m_uiNpE;
            unsigned int proc_id;

            AsyncExchangeContex ctx(vec);
            MPI_Comm commActive=this->getMPICommunicator();

            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz*dof));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p = 0 ;  recv_p <  m_uiElementRecvProcList.size(); recv_p++)
                {
                    proc_id=m_uiElementRecvProcList[recv_p];
                    MPI_Request* req=new MPI_Request();
                    par::Mpi_Irecv((recvB+ dof*m_uiNpE*m_uiRecvEleOffset[proc_id]), dof * m_uiNpE*m_uiRecvEleCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*dof*sendBSz);
                sendB=(T*)ctx.getSendBuffer();

                for(unsigned int send_p = 0; send_p < m_uiElementSendProcList.size(); send_p++) {
                    proc_id=m_uiElementSendProcList[send_p];

                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiSendEleOffset[proc_id]; k < (m_uiSendEleOffset[proc_id] + m_uiSendEleCount[proc_id]); k++)
                        {
                            for(unsigned int node =0; node < m_uiNpE; node++)
                              sendB[dof*m_uiNpE*(m_uiSendEleOffset[proc_id]) + (var*m_uiNpE*m_uiSendEleCount[proc_id]) + (k-m_uiSendEleOffset[proc_id])*m_uiNpE + node ] = (vec+var*m_uiNumTotalElements*m_uiNpE)[ (m_uiElementLocalBegin +  m_uiScatterMapElementRound1[k])*m_uiNpE + node];
                        }
                    }
                }

                
                // active send procs
                for(unsigned int send_p = 0; send_p < m_uiElementSendProcList.size(); send_p++) 
                {
                    proc_id=m_uiElementSendProcList[send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+dof*m_uiNpE*m_uiSendEleOffset[proc_id],dof*m_uiNpE*m_uiSendEleCount[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }


            }


            m_uiCommTag++;
            m_uiMPIContexts.push_back(ctx);
        
        }

        return;


    }
    
    template<typename T>
    void Mesh::readFromGhostEndEleDGVec(T* vec, unsigned int dof)
    {
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            const unsigned int activeNpes=m_uiActiveNpes;
            const unsigned int sendBSz = (m_uiSendEleOffset[activeNpes-1] + m_uiSendEleCount[activeNpes-1])*m_uiNpE;
            const unsigned int recvBSz = (m_uiRecvEleOffset[activeNpes-1] + m_uiRecvEleCount[activeNpes-1])*m_uiNpE;
            unsigned int proc_id;

            
            int ctxIndex=-1;
            for(unsigned int i=0;i<m_uiMPIContexts.size();i++)
            {
                if(m_uiMPIContexts[i].getBuffer()==vec)
                {
                    ctxIndex=i;
                    break;
                }

            }

            if(ctxIndex==-1)
            {
                std::cout<<"rank: "<<m_uiActiveRank<<" async ctx not found for vec: "<<&vec<<" in async comm end: "<<__LINE__<<std::endl;
                MPI_Abort(m_uiCommActive,0);
            }

            assert(m_uiMPIContexts[ctxIndex].getBuffer()==vec);

            MPI_Status status;
            // need to wait for the commns to finish ...
            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++) {
                MPI_Wait(m_uiMPIContexts[ctxIndex].getRequestList()[i], &status);
            }

            if(recvBSz)
            {
                // copy the recv data to the vec
                recvB=(T*)m_uiMPIContexts[ctxIndex].getRecvBuffer();

                for(unsigned int recv_p = 0 ; recv_p < m_uiElementRecvProcList.size();recv_p++){

                    proc_id = m_uiElementRecvProcList[recv_p];
                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = m_uiRecvEleOffset[proc_id]; k < (m_uiRecvEleOffset[proc_id] + m_uiRecvEleCount[proc_id]); k++)
                            for(unsigned int node =0; node < m_uiNpE; node ++)
                                (vec+var*m_uiNumTotalElements*m_uiNpE)[m_uiGhostElementRound1Index[k]*m_uiNpE+ node ] = recvB[dof*(m_uiRecvEleOffset[proc_id]*m_uiNpE) + (var*m_uiNpE*m_uiRecvEleCount[proc_id]) + (k-m_uiRecvEleOffset[proc_id])*m_uiNpE + node];
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
    void Mesh::writeFromGhostBegin(T* vec, unsigned int dof)
    {
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            // inverse direction of the read ghost 
            const std::vector<unsigned int>& nodeSendCount=this->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeSendOffset=this->getNodalRecvOffsets();

            const std::vector<unsigned int>& nodeRecvCount=this->getNodalSendCounts();
            const std::vector<unsigned int>& nodeRecvOffset=this->getNodalSendOffsets();

            const std::vector<unsigned int>& sendProcList=this->getRecvProcList();
            const std::vector<unsigned int>& recvProcList=this->getSendProcList();

            const std::vector<unsigned int>& sendNodeSM=this->getRecvNodeSM();
            const std::vector<unsigned int>& recvNodeSM=this->getSendNodeSM();


            const unsigned int activeNpes=this->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            AsyncExchangeContex ctx(vec);
            MPI_Comm commActive=this->getMPICommunicator();


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
                            sendB[dof*(nodeSendOffset[proc_id]) + (var*nodeSendCount[proc_id])+(k-nodeSendOffset[proc_id])] = (vec+var*m_uiNumActualNodes)[sendNodeSM[k]];
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
    void Mesh::writeFromGhostEnd(T* vec, ot::GWMode mode ,unsigned int dof)
    {
        if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            // inverse direction of the read ghost 
            const std::vector<unsigned int>& nodeSendCount=this->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeSendOffset=this->getNodalRecvOffsets();

            const std::vector<unsigned int>& nodeRecvCount=this->getNodalSendCounts();
            const std::vector<unsigned int>& nodeRecvOffset=this->getNodalSendOffsets();

            const std::vector<unsigned int>& sendProcList=this->getRecvProcList();
            const std::vector<unsigned int>& recvProcList=this->getSendProcList();

            const std::vector<unsigned int>& sendNodeSM=this->getRecvNodeSM();
            const std::vector<unsigned int>& recvNodeSM=this->getSendNodeSM();


            const unsigned int activeNpes=this->getMPICommSize();

            const unsigned int sendBSz=nodeSendOffset[activeNpes-1] + nodeSendCount[activeNpes-1];
            const unsigned int recvBSz=nodeRecvOffset[activeNpes-1] + nodeRecvCount[activeNpes-1];
            unsigned int proc_id;

            int ctxIndex=-1;
            for(unsigned int i=0;i<m_uiMPIContexts.size();i++)
            {
                if(m_uiMPIContexts[i].getBuffer()==vec)
                {
                    ctxIndex=i;
                    break;
                }

            }

            if(ctxIndex==-1)
            {
                std::cout<<"rank: "<<m_uiActiveRank<<" async ctx not found for vec: "<<&vec<<" in async comm end: "<<__LINE__<<std::endl;
                MPI_Abort(m_uiCommActive,0);
            }

            MPI_Status status;
            // need to wait for the commns to finish ...
            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++) {
                MPI_Wait(m_uiMPIContexts[ctxIndex].getRequestList()[i], &status);
            }

            if(mode == ot::GWMode::ACCUMILATE)
            {
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
                                (vec+var*m_uiNumActualNodes)[recvNodeSM[k]]+=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
                            }
                        }

                    }
                }
            }else
            {
                assert(mode == ot::GWMode::OVERWRITE);
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
                                (vec+var*m_uiNumActualNodes)[recvNodeSM[k]]=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
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

    template<typename T>
    void Mesh::gatherFromGhostBegin(T* vec, unsigned int dof)
    {
        this->writeFromGhostBegin(vec,dof);
    }   

    template<typename T>
    void Mesh::gatherFromGhostEnd(T*vec, std::vector<std::vector<T> >& gatherV, unsigned int dof)
    {

          if(this->getMPICommSizeGlobal()==1 || (!m_uiIsActive))
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(this->isActive())
        {
            // inverse direction of the read ghost 
            const std::vector<unsigned int>& nodeSendCount=this->getNodalRecvCounts();
            const std::vector<unsigned int>& nodeSendOffset=this->getNodalRecvOffsets();

            const std::vector<unsigned int>& nodeRecvCount=this->getNodalSendCounts();
            const std::vector<unsigned int>& nodeRecvOffset=this->getNodalSendOffsets();

            const std::vector<unsigned int>& sendProcList=this->getRecvProcList();
            const std::vector<unsigned int>& recvProcList=this->getSendProcList();

            const std::vector<unsigned int>& sendNodeSM=this->getRecvNodeSM();
            const std::vector<unsigned int>& recvNodeSM=this->getSendNodeSM();


            const unsigned int activeNpes=this->getMPICommSize();

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
                gatherV.resize(m_uiNumActualNodes);
                recvB=(T*)m_uiMPIContexts[ctxIndex].getRecvBuffer();
                std::vector<unsigned int > gcounts;
                gcounts.resize(m_uiNumActualNodes,0);

                for(unsigned int p = 0; p < recvProcList.size(); p ++)
                {
                    const unsigned int proc_id = recvProcList[p];
                    for (unsigned int k = nodeRecvOffset[proc_id]; k < (nodeRecvOffset[proc_id] + nodeRecvCount[proc_id]); k++)
                        gcounts[recvNodeSM[k]]++;
                }

                for(unsigned int k=m_uiNodeLocalBegin; k < m_uiNodeLocalEnd; k++)
                {
                    if(gcounts[k]>0)
                    {
                        gatherV[k].resize(gcounts[k]*dof);
                        gcounts[k]=0;
                    }
                        
                }

                for(unsigned int p = 0; p < recvProcList.size(); p ++)
                {
                    const unsigned int proc_id = recvProcList[p];
                    for (unsigned int k = nodeRecvOffset[proc_id]; k < (nodeRecvOffset[proc_id] + nodeRecvCount[proc_id]); k++)
                    {
                        const unsigned int gsz = gatherV[recvNodeSM[k]].size();
                        for(unsigned int v =0; v < dof; v++)
                        {
                            gatherV[recvNodeSM[k]][v * gsz + gcounts[recvNodeSM[k]]] = recvB[dof*(nodeRecvOffset[proc_id]) + (v*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
                            gcounts[recvNodeSM[k]]+=1;
                        }
                        
                    }
                }


                for(unsigned int recv_p=0;recv_p<recvProcList.size();recv_p++){

                    const unsigned int proc_id=recvProcList[recv_p];
                    for(unsigned int var=0;var<dof;var++)
                    {
                        for (unsigned int k = nodeRecvOffset[proc_id]; k < (nodeRecvOffset[proc_id] + nodeRecvCount[proc_id]); k++)
                            (vec+var*m_uiNumActualNodes)[recvNodeSM[k]]=recvB[dof*(nodeRecvOffset[proc_id]) + (var*nodeRecvCount[proc_id])+(k-nodeRecvOffset[proc_id])];
                    }
                }

                

            }

            m_uiMPIContexts[ctxIndex].deAllocateSendBuffer();
            m_uiMPIContexts[ctxIndex].deAllocateRecvBuffer();

            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++)
                delete m_uiMPIContexts[ctxIndex].getRequestList()[i];

            m_uiMPIContexts[ctxIndex].getRequestList().clear();
            m_uiMPIContexts.erase(m_uiMPIContexts.begin() + ctxIndex);

        }



    }


    template <typename T,unsigned int length,unsigned int offsetCentered,unsigned int offsetBackward,unsigned int offsetForward>
    void Mesh::applyStencil(const std::vector<T> &in, std::vector<T> &out,
                            const Stencil<T, length, offsetCentered> &centered,
                            const Stencil<T, length, offsetBackward> &backward,
                            const Stencil<T, length, offsetForward> &forward)
    {


        if(!m_uiIsActive) return;

        double t_uzip;
        double t_uzip_g[3];

        double t_zip;
        double t_zip_g[3];

        double t_stencil;
        double t_stencil_g[3];


        unsigned int blkNpe_1D;
        std::vector<T> unzipVec;
        createUnZippedVector(unzipVec);

        std::vector<T> unzipVec1;
        this->createUnZippedVector(unzipVec1,0.0);

        #ifdef PROFILE_APPLY_STENCIL
                auto t1=std::chrono::high_resolution_clock::now();
        #endif
                this->unzip(&(*(in.begin())),&(*(unzipVec.begin())));
                //std::cout<<"rank: "<<m_uiActiveRank<<" unzip completed "<<std::endl;

        #ifdef PROFILE_APPLY_STENCIL
                auto t2=std::chrono::high_resolution_clock::now();
                t_uzip=std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

                par::Mpi_Reduce(&t_uzip,t_uzip_g,1,MPI_MIN,0,m_uiCommActive);
                par::Mpi_Reduce(&t_uzip,t_uzip_g+1,1,MPI_SUM,0,m_uiCommActive);
                par::Mpi_Reduce(&t_uzip,t_uzip_g+2,1,MPI_MIN,0,m_uiCommActive);
                t_uzip_g[1]=t_uzip_g[1]/(double) m_uiActiveNpes;
        #endif



        unsigned int regLev=0;
        ot::TreeNode blkNode;






        unsigned int centeredOffset=centered.getOffset();
        unsigned int backwardOffset=backward.getOffset();
        unsigned int forwardOffset=forward.getOffset();


        // all the 3 stencil directions should be in the same.
        assert(centered.getStencilDirection()==forward.getStencilDirection());
        assert(centered.getStencilDirection()==backward.getStencilDirection());
        double h=0.0;
        unsigned int lx,ly,lz,offset,paddWidth;
        #ifdef DEBUG_UNZIP_OP
        double d_min=-0.5;
        double d_max=0.5;
        std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){ return (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
        std::function<double(double,double,double)> dx_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(cos(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
        unsigned int  x,y,z,sz,regSz;

        for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++) {

            blkNode = m_uiLocalBlockList[blk].getBlockNode();
            regLev = m_uiLocalBlockList[blk].getRegularGridLev();
            lx=m_uiLocalBlockList[blk].getAllocationSzX();
            ly=m_uiLocalBlockList[blk].getAllocationSzY();
            lz=m_uiLocalBlockList[blk].getAllocationSzZ();
            offset=m_uiLocalBlockList[blk].getOffset();
            paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();
            //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));
            h = ((blkNode.maxX() - blkNode.minX())) / ((1u << (regLev - blkNode.getLevel())) * m_uiElementOrder);
            h = 1.0 / h;
            blkNpe_1D = m_uiElementOrder * (1u << (regLev - blkNode.getLevel())) + 1 + 2 * paddWidth;
            assert(blkNpe_1D > paddWidth);


            for(unsigned int k=0;k<(blkNpe_1D);k++)
                  for(unsigned int j=0;j<(blkNpe_1D);j++)
                    for(unsigned int i=0;i<(blkNpe_1D);i++)
                    {
                       assert(((blkNode.maxX()-blkNode.minX()))%((1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder)==0);
                       sz=((blkNode.maxX()-blkNode.minX()))/((1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                       regSz=1u<<(m_uiMaxDepth-regLev);

                       if( (i>paddWidth && i<(blkNpe_1D-paddWidth-1)) && (j>paddWidth && j<(blkNpe_1D-paddWidth-1)) && (k>paddWidth && k<(blkNpe_1D-paddWidth-1)))
                       {
                           x=blkNode.getX() + (i-paddWidth)*sz;
                           y=blkNode.getY() + (j-paddWidth)*sz;
                           z=blkNode.getZ() + (k-paddWidth)*sz;
                           if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                               std::cout<<" [internal node mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;

                       }

                        if((blkNode.getX()>=regSz) && (i>=0 && i<=(paddWidth)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX()-regSz + (i+paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [left ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }

                        if((blkNode.getY()>=regSz) && (j>=0 && j<=(paddWidth)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY()-regSz + (j+paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [down ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }


                        if((blkNode.getZ()>=regSz) && (k>=0 && k<=(paddWidth)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ()-regSz + (k+paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [back ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }


                        if((blkNode.maxX()+regSz<=m_uiMeshDomain_max) && (i>=(blkNpe_1D-paddWidth) && i<(blkNpe_1D)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX()-regSz + (i+paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [right ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }


                        if((blkNode.maxY()+regSz<=m_uiMeshDomain_max) && (j>=(blkNpe_1D-paddWidth) && j<(blkNpe_1D)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY() - regSz + (j+paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [up ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }

                        if((blkNode.maxZ()+regSz<=m_uiMeshDomain_max) && (k>=(blkNpe_1D-paddWidth) && k<(blkNpe_1D)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ() - regSz + (k+paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [front ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }




                    }



        }

        #endif

        #ifdef PROFILE_APPLY_STENCIL
                t1=std::chrono::high_resolution_clock::now();
        #endif



        if(centered.getStencilDirection()==StencilDirection::STENCIL_DIR_X)
        {

            for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
            {
                blkNode=m_uiLocalBlockList[blk].getBlockNode();
                regLev=m_uiLocalBlockList[blk].getRegularGridLev();

                lx=m_uiLocalBlockList[blk].getAllocationSzX();
                ly=m_uiLocalBlockList[blk].getAllocationSzY();
                lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                offset=m_uiLocalBlockList[blk].getOffset();
                paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();

                //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));
                h=((blkNode.maxX()-blkNode.minX()))/((double)(1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                h=1.0/h;
                blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*paddWidth;
                assert(blkNpe_1D>paddWidth);


                if(blkNode.minX()==m_uiMeshDomain_min)
                {

                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<2*paddWidth;i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<forward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=forward[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-forwardOffset]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=2*paddWidth;i<(blkNpe_1D-2*paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-centeredOffset]*h;

                            }

                    if(blkNode.maxX()==m_uiMeshDomain_max)
                    {
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                                for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth);i++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<backward.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-backwardOffset]*h;
                                }

                    }else
                    {
                        assert(blkNode.maxX()<m_uiMeshDomain_max);
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                                for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth);i++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<centered.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-centeredOffset]*h;

                                }
                    }



                }else if(blkNode.maxX()==m_uiMeshDomain_max)
                {

                    assert(blkNode.minX()>m_uiMeshDomain_min);
                    assert((blkNpe_1D-2*paddWidth));
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-2*paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-centeredOffset]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<backward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-backwardOffset]*h;
                            }


                }else
                {
                    assert(blkNode.minX()>m_uiMeshDomain_min && blkNode.maxX()<m_uiMeshDomain_max);
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying centered difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+(i+index-centeredOffset)]*h;

                                }

                            }

                }


            }

        }


        if(centered.getStencilDirection()==StencilDirection::STENCIL_DIR_Y)
        {

            for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
            {
                blkNode=m_uiLocalBlockList[blk].getBlockNode();
                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));

                lx=m_uiLocalBlockList[blk].getAllocationSzX();
                ly=m_uiLocalBlockList[blk].getAllocationSzY();
                lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                offset=m_uiLocalBlockList[blk].getOffset();
                paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();


                h=((blkNode.maxY()-blkNode.minY()))/((double)(1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                h=1.0/h;
                blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*paddWidth;
                assert(blkNpe_1D>paddWidth);


                if(blkNode.minY()==m_uiMeshDomain_min)
                {

                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=paddWidth;j<2*paddWidth;j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<forward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=forward[index]*unzipVec[offset+k*(ly*lx)+(j+index-forwardOffset)*(lx)+i]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=2*paddWidth;j<(blkNpe_1D-2*paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+i]*h;

                            }

                    if(blkNode.maxY()==m_uiMeshDomain_max)
                    {
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth);j++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<backward.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+(j+index-backwardOffset)*(lx)+i]*h;
                                }

                    }else
                    {
                        assert(blkNode.maxY()<m_uiMeshDomain_max);
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth);j++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<centered.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+i]*h;

                                }
                    }



                }else if(blkNode.maxY()==m_uiMeshDomain_max)
                {

                    assert(blkNode.minY()>m_uiMeshDomain_min);
                    assert((blkNpe_1D-2*paddWidth));
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-2*paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+i]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<backward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+(j+index-backwardOffset)*(lx)+i]*h;
                            }


                }else
                {
                    assert(blkNode.minY()>m_uiMeshDomain_min && blkNode.maxY()<m_uiMeshDomain_max);
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying centered difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+(i)]*h;

                                }

                            }

                }


            }

        }




        if(centered.getStencilDirection()==StencilDirection::STENCIL_DIR_Z)
        {


            for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
            {
                blkNode=m_uiLocalBlockList[blk].getBlockNode();
                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));

                lx=m_uiLocalBlockList[blk].getAllocationSzX();
                ly=m_uiLocalBlockList[blk].getAllocationSzY();
                lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                offset=m_uiLocalBlockList[blk].getOffset();
                paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();

                h=((blkNode.maxZ()-blkNode.minZ()))/((double)(1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                h=1.0/h;
                blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*paddWidth;
                assert(blkNpe_1D>paddWidth);


                if(blkNode.minZ()==m_uiMeshDomain_min)
                {

                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=paddWidth;k<2*paddWidth;k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<forward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=forward[index]*unzipVec[offset+(k+index-forwardOffset)*(ly*lx)+(j)*(lx)+i]*h;

                            }


                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=2*paddWidth;k<(blkNpe_1D-2*paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+i]*h;

                            }

                    if(blkNode.maxZ()==m_uiMeshDomain_max)
                    {
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth);k++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<backward.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+(k+index-backwardOffset)*(ly*lx)+(j)*(lx)+i]*h;
                                }

                    }else
                    {
                        assert(blkNode.maxZ()<m_uiMeshDomain_max);
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth);k++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<centered.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+i]*h;

                                }
                    }



                }else if(blkNode.maxZ()==m_uiMeshDomain_max)
                {

                    assert(blkNode.minZ()>m_uiMeshDomain_min);
                    assert((blkNpe_1D-2*paddWidth));
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=paddWidth;k<(blkNpe_1D-2*paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+i]*h;

                            }


                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<backward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+(k+index-backwardOffset)*(ly*lx)+(j)*(lx)+i]*h;
                            }


                }else
                {
                    assert(blkNode.minZ()>m_uiMeshDomain_min && blkNode.maxZ()<m_uiMeshDomain_max);
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying centered difference difference: "<<std::endl;
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+(i)]*h;

                                }

                            }

                }


            }


        }

        #ifdef PROFILE_APPLY_STENCIL
                t2=std::chrono::high_resolution_clock::now();
                t_stencil=std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

                par::Mpi_Reduce(&t_stencil,t_stencil_g,1,MPI_MIN,0,m_uiCommActive);
                par::Mpi_Reduce(&t_stencil,t_stencil_g+1,1,MPI_SUM,0,m_uiCommActive);
                par::Mpi_Reduce(&t_stencil,t_stencil_g+2,1,MPI_MIN,0,m_uiCommActive);
                t_stencil_g[1]=t_stencil_g[1]/(double) m_uiActiveNpes;

                t1=std::chrono::high_resolution_clock::now();
        #endif
                this->createVector(out);
                this->zip(&(*(unzipVec1.begin())),&(*(out.begin())));

        #ifdef PROFILE_APPLY_STENCIL
                t2=std::chrono::high_resolution_clock::now();
                t_zip=std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                par::Mpi_Reduce(&t_zip,t_zip_g,1,MPI_MIN,0,m_uiCommActive);
                par::Mpi_Reduce(&t_zip,t_zip_g+1,1,MPI_SUM,0,m_uiCommActive);
                par::Mpi_Reduce(&t_zip,t_zip_g+2,1,MPI_MIN,0,m_uiCommActive);
                t_zip_g[1]=t_zip_g[1]/(double) m_uiActiveNpes;

                if(!m_uiActiveRank)
                {
                    std::cout<<"unzip_max \t stencil_max \t zip_max "<<std::endl;
                    std::cout<<t_uzip_g[1]<<" \t "<<t_stencil_g[1]<<" \t "<<t_zip_g[1]<<std::endl;
                }

        #endif
        unzipVec1.clear();
        unzipVec.clear();

    }


    
    template <typename pKey, typename pNode>
    void Mesh::searchKeys(std::vector<pKey>& pKeys,std::vector<pNode>& pNodes)
    {


        assert(seq::test::isSorted(pNodes));

        std::vector<Key> pKeys_cpy;
        pKeys_cpy.resize(pKeys.size());

        for(unsigned int k=0;k<pKeys.size();k++)
        {
            pKeys_cpy[k]=pKeys[k];
            pKeys_cpy[k].addOwner(k);
            pKeys_cpy[k].setSearchResult(LOOK_UP_TABLE_DEFAULT);
        }

        SFC::seqSearch::SFC_treeSearch(&(*(pKeys_cpy.begin())),&(*(pNodes.begin())),0,pKeys_cpy.size(),0,pNodes.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);

        for(unsigned int k=0;k<pKeys_cpy.size();k++)
        {
            if((pKeys_cpy[k].getFlag() & OCT_FOUND)) {
                pKeys[(*(pKeys_cpy[k].getOwnerList()))[0]].setSearchResult(pKeys_cpy[k].getSearchResult());
                pKeys[k].setFlag((pKeys[k].getFlag() | OCT_FOUND));
            }

        }

        pKeys_cpy.clear();




    }


   
    template <typename T>
    bool Mesh::isReMeshUnzip(const T **unzippedVec,const unsigned int * varIds,const unsigned int numVars,std::function<double(double,double,double,double*)>wavelet_tol,double amr_coarse_fac, double coarsen_hx)
    {

        // This is the default isRMesh code that is used as refiment criteria. (if needed some complicated application specific refinement routine please have a look
        // at the waveletAMR.h(tcc) file. )

        // new wavelet code goes here.
        bool isMeshGlobalChanged = false;
        bool isMeshLocalChanged  = false;
        //std::cout<<"calling amr"<<std::endl;
        const bool includeBdy= true; // change this to false to exclude boundary from AMR. 
        std::vector<unsigned int> refine_flags;

        if(this->isActive())
        {
            RefElement* refEl = &m_uiRefEl;
            wavelet::WaveletEl* wrefEl = new wavelet::WaveletEl(refEl);
            
            const std::vector<ot::Block>& blkList = this->getLocalBlockList();
            const unsigned int eOrder = m_uiElementOrder;
            
            const unsigned int numLocalElements = m_uiNumLocalElements;
            
            refine_flags.clear();
            refine_flags.resize(numLocalElements,OCT_NO_CHANGE);
            
            std::vector<T> blkIn;
            std::vector<double> wCout;
            const ot::TreeNode* pNodes = m_uiAllElements.data();

            std::vector<double> eleWMax;
            eleWMax.resize(numLocalElements,0);

            const unsigned int eleOfst = m_uiElementLocalBegin;

            
            for(unsigned int blk=0; blk <blkList.size(); blk++)
            {   
                const unsigned int pw = blkList[blk].get1DPadWidth();
                if((eOrder>>1u) != pw)
                {
                    std::cout<<" padding width should be half the eleOrder for generic wavelet computations. "<<std::endl;
                    MPI_Abort(this->getMPICommunicator(),0);
                }

                const unsigned int nx = (2*eOrder+1);
                const unsigned int ny = (2*eOrder+1);
                const unsigned int nz = (2*eOrder+1);

                //std::cout<<"nx "<<nx<<std::endl;
        
                blkIn.resize(numVars*nx*ny*nz);
                const unsigned int isz[] = {nx,ny,nz};
                const unsigned int bflag = blkList[blk].getBlkNodeFlag();

                for(unsigned int ele =blkList[blk].getLocalElementBegin(); ele < blkList[blk].getLocalElementEnd(); ele++)
                {

                    const unsigned int pw = blkList[blk].get1DPadWidth();
                    const bool isBdyOct = this->isBoundaryOctant(ele);

                    const double oct_dx = (1u<<(m_uiMaxDepth-pNodes[ele].getLevel()))/(double(m_uiElementOrder));
                    Point oct_pt1 = Point(pNodes[ele].minX() , pNodes[ele].minY(), pNodes[ele].minZ());
                    Point oct_pt2 = Point(pNodes[ele].minX() + oct_dx , pNodes[ele].minY() + oct_dx, pNodes[ele].minZ() + oct_dx);
                    Point domain_pt1,domain_pt2,dx_domain;
                    this->octCoordToDomainCoord(oct_pt1,domain_pt1);
                    this->octCoordToDomainCoord(oct_pt2,domain_pt2);
                    dx_domain=domain_pt2-domain_pt1;
                    double hx[3] ={dx_domain.x(),dx_domain.y(),dx_domain.z()};
                    const double tol_ele = wavelet_tol(domain_pt1.x(),domain_pt1.y(),domain_pt1.z(),hx);

                    if(!includeBdy && isBdyOct)
                    {
                        // tol small enough to not refine but not to coarsen . 
                        eleWMax[ele - eleOfst] = amr_coarse_fac*tol_ele + 1e-8; 
                        continue;
                    }
                        

                    for(unsigned int v=0; v < numVars; v++)
                    {
                        const unsigned int vid = varIds[v];
                        this->getUnzipElementalNodalValues(unzippedVec[vid],blk, ele, blkIn.data() + v*(nx*ny*nz), true);
                    }

                    //eleWMax[ele - eleOfst]=wavelet::compute_element_wavelet(this,(const wavelet::WaveletEl*)&wrefEl,blkIn.data(),tol_ele,numVars,isBdyOct);
                    // compute the wavelet
                    {

                        double wMax=0.0;
                        
                        const unsigned int nx = (2*eOrder+1);
                        const unsigned int ny = (2*eOrder+1);
                        const unsigned int nz = (2*eOrder+1); 
                        assert(pw == (eOrder>>1u));

                        const unsigned int sz_per_dof = nx*ny*nz;
                        const unsigned int isz[] = {nx,ny,nz};
                        wCout.resize(sz_per_dof);

                        const unsigned int dof = numVars;
                        for(unsigned int v=0; v < dof; v++)
                        {
                            wrefEl->compute_wavelets_3D((double*)(blkIn.data()+ v*sz_per_dof),isz,wCout,isBdyOct);
                            const double l_max = (normL2(wCout.data(),wCout.size()))/sqrt(wCout.size());

                            if(wMax < l_max)
                                wMax = l_max;

                            // for early bail out. 
                            if(wMax > tol_ele)
                                break;

                        }

                        eleWMax[ele - eleOfst] = wMax;

                    }
                    
                    // if(isBdyOct)
                    //std::cout<<"ele :  "<<ele<<" eleWMax: "<<eleWMax[ele-eleOfst]<<std::endl;
                    
                }

            }

            // delete the wavelet reference element. 
            delete wrefEl;



            // mark elements for refinement first. 
            for(unsigned int ele = m_uiElementLocalBegin; ele < m_uiElementLocalEnd; ele++)
            {
                
                const double oct_dx = (1u<<(m_uiMaxDepth-pNodes[ele].getLevel()))/(double(m_uiElementOrder));
                Point oct_pt1 = Point(pNodes[ele].minX() , pNodes[ele].minY(), pNodes[ele].minZ());
                Point oct_pt2 = Point(pNodes[ele].minX() + oct_dx , pNodes[ele].minY() + oct_dx, pNodes[ele].minZ() + oct_dx);
                Point domain_pt1,domain_pt2,dx_domain;
                this->octCoordToDomainCoord(oct_pt1,domain_pt1);
                this->octCoordToDomainCoord(oct_pt2,domain_pt2);
                dx_domain=domain_pt2-domain_pt1;
                double hx[3] ={dx_domain.x(),dx_domain.y(),dx_domain.z()};
                const double tol_ele = wavelet_tol(domain_pt1.x(),domain_pt1.y(),domain_pt1.z(),hx);

                const double l_max = eleWMax[ele - eleOfst];

                if(l_max > tol_ele )
                {
                    refine_flags[(ele-eleOfst)] = OCT_SPLIT;
                    isMeshLocalChanged=true;

                }else if( l_max < amr_coarse_fac *tol_ele)
                {
                    refine_flags[ele-eleOfst] = OCT_COARSE;
                    isMeshLocalChanged=true;

                }else
                {
                    refine_flags[ele-eleOfst] = OCT_NO_CHANGE;
                } 

            }
            

            if(isMeshLocalChanged)
                isMeshLocalChanged = this->setMeshRefinementFlags(refine_flags);
                
            
        }

        
        //par::Mpi_Allreduce(&isMeshLocalChanged,&isMeshGlobalChanged,1,MPI_LOR,this->getMPIGlobalCommunicator());
        MPI_Allreduce(&isMeshLocalChanged,&isMeshGlobalChanged,1,MPI_CXX_BOOL,MPI_LOR,this->getMPIGlobalCommunicator());
        return isMeshGlobalChanged;

        // old remesh code hard coded oly for 4th order interp. for refine wavelets and 3rd order for coarsen wavelets. (not encouraged to use :) )
        #if 0
        bool isOctChange=false;
        if(m_uiIsActive)
        {
            // remove all the previously set falgs if there is any.  THIS will change the all flags to no CHANGE
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                m_uiAllElements[ele].setFlag(((OCT_NO_CHANGE<<NUM_LEVEL_BITS)|m_uiAllElements[ele].getLevel()));

            ot::TreeNode blkNode;
            unsigned int sz[3];
            double dh[3];
            unsigned int bflag,offset;
            unsigned int regLev;
            unsigned int eIndex[3];
            double *  waveletR = NULL;
            double *  waveletC = NULL;
            unsigned int num_wr =0 ,num_wc =0;

            double * wsIn = new double[m_uiNpE];
            double * wsOut = new double[m_uiNpE];
            double ** ws = new double*[2];
            ws[0] = wsIn;
            ws[1] = wsOut;

            // upper bound for the refine and coarsen wavelets.     
            waveletR = new double[64];
            num_wr = 64;
            
            waveletC = new double[64];
            num_wc = 64 ;

            const unsigned int paddWidth=3;
            unsigned int eleIndexMin=0,eleIndexMax=0;

            double l_inf;
            double x,y,z,tol;


            // first pass to identify the refined elements.
            for(unsigned blk=0;blk<m_uiLocalBlockList.size();blk++)
            {

                blkNode=m_uiLocalBlockList[blk].getBlockNode();

                sz[0]=m_uiLocalBlockList[blk].getAllocationSzX();
                sz[1]=m_uiLocalBlockList[blk].getAllocationSzY();
                sz[2]=m_uiLocalBlockList[blk].getAllocationSzZ();

                bflag=m_uiLocalBlockList[blk].getBlkNodeFlag();
                offset=m_uiLocalBlockList[blk].getOffset();

                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;

                //if(bflag!=0) continue;

                for(unsigned int ele=m_uiLocalBlockList[blk].getLocalElementBegin();ele<m_uiLocalBlockList[blk].getLocalElementEnd();ele++)
                {

                    if((m_uiAllElements[ele].getLevel()+MAXDEAPTH_LEVEL_DIFF+1)>=m_uiMaxDepth) continue;

                    x=m_uiAllElements[ele].getX();
                    y=m_uiAllElements[ele].getY();
                    z=m_uiAllElements[ele].getZ();
                    tol=wavelet_tol(x,y,z);

                    eIndex[0]=(m_uiAllElements[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    eIndex[1]=(m_uiAllElements[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    eIndex[2]=(m_uiAllElements[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                    if((bflag &(1u<<OCT_DIR_LEFT)) && eIndex[0]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_DOWN)) && eIndex[1]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_BACK)) && eIndex[2]==eleIndexMin)   continue;

                    if((bflag &(1u<<OCT_DIR_RIGHT)) && eIndex[0]==eleIndexMax)  continue;
                    if((bflag &(1u<<OCT_DIR_UP)) && eIndex[1]==eleIndexMax)     continue;
                    if((bflag &(1u<<OCT_DIR_FRONT)) && eIndex[2]==eleIndexMax)  continue;

                    for(unsigned int var=0;var<numVars;var++)
                    {

                         refine_wavelets(&unzippedVec[varIds[var]][offset],m_uiElementOrder,eIndex,paddWidth,sz,waveletR,num_wr,(double**)ws);
                        //  for(unsigned int k=0; k<4; k+=3)
                        //      for(unsigned int j=0; j<4; j+=3)
                        //       for(unsigned int i=0; i<4; i+=3)                                
                        //         waveletR[k*16 + j*4 + i] =0;

                         l_inf=normLInfty(waveletR,num_wr);
                         //l_inf = normL2(waveletR,num_wr)/num_wr;

                            // for(unsigned int k=1; k<3; k+=1)
                            //   for(unsigned int j=1; j<3; j+=1)
                            //    for(unsigned int i=1; i<3; i+=1)
                            //     std::cout<<"ref1: (i,j,k) : " << (i-1)<<" , "<<(j-1)<<" , "<<(k-1)<<": "<<waveletR[k*16 + j*4 + i]<<std::endl;
                            
                        
                       
                        // computeRefineWavelets(&unzippedVec[varIds[var]][offset],0,m_uiElementOrder,eIndex,paddWidth,sz,waveletR);
                        // l_inf=normLInfty(waveletR,NUM_REFINE_WAVELET_COEF);

                        //     for(unsigned int k=1; k<3; k+=1)
                        //       for(unsigned int j=1; j<3; j+=1)
                        //        for(unsigned int i=1; i<3; i+=1)
                        //         std::cout<<"ref2: (i,j,k) : " << (i-1)<<" , "<<(j-1)<<" , "<<(k-1)<<": "<<waveletR[(k-1)*4 + (j-1)*2 +i-1]<<std::endl;

                        if(l_inf>tol)
                        {
                            // for(unsigned int k=0;k<num_wr;k++)
                            //    std::cout<<"elem: "<<m_uiAllElements[ele]<<" wr["<<k<<"]: "<<waveletR[k]<<std::endl;
                            assert((m_uiAllElements[ele].getLevel()+MAXDEAPTH_LEVEL_DIFF+1)<m_uiMaxDepth);
                            //std::cout<<"rank: "<<m_uiActiveRank<<" element R: "<<m_uiAllElements[ele]<<" w_tol: "<<l_inf<<std::endl;
                            m_uiAllElements[ele].setFlag(((OCT_SPLIT<<NUM_LEVEL_BITS)|m_uiAllElements[ele].getLevel()));
                            assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT);
                            break; // no point checking for the other variables since this octree needs to be refined.
                        }

                    }



                }


            }

            //second pass to identify the coarsening elements.
            for(unsigned blk=0;blk<m_uiLocalBlockList.size();blk++)
            {

                blkNode=m_uiLocalBlockList[blk].getBlockNode();

                sz[0]=m_uiLocalBlockList[blk].getAllocationSzX();
                sz[1]=m_uiLocalBlockList[blk].getAllocationSzY();
                sz[2]=m_uiLocalBlockList[blk].getAllocationSzZ();

                bflag=m_uiLocalBlockList[blk].getBlkNodeFlag();
                offset=m_uiLocalBlockList[blk].getOffset();

                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;

                dh[0]=coarsen_hx*(m_uiLocalBlockList[blk].computeGridDx());
                dh[1]=coarsen_hx*(m_uiLocalBlockList[blk].computeGridDy());
                dh[2]=coarsen_hx*(m_uiLocalBlockList[blk].computeGridDz());


                if((eleIndexMax==0) || (bflag!=0)) continue; // this implies the blocks with only 1 child and boundary blocks.

                bool isEligibleCoarsen=true;
                bool isCoarsen=true;
                ot::TreeNode tmpOct;

                for(unsigned int ele=m_uiLocalBlockList[blk].getLocalElementBegin();ele<m_uiLocalBlockList[blk].getLocalElementEnd();ele+=NUM_CHILDREN)
                {

                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());

                    isEligibleCoarsen=true;
                    for(unsigned int child=0;child<NUM_CHILDREN;child++)
                    {
                        if((m_uiAllElements[ele+child].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                        {
                            isEligibleCoarsen=false;
                            break;
                        }

                    }

                    if((isEligibleCoarsen) && (m_uiAllElements[ele].getLevel()>1))
                    {
                        tmpOct=m_uiAllElements[ele].getParent();
                        x=tmpOct.getX() + (1u<<(m_uiMaxDepth-tmpOct.getLevel()-1));
                        y=tmpOct.getY() + (1u<<(m_uiMaxDepth-tmpOct.getLevel()-1));
                        z=tmpOct.getZ() + (1u<<(m_uiMaxDepth-tmpOct.getLevel()-1));
                        tol=wavelet_tol(x,y,z);
                        tmpOct=ot::TreeNode(tmpOct.getX(),tmpOct.getY(),tmpOct.getZ(),tmpOct.getLevel()+1,m_uiDim,m_uiMaxDepth);

                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            if(tmpOct==m_uiAllElements[ele+child])
                            {
                                eIndex[0]=(m_uiAllElements[ele+child].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                                eIndex[1]=(m_uiAllElements[ele+child].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                                eIndex[2]=(m_uiAllElements[ele+child].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                                break;
                            }

                        }

                        isCoarsen=true;

                        for(unsigned int var=0;var<numVars;var++)
                        {
                            coarsen_wavelets(&unzippedVec[varIds[var]][offset],m_uiElementOrder,eIndex,paddWidth,sz,waveletC,num_wc,(double**)ws);
                            //computeCoarsenWavelets(unzippedVec[varIds[var]],offset,m_uiElementOrder,eIndex,paddWidth,sz,waveletC);
                            l_inf=normLInfty(waveletC,NUM_COARSE_WAVELET_COEF);
                            //l_inf=normLInfty(waveletC,num_wc);
                            //l_inf = normL2(waveletC,num_wc)/num_wc;
                            if(l_inf>amr_coarse_fac*tol)
                            {
                                isCoarsen=false;
                                break;
                            }

                        }


                        if(isCoarsen)
                        {

                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            {
                                m_uiAllElements[ele+child].setFlag(((OCT_COARSE<<NUM_LEVEL_BITS)|m_uiAllElements[ele].getLevel()));
                                assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE);
                                //std::cout<<"rank: "<<m_uiActiveRank<<" element C: "<<m_uiAllElements[ele]<<" is coarsening "<<l_inf<<std::endl;
                            }

                        }


                    }

                }


            }

            delete [] waveletR;
            delete [] waveletC;
            delete [] wsIn;
            delete [] wsOut;
            delete [] ws;
            
            isOctChange=false;
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)!=OCT_NO_CHANGE)//if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT) // trigger remesh only when some refinement occurs (laid back remesh :)  ) //if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)!=OCT_NO_CHANGE)
                {
                    isOctChange=true;
                    break;
                }

        }

        bool isOctChanged_g;
        MPI_Allreduce(&isOctChange,&isOctChanged_g,1,MPI_CXX_BOOL,MPI_LOR,m_uiCommGlobal);
        //if(!m_uiGlobalRank) std::cout<<"is oct changed: "<<isOctChanged_g<<std::endl;
        return isOctChanged_g;
        #endif



    }

    template<typename T>
    void Mesh::getElementNodalValues(const T* vec,T* nodalValues,unsigned int elementID, bool isDGVec ) const
    {

        if(!m_uiIsActive) return;

        // handles the element get nodal values if the vec is an element DG vector. 
        if(isDGVec)
        {
            for(unsigned int node =0; node < m_uiNpE; node++)
                nodalValues[node] = vec[elementID*m_uiNpE +  node];

            return;
        }


        #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_nodalval.start();
        #endif
        std::vector<T> edgeInpIn;
        std::vector<T> edgeInpOut;

        std::vector<T> faceInpIn;
        std::vector<T> faceInpOut;

        unsigned int cnum;
        bool isHanging;

        std::vector<unsigned int > edgeIndex;
        std::vector<unsigned int > faceIndex;


        bool nodeStatus[OCT_DIR_TOTAL];
        for(unsigned int w=0;w<OCT_DIR_TOTAL;w++)
            nodeStatus[w]=false;

        edgeInpIn.resize((m_uiElementOrder+1));
        edgeInpOut.resize((m_uiElementOrder+1));

        faceInpIn.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));
        faceInpOut.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));

        for(unsigned int k=1;k<(m_uiElementOrder);k++)
            for(unsigned int j=1;j<(m_uiElementOrder);j++)
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                {
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                    assert(!(this->isNodeHanging(elementID,i,j,k))); // internal nodes cannot be hangging.
                }
        nodeStatus[OCT_DIR_INTERNAL]=true;

        // face interpolations
        // face : OCT_DIR_LEFT (1)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_LEFT,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_LEFT, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=faceInpOut[k*(m_uiElementOrder+1)+j];

            nodeStatus[OCT_DIR_LEFT_DOWN]=true;
            nodeStatus[OCT_DIR_LEFT_UP]=true;
            nodeStatus[OCT_DIR_LEFT_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_FRONT]=true;


            nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;



        }else
        {

            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int j=1;j<m_uiElementOrder;j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]];
        }


        // face : OCT_DIR_RIGHT (2)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_RIGHT,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_RIGHT, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=faceInpOut[k*(m_uiElementOrder+1)+j];

            nodeStatus[OCT_DIR_RIGHT_DOWN]=true;
            nodeStatus[OCT_DIR_RIGHT_UP]=true;
            nodeStatus[OCT_DIR_RIGHT_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_FRONT]=true;


            nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;


        }else
        {
            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int j=1;j<m_uiElementOrder;j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]];
        }



        // face : OCT_DIR_DOWN (3)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_DOWN,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_DOWN, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=faceInpOut[k*(m_uiElementOrder+1)+i];


            nodeStatus[OCT_DIR_RIGHT_DOWN]=true;
            nodeStatus[OCT_DIR_LEFT_DOWN]=true;
            nodeStatus[OCT_DIR_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_DOWN_FRONT]=true;

            nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;


        }else
        {

            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]];
        }


        // face : OCT_DIR_UP (4)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_UP,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_UP, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=faceInpOut[k*(m_uiElementOrder+1)+i];

            nodeStatus[OCT_DIR_RIGHT_UP]=true;
            nodeStatus[OCT_DIR_LEFT_UP]=true;
            nodeStatus[OCT_DIR_UP_BACK]=true;
            nodeStatus[OCT_DIR_UP_FRONT]=true;


            nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;

        }else
        {


            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]];
        }



        // face : OCT_DIR_BACK (5)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_BACK,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_BACK, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=faceInpOut[j*(m_uiElementOrder+1)+i];

            nodeStatus[OCT_DIR_LEFT_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_BACK]=true;
            nodeStatus[OCT_DIR_UP_BACK]=true;
            nodeStatus[OCT_DIR_DOWN_BACK]=true;

            nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;


        }else
        {

            for(unsigned int j=1;j<m_uiElementOrder;j++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
        }


        // face : OCT_DIR_FRONT (6)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_FRONT,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_FRONT, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=faceInpOut[j*(m_uiElementOrder+1)+i];


            nodeStatus[OCT_DIR_LEFT_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_FRONT]=true;
            nodeStatus[OCT_DIR_UP_FRONT]=true;
            nodeStatus[OCT_DIR_DOWN_FRONT]=true;

            nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;


        }else
        {
            for(unsigned int j=1;j<m_uiElementOrder;j++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
        }



        // edge: OCT_DIR_LEFT_DOWN (1)

        if((!nodeStatus[OCT_DIR_LEFT_DOWN]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_DOWN,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=edgeInpOut[k];

                nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]];
            }
        }



        // edge: OCT_DIR_LEFT_UP (2)

        if((!nodeStatus[OCT_DIR_LEFT_UP]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_UP,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_UP,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=edgeInpOut[k];

                nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
                nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]];
            }
        }


        // edge: OCT_DIR_LEFT_BACK (3)

        if((!nodeStatus[OCT_DIR_LEFT_BACK]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=edgeInpOut[j];

                nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]];
            }
        }



        // edge: OCT_DIR_LEFT_FRONT(4)

        if((!nodeStatus[OCT_DIR_LEFT_FRONT]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=edgeInpOut[j];

                nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;
                nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]];
            }
        }



        // edge: OCT_DIR_RIGHT_DOWN (5)

        if((!nodeStatus[OCT_DIR_RIGHT_DOWN]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_DOWN,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[k];

                nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }



        // edge: OCT_DIR_RIGHT_UP (6)

        if((!nodeStatus[OCT_DIR_RIGHT_UP]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_UP,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_UP,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[k];

                nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }


        // edge: OCT_DIR_RIGHT_BACK (7)

        if((!nodeStatus[OCT_DIR_RIGHT_BACK]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[j];

                nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }



        // edge: OCT_DIR_RIGHT_FRONT(8)

        if((!nodeStatus[OCT_DIR_RIGHT_FRONT]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[j];

                nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }





        // edge: OCT_DIR_DOWN_BACK (9)

        if((!nodeStatus[OCT_DIR_DOWN_BACK]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_DOWN_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_DOWN,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;

            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]];
            }
        }


        // edge: OCT_DIR_DOWN_FRONT (10)

        if(!nodeStatus[OCT_DIR_DOWN_FRONT])
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_DOWN_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;
                nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;

            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]];
            }
        }


        // edge: OCT_DIR_UP_BACK (11)

        if(!nodeStatus[OCT_DIR_UP_BACK])
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_UP_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_UP,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;

            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]];
            }
        }



        // edge: OCT_DIR_UP_FRONT (12)

        if(!nodeStatus[OCT_DIR_UP_FRONT])
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_UP_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_UP,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;


            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]];
            }
        }


        //node: OCT_DIR_LEFT_DOWN_BACK
        if((!(this->isNodeHanging(elementID,0,0,0))) || (!nodeStatus[OCT_DIR_LEFT_DOWN_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_DOWN_BACK
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,0,0)) || (!nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]];

        //node: OCT_DIR_LEFT_UP_BACK
        if(!(this->isNodeHanging(elementID,0,m_uiElementOrder,0)) || (!nodeStatus[OCT_DIR_LEFT_UP_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_UP_BACK
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,m_uiElementOrder,0)) || (!nodeStatus[OCT_DIR_RIGHT_UP_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]];


        //node: OCT_DIR_LEFT_DOWN_FRONT
        if(!(this->isNodeHanging(elementID,0,0,m_uiElementOrder))|| (!nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_DOWN_FRONT
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,0,m_uiElementOrder))|| (!nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]];

        //node: OCT_DIR_LEFT_UP_FRONT
        if(!(this->isNodeHanging(elementID,0,m_uiElementOrder,m_uiElementOrder)) || (!nodeStatus[OCT_DIR_LEFT_UP_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_UP_FRONT
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,m_uiElementOrder,m_uiElementOrder)) || (!nodeStatus[OCT_DIR_RIGHT_UP_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]];

        #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_nodalval.stop();
        #endif

    }

    template<typename T>
    void Mesh::computeElementalContribution(const T* in,T* out,unsigned int elementID) const
    {
        if(!m_uiIsActive) return;

        const unsigned int eleOrder = m_uiElementOrder;
        const unsigned int npe_1d = eleOrder + 1;
        const unsigned int npe_2d = (eleOrder + 1) * (eleOrder + 1);
        const unsigned int nPe = (eleOrder + 1) * (eleOrder + 1) * (eleOrder + 1);

        //@todo later move this to outer allocation and reuse the memeory. 
        double* qMat = new double[nPe*nPe];
        double* qTIn = new double[nPe];

        this->getElementQMat(elementID,qMat,true);
        
        for (unsigned int i=0;i<nPe; i++)
        {
            qTIn[i] = 0;

            for(unsigned int j=0;j<nPe;j++)
            {
                qTIn[i] += qMat[j*nPe + i] *in[j]; // note the transpose. 
            }
        }

        for (unsigned int i=0;i<nPe; i++)
            out[m_uiE2NMapping_CG[elementID*nPe + i]] += qTIn[i];


        delete [] qMat;
        delete [] qTIn;

        
        return;


    }

    template<typename T>
    void Mesh::interGridTransfer(std::vector<T> & vec,const ot::Mesh* pMesh, INTERGRID_TRANSFER_MODE mode)
    {

        std::vector<T> tvec;
        pMesh->createVector<T>(tvec,0);

        this->interGridTransfer(vec.data(),tvec.data(),pMesh,mode,1);

        std::swap(vec,tvec);
        tvec.clear();
        return;
       
    }


    template<typename T>
    void Mesh::interGridTransfer(T*& vec,const ot::Mesh* pMesh, INTERGRID_TRANSFER_MODE mode, unsigned int dof)
    {

        T* tVec = pMesh->createCGVector<T>(0,dof);
        this->interGridTransfer(vec,tVec,pMesh,mode,dof);

        std::swap(vec,tVec);
        delete [] tVec;
        return;
    }

    template<typename T>
    void Mesh::interGridTransfer(T* vecIn, T* vecOut, const ot::Mesh* pMesh , INTERGRID_TRANSFER_MODE mode,unsigned int dof)
    {

        MPI_Comm comm=m_uiCommGlobal;
        int rank,npes;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        std::vector<unsigned int> sendC;
        std::vector<unsigned int> recvC;

        std::vector<unsigned int> sendOfst;
        std::vector<unsigned int> recvOfst;

        sendC.resize(npes);
        recvC.resize(npes);
        sendOfst.resize(npes);
        recvOfst.resize(npes);

        
        
        this->interGridTransferSendRecvCompute(pMesh);
        const unsigned int cg_sz_old = m_uiNumActualNodes;
        const unsigned int cg_sz_new = pMesh->getDegOfFreedom();
        const ot::TreeNode* m2prime = m_uiM2Prime.data();

        // scale the elemental counts by m_uiNpE;
        for(unsigned int p=0; p < npes; p++)
        {
            sendC[p] = m_uiIGTSendC[p] * m_uiNpE;
            recvC[p] = m_uiIGTRecvC[p] * m_uiNpE;

            sendOfst[p] = m_uiIGTSendOfst[p] *  m_uiNpE;
            recvOfst[p] = m_uiIGTRecvOfst[p] *  m_uiNpE;
        }


        std::vector<T> wVec; // dg of m2prime;
        std::vector<T> nodalVals;
        nodalVals.resize(m_uiNpE);
    
        unsigned int cnum;
        bool isHanging;
    
        std::vector<double> vallchildren;
        std::vector<T> wVec_m2;
    
        vallchildren.resize((2*m_uiElementOrder + 1)*(2*m_uiElementOrder + 1)*(2*m_uiElementOrder + 1));
        wVec_m2.resize(recvOfst[npes-1]+recvC[npes-1]);

        for(unsigned int var=0; var < dof; var++)
        {
            T* vec = vecIn  + (var * cg_sz_old);
            T* out = vecOut + (var * cg_sz_new);

            if(m_uiIsActive)
            {   

                const unsigned int npes1 = this->getMPICommSize();
                const unsigned int rank1 = this->getMPIRank(); 

                const unsigned int numM2PrimeElems = m_uiM2Prime.size();
                wVec.resize(numM2PrimeElems*m_uiNpE);
                

                //std::cout<<"rank1: "<<rank1<<" m2prime: "<<m2prime.size()<<std::endl;
                
                unsigned int m2primeCount=0;
                for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                {

                    //std::cout<<" m2primeCount: "<<m2primeCount<<" wvec offset : "<<m2primeCount*m_uiNpE<< " bound:" << (m2primeCount+1)*m_uiNpE <<" wvec size : "<<wVec.size()<<std::endl;

                    if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    {
                        this->getElementNodalValues(vec,&(*(nodalVals.begin())),ele);
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            cnum=m2prime[m2primeCount+child].getMortonIndex();
                            this->parent2ChildInterpolation(&(*(nodalVals.begin())),&(*(wVec.begin()+(m2primeCount+child)*m_uiNpE)),cnum,3);
                        }

                        m2primeCount+=NUM_CHILDREN;

                    }
                    else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                    {
                        assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                        
                        if(mode == INTERGRID_TRANSFER_MODE::P2CT)
                        {

                            const unsigned int p1d = 2*m_uiElementOrder+1;
                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            {
                                
                                this->getElementNodalValues(vec,nodalVals.data(),ele+child);
                                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                                    {
                                        cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                        const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                        const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                        const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;

                                        vallchildren[kkz*p1d*p1d  + jjy*p1d + iix] = nodalVals[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]; //vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];

                                    }

                            }

                            m_uiRefEl.I3D_Children2Parent(vallchildren.data(),&wVec[m2primeCount*m_uiNpE]);


                        }else
                        {
                            assert(mode == INTERGRID_TRANSFER_MODE::INJECTION);
                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            {
                                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                                    {

                                        isHanging=this->isNodeHanging((ele+child),i,j,k);
                                        if(isHanging)
                                        {
                                            wVec[m2primeCount*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];

                                        }else
                                        {
                                            
                                            cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                            {
                                                wVec[ m2primeCount*m_uiNpE +  (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                                                
                                            }

                                        }

                                    }

                            }

                        }

                        ele+=(NUM_CHILDREN-1);
                        m2primeCount+=1;

                    }else
                    {
                        assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                        
                        this->getElementNodalValues(vec,&(*(wVec.begin()+(m2primeCount*m_uiNpE))),ele);
                        m2primeCount+=1;

                    }

                }


                if(npes1==1 && pMesh->isActive() && pMesh->getMPICommSize()==1)
                {

                    // sequential case.

                    if(numM2PrimeElems !=pMesh->getNumLocalMeshElements())
                    {
                        std::cout<<" seq::[Inter-grid Transfer error ]: Recvn DG elements: "<<numM2PrimeElems<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;
                        MPI_Abort(comm,0);
                    }
                        

                    const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

                    const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                    const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                    const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
                    const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

                    unsigned int lookUp;
                    const unsigned int eleOrder=pMesh->getElementOrder();

                    for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                    {

                        //std::cout<<"ele: "<<ele<<"data copied "<<std::endl;
                        for(unsigned int k=0;k<eleOrder+1;k++)
                            for(unsigned int j=0;j<eleOrder+1;j++)
                                for(unsigned int i=0;i<eleOrder+1;i++)
                                {
                                    if(!(pMesh->isNodeHanging(ele,i,j,k)))
                                    {
                                        lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                        if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                            out[lookUp] = wVec[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                            
                                    }

                                }

                    }

                    continue;
                
                }


            }
            
            par::Mpi_Alltoallv_sparse(&(*(wVec.begin())), (int*)sendC.data(), (int*)sendOfst.data(), &(*(wVec_m2.begin())), (int*)recvC.data(), (int*)recvOfst.data(), comm);
            if(pMesh->isActive())
            {

                const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

                const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
                const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

                unsigned int lookUp;
                const unsigned int eleOrder=pMesh->getElementOrder();

                for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                {
                    for(unsigned int k=0;k<eleOrder+1;k++)
                        for(unsigned int j=0;j<eleOrder+1;j++)
                            for(unsigned int i=0;i<eleOrder+1;i++)
                            {
                                if(!(pMesh->isNodeHanging(ele,i,j,k)))
                                {
                                    lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                    if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                        out[lookUp]=wVec_m2[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                }
                            }

                }
            
            }
            

        }
        
        return;

    }

    template <typename T>
    void Mesh::interGridTransfer_DG(T* vecIn, T* vecOut, const ot::Mesh* pMesh, unsigned int dof)
    {
        
        // Note that this is the intergrid transfer for the DG representation of the vector, 
        // In DG / octant local representation there is no hanging nodes, each octant has it's own shared nodes. 

        MPI_Comm comm=m_uiCommGlobal;
        int rank,npes;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        std::vector<unsigned int> sendC;
        std::vector<unsigned int> recvC;

        std::vector<unsigned int> sendOfst;
        std::vector<unsigned int> recvOfst;

        sendC.resize(npes);
        recvC.resize(npes);
        sendOfst.resize(npes);
        recvOfst.resize(npes);

        std::vector<T> wVec; // dg of m2prime;
        
        this->interGridTransferSendRecvCompute(pMesh);

        const unsigned int dg_sz_old = getDegOfFreedomDG();
        const unsigned int dg_sz_new = pMesh->getDegOfFreedomDG();

        const ot::TreeNode* m2prime = m_uiM2Prime.data();
        // scale the elemental counts by m_uiNpE;
        for(unsigned int p=0; p < npes; p++)
        {
            sendC[p] = m_uiIGTSendC[p] * m_uiNpE;
            recvC[p] = m_uiIGTRecvC[p] * m_uiNpE;

            sendOfst[p] = m_uiIGTSendOfst[p] *  m_uiNpE;
            recvOfst[p] = m_uiIGTRecvOfst[p] *  m_uiNpE;
        }

        std::vector<T> wVec_m2;
        wVec_m2.resize(recvOfst[npes-1]+recvC[npes-1]);

        std::vector<T> nodalVals;
        nodalVals.resize(m_uiNpE);

        unsigned int cnum;
        bool isHanging;

        std::vector<double> vallchildren;
        vallchildren.resize((2*m_uiElementOrder + 1)*(2*m_uiElementOrder + 1)*(2*m_uiElementOrder + 1));


        for(unsigned int v=0; v < dof; v++)
        {
            T* vec = vecIn  + v * dg_sz_old;
            T* out = vecOut + v * dg_sz_new;

            if(m_uiIsActive)
            {   

                const unsigned int npes1 = this->getMPICommSize();
                const unsigned int rank1 = this->getMPIRank(); 

                const unsigned int numM2PrimeElems = m_uiM2Prime.size();
                wVec.resize(numM2PrimeElems*m_uiNpE);

                //std::cout<<"rank1: "<<rank1<<" m2prime: "<<m2prime.size()<<std::endl;
                
                unsigned int m2primeCount=0;
                for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                {
                    if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    {
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            cnum=m2prime[m2primeCount+child].getMortonIndex();
                            this->parent2ChildInterpolation(vec + ele*m_uiNpE, &(*(wVec.begin()+(m2primeCount+child)*m_uiNpE)),cnum,3);
                        }

                        m2primeCount+=NUM_CHILDREN;

                    }
                    else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                    {
                        assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                        
                        // for DG we use only one mode for the coarsening, 
                        // pure injection
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                            for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                            for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                            {
                                    cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                    const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                    const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                    const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                    //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                    if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                    {
                                        wVec[ m2primeCount*m_uiNpE +  (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = vec[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                    }
                                
                                }

                        }

                        ele+=(NUM_CHILDREN-1);
                        m2primeCount+=1;

                    }else
                    {
                        assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                        
                        for(unsigned int node=0; node < m_uiNpE; node ++)
                        wVec[(m2primeCount*m_uiNpE) + node] = vec[ele*m_uiNpE + node];
                        
                        m2primeCount+=1;

                    }

                }


                if(npes1==1 && pMesh->isActive() && pMesh->getMPICommSize()==1)
                {

                    // sequential case.

                    if(numM2PrimeElems !=pMesh->getNumLocalMeshElements())
                    {
                        std::cout<<" seq::[Inter-grid Transfer error ]: Recvn DG elements: "<<numM2PrimeElems<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;
                        MPI_Abort(comm,0);
                    }
                        

                    const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

                    const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                    const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                    const unsigned int m2LocalNodeBegin = m2LocalElemBegin * m_uiNpE ;
                    const unsigned int m2LocalNodeEnd   = m2LocalElemEnd   * m_uiNpE ;

                    unsigned int lookUp;
                    const unsigned int eleOrder=pMesh->getElementOrder();

                    for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                    {
                        for(unsigned int k=0;k<eleOrder+1;k++)
                            for(unsigned int j=0;j<eleOrder+1;j++)
                                for(unsigned int i=0;i<eleOrder+1;i++)
                                {
                                    lookUp= ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i;
                                    if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                        out[lookUp]=wVec[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                }
                    }

                    continue ;
                
                }


            }

            par::Mpi_Alltoallv_sparse(&(*(wVec.begin())), (int*)sendC.data(), (int*)sendOfst.data(), &(*(wVec_m2.begin())), (int*)recvC.data(), (int*)recvOfst.data(), comm);
        
            if(pMesh->isActive())
            {

                const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                const unsigned int m2LocalNodeBegin = m2LocalElemBegin*m_uiNpE ;
                const unsigned int m2LocalNodeEnd   = m2LocalElemEnd*m_uiNpE ;

                unsigned int lookUp;
                const unsigned int eleOrder=pMesh->getElementOrder();

                for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                {
                    for(unsigned int k=0;k<eleOrder+1;k++)
                        for(unsigned int j=0;j<eleOrder+1;j++)
                            for(unsigned int i=0;i<eleOrder+1;i++)
                            {
                                lookUp = ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i;
                                
                                if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                out[lookUp]=wVec_m2[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                
                            }

                }
            

            }

        }

        return;

    }


    template<typename T>
    void Mesh::interGridTransferCellVec(T* vecIn, T* vecOut, const ot::Mesh* pMesh,unsigned int dof,INTERGRID_TRANSFER_MODE mode)
    {

        MPI_Comm comm=m_uiCommGlobal;
        int rank,npes;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        std::vector<T> wVec; // dg of m2prime;
        this->interGridTransferSendRecvCompute(pMesh);

        // currently hard coded to the cell vec copy. 
        assert(mode==INTERGRID_TRANSFER_MODE::CELLVEC_CPY);

        const unsigned int cell_sz_old = m_uiAllElements.size();
        const unsigned int cell_sz_new = pMesh->getAllElements().size();

        const ot::TreeNode* m2prime = m_uiM2Prime.data();
        
        const unsigned int * sendC = m_uiIGTSendC.data();
        const unsigned int * recvC = m_uiIGTRecvC.data();
        const unsigned int * sendOfst = m_uiIGTSendOfst.data();
        const unsigned int * recvOfst = m_uiIGTRecvOfst.data();

        std::vector<T> wVec_m2;
        wVec_m2.resize(recvOfst[npes-1]+recvC[npes-1]);

        for(unsigned int v=0; v < dof; v++)
        {
            T* vec = vecIn  + v * cell_sz_old;
            T* out = vecOut + v * cell_sz_new;

            if(m_uiIsActive)
            {   

                const unsigned int npes1 = this->getMPICommSize();
                const unsigned int rank1 = this->getMPIRank(); 

                const unsigned int numM2PrimeElems = m_uiM2Prime.size();
                wVec.resize(numM2PrimeElems);

                //std::cout<<"rank1: "<<rank1<<" m2prime: "<<m2prime.size()<<std::endl;
                
                unsigned int m2primeCount=0;
                for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                {
                    if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    {
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            wVec[m2primeCount+child] = vec[ele];

                        m2primeCount+=NUM_CHILDREN;

                    }
                    else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                    {
                        assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                        // check if the cell vector child cells agrees with the value. 
                        assert(vec[ele]==vec[ele+NUM_CHILDREN-1]);
                        wVec[m2primeCount] = vec[ele];
                        ele+=(NUM_CHILDREN-1);
                        m2primeCount+=1;

                    }else
                    {
                        assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                        wVec[m2primeCount] = vec[ele];
                        m2primeCount+=1;

                    }

                }


                if(npes1==1 && pMesh->isActive() && pMesh->getMPICommSize()==1)
                {

                    // sequential case.

                    if(numM2PrimeElems !=pMesh->getNumLocalMeshElements())
                    {
                        std::cout<<" seq::[Inter-grid Transfer error ]: Recvn DG elements: "<<numM2PrimeElems<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;
                        MPI_Abort(comm,0);
                    }
                        

                    const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

                    const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                    const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                    const unsigned int m2LocalNodeBegin = m2LocalElemBegin * m_uiNpE ;
                    const unsigned int m2LocalNodeEnd   = m2LocalElemEnd   * m_uiNpE ;

                    unsigned int lookUp;
                    const unsigned int eleOrder=pMesh->getElementOrder();

                    for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                    {
                        out[ele] = wVec[(ele-m2LocalElemBegin)];
                    }

                    continue ;
                
                }


            }

            par::Mpi_Alltoallv_sparse(&(*(wVec.begin())), (int*)sendC, (int*)sendOfst, &(*(wVec_m2.begin())), (int*)recvC, (int*)recvOfst, comm);
        
            if(pMesh->isActive())
            {

                const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                const unsigned int m2LocalNodeBegin = m2LocalElemBegin*m_uiNpE ;
                const unsigned int m2LocalNodeEnd   = m2LocalElemEnd*m_uiNpE ;

                unsigned int lookUp;
                const unsigned int eleOrder=pMesh->getElementOrder();

                for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                    out[ele] = wVec_m2[(ele-m2LocalElemBegin)];

            }

        }

        return;






    }

    template<typename T>
    void Mesh::zip(const T* unzippedVec, T* zippedVec)
    {

        if(!m_uiIsActive) return;

        ot::TreeNode blkNode;
        unsigned int ei,ej,ek;
        unsigned int regLev;
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        unsigned int lx,ly,lz,offset,paddWidth;

        for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
        {
            blkNode=m_uiLocalBlockList[blk].getBlockNode();
            regLev=m_uiLocalBlockList[blk].getRegularGridLev();

            lx=m_uiLocalBlockList[blk].getAllocationSzX();
            ly=m_uiLocalBlockList[blk].getAllocationSzY();
            lz=m_uiLocalBlockList[blk].getAllocationSzZ();
            offset=m_uiLocalBlockList[blk].getOffset();
            paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();


            for(unsigned int elem=m_uiLocalBlockList[blk].getLocalElementBegin();elem<m_uiLocalBlockList[blk].getLocalElementEnd();elem++)
            {
                ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);


                assert(pNodes[elem].getLevel()==regLev); // this is enforced by block construction


                // (1). local nodes copy. Not need to interpolate or inject values. By block construction local octants in the block has is the same level as regular grid.
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            if((m_uiE2NMapping_DG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]/m_uiNpE)==elem)
                                zippedVec[m_uiE2NMapping_CG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]]=unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)];
                        }

            }
        }

    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_DOWN_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_DOWN-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        const unsigned int ei=0;
        const unsigned int ej=0;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_DOWN;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_RIGHT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);
        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=3;
        const unsigned int cnum2=7;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;

        unsigned int bflag=blk.getBlkNodeFlag();


        //std::cout<<" lookup : "<<pNodes[lookUp]<<" blkNode: "<<blk.getBlockNode()<<std::endl;

        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                // note this might not be the cnum1 cnum2.
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;
                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }

                edgeCount+=2;


            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];

                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_UP_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();


        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_UP-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        const unsigned int ei=0;
        const unsigned int ej=blkElem_1D-1;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_UP;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_RIGHT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=0;
        const unsigned int je=paddWidth+1;
        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=1;
        const unsigned int cnum2=5;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;


        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;

                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }

                }

                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }



                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];

                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_LEFT_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=0;
        unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_RIGHT;
        const unsigned int dir5=OCT_DIR_FRONT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);




        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=5;
        const unsigned int cnum2=7;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;


        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;
                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;


            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_BACK_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=0;
        unsigned int ej=0;
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_RIGHT;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(paddWidth+1);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=1;
        const unsigned int cnum2=3;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;

                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }



                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }
                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_DOWN_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.



        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_DOWN-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=(blkElem_1D-1);
        const unsigned int ej=0;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_DOWN;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_UP;
        ot::TreeNode blkNode=blk.getBlockNode();
        unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);
        const unsigned int ib=0;
        const unsigned int ie=paddWidth+1;

        const unsigned int cnum1=2;
        const unsigned int cnum2=6;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;

        unsigned int bflag=blk.getBlkNodeFlag();


        //std::cout<<" lookup : "<<pNodes[lookUp]<<" blkNode: "<<blk.getBlockNode()<<std::endl;

        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"RIGHT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;
                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"RIGHT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_UP_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_UP-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=blkElem_1D-1;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_UP;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_UP;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=0;
        const unsigned int je=paddWidth+1;
        const unsigned int ib=0;
        const unsigned int ie=paddWidth+1;

        const unsigned int cnum1=0;
        const unsigned int cnum2=4;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"RIGHT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;
                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }

                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"RIGHT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_RIGHT_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=blkElem_1D-1;
        unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_FRONT;
        const unsigned int dir5=OCT_DIR_FRONT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth+1;

        const unsigned int cnum1=4;
        const unsigned int cnum2=6;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;
                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_BACK_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=(blkElem_1D-1);
        unsigned int ej=0;
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_RIGHT;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(paddWidth+1);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth+1;

        const unsigned int cnum1=0;
        const unsigned int cnum2=2;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;
                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                           else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_DOWN_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_DOWN;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_LEFT;
        const unsigned int dir5=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=6;
        const unsigned int cnum2=7;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"DOWN_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"DOWN_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_DOWN_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_DOWN;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth+1;

        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=2;
        const unsigned int cnum2=3;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"DOWN_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"DOWN_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.

        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_UP_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        unsigned int ei=0;
        const unsigned int ej=(blkElem_1D-1);
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_UP;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_LEFT;
        const unsigned int dir5=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=paddWidth+1;

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=4;
        const unsigned int cnum2=5;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"UP_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"UP_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_UP_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        unsigned int ei=0;
        const unsigned int ej=(blkElem_1D-1);
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_UP;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth+1;

        const unsigned int jb=0;
        const unsigned int je=paddWidth+1;

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=0;
        const unsigned int cnum2=1;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"UP_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
                if(!eleDGValid[lookUp])
                {
                    getElementNodalValues(zippedVec,lookUpVec,lookUp);
                    eleDGValid[lookUp]=true;    
                }
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }

                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"DOWN_UP_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_LEFT_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_DOWN_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }
                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_DOWN_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth+1;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_UP_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=0;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=paddWidth+1;

        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_UP_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=paddWidth+1;

        const unsigned int ib=0;
        const unsigned int ie=paddWidth+1;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_LEFT_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_DOWN_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth+1;

        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_DOWN_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=0;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth + 1;

        const unsigned int jb=(m_uiElementOrder-paddWidth);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth +1;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_UP_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=0;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth +1 ;

        const unsigned int jb=0;
        const unsigned int je=paddWidth +1;

        const unsigned int ib=(m_uiElementOrder-paddWidth);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_UP_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth+1;

        const unsigned int jb=0;
        const unsigned int je=paddWidth+1;

        const unsigned int ib=0;
        const unsigned int ie=paddWidth+1;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=lookUpVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            T* lookUpVec = &eleDGVec[lookUp*m_uiNpE];
            if(!eleDGValid[lookUp])
            {
                getElementNodalValues(zippedVec,lookUpVec,lookUp);
                eleDGValid[lookUp]=true;    
            }
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(lookUpVec,&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::blockDiagonalUnZip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        OCT_DIR_LEFT_DOWN_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_LEFT_UP_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_LEFT_BACK_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_LEFT_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_DOWN_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_UP_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_BACK_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_DOWN_BACK_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_DOWN_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_UP_BACK_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_UP_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);


    }

    template<typename T>
    void Mesh::blockVertexUnZip(const ot::Block & blk,const T* zippedVec, T* unzippedVec, T* eleDGVec, bool* eleDGValid)
    {

        OCT_DIR_LEFT_DOWN_BACK_Unzip(blk,zippedVec,unzippedVec,eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_DOWN_BACK_Unzip(blk,zippedVec,unzippedVec,eleDGVec, eleDGValid);
        OCT_DIR_LEFT_UP_BACK_Unzip(blk,zippedVec,unzippedVec,eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_UP_BACK_Unzip(blk,zippedVec,unzippedVec,eleDGVec, eleDGValid);

        OCT_DIR_LEFT_DOWN_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_DOWN_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_LEFT_UP_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
        OCT_DIR_RIGHT_UP_FRONT_Unzip(blk,zippedVec,unzippedVec, eleDGVec, eleDGValid);
    }

    template<typename T>
    void Mesh::child2ParentInjection(const T* in , T* out, unsigned int* child, unsigned int lev) const 
    {

        for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
        {
            if(child[cnum] == LOOK_UP_TABLE_DEFAULT || m_uiAllElements[child[cnum]].getLevel()!= lev ||  !m_uiIsNodalMapValid[child[cnum]]) continue;
            
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
            for(unsigned int j=0;j<m_uiElementOrder+1;j++)
            for(unsigned int i=0;i<m_uiElementOrder+1;i++)
            {
                const bool isHanging=this->isNodeHanging(child[cnum],i,j,k);
                if(isHanging)
                {
                    out[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=in[m_uiE2NMapping_CG[child[cnum]*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                }
                else
                {
                    const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                    const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                    const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                    //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                    if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                    {
                        out[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = in [m_uiE2NMapping_CG [child[cnum]*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j *(m_uiElementOrder+1)+ i ] ];
                        
                    }

                }

            }
        }

    }

    template<typename T>
    void Mesh::unzip(const T* in, T* out, const unsigned int *blkIDs, unsigned int numblks, unsigned int dof)
    {
        if( (!m_uiIsActive) || (m_uiLocalBlockList.empty())  ) return;

        ot::TreeNode blkNode;
        unsigned int ei,ej,ek; // element wise xyz coordinates.
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        unsigned int regLev;
        //unsigned int blkNpe_1D;

        unsigned int lookUp;
        unsigned int lookUp1;
        unsigned int cnum;
        unsigned int faceCnum;

        unsigned int faceNeighCnum1[4]={0,0,0,0}; // immidiate neighbors
        unsigned int faceNeighCnum2[4]={0,0,0,0}; // neighbor's neighbors


        register unsigned int nodeLookUp_CG;
        register unsigned int nodeLookUp_DG;

        std::vector<T> interpOrInjectionOut; // interpolation or injection output.
        std::vector<T> injectionInput;// input for the injection (values from all the 8 children) (This should be put in the order of the morton ordering. )
        std::vector<T> interpolationInput;

        std::vector<T> edgeInterpIn;
        std::vector<T> edgeInterpOut;

        std::vector<T> faceInterpIn;
        std::vector<T> faceInterpOut;

        std::vector<unsigned int> edgeIndex;
        std::vector<unsigned int > faceIndex;
        std::vector<unsigned int > child;
        child.resize(NUM_CHILDREN);

        interpOrInjectionOut.resize(m_uiNpE);
        interpolationInput.resize(m_uiNpE);
        //injectionInput.resize(m_uiNpE*NUM_CHILDREN);

        std::vector<T> injectionTest;
        injectionTest.resize(m_uiNpE*NUM_CHILDREN);

        edgeIndex.resize((m_uiElementOrder+1));
        faceIndex.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));

        edgeInterpIn.resize((m_uiElementOrder+1));
        edgeInterpOut.resize((m_uiElementOrder+1));

        faceInterpIn.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));
        faceInterpOut.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));

        unsigned int mid_bit=0;
        unsigned int sz;
        bool isHanging;
        unsigned int ownerID,ii_x,jj_y,kk_z;
        unsigned int eleIndexMin=0;
        unsigned int eleIndexMax=0;
        bool edgeHanging;
        bool faceHanging;

        unsigned int lx,ly,lz,offset,paddWidth;
        bool isParentValue=false;

        unsigned int fid[(NUM_CHILDREN>>1u)];
        unsigned int cid[(NUM_CHILDREN>>1u)];

        /*if(!rank) std::cout<<"begin unzip "<<std::endl;*/
        #ifdef DEBUG_UNZIP_OP
            double d_min,d_max;
            d_min=-0.5;
            d_max=0.5;
            double x,y,z;
            unsigned int x1,y1,z1;
            std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){ return (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
        #endif

        // NOTE: Be careful when you access ghost elements for padding. (You should only access the level 1 ghost elements. You should not access the level 2 ghost elements at any time. )
        paddWidth = m_uiLocalBlockList[0].get1DPadWidth();

        if((m_uiElementOrder+1)/2 < paddWidth)
        {
            std::cout<<"rank: "<<m_uiActiveRank<<" paddiging with size : "<<paddWidth<<" is too large for element order : "<<m_uiElementOrder<<std::endl;
            MPI_Abort(m_uiCommGlobal,0);
        }

       
        assert(numblks<=m_uiLocalBlockList.size());

        std::vector<T> ele_dg_vec;
        ele_dg_vec.resize(m_uiNumTotalElements*m_uiNpE,(T)0);
        bool* eleVec_valid = new bool[m_uiAllElements.size()];

        for(unsigned int v=0; v < dof; v++)
        {
            const T* zippedVec = in  + v * m_uiNumActualNodes;
            T* unzippedVec = out + v* m_uiUnZippedVecSz; 

            for(unsigned int ii=0; ii <m_uiAllElements.size();ii++)
                eleVec_valid[ii]=false;
            
            for(unsigned int b = 0; b < numblks; b++ )
            {
                const unsigned int blk = blkIDs[b];
                blkNode=m_uiLocalBlockList[blk].getBlockNode();
                assert(blkNode.maxX()<=m_uiMeshDomain_max && blkNode.minX()>=m_uiMeshDomain_min);
                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                //blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*GHOST_WIDTH;
                //std::cout<<"rank: "<<m_uiActiveRank<<" -- blkNpw_1D: "<<blkNpe_1D<<" blkNode: "<<blkNode<<" regLev: "<<regLev<<std::endl;

                sz=1u<<(m_uiMaxDepth-regLev);
                eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;
                assert(eleIndexMax>=eleIndexMin);

                lx=m_uiLocalBlockList[blk].getAllocationSzX();
                ly=m_uiLocalBlockList[blk].getAllocationSzY();
                lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                offset=m_uiLocalBlockList[blk].getOffset();
                paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();


                for(unsigned int elem=m_uiLocalBlockList[blk].getLocalElementBegin();elem<m_uiLocalBlockList[blk].getLocalElementEnd();elem++)
                {
                    ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                    //std::cout<<"blk: "<<blk<<" : "<<blkNode<<" ek: "<<(ek)<<" ej: "<<(ej)<<" ei: "<<(ei)<<" elem: "<<m_uiAllElements[elem]<<std::endl;
                    assert(pNodes[elem].getLevel()==regLev); // this is enforced by block construction
                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_internal.start();
                    #endif

                    T* lookUpElementVec = &ele_dg_vec[elem*m_uiNpE];
                    if(!eleVec_valid[elem])
                    {
                        this->getElementNodalValues(zippedVec,lookUpElementVec,elem);
                        eleVec_valid[elem]=true;
                    }
                    
                    // (1). local nodes copy. Not need to interpolate or inject values. By block construction local octants in the block has is the same level as regular grid.
                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_cpy.start();
                    #endif
                    for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                        for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                            for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_cpy.stop();
                        dendro::timer::t_unzip_sync_internal.stop();
                    #endif
                    // (2). copy the ghost layer (we only copy GHOST_WIDTH amounts of data from the zipped array )z`
                    //---------------------------------------------------------X direction padding --------------------------------------------------------------------------------------------------------------------
                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_face[0].start();
                    #endif
                    if((pNodes[elem].minX()==blkNode.minX()))
                    {
                        assert(ei==eleIndexMin);

                        lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_LEFT];
                        if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                        {

                            if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c1.start();
                                #endif
                                assert(paddWidth<(m_uiElementOrder+1));
                                
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }
                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c1.stop();
                                #endif

                            }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c2.start();
                                #endif
                                assert(pNodes[lookUp].getLevel()+1==regLev);
                                mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                                cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | ((((pNodes[elem].getX()-sz)) >>mid_bit) & 1u));
                                //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()-sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;


                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif


                                #ifdef USE_FD_INTERP_FOR_UNZIP
                                    const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_LEFT,child.data(),fid,cid);
                                    if(st > 0)
                                    {
                                        const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                        this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                        for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                        {
                                            assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_LEFT]]);
                                            assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);
                                            //assert(m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin) ] == blk);

                                            if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                                continue;

                                            this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);


                                            const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                            const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                            const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                            
                                            const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                            const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                            const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                            const unsigned int offset_fd = blk_fd.getOffset();


                                            
                                            const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                            const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                            const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);

                                            assert(paddWidth<(m_uiElementOrder+1));
                                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                            for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                                unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+(ei_fd*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                        }

                                        


                                    }
                                #else
                                    
                                    T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                    if(!eleVec_valid[lookUp])
                                    {
                                        this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                        eleVec_valid[lookUp]=true;
                                    }
                                    this->parent2ChildInterpolation(lookUpElementVec,&(*(interpOrInjectionOut.begin())),cnum);

                                    
                                    assert(paddWidth<(m_uiElementOrder+1));
                                    for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                        for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                            for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                                unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                #endif

                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop() ;
                                    dendro::timer::t_unzip_sync_f_c2.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c3.start();
                                #endif
                                assert(pNodes[lookUp].getLevel()==(regLev+1));
                                //child.resize(NUM_CHILDREN,LOOK_UP_TABLE_DEFAULT);
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                child[1]=lookUp;
                                child[3]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                                child[5]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                                child[7]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                if(m_uiElementOrder ==4 && paddWidth==3)
                                {
                                    // we need to search for the additional points. 
                                    child[0]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_LEFT];
                                    child[2]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_LEFT];
                                    child[4]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_LEFT];
                                    child[6]=m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_LEFT];

                                }else
                                {
                                    child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_LEFT];
                                    child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_LEFT];
                                    child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_LEFT];
                                    child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_LEFT];

                                }


                                this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());

                                
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c3.stop();
                                #endif


                            }

                        }

                    }

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[0].stop();
                    #endif

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[1].start();
                    #endif

                    if((pNodes[elem].maxX()==blkNode.maxX()))
                    {
                        assert(ei==eleIndexMax);
                        lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_RIGHT];
                        if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                            {
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c1.start();
                                #endif
                                assert(paddWidth<(m_uiElementOrder+1));
                                
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }
                                
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=0;i<(paddWidth+1);i++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+((ei+1)*m_uiElementOrder+paddWidth+i)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c1.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c2.start();
                                #endif
                                assert(pNodes[lookUp].getLevel()+1==regLev);
                                mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                                cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | ((((pNodes[elem].getX()+sz)) >>mid_bit) & 1u));
                                //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()+sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;

                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif

                                #ifdef USE_FD_INTERP_FOR_UNZIP
                                    const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_RIGHT,child.data(),fid,cid);
                                    if(st > 0)
                                    {
                                        const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                        this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                        for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                        {
                                            assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_RIGHT]]);
                                            assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                            if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                                continue;

                                            this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);
                                            
                                            const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                            const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                            const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                            
                                            const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                            const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                            const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                            const unsigned int offset_fd = blk_fd.getOffset();


                                            
                                            const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                            const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                            const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);


                                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                                for(unsigned int i=0;i<(paddWidth+1);i++)
                                                unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+((ei_fd+1)*m_uiElementOrder+paddWidth+i)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                        }

                                        


                                    }
                                #else
                                    
                                    T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                    if(!eleVec_valid[lookUp])
                                    {
                                        this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                        eleVec_valid[lookUp]=true;
                                    }

                                    this->parent2ChildInterpolation(lookUpElementVec,&(*(interpOrInjectionOut.begin())),cnum);
                                    assert(paddWidth<(m_uiElementOrder+1));
                                    for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                        for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                            for(unsigned int i=0;i<(paddWidth+1);i++)
                                                unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+((ei+1)*m_uiElementOrder+paddWidth+i)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                #endif



                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c2.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c3.start();
                                #endif
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                child[0]=lookUp;
                                child[2]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                                child[4]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                                child[6]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[6]!=LOOK_UP_TABLE_DEFAULT);

                                if(m_uiElementOrder ==4 && paddWidth==3)
                                {
                                    child[1]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                    child[3]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                                    child[5]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                                    child[7]=m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_RIGHT];

                                }else
                                {
                                    child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                    child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                                    child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                                    child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_RIGHT];

                                }

                                


                                this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());

                                
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=0;i<(paddWidth+1);i++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+((ei+1)*m_uiElementOrder+paddWidth+i)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c3.stop();
                                #endif


                            }

                        }

                    }

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[1].stop();
                    #endif
                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[2].start();
                    #endif

                    //--------------------------------------------------------------------------------------------------- Y Direction----------------------------------------------------------------------------------
                    if((pNodes[elem].minY()==blkNode.minY()))
                    {
                        assert(ej==0);

                        lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_DOWN];
                        if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                        {

                            if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c1.start();
                                #endif
                                assert(paddWidth<(m_uiElementOrder+1));
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c1.stop();
                                #endif

                            }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                            {
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c2.start();
                                #endif
                                assert(pNodes[lookUp].getLevel()+1==regLev);
                                mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                                cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | (((((pNodes[elem].getY()-sz)) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));

                                //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()-sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;
                                

                                //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"parent to child interpolation executed"<<std::endl;
                                assert(paddWidth<(m_uiElementOrder+1));
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif

                                #ifdef USE_FD_INTERP_FOR_UNZIP
                                const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_DOWN,child.data(),fid,cid);
                                if(st > 0)
                                {
                                    const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                    this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                    for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                    {
                                        assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_DOWN]]);
                                        assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                        if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                            continue;
                                            
                                        this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);

                                        const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                        const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                        const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                        
                                        const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                        const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                        const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                        const unsigned int offset_fd = blk_fd.getOffset();

                                        const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);
                                        
                                        for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                            for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                                for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                                    unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                    }

                                    


                                }
                                #else
                                    T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                    if(!eleVec_valid[lookUp])
                                    {
                                        this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                        eleVec_valid[lookUp]=true;
                                    }
                                    this->parent2ChildInterpolation(lookUpElementVec,&(*(interpOrInjectionOut.begin())),cnum);

                                    for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];


                                #endif


                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c2.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                            {
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c3.start();
                                #endif
                                child[2]=lookUp;
                                child[3]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                                assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                                child[6]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                                child[7]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                if(m_uiElementOrder ==4 && paddWidth==3)
                                {
                                    child[0]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_DOWN];
                                    child[1]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_DOWN];
                                    child[4]=m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_DOWN];
                                    child[5]=m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_DOWN];


                                }else
                                {
                                    child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_DOWN];
                                    child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_DOWN];
                                    child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_DOWN];
                                    child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_DOWN];

                                }

                                
                                this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());


                                //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"child to parent interpolation executed"<<std::endl;
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c3.stop();
                                #endif

                            }

                        }

                    }

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[2].stop();
                    #endif
                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[3].start();
                    #endif
                    if((pNodes[elem].maxY()==blkNode.maxY()))
                    {
                        assert(ej==(1u<<(regLev-blkNode.getLevel()))-1);
                        lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_UP];
                        if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c1.start();
                                #endif
                                assert(paddWidth<(m_uiElementOrder+1));
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int j=0;j<(paddWidth+1);j++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+((ej+1)*m_uiElementOrder+paddWidth+j)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c1.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c2.start();
                                #endif
                                assert(pNodes[lookUp].getLevel()+1==regLev);
                                mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                                cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | (((((pNodes[elem].getY()+sz)) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));
                                //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()+sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;
                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif

                                

                                #ifdef USE_FD_INTERP_FOR_UNZIP
                                const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_UP,child.data(),fid,cid);
                                if(st > 0)
                                {
                                    const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                    this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                    for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                    {
                                        assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_UP]]);
                                        assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                        if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                            continue;

                                        this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);


                                        const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                        const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                        const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                        
                                        const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                        const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                        const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                        const unsigned int offset_fd = blk_fd.getOffset();
                                        
                                        const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);
                                        
                                        for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                            for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                                for(unsigned int j=0;j<(paddWidth+1);j++)
                                                    unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+((ej_fd+1)*m_uiElementOrder+paddWidth+j)*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                    }

                                    


                                }
                                #else
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }
                                this->parent2ChildInterpolation(lookUpElementVec,&(*(interpOrInjectionOut.begin())),cnum);

                                assert(paddWidth<(m_uiElementOrder+1));
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int j=0;j<(paddWidth+1);j++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+((ej+1)*m_uiElementOrder+paddWidth+j)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #endif


                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c2.stop();
                                #endif



                            }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c3.start();
                                #endif
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                child[0]=lookUp;
                                child[1]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                                child[4]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                                child[5]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                assert(child[5]!=LOOK_UP_TABLE_DEFAULT);

                                if(m_uiElementOrder ==4 && paddWidth==3)
                                {
                                    child[2]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                    child[3]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                    child[6]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                                    child[7]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];


                                }else
                                {
                                    child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                    child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                    child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                                    child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];
                                }

                                


                                this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());

                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int j=0;j<(paddWidth+1);j++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+((ej+1)*m_uiElementOrder+paddWidth+j)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];


                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c3.stop();
                                #endif

                            }

                        }

                    }

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_face[3].stop();
                    #endif
                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_face[4].start();
                    #endif
                    //--------------------------------------------------------------------- Z direction padding. -------------------------------------------------------------------------------------------------------

                    if((pNodes[elem].minZ()==blkNode.minZ()))
                    {
                        assert(ek==0);

                        lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_BACK];
                        if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                        {

                            if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c1.start();
                                #endif
                                assert(paddWidth<(m_uiElementOrder+1));
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c1.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c2.start();
                                #endif
                                assert(pNodes[lookUp].getLevel()+1==regLev);
                                mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                                cnum=( (((((pNodes[elem].getZ()-sz)) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));
                                //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()-sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;
                                
                                //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"parent to child interpolation executed"<<std::endl;
                                assert(paddWidth<(m_uiElementOrder+1));
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif

                                

                                #ifdef USE_FD_INTERP_FOR_UNZIP
                                const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_BACK,child.data(),fid,cid);
                                if(st > 0)
                                {
                                    const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                    this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                    for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                    {
                                        assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_BACK]]);
                                        assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                        if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                            continue;
                                        
                                        this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);
                                        const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                        const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                        const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                        
                                        const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                        const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                        const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                        const unsigned int offset_fd = blk_fd.getOffset();


                                        
                                        const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);
                                        

                                        for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                            unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                    }

                                    


                                }
                                #else
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }
                                this->parent2ChildInterpolation(lookUpElementVec,&(*(interpOrInjectionOut.begin())),cnum);

                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #endif
                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c2.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c3.start();
                                #endif
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                child[4]=lookUp;
                                child[5]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                                assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                                child[6]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                                assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                                child[7]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];
                                assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                if(m_uiElementOrder ==4 && paddWidth==3)
                                {
                                    child[0]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_BACK];
                                    child[1]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_BACK];
                                    child[2]=m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_BACK];
                                    child[3]=m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_BACK];

                                }else
                                {
                                    child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_BACK];
                                    child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_BACK];
                                    child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_BACK];
                                    child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_BACK];
                                }

                                

                                this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());


                                



                                //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"child to parent interpolation executed"<<std::endl;
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c3.stop();
                                #endif

                            }

                        }

                    }

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[4].stop();
                    #endif
                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[5].start();
                    #endif
                    
                    if((pNodes[elem].maxZ()==blkNode.maxZ()))
                    {
                        assert(ek==(1u<<(regLev-blkNode.getLevel()))-1);
                        lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_FRONT];
                        if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                        {
                            if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                            {

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c1.start();
                                #endif
                                assert(paddWidth<(m_uiElementOrder+1));
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=0;k<(paddWidth+1);k++)
                                            unzippedVec[offset+((ek+1)*m_uiElementOrder+paddWidth+k)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c1.stop();
                                #endif


                            }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                            {

                                assert(pNodes[lookUp].getLevel()+1==regLev);
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c2.start();
                                #endif
                                mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                                cnum=( (((((pNodes[elem].getZ()+sz)) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));
                                //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()+sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;

                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.start();
                                #endif

                                

                                #ifdef USE_FD_INTERP_FOR_UNZIP
                                const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_FRONT,child.data(),fid,cid);
                                if(st > 0)
                                {
                                    const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                    this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                    for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                    {
                                        assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_FRONT]]);
                                        assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                        if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                            continue;
                                        
                                        this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);
                                        
                                        const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                        const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                        const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                        
                                        const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                        const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                        const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                        const unsigned int offset_fd = blk_fd.getOffset();

                                        const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);

                                        for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=0;k<(paddWidth+1);k++)
                                            unzippedVec[offset_fd+((ek_fd+1)*m_uiElementOrder+paddWidth+k)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                    }

                                    
                                }
                                #else
                                T* lookUpElementVec = &ele_dg_vec[lookUp*m_uiNpE];
                                if(!eleVec_valid[lookUp])
                                {
                                    this->getElementNodalValues(zippedVec,lookUpElementVec,lookUp);
                                    eleVec_valid[lookUp]=true;
                                }
                                this->parent2ChildInterpolation(lookUpElementVec,&(*(interpOrInjectionOut.begin())),cnum);
                                assert(paddWidth<(m_uiElementOrder+1));

                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=0;k<(paddWidth+1);k++)
                                            unzippedVec[offset+((ek+1)*m_uiElementOrder+paddWidth+k)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #endif

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c2.stop();
                                #endif



                            }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                            {
                                // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_f_c3.start();
                                #endif
                                child[0]=lookUp;
                                child[1]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                                child[2]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                                child[3]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                assert(child[3]!=LOOK_UP_TABLE_DEFAULT);

                                if(m_uiElementOrder ==4 && paddWidth==3)
                                {
                                    child[4]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                    child[5]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                    child[6]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                    child[7]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];

                                }else
                                {
                                    child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                    child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                    child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                    child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];

                                }

                                this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());


                                

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                                #endif
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                        for(unsigned int k=0;k<(paddWidth+1);k++)
                                            unzippedVec[offset+((ek+1)*m_uiElementOrder+paddWidth+k)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                    dendro::timer::t_unzip_sync_cpy.stop();
                                    dendro::timer::t_unzip_sync_f_c3.stop();
                                #endif


                            }

                        }

                    }

                    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                        dendro::timer::t_unzip_sync_face[5].stop();
                    #endif


                }


                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_edge.start();
                #endif
                    blockDiagonalUnZip(m_uiLocalBlockList[blk],zippedVec,unzippedVec,ele_dg_vec.data(),eleVec_valid);

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_edge.stop();
                #endif

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_vtex.start();
                #endif
                    blockVertexUnZip(m_uiLocalBlockList[blk],zippedVec,unzippedVec,ele_dg_vec.data(),eleVec_valid);

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_vtex.stop();
                #endif

            }

        
        }

        delete [] eleVec_valid;

        

       



        
    }

    template<typename T>
    void Mesh::unzip_scatter(const T* in, T* out,unsigned int dof)
    {
        if(!m_uiIsActive)
            return;
        
        const ot::TreeNode* pNodes   =   m_uiAllElements.data();
        const ot::Block* blkList     =   m_uiLocalBlockList.data();
        const unsigned int eOrder    =   m_uiElementOrder;
        const unsigned int nPe       =   m_uiNpE;

        const unsigned int cgSz  =  this->getDegOfFreedom();
        const unsigned int unSz  =  this->getDegOfFreedomUnZip();   

        const unsigned int* e2n  =  this->getE2NMapping().data();
        const unsigned int* e2e  =  this->getE2EMapping().data();

        
        const unsigned int dgSz  =  nPe;
        std::vector<T> dg_ele_vec;
        dg_ele_vec.resize(dof * dgSz);
        
        T * dgWVec = dg_ele_vec.data();
        T * uzWVec = out;

        std::vector<T> p2cI_all;
        p2cI_all.resize(NUM_CHILDREN * dof * nPe);
        bool p2c_interp_valid[NUM_CHILDREN];

        const double d_compar_tol=1e-10;

        std::vector<ot::TreeNode> childOct;
        childOct.reserve(NUM_CHILDREN);

        for (unsigned int ele=0; ele < m_uiNumTotalElements; ele++)
        {   
            if(m_e2b_unzip_counts[ele]==0)
                continue;

            for(unsigned int ii=0; ii < NUM_CHILDREN; ii++)
                p2c_interp_valid[ii]=false;
            
            //get the elemental_local(dg) values
            for(unsigned int v=0; v< dof; v++)
                this->getElementNodalValues(in + v * cgSz, dgWVec + v*dgSz, ele, false);

            for(unsigned int i = 0; i < m_e2b_unzip_counts[ele]; i++)
            {   
                const unsigned int e2b_offset = m_e2b_unzip_offset[ele];
                const unsigned int blk = m_e2b_unzip_map[e2b_offset + i];
                assert(blk!=LOOK_UP_TABLE_DEFAULT && blk < m_uiLocalBlockList.size());
                
                const unsigned int regLevel  =   blkList[blk].getRegularGridLev();
                const ot::TreeNode blkNode   =   blkList[blk].getBlockNode();
                const unsigned int PW        =   blkList[blk].get1DPadWidth();

                const unsigned int lx     =  blkList[blk].getAllocationSzX();
                const unsigned int ly     =  blkList[blk].getAllocationSzY();
                const unsigned int lz     =  blkList[blk].getAllocationSzZ();
                const unsigned int offset =  blkList[blk].getOffset(); 

                const unsigned int bLev  =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();

                const double  hx = (1u<<(m_uiMaxDepth-bLev))/(double)eOrder;
                const double xmin = blkNode.minX() - PW*hx; const double xmax = blkNode.maxX() + PW*hx;
                const double ymin = blkNode.minY() - PW*hx; const double ymax = blkNode.maxY() + PW*hx;
                const double zmin = blkNode.minZ() - PW*hx; const double zmax = blkNode.maxZ() + PW*hx;


                // no interpolation needed just copy. 
                if(pNodes[ele].getLevel()==bLev)
                {
                    const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/hh;
                    
                    for(unsigned int k=0; k < eOrder+1; k++)
                    {
                        double zz  = pNodes[ele].minZ() + k*hh;
                        
                        if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                        if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                        if(zz < zmin || zz > zmax) 
                            continue;
                        const int kkz = std::round((zz-zmin)*invhh);
                        assert( std::fabs(zz-zmin-kkz*hh) < d_compar_tol);
                        assert(kkz >= 0 && kkz < lz);

                        for(unsigned int j=0; j < eOrder+1; j++)
                        {   
                            double yy  = pNodes[ele].minY() + j*hh;

                            if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                            if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                            if(yy < ymin || yy > ymax) 
                                continue;
                            const int jjy = std::round((yy-ymin)*invhh);
                            //std::cout<<"yy: "<<yy<<" (ymin + hh*jjy): "<<(ymin + hh*jjy)<<std::endl;
                            assert( std::fabs(yy-ymin-jjy*hh) < d_compar_tol);
                            assert(jjy>=0 && jjy<ly);

                            for(unsigned int i=0; i < eOrder+1; i++)
                            {
                                double xx = pNodes[ele].minX() + i*hh;

                                if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                
                                if(xx < xmin || xx > xmax) 
                                    continue;
                                const int iix = std::round((xx-xmin)*invhh);
                                assert( std::fabs(xx-xmin-iix*hh) < d_compar_tol);
                                assert(iix>=0 && iix<lx);

                                //std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<kkz<<" , "<<jjy<<", "<<iix<<")"<<" of : "<<lx<<std::endl;

                                for(unsigned int v=0; v < dof; v++)
                                    uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                    

                            }
                        
                        }
                    
                    }


                }
                else if(pNodes[ele].getLevel() > bLev)
                {
                    assert((bLev+1) == pNodes[ele].getLevel());
                    const unsigned int cnum = pNodes[ele].getMortonIndex();
                    ot::TreeNode tmpParent = pNodes[ele].getParent();
                    
                    const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/(2*hh);

                    assert(eOrder>1);
                    const unsigned int cb =(eOrder%2==0) ? 0 : 1;

                    for(unsigned int k=cb; k < eOrder+1; k+=2)
                    {
                        double zz  = (pNodes[ele].minZ() + k*hh);
                        if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                        if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                        if(zz < zmin || zz > zmax) 
                            continue;
                        const int kkz = std::round((zz-zmin)*invhh);
                        assert(kkz >= 0 && kkz < lz);

                        for(unsigned int j=cb; j < eOrder+1; j+=2)
                        {   
                            double yy  = pNodes[ele].minY() + j*hh;
                            
                            if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                            if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                            if(yy < ymin || yy > ymax) 
                                continue;

                            const int jjy = std::round((yy-ymin)*invhh);
                            assert(jjy>=0 && jjy<ly);

                            for(unsigned int i=cb; i < eOrder+1; i+=2)
                            {
                                double xx = pNodes[ele].minX() + i*hh;

                                if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                
                                if(xx < xmin || xx > xmax) 
                                    continue;
                                const int iix = std::round((xx-xmin)*invhh);
                                assert(iix>=0 && iix<lx);

                                //std::cout<<"blk: "<<blk<<" blk copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;
                                for(unsigned int v=0; v < dof; v++)
                                    uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                    
                                

                            }
                        
                        }
                    
                    }


                
                }
                else
                {   
                    assert((bLev) == (pNodes[ele].getLevel()+1));
                    childOct.clear();
                    pNodes[ele].addChildren(childOct); // note this is the ordering of SFC (depends on Hilbert or Morton. )

                    for(unsigned int child = 0; child < NUM_CHILDREN; child++)
                    {

                        if( (childOct[child].maxX() < xmin || childOct[child].minX() >=xmax)  || (childOct[child].maxY() < ymin || childOct[child].minY() >=ymax) || (childOct[child].maxZ() < zmin || childOct[child].minZ() >=zmax) )
                            continue;


                        //std::cout<<"blk: "<<blk<<" blkNode: "<<blkNode<<" child: "<<child<<" child node "<<childOct[child]<<" parent : "<<pNodes[ele]<<std::endl;
                        const double hh = (1u<<(m_uiMaxDepth - childOct[child].getLevel()))/(double) eOrder;
                        const double invhh = 1.0/hh;
                        
                        const unsigned int cnum = childOct[child].getMortonIndex();
                        if(!p2c_interp_valid[cnum])
                        {
                            for(unsigned int v=0; v < dof; v++)
                                this->parent2ChildInterpolation(&dgWVec[v*dgSz], p2cI_all.data() + cnum * dof * nPe + v * nPe , cnum, m_uiDim);
                            
                            p2c_interp_valid[cnum]=true;
                        }

                        for(unsigned int v=0; v < dof; v++)
                        {
                            
                            const T* const p2cI = p2cI_all.data() + cnum * dof * nPe + v * nPe;
                            for(unsigned int k=0; k < eOrder+1; k++)
                            {
                                double zz  = childOct[child].minZ() + k*hh;
                                
                                if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                                if(fabs(zz-zmax)<d_compar_tol) zz=zmax;
                                
                                if(zz < zmin || zz > zmax) 
                                    continue;
                                const int kkz = std::round((zz-zmin)*invhh);
                                assert(kkz >= 0 && kkz < lz);

                                for(unsigned int j=0; j < eOrder+1; j++)
                                {   
                                    double yy  = childOct[child].minY() + j*hh;
                                    
                                    if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                                    if(fabs(yy-ymax)<d_compar_tol) yy=ymax;
                                    
                                    if(yy < ymin || yy > ymax) 
                                        continue;

                                    const int jjy = std::round((yy-ymin)*invhh);
                                    assert(jjy>=0 && jjy<ly);

                                    for(unsigned int i=0; i < eOrder+1; i++)
                                    {
                                        double xx = childOct[child].minX() + i*hh;
                                        
                                        if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                        if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                        
                                        if(xx < xmin || xx > xmax) 
                                            continue;
                                        const int iix = std::round((xx-xmin)*invhh);
                                        assert(iix>=0 && iix<lx);

                                        uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  p2cI[k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];

                                    }
                                
                                }
                            
                            }
                            
                        }

                    }

                }


            }

        }
        
    }

    template <typename T>
    void Mesh::unzip(const T* in, T* out, unsigned int dof)
    {
        if( (!m_uiIsActive) || (m_uiLocalBlockList.empty())  ) return;

        // std::vector<unsigned int > blkIDs;
        // blkIDs.resize(m_uiLocalBlockList.size());

        // for(unsigned int i=0; i< m_uiLocalBlockList.size(); i++)
        //     blkIDs[i] = i ; 
        // unzip all the blocks. 
        //this->unzip(in,out,blkIDs.data(),blkIDs.size(),dof);
        this->unzip_scatter(in,out,dof);
        

    }

    #if 0
    template<typename T>
    void Mesh::readSpecialPtsBegin(const T* in)
    {

        if(m_uiGlobalNpes==1)
            return;


         // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        std::vector<T> eVec;
        eVec.resize(m_uiNpE);

        if(m_uiIsActive)
        {
            const unsigned int sendBSz=m_uiSendOffsetRePt[m_uiActiveNpes-1] + m_uiSendCountRePt[m_uiActiveNpes-1];
            const unsigned int recvBSz=m_uiRecvOffsetRePt[m_uiActiveNpes-1] + m_uiRecvCountRePt[m_uiActiveNpes-1];

            AsyncExchangeContex ctx(in);
            MPI_Comm commActive= m_uiCommActive;
            unsigned int proc_id;

            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p=0;recv_p<m_uiReqRecvProcList.size();recv_p++)
                {
                    proc_id=m_uiReqRecvProcList[recv_p];
                    MPI_Request* req=new MPI_Request();
                    par::Mpi_Irecv((recvB+m_uiRecvOffsetRePt[proc_id]),m_uiRecvCountRePt[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*sendBSz);
                sendB=(T*)ctx.getSendBuffer();


                const unsigned int nx = m_uiElementOrder + 1;
                const unsigned int ny = m_uiElementOrder + 1;
                const unsigned int nz = m_uiElementOrder + 1;

                std::vector<unsigned int>* ownerList;
                unsigned int ownerID, ii_x, jj_y, kk_z;

                
                for(unsigned int i=0; i< m_uiUnzip_3pt_ele.size(); i++)
                {
                    ot::Key tmpEleKey= m_uiUnzip_3pt_ele[i];
                    assert((tmpEleKey.getFlag() & OCT_FOUND));
                    const unsigned int eleID = tmpEleKey.getSearchResult();
                    this->getElementNodalValues(in,&(*(eVec.begin())),eleID);
                    
                    const unsigned int step_sz = ((1u<< (m_uiMaxDepth - m_uiAllElements[eleID].getLevel()))/m_uiElementOrder);
                    ownerList = tmpEleKey.getOwnerList();
                    for(unsigned int w=0; w< ownerList->size();w++)
                    {
                        
                        const unsigned int ii = (m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].minX() - m_uiAllElements[eleID].minX())/(step_sz); 
                        const unsigned int jj = (m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].minY() - m_uiAllElements[eleID].minY())/(step_sz); 
                        const unsigned int kk = (m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].minZ() - m_uiAllElements[eleID].minZ())/(step_sz);
                        const std::vector<unsigned int > * ownerList1 = m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].getOwnerList();

                        for(unsigned int w1 = 0; w1 < ownerList1->size() ; w1++)
                        {
                            // if(m_uiActiveRank==1 && (*ownerList1)[w1]<18 )
                            //     std::cout<<" rank: "<<m_uiActiveRank<<" putting : "<<m_uiUnzip_3pt_recv_keys[(*ownerList)[w]]<<" to send buf loc: "<<(*ownerList1)[w1]<<std::endl;

                            sendB[(*ownerList1)[w1]] = eVec[kk * ny * nx + jj * nx + ii];
                        }
                            
                        
                    }
                }


                // active send procs
                for(unsigned int send_p=0;send_p<m_uiReqSendProcList.size();send_p++)
                {
                    proc_id=m_uiReqSendProcList[send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+m_uiSendOffsetRePt[proc_id],m_uiSendCountRePt[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }


            }

            m_uiCommTag++;
            m_uiMPIContexts.push_back(ctx);



        }

    }

    template <typename T>
    void Mesh::readSpecialPtsEnd(const T *in, T* out)
    {
        if(m_uiGlobalNpes == 1)
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiIsActive)
        {
            const unsigned int sendBSz=m_uiSendOffsetRePt[m_uiActiveNpes-1] + m_uiSendCountRePt[m_uiActiveNpes-1];
            const unsigned int recvBSz=m_uiRecvOffsetRePt[m_uiActiveNpes-1] + m_uiRecvCountRePt[m_uiActiveNpes-1];

            //std::cout<<"rank: "<<m_uiActiveRank<<" recv sz: "<<recvBSz<<std::endl;

            unsigned int proc_id;
            unsigned int ctxIndex=0;

            for(unsigned int i=0;i<m_uiMPIContexts.size();i++)
            {
                if(m_uiMPIContexts[i].getBuffer()==in)
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
                std::memcpy(out,recvB,sizeof(T)*recvBSz);
                
                // for(unsigned int i=0; i<recvBSz ; i++ )
                //     out[i] = recvB[i];
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
    #endif

    template<typename T>
    int Mesh::getFaceNeighborValues(unsigned int eleID, const T* in, T* out, T* coords, unsigned int * neighID, unsigned int face, NeighbourLevel & level) const
    {

        if(!m_uiIsActive)
            return (0);

        const unsigned int lookUp=m_uiE2EMapping[eleID*m_uiNumDirections + face];
        if(lookUp==LOOK_UP_TABLE_DEFAULT)
            return (0);

        const unsigned int l1=m_uiAllElements[eleID].getLevel();
        const unsigned int l2=m_uiAllElements[lookUp].getLevel();

        for(unsigned int i=0;i<(NUM_CHILDREN>>1);i++)
            neighID[i]=LOOK_UP_TABLE_DEFAULT;

        int num_face_neighbours = 1;
        if(l1==l2)
        {
            // both elements are in the same level.
            level = NeighbourLevel::SAME;
            neighID[0]=lookUp;
            this->getElementNodalValues(in,out,lookUp);

            // coordinate computation
            const ot::TreeNode lookUpOct=m_uiAllElements[lookUp];
            const unsigned int sz=1u<<(m_uiMaxDepth-lookUpOct.getLevel());

            for(unsigned int k=0;k< (m_uiElementOrder+1);k++ )
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 0] = lookUpOct.minX() + i * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 1] = lookUpOct.minY() + j * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 2] = lookUpOct.minZ() + k * (sz/(T)m_uiElementOrder);
                    }
        }else if(l2<l1)
        {
            level = NeighbourLevel::COARSE;
            // lookUp octant is coaser than eleID.
            neighID[0]=lookUp;
            this->getElementNodalValues(in,out + m_uiNpE ,lookUp);

            unsigned int sz1=1u<<(m_uiMaxDepth-m_uiAllElements[eleID].getLevel());

            unsigned int x = m_uiAllElements[eleID].minX();
            unsigned int y = m_uiAllElements[eleID].minY();
            unsigned int z = m_uiAllElements[eleID].minZ();



            ot::TreeNode tmpOct;
            unsigned int cnum;
            switch (face)
            {
                case OCT_DIR_LEFT:
                    tmpOct=ot::TreeNode(x-sz1,y,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_RIGHT:
                    tmpOct=ot::TreeNode(x+sz1,y,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_DOWN:
                    tmpOct=ot::TreeNode(x,y-sz1,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_UP:
                    tmpOct=ot::TreeNode(x,y+sz1,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;


                case OCT_DIR_BACK:
                    tmpOct=ot::TreeNode(x,y,z-sz1,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_FRONT:
                    tmpOct=ot::TreeNode(x,y,z+sz1,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                default:
                    std::cout<<"global rank : "<<m_uiGlobalRank<<" dir: "<<face<<" is invalid. Function : "<<__func__<<std::endl;
                    MPI_Abort(m_uiCommGlobal,0);
                    break;
            }


            this->parent2ChildInterpolation(out+m_uiNpE,out,cnum,m_uiDim);

            // coordinate computation
            const ot::TreeNode lookUpOct=tmpOct;
            const unsigned int sz=1u<<(m_uiMaxDepth-lookUpOct.getLevel());

            for(unsigned int k=0;k< (m_uiElementOrder+1);k++ )
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 0] = lookUpOct.minX() + i * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 1] = lookUpOct.minY() + j * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 2] = lookUpOct.minZ() + k * (sz/(T)m_uiElementOrder);
                    }


        }else{
            // lookUp octant is finer than eleID.
            assert(l2>l1);

            unsigned int dir, dirOp, dir1,dir2;
            num_face_neighbours = 4;
            level = NeighbourLevel::REFINE;
            switch (face)
            {

                case OCT_DIR_LEFT:

                    dir=OCT_DIR_LEFT;
                    dirOp=OCT_DIR_RIGHT;
                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_UP;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                case OCT_DIR_RIGHT:

                    dir=OCT_DIR_RIGHT;
                    dirOp=OCT_DIR_LEFT;

                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_UP;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];

                    break;

                case OCT_DIR_DOWN:

                    dir=OCT_DIR_DOWN;
                    dirOp=OCT_DIR_UP;

                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_RIGHT;


                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];

                    break;

                case OCT_DIR_UP:

                    dir=OCT_DIR_UP;
                    dirOp=OCT_DIR_DOWN;

                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_RIGHT;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                case OCT_DIR_BACK:

                    dir=OCT_DIR_BACK;
                    dirOp=OCT_DIR_FRONT;

                    dir1=OCT_DIR_UP;
                    dir2=OCT_DIR_RIGHT;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                case OCT_DIR_FRONT:

                    dir=OCT_DIR_FRONT;
                    dirOp=OCT_DIR_BACK;

                    dir1=OCT_DIR_UP;
                    dir2=OCT_DIR_RIGHT;


                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                default:
                    std::cout<<"global rank : "<<m_uiGlobalRank<<" dir: "<<face<<" is invalid. Function : "<<__func__<<std::endl;
                    MPI_Abort(m_uiCommGlobal,0);


            }

            for(unsigned int child=0;child<(NUM_CHILDREN>>1);child++)
            {

                this->getElementNodalValues(in,out+child*m_uiNpE,neighID[child]);

                // coordinate computation
                const ot::TreeNode lookUpOct=m_uiAllElements[neighID[child]];
                const unsigned int sz=1u<<(m_uiMaxDepth-lookUpOct.getLevel());

                for(unsigned int k=0;k< (m_uiElementOrder+1);k++ )
                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                        {
                            coords[ child*m_uiNpE*m_uiDim + m_uiDim * (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 0] = lookUpOct.minX() + i * (sz/(T)m_uiElementOrder);
                            coords[ child*m_uiNpE*m_uiDim + m_uiDim * (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 1] = lookUpOct.minY() + j * (sz/(T)m_uiElementOrder);
                            coords[ child*m_uiNpE*m_uiDim + m_uiDim * (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 2] = lookUpOct.minZ() + k * (sz/(T)m_uiElementOrder);
                        }
            }




        }

        return num_face_neighbours;

    }

    template<typename T>
    void Mesh::getUnzipElementalNodalValues(const T* uzipVec, unsigned int blkID, unsigned int ele, T*out, bool isPadded) const
    {
        const ot::Block block = m_uiLocalBlockList[blkID];
        ot::TreeNode blkNode = m_uiLocalBlockList[blkID].getBlockNode();
        const unsigned int eleBegin = block.getLocalElementBegin();
        const unsigned int eleEnd = block.getLocalElementEnd();

        assert(eleBegin <= ele && ele < eleEnd);
        const unsigned int regLev=block.getRegularGridLev();
        const unsigned int lx=block.getAllocationSzX();
        const unsigned int ly=block.getAllocationSzY();
        const unsigned int lz=block.getAllocationSzZ();
        const unsigned int offset=block.getOffset();
        const unsigned int paddWidth=block.get1DPadWidth();

        const unsigned int ei=(m_uiAllElements[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
        const unsigned int ej=(m_uiAllElements[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
        const unsigned int ek=(m_uiAllElements[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);
        const unsigned int eleIDMax = m_uiLocalBlockList[blkID].getElemSz1D();

        

        if(isPadded)
        {
            const unsigned int ib = ei*m_uiElementOrder;
            const unsigned int ie = ei*m_uiElementOrder + (m_uiElementOrder+1) + 2*paddWidth ;

            const unsigned int jb = ej*m_uiElementOrder;
            const unsigned int je = ej*m_uiElementOrder + (m_uiElementOrder+1) + 2*paddWidth ;

            const unsigned int kb = ek*m_uiElementOrder;
            const unsigned int ke = ek*m_uiElementOrder + (m_uiElementOrder+1) + 2*paddWidth ;

            const unsigned int en[3] = {(m_uiElementOrder+1) + 2*paddWidth , (m_uiElementOrder+1) + 2*paddWidth, (m_uiElementOrder+1) + 2*paddWidth }; 

            for(unsigned int k=kb; k< ke; k++)
             for(unsigned int j=jb; j< je; j++)
              for(unsigned int i=ib; i< ie; i++)
               out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + j*lx + i];

            

            // copy the unzip element last point to the padding region. 
            if(m_uiAllElements[ele].minX()==0)
            {
                assert(ei==0);

                for(unsigned int k=kb; k< ke; k++)
                for(unsigned int j=jb; j< je; j++)
                for(unsigned int i=ib; i< paddWidth; i++)
                    out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + j*lx + paddWidth];

            }


            if(m_uiAllElements[ele].minY()==0)
            {
                assert(ej==0);

                for(unsigned int k=kb; k< ke; k++)
                for(unsigned int j=jb; j< paddWidth; j++)
                for(unsigned int i=ib; i< ie; i++)
                    out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + paddWidth*lx + i];

            }

            if(m_uiAllElements[ele].minZ()==0)
            {
                assert(ek==0);

                for(unsigned int k=kb; k< paddWidth; k++)
                for(unsigned int j=jb; j< je; j++)
                for(unsigned int i=ib; i< ie; i++)
                    out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + paddWidth*ly*lx + j*lx + i];

            }


            if(m_uiAllElements[ele].maxX()==(1u<<m_uiMaxDepth))
            {
                assert(ei==(eleIDMax-1));

                for(unsigned int k=kb; k< ke; k++)
                for(unsigned int j=jb; j< je; j++)
                for(unsigned int i=(ie-paddWidth); i< ie; i++)
                    out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + j*lx + (ie-paddWidth-1)];

            }

            if(m_uiAllElements[ele].maxY()==(1u<<m_uiMaxDepth))
            {
                assert(ej==(eleIDMax-1));

                for(unsigned int k=kb; k< ke; k++)
                for(unsigned int j=(je-paddWidth); j< je; j++)
                for(unsigned int i=ib; i< ie; i++)
                    out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + (je-paddWidth-1)*lx + i];

            }

            if(m_uiAllElements[ele].maxZ()==(1u<<m_uiMaxDepth))
            {
                assert(ek==(eleIDMax-1));

                for(unsigned int k=(ke-paddWidth); k< ke; k++)
                for(unsigned int j=jb; j< je; j++)
                for(unsigned int i=ib; i< ie; i++)
                    out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + (ke-paddWidth-1)*ly*lx + j*lx + i];

            }

            

        }else
        {
            const unsigned int ib = ei*m_uiElementOrder + paddWidth;
            const unsigned int ie = ei*m_uiElementOrder + (m_uiElementOrder+1) ;

            const unsigned int jb = ej*m_uiElementOrder + paddWidth ;
            const unsigned int je = ej*m_uiElementOrder + (m_uiElementOrder+1) ;

            const unsigned int kb = ek*m_uiElementOrder + paddWidth ;
            const unsigned int ke = ek*m_uiElementOrder + (m_uiElementOrder+1) ;

            const unsigned int en[3] = {(m_uiElementOrder+1) , (m_uiElementOrder+1) , (m_uiElementOrder+1) };

            for(unsigned int k=kb; k< ke; k++)
             for(unsigned int j=jb; j< je; j++)
              for(unsigned int i=ib; i< ie; i++)
                out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + j*lx + i];



        }
        
        
    }


    template<typename T>
    void Mesh::getBlkBoundaryParentNodes(const T* zipVec, T* out, T*w1, T*w2, unsigned int lookUp, const unsigned int * fid, const unsigned int* cid,const unsigned int * child)
    {

        const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
        const unsigned int eorder_by2 = (m_uiElementOrder+1)>>1u;
        const unsigned int nx = m_uiElementOrder + 1;
        const unsigned int ny = m_uiElementOrder + 1;
        const unsigned int nz = m_uiElementOrder + 1;
        
        unsigned char bit[3];

        // finner elements. 
        for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
        {
            this->getElementNodalValues(zipVec,w1,child[fid[w]]);
            //std::cout<<" cnum : "<<fid[w]<<std::endl;
            bit[0] = binOp::getBit(fid[w],0);
            bit[1] = binOp::getBit(fid[w],1);
            bit[2] = binOp::getBit(fid[w],2);

            const unsigned int kb = bit[2] * eorder_by2 ;
            unsigned int ke = kb + eorder_by2 +1;

            const unsigned int jb = bit[1] * eorder_by2 ;
            unsigned int je = jb + eorder_by2 +1;

            const unsigned int ib = bit[0] * eorder_by2 ;
            unsigned int ie = ib + eorder_by2 +1;

            for(unsigned int k=0; k< nz; k+=2)
                for(unsigned int j=0; j < ny; j+=2)
                for(unsigned int i=0; i < nx; i+=2)
                out[(kb + (k>>1u))*ny*nx + (jb + (j>>1u))*nx + (ib + (i>>1u))] = w1[k*ny*nx + j*nx + i];

        }

        this->getElementNodalValues(zipVec,w1,lookUp);
        // coarser elements. 
        for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
        {
            this->parent2ChildInterpolation(w1,w2,fid[w],m_uiDim);

            //std::cout<<" cnum : "<<fid[w]<<std::endl;
            bit[0] = binOp::getBit(cid[w],0);
            bit[1] = binOp::getBit(cid[w],1);
            bit[2] = binOp::getBit(cid[w],2);

            const unsigned int kb = bit[2] * eorder_by2 ;
            unsigned int ke = kb + eorder_by2 +1;

            const unsigned int jb = bit[1] * eorder_by2 ;
            unsigned int je = jb + eorder_by2 +1;

            const unsigned int ib = bit[0] * eorder_by2 ;
            unsigned int ie = ib + eorder_by2 +1;

            for(unsigned int k=0; k< nz; k+=2)
                for(unsigned int j=0; j < ny; j+=2)
                for(unsigned int i=0; i < nx; i+=2)
                {
                    // std::cout<< " cnum : "<<fid[w]<< " left lookup value: ijk: "<<i<<j<<k<<" "<<lookUpElementVec[(kb + (k>>1u))*ny*nx + (jb + (j>>1u))*nx + (ib + (i>>1u))]<< " inject value: "<<w2[k*ny*nx + j*nx + i]<<"";
                    // printf("lookup idx (%d,%d,%d)\n",(ib + (i>>1u)),(jb + (j>>1u)), (kb + (k>>1u)));
                    out[(kb + (k>>1u))*ny*nx + (jb + (j>>1u))*nx + (ib + (i>>1u))] = w2[k*ny*nx + j*nx + i];
                }
                

        }

    }


    template<typename T>
    void Mesh::unzipDG(const T* in, T* out, const unsigned int *blkIDs, unsigned int numblks,unsigned int dof)
    {
        if(!m_uiIsActive)
            return;
        
        const ot::TreeNode* pNodes   =   m_uiAllElements.data();
        const ot::Block* blkList     =   m_uiLocalBlockList.data();
        const unsigned int eOrder    =   m_uiElementOrder;
        const unsigned int nPe       =   m_uiNpE;

        const unsigned int dgSz  =  m_uiAllElements.size() * nPe;
        const unsigned int cgSz  =  this->getDegOfFreedom();
        const unsigned int unSz  =  this->getDegOfFreedomUnZip();   

        const unsigned int* e2n  =  this->getE2NMapping().data();
        const unsigned int* e2e  =  this->getE2EMapping().data();

        const T * dgWVec = in;
        T * uzWVec = out;

        for(unsigned int bid=0; bid < numblks; bid++)
        {
            const unsigned int blk = blkIDs[bid];
            const unsigned int regLevel  =   blkList[blk].getRegularGridLev();
            const ot::TreeNode blkNode   =   blkList[blk].getBlockNode();
            const unsigned int PW        =   blkList[blk].get1DPadWidth();

            const unsigned int lx     =  blkList[blk].getAllocationSzX();
            const unsigned int ly     =  blkList[blk].getAllocationSzY();
            const unsigned int lz     =  blkList[blk].getAllocationSzZ();
            const unsigned int offset =  blkList[blk].getOffset(); 

            const unsigned int bLev  =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
            
            std::vector<unsigned int> eid;
            eid.reserve((NUM_CHILDREN+ NUM_FACES + NUM_EDGES + 1) * 4);
            this->blkUnzipElementIDs(blk,eid);
            
            // now need to copy to the block unzip/ block asyncVector 
            const double  hx = (1u<<(m_uiMaxDepth-bLev))/(double)eOrder;
            
            const double xmin = blkNode.minX() - PW*hx; const double xmax = blkNode.maxX() + PW*hx;
            const double ymin = blkNode.minY() - PW*hx; const double ymax = blkNode.maxY() + PW*hx;
            const double zmin = blkNode.minZ() - PW*hx; const double zmax = blkNode.maxZ() + PW*hx;

            std::vector<ot::TreeNode> childOct;
            childOct.reserve(NUM_CHILDREN);

            std::vector<T> p2cI;
            p2cI.resize(nPe);

            const double d_compar_tol=1e-10;

            

            for(unsigned int m=0; m < eid.size(); m++)
            {
                const unsigned int ele = eid[m];
                
                // no interpolation needed just copy. 
                if(pNodes[ele].getLevel()==bLev)
                {
                    const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/hh;
                    
                    for(unsigned int k=0; k < eOrder+1; k++)
                    {
                        double zz  = pNodes[ele].minZ() + k*hh;
                        
                        if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                        if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                        if(zz < zmin || zz > zmax) 
                            continue;
                        const int kkz = std::round((zz-zmin)*invhh);
                        assert( std::fabs(zz-zmin-kkz*hh) < d_compar_tol);
                        assert(kkz >= 0 && kkz < lz);

                        for(unsigned int j=0; j < eOrder+1; j++)
                        {   
                            double yy  = pNodes[ele].minY() + j*hh;

                            if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                            if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                            if(yy < ymin || yy > ymax) 
                                continue;
                            const int jjy = std::round((yy-ymin)*invhh);
                            //std::cout<<"yy: "<<yy<<" (ymin + hh*jjy): "<<(ymin + hh*jjy)<<std::endl;
                            assert( std::fabs(yy-ymin-jjy*hh) < d_compar_tol);
                            assert(jjy>=0 && jjy<ly);

                            for(unsigned int i=0; i < eOrder+1; i++)
                            {
                                double xx = pNodes[ele].minX() + i*hh;

                                if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                
                                if(xx < xmin || xx > xmax) 
                                    continue;
                                const int iix = std::round((xx-xmin)*invhh);
                                assert( std::fabs(xx-xmin-iix*hh) < d_compar_tol);
                                assert(iix>=0 && iix<lx);

                                //std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<kkz<<" , "<<jjy<<", "<<iix<<")"<<" of : "<<lx<<std::endl;

                                for(unsigned int v=0; v < dof; v++)
                                    uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                    

                            }
                        
                        }
                    
                    }


                }
                else if(pNodes[ele].getLevel() > bLev)
                {
                    assert((bLev+1) == pNodes[ele].getLevel());
                    const unsigned int cnum = pNodes[ele].getMortonIndex();
                    ot::TreeNode tmpParent = pNodes[ele].getParent();
                    
                    const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/(2*hh);

                    assert(eOrder>1);
                    const unsigned int cb =(eOrder%2==0) ? 0 : 1;

                    for(unsigned int k=cb; k < eOrder+1; k+=2)
                    {
                        double zz  = (pNodes[ele].minZ() + k*hh);
                        if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                        if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                        if(zz < zmin || zz > zmax) 
                            continue;
                        const int kkz = std::round((zz-zmin)*invhh);
                        assert(kkz >= 0 && kkz < lz);

                        for(unsigned int j=cb; j < eOrder+1; j+=2)
                        {   
                            double yy  = pNodes[ele].minY() + j*hh;
                            
                            if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                            if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                            if(yy < ymin || yy > ymax) 
                                continue;

                            const int jjy = std::round((yy-ymin)*invhh);
                            assert(jjy>=0 && jjy<ly);

                            for(unsigned int i=cb; i < eOrder+1; i+=2)
                            {
                                double xx = pNodes[ele].minX() + i*hh;

                                if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                
                                if(xx < xmin || xx > xmax) 
                                    continue;
                                const int iix = std::round((xx-xmin)*invhh);
                                assert(iix>=0 && iix<lx);

                                //std::cout<<"blk: "<<blk<<" blk copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;
                                for(unsigned int v=0; v < dof; v++)
                                    uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                    
                                

                            }
                        
                        }
                    
                    }


                
                }
                else
                {
                    assert((bLev) == (pNodes[ele].getLevel()+1));
                    childOct.clear();
                    pNodes[ele].addChildren(childOct); // note this is the ordering of SFC (depends on Hilbert or Morton. )

                    for(unsigned int child = 0; child < NUM_CHILDREN; child++)
                    {

                        if( (childOct[child].maxX() < xmin || childOct[child].minX() >=xmax)  || (childOct[child].maxY() < ymin || childOct[child].minY() >=ymax) || (childOct[child].maxZ() < zmin || childOct[child].minZ() >=zmax) )
                            continue;


                        //std::cout<<"blk: "<<blk<<" blkNode: "<<blkNode<<" child: "<<child<<" child node "<<childOct[child]<<" parent : "<<pNodes[ele]<<std::endl;
                        const double hh = (1u<<(m_uiMaxDepth - childOct[child].getLevel()))/(double) eOrder;
                        const double invhh = 1.0/hh;

                        for(unsigned int v=0; v < dof; v++)
                        {
                            const unsigned int cnum = childOct[child].getMortonIndex();
                            this->parent2ChildInterpolation(&dgWVec[v*dgSz + ele*nPe],p2cI.data(),cnum,m_uiDim);

                            for(unsigned int k=0; k < eOrder+1; k++)
                            {
                                double zz  = childOct[child].minZ() + k*hh;
                                
                                if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                                if(fabs(zz-zmax)<d_compar_tol) zz=zmax;
                                
                                if(zz < zmin || zz > zmax) 
                                    continue;
                                const int kkz = std::round((zz-zmin)*invhh);
                                assert(kkz >= 0 && kkz < lz);

                                for(unsigned int j=0; j < eOrder+1; j++)
                                {   
                                    double yy  = childOct[child].minY() + j*hh;
                                    
                                    if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                                    if(fabs(yy-ymax)<d_compar_tol) yy=ymax;
                                    
                                    if(yy < ymin || yy > ymax) 
                                        continue;

                                    const int jjy = std::round((yy-ymin)*invhh);
                                    assert(jjy>=0 && jjy<ly);

                                    for(unsigned int i=0; i < eOrder+1; i++)
                                    {
                                        double xx = childOct[child].minX() + i*hh;
                                        
                                        if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                        if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                        
                                        if(xx < xmin || xx > xmax) 
                                            continue;
                                        const int iix = std::round((xx-xmin)*invhh);
                                        assert(iix>=0 && iix<lx);

                                        uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  p2cI[k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];

                                    }
                                
                                }
                            
                            }
                            
                        }

                    }

                }

            }
            


            // internal copy. 
            for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
            {
                const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                const unsigned int emin = 0;
                const unsigned int emax = (1u<<(regLevel-blkNode.getLevel()))-1;

                // #pragma unroll
                // for(unsigned int v=0; v < dof; v++)
                //     std::memcpy(dgWVec + v*dgSz + elem * nPe , dgStages[rk_s-1] + v*dgSz + elem * nPe, sizeof(T)*nPe );
                
                for(unsigned int v=0; v < dof; v++)
                for(unsigned int k=0;k<(eOrder+1);k++)
                for(unsigned int j=0;j<(eOrder+1);j++)
                for(unsigned int i=0;i<(eOrder+1);i++)
                    uzWVec[v*unSz + offset + (ek*eOrder+k+PW)*(ly*lx)+(ej*eOrder+j+PW)*(lx)+(ei*eOrder+i+PW)] = dgWVec[v*dgSz + elem*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];;

            }


        }
        

    }

    template<typename T>
    void Mesh::unzipDG_scatter(const T* in, T* out,unsigned int dof)
    {
        if(!m_uiIsActive)
            return;
        
        const ot::TreeNode* pNodes   =   m_uiAllElements.data();
        const ot::Block* blkList     =   m_uiLocalBlockList.data();
        const unsigned int eOrder    =   m_uiElementOrder;
        const unsigned int nPe       =   m_uiNpE;

        const unsigned int dgSz  =  m_uiAllElements.size() * nPe;
        const unsigned int cgSz  =  this->getDegOfFreedom();
        const unsigned int unSz  =  this->getDegOfFreedomUnZip();   

        const unsigned int* e2n  =  this->getE2NMapping().data();
        const unsigned int* e2e  =  this->getE2EMapping().data();

        const T * dgWVec = in;
        T * uzWVec = out;

        std::vector<T> p2cI_all;
        p2cI_all.resize(NUM_CHILDREN * dof * nPe);
        bool p2c_interp_valid[NUM_CHILDREN];

        const double d_compar_tol=1e-10;

        std::vector<ot::TreeNode> childOct;
        childOct.reserve(NUM_CHILDREN);

        for (unsigned int ele=0; ele < m_uiNumTotalElements; ele++)
        {   
            if(m_e2b_unzip_counts[ele]==0)
                continue;
            
            for(unsigned int ii=0; ii < NUM_CHILDREN; ii++)
                p2c_interp_valid[ii]=false;

            for(unsigned int i = 0; i < m_e2b_unzip_counts[ele]; i++)
            {   
                const unsigned int e2b_offset = m_e2b_unzip_offset[ele];
                const unsigned int blk = m_e2b_unzip_map[e2b_offset + i];
                assert(blk!=LOOK_UP_TABLE_DEFAULT && blk < m_uiLocalBlockList.size());
                
                const unsigned int regLevel  =   blkList[blk].getRegularGridLev();
                const ot::TreeNode blkNode   =   blkList[blk].getBlockNode();
                const unsigned int PW        =   blkList[blk].get1DPadWidth();

                const unsigned int lx     =  blkList[blk].getAllocationSzX();
                const unsigned int ly     =  blkList[blk].getAllocationSzY();
                const unsigned int lz     =  blkList[blk].getAllocationSzZ();
                const unsigned int offset =  blkList[blk].getOffset(); 

                const unsigned int bLev  =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();

                const double  hx = (1u<<(m_uiMaxDepth-bLev))/(double)eOrder;
                const double xmin = blkNode.minX() - PW*hx; const double xmax = blkNode.maxX() + PW*hx;
                const double ymin = blkNode.minY() - PW*hx; const double ymax = blkNode.maxY() + PW*hx;
                const double zmin = blkNode.minZ() - PW*hx; const double zmax = blkNode.maxZ() + PW*hx;


                // no interpolation needed just copy. 
                if(pNodes[ele].getLevel()==bLev)
                {
                    const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/hh;
                    
                    for(unsigned int k=0; k < eOrder+1; k++)
                    {
                        double zz  = pNodes[ele].minZ() + k*hh;
                        
                        if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                        if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                        if(zz < zmin || zz > zmax) 
                            continue;
                        const int kkz = std::round((zz-zmin)*invhh);
                        assert( std::fabs(zz-zmin-kkz*hh) < d_compar_tol);
                        assert(kkz >= 0 && kkz < lz);

                        for(unsigned int j=0; j < eOrder+1; j++)
                        {   
                            double yy  = pNodes[ele].minY() + j*hh;

                            if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                            if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                            if(yy < ymin || yy > ymax) 
                                continue;
                            const int jjy = std::round((yy-ymin)*invhh);
                            //std::cout<<"yy: "<<yy<<" (ymin + hh*jjy): "<<(ymin + hh*jjy)<<std::endl;
                            assert( std::fabs(yy-ymin-jjy*hh) < d_compar_tol);
                            assert(jjy>=0 && jjy<ly);

                            for(unsigned int i=0; i < eOrder+1; i++)
                            {
                                double xx = pNodes[ele].minX() + i*hh;

                                if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                
                                if(xx < xmin || xx > xmax) 
                                    continue;
                                const int iix = std::round((xx-xmin)*invhh);
                                assert( std::fabs(xx-xmin-iix*hh) < d_compar_tol);
                                assert(iix>=0 && iix<lx);

                                //std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<kkz<<" , "<<jjy<<", "<<iix<<")"<<" of : "<<lx<<std::endl;

                                for(unsigned int v=0; v < dof; v++)
                                    uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                    

                            }
                        
                        }
                    
                    }


                }
                else if(pNodes[ele].getLevel() > bLev)
                {
                    assert((bLev+1) == pNodes[ele].getLevel());
                    const unsigned int cnum = pNodes[ele].getMortonIndex();
                    ot::TreeNode tmpParent = pNodes[ele].getParent();
                    
                    const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/(2*hh);

                    assert(eOrder>1);
                    const unsigned int cb =(eOrder%2==0) ? 0 : 1;

                    for(unsigned int k=cb; k < eOrder+1; k+=2)
                    {
                        double zz  = (pNodes[ele].minZ() + k*hh);
                        if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                        if(fabs(zz-zmax)<d_compar_tol) zz=zmax;

                        if(zz < zmin || zz > zmax) 
                            continue;
                        const int kkz = std::round((zz-zmin)*invhh);
                        assert(kkz >= 0 && kkz < lz);

                        for(unsigned int j=cb; j < eOrder+1; j+=2)
                        {   
                            double yy  = pNodes[ele].minY() + j*hh;
                            
                            if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                            if(fabs(yy-ymax)<d_compar_tol) yy=ymax;

                            if(yy < ymin || yy > ymax) 
                                continue;

                            const int jjy = std::round((yy-ymin)*invhh);
                            assert(jjy>=0 && jjy<ly);

                            for(unsigned int i=cb; i < eOrder+1; i+=2)
                            {
                                double xx = pNodes[ele].minX() + i*hh;

                                if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                
                                if(xx < xmin || xx > xmax) 
                                    continue;
                                const int iix = std::round((xx-xmin)*invhh);
                                assert(iix>=0 && iix<lx);

                                //std::cout<<"blk: "<<blk<<" blk copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;
                                for(unsigned int v=0; v < dof; v++)
                                    uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];
                                    
                                

                            }
                        
                        }
                    
                    }


                
                }
                else
                {
                    assert((bLev) == (pNodes[ele].getLevel()+1));
                    childOct.clear();
                    pNodes[ele].addChildren(childOct); // note this is the ordering of SFC (depends on Hilbert or Morton. )

                    for(unsigned int child = 0; child < NUM_CHILDREN; child++)
                    {

                        if( (childOct[child].maxX() < xmin || childOct[child].minX() >=xmax)  || (childOct[child].maxY() < ymin || childOct[child].minY() >=ymax) || (childOct[child].maxZ() < zmin || childOct[child].minZ() >=zmax) )
                            continue;


                        //std::cout<<"blk: "<<blk<<" blkNode: "<<blkNode<<" child: "<<child<<" child node "<<childOct[child]<<" parent : "<<pNodes[ele]<<std::endl;
                        const double hh = (1u<<(m_uiMaxDepth - childOct[child].getLevel()))/(double) eOrder;
                        const double invhh = 1.0/hh;

                        const unsigned int cnum = childOct[child].getMortonIndex();
                        if(!p2c_interp_valid[cnum])
                        {
                            for(unsigned int v=0; v < dof; v++)
                                this->parent2ChildInterpolation(&dgWVec[v*dgSz + ele*nPe], p2cI_all.data() + cnum * dof * nPe + v * nPe , cnum, m_uiDim);
                            
                            p2c_interp_valid[cnum]=true;
                        }
                        
                        for(unsigned int v=0; v < dof; v++)
                        {
                            const T* const p2cI = p2cI_all.data() + cnum * dof * nPe + v * nPe;
                            for(unsigned int k=0; k < eOrder+1; k++)
                            {
                                double zz  = childOct[child].minZ() + k*hh;
                                
                                if(fabs(zz-zmin)<d_compar_tol) zz=zmin;
                                if(fabs(zz-zmax)<d_compar_tol) zz=zmax;
                                
                                if(zz < zmin || zz > zmax) 
                                    continue;
                                const int kkz = std::round((zz-zmin)*invhh);
                                assert(kkz >= 0 && kkz < lz);

                                for(unsigned int j=0; j < eOrder+1; j++)
                                {   
                                    double yy  = childOct[child].minY() + j*hh;
                                    
                                    if(fabs(yy-ymin)<d_compar_tol) yy=ymin;
                                    if(fabs(yy-ymax)<d_compar_tol) yy=ymax;
                                    
                                    if(yy < ymin || yy > ymax) 
                                        continue;

                                    const int jjy = std::round((yy-ymin)*invhh);
                                    assert(jjy>=0 && jjy<ly);

                                    for(unsigned int i=0; i < eOrder+1; i++)
                                    {
                                        double xx = childOct[child].minX() + i*hh;
                                        
                                        if(fabs(xx-xmin)<d_compar_tol) xx=xmin;
                                        if(fabs(xx-xmax)<d_compar_tol) xx=xmax;
                                        
                                        if(xx < xmin || xx > xmax) 
                                            continue;
                                        const int iix = std::round((xx-xmin)*invhh);
                                        assert(iix>=0 && iix<lx);

                                        uzWVec[v*unSz + offset + kkz*lx*ly + jjy*lx + iix] =  p2cI[k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];

                                    }
                                
                                }
                            
                            }
                            
                        }

                    }

                }


            }

        }
        
    }

    template <typename T>
    void Mesh::unzipDG(const T *in, T *out, unsigned int dof)
    {
        if( (!m_uiIsActive) || (m_uiLocalBlockList.empty())  ) return;

        std::vector<unsigned int > blkIDs;
        blkIDs.resize(m_uiLocalBlockList.size());

        for(unsigned int i=0; i< m_uiLocalBlockList.size(); i++)
            blkIDs[i] = i ; 

        // unzip all the blocks. 
        this->unzipDG(in,out,blkIDs.data(),blkIDs.size(),dof);

    }

}//end of namespase ot





