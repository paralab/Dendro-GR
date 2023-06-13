/***
 * @brief: sub oda cpp file
 * @author: Hari Sundar
 * @author: Milinda Fernando
 * 
 */

#include "sub_oda.h"


ot::subDA::subDA(ot::DA* da, std::function<double ( double, double, double ) > fx_retain, double* gSize)
{
    m_da = da;
    
    m_uiCommGlobal = m_da->getGlobalComm();
    m_uiNpesGlobal = m_da->getNpesAll();
    m_uiRankActive = m_da->getRankAll();

    unsigned int lev;
    const unsigned int maxDepth = m_uiMaxDepth;
    const unsigned int allEleSize  = m_da->getLocalElemSz() + m_da->getPreGhostElementSize() + m_da->getPostGhostElementSize();
    const unsigned int allNodeSize = m_da->getLocalNodalSz() + m_da->getPreAndPostGhostNodeSize();
    const ot::Mesh * _mesh = m_da->getMesh();
    
    const unsigned int nodeLocalBegin = _mesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd   = _mesh->getNodeLocalEnd(); 

    m_uiCommTag = 0;
    m_mpiContexts.clear();

    const unsigned int nPe = m_da->getNumNodesPerElement();
    m_uiIsActive = false;

    if(!m_da->isActive())
    {
        m_uiElementPreGhostBegin=0;
        m_uiElementPreGhostEnd=0;
        m_uiElementLocalBegin=0;
        m_uiElementLocalEnd=0;
        m_uiElementPostGhostBegin=0;
        m_uiElementPostGhostEnd=0;

        m_uiNumPreGhostElements=m_uiElementPreGhostEnd-m_uiElementPreGhostBegin;
        m_uiNumLocalElements=m_uiElementLocalEnd=m_uiElementLocalBegin;
        m_uiNumPostGhostElements=m_uiElementPostGhostEnd-m_uiElementPostGhostBegin;

        m_uiNodePreGhostBegin=0;
        m_uiNodePreGhostEnd=0;
        m_uiNodeLocalBegin=0;
        m_uiNodeLocalEnd=0;
        m_uiNodePostGhostBegin=0;
        m_uiNodePostGhostEnd=0;

        

    }else
    {
        m_uiIsActive =true;
        
        m_uiCommActive = m_da->getCommActive();
        m_uiNpesActive = m_da->getNpesActive();
        m_uiRankActive = m_da->getRankActive();

        MPI_Comm commActive = m_da->getCommActive();
        unsigned int npesActive = m_da->getNpesActive();

        auto inside = [](double d){ return d < 0.0; };
        double hx, hy, hz;
        std::vector<double> dist;
        dist.resize(nPe);
        
        Point pt;

        double xFac = gSize[0]/((double)(1<<(maxDepth)));
        double yFac = gSize[1]/((double)(1<<(maxDepth)));
        double zFac = gSize[2]/((double)(1<<(maxDepth)));

        m_ucpSkipList.clear();
        m_ucpSkipList.resize(allEleSize , 0);

        m_ucpSkipNodeList.clear();
        m_ucpSkipNodeList.resize(allNodeSize, 0);

        m_uip_DA2sub_ElemMap.resize(allEleSize, 0);
        m_uip_DA2sub_NodeMap.resize(allNodeSize, 0);
        
        for (unsigned int i=0; i<allNodeSize; ++i) {
            m_ucpSkipNodeList[i] = 1;
        }
        
        double _min[3], _max[3];

        _min[0] = gSize[0]; _min[1] = gSize[1]; _min[2] = gSize[2];
        _max[0] = 0.0; _max[1] = 0.0;  _max[2] = 0.0;

        double x, y, z;
        unsigned int num_local_skip = 0;
        
        const unsigned int nodeLocalBegin=_mesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=_mesh->getNodeLocalEnd();
        unsigned int eleOrder=_mesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(_mesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(_mesh->getE2NMapping_DG().begin()));
        const ot::TreeNode * pNodes = &(*(_mesh->getAllElements().begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z,nodeLookUp_DG;

        std::vector<double> coords;
        std::vector<DendroIntL> indices;
        indices.resize(nPe);
        coords.resize(nPe*m_uiDim);
                
        for ( m_da->init<ot::DA_FLAGS::WRITABLE>(); m_da->curr() < m_da->end<ot::DA_FLAGS::WRITABLE>(); m_da->next<ot::DA_FLAGS::WRITABLE>() )
        {
            lev = m_da->getLevel(m_da->curr());
            for(unsigned int k=0; k<(eleOrder+1); k++)
                for(unsigned int j=0; j<(eleOrder+1); j++ )
                    for(unsigned int i=0; i<(eleOrder+1); i++)
                    {
                        nodeLookUp_DG=e2n_dg[(m_da->curr())*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        _mesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                        const unsigned int len=1u<<(m_uiMaxDepth-lev);
                        coords[m_uiDim*(k*(eleOrder+1)*(eleOrder+1) + j*(eleOrder+1) +i) + 0 ] = xFac*(pNodes[ownerID].getX()+ ii_x*(len/(eleOrder)));
                        coords[m_uiDim*(k*(eleOrder+1)*(eleOrder+1) + j*(eleOrder+1) +i) + 1 ] = yFac*(pNodes[ownerID].getY()+ jj_y*(len/(eleOrder)));
                        coords[m_uiDim*(k*(eleOrder+1)*(eleOrder+1) + j*(eleOrder+1) +i) + 2 ] = zFac*(pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder)));
                        x = coords[m_uiDim*(k*(eleOrder+1)*(eleOrder+1) + j*(eleOrder+1) +i) + 0 ];
                        y = coords[m_uiDim*(k*(eleOrder+1)*(eleOrder+1) + j*(eleOrder+1) +i) + 1 ];
                        z = coords[m_uiDim*(k*(eleOrder+1)*(eleOrder+1) + j*(eleOrder+1) +i) + 2 ];
                        dist[k*(eleOrder+1)*(eleOrder+1) + j*(eleOrder+1) +i] = fx_retain(x,y,z);

                    }
            if(std::all_of(dist.begin(), dist.end(), inside ))
            {
                m_ucpSkipList[m_da->curr()]=1;
            }else
            {
                    x = coords[m_uiDim*(0) + 0 ];
                    y = coords[m_uiDim*(0) + 1 ];
                    z = coords[m_uiDim*(0) + 2 ];

                if ( x < _min[0] ) _min[0] = x;
                if ( y < _min[1] ) _min[1] = y;
                if ( z < _min[2] ) _min[2] = z;

                    x = coords[m_uiDim*(nPe-1) + 0 ];
                    y = coords[m_uiDim*(nPe-1) + 1 ];
                    z = coords[m_uiDim*(nPe-1) + 2 ];

                if ( x > _max[0] ) _max[0] = x;
                if ( y > _max[1] ) _max[1] = y;
                if ( z > _max[2] ) _max[2] = z;

                m_da->getNodeIndices(&(*(indices.begin())), m_da->curr(), true );
                for(unsigned int node =0; node < nPe ; node ++)
                {
                    if(nodeLocalBegin<=indices[node]  && indices[node]<nodeLocalEnd )
                        m_ucpSkipNodeList[indices[node]] =0;
                }
                
            }
        }

        par::Mpi_Allreduce(_min, m_dMinBB , 3 , MPI_MIN , commActive);
        par::Mpi_Allreduce(_max, m_dMaxBB , 3 , MPI_MAX , commActive);

        m_da->readFromGhostBegin(&(*(m_ucpSkipNodeList.begin())),1);
        m_da->readFromGhostEnd(&(*(m_ucpSkipNodeList.begin())));
    
        unsigned int sum =0;
        for(unsigned int i=0; i < m_uip_DA2sub_ElemMap.size(); i++ )
        {   
            m_uip_DA2sub_ElemMap[i] = sum;
            if (m_ucpSkipList[i] == 0) sum++;
        }

        m_uiTotalElementSz = sum;

        sum = 0;
        for (unsigned int i=0; i<m_uip_DA2sub_NodeMap.size(); ++i) {
            m_uip_DA2sub_NodeMap[i] = sum;
            if (m_ucpSkipNodeList[i] == 0) sum++;
        }

        m_uiTotalNodalSz=sum;

        if(m_uiTotalNodalSz ==0)
        {
            std::cout<<"subDA[Error] on processor "<<m_uiRankActive<<": empty subDA of a active DA (use weighted partition to avoid this) "<<__func__<<std::endl;
        }
        
        m_uip_sub2DA_ElemMap.resize(m_uiTotalElementSz);            
        m_uip_sub2DA_NodeMap.resize(m_uiTotalNodalSz);
        
        m_uiNumPreGhostElements = 0;
        for(unsigned int e=_mesh->getElementPreGhostBegin(); e < _mesh->getElementPreGhostEnd(); e++)
        {
            if(m_ucpSkipList[e]==0)
            {
                m_uip_sub2DA_ElemMap[m_uiNumPreGhostElements] = e;
                m_uiNumPreGhostElements++;
            }
        }

        m_uiNumLocalElements = 0;
        for(unsigned int e=_mesh->getElementLocalBegin(); e < _mesh->getElementLocalEnd() ; e++ )
        {
            if(m_ucpSkipList[e] ==0)
            {
                m_uip_sub2DA_ElemMap[ (m_uiNumPreGhostElements) +  m_uiNumLocalElements ] = e ;
                m_uiNumLocalElements++;
            }
        }

        m_uiNumPostGhostElements = 0;
        for(unsigned int e = _mesh->getElementPostGhostBegin(); e < _mesh->getElementPostGhostEnd() ; e++)
        {
            if(m_ucpSkipList[e] == 0)
            {
                m_uip_sub2DA_ElemMap[(m_uiNumPreGhostElements +  m_uiNumLocalElements) +  m_uiNumPostGhostElements] =  e;
                m_uiNumPostGhostElements++;
            }
        }


        m_uiNumPreGhostNodes = 0;
        for(unsigned int e=_mesh->getNodePreGhostBegin(); e < _mesh->getNodePreGhostEnd(); e++)
        {
            if(m_ucpSkipNodeList[e] == 0 )
            {
                m_uip_sub2DA_NodeMap[m_uiNumPreGhostNodes] = e;
                m_uiNumPreGhostNodes++;
            }
        }

        m_uiNumLocalNodes = 0;
        for(unsigned int e=_mesh->getNodeLocalBegin() ; e< _mesh->getNodeLocalEnd(); e++)
        {
            if(m_ucpSkipNodeList[e] == 0 )
            {
                m_uip_sub2DA_NodeMap[(m_uiNumPreGhostNodes) + m_uiNumLocalNodes] = e ;
                m_uiNumLocalNodes++;
            }
        }

        m_uiNumPostGhostNodes = 0;    
        for(unsigned int e =_mesh->getNodePostGhostBegin(); e < _mesh->getNodePostGhostEnd(); e++)
        {
            if(m_ucpSkipNodeList[e] ==0)
            {
                m_uip_sub2DA_NodeMap[(m_uiNumPreGhostNodes + m_uiNumLocalNodes) + m_uiNumPostGhostNodes] = e;
                m_uiNumPostGhostNodes++;
            }
        }


        m_uiElementPreGhostBegin = 0;
        m_uiElementPreGhostEnd = m_uiNumPreGhostElements;

        m_uiElementLocalBegin = m_uiNumPreGhostElements;
        m_uiElementLocalEnd = m_uiNumPreGhostElements + m_uiNumLocalElements;

        m_uiElementPostGhostBegin = m_uiNumPreGhostElements + m_uiNumLocalElements;
        m_uiElementPostGhostEnd = m_uiNumPreGhostElements + m_uiNumLocalElements + m_uiNumPostGhostElements ;


        m_uiNodePreGhostBegin = 0;
        m_uiNodePreGhostEnd = m_uiNumPreGhostNodes;

        m_uiNodeLocalBegin = m_uiNumPreGhostNodes;
        m_uiNodeLocalEnd = m_uiNumPreGhostNodes + m_uiNumLocalNodes;

        m_uiNodePostGhostBegin = m_uiNumPreGhostNodes + m_uiNumLocalNodes;
        m_uiNodePostGhostEnd = m_uiNumPreGhostNodes + m_uiNumLocalNodes + m_uiNumPostGhostNodes;


        m_uiSendCounts.resize(m_uiNpesActive);
        m_uiSendOffsets.resize(m_uiNpesActive);

        m_uiRecvCounts.resize(m_uiNpesActive);
        m_uiRecvOffsets.resize(m_uiNpesActive);

        unsigned int k=0; 
        unsigned int offset=0;
        const std::vector<unsigned int>& da_sendCount = _mesh->getNodalSendCounts();
        const std::vector<unsigned int>& da_recvCount = _mesh->getNodalRecvCounts();
        const std::vector<unsigned int>& da_sendOffset = _mesh->getNodalSendOffsets();
        const std::vector<unsigned int>& da_recvOffset = _mesh->getNodalRecvOffsets();

        const std::vector<unsigned int >& spList = _mesh->getSendProcList();
        const std::vector<unsigned int >& rpList = _mesh->getRecvProcList();

        const std::vector<unsigned int >& da_sm = _mesh->getSendNodeSM();
        const std::vector<unsigned int >& da_rm = _mesh->getRecvNodeSM();

        for (unsigned int p=0; p<spList.size(); ++p) {
            unsigned int cnt=0;  
            for(unsigned int w = (da_sendOffset[spList[p]]); w < (da_sendOffset[spList[p]] + da_sendCount[spList[p]]) ; w++ )
            {
                unsigned int idx = da_sm[w];
                if (m_ucpSkipNodeList[idx] == 0 )
                {
                    m_uiSendScatterMap.push_back(m_uip_DA2sub_NodeMap[idx]);
                    cnt++;
                }
            }

            if(cnt)
              m_uiSendProcList.push_back(spList[p]);

              m_uiSendCounts[spList[p]] = cnt;
            }

        for (unsigned int p=0; p<rpList.size(); ++p) {
            unsigned int cnt=0;  
            for(unsigned int w = (da_recvOffset[spList[p]]); w < (da_recvOffset[spList[p]] + da_recvCount[spList[p]]) ; w++ )
            {
                unsigned int idx = da_rm[w];
                if (m_ucpSkipNodeList[idx] == 0 )
                {
                  m_uiRecvScatterMap.push_back(m_uip_DA2sub_NodeMap[idx]);
                  cnt++;
                }
            }
            
            if(cnt)
              m_uiRecvProcList.push_back(rpList[p]);

            m_uiRecvCounts[rpList[p]] = cnt;
        }

        m_uiSendOffsets[0] = 0;
        m_uiRecvOffsets[0] = 0;

        omp_par::scan(&(*(m_uiSendCounts.begin())), &(*(m_uiSendCounts.begin())),npesActive);
        omp_par::scan(&(*(m_uiRecvCounts.begin())), &(*(m_uiRecvOffsets.begin())),npesActive);
        

    }





        
        
}

ot::subDA::~subDA()
{
    
}

#ifdef BUILD_WITH_PETSC
    PetscErrorCode ot::subDA::petscCreateVector(Vec &local, bool isElemental, bool isGhosted, unsigned int dof) const
    {
        unsigned int sz=0;
        MPI_Comm globalComm=this->getGlobalComm();
        if(!(m_da->isActive()))
        {
            local=NULL;

        }else {

            if(isElemental)
            {
                if(isGhosted)
                    sz=dof*m_uiTotalElementSz;
                else
                    sz=dof*m_uiLocalElementSz;

            }else {

                if(isGhosted)
                    sz=dof*m_uiTotalNodalSz;
                else
                    sz=dof*m_uiLocalNodalSz;
            }

        }

        VecCreate(globalComm,&local);
        PetscErrorCode status=VecSetSizes(local,sz,PETSC_DECIDE);

        if (this->getNpesAll() > 1) {
            VecSetType(local,VECMPI);
        } else {
            VecSetType(local,VECSEQ);
        }


        return status;


    }
#endif

