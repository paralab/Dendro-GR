/**
 * @file subSM.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief sub scatter map
 * @version 0.1
 * @date 2020-01-17
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */


#include "subSM.h"

ot::SubScatterMap::SubScatterMap(const ot::Mesh* pMesh, std::bitset<ot::TreeNode::OctantFlagType::MAX_LEVEL>* stag, std::bitset<ot::TreeNode::OctantFlagType::MAX_LEVEL>* rtag)
{
    m_uiMesh=pMesh;
    m_uiSendTag = stag;
    m_uiRecvTag = rtag;
    m_uiMesh->computeMinMaxLevel(m_uiLMin , m_uiLMax);
    m_uiAsyncTag = 0;

    this->compute_L2SM();

}

ot::SubScatterMap::~SubScatterMap()
{

    if(m_uiIsAllocated)
    {

        //const unsigned int LMAX = ot::TreeNode::OctantFlagType::MAX_LEVEL;
        for(unsigned int l=m_uiLMin; l < m_uiLMax+1; l++ )
        {
            delete [] m_uiSendCount[l];  m_uiSendCount[l] = NULL;
            delete [] m_uiSendOffset[l]; m_uiSendOffset[l] = NULL;
            delete [] m_uiRecvCount[l];  m_uiRecvCount[l] = NULL;
            delete [] m_uiRecvOffset[l]; m_uiRecvOffset[l] = NULL;


            m_uiL2SendProcList[l].clear();
            m_uiL2RecvProcList[l].clear();

            delete [] m_uiL2SendSM[l]; m_uiL2SendSM[l] = NULL;
            delete [] m_uiL2RecvSM[l]; m_uiL2RecvSM[l] = NULL;

        }

        delete [] m_uiSendCount; m_uiSendCount = NULL;
        delete [] m_uiRecvCount; m_uiRecvCount = NULL;

        delete [] m_uiSendOffset; m_uiSendOffset = NULL; 
        delete [] m_uiRecvOffset; m_uiRecvOffset = NULL;

        delete [] m_uiL2SendProcList; m_uiL2SendProcList = NULL;
        delete [] m_uiL2RecvProcList; m_uiL2RecvProcList = NULL;

        delete [] m_uiL2SendSM; m_uiL2SendSM = NULL;
        delete [] m_uiL2RecvSM; m_uiL2RecvSM = NULL;

        m_uiIsAllocated = false;

    }

    
    

}

void ot::SubScatterMap::compute_L2SM()
{
    if( !(m_uiMesh->isActive()) || m_uiMesh->getMPICommSizeGlobal()==1 )
        return;

    unsigned int activeNpes = m_uiMesh->getMPICommSize();
    MPI_Comm commActive = m_uiMesh->getMPICommunicator();
    m_uiIsAllocated =true;

    const unsigned int LMAX = m_uiLMax+1 ;//ot::TreeNode::OctantFlagType::MAX_LEVEL;

    m_uiSendCount = new unsigned int * [LMAX];
    m_uiRecvCount = new unsigned int * [LMAX];

    m_uiSendOffset = new unsigned int * [LMAX];
    m_uiRecvOffset = new unsigned int * [LMAX];

    m_uiL2SendProcList = new std::vector<unsigned int > [LMAX];
    m_uiL2RecvProcList = new std::vector<unsigned int > [LMAX];

    m_uiL2SendSM = new unsigned int * [LMAX];
    m_uiL2RecvSM = new unsigned int * [LMAX];


    for(unsigned int l = m_uiLMin; l < m_uiLMax + 1 ; l++)
    {
        m_uiSendCount[l] = new unsigned int [activeNpes];
        m_uiRecvCount[l] = new unsigned int [activeNpes];

        m_uiSendOffset[l] = new unsigned int [activeNpes];
        m_uiRecvOffset[l] = new unsigned int [activeNpes];

        for(unsigned int k=0; k < activeNpes;k++ )
        {
            m_uiSendCount[l][k] = 0;
            m_uiRecvCount[l][k] = 0;
        }

    }

    //std::cout<<__LINE__<<" initial meme alloc "<<std::endl;
    const unsigned int * const sCounts = m_uiMesh->getNodalSendCounts().data();
    const unsigned int * const rCounts = m_uiMesh->getNodalRecvCounts().data();
    const unsigned int * const sOffset = m_uiMesh->getNodalSendOffsets().data();
    const unsigned int * const rOffset = m_uiMesh->getNodalRecvOffsets().data();
    const unsigned int * const sSM = m_uiMesh->getSendNodeSM().data();
    const unsigned int * const rSM = m_uiMesh->getRecvNodeSM().data();
    const std::vector<unsigned int>& sProcList = m_uiMesh->getSendProcList();
    const std::vector<unsigned int>& rProcList = m_uiMesh->getRecvProcList();


    for(unsigned int p =0; p < sProcList.size(); p++ )
    {
        const unsigned int proc_id = sProcList[p];
        for(unsigned int k = sOffset[proc_id]; k < (sOffset[proc_id] + sCounts[proc_id]); k++)
        {
            for(unsigned int l = m_uiLMin; l < m_uiLMax + 1 ; l++)
            {
                if(m_uiSendTag[k][l]==1)
                    m_uiSendCount[l][proc_id]++;
            }
        }
    }

    for(unsigned int p =0; p < rProcList.size(); p++ )
    {
        const unsigned int proc_id = rProcList[p];
        for(unsigned int k = rOffset[proc_id]; k < (rOffset[proc_id] + rCounts[proc_id]); k++)
        {
            for(unsigned int l = m_uiLMin; l < m_uiLMax + 1 ; l++)
            {
                if(m_uiRecvTag[k][l]==1)
                    m_uiRecvCount[l][proc_id]++;
            }
        }
    }

    #ifdef __DEBUG_SUB_SM_CODE__
        
        std::vector<unsigned int> rCts; rCts.resize(activeNpes);
        for(unsigned int l= m_uiLMin; l < m_uiLMax +1; l++  )
        {
            par::Mpi_Alltoall(m_uiSendCount[l],rCts.data(),1,commActive);
            for(unsigned int p =0; p< activeNpes; p++)
                if(rCts[p] != m_uiRecvCount[l][p] )
                {
                    std::cout<<" rank : "<<  m_uiMesh->getMPIRank() <<" p : "<<p <<" recvCounts: "<<m_uiRecvCount[l][p]<<" rCnt "<<rCts[p]<<"\n";
                    std::ostringstream oss;
                    oss <<__func__<<" inconsistent send recv counts for level " << l << " encountered. ";
                    print_error(oss.str().c_str());
                    MPI_Abort(commActive,0);
                }

        }
        rCts.clear();
        for(unsigned int l = m_uiLMin; l < m_uiLMax + 1 ; l++)
        {
            io::dump_array(m_uiSendCount[l],activeNpes,"sendCnt",m_uiMesh->getMPICommunicator());
            io::dump_array(m_uiRecvCount[l],activeNpes,"recvCnt",m_uiMesh->getMPICommunicator());
        }

    #endif


    
    
    //std::cout<<"send recv counts computed. \n";

    for(unsigned int l= m_uiLMin; l < m_uiLMax +1; l++  )
    {

        m_uiSendOffset[l][0]=0;
        m_uiRecvOffset[l][0]=0;

        omp_par::scan(m_uiSendCount[l],m_uiSendOffset[l],activeNpes);
        omp_par::scan(m_uiRecvCount[l],m_uiRecvOffset[l],activeNpes);

        for(unsigned int p =0; p < activeNpes; p++)
        {
            if( m_uiSendCount[l][p] >0 )
                m_uiL2SendProcList[l].push_back(p);

            if( m_uiRecvCount[l][p] >0 )
                m_uiL2RecvProcList[l].push_back(p);

        }

        const unsigned int sendBufSz = (m_uiSendOffset[l][activeNpes-1] + m_uiSendCount[l][activeNpes-1]);
        const unsigned int recvBufSz = (m_uiRecvOffset[l][activeNpes-1] + m_uiRecvCount[l][activeNpes-1]);
        
        m_uiL2SendSM[l] = new unsigned int [sendBufSz];
        m_uiL2RecvSM[l] = new unsigned int [recvBufSz];


    }

    //std::cout<<" l2 SM allocated\n";

    std::vector<unsigned int> scounts;
    std::vector<unsigned int> rcounts;

    scounts.resize(LMAX,0);
    rcounts.resize(LMAX,0);

    for(unsigned int p =0; p < sProcList.size(); p++ )
    {
        const unsigned int proc_id = sProcList[p];
        for(unsigned int k = sOffset[proc_id]; k < (sOffset[proc_id] + sCounts[proc_id]); k++)
        {
            for(unsigned int l = m_uiLMin; l < m_uiLMax + 1 ; l++)
            {
                if(m_uiSendTag[k][l]==1)
                {
                    m_uiL2SendSM[l][scounts[l]] = sSM[k];
                    scounts[l]++;
                }
                    
                    
            }
        }
    }

    for(unsigned int p =0; p < rProcList.size(); p++ )
    {
        const unsigned int proc_id = rProcList[p];
        for(unsigned int k = rOffset[proc_id]; k < (rOffset[proc_id] + rCounts[proc_id]); k++)
        {
            for(unsigned int l = m_uiLMin; l < m_uiLMax + 1 ; l++)
            {
                if(m_uiRecvTag[k][l]==1)
                {
                    m_uiL2RecvSM[l][rcounts[l]] = rSM[k];
                    rcounts[l]++;
                }
                    
            }
        }
    }

    //std::cout<<"sm build complete \n"<<std::endl;

    scounts.clear();
    rcounts.clear();

    return ;


}
