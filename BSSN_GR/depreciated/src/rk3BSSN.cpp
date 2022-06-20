//
// Created by milinda on 12/6/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief rk3 solver for bssn
*/
//


//
// Created by milinda on 12/1/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "../include/rk3BSSN.h"

namespace ode
{
    namespace solver
    {

        RK3_BSSN::RK3_BSSN(ot::Mesh *pMesh, double pTBegin, double pTEnd,double pTh): RK(pMesh,pTBegin,pTEnd,pTh)
        {

            // allocate memory for the variables.
            m_uiVar=new double*[bssn::BSSN_NUM_VARS];
            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiVar[index]=m_uiMesh->createVector<double>();

            m_uiPrevVar=new double*[bssn::BSSN_NUM_VARS];
            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiPrevVar[index]=m_uiMesh->createVector<double>();

            m_uiVarIm=new double*[bssn::BSSN_NUM_VARS];
            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiVarIm[index]=m_uiMesh->createVector<double>();


            m_uiStage=new double**[bssn::BSSN_RK3_STAGES];
            for(unsigned int stage=0;stage<bssn::BSSN_RK3_STAGES;stage++)
            {
                m_uiStage[stage]=new double*[bssn::BSSN_NUM_VARS];
                for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                    m_uiStage[stage][index]=m_uiMesh->createVector<double>();
            }

            m_uiUnzipVar=new double*[bssn::BSSN_NUM_VARS];
            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiUnzipVar[index]=m_uiMesh->createUnZippedVector<double>();


            m_uiUnzipVarRHS=new double*[bssn::BSSN_NUM_VARS];
            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiUnzipVarRHS[index]=m_uiMesh->createUnZippedVector<double>();



            // allocate memory for the constraint variables.
            m_uiConstraintVars=new double*[bssn::BSSN_CONSTRAINT_NUM_VARS];
            for(unsigned int index=0;index<bssn::BSSN_CONSTRAINT_NUM_VARS;index++)
                m_uiConstraintVars[index]=m_uiMesh->createVector<double>();


            m_uiUnzipConstraintVars=new double*[bssn::BSSN_CONSTRAINT_NUM_VARS];
            for(unsigned int index=0;index<bssn::BSSN_CONSTRAINT_NUM_VARS;index++)
                m_uiUnzipConstraintVars[index]=m_uiMesh->createUnZippedVector<double>();


            // mpi communication
            m_uiSendNodeBuf=new double*[bssn::BSSN_ASYNC_COMM_K];
            m_uiRecvNodeBuf=new double*[bssn::BSSN_ASYNC_COMM_K];

            m_uiSendReqs=new MPI_Request*[bssn::BSSN_ASYNC_COMM_K];
            m_uiRecvReqs=new MPI_Request*[bssn::BSSN_ASYNC_COMM_K];
            m_uiSendSts=new MPI_Status*[bssn::BSSN_ASYNC_COMM_K];
            m_uiRecvSts=new MPI_Status*[bssn::BSSN_ASYNC_COMM_K];

            for(unsigned int index=0;index<bssn::BSSN_ASYNC_COMM_K;index++)
            {
                m_uiSendNodeBuf[index]=NULL;
                m_uiRecvNodeBuf[index]=NULL;

                m_uiSendReqs[index]=NULL;
                m_uiRecvReqs[index]=NULL;
                m_uiSendSts[index]=NULL;
                m_uiRecvSts[index]=NULL;
            }

            if(m_uiMesh->isActive())
            {
                // allocate mpi comm. reqs and status
                for(unsigned int index=0;index<bssn::BSSN_ASYNC_COMM_K;index++)
                {
                    if(m_uiMesh->getGhostExcgTotalSendNodeCount()!=0) m_uiSendNodeBuf[index]=new double[m_uiMesh->getGhostExcgTotalSendNodeCount()];
                    if(m_uiMesh->getGhostExcgTotalRecvNodeCount()!=0) m_uiRecvNodeBuf[index]=new double[m_uiMesh->getGhostExcgTotalRecvNodeCount()];

                    if(m_uiMesh->getSendProcListSize()!=0)
                    {
                        m_uiSendReqs[index]=new MPI_Request[m_uiMesh->getSendProcListSize()];
                        m_uiSendSts[index]=new MPI_Status[m_uiMesh->getSendProcListSize()];
                    }

                    if(m_uiMesh->getRecvProcListSize()!=0)
                    {
                        m_uiRecvReqs[index]=new MPI_Request[m_uiMesh->getRecvProcListSize()];
                        m_uiRecvSts[index]=new MPI_Status[m_uiMesh->getRecvProcListSize()];
                    }


                }

            }


        }

        RK3_BSSN::~RK3_BSSN()
        {
            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
            {
                delete [] m_uiVar[index];
                delete [] m_uiPrevVar[index];
                delete [] m_uiVarIm[index];
                delete [] m_uiUnzipVar[index];
                delete [] m_uiUnzipVarRHS[index];


            }


            delete [] m_uiVar;
            delete [] m_uiPrevVar;
            delete [] m_uiVarIm;
            delete [] m_uiUnzipVar;
            delete [] m_uiUnzipVarRHS;

            for(unsigned int stage=0;stage<bssn::BSSN_RK3_STAGES;stage++)
                for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                    delete [] m_uiStage[stage][index];

            for(unsigned int stage=0;stage<bssn::BSSN_RK3_STAGES;stage++)
                delete [] m_uiStage[stage];

            delete [] m_uiStage;


            // deallocate memory for the constraint variables.
            for(unsigned int index=0;index<bssn::BSSN_CONSTRAINT_NUM_VARS;index++)
            {
                delete [] m_uiConstraintVars[index];
                delete [] m_uiUnzipConstraintVars[index];
            }


            delete [] m_uiConstraintVars;
            delete [] m_uiUnzipConstraintVars;


            // mpi communication
            for(unsigned int index=0;index<bssn::BSSN_ASYNC_COMM_K;index++)
            {
                delete [] m_uiSendNodeBuf[index];
                delete [] m_uiRecvNodeBuf[index];

                delete [] m_uiSendReqs[index];
                delete [] m_uiRecvReqs[index];

                delete [] m_uiSendSts[index];
                delete [] m_uiRecvSts[index];

            }

            delete [] m_uiSendNodeBuf;
            delete [] m_uiRecvNodeBuf;

            delete [] m_uiSendReqs;
            delete [] m_uiSendSts;
            delete [] m_uiRecvReqs;
            delete [] m_uiRecvSts;



        }

        void RK3_BSSN::applyInitialConditions(double** zipIn)
        {
            unsigned int nodeLookUp_CG;
            unsigned int nodeLookUp_DG;
            unsigned int x,y,z,len;
            const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
            unsigned int ownerID,ii_x,jj_y,kk_z;
            unsigned int eleOrder=m_uiMesh->getElementOrder();
            const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
            const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
            const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
            const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();


            double* var=new double[bssn::BSSN_NUM_VARS];

            double mp, mm, mp_adm, mm_adm, E, J1, J2, J3;

            for(unsigned int elem=m_uiMesh->getElementLocalBegin();elem<m_uiMesh->getElementLocalEnd();elem++)
            {


                for(unsigned int k=0;k<(eleOrder+1);k++)
                    for(unsigned int j=0;j<(eleOrder+1);j++ )
                        for(unsigned int i=0;i<(eleOrder+1);i++)
                        {
                            nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                            {
                                nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                                len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                                x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                                y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                                z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));
                                assert(len%eleOrder==0);
                                if (bssn::BSSN_ID_TYPE == 0) {
                                    TwoPunctures((double)x,(double)y,(double)z,var,
                                                        &mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
                                }
                                else if (bssn::BSSN_ID_TYPE == 1) {
                                    bssn::punctureData((double)x,(double)y,(double)z,var);
                                } 
				                else if (bssn::BSSN_ID_TYPE == 2) {
                               	    bssn::KerrSchildData((double)x,(double)y,(double)z,var);
                                }
				                else if (bssn::BSSN_ID_TYPE == 3) {
                               	    bssn::noiseData((double)x,(double)y,(double)z,var);
                                }
				                else {
				                    std::cout<<"Unknown ID type"<<std::endl;
				                }
                                for(unsigned int v=0;v<bssn::BSSN_NUM_VARS;v++)
                                    zipIn[v][nodeLookUp_CG]=var[v];


                            }

                        }

            }
            
            
            for(unsigned int node=m_uiMesh->getNodeLocalBegin();node<m_uiMesh->getNodeLocalEnd();node++)
                enforce_bssn_constraints(zipIn,node);
            

            delete [] var;

        }

        void RK3_BSSN::initialGridConverge()
        {

            applyInitialConditions(m_uiPrevVar);

            bool isRefine=false;
            unsigned int oldElements,oldElements_g;
            unsigned int newElements,newElements_g;

            // refine based on all the variables
            const unsigned int refineNumVars=bssn::BSSN_NUM_REFINE_VARS;
            unsigned int refineVarIds[refineNumVars];
            for(unsigned int vIndex=0;vIndex<refineNumVars;vIndex++)
                refineVarIds[vIndex]=bssn::BSSN_REFINE_VARIABLE_INDICES[vIndex];

            double wTol=bssn::BSSN_WAVELET_TOL;
            std::function<double(double,double,double)> waveletTolFunc =[wTol](double x,double y, double z){ return bssn::computeWTol(x,y,z,wTol);};
            unsigned int iterCount=1;
            do
            {

#ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
#else
                performGhostExchangeVars(m_uiPrevVar);
                unzipVars(m_uiPrevVar,m_uiUnzipVar);
#endif


                if(bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
                    isRefine=false;
                else
                    isRefine=m_uiMesh->isReMeshUnzip((const double **)m_uiUnzipVar,refineVarIds,refineNumVars,waveletTolFunc,bssn::BSSN_DENDRO_AMR_FAC);

                if(isRefine)
                {
                    ot::Mesh* newMesh=m_uiMesh->ReMesh(bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);

                    oldElements=m_uiMesh->getNumLocalMeshElements();
                    newElements=newMesh->getNumLocalMeshElements();

                    par::Mpi_Allreduce(&oldElements,&oldElements_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());
                    par::Mpi_Allreduce(&newElements,&newElements_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());

                    if(!(m_uiMesh->getMPIRankGlobal()))std::cout<<"initial grid iteration : "<<iterCount<<" old mesh: "<<oldElements_g<<" new mesh: "<<newElements_g<<std::endl;


                    // performs the inter-grid transfer
                    intergridTransferVars(m_uiPrevVar,newMesh);

                    for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                    {
                        delete [] m_uiVar[index];
                        delete [] m_uiVarIm[index];
                        delete [] m_uiUnzipVar[index];
                        delete [] m_uiUnzipVarRHS[index];

                        m_uiVar[index]=NULL;
                        m_uiVarIm[index]=NULL;
                        m_uiUnzipVar[index]=NULL;
                        m_uiUnzipVarRHS[index]=NULL;

                        m_uiVar[index]=newMesh->createVector<double>();
                        m_uiVarIm[index]=newMesh->createVector<double>();
                        m_uiUnzipVar[index]=newMesh->createUnZippedVector<double>();
                        m_uiUnzipVarRHS[index]=newMesh->createUnZippedVector<double>();


                    }

                    for(unsigned int stage=0;stage<bssn::BSSN_RK4_STAGES;stage++)
                        for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        {
                            delete [] m_uiStage[stage][index];
                            m_uiStage[stage][index]=NULL;
                            m_uiStage[stage][index]=newMesh->createVector<double>();
                        }

                    // deallocate constraint vars allocate them for the new mesh.
                    for(unsigned int index=0;index<bssn::BSSN_CONSTRAINT_NUM_VARS;index++)
                    {
                        delete [] m_uiConstraintVars[index];
                        delete [] m_uiUnzipConstraintVars[index];

                        m_uiConstraintVars[index]=newMesh->createVector<double>();
                        m_uiUnzipConstraintVars[index]=newMesh->createUnZippedVector<double>();

                    }




                    std::swap(newMesh,m_uiMesh);
                    delete newMesh;

#ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                    // reallocates mpi resources for the the new mesh. (this will deallocate the old resources)
                    reallocateMPIResources();
#endif

                    if(m_uiMesh->isActive())
                    {
                        double l_min=vecMin(m_uiPrevVar[bssn::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                        double l_max=vecMax(m_uiPrevVar[bssn::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                        if(!(m_uiMesh->getMPIRank())) {
                            std::cout << "transfer completed:    ||VAR::U_ALPHA|| (min, max) : ("<<l_min<<", "<<l_max<<" ) "<<std::endl;
                        }

                    }

                    iterCount+=1;

                }

                applyInitialConditions(m_uiPrevVar);

            }while(isRefine && (newElements_g!=oldElements_g));


        }

        void RK3_BSSN::reallocateMPIResources()
        {
            for(unsigned int index=0;index<bssn::BSSN_ASYNC_COMM_K;index++)
            {
                delete [] m_uiSendNodeBuf[index];
                delete [] m_uiRecvNodeBuf[index];

                delete [] m_uiSendReqs[index];
                delete [] m_uiRecvReqs[index];

                delete [] m_uiSendSts[index];
                delete [] m_uiRecvSts[index];

            }

            for(unsigned int index=0;index<bssn::BSSN_ASYNC_COMM_K;index++)
            {
                m_uiSendNodeBuf[index]=NULL;
                m_uiRecvNodeBuf[index]=NULL;

                m_uiSendReqs[index]=NULL;
                m_uiRecvReqs[index]=NULL;
                m_uiSendSts[index]=NULL;
                m_uiRecvSts[index]=NULL;
            }

            if(m_uiMesh->isActive())
            {
                // allocate mpi comm. reqs and status
                for(unsigned int index=0;index<bssn::BSSN_ASYNC_COMM_K;index++)
                {
                    if(m_uiMesh->getGhostExcgTotalSendNodeCount()!=0) m_uiSendNodeBuf[index]=new double[m_uiMesh->getGhostExcgTotalSendNodeCount()];
                    if(m_uiMesh->getGhostExcgTotalRecvNodeCount()!=0) m_uiRecvNodeBuf[index]=new double[m_uiMesh->getGhostExcgTotalRecvNodeCount()];

                    if(m_uiMesh->getSendProcListSize()!=0)
                    {
                        m_uiSendReqs[index]=new MPI_Request[m_uiMesh->getSendProcListSize()];
                        m_uiSendSts[index]=new MPI_Status[m_uiMesh->getSendProcListSize()];
                    }

                    if(m_uiMesh->getRecvProcListSize()!=0)
                    {
                        m_uiRecvReqs[index]=new MPI_Request[m_uiMesh->getRecvProcListSize()];
                        m_uiRecvSts[index]=new MPI_Status[m_uiMesh->getRecvProcListSize()];
                    }


                }

            }


        }

        void RK3_BSSN::writeToVTU(double **evolZipVarIn, double ** constrZipVarIn, unsigned int numEvolVars,unsigned int numConstVars,const unsigned int * evolVarIndices, const unsigned int * constVarIndices)
        {
            bssn::timer::t_ioVtu.start();

            std::vector<std::string> pDataNames;
            double *pData[(numConstVars+numEvolVars)];

            for(unsigned int i=0;i<numEvolVars;i++)
            {
                pDataNames.push_back(std::string(bssn::BSSN_VAR_NAMES[evolVarIndices[i]]));
                pData[i]=evolZipVarIn[evolVarIndices[i]];
            }


            for(unsigned int i=0;i<numConstVars;i++)
            {
                pDataNames.push_back(std::string(bssn::BSSN_CONSTRAINT_VAR_NAMES[constVarIndices[i]]));
                pData[numEvolVars+i]=constrZipVarIn[constVarIndices[i]];
            }

            std::vector<char*> pDataNames_char;
            pDataNames_char.reserve(pDataNames.size());

            for(unsigned int  i = 0; i < pDataNames.size(); i++)
                pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

            const char * fDataNames[]={"Time","Cycle"};
            const double fData[]={m_uiCurrentTime,m_uiCurrentStep};

            char fPrefix[256];
            sprintf(fPrefix,"%s_%d",bssn::BSSN_VTU_FILE_PREFIX.c_str(),m_uiCurrentStep);

            io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);

            bssn::timer::t_ioVtu.stop();

        }

        void RK3_BSSN::performGhostExchangeVars(double** zipIn)
        {

            bssn::timer::t_ghostEx_sync.start();

            for(unsigned int v=0;v<bssn::BSSN_NUM_VARS;v++)
                m_uiMesh->performGhostExchange(zipIn[v]);

            bssn::timer::t_ghostEx_sync.stop();

        }

        void RK3_BSSN::intergridTransferVars(double**& zipIn, const ot::Mesh* pnewMesh)
        {
            bssn::timer::t_gridTransfer.start();

            for(unsigned int v=0;v<bssn::BSSN_NUM_VARS;v++)
                m_uiMesh->interGridTransfer(zipIn[v],pnewMesh);

            bssn::timer::t_gridTransfer.stop();

        }

        void RK3_BSSN::unzipVars(double** zipIn , double** uzipOut)
        {
            bssn::timer::t_unzip_sync.start();

            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiMesh->unzip(zipIn[index],uzipOut[index]);

            bssn::timer::t_unzip_sync.stop();

        }

        void RK3_BSSN::unzipVars_async(double ** zipIn , double **uzipOut)
        {

            /*bssn::timer::t_unzip_async.start();
            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiMesh->unzip_async(zipIn[index],uzipOut[index]);
            bssn::timer::t_unzip_async.stop();*/
            bssn::timer::t_unzip_async.start();

            for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var+=bssn::BSSN_ASYNC_COMM_K){

                for(unsigned int i=0;(i<bssn::BSSN_ASYNC_COMM_K);i++)
                    m_uiMesh->ghostExchangeStart(zipIn[var+i],m_uiSendNodeBuf[i],m_uiRecvNodeBuf[i],m_uiSendReqs[i],m_uiRecvReqs[i]);

                for(unsigned int i=0;(i<bssn::BSSN_ASYNC_COMM_K);i++)
                {
                    m_uiMesh->ghostExchangeRecvSync(zipIn[var + i], m_uiRecvNodeBuf[i],m_uiRecvReqs[i], m_uiRecvSts[i]);
                    m_uiMesh->unzip(zipIn[var+i],uzipOut[var+i]);
                }

                for(unsigned int i=0;(i<bssn::BSSN_ASYNC_COMM_K);i++)
                    m_uiMesh->ghostExchangeSendSync(m_uiSendReqs[i], m_uiSendSts[i]);

            }

            bssn::timer::t_unzip_async.stop();


        }


        void RK3_BSSN::zipVars(double** uzipIn , double** zipOut)
        {
            bssn::timer::t_zip.start();

            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                m_uiMesh->zip(uzipIn[index],zipOut[index]);

            bssn::timer::t_zip.stop();

        }


        void RK3_BSSN::applyBoundaryConditions()
        {

        }

        void RK3_BSSN::performSingleIteration()
        {

            bool repeatStep;
            double n_inf_max=0.0;
            double n_inf_max_g=0;

            do
            {

                repeatStep=false;
                n_inf_max=0;
                n_inf_max_g=0;

                if(m_uiMesh->isActive())
                {
                    double current_t=m_uiCurrentTime;
                    double current_t_adv=current_t;

#ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                    unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
#else
                    //1. perform ghost exchange.
                    performGhostExchangeVars(m_uiPrevVar);
                    //2. unzip all the variables.
                    unzipVars(m_uiPrevVar,m_uiUnzipVar);
#endif


                    int rank =m_uiMesh->getMPIRank();



                    const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
                    const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();

                    const std::vector<ot::Block>& blkList=m_uiMesh->getLocalBlockList();

                    unsigned int offset;
                    double ptmin[3], ptmax[3];
                    unsigned int sz[3];
                    unsigned int bflag;
                    double dx,dy,dz;
                    const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
                    const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);


                    for(unsigned int stage=0;stage<(bssn::BSSN_RK3_STAGES-1);stage++)
                    {


#ifdef DEBUG_RK_SOLVER
                        if(!rank)std::cout<<" stage: "<<stage<<" begin: "<<std::endl;
                    for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);
#endif

#ifdef BSSN_ENABLE_CUDA
                        cuda::BSSNComputeParams bssnParams;
                        bssnParams.BSSN_LAMBDA[0]=bssn::BSSN_LAMBDA[0];
                        bssnParams.BSSN_LAMBDA[1]=bssn::BSSN_LAMBDA[1];
                        bssnParams.BSSN_LAMBDA[2]=bssn::BSSN_LAMBDA[2];
                        bssnParams.BSSN_LAMBDA[3]=bssn::BSSN_LAMBDA[3];

                        bssnParams.BSSN_LAMBDA_F[0]=bssn::BSSN_LAMBDA_F[0];
                        bssnParams.BSSN_LAMBDA_F[1]=bssn::BSSN_LAMBDA_F[1];

                        bssnParams.BSSN_ETA_POWER[0]=bssn::BSSN_ETA_POWER[0];
                        bssnParams.BSSN_ETA_POWER[1]=bssn::BSSN_ETA_POWER[1];

                        bssnParams.ETA_R0=bssn::ETA_R0;
                        bssnParams.ETA_CONST=bssn::ETA_CONST;
                        bssnParams.ETA_DAMPING=bssn::ETA_DAMPING;
                        bssnParams.ETA_DAMPING_EXP=bssn::ETA_DAMPING_EXP;
                        bssnParams.KO_DISS_SIGMA=bssn::KO_DISS_SIGMA;

                        dim3 threadBlock(16,16,1);
                        cuda::computeRHS(m_uiUnzipVarRHS,(const double **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size(),(const cuda::BSSNComputeParams*) &bssnParams,threadBlock,pt_min,pt_max,1);
#else

                        for(unsigned int blk=0;blk<blkList.size();blk++)
                        {
                            offset=blkList[blk].getOffset();
                            sz[0]=blkList[blk].getAllocationSzX();
                            sz[1]=blkList[blk].getAllocationSzY();
                            sz[2]=blkList[blk].getAllocationSzZ();

                            bflag=blkList[blk].getBlkNodeFlag();

                            dx=blkList[blk].computeDx(pt_min,pt_max);
                            dy=blkList[blk].computeDy(pt_min,pt_max);
                            dz=blkList[blk].computeDz(pt_min,pt_max);

                            ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
                            ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
                            ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

                            ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
                            ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
                            ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

#ifdef BSSN_RHS_STAGED_COMP

                            bssnrhs_sep(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
#else
                            bssnrhs(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
#endif


                        }
#endif

#ifdef DEBUG_RK_SOLVER
                        if(!rank)std::cout<<" stage: "<<stage<<" af rhs UNZIP RHS TEST:"<<std::endl;
                    for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        ot::test::isUnzipInternalNaN(m_uiMesh,m_uiUnzipVarRHS[index]);
#endif



                        zipVars(m_uiUnzipVarRHS,m_uiStage[stage]);


#ifdef DEBUG_RK_SOLVER
                        for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        if(seq::test::isNAN(m_uiStage[stage][index]+m_uiMesh->getNodeLocalBegin(),m_uiMesh->getNumLocalMeshNodes()))
                            std::cout<<" var: "<<index<<" contains nan af zip  stage: "<<stage<<std::endl;
#endif

                        /*for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                             for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                                     m_uiStage[stage][index][node]*=m_uiT_h;*/



                        for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                        {
                            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                            {
                                m_uiVarIm[index][node]=m_uiPrevVar[index][node];
                                for(unsigned int s=0;s<(stage+1);s++) {
                                    //if(!rank && index==0 && node==0) std::cout<<"rk stage: "<<stage<<" im coef: "<<s<<" value: "<<RK3_U[stage+1][s]<<std::endl;
                                    m_uiVarIm[index][node] += (RK3_U[stage + 1][s] * m_uiT_h * m_uiStage[s][index][node]);
                                }

                            }
                            enforce_bssn_constraints(m_uiVarIm, node);
                        }


                        current_t_adv=current_t+RK3_T[stage+1]*m_uiT_h;

#ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                        unzipVars_async(m_uiVarIm,m_uiUnzipVar);
#else
                        performGhostExchangeVars(m_uiVarIm);
                        unzipVars(m_uiVarIm,m_uiUnzipVar);
#endif

                    }

                    current_t_adv=current_t+RK3_T[(bssn::BSSN_RK3_STAGES-1)]*m_uiT_h;


#ifdef DEBUG_RK_SOLVER
                    if(!rank)std::cout<<" stage: "<<(bssn::BSSN_RK3_STAGES-1)<<" begin: "<<std::endl;

                for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                    ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);
#endif


#ifdef BSSN_ENABLE_CUDA
                        cuda::BSSNComputeParams bssnParams;
                        bssnParams.BSSN_LAMBDA[0]=bssn::BSSN_LAMBDA[0];
                        bssnParams.BSSN_LAMBDA[1]=bssn::BSSN_LAMBDA[1];
                        bssnParams.BSSN_LAMBDA[2]=bssn::BSSN_LAMBDA[2];
                        bssnParams.BSSN_LAMBDA[3]=bssn::BSSN_LAMBDA[3];

                        bssnParams.BSSN_LAMBDA_F[0]=bssn::BSSN_LAMBDA_F[0];
                        bssnParams.BSSN_LAMBDA_F[1]=bssn::BSSN_LAMBDA_F[1];

                        bssnParams.BSSN_ETA_POWER[0]=bssn::BSSN_ETA_POWER[0];
                        bssnParams.BSSN_ETA_POWER[1]=bssn::BSSN_ETA_POWER[1];

                        bssnParams.ETA_R0=bssn::ETA_R0;
                        bssnParams.ETA_CONST=bssn::ETA_CONST;
                        bssnParams.ETA_DAMPING=bssn::ETA_DAMPING;
                        bssnParams.ETA_DAMPING_EXP=bssn::ETA_DAMPING_EXP;
                        bssnParams.KO_DISS_SIGMA=bssn::KO_DISS_SIGMA;

                        dim3 threadBlock(16,16,1);
                        cuda::computeRHS(m_uiUnzipVarRHS,(const double **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size(),(const cuda::BSSNComputeParams*) &bssnParams,threadBlock,pt_min,pt_max,1);
#else

                    for(unsigned int blk=0;blk<blkList.size();blk++)
                    {

                        offset=blkList[blk].getOffset();
                        sz[0]=blkList[blk].getAllocationSzX();
                        sz[1]=blkList[blk].getAllocationSzY();
                        sz[2]=blkList[blk].getAllocationSzZ();

                        bflag=blkList[blk].getBlkNodeFlag();

                        dx=blkList[blk].computeDx(pt_min,pt_max);
                        dy=blkList[blk].computeDy(pt_min,pt_max);
                        dz=blkList[blk].computeDz(pt_min,pt_max);

                        ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
                        ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
                        ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

                        ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
                        ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
                        ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;


#ifdef BSSN_RHS_STAGED_COMP

                        bssnrhs_sep(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
#else
                        bssnrhs(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
#endif


                    }

#endif

#ifdef DEBUG_RK_SOLVER
                    if(!rank)std::cout<<" stage: "<<(bssn::BSSN_RK3_STAGES-1)<<" af rhs UNZIP RHS TEST:"<<std::endl;
                    for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        ot::test::isUnzipInternalNaN(m_uiMesh,m_uiUnzipVarRHS[index]);
#endif

                    zipVars(m_uiUnzipVarRHS,m_uiStage[(bssn::BSSN_RK3_STAGES-1)]);

                    /*for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                            m_uiStage[(bssn::BSSN_RK3_STAGES-1)][index][node]*=m_uiT_h;*/



                    for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                    {

                        for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        {
                            m_uiVar[index][node]=m_uiPrevVar[index][node];
                            for(unsigned int s=0;s<(bssn::BSSN_RK3_STAGES);s++)
                            {
                                m_uiVar[index][node]+=(RK3_C[s]*m_uiT_h*m_uiStage[s][index][node]);
                            }


                        }
                        enforce_bssn_constraints(m_uiVar, node);

                    }



                }

                m_uiMesh->waitAll();


            }while(repeatStep);


            m_uiCurrentStep++;
            m_uiCurrentTime+=m_uiT_h;


        }


        void RK3_BSSN::rkSolve()
        {

            if(m_uiCurrentStep==0)
            {
                applyInitialConditions(m_uiPrevVar);
            }

            bool isRefine=true;
            unsigned int oldElements,oldElements_g;
            unsigned int newElements,newElements_g;

            // refine based on all the variables
            const unsigned int refineNumVars=bssn::BSSN_NUM_REFINE_VARS;
            unsigned int refineVarIds[refineNumVars];
            for(unsigned int vIndex=0;vIndex<refineNumVars;vIndex++)
                refineVarIds[vIndex]=bssn::BSSN_REFINE_VARIABLE_INDICES[vIndex];

            double wTol=bssn::BSSN_WAVELET_TOL;
            std::function<double(double,double,double)> waveletTolFunc =[wTol](double x,double y, double z){ return bssn::computeWTol(x,y,z,wTol);};


            double l_min,l_max;
            for(double t=m_uiCurrentTime;t<m_uiTimeEnd;t=t+m_uiT_h)
            {

                // checkpoint the previous solution value before going to the next step.
                bssn::timer::t_ioCheckPoint.start();
                if((m_uiMesh->isActive()) && (m_uiCurrentStep%bssn::BSSN_CHECKPT_FREQ)==0)
                    storeCheckPoint(bssn::BSSN_CHKPT_FILE_PREFIX.c_str());
                bssn::timer::t_ioCheckPoint.stop();

                // write sol to vtu.
                if((m_uiMesh->isActive()) && (m_uiCurrentStep%bssn::BSSN_TIME_STEP_OUTPUT_FREQ)==0)
                {

                    l_min=vecMin(m_uiPrevVar[bssn::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                    l_max=vecMax(m_uiPrevVar[bssn::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                    if(!m_uiMesh->getMPIRank()) {
                        std::cout << "executing step: " << m_uiCurrentStep << " dt: " << m_uiT_h << " rk_time : "<< m_uiCurrentTime << std::endl;
                        std::cout << "\t ||VAR::U_ALPHA|| (min, max) : ("<<l_min<<", "<<l_max<<" ) "<<std::endl;
                    }

                    bssn::timer::profileInfoIntermediate(bssn::BSSN_PROFILE_FILE_PREFIX.c_str(),m_uiMesh,m_uiCurrentStep);

                }

                if((m_uiCurrentStep%bssn::BSSN_TIME_STEP_OUTPUT_FREQ)==0)
                    bssn::timer::resetSnapshot();


                if((m_uiCurrentStep%bssn::BSSN_REMESH_TEST_FREQ)==0)
                {

#ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                    unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
#else
                    performGhostExchangeVars(m_uiPrevVar);
                    //isRefine=m_uiMesh->isReMesh((const double **)m_uiPrevVar,refineVarIds,refineNumVars,bssn::BSSN_WAVELET_TOL);
                    unzipVars(m_uiPrevVar,m_uiUnzipVar);
#endif
                    
//                     char fPrefix[256];
//                     sprintf(fPrefix,"%s_wavelets_%d",bssn::BSSN_VTU_FILE_PREFIX.c_str(),m_uiCurrentStep);
//                     io::vtk::waveletsToVTU(m_uiMesh,(const double **)m_uiPrevVar,(const double**)m_uiUnzipVar,refineVarIds,refineNumVars,(const char*)fPrefix);

#ifdef DEBUG_RK_SOLVER
                    if(m_uiMesh->isActive())
                    {
                        if(!m_uiMesh->getMPIRank())std::cout<<" isRemesh Unzip : "<<std::endl;
                        for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                            ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);

                    }

#endif
                    bssn::timer::t_isReMesh.start();
                    if(bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
                        isRefine=false;
                    else
                        isRefine=m_uiMesh->isReMeshUnzip((const double **)m_uiUnzipVar,refineVarIds,refineNumVars,waveletTolFunc,bssn::BSSN_DENDRO_AMR_FAC);
                    bssn::timer::t_isReMesh.stop();

                    if(isRefine)
                    {


#ifdef DEBUG_IS_REMESH
			unsigned int rank=m_uiMesh->getMPIRankGlobal();
			MPI_Comm globalComm=m_uiMesh->getMPIGlobalCommunicator();
                        std::vector<ot::TreeNode> unChanged;
                        std::vector<ot::TreeNode> refined;
                        std::vector<ot::TreeNode> coarsened;
                        std::vector<ot::TreeNode> localBlocks;

                        const ot::Block* blkList=&(*(m_uiMesh->getLocalBlockList().begin()));
                        for(unsigned int ele=0;ele<m_uiMesh->getLocalBlockList().size();ele++)
                        {
                            localBlocks.push_back(blkList[ele].getBlockNode());
                        }


                        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
                        for(unsigned int ele=m_uiMesh->getElementLocalBegin();ele<m_uiMesh->getElementLocalEnd();ele++)
                        {
                            if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE)
                            {
                                unChanged.push_back(pNodes[ele]);
                            }else if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                            {
                                refined.push_back(pNodes[ele]);
                            }else
                            {
                                assert((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE);
                                coarsened.push_back(pNodes[ele]);
                            }
                        }

                        char fN1[256];
                        char fN2[256];
                        char fN3[256];
                        char fN4[256];

                        sprintf(fN1,"unchanged_%d",m_uiCurrentStep);
                        sprintf(fN2,"refined_%d",m_uiCurrentStep);
                        sprintf(fN3,"coarsend_%d",m_uiCurrentStep);
                        sprintf(fN4,"blocks_%d",m_uiCurrentStep);

                        DendroIntL localSz=unChanged.size();
                        DendroIntL globalSz;
                        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,globalComm);
                        if(!rank) std::cout<<" total unchanged: "<<globalSz<<std::endl;

                        localSz=refined.size();
                        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,globalComm);
                        if(!rank) std::cout<<" total refined: "<<globalSz<<std::endl;


                        localSz=coarsened.size();
                        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,globalComm);
                        if(!rank) std::cout<<" total coarsend: "<<globalSz<<std::endl;


                        io::vtk::oct2vtu(&(*(unChanged.begin())),unChanged.size(),fN1,globalComm);
                        io::vtk::oct2vtu(&(*(refined.begin())),refined.size(),fN2,globalComm);
                        io::vtk::oct2vtu(&(*(coarsened.begin())),coarsened.size(),fN3,globalComm);
                        io::vtk::oct2vtu(&(*(localBlocks.begin())),localBlocks.size(),fN4,globalComm);

#endif
                        bssn::timer::t_mesh.start();
                        ot::Mesh* newMesh=m_uiMesh->ReMesh(bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
                        bssn::timer::t_mesh.stop();

                        oldElements=m_uiMesh->getNumLocalMeshElements();
                        newElements=newMesh->getNumLocalMeshElements();

                        par::Mpi_Reduce(&oldElements,&oldElements_g,1,MPI_SUM,0,m_uiMesh->getMPIGlobalCommunicator());
                        par::Mpi_Reduce(&newElements,&newElements_g,1,MPI_SUM,0,newMesh->getMPIGlobalCommunicator());

                        if(!(m_uiMesh->getMPIRankGlobal()))std::cout<<"step : "<<m_uiCurrentStep<<" time : "<<m_uiCurrentTime<<" old mesh: "<<oldElements_g<<" new mesh: "<<newElements_g<<std::endl;


                        // performs the inter-grid transfer
                        intergridTransferVars(m_uiPrevVar,newMesh);

                        for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                        {
                            delete [] m_uiVar[index];
                            delete [] m_uiVarIm[index];
                            delete [] m_uiUnzipVar[index];
                            delete [] m_uiUnzipVarRHS[index];

                            m_uiVar[index]=NULL;
                            m_uiVarIm[index]=NULL;
                            m_uiUnzipVar[index]=NULL;
                            m_uiUnzipVarRHS[index]=NULL;

                            m_uiVar[index]=newMesh->createVector<double>();
                            m_uiVarIm[index]=newMesh->createVector<double>();
                            m_uiUnzipVar[index]=newMesh->createUnZippedVector<double>();
                            m_uiUnzipVarRHS[index]=newMesh->createUnZippedVector<double>();


                        }

                        for(unsigned int stage=0;stage<bssn::BSSN_RK3_STAGES;stage++)
                            for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                            {
                                delete [] m_uiStage[stage][index];
                                m_uiStage[stage][index]=NULL;
                                m_uiStage[stage][index]=newMesh->createVector<double>();
                            }

                        // deallocate constraint vars allocate them for the new mesh.
                        for(unsigned int index=0;index<bssn::BSSN_CONSTRAINT_NUM_VARS;index++)
                        {
                            delete [] m_uiConstraintVars[index];
                            delete [] m_uiUnzipConstraintVars[index];

                            m_uiConstraintVars[index]=newMesh->createVector<double>();
                            m_uiUnzipConstraintVars[index]=newMesh->createUnZippedVector<double>();

                        }




                        std::swap(newMesh,m_uiMesh);
                        delete newMesh;

#ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                        // reallocates mpi resources for the the new mesh. (this will deallocate the old resources)
                        reallocateMPIResources();
#endif

                        if(m_uiMesh->isActive())
                        {
                            l_min=vecMin(m_uiPrevVar[bssn::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                            l_max=vecMax(m_uiPrevVar[bssn::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                            if(!(m_uiMesh->getMPIRank())) {
                                std::cout << "transfer completed:    ||VAR::U_ALPHA|| (min, max) : ("<<l_min<<", "<<l_max<<" ) "<<std::endl;
                            }

                        }



                    }

                }




                if((m_uiCurrentStep%bssn::BSSN_IO_OUTPUT_FREQ)==0)
                {

#ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                    unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
#else
                    performGhostExchangeVars(m_uiPrevVar);
                    unzipVars(m_uiPrevVar,m_uiUnzipVar);
#endif


#ifdef BSSN_COMPUTE_CONSTRAINTS

                    const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();

                    unsigned int offset;
                    double ptmin[3], ptmax[3];
                    unsigned int sz[3];
                    unsigned int bflag;
                    double dx,dy,dz;
                    const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
                    const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);

                    for(unsigned int blk=0;blk<blkList.size();blk++)
                    {
                        offset=blkList[blk].getOffset();
                        sz[0]=blkList[blk].getAllocationSzX();
                        sz[1]=blkList[blk].getAllocationSzY();
                        sz[2]=blkList[blk].getAllocationSzZ();

                        bflag=blkList[blk].getBlkNodeFlag();

                        dx=blkList[blk].computeDx(pt_min,pt_max);
                        dy=blkList[blk].computeDy(pt_min,pt_max);
                        dz=blkList[blk].computeDz(pt_min,pt_max);

                        ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
                        ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
                        ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

                        ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
                        ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
                        ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

                        physical_constraints(m_uiUnzipConstraintVars, (const double **) m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
                    }

                    /*double consVecMin[bssn::BSSN_CONSTRAINT_NUM_VARS];
                    double consVecMax[bssn::BSSN_CONSTRAINT_NUM_VARS];*/
                    double constraintMaskedL2[bssn::BSSN_CONSTRAINT_NUM_VARS];
                    for(unsigned int index=0;index<bssn::BSSN_CONSTRAINT_NUM_VARS;index++)
                    {
                        m_uiMesh->zip(m_uiUnzipConstraintVars[index],m_uiConstraintVars[index]);
                        m_uiMesh->performGhostExchange(m_uiConstraintVars[index]);
                    }

                    bssn::extractConstraints(m_uiMesh,(const double **)m_uiConstraintVars,m_uiCurrentStep);
#ifndef BSSN_KERR_SCHILD_TEST
                    bssn::extractGravitationlWaves(m_uiMesh,(const double **)m_uiConstraintVars,50,m_uiCurrentStep);
#endif

#endif

#ifdef BSSN_ENABLE_VTU_OUTPUT
                    writeToVTU(m_uiPrevVar,m_uiConstraintVars,bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT,bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT,bssn::BSSN_VTU_OUTPUT_EVOL_INDICES,bssn::BSSN_VTU_OUTPUT_CONST_INDICES);
#endif

                }



                bssn::timer::t_rkStep.start();

                performSingleIteration();

                bssn::timer::t_rkStep.stop();



                std::swap(m_uiVar,m_uiPrevVar);

                //if(m_uiCurrentStep==1) break;

            }


        }


        void RK3_BSSN::storeCheckPoint(const char * fNamePrefix)
        {


            if(m_uiMesh->isActive())
            {
                unsigned int cpIndex;
                (m_uiCurrentStep%(2*bssn::BSSN_CHECKPT_FREQ)==0) ? cpIndex=0 : cpIndex=1; // to support alternate file writing.
                unsigned int rank=m_uiMesh->getMPIRank();
                unsigned int npes=m_uiMesh->getMPICommSize();

                char fName[256];
                const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
                sprintf(fName,"%s_octree_%d_%d.oct",fNamePrefix,cpIndex,rank);
                io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

                unsigned int numVars=bssn::BSSN_NUM_VARS;
                const char ** varNames=bssn::BSSN_VAR_NAMES;

                /*for(unsigned int i=0;i<numVars;i++)
                {
                    sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,rank);
                    io::checkpoint::writeVecToFile(fName,m_uiMesh,m_uiPrevVar[i]);
                }*/

                sprintf(fName,"%s_%d_%d.var",fNamePrefix,cpIndex,rank);
                io::checkpoint::writeVecToFile(fName,m_uiMesh,(const double **)m_uiPrevVar,bssn::BSSN_NUM_VARS);

                if(!rank)
                {
                    sprintf(fName,"%s_step_%d.cp",fNamePrefix,cpIndex);
                    std::cout<<"writing : "<<fName<<std::endl;
                    std::ofstream outfile(fName);
                    if(!outfile) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                    json checkPoint;
                    checkPoint["DENDRO_RK45_TIME_BEGIN"]=m_uiTimeBegin;
                    checkPoint["DENDRO_RK45_TIME_END"]=m_uiTimeEnd;
                    checkPoint["DENDRO_RK45_ELEMENT_ORDER"]=m_uiOrder;

                    checkPoint["DENDRO_RK45_TIME_CURRENT"]=m_uiCurrentTime;
                    checkPoint["DENDRO_RK45_STEP_CURRENT"]=m_uiCurrentStep;
                    checkPoint["DENDRO_RK45_TIME_STEP_SIZE"]=m_uiT_h;
                    checkPoint["DENDRO_RK45_LAST_IO_TIME"]=m_uiCurrentTime;

                    checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"]=bssn::BSSN_WAVELET_TOL;
                    checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"]=bssn::BSSN_LOAD_IMB_TOL;
                    checkPoint["DENDRO_RK45_NUM_VARS"]=numVars; // number of variables to restore.
                    checkPoint["DENDRO_RK45_ACTIVE_COMM_SZ"]=m_uiMesh->getMPICommSize(); // (note that rank 0 is always active).

                    outfile<<std::setw(4)<<checkPoint<<std::endl;
                    outfile.close();

                }


            }



        }


        void RK3_BSSN::restoreCheckPoint(const char * fNamePrefix,MPI_Comm comm)
        {
            unsigned int numVars=0;
            std::vector<ot::TreeNode> octree;
            json checkPoint;

            int rank;
            int npes;
            m_uiComm=comm;
            MPI_Comm_rank(m_uiComm,&rank);
            MPI_Comm_size(m_uiComm,&npes);

            unsigned int activeCommSz;

            char fName[256];
            unsigned int restoreStatus=0;
            unsigned int restoreStatusGlobal=0; // 0 indicates successfully restorable.

            ot::Mesh* newMesh;
            unsigned int restoreStep[2];
            restoreStep[0]=0;
            restoreStep[1]=0;

            unsigned int restoreFileIndex=0;

            for(unsigned int cpIndex=0;cpIndex<2;cpIndex++) {

                restoreStatus=0;

                if(!rank)
                {
                    sprintf(fName,"%s_step_%d.cp",fNamePrefix,cpIndex);
                    std::ifstream infile(fName);
                    if(!infile) {std::cout<<fName<<" file open failed "<<std::endl; restoreStatus=1;}


                    if(restoreStatus==0)
                    {
                        infile>>checkPoint;
                        m_uiTimeBegin=checkPoint["DENDRO_RK45_TIME_BEGIN"];
                        m_uiTimeEnd=checkPoint["DENDRO_RK45_TIME_END"];

                        m_uiOrder=checkPoint["DENDRO_RK45_ELEMENT_ORDER"];
                        m_uiCurrentTime=checkPoint["DENDRO_RK45_TIME_CURRENT"];
                        m_uiCurrentStep=checkPoint["DENDRO_RK45_STEP_CURRENT"];
                        m_uiT_h=checkPoint["DENDRO_RK45_TIME_STEP_SIZE"];


                        bssn::BSSN_WAVELET_TOL=checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"];
                        bssn::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"];
                        numVars=checkPoint["DENDRO_RK45_NUM_VARS"];
                        activeCommSz=checkPoint["DENDRO_RK45_ACTIVE_COMM_SZ"];

                        restoreStep[cpIndex]=m_uiCurrentStep;

                    }
                }

            }

            if(!rank)
            {
                if(restoreStep[0]<restoreStep[1])
                    restoreFileIndex=1;
                else
                    restoreFileIndex=0;
            }

            par::Mpi_Bcast(&restoreFileIndex,1,0,m_uiComm);

            for(unsigned int cpIndex=restoreFileIndex;cpIndex<2;cpIndex++)
            {
                restoreStatus=0;
                octree.clear();
                if(!rank) std::cout<<" Trying to restore from checkpoint index : "<<cpIndex<<std::endl;

                if(!rank)
                {
                    sprintf(fName,"%s_step_%d.cp",fNamePrefix,cpIndex);
                    std::ifstream infile(fName);
                    if(!infile) {std::cout<<fName<<" file open failed "<<std::endl; restoreStatus=1;}


                    if(restoreStatus==0)
                    {
                        infile>>checkPoint;
                        m_uiTimeBegin=checkPoint["DENDRO_RK45_TIME_BEGIN"];
                        m_uiTimeEnd=checkPoint["DENDRO_RK45_TIME_END"];

                        m_uiOrder=checkPoint["DENDRO_RK45_ELEMENT_ORDER"];
                        m_uiCurrentTime=checkPoint["DENDRO_RK45_TIME_CURRENT"];
                        m_uiCurrentStep=checkPoint["DENDRO_RK45_STEP_CURRENT"];
                        m_uiT_h=checkPoint["DENDRO_RK45_TIME_STEP_SIZE"];


                        bssn::BSSN_WAVELET_TOL=checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"];
                        bssn::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"];
                        numVars=checkPoint["DENDRO_RK45_NUM_VARS"];
                        activeCommSz=checkPoint["DENDRO_RK45_ACTIVE_COMM_SZ"];
                    }


                }

                par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,m_uiComm);
                if(restoreStatusGlobal==1) continue;


                par::Mpi_Bcast(&m_uiTimeBegin,1,0,comm);
                par::Mpi_Bcast(&m_uiTimeEnd,1,0,comm);

                par::Mpi_Bcast(&m_uiCurrentTime,1,0,comm);
                par::Mpi_Bcast(&m_uiCurrentStep,1,0,comm);

                par::Mpi_Bcast(&m_uiT_h,1,0,comm);


                par::Mpi_Bcast(& bssn::BSSN_WAVELET_TOL,1,0,comm);
                par::Mpi_Bcast(&bssn::BSSN_LOAD_IMB_TOL,1,0,comm);

                par::Mpi_Bcast(&numVars,1,0,comm);
                par::Mpi_Bcast(&m_uiOrder,1,0,comm);
                par::Mpi_Bcast(&m_uiT_h,1,0,comm);

                par::Mpi_Bcast(&activeCommSz,1,0,comm);

                if(activeCommSz>npes)
                {
                    if(!rank)std::cout<<" [Error] : checkpoint file written from  a larger communicator than the current global comm. (i.e. communicator shrinking not allowed in the restore step. )"<<std::endl;
                    exit(0);
                }



                bool isActive=(rank<activeCommSz);

                MPI_Comm newComm;
                par::splitComm2way(isActive,&newComm,m_uiComm);

                if(isActive) {

                    int activeRank;
                    int activeNpes;

                    MPI_Comm_rank(newComm, &activeRank);
                    MPI_Comm_size(newComm, &activeNpes);
                    assert(activeNpes == activeCommSz);

                    sprintf(fName, "%s_octree_%d_%d.oct", fNamePrefix,cpIndex,activeRank);
                    restoreStatus=io::checkpoint::readOctFromFile(fName, octree);
                    assert(par::test::isUniqueAndSorted(octree, newComm));

                }

                par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,m_uiComm);
                if(restoreStatusGlobal==1) {

                    if(!rank) std::cout<<"[Error]: octree (*.oct) restore file currupted "<<std::endl;
                    continue;
                }

                newMesh=new ot::Mesh(octree,1,m_uiOrder,activeCommSz,m_uiComm);

                for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                {
                    delete [] m_uiPrevVar[index];
                    delete [] m_uiVar[index];
                    delete [] m_uiVarIm[index];
                    delete [] m_uiUnzipVar[index];
                    delete [] m_uiUnzipVarRHS[index];

                    m_uiPrevVar[index]=newMesh->createVector<double>();
                    m_uiVar[index]=newMesh->createVector<double>();
                    m_uiVarIm[index]=newMesh->createVector<double>();
                    m_uiUnzipVar[index]=newMesh->createUnZippedVector<double>();
                    m_uiUnzipVarRHS[index]=newMesh->createUnZippedVector<double>();

                }

                for(unsigned int stage=0;stage<bssn::BSSN_RK3_STAGES;stage++)
                    for(unsigned int index=0;index<bssn::BSSN_NUM_VARS;index++)
                    {
                        delete [] m_uiStage[stage][index];
                        m_uiStage[stage][index]=newMesh->createVector<double>();
                    }

                // deallocate constraint vars allocate them for the new mesh.
                for(unsigned int index=0;index<bssn::BSSN_CONSTRAINT_NUM_VARS;index++)
                {
                    delete [] m_uiConstraintVars[index];
                    delete [] m_uiUnzipConstraintVars[index];

                    m_uiConstraintVars[index]=newMesh->createVector<double>();
                    m_uiUnzipConstraintVars[index]=newMesh->createUnZippedVector<double>();

                }

                const char ** varNames=bssn::BSSN_VAR_NAMES;


                if(isActive) {

                    int activeRank;
                    int activeNpes;

                    MPI_Comm_rank(newComm, &activeRank);
                    MPI_Comm_size(newComm, &activeNpes);
                    assert(activeNpes == activeCommSz);

                    /*for(unsigned int i=0;i<numVars;i++)
                    {
                        sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,activeRank);
                        restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,m_uiPrevVar[i]);
                        if(restoreStatus==1) break;

                    }*/

                    sprintf(fName,"%s_%d_%d.var",fNamePrefix,cpIndex,activeRank);
                    restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,m_uiPrevVar,bssn::BSSN_NUM_VARS);


                }
                MPI_Comm_free(&newComm);
                par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,m_uiComm);
                if(restoreStatusGlobal==1) {

                    if(!rank) std::cout<<"[Error]: varible (*.var) restore file currupted "<<std::endl;
                    continue;
                }

                std::swap(m_uiMesh,newMesh);
                delete newMesh;
                reallocateMPIResources();
                if(restoreStatusGlobal==0) break;

            }

            if(restoreStatusGlobal==1){
                std::cout<<"rank: "<<rank<<"[Error]: rk solver restore error "<<std::endl;
                exit(0);
            }

            unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
            unsigned int totalElems;

            par::Mpi_Allreduce(&localSz,&totalElems,1,MPI_SUM,m_uiComm);


            if(!rank) std::cout<<" checkpoint at step : "<<m_uiCurrentStep<<"active Comm. sz: "<<activeCommSz<<" restore successful: "<<" restored mesh size: "<<totalElems<<std::endl;
            return;


         }



    }// end of namespace solver
} // end of namespace ode
