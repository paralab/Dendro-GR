//
// Created by milinda on 5/24/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains specific example of solving 3D transport equation using RK45 method
*/
//

#include "rkTransport.h"

namespace ode
{
    namespace solver
    {

        RK45Transport::RK45Transport(ot::Mesh *pMesh, double pTBegin, double pTEnd,double pTh): RK(pMesh,pTBegin,pTEnd,pTh)
        {
            m_uiU.clear();
            m_uiPrevU.clear();
            m_uiLoadImbTol=0.1;
            m_uiSplitterFix=2;
            m_uiLastIOTimeValue=pTBegin;

            if(m_uiOrder==4)
            {
                m_uiDxOrder4_centered=Stencil<double,5,2>(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_X);
                m_uiDxOrder4_backward=Stencil<double,5,4>(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_X);
                m_uiDxOrder4_forward=Stencil<double,5,0>(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_X);

                m_uiDyOrder4_centered=Stencil<double,5,2>(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_Y);
                m_uiDyOrder4_backward=Stencil<double,5,4>(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_Y);
                m_uiDyOrder4_forward=Stencil<double,5,0>(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_Y);

                m_uiDzOrder4_centered=Stencil<double,5,2>(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_Z);
                m_uiDzOrder4_backward=Stencil<double,5,4>(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_Z);
                m_uiDzOrder4_forward=Stencil<double,5,0>(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_Z);
            }

        }

        RK45Transport::~RK45Transport()
        {

        }

        void RK45Transport::setParameters(double *b, std::function<double(double, double, double)> g,std::function<double (double,double,double,double)>f,const char * fprefix,double wltTol )
        {

            m_uiWaveletTol=wltTol;
            m_uiCoef_b[0]=b[0];
            m_uiCoef_b[1]=b[1];
            m_uiCoef_b[2]=b[2];

            //std::cout<<"vec b: "<<m_uiCoef_b[0]<<" ,"<<m_uiCoef_b[1]<<" , "<<m_uiCoef_b[2]<<std::endl;

            m_uiFunc_g=g;
            m_uiFunc_f=f;


            if(fprefix==NULL)
            {
                m_uiFilePrefix=(char*)"sol_step";
            }else
                m_uiFilePrefix=(char *)fprefix;

            m_uiMesh->createVector(m_uiU);
            m_uiMesh->createVector(m_uiPrevU);
            m_uiMesh->createVector(m_uiDxU);
            m_uiMesh->createVector(m_uiDyU);
            m_uiMesh->createVector(m_uiDzU);

            m_uiMesh->createVector(m_uiStage_1);
            m_uiMesh->createVector(m_uiStage_2);
            m_uiMesh->createVector(m_uiStage_3);
            m_uiMesh->createVector(m_uiStage_4);
            m_uiMesh->createVector(m_uiStage_5);
            m_uiMesh->createVector(m_uiStage_6);

            m_uiMesh->createVector(m_uiDxU); // partial x derivative
            m_uiMesh->createVector(m_uiDyU); // partial y derivative
            m_uiMesh->createVector(m_uiDzU); // partial z derivative

            m_uiMesh->createVector(m_uiVecIm);

            m_uiMesh->createUnZippedVector(m_uiUZipVecIn,0.0);
            m_uiMesh->createUnZippedVector(m_uiUZipVecOut,0.0);


        }

        void  RK45Transport::setLDImbParameters(double ldTol,unsigned int sfK)
        {
            m_uiLoadImbTol=ldTol;
            m_uiSplitterFix=sfK;
        }

        void RK45Transport::applyInitialConditions()
        {
            m_uiMesh->createVector(m_uiPrevU,m_uiFunc_g);
            m_uiMesh->createVector(m_uiU,m_uiFunc_g);

        }

        void RK45Transport::performSingleIteration()
        {

            double rk45_linf=0;
            double rk45_inf_g=0;
            bool repeatStep;
            do
            {

                repeatStep=false;
                rk45_linf=0;
                rk45_inf_g=0;

                if(m_uiMesh->isActive())
                {

                    double current_t=m_uiCurrentTime;
                    double current_t_adv=current_t;


                    m_uiMesh->performGhostExchange(m_uiPrevU);
                    m_uiMesh->unzip(&(*(m_uiPrevU.begin())),&(*(m_uiUZipVecIn.begin())));

                    grad(m_uiMesh,0,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDxU.begin())));

                    grad(m_uiMesh,1,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDyU.begin())));

                    grad(m_uiMesh,2,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDzU.begin())));


                    m_uiFunc_f1=[current_t](const double x,const double y,const double z){return 0;};
                    m_uiMesh->createVector(m_uiRHS,m_uiFunc_f1);

                    unsigned int nLocalBegin=m_uiMesh->getNodeLocalBegin();
                    unsigned int nLocalEnd=m_uiMesh->getNodeLocalEnd();

                    for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
                    {
                        m_uiStage_1[k]=m_uiT_h*(m_uiRHS[k]-m_uiCoef_b[0]*m_uiDxU[k]-m_uiCoef_b[1]*m_uiDyU[k]-m_uiCoef_b[2]*m_uiDzU[k]);
                        m_uiVecIm[k]=m_uiPrevU[k]+RK45_STAGE_2_U*m_uiStage_1[k];
                    }


                    current_t_adv=current_t+RK45_STAGE_2_T*m_uiT_h;
                    //m_uiFunc_f1=[current_t_adv](const double x,const double y,const double z){return m_uiFunc_f(x,y,z,current_t_adv);};
                    m_uiMesh->createVector(m_uiRHS,m_uiFunc_f1);

                    applyBoundaryConditions(current_t_adv,&(*(m_uiVecIm.begin())));
                    m_uiMesh->performGhostExchange(m_uiVecIm);

                    m_uiMesh->unzip(&(*(m_uiVecIm.begin())),&(*(m_uiUZipVecIn.begin())));

                    grad(m_uiMesh,0,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDxU.begin())));

                    grad(m_uiMesh,1,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDyU.begin())));

                    grad(m_uiMesh,2,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDzU.begin())));

                    for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
                    {
                        m_uiStage_2[k]=m_uiT_h*(m_uiRHS[k]-m_uiCoef_b[0]*m_uiDxU[k]-m_uiCoef_b[1]*m_uiDyU[k]-m_uiCoef_b[2]*m_uiDzU[k]);
                        m_uiVecIm[k]=m_uiPrevU[k]+ (RK45_STAGE_3_U1)*m_uiStage_1[k]+(RK45_STAGE_3_U2)*m_uiStage_2[k];
                    }



                    current_t_adv=current_t+RK45_STAGE_3_T*m_uiT_h;
                    //m_uiFunc_f1=[current_t_adv](const double x,const double y,const double z){return m_uiFunc_f(x,y,z,current_t_adv);};
                    m_uiMesh->createVector(m_uiRHS,m_uiFunc_f1);
                    applyBoundaryConditions(current_t_adv,&(*(m_uiVecIm.begin())));
                    m_uiMesh->performGhostExchange(m_uiVecIm);

                    m_uiMesh->unzip(&(*(m_uiVecIm.begin())),&(*(m_uiUZipVecIn.begin())));

                    grad(m_uiMesh,0,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDxU.begin())));

                    grad(m_uiMesh,1,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDyU.begin())));

                    grad(m_uiMesh,2,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDzU.begin())));

                    for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
                    {
                        m_uiStage_3[k]=m_uiT_h*(m_uiRHS[k]-m_uiCoef_b[0]*m_uiDxU[k]-m_uiCoef_b[1]*m_uiDyU[k]-m_uiCoef_b[2]*m_uiDzU[k]);
                        m_uiVecIm[k]=m_uiPrevU[k]+ (RK45_STAGE_4_U1)*m_uiStage_1[k]+(RK45_STAGE_4_U2)*m_uiStage_2[k]+(RK45_STAGE_4_U3)*m_uiStage_3[k];
                    }



                    current_t_adv=current_t+RK45_STAGE_4_T*m_uiT_h;
                    //m_uiFunc_f1=[current_t_adv](const double x,const double y,const double z){return m_uiFunc_f(x,y,z,current_t_adv);};
                    m_uiMesh->createVector(m_uiRHS,m_uiFunc_f1);
                    applyBoundaryConditions(current_t_adv,&(*(m_uiVecIm.begin())));
                    m_uiMesh->performGhostExchange(m_uiVecIm);

                    m_uiMesh->unzip(&(*(m_uiVecIm.begin())),&(*(m_uiUZipVecIn.begin())));

                    grad(m_uiMesh,0,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDxU.begin())));

                    grad(m_uiMesh,1,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDyU.begin())));

                    grad(m_uiMesh,2,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDzU.begin())));

                    for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
                    {
                        m_uiStage_4[k]=m_uiT_h*(m_uiRHS[k]-m_uiCoef_b[0]*m_uiDxU[k]-m_uiCoef_b[1]*m_uiDyU[k]-m_uiCoef_b[2]*m_uiDzU[k]);
                        m_uiVecIm[k]=m_uiPrevU[k]+ (RK45_STAGE_5_U1)*m_uiStage_1[k]+(RK45_STAGE_5_U2)*m_uiStage_2[k]+(RK45_STAGE_5_U3)*m_uiStage_3[k]+(RK45_STAGE_5_U4)*m_uiStage_4[k];
                    }



                    current_t_adv=current_t+RK45_STAGE_5_T*m_uiT_h;
                    //m_uiFunc_f1=[current_t_adv](const double x,const double y,const double z){return m_uiFunc_f(x,y,z,current_t_adv);};
                    m_uiMesh->createVector(m_uiRHS,m_uiFunc_f1);
                    applyBoundaryConditions(current_t_adv,&(*(m_uiVecIm.begin())));
                    m_uiMesh->performGhostExchange(m_uiVecIm);

                    m_uiMesh->unzip(&(*(m_uiVecIm.begin())),&(*(m_uiUZipVecIn.begin())));

                    grad(m_uiMesh,0,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDxU.begin())));

                    grad(m_uiMesh,1,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDyU.begin())));

                    grad(m_uiMesh,2,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDzU.begin())));

                    for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
                    {
                        m_uiStage_5[k]=m_uiT_h*(m_uiRHS[k]-m_uiCoef_b[0]*m_uiDxU[k]-m_uiCoef_b[1]*m_uiDyU[k]-m_uiCoef_b[2]*m_uiDzU[k]);
                        m_uiVecIm[k]=m_uiPrevU[k]+ (RK45_STAGE_6_U1)*m_uiStage_1[k]+(RK45_STAGE_6_U2)*m_uiStage_2[k]+(RK45_STAGE_6_U3)*m_uiStage_3[k]+(RK45_STAGE_6_U4)*m_uiStage_4[k]+(RK45_STAGE_6_U5)*m_uiStage_5[k];
                    }



                    current_t_adv=current_t+RK45_STAGE_6_T*m_uiT_h;
                    //m_uiFunc_f1=[current_t_adv](const double x,const double y,const double z){return m_uiFunc_f(x,y,z,current_t_adv);};
                    m_uiMesh->createVector(m_uiRHS,m_uiFunc_f1);
                    applyBoundaryConditions(current_t_adv,&(*(m_uiVecIm.begin())));
                    m_uiMesh->performGhostExchange(m_uiVecIm);


                    m_uiMesh->unzip(&(*(m_uiVecIm.begin())),&(*(m_uiUZipVecIn.begin())));

                    grad(m_uiMesh,0,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDxU.begin())));

                    grad(m_uiMesh,1,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDyU.begin())));

                    grad(m_uiMesh,2,&(*(m_uiUZipVecIn.begin())),&(*(m_uiUZipVecOut.begin())));
                    m_uiMesh->zip(&(*(m_uiUZipVecOut.begin())),&(*(m_uiDzU.begin())));

                    for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
                    {
                        m_uiStage_6[k]=m_uiT_h*(m_uiRHS[k]-m_uiCoef_b[0]*m_uiDxU[k]-m_uiCoef_b[1]*m_uiDyU[k]-m_uiCoef_b[2]*m_uiDzU[k]);
                        m_uiU[k]=m_uiPrevU[k]+RK45_STAGE_1_COEF*m_uiStage_1[k]+RK45_STAGE_3_COEF*m_uiStage_3[k]+RK45_STAGE_4_COEF*m_uiStage_4[k]+RK45_STAGE_5_COEF*m_uiStage_5[k]+RK45_STAGE_6_COEF*m_uiStage_6[k];
                        m_uiVecIm[k]=m_uiPrevU[k]+(25.0/216.0)*m_uiStage_1[k]+(1408.0/2565.0)*m_uiStage_3[k]+(2197.0/4101.0)*m_uiStage_4[k]+(-1.0/5.0)*m_uiStage_5[k];

                    }


                    applyBoundaryConditions(m_uiCurrentTime+m_uiT_h,&(*(m_uiU.begin())));
                    applyBoundaryConditions(m_uiCurrentTime+m_uiT_h,&(*(m_uiVecIm.begin())));


                    rk45_linf=normLInfty(&(*(m_uiU.begin()+nLocalBegin)),(&*(m_uiVecIm.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),m_uiMesh->getMPICommunicator());


                }

                par::Mpi_Allreduce(&rk45_linf,&rk45_inf_g,1,MPI_MAX,m_uiComm);
                rk45_linf=rk45_inf_g;

                const double rk_tol=1e-3;

                if(rk45_linf>rk_tol)
                {
                    repeatStep=true;
                    m_uiT_h=0.8*m_uiT_h*(pow(fabs(rk_tol/rk45_linf),0.25));
                    if(!m_uiMesh->getMPIRank()) std::cout<<" repeat : "<<m_uiCurrentStep<<" with : "<<m_uiT_h<<std::endl;

                }else
                {
                    repeatStep=false;
                    m_uiT_h=0.8*m_uiT_h*(pow(fabs(rk_tol/rk45_linf),0.20));
                }


            }while(repeatStep);

            m_uiCurrentStep++;
            m_uiCurrentTime+=m_uiT_h;




        }



        void RK45Transport::applyBoundaryConditions(){}



        void RK45Transport::rkSolve()
        {
            if(m_uiCurrentStep==0)applyInitialConditions();

            bool isRefine=false;
            const ot::TreeNode * pNodes;
            unsigned int numRefinedElements=0;
            unsigned int numCoarsenElements=0;
            unsigned int numNotChangedElements=0;
            const unsigned int refineNumVars=1;
            unsigned int refineVarIds[refineNumVars];
            refineVarIds[0]=0;

            const double* varPtr []={&(*(m_uiPrevU.begin()))};
            double wTol=m_uiWaveletTol;
            std::function<double(double,double,double,double* hx)> waveletTolFunc =[wTol](double x,double y, double z,double *hx){ return wTol;};

            for(double t=m_uiCurrentTime;t<m_uiTimeEnd;t=t+m_uiT_h)
            {


                if((m_uiCurrentStep%IO_OUTPUT_FREQ)==0)
                {

                    if(!(m_uiMesh->getMPIRankGlobal()))
                        std::cout<<"executing step: "<<m_uiCurrentStep<<" dt: "<<m_uiT_h<<" rk_time : "<<m_uiCurrentTime<<std::endl;

                    char fPrefix[256];
                    sprintf(fPrefix,"%s_%d",m_uiFilePrefix,m_uiCurrentStep);

                    const char * fDataNames[]={"Time","Cyvle"};
                    const double fData[]={m_uiCurrentTime,(double)m_uiCurrentStep};
                    m_uiMesh->performGhostExchange(m_uiPrevU);
                    const char * pDataNames[]={"U(x,t)"};
                    const double * pData[]={&(*(m_uiPrevU.begin()))};
                    if(m_uiCurrentStep%IO_OUTPUT_FREQ==0)/*if((m_uiCurrentStep==0) || fabs((m_uiCurrentTime)-m_uiLastIOTimeValue)>IO_OUTPUT_GAP)*/
                    {
                        m_uiLastIOTimeValue=m_uiCurrentTime;
                        io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,1,pDataNames,pData);

                    }

                    //storeCheckPoint(m_uiFilePrefix);

                }

                if(m_uiCurrentStep%REMESH_TEST_FREQ==0)
                {

                    m_uiRefinedOctIDs.clear();
                    m_uiCoarsenOctIDs.clear();
                    m_uiMesh->performGhostExchange(m_uiPrevU);
                    m_uiMesh->unzip(&(*(m_uiPrevU.begin())),&(*(m_uiUZipVecIn.begin())));
                    const double* unzipVarPtr []={&(*(m_uiUZipVecIn.begin()))};
                    //isRefine=m_uiMesh->isReMesh(varPtr,refineVarIds,refineNumVars,m_uiWaveletTol);
                    isRefine=m_uiMesh->isReMeshUnzip((const double **)unzipVarPtr,refineVarIds,refineNumVars,waveletTolFunc);

                    if(isRefine)
                    {


                        /*pNodes=&(*(m_uiMesh->getAllElements().begin()));
                        numRefinedElements=0;
                        numCoarsenElements=0;
                        numNotChangedElements=0;
                        for(unsigned int ele=m_uiMesh->getElementLocalBegin();ele<m_uiMesh->getElementLocalEnd();ele++)
                        {
                            if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                            {
                                numRefinedElements++;
                            }
                            else if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                            {
                                numCoarsenElements++;
                                assert(pNodes[ele].getParent()==pNodes[ele+NUM_CHILDREN-1].getParent());
                                ele=ele+NUM_CHILDREN-1;

                            }else
                            {
                                assert((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                                numNotChangedElements++;
                            }

                        }

                        std::cout<<"Interation: "<<m_uiCurrentStep<<" rank: "<<m_uiMesh->getMPIRank()<<" isRefine: "<<isRefine<<" # split elements: "<<numRefinedElements<<" # coarsen elements: "<<numCoarsenElements<<" # not changed: "<<numNotChangedElements<<std::endl;*/

                        ot::Mesh* newMesh=m_uiMesh->ReMesh(m_uiLoadImbTol,m_uiSplitterFix);
                        unsigned int oldElements_g;
                        unsigned int newElements_g;
                        unsigned int oldElements=m_uiMesh->getNumLocalMeshElements();
                        unsigned int newElements=newMesh->getNumLocalMeshElements();

                        par::Mpi_Reduce(&oldElements,&oldElements_g,1,MPI_SUM,0,m_uiMesh->getMPIGlobalCommunicator());
                        par::Mpi_Reduce(&newElements,&newElements_g,1,MPI_SUM,0,newMesh->getMPIGlobalCommunicator());

                        if(!(m_uiMesh->getMPIRankGlobal()))std::cout<<"step : "<<m_uiCurrentStep<<" time : "<<m_uiCurrentTime<<" old mesh: "<<oldElements_g<<" new mesh: "<<newElements_g<<std::endl;


                        m_uiMesh->interGridTransfer(m_uiPrevU,newMesh);

                        newMesh->createVector(m_uiU);
                        newMesh->createVector(m_uiDxU);
                        newMesh->createVector(m_uiDyU);
                        newMesh->createVector(m_uiDzU);

                        newMesh->createVector(m_uiStage_1);
                        newMesh->createVector(m_uiStage_2);
                        newMesh->createVector(m_uiStage_3);
                        newMesh->createVector(m_uiStage_4);
                        newMesh->createVector(m_uiStage_5);
                        newMesh->createVector(m_uiStage_6);
                        newMesh->createVector(m_uiVecIm);

                        newMesh->createUnZippedVector(m_uiUZipVecIn,0.0);
                        newMesh->createUnZippedVector(m_uiUZipVecOut,0.0);

                        std::swap(m_uiMesh,newMesh);
                        delete newMesh;

                        if(m_uiMesh->getMPIRankGlobal()==0)std::cout<<"transfer completed"<<std::endl;


                    }
                }

                performSingleIteration();
                std::swap(m_uiU,m_uiPrevU);


            }

        }

        int RK45Transport::storeCheckPoint(const char * fNamePrefix)
        {

            unsigned int rank=m_uiMesh->getMPIRankGlobal();
            unsigned int npes=m_uiMesh->getMPICommSizeGlobal();
            char fName[256];
            const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
            sprintf(fName,"%s_octree_%d_%d_%d.oct",fNamePrefix,m_uiCurrentStep,rank,npes);
            io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

            unsigned int numVars=1;
            const char * varNames[]={"U"};
            const double * varPtr[]={&(*(m_uiPrevU.begin()))};

            for(unsigned int i=0;i<numVars;i++)
            {
                sprintf(fName,"%s_%s_%d_%d_%d.var",fNamePrefix,varNames[i],m_uiCurrentStep,rank,npes);
                io::checkpoint::writeVecToFile(fName,m_uiMesh,varPtr[i]);
            }


            if(!rank)
            {

                sprintf(fName,"%s_step_%d.cp",fNamePrefix,m_uiCurrentStep);
                std::ofstream outfile(fName);
                if(!outfile) {std::cout<<fName<<" file open failed "<<std::endl; return 1;}

                json checkPoint;
                checkPoint["DENDRO_RK45_TIME_BEGIN"]=m_uiTimeBegin;
                checkPoint["DENDRO_RK45_TIME_END"]=m_uiTimeEnd;

                checkPoint["DENDRO_RK45_ELEMENT_ORDER"]=m_uiOrder;

                checkPoint["DENDRO_RK45_TIME_CURRENT"]=m_uiCurrentTime;
                checkPoint["DENDRO_RK45_STEP_CURRENT"]=m_uiCurrentStep;
                checkPoint["DENDRO_RK45_TIME_STEP_SIZE"]=m_uiT_h;
                checkPoint["DENDRO_RK45_LAST_IO_TIME"]=m_uiLastIOTimeValue;

                checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"]=m_uiWaveletTol;
                checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"]=m_uiLoadImbTol;
                checkPoint["DENDRO_RK45_NUM_VARS"]=numVars; // number of variables to restore.

                outfile<<std::setw(4)<<checkPoint<<std::endl;
                outfile.close();

            }

            return 0;

        }


        int RK45Transport::restoreRK45Solver(const char * fNamePrefix, const unsigned int step,const MPI_Comm comm)
        {

            unsigned int numVars=0;
            std::vector<ot::TreeNode> octree;
            json checkPoint;

            int rank;
            int npes;
            m_uiComm=comm;
            MPI_Comm_rank(m_uiComm,&rank);
            MPI_Comm_size(m_uiComm,&npes);

            char fName[256];

            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp",fNamePrefix,step);
                std::ifstream infile(fName);
                if(!infile) {std::cout<<fName<<" file open failed "<<std::endl; return 1;}
                infile>>checkPoint;

                m_uiTimeBegin=checkPoint["DENDRO_RK45_TIME_BEGIN"];
                m_uiTimeEnd=checkPoint["DENDRO_RK45_TIME_END"];

                m_uiOrder=checkPoint["DENDRO_RK45_ELEMENT_ORDER"];
                m_uiCurrentTime=checkPoint["DENDRO_RK45_TIME_CURRENT"];
                m_uiCurrentStep=checkPoint["DENDRO_RK45_STEP_CURRENT"];
                m_uiT_h=checkPoint["DENDRO_RK45_TIME_STEP_SIZE"];
                m_uiLastIOTimeValue=checkPoint["DENDRO_RK45_LAST_IO_TIME"];

                m_uiWaveletTol=checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"];
                m_uiLoadImbTol=checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"];
                numVars=checkPoint["DENDRO_RK45_NUM_VARS"];
            }

            par::Mpi_Bcast(&m_uiTimeBegin,1,0,comm);
            par::Mpi_Bcast(&m_uiTimeEnd,1,0,comm);

            par::Mpi_Bcast(&m_uiCurrentTime,1,0,comm);
            par::Mpi_Bcast(&m_uiCurrentStep,1,0,comm);

            par::Mpi_Bcast(&m_uiT_h,1,0,comm);
            par::Mpi_Bcast(&m_uiLastIOTimeValue,1,0,comm);

            par::Mpi_Bcast(&m_uiWaveletTol,1,0,comm);
            par::Mpi_Bcast(&m_uiLoadImbTol,1,0,comm);

            par::Mpi_Bcast(&numVars,1,0,comm);
            par::Mpi_Bcast(&m_uiOrder,1,0,comm);
            par::Mpi_Bcast(&m_uiT_h,1,0,comm);


            sprintf(fName,"%s_octree_%d_%d_%d.oct",fNamePrefix,m_uiCurrentStep,rank,npes);
            io::checkpoint::readOctFromFile(fName,octree);
            assert(par::test::isUniqueAndSorted(octree,m_uiComm));

            ot::Mesh* newMesh=new ot::Mesh(octree,1,m_uiOrder,m_uiComm);


            newMesh->createVector(m_uiU);
            newMesh->createVector(m_uiPrevU);

            newMesh->createVector(m_uiStage_1);
            newMesh->createVector(m_uiStage_2);
            newMesh->createVector(m_uiStage_3);
            newMesh->createVector(m_uiStage_4);
            newMesh->createVector(m_uiStage_5);
            newMesh->createVector(m_uiStage_6);

            newMesh->createVector(m_uiRHS);
            newMesh->createVector(m_uiDxU); // partial x derivative
            newMesh->createVector(m_uiDyU); // partial y derivative
            newMesh->createVector(m_uiDzU); // partial z derivative
            newMesh->createVector(m_uiVecIm);

            const char * varNames[]={"U"};
            double * varPtr[]={&(*(m_uiPrevU.begin()))};


            for(unsigned int i=0;i<numVars;i++)
            {
                sprintf(fName,"%s_%s_%d_%d_%d.var",fNamePrefix,varNames[i],m_uiCurrentStep,rank,npes);
                io::checkpoint::readVecFromFile(fName,newMesh,varPtr[i]);
            }


            std::swap(m_uiMesh,newMesh);
            delete newMesh;

            unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
            unsigned int totalElems;

            par::Mpi_Reduce(&localSz,&totalElems,1,MPI_SUM,0,m_uiComm);


            if(!rank) std::cout<<" checkpoint at step : "<<step<<" restore successful "<<" restored mesh size: "<<totalElems<<std::endl;
            return 0;


        }


    } // end of namespace solver

}// end of namespace ode
