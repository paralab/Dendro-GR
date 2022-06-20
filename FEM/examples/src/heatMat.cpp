//
// Created by milinda on 11/21/18.
//

#include "heatMat.h"

HeatEq::HeatMat::HeatMat(ot::DA* da,unsigned int dof) : feMatrix(da,dof)
{
    const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
    imV1=new double[nPe];
    imV2=new double[nPe];

    Qx=new double[nPe];
    Qy=new double[nPe];
    Qz=new double[nPe];

}

HeatEq::HeatMat::~HeatMat()
{

    delete [] imV1;
    delete [] imV2;

    delete [] Qx;
    delete [] Qy;
    delete [] Qz;

    imV1=NULL;
    imV2=NULL;

    Qx=NULL;
    Qy=NULL;
    Qz=NULL;


}

void HeatEq::HeatMat::elementalMatVec(const VECType* in,VECType* out, double*coords,double scale)
{

    const RefElement* refEl=m_uiOctDA->getReferenceElement();

    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * DgT=refEl->getDgT1d();
    const double * W1d=refEl->getWgq();

    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
    const unsigned int nrp=eleOrder+1;

    Point eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
    Point eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);

    const double refElSz=refEl->getElementSz();
    //x derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Dg,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qx);

    //y derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Dg,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qy);

    //z derivative
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Dg,imV2,Qz);


    const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
    const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
    const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());


    const double Jx = 1.0/(refElSz/(double (szX)));
    const double Jy = 1.0/(refElSz/(double (szY)));
    const double Jz = 1.0/(refElSz/(double (szZ)));

    //std::cout<<"Stifness:  elem: "<<elem<<" ele Sz: "<<(elem.maxX()-elem.minX())<<" szX: "<<szX<<" Jx: "<<Jx<<" J: "<<(Jx*Jy*Jz)<<std::endl;

    for(unsigned int k=0;k<(eleOrder+1);k++)
        for(unsigned int j=0;j<(eleOrder+1);j++)
            for(unsigned int i=0;i<(eleOrder+1);i++)
            {
                Qx[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jy*Jz)/Jx)*W1d[i]*W1d[j]*W1d[k]);
                Qy[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jz)/Jy)*W1d[i]*W1d[j]*W1d[k]);
                Qz[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jy)/Jz)*W1d[i]*W1d[j]*W1d[k]);
            }



    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,DgT,Qx,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qx);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qy,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,DgT,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qy);

    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qz,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,DgT,imV2,Qz);

    for(unsigned int i=0;i<nPe;i++)
        out[i]=Qx[i]+Qy[i]+Qz[i];
}

bool HeatEq::HeatMat::preMatVec(const VECType* in,VECType* out,double scale)
{
    // apply boundary conditions.
    std::vector<unsigned int> bdyIndex;
    std::vector<double> bdyCoords;

    m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;



}

bool HeatEq::HeatMat::postMatVec(const VECType* in,VECType* out,double scale) {

    // apply boundary conditions.
    std::vector<unsigned int> bdyIndex;
    std::vector<double> bdyCoords;

    m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;
    
    return true;

}


double HeatEq::HeatMat::gridX_to_X(double x)
{
    double Rg_x=((1u<<m_uiMaxDepth)-0);
    return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
}

double HeatEq::HeatMat::gridY_to_Y(double y)
{
    double Rg_y=((1u<<m_uiMaxDepth)-0);
    return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
}


double HeatEq::HeatMat::gridZ_to_Z(double z)
{
    double Rg_z=((1u<<m_uiMaxDepth)-0);
    return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
}

int HeatEq::HeatMat::cgSolve(double * x ,double * b,int max_iter, double& tol,unsigned int var)
{
    double resid,alpha,beta,rho,rho_1;
    int status=1; // 0 indicates it has solved the system within the specified max_iter, 1 otherwise.

    const unsigned int local_dof=m_uiOctDA->getLocalNodalSz();

    MPI_Comm globalComm=m_uiOctDA->getGlobalComm();

    if(m_uiOctDA->isActive())
    {

        int activeRank=m_uiOctDA->getRankActive();
        int activeNpes=m_uiOctDA->getNpesActive();

        MPI_Comm activeComm=m_uiOctDA->getCommActive();

        double* p;
        double* z;
        double* q;
        double* Ax;
        double* Ap;
        double* r0;
        double* r1;

        m_uiOctDA->createVector(p);
        m_uiOctDA->createVector(z);
        m_uiOctDA->createVector(q);

        m_uiOctDA->createVector(Ax);
        m_uiOctDA->createVector(Ap);
        m_uiOctDA->createVector(r0);
        m_uiOctDA->createVector(r1);

        double normb = normLInfty(b,local_dof,activeComm);
        par::Mpi_Bcast(&normb,1,0,activeComm);

        if(!activeRank)
            std::cout<<"normb = "<<normb<<std::endl;

        matVec(x,Ax);

        /*char fPrefix[256];
        sprintf(fPrefix,"%s_%d","cg",0);
        const char * varNames[]={"U"};
        const double * var[]={Ax};
        io::vtk::mesh2vtuFine(mesh,fPrefix,0,NULL,NULL,1,varNames,var);
*/
        for(unsigned int i=0;i<local_dof;i++)
        {
            r0[i]=b[i]-Ax[i];
            p[i]=r0[i];
        }


        if (normb == 0.0)
            normb = 1;

        double normr=normLInfty(r0,local_dof,activeComm);
        par::Mpi_Bcast(&normr,1,0,activeComm);
        if(!activeRank) std::cout<<"initial residual : "<<(normr/normb)<<std::endl;

        if ((resid = normr / normb) <= tol) {
            tol = resid;
            max_iter = 0;

            m_uiOctDA->destroyVector(p);
            m_uiOctDA->destroyVector(z);
            m_uiOctDA->destroyVector(q);

            m_uiOctDA->destroyVector(Ax);
            m_uiOctDA->destroyVector(Ap);
            m_uiOctDA->destroyVector(r0);
            m_uiOctDA->destroyVector(r1);

            status=0;
        }

        if(status!=0)
        {

            for(unsigned int i=1;i<=max_iter;i++)
            {

                matVec(p,Ap);

                alpha=(dot(r0,r0,local_dof,activeComm)/dot(p,Ap,local_dof,activeComm));
                par::Mpi_Bcast(&alpha,1,0,activeComm);

                //if(!activeRank) std::cout<<"rank: " <<activeRank<<" alpha: "<<alpha<<std::endl;
                for(unsigned int e=0;e<local_dof;e++)
                {
                    x[e]+=alpha*p[e];
                    r1[e]=r0[e]-alpha*Ap[e];
                }

                normr=normLInfty(r1,local_dof,activeComm);
                par::Mpi_Bcast(&normr,1,0,activeComm);

                if((!activeRank) && (i%10==0)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;

                if ((resid = normr / normb) <= tol) {

                    if((!activeRank)) std::cout<<" iteration : "<<i<<" residual : "<<resid<<std::endl;
                    tol = resid;
                    m_uiOctDA->destroyVector(p);
                    m_uiOctDA->destroyVector(z);
                    m_uiOctDA->destroyVector(q);

                    m_uiOctDA->destroyVector(Ax);
                    m_uiOctDA->destroyVector(Ap);
                    m_uiOctDA->destroyVector(r0);
                    m_uiOctDA->destroyVector(r1);

                    status=0;
                    break;
                }

                beta=(dot(r1,r1,local_dof,activeComm)/dot(r0,r0,local_dof,activeComm));
                par::Mpi_Bcast(&beta,1,0,activeComm);

                //if(!activeRank) std::cout<<"<r_1,r_1> : "<<dot(r1+nodeLocalBegin,r1+nodeLocalBegin,local_dof,activeComm)<<" <r_0,r_0>: "<<dot(r0+nodeLocalBegin,r0+nodeLocalBegin,local_dof,activeComm)<<" beta "<<beta<<std::endl;



                for(unsigned int e=0;e<local_dof;e++)
                {
                    p[e]=r1[e]+beta*p[e];
                    r0[e]=r1[e];
                }


            }

            if(status!=0)
            {
                tol = resid;
                m_uiOctDA->destroyVector(p);
                m_uiOctDA->destroyVector(z);
                m_uiOctDA->destroyVector(q);

                m_uiOctDA->destroyVector(Ax);
                m_uiOctDA->destroyVector(Ap);
                m_uiOctDA->destroyVector(r0);
                m_uiOctDA->destroyVector(r1);
                status=1;

            }



        }


    }


    // bcast act as a barrier for active and inactive meshes.
    par::Mpi_Bcast(&tol,1,0,globalComm);
    return status;
}
