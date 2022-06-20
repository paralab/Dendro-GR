//
// Created by milinda on 11/22/18.
//

#include "heatVec.h"

HeatEq::HeatVec::HeatVec(ot::DA* da,unsigned int dof) : feVector(da,dof)
{
    const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
    imV1=new double[nPe];
    imV2=new double[nPe];


}

HeatEq::HeatVec::~HeatVec()
{
    delete [] imV1;
    delete [] imV2;

    imV1=NULL;
    imV2=NULL;

}

void HeatEq::HeatVec::elementalComputVec(const VECType* in,VECType* out, double*coords,double scale)
{

    const RefElement* refEl=m_uiOctDA->getReferenceElement();
    const double * Q1d=refEl->getQ1d();
    const double * QT1d=refEl->getQT1d();
    const double * Dg=refEl->getDg1d();
    const double * W1d=refEl->getWgq();

    const unsigned int eleOrder=refEl->getOrder();
    const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
    const unsigned int nrp=eleOrder+1;

    Point eleMin(coords[0*m_uiDim+0],coords[0*m_uiDim+1],coords[0*m_uiDim+2]);
    Point eleMax(coords[(nPe-1)*m_uiDim+0],coords[(nPe-1)*m_uiDim+1],coords[(nPe-1)*m_uiDim+2]);


    const double refElSz=refEl->getElementSz();

    // interpolate to quadrature points.
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,out);



    const double szX=gridX_to_X(eleMax.x())-gridX_to_X(eleMin.x());
    const double szY=gridY_to_Y(eleMax.y())-gridY_to_Y(eleMin.y());
    const double szZ=gridZ_to_Z(eleMax.z())-gridZ_to_Z(eleMin.z());


    const double Jx = 1.0/(refElSz/(double (szX)));
    const double Jy = 1.0/(refElSz/(double (szY)));
    const double Jz = 1.0/(refElSz/(double (szZ)));


    //std::cout<<"Mass:  elem: "<<elem<<" ele Sz: "<<(elem.maxX()-elem.minX())<<" szX: "<<szX<<" Jx: "<<Jx<<" J: "<<(Jx*Jy*Jz)<<std::endl;

    for(unsigned int k=0;k<(eleOrder+1);k++)
        for(unsigned int j=0;j<(eleOrder+1);j++)
            for(unsigned int i=0;i<(eleOrder+1);i++)
                out[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=(Jx*Jy*Jz*W1d[i]*W1d[j]*W1d[k]);


    // apply transpose operator
    DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,out,imV1);
    DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
    DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,out);
}




bool HeatEq::HeatVec::preComputeVec(const VECType* in,VECType* out, double scale)
{

    // apply boundary conditions.
    std::vector<unsigned int> bdyIndex;
    std::vector<double> bdyCoords;

    m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;

    return true;


}

bool HeatEq::HeatVec::postComputeVec(const VECType* in,VECType* out, double scale) {

    // apply boundary conditions.
    std::vector<unsigned int> bdyIndex;
    std::vector<double> bdyCoords;

    m_uiOctDA->getOctreeBoundaryNodeIndices(bdyIndex,bdyCoords);

    for(unsigned int i=0;i<bdyIndex.size();i++)
        out[bdyIndex[i]]=0.0;
    
    return true;


}


double HeatEq::HeatVec::gridX_to_X(double x)
{
    double Rg_x=((1u<<m_uiMaxDepth)-0);
    return (((x)/(Rg_x))*((m_uiPtMax.x()-m_uiPtMin.x()))+m_uiPtMin.x());
}

double HeatEq::HeatVec::gridY_to_Y(double y)
{
    double Rg_y=((1u<<m_uiMaxDepth)-0);
    return (((y)/(Rg_y))*((m_uiPtMax.y()-m_uiPtMin.y()))+m_uiPtMin.y());
}


double HeatEq::HeatVec::gridZ_to_Z(double z)
{
    double Rg_z=((1u<<m_uiMaxDepth)-0);
    return (((z)/(Rg_z))*((m_uiPtMax.z()-m_uiPtMin.z()))+m_uiPtMin.z());
}
