//
// Created by milinda on 10/30/18.
//

/**
 * @brief class that derived from abstract class feMat
 * RHS computation of the weak formulation
 * */

#ifndef DENDRO_5_0_FEVECTOR_H
#define DENDRO_5_0_FEVECTOR_H

#include "feVec.h"

template <typename T>
class feVector : public feVec {

protected:

    /**@brief number of unknowns */
    unsigned int m_uiDof;

    /**@brief element nodal vec in */
    VECType * m_uiEleVecIn;

    /***@brief element nodal vecOut */
    VECType * m_uiEleVecOut;

    /** elemental coordinates */
    double * m_uiEleCoords;


public:
    /**
     * @brief constructs an FEM stiffness matrix class.
     * @param[in] da: octree DA
     * */
    feVector(ot::DA* da,unsigned int dof=1);

    ~feVector();

    /**
     * @brief Evaluates the RHS of the PDE at specific points (for example evaluation at the quadrature points)
     * @param [out] out : function evaluated at specific points.
     * */
    //virtual void evalVec(VECType* out,double scale=1.0);


    /**
     * @brief Evaluates the right hand side of the weak formulations.
     * Typically the mass matrix multiplied with the load function.
     * @param [in] in: Input vector (f)
     * @param [out] out : Output vector (Mf)
     * */
    virtual void computeVec(const VECType* in,VECType* out,double scale=1.0);


    /**@brief evalVec for the elemental vec*/
    //virtual void elementalEvalVec(VECType* out,double scale=1.0)=0;

    /**@brief elemental compute vec which evaluate the elemental RHS of the weak formulation
     * */
    virtual void elementalComputVec(const VECType* in,VECType* out,double* coords=NULL,double scale=1.0)=0;

#ifdef BUILD_WITH_PETSC

    /**
     * @brief Evaluates the right hand side of the weak formulations.
     * Typically the mass matrix multiplied with the load function.
     * @param [in] in: Input vector (f)
     * @param [out] out : Output vector (Mf)
     * */
    virtual void computeVec(const Vec& in,Vec& out,double scale=1.0);

#endif


    T& asLeaf() { return static_cast<T&>(*this);}

    bool preComputeVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().preComputeVec(in,out,scale);
    }

    bool postComputeVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().postComputeVec(in,out,scale);
    }

    bool preEvalVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().preEvalVec(in,out,scale);
    }

    bool postEvalVec(const VECType* in, VECType* out,double scale=1.0) {
        return asLeaf().postEvalVec(in,out,scale);
    }

};

template <typename T>
feVector<T>::feVector(ot::DA *da,unsigned int dof) : feVec(da)
{
    m_uiDof=dof;
    const unsigned int nPe=m_uiOctDA->getNumNodesPerElement();
    m_uiEleVecIn = new  VECType[m_uiDof*nPe];
    m_uiEleVecOut = new VECType[m_uiDof*nPe];

    m_uiEleCoords= new double[m_uiDim*nPe];
}

template <typename T>
feVector<T>::~feVector()
{
    delete [] m_uiEleVecIn;
    delete [] m_uiEleVecOut;

    delete [] m_uiEleCoords;
}


template <typename T>
void feVector<T>::computeVec(const VECType* in,VECType* out,double scale)
{


    VECType* _in=NULL;
    VECType* _out=NULL;

    if(!(m_uiOctDA->isActive()))
        return;

    preComputeVec(in,out,scale);

    m_uiOctDA->nodalVecToGhostedNodal(in,_in,false,m_uiDof);
    m_uiOctDA->createVector(_out,false,true,m_uiDof);

    VECType * val=new VECType[m_uiDof];
    for(unsigned int var=0;var<m_uiDof;var++)
        val[var]=(VECType)0;

    m_uiOctDA->setVectorByScalar(_out,val,false,true,m_uiDof);

    delete [] val;


    const ot::Mesh * pMesh = m_uiOctDA->getMesh();
    const unsigned int nPe = pMesh->getNumNodesPerElement();
    double * qMat = new double[nPe*nPe];
    const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
    const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
    
    m_uiOctDA->readFromGhostBegin(_in, m_uiDof);
    const unsigned int totalNodalSize = m_uiOctDA->getTotalNodalSz();
    for (m_uiOctDA->init<ot::DA_FLAGS::INDEPENDENT>(); m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::INDEPENDENT>(); m_uiOctDA->next<ot::DA_FLAGS::INDEPENDENT>()) {

        m_uiOctDA->getElementNodalValues(_in, m_uiEleVecIn, m_uiOctDA->curr(), m_uiDof);
        const unsigned int currentId = m_uiOctDA->curr();
        pMesh->getElementQMat(m_uiOctDA->curr(),qMat,true);
        m_uiOctDA->getElementalCoords(m_uiOctDA->curr(), m_uiEleCoords);
        elementalComputVec(m_uiEleVecIn, m_uiEleVecOut, m_uiEleCoords, scale);
        for (unsigned int dof = 0; dof < m_uiDof; dof++) {
          for (unsigned int i = 0; i < nPe; i++) {
            VECType ss = (VECType) 0;
            for (unsigned int j = 0; j < nPe; j++) {
              ss += qMat[j * nPe + i] * m_uiEleVecOut[dof*nPe + j];
            }
            _out[dof*totalNodalSize + e2n_cg[currentId * nPe + i]] += (VECType) ss;

          }
        }

    }

    m_uiOctDA->readFromGhostEnd(_in, m_uiDof);

    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd = pMesh -> getElementLocalEnd();

    for (m_uiOctDA->init<ot::DA_FLAGS::W_DEPENDENT>(); m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::W_DEPENDENT>(); m_uiOctDA->next<ot::DA_FLAGS::W_DEPENDENT>()) {

        // temporary fix to skip ghost writable.     
        if( m_uiOctDA->curr()< eleLocalBegin || m_uiOctDA->curr()>=eleLocalEnd )
            continue;

        m_uiOctDA->getElementNodalValues(_in, m_uiEleVecIn, m_uiOctDA->curr(), m_uiDof);
        const unsigned int currentId = m_uiOctDA->curr();
        pMesh->getElementQMat(m_uiOctDA->curr(),qMat,true);
        m_uiOctDA->getElementalCoords(m_uiOctDA->curr(), m_uiEleCoords);
        elementalComputVec(m_uiEleVecIn, m_uiEleVecOut, m_uiEleCoords, scale);

        for (unsigned int dof = 0; dof < m_uiDof; dof++) {
          for (unsigned int i = 0; i < nPe; i++) {
            VECType ss = (VECType) 0;
            for (unsigned int j = 0; j < nPe; j++) {
              ss += qMat[j * nPe + i] * m_uiEleVecOut[dof*nPe + j];
            }
            _out[dof*totalNodalSize + e2n_cg[currentId * nPe + i]] += (VECType) ss;

          }
        }

    }


    delete [] qMat;

    m_uiOctDA->writeToGhostsBegin(_out,m_uiDof);
    m_uiOctDA->writeToGhostsEnd(_out,ot::DA_FLAGS::WriteMode::ADD_VALUES,m_uiDof);

    m_uiOctDA->ghostedNodalToNodalVec(_out,out,true,m_uiDof);

    m_uiOctDA->destroyVector(_in);
    m_uiOctDA->destroyVector(_out);

    postComputeVec(in,out,scale);

    return;


}


#ifdef BUILD_WITH_PETSC


template <typename T>
void feVector<T>::computeVec(const Vec &in, Vec &out, double scale)
{
    const PetscScalar * inArry=NULL;
    PetscScalar * outArry=NULL;

    VecGetArrayRead(in,&inArry);
    VecGetArray(out,&outArry);

    computeVec(inArry,outArry,scale);

    VecRestoreArrayRead(in,&inArry);
    VecRestoreArray(out,&outArry);
}
#endif


#endif //DENDRO_5_0_FEVECTOR_H
