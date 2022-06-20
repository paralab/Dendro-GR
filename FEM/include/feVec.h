//
// Created by milinda on 10/30/18.
//
/**
 * @file feVec.h
 * @brief feVec.h abstract interface based on Dendro 4. feVec contains an abstract interface to
 * perform FEM accumilations for a given vector.
 * @author Milinda Fernando
 *
 * **/
#ifndef DENDRO_5_0_FEVEC_H
#define DENDRO_5_0_FEVEC_H

#include "oda.h"

#ifdef BUILD_WITH_PETSC
#include "petscdmda.h"
#endif

class feVec {

protected:
    /**@brief: pointer to OCT DA*/
    ot::DA* m_uiOctDA;

    /**@brief: type of the DA*/
    ot::DAType m_uiDaType;

    /**@brief problem domain min point*/
    Point m_uiPtMin;

    /**@brief problem domain max point*/
    Point m_uiPtMax;

#ifdef BUILD_WITH_PETSC
    /**@brief: petsc DM*/
    DM m_uiPETSC_DA;
#endif

public:

    /**@brief: feVec constructor
     * @par[in] daType: type of the DA
     * */
    feVec(ot::DA* da)
    {
        m_uiOctDA=da;
    }

    /**@brief deconstructor*/
    ~feVec()
    {

    }
    /**
     * @brief Evaluates the RHS of the PDE at specific points (for example evaluation at the quadrature points)
     * @param [out] out : function evaluated at specific points.
     * */
    //virtual void evalVec(VECType* out,double scale=1.0)=0;

    /**
     * @brief Evaluates the right hand side of the weak formulations.
     * Typically the mass matrix multiplied with the load function.
     * @param [in] in: Input vector (f)
     * @param [out] out : Output vector (Mf)
     * */
    virtual void computeVec(const VECType* in,VECType* out,double scale=1.0)=0;


    /**@brief set the problem dimension*/
    inline void setProblemDimensions(const Point& pt_min, const Point& pt_max)
    {
        m_uiPtMin=pt_min;
        m_uiPtMax=pt_max;
    }

    virtual void setPlaceholder(const double * v) { assert(false); }

#ifdef BUILD_WITH_PETSC
// all PETSC function should go here.

    /**
     * @brief Evaluates the right hand side of the weak formulations.
     * Typically the mass matrix multiplied with the load function.
     * @param [in] in: Input vector (f)
     * @param [out] out : Output vector (Mf)
     * */
    virtual void computeVec(const Vec& in,Vec& out,double scale=1.0)=0;

#endif


};

#endif //DENDRO_5_0_FEVEC_H
