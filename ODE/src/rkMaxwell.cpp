//
// Created by milinda on 1/30/17.
//

/**
 * @author Milinda Fernando
 * @School of Computing, University of Utah
 * @brief Constains the implementation of Maxwell equations, with the stability properties in numerical relativity.
 *
 * Reference : https://arxiv.org/abs/gr-qc/0201051
 *
 * With the numerical stability properties we can rewrite the Maxwell's equations as follows.
 * $\partial_t A_i =-E_i -D_i \varphi $
 * $\partial_t E_i =-D^jD_jA_i + D_i\gamma  -4\piJ_i $
 * $\partial_t \gamma=-4\pi\ \rho _e -D^jD_j \varphi $
 * $\partial_t  \varphi =-D_iA^i -4\pi\rho_e$
 *
 * */

#include <rkMaxwell.h>

namespace ode
{

    namespace solver
    {

        RK45Maxwell::RK45Maxwell(ot::Mesh *pMesh,double pTBegin, double pTEnd,double pTh):RK (pMesh,pTBegin,pTEnd,pTh)
        {

            m_uiMesh->createVector(m_uiA_x);
            m_uiMesh->createVector(m_uiA_y);
            m_uiMesh->createVector(m_uiA_z);

            m_uiMesh->createVector(m_uiE_x);
            m_uiMesh->createVector(m_uiE_y);
            m_uiMesh->createVector(m_uiE_z);

            m_uiMesh->createVector(m_uiGamma);
            m_uiMesh->createVector(m_uiPhi);

        }


        void RK45Maxwell::performSingleInteration()
        {
            unsigned int neighList[NUM_CHILDREN];
            for(m_uiMesh->init<ot::WaveletDA::INDEPENDENT>();m_uiMesh->nextAvailable<ot::WaveletDA::INDEPENDENT>();m_uiMesh->next<ot::WaveletDA::INDEPENDENT>())
            {
                m_uiMesh->currentElementNeighbourIndexList(neighList);


            }


        }

        void RK45Maxwell::applyBoundaryConditions()
        {

            // Write code to apply the boundary conditions.


        }

        RK45Maxwell::~RK45Maxwell()
        {

            if(m_uiA_x!=NULL) delete [] m_uiA_x;
            if(m_uiA_y!=NULL) delete [] m_uiA_y;
            if(m_uiA_z!=NULL) delete [] m_uiA_z;

            if(m_uiE_x!=NULL) delete [] m_uiE_x;
            if(m_uiE_y!=NULL) delete [] m_uiE_y;
            if(m_uiE_z!=NULL) delete [] m_uiE_z;

            if(m_uiPhi!=NULL) delete [] m_uiPhi;
            if(m_uiGamma!=NULL) delete [] m_uiGamma;


        }



    }// end of namespace solver

}//end of namespace ode