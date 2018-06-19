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


#ifndef SFCSORTBENCH_RKMAXWELL_H
#define SFCSORTBENCH_RKMAXWELL_H

#include "rk.h"

namespace ode
{
    namespace solver
    {

        class RK45Maxwell : public RK
        {

            private :
                /// declaration of all the vector and scalar parameters specified in the above PDE system
                /** x component of the vector potential A where $B=\delta \times A$*/
                double * m_uiA_x;
                /** x component of the vector potential A where $B=\delta \times A$*/
                double * m_uiA_y;
                /** x component of the vector potential A where $B=\delta \times A$*/
                double * m_uiA_z;


                /** x component of the electric field E*/
                double * m_uiE_x;
                /** y component of the electric field E*/
                double * m_uiE_y;
                /** z component of the electric field E*/
                double * m_uiE_z;

                /**$\gamma=\partial_x A_x + \partial_y A_y + \partial_z A_z $ */
                double * m_uiGamma;

                /**scalar potential $\varphi$*/
                double * m_uiPhi;


            public :
            /**Default constructor*/
            RK45Maxwell(ot::Mesh * pMesh,double pTBegin,double pTEnd,double pTh);

            /**Destructor to release the memory for the allocated variables*/
            ~RK45Maxwell();

            /**Implementation of the base class time step function*/
            void performSingleInteration();

            /**Implementation of the base class function to apply boundary conditions. */
            void applyBoundaryConditions();



        };


    } // end of namespace solver

}//end of namespace ode


#endif //SFCSORTBENCH_RKMAXWELL_H
