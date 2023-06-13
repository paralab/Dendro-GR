//
// Created by milinda on 5/24/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains specific example of solving 3D transport equation using RK45 method
 *
 * $\partial_t U + b. grad(U)=f(x,t) $
 * $U(x,0)=g(x)=sin(x_i)sin(x_j)sin(x_k)$
 *
*/
//

#ifndef SFCSORTBENCH_RKTRANSPORT_H
#define SFCSORTBENCH_RKTRANSPORT_H

#include "rk.h"
#include "fdCoefficient.h"
#include "oct2vtk.h"
#include "checkPoint.h"
#include "rkTransportUtils.h"
#include <string>
#include <iostream>



#define IO_OUTPUT_FREQ 10
#define REMESH_TEST_FREQ 10  // test for remeshing every 10 steps.
#define IO_OUTPUT_GAP 1 // output IO after each 1 second of traversal.




namespace ode
{
    namespace solver
    {

        class RK45Transport : public RK
        {

        private :
            /** current values of  U(x,t)*/
            std::vector<double> m_uiU;

            /** previous values of  U(x,t)*/
            std::vector<double> m_uiPrevU;

            /** intermidiate vector to comput stage values*/
            std::vector<double> m_uiVecIm;

            /** stage 1- value vector of RK45 method*/
            std::vector<double> m_uiStage_1;

            /** stage 2- value vector of RK45 method*/
            std::vector<double> m_uiStage_2;

            /** stage 3- value vector of RK45 method*/
            std::vector<double> m_uiStage_3;

            /** stage 4- value vector of RK45 method*/
            std::vector<double> m_uiStage_4;

            /** stage 5- value vector of RK45 method*/
            std::vector<double> m_uiStage_5;

            /** stage 6- value vector of RK45 method*/
            std::vector<double> m_uiStage_6;

            /** variable to store the rhs of the function. */
            std::vector<double> m_uiRHS;


            /**partial derivative respect to x direction*/
            std::vector<double>m_uiDxU;
            /**partial derivative respect to x direction*/
            std::vector<double>m_uiDyU;
            /**partial derivative respect to x direction*/
            std::vector<double>m_uiDzU;


            /**unZip vector In */
            std::vector<double> m_uiUZipVecIn;

            /**unZip vector Out */
            std::vector<double> m_uiUZipVecOut;



            /** x derivative 4th order stencils*/
            Stencil<double,5,2> m_uiDxOrder4_centered;
            Stencil<double,5,4> m_uiDxOrder4_backward;
            Stencil<double,5,0> m_uiDxOrder4_forward;

            /** x derivative 4th order stencils*/
            Stencil<double,5,2> m_uiDyOrder4_centered;
            Stencil<double,5,4> m_uiDyOrder4_backward;
            Stencil<double,5,0> m_uiDyOrder4_forward;

            /** x derivative 4th order stencils*/
            Stencil<double,5,2> m_uiDzOrder4_centered;
            Stencil<double,5,4> m_uiDzOrder4_backward;
            Stencil<double,5,0> m_uiDzOrder4_forward;


            /** constant vector */
            double m_uiCoef_b[3];

            /**initial condition (value of U(x,0))*/
            std::function<double(double,double,double)> m_uiFunc_g;

            /**right hand side of the transport equation*/
            std::function<double(double,double,double,double)> m_uiFunc_f;

            /** right hand side function passing time step as a constant. */
            std::function<double(double,double,double)> m_uiFunc_f1;

            /** prefix name to output the vtu files.*/
            char* m_uiFilePrefix;

            /** wavelet tolerance used for refinement*/
            double m_uiWaveletTol;

            /**load imbalance tolerance used for flexible partioning (defaule value set to 0.1) use setLDImbParameters function to change this*/
            double m_uiLoadImbTol;

            /** spliiter fix value*/
            unsigned int m_uiSplitterFix;

            /** last io time value. */
            double m_uiLastIOTimeValue;



        public :
            /**@brief Default constructor
             * @param[in] pMesh: mesh that we will use to solve transport equation.
             * @param[in] pTBegin: begin time of time step.
             * @param[in] pTEnd: end time of time step.
             * */
            RK45Transport(ot::Mesh * pMesh,double pTBegin,double pTEnd,double pTh);

            /**@breif: Destructor to release the memory for the allocated variables
             * */
            ~RK45Transport();

            /**@brief: set parameters for transport equation
             * @param [in] constant vector b as given in transport equation. This should be an array size of 3.
             * @param [in] g: initial value when t=pTBegin.
             * */
            void setParameters(double *b, std::function<double(double,double,double)>g , std::function<double (double,double,double,double)> f ,const char * fprefix,double wltTol );

            /**
             * @brief Used to set the load balancing parameters for the solver
             * @param [in] ldTol: load imbalance toelrance
             * @param [in] sfK: splitter fix value
             * */
            void setLDImbParameters(double ldTol=0.1,unsigned int sfK=2);

            /** apply intial conditions*/
            void applyInitialConditions();

            /**@brief: Implementation of the base class time step function*/
            void performSingleIteration();

            /**@brief: Implementation of the base class function to apply boundary conditions. */
            void applyBoundaryConditions();

            /**@brief: apply boundary conditions for the each stage of rk45*/
            template <typename T>
            void applyBoundaryConditions(const double time,T* vec);

            /**@brief: starts the rk-45 solver. */
            void rkSolve();

            /**@brief: stores the all the variables that is required to restore the rk45 solver at a given stage.
             * @param[in] fNamePrefix: checkpoint file pre-fix name.
             * */
            int storeCheckPoint(const char * fNamePrefix);

            /** @brief: restore rk45 solver from a given checkpoint. This will overwrite the parameters given in the original constructor
             *  @param[in]fNamePrefix: checkpoint file pre-fix name.
             *  @param[in]step: step number which needs to be restored.
             *  @param[in]comm: MPI communicator.
             * */
            int restoreRK45Solver(const char * fNamePrefix, const unsigned int step,const MPI_Comm comm);



        };




    } // end of namespace solver

}//end of namespace ode




namespace ode
{
    namespace solver
    {

        template <typename T>
        void RK45Transport::applyBoundaryConditions(const double time,T*vec)
        {
            // enforce boundary conditions to g(x-bt)
            unsigned int x,y,z,sz;
            const ot::TreeNode* pNodes=&(*(m_uiMesh->getAllElements().begin()));

            const double grid_min=0;
            const double grid_max=1u<<(m_uiMaxDepth);
            const double d_min=-M_PI;
            const double d_max=M_PI;
            const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
            const unsigned int eleOrder=m_uiMesh->getElementOrder();
            const unsigned int *e2n=&(*(m_uiMesh->getE2NMapping().begin()));
            const unsigned int *e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
            double x1,y1,z1;
            double value=1.0;
            unsigned int owner,ii_x,jj_y,kk_z;

            for(unsigned int ele=m_uiMesh->getElementLocalBegin();ele<m_uiMesh->getElementLocalEnd();ele++)
            {

                if((pNodes[ele].minX()==grid_min)  || (pNodes[ele].minY()==grid_min)  || ((pNodes[ele].minZ()==grid_min) ))
                { // current element is a boundary element.

                    for(unsigned int k=0;k<(eleOrder+1);k++)
                        for(unsigned int j=0;j<(eleOrder+1);j++)
                            for(unsigned int i=0;i<(eleOrder+1);i++)
                            {
                                if((m_uiMesh->isNodeLocal(ele,i,j,k)))
                                {
                                    m_uiMesh->dg2eijk(e2n_dg[ele*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i],owner,ii_x,jj_y,kk_z);
                                    sz=1u<<(m_uiMaxDepth-pNodes[owner].getLevel());
                                    x=pNodes[owner].getX()+ ii_x*(sz/eleOrder);
                                    y=pNodes[owner].getY()+ jj_y*(sz/eleOrder);
                                    z=pNodes[owner].getZ()+ kk_z*(sz/eleOrder);
                                    if(((x==grid_min)|| (y==grid_min) ||  (z==grid_min) ))
                                    {
                                        assert(!m_uiMesh->isNodeHanging(owner,ii_x,jj_y,kk_z));
                                        x1=((x-grid_min)/(grid_max-grid_min))*(d_max-d_min)+d_min;
                                        y1=((y-grid_min)/(grid_max-grid_min))*(d_max-d_min)+d_min;
                                        z1=((z-grid_min)/(grid_max-grid_min))*(d_max-d_min)+d_min;

                                        /*((x1-m_uiCoef_b[0]*m_uiCurrentTime)<d_min)? value=value*(-1)*sin(2*d_min-(x1-m_uiCoef_b[0]*m_uiCurrentTime)) :value=value*sin((x1-m_uiCoef_b[0]*m_uiCurrentTime));
                                        ((y1-m_uiCoef_b[1]*m_uiCurrentTime)<d_min)? value=value*(-1)*sin(2*d_min-(y1-m_uiCoef_b[1]*m_uiCurrentTime)) :value=value*sin((y1-m_uiCoef_b[1]*m_uiCurrentTime));
                                        ((z1-m_uiCoef_b[2]*m_uiCurrentTime)<d_min)? value=value*(-1)*sin(2*d_min-(z1-m_uiCoef_b[2]*m_uiCurrentTime)) :value=value*sin((z1-m_uiCoef_b[2]*m_uiCurrentTime));*/

                                        value=sin(x1-m_uiCoef_b[0]*time)*sin(y1-m_uiCoef_b[1]*time)*sin(z1-m_uiCoef_b[2]*time);
                                        //if(!m_uiMesh->getMPIRank() && (fabs(value- m_uiU[e2n[ele*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]])>1e-2)) std::cout<<" rank: "<<m_uiMesh->getMPIRank()<<"bdy value: "<<value<<" computed value: "<< m_uiU[e2n[ele*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]]<<" i,j,k: ("<<i<<" , "<<j<<" ,"<<k<<" ) "<<std::endl;

                                        vec[e2n[ele*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]]=value;//m_uiFunc_g((x-m_uiCoef_b[0]*m_uiCurrentTime),(y-m_uiCoef_b[1]*m_uiCurrentTime),(z-m_uiCoef_b[2]*m_uiCurrentTime));
                                    }
                                }

                            }

                }
            }

        }


    } // end of namespace solver

} // end of namespace ode



#endif //SFCSORTBENCH_RKTRANSPORT_H
