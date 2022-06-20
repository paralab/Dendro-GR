//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Runge-Kutta 45 (RK45) method implementation for BSSN formulation of GR equations.
*/
//

#ifndef SFCSORTBENCH_RKBSSN_H
#define SFCSORTBENCH_RKBSSN_H


#include "rk.h"
#include "fdCoefficient.h"
#include "oct2vtk.h"
#include "checkPoint.h"
#include "mesh.h"
#include <string>
#include <iostream>
#include "grUtils.h"
#include "parameters.h"
#include "bssn_constraints.h"
#include "gr.h"
#include "rhs.h"
#include "physcon.h"
#include "psi4.h"
#include "test/meshTestUtils.h"
#include "TwoPunctures.h"

#ifdef BSSN_ENABLE_CUDA
    #include "rhs_cuda.cuh"
    #include "params_cu.h"
    #include "profile_gpu.h"
#endif


static const double RK_5_C[]={16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430, -9.0/50.0 ,2.0/55.0};
static const double RK_4_C[]={25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4101.0, -1.0/5.0};

static const double RK_T[]={1.0,1.0/4.0,3.0/8.0,12.0/13.0,1.0,0.5};

static const double RK_U[][6]={{0.0,0.0,0.0,0.0,0.0,0.0},
                               {1.0/4.0,0.0,0.0,0.0,0.0,0.0},
                               {3.0/32.0,9.0/32.0,0.0,0.0,0.0,0.0},
                               {1932.0/2197.0,-7200.0/2197.0,7296.0/2197.0,0.0,0.0,0.0},
                               {439.0/216.0,-8.0,3680.0/513.0,-845.0/4104.0,0.0,0.0},
                               {-8.0/27.0,2.0,-3544.0/2565.0,1859.0/4104.0,-11.0/40.0,0.0}};

namespace ode
{
   namespace solver
   {

           class RK45_BSSN : public RK
           {

           private:

               // variables for BSSN formulation.

               /**@brief: list of pointers to the variable set indexed by enum VAR */
               double ** m_uiVar;

               /**@brief: previous time step solution*/
               double ** m_uiPrevVar;

               /**@brief: intermidiate variable for RK*/
               double ** m_uiVarIm;

               /**@brief list of pointers to unzip version of the variables **/
               double **m_uiUnzipVar;

               /**@brief unzip rhs for each variable.*/
               double **m_uiUnzipVarRHS;

               /** stage - value vector of RK45 method*/
               double *** m_uiStage;

               /** @brief zipped version physical constrint equations*/
               double ** m_uiConstraintVars;

               /** @brief unzip physical constrint vars*/
               double ** m_uiUnzipConstraintVars;

               /**Send node bufferes for async communication*/
               double ** m_uiSendNodeBuf;

               /**recv node bufferes for async communciation*/
               double ** m_uiRecvNodeBuf;

               /**@brief mpi send reqs for evolution vars*/
               MPI_Request ** m_uiSendReqs;

               /**@brief mpi recv reqs for evolution vars*/
               MPI_Request ** m_uiRecvReqs;

               /**@brief mpi send status to sync on sends*/
               MPI_Status ** m_uiSendSts;

               /**@brief mpi recv status to sync on recv*/
               MPI_Status ** m_uiRecvSts;





           public:
               /**
                * @brief default constructor
                * @param[in] pMesh : pointer to the mesh.
                * @param[in] pTBegin: RK45 time begin
                * @param[in] pTEnd: RK45 time end
                * @param[in] pTh: times step size.
                * * */
               RK45_BSSN(ot::Mesh *pMesh, double pTBegin, double pTEnd,double pTh);

               /**@brief default destructor*/
               ~RK45_BSSN();

               /** @brief: read parameters related to BSSN simulation and store them in static variables defined in parameters.h*/
               void readConfigFile(const char * fName);

               /**@brief: starts the rk-45 solver. */
               void rkSolve();

               /** @brief: restore rk45 solver from a given checkpoint. This will overwrite the parameters given in the original constructor
                *  @param[in]fNamePrefix: checkpoint file pre-fix name.
                *  @param[in]step: step number which needs to be restored.
                *  @param[in]comm: MPI communicator.
                * */
               void restoreCheckPoint(const char * fNamePrefix,MPI_Comm comm);

           private:
               /** apply intial conditions*/
               void applyInitialConditions(double ** zipIn);

               /** refine based on the initial grid until it converges, */
               void initialGridConverge();

               /** reallocates mpi resources if the mesh is changed, (need to be called during refmesing)*/
               void reallocateMPIResources();

               /** @brief: perform ghost exchange for all vars*/
               void performGhostExchangeVars(double** zipIn);

               /**@brief: performs the intergrid transfer*/
               void intergridTransferVars(double **& zipIn, const ot::Mesh* pnewMesh);

               /**@brief unzip all the vars specified in VARS*/
               void unzipVars(double ** zipIn , double **uzipOut);

               /**@brief unzip all the vars specified in VARS*/
               void unzipVars_async(double ** zipIn , double **uzipOut);

               /**@brief zip all the variables specified in VARS*/
               void zipVars(double** uzipIn , double** zipOut);

               /**@brief write the solution to vtu file. */
               void writeToVTU(double **evolZipVarIn, double ** constrZipVarIn, unsigned int numEvolVars,unsigned int numConstVars,const unsigned int * evolVarIndices, const unsigned int * constVarIndices);

               /**@brief: Implementation of the base class time step function*/
               void performSingleIteration();

               /**@brief: Implementation of the base class function to apply boundary conditions. */
               void applyBoundaryConditions();

               /**@brief: stores the all the variables that is required to restore the rk45 solver at a given stage.
                * @param[in] fNamePrefix: checkpoint file pre-fix name.
                * */
               void storeCheckPoint(const char * fNamePrefix);


           };

    } // end of namespace solver
}// end of namespace ode



#endif //SFCSORTBENCH_RKBSSN_H
