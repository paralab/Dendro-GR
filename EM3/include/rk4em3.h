//
// Created by milinda on 12/1/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief rk4 solver for em3 equations.
*/
//

#ifndef SFCSORTBENCH_RK4EM3_H
#define SFCSORTBENCH_RK4EM3_H


#include "rk.h"
#include "fdCoefficient.h"
#include "oct2vtk.h"
#include "checkPoint.h"
#include "mesh.h"
#include <string>
#include <iostream>
#include "em3Utils.h"
#include "parameters.h"
#include "em3.h"
#include "rhs.h"
#include "meshTestUtils.h"
#include "mathMeshUtils.h"
#include "physcon.h"

static const double RK4_C[]={1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};

static const double RK4_T[]={0,1.0/2.0,1.0/2.0,1.0};
static const double RK4_U[]={0.0,1.0/2.0,1.0/2.0,1.0};

namespace ode
{
    namespace solver
    {

        class RK4_EM3 : public RK
        {

        private:

            // variables for EM3 formulation.

            /**@brief: list of pointers to the variable set indexed by enum VAR */
            double ** m_uiVar;

            /**@brief: previous time step solution*/
            double ** m_uiPrevVar;

            /**@brief: intermediate variable for RK*/
            double ** m_uiVarIm;

            /**@brief list of pointers to unzip version of the variables **/
            double **m_uiUnzipVar;

            /**@brief unzip rhs for each variable.*/
            double **m_uiUnzipVarRHS;

            /** stage - value vector of RK45 method*/
            double *** m_uiStage;

            /** @brief zipped version physical constraint equations*/
            double ** m_uiConstraintVars;

            /** @brief unzip physical constrint vars*/
            double ** m_uiUnzipConstraintVars;

            /**Send node buffers for async (send) communication*/
            double ** m_uiSendNodeBuf;

            /**recv node buffers for async (receive) communciation*/
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
            RK4_EM3(ot::Mesh *pMesh, double pTBegin, double pTEnd,double pTh);

            /**@brief default destructor*/
            ~RK4_EM3();

            /** @brief: read parameters related to EM3 simulation and store them in static variables defined in parameters.h*/
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

            /**performs initial grid convergence until the mesh converges to initial data. */
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
            //void writeToVTU(double **evolZipVarIn, double ** constrZipVarIn, unsigned int numEvolVars,unsigned int numConstVars,const unsigned int * evolVarIndices, const unsigned int * constVarIndices);
            void writeToVTU(double **evolZipVarIn, double ** constrZipVarIn, unsigned int numEvolVars,unsigned int numConstVars,const unsigned int * evolVarIndices, const unsigned int * constVarIndices, bool zslice);

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





#endif //SFCSORTBENCH_RK4EM3_H
