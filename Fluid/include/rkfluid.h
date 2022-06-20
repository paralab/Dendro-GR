/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief rk4 solver for relativistic fluid equations.
*/
//

#ifndef DENDRO_5_FLUID_RK_FLUID_H
#define DENDRO_5_FLUID_RK_FLUID_H


#include "rk.h"
#include "fdCoefficient.h"
#include "oct2vtk.h"
#include "checkPoint.h"
#include "mesh.h"
#include <string>
#include <iostream>
#include "fluidUtils.h"
#include "fluidMath.h"
#include "parameters.h"
#include "fluid.h"
#include "rhs.h"
#include "meshTestUtils.h"


namespace ode
{
    namespace solver
    {

        // coefficients for RK 3 solver ============================
        static const double RK3_C[]={1.0/6.0,1.0/6.0,2.0/3.0};
        static const double RK3_T[]={1.0,1.0,1.0};
        static const double RK3_U[][3]={{0.0,0.0,0.0},
                                        {1.0,0.0,0.0},
                                        {1.0/4.0,1.0/4.0,0.0}};
        //===========================================================

        // coefficients for RK 4 solver ==============================

        static const double RK4_C[]={1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};

        static const double RK4_T[]={0,1.0/2.0,1.0/2.0,1.0};
        static const double RK4_U[]={0.0,1.0/2.0,1.0/2.0,1.0};
        //=============================================================


        // coefficients for RK 45 solver ==============================
        static const double RK_5_C[]={16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430, -9.0/50.0 ,2.0/55.0};
        static const double RK_4_C[]={25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4101.0, -1.0/5.0};

        static const double RK_T[]={1.0,1.0/4.0,3.0/8.0,12.0/13.0,1.0,0.5};

        static const double RK_U[][6]={{0.0,0.0,0.0,0.0,0.0,0.0},
                                       {1.0/4.0,0.0,0.0,0.0,0.0,0.0},
                                       {3.0/32.0,9.0/32.0,0.0,0.0,0.0,0.0},
                                       {1932.0/2197.0,-7200.0/2197.0,7296.0/2197.0,0.0,0.0,0.0},
                                       {439.0/216.0,-8.0,3680.0/513.0,-845.0/4104.0,0.0,0.0},
                                       {-8.0/27.0,2.0,-3544.0/2565.0,1859.0/4104.0,-11.0/40.0,0.0}};



        class RK_FLUID : public RK
        {

        private:
            
            /**@brief: type of the rk solver*/
            RKType m_uiRKType;
            
            /**@brief:number of rk stages*/
            unsigned int m_uiNumRKStages;

            /**@brief: list of pointers to the variable set indexed by enum VAR */
            DendroScalar ** m_uiEvolVar;

            /**@brief: previous time step solution*/
            DendroScalar ** m_uiEvolPrevVar;
            
            /**@brief: intermidiate variable for RK*/
            DendroScalar ** m_uiEvolVarIm;

            /**@brief: list of pointers to the primitive variable set.*/
            DendroScalar ** m_uiPrimVar;
            
            /**@brief: previous time step solution for primitive variables.*/
            DendroScalar ** m_uiPrimPrevVar;

            /**@brief: intermediate primitive var values*/
            DendroScalar ** m_uiPrimVarIm;

            /**@brief unzip vector for primitive variables. */
            DendroScalar ** m_uiPrimUnzipVar;

            /**@brief list of pointers to unzip version of the variables **/
            DendroScalar **m_uiEvolUnzipVar;

            /**@brief unzip rhs for each variable.*/
            DendroScalar **m_uiEvolUnzipVarRHS;

            /**@brief stage - value vector of RK method*/
            DendroScalar *** m_uiEvolStage;

            /**Send node bufferes for async communication*/
            DendroScalar ** m_uiSendNodeBuf;

            /**recv node bufferes for async communciation*/
            DendroScalar ** m_uiRecvNodeBuf;

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
            RK_FLUID(ot::Mesh *pMesh, double pTBegin, double pTEnd,double pTh,RKType rkType);

            /**@brief default destructor*/
            ~RK_FLUID();

            /** @brief: read parameters related to FLUID simulation and store them in static variables defined in parameters.h*/
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
            /** apply initial conditions*/
            void applyInitialConditions(DendroScalar ** zipIn, DendroScalar ** zipPrimIn);
            
            /** perform initial convergence for grid and reapply initial conditions*/
            void initialGridConverge();

            /** reallocates mpi resources if the mesh is changed, (need to be called during refmesing)*/
            void reallocateMPIResources();

            /** @brief: perform ghost exchange for all vars*/
            void performGhostExchangeVars(DendroScalar** zipIn,unsigned int numVars);

            /**@brief: performs the intergrid transfer (generic version)*/
            void intergridTransferVars(DendroScalar **& zipIn, const ot::Mesh* pnewMesh,unsigned int numVars);

            /**@brief: performs the intergrid transfer for primitive variables (this is a special case to enforce 4 to 3 velocity transform)*/
            void intergridTransferPrimitiveVars(DendroScalar **& zipPrimIn, ot::Mesh* pnewMesh);

            /**@brief unzip all the vars specified in VARS*/
            void unzipVars(DendroScalar ** zipIn , DendroScalar **uzipOut,unsigned int numVars);

            /**@brief unzip all the vars specified in VARS (generic version)*/
            void unzipVars_async(DendroScalar ** zipIn , DendroScalar **uzipOut,unsigned int numVars);

            /**@brief special function to unzip primitive vars with 3,4 velocity vector transformations*/
            void unzipPrimVars(DendroScalar ** zipIn, DendroScalar** unZipOut);

            /**@brief zip all the variables specified in VARS*/
            void zipVars(DendroScalar** uzipIn , DendroScalar** zipOut, unsigned int numVars);

            /**@brief write the solution to vtu file. */
            void writeToVTU(DendroScalar **evolZipVarIn, DendroScalar** primVars, DendroScalar ** constrZipVarIn, unsigned int numEvolVars,unsigned int numPrimVars,unsigned int numConstVars,const unsigned int * evolVarIndices,const unsigned int* primVarIndices, const unsigned int * constVarIndices, bool zslice);

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





#endif //SFCSORTBENCH_RK4FLUID_H
