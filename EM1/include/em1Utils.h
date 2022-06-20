//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for EM1 simulation.
*/


#pragma once 

#include <iostream>
#include "point.h"
#include "parameters.h"
#include "mesh.h"
#include "block.h"
#include "parUtils.h"
#include "json.hpp"
#include "dendroProfileParams.h"
#include "profile_params.h"



#define Rx (em1::EM1_COMPD_MAX[0]-em1::EM1_COMPD_MIN[0])
#define Ry (em1::EM1_COMPD_MAX[1]-em1::EM1_COMPD_MIN[1])
#define Rz (em1::EM1_COMPD_MAX[2]-em1::EM1_COMPD_MIN[2])

#define RgX (em1::EM1_OCTREE_MAX[0]-em1::EM1_OCTREE_MIN[0])
#define RgY (em1::EM1_OCTREE_MAX[1]-em1::EM1_OCTREE_MIN[1])
#define RgZ (em1::EM1_OCTREE_MAX[2]-em1::EM1_OCTREE_MIN[2])

#define GRIDX_TO_X(xg) (((Rx/RgX)*(xg-em1::EM1_OCTREE_MIN[0]))+em1::EM1_COMPD_MIN[0])
#define GRIDY_TO_Y(yg) (((Ry/RgY)*(yg-em1::EM1_OCTREE_MIN[1]))+em1::EM1_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg) (((Rz/RgZ)*(zg-em1::EM1_OCTREE_MIN[2]))+em1::EM1_COMPD_MIN[2])

#define X_TO_GRIDX(xc) (((RgX/Rx)*(xc-em1::EM1_COMPD_MIN[0]))+em1::EM1_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc) (((RgY/Ry)*(yc-em1::EM1_COMPD_MIN[1]))+em1::EM1_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc) (((RgZ/Rz)*(zc-em1::EM1_COMPD_MIN[2]))+em1::EM1_OCTREE_MIN[2])

// type of rk method 
enum RKType{RK3,RK4,RK45};  


using json = nlohmann::json;
namespace em1
{
    /**
     * @brief: Read the parameter file and initialize the variables in parameters.h file.
     * @param[in] fName: file name
     * @param[in] comm: MPI communicator.
     * */
    void readParamFile(const char * fName,MPI_Comm comm);

    /**
     * @brief dump read parameters from a specified rank
     * 
     * @param sout: ostream to dump
     * @param root: rank to dump from. 
     * @param comm: mpi communicator. 
     */
    void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm);

    /**
     * @brief Initialize all the variables for a given point in space.
     * @param [in] coord: coordinates of the point.
     * @param [out] var: pointer to the list of variables, computed. var size should be (VAR::U_SYMAT5+1)
     * @note This function is taken from the old single core em1 version.
     **/
    void initData(const double xx1,const double yy1,const double zz1, double *var);

    /**
     * @brief: analytical solution based on the d'Alembert's formula.
     *
     * */
    void analyticalSol(const double xx1,const double yy1,const double zz1,const double t, double *var);

    /**
     * @brief: Generates block adaptive octree for the given binary blockhole problem.
     * @param[out] tmpNodes: created octree tmpNodes
     * @param[in] pt_min: block min point
     * @param[in] pt_max: block max point
     * @param[in] regLev: regular grid level
     * @param[in] maxDepth: maximum refinement level. 
     * @param[in] comm: MPI communicator. 
     * */
    void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm);

    /**
     * @brief Compute the wavelet tolerance as a function of space. 
     * 
     * @param x : x coord. 
     * @param y : y coord
     * @param z : z coord
     * @param tol_min : min. tolerance value. 
     * @return double 
     */
    double computeWTol(double x,double y,double z,double* hx);

    
    /**
     * @brief force refinement at the pulse. 
    */
    bool isRemeshForce(const ot::Mesh* pMesh, const double ** unzipVec, unsigned int vIndex, double refine_th, double coarsen_th, bool isOverwrite);

    /**
     * @breif: Compute L2 constraint norms. 
     */
    template <typename T>
    double computeConstraintL2Norm(const T* constraintVec, const T* maskVec, unsigned int lbegin, unsigned int lend,MPI_Comm comm);

    /**
        * @breif: Compute L2 constraint norms. 
        */
    template <typename T>
    double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec,const T* maskVector,T maskthreshoold);

    /**
        * @breif write constraints to a file. 
        */
    template<typename T>
    double extractConstraints(const ot::Mesh* mesh,const T** constraintVar,const T* maskVec, double maskthreshoold,unsigned int timestep);


    void writeBLockToBinary(const double **unzipVarsRHS,unsigned int offset,const double *pmin, const double *pmax,double* bxMin,double * bxMax, const unsigned int *sz,unsigned int blkSz,double dxFactor,const char* fprefix);


    /**
        * @brief Performs elemental based artificial dissipation based on DG book by Hesthavan
        * @param[in] mesh: const pointer to the mesh. 
        * @param[in/out] zipVars: em1 variables to which dissipation needs to be appied.   
        * @param[in] numVars: number of variables
        * @param[in] nc: cutoff order should be 0,4 since we are using 4th order elements (during interpolations etc). 
        * @param[in] s: dissipation parameter s , 
        * @note nc, s increasing, will lower the dissipation. 
        */
    template<typename T>
    void artificial_dissipation(ot::Mesh * pMesh , T** zipVars, unsigned int numVars, unsigned int nc, unsigned int s,bool isGhostEx=false);


    unsigned int getOctantWeight(const ot::TreeNode* pNode);




}// end of namespace em1



namespace em1
{

    namespace timer
    {

        /**@brief initialize all the flop counters. */
        void initFlops();

        /**@brief clears the snapshot counter for time profiler variables*/
        void resetSnapshot();


        /**@brief reduce min mean max.
         * @param [in] stat: local time
         * @param [out] stat_g 0-min, 1-mean 2-max
        * */
       template<typename T>
       void computeOverallStats(T *stat, T *stat_g, MPI_Comm comm)
       {
           int rank,npes;
           MPI_Comm_size(comm,&npes);
           MPI_Comm_rank(comm,&rank);

           par::Mpi_Reduce(stat,stat_g,1,MPI_MIN,0,comm);
           par::Mpi_Reduce(stat,stat_g+1,1,MPI_SUM,0,comm);
           par::Mpi_Reduce(stat,stat_g+2,1,MPI_MAX,0,comm);
           stat_g[1]/=(npes);

       }


        /** @breif : printout the profile parameters. */
        void profileInfo(const char* filePrefix,const ot::Mesh* pMesh);

        /** @breif : printout the profile parameters (intermediate profile information). */
        void profileInfoIntermediate(const char* filePrefix,const ot::Mesh* pMesh,const unsigned int currentStep);


    }


}

