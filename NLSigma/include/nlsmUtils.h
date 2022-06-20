//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for NLSM simulation.
*/


#ifndef SFCSORTBENCH_GRUTILS_H
#define SFCSORTBENCH_GRUTILS_H

#include "point.h"
#include "parameters.h"
#include "mesh.h"
#include "block.h"
#include "parUtils.h"
#include "json.hpp"
#include "dendroProfileParams.h"
#include "profile_params.h"


#define Rx (nlsm::NLSM_COMPD_MAX[0]-nlsm::NLSM_COMPD_MIN[0])
#define Ry (nlsm::NLSM_COMPD_MAX[1]-nlsm::NLSM_COMPD_MIN[1])
#define Rz (nlsm::NLSM_COMPD_MAX[2]-nlsm::NLSM_COMPD_MIN[2])

#define RgX (nlsm::NLSM_OCTREE_MAX[0]-nlsm::NLSM_OCTREE_MIN[0])
#define RgY (nlsm::NLSM_OCTREE_MAX[1]-nlsm::NLSM_OCTREE_MIN[1])
#define RgZ (nlsm::NLSM_OCTREE_MAX[2]-nlsm::NLSM_OCTREE_MIN[2])

#define GRIDX_TO_X(xg) (((Rx/RgX)*(xg-nlsm::NLSM_OCTREE_MIN[0]))+nlsm::NLSM_COMPD_MIN[0])
#define GRIDY_TO_Y(yg) (((Ry/RgY)*(yg-nlsm::NLSM_OCTREE_MIN[1]))+nlsm::NLSM_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg) (((Rz/RgZ)*(zg-nlsm::NLSM_OCTREE_MIN[2]))+nlsm::NLSM_COMPD_MIN[2])

#define X_TO_GRIDX(xc) (((RgX/Rx)*(xc-nlsm::NLSM_COMPD_MIN[0]))+nlsm::NLSM_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc) (((RgY/Ry)*(yc-nlsm::NLSM_COMPD_MIN[1]))+nlsm::NLSM_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc) (((RgZ/Rz)*(zc-nlsm::NLSM_COMPD_MIN[2]))+nlsm::NLSM_OCTREE_MIN[2])


using json = nlohmann::json;
namespace nlsm
{

/**
 * @brief internal variables needed for rk update.
 * */


 /**
  * @brief: Read the parameter file and initialize the variables in parameters.h file.
  * @param[in] fName: file name
  * @param[in] comm: MPI communicator.
  * */
  void readParamFile(const char * fName,MPI_Comm comm);

  /**
   * @brief Dump the read parameters to verify. 
   * 
   * @param sout stream to output read parameters. 
   * @param root rank to output read parameters
   * @param comm MPI communicator. 
   */
  void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm);

/**
 * @brief Initialize all the variables for a given point in space.
 * @param [in] coord: coordinates of the point.
 * @param [out] var: pointer to the list of variables, computed. var size should be (VAR::U_SYMAT5+1)
 * @note This function is taken from the old single core nlsm version.
 **/

 // Initial data
 void initData(const double xx1,const double yy1,const double zz1, double *var);

 /**
  * analytical solution based on the d'Alembert's formular.
  *
  * */
 void analyticalSol(const double xx1,const double yy1,const double zz1,const double t, double *var);

 /**
  * @brief: Generates block adaptive octree for the given binary blockhole problem.
  *
  * */

  void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm);

  /**
   * @brief wavelet tolerance as a function of space.
   * */
  double computeWTol(double x,double y,double z,double* hx);

  /**
   * @brief force refinement at the pulse. 
   */
  bool isRemeshForce(const ot::Mesh* pMesh, const double ** unzipVec, unsigned int vIndex, double refine_th, double coarsen_th, bool isOverwrite);

  unsigned int getEleWeight(const ot::TreeNode* pNode);

}// end of namespace nlsm



namespace nlsm
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

#endif //SFCSORTBENCH_GRUTILS_H
