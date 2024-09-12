//
// Created by milinda on 7/25/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Contains utility functions for FLUID simulation.
 */
#ifndef DENDRO_5_FLUID_FLUIDUTILS_H
#define DENDRO_5_FLUID_FLUIDUTILS_H

#include <waveletAMR.h>

#include "block.h"
#include "dendroProfileParams.h"
#include "fluidMath.h"
#include "json.hpp"
#include "mesh.h"
#include "parUtils.h"
#include "parameters.h"
#include "point.h"
#include "profile_params.h"

#define Rx (fluid::FLUID_COMPD_MAX[0] - fluid::FLUID_COMPD_MIN[0])
#define Ry (fluid::FLUID_COMPD_MAX[1] - fluid::FLUID_COMPD_MIN[1])
#define Rz (fluid::FLUID_COMPD_MAX[2] - fluid::FLUID_COMPD_MIN[2])

#define RgX (fluid::FLUID_OCTREE_MAX[0] - fluid::FLUID_OCTREE_MIN[0])
#define RgY (fluid::FLUID_OCTREE_MAX[1] - fluid::FLUID_OCTREE_MIN[1])
#define RgZ (fluid::FLUID_OCTREE_MAX[2] - fluid::FLUID_OCTREE_MIN[2])

#define GRIDX_TO_X(xg)                                  \
    (((Rx / RgX) * (xg - fluid::FLUID_OCTREE_MIN[0])) + \
     fluid::FLUID_COMPD_MIN[0])
#define GRIDY_TO_Y(yg)                                  \
    (((Ry / RgY) * (yg - fluid::FLUID_OCTREE_MIN[1])) + \
     fluid::FLUID_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg)                                  \
    (((Rz / RgZ) * (zg - fluid::FLUID_OCTREE_MIN[2])) + \
     fluid::FLUID_COMPD_MIN[2])

#define X_TO_GRIDX(xc)                                 \
    (((RgX / Rx) * (xc - fluid::FLUID_COMPD_MIN[0])) + \
     fluid::FLUID_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc)                                 \
    (((RgY / Ry) * (yc - fluid::FLUID_COMPD_MIN[1])) + \
     fluid::FLUID_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc)                                 \
    (((RgZ / Rz) * (zc - fluid::FLUID_COMPD_MIN[2])) + \
     fluid::FLUID_OCTREE_MIN[2])

using json = nlohmann::json;
namespace fluid {
/**
 * @brief These variable indexes are based on the variables defined in rkFLUID.h
 * */
enum PVAR { V_RHO = 0, V_VX, V_VY, V_VZ, V_P, V_U1, V_U2, V_U3 };
enum VAR { U_D = 0, U_SX, U_SY, U_SZ, U_TAU };

static const char* FLUID_PRIM_VAR_NAMES[] = {"V_RHO", "V_VX", "V_VY", "V_VZ",
                                             "V_P",   "V_U1", "V_U2", "V_U3"};
static const char* FLUID_EVOL_VAR_NAMES[] = {"U_D", "U_SX", "U_SY", "U_SZ",
                                             "U_TAU"};
static const char* FLUID_CONS_VAR_NAMES[] = {"C_HAM", "C_MOM0", "C_MOM1",
                                             "C_MOM2"};

static const unsigned int DIR_X           = 0;
static const unsigned int DIR_Y           = 1;
static const unsigned int DIR_Z           = 2;

/**
 * @brief: Read the parameter file and initialize the variables in parameters.h
 * file.
 * @param[in] fName: file name
 * @param[in] comm: MPI communicator.
 * */
void readParamFile(const char* fName, MPI_Comm comm);

/**
 * @brief Initialize all the variables for a given point in space.
 * @param [in] coord: coordinates of the point.
 * @param [out] var: pointer to the list of variables, computed. var size should
 *be (VAR::U_SYMAT5+1)
 * @note This function is taken from the old single core fluid version.
 **/
void initData(const DendroScalar xx1, const DendroScalar yy1,
              const DendroScalar zz1, DendroScalar* var);

/**
 * @brief: Generates block adaptive octree for the given binary blackhole
 * problem.
 * @param[in/out] tmpNodes: ot::TreeNodes for the block adaptive grids.
 * @param[in] pt_min: min point bound for the regular block.
 * @param[in] pt_max: max point bound for the regular block.
 * @param[in] regLev: regular grid level.
 * @param[in] maxDepth: max dep of the octree
 * @param[in] comm: MPI communicator
 *
 * */
void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,
                         const Point& pt_min, const Point& pt_max,
                         const unsigned int regLev, const unsigned int maxDepth,
                         MPI_Comm comm);

/**
 * @brief wavelet tolerance as a function of space.
 * */
double computeWTol(double x, double y, double z, double tol_min);

/**
 * @brief Compute the wavelet coefficient for a specified element (unzip) vector
 * @tparam T : data type
 * @param pMesh : pointer to the mesh object
 * @param wRefEl : pointer to the wavelet reference element
 * @param wavelet_tol : wavelet tolerance function
 * @param dof : degrees of freedom to consider computing wcoefficient
 * @return double
 */
template <typename T>
double computeFluidWavelet(const ot::Mesh* pMesh,
                           const wavelet::WaveletEl* wRefEl, const T* eleUnzip,
                           double* maxima, double wavelet_tol,
                           unsigned int dof = 1, bool isBdyEle = false);

/**
 * @brief: compute the wavelet refinement flags for the fluid.
 *
 * @tparam T vector data type
 * @param pMesh: pointer to the mesh data structure
 * @param refine_flags: computed refine flags (pass to set_flags in mesh class
 * later. )
 * @param unzippedVec : unzipped vector
 * @param varIds : variable id list to consider during the remesh.
 * @param numVars : number of variables to consider
 * @param wavelet_tol : wavelet tolerance function
 * @param amr_coarse_fac : coarsening factor, coarsen if wCoefficient <
 * amr_coarse_fac*wavelet_tol(ele)
 *
 * @return true when mesh needs to be updated.
 * @return false
 */
template <typename T>
bool computeFluidRemeshFlags(
    const ot::Mesh* pMesh, std::vector<unsigned int>& refine_flags,
    const T** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double)> wavelet_tol,
    double amr_coarse_fac = DENDRO_AMR_COARSEN_FAC, bool includeBdy = true);

}  // end of namespace fluid

namespace fluid {

namespace timer {

/**@brief initialize all the flop counters. */
void initFlops();

/**@brief clears the snapshot counter for time profiler variables*/
void resetSnapshot();

/**@brief reduce min mean max.
 * @param [in] stat: local time
 * @param [out] stat_g 0-min, 1-mean 2-max
 * */
template <typename T>
void computeOverallStats(T* stat, T* stat_g, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    par::Mpi_Reduce(stat, stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 2, 1, MPI_MAX, 0, comm);
    stat_g[1] /= (npes);
}

/** @breif : printout the profile parameters. */
void profileInfo(const char* filePrefix, const ot::Mesh* pMesh);

/** @breif : printout the profile parameters (intermediate profile information).
 */
void profileInfoIntermediate(const char* filePrefix, const ot::Mesh* pMesh,
                             const unsigned int currentStep);

}  // namespace timer

}  // namespace fluid

namespace util {
// Assumes a plane in z, a line in y, and a point in x.
inline unsigned int index_xm(unsigned int i, unsigned int j, unsigned int k,
                             unsigned int np, unsigned int nl,
                             unsigned int nx) {
    return i + nx * (j + k * nl);
}

// Assumes a plane in x, a line in z, and a point in y.
inline unsigned int index_ym(unsigned int i, unsigned int j, unsigned int k,
                             unsigned int np, unsigned int nl,
                             unsigned int nx) {
    return k + np * (i + j * nx);
}

// Assumes a plane in y, a line in x, and a point in z.
inline unsigned int index_zm(unsigned int i, unsigned int j, unsigned int k,
                             unsigned int np, unsigned int nl,
                             unsigned int nx) {
    return j + nl * (k + i * np);
}

inline void cal_pos_xm(unsigned int i, unsigned int j, unsigned int k,
                       double dx, double dl, double dp, double sx, double sl,
                       double sp, double* pos) {
    pos[0] = sx + (i - 0.5) * dx;
    pos[1] = sl + j * dl;
    pos[2] = sp + k * dp;
}

inline void cal_pos_ym(unsigned int i, unsigned int j, unsigned int k,
                       double dx, double dl, double dp, double sx, double sl,
                       double sp, double* pos) {
    pos[1] = sx + (i - 0.5) * dx;
    pos[2] = sl + j * dl;
    pos[0] = sp + k * dp;
}

inline void cal_pos_zm(unsigned int i, unsigned int j, unsigned int k,
                       double dx, double dl, double dp, double sx, double sl,
                       double sp, double* pos) {
    pos[2] = sx + (i - 0.5) * dx;
    pos[0] = sl + j * dl;
    pos[1] = sp + k * dp;
}
}  // namespace util

#include "fluidUtils.tcc"

#endif  // DENDRO_5_FLUID_FLUIDUTILS_H
