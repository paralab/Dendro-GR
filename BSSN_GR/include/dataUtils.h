//
// Created by milinda on 1/16/19.
//

/**
 * @author Milinda Fernando
 * School of Computing univerity of Utah
 * @brief Contain utility functions to perform processing simulation data and
 * post processing.
 *
 */

#ifndef DENDRO_5_0_DATAUTILS_H
#define DENDRO_5_0_DATAUTILS_H

#include "TreeNode.h"
#include "grDef.h"
#include "mesh.h"
#include "parameters.h"
#include "point.h"

namespace bssn {

/**
 * @brief performs the BH coordinates extraction
 * @param[in] pMesh: current mesh data structure
 * @param[in] var: considering variable for BH extraction
 * @param[in] tolerance: tolerance value for coords extraction
 * @param[in] ptIn: previous known location of the BHs.
 * @param[in] numPt: number of points
 * @param[out] ptOut: new locations of the BHs. (only rank =0 of active comm
 * will have the coords)
 * */
void extractBHCoords(const ot::Mesh* pMesh, const DendroScalar* var,
                     double tolerance, const Point* ptIn, unsigned int numPt,
                     Point* ptOut);

/**
 *@brief write the extracted BH coordinates to a file.
 *@param[in] ptLocs: point locations
 *@param[in] numPt: number of points
 *@param[in] timestep: time step
 *@param[in] time: current time.
 */
void writeBHCoordinates(const ot::Mesh* pMesh, const Point* ptLocs,
                        unsigned int numPt, unsigned int timestep, double time);

/**
 * @brief refine based on the black hole location locations.
 * @param[in] pMesh : pointer to the mesh.
 * @param[in] bhLoc : BH location
 */
bool isRemeshBH(ot::Mesh* pMesh, const Point* bhLoc);

/**
 * @brief refine only based on the alpha variable event horizon.
 *
 * @param pMesh : pointer to the mesh
 * @param unzipVec : unzip vars.
 * @param refine_th : refine tol for alpha
 * @param coarsen_th : coarsend threshold for alpha
 * @return true : is mesh need to be changed
 * @return false : otherwise.
 */
bool isRemeshEH(ot::Mesh* pMesh, const double** unzipVec, unsigned int vIndex,
                double refine_th, double coarsen_th, bool isOverwrite = true);

/**
 * @brief refine only based on the wavelets
 * @param pMesh : pointer to the mesh
 * @param unzippedVec : unzip vars.
 * @param varIds : refinement variable ids.
 * @param numVars : number of varIds
 * @param wavelet_tol: wavelet tolerance function
 * @param amr_coarsen_fac: AMR coarsening safety factor , coarsen iff w_c <
 * amr_coarsen_fac * wavelet_tol(x,y,z)
 */
bool isReMeshWAMR(
    ot::Mesh* pMesh, const double** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double, double*)> wavelet_tol,
    double amr_coarse_fac);



/**
 * @brief add refinement based on wavelets
 * @param pMesh : pointer to the mesh
 * @param unzippedVec : unzip vars.
 * @param varIds : refinement variable ids.
 * @param numVars : number of varIds
 * @param wavelet_tol: wavelet tolerance function
 * @param amr_coarsen_fac: AMR coarsening safety factor, 
 * coarsen iff w_c < amr_coarsen_fac * wavelet_tol(x,y,z)
 * @param relative_WAMR : toggle using relative wavelet values 
 */
bool addRemeshWAMR(
    ot::Mesh* pMesh, const double** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double, double*)> wavelet_tol,
    double amr_coarse_fac, bool relative_WAMR = false);


/**
 * @brief refine ratially based on the BH locations and AMR_R.
 * @param pMesh pointer to the mesh object.
 */
bool isReMeshBHRadial(ot::Mesh* pMesh);

}  // end of namespace bssn



#endif  // DENDRO_5_0_DATAUTILS_H
