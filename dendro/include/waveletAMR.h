/**
 * @file waveletAMR.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Contains utility functions to perform wavelet AMR. 
 * @version 0.1
 * @date 2020-06-13
 * 
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include "mesh.h"
#include "waveletRefEl.h"
namespace wavelet
{   

    /**
     * @brief Computes the wavelet coefficient for a specified element (unzip) vector
     * @tparam T : data type.
     * @param pMesh : pointer to the mesh object. 
     * @param wRefEl : pointer to the wavelet reference element. 
     * @param eleUnzip : element unzip vector()
     * @param wavelet_tol : wavelet tolerance function
     * @param dof :degrees of freedoms to consider computing wcoefficient. 
     * @return double 
     */
    template<typename T>
    double compute_element_wavelet(const ot::Mesh* pMesh, const WaveletEl* wRefEl, const T*eleUnzip, double wavelet_tol, unsigned int dof=1,bool isBdyEle=false);
    

    /**
     * @brief Compute the remesh flags of the entire mesh looping over local partition. 
     * 
     * @tparam T vector data type
     * @param pMesh : pointer to the mesh data strucutre
     * @param refine_flags : computed refine flags (pass to set_flags in mesh class later. )
     * @param unzippedVec : unzipped vector
     * @param varIds : variable id list to consider during the remesh. 
     * @param numVars : number of variables to consider. 
     * @param wavelet_tol : wavelet tolerance function
     * @param amr_coarse_fac : coarsening factor, coarsen if wCoefficient < amr_coarse_fac*wavelet_tol(ele)
     * @return true when mesh needs to updated. 
     * @return false 
     */
    template<typename T>
    bool compute_wavelet_remesh_flags(const ot::Mesh* pMesh, std::vector<unsigned int>& refine_flags, const T**unzippedVec, const unsigned int *varIds, const unsigned int numVars, std::function<double(double, double, double)> wavelet_tol, double amr_coarse_fac = DENDRO_AMR_COARSEN_FAC, bool includeBdy=false);



    
    

}// end of namespace wavelet. 

#include "waveletAMR.tcc"