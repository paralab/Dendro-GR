//
// Created by milinda on 12/20/18.
//

/**
 * @brief contains ODA utility functions,
 * @author Milinda Fernando, School of Computing, University of Utah
 * */

#ifndef DENDRO_5_0_ODAUTILS_H
#define DENDRO_5_0_ODAUTILS_H

#include "binUtils.h"
#include <iostream>
#include <vector>
#include "mesh.h"
#include "parUtils.h"



namespace ot
{
    /**
     * @brief flag the elements for the ODA construction
     * @param [in] pMesh : ot::Mesh
     * @param [out] flagList: flaged output
     */
    void computeODAFlags(const ot::Mesh* pMesh, std::vector<unsigned int>& flagList);

     /**
      * @brief computes the local to global nodal map.
      * @param[in] pMesh: ot::Mesh generated form 2:1 balanced octree
      * @param[out] map: constructed local to global nodal map.
      * */
    void computeLocalToGlobalNodalMap(const ot::Mesh* pMesh,std::vector<DendroIntL>& map, DendroIntL& globalNodeSz, std::vector<DendroIntL>& nodalOffset);


} // end of namespace ot

#endif //DENDRO_5_0_ODAUTILS_H
