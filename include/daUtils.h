//
// Created by milinda on 4/27/18.
//
/**
 * @author Milinda Fernando
 * School of Computing, University of Utah.
 * @brief contains utilitity functions of mesh based distributed array(da)
 * */

#ifndef SFCSORTBENCH_DAUTILS_H
#define SFCSORTBENCH_DAUTILS_H

#include "mesh.h"
#include <iostream>


namespace ot
{
    namespace da
    {

        /**
         * @brief interpolates input vector and a local element to the given coordinate which lies in the specified octant
         * @param[in] mesh : input mesh
         * @param[in] in : input vector (zipped) version defined on the mesh.
         * @param[in] coord: coordinate (size 3 vector)
         * @param[in] elementID: element ID (note that the coord should lie in the elementID)
         * @param[in] interpOrder: Lagrange interpolation order.
         * @return interpolated value based on Lagrange polynomials.
         *
         * */
        template <typename T>
        T lagrangeInterpElementToCoord(const ot::Mesh * mesh,const T* in,double * coord,unsigned int elementID,unsigned int interpOrder);


        /**
         * @brief interpolates a given input vector to given coordinate values.
         * @param[in] mesh : input mesh
         * @param[in] in : input vector (zipped) version defined on the mesh.
         * @param[in] coords: pointer to the list of coordinates (specified in the global coords, coords that are not in the local parition will be ignored)
         * @param[in] length: input length of the coordinates
         * @param[out] out: interpolated values.
         * @param[out] out_size: output size (based on how many coords in local partition of the octree)
         *
         * */
        template<typename T, typename CoordT>
        void interpolateToCoords(const ot::Mesh * mesh, const T* in, const CoordT* coords,unsigned int length,T* out,std::vector<unsigned int >& validIndices);

    } // end of namespace da
}// end of namespace ot


#include "daUtils.tcc"

#endif //SFCSORTBENCH_DAUTILS_H
