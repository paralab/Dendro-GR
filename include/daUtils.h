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
         * @param[in] domain_coord: coordinate in domain basis
         * @param[in] pt_min: min element coord in domain coordinates
         * @param[in] pt_max: max element coord in domain coordinates
         * @param[in] elementID: element ID (note that the coord should lie in the elementID)
         * @param[in] interpOrder: Lagrange interpolation order.
         * @return interpolated value based on Lagrange polynomials.
         * */
        template <typename T>
        T lagrangeInterpElementToCoord(const ot::Mesh * mesh,const T* in, double * domain_coord, const Point& pt_min, const Point&  pt_max , unsigned int elementID,unsigned int interpOrder);


        /**
         * @brief Lagrange linear interp, function computes the linear order simplex containting the pt, and then use the linear lagrange interp on the element to interpolate.
         * @param[in] mesh : input mesh
         * @param[in] in : input vector (zipped) version defined on the mesh.
         * @param[in] domain_coord: coordinate in domain basis
         * @param[in] pt_min: min element coord in domain coordinates
         * @param[in] pt_max: max element coord in domain coordinates
         * @param[in] elementID: element ID (note that the coord should lie in the elementID)
         * @param[in] interpOrder: Lagrange interpolation order.
         * @return interpolated value based on Lagrange polynomials.
         * */

        template <typename T>
        T linear_lagrange(const ot::Mesh * mesh,const T* in, double * domain_coord, const Point& pt_min, const Point&  pt_max , unsigned int elementID);


        /**
         * @brief interpolates a given input vector to given coordinate values.
         * @assumption: This routine assumes that, it transforms the octree cube domain to the physical domain of cube.  
         * @param[in] mesh : input mesh
         * @param[in] in : input vector (zipped) version defined on the mesh.
         * @param[in] domain_coords: pointer to the list of coordinates in the domain reference.  (specified in the global coords, coords that are not in the local parition will be ignored)
         * @param[in] length: input length of the coordinates
         * @param[in] grid_limit : size 2 array of points for min and max point on the grid.
         * @param[in] domain_limit : size 2 array of points for min and max point on the domain.
         * @param[out] out: interpolated values.
         * @param[out] valid_index: indices of out that are found in the current rank
         * */
        template<typename T, typename CoordT>
        void interpolateToCoords(const ot::Mesh * mesh, const T* in, const CoordT* domain_coords, unsigned int length, const Point* const grid_limit, const Point* const domain_limit, T* out,std::vector<unsigned int >& validIndices);

        /**
         * @brief interpolates a given input vector to given coordinate values.
         * @assumption: This routine assumes that, it transforms the octree cube domain to the physical domain of cube.  
         * @param[in] mesh : input mesh
         * @param[in] in : input vector (zipped) version defined on the mesh.
         * @param[in] domain_coords: pointer to the list of coordinates in the domain reference.  (specified in the global coords, coords that are not in the local parition will be ignored)
         * @param[in] length: input length of the coordinates
         * @param[in] grid_limit : size 2 array of points for min and max point on the grid.
         * @param[in] domain_limit : size 2 array of points for min and max point on the domain.
         * @param[out] out: interpolated values (only need to allocated in the root).
         * @param[in] root: root rank to gather the nodes. 
         * */
        template<typename T, typename CoordT>
        void interpolateToCoordsAndGather(const ot::Mesh * mesh, const T* in, const CoordT* domain_coords, unsigned int length, const Point* const grid_limit, const Point* const domain_limit, T* out,unsigned int root,unsigned int dof);

    } // end of namespace da
}// end of namespace ot


#include "daUtils.tcc"

#endif //SFCSORTBENCH_DAUTILS_H
