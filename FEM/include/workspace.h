//
// Created by milinda on 12/29/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief some workspcae memory allocation to compute element matvec
*/
//

#ifndef SFCSORTBENCH_WORKSPACE_H
#define SFCSORTBENCH_WORKSPACE_H

#include "point.h"

#define JACOBIAN_3D_SIZE 9
//
//        [0,1,2]
// J3d=   [3,4,5]
//        [6,7,8]
//

namespace fem
{
    namespace domain
    {

        // octree grid points min max
        extern Point grid_min;
        extern Point grid_max;

        //domain min max
        extern Point domain_min;
        extern Point domain_max;

        extern double Rg_x;
        extern double Rg_y;
        extern double Rg_z;

        extern double Rd_x;
        extern double Rd_y;
        extern double Rd_z;


        inline double gridX_to_X(double xg){return (((Rd_x/Rg_x)*(xg-fem::domain::grid_min.x()))+fem::domain::domain_min.x());};
        inline double gridY_to_Y(double yg){return (((Rd_y/Rg_y)*(yg-fem::domain::grid_min.y()))+fem::domain::domain_min.y());};
        inline double gridZ_to_Z(double zg){return (((Rd_z/Rg_z)*(zg-fem::domain::grid_min.z()))+fem::domain::domain_min.z());};

        inline double X_to_gridX(double xc){return (((Rg_x/Rd_x)*(xc-fem::domain::domain_min.x()))+fem::domain::grid_min.x());};
        inline double Y_to_gridY(double yc){return (((Rg_y/Rd_y)*(yc-fem::domain::domain_min.y()))+fem::domain::grid_min.y());};
        inline double Z_to_gridZ(double zc){return (((Rg_z/Rd_z)*(zc-fem::domain::domain_min.z()))+fem::domain::grid_min.z());};


    } // end of namespace domain.
} // end of namespace fem


namespace fem
{

    namespace operators
    {
        namespace poisson
        {

            extern double *Qx;
            extern double *Qy;
            extern double *Qz;

            extern double *imV1;
            extern double *imV2;
            extern double *imV3;

            extern double * pts_x;
            extern double * pts_y;
            extern double * pts_z;

            extern double ** jdXbydR;
            extern double ** jdRbydX;
            extern double *  jfactor;



            /**
             * @brief: This function to allocate workspace vars (Note that this will allocate all the memory for extern vars).
             *
             * */
            void allocateWorkSpace(const unsigned int n);

            /**
             * @brief: deallocates the memory of extern vars.
             * */
            void deallocateWorkSpace();


        } // end of nampspace possoin

    }// end of namespace operators

}// end of namespace fem




#endif //SFCSORTBENCH_WORKSPACE_H
