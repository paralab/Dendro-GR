//
// Created by milinda on 8/9/18.
//
/**
 * @brief Contians cuda derivs for bssn computation
 *
 * */
#ifndef SFCSORTBENCH_DERIVS_H
#define SFCSORTBENCH_DERIVS_H



#include "cuda_runtime.h"
#include <stdio.h>
#include <iostream>

// shared mem. index
#define IDX_L(i,j,k) ( (i) + tile_sz[0] * ( (j) + tile_sz[1] * (k) ) )
// gloabal index.
#define IDX_G(i,j,k) ( (i) + sz[0] * ( (j) + sz[1] * (k) ) )

namespace cuda
{


        // 1D stencil operations.
        __device__ void _RSWS_deriv42_x_1D(double * const  Dxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);
        __device__ void _RSWS_deriv42_y_1D(double * const  Dyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);
        __device__ void _RSWS_deriv42_z_1D(double * const  Dzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);

        // _RSWS_ derivates for shared mem read shared mem write in GPUs.
        // first derivatives
        __device__ void _RSWS_deriv42_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);
        __device__ void _RSWS_deriv42_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);
        __device__ void _RSWS_deriv42_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw ,unsigned bflag);

        // advective derivatives
        __device__ void _RSWS_deriv42adv_z(double * const  Dzu, const double * const  u,const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const bool * const betaz, unsigned int pw,unsigned bflag);
        __device__ void _RSWS_deriv42adv_y(double * const  Dyu, const double * const  u,const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const bool * const betay, unsigned int pw,unsigned bflag);
        __device__ void _RSWS_deriv42adv_x(double * const  Dxu, const double * const  u,const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, const bool * const betax, unsigned int pw,unsigned bflag);

        // second order derivatives
        __device__ void _RSWS_deriv42_zz(double * const  Du, const double * const  u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw ,unsigned bflag);
        __device__ void _RSWS_deriv42_yy(double * const  Du, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw , unsigned bflag);
        __device__ void _RSWS_deriv42_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw, unsigned bflag);


        // Kriess-Oliger derivatives
        __device__ void _RSWS_ko_deriv42_z(double * const Du, const double * const u, const double dz, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz, unsigned int pw,unsigned bflag);
        __device__ void _RSWS_ko_deriv42_y(double * const  Du, const double * const  u, const double dy, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag);
        __device__ void _RSWS_ko_deriv42_x(double * const  Du, const double * const  u, const double dx, const unsigned int* ijk_lm, const unsigned int * sz, const unsigned int* tile_sz,unsigned int pw, unsigned bflag);


}// end of namespace cuda
#endif //SFCSORTBENCH_DERIVS_H


