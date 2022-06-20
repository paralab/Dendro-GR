#pragma once
#include <cmath>
#include "TreeNode.h"
#include "parameters.h"

#define IDX(i, j, k) ((i) + nx * ((j) + ny * (k)))


void deriv42_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv42_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv42_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv42_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

/**
 * @brief Note that we are computing the mixed derivatives, after applying the stencils composition
 * @param DxDyu : mixed deriv. 
 * @param Dx : Dx(u)
 * @param u : input vector
 * @param dx : dx, 
 * @param dy : dy
 * @param sz : size of the u in 3D
 * @param bflag : boundary flags
 */

void deriv42_xy(double *const DxDyu, double* const Dxu, const double *const u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_xz(double *const DxDzu, double* const Dxu, const double *const u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv42_yz(double *const DyDzu, double* const Dyu, const double *const u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);

void ko_deriv42_x(double *const Du, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void ko_deriv42_y(double *const Du, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void ko_deriv42_z(double *const Du, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv22_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv22_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv22_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);


// -- (centered) 6th order derivs (in the interior) with (shifted) 4th order derivs near and on the boundaries. 
void deriv644_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv644_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv644_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv644_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv644_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv644_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv644_xy(double *const DxDyu, double* const Dxu, const double *const u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv644_xz(double *const DxDzu, double* const Dxu, const double *const u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv644_yz(double *const DyDzu, double* const Dyu, const double *const u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);


// -- (centered) 6th order derivs (in the interior) reducing to centered 4th order and centered 2nd order derivs near boundaries. On the boundary, we use a shifted second order derivative approximation.  
void deriv642_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv642_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv642_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv642_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv642_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv642_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv642_xy(double *const DxDyu, double* const Dxu, const double *const u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv642_xz(double *const DxDzu, double* const Dxu, const double *const u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv642_yz(double *const DyDzu, double* const Dyu, const double *const u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);

void ko_deriv64_x(double *const Du, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void ko_deriv64_y(double *const Du, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void ko_deriv64_z(double *const Du, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);


// -- (centered) 8th order derivs (in the interior) reducing to centered 6th order, centered 4th order and centered 2nd order derivs near boundaries. On the boundary, we use a shifted second order derivative approximation.  
void deriv8642_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8642_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8642_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8642_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8642_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8642_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8642_xy(double *const DxDyu, double* const Dxu, const double *const u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8642_xz(double *const DxDzu, double* const Dxu, const double *const u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv8642_yz(double *const DyDzu, double* const Dyu, const double *const u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);


// -- (centered) 8th order derivs (in the interior) reducing to centered 6th order, shifted 6th order and shifted 4th order derivs near boundaries. On the boundary, we use a shifted 4th order derivative approximation.  
void deriv8664_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8664_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8664_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8664_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8664_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8664_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8664_xy(double *const DxDyu, double* const Dxu, const double *const u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8664_xz(double *const DxDzu, double* const Dxu, const double *const u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv8664_yz(double *const DyDzu, double* const Dyu, const double *const u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);


// -- (centered) 8th order derivs (in the interior) reducing to centered 6th order and shifted 6th order near and on boundaries. 
void deriv8666_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8666_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8666_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8666_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8666_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8666_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8666_xy(double *const DxDyu, double* const Dxu, const double *const u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8666_xz(double *const DxDzu, double* const Dxu, const double *const u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv8666_yz(double *const DyDzu, double* const Dyu, const double *const u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);

#ifdef EM4_USE_4TH_ORDER_DERIVS
    #define deriv_x deriv42_x
    #define deriv_y deriv42_y
    #define deriv_z deriv42_z

    #define deriv_xx deriv42_xx
    #define deriv_yy deriv42_yy
    #define deriv_zz deriv42_zz

    #define deriv_xy deriv42_xy
    #define deriv_xz deriv42_xz
    #define deriv_yz deriv42_yz

    // #define adv_deriv_x deriv42adv_x
    // #define adv_deriv_y deriv42adv_y
    // #define adv_deriv_z deriv42adv_z

    #define ko_deriv_x ko_deriv42_x
    #define ko_deriv_y ko_deriv42_y
    #define ko_deriv_z ko_deriv42_z
#endif


#ifdef EM4_USE_6TH_ORDER_DERIVS
    #define deriv_x deriv644_x
    #define deriv_y deriv644_y
    #define deriv_z deriv644_z

    #define deriv_xx deriv644_xx
    #define deriv_yy deriv644_yy
    #define deriv_zz deriv644_zz

    #define deriv_xy deriv644_xy
    #define deriv_xz deriv644_xz
    #define deriv_yz deriv644_yz

    #define ko_deriv_x ko_deriv42_x
    #define ko_deriv_y ko_deriv42_y
    #define ko_deriv_z ko_deriv42_z
#endif


#ifdef EM4_USE_8TH_ORDER_DERIVS
    #define deriv_x deriv8666_x
    #define deriv_y deriv8666_y
    #define deriv_z deriv8666_z

    #define deriv_xx deriv8666_xx
    #define deriv_yy deriv8666_yy
    #define deriv_zz deriv8666_zz

    #define deriv_xy deriv8666_xy
    #define deriv_xz deriv8666_xz
    #define deriv_yz deriv8666_yz

    #define ko_deriv_x ko_deriv64_x
    #define ko_deriv_y ko_deriv64_y
    #define ko_deriv_z ko_deriv64_z
#endif
