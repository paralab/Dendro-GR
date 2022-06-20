#pragma once
#include <cmath>
#include "TreeNode.h"

#define IDX(i, j, k) ((i) + nx * ((j) + ny * (k)))

// avx vector length
#define __DERIV_AVX_SIMD_LEN__ 16
#define __RHS_AVX_SIMD_LEN__ 16
#define __DEFAULT_AVX_SIMD_LEN 16


void deriv22_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv22_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv22_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv42_x(double *const Dxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv42_y(double *const Dyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_z(double *const Dzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv42_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv42_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv42adv_x(double * const  Dxu, const double * const  u,const double dx, const unsigned int *sz, const double * const betax, unsigned bflag);
void deriv42adv_y(double * const  Dyu, const double * const  u,const double dy, const unsigned int *sz, const double * const betay, unsigned bflag);
void deriv42adv_z(double * const  Dzu, const double * const  u,const double dz, const unsigned int *sz, const double * const betaz, unsigned bflag);

void deriv42_xy(double *const DxDyu, double* const Dxu, const double *const u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv42_xz(double *const DxDzu, double* const Dxu, const double *const u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv42_yz(double *const DyDzu, double* const Dyu, const double *const u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);

void ko_deriv42_x(double *const Du, const double *const u, const double dx, const unsigned int *sz, unsigned bflag);
void ko_deriv42_y(double *const Du, const double *const u, const double dy, const unsigned int *sz, unsigned bflag);
void ko_deriv42_z(double *const Du, const double *const u, const double dz, const unsigned int *sz, unsigned bflag);


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


void deriv642_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv642_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv642_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv642_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int *sz, unsigned bflag);
void deriv642_yy(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv642_zz(double * const  Du, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv642_xy(double * const  DxDyu, double * const Dxu, const double * const  u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv642_xz(double * const  DxDzu, double * const Dxu, const double * const  u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv642_yz(double * const  DyDzu, double * const Dyu, const double * const  u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);

void ko_deriv64_x(double * const  Du, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);
void ko_deriv64_y(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void ko_deriv64_z(double * const Du, const double * const u, const double dz, const unsigned *sz, unsigned bflag);


// 8th order derivatives. 
void deriv8642_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8642_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8642_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8642_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int *sz, unsigned bflag);
void deriv8642_yy(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8642_zz(double * const  Du, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8642_xy(double * const  DxDyu, double * const Dxu, const double * const  u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8642_xz(double * const  DxDzu, double * const Dxu, const double * const  u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv8642_yz(double * const  DyDzu, double * const Dyu, const double * const  u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8664_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8664_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8664_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8664_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int *sz, unsigned bflag);
void deriv8664_yy(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8664_zz(double * const  Du, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8664_xy(double * const  DxDyu, double * const Dxu, const double * const  u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8664_xz(double * const  DxDzu, double * const Dxu, const double * const  u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv8664_yz(double * const  DyDzu, double * const Dyu, const double * const  u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8666_x(double * const  Dxu, const double * const  u, const double dx, const unsigned int *sz, unsigned bflag);
void deriv8666_y(double * const  Dyu, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8666_z(double * const  Dzu, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8666_xx(double * const  DxDxu, const double * const  u,const double dx, const unsigned int *sz, unsigned bflag);
void deriv8666_yy(double * const  Du, const double * const  u, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8666_zz(double * const  Du, const double * const  u, const double dz, const unsigned int *sz, unsigned bflag);

void deriv8666_xy(double * const  DxDyu, double * const Dxu, const double * const  u, const double dx, const double dy, const unsigned int *sz, unsigned bflag);
void deriv8666_xz(double * const  DxDzu, double * const Dxu, const double * const  u, const double dx, const double dz, const unsigned int *sz, unsigned bflag);
void deriv8666_yz(double * const  DyDzu, double * const Dyu, const double * const  u, const double dy, const double dz, const unsigned int *sz, unsigned bflag);


#ifdef NLSM_USE_4TH_ORDER_DERIVS
    #define deriv_x deriv42_x
    #define deriv_y deriv42_y
    #define deriv_z deriv42_z

    #define deriv_xx deriv42_xx
    #define deriv_yy deriv42_yy
    #define deriv_zz deriv42_zz

    #define ko_deriv_x ko_deriv42_x
    #define ko_deriv_y ko_deriv42_y
    #define ko_deriv_z ko_deriv42_z
#endif


#ifdef NLSM_USE_6TH_ORDER_DERIVS
    #define deriv_x deriv644_x
    #define deriv_y deriv644_y
    #define deriv_z deriv644_z

    #define deriv_xx deriv644_xx
    #define deriv_yy deriv644_yy
    #define deriv_zz deriv644_zz

    #define ko_deriv_x ko_deriv42_x
    #define ko_deriv_y ko_deriv42_y
    #define ko_deriv_z ko_deriv42_z
#endif

#ifdef NLSM_USE_8TH_ORDER_DERIVS
    #define deriv_x deriv8642_x
    #define deriv_y deriv8642_y
    #define deriv_z deriv8642_z

    #define deriv_xx deriv8642_xx
    #define deriv_yy deriv8642_yy
    #define deriv_zz deriv8642_zz

    #define ko_deriv_x ko_deriv64_x
    #define ko_deriv_y ko_deriv64_y
    #define ko_deriv_z ko_deriv64_z
#endif
