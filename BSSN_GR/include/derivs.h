#ifndef _DERVS_H
#define _DERVS_H

#include <cmath>

#include "TreeNode.h"

#define IDX(i, j, k) ((i) + nx * ((j) + ny * (k)))

// stores the padding width.
extern unsigned int DERIV_PW;

// avx vector length
#define __DERIV_AVX_SIMD_LEN__ 16
#define __RHS_AVX_SIMD_LEN__   16
#define __DEFAULT_AVX_SIMD_LEN 16

void deriv42_x(double *const Dxu, const double *const u, const double dx,
               const unsigned int *sz, unsigned bflag);
void deriv42_y(double *const Dyu, const double *const u, const double dy,
               const unsigned int *sz, unsigned bflag);
void deriv42_z(double *const Dzu, const double *const u, const double dz,
               const unsigned int *sz, unsigned bflag);

void deriv42adv_x(double *const Dxu, const double *const u, const double dx,
                  const unsigned int *sz, const double *const betax,
                  unsigned bflag);
void deriv42adv_y(double *const Dyu, const double *const u, const double dy,
                  const unsigned int *sz, const double *const betay,
                  unsigned bflag);
void deriv42adv_z(double *const Dzu, const double *const u, const double dz,
                  const unsigned int *sz, const double *const betaz,
                  unsigned bflag);

void deriv42_xx(double *const DxDxu, const double *const u, const double dx,
                const unsigned int *sz, unsigned bflag);
void deriv42_yy(double *const Du, const double *const u, const double dy,
                const unsigned int *sz, unsigned bflag);
void deriv42_zz(double *const Du, const double *const u, const double dz,
                const unsigned int *sz, unsigned bflag);

void ko_deriv42_x(double *const Du, const double *const u, const double dx,
                  const unsigned int *sz, unsigned bflag);
void ko_deriv42_y(double *const Du, const double *const u, const double dy,
                  const unsigned int *sz, unsigned bflag);
void ko_deriv42_z(double *const Du, const double *const u, const double dz,
                  const unsigned int *sz, unsigned bflag);

void ko_pw4_deriv42_x(double *const Du, const double *const u, const double dx,
                      const unsigned int *sz, unsigned bflag);
void ko_pw4_deriv42_y(double *const Du, const double *const u, const double dy,
                      const unsigned int *sz, unsigned bflag);
void ko_pw4_deriv42_z(double *const Du, const double *const u, const double dz,
                      const unsigned int *sz, unsigned bflag);

void disstvb3_x(double *const Du, const double *const u,
                const double *const lam, const double dx,
                const unsigned int *sz, unsigned bflag);
void disstvb3_y(double *const Du, const double *const u,
                const double *const lam, const double dx,
                const unsigned int *sz, unsigned bflag);
void disstvb3_z(double *const Du, const double *const u,
                const double *const lam, const double dx,
                const unsigned int *sz, unsigned bflag);

void disstvb5_x(double *const Du, const double *const u,
                const double *const lam, const double dx,
                const unsigned int *sz, unsigned bflag);
void disstvb5_y(double *const Du, const double *const u,
                const double *const lam, const double dx,
                const unsigned int *sz, unsigned bflag);
void disstvb5_z(double *const Du, const double *const u,
                const double *const lam, const double dx,
                const unsigned int *sz, unsigned bflag);

// 6th order derivatives.
void deriv644_x(double *const Dxu, const double *const u, const double dx,
                const unsigned int *sz, unsigned bflag);
void deriv644_y(double *const Dyu, const double *const u, const double dy,
                const unsigned int *sz, unsigned bflag);
void deriv644_z(double *const Dzu, const double *const u, const double dz,
                const unsigned int *sz, unsigned bflag);

void deriv644_xx(double *const DxDxu, const double *const u, const double dx,
                 const unsigned int *sz, unsigned bflag);
void deriv644_yy(double *const Du, const double *const u, const double dy,
                 const unsigned int *sz, unsigned bflag);
void deriv644_zz(double *const Du, const double *const u, const double dz,
                 const unsigned int *sz, unsigned bflag);

void deriv644adv_x(double *const Dxu, const double *const u, const double dx,
                   const unsigned int *sz, const double *const betax,
                   unsigned bflag);
void deriv644adv_y(double *const Dyu, const double *const u, const double dy,
                   const unsigned int *sz, const double *const betay,
                   unsigned bflag);
void deriv644adv_z(double *const Dzu, const double *const u, const double dz,
                   const unsigned int *sz, const double *const betaz,
                   unsigned bflag);

void deriv642_x(double *const Dxu, const double *const u, const double dx,
                const unsigned int *sz, unsigned bflag);
void deriv642_y(double *const Dyu, const double *const u, const double dy,
                const unsigned int *sz, unsigned bflag);
void deriv642_z(double *const Dzu, const double *const u, const double dz,
                const unsigned int *sz, unsigned bflag);

void deriv642_xx(double *const DxDxu, const double *const u, const double dx,
                 const unsigned int *sz, unsigned bflag);
void deriv642_yy(double *const Du, const double *const u, const double dy,
                 const unsigned int *sz, unsigned bflag);
void deriv642_zz(double *const Du, const double *const u, const double dz,
                 const unsigned int *sz, unsigned bflag);

void deriv642adv_x(double *const Dxu, const double *const u, const double dx,
                   const unsigned int *sz, const double *const betax,
                   unsigned bflag);
void deriv642adv_y(double *const Dyu, const double *const u, const double dy,
                   const unsigned int *sz, const double *const betay,
                   unsigned bflag);
void deriv642adv_z(double *const Dzu, const double *const u, const double dz,
                   const unsigned int *sz, const double *const betaz,
                   unsigned bflag);

void ko_deriv64_x(double *const Du, const double *const u, const double dx,
                  const unsigned int *sz, unsigned bflag);
void ko_deriv64_y(double *const Du, const double *const u, const double dy,
                  const unsigned int *sz, unsigned bflag);
void ko_deriv64_z(double *const Du, const double *const u, const double dz,
                  const unsigned *sz, unsigned bflag);

// 8th order derivatives.
void deriv8642_x(double *const Dxu, const double *const u, const double dx,
                 const unsigned int *sz, unsigned bflag);
void deriv8642_y(double *const Dyu, const double *const u, const double dy,
                 const unsigned int *sz, unsigned bflag);
void deriv8642_z(double *const Dzu, const double *const u, const double dz,
                 const unsigned int *sz, unsigned bflag);

void deriv8642_xx(double *const DxDxu, const double *const u, const double dx,
                  const unsigned int *sz, unsigned bflag);
void deriv8642_yy(double *const Du, const double *const u, const double dy,
                  const unsigned int *sz, unsigned bflag);
void deriv8642_zz(double *const Du, const double *const u, const double dz,
                  const unsigned int *sz, unsigned bflag);

void deriv8642adv_x(double *const Dxu, const double *const u, const double dx,
                    const unsigned int *sz, const double *const betax,
                    unsigned bflag);
void deriv8642adv_y(double *const Dyu, const double *const u, const double dy,
                    const unsigned int *sz, const double *const betay,
                    unsigned bflag);
void deriv8642adv_z(double *const Dzu, const double *const u, const double dz,
                    const unsigned int *sz, const double *const betaz,
                    unsigned bflag);

void deriv8664_x(double *const Dxu, const double *const u, const double dx,
                 const unsigned int *sz, unsigned bflag);
void deriv8664_y(double *const Dyu, const double *const u, const double dy,
                 const unsigned int *sz, unsigned bflag);
void deriv8664_z(double *const Dzu, const double *const u, const double dz,
                 const unsigned int *sz, unsigned bflag);

void deriv8664_xx(double *const DxDxu, const double *const u, const double dx,
                  const unsigned int *sz, unsigned bflag);
void deriv8664_yy(double *const Du, const double *const u, const double dy,
                  const unsigned int *sz, unsigned bflag);
void deriv8664_zz(double *const Du, const double *const u, const double dz,
                  const unsigned int *sz, unsigned bflag);

void deriv8644adv_x(double *const Dxu, const double *const u, const double dx,
                    const unsigned int *sz, const double *const betax,
                    unsigned bflag);
void deriv8644adv_y(double *const Dyu, const double *const u, const double dy,
                    const unsigned int *sz, const double *const betay,
                    unsigned bflag);
void deriv8644adv_z(double *const Dzu, const double *const u, const double dz,
                    const unsigned int *sz, const double *const betaz,
                    unsigned bflag);

void deriv8666_x(double *const Dxu, const double *const u, const double dx,
                 const unsigned int *sz, unsigned bflag);
void deriv8666_y(double *const Dyu, const double *const u, const double dy,
                 const unsigned int *sz, unsigned bflag);
void deriv8666_z(double *const Dzu, const double *const u, const double dz,
                 const unsigned int *sz, unsigned bflag);

void deriv8666_xx(double *const DxDxu, const double *const u, const double dx,
                  const unsigned int *sz, unsigned bflag);
void deriv8666_yy(double *const Du, const double *const u, const double dy,
                  const unsigned int *sz, unsigned bflag);
void deriv8666_zz(double *const Du, const double *const u, const double dz,
                  const unsigned int *sz, unsigned bflag);

void deriv8666adv_x(double *const Dxu, const double *const u, const double dx,
                    const unsigned int *sz, const double *const betax,
                    unsigned bflag);
void deriv8666adv_y(double *const Dyu, const double *const u, const double dy,
                    const unsigned int *sz, const double *const betay,
                    unsigned bflag);
void deriv8666adv_z(double *const Dzu, const double *const u, const double dz,
                    const unsigned int *sz, const double *const betaz,
                    unsigned bflag);

// these are the derivs that will be used in the BSSN CODE based on the FD
// derivative order.
#ifdef BSSN_USE_4TH_ORDER_DERIVS
#define deriv_x deriv42_x
#define deriv_y deriv42_y
#define deriv_z deriv42_z

#define deriv_xx deriv42_xx
#define deriv_yy deriv42_yy
#define deriv_zz deriv42_zz

#define adv_deriv_x deriv42adv_x
#define adv_deriv_y deriv42adv_y
#define adv_deriv_z deriv42adv_z

#define ko_deriv_x ko_deriv42_x
#define ko_deriv_y ko_deriv42_y
#define ko_deriv_z ko_deriv42_z
#endif

#ifdef BSSN_USE_6TH_ORDER_DERIVS
#define deriv_x deriv644_x
#define deriv_y deriv644_y
#define deriv_z deriv644_z

#define deriv_xx deriv644_xx
#define deriv_yy deriv644_yy
#define deriv_zz deriv644_zz

#define adv_deriv_x deriv644adv_x
#define adv_deriv_y deriv644adv_y
#define adv_deriv_z deriv644adv_z

#define ko_deriv_x ko_deriv42_x
#define ko_deriv_y ko_deriv42_y
#define ko_deriv_z ko_deriv42_z
#endif

#ifdef BSSN_USE_8TH_ORDER_DERIVS
#define deriv_x deriv8642_x
#define deriv_y deriv8642_y
#define deriv_z deriv8642_z

#define deriv_xx deriv8642_xx
#define deriv_yy deriv8642_yy
#define deriv_zz deriv8642_zz

#define ko_deriv_x ko_pw4_deriv42_x
#define ko_deriv_y ko_pw4_deriv42_y
#define ko_deriv_z ko_pw4_deriv42_z
#endif

#endif
