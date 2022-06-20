#ifndef FLUID_HLLE_H_
#define FLUID_HLLE_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cassert>

#include "fluid.h"

namespace hlle{
  void rhshlle(double *dtu[], const double *u[], const double *v[], const double *pmin,
               const double *pmax, const unsigned int *sz);
  void cal_nflux(double *dtu[], const double *u[], const double *v[], double *ul1[],
                 double *ur1[], double *vl1[], double *vr1[], double *F1[],
                 double *v1[], const unsigned int dir,
                 const double *pmin, const double *pmax,
                 const unsigned int *sz);
  void cal_source(double *dtu[], const double *u[], const double *v[], const double *pmin,
                  const double *pmax, const unsigned int *sz);
  void lf_flux(double *F[], double *ul1[], double *ur1[],
               double *vl1[], double *vr1[],
               const int dir, const double pos_in[], const double dx,
               const unsigned int *sz);
  void hlle_flux(double *F[], double *ul1[], double *ur1[],
                 double *vl1[], double *vr1[],
                 const int dir, const double pos_in[], const double dx,
                 const unsigned int *sz);
  void flux_x(double f[], double u[], double v[],
              double sdetg, double gd[3][3], double gu[3][3]);
  void flux_y(double f[], double u[], double v[],
              double sdetg, double gd[3][3], double gu[3][3]);
  void flux_z(double f[], double u[], double v[],
              double sdetg, double gd[3][3], double gu[3][3]);
  int cal_max_speeds(double *bp, double *bm, double vl[], double vr[],
                     double gd[3][3], double gu[3][3],
                     const int dir);
  void default_speeds(double *bp, double *bm, double gu[3][3], const int dir);

  double cal_cs(double rho, double p, double gamma);

  /**
   * For debugging purposes only, make sure that the padding region of the righthand side
   * is always zero.
   */
  bool check_padding(double *dtu[], const unsigned int *sz);
};

#endif
