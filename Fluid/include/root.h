#ifndef ROOT_H_
#define ROOT_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

namespace froot{
  int rbrac(int (*func)(double *, double, double *), double *x1, double *x2,
            double *par, int trace);
  int bisec(int (*func)(double *, double, double *),
            double *rtb, double x1, double x2,
            double tol, double *fpar, int trace);
  int rtsafe(int (*func)(double *, double *, double , double *),
             double *rts, double x1, double x2, double tol,
             double *fpar, int trace);
  int zbrent(int(*func)(double *, double, double *),
             double *rts, double x1, double x2, double tol, double *fpar);
  void check_finite_1d(double *f, int *nx, int *rc);
  void quiet_check_finite_1d(double *f, int *nx, int *rc);
}

#endif
