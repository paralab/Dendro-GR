/**
 * @file trumpet.h
 * @brief Trumpet solution functions.
 * @version 0.1
 * @date 2021-07-25
 *
 * @copyright Copyright (c) 2021
 *
 */

#pragma once
#include <cmath>
#include <iostream>

#include "grDef.h"
#include "parameters.h"
namespace trumpet_data {

void trumpetData(const double xx1, const double yy1, const double zz1,
                 double *var);
void bndcnd(double h, double &x, double y[], double dydx[]);
void derivs(double x, double y[], double dydx[]);
void hunt(double xx[], int n, double x, int *jlo);
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
          double yerr[], void (*derivs)(double, double[], double[]));
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs)(double, double[], double[]));
void odeint(double ystart[], int nvar, double x1, double x2, double eps,
            double h1, double hmin, int *nok, int *nbad,
            void (*derivs)(double, double[], double[]),
            void (*rkqs)(double[], double[], int, double *, double, double,
                         double[], double *, double *,
                         void (*)(double, double[], double[])),
            int kount);
double interpolation3(double xp[], double yp[], int np, double xb,
                      int *n_nearest_pt);
double interpolation4(double xp[], double yp[], int np, double xb,
                      int *n_nearest_pt);

}  // end of namespace trumpet_data
