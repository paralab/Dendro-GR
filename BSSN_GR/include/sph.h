/**
 * @file spherical_harmonics.h
 * @brief spherical harmonics 
 * @version 0.1
 * @date 2023-11-08
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#pragma once
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "gsl/gsl_specfunc.h"
#include "gsl/gsl_integration.h"

/**
 * @brief real spherical harmonics function
 * 
 * @param l polar mode
 * @param m azimuthal mode
 * @param theta polar angle 
 * @param phi  azimuthal angle
 * @return double 
 */
double real_spherical_harmonic(int l, int m, double theta, double phi);

/**
 * @brief real spherical harmonics function derivative w.r.t. theta
 * 
 * @param l polar mode
 * @param m azimuthal mode
 * @param theta polar angle 
 * @param phi  azimuthal angle
 * @param dorder derivative order (currently supported only for 1)
 * @return double 
 */
double real_spherical_harmonic_d_theta(int l, int m, double theta, double phi, int dorder);

/**
 * @brief real spherical harmonics function derivative w.r.t. phi
 * @param l polar mode
 * @param m azimuthal mode
 * @param theta polar angle 
 * @param phi  azimuthal angle
 * @param dorder derivative order (currently supported only for 1)
 * @return double 
 */
double real_spherical_harmonic_d_phi(int l, int m, double theta, double phi, int dorder);


inline double real_spherical_harmonic(int l, int m, double theta, double phi)
{
    const int abs_m = abs(m);
    assert(abs_m<=l);
    const int csphase  = (abs_m%2==0) ? 1 : -1;
    if(m==0)
        return gsl_sf_legendre_sphPlm(l,m, cos(theta));
    else if (m>0)
        return csphase * sqrt(2) * gsl_sf_legendre_sphPlm(l,m, cos(theta)) * cos(m * phi);
    else
        return csphase * sqrt(2) * gsl_sf_legendre_sphPlm(l, abs_m, cos(theta)) * sin( abs_m * phi);
}


inline double real_spherical_harmonic_d_theta(int l, int m, double theta, double phi, int dorder)
{
    const int abs_m = abs(m);
    assert(abs_m<=l);
    assert(dorder==1);

    if(l==0)
        return 0;
        
    const int csphase  = (abs_m%2==0) ? 1 : -1;
    long f1 = 1;
    for (int i=1; i<=(l-abs_m); ++i)
        f1 *= i;
    
    long f2 = 1;
    for (int i=1; i<=(l+abs_m); ++i)
        f2 *= i;

    const double t1   = cos(theta);
    const double u_lm = sqrt((2 * l + 1) * f1 / (4 * M_PI * f2));
    
    const double f3 = (abs_m+1 <=l) ? 0.5 * gsl_sf_legendre_Plm(l, abs_m + 1, t1) : 0 ;
    const double f4 = -0.5 * (l + abs_m) * (l-abs_m +1) * gsl_sf_legendre_Plm(l, abs(abs_m-1), t1);
    
    const double Plm_deriv = f4 + f3;
    
    if(m==0)
        return u_lm * Plm_deriv;
    else if (m>0)
        return csphase * sqrt(2) * u_lm * Plm_deriv * cos(m * phi);
    else
        return csphase * sqrt(2) * u_lm * Plm_deriv * sin( abs_m * phi);

}


inline double real_spherical_harmonic_d_phi(int l, int m, double theta, double phi, int dorder)
{
    const int abs_m = abs(m);
    assert(abs_m<=l);
    assert(dorder==1);
    const int csphase  = (abs_m%2==0) ? 1 : -1;

    if(m==0)
        return 0;
    else if (m>0)
        return csphase * sqrt(2) * gsl_sf_legendre_sphPlm(l,m, cos(theta)) * (-sin(m * phi) * m);
    else
        return csphase * sqrt(2) * gsl_sf_legendre_sphPlm(l, abs_m, cos(theta)) * (cos( abs_m * phi) * abs_m);

}