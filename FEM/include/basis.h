//
// Created by milinda on 1/12/17.
//

/**
 *
 * @author Milinda Fernando,
 * School of Computing, University of Utah.
 * @breif Contains implementation of Jacobi polynomials, and Jacobi Gauss-Lobatto quadrature  and Gauss quadrature.
 * All the basis functions that is hoped to be used in for interpolation and integration will be coded in the basis.h and basis.cpp file.
 *
 * */


#ifndef SFCSORTBENCH_JACOBIPOLYNOMIAL_H
#define SFCSORTBENCH_JACOBIPOLYNOMIAL_H


#include "assert.h"
#include <iostream>
#include <math.h>
#include <cmath>
#include "cstring"
#ifdef WITH_BLAS_LAPACK
#include "lapac.h"
#endif

#ifndef WITH_BLAS_LAPACK
template <typename T>
void printArray_1D(T *a, int length)
{
    for (int i = 0; i < length; i++) { std::cout<<a[i]<<" "; }
    std::cout<<std::endl;
}


template <typename T>
void printArray_2D(T *a, int length1,int length2)
{
    for (int i = 0; i < length1; i++) {
        for (int j = 0; j < length2; j++) {
            std::cout << a[i * length2 + j] << " ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

#endif



namespace basis {

/** Evaluate \c N'th order Jacobi Polynomial of type (\c alpha, \c beta ) > -1.
 *
 * Also (\c alpha + \c beta <> -1).
 *
 * This is C reimplemetation of the function in JacobiP from Nodal
 * Discontinuous Galerkin Methods by Jan S. Hesthaven and Tim Warburton.
 *
 * \param [in]  alpha Jacobi polynomial parameter
 * \param [in]  beta  Jacobi polynomial parameter
 * \param [in]  N     Order of the polynomial
 * \param [in]  np    size of x and p.
 * \param [in]  x     Row vector of locations to evaluate the polynomial.
 * \param [out] p     Row vector of the evaluated \a N'th order polynomial
 *
 * \note The matrix \a P is assumed to have the correct size of (1, len(x)).
 *       Also then polynomials are normalized to be orthonormal.
 */
    void jacobip(double alpha, double beta, unsigned int N,
                 double *x, double *p , unsigned int np);

/** Evaluate  the derivative of the \c N'th order Jacobi Polynomial of
 * type (\c alpha, \c beta ) > -1.
 *
 * This is C reimplemetation of the function in JacobiP from Nodal
 * Discontinuous Galerkin Methods by Jan S. Hesthaven and Tim Warburton.
 *
 * \param [in]  alpha Jacobi polynomial parameter
 * \param [in]  beta  Jacobi polynomial parameter
 * \param [in]  N     Order of the polynomial
 * \param [in]  np     size of x and p.
 * \param [in]  x     Row vector of locations to evaluate the polynomial.
 * \param [out] dp    Row vector of the evaluated derivative of the
 *                    \a N'th order polynomial
 *
 * \note The matrix \a P is assumed to have the correct size of (1, len(x)).
 *       Also then polynomials are normalized to be orthonormal.
 */
    void gradjacobip(double alpha, double beta, int N,double *x, double *dp, unsigned int np);

/** Compute the \c N'th order Gauss Lobatto quadrature points and weights.
 *
 * Where (\c alpha, \c beta ) > -1 and (\c alpha + \c beta <> -1).
 *
 * \param [in]  alpha Jacobi polynomial parameter
 * \param [in]  beta  Jacobi polynomial parameter
 * \param [in]  N     Order of the polynomial
 * \param [out] x     Row vector of Gauss Lobatto node locations
 * \param [out] w     Row vector of Gauss Lobatto weights
 *
 * \note The matrices \a x and \a w are assumed to have the correct size
 * of (1, N+1).  Also the node location will be between [-1,1].
 */
    void jacobiglq(double alpha, double beta, int N,
                   double *x, double *w);

/** Compute the \c N'th order Gauss quadrature points and weights.
 *
 * Where (\c alpha, \c beta ) > -1.
 *
 * This is C reimplemetation of the function in JacobiGQ from Nodal
 * Discontinuous Galerkin Methods by Jan S. Hesthaven and Tim Warburton.
 *
 * \param [in]  alpha Jacobi polynomial parameter
 * \param [in]  beta  Jacobi polynomial parameter
 * \param [in]  N     Order of the polynomial
 * \param [out] x     Row vector of Gauss node locations
 * \param [out] w     Row vector of Gauss weights
 *
 * \note The matrices \a x and \a w are assumed to have the correct size
 * of (1, N+1).  Also the node locaiton will be between [-1,1].
 */
    void jacobigq(double alpha, double beta, int N,
                  double *x, double *w);



}

#endif //SFCSORTBENCH_JACOBIPOLYNOMIAL_H
