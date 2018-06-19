//
// Created by milinda on 1/19/17.
//

/**
 *
 * @author Milinda Fernando
 * School of Computing University of Utah.
 * @brief Constains lapack routines such as linear system solve, eigen solve to build the interpolation matrices.
 *
 *
 * */


#ifndef SFCSORTBENCH_LAPAC_H
#define SFCSORTBENCH_LAPAC_H

#include "lapacke.h"
#include <cstring>
#include <iostream>
namespace lapack
{

/**
 *  @brief: Wrapper for LAPACK DGESV solver for AX=B. Parameters are given below.
 *  @param[in] n : number of rows or columns of linear system
 *  @param[in] nrhs: number of right hand sides.
 *  @param[in] A: matrix A
 *  @param[in] lda: leading dimention of the array A
 *  @param[in] B: matrix B
 *  @param[out] X:  matrix X (solution)
 *  @param[in]  ldb:  leading dimention of B
 *  @param[out] info:  returns the status of the solve.
 */
inline void lapack_DGESV(unsigned int n , unsigned int nrhs, const double * A, unsigned int lda, double * B, double * X, unsigned int ldb, unsigned int info)
    {
        lapack_int * ipiv = (lapack_int *)malloc(n*sizeof(lapack_int)) ;
        memcpy(X,B,sizeof(double)*n*nrhs);

        double * L=new double[n*n];
        memcpy(L,A,sizeof(double)*n*n);

        info = LAPACKE_dgesv( LAPACK_ROW_MAJOR, n, nrhs, L, lda, ipiv,
                              X,ldb);

        if( info > 0 ) {
            printf( "The diagonal element of the triangular factor of A,\n" );
            printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
            printf( "the solution could not be computed.\n" );

        }
        if (info <0) {printf(" lapack linear solve failed. \n"); }
        delete [] L;
        free(ipiv);

        return ;
    }

/**
 *  @brief: Wrapper for LAPACK DGESV compute eigen values of a square matrix of A. Parameters are given below.
 *  @param[in] n : number of rows or columns of linear system
 *  @param[in] A: matrix A
 *  @param[in] lda: leading dimention of the array A
 *  @param[out] wr: real part of eigen values
 *  @param[out] vs eigen vectors
 *  @param[out] info:  returns the status of the solve.
 */


inline void lapack_DSYEV(unsigned int n, const double * A, unsigned int lda, double * wr,double * vs,unsigned int info )
{
    memcpy(vs,A,sizeof(double)*n*n);

    info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',n,vs,lda,wr);

    /*for(unsigned int i=0;i<n;i++)
        std::cout<<" i: "<<i<<" eig : "<<wr[i]<<std::endl;

    for(unsigned int i=0;i<n;i++)
    {
        for(unsigned int j=0;j<n;j++)
        {
            std::cout<<vs[i*(n)+j]<<" ";
        }

        std::cout<<std::endl;
    }*/

    if(info!=0) std::cout<<"lapack eigen solve failed. "<<std::endl;
    return;

}






}// end of namespace


#endif //SFCSORTBENCH_LAPAC_H
