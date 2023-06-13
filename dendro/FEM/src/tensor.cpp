//
// Created by milinda on 1/15/17.
//


/**
 *
 * @author Milinda Fernando
 * @breif contains the utilities for tensor kronecker products for interpolations.
 *
 * */

#include "tensor.h"


/**
 * Along the Z axis
 * */
void DENDRO_TENSOR_AIIX_APPLY_ELEM (const int M, const double*  A, const double*  X, double*  Y) {
    int i, j, k;
    double d, e;
    for (i = 0; i < M; ++i) {
        d = A[i];
        for (j = 0; j < M * M; ++j) {
            Y[i * M * M + j] = d * X[j];
        }
        for (k = 1; k < M; ++k) {
            d = A[i + k * M];
            for (j = 0; j < M * M; ++j) {
                e = d * X[M * M * k + j];
                Y[i * M * M + j] += e;
            }
        }
    }
}

/**
 * Along the X axis
 * */

void DENDRO_TENSOR_IIAX_APPLY_ELEM(const int M, const double*  A, const double*  X, double*  Y)
{


    int i, j, k;
    double e;
    for (i = 0; i < M * M; i++) {
       // _mm_prefetch( (const char *)(X + (i+10)*M), 2);
        //_mm_prefetch( (const char *)(Y + (i+10)*M), 2);
        for (j=0; j<M; ++j) {
            e = 0;
            for (k=0; k<M; ++k) {
                e += X[i*M + k] * A[k*M + j];
            }
            Y[i*M + j] = e;
        }
    }

}

/**
 * Along the X axis. (in face interpolations. )
 * */

void DENDRO_TENSOR_IAX_APPLY_ELEM_2D(const int M, const double*  A, const double*  X, double*  Y)
{
    int i, j, k;
    double e;
    for (i = 0; i < M; i++) {
        // _mm_prefetch( (const char *)(X + (i+10)*M), 2);
        //_mm_prefetch( (const char *)(Y + (i+10)*M), 2);
        for (j=0; j<M; ++j) {
            e = 0;
            for (k=0; k<M; ++k) {
                e += X[i*M + k] * A[k*M + j];
            }
            Y[i*M + j] = e;
        }
    }

}


/**
 * Along the Y axis
 * */
void DENDRO_TENSOR_IAIX_APPLY_ELEM (const int M, const double*  A, const double*  X, double*  Y)
{
    int i, j, k, ib;
    double d, e;
    for (ib = 0; ib < M; ++ib) {
        for (i = 0; i < M; ++i) {
            d = A[i];
            for (j = 0; j < M; ++j) {
                Y[ib * M * M + i * M + j] =
                        d * X[ib * M * M + j];
            }
            for (k = 1; k < M; ++k) {
                d = A[i + k * M];
                for (j = 0; j < M; ++j) {
                    e = d * X[ib * M * M + M * k + j];
                    Y[ib * M * M + i * M + j] += e;
                }
            }
        }
    }
}

/**
 * Along the Y axis for (2D face interpolations. )
 * */
void DENDRO_TENSOR_AIX_APPLY_ELEM_2D (const int M, const double*  A, const double*  X, double*  Y)
{
    int i, j, k, ib=0;
    double d, e;
    //for (ib = 0; ib < M; ++ib) {
        for (i = 0; i < M; ++i) {
            d = A[i];
            for (j = 0; j < M; ++j) {
                Y[ib * M * M + i * M + j] =
                        d * X[ib * M * M + j];
            }
            for (k = 1; k < M; ++k) {
                d = A[i + k * M];
                for (j = 0; j < M; ++j) {
                    e = d * X[ib * M * M + M * k + j];
                    Y[ib * M * M + i * M + j] += e;
                }
            }
        }
    //}
}
