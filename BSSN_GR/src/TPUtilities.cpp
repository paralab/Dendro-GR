/* TwoPunctures:  File  "utilities.c"*/

#include "TPUtilities.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

unsigned int TP_OMP_THREADS = 16;
MPI_Comm TP_MPI_COMM        = MPI_COMM_WORLD;
unsigned int TP_SOLVE_VALID = 0;

/*---------------------------------------------------------------------------*/
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int *retval;

    retval = (int *)malloc(sizeof(int) * (nh - nl + 1));
    if (retval == NULL) CCTK_ERROR("allocation failure in ivector()");

    return retval - nl;
}

/*---------------------------------------------------------------------------*/
CCTK_REAL *dvector(long nl, long nh)
/* allocate a CCTK_REAL vector with subscript range v[nl..nh] */
{
    CCTK_REAL *retval;

    retval = (CCTK_REAL *)malloc(sizeof(CCTK_REAL) * (nh - nl + 1));
    if (retval == NULL) CCTK_ERROR("allocation failure in dvector()");

    return retval - nl;
}

/*---------------------------------------------------------------------------*/
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int **retval;

    retval = (int **)malloc(sizeof(int *) * (nrh - nrl + 1));
    if (retval == NULL) CCTK_ERROR("allocation failure (1) in imatrix()");

    /* get all memory for the matrix in on chunk */
    retval[0] = (int *)malloc(sizeof(int) * (nrh - nrl + 1) * (nch - ncl + 1));
    if (retval[0] == NULL) CCTK_ERROR("allocation failure (2) in imatrix()");

    /* apply column and row offsets */
    retval[0] -= ncl;
    retval -= nrl;

    /* slice chunk into rows */
    long width = (nch - ncl + 1);
    for (long i = nrl + 1; i <= nrh; i++) retval[i] = retval[i - 1] + width;
    assert(retval[nrh] - retval[nrl] == (nrh - nrl) * width);

    return retval;
}

/*---------------------------------------------------------------------------*/
CCTK_REAL **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    CCTK_REAL **retval;

    retval = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *) * (nrh - nrl + 1));
    if (retval == NULL) CCTK_ERROR("allocation failure (1) in dmatrix()");

    /* get all memory for the matrix in on chunk */
    retval[0] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL) * (nrh - nrl + 1) *
                                    (nch - ncl + 1));
    if (retval[0] == NULL) CCTK_ERROR("allocation failure (2) in dmatrix()");

    /* apply column and row offsets */
    retval[0] -= ncl;
    retval -= nrl;

    /* slice chunk into rows */
    long width = (nch - ncl + 1);
    for (long i = nrl + 1; i <= nrh; i++) retval[i] = retval[i - 1] + width;
    assert(retval[nrh] - retval[nrl] == (nrh - nrl) * width);

    return retval;
}

/*---------------------------------------------------------------------------*/
CCTK_REAL ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl,
                      long ndh)
/* allocate a CCTK_REAL 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    CCTK_REAL ***retval;

    /* get memory for index structures */
    retval = (CCTK_REAL ***)malloc(sizeof(CCTK_REAL **) * (nrh - nrl + 1));
    if (retval == NULL) CCTK_ERROR("allocation failure (1) in d3tensor()");

    retval[0] = (CCTK_REAL **)malloc(sizeof(CCTK_REAL *) * (nrh - nrl + 1) *
                                     (nch - ncl + 1));
    if (retval[0] == NULL) CCTK_ERROR("allocation failure (2) in d3tensor()");

    /* get all memory for the tensor in on chunk */
    retval[0][0] = (CCTK_REAL *)malloc(sizeof(CCTK_REAL) * (nrh - nrl + 1) *
                                       (nch - ncl + 1) * (ndh - ndl + 1));
    if (retval[0][0] == NULL)
        CCTK_ERROR("allocation failure (3) in d3tensor()");

    /* apply all offsets */
    retval[0][0] -= ndl;
    retval[0] -= ncl;
    retval -= nrl;

    /* slice chunk into rows and columns */
    long width = (nch - ncl + 1);
    long depth = (ndh - ndl + 1);
    for (long j = ncl + 1; j <= nch; j++) { /* first row of columns */
        retval[nrl][j] = retval[nrl][j - 1] + depth;
    }
    assert(retval[nrl][nch] - retval[nrl][ncl] == (nch - ncl) * depth);
    for (long i = nrl + 1; i <= nrh; i++) {
        retval[i] = retval[i - 1] + width;
        retval[i][ncl] =
            retval[i - 1][ncl] + width * depth; /* first cell in column */
        for (long j = ncl + 1; j <= nch; j++) {
            retval[i][j] = retval[i][j - 1] + depth;
        }
        assert(retval[i][nch] - retval[i][ncl] == (nch - ncl) * depth);
    }
    assert(retval[nrh] - retval[nrl] == (nrh - nrl) * width);
    assert(&retval[nrh][nch][ndh] - &retval[nrl][ncl][ndl] ==
           (nrh - nrl + 1) * (nch - ncl + 1) * (ndh - ndl + 1) - 1);

    return retval;
}

/*--------------------------------------------------------------------------*/
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free(v + nl);
}

/*--------------------------------------------------------------------------*/
void free_dvector(CCTK_REAL *v, long nl, long nh)
/* free an double vector allocated with dvector() */
{
    free(v + nl);
}

/*--------------------------------------------------------------------------*/
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
    free(m[nrl] + ncl);
    free(m + nrl);
}

/*--------------------------------------------------------------------------*/
void free_dmatrix(CCTK_REAL **m, long nrl, long nrh, long ncl, long nch)
/* free a CCTK_REAL matrix allocated by dmatrix() */
{
    free(m[nrl] + ncl);
    free(m + nrl);
}

/*--------------------------------------------------------------------------*/
void free_d3tensor(CCTK_REAL ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh)
/* free a CCTK_REAL d3tensor allocated by d3tensor() */
{
    free(t[nrl][ncl] + ndl);
    free(t[nrl] + ncl);
    free(t + nrl);
}

/*--------------------------------------------------------------------------*/
int minimum2(int i, int j) {
    int result = i;
    if (j < result) result = j;
    return result;
}

/*-------------------------------------------------------------------------*/
int minimum3(int i, int j, int k) {
    int result = i;
    if (j < result) result = j;
    if (k < result) result = k;
    return result;
}

/*--------------------------------------------------------------------------*/
int maximum2(int i, int j) {
    int result = i;
    if (j > result) result = j;
    return result;
}

/*--------------------------------------------------------------------------*/
int maximum3(int i, int j, int k) {
    int result = i;
    if (j > result) result = j;
    if (k > result) result = k;
    return result;
}

/*--------------------------------------------------------------------------*/
int pow_int(int mantisse, int exponent) {
    int i, result = 1;

    for (i = 1; i <= exponent; i++) result *= mantisse;

    return result;
}

/*--------------------------------------------------------------------------*/
void chebft_Zeros(CCTK_REAL u[], int n, int inv)
/* eq. 5.8.7 and 5.8.8 at x = (5.8.4) of 2nd edition C++ NR */
{
    int k, j, isignum;
    CCTK_REAL fac, sum, Pion, *c;

    c    = dvector(0, n);
    Pion = Pi / n;
    if (inv == 0) {
        fac     = 2.0 / n;
        isignum = 1;
        for (j = 0; j < n; j++) {
            sum = 0.0;
            for (k = 0; k < n; k++) sum += u[k] * cos(Pion * j * (k + 0.5));
            c[j]    = fac * sum * isignum;
            isignum = -isignum;
        }
    } else {
        for (j = 0; j < n; j++) {
            sum     = -0.5 * u[0];
            isignum = 1;
            for (k = 0; k < n; k++) {
                sum += u[k] * cos(Pion * (j + 0.5) * k) * isignum;
                isignum = -isignum;
            }
            c[j] = sum;
        }
    }
    for (j = 0; j < n; j++)
#if 0
    if (fabs(c[j]) < 5.e-16)
      u[j] = 0.0;
    else
#endif
        u[j] = c[j];
    free_dvector(c, 0, n);
}

/* --------------------------------------------------------------------------*/

void chebft_Extremes(CCTK_REAL u[], int n, int inv)
/* eq. 5.8.7 and 5.8.8 at x = (5.8.5) of 2nd edition C++ NR */
{
    int k, j, isignum, N = n - 1;
    CCTK_REAL fac, sum, PioN, *c;

    c    = dvector(0, N);
    PioN = Pi / N;
    if (inv == 0) {
        fac     = 2.0 / N;
        isignum = 1;
        for (j = 0; j < n; j++) {
            sum = 0.5 * (u[0] + u[N] * isignum);
            for (k = 1; k < N; k++) sum += u[k] * cos(PioN * j * k);
            c[j]    = fac * sum * isignum;
            isignum = -isignum;
        }
        c[N] = 0.5 * c[N];
    } else {
        for (j = 0; j < n; j++) {
            sum     = -0.5 * u[0];
            isignum = 1;
            for (k = 0; k < n; k++) {
                sum += u[k] * cos(PioN * j * k) * isignum;
                isignum = -isignum;
            }
            c[j] = sum;
        }
    }
    for (j = 0; j < n; j++) u[j] = c[j];
    free_dvector(c, 0, N);
}

/* --------------------------------------------------------------------------*/

void chder(CCTK_REAL *c, CCTK_REAL *cder, int n) {
    int j;

    cder[n]     = 0.0;
    cder[n - 1] = 0.0;
    for (j = n - 2; j >= 0; j--) cder[j] = cder[j + 2] + 2 * (j + 1) * c[j + 1];
}

/* --------------------------------------------------------------------------*/
CCTK_REAL
chebev(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c[], int m, CCTK_REAL x)
/* eq. 5.8.11 of C++ NR (2nd ed) */
{
    int j;
    CCTK_REAL djp2, djp1, dj; /* d_{j+2}, d_{j+1} and d_j */
    CCTK_REAL y;

    /* rescale input to lie within [-1,1] */
    y  = 2 * (x - 0.5 * (b + a)) / (b - a);

    dj = djp1 = 0;
    for (j = m - 1; j >= 1; j--) {
        /* advance the coefficients */
        djp2 = djp1;
        djp1 = dj;
        dj   = 2 * y * djp1 - djp2 + c[j];
    }

    return y * dj - djp1 + 0.5 * c[0];
}

/* --------------------------------------------------------------------------*/
void fourft(CCTK_REAL *u, int N, int inv)
/* a (slow) Fourier transform, seems to be just eq. 12.1.6 and 12.1.9 of C++ NR
   (2nd ed) */
{
    int l, k, iy, M;
    CCTK_REAL x, x1, fac, Pi_fac, *a, *b;

    M = N / 2;
    a = dvector(0, M);
    b = dvector(
        1, M); /* Actually: b=vector(1,M-1) but this is problematic if M=1*/
    fac    = 1. / M;
    Pi_fac = Pi * fac;
    if (inv == 0) {
        for (l = 0; l <= M; l++) {
            a[l] = 0;
            if (l > 0 && l < M) b[l] = 0;
            x1 = Pi_fac * l;
            for (k = 0; k < N; k++) {
                x = x1 * k;
                a[l] += fac * u[k] * cos(x);
                if (l > 0 && l < M) b[l] += fac * u[k] * sin(x);
            }
        }
        u[0] = a[0];
        u[M] = a[M];
        for (l = 1; l < M; l++) {
            u[l]     = a[l];
            u[l + M] = b[l];
        }
    } else {
        a[0] = u[0];
        a[M] = u[M];
        for (l = 1; l < M; l++) {
            a[l] = u[l];
            b[l] = u[M + l];
        }
        iy = 1;
        for (k = 0; k < N; k++) {
            u[k] = 0.5 * (a[0] + a[M] * iy);
            x1   = Pi_fac * k;
            for (l = 1; l < M; l++) {
                x = x1 * l;
                u[k] += a[l] * cos(x) + b[l] * sin(x);
            }
            iy = -iy;
        }
    }
    free_dvector(a, 0, M);
    free_dvector(b, 1, M);
}

/* -----------------------------------------*/
void fourder(CCTK_REAL u[], CCTK_REAL du[], int N) {
    int l, M, lpM;

    M     = N / 2;
    du[0] = 0.;
    du[M] = 0.;
    for (l = 1; l < M; l++) {
        lpM     = l + M;
        du[l]   = u[lpM] * l;
        du[lpM] = -u[l] * l;
    }
}

/* -----------------------------------------*/
void fourder2(CCTK_REAL u[], CCTK_REAL d2u[], int N) {
    int l, l2, M, lpM;

    d2u[0] = 0.;
    M      = N / 2;
    for (l = 1; l <= M; l++) {
        l2     = l * l;
        lpM    = l + M;
        d2u[l] = -u[l] * l2;
        if (l < M) d2u[lpM] = -u[lpM] * l2;
    }
}

/* ----------------------------------------- */
CCTK_REAL
fourev(CCTK_REAL *u, int N, CCTK_REAL x) {
    int l, M = N / 2;
    CCTK_REAL xl, result;

    result = 0.5 * (u[0] + u[M] * cos(x * M));
    for (l = 1; l < M; l++) {
        xl = x * l;
        result += u[l] * cos(xl) + u[M + l] * sin(xl);
    }
    return result;
}

/* ------------------------------------------------------------------------*/
CCTK_REAL
norm1(CCTK_REAL *v, int n) {
    int i;
    CCTK_REAL result = -1;

    for (i = 0; i < n; i++)
        if (fabs(v[i]) > result) result = fabs(v[i]);

    return result;
}

/* -------------------------------------------------------------------------*/
CCTK_REAL
norm2(CCTK_REAL *v, int n) {
    int i;
    CCTK_REAL result = 0;

    for (i = 0; i < n; i++) result += v[i] * v[i];

    return sqrt(result);
}

/* -------------------------------------------------------------------------*/
CCTK_REAL
scalarproduct(CCTK_REAL *v, CCTK_REAL *w, int n) {
    int i;
    CCTK_REAL result = 0;

    for (i = 0; i < n; i++) result += v[i] * w[i];

    return result;
}

/* -------------------------------------------------------------------------*/
