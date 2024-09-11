/* TwoPunctures:  File  "utilities.h"*/
#pragma once
#include <math.h>

#include "mpi.h"
#include "parameters.h"

static void code_exit(const char *s) {
    printf("%s\n", s);
    exit(-1);
}

#define CCTK_REAL  double
#define CCTK_INT   int
#define CCTK_ERROR code_exit

#define Pi  3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164 /* Pi/2*/
#define Piq 0.78539816339744830961566084582 /* Pi/4*/

#define TINY 1.0e-20
#define SWAP(a, b)   \
    {                \
        temp = (a);  \
        (a)  = (b);  \
        (b)  = temp; \
    }

#define nrerror       TP_nrerror
#define ivector       TP_ivector
#define dvector       TP_dvector
#define imatrix       TP_imatrix
#define dmatrix       TP_dmatrix
#define d3tensor      TP_d3tensor
#define free_ivector  TP_free_ivector
#define free_dvector  TP_free_dvector
#define free_imatrix  TP_free_imatrix
#define free_dmatrix  TP_free_dmatrix
#define free_d3tensor TP_free_d3tensor

void nrerror(char error_text[]);
int *ivector(long nl, long nh);
CCTK_REAL *dvector(long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
CCTK_REAL **dmatrix(long nrl, long nrh, long ncl, long nch);
CCTK_REAL ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl,
                      long ndh);
void free_ivector(int *v, long nl, long nh);
void free_dvector(CCTK_REAL *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(CCTK_REAL **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(CCTK_REAL ***t, long nrl, long nrh, long ncl, long nch,
                   long ndl, long ndh);

int minimum2(int i, int j);
int minimum3(int i, int j, int k);
int maximum2(int i, int j);
int maximum3(int i, int j, int k);
int pow_int(int mantisse, int exponent);

void chebft_Zeros(CCTK_REAL u[], int n, int inv);
void chebft_Extremes(CCTK_REAL u[], int n, int inv);
void chder(CCTK_REAL *c, CCTK_REAL *cder, int n);
CCTK_REAL chebev(CCTK_REAL a, CCTK_REAL b, CCTK_REAL c[], int m, CCTK_REAL x);
void fourft(CCTK_REAL *u, int N, int inv);
void fourder(CCTK_REAL u[], CCTK_REAL du[], int N);
void fourder2(CCTK_REAL u[], CCTK_REAL d2u[], int N);
CCTK_REAL fourev(CCTK_REAL *u, int N, CCTK_REAL x);

CCTK_REAL norm1(CCTK_REAL *v, int n);
CCTK_REAL norm2(CCTK_REAL *v, int n);
CCTK_REAL scalarproduct(CCTK_REAL *v, CCTK_REAL *w, int n);

extern unsigned int TP_OMP_THREADS;
extern MPI_Comm TP_MPI_COMM;
extern unsigned int TP_SOLVE_VALID;
