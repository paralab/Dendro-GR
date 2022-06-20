/* TwoPunctures:  File  "TwoPunctures.h"*/
#ifndef SFCSORTBENCH_TWOPUNCTURES_H
#define SFCSORTBENCH_TWOPUNCTURES_H


#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

#define TAYLOR_EXPANSION 0
#define EVALUATION 1

#define LAPSE_ANTISYMMETRIC 0
#define LAPSE_AVERAGED 1
#define LAPSE_PSIN 2
#define LAPSE_BROWNSVILLE 3

#define CS_FACTOR 0
#define CS_FACTOR_DERIVS 1
#define CS_FACTOR_SECOND_DERIVS 2

#define MT_STATIC_CONFORMAL 0
#define MT_STANDARD 1

#include "TPUtilities.h"

typedef struct DERIVS
{
  CCTK_REAL *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;

/*
Files of "TwoPunctures":
	TwoPunctures.c
	FuncAndJacobian.c
	CoordTransf.c
	Equations.c
	Newton.c
	utilities.c (see utilities.h)
**************************

*/

void TwoPunctures(const double xx1, const double yy1, const double zz1, double *vars,  double *mp, double *mm, double *mp_adm, double *mm_adm, double *E, double *J1, double *J2, double *J3);

/* Routines in  "TwoPunctures.c"*/
CCTK_REAL TestSolution (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R, CCTK_REAL phi);
void TestVector_w (CCTK_REAL *par, int nvar, int n1, int n2, int n3, CCTK_REAL *w);

/* Routines in  "FuncAndJacobian.c"*/
int Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void allocate_derivs (derivs * v, int n);
void free_derivs (derivs * v, int n);
void Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs v);
void F_of_v (int nvar, int n1, int n2, int n3, derivs v,CCTK_REAL *F, derivs u);
void J_times_dv (int nvar, int n1, int n2, int n3, derivs dv, CCTK_REAL *Jdv, derivs u);
void JFD_times_dv (int i, int j, int k, int nvar, int n1, int n2, int n3, derivs dv, derivs u, CCTK_REAL *values);
void SetMatrix_JFD (int nvar, int n1, int n2, int n3, derivs u, int *ncols, int **cols, CCTK_REAL **Matrix);
CCTK_REAL PunctEvalAtArbitPosition (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi, int nvar, int n1, int n2, int n3);
void calculate_derivs (int i, int j, int k, int ivar, int nvar, int n1,int n2, int n3, derivs v, derivs vv);
CCTK_REAL interpol (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c, derivs v); 
CCTK_REAL PunctTaylorExpandAtArbitPosition (int ivar, int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
CCTK_REAL PunctIntPolAtArbitPosition (int ivar, int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
void SpecCoef(int n1, int n2, int n3, int ivar, CCTK_REAL *v, CCTK_REAL *cf);
CCTK_REAL PunctEvalAtArbitPositionFast (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi, int nvar, int n1, int n2, int n3);
CCTK_REAL PunctIntPolAtArbitPositionFast (int ivar, int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);

/* Routines in  "CoordTransf.c"*/
void AB_To_XR (int nvar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL *X,CCTK_REAL *R, derivs U);
void C_To_c (int nvar, CCTK_REAL X, CCTK_REAL R, CCTK_REAL *x,CCTK_REAL *r, derivs U);
void rx3_To_xyz (int nvar, CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi, CCTK_REAL *y,CCTK_REAL *z, derivs U);

/* Routines in  "Equations.c"*/
CCTK_REAL BY_KKofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
void BY_Aijofxyz (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL Aij[3][3]);
void NonLinEquations (CCTK_REAL rho_adm,CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R, CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi, CCTK_REAL y, CCTK_REAL z, derivs U, CCTK_REAL *values);
void LinEquations (CCTK_REAL A, CCTK_REAL B, CCTK_REAL X, CCTK_REAL R, CCTK_REAL x, CCTK_REAL r, CCTK_REAL phi, CCTK_REAL y, CCTK_REAL z, derivs dU, derivs U, CCTK_REAL *values);


/* Routines in  "Newton.c"*/
void TestRelax ( int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL *dv);
void Newton ( int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL tol, int itmax);

/**@brief write the tp solve to a file*/
void TPStore(double *mp, double *mm, double *mp_adm, double *mm_adm, double *E, double *J1, double *J2, double *J3, const char* fprefix);
/**@brief restore the tp solve from file. */
void TPRestore(CCTK_REAL*& F, derivs& u, derivs& v, derivs& cf_v, const char* fprefix,bool mpi_bcast=true);


/* 
 27: -1.325691774825335e-03
 37: -1.325691778944117e-03
 47: -1.325691778942711e-03
 
 17: -1.510625972641537e-03
 21: -1.511443006977708e-03
 27: -1.511440785153687e-03
 37: -1.511440809549005e-03
 39: -1.511440809597588e-03
 */
#endif
