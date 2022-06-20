/* TwoPunctures:  File  "FuncAndJacobian.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "TPUtilities.h"
#include "TwoPunctures.h"

#define FAC sin(al)*sin(be)*sin(al)*sin(be)*sin(al)*sin(be)
/*#define FAC sin(al)*sin(be)*sin(al)*sin(be)*/
/*#define FAC 1*/

static inline CCTK_REAL min (CCTK_REAL const x, CCTK_REAL const y)
{
  return x<y ? x : y;
}

/* --------------------------------------------------------------------------*/
int
Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3)
{
  int i1 = i, j1 = j, k1 = k;

  if (i1 < 0)
    i1 = -(i1 + 1);
  if (i1 >= n1)
    i1 = 2 * n1 - (i1 + 1);

  if (j1 < 0)
    j1 = -(j1 + 1);
  if (j1 >= n2)
    j1 = 2 * n2 - (j1 + 1);

  if (k1 < 0)
    k1 = k1 + n3;
  if (k1 >= n3)
    k1 = k1 - n3;

  return ivar + nvar * (i1 + n1 * (j1 + n2 * k1));
}

/* --------------------------------------------------------------------------*/
void
allocate_derivs (derivs * v, int n)
{
  int m = n - 1;
  (*v).d0 = dvector (0, m);
  (*v).d1 = dvector (0, m);
  (*v).d2 = dvector (0, m);
  (*v).d3 = dvector (0, m);
  (*v).d11 = dvector (0, m);
  (*v).d12 = dvector (0, m);
  (*v).d13 = dvector (0, m);
  (*v).d22 = dvector (0, m);
  (*v).d23 = dvector (0, m);
  (*v).d33 = dvector (0, m);
}

/* --------------------------------------------------------------------------*/
void
free_derivs (derivs * v, int n)
{
  int m = n - 1;
  free_dvector ((*v).d0, 0, m);
  free_dvector ((*v).d1, 0, m);
  free_dvector ((*v).d2, 0, m);
  free_dvector ((*v).d3, 0, m);
  free_dvector ((*v).d11, 0, m);
  free_dvector ((*v).d12, 0, m);
  free_dvector ((*v).d13, 0, m);
  free_dvector ((*v).d22, 0, m);
  free_dvector ((*v).d23, 0, m);
  free_dvector ((*v).d33, 0, m);
}

/* --------------------------------------------------------------------------*/
void
Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs v)
{
  int i, j, k, ivar, N, *indx;
  CCTK_REAL *p, *dp, *d2p, *q, *dq, *r, *dr;

  N = maximum3 (n1, n2, n3);
  p = dvector (0, N);
  dp = dvector (0, N);
  d2p = dvector (0, N);
  q = dvector (0, N);
  dq = dvector (0, N);
  r = dvector (0, N);
  dr = dvector (0, N);
  indx = ivector (0, N);

  for (ivar = 0; ivar < nvar; ivar++)
  {
    for (k = 0; k < n3; k++)
    {				/* Calculation of Derivatives w.r.t. A-Dir. */
      for (j = 0; j < n2; j++)
      {				/* (Chebyshev_Zeros)*/
	for (i = 0; i < n1; i++)
	{
	  indx[i] = Index (ivar, i, j, k, nvar, n1, n2, n3);
	  p[i] = v.d0[indx[i]];
	}
	chebft_Zeros (p, n1, 0);
	chder (p, dp, n1);
	chder (dp, d2p, n1);
	chebft_Zeros (dp, n1, 1);
	chebft_Zeros (d2p, n1, 1);
	for (i = 0; i < n1; i++)
	{
	  v.d1[indx[i]] = dp[i];
	  v.d11[indx[i]] = d2p[i];
	}
      }
    }
    for (k = 0; k < n3; k++)
    {				/* Calculation of Derivatives w.r.t. B-Dir. */
      for (i = 0; i < n1; i++)
      {				/* (Chebyshev_Zeros)*/
	for (j = 0; j < n2; j++)
	{
	  indx[j] = Index (ivar, i, j, k, nvar, n1, n2, n3);
	  p[j] = v.d0[indx[j]];
	  q[j] = v.d1[indx[j]];
	}
	chebft_Zeros (p, n2, 0);
	chebft_Zeros (q, n2, 0);
	chder (p, dp, n2);
	chder (dp, d2p, n2);
	chder (q, dq, n2);
	chebft_Zeros (dp, n2, 1);
	chebft_Zeros (d2p, n2, 1);
	chebft_Zeros (dq, n2, 1);
	for (j = 0; j < n2; j++)
	{
	  v.d2[indx[j]] = dp[j];
	  v.d22[indx[j]] = d2p[j];
	  v.d12[indx[j]] = dq[j];
	}
      }
    }
    for (i = 0; i < n1; i++)
    {				/* Calculation of Derivatives w.r.t. phi-Dir. (Fourier)*/
      for (j = 0; j < n2; j++)
      {
	for (k = 0; k < n3; k++)
	{
	  indx[k] = Index (ivar, i, j, k, nvar, n1, n2, n3);
	  p[k] = v.d0[indx[k]];
	  q[k] = v.d1[indx[k]];
	  r[k] = v.d2[indx[k]];
	}
	fourft (p, n3, 0);
	fourder (p, dp, n3);
	fourder2 (p, d2p, n3);
	fourft (dp, n3, 1);
	fourft (d2p, n3, 1);
	fourft (q, n3, 0);
	fourder (q, dq, n3);
	fourft (dq, n3, 1);
	fourft (r, n3, 0);
	fourder (r, dr, n3);
	fourft (dr, n3, 1);
	for (k = 0; k < n3; k++)
	{
	  v.d3[indx[k]] = dp[k];
	  v.d33[indx[k]] = d2p[k];
	  v.d13[indx[k]] = dq[k];
	  v.d23[indx[k]] = dr[k];
	}
      }
    }
  }
  free_dvector (p, 0, N);
  free_dvector (dp, 0, N);
  free_dvector (d2p, 0, N);
  free_dvector (q, 0, N);
  free_dvector (dq, 0, N);
  free_dvector (r, 0, N);
  free_dvector (dr, 0, N);
  free_ivector (indx, 0, N);
}

/* --------------------------------------------------------------------------*/
void F_of_v (int nvar, int n1, int n2, int n3, derivs v, CCTK_REAL *F, derivs u)
{
  /*      Calculates the left hand sides of the non-linear equations F_m(v_n)=0*/
  /*      and the function u (u.d0[]) as well as its derivatives*/
  /*      (u.d1[], u.d2[], u.d3[], u.d11[], u.d12[], u.d13[], u.d22[], u.d23[], u.d33[])*/
  /*      at interior points and at the boundaries "+/-"*/
  
  
  CCTK_REAL *sources;
  sources= (CCTK_REAL *) calloc(n1*n2*n3, sizeof(CCTK_REAL));
  if (TPID::use_sources)
  {
    CCTK_REAL *s_x, *s_y, *s_z;
    
    s_x    = (CCTK_REAL *) calloc(n1*n2*n3, sizeof(CCTK_REAL));
    s_y    = (CCTK_REAL *) calloc(n1*n2*n3, sizeof(CCTK_REAL));
    s_z    = (CCTK_REAL *) calloc(n1*n2*n3, sizeof(CCTK_REAL));
    
    #pragma omp parallel for num_threads(TP_OMP_THREADS) collapse(3)
    for (int i = 0; i < n1; i++)
      for (int j = 0; j < n2; j++)
        for (int k = 0; k < n3; k++)
        {

          CCTK_REAL *values;
          derivs U;
          values = dvector (0, nvar - 1);
          allocate_derivs (&U, nvar);

          CCTK_REAL al, be, A, B, X, R, x, r, phi, y, z, Am1;
          CCTK_INT i3D = Index(0,i,j,k,1,n1,n2,n3);

          al = Pih * (2 * i + 1) / n1;
          A = -cos (al);
          be = Pih * (2 * j + 1) / n2;
          B = -cos (be);
          phi = 2. * Pi * k / n3;

          Am1 = A - 1;
          for (int ivar = 0; ivar < nvar; ivar++)
          {
            int indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
            U.d0[ivar] = Am1 * v.d0[indx];        /* U*/
            U.d1[ivar] = v.d0[indx] + Am1 * v.d1[indx];        /* U_A*/
            U.d2[ivar] = Am1 * v.d2[indx];        /* U_B*/
            U.d3[ivar] = Am1 * v.d3[indx];        /* U_3*/
            U.d11[ivar] = 2 * v.d1[indx] + Am1 * v.d11[indx];        /* U_AA*/
            U.d12[ivar] = v.d2[indx] + Am1 * v.d12[indx];        /* U_AB*/
            U.d13[ivar] = v.d3[indx] + Am1 * v.d13[indx];        /* U_AB*/
            U.d22[ivar] = Am1 * v.d22[indx];        /* U_BB*/
            U.d23[ivar] = Am1 * v.d23[indx];        /* U_B3*/
            U.d33[ivar] = Am1 * v.d33[indx];        /* U_33*/
          }
          /* Calculation of (X,R) and*/
          /* (U_X, U_R, U_3, U_XX, U_XR, U_X3, U_RR, U_R3, U_33)*/
          AB_To_XR (nvar, A, B, &X, &R, U);
          /* Calculation of (x,r) and*/
          /* (U, U_x, U_r, U_3, U_xx, U_xr, U_x3, U_rr, U_r3, U_33)*/
          C_To_c (nvar, X, R, &(s_x[i3D]), &r, U);
          /* Calculation of (y,z) and*/
          /* (U, U_x, U_y, U_z, U_xx, U_xy, U_xz, U_yy, U_yz, U_zz)*/
          rx3_To_xyz (nvar, s_x[i3D], r, phi, &(s_y[i3D]), &(s_z[i3D]), U);

          free_dvector (values, 0, nvar - 1);
          free_derivs (&U, nvar);

        }
    printf("This code not set up for sources.\n");
    //Set_Rho_ADM(cctkGH, n1*n2*n3, sources, s_x, s_y, s_z);
    free(s_z);
    free(s_y);
    free(s_x);
  }
  else
  {
    #pragma omp parallel for num_threads(TP_OMP_THREADS) collapse(3)
    for (int i = 0; i < n1; i++)
      for (int j = 0; j < n2; j++)
        for (int k = 0; k < n3; k++)
          sources[Index(0,i,j,k,1,n1,n2,n3)]=0.0;

  }
    

  Derivatives_AB3 (nvar, n1, n2, n3, v);
  CCTK_REAL psi, psi2, psi4, psi7, r_plus, r_minus;
  FILE *debugfile = NULL;
  //if (do_residuum_debug_output && CCTK_MyProc(cctkGH) == 0)
  if (TPID::do_residuum_debug_output == 0)
  {
    debugfile = fopen("res.dat", "w");
    assert(debugfile);
  }
  
  // Note: You cannot omp parallel collapse this loop there is some data dep. that I am not sure of. (I think some deriv computations going on. ) - Milinda. 
  #pragma omp parallel for num_threads(TP_OMP_THREADS) collapse(3)
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      for (int k = 0; k < n3; k++)
      {
        CCTK_REAL al, be, A, B, X, R, x, r, phi, y, z, Am1;
        al = Pih * (2 * i + 1) / n1;
        A = -cos (al);
        be = Pih * (2 * j + 1) / n2;
        B = -cos (be);
        phi = 2. * Pi * k / n3;

        Am1 = A - 1;

        CCTK_REAL *values;
        derivs U;
        values = dvector (0, nvar - 1);
        allocate_derivs (&U, nvar);

        for (int ivar = 0; ivar < nvar; ivar++)
        {
          int indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
          U.d0[ivar] = Am1 * v.d0[indx];        /* U*/
          U.d1[ivar] = v.d0[indx] + Am1 * v.d1[indx];        /* U_A*/
          U.d2[ivar] = Am1 * v.d2[indx];        /* U_B*/
          U.d3[ivar] = Am1 * v.d3[indx];        /* U_3*/
          U.d11[ivar] = 2 * v.d1[indx] + Am1 * v.d11[indx];        /* U_AA*/
          U.d12[ivar] = v.d2[indx] + Am1 * v.d12[indx];        /* U_AB*/
          U.d13[ivar] = v.d3[indx] + Am1 * v.d13[indx];        /* U_AB*/
          U.d22[ivar] = Am1 * v.d22[indx];        /* U_BB*/
          U.d23[ivar] = Am1 * v.d23[indx];        /* U_B3*/
          U.d33[ivar] = Am1 * v.d33[indx];        /* U_33*/
        }
        /* Calculation of (X,R) and*/
        /* (U_X, U_R, U_3, U_XX, U_XR, U_X3, U_RR, U_R3, U_33)*/
        AB_To_XR (nvar, A, B, &X, &R, U);
        /* Calculation of (x,r) and*/
        /* (U, U_x, U_r, U_3, U_xx, U_xr, U_x3, U_rr, U_r3, U_33)*/
        C_To_c (nvar, X, R, &x, &r, U);
        /* Calculation of (y,z) and*/
        /* (U, U_x, U_y, U_z, U_xx, U_xy, U_xz, U_yy, U_yz, U_zz)*/
        rx3_To_xyz (nvar, x, r, phi, &y, &z, U);
        NonLinEquations (sources[Index(0,i,j,k,1,n1,n2,n3)],
                         A, B, X, R, x, r, phi, y, z, U, values);
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          int indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
          F[indx] = values[ivar] * FAC;
          /* if ((i<5) && ((j<5) || (j>n2-5)))*/
          /*     F[indx] = 0.0;*/
          u.d0[indx] = U.d0[ivar];        /*  U*/
          u.d1[indx] = U.d1[ivar];        /*      U_x*/
          u.d2[indx] = U.d2[ivar];        /*      U_y*/
          u.d3[indx] = U.d3[ivar];        /*      U_z*/
          u.d11[indx] = U.d11[ivar];        /*      U_xx*/
          u.d12[indx] = U.d12[ivar];        /*      U_xy*/
          u.d13[indx] = U.d13[ivar];        /*      U_xz*/
          u.d22[indx] = U.d22[ivar];        /*      U_yy*/
          u.d23[indx] = U.d23[ivar];        /*      U_yz*/
          u.d33[indx] = U.d33[ivar];        /*      U_zz*/
        }
        
        free_dvector (values, 0, nvar - 1);
        free_derivs (&U, nvar);

      }

  free(sources);
  return;
  
}

/* --------------------------------------------------------------------------*/
void J_times_dv (int nvar, int n1, int n2, int n3, derivs dv, CCTK_REAL *Jdv, derivs u)
{
  /* Calculates the left hand sides of the non-linear equations F_m(v_n)=0*/
  /* and the function u (u.d0[]) as well as its derivatives*/
  /* (u.d1[], u.d2[], u.d3[], u.d11[], u.d12[], u.d13[], u.d22[], u.d23[], u.d33[])*/
  /* at interior points and at the boundaries "+/-"*/

  
  Derivatives_AB3 (nvar, n1, n2, n3, dv);

  #pragma omp parallel for num_threads(TP_OMP_THREADS) collapse(3)
  for (int i = 0; i < n1; i++)
  {
    for (int j = 0; j < n2; j++)
    {
      for (int k = 0; k < n3; k++)
      {
        CCTK_REAL al, be, A, B, X, R, x, r, phi, y, z, Am1, *values;
        values = dvector (0, nvar - 1);
        derivs dU, U;
        allocate_derivs (&dU, nvar);
        allocate_derivs (&U, nvar);

        al = Pih * (2 * i + 1) / n1;
        A = -cos (al);
        be = Pih * (2 * j + 1) / n2;
        B = -cos (be);
        phi = 2. * Pi * k / n3;

        Am1 = A - 1;
	      for (int ivar = 0; ivar < nvar; ivar++)
	      {
	        int indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
	        dU.d0[ivar] = Am1 * dv.d0[indx];	/* dU*/
	        dU.d1[ivar] = dv.d0[indx] + Am1 * dv.d1[indx];	/* dU_A*/
	        dU.d2[ivar] = Am1 * dv.d2[indx];	/* dU_B*/
	        dU.d3[ivar] = Am1 * dv.d3[indx];	/* dU_3*/
	        dU.d11[ivar] = 2 * dv.d1[indx] + Am1 * dv.d11[indx];	/* dU_AA*/
	        dU.d12[ivar] = dv.d2[indx] + Am1 * dv.d12[indx];	/* dU_AB*/
	        dU.d13[ivar] = dv.d3[indx] + Am1 * dv.d13[indx];	/* dU_AB*/
	        dU.d22[ivar] = Am1 * dv.d22[indx];	/* dU_BB*/
	        dU.d23[ivar] = Am1 * dv.d23[indx];	/* dU_B3*/
	        dU.d33[ivar] = Am1 * dv.d33[indx];	/* dU_33*/
	        U.d0[ivar] = u.d0[indx];	/* U   */
	        U.d1[ivar] = u.d1[indx];	/* U_x*/
	        U.d2[ivar] = u.d2[indx];	/* U_y*/
	        U.d3[ivar] = u.d3[indx];	/* U_z*/
	        U.d11[ivar] = u.d11[indx];	/* U_xx*/
	        U.d12[ivar] = u.d12[indx];	/* U_xy*/
	        U.d13[ivar] = u.d13[indx];	/* U_xz*/
	        U.d22[ivar] = u.d22[indx];	/* U_yy*/
	        U.d23[ivar] = u.d23[indx];	/* U_yz*/
	        U.d33[ivar] = u.d33[indx];	/* U_zz*/
	      }
	      /* Calculation of (X,R) and*/
	      /* (dU_X, dU_R, dU_3, dU_XX, dU_XR, dU_X3, dU_RR, dU_R3, dU_33)*/
	      AB_To_XR (nvar, A, B, &X, &R, dU);
	      /* Calculation of (x,r) and*/
	      /* (dU, dU_x, dU_r, dU_3, dU_xx, dU_xr, dU_x3, dU_rr, dU_r3, dU_33)*/
	      C_To_c (nvar, X, R, &x, &r, dU);
	      /* Calculation of (y,z) and*/
	      /* (dU, dU_x, dU_y, dU_z, dU_xx, dU_xy, dU_xz, dU_yy, dU_yz, dU_zz)*/
	      rx3_To_xyz (nvar, x, r, phi, &y, &z, dU);
	      LinEquations (A, B, X, R, x, r, phi, y, z, dU, U, values);
        for (int ivar = 0; ivar < nvar; ivar++)
        {
          int indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
          Jdv[indx] = values[ivar] * FAC;
        }

        free_dvector (values, 0, nvar - 1);
        free_derivs (&dU, nvar);
        free_derivs (&U, nvar);

      }
    }

  }

}

/* --------------------------------------------------------------------------*/
void JFD_times_dv (int i, int j, int k, int nvar, int n1, int n2, int n3, derivs dv, derivs u, CCTK_REAL *values)
{
  /* Calculates rows of the vector 'J(FD)*dv'.*/
  /* First row to be calculated: row = Index(0,      i, j, k; nvar, n1, n2, n3)*/
  /* Last  row to be calculated: row = Index(nvar-1, i, j, k; nvar, n1, n2, n3)*/
  /* These rows are stored in the vector JFDdv[0] ... JFDdv[nvar-1].*/
  int ivar, indx;
  CCTK_REAL al, be, A, B, X, R, x, r, phi, y, z, Am1;
  CCTK_REAL sin_al, sin_al_i1, sin_al_i2, sin_al_i3, cos_al;
  CCTK_REAL sin_be, sin_be_i1, sin_be_i2, sin_be_i3, cos_be;
  CCTK_REAL dV0, dV1, dV2, dV3, dV11, dV12, dV13, dV22, dV23, dV33,
    ha, ga, ga2, hb, gb, gb2, hp, gp, gp2, gagb, gagp, gbgp;
  derivs dU, U;

  allocate_derivs (&dU, nvar);
  allocate_derivs (&U, nvar);

  if (k < 0)
    k = k + n3;
  if (k >= n3)
    k = k - n3;

  ha = Pi / n1;			/* ha: Stepsize with respect to (al)*/
  al = ha * (i + 0.5);
  A = -cos (al);
  ga = 1 / ha;
  ga2 = ga * ga;

  hb = Pi / n2;			/* hb: Stepsize with respect to (be)*/
  be = hb * (j + 0.5);
  B = -cos (be);
  gb = 1 / hb;
  gb2 = gb * gb;
  gagb = ga * gb;

  hp = 2 * Pi / n3;		/* hp: Stepsize with respect to (phi)*/
  phi = hp * k;
  gp = 1 / hp;
  gp2 = gp * gp;
  gagp = ga * gp;
  gbgp = gb * gp;


  sin_al = sin (al);
  sin_be = sin (be);
  sin_al_i1 = 1 / sin_al;
  sin_be_i1 = 1 / sin_be;
  sin_al_i2 = sin_al_i1 * sin_al_i1;
  sin_be_i2 = sin_be_i1 * sin_be_i1;
  sin_al_i3 = sin_al_i1 * sin_al_i2;
  sin_be_i3 = sin_be_i1 * sin_be_i2;
  cos_al = -A;
  cos_be = -B;

  Am1 = A - 1;
  for (ivar = 0; ivar < nvar; ivar++)
  {
    int iccc = Index (ivar, i, j, k, nvar, n1, n2, n3),
      ipcc = Index (ivar, i + 1, j, k, nvar, n1, n2, n3),
      imcc = Index (ivar, i - 1, j, k, nvar, n1, n2, n3),
      icpc = Index (ivar, i, j + 1, k, nvar, n1, n2, n3),
      icmc = Index (ivar, i, j - 1, k, nvar, n1, n2, n3),
      iccp = Index (ivar, i, j, k + 1, nvar, n1, n2, n3),
      iccm = Index (ivar, i, j, k - 1, nvar, n1, n2, n3),
      icpp = Index (ivar, i, j + 1, k + 1, nvar, n1, n2, n3),
      icmp = Index (ivar, i, j - 1, k + 1, nvar, n1, n2, n3),
      icpm = Index (ivar, i, j + 1, k - 1, nvar, n1, n2, n3),
      icmm = Index (ivar, i, j - 1, k - 1, nvar, n1, n2, n3),
      ipcp = Index (ivar, i + 1, j, k + 1, nvar, n1, n2, n3),
      imcp = Index (ivar, i - 1, j, k + 1, nvar, n1, n2, n3),
      ipcm = Index (ivar, i + 1, j, k - 1, nvar, n1, n2, n3),
      imcm = Index (ivar, i - 1, j, k - 1, nvar, n1, n2, n3),
      ippc = Index (ivar, i + 1, j + 1, k, nvar, n1, n2, n3),
      impc = Index (ivar, i - 1, j + 1, k, nvar, n1, n2, n3),
      ipmc = Index (ivar, i + 1, j - 1, k, nvar, n1, n2, n3),
      immc = Index (ivar, i - 1, j - 1, k, nvar, n1, n2, n3);
    /* Derivatives of (dv) w.r.t. (al,be,phi):*/
    dV0 = dv.d0[iccc];
    dV1 = 0.5 * ga * (dv.d0[ipcc] - dv.d0[imcc]);
    dV2 = 0.5 * gb * (dv.d0[icpc] - dv.d0[icmc]);
    dV3 = 0.5 * gp * (dv.d0[iccp] - dv.d0[iccm]);
    dV11 = ga2 * (dv.d0[ipcc] + dv.d0[imcc] - 2 * dv.d0[iccc]);
    dV22 = gb2 * (dv.d0[icpc] + dv.d0[icmc] - 2 * dv.d0[iccc]);
    dV33 = gp2 * (dv.d0[iccp] + dv.d0[iccm] - 2 * dv.d0[iccc]);
    dV12 =
      0.25 * gagb * (dv.d0[ippc] - dv.d0[ipmc] + dv.d0[immc] - dv.d0[impc]);
    dV13 =
      0.25 * gagp * (dv.d0[ipcp] - dv.d0[imcp] + dv.d0[imcm] - dv.d0[ipcm]);
    dV23 =
      0.25 * gbgp * (dv.d0[icpp] - dv.d0[icpm] + dv.d0[icmm] - dv.d0[icmp]);
    /* Derivatives of (dv) w.r.t. (A,B,phi):*/
    dV11 = sin_al_i3 * (sin_al * dV11 - cos_al * dV1);
    dV12 = sin_al_i1 * sin_be_i1 * dV12;
    dV13 = sin_al_i1 * dV13;
    dV22 = sin_be_i3 * (sin_be * dV22 - cos_be * dV2);
    dV23 = sin_be_i1 * dV23;
    dV1 = sin_al_i1 * dV1;
    dV2 = sin_be_i1 * dV2;
    /* Derivatives of (dU) w.r.t. (A,B,phi):*/
    dU.d0[ivar] = Am1 * dV0;
    dU.d1[ivar] = dV0 + Am1 * dV1;
    dU.d2[ivar] = Am1 * dV2;
    dU.d3[ivar] = Am1 * dV3;
    dU.d11[ivar] = 2 * dV1 + Am1 * dV11;
    dU.d12[ivar] = dV2 + Am1 * dV12;
    dU.d13[ivar] = dV3 + Am1 * dV13;
    dU.d22[ivar] = Am1 * dV22;
    dU.d23[ivar] = Am1 * dV23;
    dU.d33[ivar] = Am1 * dV33;

    indx = Index (ivar, i, j, k, nvar, n1, n2, n3);
    U.d0[ivar] = u.d0[indx];	/* U   */
    U.d1[ivar] = u.d1[indx];	/* U_x*/
    U.d2[ivar] = u.d2[indx];	/* U_y*/
    U.d3[ivar] = u.d3[indx];	/* U_z*/
    U.d11[ivar] = u.d11[indx];	/* U_xx*/
    U.d12[ivar] = u.d12[indx];	/* U_xy*/
    U.d13[ivar] = u.d13[indx];	/* U_xz*/
    U.d22[ivar] = u.d22[indx];	/* U_yy*/
    U.d23[ivar] = u.d23[indx];	/* U_yz*/
    U.d33[ivar] = u.d33[indx];	/* U_zz*/
  }
  /* Calculation of (X,R) and*/
  /* (dU_X, dU_R, dU_3, dU_XX, dU_XR, dU_X3, dU_RR, dU_R3, dU_33)*/
  AB_To_XR (nvar, A, B, &X, &R, dU);
  /* Calculation of (x,r) and*/
  /* (dU, dU_x, dU_r, dU_3, dU_xx, dU_xr, dU_x3, dU_rr, dU_r3, dU_33)*/
  C_To_c (nvar, X, R, &x, &r, dU);
  /* Calculation of (y,z) and*/
  /* (dU, dU_x, dU_y, dU_z, dU_xx, dU_xy, dU_xz, dU_yy, dU_yz, dU_zz)*/
  rx3_To_xyz (nvar, x, r, phi, &y, &z, dU);
  LinEquations (A, B, X, R, x, r, phi, y, z, dU, U, values);
  for (ivar = 0; ivar < nvar; ivar++)
    values[ivar] *= FAC;

  free_derivs (&dU, nvar);
  free_derivs (&U, nvar);
}

/* --------------------------------------------------------------------------*/
void SetMatrix_JFD (int nvar, int n1, int n2, int n3, derivs u, int *ncols, int **cols, CCTK_REAL **Matrix)
{
  int N1, N2, N3;
  int ntotal = nvar * n1 * n2 * n3;
  CCTK_REAL *values;
  derivs dv;
  values = dvector (0, nvar - 1);
  allocate_derivs (&dv, ntotal);
  
  N1 = n1 - 1;
  N2 = n2 - 1;
  N3 = n3 - 1;

  #pragma omp parallel for num_threads(TP_OMP_THREADS) collapse(4)
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < n2; j++)
      for (int k = 0; k < n3; k++)
    	for (int ivar = 0; ivar < nvar; ivar++)
	    {
	      int row = Index (ivar, i, j, k, nvar, n1, n2, n3);
	      ncols[row] = 0;
	      dv.d0[row] = 0;
	    }

  for (int i = 0; i < n1; i++)
  {
    for (int j = 0; j < n2; j++)
    {
      for (int k = 0; k < n3; k++)
      {
	      for (int ivar = 0; ivar < nvar; ivar++)
	      {
          int column = Index (ivar, i, j, k, nvar, n1, n2, n3);
          dv.d0[column] = 1;

          int i_0 = maximum2 (0, i - 1);
          int i_1 = minimum2 (N1, i + 1);
          int j_0 = maximum2 (0, j - 1);
          int j_1 = minimum2 (N2, j + 1);
          int k_0 = k - 1;
          int k_1 = k + 1;
          /*i_0 = 0;
					i_1 = N1;
					j_0 = 0;
					j_1 = N2;
					k_0 = 0;
					k_1 = N3;*/
          
          for (int i1 = i_0; i1 <= i_1; i1++)
          {
            for (int j1 = j_0; j1 <= j_1; j1++)
            {
              for (int k1 = k_0; k1 <= k_1; k1++)
              {
                JFD_times_dv (i1, j1, k1, nvar, n1, n2, n3,dv, u, values);
                for (int ivar1 = 0; ivar1 < nvar; ivar1++)
                {
                  if (values[ivar1] != 0)
                  {
                    int row = Index (ivar1, i1, j1, k1, nvar, n1, n2, n3);
                    int mcol = ncols[row];
                    cols[row][mcol] = column;
                    Matrix[row][mcol] = values[ivar1];
                    ncols[row] += 1;
                  }
                }
              }
            }
          }

	        dv.d0[column] = 0;
	      }
      }
    }
  }
  free_derivs (&dv, ntotal);
  free_dvector (values, 0, nvar - 1);
}

/* --------------------------------------------------------------------------*/
/* Calculates the value of v at an arbitrary position (A,B,phi)*/
CCTK_REAL PunctEvalAtArbitPosition (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi, int nvar, int n1, int n2, int n3)
{
  int i, j, k, N;
  CCTK_REAL *p, *values1, **values2, result;

  N = maximum3 (n1, n2, n3);
  p = dvector (0, N);
  values1 = dvector (0, N);
  values2 = dmatrix (0, N, 0, N);

  for (k = 0; k < n3; k++)
  {
    for (j = 0; j < n2; j++)
    {
      for (i = 0; i < n1; i++)
	p[i] = v[ivar + nvar * (i + n1 * (j + n2 * k))];
      chebft_Zeros (p, n1, 0);
      values2[j][k] = chebev (-1, 1, p, n1, A);
    }
  }

  for (k = 0; k < n3; k++)
  {
    for (j = 0; j < n2; j++)
      p[j] = values2[j][k];
    chebft_Zeros (p, n2, 0);
    values1[k] = chebev (-1, 1, p, n2, B);
  }

  fourft (values1, n3, 0);
  result = fourev (values1, n3, phi);

  free_dvector (p, 0, N);
  free_dvector (values1, 0, N);
  free_dmatrix (values2, 0, N, 0, N);

  return result;
}

/* --------------------------------------------------------------------------*/
void
calculate_derivs (int i, int j, int k, int ivar, int nvar, int n1, int n2,
		  int n3, derivs v, derivs vv)
{
  CCTK_REAL al = Pih * (2 * i + 1) / n1, be = Pih * (2 * j + 1) / n2,
    sin_al = sin (al), sin2_al = sin_al * sin_al, cos_al = cos (al),
    sin_be = sin (be), sin2_be = sin_be * sin_be, cos_be = cos (be);

  vv.d0[0] = v.d0[Index (ivar, i, j, k, nvar, n1, n2, n3)];
  vv.d1[0] = v.d1[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_al;
  vv.d2[0] = v.d2[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_be;
  vv.d3[0] = v.d3[Index (ivar, i, j, k, nvar, n1, n2, n3)];
  vv.d11[0] = v.d11[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin2_al
    + v.d1[Index (ivar, i, j, k, nvar, n1, n2, n3)] * cos_al;
  vv.d12[0] =
    v.d12[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_al * sin_be;
  vv.d13[0] = v.d13[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_al;
  vv.d22[0] = v.d22[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin2_be
    + v.d2[Index (ivar, i, j, k, nvar, n1, n2, n3)] * cos_be;
  vv.d23[0] = v.d23[Index (ivar, i, j, k, nvar, n1, n2, n3)] * sin_be;
  vv.d33[0] = v.d33[Index (ivar, i, j, k, nvar, n1, n2, n3)];
}

/* --------------------------------------------------------------------------*/
CCTK_REAL
interpol (CCTK_REAL a, CCTK_REAL b, CCTK_REAL c, derivs v)
{
  return v.d0[0]
    + a * v.d1[0] + b * v.d2[0] + c * v.d3[0]
    + 0.5 * a * a * v.d11[0] + a * b * v.d12[0] + a * c * v.d13[0]
    + 0.5 * b * b * v.d22[0] + b * c * v.d23[0] + 0.5 * c * c * v.d33[0];
}

/* --------------------------------------------------------------------------*/
static CCTK_REAL
clamp_pm_one (CCTK_REAL val)
{
  return val < -1 ? -1 : val > 1 ? 1 : val;
}

/* --------------------------------------------------------------------------*/
/* Calculates the value of v at an arbitrary position (x,y,z)*/
CCTK_REAL
PunctTaylorExpandAtArbitPosition (int ivar, int nvar, int n1,
                                  int n2, int n3, derivs v, CCTK_REAL x, CCTK_REAL y,
                                  CCTK_REAL z)
{
  CCTK_REAL xs, ys, zs, rs2, phi, X, R, A, B, al, be, aux1, aux2, a, b, c,
    result, Ui;
  int i, j, k;
  derivs vv;
  allocate_derivs (&vv, 1);

  xs = x / TPID::par_b;
  ys = y / TPID::par_b;
  zs = z / TPID::par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)
    phi += 2 * Pi;

  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));

  /* Note: Range of R = asin(Q) is [0,pi] for Q in [0,1] */
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0)
    R = Pi - R;

  A = clamp_pm_one( 2 * tanh (0.5 * X) - 1 );

  /* Note: Range of R/2 - pi/4 is [ -pi/4, pi/4 ] and so range of tan
   * is [-1,1], for R in [0,pi]. */
  B = clamp_pm_one( tan (0.5 * R - Piq) );
  al = Pi - acos (A);
  be = Pi - acos (B);

  i = rint (al * n1 / Pi - 0.5);
  j = rint (be * n2 / Pi - 0.5);
  k = rint (0.5 * phi * n3 / Pi);

  a = al - Pi * (i + 0.5) / n1;
  b = be - Pi * (j + 0.5) / n2;
  c = phi - 2 * Pi * k / n3;

  calculate_derivs (i, j, k, ivar, nvar, n1, n2, n3, v, vv);
  result = interpol (a, b, c, vv);
  free_derivs (&vv, 1);

  Ui = (A - 1) * result;

  assert( std::isfinite( Ui ) );

  return Ui;
}

/* --------------------------------------------------------------------------*/
/* Calculates the value of v at an arbitrary position (x,y,z)*/
CCTK_REAL
PunctIntPolAtArbitPosition (int ivar, int nvar, int n1,
			    int n2, int n3, derivs v, CCTK_REAL x, CCTK_REAL y,
			    CCTK_REAL z)
{
  CCTK_REAL xs, ys, zs, rs2, phi, X, R, A, B, aux1, aux2, result, Ui;

  xs = x / TPID::par_b;
  ys = y / TPID::par_b;
  zs = z / TPID::par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)
    phi += 2 * Pi;

  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0)
    R = Pi - R;

  A = 2 * tanh (0.5 * X) - 1;
  B = tan (0.5 * R - Piq);

  result = PunctEvalAtArbitPosition (v.d0, ivar, A, B, phi, nvar, n1, n2, n3);

  Ui = (A - 1) * result;

  assert( std::isfinite( Ui ) );

  return Ui;
}


//////////////////////////////////////////////////////
/// Fast Spectral Interpolation Routine Stuff
//////////////////////////////////////////////////////


/* Calculates the value of v at an arbitrary position (A,B,phi)* using the fast routine */
CCTK_REAL
PunctEvalAtArbitPositionFast (CCTK_REAL *v, int ivar, CCTK_REAL A, CCTK_REAL B, CCTK_REAL phi, int nvar, int n1, int n2, int n3)
{
  int i, j, k, N;
  CCTK_REAL *p, *values1, **values2, result;
  // VASILIS: Nothing should be changed in this routine. This is used by PunctIntPolAtArbitPositionFast

  N = maximum3 (n1, n2, n3);

  p = dvector (0, N);
  values1 = dvector (0, N);
  values2 = dmatrix (0, N, 0, N);

  for (k = 0; k < n3; k++)
  {
    for (j = 0; j < n2; j++)
    {
      for (i = 0; i < n1; i++) p[i] = v[ivar + nvar * (i + n1 * (j + n2 * k))];
      //      chebft_Zeros (p, n1, 0);
      values2[j][k] = chebev (-1, 1, p, n1, A);
    }
  }

  for (k = 0; k < n3; k++)
  {
    for (j = 0; j < n2; j++) p[j] = values2[j][k];
    //    chebft_Zeros (p, n2, 0);
    values1[k] = chebev (-1, 1, p, n2, B);
  }

  //  fourft (values1, n3, 0);
  result = fourev (values1, n3, phi);

  free_dvector (p, 0, N);
  free_dvector (values1, 0, N);
  free_dmatrix (values2, 0, N, 0, N);

  return result;
  //  */
  //  return 0.;
}


// --------------------------------------------------------------------------*/
// Calculates the value of v at an arbitrary position (x,y,z) if the spectral coefficients are known //
/* --------------------------------------------------------------------------*/
CCTK_REAL
PunctIntPolAtArbitPositionFast (int ivar, int nvar, int n1,
			    int n2, int n3, derivs v, CCTK_REAL x, CCTK_REAL y,
			    CCTK_REAL z)
{
  CCTK_REAL xs, ys, zs, rs2, phi, X, R, A, B, aux1, aux2, result, Ui;
  // VASILIS: Here the struct derivs v refers to the spectral coeffiecients of variable v not the variable v itself

  xs = x / TPID::par_b;
  ys = y / TPID::par_b;
  zs = z / TPID::par_b;
  rs2 = ys * ys + zs * zs;
  phi = atan2 (z, y);
  if (phi < 0)
    phi += 2 * Pi;

  aux1 = 0.5 * (xs * xs + rs2 - 1);
  aux2 = sqrt (aux1 * aux1 + rs2);
  X = asinh (sqrt (aux1 + aux2));
  R = asin (min(1.0, sqrt (-aux1 + aux2)));
  if (x < 0)
    R = Pi - R;

  A = 2 * tanh (0.5 * X) - 1;
  B = tan (0.5 * R - Piq);

  result = PunctEvalAtArbitPositionFast (v.d0, ivar, A, B, phi, nvar, n1, n2, n3);

  Ui = (A - 1) * result;

  return Ui;
}

// Evaluates the spectral expansion coefficients of v
void SpecCoef(int n1, int n2, int n3, int ivar, CCTK_REAL *v, CCTK_REAL *cf)
{
  // VASILIS: Here v is a pointer to the values of the variable v at the collocation points and cf_v a pointer to the spectral coefficients that this routine calculates

	int i, j, k, N, n, l;
	CCTK_REAL *p, ***values3, ***values4;

	N=maximum3(n1,n2,n3);
	p=dvector(0,N);
	values3=d3tensor(0,n1,0,n2,0,n3);
	values4=d3tensor(0,n1,0,n2,0,n3);



	      // Caclulate values3[n,j,k] = a_n^{j,k} = (sum_i^(n1-1) f(A_i,B_j,phi_k) Tn(-A_i))/k_n , k_n = N/2 or N
	      for(k=0;k<n3;k++) {
		for(j=0;j<n2;j++) {

		  for(i=0;i<n1;i++) p[i]=v[ivar + (i + n1 * (j + n2 * k))];

		  chebft_Zeros(p,n1,0);
		  for (n=0;n<n1;n++)	{
		    values3[n][j][k] = p[n];
		  }
		}
	      }

	      // Caclulate values4[n,l,k] = a_{n,l}^{k} = (sum_j^(n2-1) a_n^{j,k} Tn(B_j))/k_l , k_l = N/2 or N

	      for (n = 0; n < n1; n++){
		for(k=0;k<n3;k++) {
		  for(j=0;j<n2;j++) p[j]=values3[n][j][k];
		  chebft_Zeros(p,n2,0);
		  for (l = 0; l < n2; l++){
		  values4[n][l][k] = p[l];
		  }
		}
	      }

	      // Caclulate coefficients  a_{n,l,m} = (sum_k^(n3-1) a_{n,m}^{k} fourier(phi_k))/k_m , k_m = N/2 or N
	      for (i = 0; i < n1; i++){
		for (j = 0; j < n2; j++){
		  for(k=0;k<n3;k++) p[k]=values4[i][j][k];
		  fourft(p,n3,0);
		  for (k = 0; k<n3; k++){
		    cf[ivar + (i + n1 * (j + n2 * k))] = p[k];
		  }
		}
	      }

	      free_dvector(p,0,N);
	      free_d3tensor(values3,0,n1,0,n2,0,n3);
	      free_d3tensor(values4,0,n1,0,n2,0,n3);

}
