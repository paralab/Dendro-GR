/* TwoPunctures:  File  "TwoPunctures.c"*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <ctype.h>
#include "parameters.h"
#include "grUtils.h"
#include "TwoPunctures.h"

using namespace bssn;

inline double EXTEND(double M, double r) {
  return ( M * (3./8 * pow(r, 4) / pow(TPID::TP_Extend_Radius, 5) -
                 5./4 * pow(r, 2) / pow(TPID::TP_Extend_Radius, 3) +
                 15./8 / TPID::TP_Extend_Radius));
}

/* -------------------------------------------------------------------*/
void TwoPunctures(const double xx1, const double yy1, const double zz1, double *vars,double *mp, double *mm, double *mp_adm, double *mm_adm,double *E, double *J1, double *J2, double *J3)
{

  * mp = TPID::par_m_plus;
  * mm = TPID::par_m_minus;

  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;

  int rank, npes;
  MPI_Comm_size(TP_MPI_COMM,&npes);
  MPI_Comm_rank(TP_MPI_COMM,&rank);

  int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;

  int const nvar = 1, n1 = TPID::npoints_A, n2 = TPID::npoints_B, n3 = TPID::npoints_phi;

  int const ntotal = n1 * n2 * n3 * nvar;
  
  #if 1
    int percent10 = 0;
  #endif
  static CCTK_REAL *F = NULL;
  static derivs u, v, cf_v;
  CCTK_REAL admMass;
  
  if(!F)
    TPRestore(F,u,v,cf_v,TPID::FILE_PREFIX.c_str());

  if (TPID::grid_setup_method == TAYLOR_EXPANSION)
  {
    gsm = GSM_Taylor_expansion;
  }
  else if (TPID::grid_setup_method == EVALUATION)
  {
    gsm = GSM_evaluation;
  }
  else
  {
    if(!rank)
      printf("internal error. unknown grid_setup_method = %d\n", TPID::grid_setup_method);
    
    exit(-1);
  }

  antisymmetric_lapse = (TPID::initial_lapse == LAPSE_ANTISYMMETRIC) ? 1 : 0;
  averaged_lapse = (TPID::initial_lapse == LAPSE_AVERAGED) ? 1 : 0;
	pmn_lapse = (TPID::initial_lapse == LAPSE_PSIN) ? 1 : 0 ;
  if (pmn_lapse)
  {
    /*if(!rank)
      printf("Setting initial lapse to psi ^ %f profile.\n", TPID::initial_lapse_psi_exponent);*/
  }
		
  
  brownsville_lapse = (TPID::initial_lapse == LAPSE_BROWNSVILLE) ? 1 : 0;
  
  if (brownsville_lapse)
  {
    /*if(!rank)
     printf( "Setting initial lapse to a Brownsville-style profile with exp %f.", TPID::initial_lapse_psi_exponent);*/
  }
    
  
  /*if(!rank)
    printf("Interpolating result.\n");*/

 /* Apparently in Cactus, you can choose to save one of the following
  * to 3D arrays:
  *    psi
  *    psi + first derivatives
  *    psi + first derivatives + second derivatives.
  *
  * For now I only provide the option to save psi.
  */

  int metric_type = MT_STANDARD;
  int conformal_storage = CS_FACTOR;
  int conformal_state;

  if (metric_type == MT_STATIC_CONFORMAL) {
    if (conformal_storage == CS_FACTOR) {
      conformal_state = 1;
    } else if (conformal_storage == CS_FACTOR_DERIVS) {
      conformal_state = 2;
    } else if (conformal_storage == CS_FACTOR_SECOND_DERIVS) {
      conformal_state = 3;
    }
  } else {
    conformal_state = 0;
  }

  CCTK_REAL xx, yy, zz;
  const double x= xx1; 
  const double y= yy1;
  const double z= zz1; 

  xx = x - TPID::center_offset[0];
  yy = y - TPID::center_offset[1];
  zz = z - TPID::center_offset[2];

  CCTK_REAL r_plus = sqrt(pow(xx - TPID::par_b, 2) + pow(yy, 2) + pow(zz, 2));
  CCTK_REAL r_minus = sqrt(pow(xx + TPID::par_b, 2) + pow(yy, 2) + pow(zz, 2));

  CCTK_REAL U;
  switch (gsm)
  {
    case GSM_Taylor_expansion:
      U = PunctTaylorExpandAtArbitPosition(0, nvar, n1, n2, n3, v, xx, yy, zz);
      break;
    case GSM_evaluation:
      U = PunctIntPolAtArbitPositionFast(0, nvar, n1, n2, n3, cf_v, xx, yy, zz);
      break;
    default:
      assert (0);
  }
  r_plus = pow (pow (r_plus, 4) + pow (TPID::TP_epsilon, 4), 0.25);
  r_minus = pow (pow (r_minus, 4) + pow (TPID::TP_epsilon, 4), 0.25);
  if (r_plus < TPID::TP_Tiny)
      r_plus = TPID::TP_Tiny;
  if (r_minus < TPID::TP_Tiny)
      r_minus = TPID::TP_Tiny;
  CCTK_REAL psi1 = 1 + 0.5 * *mp / r_plus + 0.5 * *mm / r_minus + U;
  if (r_plus < TPID::TP_Extend_Radius) {
    psi1 = 1
       + 0.5 * EXTEND(*mp,r_plus)
       + 0.5 * *mm / r_minus + U;
  }
  if (r_minus < TPID::TP_Extend_Radius) {
    psi1 = 1
       + 0.5 * EXTEND(*mm,r_minus)
       + 0.5 * *mp / r_plus + U;
  }
  CCTK_REAL static_psi = 1;

  CCTK_REAL Aij[3][3];
  BY_Aijofxyz (xx, yy, zz, Aij);

  CCTK_REAL old_alp=1.0;
  if (TPID::multiply_old_lapse)
      old_alp = vars[VAR::U_ALPHA];

  if ((conformal_state > 0) || (pmn_lapse) || (brownsville_lapse)) {

    CCTK_REAL xp, yp, zp, rp, ir;
    CCTK_REAL s1, s3, s5;
    CCTK_REAL p, px, py, pz, pxx, pxy, pxz, pyy, pyz, pzz;
    p = 1.0;
    px = py = pz = 0.0;
    pxx = pxy = pxz = 0.0;
    pyy = pyz = pzz = 0.0;

    /* first puncture */
    xp = xx - TPID::par_b;
    yp = yy;
    zp = zz;
    rp = sqrt (xp*xp + yp*yp + zp*zp);
    rp = pow (pow (rp, 4) + pow (TPID::TP_epsilon, 4), 0.25);
    if (rp < TPID::TP_Tiny)
      rp = TPID::TP_Tiny;

    ir = 1.0/rp;

    if (rp < TPID::TP_Extend_Radius) {
      ir = EXTEND(1., rp);
    }

    s1 = 0.5* *mp *ir;
    s3 = -s1*ir*ir;
    s5 = -3.0*s3*ir*ir;

    p += s1;

    px += xp*s3;
    py += yp*s3;
    pz += zp*s3;

    pxx += xp*xp*s5 + s3;
    pxy += xp*yp*s5;
    pxz += xp*zp*s5;
    pyy += yp*yp*s5 + s3;
    pyz += yp*zp*s5;
    pzz += zp*zp*s5 + s3;

    /* second puncture */
    xp = xx + TPID::par_b;
    yp = yy;
    zp = zz;
    rp = sqrt (xp*xp + yp*yp + zp*zp);
    rp = pow (pow (rp, 4) + pow (TPID::TP_epsilon, 4), 0.25);
    if (rp < TPID::TP_Tiny)
        rp = TPID::TP_Tiny;
    ir = 1.0/rp;

    if (rp < TPID::TP_Extend_Radius) {
      ir = EXTEND(1., rp);
    }

    s1 = 0.5* *mm *ir;
    s3 = -s1*ir*ir;
    s5 = -3.0*s3*ir*ir;

    p += s1;

    px += xp*s3;
    py += yp*s3;
    pz += zp*s3;

    pxx += xp*xp*s5 + s3;
    pxy += xp*yp*s5;
    pxz += xp*zp*s5;
    pyy += yp*yp*s5 + s3;
    pyz += yp*zp*s5;
    pzz += zp*zp*s5 + s3;

    if (conformal_state >= 1) {
      static_psi = p;
    }
    if (conformal_state >= 2 && !rank) {
      printf("Code doesn't yet work for conformal_state == 2.\n");
    }
    if (conformal_state >= 3 && !rank) {
      printf("Code doesn't yet work for conformal_state == 3.\n");
    }

    if (pmn_lapse)
      vars[VAR::U_ALPHA] = pow(p, TPID::initial_lapse_psi_exponent);
    if (brownsville_lapse)
      vars[VAR::U_ALPHA] = 2.0/(1.0+pow(p, TPID::initial_lapse_psi_exponent));

  } /* if conformal-state > 0 */

  if (antisymmetric_lapse || averaged_lapse) {
    vars[VAR::U_ALPHA] =
      ((1.0 -0.5* *mp /r_plus -0.5* *mm/r_minus)
      /(1.0 +0.5* *mp /r_plus +0.5* *mm/r_minus));

    if (r_plus < TPID::TP_Extend_Radius) {
      vars[VAR::U_ALPHA] =
        ((1.0 -0.5*EXTEND(*mp, r_plus) -0.5* *mm/r_minus)
        /(1.0 +0.5*EXTEND(*mp, r_plus) +0.5* *mm/r_minus));
    }
    if (r_minus < TPID::TP_Extend_Radius) {
      vars[VAR::U_ALPHA] =
        ((1.0 -0.5*EXTEND(*mm, r_minus) -0.5* *mp/r_plus)
        /(1.0 +0.5*EXTEND(*mp, r_minus) +0.5* *mp/r_plus));
    }

    if (averaged_lapse) {
      vars[VAR::U_ALPHA] = 0.5 * (1.0 + vars[VAR::U_ALPHA]);
    }
  }

  double gd[3][3];
  gd[0][0] = pow (psi1 / static_psi, 4);
  gd[0][1] = 0.0;
  gd[0][2] = 0.0;
  gd[1][0] = 0.0;
  gd[1][1] = pow (psi1 / static_psi, 4);
  gd[1][2] = 0.0;
  gd[2][0] = 0.0;
  gd[2][1] = 0.0;
  gd[2][2] = pow (psi1 / static_psi, 4);

  double Kd[3][3];
  Kd[0][0] = Aij[0][0] / pow(psi1, 2);
  Kd[0][1] = Aij[0][1] / pow(psi1, 2);
  Kd[0][2] = Aij[0][2] / pow(psi1, 2);
  Kd[1][0] = Kd[0][1];
  Kd[1][1] = Aij[1][1] / pow(psi1, 2);
  Kd[1][2] = Aij[1][2] / pow(psi1, 2);
  Kd[2][0] = Kd[0][2];
  Kd[2][1] = Kd[1][2];
  Kd[2][2] = Aij[2][2] / pow(psi1, 2);

  double gu[3][3], gtd[3][3], Atd[3][3];
  double detgd, idetgd, trK;
  double chi;
  double t1, t2, t4, t6, t7, t9, t10, t12, t16;

  #include "adm2bssn.h"

  vars[VAR::U_SYMGT0] = gtd[0][0];
  vars[VAR::U_SYMGT1] = gtd[0][1];
  vars[VAR::U_SYMGT2] = gtd[0][2];
  vars[VAR::U_SYMGT3] = gtd[1][1];
  vars[VAR::U_SYMGT4] = gtd[1][2];
  vars[VAR::U_SYMGT5] = gtd[2][2];

  vars[VAR::U_SYMAT0] = Atd[0][0];
  vars[VAR::U_SYMAT1] = Atd[0][1];
  vars[VAR::U_SYMAT2] = Atd[0][2];
  vars[VAR::U_SYMAT3] = Atd[1][1];
  vars[VAR::U_SYMAT4] = Atd[1][2];
  vars[VAR::U_SYMAT5] = Atd[2][2];

  vars[VAR::U_K] = trK;
  vars[VAR::U_CHI] = chi;

  vars[VAR::U_BETA0] = 0.0;
  vars[VAR::U_BETA1] = 0.0;
  vars[VAR::U_BETA2] = 0.0;

  vars[VAR::U_GT0] = 0.0;
  vars[VAR::U_GT1] = 0.0;
  vars[VAR::U_GT2] = 0.0;

  vars[VAR::U_B0] = 0.0;
  vars[VAR::U_B1] = 0.0;
  vars[VAR::U_B2] = 0.0;

  if (TPID::multiply_old_lapse)
    vars[VAR::U_ALPHA] *= old_alp;

  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal - 1);
    free_derivs (&u, ntotal);
    free_derivs (&v, ntotal);
    free_derivs (&cf_v, ntotal);
  }
}

void TPStore(double *mp, double *mm, double *mp_adm, double *mm_adm, double *E, double *J1, double *J2, double *J3, const char* fprefix)
{

  * mp = TPID::par_m_plus;
  * mm = TPID::par_m_minus;
  enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
  enum GRID_SETUP_METHOD gsm;
  int antisymmetric_lapse, averaged_lapse, pmn_lapse, brownsville_lapse;
  int const nvar = 1, n1 = TPID::npoints_A, n2 = TPID::npoints_B, n3 = TPID::npoints_phi;
  int const ntotal = n1 * n2 * n3 * nvar;

  #if 1
    int percent10 = 0;
  #endif
  
  static CCTK_REAL *F = NULL;
  static derivs u, v, cf_v;
  CCTK_REAL admMass;

  if (!F) {
    CCTK_REAL up, um;
    /* Solve only when called for the first time */
    F = dvector (0, ntotal - 1);
    allocate_derivs (&u, ntotal);
    allocate_derivs (&v, ntotal);
    allocate_derivs (&cf_v, ntotal);

    printf("INFO: b = %g\n", TPID::par_b);

    /* initialise to 0 */
    #pragma omp parallel for num_threads(TP_OMP_THREADS)
    for (int j = 0; j < ntotal; j++)
    {
      cf_v.d0[j] = 0.0;
      cf_v.d1[j] = 0.0;
      cf_v.d2[j] = 0.0;
      cf_v.d3[j] = 0.0;
      cf_v.d11[j] = 0.0;
      cf_v.d12[j] = 0.0;
      cf_v.d13[j] = 0.0;
      cf_v.d22[j] = 0.0;
      cf_v.d23[j] = 0.0;
      cf_v.d33[j] = 0.0;
      v.d0[j] = 0.0;
      v.d1[j] = 0.0;
      v.d2[j] = 0.0;
      v.d3[j] = 0.0;
      v.d11[j] = 0.0;
      v.d12[j] = 0.0;
      v.d13[j] = 0.0;
      v.d22[j] = 0.0;
      v.d23[j] = 0.0;
      v.d33[j] = 0.0;
    }
    
    /* If bare masses are not given, iteratively solve for them given the
       target ADM masses target_M_plus and target_M_minus and with initial
       guesses given by par_m_plus and par_m_minus. */
    if(!(TPID::give_bare_mass)) {
      CCTK_REAL tmp, mp_adm_err, mm_adm_err;
      char valbuf[100];

      CCTK_REAL M_p = TPID::target_M_plus;
      CCTK_REAL M_m = TPID::target_M_minus;

      printf("Attempting to find bare masses.\n");
      printf("Target ADM masses: M_p=%g and M_m=%g\n",
                  (double) M_p, (double) M_m);
      printf("ADM mass tolerance: %g\n", (double) TPID::adm_tol);

      /* Loop until both ADM masses are within adm_tol of their target */
      do {
        printf("Bare masses: mp=%.15g, mm=%.15g\n",
                    (double)*mp, (double)*mm);
        Newton (nvar, n1, n2, n3, v, TPID::Newton_tol, 1);

        F_of_v (nvar, n1, n2, n3, v, F, u);

        up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, TPID::par_b, 0., 0.);
        um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-TPID::par_b, 0., 0.);

        /* Calculate the ADM masses from the current bare mass guess */
        *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * TPID::par_b);
        *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * TPID::par_b);

        /* Check how far the current ADM masses are from the target */
        mp_adm_err = fabs(M_p-*mp_adm);
        mm_adm_err = fabs(M_m-*mm_adm);
        printf("ADM mass error: M_p_err=%.15g, M_m_err=%.15g\n",
                    mp_adm_err, mm_adm_err);

        /* Invert the ADM mass equation and update the bare mass guess so that
           it gives the correct target ADM masses */
        tmp = -4*TPID::par_b*( 1 + um + up + um*up ) +
                sqrt(16*TPID::par_b*M_m*(1 + um)*(1 + up) +
                  pow(-M_m + M_p + 4*TPID::par_b*(1 + um)*(1 + up),2));
        *mp = (tmp + M_p - M_m)/(2.*(1 + up));
        *mm = (tmp - M_p + M_m)/(2.*(1 + um));

        /* Set the par_m_plus and par_m_minus parameters */
        /*
        sprintf (valbuf,"%.17g", (double) *mp);
        CCTK_ParameterSet ("par_m_plus", "TwoPunctures", valbuf);
        sprintf (valbuf,"%.17g", (double) *mm);
        CCTK_ParameterSet ("par_m_minus", "TwoPunctures", valbuf);
        */

        TPID::par_m_plus = *mp;
        TPID::par_m_minus = *mm;

      } while ( (mp_adm_err > TPID::adm_tol) ||
                (mm_adm_err > TPID::adm_tol) );

      printf("Found bare masses.");
    }

    Newton (nvar, n1, n2, n3, v, TPID::Newton_tol, TPID::Newton_maxit);

    F_of_v (nvar, n1, n2, n3, v, F, u);

    SpecCoef(n1, n2, n3, 0, v.d0, cf_v.d0);

    printf("The two puncture masses are mp=%.17g and mm=%.17g\n", *mp, *mm);

    up = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v, TPID::par_b, 0., 0.);
    um = PunctIntPolAtArbitPosition(0, nvar, n1, n2, n3, v,-TPID::par_b, 0., 0.);

    /* Calculate the ADM masses from the current bare mass guess */
    *mp_adm = (1 + up) * *mp + *mp * *mm / (4. * TPID::par_b);
    *mm_adm = (1 + um) * *mm + *mp * *mm / (4. * TPID::par_b);

    printf("Puncture 1 ADM mass is %g\n", *mp_adm);
    printf("Puncture 2 ADM mass is %g\n", *mm_adm);

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
    admMass = (*mp + *mm
               - 4*TPID::par_b*PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, nvar, n1, n2, n3));
    printf("The total ADM mass is %g\n", admMass);
    *E = admMass;

    /*
      Run this in Mathematica (version 8 or later) with
        math -script <file>

      Needs["SymbolicC`"];
      co = Table["center_offset[" <> ToString[i] <> "]", {i, 0, 2}];
      r1 = co + {"par_b", 0, 0};
      r2 = co + {-"par_b", 0, 0};
      {p1, p2} = Table["par_P_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];
      {s1, s2} = Table["par_S_" <> bh <> "[" <> ToString[i] <> "]", {bh, {"plus", "minus"}}, {i, 0, 2}];

      J = Cross[r1, p1] + Cross[r2, p2] + s1 + s2;

      JVar = Table["*J" <> ToString[i], {i, 1, 3}];
      Print[OutputForm@StringReplace[
        ToCCodeString@MapThread[CAssign[#1, CExpression[#2]] &, {JVar, J}],
        "\"" -> ""]];
     */

    *J1 = -(TPID::center_offset[2]*TPID::par_P_minus[1]) + TPID::center_offset[1]*TPID::par_P_minus[2] - TPID::center_offset[2]*TPID::par_P_plus[1] + TPID::center_offset[1]*TPID::par_P_plus[2] + TPID::par_S_minus[0] + TPID::par_S_plus[0];
    *J2 = TPID::center_offset[2]*TPID::par_P_minus[0] - TPID::center_offset[0]*TPID::par_P_minus[2] + TPID::par_b*TPID::par_P_minus[2] + TPID::center_offset[2]*TPID::par_P_plus[0] - TPID::center_offset[0]*TPID::par_P_plus[2] - TPID::par_b*TPID::par_P_plus[2] + TPID::par_S_minus[1] + TPID::par_S_plus[1];
    *J3 = -(TPID::center_offset[1]*TPID::par_P_minus[0]) + TPID::center_offset[0]*TPID::par_P_minus[1] - TPID::par_b*TPID::par_P_minus[1] - TPID::center_offset[1]*TPID::par_P_plus[0] + TPID::center_offset[0]*TPID::par_P_plus[1] + TPID::par_b*TPID::par_P_plus[1] + TPID::par_S_minus[2] + TPID::par_S_plus[2];
  }

  char fName[300];
  sprintf(fName,"%s_tpid_sol.bin",fprefix);
  FILE *write_ptr;
  write_ptr = fopen(fName,"wb");  // w for write, b for binary
  if(write_ptr==NULL)
  {
    printf("tpid solver data write failed\n");
    return;
  }
    
  fwrite(&ntotal,sizeof(int),1,write_ptr); 
  fwrite(F,sizeof(double),ntotal-1,write_ptr); 
  
  fwrite(u.d0,sizeof(double),ntotal,write_ptr); 
  fwrite(u.d1,sizeof(double),ntotal,write_ptr);
  fwrite(u.d2,sizeof(double),ntotal,write_ptr);
  fwrite(u.d3,sizeof(double),ntotal,write_ptr);
  fwrite(u.d11,sizeof(double),ntotal,write_ptr);
  fwrite(u.d12,sizeof(double),ntotal,write_ptr);
  fwrite(u.d13,sizeof(double),ntotal,write_ptr);
  fwrite(u.d22,sizeof(double),ntotal,write_ptr);
  fwrite(u.d23,sizeof(double),ntotal,write_ptr);
  fwrite(u.d33,sizeof(double),ntotal,write_ptr);


  fwrite(v.d0,sizeof(double),ntotal,write_ptr); 
  fwrite(v.d1,sizeof(double),ntotal,write_ptr);
  fwrite(v.d2,sizeof(double),ntotal,write_ptr);
  fwrite(v.d3,sizeof(double),ntotal,write_ptr);
  fwrite(v.d11,sizeof(double),ntotal,write_ptr);
  fwrite(v.d12,sizeof(double),ntotal,write_ptr);
  fwrite(v.d13,sizeof(double),ntotal,write_ptr);
  fwrite(v.d22,sizeof(double),ntotal,write_ptr);
  fwrite(v.d23,sizeof(double),ntotal,write_ptr);
  fwrite(v.d33,sizeof(double),ntotal,write_ptr);


  fwrite(cf_v.d0,sizeof(double),ntotal,write_ptr); 
  fwrite(cf_v.d1,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d2,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d3,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d11,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d12,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d13,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d22,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d23,sizeof(double),ntotal,write_ptr);
  fwrite(cf_v.d33,sizeof(double),ntotal,write_ptr);
  fclose(write_ptr);

  // for(unsigned int w=0; w< ntotal-1; w++)
  //   printf("f[%d]: %f \n",w,F[w]);

  // for(unsigned int w=0; w< ntotal; w++)
  //   printf("store: u.d0[%d]: %f \n",w,u.d0[w]);

  free_dvector (F, 0, ntotal - 1);
  free_derivs (&u, ntotal);
  free_derivs (&v, ntotal);
  free_derivs (&cf_v, ntotal);

  return;


}

void TPRestore(CCTK_REAL*& F, derivs& u, derivs& v, derivs& cf_v, const char* fprefix,bool mpi_bcast)
{

  int rank,npes;
  MPI_Comm_rank(TP_MPI_COMM,&rank);
  MPI_Comm_size(TP_MPI_COMM,&npes);

  char fName[300];
  sprintf(fName,"%s_tpid_sol.bin",fprefix);

  int ntotal;
  FILE *read_ptr;
  size_t fr_st=0;

  read_ptr = fopen(fName,"rb");  // w for write, b for binary
  if(read_ptr==NULL)
  {
    printf("tpid solver data read failed\n");
    return;
  }

  MPI_Barrier(TP_MPI_COMM);
  if(!rank)
   std::cout<<"TPID data read begin: "<<std::endl;

  fr_st = fread(&ntotal,sizeof(int),1,read_ptr); 

  F = dvector (0, ntotal - 1);
  allocate_derivs (&u, ntotal);
  allocate_derivs (&v, ntotal);
  allocate_derivs (&cf_v, ntotal);
  
  fr_st=fread(F,sizeof(double),ntotal-1,read_ptr); 
  fr_st=fread(u.d0,sizeof(double),ntotal,read_ptr); 
  fr_st=fread(u.d1,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d2,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d3,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d11,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d12,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d13,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d22,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d23,sizeof(double),ntotal,read_ptr);
  fr_st=fread(u.d33,sizeof(double),ntotal,read_ptr);

  fr_st=fread(v.d0,sizeof(double),ntotal,read_ptr); 
  fr_st=fread(v.d1,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d2,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d3,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d11,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d12,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d13,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d22,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d23,sizeof(double),ntotal,read_ptr);
  fr_st=fread(v.d33,sizeof(double),ntotal,read_ptr);

  fr_st=fread(cf_v.d0,sizeof(double),ntotal,read_ptr); 
  fr_st=fread(cf_v.d1,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d2,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d3,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d11,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d12,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d13,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d22,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d23,sizeof(double),ntotal,read_ptr);
  fr_st=fread(cf_v.d33,sizeof(double),ntotal,read_ptr);
  
  MPI_Barrier(TP_MPI_COMM);
  if(!rank)
   std::cout<<"TPID solver successfully restored: "<<std::endl;

  return;


}


