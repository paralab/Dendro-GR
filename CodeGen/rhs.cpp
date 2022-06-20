#include "rhs.h"
#include "gr.h"

using namespace std;
using namespace bssn;

/*----------------------------------------------------------------------;
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void bssnrhs(double **unzipVarsRHS, const double **uZipVars,
             const unsigned int& offset,
             const double *pmin, const double *pmax, const unsigned int *sz,
             const unsigned int& bflag)
{



  const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
  const double *chi = &uZipVars[VAR::U_CHI][offset];
  const double *K = &uZipVars[VAR::U_K][offset];
  const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
  const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
  const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
  const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
  const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
  const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
  const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
  const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
  const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
  const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
  const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
  const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
  const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
  const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
  const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
  const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
  const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
  const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
  const double *B0 = &uZipVars[VAR::U_B0][offset];
  const double *B1 = &uZipVars[VAR::U_B1][offset];
  const double *B2 = &uZipVars[VAR::U_B2][offset];

  double *a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
  double *chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
  double *K_rhs = &unzipVarsRHS[VAR::U_K][offset];
  double *gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
  double *gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
  double *gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
  double *gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
  double *gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
  double *gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
  double *b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
  double *b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
  double *b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
  double *At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
  double *At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
  double *At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
  double *At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
  double *At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
  double *At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
  double *Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
  double *Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
  double *Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
  double *B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
  double *B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
  double *B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

  const unsigned int nx = sz[0];
  const unsigned int ny = sz[1];
  const unsigned int nz = sz[2];

  double hx = (pmax[0] - pmin[0]) / (nx - 1);
  double hy = (pmax[1] - pmin[1]) / (ny - 1);
  double hz = (pmax[2] - pmin[2]) / (nz - 1);

  const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                  BSSN_LAMBDA[2], BSSN_LAMBDA[3]};
  const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};



  int idx[3];

  unsigned int n = sz[0]*sz[1]*sz[2];

#if 0
  double vars[24];
  pmin[0] = 0.0; pmin[1] = 0.0; pmin[2] = 0.0;
  pmax[0] = 2.0; pmax[1] = 2.0; pmax[2] = 2.0;

  hx = (pmax[0] - pmin[0]) / (nx - 1);
  hy = (pmax[1] - pmin[1]) / (ny - 1);
  hz = (pmax[2] - pmin[2]) / (nz - 1);

  for (unsigned int k = 0; k < nz; k++) {
    double z = pmin[2] + k*hz;
    for (unsigned int j = 0; j < ny; j++) {
      double y = pmin[1] + j*hy;
      for (unsigned int i = 0; i < nx; i++) {
        double x = pmin[0] + i*hx;
        int pp = i + nx*(j + k*ny);
        fake_initial_data(x, y, z, vars);
        for (unsigned int m = 0; m < 24; m++) {
          uZipVars[m][offset+pp] = vars[m];
        }
      }
    }
  }
#endif

bssn::timer::t_deriv.start();


/*[[[cog
import cog
import bssnDerivs as bssnDerivs

for deriv in bssnDerivs.FUNC_D_I:
    cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")
    
for deriv in bssnDerivs.FUNC_D_IJ:
    cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

for deriv in bssnDerivs.FUNC_AD_I:
    cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")   
    
    
for var in bssnDerivs.D:
    cog.outl("\t deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
    cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
    cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))
    
    if var in bssnDerivs.DD:
        cog.outl("\t deriv_xx(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[0] + var ,var))
        cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[1] + var , bssnDerivs.PREFIX_D[0] + var ))
        cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[2] + var , bssnDerivs.PREFIX_D[0] + var ))
        
        cog.outl("\t deriv_yy(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[3] + var ,var))
        cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[4] + var , bssnDerivs.PREFIX_D[1] + var))
        
        cog.outl("\t deriv_zz(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[5] + var ,var))
        
    if var in bssnDerivs.AD:
        cog.outl("\t adv_deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_AD[0] + var ,var))
        cog.outl("\t adv_deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_AD[1] + var ,var))
        cog.outl("\t adv_deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_AD[2] + var ,var))
        
    
    
]]]*/
double* grad_0_alpha = (double*)malloc(sizeof(double)*n);
double* grad_1_alpha = (double*)malloc(sizeof(double)*n);
double* grad_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad_0_beta0 = (double*)malloc(sizeof(double)*n);
double* grad_1_beta0 = (double*)malloc(sizeof(double)*n);
double* grad_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad_0_beta1 = (double*)malloc(sizeof(double)*n);
double* grad_1_beta1 = (double*)malloc(sizeof(double)*n);
double* grad_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad_0_beta2 = (double*)malloc(sizeof(double)*n);
double* grad_1_beta2 = (double*)malloc(sizeof(double)*n);
double* grad_2_beta2 = (double*)malloc(sizeof(double)*n);
double* grad_0_B0 = (double*)malloc(sizeof(double)*n);
double* grad_1_B0 = (double*)malloc(sizeof(double)*n);
double* grad_2_B0 = (double*)malloc(sizeof(double)*n);
double* grad_0_B1 = (double*)malloc(sizeof(double)*n);
double* grad_1_B1 = (double*)malloc(sizeof(double)*n);
double* grad_2_B1 = (double*)malloc(sizeof(double)*n);
double* grad_0_B2 = (double*)malloc(sizeof(double)*n);
double* grad_1_B2 = (double*)malloc(sizeof(double)*n);
double* grad_2_B2 = (double*)malloc(sizeof(double)*n);
double* grad_0_chi = (double*)malloc(sizeof(double)*n);
double* grad_1_chi = (double*)malloc(sizeof(double)*n);
double* grad_2_chi = (double*)malloc(sizeof(double)*n);
double* grad_0_Gt0 = (double*)malloc(sizeof(double)*n);
double* grad_1_Gt0 = (double*)malloc(sizeof(double)*n);
double* grad_2_Gt0 = (double*)malloc(sizeof(double)*n);
double* grad_0_Gt1 = (double*)malloc(sizeof(double)*n);
double* grad_1_Gt1 = (double*)malloc(sizeof(double)*n);
double* grad_2_Gt1 = (double*)malloc(sizeof(double)*n);
double* grad_0_Gt2 = (double*)malloc(sizeof(double)*n);
double* grad_1_Gt2 = (double*)malloc(sizeof(double)*n);
double* grad_2_Gt2 = (double*)malloc(sizeof(double)*n);
double* grad_0_K = (double*)malloc(sizeof(double)*n);
double* grad_1_K = (double*)malloc(sizeof(double)*n);
double* grad_2_K = (double*)malloc(sizeof(double)*n);
double* grad_0_gt0 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt0 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt1 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt1 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt2 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt2 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt3 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt3 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt4 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt4 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt5 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt5 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad_0_At0 = (double*)malloc(sizeof(double)*n);
double* grad_1_At0 = (double*)malloc(sizeof(double)*n);
double* grad_2_At0 = (double*)malloc(sizeof(double)*n);
double* grad_0_At1 = (double*)malloc(sizeof(double)*n);
double* grad_1_At1 = (double*)malloc(sizeof(double)*n);
double* grad_2_At1 = (double*)malloc(sizeof(double)*n);
double* grad_0_At2 = (double*)malloc(sizeof(double)*n);
double* grad_1_At2 = (double*)malloc(sizeof(double)*n);
double* grad_2_At2 = (double*)malloc(sizeof(double)*n);
double* grad_0_At3 = (double*)malloc(sizeof(double)*n);
double* grad_1_At3 = (double*)malloc(sizeof(double)*n);
double* grad_2_At3 = (double*)malloc(sizeof(double)*n);
double* grad_0_At4 = (double*)malloc(sizeof(double)*n);
double* grad_1_At4 = (double*)malloc(sizeof(double)*n);
double* grad_2_At4 = (double*)malloc(sizeof(double)*n);
double* grad_0_At5 = (double*)malloc(sizeof(double)*n);
double* grad_1_At5 = (double*)malloc(sizeof(double)*n);
double* grad_2_At5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_chi = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_chi = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_chi = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_chi = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_chi = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_chi = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt3 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt3 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt3 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt4 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt4 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt4 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt5 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt5 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt5 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At3 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At3 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At3 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At4 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At4 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At4 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At5 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At5 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At5 = (double*)malloc(sizeof(double)*n);
double* agrad_0_alpha = (double*)malloc(sizeof(double)*n);
double* agrad_1_alpha = (double*)malloc(sizeof(double)*n);
double* agrad_2_alpha = (double*)malloc(sizeof(double)*n);
double* agrad_0_beta0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_beta0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_beta0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_beta1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_beta1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_beta1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_chi = (double*)malloc(sizeof(double)*n);
double* agrad_1_chi = (double*)malloc(sizeof(double)*n);
double* agrad_2_chi = (double*)malloc(sizeof(double)*n);
double* agrad_0_Gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_Gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_Gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_Gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_Gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_Gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_Gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_Gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_Gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_K = (double*)malloc(sizeof(double)*n);
double* agrad_1_K = (double*)malloc(sizeof(double)*n);
double* agrad_2_K = (double*)malloc(sizeof(double)*n);
double* agrad_0_B0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_B0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_B0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_B1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_B1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_B1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_B2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_B2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_B2 = (double*)malloc(sizeof(double)*n);
deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
deriv_xx(grad2_0_0_alpha, alpha, hx, sz, bflag);
deriv_y(grad2_0_1_alpha, grad_0_alpha, hx, sz, bflag);
deriv_z(grad2_0_2_alpha, grad_0_alpha, hx, sz, bflag);
deriv_yy(grad2_1_1_alpha, alpha, hx, sz, bflag);
deriv_z(grad2_1_2_alpha, grad_1_alpha, hx, sz, bflag);
deriv_zz(grad2_2_2_alpha, alpha, hx, sz, bflag);
adv_deriv_x(agrad_0_alpha, alpha, hx, sz, bflag);
adv_deriv_y(agrad_1_alpha, alpha, hx, sz, bflag);
adv_deriv_z(agrad_2_alpha, alpha, hx, sz, bflag);
deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
deriv_xx(grad2_0_0_beta0, beta0, hx, sz, bflag);
deriv_y(grad2_0_1_beta0, grad_0_beta0, hx, sz, bflag);
deriv_z(grad2_0_2_beta0, grad_0_beta0, hx, sz, bflag);
deriv_yy(grad2_1_1_beta0, beta0, hx, sz, bflag);
deriv_z(grad2_1_2_beta0, grad_1_beta0, hx, sz, bflag);
deriv_zz(grad2_2_2_beta0, beta0, hx, sz, bflag);
adv_deriv_x(agrad_0_beta0, beta0, hx, sz, bflag);
adv_deriv_y(agrad_1_beta0, beta0, hx, sz, bflag);
adv_deriv_z(agrad_2_beta0, beta0, hx, sz, bflag);
deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
deriv_xx(grad2_0_0_beta1, beta1, hx, sz, bflag);
deriv_y(grad2_0_1_beta1, grad_0_beta1, hx, sz, bflag);
deriv_z(grad2_0_2_beta1, grad_0_beta1, hx, sz, bflag);
deriv_yy(grad2_1_1_beta1, beta1, hx, sz, bflag);
deriv_z(grad2_1_2_beta1, grad_1_beta1, hx, sz, bflag);
deriv_zz(grad2_2_2_beta1, beta1, hx, sz, bflag);
adv_deriv_x(agrad_0_beta1, beta1, hx, sz, bflag);
adv_deriv_y(agrad_1_beta1, beta1, hx, sz, bflag);
adv_deriv_z(agrad_2_beta1, beta1, hx, sz, bflag);
deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
deriv_xx(grad2_0_0_beta2, beta2, hx, sz, bflag);
deriv_y(grad2_0_1_beta2, grad_0_beta2, hx, sz, bflag);
deriv_z(grad2_0_2_beta2, grad_0_beta2, hx, sz, bflag);
deriv_yy(grad2_1_1_beta2, beta2, hx, sz, bflag);
deriv_z(grad2_1_2_beta2, grad_1_beta2, hx, sz, bflag);
deriv_zz(grad2_2_2_beta2, beta2, hx, sz, bflag);
adv_deriv_x(agrad_0_beta2, beta2, hx, sz, bflag);
adv_deriv_y(agrad_1_beta2, beta2, hx, sz, bflag);
adv_deriv_z(agrad_2_beta2, beta2, hx, sz, bflag);
deriv_x(grad_0_B0, B0, hx, sz, bflag);
deriv_y(grad_1_B0, B0, hx, sz, bflag);
deriv_z(grad_2_B0, B0, hx, sz, bflag);
adv_deriv_x(agrad_0_B0, B0, hx, sz, bflag);
adv_deriv_y(agrad_1_B0, B0, hx, sz, bflag);
adv_deriv_z(agrad_2_B0, B0, hx, sz, bflag);
deriv_x(grad_0_B1, B1, hx, sz, bflag);
deriv_y(grad_1_B1, B1, hx, sz, bflag);
deriv_z(grad_2_B1, B1, hx, sz, bflag);
adv_deriv_x(agrad_0_B1, B1, hx, sz, bflag);
adv_deriv_y(agrad_1_B1, B1, hx, sz, bflag);
adv_deriv_z(agrad_2_B1, B1, hx, sz, bflag);
deriv_x(grad_0_B2, B2, hx, sz, bflag);
deriv_y(grad_1_B2, B2, hx, sz, bflag);
deriv_z(grad_2_B2, B2, hx, sz, bflag);
adv_deriv_x(agrad_0_B2, B2, hx, sz, bflag);
adv_deriv_y(agrad_1_B2, B2, hx, sz, bflag);
adv_deriv_z(agrad_2_B2, B2, hx, sz, bflag);
deriv_x(grad_0_chi, chi, hx, sz, bflag);
deriv_y(grad_1_chi, chi, hx, sz, bflag);
deriv_z(grad_2_chi, chi, hx, sz, bflag);
deriv_xx(grad2_0_0_chi, chi, hx, sz, bflag);
deriv_y(grad2_0_1_chi, grad_0_chi, hx, sz, bflag);
deriv_z(grad2_0_2_chi, grad_0_chi, hx, sz, bflag);
deriv_yy(grad2_1_1_chi, chi, hx, sz, bflag);
deriv_z(grad2_1_2_chi, grad_1_chi, hx, sz, bflag);
deriv_zz(grad2_2_2_chi, chi, hx, sz, bflag);
adv_deriv_x(agrad_0_chi, chi, hx, sz, bflag);
adv_deriv_y(agrad_1_chi, chi, hx, sz, bflag);
adv_deriv_z(agrad_2_chi, chi, hx, sz, bflag);
deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
adv_deriv_x(agrad_0_Gt0, Gt0, hx, sz, bflag);
adv_deriv_y(agrad_1_Gt0, Gt0, hx, sz, bflag);
adv_deriv_z(agrad_2_Gt0, Gt0, hx, sz, bflag);
deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
adv_deriv_x(agrad_0_Gt1, Gt1, hx, sz, bflag);
adv_deriv_y(agrad_1_Gt1, Gt1, hx, sz, bflag);
adv_deriv_z(agrad_2_Gt1, Gt1, hx, sz, bflag);
deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
adv_deriv_x(agrad_0_Gt2, Gt2, hx, sz, bflag);
adv_deriv_y(agrad_1_Gt2, Gt2, hx, sz, bflag);
adv_deriv_z(agrad_2_Gt2, Gt2, hx, sz, bflag);
deriv_x(grad_0_K, K, hx, sz, bflag);
deriv_y(grad_1_K, K, hx, sz, bflag);
deriv_z(grad_2_K, K, hx, sz, bflag);
adv_deriv_x(agrad_0_K, K, hx, sz, bflag);
adv_deriv_y(agrad_1_K, K, hx, sz, bflag);
adv_deriv_z(agrad_2_K, K, hx, sz, bflag);
deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
deriv_xx(grad2_0_0_gt0, gt0, hx, sz, bflag);
deriv_y(grad2_0_1_gt0, grad_0_gt0, hx, sz, bflag);
deriv_z(grad2_0_2_gt0, grad_0_gt0, hx, sz, bflag);
deriv_yy(grad2_1_1_gt0, gt0, hx, sz, bflag);
deriv_z(grad2_1_2_gt0, grad_1_gt0, hx, sz, bflag);
deriv_zz(grad2_2_2_gt0, gt0, hx, sz, bflag);
adv_deriv_x(agrad_0_gt0, gt0, hx, sz, bflag);
adv_deriv_y(agrad_1_gt0, gt0, hx, sz, bflag);
adv_deriv_z(agrad_2_gt0, gt0, hx, sz, bflag);
deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
deriv_xx(grad2_0_0_gt1, gt1, hx, sz, bflag);
deriv_y(grad2_0_1_gt1, grad_0_gt1, hx, sz, bflag);
deriv_z(grad2_0_2_gt1, grad_0_gt1, hx, sz, bflag);
deriv_yy(grad2_1_1_gt1, gt1, hx, sz, bflag);
deriv_z(grad2_1_2_gt1, grad_1_gt1, hx, sz, bflag);
deriv_zz(grad2_2_2_gt1, gt1, hx, sz, bflag);
adv_deriv_x(agrad_0_gt1, gt1, hx, sz, bflag);
adv_deriv_y(agrad_1_gt1, gt1, hx, sz, bflag);
adv_deriv_z(agrad_2_gt1, gt1, hx, sz, bflag);
deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
deriv_xx(grad2_0_0_gt2, gt2, hx, sz, bflag);
deriv_y(grad2_0_1_gt2, grad_0_gt2, hx, sz, bflag);
deriv_z(grad2_0_2_gt2, grad_0_gt2, hx, sz, bflag);
deriv_yy(grad2_1_1_gt2, gt2, hx, sz, bflag);
deriv_z(grad2_1_2_gt2, grad_1_gt2, hx, sz, bflag);
deriv_zz(grad2_2_2_gt2, gt2, hx, sz, bflag);
adv_deriv_x(agrad_0_gt2, gt2, hx, sz, bflag);
adv_deriv_y(agrad_1_gt2, gt2, hx, sz, bflag);
adv_deriv_z(agrad_2_gt2, gt2, hx, sz, bflag);
deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
deriv_xx(grad2_0_0_gt3, gt3, hx, sz, bflag);
deriv_y(grad2_0_1_gt3, grad_0_gt3, hx, sz, bflag);
deriv_z(grad2_0_2_gt3, grad_0_gt3, hx, sz, bflag);
deriv_yy(grad2_1_1_gt3, gt3, hx, sz, bflag);
deriv_z(grad2_1_2_gt3, grad_1_gt3, hx, sz, bflag);
deriv_zz(grad2_2_2_gt3, gt3, hx, sz, bflag);
adv_deriv_x(agrad_0_gt3, gt3, hx, sz, bflag);
adv_deriv_y(agrad_1_gt3, gt3, hx, sz, bflag);
adv_deriv_z(agrad_2_gt3, gt3, hx, sz, bflag);
deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
deriv_xx(grad2_0_0_gt4, gt4, hx, sz, bflag);
deriv_y(grad2_0_1_gt4, grad_0_gt4, hx, sz, bflag);
deriv_z(grad2_0_2_gt4, grad_0_gt4, hx, sz, bflag);
deriv_yy(grad2_1_1_gt4, gt4, hx, sz, bflag);
deriv_z(grad2_1_2_gt4, grad_1_gt4, hx, sz, bflag);
deriv_zz(grad2_2_2_gt4, gt4, hx, sz, bflag);
adv_deriv_x(agrad_0_gt4, gt4, hx, sz, bflag);
adv_deriv_y(agrad_1_gt4, gt4, hx, sz, bflag);
adv_deriv_z(agrad_2_gt4, gt4, hx, sz, bflag);
deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
deriv_xx(grad2_0_0_gt5, gt5, hx, sz, bflag);
deriv_y(grad2_0_1_gt5, grad_0_gt5, hx, sz, bflag);
deriv_z(grad2_0_2_gt5, grad_0_gt5, hx, sz, bflag);
deriv_yy(grad2_1_1_gt5, gt5, hx, sz, bflag);
deriv_z(grad2_1_2_gt5, grad_1_gt5, hx, sz, bflag);
deriv_zz(grad2_2_2_gt5, gt5, hx, sz, bflag);
adv_deriv_x(agrad_0_gt5, gt5, hx, sz, bflag);
adv_deriv_y(agrad_1_gt5, gt5, hx, sz, bflag);
adv_deriv_z(agrad_2_gt5, gt5, hx, sz, bflag);
deriv_x(grad_0_At0, At0, hx, sz, bflag);
deriv_y(grad_1_At0, At0, hx, sz, bflag);
deriv_z(grad_2_At0, At0, hx, sz, bflag);
adv_deriv_x(agrad_0_At0, At0, hx, sz, bflag);
adv_deriv_y(agrad_1_At0, At0, hx, sz, bflag);
adv_deriv_z(agrad_2_At0, At0, hx, sz, bflag);
deriv_x(grad_0_At1, At1, hx, sz, bflag);
deriv_y(grad_1_At1, At1, hx, sz, bflag);
deriv_z(grad_2_At1, At1, hx, sz, bflag);
adv_deriv_x(agrad_0_At1, At1, hx, sz, bflag);
adv_deriv_y(agrad_1_At1, At1, hx, sz, bflag);
adv_deriv_z(agrad_2_At1, At1, hx, sz, bflag);
deriv_x(grad_0_At2, At2, hx, sz, bflag);
deriv_y(grad_1_At2, At2, hx, sz, bflag);
deriv_z(grad_2_At2, At2, hx, sz, bflag);
adv_deriv_x(agrad_0_At2, At2, hx, sz, bflag);
adv_deriv_y(agrad_1_At2, At2, hx, sz, bflag);
adv_deriv_z(agrad_2_At2, At2, hx, sz, bflag);
deriv_x(grad_0_At3, At3, hx, sz, bflag);
deriv_y(grad_1_At3, At3, hx, sz, bflag);
deriv_z(grad_2_At3, At3, hx, sz, bflag);
adv_deriv_x(agrad_0_At3, At3, hx, sz, bflag);
adv_deriv_y(agrad_1_At3, At3, hx, sz, bflag);
adv_deriv_z(agrad_2_At3, At3, hx, sz, bflag);
deriv_x(grad_0_At4, At4, hx, sz, bflag);
deriv_y(grad_1_At4, At4, hx, sz, bflag);
deriv_z(grad_2_At4, At4, hx, sz, bflag);
adv_deriv_x(agrad_0_At4, At4, hx, sz, bflag);
adv_deriv_y(agrad_1_At4, At4, hx, sz, bflag);
adv_deriv_z(agrad_2_At4, At4, hx, sz, bflag);
deriv_x(grad_0_At5, At5, hx, sz, bflag);
deriv_y(grad_1_At5, At5, hx, sz, bflag);
deriv_z(grad_2_At5, At5, hx, sz, bflag);
adv_deriv_x(agrad_0_At5, At5, hx, sz, bflag);
adv_deriv_y(agrad_1_At5, At5, hx, sz, bflag);
adv_deriv_z(agrad_2_At5, At5, hx, sz, bflag);
//[[[end]]]

bssn::timer::t_deriv.stop();

  register double x;
  register double y;
  register double z;
  register unsigned int pp;

  double r_coord;
  double eta;

  //cout << "begin loop" << endl;
  for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }


bssn::timer::t_rhs.start();

#ifdef USE_ETA_FUNC
/*[[[cog

import dendro
import bssn

outs = [bssn.a_rhs, bssn.b_rhs, bssn.gt_rhs, bssn.chi_rhs, bssn.At_rhs, bssn.K_rhs, bssn.Gt_rhs, bssn.B_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
#outs=outs[0:3]
#vnames=vnames[0:3]
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 668203 
// Dendro: printing temp variables
double DENDRO_0 = 2*alpha[pp];
double DENDRO_1 = (3.0/4.0)*alpha[pp]*lambda_f[1] + (3.0/4.0)*lambda_f[0];
double DENDRO_2 = grad_0_beta0[pp];
double DENDRO_3 = (4.0/3.0)*DENDRO_2;
double DENDRO_4 = grad_1_beta1[pp];
double DENDRO_5 = (2.0/3.0)*gt0[pp];
double DENDRO_6 = grad_2_beta2[pp];
double DENDRO_7 = grad_0_beta1[pp];
double DENDRO_8 = 2*gt1[pp];
double DENDRO_9 = grad_0_beta2[pp];
double DENDRO_10 = 2*gt2[pp];
double DENDRO_11 = grad_1_beta0[pp];
double DENDRO_12 = grad_1_beta2[pp];
double DENDRO_13 = (1.0/3.0)*gt1[pp];
double DENDRO_14 = (2.0/3.0)*DENDRO_6;
double DENDRO_15 = grad_2_beta0[pp];
double DENDRO_16 = grad_2_beta1[pp];
double DENDRO_17 = (1.0/3.0)*gt2[pp];
double DENDRO_18 = (2.0/3.0)*DENDRO_4;
double DENDRO_19 = (2.0/3.0)*DENDRO_2;
double DENDRO_20 = (4.0/3.0)*DENDRO_4;
double DENDRO_21 = 2*gt4[pp];
double DENDRO_22 = (1.0/3.0)*gt4[pp];
double DENDRO_23 = (4.0/3.0)*DENDRO_6;
double DENDRO_24 = (2.0/3.0)*chi[pp];
double DENDRO_25 = DENDRO_2 + DENDRO_4 + DENDRO_6;
double DENDRO_26 = 2*At1[pp];
double DENDRO_27 = 2*At2[pp];
double DENDRO_28 = gt2[pp]*gt4[pp];
double DENDRO_29 = -DENDRO_28 + gt1[pp]*gt5[pp];
double DENDRO_30 = At0[pp]*DENDRO_29;
double DENDRO_31 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
double DENDRO_32 = At2[pp]*DENDRO_31;
double DENDRO_33 = gt0[pp]*gt5[pp];
double DENDRO_34 = pow(gt2[pp], 2);
double DENDRO_35 = DENDRO_33 - DENDRO_34;
double DENDRO_36 = At1[pp]*DENDRO_35;
double DENDRO_37 = pow(gt4[pp], 2);
double DENDRO_38 = pow(gt1[pp], 2);
double DENDRO_39 = gt3[pp]*gt5[pp];
double DENDRO_40 = -DENDRO_28*DENDRO_8 + DENDRO_34*gt3[pp] + DENDRO_37*gt0[pp] + DENDRO_38*gt5[pp] - DENDRO_39*gt0[pp];
double DENDRO_41 = 1.0/DENDRO_40;
double DENDRO_42 = DENDRO_26*DENDRO_41;
double DENDRO_43 = At1[pp]*DENDRO_29;
double DENDRO_44 = gt1[pp]*gt4[pp];
double DENDRO_45 = gt2[pp]*gt3[pp];
double DENDRO_46 = DENDRO_44 - DENDRO_45;
double DENDRO_47 = At2[pp]*DENDRO_46;
double DENDRO_48 = -DENDRO_47;
double DENDRO_49 = -DENDRO_37 + DENDRO_39;
double DENDRO_50 = At0[pp]*DENDRO_49;
double DENDRO_51 = 2*DENDRO_41;
double DENDRO_52 = At0[pp]*DENDRO_51;
double DENDRO_53 = At1[pp]*DENDRO_31;
double DENDRO_54 = gt0[pp]*gt3[pp];
double DENDRO_55 = -DENDRO_38 + DENDRO_54;
double DENDRO_56 = At2[pp]*DENDRO_55;
double DENDRO_57 = DENDRO_27*DENDRO_41;
double DENDRO_58 = grad2_0_0_alpha[pp];
double DENDRO_59 = grad_1_chi[pp];
double DENDRO_60 = grad_2_chi[pp];
double DENDRO_61 = grad_0_chi[pp];
double DENDRO_62 = DENDRO_29*DENDRO_61 + DENDRO_31*DENDRO_60;
double DENDRO_63 = -DENDRO_35*DENDRO_59 + DENDRO_62;
double DENDRO_64 = 1.0/chi[pp];
double DENDRO_65 = 0.5*DENDRO_64;
double DENDRO_66 = DENDRO_65*gt0[pp];
double DENDRO_67 = grad_0_gt1[pp];
double DENDRO_68 = 1.0*DENDRO_67;
double DENDRO_69 = grad_1_gt0[pp];
double DENDRO_70 = 0.5*DENDRO_69;
double DENDRO_71 = DENDRO_68 - DENDRO_70;
double DENDRO_72 = grad_0_gt2[pp];
double DENDRO_73 = 1.0*DENDRO_72;
double DENDRO_74 = grad_2_gt0[pp];
double DENDRO_75 = 0.5*DENDRO_74;
double DENDRO_76 = DENDRO_73 - DENDRO_75;
double DENDRO_77 = grad_0_gt0[pp];
double DENDRO_78 = 0.5*DENDRO_77;
double DENDRO_79 = DENDRO_29*DENDRO_78 + DENDRO_31*DENDRO_76;
double DENDRO_80 = -DENDRO_35*DENDRO_71 + DENDRO_79;
double DENDRO_81 = DENDRO_63*DENDRO_66 + DENDRO_80;
double DENDRO_82 = grad_1_alpha[pp];
double DENDRO_83 = 12*DENDRO_82;
double DENDRO_84 = DENDRO_41*DENDRO_83;
double DENDRO_85 = DENDRO_31*DENDRO_59;
double DENDRO_86 = DENDRO_46*DENDRO_61;
double DENDRO_87 = DENDRO_55*DENDRO_60;
double DENDRO_88 = DENDRO_85 - DENDRO_86 - DENDRO_87;
double DENDRO_89 = DENDRO_66*DENDRO_88;
double DENDRO_90 = DENDRO_46*DENDRO_78;
double DENDRO_91 = DENDRO_55*DENDRO_76;
double DENDRO_92 = DENDRO_31*DENDRO_71;
double DENDRO_93 = DENDRO_90 + DENDRO_91 - DENDRO_92;
double DENDRO_94 = grad_2_alpha[pp];
double DENDRO_95 = 12*DENDRO_94;
double DENDRO_96 = DENDRO_41*DENDRO_95;
double DENDRO_97 = DENDRO_49*DENDRO_78;
double DENDRO_98 = DENDRO_41*DENDRO_97;
double DENDRO_99 = DENDRO_46*DENDRO_76;
double DENDRO_100 = DENDRO_41*DENDRO_99;
double DENDRO_101 = DENDRO_29*DENDRO_71;
double DENDRO_102 = DENDRO_101*DENDRO_41;
double DENDRO_103 = 1.0*DENDRO_61;
double DENDRO_104 = DENDRO_41*gt0[pp];
double DENDRO_105 = -DENDRO_44 + DENDRO_45;
double DENDRO_106 = DENDRO_29*DENDRO_59;
double DENDRO_107 = DENDRO_37 - DENDRO_39;
double DENDRO_108 = DENDRO_105*DENDRO_60 + DENDRO_106 + DENDRO_107*DENDRO_61;
double DENDRO_109 = 0.5*DENDRO_108;
double DENDRO_110 = DENDRO_104*DENDRO_109;
double DENDRO_111 = DENDRO_64*(-DENDRO_103 + DENDRO_110);
double DENDRO_112 = grad_0_alpha[pp];
double DENDRO_113 = 12*DENDRO_112;
double DENDRO_114 = grad2_0_0_chi[pp];
double DENDRO_115 = -DENDRO_114;
double DENDRO_116 = -DENDRO_90 - DENDRO_91 + DENDRO_92;
double DENDRO_117 = DENDRO_41*DENDRO_60;
double DENDRO_118 = DENDRO_41*DENDRO_59;
double DENDRO_119 = DENDRO_101 - DENDRO_97 - DENDRO_99;
double DENDRO_120 = DENDRO_41*DENDRO_61;
double DENDRO_121 = 2*DENDRO_64;
double DENDRO_122 = grad_0_gt3[pp];
double DENDRO_123 = 0.5*DENDRO_122;
double DENDRO_124 = grad_1_gt1[pp];
double DENDRO_125 = 1.0*DENDRO_124;
double DENDRO_126 = -DENDRO_125;
double DENDRO_127 = DENDRO_123 + DENDRO_126;
double DENDRO_128 = grad_0_gt4[pp];
double DENDRO_129 = grad_2_gt1[pp];
double DENDRO_130 = grad_1_gt2[pp];
double DENDRO_131 = DENDRO_128 + DENDRO_129 - DENDRO_130;
double DENDRO_132 = grad_0_gt5[pp];
double DENDRO_133 = DENDRO_132*DENDRO_31 + DENDRO_29*DENDRO_74;
double DENDRO_134 = -DENDRO_131*DENDRO_35 + DENDRO_133;
double DENDRO_135 = DENDRO_127*DENDRO_134;
double DENDRO_136 = DENDRO_122*DENDRO_35;
double DENDRO_137 = DENDRO_128 - DENDRO_129 + DENDRO_130;
double DENDRO_138 = DENDRO_137*DENDRO_31;
double DENDRO_139 = DENDRO_29*DENDRO_69;
double DENDRO_140 = DENDRO_138 + DENDRO_139;
double DENDRO_141 = -DENDRO_136 + DENDRO_140;
double DENDRO_142 = DENDRO_131*DENDRO_141;
double DENDRO_143 = 0.25*DENDRO_142;
double DENDRO_144 = pow(DENDRO_40, -2);
double DENDRO_145 = DENDRO_144*DENDRO_31;
double DENDRO_146 = 4*DENDRO_145;
double DENDRO_147 = 0.5*DENDRO_132;
double DENDRO_148 = grad_2_gt2[pp];
double DENDRO_149 = 1.0*DENDRO_148;
double DENDRO_150 = -DENDRO_149;
double DENDRO_151 = DENDRO_147 + DENDRO_150;
double DENDRO_152 = DENDRO_122*DENDRO_31;
double DENDRO_153 = DENDRO_46*DENDRO_69;
double DENDRO_154 = DENDRO_137*DENDRO_55;
double DENDRO_155 = DENDRO_152 - DENDRO_153 - DENDRO_154;
double DENDRO_156 = DENDRO_151*DENDRO_155;
double DENDRO_157 = DENDRO_46*DENDRO_74;
double DENDRO_158 = DENDRO_132*DENDRO_55;
double DENDRO_159 = DENDRO_131*DENDRO_31;
double DENDRO_160 = -DENDRO_157 - DENDRO_158 + DENDRO_159;
double DENDRO_161 = DENDRO_137*DENDRO_160;
double DENDRO_162 = 0.25*DENDRO_161;
double DENDRO_163 = -DENDRO_162;
double DENDRO_164 = -DENDRO_152 + DENDRO_153 + DENDRO_154;
double DENDRO_165 = 0.25*DENDRO_132;
double DENDRO_166 = DENDRO_164*DENDRO_165;
double DENDRO_167 = -DENDRO_128 + DENDRO_129 + DENDRO_130;
double DENDRO_168 = DENDRO_157 + DENDRO_158 - DENDRO_159;
double DENDRO_169 = 0.5*DENDRO_168;
double DENDRO_170 = DENDRO_122*DENDRO_29;
double DENDRO_171 = DENDRO_49*DENDRO_69;
double DENDRO_172 = DENDRO_137*DENDRO_46;
double DENDRO_173 = -DENDRO_170 + DENDRO_171 + DENDRO_172;
double DENDRO_174 = 0.25*DENDRO_74;
double DENDRO_175 = DENDRO_173*DENDRO_174;
double DENDRO_176 = DENDRO_132*DENDRO_46;
double DENDRO_177 = DENDRO_49*DENDRO_74;
double DENDRO_178 = DENDRO_131*DENDRO_29;
double DENDRO_179 = DENDRO_176 + DENDRO_177 - DENDRO_178;
double DENDRO_180 = 0.25*DENDRO_69;
double DENDRO_181 = DENDRO_179*DENDRO_180;
double DENDRO_182 = 2*DENDRO_151;
double DENDRO_183 = DENDRO_144*DENDRO_46;
double DENDRO_184 = 4*DENDRO_183;
double DENDRO_185 = DENDRO_173*DENDRO_77;
double DENDRO_186 = 0.25*DENDRO_185;
double DENDRO_187 = -DENDRO_101 + DENDRO_97 + DENDRO_99;
double DENDRO_188 = DENDRO_187*DENDRO_69;
double DENDRO_189 = DENDRO_144*DENDRO_29;
double DENDRO_190 = 4*DENDRO_189;
double DENDRO_191 = 0.5*DENDRO_76;
double DENDRO_192 = DENDRO_179*DENDRO_77;
double DENDRO_193 = 0.25*DENDRO_192;
double DENDRO_194 = DENDRO_187*DENDRO_74;
double DENDRO_195 = DENDRO_168*DENDRO_74;
double DENDRO_196 = DENDRO_132*DENDRO_93;
double DENDRO_197 = 2.0*DENDRO_183;
double DENDRO_198 = DENDRO_164*DENDRO_74;
double DENDRO_199 = DENDRO_137*DENDRO_93;
double DENDRO_200 = 2.0*DENDRO_189;
double DENDRO_201 = DENDRO_134*DENDRO_46;
double DENDRO_202 = grad_2_gt3[pp];
double DENDRO_203 = DENDRO_202*DENDRO_35;
double DENDRO_204 = grad_1_gt5[pp];
double DENDRO_205 = DENDRO_204*DENDRO_31;
double DENDRO_206 = DENDRO_167*DENDRO_29;
double DENDRO_207 = DENDRO_205 + DENDRO_206;
double DENDRO_208 = -DENDRO_203 + DENDRO_207;
double DENDRO_209 = DENDRO_208*DENDRO_31;
double DENDRO_210 = DENDRO_141*DENDRO_29;
double DENDRO_211 = grad_2_gt5[pp];
double DENDRO_212 = 0.5*DENDRO_211;
double DENDRO_213 = DENDRO_212*DENDRO_31;
double DENDRO_214 = 0.5*DENDRO_204;
double DENDRO_215 = grad_2_gt4[pp];
double DENDRO_216 = 1.0*DENDRO_215;
double DENDRO_217 = -DENDRO_216;
double DENDRO_218 = DENDRO_214 + DENDRO_217;
double DENDRO_219 = -DENDRO_151*DENDRO_29 + DENDRO_213 + DENDRO_218*DENDRO_35;
double DENDRO_220 = DENDRO_219*DENDRO_55;
double DENDRO_221 = grad_1_gt3[pp];
double DENDRO_222 = 0.5*DENDRO_221;
double DENDRO_223 = DENDRO_222*DENDRO_35;
double DENDRO_224 = grad_1_gt4[pp];
double DENDRO_225 = 1.0*DENDRO_224;
double DENDRO_226 = 0.5*DENDRO_202;
double DENDRO_227 = DENDRO_225 - DENDRO_226;
double DENDRO_228 = DENDRO_227*DENDRO_31;
double DENDRO_229 = DENDRO_127*DENDRO_29;
double DENDRO_230 = -DENDRO_223 + DENDRO_228 - DENDRO_229;
double DENDRO_231 = DENDRO_230*DENDRO_35;
double DENDRO_232 = DENDRO_49*DENDRO_80;
double DENDRO_233 = DENDRO_144*(DENDRO_201 - 1.0*DENDRO_209 - 1.0*DENDRO_210 + DENDRO_220 + DENDRO_231 + DENDRO_232);
double DENDRO_234 = 2.0*DENDRO_69;
double DENDRO_235 = -DENDRO_176 - DENDRO_177 + DENDRO_178;
double DENDRO_236 = DENDRO_235*DENDRO_46;
double DENDRO_237 = DENDRO_202*DENDRO_29;
double DENDRO_238 = DENDRO_204*DENDRO_46;
double DENDRO_239 = DENDRO_167*DENDRO_49;
double DENDRO_240 = DENDRO_237 - DENDRO_238 - DENDRO_239;
double DENDRO_241 = DENDRO_240*DENDRO_31;
double DENDRO_242 = DENDRO_170 - DENDRO_171 - DENDRO_172;
double DENDRO_243 = DENDRO_242*DENDRO_29;
double DENDRO_244 = DENDRO_212*DENDRO_46;
double DENDRO_245 = DENDRO_218*DENDRO_29;
double DENDRO_246 = DENDRO_151*DENDRO_49;
double DENDRO_247 = -DENDRO_244 - DENDRO_245 + DENDRO_246;
double DENDRO_248 = DENDRO_247*DENDRO_55;
double DENDRO_249 = DENDRO_222*DENDRO_29;
double DENDRO_250 = DENDRO_127*DENDRO_49 - DENDRO_227*DENDRO_46 + DENDRO_249;
double DENDRO_251 = DENDRO_250*DENDRO_35;
double DENDRO_252 = DENDRO_119*DENDRO_49;
double DENDRO_253 = DENDRO_144*(DENDRO_236 - 1.0*DENDRO_241 - 1.0*DENDRO_243 + DENDRO_248 + DENDRO_251 + DENDRO_252);
double DENDRO_254 = 2.0*DENDRO_77;
double DENDRO_255 = DENDRO_160*DENDRO_46;
double DENDRO_256 = DENDRO_202*DENDRO_31;
double DENDRO_257 = DENDRO_204*DENDRO_55;
double DENDRO_258 = DENDRO_167*DENDRO_46;
double DENDRO_259 = DENDRO_256 - DENDRO_257 - DENDRO_258;
double DENDRO_260 = DENDRO_259*DENDRO_31;
double DENDRO_261 = DENDRO_155*DENDRO_29;
double DENDRO_262 = DENDRO_212*DENDRO_55;
double DENDRO_263 = DENDRO_151*DENDRO_46;
double DENDRO_264 = DENDRO_218*DENDRO_31;
double DENDRO_265 = -DENDRO_262 + DENDRO_263 - DENDRO_264;
double DENDRO_266 = DENDRO_265*DENDRO_55;
double DENDRO_267 = DENDRO_222*DENDRO_31;
double DENDRO_268 = DENDRO_127*DENDRO_46 - DENDRO_227*DENDRO_55 + DENDRO_267;
double DENDRO_269 = DENDRO_268*DENDRO_35;
double DENDRO_270 = DENDRO_116*DENDRO_49;
double DENDRO_271 = DENDRO_144*(DENDRO_255 - 1.0*DENDRO_260 - 1.0*DENDRO_261 + DENDRO_266 + DENDRO_269 + DENDRO_270);
double DENDRO_272 = 2.0*DENDRO_74;
double DENDRO_273 = grad2_2_2_chi[pp];
double DENDRO_274 = pow(DENDRO_60, 2);
double DENDRO_275 = 3*DENDRO_64;
double DENDRO_276 = DENDRO_55*(2*DENDRO_273 - DENDRO_274*DENDRO_275);
double DENDRO_277 = grad2_1_1_chi[pp];
double DENDRO_278 = pow(DENDRO_59, 2);
double DENDRO_279 = DENDRO_35*(-DENDRO_275*DENDRO_278 + 2*DENDRO_277);
double DENDRO_280 = pow(DENDRO_61, 2);
double DENDRO_281 = DENDRO_49*(2*DENDRO_114 - DENDRO_275*DENDRO_280);
double DENDRO_282 = grad2_1_2_chi[pp];
double DENDRO_283 = DENDRO_275*DENDRO_60;
double DENDRO_284 = 2*DENDRO_31;
double DENDRO_285 = DENDRO_284*(2*DENDRO_282 - DENDRO_283*DENDRO_59);
double DENDRO_286 = grad2_0_2_chi[pp];
double DENDRO_287 = 2*DENDRO_46;
double DENDRO_288 = DENDRO_287*(-DENDRO_283*DENDRO_61 + 2*DENDRO_286);
double DENDRO_289 = grad2_0_1_chi[pp];
double DENDRO_290 = DENDRO_59*DENDRO_61;
double DENDRO_291 = 2*DENDRO_29;
double DENDRO_292 = DENDRO_291*(-DENDRO_275*DENDRO_290 + 2*DENDRO_289);
double DENDRO_293 = -1.0*DENDRO_201 + DENDRO_209 + DENDRO_210 - DENDRO_220 - DENDRO_231 - DENDRO_232;
double DENDRO_294 = 2*DENDRO_118*DENDRO_293;
double DENDRO_295 = -1.0*DENDRO_255 + DENDRO_260 + DENDRO_261 - DENDRO_266 - DENDRO_269 - DENDRO_270;
double DENDRO_296 = 2*DENDRO_117*DENDRO_295;
double DENDRO_297 = -1.0*DENDRO_236 + DENDRO_241 + DENDRO_243 - DENDRO_248 - DENDRO_251 - DENDRO_252;
double DENDRO_298 = 2*DENDRO_120*DENDRO_297;
double DENDRO_299 = DENDRO_276 + DENDRO_279 + DENDRO_281 - DENDRO_285 + DENDRO_288 - DENDRO_292 + DENDRO_294 + DENDRO_296 + DENDRO_298;
double DENDRO_300 = DENDRO_104*DENDRO_64;
double DENDRO_301 = DENDRO_144*DENDRO_55;
double DENDRO_302 = 4*DENDRO_301;
double DENDRO_303 = DENDRO_160*DENDRO_302;
double DENDRO_304 = 0.25*DENDRO_122;
double DENDRO_305 = DENDRO_144*DENDRO_35;
double DENDRO_306 = 4*DENDRO_305;
double DENDRO_307 = DENDRO_141*DENDRO_306;
double DENDRO_308 = 0.25*DENDRO_128;
double DENDRO_309 = -DENDRO_308;
double DENDRO_310 = 0.75*DENDRO_130;
double DENDRO_311 = 0.25*DENDRO_129;
double DENDRO_312 = DENDRO_306*(DENDRO_309 + DENDRO_310 + DENDRO_311);
double DENDRO_313 = DENDRO_73 + DENDRO_75;
double DENDRO_314 = DENDRO_144*DENDRO_49;
double DENDRO_315 = 4*DENDRO_314;
double DENDRO_316 = DENDRO_179*DENDRO_74;
double DENDRO_317 = 3.0*DENDRO_301;
double DENDRO_318 = DENDRO_173*DENDRO_69;
double DENDRO_319 = 3.0*DENDRO_305;
double DENDRO_320 = 6.0*DENDRO_77;
double DENDRO_321 = pow(chi[pp], -2);
double DENDRO_322 = grad_0_Gt0[pp];
double DENDRO_323 = grad_0_Gt1[pp];
double DENDRO_324 = 4*gt1[pp];
double DENDRO_325 = grad_0_Gt2[pp];
double DENDRO_326 = 4*gt2[pp];
double DENDRO_327 = 0.5*DENDRO_71;
double DENDRO_328 = DENDRO_134*DENDRO_327;
double DENDRO_329 = DENDRO_41*DENDRO_46;
double DENDRO_330 = 4*DENDRO_329;
double DENDRO_331 = DENDRO_134*DENDRO_304;
double DENDRO_332 = 0.5*DENDRO_167;
double DENDRO_333 = DENDRO_141*DENDRO_327;
double DENDRO_334 = DENDRO_41*DENDRO_55;
double DENDRO_335 = 2.0*DENDRO_334;
double DENDRO_336 = DENDRO_35*DENDRO_41;
double DENDRO_337 = 2.0*DENDRO_336;
double DENDRO_338 = DENDRO_41*DENDRO_49;
double DENDRO_339 = 2.0*DENDRO_338;
double DENDRO_340 = DENDRO_141*DENDRO_69;
double DENDRO_341 = DENDRO_122*DENDRO_80;
double DENDRO_342 = DENDRO_134*DENDRO_69;
double DENDRO_343 = DENDRO_131*DENDRO_80;
double DENDRO_344 = 4.0*DENDRO_41;
double DENDRO_345 = DENDRO_31*DENDRO_344;
double DENDRO_346 = DENDRO_29*DENDRO_344;
double DENDRO_347 = 0.25*DENDRO_130;
double DENDRO_348 = 0.75*DENDRO_129;
double DENDRO_349 = 4*DENDRO_144;
double DENDRO_350 = -DENDRO_134*DENDRO_302*(DENDRO_309 + DENDRO_347 + DENDRO_348) + DENDRO_146*(DENDRO_141*DENDRO_332 + DENDRO_331) - DENDRO_184*(DENDRO_167*DENDRO_80 + DENDRO_328) + DENDRO_190*(-2*DENDRO_127*DENDRO_80 + DENDRO_333) - DENDRO_197*(DENDRO_342 + DENDRO_343) + DENDRO_200*(DENDRO_340 + DENDRO_341) - DENDRO_232*DENDRO_349*(DENDRO_68 + DENDRO_70) - DENDRO_280*DENDRO_321 + 4*DENDRO_322*gt0[pp] + DENDRO_323*DENDRO_324 + DENDRO_325*DENDRO_326 + DENDRO_330*grad2_0_2_gt0[pp] + DENDRO_335*grad2_2_2_gt0[pp] + DENDRO_337*grad2_1_1_gt0[pp] + DENDRO_339*grad2_0_0_gt0[pp] - DENDRO_345*grad2_1_2_gt0[pp] - DENDRO_346*grad2_0_1_gt0[pp];
double DENDRO_351 = 3*alpha[pp];
double DENDRO_352 = grad2_2_2_alpha[pp];
double DENDRO_353 = DENDRO_65*gt5[pp];
double DENDRO_354 = DENDRO_219 + DENDRO_353*DENDRO_63;
double DENDRO_355 = 4*DENDRO_82;
double DENDRO_356 = DENDRO_355*DENDRO_41;
double DENDRO_357 = DENDRO_46*DENDRO_60;
double DENDRO_358 = DENDRO_49*DENDRO_61;
double DENDRO_359 = DENDRO_106 - DENDRO_357 - DENDRO_358;
double DENDRO_360 = DENDRO_353*DENDRO_359;
double DENDRO_361 = 4*DENDRO_112;
double DENDRO_362 = DENDRO_361*DENDRO_41;
double DENDRO_363 = DENDRO_262*DENDRO_41;
double DENDRO_364 = DENDRO_263*DENDRO_41;
double DENDRO_365 = DENDRO_264*DENDRO_41;
double DENDRO_366 = 1.0*DENDRO_60;
double DENDRO_367 = DENDRO_41*gt5[pp];
double DENDRO_368 = DENDRO_38 - DENDRO_54;
double DENDRO_369 = DENDRO_105*DENDRO_61 + DENDRO_368*DENDRO_60 + DENDRO_85;
double DENDRO_370 = 0.5*DENDRO_369;
double DENDRO_371 = DENDRO_367*DENDRO_370;
double DENDRO_372 = DENDRO_64*(-DENDRO_366 + DENDRO_371);
double DENDRO_373 = 4*DENDRO_94;
double DENDRO_374 = -DENDRO_273;
double DENDRO_375 = DENDRO_212*DENDRO_368;
double DENDRO_376 = -DENDRO_214;
double DENDRO_377 = DENDRO_216 + DENDRO_376;
double DENDRO_378 = DENDRO_31*DENDRO_377;
double DENDRO_379 = -DENDRO_147;
double DENDRO_380 = DENDRO_149 + DENDRO_379;
double DENDRO_381 = DENDRO_105*DENDRO_380;
double DENDRO_382 = -DENDRO_33 + DENDRO_34;
double DENDRO_383 = DENDRO_213 + DENDRO_29*DENDRO_380 + DENDRO_377*DENDRO_382;
double DENDRO_384 = DENDRO_105*DENDRO_212 + DENDRO_107*DENDRO_380 + DENDRO_29*DENDRO_377;
double DENDRO_385 = 0.5*DENDRO_218;
double DENDRO_386 = DENDRO_134*DENDRO_385;
double DENDRO_387 = -DENDRO_386;
double DENDRO_388 = DENDRO_137*DENDRO_219;
double DENDRO_389 = 0.5*DENDRO_151;
double DENDRO_390 = -DENDRO_235*DENDRO_389;
double DENDRO_391 = 2*DENDRO_76;
double DENDRO_392 = DENDRO_160*DENDRO_211;
double DENDRO_393 = 0.25*DENDRO_392;
double DENDRO_394 = DENDRO_132*DENDRO_265;
double DENDRO_395 = DENDRO_240*DENDRO_389;
double DENDRO_396 = -DENDRO_395;
double DENDRO_397 = DENDRO_137*DENDRO_247;
double DENDRO_398 = DENDRO_211*DENDRO_259;
double DENDRO_399 = 0.25*DENDRO_398;
double DENDRO_400 = DENDRO_204*DENDRO_265;
double DENDRO_401 = DENDRO_240*DENDRO_76;
double DENDRO_402 = DENDRO_167*DENDRO_235;
double DENDRO_403 = 0.25*DENDRO_402;
double DENDRO_404 = 0.25*DENDRO_204;
double DENDRO_405 = DENDRO_160*DENDRO_404;
double DENDRO_406 = DENDRO_165*DENDRO_259;
double DENDRO_407 = DENDRO_174*DENDRO_240;
double DENDRO_408 = 0.5*DENDRO_137;
double DENDRO_409 = DENDRO_167*DENDRO_247;
double DENDRO_410 = 2.0*DENDRO_145;
double DENDRO_411 = DENDRO_144*DENDRO_293;
double DENDRO_412 = 2.0*DENDRO_411;
double DENDRO_413 = DENDRO_144*DENDRO_295;
double DENDRO_414 = 2.0*DENDRO_413;
double DENDRO_415 = DENDRO_144*DENDRO_297;
double DENDRO_416 = 2.0*DENDRO_415;
double DENDRO_417 = DENDRO_247*DENDRO_74;
double DENDRO_418 = -DENDRO_276 - DENDRO_279 - DENDRO_281 + DENDRO_285 - DENDRO_288 + DENDRO_292 - DENDRO_294 - DENDRO_296 - DENDRO_298;
double DENDRO_419 = DENDRO_418*DENDRO_64;
double DENDRO_420 = DENDRO_248*DENDRO_349;
double DENDRO_421 = DENDRO_220*DENDRO_349;
double DENDRO_422 = -DENDRO_311;
double DENDRO_423 = DENDRO_306*(DENDRO_308 + DENDRO_310 + DENDRO_422);
double DENDRO_424 = DENDRO_315*(-DENDRO_174 + DENDRO_73);
double DENDRO_425 = DENDRO_204*DENDRO_259;
double DENDRO_426 = DENDRO_132*DENDRO_160;
double DENDRO_427 = 3.0*DENDRO_314;
double DENDRO_428 = 6.0*DENDRO_144;
double DENDRO_429 = grad_2_Gt0[pp];
double DENDRO_430 = grad_2_Gt1[pp];
double DENDRO_431 = 4*gt4[pp];
double DENDRO_432 = grad_2_Gt2[pp];
double DENDRO_433 = -DENDRO_208*DENDRO_385;
double DENDRO_434 = DENDRO_134*DENDRO_227;
double DENDRO_435 = DENDRO_131*DENDRO_208;
double DENDRO_436 = 0.25*DENDRO_435;
double DENDRO_437 = 0.25*DENDRO_202;
double DENDRO_438 = DENDRO_134*DENDRO_437;
double DENDRO_439 = DENDRO_204*DENDRO_208;
double DENDRO_440 = DENDRO_202*DENDRO_219;
double DENDRO_441 = DENDRO_134*DENDRO_204;
double DENDRO_442 = DENDRO_131*DENDRO_219;
double DENDRO_443 = 0.75*DENDRO_128;
double DENDRO_444 = -DENDRO_134*DENDRO_315*(DENDRO_347 + DENDRO_422 + DENDRO_443) + DENDRO_146*(2*DENDRO_219*DENDRO_227 + DENDRO_433) + DENDRO_190*(DENDRO_434 + DENDRO_436) + DENDRO_190*(DENDRO_208*DENDRO_408 + DENDRO_438) - DENDRO_197*(DENDRO_441 + DENDRO_442) - DENDRO_208*DENDRO_306*(DENDRO_225 - DENDRO_437) - DENDRO_274*DENDRO_321 + DENDRO_326*DENDRO_429 + DENDRO_330*grad2_0_2_gt5[pp] + DENDRO_335*grad2_2_2_gt5[pp] + DENDRO_337*grad2_1_1_gt5[pp] + DENDRO_339*grad2_0_0_gt5[pp] - DENDRO_345*grad2_1_2_gt5[pp] - DENDRO_346*grad2_0_1_gt5[pp] + DENDRO_410*(DENDRO_439 + DENDRO_440) + DENDRO_430*DENDRO_431 + 4*DENDRO_432*gt5[pp];
double DENDRO_445 = grad2_1_1_alpha[pp];
double DENDRO_446 = DENDRO_65*gt3[pp];
double DENDRO_447 = DENDRO_373*DENDRO_41;
double DENDRO_448 = 1.0*DENDRO_59;
double DENDRO_449 = -DENDRO_448;
double DENDRO_450 = DENDRO_382*DENDRO_59 + DENDRO_62;
double DENDRO_451 = DENDRO_41*gt3[pp];
double DENDRO_452 = 0.5*DENDRO_451;
double DENDRO_453 = DENDRO_450*DENDRO_452;
double DENDRO_454 = DENDRO_228*DENDRO_41;
double DENDRO_455 = -DENDRO_223*DENDRO_41 - DENDRO_229*DENDRO_41 + DENDRO_454;
double DENDRO_456 = -DENDRO_277;
double DENDRO_457 = -DENDRO_123;
double DENDRO_458 = DENDRO_125 + DENDRO_457;
double DENDRO_459 = DENDRO_105*DENDRO_458 + DENDRO_227*DENDRO_368 + DENDRO_267;
double DENDRO_460 = DENDRO_222*DENDRO_382;
double DENDRO_461 = DENDRO_29*DENDRO_458;
double DENDRO_462 = DENDRO_105*DENDRO_227 + DENDRO_107*DENDRO_458 + DENDRO_249;
double DENDRO_463 = DENDRO_240*DENDRO_71;
double DENDRO_464 = DENDRO_167*DENDRO_242;
double DENDRO_465 = 0.25*DENDRO_464;
double DENDRO_466 = DENDRO_155*DENDRO_218;
double DENDRO_467 = DENDRO_137*DENDRO_259;
double DENDRO_468 = 0.25*DENDRO_467;
double DENDRO_469 = DENDRO_155*DENDRO_404;
double DENDRO_470 = 0.5*DENDRO_131;
double DENDRO_471 = DENDRO_180*DENDRO_240;
double DENDRO_472 = 0.5*DENDRO_127;
double DENDRO_473 = DENDRO_240*DENDRO_472;
double DENDRO_474 = -DENDRO_473;
double DENDRO_475 = DENDRO_131*DENDRO_250;
double DENDRO_476 = 0.5*DENDRO_227;
double DENDRO_477 = DENDRO_259*DENDRO_476;
double DENDRO_478 = 2*DENDRO_218*DENDRO_268;
double DENDRO_479 = DENDRO_208*DENDRO_221;
double DENDRO_480 = 0.25*DENDRO_479;
double DENDRO_481 = DENDRO_202*DENDRO_230;
double DENDRO_482 = DENDRO_155*DENDRO_476;
double DENDRO_483 = DENDRO_131*DENDRO_268;
double DENDRO_484 = -DENDRO_242*DENDRO_472;
double DENDRO_485 = 2*DENDRO_250*DENDRO_71;
double DENDRO_486 = DENDRO_141*DENDRO_221;
double DENDRO_487 = 0.25*DENDRO_486;
double DENDRO_488 = DENDRO_122*DENDRO_230;
double DENDRO_489 = DENDRO_167*DENDRO_250;
double DENDRO_490 = DENDRO_204*DENDRO_268;
double DENDRO_491 = DENDRO_137*DENDRO_268;
double DENDRO_492 = DENDRO_250*DENDRO_69;
double DENDRO_493 = DENDRO_259*DENDRO_302;
double DENDRO_494 = -DENDRO_347;
double DENDRO_495 = DENDRO_302*(DENDRO_308 + DENDRO_348 + DENDRO_494);
double DENDRO_496 = DENDRO_251*DENDRO_349;
double DENDRO_497 = DENDRO_315*(-DENDRO_180 + DENDRO_68);
double DENDRO_498 = DENDRO_315*(DENDRO_311 + DENDRO_443 + DENDRO_494);
double DENDRO_499 = grad_1_Gt0[pp];
double DENDRO_500 = grad_1_Gt1[pp];
double DENDRO_501 = grad_1_Gt2[pp];
double DENDRO_502 = DENDRO_141*DENDRO_437;
double DENDRO_503 = DENDRO_208*DENDRO_304;
double DENDRO_504 = DENDRO_202*DENDRO_208;
double DENDRO_505 = DENDRO_122*DENDRO_141;
double DENDRO_506 = -DENDRO_184*(DENDRO_123*DENDRO_208 + DENDRO_502) - DENDRO_184*(DENDRO_141*DENDRO_226 + DENDRO_503) - DENDRO_269*DENDRO_349*(DENDRO_225 + DENDRO_226) - DENDRO_278*DENDRO_321 - DENDRO_317*DENDRO_504 + DENDRO_324*DENDRO_499 + DENDRO_330*grad2_0_2_gt3[pp] + DENDRO_335*grad2_2_2_gt3[pp] + DENDRO_337*grad2_1_1_gt3[pp] + DENDRO_339*grad2_0_0_gt3[pp] - DENDRO_345*grad2_1_2_gt3[pp] - DENDRO_346*grad2_0_1_gt3[pp] - DENDRO_427*DENDRO_505 + DENDRO_431*DENDRO_501 + 4*DENDRO_500*gt3[pp];
double DENDRO_507 = DENDRO_105*DENDRO_78 + DENDRO_368*DENDRO_76 + DENDRO_92;
double DENDRO_508 = DENDRO_382*DENDRO_71 + DENDRO_79;
double DENDRO_509 = DENDRO_107*DENDRO_78;
double DENDRO_510 = DENDRO_105*DENDRO_76;
double DENDRO_511 = DENDRO_160*DENDRO_191;
double DENDRO_512 = DENDRO_235*DENDRO_77;
double DENDRO_513 = 0.25*DENDRO_512;
double DENDRO_514 = DENDRO_119*DENDRO_74;
double DENDRO_515 = DENDRO_155*DENDRO_165;
double DENDRO_516 = DENDRO_180*DENDRO_235;
double DENDRO_517 = DENDRO_174*DENDRO_242;
double DENDRO_518 = DENDRO_155*DENDRO_191;
double DENDRO_519 = DENDRO_242*DENDRO_77;
double DENDRO_520 = 0.25*DENDRO_519;
double DENDRO_521 = DENDRO_119*DENDRO_69;
double DENDRO_522 = DENDRO_116*DENDRO_137;
double DENDRO_523 = DENDRO_116*DENDRO_132;
double DENDRO_524 = DENDRO_235*DENDRO_74;
double DENDRO_525 = DENDRO_242*DENDRO_69;
double DENDRO_526 = DENDRO_259*DENDRO_74;
double DENDRO_527 = DENDRO_132*DENDRO_155;
double DENDRO_528 = DENDRO_161 + DENDRO_527;
double DENDRO_529 = DENDRO_208*DENDRO_69;
double DENDRO_530 = DENDRO_137*DENDRO_141;
double DENDRO_531 = DENDRO_122*DENDRO_134;
double DENDRO_532 = DENDRO_530 + DENDRO_531;
double DENDRO_533 = DENDRO_165*DENDRO_235;
double DENDRO_534 = DENDRO_151*DENDRO_265;
double DENDRO_535 = -DENDRO_534;
double DENDRO_536 = DENDRO_219*DENDRO_332 + 0.25*DENDRO_441;
double DENDRO_537 = DENDRO_137*DENDRO_242;
double DENDRO_538 = 0.25*DENDRO_537;
double DENDRO_539 = DENDRO_141*DENDRO_476;
double DENDRO_540 = -DENDRO_208*DENDRO_472;
double DENDRO_541 = DENDRO_539 + DENDRO_540;
double DENDRO_542 = DENDRO_119*DENDRO_76;
double DENDRO_543 = DENDRO_160*DENDRO_174;
double DENDRO_544 = 0.25*DENDRO_342 + DENDRO_408*DENDRO_80;
double DENDRO_545 = 0.25*DENDRO_134;
double DENDRO_546 = DENDRO_167*DENDRO_545 + DENDRO_214*DENDRO_80;
double DENDRO_547 = DENDRO_116*DENDRO_212;
double DENDRO_548 = -DENDRO_160*DENDRO_389;
double DENDRO_549 = DENDRO_131*DENDRO_545;
double DENDRO_550 = DENDRO_137*DENDRO_545 + DENDRO_219*DENDRO_70;
double DENDRO_551 = DENDRO_174*DENDRO_235;
double DENDRO_552 = DENDRO_247*DENDRO_78;
double DENDRO_553 = DENDRO_191*DENDRO_235 + DENDRO_552;
double DENDRO_554 = DENDRO_127*DENDRO_219;
double DENDRO_555 = 0.5*DENDRO_434;
double DENDRO_556 = -DENDRO_554 + DENDRO_555;
double DENDRO_557 = DENDRO_167*DENDRO_208;
double DENDRO_558 = 0.25*DENDRO_557;
double DENDRO_559 = DENDRO_123*DENDRO_219;
double DENDRO_560 = DENDRO_141*DENDRO_404 + DENDRO_559;
double DENDRO_561 = 0.25*DENDRO_137;
double DENDRO_562 = DENDRO_247*DENDRO_70;
double DENDRO_563 = DENDRO_235*DENDRO_561 + DENDRO_562;
double DENDRO_564 = DENDRO_265*DENDRO_332;
double DENDRO_565 = DENDRO_405 + DENDRO_406;
double DENDRO_566 = DENDRO_259*DENDRO_389;
double DENDRO_567 = -DENDRO_566;
double DENDRO_568 = 0.25*DENDRO_155;
double DENDRO_569 = DENDRO_211*DENDRO_568;
double DENDRO_570 = DENDRO_265*DENDRO_408 + DENDRO_569;
double DENDRO_571 = DENDRO_165*DENDRO_242 + DENDRO_562;
double DENDRO_572 = DENDRO_208*DENDRO_327;
double DENDRO_573 = -0.5*DENDRO_135 + DENDRO_227*DENDRO_80;
double DENDRO_574 = DENDRO_116*DENDRO_214;
double DENDRO_575 = DENDRO_191*DENDRO_259;
double DENDRO_576 = 0.25*DENDRO_167;
double DENDRO_577 = DENDRO_160*DENDRO_576;
double DENDRO_578 = 0.25*DENDRO_240;
double DENDRO_579 = DENDRO_578*DENDRO_77;
double DENDRO_580 = DENDRO_191*DENDRO_242;
double DENDRO_581 = DENDRO_579 + DENDRO_580;
double DENDRO_582 = DENDRO_119*DENDRO_408 + DENDRO_516;
double DENDRO_583 = DENDRO_167*DENDRO_259;
double DENDRO_584 = DENDRO_155*DENDRO_204;
double DENDRO_585 = DENDRO_467 + DENDRO_584;
double DENDRO_586 = 1.0*DENDRO_305;
double DENDRO_587 = -DENDRO_286;
double DENDRO_588 = DENDRO_132*DENDRO_368;
double DENDRO_589 = DENDRO_105*DENDRO_74;
double DENDRO_590 = 0.5*DENDRO_117;
double DENDRO_591 = DENDRO_131*DENDRO_382 + DENDRO_133;
double DENDRO_592 = 0.5*DENDRO_118;
double DENDRO_593 = DENDRO_105*DENDRO_132;
double DENDRO_594 = DENDRO_107*DENDRO_74;
double DENDRO_595 = 0.5*DENDRO_120;
double DENDRO_596 = 2.0*DENDRO_429;
double DENDRO_597 = DENDRO_596*gt0[pp];
double DENDRO_598 = 2.0*DENDRO_430;
double DENDRO_599 = DENDRO_598*gt1[pp];
double DENDRO_600 = 2.0*gt2[pp];
double DENDRO_601 = DENDRO_322*DENDRO_600;
double DENDRO_602 = DENDRO_432*DENDRO_600;
double DENDRO_603 = 2.0*gt4[pp];
double DENDRO_604 = DENDRO_323*DENDRO_603;
double DENDRO_605 = 2.0*gt5[pp];
double DENDRO_606 = DENDRO_325*DENDRO_605;
double DENDRO_607 = DENDRO_321*DENDRO_60;
double DENDRO_608 = -DENDRO_607*DENDRO_61;
double DENDRO_609 = DENDRO_330*grad2_0_2_gt2[pp];
double DENDRO_610 = DENDRO_335*grad2_2_2_gt2[pp];
double DENDRO_611 = DENDRO_337*grad2_1_1_gt2[pp];
double DENDRO_612 = DENDRO_339*grad2_0_0_gt2[pp];
double DENDRO_613 = -DENDRO_345*grad2_1_2_gt2[pp];
double DENDRO_614 = -DENDRO_346*grad2_0_1_gt2[pp];
double DENDRO_615 = DENDRO_41*gt2[pp];
double DENDRO_616 = -DENDRO_121*(DENDRO_587 + DENDRO_590*(DENDRO_159 + DENDRO_588 + DENDRO_589) + DENDRO_591*DENDRO_592 + DENDRO_595*(DENDRO_178 + DENDRO_593 + DENDRO_594)) + DENDRO_130*DENDRO_412 + DENDRO_148*DENDRO_414 + DENDRO_416*DENDRO_72 + DENDRO_419*DENDRO_615 + DENDRO_597 + DENDRO_599 + DENDRO_601 + DENDRO_602 + DENDRO_604 + DENDRO_606 + DENDRO_608 + DENDRO_609 + DENDRO_610 + DENDRO_611 + DENDRO_612 + DENDRO_613 + DENDRO_614;
double DENDRO_617 = grad2_0_2_alpha[pp];
double DENDRO_618 = DENDRO_157*DENDRO_41;
double DENDRO_619 = DENDRO_158*DENDRO_41;
double DENDRO_620 = DENDRO_159*DENDRO_41;
double DENDRO_621 = -DENDRO_61;
double DENDRO_622 = DENDRO_369*DENDRO_615;
double DENDRO_623 = DENDRO_64*(DENDRO_621 + DENDRO_622);
double DENDRO_624 = 2.0*DENDRO_94;
double DENDRO_625 = DENDRO_176*DENDRO_41;
double DENDRO_626 = DENDRO_177*DENDRO_41;
double DENDRO_627 = DENDRO_178*DENDRO_41;
double DENDRO_628 = -DENDRO_60;
double DENDRO_629 = DENDRO_108*DENDRO_615;
double DENDRO_630 = DENDRO_64*(DENDRO_628 + DENDRO_629);
double DENDRO_631 = 2.0*DENDRO_112;
double DENDRO_632 = 2.0*DENDRO_82;
double DENDRO_633 = DENDRO_64*gt2[pp];
double DENDRO_634 = DENDRO_41*(DENDRO_134 + DENDRO_63*DENDRO_633);
double DENDRO_635 = -4*DENDRO_617 + DENDRO_624*(-DENDRO_618 - DENDRO_619 + DENDRO_620 + DENDRO_623) + DENDRO_631*(-DENDRO_625 - DENDRO_626 + DENDRO_627 + DENDRO_630) + DENDRO_632*DENDRO_634;
double DENDRO_636 = DENDRO_240*DENDRO_74;
double DENDRO_637 = DENDRO_132*DENDRO_242 + DENDRO_636;
double DENDRO_638 = 0.5*DENDRO_392;
double DENDRO_639 = 0.25*DENDRO_583;
double DENDRO_640 = DENDRO_160*DENDRO_165 + DENDRO_547;
double DENDRO_641 = DENDRO_141*DENDRO_385;
double DENDRO_642 = -DENDRO_641;
double DENDRO_643 = -DENDRO_242*DENDRO_389;
double DENDRO_644 = DENDRO_405 + DENDRO_569;
double DENDRO_645 = DENDRO_515 + DENDRO_574;
double DENDRO_646 = DENDRO_174*DENDRO_259;
double DENDRO_647 = DENDRO_119*DENDRO_332;
double DENDRO_648 = DENDRO_240*DENDRO_69;
double DENDRO_649 = DENDRO_537 + DENDRO_648;
double DENDRO_650 = DENDRO_134*DENDRO_202;
double DENDRO_651 = DENDRO_141*DENDRO_204 + DENDRO_650;
double DENDRO_652 = 0.25*DENDRO_530;
double DENDRO_653 = DENDRO_226*DENDRO_80;
double DENDRO_654 = DENDRO_180*DENDRO_208 + DENDRO_653;
double DENDRO_655 = DENDRO_145*(DENDRO_557 + DENDRO_651) - DENDRO_184*(DENDRO_546 + DENDRO_549) - DENDRO_184*(-DENDRO_218*DENDRO_80 + DENDRO_550) + DENDRO_190*(DENDRO_143 + DENDRO_573) + DENDRO_190*(DENDRO_652 + DENDRO_654) - DENDRO_302*(DENDRO_387 + DENDRO_536) - DENDRO_306*(DENDRO_502 + DENDRO_541) - DENDRO_315*(0.5*DENDRO_343 + DENDRO_544);
double DENDRO_656 = DENDRO_160*DENDRO_202;
double DENDRO_657 = DENDRO_122*DENDRO_235;
double DENDRO_658 = 0.25*DENDRO_439;
double DENDRO_659 = DENDRO_218*DENDRO_265;
double DENDRO_660 = -DENDRO_659;
double DENDRO_661 = DENDRO_165*DENDRO_240 + DENDRO_247*DENDRO_470;
double DENDRO_662 = DENDRO_227*DENDRO_230;
double DENDRO_663 = DENDRO_259*DENDRO_437;
double DENDRO_664 = DENDRO_250*DENDRO_408;
double DENDRO_665 = DENDRO_240*DENDRO_304 + DENDRO_664;
double DENDRO_666 = DENDRO_235*DENDRO_327;
double DENDRO_667 = DENDRO_580 + DENDRO_666;
double DENDRO_668 = DENDRO_247*DENDRO_71 + 0.5*DENDRO_401;
double DENDRO_669 = DENDRO_131*DENDRO_235;
double DENDRO_670 = 0.25*DENDRO_669;
double DENDRO_671 = DENDRO_208*DENDRO_561 + DENDRO_559;
double DENDRO_672 = DENDRO_265*DENDRO_470;
double DENDRO_673 = DENDRO_160*DENDRO_385;
double DENDRO_674 = -DENDRO_673;
double DENDRO_675 = DENDRO_212*DENDRO_268;
double DENDRO_676 = -DENDRO_259*DENDRO_385 + DENDRO_675;
double DENDRO_677 = DENDRO_147*DENDRO_250;
double DENDRO_678 = DENDRO_131*DENDRO_578 + DENDRO_677;
double DENDRO_679 = DENDRO_240*DENDRO_576;
double DENDRO_680 = DENDRO_123*DENDRO_247 + DENDRO_240*DENDRO_561;
double DENDRO_681 = DENDRO_208*DENDRO_476;
double DENDRO_682 = DENDRO_219*DENDRO_222;
double DENDRO_683 = DENDRO_208*DENDRO_437 + DENDRO_682;
double DENDRO_684 = -DENDRO_235*DENDRO_472;
double DENDRO_685 = DENDRO_250*DENDRO_76;
double DENDRO_686 = 0.5*DENDRO_463 + DENDRO_685;
double DENDRO_687 = DENDRO_160*DENDRO_476;
double DENDRO_688 = 0.25*DENDRO_131;
double DENDRO_689 = DENDRO_147*DENDRO_268;
double DENDRO_690 = DENDRO_259*DENDRO_688 + DENDRO_689;
double DENDRO_691 = DENDRO_221*DENDRO_545;
double DENDRO_692 = DENDRO_503 + DENDRO_691;
double DENDRO_693 = DENDRO_230*DENDRO_408;
double DENDRO_694 = DENDRO_131*DENDRO_160;
double DENDRO_695 = 1.0*DENDRO_314;
double DENDRO_696 = -DENDRO_282;
double DENDRO_697 = DENDRO_204*DENDRO_368;
double DENDRO_698 = DENDRO_105*DENDRO_167;
double DENDRO_699 = DENDRO_202*DENDRO_382;
double DENDRO_700 = DENDRO_105*DENDRO_204 + DENDRO_107*DENDRO_167 + DENDRO_237;
double DENDRO_701 = DENDRO_41*gt4[pp];
double DENDRO_702 = DENDRO_330*grad2_0_2_gt4[pp] + DENDRO_335*grad2_2_2_gt4[pp] + DENDRO_337*grad2_1_1_gt4[pp] + DENDRO_339*grad2_0_0_gt4[pp] - DENDRO_345*grad2_1_2_gt4[pp] - DENDRO_346*grad2_0_1_gt4[pp] + DENDRO_432*DENDRO_603 + DENDRO_499*DENDRO_600 + DENDRO_500*DENDRO_603 + DENDRO_501*DENDRO_605 - DENDRO_59*DENDRO_607 + DENDRO_596*gt1[pp] + DENDRO_598*gt3[pp];
double DENDRO_703 = -DENDRO_121*(DENDRO_590*(DENDRO_256 + DENDRO_697 + DENDRO_698) + DENDRO_592*(DENDRO_207 + DENDRO_699) + DENDRO_595*DENDRO_700 + DENDRO_696) + DENDRO_128*DENDRO_416 + DENDRO_215*DENDRO_414 + DENDRO_224*DENDRO_412 + DENDRO_419*DENDRO_701 + DENDRO_702;
double DENDRO_704 = grad2_1_2_alpha[pp];
double DENDRO_705 = DENDRO_256*DENDRO_41;
double DENDRO_706 = DENDRO_257*DENDRO_41;
double DENDRO_707 = DENDRO_258*DENDRO_41;
double DENDRO_708 = -DENDRO_59;
double DENDRO_709 = DENDRO_369*DENDRO_701;
double DENDRO_710 = DENDRO_64*(DENDRO_708 + DENDRO_709);
double DENDRO_711 = DENDRO_450*DENDRO_701;
double DENDRO_712 = DENDRO_205*DENDRO_41 + DENDRO_206*DENDRO_41;
double DENDRO_713 = -DENDRO_203*DENDRO_41 + DENDRO_712;
double DENDRO_714 = DENDRO_64*gt4[pp];
double DENDRO_715 = DENDRO_359*DENDRO_714;
double DENDRO_716 = DENDRO_41*DENDRO_631*(DENDRO_240 + DENDRO_715) + DENDRO_624*(DENDRO_705 - DENDRO_706 - DENDRO_707 + DENDRO_710) + DENDRO_632*(DENDRO_64*(DENDRO_628 + DENDRO_711) + DENDRO_713) - 4*DENDRO_704;
double DENDRO_717 = 0.5*DENDRO_398;
double DENDRO_718 = 1.0*DENDRO_490;
double DENDRO_719 = 0.5*DENDRO_489;
double DENDRO_720 = 0.25*DENDRO_694;
double DENDRO_721 = DENDRO_406 + DENDRO_569;
double DENDRO_722 = DENDRO_151*DENDRO_250;
double DENDRO_723 = DENDRO_681 + DENDRO_682;
double DENDRO_724 = DENDRO_259*DENDRO_404;
double DENDRO_725 = DENDRO_250*DENDRO_75;
double DENDRO_726 = DENDRO_235*DENDRO_304 + DENDRO_725;
double DENDRO_727 = DENDRO_160*DENDRO_437 + DENDRO_689;
double DENDRO_728 = DENDRO_502 + DENDRO_503;
double DENDRO_729 = DENDRO_230*DENDRO_470 + DENDRO_691;
double DENDRO_730 = 1.0*DENDRO_183;
double DENDRO_731 = -DENDRO_184*(DENDRO_642 + DENDRO_671) - DENDRO_302*(DENDRO_219*DENDRO_226 + DENDRO_433 + DENDRO_658) - DENDRO_695*(DENDRO_142 + DENDRO_532) - DENDRO_730*(DENDRO_435 + DENDRO_651);
double DENDRO_732 = DENDRO_567 + DENDRO_674;
double DENDRO_733 = 0.5*DENDRO_486;
double DENDRO_734 = DENDRO_127*DENDRO_230;
double DENDRO_735 = -DENDRO_734;
double DENDRO_736 = DENDRO_268*DENDRO_332;
double DENDRO_737 = DENDRO_155*DENDRO_437 + DENDRO_736;
double DENDRO_738 = DENDRO_242*DENDRO_304;
double DENDRO_739 = DENDRO_250*DENDRO_70;
double DENDRO_740 = DENDRO_119*DENDRO_71;
double DENDRO_741 = DENDRO_116*DENDRO_470 + DENDRO_155*DENDRO_174;
double DENDRO_742 = DENDRO_116*DENDRO_218;
double DENDRO_743 = 0.5*DENDRO_156;
double DENDRO_744 = -DENDRO_742 - DENDRO_743;
double DENDRO_745 = DENDRO_119*DENDRO_470 + DENDRO_517;
double DENDRO_746 = DENDRO_579 + DENDRO_666;
double DENDRO_747 = DENDRO_151*DENDRO_268;
double DENDRO_748 = 0.5*DENDRO_466;
double DENDRO_749 = -DENDRO_747 - DENDRO_748;
double DENDRO_750 = DENDRO_242*DENDRO_688 + DENDRO_725;
double DENDRO_751 = DENDRO_230*DENDRO_332;
double DENDRO_752 = DENDRO_502 + DENDRO_691;
double DENDRO_753 = DENDRO_268*DENDRO_75;
double DENDRO_754 = DENDRO_131*DENDRO_568 + DENDRO_753;
double DENDRO_755 = DENDRO_250*DENDRO_78;
double DENDRO_756 = DENDRO_242*DENDRO_327 + DENDRO_755;
double DENDRO_757 = DENDRO_155*DENDRO_561;
double DENDRO_758 = DENDRO_116*DENDRO_226 + DENDRO_167*DENDRO_568;
double DENDRO_759 = DENDRO_222*DENDRO_80;
double DENDRO_760 = DENDRO_141*DENDRO_304 + DENDRO_759;
double DENDRO_761 = 1.0*DENDRO_301;
double DENDRO_762 = -DENDRO_289;
double DENDRO_763 = DENDRO_105*DENDRO_69 + DENDRO_137*DENDRO_368 + DENDRO_152;
double DENDRO_764 = DENDRO_122*DENDRO_382;
double DENDRO_765 = DENDRO_107*DENDRO_69;
double DENDRO_766 = DENDRO_105*DENDRO_137;
double DENDRO_767 = 2.0*DENDRO_499*gt0[pp];
double DENDRO_768 = 2.0*gt1[pp];
double DENDRO_769 = DENDRO_322*DENDRO_768;
double DENDRO_770 = DENDRO_500*DENDRO_768;
double DENDRO_771 = DENDRO_501*DENDRO_600;
double DENDRO_772 = 2.0*DENDRO_323*gt3[pp];
double DENDRO_773 = DENDRO_325*DENDRO_603;
double DENDRO_774 = -DENDRO_290*DENDRO_321;
double DENDRO_775 = DENDRO_330*grad2_0_2_gt1[pp];
double DENDRO_776 = DENDRO_335*grad2_2_2_gt1[pp];
double DENDRO_777 = DENDRO_337*grad2_1_1_gt1[pp];
double DENDRO_778 = DENDRO_339*grad2_0_0_gt1[pp];
double DENDRO_779 = -DENDRO_345*grad2_1_2_gt1[pp];
double DENDRO_780 = -DENDRO_346*grad2_0_1_gt1[pp];
double DENDRO_781 = DENDRO_41*gt1[pp];
double DENDRO_782 = -DENDRO_121*(DENDRO_590*DENDRO_763 + DENDRO_592*(DENDRO_140 + DENDRO_764) + DENDRO_595*(DENDRO_170 + DENDRO_765 + DENDRO_766) + DENDRO_762) + DENDRO_124*DENDRO_412 + DENDRO_129*DENDRO_414 + DENDRO_416*DENDRO_67 + DENDRO_419*DENDRO_781 + DENDRO_767 + DENDRO_769 + DENDRO_770 + DENDRO_771 + DENDRO_772 + DENDRO_773 + DENDRO_774 + DENDRO_775 + DENDRO_776 + DENDRO_777 + DENDRO_778 + DENDRO_779 + DENDRO_780;
double DENDRO_783 = 0.25*DENDRO_340;
double DENDRO_784 = DENDRO_141*DENDRO_576 + DENDRO_653;
double DENDRO_785 = -DENDRO_141*DENDRO_472;
double DENDRO_786 = DENDRO_146*(DENDRO_540 + DENDRO_752) - DENDRO_184*(DENDRO_331 + DENDRO_654) - DENDRO_184*(DENDRO_331 + DENDRO_784) + DENDRO_190*(DENDRO_760 + DENDRO_785) - DENDRO_302*(DENDRO_134*DENDRO_226 + DENDRO_558) - DENDRO_315*(1.0*DENDRO_341 + DENDRO_783);
double DENDRO_787 = grad2_0_1_alpha[pp];
double DENDRO_788 = DENDRO_450*DENDRO_781;
double DENDRO_789 = DENDRO_138*DENDRO_41 + DENDRO_139*DENDRO_41;
double DENDRO_790 = -DENDRO_136*DENDRO_41 + DENDRO_789;
double DENDRO_791 = DENDRO_170*DENDRO_41;
double DENDRO_792 = DENDRO_171*DENDRO_41;
double DENDRO_793 = DENDRO_172*DENDRO_41;
double DENDRO_794 = DENDRO_108*DENDRO_781;
double DENDRO_795 = DENDRO_64*(DENDRO_708 + DENDRO_794);
double DENDRO_796 = DENDRO_64*gt1[pp];
double DENDRO_797 = DENDRO_796*DENDRO_88;
double DENDRO_798 = DENDRO_41*DENDRO_624*(DENDRO_155 + DENDRO_797) + DENDRO_631*(DENDRO_791 - DENDRO_792 - DENDRO_793 + DENDRO_795) + DENDRO_632*(DENDRO_64*(DENDRO_621 + DENDRO_788) + DENDRO_790) - 4*DENDRO_787;
double DENDRO_799 = DENDRO_180*DENDRO_242;
double DENDRO_800 = -DENDRO_29*(DENDRO_798 + alpha[pp]*(DENDRO_145*(DENDRO_464 + DENDRO_648 + DENDRO_657) + DENDRO_145*(DENDRO_583 + DENDRO_584 + DENDRO_656) + DENDRO_146*(DENDRO_684 + DENDRO_750) + DENDRO_146*(DENDRO_687 + DENDRO_749) + DENDRO_146*(DENDRO_751 + DENDRO_752) - DENDRO_184*(DENDRO_162 + DENDRO_744) - DENDRO_184*(DENDRO_516 + DENDRO_745) - DENDRO_184*(DENDRO_647 + DENDRO_746) - DENDRO_184*(DENDRO_574 + DENDRO_646 + DENDRO_720) + DENDRO_190*(DENDRO_757 + DENDRO_758) + DENDRO_190*(DENDRO_116*DENDRO_227 + DENDRO_754) + DENDRO_190*(-DENDRO_119*DENDRO_127 + DENDRO_756) + DENDRO_190*(DENDRO_230*DENDRO_70 + DENDRO_760) + DENDRO_200*(DENDRO_119*DENDRO_122 + DENDRO_525) - DENDRO_302*(DENDRO_405 + DENDRO_732) - DENDRO_306*(DENDRO_482 + DENDRO_737) - DENDRO_306*(DENDRO_733 + DENDRO_735) - DENDRO_306*(DENDRO_484 + DENDRO_738 + DENDRO_739) - DENDRO_315*(0.5*DENDRO_522 + DENDRO_741) - DENDRO_315*(DENDRO_119*DENDRO_70 + DENDRO_520 + DENDRO_740) - DENDRO_761*(DENDRO_402 + DENDRO_636 + DENDRO_669) + DENDRO_782 + DENDRO_786)) - DENDRO_29*(DENDRO_798 + alpha[pp]*(DENDRO_146*(DENDRO_468 + DENDRO_749) + DENDRO_146*(DENDRO_471 + DENDRO_726) + DENDRO_146*(DENDRO_471 + DENDRO_750) + DENDRO_146*(DENDRO_540 + DENDRO_729) + DENDRO_146*(DENDRO_639 + DENDRO_727) + DENDRO_146*(DENDRO_728 + DENDRO_751) - DENDRO_184*(DENDRO_517 + DENDRO_746) - DENDRO_184*(DENDRO_572 + DENDRO_784) - DENDRO_184*(DENDRO_575 + DENDRO_744) - DENDRO_184*(DENDRO_579 + DENDRO_745) + DENDRO_190*(DENDRO_754 + DENDRO_757) + DENDRO_190*(DENDRO_756 + DENDRO_799) + DENDRO_190*(DENDRO_268*DENDRO_76 + DENDRO_758) + DENDRO_190*(DENDRO_119*DENDRO_123 + DENDRO_755 + DENDRO_799) + DENDRO_190*(DENDRO_230*DENDRO_71 + DENDRO_759 + DENDRO_785) + DENDRO_200*(DENDRO_230*DENDRO_69 + DENDRO_505) - DENDRO_302*(DENDRO_406 + DENDRO_732) - DENDRO_302*(DENDRO_240*DENDRO_75 + DENDRO_670) - DENDRO_306*(0.5*DENDRO_491 + DENDRO_737) - DENDRO_306*(1.0*DENDRO_492 + DENDRO_738) - DENDRO_306*(DENDRO_123*DENDRO_230 + DENDRO_487 + DENDRO_735) - DENDRO_315*(DENDRO_518 + DENDRO_741) - DENDRO_315*(0.5*DENDRO_519 + DENDRO_740) - DENDRO_315*(DENDRO_123*DENDRO_80 + DENDRO_333 + DENDRO_783) - DENDRO_730*(DENDRO_142 + DENDRO_529 + DENDRO_531) - DENDRO_730*(DENDRO_526 + DENDRO_527 + DENDRO_694) - DENDRO_761*(DENDRO_435 + DENDRO_557 + DENDRO_650) + DENDRO_782)) - DENDRO_31*(DENDRO_716 + alpha[pp]*(DENDRO_146*(DENDRO_676 + DENDRO_724) + DENDRO_146*(DENDRO_678 + DENDRO_679) + DENDRO_146*(DENDRO_680 - DENDRO_722) + DENDRO_146*(-DENDRO_218*DENDRO_230 + DENDRO_723) + DENDRO_146*(DENDRO_226*DENDRO_265 + DENDRO_675 + DENDRO_724) - DENDRO_184*(DENDRO_643 + DENDRO_668) - DENDRO_184*(DENDRO_672 + DENDRO_721) - DENDRO_184*(DENDRO_674 + DENDRO_721) + DENDRO_190*(DENDRO_465 + DENDRO_686) + DENDRO_190*(DENDRO_469 + DENDRO_690) + DENDRO_190*(DENDRO_469 + DENDRO_727) + DENDRO_190*(DENDRO_538 + DENDRO_726) + DENDRO_190*(DENDRO_539 + DENDRO_729) + DENDRO_190*(DENDRO_693 + DENDRO_728) - DENDRO_302*(DENDRO_396 + DENDRO_661) - DENDRO_302*(DENDRO_660 + DENDRO_717) - DENDRO_306*(DENDRO_663 + DENDRO_718) - DENDRO_306*(DENDRO_665 + DENDRO_719) - DENDRO_306*(DENDRO_226*DENDRO_230 + DENDRO_480 + DENDRO_662) - DENDRO_315*(DENDRO_517 + DENDRO_667) - DENDRO_315*(DENDRO_147*DENDRO_155 + DENDRO_720) + DENDRO_410*(DENDRO_204*DENDRO_230 + DENDRO_504) + DENDRO_703 - DENDRO_730*(DENDRO_637 + DENDRO_669) + DENDRO_731)) - DENDRO_31*(DENDRO_716 + alpha[pp]*(DENDRO_146*(DENDRO_679 + DENDRO_680) + DENDRO_146*(DENDRO_681 + DENDRO_683) + DENDRO_146*(-DENDRO_127*DENDRO_247 + DENDRO_678) + DENDRO_146*(DENDRO_214*DENDRO_230 + DENDRO_683) + DENDRO_146*(DENDRO_227*DENDRO_265 + DENDRO_676) - DENDRO_184*(DENDRO_403 + DENDRO_668) - DENDRO_184*(DENDRO_438 + DENDRO_560) - DENDRO_184*(DENDRO_438 + DENDRO_671) - DENDRO_184*(DENDRO_565 + DENDRO_672) - DENDRO_184*(DENDRO_570 + DENDRO_674) - DENDRO_184*(DENDRO_571 + DENDRO_670) + DENDRO_189*(DENDRO_585 + DENDRO_656) + DENDRO_189*(DENDRO_649 + DENDRO_657) + DENDRO_190*(DENDRO_539 + DENDRO_692) + DENDRO_190*(DENDRO_684 + DENDRO_686) + DENDRO_190*(DENDRO_687 + DENDRO_690) + DENDRO_190*(DENDRO_692 + DENDRO_693) - DENDRO_302*(0.5*DENDRO_409 + DENDRO_661) - DENDRO_302*(1.0*DENDRO_440 + DENDRO_658) - DENDRO_302*(DENDRO_214*DENDRO_265 + DENDRO_399 + DENDRO_660) - DENDRO_306*(DENDRO_474 + DENDRO_665) - DENDRO_306*(0.5*DENDRO_479 + DENDRO_662) - DENDRO_306*(DENDRO_214*DENDRO_268 + DENDRO_477 + DENDRO_663) - DENDRO_315*(DENDRO_516 + DENDRO_667) - DENDRO_315*(DENDRO_123*DENDRO_134 + DENDRO_652) + DENDRO_410*(DENDRO_202*DENDRO_265 + DENDRO_425) - DENDRO_695*(DENDRO_528 + DENDRO_694) + DENDRO_703)) + DENDRO_35*(DENDRO_355*(DENDRO_455 + DENDRO_64*(DENDRO_449 + DENDRO_453)) + DENDRO_362*(DENDRO_250 + DENDRO_359*DENDRO_446) - 4*DENDRO_445 + DENDRO_447*(DENDRO_268 + DENDRO_446*DENDRO_88) + alpha[pp]*(-DENDRO_121*(DENDRO_117*DENDRO_459 + DENDRO_118*(DENDRO_228 + DENDRO_460 + DENDRO_461) + DENDRO_120*DENDRO_462 + DENDRO_456) + DENDRO_122*DENDRO_416 + DENDRO_146*(DENDRO_474 + DENDRO_475) + DENDRO_146*(DENDRO_477 - DENDRO_478) + DENDRO_146*(DENDRO_480 + 1.0*DENDRO_481) - DENDRO_155*DENDRO_498 - DENDRO_184*(DENDRO_463 + DENDRO_465) - DENDRO_184*(-1.0*DENDRO_466 + DENDRO_468) - DENDRO_184*(DENDRO_242*DENDRO_470 + DENDRO_471) - DENDRO_184*(DENDRO_259*DENDRO_470 + DENDRO_469) + DENDRO_190*(DENDRO_482 + DENDRO_483) + DENDRO_190*(DENDRO_484 + DENDRO_485) + DENDRO_190*(DENDRO_487 + 1.0*DENDRO_488) + DENDRO_200*(DENDRO_486 + DENDRO_488) + DENDRO_200*(DENDRO_122*DENDRO_242 + DENDRO_492) + DENDRO_200*(DENDRO_155*DENDRO_202 + DENDRO_491) + DENDRO_202*DENDRO_414 - DENDRO_221*DENDRO_231*DENDRO_428 + DENDRO_221*DENDRO_412 - DENDRO_240*DENDRO_495 - DENDRO_242*DENDRO_497 + DENDRO_410*(DENDRO_479 + DENDRO_481) + DENDRO_410*(DENDRO_122*DENDRO_240 + DENDRO_489) + DENDRO_410*(DENDRO_202*DENDRO_259 + DENDRO_490) + DENDRO_419*DENDRO_451 - DENDRO_493*(DENDRO_216 - DENDRO_404) - DENDRO_496*(DENDRO_123 + DENDRO_125) + DENDRO_506)) + DENDRO_46*(DENDRO_635 + alpha[pp]*(DENDRO_145*(DENDRO_402 + DENDRO_637) + DENDRO_146*(DENDRO_556 + DENDRO_642) + DENDRO_146*(DENDRO_563 + DENDRO_643) + DENDRO_146*(DENDRO_564 + DENDRO_644) + DENDRO_146*(DENDRO_567 + DENDRO_644) - DENDRO_184*(DENDRO_548 + DENDRO_640) - DENDRO_184*(-DENDRO_119*DENDRO_151 + DENDRO_553) - DENDRO_184*(DENDRO_265*DENDRO_75 + DENDRO_640) + DENDRO_190*(DENDRO_517 + DENDRO_582) + DENDRO_190*(DENDRO_577 + DENDRO_645) + DENDRO_190*(DENDRO_581 + DENDRO_647) + DENDRO_190*(DENDRO_645 + DENDRO_646) - DENDRO_197*(DENDRO_119*DENDRO_132 + DENDRO_524) - DENDRO_302*(DENDRO_535 + DENDRO_638) - DENDRO_302*(DENDRO_247*DENDRO_75 + DENDRO_390 + DENDRO_533) - DENDRO_306*(DENDRO_155*DENDRO_214 + DENDRO_639) - DENDRO_315*(1.0*DENDRO_523 + DENDRO_543) - DENDRO_315*(DENDRO_119*DENDRO_75 + DENDRO_513 + DENDRO_542) - DENDRO_586*(DENDRO_464 + DENDRO_649) + DENDRO_616 + DENDRO_655)) + DENDRO_46*(DENDRO_635 + alpha[pp]*(DENDRO_146*(DENDRO_407 + DENDRO_563) + DENDRO_146*(DENDRO_407 + DENDRO_571) + DENDRO_146*(DENDRO_436 + DENDRO_556) + DENDRO_146*(DENDRO_558 + DENDRO_560) + DENDRO_146*(DENDRO_564 + DENDRO_565) + DENDRO_146*(DENDRO_567 + DENDRO_570) - DENDRO_184*(DENDRO_549 + DENDRO_550) - DENDRO_184*(DENDRO_551 + DENDRO_553) - DENDRO_184*(DENDRO_219*DENDRO_71 + DENDRO_546) - DENDRO_184*(DENDRO_119*DENDRO_147 + DENDRO_551 + DENDRO_552) - DENDRO_184*(DENDRO_265*DENDRO_76 + DENDRO_547 + DENDRO_548) + DENDRO_189*(DENDRO_526 + DENDRO_528) + DENDRO_189*(DENDRO_529 + DENDRO_532) + DENDRO_190*(DENDRO_516 + DENDRO_581) + DENDRO_190*(DENDRO_572 + DENDRO_573) + DENDRO_190*(DENDRO_579 + DENDRO_582) + DENDRO_190*(DENDRO_574 + DENDRO_575 + DENDRO_577) - DENDRO_197*(DENDRO_265*DENDRO_74 + DENDRO_426) - DENDRO_302*(1.0*DENDRO_417 + DENDRO_533) - DENDRO_302*(0.5*DENDRO_442 + DENDRO_536) - DENDRO_302*(DENDRO_147*DENDRO_265 + DENDRO_393 + DENDRO_535) - DENDRO_306*(DENDRO_503 + DENDRO_541) - DENDRO_306*(DENDRO_240*DENDRO_70 + DENDRO_538) - DENDRO_315*(DENDRO_328 + DENDRO_544) - DENDRO_315*(0.5*DENDRO_512 + DENDRO_542) - DENDRO_315*(DENDRO_116*DENDRO_147 + DENDRO_511 + DENDRO_543) - DENDRO_586*(DENDRO_583 + DENDRO_585) + DENDRO_616)) + DENDRO_49*(DENDRO_356*DENDRO_81 + DENDRO_361*(-DENDRO_100 + DENDRO_102 + DENDRO_111 - DENDRO_98) + DENDRO_447*(DENDRO_116 + DENDRO_89) - 4*DENDRO_58 + alpha[pp]*(-DENDRO_121*(DENDRO_115 + DENDRO_117*DENDRO_507 + DENDRO_118*DENDRO_508 + DENDRO_120*(DENDRO_101 + DENDRO_509 + DENDRO_510)) - DENDRO_144*DENDRO_252*DENDRO_320 + DENDRO_146*(-1.0*DENDRO_135 + DENDRO_143) + DENDRO_146*(-1.0*DENDRO_156 + DENDRO_162) + DENDRO_146*(DENDRO_160*DENDRO_332 + DENDRO_515) + DENDRO_146*(DENDRO_235*DENDRO_70 + DENDRO_517) + DENDRO_146*(DENDRO_242*DENDRO_75 + DENDRO_516) - DENDRO_155*DENDRO_312 - DENDRO_184*(DENDRO_513 + 1.0*DENDRO_514) - DENDRO_184*(-DENDRO_116*DENDRO_182 + DENDRO_511) + DENDRO_190*(DENDRO_520 + 1.0*DENDRO_521) + DENDRO_190*(DENDRO_116*DENDRO_167 + DENDRO_518) - DENDRO_197*(DENDRO_512 + DENDRO_514) - DENDRO_197*(DENDRO_160*DENDRO_74 + DENDRO_523) + DENDRO_200*(DENDRO_519 + DENDRO_521) + DENDRO_200*(DENDRO_155*DENDRO_74 + DENDRO_522) + DENDRO_234*DENDRO_411 + DENDRO_254*DENDRO_415 - DENDRO_270*DENDRO_313*DENDRO_349 + DENDRO_272*DENDRO_413 + DENDRO_300*DENDRO_418 - DENDRO_303*(DENDRO_149 - DENDRO_165) - DENDRO_307*(DENDRO_125 - DENDRO_304) - DENDRO_317*DENDRO_524 - DENDRO_319*DENDRO_525 + DENDRO_350)) + DENDRO_55*(-4*DENDRO_352 + DENDRO_354*DENDRO_356 + DENDRO_362*(DENDRO_247 + DENDRO_360) + DENDRO_373*(-DENDRO_363 + DENDRO_364 - DENDRO_365 + DENDRO_372) + alpha[pp]*(-DENDRO_121*(DENDRO_117*(DENDRO_375 + DENDRO_378 + DENDRO_381) + DENDRO_118*DENDRO_383 + DENDRO_120*DENDRO_384 + DENDRO_374) + DENDRO_132*DENDRO_416 + DENDRO_146*(DENDRO_396 + DENDRO_397) + DENDRO_146*(DENDRO_399 + 1.0*DENDRO_400) - DENDRO_184*(DENDRO_387 + DENDRO_388) - DENDRO_184*(DENDRO_393 + 1.0*DENDRO_394) - DENDRO_184*(DENDRO_247*DENDRO_391 + DENDRO_390) + DENDRO_190*(DENDRO_401 + DENDRO_403) + DENDRO_190*(DENDRO_147*DENDRO_259 + DENDRO_405) + DENDRO_190*(DENDRO_160*DENDRO_214 + DENDRO_406) + DENDRO_190*(DENDRO_235*DENDRO_408 + DENDRO_407) - DENDRO_197*(DENDRO_392 + DENDRO_394) - DENDRO_197*(DENDRO_132*DENDRO_235 + DENDRO_417) + DENDRO_204*DENDRO_412 - DENDRO_211*DENDRO_266*DENDRO_428 + DENDRO_211*DENDRO_414 - DENDRO_235*DENDRO_424 - DENDRO_240*DENDRO_423 - DENDRO_319*DENDRO_425 + DENDRO_367*DENDRO_419 + DENDRO_410*(DENDRO_398 + DENDRO_400) + DENDRO_410*(DENDRO_132*DENDRO_240 + DENDRO_409) - DENDRO_420*(DENDRO_147 + DENDRO_149) - DENDRO_421*(DENDRO_214 + DENDRO_216) - DENDRO_426*DENDRO_427 + DENDRO_444));
double DENDRO_801 = DENDRO_41*DENDRO_800;
double DENDRO_802 = (1.0/12.0)*chi[pp];
double DENDRO_803 = (1.0/3.0)*At1[pp];
double DENDRO_804 = At4[pp]*DENDRO_31;
double DENDRO_805 = At3[pp]*DENDRO_35;
double DENDRO_806 = DENDRO_43 + DENDRO_804 - DENDRO_805;
double DENDRO_807 = At3[pp]*DENDRO_29;
double DENDRO_808 = At4[pp]*DENDRO_46;
double DENDRO_809 = -At1[pp]*DENDRO_49 + DENDRO_807 - DENDRO_808;
double DENDRO_810 = At4[pp]*DENDRO_55;
double DENDRO_811 = -At1[pp]*DENDRO_46 + At3[pp]*DENDRO_31 - DENDRO_810;
double DENDRO_812 = 6.0*DENDRO_94;
double DENDRO_813 = 6.0*DENDRO_112;
double DENDRO_814 = 6.0*DENDRO_82;
double DENDRO_815 = DENDRO_131*DENDRO_179;
double DENDRO_816 = DENDRO_167*DENDRO_179;
double DENDRO_817 = -DENDRO_237 + DENDRO_238 + DENDRO_239;
double DENDRO_818 = DENDRO_74*DENDRO_817;
double DENDRO_819 = DENDRO_816 + DENDRO_818;
double DENDRO_820 = DENDRO_173*DENDRO_472;
double DENDRO_821 = -DENDRO_164*DENDRO_476;
double DENDRO_822 = DENDRO_223 - DENDRO_228 + DENDRO_229;
double DENDRO_823 = DENDRO_168*DENDRO_688;
double DENDRO_824 = -DENDRO_256 + DENDRO_257 + DENDRO_258;
double DENDRO_825 = DENDRO_214*DENDRO_93;
double DENDRO_826 = DENDRO_174*DENDRO_824 + DENDRO_825;
double DENDRO_827 = DENDRO_179*DENDRO_327;
double DENDRO_828 = DENDRO_187*DENDRO_332 + 0.25*DENDRO_77*DENDRO_817;
double DENDRO_829 = DENDRO_175 + DENDRO_181;
double DENDRO_830 = DENDRO_167*DENDRO_824;
double DENDRO_831 = 1.0*DENDRO_145;
double DENDRO_832 = DENDRO_167*DENDRO_173;
double DENDRO_833 = DENDRO_69*DENDRO_817 + DENDRO_832;
double DENDRO_834 = 2.0*DENDRO_233;
double DENDRO_835 = 2.0*DENDRO_253;
double DENDRO_836 = 2.0*DENDRO_271;
double DENDRO_837 = DENDRO_299*DENDRO_64;
double DENDRO_838 = (1.0/3.0)*At2[pp];
double DENDRO_839 = At5[pp]*DENDRO_31;
double DENDRO_840 = At2[pp]*DENDRO_29 - At4[pp]*DENDRO_35 + DENDRO_839;
double DENDRO_841 = At5[pp]*DENDRO_46;
double DENDRO_842 = -At2[pp]*DENDRO_49 + At4[pp]*DENDRO_29 - DENDRO_841;
double DENDRO_843 = -At5[pp]*DENDRO_55 + DENDRO_48 + DENDRO_804;
double DENDRO_844 = DENDRO_137*DENDRO_173;
double DENDRO_845 = DENDRO_179*DENDRO_389;
double DENDRO_846 = DENDRO_244 + DENDRO_245 - DENDRO_246;
double DENDRO_847 = 0.25*DENDRO_164*DENDRO_211;
double DENDRO_848 = DENDRO_168*DENDRO_404;
double DENDRO_849 = DENDRO_262 - DENDRO_263 + DENDRO_264;
double DENDRO_850 = DENDRO_165*DENDRO_168;
double DENDRO_851 = DENDRO_212*DENDRO_93;
double DENDRO_852 = DENDRO_173*DENDRO_191;
double DENDRO_853 = -DENDRO_847;
double DENDRO_854 = DENDRO_173*DENDRO_389;
double DENDRO_855 = DENDRO_132*DENDRO_173;
double DENDRO_856 = 2*At4[pp];
double DENDRO_857 = At3[pp]*DENDRO_51;
double DENDRO_858 = DENDRO_41*DENDRO_856;
double DENDRO_859 = DENDRO_113*DENDRO_41;
double DENDRO_860 = DENDRO_202*DENDRO_822;
double DENDRO_861 = 1.0*DENDRO_817;
double DENDRO_862 = 0.25*DENDRO_832;
double DENDRO_863 = DENDRO_164*DENDRO_404;
double DENDRO_864 = DENDRO_122*DENDRO_822;
double DENDRO_865 = (1.0/3.0)*At4[pp];
double DENDRO_866 = DENDRO_165*DENDRO_824;
double DENDRO_867 = -DENDRO_404*DENDRO_824 + DENDRO_675;
double DENDRO_868 = DENDRO_689 - DENDRO_863;
double DENDRO_869 = DENDRO_211*DENDRO_824;
double DENDRO_870 = DENDRO_204*DENDRO_849;
double DENDRO_871 = DENDRO_168*DENDRO_211;
double DENDRO_872 = DENDRO_132*DENDRO_849;
double DENDRO_873 = DENDRO_41*DENDRO_82;
double DENDRO_874 = DENDRO_109*DENDRO_64;
double DENDRO_875 = DENDRO_112*DENDRO_41;
double DENDRO_876 = DENDRO_370*DENDRO_64;
double DENDRO_877 = DENDRO_41*DENDRO_94;
double DENDRO_878 = 0.5*DENDRO_94;
double DENDRO_879 = 0.5*DENDRO_82;
double DENDRO_880 = DENDRO_284*DENDRO_41;
double DENDRO_881 = 0.5*DENDRO_112;
double DENDRO_882 = DENDRO_287*DENDRO_41;
double DENDRO_883 = DENDRO_291*DENDRO_41;
double DENDRO_884 = pow(DENDRO_29, 2);
double DENDRO_885 = pow(DENDRO_46, 2);
double DENDRO_886 = 2*DENDRO_49;
double DENDRO_887 = At0[pp]*pow(DENDRO_49, 2) + At3[pp]*DENDRO_884 + At5[pp]*DENDRO_885 - DENDRO_291*DENDRO_808 - DENDRO_43*DENDRO_886 + DENDRO_47*DENDRO_886;
double DENDRO_888 = 3*DENDRO_144;
double DENDRO_889 = pow(DENDRO_31, 2);
double DENDRO_890 = 2*DENDRO_35;
double DENDRO_891 = At0[pp]*DENDRO_884 + At3[pp]*pow(DENDRO_35, 2) + At5[pp]*DENDRO_889 + DENDRO_291*DENDRO_32 - DENDRO_43*DENDRO_890 - DENDRO_804*DENDRO_890;
double DENDRO_892 = 2*DENDRO_55;
double DENDRO_893 = At0[pp]*DENDRO_885 + At3[pp]*DENDRO_889 + At5[pp]*pow(DENDRO_55, 2) - DENDRO_287*DENDRO_53 + DENDRO_47*DENDRO_892 - DENDRO_804*DENDRO_892;
double DENDRO_894 = At2[pp]*DENDRO_885 - DENDRO_29*DENDRO_810 + DENDRO_31*DENDRO_807 - DENDRO_31*DENDRO_808 - DENDRO_43*DENDRO_46 + DENDRO_46*DENDRO_50 - DENDRO_49*DENDRO_53 + DENDRO_49*DENDRO_56 + DENDRO_55*DENDRO_841;
double DENDRO_895 = 6*DENDRO_144;
double DENDRO_896 = At1[pp]*DENDRO_884;
double DENDRO_897 = DENDRO_29*DENDRO_804;
double DENDRO_898 = DENDRO_29*DENDRO_47;
double DENDRO_899 = DENDRO_35*DENDRO_808;
double DENDRO_900 = DENDRO_31*DENDRO_841;
double DENDRO_901 = DENDRO_29*DENDRO_50;
double DENDRO_902 = DENDRO_36*DENDRO_49;
double DENDRO_903 = DENDRO_32*DENDRO_49;
double DENDRO_904 = DENDRO_35*DENDRO_807;
double DENDRO_905 = DENDRO_896 + DENDRO_897 - DENDRO_898 + DENDRO_899 - DENDRO_900 - DENDRO_901 + DENDRO_902 - DENDRO_903 - DENDRO_904;
double DENDRO_906 = At4[pp]*DENDRO_889;
double DENDRO_907 = DENDRO_31*DENDRO_43;
double DENDRO_908 = DENDRO_30*DENDRO_46;
double DENDRO_909 = DENDRO_36*DENDRO_46;
double DENDRO_910 = DENDRO_31*DENDRO_47;
double DENDRO_911 = DENDRO_29*DENDRO_56;
double DENDRO_912 = DENDRO_31*DENDRO_805;
double DENDRO_913 = DENDRO_35*DENDRO_810;
double DENDRO_914 = DENDRO_55*DENDRO_839;
double DENDRO_915 = DENDRO_906 + DENDRO_907 - DENDRO_908 + DENDRO_909 - DENDRO_910 - DENDRO_911 - DENDRO_912 + DENDRO_913 - DENDRO_914;
double DENDRO_916 = (1.0/3.0)*alpha[pp];
double DENDRO_917 = 2*DENDRO_144;
double DENDRO_918 = DENDRO_112*DENDRO_917;
double DENDRO_919 = DENDRO_82*DENDRO_917;
double DENDRO_920 = DENDRO_917*DENDRO_94;
double DENDRO_921 = grad2_0_2_beta0[pp];
double DENDRO_922 = (7.0/3.0)*DENDRO_329;
double DENDRO_923 = grad2_0_0_beta0[pp];
double DENDRO_924 = 4*grad_0_K[pp];
double DENDRO_925 = 9*DENDRO_64;
double DENDRO_926 = DENDRO_120*DENDRO_925;
double DENDRO_927 = DENDRO_41*DENDRO_916;
double DENDRO_928 = 4*grad_2_K[pp];
double DENDRO_929 = DENDRO_117*DENDRO_925;
double DENDRO_930 = 4*grad_1_K[pp];
double DENDRO_931 = DENDRO_118*DENDRO_925;
double DENDRO_932 = grad2_0_1_beta1[pp];
double DENDRO_933 = (1.0/3.0)*DENDRO_338;
double DENDRO_934 = grad2_0_2_beta2[pp];
double DENDRO_935 = grad2_1_2_beta1[pp];
double DENDRO_936 = (1.0/3.0)*DENDRO_329;
double DENDRO_937 = grad2_2_2_beta2[pp];
double DENDRO_938 = grad2_1_1_beta1[pp];
double DENDRO_939 = DENDRO_29*DENDRO_41;
double DENDRO_940 = (1.0/3.0)*DENDRO_939;
double DENDRO_941 = grad2_1_2_beta2[pp];
double DENDRO_942 = (2.0/3.0)*DENDRO_25;
double DENDRO_943 = grad2_0_1_beta0[pp];
double DENDRO_944 = (7.0/3.0)*DENDRO_939;
double DENDRO_945 = pow(DENDRO_40, -3);
double DENDRO_946 = DENDRO_0*DENDRO_945;
double DENDRO_947 = DENDRO_891*DENDRO_946;
double DENDRO_948 = DENDRO_893*DENDRO_946;
double DENDRO_949 = DENDRO_887*DENDRO_946;
double DENDRO_950 = 2.0*DENDRO_945*alpha[pp];
double DENDRO_951 = DENDRO_894*DENDRO_950;
double DENDRO_952 = DENDRO_905*DENDRO_950;
double DENDRO_953 = DENDRO_915*DENDRO_950;
double DENDRO_954 = beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp];
double DENDRO_955 = -DENDRO_11*DENDRO_411 + DENDRO_119*DENDRO_949 - DENDRO_15*DENDRO_413 - DENDRO_2*DENDRO_415 + DENDRO_235*DENDRO_951 + DENDRO_240*DENDRO_953 + DENDRO_242*DENDRO_952 + DENDRO_247*DENDRO_948 + DENDRO_250*DENDRO_947 - DENDRO_334*grad2_2_2_beta0[pp] - DENDRO_336*grad2_1_1_beta0[pp] - 4.0/3.0*DENDRO_338*DENDRO_923 + DENDRO_415*DENDRO_942 + DENDRO_880*grad2_1_2_beta0[pp] - DENDRO_887*DENDRO_918 - DENDRO_894*DENDRO_920 - DENDRO_905*DENDRO_919 - DENDRO_921*DENDRO_922 - DENDRO_927*(DENDRO_29*DENDRO_930 + DENDRO_905*DENDRO_931) - DENDRO_927*(-DENDRO_46*DENDRO_928 + DENDRO_894*DENDRO_929) - DENDRO_927*(-DENDRO_49*DENDRO_924 + DENDRO_887*DENDRO_926) - DENDRO_932*DENDRO_933 - DENDRO_933*DENDRO_934 - DENDRO_935*DENDRO_936 - DENDRO_936*DENDRO_937 + DENDRO_938*DENDRO_940 + DENDRO_940*DENDRO_941 + DENDRO_943*DENDRO_944 + DENDRO_954;
double DENDRO_956 = -DENDRO_896 - DENDRO_897 + DENDRO_898 - DENDRO_899 + DENDRO_900 + DENDRO_901 - DENDRO_902 + DENDRO_903 + DENDRO_904;
double DENDRO_957 = -DENDRO_906 - DENDRO_907 + DENDRO_908 - DENDRO_909 + DENDRO_910 + DENDRO_911 + DENDRO_912 - DENDRO_913 + DENDRO_914;
double DENDRO_958 = DENDRO_31*DENDRO_928;
double DENDRO_959 = DENDRO_29*DENDRO_924;
double DENDRO_960 = DENDRO_35*DENDRO_930;
double DENDRO_961 = DENDRO_891*DENDRO_931;
double DENDRO_962 = (1.0/3.0)*DENDRO_336;
double DENDRO_963 = DENDRO_31*DENDRO_41;
double DENDRO_964 = (1.0/3.0)*DENDRO_921;
double DENDRO_965 = (1.0/3.0)*DENDRO_963;
double DENDRO_966 = (7.0/3.0)*DENDRO_963;
double DENDRO_967 = beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp];
double DENDRO_968 = DENDRO_134*DENDRO_951 + DENDRO_219*DENDRO_948 - DENDRO_334*grad2_2_2_beta1[pp] - 4.0/3.0*DENDRO_336*DENDRO_938 - DENDRO_338*grad2_0_0_beta1[pp] + DENDRO_80*DENDRO_949 - DENDRO_882*grad2_0_2_beta1[pp] - DENDRO_891*DENDRO_919 + DENDRO_923*DENDRO_940 + DENDRO_932*DENDRO_944 + DENDRO_934*DENDRO_940 + DENDRO_935*DENDRO_966 + DENDRO_937*DENDRO_965 - DENDRO_941*DENDRO_962 - DENDRO_943*DENDRO_962 + DENDRO_963*DENDRO_964 + DENDRO_967;
double DENDRO_969 = beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp];
double DENDRO_970 = DENDRO_116*DENDRO_949 - DENDRO_12*DENDRO_411 + DENDRO_155*DENDRO_952 + DENDRO_160*DENDRO_951 + DENDRO_259*DENDRO_953 + DENDRO_265*DENDRO_948 + DENDRO_268*DENDRO_947 - 1.0/3.0*DENDRO_334*DENDRO_935 - 4.0/3.0*DENDRO_334*DENDRO_937 - DENDRO_334*DENDRO_964 - DENDRO_336*grad2_1_1_beta2[pp] - DENDRO_338*grad2_0_0_beta2[pp] - DENDRO_413*DENDRO_6 + DENDRO_413*DENDRO_942 - DENDRO_415*DENDRO_9 + DENDRO_883*grad2_0_1_beta2[pp] - DENDRO_893*DENDRO_920 - DENDRO_894*DENDRO_918 - DENDRO_915*DENDRO_919 - DENDRO_922*DENDRO_934 - DENDRO_923*DENDRO_936 - DENDRO_927*(DENDRO_31*DENDRO_930 + DENDRO_915*DENDRO_931) - DENDRO_927*(-DENDRO_46*DENDRO_924 + DENDRO_894*DENDRO_926) - DENDRO_927*(-DENDRO_55*DENDRO_928 + DENDRO_893*DENDRO_929) - DENDRO_932*DENDRO_936 + DENDRO_938*DENDRO_965 + DENDRO_941*DENDRO_966 + DENDRO_943*DENDRO_965 + DENDRO_969;
double DENDRO_971 = 2*DENDRO_61;
double DENDRO_972 = BSSN_ETA_R0*sqrt(DENDRO_41*(DENDRO_106*DENDRO_971 - DENDRO_274*DENDRO_55 - DENDRO_278*DENDRO_35 - DENDRO_280*DENDRO_49 - DENDRO_357*DENDRO_971 + 2*DENDRO_60*DENDRO_85))*pow(-pow(chi[pp], BSSN_ETA_POWER[0]) + 1, -BSSN_ETA_POWER[1]);

// Dendro: printing variables
//--
a_rhs[pp] = -DENDRO_0*K[pp] + lambda[0]*(beta0[pp]*agrad_0_alpha[pp] + beta1[pp]*agrad_1_alpha[pp] + beta2[pp]*agrad_2_alpha[pp]);
//--
b_rhs0[pp] = B0[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta0[pp] + beta1[pp]*agrad_1_beta0[pp] + beta2[pp]*agrad_2_beta0[pp]);
//--
b_rhs1[pp] = B1[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta1[pp] + beta1[pp]*agrad_1_beta1[pp] + beta2[pp]*agrad_2_beta1[pp]);
//--
b_rhs2[pp] = B2[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta2[pp] + beta1[pp]*agrad_1_beta2[pp] + beta2[pp]*agrad_2_beta2[pp]);
//--
gt_rhs00[pp] = -At0[pp]*DENDRO_0 + DENDRO_10*DENDRO_9 + DENDRO_3*gt0[pp] - DENDRO_4*DENDRO_5 - DENDRO_5*DENDRO_6 + DENDRO_7*DENDRO_8 + beta0[pp]*agrad_0_gt0[pp] + beta1[pp]*agrad_1_gt0[pp] + beta2[pp]*agrad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp]*DENDRO_0 + DENDRO_11*gt0[pp] + DENDRO_12*gt2[pp] + DENDRO_13*DENDRO_2 + DENDRO_13*DENDRO_4 - DENDRO_14*gt1[pp] + DENDRO_7*gt3[pp] + DENDRO_9*gt4[pp] + beta0[pp]*agrad_0_gt1[pp] + beta1[pp]*agrad_1_gt1[pp] + beta2[pp]*agrad_2_gt1[pp];
//--
gt_rhs02[pp] = -At2[pp]*DENDRO_0 + DENDRO_15*gt0[pp] + DENDRO_16*gt1[pp] + DENDRO_17*DENDRO_2 + DENDRO_17*DENDRO_6 - DENDRO_18*gt2[pp] + DENDRO_7*gt4[pp] + DENDRO_9*gt5[pp] + beta0[pp]*agrad_0_gt2[pp] + beta1[pp]*agrad_1_gt2[pp] + beta2[pp]*agrad_2_gt2[pp];
//--
gt_rhs11[pp] = -At3[pp]*DENDRO_0 + DENDRO_11*DENDRO_8 + DENDRO_12*DENDRO_21 - DENDRO_14*gt3[pp] - DENDRO_19*gt3[pp] + DENDRO_20*gt3[pp] + beta0[pp]*agrad_0_gt3[pp] + beta1[pp]*agrad_1_gt3[pp] + beta2[pp]*agrad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp]*DENDRO_0 + DENDRO_11*gt2[pp] + DENDRO_12*gt5[pp] + DENDRO_15*gt1[pp] + DENDRO_16*gt3[pp] - DENDRO_19*gt4[pp] + DENDRO_22*DENDRO_4 + DENDRO_22*DENDRO_6 + beta0[pp]*agrad_0_gt4[pp] + beta1[pp]*agrad_1_gt4[pp] + beta2[pp]*agrad_2_gt4[pp];
//--
gt_rhs22[pp] = -At5[pp]*DENDRO_0 + DENDRO_10*DENDRO_15 + DENDRO_16*DENDRO_21 - DENDRO_18*gt5[pp] - DENDRO_19*gt5[pp] + DENDRO_23*gt5[pp] + beta0[pp]*agrad_0_gt5[pp] + beta1[pp]*agrad_1_gt5[pp] + beta2[pp]*agrad_2_gt5[pp];
//--
chi_rhs[pp] = -DENDRO_24*DENDRO_25 + DENDRO_24*K[pp]*alpha[pp] + beta0[pp]*agrad_0_chi[pp] + beta1[pp]*agrad_1_chi[pp] + beta2[pp]*agrad_2_chi[pp];
//--
At_rhs00[pp] = -At0[pp]*DENDRO_14 - At0[pp]*DENDRO_18 + At0[pp]*DENDRO_3 + DENDRO_26*DENDRO_7 + DENDRO_27*DENDRO_9 + DENDRO_802*(-DENDRO_113*(DENDRO_100 - DENDRO_102 - DENDRO_111 + DENDRO_98) + DENDRO_351*(-DENDRO_121*(DENDRO_115 + DENDRO_116*DENDRO_117 + DENDRO_118*DENDRO_80 + DENDRO_119*DENDRO_120) - DENDRO_146*(DENDRO_135 - DENDRO_143) - DENDRO_146*(DENDRO_156 + DENDRO_163) - DENDRO_146*(DENDRO_166 + DENDRO_167*DENDRO_169) - DENDRO_146*(DENDRO_175 + DENDRO_179*DENDRO_70) - DENDRO_146*(DENDRO_173*DENDRO_75 + DENDRO_181) + DENDRO_164*DENDRO_312 + DENDRO_184*(DENDRO_193 + 1.0*DENDRO_194) - DENDRO_184*(-DENDRO_169*DENDRO_76 + DENDRO_182*DENDRO_93) + DENDRO_187*DENDRO_314*DENDRO_320 - DENDRO_190*(DENDRO_186 + 1.0*DENDRO_188) - DENDRO_190*(DENDRO_164*DENDRO_191 + 1.0*DENDRO_167*DENDRO_93) + DENDRO_197*(DENDRO_192 + DENDRO_194) + DENDRO_197*(DENDRO_195 + DENDRO_196) - DENDRO_200*(DENDRO_185 + DENDRO_188) - DENDRO_200*(DENDRO_198 + DENDRO_199) - DENDRO_233*DENDRO_234 - DENDRO_253*DENDRO_254 - DENDRO_271*DENDRO_272 - DENDRO_299*DENDRO_300 + DENDRO_303*(DENDRO_150 + DENDRO_165) + DENDRO_307*(DENDRO_126 + DENDRO_304) + DENDRO_313*DENDRO_315*DENDRO_93 + DENDRO_316*DENDRO_317 + DENDRO_318*DENDRO_319 + DENDRO_350) - 12*DENDRO_58 + DENDRO_801*gt0[pp] + DENDRO_81*DENDRO_84 - DENDRO_96*(-DENDRO_89 + DENDRO_93)) - alpha[pp]*(-At0[pp]*K[pp] + DENDRO_42*(DENDRO_30 + DENDRO_32 - DENDRO_36) + DENDRO_52*(DENDRO_43 + DENDRO_48 - DENDRO_50) + DENDRO_57*(-At0[pp]*DENDRO_46 + DENDRO_53 - DENDRO_56)) + beta0[pp]*agrad_0_At0[pp] + beta1[pp]*agrad_1_At0[pp] + beta2[pp]*agrad_2_At0[pp];
//--
At_rhs01[pp] = At0[pp]*DENDRO_11 - At1[pp]*DENDRO_14 + At2[pp]*DENDRO_12 + At3[pp]*DENDRO_7 + At4[pp]*DENDRO_9 + DENDRO_2*DENDRO_803 + DENDRO_4*DENDRO_803 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_141*DENDRO_592 + DENDRO_155*DENDRO_590 + DENDRO_242*DENDRO_595 + DENDRO_762) - DENDRO_124*DENDRO_834 - DENDRO_129*DENDRO_836 + DENDRO_146*(-DENDRO_332*DENDRO_822 + DENDRO_752) - DENDRO_146*(-DENDRO_687 + DENDRO_747 + DENDRO_748) + DENDRO_146*(-DENDRO_173*DENDRO_688 + DENDRO_179*DENDRO_472 + DENDRO_725) + DENDRO_184*(DENDRO_823 + DENDRO_826) + DENDRO_184*(DENDRO_827 + DENDRO_828) + DENDRO_184*(DENDRO_187*DENDRO_470 + DENDRO_829) + DENDRO_184*(DENDRO_163 + DENDRO_742 + DENDRO_743) + DENDRO_190*(-DENDRO_70*DENDRO_822 + DENDRO_760) + DENDRO_190*(DENDRO_127*DENDRO_187 - DENDRO_173*DENDRO_327 + DENDRO_755) - DENDRO_190*(DENDRO_164*DENDRO_561 + DENDRO_164*DENDRO_576 + DENDRO_226*DENDRO_93) + DENDRO_190*(-DENDRO_164*DENDRO_688 - DENDRO_227*DENDRO_93 + DENDRO_753) - DENDRO_200*(DENDRO_122*DENDRO_187 + DENDRO_318) + DENDRO_301*(DENDRO_815 + DENDRO_819) + DENDRO_302*(-DENDRO_405 + DENDRO_566 + DENDRO_673) + DENDRO_306*(-DENDRO_733 + DENDRO_734) - DENDRO_306*(-DENDRO_164*DENDRO_437 + DENDRO_736 + DENDRO_821) - DENDRO_306*(-DENDRO_173*DENDRO_304 + DENDRO_739 + DENDRO_820) + DENDRO_315*(DENDRO_186 + DENDRO_187*DENDRO_70 + DENDRO_187*DENDRO_71) + DENDRO_315*(0.25*DENDRO_198 + 0.5*DENDRO_199 + DENDRO_470*DENDRO_93) - DENDRO_67*DENDRO_835 + DENDRO_767 + DENDRO_769 + DENDRO_770 + DENDRO_771 + DENDRO_772 + DENDRO_773 + DENDRO_774 + DENDRO_775 + DENDRO_776 + DENDRO_777 + DENDRO_778 + DENDRO_779 + DENDRO_780 - DENDRO_781*DENDRO_837 + DENDRO_786 - DENDRO_831*(DENDRO_122*DENDRO_179 + DENDRO_833) - DENDRO_831*(DENDRO_164*DENDRO_204 + DENDRO_168*DENDRO_202 + DENDRO_830)) - DENDRO_41*DENDRO_812*(DENDRO_164 - DENDRO_797) + DENDRO_781*DENDRO_800 - 12*DENDRO_787 - DENDRO_813*(-DENDRO_791 + DENDRO_792 + DENDRO_793 - DENDRO_795) + DENDRO_814*(DENDRO_64*(DENDRO_621 + DENDRO_63*DENDRO_781) + DENDRO_790)) - alpha[pp]*(-At1[pp]*K[pp] + DENDRO_42*DENDRO_806 + DENDRO_52*DENDRO_809 + DENDRO_57*DENDRO_811) + beta0[pp]*agrad_0_At1[pp] + beta1[pp]*agrad_1_At1[pp] + beta2[pp]*agrad_2_At1[pp];
//--
At_rhs02[pp] = At0[pp]*DENDRO_15 + At1[pp]*DENDRO_16 - At2[pp]*DENDRO_18 + At4[pp]*DENDRO_7 + At5[pp]*DENDRO_9 + DENDRO_2*DENDRO_838 + DENDRO_6*DENDRO_838 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_134*DENDRO_592 + DENDRO_160*DENDRO_590 + DENDRO_235*DENDRO_595 + DENDRO_587) - DENDRO_130*DENDRO_834 - DENDRO_146*(DENDRO_554 - DENDRO_555 + DENDRO_641) + DENDRO_146*(-DENDRO_179*DENDRO_561 - DENDRO_70*DENDRO_846 + DENDRO_854) - DENDRO_146*(DENDRO_332*DENDRO_849 + DENDRO_847 + DENDRO_848) + DENDRO_146*(DENDRO_389*DENDRO_824 - DENDRO_848 + DENDRO_853) - DENDRO_148*DENDRO_836 - DENDRO_184*(DENDRO_151*DENDRO_169 - DENDRO_850 - DENDRO_851) - DENDRO_184*(DENDRO_151*DENDRO_187 - DENDRO_179*DENDRO_191 - DENDRO_78*DENDRO_846) + DENDRO_184*(DENDRO_75*DENDRO_849 + DENDRO_850 + DENDRO_851) - DENDRO_190*(DENDRO_166 + DENDRO_826) - DENDRO_190*(DENDRO_828 + DENDRO_852) - DENDRO_190*(DENDRO_187*DENDRO_408 + DENDRO_829) - DENDRO_190*(DENDRO_166 + DENDRO_168*DENDRO_576 + DENDRO_825) + DENDRO_197*(DENDRO_132*DENDRO_187 + DENDRO_316) + DENDRO_302*(DENDRO_534 - DENDRO_638) - DENDRO_302*(-DENDRO_165*DENDRO_179 - DENDRO_75*DENDRO_846 + DENDRO_845) + DENDRO_305*(DENDRO_833 + DENDRO_844) + DENDRO_306*(DENDRO_164*DENDRO_214 + 0.25*DENDRO_830) + DENDRO_315*(0.25*DENDRO_195 + 1.0*DENDRO_196) + DENDRO_315*(DENDRO_187*DENDRO_75 + DENDRO_187*DENDRO_76 + DENDRO_193) + DENDRO_597 + DENDRO_599 + DENDRO_601 + DENDRO_602 + DENDRO_604 + DENDRO_606 + DENDRO_608 + DENDRO_609 + DENDRO_610 + DENDRO_611 + DENDRO_612 + DENDRO_613 + DENDRO_614 - DENDRO_615*DENDRO_837 + DENDRO_655 - DENDRO_72*DENDRO_835 - DENDRO_831*(DENDRO_819 + DENDRO_855)) + DENDRO_615*DENDRO_800 - 12*DENDRO_617 + DENDRO_634*DENDRO_814 - DENDRO_812*(DENDRO_618 + DENDRO_619 - DENDRO_620 - DENDRO_623) - DENDRO_813*(DENDRO_625 + DENDRO_626 - DENDRO_627 - DENDRO_630)) - alpha[pp]*(-At2[pp]*K[pp] + DENDRO_42*DENDRO_840 + DENDRO_52*DENDRO_842 + DENDRO_57*DENDRO_843) + beta0[pp]*agrad_0_At2[pp] + beta1[pp]*agrad_1_At2[pp] + beta2[pp]*agrad_2_At2[pp];
//--
At_rhs11[pp] = -At3[pp]*DENDRO_14 - At3[pp]*DENDRO_19 + At3[pp]*DENDRO_20 + DENDRO_11*DENDRO_26 + DENDRO_12*DENDRO_856 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_117*DENDRO_268 + DENDRO_118*DENDRO_230 + DENDRO_120*DENDRO_250 + DENDRO_456) - DENDRO_122*DENDRO_835 - DENDRO_146*(DENDRO_473 - 1.0*DENDRO_475) + DENDRO_146*(DENDRO_480 - 1.0*DENDRO_860) - DENDRO_146*(DENDRO_476*DENDRO_824 + DENDRO_478) + DENDRO_164*DENDRO_498 + DENDRO_173*DENDRO_497 + DENDRO_184*(DENDRO_466 - DENDRO_468) + DENDRO_184*(DENDRO_173*DENDRO_470 + DENDRO_180*DENDRO_817) + DENDRO_184*(DENDRO_470*DENDRO_824 + DENDRO_863) + DENDRO_184*(DENDRO_71*DENDRO_861 + DENDRO_862) + DENDRO_190*(DENDRO_483 + DENDRO_821) + DENDRO_190*(DENDRO_485 + DENDRO_820) + DENDRO_190*(DENDRO_487 - 1.0*DENDRO_864) + DENDRO_200*(DENDRO_486 - DENDRO_864) + DENDRO_200*(-DENDRO_122*DENDRO_173 + DENDRO_492) + DENDRO_200*(-DENDRO_164*DENDRO_202 + DENDRO_491) - DENDRO_202*DENDRO_836 + 6.0*DENDRO_221*DENDRO_305*DENDRO_822 - DENDRO_221*DENDRO_834 + DENDRO_410*(DENDRO_479 - DENDRO_860) + DENDRO_410*(-DENDRO_122*DENDRO_817 + DENDRO_489) + DENDRO_410*(-DENDRO_202*DENDRO_824 + DENDRO_490) - DENDRO_451*DENDRO_837 + DENDRO_493*(DENDRO_217 + DENDRO_404) + DENDRO_495*DENDRO_817 + DENDRO_496*(DENDRO_126 + DENDRO_457) + DENDRO_506) - 12*DENDRO_445 + DENDRO_801*gt3[pp] + DENDRO_83*(DENDRO_455 + DENDRO_64*(DENDRO_449 + DENDRO_452*DENDRO_63)) + DENDRO_859*(DENDRO_250 - DENDRO_446*(-DENDRO_106 + DENDRO_357 + DENDRO_358)) + DENDRO_96*(DENDRO_268 - DENDRO_446*(-DENDRO_85 + DENDRO_86 + DENDRO_87))) - alpha[pp]*(-At3[pp]*K[pp] + DENDRO_42*DENDRO_809 + DENDRO_806*DENDRO_857 + DENDRO_811*DENDRO_858) + beta0[pp]*agrad_0_At3[pp] + beta1[pp]*agrad_1_At3[pp] + beta2[pp]*agrad_2_At3[pp];
//--
At_rhs12[pp] = At1[pp]*DENDRO_15 + At2[pp]*DENDRO_11 + At3[pp]*DENDRO_16 - At4[pp]*DENDRO_19 + At5[pp]*DENDRO_12 + DENDRO_4*DENDRO_865 + DENDRO_6*DENDRO_865 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_208*DENDRO_592 + DENDRO_240*DENDRO_595 + DENDRO_259*DENDRO_590 + DENDRO_696) - DENDRO_128*DENDRO_835 + DENDRO_146*(DENDRO_218*DENDRO_822 + DENDRO_723) + DENDRO_146*(-DENDRO_226*DENDRO_849 + DENDRO_867) + DENDRO_146*(DENDRO_385*DENDRO_824 + DENDRO_867) - DENDRO_146*(DENDRO_123*DENDRO_846 + DENDRO_561*DENDRO_817 + DENDRO_722) + DENDRO_146*(-DENDRO_576*DENDRO_817 + DENDRO_677 - DENDRO_688*DENDRO_817) + DENDRO_183*(DENDRO_815 + DENDRO_818 + DENDRO_855) - DENDRO_184*(DENDRO_169*DENDRO_218 + DENDRO_853 - DENDRO_866) - DENDRO_184*(-DENDRO_191*DENDRO_817 - DENDRO_71*DENDRO_846 + DENDRO_854) + DENDRO_184*(DENDRO_470*DENDRO_849 + DENDRO_847 + DENDRO_866) + DENDRO_190*(-DENDRO_168*DENDRO_437 + DENDRO_868) + DENDRO_190*(-DENDRO_408*DENDRO_822 + DENDRO_728) + DENDRO_190*(-DENDRO_688*DENDRO_824 + DENDRO_868) + DENDRO_190*(-DENDRO_179*DENDRO_304 + DENDRO_725 - 0.25*DENDRO_844) + DENDRO_190*(-DENDRO_327*DENDRO_817 + DENDRO_685 - DENDRO_862) + DENDRO_190*(-DENDRO_470*DENDRO_822 + DENDRO_539 + DENDRO_691) - DENDRO_215*DENDRO_836 - DENDRO_224*DENDRO_834 + DENDRO_302*(DENDRO_659 - DENDRO_717) - DENDRO_302*(-DENDRO_165*DENDRO_817 + DENDRO_389*DENDRO_817 - DENDRO_470*DENDRO_846) - DENDRO_306*(-DENDRO_437*DENDRO_824 + DENDRO_718) - DENDRO_306*(-DENDRO_226*DENDRO_822 - DENDRO_227*DENDRO_822 + DENDRO_480) - DENDRO_306*(-DENDRO_304*DENDRO_817 + DENDRO_664 + DENDRO_719) + DENDRO_315*(DENDRO_147*DENDRO_164 + DENDRO_823) + DENDRO_315*(DENDRO_175 + DENDRO_827 + DENDRO_852) + DENDRO_410*(-DENDRO_204*DENDRO_822 + DENDRO_504) - DENDRO_701*DENDRO_837 + DENDRO_702 + DENDRO_731) - DENDRO_41*DENDRO_813*(-DENDRO_715 + DENDRO_817) + DENDRO_701*DENDRO_800 - 12*DENDRO_704 - DENDRO_812*(-DENDRO_705 + DENDRO_706 + DENDRO_707 - DENDRO_710) + DENDRO_814*(DENDRO_64*(DENDRO_628 + DENDRO_63*DENDRO_701) + DENDRO_713)) - alpha[pp]*(-At4[pp]*K[pp] + DENDRO_42*DENDRO_842 + DENDRO_840*DENDRO_857 + DENDRO_843*DENDRO_858) + beta0[pp]*agrad_0_At4[pp] + beta1[pp]*agrad_1_At4[pp] + beta2[pp]*agrad_2_At4[pp];
//--
At_rhs22[pp] = -At5[pp]*DENDRO_18 - At5[pp]*DENDRO_19 + At5[pp]*DENDRO_23 + DENDRO_15*DENDRO_27 + DENDRO_16*DENDRO_856 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_117*DENDRO_265 + DENDRO_118*DENDRO_219 + DENDRO_120*DENDRO_247 + DENDRO_374) + DENDRO_132*DENDRO_168*DENDRO_427 - DENDRO_132*DENDRO_835 - DENDRO_146*(DENDRO_395 - 1.0*DENDRO_397) - DENDRO_146*(0.25*DENDRO_869 + 1.0*DENDRO_870) + DENDRO_179*DENDRO_424 + DENDRO_184*(DENDRO_386 - 1.0*DENDRO_388) + DENDRO_184*(0.25*DENDRO_871 + 1.0*DENDRO_872) - DENDRO_184*(-DENDRO_391*DENDRO_846 + DENDRO_845) - DENDRO_190*(DENDRO_147*DENDRO_824 + DENDRO_848) - DENDRO_190*(DENDRO_168*DENDRO_214 + DENDRO_866) - DENDRO_190*(DENDRO_174*DENDRO_817 + DENDRO_179*DENDRO_408) - DENDRO_190*(DENDRO_76*DENDRO_861 + 0.25*DENDRO_816) + DENDRO_197*(DENDRO_871 + DENDRO_872) + DENDRO_197*(DENDRO_132*DENDRO_179 + DENDRO_74*DENDRO_846) + DENDRO_204*DENDRO_319*DENDRO_824 - DENDRO_204*DENDRO_834 + 6.0*DENDRO_211*DENDRO_301*DENDRO_849 - DENDRO_211*DENDRO_836 - DENDRO_367*DENDRO_837 - DENDRO_410*(DENDRO_869 + DENDRO_870) - DENDRO_410*(DENDRO_132*DENDRO_817 + DENDRO_167*DENDRO_846) + DENDRO_420*(DENDRO_150 + DENDRO_379) + DENDRO_421*(DENDRO_217 + DENDRO_376) + DENDRO_423*DENDRO_817 + DENDRO_444) - 12*DENDRO_352 + DENDRO_354*DENDRO_84 + DENDRO_801*gt5[pp] - DENDRO_859*(-DENDRO_360 + DENDRO_846) - DENDRO_95*(DENDRO_363 - DENDRO_364 + DENDRO_365 - DENDRO_372)) - alpha[pp]*(At5[pp]*DENDRO_51*DENDRO_843 - At5[pp]*K[pp] + DENDRO_57*DENDRO_842 + DENDRO_840*DENDRO_858) + beta0[pp]*agrad_0_At5[pp] + beta1[pp]*agrad_1_At5[pp] + beta2[pp]*agrad_2_At5[pp];
//--
K_rhs[pp] = -DENDRO_334*chi[pp]*(-DENDRO_352 + DENDRO_873*(DENDRO_353*DENDRO_450 + DENDRO_383) + DENDRO_875*(DENDRO_384 + DENDRO_874*gt5[pp]) + DENDRO_94*(DENDRO_375*DENDRO_41 + DENDRO_378*DENDRO_41 + DENDRO_381*DENDRO_41 - DENDRO_64*(DENDRO_366 - DENDRO_371))) - DENDRO_336*chi[pp]*(-DENDRO_445 + DENDRO_82*(DENDRO_41*DENDRO_460 + DENDRO_41*DENDRO_461 + DENDRO_454 - DENDRO_64*(DENDRO_448 - DENDRO_453)) + DENDRO_875*(DENDRO_462 + DENDRO_874*gt3[pp]) + DENDRO_877*(DENDRO_459 + DENDRO_876*gt3[pp])) - DENDRO_338*chi[pp]*(DENDRO_112*(DENDRO_102 + DENDRO_41*DENDRO_509 + DENDRO_41*DENDRO_510 - DENDRO_64*(DENDRO_103 - DENDRO_110)) - DENDRO_58 + DENDRO_873*(DENDRO_450*DENDRO_66 + DENDRO_508) + DENDRO_877*(DENDRO_507 + DENDRO_876*gt0[pp])) + DENDRO_880*chi[pp]*(-DENDRO_704 + 0.5*DENDRO_875*(DENDRO_108*DENDRO_714 + DENDRO_700) + DENDRO_878*(DENDRO_41*DENDRO_697 + DENDRO_41*DENDRO_698 - DENDRO_64*(DENDRO_59 - DENDRO_709) + DENDRO_705) + DENDRO_879*(DENDRO_41*DENDRO_699 - DENDRO_64*(DENDRO_60 - DENDRO_711) + DENDRO_712)) - DENDRO_882*chi[pp]*(-DENDRO_617 + 0.5*DENDRO_873*(DENDRO_450*DENDRO_633 + DENDRO_591) + DENDRO_878*(DENDRO_41*DENDRO_588 + DENDRO_41*DENDRO_589 + DENDRO_620 - DENDRO_64*(DENDRO_61 - DENDRO_622)) + DENDRO_881*(DENDRO_41*DENDRO_593 + DENDRO_41*DENDRO_594 + DENDRO_627 - DENDRO_64*(DENDRO_60 - DENDRO_629))) + DENDRO_883*chi[pp]*(-DENDRO_787 + 0.5*DENDRO_877*(DENDRO_369*DENDRO_796 + DENDRO_763) + DENDRO_879*(DENDRO_41*DENDRO_764 - DENDRO_64*(DENDRO_61 - DENDRO_788) + DENDRO_789) + DENDRO_881*(DENDRO_41*DENDRO_765 + DENDRO_41*DENDRO_766 - DENDRO_64*(DENDRO_59 - DENDRO_794) + DENDRO_791)) + DENDRO_916*(At0[pp]*DENDRO_887*DENDRO_888 + At1[pp]*DENDRO_895*DENDRO_905 + At2[pp]*DENDRO_894*DENDRO_895 + At3[pp]*DENDRO_888*DENDRO_891 + At4[pp]*DENDRO_895*DENDRO_915 + At5[pp]*DENDRO_888*DENDRO_893 + pow(K[pp], 2)) + beta0[pp]*agrad_0_K[pp] + beta1[pp]*agrad_1_K[pp] + beta2[pp]*agrad_2_K[pp];
//--
Gt_rhs0[pp] = DENDRO_955;
//--
Gt_rhs1[pp] = -DENDRO_141*DENDRO_950*DENDRO_956 + DENDRO_16*DENDRO_271 - DENDRO_208*DENDRO_950*DENDRO_957 + DENDRO_233*DENDRO_4 - DENDRO_233*DENDRO_942 + DENDRO_253*DENDRO_7 - DENDRO_822*DENDRO_947 + DENDRO_918*DENDRO_956 + DENDRO_920*DENDRO_957 + DENDRO_927*(DENDRO_960 - DENDRO_961) - DENDRO_927*(-DENDRO_926*DENDRO_956 + DENDRO_959) - DENDRO_927*(-DENDRO_929*DENDRO_957 + DENDRO_958) + DENDRO_968;
//--
Gt_rhs2[pp] = DENDRO_970;
//--
B_rhs0[pp] = -B0[pp]*DENDRO_972 - DENDRO_954*lambda[3] + DENDRO_955 + lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] + beta2[pp]*agrad_2_B0[pp]);
//--
B_rhs1[pp] = -B1[pp]*DENDRO_972 + DENDRO_141*DENDRO_952 - DENDRO_16*DENDRO_413 + DENDRO_208*DENDRO_953 + DENDRO_230*DENDRO_947 - DENDRO_4*DENDRO_411 + DENDRO_411*DENDRO_942 - DENDRO_415*DENDRO_7 - DENDRO_905*DENDRO_918 - DENDRO_915*DENDRO_920 - DENDRO_927*(-DENDRO_960 + DENDRO_961) - DENDRO_927*(DENDRO_905*DENDRO_926 + DENDRO_959) - DENDRO_927*(DENDRO_915*DENDRO_929 + DENDRO_958) - DENDRO_967*lambda[3] + DENDRO_968 + lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] + beta2[pp]*agrad_2_B1[pp]);
//--
B_rhs2[pp] = -B2[pp]*DENDRO_972 - DENDRO_969*lambda[3] + DENDRO_970 + lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] + beta2[pp]*agrad_2_B2[pp]);
// Dendro: reduced ops: 4158
// Dendro: }}} 
//[[[end]]]
#else
/*[[[cog

import dendro
import bssn

# modify to use the eta const value. 
bssn.eta_func=bssn.eta
bssn.B_rhs = [bssn.Gt_rhs[i] - bssn.eta_func * bssn.B[i] +
         bssn.l3 * dendro.vec_j_ad_j(bssn.b, bssn.B[i]) -
         bssn.l4 * dendro.vec_j_ad_j(bssn.b, bssn.Gt[i])
         for i in dendro.e_i]

outs = [bssn.a_rhs, bssn.b_rhs, bssn.gt_rhs, bssn.chi_rhs, bssn.At_rhs, bssn.K_rhs, bssn.Gt_rhs, bssn.B_rhs]
vnames = ['a_rhs', 'b_rhs', 'gt_rhs', 'chi_rhs', 'At_rhs', 'K_rhs', 'Gt_rhs', 'B_rhs']
#outs=outs[0:3]
#vnames=vnames[0:3]
dendro.generate_cpu(outs, vnames, '[pp]')


]]]*/
// Dendro: {{{ 
// Dendro: original ops: 667753 
// Dendro: printing temp variables
double DENDRO_0 = 2*alpha[pp];
double DENDRO_1 = (3.0/4.0)*alpha[pp]*lambda_f[1] + (3.0/4.0)*lambda_f[0];
double DENDRO_2 = grad_0_beta0[pp];
double DENDRO_3 = (4.0/3.0)*DENDRO_2;
double DENDRO_4 = grad_1_beta1[pp];
double DENDRO_5 = (2.0/3.0)*gt0[pp];
double DENDRO_6 = grad_2_beta2[pp];
double DENDRO_7 = grad_0_beta1[pp];
double DENDRO_8 = 2*gt1[pp];
double DENDRO_9 = grad_0_beta2[pp];
double DENDRO_10 = 2*gt2[pp];
double DENDRO_11 = grad_1_beta0[pp];
double DENDRO_12 = grad_1_beta2[pp];
double DENDRO_13 = (1.0/3.0)*gt1[pp];
double DENDRO_14 = (2.0/3.0)*DENDRO_6;
double DENDRO_15 = grad_2_beta0[pp];
double DENDRO_16 = grad_2_beta1[pp];
double DENDRO_17 = (1.0/3.0)*gt2[pp];
double DENDRO_18 = (2.0/3.0)*DENDRO_4;
double DENDRO_19 = (2.0/3.0)*DENDRO_2;
double DENDRO_20 = (4.0/3.0)*DENDRO_4;
double DENDRO_21 = 2*gt4[pp];
double DENDRO_22 = (1.0/3.0)*gt4[pp];
double DENDRO_23 = (4.0/3.0)*DENDRO_6;
double DENDRO_24 = (2.0/3.0)*chi[pp];
double DENDRO_25 = DENDRO_2 + DENDRO_4 + DENDRO_6;
double DENDRO_26 = 2*At1[pp];
double DENDRO_27 = 2*At2[pp];
double DENDRO_28 = gt2[pp]*gt4[pp];
double DENDRO_29 = -DENDRO_28 + gt1[pp]*gt5[pp];
double DENDRO_30 = At0[pp]*DENDRO_29;
double DENDRO_31 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
double DENDRO_32 = At2[pp]*DENDRO_31;
double DENDRO_33 = gt0[pp]*gt5[pp];
double DENDRO_34 = pow(gt2[pp], 2);
double DENDRO_35 = DENDRO_33 - DENDRO_34;
double DENDRO_36 = At1[pp]*DENDRO_35;
double DENDRO_37 = pow(gt4[pp], 2);
double DENDRO_38 = pow(gt1[pp], 2);
double DENDRO_39 = gt3[pp]*gt5[pp];
double DENDRO_40 = -DENDRO_28*DENDRO_8 + DENDRO_34*gt3[pp] + DENDRO_37*gt0[pp] + DENDRO_38*gt5[pp] - DENDRO_39*gt0[pp];
double DENDRO_41 = 1.0/DENDRO_40;
double DENDRO_42 = DENDRO_26*DENDRO_41;
double DENDRO_43 = At1[pp]*DENDRO_29;
double DENDRO_44 = gt1[pp]*gt4[pp];
double DENDRO_45 = gt2[pp]*gt3[pp];
double DENDRO_46 = DENDRO_44 - DENDRO_45;
double DENDRO_47 = At2[pp]*DENDRO_46;
double DENDRO_48 = -DENDRO_47;
double DENDRO_49 = -DENDRO_37 + DENDRO_39;
double DENDRO_50 = At0[pp]*DENDRO_49;
double DENDRO_51 = 2*DENDRO_41;
double DENDRO_52 = At0[pp]*DENDRO_51;
double DENDRO_53 = At1[pp]*DENDRO_31;
double DENDRO_54 = gt0[pp]*gt3[pp];
double DENDRO_55 = -DENDRO_38 + DENDRO_54;
double DENDRO_56 = At2[pp]*DENDRO_55;
double DENDRO_57 = DENDRO_27*DENDRO_41;
double DENDRO_58 = grad2_0_0_alpha[pp];
double DENDRO_59 = grad_1_chi[pp];
double DENDRO_60 = grad_2_chi[pp];
double DENDRO_61 = grad_0_chi[pp];
double DENDRO_62 = DENDRO_29*DENDRO_61 + DENDRO_31*DENDRO_60;
double DENDRO_63 = -DENDRO_35*DENDRO_59 + DENDRO_62;
double DENDRO_64 = 1.0/chi[pp];
double DENDRO_65 = 0.5*DENDRO_64;
double DENDRO_66 = DENDRO_65*gt0[pp];
double DENDRO_67 = grad_0_gt1[pp];
double DENDRO_68 = 1.0*DENDRO_67;
double DENDRO_69 = grad_1_gt0[pp];
double DENDRO_70 = 0.5*DENDRO_69;
double DENDRO_71 = DENDRO_68 - DENDRO_70;
double DENDRO_72 = grad_0_gt2[pp];
double DENDRO_73 = 1.0*DENDRO_72;
double DENDRO_74 = grad_2_gt0[pp];
double DENDRO_75 = 0.5*DENDRO_74;
double DENDRO_76 = DENDRO_73 - DENDRO_75;
double DENDRO_77 = grad_0_gt0[pp];
double DENDRO_78 = 0.5*DENDRO_77;
double DENDRO_79 = DENDRO_29*DENDRO_78 + DENDRO_31*DENDRO_76;
double DENDRO_80 = -DENDRO_35*DENDRO_71 + DENDRO_79;
double DENDRO_81 = DENDRO_63*DENDRO_66 + DENDRO_80;
double DENDRO_82 = grad_1_alpha[pp];
double DENDRO_83 = 12*DENDRO_82;
double DENDRO_84 = DENDRO_41*DENDRO_83;
double DENDRO_85 = DENDRO_31*DENDRO_59;
double DENDRO_86 = DENDRO_46*DENDRO_61;
double DENDRO_87 = DENDRO_55*DENDRO_60;
double DENDRO_88 = DENDRO_85 - DENDRO_86 - DENDRO_87;
double DENDRO_89 = DENDRO_66*DENDRO_88;
double DENDRO_90 = DENDRO_46*DENDRO_78;
double DENDRO_91 = DENDRO_55*DENDRO_76;
double DENDRO_92 = DENDRO_31*DENDRO_71;
double DENDRO_93 = DENDRO_90 + DENDRO_91 - DENDRO_92;
double DENDRO_94 = grad_2_alpha[pp];
double DENDRO_95 = 12*DENDRO_94;
double DENDRO_96 = DENDRO_41*DENDRO_95;
double DENDRO_97 = DENDRO_49*DENDRO_78;
double DENDRO_98 = DENDRO_41*DENDRO_97;
double DENDRO_99 = DENDRO_46*DENDRO_76;
double DENDRO_100 = DENDRO_41*DENDRO_99;
double DENDRO_101 = DENDRO_29*DENDRO_71;
double DENDRO_102 = DENDRO_101*DENDRO_41;
double DENDRO_103 = 1.0*DENDRO_61;
double DENDRO_104 = DENDRO_41*gt0[pp];
double DENDRO_105 = -DENDRO_44 + DENDRO_45;
double DENDRO_106 = DENDRO_29*DENDRO_59;
double DENDRO_107 = DENDRO_37 - DENDRO_39;
double DENDRO_108 = DENDRO_105*DENDRO_60 + DENDRO_106 + DENDRO_107*DENDRO_61;
double DENDRO_109 = 0.5*DENDRO_108;
double DENDRO_110 = DENDRO_104*DENDRO_109;
double DENDRO_111 = DENDRO_64*(-DENDRO_103 + DENDRO_110);
double DENDRO_112 = grad_0_alpha[pp];
double DENDRO_113 = 12*DENDRO_112;
double DENDRO_114 = grad2_0_0_chi[pp];
double DENDRO_115 = -DENDRO_114;
double DENDRO_116 = -DENDRO_90 - DENDRO_91 + DENDRO_92;
double DENDRO_117 = DENDRO_41*DENDRO_60;
double DENDRO_118 = DENDRO_41*DENDRO_59;
double DENDRO_119 = DENDRO_101 - DENDRO_97 - DENDRO_99;
double DENDRO_120 = DENDRO_41*DENDRO_61;
double DENDRO_121 = 2*DENDRO_64;
double DENDRO_122 = grad_0_gt3[pp];
double DENDRO_123 = 0.5*DENDRO_122;
double DENDRO_124 = grad_1_gt1[pp];
double DENDRO_125 = 1.0*DENDRO_124;
double DENDRO_126 = -DENDRO_125;
double DENDRO_127 = DENDRO_123 + DENDRO_126;
double DENDRO_128 = grad_0_gt4[pp];
double DENDRO_129 = grad_2_gt1[pp];
double DENDRO_130 = grad_1_gt2[pp];
double DENDRO_131 = DENDRO_128 + DENDRO_129 - DENDRO_130;
double DENDRO_132 = grad_0_gt5[pp];
double DENDRO_133 = DENDRO_132*DENDRO_31 + DENDRO_29*DENDRO_74;
double DENDRO_134 = -DENDRO_131*DENDRO_35 + DENDRO_133;
double DENDRO_135 = DENDRO_127*DENDRO_134;
double DENDRO_136 = DENDRO_122*DENDRO_35;
double DENDRO_137 = DENDRO_128 - DENDRO_129 + DENDRO_130;
double DENDRO_138 = DENDRO_137*DENDRO_31;
double DENDRO_139 = DENDRO_29*DENDRO_69;
double DENDRO_140 = DENDRO_138 + DENDRO_139;
double DENDRO_141 = -DENDRO_136 + DENDRO_140;
double DENDRO_142 = DENDRO_131*DENDRO_141;
double DENDRO_143 = 0.25*DENDRO_142;
double DENDRO_144 = pow(DENDRO_40, -2);
double DENDRO_145 = DENDRO_144*DENDRO_31;
double DENDRO_146 = 4*DENDRO_145;
double DENDRO_147 = 0.5*DENDRO_132;
double DENDRO_148 = grad_2_gt2[pp];
double DENDRO_149 = 1.0*DENDRO_148;
double DENDRO_150 = -DENDRO_149;
double DENDRO_151 = DENDRO_147 + DENDRO_150;
double DENDRO_152 = DENDRO_122*DENDRO_31;
double DENDRO_153 = DENDRO_46*DENDRO_69;
double DENDRO_154 = DENDRO_137*DENDRO_55;
double DENDRO_155 = DENDRO_152 - DENDRO_153 - DENDRO_154;
double DENDRO_156 = DENDRO_151*DENDRO_155;
double DENDRO_157 = DENDRO_46*DENDRO_74;
double DENDRO_158 = DENDRO_132*DENDRO_55;
double DENDRO_159 = DENDRO_131*DENDRO_31;
double DENDRO_160 = -DENDRO_157 - DENDRO_158 + DENDRO_159;
double DENDRO_161 = DENDRO_137*DENDRO_160;
double DENDRO_162 = 0.25*DENDRO_161;
double DENDRO_163 = -DENDRO_162;
double DENDRO_164 = -DENDRO_152 + DENDRO_153 + DENDRO_154;
double DENDRO_165 = 0.25*DENDRO_132;
double DENDRO_166 = DENDRO_164*DENDRO_165;
double DENDRO_167 = -DENDRO_128 + DENDRO_129 + DENDRO_130;
double DENDRO_168 = DENDRO_157 + DENDRO_158 - DENDRO_159;
double DENDRO_169 = 0.5*DENDRO_168;
double DENDRO_170 = DENDRO_122*DENDRO_29;
double DENDRO_171 = DENDRO_49*DENDRO_69;
double DENDRO_172 = DENDRO_137*DENDRO_46;
double DENDRO_173 = -DENDRO_170 + DENDRO_171 + DENDRO_172;
double DENDRO_174 = 0.25*DENDRO_74;
double DENDRO_175 = DENDRO_173*DENDRO_174;
double DENDRO_176 = DENDRO_132*DENDRO_46;
double DENDRO_177 = DENDRO_49*DENDRO_74;
double DENDRO_178 = DENDRO_131*DENDRO_29;
double DENDRO_179 = DENDRO_176 + DENDRO_177 - DENDRO_178;
double DENDRO_180 = 0.25*DENDRO_69;
double DENDRO_181 = DENDRO_179*DENDRO_180;
double DENDRO_182 = 2*DENDRO_151;
double DENDRO_183 = DENDRO_144*DENDRO_46;
double DENDRO_184 = 4*DENDRO_183;
double DENDRO_185 = DENDRO_173*DENDRO_77;
double DENDRO_186 = 0.25*DENDRO_185;
double DENDRO_187 = -DENDRO_101 + DENDRO_97 + DENDRO_99;
double DENDRO_188 = DENDRO_187*DENDRO_69;
double DENDRO_189 = DENDRO_144*DENDRO_29;
double DENDRO_190 = 4*DENDRO_189;
double DENDRO_191 = 0.5*DENDRO_76;
double DENDRO_192 = DENDRO_179*DENDRO_77;
double DENDRO_193 = 0.25*DENDRO_192;
double DENDRO_194 = DENDRO_187*DENDRO_74;
double DENDRO_195 = DENDRO_168*DENDRO_74;
double DENDRO_196 = DENDRO_132*DENDRO_93;
double DENDRO_197 = 2.0*DENDRO_183;
double DENDRO_198 = DENDRO_164*DENDRO_74;
double DENDRO_199 = DENDRO_137*DENDRO_93;
double DENDRO_200 = 2.0*DENDRO_189;
double DENDRO_201 = DENDRO_134*DENDRO_46;
double DENDRO_202 = grad_2_gt3[pp];
double DENDRO_203 = DENDRO_202*DENDRO_35;
double DENDRO_204 = grad_1_gt5[pp];
double DENDRO_205 = DENDRO_204*DENDRO_31;
double DENDRO_206 = DENDRO_167*DENDRO_29;
double DENDRO_207 = DENDRO_205 + DENDRO_206;
double DENDRO_208 = -DENDRO_203 + DENDRO_207;
double DENDRO_209 = DENDRO_208*DENDRO_31;
double DENDRO_210 = DENDRO_141*DENDRO_29;
double DENDRO_211 = grad_2_gt5[pp];
double DENDRO_212 = 0.5*DENDRO_211;
double DENDRO_213 = DENDRO_212*DENDRO_31;
double DENDRO_214 = 0.5*DENDRO_204;
double DENDRO_215 = grad_2_gt4[pp];
double DENDRO_216 = 1.0*DENDRO_215;
double DENDRO_217 = -DENDRO_216;
double DENDRO_218 = DENDRO_214 + DENDRO_217;
double DENDRO_219 = -DENDRO_151*DENDRO_29 + DENDRO_213 + DENDRO_218*DENDRO_35;
double DENDRO_220 = DENDRO_219*DENDRO_55;
double DENDRO_221 = grad_1_gt3[pp];
double DENDRO_222 = 0.5*DENDRO_221;
double DENDRO_223 = DENDRO_222*DENDRO_35;
double DENDRO_224 = grad_1_gt4[pp];
double DENDRO_225 = 1.0*DENDRO_224;
double DENDRO_226 = 0.5*DENDRO_202;
double DENDRO_227 = DENDRO_225 - DENDRO_226;
double DENDRO_228 = DENDRO_227*DENDRO_31;
double DENDRO_229 = DENDRO_127*DENDRO_29;
double DENDRO_230 = -DENDRO_223 + DENDRO_228 - DENDRO_229;
double DENDRO_231 = DENDRO_230*DENDRO_35;
double DENDRO_232 = DENDRO_49*DENDRO_80;
double DENDRO_233 = DENDRO_144*(DENDRO_201 - 1.0*DENDRO_209 - 1.0*DENDRO_210 + DENDRO_220 + DENDRO_231 + DENDRO_232);
double DENDRO_234 = 2.0*DENDRO_69;
double DENDRO_235 = -DENDRO_176 - DENDRO_177 + DENDRO_178;
double DENDRO_236 = DENDRO_235*DENDRO_46;
double DENDRO_237 = DENDRO_202*DENDRO_29;
double DENDRO_238 = DENDRO_204*DENDRO_46;
double DENDRO_239 = DENDRO_167*DENDRO_49;
double DENDRO_240 = DENDRO_237 - DENDRO_238 - DENDRO_239;
double DENDRO_241 = DENDRO_240*DENDRO_31;
double DENDRO_242 = DENDRO_170 - DENDRO_171 - DENDRO_172;
double DENDRO_243 = DENDRO_242*DENDRO_29;
double DENDRO_244 = DENDRO_212*DENDRO_46;
double DENDRO_245 = DENDRO_218*DENDRO_29;
double DENDRO_246 = DENDRO_151*DENDRO_49;
double DENDRO_247 = -DENDRO_244 - DENDRO_245 + DENDRO_246;
double DENDRO_248 = DENDRO_247*DENDRO_55;
double DENDRO_249 = DENDRO_222*DENDRO_29;
double DENDRO_250 = DENDRO_127*DENDRO_49 - DENDRO_227*DENDRO_46 + DENDRO_249;
double DENDRO_251 = DENDRO_250*DENDRO_35;
double DENDRO_252 = DENDRO_119*DENDRO_49;
double DENDRO_253 = DENDRO_144*(DENDRO_236 - 1.0*DENDRO_241 - 1.0*DENDRO_243 + DENDRO_248 + DENDRO_251 + DENDRO_252);
double DENDRO_254 = 2.0*DENDRO_77;
double DENDRO_255 = DENDRO_160*DENDRO_46;
double DENDRO_256 = DENDRO_202*DENDRO_31;
double DENDRO_257 = DENDRO_204*DENDRO_55;
double DENDRO_258 = DENDRO_167*DENDRO_46;
double DENDRO_259 = DENDRO_256 - DENDRO_257 - DENDRO_258;
double DENDRO_260 = DENDRO_259*DENDRO_31;
double DENDRO_261 = DENDRO_155*DENDRO_29;
double DENDRO_262 = DENDRO_212*DENDRO_55;
double DENDRO_263 = DENDRO_151*DENDRO_46;
double DENDRO_264 = DENDRO_218*DENDRO_31;
double DENDRO_265 = -DENDRO_262 + DENDRO_263 - DENDRO_264;
double DENDRO_266 = DENDRO_265*DENDRO_55;
double DENDRO_267 = DENDRO_222*DENDRO_31;
double DENDRO_268 = DENDRO_127*DENDRO_46 - DENDRO_227*DENDRO_55 + DENDRO_267;
double DENDRO_269 = DENDRO_268*DENDRO_35;
double DENDRO_270 = DENDRO_116*DENDRO_49;
double DENDRO_271 = DENDRO_144*(DENDRO_255 - 1.0*DENDRO_260 - 1.0*DENDRO_261 + DENDRO_266 + DENDRO_269 + DENDRO_270);
double DENDRO_272 = 2.0*DENDRO_74;
double DENDRO_273 = grad2_2_2_chi[pp];
double DENDRO_274 = pow(DENDRO_60, 2);
double DENDRO_275 = 3*DENDRO_64;
double DENDRO_276 = DENDRO_55*(2*DENDRO_273 - DENDRO_274*DENDRO_275);
double DENDRO_277 = grad2_1_1_chi[pp];
double DENDRO_278 = pow(DENDRO_59, 2);
double DENDRO_279 = DENDRO_35*(-DENDRO_275*DENDRO_278 + 2*DENDRO_277);
double DENDRO_280 = pow(DENDRO_61, 2);
double DENDRO_281 = DENDRO_49*(2*DENDRO_114 - DENDRO_275*DENDRO_280);
double DENDRO_282 = grad2_1_2_chi[pp];
double DENDRO_283 = DENDRO_275*DENDRO_60;
double DENDRO_284 = 2*DENDRO_31;
double DENDRO_285 = DENDRO_284*(2*DENDRO_282 - DENDRO_283*DENDRO_59);
double DENDRO_286 = grad2_0_2_chi[pp];
double DENDRO_287 = 2*DENDRO_46;
double DENDRO_288 = DENDRO_287*(-DENDRO_283*DENDRO_61 + 2*DENDRO_286);
double DENDRO_289 = grad2_0_1_chi[pp];
double DENDRO_290 = DENDRO_59*DENDRO_61;
double DENDRO_291 = 2*DENDRO_29;
double DENDRO_292 = DENDRO_291*(-DENDRO_275*DENDRO_290 + 2*DENDRO_289);
double DENDRO_293 = -1.0*DENDRO_201 + DENDRO_209 + DENDRO_210 - DENDRO_220 - DENDRO_231 - DENDRO_232;
double DENDRO_294 = 2*DENDRO_118*DENDRO_293;
double DENDRO_295 = -1.0*DENDRO_255 + DENDRO_260 + DENDRO_261 - DENDRO_266 - DENDRO_269 - DENDRO_270;
double DENDRO_296 = 2*DENDRO_117*DENDRO_295;
double DENDRO_297 = -1.0*DENDRO_236 + DENDRO_241 + DENDRO_243 - DENDRO_248 - DENDRO_251 - DENDRO_252;
double DENDRO_298 = 2*DENDRO_120*DENDRO_297;
double DENDRO_299 = DENDRO_276 + DENDRO_279 + DENDRO_281 - DENDRO_285 + DENDRO_288 - DENDRO_292 + DENDRO_294 + DENDRO_296 + DENDRO_298;
double DENDRO_300 = DENDRO_104*DENDRO_64;
double DENDRO_301 = DENDRO_144*DENDRO_55;
double DENDRO_302 = 4*DENDRO_301;
double DENDRO_303 = DENDRO_160*DENDRO_302;
double DENDRO_304 = 0.25*DENDRO_122;
double DENDRO_305 = DENDRO_144*DENDRO_35;
double DENDRO_306 = 4*DENDRO_305;
double DENDRO_307 = DENDRO_141*DENDRO_306;
double DENDRO_308 = 0.25*DENDRO_128;
double DENDRO_309 = -DENDRO_308;
double DENDRO_310 = 0.75*DENDRO_130;
double DENDRO_311 = 0.25*DENDRO_129;
double DENDRO_312 = DENDRO_306*(DENDRO_309 + DENDRO_310 + DENDRO_311);
double DENDRO_313 = DENDRO_73 + DENDRO_75;
double DENDRO_314 = DENDRO_144*DENDRO_49;
double DENDRO_315 = 4*DENDRO_314;
double DENDRO_316 = DENDRO_179*DENDRO_74;
double DENDRO_317 = 3.0*DENDRO_301;
double DENDRO_318 = DENDRO_173*DENDRO_69;
double DENDRO_319 = 3.0*DENDRO_305;
double DENDRO_320 = 6.0*DENDRO_77;
double DENDRO_321 = pow(chi[pp], -2);
double DENDRO_322 = grad_0_Gt0[pp];
double DENDRO_323 = grad_0_Gt1[pp];
double DENDRO_324 = 4*gt1[pp];
double DENDRO_325 = grad_0_Gt2[pp];
double DENDRO_326 = 4*gt2[pp];
double DENDRO_327 = 0.5*DENDRO_71;
double DENDRO_328 = DENDRO_134*DENDRO_327;
double DENDRO_329 = DENDRO_41*DENDRO_46;
double DENDRO_330 = 4*DENDRO_329;
double DENDRO_331 = DENDRO_134*DENDRO_304;
double DENDRO_332 = 0.5*DENDRO_167;
double DENDRO_333 = DENDRO_141*DENDRO_327;
double DENDRO_334 = DENDRO_41*DENDRO_55;
double DENDRO_335 = 2.0*DENDRO_334;
double DENDRO_336 = DENDRO_35*DENDRO_41;
double DENDRO_337 = 2.0*DENDRO_336;
double DENDRO_338 = DENDRO_41*DENDRO_49;
double DENDRO_339 = 2.0*DENDRO_338;
double DENDRO_340 = DENDRO_141*DENDRO_69;
double DENDRO_341 = DENDRO_122*DENDRO_80;
double DENDRO_342 = DENDRO_134*DENDRO_69;
double DENDRO_343 = DENDRO_131*DENDRO_80;
double DENDRO_344 = 4.0*DENDRO_41;
double DENDRO_345 = DENDRO_31*DENDRO_344;
double DENDRO_346 = DENDRO_29*DENDRO_344;
double DENDRO_347 = 0.25*DENDRO_130;
double DENDRO_348 = 0.75*DENDRO_129;
double DENDRO_349 = 4*DENDRO_144;
double DENDRO_350 = -DENDRO_134*DENDRO_302*(DENDRO_309 + DENDRO_347 + DENDRO_348) + DENDRO_146*(DENDRO_141*DENDRO_332 + DENDRO_331) - DENDRO_184*(DENDRO_167*DENDRO_80 + DENDRO_328) + DENDRO_190*(-2*DENDRO_127*DENDRO_80 + DENDRO_333) - DENDRO_197*(DENDRO_342 + DENDRO_343) + DENDRO_200*(DENDRO_340 + DENDRO_341) - DENDRO_232*DENDRO_349*(DENDRO_68 + DENDRO_70) - DENDRO_280*DENDRO_321 + 4*DENDRO_322*gt0[pp] + DENDRO_323*DENDRO_324 + DENDRO_325*DENDRO_326 + DENDRO_330*grad2_0_2_gt0[pp] + DENDRO_335*grad2_2_2_gt0[pp] + DENDRO_337*grad2_1_1_gt0[pp] + DENDRO_339*grad2_0_0_gt0[pp] - DENDRO_345*grad2_1_2_gt0[pp] - DENDRO_346*grad2_0_1_gt0[pp];
double DENDRO_351 = 3*alpha[pp];
double DENDRO_352 = grad2_2_2_alpha[pp];
double DENDRO_353 = DENDRO_65*gt5[pp];
double DENDRO_354 = DENDRO_219 + DENDRO_353*DENDRO_63;
double DENDRO_355 = 4*DENDRO_82;
double DENDRO_356 = DENDRO_355*DENDRO_41;
double DENDRO_357 = DENDRO_46*DENDRO_60;
double DENDRO_358 = DENDRO_49*DENDRO_61;
double DENDRO_359 = DENDRO_106 - DENDRO_357 - DENDRO_358;
double DENDRO_360 = DENDRO_353*DENDRO_359;
double DENDRO_361 = 4*DENDRO_112;
double DENDRO_362 = DENDRO_361*DENDRO_41;
double DENDRO_363 = DENDRO_262*DENDRO_41;
double DENDRO_364 = DENDRO_263*DENDRO_41;
double DENDRO_365 = DENDRO_264*DENDRO_41;
double DENDRO_366 = 1.0*DENDRO_60;
double DENDRO_367 = DENDRO_41*gt5[pp];
double DENDRO_368 = DENDRO_38 - DENDRO_54;
double DENDRO_369 = DENDRO_105*DENDRO_61 + DENDRO_368*DENDRO_60 + DENDRO_85;
double DENDRO_370 = 0.5*DENDRO_369;
double DENDRO_371 = DENDRO_367*DENDRO_370;
double DENDRO_372 = DENDRO_64*(-DENDRO_366 + DENDRO_371);
double DENDRO_373 = 4*DENDRO_94;
double DENDRO_374 = -DENDRO_273;
double DENDRO_375 = DENDRO_212*DENDRO_368;
double DENDRO_376 = -DENDRO_214;
double DENDRO_377 = DENDRO_216 + DENDRO_376;
double DENDRO_378 = DENDRO_31*DENDRO_377;
double DENDRO_379 = -DENDRO_147;
double DENDRO_380 = DENDRO_149 + DENDRO_379;
double DENDRO_381 = DENDRO_105*DENDRO_380;
double DENDRO_382 = -DENDRO_33 + DENDRO_34;
double DENDRO_383 = DENDRO_213 + DENDRO_29*DENDRO_380 + DENDRO_377*DENDRO_382;
double DENDRO_384 = DENDRO_105*DENDRO_212 + DENDRO_107*DENDRO_380 + DENDRO_29*DENDRO_377;
double DENDRO_385 = 0.5*DENDRO_218;
double DENDRO_386 = DENDRO_134*DENDRO_385;
double DENDRO_387 = -DENDRO_386;
double DENDRO_388 = DENDRO_137*DENDRO_219;
double DENDRO_389 = 0.5*DENDRO_151;
double DENDRO_390 = -DENDRO_235*DENDRO_389;
double DENDRO_391 = 2*DENDRO_76;
double DENDRO_392 = DENDRO_160*DENDRO_211;
double DENDRO_393 = 0.25*DENDRO_392;
double DENDRO_394 = DENDRO_132*DENDRO_265;
double DENDRO_395 = DENDRO_240*DENDRO_389;
double DENDRO_396 = -DENDRO_395;
double DENDRO_397 = DENDRO_137*DENDRO_247;
double DENDRO_398 = DENDRO_211*DENDRO_259;
double DENDRO_399 = 0.25*DENDRO_398;
double DENDRO_400 = DENDRO_204*DENDRO_265;
double DENDRO_401 = DENDRO_240*DENDRO_76;
double DENDRO_402 = DENDRO_167*DENDRO_235;
double DENDRO_403 = 0.25*DENDRO_402;
double DENDRO_404 = 0.25*DENDRO_204;
double DENDRO_405 = DENDRO_160*DENDRO_404;
double DENDRO_406 = DENDRO_165*DENDRO_259;
double DENDRO_407 = DENDRO_174*DENDRO_240;
double DENDRO_408 = 0.5*DENDRO_137;
double DENDRO_409 = DENDRO_167*DENDRO_247;
double DENDRO_410 = 2.0*DENDRO_145;
double DENDRO_411 = DENDRO_144*DENDRO_293;
double DENDRO_412 = 2.0*DENDRO_411;
double DENDRO_413 = DENDRO_144*DENDRO_295;
double DENDRO_414 = 2.0*DENDRO_413;
double DENDRO_415 = DENDRO_144*DENDRO_297;
double DENDRO_416 = 2.0*DENDRO_415;
double DENDRO_417 = DENDRO_247*DENDRO_74;
double DENDRO_418 = -DENDRO_276 - DENDRO_279 - DENDRO_281 + DENDRO_285 - DENDRO_288 + DENDRO_292 - DENDRO_294 - DENDRO_296 - DENDRO_298;
double DENDRO_419 = DENDRO_418*DENDRO_64;
double DENDRO_420 = DENDRO_248*DENDRO_349;
double DENDRO_421 = DENDRO_220*DENDRO_349;
double DENDRO_422 = -DENDRO_311;
double DENDRO_423 = DENDRO_306*(DENDRO_308 + DENDRO_310 + DENDRO_422);
double DENDRO_424 = DENDRO_315*(-DENDRO_174 + DENDRO_73);
double DENDRO_425 = DENDRO_204*DENDRO_259;
double DENDRO_426 = DENDRO_132*DENDRO_160;
double DENDRO_427 = 3.0*DENDRO_314;
double DENDRO_428 = 6.0*DENDRO_144;
double DENDRO_429 = grad_2_Gt0[pp];
double DENDRO_430 = grad_2_Gt1[pp];
double DENDRO_431 = 4*gt4[pp];
double DENDRO_432 = grad_2_Gt2[pp];
double DENDRO_433 = -DENDRO_208*DENDRO_385;
double DENDRO_434 = DENDRO_134*DENDRO_227;
double DENDRO_435 = DENDRO_131*DENDRO_208;
double DENDRO_436 = 0.25*DENDRO_435;
double DENDRO_437 = 0.25*DENDRO_202;
double DENDRO_438 = DENDRO_134*DENDRO_437;
double DENDRO_439 = DENDRO_204*DENDRO_208;
double DENDRO_440 = DENDRO_202*DENDRO_219;
double DENDRO_441 = DENDRO_134*DENDRO_204;
double DENDRO_442 = DENDRO_131*DENDRO_219;
double DENDRO_443 = 0.75*DENDRO_128;
double DENDRO_444 = -DENDRO_134*DENDRO_315*(DENDRO_347 + DENDRO_422 + DENDRO_443) + DENDRO_146*(2*DENDRO_219*DENDRO_227 + DENDRO_433) + DENDRO_190*(DENDRO_434 + DENDRO_436) + DENDRO_190*(DENDRO_208*DENDRO_408 + DENDRO_438) - DENDRO_197*(DENDRO_441 + DENDRO_442) - DENDRO_208*DENDRO_306*(DENDRO_225 - DENDRO_437) - DENDRO_274*DENDRO_321 + DENDRO_326*DENDRO_429 + DENDRO_330*grad2_0_2_gt5[pp] + DENDRO_335*grad2_2_2_gt5[pp] + DENDRO_337*grad2_1_1_gt5[pp] + DENDRO_339*grad2_0_0_gt5[pp] - DENDRO_345*grad2_1_2_gt5[pp] - DENDRO_346*grad2_0_1_gt5[pp] + DENDRO_410*(DENDRO_439 + DENDRO_440) + DENDRO_430*DENDRO_431 + 4*DENDRO_432*gt5[pp];
double DENDRO_445 = grad2_1_1_alpha[pp];
double DENDRO_446 = DENDRO_65*gt3[pp];
double DENDRO_447 = DENDRO_373*DENDRO_41;
double DENDRO_448 = 1.0*DENDRO_59;
double DENDRO_449 = -DENDRO_448;
double DENDRO_450 = DENDRO_382*DENDRO_59 + DENDRO_62;
double DENDRO_451 = DENDRO_41*gt3[pp];
double DENDRO_452 = 0.5*DENDRO_451;
double DENDRO_453 = DENDRO_450*DENDRO_452;
double DENDRO_454 = DENDRO_228*DENDRO_41;
double DENDRO_455 = -DENDRO_223*DENDRO_41 - DENDRO_229*DENDRO_41 + DENDRO_454;
double DENDRO_456 = -DENDRO_277;
double DENDRO_457 = -DENDRO_123;
double DENDRO_458 = DENDRO_125 + DENDRO_457;
double DENDRO_459 = DENDRO_105*DENDRO_458 + DENDRO_227*DENDRO_368 + DENDRO_267;
double DENDRO_460 = DENDRO_222*DENDRO_382;
double DENDRO_461 = DENDRO_29*DENDRO_458;
double DENDRO_462 = DENDRO_105*DENDRO_227 + DENDRO_107*DENDRO_458 + DENDRO_249;
double DENDRO_463 = DENDRO_240*DENDRO_71;
double DENDRO_464 = DENDRO_167*DENDRO_242;
double DENDRO_465 = 0.25*DENDRO_464;
double DENDRO_466 = DENDRO_155*DENDRO_218;
double DENDRO_467 = DENDRO_137*DENDRO_259;
double DENDRO_468 = 0.25*DENDRO_467;
double DENDRO_469 = DENDRO_155*DENDRO_404;
double DENDRO_470 = 0.5*DENDRO_131;
double DENDRO_471 = DENDRO_180*DENDRO_240;
double DENDRO_472 = 0.5*DENDRO_127;
double DENDRO_473 = DENDRO_240*DENDRO_472;
double DENDRO_474 = -DENDRO_473;
double DENDRO_475 = DENDRO_131*DENDRO_250;
double DENDRO_476 = 0.5*DENDRO_227;
double DENDRO_477 = DENDRO_259*DENDRO_476;
double DENDRO_478 = 2*DENDRO_218*DENDRO_268;
double DENDRO_479 = DENDRO_208*DENDRO_221;
double DENDRO_480 = 0.25*DENDRO_479;
double DENDRO_481 = DENDRO_202*DENDRO_230;
double DENDRO_482 = DENDRO_155*DENDRO_476;
double DENDRO_483 = DENDRO_131*DENDRO_268;
double DENDRO_484 = -DENDRO_242*DENDRO_472;
double DENDRO_485 = 2*DENDRO_250*DENDRO_71;
double DENDRO_486 = DENDRO_141*DENDRO_221;
double DENDRO_487 = 0.25*DENDRO_486;
double DENDRO_488 = DENDRO_122*DENDRO_230;
double DENDRO_489 = DENDRO_167*DENDRO_250;
double DENDRO_490 = DENDRO_204*DENDRO_268;
double DENDRO_491 = DENDRO_137*DENDRO_268;
double DENDRO_492 = DENDRO_250*DENDRO_69;
double DENDRO_493 = DENDRO_259*DENDRO_302;
double DENDRO_494 = -DENDRO_347;
double DENDRO_495 = DENDRO_302*(DENDRO_308 + DENDRO_348 + DENDRO_494);
double DENDRO_496 = DENDRO_251*DENDRO_349;
double DENDRO_497 = DENDRO_315*(-DENDRO_180 + DENDRO_68);
double DENDRO_498 = DENDRO_315*(DENDRO_311 + DENDRO_443 + DENDRO_494);
double DENDRO_499 = grad_1_Gt0[pp];
double DENDRO_500 = grad_1_Gt1[pp];
double DENDRO_501 = grad_1_Gt2[pp];
double DENDRO_502 = DENDRO_141*DENDRO_437;
double DENDRO_503 = DENDRO_208*DENDRO_304;
double DENDRO_504 = DENDRO_202*DENDRO_208;
double DENDRO_505 = DENDRO_122*DENDRO_141;
double DENDRO_506 = -DENDRO_184*(DENDRO_123*DENDRO_208 + DENDRO_502) - DENDRO_184*(DENDRO_141*DENDRO_226 + DENDRO_503) - DENDRO_269*DENDRO_349*(DENDRO_225 + DENDRO_226) - DENDRO_278*DENDRO_321 - DENDRO_317*DENDRO_504 + DENDRO_324*DENDRO_499 + DENDRO_330*grad2_0_2_gt3[pp] + DENDRO_335*grad2_2_2_gt3[pp] + DENDRO_337*grad2_1_1_gt3[pp] + DENDRO_339*grad2_0_0_gt3[pp] - DENDRO_345*grad2_1_2_gt3[pp] - DENDRO_346*grad2_0_1_gt3[pp] - DENDRO_427*DENDRO_505 + DENDRO_431*DENDRO_501 + 4*DENDRO_500*gt3[pp];
double DENDRO_507 = DENDRO_105*DENDRO_78 + DENDRO_368*DENDRO_76 + DENDRO_92;
double DENDRO_508 = DENDRO_382*DENDRO_71 + DENDRO_79;
double DENDRO_509 = DENDRO_107*DENDRO_78;
double DENDRO_510 = DENDRO_105*DENDRO_76;
double DENDRO_511 = DENDRO_160*DENDRO_191;
double DENDRO_512 = DENDRO_235*DENDRO_77;
double DENDRO_513 = 0.25*DENDRO_512;
double DENDRO_514 = DENDRO_119*DENDRO_74;
double DENDRO_515 = DENDRO_155*DENDRO_165;
double DENDRO_516 = DENDRO_180*DENDRO_235;
double DENDRO_517 = DENDRO_174*DENDRO_242;
double DENDRO_518 = DENDRO_155*DENDRO_191;
double DENDRO_519 = DENDRO_242*DENDRO_77;
double DENDRO_520 = 0.25*DENDRO_519;
double DENDRO_521 = DENDRO_119*DENDRO_69;
double DENDRO_522 = DENDRO_116*DENDRO_137;
double DENDRO_523 = DENDRO_116*DENDRO_132;
double DENDRO_524 = DENDRO_235*DENDRO_74;
double DENDRO_525 = DENDRO_242*DENDRO_69;
double DENDRO_526 = DENDRO_259*DENDRO_74;
double DENDRO_527 = DENDRO_132*DENDRO_155;
double DENDRO_528 = DENDRO_161 + DENDRO_527;
double DENDRO_529 = DENDRO_208*DENDRO_69;
double DENDRO_530 = DENDRO_137*DENDRO_141;
double DENDRO_531 = DENDRO_122*DENDRO_134;
double DENDRO_532 = DENDRO_530 + DENDRO_531;
double DENDRO_533 = DENDRO_165*DENDRO_235;
double DENDRO_534 = DENDRO_151*DENDRO_265;
double DENDRO_535 = -DENDRO_534;
double DENDRO_536 = DENDRO_219*DENDRO_332 + 0.25*DENDRO_441;
double DENDRO_537 = DENDRO_137*DENDRO_242;
double DENDRO_538 = 0.25*DENDRO_537;
double DENDRO_539 = DENDRO_141*DENDRO_476;
double DENDRO_540 = -DENDRO_208*DENDRO_472;
double DENDRO_541 = DENDRO_539 + DENDRO_540;
double DENDRO_542 = DENDRO_119*DENDRO_76;
double DENDRO_543 = DENDRO_160*DENDRO_174;
double DENDRO_544 = 0.25*DENDRO_342 + DENDRO_408*DENDRO_80;
double DENDRO_545 = 0.25*DENDRO_134;
double DENDRO_546 = DENDRO_167*DENDRO_545 + DENDRO_214*DENDRO_80;
double DENDRO_547 = DENDRO_116*DENDRO_212;
double DENDRO_548 = -DENDRO_160*DENDRO_389;
double DENDRO_549 = DENDRO_131*DENDRO_545;
double DENDRO_550 = DENDRO_137*DENDRO_545 + DENDRO_219*DENDRO_70;
double DENDRO_551 = DENDRO_174*DENDRO_235;
double DENDRO_552 = DENDRO_247*DENDRO_78;
double DENDRO_553 = DENDRO_191*DENDRO_235 + DENDRO_552;
double DENDRO_554 = DENDRO_127*DENDRO_219;
double DENDRO_555 = 0.5*DENDRO_434;
double DENDRO_556 = -DENDRO_554 + DENDRO_555;
double DENDRO_557 = DENDRO_167*DENDRO_208;
double DENDRO_558 = 0.25*DENDRO_557;
double DENDRO_559 = DENDRO_123*DENDRO_219;
double DENDRO_560 = DENDRO_141*DENDRO_404 + DENDRO_559;
double DENDRO_561 = 0.25*DENDRO_137;
double DENDRO_562 = DENDRO_247*DENDRO_70;
double DENDRO_563 = DENDRO_235*DENDRO_561 + DENDRO_562;
double DENDRO_564 = DENDRO_265*DENDRO_332;
double DENDRO_565 = DENDRO_405 + DENDRO_406;
double DENDRO_566 = DENDRO_259*DENDRO_389;
double DENDRO_567 = -DENDRO_566;
double DENDRO_568 = 0.25*DENDRO_155;
double DENDRO_569 = DENDRO_211*DENDRO_568;
double DENDRO_570 = DENDRO_265*DENDRO_408 + DENDRO_569;
double DENDRO_571 = DENDRO_165*DENDRO_242 + DENDRO_562;
double DENDRO_572 = DENDRO_208*DENDRO_327;
double DENDRO_573 = -0.5*DENDRO_135 + DENDRO_227*DENDRO_80;
double DENDRO_574 = DENDRO_116*DENDRO_214;
double DENDRO_575 = DENDRO_191*DENDRO_259;
double DENDRO_576 = 0.25*DENDRO_167;
double DENDRO_577 = DENDRO_160*DENDRO_576;
double DENDRO_578 = 0.25*DENDRO_240;
double DENDRO_579 = DENDRO_578*DENDRO_77;
double DENDRO_580 = DENDRO_191*DENDRO_242;
double DENDRO_581 = DENDRO_579 + DENDRO_580;
double DENDRO_582 = DENDRO_119*DENDRO_408 + DENDRO_516;
double DENDRO_583 = DENDRO_167*DENDRO_259;
double DENDRO_584 = DENDRO_155*DENDRO_204;
double DENDRO_585 = DENDRO_467 + DENDRO_584;
double DENDRO_586 = 1.0*DENDRO_305;
double DENDRO_587 = -DENDRO_286;
double DENDRO_588 = DENDRO_132*DENDRO_368;
double DENDRO_589 = DENDRO_105*DENDRO_74;
double DENDRO_590 = 0.5*DENDRO_117;
double DENDRO_591 = DENDRO_131*DENDRO_382 + DENDRO_133;
double DENDRO_592 = 0.5*DENDRO_118;
double DENDRO_593 = DENDRO_105*DENDRO_132;
double DENDRO_594 = DENDRO_107*DENDRO_74;
double DENDRO_595 = 0.5*DENDRO_120;
double DENDRO_596 = 2.0*DENDRO_429;
double DENDRO_597 = DENDRO_596*gt0[pp];
double DENDRO_598 = 2.0*DENDRO_430;
double DENDRO_599 = DENDRO_598*gt1[pp];
double DENDRO_600 = 2.0*gt2[pp];
double DENDRO_601 = DENDRO_322*DENDRO_600;
double DENDRO_602 = DENDRO_432*DENDRO_600;
double DENDRO_603 = 2.0*gt4[pp];
double DENDRO_604 = DENDRO_323*DENDRO_603;
double DENDRO_605 = 2.0*gt5[pp];
double DENDRO_606 = DENDRO_325*DENDRO_605;
double DENDRO_607 = DENDRO_321*DENDRO_60;
double DENDRO_608 = -DENDRO_607*DENDRO_61;
double DENDRO_609 = DENDRO_330*grad2_0_2_gt2[pp];
double DENDRO_610 = DENDRO_335*grad2_2_2_gt2[pp];
double DENDRO_611 = DENDRO_337*grad2_1_1_gt2[pp];
double DENDRO_612 = DENDRO_339*grad2_0_0_gt2[pp];
double DENDRO_613 = -DENDRO_345*grad2_1_2_gt2[pp];
double DENDRO_614 = -DENDRO_346*grad2_0_1_gt2[pp];
double DENDRO_615 = DENDRO_41*gt2[pp];
double DENDRO_616 = -DENDRO_121*(DENDRO_587 + DENDRO_590*(DENDRO_159 + DENDRO_588 + DENDRO_589) + DENDRO_591*DENDRO_592 + DENDRO_595*(DENDRO_178 + DENDRO_593 + DENDRO_594)) + DENDRO_130*DENDRO_412 + DENDRO_148*DENDRO_414 + DENDRO_416*DENDRO_72 + DENDRO_419*DENDRO_615 + DENDRO_597 + DENDRO_599 + DENDRO_601 + DENDRO_602 + DENDRO_604 + DENDRO_606 + DENDRO_608 + DENDRO_609 + DENDRO_610 + DENDRO_611 + DENDRO_612 + DENDRO_613 + DENDRO_614;
double DENDRO_617 = grad2_0_2_alpha[pp];
double DENDRO_618 = DENDRO_157*DENDRO_41;
double DENDRO_619 = DENDRO_158*DENDRO_41;
double DENDRO_620 = DENDRO_159*DENDRO_41;
double DENDRO_621 = -DENDRO_61;
double DENDRO_622 = DENDRO_369*DENDRO_615;
double DENDRO_623 = DENDRO_64*(DENDRO_621 + DENDRO_622);
double DENDRO_624 = 2.0*DENDRO_94;
double DENDRO_625 = DENDRO_176*DENDRO_41;
double DENDRO_626 = DENDRO_177*DENDRO_41;
double DENDRO_627 = DENDRO_178*DENDRO_41;
double DENDRO_628 = -DENDRO_60;
double DENDRO_629 = DENDRO_108*DENDRO_615;
double DENDRO_630 = DENDRO_64*(DENDRO_628 + DENDRO_629);
double DENDRO_631 = 2.0*DENDRO_112;
double DENDRO_632 = 2.0*DENDRO_82;
double DENDRO_633 = DENDRO_64*gt2[pp];
double DENDRO_634 = DENDRO_41*(DENDRO_134 + DENDRO_63*DENDRO_633);
double DENDRO_635 = -4*DENDRO_617 + DENDRO_624*(-DENDRO_618 - DENDRO_619 + DENDRO_620 + DENDRO_623) + DENDRO_631*(-DENDRO_625 - DENDRO_626 + DENDRO_627 + DENDRO_630) + DENDRO_632*DENDRO_634;
double DENDRO_636 = DENDRO_240*DENDRO_74;
double DENDRO_637 = DENDRO_132*DENDRO_242 + DENDRO_636;
double DENDRO_638 = 0.5*DENDRO_392;
double DENDRO_639 = 0.25*DENDRO_583;
double DENDRO_640 = DENDRO_160*DENDRO_165 + DENDRO_547;
double DENDRO_641 = DENDRO_141*DENDRO_385;
double DENDRO_642 = -DENDRO_641;
double DENDRO_643 = -DENDRO_242*DENDRO_389;
double DENDRO_644 = DENDRO_405 + DENDRO_569;
double DENDRO_645 = DENDRO_515 + DENDRO_574;
double DENDRO_646 = DENDRO_174*DENDRO_259;
double DENDRO_647 = DENDRO_119*DENDRO_332;
double DENDRO_648 = DENDRO_240*DENDRO_69;
double DENDRO_649 = DENDRO_537 + DENDRO_648;
double DENDRO_650 = DENDRO_134*DENDRO_202;
double DENDRO_651 = DENDRO_141*DENDRO_204 + DENDRO_650;
double DENDRO_652 = 0.25*DENDRO_530;
double DENDRO_653 = DENDRO_226*DENDRO_80;
double DENDRO_654 = DENDRO_180*DENDRO_208 + DENDRO_653;
double DENDRO_655 = DENDRO_145*(DENDRO_557 + DENDRO_651) - DENDRO_184*(DENDRO_546 + DENDRO_549) - DENDRO_184*(-DENDRO_218*DENDRO_80 + DENDRO_550) + DENDRO_190*(DENDRO_143 + DENDRO_573) + DENDRO_190*(DENDRO_652 + DENDRO_654) - DENDRO_302*(DENDRO_387 + DENDRO_536) - DENDRO_306*(DENDRO_502 + DENDRO_541) - DENDRO_315*(0.5*DENDRO_343 + DENDRO_544);
double DENDRO_656 = DENDRO_160*DENDRO_202;
double DENDRO_657 = DENDRO_122*DENDRO_235;
double DENDRO_658 = 0.25*DENDRO_439;
double DENDRO_659 = DENDRO_218*DENDRO_265;
double DENDRO_660 = -DENDRO_659;
double DENDRO_661 = DENDRO_165*DENDRO_240 + DENDRO_247*DENDRO_470;
double DENDRO_662 = DENDRO_227*DENDRO_230;
double DENDRO_663 = DENDRO_259*DENDRO_437;
double DENDRO_664 = DENDRO_250*DENDRO_408;
double DENDRO_665 = DENDRO_240*DENDRO_304 + DENDRO_664;
double DENDRO_666 = DENDRO_235*DENDRO_327;
double DENDRO_667 = DENDRO_580 + DENDRO_666;
double DENDRO_668 = DENDRO_247*DENDRO_71 + 0.5*DENDRO_401;
double DENDRO_669 = DENDRO_131*DENDRO_235;
double DENDRO_670 = 0.25*DENDRO_669;
double DENDRO_671 = DENDRO_208*DENDRO_561 + DENDRO_559;
double DENDRO_672 = DENDRO_265*DENDRO_470;
double DENDRO_673 = DENDRO_160*DENDRO_385;
double DENDRO_674 = -DENDRO_673;
double DENDRO_675 = DENDRO_212*DENDRO_268;
double DENDRO_676 = -DENDRO_259*DENDRO_385 + DENDRO_675;
double DENDRO_677 = DENDRO_147*DENDRO_250;
double DENDRO_678 = DENDRO_131*DENDRO_578 + DENDRO_677;
double DENDRO_679 = DENDRO_240*DENDRO_576;
double DENDRO_680 = DENDRO_123*DENDRO_247 + DENDRO_240*DENDRO_561;
double DENDRO_681 = DENDRO_208*DENDRO_476;
double DENDRO_682 = DENDRO_219*DENDRO_222;
double DENDRO_683 = DENDRO_208*DENDRO_437 + DENDRO_682;
double DENDRO_684 = -DENDRO_235*DENDRO_472;
double DENDRO_685 = DENDRO_250*DENDRO_76;
double DENDRO_686 = 0.5*DENDRO_463 + DENDRO_685;
double DENDRO_687 = DENDRO_160*DENDRO_476;
double DENDRO_688 = 0.25*DENDRO_131;
double DENDRO_689 = DENDRO_147*DENDRO_268;
double DENDRO_690 = DENDRO_259*DENDRO_688 + DENDRO_689;
double DENDRO_691 = DENDRO_221*DENDRO_545;
double DENDRO_692 = DENDRO_503 + DENDRO_691;
double DENDRO_693 = DENDRO_230*DENDRO_408;
double DENDRO_694 = DENDRO_131*DENDRO_160;
double DENDRO_695 = 1.0*DENDRO_314;
double DENDRO_696 = -DENDRO_282;
double DENDRO_697 = DENDRO_204*DENDRO_368;
double DENDRO_698 = DENDRO_105*DENDRO_167;
double DENDRO_699 = DENDRO_202*DENDRO_382;
double DENDRO_700 = DENDRO_105*DENDRO_204 + DENDRO_107*DENDRO_167 + DENDRO_237;
double DENDRO_701 = DENDRO_41*gt4[pp];
double DENDRO_702 = DENDRO_330*grad2_0_2_gt4[pp] + DENDRO_335*grad2_2_2_gt4[pp] + DENDRO_337*grad2_1_1_gt4[pp] + DENDRO_339*grad2_0_0_gt4[pp] - DENDRO_345*grad2_1_2_gt4[pp] - DENDRO_346*grad2_0_1_gt4[pp] + DENDRO_432*DENDRO_603 + DENDRO_499*DENDRO_600 + DENDRO_500*DENDRO_603 + DENDRO_501*DENDRO_605 - DENDRO_59*DENDRO_607 + DENDRO_596*gt1[pp] + DENDRO_598*gt3[pp];
double DENDRO_703 = -DENDRO_121*(DENDRO_590*(DENDRO_256 + DENDRO_697 + DENDRO_698) + DENDRO_592*(DENDRO_207 + DENDRO_699) + DENDRO_595*DENDRO_700 + DENDRO_696) + DENDRO_128*DENDRO_416 + DENDRO_215*DENDRO_414 + DENDRO_224*DENDRO_412 + DENDRO_419*DENDRO_701 + DENDRO_702;
double DENDRO_704 = grad2_1_2_alpha[pp];
double DENDRO_705 = DENDRO_256*DENDRO_41;
double DENDRO_706 = DENDRO_257*DENDRO_41;
double DENDRO_707 = DENDRO_258*DENDRO_41;
double DENDRO_708 = -DENDRO_59;
double DENDRO_709 = DENDRO_369*DENDRO_701;
double DENDRO_710 = DENDRO_64*(DENDRO_708 + DENDRO_709);
double DENDRO_711 = DENDRO_450*DENDRO_701;
double DENDRO_712 = DENDRO_205*DENDRO_41 + DENDRO_206*DENDRO_41;
double DENDRO_713 = -DENDRO_203*DENDRO_41 + DENDRO_712;
double DENDRO_714 = DENDRO_64*gt4[pp];
double DENDRO_715 = DENDRO_359*DENDRO_714;
double DENDRO_716 = DENDRO_41*DENDRO_631*(DENDRO_240 + DENDRO_715) + DENDRO_624*(DENDRO_705 - DENDRO_706 - DENDRO_707 + DENDRO_710) + DENDRO_632*(DENDRO_64*(DENDRO_628 + DENDRO_711) + DENDRO_713) - 4*DENDRO_704;
double DENDRO_717 = 0.5*DENDRO_398;
double DENDRO_718 = 1.0*DENDRO_490;
double DENDRO_719 = 0.5*DENDRO_489;
double DENDRO_720 = 0.25*DENDRO_694;
double DENDRO_721 = DENDRO_406 + DENDRO_569;
double DENDRO_722 = DENDRO_151*DENDRO_250;
double DENDRO_723 = DENDRO_681 + DENDRO_682;
double DENDRO_724 = DENDRO_259*DENDRO_404;
double DENDRO_725 = DENDRO_250*DENDRO_75;
double DENDRO_726 = DENDRO_235*DENDRO_304 + DENDRO_725;
double DENDRO_727 = DENDRO_160*DENDRO_437 + DENDRO_689;
double DENDRO_728 = DENDRO_502 + DENDRO_503;
double DENDRO_729 = DENDRO_230*DENDRO_470 + DENDRO_691;
double DENDRO_730 = 1.0*DENDRO_183;
double DENDRO_731 = -DENDRO_184*(DENDRO_642 + DENDRO_671) - DENDRO_302*(DENDRO_219*DENDRO_226 + DENDRO_433 + DENDRO_658) - DENDRO_695*(DENDRO_142 + DENDRO_532) - DENDRO_730*(DENDRO_435 + DENDRO_651);
double DENDRO_732 = DENDRO_567 + DENDRO_674;
double DENDRO_733 = 0.5*DENDRO_486;
double DENDRO_734 = DENDRO_127*DENDRO_230;
double DENDRO_735 = -DENDRO_734;
double DENDRO_736 = DENDRO_268*DENDRO_332;
double DENDRO_737 = DENDRO_155*DENDRO_437 + DENDRO_736;
double DENDRO_738 = DENDRO_242*DENDRO_304;
double DENDRO_739 = DENDRO_250*DENDRO_70;
double DENDRO_740 = DENDRO_119*DENDRO_71;
double DENDRO_741 = DENDRO_116*DENDRO_470 + DENDRO_155*DENDRO_174;
double DENDRO_742 = DENDRO_116*DENDRO_218;
double DENDRO_743 = 0.5*DENDRO_156;
double DENDRO_744 = -DENDRO_742 - DENDRO_743;
double DENDRO_745 = DENDRO_119*DENDRO_470 + DENDRO_517;
double DENDRO_746 = DENDRO_579 + DENDRO_666;
double DENDRO_747 = DENDRO_151*DENDRO_268;
double DENDRO_748 = 0.5*DENDRO_466;
double DENDRO_749 = -DENDRO_747 - DENDRO_748;
double DENDRO_750 = DENDRO_242*DENDRO_688 + DENDRO_725;
double DENDRO_751 = DENDRO_230*DENDRO_332;
double DENDRO_752 = DENDRO_502 + DENDRO_691;
double DENDRO_753 = DENDRO_268*DENDRO_75;
double DENDRO_754 = DENDRO_131*DENDRO_568 + DENDRO_753;
double DENDRO_755 = DENDRO_250*DENDRO_78;
double DENDRO_756 = DENDRO_242*DENDRO_327 + DENDRO_755;
double DENDRO_757 = DENDRO_155*DENDRO_561;
double DENDRO_758 = DENDRO_116*DENDRO_226 + DENDRO_167*DENDRO_568;
double DENDRO_759 = DENDRO_222*DENDRO_80;
double DENDRO_760 = DENDRO_141*DENDRO_304 + DENDRO_759;
double DENDRO_761 = 1.0*DENDRO_301;
double DENDRO_762 = -DENDRO_289;
double DENDRO_763 = DENDRO_105*DENDRO_69 + DENDRO_137*DENDRO_368 + DENDRO_152;
double DENDRO_764 = DENDRO_122*DENDRO_382;
double DENDRO_765 = DENDRO_107*DENDRO_69;
double DENDRO_766 = DENDRO_105*DENDRO_137;
double DENDRO_767 = 2.0*DENDRO_499*gt0[pp];
double DENDRO_768 = 2.0*gt1[pp];
double DENDRO_769 = DENDRO_322*DENDRO_768;
double DENDRO_770 = DENDRO_500*DENDRO_768;
double DENDRO_771 = DENDRO_501*DENDRO_600;
double DENDRO_772 = 2.0*DENDRO_323*gt3[pp];
double DENDRO_773 = DENDRO_325*DENDRO_603;
double DENDRO_774 = -DENDRO_290*DENDRO_321;
double DENDRO_775 = DENDRO_330*grad2_0_2_gt1[pp];
double DENDRO_776 = DENDRO_335*grad2_2_2_gt1[pp];
double DENDRO_777 = DENDRO_337*grad2_1_1_gt1[pp];
double DENDRO_778 = DENDRO_339*grad2_0_0_gt1[pp];
double DENDRO_779 = -DENDRO_345*grad2_1_2_gt1[pp];
double DENDRO_780 = -DENDRO_346*grad2_0_1_gt1[pp];
double DENDRO_781 = DENDRO_41*gt1[pp];
double DENDRO_782 = -DENDRO_121*(DENDRO_590*DENDRO_763 + DENDRO_592*(DENDRO_140 + DENDRO_764) + DENDRO_595*(DENDRO_170 + DENDRO_765 + DENDRO_766) + DENDRO_762) + DENDRO_124*DENDRO_412 + DENDRO_129*DENDRO_414 + DENDRO_416*DENDRO_67 + DENDRO_419*DENDRO_781 + DENDRO_767 + DENDRO_769 + DENDRO_770 + DENDRO_771 + DENDRO_772 + DENDRO_773 + DENDRO_774 + DENDRO_775 + DENDRO_776 + DENDRO_777 + DENDRO_778 + DENDRO_779 + DENDRO_780;
double DENDRO_783 = 0.25*DENDRO_340;
double DENDRO_784 = DENDRO_141*DENDRO_576 + DENDRO_653;
double DENDRO_785 = -DENDRO_141*DENDRO_472;
double DENDRO_786 = DENDRO_146*(DENDRO_540 + DENDRO_752) - DENDRO_184*(DENDRO_331 + DENDRO_654) - DENDRO_184*(DENDRO_331 + DENDRO_784) + DENDRO_190*(DENDRO_760 + DENDRO_785) - DENDRO_302*(DENDRO_134*DENDRO_226 + DENDRO_558) - DENDRO_315*(1.0*DENDRO_341 + DENDRO_783);
double DENDRO_787 = grad2_0_1_alpha[pp];
double DENDRO_788 = DENDRO_450*DENDRO_781;
double DENDRO_789 = DENDRO_138*DENDRO_41 + DENDRO_139*DENDRO_41;
double DENDRO_790 = -DENDRO_136*DENDRO_41 + DENDRO_789;
double DENDRO_791 = DENDRO_170*DENDRO_41;
double DENDRO_792 = DENDRO_171*DENDRO_41;
double DENDRO_793 = DENDRO_172*DENDRO_41;
double DENDRO_794 = DENDRO_108*DENDRO_781;
double DENDRO_795 = DENDRO_64*(DENDRO_708 + DENDRO_794);
double DENDRO_796 = DENDRO_64*gt1[pp];
double DENDRO_797 = DENDRO_796*DENDRO_88;
double DENDRO_798 = DENDRO_41*DENDRO_624*(DENDRO_155 + DENDRO_797) + DENDRO_631*(DENDRO_791 - DENDRO_792 - DENDRO_793 + DENDRO_795) + DENDRO_632*(DENDRO_64*(DENDRO_621 + DENDRO_788) + DENDRO_790) - 4*DENDRO_787;
double DENDRO_799 = DENDRO_180*DENDRO_242;
double DENDRO_800 = -DENDRO_29*(DENDRO_798 + alpha[pp]*(DENDRO_145*(DENDRO_464 + DENDRO_648 + DENDRO_657) + DENDRO_145*(DENDRO_583 + DENDRO_584 + DENDRO_656) + DENDRO_146*(DENDRO_684 + DENDRO_750) + DENDRO_146*(DENDRO_687 + DENDRO_749) + DENDRO_146*(DENDRO_751 + DENDRO_752) - DENDRO_184*(DENDRO_162 + DENDRO_744) - DENDRO_184*(DENDRO_516 + DENDRO_745) - DENDRO_184*(DENDRO_647 + DENDRO_746) - DENDRO_184*(DENDRO_574 + DENDRO_646 + DENDRO_720) + DENDRO_190*(DENDRO_757 + DENDRO_758) + DENDRO_190*(DENDRO_116*DENDRO_227 + DENDRO_754) + DENDRO_190*(-DENDRO_119*DENDRO_127 + DENDRO_756) + DENDRO_190*(DENDRO_230*DENDRO_70 + DENDRO_760) + DENDRO_200*(DENDRO_119*DENDRO_122 + DENDRO_525) - DENDRO_302*(DENDRO_405 + DENDRO_732) - DENDRO_306*(DENDRO_482 + DENDRO_737) - DENDRO_306*(DENDRO_733 + DENDRO_735) - DENDRO_306*(DENDRO_484 + DENDRO_738 + DENDRO_739) - DENDRO_315*(0.5*DENDRO_522 + DENDRO_741) - DENDRO_315*(DENDRO_119*DENDRO_70 + DENDRO_520 + DENDRO_740) - DENDRO_761*(DENDRO_402 + DENDRO_636 + DENDRO_669) + DENDRO_782 + DENDRO_786)) - DENDRO_29*(DENDRO_798 + alpha[pp]*(DENDRO_146*(DENDRO_468 + DENDRO_749) + DENDRO_146*(DENDRO_471 + DENDRO_726) + DENDRO_146*(DENDRO_471 + DENDRO_750) + DENDRO_146*(DENDRO_540 + DENDRO_729) + DENDRO_146*(DENDRO_639 + DENDRO_727) + DENDRO_146*(DENDRO_728 + DENDRO_751) - DENDRO_184*(DENDRO_517 + DENDRO_746) - DENDRO_184*(DENDRO_572 + DENDRO_784) - DENDRO_184*(DENDRO_575 + DENDRO_744) - DENDRO_184*(DENDRO_579 + DENDRO_745) + DENDRO_190*(DENDRO_754 + DENDRO_757) + DENDRO_190*(DENDRO_756 + DENDRO_799) + DENDRO_190*(DENDRO_268*DENDRO_76 + DENDRO_758) + DENDRO_190*(DENDRO_119*DENDRO_123 + DENDRO_755 + DENDRO_799) + DENDRO_190*(DENDRO_230*DENDRO_71 + DENDRO_759 + DENDRO_785) + DENDRO_200*(DENDRO_230*DENDRO_69 + DENDRO_505) - DENDRO_302*(DENDRO_406 + DENDRO_732) - DENDRO_302*(DENDRO_240*DENDRO_75 + DENDRO_670) - DENDRO_306*(0.5*DENDRO_491 + DENDRO_737) - DENDRO_306*(1.0*DENDRO_492 + DENDRO_738) - DENDRO_306*(DENDRO_123*DENDRO_230 + DENDRO_487 + DENDRO_735) - DENDRO_315*(DENDRO_518 + DENDRO_741) - DENDRO_315*(0.5*DENDRO_519 + DENDRO_740) - DENDRO_315*(DENDRO_123*DENDRO_80 + DENDRO_333 + DENDRO_783) - DENDRO_730*(DENDRO_142 + DENDRO_529 + DENDRO_531) - DENDRO_730*(DENDRO_526 + DENDRO_527 + DENDRO_694) - DENDRO_761*(DENDRO_435 + DENDRO_557 + DENDRO_650) + DENDRO_782)) - DENDRO_31*(DENDRO_716 + alpha[pp]*(DENDRO_146*(DENDRO_676 + DENDRO_724) + DENDRO_146*(DENDRO_678 + DENDRO_679) + DENDRO_146*(DENDRO_680 - DENDRO_722) + DENDRO_146*(-DENDRO_218*DENDRO_230 + DENDRO_723) + DENDRO_146*(DENDRO_226*DENDRO_265 + DENDRO_675 + DENDRO_724) - DENDRO_184*(DENDRO_643 + DENDRO_668) - DENDRO_184*(DENDRO_672 + DENDRO_721) - DENDRO_184*(DENDRO_674 + DENDRO_721) + DENDRO_190*(DENDRO_465 + DENDRO_686) + DENDRO_190*(DENDRO_469 + DENDRO_690) + DENDRO_190*(DENDRO_469 + DENDRO_727) + DENDRO_190*(DENDRO_538 + DENDRO_726) + DENDRO_190*(DENDRO_539 + DENDRO_729) + DENDRO_190*(DENDRO_693 + DENDRO_728) - DENDRO_302*(DENDRO_396 + DENDRO_661) - DENDRO_302*(DENDRO_660 + DENDRO_717) - DENDRO_306*(DENDRO_663 + DENDRO_718) - DENDRO_306*(DENDRO_665 + DENDRO_719) - DENDRO_306*(DENDRO_226*DENDRO_230 + DENDRO_480 + DENDRO_662) - DENDRO_315*(DENDRO_517 + DENDRO_667) - DENDRO_315*(DENDRO_147*DENDRO_155 + DENDRO_720) + DENDRO_410*(DENDRO_204*DENDRO_230 + DENDRO_504) + DENDRO_703 - DENDRO_730*(DENDRO_637 + DENDRO_669) + DENDRO_731)) - DENDRO_31*(DENDRO_716 + alpha[pp]*(DENDRO_146*(DENDRO_679 + DENDRO_680) + DENDRO_146*(DENDRO_681 + DENDRO_683) + DENDRO_146*(-DENDRO_127*DENDRO_247 + DENDRO_678) + DENDRO_146*(DENDRO_214*DENDRO_230 + DENDRO_683) + DENDRO_146*(DENDRO_227*DENDRO_265 + DENDRO_676) - DENDRO_184*(DENDRO_403 + DENDRO_668) - DENDRO_184*(DENDRO_438 + DENDRO_560) - DENDRO_184*(DENDRO_438 + DENDRO_671) - DENDRO_184*(DENDRO_565 + DENDRO_672) - DENDRO_184*(DENDRO_570 + DENDRO_674) - DENDRO_184*(DENDRO_571 + DENDRO_670) + DENDRO_189*(DENDRO_585 + DENDRO_656) + DENDRO_189*(DENDRO_649 + DENDRO_657) + DENDRO_190*(DENDRO_539 + DENDRO_692) + DENDRO_190*(DENDRO_684 + DENDRO_686) + DENDRO_190*(DENDRO_687 + DENDRO_690) + DENDRO_190*(DENDRO_692 + DENDRO_693) - DENDRO_302*(0.5*DENDRO_409 + DENDRO_661) - DENDRO_302*(1.0*DENDRO_440 + DENDRO_658) - DENDRO_302*(DENDRO_214*DENDRO_265 + DENDRO_399 + DENDRO_660) - DENDRO_306*(DENDRO_474 + DENDRO_665) - DENDRO_306*(0.5*DENDRO_479 + DENDRO_662) - DENDRO_306*(DENDRO_214*DENDRO_268 + DENDRO_477 + DENDRO_663) - DENDRO_315*(DENDRO_516 + DENDRO_667) - DENDRO_315*(DENDRO_123*DENDRO_134 + DENDRO_652) + DENDRO_410*(DENDRO_202*DENDRO_265 + DENDRO_425) - DENDRO_695*(DENDRO_528 + DENDRO_694) + DENDRO_703)) + DENDRO_35*(DENDRO_355*(DENDRO_455 + DENDRO_64*(DENDRO_449 + DENDRO_453)) + DENDRO_362*(DENDRO_250 + DENDRO_359*DENDRO_446) - 4*DENDRO_445 + DENDRO_447*(DENDRO_268 + DENDRO_446*DENDRO_88) + alpha[pp]*(-DENDRO_121*(DENDRO_117*DENDRO_459 + DENDRO_118*(DENDRO_228 + DENDRO_460 + DENDRO_461) + DENDRO_120*DENDRO_462 + DENDRO_456) + DENDRO_122*DENDRO_416 + DENDRO_146*(DENDRO_474 + DENDRO_475) + DENDRO_146*(DENDRO_477 - DENDRO_478) + DENDRO_146*(DENDRO_480 + 1.0*DENDRO_481) - DENDRO_155*DENDRO_498 - DENDRO_184*(DENDRO_463 + DENDRO_465) - DENDRO_184*(-1.0*DENDRO_466 + DENDRO_468) - DENDRO_184*(DENDRO_242*DENDRO_470 + DENDRO_471) - DENDRO_184*(DENDRO_259*DENDRO_470 + DENDRO_469) + DENDRO_190*(DENDRO_482 + DENDRO_483) + DENDRO_190*(DENDRO_484 + DENDRO_485) + DENDRO_190*(DENDRO_487 + 1.0*DENDRO_488) + DENDRO_200*(DENDRO_486 + DENDRO_488) + DENDRO_200*(DENDRO_122*DENDRO_242 + DENDRO_492) + DENDRO_200*(DENDRO_155*DENDRO_202 + DENDRO_491) + DENDRO_202*DENDRO_414 - DENDRO_221*DENDRO_231*DENDRO_428 + DENDRO_221*DENDRO_412 - DENDRO_240*DENDRO_495 - DENDRO_242*DENDRO_497 + DENDRO_410*(DENDRO_479 + DENDRO_481) + DENDRO_410*(DENDRO_122*DENDRO_240 + DENDRO_489) + DENDRO_410*(DENDRO_202*DENDRO_259 + DENDRO_490) + DENDRO_419*DENDRO_451 - DENDRO_493*(DENDRO_216 - DENDRO_404) - DENDRO_496*(DENDRO_123 + DENDRO_125) + DENDRO_506)) + DENDRO_46*(DENDRO_635 + alpha[pp]*(DENDRO_145*(DENDRO_402 + DENDRO_637) + DENDRO_146*(DENDRO_556 + DENDRO_642) + DENDRO_146*(DENDRO_563 + DENDRO_643) + DENDRO_146*(DENDRO_564 + DENDRO_644) + DENDRO_146*(DENDRO_567 + DENDRO_644) - DENDRO_184*(DENDRO_548 + DENDRO_640) - DENDRO_184*(-DENDRO_119*DENDRO_151 + DENDRO_553) - DENDRO_184*(DENDRO_265*DENDRO_75 + DENDRO_640) + DENDRO_190*(DENDRO_517 + DENDRO_582) + DENDRO_190*(DENDRO_577 + DENDRO_645) + DENDRO_190*(DENDRO_581 + DENDRO_647) + DENDRO_190*(DENDRO_645 + DENDRO_646) - DENDRO_197*(DENDRO_119*DENDRO_132 + DENDRO_524) - DENDRO_302*(DENDRO_535 + DENDRO_638) - DENDRO_302*(DENDRO_247*DENDRO_75 + DENDRO_390 + DENDRO_533) - DENDRO_306*(DENDRO_155*DENDRO_214 + DENDRO_639) - DENDRO_315*(1.0*DENDRO_523 + DENDRO_543) - DENDRO_315*(DENDRO_119*DENDRO_75 + DENDRO_513 + DENDRO_542) - DENDRO_586*(DENDRO_464 + DENDRO_649) + DENDRO_616 + DENDRO_655)) + DENDRO_46*(DENDRO_635 + alpha[pp]*(DENDRO_146*(DENDRO_407 + DENDRO_563) + DENDRO_146*(DENDRO_407 + DENDRO_571) + DENDRO_146*(DENDRO_436 + DENDRO_556) + DENDRO_146*(DENDRO_558 + DENDRO_560) + DENDRO_146*(DENDRO_564 + DENDRO_565) + DENDRO_146*(DENDRO_567 + DENDRO_570) - DENDRO_184*(DENDRO_549 + DENDRO_550) - DENDRO_184*(DENDRO_551 + DENDRO_553) - DENDRO_184*(DENDRO_219*DENDRO_71 + DENDRO_546) - DENDRO_184*(DENDRO_119*DENDRO_147 + DENDRO_551 + DENDRO_552) - DENDRO_184*(DENDRO_265*DENDRO_76 + DENDRO_547 + DENDRO_548) + DENDRO_189*(DENDRO_526 + DENDRO_528) + DENDRO_189*(DENDRO_529 + DENDRO_532) + DENDRO_190*(DENDRO_516 + DENDRO_581) + DENDRO_190*(DENDRO_572 + DENDRO_573) + DENDRO_190*(DENDRO_579 + DENDRO_582) + DENDRO_190*(DENDRO_574 + DENDRO_575 + DENDRO_577) - DENDRO_197*(DENDRO_265*DENDRO_74 + DENDRO_426) - DENDRO_302*(1.0*DENDRO_417 + DENDRO_533) - DENDRO_302*(0.5*DENDRO_442 + DENDRO_536) - DENDRO_302*(DENDRO_147*DENDRO_265 + DENDRO_393 + DENDRO_535) - DENDRO_306*(DENDRO_503 + DENDRO_541) - DENDRO_306*(DENDRO_240*DENDRO_70 + DENDRO_538) - DENDRO_315*(DENDRO_328 + DENDRO_544) - DENDRO_315*(0.5*DENDRO_512 + DENDRO_542) - DENDRO_315*(DENDRO_116*DENDRO_147 + DENDRO_511 + DENDRO_543) - DENDRO_586*(DENDRO_583 + DENDRO_585) + DENDRO_616)) + DENDRO_49*(DENDRO_356*DENDRO_81 + DENDRO_361*(-DENDRO_100 + DENDRO_102 + DENDRO_111 - DENDRO_98) + DENDRO_447*(DENDRO_116 + DENDRO_89) - 4*DENDRO_58 + alpha[pp]*(-DENDRO_121*(DENDRO_115 + DENDRO_117*DENDRO_507 + DENDRO_118*DENDRO_508 + DENDRO_120*(DENDRO_101 + DENDRO_509 + DENDRO_510)) - DENDRO_144*DENDRO_252*DENDRO_320 + DENDRO_146*(-1.0*DENDRO_135 + DENDRO_143) + DENDRO_146*(-1.0*DENDRO_156 + DENDRO_162) + DENDRO_146*(DENDRO_160*DENDRO_332 + DENDRO_515) + DENDRO_146*(DENDRO_235*DENDRO_70 + DENDRO_517) + DENDRO_146*(DENDRO_242*DENDRO_75 + DENDRO_516) - DENDRO_155*DENDRO_312 - DENDRO_184*(DENDRO_513 + 1.0*DENDRO_514) - DENDRO_184*(-DENDRO_116*DENDRO_182 + DENDRO_511) + DENDRO_190*(DENDRO_520 + 1.0*DENDRO_521) + DENDRO_190*(DENDRO_116*DENDRO_167 + DENDRO_518) - DENDRO_197*(DENDRO_512 + DENDRO_514) - DENDRO_197*(DENDRO_160*DENDRO_74 + DENDRO_523) + DENDRO_200*(DENDRO_519 + DENDRO_521) + DENDRO_200*(DENDRO_155*DENDRO_74 + DENDRO_522) + DENDRO_234*DENDRO_411 + DENDRO_254*DENDRO_415 - DENDRO_270*DENDRO_313*DENDRO_349 + DENDRO_272*DENDRO_413 + DENDRO_300*DENDRO_418 - DENDRO_303*(DENDRO_149 - DENDRO_165) - DENDRO_307*(DENDRO_125 - DENDRO_304) - DENDRO_317*DENDRO_524 - DENDRO_319*DENDRO_525 + DENDRO_350)) + DENDRO_55*(-4*DENDRO_352 + DENDRO_354*DENDRO_356 + DENDRO_362*(DENDRO_247 + DENDRO_360) + DENDRO_373*(-DENDRO_363 + DENDRO_364 - DENDRO_365 + DENDRO_372) + alpha[pp]*(-DENDRO_121*(DENDRO_117*(DENDRO_375 + DENDRO_378 + DENDRO_381) + DENDRO_118*DENDRO_383 + DENDRO_120*DENDRO_384 + DENDRO_374) + DENDRO_132*DENDRO_416 + DENDRO_146*(DENDRO_396 + DENDRO_397) + DENDRO_146*(DENDRO_399 + 1.0*DENDRO_400) - DENDRO_184*(DENDRO_387 + DENDRO_388) - DENDRO_184*(DENDRO_393 + 1.0*DENDRO_394) - DENDRO_184*(DENDRO_247*DENDRO_391 + DENDRO_390) + DENDRO_190*(DENDRO_401 + DENDRO_403) + DENDRO_190*(DENDRO_147*DENDRO_259 + DENDRO_405) + DENDRO_190*(DENDRO_160*DENDRO_214 + DENDRO_406) + DENDRO_190*(DENDRO_235*DENDRO_408 + DENDRO_407) - DENDRO_197*(DENDRO_392 + DENDRO_394) - DENDRO_197*(DENDRO_132*DENDRO_235 + DENDRO_417) + DENDRO_204*DENDRO_412 - DENDRO_211*DENDRO_266*DENDRO_428 + DENDRO_211*DENDRO_414 - DENDRO_235*DENDRO_424 - DENDRO_240*DENDRO_423 - DENDRO_319*DENDRO_425 + DENDRO_367*DENDRO_419 + DENDRO_410*(DENDRO_398 + DENDRO_400) + DENDRO_410*(DENDRO_132*DENDRO_240 + DENDRO_409) - DENDRO_420*(DENDRO_147 + DENDRO_149) - DENDRO_421*(DENDRO_214 + DENDRO_216) - DENDRO_426*DENDRO_427 + DENDRO_444));
double DENDRO_801 = DENDRO_41*DENDRO_800;
double DENDRO_802 = (1.0/12.0)*chi[pp];
double DENDRO_803 = (1.0/3.0)*At1[pp];
double DENDRO_804 = At4[pp]*DENDRO_31;
double DENDRO_805 = At3[pp]*DENDRO_35;
double DENDRO_806 = DENDRO_43 + DENDRO_804 - DENDRO_805;
double DENDRO_807 = At3[pp]*DENDRO_29;
double DENDRO_808 = At4[pp]*DENDRO_46;
double DENDRO_809 = -At1[pp]*DENDRO_49 + DENDRO_807 - DENDRO_808;
double DENDRO_810 = At4[pp]*DENDRO_55;
double DENDRO_811 = -At1[pp]*DENDRO_46 + At3[pp]*DENDRO_31 - DENDRO_810;
double DENDRO_812 = 6.0*DENDRO_94;
double DENDRO_813 = 6.0*DENDRO_112;
double DENDRO_814 = 6.0*DENDRO_82;
double DENDRO_815 = DENDRO_131*DENDRO_179;
double DENDRO_816 = DENDRO_167*DENDRO_179;
double DENDRO_817 = -DENDRO_237 + DENDRO_238 + DENDRO_239;
double DENDRO_818 = DENDRO_74*DENDRO_817;
double DENDRO_819 = DENDRO_816 + DENDRO_818;
double DENDRO_820 = DENDRO_173*DENDRO_472;
double DENDRO_821 = -DENDRO_164*DENDRO_476;
double DENDRO_822 = DENDRO_223 - DENDRO_228 + DENDRO_229;
double DENDRO_823 = DENDRO_168*DENDRO_688;
double DENDRO_824 = -DENDRO_256 + DENDRO_257 + DENDRO_258;
double DENDRO_825 = DENDRO_214*DENDRO_93;
double DENDRO_826 = DENDRO_174*DENDRO_824 + DENDRO_825;
double DENDRO_827 = DENDRO_179*DENDRO_327;
double DENDRO_828 = DENDRO_187*DENDRO_332 + 0.25*DENDRO_77*DENDRO_817;
double DENDRO_829 = DENDRO_175 + DENDRO_181;
double DENDRO_830 = DENDRO_167*DENDRO_824;
double DENDRO_831 = 1.0*DENDRO_145;
double DENDRO_832 = DENDRO_167*DENDRO_173;
double DENDRO_833 = DENDRO_69*DENDRO_817 + DENDRO_832;
double DENDRO_834 = 2.0*DENDRO_233;
double DENDRO_835 = 2.0*DENDRO_253;
double DENDRO_836 = 2.0*DENDRO_271;
double DENDRO_837 = DENDRO_299*DENDRO_64;
double DENDRO_838 = (1.0/3.0)*At2[pp];
double DENDRO_839 = At5[pp]*DENDRO_31;
double DENDRO_840 = At2[pp]*DENDRO_29 - At4[pp]*DENDRO_35 + DENDRO_839;
double DENDRO_841 = At5[pp]*DENDRO_46;
double DENDRO_842 = -At2[pp]*DENDRO_49 + At4[pp]*DENDRO_29 - DENDRO_841;
double DENDRO_843 = -At5[pp]*DENDRO_55 + DENDRO_48 + DENDRO_804;
double DENDRO_844 = DENDRO_137*DENDRO_173;
double DENDRO_845 = DENDRO_179*DENDRO_389;
double DENDRO_846 = DENDRO_244 + DENDRO_245 - DENDRO_246;
double DENDRO_847 = 0.25*DENDRO_164*DENDRO_211;
double DENDRO_848 = DENDRO_168*DENDRO_404;
double DENDRO_849 = DENDRO_262 - DENDRO_263 + DENDRO_264;
double DENDRO_850 = DENDRO_165*DENDRO_168;
double DENDRO_851 = DENDRO_212*DENDRO_93;
double DENDRO_852 = DENDRO_173*DENDRO_191;
double DENDRO_853 = -DENDRO_847;
double DENDRO_854 = DENDRO_173*DENDRO_389;
double DENDRO_855 = DENDRO_132*DENDRO_173;
double DENDRO_856 = 2*At4[pp];
double DENDRO_857 = At3[pp]*DENDRO_51;
double DENDRO_858 = DENDRO_41*DENDRO_856;
double DENDRO_859 = DENDRO_113*DENDRO_41;
double DENDRO_860 = DENDRO_202*DENDRO_822;
double DENDRO_861 = 1.0*DENDRO_817;
double DENDRO_862 = 0.25*DENDRO_832;
double DENDRO_863 = DENDRO_164*DENDRO_404;
double DENDRO_864 = DENDRO_122*DENDRO_822;
double DENDRO_865 = (1.0/3.0)*At4[pp];
double DENDRO_866 = DENDRO_165*DENDRO_824;
double DENDRO_867 = -DENDRO_404*DENDRO_824 + DENDRO_675;
double DENDRO_868 = DENDRO_689 - DENDRO_863;
double DENDRO_869 = DENDRO_211*DENDRO_824;
double DENDRO_870 = DENDRO_204*DENDRO_849;
double DENDRO_871 = DENDRO_168*DENDRO_211;
double DENDRO_872 = DENDRO_132*DENDRO_849;
double DENDRO_873 = DENDRO_41*DENDRO_82;
double DENDRO_874 = DENDRO_109*DENDRO_64;
double DENDRO_875 = DENDRO_112*DENDRO_41;
double DENDRO_876 = DENDRO_370*DENDRO_64;
double DENDRO_877 = DENDRO_41*DENDRO_94;
double DENDRO_878 = 0.5*DENDRO_94;
double DENDRO_879 = 0.5*DENDRO_82;
double DENDRO_880 = DENDRO_284*DENDRO_41;
double DENDRO_881 = 0.5*DENDRO_112;
double DENDRO_882 = DENDRO_287*DENDRO_41;
double DENDRO_883 = DENDRO_291*DENDRO_41;
double DENDRO_884 = pow(DENDRO_29, 2);
double DENDRO_885 = pow(DENDRO_46, 2);
double DENDRO_886 = 2*DENDRO_49;
double DENDRO_887 = At0[pp]*pow(DENDRO_49, 2) + At3[pp]*DENDRO_884 + At5[pp]*DENDRO_885 - DENDRO_291*DENDRO_808 - DENDRO_43*DENDRO_886 + DENDRO_47*DENDRO_886;
double DENDRO_888 = 3*DENDRO_144;
double DENDRO_889 = pow(DENDRO_31, 2);
double DENDRO_890 = 2*DENDRO_35;
double DENDRO_891 = At0[pp]*DENDRO_884 + At3[pp]*pow(DENDRO_35, 2) + At5[pp]*DENDRO_889 + DENDRO_291*DENDRO_32 - DENDRO_43*DENDRO_890 - DENDRO_804*DENDRO_890;
double DENDRO_892 = 2*DENDRO_55;
double DENDRO_893 = At0[pp]*DENDRO_885 + At3[pp]*DENDRO_889 + At5[pp]*pow(DENDRO_55, 2) - DENDRO_287*DENDRO_53 + DENDRO_47*DENDRO_892 - DENDRO_804*DENDRO_892;
double DENDRO_894 = At2[pp]*DENDRO_885 - DENDRO_29*DENDRO_810 + DENDRO_31*DENDRO_807 - DENDRO_31*DENDRO_808 - DENDRO_43*DENDRO_46 + DENDRO_46*DENDRO_50 - DENDRO_49*DENDRO_53 + DENDRO_49*DENDRO_56 + DENDRO_55*DENDRO_841;
double DENDRO_895 = 6*DENDRO_144;
double DENDRO_896 = At1[pp]*DENDRO_884;
double DENDRO_897 = DENDRO_29*DENDRO_804;
double DENDRO_898 = DENDRO_29*DENDRO_47;
double DENDRO_899 = DENDRO_35*DENDRO_808;
double DENDRO_900 = DENDRO_31*DENDRO_841;
double DENDRO_901 = DENDRO_29*DENDRO_50;
double DENDRO_902 = DENDRO_36*DENDRO_49;
double DENDRO_903 = DENDRO_32*DENDRO_49;
double DENDRO_904 = DENDRO_35*DENDRO_807;
double DENDRO_905 = DENDRO_896 + DENDRO_897 - DENDRO_898 + DENDRO_899 - DENDRO_900 - DENDRO_901 + DENDRO_902 - DENDRO_903 - DENDRO_904;
double DENDRO_906 = At4[pp]*DENDRO_889;
double DENDRO_907 = DENDRO_31*DENDRO_43;
double DENDRO_908 = DENDRO_30*DENDRO_46;
double DENDRO_909 = DENDRO_36*DENDRO_46;
double DENDRO_910 = DENDRO_31*DENDRO_47;
double DENDRO_911 = DENDRO_29*DENDRO_56;
double DENDRO_912 = DENDRO_31*DENDRO_805;
double DENDRO_913 = DENDRO_35*DENDRO_810;
double DENDRO_914 = DENDRO_55*DENDRO_839;
double DENDRO_915 = DENDRO_906 + DENDRO_907 - DENDRO_908 + DENDRO_909 - DENDRO_910 - DENDRO_911 - DENDRO_912 + DENDRO_913 - DENDRO_914;
double DENDRO_916 = (1.0/3.0)*alpha[pp];
double DENDRO_917 = 2*DENDRO_144;
double DENDRO_918 = DENDRO_112*DENDRO_917;
double DENDRO_919 = DENDRO_82*DENDRO_917;
double DENDRO_920 = DENDRO_917*DENDRO_94;
double DENDRO_921 = grad2_0_2_beta0[pp];
double DENDRO_922 = (7.0/3.0)*DENDRO_329;
double DENDRO_923 = grad2_0_0_beta0[pp];
double DENDRO_924 = 4*grad_0_K[pp];
double DENDRO_925 = 9*DENDRO_64;
double DENDRO_926 = DENDRO_120*DENDRO_925;
double DENDRO_927 = DENDRO_41*DENDRO_916;
double DENDRO_928 = 4*grad_2_K[pp];
double DENDRO_929 = DENDRO_117*DENDRO_925;
double DENDRO_930 = 4*grad_1_K[pp];
double DENDRO_931 = DENDRO_118*DENDRO_925;
double DENDRO_932 = grad2_0_1_beta1[pp];
double DENDRO_933 = (1.0/3.0)*DENDRO_338;
double DENDRO_934 = grad2_0_2_beta2[pp];
double DENDRO_935 = grad2_1_2_beta1[pp];
double DENDRO_936 = (1.0/3.0)*DENDRO_329;
double DENDRO_937 = grad2_2_2_beta2[pp];
double DENDRO_938 = grad2_1_1_beta1[pp];
double DENDRO_939 = DENDRO_29*DENDRO_41;
double DENDRO_940 = (1.0/3.0)*DENDRO_939;
double DENDRO_941 = grad2_1_2_beta2[pp];
double DENDRO_942 = (2.0/3.0)*DENDRO_25;
double DENDRO_943 = grad2_0_1_beta0[pp];
double DENDRO_944 = (7.0/3.0)*DENDRO_939;
double DENDRO_945 = pow(DENDRO_40, -3);
double DENDRO_946 = DENDRO_0*DENDRO_945;
double DENDRO_947 = DENDRO_891*DENDRO_946;
double DENDRO_948 = DENDRO_893*DENDRO_946;
double DENDRO_949 = DENDRO_887*DENDRO_946;
double DENDRO_950 = 2.0*DENDRO_945*alpha[pp];
double DENDRO_951 = DENDRO_894*DENDRO_950;
double DENDRO_952 = DENDRO_905*DENDRO_950;
double DENDRO_953 = DENDRO_915*DENDRO_950;
double DENDRO_954 = beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp];
double DENDRO_955 = -DENDRO_11*DENDRO_411 + DENDRO_119*DENDRO_949 - DENDRO_15*DENDRO_413 - DENDRO_2*DENDRO_415 + DENDRO_235*DENDRO_951 + DENDRO_240*DENDRO_953 + DENDRO_242*DENDRO_952 + DENDRO_247*DENDRO_948 + DENDRO_250*DENDRO_947 - DENDRO_334*grad2_2_2_beta0[pp] - DENDRO_336*grad2_1_1_beta0[pp] - 4.0/3.0*DENDRO_338*DENDRO_923 + DENDRO_415*DENDRO_942 + DENDRO_880*grad2_1_2_beta0[pp] - DENDRO_887*DENDRO_918 - DENDRO_894*DENDRO_920 - DENDRO_905*DENDRO_919 - DENDRO_921*DENDRO_922 - DENDRO_927*(DENDRO_29*DENDRO_930 + DENDRO_905*DENDRO_931) - DENDRO_927*(-DENDRO_46*DENDRO_928 + DENDRO_894*DENDRO_929) - DENDRO_927*(-DENDRO_49*DENDRO_924 + DENDRO_887*DENDRO_926) - DENDRO_932*DENDRO_933 - DENDRO_933*DENDRO_934 - DENDRO_935*DENDRO_936 - DENDRO_936*DENDRO_937 + DENDRO_938*DENDRO_940 + DENDRO_940*DENDRO_941 + DENDRO_943*DENDRO_944 + DENDRO_954;
double DENDRO_956 = -DENDRO_896 - DENDRO_897 + DENDRO_898 - DENDRO_899 + DENDRO_900 + DENDRO_901 - DENDRO_902 + DENDRO_903 + DENDRO_904;
double DENDRO_957 = -DENDRO_906 - DENDRO_907 + DENDRO_908 - DENDRO_909 + DENDRO_910 + DENDRO_911 + DENDRO_912 - DENDRO_913 + DENDRO_914;
double DENDRO_958 = DENDRO_31*DENDRO_928;
double DENDRO_959 = DENDRO_29*DENDRO_924;
double DENDRO_960 = DENDRO_35*DENDRO_930;
double DENDRO_961 = DENDRO_891*DENDRO_931;
double DENDRO_962 = (1.0/3.0)*DENDRO_336;
double DENDRO_963 = DENDRO_31*DENDRO_41;
double DENDRO_964 = (1.0/3.0)*DENDRO_921;
double DENDRO_965 = (1.0/3.0)*DENDRO_963;
double DENDRO_966 = (7.0/3.0)*DENDRO_963;
double DENDRO_967 = beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp];
double DENDRO_968 = DENDRO_134*DENDRO_951 + DENDRO_219*DENDRO_948 - DENDRO_334*grad2_2_2_beta1[pp] - 4.0/3.0*DENDRO_336*DENDRO_938 - DENDRO_338*grad2_0_0_beta1[pp] + DENDRO_80*DENDRO_949 - DENDRO_882*grad2_0_2_beta1[pp] - DENDRO_891*DENDRO_919 + DENDRO_923*DENDRO_940 + DENDRO_932*DENDRO_944 + DENDRO_934*DENDRO_940 + DENDRO_935*DENDRO_966 + DENDRO_937*DENDRO_965 - DENDRO_941*DENDRO_962 - DENDRO_943*DENDRO_962 + DENDRO_963*DENDRO_964 + DENDRO_967;
double DENDRO_969 = beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp];
double DENDRO_970 = DENDRO_116*DENDRO_949 - DENDRO_12*DENDRO_411 + DENDRO_155*DENDRO_952 + DENDRO_160*DENDRO_951 + DENDRO_259*DENDRO_953 + DENDRO_265*DENDRO_948 + DENDRO_268*DENDRO_947 - 1.0/3.0*DENDRO_334*DENDRO_935 - 4.0/3.0*DENDRO_334*DENDRO_937 - DENDRO_334*DENDRO_964 - DENDRO_336*grad2_1_1_beta2[pp] - DENDRO_338*grad2_0_0_beta2[pp] - DENDRO_413*DENDRO_6 + DENDRO_413*DENDRO_942 - DENDRO_415*DENDRO_9 + DENDRO_883*grad2_0_1_beta2[pp] - DENDRO_893*DENDRO_920 - DENDRO_894*DENDRO_918 - DENDRO_915*DENDRO_919 - DENDRO_922*DENDRO_934 - DENDRO_923*DENDRO_936 - DENDRO_927*(DENDRO_31*DENDRO_930 + DENDRO_915*DENDRO_931) - DENDRO_927*(-DENDRO_46*DENDRO_924 + DENDRO_894*DENDRO_926) - DENDRO_927*(-DENDRO_55*DENDRO_928 + DENDRO_893*DENDRO_929) - DENDRO_932*DENDRO_936 + DENDRO_938*DENDRO_965 + DENDRO_941*DENDRO_966 + DENDRO_943*DENDRO_965 + DENDRO_969;

// Dendro: printing variables
//--
a_rhs[pp] = -DENDRO_0*K[pp] + lambda[0]*(beta0[pp]*agrad_0_alpha[pp] + beta1[pp]*agrad_1_alpha[pp] + beta2[pp]*agrad_2_alpha[pp]);
//--
b_rhs0[pp] = B0[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta0[pp] + beta1[pp]*agrad_1_beta0[pp] + beta2[pp]*agrad_2_beta0[pp]);
//--
b_rhs1[pp] = B1[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta1[pp] + beta1[pp]*agrad_1_beta1[pp] + beta2[pp]*agrad_2_beta1[pp]);
//--
b_rhs2[pp] = B2[pp]*DENDRO_1 + lambda[1]*(beta0[pp]*agrad_0_beta2[pp] + beta1[pp]*agrad_1_beta2[pp] + beta2[pp]*agrad_2_beta2[pp]);
//--
gt_rhs00[pp] = -At0[pp]*DENDRO_0 + DENDRO_10*DENDRO_9 + DENDRO_3*gt0[pp] - DENDRO_4*DENDRO_5 - DENDRO_5*DENDRO_6 + DENDRO_7*DENDRO_8 + beta0[pp]*agrad_0_gt0[pp] + beta1[pp]*agrad_1_gt0[pp] + beta2[pp]*agrad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp]*DENDRO_0 + DENDRO_11*gt0[pp] + DENDRO_12*gt2[pp] + DENDRO_13*DENDRO_2 + DENDRO_13*DENDRO_4 - DENDRO_14*gt1[pp] + DENDRO_7*gt3[pp] + DENDRO_9*gt4[pp] + beta0[pp]*agrad_0_gt1[pp] + beta1[pp]*agrad_1_gt1[pp] + beta2[pp]*agrad_2_gt1[pp];
//--
gt_rhs02[pp] = -At2[pp]*DENDRO_0 + DENDRO_15*gt0[pp] + DENDRO_16*gt1[pp] + DENDRO_17*DENDRO_2 + DENDRO_17*DENDRO_6 - DENDRO_18*gt2[pp] + DENDRO_7*gt4[pp] + DENDRO_9*gt5[pp] + beta0[pp]*agrad_0_gt2[pp] + beta1[pp]*agrad_1_gt2[pp] + beta2[pp]*agrad_2_gt2[pp];
//--
gt_rhs11[pp] = -At3[pp]*DENDRO_0 + DENDRO_11*DENDRO_8 + DENDRO_12*DENDRO_21 - DENDRO_14*gt3[pp] - DENDRO_19*gt3[pp] + DENDRO_20*gt3[pp] + beta0[pp]*agrad_0_gt3[pp] + beta1[pp]*agrad_1_gt3[pp] + beta2[pp]*agrad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp]*DENDRO_0 + DENDRO_11*gt2[pp] + DENDRO_12*gt5[pp] + DENDRO_15*gt1[pp] + DENDRO_16*gt3[pp] - DENDRO_19*gt4[pp] + DENDRO_22*DENDRO_4 + DENDRO_22*DENDRO_6 + beta0[pp]*agrad_0_gt4[pp] + beta1[pp]*agrad_1_gt4[pp] + beta2[pp]*agrad_2_gt4[pp];
//--
gt_rhs22[pp] = -At5[pp]*DENDRO_0 + DENDRO_10*DENDRO_15 + DENDRO_16*DENDRO_21 - DENDRO_18*gt5[pp] - DENDRO_19*gt5[pp] + DENDRO_23*gt5[pp] + beta0[pp]*agrad_0_gt5[pp] + beta1[pp]*agrad_1_gt5[pp] + beta2[pp]*agrad_2_gt5[pp];
//--
chi_rhs[pp] = -DENDRO_24*DENDRO_25 + DENDRO_24*K[pp]*alpha[pp] + beta0[pp]*agrad_0_chi[pp] + beta1[pp]*agrad_1_chi[pp] + beta2[pp]*agrad_2_chi[pp];
//--
At_rhs00[pp] = -At0[pp]*DENDRO_14 - At0[pp]*DENDRO_18 + At0[pp]*DENDRO_3 + DENDRO_26*DENDRO_7 + DENDRO_27*DENDRO_9 + DENDRO_802*(-DENDRO_113*(DENDRO_100 - DENDRO_102 - DENDRO_111 + DENDRO_98) + DENDRO_351*(-DENDRO_121*(DENDRO_115 + DENDRO_116*DENDRO_117 + DENDRO_118*DENDRO_80 + DENDRO_119*DENDRO_120) - DENDRO_146*(DENDRO_135 - DENDRO_143) - DENDRO_146*(DENDRO_156 + DENDRO_163) - DENDRO_146*(DENDRO_166 + DENDRO_167*DENDRO_169) - DENDRO_146*(DENDRO_175 + DENDRO_179*DENDRO_70) - DENDRO_146*(DENDRO_173*DENDRO_75 + DENDRO_181) + DENDRO_164*DENDRO_312 + DENDRO_184*(DENDRO_193 + 1.0*DENDRO_194) - DENDRO_184*(-DENDRO_169*DENDRO_76 + DENDRO_182*DENDRO_93) + DENDRO_187*DENDRO_314*DENDRO_320 - DENDRO_190*(DENDRO_186 + 1.0*DENDRO_188) - DENDRO_190*(DENDRO_164*DENDRO_191 + 1.0*DENDRO_167*DENDRO_93) + DENDRO_197*(DENDRO_192 + DENDRO_194) + DENDRO_197*(DENDRO_195 + DENDRO_196) - DENDRO_200*(DENDRO_185 + DENDRO_188) - DENDRO_200*(DENDRO_198 + DENDRO_199) - DENDRO_233*DENDRO_234 - DENDRO_253*DENDRO_254 - DENDRO_271*DENDRO_272 - DENDRO_299*DENDRO_300 + DENDRO_303*(DENDRO_150 + DENDRO_165) + DENDRO_307*(DENDRO_126 + DENDRO_304) + DENDRO_313*DENDRO_315*DENDRO_93 + DENDRO_316*DENDRO_317 + DENDRO_318*DENDRO_319 + DENDRO_350) - 12*DENDRO_58 + DENDRO_801*gt0[pp] + DENDRO_81*DENDRO_84 - DENDRO_96*(-DENDRO_89 + DENDRO_93)) - alpha[pp]*(-At0[pp]*K[pp] + DENDRO_42*(DENDRO_30 + DENDRO_32 - DENDRO_36) + DENDRO_52*(DENDRO_43 + DENDRO_48 - DENDRO_50) + DENDRO_57*(-At0[pp]*DENDRO_46 + DENDRO_53 - DENDRO_56)) + beta0[pp]*agrad_0_At0[pp] + beta1[pp]*agrad_1_At0[pp] + beta2[pp]*agrad_2_At0[pp];
//--
At_rhs01[pp] = At0[pp]*DENDRO_11 - At1[pp]*DENDRO_14 + At2[pp]*DENDRO_12 + At3[pp]*DENDRO_7 + At4[pp]*DENDRO_9 + DENDRO_2*DENDRO_803 + DENDRO_4*DENDRO_803 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_141*DENDRO_592 + DENDRO_155*DENDRO_590 + DENDRO_242*DENDRO_595 + DENDRO_762) - DENDRO_124*DENDRO_834 - DENDRO_129*DENDRO_836 + DENDRO_146*(-DENDRO_332*DENDRO_822 + DENDRO_752) - DENDRO_146*(-DENDRO_687 + DENDRO_747 + DENDRO_748) + DENDRO_146*(-DENDRO_173*DENDRO_688 + DENDRO_179*DENDRO_472 + DENDRO_725) + DENDRO_184*(DENDRO_823 + DENDRO_826) + DENDRO_184*(DENDRO_827 + DENDRO_828) + DENDRO_184*(DENDRO_187*DENDRO_470 + DENDRO_829) + DENDRO_184*(DENDRO_163 + DENDRO_742 + DENDRO_743) + DENDRO_190*(-DENDRO_70*DENDRO_822 + DENDRO_760) + DENDRO_190*(DENDRO_127*DENDRO_187 - DENDRO_173*DENDRO_327 + DENDRO_755) - DENDRO_190*(DENDRO_164*DENDRO_561 + DENDRO_164*DENDRO_576 + DENDRO_226*DENDRO_93) + DENDRO_190*(-DENDRO_164*DENDRO_688 - DENDRO_227*DENDRO_93 + DENDRO_753) - DENDRO_200*(DENDRO_122*DENDRO_187 + DENDRO_318) + DENDRO_301*(DENDRO_815 + DENDRO_819) + DENDRO_302*(-DENDRO_405 + DENDRO_566 + DENDRO_673) + DENDRO_306*(-DENDRO_733 + DENDRO_734) - DENDRO_306*(-DENDRO_164*DENDRO_437 + DENDRO_736 + DENDRO_821) - DENDRO_306*(-DENDRO_173*DENDRO_304 + DENDRO_739 + DENDRO_820) + DENDRO_315*(DENDRO_186 + DENDRO_187*DENDRO_70 + DENDRO_187*DENDRO_71) + DENDRO_315*(0.25*DENDRO_198 + 0.5*DENDRO_199 + DENDRO_470*DENDRO_93) - DENDRO_67*DENDRO_835 + DENDRO_767 + DENDRO_769 + DENDRO_770 + DENDRO_771 + DENDRO_772 + DENDRO_773 + DENDRO_774 + DENDRO_775 + DENDRO_776 + DENDRO_777 + DENDRO_778 + DENDRO_779 + DENDRO_780 - DENDRO_781*DENDRO_837 + DENDRO_786 - DENDRO_831*(DENDRO_122*DENDRO_179 + DENDRO_833) - DENDRO_831*(DENDRO_164*DENDRO_204 + DENDRO_168*DENDRO_202 + DENDRO_830)) - DENDRO_41*DENDRO_812*(DENDRO_164 - DENDRO_797) + DENDRO_781*DENDRO_800 - 12*DENDRO_787 - DENDRO_813*(-DENDRO_791 + DENDRO_792 + DENDRO_793 - DENDRO_795) + DENDRO_814*(DENDRO_64*(DENDRO_621 + DENDRO_63*DENDRO_781) + DENDRO_790)) - alpha[pp]*(-At1[pp]*K[pp] + DENDRO_42*DENDRO_806 + DENDRO_52*DENDRO_809 + DENDRO_57*DENDRO_811) + beta0[pp]*agrad_0_At1[pp] + beta1[pp]*agrad_1_At1[pp] + beta2[pp]*agrad_2_At1[pp];
//--
At_rhs02[pp] = At0[pp]*DENDRO_15 + At1[pp]*DENDRO_16 - At2[pp]*DENDRO_18 + At4[pp]*DENDRO_7 + At5[pp]*DENDRO_9 + DENDRO_2*DENDRO_838 + DENDRO_6*DENDRO_838 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_134*DENDRO_592 + DENDRO_160*DENDRO_590 + DENDRO_235*DENDRO_595 + DENDRO_587) - DENDRO_130*DENDRO_834 - DENDRO_146*(DENDRO_554 - DENDRO_555 + DENDRO_641) + DENDRO_146*(-DENDRO_179*DENDRO_561 - DENDRO_70*DENDRO_846 + DENDRO_854) - DENDRO_146*(DENDRO_332*DENDRO_849 + DENDRO_847 + DENDRO_848) + DENDRO_146*(DENDRO_389*DENDRO_824 - DENDRO_848 + DENDRO_853) - DENDRO_148*DENDRO_836 - DENDRO_184*(DENDRO_151*DENDRO_169 - DENDRO_850 - DENDRO_851) - DENDRO_184*(DENDRO_151*DENDRO_187 - DENDRO_179*DENDRO_191 - DENDRO_78*DENDRO_846) + DENDRO_184*(DENDRO_75*DENDRO_849 + DENDRO_850 + DENDRO_851) - DENDRO_190*(DENDRO_166 + DENDRO_826) - DENDRO_190*(DENDRO_828 + DENDRO_852) - DENDRO_190*(DENDRO_187*DENDRO_408 + DENDRO_829) - DENDRO_190*(DENDRO_166 + DENDRO_168*DENDRO_576 + DENDRO_825) + DENDRO_197*(DENDRO_132*DENDRO_187 + DENDRO_316) + DENDRO_302*(DENDRO_534 - DENDRO_638) - DENDRO_302*(-DENDRO_165*DENDRO_179 - DENDRO_75*DENDRO_846 + DENDRO_845) + DENDRO_305*(DENDRO_833 + DENDRO_844) + DENDRO_306*(DENDRO_164*DENDRO_214 + 0.25*DENDRO_830) + DENDRO_315*(0.25*DENDRO_195 + 1.0*DENDRO_196) + DENDRO_315*(DENDRO_187*DENDRO_75 + DENDRO_187*DENDRO_76 + DENDRO_193) + DENDRO_597 + DENDRO_599 + DENDRO_601 + DENDRO_602 + DENDRO_604 + DENDRO_606 + DENDRO_608 + DENDRO_609 + DENDRO_610 + DENDRO_611 + DENDRO_612 + DENDRO_613 + DENDRO_614 - DENDRO_615*DENDRO_837 + DENDRO_655 - DENDRO_72*DENDRO_835 - DENDRO_831*(DENDRO_819 + DENDRO_855)) + DENDRO_615*DENDRO_800 - 12*DENDRO_617 + DENDRO_634*DENDRO_814 - DENDRO_812*(DENDRO_618 + DENDRO_619 - DENDRO_620 - DENDRO_623) - DENDRO_813*(DENDRO_625 + DENDRO_626 - DENDRO_627 - DENDRO_630)) - alpha[pp]*(-At2[pp]*K[pp] + DENDRO_42*DENDRO_840 + DENDRO_52*DENDRO_842 + DENDRO_57*DENDRO_843) + beta0[pp]*agrad_0_At2[pp] + beta1[pp]*agrad_1_At2[pp] + beta2[pp]*agrad_2_At2[pp];
//--
At_rhs11[pp] = -At3[pp]*DENDRO_14 - At3[pp]*DENDRO_19 + At3[pp]*DENDRO_20 + DENDRO_11*DENDRO_26 + DENDRO_12*DENDRO_856 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_117*DENDRO_268 + DENDRO_118*DENDRO_230 + DENDRO_120*DENDRO_250 + DENDRO_456) - DENDRO_122*DENDRO_835 - DENDRO_146*(DENDRO_473 - 1.0*DENDRO_475) + DENDRO_146*(DENDRO_480 - 1.0*DENDRO_860) - DENDRO_146*(DENDRO_476*DENDRO_824 + DENDRO_478) + DENDRO_164*DENDRO_498 + DENDRO_173*DENDRO_497 + DENDRO_184*(DENDRO_466 - DENDRO_468) + DENDRO_184*(DENDRO_173*DENDRO_470 + DENDRO_180*DENDRO_817) + DENDRO_184*(DENDRO_470*DENDRO_824 + DENDRO_863) + DENDRO_184*(DENDRO_71*DENDRO_861 + DENDRO_862) + DENDRO_190*(DENDRO_483 + DENDRO_821) + DENDRO_190*(DENDRO_485 + DENDRO_820) + DENDRO_190*(DENDRO_487 - 1.0*DENDRO_864) + DENDRO_200*(DENDRO_486 - DENDRO_864) + DENDRO_200*(-DENDRO_122*DENDRO_173 + DENDRO_492) + DENDRO_200*(-DENDRO_164*DENDRO_202 + DENDRO_491) - DENDRO_202*DENDRO_836 + 6.0*DENDRO_221*DENDRO_305*DENDRO_822 - DENDRO_221*DENDRO_834 + DENDRO_410*(DENDRO_479 - DENDRO_860) + DENDRO_410*(-DENDRO_122*DENDRO_817 + DENDRO_489) + DENDRO_410*(-DENDRO_202*DENDRO_824 + DENDRO_490) - DENDRO_451*DENDRO_837 + DENDRO_493*(DENDRO_217 + DENDRO_404) + DENDRO_495*DENDRO_817 + DENDRO_496*(DENDRO_126 + DENDRO_457) + DENDRO_506) - 12*DENDRO_445 + DENDRO_801*gt3[pp] + DENDRO_83*(DENDRO_455 + DENDRO_64*(DENDRO_449 + DENDRO_452*DENDRO_63)) + DENDRO_859*(DENDRO_250 - DENDRO_446*(-DENDRO_106 + DENDRO_357 + DENDRO_358)) + DENDRO_96*(DENDRO_268 - DENDRO_446*(-DENDRO_85 + DENDRO_86 + DENDRO_87))) - alpha[pp]*(-At3[pp]*K[pp] + DENDRO_42*DENDRO_809 + DENDRO_806*DENDRO_857 + DENDRO_811*DENDRO_858) + beta0[pp]*agrad_0_At3[pp] + beta1[pp]*agrad_1_At3[pp] + beta2[pp]*agrad_2_At3[pp];
//--
At_rhs12[pp] = At1[pp]*DENDRO_15 + At2[pp]*DENDRO_11 + At3[pp]*DENDRO_16 - At4[pp]*DENDRO_19 + At5[pp]*DENDRO_12 + DENDRO_4*DENDRO_865 + DENDRO_6*DENDRO_865 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_208*DENDRO_592 + DENDRO_240*DENDRO_595 + DENDRO_259*DENDRO_590 + DENDRO_696) - DENDRO_128*DENDRO_835 + DENDRO_146*(DENDRO_218*DENDRO_822 + DENDRO_723) + DENDRO_146*(-DENDRO_226*DENDRO_849 + DENDRO_867) + DENDRO_146*(DENDRO_385*DENDRO_824 + DENDRO_867) - DENDRO_146*(DENDRO_123*DENDRO_846 + DENDRO_561*DENDRO_817 + DENDRO_722) + DENDRO_146*(-DENDRO_576*DENDRO_817 + DENDRO_677 - DENDRO_688*DENDRO_817) + DENDRO_183*(DENDRO_815 + DENDRO_818 + DENDRO_855) - DENDRO_184*(DENDRO_169*DENDRO_218 + DENDRO_853 - DENDRO_866) - DENDRO_184*(-DENDRO_191*DENDRO_817 - DENDRO_71*DENDRO_846 + DENDRO_854) + DENDRO_184*(DENDRO_470*DENDRO_849 + DENDRO_847 + DENDRO_866) + DENDRO_190*(-DENDRO_168*DENDRO_437 + DENDRO_868) + DENDRO_190*(-DENDRO_408*DENDRO_822 + DENDRO_728) + DENDRO_190*(-DENDRO_688*DENDRO_824 + DENDRO_868) + DENDRO_190*(-DENDRO_179*DENDRO_304 + DENDRO_725 - 0.25*DENDRO_844) + DENDRO_190*(-DENDRO_327*DENDRO_817 + DENDRO_685 - DENDRO_862) + DENDRO_190*(-DENDRO_470*DENDRO_822 + DENDRO_539 + DENDRO_691) - DENDRO_215*DENDRO_836 - DENDRO_224*DENDRO_834 + DENDRO_302*(DENDRO_659 - DENDRO_717) - DENDRO_302*(-DENDRO_165*DENDRO_817 + DENDRO_389*DENDRO_817 - DENDRO_470*DENDRO_846) - DENDRO_306*(-DENDRO_437*DENDRO_824 + DENDRO_718) - DENDRO_306*(-DENDRO_226*DENDRO_822 - DENDRO_227*DENDRO_822 + DENDRO_480) - DENDRO_306*(-DENDRO_304*DENDRO_817 + DENDRO_664 + DENDRO_719) + DENDRO_315*(DENDRO_147*DENDRO_164 + DENDRO_823) + DENDRO_315*(DENDRO_175 + DENDRO_827 + DENDRO_852) + DENDRO_410*(-DENDRO_204*DENDRO_822 + DENDRO_504) - DENDRO_701*DENDRO_837 + DENDRO_702 + DENDRO_731) - DENDRO_41*DENDRO_813*(-DENDRO_715 + DENDRO_817) + DENDRO_701*DENDRO_800 - 12*DENDRO_704 - DENDRO_812*(-DENDRO_705 + DENDRO_706 + DENDRO_707 - DENDRO_710) + DENDRO_814*(DENDRO_64*(DENDRO_628 + DENDRO_63*DENDRO_701) + DENDRO_713)) - alpha[pp]*(-At4[pp]*K[pp] + DENDRO_42*DENDRO_842 + DENDRO_840*DENDRO_857 + DENDRO_843*DENDRO_858) + beta0[pp]*agrad_0_At4[pp] + beta1[pp]*agrad_1_At4[pp] + beta2[pp]*agrad_2_At4[pp];
//--
At_rhs22[pp] = -At5[pp]*DENDRO_18 - At5[pp]*DENDRO_19 + At5[pp]*DENDRO_23 + DENDRO_15*DENDRO_27 + DENDRO_16*DENDRO_856 + DENDRO_802*(DENDRO_351*(-DENDRO_121*(DENDRO_117*DENDRO_265 + DENDRO_118*DENDRO_219 + DENDRO_120*DENDRO_247 + DENDRO_374) + DENDRO_132*DENDRO_168*DENDRO_427 - DENDRO_132*DENDRO_835 - DENDRO_146*(DENDRO_395 - 1.0*DENDRO_397) - DENDRO_146*(0.25*DENDRO_869 + 1.0*DENDRO_870) + DENDRO_179*DENDRO_424 + DENDRO_184*(DENDRO_386 - 1.0*DENDRO_388) + DENDRO_184*(0.25*DENDRO_871 + 1.0*DENDRO_872) - DENDRO_184*(-DENDRO_391*DENDRO_846 + DENDRO_845) - DENDRO_190*(DENDRO_147*DENDRO_824 + DENDRO_848) - DENDRO_190*(DENDRO_168*DENDRO_214 + DENDRO_866) - DENDRO_190*(DENDRO_174*DENDRO_817 + DENDRO_179*DENDRO_408) - DENDRO_190*(DENDRO_76*DENDRO_861 + 0.25*DENDRO_816) + DENDRO_197*(DENDRO_871 + DENDRO_872) + DENDRO_197*(DENDRO_132*DENDRO_179 + DENDRO_74*DENDRO_846) + DENDRO_204*DENDRO_319*DENDRO_824 - DENDRO_204*DENDRO_834 + 6.0*DENDRO_211*DENDRO_301*DENDRO_849 - DENDRO_211*DENDRO_836 - DENDRO_367*DENDRO_837 - DENDRO_410*(DENDRO_869 + DENDRO_870) - DENDRO_410*(DENDRO_132*DENDRO_817 + DENDRO_167*DENDRO_846) + DENDRO_420*(DENDRO_150 + DENDRO_379) + DENDRO_421*(DENDRO_217 + DENDRO_376) + DENDRO_423*DENDRO_817 + DENDRO_444) - 12*DENDRO_352 + DENDRO_354*DENDRO_84 + DENDRO_801*gt5[pp] - DENDRO_859*(-DENDRO_360 + DENDRO_846) - DENDRO_95*(DENDRO_363 - DENDRO_364 + DENDRO_365 - DENDRO_372)) - alpha[pp]*(At5[pp]*DENDRO_51*DENDRO_843 - At5[pp]*K[pp] + DENDRO_57*DENDRO_842 + DENDRO_840*DENDRO_858) + beta0[pp]*agrad_0_At5[pp] + beta1[pp]*agrad_1_At5[pp] + beta2[pp]*agrad_2_At5[pp];
//--
K_rhs[pp] = -DENDRO_334*chi[pp]*(-DENDRO_352 + DENDRO_873*(DENDRO_353*DENDRO_450 + DENDRO_383) + DENDRO_875*(DENDRO_384 + DENDRO_874*gt5[pp]) + DENDRO_94*(DENDRO_375*DENDRO_41 + DENDRO_378*DENDRO_41 + DENDRO_381*DENDRO_41 - DENDRO_64*(DENDRO_366 - DENDRO_371))) - DENDRO_336*chi[pp]*(-DENDRO_445 + DENDRO_82*(DENDRO_41*DENDRO_460 + DENDRO_41*DENDRO_461 + DENDRO_454 - DENDRO_64*(DENDRO_448 - DENDRO_453)) + DENDRO_875*(DENDRO_462 + DENDRO_874*gt3[pp]) + DENDRO_877*(DENDRO_459 + DENDRO_876*gt3[pp])) - DENDRO_338*chi[pp]*(DENDRO_112*(DENDRO_102 + DENDRO_41*DENDRO_509 + DENDRO_41*DENDRO_510 - DENDRO_64*(DENDRO_103 - DENDRO_110)) - DENDRO_58 + DENDRO_873*(DENDRO_450*DENDRO_66 + DENDRO_508) + DENDRO_877*(DENDRO_507 + DENDRO_876*gt0[pp])) + DENDRO_880*chi[pp]*(-DENDRO_704 + 0.5*DENDRO_875*(DENDRO_108*DENDRO_714 + DENDRO_700) + DENDRO_878*(DENDRO_41*DENDRO_697 + DENDRO_41*DENDRO_698 - DENDRO_64*(DENDRO_59 - DENDRO_709) + DENDRO_705) + DENDRO_879*(DENDRO_41*DENDRO_699 - DENDRO_64*(DENDRO_60 - DENDRO_711) + DENDRO_712)) - DENDRO_882*chi[pp]*(-DENDRO_617 + 0.5*DENDRO_873*(DENDRO_450*DENDRO_633 + DENDRO_591) + DENDRO_878*(DENDRO_41*DENDRO_588 + DENDRO_41*DENDRO_589 + DENDRO_620 - DENDRO_64*(DENDRO_61 - DENDRO_622)) + DENDRO_881*(DENDRO_41*DENDRO_593 + DENDRO_41*DENDRO_594 + DENDRO_627 - DENDRO_64*(DENDRO_60 - DENDRO_629))) + DENDRO_883*chi[pp]*(-DENDRO_787 + 0.5*DENDRO_877*(DENDRO_369*DENDRO_796 + DENDRO_763) + DENDRO_879*(DENDRO_41*DENDRO_764 - DENDRO_64*(DENDRO_61 - DENDRO_788) + DENDRO_789) + DENDRO_881*(DENDRO_41*DENDRO_765 + DENDRO_41*DENDRO_766 - DENDRO_64*(DENDRO_59 - DENDRO_794) + DENDRO_791)) + DENDRO_916*(At0[pp]*DENDRO_887*DENDRO_888 + At1[pp]*DENDRO_895*DENDRO_905 + At2[pp]*DENDRO_894*DENDRO_895 + At3[pp]*DENDRO_888*DENDRO_891 + At4[pp]*DENDRO_895*DENDRO_915 + At5[pp]*DENDRO_888*DENDRO_893 + pow(K[pp], 2)) + beta0[pp]*agrad_0_K[pp] + beta1[pp]*agrad_1_K[pp] + beta2[pp]*agrad_2_K[pp];
//--
Gt_rhs0[pp] = DENDRO_955;
//--
Gt_rhs1[pp] = -DENDRO_141*DENDRO_950*DENDRO_956 + DENDRO_16*DENDRO_271 - DENDRO_208*DENDRO_950*DENDRO_957 + DENDRO_233*DENDRO_4 - DENDRO_233*DENDRO_942 + DENDRO_253*DENDRO_7 - DENDRO_822*DENDRO_947 + DENDRO_918*DENDRO_956 + DENDRO_920*DENDRO_957 + DENDRO_927*(DENDRO_960 - DENDRO_961) - DENDRO_927*(-DENDRO_926*DENDRO_956 + DENDRO_959) - DENDRO_927*(-DENDRO_929*DENDRO_957 + DENDRO_958) + DENDRO_968;
//--
Gt_rhs2[pp] = DENDRO_970;
//--
B_rhs0[pp] = -B0[pp]*eta - DENDRO_954*lambda[3] + DENDRO_955 + lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] + beta2[pp]*agrad_2_B0[pp]);
//--
B_rhs1[pp] = -B1[pp]*eta + DENDRO_141*DENDRO_952 - DENDRO_16*DENDRO_413 + DENDRO_208*DENDRO_953 + DENDRO_230*DENDRO_947 - DENDRO_4*DENDRO_411 + DENDRO_411*DENDRO_942 - DENDRO_415*DENDRO_7 - DENDRO_905*DENDRO_918 - DENDRO_915*DENDRO_920 - DENDRO_927*(-DENDRO_960 + DENDRO_961) - DENDRO_927*(DENDRO_905*DENDRO_926 + DENDRO_959) - DENDRO_927*(DENDRO_915*DENDRO_929 + DENDRO_958) - DENDRO_967*lambda[3] + DENDRO_968 + lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] + beta2[pp]*agrad_2_B1[pp]);
//--
B_rhs2[pp] = -B2[pp]*eta - DENDRO_969*lambda[3] + DENDRO_970 + lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] + beta2[pp]*agrad_2_B2[pp]);
// Dendro: reduced ops: 4137
// Dendro: }}} 
//[[[end]]]
#endif



bssn::timer::t_rhs.stop();

       /* debugging */
        unsigned int qi = 46 - 1;
        unsigned int qj = 10 - 1;
        unsigned int qk = 60 - 1;
        unsigned int qidx = qi + nx*(qj + ny*qk);
        if (0 && qidx == pp) {
          std::cout << ".... end OPTIMIZED debug stuff..." << std::endl;
        }

      }
    }
  }


  if (bflag != 0) {

    bssn::timer::t_bdyc.start();
    
    
#ifdef BSSN_KERR_SCHILD_TEST

        freeze_bcs(a_rhs, sz, bflag);
        freeze_bcs(chi_rhs, sz, bflag);
        freeze_bcs(K_rhs, sz, bflag);
        freeze_bcs(b_rhs0, sz, bflag);
        freeze_bcs(b_rhs1, sz, bflag);
        freeze_bcs(b_rhs2, sz, bflag);
        freeze_bcs(Gt_rhs0, sz, bflag);
        freeze_bcs(Gt_rhs1, sz, bflag);
        freeze_bcs(Gt_rhs2, sz, bflag);
        freeze_bcs(B_rhs0, sz, bflag);
        freeze_bcs(B_rhs1, sz, bflag);
        freeze_bcs(B_rhs2, sz, bflag);
        freeze_bcs(At_rhs00, sz, bflag);
        freeze_bcs(At_rhs01, sz, bflag);
        freeze_bcs(At_rhs02, sz, bflag);
        freeze_bcs(At_rhs11, sz, bflag);
        freeze_bcs(At_rhs12, sz, bflag);
        freeze_bcs(At_rhs22, sz, bflag);
        freeze_bcs(gt_rhs00, sz, bflag);
        freeze_bcs(gt_rhs01, sz, bflag);
        freeze_bcs(gt_rhs02, sz, bflag);
        freeze_bcs(gt_rhs11, sz, bflag);
        freeze_bcs(gt_rhs12, sz, bflag);
        freeze_bcs(gt_rhs22, sz, bflag);

#else
        
        bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin, pmax,
             1.0, 1.0, sz, bflag);
        bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                1.0, 1.0, sz, bflag);
        bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax,
                1.0, 0.0, sz, bflag);

        bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin, pmax,
                1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin, pmax,
                1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin, pmax,
                1.0, 0.0, sz, bflag);

        bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                2.0, 0.0, sz, bflag);

        bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax,
                1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax,
                1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax,
                1.0, 0.0, sz, bflag);

        bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                2.0, 0.0, sz, bflag);

        bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                1.0, 1.0, sz, bflag);
        
#endif

    

    bssn::timer::t_bdyc.stop();
  }


bssn::timer::t_deriv.start();
// ko derives. 
/*[[[cog
import cog
import bssnDerivs as bssnDerivs
    
for var in bssnDerivs.KO:
    cog.outl("\t ko_deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
    cog.outl("\t ko_deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
    cog.outl("\t ko_deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))    
    
]]]*/
ko_deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
ko_deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
ko_deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
ko_deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
ko_deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
ko_deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
ko_deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
ko_deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
ko_deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
ko_deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
ko_deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
ko_deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
ko_deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
ko_deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
ko_deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
ko_deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
ko_deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
ko_deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
ko_deriv_x(grad_0_At0, At0, hx, sz, bflag);
ko_deriv_y(grad_1_At0, At0, hx, sz, bflag);
ko_deriv_z(grad_2_At0, At0, hx, sz, bflag);
ko_deriv_x(grad_0_At1, At1, hx, sz, bflag);
ko_deriv_y(grad_1_At1, At1, hx, sz, bflag);
ko_deriv_z(grad_2_At1, At1, hx, sz, bflag);
ko_deriv_x(grad_0_At2, At2, hx, sz, bflag);
ko_deriv_y(grad_1_At2, At2, hx, sz, bflag);
ko_deriv_z(grad_2_At2, At2, hx, sz, bflag);
ko_deriv_x(grad_0_At3, At3, hx, sz, bflag);
ko_deriv_y(grad_1_At3, At3, hx, sz, bflag);
ko_deriv_z(grad_2_At3, At3, hx, sz, bflag);
ko_deriv_x(grad_0_At4, At4, hx, sz, bflag);
ko_deriv_y(grad_1_At4, At4, hx, sz, bflag);
ko_deriv_z(grad_2_At4, At4, hx, sz, bflag);
ko_deriv_x(grad_0_At5, At5, hx, sz, bflag);
ko_deriv_y(grad_1_At5, At5, hx, sz, bflag);
ko_deriv_z(grad_2_At5, At5, hx, sz, bflag);
ko_deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
ko_deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
ko_deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
ko_deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
ko_deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
ko_deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
ko_deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
ko_deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
ko_deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
ko_deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
ko_deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
ko_deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
ko_deriv_x(grad_0_chi, chi, hx, sz, bflag);
ko_deriv_y(grad_1_chi, chi, hx, sz, bflag);
ko_deriv_z(grad_2_chi, chi, hx, sz, bflag);
ko_deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
ko_deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
ko_deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
ko_deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
ko_deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
ko_deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
ko_deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
ko_deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
ko_deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
ko_deriv_x(grad_0_K, K, hx, sz, bflag);
ko_deriv_y(grad_1_K, K, hx, sz, bflag);
ko_deriv_z(grad_2_K, K, hx, sz, bflag);
ko_deriv_x(grad_0_B0, B0, hx, sz, bflag);
ko_deriv_y(grad_1_B0, B0, hx, sz, bflag);
ko_deriv_z(grad_2_B0, B0, hx, sz, bflag);
ko_deriv_x(grad_0_B1, B1, hx, sz, bflag);
ko_deriv_y(grad_1_B1, B1, hx, sz, bflag);
ko_deriv_z(grad_2_B1, B1, hx, sz, bflag);
ko_deriv_x(grad_0_B2, B2, hx, sz, bflag);
ko_deriv_y(grad_1_B2, B2, hx, sz, bflag);
ko_deriv_z(grad_2_B2, B2, hx, sz, bflag);
//[[[end]]]

bssn::timer::t_deriv.stop();

bssn::timer::t_rhs.start();

  const  double sigma = KO_DISS_SIGMA;


  for (unsigned int k = 3; k < nz-3; k++) {
    for (unsigned int j = 3; j < ny-3; j++) {
      for (unsigned int i = 3; i < nx-3; i++) {
        pp = i + nx*(j + ny*k);

        a_rhs[pp]  += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
        b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
        b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
        b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

        gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
        gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
        gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
        gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
        gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
        gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

        chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

        At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
        At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
        At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
        At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
        At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
        At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

        K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

        Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
        Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
        Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

        B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
        B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
        B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
      }
    }
  }

  bssn::timer::t_rhs.stop();



bssn::timer::t_deriv.start();
/*[[[cog
import cog
import bssnDerivs as bssnDerivs

for deriv in bssnDerivs.FUNC_D_I:
    cog.outl("\t free(%s);" %(deriv))
    
for deriv in bssnDerivs.FUNC_D_IJ:
    cog.outl("\t free(%s);" %(deriv))

for deriv in bssnDerivs.FUNC_AD_I:
    cog.outl("\t free(%s);" %(deriv))   
    
]]]*/
free(grad_0_alpha);
free(grad_1_alpha);
free(grad_2_alpha);
free(grad_0_beta0);
free(grad_1_beta0);
free(grad_2_beta0);
free(grad_0_beta1);
free(grad_1_beta1);
free(grad_2_beta1);
free(grad_0_beta2);
free(grad_1_beta2);
free(grad_2_beta2);
free(grad_0_B0);
free(grad_1_B0);
free(grad_2_B0);
free(grad_0_B1);
free(grad_1_B1);
free(grad_2_B1);
free(grad_0_B2);
free(grad_1_B2);
free(grad_2_B2);
free(grad_0_chi);
free(grad_1_chi);
free(grad_2_chi);
free(grad_0_Gt0);
free(grad_1_Gt0);
free(grad_2_Gt0);
free(grad_0_Gt1);
free(grad_1_Gt1);
free(grad_2_Gt1);
free(grad_0_Gt2);
free(grad_1_Gt2);
free(grad_2_Gt2);
free(grad_0_K);
free(grad_1_K);
free(grad_2_K);
free(grad_0_gt0);
free(grad_1_gt0);
free(grad_2_gt0);
free(grad_0_gt1);
free(grad_1_gt1);
free(grad_2_gt1);
free(grad_0_gt2);
free(grad_1_gt2);
free(grad_2_gt2);
free(grad_0_gt3);
free(grad_1_gt3);
free(grad_2_gt3);
free(grad_0_gt4);
free(grad_1_gt4);
free(grad_2_gt4);
free(grad_0_gt5);
free(grad_1_gt5);
free(grad_2_gt5);
free(grad_0_At0);
free(grad_1_At0);
free(grad_2_At0);
free(grad_0_At1);
free(grad_1_At1);
free(grad_2_At1);
free(grad_0_At2);
free(grad_1_At2);
free(grad_2_At2);
free(grad_0_At3);
free(grad_1_At3);
free(grad_2_At3);
free(grad_0_At4);
free(grad_1_At4);
free(grad_2_At4);
free(grad_0_At5);
free(grad_1_At5);
free(grad_2_At5);
free(grad2_0_0_gt0);
free(grad2_0_1_gt0);
free(grad2_0_2_gt0);
free(grad2_1_1_gt0);
free(grad2_1_2_gt0);
free(grad2_2_2_gt0);
free(grad2_0_0_gt1);
free(grad2_0_1_gt1);
free(grad2_0_2_gt1);
free(grad2_1_1_gt1);
free(grad2_1_2_gt1);
free(grad2_2_2_gt1);
free(grad2_0_0_gt2);
free(grad2_0_1_gt2);
free(grad2_0_2_gt2);
free(grad2_1_1_gt2);
free(grad2_1_2_gt2);
free(grad2_2_2_gt2);
free(grad2_0_0_gt3);
free(grad2_0_1_gt3);
free(grad2_0_2_gt3);
free(grad2_1_1_gt3);
free(grad2_1_2_gt3);
free(grad2_2_2_gt3);
free(grad2_0_0_gt4);
free(grad2_0_1_gt4);
free(grad2_0_2_gt4);
free(grad2_1_1_gt4);
free(grad2_1_2_gt4);
free(grad2_2_2_gt4);
free(grad2_0_0_gt5);
free(grad2_0_1_gt5);
free(grad2_0_2_gt5);
free(grad2_1_1_gt5);
free(grad2_1_2_gt5);
free(grad2_2_2_gt5);
free(grad2_0_0_chi);
free(grad2_0_1_chi);
free(grad2_0_2_chi);
free(grad2_1_1_chi);
free(grad2_1_2_chi);
free(grad2_2_2_chi);
free(grad2_0_0_alpha);
free(grad2_0_1_alpha);
free(grad2_0_2_alpha);
free(grad2_1_1_alpha);
free(grad2_1_2_alpha);
free(grad2_2_2_alpha);
free(grad2_0_0_beta0);
free(grad2_0_1_beta0);
free(grad2_0_2_beta0);
free(grad2_1_1_beta0);
free(grad2_1_2_beta0);
free(grad2_2_2_beta0);
free(grad2_0_0_beta1);
free(grad2_0_1_beta1);
free(grad2_0_2_beta1);
free(grad2_1_1_beta1);
free(grad2_1_2_beta1);
free(grad2_2_2_beta1);
free(grad2_0_0_beta2);
free(grad2_0_1_beta2);
free(grad2_0_2_beta2);
free(grad2_1_1_beta2);
free(grad2_1_2_beta2);
free(grad2_2_2_beta2);
free(agrad_0_gt0);
free(agrad_1_gt0);
free(agrad_2_gt0);
free(agrad_0_gt1);
free(agrad_1_gt1);
free(agrad_2_gt1);
free(agrad_0_gt2);
free(agrad_1_gt2);
free(agrad_2_gt2);
free(agrad_0_gt3);
free(agrad_1_gt3);
free(agrad_2_gt3);
free(agrad_0_gt4);
free(agrad_1_gt4);
free(agrad_2_gt4);
free(agrad_0_gt5);
free(agrad_1_gt5);
free(agrad_2_gt5);
free(agrad_0_At0);
free(agrad_1_At0);
free(agrad_2_At0);
free(agrad_0_At1);
free(agrad_1_At1);
free(agrad_2_At1);
free(agrad_0_At2);
free(agrad_1_At2);
free(agrad_2_At2);
free(agrad_0_At3);
free(agrad_1_At3);
free(agrad_2_At3);
free(agrad_0_At4);
free(agrad_1_At4);
free(agrad_2_At4);
free(agrad_0_At5);
free(agrad_1_At5);
free(agrad_2_At5);
free(agrad_0_alpha);
free(agrad_1_alpha);
free(agrad_2_alpha);
free(agrad_0_beta0);
free(agrad_1_beta0);
free(agrad_2_beta0);
free(agrad_0_beta1);
free(agrad_1_beta1);
free(agrad_2_beta1);
free(agrad_0_beta2);
free(agrad_1_beta2);
free(agrad_2_beta2);
free(agrad_0_chi);
free(agrad_1_chi);
free(agrad_2_chi);
free(agrad_0_Gt0);
free(agrad_1_Gt0);
free(agrad_2_Gt0);
free(agrad_0_Gt1);
free(agrad_1_Gt1);
free(agrad_2_Gt1);
free(agrad_0_Gt2);
free(agrad_1_Gt2);
free(agrad_2_Gt2);
free(agrad_0_K);
free(agrad_1_K);
free(agrad_2_K);
free(agrad_0_B0);
free(agrad_1_B0);
free(agrad_2_B0);
free(agrad_0_B1);
free(agrad_1_B1);
free(agrad_2_B1);
free(agrad_0_B2);
free(agrad_1_B2);
free(agrad_2_B2);
//[[[end]]]
bssn::timer::t_deriv.stop();

#if 0
  for (unsigned int m = 0; m < 24; m++) {
    std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
  }
#endif



}



void bssnrhs_sep(double **unzipVarsRHS, const double **uZipVars,
             const unsigned int& offset,
             const double *pmin, const double *pmax, const unsigned int *sz,
             const unsigned int& bflag)
{



    const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
    const double *chi = &uZipVars[VAR::U_CHI][offset];
    const double *K = &uZipVars[VAR::U_K][offset];
    const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
    const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
    const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
    const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
    const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
    const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
    const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
    const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
    const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
    const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
    const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
    const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
    const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
    const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
    const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
    const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
    const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
    const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
    const double *B0 = &uZipVars[VAR::U_B0][offset];
    const double *B1 = &uZipVars[VAR::U_B1][offset];
    const double *B2 = &uZipVars[VAR::U_B2][offset];

    double *a_rhs = &unzipVarsRHS[VAR::U_ALPHA][offset];
    double *chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    double *K_rhs = &unzipVarsRHS[VAR::U_K][offset];
    double *gt_rhs00 = &unzipVarsRHS[VAR::U_SYMGT0][offset];
    double *gt_rhs01 = &unzipVarsRHS[VAR::U_SYMGT1][offset];
    double *gt_rhs02 = &unzipVarsRHS[VAR::U_SYMGT2][offset];
    double *gt_rhs11 = &unzipVarsRHS[VAR::U_SYMGT3][offset];
    double *gt_rhs12 = &unzipVarsRHS[VAR::U_SYMGT4][offset];
    double *gt_rhs22 = &unzipVarsRHS[VAR::U_SYMGT5][offset];
    double *b_rhs0 = &unzipVarsRHS[VAR::U_BETA0][offset];
    double *b_rhs1 = &unzipVarsRHS[VAR::U_BETA1][offset];
    double *b_rhs2 = &unzipVarsRHS[VAR::U_BETA2][offset];
    double *At_rhs00 = &unzipVarsRHS[VAR::U_SYMAT0][offset];
    double *At_rhs01 = &unzipVarsRHS[VAR::U_SYMAT1][offset];
    double *At_rhs02 = &unzipVarsRHS[VAR::U_SYMAT2][offset];
    double *At_rhs11 = &unzipVarsRHS[VAR::U_SYMAT3][offset];
    double *At_rhs12 = &unzipVarsRHS[VAR::U_SYMAT4][offset];
    double *At_rhs22 = &unzipVarsRHS[VAR::U_SYMAT5][offset];
    double *Gt_rhs0 = &unzipVarsRHS[VAR::U_GT0][offset];
    double *Gt_rhs1 = &unzipVarsRHS[VAR::U_GT1][offset];
    double *Gt_rhs2 = &unzipVarsRHS[VAR::U_GT2][offset];
    double *B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
    double *B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
    double *B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx = (pmax[0] - pmin[0]) / (nx - 1);
    double hy = (pmax[1] - pmin[1]) / (ny - 1);
    double hz = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int lambda[4] = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                    BSSN_LAMBDA[2], BSSN_LAMBDA[3]};
    const double lambda_f[2] = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};
    const double eta_power[2] = {BSSN_ETA_POWER[0], BSSN_ETA_POWER[1]};



    int idx[3];

    unsigned int n = sz[0]*sz[1]*sz[2];

    double * CalGt0 =new double[n];
    double * CalGt1 =new double[n];
    double * CalGt2 =new double[n];

    double *Gt_rhs_s1_0 =new double [n];
    double *Gt_rhs_s1_1 =new double [n];
    double *Gt_rhs_s1_2 =new double [n];

    double *Gt_rhs_s2_0 =new double [n];
    double *Gt_rhs_s2_1 =new double [n];
    double *Gt_rhs_s2_2 =new double [n];

    double *Gt_rhs_s3_0 =new double [n];
    double *Gt_rhs_s3_1 =new double [n];
    double *Gt_rhs_s3_2 =new double [n];

    double *Gt_rhs_s4_0 =new double [n];
    double *Gt_rhs_s4_1 =new double [n];
    double *Gt_rhs_s4_2 =new double [n];

    double *Gt_rhs_s5_0 =new double [n];
    double *Gt_rhs_s5_1 =new double [n];
    double *Gt_rhs_s5_2 =new double [n];

    double *Gt_rhs_s6_0 =new double [n];
    double *Gt_rhs_s6_1 =new double [n];
    double *Gt_rhs_s6_2 =new double [n];

    double *Gt_rhs_s7_0 =new double [n];
    double *Gt_rhs_s7_1 =new double [n];
    double *Gt_rhs_s7_2 =new double [n];



    bssn::timer::t_deriv.start();

/*[[[cog
import cog
import bssnDerivs as bssnDerivs

for deriv in bssnDerivs.FUNC_D_I:
    cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")
    
for deriv in bssnDerivs.FUNC_D_IJ:
    cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")

for deriv in bssnDerivs.FUNC_AD_I:
    cog.outl("\t double* "+deriv+" = (double*)malloc(sizeof(double)*n);")   
    
    
for var in bssnDerivs.D:
    cog.outl("\t deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
    cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
    cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))
    
    if var in bssnDerivs.DD:
        cog.outl("\t deriv_xx(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[0] + var ,var))
        cog.outl("\t deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[1] + var , bssnDerivs.PREFIX_D[0] + var ))
        cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[2] + var , bssnDerivs.PREFIX_D[0] + var ))
        
        cog.outl("\t deriv_yy(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[3] + var ,var))
        cog.outl("\t deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[4] + var , bssnDerivs.PREFIX_D[1] + var))
        
        cog.outl("\t deriv_zz(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_DD[5] + var ,var))
        
    if var in bssnDerivs.AD:
        cog.outl("\t adv_deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_AD[0] + var ,var))
        cog.outl("\t adv_deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_AD[1] + var ,var))
        cog.outl("\t adv_deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_AD[2] + var ,var))
        
    
    
]]]*/
double* grad_0_alpha = (double*)malloc(sizeof(double)*n);
double* grad_1_alpha = (double*)malloc(sizeof(double)*n);
double* grad_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad_0_beta0 = (double*)malloc(sizeof(double)*n);
double* grad_1_beta0 = (double*)malloc(sizeof(double)*n);
double* grad_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad_0_beta1 = (double*)malloc(sizeof(double)*n);
double* grad_1_beta1 = (double*)malloc(sizeof(double)*n);
double* grad_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad_0_beta2 = (double*)malloc(sizeof(double)*n);
double* grad_1_beta2 = (double*)malloc(sizeof(double)*n);
double* grad_2_beta2 = (double*)malloc(sizeof(double)*n);
double* grad_0_B0 = (double*)malloc(sizeof(double)*n);
double* grad_1_B0 = (double*)malloc(sizeof(double)*n);
double* grad_2_B0 = (double*)malloc(sizeof(double)*n);
double* grad_0_B1 = (double*)malloc(sizeof(double)*n);
double* grad_1_B1 = (double*)malloc(sizeof(double)*n);
double* grad_2_B1 = (double*)malloc(sizeof(double)*n);
double* grad_0_B2 = (double*)malloc(sizeof(double)*n);
double* grad_1_B2 = (double*)malloc(sizeof(double)*n);
double* grad_2_B2 = (double*)malloc(sizeof(double)*n);
double* grad_0_chi = (double*)malloc(sizeof(double)*n);
double* grad_1_chi = (double*)malloc(sizeof(double)*n);
double* grad_2_chi = (double*)malloc(sizeof(double)*n);
double* grad_0_Gt0 = (double*)malloc(sizeof(double)*n);
double* grad_1_Gt0 = (double*)malloc(sizeof(double)*n);
double* grad_2_Gt0 = (double*)malloc(sizeof(double)*n);
double* grad_0_Gt1 = (double*)malloc(sizeof(double)*n);
double* grad_1_Gt1 = (double*)malloc(sizeof(double)*n);
double* grad_2_Gt1 = (double*)malloc(sizeof(double)*n);
double* grad_0_Gt2 = (double*)malloc(sizeof(double)*n);
double* grad_1_Gt2 = (double*)malloc(sizeof(double)*n);
double* grad_2_Gt2 = (double*)malloc(sizeof(double)*n);
double* grad_0_K = (double*)malloc(sizeof(double)*n);
double* grad_1_K = (double*)malloc(sizeof(double)*n);
double* grad_2_K = (double*)malloc(sizeof(double)*n);
double* grad_0_gt0 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt0 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt1 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt1 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt2 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt2 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt3 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt3 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt4 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt4 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad_0_gt5 = (double*)malloc(sizeof(double)*n);
double* grad_1_gt5 = (double*)malloc(sizeof(double)*n);
double* grad_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad_0_At0 = (double*)malloc(sizeof(double)*n);
double* grad_1_At0 = (double*)malloc(sizeof(double)*n);
double* grad_2_At0 = (double*)malloc(sizeof(double)*n);
double* grad_0_At1 = (double*)malloc(sizeof(double)*n);
double* grad_1_At1 = (double*)malloc(sizeof(double)*n);
double* grad_2_At1 = (double*)malloc(sizeof(double)*n);
double* grad_0_At2 = (double*)malloc(sizeof(double)*n);
double* grad_1_At2 = (double*)malloc(sizeof(double)*n);
double* grad_2_At2 = (double*)malloc(sizeof(double)*n);
double* grad_0_At3 = (double*)malloc(sizeof(double)*n);
double* grad_1_At3 = (double*)malloc(sizeof(double)*n);
double* grad_2_At3 = (double*)malloc(sizeof(double)*n);
double* grad_0_At4 = (double*)malloc(sizeof(double)*n);
double* grad_1_At4 = (double*)malloc(sizeof(double)*n);
double* grad_2_At4 = (double*)malloc(sizeof(double)*n);
double* grad_0_At5 = (double*)malloc(sizeof(double)*n);
double* grad_1_At5 = (double*)malloc(sizeof(double)*n);
double* grad_2_At5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt3 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt4 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_gt5 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_chi = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_chi = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_chi = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_chi = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_chi = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_chi = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_alpha = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_beta0 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_beta1 = (double*)malloc(sizeof(double)*n);
double* grad2_0_0_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_1_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_0_2_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_1_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_1_2_beta2 = (double*)malloc(sizeof(double)*n);
double* grad2_2_2_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt3 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt3 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt3 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt4 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt4 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt4 = (double*)malloc(sizeof(double)*n);
double* agrad_0_gt5 = (double*)malloc(sizeof(double)*n);
double* agrad_1_gt5 = (double*)malloc(sizeof(double)*n);
double* agrad_2_gt5 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At3 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At3 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At3 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At4 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At4 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At4 = (double*)malloc(sizeof(double)*n);
double* agrad_0_At5 = (double*)malloc(sizeof(double)*n);
double* agrad_1_At5 = (double*)malloc(sizeof(double)*n);
double* agrad_2_At5 = (double*)malloc(sizeof(double)*n);
double* agrad_0_alpha = (double*)malloc(sizeof(double)*n);
double* agrad_1_alpha = (double*)malloc(sizeof(double)*n);
double* agrad_2_alpha = (double*)malloc(sizeof(double)*n);
double* agrad_0_beta0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_beta0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_beta0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_beta1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_beta1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_beta1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_beta2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_chi = (double*)malloc(sizeof(double)*n);
double* agrad_1_chi = (double*)malloc(sizeof(double)*n);
double* agrad_2_chi = (double*)malloc(sizeof(double)*n);
double* agrad_0_Gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_Gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_Gt0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_Gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_Gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_Gt1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_Gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_Gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_Gt2 = (double*)malloc(sizeof(double)*n);
double* agrad_0_K = (double*)malloc(sizeof(double)*n);
double* agrad_1_K = (double*)malloc(sizeof(double)*n);
double* agrad_2_K = (double*)malloc(sizeof(double)*n);
double* agrad_0_B0 = (double*)malloc(sizeof(double)*n);
double* agrad_1_B0 = (double*)malloc(sizeof(double)*n);
double* agrad_2_B0 = (double*)malloc(sizeof(double)*n);
double* agrad_0_B1 = (double*)malloc(sizeof(double)*n);
double* agrad_1_B1 = (double*)malloc(sizeof(double)*n);
double* agrad_2_B1 = (double*)malloc(sizeof(double)*n);
double* agrad_0_B2 = (double*)malloc(sizeof(double)*n);
double* agrad_1_B2 = (double*)malloc(sizeof(double)*n);
double* agrad_2_B2 = (double*)malloc(sizeof(double)*n);
deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
deriv_xx(grad2_0_0_alpha, alpha, hx, sz, bflag);
deriv_y(grad2_0_1_alpha, grad_0_alpha, hx, sz, bflag);
deriv_z(grad2_0_2_alpha, grad_0_alpha, hx, sz, bflag);
deriv_yy(grad2_1_1_alpha, alpha, hx, sz, bflag);
deriv_z(grad2_1_2_alpha, grad_1_alpha, hx, sz, bflag);
deriv_zz(grad2_2_2_alpha, alpha, hx, sz, bflag);
adv_deriv_x(agrad_0_alpha, alpha, hx, sz, bflag);
adv_deriv_y(agrad_1_alpha, alpha, hx, sz, bflag);
adv_deriv_z(agrad_2_alpha, alpha, hx, sz, bflag);
deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
deriv_xx(grad2_0_0_beta0, beta0, hx, sz, bflag);
deriv_y(grad2_0_1_beta0, grad_0_beta0, hx, sz, bflag);
deriv_z(grad2_0_2_beta0, grad_0_beta0, hx, sz, bflag);
deriv_yy(grad2_1_1_beta0, beta0, hx, sz, bflag);
deriv_z(grad2_1_2_beta0, grad_1_beta0, hx, sz, bflag);
deriv_zz(grad2_2_2_beta0, beta0, hx, sz, bflag);
adv_deriv_x(agrad_0_beta0, beta0, hx, sz, bflag);
adv_deriv_y(agrad_1_beta0, beta0, hx, sz, bflag);
adv_deriv_z(agrad_2_beta0, beta0, hx, sz, bflag);
deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
deriv_xx(grad2_0_0_beta1, beta1, hx, sz, bflag);
deriv_y(grad2_0_1_beta1, grad_0_beta1, hx, sz, bflag);
deriv_z(grad2_0_2_beta1, grad_0_beta1, hx, sz, bflag);
deriv_yy(grad2_1_1_beta1, beta1, hx, sz, bflag);
deriv_z(grad2_1_2_beta1, grad_1_beta1, hx, sz, bflag);
deriv_zz(grad2_2_2_beta1, beta1, hx, sz, bflag);
adv_deriv_x(agrad_0_beta1, beta1, hx, sz, bflag);
adv_deriv_y(agrad_1_beta1, beta1, hx, sz, bflag);
adv_deriv_z(agrad_2_beta1, beta1, hx, sz, bflag);
deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
deriv_xx(grad2_0_0_beta2, beta2, hx, sz, bflag);
deriv_y(grad2_0_1_beta2, grad_0_beta2, hx, sz, bflag);
deriv_z(grad2_0_2_beta2, grad_0_beta2, hx, sz, bflag);
deriv_yy(grad2_1_1_beta2, beta2, hx, sz, bflag);
deriv_z(grad2_1_2_beta2, grad_1_beta2, hx, sz, bflag);
deriv_zz(grad2_2_2_beta2, beta2, hx, sz, bflag);
adv_deriv_x(agrad_0_beta2, beta2, hx, sz, bflag);
adv_deriv_y(agrad_1_beta2, beta2, hx, sz, bflag);
adv_deriv_z(agrad_2_beta2, beta2, hx, sz, bflag);
deriv_x(grad_0_B0, B0, hx, sz, bflag);
deriv_y(grad_1_B0, B0, hx, sz, bflag);
deriv_z(grad_2_B0, B0, hx, sz, bflag);
adv_deriv_x(agrad_0_B0, B0, hx, sz, bflag);
adv_deriv_y(agrad_1_B0, B0, hx, sz, bflag);
adv_deriv_z(agrad_2_B0, B0, hx, sz, bflag);
deriv_x(grad_0_B1, B1, hx, sz, bflag);
deriv_y(grad_1_B1, B1, hx, sz, bflag);
deriv_z(grad_2_B1, B1, hx, sz, bflag);
adv_deriv_x(agrad_0_B1, B1, hx, sz, bflag);
adv_deriv_y(agrad_1_B1, B1, hx, sz, bflag);
adv_deriv_z(agrad_2_B1, B1, hx, sz, bflag);
deriv_x(grad_0_B2, B2, hx, sz, bflag);
deriv_y(grad_1_B2, B2, hx, sz, bflag);
deriv_z(grad_2_B2, B2, hx, sz, bflag);
adv_deriv_x(agrad_0_B2, B2, hx, sz, bflag);
adv_deriv_y(agrad_1_B2, B2, hx, sz, bflag);
adv_deriv_z(agrad_2_B2, B2, hx, sz, bflag);
deriv_x(grad_0_chi, chi, hx, sz, bflag);
deriv_y(grad_1_chi, chi, hx, sz, bflag);
deriv_z(grad_2_chi, chi, hx, sz, bflag);
deriv_xx(grad2_0_0_chi, chi, hx, sz, bflag);
deriv_y(grad2_0_1_chi, grad_0_chi, hx, sz, bflag);
deriv_z(grad2_0_2_chi, grad_0_chi, hx, sz, bflag);
deriv_yy(grad2_1_1_chi, chi, hx, sz, bflag);
deriv_z(grad2_1_2_chi, grad_1_chi, hx, sz, bflag);
deriv_zz(grad2_2_2_chi, chi, hx, sz, bflag);
adv_deriv_x(agrad_0_chi, chi, hx, sz, bflag);
adv_deriv_y(agrad_1_chi, chi, hx, sz, bflag);
adv_deriv_z(agrad_2_chi, chi, hx, sz, bflag);
deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
adv_deriv_x(agrad_0_Gt0, Gt0, hx, sz, bflag);
adv_deriv_y(agrad_1_Gt0, Gt0, hx, sz, bflag);
adv_deriv_z(agrad_2_Gt0, Gt0, hx, sz, bflag);
deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
adv_deriv_x(agrad_0_Gt1, Gt1, hx, sz, bflag);
adv_deriv_y(agrad_1_Gt1, Gt1, hx, sz, bflag);
adv_deriv_z(agrad_2_Gt1, Gt1, hx, sz, bflag);
deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
adv_deriv_x(agrad_0_Gt2, Gt2, hx, sz, bflag);
adv_deriv_y(agrad_1_Gt2, Gt2, hx, sz, bflag);
adv_deriv_z(agrad_2_Gt2, Gt2, hx, sz, bflag);
deriv_x(grad_0_K, K, hx, sz, bflag);
deriv_y(grad_1_K, K, hx, sz, bflag);
deriv_z(grad_2_K, K, hx, sz, bflag);
adv_deriv_x(agrad_0_K, K, hx, sz, bflag);
adv_deriv_y(agrad_1_K, K, hx, sz, bflag);
adv_deriv_z(agrad_2_K, K, hx, sz, bflag);
deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
deriv_xx(grad2_0_0_gt0, gt0, hx, sz, bflag);
deriv_y(grad2_0_1_gt0, grad_0_gt0, hx, sz, bflag);
deriv_z(grad2_0_2_gt0, grad_0_gt0, hx, sz, bflag);
deriv_yy(grad2_1_1_gt0, gt0, hx, sz, bflag);
deriv_z(grad2_1_2_gt0, grad_1_gt0, hx, sz, bflag);
deriv_zz(grad2_2_2_gt0, gt0, hx, sz, bflag);
adv_deriv_x(agrad_0_gt0, gt0, hx, sz, bflag);
adv_deriv_y(agrad_1_gt0, gt0, hx, sz, bflag);
adv_deriv_z(agrad_2_gt0, gt0, hx, sz, bflag);
deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
deriv_xx(grad2_0_0_gt1, gt1, hx, sz, bflag);
deriv_y(grad2_0_1_gt1, grad_0_gt1, hx, sz, bflag);
deriv_z(grad2_0_2_gt1, grad_0_gt1, hx, sz, bflag);
deriv_yy(grad2_1_1_gt1, gt1, hx, sz, bflag);
deriv_z(grad2_1_2_gt1, grad_1_gt1, hx, sz, bflag);
deriv_zz(grad2_2_2_gt1, gt1, hx, sz, bflag);
adv_deriv_x(agrad_0_gt1, gt1, hx, sz, bflag);
adv_deriv_y(agrad_1_gt1, gt1, hx, sz, bflag);
adv_deriv_z(agrad_2_gt1, gt1, hx, sz, bflag);
deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
deriv_xx(grad2_0_0_gt2, gt2, hx, sz, bflag);
deriv_y(grad2_0_1_gt2, grad_0_gt2, hx, sz, bflag);
deriv_z(grad2_0_2_gt2, grad_0_gt2, hx, sz, bflag);
deriv_yy(grad2_1_1_gt2, gt2, hx, sz, bflag);
deriv_z(grad2_1_2_gt2, grad_1_gt2, hx, sz, bflag);
deriv_zz(grad2_2_2_gt2, gt2, hx, sz, bflag);
adv_deriv_x(agrad_0_gt2, gt2, hx, sz, bflag);
adv_deriv_y(agrad_1_gt2, gt2, hx, sz, bflag);
adv_deriv_z(agrad_2_gt2, gt2, hx, sz, bflag);
deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
deriv_xx(grad2_0_0_gt3, gt3, hx, sz, bflag);
deriv_y(grad2_0_1_gt3, grad_0_gt3, hx, sz, bflag);
deriv_z(grad2_0_2_gt3, grad_0_gt3, hx, sz, bflag);
deriv_yy(grad2_1_1_gt3, gt3, hx, sz, bflag);
deriv_z(grad2_1_2_gt3, grad_1_gt3, hx, sz, bflag);
deriv_zz(grad2_2_2_gt3, gt3, hx, sz, bflag);
adv_deriv_x(agrad_0_gt3, gt3, hx, sz, bflag);
adv_deriv_y(agrad_1_gt3, gt3, hx, sz, bflag);
adv_deriv_z(agrad_2_gt3, gt3, hx, sz, bflag);
deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
deriv_xx(grad2_0_0_gt4, gt4, hx, sz, bflag);
deriv_y(grad2_0_1_gt4, grad_0_gt4, hx, sz, bflag);
deriv_z(grad2_0_2_gt4, grad_0_gt4, hx, sz, bflag);
deriv_yy(grad2_1_1_gt4, gt4, hx, sz, bflag);
deriv_z(grad2_1_2_gt4, grad_1_gt4, hx, sz, bflag);
deriv_zz(grad2_2_2_gt4, gt4, hx, sz, bflag);
adv_deriv_x(agrad_0_gt4, gt4, hx, sz, bflag);
adv_deriv_y(agrad_1_gt4, gt4, hx, sz, bflag);
adv_deriv_z(agrad_2_gt4, gt4, hx, sz, bflag);
deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
deriv_xx(grad2_0_0_gt5, gt5, hx, sz, bflag);
deriv_y(grad2_0_1_gt5, grad_0_gt5, hx, sz, bflag);
deriv_z(grad2_0_2_gt5, grad_0_gt5, hx, sz, bflag);
deriv_yy(grad2_1_1_gt5, gt5, hx, sz, bflag);
deriv_z(grad2_1_2_gt5, grad_1_gt5, hx, sz, bflag);
deriv_zz(grad2_2_2_gt5, gt5, hx, sz, bflag);
adv_deriv_x(agrad_0_gt5, gt5, hx, sz, bflag);
adv_deriv_y(agrad_1_gt5, gt5, hx, sz, bflag);
adv_deriv_z(agrad_2_gt5, gt5, hx, sz, bflag);
deriv_x(grad_0_At0, At0, hx, sz, bflag);
deriv_y(grad_1_At0, At0, hx, sz, bflag);
deriv_z(grad_2_At0, At0, hx, sz, bflag);
adv_deriv_x(agrad_0_At0, At0, hx, sz, bflag);
adv_deriv_y(agrad_1_At0, At0, hx, sz, bflag);
adv_deriv_z(agrad_2_At0, At0, hx, sz, bflag);
deriv_x(grad_0_At1, At1, hx, sz, bflag);
deriv_y(grad_1_At1, At1, hx, sz, bflag);
deriv_z(grad_2_At1, At1, hx, sz, bflag);
adv_deriv_x(agrad_0_At1, At1, hx, sz, bflag);
adv_deriv_y(agrad_1_At1, At1, hx, sz, bflag);
adv_deriv_z(agrad_2_At1, At1, hx, sz, bflag);
deriv_x(grad_0_At2, At2, hx, sz, bflag);
deriv_y(grad_1_At2, At2, hx, sz, bflag);
deriv_z(grad_2_At2, At2, hx, sz, bflag);
adv_deriv_x(agrad_0_At2, At2, hx, sz, bflag);
adv_deriv_y(agrad_1_At2, At2, hx, sz, bflag);
adv_deriv_z(agrad_2_At2, At2, hx, sz, bflag);
deriv_x(grad_0_At3, At3, hx, sz, bflag);
deriv_y(grad_1_At3, At3, hx, sz, bflag);
deriv_z(grad_2_At3, At3, hx, sz, bflag);
adv_deriv_x(agrad_0_At3, At3, hx, sz, bflag);
adv_deriv_y(agrad_1_At3, At3, hx, sz, bflag);
adv_deriv_z(agrad_2_At3, At3, hx, sz, bflag);
deriv_x(grad_0_At4, At4, hx, sz, bflag);
deriv_y(grad_1_At4, At4, hx, sz, bflag);
deriv_z(grad_2_At4, At4, hx, sz, bflag);
adv_deriv_x(agrad_0_At4, At4, hx, sz, bflag);
adv_deriv_y(agrad_1_At4, At4, hx, sz, bflag);
adv_deriv_z(agrad_2_At4, At4, hx, sz, bflag);
deriv_x(grad_0_At5, At5, hx, sz, bflag);
deriv_y(grad_1_At5, At5, hx, sz, bflag);
deriv_z(grad_2_At5, At5, hx, sz, bflag);
adv_deriv_x(agrad_0_At5, At5, hx, sz, bflag);
adv_deriv_y(agrad_1_At5, At5, hx, sz, bflag);
adv_deriv_z(agrad_2_At5, At5, hx, sz, bflag);
//[[[end]]]


    bssn::timer::t_deriv.stop();

    register double x;
    register double y;
    register double z;
    register unsigned int pp;

    double r_coord;
    double eta;

    //cout << "begin loop" << endl;
    bssn::timer::t_rhs.start();
    
    bssn::timer::t_rhs_a.start();
    
      for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
         
// a_rhs

/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.a_rhs]
vnames = ['a_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 12 
// Dendro: printing temp variables

// Dendro: printing variables
//--
a_rhs[pp] = -2*K[pp]*alpha[pp] + lambda[0]*(beta0[pp]*agrad_0_alpha[pp] + beta1[pp]*agrad_1_alpha[pp] + beta2[pp]*agrad_2_alpha[pp]);
// Dendro: reduced ops: 12
// Dendro: }}} 
//[[[end]]]
         
      }
     }
    }
    
    
    bssn::timer::t_rhs_a.stop();

    bssn::timer::t_rhs_b.start();
// b_rhs

    
 for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
    
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.b_rhs]
vnames = ['b_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 51 
// Dendro: printing temp variables
double DENDRO_0 = (3.0/4.0)*alpha[pp]*lambda_f[1] + (3.0/4.0)*lambda_f[0];

// Dendro: printing variables
//--
b_rhs0[pp] = B0[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta0[pp] + beta1[pp]*agrad_1_beta0[pp] + beta2[pp]*agrad_2_beta0[pp]);
//--
b_rhs1[pp] = B1[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta1[pp] + beta1[pp]*agrad_1_beta1[pp] + beta2[pp]*agrad_2_beta1[pp]);
//--
b_rhs2[pp] = B2[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta2[pp] + beta1[pp]*agrad_1_beta2[pp] + beta2[pp]*agrad_2_beta2[pp]);
// Dendro: reduced ops: 39
// Dendro: }}} 
//[[[end]]]
      }
    }
 }
    bssn::timer::t_rhs_b.stop();

    bssn::timer::t_rhs_gt.start();

//gt_rhs
    
 for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.gt_rhs]
vnames = ['gt_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 210 
// Dendro: printing temp variables
double DENDRO_0 = 2*alpha[pp];
double DENDRO_1 = grad_0_beta0[pp];
double DENDRO_2 = grad_1_beta1[pp];
double DENDRO_3 = (2.0/3.0)*gt0[pp];
double DENDRO_4 = grad_2_beta2[pp];
double DENDRO_5 = grad_0_beta1[pp];
double DENDRO_6 = 2*gt1[pp];
double DENDRO_7 = grad_0_beta2[pp];
double DENDRO_8 = 2*gt2[pp];
double DENDRO_9 = grad_1_beta0[pp];
double DENDRO_10 = grad_1_beta2[pp];
double DENDRO_11 = (1.0/3.0)*gt1[pp];
double DENDRO_12 = (2.0/3.0)*DENDRO_4;
double DENDRO_13 = grad_2_beta0[pp];
double DENDRO_14 = grad_2_beta1[pp];
double DENDRO_15 = (1.0/3.0)*gt2[pp];
double DENDRO_16 = (2.0/3.0)*DENDRO_2;
double DENDRO_17 = (2.0/3.0)*DENDRO_1;
double DENDRO_18 = 2*gt4[pp];
double DENDRO_19 = (1.0/3.0)*gt4[pp];

// Dendro: printing variables
//--
gt_rhs00[pp] = -At0[pp]*DENDRO_0 + (4.0/3.0)*DENDRO_1*gt0[pp] - DENDRO_2*DENDRO_3 - DENDRO_3*DENDRO_4 + DENDRO_5*DENDRO_6 + DENDRO_7*DENDRO_8 + beta0[pp]*agrad_0_gt0[pp] + beta1[pp]*agrad_1_gt0[pp] + beta2[pp]*agrad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp]*DENDRO_0 + DENDRO_1*DENDRO_11 + DENDRO_10*gt2[pp] + DENDRO_11*DENDRO_2 - DENDRO_12*gt1[pp] + DENDRO_5*gt3[pp] + DENDRO_7*gt4[pp] + DENDRO_9*gt0[pp] + beta0[pp]*agrad_0_gt1[pp] + beta1[pp]*agrad_1_gt1[pp] + beta2[pp]*agrad_2_gt1[pp];
//--
gt_rhs02[pp] = -At2[pp]*DENDRO_0 + DENDRO_1*DENDRO_15 + DENDRO_13*gt0[pp] + DENDRO_14*gt1[pp] + DENDRO_15*DENDRO_4 - DENDRO_16*gt2[pp] + DENDRO_5*gt4[pp] + DENDRO_7*gt5[pp] + beta0[pp]*agrad_0_gt2[pp] + beta1[pp]*agrad_1_gt2[pp] + beta2[pp]*agrad_2_gt2[pp];
//--
gt_rhs11[pp] = -At3[pp]*DENDRO_0 + DENDRO_10*DENDRO_18 - DENDRO_12*gt3[pp] - DENDRO_17*gt3[pp] + (4.0/3.0)*DENDRO_2*gt3[pp] + DENDRO_6*DENDRO_9 + beta0[pp]*agrad_0_gt3[pp] + beta1[pp]*agrad_1_gt3[pp] + beta2[pp]*agrad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp]*DENDRO_0 + DENDRO_10*gt5[pp] + DENDRO_13*gt1[pp] + DENDRO_14*gt3[pp] - DENDRO_17*gt4[pp] + DENDRO_19*DENDRO_2 + DENDRO_19*DENDRO_4 + DENDRO_9*gt2[pp] + beta0[pp]*agrad_0_gt4[pp] + beta1[pp]*agrad_1_gt4[pp] + beta2[pp]*agrad_2_gt4[pp];
//--
gt_rhs22[pp] = -At5[pp]*DENDRO_0 + DENDRO_13*DENDRO_8 + DENDRO_14*DENDRO_18 - DENDRO_16*gt5[pp] - DENDRO_17*gt5[pp] + (4.0/3.0)*DENDRO_4*gt5[pp] + beta0[pp]*agrad_0_gt5[pp] + beta1[pp]*agrad_1_gt5[pp] + beta2[pp]*agrad_2_gt5[pp];
// Dendro: reduced ops: 162
// Dendro: }}} 
//[[[end]]]
      }
    }
 }

    bssn::timer::t_rhs_gt.stop();

    bssn::timer::t_rhs_chi.start();
    // chi_rhs
    for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.chi_rhs]
vnames = ['chi_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 22 
// Dendro: printing temp variables
double DENDRO_0 = (2.0/3.0)*chi[pp];

// Dendro: printing variables
//--
chi_rhs[pp] = DENDRO_0*K[pp]*alpha[pp] - DENDRO_0*(grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp]) + beta0[pp]*agrad_0_chi[pp] + beta1[pp]*agrad_1_chi[pp] + beta2[pp]*agrad_2_chi[pp];
// Dendro: reduced ops: 20
// Dendro: }}} 
//[[[end]]]
      }
     }
    }

    bssn::timer::t_rhs_chi.stop();

    bssn::timer::t_rhs_At.start();
    // At_rhs
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.At_rhs]
vnames = ['At_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 630012 
// Dendro: printing temp variables
double DENDRO_0 = grad_0_beta0[pp];
double DENDRO_1 = grad_1_beta1[pp];
double DENDRO_2 = (2.0/3.0)*At0[pp];
double DENDRO_3 = grad_2_beta2[pp];
double DENDRO_4 = grad_0_beta1[pp];
double DENDRO_5 = 2*At1[pp];
double DENDRO_6 = grad_0_beta2[pp];
double DENDRO_7 = 2*At2[pp];
double DENDRO_8 = gt2[pp]*gt4[pp];
double DENDRO_9 = -DENDRO_8 + gt1[pp]*gt5[pp];
double DENDRO_10 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
double DENDRO_11 = gt0[pp]*gt5[pp];
double DENDRO_12 = pow(gt2[pp], 2);
double DENDRO_13 = DENDRO_11 - DENDRO_12;
double DENDRO_14 = pow(gt4[pp], 2);
double DENDRO_15 = pow(gt1[pp], 2);
double DENDRO_16 = gt3[pp]*gt5[pp];
double DENDRO_17 = DENDRO_12*gt3[pp] + DENDRO_14*gt0[pp] + DENDRO_15*gt5[pp] - DENDRO_16*gt0[pp] - 2*DENDRO_8*gt1[pp];
double DENDRO_18 = 1.0/DENDRO_17;
double DENDRO_19 = DENDRO_18*DENDRO_5;
double DENDRO_20 = At1[pp]*DENDRO_9;
double DENDRO_21 = gt1[pp]*gt4[pp];
double DENDRO_22 = gt2[pp]*gt3[pp];
double DENDRO_23 = DENDRO_21 - DENDRO_22;
double DENDRO_24 = -At2[pp]*DENDRO_23;
double DENDRO_25 = -DENDRO_14 + DENDRO_16;
double DENDRO_26 = 2*DENDRO_18;
double DENDRO_27 = At0[pp]*DENDRO_26;
double DENDRO_28 = gt0[pp]*gt3[pp];
double DENDRO_29 = -DENDRO_15 + DENDRO_28;
double DENDRO_30 = DENDRO_18*DENDRO_7;
double DENDRO_31 = grad2_0_0_alpha[pp];
double DENDRO_32 = 1.0/chi[pp];
double DENDRO_33 = grad_1_chi[pp];
double DENDRO_34 = grad_2_chi[pp];
double DENDRO_35 = grad_0_chi[pp];
double DENDRO_36 = DENDRO_10*DENDRO_34 + DENDRO_35*DENDRO_9;
double DENDRO_37 = -DENDRO_13*DENDRO_33 + DENDRO_36;
double DENDRO_38 = DENDRO_32*DENDRO_37;
double DENDRO_39 = 0.5*gt0[pp];
double DENDRO_40 = grad_0_gt1[pp];
double DENDRO_41 = 1.0*DENDRO_40;
double DENDRO_42 = grad_1_gt0[pp];
double DENDRO_43 = 0.5*DENDRO_42;
double DENDRO_44 = DENDRO_41 - DENDRO_43;
double DENDRO_45 = grad_0_gt2[pp];
double DENDRO_46 = 1.0*DENDRO_45;
double DENDRO_47 = grad_2_gt0[pp];
double DENDRO_48 = 0.5*DENDRO_47;
double DENDRO_49 = DENDRO_46 - DENDRO_48;
double DENDRO_50 = grad_0_gt0[pp];
double DENDRO_51 = 0.5*DENDRO_50;
double DENDRO_52 = DENDRO_10*DENDRO_49 + DENDRO_51*DENDRO_9;
double DENDRO_53 = -DENDRO_13*DENDRO_44 + DENDRO_52;
double DENDRO_54 = DENDRO_38*DENDRO_39 + DENDRO_53;
double DENDRO_55 = grad_1_alpha[pp];
double DENDRO_56 = 12*DENDRO_55;
double DENDRO_57 = DENDRO_18*DENDRO_56;
double DENDRO_58 = DENDRO_10*DENDRO_33;
double DENDRO_59 = DENDRO_23*DENDRO_35;
double DENDRO_60 = DENDRO_29*DENDRO_34;
double DENDRO_61 = DENDRO_32*(DENDRO_58 - DENDRO_59 - DENDRO_60);
double DENDRO_62 = DENDRO_39*DENDRO_61;
double DENDRO_63 = DENDRO_23*DENDRO_51;
double DENDRO_64 = DENDRO_29*DENDRO_49;
double DENDRO_65 = DENDRO_10*DENDRO_44;
double DENDRO_66 = DENDRO_63 + DENDRO_64 - DENDRO_65;
double DENDRO_67 = grad_2_alpha[pp];
double DENDRO_68 = 12*DENDRO_67;
double DENDRO_69 = DENDRO_18*DENDRO_68;
double DENDRO_70 = DENDRO_25*DENDRO_51;
double DENDRO_71 = DENDRO_18*DENDRO_70;
double DENDRO_72 = DENDRO_23*DENDRO_49;
double DENDRO_73 = DENDRO_18*DENDRO_72;
double DENDRO_74 = DENDRO_44*DENDRO_9;
double DENDRO_75 = DENDRO_18*DENDRO_74;
double DENDRO_76 = -DENDRO_21 + DENDRO_22;
double DENDRO_77 = DENDRO_33*DENDRO_9;
double DENDRO_78 = DENDRO_14 - DENDRO_16;
double DENDRO_79 = DENDRO_34*DENDRO_76 + DENDRO_35*DENDRO_78 + DENDRO_77;
double DENDRO_80 = DENDRO_18*gt0[pp];
double DENDRO_81 = DENDRO_32*(-1.0*DENDRO_35 + 0.5*DENDRO_79*DENDRO_80);
double DENDRO_82 = grad_0_alpha[pp];
double DENDRO_83 = 12*DENDRO_82;
double DENDRO_84 = grad2_0_0_chi[pp];
double DENDRO_85 = -DENDRO_84;
double DENDRO_86 = -DENDRO_63 - DENDRO_64 + DENDRO_65;
double DENDRO_87 = DENDRO_18*DENDRO_34;
double DENDRO_88 = DENDRO_18*DENDRO_33;
double DENDRO_89 = -DENDRO_70 - DENDRO_72 + DENDRO_74;
double DENDRO_90 = DENDRO_18*DENDRO_35;
double DENDRO_91 = 2*DENDRO_32;
double DENDRO_92 = grad_0_gt3[pp];
double DENDRO_93 = 0.5*DENDRO_92;
double DENDRO_94 = grad_1_gt1[pp];
double DENDRO_95 = 1.0*DENDRO_94;
double DENDRO_96 = -DENDRO_95;
double DENDRO_97 = DENDRO_93 + DENDRO_96;
double DENDRO_98 = grad_0_gt4[pp];
double DENDRO_99 = grad_2_gt1[pp];
double DENDRO_100 = grad_1_gt2[pp];
double DENDRO_101 = -DENDRO_100 + DENDRO_98 + DENDRO_99;
double DENDRO_102 = grad_0_gt5[pp];
double DENDRO_103 = DENDRO_10*DENDRO_102 + DENDRO_47*DENDRO_9;
double DENDRO_104 = -DENDRO_101*DENDRO_13 + DENDRO_103;
double DENDRO_105 = DENDRO_104*DENDRO_97;
double DENDRO_106 = DENDRO_13*DENDRO_92;
double DENDRO_107 = DENDRO_100 + DENDRO_98 - DENDRO_99;
double DENDRO_108 = DENDRO_10*DENDRO_107;
double DENDRO_109 = DENDRO_42*DENDRO_9;
double DENDRO_110 = DENDRO_108 + DENDRO_109;
double DENDRO_111 = -DENDRO_106 + DENDRO_110;
double DENDRO_112 = DENDRO_101*DENDRO_111;
double DENDRO_113 = 0.25*DENDRO_112;
double DENDRO_114 = pow(DENDRO_17, -2);
double DENDRO_115 = DENDRO_10*DENDRO_114;
double DENDRO_116 = 4*DENDRO_115;
double DENDRO_117 = 0.5*DENDRO_102;
double DENDRO_118 = grad_2_gt2[pp];
double DENDRO_119 = 1.0*DENDRO_118;
double DENDRO_120 = -DENDRO_119;
double DENDRO_121 = DENDRO_117 + DENDRO_120;
double DENDRO_122 = DENDRO_10*DENDRO_92;
double DENDRO_123 = DENDRO_23*DENDRO_42;
double DENDRO_124 = DENDRO_107*DENDRO_29;
double DENDRO_125 = DENDRO_122 - DENDRO_123 - DENDRO_124;
double DENDRO_126 = DENDRO_121*DENDRO_125;
double DENDRO_127 = DENDRO_23*DENDRO_47;
double DENDRO_128 = DENDRO_102*DENDRO_29;
double DENDRO_129 = DENDRO_10*DENDRO_101;
double DENDRO_130 = -DENDRO_127 - DENDRO_128 + DENDRO_129;
double DENDRO_131 = DENDRO_107*DENDRO_130;
double DENDRO_132 = 0.25*DENDRO_131;
double DENDRO_133 = -DENDRO_132;
double DENDRO_134 = -DENDRO_122 + DENDRO_123 + DENDRO_124;
double DENDRO_135 = 0.25*DENDRO_102;
double DENDRO_136 = DENDRO_134*DENDRO_135;
double DENDRO_137 = DENDRO_100 - DENDRO_98 + DENDRO_99;
double DENDRO_138 = DENDRO_127 + DENDRO_128 - DENDRO_129;
double DENDRO_139 = 0.5*DENDRO_138;
double DENDRO_140 = DENDRO_9*DENDRO_92;
double DENDRO_141 = DENDRO_25*DENDRO_42;
double DENDRO_142 = DENDRO_107*DENDRO_23;
double DENDRO_143 = -DENDRO_140 + DENDRO_141 + DENDRO_142;
double DENDRO_144 = 0.25*DENDRO_47;
double DENDRO_145 = DENDRO_143*DENDRO_144;
double DENDRO_146 = DENDRO_102*DENDRO_23;
double DENDRO_147 = DENDRO_25*DENDRO_47;
double DENDRO_148 = DENDRO_101*DENDRO_9;
double DENDRO_149 = DENDRO_146 + DENDRO_147 - DENDRO_148;
double DENDRO_150 = 0.25*DENDRO_42;
double DENDRO_151 = DENDRO_149*DENDRO_150;
double DENDRO_152 = 2*DENDRO_121;
double DENDRO_153 = DENDRO_114*DENDRO_23;
double DENDRO_154 = 4*DENDRO_153;
double DENDRO_155 = DENDRO_143*DENDRO_50;
double DENDRO_156 = 0.25*DENDRO_155;
double DENDRO_157 = DENDRO_70 + DENDRO_72 - DENDRO_74;
double DENDRO_158 = DENDRO_157*DENDRO_42;
double DENDRO_159 = DENDRO_114*DENDRO_9;
double DENDRO_160 = 4*DENDRO_159;
double DENDRO_161 = 0.5*DENDRO_49;
double DENDRO_162 = DENDRO_149*DENDRO_50;
double DENDRO_163 = 0.25*DENDRO_162;
double DENDRO_164 = DENDRO_157*DENDRO_47;
double DENDRO_165 = DENDRO_138*DENDRO_47;
double DENDRO_166 = DENDRO_102*DENDRO_66;
double DENDRO_167 = 2.0*DENDRO_153;
double DENDRO_168 = DENDRO_134*DENDRO_47;
double DENDRO_169 = DENDRO_107*DENDRO_66;
double DENDRO_170 = 2.0*DENDRO_159;
double DENDRO_171 = DENDRO_104*DENDRO_23;
double DENDRO_172 = grad_2_gt3[pp];
double DENDRO_173 = DENDRO_13*DENDRO_172;
double DENDRO_174 = grad_1_gt5[pp];
double DENDRO_175 = DENDRO_10*DENDRO_174;
double DENDRO_176 = DENDRO_137*DENDRO_9;
double DENDRO_177 = DENDRO_175 + DENDRO_176;
double DENDRO_178 = -DENDRO_173 + DENDRO_177;
double DENDRO_179 = DENDRO_10*DENDRO_178;
double DENDRO_180 = DENDRO_111*DENDRO_9;
double DENDRO_181 = grad_2_gt5[pp];
double DENDRO_182 = 0.5*DENDRO_181;
double DENDRO_183 = DENDRO_10*DENDRO_182;
double DENDRO_184 = 0.5*DENDRO_174;
double DENDRO_185 = grad_2_gt4[pp];
double DENDRO_186 = 1.0*DENDRO_185;
double DENDRO_187 = -DENDRO_186;
double DENDRO_188 = DENDRO_184 + DENDRO_187;
double DENDRO_189 = -DENDRO_121*DENDRO_9 + DENDRO_13*DENDRO_188 + DENDRO_183;
double DENDRO_190 = DENDRO_189*DENDRO_29;
double DENDRO_191 = grad_1_gt3[pp];
double DENDRO_192 = 0.5*DENDRO_191;
double DENDRO_193 = DENDRO_13*DENDRO_192;
double DENDRO_194 = grad_1_gt4[pp];
double DENDRO_195 = 1.0*DENDRO_194;
double DENDRO_196 = 0.5*DENDRO_172;
double DENDRO_197 = DENDRO_195 - DENDRO_196;
double DENDRO_198 = DENDRO_10*DENDRO_197;
double DENDRO_199 = DENDRO_9*DENDRO_97;
double DENDRO_200 = -DENDRO_193 + DENDRO_198 - DENDRO_199;
double DENDRO_201 = DENDRO_13*DENDRO_200;
double DENDRO_202 = DENDRO_25*DENDRO_53;
double DENDRO_203 = DENDRO_171 - 1.0*DENDRO_179 - 1.0*DENDRO_180 + DENDRO_190 + DENDRO_201 + DENDRO_202;
double DENDRO_204 = 2.0*DENDRO_114;
double DENDRO_205 = DENDRO_204*DENDRO_42;
double DENDRO_206 = -DENDRO_146 - DENDRO_147 + DENDRO_148;
double DENDRO_207 = DENDRO_206*DENDRO_23;
double DENDRO_208 = DENDRO_172*DENDRO_9;
double DENDRO_209 = DENDRO_174*DENDRO_23;
double DENDRO_210 = DENDRO_137*DENDRO_25;
double DENDRO_211 = DENDRO_208 - DENDRO_209 - DENDRO_210;
double DENDRO_212 = DENDRO_10*DENDRO_211;
double DENDRO_213 = DENDRO_140 - DENDRO_141 - DENDRO_142;
double DENDRO_214 = DENDRO_213*DENDRO_9;
double DENDRO_215 = DENDRO_182*DENDRO_23;
double DENDRO_216 = DENDRO_188*DENDRO_9;
double DENDRO_217 = DENDRO_121*DENDRO_25;
double DENDRO_218 = -DENDRO_215 - DENDRO_216 + DENDRO_217;
double DENDRO_219 = DENDRO_218*DENDRO_29;
double DENDRO_220 = DENDRO_192*DENDRO_9;
double DENDRO_221 = -DENDRO_197*DENDRO_23 + DENDRO_220 + DENDRO_25*DENDRO_97;
double DENDRO_222 = DENDRO_13*DENDRO_221;
double DENDRO_223 = DENDRO_25*DENDRO_89;
double DENDRO_224 = DENDRO_207 - 1.0*DENDRO_212 - 1.0*DENDRO_214 + DENDRO_219 + DENDRO_222 + DENDRO_223;
double DENDRO_225 = DENDRO_204*DENDRO_50;
double DENDRO_226 = DENDRO_130*DENDRO_23;
double DENDRO_227 = DENDRO_10*DENDRO_172;
double DENDRO_228 = DENDRO_174*DENDRO_29;
double DENDRO_229 = DENDRO_137*DENDRO_23;
double DENDRO_230 = DENDRO_227 - DENDRO_228 - DENDRO_229;
double DENDRO_231 = DENDRO_10*DENDRO_230;
double DENDRO_232 = DENDRO_125*DENDRO_9;
double DENDRO_233 = DENDRO_182*DENDRO_29;
double DENDRO_234 = DENDRO_121*DENDRO_23;
double DENDRO_235 = DENDRO_10*DENDRO_188;
double DENDRO_236 = -DENDRO_233 + DENDRO_234 - DENDRO_235;
double DENDRO_237 = DENDRO_236*DENDRO_29;
double DENDRO_238 = DENDRO_10*DENDRO_192;
double DENDRO_239 = -DENDRO_197*DENDRO_29 + DENDRO_23*DENDRO_97 + DENDRO_238;
double DENDRO_240 = DENDRO_13*DENDRO_239;
double DENDRO_241 = DENDRO_25*DENDRO_86;
double DENDRO_242 = DENDRO_226 - 1.0*DENDRO_231 - 1.0*DENDRO_232 + DENDRO_237 + DENDRO_240 + DENDRO_241;
double DENDRO_243 = DENDRO_204*DENDRO_47;
double DENDRO_244 = grad2_2_2_chi[pp];
double DENDRO_245 = pow(DENDRO_34, 2);
double DENDRO_246 = 3*DENDRO_32;
double DENDRO_247 = DENDRO_29*(2*DENDRO_244 - DENDRO_245*DENDRO_246);
double DENDRO_248 = grad2_1_1_chi[pp];
double DENDRO_249 = pow(DENDRO_33, 2);
double DENDRO_250 = DENDRO_13*(-DENDRO_246*DENDRO_249 + 2*DENDRO_248);
double DENDRO_251 = pow(DENDRO_35, 2);
double DENDRO_252 = DENDRO_25*(-DENDRO_246*DENDRO_251 + 2*DENDRO_84);
double DENDRO_253 = grad2_1_2_chi[pp];
double DENDRO_254 = DENDRO_246*DENDRO_34;
double DENDRO_255 = 2*DENDRO_10*(2*DENDRO_253 - DENDRO_254*DENDRO_33);
double DENDRO_256 = grad2_0_2_chi[pp];
double DENDRO_257 = 2*DENDRO_23*(-DENDRO_254*DENDRO_35 + 2*DENDRO_256);
double DENDRO_258 = grad2_0_1_chi[pp];
double DENDRO_259 = DENDRO_33*DENDRO_35;
double DENDRO_260 = 2*DENDRO_9*(-DENDRO_246*DENDRO_259 + 2*DENDRO_258);
double DENDRO_261 = -1.0*DENDRO_171 + DENDRO_179 + DENDRO_180 - DENDRO_190 - DENDRO_201 - DENDRO_202;
double DENDRO_262 = 2*DENDRO_261*DENDRO_88;
double DENDRO_263 = -1.0*DENDRO_226 + DENDRO_231 + DENDRO_232 - DENDRO_237 - DENDRO_240 - DENDRO_241;
double DENDRO_264 = 2*DENDRO_263*DENDRO_87;
double DENDRO_265 = -1.0*DENDRO_207 + DENDRO_212 + DENDRO_214 - DENDRO_219 - DENDRO_222 - DENDRO_223;
double DENDRO_266 = 2*DENDRO_265*DENDRO_90;
double DENDRO_267 = DENDRO_247 + DENDRO_250 + DENDRO_252 - DENDRO_255 + DENDRO_257 - DENDRO_260 + DENDRO_262 + DENDRO_264 + DENDRO_266;
double DENDRO_268 = DENDRO_32*DENDRO_80;
double DENDRO_269 = DENDRO_114*DENDRO_29;
double DENDRO_270 = 4*DENDRO_269;
double DENDRO_271 = DENDRO_130*DENDRO_270;
double DENDRO_272 = 0.25*DENDRO_92;
double DENDRO_273 = DENDRO_114*DENDRO_13;
double DENDRO_274 = 4*DENDRO_273;
double DENDRO_275 = DENDRO_111*DENDRO_274;
double DENDRO_276 = 0.25*DENDRO_98;
double DENDRO_277 = -DENDRO_276;
double DENDRO_278 = 0.75*DENDRO_100;
double DENDRO_279 = 0.25*DENDRO_99;
double DENDRO_280 = DENDRO_274*(DENDRO_277 + DENDRO_278 + DENDRO_279);
double DENDRO_281 = DENDRO_46 + DENDRO_48;
double DENDRO_282 = DENDRO_114*DENDRO_25;
double DENDRO_283 = 4*DENDRO_282;
double DENDRO_284 = DENDRO_149*DENDRO_47;
double DENDRO_285 = 3.0*DENDRO_269;
double DENDRO_286 = DENDRO_143*DENDRO_42;
double DENDRO_287 = 3.0*DENDRO_273;
double DENDRO_288 = 6.0*DENDRO_50;
double DENDRO_289 = pow(chi[pp], -2);
double DENDRO_290 = grad_0_Gt0[pp];
double DENDRO_291 = grad_0_Gt1[pp];
double DENDRO_292 = 4*gt1[pp];
double DENDRO_293 = grad_0_Gt2[pp];
double DENDRO_294 = 4*gt2[pp];
double DENDRO_295 = 0.5*DENDRO_44;
double DENDRO_296 = DENDRO_104*DENDRO_295;
double DENDRO_297 = 4*DENDRO_18*DENDRO_23;
double DENDRO_298 = DENDRO_104*DENDRO_272;
double DENDRO_299 = 0.5*DENDRO_137;
double DENDRO_300 = DENDRO_111*DENDRO_295;
double DENDRO_301 = 2.0*DENDRO_18;
double DENDRO_302 = DENDRO_29*DENDRO_301;
double DENDRO_303 = DENDRO_13*DENDRO_301;
double DENDRO_304 = DENDRO_25*DENDRO_301;
double DENDRO_305 = DENDRO_111*DENDRO_42;
double DENDRO_306 = DENDRO_53*DENDRO_92;
double DENDRO_307 = DENDRO_104*DENDRO_42;
double DENDRO_308 = DENDRO_101*DENDRO_53;
double DENDRO_309 = 4.0*DENDRO_18;
double DENDRO_310 = DENDRO_10*DENDRO_309;
double DENDRO_311 = DENDRO_309*DENDRO_9;
double DENDRO_312 = 0.25*DENDRO_100;
double DENDRO_313 = 0.75*DENDRO_99;
double DENDRO_314 = 4*DENDRO_114;
double DENDRO_315 = -DENDRO_104*DENDRO_270*(DENDRO_277 + DENDRO_312 + DENDRO_313) + DENDRO_116*(DENDRO_111*DENDRO_299 + DENDRO_298) - DENDRO_154*(DENDRO_137*DENDRO_53 + DENDRO_296) + DENDRO_160*(DENDRO_300 - 2*DENDRO_53*DENDRO_97) - DENDRO_167*(DENDRO_307 + DENDRO_308) + DENDRO_170*(DENDRO_305 + DENDRO_306) - DENDRO_202*DENDRO_314*(DENDRO_41 + DENDRO_43) - DENDRO_251*DENDRO_289 + 4*DENDRO_290*gt0[pp] + DENDRO_291*DENDRO_292 + DENDRO_293*DENDRO_294 + DENDRO_297*grad2_0_2_gt0[pp] + DENDRO_302*grad2_2_2_gt0[pp] + DENDRO_303*grad2_1_1_gt0[pp] + DENDRO_304*grad2_0_0_gt0[pp] - DENDRO_310*grad2_1_2_gt0[pp] - DENDRO_311*grad2_0_1_gt0[pp];
double DENDRO_316 = 3*alpha[pp];
double DENDRO_317 = grad2_2_2_alpha[pp];
double DENDRO_318 = 0.5*gt5[pp];
double DENDRO_319 = DENDRO_189 + DENDRO_318*DENDRO_38;
double DENDRO_320 = 4*DENDRO_55;
double DENDRO_321 = DENDRO_18*DENDRO_320;
double DENDRO_322 = DENDRO_23*DENDRO_34;
double DENDRO_323 = DENDRO_25*DENDRO_35;
double DENDRO_324 = DENDRO_32*(-DENDRO_322 - DENDRO_323 + DENDRO_77);
double DENDRO_325 = DENDRO_318*DENDRO_324;
double DENDRO_326 = 4*DENDRO_82;
double DENDRO_327 = DENDRO_18*DENDRO_326;
double DENDRO_328 = DENDRO_18*DENDRO_233;
double DENDRO_329 = DENDRO_18*DENDRO_234;
double DENDRO_330 = DENDRO_18*DENDRO_235;
double DENDRO_331 = DENDRO_15 - DENDRO_28;
double DENDRO_332 = DENDRO_331*DENDRO_34 + DENDRO_35*DENDRO_76 + DENDRO_58;
double DENDRO_333 = DENDRO_18*gt5[pp];
double DENDRO_334 = DENDRO_32*(0.5*DENDRO_332*DENDRO_333 - 1.0*DENDRO_34);
double DENDRO_335 = 4*DENDRO_67;
double DENDRO_336 = -DENDRO_244;
double DENDRO_337 = -DENDRO_184;
double DENDRO_338 = DENDRO_186 + DENDRO_337;
double DENDRO_339 = -DENDRO_117;
double DENDRO_340 = DENDRO_119 + DENDRO_339;
double DENDRO_341 = -DENDRO_11 + DENDRO_12;
double DENDRO_342 = 0.5*DENDRO_188;
double DENDRO_343 = DENDRO_104*DENDRO_342;
double DENDRO_344 = -DENDRO_343;
double DENDRO_345 = DENDRO_107*DENDRO_189;
double DENDRO_346 = 0.5*DENDRO_121;
double DENDRO_347 = -DENDRO_206*DENDRO_346;
double DENDRO_348 = 2*DENDRO_49;
double DENDRO_349 = DENDRO_130*DENDRO_181;
double DENDRO_350 = 0.25*DENDRO_349;
double DENDRO_351 = DENDRO_102*DENDRO_236;
double DENDRO_352 = DENDRO_211*DENDRO_346;
double DENDRO_353 = -DENDRO_352;
double DENDRO_354 = DENDRO_107*DENDRO_218;
double DENDRO_355 = DENDRO_181*DENDRO_230;
double DENDRO_356 = 0.25*DENDRO_355;
double DENDRO_357 = DENDRO_174*DENDRO_236;
double DENDRO_358 = DENDRO_211*DENDRO_49;
double DENDRO_359 = DENDRO_137*DENDRO_206;
double DENDRO_360 = 0.25*DENDRO_359;
double DENDRO_361 = 0.25*DENDRO_174;
double DENDRO_362 = DENDRO_130*DENDRO_361;
double DENDRO_363 = DENDRO_135*DENDRO_230;
double DENDRO_364 = DENDRO_144*DENDRO_211;
double DENDRO_365 = 0.5*DENDRO_107;
double DENDRO_366 = DENDRO_137*DENDRO_218;
double DENDRO_367 = 2.0*DENDRO_115;
double DENDRO_368 = DENDRO_204*DENDRO_261;
double DENDRO_369 = DENDRO_204*DENDRO_263;
double DENDRO_370 = DENDRO_204*DENDRO_265;
double DENDRO_371 = DENDRO_218*DENDRO_47;
double DENDRO_372 = -DENDRO_247 - DENDRO_250 - DENDRO_252 + DENDRO_255 - DENDRO_257 + DENDRO_260 - DENDRO_262 - DENDRO_264 - DENDRO_266;
double DENDRO_373 = DENDRO_32*DENDRO_372;
double DENDRO_374 = DENDRO_219*DENDRO_314;
double DENDRO_375 = DENDRO_190*DENDRO_314;
double DENDRO_376 = -DENDRO_279;
double DENDRO_377 = DENDRO_274*(DENDRO_276 + DENDRO_278 + DENDRO_376);
double DENDRO_378 = DENDRO_283*(-DENDRO_144 + DENDRO_46);
double DENDRO_379 = DENDRO_174*DENDRO_230;
double DENDRO_380 = DENDRO_102*DENDRO_130;
double DENDRO_381 = 3.0*DENDRO_282;
double DENDRO_382 = 6.0*DENDRO_114;
double DENDRO_383 = grad_2_Gt0[pp];
double DENDRO_384 = grad_2_Gt1[pp];
double DENDRO_385 = 4*gt4[pp];
double DENDRO_386 = grad_2_Gt2[pp];
double DENDRO_387 = -DENDRO_178*DENDRO_342;
double DENDRO_388 = DENDRO_104*DENDRO_197;
double DENDRO_389 = DENDRO_101*DENDRO_178;
double DENDRO_390 = 0.25*DENDRO_389;
double DENDRO_391 = 0.25*DENDRO_172;
double DENDRO_392 = DENDRO_104*DENDRO_391;
double DENDRO_393 = DENDRO_174*DENDRO_178;
double DENDRO_394 = DENDRO_172*DENDRO_189;
double DENDRO_395 = DENDRO_104*DENDRO_174;
double DENDRO_396 = DENDRO_101*DENDRO_189;
double DENDRO_397 = 0.75*DENDRO_98;
double DENDRO_398 = -DENDRO_104*DENDRO_283*(DENDRO_312 + DENDRO_376 + DENDRO_397) + DENDRO_116*(2*DENDRO_189*DENDRO_197 + DENDRO_387) + DENDRO_160*(DENDRO_388 + DENDRO_390) + DENDRO_160*(DENDRO_178*DENDRO_365 + DENDRO_392) - DENDRO_167*(DENDRO_395 + DENDRO_396) - DENDRO_178*DENDRO_274*(DENDRO_195 - DENDRO_391) - DENDRO_245*DENDRO_289 + DENDRO_294*DENDRO_383 + DENDRO_297*grad2_0_2_gt5[pp] + DENDRO_302*grad2_2_2_gt5[pp] + DENDRO_303*grad2_1_1_gt5[pp] + DENDRO_304*grad2_0_0_gt5[pp] - DENDRO_310*grad2_1_2_gt5[pp] - DENDRO_311*grad2_0_1_gt5[pp] + DENDRO_367*(DENDRO_393 + DENDRO_394) + DENDRO_384*DENDRO_385 + 4*DENDRO_386*gt5[pp];
double DENDRO_399 = grad2_1_1_alpha[pp];
double DENDRO_400 = 0.5*gt3[pp];
double DENDRO_401 = DENDRO_18*DENDRO_335;
double DENDRO_402 = -1.0*DENDRO_33;
double DENDRO_403 = DENDRO_33*DENDRO_341 + DENDRO_36;
double DENDRO_404 = DENDRO_18*gt3[pp];
double DENDRO_405 = 0.5*DENDRO_404;
double DENDRO_406 = -DENDRO_18*DENDRO_193 + DENDRO_18*DENDRO_198 - DENDRO_18*DENDRO_199;
double DENDRO_407 = -DENDRO_248;
double DENDRO_408 = -DENDRO_93;
double DENDRO_409 = DENDRO_408 + DENDRO_95;
double DENDRO_410 = DENDRO_211*DENDRO_44;
double DENDRO_411 = DENDRO_137*DENDRO_213;
double DENDRO_412 = 0.25*DENDRO_411;
double DENDRO_413 = DENDRO_125*DENDRO_188;
double DENDRO_414 = DENDRO_107*DENDRO_230;
double DENDRO_415 = 0.25*DENDRO_414;
double DENDRO_416 = DENDRO_125*DENDRO_361;
double DENDRO_417 = 0.5*DENDRO_101;
double DENDRO_418 = DENDRO_150*DENDRO_211;
double DENDRO_419 = 0.5*DENDRO_97;
double DENDRO_420 = DENDRO_211*DENDRO_419;
double DENDRO_421 = -DENDRO_420;
double DENDRO_422 = DENDRO_101*DENDRO_221;
double DENDRO_423 = 0.5*DENDRO_197;
double DENDRO_424 = DENDRO_230*DENDRO_423;
double DENDRO_425 = 2*DENDRO_188*DENDRO_239;
double DENDRO_426 = DENDRO_178*DENDRO_191;
double DENDRO_427 = 0.25*DENDRO_426;
double DENDRO_428 = DENDRO_172*DENDRO_200;
double DENDRO_429 = DENDRO_125*DENDRO_423;
double DENDRO_430 = DENDRO_101*DENDRO_239;
double DENDRO_431 = -DENDRO_213*DENDRO_419;
double DENDRO_432 = 2*DENDRO_221*DENDRO_44;
double DENDRO_433 = DENDRO_111*DENDRO_191;
double DENDRO_434 = 0.25*DENDRO_433;
double DENDRO_435 = DENDRO_200*DENDRO_92;
double DENDRO_436 = DENDRO_137*DENDRO_221;
double DENDRO_437 = DENDRO_174*DENDRO_239;
double DENDRO_438 = DENDRO_107*DENDRO_239;
double DENDRO_439 = DENDRO_221*DENDRO_42;
double DENDRO_440 = DENDRO_230*DENDRO_270;
double DENDRO_441 = -DENDRO_312;
double DENDRO_442 = DENDRO_270*(DENDRO_276 + DENDRO_313 + DENDRO_441);
double DENDRO_443 = DENDRO_222*DENDRO_314;
double DENDRO_444 = DENDRO_283*(-DENDRO_150 + DENDRO_41);
double DENDRO_445 = DENDRO_283*(DENDRO_279 + DENDRO_397 + DENDRO_441);
double DENDRO_446 = grad_1_Gt0[pp];
double DENDRO_447 = grad_1_Gt1[pp];
double DENDRO_448 = grad_1_Gt2[pp];
double DENDRO_449 = DENDRO_111*DENDRO_391;
double DENDRO_450 = DENDRO_178*DENDRO_272;
double DENDRO_451 = DENDRO_172*DENDRO_178;
double DENDRO_452 = DENDRO_111*DENDRO_92;
double DENDRO_453 = -DENDRO_154*(DENDRO_111*DENDRO_196 + DENDRO_450) - DENDRO_154*(DENDRO_178*DENDRO_93 + DENDRO_449) - DENDRO_240*DENDRO_314*(DENDRO_195 + DENDRO_196) - DENDRO_249*DENDRO_289 - DENDRO_285*DENDRO_451 + DENDRO_292*DENDRO_446 + DENDRO_297*grad2_0_2_gt3[pp] + DENDRO_302*grad2_2_2_gt3[pp] + DENDRO_303*grad2_1_1_gt3[pp] + DENDRO_304*grad2_0_0_gt3[pp] - DENDRO_310*grad2_1_2_gt3[pp] - DENDRO_311*grad2_0_1_gt3[pp] - DENDRO_381*DENDRO_452 + DENDRO_385*DENDRO_448 + 4*DENDRO_447*gt3[pp];
double DENDRO_454 = DENDRO_130*DENDRO_161;
double DENDRO_455 = DENDRO_206*DENDRO_50;
double DENDRO_456 = 0.25*DENDRO_455;
double DENDRO_457 = DENDRO_47*DENDRO_89;
double DENDRO_458 = DENDRO_125*DENDRO_135;
double DENDRO_459 = DENDRO_150*DENDRO_206;
double DENDRO_460 = DENDRO_144*DENDRO_213;
double DENDRO_461 = DENDRO_125*DENDRO_161;
double DENDRO_462 = DENDRO_213*DENDRO_50;
double DENDRO_463 = 0.25*DENDRO_462;
double DENDRO_464 = DENDRO_42*DENDRO_89;
double DENDRO_465 = DENDRO_107*DENDRO_86;
double DENDRO_466 = DENDRO_102*DENDRO_86;
double DENDRO_467 = DENDRO_206*DENDRO_47;
double DENDRO_468 = DENDRO_213*DENDRO_42;
double DENDRO_469 = DENDRO_230*DENDRO_47;
double DENDRO_470 = DENDRO_102*DENDRO_125;
double DENDRO_471 = DENDRO_131 + DENDRO_470;
double DENDRO_472 = DENDRO_178*DENDRO_42;
double DENDRO_473 = DENDRO_107*DENDRO_111;
double DENDRO_474 = DENDRO_104*DENDRO_92;
double DENDRO_475 = DENDRO_473 + DENDRO_474;
double DENDRO_476 = DENDRO_135*DENDRO_206;
double DENDRO_477 = DENDRO_121*DENDRO_236;
double DENDRO_478 = -DENDRO_477;
double DENDRO_479 = DENDRO_189*DENDRO_299 + 0.25*DENDRO_395;
double DENDRO_480 = DENDRO_107*DENDRO_213;
double DENDRO_481 = 0.25*DENDRO_480;
double DENDRO_482 = DENDRO_111*DENDRO_423;
double DENDRO_483 = -DENDRO_178*DENDRO_419;
double DENDRO_484 = DENDRO_482 + DENDRO_483;
double DENDRO_485 = DENDRO_49*DENDRO_89;
double DENDRO_486 = DENDRO_130*DENDRO_144;
double DENDRO_487 = 0.25*DENDRO_307 + DENDRO_365*DENDRO_53;
double DENDRO_488 = 0.25*DENDRO_104;
double DENDRO_489 = DENDRO_137*DENDRO_488 + DENDRO_184*DENDRO_53;
double DENDRO_490 = DENDRO_182*DENDRO_86;
double DENDRO_491 = -DENDRO_130*DENDRO_346;
double DENDRO_492 = DENDRO_101*DENDRO_488;
double DENDRO_493 = DENDRO_107*DENDRO_488 + DENDRO_189*DENDRO_43;
double DENDRO_494 = DENDRO_144*DENDRO_206;
double DENDRO_495 = DENDRO_218*DENDRO_51;
double DENDRO_496 = DENDRO_161*DENDRO_206 + DENDRO_495;
double DENDRO_497 = DENDRO_189*DENDRO_97;
double DENDRO_498 = 0.5*DENDRO_388;
double DENDRO_499 = -DENDRO_497 + DENDRO_498;
double DENDRO_500 = DENDRO_137*DENDRO_178;
double DENDRO_501 = 0.25*DENDRO_500;
double DENDRO_502 = DENDRO_189*DENDRO_93;
double DENDRO_503 = DENDRO_111*DENDRO_361 + DENDRO_502;
double DENDRO_504 = 0.25*DENDRO_107;
double DENDRO_505 = DENDRO_218*DENDRO_43;
double DENDRO_506 = DENDRO_206*DENDRO_504 + DENDRO_505;
double DENDRO_507 = DENDRO_236*DENDRO_299;
double DENDRO_508 = DENDRO_362 + DENDRO_363;
double DENDRO_509 = DENDRO_230*DENDRO_346;
double DENDRO_510 = -DENDRO_509;
double DENDRO_511 = 0.25*DENDRO_125;
double DENDRO_512 = DENDRO_181*DENDRO_511;
double DENDRO_513 = DENDRO_236*DENDRO_365 + DENDRO_512;
double DENDRO_514 = DENDRO_135*DENDRO_213 + DENDRO_505;
double DENDRO_515 = DENDRO_178*DENDRO_295;
double DENDRO_516 = -0.5*DENDRO_105 + DENDRO_197*DENDRO_53;
double DENDRO_517 = DENDRO_184*DENDRO_86;
double DENDRO_518 = DENDRO_161*DENDRO_230;
double DENDRO_519 = 0.25*DENDRO_137;
double DENDRO_520 = DENDRO_130*DENDRO_519;
double DENDRO_521 = 0.25*DENDRO_211;
double DENDRO_522 = DENDRO_50*DENDRO_521;
double DENDRO_523 = DENDRO_161*DENDRO_213;
double DENDRO_524 = DENDRO_522 + DENDRO_523;
double DENDRO_525 = DENDRO_365*DENDRO_89 + DENDRO_459;
double DENDRO_526 = DENDRO_137*DENDRO_230;
double DENDRO_527 = DENDRO_125*DENDRO_174;
double DENDRO_528 = DENDRO_414 + DENDRO_527;
double DENDRO_529 = 1.0*DENDRO_273;
double DENDRO_530 = -DENDRO_256;
double DENDRO_531 = 0.5*DENDRO_87;
double DENDRO_532 = 0.5*DENDRO_88;
double DENDRO_533 = 0.5*DENDRO_90;
double DENDRO_534 = 2.0*DENDRO_383;
double DENDRO_535 = DENDRO_534*gt0[pp];
double DENDRO_536 = 2.0*DENDRO_384;
double DENDRO_537 = DENDRO_536*gt1[pp];
double DENDRO_538 = 2.0*gt2[pp];
double DENDRO_539 = DENDRO_290*DENDRO_538;
double DENDRO_540 = DENDRO_386*DENDRO_538;
double DENDRO_541 = 2.0*gt4[pp];
double DENDRO_542 = DENDRO_291*DENDRO_541;
double DENDRO_543 = 2.0*gt5[pp];
double DENDRO_544 = DENDRO_293*DENDRO_543;
double DENDRO_545 = DENDRO_289*DENDRO_34;
double DENDRO_546 = -DENDRO_35*DENDRO_545;
double DENDRO_547 = DENDRO_297*grad2_0_2_gt2[pp];
double DENDRO_548 = DENDRO_302*grad2_2_2_gt2[pp];
double DENDRO_549 = DENDRO_303*grad2_1_1_gt2[pp];
double DENDRO_550 = DENDRO_304*grad2_0_0_gt2[pp];
double DENDRO_551 = -DENDRO_310*grad2_1_2_gt2[pp];
double DENDRO_552 = -DENDRO_311*grad2_0_1_gt2[pp];
double DENDRO_553 = DENDRO_18*gt2[pp];
double DENDRO_554 = DENDRO_100*DENDRO_368 + DENDRO_118*DENDRO_369 + DENDRO_370*DENDRO_45 + DENDRO_373*DENDRO_553 + DENDRO_535 + DENDRO_537 + DENDRO_539 + DENDRO_540 + DENDRO_542 + DENDRO_544 + DENDRO_546 + DENDRO_547 + DENDRO_548 + DENDRO_549 + DENDRO_550 + DENDRO_551 + DENDRO_552 - DENDRO_91*(DENDRO_530 + DENDRO_531*(DENDRO_102*DENDRO_331 + DENDRO_129 + DENDRO_47*DENDRO_76) + DENDRO_532*(DENDRO_101*DENDRO_341 + DENDRO_103) + DENDRO_533*(DENDRO_102*DENDRO_76 + DENDRO_148 + DENDRO_47*DENDRO_78));
double DENDRO_555 = grad2_0_2_alpha[pp];
double DENDRO_556 = DENDRO_127*DENDRO_18;
double DENDRO_557 = DENDRO_128*DENDRO_18;
double DENDRO_558 = DENDRO_129*DENDRO_18;
double DENDRO_559 = -DENDRO_35;
double DENDRO_560 = DENDRO_32*(DENDRO_332*DENDRO_553 + DENDRO_559);
double DENDRO_561 = 2.0*DENDRO_67;
double DENDRO_562 = DENDRO_146*DENDRO_18;
double DENDRO_563 = DENDRO_147*DENDRO_18;
double DENDRO_564 = DENDRO_148*DENDRO_18;
double DENDRO_565 = -DENDRO_34;
double DENDRO_566 = DENDRO_32*(DENDRO_553*DENDRO_79 + DENDRO_565);
double DENDRO_567 = 2.0*DENDRO_82;
double DENDRO_568 = 2.0*DENDRO_55;
double DENDRO_569 = DENDRO_18*(DENDRO_104 + DENDRO_38*gt2[pp]);
double DENDRO_570 = -4*DENDRO_555 + DENDRO_561*(-DENDRO_556 - DENDRO_557 + DENDRO_558 + DENDRO_560) + DENDRO_567*(-DENDRO_562 - DENDRO_563 + DENDRO_564 + DENDRO_566) + DENDRO_568*DENDRO_569;
double DENDRO_571 = DENDRO_211*DENDRO_47;
double DENDRO_572 = DENDRO_102*DENDRO_213 + DENDRO_571;
double DENDRO_573 = 0.5*DENDRO_349;
double DENDRO_574 = 0.25*DENDRO_526;
double DENDRO_575 = DENDRO_130*DENDRO_135 + DENDRO_490;
double DENDRO_576 = DENDRO_111*DENDRO_342;
double DENDRO_577 = -DENDRO_576;
double DENDRO_578 = -DENDRO_213*DENDRO_346;
double DENDRO_579 = DENDRO_362 + DENDRO_512;
double DENDRO_580 = DENDRO_458 + DENDRO_517;
double DENDRO_581 = DENDRO_144*DENDRO_230;
double DENDRO_582 = DENDRO_299*DENDRO_89;
double DENDRO_583 = DENDRO_211*DENDRO_42;
double DENDRO_584 = DENDRO_480 + DENDRO_583;
double DENDRO_585 = DENDRO_104*DENDRO_172;
double DENDRO_586 = DENDRO_111*DENDRO_174 + DENDRO_585;
double DENDRO_587 = 0.25*DENDRO_473;
double DENDRO_588 = DENDRO_196*DENDRO_53;
double DENDRO_589 = DENDRO_150*DENDRO_178 + DENDRO_588;
double DENDRO_590 = DENDRO_115*(DENDRO_500 + DENDRO_586) - DENDRO_154*(DENDRO_489 + DENDRO_492) - DENDRO_154*(-DENDRO_188*DENDRO_53 + DENDRO_493) + DENDRO_160*(DENDRO_113 + DENDRO_516) + DENDRO_160*(DENDRO_587 + DENDRO_589) - DENDRO_270*(DENDRO_344 + DENDRO_479) - DENDRO_274*(DENDRO_449 + DENDRO_484) - DENDRO_283*(0.5*DENDRO_308 + DENDRO_487);
double DENDRO_591 = DENDRO_130*DENDRO_172;
double DENDRO_592 = DENDRO_206*DENDRO_92;
double DENDRO_593 = 0.25*DENDRO_393;
double DENDRO_594 = DENDRO_188*DENDRO_236;
double DENDRO_595 = -DENDRO_594;
double DENDRO_596 = DENDRO_135*DENDRO_211 + DENDRO_218*DENDRO_417;
double DENDRO_597 = DENDRO_197*DENDRO_200;
double DENDRO_598 = DENDRO_230*DENDRO_391;
double DENDRO_599 = DENDRO_221*DENDRO_365;
double DENDRO_600 = DENDRO_211*DENDRO_272 + DENDRO_599;
double DENDRO_601 = DENDRO_206*DENDRO_295;
double DENDRO_602 = DENDRO_523 + DENDRO_601;
double DENDRO_603 = DENDRO_218*DENDRO_44 + 0.5*DENDRO_358;
double DENDRO_604 = DENDRO_101*DENDRO_206;
double DENDRO_605 = 0.25*DENDRO_604;
double DENDRO_606 = DENDRO_178*DENDRO_504 + DENDRO_502;
double DENDRO_607 = DENDRO_236*DENDRO_417;
double DENDRO_608 = DENDRO_130*DENDRO_342;
double DENDRO_609 = -DENDRO_608;
double DENDRO_610 = DENDRO_182*DENDRO_239;
double DENDRO_611 = -DENDRO_230*DENDRO_342 + DENDRO_610;
double DENDRO_612 = DENDRO_117*DENDRO_221;
double DENDRO_613 = DENDRO_101*DENDRO_521 + DENDRO_612;
double DENDRO_614 = DENDRO_211*DENDRO_519;
double DENDRO_615 = DENDRO_211*DENDRO_504 + DENDRO_218*DENDRO_93;
double DENDRO_616 = DENDRO_178*DENDRO_423;
double DENDRO_617 = DENDRO_189*DENDRO_192;
double DENDRO_618 = DENDRO_178*DENDRO_391 + DENDRO_617;
double DENDRO_619 = -DENDRO_206*DENDRO_419;
double DENDRO_620 = DENDRO_221*DENDRO_49;
double DENDRO_621 = 0.5*DENDRO_410 + DENDRO_620;
double DENDRO_622 = DENDRO_130*DENDRO_423;
double DENDRO_623 = 0.25*DENDRO_101;
double DENDRO_624 = DENDRO_117*DENDRO_239;
double DENDRO_625 = DENDRO_230*DENDRO_623 + DENDRO_624;
double DENDRO_626 = DENDRO_191*DENDRO_488;
double DENDRO_627 = DENDRO_450 + DENDRO_626;
double DENDRO_628 = DENDRO_200*DENDRO_365;
double DENDRO_629 = DENDRO_101*DENDRO_130;
double DENDRO_630 = 1.0*DENDRO_282;
double DENDRO_631 = -DENDRO_253;
double DENDRO_632 = DENDRO_18*gt4[pp];
double DENDRO_633 = DENDRO_297*grad2_0_2_gt4[pp] + DENDRO_302*grad2_2_2_gt4[pp] + DENDRO_303*grad2_1_1_gt4[pp] + DENDRO_304*grad2_0_0_gt4[pp] - DENDRO_310*grad2_1_2_gt4[pp] - DENDRO_311*grad2_0_1_gt4[pp] - DENDRO_33*DENDRO_545 + DENDRO_386*DENDRO_541 + DENDRO_446*DENDRO_538 + DENDRO_447*DENDRO_541 + DENDRO_448*DENDRO_543 + DENDRO_534*gt1[pp] + DENDRO_536*gt3[pp];
double DENDRO_634 = DENDRO_185*DENDRO_369 + DENDRO_194*DENDRO_368 + DENDRO_370*DENDRO_98 + DENDRO_373*DENDRO_632 + DENDRO_633 - DENDRO_91*(DENDRO_531*(DENDRO_137*DENDRO_76 + DENDRO_174*DENDRO_331 + DENDRO_227) + DENDRO_532*(DENDRO_172*DENDRO_341 + DENDRO_177) + DENDRO_533*(DENDRO_137*DENDRO_78 + DENDRO_174*DENDRO_76 + DENDRO_208) + DENDRO_631);
double DENDRO_635 = grad2_1_2_alpha[pp];
double DENDRO_636 = DENDRO_18*DENDRO_227;
double DENDRO_637 = DENDRO_18*DENDRO_228;
double DENDRO_638 = DENDRO_18*DENDRO_229;
double DENDRO_639 = -DENDRO_33;
double DENDRO_640 = DENDRO_32*(DENDRO_332*DENDRO_632 + DENDRO_639);
double DENDRO_641 = -DENDRO_173*DENDRO_18 + DENDRO_175*DENDRO_18 + DENDRO_176*DENDRO_18;
double DENDRO_642 = DENDRO_324*gt4[pp];
double DENDRO_643 = DENDRO_18*DENDRO_567*(DENDRO_211 + DENDRO_642) + DENDRO_561*(DENDRO_636 - DENDRO_637 - DENDRO_638 + DENDRO_640) + DENDRO_568*(DENDRO_32*(DENDRO_403*DENDRO_632 + DENDRO_565) + DENDRO_641) - 4*DENDRO_635;
double DENDRO_644 = 0.5*DENDRO_355;
double DENDRO_645 = 1.0*DENDRO_437;
double DENDRO_646 = 0.5*DENDRO_436;
double DENDRO_647 = 0.25*DENDRO_629;
double DENDRO_648 = DENDRO_363 + DENDRO_512;
double DENDRO_649 = DENDRO_121*DENDRO_221;
double DENDRO_650 = DENDRO_616 + DENDRO_617;
double DENDRO_651 = DENDRO_230*DENDRO_361;
double DENDRO_652 = DENDRO_221*DENDRO_48;
double DENDRO_653 = DENDRO_206*DENDRO_272 + DENDRO_652;
double DENDRO_654 = DENDRO_130*DENDRO_391 + DENDRO_624;
double DENDRO_655 = DENDRO_449 + DENDRO_450;
double DENDRO_656 = DENDRO_200*DENDRO_417 + DENDRO_626;
double DENDRO_657 = 1.0*DENDRO_153;
double DENDRO_658 = -DENDRO_154*(DENDRO_577 + DENDRO_606) - DENDRO_270*(DENDRO_189*DENDRO_196 + DENDRO_387 + DENDRO_593) - DENDRO_630*(DENDRO_112 + DENDRO_475) - DENDRO_657*(DENDRO_389 + DENDRO_586);
double DENDRO_659 = DENDRO_510 + DENDRO_609;
double DENDRO_660 = 0.5*DENDRO_433;
double DENDRO_661 = DENDRO_200*DENDRO_97;
double DENDRO_662 = -DENDRO_661;
double DENDRO_663 = DENDRO_239*DENDRO_299;
double DENDRO_664 = DENDRO_125*DENDRO_391 + DENDRO_663;
double DENDRO_665 = DENDRO_213*DENDRO_272;
double DENDRO_666 = DENDRO_221*DENDRO_43;
double DENDRO_667 = DENDRO_44*DENDRO_89;
double DENDRO_668 = DENDRO_125*DENDRO_144 + DENDRO_417*DENDRO_86;
double DENDRO_669 = DENDRO_188*DENDRO_86;
double DENDRO_670 = 0.5*DENDRO_126;
double DENDRO_671 = -DENDRO_669 - DENDRO_670;
double DENDRO_672 = DENDRO_417*DENDRO_89 + DENDRO_460;
double DENDRO_673 = DENDRO_522 + DENDRO_601;
double DENDRO_674 = DENDRO_121*DENDRO_239;
double DENDRO_675 = 0.5*DENDRO_413;
double DENDRO_676 = -DENDRO_674 - DENDRO_675;
double DENDRO_677 = DENDRO_213*DENDRO_623 + DENDRO_652;
double DENDRO_678 = DENDRO_200*DENDRO_299;
double DENDRO_679 = DENDRO_449 + DENDRO_626;
double DENDRO_680 = DENDRO_239*DENDRO_48;
double DENDRO_681 = DENDRO_101*DENDRO_511 + DENDRO_680;
double DENDRO_682 = DENDRO_221*DENDRO_51;
double DENDRO_683 = DENDRO_213*DENDRO_295 + DENDRO_682;
double DENDRO_684 = DENDRO_125*DENDRO_504;
double DENDRO_685 = DENDRO_137*DENDRO_511 + DENDRO_196*DENDRO_86;
double DENDRO_686 = DENDRO_192*DENDRO_53;
double DENDRO_687 = DENDRO_111*DENDRO_272 + DENDRO_686;
double DENDRO_688 = 1.0*DENDRO_269;
double DENDRO_689 = -DENDRO_258;
double DENDRO_690 = 2.0*DENDRO_446*gt0[pp];
double DENDRO_691 = 2.0*gt1[pp];
double DENDRO_692 = DENDRO_290*DENDRO_691;
double DENDRO_693 = DENDRO_447*DENDRO_691;
double DENDRO_694 = DENDRO_448*DENDRO_538;
double DENDRO_695 = 2.0*DENDRO_291*gt3[pp];
double DENDRO_696 = DENDRO_293*DENDRO_541;
double DENDRO_697 = -DENDRO_259*DENDRO_289;
double DENDRO_698 = DENDRO_297*grad2_0_2_gt1[pp];
double DENDRO_699 = DENDRO_302*grad2_2_2_gt1[pp];
double DENDRO_700 = DENDRO_303*grad2_1_1_gt1[pp];
double DENDRO_701 = DENDRO_304*grad2_0_0_gt1[pp];
double DENDRO_702 = -DENDRO_310*grad2_1_2_gt1[pp];
double DENDRO_703 = -DENDRO_311*grad2_0_1_gt1[pp];
double DENDRO_704 = DENDRO_18*gt1[pp];
double DENDRO_705 = DENDRO_368*DENDRO_94 + DENDRO_369*DENDRO_99 + DENDRO_370*DENDRO_40 + DENDRO_373*DENDRO_704 + DENDRO_690 + DENDRO_692 + DENDRO_693 + DENDRO_694 + DENDRO_695 + DENDRO_696 + DENDRO_697 + DENDRO_698 + DENDRO_699 + DENDRO_700 + DENDRO_701 + DENDRO_702 + DENDRO_703 - DENDRO_91*(DENDRO_531*(DENDRO_107*DENDRO_331 + DENDRO_122 + DENDRO_42*DENDRO_76) + DENDRO_532*(DENDRO_110 + DENDRO_341*DENDRO_92) + DENDRO_533*(DENDRO_107*DENDRO_76 + DENDRO_140 + DENDRO_42*DENDRO_78) + DENDRO_689);
double DENDRO_706 = 0.25*DENDRO_305;
double DENDRO_707 = DENDRO_111*DENDRO_519 + DENDRO_588;
double DENDRO_708 = -DENDRO_111*DENDRO_419;
double DENDRO_709 = DENDRO_116*(DENDRO_483 + DENDRO_679) - DENDRO_154*(DENDRO_298 + DENDRO_589) - DENDRO_154*(DENDRO_298 + DENDRO_707) + DENDRO_160*(DENDRO_687 + DENDRO_708) - DENDRO_270*(DENDRO_104*DENDRO_196 + DENDRO_501) - DENDRO_283*(1.0*DENDRO_306 + DENDRO_706);
double DENDRO_710 = grad2_0_1_alpha[pp];
double DENDRO_711 = -DENDRO_106*DENDRO_18 + DENDRO_108*DENDRO_18 + DENDRO_109*DENDRO_18;
double DENDRO_712 = DENDRO_140*DENDRO_18;
double DENDRO_713 = DENDRO_141*DENDRO_18;
double DENDRO_714 = DENDRO_142*DENDRO_18;
double DENDRO_715 = DENDRO_32*(DENDRO_639 + DENDRO_704*DENDRO_79);
double DENDRO_716 = DENDRO_61*gt1[pp];
double DENDRO_717 = DENDRO_18*DENDRO_561*(DENDRO_125 + DENDRO_716) + DENDRO_567*(DENDRO_712 - DENDRO_713 - DENDRO_714 + DENDRO_715) + DENDRO_568*(DENDRO_32*(DENDRO_403*DENDRO_704 + DENDRO_559) + DENDRO_711) - 4*DENDRO_710;
double DENDRO_718 = DENDRO_150*DENDRO_213;
double DENDRO_719 = -DENDRO_10*(DENDRO_643 + alpha[pp]*(DENDRO_116*(DENDRO_611 + DENDRO_651) + DENDRO_116*(DENDRO_613 + DENDRO_614) + DENDRO_116*(DENDRO_615 - DENDRO_649) + DENDRO_116*(-DENDRO_188*DENDRO_200 + DENDRO_650) + DENDRO_116*(DENDRO_196*DENDRO_236 + DENDRO_610 + DENDRO_651) - DENDRO_154*(DENDRO_578 + DENDRO_603) - DENDRO_154*(DENDRO_607 + DENDRO_648) - DENDRO_154*(DENDRO_609 + DENDRO_648) + DENDRO_160*(DENDRO_412 + DENDRO_621) + DENDRO_160*(DENDRO_416 + DENDRO_625) + DENDRO_160*(DENDRO_416 + DENDRO_654) + DENDRO_160*(DENDRO_481 + DENDRO_653) + DENDRO_160*(DENDRO_482 + DENDRO_656) + DENDRO_160*(DENDRO_628 + DENDRO_655) - DENDRO_270*(DENDRO_353 + DENDRO_596) - DENDRO_270*(DENDRO_595 + DENDRO_644) - DENDRO_274*(DENDRO_598 + DENDRO_645) - DENDRO_274*(DENDRO_600 + DENDRO_646) - DENDRO_274*(DENDRO_196*DENDRO_200 + DENDRO_427 + DENDRO_597) - DENDRO_283*(DENDRO_460 + DENDRO_602) - DENDRO_283*(DENDRO_117*DENDRO_125 + DENDRO_647) + DENDRO_367*(DENDRO_174*DENDRO_200 + DENDRO_451) + DENDRO_634 - DENDRO_657*(DENDRO_572 + DENDRO_604) + DENDRO_658)) - DENDRO_10*(DENDRO_643 + alpha[pp]*(DENDRO_116*(DENDRO_614 + DENDRO_615) + DENDRO_116*(DENDRO_616 + DENDRO_618) + DENDRO_116*(DENDRO_184*DENDRO_200 + DENDRO_618) + DENDRO_116*(DENDRO_197*DENDRO_236 + DENDRO_611) + DENDRO_116*(-DENDRO_218*DENDRO_97 + DENDRO_613) - DENDRO_154*(DENDRO_360 + DENDRO_603) - DENDRO_154*(DENDRO_392 + DENDRO_503) - DENDRO_154*(DENDRO_392 + DENDRO_606) - DENDRO_154*(DENDRO_508 + DENDRO_607) - DENDRO_154*(DENDRO_513 + DENDRO_609) - DENDRO_154*(DENDRO_514 + DENDRO_605) + DENDRO_159*(DENDRO_528 + DENDRO_591) + DENDRO_159*(DENDRO_584 + DENDRO_592) + DENDRO_160*(DENDRO_482 + DENDRO_627) + DENDRO_160*(DENDRO_619 + DENDRO_621) + DENDRO_160*(DENDRO_622 + DENDRO_625) + DENDRO_160*(DENDRO_627 + DENDRO_628) - DENDRO_270*(0.5*DENDRO_366 + DENDRO_596) - DENDRO_270*(1.0*DENDRO_394 + DENDRO_593) - DENDRO_270*(DENDRO_184*DENDRO_236 + DENDRO_356 + DENDRO_595) - DENDRO_274*(DENDRO_421 + DENDRO_600) - DENDRO_274*(0.5*DENDRO_426 + DENDRO_597) - DENDRO_274*(DENDRO_184*DENDRO_239 + DENDRO_424 + DENDRO_598) - DENDRO_283*(DENDRO_459 + DENDRO_602) - DENDRO_283*(DENDRO_104*DENDRO_93 + DENDRO_587) + DENDRO_367*(DENDRO_172*DENDRO_236 + DENDRO_379) - DENDRO_630*(DENDRO_471 + DENDRO_629) + DENDRO_634)) + DENDRO_13*(DENDRO_320*(DENDRO_32*(DENDRO_402 + DENDRO_403*DENDRO_405) + DENDRO_406) + DENDRO_327*(DENDRO_221 + DENDRO_324*DENDRO_400) - 4*DENDRO_399 + DENDRO_401*(DENDRO_239 + DENDRO_400*DENDRO_61) + alpha[pp]*(DENDRO_116*(DENDRO_421 + DENDRO_422) + DENDRO_116*(DENDRO_424 - DENDRO_425) + DENDRO_116*(DENDRO_427 + 1.0*DENDRO_428) - DENDRO_125*DENDRO_445 - DENDRO_154*(DENDRO_410 + DENDRO_412) - DENDRO_154*(-1.0*DENDRO_413 + DENDRO_415) - DENDRO_154*(DENDRO_213*DENDRO_417 + DENDRO_418) - DENDRO_154*(DENDRO_230*DENDRO_417 + DENDRO_416) + DENDRO_160*(DENDRO_429 + DENDRO_430) + DENDRO_160*(DENDRO_431 + DENDRO_432) + DENDRO_160*(DENDRO_434 + 1.0*DENDRO_435) + DENDRO_170*(DENDRO_433 + DENDRO_435) + DENDRO_170*(DENDRO_125*DENDRO_172 + DENDRO_438) + DENDRO_170*(DENDRO_213*DENDRO_92 + DENDRO_439) + DENDRO_172*DENDRO_369 - DENDRO_191*DENDRO_201*DENDRO_382 + DENDRO_191*DENDRO_368 - DENDRO_211*DENDRO_442 - DENDRO_213*DENDRO_444 + DENDRO_367*(DENDRO_426 + DENDRO_428) + DENDRO_367*(DENDRO_172*DENDRO_230 + DENDRO_437) + DENDRO_367*(DENDRO_211*DENDRO_92 + DENDRO_436) + DENDRO_370*DENDRO_92 + DENDRO_373*DENDRO_404 - DENDRO_440*(DENDRO_186 - DENDRO_361) - DENDRO_443*(DENDRO_93 + DENDRO_95) + DENDRO_453 - DENDRO_91*(DENDRO_407 + DENDRO_87*(DENDRO_197*DENDRO_331 + DENDRO_238 + DENDRO_409*DENDRO_76) + DENDRO_88*(DENDRO_192*DENDRO_341 + DENDRO_198 + DENDRO_409*DENDRO_9) + DENDRO_90*(DENDRO_197*DENDRO_76 + DENDRO_220 + DENDRO_409*DENDRO_78)))) + DENDRO_23*(DENDRO_570 + alpha[pp]*(DENDRO_115*(DENDRO_359 + DENDRO_572) + DENDRO_116*(DENDRO_499 + DENDRO_577) + DENDRO_116*(DENDRO_506 + DENDRO_578) + DENDRO_116*(DENDRO_507 + DENDRO_579) + DENDRO_116*(DENDRO_510 + DENDRO_579) - DENDRO_154*(DENDRO_491 + DENDRO_575) - DENDRO_154*(-DENDRO_121*DENDRO_89 + DENDRO_496) - DENDRO_154*(DENDRO_236*DENDRO_48 + DENDRO_575) + DENDRO_160*(DENDRO_460 + DENDRO_525) + DENDRO_160*(DENDRO_520 + DENDRO_580) + DENDRO_160*(DENDRO_524 + DENDRO_582) + DENDRO_160*(DENDRO_580 + DENDRO_581) - DENDRO_167*(DENDRO_102*DENDRO_89 + DENDRO_467) - DENDRO_270*(DENDRO_478 + DENDRO_573) - DENDRO_270*(DENDRO_218*DENDRO_48 + DENDRO_347 + DENDRO_476) - DENDRO_274*(DENDRO_125*DENDRO_184 + DENDRO_574) - DENDRO_283*(1.0*DENDRO_466 + DENDRO_486) - DENDRO_283*(DENDRO_456 + DENDRO_48*DENDRO_89 + DENDRO_485) - DENDRO_529*(DENDRO_411 + DENDRO_584) + DENDRO_554 + DENDRO_590)) + DENDRO_23*(DENDRO_570 + alpha[pp]*(DENDRO_116*(DENDRO_364 + DENDRO_506) + DENDRO_116*(DENDRO_364 + DENDRO_514) + DENDRO_116*(DENDRO_390 + DENDRO_499) + DENDRO_116*(DENDRO_501 + DENDRO_503) + DENDRO_116*(DENDRO_507 + DENDRO_508) + DENDRO_116*(DENDRO_510 + DENDRO_513) - DENDRO_154*(DENDRO_492 + DENDRO_493) - DENDRO_154*(DENDRO_494 + DENDRO_496) - DENDRO_154*(DENDRO_189*DENDRO_44 + DENDRO_489) - DENDRO_154*(DENDRO_117*DENDRO_89 + DENDRO_494 + DENDRO_495) - DENDRO_154*(DENDRO_236*DENDRO_49 + DENDRO_490 + DENDRO_491) + DENDRO_159*(DENDRO_469 + DENDRO_471) + DENDRO_159*(DENDRO_472 + DENDRO_475) + DENDRO_160*(DENDRO_459 + DENDRO_524) + DENDRO_160*(DENDRO_515 + DENDRO_516) + DENDRO_160*(DENDRO_522 + DENDRO_525) + DENDRO_160*(DENDRO_517 + DENDRO_518 + DENDRO_520) - DENDRO_167*(DENDRO_236*DENDRO_47 + DENDRO_380) - DENDRO_270*(1.0*DENDRO_371 + DENDRO_476) - DENDRO_270*(0.5*DENDRO_396 + DENDRO_479) - DENDRO_270*(DENDRO_117*DENDRO_236 + DENDRO_350 + DENDRO_478) - DENDRO_274*(DENDRO_450 + DENDRO_484) - DENDRO_274*(DENDRO_211*DENDRO_43 + DENDRO_481) - DENDRO_283*(DENDRO_296 + DENDRO_487) - DENDRO_283*(0.5*DENDRO_455 + DENDRO_485) - DENDRO_283*(DENDRO_117*DENDRO_86 + DENDRO_454 + DENDRO_486) - DENDRO_529*(DENDRO_526 + DENDRO_528) + DENDRO_554)) + DENDRO_25*(-4*DENDRO_31 + DENDRO_321*DENDRO_54 + DENDRO_326*(-DENDRO_71 - DENDRO_73 + DENDRO_75 + DENDRO_81) + DENDRO_401*(DENDRO_62 + DENDRO_86) + alpha[pp]*(-DENDRO_114*DENDRO_223*DENDRO_288 + DENDRO_116*(-1.0*DENDRO_105 + DENDRO_113) + DENDRO_116*(-1.0*DENDRO_126 + DENDRO_132) + DENDRO_116*(DENDRO_130*DENDRO_299 + DENDRO_458) + DENDRO_116*(DENDRO_206*DENDRO_43 + DENDRO_460) + DENDRO_116*(DENDRO_213*DENDRO_48 + DENDRO_459) - DENDRO_125*DENDRO_280 - DENDRO_154*(DENDRO_456 + 1.0*DENDRO_457) - DENDRO_154*(-DENDRO_152*DENDRO_86 + DENDRO_454) + DENDRO_160*(DENDRO_463 + 1.0*DENDRO_464) + DENDRO_160*(DENDRO_137*DENDRO_86 + DENDRO_461) - DENDRO_167*(DENDRO_455 + DENDRO_457) - DENDRO_167*(DENDRO_130*DENDRO_47 + DENDRO_466) + DENDRO_170*(DENDRO_462 + DENDRO_464) + DENDRO_170*(DENDRO_125*DENDRO_47 + DENDRO_465) + DENDRO_205*DENDRO_261 + DENDRO_225*DENDRO_265 - DENDRO_241*DENDRO_281*DENDRO_314 + DENDRO_243*DENDRO_263 + DENDRO_268*DENDRO_372 - DENDRO_271*(DENDRO_119 - DENDRO_135) - DENDRO_275*(-DENDRO_272 + DENDRO_95) - DENDRO_285*DENDRO_467 - DENDRO_287*DENDRO_468 + DENDRO_315 - DENDRO_91*(DENDRO_85 + DENDRO_87*(DENDRO_331*DENDRO_49 + DENDRO_51*DENDRO_76 + DENDRO_65) + DENDRO_88*(DENDRO_341*DENDRO_44 + DENDRO_52) + DENDRO_90*(DENDRO_49*DENDRO_76 + DENDRO_51*DENDRO_78 + DENDRO_74)))) + DENDRO_29*(-4*DENDRO_317 + DENDRO_319*DENDRO_321 + DENDRO_327*(DENDRO_218 + DENDRO_325) + DENDRO_335*(-DENDRO_328 + DENDRO_329 - DENDRO_330 + DENDRO_334) + alpha[pp]*(DENDRO_102*DENDRO_370 + DENDRO_116*(DENDRO_353 + DENDRO_354) + DENDRO_116*(DENDRO_356 + 1.0*DENDRO_357) - DENDRO_154*(DENDRO_344 + DENDRO_345) - DENDRO_154*(DENDRO_350 + 1.0*DENDRO_351) - DENDRO_154*(DENDRO_218*DENDRO_348 + DENDRO_347) + DENDRO_160*(DENDRO_358 + DENDRO_360) + DENDRO_160*(DENDRO_117*DENDRO_230 + DENDRO_362) + DENDRO_160*(DENDRO_130*DENDRO_184 + DENDRO_363) + DENDRO_160*(DENDRO_206*DENDRO_365 + DENDRO_364) - DENDRO_167*(DENDRO_349 + DENDRO_351) - DENDRO_167*(DENDRO_102*DENDRO_206 + DENDRO_371) + DENDRO_174*DENDRO_368 - DENDRO_181*DENDRO_237*DENDRO_382 + DENDRO_181*DENDRO_369 - DENDRO_206*DENDRO_378 - DENDRO_211*DENDRO_377 - DENDRO_287*DENDRO_379 + DENDRO_333*DENDRO_373 + DENDRO_367*(DENDRO_355 + DENDRO_357) + DENDRO_367*(DENDRO_102*DENDRO_211 + DENDRO_366) - DENDRO_374*(DENDRO_117 + DENDRO_119) - DENDRO_375*(DENDRO_184 + DENDRO_186) - DENDRO_380*DENDRO_381 + DENDRO_398 - DENDRO_91*(DENDRO_336 + DENDRO_87*(DENDRO_10*DENDRO_338 + DENDRO_182*DENDRO_331 + DENDRO_340*DENDRO_76) + DENDRO_88*(DENDRO_183 + DENDRO_338*DENDRO_341 + DENDRO_340*DENDRO_9) + DENDRO_90*(DENDRO_182*DENDRO_76 + DENDRO_338*DENDRO_9 + DENDRO_340*DENDRO_78)))) - DENDRO_9*(DENDRO_717 + alpha[pp]*(DENDRO_115*(DENDRO_411 + DENDRO_583 + DENDRO_592) + DENDRO_115*(DENDRO_526 + DENDRO_527 + DENDRO_591) + DENDRO_116*(DENDRO_619 + DENDRO_677) + DENDRO_116*(DENDRO_622 + DENDRO_676) + DENDRO_116*(DENDRO_678 + DENDRO_679) - DENDRO_154*(DENDRO_132 + DENDRO_671) - DENDRO_154*(DENDRO_459 + DENDRO_672) - DENDRO_154*(DENDRO_582 + DENDRO_673) - DENDRO_154*(DENDRO_517 + DENDRO_581 + DENDRO_647) + DENDRO_160*(DENDRO_683 - DENDRO_89*DENDRO_97) + DENDRO_160*(DENDRO_684 + DENDRO_685) + DENDRO_160*(DENDRO_197*DENDRO_86 + DENDRO_681) + DENDRO_160*(DENDRO_200*DENDRO_43 + DENDRO_687) + DENDRO_170*(DENDRO_468 + DENDRO_89*DENDRO_92) - DENDRO_270*(DENDRO_362 + DENDRO_659) - DENDRO_274*(DENDRO_429 + DENDRO_664) - DENDRO_274*(DENDRO_660 + DENDRO_662) - DENDRO_274*(DENDRO_431 + DENDRO_665 + DENDRO_666) - DENDRO_283*(0.5*DENDRO_465 + DENDRO_668) - DENDRO_283*(DENDRO_43*DENDRO_89 + DENDRO_463 + DENDRO_667) - DENDRO_688*(DENDRO_359 + DENDRO_571 + DENDRO_604) + DENDRO_705 + DENDRO_709)) - DENDRO_9*(DENDRO_717 + alpha[pp]*(DENDRO_116*(DENDRO_415 + DENDRO_676) + DENDRO_116*(DENDRO_418 + DENDRO_653) + DENDRO_116*(DENDRO_418 + DENDRO_677) + DENDRO_116*(DENDRO_483 + DENDRO_656) + DENDRO_116*(DENDRO_574 + DENDRO_654) + DENDRO_116*(DENDRO_655 + DENDRO_678) - DENDRO_154*(DENDRO_460 + DENDRO_673) - DENDRO_154*(DENDRO_515 + DENDRO_707) - DENDRO_154*(DENDRO_518 + DENDRO_671) - DENDRO_154*(DENDRO_522 + DENDRO_672) + DENDRO_160*(DENDRO_681 + DENDRO_684) + DENDRO_160*(DENDRO_683 + DENDRO_718) + DENDRO_160*(DENDRO_239*DENDRO_49 + DENDRO_685) + DENDRO_160*(DENDRO_682 + DENDRO_718 + DENDRO_89*DENDRO_93) + DENDRO_160*(DENDRO_200*DENDRO_44 + DENDRO_686 + DENDRO_708) + DENDRO_170*(DENDRO_200*DENDRO_42 + DENDRO_452) - DENDRO_270*(DENDRO_363 + DENDRO_659) - DENDRO_270*(DENDRO_211*DENDRO_48 + DENDRO_605) - DENDRO_274*(0.5*DENDRO_438 + DENDRO_664) - DENDRO_274*(1.0*DENDRO_439 + DENDRO_665) - DENDRO_274*(DENDRO_200*DENDRO_93 + DENDRO_434 + DENDRO_662) - DENDRO_283*(DENDRO_461 + DENDRO_668) - DENDRO_283*(0.5*DENDRO_462 + DENDRO_667) - DENDRO_283*(DENDRO_300 + DENDRO_53*DENDRO_93 + DENDRO_706) - DENDRO_657*(DENDRO_112 + DENDRO_472 + DENDRO_474) - DENDRO_657*(DENDRO_469 + DENDRO_470 + DENDRO_629) - DENDRO_688*(DENDRO_389 + DENDRO_500 + DENDRO_585) + DENDRO_705));
double DENDRO_720 = DENDRO_18*DENDRO_719;
double DENDRO_721 = (1.0/12.0)*chi[pp];
double DENDRO_722 = grad_1_beta0[pp];
double DENDRO_723 = grad_1_beta2[pp];
double DENDRO_724 = (1.0/3.0)*At1[pp];
double DENDRO_725 = (2.0/3.0)*DENDRO_3;
double DENDRO_726 = At4[pp]*DENDRO_10;
double DENDRO_727 = -At3[pp]*DENDRO_13 + DENDRO_20 + DENDRO_726;
double DENDRO_728 = -At1[pp]*DENDRO_25 + At3[pp]*DENDRO_9 - At4[pp]*DENDRO_23;
double DENDRO_729 = -At1[pp]*DENDRO_23 + At3[pp]*DENDRO_10 - At4[pp]*DENDRO_29;
double DENDRO_730 = 6.0*DENDRO_67;
double DENDRO_731 = 6.0*DENDRO_82;
double DENDRO_732 = 6.0*DENDRO_55;
double DENDRO_733 = DENDRO_101*DENDRO_149;
double DENDRO_734 = DENDRO_137*DENDRO_149;
double DENDRO_735 = -DENDRO_208 + DENDRO_209 + DENDRO_210;
double DENDRO_736 = DENDRO_47*DENDRO_735;
double DENDRO_737 = DENDRO_734 + DENDRO_736;
double DENDRO_738 = DENDRO_143*DENDRO_419;
double DENDRO_739 = -DENDRO_134*DENDRO_423;
double DENDRO_740 = DENDRO_193 - DENDRO_198 + DENDRO_199;
double DENDRO_741 = DENDRO_138*DENDRO_623;
double DENDRO_742 = -DENDRO_227 + DENDRO_228 + DENDRO_229;
double DENDRO_743 = DENDRO_184*DENDRO_66;
double DENDRO_744 = DENDRO_144*DENDRO_742 + DENDRO_743;
double DENDRO_745 = DENDRO_149*DENDRO_295;
double DENDRO_746 = DENDRO_157*DENDRO_299 + 0.25*DENDRO_50*DENDRO_735;
double DENDRO_747 = DENDRO_145 + DENDRO_151;
double DENDRO_748 = DENDRO_137*DENDRO_742;
double DENDRO_749 = 1.0*DENDRO_115;
double DENDRO_750 = DENDRO_137*DENDRO_143;
double DENDRO_751 = DENDRO_42*DENDRO_735 + DENDRO_750;
double DENDRO_752 = DENDRO_203*DENDRO_204;
double DENDRO_753 = DENDRO_204*DENDRO_224;
double DENDRO_754 = DENDRO_204*DENDRO_242;
double DENDRO_755 = DENDRO_267*DENDRO_32;
double DENDRO_756 = grad_2_beta0[pp];
double DENDRO_757 = grad_2_beta1[pp];
double DENDRO_758 = (1.0/3.0)*At2[pp];
double DENDRO_759 = (2.0/3.0)*DENDRO_1;
double DENDRO_760 = At2[pp]*DENDRO_9 - At4[pp]*DENDRO_13 + At5[pp]*DENDRO_10;
double DENDRO_761 = -At2[pp]*DENDRO_25 + At4[pp]*DENDRO_9 - At5[pp]*DENDRO_23;
double DENDRO_762 = -At5[pp]*DENDRO_29 + DENDRO_24 + DENDRO_726;
double DENDRO_763 = DENDRO_107*DENDRO_143;
double DENDRO_764 = DENDRO_149*DENDRO_346;
double DENDRO_765 = DENDRO_215 + DENDRO_216 - DENDRO_217;
double DENDRO_766 = 0.25*DENDRO_134*DENDRO_181;
double DENDRO_767 = DENDRO_138*DENDRO_361;
double DENDRO_768 = DENDRO_233 - DENDRO_234 + DENDRO_235;
double DENDRO_769 = DENDRO_135*DENDRO_138;
double DENDRO_770 = DENDRO_182*DENDRO_66;
double DENDRO_771 = DENDRO_143*DENDRO_161;
double DENDRO_772 = -DENDRO_766;
double DENDRO_773 = DENDRO_143*DENDRO_346;
double DENDRO_774 = DENDRO_102*DENDRO_143;
double DENDRO_775 = (2.0/3.0)*DENDRO_0;
double DENDRO_776 = 2*At4[pp];
double DENDRO_777 = At3[pp]*DENDRO_26;
double DENDRO_778 = DENDRO_18*DENDRO_776;
double DENDRO_779 = DENDRO_32*DENDRO_400;
double DENDRO_780 = DENDRO_18*DENDRO_83;
double DENDRO_781 = DENDRO_172*DENDRO_740;
double DENDRO_782 = 1.0*DENDRO_735;
double DENDRO_783 = 0.25*DENDRO_750;
double DENDRO_784 = DENDRO_134*DENDRO_361;
double DENDRO_785 = DENDRO_740*DENDRO_92;
double DENDRO_786 = (1.0/3.0)*At4[pp];
double DENDRO_787 = DENDRO_135*DENDRO_742;
double DENDRO_788 = -DENDRO_361*DENDRO_742 + DENDRO_610;
double DENDRO_789 = DENDRO_624 - DENDRO_784;
double DENDRO_790 = DENDRO_181*DENDRO_742;
double DENDRO_791 = DENDRO_174*DENDRO_768;
double DENDRO_792 = DENDRO_138*DENDRO_181;
double DENDRO_793 = DENDRO_102*DENDRO_768;

// Dendro: printing variables
//--
At_rhs00[pp] = (4.0/3.0)*At0[pp]*DENDRO_0 - DENDRO_1*DENDRO_2 - DENDRO_2*DENDRO_3 + DENDRO_4*DENDRO_5 + DENDRO_6*DENDRO_7 + DENDRO_721*(-12*DENDRO_31 + DENDRO_316*(-DENDRO_116*(DENDRO_105 - DENDRO_113) - DENDRO_116*(DENDRO_126 + DENDRO_133) - DENDRO_116*(DENDRO_136 + DENDRO_137*DENDRO_139) - DENDRO_116*(DENDRO_145 + DENDRO_149*DENDRO_43) - DENDRO_116*(DENDRO_143*DENDRO_48 + DENDRO_151) + DENDRO_134*DENDRO_280 + DENDRO_154*(DENDRO_163 + 1.0*DENDRO_164) - DENDRO_154*(-DENDRO_139*DENDRO_49 + DENDRO_152*DENDRO_66) + DENDRO_157*DENDRO_282*DENDRO_288 - DENDRO_160*(DENDRO_156 + 1.0*DENDRO_158) - DENDRO_160*(DENDRO_134*DENDRO_161 + 1.0*DENDRO_137*DENDRO_66) + DENDRO_167*(DENDRO_162 + DENDRO_164) + DENDRO_167*(DENDRO_165 + DENDRO_166) - DENDRO_170*(DENDRO_155 + DENDRO_158) - DENDRO_170*(DENDRO_168 + DENDRO_169) - DENDRO_203*DENDRO_205 - DENDRO_224*DENDRO_225 - DENDRO_242*DENDRO_243 - DENDRO_267*DENDRO_268 + DENDRO_271*(DENDRO_120 + DENDRO_135) + DENDRO_275*(DENDRO_272 + DENDRO_96) + DENDRO_281*DENDRO_283*DENDRO_66 + DENDRO_284*DENDRO_285 + DENDRO_286*DENDRO_287 + DENDRO_315 - DENDRO_91*(DENDRO_53*DENDRO_88 + DENDRO_85 + DENDRO_86*DENDRO_87 + DENDRO_89*DENDRO_90)) + DENDRO_54*DENDRO_57 - DENDRO_69*(-DENDRO_62 + DENDRO_66) + DENDRO_720*gt0[pp] - DENDRO_83*(DENDRO_71 + DENDRO_73 - DENDRO_75 - DENDRO_81)) - alpha[pp]*(-At0[pp]*K[pp] + DENDRO_19*(At0[pp]*DENDRO_9 - At1[pp]*DENDRO_13 + At2[pp]*DENDRO_10) + DENDRO_27*(-At0[pp]*DENDRO_25 + DENDRO_20 + DENDRO_24) + DENDRO_30*(-At0[pp]*DENDRO_23 + At1[pp]*DENDRO_10 - At2[pp]*DENDRO_29)) + beta0[pp]*agrad_0_At0[pp] + beta1[pp]*agrad_1_At0[pp] + beta2[pp]*agrad_2_At0[pp];
//--
At_rhs01[pp] = At0[pp]*DENDRO_722 - At1[pp]*DENDRO_725 + At2[pp]*DENDRO_723 + At3[pp]*DENDRO_4 + At4[pp]*DENDRO_6 + DENDRO_0*DENDRO_724 + DENDRO_1*DENDRO_724 + DENDRO_721*(-DENDRO_18*DENDRO_730*(DENDRO_134 - DENDRO_716) + DENDRO_316*(DENDRO_116*(-DENDRO_299*DENDRO_740 + DENDRO_679) - DENDRO_116*(-DENDRO_622 + DENDRO_674 + DENDRO_675) + DENDRO_116*(-DENDRO_143*DENDRO_623 + DENDRO_149*DENDRO_419 + DENDRO_652) + DENDRO_154*(DENDRO_741 + DENDRO_744) + DENDRO_154*(DENDRO_745 + DENDRO_746) + DENDRO_154*(DENDRO_157*DENDRO_417 + DENDRO_747) + DENDRO_154*(DENDRO_133 + DENDRO_669 + DENDRO_670) + DENDRO_160*(-DENDRO_43*DENDRO_740 + DENDRO_687) - DENDRO_160*(DENDRO_134*DENDRO_504 + DENDRO_134*DENDRO_519 + DENDRO_196*DENDRO_66) + DENDRO_160*(-DENDRO_134*DENDRO_623 - DENDRO_197*DENDRO_66 + DENDRO_680) + DENDRO_160*(-DENDRO_143*DENDRO_295 + DENDRO_157*DENDRO_97 + DENDRO_682) - DENDRO_170*(DENDRO_157*DENDRO_92 + DENDRO_286) + DENDRO_269*(DENDRO_733 + DENDRO_737) + DENDRO_270*(-DENDRO_362 + DENDRO_509 + DENDRO_608) + DENDRO_274*(-DENDRO_660 + DENDRO_661) - DENDRO_274*(-DENDRO_134*DENDRO_391 + DENDRO_663 + DENDRO_739) - DENDRO_274*(-DENDRO_143*DENDRO_272 + DENDRO_666 + DENDRO_738) + DENDRO_283*(DENDRO_156 + DENDRO_157*DENDRO_43 + DENDRO_157*DENDRO_44) + DENDRO_283*(0.25*DENDRO_168 + 0.5*DENDRO_169 + DENDRO_417*DENDRO_66) - DENDRO_40*DENDRO_753 + DENDRO_690 + DENDRO_692 + DENDRO_693 + DENDRO_694 + DENDRO_695 + DENDRO_696 + DENDRO_697 + DENDRO_698 + DENDRO_699 + DENDRO_700 + DENDRO_701 + DENDRO_702 + DENDRO_703 - DENDRO_704*DENDRO_755 + DENDRO_709 - DENDRO_749*(DENDRO_149*DENDRO_92 + DENDRO_751) - DENDRO_749*(DENDRO_134*DENDRO_174 + DENDRO_138*DENDRO_172 + DENDRO_748) - DENDRO_752*DENDRO_94 - DENDRO_754*DENDRO_99 - DENDRO_91*(DENDRO_111*DENDRO_532 + DENDRO_125*DENDRO_531 + DENDRO_213*DENDRO_533 + DENDRO_689)) + DENDRO_704*DENDRO_719 - 12*DENDRO_710 - DENDRO_731*(-DENDRO_712 + DENDRO_713 + DENDRO_714 - DENDRO_715) + DENDRO_732*(DENDRO_32*(DENDRO_37*DENDRO_704 + DENDRO_559) + DENDRO_711)) - alpha[pp]*(-At1[pp]*K[pp] + DENDRO_19*DENDRO_727 + DENDRO_27*DENDRO_728 + DENDRO_30*DENDRO_729) + beta0[pp]*agrad_0_At1[pp] + beta1[pp]*agrad_1_At1[pp] + beta2[pp]*agrad_2_At1[pp];
//--
At_rhs02[pp] = At0[pp]*DENDRO_756 + At1[pp]*DENDRO_757 - At2[pp]*DENDRO_759 + At4[pp]*DENDRO_4 + At5[pp]*DENDRO_6 + DENDRO_0*DENDRO_758 + DENDRO_3*DENDRO_758 + DENDRO_721*(DENDRO_316*(-DENDRO_100*DENDRO_752 - DENDRO_116*(DENDRO_497 - DENDRO_498 + DENDRO_576) + DENDRO_116*(-DENDRO_149*DENDRO_504 - DENDRO_43*DENDRO_765 + DENDRO_773) - DENDRO_116*(DENDRO_299*DENDRO_768 + DENDRO_766 + DENDRO_767) + DENDRO_116*(DENDRO_346*DENDRO_742 - DENDRO_767 + DENDRO_772) - DENDRO_118*DENDRO_754 - DENDRO_154*(DENDRO_121*DENDRO_139 - DENDRO_769 - DENDRO_770) - DENDRO_154*(DENDRO_121*DENDRO_157 - DENDRO_149*DENDRO_161 - DENDRO_51*DENDRO_765) + DENDRO_154*(DENDRO_48*DENDRO_768 + DENDRO_769 + DENDRO_770) - DENDRO_160*(DENDRO_136 + DENDRO_744) - DENDRO_160*(DENDRO_746 + DENDRO_771) - DENDRO_160*(DENDRO_157*DENDRO_365 + DENDRO_747) - DENDRO_160*(DENDRO_136 + DENDRO_138*DENDRO_519 + DENDRO_743) + DENDRO_167*(DENDRO_102*DENDRO_157 + DENDRO_284) + DENDRO_270*(DENDRO_477 - DENDRO_573) - DENDRO_270*(-DENDRO_135*DENDRO_149 - DENDRO_48*DENDRO_765 + DENDRO_764) + DENDRO_273*(DENDRO_751 + DENDRO_763) + DENDRO_274*(DENDRO_134*DENDRO_184 + 0.25*DENDRO_748) + DENDRO_283*(0.25*DENDRO_165 + 1.0*DENDRO_166) + DENDRO_283*(DENDRO_157*DENDRO_48 + DENDRO_157*DENDRO_49 + DENDRO_163) - DENDRO_45*DENDRO_753 + DENDRO_535 + DENDRO_537 + DENDRO_539 + DENDRO_540 + DENDRO_542 + DENDRO_544 + DENDRO_546 + DENDRO_547 + DENDRO_548 + DENDRO_549 + DENDRO_550 + DENDRO_551 + DENDRO_552 - DENDRO_553*DENDRO_755 + DENDRO_590 - DENDRO_749*(DENDRO_737 + DENDRO_774) - DENDRO_91*(DENDRO_104*DENDRO_532 + DENDRO_130*DENDRO_531 + DENDRO_206*DENDRO_533 + DENDRO_530)) + DENDRO_553*DENDRO_719 - 12*DENDRO_555 + DENDRO_569*DENDRO_732 - DENDRO_730*(DENDRO_556 + DENDRO_557 - DENDRO_558 - DENDRO_560) - DENDRO_731*(DENDRO_562 + DENDRO_563 - DENDRO_564 - DENDRO_566)) - alpha[pp]*(-At2[pp]*K[pp] + DENDRO_19*DENDRO_760 + DENDRO_27*DENDRO_761 + DENDRO_30*DENDRO_762) + beta0[pp]*agrad_0_At2[pp] + beta1[pp]*agrad_1_At2[pp] + beta2[pp]*agrad_2_At2[pp];
//--
At_rhs11[pp] = (4.0/3.0)*At3[pp]*DENDRO_1 - At3[pp]*DENDRO_725 - At3[pp]*DENDRO_775 + DENDRO_5*DENDRO_722 + DENDRO_721*(DENDRO_316*(-DENDRO_116*(DENDRO_420 - 1.0*DENDRO_422) + DENDRO_116*(DENDRO_427 - 1.0*DENDRO_781) - DENDRO_116*(DENDRO_423*DENDRO_742 + DENDRO_425) + DENDRO_134*DENDRO_445 + DENDRO_143*DENDRO_444 + DENDRO_154*(DENDRO_413 - DENDRO_415) + DENDRO_154*(DENDRO_143*DENDRO_417 + DENDRO_150*DENDRO_735) + DENDRO_154*(DENDRO_417*DENDRO_742 + DENDRO_784) + DENDRO_154*(DENDRO_44*DENDRO_782 + DENDRO_783) + DENDRO_160*(DENDRO_430 + DENDRO_739) + DENDRO_160*(DENDRO_432 + DENDRO_738) + DENDRO_160*(DENDRO_434 - 1.0*DENDRO_785) + DENDRO_170*(DENDRO_433 - DENDRO_785) + DENDRO_170*(-DENDRO_134*DENDRO_172 + DENDRO_438) + DENDRO_170*(-DENDRO_143*DENDRO_92 + DENDRO_439) - DENDRO_172*DENDRO_754 + 6.0*DENDRO_191*DENDRO_273*DENDRO_740 - DENDRO_191*DENDRO_752 + DENDRO_367*(DENDRO_426 - DENDRO_781) + DENDRO_367*(DENDRO_436 - DENDRO_735*DENDRO_92) + DENDRO_367*(-DENDRO_172*DENDRO_742 + DENDRO_437) - DENDRO_404*DENDRO_755 + DENDRO_440*(DENDRO_187 + DENDRO_361) + DENDRO_442*DENDRO_735 + DENDRO_443*(DENDRO_408 + DENDRO_96) + DENDRO_453 - DENDRO_753*DENDRO_92 - DENDRO_91*(DENDRO_200*DENDRO_88 + DENDRO_221*DENDRO_90 + DENDRO_239*DENDRO_87 + DENDRO_407)) - 12*DENDRO_399 + DENDRO_56*(DENDRO_32*(DENDRO_37*DENDRO_405 + DENDRO_402) + DENDRO_406) + DENDRO_69*(DENDRO_239 - DENDRO_779*(-DENDRO_58 + DENDRO_59 + DENDRO_60)) + DENDRO_720*gt3[pp] + DENDRO_780*(DENDRO_221 - DENDRO_779*(DENDRO_322 + DENDRO_323 - DENDRO_77))) + DENDRO_723*DENDRO_776 - alpha[pp]*(-At3[pp]*K[pp] + DENDRO_19*DENDRO_728 + DENDRO_727*DENDRO_777 + DENDRO_729*DENDRO_778) + beta0[pp]*agrad_0_At3[pp] + beta1[pp]*agrad_1_At3[pp] + beta2[pp]*agrad_2_At3[pp];
//--
At_rhs12[pp] = At1[pp]*DENDRO_756 + At2[pp]*DENDRO_722 + At3[pp]*DENDRO_757 - At4[pp]*DENDRO_775 + At5[pp]*DENDRO_723 + DENDRO_1*DENDRO_786 + DENDRO_3*DENDRO_786 + DENDRO_721*(-DENDRO_18*DENDRO_731*(-DENDRO_642 + DENDRO_735) + DENDRO_316*(DENDRO_116*(DENDRO_188*DENDRO_740 + DENDRO_650) + DENDRO_116*(-DENDRO_196*DENDRO_768 + DENDRO_788) + DENDRO_116*(DENDRO_342*DENDRO_742 + DENDRO_788) - DENDRO_116*(DENDRO_504*DENDRO_735 + DENDRO_649 + DENDRO_765*DENDRO_93) + DENDRO_116*(-DENDRO_519*DENDRO_735 + DENDRO_612 - DENDRO_623*DENDRO_735) + DENDRO_153*(DENDRO_733 + DENDRO_736 + DENDRO_774) - DENDRO_154*(DENDRO_139*DENDRO_188 + DENDRO_772 - DENDRO_787) - DENDRO_154*(-DENDRO_161*DENDRO_735 - DENDRO_44*DENDRO_765 + DENDRO_773) + DENDRO_154*(DENDRO_417*DENDRO_768 + DENDRO_766 + DENDRO_787) + DENDRO_160*(-DENDRO_138*DENDRO_391 + DENDRO_789) + DENDRO_160*(-DENDRO_365*DENDRO_740 + DENDRO_655) + DENDRO_160*(-DENDRO_623*DENDRO_742 + DENDRO_789) + DENDRO_160*(-DENDRO_149*DENDRO_272 + DENDRO_652 - 0.25*DENDRO_763) + DENDRO_160*(-DENDRO_295*DENDRO_735 + DENDRO_620 - DENDRO_783) + DENDRO_160*(-DENDRO_417*DENDRO_740 + DENDRO_482 + DENDRO_626) - DENDRO_185*DENDRO_754 - DENDRO_194*DENDRO_752 + DENDRO_270*(DENDRO_594 - DENDRO_644) - DENDRO_270*(-DENDRO_135*DENDRO_735 + DENDRO_346*DENDRO_735 - DENDRO_417*DENDRO_765) - DENDRO_274*(-DENDRO_391*DENDRO_742 + DENDRO_645) - DENDRO_274*(-DENDRO_196*DENDRO_740 - DENDRO_197*DENDRO_740 + DENDRO_427) - DENDRO_274*(-DENDRO_272*DENDRO_735 + DENDRO_599 + DENDRO_646) + DENDRO_283*(DENDRO_117*DENDRO_134 + DENDRO_741) + DENDRO_283*(DENDRO_145 + DENDRO_745 + DENDRO_771) + DENDRO_367*(-DENDRO_174*DENDRO_740 + DENDRO_451) - DENDRO_632*DENDRO_755 + DENDRO_633 + DENDRO_658 - DENDRO_753*DENDRO_98 - DENDRO_91*(DENDRO_178*DENDRO_532 + DENDRO_211*DENDRO_533 + DENDRO_230*DENDRO_531 + DENDRO_631)) + DENDRO_632*DENDRO_719 - 12*DENDRO_635 - DENDRO_730*(-DENDRO_636 + DENDRO_637 + DENDRO_638 - DENDRO_640) + DENDRO_732*(DENDRO_32*(DENDRO_37*DENDRO_632 + DENDRO_565) + DENDRO_641)) - alpha[pp]*(-At4[pp]*K[pp] + DENDRO_19*DENDRO_761 + DENDRO_760*DENDRO_777 + DENDRO_762*DENDRO_778) + beta0[pp]*agrad_0_At4[pp] + beta1[pp]*agrad_1_At4[pp] + beta2[pp]*agrad_2_At4[pp];
//--
At_rhs22[pp] = (4.0/3.0)*At5[pp]*DENDRO_3 - At5[pp]*DENDRO_759 - At5[pp]*DENDRO_775 + DENDRO_7*DENDRO_756 + DENDRO_721*(DENDRO_316*(DENDRO_102*DENDRO_138*DENDRO_381 - DENDRO_102*DENDRO_753 - DENDRO_116*(DENDRO_352 - 1.0*DENDRO_354) - DENDRO_116*(0.25*DENDRO_790 + 1.0*DENDRO_791) + DENDRO_149*DENDRO_378 + DENDRO_154*(DENDRO_343 - 1.0*DENDRO_345) + DENDRO_154*(0.25*DENDRO_792 + 1.0*DENDRO_793) - DENDRO_154*(-DENDRO_348*DENDRO_765 + DENDRO_764) - DENDRO_160*(DENDRO_117*DENDRO_742 + DENDRO_767) - DENDRO_160*(DENDRO_138*DENDRO_184 + DENDRO_787) - DENDRO_160*(DENDRO_144*DENDRO_735 + DENDRO_149*DENDRO_365) - DENDRO_160*(DENDRO_49*DENDRO_782 + 0.25*DENDRO_734) + DENDRO_167*(DENDRO_792 + DENDRO_793) + DENDRO_167*(DENDRO_102*DENDRO_149 + DENDRO_47*DENDRO_765) + DENDRO_174*DENDRO_287*DENDRO_742 - DENDRO_174*DENDRO_752 + 6.0*DENDRO_181*DENDRO_269*DENDRO_768 - DENDRO_181*DENDRO_754 - DENDRO_333*DENDRO_755 - DENDRO_367*(DENDRO_790 + DENDRO_791) - DENDRO_367*(DENDRO_102*DENDRO_735 + DENDRO_137*DENDRO_765) + DENDRO_374*(DENDRO_120 + DENDRO_339) + DENDRO_375*(DENDRO_187 + DENDRO_337) + DENDRO_377*DENDRO_735 + DENDRO_398 - DENDRO_91*(DENDRO_189*DENDRO_88 + DENDRO_218*DENDRO_90 + DENDRO_236*DENDRO_87 + DENDRO_336)) - 12*DENDRO_317 + DENDRO_319*DENDRO_57 - DENDRO_68*(DENDRO_328 - DENDRO_329 + DENDRO_330 - DENDRO_334) + DENDRO_720*gt5[pp] - DENDRO_780*(-DENDRO_325 + DENDRO_765)) + DENDRO_757*DENDRO_776 - alpha[pp]*(At5[pp]*DENDRO_26*DENDRO_762 - At5[pp]*K[pp] + DENDRO_30*DENDRO_761 + DENDRO_760*DENDRO_778) + beta0[pp]*agrad_0_At5[pp] + beta1[pp]*agrad_1_At5[pp] + beta2[pp]*agrad_2_At5[pp];
// Dendro: reduced ops: 3279
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    bssn::timer::t_rhs_At.stop();

    bssn::timer::t_rhs_K.start();
    //K_rhs
    for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.K_rhs]
vnames = ['K_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 3960 
// Dendro: printing temp variables
double DENDRO_0 = gt0[pp]*gt3[pp];
double DENDRO_1 = pow(gt1[pp], 2);
double DENDRO_2 = DENDRO_0 - DENDRO_1;
double DENDRO_3 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
double DENDRO_4 = 0.5*grad_2_gt5[pp];
double DENDRO_5 = pow(gt2[pp], 2);
double DENDRO_6 = gt0[pp]*gt5[pp];
double DENDRO_7 = DENDRO_5 - DENDRO_6;
double DENDRO_8 = grad_1_gt5[pp];
double DENDRO_9 = -0.5*DENDRO_8 + 1.0*grad_2_gt4[pp];
double DENDRO_10 = gt2[pp]*gt4[pp];
double DENDRO_11 = -DENDRO_10 + gt1[pp]*gt5[pp];
double DENDRO_12 = grad_0_gt5[pp];
double DENDRO_13 = -0.5*DENDRO_12 + 1.0*grad_2_gt2[pp];
double DENDRO_14 = 1.0/chi[pp];
double DENDRO_15 = grad_2_chi[pp];
double DENDRO_16 = grad_1_chi[pp];
double DENDRO_17 = grad_0_chi[pp];
double DENDRO_18 = DENDRO_11*DENDRO_17 + DENDRO_15*DENDRO_3 + DENDRO_16*DENDRO_7;
double DENDRO_19 = DENDRO_14*DENDRO_18;
double DENDRO_20 = 0.5*gt5[pp];
double DENDRO_21 = grad_1_alpha[pp];
double DENDRO_22 = pow(gt4[pp], 2);
double DENDRO_23 = gt3[pp]*gt5[pp];
double DENDRO_24 = DENDRO_1*gt5[pp] - 2*DENDRO_10*gt1[pp] + DENDRO_22*gt0[pp] - DENDRO_23*gt0[pp] + DENDRO_5*gt3[pp];
double DENDRO_25 = 1.0/DENDRO_24;
double DENDRO_26 = DENDRO_21*DENDRO_25;
double DENDRO_27 = gt2[pp]*gt3[pp];
double DENDRO_28 = gt1[pp]*gt4[pp];
double DENDRO_29 = DENDRO_27 - DENDRO_28;
double DENDRO_30 = DENDRO_22 - DENDRO_23;
double DENDRO_31 = DENDRO_11*DENDRO_16 + DENDRO_15*DENDRO_29 + DENDRO_17*DENDRO_30;
double DENDRO_32 = DENDRO_14*DENDRO_31;
double DENDRO_33 = grad_0_alpha[pp];
double DENDRO_34 = DENDRO_25*DENDRO_33;
double DENDRO_35 = grad_2_alpha[pp];
double DENDRO_36 = -DENDRO_0 + DENDRO_1;
double DENDRO_37 = DENDRO_25*DENDRO_36;
double DENDRO_38 = DENDRO_25*DENDRO_3;
double DENDRO_39 = DENDRO_25*DENDRO_29;
double DENDRO_40 = DENDRO_15*DENDRO_36 + DENDRO_16*DENDRO_3 + DENDRO_17*DENDRO_29;
double DENDRO_41 = DENDRO_25*chi[pp];
double DENDRO_42 = -DENDRO_5 + DENDRO_6;
double DENDRO_43 = 0.5*grad_1_gt3[pp];
double DENDRO_44 = grad_2_gt3[pp];
double DENDRO_45 = -0.5*DENDRO_44 + 1.0*grad_1_gt4[pp];
double DENDRO_46 = grad_0_gt3[pp];
double DENDRO_47 = -0.5*DENDRO_46 + 1.0*grad_1_gt1[pp];
double DENDRO_48 = DENDRO_14*DENDRO_40;
double DENDRO_49 = 0.5*gt3[pp];
double DENDRO_50 = DENDRO_25*DENDRO_35;
double DENDRO_51 = DENDRO_25*DENDRO_7;
double DENDRO_52 = DENDRO_11*DENDRO_25;
double DENDRO_53 = -DENDRO_22 + DENDRO_23;
double DENDRO_54 = 0.5*grad_0_gt0[pp];
double DENDRO_55 = grad_2_gt0[pp];
double DENDRO_56 = -0.5*DENDRO_55 + 1.0*grad_0_gt2[pp];
double DENDRO_57 = grad_1_gt0[pp];
double DENDRO_58 = -0.5*DENDRO_57 + 1.0*grad_0_gt1[pp];
double DENDRO_59 = 0.5*gt0[pp];
double DENDRO_60 = DENDRO_25*DENDRO_30;
double DENDRO_61 = grad_1_gt2[pp];
double DENDRO_62 = grad_2_gt1[pp];
double DENDRO_63 = grad_0_gt4[pp];
double DENDRO_64 = DENDRO_61 + DENDRO_62 - DENDRO_63;
double DENDRO_65 = DENDRO_25*gt4[pp];
double DENDRO_66 = 0.5*DENDRO_35;
double DENDRO_67 = 0.5*DENDRO_21;
double DENDRO_68 = 2*chi[pp];
double DENDRO_69 = -DENDRO_61 + DENDRO_62 + DENDRO_63;
double DENDRO_70 = DENDRO_25*gt2[pp];
double DENDRO_71 = 0.5*DENDRO_33;
double DENDRO_72 = -DENDRO_27 + DENDRO_28;
double DENDRO_73 = 2*DENDRO_72;
double DENDRO_74 = DENDRO_61 - DENDRO_62 + DENDRO_63;
double DENDRO_75 = DENDRO_25*gt1[pp];
double DENDRO_76 = pow(DENDRO_11, 2);
double DENDRO_77 = pow(DENDRO_72, 2);
double DENDRO_78 = At4[pp]*DENDRO_72;
double DENDRO_79 = 2*DENDRO_11;
double DENDRO_80 = At1[pp]*DENDRO_53;
double DENDRO_81 = At2[pp]*DENDRO_53;
double DENDRO_82 = pow(DENDRO_24, -2);
double DENDRO_83 = 3*DENDRO_82;
double DENDRO_84 = pow(DENDRO_3, 2);
double DENDRO_85 = DENDRO_11*DENDRO_3;
double DENDRO_86 = At1[pp]*DENDRO_42;
double DENDRO_87 = 2*DENDRO_3;
double DENDRO_88 = DENDRO_3*DENDRO_72;
double DENDRO_89 = At2[pp]*DENDRO_2;
double DENDRO_90 = At4[pp]*DENDRO_2;
double DENDRO_91 = At0[pp]*DENDRO_53;
double DENDRO_92 = DENDRO_11*DENDRO_72;
double DENDRO_93 = At5[pp]*DENDRO_72;
double DENDRO_94 = 6*DENDRO_82;
double DENDRO_95 = At3[pp]*DENDRO_42;

// Dendro: printing variables
//--
K_rhs[pp] = -DENDRO_2*DENDRO_41*(DENDRO_26*(DENDRO_11*DENDRO_13 + DENDRO_19*DENDRO_20 + DENDRO_3*DENDRO_4 + DENDRO_7*DENDRO_9) + DENDRO_34*(DENDRO_11*DENDRO_9 + DENDRO_13*DENDRO_30 + DENDRO_20*DENDRO_32 + DENDRO_29*DENDRO_4) + DENDRO_35*(DENDRO_13*DENDRO_39 - DENDRO_14*(1.0*DENDRO_15 - DENDRO_20*DENDRO_25*DENDRO_40) + DENDRO_37*DENDRO_4 + DENDRO_38*DENDRO_9) - grad2_2_2_alpha[pp]) + DENDRO_38*DENDRO_68*(0.5*DENDRO_34*(DENDRO_11*DENDRO_44 + DENDRO_29*DENDRO_8 + DENDRO_30*DENDRO_64 + DENDRO_32*gt4[pp]) + DENDRO_66*(-DENDRO_14*(DENDRO_16 - DENDRO_40*DENDRO_65) + DENDRO_37*DENDRO_8 + DENDRO_38*DENDRO_44 + DENDRO_39*DENDRO_64) + DENDRO_67*(-DENDRO_14*(DENDRO_15 - DENDRO_18*DENDRO_65) + DENDRO_38*DENDRO_8 + DENDRO_44*DENDRO_51 + DENDRO_52*DENDRO_64) - grad2_1_2_alpha[pp]) - DENDRO_41*DENDRO_42*(DENDRO_21*(-DENDRO_14*(1.0*DENDRO_16 - DENDRO_18*DENDRO_25*DENDRO_49) + DENDRO_38*DENDRO_45 + DENDRO_43*DENDRO_51 + DENDRO_47*DENDRO_52) + DENDRO_34*(DENDRO_11*DENDRO_43 + DENDRO_29*DENDRO_45 + DENDRO_30*DENDRO_47 + DENDRO_32*DENDRO_49) + DENDRO_50*(DENDRO_29*DENDRO_47 + DENDRO_3*DENDRO_43 + DENDRO_36*DENDRO_45 + DENDRO_48*DENDRO_49) - grad2_1_1_alpha[pp]) - DENDRO_41*DENDRO_53*(DENDRO_26*(DENDRO_11*DENDRO_54 + DENDRO_19*DENDRO_59 + DENDRO_3*DENDRO_56 + DENDRO_58*DENDRO_7) + DENDRO_33*(-DENDRO_14*(1.0*DENDRO_17 - DENDRO_25*DENDRO_31*DENDRO_59) + DENDRO_39*DENDRO_56 + DENDRO_52*DENDRO_58 + DENDRO_54*DENDRO_60) + DENDRO_50*(DENDRO_29*DENDRO_54 + DENDRO_3*DENDRO_58 + DENDRO_36*DENDRO_56 + DENDRO_48*DENDRO_59) - grad2_0_0_alpha[pp]) - DENDRO_41*DENDRO_73*(0.5*DENDRO_26*(DENDRO_11*DENDRO_55 + DENDRO_12*DENDRO_3 + DENDRO_19*gt2[pp] + DENDRO_69*DENDRO_7) + DENDRO_66*(DENDRO_12*DENDRO_37 - DENDRO_14*(DENDRO_17 - DENDRO_40*DENDRO_70) + DENDRO_38*DENDRO_69 + DENDRO_39*DENDRO_55) + DENDRO_71*(DENDRO_12*DENDRO_39 - DENDRO_14*(DENDRO_15 - DENDRO_31*DENDRO_70) + DENDRO_52*DENDRO_69 + DENDRO_55*DENDRO_60) - grad2_0_2_alpha[pp]) + DENDRO_52*DENDRO_68*(0.5*DENDRO_50*(DENDRO_29*DENDRO_57 + DENDRO_3*DENDRO_46 + DENDRO_36*DENDRO_74 + DENDRO_48*gt1[pp]) + DENDRO_67*(-DENDRO_14*(DENDRO_17 - DENDRO_18*DENDRO_75) + DENDRO_38*DENDRO_74 + DENDRO_46*DENDRO_51 + DENDRO_52*DENDRO_57) + DENDRO_71*(-DENDRO_14*(DENDRO_16 - DENDRO_31*DENDRO_75) + DENDRO_39*DENDRO_74 + DENDRO_46*DENDRO_52 + DENDRO_57*DENDRO_60) - grad2_0_1_alpha[pp]) + (1.0/3.0)*alpha[pp]*(At0[pp]*DENDRO_83*(At0[pp]*pow(DENDRO_53, 2) + At3[pp]*DENDRO_76 + At5[pp]*DENDRO_77 + DENDRO_73*DENDRO_81 - DENDRO_78*DENDRO_79 - DENDRO_79*DENDRO_80) + At1[pp]*DENDRO_94*(At1[pp]*DENDRO_76 - At2[pp]*DENDRO_92 + At4[pp]*DENDRO_85 - DENDRO_11*DENDRO_91 - DENDRO_11*DENDRO_95 - DENDRO_3*DENDRO_81 - DENDRO_3*DENDRO_93 + DENDRO_42*DENDRO_78 + DENDRO_42*DENDRO_80) + At2[pp]*DENDRO_94*(-At1[pp]*DENDRO_92 + At2[pp]*DENDRO_77 + At3[pp]*DENDRO_85 - DENDRO_11*DENDRO_90 + DENDRO_2*DENDRO_81 + DENDRO_2*DENDRO_93 - DENDRO_3*DENDRO_78 - DENDRO_3*DENDRO_80 + DENDRO_72*DENDRO_91) + At3[pp]*DENDRO_83*(At0[pp]*DENDRO_76 + 2*At2[pp]*DENDRO_85 + At3[pp]*pow(DENDRO_42, 2) - At4[pp]*DENDRO_42*DENDRO_87 + At5[pp]*DENDRO_84 - DENDRO_79*DENDRO_86) + At4[pp]*DENDRO_94*(-At0[pp]*DENDRO_92 + At1[pp]*DENDRO_85 - At2[pp]*DENDRO_88 + At4[pp]*DENDRO_84 - At5[pp]*DENDRO_2*DENDRO_3 - DENDRO_11*DENDRO_89 - DENDRO_3*DENDRO_95 + DENDRO_42*DENDRO_90 + DENDRO_72*DENDRO_86) + At5[pp]*DENDRO_83*(At0[pp]*DENDRO_77 - 2*At1[pp]*DENDRO_88 + At3[pp]*DENDRO_84 + At5[pp]*pow(DENDRO_2, 2) + DENDRO_73*DENDRO_89 - DENDRO_87*DENDRO_90) + pow(K[pp], 2)) + beta0[pp]*agrad_0_K[pp] + beta1[pp]*agrad_1_K[pp] + beta2[pp]*agrad_2_K[pp];
// Dendro: reduced ops: 479
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    bssn::timer::t_rhs_K.stop();

    bssn::timer::t_rhs_Gt.start();

    //#include "CalGt.cpp"
    //CalGt
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.CalGt]
vnames = ['CalGt']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 1716 
// Dendro: printing temp variables
double DENDRO_0 = pow(gt4[pp], 2);
double DENDRO_1 = pow(gt1[pp], 2);
double DENDRO_2 = pow(gt2[pp], 2);
double DENDRO_3 = gt3[pp]*gt5[pp];
double DENDRO_4 = gt1[pp]*gt2[pp];
double DENDRO_5 = pow(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt0[pp] - 2*DENDRO_4*gt4[pp], -2);
double DENDRO_6 = -DENDRO_4 + gt0[pp]*gt4[pp];
double DENDRO_7 = grad_2_gt3[pp];
double DENDRO_8 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
double DENDRO_9 = grad_1_gt5[pp];
double DENDRO_10 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
double DENDRO_11 = -DENDRO_0 + DENDRO_3;
double DENDRO_12 = grad_1_gt2[pp];
double DENDRO_13 = grad_2_gt1[pp];
double DENDRO_14 = grad_0_gt4[pp];
double DENDRO_15 = DENDRO_12 + DENDRO_13 - DENDRO_14;
double DENDRO_16 = grad_0_gt3[pp];
double DENDRO_17 = grad_1_gt0[pp];
double DENDRO_18 = DENDRO_12 - DENDRO_13 + DENDRO_14;
double DENDRO_19 = grad_0_gt5[pp];
double DENDRO_20 = grad_2_gt0[pp];
double DENDRO_21 = -DENDRO_12 + DENDRO_13 + DENDRO_14;
double DENDRO_22 = 1.0*DENDRO_10;
double DENDRO_23 = -DENDRO_1 + gt0[pp]*gt3[pp];
double DENDRO_24 = 0.5*grad_2_gt5[pp];
double DENDRO_25 = 0.5*DENDRO_9 - 1.0*grad_2_gt4[pp];
double DENDRO_26 = 0.5*DENDRO_19 - 1.0*grad_2_gt2[pp];
double DENDRO_27 = -DENDRO_2 + gt0[pp]*gt5[pp];
double DENDRO_28 = 0.5*grad_1_gt3[pp];
double DENDRO_29 = -0.5*DENDRO_7 + 1.0*grad_1_gt4[pp];
double DENDRO_30 = 0.5*DENDRO_16 - 1.0*grad_1_gt1[pp];
double DENDRO_31 = 0.5*grad_0_gt0[pp];
double DENDRO_32 = -0.5*DENDRO_17 + 1.0*grad_0_gt1[pp];
double DENDRO_33 = -0.5*DENDRO_20 + 1.0*grad_0_gt2[pp];

// Dendro: printing variables
//--
CalGt0[pp] = DENDRO_5*(-DENDRO_11*(-DENDRO_10*DENDRO_33 - DENDRO_11*DENDRO_31 + DENDRO_32*DENDRO_8) - DENDRO_22*(-DENDRO_10*DENDRO_19 - DENDRO_11*DENDRO_20 + DENDRO_21*DENDRO_8) - DENDRO_23*(-DENDRO_10*DENDRO_24 + DENDRO_11*DENDRO_26 - DENDRO_25*DENDRO_8) - DENDRO_27*(-DENDRO_10*DENDRO_29 + DENDRO_11*DENDRO_30 + DENDRO_28*DENDRO_8) + DENDRO_6*(-DENDRO_10*DENDRO_9 - DENDRO_11*DENDRO_15 + DENDRO_7*DENDRO_8) + DENDRO_8*(-DENDRO_10*DENDRO_18 - DENDRO_11*DENDRO_17 + DENDRO_16*DENDRO_8));
//--
CalGt1[pp] = DENDRO_5*(-DENDRO_11*(-DENDRO_27*DENDRO_32 + DENDRO_31*DENDRO_8 + DENDRO_33*DENDRO_6) - DENDRO_22*(DENDRO_19*DENDRO_6 + DENDRO_20*DENDRO_8 - DENDRO_21*DENDRO_27) - DENDRO_23*(DENDRO_24*DENDRO_6 + DENDRO_25*DENDRO_27 - DENDRO_26*DENDRO_8) - DENDRO_27*(-DENDRO_27*DENDRO_28 + DENDRO_29*DENDRO_6 - DENDRO_30*DENDRO_8) + DENDRO_6*(DENDRO_15*DENDRO_8 - DENDRO_27*DENDRO_7 + DENDRO_6*DENDRO_9) + DENDRO_8*(-DENDRO_16*DENDRO_27 + DENDRO_17*DENDRO_8 + DENDRO_18*DENDRO_6));
//--
CalGt2[pp] = DENDRO_5*(-DENDRO_11*(-DENDRO_10*DENDRO_31 - DENDRO_23*DENDRO_33 + DENDRO_32*DENDRO_6) - DENDRO_22*(-DENDRO_10*DENDRO_20 - DENDRO_19*DENDRO_23 + DENDRO_21*DENDRO_6) - DENDRO_23*(DENDRO_10*DENDRO_26 - DENDRO_23*DENDRO_24 - DENDRO_25*DENDRO_6) - DENDRO_27*(DENDRO_10*DENDRO_30 - DENDRO_23*DENDRO_29 + DENDRO_28*DENDRO_6) + DENDRO_6*(-DENDRO_10*DENDRO_15 - DENDRO_23*DENDRO_9 + DENDRO_6*DENDRO_7) + DENDRO_8*(-DENDRO_10*DENDRO_17 + DENDRO_16*DENDRO_6 - DENDRO_18*DENDRO_23));
// Dendro: reduced ops: 202
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    

    //#include "Gt_rhs_s1_.cpp"
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs_s1]
vnames = ['Gt_rhs_s1']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 24 
// Dendro: printing temp variables

// Dendro: printing variables
//--
Gt_rhs_s10[pp] = beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp];
//--
Gt_rhs_s11[pp] = beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp];
//--
Gt_rhs_s12[pp] = beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp];
// Dendro: reduced ops: 24
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    
    //#include "Gt_rhs_s2_.cpp"
    
    for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs_s2]
vnames = ['Gt_rhs_s2']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 24 
// Dendro: printing temp variables

// Dendro: printing variables
//--
Gt_rhs_s20[pp] = CalGt0[pp]*grad_0_beta0[pp] + CalGt1[pp]*grad_1_beta0[pp] + CalGt2[pp]*grad_2_beta0[pp];
//--
Gt_rhs_s21[pp] = CalGt0[pp]*grad_0_beta1[pp] + CalGt1[pp]*grad_1_beta1[pp] + CalGt2[pp]*grad_2_beta1[pp];
//--
Gt_rhs_s22[pp] = CalGt0[pp]*grad_0_beta2[pp] + CalGt1[pp]*grad_1_beta2[pp] + CalGt2[pp]*grad_2_beta2[pp];
// Dendro: reduced ops: 24
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    
    //#include "Gt_rhs_s3_.cpp"
    
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs_s3]
vnames = ['Gt_rhs_s3']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 18 
// Dendro: printing temp variables
double DENDRO_0 = grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp];

// Dendro: printing variables
//--
Gt_rhs_s30[pp] = CalGt0[pp]*DENDRO_0;
//--
Gt_rhs_s31[pp] = CalGt1[pp]*DENDRO_0;
//--
Gt_rhs_s32[pp] = CalGt2[pp]*DENDRO_0;
// Dendro: reduced ops: 8
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    
    //#include "Gt_rhs_s4_.cpp"
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs_s4]
vnames = ['Gt_rhs_s4']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 909 
// Dendro: printing temp variables
double DENDRO_0 = pow(gt4[pp], 2);
double DENDRO_1 = pow(gt1[pp], 2);
double DENDRO_2 = pow(gt2[pp], 2);
double DENDRO_3 = gt0[pp]*gt3[pp];
double DENDRO_4 = gt1[pp]*gt2[pp];
double DENDRO_5 = 1.0/(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt5[pp] - 2*DENDRO_4*gt4[pp]);
double DENDRO_6 = -DENDRO_4 + gt0[pp]*gt4[pp];
double DENDRO_7 = grad2_0_2_beta0[pp];
double DENDRO_8 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
double DENDRO_9 = (7.0/3.0)*DENDRO_8;
double DENDRO_10 = grad2_1_2_beta1[pp];
double DENDRO_11 = (1.0/3.0)*DENDRO_8;
double DENDRO_12 = grad2_2_2_beta2[pp];
double DENDRO_13 = grad2_0_1_beta0[pp];
double DENDRO_14 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
double DENDRO_15 = (7.0/3.0)*DENDRO_14;
double DENDRO_16 = grad2_1_1_beta1[pp];
double DENDRO_17 = (1.0/3.0)*DENDRO_14;
double DENDRO_18 = grad2_1_2_beta2[pp];
double DENDRO_19 = -DENDRO_1 + DENDRO_3;
double DENDRO_20 = -DENDRO_2 + gt0[pp]*gt5[pp];
double DENDRO_21 = grad2_0_0_beta0[pp];
double DENDRO_22 = -DENDRO_0 + gt3[pp]*gt5[pp];
double DENDRO_23 = grad2_0_1_beta1[pp];
double DENDRO_24 = (1.0/3.0)*DENDRO_22;
double DENDRO_25 = grad2_0_2_beta2[pp];
double DENDRO_26 = (1.0/3.0)*DENDRO_6;
double DENDRO_27 = (7.0/3.0)*DENDRO_6;
double DENDRO_28 = (1.0/3.0)*DENDRO_20;
double DENDRO_29 = (1.0/3.0)*DENDRO_19;

// Dendro: printing variables
//--
Gt_rhs_s40[pp] = DENDRO_5*(-DENDRO_10*DENDRO_11 - DENDRO_11*DENDRO_12 + DENDRO_13*DENDRO_15 + DENDRO_16*DENDRO_17 + DENDRO_17*DENDRO_18 - DENDRO_19*grad2_2_2_beta0[pp] - DENDRO_20*grad2_1_1_beta0[pp] - 4.0/3.0*DENDRO_21*DENDRO_22 - DENDRO_23*DENDRO_24 - DENDRO_24*DENDRO_25 + 2*DENDRO_6*grad2_1_2_beta0[pp] - DENDRO_7*DENDRO_9);
//--
Gt_rhs_s41[pp] = DENDRO_5*(DENDRO_10*DENDRO_27 + DENDRO_12*DENDRO_26 - DENDRO_13*DENDRO_28 + DENDRO_15*DENDRO_23 - 4.0/3.0*DENDRO_16*DENDRO_20 + DENDRO_17*DENDRO_21 + DENDRO_17*DENDRO_25 - DENDRO_18*DENDRO_28 - DENDRO_19*grad2_2_2_beta1[pp] - DENDRO_22*grad2_0_0_beta1[pp] + DENDRO_26*DENDRO_7 - 2*DENDRO_8*grad2_0_2_beta1[pp]);
//--
Gt_rhs_s42[pp] = DENDRO_5*(-DENDRO_10*DENDRO_29 - DENDRO_11*DENDRO_21 - DENDRO_11*DENDRO_23 - 4.0/3.0*DENDRO_12*DENDRO_19 + DENDRO_13*DENDRO_26 + 2*DENDRO_14*grad2_0_1_beta2[pp] + DENDRO_16*DENDRO_26 + DENDRO_18*DENDRO_27 - DENDRO_20*grad2_1_1_beta2[pp] - DENDRO_22*grad2_0_0_beta2[pp] - DENDRO_25*DENDRO_9 - DENDRO_29*DENDRO_7);
// Dendro: reduced ops: 140
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    
    //#include "Gt_rhs_s5_.cpp"
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs_s5]
vnames = ['Gt_rhs_s5']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 1914 
// Dendro: printing temp variables
double DENDRO_0 = grad_0_alpha[pp];
double DENDRO_1 = gt2[pp]*gt4[pp];
double DENDRO_2 = -DENDRO_1 + gt1[pp]*gt5[pp];
double DENDRO_3 = pow(DENDRO_2, 2);
double DENDRO_4 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
double DENDRO_5 = pow(DENDRO_4, 2);
double DENDRO_6 = gt3[pp]*gt5[pp];
double DENDRO_7 = pow(gt4[pp], 2);
double DENDRO_8 = DENDRO_6 - DENDRO_7;
double DENDRO_9 = DENDRO_2*DENDRO_4;
double DENDRO_10 = 2*At4[pp];
double DENDRO_11 = At1[pp]*DENDRO_8;
double DENDRO_12 = 2*DENDRO_2;
double DENDRO_13 = At2[pp]*DENDRO_8;
double DENDRO_14 = 2*DENDRO_4;
double DENDRO_15 = grad_2_alpha[pp];
double DENDRO_16 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
double DENDRO_17 = At3[pp]*DENDRO_2;
double DENDRO_18 = At0[pp]*DENDRO_8;
double DENDRO_19 = At4[pp]*DENDRO_4;
double DENDRO_20 = pow(gt1[pp], 2);
double DENDRO_21 = -DENDRO_20 + gt0[pp]*gt3[pp];
double DENDRO_22 = At5[pp]*DENDRO_4;
double DENDRO_23 = At4[pp]*DENDRO_2;
double DENDRO_24 = -At1[pp]*DENDRO_9 + At2[pp]*DENDRO_5 - DENDRO_11*DENDRO_16 + DENDRO_13*DENDRO_21 + DENDRO_16*DENDRO_17 - DENDRO_16*DENDRO_19 + DENDRO_18*DENDRO_4 + DENDRO_21*DENDRO_22 - DENDRO_21*DENDRO_23;
double DENDRO_25 = grad_1_alpha[pp];
double DENDRO_26 = pow(gt2[pp], 2);
double DENDRO_27 = -DENDRO_26 + gt0[pp]*gt5[pp];
double DENDRO_28 = -At1[pp]*DENDRO_3 + At2[pp]*DENDRO_9 - DENDRO_11*DENDRO_27 + DENDRO_13*DENDRO_16 + DENDRO_16*DENDRO_22 - DENDRO_16*DENDRO_23 + DENDRO_17*DENDRO_27 + DENDRO_18*DENDRO_2 - DENDRO_19*DENDRO_27;
double DENDRO_29 = 2/pow(-2*DENDRO_1*gt1[pp] + DENDRO_20*gt5[pp] + DENDRO_26*gt3[pp] - DENDRO_6*gt0[pp] + DENDRO_7*gt0[pp], 2);
double DENDRO_30 = pow(DENDRO_16, 2);
double DENDRO_31 = At2[pp]*DENDRO_2;
double DENDRO_32 = At1[pp]*DENDRO_27;
double DENDRO_33 = DENDRO_16*DENDRO_27;
double DENDRO_34 = DENDRO_16*DENDRO_4;
double DENDRO_35 = DENDRO_16*DENDRO_21;
double DENDRO_36 = At0[pp]*DENDRO_9 - At1[pp]*DENDRO_16*DENDRO_2 + At2[pp]*DENDRO_34 + At3[pp]*DENDRO_33 - At4[pp]*DENDRO_21*DENDRO_27 - At4[pp]*DENDRO_30 + At5[pp]*DENDRO_35 + DENDRO_21*DENDRO_31 - DENDRO_32*DENDRO_4;

// Dendro: printing variables
//--
Gt_rhs_s50[pp] = DENDRO_29*(DENDRO_0*(At0[pp]*pow(DENDRO_8, 2) + At3[pp]*DENDRO_3 + At5[pp]*DENDRO_5 - DENDRO_10*DENDRO_9 - DENDRO_11*DENDRO_12 + DENDRO_13*DENDRO_14) + DENDRO_15*DENDRO_24 - DENDRO_25*DENDRO_28);
//--
Gt_rhs_s51[pp] = DENDRO_29*(-DENDRO_0*DENDRO_28 - DENDRO_15*DENDRO_36 + DENDRO_25*(At0[pp]*DENDRO_3 + At3[pp]*pow(DENDRO_27, 2) + At5[pp]*DENDRO_30 - DENDRO_10*DENDRO_33 - DENDRO_12*DENDRO_32 + 2*DENDRO_16*DENDRO_31));
//--
Gt_rhs_s52[pp] = DENDRO_29*(DENDRO_0*DENDRO_24 + DENDRO_15*(At0[pp]*DENDRO_5 - 2*At1[pp]*DENDRO_34 + At2[pp]*DENDRO_14*DENDRO_21 + At3[pp]*DENDRO_30 + At5[pp]*pow(DENDRO_21, 2) - DENDRO_10*DENDRO_35) - DENDRO_25*DENDRO_36);
// Dendro: reduced ops: 162
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    
    //#include "Gt_rhs_s6_.cpp"
    
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs_s6]
vnames = ['Gt_rhs_s6']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 4812 
// Dendro: printing temp variables
double DENDRO_0 = grad_1_gt3[pp];
double DENDRO_1 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
double DENDRO_2 = 0.5*DENDRO_1;
double DENDRO_3 = gt1[pp]*gt4[pp];
double DENDRO_4 = DENDRO_3 - gt2[pp]*gt3[pp];
double DENDRO_5 = grad_2_gt3[pp];
double DENDRO_6 = -0.5*DENDRO_5 + 1.0*grad_1_gt4[pp];
double DENDRO_7 = gt3[pp]*gt5[pp];
double DENDRO_8 = pow(gt4[pp], 2);
double DENDRO_9 = DENDRO_7 - DENDRO_8;
double DENDRO_10 = grad_0_gt3[pp];
double DENDRO_11 = 0.5*DENDRO_10 - 1.0*grad_1_gt1[pp];
double DENDRO_12 = pow(DENDRO_1, 2);
double DENDRO_13 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
double DENDRO_14 = pow(DENDRO_13, 2);
double DENDRO_15 = pow(gt2[pp], 2);
double DENDRO_16 = -DENDRO_15 + gt0[pp]*gt5[pp];
double DENDRO_17 = At2[pp]*DENDRO_1;
double DENDRO_18 = 2*DENDRO_13;
double DENDRO_19 = At1[pp]*DENDRO_16;
double DENDRO_20 = 2*DENDRO_1;
double DENDRO_21 = At4[pp]*DENDRO_16;
double DENDRO_22 = At0[pp]*DENDRO_12 + At3[pp]*pow(DENDRO_16, 2) + At5[pp]*DENDRO_14 + DENDRO_17*DENDRO_18 - DENDRO_18*DENDRO_21 - DENDRO_19*DENDRO_20;
double DENDRO_23 = grad_0_gt0[pp];
double DENDRO_24 = 0.5*DENDRO_23;
double DENDRO_25 = grad_2_gt0[pp];
double DENDRO_26 = -0.5*DENDRO_25 + 1.0*grad_0_gt2[pp];
double DENDRO_27 = grad_1_gt0[pp];
double DENDRO_28 = -0.5*DENDRO_27 + 1.0*grad_0_gt1[pp];
double DENDRO_29 = pow(DENDRO_4, 2);
double DENDRO_30 = At4[pp]*DENDRO_4;
double DENDRO_31 = At1[pp]*DENDRO_1;
double DENDRO_32 = DENDRO_4*DENDRO_9;
double DENDRO_33 = At0[pp]*pow(DENDRO_9, 2) + 2*At2[pp]*DENDRO_32 + At3[pp]*DENDRO_12 + At5[pp]*DENDRO_29 - DENDRO_20*DENDRO_30 - 2*DENDRO_31*DENDRO_9;
double DENDRO_34 = 0.5*grad_2_gt5[pp];
double DENDRO_35 = grad_1_gt5[pp];
double DENDRO_36 = 0.5*DENDRO_35 - 1.0*grad_2_gt4[pp];
double DENDRO_37 = grad_0_gt5[pp];
double DENDRO_38 = 0.5*DENDRO_37 - 1.0*grad_2_gt2[pp];
double DENDRO_39 = pow(gt1[pp], 2);
double DENDRO_40 = -DENDRO_39 + gt0[pp]*gt3[pp];
double DENDRO_41 = At1[pp]*DENDRO_13;
double DENDRO_42 = 2*DENDRO_4;
double DENDRO_43 = At2[pp]*DENDRO_40;
double DENDRO_44 = At4[pp]*DENDRO_40;
double DENDRO_45 = At0[pp]*DENDRO_29 + At3[pp]*DENDRO_14 + At5[pp]*pow(DENDRO_40, 2) - DENDRO_18*DENDRO_44 - DENDRO_41*DENDRO_42 + DENDRO_42*DENDRO_43;
double DENDRO_46 = grad_0_gt4[pp];
double DENDRO_47 = grad_1_gt2[pp];
double DENDRO_48 = grad_2_gt1[pp];
double DENDRO_49 = DENDRO_46 + DENDRO_47 - DENDRO_48;
double DENDRO_50 = At5[pp]*DENDRO_13;
double DENDRO_51 = At0[pp]*DENDRO_1;
double DENDRO_52 = At2[pp]*DENDRO_13;
double DENDRO_53 = At3[pp]*DENDRO_16;
double DENDRO_54 = DENDRO_1*DENDRO_13;
double DENDRO_55 = -At1[pp]*DENDRO_12 - At4[pp]*DENDRO_54 + DENDRO_1*DENDRO_53 + DENDRO_17*DENDRO_4 - DENDRO_19*DENDRO_9 - DENDRO_21*DENDRO_4 + DENDRO_4*DENDRO_50 + DENDRO_51*DENDRO_9 + DENDRO_52*DENDRO_9;
double DENDRO_56 = -DENDRO_46 + DENDRO_47 + DENDRO_48;
double DENDRO_57 = -At1[pp]*DENDRO_54 - At4[pp]*DENDRO_14 + DENDRO_13*DENDRO_53 + DENDRO_17*DENDRO_40 - DENDRO_19*DENDRO_4 - DENDRO_21*DENDRO_40 + DENDRO_4*DENDRO_51 + DENDRO_4*DENDRO_52 + DENDRO_40*DENDRO_50;
double DENDRO_58 = DENDRO_46 - DENDRO_47 + DENDRO_48;
double DENDRO_59 = At0[pp]*DENDRO_32 + At2[pp]*DENDRO_29 + At3[pp]*DENDRO_54 + At5[pp]*DENDRO_4*DENDRO_40 - DENDRO_1*DENDRO_44 - DENDRO_13*DENDRO_30 - DENDRO_31*DENDRO_4 - DENDRO_41*DENDRO_9 + DENDRO_43*DENDRO_9;
double DENDRO_60 = 1.0*DENDRO_59;
double DENDRO_61 = 2*alpha[pp]/pow(DENDRO_15*gt3[pp] - 2*DENDRO_3*gt2[pp] + DENDRO_39*gt5[pp] - DENDRO_7*gt0[pp] + DENDRO_8*gt0[pp], 3);
double DENDRO_62 = 0.5*DENDRO_0;

// Dendro: printing variables
//--
Gt_rhs_s60[pp] = DENDRO_61*(DENDRO_22*(DENDRO_0*DENDRO_2 + DENDRO_11*DENDRO_9 - DENDRO_4*DENDRO_6) - DENDRO_33*(-DENDRO_1*DENDRO_28 + DENDRO_24*DENDRO_9 + DENDRO_26*DENDRO_4) - DENDRO_45*(DENDRO_1*DENDRO_36 + DENDRO_34*DENDRO_4 - DENDRO_38*DENDRO_9) + DENDRO_55*(-DENDRO_1*DENDRO_10 + DENDRO_27*DENDRO_9 + DENDRO_4*DENDRO_49) + DENDRO_57*(-DENDRO_1*DENDRO_5 + DENDRO_35*DENDRO_4 + DENDRO_56*DENDRO_9) - DENDRO_60*(-DENDRO_1*DENDRO_58 + DENDRO_25*DENDRO_9 + DENDRO_37*DENDRO_4));
//--
Gt_rhs_s61[pp] = DENDRO_61*(-DENDRO_22*(DENDRO_1*DENDRO_11 - DENDRO_13*DENDRO_6 + DENDRO_16*DENDRO_62) + DENDRO_33*(DENDRO_13*DENDRO_26 - DENDRO_16*DENDRO_28 + DENDRO_2*DENDRO_23) + DENDRO_45*(-DENDRO_1*DENDRO_38 + DENDRO_13*DENDRO_34 + DENDRO_16*DENDRO_36) - 1.0*DENDRO_55*(DENDRO_1*DENDRO_27 - DENDRO_10*DENDRO_16 + DENDRO_13*DENDRO_49) - 1.0*DENDRO_57*(DENDRO_1*DENDRO_56 + DENDRO_13*DENDRO_35 - DENDRO_16*DENDRO_5) + DENDRO_59*(DENDRO_1*DENDRO_25 + DENDRO_13*DENDRO_37 - DENDRO_16*DENDRO_58));
//--
Gt_rhs_s62[pp] = DENDRO_61*(DENDRO_22*(DENDRO_11*DENDRO_4 + DENDRO_13*DENDRO_62 - DENDRO_40*DENDRO_6) - DENDRO_33*(-DENDRO_13*DENDRO_28 + DENDRO_24*DENDRO_4 + DENDRO_26*DENDRO_40) - DENDRO_45*(DENDRO_13*DENDRO_36 + DENDRO_34*DENDRO_40 - DENDRO_38*DENDRO_4) + DENDRO_55*(-DENDRO_10*DENDRO_13 + DENDRO_27*DENDRO_4 + DENDRO_40*DENDRO_49) + DENDRO_57*(-DENDRO_13*DENDRO_5 + DENDRO_35*DENDRO_40 + DENDRO_4*DENDRO_56) - DENDRO_60*(-DENDRO_13*DENDRO_58 + DENDRO_25*DENDRO_4 + DENDRO_37*DENDRO_40));
// Dendro: reduced ops: 316
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    
    //#include "Gt_rhs_s7_.cpp"
    
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs_s7]
vnames = ['Gt_rhs_s7']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 2121 
// Dendro: printing temp variables
double DENDRO_0 = gt1[pp]*gt4[pp];
double DENDRO_1 = DENDRO_0 - gt2[pp]*gt3[pp];
double DENDRO_2 = (4.0/3.0)*grad_2_K[pp];
double DENDRO_3 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
double DENDRO_4 = (4.0/3.0)*grad_1_K[pp];
double DENDRO_5 = gt3[pp]*gt5[pp];
double DENDRO_6 = pow(gt4[pp], 2);
double DENDRO_7 = DENDRO_5 - DENDRO_6;
double DENDRO_8 = (4.0/3.0)*grad_0_K[pp];
double DENDRO_9 = pow(DENDRO_3, 2);
double DENDRO_10 = pow(DENDRO_1, 2);
double DENDRO_11 = DENDRO_1*DENDRO_3;
double DENDRO_12 = 2*At4[pp];
double DENDRO_13 = At1[pp]*DENDRO_7;
double DENDRO_14 = 2*DENDRO_3;
double DENDRO_15 = At2[pp]*DENDRO_7;
double DENDRO_16 = 2*DENDRO_1;
double DENDRO_17 = grad_0_chi[pp];
double DENDRO_18 = pow(gt1[pp], 2);
double DENDRO_19 = pow(gt2[pp], 2);
double DENDRO_20 = 1.0/(-2*DENDRO_0*gt2[pp] + DENDRO_18*gt5[pp] + DENDRO_19*gt3[pp] - DENDRO_5*gt0[pp] + DENDRO_6*gt0[pp]);
double DENDRO_21 = 3*DENDRO_20/chi[pp];
double DENDRO_22 = DENDRO_17*DENDRO_21;
double DENDRO_23 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
double DENDRO_24 = At3[pp]*DENDRO_3;
double DENDRO_25 = At0[pp]*DENDRO_7;
double DENDRO_26 = At4[pp]*DENDRO_1;
double DENDRO_27 = -DENDRO_18 + gt0[pp]*gt3[pp];
double DENDRO_28 = At5[pp]*DENDRO_1;
double DENDRO_29 = At4[pp]*DENDRO_3;
double DENDRO_30 = -At1[pp]*DENDRO_11 + At2[pp]*DENDRO_10 + DENDRO_1*DENDRO_25 - DENDRO_13*DENDRO_23 + DENDRO_15*DENDRO_27 + DENDRO_23*DENDRO_24 - DENDRO_23*DENDRO_26 + DENDRO_27*DENDRO_28 - DENDRO_27*DENDRO_29;
double DENDRO_31 = DENDRO_21*grad_2_chi[pp];
double DENDRO_32 = grad_1_chi[pp];
double DENDRO_33 = -DENDRO_19 + gt0[pp]*gt5[pp];
double DENDRO_34 = DENDRO_21*(-At1[pp]*DENDRO_9 + At2[pp]*DENDRO_11 - DENDRO_13*DENDRO_33 + DENDRO_15*DENDRO_23 + DENDRO_23*DENDRO_28 - DENDRO_23*DENDRO_29 + DENDRO_24*DENDRO_33 + DENDRO_25*DENDRO_3 - DENDRO_26*DENDRO_33);
double DENDRO_35 = DENDRO_20*alpha[pp];
double DENDRO_36 = pow(DENDRO_23, 2);
double DENDRO_37 = At2[pp]*DENDRO_3;
double DENDRO_38 = At1[pp]*DENDRO_33;
double DENDRO_39 = DENDRO_23*DENDRO_33;
double DENDRO_40 = DENDRO_21*DENDRO_32;
double DENDRO_41 = DENDRO_1*DENDRO_23;
double DENDRO_42 = DENDRO_23*DENDRO_27;
double DENDRO_43 = At0[pp]*DENDRO_11 - At1[pp]*DENDRO_23*DENDRO_3 + At2[pp]*DENDRO_41 + At3[pp]*DENDRO_39 - At4[pp]*DENDRO_27*DENDRO_33 - At4[pp]*DENDRO_36 + At5[pp]*DENDRO_42 - DENDRO_1*DENDRO_38 + DENDRO_27*DENDRO_37;

// Dendro: printing variables
//--
Gt_rhs_s70[pp] = DENDRO_35*(-DENDRO_1*DENDRO_2 + DENDRO_22*(At0[pp]*pow(DENDRO_7, 2) + At3[pp]*DENDRO_9 + At5[pp]*DENDRO_10 - DENDRO_11*DENDRO_12 - DENDRO_13*DENDRO_14 + DENDRO_15*DENDRO_16) + DENDRO_3*DENDRO_4 + DENDRO_30*DENDRO_31 - DENDRO_32*DENDRO_34 - DENDRO_7*DENDRO_8);
//--
Gt_rhs_s71[pp] = DENDRO_35*(-DENDRO_17*DENDRO_34 + DENDRO_2*DENDRO_23 + DENDRO_3*DENDRO_8 - DENDRO_31*DENDRO_43 - DENDRO_33*DENDRO_4 + DENDRO_40*(At0[pp]*DENDRO_9 + At3[pp]*pow(DENDRO_33, 2) + At5[pp]*DENDRO_36 - DENDRO_12*DENDRO_39 - DENDRO_14*DENDRO_38 + 2*DENDRO_23*DENDRO_37));
//--
Gt_rhs_s72[pp] = DENDRO_35*(-DENDRO_1*DENDRO_8 - DENDRO_2*DENDRO_27 + DENDRO_22*DENDRO_30 + DENDRO_23*DENDRO_4 + DENDRO_31*(At0[pp]*DENDRO_10 - 2*At1[pp]*DENDRO_41 + At2[pp]*DENDRO_16*DENDRO_27 + At3[pp]*DENDRO_36 + At5[pp]*pow(DENDRO_27, 2) - DENDRO_12*DENDRO_42) - DENDRO_40*DENDRO_43);
// Dendro: reduced ops: 195
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    //#include "Gt_rhs.cpp"

for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.Gt_rhs]
vnames = ['Gt_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 24 
// Dendro: printing temp variables

// Dendro: printing variables
//--
Gt_rhs0[pp] = Gt_rhs_s1_0[pp] - Gt_rhs_s2_0[pp] + (2.0/3.0)*Gt_rhs_s3_0[pp] + Gt_rhs_s4_0[pp] - Gt_rhs_s5_0[pp] + Gt_rhs_s6_0[pp] - Gt_rhs_s7_0[pp];
//--
Gt_rhs1[pp] = Gt_rhs_s1_1[pp] - Gt_rhs_s2_1[pp] + (2.0/3.0)*Gt_rhs_s3_1[pp] + Gt_rhs_s4_1[pp] - Gt_rhs_s5_1[pp] + Gt_rhs_s6_1[pp] - Gt_rhs_s7_1[pp];
//--
Gt_rhs2[pp] = Gt_rhs_s1_2[pp] - Gt_rhs_s2_2[pp] + (2.0/3.0)*Gt_rhs_s3_2[pp] + Gt_rhs_s4_2[pp] - Gt_rhs_s5_2[pp] + Gt_rhs_s6_2[pp] - Gt_rhs_s7_2[pp];
// Dendro: reduced ops: 24
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
    
    
    bssn::timer::t_rhs_Gt.stop();

    bssn::timer::t_rhs_B.start();
#ifdef USE_ETA_FUNC
//#include "B_rhs_eta_func.cpp"
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

outs = [bssn_stages.B_rhs]
vnames = ['B_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 516 
// Dendro: printing temp variables
double DENDRO_0 = pow(gt4[pp], 2);
double DENDRO_1 = pow(gt1[pp], 2);
double DENDRO_2 = pow(gt2[pp], 2);
double DENDRO_3 = gt3[pp]*gt5[pp];
double DENDRO_4 = gt1[pp]*gt4[pp];
double DENDRO_5 = grad_2_chi[pp];
double DENDRO_6 = grad_1_chi[pp];
double DENDRO_7 = grad_0_chi[pp];
double DENDRO_8 = 2*DENDRO_5;
double DENDRO_9 = BSSN_ETA_R0*sqrt((-pow(DENDRO_5, 2)*(-DENDRO_1 + gt0[pp]*gt3[pp]) - pow(DENDRO_6, 2)*(-DENDRO_2 + gt0[pp]*gt5[pp]) + 2*DENDRO_6*DENDRO_7*(gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp]) + DENDRO_6*DENDRO_8*(gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp]) - pow(DENDRO_7, 2)*(-DENDRO_0 + DENDRO_3) - DENDRO_7*DENDRO_8*(DENDRO_4 - gt2[pp]*gt3[pp]))/(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt0[pp] - 2*DENDRO_4*gt2[pp]))*pow(-pow(chi[pp], BSSN_ETA_POWER[0]) + 1, -BSSN_ETA_POWER[1]);

// Dendro: printing variables
//--
B_rhs0[pp] = -B0[pp]*DENDRO_9 + Gt_rhs0[pp] + lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] + beta2[pp]*agrad_2_B0[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp]);
//--
B_rhs1[pp] = -B1[pp]*DENDRO_9 + Gt_rhs1[pp] + lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] + beta2[pp]*agrad_2_B1[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp]);
//--
B_rhs2[pp] = -B2[pp]*DENDRO_9 + Gt_rhs2[pp] + lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] + beta2[pp]*agrad_2_B2[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp]);
// Dendro: reduced ops: 124
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
#else

//#include "B_rhs_eta_const.cpp"
for (unsigned int k = 3; k < nz-3; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = 3; j < ny-3; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = 3; i < nx-3; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         r_coord = sqrt(x*x + y*y + z*z);
         eta=ETA_CONST;
         if (r_coord >= ETA_R0) {
          eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
         }
/*[[[cog

import dendro
import bssn_stages

bssn_stages.eta_func=bssn_stages.eta
bssn_stages.B_rhs = [bssn_stages.Gt_rhs[i] - bssn_stages.eta_func * bssn_stages.B[i] +
         bssn_stages.l3 * dendro.vec_j_ad_j(bssn_stages.b, bssn_stages.B[i]) -
         bssn_stages.l4 * dendro.vec_j_ad_j(bssn_stages.b, bssn_stages.Gt[i])
         for i in dendro.e_i]

outs = [bssn_stages.B_rhs]
vnames = ['B_rhs']
dendro.generate_cpu(outs, vnames, '[pp]')

]]]*/
// Dendro: {{{ 
// Dendro: original ops: 90 
// Dendro: printing temp variables

// Dendro: printing variables
//--
B_rhs0[pp] = -B0[pp]*eta + Gt_rhs_s1_0[pp] - Gt_rhs_s2_0[pp] + (2.0/3.0)*Gt_rhs_s3_0[pp] + Gt_rhs_s4_0[pp] - Gt_rhs_s5_0[pp] + Gt_rhs_s6_0[pp] - Gt_rhs_s7_0[pp] + lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] + beta2[pp]*agrad_2_B0[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt0[pp] + beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp]);
//--
B_rhs1[pp] = -B1[pp]*eta + Gt_rhs_s1_1[pp] - Gt_rhs_s2_1[pp] + (2.0/3.0)*Gt_rhs_s3_1[pp] + Gt_rhs_s4_1[pp] - Gt_rhs_s5_1[pp] + Gt_rhs_s6_1[pp] - Gt_rhs_s7_1[pp] + lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] + beta2[pp]*agrad_2_B1[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt1[pp] + beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp]);
//--
B_rhs2[pp] = -B2[pp]*eta + Gt_rhs_s1_2[pp] - Gt_rhs_s2_2[pp] + (2.0/3.0)*Gt_rhs_s3_2[pp] + Gt_rhs_s4_2[pp] - Gt_rhs_s5_2[pp] + Gt_rhs_s6_2[pp] - Gt_rhs_s7_2[pp] + lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] + beta2[pp]*agrad_2_B2[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt2[pp] + beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp]);
// Dendro: reduced ops: 90
// Dendro: }}} 
//[[[end]]]
      }
     }
    }
#endif

    bssn::timer::t_rhs_B.stop();
    
    bssn::timer::t_rhs.stop();


    delete [] CalGt0;
    delete [] CalGt1;
    delete [] CalGt2;

    delete [] Gt_rhs_s1_0;
    delete [] Gt_rhs_s1_1;
    delete [] Gt_rhs_s1_2;

    delete [] Gt_rhs_s2_0;
    delete [] Gt_rhs_s2_1;
    delete [] Gt_rhs_s2_2;

    delete [] Gt_rhs_s3_0;
    delete [] Gt_rhs_s3_1;
    delete [] Gt_rhs_s3_2;

    delete [] Gt_rhs_s4_0;
    delete [] Gt_rhs_s4_1;
    delete [] Gt_rhs_s4_2;

    delete [] Gt_rhs_s5_0;
    delete [] Gt_rhs_s5_1;
    delete [] Gt_rhs_s5_2;

    delete [] Gt_rhs_s6_0;
    delete [] Gt_rhs_s6_1;
    delete [] Gt_rhs_s6_2;

    delete [] Gt_rhs_s7_0;
    delete [] Gt_rhs_s7_1;
    delete [] Gt_rhs_s7_2;



    if (bflag != 0) {

        bssn::timer::t_bdyc.start();

#ifdef BSSN_KERR_SCHILD_TEST

        freeze_bcs(a_rhs, sz, bflag);
        freeze_bcs(chi_rhs, sz, bflag);
        freeze_bcs(K_rhs, sz, bflag);
        freeze_bcs(b_rhs0, sz, bflag);
        freeze_bcs(b_rhs1, sz, bflag);
        freeze_bcs(b_rhs2, sz, bflag);
        freeze_bcs(Gt_rhs0, sz, bflag);
        freeze_bcs(Gt_rhs1, sz, bflag);
        freeze_bcs(Gt_rhs2, sz, bflag);
        freeze_bcs(B_rhs0, sz, bflag);
        freeze_bcs(B_rhs1, sz, bflag);
        freeze_bcs(B_rhs2, sz, bflag);
        freeze_bcs(At_rhs00, sz, bflag);
        freeze_bcs(At_rhs01, sz, bflag);
        freeze_bcs(At_rhs02, sz, bflag);
        freeze_bcs(At_rhs11, sz, bflag);
        freeze_bcs(At_rhs12, sz, bflag);
        freeze_bcs(At_rhs22, sz, bflag);
        freeze_bcs(gt_rhs00, sz, bflag);
        freeze_bcs(gt_rhs01, sz, bflag);
        freeze_bcs(gt_rhs02, sz, bflag);
        freeze_bcs(gt_rhs11, sz, bflag);
        freeze_bcs(gt_rhs12, sz, bflag);
        freeze_bcs(gt_rhs22, sz, bflag);

#else
        bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax,
                 1.0, 0.0, sz, bflag);

        bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                 2.0, 0.0, sz, bflag);

        bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                 1.0, 1.0, sz, bflag);
#endif
        bssn::timer::t_bdyc.stop();
    }


    bssn::timer::t_deriv.start();
    
// ko derivs

/*[[[cog
import cog
import bssnDerivs as bssnDerivs
    
for var in bssnDerivs.KO:
    cog.outl("\t ko_deriv_x(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[0] + var ,var))
    cog.outl("\t ko_deriv_y(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[1] + var ,var))
    cog.outl("\t ko_deriv_z(%s, %s, hx, sz, bflag);" %(bssnDerivs.PREFIX_D[2] + var ,var))    
    
]]]*/
ko_deriv_x(grad_0_gt0, gt0, hx, sz, bflag);
ko_deriv_y(grad_1_gt0, gt0, hx, sz, bflag);
ko_deriv_z(grad_2_gt0, gt0, hx, sz, bflag);
ko_deriv_x(grad_0_gt1, gt1, hx, sz, bflag);
ko_deriv_y(grad_1_gt1, gt1, hx, sz, bflag);
ko_deriv_z(grad_2_gt1, gt1, hx, sz, bflag);
ko_deriv_x(grad_0_gt2, gt2, hx, sz, bflag);
ko_deriv_y(grad_1_gt2, gt2, hx, sz, bflag);
ko_deriv_z(grad_2_gt2, gt2, hx, sz, bflag);
ko_deriv_x(grad_0_gt3, gt3, hx, sz, bflag);
ko_deriv_y(grad_1_gt3, gt3, hx, sz, bflag);
ko_deriv_z(grad_2_gt3, gt3, hx, sz, bflag);
ko_deriv_x(grad_0_gt4, gt4, hx, sz, bflag);
ko_deriv_y(grad_1_gt4, gt4, hx, sz, bflag);
ko_deriv_z(grad_2_gt4, gt4, hx, sz, bflag);
ko_deriv_x(grad_0_gt5, gt5, hx, sz, bflag);
ko_deriv_y(grad_1_gt5, gt5, hx, sz, bflag);
ko_deriv_z(grad_2_gt5, gt5, hx, sz, bflag);
ko_deriv_x(grad_0_At0, At0, hx, sz, bflag);
ko_deriv_y(grad_1_At0, At0, hx, sz, bflag);
ko_deriv_z(grad_2_At0, At0, hx, sz, bflag);
ko_deriv_x(grad_0_At1, At1, hx, sz, bflag);
ko_deriv_y(grad_1_At1, At1, hx, sz, bflag);
ko_deriv_z(grad_2_At1, At1, hx, sz, bflag);
ko_deriv_x(grad_0_At2, At2, hx, sz, bflag);
ko_deriv_y(grad_1_At2, At2, hx, sz, bflag);
ko_deriv_z(grad_2_At2, At2, hx, sz, bflag);
ko_deriv_x(grad_0_At3, At3, hx, sz, bflag);
ko_deriv_y(grad_1_At3, At3, hx, sz, bflag);
ko_deriv_z(grad_2_At3, At3, hx, sz, bflag);
ko_deriv_x(grad_0_At4, At4, hx, sz, bflag);
ko_deriv_y(grad_1_At4, At4, hx, sz, bflag);
ko_deriv_z(grad_2_At4, At4, hx, sz, bflag);
ko_deriv_x(grad_0_At5, At5, hx, sz, bflag);
ko_deriv_y(grad_1_At5, At5, hx, sz, bflag);
ko_deriv_z(grad_2_At5, At5, hx, sz, bflag);
ko_deriv_x(grad_0_alpha, alpha, hx, sz, bflag);
ko_deriv_y(grad_1_alpha, alpha, hx, sz, bflag);
ko_deriv_z(grad_2_alpha, alpha, hx, sz, bflag);
ko_deriv_x(grad_0_beta0, beta0, hx, sz, bflag);
ko_deriv_y(grad_1_beta0, beta0, hx, sz, bflag);
ko_deriv_z(grad_2_beta0, beta0, hx, sz, bflag);
ko_deriv_x(grad_0_beta1, beta1, hx, sz, bflag);
ko_deriv_y(grad_1_beta1, beta1, hx, sz, bflag);
ko_deriv_z(grad_2_beta1, beta1, hx, sz, bflag);
ko_deriv_x(grad_0_beta2, beta2, hx, sz, bflag);
ko_deriv_y(grad_1_beta2, beta2, hx, sz, bflag);
ko_deriv_z(grad_2_beta2, beta2, hx, sz, bflag);
ko_deriv_x(grad_0_chi, chi, hx, sz, bflag);
ko_deriv_y(grad_1_chi, chi, hx, sz, bflag);
ko_deriv_z(grad_2_chi, chi, hx, sz, bflag);
ko_deriv_x(grad_0_Gt0, Gt0, hx, sz, bflag);
ko_deriv_y(grad_1_Gt0, Gt0, hx, sz, bflag);
ko_deriv_z(grad_2_Gt0, Gt0, hx, sz, bflag);
ko_deriv_x(grad_0_Gt1, Gt1, hx, sz, bflag);
ko_deriv_y(grad_1_Gt1, Gt1, hx, sz, bflag);
ko_deriv_z(grad_2_Gt1, Gt1, hx, sz, bflag);
ko_deriv_x(grad_0_Gt2, Gt2, hx, sz, bflag);
ko_deriv_y(grad_1_Gt2, Gt2, hx, sz, bflag);
ko_deriv_z(grad_2_Gt2, Gt2, hx, sz, bflag);
ko_deriv_x(grad_0_K, K, hx, sz, bflag);
ko_deriv_y(grad_1_K, K, hx, sz, bflag);
ko_deriv_z(grad_2_K, K, hx, sz, bflag);
ko_deriv_x(grad_0_B0, B0, hx, sz, bflag);
ko_deriv_y(grad_1_B0, B0, hx, sz, bflag);
ko_deriv_z(grad_2_B0, B0, hx, sz, bflag);
ko_deriv_x(grad_0_B1, B1, hx, sz, bflag);
ko_deriv_y(grad_1_B1, B1, hx, sz, bflag);
ko_deriv_z(grad_2_B1, B1, hx, sz, bflag);
ko_deriv_x(grad_0_B2, B2, hx, sz, bflag);
ko_deriv_y(grad_1_B2, B2, hx, sz, bflag);
ko_deriv_z(grad_2_B2, B2, hx, sz, bflag);
//[[[end]]]


    bssn::timer::t_deriv.stop();

    bssn::timer::t_rhs.start();

    const  double sigma = KO_DISS_SIGMA;


    for (unsigned int k = 3; k < nz-3; k++) {
        for (unsigned int j = 3; j < ny-3; j++) {
            for (unsigned int i = 3; i < nx-3; i++) {
                pp = i + nx*(j + ny*k);

                a_rhs[pp]  += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
                b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
                b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
                b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

                gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
            }
        }
    }

    bssn::timer::t_rhs.stop();

    //bssn::writeBLockToBinary((const double**)unzipVarsRHS,offset,pmin,pmax,bxMin,bxMax,sz,15,1.0,"af_ko");

    bssn::timer::t_deriv.start();
    //dealloc derivative variables. 
/*[[[cog
import cog
import bssnDerivs as bssnDerivs

for deriv in bssnDerivs.FUNC_D_I:
    cog.outl("\t free(%s);" %(deriv))
    
for deriv in bssnDerivs.FUNC_D_IJ:
    cog.outl("\t free(%s);" %(deriv))

for deriv in bssnDerivs.FUNC_AD_I:
    cog.outl("\t free(%s);" %(deriv))   
    
]]]*/
free(grad_0_alpha);
free(grad_1_alpha);
free(grad_2_alpha);
free(grad_0_beta0);
free(grad_1_beta0);
free(grad_2_beta0);
free(grad_0_beta1);
free(grad_1_beta1);
free(grad_2_beta1);
free(grad_0_beta2);
free(grad_1_beta2);
free(grad_2_beta2);
free(grad_0_B0);
free(grad_1_B0);
free(grad_2_B0);
free(grad_0_B1);
free(grad_1_B1);
free(grad_2_B1);
free(grad_0_B2);
free(grad_1_B2);
free(grad_2_B2);
free(grad_0_chi);
free(grad_1_chi);
free(grad_2_chi);
free(grad_0_Gt0);
free(grad_1_Gt0);
free(grad_2_Gt0);
free(grad_0_Gt1);
free(grad_1_Gt1);
free(grad_2_Gt1);
free(grad_0_Gt2);
free(grad_1_Gt2);
free(grad_2_Gt2);
free(grad_0_K);
free(grad_1_K);
free(grad_2_K);
free(grad_0_gt0);
free(grad_1_gt0);
free(grad_2_gt0);
free(grad_0_gt1);
free(grad_1_gt1);
free(grad_2_gt1);
free(grad_0_gt2);
free(grad_1_gt2);
free(grad_2_gt2);
free(grad_0_gt3);
free(grad_1_gt3);
free(grad_2_gt3);
free(grad_0_gt4);
free(grad_1_gt4);
free(grad_2_gt4);
free(grad_0_gt5);
free(grad_1_gt5);
free(grad_2_gt5);
free(grad_0_At0);
free(grad_1_At0);
free(grad_2_At0);
free(grad_0_At1);
free(grad_1_At1);
free(grad_2_At1);
free(grad_0_At2);
free(grad_1_At2);
free(grad_2_At2);
free(grad_0_At3);
free(grad_1_At3);
free(grad_2_At3);
free(grad_0_At4);
free(grad_1_At4);
free(grad_2_At4);
free(grad_0_At5);
free(grad_1_At5);
free(grad_2_At5);
free(grad2_0_0_gt0);
free(grad2_0_1_gt0);
free(grad2_0_2_gt0);
free(grad2_1_1_gt0);
free(grad2_1_2_gt0);
free(grad2_2_2_gt0);
free(grad2_0_0_gt1);
free(grad2_0_1_gt1);
free(grad2_0_2_gt1);
free(grad2_1_1_gt1);
free(grad2_1_2_gt1);
free(grad2_2_2_gt1);
free(grad2_0_0_gt2);
free(grad2_0_1_gt2);
free(grad2_0_2_gt2);
free(grad2_1_1_gt2);
free(grad2_1_2_gt2);
free(grad2_2_2_gt2);
free(grad2_0_0_gt3);
free(grad2_0_1_gt3);
free(grad2_0_2_gt3);
free(grad2_1_1_gt3);
free(grad2_1_2_gt3);
free(grad2_2_2_gt3);
free(grad2_0_0_gt4);
free(grad2_0_1_gt4);
free(grad2_0_2_gt4);
free(grad2_1_1_gt4);
free(grad2_1_2_gt4);
free(grad2_2_2_gt4);
free(grad2_0_0_gt5);
free(grad2_0_1_gt5);
free(grad2_0_2_gt5);
free(grad2_1_1_gt5);
free(grad2_1_2_gt5);
free(grad2_2_2_gt5);
free(grad2_0_0_chi);
free(grad2_0_1_chi);
free(grad2_0_2_chi);
free(grad2_1_1_chi);
free(grad2_1_2_chi);
free(grad2_2_2_chi);
free(grad2_0_0_alpha);
free(grad2_0_1_alpha);
free(grad2_0_2_alpha);
free(grad2_1_1_alpha);
free(grad2_1_2_alpha);
free(grad2_2_2_alpha);
free(grad2_0_0_beta0);
free(grad2_0_1_beta0);
free(grad2_0_2_beta0);
free(grad2_1_1_beta0);
free(grad2_1_2_beta0);
free(grad2_2_2_beta0);
free(grad2_0_0_beta1);
free(grad2_0_1_beta1);
free(grad2_0_2_beta1);
free(grad2_1_1_beta1);
free(grad2_1_2_beta1);
free(grad2_2_2_beta1);
free(grad2_0_0_beta2);
free(grad2_0_1_beta2);
free(grad2_0_2_beta2);
free(grad2_1_1_beta2);
free(grad2_1_2_beta2);
free(grad2_2_2_beta2);
free(agrad_0_gt0);
free(agrad_1_gt0);
free(agrad_2_gt0);
free(agrad_0_gt1);
free(agrad_1_gt1);
free(agrad_2_gt1);
free(agrad_0_gt2);
free(agrad_1_gt2);
free(agrad_2_gt2);
free(agrad_0_gt3);
free(agrad_1_gt3);
free(agrad_2_gt3);
free(agrad_0_gt4);
free(agrad_1_gt4);
free(agrad_2_gt4);
free(agrad_0_gt5);
free(agrad_1_gt5);
free(agrad_2_gt5);
free(agrad_0_At0);
free(agrad_1_At0);
free(agrad_2_At0);
free(agrad_0_At1);
free(agrad_1_At1);
free(agrad_2_At1);
free(agrad_0_At2);
free(agrad_1_At2);
free(agrad_2_At2);
free(agrad_0_At3);
free(agrad_1_At3);
free(agrad_2_At3);
free(agrad_0_At4);
free(agrad_1_At4);
free(agrad_2_At4);
free(agrad_0_At5);
free(agrad_1_At5);
free(agrad_2_At5);
free(agrad_0_alpha);
free(agrad_1_alpha);
free(agrad_2_alpha);
free(agrad_0_beta0);
free(agrad_1_beta0);
free(agrad_2_beta0);
free(agrad_0_beta1);
free(agrad_1_beta1);
free(agrad_2_beta1);
free(agrad_0_beta2);
free(agrad_1_beta2);
free(agrad_2_beta2);
free(agrad_0_chi);
free(agrad_1_chi);
free(agrad_2_chi);
free(agrad_0_Gt0);
free(agrad_1_Gt0);
free(agrad_2_Gt0);
free(agrad_0_Gt1);
free(agrad_1_Gt1);
free(agrad_2_Gt1);
free(agrad_0_Gt2);
free(agrad_1_Gt2);
free(agrad_2_Gt2);
free(agrad_0_K);
free(agrad_1_K);
free(agrad_2_K);
free(agrad_0_B0);
free(agrad_1_B0);
free(agrad_2_B0);
free(agrad_0_B1);
free(agrad_1_B1);
free(agrad_2_B1);
free(agrad_0_B2);
free(agrad_1_B2);
free(agrad_2_B2);
//[[[end]]]
    bssn::timer::t_deriv.stop();

#if 0
    for (unsigned int m = 0; m < 24; m++) {
    std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
  }
#endif



}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void bssn_bcs(double *f_rhs, const double *f,
              const double *dxf, const double *dyf, const double *dzf,
              const double *pmin, const double *pmax,
              const double f_falloff, const double f_asymptotic,
              const unsigned int *sz, const unsigned int &bflag)
{

  const unsigned int nx = sz[0];
  const unsigned int ny = sz[1];
  const unsigned int nz = sz[2];

  double hx = (pmax[0] - pmin[0]) / (nx - 1);
  double hy = (pmax[1] - pmin[1]) / (ny - 1);
  double hz = (pmax[2] - pmin[2]) / (nz - 1);

  unsigned int ib = 3;
  unsigned int jb = 3;
  unsigned int kb = 3;
  unsigned int ie = sz[0]-3;
  unsigned int je = sz[1]-3;
  unsigned int ke = sz[2]-3;

  double x,y,z;
  unsigned int pp;
  double inv_r;

  //std::cout<<"boundary bssnrhs: size [ "<<nx<<", "<<ny<<", "<<nz<<" ]"<<std::endl;
  //std::cout<<"boundary bssnrhs: pmin [ "<<pmin[0]<<", "<<pmin[1]<<", "<<pmin[2]<<" ]"<<std::endl;
  //std::cout<<"boundary bssnrhs: pmax [ "<<pmax[0]<<", "<<pmax[1]<<", "<<pmax[2]<<" ]"<<std::endl;

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    double x = pmin[0] + ib*hx;
    for (unsigned int k = kb; k < ke; k++) {
       z = pmin[2] + k*hz;
      for (unsigned int j = jb; j < je; j++) {
         y = pmin[1] + j*hy;
         pp = IDX(ib,j,k);
         inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        f_rhs[pp] = -  inv_r * (
                         x * dxf[pp]
                       + y * dyf[pp]
                       + z * dzf[pp]
                       + f_falloff * (   f[pp] - f_asymptotic ) );

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
     x = pmin[0] + (ie-1)*hx;
    for (unsigned int k = kb; k < ke; k++) {
       z = pmin[2] + k*hz;
      for (unsigned int j = jb; j < je; j++) {
         y = pmin[1] + j*hy;
         pp = IDX((ie-1),j,k);
         inv_r = 1.0 / sqrt(x*x + y*y + z*z);

        f_rhs[pp] = -  inv_r * (
                         x * dxf[pp]
                       + y * dyf[pp]
                       + z * dzf[pp]
                       + f_falloff * (   f[pp] - f_asymptotic ) );

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
     y = pmin[1] + jb*hy;
    for (unsigned int k = kb; k < ke; k++) {
       z = pmin[2] + k*hz;
      for (unsigned int i = ib; i < ie; i++) {
         x = pmin[0] + i*hx;
         inv_r = 1.0 / sqrt(x*x + y*y + z*z);
         pp = IDX(i,jb,k);

        f_rhs[pp] = -  inv_r * (
                         x * dxf[pp]
                       + y * dyf[pp]
                       + z * dzf[pp]
                       + f_falloff * (   f[pp] - f_asymptotic ) );

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
     y = pmin[1] + (je-1)*hy;
    for (unsigned int k = kb; k < ke; k++) {
       z = pmin[2] + k*hz;
      for (unsigned int i = ib; i < ie; i++) {
         x = pmin[0] + i*hx;
         inv_r = 1.0 / sqrt(x*x + y*y + z*z);
         pp = IDX(i,(je-1),k);

        f_rhs[pp] = -  inv_r * (
                         x * dxf[pp]
                       + y * dyf[pp]
                       + z * dzf[pp]
                       + f_falloff * (   f[pp] - f_asymptotic ) );

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
     z = pmin[2] + kb*hz;
    for (unsigned int j = jb; j < je; j++) {
       y = pmin[1] + j*hy;
      for (unsigned int i = ib; i < ie; i++) {
         x = pmin[0] + i*hx;
         inv_r = 1.0 / sqrt(x*x + y*y + z*z);
         pp = IDX(i,j,kb);

        f_rhs[pp] = - inv_r * (
                         x * dxf[pp]
                       + y * dyf[pp]
                       + z * dzf[pp]
                       + f_falloff * (   f[pp] - f_asymptotic ) );

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    z = pmin[2] + (ke-1)*hz;
    for (unsigned int j = jb; j < je; j++) {
      y = pmin[1] + j*hy;
      for (unsigned int i = ib; i < ie; i++) {
        x = pmin[0] + i*hx;
        inv_r = 1.0 / sqrt(x*x + y*y + z*z);
        pp = IDX(i,j,(ke-1));

        f_rhs[pp] = - inv_r * (
                         x * dxf[pp]
                       + y * dyf[pp]
                       + z * dzf[pp]
                       + f_falloff * (   f[pp] - f_asymptotic ) );

      }
    }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void freeze_bcs(double *f_rhs, const unsigned int *sz, const unsigned int &bflag)
{

  const unsigned int nx = sz[0];
  const unsigned int ny = sz[1];
  const unsigned int nz = sz[2];

  unsigned int ib = 3;
  unsigned int jb = 3;
  unsigned int kb = 3;
  unsigned int ie = sz[0]-3;
  unsigned int je = sz[1]-3;
  unsigned int ke = sz[2]-3;

  unsigned int pp;

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (unsigned int k = kb; k < ke; k++) {
      for (unsigned int j = jb; j < je; j++) {
        pp = IDX(ib,j,k);
        f_rhs[pp] = 0.0;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (unsigned int k = kb; k < ke; k++) {
      for (unsigned int j = jb; j < je; j++) {
        pp = IDX((ie-1),j,k);
        f_rhs[pp] = 0.0;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (unsigned int k = kb; k < ke; k++) {
      for (unsigned int i = ib; i < ie; i++) {
        pp = IDX(i,jb,k);
        f_rhs[pp] = 0.0;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (unsigned int k = kb; k < ke; k++) {
      for (unsigned int i = ib; i < ie; i++) {
        pp = IDX(i,(je-1),k);
        f_rhs[pp] = 0.0;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (unsigned int j = jb; j < je; j++) {
      for (unsigned int i = ib; i < ie; i++) {
        pp = IDX(i,j,kb);
        f_rhs[pp] = 0.0;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (unsigned int j = jb; j < je; j++) {
      for (unsigned int i = ib; i < ie; i++) {
        pp = IDX(i,j,(ke-1));
        f_rhs[pp] = 0.0;
      }
    }
  }

}

#if 0
/*--------------------------------------------------------------
 * Kerr-Schild data
 *--------------------------------------------------------------*/

void ks_initial_data(double x, double y, double z, double *u)
{

  u[VAR::U_ALPHA] = 0.0;
  u[VAR::U_BETA0] = 0.0;
  u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
  u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

  u[VAR::U_B0] = 0.0;
  u[VAR::U_B1] = 0.0;
  u[VAR::U_B2] = 0.0;

  u[VAR::U_GT0] = Gamt_1;
  u[VAR::U_GT1] = Gamt_2;
  u[VAR::U_GT2] = Gamt_3;

  u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));

  u[VAR::U_SYMGT0] = 1.00+0.2*sin(x+z)*cos(y);
  u[VAR::U_SYMGT3] = 1.00+0.2*cos(y)*cos(z+ x);
  u[VAR::U_SYMGT5] = 1.00 / ( u[VAR::U_SYMGT0] + u[VAR::U_SYMGT3]);
  u[VAR::U_SYMGT1] = 0.7*cos(x*x + y*y);
  u[VAR::U_SYMGT2] = 0.3*sin(z)*cos(x);
  u[VAR::U_SYMGT4] = -0.5*sin(x*x)*cos(y)*cos(z);

  u[VAR::U_K] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
           +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
           +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
               *exp(-4.0*cos(x)*sin(y))*cos(z);

  u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
  u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
  u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
  u[VAR::U_SYMAT3] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
  u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
  u[VAR::U_SYMAT5] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));

}
#endif
/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
