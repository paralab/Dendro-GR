#include "physcon.h"
#include "gr.h"

using namespace bssn;


/*----------------------------------------------------------------------
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void physical_constraints(double **uZipConVars, const double **uZipVars,
                       const unsigned int& offset,
                       const double *pmin, const double *pmax,
                       const unsigned int *sz, const unsigned int& bflag)
{

  const unsigned int nx = sz[0];
  const unsigned int ny = sz[1];
  const unsigned int nz = sz[2];
  const unsigned int n = nx * ny * nz;

  const double hx = (pmax[0] - pmin[0]) / (nx - 1);
  const double hy = (pmax[1] - pmin[1]) / (ny - 1);
  const double hz = (pmax[2] - pmin[2]) / (nz - 1);

  double * const ham = &uZipConVars[VAR_CONSTRAINT::C_HAM][offset];
  double * const mom0 = &uZipConVars[VAR_CONSTRAINT::C_MOM0][offset];
  double * const mom1 = &uZipConVars[VAR_CONSTRAINT::C_MOM1][offset];
  double * const mom2 = &uZipConVars[VAR_CONSTRAINT::C_MOM2][offset];
  double * const psi4_real = &uZipConVars[VAR_CONSTRAINT::C_PSI4_REAL][offset];
  double * const psi4_img = &uZipConVars[VAR_CONSTRAINT::C_PSI4_IMG][offset];

  const double * const alpha = &uZipVars[VAR::U_ALPHA][offset];
  const double * const chi = &uZipVars[VAR::U_CHI][offset];
  const double * const K = &uZipVars[VAR::U_K][offset];
  const double * const gt0 = &uZipVars[VAR::U_SYMGT0][offset];
  const double * const gt1 = &uZipVars[VAR::U_SYMGT1][offset];
  const double * const gt2 = &uZipVars[VAR::U_SYMGT2][offset];
  const double * const gt3 = &uZipVars[VAR::U_SYMGT3][offset];
  const double * const gt4 = &uZipVars[VAR::U_SYMGT4][offset];
  const double * const gt5 = &uZipVars[VAR::U_SYMGT5][offset];
  const double * const beta0 = &uZipVars[VAR::U_BETA0][offset];
  const double * const beta1 = &uZipVars[VAR::U_BETA1][offset];
  const double * const beta2 = &uZipVars[VAR::U_BETA2][offset];
  const double * const At0 = &uZipVars[VAR::U_SYMAT0][offset];
  const double * const At1 = &uZipVars[VAR::U_SYMAT1][offset];
  const double * const At2 = &uZipVars[VAR::U_SYMAT2][offset];
  const double * const At3 = &uZipVars[VAR::U_SYMAT3][offset];
  const double * const At4 = &uZipVars[VAR::U_SYMAT4][offset];
  const double * const At5 = &uZipVars[VAR::U_SYMAT5][offset];
  const double * const Gt0 = &uZipVars[VAR::U_GT0][offset];
  const double * const Gt1 = &uZipVars[VAR::U_GT1][offset];
  const double * const Gt2 = &uZipVars[VAR::U_GT2][offset];
  const double * const B0 = &uZipVars[VAR::U_B0][offset];
  const double * const B1 = &uZipVars[VAR::U_B1][offset];
  const double * const B2 = &uZipVars[VAR::U_B2][offset];
  const unsigned int PW=bssn::BSSN_PADDING_WIDTH;

  const unsigned int BLK_SZ=n;
  double * const deriv_base = bssn::BSSN_DERIV_WORKSPACE;
  #include "bssnrhs_evar_derivs.h"
  #include "constraint_derivs.h"

  // enforce hamiltonian and momentum constraints
  for (unsigned int k = PW; k < nz - PW; k++) {
    for (unsigned int j = PW; j < ny - PW; j++) {
      #ifdef BSSN_ENABLE_AVX
        #ifdef __INTEL_COMPILER
        #pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
        #pragma ivdep
        #endif
      #endif
      for (unsigned int i = PW; i < nx - PW; i++) {
        const double x = pmin[0] + i*hx;
        const double y = pmin[1] + j*hy;
        const double z = pmin[2] + k*hz;
        const unsigned int pp = i + nx * (j + ny * k);
        
        #include "physconeqs.cpp"

        if(fabs(x) <=1e-7 && fabs(y)<=1e-7)
        {
          psi4_real[pp]=0.0;
          psi4_img[pp]=0.0;
        }

      }
    }
  }


}
