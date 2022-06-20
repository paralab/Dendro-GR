#include "physcon.h"
#include "em2.h"

using namespace em2;


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

  double PI = 3.14159265358979323846 ;

  double hx = (pmax[0] - pmin[0]) / (nx - 1);
  double hy = (pmax[1] - pmin[1]) / (ny - 1);
  double hz = (pmax[2] - pmin[2]) / (nz - 1);

  double *divA = &uZipConVars[VAR_CONSTRAINT::C_DIVA][offset];
  double *divE = &uZipConVars[VAR_CONSTRAINT::C_DIVE][offset];

  const double *A0 = &uZipVars[VAR::U_A0][offset];
  const double *A1 = &uZipVars[VAR::U_A1][offset];
  const double *A2 = &uZipVars[VAR::U_A2][offset];
  
  const double *E0 = &uZipVars[VAR::U_E0][offset];
  const double *E1 = &uZipVars[VAR::U_E1][offset];
  const double *E2 = &uZipVars[VAR::U_E2][offset];

  const double *Gamma = &uZipVars[VAR::U_GAMMA][offset];

  double *grad_0_A0 = new double [n] ; 
  double *grad_1_A1 = new double [n] ; 
  double *grad_2_A2 = new double [n] ; 
  
  double *grad_0_E0 = new double [n] ; 
  double *grad_1_E1 = new double [n] ; 
  double *grad_2_E2 = new double [n] ; 

  deriv_x(grad_0_A0, A0, hx, sz, bflag); 
  deriv_y(grad_1_A1, A1, hy, sz, bflag); 
  deriv_z(grad_2_A2, A2, hz, sz, bflag); 

  deriv_x(grad_0_E0, E0, hx, sz, bflag); 
  deriv_y(grad_1_E1, E1, hy, sz, bflag); 
  deriv_z(grad_2_E2, E2, hz, sz, bflag); 


  double *rho_e = new double [n] ; 
  const unsigned int PW = em2::EM2_PADDING_WIDTH; 

  for ( unsigned int m=0; m<n ; m++ )
  { 
      rho_e[m] = 0.0 ; 
  }


  for (unsigned int k = PW; k < nz - PW; k++) {
    for (unsigned int j = PW; j < ny - PW; j++) {
      for (unsigned int i = PW; i < nx - PW; i++) {
        
        const double x = pmin[0] + i*hx;
        const double y = pmin[1] + j*hy;
        const double z = pmin[2] + k*hz;

        unsigned int pp = i + nx * (j + ny * k);
        #include "physcon_eqs.cpp"

      }
    }
  }

  delete [] grad_0_A0 ; 
  delete [] grad_1_A1 ; 
  delete [] grad_2_A2 ; 
  delete [] grad_0_E0 ; 
  delete [] grad_1_E1 ; 
  delete [] grad_2_E2 ; 
  delete [] rho_e ; 


}
