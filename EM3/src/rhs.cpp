#include "rhs.h"

using namespace std;
using namespace em3;

/*----------------------------------------------------------------------;
 *
 * RHS for NOT the non-linear sigma model
 *
 *----------------------------------------------------------------------*/
void em3rhs(double **unzipVarsRHS, const double **uZipVars,
             const unsigned int& offset,
             const double *pmin, const double *pmax, const unsigned int *sz,
             const unsigned int& bflag)
{

  const double *E0 = &uZipVars[VAR::U_E0][offset];
  const double *E1 = &uZipVars[VAR::U_E1][offset];
  const double *E2 = &uZipVars[VAR::U_E2][offset];
  const double *B0 = &uZipVars[VAR::U_B0][offset];
  const double *B1 = &uZipVars[VAR::U_B1][offset];
  const double *B2 = &uZipVars[VAR::U_B2][offset];

  double *E_rhs0 = &unzipVarsRHS[VAR::U_E0][offset];
  double *E_rhs1 = &unzipVarsRHS[VAR::U_E1][offset];
  double *E_rhs2 = &unzipVarsRHS[VAR::U_E2][offset];
  double *B_rhs0 = &unzipVarsRHS[VAR::U_B0][offset];
  double *B_rhs1 = &unzipVarsRHS[VAR::U_B1][offset];
  double *B_rhs2 = &unzipVarsRHS[VAR::U_B2][offset];

  mem::memory_pool<DendroScalar>* __mem_pool = &EM3_MEM_POOL;
  
  const unsigned int nx = sz[0];
  const unsigned int ny = sz[1];
  const unsigned int nz = sz[2];

  double PI = 3.14159265358979323846 ;  

  double hx = (pmax[0] - pmin[0]) / (nx - 1);
  double hy = (pmax[1] - pmin[1]) / (ny - 1);
  double hz = (pmax[2] - pmin[2]) / (nz - 1);

  int idx[3];

  unsigned int n = sz[0]*sz[1]*sz[2];

  em3::timer::t_deriv.start();

  double *rho_e = __mem_pool->allocate(n);

  double *J0 = __mem_pool->allocate(n);
  double *J1 = __mem_pool->allocate(n);
  double *J2 = __mem_pool->allocate(n);

  for(unsigned int m=0; m < n ; m++)
  {
    
    rho_e[m] = 0.0; 
    
    J0[m] = 0.0; 
    J1[m] = 0.0;
    J2[m] = 0.0;
    
  }


  double *grad_0_E0 = __mem_pool->allocate(n);
  double *grad_1_E0 = __mem_pool->allocate(n);
  double *grad_2_E0 = __mem_pool->allocate(n);

  double *grad_0_E1 = __mem_pool->allocate(n);
  double *grad_1_E1 = __mem_pool->allocate(n);
  double *grad_2_E1 = __mem_pool->allocate(n);

  double *grad_0_E2 = __mem_pool->allocate(n);
  double *grad_1_E2 = __mem_pool->allocate(n);
  double *grad_2_E2 = __mem_pool->allocate(n);

  double *grad_0_B0 = __mem_pool->allocate(n);
  double *grad_1_B0 = __mem_pool->allocate(n);
  double *grad_2_B0 = __mem_pool->allocate(n);

  double *grad_0_B1 = __mem_pool->allocate(n);
  double *grad_1_B1 = __mem_pool->allocate(n);
  double *grad_2_B1 = __mem_pool->allocate(n);

  double *grad_0_B2 = __mem_pool->allocate(n);
  double *grad_1_B2 = __mem_pool->allocate(n);
  double *grad_2_B2 = __mem_pool->allocate(n);

  deriv_x(grad_0_E0, E0, hx, sz, bflag);
  deriv_y(grad_1_E0, E0, hy, sz, bflag);
  deriv_z(grad_2_E0, E0, hz, sz, bflag);

  deriv_x(grad_0_E1, E1, hx, sz, bflag);
  deriv_y(grad_1_E1, E1, hy, sz, bflag);
  deriv_z(grad_2_E1, E1, hz, sz, bflag);

  deriv_x(grad_0_E2, E2, hx, sz, bflag);
  deriv_y(grad_1_E2, E2, hy, sz, bflag);
  deriv_z(grad_2_E2, E2, hz, sz, bflag);

  deriv_x(grad_0_B0, B0, hx, sz, bflag);
  deriv_y(grad_1_B0, B0, hy, sz, bflag);
  deriv_z(grad_2_B0, B0, hz, sz, bflag);

  deriv_x(grad_0_B1, B1, hx, sz, bflag);
  deriv_y(grad_1_B1, B1, hy, sz, bflag);
  deriv_z(grad_2_B1, B1, hz, sz, bflag);

  deriv_x(grad_0_B2, B2, hx, sz, bflag);
  deriv_y(grad_1_B2, B2, hy, sz, bflag);
  deriv_z(grad_2_B2, B2, hz, sz, bflag);



  em3::timer::t_deriv.stop();

  register double x;
  register double y;
  register double z;
  register unsigned int pp;

  double r;
  double eta;
  const unsigned int PW = em3::EM3_PADDING_WIDTH;

  //cout << "begin loop" << endl;
  for (unsigned int k = PW ; k < nz-PW; k++) {
      z = pmin[2] + k*hz;

    for (unsigned int j = PW; j < ny-PW; j++) {
       y = pmin[1] + j*hy;

      for (unsigned int i = PW; i < nx-PW; i++) {
         x = pmin[0] + i*hx;
         pp = i + nx*(j + ny*k);
         //r= sqrt(x*x + y*y + z*z);

          em3::timer::t_rhs.start();
            #include "em3_eqs.cpp"
          em3::timer::t_rhs.stop();
         
      }
    }
  }

    if (bflag != 0) {

      em3::timer::t_bdyc.start();

      
      em3_bcs(E_rhs0, E0, grad_0_E0, grad_1_E0, grad_2_E0, pmin, pmax,
              2.0, 0.0, sz, bflag);
      em3_bcs(E_rhs1, E1, grad_0_E1, grad_1_E1, grad_2_E1, pmin, pmax,
              2.0, 0.0, sz, bflag);
      em3_bcs(E_rhs2, E2, grad_0_E2, grad_1_E2, grad_2_E2, pmin, pmax,
              2.0, 0.0, sz, bflag);
      
      em3_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax,
              2.0, 0.0, sz, bflag);
      em3_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax,
              2.0, 0.0, sz, bflag);
      em3_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax,
              2.0, 0.0, sz, bflag);
      
      em3::timer::t_bdyc.stop();
    }


    em3::timer::t_deriv.start();

    ko_deriv_x(grad_0_E0, E0, hx, sz, bflag);
    ko_deriv_y(grad_1_E0, E0, hy, sz, bflag);
    ko_deriv_z(grad_2_E0, E0, hz, sz, bflag);

    ko_deriv_x(grad_0_E1, E1, hx, sz, bflag);
    ko_deriv_y(grad_1_E1, E1, hy, sz, bflag);
    ko_deriv_z(grad_2_E1, E1, hz, sz, bflag);

    ko_deriv_x(grad_0_E2, E2, hx, sz, bflag);
    ko_deriv_y(grad_1_E2, E2, hy, sz, bflag);
    ko_deriv_z(grad_2_E2, E2, hz, sz, bflag);

    ko_deriv_x(grad_0_B0, B0, hx, sz, bflag);
    ko_deriv_y(grad_1_B0, B0, hy, sz, bflag);
    ko_deriv_z(grad_2_B0, B0, hz, sz, bflag);

    ko_deriv_x(grad_0_B1, B1, hx, sz, bflag);
    ko_deriv_y(grad_1_B1, B1, hy, sz, bflag);
    ko_deriv_z(grad_2_B1, B1, hz, sz, bflag);

    ko_deriv_x(grad_0_B2, B2, hx, sz, bflag);
    ko_deriv_y(grad_1_B2, B2, hy, sz, bflag);
    ko_deriv_z(grad_2_B2, B2, hz, sz, bflag);

    em3::timer::t_deriv.stop();


    em3::timer::t_rhs.start();

    const  double sigma = KO_DISS_SIGMA;

      for (unsigned int k = PW; k < nz-PW; k++) {
        for (unsigned int j = PW; j < ny-PW; j++) {
          for (unsigned int i = PW; i < nx-PW; i++) {
            pp = i + nx*(j + ny*k);
    
            E_rhs0[pp] += sigma*(grad_0_E0[pp] + grad_1_E0[pp] + grad_2_E0[pp]);
            E_rhs1[pp] += sigma*(grad_0_E1[pp] + grad_1_E1[pp] + grad_2_E1[pp]);
            E_rhs2[pp] += sigma*(grad_0_E2[pp] + grad_1_E2[pp] + grad_2_E2[pp]);
            
            B_rhs0[pp] += sigma*(grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
            B_rhs1[pp] += sigma*(grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
            B_rhs2[pp] += sigma*(grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
    
          }
        }
      }

    em3::timer::t_rhs.stop();


    em3::timer::t_deriv.start();

    __mem_pool->free(grad_0_E0);
    __mem_pool->free(grad_1_E0);
    __mem_pool->free(grad_2_E0);

    __mem_pool->free(grad_0_E1);
    __mem_pool->free(grad_1_E1);
    __mem_pool->free(grad_2_E1);

    __mem_pool->free(grad_0_E2);
    __mem_pool->free(grad_1_E2);
    __mem_pool->free(grad_2_E2);

    __mem_pool->free(grad_0_B0);
    __mem_pool->free(grad_1_B0);
    __mem_pool->free(grad_2_B0);

    __mem_pool->free(grad_0_B1);
    __mem_pool->free(grad_1_B1);
    __mem_pool->free(grad_2_B1);

    __mem_pool->free(grad_0_B2);
    __mem_pool->free(grad_1_B2);
    __mem_pool->free(grad_2_B2);

    __mem_pool->free(J0);
    __mem_pool->free(J1);
    __mem_pool->free(J2);
    
    __mem_pool->free(rho_e);
    
    em3::timer::t_deriv.stop();

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
void em3_bcs(double *f_rhs, const double *f,
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

  const unsigned int PW = em3::EM3_PADDING_WIDTH;

  unsigned int ib = PW;
  unsigned int jb = PW;
  unsigned int kb = PW;
  unsigned int ie = sz[0]-PW;
  unsigned int je = sz[1]-PW;
  unsigned int ke = sz[2]-PW;

  double x,y,z;
  unsigned int pp;
  double inv_r;

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    double x = pmin[0] + ib*hx;
    for (unsigned int k = kb; k < ke; k++) {
       z = pmin[2] + k*hz;
      for (unsigned int j = jb; j < je; j++) {
         y = pmin[1] + j*hy;
         pp = IDX(ib,j,k);
         inv_r = 1.0 / sqrt(x*x + y*y + z*z);

#ifdef EM3_DIRICHLET_BDY
        f_rhs[pp] = 0.0;
#else
        f_rhs[pp] = -  inv_r * (   x * dxf[pp]
                                 + y * dyf[pp]
                                 + z * dzf[pp]
                                 + f_falloff * ( f[pp] - f_asymptotic ) 
                               );
#endif

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

#ifdef EM3_DIRICHLET_BDY
        f_rhs[pp] = 0.0;
#else
        f_rhs[pp] = - inv_r * (   x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * ( f[pp] - f_asymptotic ) 
                              ) ;
#endif

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

#ifdef EM3_DIRICHLET_BDY
        f_rhs[pp] = 0.0;
#else
        f_rhs[pp] = - inv_r * (   x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * ( f[pp] - f_asymptotic ) 
                              ) ;
#endif

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

#ifdef EM3_DIRICHLET_BDY
        f_rhs[pp] = 0.0;
#else
        f_rhs[pp] = - inv_r * (   x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * ( f[pp] - f_asymptotic ) 
                              ) ;
#endif

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

#ifdef EM3_DIRICHLET_BDY
        f_rhs[pp] = 0.0;
#else
        f_rhs[pp] = - inv_r * (   x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * ( f[pp] - f_asymptotic ) 
                              ) ;
#endif

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
        
#ifdef EM3_DIRICHLET_BDY
        f_rhs[pp] = 0.0;
#else
        f_rhs[pp] = - inv_r * (   x * dxf[pp]
                                + y * dyf[pp]
                                + z * dzf[pp]
                                + f_falloff * ( f[pp] - f_asymptotic ) 
                              ) ;
#endif
        

      }
    }
  }

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
