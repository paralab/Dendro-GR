#include <cmath>
#include <iostream>
#include "derivs.h"


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv22_xx(double *const DxDxu, const double *const u, const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx_sqrd = 1.0/(dx*dx);
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);

        DxDxu[pp] = (          u[pp-1] 
                        -2.0 * u[pp]
                         +     u[pp+1]
                        
                    ) * idx_sqrd;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
    
        DxDxu[IDX(3,j,k)] = (         u[IDX(3,j,k)]
                              - 2.0 * u[IDX(4,j,k)]
                              +       u[IDX(5,j,k)]
                            ) * idx_sqrd;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {

        DxDxu[IDX(ie-1,j,k)] = (       u[IDX(ie-3,j,k)]
                               - 2.0 * u[IDX(ie-2,j,k)]
                               +       u[IDX(ie-1,j,k)]
                              ) * idx_sqrd;
      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif
}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv22_yy(double *const DyDyu, const double *const u, const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy_sqrd = 1.0/(dy*dy);
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        DyDyu[pp] = (  u[pp-nx] 
                    - 2.0 * u[pp] 
                    + u[pp+nx] 
                 ) * idy_sqrd;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        
        DyDyu[IDX(i,3,k)] = (         u[IDX(i,3,k)]
                           - 2.0 * u[IDX(i,4,k)]
                           +       u[IDX(i,5,k)]
                        ) * idy_sqrd;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {

    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        DyDyu[IDX(i,je-1,k)] = (      u[IDX(i,je-3,k)]
                           - 2.0 * u[IDX(i,je-2,k)]
                           +       u[IDX(i,je-1,k)]
                          ) * idy_sqrd;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}




/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv22_zz(double *const DzDzu, const double *const u, const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz_sqrd = 1.0/(dz*dz);
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        DzDzu[pp] = (  u[pp-n] - 2.0 * u[pp] + u[pp+n] ) * idz_sqrd;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        

        DzDzu[IDX(i,j,3)] = (         u[IDX(i,j,3)]
                           - 2.0 * u[IDX(i,j,4)]
                           +       u[IDX(i,j,5)]
                        ) * idz_sqrd;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {

        DzDzu[IDX(i,j,ke-1)] = (      u[IDX(i,j,ke-3)]
                           - 2.0 * u[IDX(i,j,ke-2)]
                           +       u[IDX(i,j,ke-1)]
                          ) * idz_sqrd;
      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}




/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_x(double * const  Dxu, const double * const  u,
               const double dx, const unsigned int *sz, unsigned bflag)
{

  const double idx = 1.0/dx;
  const double idx_by_2 = 0.5 * idx;
  const double idx_by_12 = idx / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 1;
  const int kb = 1;
  const int ie = sz[0]-3;
  const int je = sz[1]-1;
  const int ke = sz[2]-1;
  
  const int n=1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        Dxu[pp] = (u[pp-2] -8.0*u[pp-1] + 8.0*u[pp+1] - u[pp+2] ) * idx_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        Dxu[IDX(3,j,k)] = ( -  3.0 * u[IDX(3,j,k)]
                            +  4.0 * u[IDX(4,j,k)]
                            -        u[IDX(5,j,k)]
                          ) * idx_by_2;
        Dxu[IDX(4,j,k)] = ( - u[IDX(3,j,k)]
                            + u[IDX(5,j,k)]
                          ) * idx_by_2;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        Dxu[IDX(ie-2,j,k)] = ( - u[IDX(ie-3,j,k)]
                               + u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

        Dxu[IDX(ie-1,j,k)] = (        u[IDX(ie-3,j,k)]
                              - 4.0 * u[IDX(ie-2,j,k)]
                              + 3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif


}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_y(double * const  Dyu, const double * const  u,
               const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0/dy;
  const double idy_by_2 = 0.50 * idy;
  const double idy_by_12 = idy / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 1;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-1;

    const int n=nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        Dyu[pp] = (u[pp-2*nx] - 8.0*u[pp-nx] + 8.0*u[pp+nx] - u[pp+2*nx])*idy_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i, 3,k)] = ( - 3.0 * u[IDX(i,3,k)]
                            +  4.0 * u[IDX(i,4,k)]
                            -        u[IDX(i,5,k)]
                          ) * idy_by_2;

        Dyu[IDX(i,4,k)] = ( - u[IDX(i,3,k)]
                            + u[IDX(i,5,k)]
                          ) * idy_by_2;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i,je-2,k)] = ( - u[IDX(i,je-3,k)]
                               + u[IDX(i,je-1,k)]
                             ) * idy_by_2;

        Dyu[IDX(i,je-1,k)] = (        u[IDX(i,je-3,k)]
                              - 4.0 * u[IDX(i,je-2,k)]
                              + 3.0 * u[IDX(i,je-1,k)]
                          ) * idy_by_2;
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_z(double * const  Dzu, const double * const  u,
               const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0/dz;
  const double idz_by_2 = 0.50 * idz;
  const double idz_by_12 = idz / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-3;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        Dzu[pp] = (u[pp-2*n] - 8.0*u[pp-n] + 8.0*u[pp+n] - u[pp+2*n]) * idz_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i, j, 3)] = ( -  3.0 *  u[IDX(i,j,3)]
                              +  4.0 * u[IDX(i,j,4)]
                              -        u[IDX(i,j,5)]
                            ) * idz_by_2;

        Dzu[IDX(i,j,4)] = ( - u[IDX(i,j,3)]
                            + u[IDX(i,j,5)]
                          ) * idz_by_2;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i,j,ke-2)] = ( - u[IDX(i,j,ke-3)]
                               + u[IDX(i,j,ke-1)]
                             ) * idz_by_2;

        Dzu[IDX(i,j,ke-1)] = (        u[IDX(i,j,ke-3)]
                              - 4.0 * u[IDX(i,j,ke-2)]
                              + 3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_2;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_xx(double * const  DxDxu, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  const double idx_sqrd = 1.0/(dx*dx);
  const double idx_sqrd_by_12 = idx_sqrd / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        DxDxu[pp] = (   -        u[pp-2]
                        + 16.0 * u[pp-1]
                        - 30.0 * u[pp  ]
                        + 16.0 * u[pp+1]
                        -        u[pp+2]
                    ) * idx_sqrd_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        DxDxu[IDX(3,j,k)] = (   2.0 * u[IDX(3,j,k)]
                              - 5.0 * u[IDX(4,j,k)]
                              + 4.0 * u[IDX(5,j,k)]
                              -       u[IDX(6,j,k)]
                            ) * idx_sqrd;

        DxDxu[IDX(4,j,k)] = (         u[IDX(3,j,k)]
                              - 2.0 * u[IDX(4,j,k)]
                              +       u[IDX(5,j,k)]
                            ) * idx_sqrd;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        DxDxu[IDX(ie-2,j,k)] = (         u[IDX(ie-3,j,k)]
                                 - 2.0 * u[IDX(ie-2,j,k)]
                                 +       u[IDX(ie-1,j,k)]
                               ) * idx_sqrd;

        DxDxu[IDX(ie-1,j,k)] = ( -       u[IDX(ie-4,j,k)]
                                 + 4.0 * u[IDX(ie-3,j,k)]
                                 - 5.0 * u[IDX(ie-2,j,k)]
                                 + 2.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_yy(double * const  DyDyu, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_12 = idy_sqrd / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        DyDyu[pp] = ( -        u[pp-2*nx] 
                      + 16.0 * u[pp-  nx] 
                      - 30.0 * u[pp     ]
                      + 16.0 * u[pp+  nx] 
                      -        u[pp+2*nx]
                    ) * idy_sqrd_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        DyDyu[IDX(i,3,k)] = (   2.0 * u[IDX(i,3,k)]
                              - 5.0 * u[IDX(i,4,k)]
                              + 4.0 * u[IDX(i,5,k)]
                              -       u[IDX(i,6,k)]
                            ) * idy_sqrd;

        DyDyu[IDX(i,4,k)] = (         u[IDX(i,3,k)]
                              - 2.0 * u[IDX(i,4,k)]
                              +       u[IDX(i,5,k)]
                            ) * idy_sqrd;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        DyDyu[IDX(i,je-2,k)] = (         u[IDX(i,je-3,k)]
                                 - 2.0 * u[IDX(i,je-2,k)]
                                 +       u[IDX(i,je-1,k)]
                               ) * idy_sqrd;

        DyDyu[IDX(i,je-1,k)] = ( -       u[IDX(i,je-4,k)]
                                 + 4.0 * u[IDX(i,je-3,k)]
                                 - 5.0 * u[IDX(i,je-2,k)]
                                 + 2.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_zz(double * const  DzDzu, const double * const  u,
                const double dz, const unsigned int *sz, unsigned bflag)
{

  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_12 = idz_sqrd / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        DzDzu[pp] = ( -        u[pp-2*n] 
                      + 16.0 * u[pp-  n] 
                      - 30.0 * u[pp    ]
                      + 16.0 * u[pp+  n] 
                      -        u[pp+2*n] 
                    ) * idz_sqrd_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        DzDzu[IDX(i,j,3)] = (   2.0 * u[IDX(i,j,3)]
                              - 5.0 * u[IDX(i,j,4)]
                              + 4.0 * u[IDX(i,j,5)]
                              -       u[IDX(i,j,6)]
                            ) * idz_sqrd;

        DzDzu[IDX(i,j,4)] = (         u[IDX(i,j,3)]
                              - 2.0 * u[IDX(i,j,4)]
                              +       u[IDX(i,j,5)]
                            ) * idz_sqrd;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        DzDzu[IDX(i,j,ke-2)] = (         u[IDX(i,j,ke-3)]
                                 - 2.0 * u[IDX(i,j,ke-2)]
                                 +       u[IDX(i,j,ke-1)]
                               ) * idz_sqrd;

        DzDzu[IDX(i,j,ke-1)] = ( -       u[IDX(i,j,ke-4)]
                                 + 4.0 * u[IDX(i,j,ke-3)]
                                 - 5.0 * u[IDX(i,j,ke-2)]
                                 + 2.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_x(double * const  Dxu, const double * const  u,
                  const double dx, const unsigned int *sz,
                  const double * const betax, unsigned bflag)
{

  const double idx = 1.0/dx;
  const double idx_by_2 = 0.50 * idx;
  const double idx_by_12 = idx / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if (betax[pp] > 0.0 ) {
          Dxu[pp] = ( -  3.0 * u[pp-1]
                      - 10.0 * u[pp]
                      + 18.0 * u[pp+1]
                      -  6.0 * u[pp+2]
                      +        u[pp+3]
                    ) * idx_by_12;
        }
        else {
          Dxu[pp] = ( -        u[pp-3]
                      +  6.0 * u[pp-2]
                      - 18.0 * u[pp-1]
                      + 10.0 * u[pp]
                      +  3.0 * u[pp+1]
                    ) * idx_by_12;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Dxu[IDX(3,j,k)] = ( -  3.0 * u[IDX(3,j,k)]
                            +  4.0 * u[IDX(4,j,k)]
                            -        u[IDX(5,j,k)]
                          ) * idx_by_2;

        if (betax[IDX(4,j,k)] > 0.0) {
          Dxu[IDX(4,j,k)] = ( -  3.0 * u[IDX(4,j,k)]
                              +  4.0 * u[IDX(5,j,k)]
                              -        u[IDX(6,j,k)]
                            ) * idx_by_2;
        }
        else {
          Dxu[IDX(4,j,k)] = ( -         u[IDX(3,j,k)]
                               +        u[IDX(5,j,k)]
                            ) * idx_by_2;
        }

        if (betax[IDX(5,j,k)] > 0.0 ) {
          Dxu[IDX(5,j,k)] = (-  3.0 * u[IDX(4,j,k)]
                             - 10.0 * u[IDX(5,j,k)]
                             + 18.0 * u[IDX(6,j,k)]
                             -  6.0 * u[IDX(7,j,k)]
                             +        u[IDX(8,j,k)]
                           ) * idx_by_12;
        }
        else {
          Dxu[IDX(5,j,k)] = (           u[IDX(3,j,k)]
                               -  4.0 * u[IDX(4,j,k)]
                               +  3.0 * u[IDX(5,j,k)]
                            ) * idx_by_2;
        }

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        if ( betax[IDX(ie-3,j,k)] > 0.0 ) {
          Dxu[IDX(ie-3,j,k)] = (  - 3.0 * u[IDX(ie-3,j,k)]
                                  + 4.0 * u[IDX(ie-2,j,k)]
                                  -       u[IDX(ie-1,j,k)]
                               ) * idx_by_2;
        }
        else {
          Dxu[IDX(ie-3,j,k)] = ( -   u[IDX(ie-6,j,k)]
                            +  6.0 * u[IDX(ie-5,j,k)]
                            - 18.0 * u[IDX(ie-4,j,k)]
                            + 10.0 * u[IDX(ie-3  ,j,k)]
                            +  3.0 * u[IDX(ie-2,j,k)]
                          ) * idx_by_12;
        }

        if (betax[IDX(ie-2,j,k)] > 0.0 ) {
          Dxu[IDX(ie-2,j,k)] = (  -  u[IDX(ie-3,j,k)]
                                  +  u[IDX(ie-1,j,k)]
                               ) * idx_by_2;
        }
        else {
          Dxu[IDX(ie-2,j,k)] = (     u[IDX(ie-4,j,k)]
                             - 4.0 * u[IDX(ie-3,j,k)]
                             + 3.0 * u[IDX(ie-2,j,k)]
                               ) * idx_by_2;
        }

        Dxu[IDX(ie-1,j,k)] = (          u[IDX(ie-3,j,k)]
                                - 4.0 * u[IDX(ie-2,j,k)]
                                + 3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_y(double * const  Dyu, const double * const  u,
                  const double dy, const unsigned int *sz,
                  const double * const betay, unsigned bflag)
{

  const double idy = 1.0/dy;
  const double idy_by_2 = 0.50 * idy;
  const double idy_by_12 = idy / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        if (betay[pp] > 0.0 ) {
          Dyu[pp] = ( -  3.0 * u[pp-nx]
                      - 10.0 * u[pp]
                      + 18.0 * u[pp+nx]
                      -  6.0 * u[pp+2*nx]
                      +        u[pp+3*nx]
                    ) * idy_by_12;
        }
        else {
          Dyu[pp] = ( - u[pp-3*nx]
                      +  6.0 * u[pp-2*nx]
                      - 18.0 * u[pp-nx]
                      + 10.0 * u[pp]
                      +  3.0 * u[pp+nx]
                    ) * idy_by_12;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i,3,k)] = ( -  3.0 * u[IDX(i,3,k)]
                            +  4.0 * u[IDX(i,4,k)]
                            -        u[IDX(i,5,k)]
                          ) * idy_by_2;

        if (betay[IDX(i,4,k)] > 0.0) {
          Dyu[IDX(i,4,k)] = ( -  3.0 * u[IDX(i,4,k)]
                              +  4.0 * u[IDX(i,5,k)]
                              -        u[IDX(i,6,k)]
                            ) * idy_by_2;
        }
        else {
          Dyu[IDX(i,4,k)] = ( -         u[IDX(i,3,k)]
                               +        u[IDX(i,5,k)]
                            ) * idy_by_2;
        }

        if (betay[IDX(i,5,k)] > 0.0 ) {
          Dyu[IDX(i,5,k)] = ( -  3.0 * u[IDX(i,4,k)]
                              - 10.0 * u[IDX(i,5,k)]
                              + 18.0 * u[IDX(i,6,k)]
                              -  6.0 * u[IDX(i,7,k)]
                             +         u[IDX(i,8,k)]
                           ) * idy_by_12;
        }
        else {
          Dyu[IDX(i,5,k)] = (           u[IDX(i,3,k)]
                               -  4.0 * u[IDX(i,4,k)]
                               +  3.0 * u[IDX(i,5,k)]
                            ) * idy_by_2;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      for (int i = ib; i < ie; i++) {
        if ( betay[IDX(i,je-3,k)] > 0.0 ) {
          Dyu[IDX(i,je-3,k)] = (  - 3.0 * u[IDX(i,je-3,k)]
                                + 4.0 * u[IDX(i,je-2,k)]
                                -       u[IDX(i,je-1,k)]
                             ) * idy_by_2;
        }
        else {
          Dyu[IDX(i,je-3,k)] = ( -   u[IDX(i,je-6,k)]
                            +  6.0 * u[IDX(i,je-5,k)]
                            - 18.0 * u[IDX(i,je-4,k)]
                            + 10.0 * u[IDX(i,je-3,k)]
                            +  3.0 * u[IDX(i,je-2,k)]
                          ) * idy_by_12;
        }

        if (betay[IDX(i,je-2,k)] > 0.0 ) {
          Dyu[IDX(i,je-2,k)] = (  -  u[IDX(i,je-3,k)]
                                  +  u[IDX(i,je-1,k)]
                               ) * idy_by_2;
        }
        else {
          Dyu[IDX(i,je-2,k)] = (     u[IDX(i,je-4,k)]
                             - 4.0 * u[IDX(i,je-3,k)]
                             + 3.0 * u[IDX(i,je-2,k)]
                               ) * idy_by_2;
        }

        Dyu[IDX(i,je-1,k)] = (          u[IDX(i,je-3,k)]
                                - 4.0 * u[IDX(i,je-2,k)]
                                + 3.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_2;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42adv_z(double * const  Dzu, const double * const  u,
                  const double dz, const unsigned int *sz,
                  const double * const betaz, unsigned bflag)
{

  const double idz = 1.0/dz;
  const double idz_by_2 = 0.50 * idz;
  const double idz_by_12 = idz / 12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        if (betaz[pp] > 0.0 ) {
          Dzu[pp] = ( -  3.0 * u[pp-n]
                      - 10.0 * u[pp]
                      + 18.0 * u[pp+n]
                      -  6.0 * u[pp+2*n]
                      +        u[pp+3*n]
                    ) * idz_by_12;
        }
        else {
          Dzu[pp] = ( -        u[pp-3*n]
                      +  6.0 * u[pp-2*n]
                      - 18.0 * u[pp-n]
                      + 10.0 * u[pp]
                      +  3.0 * u[pp+n]
                    ) * idz_by_12;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i,j,3)] = ( -  3.0 * u[IDX(i,j,3)]
                            +  4.0 * u[IDX(i,j,4)]
                            -        u[IDX(i,j,5)]
                          ) * idz_by_2;

        if (betaz[IDX(i,j,4)] > 0.0) {
          Dzu[IDX(i,j,4)] = ( -  3.0 * u[IDX(i,j,4)]
                              +  4.0 * u[IDX(i,j,5)]
                              -        u[IDX(i,j,6)]
                            ) * idz_by_2;
        }
        else {
          Dzu[IDX(i,j,4)] = ( -         u[IDX(i,j,3)]
                               +        u[IDX(i,j,5)]
                            ) * idz_by_2;
        }

        if (betaz[IDX(i,j,5)] > 0.0 ) {
          Dzu[IDX(i,j,5)] = ( -  3.0 * u[IDX(i,j,4)]
                              - 10.0 * u[IDX(i,j,5)]
                              + 18.0 * u[IDX(i,j,6)]
                              -  6.0 * u[IDX(i,j,7)]
                             +         u[IDX(i,j,8)]
                           ) * idz_by_12;
        }
        else {
          Dzu[IDX(i,j,5)] = (           u[IDX(i,j,3)]
                               -  4.0 * u[IDX(i,j,4)]
                               +  3.0 * u[IDX(i,j,5)]
                            ) * idz_by_2;
        }
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        if ( betaz[IDX(i,j,ke-3)] > 0.0 ) {
          Dzu[IDX(i,j,ke-3)] = (  - 3.0 * u[IDX(i,j,ke-3)]
                                  + 4.0 * u[IDX(i,j,ke-2)]
                                  -       u[IDX(i,j,ke-1)]
                               ) * idz_by_2;
        }
        else {
          Dzu[IDX(i,j,ke-3)] = ( -        u[IDX(i,j,ke-6)]
                                 +  6.0 * u[IDX(i,j,ke-5)]
                                 - 18.0 * u[IDX(i,j,ke-4)]
                                 + 10.0 * u[IDX(i,j,ke-3)]
                                 +  3.0 * u[IDX(i,j,ke-2)]
                               ) * idz_by_12;
        }

        if (betaz[IDX(i,j,ke-2)] > 0.0 ) {
          Dzu[IDX(i,j,ke-2)] = (  -  u[IDX(i,j,ke-3)]
                                  +  u[IDX(i,j,ke-1)]
                               ) * idz_by_2;
        }
        else {
          Dzu[IDX(i,j,ke-2)] = (         u[IDX(i,j,ke-4)]
                                 - 4.0 * u[IDX(i,j,ke-3)]
                                 + 3.0 * u[IDX(i,j,ke-2)]
                               ) * idz_by_2;
        }

        Dzu[IDX(i,j,ke-1)] = (          u[IDX(i,j,ke-3)]
                                - 4.0 * u[IDX(i,j,ke-2)]
                                + 3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_2;
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv42_xy( double * const DxDyu, double* const Dxu, 
                 const double * const u, const double dx, 
                 const double dy, const unsigned int *sz, unsigned bflag )
{
    deriv42_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv42_y( DxDyu, Dxu, dy, sz, bflag ) ; 
} 


void deriv42_xz( double * const DxDzu, double* const Dxu, const double * const u, const double dx, 
                 const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv42_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv42_z( DxDzu, Dxu, dz, sz, bflag ) ; 
} 


void deriv42_yz( double * const DyDzu, double* const Dyu, const double * const u, const double dy, 
                 const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv42_y( Dyu  ,   u, dy, sz, bflag ) ; 
    deriv42_z( DyDzu, Dyu, dz, sz, bflag ) ; 
} 

/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv42_x(double * const  Du, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_6_dx = -1.0 / 64.0 / dx;

  double smr3=59.0/48.0*64*dx;
  double smr2=43.0/48.0*64*dx;
  double smr1=49.0/48.0*64*dx;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];

  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;
  
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
       for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_6_dx *
                         (
                         -      u[pp-3]
                         +  6.0*u[pp-2]
                         - 15.0*u[pp-1]
                         + 20.0*u[pp]
                         - 15.0*u[pp+1]
                         +  6.0*u[pp+2]
                         -      u[pp+3]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        Du[IDX(3,j,k)] =  (      u[IDX(6,j,k)]
                           - 3.0*u[IDX(5,j,k)]
                           + 3.0*u[IDX(4,j,k)]
                           -     u[IDX(3,j,k)]
                          )/smr3;
        Du[IDX(4,j,k)] =  (
                                 u[IDX(7,j,k)]
                          -  6.0*u[IDX(6,j,k)]
                          + 12.0*u[IDX(5,j,k)]
                          - 10.0*u[IDX(4,j,k)]
                          +  3.0*u[IDX(3,j,k)]
                          )/smr2;
        Du[IDX(5,j,k)] =  (
                                 u[IDX(8,j,k)]
                          -  6.0*u[IDX(7,j,k)]
                          + 15.0*u[IDX(6,j,k)]
                          - 19.0*u[IDX(5,j,k)]
                          + 12.0*u[IDX(4,j,k)]
                          -  3.0*u[IDX(3,j,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
         Du[IDX(ie-3,j,k)] = (
                                 u[IDX(ie-6,j,k)]
                          -  6.0*u[IDX(ie-5,j,k)]
                          + 15.0*u[IDX(ie-4,j,k)]
                          - 19.0*u[IDX(ie-3,j,k)]
                          + 12.0*u[IDX(ie-2,j,k)]
                          -  3.0*u[IDX(ie-1,j,k)]
                           )/spr1;

         Du[IDX(ie-2,j,k)] = (
                                 u[IDX(ie-5,j,k)]
                          -  6.0*u[IDX(ie-4,j,k)]
                          + 12.0*u[IDX(ie-3,j,k)]
                          - 10.0*u[IDX(ie-2,j,k)]
                          +  3.0*u[IDX(ie-1,j,k)]
                           )/spr2;

         Du[IDX(ie-1,j,k)] = (
                                 u[IDX(ie-4,j,k)]
                          -  3.0*u[IDX(ie-3,j,k)]
                          +  3.0*u[IDX(ie-2,j,k)]
                          -      u[IDX(ie-1,j,k)]
                           )/spr3;
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv42_y(double * const  Du, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_6_dy = -1.0 / 64.0 / dy;

  double smr3=59.0/48.0*64*dy;
  double smr2=43.0/48.0*64*dy;
  double smr1=49.0/48.0*64*dy;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;
  
  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
       for (int j = jb; j < je; j++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_6_dy *
                         (
                         -      u[pp-3*nx]
                         +  6.0*u[pp-2*nx]
                         - 15.0*u[pp-nx]
                         + 20.0*u[pp]
                         - 15.0*u[pp+nx]
                         +  6.0*u[pp+2*nx]
                         -      u[pp+3*nx]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,3,k)] =  (      u[IDX(i,6,k)]
                           - 3.0*u[IDX(i,5,k)]
                           + 3.0*u[IDX(i,4,k)]
                           -     u[IDX(i,3,k)]
                          )/smr3;
        Du[IDX(i,4,k)] =  (
                                 u[IDX(i,7,k)]
                          -  6.0*u[IDX(i,6,k)]
                          + 12.0*u[IDX(i,5,k)]
                          - 10.0*u[IDX(i,4,k)]
                          +  3.0*u[IDX(i,3,k)]
                          )/smr2;
        Du[IDX(i,5,k)] =  (
                                 u[IDX(i,8,k)]
                          -  6.0*u[IDX(i,7,k)]
                          + 15.0*u[IDX(i,6,k)]
                          - 19.0*u[IDX(i,5,k)]
                          + 12.0*u[IDX(i,4,k)]
                          -  3.0*u[IDX(i,3,k)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,je-3,k)] = (
                                 u[IDX(i,je-6,k)]
                          -  6.0*u[IDX(i,je-5,k)]
                          + 15.0*u[IDX(i,je-4,k)]
                          - 19.0*u[IDX(i,je-3,k)]
                          + 12.0*u[IDX(i,je-2,k)]
                          -  3.0*u[IDX(i,je-1,k)]
                           )/spr1;

       Du[IDX(i,je-2,k)] = (
                                 u[IDX(i,je-5,k)]
                          -  6.0*u[IDX(i,je-4,k)]
                          + 12.0*u[IDX(i,je-3,k)]
                          - 10.0*u[IDX(i,je-2,k)]
                          +  3.0*u[IDX(i,je-1,k)]
                           )/spr2;

       Du[IDX(i,je-1,k)] = (
                                 u[IDX(i,je-4,k)]
                          -  3.0*u[IDX(i,je-3,k)]
                          +  3.0*u[IDX(i,je-2,k)]
                          -      u[IDX(i,je-1,k)]
                           )/spr3;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv42_z(double * const  Du, const double * const  u,
                const double dz, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_6_dz = -1.0 / 64.0 / dz;

  double smr3=59.0/48.0*64*dz;
  double smr2=43.0/48.0*64*dz;
  double smr1=49.0/48.0*64*dz;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];

  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;
  
  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
       for (int k = kb; k < ke; k++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_6_dz *
                         (
                         -      u[pp-3*n]
                         +  6.0*u[pp-2*n]
                         - 15.0*u[pp-n]
                         + 20.0*u[pp]
                         - 15.0*u[pp+n]
                         +  6.0*u[pp+2*n]
                         -      u[pp+3*n]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,j,3)] =  (      u[IDX(i,j,6)]
                           - 3.0*u[IDX(i,j,5)]
                           + 3.0*u[IDX(i,j,4)]
                           -     u[IDX(i,j,3)]
                          )/smr3;
        Du[IDX(i,j,4)] =  (
                                 u[IDX(i,j,7)]
                          -  6.0*u[IDX(i,j,6)]
                          + 12.0*u[IDX(i,j,5)]
                          - 10.0*u[IDX(i,j,4)]
                          +  3.0*u[IDX(i,j,3)]
                          )/smr2;
        Du[IDX(i,j,5)] =  (
                                 u[IDX(i,j,8)]
                          -  6.0*u[IDX(i,j,7)]
                          + 15.0*u[IDX(i,j,6)]
                          - 19.0*u[IDX(i,j,5)]
                          + 12.0*u[IDX(i,j,4)]
                          -  3.0*u[IDX(i,j,3)]
                          )/smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,j,ke-3)] = (
                                 u[IDX(i,j,ke-6)]
                          -  6.0*u[IDX(i,j,ke-5)]
                          + 15.0*u[IDX(i,j,ke-4)]
                          - 19.0*u[IDX(i,j,ke-3)]
                          + 12.0*u[IDX(i,j,ke-2)]
                          -  3.0*u[IDX(i,j,ke-1)]
                           )/spr1;

         Du[IDX(i,j,ke-2)] = (
                                 u[IDX(i,j,ke-5)]
                          -  6.0*u[IDX(i,j,ke-4)]
                          + 12.0*u[IDX(i,j,ke-3)]
                          - 10.0*u[IDX(i,j,ke-2)]
                          +  3.0*u[IDX(i,j,ke-1)]
                           )/spr2;

         Du[IDX(i,j,ke-1)] = (
                                 u[IDX(i,j,ke-4)]
                          -  3.0*u[IDX(i,j,ke-3)]
                          +  3.0*u[IDX(i,j,ke-2)]
                          -      u[IDX(i,j,ke-1)]
                           )/spr3;
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv64_x(double * const  Du, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_8_dx = - 1.0 / 256.0 / dx;

  double smr4=17.0/48.0*256*dx;
  double smr3=59.0/48.0*256*dx;
  double smr2=43.0/48.0*256*dx;
  double smr1=49.0/48.0*256*dx;
  double spr4=smr4;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;



  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_8_dx *
                         (          u[pp-4]
                           -  8.0 * u[pp-3]
                           + 28.0 * u[pp-2]
                           - 56.0 * u[pp-1]
                           + 70.0 * u[pp  ]
                           - 56.0 * u[pp+1]
                           + 28.0 * u[pp+2]
                           -  8.0 * u[pp+3]
                           +        u[pp+4]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        Du[IDX(4,j,k)] =  ( -       u[IDX(4,j,k)]
                            + 4.0 * u[IDX(5,j,k)]
                            - 6.0 * u[IDX(6,j,k)]
                            + 4.0 * u[IDX(7,j,k)]
                            -       u[IDX(8,j,k)]
                          ) / smr4;
        
        Du[IDX(5,j,k)] =  (    3.0 * u[IDX(4,j,k)]
                            - 11.0 * u[IDX(5,j,k)]
                            + 15.0 * u[IDX(6,j,k)]
                            -  9.0 * u[IDX(7,j,k)]
                            +  2.0 * u[IDX(8,j,k)]
                          ) / smr3;
        
        Du[IDX(6,j,k)] =  ( - 3.0 * u[IDX(4,j,k)]
                            + 9.0 * u[IDX(5,j,k)]
                            - 8.0 * u[IDX(6,j,k)]
                            + 3.0 * u[IDX(8,j,k)]
                            -       u[IDX(9,j,k)]
                          ) / smr2;
        
        Du[IDX(7,j,k)] =  (          u[IDX( 4,j,k)]
                            -        u[IDX( 5,j,k)]
                            -  6.0 * u[IDX( 6,j,k)]
                            + 15.0 * u[IDX( 7,j,k)]
                            - 14.0 * u[IDX( 8,j,k)]
                            +  6.0 * u[IDX( 9,j,k)]
                            -        u[IDX(10,j,k)]
                          ) / smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
         Du[IDX(ie-4,j,k)] = ( -        u[IDX(ie-7,j,k)]
                               +  6.0 * u[IDX(ie-6,j,k)]
                               - 14.0 * u[IDX(ie-5,j,k)]
                               + 15.0 * u[IDX(ie-4,j,k)]
                               -  6.0 * u[IDX(ie-3,j,k)]
                               -        u[IDX(ie-2,j,k)]
                               +        u[IDX(ie-1,j,k)]
                             ) / spr1;

         Du[IDX(ie-3,j,k)] = ( -        u[IDX(ie-6,j,k)]
                               +  3.0 * u[IDX(ie-5,j,k)]
                               -  8.0 * u[IDX(ie-3,j,k)]
                               +  9.0 * u[IDX(ie-2,j,k)]
                               -  3.0 * u[IDX(ie-1,j,k)]
                             ) / spr2;

         Du[IDX(ie-2,j,k)] = (    2.0 * u[IDX(ie-5,j,k)]
                               -  9.0 * u[IDX(ie-4,j,k)]
                               + 15.0 * u[IDX(ie-3,j,k)]
                               - 11.0 * u[IDX(ie-2,j,k)]
                               +  3.0 * u[IDX(ie-1,j,k)]
                             ) / spr3;

         Du[IDX(ie-1,j,k)] = ( -       u[IDX(ie-5,j,k)]
                               + 4.0 * u[IDX(ie-4,j,k)]
                               - 6.0 * u[IDX(ie-3,j,k)]
                               + 4.0 * u[IDX(ie-2,j,k)]
                               -       u[IDX(ie-1,j,k)]
                             ) / spr4;
      
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}


/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv64_y(double * const  Du, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_8_dy = - 1.0 / 256.0 / dy;

  double smr4=17.0/48.0*256*dy;
  double smr3=59.0/48.0*256*dy;
  double smr2=43.0/48.0*256*dy;
  double smr1=49.0/48.0*256*dy;
  double spr4=smr4;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;
  
  const int n = nx ; 

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_8_dy *
                         (          u[pp-4*n]
                           -  8.0 * u[pp-3*n]
                           + 28.0 * u[pp-2*n]
                           - 56.0 * u[pp-  n]
                           + 70.0 * u[pp  ]
                           - 56.0 * u[pp+  n]
                           + 28.0 * u[pp+2*n]
                           -  8.0 * u[pp+3*n]
                           +        u[pp+4*n]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,4,k)] =  ( -       u[IDX(i,4,k)]
                            + 4.0 * u[IDX(i,5,k)]
                            - 6.0 * u[IDX(i,6,k)]
                            + 4.0 * u[IDX(i,7,k)]
                            -       u[IDX(i,8,k)]
                          ) / smr4;
        
        Du[IDX(i,5,k)] =  (    3.0 * u[IDX(i,4,k)]
                            - 11.0 * u[IDX(i,5,k)]
                            + 15.0 * u[IDX(i,6,k)]
                            -  9.0 * u[IDX(i,7,k)]
                            +  2.0 * u[IDX(i,8,k)]
                          ) / smr3;
        
        Du[IDX(i,6,k)] =  ( - 3.0 * u[IDX(i,4,k)]
                            + 9.0 * u[IDX(i,5,k)]
                            - 8.0 * u[IDX(i,6,k)]
                            + 3.0 * u[IDX(i,8,k)]
                            -       u[IDX(i,9,k)]
                          ) / smr2;
        
        Du[IDX(i,7,k)] =  (          u[IDX(i, 4,k)]
                            -        u[IDX(i, 5,k)]
                            -  6.0 * u[IDX(i, 6,k)]
                            + 15.0 * u[IDX(i, 7,k)]
                            - 14.0 * u[IDX(i, 8,k)]
                            +  6.0 * u[IDX(i, 9,k)]
                            -        u[IDX(i,10,k)]
                          ) / smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,je-4,k)] = ( -        u[IDX(i,je-7,k)]
                               +  6.0 * u[IDX(i,je-6,k)]
                               - 14.0 * u[IDX(i,je-5,k)]
                               + 15.0 * u[IDX(i,je-4,k)]
                               -  6.0 * u[IDX(i,je-3,k)]
                               -        u[IDX(i,je-2,k)]
                               +        u[IDX(i,je-1,k)]
                             ) / spr1;

         Du[IDX(i,je-3,k)] = ( -        u[IDX(i,je-6,k)]
                               +  3.0 * u[IDX(i,je-5,k)]
                               -  8.0 * u[IDX(i,je-3,k)]
                               +  9.0 * u[IDX(i,je-2,k)]
                               -  3.0 * u[IDX(i,je-1,k)]
                             ) / spr2;

         Du[IDX(i,je-2,k)] = (    2.0 * u[IDX(i,je-5,k)]
                               -  9.0 * u[IDX(i,je-4,k)]
                               + 15.0 * u[IDX(i,je-3,k)]
                               - 11.0 * u[IDX(i,je-2,k)]
                               +  3.0 * u[IDX(i,je-1,k)]
                             ) / spr3;

         Du[IDX(i,je-1,k)] = ( -       u[IDX(i,je-5,k)]
                               + 4.0 * u[IDX(i,je-4,k)]
                               - 6.0 * u[IDX(i,je-3,k)]
                               + 4.0 * u[IDX(i,je-2,k)]
                               -       u[IDX(i,je-1,k)]
                             ) / spr4;
      
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}


/*----------------------------------------------------------------------
 *
 *
 *
 *----------------------------------------------------------------------*/
void ko_deriv64_z(double * const  Du, const double * const  u,
                const double dz, const unsigned int *sz, unsigned bflag)
{

  double pre_factor_8_dz = - 1.0 / 256.0 / dz;

  double smr4=17.0/48.0*256*dz;
  double smr3=59.0/48.0*256*dz;
  double smr2=43.0/48.0*256*dz;
  double smr1=49.0/48.0*256*dz;
  double spr4=smr4;
  double spr3=smr3;
  double spr2=smr2;
  double spr1=smr1;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;
  
  const int n = nx * ny ; 

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
          const int pp = IDX(i,j,k);
          Du[pp] = pre_factor_8_dz *
                         (          u[pp-4*n]
                           -  8.0 * u[pp-3*n]
                           + 28.0 * u[pp-2*n]
                           - 56.0 * u[pp-  n]
                           + 70.0 * u[pp  ]
                           - 56.0 * u[pp+  n]
                           + 28.0 * u[pp+2*n]
                           -  8.0 * u[pp+3*n]
                           +        u[pp+4*n]
                         );
       }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Du[IDX(i,j,4)] =  ( -       u[IDX(i,j,4)]
                            + 4.0 * u[IDX(i,j,5)]
                            - 6.0 * u[IDX(i,j,6)]
                            + 4.0 * u[IDX(i,j,7)]
                            -       u[IDX(i,j,8)]
                          ) / smr4;
        
        Du[IDX(i,j,5)] =  (    3.0 * u[IDX(i,j,4)]
                            - 11.0 * u[IDX(i,j,5)]
                            + 15.0 * u[IDX(i,j,6)]
                            -  9.0 * u[IDX(i,j,7)]
                            +  2.0 * u[IDX(i,j,8)]
                          ) / smr3;
        
        Du[IDX(i,j,6)] =  ( - 3.0 * u[IDX(i,j,4)]
                            + 9.0 * u[IDX(i,j,5)]
                            - 8.0 * u[IDX(i,j,6)]
                            + 3.0 * u[IDX(i,j,8)]
                            -       u[IDX(i,j,9)]
                          ) / smr2;
        
        Du[IDX(i,j,7)] =  (          u[IDX(i,j, 4)]
                            -        u[IDX(i,j, 5)]
                            -  6.0 * u[IDX(i,j, 6)]
                            + 15.0 * u[IDX(i,j, 7)]
                            - 14.0 * u[IDX(i,j, 8)]
                            +  6.0 * u[IDX(i,j, 9)]
                            -        u[IDX(i,j,10)]
                          ) / smr1;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
         Du[IDX(i,j,ke-4)] = ( -        u[IDX(i,j,ke-7)]
                               +  6.0 * u[IDX(i,j,ke-6)]
                               - 14.0 * u[IDX(i,j,ke-5)]
                               + 15.0 * u[IDX(i,j,ke-4)]
                               -  6.0 * u[IDX(i,j,ke-3)]
                               -        u[IDX(i,j,ke-2)]
                               +        u[IDX(i,j,ke-1)]
                             ) / spr1;

         Du[IDX(i,j,ke-3)] = ( -        u[IDX(i,j,ke-6)]
                               +  3.0 * u[IDX(i,j,ke-5)]
                               -  8.0 * u[IDX(i,j,ke-3)]
                               +  9.0 * u[IDX(i,j,ke-2)]
                               -  3.0 * u[IDX(i,j,ke-1)]
                             ) / spr2;

         Du[IDX(i,j,ke-2)] = (    2.0 * u[IDX(i,j,ke-5)]
                               -  9.0 * u[IDX(i,j,ke-4)]
                               + 15.0 * u[IDX(i,j,ke-3)]
                               - 11.0 * u[IDX(i,j,ke-2)]
                               +  3.0 * u[IDX(i,j,ke-1)]
                             ) / spr3;

         Du[IDX(i,j,ke-1)] = ( -       u[IDX(i,j,ke-5)]
                               + 4.0 * u[IDX(i,j,ke-4)]
                               - 6.0 * u[IDX(i,j,ke-3)]
                               + 4.0 * u[IDX(i,j,ke-2)]
                               -       u[IDX(i,j,ke-1)]
                             ) / spr4;
      
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        if(std::isnan(Du[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
      }
    }
  }
#endif

}





/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv644_x(double * const  Dxu, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx = 1.0 / dx;
  const double idx_by_12 = idx / 12.0;
  const double idx_by_60 = idx / 60.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 0;
  const int kb = 0;
  const int ie = sz[0]-3;
  const int je = sz[1]-0;
  const int ke = sz[2]-0;
  
  const int n = 1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        Dxu[pp] = ( -        u[pp-3] 
                    +  9.0 * u[pp-2]
                    - 45.0 * u[pp-1]
                    + 45.0 * u[pp+1]
                    -  9.0 * u[pp+2]
                    +        u[pp+3] ) * idx_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        Dxu[IDX(3,j,k)] = ( - 25.0 * u[IDX(3,j,k)]
                            + 48.0 * u[IDX(4,j,k)]
                            - 36.0 * u[IDX(5,j,k)]
                            + 16.0 * u[IDX(6,j,k)]
                            -  3.0 * u[IDX(7,j,k)]
                          ) * idx_by_12;
        
        Dxu[IDX(4,j,k)] = ( -  3.0 * u[IDX(3,j,k)]
                            - 10.0 * u[IDX(4,j,k)]
                            + 18.0 * u[IDX(5,j,k)]
                            -  6.0 * u[IDX(6,j,k)]
                            +        u[IDX(7,j,k)]
                          ) * idx_by_12;

        Dxu[IDX(5,j,k)] = ( +       u[IDX(3,j,k)]
                            - 8.0 * u[IDX(4,j,k)]
                            + 8.0 * u[IDX(6,j,k)]
                            -       u[IDX(7,j,k)]
                          ) * idx_by_12;
        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {

        Dxu[IDX(ie-3,j,k)] = ( +       u[IDX(ie-5,j,k)]
                               - 8.0 * u[IDX(ie-4,j,k)]
                               + 8.0 * u[IDX(ie-2,j,k)]
                               -       u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

        Dxu[IDX(ie-2,j,k)] = ( -        u[IDX(ie-5,j,k)]
                               +  6.0 * u[IDX(ie-4,j,k)]
                               - 18.0 * u[IDX(ie-3,j,k)]
                               + 10.0 * u[IDX(ie-2,j,k)]
                               +  3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;
                             
        Dxu[IDX(ie-1,j,k)] = (    3.0 * u[IDX(ie-5,j,k)]
                               - 16.0 * u[IDX(ie-4,j,k)]
                               + 36.0 * u[IDX(ie-3,j,k)]
                               - 48.0 * u[IDX(ie-2,j,k)]
                               + 25.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif


}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv644_y(double * const  Dyu, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0 / dy;
  const double idy_by_12 = idy / 12.0;
  const double idy_by_60 = idy / 60.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 0;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-0;
  
  const int n = nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        Dyu[pp] = ( -        u[pp-3*nx] 
                    +  9.0 * u[pp-2*nx]
                    - 45.0 * u[pp-  nx]
                    + 45.0 * u[pp+  nx]
                    -  9.0 * u[pp+2*nx]
                    +        u[pp+3*nx]) * idy_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Dyu[IDX(i,3,k)] =  ( - 25.0 * u[IDX(i,3,k)]
                             + 48.0 * u[IDX(i,4,k)]
                             - 36.0 * u[IDX(i,5,k)]
                             + 16.0 * u[IDX(i,6,k)]
                             -  3.0 * u[IDX(i,7,k)]
                           ) * idy_by_12;

        Dyu[IDX(i,4,k)] = ( -  3.0 * u[IDX(i,3,k)]
                            - 10.0 * u[IDX(i,4,k)]
                            + 18.0 * u[IDX(i,5,k)]
                            -  6.0 * u[IDX(i,6,k)]
                            +        u[IDX(i,7,k)]
                          ) * idy_by_12;

        Dyu[IDX(i,5,k)] = (         u[IDX(i,3,k)]
                            - 8.0 * u[IDX(i,4,k)]
                            + 8.0 * u[IDX(i,6,k)]
                            -       u[IDX(i,7,k)]
                          ) * idy_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        

        Dyu[IDX(i,je-3,k)] = (         u[IDX(i,je-5,k)]
                               - 8.0 * u[IDX(i,je-4,k)]
                               + 8.0 * u[IDX(i,je-2,k)]
                               -       u[IDX(i,je-1,k)]
                             ) * idy_by_12;


        Dyu[IDX(i,je-2,k)] = ( -        u[IDX(i,je-5,k)]
                               +  6.0 * u[IDX(i,je-4,k)]
                               - 18.0 * u[IDX(i,je-3,k)]
                               + 10.0 * u[IDX(i,je-2,k)]
                               +  3.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;
                             

        Dyu[IDX(i,je-1,k)] = ( +  3.0 * u[IDX(i,je-5,k)]
                               - 16.0 * u[IDX(i,je-4,k)]
                               + 36.0 * u[IDX(i,je-3,k)]
                               - 48.0 * u[IDX(i,je-2,k)]
                               + 25.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
//void deriv64_z(double * const  Dzu, const double * const  u,
void deriv644_z(double * const  Dzu, const double * const  u,
                const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0 / dz;
  const double idz_by_12 = idz / 12.0;
  const double idz_by_60 = idz / 60.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-3;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        Dzu[pp] = ( -        u[pp-3*n] 
                    +  9.0 * u[pp-2*n]
                    - 45.0 * u[pp-  n]
                    + 45.0 * u[pp+  n]
                    -  9.0 * u[pp+2*n]
                    +        u[pp+3*n]) * idz_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i, j,3)] = ( - 25.0 * u[IDX(i,j,3)]
                             + 48.0 * u[IDX(i,j,4)]
                             - 36.0 * u[IDX(i,j,5)]
                             + 16.0 * u[IDX(i,j,6)]
                             -  3.0 * u[IDX(i,j,7)]
                          ) * idz_by_12;

        Dzu[IDX(i,j,4)] = ( -  3.0 * u[IDX(i,j,3)]
                            - 10.0 * u[IDX(i,j,4)]
                            + 18.0 * u[IDX(i,j,5)]
                            -  6.0 * u[IDX(i,j,6)]
                            +        u[IDX(i,j,7)]
                          ) * idz_by_12;

        Dzu[IDX(i,j,5)] = (         u[IDX(i,j,3)]
                            - 8.0 * u[IDX(i,j,4)]
                            + 8.0 * u[IDX(i,j,6)]
                            -       u[IDX(i,j,7)]
                          ) * idz_by_12;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        Dzu[IDX(i,j,ke-3)] = (         u[IDX(i,j,ke-5)]
                               - 8.0 * u[IDX(i,j,ke-4)]
                               + 8.0 * u[IDX(i,j,ke-2)]
                               -       u[IDX(i,j,ke-1)]
                             ) * idz_by_12;

        Dzu[IDX(i,j,ke-2)] = ( -        u[IDX(i,j,ke-5)]
                               +  6.0 * u[IDX(i,j,ke-4)]
                               - 18.0 * u[IDX(i,j,ke-3)]
                               + 10.0 * u[IDX(i,j,ke-2)]
                               +  3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_12;
                             
        Dzu[IDX(i,j,ke-1)] = (    3.0 * u[IDX(i,j,ke-5)]
                               - 16.0 * u[IDX(i,j,ke-4)]
                               + 36.0 * u[IDX(i,j,ke-3)]
                               - 48.0 * u[IDX(i,j,ke-2)]
                               + 25.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_12;

      }
    }
  }

  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv644_xx(double * const  DxDxu, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{

  const double idx_sqrd = 1.0 / (dx*dx);
  const double idx_sqrd_by_180 = idx_sqrd / 180.0;
  const double idx_sqrd_by_12  = idx_sqrd / 12.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);

        DxDxu[pp] = (     2.0 * u[pp-3]
                      -  27.0 * u[pp-2]
                      + 270.0 * u[pp-1]
                      - 490.0 * u[pp  ]
                      + 270.0 * u[pp+1]
                      -  27.0 * u[pp+2]
                      +   2.0 * u[pp+3]
                    ) * idx_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        // The above two should be replaced by 4th order approximations: 
        DxDxu[IDX(3,j,k)] = (    45.0 * u[IDX(3,j,k)]
                              - 154.0 * u[IDX(4,j,k)]
                              + 214.0 * u[IDX(5,j,k)]
                              - 156.0 * u[IDX(6,j,k)]
                              +  61.0 * u[IDX(7,j,k)]
                              -  10.0 * u[IDX(8,j,k)]
                            ) * idx_sqrd_by_12;
          
        DxDxu[IDX(4,j,k)] = (    10.0 * u[IDX(3,j,k)]
                              -  15.0 * u[IDX(4,j,k)]
                              -   4.0 * u[IDX(5,j,k)]
                              +  14.0 * u[IDX(6,j,k)]
                              -   6.0 * u[IDX(7,j,k)]
                              +         u[IDX(8,j,k)]
                            ) * idx_sqrd_by_12;
        
        DxDxu[IDX(5,j,k)] = ( -         u[IDX(3,j,k)]
                              +  16.0 * u[IDX(4,j,k)]
                              -  30.0 * u[IDX(5,j,k)]
                              +  16.0 * u[IDX(6,j,k)]
                              -         u[IDX(7,j,k)]
                            ) * idx_sqrd_by_12;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        DxDxu[IDX(ie-3,j,k)] = ( -         u[IDX(ie-5,j,k)]
                                 +  16.0 * u[IDX(ie-4,j,k)]
                                 -  30.0 * u[IDX(ie-3,j,k)]
                                 +  16.0 * u[IDX(ie-2,j,k)]
                                 -         u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;
        
        // The above two should be replaced by 4th order approximations: 
        DxDxu[IDX(ie-2,j,k)] = (           u[IDX(ie-6,j,k)]
                                 -   6.0 * u[IDX(ie-5,j,k)]
                                 +  14.0 * u[IDX(ie-4,j,k)]
                                 -   4.0 * u[IDX(ie-3,j,k)]
                                 -  15.0 * u[IDX(ie-2,j,k)]
                                 +  10.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;
         
        DxDxu[IDX(ie-1,j,k)] = ( -  10.0 * u[IDX(ie-6,j,k)]
                                 +  61.0 * u[IDX(ie-5,j,k)]
                                 - 156.0 * u[IDX(ie-4,j,k)]
                                 + 214.0 * u[IDX(ie-3,j,k)]
                                 - 154.0 * u[IDX(ie-2,j,k)]
                                 +  45.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv644_yy(double * const  DyDyu, const double * const  u,
                 const double dy, const unsigned int *sz, unsigned bflag)
{

  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_180 = idy_sqrd / 180.0;
  const double idy_sqrd_by_12  = idy_sqrd /  12.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        DyDyu[pp] = (     2.0 * u[pp-3*nx]
                      -  27.0 * u[pp-2*nx]
                      + 270.0 * u[pp-  nx]
                      - 490.0 * u[pp     ]
                      + 270.0 * u[pp+  nx]
                      -  27.0 * u[pp+2*nx]
                      +   2.0 * u[pp+3*nx]
                    ) * idy_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        // The above two should be replaced by 4th order approximations: 
        DyDyu[IDX(i,3,k)] = (    45.0 * u[IDX(i,3,k)]
                              - 154.0 * u[IDX(i,4,k)]
                              + 214.0 * u[IDX(i,5,k)]
                              - 156.0 * u[IDX(i,6,k)]
                              +  61.0 * u[IDX(i,7,k)]
                              -  10.0 * u[IDX(i,8,k)]
                            ) * idy_sqrd_by_12;
        
        DyDyu[IDX(i,4,k)] = (   10.0 * u[IDX(i,3,k)]
                              - 15.0 * u[IDX(i,4,k)]
                              -  4.0 * u[IDX(i,5,k)]
                              + 14.0 * u[IDX(i,6,k)]
                              -  6.0 * u[IDX(i,7,k)]
                              +        u[IDX(i,8,k)]
                            ) * idy_sqrd_by_12;

        DyDyu[IDX(i,5,k)] = ( -        u[IDX(i,3,k)]
                              + 16.0 * u[IDX(i,4,k)]
                              - 30.0 * u[IDX(i,5,k)]
                              + 16.0 * u[IDX(i,6,k)]
                              -        u[IDX(i,7,k)]
                            ) * idy_sqrd_by_12;


        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        DyDyu[IDX(i,je-3,k)] = ( -        u[IDX(i,je-5,k)]
                                 + 16.0 * u[IDX(i,je-4,k)]
                                 - 30.0 * u[IDX(i,je-3,k)]
                                 + 16.0 * u[IDX(i,je-2,k)]
                                 -        u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

        // The above two should be replaced by 4th order approximations: 
        DyDyu[IDX(i,je-2,k)] = (          u[IDX(i,je-6,k)]
                                 -  6.0 * u[IDX(i,je-5,k)]
                                 + 14.0 * u[IDX(i,je-4,k)]
                                 -  4.0 * u[IDX(i,je-3,k)]
                                 - 15.0 * u[IDX(i,je-2,k)]
                                 + 10.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

        DyDyu[IDX(i,je-1,k)] = ( -  10.0 * u[IDX(i,je-6,k)]
                                 +  61.0 * u[IDX(i,je-5,k)]
                                 - 156.0 * u[IDX(i,je-4,k)]
                                 + 214.0 * u[IDX(i,je-3,k)]
                                 - 154.0 * u[IDX(i,je-2,k)]
                                 +  45.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
//void deriv64_zz(double * const  DzDzu, const double * const  u,
void deriv644_zz(double * const  DzDzu, const double * const  u,
                 const double dz, const unsigned int *sz, unsigned bflag)
{

  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_180 = idz_sqrd / 180.0;
  const double idz_sqrd_by_12  = idz_sqrd /  12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        DzDzu[pp] = (     2.0 * u[pp-3*n]
                      -  27.0 * u[pp-2*n]
                      + 270.0 * u[pp-  n]
                      - 490.0 * u[pp    ]
                      + 270.0 * u[pp+  n]
                      -  27.0 * u[pp+2*n]
                      +   2.0 * u[pp+3*n]
                    ) * idz_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        // The above two should be replaced by 4th order approximations: 
        DzDzu[IDX(i,j,3)] = (    45.0 * u[IDX(i,j,3)]
                              - 154.0 * u[IDX(i,j,4)]
                              + 214.0 * u[IDX(i,j,5)]
                              - 156.0 * u[IDX(i,j,6)]
                              +  61.0 * u[IDX(i,j,7)]
                              -  10.0 * u[IDX(i,j,8)]
                            ) * idz_sqrd_by_12;
        
        DzDzu[IDX(i,j,4)] = (    10.0 * u[IDX(i,j,3)]
                              -  15.0 * u[IDX(i,j,4)]
                              -   4.0 * u[IDX(i,j,5)]
                              +  14.0 * u[IDX(i,j,6)]
                              -   6.0 * u[IDX(i,j,7)]
                              +         u[IDX(i,j,8)]
                            ) * idz_sqrd_by_12;

        DzDzu[IDX(i,j,5)] = ( -         u[IDX(i,j,3)]
                              +  16.0 * u[IDX(i,j,4)]
                              -  30.0 * u[IDX(i,j,5)]
                              +  16.0 * u[IDX(i,j,6)]
                              -         u[IDX(i,j,7)]
                            ) * idz_sqrd_by_12;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        DzDzu[IDX(i,j,ke-3)] = ( -         u[IDX(i,j,ke-5)]
                                 +  16.0 * u[IDX(i,j,ke-4)]
                                 -  30.0 * u[IDX(i,j,ke-3)]
                                 +  16.0 * u[IDX(i,j,ke-2)]
                                 -         u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;
        // The above two should be replaced by 4th order approximations: 
        DzDzu[IDX(i,j,ke-2)] = (           u[IDX(i,j,ke-6)]
                                 -   6.0 * u[IDX(i,j,ke-5)]
                                 +  14.0 * u[IDX(i,j,ke-4)]
                                 -   4.0 * u[IDX(i,j,ke-3)]
                                 -  15.0 * u[IDX(i,j,ke-2)]
                                 +  10.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;
        
        DzDzu[IDX(i,j,ke-1)] = ( -  10.0 * u[IDX(i,j,ke-6)]
                                 +  61.0 * u[IDX(i,j,ke-5)]
                                 - 156.0 * u[IDX(i,j,ke-4)]
                                 + 214.0 * u[IDX(i,j,ke-3)]
                                 - 154.0 * u[IDX(i,j,ke-2)]
                                 +  45.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv644_xy( double * const DxDyu, double* const Dxu, 
                  const double * const u, const double dx, 
                  const double dy, const unsigned int *sz, unsigned bflag )
{
     deriv644_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv644_y( DxDyu, Dxu, dy, sz, bflag ) ; 
} 

void deriv644_xz( double * const DxDzu, double* const Dxu, 
                  const double * const u, const double dx, 
                  const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv644_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv644_z( DxDzu, Dxu, dz, sz, bflag ) ; 
} 


void deriv644_yz( double * const DyDzu, double* const Dyu, 
                  const double * const u, const double dy, 
                  const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv644_y( Dyu  ,   u, dy, sz, bflag ) ; 
    deriv644_z( DyDzu, Dyu, dz, sz, bflag ) ; 
} 




/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv642_x(double * const  Dxu, const double * const  u,
                const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx = 1.0 / dx;
  const double idx_by_2 = 0.5 * idx;
  const double idx_by_12 = idx / 12.0;
  const double idx_by_60 = idx / 60.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 0;
  const int kb = 0;
  const int ie = sz[0]-3;
  const int je = sz[1]-0;
  const int ke = sz[2]-0;
  
  const int n = 1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        Dxu[pp] = ( -        u[pp-3] 
                    +  9.0 * u[pp-2]
                    - 45.0 * u[pp-1]
                    + 45.0 * u[pp+1]
                    -  9.0 * u[pp+2]
                    +        u[pp+3] ) * idx_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
       
        // This is a (totally) shifted second order stencil.   
        Dxu[IDX(3,j,k)] = ( - 3.0 * u[IDX(3,j,k)]
                            + 4.0 * u[IDX(4,j,k)]
                            -       u[IDX(5,j,k)]
                          ) * idx_by_2;
        
        // This is a centered second order stencil.   
        Dxu[IDX(4,j,k)] = ( - u[IDX(3,j,k)]
                            + u[IDX(5,j,k)]
                          ) * idx_by_2;

        // This is a centered fourth order stencil.   
        Dxu[IDX(5,j,k)] = (         u[IDX(3,j,k)]
                            - 8.0 * u[IDX(4,j,k)]
                            + 8.0 * u[IDX(6,j,k)]
                            -       u[IDX(7,j,k)]
                          ) * idx_by_12;
        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered fourth order stencil.   
        Dxu[IDX(ie-3,j,k)] = (         u[IDX(ie-5,j,k)]
                               - 8.0 * u[IDX(ie-4,j,k)]
                               + 8.0 * u[IDX(ie-2,j,k)]
                               -       u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

        // This is a centered second order stencil.   
        Dxu[IDX(ie-2,j,k)] = ( - u[IDX(ie-3,j,k)]
                               + u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

        // This is a (totally) shifted second order stencil.   
        Dxu[IDX(ie-1,j,k)] = (         u[IDX(ie-3,j,k)]
                               - 4.0 * u[IDX(ie-2,j,k)]
                               + 3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif


}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv642_y(double * const  Dyu, const double * const  u,
                const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0 / dy;
  const double idy_by_2 = 0.5 * idy;
  const double idy_by_12 = idy / 12.0;
  const double idy_by_60 = idy / 60.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 0;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-0;
  
  const int n = nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        Dyu[pp] = ( -        u[pp-3*n] 
                    +  9.0 * u[pp-2*n]
                    - 45.0 * u[pp-  n]
                    + 45.0 * u[pp+  n]
                    -  9.0 * u[pp+2*n]
                    +        u[pp+3*n] ) * idy_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted second order stencil.   
        Dyu[IDX(i,3,k)] = ( - 3.0 * u[IDX(i,3,k)]
                            + 4.0 * u[IDX(i,4,k)]
                            -       u[IDX(i,5,k)]
                          ) * idy_by_2;

        // This is a centered second order stencil.   
        Dyu[IDX(i,4,k)] = ( - u[IDX(i,3,k)]
                            + u[IDX(i,5,k)]
                          ) * idy_by_2;

        // This is a centered fourth order stencil.   
        Dyu[IDX(i,5,k)] = (         u[IDX(i,3,k)]
                            - 8.0 * u[IDX(i,4,k)]
                            + 8.0 * u[IDX(i,6,k)]
                            -       u[IDX(i,7,k)]
                          ) * idy_by_12;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered fourth order stencil.   
        Dyu[IDX(i,je-3,k)] = (         u[IDX(i,je-5,k)]
                               - 8.0 * u[IDX(i,je-4,k)]
                               + 8.0 * u[IDX(i,je-2,k)]
                               -       u[IDX(i,je-1,k)]
                             ) * idy_by_12;

        // This is a centered second order stencil.   
        Dyu[IDX(i,je-2,k)] = ( - u[IDX(i,je-3,k)]
                               + u[IDX(i,je-1,k)]
                             ) * idy_by_2;
                             
        // This is a (totally) shifted second order stencil.   
        Dyu[IDX(i,je-1,k)] = (         u[IDX(i,je-3,k)]
                               - 4.0 * u[IDX(i,je-2,k)]
                               + 3.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_2;
      
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv642_z(double * const  Dzu, const double * const  u,
                const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0 / dz;
  const double idz_by_2 = 0.5 * idz;
  const double idz_by_12 = idz / 12.0;
  const double idz_by_60 = idz / 60.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0]-3;
  const int je = sz[1]-3;
  const int ke = sz[2]-3;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        Dzu[pp] = ( -        u[pp-3*n] 
                    +  9.0 * u[pp-2*n]
                    - 45.0 * u[pp-  n]
                    + 45.0 * u[pp+  n]
                    -  9.0 * u[pp+2*n]
                    +        u[pp+3*n] ) * idz_by_60;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted second order stencil.   
        Dzu[IDX(i,j,3)] = ( - 3.0 * u[IDX(i,j,3)]
                            + 4.0 * u[IDX(i,j,4)]
                            -       u[IDX(i,j,5)]
                          ) * idz_by_2;

        // This is a centered second order stencil.   
        Dzu[IDX(i,j,4)] = ( - u[IDX(i,j,3)]
                            + u[IDX(i,j,5)]
                          ) * idz_by_2;

        // This is a centered fourth order stencil.   
        Dzu[IDX(i,j,5)] = (         u[IDX(i,j,3)]
                            - 8.0 * u[IDX(i,j,4)]
                            + 8.0 * u[IDX(i,j,6)]
                            -       u[IDX(i,j,7)]
                          ) * idz_by_12;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered fourth order stencil.   
        Dzu[IDX(i,j,ke-3)] = (         u[IDX(i,j,ke-5)]
                               - 8.0 * u[IDX(i,j,ke-4)]
                               + 8.0 * u[IDX(i,j,ke-2)]
                               -       u[IDX(i,j,ke-1)]
                             ) * idz_by_12;

        // This is a centered second order stencil.   
        Dzu[IDX(i,j,ke-2)] = ( - u[IDX(i,j,ke-3)]
                               + u[IDX(i,j,ke-1)]
                             ) * idz_by_2;
                             
        // This is a (totally) shifted second order stencil.   
        Dzu[IDX(i,j,ke-1)] = (         u[IDX(i,j,ke-3)]
                               - 4.0 * u[IDX(i,j,ke-2)]
                               + 3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_2;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv642_xx(double * const  DxDxu, const double * const  u,
                 const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx_sqrd = 1.0 / (dx*dx);
  const double idx_sqrd_by_180 = idx_sqrd / 180.0;
  const double idx_sqrd_by_12  = idx_sqrd / 12.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);

        DxDxu[pp] = (     2.0 * u[pp-3]
                      -  27.0 * u[pp-2]
                      + 270.0 * u[pp-1]
                      - 490.0 * u[pp  ]
                      + 270.0 * u[pp+1]
                      -  27.0 * u[pp+2]
                      +   2.0 * u[pp+3]
                    ) * idx_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a (totally) shifted second order stencil.   
        DxDxu[IDX(3,j,k)] = (   2.0 * u[IDX(3,j,k)]
                              - 5.0 * u[IDX(4,j,k)]
                              + 4.0 * u[IDX(5,j,k)]
                              -       u[IDX(6,j,k)]
                            ) * idx_sqrd;
          
        // This is a centered second order stencil.   
        DxDxu[IDX(4,j,k)] = (         u[IDX(3,j,k)]
                              - 2.0 * u[IDX(4,j,k)]
                              +       u[IDX(5,j,k)]
                            ) * idx_sqrd;
        
        // This is a centered fourth order stencil.   
        DxDxu[IDX(5,j,k)] = ( -         u[IDX(3,j,k)]
                              +  16.0 * u[IDX(4,j,k)]
                              -  30.0 * u[IDX(5,j,k)]
                              +  16.0 * u[IDX(6,j,k)]
                              -         u[IDX(7,j,k)]
                            ) * idx_sqrd_by_12;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered fourth order stencil.   
        DxDxu[IDX(ie-3,j,k)] = ( -        u[IDX(ie-5,j,k)]
                                 + 16.0 * u[IDX(ie-4,j,k)]
                                 - 30.0 * u[IDX(ie-3,j,k)]
                                 + 16.0 * u[IDX(ie-2,j,k)]
                                 -        u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;
        
        // This is a centered second order stencil.   
        DxDxu[IDX(ie-2,j,k)] = (         u[IDX(ie-3,j,k)]
                                 - 2.0 * u[IDX(ie-2,j,k)]
                                 +       u[IDX(ie-1,j,k)]
                               ) * idx_sqrd;
         
        // This is a (totally) shifted second order stencil.   
        DxDxu[IDX(ie-1,j,k)] = ( -       u[IDX(ie-4,j,k)]
                                 + 4.0 * u[IDX(ie-3,j,k)]
                                 - 5.0 * u[IDX(ie-2,j,k)]
                                 + 2.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv642_yy(double * const  DyDyu, const double * const  u,
                 const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_180 = idy_sqrd / 180.0;
  const double idy_sqrd_by_12  = idy_sqrd /  12.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
        
        DyDyu[pp] = (     2.0 * u[pp-3*nx]
                      -  27.0 * u[pp-2*nx]
                      + 270.0 * u[pp-  nx]
                      - 490.0 * u[pp     ]
                      + 270.0 * u[pp+  nx]
                      -  27.0 * u[pp+2*nx]
                      +   2.0 * u[pp+3*nx]
                    ) * idy_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a (totally) shifted second order stencil.   
        DyDyu[IDX(i,3,k)] = (   2.0 * u[IDX(i,3,k)]
                              - 5.0 * u[IDX(i,4,k)]
                              + 4.0 * u[IDX(i,5,k)]
                              -       u[IDX(i,6,k)]
                            ) * idy_sqrd;
        
        // This is a centered second order stencil.   
        DyDyu[IDX(i,4,k)] = (         u[IDX(i,3,k)]
                              - 2.0 * u[IDX(i,4,k)]
                              +       u[IDX(i,5,k)]
                            ) * idy_sqrd;

        // This is a centered fourth order stencil.   
        DyDyu[IDX(i,5,k)] = ( -         u[IDX(i,3,k)]
                              +  16.0 * u[IDX(i,4,k)]
                              -  30.0 * u[IDX(i,5,k)]
                              +  16.0 * u[IDX(i,6,k)]
                              -         u[IDX(i,7,k)]
                            ) * idy_sqrd_by_12;
        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered fourth order stencil.   
        DyDyu[IDX(i,je-3,k)] = ( -        u[IDX(i,je-5,k)]
                                 + 16.0 * u[IDX(i,je-4,k)]
                                 - 30.0 * u[IDX(i,je-3,k)]
                                 + 16.0 * u[IDX(i,je-2,k)]
                                 -        u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

        // This is a centered second order stencil.   
        DyDyu[IDX(i,je-2,k)] = (         u[IDX(i,je-3,k)]
                                 - 2.0 * u[IDX(i,je-2,k)]
                                 +       u[IDX(i,je-1,k)]
                               ) * idy_sqrd; 

        // This is a (totally) shifted second order stencil.   
        DyDyu[IDX(i,je-1,k)] = ( -       u[IDX(i,je-4,k)]
                                 + 4.0 * u[IDX(i,je-3,k)]
                                 - 5.0 * u[IDX(i,je-2,k)]
                                 + 2.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv642_zz(double * const  DzDzu, const double * const  u,
                 const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_180 = idz_sqrd / 180.0;
  const double idz_sqrd_by_12  = idz_sqrd /  12.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 3;
  const int jb = 3;
  const int kb = 3;
  const int ie = sz[0] - 3;
  const int je = sz[1] - 3;
  const int ke = sz[2] - 3;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        
        DzDzu[pp] = (     2.0 * u[pp-3*n]
                      -  27.0 * u[pp-2*n]
                      + 270.0 * u[pp-  n]
                      - 490.0 * u[pp    ]
                      + 270.0 * u[pp+  n]
                      -  27.0 * u[pp+2*n]
                      +   2.0 * u[pp+3*n]
                    ) * idz_sqrd_by_180;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a (totally) shifted second order stencil.   
        DzDzu[IDX(i,j,3)] = (   2.0 * u[IDX(i,j,3)]
                              - 5.0 * u[IDX(i,j,4)]
                              + 4.0 * u[IDX(i,j,5)]
                              -       u[IDX(i,j,6)]
                            ) * idz_sqrd;
        
        // This is a centered second order stencil.   
        DzDzu[IDX(i,j,4)] = (         u[IDX(i,j,3)]
                              - 2.0 * u[IDX(i,j,4)]
                              +       u[IDX(i,j,5)]
                            ) * idz_sqrd;

        // This is a centered fourth order stencil.   
        DzDzu[IDX(i,j,5)] = ( -         u[IDX(i,j,3)]
                              +  16.0 * u[IDX(i,j,4)]
                              -  30.0 * u[IDX(i,j,5)]
                              +  16.0 * u[IDX(i,j,6)]
                              -         u[IDX(i,j,7)]
                            ) * idz_sqrd_by_12;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered fourth order stencil.   
        DzDzu[IDX(i,j,ke-3)] = ( -        u[IDX(i,j,ke-5)]
                                 + 16.0 * u[IDX(i,j,ke-4)]
                                 - 30.0 * u[IDX(i,j,ke-3)]
                                 + 16.0 * u[IDX(i,j,ke-2)]
                                 -        u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;

        // This is a centered second order stencil.   
        DzDzu[IDX(i,j,ke-2)] = (         u[IDX(i,j,ke-3)]
                                 - 2.0 * u[IDX(i,j,ke-2)]
                                 +       u[IDX(i,j,ke-1)]
                               ) * idz_sqrd;
        
        // This is a (totally) shifted second order stencil.   
        DzDzu[IDX(i,j,ke-1)] = ( -       u[IDX(i,j,ke-4)]
                                 + 4.0 * u[IDX(i,j,ke-3)]
                                 - 5.0 * u[IDX(i,j,ke-2)]
                                 + 2.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv642_xy( double * const DxDyu, double* const Dxu, 
                  const double * const u, const double dx, 
                  const double dy, const unsigned int *sz, unsigned bflag )
{
    deriv642_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv642_y( DxDyu, Dxu, dy, sz, bflag ) ; 
} 


void deriv642_xz( double * const DxDzu, double* const Dxu, 
                  const double * const u, const double dx, 
                  const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv642_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv642_z( DxDzu, Dxu, dz, sz, bflag ) ; 
} 


void deriv642_yz( double * const DyDzu, double* const Dyu, 
                  const double * const u, const double dy, 
                  const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv642_y( Dyu  ,   u, dy, sz, bflag ) ; 
    deriv642_z( DyDzu, Dyu, dz, sz, bflag ) ; 
} 



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8642_x(double * const  Dxu, const double * const  u,
                 const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx = 1.0 / dx;
  const double idx_by_2 = 0.5 * idx;
  const double idx_by_12 = idx / 12.0;
  const double idx_by_60 = idx / 60.0;
  const double idx_by_2520 = idx / 2520.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 0;
  const int kb = 0;
  const int ie = sz[0]-4;
  const int je = sz[1]-0;
  const int ke = sz[2]-0;
  
  const int n = 1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        
        Dxu[pp] = (      9.0 * u[pp-4] 
                    -   96.0 * u[pp-3]
                    +  504.0 * u[pp-2]
                    - 2016.0 * u[pp-1]
                    + 2016.0 * u[pp+1]
                    -  504.0 * u[pp+2]
                    +   96.0 * u[pp+3]
                    -    9.0 * u[pp+4] ) * idx_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
       
        // This is a (totally) shifted second order stencil.   
        Dxu[IDX(4,j,k)] = ( - 3.0 * u[IDX(4,j,k)]
                            + 4.0 * u[IDX(5,j,k)]
                            -       u[IDX(6,j,k)]
                          ) * idx_by_2;
        
        // This is a centered second order stencil.   
        Dxu[IDX(5,j,k)] = ( - u[IDX(4,j,k)]
                            + u[IDX(6,j,k)]
                          ) * idx_by_2;

        // This is a centered fourth order stencil.   
        Dxu[IDX(6,j,k)] = (         u[IDX(4,j,k)]
                            - 8.0 * u[IDX(5,j,k)]
                            + 8.0 * u[IDX(7,j,k)]
                            -       u[IDX(8,j,k)]
                          ) * idx_by_12;

        // This is a centered sixth order stencil.   
        Dxu[IDX(7,j,k)] = ( -        u[IDX( 4,j,k)]
                            +  9.0 * u[IDX( 5,j,k)]
                            - 45.0 * u[IDX( 6,j,k)]
                            + 45.0 * u[IDX( 8,j,k)]
                            -  9.0 * u[IDX( 9,j,k)]
                            +        u[IDX(10,j,k)]
                          ) * idx_by_60;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered sixth order stencil.   
        Dxu[IDX(ie-4,j,k)] = ( -        u[IDX(ie-7,j,k)]
                               +  9.0 * u[IDX(ie-6,j,k)]
                               - 45.0 * u[IDX(ie-5,j,k)]
                               + 45.0 * u[IDX(ie-3,j,k)]
                               -  9.0 * u[IDX(ie-2,j,k)]
                               +        u[IDX(ie-1,j,k)]
                             ) * idx_by_60;

        // This is a centered fourth order stencil.   
        Dxu[IDX(ie-3,j,k)] = (         u[IDX(ie-5,j,k)]
                               - 8.0 * u[IDX(ie-4,j,k)]
                               + 8.0 * u[IDX(ie-2,j,k)]
                               -       u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

        // This is a centered second order stencil.   
        Dxu[IDX(ie-2,j,k)] = ( - u[IDX(ie-3,j,k)]
                               + u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

        // This is a (totally) shifted second order stencil.   
        Dxu[IDX(ie-1,j,k)] = (         u[IDX(ie-3,j,k)]
                               - 4.0 * u[IDX(ie-2,j,k)]
                               + 3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_2;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif


}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8642_y(double * const  Dyu, const double * const  u,
                 const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0 / dy;
  const double idy_by_2 = 0.5 * idy;
  const double idy_by_12 = idy / 12.0;
  const double idy_by_60 = idy / 60.0;
  const double idy_by_2520 = idy / 2520.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 0;
  const int ie = sz[0]-4;
  const int je = sz[1]-4;
  const int ke = sz[2]-0;
  
  const int n = nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);

        Dyu[pp] = (      9.0 * u[pp-4*n] 
                    -   96.0 * u[pp-3*n]
                    +  504.0 * u[pp-2*n]
                    - 2016.0 * u[pp-  n]
                    + 2016.0 * u[pp+  n]
                    -  504.0 * u[pp+2*n]
                    +   96.0 * u[pp+3*n]
                    -    9.0 * u[pp+4*n] ) * idy_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted second order stencil.   
        Dyu[IDX(i,4,k)] = ( - 3.0 * u[IDX(i,4,k)]
                            + 4.0 * u[IDX(i,5,k)]
                            -       u[IDX(i,6,k)]
                          ) * idy_by_2;

        // This is a centered second order stencil.   
        Dyu[IDX(i,5,k)] = ( - u[IDX(i,4,k)]
                            + u[IDX(i,6,k)]
                          ) * idy_by_2;

        // This is a centered fourth order stencil.   
        Dyu[IDX(i,6,k)] = (         u[IDX(i,4,k)]
                            - 8.0 * u[IDX(i,5,k)]
                            + 8.0 * u[IDX(i,7,k)]
                            -       u[IDX(i,8,k)]
                          ) * idy_by_12;

        // This is a centered sixth order stencil.   
        Dyu[IDX(i,7,k)] = ( -        u[IDX(i, 4,k)]
                            +  9.0 * u[IDX(i, 5,k)]
                            - 45.0 * u[IDX(i, 6,k)]
                            + 45.0 * u[IDX(i, 8,k)]
                            -  9.0 * u[IDX(i, 9,k)]
                            +        u[IDX(i,10,k)]
                          ) * idy_by_60;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered sixth order stencil.   
        Dyu[IDX(i,je-4,k)] = ( -        u[IDX(i,je-7,k)]
                               +  9.0 * u[IDX(i,je-6,k)]
                               - 45.0 * u[IDX(i,je-5,k)]
                               + 45.0 * u[IDX(i,je-3,k)]
                               -  9.0 * u[IDX(i,je-2,k)]
                               +        u[IDX(i,je-1,k)]
                             ) * idy_by_60;

        // This is a centered fourth order stencil.   
        Dyu[IDX(i,je-3,k)] = (         u[IDX(i,je-5,k)]
                               - 8.0 * u[IDX(i,je-4,k)]
                               + 8.0 * u[IDX(i,je-2,k)]
                               -       u[IDX(i,je-1,k)]
                             ) * idy_by_12;

        // This is a centered second order stencil.   
        Dyu[IDX(i,je-2,k)] = ( - u[IDX(i,je-3,k)]
                               + u[IDX(i,je-1,k)]
                             ) * idy_by_2;
                             
        // This is a (totally) shifted second order stencil.   
        Dyu[IDX(i,je-1,k)] = (         u[IDX(i,je-3,k)]
                               - 4.0 * u[IDX(i,je-2,k)]
                               + 3.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_2;
      
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8642_z(double * const  Dzu, const double * const  u,
                 const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0 / dz;
  const double idz_by_2 = 0.5 * idz;
  const double idz_by_12 = idz / 12.0;
  const double idz_by_60 = idz / 60.0;
  const double idz_by_2520 = idz / 2520.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0]-4;
  const int je = sz[1]-4;
  const int ke = sz[2]-4;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);

        Dzu[pp] = (      9.0 * u[pp-4*n] 
                    -   96.0 * u[pp-3*n]
                    +  504.0 * u[pp-2*n]
                    - 2016.0 * u[pp-  n]
                    + 2016.0 * u[pp+  n]
                    -  504.0 * u[pp+2*n]
                    +   96.0 * u[pp+3*n]
                    -    9.0 * u[pp+4*n] ) * idz_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted second order stencil.   
        Dzu[IDX(i,j,4)] = ( - 3.0 * u[IDX(i,j,4)]
                            + 4.0 * u[IDX(i,j,5)]
                            -       u[IDX(i,j,6)]
                          ) * idz_by_2;

        // This is a centered second order stencil.   
        Dzu[IDX(i,j,5)] = ( - u[IDX(i,j,4)]
                            + u[IDX(i,j,6)]
                          ) * idz_by_2;

        // This is a centered fourth order stencil.   
        Dzu[IDX(i,j,6)] = (         u[IDX(i,j,4)]
                            - 8.0 * u[IDX(i,j,5)]
                            + 8.0 * u[IDX(i,j,7)]
                            -       u[IDX(i,j,8)]
                          ) * idz_by_12;

        // This is a centered sixth order stencil.   
        Dzu[IDX(i,j,7)] = ( -        u[IDX(i,j, 4)]
                            +  9.0 * u[IDX(i,j, 5)]
                            - 45.0 * u[IDX(i,j, 6)]
                            + 45.0 * u[IDX(i,j, 8)]
                            -  9.0 * u[IDX(i,j, 9)]
                            +        u[IDX(i,j,10)]
                          ) * idz_by_60;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered sixth order stencil.   
        Dzu[IDX(i,j,ke-4)] = ( -        u[IDX(i,j,ke-7)]
                               +  9.0 * u[IDX(i,j,ke-6)]
                               - 45.0 * u[IDX(i,j,ke-5)]
                               + 45.0 * u[IDX(i,j,ke-3)]
                               -  9.0 * u[IDX(i,j,ke-2)]
                               +        u[IDX(i,j,ke-1)]
                             ) * idz_by_60;

        // This is a centered fourth order stencil.   
        Dzu[IDX(i,j,ke-3)] = (         u[IDX(i,j,ke-5)]
                               - 8.0 * u[IDX(i,j,ke-4)]
                               + 8.0 * u[IDX(i,j,ke-2)]
                               -       u[IDX(i,j,ke-1)]
                             ) * idz_by_12;

        // This is a centered second order stencil.   
        Dzu[IDX(i,j,ke-2)] = ( - u[IDX(i,j,ke-3)]
                               + u[IDX(i,j,ke-1)]
                             ) * idz_by_2;
                             
        // This is a (totally) shifted second order stencil.   
        Dzu[IDX(i,j,ke-1)] = (         u[IDX(i,j,ke-3)]
                               - 4.0 * u[IDX(i,j,ke-2)]
                               + 3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_2;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8642_xx(double * const  DxDxu, const double * const  u,
                  const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx_sqrd = 1.0 / (dx*dx);
  const double idx_sqrd_by_12  = idx_sqrd / 12.0;
  const double idx_sqrd_by_180 = idx_sqrd / 180.0;
  const double idx_sqrd_by_5040 = idx_sqrd / 5040.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);

        DxDxu[pp] = ( -     9.0 * u[pp-4]
                      +   128.0 * u[pp-3]
                      -  1008.0 * u[pp-2]
                      +  8064.0 * u[pp-1]
                      - 14350.0 * u[pp  ]
                      +  8064.0 * u[pp+1]
                      -  1008.0 * u[pp+2]
                      +   128.0 * u[pp+3]
                      -     9.0 * u[pp+4]
                    ) * idx_sqrd_by_5040;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {

        // This is a (totally) shifted second order stencil.   
        DxDxu[IDX(4,j,k)] = (   2.0 * u[IDX(4,j,k)]
                              - 5.0 * u[IDX(5,j,k)]
                              + 4.0 * u[IDX(6,j,k)]
                              -       u[IDX(7,j,k)]
                            ) * idx_sqrd;
          
        // This is a centered second order stencil.   
        DxDxu[IDX(5,j,k)] = (         u[IDX(4,j,k)]
                              - 2.0 * u[IDX(5,j,k)]
                              +       u[IDX(6,j,k)]
                            ) * idx_sqrd;
        
        // This is a centered fourth order stencil.   
        DxDxu[IDX(6,j,k)] = ( -         u[IDX(4,j,k)]
                              +  16.0 * u[IDX(5,j,k)]
                              -  30.0 * u[IDX(6,j,k)]
                              +  16.0 * u[IDX(7,j,k)]
                              -         u[IDX(8,j,k)]
                            ) * idx_sqrd_by_12;
        
        // This is a centered sixth order stencil.   
        DxDxu[IDX(7,j,k)] = (      2.0 * u[IDX( 4,j,k)]
                              -   27.0 * u[IDX( 5,j,k)]
                              +  270.0 * u[IDX( 6,j,k)]
                              -  490.0 * u[IDX( 7,j,k)]
                              +  270.0 * u[IDX( 8,j,k)]
                              -   27.0 * u[IDX( 9,j,k)]
                              +    2.0 * u[IDX(10,j,k)]
                            ) * idx_sqrd_by_180;
     
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered sixth order stencil.   
        DxDxu[IDX(ie-4,j,k)] = (     2.0 * u[IDX(ie-7,j,k)]
                                 -  27.0 * u[IDX(ie-6,j,k)]
                                 + 270.0 * u[IDX(ie-5,j,k)]
                                 - 490.0 * u[IDX(ie-4,j,k)]
                                 + 270.0 * u[IDX(ie-3,j,k)]
                                 -  27.0 * u[IDX(ie-2,j,k)]
                                 +   2.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_180;
        
        // This is a centered fourth order stencil.   
        DxDxu[IDX(ie-3,j,k)] = ( -        u[IDX(ie-5,j,k)]
                                 + 16.0 * u[IDX(ie-4,j,k)]
                                 - 30.0 * u[IDX(ie-3,j,k)]
                                 + 16.0 * u[IDX(ie-2,j,k)]
                                 -        u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;
        
        // This is a centered second order stencil.   
        DxDxu[IDX(ie-2,j,k)] = (         u[IDX(ie-3,j,k)]
                                 - 2.0 * u[IDX(ie-2,j,k)]
                                 +       u[IDX(ie-1,j,k)]
                               ) * idx_sqrd;
         
        // This is a (totally) shifted second order stencil.   
        DxDxu[IDX(ie-1,j,k)] = ( -       u[IDX(ie-4,j,k)]
                                 + 4.0 * u[IDX(ie-3,j,k)]
                                 - 5.0 * u[IDX(ie-2,j,k)]
                                 + 2.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8642_yy(double * const  DyDyu, const double * const  u,
                  const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_12  = idy_sqrd /  12.0;
  const double idy_sqrd_by_180 = idy_sqrd / 180.0;
  const double idy_sqrd_by_5040 = idy_sqrd / 5040.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
 
        DyDyu[pp] = ( -     9.0 * u[pp-4*nx]
                      +   128.0 * u[pp-3*nx]
                      -  1008.0 * u[pp-2*nx]
                      +  8064.0 * u[pp-  nx]
                      - 14350.0 * u[pp     ]
                      +  8064.0 * u[pp+  nx]
                      -  1008.0 * u[pp+2*nx]
                      +   128.0 * u[pp+3*nx]
                      -     9.0 * u[pp+4*nx]
                    ) * idy_sqrd_by_5040;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted second order stencil.   
        DyDyu[IDX(i,4,k)] = (   2.0 * u[IDX(i,4,k)]
                              - 5.0 * u[IDX(i,5,k)]
                              + 4.0 * u[IDX(i,6,k)]
                              -       u[IDX(i,7,k)]
                            ) * idy_sqrd;
        
        // This is a centered second order stencil.   
        DyDyu[IDX(i,5,k)] = (         u[IDX(i,4,k)]
                              - 2.0 * u[IDX(i,5,k)]
                              +       u[IDX(i,6,k)]
                            ) * idy_sqrd;

        // This is a centered fourth order stencil.   
        DyDyu[IDX(i,6,k)] = ( -         u[IDX(i,4,k)]
                              +  16.0 * u[IDX(i,5,k)]
                              -  30.0 * u[IDX(i,6,k)]
                              +  16.0 * u[IDX(i,7,k)]
                              -         u[IDX(i,8,k)]
                            ) * idy_sqrd_by_12;
        
        // This is a centered sixth order stencil.   
        DyDyu[IDX(i,7,k)] = (      2.0 * u[IDX(i, 4,k)]
                              -   27.0 * u[IDX(i, 5,k)]
                              +  270.0 * u[IDX(i, 6,k)]
                              -  490.0 * u[IDX(i, 7,k)]
                              +  270.0 * u[IDX(i, 8,k)]
                              -   27.0 * u[IDX(i, 9,k)]
                              +    2.0 * u[IDX(i,10,k)]
                            ) * idy_sqrd_by_180;
        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a centered sixth order stencil.   
        DyDyu[IDX(i,je-4,k)] = (     2.0 * u[IDX(i,je-7,k)]
                                 -  27.0 * u[IDX(i,je-6,k)]
                                 + 270.0 * u[IDX(i,je-5,k)]
                                 - 490.0 * u[IDX(i,je-4,k)]
                                 + 270.0 * u[IDX(i,je-3,k)]
                                 -  27.0 * u[IDX(i,je-2,k)]
                                 +   2.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_180;
        
        // This is a centered fourth order stencil.   
        DyDyu[IDX(i,je-3,k)] = ( -        u[IDX(i,je-5,k)]
                                 + 16.0 * u[IDX(i,je-4,k)]
                                 - 30.0 * u[IDX(i,je-3,k)]
                                 + 16.0 * u[IDX(i,je-2,k)]
                                 -        u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

        // This is a centered second order stencil.   
        DyDyu[IDX(i,je-2,k)] = (         u[IDX(i,je-3,k)]
                                 - 2.0 * u[IDX(i,je-2,k)]
                                 +       u[IDX(i,je-1,k)]
                               ) * idy_sqrd; 

        // This is a (totally) shifted second order stencil.   
        DyDyu[IDX(i,je-1,k)] = ( -       u[IDX(i,je-4,k)]
                                 + 4.0 * u[IDX(i,je-3,k)]
                                 - 5.0 * u[IDX(i,je-2,k)]
                                 + 2.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8642_zz(double * const  DzDzu, const double * const  u,
                  const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_12  = idz_sqrd /  12.0;
  const double idz_sqrd_by_180 = idz_sqrd / 180.0;
  const double idz_sqrd_by_5040 = idz_sqrd / 5040.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);

        DzDzu[pp] = ( -     9.0 * u[pp-4*n]
                      +   128.0 * u[pp-3*n]
                      -  1008.0 * u[pp-2*n]
                      +  8064.0 * u[pp-  n]
                      - 14350.0 * u[pp    ]
                      +  8064.0 * u[pp+  n]
                      -  1008.0 * u[pp+2*n]
                      +   128.0 * u[pp+3*n]
                      -     9.0 * u[pp+4*n]
                    ) * idz_sqrd_by_5040;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a (totally) shifted second order stencil.   
        DzDzu[IDX(i,j,4)] = (   2.0 * u[IDX(i,j,4)]
                              - 5.0 * u[IDX(i,j,5)]
                              + 4.0 * u[IDX(i,j,6)]
                              -       u[IDX(i,j,7)]
                            ) * idz_sqrd;
        
        // This is a centered second order stencil.   
        DzDzu[IDX(i,j,5)] = (         u[IDX(i,j,4)]
                              - 2.0 * u[IDX(i,j,5)]
                              +       u[IDX(i,j,6)]
                            ) * idz_sqrd;

        // This is a centered fourth order stencil.   
        DzDzu[IDX(i,j,6)] = ( -         u[IDX(i,j,4)]
                              +  16.0 * u[IDX(i,j,5)]
                              -  30.0 * u[IDX(i,j,6)]
                              +  16.0 * u[IDX(i,j,7)]
                              -         u[IDX(i,j,8)]
                            ) * idz_sqrd_by_12;
        
        // This is a centered sixth order stencil.   
        DzDzu[IDX(i,j,7)] = (      2.0 * u[IDX(i,j, 4)]
                              -   27.0 * u[IDX(i,j, 5)]
                              +  270.0 * u[IDX(i,j, 6)]
                              -  490.0 * u[IDX(i,j, 7)]
                              +  270.0 * u[IDX(i,j, 8)]
                              -   27.0 * u[IDX(i,j, 9)]
                              +    2.0 * u[IDX(i,j,10)]
                            ) * idz_sqrd_by_180;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered sixth order stencil.   
        DzDzu[IDX(i,j,ke-4)] = (     2.0 * u[IDX(i,j,ke-7)]
                                 -  27.0 * u[IDX(i,j,ke-6)]
                                 + 270.0 * u[IDX(i,j,ke-5)]
                                 - 490.0 * u[IDX(i,j,ke-4)]
                                 + 270.0 * u[IDX(i,j,ke-3)]
                                 -  27.0 * u[IDX(i,j,ke-2)]
                                 +   2.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_180;
        
        // This is a centered fourth order stencil.   
        DzDzu[IDX(i,j,ke-3)] = ( -        u[IDX(i,j,ke-5)]
                                 + 16.0 * u[IDX(i,j,ke-4)]
                                 - 30.0 * u[IDX(i,j,ke-3)]
                                 + 16.0 * u[IDX(i,j,ke-2)]
                                 -        u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;

        // This is a centered second order stencil.   
        DzDzu[IDX(i,j,ke-2)] = (         u[IDX(i,j,ke-3)]
                                 - 2.0 * u[IDX(i,j,ke-2)]
                                 +       u[IDX(i,j,ke-1)]
                               ) * idz_sqrd;
        
        // This is a (totally) shifted second order stencil.   
        DzDzu[IDX(i,j,ke-1)] = ( -       u[IDX(i,j,ke-4)]
                                 + 4.0 * u[IDX(i,j,ke-3)]
                                 - 5.0 * u[IDX(i,j,ke-2)]
                                 + 2.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd;

      }
    }
  }

  
  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8642_xy( double * const DxDyu, double* const Dxu, 
                   const double * const u, const double dx, 
                   const double dy, const unsigned int *sz, unsigned bflag )
{
    deriv8642_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv8642_y( DxDyu, Dxu, dy, sz, bflag ) ; 
} 


void deriv8642_xz( double * const DxDzu, double* const Dxu, 
                   const double * const u, const double dx, 
                   const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv8642_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv8642_z( DxDzu, Dxu, dz, sz, bflag ) ; 
} 


void deriv8642_yz( double * const DyDzu, double* const Dyu, 
                   const double * const u, const double dy, 
                   const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv8642_y( Dyu  ,   u, dy, sz, bflag ) ; 
    deriv8642_z( DyDzu, Dyu, dz, sz, bflag ) ; 
} 




/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8664_x(double * const  Dxu, const double * const  u,
                 const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx = 1.0 / dx;
  const double idx_by_2 = 0.5 * idx;
  const double idx_by_12 = idx / 12.0;
  const double idx_by_60 = idx / 60.0;
  const double idx_by_2520 = idx / 2520.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 0;
  const int kb = 0;
  const int ie = sz[0]-4;
  const int je = sz[1]-0;
  const int ke = sz[2]-0;
  
  const int n = 1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        
        Dxu[pp] = (      9.0 * u[pp-4] 
                    -   96.0 * u[pp-3]
                    +  504.0 * u[pp-2]
                    - 2016.0 * u[pp-1]
                    + 2016.0 * u[pp+1]
                    -  504.0 * u[pp+2]
                    +   96.0 * u[pp+3]
                    -    9.0 * u[pp+4] ) * idx_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
       
        // This is a (totally) shifted fourth order stencil.   
        Dxu[IDX(4,j,k)] = ( - 25.0 * u[IDX(4,j,k)]
                            + 48.0 * u[IDX(5,j,k)]
                            - 36.0 * u[IDX(6,j,k)]
                            + 16.0 * u[IDX(7,j,k)]
                            -  3.0 * u[IDX(8,j,k)]
                          ) * idx_by_12;
        
        // This is a (partially) shifted fourth order stencil.   
        Dxu[IDX(5,j,k)] = ( -  3.0 * u[IDX(4,j,k)]
                            - 10.0 * u[IDX(5,j,k)]
                            + 18.0 * u[IDX(6,j,k)]
                            -  6.0 * u[IDX(7,j,k)]
                            +        u[IDX(8,j,k)]
                          ) * idx_by_12;

        // This is a shifted sixth order stencil.   
        Dxu[IDX(6,j,k)] = (    2.0 * u[IDX( 4,j,k)]
                            - 24.0 * u[IDX( 5,j,k)]
                            - 35.0 * u[IDX( 6,j,k)]
                            + 80.0 * u[IDX( 7,j,k)]
                            - 30.0 * u[IDX( 8,j,k)]
                            +  8.0 * u[IDX( 9,j,k)]
                            -        u[IDX(10,j,k)]
                          ) * idx_by_60;

        // This is a centered sixth order stencil.   
        Dxu[IDX(7,j,k)] = ( -        u[IDX( 4,j,k)]
                            +  9.0 * u[IDX( 5,j,k)]
                            - 45.0 * u[IDX( 6,j,k)]
                            + 45.0 * u[IDX( 8,j,k)]
                            -  9.0 * u[IDX( 9,j,k)]
                            +        u[IDX(10,j,k)]
                          ) * idx_by_60;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered sixth order stencil.   
        Dxu[IDX(ie-4,j,k)] = ( -        u[IDX(ie-7,j,k)]
                               +  9.0 * u[IDX(ie-6,j,k)]
                               - 45.0 * u[IDX(ie-5,j,k)]
                               + 45.0 * u[IDX(ie-3,j,k)]
                               -  9.0 * u[IDX(ie-2,j,k)]
                               +        u[IDX(ie-1,j,k)]
                             ) * idx_by_60;

        // This is a shifted sixth order stencil.   
        Dxu[IDX(ie-3,j,k)] = (          u[IDX(ie-7,j,k)]
                               -  8.0 * u[IDX(ie-6,j,k)]
                               + 30.0 * u[IDX(ie-5,j,k)]
                               - 80.0 * u[IDX(ie-4,j,k)]
                               + 35.0 * u[IDX(ie-3,j,k)]
                               + 24.0 * u[IDX(ie-2,j,k)]
                               -  2.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_60;

        // This is a (partially) shifted fourth order stencil.   
        Dxu[IDX(ie-2,j,k)] = ( -        u[IDX(ie-5,j,k)]
                               +  6.0 * u[IDX(ie-4,j,k)]
                               - 18.0 * u[IDX(ie-3,j,k)]
                               + 10.0 * u[IDX(ie-2,j,k)]
                               +  3.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

        // This is a (totally) shifted fourth order stencil.   
        Dxu[IDX(ie-1,j,k)] = (    3.0 * u[IDX(ie-5,j,k)]
                               - 16.0 * u[IDX(ie-4,j,k)]
                               + 36.0 * u[IDX(ie-3,j,k)]
                               - 48.0 * u[IDX(ie-2,j,k)]
                               + 25.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif


}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8664_y(double * const  Dyu, const double * const  u,
                 const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0 / dy;
  const double idy_by_2 = 0.5 * idy;
  const double idy_by_12 = idy / 12.0;
  const double idy_by_60 = idy / 60.0;
  const double idy_by_2520 = idy / 2520.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 0;
  const int ie = sz[0]-4;
  const int je = sz[1]-4;
  const int ke = sz[2]-0;
  
  const int n = nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);

        Dyu[pp] = (      9.0 * u[pp-4*n] 
                    -   96.0 * u[pp-3*n]
                    +  504.0 * u[pp-2*n]
                    - 2016.0 * u[pp-  n]
                    + 2016.0 * u[pp+  n]
                    -  504.0 * u[pp+2*n]
                    +   96.0 * u[pp+3*n]
                    -    9.0 * u[pp+4*n] ) * idy_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted fourth order stencil.   
        Dyu[IDX(i,4,k)] = ( - 25.0 * u[IDX(i,4,k)]
                            + 48.0 * u[IDX(i,5,k)]
                            - 36.0 * u[IDX(i,6,k)]
                            + 16.0 * u[IDX(i,7,k)]
                            -  3.0 * u[IDX(i,8,k)]
                          ) * idy_by_12;

        // This is a (partially) shifted fourth order stencil.   
        Dyu[IDX(i,5,k)] = ( -  3.0 * u[IDX(i,4,k)]
                            - 10.0 * u[IDX(i,5,k)]
                            + 18.0 * u[IDX(i,6,k)]
                            -  6.0 * u[IDX(i,7,k)]
                            +        u[IDX(i,8,k)]
                          ) * idy_by_12;

        // This is a shifted sixth order stencil.   
        Dyu[IDX(i,6,k)] = (    2.0 * u[IDX(i, 4,k)]
                            - 24.0 * u[IDX(i, 5,k)]
                            - 35.0 * u[IDX(i, 6,k)]
                            + 80.0 * u[IDX(i, 7,k)]
                            - 30.0 * u[IDX(i, 8,k)]
                            +  8.0 * u[IDX(i, 9,k)]
                            -        u[IDX(i,10,k)]
                          ) * idy_by_60;

        // This is a centered sixth order stencil.   
        Dyu[IDX(i,7,k)] = ( -        u[IDX(i, 4,k)]
                            +  9.0 * u[IDX(i, 5,k)]
                            - 45.0 * u[IDX(i, 6,k)]
                            + 45.0 * u[IDX(i, 8,k)]
                            -  9.0 * u[IDX(i, 9,k)]
                            +        u[IDX(i,10,k)]
                          ) * idy_by_60;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a centered sixth order stencil.   
        Dyu[IDX(i,je-4,k)] = ( -        u[IDX(i,je-7,k)]
                               +  9.0 * u[IDX(i,je-6,k)]
                               - 45.0 * u[IDX(i,je-5,k)]
                               + 45.0 * u[IDX(i,je-3,k)]
                               -  9.0 * u[IDX(i,je-2,k)]
                               +        u[IDX(i,je-1,k)]
                             ) * idy_by_60;

        // This is a shifted sixth order stencil.   
        Dyu[IDX(i,je-3,k)] = (          u[IDX(i,je-7,k)]
                               -  8.0 * u[IDX(i,je-6,k)]
                               + 30.0 * u[IDX(i,je-5,k)]
                               - 80.0 * u[IDX(i,je-4,k)]
                               + 35.0 * u[IDX(i,je-3,k)]
                               + 24.0 * u[IDX(i,je-2,k)]
                               -  2.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_60;

        // This is a (partially) shifted fourth order stencil.   
        Dyu[IDX(i,je-2,k)] = ( -        u[IDX(i,je-5,k)]
                               +  6.0 * u[IDX(i,je-4,k)]
                               - 18.0 * u[IDX(i,je-3,k)]
                               + 10.0 * u[IDX(i,je-2,k)]
                               +  3.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;
                             
        // This is a (totally) shifted fourth order stencil.   
        Dyu[IDX(i,je-1,k)] = (    3.0 * u[IDX(i,je-5,k)]
                               - 16.0 * u[IDX(i,je-4,k)]
                               + 36.0 * u[IDX(i,je-3,k)]
                               - 48.0 * u[IDX(i,je-2,k)]
                               + 25.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;
      
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8664_z(double * const  Dzu, const double * const  u,
                 const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0 / dz;
  const double idz_by_2 = 0.5 * idz;
  const double idz_by_12 = idz / 12.0;
  const double idz_by_60 = idz / 60.0;
  const double idz_by_2520 = idz / 2520.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0]-4;
  const int je = sz[1]-4;
  const int ke = sz[2]-4;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);

        Dzu[pp] = (      9.0 * u[pp-4*n] 
                    -   96.0 * u[pp-3*n]
                    +  504.0 * u[pp-2*n]
                    - 2016.0 * u[pp-  n]
                    + 2016.0 * u[pp+  n]
                    -  504.0 * u[pp+2*n]
                    +   96.0 * u[pp+3*n]
                    -    9.0 * u[pp+4*n] ) * idz_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted fourth order stencil.   
        Dzu[IDX(i,j,4)] = ( - 25.0 * u[IDX(i,j,4)]
                            + 48.0 * u[IDX(i,j,5)]
                            - 36.0 * u[IDX(i,j,6)]
                            + 16.0 * u[IDX(i,j,7)]
                            -  3.0 * u[IDX(i,j,8)]
                          ) * idz_by_12;

        // This is a (partially) shifted fourth order stencil.   
        Dzu[IDX(i,j,5)] = ( -  3.0 * u[IDX(i,j,4)]
                            - 10.0 * u[IDX(i,j,5)]
                            + 18.0 * u[IDX(i,j,6)]
                            -  6.0 * u[IDX(i,j,7)]
                            +        u[IDX(i,j,8)]
                          ) * idz_by_12;

        // This is a shifted sixth order stencil.   
        Dzu[IDX(i,j,6)] = (    2.0 * u[IDX(i,j, 4)]
                            - 24.0 * u[IDX(i,j, 5)]
                            - 35.0 * u[IDX(i,j, 6)]
                            + 80.0 * u[IDX(i,j, 7)]
                            - 30.0 * u[IDX(i,j, 8)]
                            +  8.0 * u[IDX(i,j, 9)]
                            -        u[IDX(i,j,10)]
                          ) * idz_by_60;

        // This is a centered sixth order stencil.   
        Dzu[IDX(i,j,7)] = ( -        u[IDX(i,j, 4)]
                            +  9.0 * u[IDX(i,j, 5)]
                            - 45.0 * u[IDX(i,j, 6)]
                            + 45.0 * u[IDX(i,j, 8)]
                            -  9.0 * u[IDX(i,j, 9)]
                            +        u[IDX(i,j,10)]
                          ) * idz_by_60;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a centered sixth order stencil.   
        Dzu[IDX(i,j,ke-4)] = ( -        u[IDX(i,j,ke-7)]
                               +  9.0 * u[IDX(i,j,ke-6)]
                               - 45.0 * u[IDX(i,j,ke-5)]
                               + 45.0 * u[IDX(i,j,ke-3)]
                               -  9.0 * u[IDX(i,j,ke-2)]
                               +        u[IDX(i,j,ke-1)]
                             ) * idz_by_60;

        // This is a shifted sixth order stencil.   
        Dzu[IDX(i,j,ke-3)] = (          u[IDX(i,j,ke-7)]
                               -  8.0 * u[IDX(i,j,ke-6)]
                               + 30.0 * u[IDX(i,j,ke-5)]
                               - 80.0 * u[IDX(i,j,ke-4)]
                               + 35.0 * u[IDX(i,j,ke-3)]
                               + 24.0 * u[IDX(i,j,ke-2)]
                               -  2.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_60;

        // This is a (partially) shifted fourth order stencil.   
        Dzu[IDX(i,j,ke-2)] = ( -        u[IDX(i,j,ke-5)]
                               +  6.0 * u[IDX(i,j,ke-4)]
                               - 18.0 * u[IDX(i,j,ke-3)]
                               + 10.0 * u[IDX(i,j,ke-2)]
                               +  3.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_12;
                             
        // This is a (totally) shifted fourth order stencil.   
        Dzu[IDX(i,j,ke-1)] = (    3.0 * u[IDX(i,j,ke-5)]
                               - 16.0 * u[IDX(i,j,ke-4)]
                               + 36.0 * u[IDX(i,j,ke-3)]
                               - 48.0 * u[IDX(i,j,ke-2)]
                               + 25.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8664_xx(double * const  DxDxu, const double * const  u,
                  const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx_sqrd = 1.0 / (dx*dx);
  const double idx_sqrd_by_12  = idx_sqrd / 12.0;
  const double idx_sqrd_by_180 = idx_sqrd / 180.0;
  const double idx_sqrd_by_5040 = idx_sqrd / 5040.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);

        DxDxu[pp] = ( -     9.0 * u[pp-4]
                      +   128.0 * u[pp-3]
                      -  1008.0 * u[pp-2]
                      +  8064.0 * u[pp-1]
                      - 14350.0 * u[pp  ]
                      +  8064.0 * u[pp+1]
                      -  1008.0 * u[pp+2]
                      +   128.0 * u[pp+3]
                      -     9.0 * u[pp+4]
                    ) * idx_sqrd_by_5040;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a (totally) shifted fourth order stencil.   
        DxDxu[IDX(4,j,k)] = (    45.0 * u[IDX(4,j,k)]
                              - 154.0 * u[IDX(5,j,k)]
                              + 214.0 * u[IDX(6,j,k)]
                              - 156.0 * u[IDX(7,j,k)]
                              +  61.0 * u[IDX(8,j,k)]
                              -  10.0 * u[IDX(9,j,k)]
                            ) * idx_sqrd_by_12;
          
        // This is a (partially) shifted fourth order stencil.   
        DxDxu[IDX(5,j,k)] = (   10.0 * u[IDX(4,j,k)]
                              - 15.0 * u[IDX(5,j,k)]
                              -  4.0 * u[IDX(6,j,k)]
                              + 14.0 * u[IDX(7,j,k)]
                              -  6.0 * u[IDX(8,j,k)]
                              +        u[IDX(9,j,k)]
                            ) * idx_sqrd_by_12;
        
        // This is a shifted sixth order stencil.   
        DxDxu[IDX(6,j,k)] = ( -  11.0 * u[IDX( 4,j,k)]
                              + 214.0 * u[IDX( 5,j,k)]
                              - 378.0 * u[IDX( 6,j,k)]
                              + 130.0 * u[IDX( 7,j,k)]
                              +  85.0 * u[IDX( 8,j,k)]
                              -  54.0 * u[IDX( 9,j,k)]
                              +  16.0 * u[IDX(10,j,k)]
                              -   2.0 * u[IDX(11,j,k)]
                            ) * idx_sqrd_by_180;

        // This is a centered sixth order stencil.   
        DxDxu[IDX(7,j,k)] = (      2.0 * u[IDX( 4,j,k)]
                              -   27.0 * u[IDX( 5,j,k)]
                              +  270.0 * u[IDX( 6,j,k)]
                              -  490.0 * u[IDX( 7,j,k)]
                              +  270.0 * u[IDX( 8,j,k)]
                              -   27.0 * u[IDX( 9,j,k)]
                              +    2.0 * u[IDX(10,j,k)]
                            ) * idx_sqrd_by_180;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered sixth order stencil.   
        DxDxu[IDX(ie-4,j,k)] = (     2.0 * u[IDX(ie-7,j,k)]
                                 -  27.0 * u[IDX(ie-6,j,k)]
                                 + 270.0 * u[IDX(ie-5,j,k)]
                                 - 490.0 * u[IDX(ie-4,j,k)]
                                 + 270.0 * u[IDX(ie-3,j,k)]
                                 -  27.0 * u[IDX(ie-2,j,k)]
                                 +   2.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_180;
        
        // This is a shifted sixth order stencil.   
        DxDxu[IDX(ie-3,j,k)] = ( -   2.0 * u[IDX(ie-8,j,k)]
                                 +  16.0 * u[IDX(ie-7,j,k)]
                                 -  54.0 * u[IDX(ie-6,j,k)]
                                 +  85.0 * u[IDX(ie-5,j,k)]
                                 + 130.0 * u[IDX(ie-4,j,k)]
                                 - 378.0 * u[IDX(ie-3,j,k)]
                                 + 214.0 * u[IDX(ie-2,j,k)]
                                 -  11.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_180;
        
        // This is a (partially) shifted fourth order stencil.   
        DxDxu[IDX(ie-2,j,k)] = (          u[IDX(ie-6,j,k)]
                                 -  6.0 * u[IDX(ie-5,j,k)]
                                 + 14.0 * u[IDX(ie-4,j,k)]
                                 -  4.0 * u[IDX(ie-3,j,k)]
                                 - 15.0 * u[IDX(ie-2,j,k)]
                                 + 10.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;
         
        // XThis is a (totally) shifted fourth order stencil.   
        DxDxu[IDX(ie-1,j,k)] = ( -  10.0 * u[IDX(ie-6,j,k)]
                                 +  61.0 * u[IDX(ie-5,j,k)]
                                 - 156.0 * u[IDX(ie-4,j,k)]
                                 + 214.0 * u[IDX(ie-3,j,k)]
                                 - 154.0 * u[IDX(ie-2,j,k)]
                                 +  45.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8664_yy(double * const  DyDyu, const double * const  u,
                  const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_12  = idy_sqrd /  12.0;
  const double idy_sqrd_by_180 = idy_sqrd / 180.0;
  const double idy_sqrd_by_5040 = idy_sqrd / 5040.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
 
        DyDyu[pp] = ( -     9.0 * u[pp-4*nx]
                      +   128.0 * u[pp-3*nx]
                      -  1008.0 * u[pp-2*nx]
                      +  8064.0 * u[pp-  nx]
                      - 14350.0 * u[pp     ]
                      +  8064.0 * u[pp+  nx]
                      -  1008.0 * u[pp+2*nx]
                      +   128.0 * u[pp+3*nx]
                      -     9.0 * u[pp+4*nx]
                    ) * idy_sqrd_by_5040;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted fourth order stencil.   
        DyDyu[IDX(i,4,k)] = (    45.0 * u[IDX(i,4,k)]
                              - 154.0 * u[IDX(i,5,k)]
                              + 214.0 * u[IDX(i,6,k)]
                              - 156.0 * u[IDX(i,7,k)]
                              +  61.0 * u[IDX(i,8,k)]
                              -  10.0 * u[IDX(i,9,k)]
                            ) * idy_sqrd_by_12;
        
        // This is a (partially) shifted fourth order stencil.   
        DyDyu[IDX(i,5,k)] = (   10.0 * u[IDX(i,4,k)]
                              - 15.0 * u[IDX(i,5,k)]
                              -  4.0 * u[IDX(i,6,k)]
                              + 14.0 * u[IDX(i,7,k)]
                              -  6.0 * u[IDX(i,8,k)]
                              +        u[IDX(i,9,k)]
                            ) * idy_sqrd_by_12;

        // This is a shifted sixth order stencil.   
        DyDyu[IDX(i,6,k)] = ( -  11.0 * u[IDX(i, 4,k)]
                              + 214.0 * u[IDX(i, 5,k)]
                              - 378.0 * u[IDX(i, 6,k)]
                              + 130.0 * u[IDX(i, 7,k)]
                              +  85.0 * u[IDX(i, 8,k)]
                              -  54.0 * u[IDX(i, 9,k)]
                              +  16.0 * u[IDX(i,10,k)]
                              -   2.0 * u[IDX(i,11,k)]
                            ) * idy_sqrd_by_180;
        
        // This is a centered sixth order stencil.   
        DyDyu[IDX(i,7,k)] = (      2.0 * u[IDX(i, 4,k)]
                              -   27.0 * u[IDX(i, 5,k)]
                              +  270.0 * u[IDX(i, 6,k)]
                              -  490.0 * u[IDX(i, 7,k)]
                              +  270.0 * u[IDX(i, 8,k)]
                              -   27.0 * u[IDX(i, 9,k)]
                              +    2.0 * u[IDX(i,10,k)]
                            ) * idy_sqrd_by_180;
        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a centered sixth order stencil.   
        DyDyu[IDX(i,je-4,k)] = (     2.0 * u[IDX(i,je-7,k)]
                                 -  27.0 * u[IDX(i,je-6,k)]
                                 + 270.0 * u[IDX(i,je-5,k)]
                                 - 490.0 * u[IDX(i,je-4,k)]
                                 + 270.0 * u[IDX(i,je-3,k)]
                                 -  27.0 * u[IDX(i,je-2,k)]
                                 +   2.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_180;

        // This is a shifted sixth order stencil.   
        DyDyu[IDX(i,je-3,k)] = ( -   2.0 * u[IDX(i,je-8,k)]
                                 +  16.0 * u[IDX(i,je-7,k)]
                                 -  54.0 * u[IDX(i,je-6,k)]
                                 +  85.0 * u[IDX(i,je-5,k)]
                                 + 130.0 * u[IDX(i,je-4,k)]
                                 - 378.0 * u[IDX(i,je-3,k)]
                                 + 214.0 * u[IDX(i,je-2,k)]
                                 -  11.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_180;

        // This is a (partially) shifted fourth order stencil.   
        DyDyu[IDX(i,je-2,k)] = (          u[IDX(i,je-6,k)]
                                 -  6.0 * u[IDX(i,je-5,k)]
                                 + 14.0 * u[IDX(i,je-4,k)]
                                 -  4.0 * u[IDX(i,je-3,k)]
                                 - 15.0 * u[IDX(i,je-2,k)]
                                 + 10.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12; 

        // XThis is a (totally) shifted fourth order stencil.   
        DyDyu[IDX(i,je-1,k)] = ( -  10.0 * u[IDX(i,je-6,k)]
                                 +  61.0 * u[IDX(i,je-5,k)]
                                 - 156.0 * u[IDX(i,je-4,k)]
                                 + 254.0 * u[IDX(i,je-3,k)]
                                 - 154.0 * u[IDX(i,je-2,k)]
                                 +  45.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8664_zz(double * const  DzDzu, const double * const  u,
                  const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_12  = idz_sqrd /  12.0;
  const double idz_sqrd_by_180 = idz_sqrd / 180.0;
  const double idz_sqrd_by_5040 = idz_sqrd / 5040.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);
        DzDzu[pp] = ( -     9.0 * u[pp-4*n]
                      +   128.0 * u[pp-3*n]
                      -  1008.0 * u[pp-2*n]
                      +  8064.0 * u[pp-  n]
                      - 14350.0 * u[pp    ]
                      +  8064.0 * u[pp+  n]
                      -  1008.0 * u[pp+2*n]
                      +   128.0 * u[pp+3*n]
                      -     9.0 * u[pp+4*n]
                    ) * idz_sqrd_by_5040;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted fourth order stencil.   
        DzDzu[IDX(i,j,4)] = (    45.0 * u[IDX(i,j,4)]
                              - 154.0 * u[IDX(i,j,5)]
                              + 214.0 * u[IDX(i,j,6)]
                              - 156.0 * u[IDX(i,j,7)]
                              +  61.0 * u[IDX(i,j,8)]
                              -  10.0 * u[IDX(i,j,9)]
                            ) * idz_sqrd_by_12;
        
        // This is a (partially) shifted fourth order stencil.   
        DzDzu[IDX(i,j,5)] = (   10.0 * u[IDX(i,j,4)]
                              - 15.0 * u[IDX(i,j,5)]
                              -  4.0 * u[IDX(i,j,6)]
                              + 14.0 * u[IDX(i,j,7)]
                              -  6.0 * u[IDX(i,j,8)]
                              +        u[IDX(i,j,9)]
                            ) * idz_sqrd_by_12;

        // This is a shifted sixth order stencil.   
        DzDzu[IDX(i,j,6)] = ( -  11.0 * u[IDX(i,j, 4)]
                              + 214.0 * u[IDX(i,j, 5)]
                              - 378.0 * u[IDX(i,j, 6)]
                              + 130.0 * u[IDX(i,j, 7)]
                              +  85.0 * u[IDX(i,j, 8)]
                              -  54.0 * u[IDX(i,j, 9)]
                              +  16.0 * u[IDX(i,j,10)]
                              -   2.0 * u[IDX(i,j,11)]
                            ) * idz_sqrd_by_180;
        
        // This is a centered sixth order stencil.   
        DzDzu[IDX(i,j,7)] = (      2.0 * u[IDX(i,j, 4)]
                              -   27.0 * u[IDX(i,j, 5)]
                              +  270.0 * u[IDX(i,j, 6)]
                              -  490.0 * u[IDX(i,j, 7)]
                              +  270.0 * u[IDX(i,j, 8)]
                              -   27.0 * u[IDX(i,j, 9)]
                              +    2.0 * u[IDX(i,j,10)]
                            ) * idz_sqrd_by_180;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered sixth order stencil.   
        DzDzu[IDX(i,j,ke-4)] = (     2.0 * u[IDX(i,j,ke-7)]
                                 -  27.0 * u[IDX(i,j,ke-6)]
                                 + 270.0 * u[IDX(i,j,ke-5)]
                                 - 490.0 * u[IDX(i,j,ke-4)]
                                 + 270.0 * u[IDX(i,j,ke-3)]
                                 -  27.0 * u[IDX(i,j,ke-2)]
                                 +   2.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_180;

        // This is a shifted sixth order stencil.   
        DzDzu[IDX(i,j,ke-3)] = ( -   2.0 * u[IDX(i,j,ke-8)]
                                 +  16.0 * u[IDX(i,j,ke-7)]
                                 -  54.0 * u[IDX(i,j,ke-6)]
                                 +  85.0 * u[IDX(i,j,ke-5)]
                                 + 130.0 * u[IDX(i,j,ke-4)]
                                 - 378.0 * u[IDX(i,j,ke-3)]
                                 + 214.0 * u[IDX(i,j,ke-2)]
                                 -  11.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_180;

        // This is a (partially) shifted fourth order stencil.   
        DzDzu[IDX(i,j,ke-2)] = (          u[IDX(i,j,ke-6)]
                                 -  6.0 * u[IDX(i,j,ke-5)]
                                 + 14.0 * u[IDX(i,j,ke-4)]
                                 -  4.0 * u[IDX(i,j,ke-3)]
                                 - 15.0 * u[IDX(i,j,ke-2)]
                                 + 10.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;
        
        // XThis is a (totally) shifted fourth order stencil.   
        DzDzu[IDX(i,j,ke-1)] = ( -  10.0 * u[IDX(i,j,ke-6)]
                                 +  61.0 * u[IDX(i,j,ke-5)]
                                 - 156.0 * u[IDX(i,j,ke-4)]
                                 + 214.0 * u[IDX(i,j,ke-3)]
                                 - 154.0 * u[IDX(i,j,ke-2)]
                                 +  45.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_12;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8664_xy( double * const DxDyu, double* const Dxu, 
                   const double * const u, const double dx, 
                   const double dy, const unsigned int *sz, unsigned bflag )
{
    deriv8664_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv8664_y( DxDyu, Dxu, dy, sz, bflag ) ; 
} 


void deriv8664_xz( double * const DxDzu, double* const Dxu, 
                   const double * const u, const double dx, 
                   const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv8664_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv8664_z( DxDzu, Dxu, dz, sz, bflag ) ; 
} 


void deriv8664_yz( double * const DyDzu, double* const Dyu, 
                   const double * const u, const double dy, 
                   const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv8664_y( Dyu  ,   u, dy, sz, bflag ) ; 
    deriv8664_z( DyDzu, Dyu, dz, sz, bflag ) ; 
} 


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8666_x(double * const  Dxu, const double * const  u,
                 const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx = 1.0 / dx;
  const double idx_by_2 = 0.5 * idx;
  const double idx_by_12 = idx / 12.0;
  const double idx_by_60 = idx / 60.0;
  const double idx_by_2520 = idx / 2520.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 0;
  const int kb = 0;
  const int ie = sz[0]-4;
  const int je = sz[1]-0;
  const int ke = sz[2]-0;
  
  const int n = 1;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);
        
        Dxu[pp] = (      9.0 * u[pp-4] 
                    -   96.0 * u[pp-3]
                    +  504.0 * u[pp-2]
                    - 2016.0 * u[pp-1]
                    + 2016.0 * u[pp+1]
                    -  504.0 * u[pp+2]
                    +   96.0 * u[pp+3]
                    -    9.0 * u[pp+4] ) * idx_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
       
        // This is a (totally) shifted sixth order stencil.   
        Dxu[IDX(4,j,k)] = ( - 147.0 * u[IDX( 4,j,k)]
                            + 360.0 * u[IDX( 5,j,k)]
                            - 450.0 * u[IDX( 6,j,k)]
                            + 400.0 * u[IDX( 7,j,k)]
                            - 225.0 * u[IDX( 8,j,k)]
                            +  72.0 * u[IDX( 9,j,k)]
                            -  10.0 * u[IDX(10,j,k)]
                          ) * idx_by_60;
        
        // This is a shifted sixth order stencil.   
        Dxu[IDX(5,j,k)] = ( -  10.0 * u[IDX( 4,j,k)]
                            -  77.0 * u[IDX( 5,j,k)]
                            + 150.0 * u[IDX( 6,j,k)]
                            - 100.0 * u[IDX( 7,j,k)]
                            +  50.0 * u[IDX( 8,j,k)]
                            -  15.0 * u[IDX( 9,j,k)]
                            +   2.0 * u[IDX(10,j,k)]
                          ) * idx_by_60;

        // This is a shifted sixth order stencil.   
        Dxu[IDX(6,j,k)] = (    2.0 * u[IDX( 4,j,k)]
                            - 24.0 * u[IDX( 5,j,k)]
                            - 35.0 * u[IDX( 6,j,k)]
                            + 80.0 * u[IDX( 7,j,k)]
                            - 30.0 * u[IDX( 8,j,k)]
                            +  8.0 * u[IDX( 9,j,k)]
                            -        u[IDX(10,j,k)]
                          ) * idx_by_60;

        // This is a centered sixth order stencil.   
        Dxu[IDX(7,j,k)] = ( -        u[IDX( 4,j,k)]
                            +  9.0 * u[IDX( 5,j,k)]
                            - 45.0 * u[IDX( 6,j,k)]
                            + 45.0 * u[IDX( 8,j,k)]
                            -  9.0 * u[IDX( 9,j,k)]
                            +        u[IDX(10,j,k)]
                          ) * idx_by_60;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered sixth order stencil.   
        Dxu[IDX(ie-4,j,k)] = ( -        u[IDX(ie-7,j,k)]
                               +  9.0 * u[IDX(ie-6,j,k)]
                               - 45.0 * u[IDX(ie-5,j,k)]
                               + 45.0 * u[IDX(ie-3,j,k)]
                               -  9.0 * u[IDX(ie-2,j,k)]
                               +        u[IDX(ie-1,j,k)]
                             ) * idx_by_60;

        // This is a shifted sixth order stencil.   
        Dxu[IDX(ie-3,j,k)] = (          u[IDX(ie-7,j,k)]
                               -  8.0 * u[IDX(ie-6,j,k)]
                               + 30.0 * u[IDX(ie-5,j,k)]
                               - 80.0 * u[IDX(ie-4,j,k)]
                               + 35.0 * u[IDX(ie-3,j,k)]
                               + 24.0 * u[IDX(ie-2,j,k)]
                               -  2.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_60;

        // This is a shifted sixth order stencil.   
        Dxu[IDX(ie-2,j,k)] = ( -   2.0 * u[IDX(ie-7,j,k)]
                               +  15.0 * u[IDX(ie-6,j,k)]
                               -  50.0 * u[IDX(ie-5,j,k)]
                               + 100.0 * u[IDX(ie-4,j,k)]
                               - 150.0 * u[IDX(ie-3,j,k)]
                               +  77.0 * u[IDX(ie-2,j,k)]
                               +  10.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_60;

        // This is a shifted sixth order stencil.   
        Dxu[IDX(ie-1,j,k)] = (    10.0 * u[IDX(ie-7,j,k)]
                               -  72.0 * u[IDX(ie-6,j,k)]
                               + 225.0 * u[IDX(ie-5,j,k)]
                               - 400.0 * u[IDX(ie-4,j,k)]
                               + 450.0 * u[IDX(ie-3,j,k)]
                               - 360.0 * u[IDX(ie-2,j,k)]
                               + 147.0 * u[IDX(ie-1,j,k)]
                             ) * idx_by_60;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < sz[2]-4; k++) {
      for (int j = jb; j < sz[1]-4; j++) {
        for (int i = ib; i < sz[0]-4; i++) {
          const int pp = IDX(i,j,k);
          if(isnan(Dxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif


}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8666_y(double * const  Dyu, const double * const  u,
                 const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy = 1.0 / dy;
  const double idy_by_2 = 0.5 * idy;
  const double idy_by_12 = idy / 12.0;
  const double idy_by_60 = idy / 60.0;
  const double idy_by_2520 = idy / 2520.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 0;
  const int ie = sz[0]-4;
  const int je = sz[1]-4;
  const int ke = sz[2]-0;
  
  const int n = nx;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);

        Dyu[pp] = (      9.0 * u[pp-4*n] 
                    -   96.0 * u[pp-3*n]
                    +  504.0 * u[pp-2*n]
                    - 2016.0 * u[pp-  n]
                    + 2016.0 * u[pp+  n]
                    -  504.0 * u[pp+2*n]
                    +   96.0 * u[pp+3*n]
                    -    9.0 * u[pp+4*n] ) * idy_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted sixth order stencil.   
        Dyu[IDX(i,4,k)] = ( - 147.0 * u[IDX(i, 4,k)]
                            + 360.0 * u[IDX(i, 5,k)]
                            - 450.0 * u[IDX(i, 6,k)]
                            + 400.0 * u[IDX(i, 7,k)]
                            - 225.0 * u[IDX(i, 8,k)]
                            +  72.0 * u[IDX(i, 9,k)]
                            -  10.0 * u[IDX(i,10,k)]
                          ) * idy_by_60;

        // This is a shifted sixth order stencil.   
         Dyu[IDX(i,5,k)] = ( - 10.0 * u[IDX(i, 4,k)]
                            -  77.0 * u[IDX(i, 5,k)]
                            + 150.0 * u[IDX(i, 6,k)]
                            - 100.0 * u[IDX(i, 7,k)]
                            +  50.0 * u[IDX(i, 8,k)]
                            -  15.0 * u[IDX(i, 9,k)]
                            +   2.0 * u[IDX(i,10,k)]
                          ) * idy_by_60;

        // This is a shifted sixth order stencil.   
        Dyu[IDX(i,6,k)] = (    2.0 * u[IDX(i, 4,k)]
                            - 24.0 * u[IDX(i, 5,k)]
                            - 35.0 * u[IDX(i, 6,k)]
                            + 80.0 * u[IDX(i, 7,k)]
                            - 30.0 * u[IDX(i, 8,k)]
                            +  8.0 * u[IDX(i, 9,k)]
                            -        u[IDX(i,10,k)]
                          ) * idy_by_60;

        // This is a centered sixth order stencil.   
        Dyu[IDX(i,7,k)] = ( -        u[IDX(i, 4,k)]
                            +  9.0 * u[IDX(i, 5,k)]
                            - 45.0 * u[IDX(i, 6,k)]
                            + 45.0 * u[IDX(i, 8,k)]
                            -  9.0 * u[IDX(i, 9,k)]
                            +        u[IDX(i,10,k)]
                          ) * idy_by_60;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a centered sixth order stencil.   
        Dyu[IDX(i,je-4,k)] = ( -        u[IDX(i,je-7,k)]
                               +  9.0 * u[IDX(i,je-6,k)]
                               - 45.0 * u[IDX(i,je-5,k)]
                               + 45.0 * u[IDX(i,je-3,k)]
                               -  9.0 * u[IDX(i,je-2,k)]
                               +        u[IDX(i,je-1,k)]
                             ) * idy_by_60;

        // This is a shifted sixth order stencil.   
        Dyu[IDX(i,je-3,k)] = (          u[IDX(i,je-7,k)]
                               -  8.0 * u[IDX(i,je-6,k)]
                               + 30.0 * u[IDX(i,je-5,k)]
                               - 80.0 * u[IDX(i,je-4,k)]
                               + 35.0 * u[IDX(i,je-3,k)]
                               + 24.0 * u[IDX(i,je-2,k)]
                               -  2.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_60;

        // This is a shifted sixth order stencil.   
        Dyu[IDX(i,je-2,k)] = ( -   2.0 * u[IDX(i,je-7,k)]
                               +  15.0 * u[IDX(i,je-6,k)]
                               -  50.0 * u[IDX(i,je-5,k)]
                               + 100.0 * u[IDX(i,je-4,k)]
                               - 150.0 * u[IDX(i,je-3,k)]
                               +  77.0 * u[IDX(i,je-2,k)]
                               +  10.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_60;
                             
        // This is a (totally) shifted sixth order stencil.   
        Dyu[IDX(i,je-1,k)] = (    10.0 * u[IDX(i,je-7,k)]
                               -  72.0 * u[IDX(i,je-6,k)]
                               + 225.0 * u[IDX(i,je-5,k)]
                               - 400.0 * u[IDX(i,je-4,k)]
                               + 450.0 * u[IDX(i,je-3,k)]
                               - 360.0 * u[IDX(i,je-2,k)]
                               + 147.0 * u[IDX(i,je-1,k)]
                             ) * idy_by_12;
      
      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < sz[2]-4; k++) {
      for (int j = jb; j < sz[1]-4; j++) {
        for (int i = ib; i < sz[0]-4; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8666_z(double * const  Dzu, const double * const  u,
                 const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz = 1.0 / dz;
  const double idz_by_2 = 0.5 * idz;
  const double idz_by_12 = idz / 12.0;
  const double idz_by_60 = idz / 60.0;
  const double idz_by_2520 = idz / 2520.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0]-4;
  const int je = sz[1]-4;
  const int ke = sz[2]-4;

  const int n = nx*ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);

        Dzu[pp] = (      9.0 * u[pp-4*n] 
                    -   96.0 * u[pp-3*n]
                    +  504.0 * u[pp-2*n]
                    - 2016.0 * u[pp-  n]
                    + 2016.0 * u[pp+  n]
                    -  504.0 * u[pp+2*n]
                    +   96.0 * u[pp+3*n]
                    -    9.0 * u[pp+4*n] ) * idz_by_2520;
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted sixth order stencil.   
        Dzu[IDX(i,j,4)] = ( - 147.0 * u[IDX(i,j, 4)]
                            + 360.0 * u[IDX(i,j, 5)]
                            - 450.0 * u[IDX(i,j, 6)]
                            + 400.0 * u[IDX(i,j, 7)]
                            - 225.0 * u[IDX(i,j, 8)]
                            +  72.0 * u[IDX(i,j, 9)]
                            -  10.0 * u[IDX(i,j,10)]
                          ) * idz_by_60;

        // This is a shifted sixth order stencil.   
        Dzu[IDX(i,j,5)] = ( -  10.0 * u[IDX(i,j, 4)]
                            -  77.0 * u[IDX(i,j, 5)]
                            + 150.0 * u[IDX(i,j, 6)]
                            - 100.0 * u[IDX(i,j, 7)]
                            +  50.0 * u[IDX(i,j, 8)]
                            -  15.0 * u[IDX(i,j, 9)]
                            +   2.0 * u[IDX(i,j,10)]
                          ) * idz_by_60;

        // This is a shifted sixth order stencil.   
        Dzu[IDX(i,j,6)] = (    2.0 * u[IDX(i,j, 4)]
                            - 24.0 * u[IDX(i,j, 5)]
                            - 35.0 * u[IDX(i,j, 6)]
                            + 80.0 * u[IDX(i,j, 7)]
                            - 30.0 * u[IDX(i,j, 8)]
                            +  8.0 * u[IDX(i,j, 9)]
                            -        u[IDX(i,j,10)]
                          ) * idz_by_60;

        // This is a centered sixth order stencil.   
        Dzu[IDX(i,j,7)] = ( -        u[IDX(i,j, 4)]
                            +  9.0 * u[IDX(i,j, 5)]
                            - 45.0 * u[IDX(i,j, 6)]
                            + 45.0 * u[IDX(i,j, 8)]
                            -  9.0 * u[IDX(i,j, 9)]
                            +        u[IDX(i,j,10)]
                          ) * idz_by_60;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a centered sixth order stencil.   
        Dzu[IDX(i,j,ke-4)] = ( -        u[IDX(i,j,ke-7)]
                               +  9.0 * u[IDX(i,j,ke-6)]
                               - 45.0 * u[IDX(i,j,ke-5)]
                               + 45.0 * u[IDX(i,j,ke-3)]
                               -  9.0 * u[IDX(i,j,ke-2)]
                               +        u[IDX(i,j,ke-1)]
                             ) * idz_by_60;

        // This is a shifted sixth order stencil.   
        Dzu[IDX(i,j,ke-3)] = (          u[IDX(i,j,ke-7)]
                               -  8.0 * u[IDX(i,j,ke-6)]
                               + 30.0 * u[IDX(i,j,ke-5)]
                               - 80.0 * u[IDX(i,j,ke-4)]
                               + 35.0 * u[IDX(i,j,ke-3)]
                               + 24.0 * u[IDX(i,j,ke-2)]
                               -  2.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_60;

        // This is a (partially) shifted sixth order stencil.   
        Dzu[IDX(i,j,ke-2)] = ( -   2.0 * u[IDX(i,j,ke-7)]
                               +  15.0 * u[IDX(i,j,ke-6)]
                               -  50.0 * u[IDX(i,j,ke-5)]
                               + 100.0 * u[IDX(i,j,ke-4)]
                               - 150.0 * u[IDX(i,j,ke-3)]
                               +  77.0 * u[IDX(i,j,ke-2)]
                               +  10.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_60;
                             
        // This is a (totally) shifted sixth order stencil.   
        Dzu[IDX(i,j,ke-1)] = (    10.0 * u[IDX(i,j,ke-7)]
                               -  72.0 * u[IDX(i,j,ke-6)]
                               + 225.0 * u[IDX(i,j,ke-5)]
                               - 400.0 * u[IDX(i,j,ke-4)]
                               + 450.0 * u[IDX(i,j,ke-3)]
                               - 360.0 * u[IDX(i,j,ke-2)]
                               + 147.0 * u[IDX(i,j,ke-1)]
                             ) * idz_by_60;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < sz[2]-4; k++) {
      for (int j = jb; j < sz[1]-4; j++) {
        for (int i = ib; i < sz[0]-4; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(Dzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8666_xx(double * const  DxDxu, const double * const  u,
                  const double dx, const unsigned int *sz, unsigned bflag)
{
  const double idx_sqrd = 1.0 / (dx*dx);
  const double idx_sqrd_by_12  = idx_sqrd / 12.0;
  const double idx_sqrd_by_180 = idx_sqrd / 180.0;
  const double idx_sqrd_by_5040 = idx_sqrd / 5040.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  for (int k = kb; k < ke; k++) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        const int pp = IDX(i,j,k);

        DxDxu[pp] = ( -     9.0 * u[pp-4]
                      +   128.0 * u[pp-3]
                      -  1008.0 * u[pp-2]
                      +  8064.0 * u[pp-1]
                      - 14350.0 * u[pp  ]
                      +  8064.0 * u[pp+1]
                      -  1008.0 * u[pp+2]
                      +   128.0 * u[pp+3]
                      -     9.0 * u[pp+4]
                    ) * idx_sqrd_by_5040;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_LEFT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a (totally) shifted sixth order stencil.   
        DxDxu[IDX(4,j,k)] = (    938.0 * u[IDX( 4,j,k)]
                              - 4014.0 * u[IDX( 5,j,k)]
                              + 7911.0 * u[IDX( 6,j,k)]
                              - 9490.0 * u[IDX( 7,j,k)]
                              + 7380.0 * u[IDX( 8,j,k)]
                              - 3618.0 * u[IDX( 9,j,k)]
                              + 1019.0 * u[IDX(10,j,k)]
                              -  126.0 * u[IDX(11,j,k)]
                            ) * idx_sqrd_by_180;
          
        // This is a (partially) shifted sixth order stencil.   
        DxDxu[IDX(5,j,k)] = (   126.0 * u[IDX( 4,j,k)]
                              -  70.0 * u[IDX( 5,j,k)]
                              - 486.0 * u[IDX( 6,j,k)]
                              + 855.0 * u[IDX( 7,j,k)]
                              - 670.0 * u[IDX( 8,j,k)]
                              + 324.0 * u[IDX( 9,j,k)]
                              -  90.0 * u[IDX(10,j,k)]
                              +  11.0 * u[IDX(11,j,k)]
                            ) * idx_sqrd_by_180;

        // This is a shifted sixth order stencil.   
        DxDxu[IDX(6,j,k)] = ( -  11.0 * u[IDX( 4,j,k)]
                              + 214.0 * u[IDX( 5,j,k)]
                              - 378.0 * u[IDX( 6,j,k)]
                              + 130.0 * u[IDX( 7,j,k)]
                              +  85.0 * u[IDX( 8,j,k)]
                              -  54.0 * u[IDX( 9,j,k)]
                              +  16.0 * u[IDX(10,j,k)]
                              -   2.0 * u[IDX(11,j,k)]
                            ) * idx_sqrd_by_180;

        // This is a centered sixth order stencil.   
        DxDxu[IDX(7,j,k)] = (      2.0 * u[IDX( 4,j,k)]
                              -   27.0 * u[IDX( 5,j,k)]
                              +  270.0 * u[IDX( 6,j,k)]
                              -  490.0 * u[IDX( 7,j,k)]
                              +  270.0 * u[IDX( 8,j,k)]
                              -   27.0 * u[IDX( 9,j,k)]
                              +    2.0 * u[IDX(10,j,k)]
                            ) * idx_sqrd_by_180;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_RIGHT)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {

        // This is a centered sixth order stencil.   
        DxDxu[IDX(ie-4,j,k)] = (     2.0 * u[IDX(ie-7,j,k)]
                                 -  27.0 * u[IDX(ie-6,j,k)]
                                 + 270.0 * u[IDX(ie-5,j,k)]
                                 - 490.0 * u[IDX(ie-4,j,k)]
                                 + 270.0 * u[IDX(ie-3,j,k)]
                                 -  27.0 * u[IDX(ie-2,j,k)]
                                 +   2.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_180;
        
        // This is a shifted sixth order stencil.   
        DxDxu[IDX(ie-3,j,k)] = ( -   2.0 * u[IDX(ie-8,j,k)]
                                 +  16.0 * u[IDX(ie-7,j,k)]
                                 -  54.0 * u[IDX(ie-6,j,k)]
                                 +  85.0 * u[IDX(ie-5,j,k)]
                                 + 130.0 * u[IDX(ie-4,j,k)]
                                 - 378.0 * u[IDX(ie-3,j,k)]
                                 + 214.0 * u[IDX(ie-2,j,k)]
                                 -  11.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_180;
 
        // This is a (partially) shifted sixth order stencil.   
        DxDxu[IDX(ie-2,j,k)] = (    11.0 * u[IDX(ie-8,j,k)]
                                 -  90.0 * u[IDX(ie-7,j,k)]
                                 + 324.0 * u[IDX(ie-6,j,k)]
                                 - 670.0 * u[IDX(ie-5,j,k)]
                                 + 855.0 * u[IDX(ie-4,j,k)]
                                 - 486.0 * u[IDX(ie-3,j,k)]
                                 -  70.0 * u[IDX(ie-2,j,k)]
                                 + 126.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_180;
         
        // XThis is a (totally) shifted sixth order stencil.   
        DxDxu[IDX(ie-1,j,k)] = ( -  126.0 * u[IDX(ie-8,j,k)]
                                 + 1019.0 * u[IDX(ie-7,j,k)]
                                 - 3618.0 * u[IDX(ie-6,j,k)]
                                 + 7380.0 * u[IDX(ie-5,j,k)]
                                 - 9490.0 * u[IDX(ie-4,j,k)]
                                 + 7911.0 * u[IDX(ie-3,j,k)]
                                 - 4014.0 * u[IDX(ie-2,j,k)]
                                 +  938.0 * u[IDX(ie-1,j,k)]
                               ) * idx_sqrd_by_180;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DxDxu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8666_yy(double * const  DyDyu, const double * const  u,
                  const double dy, const unsigned int *sz, unsigned bflag)
{
  const double idy_sqrd = 1.0/(dy*dy);
  const double idy_sqrd_by_12  = idy_sqrd /  12.0;
  const double idy_sqrd_by_180 = idy_sqrd / 180.0;
  const double idy_sqrd_by_5040 = idy_sqrd / 5040.0;
  
  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  for (int k = kb; k < ke; k++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int j = jb; j < je; j++) {
        const int pp = IDX(i,j,k);
 
        DyDyu[pp] = ( -     9.0 * u[pp-4*nx]
                      +   128.0 * u[pp-3*nx]
                      -  1008.0 * u[pp-2*nx]
                      +  8064.0 * u[pp-  nx]
                      - 14350.0 * u[pp     ]
                      +  8064.0 * u[pp+  nx]
                      -  1008.0 * u[pp+2*nx]
                      +   128.0 * u[pp+3*nx]
                      -     9.0 * u[pp+4*nx]
                    ) * idy_sqrd_by_5040;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_DOWN)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted sixth order stencil.   
        DyDyu[IDX(i,4,k)] = (    938.0 * u[IDX(i, 4,k)]
                              - 4014.0 * u[IDX(i, 5,k)]
                              + 7911.0 * u[IDX(i, 6,k)]
                              - 9490.0 * u[IDX(i, 7,k)]
                              + 7380.0 * u[IDX(i, 8,k)]
                              - 3618.0 * u[IDX(i, 9,k)]
                              + 1019.0 * u[IDX(i,10,k)]
                              -  126.0 * u[IDX(i,11,k)]
                            ) * idy_sqrd_by_180;
        
        // This is a (partially) shifted sixth order stencil.   
        DyDyu[IDX(i,5,k)] = (   126.0 * u[IDX(i, 4,k)]
                              -  70.0 * u[IDX(i, 5,k)]
                              - 486.0 * u[IDX(i, 6,k)]
                              + 855.0 * u[IDX(i, 7,k)]
                              - 670.0 * u[IDX(i, 8,k)]
                              + 324.0 * u[IDX(i, 9,k)]
                              -  90.0 * u[IDX(i,10,k)]
                              +  11.0 * u[IDX(i,11,k)]
                            ) * idy_sqrd_by_180;

        // This is a shifted sixth order stencil.   
        DyDyu[IDX(i,6,k)] = ( -  11.0 * u[IDX(i, 4,k)]
                              + 214.0 * u[IDX(i, 5,k)]
                              - 378.0 * u[IDX(i, 6,k)]
                              + 130.0 * u[IDX(i, 7,k)]
                              +  85.0 * u[IDX(i, 8,k)]
                              -  54.0 * u[IDX(i, 9,k)]
                              +  16.0 * u[IDX(i,10,k)]
                              -   2.0 * u[IDX(i,11,k)]
                            ) * idy_sqrd_by_180;
        
        // This is a centered sixth order stencil.   
        DyDyu[IDX(i,7,k)] = (      2.0 * u[IDX(i, 4,k)]
                              -   27.0 * u[IDX(i, 5,k)]
                              +  270.0 * u[IDX(i, 6,k)]
                              -  490.0 * u[IDX(i, 7,k)]
                              +  270.0 * u[IDX(i, 8,k)]
                              -   27.0 * u[IDX(i, 9,k)]
                              +    2.0 * u[IDX(i,10,k)]
                            ) * idy_sqrd_by_180;
        
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_UP)) {
    for (int k = kb; k < ke; k++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a centered sixth order stencil.   
        DyDyu[IDX(i,je-4,k)] = (     2.0 * u[IDX(i,je-7,k)]
                                 -  27.0 * u[IDX(i,je-6,k)]
                                 + 270.0 * u[IDX(i,je-5,k)]
                                 - 490.0 * u[IDX(i,je-4,k)]
                                 + 270.0 * u[IDX(i,je-3,k)]
                                 -  27.0 * u[IDX(i,je-2,k)]
                                 +   2.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_180;

        // This is a shifted sixth order stencil.   
        DyDyu[IDX(i,je-3,k)] = ( -   2.0 * u[IDX(i,je-8,k)]
                                 +  16.0 * u[IDX(i,je-7,k)]
                                 -  54.0 * u[IDX(i,je-6,k)]
                                 +  85.0 * u[IDX(i,je-5,k)]
                                 + 130.0 * u[IDX(i,je-4,k)]
                                 - 378.0 * u[IDX(i,je-3,k)]
                                 + 214.0 * u[IDX(i,je-2,k)]
                                 -  11.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_180;

        // This is a (partially) shifted sixth order stencil.   
        DyDyu[IDX(i,je-2,k)] = (    11.0 * u[IDX(i,je-8,k)]
                                 -  90.0 * u[IDX(i,je-7,k)]
                                 + 324.0 * u[IDX(i,je-6,k)]
                                 - 670.0 * u[IDX(i,je-5,k)]
                                 + 855.0 * u[IDX(i,je-4,k)]
                                 - 486.0 * u[IDX(i,je-3,k)]
                                 -  70.0 * u[IDX(i,je-2,k)]
                                 + 126.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_180; 

        // XThis is a (totally) shifted sixth order stencil.   
        DyDyu[IDX(i,je-1,k)] = ( -  126.0 * u[IDX(i,je-8,k)]
                                 + 1019.0 * u[IDX(i,je-7,k)]
                                 - 3618.0 * u[IDX(i,je-6,k)]
                                 + 7380.0 * u[IDX(i,je-5,k)]
                                 - 9490.0 * u[IDX(i,je-4,k)]
                                 + 7911.0 * u[IDX(i,je-3,k)]
                                 - 4014.0 * u[IDX(i,je-2,k)]
                                 +  938.0 * u[IDX(i,je-1,k)]
                               ) * idy_sqrd_by_180;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DyDyu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}



/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8666_zz(double * const  DzDzu, const double * const  u,
                  const double dz, const unsigned int *sz, unsigned bflag)
{
  const double idz_sqrd = 1.0/(dz*dz);
  const double idz_sqrd_by_12  = idz_sqrd /  12.0;
  const double idz_sqrd_by_180 = idz_sqrd / 180.0;
  const double idz_sqrd_by_5040 = idz_sqrd / 5040.0;

  const int nx = sz[0];
  const int ny = sz[1];
  const int nz = sz[2];
  const int ib = 4;
  const int jb = 4;
  const int kb = 4;
  const int ie = sz[0] - 4;
  const int je = sz[1] - 4;
  const int ke = sz[2] - 4;

  const int n = nx * ny;

  for (int j = jb; j < je; j++) {
    for (int i = ib; i < ie; i++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int k = kb; k < ke; k++) {
        const int pp = IDX(i,j,k);

        DzDzu[pp] = ( -     9.0 * u[pp-4*n]
                      +   128.0 * u[pp-3*n]
                      -  1008.0 * u[pp-2*n]
                      +  8064.0 * u[pp-  n]
                      - 14350.0 * u[pp    ]
                      +  8064.0 * u[pp+  n]
                      -  1008.0 * u[pp+2*n]
                      +   128.0 * u[pp+3*n]
                      -     9.0 * u[pp+4*n]
                    ) * idz_sqrd_by_5040;
      
      }
    }
  }

  if (bflag & (1u<<OCT_DIR_BACK)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {

        // This is a (totally) shifted sixth order stencil.   
        DzDzu[IDX(i,j,4)] = (    938.0 * u[IDX(i,j, 4)]
                              - 4014.0 * u[IDX(i,j, 5)]
                              + 7911.0 * u[IDX(i,j, 6)]
                              - 9490.0 * u[IDX(i,j, 7)]
                              + 7380.0 * u[IDX(i,j, 8)]
                              - 3618.0 * u[IDX(i,j, 9)]
                              - 1019.0 * u[IDX(i,j,10)]
                              -  126.0 * u[IDX(i,j,11)]
                            ) * idz_sqrd_by_180;
        
        // This is a (partially) shifted sixth order stencil.   
        DzDzu[IDX(i,j,5)] = (   126.0 * u[IDX(i,j, 4)]
                              -  70.0 * u[IDX(i,j, 5)]
                              - 486.0 * u[IDX(i,j, 6)]
                              + 855.0 * u[IDX(i,j, 7)]
                              - 670.0 * u[IDX(i,j, 8)]
                              + 324.0 * u[IDX(i,j, 9)]
                              -  90.0 * u[IDX(i,j,10)]
                              +  11.0 * u[IDX(i,j,11)]
                            ) * idz_sqrd_by_180;

        // This is a shifted sixth order stencil.   
        DzDzu[IDX(i,j,6)] = ( -  11.0 * u[IDX(i,j, 4)]
                              + 214.0 * u[IDX(i,j, 5)]
                              - 378.0 * u[IDX(i,j, 6)]
                              + 130.0 * u[IDX(i,j, 7)]
                              +  85.0 * u[IDX(i,j, 8)]
                              -  54.0 * u[IDX(i,j, 9)]
                              +  16.0 * u[IDX(i,j,10)]
                              -   2.0 * u[IDX(i,j,11)]
                            ) * idz_sqrd_by_180;
        
        // This is a centered sixth order stencil.   
        DzDzu[IDX(i,j,7)] = (      2.0 * u[IDX(i,j, 4)]
                              -   27.0 * u[IDX(i,j, 5)]
                              +  270.0 * u[IDX(i,j, 6)]
                              -  490.0 * u[IDX(i,j, 7)]
                              +  270.0 * u[IDX(i,j, 8)]
                              -   27.0 * u[IDX(i,j, 9)]
                              +    2.0 * u[IDX(i,j,10)]
                            ) * idz_sqrd_by_180;

      }
    }
  }

  if (bflag & (1u<<OCT_DIR_FRONT)) {
    for (int j = jb; j < je; j++) {
      #ifdef DERIV_ENABLE_AVX
        #ifdef __INTEL_COMPILER
          #pragma vector vectorlength(__DERIV_AVX_SIMD_LEN__) vecremainder
          #pragma ivdep
        #endif
      #endif
      for (int i = ib; i < ie; i++) {
        
        // This is a centered sixth order stencil.   
        DzDzu[IDX(i,j,ke-4)] = (     2.0 * u[IDX(i,j,ke-7)]
                                 -  27.0 * u[IDX(i,j,ke-6)]
                                 + 270.0 * u[IDX(i,j,ke-5)]
                                 - 490.0 * u[IDX(i,j,ke-4)]
                                 + 270.0 * u[IDX(i,j,ke-3)]
                                 -  27.0 * u[IDX(i,j,ke-2)]
                                 +   2.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_180;

        // This is a shifted sixth order stencil.   
        DzDzu[IDX(i,j,ke-3)] = ( -   2.0 * u[IDX(i,j,ke-8)]
                                 +  16.0 * u[IDX(i,j,ke-7)]
                                 -  54.0 * u[IDX(i,j,ke-6)]
                                 +  85.0 * u[IDX(i,j,ke-5)]
                                 + 130.0 * u[IDX(i,j,ke-4)]
                                 - 378.0 * u[IDX(i,j,ke-3)]
                                 + 214.0 * u[IDX(i,j,ke-2)]
                                 -  11.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_180;

        // This is a (partially) shifted sixth order stencil.   
        DzDzu[IDX(i,j,ke-2)] = (    11.0 * u[IDX(i,j,ke-8)]
                                 -  90.0 * u[IDX(i,j,ke-7)]
                                 - 324.0 * u[IDX(i,j,ke-6)]
                                 - 670.0 * u[IDX(i,j,ke-5)]
                                 + 855.0 * u[IDX(i,j,ke-4)]
                                 - 486.0 * u[IDX(i,j,ke-3)]
                                 -  70.0 * u[IDX(i,j,ke-2)]
                                 + 126.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_180;
        
        // XThis is a (totally) shifted sixth order stencil.   
        DzDzu[IDX(i,j,ke-1)] = ( -  126.0 * u[IDX(i,j,ke-8)]
                                 + 1019.0 * u[IDX(i,j,ke-7)]
                                 - 3618.0 * u[IDX(i,j,ke-6)]
                                 + 7380.0 * u[IDX(i,j,ke-5)]
                                 - 9490.0 * u[IDX(i,j,ke-4)]
                                 + 7911.0 * u[IDX(i,j,ke-3)]
                                 - 4014.0 * u[IDX(i,j,ke-2)]
                                 +  938.0 * u[IDX(i,j,ke-1)]
                               ) * idz_sqrd_by_180;

      }
    }
  }


  #ifdef DEBUG_DERIVS_COMP
  #pragma message("DEBUG_DERIVS_COMP: ON")
    for (int k = kb; k < ke; k++) {
      for (int j = jb; j < je; j++) {
        for (int i = ib; i < ie; i++) {
          const int pp = IDX(i,j,k);
          if(std::isnan(DzDzu[pp])) std::cout<<"NAN detected function "<<__func__<<" file: "<<__FILE__<<" line: "<<__LINE__<<std::endl;
        }
      }
    }
  #endif

}


/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void deriv8666_xy( double * const DxDyu, double* const Dxu, 
                   const double * const u, const double dx, 
                   const double dy, const unsigned int *sz, unsigned bflag )
{
    deriv8666_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv8666_y( DxDyu, Dxu, dy, sz, bflag ) ; 
} 


void deriv8666_xz( double * const DxDzu, double* const Dxu, 
                   const double * const u, const double dx, 
                   const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv8666_x( Dxu  ,   u, dx, sz, bflag ) ; 
    deriv8666_z( DxDzu, Dxu, dz, sz, bflag ) ; 
} 


void deriv8666_yz( double * const DyDzu, double* const Dyu, 
                   const double * const u, const double dy, 
                   const double dz, const unsigned int *sz, unsigned bflag )
{
    deriv8666_y( Dyu  ,   u, dy, sz, bflag ) ; 
    deriv8666_z( DyDzu, Dyu, dz, sz, bflag ) ; 
} 



