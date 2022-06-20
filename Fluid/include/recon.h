#ifndef RECON_H
#define RECON_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "fluid.h"
#include "fluidUtils.h"

#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

#define MINMOD4(w,x,y,z) \
    (0.125*( copysign(1.0,w)+copysign(1.0,x))*fabs((copysign(1.0,w)+copysign(1.0,y)) * (copysign(1.0,w)+copysign(1.0,z))) * fmin(fabs(w),fmin(fabs(x),fmin(fabs(y),fabs(z)))))

namespace recon{
  namespace none{
    void reconstruct(const int n, const double* const RESTRICT v,
                     double* const RESTRICT vl, double* const RESTRICT vr);
  }
  
  namespace minmod{
    void reconstruct(const int n, const double* const RESTRICT v,
                     double* const RESTRICT vl, double* const RESTRICT vr);

    void reconstructPt(const int n, const int i, const double * const RESTRICT v,
                       double* const RESTRICT vl, double* const RESTRICT vr);
  }
  
  namespace weno3{
    void reconstruct(const int n, const double* const RESTRICT v,
                     double* const RESTRICT vl, double* const RESTRICT vr);
  }

  namespace weno5{
    void reconstruct(const int n, const double* const RESTRICT v,
                     double* const RESTRICT vl, double* const RESTRICT vr);
  }

  namespace mp5{
    void reconstruct(const int n, const double* const RESTRICT v,
                     double* const RESTRICT vl, double* const RESTRICT vr);
    static inline double MP5(const double am2, const double am1, const double a,
                             const double ap1, const double ap2, const double anorm,
                             const double mp5_eps, const double mp5_alpha){
      const double vl = (2.0*am2 - 13.0*am1 + 47.0*a + 27.0*ap1 - 3.0*ap2) / 60.0;
      const double x = ap1-a;
      const double y = mp5_alpha*(a-am1);
      const double vmp = a + (0.5*(copysign(1.0,x) + copysign(1.0,y)) * fmin(fabs(x),fabs(y)));
      if((vl-a)*(vl-vmp) <= mp5_eps*anorm){
        return vl;
      }
      else{
        const double djm1 = am2 - 2.0*am1 + a;
        const double dj   = am1 - 2.0*a + ap1;
        const double djp1 = a - 2.0*ap1 + ap2;
        const double dm4jph = MINMOD4(4.0*dj-djp1, 4.0*djp1-dj, dj, djp1);
        const double dm4jmh = MINMOD4(4.0*dj-djm1, 4.0*djm1-dj, dj, djm1);
        const double vul = a + mp5_alpha*(a-am1);
        const double vav = 0.5*(a+ap1);
        const double vmd = vav - 0.5*dm4jph;
        const double vlc = a + 0.5*(a-am1) + 4.0/3.0*dm4jmh;
        const double vmin = fmax(fmin(a,fmin(ap1,vmd)), fmin(a, fmin(vul,vlc)));
        const double vmax = fmin(fmax(a,fmax(ap1,vmd)), fmax(a, fmax(vul,vlc)));
        return vl + (0.5*(copysign(1.0,vmin-vl) + copysign(1.0,vmax-vl)) * fmin(fabs(vmin-vl),fabs(vmax-vl)));
      }
    }
  }

  void reconstructvars(double *vl[], double *vr[], double *v1[], double *F1[],
                       const int dir, const double pos_in[], const double dx,
                       const unsigned int *sz);

  void recon_help(double *vl[], double *vr[], double *v1[], double *wv[],
                   const int dir, const double pos_in[], const double dx,
                   const unsigned int *sz);

  int cal_wv(double *wv[], double *v[], const int dir, const double pos_in[],
             const double dx, const unsigned int *sz);
}

#endif
