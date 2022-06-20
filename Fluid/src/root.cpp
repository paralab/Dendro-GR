#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "mpi.h"
#include "fluid.h"
#include "root.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#define LTRACE 0


// rbrac {{{
/*-------------------------------------------------------------------------
 *
 *  Routine from numerical recipies.  Attempt to bracket the roots of
 *  a function.  The function is in the routine func(double *par)
 *
 *-------------------------------------------------------------------------*/
int froot::rbrac(int (*func)(double *, double , double *), double *x1, double *x2,
           double *par, int trace)
{

  int    rc, j;
  double f1, f2;
  const double FACTOR = 1.5;
  const double MAXIT  = 50;

  if (*x1 == *x2) {
    if (trace) printf("RBRAC.  Initially x1==x2.  x1=%g, x2=%g\n",*x1,*x2);
    return -1;
  }

  rc = (*func)(&f1, *x1, par);
  if (rc < 0) {
    if (trace) printf("RBRAC: x1 out of range, x1=%g, f1=%g\n",*x1,f1);
    return -1;
  }
  rc = (*func)(&f2, *x2, par);
  if (rc < 0) {
    if (trace) printf("RBRAC: x2 out of range, x2=%g, f2=%g\n",*x2,f2);
    return -1;
  }
  for (j = 0; j < MAXIT; j++) {
    if (trace) {
      printf("RBRAC: Loop x1=%g, f1=%g, x2=%g, f2=%g\n",*x1,f1,*x2,f2);
    }
    if (f1*f2 < 0.0) {
      return 1;
    }
    if (fabs(f1) < fabs(f2)) {
      *x1 += FACTOR*(*x1 - *x2);
      rc = (*func)(&f1, *x1, par);
      if (rc < 0) {
        if (trace) printf("RBRAC (loop): x1 out of range, x1=%g, f1=%g\n",
             *x1,f1);
        return -1;
      }
    }
    else {
      *x2 += FACTOR*(*x2 - *x1);
      rc = (*func)(&f2, *x2, par);
      if (rc < 0) {
        if (trace) printf("RBRAC (loop): x2 out of range, x2=%g, f2=%g\n",
             *x2,f2);
        return -1;
      }
    }
  }

  if (trace) {
    printf("RBRAC:  too many iterations\n");
  }
  return -2;
}
// }}}

// bisec {{{
int froot::bisec(int (*func)(double *, double , double *),
          double *rtb, double x1, double x2,
          double tol, double *fpar, int trace)
{
  int rc, j;
  double dx, f, fm, xm;
  const int MAXIT = 100;

  rc = (*func)(&f, x1, fpar);
  if (rc < 0) {
    if (trace) {
      printf("BISEC: x1 out of range, x1 = %f\n",x1);
    }
    return -1;
  }
  rc = (*func)(&fm, x2, fpar);
  if (rc < 0) {
    if (trace) {
      printf("BISEC: x2 out of range, x2 = %f\n",x2);
    }
    return -1;
  }

  if (f*fm >= 0.0) {
    printf("BISEC:  Error, root has not been bracketed\n");
    printf("        x1=%f, f1=%f, x2=%f, f2=%f\n",x1,f,x2,fm);
    return -1;
  }

  if (f < 0.0) {
    *rtb = x1;
    dx = x2 - x1;
  }
  else {
    *rtb = x2;
    dx = x1 - x2;
  }
  for (j = 1; j < MAXIT; j++) {
    dx *= 0.5;
    xm = *rtb + dx;
    rc = (*func)(&fm, xm, fpar);
    if (rc < 0) {
      if (trace) {
        printf("BISEC (loop): xm out of range, xm = %f\n",xm);
      }
      return -1;
    }
    if (fm <= 0.0) *rtb = xm;
    if (fabs(dx) < tol || fm == 0.0) return 1;
  }

  if (trace) {
    printf("BISEC:  did not find root\n");
  }
  return -1;

}
// }}}

// rtsafe {{{
/*-------------------------------------------------------------------------
 *
 *  func - returns f(x0), df/dx(x0).  The arguments are
 *         func(f, df, x0, fpar)
 *
 *-------------------------------------------------------------------------*/
int froot::rtsafe(int (*func)(double *, double *, double , double *),
              double *rts, double x1, double x2, double tol,
              double *fpar,int trace)
{

  int j, rc;
  double df, dx, dxold, f, fh, fl;
  double tmp, xh, xl;

  const double MAXIT  = 50;

  rc = (*func)(&fl, &df, x1, fpar);
  if (rc != 1) return -1;
  rc = (*func)(&fh, &df, x2, fpar);
  if (rc != 1) return -1;

  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
    if ( trace ) {
      printf("RTSAFE: Error.  Root must be bracketed.\n");
      printf("x1 = %g, f1 = %g\n",x1,fl);
      printf("x2 = %g, f2 = %g\n",x2,fh);
    }
    return -1;
  }

  if (fl == 0.0) {
    *rts = x1;
    return 1;
  }
  if (fh == 0.0) {
    *rts = x2;
    return 1;
  }

  if (fl < 0.0) {
    xl = x1;
    xh = x2;
  }
  else {
    xl = x2;
    xh = x1;
  }

  *rts = 0.5*(x1 + x2);
  dxold = fabs(x2 - x1);
  dx = dxold;
  rc = (*func)(&f, &df, *rts, fpar);
  if (trace) {
    printf("RTSAFE: rts=%g, f=%g, df=%g\n",*rts, f, df);
  }

  if (rc != 1) return -1;

  for (j=0; j < MAXIT; j++) {
    if ((((*rts-xh)*df-f)*((*rts-xl)*df-f) > 0.0) ||
          (fabs(2.0*f) > fabs(dxold*df))) {
      /* Bisect if Newton solve would give root out of range, or if the
       * solver isn't converging quickly.
       */
      dxold = dx;
      dx = 0.5*(xh-xl);
      *rts = xl + dx;
      if (xl == *rts) return 1;
    }
    else {
      /* Newton solve */
      dxold = dx;
      dx = f/df;
      tmp = *rts;
      *rts -= dx;
      if (tmp == *rts) return 1;
    }
    if (fabs(dx) < tol) return 1;

    rc = (*func)(&f, &df, *rts, fpar);
    if (trace) {
      printf("RTSAFE: rts=%g, f=%g, df=%g\n",*rts, f, df);
    }
    if (rc != 1) return -1;

    if (f < 0.0)
      xl = *rts;
    else
      xh = *rts;
  }

  if ( trace ) {
    printf("RTSAFE:  No solution.\n");
  }
  return -1;

}
// }}}

// zbrent {{{
/*-------------------------------------------------------------------------
 *
 *  func - returns f(x0).  The arguments are func(f, x0, fpar)
 *
 *-------------------------------------------------------------------------*/
int froot::zbrent(int(*func)(double *, double , double *),
              double *rts, double x1, double x2, double tol, double *fpar)
{
  const int ltrace  = 0;
  const int ltrace2 = 1;
  const int ITMAX = 100;
  const double EPS = 1.0e-15;
  int iter;
  int rc;

  double a=x1,b=x2,c=x2,d,e,min1,min2;
  double fa, fb, fc;
  double p, q, r, s, tol1, xm;

  if (ltrace2) printf(">>> begin zbrent\n");

  rc = (*func)(&fa, a, fpar);
  rc = (*func)(&fb, b, fpar);

  if (ltrace) printf(">>> rc %d\n",rc);
  // initialize d,e to avoid a warning
  d = -999.e9;
  e = -999.e9;

#if 0
  if ((fa * fb) > 0.0) {
    if (fa < 0.0 && fb < 0.0) {
      for (int k = 0; k < 30; k++) {
        b *= 1.1;
        rc = (*func)(&fb, b, fpar);
        if (fb < 0.0) {
          a = b;
          fa = fb;
        }
        else {
          break;
        }
      }
      if (fb < 0.0) {
        printf("ZBRENT >>> failed to reset bounds. x2=%21.15e, b=%21.15e\n",x2,b);
      }
      else {
        if (ltrace2) printf("ZBRENT >>> bounds sucessfully reset. x2=%21.15e, b=%21.15e\n",x2,b);
      }
    }
    else {
      for (int k = 0; k < 30; k++) {
        a *= 0.9;
        rc = (*func)(&fa, a, fpar);
        if (fa > 0.0) {
          b = a;
          fb = fa;
        }
        else {
          break;
        }
      }
      if (fa > 0.0) {
        printf("ZBRENT >>> failed to reset bounds. x1=%21.15e, a=%21.15e\n",x1,a);
      }
      else {
        printf("ZBRENT >>> bounds sucessfully reset. x1=%21.15e, a=%21.15e\n",x1,a);
      }

    }
  }
#endif

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    printf("Root must be bracketed in zbrent.  fa=%g, fb=%g\n",fa,fb);
    return -1;
  }

  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if (LTRACE) {
      printf("a=%g, b=%g, c=%g, fa=%g, fb=%g, fc=%g\n",a,b,c,fa,fb,fc);
    }
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) {
      if (fb >= 0.0) {
        printf("sending root b. fb=%g, fc=%g, b=%21.15e, c=%21.15e\n",fb,fc,b,c);
        *rts = b;
      }
      else if (fc > 0.0) {
        printf("sending root c. fb=%g, fc=%g, b=%21.15e, c=%21.15e\n",fb,fc,b,c);
        *rts = c;
      }
      else if (fa > 0.0) {
        printf("sending root a. fb=%g, fa=%g, b=%21.15e, a=%21.15e\n",fb,fa,b,a);
        *rts = a;
      }
      else {
        *rts = b;
        printf("ZBRENT. Failed to distinguish the proper root. Should not happen.\n");
        exit(2);
      }
      //*rts = b;
      return 1;
    }
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
        p=2.0*xm*s;
        q=1.0-s;
      } else {
        q=fa/fc;
        r=fb/fc;
        p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
        e=d;
        d=p/q;
      } else {
        d=xm;
        e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    rc=(*func)(&fb, b, fpar);
  }
  printf("ZBRENT >>> Maximum number of iterations exceeded in zbrent\n");
  printf("           a = %21.15e, fa = %21.15e\n",a,fa);
  printf("           b = %21.15e, fb = %21.15e\n",b,fb);
  printf("           c = %21.15e, fc = %21.15e\n",c,fc);
  return -1;

}
// }}}

// check_finite_1d {{{
void froot::check_finite_1d(double *f, int *nx, int *rc)
{
  int i;
  int n = *nx;

  *rc = 1;
  for (i = 0; i < n; i++) {
    if (std::isfinite(f[i]) == 0) {
      printf("### check_finite_1d: error.  i=%d, f=%g\n",i,f[i]);
      *rc = 0;
    }
  }

}
// }}}

// check_finite_1d_ {{{
/*void froot::check_finite_1d_(double *f, int *nx, int *rc)
{
  int i;
  int n = *nx;

  *rc = 1;
  for (i = 0; i < n; i++) {
    if (isfinite(f[i]) == 0) {
      printf("### check_finite_1d: error.  i=%d, f=%g\n",i,f[i]);
      *rc = 0;
    }
  }

}*/
// }}}

// quiet_check_finite_1d {{{
void froot::quiet_check_finite_1d(double *f, int *nx, int *rc)
{
  int i;
  int n = *nx;

  *rc = 1;
  for (i = 0; i < n; i++) {
    if (std::isfinite(f[i]) == 0) {
      *rc = 0;
    }
  }

}
// }}}
