#include "trumpet.h"

namespace trumpet_data {

const unsigned int KMAX = 5000;
double dxsav;
static double *xp;
static double *yp[2];

void trumpetData(const double xx1, const double yy1, const double zz1,
                 double *var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    double eps      = 1.e-10;  // tolerance on the integrator
    double h1       = -1.e-5;  // starting step off the critical point
    double hmin     = 1.e-25;  // min allowed stepsize for adaptive integrator

    double bh_mass  = 7.5;
    double alpha_c  = 0.16227766016837933200;
    double C_sq     = 1.5543095902183302040;
    double bigR_c   = 1.5405694150420948330 * bh_mass;
    double r_c      = 0.30405997036 * bh_mass;

    static bool firstcall = true;

    const int neq         = 2;
    static double *alpha0;
    static double *bigR0;
    static double *r_iso;

    static int np;

    const double third = 1.0 / 3.0;

    if (firstcall == true) {
        // solve ODE system
        std::cout << "trumpetData:  first call. Solve the ODE" << std::endl;

        xp          = new double[KMAX];
        yp[0]       = new double[KMAX];
        yp[1]       = new double[KMAX];
        alpha0      = new double[KMAX];
        bigR0       = new double[KMAX];
        r_iso       = new double[KMAX];

        // integrate inward from r_C to r=0.
        double rmin = 0.0;     // min value for inward r integration
        double rmax = 1000.0;  // max value for outward r integration
        double rstart, ystart[2], dydr[2];
        bndcnd(h1, rstart, ystart, dydr);
        double rend = rmin + fabs(h1);
        int nok, nbad;
        int kount;
        odeint(ystart, neq, rstart, rend, eps, h1, hmin, &nok, &nbad, derivs,
               rkqs, kount);

        int kountin = kount;
        for (unsigned int i = 0; i < kountin; i++) {
            r_iso[i]  = xp[kountin - 1 - i];
            alpha0[i] = yp[0][kountin - 1 - i];
            bigR0[i]  = yp[1][kountin - 1 - i];
            // std::cout<<"<in> xp = "<<xp[i]<<", i="<<i<<", alpha0 =
            // "<<yp[0][i]<<", bigR0 = "<<yp[1][i]<<std::endl;
        }

        // integrate outwards from r_c to large r
        h1   = fabs(h1);
        rend = rmax;

        std::cout << "integrating outwards ... " << std::endl;
        bndcnd(h1, rstart, ystart, dydr);
        odeint(ystart, neq, rstart, rend, eps, h1, hmin, &nok, &nbad, derivs,
               rkqs, kount);

        int kountout = kount;
        for (unsigned int i = 0; i < kountout; i++) {
            r_iso[i + kountin]  = xp[i];
            alpha0[i + kountin] = yp[0][i];
            bigR0[i + kountin]  = yp[1][i];
            // std::cout<<"<out> xp = "<<xp[i]<<", alpha0 = "<<yp[0][i]<<",
            // bigR0 = "<<yp[1][i]<<std::endl;
        }
        np        = kountin + kountout;

        firstcall = false;
    }

    int nrby_indx = np / 2;
    double ax_eps = 1.0e-5;

    double tenh1  = 10.0 * fabs(h1);
    double rbar   = sqrt(xx * xx + yy * yy + zz * zz);
    if (fabs(xx) < tenh1 && fabs(yy) < tenh1 && fabs(zz) < tenh1) {
        rbar = tenh1;
    }
    double alpha = interpolation4(r_iso, alpha0, np, rbar, &nrby_indx);
    double bigR  = interpolation4(r_iso, bigR0, np, rbar, &nrby_indx);

    if (fabs(xx) < ax_eps && fabs(yy) < ax_eps && fabs(zz) < ax_eps) {
        rbar = sqrt(xx * xx + yy * yy + zz * zz + ax_eps * ax_eps);
    }

    double f0            = 1.0 - 2.0 * bh_mass / bigR;

    double Rsq_dalpha_dR = 4.0 * (alpha * alpha - f0 - 0.5 * bh_mass / bigR) /
                           (alpha * alpha - 2.0 * alpha - f0) * bigR;

    double tmp_sqrt = sqrt(alpha * alpha - f0);

    double tmp_chi  = (rbar * rbar) / (bigR * bigR);

    double tmp_trK = 2.0 * tmp_sqrt / bigR + (alpha * Rsq_dalpha_dR - bh_mass) /
                                                 (tmp_sqrt * bigR * bigR);

    double tmp_beta = tmp_sqrt / bigR;

    double tmp_Atilde =
        ((alpha * Rsq_dalpha_dR - bh_mass) / tmp_sqrt - bigR * tmp_sqrt) /
        bigR / bigR;

    var[bssn::VAR::U_ALPHA] = alpha;
    var[bssn::VAR::U_ALPHA] =
        std::max(var[bssn::VAR::U_ALPHA], bssn::CHI_FLOOR);

    var[bssn::VAR::U_CHI] = tmp_chi;

    if (var[bssn::VAR::U_CHI] < bssn::CHI_FLOOR)
        var[bssn::VAR::U_CHI] = bssn::CHI_FLOOR;

    var[bssn::VAR::U_K]      = tmp_trK;

    var[bssn::VAR::U_BETA0]  = xx * tmp_beta;
    var[bssn::VAR::U_BETA1]  = yy * tmp_beta;
    var[bssn::VAR::U_BETA2]  = zz * tmp_beta;

    var[bssn::VAR::U_GT0]    = 0.0;
    var[bssn::VAR::U_GT1]    = 0.0;
    var[bssn::VAR::U_GT2]    = 0.0;

    var[bssn::VAR::U_B0]     = 0.0;
    var[bssn::VAR::U_B1]     = 0.0;
    var[bssn::VAR::U_B2]     = 0.0;

    var[bssn::VAR::U_SYMGT0] = 1.0;  // XX
    var[bssn::VAR::U_SYMGT1] = 0.0;  // XY
    var[bssn::VAR::U_SYMGT2] = 0.0;  // XZ
    var[bssn::VAR::U_SYMGT3] = 1.0;  // YY
    var[bssn::VAR::U_SYMGT4] = 0.0;  // YZ
    var[bssn::VAR::U_SYMGT5] = 1.0;  // ZZ

    var[bssn::VAR::U_SYMAT0] = tmp_Atilde * ((xx / rbar) * (xx / rbar) - third);
    var[bssn::VAR::U_SYMAT1] = tmp_Atilde * xx * yy / (rbar * rbar);
    var[bssn::VAR::U_SYMAT2] = tmp_Atilde * xx * zz / (rbar * rbar);
    var[bssn::VAR::U_SYMAT3] = tmp_Atilde * ((yy / rbar) * (yy / rbar) - third);
    var[bssn::VAR::U_SYMAT4] = tmp_Atilde * yy * zz / (rbar * rbar);
    var[bssn::VAR::U_SYMAT5] = tmp_Atilde * ((zz / rbar) * (zz / rbar) - third);
}

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
void bndcnd(double h, double &x, double y[], double dydx[]) {
    double bh_mass   = 7.5;
    double alpha_c   = 0.16227766016837933200;
    double bigR_c    = 1.5405694150420948330 * bh_mass;
    double r_c       = 0.30405997036 * bh_mass;

    double alpha     = alpha_c;
    double bigR      = bigR_c;
    double r         = r_c;

    double bigRsq    = bigR * bigR;
    double alphasq   = alpha * alpha;

    double dbigR_dr  = alpha_c * bigR_c / r_c;

    double df0_dr    = 2.0 * bh_mass / bigRsq * dbigR_dr;

    double tmp1      = 1.5 * bh_mass / bigRsq * dbigR_dr;

    double dalpha_dr = 0.25 * (r * df0_dr + 8.0 * alphasq) / (alpha - 1.0) -
                       0.25 *
                           sqrt(pow((r * df0_dr + 8.0 * alphasq), 2) +
                                32.0 * alpha * (1.0 - alpha) * r * tmp1) /
                           (alpha - 1.0);

    dalpha_dr = dalpha_dr / r;

    double ddbigR_drdr =
        (dalpha_dr * bigR + alpha * dbigR_dr) / r - alpha * bigR / (r * r);

    double ddf0_drdr = -4.0 * bh_mass * dbigR_dr * dbigR_dr / (bigR * bigRsq) +
                       2.0 * bh_mass * ddbigR_drdr / bigRsq;

    double tmp2 = -3.00 * bh_mass * dbigR_dr * dbigR_dr / (bigR * bigRsq) +
                  1.50 * bh_mass * ddbigR_drdr / bigRsq;

    double ddalpha_drdr =
        -2.0 * pow((r * dalpha_dr), 3) -
        3.0 * r * dalpha_dr *
            (2.0 * r * dalpha_dr * (alpha - 1.0) - r * df0_dr) +
        (r * dalpha_dr) * r * r * ddf0_drdr +
        (8.0 * r * dalpha_dr + 4.0 * alpha) *
            (2.0 * alpha * r * dalpha_dr - r * tmp1) +
        4.0 * alpha * (2.0 * pow((r * dalpha_dr), 2) - r * r * tmp2);

    ddalpha_drdr =
        ddalpha_drdr / (6.0 * pow(r, 3) * (alpha - 1.0) * dalpha_dr -
                        2.0 * pow(r, 3) * df0_dr - 8.0 * r * r * alphasq);

    x       = r + h;

    y[0]    = alpha + h * dalpha_dr + 0.5 * h * h * ddalpha_drdr;
    y[1]    = bigR + h * dbigR_dr + 0.5 * h * h * ddbigR_drdr;

    dydx[0] = dalpha_dr + h * ddalpha_drdr;
    dydx[1] = dbigR_dr + h * ddbigR_drdr;
}

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
void derivs(double x, double y[], double dydx[]) {
    double alpha            = y[0];
    double bigR             = y[1];
    double r                = x;

    double bh_mass          = 7.5;

    double alpha_sq         = alpha * alpha;

    double alpha_sq_minus_1 = alpha_sq - 1.0;
    double M_over_R         = bh_mass / bigR;

    double dalpha_dr = 4.0 * alpha * (alpha_sq_minus_1 + 1.5 * M_over_R) /
                       (alpha_sq_minus_1 - 2.0 * alpha + 2.0 * M_over_R) / r;

    double dbigR_dr = alpha * bigR / r;

    dydx[0]         = dalpha_dr;
    dydx[1]         = dbigR_dr;
}

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
void hunt(double xx[], int n, double x, int *jlo) {
    unsigned long jm, jhi, inc;
    int ascnd;

    ascnd = (xx[n - 1] > xx[0]);
    // if (*jlo < 0 || *jlo > n)
    if (*jlo < 0 || *jlo > n - 1) {
        //*jlo=0;
        *jlo = -1;
        // jhi=n-1;
        jhi  = n;
    } else {
        inc = 1;
        if (x >= xx[*jlo] == ascnd) {
            if (*jlo == n - 1) return;
            jhi = (*jlo) + 1;
            while (x >= xx[jhi] == ascnd) {
                *jlo = jhi;
                inc += inc;
                jhi = (*jlo) + inc;
                if (jhi > n - 1) {
                    jhi = n;
                    break;
                }
            }
        } else {
            if (*jlo == 0) {
                *jlo = -1;
                return;
            }
            jhi = (*jlo)--;
            while (x < xx[*jlo] == ascnd) {
                jhi = (*jlo);
                inc <<= 1;
                if (inc >= jhi) {
                    *jlo = -1;
                    break;
                } else
                    *jlo = jhi - inc;
            }
        }
    }
    while (jhi - (*jlo) != 1) {
        jm = (jhi + (*jlo)) >> 1;
        if (x > xx[jm] == ascnd)
            *jlo = jm;
        else
            jhi = jm;
    }
} /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
void rkck(double y[], double dydx[], int n, double x, double h, double yout[],
          double yerr[], void (*derivs)(double, double[], double[])) {
    int i;
    static double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875, b21 = 0.2,
                  b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9,
                  b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
                  b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0,
                  b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0,
                  b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
                  c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0,
                  c6 = 512.0 / 1771.0, dc5 = -277.0 / 14336.0;
    double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
           dc4 = c4 - 13525.0 / 55296.0, dc6 = c6 - 0.25;

    double *ak2   = new double[n];
    double *ak3   = new double[n];
    double *ak4   = new double[n];
    double *ak5   = new double[n];
    double *ak6   = new double[n];
    double *ytemp = new double[n];

    for (i = 0; i < n; i++) ytemp[i] = y[i] + b21 * h * dydx[i];
    (*derivs)(x + a2 * h, ytemp, ak2);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    (*derivs)(x + a3 * h, ytemp, ak3);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    (*derivs)(x + a4 * h, ytemp, ak4);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] +
                               b54 * ak4[i]);
    (*derivs)(x + a5 * h, ytemp, ak5);
    for (i = 0; i < n; i++)
        ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] +
                               b64 * ak4[i] + b65 * ak5[i]);
    (*derivs)(x + a6 * h, ytemp, ak6);
    for (i = 0; i < n; i++)
        yout[i] =
            y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
    for (i = 0; i < n; i++)
        yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] +
                       dc5 * ak5[i] + dc6 * ak6[i]);
    delete[] ytemp;
    delete[] ak6;
    delete[] ak5;
    delete[] ak4;
    delete[] ak3;
    delete[] ak2;
} /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
void rkqs(double y[], double dydx[], int n, double *x, double htry, double eps,
          double yscal[], double *hdid, double *hnext,
          void (*derivs)(double, double[], double[])) {
    const double SAFETY = 0.9;
    const double PGROW  = -0.2;
    const double PSHRNK = -0.25;
    const double ERRCON = 1.89e-4;

    void rkck(double y[], double dydx[], int n, double x, double h,
              double yout[], double yerr[],
              void (*derivs)(double, double[], double[]));
    int i;
    double errmax, h, xnew;

    double *yerr  = new double[n];
    double *ytemp = new double[n];
    h             = htry;
    for (;;) {
        rkck(y, dydx, n, *x, h, ytemp, yerr, derivs);
        errmax = 0.0;
        for (i = 0; i < n; i++)
            errmax = std::max(errmax, fabs(yerr[i] / yscal[i]));
        errmax /= eps;
        if (errmax > 1.0) {
            h = SAFETY * h * pow(errmax, PSHRNK);
            if (h < 0.1 * h) h *= 0.1;
            xnew = (*x) + h;
            if (xnew == *x)
                std::cerr << "stepsize underflow in rkqs" << std::endl;
            continue;
        } else {
            if (errmax > ERRCON)
                *hnext = SAFETY * h * pow(errmax, PGROW);
            else
                *hnext = 5.0 * h;
            *x += (*hdid = h);
            for (i = 0; i < n; i++) y[i] = ytemp[i];
            break;
        }
    }
    delete[] ytemp;
    delete[] yerr;
} /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
void odeint(double ystart[], int nvar, double x1, double x2, double eps,
            double h1, double hmin, int *nok, int *nbad,
            void (*derivs)(double, double[], double[]),
            void (*rkqs)(double[], double[], int, double *, double, double,
                         double[], double *, double *,
                         void (*)(double, double[], double[])),
            int kount) {
    int nstp, i;
    double xsav, x, hnext, hdid, h;
    const int MAXSTP  = 10000;
    const double TINY = 1.0e-30;

    double *yscal     = new double[nvar];
    double *y         = new double[nvar];
    double *dydx      = new double[nvar];
    x                 = x1;
    h                 = copysign(h1, x2 - x1);
    *nok = (*nbad) = kount = 0;
    for (i = 0; i < nvar; i++) y[i] = ystart[i];
    if (KMAX > 0) xsav = x - dxsav * 2.0;
    for (nstp = 0; nstp < MAXSTP; nstp++) {
        // std::cout<<"odeint: nstp="<<nstp<<", kount="<<kount<<std::endl;
        (*derivs)(x, y, dydx);
        for (i = 0; i < nvar; i++)
            yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
        if (KMAX > 0 && kount < KMAX - 1 && fabs(x - xsav) > fabs(dxsav)) {
            xp[kount] = x;
            for (i = 0; i < nvar; i++) yp[i][kount] = y[i];
            xsav = x;
            kount++;
        }
        if ((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x;
        (*rkqs)(y, dydx, nvar, &x, h, eps, yscal, &hdid, &hnext, derivs);
        if (hdid == h)
            ++(nok);
        else
            ++(nbad);
        if ((x - x2) * (x2 - x1) >= 0.0) {
            for (i = 0; i < nvar; i++) ystart[i] = y[i];
            if (KMAX) {
                xp[kount] = x;
                for (i = 0; i < nvar; i++) yp[i][kount] = y[i];
                kount++;
            }
            delete[] dydx;
            delete[] y;
            delete[] yscal;
            return;
        }
        if (fabs(hnext) <= hmin)
            std::cerr << "Step size too small in odeint" << std::endl;
        h = hnext;
    }
    std::cerr << "Too many steps in routine odeint" << std::endl;
}
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
double interpolation3(double xp[], double yp[], int np, double xb,
                      int *n_nearest_pt) {
    int k;     /* index of 1st point */
    int m = 4; /* degree of interpolation */
    double y;  /* intermediate value */

    hunt(xp, np, xb, n_nearest_pt);

    k = std::min(std::max((*n_nearest_pt) - (m - 1) / 2, 1), np + 1 - m);

    double DBL_EPSILON = 1.e-12;

    if (xb == xp[k] || xb == xp[k + 1] || xb == xp[k + 2] || xb == xp[k + 3]) {
        xb += 3.0 * DBL_EPSILON;
    }

    y = (xb - xp[k + 1]) * (xb - xp[k + 2]) * (xb - xp[k + 3]) * yp[k] /
            ((xp[k] - xp[k + 1]) * (xp[k] - xp[k + 2]) * (xp[k] - xp[k + 3]))

        + (xb - xp[k]) * (xb - xp[k + 2]) * (xb - xp[k + 3]) * yp[k + 1] /
              ((xp[k + 1] - xp[k]) * (xp[k + 1] - xp[k + 2]) *
               (xp[k + 1] - xp[k + 3]))

        + (xb - xp[k]) * (xb - xp[k + 1]) * (xb - xp[k + 3]) * yp[k + 2] /
              ((xp[k + 2] - xp[k]) * (xp[k + 2] - xp[k + 1]) *
               (xp[k + 2] - xp[k + 3]))

        + (xb - xp[k]) * (xb - xp[k + 1]) * (xb - xp[k + 2]) * yp[k + 3] /
              ((xp[k + 3] - xp[k]) * (xp[k + 3] - xp[k + 1]) *
               (xp[k + 3] - xp[k + 2]));

    return (y);
}

/*---------------------------------------------------------------------
 *
 *
 *
 *---------------------------------------------------------------------*/
double interpolation4(double xx[], double yy[], int np, double xb,
                      int *n_nearest_pt) {
    int k;     /* index of 1st point */
    int m = 5; /* degree of interpolation */
    double y;  /* intermediate value */

    hunt(xx, np, xb, n_nearest_pt);

    k = std::min(std::max((*n_nearest_pt) - (m - 1) / 2, 1), np + 1 - m);

#if 0
        if ( fabs( xb - 1.0 ) < 1.e-13 )
        { // xb is 1.0 so just return with the corresponding y value -- no interp necessary 
            y = yy[np] ;
            return(y) ;
        }
        if ( fabs( xb ) < 1.e-13 )
        { // xb is zero so just return with the corresponding y value -- no interp necessary 
            y = yy[0] ;
            return(y) ;
        }
#endif

    double DBL_EPSILON = 1.e-12;

    if (xb == xx[k] || xb == xx[k + 1] || xb == xx[k + 2] || xb == xx[k + 3] ||
        xb == xx[k + 4]) {
        xb += 3.0 * DBL_EPSILON;
    }

    double xtmp0   = xb - xx[k];
    double xtmp1   = xb - xx[k + 1];
    double xtmp2   = xb - xx[k + 2];
    double xtmp3   = xb - xx[k + 3];
    double xtmp4   = xb - xx[k + 4];
    double xdiff01 = xx[k] - xx[k + 1];
    double xdiff02 = xx[k] - xx[k + 2];
    double xdiff03 = xx[k] - xx[k + 3];
    double xdiff04 = xx[k] - xx[k + 4];
    double xdiff12 = xx[k + 1] - xx[k + 2];
    double xdiff13 = xx[k + 1] - xx[k + 3];
    double xdiff14 = xx[k + 1] - xx[k + 4];
    double xdiff23 = xx[k + 2] - xx[k + 3];
    double xdiff24 = xx[k + 2] - xx[k + 4];
    double xdiff34 = xx[k + 3] - xx[k + 4];

    y              = xtmp1 * xtmp2 * xtmp3 * xtmp4 * yy[k] /
            (xdiff01 * xdiff02 * xdiff03 * xdiff04)

        - xtmp0 * xtmp2 * xtmp3 * xtmp4 * yy[k + 1] /
              (xdiff01 * xdiff12 * xdiff13 * xdiff14)

        + xtmp0 * xtmp1 * xtmp3 * xtmp4 * yy[k + 2] /
              (xdiff02 * xdiff12 * xdiff23 * xdiff24)

        - xtmp0 * xtmp1 * xtmp2 * xtmp4 * yy[k + 3] /
              (xdiff03 * xdiff13 * xdiff23 * xdiff34) +
        xtmp0 * xtmp1 * xtmp2 * xtmp3 * yy[k + 4] /
            (xdiff04 * xdiff14 * xdiff24 * xdiff34);

    return (y);
}

}  // end of namespace trumpet_data
