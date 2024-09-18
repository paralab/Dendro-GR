#include "recon.h"

#include "geom.h"
#include "parameters.h"

using namespace fluid;
using namespace recon;

// none::reconstruct {{{
void recon::none::reconstruct(const int n, const double* const RESTRICT v,
                              double* const RESTRICT vl,
                              double* const RESTRICT vr) {
    vl[0] = v[0];
    for (int i = 0; i < n - 1; i++) {
        vl[i + 1] = v[i];
        vr[i]     = v[i];
    }
    vr[n - 1] = v[n - 1];
}
// }}}

// minmod::reconstruct {{{
void recon::minmod::reconstruct(const int n, const double* const RESTRICT v,
                                double* const RESTRICT vl,
                                double* const RESTRICT vr) {
    vr[0] = v[0];
    vl[0] = v[0];

    for (int i = 1; i < n - 1; i++) {
        double df1 = (v[i] - v[i - 1]);
        double df2 = (v[i + 1] - v[i]);

        double slp = 0.5 * (copysign(1.0, df1) + copysign(1.0, df2)) *
                     std::min(fabs(df1), fabs(df2));
        vl[i + 1] = v[i] + 0.5 * slp;
        vr[i]     = v[i] - 0.5 * slp;
    }
    vl[n - 1] = v[n - 1];
    vr[n - 1] = v[n - 1];
}
// }}}

// minmod::reconstructPt {{{
void recon::minmod::reconstructPt(const int n, const int i,
                                  const double* const RESTRICT v,
                                  double* const RESTRICT vl,
                                  double* const RESTRICT vr) {
    if (i == 0) {
        vl[0] = v[0];
        vr[0] = v[0];

        if (i < n - 1) {
            vl[1] = v[0];
        }
    } else if (i > 0 && i < n - 1) {
        double df1 = v[i] - v[i - 1];
        double df2 = v[i + 1] - v[i];
        double slp = 0.5 * (copysign(1.0, df1) + copysign(1.0, df2)) *
                     std::min(fabs(df1), fabs(df2));
        vl[i + 1] = v[i] + 0.5 * slp;
        vr[i]     = v[i] - 0.5 * slp;
    } else if (i == n - 1) {
        vl[i] = v[i];
        vr[i] = v[i];
    }
}
// }}}

// weno3::reconstruct {{{
// FIXME: There's probably a bug somewhere in here; it doesn't seem to
// be as stable as it should be.
void recon::weno3::reconstruct(const int n, const double* const RESTRICT u,
                               double* const RESTRICT ul,
                               double* const RESTRICT ur) {
    // WENO3 uses three different interpolations for the reconstructed point:
    // one for a three-point stencil on the left, one for a three-point stencil
    // on the right, and one for a centered stencil. It then weights them
    // appropriately to achieve the most accurate reconstruction based on the
    // behavior of the points in the stencil.

    // Here are the conventions used:
    // We determine what is left and right in terms of their position
    // relative to the interface at u[i-1/2]. Therefore ur[i] is on the right
    // of this interface( so in the same cell as u[i], but located on the
    // left), and ul[i] is on the left of this interface (in the same cell as
    // u[i-1], but to its right).
    // These are the interpolations used:
    //
    // u1[i-1/2] = 1/3*u[i  ] + 5/6*u[i-1] - 1/6 *u[i-2]
    // u2[i-1/2] = 1/3*u[i-1] + 5/6*u[i  ] - 1/6 *u[i+1]
    // u3[i-1/2] = 1/3*u[i+2] - 7/6*u[i+1] + 11/6*u[i  ]
    //
    // u1[i+1/2] = 1/3*u[i-2] - 7/6*u[i-1] + 11/6*u[i  ]
    // u2[i+1/2] = 1/3*u[i+1] + 5/6*u[i  ] - 1/6 *u[i-1]
    // u3[i+1/2] = 1/3*u[i  ] + 5/6*u[i+1] - 1/6 *u[i+2]

    // Polynomial weights
    double pw[] = {1.0 / 6.0, 5.0 / 6.0, 1.0 / 3.0, 7.0 / 6.0, 11.0 / 6.0};
    // Linear weights
    double lw[] = {1.0 / 10.0, 3.0 / 5.0, 3.0 / 10.0};

    double fac  = 13.0 / 12.0;
    // Interpolated points.
    double u1, u2, u3;
    // Smoothness indicators.
    double b1, b2, b3;
    // Nonlinear weights
    double w1, w2, w3;
    double wt1, wt2, wt3;

    for (int i = 2; i < n - 2; i++) {
        // Start with the right side.

        // Construct the interpolations.
        u1 = pw[2] * u[i] + pw[1] * u[i - 1] - pw[0] * u[i - 2];
        u2 = pw[2] * u[i - 1] + pw[1] * u[i] - pw[0] * u[i + 1];
        u3 = pw[2] * u[i + 2] - pw[3] * u[i + 1] + pw[4] * u[i];

        // Construct the smoothness indicators.
        b1 = fac * (u[i - 2] - 2.0 * u[i - 1] + u[i]) *
                 (u[i - 2] - 2.0 * u[i - 1] + u[i]) +
             0.25 * (u[i - 2] - 4.0 * u[i - 1] + 3.0 * u[i]) *
                 (u[i - 2] - 4.0 * u[i - 1] + 3.0 * u[i]);
        b2 = fac * (u[i - 1] - 2.0 * u[i] + u[i + 1]) *
                 (u[i - 1] - 2.0 * u[i] + u[i + 1]) +
             0.25 * (u[i - 1] - u[i + 1]) * (u[i - 1] - u[i + 1]);
        b3 = fac * (u[i] - 2.0 * u[i + 1] + u[i + 2]) *
                 (u[i] - 2.0 * u[i + 1] + u[i + 2]) +
             0.25 * (u[i + 2] - 4.0 * u[i + 1] + 3.0 * u[i]) *
                 (u[i + 2] - 4.0 * u[i + 1] + 3.0 * u[i]);

        // Construct the nonlinear weights.
        double tau5 = fabs(b1 - b3);
        wt1 = lw[0] * (1.0 + (tau5 / (1.0e-15 + b1)) * (tau5 / (1.0e-15 + b1)));
        wt2 = lw[1] * (1.0 + (tau5 / (1.0e-15 + b2)) * (tau5 / (1.0e-15 + b2)));
        wt3 = lw[2] * (1.0 + (tau5 / (1.0e-15 + b3)) * (tau5 / (1.0e-15 + b3)));
        double isum = 1.0 / (wt1 + wt2 + wt3);
        w1          = wt1 * isum;
        w2          = wt2 * isum;
        w3          = wt3 * isum;

        // Now construct the final reconstructed value as a linear combination
        // of the weights.
        ur[i]       = w1 * u1 + w2 * u2 + w3 * u3;

        // Now we need to reconstruct ul[i+1].
        u1          = pw[2] * u[i - 2] - pw[3] * u[i - 1] + pw[4] * u[i];
        u2          = pw[2] * u[i + 1] + pw[1] * u[i] - pw[0] * u[i - 1];
        u3          = pw[2] * u[i] + pw[1] * u[i + 1] - pw[0] * u[i + 2];

        // Our weights and smoothness indicators haven't changed, so all we need
        // to do is reconstruct ul based on these new interpolations.
        ul[i + 1]   = w1 * u1 + w2 * u2 + w3 * u3;
    }

    // We can't reconstruct points at the edges of the grid. Just assume
    // no reconstruction here.
    ul[0]     = u[0];
    ur[0]     = u[0];
    ul[1]     = u[0];
    ur[1]     = u[1];
    ul[2]     = u[1];

    ul[n - 1] = u[n - 2];
    ur[n - 1] = u[n - 1];
    ur[n - 2] = u[n - 2];
}
// }}}

// weno5::reconstruct {{{
const double weno_eps   = 1.0e-20;
int do_wenoz            = 1;
int do_adaptive_epsilon = 0;
/**
   WENO5 reconstruction operator.
   Supports standard WENO5 (with and without adaptive epsilon), and WENO-Z.
*/
/*********************************************************************************
 *
 *    WENO5 reconstruction
 *
 *    DWN - Unfortunately, not everyone uses my notation. aplus is the
 *          reconstructed state on the left side of the boundary at i+1/2.
 *          aminus is the state on the right side of the boundary at i-1/2.
 *          FIXME: We may need to adjust the index of aplus by one?
 *
 ********************************************************************************/
void recon::weno5::reconstruct(const int nx, const double* const RESTRICT a,
                               double* const RESTRICT aplus,
                               double* const RESTRICT aminus) {
#define A(i_)      (a[ijk[i_]])
#define Aplus(i_)  (aplus[ijk[i_]])
#define Aminus(i_) (aminus[ijk[i_]])
#define SQR(x)     ((x) * (x))

    for (unsigned int i = 2; i < nx - 2; ++i) {
        const unsigned int ijk[5] = {i - 2, i - 1, i, i + 1, i + 2};
        assert(!(do_wenoz && do_adaptive_epsilon) &&
               "Adaptive_epsilon not supported for WENO-Z");

        if (do_wenoz) {
            static const double weno_coeffs[3][5] = {
                {2.0 / 6.0, -7.0 / 6.0, 11.0 / 6.0, 0, 0},
                {0, -1.0 / 6.0, 5.0 / 6.0, 2.0 / 6.0, 0},
                {0, 0, 2.0 / 6.0, 5.0 / 6.0, -1.0 / 6.0}};

            const double beta1 =
                13.0 / 12.0 * (A(0) - 2.0 * A(1) + A(2)) *
                    (A(0) - 2.0 * A(1) + A(2)) +
                1.0 / 4.0 * SQR(A(0) - 4.0 * A(1) + 3.0 * A(2));
            const double beta2 = 13.0 / 12.0 * (A(1) - 2.0 * A(2) + A(3)) *
                                     (A(1) - 2.0 * A(2) + A(3)) +
                                 1.0 / 4.0 * SQR(A(1) - A(3));
            const double beta3 =
                13.0 / 12.0 * (A(2) - 2.0 * A(3) + A(4)) *
                    (A(2) - 2.0 * A(3) + A(4)) +
                1.0 / 4.0 * SQR(3.0 * A(2) - 4.0 * A(3) + A(4));

            // Compute weights according to WENO-Z algorithm.
            const double wbarplus1 =
                1.0 / 10.0 * (1.0 + fabs(beta1 - beta3) / (weno_eps + beta1));
            const double wbarplus2 =
                3.0 / 5.0 * (1.0 + fabs(beta1 - beta3) / (weno_eps + beta2));
            const double wbarplus3 =
                3.0 / 10.0 * (1.0 + fabs(beta1 - beta3) / (weno_eps + beta3));

            const double wbarminus1 =
                3.0 / 10.0 * (1.0 + fabs(beta1 - beta3) / (weno_eps + beta1));
            const double wbarminus2 =
                3.0 / 5.0 * (1.0 + fabs(beta1 - beta3) / (weno_eps + beta2));
            const double wbarminus3 =
                1.0 / 10.0 * (1.0 + fabs(beta1 - beta3) / (weno_eps + beta3));

            const double iwbarplussum =
                1.0 / (wbarplus1 + wbarplus2 + wbarplus3);

            const double wplus1 = wbarplus1 * iwbarplussum;
            const double wplus2 = wbarplus2 * iwbarplussum;
            const double wplus3 = wbarplus3 * iwbarplussum;

            const double iwbarminussum =
                1.0 / (wbarminus1 + wbarminus2 + wbarminus3);

            const double wminus1 = wbarminus1 * iwbarminussum;
            const double wminus2 = wbarminus2 * iwbarminussum;
            const double wminus3 = wbarminus3 * iwbarminussum;

            // Calculate the reconstruction
            Aplus(3)             = 0;
            Aminus(2)            = 0;
            for (int j = 0; j < 5; ++j) {
                Aplus(3) +=
                    (wplus1 * weno_coeffs[0][j] + wplus2 * weno_coeffs[1][j] +
                     wplus3 * weno_coeffs[2][j]) *
                    A(j);
                Aminus(2) += (wminus1 * weno_coeffs[2][4 - j] +
                              wminus2 * weno_coeffs[1][4 - j] +
                              wminus3 * weno_coeffs[0][4 - j]) *
                             A(j);
            }
        } else {
            static const double beta_shu[3][6] = {
                {4.0 / 3.0, -19.0 / 3.0, 25.0 / 3.0, 11.0 / 3.0, -31.0 / 3.0,
                 10.0 / 3.0},
                {4.0 / 3.0, -13.0 / 3.0, 13.0 / 3.0, 5.0 / 3.0, -13.0 / 3.0,
                 4.0 / 3.0},
                {10.0 / 3.0, -31.0 / 3.0, 25.0 / 3.0, 11.0 / 3.0, -19.0 / 3.0,
                 4.0 / 3.0}};
            static const double weno_coeffs[3][5] = {
                {3.0 / 8.0, -5.0 / 4.0, 15.0 / 8.0, 0, 0},
                {0, -1.0 / 8.0, 3.0 / 4.0, 3.0 / 8.0, 0},
                {0, 0, 3.0 / 8.0, 3.0 / 4.0, -1.0 / 8.0}};

            // Compute smoothness indicators
            // This is from Tchekhovskoy et al 2007 (WHAM code paper).
            double beta1 =
                beta_shu[0][0] * (A(0) * A(0)) + beta_shu[0][1] * A(0) * A(1) +
                beta_shu[0][2] * A(1) * A(1) + beta_shu[0][3] * A(0) * A(2) +
                beta_shu[0][4] * A(1) * A(2) + beta_shu[0][5] * A(2) * A(2);

            double beta2 =
                beta_shu[1][0] * A(1) * A(1) + beta_shu[1][1] * A(1) * A(2) +
                beta_shu[1][2] * A(2) * A(2) + beta_shu[1][3] * A(1) * A(3) +
                beta_shu[1][4] * A(2) * A(3) + beta_shu[1][5] * A(3) * A(3);

            double beta3 =
                beta_shu[2][0] * A(2) * A(2) + beta_shu[2][1] * A(2) * A(3) +
                beta_shu[2][2] * A(3) * A(3) + beta_shu[2][3] * A(2) * A(4) +
                beta_shu[2][4] * A(3) * A(4) + beta_shu[2][5] * A(4) * A(4);

            if (do_adaptive_epsilon) {
                const double vnorm = (A(0) * A(0) + A(1) * A(1) + A(2) * A(2) +
                                      A(3) * A(3) + A(4) * A(4));

                beta1 += 100.0 * weno_eps * (vnorm + 1.0);
                beta2 += 100.0 * weno_eps * (vnorm + 1.0);
                beta3 += 100.0 * weno_eps * (vnorm + 1.0);

                const double ibetanorm = 1.0 / (beta1 + beta2 + beta3);

                beta1 *= ibetanorm;
                beta2 *= ibetanorm;
                beta3 *= ibetanorm;
            }

            const double wbarplus1 =
                1.0 / 16.0 / ((weno_eps + beta1) * (weno_eps + beta1));
            const double wbarplus2 =
                5.0 / 8.0 / ((weno_eps + beta2) * (weno_eps + beta2));
            const double wbarplus3 =
                5.0 / 16.0 / ((weno_eps + beta3) * (weno_eps + beta3));

            const double iwbarplussum =
                1.0 / (wbarplus1 + wbarplus2 + wbarplus3);

            const double wplus1 = wbarplus1 * iwbarplussum;
            const double wplus2 = wbarplus2 * iwbarplussum;
            const double wplus3 = wbarplus3 * iwbarplussum;

            const double wbarminus1 =
                5.0 / 16.0 / ((weno_eps + beta1) * (weno_eps + beta1));
            const double wbarminus2 =
                5.0 / 8.0 / ((weno_eps + beta2) * (weno_eps + beta2));
            const double wbarminus3 =
                1.0 / 16.0 / ((weno_eps + beta3) * (weno_eps + beta3));

            const double iwbarminussum =
                1.0 / (wbarminus1 + wbarminus2 + wbarminus3);

            const double wminus1 = wbarminus1 * iwbarminussum;
            const double wminus2 = wbarminus2 * iwbarminussum;
            const double wminus3 = wbarminus3 * iwbarminussum;

            // Calculate the reconstruction
            Aplus(3)             = 0;
            Aminus(2)            = 0;
            for (unsigned int j = 0; j < 5; ++j) {
                Aplus(3) +=
                    (wplus1 * weno_coeffs[0][j] + wplus2 * weno_coeffs[1][j] +
                     wplus3 * weno_coeffs[2][j]) *
                    A(j);
                Aminus(2) += (wminus1 * weno_coeffs[2][4 - j] +
                              wminus2 * weno_coeffs[1][4 - j] +
                              wminus3 * weno_coeffs[0][4 - j]) *
                             A(j);
            }
        }
    }

    aplus[0]       = a[0];
    aplus[1]       = a[0];
    aplus[2]       = a[1];
    aminus[0]      = a[0];
    aminus[1]      = a[1];

    aminus[nx - 2] = a[nx - 2];
    aminus[nx - 1] = a[nx - 1];
    aplus[nx - 1]  = a[nx - 2];
}
// }}}

// mp5::reconstruct {{{
int do_MP5_adaptive_epsilon = 1;
/*********************************************************************************
 *
 *    MP5 reconstruction
 *
 *    DWN - Unfortunately, not everyone uses my notation. aplus is the
 *          reconstructed state on the left side of the boundary at i+1/2.
 *          aminus is the state on the right side of the boundary at i-1/2.
 *          FIXME: We may need to adjust the index of aplus by one?
 *
 ********************************************************************************/
void recon::mp5::reconstruct(const int nx, const double* const RESTRICT a,
                             double* const RESTRICT aplus,
                             double* const RESTRICT aminus) {
#define A(i_)      (a[ijk[i_]])
#define Aplus(i_)  (aplus[ijk[i_]])
#define Aminus(i_) (aminus[ijk[i_]])

    const double mp5_alpha = 4.0;
    const double mp5_eps   = 1.0e-10;

    for (int i = 2; i < nx - 2; ++i) {
        const int ijk[5] = {i - 2, i - 1, i, i + 1, i + 2};

        if (!do_MP5_adaptive_epsilon) {
            Aplus(3) =
                MP5(A(0), A(1), A(2), A(3), A(4), 1.0, mp5_eps, mp5_alpha);
            Aminus(2) =
                MP5(A(4), A(3), A(2), A(1), A(0), 1.0, mp5_eps, mp5_alpha);
        } else {
            const double anorm = sqrt(A(0) * A(0) + A(1) * A(1) + A(2) * A(2) +
                                      A(3) * A(3) + A(4) * A(4));

            Aplus(3) =
                MP5(A(0), A(1), A(2), A(3), A(4), anorm, mp5_eps, mp5_alpha);
            Aminus(2) =
                MP5(A(4), A(3), A(2), A(1), A(0), anorm, mp5_eps, mp5_alpha);
        }
    }

    aplus[0]       = a[0];
    aplus[1]       = a[0];
    aplus[2]       = a[1];
    aminus[0]      = a[0];
    aminus[1]      = a[1];

    aminus[nx - 2] = a[nx - 2];
    aminus[nx - 1] = a[nx - 1];
    aplus[nx - 1]  = a[nx - 2];
}
// }}}

// reconstructvars {{{
void recon::reconstructvars(double* vl[], double* vr[], double* v1[],
                            double* wv[], const int dir, const double pos_in[],
                            const double dx, const unsigned int* sz) {
#define NBUFFER      3
#define STENCIL_SIZE 3

    // memory for reconstruction
    /*static int *mask;
    static double *s1[FLUID_NUM_PRIM_VARS], *sl[FLUID_NUM_PRIM_VARS],
    *sr[FLUID_NUM_PRIM_VARS], *swv[3], *sxi[3]; int sshp[3]; int xdir = (dir +
    1) % 3; int ydir = (dir + 2) % 3; sshp[dir] = STENCIL_SIZE; sshp[xdir] = 1;
    sshp[ydir] = 1;*/

    // static int first_call = 1;
    /*unsigned int nx = sz[0];
    unsigned int ny = sz[1];
    unsigned int nz = sz[2];
    unsigned int dim = max(nx,max(ny,nz));*/
    unsigned int num = FLUID_NUM_PRIM_VARS;

    /*static unsigned int size = 2*dim;
    if(first_call == 1){
      mask = (int *)malloc(size*sizeof(int));
      for(int m = 0; m < dim; m++){
        s1[m] = (double*)malloc(STENCIL_SIZE*sizeof(double));
        sl[m] = (double*)malloc(STENCIL_SIZE*sizeof(double));
        sr[m] = (double*)malloc(STENCIL_SIZE*sizeof(double));
      }
      for(int m = 0; m < 3; m++){
        swv[m] = (double*)malloc(STENCIL_SIZE * sizeof(double));
        sxi[m] = (double*)malloc(STENCIL_SIZE * sizeof(double));
      }
      first_call = 0;
    }

    // If the dimension size exceeds the allocated memory, we need to reallocate
    the mask. if(dim > size){ size = dim * 2; free(mask); mask = (int
    *)malloc(size*sizeof(int));
    }*/

    int n            = sz[dir];
    // double xi, dx;

    /*if(dir == 0){
      n = nx;
      xi = pmin[0];
      dx = (pmax[0] - pmin[0])/(n-1);
    }
    else if(dir == 1){
      n = ny;
      xi = pmin[1];
      dx = (pmax[1] - pmin[1])/(n-1);
    }
    else if(dir == 2){
      n = nz;
      xi = pmin[2];
      dx = (pmax[2] - pmin[2])/(n-1);
    }*/

    /*for(int k = 0; k < n; k++){
      if(v1[PVAR::V_RHO][k] < 1.0e-13){
        printf("...reconstructvars [%d]: rho[%d] =
    %g\n",dir,k,v1[PVAR::V_RHO][k]);
      }
    }*/

    recon_help(vl, vr, v1, wv, dir, pos_in, dx, sz);

    /* Apply floor to reconstructed variables */
    const double vacuum_reset     = FLUID_VACUUM_RESET;
    const double vacuum_tau_reset = FLUID_VACUUM_TAU_RESET;

#define USE_LOWER_ORDER_RECON 1
#if USE_LOWER_ORDER_RECON == 1
    // Lower order recon based on pointwise methods {{{
    void (*fallback)(const int n, const int i, const double* const RESTRICT v,
                     double* const RESTRICT vl, double* const RESTRICT vr) =
        &minmod::reconstructPt;

    for (int i = 0; i < n; i++) {
        if (vl[PVAR::V_RHO][i] < vacuum_reset) {
            fallback(n, (i > 0) ? i - 1 : 0, v1[PVAR::V_RHO], vl[PVAR::V_RHO],
                     vr[PVAR::V_RHO]);
        }
        if (vr[PVAR::V_RHO][i] < vacuum_reset) {
            fallback(n, i, v1[PVAR::V_RHO], vl[PVAR::V_RHO], vr[PVAR::V_RHO]);
        }
        if (vl[PVAR::V_P][i] < vacuum_tau_reset) {
            fallback(n, (i > 0) ? i - 1 : 0, v1[PVAR::V_P], vl[PVAR::V_P],
                     vr[PVAR::V_P]);
        }
        if (vr[PVAR::V_P][i] < vacuum_tau_reset) {
            fallback(n, i, v1[PVAR::V_P], vl[PVAR::V_P], vr[PVAR::V_P]);
        }
    }
// }}}
#else
    // No lower order recon, just floor {{{
    for (int i = 0; i < n; i++) {
        if (vl[PVAR::V_RHO][i] < vacuum_reset) {
            vl[PVAR::V_RHO][i] = vacuum_reset;
            vl[PVAR::V_VX][i]  = 0.0;
            vl[PVAR::V_VY][i]  = 0.0;
            vl[PVAR::V_VZ][i]  = 0.0;
        }
        if (vr[PVAR::V_RHO][i] < vacuum_reset) {
            vr[PVAR::V_RHO][i] = vacuum_reset;
            vr[PVAR::V_VX][i]  = 0.0;
            vr[PVAR::V_VY][i]  = 0.0;
            vr[PVAR::V_VZ][i]  = 0.0;
        }
        // vl[PVAR::V_RHO][i] = fmax(vl[PVAR::V_RHO][i], vacuum_reset);
        // vr[PVAR::V_RHO][i] = fmax(vr[PVAR::V_RHO][i], vacuum_reset);
        vl[PVAR::V_P][i] = fmax(vl[PVAR::V_P][i], vacuum_tau_reset);
        vr[PVAR::V_P][i] = fmax(vr[PVAR::V_P][i], vacuum_tau_reset);
    }
// }}}
#endif
}
// }}}

// recon_help {{{
void recon::recon_help(double* vl[], double* vr[], double* v1[], double* wv[],
                       const int dir, const double pos_in[], const double dx,
                       const unsigned int* sz) {
    void (*recon)(const int n, const double* const RESTRICT v,
                  double* const RESTRICT vl, double* const RESTRICT vr) =
        weno3::reconstruct;

    if (FLUID_RECON_METHOD == 4) {
        recon = &mp5::reconstruct;
    } else if (FLUID_RECON_METHOD == 3) {
        recon = &weno5::reconstruct;
    } else if (FLUID_RECON_METHOD == 2) {
        recon = &weno3::reconstruct;
    } else if (FLUID_RECON_METHOD == 1) {
        recon = &minmod::reconstruct;
    } else if (FLUID_RECON_METHOD == 0) {
        recon = &none::reconstruct;
    } else {
        printf("recon_help: Unknown reconstruction method.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }

    double pos[3];
    int n;
    double xm;

    switch (dir) {
        case 0:
            n  = sz[0];
            xm = pos_in[0];
            break;
        case 1:
            n  = sz[1];
            xm = pos_in[1];
            break;
        case 2:
            n  = sz[2];
            xm = pos_in[2];
            break;
        default:
            fprintf(stderr, "recon_help:  Unknown direction.\n");
            exit(2);
    }

    recon(n, v1[PVAR::V_RHO], vl[PVAR::V_RHO], vr[PVAR::V_RHO]);

    recon(n, v1[PVAR::V_P], vl[PVAR::V_P], vr[PVAR::V_P]);

    // Reconstruct the four-velocity.
    cal_wv(wv, v1, dir, pos_in, dx, sz);

    recon(n, wv[0], vl[PVAR::V_VX], vr[PVAR::V_VX]);
    recon(n, wv[1], vl[PVAR::V_VY], vr[PVAR::V_VY]);
    recon(n, wv[2], vl[PVAR::V_VZ], vr[PVAR::V_VZ]);

    for (int i = 0; i < n; i++) {
        pos[dir] = xm + dx * ((double)i - 0.5);
        double wvl_pt[3], wvr_pt[3];

        wvl_pt[0] = vl[PVAR::V_VX][i];
        wvl_pt[1] = vl[PVAR::V_VY][i];
        wvl_pt[2] = vl[PVAR::V_VZ][i];

        wvr_pt[0] = vr[PVAR::V_VX][i];
        wvr_pt[1] = vr[PVAR::V_VY][i];
        wvr_pt[2] = vr[PVAR::V_VZ][i];

        double gd[3][3], detg;
        metric(gd, &detg, pos);

        double wl         = sqrt(1.0 + square_vector(wvl_pt, gd));
        double wr         = sqrt(1.0 + square_vector(wvr_pt, gd));
        wl                = fmax(1.0, wl);
        wr                = fmax(1.0, wr);

        vl[PVAR::V_VX][i] = wvl_pt[0] / wl;
        vl[PVAR::V_VY][i] = wvl_pt[1] / wl;
        vl[PVAR::V_VZ][i] = wvl_pt[2] / wl;

        vr[PVAR::V_VX][i] = wvr_pt[0] / wr;
        vr[PVAR::V_VY][i] = wvr_pt[1] / wr;
        vr[PVAR::V_VZ][i] = wvr_pt[2] / wr;
    }

    pos[0] = pos_in[0];
    pos[1] = pos_in[1];
    pos[2] = pos_in[2];
}
// }}}

// cal_wv {{{
int recon::cal_wv(double* wv[], double* v[], const int dir,
                  const double pos_in[], const double dx,
                  const unsigned int* sz) {
    double pos[3];
    double n;

    switch (dir) {
        case 0:
            pos[1] = pos_in[1];
            pos[2] = pos_in[2];
            n      = sz[0];
            break;
        case 1:
            pos[0] = pos_in[0];
            pos[2] = pos_in[2];
            n      = sz[1];
            break;
        case 2:
            pos[0] = pos_in[0];
            pos[1] = pos_in[1];
            n      = sz[2];
            break;
    }

    for (int i = 0; i < n; i++) {
        pos[dir] = pos_in[dir] + dx * (i - 0.5);

        double gd[3][3], detg;
        metric(gd, &detg, pos);

        double vu[3];
        vu[0]      = v[PVAR::V_VX][i];
        vu[1]      = v[PVAR::V_VY][i];
        vu[2]      = v[PVAR::V_VZ][i];
        double vsq = square_vector(vu, gd);

        // check to account for round-off errors
        double t1  = vsq - 1.0;
        while (t1 > 0.0 && t1 < 1.0e-15) {
            for (int m = 0; m < 3; m++) {
                if (vu[m] < 0.0) {
                    vu[m] += 5.0e-16;
                } else {
                    vu[m] -= 5.0e-16;
                }
            }
            vsq = square_vector(vu, gd);
            t1  = vsq - 1.0;
        }

        // potential problem.
        if (vsq >= 1.0) {
            printf(
                "## cal_wv: before reconstruction, vsq >= 1.0 vsq = %25.20e.",
                vsq);
            printf(" v = ( %g, %g, %g )", vu[0], vu[1], vu[2]);
            printf(" g = diag( %g, %g, %g )\n", gd[0][0], gd[1][1], gd[2][2]);
            double t1 = sqrt(1.0 / vsq - 6.0e-16);
            for (int m = 0; m < 3; m++) {
                vu[m] *= t1;
            }
            vsq = square_vector(vu, gd);
            if (vsq >= 1.0) {
                printf("#### cal_wv: problem with vsq.\n");
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }

        double W = fmax(1.0 / sqrt(1.0 - vsq), 1.0);

        wv[0][i] = W * vu[0];
        wv[1][i] = W * vu[1];
        wv[2][i] = W * vu[2];
    }
}
// }}}
