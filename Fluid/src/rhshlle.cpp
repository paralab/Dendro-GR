#include "rhshlle.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <new>

#include "derivs.h"
#include "fluid.h"
#include "fluidMath.h"
#include "fluidUtils.h"
#include "geom.h"
#include "recon.h"

using namespace std;
using namespace fluid;
using namespace hlle;

// #define DEBUG_FLUX
// #define LF_FLUX

// hllerhs {{{
void hlle::rhshlle(double *dtu[], const double *u[], const double *v[],
                   const double *pmin, const double *pmax,
                   const unsigned int *sz) {
    // Get some variables we're sure to use frequently.
    unsigned int nx         = sz[DIR_X];
    unsigned int ny         = sz[DIR_Y];
    unsigned int nz         = sz[DIR_Z];
    unsigned int dim        = max(nx, max(ny, nz));
    const unsigned int numu = FLUID_NUM_EVOL_VARS;
// We need to subtract off the spatial components of the four-velocity.
#pragma message( \
    "Primitive variable count hard-coded. Works now, but could break in the future.")
    const unsigned int nump = FLUID_NUM_PRIM_VARS - 3;

    // We need left and right side grids. We overestimate the
    // size of the grid needed so we don't have to reallocate
    // as frequently.
    static double *F1[numu], *ul1[numu], *ur1[numu];
    static double *v1[nump], *vl1[nump], *vr1[nump];
    static unsigned int size = 2 * dim;

    static int first_call    = 1;
    if (first_call == 1) {
        for (int m = 0; m < numu; m++) {
            F1[m]  = new double[size];
            ul1[m] = new double[size];
            ur1[m] = new double[size];

            if (F1[m] == nullptr || ul1[m] == nullptr || ur1[m] == nullptr) {
                printf("hllerhs:    Out of memory. Die here.\n");
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
        for (int m = 0; m < nump; m++) {
            v1[m]  = new double[size];
            vl1[m] = new double[size];
            vr1[m] = new double[size];
            if (v1[m] == nullptr || vl1[m] == nullptr || vr1[m] == nullptr) {
                printf("hllerhs:    Out of memory. Die here.\n");
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
        first_call = 0;
    }

    // If the grid is too big to fit the memory we've allocated,
    // we need to allocate more.
    if (dim > size) {
        size = dim * 2;

        for (int m = 0; m < numu; m++) {
            delete[] F1[m];
            delete[] ul1[m];
            delete[] ur1[m];
            delete[] v1[m];
            delete[] vl1[m];
            delete[] vr1[m];

            F1[m]  = new double[size];
            ul1[m] = new double[size];
            ur1[m] = new double[size];
            v1[m]  = new double[size];
            vl1[m] = new double[size];
            vr1[m] = new double[size];

            if (F1[m] == nullptr || ul1[m] == nullptr || ur1[m] == nullptr) {
                printf("hllerhs:    Out of memory. Die here.\n");
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
        for (int m = 0; m < nump; m++) {
            v1[m]  = new double[size];
            vl1[m] = new double[size];
            vr1[m] = new double[size];
            if (v1[m] == nullptr || vl1[m] == nullptr || vr1[m] == nullptr) {
                printf("hllerhs:    Out of memory. Die here.\n");
                MPI_Abort(MPI_COMM_WORLD, 0);
            }
        }
    }

    // Calculate the flux in all three directions.
    cal_nflux(dtu, u, v, ul1, ur1, vl1, vr1, F1, v1, DIR_X, pmin, pmax, sz);
#ifdef DEBUG_FLUX
    if (!check_padding(dtu, sz)) {
        printf(
            "An error occurred in the padding points for the righthand side in "
            "DIR_X!");
    }
#endif
    cal_nflux(dtu, u, v, ul1, ur1, vl1, vr1, F1, v1, DIR_Y, pmin, pmax, sz);
#ifdef DEBUG_FLUX
    if (!check_padding(dtu, sz)) {
        printf(
            "An error occurred in the padding points for the righthand side in "
            "DIR_Y!");
    }
#endif
    cal_nflux(dtu, u, v, ul1, ur1, vl1, vr1, F1, v1, DIR_Z, pmin, pmax, sz);
#ifdef DEBUG_FLUX
    if (!check_padding(dtu, sz)) {
        printf(
            "An error occurred in the padding points for the righthand side in "
            "DIR_Z!");
    }
#endif

    // If we're in cylindrical coordinates, throw in a source term.
    if (FLUID_COORDS == 1) {
        cal_source(dtu, u, v, pmin, pmax, sz);
    }

    // Check that the right-hand side is physical. Based on the inequality
    //        S^2 < (tau + D)^2 - D^2,
    // we can derive a relationship
    //        |S|d|S|/dt < (tau + D)(dtau/dt + dD/dt) - DdD/dt
    // which must hold.
    // This needs to be generalized for alternative EOSes; assumes an ideal gas.
    // FIXME: hardcoded for 3 ghost points.
    /*const unsigned int nb = 3;
    for(unsigned int i = nb; i < sz[0] - nb; i++){
      for(unsigned int j = nb; j < sz[1] - nb; j++){
        for(unsigned int k = nb; k < sz[2] - nb; k++){
          unsigned int pp = i + sz[0]*(j + k*sz[1]);
          // Collect the conserved variables and their righthand sides.
          double D = u[fluid::VAR::U_D][pp];
          double Sd[3] =
    {u[fluid::VAR::U_SX][pp],u[fluid::VAR::U_SY][pp],u[fluid::VAR::U_SZ][pp]};
          double tau = u[fluid::VAR::U_TAU][pp];
          double dtD = dtu[fluid::VAR::U_D][pp];
          double dtSd[3] =
    {dtu[fluid::VAR::U_SX][pp],dtu[fluid::VAR::U_SY][pp],dtu[fluid::VAR::U_SZ][pp]};
          double dttau = dtu[fluid::VAR::U_TAU][pp];
          // Calculate the position so we can calculate the metric.
          double pos[3] = {0};
          double dx = (pmax[0] - pmin[0])/(sz[0] - 1);
          double dy = (pmax[1] - pmin[1])/(sz[1] - 1);
          double dz = (pmax[2] - pmin[2])/(sz[2] - 1);
          pos[0] = pmin[0] + i*dx;
          pos[1] = pmin[1] + j*dy;
          pos[2] = pmin[2] + k*dz;

          double gu[3][3], gd[3][3], detg;
          metric_vars(gu, &detg, gd, pos);

          // Get the absolute value of S as well as the up components of dtS and
    S.
          //double S = sqrt(square_form(Sd, gu));
          double Su[3];
          double dtSu[3];
          form_raise_index(Su, Sd, gu);
          form_raise_index(dtSu, dtSd, gu);
          double SdtS = 0.5*(dtSd[0]*Su[0] + dtSd[1]*Su[1] + dtSd[2]*Su[2] +
                        dtSu[0]*Sd[0] + dtSu[1]*Sd[1] + dtSu[2]*Sd[2]);
          double SdtSmax = dtD*tau + dttau*(D + tau);
          if(SdtSmax < SdtS){
            printf("rhshlle:  the righthand side is unphysical!\n");
            printf("rhshlle:  SdtS = %25.20e\n",SdtS);
            printf("rhshlle:  SdtSmax = %25.20e\n",SdtSmax);
          }
        }
      }
    }*/
}
// }}}

// cal_nflux {{{
void hlle::cal_nflux(double *dtu[], const double *u[], const double *v[],
                     double *ul1[], double *ur1[], double *vl1[], double *vr1[],
                     double *F1[], double *v1[], const unsigned int dir,
                     const double *pmin, const double *pmax,
                     const unsigned int *sz) {
    // We need to determine a loop ordering based on the direction
    // we're going in.
    // int perm[3];
    unsigned int ip, il, ix;  // direction index
    unsigned int np, nl, nx;  // Grid dimensions
    double dp, dl, dx;        // Grid increments
    double sp, sl, sx;        // Grid starting points
    double xp, xl, xx;

    double nbp = 3;  // Padding points.

    // An index calculating function.
    unsigned int (*index)(unsigned int, unsigned int, unsigned int,
                          unsigned int, unsigned int, unsigned int);
    // A position calculating function.
    void (*cal_pos)(unsigned int, unsigned int, unsigned int, double, double,
                    double, double, double, double, double *);
    void (*deriv)(double, double, double, double *);
    // Set up direction-dependent stuff {{{
    switch (dir) {
        case DIR_X:
            // Loop order should be a plane in z, a line in y, and a point in x.
            ip      = DIR_Z;
            il      = DIR_Y;
            ix      = DIR_X;

            // Define our indexing function.
            index   = &util::index_xm;
            // Calculate our position.
            cal_pos = &util::cal_pos_xm;
            // Define our derivative.
            deriv   = &deriv_xm;
            break;
        case DIR_Y:
            // Loop order should be a plane in x, a line in z, and a point in y.
            ip      = DIR_X;
            il      = DIR_Z;
            ix      = DIR_Y;

            // Define our indexing function.
            index   = &util::index_ym;
            // Calculate our position.
            cal_pos = &util::cal_pos_ym;
            // Define our derivative.
            deriv   = &deriv_yzm;
            break;
        case DIR_Z:
            // Loop order should be a plane in y, a line in x, and a point in z.
            ip      = DIR_Y;
            il      = DIR_X;
            ix      = DIR_Z;

            // Define our indexing function.
            index   = &util::index_zm;
            // Calculate our position.
            cal_pos = &util::cal_pos_zm;
            // Define our derivative.
            deriv   = &deriv_yzm;
            break;
        default:
            printf(
                "cal_nflux: dir=%d does not correspond to a valid direction.\n",
                dir);
            MPI_Abort(MPI_COMM_WORLD, 0);
            break;
    }
    // }}}

    // Get the grid dimensions
    np = sz[ip];
    nl = sz[il];
    nx = sz[ix];

    // Calculate the grid increments
    dp = (pmax[ip] - pmin[ip]) / (np - 1);
    dl = (pmax[il] - pmin[il]) / (nl - 1);
    dx = (pmax[ix] - pmin[ix]) / (nx - 1);

    // Get the grid starting points.
    sp = pmin[ip];
    sl = pmin[il];
    sx = pmin[ix];

    // Loop through the planes
    for (unsigned int k = nbp; k < np - nbp; k++) {
        // Calculate the cell-centered coordinate of the plane.
        xp = sp + k * dp;

        for (unsigned int j = nbp; j < nl - nbp; j++) {
            // Calculate the cell-centered coordinate of the line.
            xl = sl + j * dl;

            // Extract constant lines, reconstruct, and calculate
            // the numerical flux.
            for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
                for (unsigned int i = 0; i < nx; i++) {
                    v1[m][i] = v[m][index(i, j, k, np, nl, nx)];
                }
                // I think this is not needed. (the behaviour is different with
                // this and without this, hence I reverted it to the previous
                // version. @jacob can you please make sure this is correct. )
                // v1[m][0] = v1[m][1];
                // v1[m][nx-1] = v1[m][nx-2];
            }

            double pos[3];
            cal_pos(0, j, k, dx, dl, dp, sx, sl, sp, pos);
            recon::reconstructvars(vl1, vr1, v1, F1, dir, pos, dx, sz);

            // Calculate the conserved variables at the cell interfaces
            for (unsigned int i = 0; i < nx; i++) {
                double ulpt[FLUID_NUM_EVOL_VARS], urpt[FLUID_NUM_EVOL_VARS],
                    vlpt[FLUID_NUM_EVOL_VARS], vrpt[FLUID_NUM_EVOL_VARS];
                cal_pos(i, j, k, dx, dl, dp, sx, sl, sp, pos);

                for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
                    vlpt[m] = vl1[m][i];
                    vrpt[m] = vr1[m][i];
                }

                double gd[3][3], detg;
                metric(gd, &detg, pos);
                double sdetg = sqrt(detg);

                fluid_math::prim_to_con_pt(vlpt, ulpt, gd, sdetg);
                fluid_math::prim_to_con_pt(vrpt, urpt, gd, sdetg);

                for (int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
                    ul1[m][i] = ulpt[m];
                    ur1[m][i] = urpt[m];
                }
            }

            // Calculate the flux.
            pos[dir] = pmin[dir];
#ifdef LF_FLUX
            lf_flux(F1, ul1, ur1, vl1, vr1, dir, pos, dx, sz);
#else
            hlle_flux(F1, ul1, ur1, vl1, vr1, dir, pos, dx, sz);
#endif

            // Calculate the spatial derivative
            for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
                for (unsigned int i = nbp; i < nx - nbp; i++) {
                    unsigned int pp = index(i, j, k, np, nl, nx);
                    double bd       = dtu[m][pp];
                    deriv(F1[m][i + 1], F1[m][i], dx, &bd);
                    dtu[m][pp] = bd;
                }
            }
        }
    }
}
// }}}

// cal_source {{{
void hlle::cal_source(double *dtu[], const double *u[], const double *v[],
                      const double *pmin, const double *pmax,
                      const unsigned int *sz) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

#pragma message("FIXME: ghost point width hardcoded to 3.")
    const int ghosts = 3;
    const int ai     = ghosts - 1;

    for (unsigned int k = 0; k < nz; k++) {
        for (unsigned int j = 0; j < ny; j++) {
            for (unsigned int i = 0; i < nx; i++) {
                unsigned int pp = i + nx * (j + ny * k);

                // Modify the momentum.
                dtu[VAR::U_SX][pp] += v[PVAR::V_P][pp];
            }
        }
    }
}
// }}}

// lf_flux {{{
void hlle::lf_flux(double *F[], double *ul1[], double *ur1[], double *vl1[],
                   double *vr1[], const int dir, const double pos_in[],
                   const double dx, const unsigned int *sz) {
    int n;
    double xm;
    double pos[3];
    void (*flux)(double[], double[], double[], double, double[3][3],
                 double[3][3]);
    switch (dir) {
        case 0:
            n    = sz[0];
            xm   = pos_in[0];
            flux = &flux_x;
            break;
        case 1:
            n    = sz[1];
            xm   = pos_in[1];
            flux = &flux_y;
            break;
        case 2:
            n    = sz[2];
            xm   = pos_in[2];
            flux = &flux_z;
            break;
    }

    pos[0] = pos_in[0];
    pos[1] = pos_in[1];
    pos[2] = pos_in[2];

    for (int i = 0; i < n; i++) {
        pos[dir] = xm + (i - 0.5) * dx;

        double ul[FLUID_NUM_EVOL_VARS], ur[FLUID_NUM_EVOL_VARS],
            fl[FLUID_NUM_EVOL_VARS], fr[FLUID_NUM_EVOL_VARS],
            vl[FLUID_NUM_EVOL_VARS], vr[FLUID_NUM_EVOL_VARS];
        double gd[3][3], gu[3][3], detg;

        metric_vars(gu, &detg, gd, pos);
        double sdetg = sqrt(detg);

        for (int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            ul[m] = ul1[m][i];
            ur[m] = ur1[m][i];
        }
        for (int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            vl[m] = vl1[m][i];
            vr[m] = vr1[m][i];
        }

        flux(fl, ul, vl, sdetg, gd, gu);
        flux(fr, ur, vr, sdetg, gd, gu);

        double bp, bm;
        // int rc = cal_max_speeds(&bp, &bm, vl, vr, gd, gu, dir);
        // double alpha = fmax(fabs(bp), fabs(bm));
        double alpha = 1.0;
        for (unsigned int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            // Not quite correct, but it should be stable enough at least to
            // let us check for bugs.
            F[m][i] = 0.5 * ((fl[m] + fr[m]) - alpha * (ur[m] - ul[m]));
        }
    }
}
// }}}

// hlle_flux {{{
void hlle::hlle_flux(double *F[], double *ul1[], double *ur1[], double *vl1[],
                     double *vr1[], const int dir, const double pos_in[],
                     const double dx, const unsigned int *sz) {
    int n;
    double xm;
    double pos[3];
    void (*flux)(double[], double[], double[], double, double[3][3],
                 double[3][3]);
    switch (dir) {
        case 0:
            n    = sz[0];
            xm   = pos_in[0];
            flux = &flux_x;
            break;
        case 1:
            n    = sz[1];
            xm   = pos_in[1];
            flux = &flux_y;
            break;
        case 2:
            n    = sz[2];
            xm   = pos_in[2];
            flux = &flux_z;
            break;
    }

    pos[0] = pos_in[0];
    pos[1] = pos_in[1];
    pos[2] = pos_in[2];

    for (int i = 0; i < n; i++) {
        pos[dir] = xm + (i - 0.5) * dx;

        double ul[FLUID_NUM_EVOL_VARS], ur[FLUID_NUM_EVOL_VARS],
            fl[FLUID_NUM_EVOL_VARS], fr[FLUID_NUM_EVOL_VARS],
            vl[FLUID_NUM_EVOL_VARS], vr[FLUID_NUM_EVOL_VARS];
        double gd[3][3], gu[3][3], detg;

        metric_vars(gu, &detg, gd, pos);
        double sdetg = sqrt(detg);

        for (int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            ul[m] = ul1[m][i];
            ur[m] = ur1[m][i];
        }
        for (int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            vl[m] = vl1[m][i];
            vr[m] = vr1[m][i];
        }

        flux(fl, ul, vl, sdetg, gd, gu);
        flux(fr, ur, vr, sdetg, gd, gu);

        double bp, bm;

        int rc = cal_max_speeds(&bp, &bm, vl, vr, gd, gu, dir);

        if (rc < 0) default_speeds(&bp, &bm, gu, dir);

        for (int m = 0; m < FLUID_NUM_EVOL_VARS; m++) {
            F[m][i] = ((bp * fl[m] - bm * fr[m]) + bp * bm * (ur[m] - ul[m])) /
                      (bp - bm);
        }
    }
}
// }}}

// flux_x {{{
void hlle::flux_x(double f[], double u[], double v[], double sdetg,
                  double gd[3][3], double gu[3][3]) {
    const double Alpha   = 1.0;
    const double Beta[3] = {0.0, 0.0, 0.0};

    // if( sdetg < 1.0e-7){
    if (sdetg < 1.0e-15) {
        f[VAR::U_D]   = 0.0;
        f[VAR::U_SX]  = 0.0;
        f[VAR::U_SY]  = 0.0;
        f[VAR::U_SZ]  = 0.0;
        f[VAR::U_TAU] = 0.0;
        return;
    }

    double D, Tau, P;
    double Sd[3], Su[3], vd[3], vu[3];
    // double isdetg = 1.0 / (sdetg + 1.0e-15);

    /*D     = isdetg * u[VAR::U_D  ];
    Sd[0] = isdetg * u[VAR::U_SX ];
    Sd[1] = isdetg * u[VAR::U_SY ];
    Sd[2] = isdetg * u[VAR::U_SZ ];
    Tau   = isdetg * u[VAR::U_TAU];*/
    D     = u[VAR::U_D];
    Sd[0] = u[VAR::U_SX];
    Sd[1] = u[VAR::U_SY];
    Sd[2] = u[VAR::U_SZ];
    Tau   = u[VAR::U_TAU];

    vu[0] = v[PVAR::V_VX];
    vu[1] = v[PVAR::V_VY];
    vu[2] = v[PVAR::V_VZ];
    P     = v[PVAR::V_P];

    // double vsq = square_vector(vu, gd);
    vector_lower_index(vd, vu, gd);
    form_raise_index(Su, Sd, gu);

    // double Wsq = 1.0/(1.0 - vsq);

    /*f[VAR::U_D  ] = Alpha*sdetg*D*(vu[0] - Beta[0]/Alpha);
    f[VAR::U_SX ] = Alpha*sdetg*(Sd[0]*vu[0] + P - Beta[0]/Alpha*Sd[0]);
    f[VAR::U_SY ] = Alpha*sdetg*(Sd[1]*vu[0] - Beta[0]/Alpha*Sd[1]);
    f[VAR::U_SZ ] = Alpha*sdetg*(Sd[2]*vu[0] - Beta[0]/Alpha*Sd[2]);
    f[VAR::U_TAU] = Alpha*sdetg*(Su[0] - D*vu[0]) - sdetg*Beta[0]*Tau;*/
    f[VAR::U_D] = Alpha * D * (vu[0] - Beta[0] / Alpha);
    f[VAR::U_SX] =
        Alpha * (Sd[0] * vu[0] + sdetg * P - Beta[0] / Alpha * Sd[0]);
    f[VAR::U_SY]  = Alpha * (Sd[1] * vu[0] - Beta[0] / Alpha * Sd[1]);
    f[VAR::U_SZ]  = Alpha * (Sd[2] * vu[0] - Beta[0] / Alpha * Sd[2]);
    f[VAR::U_TAU] = Alpha * (Su[0] - D * vu[0]) - Beta[0] * Tau;
}
// }}}

// flux_y {{{
void hlle::flux_y(double f[], double u[], double v[], double sdetg,
                  double gd[3][3], double gu[3][3]) {
    const double Alpha   = 1.0;
    const double Beta[3] = {0.0, 0.0, 0.0};

    // if( sdetg < 1.0e-7){
    if (sdetg < 1.0e-15) {
        f[VAR::U_D]   = 0.0;
        f[VAR::U_SX]  = 0.0;
        f[VAR::U_SY]  = 0.0;
        f[VAR::U_SZ]  = 0.0;
        f[VAR::U_TAU] = 0.0;
        return;
    }

    double D, Tau, P;
    double Sd[3], Su[3], vd[3], vu[3];
    // double isdetg = 1.0 / (sdetg + 1.0e-15);

    /*D     = isdetg * u[VAR::U_D  ];
    Sd[0] = isdetg * u[VAR::U_SX ];
    Sd[1] = isdetg * u[VAR::U_SY ];
    Sd[2] = isdetg * u[VAR::U_SZ ];
    Tau   = isdetg * u[VAR::U_TAU];*/
    D     = u[VAR::U_D];
    Sd[0] = u[VAR::U_SX];
    Sd[1] = u[VAR::U_SY];
    Sd[2] = u[VAR::U_SZ];
    Tau   = u[VAR::U_TAU];

    vu[0] = v[PVAR::V_VX];
    vu[1] = v[PVAR::V_VY];
    vu[2] = v[PVAR::V_VZ];
    P     = v[PVAR::V_P];

    // double vsq = square_vector(vu, gd);
    vector_lower_index(vd, vu, gd);
    form_raise_index(Su, Sd, gu);

    // double Wsq = 1.0/(1.0 - vsq);

    /*f[VAR::U_D  ] = Alpha*sdetg*D*(vu[1] - Beta[1]/Alpha);
    f[VAR::U_SX ] = Alpha*sdetg*(Sd[0]*vu[1] - Beta[1]/Alpha*Sd[0]);
    f[VAR::U_SY ] = Alpha*sdetg*(Sd[1]*vu[1] + P - Beta[1]/Alpha*Sd[1]);
    f[VAR::U_SZ ] = Alpha*sdetg*(Sd[2]*vu[1] - Beta[1]/Alpha*Sd[2]);
    f[VAR::U_TAU] = Alpha*sdetg*(Su[1] - D*vu[1]) - sdetg*Beta[1]*Tau;*/
    f[VAR::U_D]  = Alpha * D * (vu[1] - Beta[1] / Alpha);
    f[VAR::U_SX] = Alpha * (Sd[0] * vu[1] - Beta[1] / Alpha * Sd[0]);
    f[VAR::U_SY] =
        Alpha * (Sd[1] * vu[1] + sdetg * P - Beta[1] / Alpha * Sd[1]);
    f[VAR::U_SZ]  = Alpha * (Sd[2] * vu[1] - Beta[1] / Alpha * Sd[2]);
    f[VAR::U_TAU] = Alpha * (Su[1] - D * vu[1]) - Beta[1] * Tau;
}
// }}}

// flux_z {{{
void hlle::flux_z(double f[], double u[], double v[], double sdetg,
                  double gd[3][3], double gu[3][3]) {
    const double Alpha   = 1.0;
    const double Beta[3] = {0.0, 0.0, 0.0};

    // if( sdetg < 1.0e-7){
    if (sdetg < 1.0e-15) {
        f[VAR::U_D]   = 0.0;
        f[VAR::U_SX]  = 0.0;
        f[VAR::U_SY]  = 0.0;
        f[VAR::U_SZ]  = 0.0;
        f[VAR::U_TAU] = 0.0;
        return;
    }

    double D, Tau, P;
    double Sd[3], Su[3], vd[3], vu[3];
    // double isdetg = 1.0 / (sdetg + 1.0e-15);

    /*D     = isdetg * u[VAR::U_D  ];
    Sd[0] = isdetg * u[VAR::U_SX ];
    Sd[1] = isdetg * u[VAR::U_SY ];
    Sd[2] = isdetg * u[VAR::U_SZ ];
    Tau   = isdetg * u[VAR::U_TAU];*/
    D     = u[VAR::U_D];
    Sd[0] = u[VAR::U_SX];
    Sd[1] = u[VAR::U_SY];
    Sd[2] = u[VAR::U_SZ];
    Tau   = u[VAR::U_TAU];

    vu[0] = v[PVAR::V_VX];
    vu[1] = v[PVAR::V_VY];
    vu[2] = v[PVAR::V_VZ];
    P     = v[PVAR::V_P];

    // double vsq = square_vector(vu, gd);
    vector_lower_index(vd, vu, gd);
    form_raise_index(Su, Sd, gu);

    // double Wsq = 1.0/(1.0 - vsq);

    /*f[VAR::U_D  ] = Alpha*sdetg*D*(vu[2] - Beta[2]/Alpha);
    f[VAR::U_SX ] = Alpha*sdetg*(Sd[0]*vu[2] - Beta[2]/Alpha*Sd[0]);
    f[VAR::U_SY ] = Alpha*sdetg*(Sd[1]*vu[2] - Beta[2]/Alpha*Sd[1]);
    f[VAR::U_SZ ] = Alpha*sdetg*(Sd[2]*vu[2] + P - Beta[2]/Alpha*Sd[2]);
    f[VAR::U_TAU] = Alpha*sdetg*(Su[2] - D*vu[2]) - sdetg*Beta[2]*Tau;*/
    f[VAR::U_D]  = Alpha * D * (vu[2] - Beta[2] / Alpha);
    f[VAR::U_SX] = Alpha * (Sd[0] * vu[2] - Beta[2] / Alpha * Sd[0]);
    f[VAR::U_SY] = Alpha * (Sd[1] * vu[2] - Beta[2] / Alpha * Sd[1]);
    f[VAR::U_SZ] =
        Alpha * (Sd[2] * vu[2] + sdetg * P - Beta[2] / Alpha * Sd[2]);
    f[VAR::U_TAU] = Alpha * (Su[2] - D * vu[2]) - Beta[2] * Tau;
}
// }}}

// cal_max_speeds {{{
int hlle::cal_max_speeds(double *bp, double *bm, double vl[], double vr[],
                         double gd[3][3], double gu[3][3], const int dir) {
    const double Alpha   = 1.0;
    const double Beta[3] = {0.0, 0.0, 0.0};
    const double gamma   = FLUID_GAMMA;

    double vlu[3], vru[3];

    double rhol  = vl[PVAR::V_RHO];
    double Pl    = vl[PVAR::V_P];
    vlu[0]       = vl[PVAR::V_VX];
    vlu[1]       = vl[PVAR::V_VY];
    vlu[2]       = vl[PVAR::V_VZ];

    double rhor  = vr[PVAR::V_RHO];
    double Pr    = vr[PVAR::V_P];
    vru[0]       = vr[PVAR::V_VX];
    vru[1]       = vr[PVAR::V_VY];
    vru[2]       = vr[PVAR::V_VZ];

    double csl   = cal_cs(rhol, Pl, gamma);
    double csr   = cal_cs(rhor, Pr, gamma);

    double cslsq = csl * csl;
    double csrsq = csr * csr;

    double hkk   = gu[dir][dir];
    double betak = Beta[dir];
    double vlk   = vlu[dir];
    double vrk   = vru[dir];

    double vlsq  = square_vector(vlu, gd);
    double vrsq  = square_vector(vru, gd);

    assert(vlsq < 1.0);
    assert(vrsq < 1.0);

    double dis;
    dis = cslsq * (1.0 - vlsq) *
          (hkk * (1.0 - vlsq * cslsq) - vlk * vlk * (1.0 - cslsq));
    if (dis < 0.0) {
        return -1;
    }
    double Lpl =
        Alpha / (1.0 - vlsq * cslsq) * (vlk * (1.0 - cslsq) + sqrt(dis)) -
        betak;
    double Lml =
        Alpha / (1.0 - vlsq * cslsq) * (vlk * (1.0 - cslsq) - sqrt(dis)) -
        betak;

    dis = csrsq * (1.0 - vrsq) *
          (hkk * (1.0 - vrsq * csrsq) - vrk * vrk * (1.0 - csrsq));
    if (dis < 0.0) {
        return -1;
    }
    double Lpr =
        Alpha / (1.0 - vrsq * csrsq) * (vrk * (1.0 - csrsq) + sqrt(dis)) -
        betak;
    double Lmr =
        Alpha / (1.0 - vrsq * csrsq) * (vrk * (1.0 - csrsq) - sqrt(dis)) -
        betak;

    double t1 = fmax(0.0, Lpl);
    *bp       = fmax(t1, Lpr);

    t1        = fmin(0.0, Lml);
    *bm       = fmin(t1, Lmr);

    return 1;
}
// }}}

// default_speeds {{{
void hlle::default_speeds(double *bp, double *bm, double gu[3][3],
                          const int dir) {
    const double Alpha   = 1.0;
    const double Beta[3] = {0.0, 0.0, 0.0};

    *bp                  = Alpha * sqrt(gu[dir][dir]) - Beta[dir];
    *bm                  = -Alpha * sqrt(gu[dir][dir]) - Beta[dir];
}
// }}}

// cal_cs {{{
double hlle::cal_cs(double rho, double p, double gamma) {
    const double factor = 1.0e3;
    double y, cs;

    if (rho * factor < p) {
        y  = (gamma - 1.0) * rho / (gamma * p);
        cs = sqrt(gamma - 1.0) *
             (1.0 + (-0.5 + (0.375 + (-0.3125 + 0.2734375 * y) * y) * y) * y);
    } else {
        cs =
            sqrt(gamma * (gamma - 1.0) * p / ((gamma - 1.0) * rho + gamma * p));
    }

    return cs;
}
// }}}

// check_padding {{{
bool hlle::check_padding(double *dtu[], const unsigned int *sz) {
    unsigned int nb = 3;
    for (unsigned int m = 0; m < 5; m++) {
        for (unsigned int k = 0; k < sz[2]; k++) {
            for (unsigned int j = 0; j < sz[1]; j++) {
                for (unsigned int i = 0; i < sz[0]; i++) {
                    unsigned int pp = i + sz[0] * (j + sz[1] * k);
                    if (i >= nb && i < sz[0] - nb && j >= nb &&
                        j < sz[1] - nb && k >= nb && k < sz[2] - nb) {
                        continue;
                    }
                    if (dtu[m][pp] != 0.0) {
                        printf("Padding error: dtu[%d][%d,%d,%d] = %g != 0!\n",
                               m, i, j, k, dtu[m][pp]);
                        return false;
                    }
                }
            }
        }

        return true;
    }
}

// }}}
