#include "rhs.h"

using namespace std;
using namespace em1;

/*----------------------------------------------------------------------;
 *
 * RHS for non-linear sigma model
 *
 *----------------------------------------------------------------------*/
void em1rhs(double **unzipVarsRHS, const double **uZipVars,
            const unsigned int &offset, const double *pmin, const double *pmax,
            const unsigned int *sz, const unsigned int &bflag) {
    const double *psi     = &uZipVars[VAR::U_PSI][offset];
    const double *E0      = &uZipVars[VAR::U_E0][offset];
    const double *E1      = &uZipVars[VAR::U_E1][offset];
    const double *E2      = &uZipVars[VAR::U_E2][offset];
    const double *A0      = &uZipVars[VAR::U_A0][offset];
    const double *A1      = &uZipVars[VAR::U_A1][offset];
    const double *A2      = &uZipVars[VAR::U_A2][offset];

    double *psi_rhs       = &unzipVarsRHS[VAR::U_PSI][offset];
    double *E_rhs0        = &unzipVarsRHS[VAR::U_E0][offset];
    double *E_rhs1        = &unzipVarsRHS[VAR::U_E1][offset];
    double *E_rhs2        = &unzipVarsRHS[VAR::U_E2][offset];
    double *A_rhs0        = &unzipVarsRHS[VAR::U_A0][offset];
    double *A_rhs1        = &unzipVarsRHS[VAR::U_A1][offset];
    double *A_rhs2        = &unzipVarsRHS[VAR::U_A2][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    mem::memory_pool<DendroScalar> *__mem_pool = &EM1_MEM_POOL;

    double PI                                  = 3.14159265358979323846;

    double hx                                  = (pmax[0] - pmin[0]) / (nx - 1);
    double hy                                  = (pmax[1] - pmin[1]) / (ny - 1);
    double hz                                  = (pmax[2] - pmin[2]) / (nz - 1);

    int idx[3];

    unsigned int n = sz[0] * sz[1] * sz[2];

    em1::timer::t_deriv.start();

    double *rho_e = __mem_pool->allocate(n);

    double *J0    = __mem_pool->allocate(n);
    double *J1    = __mem_pool->allocate(n);
    double *J2    = __mem_pool->allocate(n);

    for (unsigned int m = 0; m < n; m++) {
        rho_e[m] = 0.0;

        J0[m]    = 0.0;
        J1[m]    = 0.0;
        J2[m]    = 0.0;
    }

    double *grad_0_E0    = __mem_pool->allocate(n);
    double *grad_1_E0    = __mem_pool->allocate(n);
    double *grad_2_E0    = __mem_pool->allocate(n);

    double *grad_0_E1    = __mem_pool->allocate(n);
    double *grad_1_E1    = __mem_pool->allocate(n);
    double *grad_2_E1    = __mem_pool->allocate(n);

    double *grad_0_E2    = __mem_pool->allocate(n);
    double *grad_1_E2    = __mem_pool->allocate(n);
    double *grad_2_E2    = __mem_pool->allocate(n);

    double *grad_0_A0    = __mem_pool->allocate(n);
    double *grad_1_A0    = __mem_pool->allocate(n);
    double *grad_2_A0    = __mem_pool->allocate(n);

    double *grad_0_A1    = __mem_pool->allocate(n);
    double *grad_1_A1    = __mem_pool->allocate(n);
    double *grad_2_A1    = __mem_pool->allocate(n);

    double *grad_0_A2    = __mem_pool->allocate(n);
    double *grad_1_A2    = __mem_pool->allocate(n);
    double *grad_2_A2    = __mem_pool->allocate(n);

    double *grad_0_psi   = __mem_pool->allocate(n);
    double *grad_1_psi   = __mem_pool->allocate(n);
    double *grad_2_psi   = __mem_pool->allocate(n);

    double *grad2_0_1_A0 = __mem_pool->allocate(n);
    double *grad2_0_2_A0 = __mem_pool->allocate(n);
    double *grad2_1_1_A0 = __mem_pool->allocate(n);
    double *grad2_2_2_A0 = __mem_pool->allocate(n);

    double *grad2_0_0_A1 = __mem_pool->allocate(n);
    double *grad2_0_1_A1 = __mem_pool->allocate(n);
    double *grad2_1_2_A1 = __mem_pool->allocate(n);
    double *grad2_2_2_A1 = __mem_pool->allocate(n);

    double *grad2_0_0_A2 = __mem_pool->allocate(n);
    double *grad2_0_2_A2 = __mem_pool->allocate(n);
    double *grad2_1_1_A2 = __mem_pool->allocate(n);
    double *grad2_1_2_A2 = __mem_pool->allocate(n);

    deriv_x(grad_0_psi, psi, hx, sz, bflag);
    deriv_y(grad_1_psi, psi, hy, sz, bflag);
    deriv_z(grad_2_psi, psi, hz, sz, bflag);

    // derivs for A0.
    deriv_x(grad_0_A0, A0, hx, sz, bflag);            // dx  (A0)
    deriv_y(grad2_0_1_A0, grad_0_A0, hy, sz, bflag);  // dxy (A0)
    deriv_z(grad2_0_2_A0, grad_0_A0, hz, sz, bflag);  // dxz (A0)
    deriv_yy(grad2_1_1_A0, A0, hy, sz, bflag);
    deriv_zz(grad2_2_2_A0, A0, hz, sz, bflag);

    // derivs for A1
    deriv_x(grad_0_A1, A1, hx, sz, bflag);            // dx(A1)
    deriv_y(grad2_0_1_A1, grad_0_A1, hy, sz, bflag);  // dxy(A1)

    deriv_y(grad_1_A1, A1, hy, sz, bflag);            // dy(A1)
    deriv_z(grad2_1_2_A1, grad_1_A1, hz, sz, bflag);  // dyz(A1)
    deriv_xx(grad2_0_0_A1, A1, hx, sz, bflag);
    deriv_zz(grad2_2_2_A1, A1, hz, sz, bflag);

    // derivs for A2
    deriv_x(grad_0_A2, A2, hx, sz, bflag);            // dx(A2)
    deriv_z(grad2_0_2_A2, grad_0_A2, hz, sz, bflag);  // dxz(A2)

    deriv_y(grad_1_A2, A2, hy, sz, bflag);            // dy(A2)
    deriv_z(grad2_1_2_A2, grad_1_A2, hz, sz, bflag);  // dyz(A2)

    deriv_xx(grad2_0_0_A2, A2, hx, sz, bflag);
    deriv_yy(grad2_1_1_A2, A2, hy, sz, bflag);

    deriv_z(grad_2_A2, A2, hz, sz, bflag);

    em1::timer::t_deriv.stop();

    double r;
    double eta;
    const unsigned int PW = em1::EM1_PADDING_WIDTH;
    // cout << "begin loop" << endl;
    for (unsigned int k = PW; k < nz - PW; k++) {
        for (unsigned int j = PW; j < ny - PW; j++) {
            for (unsigned int i = PW; i < nx - PW; i++) {
                const double x        = pmin[0] + i * hx;
                const double y        = pmin[1] + j * hy;
                const double z        = pmin[2] + k * hz;
                const unsigned int pp = i + nx * (j + ny * k);
                // r= sqrt(x*x + y*y + z*z);

                em1::timer::t_rhs.start();
#include "em1_eqs.cpp"
                em1::timer::t_rhs.stop();
            }
        }
    }

    if (bflag != 0) {
        em1::timer::t_bdyc.start();

        // Some of these are redundant ...
        // deriv_x(grad_0_psi, psi, hx, sz, bflag);
        // deriv_y(grad_1_psi, psi, hy, sz, bflag);
        // deriv_z(grad_2_psi, psi, hz, sz, bflag);

        deriv_x(grad_0_E0, E0, hx, sz, bflag);
        deriv_y(grad_1_E0, E0, hy, sz, bflag);
        deriv_z(grad_2_E0, E0, hz, sz, bflag);

        deriv_x(grad_0_E1, E1, hx, sz, bflag);
        deriv_y(grad_1_E1, E1, hy, sz, bflag);
        deriv_z(grad_2_E1, E1, hz, sz, bflag);

        deriv_x(grad_0_E2, E2, hx, sz, bflag);
        deriv_y(grad_1_E2, E2, hy, sz, bflag);
        deriv_z(grad_2_E2, E2, hz, sz, bflag);

        // deriv_x(grad_0_A0, A0, hx, sz, bflag);
        deriv_y(grad_1_A0, A0, hy, sz, bflag);
        deriv_z(grad_2_A0, A0, hz, sz, bflag);

        // deriv_x(grad_0_A1, A1, hx, sz, bflag);
        // deriv_y(grad_1_A1, A1, hy, sz, bflag);
        deriv_z(grad_2_A1, A1, hz, sz, bflag);

        // deriv_x(grad_0_A2, A2, hx, sz, bflag);
        // deriv_y(grad_1_A2, A2, hy, sz, bflag);
        // deriv_z(grad_2_A2, A2, hz, sz, bflag);

        em1_bcs(psi_rhs, psi, grad_0_psi, grad_1_psi, grad_2_psi, pmin, pmax,
                1.0, 0.0, sz, bflag);

        em1_bcs(E_rhs0, E0, grad_0_E0, grad_1_E0, grad_2_E0, pmin, pmax, 2.0,
                0.0, sz, bflag);
        em1_bcs(E_rhs1, E1, grad_0_E1, grad_1_E1, grad_2_E1, pmin, pmax, 2.0,
                0.0, sz, bflag);
        em1_bcs(E_rhs2, E2, grad_0_E2, grad_1_E2, grad_2_E2, pmin, pmax, 2.0,
                0.0, sz, bflag);

        em1_bcs(A_rhs0, A0, grad_0_A0, grad_1_A0, grad_2_A0, pmin, pmax, 1.0,
                0.0, sz, bflag);
        em1_bcs(A_rhs1, A1, grad_0_A1, grad_1_A1, grad_2_A1, pmin, pmax, 1.0,
                0.0, sz, bflag);
        em1_bcs(A_rhs2, A2, grad_0_A2, grad_1_A2, grad_2_A2, pmin, pmax, 1.0,
                0.0, sz, bflag);

        em1::timer::t_bdyc.stop();
    }

    em1::timer::t_deriv.start();

    ko_deriv_x(grad_0_psi, psi, hx, sz, bflag);
    ko_deriv_y(grad_1_psi, psi, hy, sz, bflag);
    ko_deriv_z(grad_2_psi, psi, hz, sz, bflag);

    ko_deriv_x(grad_0_E0, E0, hx, sz, bflag);
    ko_deriv_y(grad_1_E0, E0, hy, sz, bflag);
    ko_deriv_z(grad_2_E0, E0, hz, sz, bflag);

    ko_deriv_x(grad_0_E1, E1, hx, sz, bflag);
    ko_deriv_y(grad_1_E1, E1, hy, sz, bflag);
    ko_deriv_z(grad_2_E1, E1, hz, sz, bflag);

    ko_deriv_x(grad_0_E2, E2, hx, sz, bflag);
    ko_deriv_y(grad_1_E2, E2, hy, sz, bflag);
    ko_deriv_z(grad_2_E2, E2, hz, sz, bflag);

    ko_deriv_x(grad_0_A0, A0, hx, sz, bflag);
    ko_deriv_y(grad_1_A0, A0, hy, sz, bflag);
    ko_deriv_z(grad_2_A0, A0, hz, sz, bflag);

    ko_deriv_x(grad_0_A1, A1, hx, sz, bflag);
    ko_deriv_y(grad_1_A1, A1, hy, sz, bflag);
    ko_deriv_z(grad_2_A1, A1, hz, sz, bflag);

    ko_deriv_x(grad_0_A2, A2, hx, sz, bflag);
    ko_deriv_y(grad_1_A2, A2, hy, sz, bflag);
    ko_deriv_z(grad_2_A2, A2, hz, sz, bflag);

    em1::timer::t_deriv.stop();

    em1::timer::t_rhs.start();

    const double sigma = KO_DISS_SIGMA;

    for (unsigned int k = PW; k < nz - PW; k++) {
        for (unsigned int j = PW; j < ny - PW; j++) {
            for (unsigned int i = PW; i < nx - PW; i++) {
                const unsigned int pp = i + nx * (j + ny * k);

                psi_rhs[pp] +=
                    sigma * (grad_0_psi[pp] + grad_1_psi[pp] + grad_2_psi[pp]);

                E_rhs0[pp] +=
                    sigma * (grad_0_E0[pp] + grad_1_E0[pp] + grad_2_E0[pp]);
                E_rhs1[pp] +=
                    sigma * (grad_0_E1[pp] + grad_1_E1[pp] + grad_2_E1[pp]);
                E_rhs2[pp] +=
                    sigma * (grad_0_E2[pp] + grad_1_E2[pp] + grad_2_E2[pp]);

                A_rhs0[pp] +=
                    sigma * (grad_0_A0[pp] + grad_1_A0[pp] + grad_2_A0[pp]);
                A_rhs1[pp] +=
                    sigma * (grad_0_A1[pp] + grad_1_A1[pp] + grad_2_A1[pp]);
                A_rhs2[pp] +=
                    sigma * (grad_0_A2[pp] + grad_1_A2[pp] + grad_2_A2[pp]);
            }
        }
    }

    em1::timer::t_rhs.stop();

    em1::timer::t_deriv.start();

    __mem_pool->free(grad2_0_1_A0);
    __mem_pool->free(grad2_0_2_A0);
    __mem_pool->free(grad2_1_1_A0);
    __mem_pool->free(grad2_2_2_A0);

    __mem_pool->free(grad2_0_0_A1);
    __mem_pool->free(grad2_0_1_A1);
    __mem_pool->free(grad2_1_2_A1);
    __mem_pool->free(grad2_2_2_A1);

    __mem_pool->free(grad2_0_0_A2);
    __mem_pool->free(grad2_0_2_A2);
    __mem_pool->free(grad2_1_1_A2);
    __mem_pool->free(grad2_1_2_A2);

    __mem_pool->free(grad_0_psi);
    __mem_pool->free(grad_1_psi);
    __mem_pool->free(grad_2_psi);

    __mem_pool->free(grad_0_E0);
    __mem_pool->free(grad_1_E0);
    __mem_pool->free(grad_2_E0);

    __mem_pool->free(grad_0_E1);
    __mem_pool->free(grad_1_E1);
    __mem_pool->free(grad_2_E1);

    __mem_pool->free(grad_0_E2);
    __mem_pool->free(grad_1_E2);
    __mem_pool->free(grad_2_E2);

    __mem_pool->free(grad_0_A0);
    __mem_pool->free(grad_1_A0);
    __mem_pool->free(grad_2_A0);

    __mem_pool->free(grad_0_A1);
    __mem_pool->free(grad_1_A1);
    __mem_pool->free(grad_2_A1);

    __mem_pool->free(grad_0_A2);
    __mem_pool->free(grad_1_A2);
    __mem_pool->free(grad_2_A2);

    __mem_pool->free(J0);
    __mem_pool->free(J1);
    __mem_pool->free(J2);

    __mem_pool->free(rho_e);

    em1::timer::t_deriv.stop();

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
void em1_bcs(double *f_rhs, const double *f, const double *dxf,
             const double *dyf, const double *dzf, const double *pmin,
             const double *pmax, const double f_falloff,
             const double f_asymptotic, const unsigned int *sz,
             const unsigned int &bflag) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx             = (pmax[0] - pmin[0]) / (nx - 1);
    double hy             = (pmax[1] - pmin[1]) / (ny - 1);
    double hz             = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int PW = em1::EM1_PADDING_WIDTH;
    unsigned int ib       = PW;
    unsigned int jb       = PW;
    unsigned int kb       = PW;
    unsigned int ie       = sz[0] - PW;
    unsigned int je       = sz[1] - PW;
    unsigned int ke       = sz[2] - PW;

    double x, y, z;
    unsigned int pp;
    double inv_r;

    if (bflag & (1u << OCT_DIR_LEFT)) {
        double x = pmin[0] + ib * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y     = pmin[1] + j * hy;
                pp    = IDX(ib, j, k);
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);

#ifdef EM1_DIRICHLET_BDY
                f_rhs[pp] = 0.0;
#else
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
#endif
            }
        }
    }

    if (bflag & (1u << OCT_DIR_RIGHT)) {
        x = pmin[0] + (ie - 1) * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y     = pmin[1] + j * hy;
                pp    = IDX((ie - 1), j, k);
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);

#ifdef EM1_DIRICHLET_BDY
                f_rhs[pp] = 0.0;
#else
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
#endif
            }
        }
    }

    if (bflag & (1u << OCT_DIR_DOWN)) {
        y = pmin[1] + jb * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x     = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp    = IDX(i, jb, k);

#ifdef EM1_DIRICHLET_BDY
                f_rhs[pp] = 0.0;
#else
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
#endif
            }
        }
    }

    if (bflag & (1u << OCT_DIR_UP)) {
        y = pmin[1] + (je - 1) * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x     = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp    = IDX(i, (je - 1), k);

#ifdef EM1_DIRICHLET_BDY
                f_rhs[pp] = 0.0;
#else
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
#endif
            }
        }
    }

    if (bflag & (1u << OCT_DIR_BACK)) {
        z = pmin[2] + kb * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x     = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp    = IDX(i, j, kb);

#ifdef EM1_DIRICHLET_BDY
                f_rhs[pp] = 0.0;
#else
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
#endif
            }
        }
    }

    if (bflag & (1u << OCT_DIR_FRONT)) {
        z = pmin[2] + (ke - 1) * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x     = pmin[0] + i * hx;
                inv_r = 1.0 / sqrt(x * x + y * y + z * z);
                pp    = IDX(i, j, (ke - 1));

#ifdef EM1_DIRICHLET_BDY
                f_rhs[pp] = 0.0;
#else
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
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
