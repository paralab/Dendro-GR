#include "rhs.h"

using namespace std;
using namespace maxwell;

/*----------------------------------------------------------------------;
 *
 * RHS for non-linear sigma model
 *
 *----------------------------------------------------------------------*/
void maxwellRhs(double t, double **unzipVarsRHS, const double **uZipVars,
                const unsigned int &offset, const double *pmin,
                const double *pmax, const unsigned int *sz,
                const unsigned int &bflag) {
    // const double *chi = &uZipVars[VAR::U_CHI][offset];
    // const double *phi = &uZipVars[VAR::U_PHI][offset];

    const double *AX      = &uZipVars[VAR::U_AX][offset];
    const double *AY      = &uZipVars[VAR::U_AY][offset];
    const double *AZ      = &uZipVars[VAR::U_AZ][offset];

    const double *EX      = &uZipVars[VAR::U_EX][offset];
    const double *EY      = &uZipVars[VAR::U_EY][offset];
    const double *EZ      = &uZipVars[VAR::U_EZ][offset];

    const double *GAM     = &uZipVars[VAR::U_GAM][offset];
    const double *PSI     = &uZipVars[VAR::U_PSI][offset];

    // double *chi_rhs = &unzipVarsRHS[VAR::U_CHI][offset];
    // double *phi_rhs = &unzipVarsRHS[VAR::U_PHI][offset];

    double *AX_rhs        = &unzipVarsRHS[VAR::U_AX][offset];
    double *AY_rhs        = &unzipVarsRHS[VAR::U_AY][offset];
    double *AZ_rhs        = &unzipVarsRHS[VAR::U_AZ][offset];

    double *EX_rhs        = &unzipVarsRHS[VAR::U_EX][offset];
    double *EY_rhs        = &unzipVarsRHS[VAR::U_EY][offset];
    double *EZ_rhs        = &unzipVarsRHS[VAR::U_EZ][offset];

    double *GAM_rhs       = &unzipVarsRHS[VAR::U_GAM][offset];
    double *PSI_rhs       = &unzipVarsRHS[VAR::U_PSI][offset];

    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    double hx             = (pmax[0] - pmin[0]) / (nx - 1);
    double hy             = (pmax[1] - pmin[1]) / (ny - 1);
    double hz             = (pmax[2] - pmin[2]) / (nz - 1);

    int idx[3];

    unsigned int n = sz[0] * sz[1] * sz[2];

    maxwell::timer::t_deriv.start();

    /*******NLSM*****
    double *grad_0_chi = new double [n];
    double *grad_1_chi = new double [n];
    double *grad_2_chi = new double [n];

    double *grad_0_phi = new double [n];
    double *grad_1_phi = new double [n];
    double *grad_2_phi = new double [n];

    double *grad2_0_0_chi = new double [n];
    double *grad2_1_1_chi = new double [n];
    double *grad2_2_2_chi = new double [n];

    deriv_xx(grad2_0_0_chi, chi, hx, sz, bflag);
    deriv_yy(grad2_1_1_chi, chi, hy, sz, bflag);
    deriv_zz(grad2_2_2_chi, chi, hz, sz, bflag);
    *******************/

    double *grad_0_AX     = new double[n];
    double *grad_1_AX     = new double[n];
    double *grad_2_AX     = new double[n];

    double *grad_0_AY     = new double[n];
    double *grad_1_AY     = new double[n];
    double *grad_2_AY     = new double[n];

    double *grad_0_AZ     = new double[n];
    double *grad_1_AZ     = new double[n];
    double *grad_2_AZ     = new double[n];

    double *grad_0_EX     = new double[n];
    double *grad_1_EX     = new double[n];
    double *grad_2_EX     = new double[n];

    double *grad_0_EY     = new double[n];
    double *grad_1_EY     = new double[n];
    double *grad_2_EY     = new double[n];

    double *grad_0_EZ     = new double[n];
    double *grad_1_EZ     = new double[n];
    double *grad_2_EZ     = new double[n];

    double *grad_0_GAM    = new double[n];
    double *grad_1_GAM    = new double[n];
    double *grad_2_GAM    = new double[n];

    double *grad_0_PSI    = new double[n];
    double *grad_1_PSI    = new double[n];
    double *grad_2_PSI    = new double[n];

    double *grad2_0_0_AX  = new double[n];
    double *grad2_1_1_AX  = new double[n];
    double *grad2_2_2_AX  = new double[n];

    double *grad2_0_0_AY  = new double[n];
    double *grad2_1_1_AY  = new double[n];
    double *grad2_2_2_AY  = new double[n];

    double *grad2_0_0_AZ  = new double[n];
    double *grad2_1_1_AZ  = new double[n];
    double *grad2_2_2_AZ  = new double[n];

    double *grad2_0_0_PSI = new double[n];
    double *grad2_1_1_PSI = new double[n];
    double *grad2_2_2_PSI = new double[n];

    deriv_x(grad_0_AX, AX, hx, sz, bflag);
    deriv_y(grad_1_AX, AX, hy, sz, bflag);
    deriv_z(grad_2_AX, AX, hz, sz, bflag);

    deriv_x(grad_0_AY, AY, hx, sz, bflag);
    deriv_y(grad_1_AY, AY, hy, sz, bflag);
    deriv_z(grad_2_AY, AY, hz, sz, bflag);

    deriv_x(grad_0_AZ, AZ, hx, sz, bflag);
    deriv_y(grad_1_AZ, AZ, hy, sz, bflag);
    deriv_z(grad_2_AZ, AZ, hz, sz, bflag);

    deriv_x(grad_0_EX, EX, hx, sz, bflag);
    deriv_y(grad_1_EX, EX, hy, sz, bflag);
    deriv_z(grad_2_EX, EX, hz, sz, bflag);

    deriv_x(grad_0_EY, EY, hx, sz, bflag);
    deriv_y(grad_1_EY, EY, hy, sz, bflag);
    deriv_z(grad_2_EY, EY, hz, sz, bflag);

    deriv_x(grad_0_EZ, EZ, hx, sz, bflag);
    deriv_y(grad_1_EZ, EZ, hy, sz, bflag);
    deriv_z(grad_2_EZ, EZ, hz, sz, bflag);

    deriv_x(grad_0_GAM, GAM, hx, sz, bflag);
    deriv_y(grad_1_GAM, GAM, hy, sz, bflag);
    deriv_z(grad_2_GAM, GAM, hz, sz, bflag);

    deriv_x(grad_0_PSI, PSI, hx, sz, bflag);
    deriv_y(grad_1_PSI, PSI, hy, sz, bflag);
    deriv_z(grad_2_PSI, PSI, hz, sz, bflag);

    deriv_xx(grad2_0_0_AX, AX, hx, sz, bflag);
    deriv_yy(grad2_1_1_AX, AX, hy, sz, bflag);
    deriv_zz(grad2_2_2_AX, AX, hz, sz, bflag);

    deriv_xx(grad2_0_0_AY, AY, hx, sz, bflag);
    deriv_yy(grad2_1_1_AY, AY, hy, sz, bflag);
    deriv_zz(grad2_2_2_AY, AY, hz, sz, bflag);

    deriv_xx(grad2_0_0_AZ, AZ, hx, sz, bflag);
    deriv_yy(grad2_1_1_AZ, AZ, hy, sz, bflag);
    deriv_zz(grad2_2_2_AZ, AZ, hz, sz, bflag);

    deriv_xx(grad2_0_0_PSI, PSI, hx, sz, bflag);
    deriv_yy(grad2_1_1_PSI, PSI, hy, sz, bflag);
    deriv_zz(grad2_2_2_PSI, PSI, hz, sz, bflag);

    maxwell::timer::t_deriv.stop();

    DendroRegister double x;
    DendroRegister double y;
    DendroRegister double z;
    DendroRegister unsigned int pp;
    // DAVE'S CODE STARTS HERE
    double *jx = new double[n];
    for (int counter = 0; counter < n; counter++) {
        jx[counter] = 0;
    }
    bool antennaSource = false;
    if (antennaSource == true) {
        for (unsigned int k = 3; k < nz - 3; k++) {
            z = pmin[2] + k * hz;

            for (unsigned int j = 3; j < ny - 3; j++) {
                y = pmin[1] + j * hy;

                for (unsigned int i = 3; i < nx - 3; i++) {
                    x         = pmin[0] + i * hx;
                    double rr = sqrt(x * x + y * y);
                    double kk = 3.14 / 5;
                    pp        = i + nx * (j + ny * k);
                    if ((x * x < 1e-8) && (y < 0.0) && (z * z < 1e-8) &&
                        (y > -5.0)) {
                        jx[pp] = sin(0.5 * t - rr * kk);
                    } else if (((y - x) * (y - x) < 1e-10) && (z * z < 1e-10) &&
                               (x > 0) && (y > 0) && (x < 3) && (y < 3)) {
                        jx[pp] = sin(0.5 * t - rr * kk);
                    } else if (((y + x) * (y + x) < 1e-10) && (z * z < 1e-10) &&
                               (x < 0) && (y > 0) && (y < 3)) {
                        jx[pp] = sin(0.5 * t - rr * kk);
                    } else {
                        jx[pp] = 0;
                    }
                }
            }
        }

    }  // End of antenna if statement

    // DAVE'S CODE ENDS HERE

    // cout << "begin loop" << endl;
    for (unsigned int k = 3; k < nz - 3; k++) {
        z = pmin[2] + k * hz;

        for (unsigned int j = 3; j < ny - 3; j++) {
            y = pmin[1] + j * hy;

            for (unsigned int i = 3; i < nx - 3; i++) {
                x          = pmin[0] + i * hx;
                pp         = i + nx * (j + ny * k);

                AX_rhs[pp] = -EX[pp] - grad_0_PSI[pp];
                AY_rhs[pp] = -EY[pp] - grad_1_PSI[pp];
                AZ_rhs[pp] = -EZ[pp] - grad_2_PSI[pp];

                EX_rhs[pp] =
                    grad_0_GAM[pp] -
                    (grad2_0_0_AX[pp] + grad2_1_1_AX[pp] + grad2_2_2_AX[pp]) -
                    4 * 3.14 * jx[pp];
                EY_rhs[pp] =
                    grad_1_GAM[pp] -
                    (grad2_0_0_AY[pp] + grad2_1_1_AY[pp] + grad2_2_2_AY[pp]);
                EZ_rhs[pp] =
                    grad_2_GAM[pp] -
                    (grad2_0_0_AZ[pp] + grad2_1_1_AZ[pp] + grad2_2_2_AZ[pp]);

                GAM_rhs[pp] = -(grad2_0_0_PSI[pp] + grad2_1_1_PSI[pp] +
                                grad2_2_2_PSI[pp]);

                PSI_rhs[pp] = -GAM[pp];
            }
        }
    }

    if (bflag != 0) {
        maxwell::timer::t_bdyc.start();

        maxwell_bcs(AX_rhs, AX, grad_0_AX, grad_1_AX, grad_2_AX, pmin, pmax,
                    1.0, 0.0, sz, bflag);
        maxwell_bcs(AY_rhs, AY, grad_0_AY, grad_1_AY, grad_2_AY, pmin, pmax,
                    1.0, 0.0, sz, bflag);
        maxwell_bcs(AZ_rhs, AZ, grad_0_AZ, grad_1_AZ, grad_2_AZ, pmin, pmax,
                    1.0, 0.0, sz, bflag);

        maxwell_bcs(EX_rhs, EX, grad_0_EX, grad_1_EX, grad_2_EX, pmin, pmax,
                    2.0, 0.0, sz, bflag);
        maxwell_bcs(EY_rhs, EY, grad_0_EY, grad_1_EY, grad_2_EY, pmin, pmax,
                    2.0, 0.0, sz, bflag);
        maxwell_bcs(EZ_rhs, EZ, grad_0_EZ, grad_1_EZ, grad_2_EZ, pmin, pmax,
                    2.0, 0.0, sz, bflag);

        maxwell_bcs(GAM_rhs, GAM, grad_0_GAM, grad_1_GAM, grad_2_GAM, pmin,
                    pmax, 2.0, 0.0, sz, bflag);

        maxwell_bcs(PSI_rhs, PSI, grad_0_PSI, grad_1_PSI, grad_2_PSI, pmin,
                    pmax, 1.0, 0.0, sz, bflag);

        maxwell::timer::t_bdyc.stop();
    }

    /****8
    maxwell::timer::t_deriv.start();


    ko_deriv_x(grad_0_chi, chi, hx, sz, bflag);
    ko_deriv_y(grad_1_chi, chi, hy, sz, bflag);
    ko_deriv_z(grad_2_chi, chi, hz, sz, bflag);

    ko_deriv_x(grad_0_phi, phi, hx, sz, bflag);
    ko_deriv_y(grad_1_phi, phi, hy, sz, bflag);
    ko_deriv_z(grad_2_phi, phi, hz, sz, bflag);
    maxwell::timer::t_deriv.stop();

    maxwell::timer::t_rhs.start();

      const  double sigma = KO_DISS_SIGMA;


      for (unsigned int k = 3; k < nz-3; k++) {
        for (unsigned int j = 3; j < ny-3; j++) {
          for (unsigned int i = 3; i < nx-3; i++) {
            pp = i + nx*(j + ny*k);

            chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] +
    grad_2_chi[pp]); phi_rhs[pp]  += sigma * (grad_0_phi[pp] + grad_1_phi[pp] +
    grad_2_phi[pp]);

          }
        }
      }

      maxwell::timer::t_rhs.stop();
    ************/

    maxwell::timer::t_deriv.start();

    delete[] jx;

    delete[] grad_0_AX;
    delete[] grad_1_AX;
    delete[] grad_2_AX;

    delete[] grad_0_AY;
    delete[] grad_1_AY;
    delete[] grad_2_AY;

    delete[] grad_0_AZ;
    delete[] grad_1_AZ;
    delete[] grad_2_AZ;

    delete[] grad_0_EX;
    delete[] grad_1_EX;
    delete[] grad_2_EX;

    delete[] grad_0_EY;
    delete[] grad_1_EY;
    delete[] grad_2_EY;

    delete[] grad_0_EZ;
    delete[] grad_1_EZ;
    delete[] grad_2_EZ;

    delete[] grad_0_GAM;
    delete[] grad_1_GAM;
    delete[] grad_2_GAM;

    delete[] grad_0_PSI;
    delete[] grad_1_PSI;
    delete[] grad_2_PSI;

    delete[] grad2_0_0_AX;
    delete[] grad2_1_1_AX;
    delete[] grad2_2_2_AX;

    delete[] grad2_0_0_AY;
    delete[] grad2_1_1_AY;
    delete[] grad2_2_2_AY;

    delete[] grad2_0_0_AZ;
    delete[] grad2_1_1_AZ;
    delete[] grad2_2_2_AZ;

    delete[] grad2_0_0_PSI;
    delete[] grad2_1_1_PSI;
    delete[] grad2_2_2_PSI;

    maxwell::timer::t_deriv.stop();

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
void maxwell_bcs(double *f_rhs, const double *f, const double *dxf,
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

    unsigned int ib       = 3;
    unsigned int jb       = 3;
    unsigned int kb       = 3;
    unsigned int ie       = sz[0] - 3;
    unsigned int je       = sz[1] - 3;
    unsigned int ke       = sz[2] - 3;

    double x, y, z;
    unsigned int pp;
    double inv_r;

    if (bflag & (1u << OCT_DIR_LEFT)) {
        double x = pmin[0] + ib * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y         = pmin[1] + j * hy;
                pp        = IDX(ib, j, k);
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_RIGHT)) {
        x = pmin[0] + (ie - 1) * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y         = pmin[1] + j * hy;
                pp        = IDX((ie - 1), j, k);
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_DOWN)) {
        y = pmin[1] + jb * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, jb, k);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_UP)) {
        y = pmin[1] + (je - 1) * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, (je - 1), k);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_BACK)) {
        z = pmin[2] + kb * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, j, kb);

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }

    if (bflag & (1u << OCT_DIR_FRONT)) {
        z = pmin[2] + (ke - 1) * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, j, (ke - 1));

                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
