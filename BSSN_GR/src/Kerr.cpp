double mass   = bssn::BH1.getBHMass();
double a_spin = bssn::BH1.getBHSpin();
// double mass = emda::EMDA_SINGLE_BH_MASS ;
// double a_spin = emda::EMDA_SINGLE_BH_SPIN ;
double rbar;
double cyl_rho_bar;
double xbar            = xx;
double ybar            = yy;
double zbar            = zz;

double cyl_rho_bar_sqr = xx * xx + yy * yy;
double rbar_sqr        = cyl_rho_bar_sqr + zz * zz;

double cos_theta;
double sin_theta;

if (fabs(xbar) > 1.0e-4 || fabs(ybar) > 1.0e-4)  // not on axis
{
    cyl_rho_bar = sqrt(cyl_rho_bar_sqr);
    rbar        = sqrt(rbar_sqr);
    cos_theta   = zbar / rbar;
    sin_theta   = cyl_rho_bar / rbar;
} else  // on axis
{
    cyl_rho_bar = 1.0e-10;

    if (fabs(zbar) > 1.0e-4)  // on axis but not at the origin
    {
        rbar      = sqrt(rbar_sqr);
        cos_theta = zbar / rbar;
        sin_theta = 0.0;  // cyl_rho_bar / rbar ;
    } else                // at the origin
    {
        rbar      = 1.0e-10;
        cos_theta = 0.0;
        sin_theta = 0.0;
    }
}

double r_boylin;
double f_0;
double rho_bar_sqr;
double Delta_bar;
double Sigma_bar;
double CC;
double beta_phi_up;

f_0 = (rbar + 0.5 * mass) * (rbar + 0.5 * mass) - 0.25 * a_spin * a_spin;
rho_bar_sqr = f_0 * f_0 + a_spin * a_spin * zbar * zbar;
Delta_bar   = f_0 * f_0 - 2 * mass * rbar * f_0 + a_spin * a_spin * rbar * rbar;
Sigma_bar   = pow(f_0 * f_0 + a_spin * a_spin * rbar * rbar, 2) -
            a_spin * a_spin * (xbar * xbar + ybar * ybar) * Delta_bar;
CC                = rho_bar_sqr * rho_bar_sqr / Sigma_bar;
beta_phi_up       = -a_spin * 2 * mass * f_0 * pow(rbar, 3) / Sigma_bar;

r_boylin          = f_0 / rbar;

var[VAR::U_ALPHA] = sqrt(Delta_bar * rho_bar_sqr / Sigma_bar);
var[VAR::U_CHI]   = pow(rbar, 4) / pow(Sigma_bar * rho_bar_sqr, 1.0 / 3.0);
// var[VAR::U_CHI] = pow(rbar/(rbar+0.5*mass), 4);
var[VAR::U_K]     = 0.0;

var[VAR::U_BETA0] = -ybar * beta_phi_up;  // shift_x
var[VAR::U_BETA1] = xbar * beta_phi_up;   // shift_y
var[VAR::U_BETA2] = 0.0;                  // shift_z

// var[VAR::U_CAP_GT0] = 0.0 ;  // Gt_x
// var[VAR::U_CAP_GT1] = 0.0 ;  // Gt_y
// var[VAR::U_CAP_GT2] = 0.0 ;  // Gt_z

var[VAR::U_GT0] =
    pow(a_spin, 2) * xbar * rho_bar_sqr * pow(CC, -4.0 / 3.0) / 3.0 /
    pow(Sigma_bar, 2) *
    (2 * pow(sin_theta, 2) * sqrt(Delta_bar) *
         (rho_bar_sqr * (f_0 - mass * rbar) +
          2 * f_0 * rbar * (2 * mass * f_0)) +
     2 * pow(cos_theta, 2) * (rho_bar_sqr * Delta_bar - 2 * Sigma_bar) -
     3 * rho_bar_sqr * (rho_bar_sqr + 2 * mass * f_0 * rbar));

var[VAR::U_GT1] =
    pow(a_spin, 2) * ybar * rho_bar_sqr * pow(CC, -4.0 / 3.0) / 3.0 /
    pow(Sigma_bar, 2) *
    (2 * pow(sin_theta, 2) * sqrt(Delta_bar) *
         (rho_bar_sqr * (f_0 - mass * rbar) +
          2 * f_0 * rbar * (2 * mass * f_0)) +
     2 * pow(cos_theta, 2) * (rho_bar_sqr * Delta_bar - 2 * Sigma_bar) -
     3 * rho_bar_sqr * (rho_bar_sqr + 2 * mass * f_0 * rbar));

var[VAR::U_GT2] = 2 * pow(a_spin, 2) * zbar * pow(sin_theta, 2) * rho_bar_sqr *
                  pow(CC, -4.0 / 3.0) / 3.0 / pow(Sigma_bar, 2) *
                  (sqrt(Delta_bar) * (rho_bar_sqr * (f_0 - mass * rbar) +
                                      2 * f_0 * rbar * (2 * mass * f_0)) -
                   (rho_bar_sqr * Delta_bar - 2 * Sigma_bar));

var[VAR::U_B0] = 0.0;  // gaugeB0
var[VAR::U_B1] = 0.0;  // gaugeB1
var[VAR::U_B2] = 0.0;  // gaugeB2
// var[VAR::U_GAUGEB0] = Sigma_bar ;  // gaugeB0
// var[VAR::U_GAUGEB1] = Delta_bar ;  // gaugeB1
// var[VAR::U_GAUGEB2] = CC ;  // gaugeB2

// var[VAR::U_GT00] = pow(CC, -2.0/3.0) * (xbar*xbar*CC + ybar*ybar) /
// cyl_rho_bar ;
var[VAR::U_SYMGT0] =
    pow(CC, -2.0 / 3.0) * (CC - (CC - 1.0) * ybar * ybar / cyl_rho_bar);
// var[VAR::U_GT00] = 1.0 ;
var[VAR::U_SYMGT1] =
    xbar * ybar * pow(CC, -2.0 / 3.0) * (CC - 1.0) / cyl_rho_bar;
// var[VAR::U_GT01] = 0.0 ;
var[VAR::U_SYMGT2] = 0.0;
// var[VAR::U_GT11] = pow(CC, -2.0/3.0) * (ybar*ybar*CC + xbar*xbar) /
// cyl_rho_bar ;
var[VAR::U_SYMGT3] =
    pow(CC, -2.0 / 3.0) * (CC - (CC - 1.0) * xbar * xbar / cyl_rho_bar);
// var[VAR::U_GT11] = 1.0 ;
var[VAR::U_SYMGT4] = 0.0;
var[VAR::U_SYMGT5] = pow(CC, 1.0 / 3.0);

var[VAR::U_SYMAT0] =
    2 * a_spin * xbar * ybar * rbar / rho_bar_sqr /
    pow(Sigma_bar * rho_bar_sqr, 5.0 / 6.0) *
    (mass * rho_bar_sqr * (pow(f_0, 2) + pow(a_spin * rbar, 2)) -
     f_0 * (2 * mass * f_0) *
         (rho_bar_sqr + pow(f_0, 2) + pow(a_spin * rbar, 2)) +
     pow(a_spin * zbar, 2) * sqrt(Delta_bar) * 2 * mass * f_0);

var[VAR::U_SYMAT1] =
    -a_spin * (xbar * xbar - ybar * ybar) * rbar / rho_bar_sqr /
    pow(Sigma_bar * rho_bar_sqr, 5.0 / 6.0) *
    (mass * rho_bar_sqr * (pow(f_0, 2) + pow(a_spin * rbar, 2)) -
     f_0 * (2 * mass * f_0) *
         (rho_bar_sqr + pow(f_0, 2) + pow(a_spin * rbar, 2)) +
     pow(a_spin * zbar, 2) * sqrt(Delta_bar) * 2 * mass * f_0);

var[VAR::U_SYMAT2] =
    a_spin * ybar * zbar * rbar / rho_bar_sqr /
    pow(Sigma_bar * rho_bar_sqr, 5.0 / 6.0) *
    (mass * rho_bar_sqr * (pow(f_0, 2) + pow(a_spin * rbar, 2)) -
     f_0 * (2 * mass * f_0) *
         (rho_bar_sqr + pow(f_0, 2) + pow(a_spin * rbar, 2)) -
     pow(a_spin, 2) * (xbar * xbar + ybar * ybar) * sqrt(Delta_bar) * 2 * mass *
         f_0);

var[VAR::U_SYMAT3] = -var[VAR::U_SYMAT0];

var[VAR::U_SYMAT3] =
    -a_spin * xbar * zbar * rbar / rho_bar_sqr /
    pow(Sigma_bar * rho_bar_sqr, 5.0 / 6.0) *
    (mass * rho_bar_sqr * (pow(f_0, 2) + pow(a_spin * rbar, 2)) -
     f_0 * (2 * mass * f_0) *
         (rho_bar_sqr + pow(f_0, 2) + pow(a_spin * rbar, 2)) -
     pow(a_spin, 2) * (xbar * xbar + ybar * ybar) * sqrt(Delta_bar) * 2 * mass *
         f_0);

var[VAR::U_SYMAT5] = 0.0;

// var[VAR::U_DILATONPHI] = 0.0 ;  // the dilaton
// var[VAR::U_CAPITALPI]  = 0.0 ;  // time deriv of dilaton
// //var[VAR::U_DILATONPHI] = sin_theta ;  // the dilaton
// //var[VAR::U_CAPITALPI]  = cos_theta;  // time deriv of dilaton

// //var[VAR::U_KAPPA]     = Delta_bar ;   //  the axion
// //var[VAR::U_CAPITALXI] = f_0 ;   // time deriv of axion
// var[VAR::U_KAPPA]     = 0.0 ;   //  the axion
// var[VAR::U_CAPITALXI] = 0.0 ;   // time deriv of axion

// var[VAR::U_DAMPINGPSI] = 0.0 ;  // dampingPsi
// var[VAR::U_DAMPINGPHI] = 0.0 ;  // dampingPhi
// //var[VAR::U_DAMPINGPSI] = beta_phi_up ;  // dampingPsi
// //var[VAR::U_DAMPINGPHI] = 0.0 ;  // dampingPhi

// var[VAR::U_PERPE0] = 0.0 ;  // E0
// var[VAR::U_PERPE1] = 0.0 ;  // E1
// var[VAR::U_PERPE2] = 0.0 ;  // E2
// var[VAR::U_PERPE0] = 1.0 ;  // E0
// var[VAR::U_PERPE1] = xbar*xbar + ybar*ybar ;  // E1
// var[VAR::U_PERPE2] = CC ;  // E2
// var[VAR::U_PERPE0] = sin_theta ;
// var[VAR::U_PERPE1] = cos_theta ;
// var[VAR::U_PERPE2] = CC ;  // E2

// var[VAR::U_PERPB0] = 0.0 ;  // B0
// var[VAR::U_PERPB1] = 0.0 ;  // B1
// var[VAR::U_PERPB2] = 0.0 ;  // B2
// var[VAR::U_PERPB0] = sin_theta ;  // B0
// var[VAR::U_PERPB1] = cos_theta ;  // B1
// var[VAR::U_PERPB2] = 0.0 ;  // B2
// var[VAR::U_PERPB0] = rbar_sqr ;  // B0
// var[VAR::U_PERPB1] = rbar ;  // B1
// var[VAR::U_PERPB2] = zbar/rbar ;  // B2
