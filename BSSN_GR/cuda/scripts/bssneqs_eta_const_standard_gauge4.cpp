double DENDRO_igt0;
double DENDRO_igt1;
double DENDRO_igt2;
double DENDRO_igt3;
double DENDRO_igt4;
double DENDRO_igt5;
double DENDRO_C1_k0_0;
double DENDRO_C1_k0_1;
double DENDRO_C1_k0_2;
double DENDRO_C1_k0_3;
double DENDRO_C1_k0_4;
double DENDRO_C1_k0_5;
double DENDRO_C1_k1_0;
double DENDRO_C1_k1_1;
double DENDRO_C1_k1_2;
double DENDRO_C1_k1_3;
double DENDRO_C1_k1_4;
double DENDRO_C1_k1_5;
double DENDRO_C1_k2_0;
double DENDRO_C1_k2_1;
double DENDRO_C1_k2_2;
double DENDRO_C1_k2_3;
double DENDRO_C1_k2_4;
double DENDRO_C1_k2_5;
double DENDRO_C2_k0_0;
double DENDRO_C2_k0_1;
double DENDRO_C2_k0_2;
double DENDRO_C2_k0_3;
double DENDRO_C2_k0_4;
double DENDRO_C2_k0_5;
double DENDRO_C2_k1_0;
double DENDRO_C2_k1_1;
double DENDRO_C2_k1_2;
double DENDRO_C2_k1_3;
double DENDRO_C2_k1_4;
double DENDRO_C2_k1_5;
double DENDRO_C2_k2_0;
double DENDRO_C2_k2_1;
double DENDRO_C2_k2_2;
double DENDRO_C2_k2_3;
double DENDRO_C2_k2_4;
double DENDRO_C2_k2_5;
double DENDRO_C3_k0_0;
double DENDRO_C3_k0_1;
double DENDRO_C3_k0_2;
double DENDRO_C3_k0_3;
double DENDRO_C3_k0_4;
double DENDRO_C3_k0_5;
double DENDRO_C3_k1_0;
double DENDRO_C3_k1_1;
double DENDRO_C3_k1_2;
double DENDRO_C3_k1_3;
double DENDRO_C3_k1_4;
double DENDRO_C3_k1_5;
double DENDRO_C3_k2_0;
double DENDRO_C3_k2_1;
double DENDRO_C3_k2_2;
double DENDRO_C3_k2_3;
double DENDRO_C3_k2_4;
double DENDRO_C3_k2_5;
double DENDRO_RIJ0;
double DENDRO_RIJ1;
double DENDRO_RIJ2;
double DENDRO_RIJ3;
double DENDRO_RIJ4;
double DENDRO_RIJ5;
double DENDRO_Gtk0;
double DENDRO_Gtk1;
double DENDRO_Gtk2;


{
// Dendro: {{{ 
// Dendro: original ops: 9 
// Dendro: printing temp variables

// Dendro: printing variables
//--
a_rhs[pp] = -2*K[pp]*alpha[pp] + lambda[0]*((deriv_evars->grad_0_alpha[d_pp])*beta0[pp] + (deriv_evars->grad_1_alpha[d_pp])*beta1[pp] + (deriv_evars->grad_2_alpha[d_pp])*beta2[pp]);
// Dendro: reduced ops: 9
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
b_rhs0[pp] = (3.0/4.0)*B0[pp]*(alpha[pp]*lambda_f[1] + lambda_f[0]) + lambda[1]*((deriv_evars->grad_0_beta0[d_pp])*beta0[pp] + (deriv_evars->grad_1_beta0[d_pp])*beta1[pp] + (deriv_evars->grad_2_beta0[d_pp])*beta2[pp]);
// Dendro: reduced ops: 12
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
b_rhs1[pp] = (3.0/4.0)*B1[pp]*(alpha[pp]*lambda_f[1] + lambda_f[0]) + lambda[1]*((deriv_evars->grad_0_beta1[d_pp])*beta0[pp] + (deriv_evars->grad_1_beta1[d_pp])*beta1[pp] + (deriv_evars->grad_2_beta1[d_pp])*beta2[pp]);
// Dendro: reduced ops: 12
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
b_rhs2[pp] = (3.0/4.0)*B2[pp]*(alpha[pp]*lambda_f[1] + lambda_f[0]) + lambda[1]*((deriv_evars->grad_0_beta2[d_pp])*beta0[pp] + (deriv_evars->grad_1_beta2[d_pp])*beta1[pp] + (deriv_evars->grad_2_beta2[d_pp])*beta2[pp]);
// Dendro: reduced ops: 12
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 26 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*gt0[pp];

// Dendro: printing variables
//--
gt_rhs00[pp] = (4.0/3.0)*(deriv_evars->grad_0_beta0[d_pp])*gt0[pp] + 2*(deriv_evars->grad_0_beta1[d_pp])*gt1[pp] + 2*(deriv_evars->grad_0_beta2[d_pp])*gt2[pp] + (deriv_evars->grad_0_gt0[d_pp])*beta0[pp] - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_gt0[d_pp])*beta1[pp] - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 + (deriv_evars->grad_2_gt0[d_pp])*beta2[pp] - 2*At0[pp]*alpha[pp];
// Dendro: reduced ops: 24
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 26 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*gt1[pp];

// Dendro: printing variables
//--
gt_rhs01[pp] = (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta1[d_pp])*gt3[pp] + (deriv_evars->grad_0_beta2[d_pp])*gt4[pp] + (deriv_evars->grad_0_gt1[d_pp])*beta0[pp] + (deriv_evars->grad_1_beta0[d_pp])*gt0[pp] + (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_beta2[d_pp])*gt2[pp] + (deriv_evars->grad_1_gt1[d_pp])*beta1[pp] - 2.0/3.0*(deriv_evars->grad_2_beta2[d_pp])*gt1[pp] + (deriv_evars->grad_2_gt1[d_pp])*beta2[pp] - 2*At1[pp]*alpha[pp];
// Dendro: reduced ops: 25
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 26 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*gt2[pp];

// Dendro: printing variables
//--
gt_rhs02[pp] = (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta1[d_pp])*gt4[pp] + (deriv_evars->grad_0_beta2[d_pp])*gt5[pp] + (deriv_evars->grad_0_gt2[d_pp])*beta0[pp] - 2.0/3.0*(deriv_evars->grad_1_beta1[d_pp])*gt2[pp] + (deriv_evars->grad_1_gt2[d_pp])*beta1[pp] + (deriv_evars->grad_2_beta0[d_pp])*gt0[pp] + (deriv_evars->grad_2_beta1[d_pp])*gt1[pp] + (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 + (deriv_evars->grad_2_gt2[d_pp])*beta2[pp] - 2*At2[pp]*alpha[pp];
// Dendro: reduced ops: 25
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 26 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*gt3[pp];

// Dendro: printing variables
//--
gt_rhs11[pp] = -(deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_gt3[d_pp])*beta0[pp] + 2*(deriv_evars->grad_1_beta0[d_pp])*gt1[pp] + (4.0/3.0)*(deriv_evars->grad_1_beta1[d_pp])*gt3[pp] + 2*(deriv_evars->grad_1_beta2[d_pp])*gt4[pp] + (deriv_evars->grad_1_gt3[d_pp])*beta1[pp] - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 + (deriv_evars->grad_2_gt3[d_pp])*beta2[pp] - 2*At3[pp]*alpha[pp];
// Dendro: reduced ops: 24
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 26 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*gt4[pp];

// Dendro: printing variables
//--
gt_rhs12[pp] = -2.0/3.0*(deriv_evars->grad_0_beta0[d_pp])*gt4[pp] + (deriv_evars->grad_0_gt4[d_pp])*beta0[pp] + (deriv_evars->grad_1_beta0[d_pp])*gt2[pp] + (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_beta2[d_pp])*gt5[pp] + (deriv_evars->grad_1_gt4[d_pp])*beta1[pp] + (deriv_evars->grad_2_beta0[d_pp])*gt1[pp] + (deriv_evars->grad_2_beta1[d_pp])*gt3[pp] + (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 + (deriv_evars->grad_2_gt4[d_pp])*beta2[pp] - 2*At4[pp]*alpha[pp];
// Dendro: reduced ops: 25
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 26 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*gt5[pp];

// Dendro: printing variables
//--
gt_rhs22[pp] = -(deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_gt5[d_pp])*beta0[pp] - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_gt5[d_pp])*beta1[pp] + 2*(deriv_evars->grad_2_beta0[d_pp])*gt2[pp] + 2*(deriv_evars->grad_2_beta1[d_pp])*gt4[pp] + (4.0/3.0)*(deriv_evars->grad_2_beta2[d_pp])*gt5[pp] + (deriv_evars->grad_2_gt5[d_pp])*beta2[pp] - 2*At5[pp]*alpha[pp];
// Dendro: reduced ops: 24
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 16 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*chi[pp];

// Dendro: printing variables
//--
chi_rhs[pp] = (deriv_evars->grad_0_chi[d_pp])*beta0[pp] + (deriv_evars->grad_1_chi[d_pp])*beta1[pp] + (deriv_evars->grad_2_chi[d_pp])*beta2[pp] + DENDRO_0*K[pp]*alpha[pp] - DENDRO_0*((deriv_evars->grad_0_beta0[d_pp]) + (deriv_evars->grad_1_beta1[d_pp]) + (deriv_evars->grad_2_beta2[d_pp]));
// Dendro: reduced ops: 14
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 114 
// Dendro: printing temp variables
const double DENDRO_0 = gt3[pp]*gt5[pp];
const double DENDRO_1 = pow(gt4[pp], 2);
const double DENDRO_2 = pow(gt1[pp], 2);
const double DENDRO_3 = pow(gt2[pp], 2);
const double DENDRO_4 = gt2[pp]*gt4[pp];
const double DENDRO_5 = 1.0/(-DENDRO_0*gt0[pp] + DENDRO_1*gt0[pp] + DENDRO_2*gt5[pp] + DENDRO_3*gt3[pp] - 2*DENDRO_4*gt1[pp]);

// Dendro: printing variables
//--
DENDRO_igt0 = -DENDRO_5*(DENDRO_0 - DENDRO_1);
//--
DENDRO_igt1 = DENDRO_5*(-DENDRO_4 + gt1[pp]*gt5[pp]);
//--
DENDRO_igt2 = -DENDRO_5*(gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp]);
//--
DENDRO_igt3 = -DENDRO_5*(-DENDRO_3 + gt0[pp]*gt5[pp]);
//--
DENDRO_igt4 = DENDRO_5*(gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp]);
//--
DENDRO_igt5 = -DENDRO_5*(-DENDRO_2 + gt0[pp]*gt3[pp]);
// Dendro: reduced ops: 39
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
DENDRO_C1_k0_0 = 0.5*(deriv_evars->grad_0_gt0[d_pp]);
//--
DENDRO_C1_k0_1 = 0.5*(deriv_evars->grad_1_gt0[d_pp]);
//--
DENDRO_C1_k0_2 = 0.5*(deriv_evars->grad_2_gt0[d_pp]);
//--
DENDRO_C1_k0_3 = -0.5*(deriv_evars->grad_0_gt3[d_pp]) + 1.0*(deriv_evars->grad_1_gt1[d_pp]);
//--
DENDRO_C1_k0_4 = 0.5*(-(deriv_evars->grad_0_gt4[d_pp]) + (deriv_evars->grad_1_gt2[d_pp]) + (deriv_evars->grad_2_gt1[d_pp]));
//--
DENDRO_C1_k0_5 = -0.5*(deriv_evars->grad_0_gt5[d_pp]) + 1.0*(deriv_evars->grad_2_gt2[d_pp]);
// Dendro: reduced ops: 12
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
DENDRO_C1_k1_0 = 1.0*(deriv_evars->grad_0_gt1[d_pp]) - 0.5*(deriv_evars->grad_1_gt0[d_pp]);
//--
DENDRO_C1_k1_1 = 0.5*(deriv_evars->grad_0_gt3[d_pp]);
//--
DENDRO_C1_k1_2 = 0.5*((deriv_evars->grad_0_gt4[d_pp]) - (deriv_evars->grad_1_gt2[d_pp]) + (deriv_evars->grad_2_gt1[d_pp]));
//--
DENDRO_C1_k1_3 = 0.5*(deriv_evars->grad_1_gt3[d_pp]);
//--
DENDRO_C1_k1_4 = 0.5*(deriv_evars->grad_2_gt3[d_pp]);
//--
DENDRO_C1_k1_5 = -0.5*(deriv_evars->grad_1_gt5[d_pp]) + 1.0*(deriv_evars->grad_2_gt4[d_pp]);
// Dendro: reduced ops: 12
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
DENDRO_C1_k2_0 = 1.0*(deriv_evars->grad_0_gt2[d_pp]) - 0.5*(deriv_evars->grad_2_gt0[d_pp]);
//--
DENDRO_C1_k2_1 = 0.5*((deriv_evars->grad_0_gt4[d_pp]) + (deriv_evars->grad_1_gt2[d_pp]) - (deriv_evars->grad_2_gt1[d_pp]));
//--
DENDRO_C1_k2_2 = 0.5*(deriv_evars->grad_0_gt5[d_pp]);
//--
DENDRO_C1_k2_3 = 1.0*(deriv_evars->grad_1_gt4[d_pp]) - 0.5*(deriv_evars->grad_2_gt3[d_pp]);
//--
DENDRO_C1_k2_4 = 0.5*(deriv_evars->grad_1_gt5[d_pp]);
//--
DENDRO_C1_k2_5 = 0.5*(deriv_evars->grad_2_gt5[d_pp]);
// Dendro: reduced ops: 12
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 30 
// Dendro: printing temp variables

// Dendro: printing variables
//--
DENDRO_C2_k0_0 = DENDRO_C1_k0_0*DENDRO_igt0 + DENDRO_C1_k1_0*DENDRO_igt1 + DENDRO_C1_k2_0*DENDRO_igt2;
//--
DENDRO_C2_k0_1 = DENDRO_C1_k0_1*DENDRO_igt0 + DENDRO_C1_k1_1*DENDRO_igt1 + DENDRO_C1_k2_1*DENDRO_igt2;
//--
DENDRO_C2_k0_2 = DENDRO_C1_k0_2*DENDRO_igt0 + DENDRO_C1_k1_2*DENDRO_igt1 + DENDRO_C1_k2_2*DENDRO_igt2;
//--
DENDRO_C2_k0_3 = DENDRO_C1_k0_3*DENDRO_igt0 + DENDRO_C1_k1_3*DENDRO_igt1 + DENDRO_C1_k2_3*DENDRO_igt2;
//--
DENDRO_C2_k0_4 = DENDRO_C1_k0_4*DENDRO_igt0 + DENDRO_C1_k1_4*DENDRO_igt1 + DENDRO_C1_k2_4*DENDRO_igt2;
//--
DENDRO_C2_k0_5 = DENDRO_C1_k0_5*DENDRO_igt0 + DENDRO_C1_k1_5*DENDRO_igt1 + DENDRO_C1_k2_5*DENDRO_igt2;
// Dendro: reduced ops: 30
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 30 
// Dendro: printing temp variables

// Dendro: printing variables
//--
DENDRO_C2_k1_0 = DENDRO_C1_k0_0*DENDRO_igt1 + DENDRO_C1_k1_0*DENDRO_igt3 + DENDRO_C1_k2_0*DENDRO_igt4;
//--
DENDRO_C2_k1_1 = DENDRO_C1_k0_1*DENDRO_igt1 + DENDRO_C1_k1_1*DENDRO_igt3 + DENDRO_C1_k2_1*DENDRO_igt4;
//--
DENDRO_C2_k1_2 = DENDRO_C1_k0_2*DENDRO_igt1 + DENDRO_C1_k1_2*DENDRO_igt3 + DENDRO_C1_k2_2*DENDRO_igt4;
//--
DENDRO_C2_k1_3 = DENDRO_C1_k0_3*DENDRO_igt1 + DENDRO_C1_k1_3*DENDRO_igt3 + DENDRO_C1_k2_3*DENDRO_igt4;
//--
DENDRO_C2_k1_4 = DENDRO_C1_k0_4*DENDRO_igt1 + DENDRO_C1_k1_4*DENDRO_igt3 + DENDRO_C1_k2_4*DENDRO_igt4;
//--
DENDRO_C2_k1_5 = DENDRO_C1_k0_5*DENDRO_igt1 + DENDRO_C1_k1_5*DENDRO_igt3 + DENDRO_C1_k2_5*DENDRO_igt4;
// Dendro: reduced ops: 30
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 30 
// Dendro: printing temp variables

// Dendro: printing variables
//--
DENDRO_C2_k2_0 = DENDRO_C1_k0_0*DENDRO_igt2 + DENDRO_C1_k1_0*DENDRO_igt4 + DENDRO_C1_k2_0*DENDRO_igt5;
//--
DENDRO_C2_k2_1 = DENDRO_C1_k0_1*DENDRO_igt2 + DENDRO_C1_k1_1*DENDRO_igt4 + DENDRO_C1_k2_1*DENDRO_igt5;
//--
DENDRO_C2_k2_2 = DENDRO_C1_k0_2*DENDRO_igt2 + DENDRO_C1_k1_2*DENDRO_igt4 + DENDRO_C1_k2_2*DENDRO_igt5;
//--
DENDRO_C2_k2_3 = DENDRO_C1_k0_3*DENDRO_igt2 + DENDRO_C1_k1_3*DENDRO_igt4 + DENDRO_C1_k2_3*DENDRO_igt5;
//--
DENDRO_C2_k2_4 = DENDRO_C1_k0_4*DENDRO_igt2 + DENDRO_C1_k1_4*DENDRO_igt4 + DENDRO_C1_k2_4*DENDRO_igt5;
//--
DENDRO_C2_k2_5 = DENDRO_C1_k0_5*DENDRO_igt2 + DENDRO_C1_k1_5*DENDRO_igt4 + DENDRO_C1_k2_5*DENDRO_igt5;
// Dendro: reduced ops: 30
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 60 
// Dendro: printing temp variables
const double DENDRO_0 = 1.0/chi[pp];
const double DENDRO_1 = 0.5*(deriv_evars->grad_0_chi[d_pp])*DENDRO_igt0 + 0.5*(deriv_evars->grad_1_chi[d_pp])*DENDRO_igt1 + 0.5*(deriv_evars->grad_2_chi[d_pp])*DENDRO_igt2;
const double DENDRO_2 = DENDRO_0*DENDRO_1;

// Dendro: printing variables
//--
DENDRO_C3_k0_0 = -DENDRO_0*(1.0*(deriv_evars->grad_0_chi[d_pp]) - DENDRO_1*gt0[pp]) + DENDRO_C2_k0_0;
//--
DENDRO_C3_k0_1 = -DENDRO_0*(0.5*(deriv_evars->grad_1_chi[d_pp]) - DENDRO_1*gt1[pp]) + DENDRO_C2_k0_1;
//--
DENDRO_C3_k0_2 = -DENDRO_0*(0.5*(deriv_evars->grad_2_chi[d_pp]) - DENDRO_1*gt2[pp]) + DENDRO_C2_k0_2;
//--
DENDRO_C3_k0_3 = DENDRO_2*gt3[pp] + DENDRO_C2_k0_3;
//--
DENDRO_C3_k0_4 = DENDRO_2*gt4[pp] + DENDRO_C2_k0_4;
//--
DENDRO_C3_k0_5 = DENDRO_2*gt5[pp] + DENDRO_C2_k0_5;
// Dendro: reduced ops: 31
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 60 
// Dendro: printing temp variables
const double DENDRO_0 = 1.0/chi[pp];
const double DENDRO_1 = 0.5*(deriv_evars->grad_0_chi[d_pp])*DENDRO_igt1 + 0.5*(deriv_evars->grad_1_chi[d_pp])*DENDRO_igt3 + 0.5*(deriv_evars->grad_2_chi[d_pp])*DENDRO_igt4;
const double DENDRO_2 = DENDRO_0*DENDRO_1;

// Dendro: printing variables
//--
DENDRO_C3_k1_0 = DENDRO_2*gt0[pp] + DENDRO_C2_k1_0;
//--
DENDRO_C3_k1_1 = -DENDRO_0*(0.5*(deriv_evars->grad_0_chi[d_pp]) - DENDRO_1*gt1[pp]) + DENDRO_C2_k1_1;
//--
DENDRO_C3_k1_2 = DENDRO_2*gt2[pp] + DENDRO_C2_k1_2;
//--
DENDRO_C3_k1_3 = -DENDRO_0*(1.0*(deriv_evars->grad_1_chi[d_pp]) - DENDRO_1*gt3[pp]) + DENDRO_C2_k1_3;
//--
DENDRO_C3_k1_4 = -DENDRO_0*(0.5*(deriv_evars->grad_2_chi[d_pp]) - DENDRO_1*gt4[pp]) + DENDRO_C2_k1_4;
//--
DENDRO_C3_k1_5 = DENDRO_2*gt5[pp] + DENDRO_C2_k1_5;
// Dendro: reduced ops: 31
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 60 
// Dendro: printing temp variables
const double DENDRO_0 = 1.0/chi[pp];
const double DENDRO_1 = 0.5*(deriv_evars->grad_0_chi[d_pp])*DENDRO_igt2 + 0.5*(deriv_evars->grad_1_chi[d_pp])*DENDRO_igt4 + 0.5*(deriv_evars->grad_2_chi[d_pp])*DENDRO_igt5;
const double DENDRO_2 = DENDRO_0*DENDRO_1;

// Dendro: printing variables
//--
DENDRO_C3_k2_0 = DENDRO_2*gt0[pp] + DENDRO_C2_k2_0;
//--
DENDRO_C3_k2_1 = DENDRO_2*gt1[pp] + DENDRO_C2_k2_1;
//--
DENDRO_C3_k2_2 = -DENDRO_0*(0.5*(deriv_evars->grad_0_chi[d_pp]) - DENDRO_1*gt2[pp]) + DENDRO_C2_k2_2;
//--
DENDRO_C3_k2_3 = DENDRO_2*gt3[pp] + DENDRO_C2_k2_3;
//--
DENDRO_C3_k2_4 = -DENDRO_0*(0.5*(deriv_evars->grad_1_chi[d_pp]) - DENDRO_1*gt4[pp]) + DENDRO_C2_k2_4;
//--
DENDRO_C3_k2_5 = -DENDRO_0*(1.0*(deriv_evars->grad_2_chi[d_pp]) - DENDRO_1*gt5[pp]) + DENDRO_C2_k2_5;
// Dendro: reduced ops: 31
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 363 
// Dendro: printing temp variables
const double DENDRO_0 = pow(chi[pp], 2);
const double DENDRO_1 = pow((deriv_evars->grad_0_chi[d_pp]), 2);
const double DENDRO_2 = 2*chi[pp];
const double DENDRO_3 = 3*(deriv_evars->grad_0_chi[d_pp]);
const double DENDRO_4 = 2*DENDRO_igt1;
const double DENDRO_5 = 2*DENDRO_igt2;
const double DENDRO_6 = 2*DENDRO_igt4;
const double DENDRO_7 = DENDRO_C2_k0_0*DENDRO_igt0;
const double DENDRO_8 = DENDRO_4*DENDRO_C2_k0_1 + DENDRO_5*DENDRO_C2_k0_2 + DENDRO_6*DENDRO_C2_k0_4 + DENDRO_7 + DENDRO_C2_k0_3*DENDRO_igt3 + DENDRO_C2_k0_5*DENDRO_igt5;
const double DENDRO_9 = DENDRO_C2_k1_0*DENDRO_igt0;
const double DENDRO_10 = DENDRO_4*DENDRO_C2_k1_1 + DENDRO_5*DENDRO_C2_k1_2 + DENDRO_6*DENDRO_C2_k1_4 + DENDRO_9 + DENDRO_C2_k1_3*DENDRO_igt3 + DENDRO_C2_k1_5*DENDRO_igt5;
const double DENDRO_11 = DENDRO_C2_k2_0*DENDRO_igt0;
const double DENDRO_12 = DENDRO_11 + DENDRO_4*DENDRO_C2_k2_1 + DENDRO_5*DENDRO_C2_k2_2 + DENDRO_6*DENDRO_C2_k2_4 + DENDRO_C2_k2_3*DENDRO_igt3 + DENDRO_C2_k2_5*DENDRO_igt5;
const double DENDRO_13 = 2*DENDRO_C1_k0_1;
const double DENDRO_14 = 2*DENDRO_C1_k0_3;
const double DENDRO_15 = 4*DENDRO_igt3;
const double DENDRO_16 = 2*DENDRO_C1_k0_4;
const double DENDRO_17 = 4*DENDRO_igt5;
const double DENDRO_18 = 2*DENDRO_C1_k0_2;
const double DENDRO_19 = 2*DENDRO_C1_k0_5;
const double DENDRO_20 = DENDRO_C1_k0_0*DENDRO_C2_k0_1;
const double DENDRO_21 = DENDRO_C1_k0_1*DENDRO_C2_k0_0;
const double DENDRO_22 = 4*DENDRO_igt1;
const double DENDRO_23 = DENDRO_C1_k0_0*DENDRO_C2_k0_2;
const double DENDRO_24 = DENDRO_C1_k0_2*DENDRO_C2_k0_0;
const double DENDRO_25 = 4*DENDRO_igt2;
const double DENDRO_26 = DENDRO_C1_k0_1*DENDRO_C2_k0_2;
const double DENDRO_27 = DENDRO_C1_k0_2*DENDRO_C2_k0_1;
const double DENDRO_28 = 4*DENDRO_igt4;

// Dendro: printing variables
//--
DENDRO_RIJ0 = -1.0/4.0*(-DENDRO_0*(4.0*(deriv_evars->grad_0_Gt0[d_pp])*gt0[pp] + 4.0*(deriv_evars->grad_0_Gt1[d_pp])*gt1[pp] + 4.0*(deriv_evars->grad_0_Gt2[d_pp])*gt2[pp] + 4.0*DENDRO_10*DENDRO_C1_k0_1 + 4*DENDRO_11*(DENDRO_18 + DENDRO_C1_k2_0) + 4.0*DENDRO_12*DENDRO_C1_k0_2 + DENDRO_15*DENDRO_C2_k1_1*(DENDRO_14 + DENDRO_C1_k1_1) + DENDRO_15*DENDRO_C2_k2_1*(DENDRO_16 + DENDRO_C1_k2_1) + DENDRO_17*DENDRO_C2_k1_2*(DENDRO_16 + DENDRO_C1_k1_2) + DENDRO_17*DENDRO_C2_k2_2*(DENDRO_19 + DENDRO_C1_k2_2) + DENDRO_22*(DENDRO_20 + 2*DENDRO_21) + DENDRO_22*(2*DENDRO_20 + DENDRO_21) + DENDRO_22*(DENDRO_13*DENDRO_C2_k1_1 + DENDRO_C1_k1_1*DENDRO_C2_k1_0) + DENDRO_22*(DENDRO_14*DENDRO_C2_k1_0 + DENDRO_C1_k1_0*DENDRO_C2_k1_1) + DENDRO_22*(DENDRO_16*DENDRO_C2_k2_0 + DENDRO_C1_k2_0*DENDRO_C2_k2_1) + DENDRO_22*(DENDRO_18*DENDRO_C2_k2_1 + DENDRO_C1_k2_1*DENDRO_C2_k2_0) + DENDRO_25*(DENDRO_23 + 2*DENDRO_24) + DENDRO_25*(2*DENDRO_23 + DENDRO_24) + DENDRO_25*(DENDRO_13*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_0) + DENDRO_25*(DENDRO_16*DENDRO_C2_k1_0 + DENDRO_C1_k1_0*DENDRO_C2_k1_2) + DENDRO_25*(DENDRO_18*DENDRO_C2_k2_2 + DENDRO_C1_k2_2*DENDRO_C2_k2_0) + DENDRO_25*(DENDRO_19*DENDRO_C2_k2_0 + DENDRO_C1_k2_0*DENDRO_C2_k2_2) + DENDRO_28*(DENDRO_26 + 2*DENDRO_27) + DENDRO_28*(2*DENDRO_26 + DENDRO_27) + DENDRO_28*(DENDRO_14*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_1) + DENDRO_28*(DENDRO_16*DENDRO_C2_k1_1 + DENDRO_C1_k1_1*DENDRO_C2_k1_2) + DENDRO_28*(DENDRO_16*DENDRO_C2_k2_2 + DENDRO_C1_k2_2*DENDRO_C2_k2_1) + DENDRO_28*(DENDRO_19*DENDRO_C2_k2_1 + DENDRO_C1_k2_1*DENDRO_C2_k2_2) + 12*DENDRO_7*DENDRO_C1_k0_0 + 4.0*DENDRO_8*DENDRO_C1_k0_0 + 4*DENDRO_9*(DENDRO_13 + DENDRO_C1_k1_0) + 12*DENDRO_C1_k0_1*DENDRO_C2_k0_1*DENDRO_igt3 + 12*DENDRO_C1_k0_2*DENDRO_C2_k0_2*DENDRO_igt5 - 2.0*DENDRO_igt0*deriv_evars->grad2_0_0_gt0[d_pp] - 4.0*DENDRO_igt1*deriv_evars->grad2_0_1_gt0[d_pp] - 4.0*DENDRO_igt2*deriv_evars->grad2_0_2_gt0[d_pp] - 2.0*DENDRO_igt3*deriv_evars->grad2_1_1_gt0[d_pp] - 4.0*DENDRO_igt4*deriv_evars->grad2_1_2_gt0[d_pp] - 2.0*DENDRO_igt5*deriv_evars->grad2_2_2_gt0[d_pp]) + DENDRO_1 + DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_C2_k0_0 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_C2_k1_0 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_C2_k2_0 - deriv_evars->grad2_0_0_chi[d_pp]) + gt0[pp]*(DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_8 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_10 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_12) + DENDRO_4*((deriv_evars->grad_1_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_1_chi[d_pp]) + DENDRO_5*((deriv_evars->grad_2_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_2_chi[d_pp]) + DENDRO_6*(3*(deriv_evars->grad_1_chi[d_pp])*(deriv_evars->grad_2_chi[d_pp]) - DENDRO_2*deriv_evars->grad2_1_2_chi[d_pp]) + DENDRO_igt0*(3*DENDRO_1 - DENDRO_2*deriv_evars->grad2_0_0_chi[d_pp]) + DENDRO_igt3*(3*pow((deriv_evars->grad_1_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_1_1_chi[d_pp]) + DENDRO_igt5*(3*pow((deriv_evars->grad_2_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_2_2_chi[d_pp])))/DENDRO_0;
// Dendro: reduced ops: 264
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 413 
// Dendro: printing temp variables
const double DENDRO_0 = pow(chi[pp], 2);
const double DENDRO_1 = (deriv_evars->grad_0_chi[d_pp])*(deriv_evars->grad_1_chi[d_pp]);
const double DENDRO_2 = 2*chi[pp];
const double DENDRO_3 = 2*DENDRO_igt1;
const double DENDRO_4 = 3*(deriv_evars->grad_2_chi[d_pp]);
const double DENDRO_5 = 2*DENDRO_igt2;
const double DENDRO_6 = 2*DENDRO_igt4;
const double DENDRO_7 = DENDRO_3*DENDRO_C2_k0_1 + DENDRO_5*DENDRO_C2_k0_2 + DENDRO_6*DENDRO_C2_k0_4 + DENDRO_C2_k0_0*DENDRO_igt0 + DENDRO_C2_k0_3*DENDRO_igt3 + DENDRO_C2_k0_5*DENDRO_igt5;
const double DENDRO_8 = DENDRO_3*DENDRO_C2_k1_1 + DENDRO_5*DENDRO_C2_k1_2 + DENDRO_6*DENDRO_C2_k1_4 + DENDRO_C2_k1_0*DENDRO_igt0 + DENDRO_C2_k1_3*DENDRO_igt3 + DENDRO_C2_k1_5*DENDRO_igt5;
const double DENDRO_9 = DENDRO_3*DENDRO_C2_k2_1 + DENDRO_5*DENDRO_C2_k2_2 + DENDRO_6*DENDRO_C2_k2_4 + DENDRO_C2_k2_0*DENDRO_igt0 + DENDRO_C2_k2_3*DENDRO_igt3 + DENDRO_C2_k2_5*DENDRO_igt5;
const double DENDRO_10 = 2.0*gt1[pp];
const double DENDRO_11 = 4*DENDRO_igt0;
const double DENDRO_12 = 4*DENDRO_igt1;
const double DENDRO_13 = 4*DENDRO_igt3;
const double DENDRO_14 = 4*DENDRO_igt5;
const double DENDRO_15 = DENDRO_C1_k1_1*DENDRO_C2_k1_1 + DENDRO_C1_k1_3*DENDRO_C2_k1_0;
const double DENDRO_16 = 4*DENDRO_igt2;
const double DENDRO_17 = DENDRO_C1_k1_1*DENDRO_C2_k1_2 + DENDRO_C1_k1_4*DENDRO_C2_k1_0;
const double DENDRO_18 = 4*DENDRO_igt4;
const double DENDRO_19 = DENDRO_C1_k1_3*DENDRO_C2_k1_2 + DENDRO_C1_k1_4*DENDRO_C2_k1_1;

// Dendro: printing variables
//--
DENDRO_RIJ1 = -1.0/4.0*(-DENDRO_0*((deriv_evars->grad_0_Gt0[d_pp])*DENDRO_10 + 2.0*(deriv_evars->grad_0_Gt1[d_pp])*gt3[pp] + 2.0*(deriv_evars->grad_0_Gt2[d_pp])*gt4[pp] + 2.0*(deriv_evars->grad_1_Gt0[d_pp])*gt0[pp] + (deriv_evars->grad_1_Gt1[d_pp])*DENDRO_10 + 2.0*(deriv_evars->grad_1_Gt2[d_pp])*gt2[pp] + DENDRO_11*(DENDRO_C1_k0_1*DENDRO_C2_k1_1 + 2*DENDRO_C1_k1_1*DENDRO_C2_k1_0) + DENDRO_11*(DENDRO_C1_k0_0*DENDRO_C2_k0_1 + DENDRO_C1_k0_1*DENDRO_C2_k0_0 + DENDRO_C1_k1_0*DENDRO_C2_k0_0) + DENDRO_11*(DENDRO_C1_k0_2*DENDRO_C2_k2_1 + DENDRO_C1_k1_2*DENDRO_C2_k2_0 + DENDRO_C1_k2_1*DENDRO_C2_k2_0) + DENDRO_12*(DENDRO_15 + DENDRO_C1_k0_1*DENDRO_C2_k1_3) + DENDRO_12*(DENDRO_15 + DENDRO_C1_k0_3*DENDRO_C2_k1_1) + DENDRO_12*(2*DENDRO_C1_k0_1*DENDRO_C2_k0_1 + DENDRO_C1_k1_1*DENDRO_C2_k0_0) + DENDRO_12*(DENDRO_C1_k0_0*DENDRO_C2_k0_3 + DENDRO_C1_k0_3*DENDRO_C2_k0_0 + DENDRO_C1_k1_0*DENDRO_C2_k0_1) + DENDRO_12*(DENDRO_C1_k0_2*DENDRO_C2_k2_3 + DENDRO_C1_k1_2*DENDRO_C2_k2_1 + DENDRO_C1_k2_3*DENDRO_C2_k2_0) + DENDRO_12*(DENDRO_C1_k0_4*DENDRO_C2_k2_1 + DENDRO_C1_k1_4*DENDRO_C2_k2_0 + DENDRO_C1_k2_1*DENDRO_C2_k2_1) + DENDRO_13*(DENDRO_C1_k0_3*DENDRO_C2_k1_3 + 2*DENDRO_C1_k1_3*DENDRO_C2_k1_1) + DENDRO_13*(DENDRO_C1_k0_1*DENDRO_C2_k0_3 + DENDRO_C1_k0_3*DENDRO_C2_k0_1 + DENDRO_C1_k1_1*DENDRO_C2_k0_1) + DENDRO_13*(DENDRO_C1_k0_4*DENDRO_C2_k2_3 + DENDRO_C1_k1_4*DENDRO_C2_k2_1 + DENDRO_C1_k2_3*DENDRO_C2_k2_1) + DENDRO_14*(DENDRO_C1_k0_4*DENDRO_C2_k1_4 + 2*DENDRO_C1_k1_4*DENDRO_C2_k1_2) + DENDRO_14*(DENDRO_C1_k0_2*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_2 + DENDRO_C1_k1_2*DENDRO_C2_k0_2) + DENDRO_14*(DENDRO_C1_k0_5*DENDRO_C2_k2_4 + DENDRO_C1_k1_5*DENDRO_C2_k2_2 + DENDRO_C1_k2_4*DENDRO_C2_k2_2) + DENDRO_16*(DENDRO_17 + DENDRO_C1_k0_1*DENDRO_C2_k1_4) + DENDRO_16*(DENDRO_17 + DENDRO_C1_k0_4*DENDRO_C2_k1_1) + DENDRO_16*(DENDRO_C1_k0_0*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_0 + DENDRO_C1_k1_0*DENDRO_C2_k0_2) + DENDRO_16*(DENDRO_C1_k0_1*DENDRO_C2_k0_2 + DENDRO_C1_k0_2*DENDRO_C2_k0_1 + DENDRO_C1_k1_2*DENDRO_C2_k0_0) + DENDRO_16*(DENDRO_C1_k0_2*DENDRO_C2_k2_4 + DENDRO_C1_k1_2*DENDRO_C2_k2_2 + DENDRO_C1_k2_4*DENDRO_C2_k2_0) + DENDRO_16*(DENDRO_C1_k0_5*DENDRO_C2_k2_1 + DENDRO_C1_k1_5*DENDRO_C2_k2_0 + DENDRO_C1_k2_1*DENDRO_C2_k2_2) + DENDRO_18*(DENDRO_19 + DENDRO_C1_k0_3*DENDRO_C2_k1_4) + DENDRO_18*(DENDRO_19 + DENDRO_C1_k0_4*DENDRO_C2_k1_3) + DENDRO_18*(DENDRO_C1_k0_1*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_1 + DENDRO_C1_k1_1*DENDRO_C2_k0_2) + DENDRO_18*(DENDRO_C1_k0_2*DENDRO_C2_k0_3 + DENDRO_C1_k0_3*DENDRO_C2_k0_2 + DENDRO_C1_k1_2*DENDRO_C2_k0_1) + DENDRO_18*(DENDRO_C1_k0_4*DENDRO_C2_k2_4 + DENDRO_C1_k1_4*DENDRO_C2_k2_2 + DENDRO_C1_k2_4*DENDRO_C2_k2_1) + DENDRO_18*(DENDRO_C1_k0_5*DENDRO_C2_k2_3 + DENDRO_C1_k1_5*DENDRO_C2_k2_1 + DENDRO_C1_k2_3*DENDRO_C2_k2_2) + 2.0*DENDRO_7*(DENDRO_C1_k0_1 + DENDRO_C1_k1_0) + 2.0*DENDRO_8*(DENDRO_C1_k0_3 + DENDRO_C1_k1_1) + 2.0*DENDRO_9*(DENDRO_C1_k0_4 + DENDRO_C1_k1_2) - 2.0*DENDRO_igt0*deriv_evars->grad2_0_0_gt1[d_pp] - 4.0*DENDRO_igt1*deriv_evars->grad2_0_1_gt1[d_pp] - 4.0*DENDRO_igt2*deriv_evars->grad2_0_2_gt1[d_pp] - 2.0*DENDRO_igt3*deriv_evars->grad2_1_1_gt1[d_pp] - 4.0*DENDRO_igt4*deriv_evars->grad2_1_2_gt1[d_pp] - 2.0*DENDRO_igt5*deriv_evars->grad2_2_2_gt1[d_pp]) + DENDRO_1 + DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_C2_k0_1 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_C2_k1_1 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_C2_k2_1 - deriv_evars->grad2_0_1_chi[d_pp]) + gt1[pp]*(DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_7 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_8 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_9) + DENDRO_3*(3*DENDRO_1 - DENDRO_2*deriv_evars->grad2_0_1_chi[d_pp]) + DENDRO_5*((deriv_evars->grad_0_chi[d_pp])*DENDRO_4 - DENDRO_2*deriv_evars->grad2_0_2_chi[d_pp]) + DENDRO_6*((deriv_evars->grad_1_chi[d_pp])*DENDRO_4 - DENDRO_2*deriv_evars->grad2_1_2_chi[d_pp]) + DENDRO_igt0*(3*pow((deriv_evars->grad_0_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_0_0_chi[d_pp]) + DENDRO_igt3*(3*pow((deriv_evars->grad_1_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_1_1_chi[d_pp]) + DENDRO_igt5*(3*pow((deriv_evars->grad_2_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_2_2_chi[d_pp])))/DENDRO_0;
// Dendro: reduced ops: 322
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 413 
// Dendro: printing temp variables
const double DENDRO_0 = pow(chi[pp], 2);
const double DENDRO_1 = (deriv_evars->grad_0_chi[d_pp])*(deriv_evars->grad_2_chi[d_pp]);
const double DENDRO_2 = 2*chi[pp];
const double DENDRO_3 = 3*(deriv_evars->grad_1_chi[d_pp]);
const double DENDRO_4 = 2*DENDRO_igt1;
const double DENDRO_5 = 2*DENDRO_igt2;
const double DENDRO_6 = 2*DENDRO_igt4;
const double DENDRO_7 = DENDRO_4*DENDRO_C2_k0_1 + DENDRO_5*DENDRO_C2_k0_2 + DENDRO_6*DENDRO_C2_k0_4 + DENDRO_C2_k0_0*DENDRO_igt0 + DENDRO_C2_k0_3*DENDRO_igt3 + DENDRO_C2_k0_5*DENDRO_igt5;
const double DENDRO_8 = DENDRO_4*DENDRO_C2_k1_1 + DENDRO_5*DENDRO_C2_k1_2 + DENDRO_6*DENDRO_C2_k1_4 + DENDRO_C2_k1_0*DENDRO_igt0 + DENDRO_C2_k1_3*DENDRO_igt3 + DENDRO_C2_k1_5*DENDRO_igt5;
const double DENDRO_9 = DENDRO_4*DENDRO_C2_k2_1 + DENDRO_5*DENDRO_C2_k2_2 + DENDRO_6*DENDRO_C2_k2_4 + DENDRO_C2_k2_0*DENDRO_igt0 + DENDRO_C2_k2_3*DENDRO_igt3 + DENDRO_C2_k2_5*DENDRO_igt5;
const double DENDRO_10 = 2.0*gt2[pp];
const double DENDRO_11 = 4*DENDRO_igt0;
const double DENDRO_12 = 4*DENDRO_igt2;
const double DENDRO_13 = 4*DENDRO_igt3;
const double DENDRO_14 = 4*DENDRO_igt5;
const double DENDRO_15 = 4*DENDRO_igt1;
const double DENDRO_16 = DENDRO_C1_k2_2*DENDRO_C2_k2_1 + DENDRO_C1_k2_4*DENDRO_C2_k2_0;
const double DENDRO_17 = DENDRO_C1_k2_2*DENDRO_C2_k2_2 + DENDRO_C1_k2_5*DENDRO_C2_k2_0;
const double DENDRO_18 = 4*DENDRO_igt4;
const double DENDRO_19 = DENDRO_C1_k2_4*DENDRO_C2_k2_2 + DENDRO_C1_k2_5*DENDRO_C2_k2_1;

// Dendro: printing variables
//--
DENDRO_RIJ2 = -1.0/4.0*(-DENDRO_0*((deriv_evars->grad_0_Gt0[d_pp])*DENDRO_10 + 2.0*(deriv_evars->grad_0_Gt1[d_pp])*gt4[pp] + 2.0*(deriv_evars->grad_0_Gt2[d_pp])*gt5[pp] + 2.0*(deriv_evars->grad_2_Gt0[d_pp])*gt0[pp] + 2.0*(deriv_evars->grad_2_Gt1[d_pp])*gt1[pp] + (deriv_evars->grad_2_Gt2[d_pp])*DENDRO_10 + DENDRO_11*(DENDRO_C1_k0_2*DENDRO_C2_k2_2 + 2*DENDRO_C1_k2_2*DENDRO_C2_k2_0) + DENDRO_11*(DENDRO_C1_k0_0*DENDRO_C2_k0_2 + DENDRO_C1_k0_2*DENDRO_C2_k0_0 + DENDRO_C1_k2_0*DENDRO_C2_k0_0) + DENDRO_11*(DENDRO_C1_k0_1*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_0 + DENDRO_C1_k2_1*DENDRO_C2_k1_0) + DENDRO_12*(DENDRO_17 + DENDRO_C1_k0_2*DENDRO_C2_k2_5) + DENDRO_12*(DENDRO_17 + DENDRO_C1_k0_5*DENDRO_C2_k2_2) + DENDRO_12*(2*DENDRO_C1_k0_2*DENDRO_C2_k0_2 + DENDRO_C1_k2_2*DENDRO_C2_k0_0) + DENDRO_12*(DENDRO_C1_k0_0*DENDRO_C2_k0_5 + DENDRO_C1_k0_5*DENDRO_C2_k0_0 + DENDRO_C1_k2_0*DENDRO_C2_k0_2) + DENDRO_12*(DENDRO_C1_k0_1*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_0 + DENDRO_C1_k2_1*DENDRO_C2_k1_2) + DENDRO_12*(DENDRO_C1_k0_4*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_2 + DENDRO_C1_k2_4*DENDRO_C2_k1_0) + DENDRO_13*(DENDRO_C1_k0_4*DENDRO_C2_k2_4 + 2*DENDRO_C1_k2_4*DENDRO_C2_k2_1) + DENDRO_13*(DENDRO_C1_k0_1*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_1 + DENDRO_C1_k2_1*DENDRO_C2_k0_1) + DENDRO_13*(DENDRO_C1_k0_3*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_1 + DENDRO_C1_k2_3*DENDRO_C2_k1_1) + DENDRO_14*(DENDRO_C1_k0_5*DENDRO_C2_k2_5 + 2*DENDRO_C1_k2_5*DENDRO_C2_k2_2) + DENDRO_14*(DENDRO_C1_k0_2*DENDRO_C2_k0_5 + DENDRO_C1_k0_5*DENDRO_C2_k0_2 + DENDRO_C1_k2_2*DENDRO_C2_k0_2) + DENDRO_14*(DENDRO_C1_k0_4*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_2 + DENDRO_C1_k2_4*DENDRO_C2_k1_2) + DENDRO_15*(DENDRO_16 + DENDRO_C1_k0_2*DENDRO_C2_k2_4) + DENDRO_15*(DENDRO_16 + DENDRO_C1_k0_4*DENDRO_C2_k2_2) + DENDRO_15*(DENDRO_C1_k0_0*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_0 + DENDRO_C1_k2_0*DENDRO_C2_k0_1) + DENDRO_15*(DENDRO_C1_k0_1*DENDRO_C2_k0_2 + DENDRO_C1_k0_2*DENDRO_C2_k0_1 + DENDRO_C1_k2_1*DENDRO_C2_k0_0) + DENDRO_15*(DENDRO_C1_k0_1*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_0 + DENDRO_C1_k2_1*DENDRO_C2_k1_1) + DENDRO_15*(DENDRO_C1_k0_3*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_1 + DENDRO_C1_k2_3*DENDRO_C2_k1_0) + DENDRO_18*(DENDRO_19 + DENDRO_C1_k0_4*DENDRO_C2_k2_5) + DENDRO_18*(DENDRO_19 + DENDRO_C1_k0_5*DENDRO_C2_k2_4) + DENDRO_18*(DENDRO_C1_k0_1*DENDRO_C2_k0_5 + DENDRO_C1_k0_5*DENDRO_C2_k0_1 + DENDRO_C1_k2_1*DENDRO_C2_k0_2) + DENDRO_18*(DENDRO_C1_k0_2*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_2 + DENDRO_C1_k2_2*DENDRO_C2_k0_1) + DENDRO_18*(DENDRO_C1_k0_3*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_1 + DENDRO_C1_k2_3*DENDRO_C2_k1_2) + DENDRO_18*(DENDRO_C1_k0_4*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_2 + DENDRO_C1_k2_4*DENDRO_C2_k1_1) + 2.0*DENDRO_7*(DENDRO_C1_k0_2 + DENDRO_C1_k2_0) + 2.0*DENDRO_8*(DENDRO_C1_k0_4 + DENDRO_C1_k2_1) + 2.0*DENDRO_9*(DENDRO_C1_k0_5 + DENDRO_C1_k2_2) - 2.0*DENDRO_igt0*deriv_evars->grad2_0_0_gt2[d_pp] - 4.0*DENDRO_igt1*deriv_evars->grad2_0_1_gt2[d_pp] - 4.0*DENDRO_igt2*deriv_evars->grad2_0_2_gt2[d_pp] - 2.0*DENDRO_igt3*deriv_evars->grad2_1_1_gt2[d_pp] - 4.0*DENDRO_igt4*deriv_evars->grad2_1_2_gt2[d_pp] - 2.0*DENDRO_igt5*deriv_evars->grad2_2_2_gt2[d_pp]) + DENDRO_1 + DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_C2_k0_2 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_C2_k1_2 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_C2_k2_2 - deriv_evars->grad2_0_2_chi[d_pp]) + gt2[pp]*(DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_7 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_8 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_9) + DENDRO_4*((deriv_evars->grad_0_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_1_chi[d_pp]) + DENDRO_5*(3*DENDRO_1 - DENDRO_2*deriv_evars->grad2_0_2_chi[d_pp]) + DENDRO_6*((deriv_evars->grad_2_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_1_2_chi[d_pp]) + DENDRO_igt0*(3*pow((deriv_evars->grad_0_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_0_0_chi[d_pp]) + DENDRO_igt3*(3*pow((deriv_evars->grad_1_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_1_1_chi[d_pp]) + DENDRO_igt5*(3*pow((deriv_evars->grad_2_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_2_2_chi[d_pp])))/DENDRO_0;
// Dendro: reduced ops: 322
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 363 
// Dendro: printing temp variables
const double DENDRO_0 = pow(chi[pp], 2);
const double DENDRO_1 = pow((deriv_evars->grad_1_chi[d_pp]), 2);
const double DENDRO_2 = 2*chi[pp];
const double DENDRO_3 = 3*(deriv_evars->grad_0_chi[d_pp]);
const double DENDRO_4 = 2*DENDRO_igt1;
const double DENDRO_5 = 2*DENDRO_igt2;
const double DENDRO_6 = 2*DENDRO_igt4;
const double DENDRO_7 = DENDRO_C2_k0_3*DENDRO_igt3;
const double DENDRO_8 = DENDRO_4*DENDRO_C2_k0_1 + DENDRO_5*DENDRO_C2_k0_2 + DENDRO_6*DENDRO_C2_k0_4 + DENDRO_7 + DENDRO_C2_k0_0*DENDRO_igt0 + DENDRO_C2_k0_5*DENDRO_igt5;
const double DENDRO_9 = DENDRO_C2_k1_3*DENDRO_igt3;
const double DENDRO_10 = DENDRO_4*DENDRO_C2_k1_1 + DENDRO_5*DENDRO_C2_k1_2 + DENDRO_6*DENDRO_C2_k1_4 + DENDRO_9 + DENDRO_C2_k1_0*DENDRO_igt0 + DENDRO_C2_k1_5*DENDRO_igt5;
const double DENDRO_11 = DENDRO_C2_k2_3*DENDRO_igt3;
const double DENDRO_12 = DENDRO_11 + DENDRO_4*DENDRO_C2_k2_1 + DENDRO_5*DENDRO_C2_k2_2 + DENDRO_6*DENDRO_C2_k2_4 + DENDRO_C2_k2_0*DENDRO_igt0 + DENDRO_C2_k2_5*DENDRO_igt5;
const double DENDRO_13 = 2*DENDRO_C1_k1_0;
const double DENDRO_14 = 4*DENDRO_igt0;
const double DENDRO_15 = 2*DENDRO_C1_k1_1;
const double DENDRO_16 = 2*DENDRO_C1_k1_2;
const double DENDRO_17 = 4*DENDRO_igt5;
const double DENDRO_18 = 2*DENDRO_C1_k1_4;
const double DENDRO_19 = 2*DENDRO_C1_k1_5;
const double DENDRO_20 = 4*DENDRO_igt1;
const double DENDRO_21 = DENDRO_C1_k1_1*DENDRO_C2_k1_3;
const double DENDRO_22 = DENDRO_C1_k1_3*DENDRO_C2_k1_1;
const double DENDRO_23 = 4*DENDRO_igt2;
const double DENDRO_24 = DENDRO_C1_k1_1*DENDRO_C2_k1_4;
const double DENDRO_25 = DENDRO_C1_k1_4*DENDRO_C2_k1_1;
const double DENDRO_26 = 4*DENDRO_igt4;
const double DENDRO_27 = DENDRO_C1_k1_3*DENDRO_C2_k1_4;
const double DENDRO_28 = DENDRO_C1_k1_4*DENDRO_C2_k1_3;

// Dendro: printing variables
//--
DENDRO_RIJ3 = -1.0/4.0*(-DENDRO_0*(4.0*(deriv_evars->grad_1_Gt0[d_pp])*gt1[pp] + 4.0*(deriv_evars->grad_1_Gt1[d_pp])*gt3[pp] + 4.0*(deriv_evars->grad_1_Gt2[d_pp])*gt4[pp] + 4.0*DENDRO_10*DENDRO_C1_k1_3 + 4*DENDRO_11*(DENDRO_18 + DENDRO_C1_k2_3) + 4.0*DENDRO_12*DENDRO_C1_k1_4 + DENDRO_14*DENDRO_C2_k0_1*(DENDRO_13 + DENDRO_C1_k0_1) + DENDRO_14*DENDRO_C2_k2_1*(DENDRO_16 + DENDRO_C1_k2_1) + DENDRO_17*DENDRO_C2_k0_4*(DENDRO_16 + DENDRO_C1_k0_4) + DENDRO_17*DENDRO_C2_k2_4*(DENDRO_19 + DENDRO_C1_k2_4) + DENDRO_20*(DENDRO_21 + 2*DENDRO_22) + DENDRO_20*(2*DENDRO_21 + DENDRO_22) + DENDRO_20*(DENDRO_13*DENDRO_C2_k0_3 + DENDRO_C1_k0_3*DENDRO_C2_k0_1) + DENDRO_20*(DENDRO_15*DENDRO_C2_k0_1 + DENDRO_C1_k0_1*DENDRO_C2_k0_3) + DENDRO_20*(DENDRO_16*DENDRO_C2_k2_3 + DENDRO_C1_k2_3*DENDRO_C2_k2_1) + DENDRO_20*(DENDRO_18*DENDRO_C2_k2_1 + DENDRO_C1_k2_1*DENDRO_C2_k2_3) + DENDRO_23*(DENDRO_24 + 2*DENDRO_25) + DENDRO_23*(2*DENDRO_24 + DENDRO_25) + DENDRO_23*(DENDRO_13*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_1) + DENDRO_23*(DENDRO_16*DENDRO_C2_k0_1 + DENDRO_C1_k0_1*DENDRO_C2_k0_4) + DENDRO_23*(DENDRO_16*DENDRO_C2_k2_4 + DENDRO_C1_k2_4*DENDRO_C2_k2_1) + DENDRO_23*(DENDRO_19*DENDRO_C2_k2_1 + DENDRO_C1_k2_1*DENDRO_C2_k2_4) + DENDRO_26*(DENDRO_27 + 2*DENDRO_28) + DENDRO_26*(2*DENDRO_27 + DENDRO_28) + DENDRO_26*(DENDRO_15*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_3) + DENDRO_26*(DENDRO_16*DENDRO_C2_k0_3 + DENDRO_C1_k0_3*DENDRO_C2_k0_4) + DENDRO_26*(DENDRO_18*DENDRO_C2_k2_4 + DENDRO_C1_k2_4*DENDRO_C2_k2_3) + DENDRO_26*(DENDRO_19*DENDRO_C2_k2_3 + DENDRO_C1_k2_3*DENDRO_C2_k2_4) + 4*DENDRO_7*(DENDRO_15 + DENDRO_C1_k0_3) + 4.0*DENDRO_8*DENDRO_C1_k1_1 + 12*DENDRO_9*DENDRO_C1_k1_3 + 12*DENDRO_C1_k1_1*DENDRO_C2_k1_1*DENDRO_igt0 + 12*DENDRO_C1_k1_4*DENDRO_C2_k1_4*DENDRO_igt5 - 2.0*DENDRO_igt0*deriv_evars->grad2_0_0_gt3[d_pp] - 4.0*DENDRO_igt1*deriv_evars->grad2_0_1_gt3[d_pp] - 4.0*DENDRO_igt2*deriv_evars->grad2_0_2_gt3[d_pp] - 2.0*DENDRO_igt3*deriv_evars->grad2_1_1_gt3[d_pp] - 4.0*DENDRO_igt4*deriv_evars->grad2_1_2_gt3[d_pp] - 2.0*DENDRO_igt5*deriv_evars->grad2_2_2_gt3[d_pp]) + DENDRO_1 + DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_C2_k0_3 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_C2_k1_3 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_C2_k2_3 - deriv_evars->grad2_1_1_chi[d_pp]) + gt3[pp]*(DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_8 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_10 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_12) + DENDRO_4*((deriv_evars->grad_1_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_1_chi[d_pp]) + DENDRO_5*((deriv_evars->grad_2_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_2_chi[d_pp]) + DENDRO_6*(3*(deriv_evars->grad_1_chi[d_pp])*(deriv_evars->grad_2_chi[d_pp]) - DENDRO_2*deriv_evars->grad2_1_2_chi[d_pp]) + DENDRO_igt0*(3*pow((deriv_evars->grad_0_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_0_0_chi[d_pp]) + DENDRO_igt3*(3*DENDRO_1 - DENDRO_2*deriv_evars->grad2_1_1_chi[d_pp]) + DENDRO_igt5*(3*pow((deriv_evars->grad_2_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_2_2_chi[d_pp])))/DENDRO_0;
// Dendro: reduced ops: 264
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 413 
// Dendro: printing temp variables
const double DENDRO_0 = pow(chi[pp], 2);
const double DENDRO_1 = (deriv_evars->grad_1_chi[d_pp])*(deriv_evars->grad_2_chi[d_pp]);
const double DENDRO_2 = 2*chi[pp];
const double DENDRO_3 = 3*(deriv_evars->grad_0_chi[d_pp]);
const double DENDRO_4 = 2*DENDRO_igt1;
const double DENDRO_5 = 2*DENDRO_igt2;
const double DENDRO_6 = 2*DENDRO_igt4;
const double DENDRO_7 = DENDRO_4*DENDRO_C2_k0_1 + DENDRO_5*DENDRO_C2_k0_2 + DENDRO_6*DENDRO_C2_k0_4 + DENDRO_C2_k0_0*DENDRO_igt0 + DENDRO_C2_k0_3*DENDRO_igt3 + DENDRO_C2_k0_5*DENDRO_igt5;
const double DENDRO_8 = DENDRO_4*DENDRO_C2_k1_1 + DENDRO_5*DENDRO_C2_k1_2 + DENDRO_6*DENDRO_C2_k1_4 + DENDRO_C2_k1_0*DENDRO_igt0 + DENDRO_C2_k1_3*DENDRO_igt3 + DENDRO_C2_k1_5*DENDRO_igt5;
const double DENDRO_9 = DENDRO_4*DENDRO_C2_k2_1 + DENDRO_5*DENDRO_C2_k2_2 + DENDRO_6*DENDRO_C2_k2_4 + DENDRO_C2_k2_0*DENDRO_igt0 + DENDRO_C2_k2_3*DENDRO_igt3 + DENDRO_C2_k2_5*DENDRO_igt5;
const double DENDRO_10 = 2.0*gt4[pp];
const double DENDRO_11 = 4*DENDRO_igt0;
const double DENDRO_12 = 4*DENDRO_igt3;
const double DENDRO_13 = 4*DENDRO_igt4;
const double DENDRO_14 = 4*DENDRO_igt5;
const double DENDRO_15 = 4*DENDRO_igt1;
const double DENDRO_16 = DENDRO_C1_k2_2*DENDRO_C2_k2_3 + DENDRO_C1_k2_4*DENDRO_C2_k2_1;
const double DENDRO_17 = 4*DENDRO_igt2;
const double DENDRO_18 = DENDRO_C1_k2_2*DENDRO_C2_k2_4 + DENDRO_C1_k2_5*DENDRO_C2_k2_1;
const double DENDRO_19 = DENDRO_C1_k2_4*DENDRO_C2_k2_4 + DENDRO_C1_k2_5*DENDRO_C2_k2_3;

// Dendro: printing variables
//--
DENDRO_RIJ4 = -1.0/4.0*(-DENDRO_0*(2.0*(deriv_evars->grad_1_Gt0[d_pp])*gt2[pp] + (deriv_evars->grad_1_Gt1[d_pp])*DENDRO_10 + 2.0*(deriv_evars->grad_1_Gt2[d_pp])*gt5[pp] + 2.0*(deriv_evars->grad_2_Gt0[d_pp])*gt1[pp] + 2.0*(deriv_evars->grad_2_Gt1[d_pp])*gt3[pp] + (deriv_evars->grad_2_Gt2[d_pp])*DENDRO_10 + DENDRO_11*(DENDRO_C1_k1_2*DENDRO_C2_k2_2 + 2*DENDRO_C1_k2_2*DENDRO_C2_k2_1) + DENDRO_11*(DENDRO_C1_k0_2*DENDRO_C2_k0_1 + DENDRO_C1_k1_0*DENDRO_C2_k0_2 + DENDRO_C1_k2_0*DENDRO_C2_k0_1) + DENDRO_11*(DENDRO_C1_k1_1*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_1 + DENDRO_C1_k2_1*DENDRO_C2_k1_1) + DENDRO_12*(DENDRO_C1_k1_4*DENDRO_C2_k2_4 + 2*DENDRO_C1_k2_4*DENDRO_C2_k2_3) + DENDRO_12*(DENDRO_C1_k0_4*DENDRO_C2_k0_3 + DENDRO_C1_k1_1*DENDRO_C2_k0_4 + DENDRO_C1_k2_1*DENDRO_C2_k0_3) + DENDRO_12*(DENDRO_C1_k1_3*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_3 + DENDRO_C1_k2_3*DENDRO_C2_k1_3) + DENDRO_13*(DENDRO_19 + DENDRO_C1_k1_4*DENDRO_C2_k2_5) + DENDRO_13*(DENDRO_19 + DENDRO_C1_k1_5*DENDRO_C2_k2_4) + DENDRO_13*(2*DENDRO_C1_k1_4*DENDRO_C2_k1_4 + DENDRO_C1_k2_4*DENDRO_C2_k1_3) + DENDRO_13*(DENDRO_C1_k0_4*DENDRO_C2_k0_4 + DENDRO_C1_k1_2*DENDRO_C2_k0_4 + DENDRO_C1_k2_2*DENDRO_C2_k0_3) + DENDRO_13*(DENDRO_C1_k0_5*DENDRO_C2_k0_3 + DENDRO_C1_k1_1*DENDRO_C2_k0_5 + DENDRO_C1_k2_1*DENDRO_C2_k0_4) + DENDRO_13*(DENDRO_C1_k1_3*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_3 + DENDRO_C1_k2_3*DENDRO_C2_k1_4) + DENDRO_14*(DENDRO_C1_k1_5*DENDRO_C2_k2_5 + 2*DENDRO_C1_k2_5*DENDRO_C2_k2_4) + DENDRO_14*(DENDRO_C1_k0_5*DENDRO_C2_k0_4 + DENDRO_C1_k1_2*DENDRO_C2_k0_5 + DENDRO_C1_k2_2*DENDRO_C2_k0_4) + DENDRO_14*(DENDRO_C1_k1_4*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_4 + DENDRO_C1_k2_4*DENDRO_C2_k1_4) + DENDRO_15*(DENDRO_16 + DENDRO_C1_k1_2*DENDRO_C2_k2_4) + DENDRO_15*(DENDRO_16 + DENDRO_C1_k1_4*DENDRO_C2_k2_2) + DENDRO_15*(DENDRO_C1_k0_2*DENDRO_C2_k0_3 + DENDRO_C1_k1_1*DENDRO_C2_k0_2 + DENDRO_C1_k2_1*DENDRO_C2_k0_1) + DENDRO_15*(DENDRO_C1_k0_4*DENDRO_C2_k0_1 + DENDRO_C1_k1_0*DENDRO_C2_k0_4 + DENDRO_C1_k2_0*DENDRO_C2_k0_3) + DENDRO_15*(DENDRO_C1_k1_1*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_1 + DENDRO_C1_k2_1*DENDRO_C2_k1_3) + DENDRO_15*(DENDRO_C1_k1_2*DENDRO_C2_k1_3 + DENDRO_C1_k1_3*DENDRO_C2_k1_2 + DENDRO_C1_k2_3*DENDRO_C2_k1_1) + DENDRO_17*(DENDRO_18 + DENDRO_C1_k1_2*DENDRO_C2_k2_5) + DENDRO_17*(DENDRO_18 + DENDRO_C1_k1_5*DENDRO_C2_k2_2) + DENDRO_17*(DENDRO_C1_k0_2*DENDRO_C2_k0_4 + DENDRO_C1_k1_2*DENDRO_C2_k0_2 + DENDRO_C1_k2_2*DENDRO_C2_k0_1) + DENDRO_17*(DENDRO_C1_k0_5*DENDRO_C2_k0_1 + DENDRO_C1_k1_0*DENDRO_C2_k0_5 + DENDRO_C1_k2_0*DENDRO_C2_k0_4) + DENDRO_17*(DENDRO_C1_k1_1*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_1 + DENDRO_C1_k2_1*DENDRO_C2_k1_4) + DENDRO_17*(DENDRO_C1_k1_2*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_2 + DENDRO_C1_k2_4*DENDRO_C2_k1_1) + 2.0*DENDRO_7*(DENDRO_C1_k1_2 + DENDRO_C1_k2_1) + 2.0*DENDRO_8*(DENDRO_C1_k1_4 + DENDRO_C1_k2_3) + 2.0*DENDRO_9*(DENDRO_C1_k1_5 + DENDRO_C1_k2_4) - 2.0*DENDRO_igt0*deriv_evars->grad2_0_0_gt4[d_pp] - 4.0*DENDRO_igt1*deriv_evars->grad2_0_1_gt4[d_pp] - 4.0*DENDRO_igt2*deriv_evars->grad2_0_2_gt4[d_pp] - 2.0*DENDRO_igt3*deriv_evars->grad2_1_1_gt4[d_pp] - 4.0*DENDRO_igt4*deriv_evars->grad2_1_2_gt4[d_pp] - 2.0*DENDRO_igt5*deriv_evars->grad2_2_2_gt4[d_pp]) + DENDRO_1 + DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_C2_k0_4 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_C2_k1_4 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_C2_k2_4 - deriv_evars->grad2_1_2_chi[d_pp]) + gt4[pp]*(DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_7 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_8 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_9) + DENDRO_4*((deriv_evars->grad_1_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_1_chi[d_pp]) + DENDRO_5*((deriv_evars->grad_2_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_2_chi[d_pp]) + DENDRO_6*(3*DENDRO_1 - DENDRO_2*deriv_evars->grad2_1_2_chi[d_pp]) + DENDRO_igt0*(3*pow((deriv_evars->grad_0_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_0_0_chi[d_pp]) + DENDRO_igt3*(3*pow((deriv_evars->grad_1_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_1_1_chi[d_pp]) + DENDRO_igt5*(3*pow((deriv_evars->grad_2_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_2_2_chi[d_pp])))/DENDRO_0;
// Dendro: reduced ops: 322
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 363 
// Dendro: printing temp variables
const double DENDRO_0 = pow(chi[pp], 2);
const double DENDRO_1 = pow((deriv_evars->grad_2_chi[d_pp]), 2);
const double DENDRO_2 = 2*chi[pp];
const double DENDRO_3 = 3*(deriv_evars->grad_0_chi[d_pp]);
const double DENDRO_4 = 2*DENDRO_igt1;
const double DENDRO_5 = 2*DENDRO_igt2;
const double DENDRO_6 = 2*DENDRO_igt4;
const double DENDRO_7 = DENDRO_C2_k0_5*DENDRO_igt5;
const double DENDRO_8 = DENDRO_4*DENDRO_C2_k0_1 + DENDRO_5*DENDRO_C2_k0_2 + DENDRO_6*DENDRO_C2_k0_4 + DENDRO_7 + DENDRO_C2_k0_0*DENDRO_igt0 + DENDRO_C2_k0_3*DENDRO_igt3;
const double DENDRO_9 = DENDRO_C2_k1_5*DENDRO_igt5;
const double DENDRO_10 = DENDRO_4*DENDRO_C2_k1_1 + DENDRO_5*DENDRO_C2_k1_2 + DENDRO_6*DENDRO_C2_k1_4 + DENDRO_9 + DENDRO_C2_k1_0*DENDRO_igt0 + DENDRO_C2_k1_3*DENDRO_igt3;
const double DENDRO_11 = DENDRO_C2_k2_5*DENDRO_igt5;
const double DENDRO_12 = DENDRO_11 + DENDRO_4*DENDRO_C2_k2_1 + DENDRO_5*DENDRO_C2_k2_2 + DENDRO_6*DENDRO_C2_k2_4 + DENDRO_C2_k2_0*DENDRO_igt0 + DENDRO_C2_k2_3*DENDRO_igt3;
const double DENDRO_13 = 2*DENDRO_C1_k2_0;
const double DENDRO_14 = 4*DENDRO_igt0;
const double DENDRO_15 = 2*DENDRO_C1_k2_1;
const double DENDRO_16 = 4*DENDRO_igt3;
const double DENDRO_17 = 2*DENDRO_C1_k2_2;
const double DENDRO_18 = 2*DENDRO_C1_k2_3;
const double DENDRO_19 = 2*DENDRO_C1_k2_4;
const double DENDRO_20 = 4*DENDRO_igt1;
const double DENDRO_21 = DENDRO_C1_k2_2*DENDRO_C2_k2_4;
const double DENDRO_22 = DENDRO_C1_k2_4*DENDRO_C2_k2_2;
const double DENDRO_23 = 4*DENDRO_igt2;
const double DENDRO_24 = DENDRO_C1_k2_2*DENDRO_C2_k2_5;
const double DENDRO_25 = DENDRO_C1_k2_5*DENDRO_C2_k2_2;
const double DENDRO_26 = 4*DENDRO_igt4;
const double DENDRO_27 = DENDRO_C1_k2_4*DENDRO_C2_k2_5;
const double DENDRO_28 = DENDRO_C1_k2_5*DENDRO_C2_k2_4;

// Dendro: printing variables
//--
DENDRO_RIJ5 = -1.0/4.0*(-DENDRO_0*(4.0*(deriv_evars->grad_2_Gt0[d_pp])*gt2[pp] + 4.0*(deriv_evars->grad_2_Gt1[d_pp])*gt4[pp] + 4.0*(deriv_evars->grad_2_Gt2[d_pp])*gt5[pp] + 4.0*DENDRO_10*DENDRO_C1_k2_4 + 12*DENDRO_11*DENDRO_C1_k2_5 + 4.0*DENDRO_12*DENDRO_C1_k2_5 + DENDRO_14*DENDRO_C2_k0_2*(DENDRO_13 + DENDRO_C1_k0_2) + DENDRO_14*DENDRO_C2_k1_2*(DENDRO_15 + DENDRO_C1_k1_2) + DENDRO_16*DENDRO_C2_k0_4*(DENDRO_15 + DENDRO_C1_k0_4) + DENDRO_16*DENDRO_C2_k1_4*(DENDRO_18 + DENDRO_C1_k1_4) + DENDRO_20*(DENDRO_21 + 2*DENDRO_22) + DENDRO_20*(2*DENDRO_21 + DENDRO_22) + DENDRO_20*(DENDRO_13*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_2) + DENDRO_20*(DENDRO_15*DENDRO_C2_k0_2 + DENDRO_C1_k0_2*DENDRO_C2_k0_4) + DENDRO_20*(DENDRO_15*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_2) + DENDRO_20*(DENDRO_18*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_4) + DENDRO_23*(DENDRO_24 + 2*DENDRO_25) + DENDRO_23*(2*DENDRO_24 + DENDRO_25) + DENDRO_23*(DENDRO_13*DENDRO_C2_k0_5 + DENDRO_C1_k0_5*DENDRO_C2_k0_2) + DENDRO_23*(DENDRO_15*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_2) + DENDRO_23*(DENDRO_17*DENDRO_C2_k0_2 + DENDRO_C1_k0_2*DENDRO_C2_k0_5) + DENDRO_23*(DENDRO_19*DENDRO_C2_k1_2 + DENDRO_C1_k1_2*DENDRO_C2_k1_5) + DENDRO_26*(DENDRO_27 + 2*DENDRO_28) + DENDRO_26*(2*DENDRO_27 + DENDRO_28) + DENDRO_26*(DENDRO_15*DENDRO_C2_k0_5 + DENDRO_C1_k0_5*DENDRO_C2_k0_4) + DENDRO_26*(DENDRO_17*DENDRO_C2_k0_4 + DENDRO_C1_k0_4*DENDRO_C2_k0_5) + DENDRO_26*(DENDRO_18*DENDRO_C2_k1_5 + DENDRO_C1_k1_5*DENDRO_C2_k1_4) + DENDRO_26*(DENDRO_19*DENDRO_C2_k1_4 + DENDRO_C1_k1_4*DENDRO_C2_k1_5) + 4*DENDRO_7*(DENDRO_17 + DENDRO_C1_k0_5) + 4.0*DENDRO_8*DENDRO_C1_k2_2 + 4*DENDRO_9*(DENDRO_19 + DENDRO_C1_k1_5) + 12*DENDRO_C1_k2_2*DENDRO_C2_k2_2*DENDRO_igt0 + 12*DENDRO_C1_k2_4*DENDRO_C2_k2_4*DENDRO_igt3 - 2.0*DENDRO_igt0*deriv_evars->grad2_0_0_gt5[d_pp] - 4.0*DENDRO_igt1*deriv_evars->grad2_0_1_gt5[d_pp] - 4.0*DENDRO_igt2*deriv_evars->grad2_0_2_gt5[d_pp] - 2.0*DENDRO_igt3*deriv_evars->grad2_1_1_gt5[d_pp] - 4.0*DENDRO_igt4*deriv_evars->grad2_1_2_gt5[d_pp] - 2.0*DENDRO_igt5*deriv_evars->grad2_2_2_gt5[d_pp]) + DENDRO_1 + DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_C2_k0_5 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_C2_k1_5 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_C2_k2_5 - deriv_evars->grad2_2_2_chi[d_pp]) + gt5[pp]*(DENDRO_2*((deriv_evars->grad_0_chi[d_pp])*DENDRO_8 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_10 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_12) + DENDRO_4*((deriv_evars->grad_1_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_1_chi[d_pp]) + DENDRO_5*((deriv_evars->grad_2_chi[d_pp])*DENDRO_3 - DENDRO_2*deriv_evars->grad2_0_2_chi[d_pp]) + DENDRO_6*(3*(deriv_evars->grad_1_chi[d_pp])*(deriv_evars->grad_2_chi[d_pp]) - DENDRO_2*deriv_evars->grad2_1_2_chi[d_pp]) + DENDRO_igt0*(3*pow((deriv_evars->grad_0_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_0_0_chi[d_pp]) + DENDRO_igt3*(3*pow((deriv_evars->grad_1_chi[d_pp]), 2) - DENDRO_2*deriv_evars->grad2_1_1_chi[d_pp]) + DENDRO_igt5*(3*DENDRO_1 - DENDRO_2*deriv_evars->grad2_2_2_chi[d_pp])))/DENDRO_0;
// Dendro: reduced ops: 264
// Dendro: }}} 

}
{
// Dendro: {{{ 
// Dendro: original ops: 42 
// Dendro: printing temp variables
const double DENDRO_0 = 2*DENDRO_igt1;
const double DENDRO_1 = 2*DENDRO_igt2;
const double DENDRO_2 = 2*DENDRO_igt4;

// Dendro: printing variables
//--
DENDRO_Gtk0 = DENDRO_0*DENDRO_C2_k0_1 + DENDRO_1*DENDRO_C2_k0_2 + DENDRO_2*DENDRO_C2_k0_4 + DENDRO_C2_k0_0*DENDRO_igt0 + DENDRO_C2_k0_3*DENDRO_igt3 + DENDRO_C2_k0_5*DENDRO_igt5;
//--
DENDRO_Gtk1 = DENDRO_0*DENDRO_C2_k1_1 + DENDRO_1*DENDRO_C2_k1_2 + DENDRO_2*DENDRO_C2_k1_4 + DENDRO_C2_k1_0*DENDRO_igt0 + DENDRO_C2_k1_3*DENDRO_igt3 + DENDRO_C2_k1_5*DENDRO_igt5;
//--
DENDRO_Gtk2 = DENDRO_0*DENDRO_C2_k2_1 + DENDRO_1*DENDRO_C2_k2_2 + DENDRO_2*DENDRO_C2_k2_4 + DENDRO_C2_k2_0*DENDRO_igt0 + DENDRO_C2_k2_3*DENDRO_igt3 + DENDRO_C2_k2_5*DENDRO_igt5;
// Dendro: reduced ops: 36
// Dendro: }}} 

}

{
// Dendro: {{{ 
// Dendro: original ops: 125 
// Dendro: printing temp variables
const double DENDRO_0 = 2*At1[pp];
const double DENDRO_1 = 2*At2[pp];
const double DENDRO_2 = (2.0/3.0)*At0[pp];
const double DENDRO_3 = (deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_0;
const double DENDRO_4 = (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_0;
const double DENDRO_5 = (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_0;
const double DENDRO_6 = DENDRO_RIJ0*alpha[pp];

// Dendro: printing variables
//--
At_rhs00[pp] = (deriv_evars->grad_0_At0[d_pp])*beta0[pp] + (4.0/3.0)*(deriv_evars->grad_0_beta0[d_pp])*At0[pp] + (deriv_evars->grad_0_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta2[d_pp])*DENDRO_1 + (deriv_evars->grad_1_At0[d_pp])*beta1[pp] - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_2 + (deriv_evars->grad_2_At0[d_pp])*beta2[pp] - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_2 - alpha[pp]*(-At0[pp]*K[pp] + 2*At0[pp]*(At0[pp]*DENDRO_igt0 + At1[pp]*DENDRO_igt1 + At2[pp]*DENDRO_igt2) + DENDRO_0*(At0[pp]*DENDRO_igt1 + At1[pp]*DENDRO_igt3 + At2[pp]*DENDRO_igt4) + DENDRO_1*(At0[pp]*DENDRO_igt2 + At1[pp]*DENDRO_igt4 + At2[pp]*DENDRO_igt5)) + (1.0/3.0)*chi[pp]*(3*DENDRO_3 + 3*DENDRO_4 + 3*DENDRO_5 + 3*DENDRO_6 - 3*deriv_evars->grad2_0_0_alpha[d_pp] - gt0[pp]*(DENDRO_igt0*(DENDRO_3 + DENDRO_4 + DENDRO_5 + DENDRO_6 - deriv_evars->grad2_0_0_alpha[d_pp]) + 2*DENDRO_igt1*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_1 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_1 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_1 + DENDRO_RIJ1*alpha[pp] - deriv_evars->grad2_0_1_alpha[d_pp]) + 2*DENDRO_igt2*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_2 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_2 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_2 + DENDRO_RIJ2*alpha[pp] - deriv_evars->grad2_0_2_alpha[d_pp]) + DENDRO_igt3*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_3 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_3 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_3 + DENDRO_RIJ3*alpha[pp] - deriv_evars->grad2_1_1_alpha[d_pp]) + 2*DENDRO_igt4*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_4 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_4 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_4 + DENDRO_RIJ4*alpha[pp] - deriv_evars->grad2_1_2_alpha[d_pp]) + DENDRO_igt5*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_5 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_5 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_5 + DENDRO_RIJ5*alpha[pp] - deriv_evars->grad2_2_2_alpha[d_pp])));
// Dendro: reduced ops: 122
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 125 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*At1[pp];
const double DENDRO_1 = (deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_1;
const double DENDRO_2 = (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_1;
const double DENDRO_3 = (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_1;
const double DENDRO_4 = DENDRO_RIJ1*alpha[pp];

// Dendro: printing variables
//--
At_rhs01[pp] = (deriv_evars->grad_0_At1[d_pp])*beta0[pp] + (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta1[d_pp])*At3[pp] + (deriv_evars->grad_0_beta2[d_pp])*At4[pp] + (deriv_evars->grad_1_At1[d_pp])*beta1[pp] + (deriv_evars->grad_1_beta0[d_pp])*At0[pp] + (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_beta2[d_pp])*At2[pp] + (deriv_evars->grad_2_At1[d_pp])*beta2[pp] - 2.0/3.0*(deriv_evars->grad_2_beta2[d_pp])*At1[pp] - alpha[pp]*(2*At0[pp]*(At1[pp]*DENDRO_igt0 + At3[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt2) - At1[pp]*K[pp] + 2*At1[pp]*(At1[pp]*DENDRO_igt1 + At3[pp]*DENDRO_igt3 + At4[pp]*DENDRO_igt4) + 2*At2[pp]*(At1[pp]*DENDRO_igt2 + At3[pp]*DENDRO_igt4 + At4[pp]*DENDRO_igt5)) + (1.0/3.0)*chi[pp]*(3*DENDRO_1 + 3*DENDRO_2 + 3*DENDRO_3 + 3*DENDRO_4 - 3*deriv_evars->grad2_0_1_alpha[d_pp] - gt1[pp]*(DENDRO_igt0*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_0 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_0 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_0 + DENDRO_RIJ0*alpha[pp] - deriv_evars->grad2_0_0_alpha[d_pp]) + 2*DENDRO_igt1*(DENDRO_1 + DENDRO_2 + DENDRO_3 + DENDRO_4 - deriv_evars->grad2_0_1_alpha[d_pp]) + 2*DENDRO_igt2*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_2 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_2 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_2 + DENDRO_RIJ2*alpha[pp] - deriv_evars->grad2_0_2_alpha[d_pp]) + DENDRO_igt3*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_3 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_3 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_3 + DENDRO_RIJ3*alpha[pp] - deriv_evars->grad2_1_1_alpha[d_pp]) + 2*DENDRO_igt4*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_4 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_4 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_4 + DENDRO_RIJ4*alpha[pp] - deriv_evars->grad2_1_2_alpha[d_pp]) + DENDRO_igt5*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_5 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_5 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_5 + DENDRO_RIJ5*alpha[pp] - deriv_evars->grad2_2_2_alpha[d_pp])));
// Dendro: reduced ops: 125
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 125 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*At2[pp];
const double DENDRO_1 = (deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_2;
const double DENDRO_2 = (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_2;
const double DENDRO_3 = (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_2;
const double DENDRO_4 = DENDRO_RIJ2*alpha[pp];

// Dendro: printing variables
//--
At_rhs02[pp] = (deriv_evars->grad_0_At2[d_pp])*beta0[pp] + (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta1[d_pp])*At4[pp] + (deriv_evars->grad_0_beta2[d_pp])*At5[pp] + (deriv_evars->grad_1_At2[d_pp])*beta1[pp] - 2.0/3.0*(deriv_evars->grad_1_beta1[d_pp])*At2[pp] + (deriv_evars->grad_2_At2[d_pp])*beta2[pp] + (deriv_evars->grad_2_beta0[d_pp])*At0[pp] + (deriv_evars->grad_2_beta1[d_pp])*At1[pp] + (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 - alpha[pp]*(2*At0[pp]*(At2[pp]*DENDRO_igt0 + At4[pp]*DENDRO_igt1 + At5[pp]*DENDRO_igt2) + 2*At1[pp]*(At2[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt3 + At5[pp]*DENDRO_igt4) - At2[pp]*K[pp] + 2*At2[pp]*(At2[pp]*DENDRO_igt2 + At4[pp]*DENDRO_igt4 + At5[pp]*DENDRO_igt5)) + (1.0/3.0)*chi[pp]*(3*DENDRO_1 + 3*DENDRO_2 + 3*DENDRO_3 + 3*DENDRO_4 - 3*deriv_evars->grad2_0_2_alpha[d_pp] - gt2[pp]*(DENDRO_igt0*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_0 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_0 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_0 + DENDRO_RIJ0*alpha[pp] - deriv_evars->grad2_0_0_alpha[d_pp]) + 2*DENDRO_igt1*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_1 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_1 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_1 + DENDRO_RIJ1*alpha[pp] - deriv_evars->grad2_0_1_alpha[d_pp]) + 2*DENDRO_igt2*(DENDRO_1 + DENDRO_2 + DENDRO_3 + DENDRO_4 - deriv_evars->grad2_0_2_alpha[d_pp]) + DENDRO_igt3*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_3 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_3 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_3 + DENDRO_RIJ3*alpha[pp] - deriv_evars->grad2_1_1_alpha[d_pp]) + 2*DENDRO_igt4*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_4 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_4 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_4 + DENDRO_RIJ4*alpha[pp] - deriv_evars->grad2_1_2_alpha[d_pp]) + DENDRO_igt5*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_5 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_5 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_5 + DENDRO_RIJ5*alpha[pp] - deriv_evars->grad2_2_2_alpha[d_pp])));
// Dendro: reduced ops: 125
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 125 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*At3[pp];
const double DENDRO_1 = 2*At1[pp];
const double DENDRO_2 = 2*At4[pp];
const double DENDRO_3 = (deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_3;
const double DENDRO_4 = (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_3;
const double DENDRO_5 = (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_3;
const double DENDRO_6 = DENDRO_RIJ3*alpha[pp];

// Dendro: printing variables
//--
At_rhs11[pp] = (deriv_evars->grad_0_At3[d_pp])*beta0[pp] - (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_1_At3[d_pp])*beta1[pp] + (deriv_evars->grad_1_beta0[d_pp])*DENDRO_1 + (4.0/3.0)*(deriv_evars->grad_1_beta1[d_pp])*At3[pp] + (deriv_evars->grad_1_beta2[d_pp])*DENDRO_2 + (deriv_evars->grad_2_At3[d_pp])*beta2[pp] - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 - alpha[pp]*(-At3[pp]*K[pp] + 2*At3[pp]*(At1[pp]*DENDRO_igt1 + At3[pp]*DENDRO_igt3 + At4[pp]*DENDRO_igt4) + DENDRO_1*(At1[pp]*DENDRO_igt0 + At3[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt2) + DENDRO_2*(At1[pp]*DENDRO_igt2 + At3[pp]*DENDRO_igt4 + At4[pp]*DENDRO_igt5)) + (1.0/3.0)*chi[pp]*(3*DENDRO_3 + 3*DENDRO_4 + 3*DENDRO_5 + 3*DENDRO_6 - 3*deriv_evars->grad2_1_1_alpha[d_pp] - gt3[pp]*(DENDRO_igt0*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_0 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_0 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_0 + DENDRO_RIJ0*alpha[pp] - deriv_evars->grad2_0_0_alpha[d_pp]) + 2*DENDRO_igt1*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_1 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_1 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_1 + DENDRO_RIJ1*alpha[pp] - deriv_evars->grad2_0_1_alpha[d_pp]) + 2*DENDRO_igt2*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_2 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_2 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_2 + DENDRO_RIJ2*alpha[pp] - deriv_evars->grad2_0_2_alpha[d_pp]) + DENDRO_igt3*(DENDRO_3 + DENDRO_4 + DENDRO_5 + DENDRO_6 - deriv_evars->grad2_1_1_alpha[d_pp]) + 2*DENDRO_igt4*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_4 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_4 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_4 + DENDRO_RIJ4*alpha[pp] - deriv_evars->grad2_1_2_alpha[d_pp]) + DENDRO_igt5*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_5 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_5 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_5 + DENDRO_RIJ5*alpha[pp] - deriv_evars->grad2_2_2_alpha[d_pp])));
// Dendro: reduced ops: 122
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 125 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*At4[pp];
const double DENDRO_1 = (deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_4;
const double DENDRO_2 = (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_4;
const double DENDRO_3 = (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_4;
const double DENDRO_4 = DENDRO_RIJ4*alpha[pp];

// Dendro: printing variables
//--
At_rhs12[pp] = (deriv_evars->grad_0_At4[d_pp])*beta0[pp] - 2.0/3.0*(deriv_evars->grad_0_beta0[d_pp])*At4[pp] + (deriv_evars->grad_1_At4[d_pp])*beta1[pp] + (deriv_evars->grad_1_beta0[d_pp])*At2[pp] + (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_beta2[d_pp])*At5[pp] + (deriv_evars->grad_2_At4[d_pp])*beta2[pp] + (deriv_evars->grad_2_beta0[d_pp])*At1[pp] + (deriv_evars->grad_2_beta1[d_pp])*At3[pp] + (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 - alpha[pp]*(2*At1[pp]*(At2[pp]*DENDRO_igt0 + At4[pp]*DENDRO_igt1 + At5[pp]*DENDRO_igt2) + 2*At3[pp]*(At2[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt3 + At5[pp]*DENDRO_igt4) - At4[pp]*K[pp] + 2*At4[pp]*(At2[pp]*DENDRO_igt2 + At4[pp]*DENDRO_igt4 + At5[pp]*DENDRO_igt5)) + (1.0/3.0)*chi[pp]*(3*DENDRO_1 + 3*DENDRO_2 + 3*DENDRO_3 + 3*DENDRO_4 - 3*deriv_evars->grad2_1_2_alpha[d_pp] - gt4[pp]*(DENDRO_igt0*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_0 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_0 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_0 + DENDRO_RIJ0*alpha[pp] - deriv_evars->grad2_0_0_alpha[d_pp]) + 2*DENDRO_igt1*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_1 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_1 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_1 + DENDRO_RIJ1*alpha[pp] - deriv_evars->grad2_0_1_alpha[d_pp]) + 2*DENDRO_igt2*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_2 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_2 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_2 + DENDRO_RIJ2*alpha[pp] - deriv_evars->grad2_0_2_alpha[d_pp]) + DENDRO_igt3*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_3 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_3 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_3 + DENDRO_RIJ3*alpha[pp] - deriv_evars->grad2_1_1_alpha[d_pp]) + 2*DENDRO_igt4*(DENDRO_1 + DENDRO_2 + DENDRO_3 + DENDRO_4 - deriv_evars->grad2_1_2_alpha[d_pp]) + DENDRO_igt5*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_5 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_5 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_5 + DENDRO_RIJ5*alpha[pp] - deriv_evars->grad2_2_2_alpha[d_pp])));
// Dendro: reduced ops: 125
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 125 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*At5[pp];
const double DENDRO_1 = 2*At2[pp];
const double DENDRO_2 = 2*At4[pp];
const double DENDRO_3 = (deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_5;
const double DENDRO_4 = (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_5;
const double DENDRO_5 = (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_5;
const double DENDRO_6 = DENDRO_RIJ5*alpha[pp];

// Dendro: printing variables
//--
At_rhs22[pp] = (deriv_evars->grad_0_At5[d_pp])*beta0[pp] - (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_1_At5[d_pp])*beta1[pp] - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_2_At5[d_pp])*beta2[pp] + (deriv_evars->grad_2_beta0[d_pp])*DENDRO_1 + (deriv_evars->grad_2_beta1[d_pp])*DENDRO_2 + (4.0/3.0)*(deriv_evars->grad_2_beta2[d_pp])*At5[pp] - alpha[pp]*(-At5[pp]*K[pp] + 2*At5[pp]*(At2[pp]*DENDRO_igt2 + At4[pp]*DENDRO_igt4 + At5[pp]*DENDRO_igt5) + DENDRO_1*(At2[pp]*DENDRO_igt0 + At4[pp]*DENDRO_igt1 + At5[pp]*DENDRO_igt2) + DENDRO_2*(At2[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt3 + At5[pp]*DENDRO_igt4)) + (1.0/3.0)*chi[pp]*(3*DENDRO_3 + 3*DENDRO_4 + 3*DENDRO_5 + 3*DENDRO_6 - 3*deriv_evars->grad2_2_2_alpha[d_pp] - gt5[pp]*(DENDRO_igt0*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_0 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_0 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_0 + DENDRO_RIJ0*alpha[pp] - deriv_evars->grad2_0_0_alpha[d_pp]) + 2*DENDRO_igt1*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_1 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_1 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_1 + DENDRO_RIJ1*alpha[pp] - deriv_evars->grad2_0_1_alpha[d_pp]) + 2*DENDRO_igt2*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_2 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_2 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_2 + DENDRO_RIJ2*alpha[pp] - deriv_evars->grad2_0_2_alpha[d_pp]) + DENDRO_igt3*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_3 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_3 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_3 + DENDRO_RIJ3*alpha[pp] - deriv_evars->grad2_1_1_alpha[d_pp]) + 2*DENDRO_igt4*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_4 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_4 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_4 + DENDRO_RIJ4*alpha[pp] - deriv_evars->grad2_1_2_alpha[d_pp]) + DENDRO_igt5*(DENDRO_3 + DENDRO_4 + DENDRO_5 + DENDRO_6 - deriv_evars->grad2_2_2_alpha[d_pp])));
// Dendro: reduced ops: 122
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 333 
// Dendro: printing temp variables
const double DENDRO_0 = pow(gt4[pp], 2);
const double DENDRO_1 = pow(gt1[pp], 2);
const double DENDRO_2 = pow(gt2[pp], 2);
const double DENDRO_3 = gt0[pp]*gt3[pp];
const double DENDRO_4 = gt1[pp]*gt4[pp];
const double DENDRO_5 = chi[pp]/(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt5[pp] - 2*DENDRO_4*gt2[pp]);
const double DENDRO_6 = 2*DENDRO_5;
const double DENDRO_7 = (DENDRO_igt1 * DENDRO_igt1);
const double DENDRO_8 = (DENDRO_igt2 * DENDRO_igt2);
const double DENDRO_9 = At1[pp]*DENDRO_igt0;
const double DENDRO_10 = At2[pp]*DENDRO_igt0;
const double DENDRO_11 = DENDRO_igt1*DENDRO_igt2;
const double DENDRO_12 = 2*At4[pp];
const double DENDRO_13 = (DENDRO_igt4 * DENDRO_igt4);
const double DENDRO_14 = At1[pp]*DENDRO_igt1;
const double DENDRO_15 = At2[pp]*DENDRO_igt1;
const double DENDRO_16 = 2*DENDRO_igt4;
const double DENDRO_17 = DENDRO_igt3*DENDRO_igt4;
const double DENDRO_18 = At1[pp]*DENDRO_igt2;
const double DENDRO_19 = At2[pp]*DENDRO_igt2;
const double DENDRO_20 = At4[pp]*DENDRO_igt5;
const double DENDRO_21 = At0[pp]*DENDRO_igt0;
const double DENDRO_22 = At3[pp]*DENDRO_igt1;
const double DENDRO_23 = At4[pp]*DENDRO_igt1;
const double DENDRO_24 = At4[pp]*DENDRO_igt2;
const double DENDRO_25 = At5[pp]*DENDRO_igt2;

// Dendro: printing variables
//--
K_rhs[pp] = (deriv_evars->grad_0_K[d_pp])*beta0[pp] + (deriv_evars->grad_1_K[d_pp])*beta1[pp] + (deriv_evars->grad_2_K[d_pp])*beta2[pp] - DENDRO_5*(-DENDRO_0 + gt3[pp]*gt5[pp])*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_0 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_0 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_0 - deriv_evars->grad2_0_0_alpha[d_pp]) - DENDRO_5*(-DENDRO_1 + DENDRO_3)*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_5 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_5 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_5 - deriv_evars->grad2_2_2_alpha[d_pp]) - DENDRO_5*(-DENDRO_2 + gt0[pp]*gt5[pp])*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_3 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_3 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_3 - deriv_evars->grad2_1_1_alpha[d_pp]) - DENDRO_6*(DENDRO_4 - gt2[pp]*gt3[pp])*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_2 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_2 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_2 - deriv_evars->grad2_0_2_alpha[d_pp]) + DENDRO_6*(gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp])*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_4 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_4 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_4 - deriv_evars->grad2_1_2_alpha[d_pp]) + DENDRO_6*(gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp])*((deriv_evars->grad_0_alpha[d_pp])*DENDRO_C3_k0_1 + (deriv_evars->grad_1_alpha[d_pp])*DENDRO_C3_k1_1 + (deriv_evars->grad_2_alpha[d_pp])*DENDRO_C3_k2_1 - deriv_evars->grad2_0_1_alpha[d_pp]) + (1.0/3.0)*alpha[pp]*(3*At0[pp]*(At0[pp]*(DENDRO_igt0 * DENDRO_igt0) + At3[pp]*DENDRO_7 + At5[pp]*DENDRO_8 + 2*DENDRO_10*DENDRO_igt2 + DENDRO_11*DENDRO_12 + 2*DENDRO_9*DENDRO_igt1) + 6*At1[pp]*(At1[pp]*DENDRO_7 + At2[pp]*DENDRO_11 + DENDRO_10*DENDRO_igt4 + DENDRO_21*DENDRO_igt1 + DENDRO_22*DENDRO_igt3 + DENDRO_23*DENDRO_igt4 + DENDRO_24*DENDRO_igt3 + DENDRO_25*DENDRO_igt4 + DENDRO_9*DENDRO_igt3) + 6*At2[pp]*(At1[pp]*DENDRO_11 + At2[pp]*DENDRO_8 + DENDRO_10*DENDRO_igt5 + DENDRO_21*DENDRO_igt2 + DENDRO_22*DENDRO_igt4 + DENDRO_23*DENDRO_igt5 + DENDRO_24*DENDRO_igt4 + DENDRO_25*DENDRO_igt5 + DENDRO_9*DENDRO_igt4) + 3*At3[pp]*(At0[pp]*DENDRO_7 + At3[pp]*(DENDRO_igt3 * DENDRO_igt3) + At5[pp]*DENDRO_13 + DENDRO_12*DENDRO_17 + 2*DENDRO_14*DENDRO_igt3 + DENDRO_15*DENDRO_16) + 6*At4[pp]*(At0[pp]*DENDRO_11 + At3[pp]*DENDRO_17 + At4[pp]*DENDRO_13 + At5[pp]*DENDRO_igt4*DENDRO_igt5 + DENDRO_14*DENDRO_igt4 + DENDRO_15*DENDRO_igt5 + DENDRO_18*DENDRO_igt3 + DENDRO_19*DENDRO_igt4 + DENDRO_20*DENDRO_igt3) + 3*At5[pp]*(At0[pp]*DENDRO_8 + At3[pp]*DENDRO_13 + At5[pp]*(DENDRO_igt5 * DENDRO_igt5) + DENDRO_16*DENDRO_18 + DENDRO_16*DENDRO_20 + 2*DENDRO_19*DENDRO_igt5) + pow(K[pp], 2));
// Dendro: reduced ops: 222
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 414 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*DENDRO_igt0;
const double DENDRO_1 = (1.0/3.0)*DENDRO_igt1;
const double DENDRO_2 = (1.0/3.0)*DENDRO_igt2;
const double DENDRO_3 = 2*DENDRO_igt4;
const double DENDRO_4 = (DENDRO_igt1 * DENDRO_igt1);
const double DENDRO_5 = (DENDRO_igt2 * DENDRO_igt2);
const double DENDRO_6 = At1[pp]*DENDRO_igt0;
const double DENDRO_7 = At2[pp]*DENDRO_igt0;
const double DENDRO_8 = DENDRO_igt1*DENDRO_igt2;
const double DENDRO_9 = At0[pp]*(DENDRO_igt0 * DENDRO_igt0) + At3[pp]*DENDRO_4 + 2*At4[pp]*DENDRO_8 + At5[pp]*DENDRO_5 + 2*DENDRO_6*DENDRO_igt1 + 2*DENDRO_7*DENDRO_igt2;
const double DENDRO_10 = 2*DENDRO_9;
const double DENDRO_11 = (DENDRO_igt4 * DENDRO_igt4);
const double DENDRO_12 = At1[pp]*DENDRO_igt1;
const double DENDRO_13 = At2[pp]*DENDRO_igt1;
const double DENDRO_14 = At4[pp]*DENDRO_igt3;
const double DENDRO_15 = 2*alpha[pp];
const double DENDRO_16 = At1[pp]*DENDRO_igt2;
const double DENDRO_17 = At2[pp]*DENDRO_igt2;
const double DENDRO_18 = At0[pp]*DENDRO_igt0;
const double DENDRO_19 = At3[pp]*DENDRO_igt1;
const double DENDRO_20 = At4[pp]*DENDRO_igt1;
const double DENDRO_21 = At4[pp]*DENDRO_igt2;
const double DENDRO_22 = At5[pp]*DENDRO_igt2;
const double DENDRO_23 = At1[pp]*DENDRO_4 + At2[pp]*DENDRO_8 + DENDRO_18*DENDRO_igt1 + DENDRO_19*DENDRO_igt3 + DENDRO_20*DENDRO_igt4 + DENDRO_21*DENDRO_igt3 + DENDRO_22*DENDRO_igt4 + DENDRO_6*DENDRO_igt3 + DENDRO_7*DENDRO_igt4;
const double DENDRO_24 = At1[pp]*DENDRO_8 + At2[pp]*DENDRO_5 + DENDRO_18*DENDRO_igt2 + DENDRO_19*DENDRO_igt4 + DENDRO_20*DENDRO_igt5 + DENDRO_21*DENDRO_igt4 + DENDRO_22*DENDRO_igt5 + DENDRO_6*DENDRO_igt4 + DENDRO_7*DENDRO_igt5;
const double DENDRO_25 = 4*alpha[pp];
const double DENDRO_26 = 9/chi[pp];
const double DENDRO_27 = (1.0/3.0)*alpha[pp];

// Dendro: printing variables
//--
Gt_rhs0[pp] = (deriv_evars->grad_0_Gt0[d_pp])*beta0[pp] - (deriv_evars->grad_0_alpha[d_pp])*DENDRO_10 - (deriv_evars->grad_0_beta0[d_pp])*DENDRO_Gtk0 + (deriv_evars->grad_1_Gt0[d_pp])*beta1[pp] - 2*(deriv_evars->grad_1_alpha[d_pp])*DENDRO_23 - (deriv_evars->grad_1_beta0[d_pp])*DENDRO_Gtk1 + (deriv_evars->grad_2_Gt0[d_pp])*beta2[pp] - 2*(deriv_evars->grad_2_alpha[d_pp])*DENDRO_24 - (deriv_evars->grad_2_beta0[d_pp])*DENDRO_Gtk2 + DENDRO_0*deriv_evars->grad2_0_1_beta1[d_pp] + DENDRO_0*deriv_evars->grad2_0_2_beta2[d_pp] + DENDRO_1*deriv_evars->grad2_1_1_beta1[d_pp] + DENDRO_1*deriv_evars->grad2_1_2_beta2[d_pp] + DENDRO_10*DENDRO_C2_k0_0*alpha[pp] + DENDRO_15*DENDRO_C2_k0_3*(At0[pp]*DENDRO_4 + At3[pp]*(DENDRO_igt3 * DENDRO_igt3) + At5[pp]*DENDRO_11 + 2*DENDRO_12*DENDRO_igt3 + DENDRO_13*DENDRO_3 + DENDRO_14*DENDRO_3) + DENDRO_15*DENDRO_C2_k0_5*(At0[pp]*DENDRO_5 + At3[pp]*DENDRO_11 + At4[pp]*DENDRO_3*DENDRO_igt5 + At5[pp]*(DENDRO_igt5 * DENDRO_igt5) + DENDRO_16*DENDRO_3 + 2*DENDRO_17*DENDRO_igt5) + DENDRO_2*deriv_evars->grad2_1_2_beta1[d_pp] + DENDRO_2*deriv_evars->grad2_2_2_beta2[d_pp] + DENDRO_23*DENDRO_25*DENDRO_C2_k0_1 + DENDRO_24*DENDRO_25*DENDRO_C2_k0_2 + DENDRO_25*DENDRO_C2_k0_4*(At0[pp]*DENDRO_8 + At3[pp]*DENDRO_igt3*DENDRO_igt4 + At4[pp]*DENDRO_11 + At5[pp]*DENDRO_igt4*DENDRO_igt5 + DENDRO_12*DENDRO_igt4 + DENDRO_13*DENDRO_igt5 + DENDRO_14*DENDRO_igt5 + DENDRO_16*DENDRO_igt3 + DENDRO_17*DENDRO_igt4) - DENDRO_27*(4*(deriv_evars->grad_0_K[d_pp])*DENDRO_igt0 + (deriv_evars->grad_0_chi[d_pp])*DENDRO_26*DENDRO_9) - DENDRO_27*(4*(deriv_evars->grad_1_K[d_pp])*DENDRO_igt1 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_23*DENDRO_26) - DENDRO_27*(4*(deriv_evars->grad_2_K[d_pp])*DENDRO_igt2 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_24*DENDRO_26) + DENDRO_3*deriv_evars->grad2_1_2_beta0[d_pp] + (2.0/3.0)*DENDRO_Gtk0*((deriv_evars->grad_0_beta0[d_pp]) + (deriv_evars->grad_1_beta1[d_pp]) + (deriv_evars->grad_2_beta2[d_pp])) + (4.0/3.0)*DENDRO_igt0*deriv_evars->grad2_0_0_beta0[d_pp] + (7.0/3.0)*DENDRO_igt1*deriv_evars->grad2_0_1_beta0[d_pp] + (7.0/3.0)*DENDRO_igt2*deriv_evars->grad2_0_2_beta0[d_pp] + DENDRO_igt3*deriv_evars->grad2_1_1_beta0[d_pp] + DENDRO_igt5*deriv_evars->grad2_2_2_beta0[d_pp];
// Dendro: reduced ops: 214
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 414 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*DENDRO_igt1;
const double DENDRO_1 = 2*DENDRO_igt2;
const double DENDRO_2 = (1.0/3.0)*DENDRO_igt3;
const double DENDRO_3 = (1.0/3.0)*DENDRO_igt4;
const double DENDRO_4 = (DENDRO_igt1 * DENDRO_igt1);
const double DENDRO_5 = (DENDRO_igt4 * DENDRO_igt4);
const double DENDRO_6 = At1[pp]*DENDRO_igt3;
const double DENDRO_7 = 2*DENDRO_igt1;
const double DENDRO_8 = At2[pp]*DENDRO_igt4;
const double DENDRO_9 = At4[pp]*DENDRO_igt3;
const double DENDRO_10 = At0[pp]*DENDRO_4 + At3[pp]*(DENDRO_igt3 * DENDRO_igt3) + At5[pp]*DENDRO_5 + DENDRO_6*DENDRO_7 + DENDRO_7*DENDRO_8 + 2*DENDRO_9*DENDRO_igt4;
const double DENDRO_11 = 2*DENDRO_10;
const double DENDRO_12 = (DENDRO_igt2 * DENDRO_igt2);
const double DENDRO_13 = At1[pp]*DENDRO_igt0;
const double DENDRO_14 = At2[pp]*DENDRO_igt0;
const double DENDRO_15 = At4[pp]*DENDRO_igt1;
const double DENDRO_16 = 2*alpha[pp];
const double DENDRO_17 = At4[pp]*DENDRO_igt4;
const double DENDRO_18 = At0[pp]*DENDRO_igt1;
const double DENDRO_19 = At2[pp]*DENDRO_igt1;
const double DENDRO_20 = At3[pp]*DENDRO_igt3;
const double DENDRO_21 = DENDRO_igt1*DENDRO_igt4;
const double DENDRO_22 = At5[pp]*DENDRO_igt4;
const double DENDRO_23 = At1[pp]*DENDRO_4 + At4[pp]*DENDRO_21 + DENDRO_18*DENDRO_igt0 + DENDRO_19*DENDRO_igt2 + DENDRO_20*DENDRO_igt1 + DENDRO_22*DENDRO_igt2 + DENDRO_6*DENDRO_igt0 + DENDRO_8*DENDRO_igt0 + DENDRO_9*DENDRO_igt2;
const double DENDRO_24 = At1[pp]*DENDRO_21 + At4[pp]*DENDRO_5 + DENDRO_18*DENDRO_igt2 + DENDRO_19*DENDRO_igt5 + DENDRO_20*DENDRO_igt4 + DENDRO_22*DENDRO_igt5 + DENDRO_6*DENDRO_igt2 + DENDRO_8*DENDRO_igt2 + DENDRO_9*DENDRO_igt5;
const double DENDRO_25 = 4*alpha[pp];
const double DENDRO_26 = 9/chi[pp];
const double DENDRO_27 = (1.0/3.0)*alpha[pp];

// Dendro: printing variables
//--
Gt_rhs1[pp] = (deriv_evars->grad_0_Gt1[d_pp])*beta0[pp] - 2*(deriv_evars->grad_0_alpha[d_pp])*DENDRO_23 - (deriv_evars->grad_0_beta1[d_pp])*DENDRO_Gtk0 + (deriv_evars->grad_1_Gt1[d_pp])*beta1[pp] - (deriv_evars->grad_1_alpha[d_pp])*DENDRO_11 - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_Gtk1 + (deriv_evars->grad_2_Gt1[d_pp])*beta2[pp] - 2*(deriv_evars->grad_2_alpha[d_pp])*DENDRO_24 - (deriv_evars->grad_2_beta1[d_pp])*DENDRO_Gtk2 + DENDRO_0*deriv_evars->grad2_0_0_beta0[d_pp] + DENDRO_0*deriv_evars->grad2_0_2_beta2[d_pp] + DENDRO_1*deriv_evars->grad2_0_2_beta1[d_pp] + DENDRO_11*DENDRO_C2_k1_3*alpha[pp] + DENDRO_16*DENDRO_C2_k1_0*(At0[pp]*(DENDRO_igt0 * DENDRO_igt0) + At3[pp]*DENDRO_4 + At5[pp]*DENDRO_12 + DENDRO_1*DENDRO_14 + DENDRO_1*DENDRO_15 + DENDRO_13*DENDRO_7) + DENDRO_16*DENDRO_C2_k1_5*(At0[pp]*DENDRO_12 + At1[pp]*DENDRO_1*DENDRO_igt4 + At2[pp]*DENDRO_1*DENDRO_igt5 + At3[pp]*DENDRO_5 + At5[pp]*(DENDRO_igt5 * DENDRO_igt5) + 2*DENDRO_17*DENDRO_igt5) + DENDRO_2*deriv_evars->grad2_0_1_beta0[d_pp] + DENDRO_2*deriv_evars->grad2_1_2_beta2[d_pp] + DENDRO_23*DENDRO_25*DENDRO_C2_k1_1 + DENDRO_24*DENDRO_25*DENDRO_C2_k1_4 + DENDRO_25*DENDRO_C2_k1_2*(At0[pp]*DENDRO_igt0*DENDRO_igt2 + At1[pp]*DENDRO_igt1*DENDRO_igt2 + At2[pp]*DENDRO_12 + At3[pp]*DENDRO_21 + At5[pp]*DENDRO_igt2*DENDRO_igt5 + DENDRO_13*DENDRO_igt4 + DENDRO_14*DENDRO_igt5 + DENDRO_15*DENDRO_igt5 + DENDRO_17*DENDRO_igt2) - DENDRO_27*(4*(deriv_evars->grad_0_K[d_pp])*DENDRO_igt1 + (deriv_evars->grad_0_chi[d_pp])*DENDRO_23*DENDRO_26) - DENDRO_27*(4*(deriv_evars->grad_1_K[d_pp])*DENDRO_igt3 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_10*DENDRO_26) - DENDRO_27*(4*(deriv_evars->grad_2_K[d_pp])*DENDRO_igt4 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_24*DENDRO_26) + DENDRO_3*deriv_evars->grad2_0_2_beta0[d_pp] + DENDRO_3*deriv_evars->grad2_2_2_beta2[d_pp] + (2.0/3.0)*DENDRO_Gtk1*((deriv_evars->grad_0_beta0[d_pp]) + (deriv_evars->grad_1_beta1[d_pp]) + (deriv_evars->grad_2_beta2[d_pp])) + DENDRO_igt0*deriv_evars->grad2_0_0_beta1[d_pp] + (7.0/3.0)*DENDRO_igt1*deriv_evars->grad2_0_1_beta1[d_pp] + (4.0/3.0)*DENDRO_igt3*deriv_evars->grad2_1_1_beta1[d_pp] + (7.0/3.0)*DENDRO_igt4*deriv_evars->grad2_1_2_beta1[d_pp] + DENDRO_igt5*deriv_evars->grad2_2_2_beta1[d_pp];
// Dendro: reduced ops: 213
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 414 
// Dendro: printing temp variables
const double DENDRO_0 = 2*DENDRO_igt1;
const double DENDRO_1 = (1.0/3.0)*DENDRO_igt2;
const double DENDRO_2 = (1.0/3.0)*DENDRO_igt4;
const double DENDRO_3 = (1.0/3.0)*DENDRO_igt5;
const double DENDRO_4 = (DENDRO_igt2 * DENDRO_igt2);
const double DENDRO_5 = (DENDRO_igt4 * DENDRO_igt4);
const double DENDRO_6 = At1[pp]*DENDRO_igt4;
const double DENDRO_7 = 2*DENDRO_igt2;
const double DENDRO_8 = At2[pp]*DENDRO_igt5;
const double DENDRO_9 = At4[pp]*DENDRO_igt5;
const double DENDRO_10 = At0[pp]*DENDRO_4 + At3[pp]*DENDRO_5 + At5[pp]*(DENDRO_igt5 * DENDRO_igt5) + DENDRO_6*DENDRO_7 + DENDRO_7*DENDRO_8 + 2*DENDRO_9*DENDRO_igt4;
const double DENDRO_11 = 2*DENDRO_10;
const double DENDRO_12 = (DENDRO_igt1 * DENDRO_igt1);
const double DENDRO_13 = At1[pp]*DENDRO_igt0;
const double DENDRO_14 = At2[pp]*DENDRO_igt0;
const double DENDRO_15 = At4[pp]*DENDRO_igt2;
const double DENDRO_16 = 2*alpha[pp];
const double DENDRO_17 = At4[pp]*DENDRO_igt4;
const double DENDRO_18 = At0[pp]*DENDRO_igt2;
const double DENDRO_19 = At1[pp]*DENDRO_igt2;
const double DENDRO_20 = At3[pp]*DENDRO_igt4;
const double DENDRO_21 = DENDRO_igt2*DENDRO_igt4;
const double DENDRO_22 = At5[pp]*DENDRO_igt5;
const double DENDRO_23 = At2[pp]*DENDRO_4 + At4[pp]*DENDRO_21 + DENDRO_18*DENDRO_igt0 + DENDRO_19*DENDRO_igt1 + DENDRO_20*DENDRO_igt1 + DENDRO_22*DENDRO_igt2 + DENDRO_6*DENDRO_igt0 + DENDRO_8*DENDRO_igt0 + DENDRO_9*DENDRO_igt1;
const double DENDRO_24 = At2[pp]*DENDRO_21 + At4[pp]*DENDRO_5 + DENDRO_18*DENDRO_igt1 + DENDRO_19*DENDRO_igt3 + DENDRO_20*DENDRO_igt3 + DENDRO_22*DENDRO_igt4 + DENDRO_6*DENDRO_igt1 + DENDRO_8*DENDRO_igt1 + DENDRO_9*DENDRO_igt3;
const double DENDRO_25 = 4*alpha[pp];
const double DENDRO_26 = 9/chi[pp];
const double DENDRO_27 = (1.0/3.0)*alpha[pp];

// Dendro: printing variables
//--
Gt_rhs2[pp] = (deriv_evars->grad_0_Gt2[d_pp])*beta0[pp] - 2*(deriv_evars->grad_0_alpha[d_pp])*DENDRO_23 - (deriv_evars->grad_0_beta2[d_pp])*DENDRO_Gtk0 + (deriv_evars->grad_1_Gt2[d_pp])*beta1[pp] - 2*(deriv_evars->grad_1_alpha[d_pp])*DENDRO_24 - (deriv_evars->grad_1_beta2[d_pp])*DENDRO_Gtk1 + (deriv_evars->grad_2_Gt2[d_pp])*beta2[pp] - (deriv_evars->grad_2_alpha[d_pp])*DENDRO_11 - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_Gtk2 + DENDRO_0*deriv_evars->grad2_0_1_beta2[d_pp] + DENDRO_1*deriv_evars->grad2_0_0_beta0[d_pp] + DENDRO_1*deriv_evars->grad2_0_1_beta1[d_pp] + DENDRO_11*DENDRO_C2_k2_5*alpha[pp] + DENDRO_16*DENDRO_C2_k2_0*(At0[pp]*(DENDRO_igt0 * DENDRO_igt0) + At3[pp]*DENDRO_12 + At5[pp]*DENDRO_4 + DENDRO_0*DENDRO_13 + DENDRO_0*DENDRO_15 + DENDRO_14*DENDRO_7) + DENDRO_16*DENDRO_C2_k2_3*(At0[pp]*DENDRO_12 + At1[pp]*DENDRO_0*DENDRO_igt3 + At2[pp]*DENDRO_0*DENDRO_igt4 + At3[pp]*(DENDRO_igt3 * DENDRO_igt3) + At5[pp]*DENDRO_5 + 2*DENDRO_17*DENDRO_igt3) + DENDRO_2*deriv_evars->grad2_0_1_beta0[d_pp] + DENDRO_2*deriv_evars->grad2_1_1_beta1[d_pp] + DENDRO_23*DENDRO_25*DENDRO_C2_k2_2 + DENDRO_24*DENDRO_25*DENDRO_C2_k2_4 + DENDRO_25*DENDRO_C2_k2_1*(At0[pp]*DENDRO_igt0*DENDRO_igt1 + At1[pp]*DENDRO_12 + At2[pp]*DENDRO_igt1*DENDRO_igt2 + At3[pp]*DENDRO_igt1*DENDRO_igt3 + At5[pp]*DENDRO_21 + DENDRO_13*DENDRO_igt3 + DENDRO_14*DENDRO_igt4 + DENDRO_15*DENDRO_igt3 + DENDRO_17*DENDRO_igt1) - DENDRO_27*(4*(deriv_evars->grad_0_K[d_pp])*DENDRO_igt2 + (deriv_evars->grad_0_chi[d_pp])*DENDRO_23*DENDRO_26) - DENDRO_27*(4*(deriv_evars->grad_1_K[d_pp])*DENDRO_igt4 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_24*DENDRO_26) - DENDRO_27*(4*(deriv_evars->grad_2_K[d_pp])*DENDRO_igt5 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_10*DENDRO_26) + DENDRO_3*deriv_evars->grad2_0_2_beta0[d_pp] + DENDRO_3*deriv_evars->grad2_1_2_beta1[d_pp] + (2.0/3.0)*DENDRO_Gtk2*((deriv_evars->grad_0_beta0[d_pp]) + (deriv_evars->grad_1_beta1[d_pp]) + (deriv_evars->grad_2_beta2[d_pp])) + DENDRO_igt0*deriv_evars->grad2_0_0_beta2[d_pp] + (7.0/3.0)*DENDRO_igt2*deriv_evars->grad2_0_2_beta2[d_pp] + DENDRO_igt3*deriv_evars->grad2_1_1_beta2[d_pp] + (7.0/3.0)*DENDRO_igt4*deriv_evars->grad2_1_2_beta2[d_pp] + (4.0/3.0)*DENDRO_igt5*deriv_evars->grad2_2_2_beta2[d_pp];
// Dendro: reduced ops: 213
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 430 
// Dendro: printing temp variables
const double DENDRO_0 = (deriv_evars->grad_0_Gt0[d_pp])*beta0[pp] + (deriv_evars->grad_1_Gt0[d_pp])*beta1[pp] + (deriv_evars->grad_2_Gt0[d_pp])*beta2[pp];
const double DENDRO_1 = (DENDRO_igt1 * DENDRO_igt1);
const double DENDRO_2 = (DENDRO_igt2 * DENDRO_igt2);
const double DENDRO_3 = At1[pp]*DENDRO_igt0;
const double DENDRO_4 = At2[pp]*DENDRO_igt0;
const double DENDRO_5 = DENDRO_igt1*DENDRO_igt2;
const double DENDRO_6 = At0[pp]*(DENDRO_igt0 * DENDRO_igt0) + At3[pp]*DENDRO_1 + 2*At4[pp]*DENDRO_5 + At5[pp]*DENDRO_2 + 2*DENDRO_3*DENDRO_igt1 + 2*DENDRO_4*DENDRO_igt2;
const double DENDRO_7 = 2*DENDRO_6;
const double DENDRO_8 = At0[pp]*DENDRO_igt0;
const double DENDRO_9 = At3[pp]*DENDRO_igt1;
const double DENDRO_10 = At4[pp]*DENDRO_igt1;
const double DENDRO_11 = At4[pp]*DENDRO_igt2;
const double DENDRO_12 = At5[pp]*DENDRO_igt2;
const double DENDRO_13 = At1[pp]*DENDRO_1 + At2[pp]*DENDRO_5 + DENDRO_10*DENDRO_igt4 + DENDRO_11*DENDRO_igt3 + DENDRO_12*DENDRO_igt4 + DENDRO_3*DENDRO_igt3 + DENDRO_4*DENDRO_igt4 + DENDRO_8*DENDRO_igt1 + DENDRO_9*DENDRO_igt3;
const double DENDRO_14 = At1[pp]*DENDRO_5 + At2[pp]*DENDRO_2 + DENDRO_10*DENDRO_igt5 + DENDRO_11*DENDRO_igt4 + DENDRO_12*DENDRO_igt5 + DENDRO_3*DENDRO_igt4 + DENDRO_4*DENDRO_igt5 + DENDRO_8*DENDRO_igt2 + DENDRO_9*DENDRO_igt4;
const double DENDRO_15 = 2*DENDRO_igt4;
const double DENDRO_16 = 9/chi[pp];
const double DENDRO_17 = (1.0/3.0)*alpha[pp];
const double DENDRO_18 = (1.0/3.0)*DENDRO_igt0;
const double DENDRO_19 = (1.0/3.0)*DENDRO_igt1;
const double DENDRO_20 = (1.0/3.0)*DENDRO_igt2;
const double DENDRO_21 = (DENDRO_igt4 * DENDRO_igt4);
const double DENDRO_22 = At1[pp]*DENDRO_igt1;
const double DENDRO_23 = At2[pp]*DENDRO_igt1;
const double DENDRO_24 = At4[pp]*DENDRO_igt3;
const double DENDRO_25 = 2*alpha[pp];
const double DENDRO_26 = At1[pp]*DENDRO_igt2;
const double DENDRO_27 = At2[pp]*DENDRO_igt2;
const double DENDRO_28 = 4*alpha[pp];

// Dendro: printing variables
//--
B_rhs0[pp] = -(deriv_evars->grad_0_alpha[d_pp])*DENDRO_7 - (deriv_evars->grad_0_beta0[d_pp])*DENDRO_Gtk0 - 2*(deriv_evars->grad_1_alpha[d_pp])*DENDRO_13 - (deriv_evars->grad_1_beta0[d_pp])*DENDRO_Gtk1 - 2*(deriv_evars->grad_2_alpha[d_pp])*DENDRO_14 - (deriv_evars->grad_2_beta0[d_pp])*DENDRO_Gtk2 - B0[pp]*eta - DENDRO_0*lambda[3] + DENDRO_0 + DENDRO_13*DENDRO_28*DENDRO_C2_k0_1 + DENDRO_14*DENDRO_28*DENDRO_C2_k0_2 + DENDRO_15*deriv_evars->grad2_1_2_beta0[d_pp] - DENDRO_17*(4*(deriv_evars->grad_0_K[d_pp])*DENDRO_igt0 + (deriv_evars->grad_0_chi[d_pp])*DENDRO_16*DENDRO_6) - DENDRO_17*(4*(deriv_evars->grad_1_K[d_pp])*DENDRO_igt1 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_13*DENDRO_16) - DENDRO_17*(4*(deriv_evars->grad_2_K[d_pp])*DENDRO_igt2 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_14*DENDRO_16) + DENDRO_18*deriv_evars->grad2_0_1_beta1[d_pp] + DENDRO_18*deriv_evars->grad2_0_2_beta2[d_pp] + DENDRO_19*deriv_evars->grad2_1_1_beta1[d_pp] + DENDRO_19*deriv_evars->grad2_1_2_beta2[d_pp] + DENDRO_20*deriv_evars->grad2_1_2_beta1[d_pp] + DENDRO_20*deriv_evars->grad2_2_2_beta2[d_pp] + DENDRO_25*DENDRO_C2_k0_3*(At0[pp]*DENDRO_1 + At3[pp]*(DENDRO_igt3 * DENDRO_igt3) + At5[pp]*DENDRO_21 + DENDRO_15*DENDRO_23 + DENDRO_15*DENDRO_24 + 2*DENDRO_22*DENDRO_igt3) + DENDRO_25*DENDRO_C2_k0_5*(At0[pp]*DENDRO_2 + At3[pp]*DENDRO_21 + At4[pp]*DENDRO_15*DENDRO_igt5 + At5[pp]*(DENDRO_igt5 * DENDRO_igt5) + DENDRO_15*DENDRO_26 + 2*DENDRO_27*DENDRO_igt5) + DENDRO_28*DENDRO_C2_k0_4*(At0[pp]*DENDRO_5 + At3[pp]*DENDRO_igt3*DENDRO_igt4 + At4[pp]*DENDRO_21 + At5[pp]*DENDRO_igt4*DENDRO_igt5 + DENDRO_22*DENDRO_igt4 + DENDRO_23*DENDRO_igt5 + DENDRO_24*DENDRO_igt5 + DENDRO_26*DENDRO_igt3 + DENDRO_27*DENDRO_igt4) + DENDRO_7*DENDRO_C2_k0_0*alpha[pp] + (2.0/3.0)*DENDRO_Gtk0*((deriv_evars->grad_0_beta0[d_pp]) + (deriv_evars->grad_1_beta1[d_pp]) + (deriv_evars->grad_2_beta2[d_pp])) + (4.0/3.0)*DENDRO_igt0*deriv_evars->grad2_0_0_beta0[d_pp] + (7.0/3.0)*DENDRO_igt1*deriv_evars->grad2_0_1_beta0[d_pp] + (7.0/3.0)*DENDRO_igt2*deriv_evars->grad2_0_2_beta0[d_pp] + DENDRO_igt3*deriv_evars->grad2_1_1_beta0[d_pp] + DENDRO_igt5*deriv_evars->grad2_2_2_beta0[d_pp] + lambda[2]*((deriv_evars->grad_0_B0[d_pp])*beta0[pp] + (deriv_evars->grad_1_B0[d_pp])*beta1[pp] + (deriv_evars->grad_2_B0[d_pp])*beta2[pp]);
// Dendro: reduced ops: 225
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 430 
// Dendro: printing temp variables
const double DENDRO_0 = (deriv_evars->grad_0_Gt1[d_pp])*beta0[pp] + (deriv_evars->grad_1_Gt1[d_pp])*beta1[pp] + (deriv_evars->grad_2_Gt1[d_pp])*beta2[pp];
const double DENDRO_1 = At0[pp]*DENDRO_igt1;
const double DENDRO_2 = At1[pp]*DENDRO_igt3;
const double DENDRO_3 = At2[pp]*DENDRO_igt4;
const double DENDRO_4 = At2[pp]*DENDRO_igt1;
const double DENDRO_5 = At3[pp]*DENDRO_igt3;
const double DENDRO_6 = DENDRO_igt1*DENDRO_igt4;
const double DENDRO_7 = At4[pp]*DENDRO_igt3;
const double DENDRO_8 = At5[pp]*DENDRO_igt4;
const double DENDRO_9 = (DENDRO_igt1 * DENDRO_igt1);
const double DENDRO_10 = At1[pp]*DENDRO_9 + At4[pp]*DENDRO_6 + DENDRO_1*DENDRO_igt0 + DENDRO_2*DENDRO_igt0 + DENDRO_3*DENDRO_igt0 + DENDRO_4*DENDRO_igt2 + DENDRO_5*DENDRO_igt1 + DENDRO_7*DENDRO_igt2 + DENDRO_8*DENDRO_igt2;
const double DENDRO_11 = (DENDRO_igt4 * DENDRO_igt4);
const double DENDRO_12 = 2*DENDRO_igt1;
const double DENDRO_13 = At0[pp]*DENDRO_9 + At3[pp]*(DENDRO_igt3 * DENDRO_igt3) + At5[pp]*DENDRO_11 + DENDRO_12*DENDRO_2 + DENDRO_12*DENDRO_3 + 2*DENDRO_7*DENDRO_igt4;
const double DENDRO_14 = 2*DENDRO_13;
const double DENDRO_15 = At1[pp]*DENDRO_6 + At4[pp]*DENDRO_11 + DENDRO_1*DENDRO_igt2 + DENDRO_2*DENDRO_igt2 + DENDRO_3*DENDRO_igt2 + DENDRO_4*DENDRO_igt5 + DENDRO_5*DENDRO_igt4 + DENDRO_7*DENDRO_igt5 + DENDRO_8*DENDRO_igt5;
const double DENDRO_16 = 2*DENDRO_igt2;
const double DENDRO_17 = 9/chi[pp];
const double DENDRO_18 = (1.0/3.0)*alpha[pp];
const double DENDRO_19 = (1.0/3.0)*DENDRO_igt1;
const double DENDRO_20 = (1.0/3.0)*DENDRO_igt3;
const double DENDRO_21 = (1.0/3.0)*DENDRO_igt4;
const double DENDRO_22 = (DENDRO_igt2 * DENDRO_igt2);
const double DENDRO_23 = At1[pp]*DENDRO_igt0;
const double DENDRO_24 = At2[pp]*DENDRO_igt0;
const double DENDRO_25 = At4[pp]*DENDRO_igt1;
const double DENDRO_26 = 2*alpha[pp];
const double DENDRO_27 = At4[pp]*DENDRO_igt4;
const double DENDRO_28 = 4*alpha[pp];

// Dendro: printing variables
//--
B_rhs1[pp] = -2*(deriv_evars->grad_0_alpha[d_pp])*DENDRO_10 - (deriv_evars->grad_0_beta1[d_pp])*DENDRO_Gtk0 - (deriv_evars->grad_1_alpha[d_pp])*DENDRO_14 - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_Gtk1 - 2*(deriv_evars->grad_2_alpha[d_pp])*DENDRO_15 - (deriv_evars->grad_2_beta1[d_pp])*DENDRO_Gtk2 - B1[pp]*eta - DENDRO_0*lambda[3] + DENDRO_0 + DENDRO_10*DENDRO_28*DENDRO_C2_k1_1 + DENDRO_14*DENDRO_C2_k1_3*alpha[pp] + DENDRO_15*DENDRO_28*DENDRO_C2_k1_4 + DENDRO_16*deriv_evars->grad2_0_2_beta1[d_pp] - DENDRO_18*(4*(deriv_evars->grad_0_K[d_pp])*DENDRO_igt1 + (deriv_evars->grad_0_chi[d_pp])*DENDRO_10*DENDRO_17) - DENDRO_18*(4*(deriv_evars->grad_1_K[d_pp])*DENDRO_igt3 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_13*DENDRO_17) - DENDRO_18*(4*(deriv_evars->grad_2_K[d_pp])*DENDRO_igt4 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_15*DENDRO_17) + DENDRO_19*deriv_evars->grad2_0_0_beta0[d_pp] + DENDRO_19*deriv_evars->grad2_0_2_beta2[d_pp] + DENDRO_20*deriv_evars->grad2_0_1_beta0[d_pp] + DENDRO_20*deriv_evars->grad2_1_2_beta2[d_pp] + DENDRO_21*deriv_evars->grad2_0_2_beta0[d_pp] + DENDRO_21*deriv_evars->grad2_2_2_beta2[d_pp] + DENDRO_26*DENDRO_C2_k1_0*(At0[pp]*(DENDRO_igt0 * DENDRO_igt0) + At3[pp]*DENDRO_9 + At5[pp]*DENDRO_22 + DENDRO_12*DENDRO_23 + DENDRO_16*DENDRO_24 + DENDRO_16*DENDRO_25) + DENDRO_26*DENDRO_C2_k1_5*(At0[pp]*DENDRO_22 + At1[pp]*DENDRO_16*DENDRO_igt4 + At2[pp]*DENDRO_16*DENDRO_igt5 + At3[pp]*DENDRO_11 + At5[pp]*(DENDRO_igt5 * DENDRO_igt5) + 2*DENDRO_27*DENDRO_igt5) + DENDRO_28*DENDRO_C2_k1_2*(At0[pp]*DENDRO_igt0*DENDRO_igt2 + At1[pp]*DENDRO_igt1*DENDRO_igt2 + At2[pp]*DENDRO_22 + At3[pp]*DENDRO_6 + At5[pp]*DENDRO_igt2*DENDRO_igt5 + DENDRO_23*DENDRO_igt4 + DENDRO_24*DENDRO_igt5 + DENDRO_25*DENDRO_igt5 + DENDRO_27*DENDRO_igt2) + (2.0/3.0)*DENDRO_Gtk1*((deriv_evars->grad_0_beta0[d_pp]) + (deriv_evars->grad_1_beta1[d_pp]) + (deriv_evars->grad_2_beta2[d_pp])) + DENDRO_igt0*deriv_evars->grad2_0_0_beta1[d_pp] + (7.0/3.0)*DENDRO_igt1*deriv_evars->grad2_0_1_beta1[d_pp] + (4.0/3.0)*DENDRO_igt3*deriv_evars->grad2_1_1_beta1[d_pp] + (7.0/3.0)*DENDRO_igt4*deriv_evars->grad2_1_2_beta1[d_pp] + DENDRO_igt5*deriv_evars->grad2_2_2_beta1[d_pp] + lambda[2]*((deriv_evars->grad_0_B1[d_pp])*beta0[pp] + (deriv_evars->grad_1_B1[d_pp])*beta1[pp] + (deriv_evars->grad_2_B1[d_pp])*beta2[pp]);
// Dendro: reduced ops: 224
// Dendro: }}} 
}
{
// Dendro: {{{ 
// Dendro: original ops: 430 
// Dendro: printing temp variables
const double DENDRO_0 = (deriv_evars->grad_0_Gt2[d_pp])*beta0[pp] + (deriv_evars->grad_1_Gt2[d_pp])*beta1[pp] + (deriv_evars->grad_2_Gt2[d_pp])*beta2[pp];
const double DENDRO_1 = At0[pp]*DENDRO_igt2;
const double DENDRO_2 = At1[pp]*DENDRO_igt4;
const double DENDRO_3 = At1[pp]*DENDRO_igt2;
const double DENDRO_4 = At2[pp]*DENDRO_igt5;
const double DENDRO_5 = At3[pp]*DENDRO_igt4;
const double DENDRO_6 = At4[pp]*DENDRO_igt5;
const double DENDRO_7 = DENDRO_igt2*DENDRO_igt4;
const double DENDRO_8 = At5[pp]*DENDRO_igt5;
const double DENDRO_9 = (DENDRO_igt2 * DENDRO_igt2);
const double DENDRO_10 = At2[pp]*DENDRO_9 + At4[pp]*DENDRO_7 + DENDRO_1*DENDRO_igt0 + DENDRO_2*DENDRO_igt0 + DENDRO_3*DENDRO_igt1 + DENDRO_4*DENDRO_igt0 + DENDRO_5*DENDRO_igt1 + DENDRO_6*DENDRO_igt1 + DENDRO_8*DENDRO_igt2;
const double DENDRO_11 = (DENDRO_igt4 * DENDRO_igt4);
const double DENDRO_12 = At2[pp]*DENDRO_7 + At4[pp]*DENDRO_11 + DENDRO_1*DENDRO_igt1 + DENDRO_2*DENDRO_igt1 + DENDRO_3*DENDRO_igt3 + DENDRO_4*DENDRO_igt1 + DENDRO_5*DENDRO_igt3 + DENDRO_6*DENDRO_igt3 + DENDRO_8*DENDRO_igt4;
const double DENDRO_13 = 2*DENDRO_igt2;
const double DENDRO_14 = At0[pp]*DENDRO_9 + At3[pp]*DENDRO_11 + At5[pp]*(DENDRO_igt5 * DENDRO_igt5) + DENDRO_13*DENDRO_2 + DENDRO_13*DENDRO_4 + 2*DENDRO_6*DENDRO_igt4;
const double DENDRO_15 = 2*DENDRO_14;
const double DENDRO_16 = 2*DENDRO_igt1;
const double DENDRO_17 = 9/chi[pp];
const double DENDRO_18 = (1.0/3.0)*alpha[pp];
const double DENDRO_19 = (1.0/3.0)*DENDRO_igt2;
const double DENDRO_20 = (1.0/3.0)*DENDRO_igt4;
const double DENDRO_21 = (1.0/3.0)*DENDRO_igt5;
const double DENDRO_22 = (DENDRO_igt1 * DENDRO_igt1);
const double DENDRO_23 = At1[pp]*DENDRO_igt0;
const double DENDRO_24 = At2[pp]*DENDRO_igt0;
const double DENDRO_25 = At4[pp]*DENDRO_igt2;
const double DENDRO_26 = 2*alpha[pp];
const double DENDRO_27 = At4[pp]*DENDRO_igt4;
const double DENDRO_28 = 4*alpha[pp];

// Dendro: printing variables
//--
B_rhs2[pp] = -2*(deriv_evars->grad_0_alpha[d_pp])*DENDRO_10 - (deriv_evars->grad_0_beta2[d_pp])*DENDRO_Gtk0 - 2*(deriv_evars->grad_1_alpha[d_pp])*DENDRO_12 - (deriv_evars->grad_1_beta2[d_pp])*DENDRO_Gtk1 - (deriv_evars->grad_2_alpha[d_pp])*DENDRO_15 - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_Gtk2 - B2[pp]*eta - DENDRO_0*lambda[3] + DENDRO_0 + DENDRO_10*DENDRO_28*DENDRO_C2_k2_2 + DENDRO_12*DENDRO_28*DENDRO_C2_k2_4 + DENDRO_15*DENDRO_C2_k2_5*alpha[pp] + DENDRO_16*deriv_evars->grad2_0_1_beta2[d_pp] - DENDRO_18*(4*(deriv_evars->grad_0_K[d_pp])*DENDRO_igt2 + (deriv_evars->grad_0_chi[d_pp])*DENDRO_10*DENDRO_17) - DENDRO_18*(4*(deriv_evars->grad_1_K[d_pp])*DENDRO_igt4 + (deriv_evars->grad_1_chi[d_pp])*DENDRO_12*DENDRO_17) - DENDRO_18*(4*(deriv_evars->grad_2_K[d_pp])*DENDRO_igt5 + (deriv_evars->grad_2_chi[d_pp])*DENDRO_14*DENDRO_17) + DENDRO_19*deriv_evars->grad2_0_0_beta0[d_pp] + DENDRO_19*deriv_evars->grad2_0_1_beta1[d_pp] + DENDRO_20*deriv_evars->grad2_0_1_beta0[d_pp] + DENDRO_20*deriv_evars->grad2_1_1_beta1[d_pp] + DENDRO_21*deriv_evars->grad2_0_2_beta0[d_pp] + DENDRO_21*deriv_evars->grad2_1_2_beta1[d_pp] + DENDRO_26*DENDRO_C2_k2_0*(At0[pp]*(DENDRO_igt0 * DENDRO_igt0) + At3[pp]*DENDRO_22 + At5[pp]*DENDRO_9 + DENDRO_13*DENDRO_24 + DENDRO_16*DENDRO_23 + DENDRO_16*DENDRO_25) + DENDRO_26*DENDRO_C2_k2_3*(At0[pp]*DENDRO_22 + At1[pp]*DENDRO_16*DENDRO_igt3 + At2[pp]*DENDRO_16*DENDRO_igt4 + At3[pp]*(DENDRO_igt3 * DENDRO_igt3) + At5[pp]*DENDRO_11 + 2*DENDRO_27*DENDRO_igt3) + DENDRO_28*DENDRO_C2_k2_1*(At0[pp]*DENDRO_igt0*DENDRO_igt1 + At1[pp]*DENDRO_22 + At2[pp]*DENDRO_igt1*DENDRO_igt2 + At3[pp]*DENDRO_igt1*DENDRO_igt3 + At5[pp]*DENDRO_7 + DENDRO_23*DENDRO_igt3 + DENDRO_24*DENDRO_igt4 + DENDRO_25*DENDRO_igt3 + DENDRO_27*DENDRO_igt1) + (2.0/3.0)*DENDRO_Gtk2*((deriv_evars->grad_0_beta0[d_pp]) + (deriv_evars->grad_1_beta1[d_pp]) + (deriv_evars->grad_2_beta2[d_pp])) + DENDRO_igt0*deriv_evars->grad2_0_0_beta2[d_pp] + (7.0/3.0)*DENDRO_igt2*deriv_evars->grad2_0_2_beta2[d_pp] + DENDRO_igt3*deriv_evars->grad2_1_1_beta2[d_pp] + (7.0/3.0)*DENDRO_igt4*deriv_evars->grad2_1_2_beta2[d_pp] + (4.0/3.0)*DENDRO_igt5*deriv_evars->grad2_2_2_beta2[d_pp] + lambda[2]*((deriv_evars->grad_0_B2[d_pp])*beta0[pp] + (deriv_evars->grad_1_B2[d_pp])*beta1[pp] + (deriv_evars->grad_2_B2[d_pp])*beta2[pp]);
// Dendro: reduced ops: 224
// Dendro: }}} 
}
