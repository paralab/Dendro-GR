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
// Dendro: original ops: 51 
// Dendro: printing temp variables
const double DENDRO_0 = 2*At1[pp];
const double DENDRO_1 = 2*At2[pp];
const double DENDRO_2 = (2.0/3.0)*At0[pp];

// Dendro: printing variables
//--
At_rhs00[pp] = (deriv_evars->grad_0_At0[d_pp])*beta0[pp] + (4.0/3.0)*(deriv_evars->grad_0_beta0[d_pp])*At0[pp] + (deriv_evars->grad_0_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta2[d_pp])*DENDRO_1 + (deriv_evars->grad_1_At0[d_pp])*beta1[pp] - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_2 + (deriv_evars->grad_2_At0[d_pp])*beta2[pp] - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_2 + At_rhs00[pp] - alpha[pp]*(-At0[pp]*K[pp] + 2*At0[pp]*(At0[pp]*DENDRO_igt0 + At1[pp]*DENDRO_igt1 + At2[pp]*DENDRO_igt2) + DENDRO_0*(At0[pp]*DENDRO_igt1 + At1[pp]*DENDRO_igt3 + At2[pp]*DENDRO_igt4) + DENDRO_1*(At0[pp]*DENDRO_igt2 + At1[pp]*DENDRO_igt4 + At2[pp]*DENDRO_igt5));
// Dendro: reduced ops: 47
// Dendro: }}} 

}

{
// Dendro: {{{ 
// Dendro: original ops: 51 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*At1[pp];

// Dendro: printing variables
//--
At_rhs01[pp] = (deriv_evars->grad_0_At1[d_pp])*beta0[pp] + (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta1[d_pp])*At3[pp] + (deriv_evars->grad_0_beta2[d_pp])*At4[pp] + (deriv_evars->grad_1_At1[d_pp])*beta1[pp] + (deriv_evars->grad_1_beta0[d_pp])*At0[pp] + (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_beta2[d_pp])*At2[pp] + (deriv_evars->grad_2_At1[d_pp])*beta2[pp] - 2.0/3.0*(deriv_evars->grad_2_beta2[d_pp])*At1[pp] + At_rhs01[pp] - alpha[pp]*(2*At0[pp]*(At1[pp]*DENDRO_igt0 + At3[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt2) - At1[pp]*K[pp] + 2*At1[pp]*(At1[pp]*DENDRO_igt1 + At3[pp]*DENDRO_igt3 + At4[pp]*DENDRO_igt4) + 2*At2[pp]*(At1[pp]*DENDRO_igt2 + At3[pp]*DENDRO_igt4 + At4[pp]*DENDRO_igt5));
// Dendro: reduced ops: 50
// Dendro: }}} 

}

{
// Dendro: {{{ 
// Dendro: original ops: 51 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*At2[pp];

// Dendro: printing variables
//--
At_rhs02[pp] = (deriv_evars->grad_0_At2[d_pp])*beta0[pp] + (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_0_beta1[d_pp])*At4[pp] + (deriv_evars->grad_0_beta2[d_pp])*At5[pp] + (deriv_evars->grad_1_At2[d_pp])*beta1[pp] - 2.0/3.0*(deriv_evars->grad_1_beta1[d_pp])*At2[pp] + (deriv_evars->grad_2_At2[d_pp])*beta2[pp] + (deriv_evars->grad_2_beta0[d_pp])*At0[pp] + (deriv_evars->grad_2_beta1[d_pp])*At1[pp] + (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 + At_rhs02[pp] - alpha[pp]*(2*At0[pp]*(At2[pp]*DENDRO_igt0 + At4[pp]*DENDRO_igt1 + At5[pp]*DENDRO_igt2) + 2*At1[pp]*(At2[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt3 + At5[pp]*DENDRO_igt4) - At2[pp]*K[pp] + 2*At2[pp]*(At2[pp]*DENDRO_igt2 + At4[pp]*DENDRO_igt4 + At5[pp]*DENDRO_igt5));
// Dendro: reduced ops: 50
// Dendro: }}} 

}

{
// Dendro: {{{ 
// Dendro: original ops: 51 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*At3[pp];
const double DENDRO_1 = 2*At1[pp];
const double DENDRO_2 = 2*At4[pp];

// Dendro: printing variables
//--
At_rhs11[pp] = (deriv_evars->grad_0_At3[d_pp])*beta0[pp] - (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_1_At3[d_pp])*beta1[pp] + (deriv_evars->grad_1_beta0[d_pp])*DENDRO_1 + (4.0/3.0)*(deriv_evars->grad_1_beta1[d_pp])*At3[pp] + (deriv_evars->grad_1_beta2[d_pp])*DENDRO_2 + (deriv_evars->grad_2_At3[d_pp])*beta2[pp] - (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 + At_rhs11[pp] - alpha[pp]*(-At3[pp]*K[pp] + 2*At3[pp]*(At1[pp]*DENDRO_igt1 + At3[pp]*DENDRO_igt3 + At4[pp]*DENDRO_igt4) + DENDRO_1*(At1[pp]*DENDRO_igt0 + At3[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt2) + DENDRO_2*(At1[pp]*DENDRO_igt2 + At3[pp]*DENDRO_igt4 + At4[pp]*DENDRO_igt5));
// Dendro: reduced ops: 47
// Dendro: }}} 

}

{
// Dendro: {{{ 
// Dendro: original ops: 51 
// Dendro: printing temp variables
const double DENDRO_0 = (1.0/3.0)*At4[pp];

// Dendro: printing variables
//--
At_rhs12[pp] = (deriv_evars->grad_0_At4[d_pp])*beta0[pp] - 2.0/3.0*(deriv_evars->grad_0_beta0[d_pp])*At4[pp] + (deriv_evars->grad_1_At4[d_pp])*beta1[pp] + (deriv_evars->grad_1_beta0[d_pp])*At2[pp] + (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_1_beta2[d_pp])*At5[pp] + (deriv_evars->grad_2_At4[d_pp])*beta2[pp] + (deriv_evars->grad_2_beta0[d_pp])*At1[pp] + (deriv_evars->grad_2_beta1[d_pp])*At3[pp] + (deriv_evars->grad_2_beta2[d_pp])*DENDRO_0 + At_rhs12[pp] - alpha[pp]*(2*At1[pp]*(At2[pp]*DENDRO_igt0 + At4[pp]*DENDRO_igt1 + At5[pp]*DENDRO_igt2) + 2*At3[pp]*(At2[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt3 + At5[pp]*DENDRO_igt4) - At4[pp]*K[pp] + 2*At4[pp]*(At2[pp]*DENDRO_igt2 + At4[pp]*DENDRO_igt4 + At5[pp]*DENDRO_igt5));
// Dendro: reduced ops: 50
// Dendro: }}} 

}

{
// Dendro: {{{ 
// Dendro: original ops: 51 
// Dendro: printing temp variables
const double DENDRO_0 = (2.0/3.0)*At5[pp];
const double DENDRO_1 = 2*At2[pp];
const double DENDRO_2 = 2*At4[pp];

// Dendro: printing variables
//--
At_rhs22[pp] = (deriv_evars->grad_0_At5[d_pp])*beta0[pp] - (deriv_evars->grad_0_beta0[d_pp])*DENDRO_0 + (deriv_evars->grad_1_At5[d_pp])*beta1[pp] - (deriv_evars->grad_1_beta1[d_pp])*DENDRO_0 + (deriv_evars->grad_2_At5[d_pp])*beta2[pp] + (deriv_evars->grad_2_beta0[d_pp])*DENDRO_1 + (deriv_evars->grad_2_beta1[d_pp])*DENDRO_2 + (4.0/3.0)*(deriv_evars->grad_2_beta2[d_pp])*At5[pp] + At_rhs22[pp] - alpha[pp]*(-At5[pp]*K[pp] + 2*At5[pp]*(At2[pp]*DENDRO_igt2 + At4[pp]*DENDRO_igt4 + At5[pp]*DENDRO_igt5) + DENDRO_1*(At2[pp]*DENDRO_igt0 + At4[pp]*DENDRO_igt1 + At5[pp]*DENDRO_igt2) + DENDRO_2*(At2[pp]*DENDRO_igt1 + At4[pp]*DENDRO_igt3 + At5[pp]*DENDRO_igt4));
// Dendro: reduced ops: 47
// Dendro: }}} 

}
