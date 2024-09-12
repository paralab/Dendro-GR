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
    const double DENDRO_0 = gt3[pp] * gt5[pp];
    const double DENDRO_1 = pow(gt4[pp], 2);
    const double DENDRO_2 = pow(gt1[pp], 2);
    const double DENDRO_3 = pow(gt2[pp], 2);
    const double DENDRO_4 = gt2[pp] * gt4[pp];
    const double DENDRO_5 =
        1.0 / (-DENDRO_0 * gt0[pp] + DENDRO_1 * gt0[pp] + DENDRO_2 * gt5[pp] +
               DENDRO_3 * gt3[pp] - 2 * DENDRO_4 * gt1[pp]);

    // Dendro: printing variables
    //--
    DENDRO_igt0 = -DENDRO_5 * (DENDRO_0 - DENDRO_1);
    //--
    DENDRO_igt1 = DENDRO_5 * (-DENDRO_4 + gt1[pp] * gt5[pp]);
    //--
    DENDRO_igt2 = -DENDRO_5 * (gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp]);
    //--
    DENDRO_igt3 = -DENDRO_5 * (-DENDRO_3 + gt0[pp] * gt5[pp]);
    //--
    DENDRO_igt4 = DENDRO_5 * (gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp]);
    //--
    DENDRO_igt5 = -DENDRO_5 * (-DENDRO_2 + gt0[pp] * gt3[pp]);
    // Dendro: reduced ops: 39
    // Dendro: }}}
}

{
    // Dendro: {{{
    // Dendro: original ops: 14
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    DENDRO_C1_k0_0 = 0.5 * (deriv_evars->grad_0_gt0[d_pp]);
    //--
    DENDRO_C1_k0_1 = 0.5 * (deriv_evars->grad_1_gt0[d_pp]);
    //--
    DENDRO_C1_k0_2 = 0.5 * (deriv_evars->grad_2_gt0[d_pp]);
    //--
    DENDRO_C1_k0_3 = -0.5 * (deriv_evars->grad_0_gt3[d_pp]) +
                     1.0 * (deriv_evars->grad_1_gt1[d_pp]);
    //--
    DENDRO_C1_k0_4 = 0.5 * (-(deriv_evars->grad_0_gt4[d_pp]) +
                            (deriv_evars->grad_1_gt2[d_pp]) +
                            (deriv_evars->grad_2_gt1[d_pp]));
    //--
    DENDRO_C1_k0_5 = -0.5 * (deriv_evars->grad_0_gt5[d_pp]) +
                     1.0 * (deriv_evars->grad_2_gt2[d_pp]);
    // Dendro: reduced ops: 12
    // Dendro: }}}
}

{
    // Dendro: {{{
    // Dendro: original ops: 14
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    DENDRO_C1_k1_0 = 1.0 * (deriv_evars->grad_0_gt1[d_pp]) -
                     0.5 * (deriv_evars->grad_1_gt0[d_pp]);
    //--
    DENDRO_C1_k1_1 = 0.5 * (deriv_evars->grad_0_gt3[d_pp]);
    //--
    DENDRO_C1_k1_2 = 0.5 * ((deriv_evars->grad_0_gt4[d_pp]) -
                            (deriv_evars->grad_1_gt2[d_pp]) +
                            (deriv_evars->grad_2_gt1[d_pp]));
    //--
    DENDRO_C1_k1_3 = 0.5 * (deriv_evars->grad_1_gt3[d_pp]);
    //--
    DENDRO_C1_k1_4 = 0.5 * (deriv_evars->grad_2_gt3[d_pp]);
    //--
    DENDRO_C1_k1_5 = -0.5 * (deriv_evars->grad_1_gt5[d_pp]) +
                     1.0 * (deriv_evars->grad_2_gt4[d_pp]);
    // Dendro: reduced ops: 12
    // Dendro: }}}
}

{
    // Dendro: {{{
    // Dendro: original ops: 14
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    DENDRO_C1_k2_0 = 1.0 * (deriv_evars->grad_0_gt2[d_pp]) -
                     0.5 * (deriv_evars->grad_2_gt0[d_pp]);
    //--
    DENDRO_C1_k2_1 = 0.5 * ((deriv_evars->grad_0_gt4[d_pp]) +
                            (deriv_evars->grad_1_gt2[d_pp]) -
                            (deriv_evars->grad_2_gt1[d_pp]));
    //--
    DENDRO_C1_k2_2 = 0.5 * (deriv_evars->grad_0_gt5[d_pp]);
    //--
    DENDRO_C1_k2_3 = 1.0 * (deriv_evars->grad_1_gt4[d_pp]) -
                     0.5 * (deriv_evars->grad_2_gt3[d_pp]);
    //--
    DENDRO_C1_k2_4 = 0.5 * (deriv_evars->grad_1_gt5[d_pp]);
    //--
    DENDRO_C1_k2_5 = 0.5 * (deriv_evars->grad_2_gt5[d_pp]);
    // Dendro: reduced ops: 12
    // Dendro: }}}
}

{
    // Dendro: {{{
    // Dendro: original ops: 30
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    DENDRO_C2_k0_0 = DENDRO_C1_k0_0 * DENDRO_igt0 +
                     DENDRO_C1_k1_0 * DENDRO_igt1 +
                     DENDRO_C1_k2_0 * DENDRO_igt2;
    //--
    DENDRO_C2_k0_1 = DENDRO_C1_k0_1 * DENDRO_igt0 +
                     DENDRO_C1_k1_1 * DENDRO_igt1 +
                     DENDRO_C1_k2_1 * DENDRO_igt2;
    //--
    DENDRO_C2_k0_2 = DENDRO_C1_k0_2 * DENDRO_igt0 +
                     DENDRO_C1_k1_2 * DENDRO_igt1 +
                     DENDRO_C1_k2_2 * DENDRO_igt2;
    //--
    DENDRO_C2_k0_3 = DENDRO_C1_k0_3 * DENDRO_igt0 +
                     DENDRO_C1_k1_3 * DENDRO_igt1 +
                     DENDRO_C1_k2_3 * DENDRO_igt2;
    //--
    DENDRO_C2_k0_4 = DENDRO_C1_k0_4 * DENDRO_igt0 +
                     DENDRO_C1_k1_4 * DENDRO_igt1 +
                     DENDRO_C1_k2_4 * DENDRO_igt2;
    //--
    DENDRO_C2_k0_5 = DENDRO_C1_k0_5 * DENDRO_igt0 +
                     DENDRO_C1_k1_5 * DENDRO_igt1 +
                     DENDRO_C1_k2_5 * DENDRO_igt2;
    // Dendro: reduced ops: 30
    // Dendro: }}}
}

{
    // Dendro: {{{
    // Dendro: original ops: 30
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    DENDRO_C2_k1_0 = DENDRO_C1_k0_0 * DENDRO_igt1 +
                     DENDRO_C1_k1_0 * DENDRO_igt3 +
                     DENDRO_C1_k2_0 * DENDRO_igt4;
    //--
    DENDRO_C2_k1_1 = DENDRO_C1_k0_1 * DENDRO_igt1 +
                     DENDRO_C1_k1_1 * DENDRO_igt3 +
                     DENDRO_C1_k2_1 * DENDRO_igt4;
    //--
    DENDRO_C2_k1_2 = DENDRO_C1_k0_2 * DENDRO_igt1 +
                     DENDRO_C1_k1_2 * DENDRO_igt3 +
                     DENDRO_C1_k2_2 * DENDRO_igt4;
    //--
    DENDRO_C2_k1_3 = DENDRO_C1_k0_3 * DENDRO_igt1 +
                     DENDRO_C1_k1_3 * DENDRO_igt3 +
                     DENDRO_C1_k2_3 * DENDRO_igt4;
    //--
    DENDRO_C2_k1_4 = DENDRO_C1_k0_4 * DENDRO_igt1 +
                     DENDRO_C1_k1_4 * DENDRO_igt3 +
                     DENDRO_C1_k2_4 * DENDRO_igt4;
    //--
    DENDRO_C2_k1_5 = DENDRO_C1_k0_5 * DENDRO_igt1 +
                     DENDRO_C1_k1_5 * DENDRO_igt3 +
                     DENDRO_C1_k2_5 * DENDRO_igt4;
    // Dendro: reduced ops: 30
    // Dendro: }}}
}

{
    // Dendro: {{{
    // Dendro: original ops: 30
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    DENDRO_C2_k2_0 = DENDRO_C1_k0_0 * DENDRO_igt2 +
                     DENDRO_C1_k1_0 * DENDRO_igt4 +
                     DENDRO_C1_k2_0 * DENDRO_igt5;
    //--
    DENDRO_C2_k2_1 = DENDRO_C1_k0_1 * DENDRO_igt2 +
                     DENDRO_C1_k1_1 * DENDRO_igt4 +
                     DENDRO_C1_k2_1 * DENDRO_igt5;
    //--
    DENDRO_C2_k2_2 = DENDRO_C1_k0_2 * DENDRO_igt2 +
                     DENDRO_C1_k1_2 * DENDRO_igt4 +
                     DENDRO_C1_k2_2 * DENDRO_igt5;
    //--
    DENDRO_C2_k2_3 = DENDRO_C1_k0_3 * DENDRO_igt2 +
                     DENDRO_C1_k1_3 * DENDRO_igt4 +
                     DENDRO_C1_k2_3 * DENDRO_igt5;
    //--
    DENDRO_C2_k2_4 = DENDRO_C1_k0_4 * DENDRO_igt2 +
                     DENDRO_C1_k1_4 * DENDRO_igt4 +
                     DENDRO_C1_k2_4 * DENDRO_igt5;
    //--
    DENDRO_C2_k2_5 = DENDRO_C1_k0_5 * DENDRO_igt2 +
                     DENDRO_C1_k1_5 * DENDRO_igt4 +
                     DENDRO_C1_k2_5 * DENDRO_igt5;
    // Dendro: reduced ops: 30
    // Dendro: }}}
}
// Dendro: {{{
// Dendro: original ops: 413
// Dendro: printing temp variables
const double DENDRO_0 = pow(chi[pp], 2);
const double DENDRO_1 =
    (deriv_evars->grad_0_chi[d_pp]) * (deriv_evars->grad_2_chi[d_pp]);
const double DENDRO_2 = 2 * chi[pp];
const double DENDRO_3 = 3 * (deriv_evars->grad_1_chi[d_pp]);
const double DENDRO_4 = 2 * DENDRO_igt1;
const double DENDRO_5 = 2 * DENDRO_igt2;
const double DENDRO_6 = 2 * DENDRO_igt4;
const double DENDRO_7 =
    DENDRO_4 * DENDRO_C2_k0_1 + DENDRO_5 * DENDRO_C2_k0_2 +
    DENDRO_6 * DENDRO_C2_k0_4 + DENDRO_C2_k0_0 * DENDRO_igt0 +
    DENDRO_C2_k0_3 * DENDRO_igt3 + DENDRO_C2_k0_5 * DENDRO_igt5;
const double DENDRO_8 =
    DENDRO_4 * DENDRO_C2_k1_1 + DENDRO_5 * DENDRO_C2_k1_2 +
    DENDRO_6 * DENDRO_C2_k1_4 + DENDRO_C2_k1_0 * DENDRO_igt0 +
    DENDRO_C2_k1_3 * DENDRO_igt3 + DENDRO_C2_k1_5 * DENDRO_igt5;
const double DENDRO_9 =
    DENDRO_4 * DENDRO_C2_k2_1 + DENDRO_5 * DENDRO_C2_k2_2 +
    DENDRO_6 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
    DENDRO_C2_k2_3 * DENDRO_igt3 + DENDRO_C2_k2_5 * DENDRO_igt5;
const double DENDRO_10 = 2.0 * gt2[pp];
const double DENDRO_11 = 4 * DENDRO_igt0;
const double DENDRO_12 = 4 * DENDRO_igt2;
const double DENDRO_13 = 4 * DENDRO_igt3;
const double DENDRO_14 = 4 * DENDRO_igt5;
const double DENDRO_15 = 4 * DENDRO_igt1;
const double DENDRO_16 =
    DENDRO_C1_k2_2 * DENDRO_C2_k2_1 + DENDRO_C1_k2_4 * DENDRO_C2_k2_0;
const double DENDRO_17 =
    DENDRO_C1_k2_2 * DENDRO_C2_k2_2 + DENDRO_C1_k2_5 * DENDRO_C2_k2_0;
const double DENDRO_18 = 4 * DENDRO_igt4;
const double DENDRO_19 =
    DENDRO_C1_k2_4 * DENDRO_C2_k2_2 + DENDRO_C1_k2_5 * DENDRO_C2_k2_1;

// Dendro: printing variables
//--
At_rhs02[pp] =
    -1.0 / 4.0 *
    (-DENDRO_0 * ((deriv_evars->grad_0_Gt0[d_pp]) * DENDRO_10 +
                  2.0 * (deriv_evars->grad_0_Gt1[d_pp]) * gt4[pp] +
                  2.0 * (deriv_evars->grad_0_Gt2[d_pp]) * gt5[pp] +
                  2.0 * (deriv_evars->grad_2_Gt0[d_pp]) * gt0[pp] +
                  2.0 * (deriv_evars->grad_2_Gt1[d_pp]) * gt1[pp] +
                  (deriv_evars->grad_2_Gt2[d_pp]) * DENDRO_10 +
                  DENDRO_11 * (DENDRO_C1_k0_2 * DENDRO_C2_k2_2 +
                               2 * DENDRO_C1_k2_2 * DENDRO_C2_k2_0) +
                  DENDRO_11 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k0_2 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_0) +
                  DENDRO_11 * (DENDRO_C1_k0_1 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_0) +
                  DENDRO_12 * (DENDRO_17 + DENDRO_C1_k0_2 * DENDRO_C2_k2_5) +
                  DENDRO_12 * (DENDRO_17 + DENDRO_C1_k0_5 * DENDRO_C2_k2_2) +
                  DENDRO_12 * (2 * DENDRO_C1_k0_2 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_0) +
                  DENDRO_12 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_2) +
                  DENDRO_12 * (DENDRO_C1_k0_1 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_2) +
                  DENDRO_12 * (DENDRO_C1_k0_4 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_0) +
                  DENDRO_13 * (DENDRO_C1_k0_4 * DENDRO_C2_k2_4 +
                               2 * DENDRO_C1_k2_4 * DENDRO_C2_k2_1) +
                  DENDRO_13 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_1) +
                  DENDRO_13 * (DENDRO_C1_k0_3 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_1) +
                  DENDRO_14 * (DENDRO_C1_k0_5 * DENDRO_C2_k2_5 +
                               2 * DENDRO_C1_k2_5 * DENDRO_C2_k2_2) +
                  DENDRO_14 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_2) +
                  DENDRO_14 * (DENDRO_C1_k0_4 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_2) +
                  DENDRO_15 * (DENDRO_16 + DENDRO_C1_k0_2 * DENDRO_C2_k2_4) +
                  DENDRO_15 * (DENDRO_16 + DENDRO_C1_k0_4 * DENDRO_C2_k2_2) +
                  DENDRO_15 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_1) +
                  DENDRO_15 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k0_2 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_0) +
                  DENDRO_15 * (DENDRO_C1_k0_1 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_1) +
                  DENDRO_15 * (DENDRO_C1_k0_3 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_0) +
                  DENDRO_18 * (DENDRO_19 + DENDRO_C1_k0_4 * DENDRO_C2_k2_5) +
                  DENDRO_18 * (DENDRO_19 + DENDRO_C1_k0_5 * DENDRO_C2_k2_4) +
                  DENDRO_18 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_2) +
                  DENDRO_18 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_1) +
                  DENDRO_18 * (DENDRO_C1_k0_3 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_2) +
                  DENDRO_18 * (DENDRO_C1_k0_4 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_1) +
                  2.0 * DENDRO_7 * (DENDRO_C1_k0_2 + DENDRO_C1_k2_0) +
                  2.0 * DENDRO_8 * (DENDRO_C1_k0_4 + DENDRO_C1_k2_1) +
                  2.0 * DENDRO_9 * (DENDRO_C1_k0_5 + DENDRO_C1_k2_2) -
                  2.0 * DENDRO_igt0 * deriv_evars->grad2_0_0_gt2[d_pp] -
                  4.0 * DENDRO_igt1 * deriv_evars->grad2_0_1_gt2[d_pp] -
                  4.0 * DENDRO_igt2 * deriv_evars->grad2_0_2_gt2[d_pp] -
                  2.0 * DENDRO_igt3 * deriv_evars->grad2_1_1_gt2[d_pp] -
                  4.0 * DENDRO_igt4 * deriv_evars->grad2_1_2_gt2[d_pp] -
                  2.0 * DENDRO_igt5 * deriv_evars->grad2_2_2_gt2[d_pp]) +
     DENDRO_1 +
     DENDRO_2 * ((deriv_evars->grad_0_chi[d_pp]) * DENDRO_C2_k0_2 +
                 (deriv_evars->grad_1_chi[d_pp]) * DENDRO_C2_k1_2 +
                 (deriv_evars->grad_2_chi[d_pp]) * DENDRO_C2_k2_2 -
                 deriv_evars->grad2_0_2_chi[d_pp]) +
     gt2[pp] * (DENDRO_2 * ((deriv_evars->grad_0_chi[d_pp]) * DENDRO_7 +
                            (deriv_evars->grad_1_chi[d_pp]) * DENDRO_8 +
                            (deriv_evars->grad_2_chi[d_pp]) * DENDRO_9) +
                DENDRO_4 * ((deriv_evars->grad_0_chi[d_pp]) * DENDRO_3 -
                            DENDRO_2 * deriv_evars->grad2_0_1_chi[d_pp]) +
                DENDRO_5 * (3 * DENDRO_1 -
                            DENDRO_2 * deriv_evars->grad2_0_2_chi[d_pp]) +
                DENDRO_6 * ((deriv_evars->grad_2_chi[d_pp]) * DENDRO_3 -
                            DENDRO_2 * deriv_evars->grad2_1_2_chi[d_pp]) +
                DENDRO_igt0 * (3 * pow((deriv_evars->grad_0_chi[d_pp]), 2) -
                               DENDRO_2 * deriv_evars->grad2_0_0_chi[d_pp]) +
                DENDRO_igt3 * (3 * pow((deriv_evars->grad_1_chi[d_pp]), 2) -
                               DENDRO_2 * deriv_evars->grad2_1_1_chi[d_pp]) +
                DENDRO_igt5 * (3 * pow((deriv_evars->grad_2_chi[d_pp]), 2) -
                               DENDRO_2 * deriv_evars->grad2_2_2_chi[d_pp]))) /
    DENDRO_0;
// Dendro: reduced ops: 322
// Dendro: }}}
