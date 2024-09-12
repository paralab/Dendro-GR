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
{
    // Dendro: {{{
    // Dendro: original ops: 60
    // Dendro: printing temp variables
    const double DENDRO_0 = 1.0 / chi[pp];
    const double DENDRO_1 =
        0.5 * (deriv_evars->grad_0_chi[d_pp]) * DENDRO_igt0 +
        0.5 * (deriv_evars->grad_1_chi[d_pp]) * DENDRO_igt1 +
        0.5 * (deriv_evars->grad_2_chi[d_pp]) * DENDRO_igt2;
    const double DENDRO_2 = DENDRO_0 * DENDRO_1;

    // Dendro: printing variables
    //--
    DENDRO_C3_k0_0        = -DENDRO_0 * (1.0 * (deriv_evars->grad_0_chi[d_pp]) -
                                  DENDRO_1 * gt0[pp]) +
                     DENDRO_C2_k0_0;
    //--
    DENDRO_C3_k0_1 = -DENDRO_0 * (0.5 * (deriv_evars->grad_1_chi[d_pp]) -
                                  DENDRO_1 * gt1[pp]) +
                     DENDRO_C2_k0_1;
    //--
    DENDRO_C3_k0_2 = -DENDRO_0 * (0.5 * (deriv_evars->grad_2_chi[d_pp]) -
                                  DENDRO_1 * gt2[pp]) +
                     DENDRO_C2_k0_2;
    //--
    DENDRO_C3_k0_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k0_3;
    //--
    DENDRO_C3_k0_4 = DENDRO_2 * gt4[pp] + DENDRO_C2_k0_4;
    //--
    DENDRO_C3_k0_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k0_5;
    // Dendro: reduced ops: 31
    // Dendro: }}}
}
{
    // Dendro: {{{
    // Dendro: original ops: 60
    // Dendro: printing temp variables
    const double DENDRO_0 = 1.0 / chi[pp];
    const double DENDRO_1 =
        0.5 * (deriv_evars->grad_0_chi[d_pp]) * DENDRO_igt1 +
        0.5 * (deriv_evars->grad_1_chi[d_pp]) * DENDRO_igt3 +
        0.5 * (deriv_evars->grad_2_chi[d_pp]) * DENDRO_igt4;
    const double DENDRO_2 = DENDRO_0 * DENDRO_1;

    // Dendro: printing variables
    //--
    DENDRO_C3_k1_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k1_0;
    //--
    DENDRO_C3_k1_1        = -DENDRO_0 * (0.5 * (deriv_evars->grad_0_chi[d_pp]) -
                                  DENDRO_1 * gt1[pp]) +
                     DENDRO_C2_k1_1;
    //--
    DENDRO_C3_k1_2 = DENDRO_2 * gt2[pp] + DENDRO_C2_k1_2;
    //--
    DENDRO_C3_k1_3 = -DENDRO_0 * (1.0 * (deriv_evars->grad_1_chi[d_pp]) -
                                  DENDRO_1 * gt3[pp]) +
                     DENDRO_C2_k1_3;
    //--
    DENDRO_C3_k1_4 = -DENDRO_0 * (0.5 * (deriv_evars->grad_2_chi[d_pp]) -
                                  DENDRO_1 * gt4[pp]) +
                     DENDRO_C2_k1_4;
    //--
    DENDRO_C3_k1_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k1_5;
    // Dendro: reduced ops: 31
    // Dendro: }}}
}
{
    // Dendro: {{{
    // Dendro: original ops: 60
    // Dendro: printing temp variables
    const double DENDRO_0 = 1.0 / chi[pp];
    const double DENDRO_1 =
        0.5 * (deriv_evars->grad_0_chi[d_pp]) * DENDRO_igt2 +
        0.5 * (deriv_evars->grad_1_chi[d_pp]) * DENDRO_igt4 +
        0.5 * (deriv_evars->grad_2_chi[d_pp]) * DENDRO_igt5;
    const double DENDRO_2 = DENDRO_0 * DENDRO_1;

    // Dendro: printing variables
    //--
    DENDRO_C3_k2_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k2_0;
    //--
    DENDRO_C3_k2_1        = DENDRO_2 * gt1[pp] + DENDRO_C2_k2_1;
    //--
    DENDRO_C3_k2_2        = -DENDRO_0 * (0.5 * (deriv_evars->grad_0_chi[d_pp]) -
                                  DENDRO_1 * gt2[pp]) +
                     DENDRO_C2_k2_2;
    //--
    DENDRO_C3_k2_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k2_3;
    //--
    DENDRO_C3_k2_4 = -DENDRO_0 * (0.5 * (deriv_evars->grad_1_chi[d_pp]) -
                                  DENDRO_1 * gt4[pp]) +
                     DENDRO_C2_k2_4;
    //--
    DENDRO_C3_k2_5 = -DENDRO_0 * (1.0 * (deriv_evars->grad_2_chi[d_pp]) -
                                  DENDRO_1 * gt5[pp]) +
                     DENDRO_C2_k2_5;
    // Dendro: reduced ops: 31
    // Dendro: }}}
}
{
    // Dendro: {{{
    // Dendro: original ops: 42
    // Dendro: printing temp variables
    const double DENDRO_0 = 2 * DENDRO_igt1;
    const double DENDRO_1 = 2 * DENDRO_igt2;
    const double DENDRO_2 = 2 * DENDRO_igt4;

    // Dendro: printing variables
    //--
    DENDRO_Gtk0 = DENDRO_0 * DENDRO_C2_k0_1 + DENDRO_1 * DENDRO_C2_k0_2 +
                  DENDRO_2 * DENDRO_C2_k0_4 + DENDRO_C2_k0_0 * DENDRO_igt0 +
                  DENDRO_C2_k0_3 * DENDRO_igt3 + DENDRO_C2_k0_5 * DENDRO_igt5;
    //--
    DENDRO_Gtk1 = DENDRO_0 * DENDRO_C2_k1_1 + DENDRO_1 * DENDRO_C2_k1_2 +
                  DENDRO_2 * DENDRO_C2_k1_4 + DENDRO_C2_k1_0 * DENDRO_igt0 +
                  DENDRO_C2_k1_3 * DENDRO_igt3 + DENDRO_C2_k1_5 * DENDRO_igt5;
    //--
    DENDRO_Gtk2 = DENDRO_0 * DENDRO_C2_k2_1 + DENDRO_1 * DENDRO_C2_k2_2 +
                  DENDRO_2 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
                  DENDRO_C2_k2_3 * DENDRO_igt3 + DENDRO_C2_k2_5 * DENDRO_igt5;
    // Dendro: reduced ops: 36
    // Dendro: }}}
}

{
    // Dendro: {{{
    // Dendro: original ops: 2532
    // Dendro: printing temp variables
    const double DENDRO_0 = (DENDRO_igt1 * DENDRO_igt1);
    const double DENDRO_1 = (DENDRO_igt2 * DENDRO_igt2);
    const double DENDRO_2 = At1[pp] * DENDRO_igt0;
    const double DENDRO_3 = 2 * DENDRO_igt1;
    const double DENDRO_4 = At2[pp] * DENDRO_igt0;
    const double DENDRO_5 = 2 * DENDRO_igt2;
    const double DENDRO_6 = DENDRO_igt1 * DENDRO_igt2;
    const double DENDRO_7 = At0[pp] * (DENDRO_igt0 * DENDRO_igt0) +
                            At3[pp] * DENDRO_0 + 2 * At4[pp] * DENDRO_6 +
                            At5[pp] * DENDRO_1 + DENDRO_2 * DENDRO_3 +
                            DENDRO_4 * DENDRO_5;
    const double DENDRO_8  = 2 * (deriv_evars->grad_0_alpha[d_pp]);
    const double DENDRO_9  = At0[pp] * DENDRO_igt0;
    const double DENDRO_10 = At3[pp] * DENDRO_igt1;
    const double DENDRO_11 = At4[pp] * DENDRO_igt1;
    const double DENDRO_12 = At4[pp] * DENDRO_igt2;
    const double DENDRO_13 = At5[pp] * DENDRO_igt2;
    const double DENDRO_14 = At1[pp] * DENDRO_0 + At2[pp] * DENDRO_6 +
                             DENDRO_10 * DENDRO_igt3 + DENDRO_11 * DENDRO_igt4 +
                             DENDRO_12 * DENDRO_igt3 + DENDRO_13 * DENDRO_igt4 +
                             DENDRO_2 * DENDRO_igt3 + DENDRO_4 * DENDRO_igt4 +
                             DENDRO_9 * DENDRO_igt1;
    const double DENDRO_15 = 2 * (deriv_evars->grad_1_alpha[d_pp]);
    const double DENDRO_16 = At1[pp] * DENDRO_6 + At2[pp] * DENDRO_1 +
                             DENDRO_10 * DENDRO_igt4 + DENDRO_11 * DENDRO_igt5 +
                             DENDRO_12 * DENDRO_igt4 + DENDRO_13 * DENDRO_igt5 +
                             DENDRO_2 * DENDRO_igt4 + DENDRO_4 * DENDRO_igt5 +
                             DENDRO_9 * DENDRO_igt2;
    const double DENDRO_17 = 2 * (deriv_evars->grad_2_alpha[d_pp]);
    const double DENDRO_18 = 2 * DENDRO_igt4;
    const double DENDRO_19 = 4 * (deriv_evars->grad_0_K[d_pp]);
    const double DENDRO_20 = 9 / chi[pp];
    const double DENDRO_21 = (deriv_evars->grad_0_chi[d_pp]) * DENDRO_20;
    const double DENDRO_22 = (1.0 / 3.0) * alpha[pp];
    const double DENDRO_23 = 4 * (deriv_evars->grad_1_K[d_pp]);
    const double DENDRO_24 = (deriv_evars->grad_1_chi[d_pp]) * DENDRO_20;
    const double DENDRO_25 = 4 * (deriv_evars->grad_2_K[d_pp]);
    const double DENDRO_26 = (deriv_evars->grad_2_chi[d_pp]) * DENDRO_20;
    const double DENDRO_27 = (1.0 / 3.0) * DENDRO_igt0;
    const double DENDRO_28 = (1.0 / 3.0) * DENDRO_igt1;
    const double DENDRO_29 = (1.0 / 3.0) * DENDRO_igt2;
    const double DENDRO_30 = (2.0 / 3.0) * (deriv_evars->grad_0_beta0[d_pp]) +
                             (2.0 / 3.0) * (deriv_evars->grad_1_beta1[d_pp]) +
                             (2.0 / 3.0) * (deriv_evars->grad_2_beta2[d_pp]);
    const double DENDRO_31 = (7.0 / 3.0) * DENDRO_igt1;
    const double DENDRO_32 = (7.0 / 3.0) * DENDRO_igt2;
    const double DENDRO_33 = 2 * alpha[pp];
    const double DENDRO_34 = DENDRO_33 * DENDRO_7;
    const double DENDRO_35 = (DENDRO_igt4 * DENDRO_igt4);
    const double DENDRO_36 = At1[pp] * DENDRO_igt1;
    const double DENDRO_37 = At2[pp] * DENDRO_igt1;
    const double DENDRO_38 = At4[pp] * DENDRO_igt3;
    const double DENDRO_39 =
        At0[pp] * DENDRO_0 + At3[pp] * (DENDRO_igt3 * DENDRO_igt3) +
        At5[pp] * DENDRO_35 + DENDRO_18 * DENDRO_37 + DENDRO_18 * DENDRO_38 +
        2 * DENDRO_36 * DENDRO_igt3;
    const double DENDRO_40 = DENDRO_33 * DENDRO_39;
    const double DENDRO_41 = At1[pp] * DENDRO_igt2;
    const double DENDRO_42 = At2[pp] * DENDRO_igt2;
    const double DENDRO_43 = At0[pp] * DENDRO_1 + At3[pp] * DENDRO_35 +
                             At4[pp] * DENDRO_18 * DENDRO_igt5 +
                             At5[pp] * (DENDRO_igt5 * DENDRO_igt5) +
                             DENDRO_18 * DENDRO_41 +
                             2 * DENDRO_42 * DENDRO_igt5;
    const double DENDRO_44 = DENDRO_33 * DENDRO_43;
    const double DENDRO_45 = 4 * alpha[pp];
    const double DENDRO_46 = DENDRO_14 * DENDRO_45;
    const double DENDRO_47 = DENDRO_16 * DENDRO_45;
    const double DENDRO_48 =
        At0[pp] * DENDRO_6 + At3[pp] * DENDRO_igt3 * DENDRO_igt4 +
        At4[pp] * DENDRO_35 + At5[pp] * DENDRO_igt4 * DENDRO_igt5 +
        DENDRO_36 * DENDRO_igt4 + DENDRO_37 * DENDRO_igt5 +
        DENDRO_38 * DENDRO_igt5 + DENDRO_41 * DENDRO_igt3 +
        DENDRO_42 * DENDRO_igt4;
    const double DENDRO_49 = DENDRO_45 * DENDRO_48;
    const double DENDRO_50 = (deriv_evars->grad_0_Gt0[d_pp]) * beta0[pp] +
                             (deriv_evars->grad_1_Gt0[d_pp]) * beta1[pp] +
                             (deriv_evars->grad_2_Gt0[d_pp]) * beta2[pp];
    const double DENDRO_51 =
        -(deriv_evars->grad_0_beta0[d_pp]) * DENDRO_Gtk0 -
        (deriv_evars->grad_1_beta0[d_pp]) * DENDRO_Gtk1 -
        (deriv_evars->grad_2_beta0[d_pp]) * DENDRO_Gtk2 -
        DENDRO_14 * DENDRO_15 - DENDRO_16 * DENDRO_17 +
        DENDRO_18 * deriv_evars->grad2_1_2_beta0[d_pp] -
        DENDRO_22 * (DENDRO_14 * DENDRO_24 + DENDRO_23 * DENDRO_igt1) -
        DENDRO_22 * (DENDRO_16 * DENDRO_26 + DENDRO_25 * DENDRO_igt2) -
        DENDRO_22 * (DENDRO_19 * DENDRO_igt0 + DENDRO_21 * DENDRO_7) +
        DENDRO_27 * deriv_evars->grad2_0_1_beta1[d_pp] +
        DENDRO_27 * deriv_evars->grad2_0_2_beta2[d_pp] +
        DENDRO_28 * deriv_evars->grad2_1_1_beta1[d_pp] +
        DENDRO_28 * deriv_evars->grad2_1_2_beta2[d_pp] +
        DENDRO_29 * deriv_evars->grad2_1_2_beta1[d_pp] +
        DENDRO_29 * deriv_evars->grad2_2_2_beta2[d_pp] +
        DENDRO_30 * DENDRO_Gtk0 +
        DENDRO_31 * deriv_evars->grad2_0_1_beta0[d_pp] +
        DENDRO_32 * deriv_evars->grad2_0_2_beta0[d_pp] +
        DENDRO_34 * DENDRO_C2_k0_0 + DENDRO_40 * DENDRO_C2_k0_3 +
        DENDRO_44 * DENDRO_C2_k0_5 + DENDRO_46 * DENDRO_C2_k0_1 +
        DENDRO_47 * DENDRO_C2_k0_2 + DENDRO_49 * DENDRO_C2_k0_4 + DENDRO_50 -
        DENDRO_7 * DENDRO_8 +
        (4.0 / 3.0) * DENDRO_igt0 * deriv_evars->grad2_0_0_beta0[d_pp] +
        DENDRO_igt3 * deriv_evars->grad2_1_1_beta0[d_pp] +
        DENDRO_igt5 * deriv_evars->grad2_2_2_beta0[d_pp];
    const double DENDRO_52 = (1.0 / 3.0) * DENDRO_igt3;
    const double DENDRO_53 = (1.0 / 3.0) * DENDRO_igt4;
    const double DENDRO_54 = (7.0 / 3.0) * DENDRO_igt4;
    const double DENDRO_55 = (deriv_evars->grad_0_Gt1[d_pp]) * beta0[pp] +
                             (deriv_evars->grad_1_Gt1[d_pp]) * beta1[pp] +
                             (deriv_evars->grad_2_Gt1[d_pp]) * beta2[pp];
    const double DENDRO_56 =
        -(deriv_evars->grad_0_beta1[d_pp]) * DENDRO_Gtk0 -
        (deriv_evars->grad_1_beta1[d_pp]) * DENDRO_Gtk1 -
        (deriv_evars->grad_2_beta1[d_pp]) * DENDRO_Gtk2 - DENDRO_14 * DENDRO_8 -
        DENDRO_15 * DENDRO_39 - DENDRO_17 * DENDRO_48 -
        DENDRO_22 * (DENDRO_14 * DENDRO_21 + DENDRO_19 * DENDRO_igt1) -
        DENDRO_22 * (DENDRO_23 * DENDRO_igt3 + DENDRO_24 * DENDRO_39) -
        DENDRO_22 * (DENDRO_25 * DENDRO_igt4 + DENDRO_26 * DENDRO_48) +
        DENDRO_28 * deriv_evars->grad2_0_0_beta0[d_pp] +
        DENDRO_28 * deriv_evars->grad2_0_2_beta2[d_pp] +
        DENDRO_30 * DENDRO_Gtk1 +
        DENDRO_31 * deriv_evars->grad2_0_1_beta1[d_pp] +
        DENDRO_34 * DENDRO_C2_k1_0 + DENDRO_40 * DENDRO_C2_k1_3 +
        DENDRO_44 * DENDRO_C2_k1_5 + DENDRO_46 * DENDRO_C2_k1_1 +
        DENDRO_47 * DENDRO_C2_k1_2 + DENDRO_49 * DENDRO_C2_k1_4 +
        DENDRO_5 * deriv_evars->grad2_0_2_beta1[d_pp] +
        DENDRO_52 * deriv_evars->grad2_0_1_beta0[d_pp] +
        DENDRO_52 * deriv_evars->grad2_1_2_beta2[d_pp] +
        DENDRO_53 * deriv_evars->grad2_0_2_beta0[d_pp] +
        DENDRO_53 * deriv_evars->grad2_2_2_beta2[d_pp] +
        DENDRO_54 * deriv_evars->grad2_1_2_beta1[d_pp] + DENDRO_55 +
        DENDRO_igt0 * deriv_evars->grad2_0_0_beta1[d_pp] +
        (4.0 / 3.0) * DENDRO_igt3 * deriv_evars->grad2_1_1_beta1[d_pp] +
        DENDRO_igt5 * deriv_evars->grad2_2_2_beta1[d_pp];
    const double DENDRO_57 = (1.0 / 3.0) * DENDRO_igt5;
    const double DENDRO_58 = (deriv_evars->grad_0_Gt2[d_pp]) * beta0[pp] +
                             (deriv_evars->grad_1_Gt2[d_pp]) * beta1[pp] +
                             (deriv_evars->grad_2_Gt2[d_pp]) * beta2[pp];
    const double DENDRO_59 =
        -(deriv_evars->grad_0_beta2[d_pp]) * DENDRO_Gtk0 -
        (deriv_evars->grad_1_beta2[d_pp]) * DENDRO_Gtk1 -
        (deriv_evars->grad_2_beta2[d_pp]) * DENDRO_Gtk2 -
        DENDRO_15 * DENDRO_48 - DENDRO_16 * DENDRO_8 - DENDRO_17 * DENDRO_43 -
        DENDRO_22 * (DENDRO_16 * DENDRO_21 + DENDRO_19 * DENDRO_igt2) -
        DENDRO_22 * (DENDRO_23 * DENDRO_igt4 + DENDRO_24 * DENDRO_48) -
        DENDRO_22 * (DENDRO_25 * DENDRO_igt5 + DENDRO_26 * DENDRO_43) +
        DENDRO_29 * deriv_evars->grad2_0_0_beta0[d_pp] +
        DENDRO_29 * deriv_evars->grad2_0_1_beta1[d_pp] +
        DENDRO_3 * deriv_evars->grad2_0_1_beta2[d_pp] +
        DENDRO_30 * DENDRO_Gtk2 +
        DENDRO_32 * deriv_evars->grad2_0_2_beta2[d_pp] +
        DENDRO_34 * DENDRO_C2_k2_0 + DENDRO_40 * DENDRO_C2_k2_3 +
        DENDRO_44 * DENDRO_C2_k2_5 + DENDRO_46 * DENDRO_C2_k2_1 +
        DENDRO_47 * DENDRO_C2_k2_2 + DENDRO_49 * DENDRO_C2_k2_4 +
        DENDRO_53 * deriv_evars->grad2_0_1_beta0[d_pp] +
        DENDRO_53 * deriv_evars->grad2_1_1_beta1[d_pp] +
        DENDRO_54 * deriv_evars->grad2_1_2_beta2[d_pp] +
        DENDRO_57 * deriv_evars->grad2_0_2_beta0[d_pp] +
        DENDRO_57 * deriv_evars->grad2_1_2_beta1[d_pp] + DENDRO_58 +
        DENDRO_igt0 * deriv_evars->grad2_0_0_beta2[d_pp] +
        DENDRO_igt3 * deriv_evars->grad2_1_1_beta2[d_pp] +
        (4.0 / 3.0) * DENDRO_igt5 * deriv_evars->grad2_2_2_beta2[d_pp];

    // Dendro: printing variables
    //--
    Gt_rhs0[pp] = DENDRO_51;
    //--
    Gt_rhs1[pp] = DENDRO_56;
    //--
    Gt_rhs2[pp] = DENDRO_59;
    //--
    B_rhs0[pp]  = -B0[pp] * eta - DENDRO_50 * lambda[3] + DENDRO_51 +
                 lambda[2] * ((deriv_evars->grad_0_B0[d_pp]) * beta0[pp] +
                              (deriv_evars->grad_1_B0[d_pp]) * beta1[pp] +
                              (deriv_evars->grad_2_B0[d_pp]) * beta2[pp]);
    //--
    B_rhs1[pp] = -B1[pp] * eta - DENDRO_55 * lambda[3] + DENDRO_56 +
                 lambda[2] * ((deriv_evars->grad_0_B1[d_pp]) * beta0[pp] +
                              (deriv_evars->grad_1_B1[d_pp]) * beta1[pp] +
                              (deriv_evars->grad_2_B1[d_pp]) * beta2[pp]);
    //--
    B_rhs2[pp] = -B2[pp] * eta - DENDRO_58 * lambda[3] + DENDRO_59 +
                 lambda[2] * ((deriv_evars->grad_0_B2[d_pp]) * beta0[pp] +
                              (deriv_evars->grad_1_B2[d_pp]) * beta1[pp] +
                              (deriv_evars->grad_2_B2[d_pp]) * beta2[pp]);
    // Dendro: reduced ops: 400
    // Dendro: }}}
}
