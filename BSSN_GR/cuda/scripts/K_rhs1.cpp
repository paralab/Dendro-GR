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
double AijAuij = pow(At0[pp], 2) * (DENDRO_igt0 * DENDRO_igt0) +
                 4 * At0[pp] * At1[pp] * DENDRO_igt0 * DENDRO_igt1 +
                 4 * At0[pp] * At2[pp] * DENDRO_igt0 * DENDRO_igt2 +
                 2 * At0[pp] * At3[pp] * (DENDRO_igt1 * DENDRO_igt1) +
                 4 * At0[pp] * At4[pp] * DENDRO_igt1 * DENDRO_igt2 +
                 2 * At0[pp] * At5[pp] * (DENDRO_igt2 * DENDRO_igt2) +
                 2 * pow(At1[pp], 2) * DENDRO_igt0 * DENDRO_igt3 +
                 2 * pow(At1[pp], 2) * (DENDRO_igt1 * DENDRO_igt1) +
                 4 * At1[pp] * At2[pp] * DENDRO_igt0 * DENDRO_igt4 +
                 4 * At1[pp] * At2[pp] * DENDRO_igt1 * DENDRO_igt2 +
                 4 * At1[pp] * At3[pp] * DENDRO_igt1 * DENDRO_igt3 +
                 4 * At1[pp] * At4[pp] * DENDRO_igt1 * DENDRO_igt4 +
                 4 * At1[pp] * At4[pp] * DENDRO_igt2 * DENDRO_igt3 +
                 4 * At1[pp] * At5[pp] * DENDRO_igt2 * DENDRO_igt4 +
                 2 * pow(At2[pp], 2) * DENDRO_igt0 * DENDRO_igt5 +
                 2 * pow(At2[pp], 2) * (DENDRO_igt2 * DENDRO_igt2) +
                 4 * At2[pp] * At3[pp] * DENDRO_igt1 * DENDRO_igt4 +
                 4 * At2[pp] * At4[pp] * DENDRO_igt1 * DENDRO_igt5 +
                 4 * At2[pp] * At4[pp] * DENDRO_igt2 * DENDRO_igt4 +
                 4 * At2[pp] * At5[pp] * DENDRO_igt2 * DENDRO_igt5 +
                 pow(At3[pp], 2) * (DENDRO_igt3 * DENDRO_igt3) +
                 4 * At3[pp] * At4[pp] * DENDRO_igt3 * DENDRO_igt4 +
                 2 * At3[pp] * At5[pp] * (DENDRO_igt4 * DENDRO_igt4) +
                 2 * pow(At4[pp], 2) * DENDRO_igt3 * DENDRO_igt5 +
                 2 * pow(At4[pp], 2) * (DENDRO_igt4 * DENDRO_igt4) +
                 4 * At4[pp] * At5[pp] * DENDRO_igt4 * DENDRO_igt5 +
                 pow(At5[pp], 2) * (DENDRO_igt5 * DENDRO_igt5);

{
    // Dendro: {{{
    // Dendro: original ops: 10
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    K_rhs[pp] = (deriv_evars->grad_0_K[d_pp]) * beta0[pp] +
                (deriv_evars->grad_1_K[d_pp]) * beta1[pp] +
                (deriv_evars->grad_2_K[d_pp]) * beta2[pp] +
                (1.0 / 3.0) * alpha[pp] * (3 * AijAuij + pow(K[pp], 2));
    // Dendro: reduced ops: 11
    // Dendro: }}}
}
