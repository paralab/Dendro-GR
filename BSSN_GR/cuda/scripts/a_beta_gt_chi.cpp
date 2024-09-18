
{
    // Dendro: {{{
    // Dendro: original ops: 223
    // Dendro: printing temp variables
    const double DENDRO_0 = 2 * alpha[pp];
    const double DENDRO_1 =
        (3.0 / 4.0) * alpha[pp] * lambda_f[1] + (3.0 / 4.0) * lambda_f[0];
    const double DENDRO_2  = 2 * gt1[pp];
    const double DENDRO_3  = 2 * gt2[pp];
    const double DENDRO_4  = (2.0 / 3.0) * gt0[pp];
    const double DENDRO_5  = (1.0 / 3.0) * gt1[pp];
    const double DENDRO_6  = (2.0 / 3.0) * (deriv_evars->grad_2_beta2[d_pp]);
    const double DENDRO_7  = (1.0 / 3.0) * gt2[pp];
    const double DENDRO_8  = (2.0 / 3.0) * (deriv_evars->grad_1_beta1[d_pp]);
    const double DENDRO_9  = (2.0 / 3.0) * (deriv_evars->grad_0_beta0[d_pp]);
    const double DENDRO_10 = 2 * gt4[pp];
    const double DENDRO_11 = (1.0 / 3.0) * gt4[pp];
    const double DENDRO_12 = (2.0 / 3.0) * chi[pp];

    // Dendro: printing variables
    //--
    a_rhs[pp]              = -DENDRO_0 * K[pp] +
                lambda[0] * ((deriv_evars->grad_0_alpha[d_pp]) * beta0[pp] +
                             (deriv_evars->grad_1_alpha[d_pp]) * beta1[pp] +
                             (deriv_evars->grad_2_alpha[d_pp]) * beta2[pp]);
    //--
    b_rhs0[pp] = B0[pp] * DENDRO_1 +
                 lambda[1] * ((deriv_evars->grad_0_beta0[d_pp]) * beta0[pp] +
                              (deriv_evars->grad_1_beta0[d_pp]) * beta1[pp] +
                              (deriv_evars->grad_2_beta0[d_pp]) * beta2[pp]);
    //--
    b_rhs1[pp] = B1[pp] * DENDRO_1 +
                 lambda[1] * ((deriv_evars->grad_0_beta1[d_pp]) * beta0[pp] +
                              (deriv_evars->grad_1_beta1[d_pp]) * beta1[pp] +
                              (deriv_evars->grad_2_beta1[d_pp]) * beta2[pp]);
    //--
    b_rhs2[pp] = B2[pp] * DENDRO_1 +
                 lambda[1] * ((deriv_evars->grad_0_beta2[d_pp]) * beta0[pp] +
                              (deriv_evars->grad_1_beta2[d_pp]) * beta1[pp] +
                              (deriv_evars->grad_2_beta2[d_pp]) * beta2[pp]);
    //--
    gt_rhs00[pp] = (4.0 / 3.0) * (deriv_evars->grad_0_beta0[d_pp]) * gt0[pp] +
                   (deriv_evars->grad_0_beta1[d_pp]) * DENDRO_2 +
                   (deriv_evars->grad_0_beta2[d_pp]) * DENDRO_3 +
                   (deriv_evars->grad_0_gt0[d_pp]) * beta0[pp] -
                   (deriv_evars->grad_1_beta1[d_pp]) * DENDRO_4 +
                   (deriv_evars->grad_1_gt0[d_pp]) * beta1[pp] -
                   (deriv_evars->grad_2_beta2[d_pp]) * DENDRO_4 +
                   (deriv_evars->grad_2_gt0[d_pp]) * beta2[pp] -
                   At0[pp] * DENDRO_0;
    //--
    gt_rhs01[pp] = (deriv_evars->grad_0_beta0[d_pp]) * DENDRO_5 +
                   (deriv_evars->grad_0_beta1[d_pp]) * gt3[pp] +
                   (deriv_evars->grad_0_beta2[d_pp]) * gt4[pp] +
                   (deriv_evars->grad_0_gt1[d_pp]) * beta0[pp] +
                   (deriv_evars->grad_1_beta0[d_pp]) * gt0[pp] +
                   (deriv_evars->grad_1_beta1[d_pp]) * DENDRO_5 +
                   (deriv_evars->grad_1_beta2[d_pp]) * gt2[pp] +
                   (deriv_evars->grad_1_gt1[d_pp]) * beta1[pp] +
                   (deriv_evars->grad_2_gt1[d_pp]) * beta2[pp] -
                   At1[pp] * DENDRO_0 - DENDRO_6 * gt1[pp];
    //--
    gt_rhs02[pp] = (deriv_evars->grad_0_beta0[d_pp]) * DENDRO_7 +
                   (deriv_evars->grad_0_beta1[d_pp]) * gt4[pp] +
                   (deriv_evars->grad_0_beta2[d_pp]) * gt5[pp] +
                   (deriv_evars->grad_0_gt2[d_pp]) * beta0[pp] +
                   (deriv_evars->grad_1_gt2[d_pp]) * beta1[pp] +
                   (deriv_evars->grad_2_beta0[d_pp]) * gt0[pp] +
                   (deriv_evars->grad_2_beta1[d_pp]) * gt1[pp] +
                   (deriv_evars->grad_2_beta2[d_pp]) * DENDRO_7 +
                   (deriv_evars->grad_2_gt2[d_pp]) * beta2[pp] -
                   At2[pp] * DENDRO_0 - DENDRO_8 * gt2[pp];
    //--
    gt_rhs11[pp] = (deriv_evars->grad_0_gt3[d_pp]) * beta0[pp] +
                   (deriv_evars->grad_1_beta0[d_pp]) * DENDRO_2 +
                   (4.0 / 3.0) * (deriv_evars->grad_1_beta1[d_pp]) * gt3[pp] +
                   (deriv_evars->grad_1_beta2[d_pp]) * DENDRO_10 +
                   (deriv_evars->grad_1_gt3[d_pp]) * beta1[pp] +
                   (deriv_evars->grad_2_gt3[d_pp]) * beta2[pp] -
                   At3[pp] * DENDRO_0 - DENDRO_6 * gt3[pp] - DENDRO_9 * gt3[pp];
    //--
    gt_rhs12[pp] = (deriv_evars->grad_0_gt4[d_pp]) * beta0[pp] +
                   (deriv_evars->grad_1_beta0[d_pp]) * gt2[pp] +
                   (deriv_evars->grad_1_beta1[d_pp]) * DENDRO_11 +
                   (deriv_evars->grad_1_beta2[d_pp]) * gt5[pp] +
                   (deriv_evars->grad_1_gt4[d_pp]) * beta1[pp] +
                   (deriv_evars->grad_2_beta0[d_pp]) * gt1[pp] +
                   (deriv_evars->grad_2_beta1[d_pp]) * gt3[pp] +
                   (deriv_evars->grad_2_beta2[d_pp]) * DENDRO_11 +
                   (deriv_evars->grad_2_gt4[d_pp]) * beta2[pp] -
                   At4[pp] * DENDRO_0 - DENDRO_9 * gt4[pp];
    //--
    gt_rhs22[pp] = (deriv_evars->grad_0_gt5[d_pp]) * beta0[pp] +
                   (deriv_evars->grad_1_gt5[d_pp]) * beta1[pp] +
                   (deriv_evars->grad_2_beta0[d_pp]) * DENDRO_3 +
                   (deriv_evars->grad_2_beta1[d_pp]) * DENDRO_10 +
                   (4.0 / 3.0) * (deriv_evars->grad_2_beta2[d_pp]) * gt5[pp] +
                   (deriv_evars->grad_2_gt5[d_pp]) * beta2[pp] -
                   At5[pp] * DENDRO_0 - DENDRO_8 * gt5[pp] - DENDRO_9 * gt5[pp];
    //--
    chi_rhs[pp] = (deriv_evars->grad_0_chi[d_pp]) * beta0[pp] +
                  (deriv_evars->grad_1_chi[d_pp]) * beta1[pp] +
                  (deriv_evars->grad_2_chi[d_pp]) * beta2[pp] +
                  DENDRO_12 * K[pp] * alpha[pp] -
                  DENDRO_12 * ((deriv_evars->grad_0_beta0[d_pp]) +
                               (deriv_evars->grad_1_beta1[d_pp]) +
                               (deriv_evars->grad_2_beta2[d_pp]));
    // Dendro: reduced ops: 187
    // Dendro: }}}
}
