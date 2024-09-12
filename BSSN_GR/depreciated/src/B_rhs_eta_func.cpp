bssn::timer::t_rhs.start();
for (unsigned int k = PW; k < nz - PW; k++) {
    z = pmin[2] + k * hz;
    for (unsigned int j = PW; j < ny - PW; j++) {
        y = pmin[1] + j * hy;
        for (unsigned int i = PW; i < nx - PW; i++) {
            x       = pmin[0] + i * hx;
            pp      = i + nx * (j + ny * k);
            r_coord = sqrt(x * x + y * y + z * z);
            eta     = ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
            }
            // Dendro: {{{
            // Dendro: original ops:  516
            // Dendro: printing temp variables
            double DENDRO_0 = pow(gt4[pp], 2);
            double DENDRO_1 = pow(gt1[pp], 2);
            double DENDRO_2 = pow(gt2[pp], 2);
            double DENDRO_3 = gt3[pp] * gt5[pp];
            double DENDRO_4 = gt1[pp] * gt4[pp];
            double DENDRO_5 = grad_2_chi[pp];
            double DENDRO_6 = grad_1_chi[pp];
            double DENDRO_7 = grad_0_chi[pp];
            double DENDRO_8 = 2 * DENDRO_5;
            double DENDRO_9 =
                BSSN_ETA_R0 *
                sqrt((-pow(DENDRO_5, 2) * (-DENDRO_1 + gt0[pp] * gt3[pp]) -
                      pow(DENDRO_6, 2) * (-DENDRO_2 + gt0[pp] * gt5[pp]) +
                      2 * DENDRO_6 * DENDRO_7 *
                          (gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp]) +
                      DENDRO_6 * DENDRO_8 *
                          (gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp]) -
                      pow(DENDRO_7, 2) * (-DENDRO_0 + DENDRO_3) -
                      DENDRO_7 * DENDRO_8 * (DENDRO_4 - gt2[pp] * gt3[pp])) /
                     (DENDRO_0 * gt0[pp] + DENDRO_1 * gt5[pp] +
                      DENDRO_2 * gt3[pp] - DENDRO_3 * gt0[pp] -
                      2 * DENDRO_4 * gt2[pp])) *
                pow(-pow(chi[pp], BSSN_ETA_POWER[0]) + 1, -BSSN_ETA_POWER[1]);
            // Dendro: printing variables
            //--
            B_rhs0[pp] = -B0[pp] * DENDRO_9 + Gt_rhs0[pp] +
                         lambda[2] * (beta0[pp] * agrad_0_B0[pp] +
                                      beta1[pp] * agrad_1_B0[pp] +
                                      beta2[pp] * agrad_2_B0[pp]) -
                         lambda[3] * (beta0[pp] * agrad_0_Gt0[pp] +
                                      beta1[pp] * agrad_1_Gt0[pp] +
                                      beta2[pp] * agrad_2_Gt0[pp]);
            //--
            B_rhs1[pp] = -B1[pp] * DENDRO_9 + Gt_rhs1[pp] +
                         lambda[2] * (beta0[pp] * agrad_0_B1[pp] +
                                      beta1[pp] * agrad_1_B1[pp] +
                                      beta2[pp] * agrad_2_B1[pp]) -
                         lambda[3] * (beta0[pp] * agrad_0_Gt1[pp] +
                                      beta1[pp] * agrad_1_Gt1[pp] +
                                      beta2[pp] * agrad_2_Gt1[pp]);
            //--
            B_rhs2[pp] = -B2[pp] * DENDRO_9 + Gt_rhs2[pp] +
                         lambda[2] * (beta0[pp] * agrad_0_B2[pp] +
                                      beta1[pp] * agrad_1_B2[pp] +
                                      beta2[pp] * agrad_2_B2[pp]) -
                         lambda[3] * (beta0[pp] * agrad_0_Gt2[pp] +
                                      beta1[pp] * agrad_1_Gt2[pp] +
                                      beta2[pp] * agrad_2_Gt2[pp]);
            // Dendro: reduced ops:  124
            // Dendro: }}}
            /* debugging */
            /*unsigned int qi = 46 - 1;
            unsigned int qj = 10 - 1;
            unsigned int qk = 60 - 1;
            unsigned int qidx = qi + nx*(qj + ny*qk);
            if (0 && qidx == pp) {
            std::cout << ".... end OPTIMIZED debug stuff..." << std::endl;
            }*/
        }
    }
}
bssn::timer::t_rhs.stop();
