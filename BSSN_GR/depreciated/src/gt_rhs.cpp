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
            // Dendro: original ops:  210
            // Dendro: printing temp variables
            double DENDRO_0  = 2 * alpha[pp];
            double DENDRO_1  = grad_0_beta0[pp];
            double DENDRO_2  = (2.0L / 3.0L) * gt0[pp];
            double DENDRO_3  = grad_1_beta1[pp];
            double DENDRO_4  = grad_2_beta2[pp];
            double DENDRO_5  = 2 * gt1[pp];
            double DENDRO_6  = grad_0_beta1[pp];
            double DENDRO_7  = 2 * gt2[pp];
            double DENDRO_8  = grad_0_beta2[pp];
            double DENDRO_9  = grad_1_beta0[pp];
            double DENDRO_10 = grad_1_beta2[pp];
            double DENDRO_11 = (1.0L / 3.0L) * gt1[pp];
            double DENDRO_12 = (2.0L / 3.0L) * DENDRO_4;
            double DENDRO_13 = grad_2_beta0[pp];
            double DENDRO_14 = grad_2_beta1[pp];
            double DENDRO_15 = (1.0L / 3.0L) * gt2[pp];
            double DENDRO_16 = (2.0L / 3.0L) * DENDRO_3;
            double DENDRO_17 = (2.0L / 3.0L) * DENDRO_1;
            double DENDRO_18 = 2 * gt4[pp];
            double DENDRO_19 = (1.0L / 3.0L) * gt4[pp];
            // Dendro: printing variables
            //--
            gt_rhs00[pp] =
                -At0[pp] * DENDRO_0 + (4.0L / 3.0L) * DENDRO_1 * gt0[pp] -
                DENDRO_2 * DENDRO_3 - DENDRO_2 * DENDRO_4 +
                DENDRO_5 * DENDRO_6 + DENDRO_7 * DENDRO_8 +
                beta0[pp] * agrad_0_gt0[pp] + beta1[pp] * agrad_1_gt0[pp] +
                beta2[pp] * agrad_2_gt0[pp];
            //--
            gt_rhs01[pp] =
                -At1[pp] * DENDRO_0 + DENDRO_1 * DENDRO_11 +
                DENDRO_10 * gt2[pp] + DENDRO_11 * DENDRO_3 -
                DENDRO_12 * gt1[pp] + DENDRO_6 * gt3[pp] + DENDRO_8 * gt4[pp] +
                DENDRO_9 * gt0[pp] + beta0[pp] * agrad_0_gt1[pp] +
                beta1[pp] * agrad_1_gt1[pp] + beta2[pp] * agrad_2_gt1[pp];
            //--
            gt_rhs02[pp] = -At2[pp] * DENDRO_0 + DENDRO_1 * DENDRO_15 +
                           DENDRO_13 * gt0[pp] + DENDRO_14 * gt1[pp] +
                           DENDRO_15 * DENDRO_4 - DENDRO_16 * gt2[pp] +
                           DENDRO_6 * gt4[pp] + DENDRO_8 * gt5[pp] +
                           beta0[pp] * agrad_0_gt2[pp] +
                           beta1[pp] * agrad_1_gt2[pp] +
                           beta2[pp] * agrad_2_gt2[pp];
            //--
            gt_rhs11[pp] = -At3[pp] * DENDRO_0 + DENDRO_10 * DENDRO_18 -
                           DENDRO_12 * gt3[pp] - DENDRO_17 * gt3[pp] +
                           (4.0L / 3.0L) * DENDRO_3 * gt3[pp] +
                           DENDRO_5 * DENDRO_9 + beta0[pp] * agrad_0_gt3[pp] +
                           beta1[pp] * agrad_1_gt3[pp] +
                           beta2[pp] * agrad_2_gt3[pp];
            //--
            gt_rhs12[pp] = -At4[pp] * DENDRO_0 + DENDRO_10 * gt5[pp] +
                           DENDRO_13 * gt1[pp] + DENDRO_14 * gt3[pp] -
                           DENDRO_17 * gt4[pp] + DENDRO_19 * DENDRO_3 +
                           DENDRO_19 * DENDRO_4 + DENDRO_9 * gt2[pp] +
                           beta0[pp] * agrad_0_gt4[pp] +
                           beta1[pp] * agrad_1_gt4[pp] +
                           beta2[pp] * agrad_2_gt4[pp];
            //--
            gt_rhs22[pp] =
                -At5[pp] * DENDRO_0 + DENDRO_13 * DENDRO_7 +
                DENDRO_14 * DENDRO_18 - DENDRO_16 * gt5[pp] -
                DENDRO_17 * gt5[pp] + (4.0L / 3.0L) * DENDRO_4 * gt5[pp] +
                beta0[pp] * agrad_0_gt5[pp] + beta1[pp] * agrad_1_gt5[pp] +
                beta2[pp] * agrad_2_gt5[pp];
            // Dendro: reduced ops:  162
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
