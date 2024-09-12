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
            // Dendro: original ops:  1914
            // Dendro: printing temp variables
            double DENDRO_0 = pow(gt4[pp], 2);
            double DENDRO_1 = pow(gt1[pp], 2);
            double DENDRO_2 = pow(gt2[pp], 2);
            double DENDRO_3 = gt0[pp] * gt3[pp];
            double DENDRO_4 = gt1[pp] * gt2[pp];
            double DENDRO_5 =
                2 / pow(DENDRO_0 * gt0[pp] + DENDRO_1 * gt5[pp] +
                            DENDRO_2 * gt3[pp] - DENDRO_3 * gt5[pp] -
                            2 * DENDRO_4 * gt4[pp],
                        2);
            double DENDRO_6  = grad_0_alpha[pp];
            double DENDRO_7  = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
            double DENDRO_8  = pow(DENDRO_7, 2);
            double DENDRO_9  = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
            double DENDRO_10 = pow(DENDRO_9, 2);
            double DENDRO_11 = -DENDRO_0 + gt3[pp] * gt5[pp];
            double DENDRO_12 = DENDRO_7 * DENDRO_9;
            double DENDRO_13 = 2 * At1[pp] * DENDRO_7;
            double DENDRO_14 = 2 * At2[pp] * DENDRO_9;
            double DENDRO_15 = grad_2_alpha[pp];
            double DENDRO_16 = -DENDRO_4 + gt0[pp] * gt4[pp];
            double DENDRO_17 = DENDRO_16 * DENDRO_7;
            double DENDRO_18 = At0[pp] * DENDRO_11;
            double DENDRO_19 = DENDRO_16 * DENDRO_9;
            double DENDRO_20 = -DENDRO_1 + DENDRO_3;
            double DENDRO_21 = At5[pp] * DENDRO_20;
            double DENDRO_22 = DENDRO_11 * DENDRO_16;
            double DENDRO_23 = DENDRO_20 * DENDRO_7;
            double DENDRO_24 = -At1[pp] * DENDRO_12 - At1[pp] * DENDRO_22 +
                               At2[pp] * DENDRO_10 +
                               At2[pp] * DENDRO_11 * DENDRO_20 +
                               At3[pp] * DENDRO_17 - At4[pp] * DENDRO_19 -
                               At4[pp] * DENDRO_23 + DENDRO_18 * DENDRO_9 +
                               DENDRO_21 * DENDRO_9;
            double DENDRO_25 = grad_1_alpha[pp];
            double DENDRO_26 = -DENDRO_2 + gt0[pp] * gt5[pp];
            double DENDRO_27 = At3[pp] * DENDRO_26;
            double DENDRO_28 = DENDRO_26 * DENDRO_9;
            double DENDRO_29 = -At1[pp] * DENDRO_11 * DENDRO_26 -
                               At1[pp] * DENDRO_8 + At2[pp] * DENDRO_12 +
                               At2[pp] * DENDRO_22 - At4[pp] * DENDRO_17 -
                               At4[pp] * DENDRO_28 + At5[pp] * DENDRO_19 +
                               DENDRO_18 * DENDRO_7 + DENDRO_27 * DENDRO_7;
            double DENDRO_30 = pow(DENDRO_16, 2);
            double DENDRO_31 = 2 * At4[pp] * DENDRO_16;
            double DENDRO_32 = At0[pp] * DENDRO_12 - At1[pp] * DENDRO_17 -
                               At1[pp] * DENDRO_28 + At2[pp] * DENDRO_19 +
                               At2[pp] * DENDRO_23 -
                               At4[pp] * DENDRO_20 * DENDRO_26 -
                               At4[pp] * DENDRO_30 + DENDRO_16 * DENDRO_21 +
                               DENDRO_16 * DENDRO_27;
            // Dendro: printing variables
            //--
            Gt_rhs_s5_0[pp] =
                DENDRO_5 *
                (DENDRO_15 * DENDRO_24 - DENDRO_25 * DENDRO_29 +
                 DENDRO_6 * (At0[pp] * pow(DENDRO_11, 2) + At3[pp] * DENDRO_8 -
                             2 * At4[pp] * DENDRO_12 + At5[pp] * DENDRO_10 -
                             DENDRO_11 * DENDRO_13 + DENDRO_11 * DENDRO_14));
            //--
            Gt_rhs_s5_1[pp] =
                DENDRO_5 *
                (-DENDRO_15 * DENDRO_32 +
                 DENDRO_25 *
                     (At0[pp] * DENDRO_8 + 2 * At2[pp] * DENDRO_17 +
                      At3[pp] * pow(DENDRO_26, 2) + At5[pp] * DENDRO_30 -
                      DENDRO_13 * DENDRO_26 - DENDRO_26 * DENDRO_31) -
                 DENDRO_29 * DENDRO_6);
            //--
            Gt_rhs_s5_2[pp] =
                DENDRO_5 *
                (DENDRO_15 *
                     (At0[pp] * DENDRO_10 - 2 * At1[pp] * DENDRO_19 +
                      At3[pp] * DENDRO_30 + At5[pp] * pow(DENDRO_20, 2) +
                      DENDRO_14 * DENDRO_20 - DENDRO_20 * DENDRO_31) +
                 DENDRO_24 * DENDRO_6 - DENDRO_25 * DENDRO_32);
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
