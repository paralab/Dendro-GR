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
            // Dendro: original ops:  1716
            // Dendro: printing temp variables
            double DENDRO_0  = pow(gt4[pp], 2);
            double DENDRO_1  = pow(gt1[pp], 2);
            double DENDRO_2  = pow(gt2[pp], 2);
            double DENDRO_3  = gt0[pp] * gt3[pp];
            double DENDRO_4  = gt1[pp] * gt2[pp];
            double DENDRO_5  = pow(DENDRO_0 * gt0[pp] + DENDRO_1 * gt5[pp] +
                                       DENDRO_2 * gt3[pp] - DENDRO_3 * gt5[pp] -
                                       2 * DENDRO_4 * gt4[pp],
                                   -2);
            double DENDRO_6  = -DENDRO_4 + gt0[pp] * gt4[pp];
            double DENDRO_7  = grad_2_gt3[pp];
            double DENDRO_8  = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
            double DENDRO_9  = grad_1_gt5[pp];
            double DENDRO_10 = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
            double DENDRO_11 = -DENDRO_0 + gt3[pp] * gt5[pp];
            double DENDRO_12 = grad_1_gt2[pp];
            double DENDRO_13 = grad_2_gt1[pp];
            double DENDRO_14 = grad_0_gt4[pp];
            double DENDRO_15 = DENDRO_12 + DENDRO_13 - DENDRO_14;
            double DENDRO_16 = grad_0_gt3[pp];
            double DENDRO_17 = grad_1_gt0[pp];
            double DENDRO_18 = DENDRO_12 - DENDRO_13 + DENDRO_14;
            double DENDRO_19 =
                1.0 * gt1[pp] * gt4[pp] - 1.0 * gt2[pp] * gt3[pp];
            double DENDRO_20 = grad_0_gt5[pp];
            double DENDRO_21 = grad_2_gt0[pp];
            double DENDRO_22 = -DENDRO_12 + DENDRO_13 + DENDRO_14;
            double DENDRO_23 = -DENDRO_1 + DENDRO_3;
            double DENDRO_24 = grad_2_gt5[pp];
            double DENDRO_25 =
                0.5 * gt1[pp] * gt4[pp] - 0.5 * gt2[pp] * gt3[pp];
            double DENDRO_26 = 0.5 * DENDRO_9 - 1.0 * grad_2_gt4[pp];
            double DENDRO_27 = 0.5 * DENDRO_20 - 1.0 * grad_2_gt2[pp];
            double DENDRO_28 = -DENDRO_2 + gt0[pp] * gt5[pp];
            double DENDRO_29 = grad_1_gt3[pp];
            double DENDRO_30 =
                0.5 * gt1[pp] * gt5[pp] - 0.5 * gt2[pp] * gt4[pp];
            double DENDRO_31 = -0.5 * DENDRO_7 + 1.0 * grad_1_gt4[pp];
            double DENDRO_32 = 0.5 * DENDRO_16 - 1.0 * grad_1_gt1[pp];
            double DENDRO_33 = grad_0_gt0[pp];
            double DENDRO_34 = -0.5 * DENDRO_17 + 1.0 * grad_0_gt1[pp];
            double DENDRO_35 = -0.5 * DENDRO_21 + 1.0 * grad_0_gt2[pp];
            double DENDRO_36 =
                0.5 * gt0[pp] * gt4[pp] - 0.5 * gt1[pp] * gt2[pp];
            // Dendro: printing variables
            //--
            CalGt0[pp] =
                DENDRO_5 *
                (-DENDRO_11 *
                     (-DENDRO_10 * DENDRO_35 - 0.5 * DENDRO_11 * DENDRO_33 +
                      DENDRO_34 * DENDRO_8) -
                 DENDRO_19 * (-DENDRO_10 * DENDRO_20 - DENDRO_11 * DENDRO_21 +
                              DENDRO_22 * DENDRO_8) -
                 DENDRO_23 * (DENDRO_11 * DENDRO_27 - DENDRO_24 * DENDRO_25 -
                              DENDRO_26 * DENDRO_8) -
                 DENDRO_28 * (-DENDRO_10 * DENDRO_31 + DENDRO_11 * DENDRO_32 +
                              DENDRO_29 * DENDRO_30) +
                 DENDRO_6 * (-DENDRO_10 * DENDRO_9 - DENDRO_11 * DENDRO_15 +
                             DENDRO_7 * DENDRO_8) +
                 DENDRO_8 * (-DENDRO_10 * DENDRO_18 - DENDRO_11 * DENDRO_17 +
                             DENDRO_16 * DENDRO_8));
            //--
            CalGt1[pp] =
                DENDRO_5 *
                (-DENDRO_11 * (-DENDRO_28 * DENDRO_34 + DENDRO_30 * DENDRO_33 +
                               DENDRO_35 * DENDRO_6) -
                 DENDRO_19 * (DENDRO_20 * DENDRO_6 + DENDRO_21 * DENDRO_8 -
                              DENDRO_22 * DENDRO_28) -
                 DENDRO_23 * (DENDRO_24 * DENDRO_36 + DENDRO_26 * DENDRO_28 -
                              DENDRO_27 * DENDRO_8) -
                 DENDRO_28 * (-0.5 * DENDRO_28 * DENDRO_29 +
                              DENDRO_31 * DENDRO_6 - DENDRO_32 * DENDRO_8) +
                 DENDRO_6 * (DENDRO_15 * DENDRO_8 - DENDRO_28 * DENDRO_7 +
                             DENDRO_6 * DENDRO_9) +
                 DENDRO_8 * (-DENDRO_16 * DENDRO_28 + DENDRO_17 * DENDRO_8 +
                             DENDRO_18 * DENDRO_6));
            //--
            CalGt2[pp] =
                DENDRO_5 *
                (-DENDRO_11 * (-DENDRO_23 * DENDRO_35 - DENDRO_25 * DENDRO_33 +
                               DENDRO_34 * DENDRO_6) -
                 DENDRO_19 * (-DENDRO_10 * DENDRO_21 - DENDRO_20 * DENDRO_23 +
                              DENDRO_22 * DENDRO_6) -
                 DENDRO_23 *
                     (DENDRO_10 * DENDRO_27 - 0.5 * DENDRO_23 * DENDRO_24 -
                      DENDRO_26 * DENDRO_6) -
                 DENDRO_28 * (DENDRO_10 * DENDRO_32 - DENDRO_23 * DENDRO_31 +
                              DENDRO_29 * DENDRO_36) +
                 DENDRO_6 * (-DENDRO_10 * DENDRO_15 - DENDRO_23 * DENDRO_9 +
                             DENDRO_6 * DENDRO_7) +
                 DENDRO_8 * (-DENDRO_10 * DENDRO_17 + DENDRO_16 * DENDRO_6 -
                             DENDRO_18 * DENDRO_23));
            // Dendro: reduced ops:  221
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
