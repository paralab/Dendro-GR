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
            // Dendro: original ops:  4812
            // Dendro: printing temp variables
            double DENDRO_0 = pow(gt4[pp], 2);
            double DENDRO_1 = pow(gt1[pp], 2);
            double DENDRO_2 = pow(gt2[pp], 2);
            double DENDRO_3 = gt0[pp] * gt3[pp];
            double DENDRO_4 = gt1[pp] * gt2[pp];
            double DENDRO_5 = 2 * alpha[pp] /
                              pow(DENDRO_0 * gt0[pp] + DENDRO_1 * gt5[pp] +
                                      DENDRO_2 * gt3[pp] - DENDRO_3 * gt5[pp] -
                                      2 * DENDRO_4 * gt4[pp],
                                  3);
            double DENDRO_6 = grad_1_gt3[pp];
            double DENDRO_7 = 0.5 * gt1[pp] * gt5[pp] - 0.5 * gt2[pp] * gt4[pp];
            double DENDRO_8 = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
            double DENDRO_9 = grad_2_gt3[pp];
            double DENDRO_10 = -0.5 * DENDRO_9 + 1.0 * grad_1_gt4[pp];
            double DENDRO_11 = -DENDRO_0 + gt3[pp] * gt5[pp];
            double DENDRO_12 = grad_0_gt3[pp];
            double DENDRO_13 = 0.5 * DENDRO_12 - 1.0 * grad_1_gt1[pp];
            double DENDRO_14 = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
            double DENDRO_15 = pow(DENDRO_14, 2);
            double DENDRO_16 = -DENDRO_4 + gt0[pp] * gt4[pp];
            double DENDRO_17 = pow(DENDRO_16, 2);
            double DENDRO_18 = -DENDRO_2 + gt0[pp] * gt5[pp];
            double DENDRO_19 = DENDRO_14 * DENDRO_16;
            double DENDRO_20 = 2 * At1[pp] * DENDRO_14;
            double DENDRO_21 = 2 * At4[pp] * DENDRO_16;
            double DENDRO_22 = At0[pp] * DENDRO_15 + 2 * At2[pp] * DENDRO_19 +
                               At3[pp] * pow(DENDRO_18, 2) +
                               At5[pp] * DENDRO_17 - DENDRO_18 * DENDRO_20 -
                               DENDRO_18 * DENDRO_21;
            double DENDRO_23 = grad_0_gt0[pp];
            double DENDRO_24 = grad_2_gt0[pp];
            double DENDRO_25 = -0.5 * DENDRO_24 + 1.0 * grad_0_gt2[pp];
            double DENDRO_26 = grad_1_gt0[pp];
            double DENDRO_27 = -0.5 * DENDRO_26 + 1.0 * grad_0_gt1[pp];
            double DENDRO_28 = pow(DENDRO_8, 2);
            double DENDRO_29 = DENDRO_14 * DENDRO_8;
            double DENDRO_30 = 2 * At2[pp] * DENDRO_8;
            double DENDRO_31 = At0[pp] * pow(DENDRO_11, 2) +
                               At3[pp] * DENDRO_15 - 2 * At4[pp] * DENDRO_29 +
                               At5[pp] * DENDRO_28 - DENDRO_11 * DENDRO_20 +
                               DENDRO_11 * DENDRO_30;
            double DENDRO_32 = grad_2_gt5[pp];
            double DENDRO_33 =
                0.5 * gt1[pp] * gt4[pp] - 0.5 * gt2[pp] * gt3[pp];
            double DENDRO_34 = grad_1_gt5[pp];
            double DENDRO_35 = 0.5 * DENDRO_34 - 1.0 * grad_2_gt4[pp];
            double DENDRO_36 = grad_0_gt5[pp];
            double DENDRO_37 = 0.5 * DENDRO_36 - 1.0 * grad_2_gt2[pp];
            double DENDRO_38 = -DENDRO_1 + DENDRO_3;
            double DENDRO_39 = DENDRO_16 * DENDRO_8;
            double DENDRO_40 = At0[pp] * DENDRO_28 - 2 * At1[pp] * DENDRO_39 +
                               At3[pp] * DENDRO_17 +
                               At5[pp] * pow(DENDRO_38, 2) -
                               DENDRO_21 * DENDRO_38 + DENDRO_30 * DENDRO_38;
            double DENDRO_41 = grad_0_gt4[pp];
            double DENDRO_42 = grad_1_gt2[pp];
            double DENDRO_43 = grad_2_gt1[pp];
            double DENDRO_44 = DENDRO_41 + DENDRO_42 - DENDRO_43;
            double DENDRO_45 = At0[pp] * DENDRO_11;
            double DENDRO_46 = DENDRO_11 * DENDRO_16;
            double DENDRO_47 = At3[pp] * DENDRO_18;
            double DENDRO_48 = DENDRO_18 * DENDRO_8;
            double DENDRO_49 = -At1[pp] * DENDRO_11 * DENDRO_18 -
                               At1[pp] * DENDRO_15 + At2[pp] * DENDRO_29 +
                               At2[pp] * DENDRO_46 - At4[pp] * DENDRO_19 -
                               At4[pp] * DENDRO_48 + At5[pp] * DENDRO_39 +
                               DENDRO_14 * DENDRO_45 + DENDRO_14 * DENDRO_47;
            double DENDRO_50 = -DENDRO_41 + DENDRO_42 + DENDRO_43;
            double DENDRO_51 = DENDRO_14 * DENDRO_38;
            double DENDRO_52 = At5[pp] * DENDRO_38;
            double DENDRO_53 = At0[pp] * DENDRO_29 - At1[pp] * DENDRO_19 -
                               At1[pp] * DENDRO_48 + At2[pp] * DENDRO_39 +
                               At2[pp] * DENDRO_51 - At4[pp] * DENDRO_17 -
                               At4[pp] * DENDRO_18 * DENDRO_38 +
                               DENDRO_16 * DENDRO_47 + DENDRO_16 * DENDRO_52;
            double DENDRO_54 = DENDRO_41 - DENDRO_42 + DENDRO_43;
            double DENDRO_55 = 1.0 * At0[pp] * DENDRO_11 * DENDRO_8 -
                               1.0 * At1[pp] * DENDRO_11 * DENDRO_16 -
                               1.0 * At1[pp] * DENDRO_14 * DENDRO_8 +
                               1.0 * At2[pp] * DENDRO_11 * DENDRO_38 +
                               1.0 * At2[pp] * DENDRO_28 +
                               1.0 * At3[pp] * DENDRO_14 * DENDRO_16 -
                               1.0 * At4[pp] * DENDRO_14 * DENDRO_38 -
                               1.0 * At4[pp] * DENDRO_16 * DENDRO_8 +
                               1.0 * At5[pp] * DENDRO_38 * DENDRO_8;
            double DENDRO_56 =
                0.5 * gt0[pp] * gt4[pp] - 0.5 * gt1[pp] * gt2[pp];
            // Dendro: printing variables
            //--
            Gt_rhs_s6_0[pp] =
                DENDRO_5 *
                (DENDRO_22 * (-DENDRO_10 * DENDRO_8 + DENDRO_11 * DENDRO_13 +
                              DENDRO_6 * DENDRO_7) -
                 DENDRO_31 * (0.5 * DENDRO_11 * DENDRO_23 -
                              DENDRO_14 * DENDRO_27 + DENDRO_25 * DENDRO_8) -
                 DENDRO_40 * (-DENDRO_11 * DENDRO_37 + DENDRO_14 * DENDRO_35 +
                              DENDRO_32 * DENDRO_33) +
                 DENDRO_49 * (DENDRO_11 * DENDRO_26 - DENDRO_12 * DENDRO_14 +
                              DENDRO_44 * DENDRO_8) +
                 DENDRO_53 * (DENDRO_11 * DENDRO_50 - DENDRO_14 * DENDRO_9 +
                              DENDRO_34 * DENDRO_8) -
                 DENDRO_55 * (DENDRO_11 * DENDRO_24 - DENDRO_14 * DENDRO_54 +
                              DENDRO_36 * DENDRO_8));
            //--
            Gt_rhs_s6_1[pp] =
                DENDRO_5 *
                (-DENDRO_22 * (-DENDRO_10 * DENDRO_16 + DENDRO_13 * DENDRO_14 +
                               0.5 * DENDRO_18 * DENDRO_6) +
                 DENDRO_31 * (DENDRO_16 * DENDRO_25 - DENDRO_18 * DENDRO_27 +
                              DENDRO_23 * DENDRO_7) +
                 DENDRO_40 * (-DENDRO_14 * DENDRO_37 + DENDRO_18 * DENDRO_35 +
                              DENDRO_32 * DENDRO_56) -
                 1.0 * DENDRO_49 *
                     (-DENDRO_12 * DENDRO_18 + DENDRO_14 * DENDRO_26 +
                      DENDRO_16 * DENDRO_44) -
                 1.0 * DENDRO_53 *
                     (DENDRO_14 * DENDRO_50 + DENDRO_16 * DENDRO_34 -
                      DENDRO_18 * DENDRO_9) +
                 (DENDRO_14 * DENDRO_24 + DENDRO_16 * DENDRO_36 -
                  DENDRO_18 * DENDRO_54) *
                     (-At1[pp] * DENDRO_29 - At1[pp] * DENDRO_46 +
                      At2[pp] * DENDRO_11 * DENDRO_38 + At2[pp] * DENDRO_28 +
                      At3[pp] * DENDRO_19 - At4[pp] * DENDRO_39 -
                      At4[pp] * DENDRO_51 + DENDRO_45 * DENDRO_8 +
                      DENDRO_52 * DENDRO_8));
            //--
            Gt_rhs_s6_2[pp] =
                DENDRO_5 *
                (DENDRO_22 * (-DENDRO_10 * DENDRO_38 + DENDRO_13 * DENDRO_8 +
                              DENDRO_56 * DENDRO_6) -
                 DENDRO_31 * (-DENDRO_16 * DENDRO_27 + DENDRO_23 * DENDRO_33 +
                              DENDRO_25 * DENDRO_38) -
                 DENDRO_40 *
                     (DENDRO_16 * DENDRO_35 + 0.5 * DENDRO_32 * DENDRO_38 -
                      DENDRO_37 * DENDRO_8) +
                 DENDRO_49 * (-DENDRO_12 * DENDRO_16 + DENDRO_26 * DENDRO_8 +
                              DENDRO_38 * DENDRO_44) +
                 DENDRO_53 * (-DENDRO_16 * DENDRO_9 + DENDRO_34 * DENDRO_38 +
                              DENDRO_50 * DENDRO_8) -
                 DENDRO_55 * (-DENDRO_16 * DENDRO_54 + DENDRO_24 * DENDRO_8 +
                              DENDRO_36 * DENDRO_38));
            // Dendro: reduced ops:  364
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
