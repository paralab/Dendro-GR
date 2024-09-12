
#ifdef USE_ETA_FUNC

bssn::timer::t_rhs.start();
for (unsigned int k = 3; k < nz - 3; k++) {
    z = pmin[2] + k * hz;
    for (unsigned int j = 3; j < ny - 3; j++) {
        y = pmin[1] + j * hy;
        for (unsigned int i = 3; i < nx - 3; i++) {
            x       = pmin[0] + i * hx;
            pp      = i + nx * (j + ny * k);
            r_coord = sqrt(x * x + y * y + z * z);
            eta     = ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
            }
            const double R0          = BSSN_ETA_R0;
            const double eta_power[] = {BSSN_ETA_POWER[0], BSSN_ETA_POWER[1]};
            // Dendro: {{{
            // Dendro: original ops:  516
            // Dendro: printing temp variables
            double DENDRO_0          = pow(gt4[pp], 2);
            double DENDRO_1          = pow(gt1[pp], 2);
            double DENDRO_2          = pow(gt2[pp], 2);
            double DENDRO_3          = gt0[pp] * gt3[pp];
            double DENDRO_4          = gt1[pp] * gt2[pp];
            double DENDRO_5          = grad_2_chi[pp];
            double DENDRO_6          = grad_1_chi[pp];
            double DENDRO_7          = grad_0_chi[pp];
            double DENDRO_8          = 2 * DENDRO_5;
            double DENDRO_9 =
                R0 *
                sqrt((-pow(DENDRO_5, 2) * (-DENDRO_1 + DENDRO_3) -
                      pow(DENDRO_6, 2) * (-DENDRO_2 + gt0[pp] * gt5[pp]) +
                      2 * DENDRO_6 * DENDRO_7 *
                          (gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp]) +
                      DENDRO_6 * DENDRO_8 * (-DENDRO_4 + gt0[pp] * gt4[pp]) -
                      pow(DENDRO_7, 2) * (-DENDRO_0 + gt3[pp] * gt5[pp]) -
                      DENDRO_7 * DENDRO_8 *
                          (gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp])) /
                     (DENDRO_0 * gt0[pp] + DENDRO_1 * gt5[pp] +
                      DENDRO_2 * gt3[pp] - DENDRO_3 * gt5[pp] -
                      2 * DENDRO_4 * gt4[pp])) *
                pow(-pow(chi[pp], eta_power[0]) + 1, -eta_power[1]);
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

#else
const double R0 = ETA_R0;
bssn::timer::t_rhs.start();
for (unsigned int k = 3; k < nz - 3; k++) {
    z = pmin[2] + k * hz;
    for (unsigned int j = 3; j < ny - 3; j++) {
        y = pmin[1] + j * hy;
        for (unsigned int i = 3; i < nx - 3; i++) {
            x       = pmin[0] + i * hx;
            pp      = i + nx * (j + ny * k);
            r_coord = sqrt(x * x + y * y + z * z);
            eta     = ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
            }

            // Dendro: {{{
            // Dendro: original ops: 17226
            // Dendro: printing temp variables
            const double DENDRO_0 = beta0[pp] * agrad_0_Gt0[pp] +
                                    beta1[pp] * agrad_1_Gt0[pp] +
                                    beta2[pp] * agrad_2_Gt0[pp];
            const double DENDRO_1 =
                2 * gt0[pp] * gt4[pp] - 2 * gt1[pp] * gt2[pp];
            const double DENDRO_2 = pow(gt4[pp], 2);
            const double DENDRO_3 = pow(gt1[pp], 2);
            const double DENDRO_4 = pow(gt2[pp], 2);
            const double DENDRO_5 = gt0[pp] * gt3[pp];
            const double DENDRO_6 = gt1[pp] * gt2[pp];
            const double DENDRO_7 = DENDRO_2 * gt0[pp] + DENDRO_3 * gt5[pp] +
                                    DENDRO_4 * gt3[pp] - DENDRO_5 * gt5[pp] -
                                    2 * DENDRO_6 * gt4[pp];
            const double DENDRO_8  = 1.0 / DENDRO_7;
            const double DENDRO_9  = grad2_0_2_beta0[pp];
            const double DENDRO_10 = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
            const double DENDRO_11 = (7.0L / 3.0L) * DENDRO_10 * DENDRO_8;
            const double DENDRO_12 = grad2_1_2_beta1[pp];
            const double DENDRO_13 = (1.0L / 3.0L) * DENDRO_10 * DENDRO_8;
            const double DENDRO_14 = grad2_2_2_beta2[pp];
            const double DENDRO_15 = grad2_0_1_beta0[pp];
            const double DENDRO_16 = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
            const double DENDRO_17 = (7.0L / 3.0L) * DENDRO_16 * DENDRO_8;
            const double DENDRO_18 = grad2_1_1_beta1[pp];
            const double DENDRO_19 = (1.0L / 3.0L) * DENDRO_16 * DENDRO_8;
            const double DENDRO_20 = grad2_1_2_beta2[pp];
            const double DENDRO_21 = -DENDRO_3 + DENDRO_5;
            const double DENDRO_22 = DENDRO_21 * DENDRO_8;
            const double DENDRO_23 = -DENDRO_4 + gt0[pp] * gt5[pp];
            const double DENDRO_24 = DENDRO_23 * DENDRO_8;
            const double DENDRO_25 = grad2_0_0_beta0[pp];
            const double DENDRO_26 = -DENDRO_2 + gt3[pp] * gt5[pp];
            const double DENDRO_27 = DENDRO_26 * DENDRO_8;
            const double DENDRO_28 = grad2_0_1_beta1[pp];
            const double DENDRO_29 = (1.0L / 3.0L) * DENDRO_26 * DENDRO_8;
            const double DENDRO_30 = grad2_0_2_beta2[pp];
            const double DENDRO_31 = pow(DENDRO_7, -2);
            const double DENDRO_32 = 2 * DENDRO_31 * grad_0_alpha[pp];
            const double DENDRO_33 = pow(DENDRO_16, 2);
            const double DENDRO_34 = pow(DENDRO_10, 2);
            const double DENDRO_35 =
                2 * gt1[pp] * gt5[pp] - 2 * gt2[pp] * gt4[pp];
            const double DENDRO_36 =
                2 * gt1[pp] * gt4[pp] - 2 * gt2[pp] * gt3[pp];
            const double DENDRO_37 =
                At0[pp] * pow(DENDRO_26, 2) - At1[pp] * DENDRO_26 * DENDRO_35 +
                At2[pp] * DENDRO_26 * DENDRO_36 + At3[pp] * DENDRO_33 -
                At4[pp] * DENDRO_10 * DENDRO_35 + At5[pp] * DENDRO_34;
            const double DENDRO_38 = grad_2_chi[pp];
            const double DENDRO_39 = grad_1_chi[pp];
            const double DENDRO_40 = grad_0_chi[pp];
            const double DENDRO_41 = 2 * DENDRO_38;
            const double DENDRO_42 = -DENDRO_6 + gt0[pp] * gt4[pp];
            const double DENDRO_43 =
                R0 *
                sqrt(DENDRO_8 * (-DENDRO_10 * DENDRO_40 * DENDRO_41 -
                                 DENDRO_21 * pow(DENDRO_38, 2) -
                                 DENDRO_23 * pow(DENDRO_39, 2) -
                                 DENDRO_26 * pow(DENDRO_40, 2) +
                                 DENDRO_35 * DENDRO_39 * DENDRO_40 +
                                 DENDRO_39 * DENDRO_41 * DENDRO_42)) *
                pow(-pow(chi[pp], eta_power[0]) + 1, -eta_power[1]);
            const double DENDRO_44 = (1.0L / 3.0L) * DENDRO_8 * alpha[pp];
            const double DENDRO_45 = grad_0_K[pp];
            const double DENDRO_46 = 1.0 / chi[pp];
            const double DENDRO_47 = 9 * DENDRO_40 * DENDRO_46 * DENDRO_8;
            const double DENDRO_48 = grad_0_gt0[pp];
            const double DENDRO_49 = grad_1_gt0[pp];
            const double DENDRO_50 = -0.5 * DENDRO_49 + 1.0 * grad_0_gt1[pp];
            const double DENDRO_51 = grad_2_gt0[pp];
            const double DENDRO_52 = -0.5 * DENDRO_51 + 1.0 * grad_0_gt2[pp];
            const double DENDRO_53 = -DENDRO_10 * DENDRO_52 +
                                     DENDRO_16 * DENDRO_50 -
                                     0.5 * DENDRO_26 * DENDRO_48;
            const double DENDRO_54 = pow(DENDRO_7, -3);
            const double DENDRO_55 = 2 * DENDRO_37 * DENDRO_54 * alpha[pp];
            const double DENDRO_56 = grad_1_gt3[pp];
            const double DENDRO_57 =
                0.5 * gt1[pp] * gt5[pp] - 0.5 * gt2[pp] * gt4[pp];
            const double DENDRO_58 = grad_2_gt3[pp];
            const double DENDRO_59 = -0.5 * DENDRO_58 + 1.0 * grad_1_gt4[pp];
            const double DENDRO_60 = grad_0_gt3[pp];
            const double DENDRO_61 = 0.5 * DENDRO_60 - 1.0 * grad_1_gt1[pp];
            const double DENDRO_62 = -DENDRO_10 * DENDRO_59 +
                                     DENDRO_26 * DENDRO_61 +
                                     DENDRO_56 * DENDRO_57;
            const double DENDRO_63 = pow(DENDRO_42, 2);
            const double DENDRO_64 = At1[pp] * DENDRO_23;
            const double DENDRO_65 =
                At0[pp] * DENDRO_33 + At2[pp] * DENDRO_35 * DENDRO_42 +
                At3[pp] * pow(DENDRO_23, 2) - At4[pp] * DENDRO_1 * DENDRO_23 +
                At5[pp] * DENDRO_63 - DENDRO_35 * DENDRO_64;
            const double DENDRO_66 = 2 * DENDRO_54 * DENDRO_65 * alpha[pp];
            const double DENDRO_67 = grad_2_gt5[pp];
            const double DENDRO_68 =
                0.5 * gt1[pp] * gt4[pp] - 0.5 * gt2[pp] * gt3[pp];
            const double DENDRO_69 = grad_1_gt5[pp];
            const double DENDRO_70 = 0.5 * DENDRO_69 - 1.0 * grad_2_gt4[pp];
            const double DENDRO_71 = grad_0_gt5[pp];
            const double DENDRO_72 = 0.5 * DENDRO_71 - 1.0 * grad_2_gt2[pp];
            const double DENDRO_73 = -DENDRO_16 * DENDRO_70 +
                                     DENDRO_26 * DENDRO_72 -
                                     DENDRO_67 * DENDRO_68;
            const double DENDRO_74 = DENDRO_10 * DENDRO_42;
            const double DENDRO_75 = At2[pp] * DENDRO_21;
            const double DENDRO_76 = At4[pp] * DENDRO_21;
            const double DENDRO_77 =
                At0[pp] * DENDRO_34 - 2 * At1[pp] * DENDRO_74 +
                At3[pp] * DENDRO_63 + At5[pp] * pow(DENDRO_21, 2) -
                DENDRO_1 * DENDRO_76 + DENDRO_36 * DENDRO_75;
            const double DENDRO_78 = 2 * DENDRO_54 * DENDRO_77 * alpha[pp];
            const double DENDRO_79 = 2 * DENDRO_31 * grad_2_alpha[pp];
            const double DENDRO_80 = DENDRO_16 * DENDRO_42;
            const double DENDRO_81 = At0[pp] * DENDRO_26;
            const double DENDRO_82 = DENDRO_10 * DENDRO_16;
            const double DENDRO_83 = At5[pp] * DENDRO_21;
            const double DENDRO_84 = DENDRO_26 * DENDRO_42;
            const double DENDRO_85 = DENDRO_16 * DENDRO_21;
            const double DENDRO_86 =
                -At1[pp] * DENDRO_82 - At1[pp] * DENDRO_84 +
                At2[pp] * DENDRO_34 + At3[pp] * DENDRO_80 -
                At4[pp] * DENDRO_74 - At4[pp] * DENDRO_85 +
                DENDRO_10 * DENDRO_81 + DENDRO_10 * DENDRO_83 +
                DENDRO_26 * DENDRO_75;
            const double DENDRO_87 = 2 * DENDRO_31 * grad_1_alpha[pp];
            const double DENDRO_88 = DENDRO_10 * DENDRO_23;
            const double DENDRO_89 = At3[pp] * DENDRO_23;
            const double DENDRO_90 = At1[pp] * DENDRO_33 - At2[pp] * DENDRO_82 -
                                     At2[pp] * DENDRO_84 + At4[pp] * DENDRO_80 +
                                     At4[pp] * DENDRO_88 - At5[pp] * DENDRO_74 -
                                     DENDRO_16 * DENDRO_81 -
                                     DENDRO_16 * DENDRO_89 +
                                     DENDRO_26 * DENDRO_64;
            const double DENDRO_91 = grad_0_gt4[pp];
            const double DENDRO_92 = grad_2_gt1[pp];
            const double DENDRO_93 = grad_1_gt2[pp];
            const double DENDRO_94 = DENDRO_91 + DENDRO_92 - DENDRO_93;
            const double DENDRO_95 = -DENDRO_10 * DENDRO_71 +
                                     DENDRO_16 * DENDRO_94 -
                                     DENDRO_26 * DENDRO_51;
            const double DENDRO_96 = 2.0 * DENDRO_54 * DENDRO_86 * alpha[pp];
            const double DENDRO_97 = grad_2_K[pp];
            const double DENDRO_98 =
                4 * gt1[pp] * gt4[pp] - 4 * gt2[pp] * gt3[pp];
            const double DENDRO_99  = 9 * DENDRO_38 * DENDRO_46 * DENDRO_8;
            const double DENDRO_100 = DENDRO_91 - DENDRO_92 + DENDRO_93;
            const double DENDRO_101 = -DENDRO_10 * DENDRO_100 +
                                      DENDRO_16 * DENDRO_60 -
                                      DENDRO_26 * DENDRO_49;
            const double DENDRO_102 = 2.0 * DENDRO_54 * DENDRO_90 * alpha[pp];
            const double DENDRO_103 = -DENDRO_91 + DENDRO_92 + DENDRO_93;
            const double DENDRO_104 = -DENDRO_10 * DENDRO_69 -
                                      DENDRO_103 * DENDRO_26 +
                                      DENDRO_16 * DENDRO_58;
            const double DENDRO_105 =
                -At0[pp] * DENDRO_82 + At1[pp] * DENDRO_80 +
                At1[pp] * DENDRO_88 - At2[pp] * DENDRO_74 -
                At2[pp] * DENDRO_85 + At4[pp] * DENDRO_63 +
                DENDRO_23 * DENDRO_76 - DENDRO_42 * DENDRO_83 -
                DENDRO_42 * DENDRO_89;
            const double DENDRO_106 = 2.0 * DENDRO_105 * DENDRO_54 * alpha[pp];
            const double DENDRO_107 = grad_1_K[pp];
            const double DENDRO_108 =
                4 * gt1[pp] * gt5[pp] - 4 * gt2[pp] * gt4[pp];
            const double DENDRO_109 = 9 * DENDRO_39 * DENDRO_46 * DENDRO_8;
            const double DENDRO_110 = DENDRO_103 * DENDRO_16 -
                                      DENDRO_23 * DENDRO_58 +
                                      DENDRO_42 * DENDRO_69;
            const double DENDRO_111 = DENDRO_100 * DENDRO_42 +
                                      DENDRO_16 * DENDRO_49 -
                                      DENDRO_23 * DENDRO_60;
            const double DENDRO_112 =
                1.0 * gt1[pp] * gt4[pp] - 1.0 * gt2[pp] * gt3[pp];
            const double DENDRO_113 = DENDRO_16 * DENDRO_51 -
                                      DENDRO_23 * DENDRO_94 +
                                      DENDRO_42 * DENDRO_71;
            const double DENDRO_114 =
                0.5 * gt0[pp] * gt4[pp] - 0.5 * gt1[pp] * gt2[pp];
            const double DENDRO_115 = DENDRO_114 * DENDRO_67 -
                                      DENDRO_16 * DENDRO_72 +
                                      DENDRO_23 * DENDRO_70;
            const double DENDRO_116 = -DENDRO_16 * DENDRO_61 -
                                      0.5 * DENDRO_23 * DENDRO_56 +
                                      DENDRO_42 * DENDRO_59;
            const double DENDRO_117 = -DENDRO_23 * DENDRO_50 +
                                      DENDRO_42 * DENDRO_52 +
                                      DENDRO_48 * DENDRO_57;
            const double DENDRO_118 =
                DENDRO_31 * (DENDRO_110 * DENDRO_42 + DENDRO_111 * DENDRO_16 -
                             DENDRO_112 * DENDRO_113 - DENDRO_115 * DENDRO_21 -
                             DENDRO_116 * DENDRO_23 - DENDRO_117 * DENDRO_26);
            const double DENDRO_119 = -DENDRO_10 * DENDRO_103 -
                                      DENDRO_21 * DENDRO_69 +
                                      DENDRO_42 * DENDRO_58;
            const double DENDRO_120 = -DENDRO_10 * DENDRO_49 -
                                      DENDRO_100 * DENDRO_21 +
                                      DENDRO_42 * DENDRO_60;
            const double DENDRO_121 = -DENDRO_10 * DENDRO_51 -
                                      DENDRO_21 * DENDRO_71 +
                                      DENDRO_42 * DENDRO_94;
            const double DENDRO_122 = DENDRO_10 * DENDRO_72 -
                                      0.5 * DENDRO_21 * DENDRO_67 -
                                      DENDRO_42 * DENDRO_70;
            const double DENDRO_123 = DENDRO_10 * DENDRO_61 +
                                      DENDRO_114 * DENDRO_56 -
                                      DENDRO_21 * DENDRO_59;
            const double DENDRO_124 = -DENDRO_21 * DENDRO_52 +
                                      DENDRO_42 * DENDRO_50 -
                                      DENDRO_48 * DENDRO_68;
            const double DENDRO_125 =
                DENDRO_31 * (-DENDRO_112 * DENDRO_121 + DENDRO_119 * DENDRO_42 +
                             DENDRO_120 * DENDRO_16 - DENDRO_122 * DENDRO_21 -
                             DENDRO_123 * DENDRO_23 - DENDRO_124 * DENDRO_26);
            const double DENDRO_126 = grad_0_beta0[pp];
            const double DENDRO_127 =
                DENDRO_31 * (DENDRO_101 * DENDRO_16 + DENDRO_104 * DENDRO_42 -
                             DENDRO_112 * DENDRO_95 - DENDRO_21 * DENDRO_73 -
                             DENDRO_23 * DENDRO_62 - DENDRO_26 * DENDRO_53);
            const double DENDRO_128 = grad_1_beta1[pp];
            const double DENDRO_129 = grad_2_beta2[pp];
            const double DENDRO_130 = (2.0L / 3.0L) * DENDRO_126 +
                                      (2.0L / 3.0L) * DENDRO_128 +
                                      (2.0L / 3.0L) * DENDRO_129;
            const double DENDRO_131 = beta0[pp] * agrad_0_Gt1[pp] +
                                      beta1[pp] * agrad_1_Gt1[pp] +
                                      beta2[pp] * agrad_2_Gt1[pp];
            const double DENDRO_132 = (1.0L / 3.0L) * DENDRO_42 * DENDRO_8;
            const double DENDRO_133 = (7.0L / 3.0L) * DENDRO_42 * DENDRO_8;
            const double DENDRO_134 = (1.0L / 3.0L) * DENDRO_23 * DENDRO_8;
            const double DENDRO_135 =
                4 * gt0[pp] * gt4[pp] - 4 * gt1[pp] * gt2[pp];
            const double DENDRO_136 = beta0[pp] * agrad_0_Gt2[pp] +
                                      beta1[pp] * agrad_1_Gt2[pp] +
                                      beta2[pp] * agrad_2_Gt2[pp];
            const double DENDRO_137 = (1.0L / 3.0L) * DENDRO_21 * DENDRO_8;
            // Dendro: printing variables

            B_rhs0[pp] =
                -B0[pp] * DENDRO_43 - DENDRO_0 * lambda[3] + DENDRO_0 +
                DENDRO_1 * DENDRO_8 * grad2_1_2_beta0[pp] +
                DENDRO_101 * DENDRO_102 + DENDRO_104 * DENDRO_106 -
                DENDRO_11 * DENDRO_9 - DENDRO_118 * grad_1_beta0[pp] -
                DENDRO_12 * DENDRO_13 - DENDRO_125 * grad_2_beta0[pp] -
                DENDRO_126 * DENDRO_127 + DENDRO_127 * DENDRO_130 -
                DENDRO_13 * DENDRO_14 + DENDRO_15 * DENDRO_17 +
                DENDRO_18 * DENDRO_19 + DENDRO_19 * DENDRO_20 -
                DENDRO_22 * grad2_2_2_beta0[pp] -
                DENDRO_24 * grad2_1_1_beta0[pp] -
                4.0L / 3.0L * DENDRO_25 * DENDRO_27 - DENDRO_28 * DENDRO_29 -
                DENDRO_29 * DENDRO_30 - DENDRO_32 * DENDRO_37 -
                DENDRO_44 * (DENDRO_107 * DENDRO_108 + DENDRO_109 * DENDRO_90) -
                DENDRO_44 *
                    (-4 * DENDRO_26 * DENDRO_45 + DENDRO_37 * DENDRO_47) -
                DENDRO_44 * (DENDRO_86 * DENDRO_99 - DENDRO_97 * DENDRO_98) +
                DENDRO_53 * DENDRO_55 + DENDRO_62 * DENDRO_66 +
                DENDRO_73 * DENDRO_78 - DENDRO_79 * DENDRO_86 -
                DENDRO_87 * DENDRO_90 + DENDRO_95 * DENDRO_96 +
                lambda[2] *
                    (beta0[pp] * agrad_0_B0[pp] + beta1[pp] * agrad_1_B0[pp] +
                     beta2[pp] * agrad_2_B0[pp]);
            B_rhs1[pp] =
                -B1[pp] * DENDRO_43 + DENDRO_102 * DENDRO_111 -
                DENDRO_105 * DENDRO_79 + DENDRO_106 * DENDRO_110 +
                DENDRO_113 * DENDRO_96 + DENDRO_115 * DENDRO_78 +
                DENDRO_116 * DENDRO_66 + DENDRO_117 * DENDRO_55 -
                DENDRO_118 * DENDRO_128 + DENDRO_118 * DENDRO_130 +
                DENDRO_12 * DENDRO_133 - DENDRO_125 * grad_2_beta1[pp] -
                DENDRO_127 * grad_0_beta1[pp] - DENDRO_131 * lambda[3] +
                DENDRO_131 + DENDRO_132 * DENDRO_14 + DENDRO_132 * DENDRO_9 -
                DENDRO_134 * DENDRO_15 - DENDRO_134 * DENDRO_20 +
                DENDRO_17 * DENDRO_28 - 4.0L / 3.0L * DENDRO_18 * DENDRO_24 +
                DENDRO_19 * DENDRO_25 + DENDRO_19 * DENDRO_30 -
                DENDRO_22 * grad2_2_2_beta1[pp] -
                DENDRO_27 * grad2_0_0_beta1[pp] - DENDRO_32 * DENDRO_90 -
                DENDRO_36 * DENDRO_8 * grad2_0_2_beta1[pp] -
                DENDRO_44 * (DENDRO_105 * DENDRO_99 + DENDRO_135 * DENDRO_97) -
                DENDRO_44 *
                    (-4 * DENDRO_107 * DENDRO_23 + DENDRO_109 * DENDRO_65) -
                DENDRO_44 * (DENDRO_108 * DENDRO_45 + DENDRO_47 * DENDRO_90) -
                DENDRO_65 * DENDRO_87 +
                lambda[2] *
                    (beta0[pp] * agrad_0_B1[pp] + beta1[pp] * agrad_1_B1[pp] +
                     beta2[pp] * agrad_2_B1[pp]);
            B_rhs2[pp] =
                -B2[pp] * DENDRO_43 + DENDRO_102 * DENDRO_120 -
                DENDRO_105 * DENDRO_87 + DENDRO_106 * DENDRO_119 -
                DENDRO_11 * DENDRO_30 - DENDRO_118 * grad_1_beta2[pp] -
                DENDRO_12 * DENDRO_137 + DENDRO_121 * DENDRO_96 +
                DENDRO_122 * DENDRO_78 + DENDRO_123 * DENDRO_66 +
                DENDRO_124 * DENDRO_55 - DENDRO_125 * DENDRO_129 +
                DENDRO_125 * DENDRO_130 - DENDRO_127 * grad_0_beta2[pp] -
                DENDRO_13 * DENDRO_25 - DENDRO_13 * DENDRO_28 +
                DENDRO_132 * DENDRO_15 + DENDRO_132 * DENDRO_18 +
                DENDRO_133 * DENDRO_20 - DENDRO_136 * lambda[3] + DENDRO_136 -
                DENDRO_137 * DENDRO_9 - 4.0L / 3.0L * DENDRO_14 * DENDRO_22 -
                DENDRO_24 * grad2_1_1_beta2[pp] -
                DENDRO_27 * grad2_0_0_beta2[pp] - DENDRO_32 * DENDRO_86 +
                DENDRO_35 * DENDRO_8 * grad2_0_1_beta2[pp] -
                DENDRO_44 *
                    (DENDRO_105 * DENDRO_109 + DENDRO_107 * DENDRO_135) -
                DENDRO_44 *
                    (-4 * DENDRO_21 * DENDRO_97 + DENDRO_77 * DENDRO_99) -
                DENDRO_44 * (-DENDRO_45 * DENDRO_98 + DENDRO_47 * DENDRO_86) -
                DENDRO_77 * DENDRO_79 +
                lambda[2] *
                    (beta0[pp] * agrad_0_B2[pp] + beta1[pp] * agrad_1_B2[pp] +
                     beta2[pp] * agrad_2_B2[pp]);
            // Dendro: reduced ops: 765
            // Dendro: }}}

            /*
            // Dendro: {{{
            // Dendro: original ops:  66
            // Dendro: printing temp variables
            // Dendro: printing variables
            //--
            B_rhs0[pp] = -B0[pp]*eta + Gt_rhs0[pp] +
            lambda[2]*(beta0[pp]*agrad_0_B0[pp] + beta1[pp]*agrad_1_B0[pp] +
            beta2[pp]*agrad_2_B0[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt0[pp] +
            beta1[pp]*agrad_1_Gt0[pp] + beta2[pp]*agrad_2_Gt0[pp]);
            //--
            B_rhs1[pp] = -B1[pp]*eta + Gt_rhs1[pp] +
            lambda[2]*(beta0[pp]*agrad_0_B1[pp] + beta1[pp]*agrad_1_B1[pp] +
            beta2[pp]*agrad_2_B1[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt1[pp] +
            beta1[pp]*agrad_1_Gt1[pp] + beta2[pp]*agrad_2_Gt1[pp]);
            //--
            B_rhs2[pp] = -B2[pp]*eta + Gt_rhs2[pp] +
            lambda[2]*(beta0[pp]*agrad_0_B2[pp] + beta1[pp]*agrad_1_B2[pp] +
            beta2[pp]*agrad_2_B2[pp]) - lambda[3]*(beta0[pp]*agrad_0_Gt2[pp] +
            beta1[pp]*agrad_1_Gt2[pp] + beta2[pp]*agrad_2_Gt2[pp]);
            // Dendro: reduced ops:  66
            // Dendro: }}}*/
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

#endif