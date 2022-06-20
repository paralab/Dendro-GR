    bssn::timer::t_rhs.start();
for (unsigned int k = PW; k < nz-PW; k++) { 
    z = pmin[2] + k*hz;
for (unsigned int j = PW; j < ny-PW; j++) { 
    y = pmin[1] + j*hy; 
for (unsigned int i = PW; i < nx-PW; i++) {
    x = pmin[0] + i*hx;
    pp = i + nx*(j + ny*k);
    r_coord = sqrt(x*x + y*y + z*z);
    eta=ETA_CONST;
    if (r_coord >= ETA_R0) {
    eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
    }
// Dendro: {{{ 
// Dendro: original ops:  3960
// Dendro: printing temp variables
double DENDRO_0 = pow(gt4[pp], 2);
double DENDRO_1 = pow(gt1[pp], 2);
double DENDRO_2 = pow(gt2[pp], 2);
double DENDRO_3 = gt0[pp]*gt3[pp];
double DENDRO_4 = gt1[pp]*gt2[pp];
double DENDRO_5 = DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt5[pp] - 2*DENDRO_4*gt4[pp];
double DENDRO_6 = 1.0/DENDRO_5;
double DENDRO_7 = DENDRO_6*chi[pp];
double DENDRO_8 = -DENDRO_1 + DENDRO_3;
double DENDRO_9 = grad_1_alpha[pp];
double DENDRO_10 = DENDRO_6*DENDRO_9;
double DENDRO_11 = grad_2_gt5[pp];
double DENDRO_12 = 0.5*gt0[pp]*gt4[pp] - 0.5*gt1[pp]*gt2[pp];
double DENDRO_13 = gt0[pp]*gt5[pp];
double DENDRO_14 = -DENDRO_13 + DENDRO_2;
double DENDRO_15 = grad_1_gt5[pp];
double DENDRO_16 = -0.5*DENDRO_15 + 1.0*grad_2_gt4[pp];
double DENDRO_17 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
double DENDRO_18 = grad_0_gt5[pp];
double DENDRO_19 = -0.5*DENDRO_18 + 1.0*grad_2_gt2[pp];
double DENDRO_20 = 0.5*gt5[pp];
double DENDRO_21 = 1.0/chi[pp];
double DENDRO_22 = grad_2_chi[pp];
double DENDRO_23 = -DENDRO_4 + gt0[pp]*gt4[pp];
double DENDRO_24 = grad_1_chi[pp];
double DENDRO_25 = grad_0_chi[pp];
double DENDRO_26 = DENDRO_14*DENDRO_24 + DENDRO_17*DENDRO_25 + DENDRO_22*DENDRO_23;
double DENDRO_27 = DENDRO_21*DENDRO_26;
double DENDRO_28 = grad_0_alpha[pp];
double DENDRO_29 = DENDRO_28*DENDRO_6;
double DENDRO_30 = -0.5*gt1[pp]*gt4[pp] + 0.5*gt2[pp]*gt3[pp];
double DENDRO_31 = gt3[pp]*gt5[pp];
double DENDRO_32 = DENDRO_0 - DENDRO_31;
double DENDRO_33 = gt2[pp]*gt3[pp];
double DENDRO_34 = gt1[pp]*gt4[pp];
double DENDRO_35 = DENDRO_33 - DENDRO_34;
double DENDRO_36 = DENDRO_17*DENDRO_24 + DENDRO_22*DENDRO_35 + DENDRO_25*DENDRO_32;
double DENDRO_37 = DENDRO_21*DENDRO_36;
double DENDRO_38 = grad_2_alpha[pp];
double DENDRO_39 = DENDRO_1 - DENDRO_3;
double DENDRO_40 = DENDRO_39*DENDRO_6;
double DENDRO_41 = DENDRO_23*DENDRO_6;
double DENDRO_42 = DENDRO_35*DENDRO_6;
double DENDRO_43 = DENDRO_22*DENDRO_39 + DENDRO_23*DENDRO_24 + DENDRO_25*DENDRO_35;
double DENDRO_44 = DENDRO_13 - DENDRO_2;
double DENDRO_45 = DENDRO_38*DENDRO_6;
double DENDRO_46 = grad_1_gt3[pp];
double DENDRO_47 = grad_2_gt3[pp];
double DENDRO_48 = -0.5*DENDRO_47 + 1.0*grad_1_gt4[pp];
double DENDRO_49 = grad_0_gt3[pp];
double DENDRO_50 = -0.5*DENDRO_49 + 1.0*grad_1_gt1[pp];
double DENDRO_51 = 0.5*gt3[pp];
double DENDRO_52 = DENDRO_21*DENDRO_43;
double DENDRO_53 = 0.5*gt1[pp]*gt5[pp] - 0.5*gt2[pp]*gt4[pp];
double DENDRO_54 = DENDRO_14*DENDRO_6;
double DENDRO_55 = DENDRO_17*DENDRO_6;
double DENDRO_56 = -DENDRO_0 + DENDRO_31;
double DENDRO_57 = grad_0_gt0[pp];
double DENDRO_58 = grad_2_gt0[pp];
double DENDRO_59 = -0.5*DENDRO_58 + 1.0*grad_0_gt2[pp];
double DENDRO_60 = grad_1_gt0[pp];
double DENDRO_61 = -0.5*DENDRO_60 + 1.0*grad_0_gt1[pp];
double DENDRO_62 = 0.5*gt0[pp];
double DENDRO_63 = DENDRO_32*DENDRO_6;
double DENDRO_64 = 2*chi[pp];
double DENDRO_65 = 0.5*DENDRO_6;
double DENDRO_66 = grad_1_gt2[pp];
double DENDRO_67 = grad_2_gt1[pp];
double DENDRO_68 = grad_0_gt4[pp];
double DENDRO_69 = DENDRO_66 + DENDRO_67 - DENDRO_68;
double DENDRO_70 = 0.5*DENDRO_38;
double DENDRO_71 = DENDRO_6*gt4[pp];
double DENDRO_72 = 0.5*DENDRO_9;
double DENDRO_73 = -DENDRO_33 + DENDRO_34;
double DENDRO_74 = -DENDRO_66 + DENDRO_67 + DENDRO_68;
double DENDRO_75 = DENDRO_6*gt2[pp];
double DENDRO_76 = 0.5*DENDRO_28;
double DENDRO_77 = DENDRO_66 - DENDRO_67 + DENDRO_68;
double DENDRO_78 = DENDRO_6*gt1[pp];
double DENDRO_79 = pow(DENDRO_5, -2);
double DENDRO_80 = 3*DENDRO_79;
double DENDRO_81 = pow(DENDRO_17, 2);
double DENDRO_82 = pow(DENDRO_73, 2);
double DENDRO_83 = DENDRO_17*DENDRO_73;
double DENDRO_84 = 2*At1[pp]*DENDRO_17;
double DENDRO_85 = 2*At2[pp]*DENDRO_73;
double DENDRO_86 = pow(DENDRO_23, 2);
double DENDRO_87 = DENDRO_17*DENDRO_23;
double DENDRO_88 = 2*At4[pp]*DENDRO_23;
double DENDRO_89 = DENDRO_23*DENDRO_73;
double DENDRO_90 = 6*DENDRO_79;
double DENDRO_91 = At0[pp]*DENDRO_56;
double DENDRO_92 = At5[pp]*DENDRO_8;
double DENDRO_93 = DENDRO_23*DENDRO_56;
double DENDRO_94 = DENDRO_17*DENDRO_8;
double DENDRO_95 = DENDRO_44*DENDRO_73;
double DENDRO_96 = At3[pp]*DENDRO_44;
// Dendro: printing variables
//--
K_rhs[pp] = DENDRO_17*DENDRO_6*DENDRO_64*(DENDRO_38*DENDRO_65*(DENDRO_23*DENDRO_49 + DENDRO_35*DENDRO_60 + DENDRO_39*DENDRO_77 + DENDRO_52*gt1[pp]) + DENDRO_72*(-DENDRO_21*(DENDRO_25 - DENDRO_26*DENDRO_78) + DENDRO_41*DENDRO_77 + DENDRO_49*DENDRO_54 + DENDRO_55*DENDRO_60) + DENDRO_76*(-DENDRO_21*(DENDRO_24 - DENDRO_36*DENDRO_78) + DENDRO_42*DENDRO_77 + DENDRO_49*DENDRO_55 + DENDRO_60*DENDRO_63) - grad2_0_1_alpha[pp]) + DENDRO_23*DENDRO_6*DENDRO_64*(DENDRO_28*DENDRO_65*(DENDRO_15*DENDRO_35 + DENDRO_17*DENDRO_47 + DENDRO_32*DENDRO_69 + DENDRO_37*gt4[pp]) + DENDRO_70*(DENDRO_15*DENDRO_40 - DENDRO_21*(DENDRO_24 - DENDRO_43*DENDRO_71) + DENDRO_41*DENDRO_47 + DENDRO_42*DENDRO_69) + DENDRO_72*(DENDRO_15*DENDRO_41 - DENDRO_21*(DENDRO_22 - DENDRO_26*DENDRO_71) + DENDRO_47*DENDRO_54 + DENDRO_55*DENDRO_69) - grad2_1_2_alpha[pp]) - DENDRO_44*DENDRO_7*(DENDRO_29*(DENDRO_32*DENDRO_50 + DENDRO_35*DENDRO_48 + DENDRO_37*DENDRO_51 + DENDRO_46*DENDRO_53) + DENDRO_45*(DENDRO_12*DENDRO_46 + DENDRO_35*DENDRO_50 + DENDRO_39*DENDRO_48 + DENDRO_51*DENDRO_52) + DENDRO_9*(-DENDRO_21*(1.0*DENDRO_24 - DENDRO_26*DENDRO_51*DENDRO_6) + DENDRO_41*DENDRO_48 + 0.5*DENDRO_46*DENDRO_54 + DENDRO_50*DENDRO_55) - grad2_1_1_alpha[pp]) - DENDRO_56*DENDRO_7*(DENDRO_10*(DENDRO_14*DENDRO_61 + DENDRO_23*DENDRO_59 + DENDRO_27*DENDRO_62 + DENDRO_53*DENDRO_57) + DENDRO_28*(-DENDRO_21*(1.0*DENDRO_25 - DENDRO_36*DENDRO_6*DENDRO_62) + DENDRO_42*DENDRO_59 + DENDRO_55*DENDRO_61 + 0.5*DENDRO_57*DENDRO_63) + DENDRO_45*(DENDRO_23*DENDRO_61 + DENDRO_30*DENDRO_57 + DENDRO_39*DENDRO_59 + DENDRO_52*DENDRO_62) - grad2_0_0_alpha[pp]) - 2*DENDRO_7*DENDRO_73*(DENDRO_65*DENDRO_9*(DENDRO_14*DENDRO_74 + DENDRO_17*DENDRO_58 + DENDRO_18*DENDRO_23 + DENDRO_27*gt2[pp]) + DENDRO_70*(DENDRO_18*DENDRO_40 - DENDRO_21*(DENDRO_25 - DENDRO_43*DENDRO_75) + DENDRO_41*DENDRO_74 + DENDRO_42*DENDRO_58) + DENDRO_76*(DENDRO_18*DENDRO_42 - DENDRO_21*(DENDRO_22 - DENDRO_36*DENDRO_75) + DENDRO_55*DENDRO_74 + DENDRO_58*DENDRO_63) - grad2_0_2_alpha[pp]) - DENDRO_7*DENDRO_8*(DENDRO_10*(DENDRO_11*DENDRO_12 + DENDRO_14*DENDRO_16 + DENDRO_17*DENDRO_19 + DENDRO_20*DENDRO_27) + DENDRO_29*(DENDRO_11*DENDRO_30 + DENDRO_16*DENDRO_17 + DENDRO_19*DENDRO_32 + DENDRO_20*DENDRO_37) + DENDRO_38*(0.5*DENDRO_11*DENDRO_40 + DENDRO_16*DENDRO_41 + DENDRO_19*DENDRO_42 - DENDRO_21*(-DENDRO_20*DENDRO_43*DENDRO_6 + 1.0*DENDRO_22)) - grad2_2_2_alpha[pp]) + (1.0L/3.0L)*alpha[pp]*(At0[pp]*DENDRO_80*(At0[pp]*pow(DENDRO_56, 2) + At3[pp]*DENDRO_81 - 2*At4[pp]*DENDRO_83 + At5[pp]*DENDRO_82 - DENDRO_56*DENDRO_84 + DENDRO_56*DENDRO_85) + At1[pp]*DENDRO_90*(At1[pp]*DENDRO_44*DENDRO_56 + At1[pp]*DENDRO_81 - At2[pp]*DENDRO_83 - At2[pp]*DENDRO_93 + At4[pp]*DENDRO_87 + At4[pp]*DENDRO_95 - At5[pp]*DENDRO_89 - DENDRO_17*DENDRO_91 - DENDRO_17*DENDRO_96) + At2[pp]*DENDRO_90*(-At1[pp]*DENDRO_83 - At1[pp]*DENDRO_93 + At2[pp]*DENDRO_56*DENDRO_8 + At2[pp]*DENDRO_82 + At3[pp]*DENDRO_87 - At4[pp]*DENDRO_89 - At4[pp]*DENDRO_94 + DENDRO_73*DENDRO_91 + DENDRO_73*DENDRO_92) + At3[pp]*DENDRO_80*(At0[pp]*DENDRO_81 + 2*At2[pp]*DENDRO_87 + At3[pp]*pow(DENDRO_44, 2) + At5[pp]*DENDRO_86 - DENDRO_44*DENDRO_84 - DENDRO_44*DENDRO_88) + At4[pp]*DENDRO_90*(-At0[pp]*DENDRO_83 + At1[pp]*DENDRO_87 + At1[pp]*DENDRO_95 - At2[pp]*DENDRO_89 - At2[pp]*DENDRO_94 + At4[pp]*DENDRO_44*DENDRO_8 + At4[pp]*DENDRO_86 - DENDRO_23*DENDRO_92 - DENDRO_23*DENDRO_96) + At5[pp]*DENDRO_80*(At0[pp]*DENDRO_82 - 2*At1[pp]*DENDRO_89 + At3[pp]*DENDRO_86 + At5[pp]*pow(DENDRO_8, 2) + DENDRO_8*DENDRO_85 - DENDRO_8*DENDRO_88) + pow(K[pp], 2)) + beta0[pp]*agrad_0_K[pp] + beta1[pp]*agrad_1_K[pp] + beta2[pp]*agrad_2_K[pp];
// Dendro: reduced ops:  500
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
