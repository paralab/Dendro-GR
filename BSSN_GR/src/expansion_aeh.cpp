// Dendro: {{{ 
// Dendro: original ops: 4034 
// Dendro: printing temp variables
const double DENDRO_0 = grad_2_F[pp];
const double DENDRO_1 = gt1[pp]*gt4[pp];
const double DENDRO_2 = gt2[pp]*gt3[pp];
const double DENDRO_3 = DENDRO_1 - DENDRO_2;
const double DENDRO_4 = grad_1_F[pp];
const double DENDRO_5 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
const double DENDRO_6 = grad_0_F[pp];
const double DENDRO_7 = gt3[pp]*gt5[pp];
const double DENDRO_8 = pow(gt4[pp], 2);
const double DENDRO_9 = DENDRO_0*DENDRO_3 - DENDRO_4*DENDRO_5 + DENDRO_6*(DENDRO_7 - DENDRO_8);
const double DENDRO_10 = gt0[pp]*gt3[pp];
const double DENDRO_11 = pow(gt1[pp], 2);
const double DENDRO_12 = gt0[pp]*gt4[pp] - gt1[pp]*gt2[pp];
const double DENDRO_13 = DENDRO_0*(DENDRO_10 - DENDRO_11) - DENDRO_12*DENDRO_4 + DENDRO_3*DENDRO_6;
const double DENDRO_14 = gt0[pp]*gt5[pp];
const double DENDRO_15 = pow(gt2[pp], 2);
const double DENDRO_16 = DENDRO_0*DENDRO_12 - DENDRO_4*(DENDRO_14 - DENDRO_15) + DENDRO_5*DENDRO_6;
const double DENDRO_17 = DENDRO_0*DENDRO_13 - DENDRO_16*DENDRO_4 + DENDRO_6*DENDRO_9;
const double DENDRO_18 = 1.0/DENDRO_17;
const double DENDRO_19 = -DENDRO_7 + DENDRO_8;
const double DENDRO_20 = 1.0/chi[pp];
const double DENDRO_21 = -DENDRO_1 + DENDRO_2;
const double DENDRO_22 = 0.5*grad_0_gt0[pp];
const double DENDRO_23 = -DENDRO_10 + DENDRO_11;
const double DENDRO_24 = grad_2_gt0[pp];
const double DENDRO_25 = -0.5*DENDRO_24 + 1.0*grad_0_gt2[pp];
const double DENDRO_26 = grad_1_gt0[pp];
const double DENDRO_27 = -0.5*DENDRO_26 + 1.0*grad_0_gt1[pp];
const double DENDRO_28 = grad_2_chi[pp];
const double DENDRO_29 = grad_1_chi[pp];
const double DENDRO_30 = grad_0_chi[pp];
const double DENDRO_31 = DENDRO_12*DENDRO_29 + DENDRO_21*DENDRO_30 + DENDRO_23*DENDRO_28;
const double DENDRO_32 = DENDRO_20*DENDRO_31;
const double DENDRO_33 = 0.5*gt0[pp];
const double DENDRO_34 = 1.0/(-2*DENDRO_1*gt2[pp] - DENDRO_10*gt5[pp] + DENDRO_11*gt5[pp] + DENDRO_15*gt3[pp] + DENDRO_8*gt0[pp]);
const double DENDRO_35 = DENDRO_0*DENDRO_34;
const double DENDRO_36 = -DENDRO_14 + DENDRO_15;
const double DENDRO_37 = DENDRO_12*DENDRO_28 + DENDRO_29*DENDRO_36 + DENDRO_30*DENDRO_5;
const double DENDRO_38 = DENDRO_20*DENDRO_37;
const double DENDRO_39 = DENDRO_34*DENDRO_4;
const double DENDRO_40 = DENDRO_19*DENDRO_34;
const double DENDRO_41 = DENDRO_21*DENDRO_34;
const double DENDRO_42 = DENDRO_34*DENDRO_5;
const double DENDRO_43 = DENDRO_19*DENDRO_30 + DENDRO_21*DENDRO_28 + DENDRO_29*DENDRO_5;
const double DENDRO_44 = DENDRO_34*chi[pp];
const double DENDRO_45 = 3/sqrt(-DENDRO_17*DENDRO_44);
const double DENDRO_46 = 0.5*grad_1_gt3[pp];
const double DENDRO_47 = grad_2_gt3[pp];
const double DENDRO_48 = -0.5*DENDRO_47 + 1.0*grad_1_gt4[pp];
const double DENDRO_49 = grad_0_gt3[pp];
const double DENDRO_50 = -0.5*DENDRO_49 + 1.0*grad_1_gt1[pp];
const double DENDRO_51 = 0.5*gt3[pp];
const double DENDRO_52 = DENDRO_20*DENDRO_43;
const double DENDRO_53 = DENDRO_34*DENDRO_6;
const double DENDRO_54 = DENDRO_34*DENDRO_36;
const double DENDRO_55 = DENDRO_12*DENDRO_34;
const double DENDRO_56 = 0.5*grad_2_gt5[pp];
const double DENDRO_57 = grad_1_gt5[pp];
const double DENDRO_58 = -0.5*DENDRO_57 + 1.0*grad_2_gt4[pp];
const double DENDRO_59 = grad_0_gt5[pp];
const double DENDRO_60 = -0.5*DENDRO_59 + 1.0*grad_2_gt2[pp];
const double DENDRO_61 = 0.5*gt5[pp];
const double DENDRO_62 = DENDRO_23*DENDRO_34;
const double DENDRO_63 = DENDRO_18*DENDRO_9;
const double DENDRO_64 = grad_0_gt4[pp];
const double DENDRO_65 = grad_2_gt1[pp];
const double DENDRO_66 = grad_1_gt2[pp];
const double DENDRO_67 = DENDRO_64 + DENDRO_65 - DENDRO_66;
const double DENDRO_68 = DENDRO_34*gt2[pp];
const double DENDRO_69 = 0.5*DENDRO_0;
const double DENDRO_70 = 0.5*DENDRO_6;
const double DENDRO_71 = DENDRO_64 - DENDRO_65 + DENDRO_66;
const double DENDRO_72 = DENDRO_34*gt1[pp];
const double DENDRO_73 = 0.5*DENDRO_4;
const double DENDRO_74 = -DENDRO_64 + DENDRO_65 + DENDRO_66;
const double DENDRO_75 = DENDRO_34*gt4[pp];

// Dendro: printing variables
//--
H[pp] = -1.0/3.0*DENDRO_44*(2*(DENDRO_12 - DENDRO_13*DENDRO_16*DENDRO_18)*(DENDRO_20*(3*At4[pp] + K[pp]*gt4[pp]) + DENDRO_45*(0.5*DENDRO_53*(DENDRO_19*DENDRO_74 + DENDRO_21*DENDRO_57 + DENDRO_47*DENDRO_5 + DENDRO_52*gt4[pp]) + DENDRO_69*(-DENDRO_20*(DENDRO_29 - DENDRO_31*DENDRO_75) + DENDRO_41*DENDRO_74 + DENDRO_47*DENDRO_55 + DENDRO_57*DENDRO_62) + DENDRO_73*(-DENDRO_20*(DENDRO_28 - DENDRO_37*DENDRO_75) + DENDRO_42*DENDRO_74 + DENDRO_47*DENDRO_54 + DENDRO_55*DENDRO_57) - grad2_1_2_F[pp])) + 2*(DENDRO_13*DENDRO_63 + DENDRO_21)*(DENDRO_20*(3*At2[pp] + K[pp]*gt2[pp]) + DENDRO_45*(0.5*DENDRO_39*(DENDRO_12*DENDRO_59 + DENDRO_24*DENDRO_5 + DENDRO_36*DENDRO_67 + DENDRO_38*gt2[pp]) + DENDRO_69*(-DENDRO_20*(DENDRO_30 - DENDRO_31*DENDRO_68) + DENDRO_24*DENDRO_41 + DENDRO_55*DENDRO_67 + DENDRO_59*DENDRO_62) + DENDRO_70*(-DENDRO_20*(DENDRO_28 - DENDRO_43*DENDRO_68) + DENDRO_24*DENDRO_40 + DENDRO_41*DENDRO_59 + DENDRO_42*DENDRO_67) - grad2_0_2_F[pp])) + ((DENDRO_13 * DENDRO_13)*DENDRO_18 + DENDRO_23)*(DENDRO_20*(3*At5[pp] + K[pp]*gt5[pp]) + DENDRO_45*(DENDRO_0*(-DENDRO_20*(1.0*DENDRO_28 - DENDRO_31*DENDRO_34*DENDRO_61) + DENDRO_41*DENDRO_60 + DENDRO_55*DENDRO_58 + DENDRO_56*DENDRO_62) + DENDRO_39*(DENDRO_12*DENDRO_56 + DENDRO_36*DENDRO_58 + DENDRO_38*DENDRO_61 + DENDRO_5*DENDRO_60) + DENDRO_53*(DENDRO_19*DENDRO_60 + DENDRO_21*DENDRO_56 + DENDRO_5*DENDRO_58 + DENDRO_52*DENDRO_61) - grad2_2_2_F[pp])) + 2*(-DENDRO_16*DENDRO_63 + DENDRO_5)*(DENDRO_20*(3*At1[pp] + K[pp]*gt1[pp]) + DENDRO_45*(0.5*DENDRO_35*(DENDRO_12*DENDRO_49 + DENDRO_21*DENDRO_26 + DENDRO_23*DENDRO_71 + DENDRO_32*gt1[pp]) + DENDRO_70*(-DENDRO_20*(DENDRO_29 - DENDRO_43*DENDRO_72) + DENDRO_26*DENDRO_40 + DENDRO_41*DENDRO_71 + DENDRO_42*DENDRO_49) + DENDRO_73*(-DENDRO_20*(DENDRO_30 - DENDRO_37*DENDRO_72) + DENDRO_26*DENDRO_42 + DENDRO_49*DENDRO_54 + DENDRO_55*DENDRO_71) - grad2_0_1_F[pp])) + ((DENDRO_16 * DENDRO_16)*DENDRO_18 + DENDRO_36)*(DENDRO_20*(3*At3[pp] + K[pp]*gt3[pp]) + DENDRO_45*(DENDRO_35*(DENDRO_12*DENDRO_46 + DENDRO_21*DENDRO_50 + DENDRO_23*DENDRO_48 + DENDRO_32*DENDRO_51) + DENDRO_4*(-DENDRO_20*(1.0*DENDRO_29 - DENDRO_34*DENDRO_37*DENDRO_51) + DENDRO_42*DENDRO_50 + DENDRO_46*DENDRO_54 + DENDRO_48*DENDRO_55) + DENDRO_53*(DENDRO_19*DENDRO_50 + DENDRO_21*DENDRO_48 + DENDRO_46*DENDRO_5 + DENDRO_51*DENDRO_52) - grad2_1_1_F[pp])) + (DENDRO_18*(DENDRO_9 * DENDRO_9) + DENDRO_19)*(DENDRO_20*(3*At0[pp] + K[pp]*gt0[pp]) + DENDRO_45*(DENDRO_35*(DENDRO_12*DENDRO_27 + DENDRO_21*DENDRO_22 + DENDRO_23*DENDRO_25 + DENDRO_32*DENDRO_33) + DENDRO_39*(DENDRO_12*DENDRO_25 + DENDRO_22*DENDRO_5 + DENDRO_27*DENDRO_36 + DENDRO_33*DENDRO_38) + DENDRO_6*(-DENDRO_20*(1.0*DENDRO_30 - DENDRO_33*DENDRO_34*DENDRO_43) + DENDRO_22*DENDRO_40 + DENDRO_25*DENDRO_41 + DENDRO_27*DENDRO_42) - grad2_0_0_F[pp])));
// Dendro: reduced ops: 413
// Dendro: }}} 
