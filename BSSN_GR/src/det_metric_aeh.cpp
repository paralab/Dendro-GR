// Dendro: {{{ 
// Dendro: original ops: 1759 
// Dendro: printing temp variables
const double DENDRO_0 = gt0[pp]*gt3[pp];
const double DENDRO_1 = 2*grad_0_A1;
const double DENDRO_2 = DENDRO_1*grad_0_A0;
const double DENDRO_3 = DENDRO_2*grad_1_A1;
const double DENDRO_4 = DENDRO_3*grad_1_A0;
const double DENDRO_5 = gt0[pp]*gt4[pp];
const double DENDRO_6 = grad_1_A0*grad_1_A2;
const double DENDRO_7 = DENDRO_2*DENDRO_6;
const double DENDRO_8 = gt1[pp]*gt2[pp];
const double DENDRO_9 = gt1[pp]*gt4[pp];
const double DENDRO_10 = 2*DENDRO_9;
const double DENDRO_11 = DENDRO_10*grad_0_A1*grad_1_A1;
const double DENDRO_12 = gt2[pp]*gt3[pp];
const double DENDRO_13 = grad_0_A0*grad_0_A2;
const double DENDRO_14 = 2*grad_1_A1;
const double DENDRO_15 = DENDRO_13*DENDRO_14*grad_1_A0;
const double DENDRO_16 = gt0[pp]*gt5[pp];
const double DENDRO_17 = 2*DENDRO_13*DENDRO_6;
const double DENDRO_18 = gt1[pp]*gt5[pp];
const double DENDRO_19 = DENDRO_14*grad_1_A2;
const double DENDRO_20 = DENDRO_13*DENDRO_19;
const double DENDRO_21 = gt2[pp]*gt4[pp];
const double DENDRO_22 = DENDRO_1*grad_0_A2;
const double DENDRO_23 = DENDRO_22*grad_1_A1;
const double DENDRO_24 = DENDRO_21*grad_1_A0;
const double DENDRO_25 = gt3[pp]*gt5[pp];
const double DENDRO_26 = DENDRO_23*grad_1_A2;
const double DENDRO_27 = (grad_0_A0 * grad_0_A0);
const double DENDRO_28 = (grad_1_A1 * grad_1_A1);
const double DENDRO_29 = DENDRO_27*DENDRO_28;
const double DENDRO_30 = (grad_1_A2 * grad_1_A2);
const double DENDRO_31 = DENDRO_27*DENDRO_30;
const double DENDRO_32 = (grad_0_A1 * grad_0_A1);
const double DENDRO_33 = (grad_1_A0 * grad_1_A0);
const double DENDRO_34 = DENDRO_32*DENDRO_33;
const double DENDRO_35 = DENDRO_30*DENDRO_32;
const double DENDRO_36 = (grad_0_A2 * grad_0_A2);
const double DENDRO_37 = DENDRO_33*DENDRO_36;
const double DENDRO_38 = DENDRO_28*DENDRO_36;
const double DENDRO_39 = pow(gt1[pp], 2);
const double DENDRO_40 = DENDRO_2*DENDRO_30;
const double DENDRO_41 = pow(gt2[pp], 2);
const double DENDRO_42 = DENDRO_13*DENDRO_28;
const double DENDRO_43 = 2*DENDRO_12;
const double DENDRO_44 = DENDRO_19*DENDRO_27;
const double DENDRO_45 = DENDRO_22*DENDRO_33;
const double DENDRO_46 = pow(gt4[pp], 2);
const double DENDRO_47 = DENDRO_32*DENDRO_6;
const double DENDRO_48 = DENDRO_14*DENDRO_36;
const double DENDRO_49 = sqrt(fabs((DENDRO_0*DENDRO_29 + DENDRO_0*DENDRO_34 - DENDRO_0*DENDRO_4 - DENDRO_10*DENDRO_42 - DENDRO_10*DENDRO_47 + DENDRO_11*grad_0_A0*grad_1_A2 + DENDRO_11*grad_0_A2*grad_1_A0 - DENDRO_12*DENDRO_23*grad_1_A0 - DENDRO_12*DENDRO_3*grad_1_A2 - DENDRO_15*DENDRO_5 + DENDRO_15*DENDRO_8 - DENDRO_16*DENDRO_17 + DENDRO_16*DENDRO_31 + DENDRO_16*DENDRO_37 + DENDRO_17*DENDRO_41 - DENDRO_18*DENDRO_20 - DENDRO_18*DENDRO_22*DENDRO_6 + DENDRO_18*DENDRO_40 + DENDRO_18*DENDRO_48*grad_1_A0 + DENDRO_20*DENDRO_21 - DENDRO_21*DENDRO_40 + DENDRO_22*DENDRO_24*grad_1_A2 - DENDRO_24*DENDRO_48 - DENDRO_25*DENDRO_26 + DENDRO_25*DENDRO_35 + DENDRO_25*DENDRO_38 + DENDRO_26*DENDRO_46 - DENDRO_29*DENDRO_39 - DENDRO_31*DENDRO_41 - DENDRO_34*DENDRO_39 - DENDRO_35*DENDRO_46 - DENDRO_37*DENDRO_41 - DENDRO_38*DENDRO_46 + DENDRO_39*DENDRO_4 + DENDRO_42*DENDRO_43 + DENDRO_43*DENDRO_47 + DENDRO_44*DENDRO_5 - DENDRO_44*DENDRO_8 + DENDRO_45*DENDRO_5 - DENDRO_45*DENDRO_8 - DENDRO_5*DENDRO_7 + DENDRO_7*DENDRO_8)/pow(chi[pp], 2)));
const double DENDRO_50 = -DENDRO_12 + DENDRO_9;
const double DENDRO_51 = DENDRO_18 - DENDRO_21;
const double DENDRO_52 = DENDRO_50*grad_2_F - DENDRO_51*grad_1_F + grad_0_F*(DENDRO_25 - DENDRO_46);
const double DENDRO_53 = 3*At2[pp] + K[pp]*gt2[pp];
const double DENDRO_54 = DENDRO_52*DENDRO_53;
const double DENDRO_55 = DENDRO_5 - DENDRO_8;
const double DENDRO_56 = DENDRO_50*grad_0_F - DENDRO_55*grad_1_F + grad_2_F*(DENDRO_0 - DENDRO_39);
const double DENDRO_57 = DENDRO_56*(3*At5[pp] + K[pp]*gt5[pp]);
const double DENDRO_58 = 3*At3[pp] + K[pp]*gt3[pp];
const double DENDRO_59 = DENDRO_51*grad_0_F + DENDRO_55*grad_2_F - grad_1_F*(DENDRO_16 - DENDRO_41);
const double DENDRO_60 = DENDRO_59*zz;
const double DENDRO_61 = 3*At4[pp] + K[pp]*gt4[pp];
const double DENDRO_62 = DENDRO_59*DENDRO_61;
const double DENDRO_63 = 3*At1[pp] + K[pp]*gt1[pp];
const double DENDRO_64 = DENDRO_52*zz;
const double DENDRO_65 = DENDRO_56*zz;
const double DENDRO_66 = 1.0/(-DENDRO_10*gt2[pp] - DENDRO_25*gt0[pp] + DENDRO_39*gt5[pp] + DENDRO_41*gt3[pp] + DENDRO_46*gt0[pp]);
const double DENDRO_67 = (1.0/3.0)*DENDRO_49*DENDRO_66/sqrt(-DENDRO_66*chi[pp]*(DENDRO_52*grad_0_F + DENDRO_56*grad_2_F - DENDRO_59*grad_1_F));
const double DENDRO_68 = 3*At0[pp] + K[pp]*gt0[pp];

// Dendro: printing variables
//--
sqrt_det_m_ab[pp] = DENDRO_49;
//--
J_x[pp] = -DENDRO_67*(DENDRO_54*yy + DENDRO_57*yy + DENDRO_58*DENDRO_60 - DENDRO_61*DENDRO_65 - DENDRO_62*yy - DENDRO_63*DENDRO_64);
//--
J_y[pp] = DENDRO_67*(-DENDRO_53*DENDRO_65 + DENDRO_54*xx + DENDRO_57*xx + DENDRO_60*DENDRO_63 - DENDRO_62*xx - DENDRO_64*DENDRO_68);
//--
J_z[pp] = -DENDRO_67*(DENDRO_52*DENDRO_63*xx - DENDRO_52*DENDRO_68*yy - DENDRO_53*DENDRO_56*yy + DENDRO_56*DENDRO_61*xx - DENDRO_58*DENDRO_59*xx + DENDRO_59*DENDRO_63*yy);
// Dendro: reduced ops: 259
// Dendro: }}} 
