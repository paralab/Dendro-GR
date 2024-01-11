// Dendro: {{{ 
// Dendro: original ops: 283 
// Dendro: printing temp variables
const double DENDRO_0 = gt0[pp]*gt3[pp];
const double DENDRO_1 = grad_0_A1*grad_1_A0;
const double DENDRO_2 = 2*grad_0_A0;
const double DENDRO_3 = DENDRO_2*grad_1_A1;
const double DENDRO_4 = DENDRO_1*DENDRO_3;
const double DENDRO_5 = 2*gt4[pp];
const double DENDRO_6 = DENDRO_5*grad_0_A0;
const double DENDRO_7 = DENDRO_6*gt0[pp];
const double DENDRO_8 = DENDRO_1*grad_1_A2;
const double DENDRO_9 = grad_0_A1*gt2[pp];
const double DENDRO_10 = 2*DENDRO_9;
const double DENDRO_11 = grad_1_A2*gt1[pp];
const double DENDRO_12 = grad_0_A0*grad_1_A0;
const double DENDRO_13 = grad_0_A1*grad_1_A1;
const double DENDRO_14 = grad_1_A1*grad_1_A2;
const double DENDRO_15 = DENDRO_10*gt3[pp];
const double DENDRO_16 = grad_0_A2*grad_1_A0*grad_1_A1;
const double DENDRO_17 = grad_0_A2*gt1[pp];
const double DENDRO_18 = grad_1_A1*gt2[pp];
const double DENDRO_19 = 2*DENDRO_18;
const double DENDRO_20 = gt0[pp]*gt5[pp];
const double DENDRO_21 = DENDRO_2*grad_0_A2;
const double DENDRO_22 = grad_1_A0*grad_1_A2;
const double DENDRO_23 = DENDRO_21*DENDRO_22;
const double DENDRO_24 = DENDRO_17*gt5[pp];
const double DENDRO_25 = grad_0_A2*grad_1_A2;
const double DENDRO_26 = DENDRO_5*grad_1_A0;
const double DENDRO_27 = gt3[pp]*gt5[pp];
const double DENDRO_28 = grad_0_A1*grad_0_A2;
const double DENDRO_29 = 2*grad_1_A1;
const double DENDRO_30 = DENDRO_28*DENDRO_29*grad_1_A2;
const double DENDRO_31 = (grad_0_A0 * grad_0_A0);
const double DENDRO_32 = (grad_1_A1 * grad_1_A1);
const double DENDRO_33 = DENDRO_31*DENDRO_32;
const double DENDRO_34 = (grad_1_A2 * grad_1_A2);
const double DENDRO_35 = DENDRO_31*DENDRO_34;
const double DENDRO_36 = (grad_0_A1 * grad_0_A1);
const double DENDRO_37 = (grad_1_A0 * grad_1_A0);
const double DENDRO_38 = DENDRO_36*DENDRO_37;
const double DENDRO_39 = DENDRO_34*DENDRO_36;
const double DENDRO_40 = (grad_0_A2 * grad_0_A2);
const double DENDRO_41 = DENDRO_37*DENDRO_40;
const double DENDRO_42 = DENDRO_32*DENDRO_40;
const double DENDRO_43 = pow(gt1[pp], 2);
const double DENDRO_44 = gt1[pp]*gt5[pp];
const double DENDRO_45 = pow(gt2[pp], 2);
const double DENDRO_46 = gt2[pp]*gt3[pp];
const double DENDRO_47 = DENDRO_5*gt0[pp];
const double DENDRO_48 = pow(gt4[pp], 2);

// Dendro: printing variables
//--
det_m_ab[pp] = (DENDRO_0*DENDRO_33 + DENDRO_0*DENDRO_38 - DENDRO_0*DENDRO_4 + DENDRO_10*DENDRO_11*DENDRO_12 - DENDRO_10*DENDRO_17*DENDRO_37 + DENDRO_11*DENDRO_13*DENDRO_6 - DENDRO_11*DENDRO_19*DENDRO_31 - DENDRO_11*DENDRO_26*DENDRO_36 + DENDRO_12*DENDRO_17*DENDRO_19 + DENDRO_13*DENDRO_17*DENDRO_26 - DENDRO_14*DENDRO_15*grad_0_A0 + DENDRO_14*DENDRO_31*DENDRO_47 - DENDRO_15*DENDRO_16 - DENDRO_16*DENDRO_7 - DENDRO_17*DENDRO_32*DENDRO_6 + DENDRO_18*DENDRO_25*DENDRO_6 - DENDRO_18*DENDRO_26*DENDRO_40 + DENDRO_2*DENDRO_34*DENDRO_44*grad_0_A1 - DENDRO_20*DENDRO_23 + DENDRO_20*DENDRO_35 + DENDRO_20*DENDRO_41 + DENDRO_21*DENDRO_32*DENDRO_46 + 2*DENDRO_22*DENDRO_36*DENDRO_46 + DENDRO_23*DENDRO_45 - DENDRO_24*DENDRO_3*grad_1_A2 - 2*DENDRO_24*DENDRO_8 + DENDRO_25*DENDRO_26*DENDRO_9 - DENDRO_27*DENDRO_30 + DENDRO_27*DENDRO_39 + DENDRO_27*DENDRO_42 + DENDRO_28*DENDRO_37*DENDRO_47 + DENDRO_29*DENDRO_40*DENDRO_44*grad_1_A0 + DENDRO_30*DENDRO_48 - DENDRO_33*DENDRO_43 - DENDRO_34*DENDRO_6*DENDRO_9 - DENDRO_35*DENDRO_45 - DENDRO_38*DENDRO_43 - DENDRO_39*DENDRO_48 + DENDRO_4*DENDRO_43 - DENDRO_41*DENDRO_45 - DENDRO_42*DENDRO_48 - DENDRO_7*DENDRO_8)/pow(chi[pp], 2);
// Dendro: reduced ops: 160
// Dendro: }}} 
