double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;double DENDRO_33;double DENDRO_34;double DENDRO_35;double DENDRO_36;double DENDRO_37;double DENDRO_38;double DENDRO_39;double DENDRO_40;double DENDRO_41;double DENDRO_42;double DENDRO_43;double DENDRO_44;double DENDRO_45;double DENDRO_46;double DENDRO_47;double DENDRO_48;double DENDRO_49;double DENDRO_50;double DENDRO_51;double DENDRO_52;double DENDRO_53;double DENDRO_54;double DENDRO_55;double DENDRO_56;double DENDRO_57;double DENDRO_58;double DENDRO_59;double DENDRO_60;double DENDRO_61;double DENDRO_62;double DENDRO_63;double DENDRO_64;double DENDRO_65;double DENDRO_66;double DENDRO_67;double DENDRO_68;double DENDRO_69;double DENDRO_70;double DENDRO_71;double DENDRO_72;double DENDRO_73;double DENDRO_74;double DENDRO_75;double DENDRO_76;double DENDRO_77;double DENDRO_78;double DENDRO_79;double DENDRO_80;double DENDRO_81;double DENDRO_82;double DENDRO_83;double DENDRO_84;double DENDRO_85;double DENDRO_86;double DENDRO_87;double DENDRO_88;double DENDRO_89;double DENDRO_90;double DENDRO_91;double DENDRO_92;double DENDRO_93;double DENDRO_94;double DENDRO_95;double DENDRO_96;double DENDRO_97;double DENDRO_98;double DENDRO_99;
// Dendro: {{{ 
// Dendro: original ops: 114 
// Dendro: printing temp variables
DENDRO_0 = gt3*gt5;
DENDRO_1 = (gt4 * gt4);
DENDRO_2 = (gt1 * gt1);
DENDRO_3 = (gt2 * gt2);
DENDRO_4 = gt2*gt4;
DENDRO_5 = 1.0/(-DENDRO_0*gt0 + DENDRO_1*gt0 + DENDRO_2*gt5 + DENDRO_3*gt3 - 2*DENDRO_4*gt1);

// Dendro: printing variables
//--
double DENDRO_igt0 = -DENDRO_5*(DENDRO_0 - DENDRO_1);
//--
double DENDRO_igt1 = DENDRO_5*(-DENDRO_4 + gt1*gt5);
//--
double DENDRO_igt2 = -DENDRO_5*(gt1*gt4 - gt2*gt3);
//--
double DENDRO_igt3 = -DENDRO_5*(-DENDRO_3 + gt0*gt5);
//--
double DENDRO_igt4 = DENDRO_5*(gt0*gt4 - gt1*gt2);
//--
double DENDRO_igt5 = -DENDRO_5*(-DENDRO_2 + gt0*gt3);
// Dendro: reduced ops: 39
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C1_k0_0 = 0.5*grad_0_gt0;
//--
const double DENDRO_C1_k0_1 = 0.5*grad_1_gt0;
//--
const double DENDRO_C1_k0_2 = 0.5*grad_2_gt0;
//--
const double DENDRO_C1_k0_3 = -0.5*grad_0_gt3 + 1.0*grad_1_gt1;
//--
const double DENDRO_C1_k0_4 = 0.5*(-grad_0_gt4 + grad_1_gt2 + grad_2_gt1);
//--
const double DENDRO_C1_k0_5 = -0.5*grad_0_gt5 + 1.0*grad_2_gt2;
// Dendro: reduced ops: 12
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C1_k1_0 = 1.0*grad_0_gt1 - 0.5*grad_1_gt0;
//--
const double DENDRO_C1_k1_1 = 0.5*grad_0_gt3;
//--
const double DENDRO_C1_k1_2 = 0.5*(grad_0_gt4 - grad_1_gt2 + grad_2_gt1);
//--
const double DENDRO_C1_k1_3 = 0.5*grad_1_gt3;
//--
const double DENDRO_C1_k1_4 = 0.5*grad_2_gt3;
//--
const double DENDRO_C1_k1_5 = -0.5*grad_1_gt5 + 1.0*grad_2_gt4;
// Dendro: reduced ops: 12
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 14 
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C1_k2_0 = 1.0*grad_0_gt2 - 0.5*grad_2_gt0;
//--
const double DENDRO_C1_k2_1 = 0.5*(grad_0_gt4 + grad_1_gt2 - grad_2_gt1);
//--
const double DENDRO_C1_k2_2 = 0.5*grad_0_gt5;
//--
const double DENDRO_C1_k2_3 = 1.0*grad_1_gt4 - 0.5*grad_2_gt3;
//--
const double DENDRO_C1_k2_4 = 0.5*grad_1_gt5;
//--
const double DENDRO_C1_k2_5 = 0.5*grad_2_gt5;
// Dendro: reduced ops: 12
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 30 
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C2_k0_0 = DENDRO_C1_k0_0*DENDRO_igt0 + DENDRO_C1_k1_0*DENDRO_igt1 + DENDRO_C1_k2_0*DENDRO_igt2;
//--
const double DENDRO_C2_k0_1 = DENDRO_C1_k0_1*DENDRO_igt0 + DENDRO_C1_k1_1*DENDRO_igt1 + DENDRO_C1_k2_1*DENDRO_igt2;
//--
const double DENDRO_C2_k0_2 = DENDRO_C1_k0_2*DENDRO_igt0 + DENDRO_C1_k1_2*DENDRO_igt1 + DENDRO_C1_k2_2*DENDRO_igt2;
//--
const double DENDRO_C2_k0_3 = DENDRO_C1_k0_3*DENDRO_igt0 + DENDRO_C1_k1_3*DENDRO_igt1 + DENDRO_C1_k2_3*DENDRO_igt2;
//--
const double DENDRO_C2_k0_4 = DENDRO_C1_k0_4*DENDRO_igt0 + DENDRO_C1_k1_4*DENDRO_igt1 + DENDRO_C1_k2_4*DENDRO_igt2;
//--
const double DENDRO_C2_k0_5 = DENDRO_C1_k0_5*DENDRO_igt0 + DENDRO_C1_k1_5*DENDRO_igt1 + DENDRO_C1_k2_5*DENDRO_igt2;
// Dendro: reduced ops: 30
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 30 
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C2_k1_0 = DENDRO_C1_k0_0*DENDRO_igt1 + DENDRO_C1_k1_0*DENDRO_igt3 + DENDRO_C1_k2_0*DENDRO_igt4;
//--
const double DENDRO_C2_k1_1 = DENDRO_C1_k0_1*DENDRO_igt1 + DENDRO_C1_k1_1*DENDRO_igt3 + DENDRO_C1_k2_1*DENDRO_igt4;
//--
const double DENDRO_C2_k1_2 = DENDRO_C1_k0_2*DENDRO_igt1 + DENDRO_C1_k1_2*DENDRO_igt3 + DENDRO_C1_k2_2*DENDRO_igt4;
//--
const double DENDRO_C2_k1_3 = DENDRO_C1_k0_3*DENDRO_igt1 + DENDRO_C1_k1_3*DENDRO_igt3 + DENDRO_C1_k2_3*DENDRO_igt4;
//--
const double DENDRO_C2_k1_4 = DENDRO_C1_k0_4*DENDRO_igt1 + DENDRO_C1_k1_4*DENDRO_igt3 + DENDRO_C1_k2_4*DENDRO_igt4;
//--
const double DENDRO_C2_k1_5 = DENDRO_C1_k0_5*DENDRO_igt1 + DENDRO_C1_k1_5*DENDRO_igt3 + DENDRO_C1_k2_5*DENDRO_igt4;
// Dendro: reduced ops: 30
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 30 
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C2_k2_0 = DENDRO_C1_k0_0*DENDRO_igt2 + DENDRO_C1_k1_0*DENDRO_igt4 + DENDRO_C1_k2_0*DENDRO_igt5;
//--
const double DENDRO_C2_k2_1 = DENDRO_C1_k0_1*DENDRO_igt2 + DENDRO_C1_k1_1*DENDRO_igt4 + DENDRO_C1_k2_1*DENDRO_igt5;
//--
const double DENDRO_C2_k2_2 = DENDRO_C1_k0_2*DENDRO_igt2 + DENDRO_C1_k1_2*DENDRO_igt4 + DENDRO_C1_k2_2*DENDRO_igt5;
//--
const double DENDRO_C2_k2_3 = DENDRO_C1_k0_3*DENDRO_igt2 + DENDRO_C1_k1_3*DENDRO_igt4 + DENDRO_C1_k2_3*DENDRO_igt5;
//--
const double DENDRO_C2_k2_4 = DENDRO_C1_k0_4*DENDRO_igt2 + DENDRO_C1_k1_4*DENDRO_igt4 + DENDRO_C1_k2_4*DENDRO_igt5;
//--
const double DENDRO_C2_k2_5 = DENDRO_C1_k0_5*DENDRO_igt2 + DENDRO_C1_k1_5*DENDRO_igt4 + DENDRO_C1_k2_5*DENDRO_igt5;
// Dendro: {{{ 
// Dendro: original ops: 60 
// Dendro: printing temp variables
DENDRO_0 = 1.0/chi;
DENDRO_1 = 0.5*DENDRO_igt0*grad_0_chi + 0.5*DENDRO_igt1*grad_1_chi + 0.5*DENDRO_igt2*grad_2_chi;
DENDRO_2 = DENDRO_0*DENDRO_1;

// Dendro: printing variables
//--
const double DENDRO_C3_k0_0 = -DENDRO_0*(-DENDRO_1*gt0 + 1.0*grad_0_chi) + DENDRO_C2_k0_0;
//--
const double DENDRO_C3_k0_1 = -DENDRO_0*(-DENDRO_1*gt1 + 0.5*grad_1_chi) + DENDRO_C2_k0_1;
//--
const double DENDRO_C3_k0_2 = -DENDRO_0*(-DENDRO_1*gt2 + 0.5*grad_2_chi) + DENDRO_C2_k0_2;
//--
const double DENDRO_C3_k0_3 = DENDRO_2*gt3 + DENDRO_C2_k0_3;
//--
const double DENDRO_C3_k0_4 = DENDRO_2*gt4 + DENDRO_C2_k0_4;
//--
const double DENDRO_C3_k0_5 = DENDRO_2*gt5 + DENDRO_C2_k0_5;
// Dendro: reduced ops: 31
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 60 
// Dendro: printing temp variables
DENDRO_0 = 1.0/chi;
DENDRO_1 = 0.5*DENDRO_igt1*grad_0_chi + 0.5*DENDRO_igt3*grad_1_chi + 0.5*DENDRO_igt4*grad_2_chi;
DENDRO_2 = DENDRO_0*DENDRO_1;

// Dendro: printing variables
//--
const double DENDRO_C3_k1_0 = DENDRO_2*gt0 + DENDRO_C2_k1_0;
//--
const double DENDRO_C3_k1_1 = -DENDRO_0*(-DENDRO_1*gt1 + 0.5*grad_0_chi) + DENDRO_C2_k1_1;
//--
const double DENDRO_C3_k1_2 = DENDRO_2*gt2 + DENDRO_C2_k1_2;
//--
const double DENDRO_C3_k1_3 = -DENDRO_0*(-DENDRO_1*gt3 + 1.0*grad_1_chi) + DENDRO_C2_k1_3;
//--
const double DENDRO_C3_k1_4 = -DENDRO_0*(-DENDRO_1*gt4 + 0.5*grad_2_chi) + DENDRO_C2_k1_4;
//--
const double DENDRO_C3_k1_5 = DENDRO_2*gt5 + DENDRO_C2_k1_5;
// Dendro: reduced ops: 31
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 60 
// Dendro: printing temp variables
DENDRO_0 = 1.0/chi;
DENDRO_1 = 0.5*DENDRO_igt2*grad_0_chi + 0.5*DENDRO_igt4*grad_1_chi + 0.5*DENDRO_igt5*grad_2_chi;
DENDRO_2 = DENDRO_0*DENDRO_1;

// Dendro: printing variables
//--
const double DENDRO_C3_k2_0 = DENDRO_2*gt0 + DENDRO_C2_k2_0;
//--
const double DENDRO_C3_k2_1 = DENDRO_2*gt1 + DENDRO_C2_k2_1;
//--
const double DENDRO_C3_k2_2 = -DENDRO_0*(-DENDRO_1*gt2 + 0.5*grad_0_chi) + DENDRO_C2_k2_2;
//--
const double DENDRO_C3_k2_3 = DENDRO_2*gt3 + DENDRO_C2_k2_3;
//--
const double DENDRO_C3_k2_4 = -DENDRO_0*(-DENDRO_1*gt4 + 0.5*grad_1_chi) + DENDRO_C2_k2_4;
//--
const double DENDRO_C3_k2_5 = -DENDRO_0*(-DENDRO_1*gt5 + 1.0*grad_2_chi) + DENDRO_C2_k2_5;
// Dendro: reduced ops: 31
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 42 
// Dendro: printing temp variables
DENDRO_0 = 2*DENDRO_igt1;
DENDRO_1 = 2*DENDRO_igt2;
DENDRO_2 = 2*DENDRO_igt4;

// Dendro: printing variables
//--
const double DENDRO_Gtk0 = DENDRO_0*DENDRO_C2_k0_1 + DENDRO_1*DENDRO_C2_k0_2 + DENDRO_2*DENDRO_C2_k0_4 + DENDRO_C2_k0_0*DENDRO_igt0 + DENDRO_C2_k0_3*DENDRO_igt3 + DENDRO_C2_k0_5*DENDRO_igt5;
//--
const double DENDRO_Gtk1 = DENDRO_0*DENDRO_C2_k1_1 + DENDRO_1*DENDRO_C2_k1_2 + DENDRO_2*DENDRO_C2_k1_4 + DENDRO_C2_k1_0*DENDRO_igt0 + DENDRO_C2_k1_3*DENDRO_igt3 + DENDRO_C2_k1_5*DENDRO_igt5;
//--
const double DENDRO_Gtk2 = DENDRO_0*DENDRO_C2_k2_1 + DENDRO_1*DENDRO_C2_k2_2 + DENDRO_2*DENDRO_C2_k2_4 + DENDRO_C2_k2_0*DENDRO_igt0 + DENDRO_C2_k2_3*DENDRO_igt3 + DENDRO_C2_k2_5*DENDRO_igt5;
// Dendro: reduced ops: 36
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 414 
// Dendro: printing temp variables
DENDRO_0 = (1.0/3.0)*DENDRO_igt0;
DENDRO_1 = (1.0/3.0)*DENDRO_igt1;
DENDRO_2 = (1.0/3.0)*DENDRO_igt2;
DENDRO_3 = 2*DENDRO_igt4;
DENDRO_4 = (DENDRO_igt1 * DENDRO_igt1);
DENDRO_5 = (DENDRO_igt2 * DENDRO_igt2);
DENDRO_6 = At1*DENDRO_igt0;
DENDRO_7 = At2*DENDRO_igt0;
DENDRO_8 = DENDRO_igt1*DENDRO_igt2;
DENDRO_9 = At0*(DENDRO_igt0 * DENDRO_igt0) + At3*DENDRO_4 + 2*At4*DENDRO_8 + At5*DENDRO_5 + 2*DENDRO_6*DENDRO_igt1 + 2*DENDRO_7*DENDRO_igt2;
DENDRO_10 = 2*DENDRO_9;
DENDRO_11 = (DENDRO_igt4 * DENDRO_igt4);
DENDRO_12 = At1*DENDRO_igt1;
DENDRO_13 = At2*DENDRO_igt1;
DENDRO_14 = At4*DENDRO_igt3;
DENDRO_15 = 2*alpha;
DENDRO_16 = At1*DENDRO_igt2;
DENDRO_17 = At2*DENDRO_igt2;
DENDRO_18 = At0*DENDRO_igt0;
DENDRO_19 = At3*DENDRO_igt1;
DENDRO_20 = At4*DENDRO_igt1;
DENDRO_21 = At4*DENDRO_igt2;
DENDRO_22 = At5*DENDRO_igt2;
DENDRO_23 = At1*DENDRO_4 + At2*DENDRO_8 + DENDRO_18*DENDRO_igt1 + DENDRO_19*DENDRO_igt3 + DENDRO_20*DENDRO_igt4 + DENDRO_21*DENDRO_igt3 + DENDRO_22*DENDRO_igt4 + DENDRO_6*DENDRO_igt3 + DENDRO_7*DENDRO_igt4;
DENDRO_24 = At1*DENDRO_8 + At2*DENDRO_5 + DENDRO_18*DENDRO_igt2 + DENDRO_19*DENDRO_igt4 + DENDRO_20*DENDRO_igt5 + DENDRO_21*DENDRO_igt4 + DENDRO_22*DENDRO_igt5 + DENDRO_6*DENDRO_igt4 + DENDRO_7*DENDRO_igt5;
DENDRO_25 = 4*alpha;
DENDRO_26 = 9/chi;
DENDRO_27 = (1.0/3.0)*alpha;

// Dendro: printing variables
//--
Gt_rhs0[pp] = DENDRO_0*grad2_0_1_beta1 + DENDRO_0*grad2_0_2_beta2 + DENDRO_1*grad2_1_1_beta1 + DENDRO_1*grad2_1_2_beta2 + DENDRO_10*DENDRO_C2_k0_0*alpha - DENDRO_10*grad_0_alpha + DENDRO_15*DENDRO_C2_k0_3*(At0*DENDRO_4 + At3*(DENDRO_igt3 * DENDRO_igt3) + At5*DENDRO_11 + 2*DENDRO_12*DENDRO_igt3 + DENDRO_13*DENDRO_3 + DENDRO_14*DENDRO_3) + DENDRO_15*DENDRO_C2_k0_5*(At0*DENDRO_5 + At3*DENDRO_11 + At4*DENDRO_3*DENDRO_igt5 + At5*(DENDRO_igt5 * DENDRO_igt5) + DENDRO_16*DENDRO_3 + 2*DENDRO_17*DENDRO_igt5) + DENDRO_2*grad2_1_2_beta1 + DENDRO_2*grad2_2_2_beta2 + DENDRO_23*DENDRO_25*DENDRO_C2_k0_1 - 2*DENDRO_23*grad_1_alpha + DENDRO_24*DENDRO_25*DENDRO_C2_k0_2 - 2*DENDRO_24*grad_2_alpha + DENDRO_25*DENDRO_C2_k0_4*(At0*DENDRO_8 + At3*DENDRO_igt3*DENDRO_igt4 + At4*DENDRO_11 + At5*DENDRO_igt4*DENDRO_igt5 + DENDRO_12*DENDRO_igt4 + DENDRO_13*DENDRO_igt5 + DENDRO_14*DENDRO_igt5 + DENDRO_16*DENDRO_igt3 + DENDRO_17*DENDRO_igt4) - DENDRO_27*(DENDRO_23*DENDRO_26*grad_1_chi + 4*DENDRO_igt1*grad_1_K) - DENDRO_27*(DENDRO_24*DENDRO_26*grad_2_chi + 4*DENDRO_igt2*grad_2_K) - DENDRO_27*(DENDRO_26*DENDRO_9*grad_0_chi + 4*DENDRO_igt0*grad_0_K) + DENDRO_3*grad2_1_2_beta0 - DENDRO_Gtk0*grad_0_beta0 + (2.0/3.0)*DENDRO_Gtk0*(grad_0_beta0 + grad_1_beta1 + grad_2_beta2) - DENDRO_Gtk1*grad_1_beta0 - DENDRO_Gtk2*grad_2_beta0 + (4.0/3.0)*DENDRO_igt0*grad2_0_0_beta0 + (7.0/3.0)*DENDRO_igt1*grad2_0_1_beta0 + (7.0/3.0)*DENDRO_igt2*grad2_0_2_beta0 + DENDRO_igt3*grad2_1_1_beta0 + DENDRO_igt5*grad2_2_2_beta0 + beta0*grad_0_Gt0 + beta1*grad_1_Gt0 + beta2*grad_2_Gt0;

return;
// Dendro: reduced ops: 225
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 430 
// Dendro: printing temp variables
DENDRO_0 = At0*DENDRO_igt1;
DENDRO_1 = At1*DENDRO_igt3;
DENDRO_2 = At2*DENDRO_igt4;
DENDRO_3 = At2*DENDRO_igt1;
DENDRO_4 = At3*DENDRO_igt3;
DENDRO_5 = DENDRO_igt1*DENDRO_igt4;
DENDRO_6 = At4*DENDRO_igt3;
DENDRO_7 = At5*DENDRO_igt4;
DENDRO_8 = (DENDRO_igt1 * DENDRO_igt1);
DENDRO_9 = At1*DENDRO_8 + At4*DENDRO_5 + DENDRO_0*DENDRO_igt0 + DENDRO_1*DENDRO_igt0 + DENDRO_2*DENDRO_igt0 + DENDRO_3*DENDRO_igt2 + DENDRO_4*DENDRO_igt1 + DENDRO_6*DENDRO_igt2 + DENDRO_7*DENDRO_igt2;
DENDRO_10 = (DENDRO_igt4 * DENDRO_igt4);
DENDRO_11 = 2*DENDRO_igt1;
DENDRO_12 = At0*DENDRO_8 + At3*(DENDRO_igt3 * DENDRO_igt3) + At5*DENDRO_10 + DENDRO_1*DENDRO_11 + DENDRO_11*DENDRO_2 + 2*DENDRO_6*DENDRO_igt4;
DENDRO_13 = 2*DENDRO_12;
DENDRO_14 = At1*DENDRO_5 + At4*DENDRO_10 + DENDRO_0*DENDRO_igt2 + DENDRO_1*DENDRO_igt2 + DENDRO_2*DENDRO_igt2 + DENDRO_3*DENDRO_igt5 + DENDRO_4*DENDRO_igt4 + DENDRO_6*DENDRO_igt5 + DENDRO_7*DENDRO_igt5;
DENDRO_15 = 2*DENDRO_igt2;
DENDRO_16 = 9/chi;
DENDRO_17 = (1.0/3.0)*alpha;
DENDRO_18 = (1.0/3.0)*DENDRO_igt1;
DENDRO_19 = (1.0/3.0)*DENDRO_igt3;
DENDRO_20 = (1.0/3.0)*DENDRO_igt4;
DENDRO_21 = (DENDRO_igt2 * DENDRO_igt2);
DENDRO_22 = At1*DENDRO_igt0;
DENDRO_23 = At2*DENDRO_igt0;
DENDRO_24 = At4*DENDRO_igt1;
DENDRO_25 = 2*alpha;
DENDRO_26 = At4*DENDRO_igt4;
DENDRO_27 = 4*alpha;
DENDRO_28 = beta0*grad_0_Gt1 + beta1*grad_1_Gt1 + beta2*grad_2_Gt1;

// Dendro: printing variables
//--
Gt_rhs1[pp] = DENDRO_13*DENDRO_C2_k1_3*alpha - DENDRO_13*grad_1_alpha + DENDRO_14*DENDRO_27*DENDRO_C2_k1_4 - 2*DENDRO_14*grad_2_alpha + DENDRO_15*grad2_0_2_beta1 - DENDRO_17*(DENDRO_12*DENDRO_16*grad_1_chi + 4*DENDRO_igt3*grad_1_K) - DENDRO_17*(DENDRO_14*DENDRO_16*grad_2_chi + 4*DENDRO_igt4*grad_2_K) - DENDRO_17*(DENDRO_16*DENDRO_9*grad_0_chi + 4*DENDRO_igt1*grad_0_K) + DENDRO_18*grad2_0_0_beta0 + DENDRO_18*grad2_0_2_beta2 + DENDRO_19*grad2_0_1_beta0 + DENDRO_19*grad2_1_2_beta2 + DENDRO_20*grad2_0_2_beta0 + DENDRO_20*grad2_2_2_beta2 + DENDRO_25*DENDRO_C2_k1_0*(At0*(DENDRO_igt0 * DENDRO_igt0) + At3*DENDRO_8 + At5*DENDRO_21 + DENDRO_11*DENDRO_22 + DENDRO_15*DENDRO_23 + DENDRO_15*DENDRO_24) + DENDRO_25*DENDRO_C2_k1_5*(At0*DENDRO_21 + At1*DENDRO_15*DENDRO_igt4 + At2*DENDRO_15*DENDRO_igt5 + At3*DENDRO_10 + At5*(DENDRO_igt5 * DENDRO_igt5) + 2*DENDRO_26*DENDRO_igt5) + DENDRO_27*DENDRO_9*DENDRO_C2_k1_1 + DENDRO_27*DENDRO_C2_k1_2*(At0*DENDRO_igt0*DENDRO_igt2 + At1*DENDRO_igt1*DENDRO_igt2 + At2*DENDRO_21 + At3*DENDRO_5 + At5*DENDRO_igt2*DENDRO_igt5 + DENDRO_22*DENDRO_igt4 + DENDRO_23*DENDRO_igt5 + DENDRO_24*DENDRO_igt5 + DENDRO_26*DENDRO_igt2) + DENDRO_28 - 2*DENDRO_9*grad_0_alpha - DENDRO_Gtk0*grad_0_beta1 - DENDRO_Gtk1*grad_1_beta1 + (2.0/3.0)*DENDRO_Gtk1*(grad_0_beta0 + grad_1_beta1 + grad_2_beta2) - DENDRO_Gtk2*grad_2_beta1 + DENDRO_igt0*grad2_0_0_beta1 + (7.0/3.0)*DENDRO_igt1*grad2_0_1_beta1 + (4.0/3.0)*DENDRO_igt3*grad2_1_1_beta1 + (7.0/3.0)*DENDRO_igt4*grad2_1_2_beta1 + DENDRO_igt5*grad2_2_2_beta1;
//--
B_rhs1[pp] = -B1*eta - DENDRO_28*lambda[3] + Gt_rhs1[pp] + lambda[2]*(beta0*grad_0_B1 + beta1*grad_1_B1 + beta2*grad_2_B1);
// Dendro: reduced ops: 224
// Dendro: }}} 
// Dendro: {{{ 
// Dendro: original ops: 430 
// Dendro: printing temp variables
DENDRO_0 = At0*DENDRO_igt2;
DENDRO_1 = At1*DENDRO_igt4;
DENDRO_2 = At1*DENDRO_igt2;
DENDRO_3 = At2*DENDRO_igt5;
DENDRO_4 = At3*DENDRO_igt4;
DENDRO_5 = At4*DENDRO_igt5;
DENDRO_6 = DENDRO_igt2*DENDRO_igt4;
DENDRO_7 = At5*DENDRO_igt5;
DENDRO_8 = (DENDRO_igt2 * DENDRO_igt2);
DENDRO_9 = At2*DENDRO_8 + At4*DENDRO_6 + DENDRO_0*DENDRO_igt0 + DENDRO_1*DENDRO_igt0 + DENDRO_2*DENDRO_igt1 + DENDRO_3*DENDRO_igt0 + DENDRO_4*DENDRO_igt1 + DENDRO_5*DENDRO_igt1 + DENDRO_7*DENDRO_igt2;
DENDRO_10 = (DENDRO_igt4 * DENDRO_igt4);
DENDRO_11 = At2*DENDRO_6 + At4*DENDRO_10 + DENDRO_0*DENDRO_igt1 + DENDRO_1*DENDRO_igt1 + DENDRO_2*DENDRO_igt3 + DENDRO_3*DENDRO_igt1 + DENDRO_4*DENDRO_igt3 + DENDRO_5*DENDRO_igt3 + DENDRO_7*DENDRO_igt4;
DENDRO_12 = 2*DENDRO_igt2;
DENDRO_13 = At0*DENDRO_8 + At3*DENDRO_10 + At5*(DENDRO_igt5 * DENDRO_igt5) + DENDRO_1*DENDRO_12 + DENDRO_12*DENDRO_3 + 2*DENDRO_5*DENDRO_igt4;
DENDRO_14 = 2*DENDRO_13;
DENDRO_15 = 2*DENDRO_igt1;
DENDRO_16 = 9/chi;
DENDRO_17 = (1.0/3.0)*alpha;
DENDRO_18 = (1.0/3.0)*DENDRO_igt2;
DENDRO_19 = (1.0/3.0)*DENDRO_igt4;
DENDRO_20 = (1.0/3.0)*DENDRO_igt5;
DENDRO_21 = (DENDRO_igt1 * DENDRO_igt1);
DENDRO_22 = At1*DENDRO_igt0;
DENDRO_23 = At2*DENDRO_igt0;
DENDRO_24 = At4*DENDRO_igt2;
DENDRO_25 = 2*alpha;
DENDRO_26 = At4*DENDRO_igt4;
DENDRO_27 = 4*alpha;
DENDRO_28 = beta0*grad_0_Gt2 + beta1*grad_1_Gt2 + beta2*grad_2_Gt2;

// Dendro: printing variables
//--
Gt_rhs2[pp] = DENDRO_11*DENDRO_27*DENDRO_C2_k2_4 - 2*DENDRO_11*grad_1_alpha + DENDRO_14*DENDRO_C2_k2_5*alpha - DENDRO_14*grad_2_alpha + DENDRO_15*grad2_0_1_beta2 - DENDRO_17*(DENDRO_11*DENDRO_16*grad_1_chi + 4*DENDRO_igt4*grad_1_K) - DENDRO_17*(DENDRO_13*DENDRO_16*grad_2_chi + 4*DENDRO_igt5*grad_2_K) - DENDRO_17*(DENDRO_16*DENDRO_9*grad_0_chi + 4*DENDRO_igt2*grad_0_K) + DENDRO_18*grad2_0_0_beta0 + DENDRO_18*grad2_0_1_beta1 + DENDRO_19*grad2_0_1_beta0 + DENDRO_19*grad2_1_1_beta1 + DENDRO_20*grad2_0_2_beta0 + DENDRO_20*grad2_1_2_beta1 + DENDRO_25*DENDRO_C2_k2_0*(At0*(DENDRO_igt0 * DENDRO_igt0) + At3*DENDRO_21 + At5*DENDRO_8 + DENDRO_12*DENDRO_23 + DENDRO_15*DENDRO_22 + DENDRO_15*DENDRO_24) + DENDRO_25*DENDRO_C2_k2_3*(At0*DENDRO_21 + At1*DENDRO_15*DENDRO_igt3 + At2*DENDRO_15*DENDRO_igt4 + At3*(DENDRO_igt3 * DENDRO_igt3) + At5*DENDRO_10 + 2*DENDRO_26*DENDRO_igt3) + DENDRO_27*DENDRO_9*DENDRO_C2_k2_2 + DENDRO_27*DENDRO_C2_k2_1*(At0*DENDRO_igt0*DENDRO_igt1 + At1*DENDRO_21 + At2*DENDRO_igt1*DENDRO_igt2 + At3*DENDRO_igt1*DENDRO_igt3 + At5*DENDRO_6 + DENDRO_22*DENDRO_igt3 + DENDRO_23*DENDRO_igt4 + DENDRO_24*DENDRO_igt3 + DENDRO_26*DENDRO_igt1) + DENDRO_28 - 2*DENDRO_9*grad_0_alpha - DENDRO_Gtk0*grad_0_beta2 - DENDRO_Gtk1*grad_1_beta2 - DENDRO_Gtk2*grad_2_beta2 + (2.0/3.0)*DENDRO_Gtk2*(grad_0_beta0 + grad_1_beta1 + grad_2_beta2) + DENDRO_igt0*grad2_0_0_beta2 + (7.0/3.0)*DENDRO_igt2*grad2_0_2_beta2 + DENDRO_igt3*grad2_1_1_beta2 + (7.0/3.0)*DENDRO_igt4*grad2_1_2_beta2 + (4.0/3.0)*DENDRO_igt5*grad2_2_2_beta2;
//--
B_rhs2[pp] = -B2*eta - DENDRO_28*lambda[3] + Gt_rhs2[pp] + lambda[2]*(beta0*grad_0_B2 + beta1*grad_1_B2 + beta2*grad_2_B2);
// Dendro: reduced ops: 224
// Dendro: }}} 

