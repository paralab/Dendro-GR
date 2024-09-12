double DENDRO_0;
double DENDRO_1;
double DENDRO_2;
double DENDRO_3;
double DENDRO_4;
double DENDRO_5;
double DENDRO_6;
double DENDRO_7;
double DENDRO_8;
double DENDRO_9;
double DENDRO_10;
double DENDRO_11;
double DENDRO_12;
double DENDRO_13;
double DENDRO_14;
double DENDRO_15;
double DENDRO_16;
double DENDRO_17;
double DENDRO_18;
double DENDRO_19;
double DENDRO_20;
double DENDRO_21;
double DENDRO_22;
double DENDRO_23;
double DENDRO_24;
double DENDRO_25;
double DENDRO_26;
double DENDRO_27;
double DENDRO_28;
double DENDRO_29;
double DENDRO_30;
double DENDRO_31;
double DENDRO_32;
double DENDRO_33;
double DENDRO_34;
double DENDRO_35;
double DENDRO_36;
double DENDRO_37;
double DENDRO_38;
double DENDRO_39;
double DENDRO_40;
double DENDRO_41;
double DENDRO_42;
double DENDRO_43;
double DENDRO_44;
double DENDRO_45;
double DENDRO_46;
double DENDRO_47;
double DENDRO_48;
double DENDRO_49;
double DENDRO_50;
double DENDRO_51;
double DENDRO_52;
double DENDRO_53;
double DENDRO_54;
double DENDRO_55;
double DENDRO_56;
double DENDRO_57;
double DENDRO_58;
double DENDRO_59;
double DENDRO_60;
double DENDRO_61;
double DENDRO_62;
double DENDRO_63;
double DENDRO_64;
double DENDRO_65;
double DENDRO_66;
double DENDRO_67;
double DENDRO_68;
double DENDRO_69;
double DENDRO_70;
double DENDRO_71;
double DENDRO_72;
double DENDRO_73;
double DENDRO_74;
double DENDRO_75;
double DENDRO_76;
double DENDRO_77;
double DENDRO_78;
double DENDRO_79;
double DENDRO_80;
double DENDRO_81;
double DENDRO_82;
double DENDRO_83;
double DENDRO_84;
double DENDRO_85;
double DENDRO_86;
double DENDRO_87;
double DENDRO_88;
double DENDRO_89;
double DENDRO_90;
double DENDRO_91;
double DENDRO_92;
double DENDRO_93;
double DENDRO_94;
double DENDRO_95;
double DENDRO_96;
double DENDRO_97;
double DENDRO_98;
double DENDRO_99;
double DENDRO_100;
double DENDRO_101;
double DENDRO_102;
double DENDRO_103;
double DENDRO_104;
double DENDRO_105;
double DENDRO_106;
double DENDRO_107;
double DENDRO_108;
double DENDRO_109;
double DENDRO_110;
double DENDRO_111;
double DENDRO_112;
double DENDRO_113;
double DENDRO_114;
double DENDRO_115;
double DENDRO_116;
double DENDRO_117;
double DENDRO_118;
double DENDRO_119;
double DENDRO_120;
double DENDRO_121;
double DENDRO_122;
double DENDRO_123;
double DENDRO_124;
double DENDRO_125;
double DENDRO_126;
double DENDRO_127;
double DENDRO_128;
double DENDRO_129;
double DENDRO_130;
double DENDRO_131;
double DENDRO_132;
double DENDRO_133;
double DENDRO_134;
double DENDRO_135;
double DENDRO_136;
double DENDRO_137;
double DENDRO_138;
double DENDRO_139;
double DENDRO_140;
double DENDRO_141;
double DENDRO_142;
double DENDRO_143;
double DENDRO_144;
double DENDRO_145;
double DENDRO_146;
double DENDRO_147;
double DENDRO_148;
double DENDRO_149;
double DENDRO_150;
double DENDRO_151;
double DENDRO_152;
double DENDRO_153;
double DENDRO_154;
// Dendro: {{{
// Dendro: original ops: 114
// Dendro: printing temp variables
DENDRO_0 = gt3 * gt5;
DENDRO_1 = pow(gt4, 2);
DENDRO_2 = pow(gt1, 2);
DENDRO_3 = pow(gt2, 2);
DENDRO_4 = gt2 * gt4;
DENDRO_5 = 1.0 / (-DENDRO_0 * gt0 + DENDRO_1 * gt0 + DENDRO_2 * gt5 +
                  DENDRO_3 * gt3 - 2 * DENDRO_4 * gt1);

// Dendro: printing variables
//--
const double DENDRO_igt0    = -DENDRO_5 * (DENDRO_0 - DENDRO_1);
//--
const double DENDRO_igt1    = DENDRO_5 * (-DENDRO_4 + gt1 * gt5);
//--
const double DENDRO_igt2    = -DENDRO_5 * (gt1 * gt4 - gt2 * gt3);
//--
const double DENDRO_igt3    = -DENDRO_5 * (-DENDRO_3 + gt0 * gt5);
//--
const double DENDRO_igt4    = DENDRO_5 * (gt0 * gt4 - gt1 * gt2);
//--
const double DENDRO_igt5    = -DENDRO_5 * (-DENDRO_2 + gt0 * gt3);
// Dendro: reduced ops: 39
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 14
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C1_k0_0 = 0.5 * grad_0_gt0;
//--
const double DENDRO_C1_k0_1 = 0.5 * grad_1_gt0;
//--
const double DENDRO_C1_k0_2 = 0.5 * grad_2_gt0;
//--
const double DENDRO_C1_k0_3 = -0.5 * grad_0_gt3 + 1.0 * grad_1_gt1;
//--
const double DENDRO_C1_k0_4 = 0.5 * (-grad_0_gt4 + grad_1_gt2 + grad_2_gt1);
//--
const double DENDRO_C1_k0_5 = -0.5 * grad_0_gt5 + 1.0 * grad_2_gt2;
// Dendro: reduced ops: 12
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 14
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C1_k1_0 = 1.0 * grad_0_gt1 - 0.5 * grad_1_gt0;
//--
const double DENDRO_C1_k1_1 = 0.5 * grad_0_gt3;
//--
const double DENDRO_C1_k1_2 = 0.5 * (grad_0_gt4 - grad_1_gt2 + grad_2_gt1);
//--
const double DENDRO_C1_k1_3 = 0.5 * grad_1_gt3;
//--
const double DENDRO_C1_k1_4 = 0.5 * grad_2_gt3;
//--
const double DENDRO_C1_k1_5 = -0.5 * grad_1_gt5 + 1.0 * grad_2_gt4;
// Dendro: reduced ops: 12
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 14
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C1_k2_0 = 1.0 * grad_0_gt2 - 0.5 * grad_2_gt0;
//--
const double DENDRO_C1_k2_1 = 0.5 * (grad_0_gt4 + grad_1_gt2 - grad_2_gt1);
//--
const double DENDRO_C1_k2_2 = 0.5 * grad_0_gt5;
//--
const double DENDRO_C1_k2_3 = 1.0 * grad_1_gt4 - 0.5 * grad_2_gt3;
//--
const double DENDRO_C1_k2_4 = 0.5 * grad_1_gt5;
//--
const double DENDRO_C1_k2_5 = 0.5 * grad_2_gt5;
// Dendro: reduced ops: 12
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 30
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C2_k0_0 = DENDRO_C1_k0_0 * DENDRO_igt0 +
                              DENDRO_C1_k1_0 * DENDRO_igt1 +
                              DENDRO_C1_k2_0 * DENDRO_igt2;
//--
const double DENDRO_C2_k0_1 = DENDRO_C1_k0_1 * DENDRO_igt0 +
                              DENDRO_C1_k1_1 * DENDRO_igt1 +
                              DENDRO_C1_k2_1 * DENDRO_igt2;
//--
const double DENDRO_C2_k0_2 = DENDRO_C1_k0_2 * DENDRO_igt0 +
                              DENDRO_C1_k1_2 * DENDRO_igt1 +
                              DENDRO_C1_k2_2 * DENDRO_igt2;
//--
const double DENDRO_C2_k0_3 = DENDRO_C1_k0_3 * DENDRO_igt0 +
                              DENDRO_C1_k1_3 * DENDRO_igt1 +
                              DENDRO_C1_k2_3 * DENDRO_igt2;
//--
const double DENDRO_C2_k0_4 = DENDRO_C1_k0_4 * DENDRO_igt0 +
                              DENDRO_C1_k1_4 * DENDRO_igt1 +
                              DENDRO_C1_k2_4 * DENDRO_igt2;
//--
const double DENDRO_C2_k0_5 = DENDRO_C1_k0_5 * DENDRO_igt0 +
                              DENDRO_C1_k1_5 * DENDRO_igt1 +
                              DENDRO_C1_k2_5 * DENDRO_igt2;
// Dendro: reduced ops: 30
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 30
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C2_k1_0 = DENDRO_C1_k0_0 * DENDRO_igt1 +
                              DENDRO_C1_k1_0 * DENDRO_igt3 +
                              DENDRO_C1_k2_0 * DENDRO_igt4;
//--
const double DENDRO_C2_k1_1 = DENDRO_C1_k0_1 * DENDRO_igt1 +
                              DENDRO_C1_k1_1 * DENDRO_igt3 +
                              DENDRO_C1_k2_1 * DENDRO_igt4;
//--
const double DENDRO_C2_k1_2 = DENDRO_C1_k0_2 * DENDRO_igt1 +
                              DENDRO_C1_k1_2 * DENDRO_igt3 +
                              DENDRO_C1_k2_2 * DENDRO_igt4;
//--
const double DENDRO_C2_k1_3 = DENDRO_C1_k0_3 * DENDRO_igt1 +
                              DENDRO_C1_k1_3 * DENDRO_igt3 +
                              DENDRO_C1_k2_3 * DENDRO_igt4;
//--
const double DENDRO_C2_k1_4 = DENDRO_C1_k0_4 * DENDRO_igt1 +
                              DENDRO_C1_k1_4 * DENDRO_igt3 +
                              DENDRO_C1_k2_4 * DENDRO_igt4;
//--
const double DENDRO_C2_k1_5 = DENDRO_C1_k0_5 * DENDRO_igt1 +
                              DENDRO_C1_k1_5 * DENDRO_igt3 +
                              DENDRO_C1_k2_5 * DENDRO_igt4;
// Dendro: reduced ops: 30
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 30
// Dendro: printing temp variables

// Dendro: printing variables
//--
const double DENDRO_C2_k2_0 = DENDRO_C1_k0_0 * DENDRO_igt2 +
                              DENDRO_C1_k1_0 * DENDRO_igt4 +
                              DENDRO_C1_k2_0 * DENDRO_igt5;
//--
const double DENDRO_C2_k2_1 = DENDRO_C1_k0_1 * DENDRO_igt2 +
                              DENDRO_C1_k1_1 * DENDRO_igt4 +
                              DENDRO_C1_k2_1 * DENDRO_igt5;
//--
const double DENDRO_C2_k2_2 = DENDRO_C1_k0_2 * DENDRO_igt2 +
                              DENDRO_C1_k1_2 * DENDRO_igt4 +
                              DENDRO_C1_k2_2 * DENDRO_igt5;
//--
const double DENDRO_C2_k2_3 = DENDRO_C1_k0_3 * DENDRO_igt2 +
                              DENDRO_C1_k1_3 * DENDRO_igt4 +
                              DENDRO_C1_k2_3 * DENDRO_igt5;
//--
const double DENDRO_C2_k2_4 = DENDRO_C1_k0_4 * DENDRO_igt2 +
                              DENDRO_C1_k1_4 * DENDRO_igt4 +
                              DENDRO_C1_k2_4 * DENDRO_igt5;
//--
const double DENDRO_C2_k2_5 = DENDRO_C1_k0_5 * DENDRO_igt2 +
                              DENDRO_C1_k1_5 * DENDRO_igt4 +
                              DENDRO_C1_k2_5 * DENDRO_igt5;
// Dendro: reduced ops: 30
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 60
// Dendro: printing temp variables
DENDRO_0 = 1.0 / chi;
DENDRO_1 = 0.5 * DENDRO_igt0 * grad_0_chi + 0.5 * DENDRO_igt1 * grad_1_chi +
           0.5 * DENDRO_igt2 * grad_2_chi;
DENDRO_2 = DENDRO_0 * DENDRO_1;

// Dendro: printing variables
//--
const double DENDRO_C3_k0_0 =
    -DENDRO_0 * (-DENDRO_1 * gt0 + 1.0 * grad_0_chi) + DENDRO_C2_k0_0;
//--
const double DENDRO_C3_k0_1 =
    -DENDRO_0 * (-DENDRO_1 * gt1 + 0.5 * grad_1_chi) + DENDRO_C2_k0_1;
//--
const double DENDRO_C3_k0_2 =
    -DENDRO_0 * (-DENDRO_1 * gt2 + 0.5 * grad_2_chi) + DENDRO_C2_k0_2;
//--
const double DENDRO_C3_k0_3 = DENDRO_2 * gt3 + DENDRO_C2_k0_3;
//--
const double DENDRO_C3_k0_4 = DENDRO_2 * gt4 + DENDRO_C2_k0_4;
//--
const double DENDRO_C3_k0_5 = DENDRO_2 * gt5 + DENDRO_C2_k0_5;
// Dendro: reduced ops: 31
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 60
// Dendro: printing temp variables
DENDRO_0                    = 1.0 / chi;
DENDRO_1 = 0.5 * DENDRO_igt1 * grad_0_chi + 0.5 * DENDRO_igt3 * grad_1_chi +
           0.5 * DENDRO_igt4 * grad_2_chi;
DENDRO_2                    = DENDRO_0 * DENDRO_1;

// Dendro: printing variables
//--
const double DENDRO_C3_k1_0 = DENDRO_2 * gt0 + DENDRO_C2_k1_0;
//--
const double DENDRO_C3_k1_1 =
    -DENDRO_0 * (-DENDRO_1 * gt1 + 0.5 * grad_0_chi) + DENDRO_C2_k1_1;
//--
const double DENDRO_C3_k1_2 = DENDRO_2 * gt2 + DENDRO_C2_k1_2;
//--
const double DENDRO_C3_k1_3 =
    -DENDRO_0 * (-DENDRO_1 * gt3 + 1.0 * grad_1_chi) + DENDRO_C2_k1_3;
//--
const double DENDRO_C3_k1_4 =
    -DENDRO_0 * (-DENDRO_1 * gt4 + 0.5 * grad_2_chi) + DENDRO_C2_k1_4;
//--
const double DENDRO_C3_k1_5 = DENDRO_2 * gt5 + DENDRO_C2_k1_5;
// Dendro: reduced ops: 31
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 60
// Dendro: printing temp variables
DENDRO_0                    = 1.0 / chi;
DENDRO_1 = 0.5 * DENDRO_igt2 * grad_0_chi + 0.5 * DENDRO_igt4 * grad_1_chi +
           0.5 * DENDRO_igt5 * grad_2_chi;
DENDRO_2                    = DENDRO_0 * DENDRO_1;

// Dendro: printing variables
//--
const double DENDRO_C3_k2_0 = DENDRO_2 * gt0 + DENDRO_C2_k2_0;
//--
const double DENDRO_C3_k2_1 = DENDRO_2 * gt1 + DENDRO_C2_k2_1;
//--
const double DENDRO_C3_k2_2 =
    -DENDRO_0 * (-DENDRO_1 * gt2 + 0.5 * grad_0_chi) + DENDRO_C2_k2_2;
//--
const double DENDRO_C3_k2_3 = DENDRO_2 * gt3 + DENDRO_C2_k2_3;
//--
const double DENDRO_C3_k2_4 =
    -DENDRO_0 * (-DENDRO_1 * gt4 + 0.5 * grad_1_chi) + DENDRO_C2_k2_4;
//--
const double DENDRO_C3_k2_5 =
    -DENDRO_0 * (-DENDRO_1 * gt5 + 1.0 * grad_2_chi) + DENDRO_C2_k2_5;
// Dendro: reduced ops: 31
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 42
// Dendro: printing temp variables
DENDRO_0 = 2 * DENDRO_igt1;
DENDRO_1 = 2 * DENDRO_igt2;
DENDRO_2 = 2 * DENDRO_igt4;

// Dendro: printing variables
//--
const double DENDRO_Gtk0 =
    DENDRO_0 * DENDRO_C2_k0_1 + DENDRO_1 * DENDRO_C2_k0_2 +
    DENDRO_2 * DENDRO_C2_k0_4 + DENDRO_C2_k0_0 * DENDRO_igt0 +
    DENDRO_C2_k0_3 * DENDRO_igt3 + DENDRO_C2_k0_5 * DENDRO_igt5;
//--
const double DENDRO_Gtk1 =
    DENDRO_0 * DENDRO_C2_k1_1 + DENDRO_1 * DENDRO_C2_k1_2 +
    DENDRO_2 * DENDRO_C2_k1_4 + DENDRO_C2_k1_0 * DENDRO_igt0 +
    DENDRO_C2_k1_3 * DENDRO_igt3 + DENDRO_C2_k1_5 * DENDRO_igt5;
//--
const double DENDRO_Gtk2 =
    DENDRO_0 * DENDRO_C2_k2_1 + DENDRO_1 * DENDRO_C2_k2_2 +
    DENDRO_2 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
    DENDRO_C2_k2_3 * DENDRO_igt3 + DENDRO_C2_k2_5 * DENDRO_igt5;
// Dendro: reduced ops: 36
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 2181
// Dendro: printing temp variables
DENDRO_0  = 0.5 * DENDRO_igt0;
DENDRO_1  = 1.0 * DENDRO_igt1;
DENDRO_2  = 1.0 * DENDRO_igt2;
DENDRO_3  = 0.5 * DENDRO_igt3;
DENDRO_4  = 1.0 * DENDRO_igt4;
DENDRO_5  = 0.5 * DENDRO_igt5;
DENDRO_6  = DENDRO_C2_k0_0 * DENDRO_igt0;
DENDRO_7  = 3 * DENDRO_igt3;
DENDRO_8  = 3 * DENDRO_igt5;
DENDRO_9  = pow(grad_0_chi, 2);
DENDRO_10 = (1.0 / 4.0) / pow(chi, 2);
DENDRO_11 = 2 * DENDRO_C1_k0_1;
DENDRO_12 = DENDRO_C2_k1_0 * DENDRO_igt0;
DENDRO_13 = 2 * DENDRO_C1_k0_3;
DENDRO_14 = 2 * DENDRO_C1_k0_4;
DENDRO_15 = 2 * DENDRO_C1_k0_2;
DENDRO_16 = DENDRO_C2_k2_0 * DENDRO_igt0;
DENDRO_17 = 2 * DENDRO_C1_k0_5;
DENDRO_18 = DENDRO_C1_k0_0 * DENDRO_C2_k0_1;
DENDRO_19 = DENDRO_C1_k0_1 * DENDRO_C2_k0_0;
DENDRO_20 = DENDRO_C1_k1_1 * DENDRO_C2_k1_0;
DENDRO_21 = DENDRO_C1_k2_1 * DENDRO_C2_k2_0;
DENDRO_22 = DENDRO_C1_k0_0 * DENDRO_C2_k0_2;
DENDRO_23 = DENDRO_C1_k0_2 * DENDRO_C2_k0_0;
DENDRO_24 = DENDRO_C1_k1_2 * DENDRO_C2_k1_0;
DENDRO_25 = DENDRO_C1_k2_2 * DENDRO_C2_k2_0;
DENDRO_26 = DENDRO_C1_k0_1 * DENDRO_C2_k0_2;
DENDRO_27 = DENDRO_C1_k0_2 * DENDRO_C2_k0_1;
DENDRO_28 = DENDRO_C1_k1_2 * DENDRO_C2_k1_1;
DENDRO_29 = DENDRO_C1_k1_1 * DENDRO_C2_k1_2;
DENDRO_30 = DENDRO_C1_k2_2 * DENDRO_C2_k2_1;
DENDRO_31 = DENDRO_C1_k2_1 * DENDRO_C2_k2_2;
DENDRO_32 = 1.0 / chi;
DENDRO_33 = (1.0 / 2.0) * DENDRO_32;
DENDRO_34 = DENDRO_C2_k0_3 * DENDRO_igt3;
DENDRO_35 = DENDRO_C2_k0_5 * DENDRO_igt5;
DENDRO_36 = 2 * DENDRO_igt1;
DENDRO_37 = 2 * DENDRO_igt2;
DENDRO_38 = 2 * DENDRO_igt4;
DENDRO_39 = DENDRO_34 + DENDRO_35 + DENDRO_36 * DENDRO_C2_k0_1 +
            DENDRO_37 * DENDRO_C2_k0_2 + DENDRO_38 * DENDRO_C2_k0_4 + DENDRO_6;
DENDRO_40 = DENDRO_C2_k1_3 * DENDRO_igt3;
DENDRO_41 = DENDRO_C2_k1_5 * DENDRO_igt5;
DENDRO_42 = DENDRO_12 + DENDRO_36 * DENDRO_C2_k1_1 +
            DENDRO_37 * DENDRO_C2_k1_2 + DENDRO_38 * DENDRO_C2_k1_4 +
            DENDRO_40 + DENDRO_41;
DENDRO_43 = DENDRO_C2_k2_3 * DENDRO_igt3;
DENDRO_44 = DENDRO_C2_k2_5 * DENDRO_igt5;
DENDRO_45 = DENDRO_16 + DENDRO_36 * DENDRO_C2_k2_1 +
            DENDRO_37 * DENDRO_C2_k2_2 + DENDRO_38 * DENDRO_C2_k2_4 +
            DENDRO_43 + DENDRO_44;
DENDRO_46 = 3 * DENDRO_32;
DENDRO_47 = DENDRO_46 * grad_0_chi;
DENDRO_48 = pow(grad_1_chi, 2);
DENDRO_49 = grad_1_chi * grad_2_chi;
DENDRO_50 = pow(grad_2_chi, 2);
DENDRO_51 = (1.0 / 4.0) * DENDRO_32 *
            (DENDRO_36 * (-DENDRO_47 * grad_1_chi + 2 * grad2_0_1_chi) +
             DENDRO_37 * (-DENDRO_47 * grad_2_chi + 2 * grad2_0_2_chi) +
             DENDRO_38 * (-DENDRO_46 * DENDRO_49 + 2 * grad2_1_2_chi) -
             2 * DENDRO_39 * grad_0_chi - 2 * DENDRO_42 * grad_1_chi -
             2 * DENDRO_45 * grad_2_chi +
             DENDRO_igt0 * (-DENDRO_46 * DENDRO_9 + 2 * grad2_0_0_chi) +
             DENDRO_igt3 * (-DENDRO_46 * DENDRO_48 + 2 * grad2_1_1_chi) +
             DENDRO_igt5 * (-DENDRO_46 * DENDRO_50 + 2 * grad2_2_2_chi));
DENDRO_52  = 0.5 * gt1;
DENDRO_53  = 0.5 * grad_0_Gt1;
DENDRO_54  = 0.5 * grad_0_Gt2;
DENDRO_55  = 0.5 * gt0;
DENDRO_56  = 0.5 * gt2;
DENDRO_57  = DENDRO_10 * grad_0_chi;
DENDRO_58  = DENDRO_C1_k1_3 * DENDRO_C2_k1_1;
DENDRO_59  = 2 * DENDRO_58;
DENDRO_60  = DENDRO_C1_k0_4 * DENDRO_C2_k1_4;
DENDRO_61  = DENDRO_C1_k1_4 * DENDRO_C2_k1_2;
DENDRO_62  = DENDRO_C1_k1_1 * DENDRO_C2_k1_1;
DENDRO_63  = DENDRO_62 + DENDRO_C1_k1_3 * DENDRO_C2_k1_0;
DENDRO_64  = DENDRO_C1_k1_0 * DENDRO_C2_k0_2;
DENDRO_65  = DENDRO_C1_k0_0 * DENDRO_C2_k0_4 + DENDRO_C1_k0_4 * DENDRO_C2_k0_0;
DENDRO_66  = DENDRO_26 + DENDRO_27;
DENDRO_67  = DENDRO_C1_k0_1 * DENDRO_C2_k1_4;
DENDRO_68  = DENDRO_C1_k1_4 * DENDRO_C2_k1_0;
DENDRO_69  = DENDRO_29 + DENDRO_68;
DENDRO_70  = DENDRO_C1_k1_2 * DENDRO_C2_k2_2;
DENDRO_71  = DENDRO_C1_k2_4 * DENDRO_C2_k2_0;
DENDRO_72  = DENDRO_71 + DENDRO_C1_k0_2 * DENDRO_C2_k2_4;
DENDRO_73  = DENDRO_C1_k0_1 * DENDRO_C2_k0_3;
DENDRO_74  = DENDRO_C1_k0_3 * DENDRO_C2_k0_1;
DENDRO_75  = DENDRO_C1_k1_1 * DENDRO_C2_k0_1;
DENDRO_76  = DENDRO_C1_k1_4 * DENDRO_C2_k2_1;
DENDRO_77  = DENDRO_C1_k2_3 * DENDRO_C2_k2_1;
DENDRO_78  = DENDRO_C1_k1_1 * DENDRO_C2_k0_2;
DENDRO_79  = DENDRO_C1_k0_1 * DENDRO_C2_k0_4;
DENDRO_80  = DENDRO_C1_k0_4 * DENDRO_C2_k0_1;
DENDRO_81  = DENDRO_79 + DENDRO_80;
DENDRO_82  = DENDRO_C1_k0_2 * DENDRO_C2_k0_3;
DENDRO_83  = DENDRO_C1_k1_2 * DENDRO_C2_k0_1;
DENDRO_84  = DENDRO_C1_k0_3 * DENDRO_C2_k1_4;
DENDRO_85  = DENDRO_C1_k1_4 * DENDRO_C2_k1_1;
DENDRO_86  = DENDRO_C1_k1_3 * DENDRO_C2_k1_2;
DENDRO_87  = DENDRO_85 + DENDRO_86;
DENDRO_88  = DENDRO_C1_k0_4 * DENDRO_C2_k2_4;
DENDRO_89  = DENDRO_C1_k2_4 * DENDRO_C2_k2_1;
DENDRO_90  = DENDRO_89 + DENDRO_C1_k1_4 * DENDRO_C2_k2_2;
DENDRO_91  = DENDRO_C1_k1_5 * DENDRO_C2_k2_1;
DENDRO_92  = DENDRO_C1_k1_2 * DENDRO_C2_k0_2;
DENDRO_93  = DENDRO_C1_k0_2 * DENDRO_C2_k0_4;
DENDRO_94  = DENDRO_C1_k0_4 * DENDRO_C2_k0_2;
DENDRO_95  = DENDRO_93 + DENDRO_94;
DENDRO_96  = DENDRO_C1_k1_5 * DENDRO_C2_k2_2;
DENDRO_97  = DENDRO_C1_k2_4 * DENDRO_C2_k2_2;
DENDRO_98  = DENDRO_97 + DENDRO_C1_k0_5 * DENDRO_C2_k2_4;
DENDRO_99  = 0.5 * DENDRO_39;
DENDRO_100 = 0.5 * DENDRO_42;
DENDRO_101 = 0.5 * DENDRO_45;
DENDRO_102 = DENDRO_C1_k2_5 * DENDRO_C2_k2_2;
DENDRO_103 = 2 * DENDRO_102;
DENDRO_104 = DENDRO_C1_k2_0 * DENDRO_C2_k0_1;
DENDRO_105 = DENDRO_C1_k2_1 * DENDRO_C2_k1_1;
DENDRO_106 = DENDRO_C1_k2_2 * DENDRO_C2_k2_2;
DENDRO_107 = DENDRO_106 + DENDRO_C1_k2_5 * DENDRO_C2_k2_0;
DENDRO_108 = DENDRO_C1_k2_1 * DENDRO_C2_k0_1;
DENDRO_109 = DENDRO_C1_k2_3 * DENDRO_C2_k1_1;
DENDRO_110 = DENDRO_C1_k0_5 * DENDRO_C2_k0_1;
DENDRO_111 = DENDRO_C1_k2_1 * DENDRO_C2_k0_2;
DENDRO_112 = DENDRO_C1_k2_2 * DENDRO_C2_k0_1;
DENDRO_113 = DENDRO_C1_k1_5 * DENDRO_C2_k1_1;
DENDRO_114 = DENDRO_C1_k2_3 * DENDRO_C2_k1_2;
DENDRO_115 = DENDRO_61 + DENDRO_C1_k2_4 * DENDRO_C2_k1_1;
DENDRO_116 = DENDRO_C1_k2_5 * DENDRO_C2_k2_1;
DENDRO_117 = DENDRO_C1_k0_2 * DENDRO_C2_k0_5;
DENDRO_118 = DENDRO_C1_k0_5 * DENDRO_C2_k0_2;
DENDRO_119 = DENDRO_C1_k2_2 * DENDRO_C2_k0_2;
DENDRO_120 = DENDRO_C1_k1_5 * DENDRO_C2_k1_2;
DENDRO_121 = DENDRO_C1_k2_4 * DENDRO_C2_k1_2;
DENDRO_122 = 3 * DENDRO_igt0;
DENDRO_123 = 2 * DENDRO_C1_k1_0;
DENDRO_124 = 2 * DENDRO_C1_k1_1;
DENDRO_125 = 2 * DENDRO_C1_k1_2;
DENDRO_126 = 2 * DENDRO_C1_k1_4;
DENDRO_127 = 2 * DENDRO_C1_k1_5;
DENDRO_128 = DENDRO_C1_k1_1 * DENDRO_C2_k1_3;
DENDRO_129 = DENDRO_C1_k1_1 * DENDRO_C2_k1_4;
DENDRO_130 = DENDRO_C1_k0_4 * DENDRO_C2_k0_3;
DENDRO_131 = DENDRO_C1_k1_3 * DENDRO_C2_k1_4;
DENDRO_132 = DENDRO_C1_k1_4 * DENDRO_C2_k1_3;
DENDRO_133 = DENDRO_C1_k2_4 * DENDRO_C2_k2_3;
DENDRO_134 = 0.5 * gt4;
DENDRO_135 = DENDRO_C1_k2_5 * DENDRO_C2_k2_4;
DENDRO_136 = 2 * DENDRO_135;
DENDRO_137 = DENDRO_C1_k2_2 * DENDRO_C2_k2_3;
DENDRO_138 = DENDRO_C1_k2_0 * DENDRO_C2_k0_4;
DENDRO_139 = DENDRO_C1_k2_1 * DENDRO_C2_k1_4;
DENDRO_140 = DENDRO_C1_k1_2 * DENDRO_C2_k1_4;
DENDRO_141 = DENDRO_C1_k2_2 * DENDRO_C2_k2_4;
DENDRO_142 = DENDRO_116 + DENDRO_141;
DENDRO_143 = DENDRO_C1_k2_4 * DENDRO_C2_k2_4;
DENDRO_144 = DENDRO_143 + DENDRO_C1_k2_5 * DENDRO_C2_k2_3;
DENDRO_145 = DENDRO_C1_k0_5 * DENDRO_C2_k0_4;
DENDRO_146 = DENDRO_C1_k2_2 * DENDRO_C2_k0_4;
DENDRO_147 = DENDRO_C1_k1_4 * DENDRO_C2_k1_5;
DENDRO_148 = DENDRO_C1_k1_5 * DENDRO_C2_k1_4;
DENDRO_149 = DENDRO_C1_k2_4 * DENDRO_C2_k1_4;
DENDRO_150 = 2 * DENDRO_C1_k2_0;
DENDRO_151 = 2 * DENDRO_C1_k2_1;
DENDRO_152 = 2 * DENDRO_C1_k2_3;
DENDRO_153 = DENDRO_C1_k2_2 * DENDRO_C2_k2_5;
DENDRO_154 = DENDRO_C1_k2_4 * DENDRO_C2_k2_5;

// Dendro: printing variables
//--
const double DENDRO_RIJ0 =
    -DENDRO_0 * grad2_0_0_gt0 - DENDRO_1 * grad2_0_1_gt0 -
    DENDRO_10 * DENDRO_9 + DENDRO_12 * (DENDRO_11 + DENDRO_C1_k1_0) +
    DENDRO_16 * (DENDRO_15 + DENDRO_C1_k2_0) - DENDRO_2 * grad2_0_2_gt0 -
    DENDRO_3 * grad2_1_1_gt0 -
    DENDRO_33 * (DENDRO_C2_k0_0 * grad_0_chi + DENDRO_C2_k1_0 * grad_1_chi +
                 DENDRO_C2_k2_0 * grad_2_chi - grad2_0_0_chi) +
    DENDRO_39 * DENDRO_C1_k0_0 - DENDRO_4 * grad2_1_2_gt0 +
    DENDRO_42 * DENDRO_C1_k0_1 + DENDRO_45 * DENDRO_C1_k0_2 -
    DENDRO_5 * grad2_2_2_gt0 + DENDRO_51 * gt0 + 3 * DENDRO_6 * DENDRO_C1_k0_0 +
    DENDRO_7 * DENDRO_C1_k0_1 * DENDRO_C2_k0_1 +
    DENDRO_8 * DENDRO_C1_k0_2 * DENDRO_C2_k0_2 +
    DENDRO_C2_k1_1 * DENDRO_igt3 * (DENDRO_13 + DENDRO_C1_k1_1) +
    DENDRO_C2_k1_2 * DENDRO_igt5 * (DENDRO_14 + DENDRO_C1_k1_2) +
    DENDRO_C2_k2_1 * DENDRO_igt3 * (DENDRO_14 + DENDRO_C1_k2_1) +
    DENDRO_C2_k2_2 * DENDRO_igt5 * (DENDRO_17 + DENDRO_C1_k2_2) +
    DENDRO_igt1 * (DENDRO_18 + 2 * DENDRO_19) +
    DENDRO_igt1 * (2 * DENDRO_18 + DENDRO_19) +
    DENDRO_igt1 * (DENDRO_11 * DENDRO_C2_k1_1 + DENDRO_20) +
    DENDRO_igt1 *
        (DENDRO_13 * DENDRO_C2_k1_0 + DENDRO_C1_k1_0 * DENDRO_C2_k1_1) +
    DENDRO_igt1 *
        (DENDRO_14 * DENDRO_C2_k2_0 + DENDRO_C1_k2_0 * DENDRO_C2_k2_1) +
    DENDRO_igt1 * (DENDRO_15 * DENDRO_C2_k2_1 + DENDRO_21) +
    DENDRO_igt2 * (DENDRO_22 + 2 * DENDRO_23) +
    DENDRO_igt2 * (2 * DENDRO_22 + DENDRO_23) +
    DENDRO_igt2 * (DENDRO_11 * DENDRO_C2_k1_2 + DENDRO_24) +
    DENDRO_igt2 *
        (DENDRO_14 * DENDRO_C2_k1_0 + DENDRO_C1_k1_0 * DENDRO_C2_k1_2) +
    DENDRO_igt2 * (DENDRO_15 * DENDRO_C2_k2_2 + DENDRO_25) +
    DENDRO_igt2 *
        (DENDRO_17 * DENDRO_C2_k2_0 + DENDRO_C1_k2_0 * DENDRO_C2_k2_2) +
    DENDRO_igt4 * (DENDRO_26 + 2 * DENDRO_27) +
    DENDRO_igt4 * (2 * DENDRO_26 + DENDRO_27) +
    DENDRO_igt4 * (DENDRO_13 * DENDRO_C2_k1_2 + DENDRO_28) +
    DENDRO_igt4 * (DENDRO_14 * DENDRO_C2_k1_1 + DENDRO_29) +
    DENDRO_igt4 * (DENDRO_14 * DENDRO_C2_k2_2 + DENDRO_30) +
    DENDRO_igt4 * (DENDRO_17 * DENDRO_C2_k2_1 + DENDRO_31) + grad_0_Gt0 * gt0 +
    grad_0_Gt1 * gt1 + grad_0_Gt2 * gt2;
//--
const double DENDRO_RIJ1 =
    -DENDRO_0 * grad2_0_0_gt1 - DENDRO_1 * grad2_0_1_gt1 +
    DENDRO_100 * (DENDRO_C1_k0_3 + DENDRO_C1_k1_1) +
    DENDRO_101 * (DENDRO_C1_k0_4 + DENDRO_C1_k1_2) - DENDRO_2 * grad2_0_2_gt1 -
    DENDRO_3 * grad2_1_1_gt1 -
    DENDRO_33 * (DENDRO_C2_k0_1 * grad_0_chi + DENDRO_C2_k1_1 * grad_1_chi +
                 DENDRO_C2_k2_1 * grad_2_chi - grad2_0_1_chi) -
    DENDRO_4 * grad2_1_2_gt1 - DENDRO_5 * grad2_2_2_gt1 + DENDRO_51 * gt1 +
    DENDRO_52 * grad_0_Gt0 + DENDRO_52 * grad_1_Gt1 + DENDRO_53 * gt3 +
    DENDRO_54 * gt4 + DENDRO_55 * grad_1_Gt0 + DENDRO_56 * grad_1_Gt2 -
    DENDRO_57 * grad_1_chi + DENDRO_99 * (DENDRO_C1_k0_1 + DENDRO_C1_k1_0) +
    DENDRO_igt0 * (2 * DENDRO_20 + DENDRO_C1_k0_1 * DENDRO_C2_k1_1) +
    DENDRO_igt0 * (DENDRO_18 + DENDRO_19 + DENDRO_C1_k1_0 * DENDRO_C2_k0_0) +
    DENDRO_igt0 * (DENDRO_21 + DENDRO_C1_k0_2 * DENDRO_C2_k2_1 +
                   DENDRO_C1_k1_2 * DENDRO_C2_k2_0) +
    DENDRO_igt1 * (DENDRO_63 + DENDRO_C1_k0_1 * DENDRO_C2_k1_3) +
    DENDRO_igt1 * (DENDRO_63 + DENDRO_C1_k0_3 * DENDRO_C2_k1_1) +
    DENDRO_igt1 *
        (DENDRO_11 * DENDRO_C2_k0_1 + DENDRO_C1_k1_1 * DENDRO_C2_k0_0) +
    DENDRO_igt1 *
        (DENDRO_C1_k0_0 * DENDRO_C2_k0_3 + DENDRO_C1_k0_3 * DENDRO_C2_k0_0 +
         DENDRO_C1_k1_0 * DENDRO_C2_k0_1) +
    DENDRO_igt1 *
        (DENDRO_C1_k0_2 * DENDRO_C2_k2_3 + DENDRO_C1_k1_2 * DENDRO_C2_k2_1 +
         DENDRO_C1_k2_3 * DENDRO_C2_k2_0) +
    DENDRO_igt1 *
        (DENDRO_C1_k0_4 * DENDRO_C2_k2_1 + DENDRO_C1_k1_4 * DENDRO_C2_k2_0 +
         DENDRO_C1_k2_1 * DENDRO_C2_k2_1) +
    DENDRO_igt2 * (DENDRO_64 + DENDRO_65) +
    DENDRO_igt2 * (DENDRO_66 + DENDRO_C1_k1_2 * DENDRO_C2_k0_0) +
    DENDRO_igt2 * (DENDRO_67 + DENDRO_69) +
    DENDRO_igt2 * (DENDRO_69 + DENDRO_C1_k0_4 * DENDRO_C2_k1_1) +
    DENDRO_igt2 * (DENDRO_70 + DENDRO_72) +
    DENDRO_igt2 * (DENDRO_31 + DENDRO_C1_k0_5 * DENDRO_C2_k2_1 +
                   DENDRO_C1_k1_5 * DENDRO_C2_k2_0) +
    DENDRO_igt3 * (DENDRO_59 + DENDRO_C1_k0_3 * DENDRO_C2_k1_3) +
    DENDRO_igt3 * (DENDRO_73 + DENDRO_74 + DENDRO_75) +
    DENDRO_igt3 * (DENDRO_76 + DENDRO_77 + DENDRO_C1_k0_4 * DENDRO_C2_k2_3) +
    DENDRO_igt4 * (DENDRO_78 + DENDRO_81) +
    DENDRO_igt4 * (DENDRO_84 + DENDRO_87) +
    DENDRO_igt4 * (DENDRO_87 + DENDRO_C1_k0_4 * DENDRO_C2_k1_3) +
    DENDRO_igt4 * (DENDRO_88 + DENDRO_90) +
    DENDRO_igt4 * (DENDRO_82 + DENDRO_83 + DENDRO_C1_k0_3 * DENDRO_C2_k0_2) +
    DENDRO_igt4 * (DENDRO_91 + DENDRO_C1_k0_5 * DENDRO_C2_k2_3 +
                   DENDRO_C1_k2_3 * DENDRO_C2_k2_2) +
    DENDRO_igt5 * (DENDRO_60 + 2 * DENDRO_61) +
    DENDRO_igt5 * (DENDRO_92 + DENDRO_95) +
    DENDRO_igt5 * (DENDRO_96 + DENDRO_98);
//--
const double DENDRO_RIJ2 =
    -DENDRO_0 * grad2_0_0_gt2 - DENDRO_1 * grad2_0_1_gt2 +
    DENDRO_100 * (DENDRO_C1_k0_4 + DENDRO_C1_k2_1) +
    DENDRO_101 * (DENDRO_C1_k0_5 + DENDRO_C1_k2_2) - DENDRO_2 * grad2_0_2_gt2 -
    DENDRO_3 * grad2_1_1_gt2 -
    DENDRO_33 * (DENDRO_C2_k0_2 * grad_0_chi + DENDRO_C2_k1_2 * grad_1_chi +
                 DENDRO_C2_k2_2 * grad_2_chi - grad2_0_2_chi) -
    DENDRO_4 * grad2_1_2_gt2 - DENDRO_5 * grad2_2_2_gt2 + DENDRO_51 * gt2 +
    DENDRO_52 * grad_2_Gt1 + DENDRO_53 * gt4 + DENDRO_54 * gt5 +
    DENDRO_55 * grad_2_Gt0 + DENDRO_56 * grad_0_Gt0 + DENDRO_56 * grad_2_Gt2 -
    DENDRO_57 * grad_2_chi + DENDRO_99 * (DENDRO_C1_k0_2 + DENDRO_C1_k2_0) +
    DENDRO_igt0 * (2 * DENDRO_25 + DENDRO_C1_k0_2 * DENDRO_C2_k2_2) +
    DENDRO_igt0 * (DENDRO_22 + DENDRO_23 + DENDRO_C1_k2_0 * DENDRO_C2_k0_0) +
    DENDRO_igt0 * (DENDRO_24 + DENDRO_C1_k0_1 * DENDRO_C2_k1_2 +
                   DENDRO_C1_k2_1 * DENDRO_C2_k1_0) +
    DENDRO_igt1 * (DENDRO_104 + DENDRO_65) +
    DENDRO_igt1 * (DENDRO_30 + DENDRO_72) +
    DENDRO_igt1 * (DENDRO_66 + DENDRO_C1_k2_1 * DENDRO_C2_k0_0) +
    DENDRO_igt1 * (DENDRO_105 + DENDRO_67 + DENDRO_68) +
    DENDRO_igt1 * (DENDRO_28 + DENDRO_C1_k0_3 * DENDRO_C2_k1_2 +
                   DENDRO_C1_k2_3 * DENDRO_C2_k1_0) +
    DENDRO_igt1 * (DENDRO_30 + DENDRO_71 + DENDRO_C1_k0_4 * DENDRO_C2_k2_2) +
    DENDRO_igt2 * (DENDRO_107 + DENDRO_C1_k0_2 * DENDRO_C2_k2_5) +
    DENDRO_igt2 * (DENDRO_107 + DENDRO_C1_k0_5 * DENDRO_C2_k2_2) +
    DENDRO_igt2 *
        (DENDRO_15 * DENDRO_C2_k0_2 + DENDRO_C1_k2_2 * DENDRO_C2_k0_0) +
    DENDRO_igt2 *
        (DENDRO_C1_k0_0 * DENDRO_C2_k0_5 + DENDRO_C1_k0_5 * DENDRO_C2_k0_0 +
         DENDRO_C1_k2_0 * DENDRO_C2_k0_2) +
    DENDRO_igt2 *
        (DENDRO_C1_k0_1 * DENDRO_C2_k1_5 + DENDRO_C1_k1_5 * DENDRO_C2_k1_0 +
         DENDRO_C1_k2_1 * DENDRO_C2_k1_2) +
    DENDRO_igt2 *
        (DENDRO_C1_k0_4 * DENDRO_C2_k1_2 + DENDRO_C1_k1_2 * DENDRO_C2_k1_2 +
         DENDRO_C1_k2_4 * DENDRO_C2_k1_0) +
    DENDRO_igt3 * (DENDRO_108 + DENDRO_81) +
    DENDRO_igt3 * (DENDRO_88 + 2 * DENDRO_89) +
    DENDRO_igt3 * (DENDRO_109 + DENDRO_84 + DENDRO_85) +
    DENDRO_igt4 * (DENDRO_112 + DENDRO_95) +
    DENDRO_igt4 * (DENDRO_115 + DENDRO_60) +
    DENDRO_igt4 * (DENDRO_116 + DENDRO_98) +
    DENDRO_igt4 * (DENDRO_110 + DENDRO_111 + DENDRO_C1_k0_1 * DENDRO_C2_k0_5) +
    DENDRO_igt4 * (DENDRO_113 + DENDRO_114 + DENDRO_C1_k0_3 * DENDRO_C2_k1_5) +
    DENDRO_igt4 * (DENDRO_116 + DENDRO_97 + DENDRO_C1_k0_4 * DENDRO_C2_k2_5) +
    DENDRO_igt5 * (DENDRO_103 + DENDRO_C1_k0_5 * DENDRO_C2_k2_5) +
    DENDRO_igt5 * (DENDRO_117 + DENDRO_118 + DENDRO_119) +
    DENDRO_igt5 * (DENDRO_120 + DENDRO_121 + DENDRO_C1_k0_4 * DENDRO_C2_k1_5);
//--
const double DENDRO_RIJ3 =
    -DENDRO_0 * grad2_0_0_gt3 - DENDRO_1 * grad2_0_1_gt3 -
    DENDRO_10 * DENDRO_48 + DENDRO_122 * DENDRO_62 - DENDRO_2 * grad2_0_2_gt3 -
    DENDRO_3 * grad2_1_1_gt3 -
    DENDRO_33 * (DENDRO_C2_k0_3 * grad_0_chi + DENDRO_C2_k1_3 * grad_1_chi +
                 DENDRO_C2_k2_3 * grad_2_chi - grad2_1_1_chi) +
    DENDRO_34 * (DENDRO_124 + DENDRO_C1_k0_3) + DENDRO_39 * DENDRO_C1_k1_1 -
    DENDRO_4 * grad2_1_2_gt3 + 3 * DENDRO_40 * DENDRO_C1_k1_3 +
    DENDRO_42 * DENDRO_C1_k1_3 + DENDRO_43 * (DENDRO_126 + DENDRO_C1_k2_3) +
    DENDRO_45 * DENDRO_C1_k1_4 - DENDRO_5 * grad2_2_2_gt3 + DENDRO_51 * gt3 +
    DENDRO_8 * DENDRO_C1_k1_4 * DENDRO_C2_k1_4 +
    DENDRO_C2_k0_1 * DENDRO_igt0 * (DENDRO_123 + DENDRO_C1_k0_1) +
    DENDRO_C2_k0_4 * DENDRO_igt5 * (DENDRO_125 + DENDRO_C1_k0_4) +
    DENDRO_C2_k2_1 * DENDRO_igt0 * (DENDRO_125 + DENDRO_C1_k2_1) +
    DENDRO_C2_k2_4 * DENDRO_igt5 * (DENDRO_127 + DENDRO_C1_k2_4) +
    DENDRO_igt1 * (DENDRO_128 + DENDRO_59) +
    DENDRO_igt1 * (2 * DENDRO_128 + DENDRO_58) +
    DENDRO_igt1 * (DENDRO_73 + 2 * DENDRO_75) +
    DENDRO_igt1 * (2 * DENDRO_76 + DENDRO_C1_k2_1 * DENDRO_C2_k2_3) +
    DENDRO_igt1 * (DENDRO_123 * DENDRO_C2_k0_3 + DENDRO_74) +
    DENDRO_igt1 * (DENDRO_125 * DENDRO_C2_k2_3 + DENDRO_77) +
    DENDRO_igt2 * (DENDRO_129 + 2 * DENDRO_85) +
    DENDRO_igt2 * (2 * DENDRO_129 + DENDRO_85) +
    DENDRO_igt2 * (DENDRO_79 + 2 * DENDRO_83) +
    DENDRO_igt2 * (2 * DENDRO_91 + DENDRO_C1_k2_1 * DENDRO_C2_k2_4) +
    DENDRO_igt2 * (DENDRO_123 * DENDRO_C2_k0_4 + DENDRO_80) +
    DENDRO_igt2 * (DENDRO_125 * DENDRO_C2_k2_4 + DENDRO_89) +
    DENDRO_igt4 * (DENDRO_131 + 2 * DENDRO_132) +
    DENDRO_igt4 * (2 * DENDRO_131 + DENDRO_132) +
    DENDRO_igt4 * (DENDRO_124 * DENDRO_C2_k0_4 + DENDRO_130) +
    DENDRO_igt4 *
        (DENDRO_125 * DENDRO_C2_k0_3 + DENDRO_C1_k0_3 * DENDRO_C2_k0_4) +
    DENDRO_igt4 * (DENDRO_126 * DENDRO_C2_k2_4 + DENDRO_133) +
    DENDRO_igt4 *
        (DENDRO_127 * DENDRO_C2_k2_3 + DENDRO_C1_k2_3 * DENDRO_C2_k2_4) +
    grad_1_Gt0 * gt1 + grad_1_Gt1 * gt3 + grad_1_Gt2 * gt4;
//--
const double DENDRO_RIJ4 =
    -DENDRO_0 * grad2_0_0_gt4 - DENDRO_1 * grad2_0_1_gt4 -
    DENDRO_10 * DENDRO_49 + DENDRO_100 * (DENDRO_C1_k1_4 + DENDRO_C1_k2_3) +
    DENDRO_101 * (DENDRO_C1_k1_5 + DENDRO_C1_k2_4) + DENDRO_134 * grad_1_Gt1 +
    DENDRO_134 * grad_2_Gt2 - DENDRO_2 * grad2_0_2_gt4 -
    DENDRO_3 * grad2_1_1_gt4 -
    DENDRO_33 * (DENDRO_C2_k0_4 * grad_0_chi + DENDRO_C2_k1_4 * grad_1_chi +
                 DENDRO_C2_k2_4 * grad_2_chi - grad2_1_2_chi) -
    DENDRO_4 * grad2_1_2_gt4 - DENDRO_5 * grad2_2_2_gt4 + DENDRO_51 * gt4 +
    DENDRO_52 * grad_2_Gt0 + DENDRO_56 * grad_1_Gt0 +
    DENDRO_99 * (DENDRO_C1_k1_2 + DENDRO_C1_k2_1) +
    DENDRO_igt0 * (2 * DENDRO_30 + DENDRO_70) +
    DENDRO_igt0 * (DENDRO_104 + DENDRO_27 + DENDRO_64) +
    DENDRO_igt0 * (DENDRO_105 + DENDRO_28 + DENDRO_29) +
    DENDRO_igt1 * (DENDRO_137 + DENDRO_90) +
    DENDRO_igt1 * (DENDRO_108 + DENDRO_78 + DENDRO_82) +
    DENDRO_igt1 * (DENDRO_109 + DENDRO_86 + DENDRO_C1_k1_2 * DENDRO_C2_k1_3) +
    DENDRO_igt1 * (DENDRO_129 + DENDRO_85 + DENDRO_C1_k2_1 * DENDRO_C2_k1_3) +
    DENDRO_igt1 * (DENDRO_137 + DENDRO_89 + DENDRO_C1_k1_2 * DENDRO_C2_k2_4) +
    DENDRO_igt1 * (DENDRO_80 + DENDRO_C1_k1_0 * DENDRO_C2_k0_4 +
                   DENDRO_C1_k2_0 * DENDRO_C2_k0_3) +
    DENDRO_igt2 * (DENDRO_115 + DENDRO_140) +
    DENDRO_igt2 * (DENDRO_142 + DENDRO_96) +
    DENDRO_igt2 * (DENDRO_142 + DENDRO_C1_k1_2 * DENDRO_C2_k2_5) +
    DENDRO_igt2 * (DENDRO_110 + DENDRO_138 + DENDRO_C1_k1_0 * DENDRO_C2_k0_5) +
    DENDRO_igt2 * (DENDRO_112 + DENDRO_92 + DENDRO_93) +
    DENDRO_igt2 * (DENDRO_113 + DENDRO_139 + DENDRO_C1_k1_1 * DENDRO_C2_k1_5) +
    DENDRO_igt3 * (2 * DENDRO_133 + DENDRO_C1_k1_4 * DENDRO_C2_k2_4) +
    DENDRO_igt3 * (DENDRO_130 + DENDRO_C1_k1_1 * DENDRO_C2_k0_4 +
                   DENDRO_C1_k2_1 * DENDRO_C2_k0_3) +
    DENDRO_igt3 * (DENDRO_131 + DENDRO_132 + DENDRO_C1_k2_3 * DENDRO_C2_k1_3) +
    DENDRO_igt4 * (DENDRO_144 + DENDRO_C1_k1_4 * DENDRO_C2_k2_5) +
    DENDRO_igt4 * (DENDRO_144 + DENDRO_C1_k1_5 * DENDRO_C2_k2_4) +
    DENDRO_igt4 *
        (DENDRO_126 * DENDRO_C2_k1_4 + DENDRO_C1_k2_4 * DENDRO_C2_k1_3) +
    DENDRO_igt4 *
        (DENDRO_C1_k0_4 * DENDRO_C2_k0_4 + DENDRO_C1_k1_2 * DENDRO_C2_k0_4 +
         DENDRO_C1_k2_2 * DENDRO_C2_k0_3) +
    DENDRO_igt4 *
        (DENDRO_C1_k0_5 * DENDRO_C2_k0_3 + DENDRO_C1_k1_1 * DENDRO_C2_k0_5 +
         DENDRO_C1_k2_1 * DENDRO_C2_k0_4) +
    DENDRO_igt4 *
        (DENDRO_C1_k1_3 * DENDRO_C2_k1_5 + DENDRO_C1_k1_5 * DENDRO_C2_k1_3 +
         DENDRO_C1_k2_3 * DENDRO_C2_k1_4) +
    DENDRO_igt5 * (DENDRO_136 + DENDRO_C1_k1_5 * DENDRO_C2_k2_5) +
    DENDRO_igt5 * (DENDRO_145 + DENDRO_146 + DENDRO_C1_k1_2 * DENDRO_C2_k0_5) +
    DENDRO_igt5 * (DENDRO_147 + DENDRO_148 + DENDRO_149) +
    0.5 * grad_1_Gt2 * gt5 + 0.5 * grad_2_Gt1 * gt3;
//--
const double DENDRO_RIJ5 =
    -DENDRO_0 * grad2_0_0_gt5 - DENDRO_1 * grad2_0_1_gt5 -
    DENDRO_10 * DENDRO_50 + DENDRO_106 * DENDRO_122 + DENDRO_143 * DENDRO_7 -
    DENDRO_2 * grad2_0_2_gt5 - DENDRO_3 * grad2_1_1_gt5 -
    DENDRO_33 * (DENDRO_C2_k0_5 * grad_0_chi + DENDRO_C2_k1_5 * grad_1_chi +
                 DENDRO_C2_k2_5 * grad_2_chi - grad2_2_2_chi) +
    DENDRO_35 * (DENDRO_C1_k0_5 + 2 * DENDRO_C1_k2_2) +
    DENDRO_39 * DENDRO_C1_k2_2 - DENDRO_4 * grad2_1_2_gt5 +
    DENDRO_41 * (DENDRO_C1_k1_5 + 2 * DENDRO_C1_k2_4) +
    DENDRO_42 * DENDRO_C1_k2_4 + 3 * DENDRO_44 * DENDRO_C1_k2_5 +
    DENDRO_45 * DENDRO_C1_k2_5 - DENDRO_5 * grad2_2_2_gt5 + DENDRO_51 * gt5 +
    DENDRO_C2_k0_2 * DENDRO_igt0 * (DENDRO_150 + DENDRO_C1_k0_2) +
    DENDRO_C2_k0_4 * DENDRO_igt3 * (DENDRO_151 + DENDRO_C1_k0_4) +
    DENDRO_C2_k1_2 * DENDRO_igt0 * (DENDRO_151 + DENDRO_C1_k1_2) +
    DENDRO_C2_k1_4 * DENDRO_igt3 * (DENDRO_152 + DENDRO_C1_k1_4) +
    DENDRO_igt1 * (2 * DENDRO_111 + DENDRO_93) +
    DENDRO_igt1 * (2 * DENDRO_114 + DENDRO_140) +
    DENDRO_igt1 * (2 * DENDRO_138 + DENDRO_94) +
    DENDRO_igt1 * (2 * DENDRO_139 + DENDRO_61) +
    DENDRO_igt1 * (DENDRO_141 + 2 * DENDRO_97) +
    DENDRO_igt1 * (2 * DENDRO_141 + DENDRO_97) +
    DENDRO_igt2 * (DENDRO_102 + 2 * DENDRO_153) +
    DENDRO_igt2 * (DENDRO_103 + DENDRO_153) +
    DENDRO_igt2 * (DENDRO_117 + 2 * DENDRO_119) +
    DENDRO_igt2 * (DENDRO_118 + DENDRO_150 * DENDRO_C2_k0_5) +
    DENDRO_igt2 * (DENDRO_120 + DENDRO_151 * DENDRO_C2_k1_5) +
    DENDRO_igt2 * (2 * DENDRO_121 + DENDRO_C1_k1_2 * DENDRO_C2_k1_5) +
    DENDRO_igt4 * (DENDRO_135 + 2 * DENDRO_154) +
    DENDRO_igt4 * (DENDRO_136 + DENDRO_154) +
    DENDRO_igt4 * (DENDRO_145 + DENDRO_151 * DENDRO_C2_k0_5) +
    DENDRO_igt4 * (2 * DENDRO_146 + DENDRO_C1_k0_4 * DENDRO_C2_k0_5) +
    DENDRO_igt4 * (DENDRO_147 + 2 * DENDRO_149) +
    DENDRO_igt4 * (DENDRO_148 + DENDRO_152 * DENDRO_C2_k1_5) +
    grad_2_Gt0 * gt2 + grad_2_Gt1 * gt4 + grad_2_Gt2 * gt5;
// Dendro: reduced ops: 1150
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 9
// Dendro: printing temp variables

// Dendro: printing variables
//--
a_rhs[pp] =
    -2 * K * alpha + lambda[0] * (beta0 * grad_0_alpha + beta1 * grad_1_alpha +
                                  beta2 * grad_2_alpha);
// Dendro: reduced ops: 9
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 42
// Dendro: printing temp variables
DENDRO_0 = (3.0 / 4.0) * alpha * lambda_f[1] + (3.0 / 4.0) * lambda_f[0];

// Dendro: printing variables
//--
b_rhs0[pp] =
    B0 * DENDRO_0 + lambda[1] * (beta0 * grad_0_beta0 + beta1 * grad_1_beta0 +
                                 beta2 * grad_2_beta0);
//--
b_rhs1[pp] =
    B1 * DENDRO_0 + lambda[1] * (beta0 * grad_0_beta1 + beta1 * grad_1_beta1 +
                                 beta2 * grad_2_beta1);
//--
b_rhs2[pp] =
    B2 * DENDRO_0 + lambda[1] * (beta0 * grad_0_beta2 + beta1 * grad_1_beta2 +
                                 beta2 * grad_2_beta2);
// Dendro: reduced ops: 30
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 156
// Dendro: printing temp variables
DENDRO_0  = 2 * alpha;
DENDRO_1  = 2 * gt1;
DENDRO_2  = 2 * gt2;
DENDRO_3  = (2.0 / 3.0) * gt0;
DENDRO_4  = (1.0 / 3.0) * gt1;
DENDRO_5  = (2.0 / 3.0) * grad_2_beta2;
DENDRO_6  = (1.0 / 3.0) * gt2;
DENDRO_7  = (2.0 / 3.0) * grad_1_beta1;
DENDRO_8  = (2.0 / 3.0) * grad_0_beta0;
DENDRO_9  = 2 * gt4;
DENDRO_10 = (1.0 / 3.0) * gt4;

// Dendro: printing variables
//--
gt_rhs00[pp] =
    -At0 * DENDRO_0 + DENDRO_1 * grad_0_beta1 + DENDRO_2 * grad_0_beta2 -
    DENDRO_3 * grad_1_beta1 - DENDRO_3 * grad_2_beta2 + beta0 * grad_0_gt0 +
    beta1 * grad_1_gt0 + beta2 * grad_2_gt0 + (4.0 / 3.0) * grad_0_beta0 * gt0;
//--
gt_rhs01[pp] = -At1 * DENDRO_0 + DENDRO_4 * grad_0_beta0 +
               DENDRO_4 * grad_1_beta1 - DENDRO_5 * gt1 + beta0 * grad_0_gt1 +
               beta1 * grad_1_gt1 + beta2 * grad_2_gt1 + grad_0_beta1 * gt3 +
               grad_0_beta2 * gt4 + grad_1_beta0 * gt0 + grad_1_beta2 * gt2;
//--
gt_rhs02[pp] = -At2 * DENDRO_0 + DENDRO_6 * grad_0_beta0 +
               DENDRO_6 * grad_2_beta2 - DENDRO_7 * gt2 + beta0 * grad_0_gt2 +
               beta1 * grad_1_gt2 + beta2 * grad_2_gt2 + grad_0_beta1 * gt4 +
               grad_0_beta2 * gt5 + grad_2_beta0 * gt0 + grad_2_beta1 * gt1;
//--
gt_rhs11[pp] = -At3 * DENDRO_0 + DENDRO_1 * grad_1_beta0 - DENDRO_5 * gt3 -
               DENDRO_8 * gt3 + DENDRO_9 * grad_1_beta2 + beta0 * grad_0_gt3 +
               beta1 * grad_1_gt3 + beta2 * grad_2_gt3 +
               (4.0 / 3.0) * grad_1_beta1 * gt3;
//--
gt_rhs12[pp] = -At4 * DENDRO_0 + DENDRO_10 * grad_1_beta1 +
               DENDRO_10 * grad_2_beta2 - DENDRO_8 * gt4 + beta0 * grad_0_gt4 +
               beta1 * grad_1_gt4 + beta2 * grad_2_gt4 + grad_1_beta0 * gt2 +
               grad_1_beta2 * gt5 + grad_2_beta0 * gt1 + grad_2_beta1 * gt3;
//--
gt_rhs22[pp] = -At5 * DENDRO_0 + DENDRO_2 * grad_2_beta0 - DENDRO_7 * gt5 -
               DENDRO_8 * gt5 + DENDRO_9 * grad_2_beta1 + beta0 * grad_0_gt5 +
               beta1 * grad_1_gt5 + beta2 * grad_2_gt5 +
               (4.0 / 3.0) * grad_2_beta2 * gt5;
// Dendro: reduced ops: 135
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 16
// Dendro: printing temp variables
DENDRO_0    = (2.0 / 3.0) * chi;

// Dendro: printing variables
//--
chi_rhs[pp] = DENDRO_0 * K * alpha -
              DENDRO_0 * (grad_0_beta0 + grad_1_beta1 + grad_2_beta2) +
              beta0 * grad_0_chi + beta1 * grad_1_chi + beta2 * grad_2_chi;
// Dendro: reduced ops: 14
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 750
// Dendro: printing temp variables
DENDRO_0  = (2.0 / 3.0) * At0;
DENDRO_1  = 2 * At1;
DENDRO_2  = 2 * At2;
DENDRO_3  = At1 * DENDRO_igt1;
DENDRO_4  = At2 * DENDRO_igt2;
DENDRO_5  = 2 * At0;
DENDRO_6  = DENDRO_C3_k0_0 * grad_0_alpha;
DENDRO_7  = DENDRO_C3_k1_0 * grad_1_alpha;
DENDRO_8  = DENDRO_C3_k2_0 * grad_2_alpha;
DENDRO_9  = DENDRO_RIJ0 * alpha;
DENDRO_10 = DENDRO_C3_k0_3 * grad_0_alpha;
DENDRO_11 = DENDRO_C3_k1_3 * grad_1_alpha;
DENDRO_12 = DENDRO_C3_k2_3 * grad_2_alpha;
DENDRO_13 = DENDRO_RIJ3 * alpha;
DENDRO_14 = DENDRO_C3_k0_5 * grad_0_alpha;
DENDRO_15 = DENDRO_C3_k1_5 * grad_1_alpha;
DENDRO_16 = DENDRO_C3_k2_5 * grad_2_alpha;
DENDRO_17 = DENDRO_RIJ5 * alpha;
DENDRO_18 = DENDRO_C3_k0_1 * grad_0_alpha;
DENDRO_19 = DENDRO_C3_k1_1 * grad_1_alpha;
DENDRO_20 = DENDRO_C3_k2_1 * grad_2_alpha;
DENDRO_21 = DENDRO_RIJ1 * alpha;
DENDRO_22 = DENDRO_C3_k0_2 * grad_0_alpha;
DENDRO_23 = DENDRO_C3_k1_2 * grad_1_alpha;
DENDRO_24 = DENDRO_C3_k2_2 * grad_2_alpha;
DENDRO_25 = DENDRO_RIJ2 * alpha;
DENDRO_26 = DENDRO_C3_k0_4 * grad_0_alpha;
DENDRO_27 = DENDRO_C3_k1_4 * grad_1_alpha;
DENDRO_28 = DENDRO_C3_k2_4 * grad_2_alpha;
DENDRO_29 = DENDRO_RIJ4 * alpha;
DENDRO_30 =
    DENDRO_igt0 *
        (DENDRO_6 + DENDRO_7 + DENDRO_8 + DENDRO_9 - grad2_0_0_alpha) +
    2 * DENDRO_igt1 *
        (DENDRO_18 + DENDRO_19 + DENDRO_20 + DENDRO_21 - grad2_0_1_alpha) +
    2 * DENDRO_igt2 *
        (DENDRO_22 + DENDRO_23 + DENDRO_24 + DENDRO_25 - grad2_0_2_alpha) +
    DENDRO_igt3 *
        (DENDRO_10 + DENDRO_11 + DENDRO_12 + DENDRO_13 - grad2_1_1_alpha) +
    2 * DENDRO_igt4 *
        (DENDRO_26 + DENDRO_27 + DENDRO_28 + DENDRO_29 - grad2_1_2_alpha) +
    DENDRO_igt5 *
        (DENDRO_14 + DENDRO_15 + DENDRO_16 + DENDRO_17 - grad2_2_2_alpha);
DENDRO_31 = (1.0 / 3.0) * chi;
DENDRO_32 = (1.0 / 3.0) * At1;
DENDRO_33 = (2.0 / 3.0) * grad_2_beta2;
DENDRO_34 = At1 * DENDRO_igt0 + At3 * DENDRO_igt1 + At4 * DENDRO_igt2;
DENDRO_35 = At4 * DENDRO_igt4;
DENDRO_36 = At3 * DENDRO_igt3 + DENDRO_3 + DENDRO_35;
DENDRO_37 = At1 * DENDRO_igt2 + At3 * DENDRO_igt4 + At4 * DENDRO_igt5;
DENDRO_38 = (1.0 / 3.0) * At2;
DENDRO_39 = (2.0 / 3.0) * grad_1_beta1;
DENDRO_40 = At2 * DENDRO_igt0 + At4 * DENDRO_igt1 + At5 * DENDRO_igt2;
DENDRO_41 = At2 * DENDRO_igt1 + At4 * DENDRO_igt3 + At5 * DENDRO_igt4;
DENDRO_42 = At5 * DENDRO_igt5 + DENDRO_35 + DENDRO_4;
DENDRO_43 = (2.0 / 3.0) * grad_0_beta0;
DENDRO_44 = 2 * At4;
DENDRO_45 = 2 * At3;
DENDRO_46 = (1.0 / 3.0) * At4;

// Dendro: printing variables
//--
At_rhs00[pp] =
    (4.0 / 3.0) * At0 * grad_0_beta0 - DENDRO_0 * grad_1_beta1 -
    DENDRO_0 * grad_2_beta2 + DENDRO_1 * grad_0_beta1 +
    DENDRO_2 * grad_0_beta2 +
    DENDRO_31 * (-DENDRO_30 * gt0 + 3 * DENDRO_6 + 3 * DENDRO_7 + 3 * DENDRO_8 +
                 3 * DENDRO_9 - 3 * grad2_0_0_alpha) -
    alpha * (-At0 * K +
             DENDRO_1 *
                 (At0 * DENDRO_igt1 + At1 * DENDRO_igt3 + At2 * DENDRO_igt4) +
             DENDRO_2 *
                 (At0 * DENDRO_igt2 + At1 * DENDRO_igt4 + At2 * DENDRO_igt5) +
             DENDRO_5 * (At0 * DENDRO_igt0 + DENDRO_3 + DENDRO_4)) +
    beta0 * grad_0_At0 + beta1 * grad_1_At0 + beta2 * grad_2_At0;
//--
At_rhs01[pp] =
    At0 * grad_1_beta0 - At1 * DENDRO_33 + At2 * grad_1_beta2 +
    At3 * grad_0_beta1 + At4 * grad_0_beta2 +
    DENDRO_31 * (3 * DENDRO_18 + 3 * DENDRO_19 + 3 * DENDRO_20 + 3 * DENDRO_21 -
                 DENDRO_30 * gt1 - 3 * grad2_0_1_alpha) +
    DENDRO_32 * grad_0_beta0 + DENDRO_32 * grad_1_beta1 -
    alpha * (-At1 * K + DENDRO_1 * DENDRO_36 + DENDRO_2 * DENDRO_37 +
             DENDRO_34 * DENDRO_5) +
    beta0 * grad_0_At1 + beta1 * grad_1_At1 + beta2 * grad_2_At1;
//--
At_rhs02[pp] =
    At0 * grad_2_beta0 + At1 * grad_2_beta1 - At2 * DENDRO_39 +
    At4 * grad_0_beta1 + At5 * grad_0_beta2 +
    DENDRO_31 * (3 * DENDRO_22 + 3 * DENDRO_23 + 3 * DENDRO_24 + 3 * DENDRO_25 -
                 DENDRO_30 * gt2 - 3 * grad2_0_2_alpha) +
    DENDRO_38 * grad_0_beta0 + DENDRO_38 * grad_2_beta2 -
    alpha * (-At2 * K + DENDRO_1 * DENDRO_41 + DENDRO_2 * DENDRO_42 +
             DENDRO_40 * DENDRO_5) +
    beta0 * grad_0_At2 + beta1 * grad_1_At2 + beta2 * grad_2_At2;
//--
At_rhs11[pp] =
    -At3 * DENDRO_33 - At3 * DENDRO_43 + (4.0 / 3.0) * At3 * grad_1_beta1 +
    DENDRO_1 * grad_1_beta0 +
    DENDRO_31 * (3 * DENDRO_10 + 3 * DENDRO_11 + 3 * DENDRO_12 + 3 * DENDRO_13 -
                 DENDRO_30 * gt3 - 3 * grad2_1_1_alpha) +
    DENDRO_44 * grad_1_beta2 -
    alpha * (-At3 * K + DENDRO_1 * DENDRO_34 + DENDRO_36 * DENDRO_45 +
             DENDRO_37 * DENDRO_44) +
    beta0 * grad_0_At3 + beta1 * grad_1_At3 + beta2 * grad_2_At3;
//--
At_rhs12[pp] =
    At1 * grad_2_beta0 + At2 * grad_1_beta0 + At3 * grad_2_beta1 -
    At4 * DENDRO_43 + At5 * grad_1_beta2 +
    DENDRO_31 * (3 * DENDRO_26 + 3 * DENDRO_27 + 3 * DENDRO_28 + 3 * DENDRO_29 -
                 DENDRO_30 * gt4 - 3 * grad2_1_2_alpha) +
    DENDRO_46 * grad_1_beta1 + DENDRO_46 * grad_2_beta2 -
    alpha * (-At4 * K + DENDRO_1 * DENDRO_40 + DENDRO_41 * DENDRO_45 +
             DENDRO_42 * DENDRO_44) +
    beta0 * grad_0_At4 + beta1 * grad_1_At4 + beta2 * grad_2_At4;
//--
At_rhs22[pp] =
    -At5 * DENDRO_39 - At5 * DENDRO_43 + (4.0 / 3.0) * At5 * grad_2_beta2 +
    DENDRO_2 * grad_2_beta0 +
    DENDRO_31 * (3 * DENDRO_14 + 3 * DENDRO_15 + 3 * DENDRO_16 + 3 * DENDRO_17 -
                 DENDRO_30 * gt5 - 3 * grad2_2_2_alpha) +
    DENDRO_44 * grad_2_beta1 -
    alpha * (2 * At5 * DENDRO_42 - At5 * K + DENDRO_2 * DENDRO_40 +
             DENDRO_41 * DENDRO_44) +
    beta0 * grad_0_At5 + beta1 * grad_1_At5 + beta2 * grad_2_At5;
// Dendro: reduced ops: 362
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 333
// Dendro: printing temp variables
DENDRO_0  = pow(gt4, 2);
DENDRO_1  = pow(gt1, 2);
DENDRO_2  = pow(gt2, 2);
DENDRO_3  = gt0 * gt3;
DENDRO_4  = gt1 * gt4;
DENDRO_5  = chi / (DENDRO_0 * gt0 + DENDRO_1 * gt5 + DENDRO_2 * gt3 -
                  DENDRO_3 * gt5 - 2 * DENDRO_4 * gt2);
DENDRO_6  = 2 * DENDRO_5;
DENDRO_7  = pow(DENDRO_igt1, 2);
DENDRO_8  = pow(DENDRO_igt2, 2);
DENDRO_9  = At1 * DENDRO_igt0;
DENDRO_10 = At2 * DENDRO_igt0;
DENDRO_11 = DENDRO_igt1 * DENDRO_igt2;
DENDRO_12 = 2 * At4;
DENDRO_13 = pow(DENDRO_igt4, 2);
DENDRO_14 = At1 * DENDRO_igt1;
DENDRO_15 = At2 * DENDRO_igt1;
DENDRO_16 = 2 * DENDRO_igt4;
DENDRO_17 = DENDRO_igt3 * DENDRO_igt4;
DENDRO_18 = At1 * DENDRO_igt2;
DENDRO_19 = At2 * DENDRO_igt2;
DENDRO_20 = At4 * DENDRO_igt5;
DENDRO_21 = At0 * DENDRO_igt0;
DENDRO_22 = At3 * DENDRO_igt1;
DENDRO_23 = At4 * DENDRO_igt1;
DENDRO_24 = At4 * DENDRO_igt2;
DENDRO_25 = At5 * DENDRO_igt2;

// Dendro: printing variables
//--
K_rhs[pp] =
    -DENDRO_5 * (-DENDRO_0 + gt3 * gt5) *
        (DENDRO_C3_k0_0 * grad_0_alpha + DENDRO_C3_k1_0 * grad_1_alpha +
         DENDRO_C3_k2_0 * grad_2_alpha - grad2_0_0_alpha) -
    DENDRO_5 * (-DENDRO_1 + DENDRO_3) *
        (DENDRO_C3_k0_5 * grad_0_alpha + DENDRO_C3_k1_5 * grad_1_alpha +
         DENDRO_C3_k2_5 * grad_2_alpha - grad2_2_2_alpha) -
    DENDRO_5 * (-DENDRO_2 + gt0 * gt5) *
        (DENDRO_C3_k0_3 * grad_0_alpha + DENDRO_C3_k1_3 * grad_1_alpha +
         DENDRO_C3_k2_3 * grad_2_alpha - grad2_1_1_alpha) -
    DENDRO_6 * (DENDRO_4 - gt2 * gt3) *
        (DENDRO_C3_k0_2 * grad_0_alpha + DENDRO_C3_k1_2 * grad_1_alpha +
         DENDRO_C3_k2_2 * grad_2_alpha - grad2_0_2_alpha) +
    DENDRO_6 * (gt0 * gt4 - gt1 * gt2) *
        (DENDRO_C3_k0_4 * grad_0_alpha + DENDRO_C3_k1_4 * grad_1_alpha +
         DENDRO_C3_k2_4 * grad_2_alpha - grad2_1_2_alpha) +
    DENDRO_6 * (gt1 * gt5 - gt2 * gt4) *
        (DENDRO_C3_k0_1 * grad_0_alpha + DENDRO_C3_k1_1 * grad_1_alpha +
         DENDRO_C3_k2_1 * grad_2_alpha - grad2_0_1_alpha) +
    (1.0 / 3.0) * alpha *
        (3 * At0 *
             (At0 * pow(DENDRO_igt0, 2) + At3 * DENDRO_7 + At5 * DENDRO_8 +
              2 * DENDRO_10 * DENDRO_igt2 + DENDRO_11 * DENDRO_12 +
              2 * DENDRO_9 * DENDRO_igt1) +
         6 * At1 *
             (At1 * DENDRO_7 + At2 * DENDRO_11 + DENDRO_10 * DENDRO_igt4 +
              DENDRO_21 * DENDRO_igt1 + DENDRO_22 * DENDRO_igt3 +
              DENDRO_23 * DENDRO_igt4 + DENDRO_24 * DENDRO_igt3 +
              DENDRO_25 * DENDRO_igt4 + DENDRO_9 * DENDRO_igt3) +
         6 * At2 *
             (At1 * DENDRO_11 + At2 * DENDRO_8 + DENDRO_10 * DENDRO_igt5 +
              DENDRO_21 * DENDRO_igt2 + DENDRO_22 * DENDRO_igt4 +
              DENDRO_23 * DENDRO_igt5 + DENDRO_24 * DENDRO_igt4 +
              DENDRO_25 * DENDRO_igt5 + DENDRO_9 * DENDRO_igt4) +
         3 * At3 *
             (At0 * DENDRO_7 + At3 * pow(DENDRO_igt3, 2) + At5 * DENDRO_13 +
              DENDRO_12 * DENDRO_17 + 2 * DENDRO_14 * DENDRO_igt3 +
              DENDRO_15 * DENDRO_16) +
         6 * At4 *
             (At0 * DENDRO_11 + At3 * DENDRO_17 + At4 * DENDRO_13 +
              At5 * DENDRO_igt4 * DENDRO_igt5 + DENDRO_14 * DENDRO_igt4 +
              DENDRO_15 * DENDRO_igt5 + DENDRO_18 * DENDRO_igt3 +
              DENDRO_19 * DENDRO_igt4 + DENDRO_20 * DENDRO_igt3) +
         3 * At5 *
             (At0 * DENDRO_8 + At3 * DENDRO_13 + At5 * pow(DENDRO_igt5, 2) +
              DENDRO_16 * DENDRO_18 + DENDRO_16 * DENDRO_20 +
              2 * DENDRO_19 * DENDRO_igt5) +
         pow(K, 2)) +
    beta0 * grad_0_K + beta1 * grad_1_K + beta2 * grad_2_K;
// Dendro: reduced ops: 222
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 1242
// Dendro: printing temp variables
DENDRO_0 = (1.0 / 3.0) * DENDRO_igt0;
DENDRO_1 = (7.0 / 3.0) * DENDRO_igt1;
DENDRO_2 = (1.0 / 3.0) * DENDRO_igt1;
DENDRO_3 = (7.0 / 3.0) * DENDRO_igt2;
DENDRO_4 = (1.0 / 3.0) * DENDRO_igt2;
DENDRO_5 = 2 * DENDRO_igt4;
DENDRO_6 = (2.0 / 3.0) * grad_0_beta0 + (2.0 / 3.0) * grad_1_beta1 +
           (2.0 / 3.0) * grad_2_beta2;
DENDRO_7  = pow(DENDRO_igt1, 2);
DENDRO_8  = pow(DENDRO_igt2, 2);
DENDRO_9  = At1 * DENDRO_igt0;
DENDRO_10 = 2 * DENDRO_igt1;
DENDRO_11 = At2 * DENDRO_igt0;
DENDRO_12 = 2 * DENDRO_igt2;
DENDRO_13 = DENDRO_igt1 * DENDRO_igt2;
DENDRO_14 = At0 * pow(DENDRO_igt0, 2) + At3 * DENDRO_7 + 2 * At4 * DENDRO_13 +
            At5 * DENDRO_8 + DENDRO_10 * DENDRO_9 + DENDRO_11 * DENDRO_12;
DENDRO_15 = 2 * grad_0_alpha;
DENDRO_16 = 2 * alpha;
DENDRO_17 = DENDRO_14 * DENDRO_16;
DENDRO_18 = pow(DENDRO_igt4, 2);
DENDRO_19 = At1 * DENDRO_igt1;
DENDRO_20 = At2 * DENDRO_igt1;
DENDRO_21 = At4 * DENDRO_igt3;
DENDRO_22 = At0 * DENDRO_7 + At3 * pow(DENDRO_igt3, 2) + At5 * DENDRO_18 +
            2 * DENDRO_19 * DENDRO_igt3 + DENDRO_20 * DENDRO_5 +
            DENDRO_21 * DENDRO_5;
DENDRO_23 = DENDRO_16 * DENDRO_22;
DENDRO_24 = At1 * DENDRO_igt2;
DENDRO_25 = At2 * DENDRO_igt2;
DENDRO_26 = At0 * DENDRO_8 + At3 * DENDRO_18 + At4 * DENDRO_5 * DENDRO_igt5 +
            At5 * pow(DENDRO_igt5, 2) + DENDRO_24 * DENDRO_5 +
            2 * DENDRO_25 * DENDRO_igt5;
DENDRO_27 = DENDRO_16 * DENDRO_26;
DENDRO_28 = At0 * DENDRO_igt0;
DENDRO_29 = At3 * DENDRO_igt1;
DENDRO_30 = At4 * DENDRO_igt1;
DENDRO_31 = At4 * DENDRO_igt2;
DENDRO_32 = At5 * DENDRO_igt2;
DENDRO_33 = At1 * DENDRO_7 + At2 * DENDRO_13 + DENDRO_11 * DENDRO_igt4 +
            DENDRO_28 * DENDRO_igt1 + DENDRO_29 * DENDRO_igt3 +
            DENDRO_30 * DENDRO_igt4 + DENDRO_31 * DENDRO_igt3 +
            DENDRO_32 * DENDRO_igt4 + DENDRO_9 * DENDRO_igt3;
DENDRO_34 = 2 * grad_1_alpha;
DENDRO_35 = At1 * DENDRO_13 + At2 * DENDRO_8 + DENDRO_11 * DENDRO_igt5 +
            DENDRO_28 * DENDRO_igt2 + DENDRO_29 * DENDRO_igt4 +
            DENDRO_30 * DENDRO_igt5 + DENDRO_31 * DENDRO_igt4 +
            DENDRO_32 * DENDRO_igt5 + DENDRO_9 * DENDRO_igt4;
DENDRO_36 = 2 * grad_2_alpha;
DENDRO_37 = 4 * alpha;
DENDRO_38 = DENDRO_33 * DENDRO_37;
DENDRO_39 = DENDRO_35 * DENDRO_37;
DENDRO_40 = At0 * DENDRO_13 + At3 * DENDRO_igt3 * DENDRO_igt4 +
            At4 * DENDRO_18 + At5 * DENDRO_igt4 * DENDRO_igt5 +
            DENDRO_19 * DENDRO_igt4 + DENDRO_20 * DENDRO_igt5 +
            DENDRO_21 * DENDRO_igt5 + DENDRO_24 * DENDRO_igt3 +
            DENDRO_25 * DENDRO_igt4;
DENDRO_41 = DENDRO_37 * DENDRO_40;
DENDRO_42 = 4 * grad_0_K;
DENDRO_43 = 9 / chi;
DENDRO_44 = DENDRO_43 * grad_0_chi;
DENDRO_45 = (1.0 / 3.0) * alpha;
DENDRO_46 = 4 * grad_1_K;
DENDRO_47 = DENDRO_43 * grad_1_chi;
DENDRO_48 = 4 * grad_2_K;
DENDRO_49 = DENDRO_43 * grad_2_chi;
DENDRO_50 = (1.0 / 3.0) * DENDRO_igt3;
DENDRO_51 = (1.0 / 3.0) * DENDRO_igt4;
DENDRO_52 = (7.0 / 3.0) * DENDRO_igt4;
DENDRO_53 = (1.0 / 3.0) * DENDRO_igt5;

// Dendro: printing variables
//--
Gt_rhs0[pp] =
    DENDRO_0 * grad2_0_1_beta1 + DENDRO_0 * grad2_0_2_beta2 +
    DENDRO_1 * grad2_0_1_beta0 - DENDRO_14 * DENDRO_15 +
    DENDRO_17 * DENDRO_C2_k0_0 + DENDRO_2 * grad2_1_1_beta1 +
    DENDRO_2 * grad2_1_2_beta2 + DENDRO_23 * DENDRO_C2_k0_3 +
    DENDRO_27 * DENDRO_C2_k0_5 + DENDRO_3 * grad2_0_2_beta0 -
    DENDRO_33 * DENDRO_34 - DENDRO_35 * DENDRO_36 + DENDRO_38 * DENDRO_C2_k0_1 +
    DENDRO_39 * DENDRO_C2_k0_2 + DENDRO_4 * grad2_1_2_beta1 +
    DENDRO_4 * grad2_2_2_beta2 + DENDRO_41 * DENDRO_C2_k0_4 -
    DENDRO_45 * (DENDRO_14 * DENDRO_44 + DENDRO_42 * DENDRO_igt0) -
    DENDRO_45 * (DENDRO_33 * DENDRO_47 + DENDRO_46 * DENDRO_igt1) -
    DENDRO_45 * (DENDRO_35 * DENDRO_49 + DENDRO_48 * DENDRO_igt2) +
    DENDRO_5 * grad2_1_2_beta0 + DENDRO_6 * DENDRO_Gtk0 -
    DENDRO_Gtk0 * grad_0_beta0 - DENDRO_Gtk1 * grad_1_beta0 -
    DENDRO_Gtk2 * grad_2_beta0 + (4.0 / 3.0) * DENDRO_igt0 * grad2_0_0_beta0 +
    DENDRO_igt3 * grad2_1_1_beta0 + DENDRO_igt5 * grad2_2_2_beta0 +
    beta0 * grad_0_Gt0 + beta1 * grad_1_Gt0 + beta2 * grad_2_Gt0;
//--
Gt_rhs1[pp] = DENDRO_1 * grad2_0_1_beta1 + DENDRO_12 * grad2_0_2_beta1 -
              DENDRO_15 * DENDRO_33 + DENDRO_17 * DENDRO_C2_k1_0 +
              DENDRO_2 * grad2_0_0_beta0 + DENDRO_2 * grad2_0_2_beta2 -
              DENDRO_22 * DENDRO_34 + DENDRO_23 * DENDRO_C2_k1_3 +
              DENDRO_27 * DENDRO_C2_k1_5 - DENDRO_36 * DENDRO_40 +
              DENDRO_38 * DENDRO_C2_k1_1 + DENDRO_39 * DENDRO_C2_k1_2 +
              DENDRO_41 * DENDRO_C2_k1_4 -
              DENDRO_45 * (DENDRO_22 * DENDRO_47 + DENDRO_46 * DENDRO_igt3) -
              DENDRO_45 * (DENDRO_33 * DENDRO_44 + DENDRO_42 * DENDRO_igt1) -
              DENDRO_45 * (DENDRO_40 * DENDRO_49 + DENDRO_48 * DENDRO_igt4) +
              DENDRO_50 * grad2_0_1_beta0 + DENDRO_50 * grad2_1_2_beta2 +
              DENDRO_51 * grad2_0_2_beta0 + DENDRO_51 * grad2_2_2_beta2 +
              DENDRO_52 * grad2_1_2_beta1 + DENDRO_6 * DENDRO_Gtk1 -
              DENDRO_Gtk0 * grad_0_beta1 - DENDRO_Gtk1 * grad_1_beta1 -
              DENDRO_Gtk2 * grad_2_beta1 + DENDRO_igt0 * grad2_0_0_beta1 +
              (4.0 / 3.0) * DENDRO_igt3 * grad2_1_1_beta1 +
              DENDRO_igt5 * grad2_2_2_beta1 + beta0 * grad_0_Gt1 +
              beta1 * grad_1_Gt1 + beta2 * grad_2_Gt1;
//--
Gt_rhs2[pp] = DENDRO_10 * grad2_0_1_beta2 - DENDRO_15 * DENDRO_35 +
              DENDRO_17 * DENDRO_C2_k2_0 + DENDRO_23 * DENDRO_C2_k2_3 -
              DENDRO_26 * DENDRO_36 + DENDRO_27 * DENDRO_C2_k2_5 +
              DENDRO_3 * grad2_0_2_beta2 - DENDRO_34 * DENDRO_40 +
              DENDRO_38 * DENDRO_C2_k2_1 + DENDRO_39 * DENDRO_C2_k2_2 +
              DENDRO_4 * grad2_0_0_beta0 + DENDRO_4 * grad2_0_1_beta1 +
              DENDRO_41 * DENDRO_C2_k2_4 -
              DENDRO_45 * (DENDRO_26 * DENDRO_49 + DENDRO_48 * DENDRO_igt5) -
              DENDRO_45 * (DENDRO_35 * DENDRO_44 + DENDRO_42 * DENDRO_igt2) -
              DENDRO_45 * (DENDRO_40 * DENDRO_47 + DENDRO_46 * DENDRO_igt4) +
              DENDRO_51 * grad2_0_1_beta0 + DENDRO_51 * grad2_1_1_beta1 +
              DENDRO_52 * grad2_1_2_beta2 + DENDRO_53 * grad2_0_2_beta0 +
              DENDRO_53 * grad2_1_2_beta1 + DENDRO_6 * DENDRO_Gtk2 -
              DENDRO_Gtk0 * grad_0_beta2 - DENDRO_Gtk1 * grad_1_beta2 -
              DENDRO_Gtk2 * grad_2_beta2 + DENDRO_igt0 * grad2_0_0_beta2 +
              DENDRO_igt3 * grad2_1_1_beta2 +
              (4.0 / 3.0) * DENDRO_igt5 * grad2_2_2_beta2 + beta0 * grad_0_Gt2 +
              beta1 * grad_1_Gt2 + beta2 * grad_2_Gt2;
// Dendro: reduced ops: 367
// Dendro: }}}
// Dendro: {{{
// Dendro: original ops: 1290
// Dendro: printing temp variables
DENDRO_0 = beta0 * grad_0_Gt0 + beta1 * grad_1_Gt0 + beta2 * grad_2_Gt0;
DENDRO_1 = pow(DENDRO_igt1, 2);
DENDRO_2 = pow(DENDRO_igt2, 2);
DENDRO_3 = At1 * DENDRO_igt0;
DENDRO_4 = 2 * DENDRO_igt1;
DENDRO_5 = At2 * DENDRO_igt0;
DENDRO_6 = 2 * DENDRO_igt2;
DENDRO_7 = DENDRO_igt1 * DENDRO_igt2;
DENDRO_8 = At0 * pow(DENDRO_igt0, 2) + At3 * DENDRO_1 + 2 * At4 * DENDRO_7 +
           At5 * DENDRO_2 + DENDRO_3 * DENDRO_4 + DENDRO_5 * DENDRO_6;
DENDRO_9  = 2 * grad_0_alpha;
DENDRO_10 = At0 * DENDRO_igt0;
DENDRO_11 = At3 * DENDRO_igt1;
DENDRO_12 = At4 * DENDRO_igt1;
DENDRO_13 = At4 * DENDRO_igt2;
DENDRO_14 = At5 * DENDRO_igt2;
DENDRO_15 = At1 * DENDRO_1 + At2 * DENDRO_7 + DENDRO_10 * DENDRO_igt1 +
            DENDRO_11 * DENDRO_igt3 + DENDRO_12 * DENDRO_igt4 +
            DENDRO_13 * DENDRO_igt3 + DENDRO_14 * DENDRO_igt4 +
            DENDRO_3 * DENDRO_igt3 + DENDRO_5 * DENDRO_igt4;
DENDRO_16 = 2 * grad_1_alpha;
DENDRO_17 = At1 * DENDRO_7 + At2 * DENDRO_2 + DENDRO_10 * DENDRO_igt2 +
            DENDRO_11 * DENDRO_igt4 + DENDRO_12 * DENDRO_igt5 +
            DENDRO_13 * DENDRO_igt4 + DENDRO_14 * DENDRO_igt5 +
            DENDRO_3 * DENDRO_igt4 + DENDRO_5 * DENDRO_igt5;
DENDRO_18 = 2 * grad_2_alpha;
DENDRO_19 = 2 * DENDRO_igt4;
DENDRO_20 = 4 * grad_0_K;
DENDRO_21 = 9 / chi;
DENDRO_22 = DENDRO_21 * grad_0_chi;
DENDRO_23 = (1.0 / 3.0) * alpha;
DENDRO_24 = 4 * grad_1_K;
DENDRO_25 = DENDRO_21 * grad_1_chi;
DENDRO_26 = 4 * grad_2_K;
DENDRO_27 = DENDRO_21 * grad_2_chi;
DENDRO_28 = (1.0 / 3.0) * DENDRO_igt0;
DENDRO_29 = (1.0 / 3.0) * DENDRO_igt1;
DENDRO_30 = (1.0 / 3.0) * DENDRO_igt2;
DENDRO_31 = (2.0 / 3.0) * grad_0_beta0 + (2.0 / 3.0) * grad_1_beta1 +
            (2.0 / 3.0) * grad_2_beta2;
DENDRO_32 = (7.0 / 3.0) * DENDRO_igt1;
DENDRO_33 = (7.0 / 3.0) * DENDRO_igt2;
DENDRO_34 = 2 * alpha;
DENDRO_35 = DENDRO_34 * DENDRO_8;
DENDRO_36 = pow(DENDRO_igt4, 2);
DENDRO_37 = At1 * DENDRO_igt1;
DENDRO_38 = At2 * DENDRO_igt1;
DENDRO_39 = At4 * DENDRO_igt3;
DENDRO_40 = At0 * DENDRO_1 + At3 * pow(DENDRO_igt3, 2) + At5 * DENDRO_36 +
            DENDRO_19 * DENDRO_38 + DENDRO_19 * DENDRO_39 +
            2 * DENDRO_37 * DENDRO_igt3;
DENDRO_41 = DENDRO_34 * DENDRO_40;
DENDRO_42 = At1 * DENDRO_igt2;
DENDRO_43 = At2 * DENDRO_igt2;
DENDRO_44 = At0 * DENDRO_2 + At3 * DENDRO_36 + At4 * DENDRO_19 * DENDRO_igt5 +
            At5 * pow(DENDRO_igt5, 2) + DENDRO_19 * DENDRO_42 +
            2 * DENDRO_43 * DENDRO_igt5;
DENDRO_45 = DENDRO_34 * DENDRO_44;
DENDRO_46 = 4 * alpha;
DENDRO_47 = DENDRO_15 * DENDRO_46;
DENDRO_48 = DENDRO_17 * DENDRO_46;
DENDRO_49 = At0 * DENDRO_7 + At3 * DENDRO_igt3 * DENDRO_igt4 + At4 * DENDRO_36 +
            At5 * DENDRO_igt4 * DENDRO_igt5 + DENDRO_37 * DENDRO_igt4 +
            DENDRO_38 * DENDRO_igt5 + DENDRO_39 * DENDRO_igt5 +
            DENDRO_42 * DENDRO_igt3 + DENDRO_43 * DENDRO_igt4;
DENDRO_50 = DENDRO_46 * DENDRO_49;
DENDRO_51 = beta0 * grad_0_Gt1 + beta1 * grad_1_Gt1 + beta2 * grad_2_Gt1;
DENDRO_52 = (1.0 / 3.0) * DENDRO_igt3;
DENDRO_53 = (1.0 / 3.0) * DENDRO_igt4;
DENDRO_54 = (7.0 / 3.0) * DENDRO_igt4;
DENDRO_55 = beta0 * grad_0_Gt2 + beta1 * grad_1_Gt2 + beta2 * grad_2_Gt2;
DENDRO_56 = (1.0 / 3.0) * DENDRO_igt5;

// Dendro: printing variables
//--
B_rhs0[pp] =
    -B0 * eta - DENDRO_0 * lambda[3] + DENDRO_0 - DENDRO_15 * DENDRO_16 -
    DENDRO_17 * DENDRO_18 + DENDRO_19 * grad2_1_2_beta0 -
    DENDRO_23 * (DENDRO_15 * DENDRO_25 + DENDRO_24 * DENDRO_igt1) -
    DENDRO_23 * (DENDRO_17 * DENDRO_27 + DENDRO_26 * DENDRO_igt2) -
    DENDRO_23 * (DENDRO_20 * DENDRO_igt0 + DENDRO_22 * DENDRO_8) +
    DENDRO_28 * grad2_0_1_beta1 + DENDRO_28 * grad2_0_2_beta2 +
    DENDRO_29 * grad2_1_1_beta1 + DENDRO_29 * grad2_1_2_beta2 +
    DENDRO_30 * grad2_1_2_beta1 + DENDRO_30 * grad2_2_2_beta2 +
    DENDRO_31 * DENDRO_Gtk0 + DENDRO_32 * grad2_0_1_beta0 +
    DENDRO_33 * grad2_0_2_beta0 + DENDRO_35 * DENDRO_C2_k0_0 +
    DENDRO_41 * DENDRO_C2_k0_3 + DENDRO_45 * DENDRO_C2_k0_5 +
    DENDRO_47 * DENDRO_C2_k0_1 + DENDRO_48 * DENDRO_C2_k0_2 +
    DENDRO_50 * DENDRO_C2_k0_4 - DENDRO_8 * DENDRO_9 -
    DENDRO_Gtk0 * grad_0_beta0 - DENDRO_Gtk1 * grad_1_beta0 -
    DENDRO_Gtk2 * grad_2_beta0 + (4.0 / 3.0) * DENDRO_igt0 * grad2_0_0_beta0 +
    DENDRO_igt3 * grad2_1_1_beta0 + DENDRO_igt5 * grad2_2_2_beta0 +
    lambda[2] * (beta0 * grad_0_B0 + beta1 * grad_1_B0 + beta2 * grad_2_B0);
//--
B_rhs1[pp] =
    -B1 * eta - DENDRO_15 * DENDRO_9 - DENDRO_16 * DENDRO_40 -
    DENDRO_18 * DENDRO_49 -
    DENDRO_23 * (DENDRO_15 * DENDRO_22 + DENDRO_20 * DENDRO_igt1) -
    DENDRO_23 * (DENDRO_24 * DENDRO_igt3 + DENDRO_25 * DENDRO_40) -
    DENDRO_23 * (DENDRO_26 * DENDRO_igt4 + DENDRO_27 * DENDRO_49) +
    DENDRO_29 * grad2_0_0_beta0 + DENDRO_29 * grad2_0_2_beta2 +
    DENDRO_31 * DENDRO_Gtk1 + DENDRO_32 * grad2_0_1_beta1 +
    DENDRO_35 * DENDRO_C2_k1_0 + DENDRO_41 * DENDRO_C2_k1_3 +
    DENDRO_45 * DENDRO_C2_k1_5 + DENDRO_47 * DENDRO_C2_k1_1 +
    DENDRO_48 * DENDRO_C2_k1_2 + DENDRO_50 * DENDRO_C2_k1_4 -
    DENDRO_51 * lambda[3] + DENDRO_51 + DENDRO_52 * grad2_0_1_beta0 +
    DENDRO_52 * grad2_1_2_beta2 + DENDRO_53 * grad2_0_2_beta0 +
    DENDRO_53 * grad2_2_2_beta2 + DENDRO_54 * grad2_1_2_beta1 +
    DENDRO_6 * grad2_0_2_beta1 - DENDRO_Gtk0 * grad_0_beta1 -
    DENDRO_Gtk1 * grad_1_beta1 - DENDRO_Gtk2 * grad_2_beta1 +
    DENDRO_igt0 * grad2_0_0_beta1 +
    (4.0 / 3.0) * DENDRO_igt3 * grad2_1_1_beta1 +
    DENDRO_igt5 * grad2_2_2_beta1 +
    lambda[2] * (beta0 * grad_0_B1 + beta1 * grad_1_B1 + beta2 * grad_2_B1);
//--
B_rhs2[pp] =
    -B2 * eta - DENDRO_16 * DENDRO_49 - DENDRO_17 * DENDRO_9 -
    DENDRO_18 * DENDRO_44 -
    DENDRO_23 * (DENDRO_17 * DENDRO_22 + DENDRO_20 * DENDRO_igt2) -
    DENDRO_23 * (DENDRO_24 * DENDRO_igt4 + DENDRO_25 * DENDRO_49) -
    DENDRO_23 * (DENDRO_26 * DENDRO_igt5 + DENDRO_27 * DENDRO_44) +
    DENDRO_30 * grad2_0_0_beta0 + DENDRO_30 * grad2_0_1_beta1 +
    DENDRO_31 * DENDRO_Gtk2 + DENDRO_33 * grad2_0_2_beta2 +
    DENDRO_35 * DENDRO_C2_k2_0 + DENDRO_4 * grad2_0_1_beta2 +
    DENDRO_41 * DENDRO_C2_k2_3 + DENDRO_45 * DENDRO_C2_k2_5 +
    DENDRO_47 * DENDRO_C2_k2_1 + DENDRO_48 * DENDRO_C2_k2_2 +
    DENDRO_50 * DENDRO_C2_k2_4 + DENDRO_53 * grad2_0_1_beta0 +
    DENDRO_53 * grad2_1_1_beta1 + DENDRO_54 * grad2_1_2_beta2 -
    DENDRO_55 * lambda[3] + DENDRO_55 + DENDRO_56 * grad2_0_2_beta0 +
    DENDRO_56 * grad2_1_2_beta1 - DENDRO_Gtk0 * grad_0_beta2 -
    DENDRO_Gtk1 * grad_1_beta2 - DENDRO_Gtk2 * grad_2_beta2 +
    DENDRO_igt0 * grad2_0_0_beta2 + DENDRO_igt3 * grad2_1_1_beta2 +
    (4.0 / 3.0) * DENDRO_igt5 * grad2_2_2_beta2 +
    lambda[2] * (beta0 * grad_0_B2 + beta1 * grad_1_B2 + beta2 * grad_2_B2);
// Dendro: reduced ops: 400
// Dendro: }}}
