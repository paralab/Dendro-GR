// Dendro: {{{
// Dendro: original ops: 6898
// Dendro: printing temp variables
const double DENDRO_0  = grad_1_gt5[pp];
const double DENDRO_1  = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
const double DENDRO_2  = grad_2_gt3[pp];
const double DENDRO_3  = gt2[pp] * gt4[pp];
const double DENDRO_4  = -DENDRO_3 + gt1[pp] * gt5[pp];
const double DENDRO_5  = pow(gt4[pp], 2);
const double DENDRO_6  = -DENDRO_5 + gt3[pp] * gt5[pp];
const double DENDRO_7  = grad_1_gt2[pp];
const double DENDRO_8  = grad_2_gt1[pp];
const double DENDRO_9  = grad_0_gt4[pp];
const double DENDRO_10 = DENDRO_7 + DENDRO_8 - DENDRO_9;
const double DENDRO_11 = gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp];
const double DENDRO_12 = pow(gt1[pp], 2);
const double DENDRO_13 = pow(gt2[pp], 2);
const double DENDRO_14 = gt0[pp] * gt3[pp];
const double DENDRO_15 = 2 * gt1[pp];
const double DENDRO_16 =
    pow(DENDRO_12 * gt5[pp] + DENDRO_13 * gt3[pp] - DENDRO_14 * gt5[pp] -
            DENDRO_15 * DENDRO_3 + DENDRO_5 * gt0[pp],
        -2);
const double DENDRO_17 = 1.0 * DENDRO_16;
const double DENDRO_18 = DENDRO_11 * DENDRO_17;
const double DENDRO_19 = grad_0_gt5[pp];
const double DENDRO_20 = grad_2_gt0[pp];
const double DENDRO_21 = -DENDRO_7 + DENDRO_8 + DENDRO_9;
const double DENDRO_22 = DENDRO_1 * DENDRO_17;
const double DENDRO_23 = grad_0_gt3[pp];
const double DENDRO_24 = grad_1_gt0[pp];
const double DENDRO_25 = DENDRO_7 - DENDRO_8 + DENDRO_9;
const double DENDRO_26 = DENDRO_17 * DENDRO_4;
const double DENDRO_27 = 0.5 * grad_2_gt5[pp];
const double DENDRO_28 = 0.5 * DENDRO_0 - 1.0 * grad_2_gt4[pp];
const double DENDRO_29 = 0.5 * DENDRO_19 - 1.0 * grad_2_gt2[pp];
const double DENDRO_30 = -DENDRO_12 + DENDRO_14;
const double DENDRO_31 = DENDRO_16 * DENDRO_30;
const double DENDRO_32 = 0.5 * grad_1_gt3[pp];
const double DENDRO_33 = -0.5 * DENDRO_2 + 1.0 * grad_1_gt4[pp];
const double DENDRO_34 = 0.5 * DENDRO_23 - 1.0 * grad_1_gt1[pp];
const double DENDRO_35 = -DENDRO_13 + gt0[pp] * gt5[pp];
const double DENDRO_36 = DENDRO_16 * DENDRO_35;
const double DENDRO_37 = 0.5 * grad_0_gt0[pp];
const double DENDRO_38 = -0.5 * DENDRO_20 + 1.0 * grad_0_gt2[pp];
const double DENDRO_39 = -0.5 * DENDRO_24 + 1.0 * grad_0_gt1[pp];
const double DENDRO_40 = DENDRO_16 * DENDRO_6;
const double DENDRO_41 =
    -DENDRO_18 *
        (-DENDRO_0 * DENDRO_1 - DENDRO_10 * DENDRO_6 + DENDRO_2 * DENDRO_4) +
    DENDRO_22 *
        (-DENDRO_1 * DENDRO_19 - DENDRO_20 * DENDRO_6 + DENDRO_21 * DENDRO_4) -
    DENDRO_26 *
        (-DENDRO_1 * DENDRO_25 + DENDRO_23 * DENDRO_4 - DENDRO_24 * DENDRO_6) +
    DENDRO_31 *
        (-DENDRO_1 * DENDRO_27 - DENDRO_28 * DENDRO_4 + DENDRO_29 * DENDRO_6) +
    DENDRO_36 *
        (-DENDRO_1 * DENDRO_33 + DENDRO_32 * DENDRO_4 + DENDRO_34 * DENDRO_6) +
    DENDRO_40 *
        (-DENDRO_1 * DENDRO_38 - DENDRO_37 * DENDRO_6 + DENDRO_39 * DENDRO_4) +
    Gt0[pp];
const double DENDRO_42 =
    -DENDRO_18 *
        (DENDRO_0 * DENDRO_11 + DENDRO_10 * DENDRO_4 - DENDRO_2 * DENDRO_35) +
    DENDRO_22 *
        (DENDRO_11 * DENDRO_19 + DENDRO_20 * DENDRO_4 - DENDRO_21 * DENDRO_35) -
    DENDRO_26 *
        (DENDRO_11 * DENDRO_25 - DENDRO_23 * DENDRO_35 + DENDRO_24 * DENDRO_4) +
    DENDRO_31 *
        (DENDRO_11 * DENDRO_27 + DENDRO_28 * DENDRO_35 - DENDRO_29 * DENDRO_4) +
    DENDRO_36 *
        (DENDRO_11 * DENDRO_33 - DENDRO_32 * DENDRO_35 - DENDRO_34 * DENDRO_4) +
    DENDRO_40 *
        (DENDRO_11 * DENDRO_38 - DENDRO_35 * DENDRO_39 + DENDRO_37 * DENDRO_4) +
    Gt1[pp];
const double DENDRO_43 =
    -DENDRO_18 *
        (-DENDRO_0 * DENDRO_30 - DENDRO_1 * DENDRO_10 + DENDRO_11 * DENDRO_2) +
    DENDRO_22 * (-DENDRO_1 * DENDRO_20 + DENDRO_11 * DENDRO_21 -
                 DENDRO_19 * DENDRO_30) -
    DENDRO_26 * (-DENDRO_1 * DENDRO_24 + DENDRO_11 * DENDRO_23 -
                 DENDRO_25 * DENDRO_30) +
    DENDRO_31 *
        (DENDRO_1 * DENDRO_29 - DENDRO_11 * DENDRO_28 - DENDRO_27 * DENDRO_30) +
    DENDRO_36 *
        (DENDRO_1 * DENDRO_34 + DENDRO_11 * DENDRO_32 - DENDRO_30 * DENDRO_33) +
    DENDRO_40 * (-DENDRO_1 * DENDRO_37 + DENDRO_11 * DENDRO_39 -
                 DENDRO_30 * DENDRO_38) +
    Gt2[pp];
const double DENDRO_44 = -DENDRO_42;
const double DENDRO_45 = -DENDRO_41;
const double DENDRO_46 = -DENDRO_43;
const double DENDRO_47 = 2 * DENDRO_46;

// Dendro: printing variables
//--
cgt0[pp]               = DENDRO_41;
//--
cgt1[pp]               = DENDRO_42;
//--
cgt2[pp]               = DENDRO_43;
//--
cgt_mag[pp] =
    sqrt(DENDRO_15 * DENDRO_44 * DENDRO_45 + (DENDRO_44 * DENDRO_44) * gt3[pp] +
         DENDRO_44 * DENDRO_47 * gt4[pp] + (DENDRO_45 * DENDRO_45) * gt0[pp] +
         DENDRO_45 * DENDRO_47 * gt2[pp] + (DENDRO_46 * DENDRO_46) * gt5[pp]);
// Dendro: reduced ops: 231
// Dendro: }}}
