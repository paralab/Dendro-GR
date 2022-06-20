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
// Dendro: original ops:  2121
// Dendro: printing temp variables
double DENDRO_0 = pow(gt4[pp], 2);
double DENDRO_1 = pow(gt1[pp], 2);
double DENDRO_2 = pow(gt2[pp], 2);
double DENDRO_3 = gt0[pp]*gt3[pp];
double DENDRO_4 = gt1[pp]*gt2[pp];
double DENDRO_5 = 1.0/(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt5[pp] - 2*DENDRO_4*gt4[pp]);
double DENDRO_6 = DENDRO_5*alpha[pp];
double DENDRO_7 = grad_2_K[pp];
double DENDRO_8 = (4.0L/3.0L)*gt1[pp]*gt4[pp] - 4.0L/3.0L*gt2[pp]*gt3[pp];
double DENDRO_9 = grad_1_K[pp];
double DENDRO_10 = (4.0L/3.0L)*gt1[pp]*gt5[pp] - 4.0L/3.0L*gt2[pp]*gt4[pp];
double DENDRO_11 = grad_0_K[pp];
double DENDRO_12 = -DENDRO_0 + gt3[pp]*gt5[pp];
double DENDRO_13 = 1.0/chi[pp];
double DENDRO_14 = 3*DENDRO_13*DENDRO_5*grad_0_chi[pp];
double DENDRO_15 = gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp];
double DENDRO_16 = pow(DENDRO_15, 2);
double DENDRO_17 = gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp];
double DENDRO_18 = pow(DENDRO_17, 2);
double DENDRO_19 = DENDRO_15*DENDRO_17;
double DENDRO_20 = 2*At1[pp]*DENDRO_15;
double DENDRO_21 = 2*At2[pp]*DENDRO_17;
double DENDRO_22 = 3*DENDRO_13*DENDRO_5*grad_2_chi[pp];
double DENDRO_23 = -DENDRO_4 + gt0[pp]*gt4[pp];
double DENDRO_24 = DENDRO_15*DENDRO_23;
double DENDRO_25 = At0[pp]*DENDRO_12;
double DENDRO_26 = DENDRO_17*DENDRO_23;
double DENDRO_27 = -DENDRO_1 + DENDRO_3;
double DENDRO_28 = At5[pp]*DENDRO_27;
double DENDRO_29 = DENDRO_12*DENDRO_23;
double DENDRO_30 = DENDRO_15*DENDRO_27;
double DENDRO_31 = -At1[pp]*DENDRO_19 - At1[pp]*DENDRO_29 + At2[pp]*DENDRO_12*DENDRO_27 + At2[pp]*DENDRO_18 + At3[pp]*DENDRO_24 - At4[pp]*DENDRO_26 - At4[pp]*DENDRO_30 + DENDRO_17*DENDRO_25 + DENDRO_17*DENDRO_28;
double DENDRO_32 = 3*DENDRO_13*DENDRO_5*grad_1_chi[pp];
double DENDRO_33 = -DENDRO_2 + gt0[pp]*gt5[pp];
double DENDRO_34 = At3[pp]*DENDRO_33;
double DENDRO_35 = DENDRO_17*DENDRO_33;
double DENDRO_36 = -At1[pp]*DENDRO_12*DENDRO_33 - At1[pp]*DENDRO_16 + At2[pp]*DENDRO_19 + At2[pp]*DENDRO_29 - At4[pp]*DENDRO_24 - At4[pp]*DENDRO_35 + At5[pp]*DENDRO_26 + DENDRO_15*DENDRO_25 + DENDRO_15*DENDRO_34;
double DENDRO_37 = (4.0L/3.0L)*gt0[pp]*gt4[pp] - 4.0L/3.0L*gt1[pp]*gt2[pp];
double DENDRO_38 = pow(DENDRO_23, 2);
double DENDRO_39 = 2*At4[pp]*DENDRO_23;
double DENDRO_40 = At0[pp]*DENDRO_19 - At1[pp]*DENDRO_24 - At1[pp]*DENDRO_35 + At2[pp]*DENDRO_26 + At2[pp]*DENDRO_30 - At4[pp]*DENDRO_27*DENDRO_33 - At4[pp]*DENDRO_38 + DENDRO_23*DENDRO_28 + DENDRO_23*DENDRO_34;
// Dendro: printing variables
//--
Gt_rhs_s7_0[pp] = DENDRO_6*(DENDRO_10*DENDRO_9 - 4.0L/3.0L*DENDRO_11*DENDRO_12 + DENDRO_14*(At0[pp]*pow(DENDRO_12, 2) + At3[pp]*DENDRO_16 - 2*At4[pp]*DENDRO_19 + At5[pp]*DENDRO_18 - DENDRO_12*DENDRO_20 + DENDRO_12*DENDRO_21) + DENDRO_22*DENDRO_31 - DENDRO_32*DENDRO_36 - DENDRO_7*DENDRO_8);
//--
Gt_rhs_s7_1[pp] = DENDRO_6*(DENDRO_10*DENDRO_11 - DENDRO_14*DENDRO_36 - DENDRO_22*DENDRO_40 + DENDRO_32*(At0[pp]*DENDRO_16 + 2*At2[pp]*DENDRO_24 + At3[pp]*pow(DENDRO_33, 2) + At5[pp]*DENDRO_38 - DENDRO_20*DENDRO_33 - DENDRO_33*DENDRO_39) - 4.0L/3.0L*DENDRO_33*DENDRO_9 + DENDRO_37*DENDRO_7);
//--
Gt_rhs_s7_2[pp] = DENDRO_6*(-DENDRO_11*DENDRO_8 + DENDRO_14*DENDRO_31 + DENDRO_22*(At0[pp]*DENDRO_18 - 2*At1[pp]*DENDRO_26 + At3[pp]*DENDRO_38 + At5[pp]*pow(DENDRO_27, 2) + DENDRO_21*DENDRO_27 - DENDRO_27*DENDRO_39) - 4.0L/3.0L*DENDRO_27*DENDRO_7 - DENDRO_32*DENDRO_40 + DENDRO_37*DENDRO_9);
// Dendro: reduced ops:  220
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
