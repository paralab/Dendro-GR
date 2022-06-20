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
// Dendro: original ops:  909
// Dendro: printing temp variables
double DENDRO_0 = pow(gt4[pp], 2);
double DENDRO_1 = pow(gt1[pp], 2);
double DENDRO_2 = pow(gt2[pp], 2);
double DENDRO_3 = gt0[pp]*gt3[pp];
double DENDRO_4 = gt1[pp]*gt2[pp];
double DENDRO_5 = 1.0/(DENDRO_0*gt0[pp] + DENDRO_1*gt5[pp] + DENDRO_2*gt3[pp] - DENDRO_3*gt5[pp] - 2*DENDRO_4*gt4[pp]);
double DENDRO_6 = grad2_0_2_beta0[pp];
double DENDRO_7 = (7.0L/3.0L)*gt1[pp]*gt4[pp] - 7.0L/3.0L*gt2[pp]*gt3[pp];
double DENDRO_8 = grad2_1_2_beta1[pp];
double DENDRO_9 = (1.0L/3.0L)*gt1[pp]*gt4[pp] - 1.0L/3.0L*gt2[pp]*gt3[pp];
double DENDRO_10 = grad2_2_2_beta2[pp];
double DENDRO_11 = grad2_0_1_beta0[pp];
double DENDRO_12 = (7.0L/3.0L)*gt1[pp]*gt5[pp] - 7.0L/3.0L*gt2[pp]*gt4[pp];
double DENDRO_13 = grad2_1_1_beta1[pp];
double DENDRO_14 = (1.0L/3.0L)*gt1[pp]*gt5[pp] - 1.0L/3.0L*gt2[pp]*gt4[pp];
double DENDRO_15 = grad2_1_2_beta2[pp];
double DENDRO_16 = -DENDRO_1 + DENDRO_3;
double DENDRO_17 = -DENDRO_2 + gt0[pp]*gt5[pp];
double DENDRO_18 = grad2_0_0_beta0[pp];
double DENDRO_19 = -DENDRO_0 + gt3[pp]*gt5[pp];
double DENDRO_20 = grad2_0_1_beta1[pp];
double DENDRO_21 = -1.0L/3.0L*DENDRO_0 + (1.0L/3.0L)*gt3[pp]*gt5[pp];
double DENDRO_22 = grad2_0_2_beta2[pp];
double DENDRO_23 = (1.0L/3.0L)*gt0[pp]*gt4[pp] - 1.0L/3.0L*gt1[pp]*gt2[pp];
double DENDRO_24 = (7.0L/3.0L)*gt0[pp]*gt4[pp] - 7.0L/3.0L*gt1[pp]*gt2[pp];
double DENDRO_25 = -1.0L/3.0L*DENDRO_2 + (1.0L/3.0L)*gt0[pp]*gt5[pp];
double DENDRO_26 = -1.0L/3.0L*DENDRO_1 + (1.0L/3.0L)*gt0[pp]*gt3[pp];
// Dendro: printing variables
//--
Gt_rhs_s4_0[pp] = DENDRO_5*(-DENDRO_10*DENDRO_9 + DENDRO_11*DENDRO_12 + DENDRO_13*DENDRO_14 + DENDRO_14*DENDRO_15 - DENDRO_16*grad2_2_2_beta0[pp] - DENDRO_17*grad2_1_1_beta0[pp] - 4.0L/3.0L*DENDRO_18*DENDRO_19 - DENDRO_20*DENDRO_21 - DENDRO_21*DENDRO_22 - DENDRO_6*DENDRO_7 - DENDRO_8*DENDRO_9 + 2*(-DENDRO_4 + gt0[pp]*gt4[pp])*grad2_1_2_beta0[pp]);
//--
Gt_rhs_s4_1[pp] = DENDRO_5*(DENDRO_10*DENDRO_23 - DENDRO_11*DENDRO_25 + DENDRO_12*DENDRO_20 - 4.0L/3.0L*DENDRO_13*DENDRO_17 + DENDRO_14*DENDRO_18 + DENDRO_14*DENDRO_22 - DENDRO_15*DENDRO_25 - DENDRO_16*grad2_2_2_beta1[pp] - DENDRO_19*grad2_0_0_beta1[pp] + DENDRO_23*DENDRO_6 + DENDRO_24*DENDRO_8 - 2*(gt1[pp]*gt4[pp] - gt2[pp]*gt3[pp])*grad2_0_2_beta1[pp]);
//--
Gt_rhs_s4_2[pp] = DENDRO_5*(-4.0L/3.0L*DENDRO_10*DENDRO_16 + DENDRO_11*DENDRO_23 + DENDRO_13*DENDRO_23 + DENDRO_15*DENDRO_24 - DENDRO_17*grad2_1_1_beta2[pp] - DENDRO_18*DENDRO_9 - DENDRO_19*grad2_0_0_beta2[pp] - DENDRO_20*DENDRO_9 - DENDRO_22*DENDRO_7 - DENDRO_26*DENDRO_6 - DENDRO_26*DENDRO_8 + 2*(gt1[pp]*gt5[pp] - gt2[pp]*gt4[pp])*grad2_0_1_beta2[pp]);
// Dendro: reduced ops:  176
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
