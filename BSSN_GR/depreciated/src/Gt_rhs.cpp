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
// Dendro: original ops:  24
// Dendro: printing temp variables
// Dendro: printing variables
//--
Gt_rhs0[pp] = Gt_rhs_s1_0[pp] - Gt_rhs_s2_0[pp] + (2.0L/3.0L)*Gt_rhs_s3_0[pp] + Gt_rhs_s4_0[pp] - Gt_rhs_s5_0[pp] + Gt_rhs_s6_0[pp] - Gt_rhs_s7_0[pp];
//--
Gt_rhs1[pp] = Gt_rhs_s1_1[pp] - Gt_rhs_s2_1[pp] + (2.0L/3.0L)*Gt_rhs_s3_1[pp] + Gt_rhs_s4_1[pp] - Gt_rhs_s5_1[pp] + Gt_rhs_s6_1[pp] - Gt_rhs_s7_1[pp];
//--
Gt_rhs2[pp] = Gt_rhs_s1_2[pp] - Gt_rhs_s2_2[pp] + (2.0L/3.0L)*Gt_rhs_s3_2[pp] + Gt_rhs_s4_2[pp] - Gt_rhs_s5_2[pp] + Gt_rhs_s6_2[pp] - Gt_rhs_s7_2[pp];
// Dendro: reduced ops:  24
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
