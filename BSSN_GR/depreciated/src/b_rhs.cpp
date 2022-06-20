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
// Dendro: original ops:  66
// Dendro: printing temp variables
double DENDRO_0 = (3.0L/4.0L)*alpha[pp]*lambda_f[1] + (3.0L/4.0L)*lambda_f[0];
// Dendro: printing variables
//--
b_rhs0[pp] = B0[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta0[pp] + beta1[pp]*agrad_1_beta0[pp] + beta2[pp]*agrad_2_beta0[pp]);
//--
b_rhs1[pp] = B1[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta1[pp] + beta1[pp]*agrad_1_beta1[pp] + beta2[pp]*agrad_2_beta1[pp]);
//--
b_rhs2[pp] = B2[pp]*DENDRO_0 + lambda[1]*(beta0[pp]*agrad_0_beta2[pp] + beta1[pp]*agrad_1_beta2[pp] + beta2[pp]*agrad_2_beta2[pp]);
//--
//b_rhs3[pp] = kograd_0_beta0[pp] + kograd_1_beta0[pp] + kograd_2_beta0[pp];
//--
//b_rhs4[pp] = kograd_0_beta1[pp] + kograd_1_beta1[pp] + kograd_2_beta1[pp];
//--
//b_rhs5[pp] = kograd_0_beta2[pp] + kograd_1_beta2[pp] + kograd_2_beta2[pp];
// Dendro: reduced ops:  54
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
