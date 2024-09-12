bssn::timer::t_rhs.start();
for (unsigned int k = PW; k < nz - PW; k++) {
    z = pmin[2] + k * hz;
    for (unsigned int j = PW; j < ny - PW; j++) {
        y = pmin[1] + j * hy;
        for (unsigned int i = PW; i < nx - PW; i++) {
            x       = pmin[0] + i * hx;
            pp      = i + nx * (j + ny * k);
            r_coord = sqrt(x * x + y * y + z * z);
            eta     = ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
            }
            // Dendro: {{{
            // Dendro: original ops:  24
            // Dendro: printing temp variables
            // Dendro: printing variables
            //--
            Gt_rhs_s1_0[pp] = beta0[pp] * agrad_0_Gt0[pp] +
                              beta1[pp] * agrad_1_Gt0[pp] +
                              beta2[pp] * agrad_2_Gt0[pp];
            //--
            Gt_rhs_s1_1[pp] = beta0[pp] * agrad_0_Gt1[pp] +
                              beta1[pp] * agrad_1_Gt1[pp] +
                              beta2[pp] * agrad_2_Gt1[pp];
            //--
            Gt_rhs_s1_2[pp] = beta0[pp] * agrad_0_Gt2[pp] +
                              beta1[pp] * agrad_1_Gt2[pp] +
                              beta2[pp] * agrad_2_Gt2[pp];
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
