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
// Dendro: original ops:  12
// Dendro: printing temp variables
// Dendro: printing variables
//--
#ifdef BSSN_KERR_SCHILD_TEST
            a_rhs[pp] = 0.0;
#else
            a_rhs[pp] = -2 * K[pp] * alpha[pp] +
                        lambda[0] * (beta0[pp] * agrad_0_alpha[pp] +
                                     beta1[pp] * agrad_1_alpha[pp] +
                                     beta2[pp] * agrad_2_alpha[pp]);
#endif
            // Dendro: reduced ops:  12
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
