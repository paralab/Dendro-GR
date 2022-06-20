// Dendro: {{{ 
// Dendro: original ops:  49
// Dendro: printing temp variables
double DENDRO_0 = 4.0*PI;

// Dendro: printing variables
//--
A_rhs0[pp] = -E0[pp] - grad_0_psi[pp];
//--
A_rhs1[pp] = -E1[pp] - grad_1_psi[pp];
//--
A_rhs2[pp] = -E2[pp] - grad_2_psi[pp];
//--
E_rhs0[pp] = -DENDRO_0*J0[pp] + grad_0_Gamma[pp] - grad2_0_0_A0[pp] - grad2_1_1_A0[pp] - grad2_2_2_A0[pp];
//--
E_rhs1[pp] = -DENDRO_0*J1[pp] + grad_1_Gamma[pp] - grad2_0_0_A1[pp] - grad2_1_1_A1[pp] - grad2_2_2_A1[pp];
//--
E_rhs2[pp] = -DENDRO_0*J2[pp] + grad_2_Gamma[pp] - grad2_0_0_A2[pp] - grad2_1_1_A2[pp] - grad2_2_2_A2[pp];
//--
Gamma_rhs[pp] = -DENDRO_0*rho_e[pp] - grad2_0_0_psi[pp] - grad2_1_1_psi[pp] - grad2_2_2_psi[pp];
//--
psi_rhs[pp] = -Gamma[pp];
// Dendro: reduced ops:  46
// Dendro: }}} 
