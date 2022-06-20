// Dendro: {{{ 
// Dendro: original ops:  56
// Dendro: printing temp variables
double DENDRO_0 = 4.0*PI;

// Dendro: printing variables
//--
B_rhs0[pp] = grad_0_Phi[pp] - grad_1_E2[pp] + grad_2_E1[pp];
//--
B_rhs1[pp] = grad_0_E2[pp] + grad_1_Phi[pp] - grad_2_E0[pp];
//--
B_rhs2[pp] = -grad_0_E1[pp] + grad_1_E0[pp] + grad_2_Phi[pp];
//--
E_rhs0[pp] = -DENDRO_0*J0[pp] - grad_0_Psi[pp] + grad_1_B2[pp] - grad_2_B1[pp];
//--
E_rhs1[pp] = -DENDRO_0*J1[pp] - grad_0_B2[pp] - grad_1_Psi[pp] + grad_2_B0[pp];
//--
E_rhs2[pp] = -DENDRO_0*J2[pp] + grad_0_B1[pp] - grad_1_B0[pp] - grad_2_Psi[pp];
//--
Psi_rhs[pp] = 4*PI*rho_e[pp] - Psi[pp]*kappa_1 - grad_0_E0[pp] - grad_1_E1[pp] - grad_2_E2[pp];
//--
Phi_rhs[pp] = -Phi[pp]*kappa_2 + grad_0_B0[pp] + grad_1_B1[pp] + grad_2_B2[pp];
// Dendro: reduced ops:  54
// Dendro: }}} 
