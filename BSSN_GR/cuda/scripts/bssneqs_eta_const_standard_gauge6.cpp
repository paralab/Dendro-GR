device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, alpha, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_alpha = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_alpha = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_alpha = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_alpha = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_alpha = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_alpha = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_alpha = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_alpha = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_alpha = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_alpha = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_alpha = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_alpha = Du[gidx];
__syncthreads();

{
    // Dendro: {{{
    // Dendro: original ops: 9
    // Dendro: printing temp variables

    // Dendro: printing variables
    //--
    a_rhs[pp] = -2 * K[pp] * alpha[pp] + lambda[0] * (beta0[pp] * grad_0_alpha +
                                                      beta1[pp] * grad_1_alpha +
                                                      beta2[pp] * grad_2_alpha);
    // Dendro: reduced ops: 9
    // Dendro: }}}
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&a_rhs[gidx], alpha[gidx], grad_0_alpha,
                                grad_1_alpha, grad_2_alpha, 1.0, 1.0, blk);
    }
    a_rhs[pp] += ko_sigma * (kograd_0_alpha + kograd_1_alpha + kograd_2_alpha);
    __syncthreads();
}

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, beta0, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_beta0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_beta0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_beta0 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_beta0 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_beta0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_beta0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_beta0 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_beta0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_beta0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_beta0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_beta0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_beta0 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, beta1, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_beta1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_beta1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_beta1 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_beta1 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_beta1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_beta1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_beta1 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_beta1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_beta1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_beta1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_beta1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_beta1 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, beta2, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_beta2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_beta2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_beta2 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_beta2 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_beta2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_beta2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_beta2 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_beta2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_beta2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_beta2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_beta2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_beta2 = Du[gidx];
__syncthreads();

{
    // Dendro: {{{
    // Dendro: original ops: 42
    // Dendro: printing temp variables
    const double DENDRO_0 =
        (3.0 / 4.0) * alpha[pp] * lambda_f[1] + (3.0 / 4.0) * lambda_f[0];

    // Dendro: printing variables
    //--
    b_rhs0[pp] = B0[pp] * DENDRO_0 + lambda[1] * (beta0[pp] * grad_0_beta0 +
                                                  beta1[pp] * grad_1_beta0 +
                                                  beta2[pp] * grad_2_beta0);
    //--
    b_rhs1[pp] = B1[pp] * DENDRO_0 + lambda[1] * (beta0[pp] * grad_0_beta1 +
                                                  beta1[pp] * grad_1_beta1 +
                                                  beta2[pp] * grad_2_beta1);
    //--
    b_rhs2[pp] = B2[pp] * DENDRO_0 + lambda[1] * (beta0[pp] * grad_0_beta2 +
                                                  beta1[pp] * grad_1_beta2 +
                                                  beta2[pp] * grad_2_beta2);
    // Dendro: reduced ops: 30
    // Dendro: }}}
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&b_rhs0[gidx], beta0[gidx], grad_0_beta0,
                                grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&b_rhs1[gidx], beta1[gidx], grad_0_beta1,
                                grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&b_rhs2[gidx], beta2[gidx], grad_0_beta2,
                                grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk);
    }
    b_rhs0[pp] += ko_sigma * (kograd_0_beta0 + kograd_1_beta0 + kograd_2_beta0);
    b_rhs1[pp] += ko_sigma * (kograd_0_beta1 + kograd_1_beta1 + kograd_2_beta1);
    b_rhs2[pp] += ko_sigma * (kograd_0_beta2 + kograd_1_beta2 + kograd_2_beta2);
    __syncthreads();
}

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, chi, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_chi = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_chi = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_chi = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_chi = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_chi = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_chi = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_chi = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_chi = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_chi = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_chi = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_chi = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_chi = Du[gidx];
__syncthreads();

{
    // Dendro: {{{
    // Dendro: original ops: 16
    // Dendro: printing temp variables
    const double DENDRO_0 = (2.0 / 3.0) * chi[pp];

    // Dendro: printing variables
    //--
    chi_rhs[pp]           = DENDRO_0 * K[pp] * alpha[pp] -
                  DENDRO_0 * (grad_0_beta0 + grad_1_beta1 + grad_2_beta2) +
                  beta0[pp] * grad_0_chi + beta1[pp] * grad_1_chi +
                  beta2[pp] * grad_2_chi;
    // Dendro: reduced ops: 14
    // Dendro: }}}
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&chi_rhs[gidx], chi[gidx], grad_0_chi,
                                grad_1_chi, grad_2_chi, 1.0, 1.0, blk);
    }
    chi_rhs[pp] += ko_sigma * (kograd_0_chi + kograd_1_chi + kograd_2_chi);
    __syncthreads();
}

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, gt0, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt0 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt0 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt0 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt0 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, gt1, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt1 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt1 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt1 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt1 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, gt2, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt2 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt2 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt2 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt2 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, gt3, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt3 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt3 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt3 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt3 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt3 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt3 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt3 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt3 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt3 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt3 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt3 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt3 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, gt4, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt4 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt4 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt4 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt4 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt4 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt4 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt4 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt4 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt4 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt4 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt4 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt4 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, gt5, blk);
device::__blk1_deriv644_xx<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt5 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt5 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt5 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt5 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt5 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt5 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(DDu, Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt5 = DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt5 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt5 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt5 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt5 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt5 = Du[gidx];
__syncthreads();

{
    // Dendro: {{{
    // Dendro: original ops: 156
    // Dendro: printing temp variables
    const double DENDRO_0  = 2 * alpha[pp];
    const double DENDRO_1  = 2 * gt1[pp];
    const double DENDRO_2  = 2 * gt2[pp];
    const double DENDRO_3  = (2.0 / 3.0) * gt0[pp];
    const double DENDRO_4  = (1.0 / 3.0) * gt1[pp];
    const double DENDRO_5  = (2.0 / 3.0) * grad_2_beta2;
    const double DENDRO_6  = (1.0 / 3.0) * gt2[pp];
    const double DENDRO_7  = (2.0 / 3.0) * grad_1_beta1;
    const double DENDRO_8  = (2.0 / 3.0) * grad_0_beta0;
    const double DENDRO_9  = 2 * gt4[pp];
    const double DENDRO_10 = (1.0 / 3.0) * gt4[pp];

    // Dendro: printing variables
    //--
    gt_rhs00[pp]           = -At0[pp] * DENDRO_0 + DENDRO_1 * grad_0_beta1 +
                   DENDRO_2 * grad_0_beta2 - DENDRO_3 * grad_1_beta1 -
                   DENDRO_3 * grad_2_beta2 + beta0[pp] * grad_0_gt0 +
                   beta1[pp] * grad_1_gt0 + beta2[pp] * grad_2_gt0 +
                   (4.0 / 3.0) * grad_0_beta0 * gt0[pp];
    //--
    gt_rhs01[pp] = -At1[pp] * DENDRO_0 + DENDRO_4 * grad_0_beta0 +
                   DENDRO_4 * grad_1_beta1 - DENDRO_5 * gt1[pp] +
                   beta0[pp] * grad_0_gt1 + beta1[pp] * grad_1_gt1 +
                   beta2[pp] * grad_2_gt1 + grad_0_beta1 * gt3[pp] +
                   grad_0_beta2 * gt4[pp] + grad_1_beta0 * gt0[pp] +
                   grad_1_beta2 * gt2[pp];
    //--
    gt_rhs02[pp] = -At2[pp] * DENDRO_0 + DENDRO_6 * grad_0_beta0 +
                   DENDRO_6 * grad_2_beta2 - DENDRO_7 * gt2[pp] +
                   beta0[pp] * grad_0_gt2 + beta1[pp] * grad_1_gt2 +
                   beta2[pp] * grad_2_gt2 + grad_0_beta1 * gt4[pp] +
                   grad_0_beta2 * gt5[pp] + grad_2_beta0 * gt0[pp] +
                   grad_2_beta1 * gt1[pp];
    //--
    gt_rhs11[pp] = -At3[pp] * DENDRO_0 + DENDRO_1 * grad_1_beta0 -
                   DENDRO_5 * gt3[pp] - DENDRO_8 * gt3[pp] +
                   DENDRO_9 * grad_1_beta2 + beta0[pp] * grad_0_gt3 +
                   beta1[pp] * grad_1_gt3 + beta2[pp] * grad_2_gt3 +
                   (4.0 / 3.0) * grad_1_beta1 * gt3[pp];
    //--
    gt_rhs12[pp] = -At4[pp] * DENDRO_0 + DENDRO_10 * grad_1_beta1 +
                   DENDRO_10 * grad_2_beta2 - DENDRO_8 * gt4[pp] +
                   beta0[pp] * grad_0_gt4 + beta1[pp] * grad_1_gt4 +
                   beta2[pp] * grad_2_gt4 + grad_1_beta0 * gt2[pp] +
                   grad_1_beta2 * gt5[pp] + grad_2_beta0 * gt1[pp] +
                   grad_2_beta1 * gt3[pp];
    //--
    gt_rhs22[pp] = -At5[pp] * DENDRO_0 + DENDRO_2 * grad_2_beta0 -
                   DENDRO_7 * gt5[pp] - DENDRO_8 * gt5[pp] +
                   DENDRO_9 * grad_2_beta1 + beta0[pp] * grad_0_gt5 +
                   beta1[pp] * grad_1_gt5 + beta2[pp] * grad_2_gt5 +
                   (4.0 / 3.0) * grad_2_beta2 * gt5[pp];
    // Dendro: reduced ops: 135
    // Dendro: }}}
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&gt_rhs00[gidx], gt0[gidx], grad_0_gt0,
                                grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk);
        radiative_bc_pt<pw, nx>(&gt_rhs01[gidx], gt1[gidx], grad_0_gt1,
                                grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&gt_rhs02[gidx], gt2[gidx], grad_0_gt2,
                                grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&gt_rhs11[gidx], gt3[gidx], grad_0_gt3,
                                grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk);
        radiative_bc_pt<pw, nx>(&gt_rhs12[gidx], gt4[gidx], grad_0_gt4,
                                grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&gt_rhs22[gidx], gt5[gidx], grad_0_gt5,
                                grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk);
    }
    gt_rhs00[pp] += ko_sigma * (kograd_0_gt0 + kograd_1_gt0 + kograd_2_gt0);
    gt_rhs01[pp] += ko_sigma * (kograd_0_gt1 + kograd_1_gt1 + kograd_2_gt1);
    gt_rhs02[pp] += ko_sigma * (kograd_0_gt2 + kograd_1_gt2 + kograd_2_gt2);
    gt_rhs11[pp] += ko_sigma * (kograd_0_gt3 + kograd_1_gt3 + kograd_2_gt3);
    gt_rhs12[pp] += ko_sigma * (kograd_0_gt4 + kograd_1_gt4 + kograd_2_gt4);
    gt_rhs22[pp] += ko_sigma * (kograd_0_gt5 + kograd_1_gt5 + kograd_2_gt5);
    __syncthreads();
}

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, K, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_K = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_K = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_K = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_K = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_K = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_K = Du[gidx];
__syncthreads();

{
    double DENDRO_igt0;
    double DENDRO_igt1;
    double DENDRO_igt2;
    double DENDRO_igt3;
    double DENDRO_igt4;
    double DENDRO_igt5;
    double DENDRO_C1_k0_0;
    double DENDRO_C1_k0_1;
    double DENDRO_C1_k0_2;
    double DENDRO_C1_k0_3;
    double DENDRO_C1_k0_4;
    double DENDRO_C1_k0_5;
    double DENDRO_C1_k1_0;
    double DENDRO_C1_k1_1;
    double DENDRO_C1_k1_2;
    double DENDRO_C1_k1_3;
    double DENDRO_C1_k1_4;
    double DENDRO_C1_k1_5;
    double DENDRO_C1_k2_0;
    double DENDRO_C1_k2_1;
    double DENDRO_C1_k2_2;
    double DENDRO_C1_k2_3;
    double DENDRO_C1_k2_4;
    double DENDRO_C1_k2_5;
    double DENDRO_C2_k0_0;
    double DENDRO_C2_k0_1;
    double DENDRO_C2_k0_2;
    double DENDRO_C2_k0_3;
    double DENDRO_C2_k0_4;
    double DENDRO_C2_k0_5;
    double DENDRO_C2_k1_0;
    double DENDRO_C2_k1_1;
    double DENDRO_C2_k1_2;
    double DENDRO_C2_k1_3;
    double DENDRO_C2_k1_4;
    double DENDRO_C2_k1_5;
    double DENDRO_C2_k2_0;
    double DENDRO_C2_k2_1;
    double DENDRO_C2_k2_2;
    double DENDRO_C2_k2_3;
    double DENDRO_C2_k2_4;
    double DENDRO_C2_k2_5;
    double DENDRO_C3_k0_0;
    double DENDRO_C3_k0_1;
    double DENDRO_C3_k0_2;
    double DENDRO_C3_k0_3;
    double DENDRO_C3_k0_4;
    double DENDRO_C3_k0_5;
    double DENDRO_C3_k1_0;
    double DENDRO_C3_k1_1;
    double DENDRO_C3_k1_2;
    double DENDRO_C3_k1_3;
    double DENDRO_C3_k1_4;
    double DENDRO_C3_k1_5;
    double DENDRO_C3_k2_0;
    double DENDRO_C3_k2_1;
    double DENDRO_C3_k2_2;
    double DENDRO_C3_k2_3;
    double DENDRO_C3_k2_4;
    double DENDRO_C3_k2_5;
    double DENDRO_RIJ0;
    double DENDRO_RIJ1;
    double DENDRO_RIJ2;
    double DENDRO_RIJ3;
    double DENDRO_RIJ4;
    double DENDRO_RIJ5;
    double DENDRO_Gtk0;
    double DENDRO_Gtk1;
    double DENDRO_Gtk2;
    {
        // Dendro: {{{
        // Dendro: original ops: 114
        // Dendro: printing temp variables
        const double DENDRO_0 = gt3[pp] * gt5[pp];
        const double DENDRO_1 = pow(gt4[pp], 2);
        const double DENDRO_2 = pow(gt1[pp], 2);
        const double DENDRO_3 = pow(gt2[pp], 2);
        const double DENDRO_4 = gt2[pp] * gt4[pp];
        const double DENDRO_5 =
            1.0 /
            (-DENDRO_0 * gt0[pp] + DENDRO_1 * gt0[pp] + DENDRO_2 * gt5[pp] +
             DENDRO_3 * gt3[pp] - 2 * DENDRO_4 * gt1[pp]);

        // Dendro: printing variables
        //--
        DENDRO_igt0 = -DENDRO_5 * (DENDRO_0 - DENDRO_1);
        //--
        DENDRO_igt1 = DENDRO_5 * (-DENDRO_4 + gt1[pp] * gt5[pp]);
        //--
        DENDRO_igt2 = -DENDRO_5 * (gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp]);
        //--
        DENDRO_igt3 = -DENDRO_5 * (-DENDRO_3 + gt0[pp] * gt5[pp]);
        //--
        DENDRO_igt4 = DENDRO_5 * (gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp]);
        //--
        DENDRO_igt5 = -DENDRO_5 * (-DENDRO_2 + gt0[pp] * gt3[pp]);
        // Dendro: reduced ops: 39
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k0_0 = 0.5 * grad_0_gt0;
        //--
        DENDRO_C1_k0_1 = 0.5 * grad_1_gt0;
        //--
        DENDRO_C1_k0_2 = 0.5 * grad_2_gt0;
        //--
        DENDRO_C1_k0_3 = -0.5 * grad_0_gt3 + 1.0 * grad_1_gt1;
        //--
        DENDRO_C1_k0_4 = 0.5 * (-grad_0_gt4 + grad_1_gt2 + grad_2_gt1);
        //--
        DENDRO_C1_k0_5 = -0.5 * grad_0_gt5 + 1.0 * grad_2_gt2;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k1_0 = 1.0 * grad_0_gt1 - 0.5 * grad_1_gt0;
        //--
        DENDRO_C1_k1_1 = 0.5 * grad_0_gt3;
        //--
        DENDRO_C1_k1_2 = 0.5 * (grad_0_gt4 - grad_1_gt2 + grad_2_gt1);
        //--
        DENDRO_C1_k1_3 = 0.5 * grad_1_gt3;
        //--
        DENDRO_C1_k1_4 = 0.5 * grad_2_gt3;
        //--
        DENDRO_C1_k1_5 = -0.5 * grad_1_gt5 + 1.0 * grad_2_gt4;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k2_0 = 1.0 * grad_0_gt2 - 0.5 * grad_2_gt0;
        //--
        DENDRO_C1_k2_1 = 0.5 * (grad_0_gt4 + grad_1_gt2 - grad_2_gt1);
        //--
        DENDRO_C1_k2_2 = 0.5 * grad_0_gt5;
        //--
        DENDRO_C1_k2_3 = 1.0 * grad_1_gt4 - 0.5 * grad_2_gt3;
        //--
        DENDRO_C1_k2_4 = 0.5 * grad_1_gt5;
        //--
        DENDRO_C1_k2_5 = 0.5 * grad_2_gt5;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k0_0 = DENDRO_C1_k0_0 * DENDRO_igt0 +
                         DENDRO_C1_k1_0 * DENDRO_igt1 +
                         DENDRO_C1_k2_0 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_1 = DENDRO_C1_k0_1 * DENDRO_igt0 +
                         DENDRO_C1_k1_1 * DENDRO_igt1 +
                         DENDRO_C1_k2_1 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_2 = DENDRO_C1_k0_2 * DENDRO_igt0 +
                         DENDRO_C1_k1_2 * DENDRO_igt1 +
                         DENDRO_C1_k2_2 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_3 = DENDRO_C1_k0_3 * DENDRO_igt0 +
                         DENDRO_C1_k1_3 * DENDRO_igt1 +
                         DENDRO_C1_k2_3 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_4 = DENDRO_C1_k0_4 * DENDRO_igt0 +
                         DENDRO_C1_k1_4 * DENDRO_igt1 +
                         DENDRO_C1_k2_4 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_5 = DENDRO_C1_k0_5 * DENDRO_igt0 +
                         DENDRO_C1_k1_5 * DENDRO_igt1 +
                         DENDRO_C1_k2_5 * DENDRO_igt2;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k1_0 = DENDRO_C1_k0_0 * DENDRO_igt1 +
                         DENDRO_C1_k1_0 * DENDRO_igt3 +
                         DENDRO_C1_k2_0 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_1 = DENDRO_C1_k0_1 * DENDRO_igt1 +
                         DENDRO_C1_k1_1 * DENDRO_igt3 +
                         DENDRO_C1_k2_1 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_2 = DENDRO_C1_k0_2 * DENDRO_igt1 +
                         DENDRO_C1_k1_2 * DENDRO_igt3 +
                         DENDRO_C1_k2_2 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_3 = DENDRO_C1_k0_3 * DENDRO_igt1 +
                         DENDRO_C1_k1_3 * DENDRO_igt3 +
                         DENDRO_C1_k2_3 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_4 = DENDRO_C1_k0_4 * DENDRO_igt1 +
                         DENDRO_C1_k1_4 * DENDRO_igt3 +
                         DENDRO_C1_k2_4 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_5 = DENDRO_C1_k0_5 * DENDRO_igt1 +
                         DENDRO_C1_k1_5 * DENDRO_igt3 +
                         DENDRO_C1_k2_5 * DENDRO_igt4;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k2_0 = DENDRO_C1_k0_0 * DENDRO_igt2 +
                         DENDRO_C1_k1_0 * DENDRO_igt4 +
                         DENDRO_C1_k2_0 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_1 = DENDRO_C1_k0_1 * DENDRO_igt2 +
                         DENDRO_C1_k1_1 * DENDRO_igt4 +
                         DENDRO_C1_k2_1 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_2 = DENDRO_C1_k0_2 * DENDRO_igt2 +
                         DENDRO_C1_k1_2 * DENDRO_igt4 +
                         DENDRO_C1_k2_2 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_3 = DENDRO_C1_k0_3 * DENDRO_igt2 +
                         DENDRO_C1_k1_3 * DENDRO_igt4 +
                         DENDRO_C1_k2_3 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_4 = DENDRO_C1_k0_4 * DENDRO_igt2 +
                         DENDRO_C1_k1_4 * DENDRO_igt4 +
                         DENDRO_C1_k2_4 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_5 = DENDRO_C1_k0_5 * DENDRO_igt2 +
                         DENDRO_C1_k1_5 * DENDRO_igt4 +
                         DENDRO_C1_k2_5 * DENDRO_igt5;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt0 * grad_0_chi +
                                0.5 * DENDRO_igt1 * grad_1_chi +
                                0.5 * DENDRO_igt2 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k0_0 = -DENDRO_0 * (-DENDRO_1 * gt0[pp] + 1.0 * grad_0_chi) +
                         DENDRO_C2_k0_0;
        //--
        DENDRO_C3_k0_1 = -DENDRO_0 * (-DENDRO_1 * gt1[pp] + 0.5 * grad_1_chi) +
                         DENDRO_C2_k0_1;
        //--
        DENDRO_C3_k0_2 = -DENDRO_0 * (-DENDRO_1 * gt2[pp] + 0.5 * grad_2_chi) +
                         DENDRO_C2_k0_2;
        //--
        DENDRO_C3_k0_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k0_3;
        //--
        DENDRO_C3_k0_4 = DENDRO_2 * gt4[pp] + DENDRO_C2_k0_4;
        //--
        DENDRO_C3_k0_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k0_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt1 * grad_0_chi +
                                0.5 * DENDRO_igt3 * grad_1_chi +
                                0.5 * DENDRO_igt4 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k1_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k1_0;
        //--
        DENDRO_C3_k1_1 = -DENDRO_0 * (-DENDRO_1 * gt1[pp] + 0.5 * grad_0_chi) +
                         DENDRO_C2_k1_1;
        //--
        DENDRO_C3_k1_2 = DENDRO_2 * gt2[pp] + DENDRO_C2_k1_2;
        //--
        DENDRO_C3_k1_3 = -DENDRO_0 * (-DENDRO_1 * gt3[pp] + 1.0 * grad_1_chi) +
                         DENDRO_C2_k1_3;
        //--
        DENDRO_C3_k1_4 = -DENDRO_0 * (-DENDRO_1 * gt4[pp] + 0.5 * grad_2_chi) +
                         DENDRO_C2_k1_4;
        //--
        DENDRO_C3_k1_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k1_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt2 * grad_0_chi +
                                0.5 * DENDRO_igt4 * grad_1_chi +
                                0.5 * DENDRO_igt5 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k2_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k2_0;
        //--
        DENDRO_C3_k2_1        = DENDRO_2 * gt1[pp] + DENDRO_C2_k2_1;
        //--
        DENDRO_C3_k2_2 = -DENDRO_0 * (-DENDRO_1 * gt2[pp] + 0.5 * grad_0_chi) +
                         DENDRO_C2_k2_2;
        //--
        DENDRO_C3_k2_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k2_3;
        //--
        DENDRO_C3_k2_4 = -DENDRO_0 * (-DENDRO_1 * gt4[pp] + 0.5 * grad_1_chi) +
                         DENDRO_C2_k2_4;
        //--
        DENDRO_C3_k2_5 = -DENDRO_0 * (-DENDRO_1 * gt5[pp] + 1.0 * grad_2_chi) +
                         DENDRO_C2_k2_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    // Dendro: {{{
    // Dendro: original ops: 333
    // Dendro: printing temp variables
    const double DENDRO_0 = pow(gt4[pp], 2);
    const double DENDRO_1 = pow(gt1[pp], 2);
    const double DENDRO_2 = pow(gt2[pp], 2);
    const double DENDRO_3 = gt0[pp] * gt3[pp];
    const double DENDRO_4 = gt1[pp] * gt4[pp];
    const double DENDRO_5 = chi[pp] / (DENDRO_0 * gt0[pp] + DENDRO_1 * gt5[pp] +
                                       DENDRO_2 * gt3[pp] - DENDRO_3 * gt5[pp] -
                                       2 * DENDRO_4 * gt2[pp]);
    const double DENDRO_6 = 2 * DENDRO_5;
    const double DENDRO_7 = (DENDRO_igt1 * DENDRO_igt1);
    const double DENDRO_8 = (DENDRO_igt2 * DENDRO_igt2);
    const double DENDRO_9 = At1[pp] * DENDRO_igt0;
    const double DENDRO_10 = At2[pp] * DENDRO_igt0;
    const double DENDRO_11 = DENDRO_igt1 * DENDRO_igt2;
    const double DENDRO_12 = 2 * At4[pp];
    const double DENDRO_13 = (DENDRO_igt4 * DENDRO_igt4);
    const double DENDRO_14 = At1[pp] * DENDRO_igt1;
    const double DENDRO_15 = At2[pp] * DENDRO_igt1;
    const double DENDRO_16 = 2 * DENDRO_igt4;
    const double DENDRO_17 = DENDRO_igt3 * DENDRO_igt4;
    const double DENDRO_18 = At1[pp] * DENDRO_igt2;
    const double DENDRO_19 = At2[pp] * DENDRO_igt2;
    const double DENDRO_20 = At4[pp] * DENDRO_igt5;
    const double DENDRO_21 = At0[pp] * DENDRO_igt0;
    const double DENDRO_22 = At3[pp] * DENDRO_igt1;
    const double DENDRO_23 = At4[pp] * DENDRO_igt1;
    const double DENDRO_24 = At4[pp] * DENDRO_igt2;
    const double DENDRO_25 = At5[pp] * DENDRO_igt2;

    // Dendro: printing variables
    //--
    K_rhs[pp] =
        -DENDRO_5 * (-DENDRO_0 + gt3[pp] * gt5[pp]) *
            (DENDRO_C3_k0_0 * grad_0_alpha + DENDRO_C3_k1_0 * grad_1_alpha +
             DENDRO_C3_k2_0 * grad_2_alpha - grad2_0_0_alpha) -
        DENDRO_5 * (-DENDRO_1 + DENDRO_3) *
            (DENDRO_C3_k0_5 * grad_0_alpha + DENDRO_C3_k1_5 * grad_1_alpha +
             DENDRO_C3_k2_5 * grad_2_alpha - grad2_2_2_alpha) -
        DENDRO_5 * (-DENDRO_2 + gt0[pp] * gt5[pp]) *
            (DENDRO_C3_k0_3 * grad_0_alpha + DENDRO_C3_k1_3 * grad_1_alpha +
             DENDRO_C3_k2_3 * grad_2_alpha - grad2_1_1_alpha) -
        DENDRO_6 * (DENDRO_4 - gt2[pp] * gt3[pp]) *
            (DENDRO_C3_k0_2 * grad_0_alpha + DENDRO_C3_k1_2 * grad_1_alpha +
             DENDRO_C3_k2_2 * grad_2_alpha - grad2_0_2_alpha) +
        DENDRO_6 * (gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp]) *
            (DENDRO_C3_k0_4 * grad_0_alpha + DENDRO_C3_k1_4 * grad_1_alpha +
             DENDRO_C3_k2_4 * grad_2_alpha - grad2_1_2_alpha) +
        DENDRO_6 * (gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp]) *
            (DENDRO_C3_k0_1 * grad_0_alpha + DENDRO_C3_k1_1 * grad_1_alpha +
             DENDRO_C3_k2_1 * grad_2_alpha - grad2_0_1_alpha) +
        (1.0 / 3.0) * alpha[pp] *
            (3 * At0[pp] *
                 (At0[pp] * (DENDRO_igt0 * DENDRO_igt0) + At3[pp] * DENDRO_7 +
                  At5[pp] * DENDRO_8 + 2 * DENDRO_10 * DENDRO_igt2 +
                  DENDRO_11 * DENDRO_12 + 2 * DENDRO_9 * DENDRO_igt1) +
             6 * At1[pp] *
                 (At1[pp] * DENDRO_7 + At2[pp] * DENDRO_11 +
                  DENDRO_10 * DENDRO_igt4 + DENDRO_21 * DENDRO_igt1 +
                  DENDRO_22 * DENDRO_igt3 + DENDRO_23 * DENDRO_igt4 +
                  DENDRO_24 * DENDRO_igt3 + DENDRO_25 * DENDRO_igt4 +
                  DENDRO_9 * DENDRO_igt3) +
             6 * At2[pp] *
                 (At1[pp] * DENDRO_11 + At2[pp] * DENDRO_8 +
                  DENDRO_10 * DENDRO_igt5 + DENDRO_21 * DENDRO_igt2 +
                  DENDRO_22 * DENDRO_igt4 + DENDRO_23 * DENDRO_igt5 +
                  DENDRO_24 * DENDRO_igt4 + DENDRO_25 * DENDRO_igt5 +
                  DENDRO_9 * DENDRO_igt4) +
             3 * At3[pp] *
                 (At0[pp] * DENDRO_7 + At3[pp] * (DENDRO_igt3 * DENDRO_igt3) +
                  At5[pp] * DENDRO_13 + DENDRO_12 * DENDRO_17 +
                  2 * DENDRO_14 * DENDRO_igt3 + DENDRO_15 * DENDRO_16) +
             6 * At4[pp] *
                 (At0[pp] * DENDRO_11 + At3[pp] * DENDRO_17 +
                  At4[pp] * DENDRO_13 + At5[pp] * DENDRO_igt4 * DENDRO_igt5 +
                  DENDRO_14 * DENDRO_igt4 + DENDRO_15 * DENDRO_igt5 +
                  DENDRO_18 * DENDRO_igt3 + DENDRO_19 * DENDRO_igt4 +
                  DENDRO_20 * DENDRO_igt3) +
             3 * At5[pp] *
                 (At0[pp] * DENDRO_8 + At3[pp] * DENDRO_13 +
                  At5[pp] * (DENDRO_igt5 * DENDRO_igt5) +
                  DENDRO_16 * DENDRO_18 + DENDRO_16 * DENDRO_20 +
                  2 * DENDRO_19 * DENDRO_igt5) +
             pow(K[pp], 2)) +
        beta0[pp] * grad_0_K + beta1[pp] * grad_1_K + beta2[pp] * grad_2_K;
    // Dendro: reduced ops: 222
    // Dendro: }}}
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&K_rhs[gidx], K[gidx], grad_0_K, grad_1_K,
                                grad_2_K, 1.0, 0.0, blk);
    }
    K_rhs[pp] += ko_sigma * (kograd_0_K + kograd_1_K + kograd_2_K);
    __syncthreads();
}

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, Gt0, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_Gt0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_Gt0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_Gt0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_Gt0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_Gt0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_Gt0 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, Gt1, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_Gt1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_Gt1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_Gt1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_Gt1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_Gt1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_Gt1 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, Gt2, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_Gt2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_Gt2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_Gt2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_Gt2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_Gt2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_Gt2 = Du[gidx];
__syncthreads();
device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, B0, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_B0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_B0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_B0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_B0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_B0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_B0 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, B1, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_B1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_B1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_B1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_B1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_B1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_B1 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, B2, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_B2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_B2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_B2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_B2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_B2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_B2 = Du[gidx];
__syncthreads();

{
    double DENDRO_igt0;
    double DENDRO_igt1;
    double DENDRO_igt2;
    double DENDRO_igt3;
    double DENDRO_igt4;
    double DENDRO_igt5;
    double DENDRO_C1_k0_0;
    double DENDRO_C1_k0_1;
    double DENDRO_C1_k0_2;
    double DENDRO_C1_k0_3;
    double DENDRO_C1_k0_4;
    double DENDRO_C1_k0_5;
    double DENDRO_C1_k1_0;
    double DENDRO_C1_k1_1;
    double DENDRO_C1_k1_2;
    double DENDRO_C1_k1_3;
    double DENDRO_C1_k1_4;
    double DENDRO_C1_k1_5;
    double DENDRO_C1_k2_0;
    double DENDRO_C1_k2_1;
    double DENDRO_C1_k2_2;
    double DENDRO_C1_k2_3;
    double DENDRO_C1_k2_4;
    double DENDRO_C1_k2_5;
    double DENDRO_C2_k0_0;
    double DENDRO_C2_k0_1;
    double DENDRO_C2_k0_2;
    double DENDRO_C2_k0_3;
    double DENDRO_C2_k0_4;
    double DENDRO_C2_k0_5;
    double DENDRO_C2_k1_0;
    double DENDRO_C2_k1_1;
    double DENDRO_C2_k1_2;
    double DENDRO_C2_k1_3;
    double DENDRO_C2_k1_4;
    double DENDRO_C2_k1_5;
    double DENDRO_C2_k2_0;
    double DENDRO_C2_k2_1;
    double DENDRO_C2_k2_2;
    double DENDRO_C2_k2_3;
    double DENDRO_C2_k2_4;
    double DENDRO_C2_k2_5;
    double DENDRO_C3_k0_0;
    double DENDRO_C3_k0_1;
    double DENDRO_C3_k0_2;
    double DENDRO_C3_k0_3;
    double DENDRO_C3_k0_4;
    double DENDRO_C3_k0_5;
    double DENDRO_C3_k1_0;
    double DENDRO_C3_k1_1;
    double DENDRO_C3_k1_2;
    double DENDRO_C3_k1_3;
    double DENDRO_C3_k1_4;
    double DENDRO_C3_k1_5;
    double DENDRO_C3_k2_0;
    double DENDRO_C3_k2_1;
    double DENDRO_C3_k2_2;
    double DENDRO_C3_k2_3;
    double DENDRO_C3_k2_4;
    double DENDRO_C3_k2_5;
    double DENDRO_RIJ0;
    double DENDRO_RIJ1;
    double DENDRO_RIJ2;
    double DENDRO_RIJ3;
    double DENDRO_RIJ4;
    double DENDRO_RIJ5;
    double DENDRO_Gtk0;
    double DENDRO_Gtk1;
    double DENDRO_Gtk2;
    {
        // Dendro: {{{
        // Dendro: original ops: 114
        // Dendro: printing temp variables
        const double DENDRO_0 = gt3[pp] * gt5[pp];
        const double DENDRO_1 = pow(gt4[pp], 2);
        const double DENDRO_2 = pow(gt1[pp], 2);
        const double DENDRO_3 = pow(gt2[pp], 2);
        const double DENDRO_4 = gt2[pp] * gt4[pp];
        const double DENDRO_5 =
            1.0 /
            (-DENDRO_0 * gt0[pp] + DENDRO_1 * gt0[pp] + DENDRO_2 * gt5[pp] +
             DENDRO_3 * gt3[pp] - 2 * DENDRO_4 * gt1[pp]);

        // Dendro: printing variables
        //--
        DENDRO_igt0 = -DENDRO_5 * (DENDRO_0 - DENDRO_1);
        //--
        DENDRO_igt1 = DENDRO_5 * (-DENDRO_4 + gt1[pp] * gt5[pp]);
        //--
        DENDRO_igt2 = -DENDRO_5 * (gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp]);
        //--
        DENDRO_igt3 = -DENDRO_5 * (-DENDRO_3 + gt0[pp] * gt5[pp]);
        //--
        DENDRO_igt4 = DENDRO_5 * (gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp]);
        //--
        DENDRO_igt5 = -DENDRO_5 * (-DENDRO_2 + gt0[pp] * gt3[pp]);
        // Dendro: reduced ops: 39
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k0_0 = 0.5 * grad_0_gt0;
        //--
        DENDRO_C1_k0_1 = 0.5 * grad_1_gt0;
        //--
        DENDRO_C1_k0_2 = 0.5 * grad_2_gt0;
        //--
        DENDRO_C1_k0_3 = -0.5 * grad_0_gt3 + 1.0 * grad_1_gt1;
        //--
        DENDRO_C1_k0_4 = 0.5 * (-grad_0_gt4 + grad_1_gt2 + grad_2_gt1);
        //--
        DENDRO_C1_k0_5 = -0.5 * grad_0_gt5 + 1.0 * grad_2_gt2;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k1_0 = 1.0 * grad_0_gt1 - 0.5 * grad_1_gt0;
        //--
        DENDRO_C1_k1_1 = 0.5 * grad_0_gt3;
        //--
        DENDRO_C1_k1_2 = 0.5 * (grad_0_gt4 - grad_1_gt2 + grad_2_gt1);
        //--
        DENDRO_C1_k1_3 = 0.5 * grad_1_gt3;
        //--
        DENDRO_C1_k1_4 = 0.5 * grad_2_gt3;
        //--
        DENDRO_C1_k1_5 = -0.5 * grad_1_gt5 + 1.0 * grad_2_gt4;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k2_0 = 1.0 * grad_0_gt2 - 0.5 * grad_2_gt0;
        //--
        DENDRO_C1_k2_1 = 0.5 * (grad_0_gt4 + grad_1_gt2 - grad_2_gt1);
        //--
        DENDRO_C1_k2_2 = 0.5 * grad_0_gt5;
        //--
        DENDRO_C1_k2_3 = 1.0 * grad_1_gt4 - 0.5 * grad_2_gt3;
        //--
        DENDRO_C1_k2_4 = 0.5 * grad_1_gt5;
        //--
        DENDRO_C1_k2_5 = 0.5 * grad_2_gt5;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k0_0 = DENDRO_C1_k0_0 * DENDRO_igt0 +
                         DENDRO_C1_k1_0 * DENDRO_igt1 +
                         DENDRO_C1_k2_0 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_1 = DENDRO_C1_k0_1 * DENDRO_igt0 +
                         DENDRO_C1_k1_1 * DENDRO_igt1 +
                         DENDRO_C1_k2_1 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_2 = DENDRO_C1_k0_2 * DENDRO_igt0 +
                         DENDRO_C1_k1_2 * DENDRO_igt1 +
                         DENDRO_C1_k2_2 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_3 = DENDRO_C1_k0_3 * DENDRO_igt0 +
                         DENDRO_C1_k1_3 * DENDRO_igt1 +
                         DENDRO_C1_k2_3 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_4 = DENDRO_C1_k0_4 * DENDRO_igt0 +
                         DENDRO_C1_k1_4 * DENDRO_igt1 +
                         DENDRO_C1_k2_4 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_5 = DENDRO_C1_k0_5 * DENDRO_igt0 +
                         DENDRO_C1_k1_5 * DENDRO_igt1 +
                         DENDRO_C1_k2_5 * DENDRO_igt2;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k1_0 = DENDRO_C1_k0_0 * DENDRO_igt1 +
                         DENDRO_C1_k1_0 * DENDRO_igt3 +
                         DENDRO_C1_k2_0 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_1 = DENDRO_C1_k0_1 * DENDRO_igt1 +
                         DENDRO_C1_k1_1 * DENDRO_igt3 +
                         DENDRO_C1_k2_1 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_2 = DENDRO_C1_k0_2 * DENDRO_igt1 +
                         DENDRO_C1_k1_2 * DENDRO_igt3 +
                         DENDRO_C1_k2_2 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_3 = DENDRO_C1_k0_3 * DENDRO_igt1 +
                         DENDRO_C1_k1_3 * DENDRO_igt3 +
                         DENDRO_C1_k2_3 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_4 = DENDRO_C1_k0_4 * DENDRO_igt1 +
                         DENDRO_C1_k1_4 * DENDRO_igt3 +
                         DENDRO_C1_k2_4 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_5 = DENDRO_C1_k0_5 * DENDRO_igt1 +
                         DENDRO_C1_k1_5 * DENDRO_igt3 +
                         DENDRO_C1_k2_5 * DENDRO_igt4;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k2_0 = DENDRO_C1_k0_0 * DENDRO_igt2 +
                         DENDRO_C1_k1_0 * DENDRO_igt4 +
                         DENDRO_C1_k2_0 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_1 = DENDRO_C1_k0_1 * DENDRO_igt2 +
                         DENDRO_C1_k1_1 * DENDRO_igt4 +
                         DENDRO_C1_k2_1 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_2 = DENDRO_C1_k0_2 * DENDRO_igt2 +
                         DENDRO_C1_k1_2 * DENDRO_igt4 +
                         DENDRO_C1_k2_2 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_3 = DENDRO_C1_k0_3 * DENDRO_igt2 +
                         DENDRO_C1_k1_3 * DENDRO_igt4 +
                         DENDRO_C1_k2_3 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_4 = DENDRO_C1_k0_4 * DENDRO_igt2 +
                         DENDRO_C1_k1_4 * DENDRO_igt4 +
                         DENDRO_C1_k2_4 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_5 = DENDRO_C1_k0_5 * DENDRO_igt2 +
                         DENDRO_C1_k1_5 * DENDRO_igt4 +
                         DENDRO_C1_k2_5 * DENDRO_igt5;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt0 * grad_0_chi +
                                0.5 * DENDRO_igt1 * grad_1_chi +
                                0.5 * DENDRO_igt2 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k0_0 = -DENDRO_0 * (-DENDRO_1 * gt0[pp] + 1.0 * grad_0_chi) +
                         DENDRO_C2_k0_0;
        //--
        DENDRO_C3_k0_1 = -DENDRO_0 * (-DENDRO_1 * gt1[pp] + 0.5 * grad_1_chi) +
                         DENDRO_C2_k0_1;
        //--
        DENDRO_C3_k0_2 = -DENDRO_0 * (-DENDRO_1 * gt2[pp] + 0.5 * grad_2_chi) +
                         DENDRO_C2_k0_2;
        //--
        DENDRO_C3_k0_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k0_3;
        //--
        DENDRO_C3_k0_4 = DENDRO_2 * gt4[pp] + DENDRO_C2_k0_4;
        //--
        DENDRO_C3_k0_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k0_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt1 * grad_0_chi +
                                0.5 * DENDRO_igt3 * grad_1_chi +
                                0.5 * DENDRO_igt4 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k1_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k1_0;
        //--
        DENDRO_C3_k1_1 = -DENDRO_0 * (-DENDRO_1 * gt1[pp] + 0.5 * grad_0_chi) +
                         DENDRO_C2_k1_1;
        //--
        DENDRO_C3_k1_2 = DENDRO_2 * gt2[pp] + DENDRO_C2_k1_2;
        //--
        DENDRO_C3_k1_3 = -DENDRO_0 * (-DENDRO_1 * gt3[pp] + 1.0 * grad_1_chi) +
                         DENDRO_C2_k1_3;
        //--
        DENDRO_C3_k1_4 = -DENDRO_0 * (-DENDRO_1 * gt4[pp] + 0.5 * grad_2_chi) +
                         DENDRO_C2_k1_4;
        //--
        DENDRO_C3_k1_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k1_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt2 * grad_0_chi +
                                0.5 * DENDRO_igt4 * grad_1_chi +
                                0.5 * DENDRO_igt5 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k2_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k2_0;
        //--
        DENDRO_C3_k2_1        = DENDRO_2 * gt1[pp] + DENDRO_C2_k2_1;
        //--
        DENDRO_C3_k2_2 = -DENDRO_0 * (-DENDRO_1 * gt2[pp] + 0.5 * grad_0_chi) +
                         DENDRO_C2_k2_2;
        //--
        DENDRO_C3_k2_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k2_3;
        //--
        DENDRO_C3_k2_4 = -DENDRO_0 * (-DENDRO_1 * gt4[pp] + 0.5 * grad_1_chi) +
                         DENDRO_C2_k2_4;
        //--
        DENDRO_C3_k2_5 = -DENDRO_0 * (-DENDRO_1 * gt5[pp] + 1.0 * grad_2_chi) +
                         DENDRO_C2_k2_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 42
        // Dendro: printing temp variables
        const double DENDRO_0 = 2 * DENDRO_igt1;
        const double DENDRO_1 = 2 * DENDRO_igt2;
        const double DENDRO_2 = 2 * DENDRO_igt4;

        // Dendro: printing variables
        //--
        DENDRO_Gtk0 = DENDRO_0 * DENDRO_C2_k0_1 + DENDRO_1 * DENDRO_C2_k0_2 +
                      DENDRO_2 * DENDRO_C2_k0_4 + DENDRO_C2_k0_0 * DENDRO_igt0 +
                      DENDRO_C2_k0_3 * DENDRO_igt3 +
                      DENDRO_C2_k0_5 * DENDRO_igt5;
        //--
        DENDRO_Gtk1 = DENDRO_0 * DENDRO_C2_k1_1 + DENDRO_1 * DENDRO_C2_k1_2 +
                      DENDRO_2 * DENDRO_C2_k1_4 + DENDRO_C2_k1_0 * DENDRO_igt0 +
                      DENDRO_C2_k1_3 * DENDRO_igt3 +
                      DENDRO_C2_k1_5 * DENDRO_igt5;
        //--
        DENDRO_Gtk2 = DENDRO_0 * DENDRO_C2_k2_1 + DENDRO_1 * DENDRO_C2_k2_2 +
                      DENDRO_2 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
                      DENDRO_C2_k2_3 * DENDRO_igt3 +
                      DENDRO_C2_k2_5 * DENDRO_igt5;
        // Dendro: reduced ops: 36
        // Dendro: }}}
    }
    // Dendro: {{{
    // Dendro: original ops: 2532
    // Dendro: printing temp variables
    const double DENDRO_0 = (DENDRO_igt1 * DENDRO_igt1);
    const double DENDRO_1 = (DENDRO_igt2 * DENDRO_igt2);
    const double DENDRO_2 = At1[pp] * DENDRO_igt0;
    const double DENDRO_3 = 2 * DENDRO_igt1;
    const double DENDRO_4 = At2[pp] * DENDRO_igt0;
    const double DENDRO_5 = 2 * DENDRO_igt2;
    const double DENDRO_6 = DENDRO_igt1 * DENDRO_igt2;
    const double DENDRO_7 = At0[pp] * (DENDRO_igt0 * DENDRO_igt0) +
                            At3[pp] * DENDRO_0 + 2 * At4[pp] * DENDRO_6 +
                            At5[pp] * DENDRO_1 + DENDRO_2 * DENDRO_3 +
                            DENDRO_4 * DENDRO_5;
    const double DENDRO_8  = 2 * grad_0_alpha;
    const double DENDRO_9  = At0[pp] * DENDRO_igt0;
    const double DENDRO_10 = At3[pp] * DENDRO_igt1;
    const double DENDRO_11 = At4[pp] * DENDRO_igt1;
    const double DENDRO_12 = At4[pp] * DENDRO_igt2;
    const double DENDRO_13 = At5[pp] * DENDRO_igt2;
    const double DENDRO_14 = At1[pp] * DENDRO_0 + At2[pp] * DENDRO_6 +
                             DENDRO_10 * DENDRO_igt3 + DENDRO_11 * DENDRO_igt4 +
                             DENDRO_12 * DENDRO_igt3 + DENDRO_13 * DENDRO_igt4 +
                             DENDRO_2 * DENDRO_igt3 + DENDRO_4 * DENDRO_igt4 +
                             DENDRO_9 * DENDRO_igt1;
    const double DENDRO_15 = 2 * grad_1_alpha;
    const double DENDRO_16 = At1[pp] * DENDRO_6 + At2[pp] * DENDRO_1 +
                             DENDRO_10 * DENDRO_igt4 + DENDRO_11 * DENDRO_igt5 +
                             DENDRO_12 * DENDRO_igt4 + DENDRO_13 * DENDRO_igt5 +
                             DENDRO_2 * DENDRO_igt4 + DENDRO_4 * DENDRO_igt5 +
                             DENDRO_9 * DENDRO_igt2;
    const double DENDRO_17 = 2 * grad_2_alpha;
    const double DENDRO_18 = 2 * DENDRO_igt4;
    const double DENDRO_19 = 4 * grad_0_K;
    const double DENDRO_20 = 9 / chi[pp];
    const double DENDRO_21 = DENDRO_20 * grad_0_chi;
    const double DENDRO_22 = (1.0 / 3.0) * alpha[pp];
    const double DENDRO_23 = 4 * grad_1_K;
    const double DENDRO_24 = DENDRO_20 * grad_1_chi;
    const double DENDRO_25 = 4 * grad_2_K;
    const double DENDRO_26 = DENDRO_20 * grad_2_chi;
    const double DENDRO_27 = (1.0 / 3.0) * DENDRO_igt0;
    const double DENDRO_28 = (1.0 / 3.0) * DENDRO_igt1;
    const double DENDRO_29 = (1.0 / 3.0) * DENDRO_igt2;
    const double DENDRO_30 = (2.0 / 3.0) * grad_0_beta0 +
                             (2.0 / 3.0) * grad_1_beta1 +
                             (2.0 / 3.0) * grad_2_beta2;
    const double DENDRO_31 = (7.0 / 3.0) * DENDRO_igt1;
    const double DENDRO_32 = (7.0 / 3.0) * DENDRO_igt2;
    const double DENDRO_33 = 2 * alpha[pp];
    const double DENDRO_34 = DENDRO_33 * DENDRO_7;
    const double DENDRO_35 = (DENDRO_igt4 * DENDRO_igt4);
    const double DENDRO_36 = At1[pp] * DENDRO_igt1;
    const double DENDRO_37 = At2[pp] * DENDRO_igt1;
    const double DENDRO_38 = At4[pp] * DENDRO_igt3;
    const double DENDRO_39 =
        At0[pp] * DENDRO_0 + At3[pp] * (DENDRO_igt3 * DENDRO_igt3) +
        At5[pp] * DENDRO_35 + DENDRO_18 * DENDRO_37 + DENDRO_18 * DENDRO_38 +
        2 * DENDRO_36 * DENDRO_igt3;
    const double DENDRO_40 = DENDRO_33 * DENDRO_39;
    const double DENDRO_41 = At1[pp] * DENDRO_igt2;
    const double DENDRO_42 = At2[pp] * DENDRO_igt2;
    const double DENDRO_43 = At0[pp] * DENDRO_1 + At3[pp] * DENDRO_35 +
                             At4[pp] * DENDRO_18 * DENDRO_igt5 +
                             At5[pp] * (DENDRO_igt5 * DENDRO_igt5) +
                             DENDRO_18 * DENDRO_41 +
                             2 * DENDRO_42 * DENDRO_igt5;
    const double DENDRO_44 = DENDRO_33 * DENDRO_43;
    const double DENDRO_45 = 4 * alpha[pp];
    const double DENDRO_46 = DENDRO_14 * DENDRO_45;
    const double DENDRO_47 = DENDRO_16 * DENDRO_45;
    const double DENDRO_48 =
        At0[pp] * DENDRO_6 + At3[pp] * DENDRO_igt3 * DENDRO_igt4 +
        At4[pp] * DENDRO_35 + At5[pp] * DENDRO_igt4 * DENDRO_igt5 +
        DENDRO_36 * DENDRO_igt4 + DENDRO_37 * DENDRO_igt5 +
        DENDRO_38 * DENDRO_igt5 + DENDRO_41 * DENDRO_igt3 +
        DENDRO_42 * DENDRO_igt4;
    const double DENDRO_49 = DENDRO_45 * DENDRO_48;
    const double DENDRO_50 = beta0[pp] * grad_0_Gt0 + beta1[pp] * grad_1_Gt0 +
                             beta2[pp] * grad_2_Gt0;
    const double DENDRO_51 =
        -DENDRO_14 * DENDRO_15 - DENDRO_16 * DENDRO_17 +
        DENDRO_18 * grad2_1_2_beta0 -
        DENDRO_22 * (DENDRO_14 * DENDRO_24 + DENDRO_23 * DENDRO_igt1) -
        DENDRO_22 * (DENDRO_16 * DENDRO_26 + DENDRO_25 * DENDRO_igt2) -
        DENDRO_22 * (DENDRO_19 * DENDRO_igt0 + DENDRO_21 * DENDRO_7) +
        DENDRO_27 * grad2_0_1_beta1 + DENDRO_27 * grad2_0_2_beta2 +
        DENDRO_28 * grad2_1_1_beta1 + DENDRO_28 * grad2_1_2_beta2 +
        DENDRO_29 * grad2_1_2_beta1 + DENDRO_29 * grad2_2_2_beta2 +
        DENDRO_30 * DENDRO_Gtk0 + DENDRO_31 * grad2_0_1_beta0 +
        DENDRO_32 * grad2_0_2_beta0 + DENDRO_34 * DENDRO_C2_k0_0 +
        DENDRO_40 * DENDRO_C2_k0_3 + DENDRO_44 * DENDRO_C2_k0_5 +
        DENDRO_46 * DENDRO_C2_k0_1 + DENDRO_47 * DENDRO_C2_k0_2 +
        DENDRO_49 * DENDRO_C2_k0_4 + DENDRO_50 - DENDRO_7 * DENDRO_8 -
        DENDRO_Gtk0 * grad_0_beta0 - DENDRO_Gtk1 * grad_1_beta0 -
        DENDRO_Gtk2 * grad_2_beta0 +
        (4.0 / 3.0) * DENDRO_igt0 * grad2_0_0_beta0 +
        DENDRO_igt3 * grad2_1_1_beta0 + DENDRO_igt5 * grad2_2_2_beta0;
    const double DENDRO_52 = (1.0 / 3.0) * DENDRO_igt3;
    const double DENDRO_53 = (1.0 / 3.0) * DENDRO_igt4;
    const double DENDRO_54 = (7.0 / 3.0) * DENDRO_igt4;
    const double DENDRO_55 = beta0[pp] * grad_0_Gt1 + beta1[pp] * grad_1_Gt1 +
                             beta2[pp] * grad_2_Gt1;
    const double DENDRO_56 =
        -DENDRO_14 * DENDRO_8 - DENDRO_15 * DENDRO_39 - DENDRO_17 * DENDRO_48 -
        DENDRO_22 * (DENDRO_14 * DENDRO_21 + DENDRO_19 * DENDRO_igt1) -
        DENDRO_22 * (DENDRO_23 * DENDRO_igt3 + DENDRO_24 * DENDRO_39) -
        DENDRO_22 * (DENDRO_25 * DENDRO_igt4 + DENDRO_26 * DENDRO_48) +
        DENDRO_28 * grad2_0_0_beta0 + DENDRO_28 * grad2_0_2_beta2 +
        DENDRO_30 * DENDRO_Gtk1 + DENDRO_31 * grad2_0_1_beta1 +
        DENDRO_34 * DENDRO_C2_k1_0 + DENDRO_40 * DENDRO_C2_k1_3 +
        DENDRO_44 * DENDRO_C2_k1_5 + DENDRO_46 * DENDRO_C2_k1_1 +
        DENDRO_47 * DENDRO_C2_k1_2 + DENDRO_49 * DENDRO_C2_k1_4 +
        DENDRO_5 * grad2_0_2_beta1 + DENDRO_52 * grad2_0_1_beta0 +
        DENDRO_52 * grad2_1_2_beta2 + DENDRO_53 * grad2_0_2_beta0 +
        DENDRO_53 * grad2_2_2_beta2 + DENDRO_54 * grad2_1_2_beta1 + DENDRO_55 -
        DENDRO_Gtk0 * grad_0_beta1 - DENDRO_Gtk1 * grad_1_beta1 -
        DENDRO_Gtk2 * grad_2_beta1 + DENDRO_igt0 * grad2_0_0_beta1 +
        (4.0 / 3.0) * DENDRO_igt3 * grad2_1_1_beta1 +
        DENDRO_igt5 * grad2_2_2_beta1;
    const double DENDRO_57 = (1.0 / 3.0) * DENDRO_igt5;
    const double DENDRO_58 = beta0[pp] * grad_0_Gt2 + beta1[pp] * grad_1_Gt2 +
                             beta2[pp] * grad_2_Gt2;
    const double DENDRO_59 =
        -DENDRO_15 * DENDRO_48 - DENDRO_16 * DENDRO_8 - DENDRO_17 * DENDRO_43 -
        DENDRO_22 * (DENDRO_16 * DENDRO_21 + DENDRO_19 * DENDRO_igt2) -
        DENDRO_22 * (DENDRO_23 * DENDRO_igt4 + DENDRO_24 * DENDRO_48) -
        DENDRO_22 * (DENDRO_25 * DENDRO_igt5 + DENDRO_26 * DENDRO_43) +
        DENDRO_29 * grad2_0_0_beta0 + DENDRO_29 * grad2_0_1_beta1 +
        DENDRO_3 * grad2_0_1_beta2 + DENDRO_30 * DENDRO_Gtk2 +
        DENDRO_32 * grad2_0_2_beta2 + DENDRO_34 * DENDRO_C2_k2_0 +
        DENDRO_40 * DENDRO_C2_k2_3 + DENDRO_44 * DENDRO_C2_k2_5 +
        DENDRO_46 * DENDRO_C2_k2_1 + DENDRO_47 * DENDRO_C2_k2_2 +
        DENDRO_49 * DENDRO_C2_k2_4 + DENDRO_53 * grad2_0_1_beta0 +
        DENDRO_53 * grad2_1_1_beta1 + DENDRO_54 * grad2_1_2_beta2 +
        DENDRO_57 * grad2_0_2_beta0 + DENDRO_57 * grad2_1_2_beta1 + DENDRO_58 -
        DENDRO_Gtk0 * grad_0_beta2 - DENDRO_Gtk1 * grad_1_beta2 -
        DENDRO_Gtk2 * grad_2_beta2 + DENDRO_igt0 * grad2_0_0_beta2 +
        DENDRO_igt3 * grad2_1_1_beta2 +
        (4.0 / 3.0) * DENDRO_igt5 * grad2_2_2_beta2;

    // Dendro: printing variables
    //--
    Gt_rhs0[pp] = DENDRO_51;
    //--
    Gt_rhs1[pp] = DENDRO_56;
    //--
    Gt_rhs2[pp] = DENDRO_59;
    //--
    B_rhs0[pp]  = -B0[pp] * eta - DENDRO_50 * lambda[3] + DENDRO_51 +
                 lambda[2] * (beta0[pp] * grad_0_B0 + beta1[pp] * grad_1_B0 +
                              beta2[pp] * grad_2_B0);
    //--
    B_rhs1[pp] = -B1[pp] * eta - DENDRO_55 * lambda[3] + DENDRO_56 +
                 lambda[2] * (beta0[pp] * grad_0_B1 + beta1[pp] * grad_1_B1 +
                              beta2[pp] * grad_2_B1);
    //--
    B_rhs2[pp] = -B2[pp] * eta - DENDRO_58 * lambda[3] + DENDRO_59 +
                 lambda[2] * (beta0[pp] * grad_0_B2 + beta1[pp] * grad_1_B2 +
                              beta2[pp] * grad_2_B2);
    // Dendro: reduced ops: 400
    // Dendro: }}}
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&Gt_rhs0[gidx], Gt0[gidx], grad_0_Gt0,
                                grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&Gt_rhs1[gidx], Gt1[gidx], grad_0_Gt1,
                                grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&Gt_rhs2[gidx], Gt2[gidx], grad_0_Gt2,
                                grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk);
    }
    Gt_rhs0[pp] += ko_sigma * (kograd_0_Gt0 + kograd_1_Gt0 + kograd_2_Gt0);
    Gt_rhs1[pp] += ko_sigma * (kograd_0_Gt1 + kograd_1_Gt1 + kograd_2_Gt1);
    Gt_rhs2[pp] += ko_sigma * (kograd_0_Gt2 + kograd_1_Gt2 + kograd_2_Gt2);
    __syncthreads();
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&B_rhs0[gidx], B0[gidx], grad_0_B0, grad_1_B0,
                                grad_2_B0, 1.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&B_rhs1[gidx], B1[gidx], grad_0_B1, grad_1_B1,
                                grad_2_B1, 1.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&B_rhs2[gidx], B2[gidx], grad_0_B2, grad_1_B2,
                                grad_2_B2, 1.0, 0.0, blk);
    }
    B_rhs0[pp] += ko_sigma * (kograd_0_B0 + kograd_1_B0 + kograd_2_B0);
    B_rhs1[pp] += ko_sigma * (kograd_0_B1 + kograd_1_B1 + kograd_2_B1);
    B_rhs2[pp] += ko_sigma * (kograd_0_B2 + kograd_1_B2 + kograd_2_B2);
    __syncthreads();
}

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, At0, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At0 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At0 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At0 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, At1, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At1 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At1 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At1 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, At2, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At2 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At2 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At2 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, At3, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At3 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At3 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At3 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At3 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At3 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At3 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, At4, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At4 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At4 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At4 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At4 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At4 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At4 = Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL, pw, nx>(su, At5, blk);
device::__blk1_deriv644_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At5 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At5 = Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At5 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_x<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At5 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At5 = Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw, pencils, pencil_sz>(Du, su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At5 = Du[gidx];
__syncthreads();

{
    double DENDRO_igt0;
    double DENDRO_igt1;
    double DENDRO_igt2;
    double DENDRO_igt3;
    double DENDRO_igt4;
    double DENDRO_igt5;
    double DENDRO_C1_k0_0;
    double DENDRO_C1_k0_1;
    double DENDRO_C1_k0_2;
    double DENDRO_C1_k0_3;
    double DENDRO_C1_k0_4;
    double DENDRO_C1_k0_5;
    double DENDRO_C1_k1_0;
    double DENDRO_C1_k1_1;
    double DENDRO_C1_k1_2;
    double DENDRO_C1_k1_3;
    double DENDRO_C1_k1_4;
    double DENDRO_C1_k1_5;
    double DENDRO_C1_k2_0;
    double DENDRO_C1_k2_1;
    double DENDRO_C1_k2_2;
    double DENDRO_C1_k2_3;
    double DENDRO_C1_k2_4;
    double DENDRO_C1_k2_5;
    double DENDRO_C2_k0_0;
    double DENDRO_C2_k0_1;
    double DENDRO_C2_k0_2;
    double DENDRO_C2_k0_3;
    double DENDRO_C2_k0_4;
    double DENDRO_C2_k0_5;
    double DENDRO_C2_k1_0;
    double DENDRO_C2_k1_1;
    double DENDRO_C2_k1_2;
    double DENDRO_C2_k1_3;
    double DENDRO_C2_k1_4;
    double DENDRO_C2_k1_5;
    double DENDRO_C2_k2_0;
    double DENDRO_C2_k2_1;
    double DENDRO_C2_k2_2;
    double DENDRO_C2_k2_3;
    double DENDRO_C2_k2_4;
    double DENDRO_C2_k2_5;
    double DENDRO_C3_k0_0;
    double DENDRO_C3_k0_1;
    double DENDRO_C3_k0_2;
    double DENDRO_C3_k0_3;
    double DENDRO_C3_k0_4;
    double DENDRO_C3_k0_5;
    double DENDRO_C3_k1_0;
    double DENDRO_C3_k1_1;
    double DENDRO_C3_k1_2;
    double DENDRO_C3_k1_3;
    double DENDRO_C3_k1_4;
    double DENDRO_C3_k1_5;
    double DENDRO_C3_k2_0;
    double DENDRO_C3_k2_1;
    double DENDRO_C3_k2_2;
    double DENDRO_C3_k2_3;
    double DENDRO_C3_k2_4;
    double DENDRO_C3_k2_5;
    double DENDRO_RIJ0;
    double DENDRO_RIJ1;
    double DENDRO_RIJ2;
    double DENDRO_RIJ3;
    double DENDRO_RIJ4;
    double DENDRO_RIJ5;
    double DENDRO_Gtk0;
    double DENDRO_Gtk1;
    double DENDRO_Gtk2;
    {
        // Dendro: {{{
        // Dendro: original ops: 114
        // Dendro: printing temp variables
        const double DENDRO_0 = gt3[pp] * gt5[pp];
        const double DENDRO_1 = pow(gt4[pp], 2);
        const double DENDRO_2 = pow(gt1[pp], 2);
        const double DENDRO_3 = pow(gt2[pp], 2);
        const double DENDRO_4 = gt2[pp] * gt4[pp];
        const double DENDRO_5 =
            1.0 /
            (-DENDRO_0 * gt0[pp] + DENDRO_1 * gt0[pp] + DENDRO_2 * gt5[pp] +
             DENDRO_3 * gt3[pp] - 2 * DENDRO_4 * gt1[pp]);

        // Dendro: printing variables
        //--
        DENDRO_igt0 = -DENDRO_5 * (DENDRO_0 - DENDRO_1);
        //--
        DENDRO_igt1 = DENDRO_5 * (-DENDRO_4 + gt1[pp] * gt5[pp]);
        //--
        DENDRO_igt2 = -DENDRO_5 * (gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp]);
        //--
        DENDRO_igt3 = -DENDRO_5 * (-DENDRO_3 + gt0[pp] * gt5[pp]);
        //--
        DENDRO_igt4 = DENDRO_5 * (gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp]);
        //--
        DENDRO_igt5 = -DENDRO_5 * (-DENDRO_2 + gt0[pp] * gt3[pp]);
        // Dendro: reduced ops: 39
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k0_0 = 0.5 * grad_0_gt0;
        //--
        DENDRO_C1_k0_1 = 0.5 * grad_1_gt0;
        //--
        DENDRO_C1_k0_2 = 0.5 * grad_2_gt0;
        //--
        DENDRO_C1_k0_3 = -0.5 * grad_0_gt3 + 1.0 * grad_1_gt1;
        //--
        DENDRO_C1_k0_4 = 0.5 * (-grad_0_gt4 + grad_1_gt2 + grad_2_gt1);
        //--
        DENDRO_C1_k0_5 = -0.5 * grad_0_gt5 + 1.0 * grad_2_gt2;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k1_0 = 1.0 * grad_0_gt1 - 0.5 * grad_1_gt0;
        //--
        DENDRO_C1_k1_1 = 0.5 * grad_0_gt3;
        //--
        DENDRO_C1_k1_2 = 0.5 * (grad_0_gt4 - grad_1_gt2 + grad_2_gt1);
        //--
        DENDRO_C1_k1_3 = 0.5 * grad_1_gt3;
        //--
        DENDRO_C1_k1_4 = 0.5 * grad_2_gt3;
        //--
        DENDRO_C1_k1_5 = -0.5 * grad_1_gt5 + 1.0 * grad_2_gt4;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 14
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C1_k2_0 = 1.0 * grad_0_gt2 - 0.5 * grad_2_gt0;
        //--
        DENDRO_C1_k2_1 = 0.5 * (grad_0_gt4 + grad_1_gt2 - grad_2_gt1);
        //--
        DENDRO_C1_k2_2 = 0.5 * grad_0_gt5;
        //--
        DENDRO_C1_k2_3 = 1.0 * grad_1_gt4 - 0.5 * grad_2_gt3;
        //--
        DENDRO_C1_k2_4 = 0.5 * grad_1_gt5;
        //--
        DENDRO_C1_k2_5 = 0.5 * grad_2_gt5;
        // Dendro: reduced ops: 12
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k0_0 = DENDRO_C1_k0_0 * DENDRO_igt0 +
                         DENDRO_C1_k1_0 * DENDRO_igt1 +
                         DENDRO_C1_k2_0 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_1 = DENDRO_C1_k0_1 * DENDRO_igt0 +
                         DENDRO_C1_k1_1 * DENDRO_igt1 +
                         DENDRO_C1_k2_1 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_2 = DENDRO_C1_k0_2 * DENDRO_igt0 +
                         DENDRO_C1_k1_2 * DENDRO_igt1 +
                         DENDRO_C1_k2_2 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_3 = DENDRO_C1_k0_3 * DENDRO_igt0 +
                         DENDRO_C1_k1_3 * DENDRO_igt1 +
                         DENDRO_C1_k2_3 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_4 = DENDRO_C1_k0_4 * DENDRO_igt0 +
                         DENDRO_C1_k1_4 * DENDRO_igt1 +
                         DENDRO_C1_k2_4 * DENDRO_igt2;
        //--
        DENDRO_C2_k0_5 = DENDRO_C1_k0_5 * DENDRO_igt0 +
                         DENDRO_C1_k1_5 * DENDRO_igt1 +
                         DENDRO_C1_k2_5 * DENDRO_igt2;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k1_0 = DENDRO_C1_k0_0 * DENDRO_igt1 +
                         DENDRO_C1_k1_0 * DENDRO_igt3 +
                         DENDRO_C1_k2_0 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_1 = DENDRO_C1_k0_1 * DENDRO_igt1 +
                         DENDRO_C1_k1_1 * DENDRO_igt3 +
                         DENDRO_C1_k2_1 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_2 = DENDRO_C1_k0_2 * DENDRO_igt1 +
                         DENDRO_C1_k1_2 * DENDRO_igt3 +
                         DENDRO_C1_k2_2 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_3 = DENDRO_C1_k0_3 * DENDRO_igt1 +
                         DENDRO_C1_k1_3 * DENDRO_igt3 +
                         DENDRO_C1_k2_3 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_4 = DENDRO_C1_k0_4 * DENDRO_igt1 +
                         DENDRO_C1_k1_4 * DENDRO_igt3 +
                         DENDRO_C1_k2_4 * DENDRO_igt4;
        //--
        DENDRO_C2_k1_5 = DENDRO_C1_k0_5 * DENDRO_igt1 +
                         DENDRO_C1_k1_5 * DENDRO_igt3 +
                         DENDRO_C1_k2_5 * DENDRO_igt4;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 30
        // Dendro: printing temp variables

        // Dendro: printing variables
        //--
        DENDRO_C2_k2_0 = DENDRO_C1_k0_0 * DENDRO_igt2 +
                         DENDRO_C1_k1_0 * DENDRO_igt4 +
                         DENDRO_C1_k2_0 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_1 = DENDRO_C1_k0_1 * DENDRO_igt2 +
                         DENDRO_C1_k1_1 * DENDRO_igt4 +
                         DENDRO_C1_k2_1 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_2 = DENDRO_C1_k0_2 * DENDRO_igt2 +
                         DENDRO_C1_k1_2 * DENDRO_igt4 +
                         DENDRO_C1_k2_2 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_3 = DENDRO_C1_k0_3 * DENDRO_igt2 +
                         DENDRO_C1_k1_3 * DENDRO_igt4 +
                         DENDRO_C1_k2_3 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_4 = DENDRO_C1_k0_4 * DENDRO_igt2 +
                         DENDRO_C1_k1_4 * DENDRO_igt4 +
                         DENDRO_C1_k2_4 * DENDRO_igt5;
        //--
        DENDRO_C2_k2_5 = DENDRO_C1_k0_5 * DENDRO_igt2 +
                         DENDRO_C1_k1_5 * DENDRO_igt4 +
                         DENDRO_C1_k2_5 * DENDRO_igt5;
        // Dendro: reduced ops: 30
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt0 * grad_0_chi +
                                0.5 * DENDRO_igt1 * grad_1_chi +
                                0.5 * DENDRO_igt2 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k0_0 = -DENDRO_0 * (-DENDRO_1 * gt0[pp] + 1.0 * grad_0_chi) +
                         DENDRO_C2_k0_0;
        //--
        DENDRO_C3_k0_1 = -DENDRO_0 * (-DENDRO_1 * gt1[pp] + 0.5 * grad_1_chi) +
                         DENDRO_C2_k0_1;
        //--
        DENDRO_C3_k0_2 = -DENDRO_0 * (-DENDRO_1 * gt2[pp] + 0.5 * grad_2_chi) +
                         DENDRO_C2_k0_2;
        //--
        DENDRO_C3_k0_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k0_3;
        //--
        DENDRO_C3_k0_4 = DENDRO_2 * gt4[pp] + DENDRO_C2_k0_4;
        //--
        DENDRO_C3_k0_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k0_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt1 * grad_0_chi +
                                0.5 * DENDRO_igt3 * grad_1_chi +
                                0.5 * DENDRO_igt4 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k1_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k1_0;
        //--
        DENDRO_C3_k1_1 = -DENDRO_0 * (-DENDRO_1 * gt1[pp] + 0.5 * grad_0_chi) +
                         DENDRO_C2_k1_1;
        //--
        DENDRO_C3_k1_2 = DENDRO_2 * gt2[pp] + DENDRO_C2_k1_2;
        //--
        DENDRO_C3_k1_3 = -DENDRO_0 * (-DENDRO_1 * gt3[pp] + 1.0 * grad_1_chi) +
                         DENDRO_C2_k1_3;
        //--
        DENDRO_C3_k1_4 = -DENDRO_0 * (-DENDRO_1 * gt4[pp] + 0.5 * grad_2_chi) +
                         DENDRO_C2_k1_4;
        //--
        DENDRO_C3_k1_5 = DENDRO_2 * gt5[pp] + DENDRO_C2_k1_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 60
        // Dendro: printing temp variables
        const double DENDRO_0 = 1.0 / chi[pp];
        const double DENDRO_1 = 0.5 * DENDRO_igt2 * grad_0_chi +
                                0.5 * DENDRO_igt4 * grad_1_chi +
                                0.5 * DENDRO_igt5 * grad_2_chi;
        const double DENDRO_2 = DENDRO_0 * DENDRO_1;

        // Dendro: printing variables
        //--
        DENDRO_C3_k2_0        = DENDRO_2 * gt0[pp] + DENDRO_C2_k2_0;
        //--
        DENDRO_C3_k2_1        = DENDRO_2 * gt1[pp] + DENDRO_C2_k2_1;
        //--
        DENDRO_C3_k2_2 = -DENDRO_0 * (-DENDRO_1 * gt2[pp] + 0.5 * grad_0_chi) +
                         DENDRO_C2_k2_2;
        //--
        DENDRO_C3_k2_3 = DENDRO_2 * gt3[pp] + DENDRO_C2_k2_3;
        //--
        DENDRO_C3_k2_4 = -DENDRO_0 * (-DENDRO_1 * gt4[pp] + 0.5 * grad_1_chi) +
                         DENDRO_C2_k2_4;
        //--
        DENDRO_C3_k2_5 = -DENDRO_0 * (-DENDRO_1 * gt5[pp] + 1.0 * grad_2_chi) +
                         DENDRO_C2_k2_5;
        // Dendro: reduced ops: 31
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 363
        // Dendro: printing temp variables
        const double DENDRO_0 = pow(chi[pp], 2);
        const double DENDRO_1 = (grad_0_chi * grad_0_chi);
        const double DENDRO_2 = 2 * chi[pp];
        const double DENDRO_3 = 3 * grad_0_chi;
        const double DENDRO_4 = 2 * DENDRO_igt1;
        const double DENDRO_5 = 2 * DENDRO_igt2;
        const double DENDRO_6 = 2 * DENDRO_igt4;
        const double DENDRO_7 = DENDRO_C2_k0_0 * DENDRO_igt0;
        const double DENDRO_8 =
            DENDRO_4 * DENDRO_C2_k0_1 + DENDRO_5 * DENDRO_C2_k0_2 +
            DENDRO_6 * DENDRO_C2_k0_4 + DENDRO_7 +
            DENDRO_C2_k0_3 * DENDRO_igt3 + DENDRO_C2_k0_5 * DENDRO_igt5;
        const double DENDRO_9 = DENDRO_C2_k1_0 * DENDRO_igt0;
        const double DENDRO_10 =
            DENDRO_4 * DENDRO_C2_k1_1 + DENDRO_5 * DENDRO_C2_k1_2 +
            DENDRO_6 * DENDRO_C2_k1_4 + DENDRO_9 +
            DENDRO_C2_k1_3 * DENDRO_igt3 + DENDRO_C2_k1_5 * DENDRO_igt5;
        const double DENDRO_11 = DENDRO_C2_k2_0 * DENDRO_igt0;
        const double DENDRO_12 =
            DENDRO_11 + DENDRO_4 * DENDRO_C2_k2_1 + DENDRO_5 * DENDRO_C2_k2_2 +
            DENDRO_6 * DENDRO_C2_k2_4 + DENDRO_C2_k2_3 * DENDRO_igt3 +
            DENDRO_C2_k2_5 * DENDRO_igt5;
        const double DENDRO_13 = 2 * DENDRO_C1_k0_1;
        const double DENDRO_14 = 2 * DENDRO_C1_k0_3;
        const double DENDRO_15 = 4 * DENDRO_igt3;
        const double DENDRO_16 = 2 * DENDRO_C1_k0_4;
        const double DENDRO_17 = 4 * DENDRO_igt5;
        const double DENDRO_18 = 2 * DENDRO_C1_k0_2;
        const double DENDRO_19 = 2 * DENDRO_C1_k0_5;
        const double DENDRO_20 = DENDRO_C1_k0_0 * DENDRO_C2_k0_1;
        const double DENDRO_21 = DENDRO_C1_k0_1 * DENDRO_C2_k0_0;
        const double DENDRO_22 = 4 * DENDRO_igt1;
        const double DENDRO_23 = DENDRO_C1_k0_0 * DENDRO_C2_k0_2;
        const double DENDRO_24 = DENDRO_C1_k0_2 * DENDRO_C2_k0_0;
        const double DENDRO_25 = 4 * DENDRO_igt2;
        const double DENDRO_26 = DENDRO_C1_k0_1 * DENDRO_C2_k0_2;
        const double DENDRO_27 = DENDRO_C1_k0_2 * DENDRO_C2_k0_1;
        const double DENDRO_28 = 4 * DENDRO_igt4;

        // Dendro: printing variables
        //--
        DENDRO_RIJ0 =
            (1.0 / 4.0) *
            (DENDRO_0 *
                 (4.0 * DENDRO_10 * DENDRO_C1_k0_1 +
                  4 * DENDRO_11 * (DENDRO_18 + DENDRO_C1_k2_0) +
                  4.0 * DENDRO_12 * DENDRO_C1_k0_2 +
                  DENDRO_15 * DENDRO_C2_k1_1 * (DENDRO_14 + DENDRO_C1_k1_1) +
                  DENDRO_15 * DENDRO_C2_k2_1 * (DENDRO_16 + DENDRO_C1_k2_1) +
                  DENDRO_17 * DENDRO_C2_k1_2 * (DENDRO_16 + DENDRO_C1_k1_2) +
                  DENDRO_17 * DENDRO_C2_k2_2 * (DENDRO_19 + DENDRO_C1_k2_2) +
                  DENDRO_22 * (DENDRO_20 + 2 * DENDRO_21) +
                  DENDRO_22 * (2 * DENDRO_20 + DENDRO_21) +
                  DENDRO_22 * (DENDRO_13 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k1_0) +
                  DENDRO_22 * (DENDRO_14 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k1_1) +
                  DENDRO_22 * (DENDRO_16 * DENDRO_C2_k2_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k2_1) +
                  DENDRO_22 * (DENDRO_18 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k2_0) +
                  DENDRO_25 * (DENDRO_23 + 2 * DENDRO_24) +
                  DENDRO_25 * (2 * DENDRO_23 + DENDRO_24) +
                  DENDRO_25 * (DENDRO_13 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_0) +
                  DENDRO_25 * (DENDRO_16 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k1_2) +
                  DENDRO_25 * (DENDRO_18 * DENDRO_C2_k2_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k2_0) +
                  DENDRO_25 * (DENDRO_19 * DENDRO_C2_k2_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k2_2) +
                  DENDRO_28 * (DENDRO_26 + 2 * DENDRO_27) +
                  DENDRO_28 * (2 * DENDRO_26 + DENDRO_27) +
                  DENDRO_28 * (DENDRO_14 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_1) +
                  DENDRO_28 * (DENDRO_16 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k1_2) +
                  DENDRO_28 * (DENDRO_16 * DENDRO_C2_k2_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k2_1) +
                  DENDRO_28 * (DENDRO_19 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k2_2) +
                  12 * DENDRO_7 * DENDRO_C1_k0_0 +
                  4.0 * DENDRO_8 * DENDRO_C1_k0_0 +
                  4 * DENDRO_9 * (DENDRO_13 + DENDRO_C1_k1_0) +
                  12 * DENDRO_C1_k0_1 * DENDRO_C2_k0_1 * DENDRO_igt3 +
                  12 * DENDRO_C1_k0_2 * DENDRO_C2_k0_2 * DENDRO_igt5 -
                  2.0 * DENDRO_igt0 * grad2_0_0_gt0 -
                  4.0 * DENDRO_igt1 * grad2_0_1_gt0 -
                  4.0 * DENDRO_igt2 * grad2_0_2_gt0 -
                  2.0 * DENDRO_igt3 * grad2_1_1_gt0 -
                  4.0 * DENDRO_igt4 * grad2_1_2_gt0 -
                  2.0 * DENDRO_igt5 * grad2_2_2_gt0 +
                  4.0 * grad_0_Gt0 * gt0[pp] + 4.0 * grad_0_Gt1 * gt1[pp] +
                  4.0 * grad_0_Gt2 * gt2[pp]) -
             DENDRO_1 -
             DENDRO_2 *
                 (DENDRO_C2_k0_0 * grad_0_chi + DENDRO_C2_k1_0 * grad_1_chi +
                  DENDRO_C2_k2_0 * grad_2_chi - grad2_0_0_chi) +
             gt0[pp] *
                 (-DENDRO_2 * (DENDRO_10 * grad_1_chi + DENDRO_12 * grad_2_chi +
                               DENDRO_8 * grad_0_chi) +
                  DENDRO_4 *
                      (DENDRO_2 * grad2_0_1_chi - DENDRO_3 * grad_1_chi) +
                  DENDRO_5 *
                      (DENDRO_2 * grad2_0_2_chi - DENDRO_3 * grad_2_chi) +
                  DENDRO_6 *
                      (DENDRO_2 * grad2_1_2_chi - 3 * grad_1_chi * grad_2_chi) +
                  DENDRO_igt0 * (-3 * DENDRO_1 + DENDRO_2 * grad2_0_0_chi) +
                  DENDRO_igt3 * (DENDRO_2 * grad2_1_1_chi -
                                 3 * (grad_1_chi * grad_1_chi)) +
                  DENDRO_igt5 * (DENDRO_2 * grad2_2_2_chi -
                                 3 * (grad_2_chi * grad_2_chi)))) /
            DENDRO_0;
        // Dendro: reduced ops: 263
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 413
        // Dendro: printing temp variables
        const double DENDRO_0 = pow(chi[pp], 2);
        const double DENDRO_1 = grad_0_chi * grad_1_chi;
        const double DENDRO_2 = 2 * chi[pp];
        const double DENDRO_3 = 2 * DENDRO_igt1;
        const double DENDRO_4 = 3 * grad_2_chi;
        const double DENDRO_5 = 2 * DENDRO_igt2;
        const double DENDRO_6 = 2 * DENDRO_igt4;
        const double DENDRO_7 =
            DENDRO_3 * DENDRO_C2_k0_1 + DENDRO_5 * DENDRO_C2_k0_2 +
            DENDRO_6 * DENDRO_C2_k0_4 + DENDRO_C2_k0_0 * DENDRO_igt0 +
            DENDRO_C2_k0_3 * DENDRO_igt3 + DENDRO_C2_k0_5 * DENDRO_igt5;
        const double DENDRO_8 =
            DENDRO_3 * DENDRO_C2_k1_1 + DENDRO_5 * DENDRO_C2_k1_2 +
            DENDRO_6 * DENDRO_C2_k1_4 + DENDRO_C2_k1_0 * DENDRO_igt0 +
            DENDRO_C2_k1_3 * DENDRO_igt3 + DENDRO_C2_k1_5 * DENDRO_igt5;
        const double DENDRO_9 =
            DENDRO_3 * DENDRO_C2_k2_1 + DENDRO_5 * DENDRO_C2_k2_2 +
            DENDRO_6 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
            DENDRO_C2_k2_3 * DENDRO_igt3 + DENDRO_C2_k2_5 * DENDRO_igt5;
        const double DENDRO_10 = 2.0 * gt1[pp];
        const double DENDRO_11 = 4 * DENDRO_igt0;
        const double DENDRO_12 = 4 * DENDRO_igt1;
        const double DENDRO_13 = 4 * DENDRO_igt3;
        const double DENDRO_14 = 4 * DENDRO_igt5;
        const double DENDRO_15 =
            DENDRO_C1_k1_1 * DENDRO_C2_k1_1 + DENDRO_C1_k1_3 * DENDRO_C2_k1_0;
        const double DENDRO_16 = 4 * DENDRO_igt2;
        const double DENDRO_17 =
            DENDRO_C1_k1_1 * DENDRO_C2_k1_2 + DENDRO_C1_k1_4 * DENDRO_C2_k1_0;
        const double DENDRO_18 = 4 * DENDRO_igt4;
        const double DENDRO_19 =
            DENDRO_C1_k1_3 * DENDRO_C2_k1_2 + DENDRO_C1_k1_4 * DENDRO_C2_k1_1;

        // Dendro: printing variables
        //--
        DENDRO_RIJ1 =
            (1.0 / 4.0) *
            (DENDRO_0 *
                 (DENDRO_10 * grad_0_Gt0 + DENDRO_10 * grad_1_Gt1 +
                  DENDRO_11 * (DENDRO_C1_k0_1 * DENDRO_C2_k1_1 +
                               2 * DENDRO_C1_k1_1 * DENDRO_C2_k1_0) +
                  DENDRO_11 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k0_1 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k0_0) +
                  DENDRO_11 * (DENDRO_C1_k0_2 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k2_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k2_0) +
                  DENDRO_12 * (DENDRO_15 + DENDRO_C1_k0_1 * DENDRO_C2_k1_3) +
                  DENDRO_12 * (DENDRO_15 + DENDRO_C1_k0_3 * DENDRO_C2_k1_1) +
                  DENDRO_12 * (2 * DENDRO_C1_k0_1 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k0_0) +
                  DENDRO_12 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k0_3 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k0_1) +
                  DENDRO_12 * (DENDRO_C1_k0_2 * DENDRO_C2_k2_3 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k2_0) +
                  DENDRO_12 * (DENDRO_C1_k0_4 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k2_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k2_1) +
                  DENDRO_13 * (DENDRO_C1_k0_3 * DENDRO_C2_k1_3 +
                               2 * DENDRO_C1_k1_3 * DENDRO_C2_k1_1) +
                  DENDRO_13 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k0_3 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k0_1) +
                  DENDRO_13 * (DENDRO_C1_k0_4 * DENDRO_C2_k2_3 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k2_1) +
                  DENDRO_14 * (DENDRO_C1_k0_4 * DENDRO_C2_k1_4 +
                               2 * DENDRO_C1_k1_4 * DENDRO_C2_k1_2) +
                  DENDRO_14 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k0_2) +
                  DENDRO_14 * (DENDRO_C1_k0_5 * DENDRO_C2_k2_4 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k2_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k2_2) +
                  DENDRO_16 * (DENDRO_17 + DENDRO_C1_k0_1 * DENDRO_C2_k1_4) +
                  DENDRO_16 * (DENDRO_17 + DENDRO_C1_k0_4 * DENDRO_C2_k1_1) +
                  DENDRO_16 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k0_2) +
                  DENDRO_16 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k0_2 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k0_0) +
                  DENDRO_16 * (DENDRO_C1_k0_2 * DENDRO_C2_k2_4 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k2_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k2_0) +
                  DENDRO_16 * (DENDRO_C1_k0_5 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k2_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k2_2) +
                  DENDRO_18 * (DENDRO_19 + DENDRO_C1_k0_3 * DENDRO_C2_k1_4) +
                  DENDRO_18 * (DENDRO_19 + DENDRO_C1_k0_4 * DENDRO_C2_k1_3) +
                  DENDRO_18 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k0_2) +
                  DENDRO_18 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k0_3 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k0_1) +
                  DENDRO_18 * (DENDRO_C1_k0_4 * DENDRO_C2_k2_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k2_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k2_1) +
                  DENDRO_18 * (DENDRO_C1_k0_5 * DENDRO_C2_k2_3 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k2_2) +
                  2.0 * DENDRO_7 * (DENDRO_C1_k0_1 + DENDRO_C1_k1_0) +
                  2.0 * DENDRO_8 * (DENDRO_C1_k0_3 + DENDRO_C1_k1_1) +
                  2.0 * DENDRO_9 * (DENDRO_C1_k0_4 + DENDRO_C1_k1_2) -
                  2.0 * DENDRO_igt0 * grad2_0_0_gt1 -
                  4.0 * DENDRO_igt1 * grad2_0_1_gt1 -
                  4.0 * DENDRO_igt2 * grad2_0_2_gt1 -
                  2.0 * DENDRO_igt3 * grad2_1_1_gt1 -
                  4.0 * DENDRO_igt4 * grad2_1_2_gt1 -
                  2.0 * DENDRO_igt5 * grad2_2_2_gt1 +
                  2.0 * grad_0_Gt1 * gt3[pp] + 2.0 * grad_0_Gt2 * gt4[pp] +
                  2.0 * grad_1_Gt0 * gt0[pp] + 2.0 * grad_1_Gt2 * gt2[pp]) -
             DENDRO_1 -
             DENDRO_2 *
                 (DENDRO_C2_k0_1 * grad_0_chi + DENDRO_C2_k1_1 * grad_1_chi +
                  DENDRO_C2_k2_1 * grad_2_chi - grad2_0_1_chi) +
             gt1[pp] *
                 (-DENDRO_2 * (DENDRO_7 * grad_0_chi + DENDRO_8 * grad_1_chi +
                               DENDRO_9 * grad_2_chi) +
                  DENDRO_3 * (-3 * DENDRO_1 + DENDRO_2 * grad2_0_1_chi) +
                  DENDRO_5 *
                      (DENDRO_2 * grad2_0_2_chi - DENDRO_4 * grad_0_chi) +
                  DENDRO_6 *
                      (DENDRO_2 * grad2_1_2_chi - DENDRO_4 * grad_1_chi) +
                  DENDRO_igt0 * (DENDRO_2 * grad2_0_0_chi -
                                 3 * (grad_0_chi * grad_0_chi)) +
                  DENDRO_igt3 * (DENDRO_2 * grad2_1_1_chi -
                                 3 * (grad_1_chi * grad_1_chi)) +
                  DENDRO_igt5 * (DENDRO_2 * grad2_2_2_chi -
                                 3 * (grad_2_chi * grad_2_chi)))) /
            DENDRO_0;
        // Dendro: reduced ops: 321
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 413
        // Dendro: printing temp variables
        const double DENDRO_0 = pow(chi[pp], 2);
        const double DENDRO_1 = grad_0_chi * grad_2_chi;
        const double DENDRO_2 = 2 * chi[pp];
        const double DENDRO_3 = 3 * grad_1_chi;
        const double DENDRO_4 = 2 * DENDRO_igt1;
        const double DENDRO_5 = 2 * DENDRO_igt2;
        const double DENDRO_6 = 2 * DENDRO_igt4;
        const double DENDRO_7 =
            DENDRO_4 * DENDRO_C2_k0_1 + DENDRO_5 * DENDRO_C2_k0_2 +
            DENDRO_6 * DENDRO_C2_k0_4 + DENDRO_C2_k0_0 * DENDRO_igt0 +
            DENDRO_C2_k0_3 * DENDRO_igt3 + DENDRO_C2_k0_5 * DENDRO_igt5;
        const double DENDRO_8 =
            DENDRO_4 * DENDRO_C2_k1_1 + DENDRO_5 * DENDRO_C2_k1_2 +
            DENDRO_6 * DENDRO_C2_k1_4 + DENDRO_C2_k1_0 * DENDRO_igt0 +
            DENDRO_C2_k1_3 * DENDRO_igt3 + DENDRO_C2_k1_5 * DENDRO_igt5;
        const double DENDRO_9 =
            DENDRO_4 * DENDRO_C2_k2_1 + DENDRO_5 * DENDRO_C2_k2_2 +
            DENDRO_6 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
            DENDRO_C2_k2_3 * DENDRO_igt3 + DENDRO_C2_k2_5 * DENDRO_igt5;
        const double DENDRO_10 = 2.0 * gt2[pp];
        const double DENDRO_11 = 4 * DENDRO_igt0;
        const double DENDRO_12 = 4 * DENDRO_igt2;
        const double DENDRO_13 = 4 * DENDRO_igt3;
        const double DENDRO_14 = 4 * DENDRO_igt5;
        const double DENDRO_15 = 4 * DENDRO_igt1;
        const double DENDRO_16 =
            DENDRO_C1_k2_2 * DENDRO_C2_k2_1 + DENDRO_C1_k2_4 * DENDRO_C2_k2_0;
        const double DENDRO_17 =
            DENDRO_C1_k2_2 * DENDRO_C2_k2_2 + DENDRO_C1_k2_5 * DENDRO_C2_k2_0;
        const double DENDRO_18 = 4 * DENDRO_igt4;
        const double DENDRO_19 =
            DENDRO_C1_k2_4 * DENDRO_C2_k2_2 + DENDRO_C1_k2_5 * DENDRO_C2_k2_1;

        // Dendro: printing variables
        //--
        DENDRO_RIJ2 =
            (1.0 / 4.0) *
            (DENDRO_0 *
                 (DENDRO_10 * grad_0_Gt0 + DENDRO_10 * grad_2_Gt2 +
                  DENDRO_11 * (DENDRO_C1_k0_2 * DENDRO_C2_k2_2 +
                               2 * DENDRO_C1_k2_2 * DENDRO_C2_k2_0) +
                  DENDRO_11 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k0_2 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_0) +
                  DENDRO_11 * (DENDRO_C1_k0_1 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_0) +
                  DENDRO_12 * (DENDRO_17 + DENDRO_C1_k0_2 * DENDRO_C2_k2_5) +
                  DENDRO_12 * (DENDRO_17 + DENDRO_C1_k0_5 * DENDRO_C2_k2_2) +
                  DENDRO_12 * (2 * DENDRO_C1_k0_2 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_0) +
                  DENDRO_12 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_2) +
                  DENDRO_12 * (DENDRO_C1_k0_1 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_2) +
                  DENDRO_12 * (DENDRO_C1_k0_4 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_0) +
                  DENDRO_13 * (DENDRO_C1_k0_4 * DENDRO_C2_k2_4 +
                               2 * DENDRO_C1_k2_4 * DENDRO_C2_k2_1) +
                  DENDRO_13 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_1) +
                  DENDRO_13 * (DENDRO_C1_k0_3 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_1) +
                  DENDRO_14 * (DENDRO_C1_k0_5 * DENDRO_C2_k2_5 +
                               2 * DENDRO_C1_k2_5 * DENDRO_C2_k2_2) +
                  DENDRO_14 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_2) +
                  DENDRO_14 * (DENDRO_C1_k0_4 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_2) +
                  DENDRO_15 * (DENDRO_16 + DENDRO_C1_k0_2 * DENDRO_C2_k2_4) +
                  DENDRO_15 * (DENDRO_16 + DENDRO_C1_k0_4 * DENDRO_C2_k2_2) +
                  DENDRO_15 * (DENDRO_C1_k0_0 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_0 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_1) +
                  DENDRO_15 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k0_2 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_0) +
                  DENDRO_15 * (DENDRO_C1_k0_1 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_0 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_1) +
                  DENDRO_15 * (DENDRO_C1_k0_3 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_0) +
                  DENDRO_18 * (DENDRO_19 + DENDRO_C1_k0_4 * DENDRO_C2_k2_5) +
                  DENDRO_18 * (DENDRO_19 + DENDRO_C1_k0_5 * DENDRO_C2_k2_4) +
                  DENDRO_18 * (DENDRO_C1_k0_1 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_2) +
                  DENDRO_18 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_1) +
                  DENDRO_18 * (DENDRO_C1_k0_3 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_2) +
                  DENDRO_18 * (DENDRO_C1_k0_4 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_1) +
                  2.0 * DENDRO_7 * (DENDRO_C1_k0_2 + DENDRO_C1_k2_0) +
                  2.0 * DENDRO_8 * (DENDRO_C1_k0_4 + DENDRO_C1_k2_1) +
                  2.0 * DENDRO_9 * (DENDRO_C1_k0_5 + DENDRO_C1_k2_2) -
                  2.0 * DENDRO_igt0 * grad2_0_0_gt2 -
                  4.0 * DENDRO_igt1 * grad2_0_1_gt2 -
                  4.0 * DENDRO_igt2 * grad2_0_2_gt2 -
                  2.0 * DENDRO_igt3 * grad2_1_1_gt2 -
                  4.0 * DENDRO_igt4 * grad2_1_2_gt2 -
                  2.0 * DENDRO_igt5 * grad2_2_2_gt2 +
                  2.0 * grad_0_Gt1 * gt4[pp] + 2.0 * grad_0_Gt2 * gt5[pp] +
                  2.0 * grad_2_Gt0 * gt0[pp] + 2.0 * grad_2_Gt1 * gt1[pp]) -
             DENDRO_1 -
             DENDRO_2 *
                 (DENDRO_C2_k0_2 * grad_0_chi + DENDRO_C2_k1_2 * grad_1_chi +
                  DENDRO_C2_k2_2 * grad_2_chi - grad2_0_2_chi) +
             gt2[pp] *
                 (-DENDRO_2 * (DENDRO_7 * grad_0_chi + DENDRO_8 * grad_1_chi +
                               DENDRO_9 * grad_2_chi) +
                  DENDRO_4 *
                      (DENDRO_2 * grad2_0_1_chi - DENDRO_3 * grad_0_chi) +
                  DENDRO_5 * (-3 * DENDRO_1 + DENDRO_2 * grad2_0_2_chi) +
                  DENDRO_6 *
                      (DENDRO_2 * grad2_1_2_chi - DENDRO_3 * grad_2_chi) +
                  DENDRO_igt0 * (DENDRO_2 * grad2_0_0_chi -
                                 3 * (grad_0_chi * grad_0_chi)) +
                  DENDRO_igt3 * (DENDRO_2 * grad2_1_1_chi -
                                 3 * (grad_1_chi * grad_1_chi)) +
                  DENDRO_igt5 * (DENDRO_2 * grad2_2_2_chi -
                                 3 * (grad_2_chi * grad_2_chi)))) /
            DENDRO_0;
        // Dendro: reduced ops: 321
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 363
        // Dendro: printing temp variables
        const double DENDRO_0 = pow(chi[pp], 2);
        const double DENDRO_1 = (grad_1_chi * grad_1_chi);
        const double DENDRO_2 = 2 * chi[pp];
        const double DENDRO_3 = 3 * grad_0_chi;
        const double DENDRO_4 = 2 * DENDRO_igt1;
        const double DENDRO_5 = 2 * DENDRO_igt2;
        const double DENDRO_6 = 2 * DENDRO_igt4;
        const double DENDRO_7 = DENDRO_C2_k0_3 * DENDRO_igt3;
        const double DENDRO_8 =
            DENDRO_4 * DENDRO_C2_k0_1 + DENDRO_5 * DENDRO_C2_k0_2 +
            DENDRO_6 * DENDRO_C2_k0_4 + DENDRO_7 +
            DENDRO_C2_k0_0 * DENDRO_igt0 + DENDRO_C2_k0_5 * DENDRO_igt5;
        const double DENDRO_9 = DENDRO_C2_k1_3 * DENDRO_igt3;
        const double DENDRO_10 =
            DENDRO_4 * DENDRO_C2_k1_1 + DENDRO_5 * DENDRO_C2_k1_2 +
            DENDRO_6 * DENDRO_C2_k1_4 + DENDRO_9 +
            DENDRO_C2_k1_0 * DENDRO_igt0 + DENDRO_C2_k1_5 * DENDRO_igt5;
        const double DENDRO_11 = DENDRO_C2_k2_3 * DENDRO_igt3;
        const double DENDRO_12 =
            DENDRO_11 + DENDRO_4 * DENDRO_C2_k2_1 + DENDRO_5 * DENDRO_C2_k2_2 +
            DENDRO_6 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
            DENDRO_C2_k2_5 * DENDRO_igt5;
        const double DENDRO_13 = 2 * DENDRO_C1_k1_0;
        const double DENDRO_14 = 4 * DENDRO_igt0;
        const double DENDRO_15 = 2 * DENDRO_C1_k1_1;
        const double DENDRO_16 = 2 * DENDRO_C1_k1_2;
        const double DENDRO_17 = 4 * DENDRO_igt5;
        const double DENDRO_18 = 2 * DENDRO_C1_k1_4;
        const double DENDRO_19 = 2 * DENDRO_C1_k1_5;
        const double DENDRO_20 = 4 * DENDRO_igt1;
        const double DENDRO_21 = DENDRO_C1_k1_1 * DENDRO_C2_k1_3;
        const double DENDRO_22 = DENDRO_C1_k1_3 * DENDRO_C2_k1_1;
        const double DENDRO_23 = 4 * DENDRO_igt2;
        const double DENDRO_24 = DENDRO_C1_k1_1 * DENDRO_C2_k1_4;
        const double DENDRO_25 = DENDRO_C1_k1_4 * DENDRO_C2_k1_1;
        const double DENDRO_26 = 4 * DENDRO_igt4;
        const double DENDRO_27 = DENDRO_C1_k1_3 * DENDRO_C2_k1_4;
        const double DENDRO_28 = DENDRO_C1_k1_4 * DENDRO_C2_k1_3;

        // Dendro: printing variables
        //--
        DENDRO_RIJ3 =
            (1.0 / 4.0) *
            (DENDRO_0 *
                 (4.0 * DENDRO_10 * DENDRO_C1_k1_3 +
                  4 * DENDRO_11 * (DENDRO_18 + DENDRO_C1_k2_3) +
                  4.0 * DENDRO_12 * DENDRO_C1_k1_4 +
                  DENDRO_14 * DENDRO_C2_k0_1 * (DENDRO_13 + DENDRO_C1_k0_1) +
                  DENDRO_14 * DENDRO_C2_k2_1 * (DENDRO_16 + DENDRO_C1_k2_1) +
                  DENDRO_17 * DENDRO_C2_k0_4 * (DENDRO_16 + DENDRO_C1_k0_4) +
                  DENDRO_17 * DENDRO_C2_k2_4 * (DENDRO_19 + DENDRO_C1_k2_4) +
                  DENDRO_20 * (DENDRO_21 + 2 * DENDRO_22) +
                  DENDRO_20 * (2 * DENDRO_21 + DENDRO_22) +
                  DENDRO_20 * (DENDRO_13 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k0_3 * DENDRO_C2_k0_1) +
                  DENDRO_20 * (DENDRO_15 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k0_1 * DENDRO_C2_k0_3) +
                  DENDRO_20 * (DENDRO_16 * DENDRO_C2_k2_3 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k2_1) +
                  DENDRO_20 * (DENDRO_18 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k2_3) +
                  DENDRO_23 * (DENDRO_24 + 2 * DENDRO_25) +
                  DENDRO_23 * (2 * DENDRO_24 + DENDRO_25) +
                  DENDRO_23 * (DENDRO_13 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_1) +
                  DENDRO_23 * (DENDRO_16 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k0_1 * DENDRO_C2_k0_4) +
                  DENDRO_23 * (DENDRO_16 * DENDRO_C2_k2_4 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k2_1) +
                  DENDRO_23 * (DENDRO_19 * DENDRO_C2_k2_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k2_4) +
                  DENDRO_26 * (DENDRO_27 + 2 * DENDRO_28) +
                  DENDRO_26 * (2 * DENDRO_27 + DENDRO_28) +
                  DENDRO_26 * (DENDRO_15 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_3) +
                  DENDRO_26 * (DENDRO_16 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k0_3 * DENDRO_C2_k0_4) +
                  DENDRO_26 * (DENDRO_18 * DENDRO_C2_k2_4 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k2_3) +
                  DENDRO_26 * (DENDRO_19 * DENDRO_C2_k2_3 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k2_4) +
                  4 * DENDRO_7 * (DENDRO_15 + DENDRO_C1_k0_3) +
                  4.0 * DENDRO_8 * DENDRO_C1_k1_1 +
                  12 * DENDRO_9 * DENDRO_C1_k1_3 +
                  12 * DENDRO_C1_k1_1 * DENDRO_C2_k1_1 * DENDRO_igt0 +
                  12 * DENDRO_C1_k1_4 * DENDRO_C2_k1_4 * DENDRO_igt5 -
                  2.0 * DENDRO_igt0 * grad2_0_0_gt3 -
                  4.0 * DENDRO_igt1 * grad2_0_1_gt3 -
                  4.0 * DENDRO_igt2 * grad2_0_2_gt3 -
                  2.0 * DENDRO_igt3 * grad2_1_1_gt3 -
                  4.0 * DENDRO_igt4 * grad2_1_2_gt3 -
                  2.0 * DENDRO_igt5 * grad2_2_2_gt3 +
                  4.0 * grad_1_Gt0 * gt1[pp] + 4.0 * grad_1_Gt1 * gt3[pp] +
                  4.0 * grad_1_Gt2 * gt4[pp]) -
             DENDRO_1 -
             DENDRO_2 *
                 (DENDRO_C2_k0_3 * grad_0_chi + DENDRO_C2_k1_3 * grad_1_chi +
                  DENDRO_C2_k2_3 * grad_2_chi - grad2_1_1_chi) +
             gt3[pp] *
                 (-DENDRO_2 * (DENDRO_10 * grad_1_chi + DENDRO_12 * grad_2_chi +
                               DENDRO_8 * grad_0_chi) +
                  DENDRO_4 *
                      (DENDRO_2 * grad2_0_1_chi - DENDRO_3 * grad_1_chi) +
                  DENDRO_5 *
                      (DENDRO_2 * grad2_0_2_chi - DENDRO_3 * grad_2_chi) +
                  DENDRO_6 *
                      (DENDRO_2 * grad2_1_2_chi - 3 * grad_1_chi * grad_2_chi) +
                  DENDRO_igt0 * (DENDRO_2 * grad2_0_0_chi -
                                 3 * (grad_0_chi * grad_0_chi)) +
                  DENDRO_igt3 * (-3 * DENDRO_1 + DENDRO_2 * grad2_1_1_chi) +
                  DENDRO_igt5 * (DENDRO_2 * grad2_2_2_chi -
                                 3 * (grad_2_chi * grad_2_chi)))) /
            DENDRO_0;
        // Dendro: reduced ops: 263
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 413
        // Dendro: printing temp variables
        const double DENDRO_0 = pow(chi[pp], 2);
        const double DENDRO_1 = grad_1_chi * grad_2_chi;
        const double DENDRO_2 = 2 * chi[pp];
        const double DENDRO_3 = 3 * grad_0_chi;
        const double DENDRO_4 = 2 * DENDRO_igt1;
        const double DENDRO_5 = 2 * DENDRO_igt2;
        const double DENDRO_6 = 2 * DENDRO_igt4;
        const double DENDRO_7 =
            DENDRO_4 * DENDRO_C2_k0_1 + DENDRO_5 * DENDRO_C2_k0_2 +
            DENDRO_6 * DENDRO_C2_k0_4 + DENDRO_C2_k0_0 * DENDRO_igt0 +
            DENDRO_C2_k0_3 * DENDRO_igt3 + DENDRO_C2_k0_5 * DENDRO_igt5;
        const double DENDRO_8 =
            DENDRO_4 * DENDRO_C2_k1_1 + DENDRO_5 * DENDRO_C2_k1_2 +
            DENDRO_6 * DENDRO_C2_k1_4 + DENDRO_C2_k1_0 * DENDRO_igt0 +
            DENDRO_C2_k1_3 * DENDRO_igt3 + DENDRO_C2_k1_5 * DENDRO_igt5;
        const double DENDRO_9 =
            DENDRO_4 * DENDRO_C2_k2_1 + DENDRO_5 * DENDRO_C2_k2_2 +
            DENDRO_6 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
            DENDRO_C2_k2_3 * DENDRO_igt3 + DENDRO_C2_k2_5 * DENDRO_igt5;
        const double DENDRO_10 = 2.0 * gt4[pp];
        const double DENDRO_11 = 4 * DENDRO_igt0;
        const double DENDRO_12 = 4 * DENDRO_igt3;
        const double DENDRO_13 = 4 * DENDRO_igt4;
        const double DENDRO_14 = 4 * DENDRO_igt5;
        const double DENDRO_15 = 4 * DENDRO_igt1;
        const double DENDRO_16 =
            DENDRO_C1_k2_2 * DENDRO_C2_k2_3 + DENDRO_C1_k2_4 * DENDRO_C2_k2_1;
        const double DENDRO_17 = 4 * DENDRO_igt2;
        const double DENDRO_18 =
            DENDRO_C1_k2_2 * DENDRO_C2_k2_4 + DENDRO_C1_k2_5 * DENDRO_C2_k2_1;
        const double DENDRO_19 =
            DENDRO_C1_k2_4 * DENDRO_C2_k2_4 + DENDRO_C1_k2_5 * DENDRO_C2_k2_3;

        // Dendro: printing variables
        //--
        DENDRO_RIJ4 =
            (1.0 / 4.0) *
            (DENDRO_0 *
                 (DENDRO_10 * grad_1_Gt1 + DENDRO_10 * grad_2_Gt2 +
                  DENDRO_11 * (DENDRO_C1_k1_2 * DENDRO_C2_k2_2 +
                               2 * DENDRO_C1_k2_2 * DENDRO_C2_k2_1) +
                  DENDRO_11 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_1) +
                  DENDRO_11 * (DENDRO_C1_k1_1 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_1) +
                  DENDRO_12 * (DENDRO_C1_k1_4 * DENDRO_C2_k2_4 +
                               2 * DENDRO_C1_k2_4 * DENDRO_C2_k2_3) +
                  DENDRO_12 * (DENDRO_C1_k0_4 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_3) +
                  DENDRO_12 * (DENDRO_C1_k1_3 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_3 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_3) +
                  DENDRO_13 * (DENDRO_19 + DENDRO_C1_k1_4 * DENDRO_C2_k2_5) +
                  DENDRO_13 * (DENDRO_19 + DENDRO_C1_k1_5 * DENDRO_C2_k2_4) +
                  DENDRO_13 * (2 * DENDRO_C1_k1_4 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_3) +
                  DENDRO_13 * (DENDRO_C1_k0_4 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_3) +
                  DENDRO_13 * (DENDRO_C1_k0_5 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_4) +
                  DENDRO_13 * (DENDRO_C1_k1_3 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_3 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_4) +
                  DENDRO_14 * (DENDRO_C1_k1_5 * DENDRO_C2_k2_5 +
                               2 * DENDRO_C1_k2_5 * DENDRO_C2_k2_4) +
                  DENDRO_14 * (DENDRO_C1_k0_5 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_4) +
                  DENDRO_14 * (DENDRO_C1_k1_4 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_4) +
                  DENDRO_15 * (DENDRO_16 + DENDRO_C1_k1_2 * DENDRO_C2_k2_4) +
                  DENDRO_15 * (DENDRO_16 + DENDRO_C1_k1_4 * DENDRO_C2_k2_2) +
                  DENDRO_15 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_3 +
                               DENDRO_C1_k1_1 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k0_1) +
                  DENDRO_15 * (DENDRO_C1_k0_4 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_3) +
                  DENDRO_15 * (DENDRO_C1_k1_1 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_3) +
                  DENDRO_15 * (DENDRO_C1_k1_2 * DENDRO_C2_k1_3 +
                               DENDRO_C1_k1_3 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_3 * DENDRO_C2_k1_1) +
                  DENDRO_17 * (DENDRO_18 + DENDRO_C1_k1_2 * DENDRO_C2_k2_5) +
                  DENDRO_17 * (DENDRO_18 + DENDRO_C1_k1_5 * DENDRO_C2_k2_2) +
                  DENDRO_17 * (DENDRO_C1_k0_2 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k2_2 * DENDRO_C2_k0_1) +
                  DENDRO_17 * (DENDRO_C1_k0_5 * DENDRO_C2_k0_1 +
                               DENDRO_C1_k1_0 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k2_0 * DENDRO_C2_k0_4) +
                  DENDRO_17 * (DENDRO_C1_k1_1 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_1 +
                               DENDRO_C1_k2_1 * DENDRO_C2_k1_4) +
                  DENDRO_17 * (DENDRO_C1_k1_2 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k2_4 * DENDRO_C2_k1_1) +
                  2.0 * DENDRO_7 * (DENDRO_C1_k1_2 + DENDRO_C1_k2_1) +
                  2.0 * DENDRO_8 * (DENDRO_C1_k1_4 + DENDRO_C1_k2_3) +
                  2.0 * DENDRO_9 * (DENDRO_C1_k1_5 + DENDRO_C1_k2_4) -
                  2.0 * DENDRO_igt0 * grad2_0_0_gt4 -
                  4.0 * DENDRO_igt1 * grad2_0_1_gt4 -
                  4.0 * DENDRO_igt2 * grad2_0_2_gt4 -
                  2.0 * DENDRO_igt3 * grad2_1_1_gt4 -
                  4.0 * DENDRO_igt4 * grad2_1_2_gt4 -
                  2.0 * DENDRO_igt5 * grad2_2_2_gt4 +
                  2.0 * grad_1_Gt0 * gt2[pp] + 2.0 * grad_1_Gt2 * gt5[pp] +
                  2.0 * grad_2_Gt0 * gt1[pp] + 2.0 * grad_2_Gt1 * gt3[pp]) -
             DENDRO_1 -
             DENDRO_2 *
                 (DENDRO_C2_k0_4 * grad_0_chi + DENDRO_C2_k1_4 * grad_1_chi +
                  DENDRO_C2_k2_4 * grad_2_chi - grad2_1_2_chi) +
             gt4[pp] *
                 (-DENDRO_2 * (DENDRO_7 * grad_0_chi + DENDRO_8 * grad_1_chi +
                               DENDRO_9 * grad_2_chi) +
                  DENDRO_4 *
                      (DENDRO_2 * grad2_0_1_chi - DENDRO_3 * grad_1_chi) +
                  DENDRO_5 *
                      (DENDRO_2 * grad2_0_2_chi - DENDRO_3 * grad_2_chi) +
                  DENDRO_6 * (-3 * DENDRO_1 + DENDRO_2 * grad2_1_2_chi) +
                  DENDRO_igt0 * (DENDRO_2 * grad2_0_0_chi -
                                 3 * (grad_0_chi * grad_0_chi)) +
                  DENDRO_igt3 * (DENDRO_2 * grad2_1_1_chi -
                                 3 * (grad_1_chi * grad_1_chi)) +
                  DENDRO_igt5 * (DENDRO_2 * grad2_2_2_chi -
                                 3 * (grad_2_chi * grad_2_chi)))) /
            DENDRO_0;
        // Dendro: reduced ops: 321
        // Dendro: }}}
    }
    {
        // Dendro: {{{
        // Dendro: original ops: 363
        // Dendro: printing temp variables
        const double DENDRO_0 = pow(chi[pp], 2);
        const double DENDRO_1 = (grad_2_chi * grad_2_chi);
        const double DENDRO_2 = 2 * chi[pp];
        const double DENDRO_3 = 3 * grad_0_chi;
        const double DENDRO_4 = 2 * DENDRO_igt1;
        const double DENDRO_5 = 2 * DENDRO_igt2;
        const double DENDRO_6 = 2 * DENDRO_igt4;
        const double DENDRO_7 = DENDRO_C2_k0_5 * DENDRO_igt5;
        const double DENDRO_8 =
            DENDRO_4 * DENDRO_C2_k0_1 + DENDRO_5 * DENDRO_C2_k0_2 +
            DENDRO_6 * DENDRO_C2_k0_4 + DENDRO_7 +
            DENDRO_C2_k0_0 * DENDRO_igt0 + DENDRO_C2_k0_3 * DENDRO_igt3;
        const double DENDRO_9 = DENDRO_C2_k1_5 * DENDRO_igt5;
        const double DENDRO_10 =
            DENDRO_4 * DENDRO_C2_k1_1 + DENDRO_5 * DENDRO_C2_k1_2 +
            DENDRO_6 * DENDRO_C2_k1_4 + DENDRO_9 +
            DENDRO_C2_k1_0 * DENDRO_igt0 + DENDRO_C2_k1_3 * DENDRO_igt3;
        const double DENDRO_11 = DENDRO_C2_k2_5 * DENDRO_igt5;
        const double DENDRO_12 =
            DENDRO_11 + DENDRO_4 * DENDRO_C2_k2_1 + DENDRO_5 * DENDRO_C2_k2_2 +
            DENDRO_6 * DENDRO_C2_k2_4 + DENDRO_C2_k2_0 * DENDRO_igt0 +
            DENDRO_C2_k2_3 * DENDRO_igt3;
        const double DENDRO_13 = 2 * DENDRO_C1_k2_0;
        const double DENDRO_14 = 4 * DENDRO_igt0;
        const double DENDRO_15 = 2 * DENDRO_C1_k2_1;
        const double DENDRO_16 = 4 * DENDRO_igt3;
        const double DENDRO_17 = 2 * DENDRO_C1_k2_2;
        const double DENDRO_18 = 2 * DENDRO_C1_k2_3;
        const double DENDRO_19 = 2 * DENDRO_C1_k2_4;
        const double DENDRO_20 = 4 * DENDRO_igt1;
        const double DENDRO_21 = DENDRO_C1_k2_2 * DENDRO_C2_k2_4;
        const double DENDRO_22 = DENDRO_C1_k2_4 * DENDRO_C2_k2_2;
        const double DENDRO_23 = 4 * DENDRO_igt2;
        const double DENDRO_24 = DENDRO_C1_k2_2 * DENDRO_C2_k2_5;
        const double DENDRO_25 = DENDRO_C1_k2_5 * DENDRO_C2_k2_2;
        const double DENDRO_26 = 4 * DENDRO_igt4;
        const double DENDRO_27 = DENDRO_C1_k2_4 * DENDRO_C2_k2_5;
        const double DENDRO_28 = DENDRO_C1_k2_5 * DENDRO_C2_k2_4;

        // Dendro: printing variables
        //--
        DENDRO_RIJ5 =
            (1.0 / 4.0) *
            (DENDRO_0 *
                 (4.0 * DENDRO_10 * DENDRO_C1_k2_4 +
                  12 * DENDRO_11 * DENDRO_C1_k2_5 +
                  4.0 * DENDRO_12 * DENDRO_C1_k2_5 +
                  DENDRO_14 * DENDRO_C2_k0_2 * (DENDRO_13 + DENDRO_C1_k0_2) +
                  DENDRO_14 * DENDRO_C2_k1_2 * (DENDRO_15 + DENDRO_C1_k1_2) +
                  DENDRO_16 * DENDRO_C2_k0_4 * (DENDRO_15 + DENDRO_C1_k0_4) +
                  DENDRO_16 * DENDRO_C2_k1_4 * (DENDRO_18 + DENDRO_C1_k1_4) +
                  DENDRO_20 * (DENDRO_21 + 2 * DENDRO_22) +
                  DENDRO_20 * (2 * DENDRO_21 + DENDRO_22) +
                  DENDRO_20 * (DENDRO_13 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_2) +
                  DENDRO_20 * (DENDRO_15 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k0_2 * DENDRO_C2_k0_4) +
                  DENDRO_20 * (DENDRO_15 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_2) +
                  DENDRO_20 * (DENDRO_18 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_4) +
                  DENDRO_23 * (DENDRO_24 + 2 * DENDRO_25) +
                  DENDRO_23 * (2 * DENDRO_24 + DENDRO_25) +
                  DENDRO_23 * (DENDRO_13 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_2) +
                  DENDRO_23 * (DENDRO_15 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_2) +
                  DENDRO_23 * (DENDRO_17 * DENDRO_C2_k0_2 +
                               DENDRO_C1_k0_2 * DENDRO_C2_k0_5) +
                  DENDRO_23 * (DENDRO_19 * DENDRO_C2_k1_2 +
                               DENDRO_C1_k1_2 * DENDRO_C2_k1_5) +
                  DENDRO_26 * (DENDRO_27 + 2 * DENDRO_28) +
                  DENDRO_26 * (2 * DENDRO_27 + DENDRO_28) +
                  DENDRO_26 * (DENDRO_15 * DENDRO_C2_k0_5 +
                               DENDRO_C1_k0_5 * DENDRO_C2_k0_4) +
                  DENDRO_26 * (DENDRO_17 * DENDRO_C2_k0_4 +
                               DENDRO_C1_k0_4 * DENDRO_C2_k0_5) +
                  DENDRO_26 * (DENDRO_18 * DENDRO_C2_k1_5 +
                               DENDRO_C1_k1_5 * DENDRO_C2_k1_4) +
                  DENDRO_26 * (DENDRO_19 * DENDRO_C2_k1_4 +
                               DENDRO_C1_k1_4 * DENDRO_C2_k1_5) +
                  4 * DENDRO_7 * (DENDRO_17 + DENDRO_C1_k0_5) +
                  4.0 * DENDRO_8 * DENDRO_C1_k2_2 +
                  4 * DENDRO_9 * (DENDRO_19 + DENDRO_C1_k1_5) +
                  12 * DENDRO_C1_k2_2 * DENDRO_C2_k2_2 * DENDRO_igt0 +
                  12 * DENDRO_C1_k2_4 * DENDRO_C2_k2_4 * DENDRO_igt3 -
                  2.0 * DENDRO_igt0 * grad2_0_0_gt5 -
                  4.0 * DENDRO_igt1 * grad2_0_1_gt5 -
                  4.0 * DENDRO_igt2 * grad2_0_2_gt5 -
                  2.0 * DENDRO_igt3 * grad2_1_1_gt5 -
                  4.0 * DENDRO_igt4 * grad2_1_2_gt5 -
                  2.0 * DENDRO_igt5 * grad2_2_2_gt5 +
                  4.0 * grad_2_Gt0 * gt2[pp] + 4.0 * grad_2_Gt1 * gt4[pp] +
                  4.0 * grad_2_Gt2 * gt5[pp]) -
             DENDRO_1 -
             DENDRO_2 *
                 (DENDRO_C2_k0_5 * grad_0_chi + DENDRO_C2_k1_5 * grad_1_chi +
                  DENDRO_C2_k2_5 * grad_2_chi - grad2_2_2_chi) +
             gt5[pp] *
                 (-DENDRO_2 * (DENDRO_10 * grad_1_chi + DENDRO_12 * grad_2_chi +
                               DENDRO_8 * grad_0_chi) +
                  DENDRO_4 *
                      (DENDRO_2 * grad2_0_1_chi - DENDRO_3 * grad_1_chi) +
                  DENDRO_5 *
                      (DENDRO_2 * grad2_0_2_chi - DENDRO_3 * grad_2_chi) +
                  DENDRO_6 *
                      (DENDRO_2 * grad2_1_2_chi - 3 * grad_1_chi * grad_2_chi) +
                  DENDRO_igt0 * (DENDRO_2 * grad2_0_0_chi -
                                 3 * (grad_0_chi * grad_0_chi)) +
                  DENDRO_igt3 * (DENDRO_2 * grad2_1_1_chi -
                                 3 * (grad_1_chi * grad_1_chi)) +
                  DENDRO_igt5 * (-3 * DENDRO_1 + DENDRO_2 * grad2_2_2_chi))) /
            DENDRO_0;
        // Dendro: reduced ops: 263
        // Dendro: }}}
    }
    // Dendro: {{{
    // Dendro: original ops: 750
    // Dendro: printing temp variables
    const double DENDRO_0  = (2.0 / 3.0) * At0[pp];
    const double DENDRO_1  = 2 * At1[pp];
    const double DENDRO_2  = 2 * At2[pp];
    const double DENDRO_3  = At1[pp] * DENDRO_igt1;
    const double DENDRO_4  = At2[pp] * DENDRO_igt2;
    const double DENDRO_5  = 2 * At0[pp];
    const double DENDRO_6  = DENDRO_C3_k0_0 * grad_0_alpha;
    const double DENDRO_7  = DENDRO_C3_k1_0 * grad_1_alpha;
    const double DENDRO_8  = DENDRO_C3_k2_0 * grad_2_alpha;
    const double DENDRO_9  = DENDRO_RIJ0 * alpha[pp];
    const double DENDRO_10 = DENDRO_C3_k0_3 * grad_0_alpha;
    const double DENDRO_11 = DENDRO_C3_k1_3 * grad_1_alpha;
    const double DENDRO_12 = DENDRO_C3_k2_3 * grad_2_alpha;
    const double DENDRO_13 = DENDRO_RIJ3 * alpha[pp];
    const double DENDRO_14 = DENDRO_C3_k0_5 * grad_0_alpha;
    const double DENDRO_15 = DENDRO_C3_k1_5 * grad_1_alpha;
    const double DENDRO_16 = DENDRO_C3_k2_5 * grad_2_alpha;
    const double DENDRO_17 = DENDRO_RIJ5 * alpha[pp];
    const double DENDRO_18 = DENDRO_C3_k0_1 * grad_0_alpha;
    const double DENDRO_19 = DENDRO_C3_k1_1 * grad_1_alpha;
    const double DENDRO_20 = DENDRO_C3_k2_1 * grad_2_alpha;
    const double DENDRO_21 = DENDRO_RIJ1 * alpha[pp];
    const double DENDRO_22 = DENDRO_C3_k0_2 * grad_0_alpha;
    const double DENDRO_23 = DENDRO_C3_k1_2 * grad_1_alpha;
    const double DENDRO_24 = DENDRO_C3_k2_2 * grad_2_alpha;
    const double DENDRO_25 = DENDRO_RIJ2 * alpha[pp];
    const double DENDRO_26 = DENDRO_C3_k0_4 * grad_0_alpha;
    const double DENDRO_27 = DENDRO_C3_k1_4 * grad_1_alpha;
    const double DENDRO_28 = DENDRO_C3_k2_4 * grad_2_alpha;
    const double DENDRO_29 = DENDRO_RIJ4 * alpha[pp];
    const double DENDRO_30 =
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
    const double DENDRO_31 = (1.0 / 3.0) * chi[pp];
    const double DENDRO_32 = (1.0 / 3.0) * At1[pp];
    const double DENDRO_33 = (2.0 / 3.0) * grad_2_beta2;
    const double DENDRO_34 =
        At1[pp] * DENDRO_igt0 + At3[pp] * DENDRO_igt1 + At4[pp] * DENDRO_igt2;
    const double DENDRO_35 = At4[pp] * DENDRO_igt4;
    const double DENDRO_36 = At3[pp] * DENDRO_igt3 + DENDRO_3 + DENDRO_35;
    const double DENDRO_37 =
        At1[pp] * DENDRO_igt2 + At3[pp] * DENDRO_igt4 + At4[pp] * DENDRO_igt5;
    const double DENDRO_38 = (1.0 / 3.0) * At2[pp];
    const double DENDRO_39 = (2.0 / 3.0) * grad_1_beta1;
    const double DENDRO_40 =
        At2[pp] * DENDRO_igt0 + At4[pp] * DENDRO_igt1 + At5[pp] * DENDRO_igt2;
    const double DENDRO_41 =
        At2[pp] * DENDRO_igt1 + At4[pp] * DENDRO_igt3 + At5[pp] * DENDRO_igt4;
    const double DENDRO_42 = At5[pp] * DENDRO_igt5 + DENDRO_35 + DENDRO_4;
    const double DENDRO_43 = (2.0 / 3.0) * grad_0_beta0;
    const double DENDRO_44 = 2 * At4[pp];
    const double DENDRO_45 = 2 * At3[pp];
    const double DENDRO_46 = (1.0 / 3.0) * At4[pp];

    // Dendro: printing variables
    //--
    At_rhs00[pp] =
        (4.0 / 3.0) * At0[pp] * grad_0_beta0 - DENDRO_0 * grad_1_beta1 -
        DENDRO_0 * grad_2_beta2 + DENDRO_1 * grad_0_beta1 +
        DENDRO_2 * grad_0_beta2 +
        DENDRO_31 * (-DENDRO_30 * gt0[pp] + 3 * DENDRO_6 + 3 * DENDRO_7 +
                     3 * DENDRO_8 + 3 * DENDRO_9 - 3 * grad2_0_0_alpha) -
        alpha[pp] * (-At0[pp] * K[pp] +
                     DENDRO_1 * (At0[pp] * DENDRO_igt1 + At1[pp] * DENDRO_igt3 +
                                 At2[pp] * DENDRO_igt4) +
                     DENDRO_2 * (At0[pp] * DENDRO_igt2 + At1[pp] * DENDRO_igt4 +
                                 At2[pp] * DENDRO_igt5) +
                     DENDRO_5 * (At0[pp] * DENDRO_igt0 + DENDRO_3 + DENDRO_4)) +
        beta0[pp] * grad_0_At0 + beta1[pp] * grad_1_At0 +
        beta2[pp] * grad_2_At0;
    //--
    At_rhs01[pp] = At0[pp] * grad_1_beta0 - At1[pp] * DENDRO_33 +
                   At2[pp] * grad_1_beta2 + At3[pp] * grad_0_beta1 +
                   At4[pp] * grad_0_beta2 +
                   DENDRO_31 * (3 * DENDRO_18 + 3 * DENDRO_19 + 3 * DENDRO_20 +
                                3 * DENDRO_21 - DENDRO_30 * gt1[pp] -
                                3 * grad2_0_1_alpha) +
                   DENDRO_32 * grad_0_beta0 + DENDRO_32 * grad_1_beta1 -
                   alpha[pp] * (-At1[pp] * K[pp] + DENDRO_1 * DENDRO_36 +
                                DENDRO_2 * DENDRO_37 + DENDRO_34 * DENDRO_5) +
                   beta0[pp] * grad_0_At1 + beta1[pp] * grad_1_At1 +
                   beta2[pp] * grad_2_At1;
    //--
    At_rhs02[pp] = At0[pp] * grad_2_beta0 + At1[pp] * grad_2_beta1 -
                   At2[pp] * DENDRO_39 + At4[pp] * grad_0_beta1 +
                   At5[pp] * grad_0_beta2 +
                   DENDRO_31 * (3 * DENDRO_22 + 3 * DENDRO_23 + 3 * DENDRO_24 +
                                3 * DENDRO_25 - DENDRO_30 * gt2[pp] -
                                3 * grad2_0_2_alpha) +
                   DENDRO_38 * grad_0_beta0 + DENDRO_38 * grad_2_beta2 -
                   alpha[pp] * (-At2[pp] * K[pp] + DENDRO_1 * DENDRO_41 +
                                DENDRO_2 * DENDRO_42 + DENDRO_40 * DENDRO_5) +
                   beta0[pp] * grad_0_At2 + beta1[pp] * grad_1_At2 +
                   beta2[pp] * grad_2_At2;
    //--
    At_rhs11[pp] = -At3[pp] * DENDRO_33 - At3[pp] * DENDRO_43 +
                   (4.0 / 3.0) * At3[pp] * grad_1_beta1 +
                   DENDRO_1 * grad_1_beta0 +
                   DENDRO_31 * (3 * DENDRO_10 + 3 * DENDRO_11 + 3 * DENDRO_12 +
                                3 * DENDRO_13 - DENDRO_30 * gt3[pp] -
                                3 * grad2_1_1_alpha) +
                   DENDRO_44 * grad_1_beta2 -
                   alpha[pp] * (-At3[pp] * K[pp] + DENDRO_1 * DENDRO_34 +
                                DENDRO_36 * DENDRO_45 + DENDRO_37 * DENDRO_44) +
                   beta0[pp] * grad_0_At3 + beta1[pp] * grad_1_At3 +
                   beta2[pp] * grad_2_At3;
    //--
    At_rhs12[pp] = At1[pp] * grad_2_beta0 + At2[pp] * grad_1_beta0 +
                   At3[pp] * grad_2_beta1 - At4[pp] * DENDRO_43 +
                   At5[pp] * grad_1_beta2 +
                   DENDRO_31 * (3 * DENDRO_26 + 3 * DENDRO_27 + 3 * DENDRO_28 +
                                3 * DENDRO_29 - DENDRO_30 * gt4[pp] -
                                3 * grad2_1_2_alpha) +
                   DENDRO_46 * grad_1_beta1 + DENDRO_46 * grad_2_beta2 -
                   alpha[pp] * (-At4[pp] * K[pp] + DENDRO_1 * DENDRO_40 +
                                DENDRO_41 * DENDRO_45 + DENDRO_42 * DENDRO_44) +
                   beta0[pp] * grad_0_At4 + beta1[pp] * grad_1_At4 +
                   beta2[pp] * grad_2_At4;
    //--
    At_rhs22[pp] = -At5[pp] * DENDRO_39 - At5[pp] * DENDRO_43 +
                   (4.0 / 3.0) * At5[pp] * grad_2_beta2 +
                   DENDRO_2 * grad_2_beta0 +
                   DENDRO_31 * (3 * DENDRO_14 + 3 * DENDRO_15 + 3 * DENDRO_16 +
                                3 * DENDRO_17 - DENDRO_30 * gt5[pp] -
                                3 * grad2_2_2_alpha) +
                   DENDRO_44 * grad_2_beta1 -
                   alpha[pp] * (2 * At5[pp] * DENDRO_42 - At5[pp] * K[pp] +
                                DENDRO_2 * DENDRO_40 + DENDRO_41 * DENDRO_44) +
                   beta0[pp] * grad_0_At5 + beta1[pp] * grad_1_At5 +
                   beta2[pp] * grad_2_At5;
    // Dendro: reduced ops: 362
    // Dendro: }}}
    if (blk->m_bflag != 0) {
        radiative_bc_pt<pw, nx>(&At_rhs00[gidx], At0[gidx], grad_0_At0,
                                grad_1_At0, grad_2_At0, 2.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&At_rhs01[gidx], At1[gidx], grad_0_At1,
                                grad_1_At1, grad_2_At1, 2.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&At_rhs02[gidx], At2[gidx], grad_0_At2,
                                grad_1_At2, grad_2_At2, 2.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&At_rhs11[gidx], At3[gidx], grad_0_At3,
                                grad_1_At3, grad_2_At3, 2.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&At_rhs12[gidx], At4[gidx], grad_0_At4,
                                grad_1_At4, grad_2_At4, 2.0, 0.0, blk);
        radiative_bc_pt<pw, nx>(&At_rhs22[gidx], At5[gidx], grad_0_At5,
                                grad_1_At5, grad_2_At5, 2.0, 0.0, blk);
    }
    At_rhs00[pp] += ko_sigma * (kograd_0_At0 + kograd_1_At0 + kograd_2_At0);
    At_rhs01[pp] += ko_sigma * (kograd_0_At1 + kograd_1_At1 + kograd_2_At1);
    At_rhs02[pp] += ko_sigma * (kograd_0_At2 + kograd_1_At2 + kograd_2_At2);
    At_rhs11[pp] += ko_sigma * (kograd_0_At3 + kograd_1_At3 + kograd_2_At3);
    At_rhs12[pp] += ko_sigma * (kograd_0_At4 + kograd_1_At4 + kograd_2_At4);
    At_rhs22[pp] += ko_sigma * (kograd_0_At5 + kograd_1_At5 + kograd_2_At5);
    __syncthreads();
}
