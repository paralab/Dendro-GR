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
    double DENDRO_0 = grad_2_alpha;
    double DENDRO_1 = beta2[pp];
    double DENDRO_2 = grad_1_alpha;
    double DENDRO_3 = beta1[pp];
    double DENDRO_4 = grad_0_alpha;
    double DENDRO_5 = beta0[pp];

    // initialize reduction for
    double DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_1;
    DENDRO_6 *= DENDRO_0;

    // initialize reduction for
    DENDRO_0 = 1;

    DENDRO_0 *= DENDRO_3;
    DENDRO_0 *= DENDRO_2;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_5;
    DENDRO_1 *= DENDRO_4;
    DENDRO_2 = alpha[pp];
    DENDRO_3 = K[pp];
    DENDRO_4 = -2;

    // initialize reduction for
    DENDRO_5 = 0;

    DENDRO_5 += DENDRO_1;
    DENDRO_5 += DENDRO_0;
    DENDRO_5 += DENDRO_6;
    DENDRO_0 = lambda[0];

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_4;
    DENDRO_1 *= DENDRO_3;
    DENDRO_1 *= DENDRO_2;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_0;
    DENDRO_2 *= DENDRO_5;

    // initialize reduction for
    DENDRO_0 = 0;

    DENDRO_0 += DENDRO_2;
    DENDRO_0 += DENDRO_1;
    a_rhs[pp] = DENDRO_0;
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
    double DENDRO_0  = grad_2_beta2;
    double DENDRO_1  = beta2[pp];
    double DENDRO_2  = grad_1_beta2;
    double DENDRO_3  = beta1[pp];
    double DENDRO_4  = grad_0_beta2;
    double DENDRO_5  = beta0[pp];
    double DENDRO_6  = lambda_f[1];
    double DENDRO_7  = alpha[pp];
    double DENDRO_8  = 3.0 / 4.0;
    double DENDRO_9  = lambda_f[0];
    double DENDRO_10 = grad_2_beta1;
    double DENDRO_11 = grad_1_beta1;
    double DENDRO_12 = grad_0_beta1;
    double DENDRO_13 = grad_2_beta0;
    double DENDRO_14 = grad_1_beta0;
    double DENDRO_15 = grad_0_beta0;

    // initialize reduction for
    double DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_1;
    DENDRO_16 *= DENDRO_0;

    // initialize reduction for
    DENDRO_0 = 1;

    DENDRO_0 *= DENDRO_3;
    DENDRO_0 *= DENDRO_2;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_5;
    DENDRO_2 *= DENDRO_4;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_8;
    DENDRO_4 *= DENDRO_7;
    DENDRO_4 *= DENDRO_6;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_8;
    DENDRO_6 *= DENDRO_9;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_1;
    DENDRO_7 *= DENDRO_10;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_3;
    DENDRO_8 *= DENDRO_11;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_5;
    DENDRO_9 *= DENDRO_12;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_1;
    DENDRO_10 *= DENDRO_13;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_3;
    DENDRO_1 *= DENDRO_14;

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_5;
    DENDRO_3 *= DENDRO_15;

    // initialize reduction for
    DENDRO_5 = 0;

    DENDRO_5 += DENDRO_2;
    DENDRO_5 += DENDRO_0;
    DENDRO_5 += DENDRO_16;
    DENDRO_0 = lambda[1];

    // initialize reduction for
    DENDRO_2 = 0;

    DENDRO_2 += DENDRO_6;
    DENDRO_2 += DENDRO_4;
    DENDRO_4 = B2[pp];

    // initialize reduction for
    DENDRO_6 = 0;

    DENDRO_6 += DENDRO_9;
    DENDRO_6 += DENDRO_8;
    DENDRO_6 += DENDRO_7;
    DENDRO_7 = B1[pp];

    // initialize reduction for
    DENDRO_8 = 0;

    DENDRO_8 += DENDRO_3;
    DENDRO_8 += DENDRO_1;
    DENDRO_8 += DENDRO_10;
    DENDRO_1 = B0[pp];

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_0;
    DENDRO_3 *= DENDRO_5;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_4;
    DENDRO_5 *= DENDRO_2;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_0;
    DENDRO_4 *= DENDRO_6;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_7;
    DENDRO_6 *= DENDRO_2;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_0;
    DENDRO_7 *= DENDRO_8;

    // initialize reduction for
    DENDRO_0 = 1;

    DENDRO_0 *= DENDRO_1;
    DENDRO_0 *= DENDRO_2;

    // initialize reduction for
    DENDRO_1 = 0;

    DENDRO_1 += DENDRO_5;
    DENDRO_1 += DENDRO_3;
    b_rhs2[pp] = DENDRO_1;

    // initialize reduction for
    DENDRO_1   = 0;

    DENDRO_1 += DENDRO_6;
    DENDRO_1 += DENDRO_4;
    b_rhs1[pp] = DENDRO_1;

    // initialize reduction for
    DENDRO_1   = 0;

    DENDRO_1 += DENDRO_0;
    DENDRO_1 += DENDRO_7;
    b_rhs0[pp] = DENDRO_1;
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
    double DENDRO_0 = grad_2_beta2;
    double DENDRO_1 = grad_1_beta1;
    double DENDRO_2 = grad_0_beta0;
    double DENDRO_3 = chi[pp];
    double DENDRO_4 = alpha[pp];
    double DENDRO_5 = K[pp];
    double DENDRO_6 = 2.0 / 3.0;

    // initialize reduction for
    double DENDRO_7 = 0;

    DENDRO_7 += DENDRO_2;
    DENDRO_7 += DENDRO_1;
    DENDRO_7 += DENDRO_0;
    DENDRO_0         = -2.0 / 3.0;
    DENDRO_1         = grad_2_chi;
    DENDRO_2         = beta2[pp];
    double DENDRO_8  = grad_1_chi;
    double DENDRO_9  = beta1[pp];
    double DENDRO_10 = grad_0_chi;
    double DENDRO_11 = beta0[pp];

    // initialize reduction for
    double DENDRO_12 = 1;

    DENDRO_12 *= DENDRO_6;
    DENDRO_12 *= DENDRO_5;
    DENDRO_12 *= DENDRO_4;
    DENDRO_12 *= DENDRO_3;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_0;
    DENDRO_4 *= DENDRO_3;
    DENDRO_4 *= DENDRO_7;

    // initialize reduction for
    DENDRO_0 = 1;

    DENDRO_0 *= DENDRO_2;
    DENDRO_0 *= DENDRO_1;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_9;
    DENDRO_1 *= DENDRO_8;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_11;
    DENDRO_2 *= DENDRO_10;

    // initialize reduction for
    DENDRO_3 = 0;

    DENDRO_3 += DENDRO_2;
    DENDRO_3 += DENDRO_1;
    DENDRO_3 += DENDRO_0;
    DENDRO_3 += DENDRO_4;
    DENDRO_3 += DENDRO_12;
    chi_rhs[pp] = DENDRO_3;
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
    double DENDRO_0  = gt5[pp];
    double DENDRO_1  = grad_2_beta2;
    double DENDRO_2  = 4.0 / 3.0;
    double DENDRO_3  = grad_1_beta1;
    double DENDRO_4  = -2.0 / 3.0;
    double DENDRO_5  = grad_0_beta0;
    double DENDRO_6  = gt4[pp];
    double DENDRO_7  = grad_2_beta1;
    double DENDRO_8  = 2;
    double DENDRO_9  = gt2[pp];
    double DENDRO_10 = grad_2_beta0;
    double DENDRO_11 = alpha[pp];
    double DENDRO_12 = At5[pp];
    double DENDRO_13 = -2;
    double DENDRO_14 = grad_2_gt5;
    double DENDRO_15 = beta2[pp];
    double DENDRO_16 = grad_1_gt5;
    double DENDRO_17 = beta1[pp];
    double DENDRO_18 = grad_0_gt5;
    double DENDRO_19 = beta0[pp];
    double DENDRO_20 = 1.0 / 3.0;
    double DENDRO_21 = At4[pp];
    double DENDRO_22 = gt3[pp];
    double DENDRO_23 = gt1[pp];
    double DENDRO_24 = grad_1_beta2;
    double DENDRO_25 = grad_1_beta0;
    double DENDRO_26 = grad_2_gt4;
    double DENDRO_27 = grad_1_gt4;
    double DENDRO_28 = grad_0_gt4;
    double DENDRO_29 = At3[pp];
    double DENDRO_30 = grad_2_gt3;
    double DENDRO_31 = grad_1_gt3;
    double DENDRO_32 = grad_0_gt3;
    double DENDRO_33 = At2[pp];
    double DENDRO_34 = gt0[pp];
    double DENDRO_35 = grad_0_beta2;
    double DENDRO_36 = grad_0_beta1;
    double DENDRO_37 = grad_2_gt2;
    double DENDRO_38 = grad_1_gt2;
    double DENDRO_39 = grad_0_gt2;
    double DENDRO_40 = At1[pp];
    double DENDRO_41 = grad_2_gt1;
    double DENDRO_42 = grad_1_gt1;
    double DENDRO_43 = grad_0_gt1;
    double DENDRO_44 = At0[pp];
    double DENDRO_45 = grad_2_gt0;
    double DENDRO_46 = grad_1_gt0;
    double DENDRO_47 = grad_0_gt0;

    // initialize reduction for
    double DENDRO_48 = 1;

    DENDRO_48 *= DENDRO_2;
    DENDRO_48 *= DENDRO_1;
    DENDRO_48 *= DENDRO_0;

    // initialize reduction for
    double DENDRO_49 = 1;

    DENDRO_49 *= DENDRO_4;
    DENDRO_49 *= DENDRO_3;
    DENDRO_49 *= DENDRO_0;

    // initialize reduction for
    double DENDRO_50 = 1;

    DENDRO_50 *= DENDRO_4;
    DENDRO_50 *= DENDRO_5;
    DENDRO_50 *= DENDRO_0;

    // initialize reduction for
    double DENDRO_51 = 1;

    DENDRO_51 *= DENDRO_8;
    DENDRO_51 *= DENDRO_7;
    DENDRO_51 *= DENDRO_6;

    // initialize reduction for
    double DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_8;
    DENDRO_52 *= DENDRO_10;
    DENDRO_52 *= DENDRO_9;

    // initialize reduction for
    double DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_13;
    DENDRO_53 *= DENDRO_12;
    DENDRO_53 *= DENDRO_11;

    // initialize reduction for
    DENDRO_12 = 1;

    DENDRO_12 *= DENDRO_15;
    DENDRO_12 *= DENDRO_14;

    // initialize reduction for
    DENDRO_14 = 1;

    DENDRO_14 *= DENDRO_17;
    DENDRO_14 *= DENDRO_16;

    // initialize reduction for
    DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_19;
    DENDRO_16 *= DENDRO_18;

    // initialize reduction for
    DENDRO_18 = 1;

    DENDRO_18 *= DENDRO_20;
    DENDRO_18 *= DENDRO_1;
    DENDRO_18 *= DENDRO_6;

    // initialize reduction for
    double DENDRO_54 = 1;

    DENDRO_54 *= DENDRO_20;
    DENDRO_54 *= DENDRO_3;
    DENDRO_54 *= DENDRO_6;

    // initialize reduction for
    double DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_4;
    DENDRO_55 *= DENDRO_5;
    DENDRO_55 *= DENDRO_6;

    // initialize reduction for
    double DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_13;
    DENDRO_56 *= DENDRO_21;
    DENDRO_56 *= DENDRO_11;

    // initialize reduction for
    DENDRO_21 = 1;

    DENDRO_21 *= DENDRO_7;
    DENDRO_21 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_57 = 1;

    DENDRO_57 *= DENDRO_10;
    DENDRO_57 *= DENDRO_23;

    // initialize reduction for
    double DENDRO_58 = 1;

    DENDRO_58 *= DENDRO_24;
    DENDRO_58 *= DENDRO_0;

    // initialize reduction for
    double DENDRO_59 = 1;

    DENDRO_59 *= DENDRO_25;
    DENDRO_59 *= DENDRO_9;

    // initialize reduction for
    double DENDRO_60 = 1;

    DENDRO_60 *= DENDRO_15;
    DENDRO_60 *= DENDRO_26;

    // initialize reduction for
    DENDRO_26 = 1;

    DENDRO_26 *= DENDRO_17;
    DENDRO_26 *= DENDRO_27;

    // initialize reduction for
    DENDRO_27 = 1;

    DENDRO_27 *= DENDRO_19;
    DENDRO_27 *= DENDRO_28;

    // initialize reduction for
    DENDRO_28 = 1;

    DENDRO_28 *= DENDRO_2;
    DENDRO_28 *= DENDRO_3;
    DENDRO_28 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_4;
    DENDRO_61 *= DENDRO_1;
    DENDRO_61 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_4;
    DENDRO_62 *= DENDRO_5;
    DENDRO_62 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_63 = 1;

    DENDRO_63 *= DENDRO_8;
    DENDRO_63 *= DENDRO_24;
    DENDRO_63 *= DENDRO_6;

    // initialize reduction for
    double DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_8;
    DENDRO_64 *= DENDRO_25;
    DENDRO_64 *= DENDRO_23;

    // initialize reduction for
    double DENDRO_65 = 1;

    DENDRO_65 *= DENDRO_13;
    DENDRO_65 *= DENDRO_29;
    DENDRO_65 *= DENDRO_11;

    // initialize reduction for
    DENDRO_29 = 1;

    DENDRO_29 *= DENDRO_15;
    DENDRO_29 *= DENDRO_30;

    // initialize reduction for
    DENDRO_30 = 1;

    DENDRO_30 *= DENDRO_17;
    DENDRO_30 *= DENDRO_31;

    // initialize reduction for
    DENDRO_31 = 1;

    DENDRO_31 *= DENDRO_19;
    DENDRO_31 *= DENDRO_32;

    // initialize reduction for
    DENDRO_32 = 1;

    DENDRO_32 *= DENDRO_20;
    DENDRO_32 *= DENDRO_1;
    DENDRO_32 *= DENDRO_9;

    // initialize reduction for
    double DENDRO_66 = 1;

    DENDRO_66 *= DENDRO_20;
    DENDRO_66 *= DENDRO_5;
    DENDRO_66 *= DENDRO_9;

    // initialize reduction for
    double DENDRO_67 = 1;

    DENDRO_67 *= DENDRO_4;
    DENDRO_67 *= DENDRO_3;
    DENDRO_67 *= DENDRO_9;

    // initialize reduction for
    double DENDRO_68 = 1;

    DENDRO_68 *= DENDRO_13;
    DENDRO_68 *= DENDRO_33;
    DENDRO_68 *= DENDRO_11;

    // initialize reduction for
    DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_7;
    DENDRO_33 *= DENDRO_23;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_10;
    DENDRO_7 *= DENDRO_34;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_35;
    DENDRO_10 *= DENDRO_0;

    // initialize reduction for
    DENDRO_0 = 1;

    DENDRO_0 *= DENDRO_36;
    DENDRO_0 *= DENDRO_6;

    // initialize reduction for
    double DENDRO_69 = 1;

    DENDRO_69 *= DENDRO_15;
    DENDRO_69 *= DENDRO_37;

    // initialize reduction for
    DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_17;
    DENDRO_37 *= DENDRO_38;

    // initialize reduction for
    DENDRO_38 = 1;

    DENDRO_38 *= DENDRO_19;
    DENDRO_38 *= DENDRO_39;

    // initialize reduction for
    DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_20;
    DENDRO_39 *= DENDRO_3;
    DENDRO_39 *= DENDRO_23;

    // initialize reduction for
    double DENDRO_70 = 1;

    DENDRO_70 *= DENDRO_20;
    DENDRO_70 *= DENDRO_5;
    DENDRO_70 *= DENDRO_23;

    // initialize reduction for
    DENDRO_20 = 1;

    DENDRO_20 *= DENDRO_4;
    DENDRO_20 *= DENDRO_1;
    DENDRO_20 *= DENDRO_23;

    // initialize reduction for
    double DENDRO_71 = 1;

    DENDRO_71 *= DENDRO_13;
    DENDRO_71 *= DENDRO_40;
    DENDRO_71 *= DENDRO_11;

    // initialize reduction for
    DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_24;
    DENDRO_40 *= DENDRO_9;

    // initialize reduction for
    DENDRO_24 = 1;

    DENDRO_24 *= DENDRO_25;
    DENDRO_24 *= DENDRO_34;

    // initialize reduction for
    DENDRO_25 = 1;

    DENDRO_25 *= DENDRO_35;
    DENDRO_25 *= DENDRO_6;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_36;
    DENDRO_6 *= DENDRO_22;

    // initialize reduction for
    DENDRO_22 = 1;

    DENDRO_22 *= DENDRO_15;
    DENDRO_22 *= DENDRO_41;

    // initialize reduction for
    DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_17;
    DENDRO_41 *= DENDRO_42;

    // initialize reduction for
    DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_19;
    DENDRO_42 *= DENDRO_43;

    // initialize reduction for
    DENDRO_43 = 1;

    DENDRO_43 *= DENDRO_2;
    DENDRO_43 *= DENDRO_5;
    DENDRO_43 *= DENDRO_34;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_4;
    DENDRO_2 *= DENDRO_1;
    DENDRO_2 *= DENDRO_34;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_4;
    DENDRO_1 *= DENDRO_3;
    DENDRO_1 *= DENDRO_34;

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_8;
    DENDRO_3 *= DENDRO_35;
    DENDRO_3 *= DENDRO_9;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_8;
    DENDRO_4 *= DENDRO_36;
    DENDRO_4 *= DENDRO_23;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_13;
    DENDRO_5 *= DENDRO_44;
    DENDRO_5 *= DENDRO_11;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_15;
    DENDRO_8 *= DENDRO_45;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_17;
    DENDRO_9 *= DENDRO_46;

    // initialize reduction for
    DENDRO_11 = 1;

    DENDRO_11 *= DENDRO_19;
    DENDRO_11 *= DENDRO_47;

    // initialize reduction for
    DENDRO_13 = 0;

    DENDRO_13 += DENDRO_16;
    DENDRO_13 += DENDRO_14;
    DENDRO_13 += DENDRO_12;
    DENDRO_13 += DENDRO_53;
    DENDRO_13 += DENDRO_52;
    DENDRO_13 += DENDRO_51;
    DENDRO_13 += DENDRO_50;
    DENDRO_13 += DENDRO_49;
    DENDRO_13 += DENDRO_48;
    gt_rhs22[pp] = DENDRO_13;

    // initialize reduction for
    DENDRO_12    = 0;

    DENDRO_12 += DENDRO_27;
    DENDRO_12 += DENDRO_26;
    DENDRO_12 += DENDRO_60;
    DENDRO_12 += DENDRO_59;
    DENDRO_12 += DENDRO_58;
    DENDRO_12 += DENDRO_57;
    DENDRO_12 += DENDRO_21;
    DENDRO_12 += DENDRO_56;
    DENDRO_12 += DENDRO_55;
    DENDRO_12 += DENDRO_54;
    DENDRO_12 += DENDRO_18;
    gt_rhs12[pp] = DENDRO_12;

    // initialize reduction for
    DENDRO_12    = 0;

    DENDRO_12 += DENDRO_31;
    DENDRO_12 += DENDRO_30;
    DENDRO_12 += DENDRO_29;
    DENDRO_12 += DENDRO_65;
    DENDRO_12 += DENDRO_64;
    DENDRO_12 += DENDRO_63;
    DENDRO_12 += DENDRO_62;
    DENDRO_12 += DENDRO_61;
    DENDRO_12 += DENDRO_28;
    gt_rhs11[pp] = DENDRO_12;

    // initialize reduction for
    DENDRO_12    = 0;

    DENDRO_12 += DENDRO_38;
    DENDRO_12 += DENDRO_37;
    DENDRO_12 += DENDRO_69;
    DENDRO_12 += DENDRO_0;
    DENDRO_12 += DENDRO_10;
    DENDRO_12 += DENDRO_7;
    DENDRO_12 += DENDRO_33;
    DENDRO_12 += DENDRO_68;
    DENDRO_12 += DENDRO_67;
    DENDRO_12 += DENDRO_66;
    DENDRO_12 += DENDRO_32;
    gt_rhs02[pp] = DENDRO_12;

    // initialize reduction for
    DENDRO_0     = 0;

    DENDRO_0 += DENDRO_42;
    DENDRO_0 += DENDRO_41;
    DENDRO_0 += DENDRO_22;
    DENDRO_0 += DENDRO_6;
    DENDRO_0 += DENDRO_25;
    DENDRO_0 += DENDRO_24;
    DENDRO_0 += DENDRO_40;
    DENDRO_0 += DENDRO_71;
    DENDRO_0 += DENDRO_20;
    DENDRO_0 += DENDRO_70;
    DENDRO_0 += DENDRO_39;
    gt_rhs01[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0     = 0;

    DENDRO_0 += DENDRO_11;
    DENDRO_0 += DENDRO_9;
    DENDRO_0 += DENDRO_8;
    DENDRO_0 += DENDRO_5;
    DENDRO_0 += DENDRO_4;
    DENDRO_0 += DENDRO_3;
    DENDRO_0 += DENDRO_1;
    DENDRO_0 += DENDRO_2;
    DENDRO_0 += DENDRO_43;
    gt_rhs00[pp] = DENDRO_0;
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
        double DENDRO_0 = 2;
        double DENDRO_1 = gt1[pp];
        double DENDRO_2 = gt2[pp];
        double DENDRO_3 = gt4[pp];
        double DENDRO_4 = -2;
        double DENDRO_5 = gt5[pp];
        double DENDRO_6 = gt3[pp];
        double DENDRO_7 = gt0[pp];
        double DENDRO_8 = -1;
        double DENDRO_9 = DENDRO_1;

        DENDRO_9 *= DENDRO_1;
        double DENDRO_10 = DENDRO_2;

        DENDRO_10 *= DENDRO_2;
        double DENDRO_11 = DENDRO_3;

        DENDRO_11 *= DENDRO_3;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_4;
        DENDRO_0 *= DENDRO_1;
        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_8;
        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_5;
        DENDRO_12 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_6;
        DENDRO_13 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_8;
        DENDRO_15 *= DENDRO_7;
        DENDRO_15 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_16 = 0;

        DENDRO_16 += DENDRO_14;
        DENDRO_16 += DENDRO_13;
        DENDRO_16 += DENDRO_12;
        DENDRO_16 += DENDRO_4;
        DENDRO_16 += DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_8;
        DENDRO_0 *= DENDRO_1;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_3;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_8;
        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_3;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_2;
        DENDRO_13 *= DENDRO_6;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_8;
        DENDRO_14 *= DENDRO_2;
        DENDRO_14 *= DENDRO_3;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_8;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_9;
        DENDRO_3 += DENDRO_15;
        DENDRO_5 = DENDRO_16;

        DENDRO_5 = 1 / DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_4;
        DENDRO_6 += DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 0;

        DENDRO_0 += DENDRO_10;
        DENDRO_0 += DENDRO_12;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_13;
        DENDRO_4 += DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 0;

        DENDRO_7 += DENDRO_2;
        DENDRO_7 += DENDRO_14;

        // initialize reduction for
        DENDRO_2 = 0;

        DENDRO_2 += DENDRO_11;
        DENDRO_2 += DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_3;
        DENDRO_igt5 = DENDRO_1;

        // initialize reduction for
        DENDRO_1    = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_6;
        DENDRO_igt4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1    = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_0;
        DENDRO_igt3 = DENDRO_1;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_4;
        DENDRO_igt2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_7;
        DENDRO_igt1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_2;
        DENDRO_igt0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_0_gt5;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_2_gt2;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_0_gt4;
        double DENDRO_5  = grad_2_gt1;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_1_gt2;
        double DENDRO_8  = grad_0_gt3;
        double DENDRO_9  = grad_1_gt1;

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;
        DENDRO_3         = grad_2_gt0;
        DENDRO_8         = grad_1_gt0;
        DENDRO_9         = grad_0_gt0;

        // initialize reduction for
        double DENDRO_11 = 0;

        DENDRO_11 += DENDRO_0;
        DENDRO_11 += DENDRO_10;
        DENDRO_C1_k0_5 = DENDRO_11;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_2;
        DENDRO_C1_k0_4 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_7;
        DENDRO_C1_k0_3 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_3;
        DENDRO_C1_k0_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_8;
        DENDRO_C1_k0_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_9;
        DENDRO_C1_k0_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_1_gt5;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_2_gt4;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_1_gt2;
        double DENDRO_5  = grad_2_gt1;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_0_gt4;
        double DENDRO_8  = grad_1_gt0;
        double DENDRO_9  = grad_0_gt1;

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;
        DENDRO_2         = grad_2_gt3;
        double DENDRO_11 = grad_1_gt3;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;
        DENDRO_7         = grad_0_gt3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_0;
        DENDRO_3 += DENDRO_10;
        DENDRO_C1_k1_5 = DENDRO_3;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_2;
        DENDRO_C1_k1_4 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_11;
        DENDRO_C1_k1_3 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_12;
        DENDRO_C1_k1_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_7;
        DENDRO_C1_k1_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_13;
        DENDRO_C1_k1_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_2_gt3;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_1_gt4;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_2_gt1;
        double DENDRO_5  = grad_1_gt2;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_0_gt4;
        double DENDRO_8  = grad_2_gt0;
        double DENDRO_9  = grad_0_gt2;
        double DENDRO_10 = grad_2_gt5;
        double DENDRO_11 = grad_1_gt5;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;
        DENDRO_2         = grad_0_gt5;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_6;
        DENDRO_3 *= DENDRO_10;
        DENDRO_C1_k2_5 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 1;

        DENDRO_3 *= DENDRO_6;
        DENDRO_3 *= DENDRO_11;
        DENDRO_C1_k2_4 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_0;
        DENDRO_3 += DENDRO_12;
        DENDRO_C1_k2_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_2;
        DENDRO_C1_k2_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_13;
        DENDRO_C1_k2_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_7;
        DENDRO_C1_k2_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = DENDRO_igt2;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt1;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt0;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k0_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k0_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k0_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k0_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k0_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k0_0 = DENDRO_1;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt3;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt1;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k1_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k1_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k1_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k1_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k1_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k1_0 = DENDRO_1;
    }
    {
        double DENDRO_0  = DENDRO_igt5;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt4;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt2;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k2_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k2_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k2_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k2_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k2_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k2_0 = DENDRO_1;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt2;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt1;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt0;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1         = gt2[pp];
        DENDRO_3         = -1;
        DENDRO_6         = gt1[pp];
        double DENDRO_7  = gt0[pp];
        double DENDRO_8  = 2;
        double DENDRO_9  = chi[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_3;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_3;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;
        DENDRO_4         = DENDRO_9;

        DENDRO_4         = 1 / DENDRO_4;
        DENDRO_3         = gt5[pp];
        DENDRO_8         = 0.5;
        DENDRO_9         = gt4[pp];
        double DENDRO_11 = gt3[pp];

        // initialize reduction for
        double DENDRO_12 = 0;

        DENDRO_12 += DENDRO_0;
        DENDRO_12 += DENDRO_10;
        DENDRO_0  = -0.5;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_2;
        DENDRO_10 += DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_8;
        DENDRO_2 *= DENDRO_3;
        DENDRO_2 *= DENDRO_4;
        DENDRO_2 *= DENDRO_5;
        DENDRO_3 = DENDRO_C2_k0_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_9;
        DENDRO_6 *= DENDRO_4;
        DENDRO_6 *= DENDRO_5;
        DENDRO_7 = DENDRO_C2_k0_4;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_8;
        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_4;
        DENDRO_9 *= DENDRO_5;
        DENDRO_5 = DENDRO_C2_k0_3;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_0;
        DENDRO_8 *= DENDRO_4;
        DENDRO_8 *= DENDRO_12;
        DENDRO_11 = DENDRO_C2_k0_2;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_0;
        DENDRO_12 *= DENDRO_4;
        DENDRO_12 *= DENDRO_10;
        DENDRO_10        = DENDRO_C2_k0_1;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_0;
        DENDRO_13 *= DENDRO_4;
        DENDRO_13 *= DENDRO_1;
        DENDRO_0 = DENDRO_C2_k0_0;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_3;
        DENDRO_1 += DENDRO_2;
        DENDRO_C3_k0_5 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_C3_k0_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_5;
        DENDRO_1 += DENDRO_9;
        DENDRO_C3_k0_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_11;
        DENDRO_1 += DENDRO_8;
        DENDRO_C3_k0_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_12;
        DENDRO_C3_k0_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_13;
        DENDRO_C3_k0_0 = DENDRO_1;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt4;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt3;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt1;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1         = gt4[pp];
        DENDRO_3         = -1;
        DENDRO_6         = gt3[pp];
        double DENDRO_7  = 2;
        double DENDRO_8  = gt1[pp];
        double DENDRO_9  = chi[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_3;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_3;
        DENDRO_2 *= DENDRO_8;
        DENDRO_2 *= DENDRO_5;
        DENDRO_7 = DENDRO_9;

        DENDRO_7 = 1 / DENDRO_7;
        DENDRO_3 = gt5[pp];
        DENDRO_8 = 0.5;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_0;
        DENDRO_9 += DENDRO_10;
        DENDRO_0  = -0.5;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_6;
        DENDRO_10 += DENDRO_1;
        DENDRO_1 = gt2[pp];

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_4;
        DENDRO_6 += DENDRO_2;
        DENDRO_2 = gt0[pp];

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_8;
        DENDRO_4 *= DENDRO_3;
        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_5;
        DENDRO_3         = DENDRO_C2_k1_5;

        // initialize reduction for
        double DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_0;
        DENDRO_11 *= DENDRO_7;
        DENDRO_11 *= DENDRO_9;
        DENDRO_9         = DENDRO_C2_k1_4;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_0;
        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_10;
        DENDRO_10        = DENDRO_C2_k1_3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_8;
        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_7;
        DENDRO_13 *= DENDRO_5;
        DENDRO_1         = DENDRO_C2_k1_2;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_0;
        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_6;
        DENDRO_0 = DENDRO_C2_k1_1;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_2;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;
        DENDRO_2 = DENDRO_C2_k1_0;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_4;
        DENDRO_C3_k1_5 = DENDRO_5;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_9;
        DENDRO_3 += DENDRO_11;
        DENDRO_C3_k1_4 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_10;
        DENDRO_3 += DENDRO_12;
        DENDRO_C3_k1_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_1;
        DENDRO_3 += DENDRO_13;
        DENDRO_C3_k1_2 = DENDRO_3;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_14;
        DENDRO_C3_k1_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_2;
        DENDRO_0 += DENDRO_6;
        DENDRO_C3_k1_0 = DENDRO_0;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt5;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt4;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt2;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1        = gt5[pp];
        DENDRO_3        = -1;
        DENDRO_6        = 2;
        double DENDRO_7 = gt4[pp];
        double DENDRO_8 = gt2[pp];

        // initialize reduction for
        double DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_3;
        DENDRO_9 *= DENDRO_1;
        DENDRO_9 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_0;
        DENDRO_0 = chi[pp];

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_3;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_3;
        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_5;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_1;
        DENDRO_8 += DENDRO_9;
        DENDRO_1 = DENDRO_0;

        DENDRO_1 = 1 / DENDRO_1;
        DENDRO_0 = -0.5;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_2;
        DENDRO_3 += DENDRO_6;
        DENDRO_2 = gt3[pp];
        DENDRO_6 = 0.5;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_4;
        DENDRO_9 += DENDRO_7;
        DENDRO_4         = gt1[pp];
        DENDRO_7         = gt0[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_0;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_8;
        DENDRO_8         = DENDRO_C2_k2_5;

        // initialize reduction for
        double DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_0;
        DENDRO_11 *= DENDRO_1;
        DENDRO_11 *= DENDRO_3;
        DENDRO_3         = DENDRO_C2_k2_4;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_6;
        DENDRO_12 *= DENDRO_2;
        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_5;
        DENDRO_2         = DENDRO_C2_k2_3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_0;
        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_9;
        DENDRO_0 = DENDRO_C2_k2_2;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_6;
        DENDRO_9 *= DENDRO_4;
        DENDRO_9 *= DENDRO_1;
        DENDRO_9 *= DENDRO_5;
        DENDRO_4         = DENDRO_C2_k2_1;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_6;
        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_1;
        DENDRO_14 *= DENDRO_5;
        DENDRO_1 = DENDRO_C2_k2_0;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_8;
        DENDRO_5 += DENDRO_10;
        DENDRO_C3_k2_5 = DENDRO_5;

        // initialize reduction for
        DENDRO_5       = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_11;
        DENDRO_C3_k2_4 = DENDRO_5;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_2;
        DENDRO_3 += DENDRO_12;
        DENDRO_C3_k2_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_2       = 0;

        DENDRO_2 += DENDRO_0;
        DENDRO_2 += DENDRO_13;
        DENDRO_C3_k2_2 = DENDRO_2;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_9;
        DENDRO_C3_k2_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_14;
        DENDRO_C3_k2_0 = DENDRO_0;
    }
    double DENDRO_0  = 2;
    double DENDRO_1  = DENDRO_igt4;
    double DENDRO_2  = DENDRO_igt2;
    double DENDRO_3  = DENDRO_igt1;
    double DENDRO_4  = DENDRO_igt5;
    double DENDRO_5  = DENDRO_igt3;
    double DENDRO_6  = DENDRO_igt0;
    double DENDRO_7  = gt1[pp];
    double DENDRO_8  = gt2[pp];
    double DENDRO_9  = gt4[pp];
    double DENDRO_10 = At5[pp];
    double DENDRO_11 = At4[pp];
    double DENDRO_12 = At3[pp];
    double DENDRO_13 = At2[pp];
    double DENDRO_14 = At1[pp];
    double DENDRO_15 = At0[pp];
    double DENDRO_16 = DENDRO_1;

    DENDRO_16 *= DENDRO_1;
    double DENDRO_17 = DENDRO_2;

    DENDRO_17 *= DENDRO_2;
    double DENDRO_18 = DENDRO_3;

    DENDRO_18 *= DENDRO_3;
    double DENDRO_19 = DENDRO_4;

    DENDRO_19 *= DENDRO_4;
    double DENDRO_20 = DENDRO_5;

    DENDRO_20 *= DENDRO_5;
    double DENDRO_21 = DENDRO_6;

    DENDRO_21 *= DENDRO_6;
    double DENDRO_22 = gt5[pp];
    double DENDRO_23 = gt3[pp];
    double DENDRO_24 = gt0[pp];
    double DENDRO_25 = DENDRO_7;

    DENDRO_25 *= DENDRO_7;
    double DENDRO_26 = -1;
    double DENDRO_27 = DENDRO_8;

    DENDRO_27 *= DENDRO_8;
    double DENDRO_28 = DENDRO_9;

    DENDRO_28 *= DENDRO_9;
    double DENDRO_29 = -2;

    // initialize reduction for
    double DENDRO_30 = 1;

    DENDRO_30 *= DENDRO_10;
    DENDRO_30 *= DENDRO_1;
    DENDRO_30 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_31 = 1;

    DENDRO_31 *= DENDRO_11;
    DENDRO_31 *= DENDRO_5;
    DENDRO_31 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_32 = 1;

    DENDRO_32 *= DENDRO_12;
    DENDRO_32 *= DENDRO_5;
    DENDRO_32 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_13;
    DENDRO_33 *= DENDRO_2;
    DENDRO_33 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_34 = 1;

    DENDRO_34 *= DENDRO_13;
    DENDRO_34 *= DENDRO_3;
    DENDRO_34 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_35 = 1;

    DENDRO_35 *= DENDRO_14;
    DENDRO_35 *= DENDRO_2;
    DENDRO_35 *= DENDRO_5;

    // initialize reduction for
    double DENDRO_36 = 1;

    DENDRO_36 *= DENDRO_14;
    DENDRO_36 *= DENDRO_3;
    DENDRO_36 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_15;
    DENDRO_37 *= DENDRO_3;
    DENDRO_37 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_38 = 1;

    DENDRO_38 *= DENDRO_11;
    DENDRO_38 *= DENDRO_16;

    // initialize reduction for
    double DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_10;
    DENDRO_39 *= DENDRO_2;
    DENDRO_39 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_11;
    DENDRO_40 *= DENDRO_2;
    DENDRO_40 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_11;
    DENDRO_41 *= DENDRO_3;
    DENDRO_41 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_12;
    DENDRO_42 *= DENDRO_3;
    DENDRO_42 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_43 = 1;

    DENDRO_43 *= DENDRO_13;
    DENDRO_43 *= DENDRO_6;
    DENDRO_43 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_44 = 1;

    DENDRO_44 *= DENDRO_14;
    DENDRO_44 *= DENDRO_3;
    DENDRO_44 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_45 = 1;

    DENDRO_45 *= DENDRO_14;
    DENDRO_45 *= DENDRO_6;
    DENDRO_45 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_46 = 1;

    DENDRO_46 *= DENDRO_15;
    DENDRO_46 *= DENDRO_6;
    DENDRO_46 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_47 = 1;

    DENDRO_47 *= DENDRO_13;
    DENDRO_47 *= DENDRO_17;

    // initialize reduction for
    double DENDRO_48 = 1;

    DENDRO_48 *= DENDRO_10;
    DENDRO_48 *= DENDRO_2;
    DENDRO_48 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_49 = 1;

    DENDRO_49 *= DENDRO_11;
    DENDRO_49 *= DENDRO_2;
    DENDRO_49 *= DENDRO_5;

    // initialize reduction for
    double DENDRO_50 = 1;

    DENDRO_50 *= DENDRO_11;
    DENDRO_50 *= DENDRO_3;
    DENDRO_50 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_51 = 1;

    DENDRO_51 *= DENDRO_12;
    DENDRO_51 *= DENDRO_3;
    DENDRO_51 *= DENDRO_5;

    // initialize reduction for
    double DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_13;
    DENDRO_52 *= DENDRO_3;
    DENDRO_52 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_13;
    DENDRO_53 *= DENDRO_6;
    DENDRO_53 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_54 = 1;

    DENDRO_54 *= DENDRO_14;
    DENDRO_54 *= DENDRO_6;
    DENDRO_54 *= DENDRO_5;

    // initialize reduction for
    double DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_15;
    DENDRO_55 *= DENDRO_6;
    DENDRO_55 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_14;
    DENDRO_56 *= DENDRO_18;

    // initialize reduction for
    double DENDRO_57 = 1;

    DENDRO_57 *= DENDRO_0;
    DENDRO_57 *= DENDRO_11;
    DENDRO_57 *= DENDRO_1;
    DENDRO_57 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_58 = 1;

    DENDRO_58 *= DENDRO_0;
    DENDRO_58 *= DENDRO_13;
    DENDRO_58 *= DENDRO_2;
    DENDRO_58 *= DENDRO_4;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_0;
    DENDRO_4 *= DENDRO_14;
    DENDRO_4 *= DENDRO_2;
    DENDRO_4 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_59 = 1;

    DENDRO_59 *= DENDRO_10;
    DENDRO_59 *= DENDRO_19;

    // initialize reduction for
    DENDRO_19 = 1;

    DENDRO_19 *= DENDRO_12;
    DENDRO_19 *= DENDRO_16;

    // initialize reduction for
    double DENDRO_60 = 1;

    DENDRO_60 *= DENDRO_15;
    DENDRO_60 *= DENDRO_17;

    // initialize reduction for
    double DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_0;
    DENDRO_61 *= DENDRO_11;
    DENDRO_61 *= DENDRO_5;
    DENDRO_61 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_0;
    DENDRO_62 *= DENDRO_13;
    DENDRO_62 *= DENDRO_3;
    DENDRO_62 *= DENDRO_1;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_0;
    DENDRO_1 *= DENDRO_14;
    DENDRO_1 *= DENDRO_3;
    DENDRO_1 *= DENDRO_5;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_10;
    DENDRO_5 *= DENDRO_16;

    // initialize reduction for
    DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_12;
    DENDRO_16 *= DENDRO_20;

    // initialize reduction for
    DENDRO_20 = 1;

    DENDRO_20 *= DENDRO_15;
    DENDRO_20 *= DENDRO_18;

    // initialize reduction for
    double DENDRO_63 = 1;

    DENDRO_63 *= DENDRO_0;
    DENDRO_63 *= DENDRO_11;
    DENDRO_63 *= DENDRO_3;
    DENDRO_63 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_0;
    DENDRO_64 *= DENDRO_13;
    DENDRO_64 *= DENDRO_6;
    DENDRO_64 *= DENDRO_2;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_0;
    DENDRO_2 *= DENDRO_14;
    DENDRO_2 *= DENDRO_6;
    DENDRO_2 *= DENDRO_3;

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_10;
    DENDRO_3 *= DENDRO_17;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_12;
    DENDRO_6 *= DENDRO_18;

    // initialize reduction for
    DENDRO_17 = 1;

    DENDRO_17 *= DENDRO_15;
    DENDRO_17 *= DENDRO_21;
    DENDRO_18        = K[pp];
    DENDRO_21        = grad_2_alpha;
    double DENDRO_65 = DENDRO_C3_k2_2;
    double DENDRO_66 = grad_1_alpha;
    double DENDRO_67 = DENDRO_C3_k1_2;
    double DENDRO_68 = grad_0_alpha;
    double DENDRO_69 = DENDRO_C3_k0_2;

    // initialize reduction for
    double DENDRO_70 = 1;

    DENDRO_70 *= DENDRO_0;
    DENDRO_70 *= DENDRO_7;
    DENDRO_70 *= DENDRO_8;
    DENDRO_70 *= DENDRO_9;

    // initialize reduction for
    double DENDRO_71 = 1;

    DENDRO_71 *= DENDRO_24;
    DENDRO_71 *= DENDRO_23;
    DENDRO_71 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_72 = 1;

    DENDRO_72 *= DENDRO_26;
    DENDRO_72 *= DENDRO_22;
    DENDRO_72 *= DENDRO_25;

    // initialize reduction for
    double DENDRO_73 = 1;

    DENDRO_73 *= DENDRO_26;
    DENDRO_73 *= DENDRO_23;
    DENDRO_73 *= DENDRO_27;

    // initialize reduction for
    double DENDRO_74 = 1;

    DENDRO_74 *= DENDRO_26;
    DENDRO_74 *= DENDRO_24;
    DENDRO_74 *= DENDRO_28;
    double DENDRO_75 = DENDRO_C3_k2_1;
    double DENDRO_76 = DENDRO_C3_k1_1;
    double DENDRO_77 = DENDRO_C3_k0_1;

    // initialize reduction for
    double DENDRO_78 = 1;

    DENDRO_78 *= DENDRO_29;
    DENDRO_78 *= DENDRO_7;
    DENDRO_78 *= DENDRO_8;
    DENDRO_78 *= DENDRO_9;

    // initialize reduction for
    double DENDRO_79 = 1;

    DENDRO_79 *= DENDRO_26;
    DENDRO_79 *= DENDRO_24;
    DENDRO_79 *= DENDRO_23;
    DENDRO_79 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_80 = 1;

    DENDRO_80 *= DENDRO_22;
    DENDRO_80 *= DENDRO_25;

    // initialize reduction for
    double DENDRO_81 = 1;

    DENDRO_81 *= DENDRO_23;
    DENDRO_81 *= DENDRO_27;

    // initialize reduction for
    double DENDRO_82 = 1;

    DENDRO_82 *= DENDRO_24;
    DENDRO_82 *= DENDRO_28;
    double DENDRO_83 = DENDRO_C3_k2_4;
    double DENDRO_84 = DENDRO_C3_k1_4;
    double DENDRO_85 = DENDRO_C3_k0_4;
    double DENDRO_86 = DENDRO_C3_k2_0;
    double DENDRO_87 = DENDRO_C3_k1_0;
    double DENDRO_88 = DENDRO_C3_k0_0;
    double DENDRO_89 = DENDRO_C3_k2_3;
    double DENDRO_90 = DENDRO_C3_k1_3;
    double DENDRO_91 = DENDRO_C3_k0_3;
    double DENDRO_92 = DENDRO_C3_k2_5;
    double DENDRO_93 = DENDRO_C3_k1_5;
    double DENDRO_94 = DENDRO_C3_k0_5;

    // initialize reduction for
    double DENDRO_95 = 0;

    DENDRO_95 += DENDRO_38;
    DENDRO_95 += DENDRO_37;
    DENDRO_95 += DENDRO_36;
    DENDRO_95 += DENDRO_35;
    DENDRO_95 += DENDRO_34;
    DENDRO_95 += DENDRO_33;
    DENDRO_95 += DENDRO_32;
    DENDRO_95 += DENDRO_31;
    DENDRO_95 += DENDRO_30;

    // initialize reduction for
    DENDRO_30 = 0;

    DENDRO_30 += DENDRO_47;
    DENDRO_30 += DENDRO_46;
    DENDRO_30 += DENDRO_45;
    DENDRO_30 += DENDRO_44;
    DENDRO_30 += DENDRO_43;
    DENDRO_30 += DENDRO_42;
    DENDRO_30 += DENDRO_41;
    DENDRO_30 += DENDRO_40;
    DENDRO_30 += DENDRO_39;

    // initialize reduction for
    DENDRO_31 = 0;

    DENDRO_31 += DENDRO_56;
    DENDRO_31 += DENDRO_55;
    DENDRO_31 += DENDRO_54;
    DENDRO_31 += DENDRO_53;
    DENDRO_31 += DENDRO_52;
    DENDRO_31 += DENDRO_51;
    DENDRO_31 += DENDRO_50;
    DENDRO_31 += DENDRO_49;
    DENDRO_31 += DENDRO_48;

    // initialize reduction for
    DENDRO_32 = 0;

    DENDRO_32 += DENDRO_60;
    DENDRO_32 += DENDRO_19;
    DENDRO_32 += DENDRO_59;
    DENDRO_32 += DENDRO_4;
    DENDRO_32 += DENDRO_58;
    DENDRO_32 += DENDRO_57;

    // initialize reduction for
    DENDRO_4 = 0;

    DENDRO_4 += DENDRO_20;
    DENDRO_4 += DENDRO_16;
    DENDRO_4 += DENDRO_5;
    DENDRO_4 += DENDRO_1;
    DENDRO_4 += DENDRO_62;
    DENDRO_4 += DENDRO_61;

    // initialize reduction for
    DENDRO_1 = 0;

    DENDRO_1 += DENDRO_17;
    DENDRO_1 += DENDRO_6;
    DENDRO_1 += DENDRO_3;
    DENDRO_1 += DENDRO_2;
    DENDRO_1 += DENDRO_64;
    DENDRO_1 += DENDRO_63;
    DENDRO_2 = DENDRO_18;

    DENDRO_2 *= DENDRO_18;
    DENDRO_3 = 1.0 / 3.0;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_26;
    DENDRO_5 *= DENDRO_65;
    DENDRO_5 *= DENDRO_21;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_26;
    DENDRO_6 *= DENDRO_67;
    DENDRO_6 *= DENDRO_66;

    // initialize reduction for
    DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_26;
    DENDRO_16 *= DENDRO_69;
    DENDRO_16 *= DENDRO_68;
    DENDRO_17 = grad2_0_2_alpha;

    // initialize reduction for
    DENDRO_18 = 1;

    DENDRO_18 *= DENDRO_26;
    DENDRO_18 *= DENDRO_8;
    DENDRO_18 *= DENDRO_23;

    // initialize reduction for
    DENDRO_19 = 1;

    DENDRO_19 *= DENDRO_7;
    DENDRO_19 *= DENDRO_9;

    // initialize reduction for
    DENDRO_20 = 0;

    DENDRO_20 += DENDRO_74;
    DENDRO_20 += DENDRO_73;
    DENDRO_20 += DENDRO_72;
    DENDRO_20 += DENDRO_71;
    DENDRO_20 += DENDRO_70;

    // initialize reduction for
    DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_26;
    DENDRO_33 *= DENDRO_75;
    DENDRO_33 *= DENDRO_21;

    // initialize reduction for
    DENDRO_34 = 1;

    DENDRO_34 *= DENDRO_26;
    DENDRO_34 *= DENDRO_76;
    DENDRO_34 *= DENDRO_66;

    // initialize reduction for
    DENDRO_35 = 1;

    DENDRO_35 *= DENDRO_26;
    DENDRO_35 *= DENDRO_77;
    DENDRO_35 *= DENDRO_68;
    DENDRO_36 = grad2_0_1_alpha;

    // initialize reduction for
    DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_26;
    DENDRO_37 *= DENDRO_8;
    DENDRO_37 *= DENDRO_9;

    // initialize reduction for
    DENDRO_38 = 1;

    DENDRO_38 *= DENDRO_7;
    DENDRO_38 *= DENDRO_22;

    // initialize reduction for
    DENDRO_39 = 0;

    DENDRO_39 += DENDRO_82;
    DENDRO_39 += DENDRO_81;
    DENDRO_39 += DENDRO_80;
    DENDRO_39 += DENDRO_79;
    DENDRO_39 += DENDRO_78;

    // initialize reduction for
    DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_26;
    DENDRO_40 *= DENDRO_83;
    DENDRO_40 *= DENDRO_21;

    // initialize reduction for
    DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_26;
    DENDRO_41 *= DENDRO_84;
    DENDRO_41 *= DENDRO_66;

    // initialize reduction for
    DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_26;
    DENDRO_42 *= DENDRO_85;
    DENDRO_42 *= DENDRO_68;
    DENDRO_43 = grad2_1_2_alpha;

    // initialize reduction for
    DENDRO_44 = 1;

    DENDRO_44 *= DENDRO_26;
    DENDRO_44 *= DENDRO_7;
    DENDRO_44 *= DENDRO_8;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_24;
    DENDRO_7 *= DENDRO_9;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_26;
    DENDRO_8 *= DENDRO_86;
    DENDRO_8 *= DENDRO_21;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_26;
    DENDRO_9 *= DENDRO_87;
    DENDRO_9 *= DENDRO_66;

    // initialize reduction for
    DENDRO_45 = 1;

    DENDRO_45 *= DENDRO_26;
    DENDRO_45 *= DENDRO_88;
    DENDRO_45 *= DENDRO_68;
    DENDRO_46 = grad2_0_0_alpha;

    // initialize reduction for
    DENDRO_47 = 1;

    DENDRO_47 *= DENDRO_23;
    DENDRO_47 *= DENDRO_22;

    // initialize reduction for
    DENDRO_48 = 1;

    DENDRO_48 *= DENDRO_26;
    DENDRO_48 *= DENDRO_28;

    // initialize reduction for
    DENDRO_28 = 1;

    DENDRO_28 *= DENDRO_26;
    DENDRO_28 *= DENDRO_89;
    DENDRO_28 *= DENDRO_21;

    // initialize reduction for
    DENDRO_49 = 1;

    DENDRO_49 *= DENDRO_26;
    DENDRO_49 *= DENDRO_90;
    DENDRO_49 *= DENDRO_66;

    // initialize reduction for
    DENDRO_50 = 1;

    DENDRO_50 *= DENDRO_26;
    DENDRO_50 *= DENDRO_91;
    DENDRO_50 *= DENDRO_68;
    DENDRO_51 = grad2_1_1_alpha;

    // initialize reduction for
    DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_24;
    DENDRO_52 *= DENDRO_22;

    // initialize reduction for
    DENDRO_22 = 1;

    DENDRO_22 *= DENDRO_26;
    DENDRO_22 *= DENDRO_27;

    // initialize reduction for
    DENDRO_27 = 1;

    DENDRO_27 *= DENDRO_26;
    DENDRO_27 *= DENDRO_92;
    DENDRO_27 *= DENDRO_21;

    // initialize reduction for
    DENDRO_21 = 1;

    DENDRO_21 *= DENDRO_26;
    DENDRO_21 *= DENDRO_93;
    DENDRO_21 *= DENDRO_66;

    // initialize reduction for
    DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_26;
    DENDRO_53 *= DENDRO_94;
    DENDRO_53 *= DENDRO_68;
    DENDRO_54 = grad2_2_2_alpha;

    // initialize reduction for
    DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_24;
    DENDRO_55 *= DENDRO_23;

    // initialize reduction for
    DENDRO_23 = 1;

    DENDRO_23 *= DENDRO_26;
    DENDRO_23 *= DENDRO_25;

    // initialize reduction for
    DENDRO_24 = 1;

    DENDRO_24 *= DENDRO_0;
    DENDRO_24 *= DENDRO_11;
    DENDRO_24 *= DENDRO_95;

    // initialize reduction for
    DENDRO_11 = 1;

    DENDRO_11 *= DENDRO_0;
    DENDRO_11 *= DENDRO_13;
    DENDRO_11 *= DENDRO_30;

    // initialize reduction for
    DENDRO_13 = 1;

    DENDRO_13 *= DENDRO_0;
    DENDRO_13 *= DENDRO_14;
    DENDRO_13 *= DENDRO_31;

    // initialize reduction for
    DENDRO_0 = 1;

    DENDRO_0 *= DENDRO_10;
    DENDRO_0 *= DENDRO_32;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_12;
    DENDRO_10 *= DENDRO_4;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_15;
    DENDRO_4 *= DENDRO_1;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_3;
    DENDRO_1 *= DENDRO_2;

    // initialize reduction for
    DENDRO_2 = 0;

    DENDRO_2 += DENDRO_17;
    DENDRO_2 += DENDRO_16;
    DENDRO_2 += DENDRO_6;
    DENDRO_2 += DENDRO_5;

    // initialize reduction for
    DENDRO_3 = 0;

    DENDRO_3 += DENDRO_19;
    DENDRO_3 += DENDRO_18;
    DENDRO_5  = DENDRO_20;

    DENDRO_5  = 1 / DENDRO_5;
    DENDRO_6  = chi[pp];

    // initialize reduction for
    DENDRO_12 = 0;

    DENDRO_12 += DENDRO_36;
    DENDRO_12 += DENDRO_35;
    DENDRO_12 += DENDRO_34;
    DENDRO_12 += DENDRO_33;

    // initialize reduction for
    DENDRO_14 = 0;

    DENDRO_14 += DENDRO_38;
    DENDRO_14 += DENDRO_37;
    DENDRO_15 = DENDRO_39;

    DENDRO_15 = 1 / DENDRO_15;

    // initialize reduction for
    DENDRO_16 = 0;

    DENDRO_16 += DENDRO_43;
    DENDRO_16 += DENDRO_42;
    DENDRO_16 += DENDRO_41;
    DENDRO_16 += DENDRO_40;

    // initialize reduction for
    DENDRO_17 = 0;

    DENDRO_17 += DENDRO_7;
    DENDRO_17 += DENDRO_44;

    // initialize reduction for
    DENDRO_7 = 0;

    DENDRO_7 += DENDRO_46;
    DENDRO_7 += DENDRO_45;
    DENDRO_7 += DENDRO_9;
    DENDRO_7 += DENDRO_8;

    // initialize reduction for
    DENDRO_8 = 0;

    DENDRO_8 += DENDRO_48;
    DENDRO_8 += DENDRO_47;

    // initialize reduction for
    DENDRO_9 = 0;

    DENDRO_9 += DENDRO_51;
    DENDRO_9 += DENDRO_50;
    DENDRO_9 += DENDRO_49;
    DENDRO_9 += DENDRO_28;

    // initialize reduction for
    DENDRO_18 = 0;

    DENDRO_18 += DENDRO_22;
    DENDRO_18 += DENDRO_52;

    // initialize reduction for
    DENDRO_19 = 0;

    DENDRO_19 += DENDRO_54;
    DENDRO_19 += DENDRO_53;
    DENDRO_19 += DENDRO_21;
    DENDRO_19 += DENDRO_27;

    // initialize reduction for
    DENDRO_20 = 0;

    DENDRO_20 += DENDRO_23;
    DENDRO_20 += DENDRO_55;
    DENDRO_21 = grad_2_K;
    DENDRO_22 = beta2[pp];
    DENDRO_23 = grad_1_K;
    DENDRO_25 = beta1[pp];
    DENDRO_27 = grad_0_K;
    DENDRO_28 = beta0[pp];

    // initialize reduction for
    DENDRO_30 = 0;

    DENDRO_30 += DENDRO_1;
    DENDRO_30 += DENDRO_4;
    DENDRO_30 += DENDRO_10;
    DENDRO_30 += DENDRO_0;
    DENDRO_30 += DENDRO_13;
    DENDRO_30 += DENDRO_11;
    DENDRO_30 += DENDRO_24;
    DENDRO_0 = alpha[pp];

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_29;
    DENDRO_1 *= DENDRO_6;
    DENDRO_1 *= DENDRO_5;
    DENDRO_1 *= DENDRO_3;
    DENDRO_1 *= DENDRO_2;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_29;
    DENDRO_2 *= DENDRO_6;
    DENDRO_2 *= DENDRO_15;
    DENDRO_2 *= DENDRO_14;
    DENDRO_2 *= DENDRO_12;

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_29;
    DENDRO_3 *= DENDRO_6;
    DENDRO_3 *= DENDRO_15;
    DENDRO_3 *= DENDRO_17;
    DENDRO_3 *= DENDRO_16;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_26;
    DENDRO_4 *= DENDRO_6;
    DENDRO_4 *= DENDRO_5;
    DENDRO_4 *= DENDRO_8;
    DENDRO_4 *= DENDRO_7;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_26;
    DENDRO_7 *= DENDRO_6;
    DENDRO_7 *= DENDRO_5;
    DENDRO_7 *= DENDRO_18;
    DENDRO_7 *= DENDRO_9;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_26;
    DENDRO_8 *= DENDRO_6;
    DENDRO_8 *= DENDRO_5;
    DENDRO_8 *= DENDRO_20;
    DENDRO_8 *= DENDRO_19;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_22;
    DENDRO_5 *= DENDRO_21;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_25;
    DENDRO_6 *= DENDRO_23;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_28;
    DENDRO_9 *= DENDRO_27;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_0;
    DENDRO_10 *= DENDRO_30;

    // initialize reduction for
    DENDRO_0 = 0;

    DENDRO_0 += DENDRO_10;
    DENDRO_0 += DENDRO_9;
    DENDRO_0 += DENDRO_6;
    DENDRO_0 += DENDRO_5;
    DENDRO_0 += DENDRO_8;
    DENDRO_0 += DENDRO_7;
    DENDRO_0 += DENDRO_4;
    DENDRO_0 += DENDRO_3;
    DENDRO_0 += DENDRO_2;
    DENDRO_0 += DENDRO_1;
    K_rhs[pp] = DENDRO_0;
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
        double DENDRO_0 = 2;
        double DENDRO_1 = gt1[pp];
        double DENDRO_2 = gt2[pp];
        double DENDRO_3 = gt4[pp];
        double DENDRO_4 = -2;
        double DENDRO_5 = gt5[pp];
        double DENDRO_6 = gt3[pp];
        double DENDRO_7 = gt0[pp];
        double DENDRO_8 = -1;
        double DENDRO_9 = DENDRO_1;

        DENDRO_9 *= DENDRO_1;
        double DENDRO_10 = DENDRO_2;

        DENDRO_10 *= DENDRO_2;
        double DENDRO_11 = DENDRO_3;

        DENDRO_11 *= DENDRO_3;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_4;
        DENDRO_0 *= DENDRO_1;
        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_8;
        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_5;
        DENDRO_12 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_6;
        DENDRO_13 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_8;
        DENDRO_15 *= DENDRO_7;
        DENDRO_15 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_16 = 0;

        DENDRO_16 += DENDRO_14;
        DENDRO_16 += DENDRO_13;
        DENDRO_16 += DENDRO_12;
        DENDRO_16 += DENDRO_4;
        DENDRO_16 += DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_8;
        DENDRO_0 *= DENDRO_1;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_3;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_8;
        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_3;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_2;
        DENDRO_13 *= DENDRO_6;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_8;
        DENDRO_14 *= DENDRO_2;
        DENDRO_14 *= DENDRO_3;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_8;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_9;
        DENDRO_3 += DENDRO_15;
        DENDRO_5 = DENDRO_16;

        DENDRO_5 = 1 / DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_4;
        DENDRO_6 += DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 0;

        DENDRO_0 += DENDRO_10;
        DENDRO_0 += DENDRO_12;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_13;
        DENDRO_4 += DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 0;

        DENDRO_7 += DENDRO_2;
        DENDRO_7 += DENDRO_14;

        // initialize reduction for
        DENDRO_2 = 0;

        DENDRO_2 += DENDRO_11;
        DENDRO_2 += DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_3;
        DENDRO_igt5 = DENDRO_1;

        // initialize reduction for
        DENDRO_1    = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_6;
        DENDRO_igt4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1    = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_0;
        DENDRO_igt3 = DENDRO_1;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_4;
        DENDRO_igt2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_7;
        DENDRO_igt1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_2;
        DENDRO_igt0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_0_gt5;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_2_gt2;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_0_gt4;
        double DENDRO_5  = grad_2_gt1;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_1_gt2;
        double DENDRO_8  = grad_0_gt3;
        double DENDRO_9  = grad_1_gt1;

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;
        DENDRO_3         = grad_2_gt0;
        DENDRO_8         = grad_1_gt0;
        DENDRO_9         = grad_0_gt0;

        // initialize reduction for
        double DENDRO_11 = 0;

        DENDRO_11 += DENDRO_0;
        DENDRO_11 += DENDRO_10;
        DENDRO_C1_k0_5 = DENDRO_11;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_2;
        DENDRO_C1_k0_4 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_7;
        DENDRO_C1_k0_3 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_3;
        DENDRO_C1_k0_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_8;
        DENDRO_C1_k0_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_9;
        DENDRO_C1_k0_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_1_gt5;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_2_gt4;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_1_gt2;
        double DENDRO_5  = grad_2_gt1;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_0_gt4;
        double DENDRO_8  = grad_1_gt0;
        double DENDRO_9  = grad_0_gt1;

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;
        DENDRO_2         = grad_2_gt3;
        double DENDRO_11 = grad_1_gt3;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;
        DENDRO_7         = grad_0_gt3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_0;
        DENDRO_3 += DENDRO_10;
        DENDRO_C1_k1_5 = DENDRO_3;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_2;
        DENDRO_C1_k1_4 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_11;
        DENDRO_C1_k1_3 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_12;
        DENDRO_C1_k1_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_7;
        DENDRO_C1_k1_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_13;
        DENDRO_C1_k1_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_2_gt3;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_1_gt4;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_2_gt1;
        double DENDRO_5  = grad_1_gt2;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_0_gt4;
        double DENDRO_8  = grad_2_gt0;
        double DENDRO_9  = grad_0_gt2;
        double DENDRO_10 = grad_2_gt5;
        double DENDRO_11 = grad_1_gt5;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;
        DENDRO_2         = grad_0_gt5;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_6;
        DENDRO_3 *= DENDRO_10;
        DENDRO_C1_k2_5 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 1;

        DENDRO_3 *= DENDRO_6;
        DENDRO_3 *= DENDRO_11;
        DENDRO_C1_k2_4 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_0;
        DENDRO_3 += DENDRO_12;
        DENDRO_C1_k2_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_2;
        DENDRO_C1_k2_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_13;
        DENDRO_C1_k2_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_7;
        DENDRO_C1_k2_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = DENDRO_igt2;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt1;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt0;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k0_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k0_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k0_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k0_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k0_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k0_0 = DENDRO_1;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt3;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt1;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k1_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k1_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k1_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k1_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k1_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k1_0 = DENDRO_1;
    }
    {
        double DENDRO_0  = DENDRO_igt5;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt4;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt2;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k2_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k2_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k2_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k2_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k2_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k2_0 = DENDRO_1;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt2;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt1;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt0;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1         = gt2[pp];
        DENDRO_3         = -1;
        DENDRO_6         = gt1[pp];
        double DENDRO_7  = gt0[pp];
        double DENDRO_8  = 2;
        double DENDRO_9  = chi[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_3;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_3;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;
        DENDRO_4         = DENDRO_9;

        DENDRO_4         = 1 / DENDRO_4;
        DENDRO_3         = gt5[pp];
        DENDRO_8         = 0.5;
        DENDRO_9         = gt4[pp];
        double DENDRO_11 = gt3[pp];

        // initialize reduction for
        double DENDRO_12 = 0;

        DENDRO_12 += DENDRO_0;
        DENDRO_12 += DENDRO_10;
        DENDRO_0  = -0.5;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_2;
        DENDRO_10 += DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_8;
        DENDRO_2 *= DENDRO_3;
        DENDRO_2 *= DENDRO_4;
        DENDRO_2 *= DENDRO_5;
        DENDRO_3 = DENDRO_C2_k0_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_9;
        DENDRO_6 *= DENDRO_4;
        DENDRO_6 *= DENDRO_5;
        DENDRO_7 = DENDRO_C2_k0_4;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_8;
        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_4;
        DENDRO_9 *= DENDRO_5;
        DENDRO_5 = DENDRO_C2_k0_3;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_0;
        DENDRO_8 *= DENDRO_4;
        DENDRO_8 *= DENDRO_12;
        DENDRO_11 = DENDRO_C2_k0_2;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_0;
        DENDRO_12 *= DENDRO_4;
        DENDRO_12 *= DENDRO_10;
        DENDRO_10        = DENDRO_C2_k0_1;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_0;
        DENDRO_13 *= DENDRO_4;
        DENDRO_13 *= DENDRO_1;
        DENDRO_0 = DENDRO_C2_k0_0;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_3;
        DENDRO_1 += DENDRO_2;
        DENDRO_C3_k0_5 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_C3_k0_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_5;
        DENDRO_1 += DENDRO_9;
        DENDRO_C3_k0_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_11;
        DENDRO_1 += DENDRO_8;
        DENDRO_C3_k0_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_12;
        DENDRO_C3_k0_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_13;
        DENDRO_C3_k0_0 = DENDRO_1;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt4;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt3;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt1;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1         = gt4[pp];
        DENDRO_3         = -1;
        DENDRO_6         = gt3[pp];
        double DENDRO_7  = 2;
        double DENDRO_8  = gt1[pp];
        double DENDRO_9  = chi[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_3;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_3;
        DENDRO_2 *= DENDRO_8;
        DENDRO_2 *= DENDRO_5;
        DENDRO_7 = DENDRO_9;

        DENDRO_7 = 1 / DENDRO_7;
        DENDRO_3 = gt5[pp];
        DENDRO_8 = 0.5;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_0;
        DENDRO_9 += DENDRO_10;
        DENDRO_0  = -0.5;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_6;
        DENDRO_10 += DENDRO_1;
        DENDRO_1 = gt2[pp];

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_4;
        DENDRO_6 += DENDRO_2;
        DENDRO_2 = gt0[pp];

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_8;
        DENDRO_4 *= DENDRO_3;
        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_5;
        DENDRO_3         = DENDRO_C2_k1_5;

        // initialize reduction for
        double DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_0;
        DENDRO_11 *= DENDRO_7;
        DENDRO_11 *= DENDRO_9;
        DENDRO_9         = DENDRO_C2_k1_4;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_0;
        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_10;
        DENDRO_10        = DENDRO_C2_k1_3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_8;
        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_7;
        DENDRO_13 *= DENDRO_5;
        DENDRO_1         = DENDRO_C2_k1_2;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_0;
        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_6;
        DENDRO_0 = DENDRO_C2_k1_1;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_2;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;
        DENDRO_2 = DENDRO_C2_k1_0;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_4;
        DENDRO_C3_k1_5 = DENDRO_5;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_9;
        DENDRO_3 += DENDRO_11;
        DENDRO_C3_k1_4 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_10;
        DENDRO_3 += DENDRO_12;
        DENDRO_C3_k1_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_1;
        DENDRO_3 += DENDRO_13;
        DENDRO_C3_k1_2 = DENDRO_3;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_14;
        DENDRO_C3_k1_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_2;
        DENDRO_0 += DENDRO_6;
        DENDRO_C3_k1_0 = DENDRO_0;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt5;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt4;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt2;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1        = gt5[pp];
        DENDRO_3        = -1;
        DENDRO_6        = 2;
        double DENDRO_7 = gt4[pp];
        double DENDRO_8 = gt2[pp];

        // initialize reduction for
        double DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_3;
        DENDRO_9 *= DENDRO_1;
        DENDRO_9 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_0;
        DENDRO_0 = chi[pp];

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_3;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_3;
        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_5;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_1;
        DENDRO_8 += DENDRO_9;
        DENDRO_1 = DENDRO_0;

        DENDRO_1 = 1 / DENDRO_1;
        DENDRO_0 = -0.5;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_2;
        DENDRO_3 += DENDRO_6;
        DENDRO_2 = gt3[pp];
        DENDRO_6 = 0.5;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_4;
        DENDRO_9 += DENDRO_7;
        DENDRO_4         = gt1[pp];
        DENDRO_7         = gt0[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_0;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_8;
        DENDRO_8         = DENDRO_C2_k2_5;

        // initialize reduction for
        double DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_0;
        DENDRO_11 *= DENDRO_1;
        DENDRO_11 *= DENDRO_3;
        DENDRO_3         = DENDRO_C2_k2_4;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_6;
        DENDRO_12 *= DENDRO_2;
        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_5;
        DENDRO_2         = DENDRO_C2_k2_3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_0;
        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_9;
        DENDRO_0 = DENDRO_C2_k2_2;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_6;
        DENDRO_9 *= DENDRO_4;
        DENDRO_9 *= DENDRO_1;
        DENDRO_9 *= DENDRO_5;
        DENDRO_4         = DENDRO_C2_k2_1;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_6;
        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_1;
        DENDRO_14 *= DENDRO_5;
        DENDRO_1 = DENDRO_C2_k2_0;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_8;
        DENDRO_5 += DENDRO_10;
        DENDRO_C3_k2_5 = DENDRO_5;

        // initialize reduction for
        DENDRO_5       = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_11;
        DENDRO_C3_k2_4 = DENDRO_5;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_2;
        DENDRO_3 += DENDRO_12;
        DENDRO_C3_k2_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_2       = 0;

        DENDRO_2 += DENDRO_0;
        DENDRO_2 += DENDRO_13;
        DENDRO_C3_k2_2 = DENDRO_2;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_9;
        DENDRO_C3_k2_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_14;
        DENDRO_C3_k2_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C2_k2_4;
        double DENDRO_2  = 2;
        double DENDRO_3  = DENDRO_igt2;
        double DENDRO_4  = DENDRO_C2_k2_2;
        double DENDRO_5  = DENDRO_igt1;
        double DENDRO_6  = DENDRO_C2_k2_1;
        double DENDRO_7  = DENDRO_igt5;
        double DENDRO_8  = DENDRO_C2_k2_5;
        double DENDRO_9  = DENDRO_igt3;
        double DENDRO_10 = DENDRO_C2_k2_3;
        double DENDRO_11 = DENDRO_igt0;
        double DENDRO_12 = DENDRO_C2_k2_0;
        double DENDRO_13 = DENDRO_C2_k1_4;
        double DENDRO_14 = DENDRO_C2_k1_2;
        double DENDRO_15 = DENDRO_C2_k1_1;
        double DENDRO_16 = DENDRO_C2_k1_5;
        double DENDRO_17 = DENDRO_C2_k1_3;
        double DENDRO_18 = DENDRO_C2_k1_0;
        double DENDRO_19 = DENDRO_C2_k0_4;
        double DENDRO_20 = DENDRO_C2_k0_2;
        double DENDRO_21 = DENDRO_C2_k0_1;
        double DENDRO_22 = DENDRO_C2_k0_5;
        double DENDRO_23 = DENDRO_C2_k0_3;
        double DENDRO_24 = DENDRO_C2_k0_0;

        // initialize reduction for
        double DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_1;
        DENDRO_25 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_2;
        DENDRO_1 *= DENDRO_4;
        DENDRO_1 *= DENDRO_3;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_2;
        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_7;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_10;
        DENDRO_8 *= DENDRO_9;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_12;
        DENDRO_10 *= DENDRO_11;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_2;
        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_0;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_2;
        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_3;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_2;
        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_5;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_7;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_9;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_11;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_2;
        DENDRO_18 *= DENDRO_19;
        DENDRO_18 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_20;
        DENDRO_0 *= DENDRO_3;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_2;
        DENDRO_3 *= DENDRO_21;
        DENDRO_3 *= DENDRO_5;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_22;
        DENDRO_2 *= DENDRO_7;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_23;
        DENDRO_5 *= DENDRO_9;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_24;
        DENDRO_7 *= DENDRO_11;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_10;
        DENDRO_9 += DENDRO_8;
        DENDRO_9 += DENDRO_6;
        DENDRO_9 += DENDRO_4;
        DENDRO_9 += DENDRO_1;
        DENDRO_9 += DENDRO_25;
        DENDRO_Gtk2 = DENDRO_9;

        // initialize reduction for
        DENDRO_1    = 0;

        DENDRO_1 += DENDRO_17;
        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_Gtk1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1    = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_5;
        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_3;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_18;
        DENDRO_Gtk0 = DENDRO_1;
    }
    double DENDRO_0  = 2;
    double DENDRO_1  = DENDRO_igt5;
    double DENDRO_2  = DENDRO_igt4;
    double DENDRO_3  = DENDRO_igt2;
    double DENDRO_4  = DENDRO_igt3;
    double DENDRO_5  = DENDRO_igt1;
    double DENDRO_6  = DENDRO_igt0;
    double DENDRO_7  = At4[pp];
    double DENDRO_8  = At2[pp];
    double DENDRO_9  = At1[pp];
    double DENDRO_10 = DENDRO_1;

    DENDRO_10 *= DENDRO_1;
    double DENDRO_11 = At5[pp];
    double DENDRO_12 = DENDRO_2;

    DENDRO_12 *= DENDRO_2;
    double DENDRO_13 = At3[pp];
    double DENDRO_14 = DENDRO_3;

    DENDRO_14 *= DENDRO_3;
    double DENDRO_15 = At0[pp];
    double DENDRO_16 = DENDRO_4;

    DENDRO_16 *= DENDRO_4;
    double DENDRO_17 = DENDRO_5;

    DENDRO_17 *= DENDRO_5;
    double DENDRO_18 = DENDRO_6;

    DENDRO_18 *= DENDRO_6;

    // initialize reduction for
    double DENDRO_19 = 1;

    DENDRO_19 *= DENDRO_0;
    DENDRO_19 *= DENDRO_7;
    DENDRO_19 *= DENDRO_2;
    DENDRO_19 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_20 = 1;

    DENDRO_20 *= DENDRO_0;
    DENDRO_20 *= DENDRO_8;
    DENDRO_20 *= DENDRO_3;
    DENDRO_20 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_21 = 1;

    DENDRO_21 *= DENDRO_0;
    DENDRO_21 *= DENDRO_9;
    DENDRO_21 *= DENDRO_3;
    DENDRO_21 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_22 = 1;

    DENDRO_22 *= DENDRO_11;
    DENDRO_22 *= DENDRO_10;

    // initialize reduction for
    double DENDRO_23 = 1;

    DENDRO_23 *= DENDRO_13;
    DENDRO_23 *= DENDRO_12;

    // initialize reduction for
    double DENDRO_24 = 1;

    DENDRO_24 *= DENDRO_15;
    DENDRO_24 *= DENDRO_14;
    double DENDRO_25 = -1;
    double DENDRO_26 = chi[pp];

    // initialize reduction for
    double DENDRO_27 = 1;

    DENDRO_27 *= DENDRO_11;
    DENDRO_27 *= DENDRO_2;
    DENDRO_27 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_28 = 1;

    DENDRO_28 *= DENDRO_7;
    DENDRO_28 *= DENDRO_4;
    DENDRO_28 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_29 = 1;

    DENDRO_29 *= DENDRO_13;
    DENDRO_29 *= DENDRO_4;
    DENDRO_29 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_30 = 1;

    DENDRO_30 *= DENDRO_8;
    DENDRO_30 *= DENDRO_3;
    DENDRO_30 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_31 = 1;

    DENDRO_31 *= DENDRO_8;
    DENDRO_31 *= DENDRO_5;
    DENDRO_31 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_32 = 1;

    DENDRO_32 *= DENDRO_9;
    DENDRO_32 *= DENDRO_3;
    DENDRO_32 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_9;
    DENDRO_33 *= DENDRO_5;
    DENDRO_33 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_34 = 1;

    DENDRO_34 *= DENDRO_15;
    DENDRO_34 *= DENDRO_5;
    DENDRO_34 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_35 = 1;

    DENDRO_35 *= DENDRO_7;
    DENDRO_35 *= DENDRO_12;

    // initialize reduction for
    double DENDRO_36 = 1;

    DENDRO_36 *= DENDRO_11;
    DENDRO_36 *= DENDRO_3;
    DENDRO_36 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_7;
    DENDRO_37 *= DENDRO_3;
    DENDRO_37 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_38 = 1;

    DENDRO_38 *= DENDRO_7;
    DENDRO_38 *= DENDRO_5;
    DENDRO_38 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_13;
    DENDRO_39 *= DENDRO_5;
    DENDRO_39 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_8;
    DENDRO_40 *= DENDRO_6;
    DENDRO_40 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_9;
    DENDRO_41 *= DENDRO_5;
    DENDRO_41 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_9;
    DENDRO_42 *= DENDRO_6;
    DENDRO_42 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_43 = 1;

    DENDRO_43 *= DENDRO_15;
    DENDRO_43 *= DENDRO_6;
    DENDRO_43 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_44 = 1;

    DENDRO_44 *= DENDRO_8;
    DENDRO_44 *= DENDRO_14;

    // initialize reduction for
    double DENDRO_45 = 1;

    DENDRO_45 *= DENDRO_0;
    DENDRO_45 *= DENDRO_7;
    DENDRO_45 *= DENDRO_4;
    DENDRO_45 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_46 = 1;

    DENDRO_46 *= DENDRO_0;
    DENDRO_46 *= DENDRO_8;
    DENDRO_46 *= DENDRO_5;
    DENDRO_46 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_47 = 1;

    DENDRO_47 *= DENDRO_0;
    DENDRO_47 *= DENDRO_9;
    DENDRO_47 *= DENDRO_5;
    DENDRO_47 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_48 = 1;

    DENDRO_48 *= DENDRO_11;
    DENDRO_48 *= DENDRO_12;

    // initialize reduction for
    double DENDRO_49 = 1;

    DENDRO_49 *= DENDRO_13;
    DENDRO_49 *= DENDRO_16;

    // initialize reduction for
    double DENDRO_50 = 1;

    DENDRO_50 *= DENDRO_15;
    DENDRO_50 *= DENDRO_17;

    // initialize reduction for
    double DENDRO_51 = 1;

    DENDRO_51 *= DENDRO_11;
    DENDRO_51 *= DENDRO_3;
    DENDRO_51 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_7;
    DENDRO_52 *= DENDRO_3;
    DENDRO_52 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_7;
    DENDRO_53 *= DENDRO_5;
    DENDRO_53 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_54 = 1;

    DENDRO_54 *= DENDRO_13;
    DENDRO_54 *= DENDRO_5;
    DENDRO_54 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_8;
    DENDRO_55 *= DENDRO_5;
    DENDRO_55 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_8;
    DENDRO_56 *= DENDRO_6;
    DENDRO_56 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_57 = 1;

    DENDRO_57 *= DENDRO_9;
    DENDRO_57 *= DENDRO_6;
    DENDRO_57 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_58 = 1;

    DENDRO_58 *= DENDRO_15;
    DENDRO_58 *= DENDRO_6;
    DENDRO_58 *= DENDRO_5;

    // initialize reduction for
    double DENDRO_59 = 1;

    DENDRO_59 *= DENDRO_9;
    DENDRO_59 *= DENDRO_17;

    // initialize reduction for
    double DENDRO_60 = 1;

    DENDRO_60 *= DENDRO_0;
    DENDRO_60 *= DENDRO_7;
    DENDRO_60 *= DENDRO_5;
    DENDRO_60 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_0;
    DENDRO_61 *= DENDRO_8;
    DENDRO_61 *= DENDRO_6;
    DENDRO_61 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_0;
    DENDRO_62 *= DENDRO_9;
    DENDRO_62 *= DENDRO_6;
    DENDRO_62 *= DENDRO_5;

    // initialize reduction for
    double DENDRO_63 = 1;

    DENDRO_63 *= DENDRO_11;
    DENDRO_63 *= DENDRO_14;

    // initialize reduction for
    double DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_13;
    DENDRO_64 *= DENDRO_17;

    // initialize reduction for
    double DENDRO_65 = 1;

    DENDRO_65 *= DENDRO_15;
    DENDRO_65 *= DENDRO_18;
    double DENDRO_66 = grad_2_Gt2;
    double DENDRO_67 = beta2[pp];
    double DENDRO_68 = grad_1_Gt2;
    double DENDRO_69 = beta1[pp];
    double DENDRO_70 = grad_0_Gt2;
    double DENDRO_71 = beta0[pp];
    double DENDRO_72 = 4;

    // initialize reduction for
    double DENDRO_73 = 0;

    DENDRO_73 += DENDRO_24;
    DENDRO_73 += DENDRO_23;
    DENDRO_73 += DENDRO_22;
    DENDRO_73 += DENDRO_21;
    DENDRO_73 += DENDRO_20;
    DENDRO_73 += DENDRO_19;
    DENDRO_19 = DENDRO_26;

    DENDRO_19 = 1 / DENDRO_19;
    DENDRO_20 = grad_2_chi;
    DENDRO_21 = 3;
    DENDRO_22 = grad_2_K;
    DENDRO_23 = 4.0 / 3.0;

    // initialize reduction for
    DENDRO_24 = 0;

    DENDRO_24 += DENDRO_35;
    DENDRO_24 += DENDRO_34;
    DENDRO_24 += DENDRO_33;
    DENDRO_24 += DENDRO_32;
    DENDRO_24 += DENDRO_31;
    DENDRO_24 += DENDRO_30;
    DENDRO_24 += DENDRO_29;
    DENDRO_24 += DENDRO_28;
    DENDRO_24 += DENDRO_27;
    DENDRO_26 = grad_1_chi;
    DENDRO_27 = grad_1_K;

    // initialize reduction for
    DENDRO_28 = 0;

    DENDRO_28 += DENDRO_44;
    DENDRO_28 += DENDRO_43;
    DENDRO_28 += DENDRO_42;
    DENDRO_28 += DENDRO_41;
    DENDRO_28 += DENDRO_40;
    DENDRO_28 += DENDRO_39;
    DENDRO_28 += DENDRO_38;
    DENDRO_28 += DENDRO_37;
    DENDRO_28 += DENDRO_36;
    DENDRO_29 = grad_0_chi;
    DENDRO_30 = grad_0_K;
    DENDRO_31 = grad_2_B2;
    DENDRO_32 = grad_1_B2;
    DENDRO_33 = grad_0_B2;
    DENDRO_34 = grad_2_Gt1;
    DENDRO_35 = grad_1_Gt1;
    DENDRO_36 = grad_0_Gt1;

    // initialize reduction for
    DENDRO_37 = 0;

    DENDRO_37 += DENDRO_50;
    DENDRO_37 += DENDRO_49;
    DENDRO_37 += DENDRO_48;
    DENDRO_37 += DENDRO_47;
    DENDRO_37 += DENDRO_46;
    DENDRO_37 += DENDRO_45;

    // initialize reduction for
    DENDRO_38 = 0;

    DENDRO_38 += DENDRO_59;
    DENDRO_38 += DENDRO_58;
    DENDRO_38 += DENDRO_57;
    DENDRO_38 += DENDRO_56;
    DENDRO_38 += DENDRO_55;
    DENDRO_38 += DENDRO_54;
    DENDRO_38 += DENDRO_53;
    DENDRO_38 += DENDRO_52;
    DENDRO_38 += DENDRO_51;
    DENDRO_39 = grad_2_B1;
    DENDRO_40 = grad_1_B1;
    DENDRO_41 = grad_0_B1;
    DENDRO_42 = grad_2_Gt0;
    DENDRO_43 = grad_1_Gt0;
    DENDRO_44 = grad_0_Gt0;

    // initialize reduction for
    DENDRO_45 = 0;

    DENDRO_45 += DENDRO_65;
    DENDRO_45 += DENDRO_64;
    DENDRO_45 += DENDRO_63;
    DENDRO_45 += DENDRO_62;
    DENDRO_45 += DENDRO_61;
    DENDRO_45 += DENDRO_60;
    DENDRO_46 = grad_2_B0;
    DENDRO_47 = grad_1_B0;
    DENDRO_48 = grad_0_B0;
    DENDRO_49 = grad_2_beta2;
    DENDRO_50 = grad_1_beta1;
    DENDRO_51 = grad_0_beta0;

    // initialize reduction for
    DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_67;
    DENDRO_52 *= DENDRO_66;

    // initialize reduction for
    DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_69;
    DENDRO_53 *= DENDRO_68;

    // initialize reduction for
    DENDRO_54 = 1;

    DENDRO_54 *= DENDRO_71;
    DENDRO_54 *= DENDRO_70;

    // initialize reduction for
    DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_72;
    DENDRO_55 *= DENDRO_7;
    DENDRO_55 *= DENDRO_2;
    DENDRO_55 *= DENDRO_1;

    // initialize reduction for
    DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_72;
    DENDRO_56 *= DENDRO_8;
    DENDRO_56 *= DENDRO_3;
    DENDRO_56 *= DENDRO_1;

    // initialize reduction for
    DENDRO_57 = 1;

    DENDRO_57 *= DENDRO_72;
    DENDRO_57 *= DENDRO_9;
    DENDRO_57 *= DENDRO_3;
    DENDRO_57 *= DENDRO_2;

    // initialize reduction for
    DENDRO_58 = 1;

    DENDRO_58 *= DENDRO_0;
    DENDRO_58 *= DENDRO_11;
    DENDRO_58 *= DENDRO_10;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_0;
    DENDRO_10 *= DENDRO_13;
    DENDRO_10 *= DENDRO_12;

    // initialize reduction for
    DENDRO_59 = 1;

    DENDRO_59 *= DENDRO_0;
    DENDRO_59 *= DENDRO_15;
    DENDRO_59 *= DENDRO_14;

    // initialize reduction for
    DENDRO_60 = 1;

    DENDRO_60 *= DENDRO_0;
    DENDRO_60 *= DENDRO_11;
    DENDRO_60 *= DENDRO_2;
    DENDRO_60 *= DENDRO_1;

    // initialize reduction for
    DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_0;
    DENDRO_61 *= DENDRO_7;
    DENDRO_61 *= DENDRO_4;
    DENDRO_61 *= DENDRO_1;

    // initialize reduction for
    DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_0;
    DENDRO_62 *= DENDRO_13;
    DENDRO_62 *= DENDRO_4;
    DENDRO_62 *= DENDRO_2;

    // initialize reduction for
    DENDRO_63 = 1;

    DENDRO_63 *= DENDRO_0;
    DENDRO_63 *= DENDRO_8;
    DENDRO_63 *= DENDRO_3;
    DENDRO_63 *= DENDRO_2;

    // initialize reduction for
    DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_0;
    DENDRO_64 *= DENDRO_8;
    DENDRO_64 *= DENDRO_5;
    DENDRO_64 *= DENDRO_1;

    // initialize reduction for
    DENDRO_65 = 1;

    DENDRO_65 *= DENDRO_0;
    DENDRO_65 *= DENDRO_9;
    DENDRO_65 *= DENDRO_3;
    DENDRO_65 *= DENDRO_4;

    // initialize reduction for
    DENDRO_66 = 1;

    DENDRO_66 *= DENDRO_0;
    DENDRO_66 *= DENDRO_9;
    DENDRO_66 *= DENDRO_5;
    DENDRO_66 *= DENDRO_2;

    // initialize reduction for
    DENDRO_68 = 1;

    DENDRO_68 *= DENDRO_0;
    DENDRO_68 *= DENDRO_15;
    DENDRO_68 *= DENDRO_5;
    DENDRO_68 *= DENDRO_3;

    // initialize reduction for
    DENDRO_70 = 1;

    DENDRO_70 *= DENDRO_0;
    DENDRO_70 *= DENDRO_7;
    DENDRO_70 *= DENDRO_12;

    // initialize reduction for
    double DENDRO_74 = 1;

    DENDRO_74 *= DENDRO_0;
    DENDRO_74 *= DENDRO_11;
    DENDRO_74 *= DENDRO_3;
    DENDRO_74 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_75 = 1;

    DENDRO_75 *= DENDRO_0;
    DENDRO_75 *= DENDRO_7;
    DENDRO_75 *= DENDRO_3;
    DENDRO_75 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_76 = 1;

    DENDRO_76 *= DENDRO_0;
    DENDRO_76 *= DENDRO_7;
    DENDRO_76 *= DENDRO_5;
    DENDRO_76 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_77 = 1;

    DENDRO_77 *= DENDRO_0;
    DENDRO_77 *= DENDRO_13;
    DENDRO_77 *= DENDRO_5;
    DENDRO_77 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_78 = 1;

    DENDRO_78 *= DENDRO_0;
    DENDRO_78 *= DENDRO_8;
    DENDRO_78 *= DENDRO_6;
    DENDRO_78 *= DENDRO_1;

    // initialize reduction for
    double DENDRO_79 = 1;

    DENDRO_79 *= DENDRO_0;
    DENDRO_79 *= DENDRO_9;
    DENDRO_79 *= DENDRO_5;
    DENDRO_79 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_80 = 1;

    DENDRO_80 *= DENDRO_0;
    DENDRO_80 *= DENDRO_9;
    DENDRO_80 *= DENDRO_6;
    DENDRO_80 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_81 = 1;

    DENDRO_81 *= DENDRO_0;
    DENDRO_81 *= DENDRO_15;
    DENDRO_81 *= DENDRO_6;
    DENDRO_81 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_82 = 1;

    DENDRO_82 *= DENDRO_0;
    DENDRO_82 *= DENDRO_8;
    DENDRO_82 *= DENDRO_14;

    // initialize reduction for
    double DENDRO_83 = 1;

    DENDRO_83 *= DENDRO_21;
    DENDRO_83 *= DENDRO_20;
    DENDRO_83 *= DENDRO_19;
    DENDRO_83 *= DENDRO_73;

    // initialize reduction for
    double DENDRO_84 = 1;

    DENDRO_84 *= DENDRO_23;
    DENDRO_84 *= DENDRO_1;
    DENDRO_84 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_85 = 1;

    DENDRO_85 *= DENDRO_21;
    DENDRO_85 *= DENDRO_26;
    DENDRO_85 *= DENDRO_19;
    DENDRO_85 *= DENDRO_24;

    // initialize reduction for
    double DENDRO_86 = 1;

    DENDRO_86 *= DENDRO_23;
    DENDRO_86 *= DENDRO_2;
    DENDRO_86 *= DENDRO_27;

    // initialize reduction for
    double DENDRO_87 = 1;

    DENDRO_87 *= DENDRO_21;
    DENDRO_87 *= DENDRO_29;
    DENDRO_87 *= DENDRO_19;
    DENDRO_87 *= DENDRO_28;

    // initialize reduction for
    double DENDRO_88 = 1;

    DENDRO_88 *= DENDRO_23;
    DENDRO_88 *= DENDRO_3;
    DENDRO_88 *= DENDRO_30;

    // initialize reduction for
    double DENDRO_89 = 1;

    DENDRO_89 *= DENDRO_67;
    DENDRO_89 *= DENDRO_31;

    // initialize reduction for
    DENDRO_31 = 1;

    DENDRO_31 *= DENDRO_69;
    DENDRO_31 *= DENDRO_32;

    // initialize reduction for
    DENDRO_32 = 1;

    DENDRO_32 *= DENDRO_71;
    DENDRO_32 *= DENDRO_33;

    // initialize reduction for
    DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_67;
    DENDRO_33 *= DENDRO_34;

    // initialize reduction for
    DENDRO_34 = 1;

    DENDRO_34 *= DENDRO_69;
    DENDRO_34 *= DENDRO_35;

    // initialize reduction for
    DENDRO_35 = 1;

    DENDRO_35 *= DENDRO_71;
    DENDRO_35 *= DENDRO_36;

    // initialize reduction for
    DENDRO_36 = 1;

    DENDRO_36 *= DENDRO_72;
    DENDRO_36 *= DENDRO_7;
    DENDRO_36 *= DENDRO_4;
    DENDRO_36 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_90 = 1;

    DENDRO_90 *= DENDRO_72;
    DENDRO_90 *= DENDRO_8;
    DENDRO_90 *= DENDRO_5;
    DENDRO_90 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_91 = 1;

    DENDRO_91 *= DENDRO_72;
    DENDRO_91 *= DENDRO_9;
    DENDRO_91 *= DENDRO_5;
    DENDRO_91 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_92 = 1;

    DENDRO_92 *= DENDRO_0;
    DENDRO_92 *= DENDRO_11;
    DENDRO_92 *= DENDRO_12;

    // initialize reduction for
    DENDRO_12 = 1;

    DENDRO_12 *= DENDRO_0;
    DENDRO_12 *= DENDRO_13;
    DENDRO_12 *= DENDRO_16;

    // initialize reduction for
    DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_0;
    DENDRO_16 *= DENDRO_15;
    DENDRO_16 *= DENDRO_17;

    // initialize reduction for
    double DENDRO_93 = 1;

    DENDRO_93 *= DENDRO_0;
    DENDRO_93 *= DENDRO_11;
    DENDRO_93 *= DENDRO_3;
    DENDRO_93 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_94 = 1;

    DENDRO_94 *= DENDRO_0;
    DENDRO_94 *= DENDRO_7;
    DENDRO_94 *= DENDRO_3;
    DENDRO_94 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_95 = 1;

    DENDRO_95 *= DENDRO_0;
    DENDRO_95 *= DENDRO_7;
    DENDRO_95 *= DENDRO_5;
    DENDRO_95 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_96 = 1;

    DENDRO_96 *= DENDRO_0;
    DENDRO_96 *= DENDRO_13;
    DENDRO_96 *= DENDRO_5;
    DENDRO_96 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_97 = 1;

    DENDRO_97 *= DENDRO_0;
    DENDRO_97 *= DENDRO_8;
    DENDRO_97 *= DENDRO_5;
    DENDRO_97 *= DENDRO_3;

    // initialize reduction for
    double DENDRO_98 = 1;

    DENDRO_98 *= DENDRO_0;
    DENDRO_98 *= DENDRO_8;
    DENDRO_98 *= DENDRO_6;
    DENDRO_98 *= DENDRO_2;

    // initialize reduction for
    double DENDRO_99 = 1;

    DENDRO_99 *= DENDRO_0;
    DENDRO_99 *= DENDRO_9;
    DENDRO_99 *= DENDRO_6;
    DENDRO_99 *= DENDRO_4;

    // initialize reduction for
    double DENDRO_100 = 1;

    DENDRO_100 *= DENDRO_0;
    DENDRO_100 *= DENDRO_15;
    DENDRO_100 *= DENDRO_6;
    DENDRO_100 *= DENDRO_5;

    // initialize reduction for
    double DENDRO_101 = 1;

    DENDRO_101 *= DENDRO_0;
    DENDRO_101 *= DENDRO_9;
    DENDRO_101 *= DENDRO_17;

    // initialize reduction for
    double DENDRO_102 = 1;

    DENDRO_102 *= DENDRO_21;
    DENDRO_102 *= DENDRO_20;
    DENDRO_102 *= DENDRO_19;
    DENDRO_102 *= DENDRO_24;

    // initialize reduction for
    double DENDRO_103 = 1;

    DENDRO_103 *= DENDRO_23;
    DENDRO_103 *= DENDRO_2;
    DENDRO_103 *= DENDRO_22;

    // initialize reduction for
    double DENDRO_104 = 1;

    DENDRO_104 *= DENDRO_21;
    DENDRO_104 *= DENDRO_26;
    DENDRO_104 *= DENDRO_19;
    DENDRO_104 *= DENDRO_37;

    // initialize reduction for
    double DENDRO_105 = 1;

    DENDRO_105 *= DENDRO_23;
    DENDRO_105 *= DENDRO_4;
    DENDRO_105 *= DENDRO_27;

    // initialize reduction for
    double DENDRO_106 = 1;

    DENDRO_106 *= DENDRO_21;
    DENDRO_106 *= DENDRO_29;
    DENDRO_106 *= DENDRO_19;
    DENDRO_106 *= DENDRO_38;

    // initialize reduction for
    double DENDRO_107 = 1;

    DENDRO_107 *= DENDRO_23;
    DENDRO_107 *= DENDRO_5;
    DENDRO_107 *= DENDRO_30;

    // initialize reduction for
    double DENDRO_108 = 1;

    DENDRO_108 *= DENDRO_67;
    DENDRO_108 *= DENDRO_39;

    // initialize reduction for
    DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_69;
    DENDRO_39 *= DENDRO_40;

    // initialize reduction for
    DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_71;
    DENDRO_40 *= DENDRO_41;

    // initialize reduction for
    DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_67;
    DENDRO_41 *= DENDRO_42;

    // initialize reduction for
    DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_69;
    DENDRO_42 *= DENDRO_43;

    // initialize reduction for
    DENDRO_43 = 1;

    DENDRO_43 *= DENDRO_71;
    DENDRO_43 *= DENDRO_44;

    // initialize reduction for
    DENDRO_44 = 1;

    DENDRO_44 *= DENDRO_72;
    DENDRO_44 *= DENDRO_7;
    DENDRO_44 *= DENDRO_5;
    DENDRO_44 *= DENDRO_3;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_72;
    DENDRO_7 *= DENDRO_8;
    DENDRO_7 *= DENDRO_6;
    DENDRO_7 *= DENDRO_3;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_72;
    DENDRO_8 *= DENDRO_9;
    DENDRO_8 *= DENDRO_6;
    DENDRO_8 *= DENDRO_5;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_0;
    DENDRO_9 *= DENDRO_11;
    DENDRO_9 *= DENDRO_14;

    // initialize reduction for
    DENDRO_11 = 1;

    DENDRO_11 *= DENDRO_0;
    DENDRO_11 *= DENDRO_13;
    DENDRO_11 *= DENDRO_17;

    // initialize reduction for
    DENDRO_13 = 1;

    DENDRO_13 *= DENDRO_0;
    DENDRO_13 *= DENDRO_15;
    DENDRO_13 *= DENDRO_18;

    // initialize reduction for
    DENDRO_14 = 1;

    DENDRO_14 *= DENDRO_21;
    DENDRO_14 *= DENDRO_20;
    DENDRO_14 *= DENDRO_19;
    DENDRO_14 *= DENDRO_28;

    // initialize reduction for
    DENDRO_15 = 1;

    DENDRO_15 *= DENDRO_23;
    DENDRO_15 *= DENDRO_3;
    DENDRO_15 *= DENDRO_22;

    // initialize reduction for
    DENDRO_17 = 1;

    DENDRO_17 *= DENDRO_21;
    DENDRO_17 *= DENDRO_26;
    DENDRO_17 *= DENDRO_19;
    DENDRO_17 *= DENDRO_38;

    // initialize reduction for
    DENDRO_18 = 1;

    DENDRO_18 *= DENDRO_23;
    DENDRO_18 *= DENDRO_5;
    DENDRO_18 *= DENDRO_27;

    // initialize reduction for
    DENDRO_20 = 1;

    DENDRO_20 *= DENDRO_21;
    DENDRO_20 *= DENDRO_29;
    DENDRO_20 *= DENDRO_19;
    DENDRO_20 *= DENDRO_45;

    // initialize reduction for
    DENDRO_19 = 1;

    DENDRO_19 *= DENDRO_23;
    DENDRO_19 *= DENDRO_6;
    DENDRO_19 *= DENDRO_30;

    // initialize reduction for
    DENDRO_21 = 1;

    DENDRO_21 *= DENDRO_67;
    DENDRO_21 *= DENDRO_46;

    // initialize reduction for
    DENDRO_22 = 1;

    DENDRO_22 *= DENDRO_69;
    DENDRO_22 *= DENDRO_47;

    // initialize reduction for
    DENDRO_26 = 1;

    DENDRO_26 *= DENDRO_71;
    DENDRO_26 *= DENDRO_48;
    DENDRO_27         = alpha[pp];
    DENDRO_29         = DENDRO_C2_k2_4;
    DENDRO_30         = DENDRO_C2_k2_2;
    DENDRO_46         = DENDRO_C2_k2_1;
    DENDRO_47         = DENDRO_C2_k2_5;
    DENDRO_48         = DENDRO_C2_k2_3;
    DENDRO_67         = DENDRO_C2_k2_0;
    DENDRO_69         = grad2_1_2_beta2;
    DENDRO_71         = 7.0 / 3.0;
    double DENDRO_109 = grad2_0_2_beta2;
    double DENDRO_110 = grad2_2_2_beta2;

    // initialize reduction for
    double DENDRO_111 = 0;

    DENDRO_111 += DENDRO_51;
    DENDRO_111 += DENDRO_50;
    DENDRO_111 += DENDRO_49;
    double DENDRO_112 = DENDRO_Gtk2;
    double DENDRO_113 = 2.0 / 3.0;
    double DENDRO_114 = grad2_1_2_beta1;
    double DENDRO_115 = 1.0 / 3.0;
    double DENDRO_116 = grad2_0_2_beta0;
    double DENDRO_117 = grad2_1_1_beta1;
    double DENDRO_118 = grad2_0_1_beta0;
    double DENDRO_119 = grad2_0_1_beta1;
    double DENDRO_120 = grad2_0_0_beta0;
    double DENDRO_121 = grad2_0_1_beta2;

    // initialize reduction for
    double DENDRO_122 = 0;

    DENDRO_122 += DENDRO_54;
    DENDRO_122 += DENDRO_53;
    DENDRO_122 += DENDRO_52;
    double DENDRO_123 = lambda[3];

    // initialize reduction for
    double DENDRO_124 = 0;

    DENDRO_124 += DENDRO_59;
    DENDRO_124 += DENDRO_10;
    DENDRO_124 += DENDRO_58;
    DENDRO_124 += DENDRO_57;
    DENDRO_124 += DENDRO_56;
    DENDRO_124 += DENDRO_55;
    DENDRO_10 = grad_2_alpha;

    // initialize reduction for
    DENDRO_55 = 0;

    DENDRO_55 += DENDRO_70;
    DENDRO_55 += DENDRO_68;
    DENDRO_55 += DENDRO_66;
    DENDRO_55 += DENDRO_65;
    DENDRO_55 += DENDRO_64;
    DENDRO_55 += DENDRO_63;
    DENDRO_55 += DENDRO_62;
    DENDRO_55 += DENDRO_61;
    DENDRO_55 += DENDRO_60;
    DENDRO_56 = grad_1_alpha;

    // initialize reduction for
    DENDRO_57 = 0;

    DENDRO_57 += DENDRO_82;
    DENDRO_57 += DENDRO_81;
    DENDRO_57 += DENDRO_80;
    DENDRO_57 += DENDRO_79;
    DENDRO_57 += DENDRO_78;
    DENDRO_57 += DENDRO_77;
    DENDRO_57 += DENDRO_76;
    DENDRO_57 += DENDRO_75;
    DENDRO_57 += DENDRO_74;
    DENDRO_58 = grad_0_alpha;

    // initialize reduction for
    DENDRO_59 = 0;

    DENDRO_59 += DENDRO_84;
    DENDRO_59 += DENDRO_83;

    // initialize reduction for
    DENDRO_60 = 0;

    DENDRO_60 += DENDRO_86;
    DENDRO_60 += DENDRO_85;

    // initialize reduction for
    DENDRO_61 = 0;

    DENDRO_61 += DENDRO_88;
    DENDRO_61 += DENDRO_87;
    DENDRO_62 = grad_1_beta2;
    DENDRO_63 = DENDRO_Gtk1;
    DENDRO_64 = grad_0_beta2;
    DENDRO_65 = DENDRO_Gtk0;
    DENDRO_66 = eta;
    DENDRO_68 = B2[pp];

    // initialize reduction for
    DENDRO_70 = 0;

    DENDRO_70 += DENDRO_32;
    DENDRO_70 += DENDRO_31;
    DENDRO_70 += DENDRO_89;
    DENDRO_31 = lambda[2];
    DENDRO_32 = grad2_1_1_beta2;
    DENDRO_74 = grad2_0_0_beta2;
    DENDRO_75 = DENDRO_C2_k1_4;
    DENDRO_76 = DENDRO_C2_k1_2;
    DENDRO_77 = DENDRO_C2_k1_1;
    DENDRO_78 = DENDRO_C2_k1_5;
    DENDRO_79 = DENDRO_C2_k1_3;
    DENDRO_80 = DENDRO_C2_k1_0;
    DENDRO_81 = grad2_0_2_beta1;

    // initialize reduction for
    DENDRO_82 = 0;

    DENDRO_82 += DENDRO_35;
    DENDRO_82 += DENDRO_34;
    DENDRO_82 += DENDRO_33;

    // initialize reduction for
    DENDRO_83 = 0;

    DENDRO_83 += DENDRO_16;
    DENDRO_83 += DENDRO_12;
    DENDRO_83 += DENDRO_92;
    DENDRO_83 += DENDRO_91;
    DENDRO_83 += DENDRO_90;
    DENDRO_83 += DENDRO_36;

    // initialize reduction for
    DENDRO_12 = 0;

    DENDRO_12 += DENDRO_101;
    DENDRO_12 += DENDRO_100;
    DENDRO_12 += DENDRO_99;
    DENDRO_12 += DENDRO_98;
    DENDRO_12 += DENDRO_97;
    DENDRO_12 += DENDRO_96;
    DENDRO_12 += DENDRO_95;
    DENDRO_12 += DENDRO_94;
    DENDRO_12 += DENDRO_93;

    // initialize reduction for
    DENDRO_16 = 0;

    DENDRO_16 += DENDRO_103;
    DENDRO_16 += DENDRO_102;

    // initialize reduction for
    DENDRO_36 = 0;

    DENDRO_36 += DENDRO_105;
    DENDRO_36 += DENDRO_104;

    // initialize reduction for
    DENDRO_84 = 0;

    DENDRO_84 += DENDRO_107;
    DENDRO_84 += DENDRO_106;
    DENDRO_85 = grad_2_beta1;
    DENDRO_86 = grad_0_beta1;
    DENDRO_87 = B1[pp];

    // initialize reduction for
    DENDRO_88 = 0;

    DENDRO_88 += DENDRO_40;
    DENDRO_88 += DENDRO_39;
    DENDRO_88 += DENDRO_108;
    DENDRO_39 = grad2_2_2_beta1;
    DENDRO_40 = grad2_0_0_beta1;
    DENDRO_89 = DENDRO_C2_k0_4;
    DENDRO_90 = DENDRO_C2_k0_2;
    DENDRO_91 = DENDRO_C2_k0_1;
    DENDRO_92 = DENDRO_C2_k0_5;
    DENDRO_93 = DENDRO_C2_k0_3;
    DENDRO_94 = DENDRO_C2_k0_0;
    DENDRO_95 = grad2_1_2_beta0;

    // initialize reduction for
    DENDRO_96 = 0;

    DENDRO_96 += DENDRO_43;
    DENDRO_96 += DENDRO_42;
    DENDRO_96 += DENDRO_41;

    // initialize reduction for
    DENDRO_97 = 0;

    DENDRO_97 += DENDRO_13;
    DENDRO_97 += DENDRO_11;
    DENDRO_97 += DENDRO_9;
    DENDRO_97 += DENDRO_8;
    DENDRO_97 += DENDRO_7;
    DENDRO_97 += DENDRO_44;

    // initialize reduction for
    DENDRO_7 = 0;

    DENDRO_7 += DENDRO_15;
    DENDRO_7 += DENDRO_14;

    // initialize reduction for
    DENDRO_8 = 0;

    DENDRO_8 += DENDRO_18;
    DENDRO_8 += DENDRO_17;

    // initialize reduction for
    DENDRO_9 = 0;

    DENDRO_9 += DENDRO_19;
    DENDRO_9 += DENDRO_20;
    DENDRO_11 = grad_2_beta0;
    DENDRO_13 = grad_1_beta0;
    DENDRO_14 = B0[pp];

    // initialize reduction for
    DENDRO_15 = 0;

    DENDRO_15 += DENDRO_26;
    DENDRO_15 += DENDRO_22;
    DENDRO_15 += DENDRO_21;
    DENDRO_17 = grad2_2_2_beta0;
    DENDRO_18 = grad2_1_1_beta0;

    // initialize reduction for
    DENDRO_19 = 1;

    DENDRO_19 *= DENDRO_72;
    DENDRO_19 *= DENDRO_29;
    DENDRO_19 *= DENDRO_27;
    DENDRO_19 *= DENDRO_24;

    // initialize reduction for
    DENDRO_20 = 1;

    DENDRO_20 *= DENDRO_72;
    DENDRO_20 *= DENDRO_30;
    DENDRO_20 *= DENDRO_27;
    DENDRO_20 *= DENDRO_28;

    // initialize reduction for
    DENDRO_21 = 1;

    DENDRO_21 *= DENDRO_72;
    DENDRO_21 *= DENDRO_46;
    DENDRO_21 *= DENDRO_27;
    DENDRO_21 *= DENDRO_38;

    // initialize reduction for
    DENDRO_22 = 1;

    DENDRO_22 *= DENDRO_0;
    DENDRO_22 *= DENDRO_47;
    DENDRO_22 *= DENDRO_27;
    DENDRO_22 *= DENDRO_73;

    // initialize reduction for
    DENDRO_26 = 1;

    DENDRO_26 *= DENDRO_0;
    DENDRO_26 *= DENDRO_48;
    DENDRO_26 *= DENDRO_27;
    DENDRO_26 *= DENDRO_37;

    // initialize reduction for
    DENDRO_29 = 1;

    DENDRO_29 *= DENDRO_0;
    DENDRO_29 *= DENDRO_67;
    DENDRO_29 *= DENDRO_27;
    DENDRO_29 *= DENDRO_45;

    // initialize reduction for
    DENDRO_30 = 1;

    DENDRO_30 *= DENDRO_71;
    DENDRO_30 *= DENDRO_2;
    DENDRO_30 *= DENDRO_69;

    // initialize reduction for
    DENDRO_44 = 1;

    DENDRO_44 *= DENDRO_71;
    DENDRO_44 *= DENDRO_3;
    DENDRO_44 *= DENDRO_109;

    // initialize reduction for
    DENDRO_46 = 1;

    DENDRO_46 *= DENDRO_23;
    DENDRO_46 *= DENDRO_1;
    DENDRO_46 *= DENDRO_110;

    // initialize reduction for
    DENDRO_47 = 1;

    DENDRO_47 *= DENDRO_113;
    DENDRO_47 *= DENDRO_112;
    DENDRO_47 *= DENDRO_111;

    // initialize reduction for
    DENDRO_48 = 1;

    DENDRO_48 *= DENDRO_115;
    DENDRO_48 *= DENDRO_1;
    DENDRO_48 *= DENDRO_114;

    // initialize reduction for
    DENDRO_67 = 1;

    DENDRO_67 *= DENDRO_115;
    DENDRO_67 *= DENDRO_1;
    DENDRO_67 *= DENDRO_116;

    // initialize reduction for
    DENDRO_98 = 1;

    DENDRO_98 *= DENDRO_115;
    DENDRO_98 *= DENDRO_2;
    DENDRO_98 *= DENDRO_117;

    // initialize reduction for
    DENDRO_99 = 1;

    DENDRO_99 *= DENDRO_115;
    DENDRO_99 *= DENDRO_2;
    DENDRO_99 *= DENDRO_118;

    // initialize reduction for
    DENDRO_100 = 1;

    DENDRO_100 *= DENDRO_115;
    DENDRO_100 *= DENDRO_3;
    DENDRO_100 *= DENDRO_119;

    // initialize reduction for
    DENDRO_101 = 1;

    DENDRO_101 *= DENDRO_115;
    DENDRO_101 *= DENDRO_3;
    DENDRO_101 *= DENDRO_120;

    // initialize reduction for
    DENDRO_102 = 1;

    DENDRO_102 *= DENDRO_0;
    DENDRO_102 *= DENDRO_5;
    DENDRO_102 *= DENDRO_121;

    // initialize reduction for
    DENDRO_103 = 1;

    DENDRO_103 *= DENDRO_25;
    DENDRO_103 *= DENDRO_123;
    DENDRO_103 *= DENDRO_122;

    // initialize reduction for
    DENDRO_104 = 1;

    DENDRO_104 *= DENDRO_25;
    DENDRO_104 *= DENDRO_10;
    DENDRO_104 *= DENDRO_124;

    // initialize reduction for
    DENDRO_105 = 1;

    DENDRO_105 *= DENDRO_25;
    DENDRO_105 *= DENDRO_56;
    DENDRO_105 *= DENDRO_55;

    // initialize reduction for
    DENDRO_106 = 1;

    DENDRO_106 *= DENDRO_25;
    DENDRO_106 *= DENDRO_58;
    DENDRO_106 *= DENDRO_57;

    // initialize reduction for
    DENDRO_107 = 1;

    DENDRO_107 *= DENDRO_25;
    DENDRO_107 *= DENDRO_27;
    DENDRO_107 *= DENDRO_59;

    // initialize reduction for
    DENDRO_59 = 1;

    DENDRO_59 *= DENDRO_25;
    DENDRO_59 *= DENDRO_27;
    DENDRO_59 *= DENDRO_60;

    // initialize reduction for
    DENDRO_60 = 1;

    DENDRO_60 *= DENDRO_25;
    DENDRO_60 *= DENDRO_27;
    DENDRO_60 *= DENDRO_61;

    // initialize reduction for
    DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_25;
    DENDRO_61 *= DENDRO_112;
    DENDRO_61 *= DENDRO_49;

    // initialize reduction for
    DENDRO_49 = 1;

    DENDRO_49 *= DENDRO_25;
    DENDRO_49 *= DENDRO_63;
    DENDRO_49 *= DENDRO_62;

    // initialize reduction for
    DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_25;
    DENDRO_62 *= DENDRO_65;
    DENDRO_62 *= DENDRO_64;

    // initialize reduction for
    DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_25;
    DENDRO_64 *= DENDRO_68;
    DENDRO_64 *= DENDRO_66;

    // initialize reduction for
    DENDRO_68 = 1;

    DENDRO_68 *= DENDRO_31;
    DENDRO_68 *= DENDRO_70;

    // initialize reduction for
    DENDRO_70 = 1;

    DENDRO_70 *= DENDRO_4;
    DENDRO_70 *= DENDRO_32;

    // initialize reduction for
    DENDRO_32 = 1;

    DENDRO_32 *= DENDRO_6;
    DENDRO_32 *= DENDRO_74;

    // initialize reduction for
    DENDRO_74 = 1;

    DENDRO_74 *= DENDRO_72;
    DENDRO_74 *= DENDRO_75;
    DENDRO_74 *= DENDRO_27;
    DENDRO_74 *= DENDRO_24;

    // initialize reduction for
    DENDRO_75 = 1;

    DENDRO_75 *= DENDRO_72;
    DENDRO_75 *= DENDRO_76;
    DENDRO_75 *= DENDRO_27;
    DENDRO_75 *= DENDRO_28;

    // initialize reduction for
    DENDRO_76 = 1;

    DENDRO_76 *= DENDRO_72;
    DENDRO_76 *= DENDRO_77;
    DENDRO_76 *= DENDRO_27;
    DENDRO_76 *= DENDRO_38;

    // initialize reduction for
    DENDRO_77 = 1;

    DENDRO_77 *= DENDRO_0;
    DENDRO_77 *= DENDRO_78;
    DENDRO_77 *= DENDRO_27;
    DENDRO_77 *= DENDRO_73;

    // initialize reduction for
    DENDRO_78 = 1;

    DENDRO_78 *= DENDRO_0;
    DENDRO_78 *= DENDRO_79;
    DENDRO_78 *= DENDRO_27;
    DENDRO_78 *= DENDRO_37;

    // initialize reduction for
    DENDRO_79 = 1;

    DENDRO_79 *= DENDRO_0;
    DENDRO_79 *= DENDRO_80;
    DENDRO_79 *= DENDRO_27;
    DENDRO_79 *= DENDRO_45;

    // initialize reduction for
    DENDRO_80 = 1;

    DENDRO_80 *= DENDRO_71;
    DENDRO_80 *= DENDRO_2;
    DENDRO_80 *= DENDRO_114;

    // initialize reduction for
    DENDRO_108 = 1;

    DENDRO_108 *= DENDRO_71;
    DENDRO_108 *= DENDRO_5;
    DENDRO_108 *= DENDRO_119;

    // initialize reduction for
    DENDRO_121 = 1;

    DENDRO_121 *= DENDRO_23;
    DENDRO_121 *= DENDRO_4;
    DENDRO_121 *= DENDRO_117;

    // initialize reduction for
    DENDRO_122 = 1;

    DENDRO_122 *= DENDRO_113;
    DENDRO_122 *= DENDRO_63;
    DENDRO_122 *= DENDRO_111;

    // initialize reduction for
    DENDRO_124 = 1;

    DENDRO_124 *= DENDRO_115;
    DENDRO_124 *= DENDRO_2;
    DENDRO_124 *= DENDRO_110;

    // initialize reduction for
    double DENDRO_125 = 1;

    DENDRO_125 *= DENDRO_115;
    DENDRO_125 *= DENDRO_2;
    DENDRO_125 *= DENDRO_116;

    // initialize reduction for
    double DENDRO_126 = 1;

    DENDRO_126 *= DENDRO_115;
    DENDRO_126 *= DENDRO_4;
    DENDRO_126 *= DENDRO_69;

    // initialize reduction for
    double DENDRO_127 = 1;

    DENDRO_127 *= DENDRO_115;
    DENDRO_127 *= DENDRO_4;
    DENDRO_127 *= DENDRO_118;

    // initialize reduction for
    double DENDRO_128 = 1;

    DENDRO_128 *= DENDRO_115;
    DENDRO_128 *= DENDRO_5;
    DENDRO_128 *= DENDRO_109;

    // initialize reduction for
    double DENDRO_129 = 1;

    DENDRO_129 *= DENDRO_115;
    DENDRO_129 *= DENDRO_5;
    DENDRO_129 *= DENDRO_120;

    // initialize reduction for
    double DENDRO_130 = 1;

    DENDRO_130 *= DENDRO_0;
    DENDRO_130 *= DENDRO_3;
    DENDRO_130 *= DENDRO_81;

    // initialize reduction for
    DENDRO_81 = 1;

    DENDRO_81 *= DENDRO_25;
    DENDRO_81 *= DENDRO_123;
    DENDRO_81 *= DENDRO_82;

    // initialize reduction for
    DENDRO_82 = 1;

    DENDRO_82 *= DENDRO_25;
    DENDRO_82 *= DENDRO_10;
    DENDRO_82 *= DENDRO_55;

    // initialize reduction for
    DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_25;
    DENDRO_55 *= DENDRO_56;
    DENDRO_55 *= DENDRO_83;

    // initialize reduction for
    DENDRO_83 = 1;

    DENDRO_83 *= DENDRO_25;
    DENDRO_83 *= DENDRO_58;
    DENDRO_83 *= DENDRO_12;

    // initialize reduction for
    double DENDRO_131 = 1;

    DENDRO_131 *= DENDRO_25;
    DENDRO_131 *= DENDRO_27;
    DENDRO_131 *= DENDRO_16;

    // initialize reduction for
    DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_25;
    DENDRO_16 *= DENDRO_27;
    DENDRO_16 *= DENDRO_36;

    // initialize reduction for
    DENDRO_36 = 1;

    DENDRO_36 *= DENDRO_25;
    DENDRO_36 *= DENDRO_27;
    DENDRO_36 *= DENDRO_84;

    // initialize reduction for
    DENDRO_84 = 1;

    DENDRO_84 *= DENDRO_25;
    DENDRO_84 *= DENDRO_112;
    DENDRO_84 *= DENDRO_85;

    // initialize reduction for
    DENDRO_85 = 1;

    DENDRO_85 *= DENDRO_25;
    DENDRO_85 *= DENDRO_63;
    DENDRO_85 *= DENDRO_50;

    // initialize reduction for
    DENDRO_50 = 1;

    DENDRO_50 *= DENDRO_25;
    DENDRO_50 *= DENDRO_65;
    DENDRO_50 *= DENDRO_86;

    // initialize reduction for
    DENDRO_86 = 1;

    DENDRO_86 *= DENDRO_25;
    DENDRO_86 *= DENDRO_87;
    DENDRO_86 *= DENDRO_66;

    // initialize reduction for
    DENDRO_87 = 1;

    DENDRO_87 *= DENDRO_31;
    DENDRO_87 *= DENDRO_88;

    // initialize reduction for
    DENDRO_88 = 1;

    DENDRO_88 *= DENDRO_1;
    DENDRO_88 *= DENDRO_39;

    // initialize reduction for
    DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_6;
    DENDRO_39 *= DENDRO_40;

    // initialize reduction for
    DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_72;
    DENDRO_40 *= DENDRO_89;
    DENDRO_40 *= DENDRO_27;
    DENDRO_40 *= DENDRO_24;

    // initialize reduction for
    DENDRO_24 = 1;

    DENDRO_24 *= DENDRO_72;
    DENDRO_24 *= DENDRO_90;
    DENDRO_24 *= DENDRO_27;
    DENDRO_24 *= DENDRO_28;

    // initialize reduction for
    DENDRO_28 = 1;

    DENDRO_28 *= DENDRO_72;
    DENDRO_28 *= DENDRO_91;
    DENDRO_28 *= DENDRO_27;
    DENDRO_28 *= DENDRO_38;

    // initialize reduction for
    DENDRO_38 = 1;

    DENDRO_38 *= DENDRO_0;
    DENDRO_38 *= DENDRO_92;
    DENDRO_38 *= DENDRO_27;
    DENDRO_38 *= DENDRO_73;

    // initialize reduction for
    DENDRO_72 = 1;

    DENDRO_72 *= DENDRO_0;
    DENDRO_72 *= DENDRO_93;
    DENDRO_72 *= DENDRO_27;
    DENDRO_72 *= DENDRO_37;

    // initialize reduction for
    DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_0;
    DENDRO_37 *= DENDRO_94;
    DENDRO_37 *= DENDRO_27;
    DENDRO_37 *= DENDRO_45;

    // initialize reduction for
    DENDRO_45 = 1;

    DENDRO_45 *= DENDRO_71;
    DENDRO_45 *= DENDRO_3;
    DENDRO_45 *= DENDRO_116;

    // initialize reduction for
    DENDRO_73 = 1;

    DENDRO_73 *= DENDRO_71;
    DENDRO_73 *= DENDRO_5;
    DENDRO_73 *= DENDRO_118;

    // initialize reduction for
    DENDRO_71 = 1;

    DENDRO_71 *= DENDRO_23;
    DENDRO_71 *= DENDRO_6;
    DENDRO_71 *= DENDRO_120;

    // initialize reduction for
    DENDRO_23 = 1;

    DENDRO_23 *= DENDRO_113;
    DENDRO_23 *= DENDRO_65;
    DENDRO_23 *= DENDRO_111;

    // initialize reduction for
    DENDRO_89 = 1;

    DENDRO_89 *= DENDRO_115;
    DENDRO_89 *= DENDRO_3;
    DENDRO_89 *= DENDRO_110;

    // initialize reduction for
    DENDRO_90 = 1;

    DENDRO_90 *= DENDRO_115;
    DENDRO_90 *= DENDRO_3;
    DENDRO_90 *= DENDRO_114;

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_115;
    DENDRO_3 *= DENDRO_5;
    DENDRO_3 *= DENDRO_69;

    // initialize reduction for
    DENDRO_69 = 1;

    DENDRO_69 *= DENDRO_115;
    DENDRO_69 *= DENDRO_5;
    DENDRO_69 *= DENDRO_117;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_115;
    DENDRO_5 *= DENDRO_6;
    DENDRO_5 *= DENDRO_109;

    // initialize reduction for
    DENDRO_91 = 1;

    DENDRO_91 *= DENDRO_115;
    DENDRO_91 *= DENDRO_6;
    DENDRO_91 *= DENDRO_119;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_0;
    DENDRO_6 *= DENDRO_2;
    DENDRO_6 *= DENDRO_95;

    // initialize reduction for
    DENDRO_0 = 1;

    DENDRO_0 *= DENDRO_25;
    DENDRO_0 *= DENDRO_123;
    DENDRO_0 *= DENDRO_96;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_25;
    DENDRO_2 *= DENDRO_10;
    DENDRO_2 *= DENDRO_57;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_25;
    DENDRO_10 *= DENDRO_56;
    DENDRO_10 *= DENDRO_12;

    // initialize reduction for
    DENDRO_12 = 1;

    DENDRO_12 *= DENDRO_25;
    DENDRO_12 *= DENDRO_58;
    DENDRO_12 *= DENDRO_97;

    // initialize reduction for
    DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_25;
    DENDRO_56 *= DENDRO_27;
    DENDRO_56 *= DENDRO_7;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_25;
    DENDRO_7 *= DENDRO_27;
    DENDRO_7 *= DENDRO_8;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_25;
    DENDRO_8 *= DENDRO_27;
    DENDRO_8 *= DENDRO_9;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_25;
    DENDRO_9 *= DENDRO_112;
    DENDRO_9 *= DENDRO_11;

    // initialize reduction for
    DENDRO_11 = 1;

    DENDRO_11 *= DENDRO_25;
    DENDRO_11 *= DENDRO_63;
    DENDRO_11 *= DENDRO_13;

    // initialize reduction for
    DENDRO_13 = 1;

    DENDRO_13 *= DENDRO_25;
    DENDRO_13 *= DENDRO_65;
    DENDRO_13 *= DENDRO_51;

    // initialize reduction for
    DENDRO_27 = 1;

    DENDRO_27 *= DENDRO_25;
    DENDRO_27 *= DENDRO_14;
    DENDRO_27 *= DENDRO_66;

    // initialize reduction for
    DENDRO_14 = 1;

    DENDRO_14 *= DENDRO_31;
    DENDRO_14 *= DENDRO_15;

    // initialize reduction for
    DENDRO_15 = 1;

    DENDRO_15 *= DENDRO_1;
    DENDRO_15 *= DENDRO_17;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_4;
    DENDRO_1 *= DENDRO_18;

    // initialize reduction for
    DENDRO_4 = 0;

    DENDRO_4 += DENDRO_32;
    DENDRO_4 += DENDRO_70;
    DENDRO_4 += DENDRO_54;
    DENDRO_4 += DENDRO_53;
    DENDRO_4 += DENDRO_52;
    DENDRO_4 += DENDRO_68;
    DENDRO_4 += DENDRO_64;
    DENDRO_4 += DENDRO_62;
    DENDRO_4 += DENDRO_49;
    DENDRO_4 += DENDRO_61;
    DENDRO_4 += DENDRO_60;
    DENDRO_4 += DENDRO_59;
    DENDRO_4 += DENDRO_107;
    DENDRO_4 += DENDRO_106;
    DENDRO_4 += DENDRO_105;
    DENDRO_4 += DENDRO_104;
    DENDRO_4 += DENDRO_103;
    DENDRO_4 += DENDRO_102;
    DENDRO_4 += DENDRO_101;
    DENDRO_4 += DENDRO_100;
    DENDRO_4 += DENDRO_99;
    DENDRO_4 += DENDRO_98;
    DENDRO_4 += DENDRO_67;
    DENDRO_4 += DENDRO_48;
    DENDRO_4 += DENDRO_47;
    DENDRO_4 += DENDRO_46;
    DENDRO_4 += DENDRO_44;
    DENDRO_4 += DENDRO_30;
    DENDRO_4 += DENDRO_29;
    DENDRO_4 += DENDRO_26;
    DENDRO_4 += DENDRO_22;
    DENDRO_4 += DENDRO_21;
    DENDRO_4 += DENDRO_20;
    DENDRO_4 += DENDRO_19;
    B_rhs2[pp] = DENDRO_4;

    // initialize reduction for
    DENDRO_4   = 0;

    DENDRO_4 += DENDRO_39;
    DENDRO_4 += DENDRO_88;
    DENDRO_4 += DENDRO_35;
    DENDRO_4 += DENDRO_34;
    DENDRO_4 += DENDRO_33;
    DENDRO_4 += DENDRO_87;
    DENDRO_4 += DENDRO_86;
    DENDRO_4 += DENDRO_50;
    DENDRO_4 += DENDRO_85;
    DENDRO_4 += DENDRO_84;
    DENDRO_4 += DENDRO_36;
    DENDRO_4 += DENDRO_16;
    DENDRO_4 += DENDRO_131;
    DENDRO_4 += DENDRO_83;
    DENDRO_4 += DENDRO_55;
    DENDRO_4 += DENDRO_82;
    DENDRO_4 += DENDRO_81;
    DENDRO_4 += DENDRO_130;
    DENDRO_4 += DENDRO_129;
    DENDRO_4 += DENDRO_128;
    DENDRO_4 += DENDRO_127;
    DENDRO_4 += DENDRO_126;
    DENDRO_4 += DENDRO_125;
    DENDRO_4 += DENDRO_124;
    DENDRO_4 += DENDRO_122;
    DENDRO_4 += DENDRO_121;
    DENDRO_4 += DENDRO_108;
    DENDRO_4 += DENDRO_80;
    DENDRO_4 += DENDRO_79;
    DENDRO_4 += DENDRO_78;
    DENDRO_4 += DENDRO_77;
    DENDRO_4 += DENDRO_76;
    DENDRO_4 += DENDRO_75;
    DENDRO_4 += DENDRO_74;
    B_rhs1[pp] = DENDRO_4;

    // initialize reduction for
    DENDRO_4   = 0;

    DENDRO_4 += DENDRO_1;
    DENDRO_4 += DENDRO_15;
    DENDRO_4 += DENDRO_43;
    DENDRO_4 += DENDRO_42;
    DENDRO_4 += DENDRO_41;
    DENDRO_4 += DENDRO_14;
    DENDRO_4 += DENDRO_27;
    DENDRO_4 += DENDRO_13;
    DENDRO_4 += DENDRO_11;
    DENDRO_4 += DENDRO_9;
    DENDRO_4 += DENDRO_8;
    DENDRO_4 += DENDRO_7;
    DENDRO_4 += DENDRO_56;
    DENDRO_4 += DENDRO_12;
    DENDRO_4 += DENDRO_10;
    DENDRO_4 += DENDRO_2;
    DENDRO_4 += DENDRO_0;
    DENDRO_4 += DENDRO_6;
    DENDRO_4 += DENDRO_91;
    DENDRO_4 += DENDRO_5;
    DENDRO_4 += DENDRO_69;
    DENDRO_4 += DENDRO_3;
    DENDRO_4 += DENDRO_90;
    DENDRO_4 += DENDRO_89;
    DENDRO_4 += DENDRO_23;
    DENDRO_4 += DENDRO_71;
    DENDRO_4 += DENDRO_73;
    DENDRO_4 += DENDRO_45;
    DENDRO_4 += DENDRO_37;
    DENDRO_4 += DENDRO_72;
    DENDRO_4 += DENDRO_38;
    DENDRO_4 += DENDRO_28;
    DENDRO_4 += DENDRO_24;
    DENDRO_4 += DENDRO_40;
    B_rhs0[pp] = DENDRO_4;

    // initialize reduction for
    DENDRO_0   = 0;

    DENDRO_0 += DENDRO_32;
    DENDRO_0 += DENDRO_70;
    DENDRO_0 += DENDRO_54;
    DENDRO_0 += DENDRO_53;
    DENDRO_0 += DENDRO_52;
    DENDRO_0 += DENDRO_62;
    DENDRO_0 += DENDRO_49;
    DENDRO_0 += DENDRO_61;
    DENDRO_0 += DENDRO_60;
    DENDRO_0 += DENDRO_59;
    DENDRO_0 += DENDRO_107;
    DENDRO_0 += DENDRO_106;
    DENDRO_0 += DENDRO_105;
    DENDRO_0 += DENDRO_104;
    DENDRO_0 += DENDRO_102;
    DENDRO_0 += DENDRO_101;
    DENDRO_0 += DENDRO_100;
    DENDRO_0 += DENDRO_99;
    DENDRO_0 += DENDRO_98;
    DENDRO_0 += DENDRO_67;
    DENDRO_0 += DENDRO_48;
    DENDRO_0 += DENDRO_47;
    DENDRO_0 += DENDRO_46;
    DENDRO_0 += DENDRO_44;
    DENDRO_0 += DENDRO_30;
    DENDRO_0 += DENDRO_29;
    DENDRO_0 += DENDRO_26;
    DENDRO_0 += DENDRO_22;
    DENDRO_0 += DENDRO_21;
    DENDRO_0 += DENDRO_20;
    DENDRO_0 += DENDRO_19;
    Gt_rhs2[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0    = 0;

    DENDRO_0 += DENDRO_39;
    DENDRO_0 += DENDRO_88;
    DENDRO_0 += DENDRO_35;
    DENDRO_0 += DENDRO_34;
    DENDRO_0 += DENDRO_33;
    DENDRO_0 += DENDRO_50;
    DENDRO_0 += DENDRO_85;
    DENDRO_0 += DENDRO_84;
    DENDRO_0 += DENDRO_36;
    DENDRO_0 += DENDRO_16;
    DENDRO_0 += DENDRO_131;
    DENDRO_0 += DENDRO_83;
    DENDRO_0 += DENDRO_55;
    DENDRO_0 += DENDRO_82;
    DENDRO_0 += DENDRO_130;
    DENDRO_0 += DENDRO_129;
    DENDRO_0 += DENDRO_128;
    DENDRO_0 += DENDRO_127;
    DENDRO_0 += DENDRO_126;
    DENDRO_0 += DENDRO_125;
    DENDRO_0 += DENDRO_124;
    DENDRO_0 += DENDRO_122;
    DENDRO_0 += DENDRO_121;
    DENDRO_0 += DENDRO_108;
    DENDRO_0 += DENDRO_80;
    DENDRO_0 += DENDRO_79;
    DENDRO_0 += DENDRO_78;
    DENDRO_0 += DENDRO_77;
    DENDRO_0 += DENDRO_76;
    DENDRO_0 += DENDRO_75;
    DENDRO_0 += DENDRO_74;
    Gt_rhs1[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0    = 0;

    DENDRO_0 += DENDRO_1;
    DENDRO_0 += DENDRO_15;
    DENDRO_0 += DENDRO_43;
    DENDRO_0 += DENDRO_42;
    DENDRO_0 += DENDRO_41;
    DENDRO_0 += DENDRO_13;
    DENDRO_0 += DENDRO_11;
    DENDRO_0 += DENDRO_9;
    DENDRO_0 += DENDRO_8;
    DENDRO_0 += DENDRO_7;
    DENDRO_0 += DENDRO_56;
    DENDRO_0 += DENDRO_12;
    DENDRO_0 += DENDRO_10;
    DENDRO_0 += DENDRO_2;
    DENDRO_0 += DENDRO_6;
    DENDRO_0 += DENDRO_91;
    DENDRO_0 += DENDRO_5;
    DENDRO_0 += DENDRO_69;
    DENDRO_0 += DENDRO_3;
    DENDRO_0 += DENDRO_90;
    DENDRO_0 += DENDRO_89;
    DENDRO_0 += DENDRO_23;
    DENDRO_0 += DENDRO_71;
    DENDRO_0 += DENDRO_73;
    DENDRO_0 += DENDRO_45;
    DENDRO_0 += DENDRO_37;
    DENDRO_0 += DENDRO_72;
    DENDRO_0 += DENDRO_38;
    DENDRO_0 += DENDRO_28;
    DENDRO_0 += DENDRO_24;
    DENDRO_0 += DENDRO_40;
    Gt_rhs0[pp] = DENDRO_0;
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
        double DENDRO_0 = 2;
        double DENDRO_1 = gt1[pp];
        double DENDRO_2 = gt2[pp];
        double DENDRO_3 = gt4[pp];
        double DENDRO_4 = -2;
        double DENDRO_5 = gt5[pp];
        double DENDRO_6 = gt3[pp];
        double DENDRO_7 = gt0[pp];
        double DENDRO_8 = -1;
        double DENDRO_9 = DENDRO_1;

        DENDRO_9 *= DENDRO_1;
        double DENDRO_10 = DENDRO_2;

        DENDRO_10 *= DENDRO_2;
        double DENDRO_11 = DENDRO_3;

        DENDRO_11 *= DENDRO_3;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_4;
        DENDRO_0 *= DENDRO_1;
        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_8;
        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_5;
        DENDRO_12 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_6;
        DENDRO_13 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_8;
        DENDRO_15 *= DENDRO_7;
        DENDRO_15 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_16 = 0;

        DENDRO_16 += DENDRO_14;
        DENDRO_16 += DENDRO_13;
        DENDRO_16 += DENDRO_12;
        DENDRO_16 += DENDRO_4;
        DENDRO_16 += DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_8;
        DENDRO_0 *= DENDRO_1;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_3;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_8;
        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_3;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_2;
        DENDRO_13 *= DENDRO_6;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_8;
        DENDRO_14 *= DENDRO_2;
        DENDRO_14 *= DENDRO_3;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_8;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_9;
        DENDRO_3 += DENDRO_15;
        DENDRO_5 = DENDRO_16;

        DENDRO_5 = 1 / DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_4;
        DENDRO_6 += DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 0;

        DENDRO_0 += DENDRO_10;
        DENDRO_0 += DENDRO_12;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_13;
        DENDRO_4 += DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 0;

        DENDRO_7 += DENDRO_2;
        DENDRO_7 += DENDRO_14;

        // initialize reduction for
        DENDRO_2 = 0;

        DENDRO_2 += DENDRO_11;
        DENDRO_2 += DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_3;
        DENDRO_igt5 = DENDRO_1;

        // initialize reduction for
        DENDRO_1    = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_6;
        DENDRO_igt4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1    = 1;

        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_0;
        DENDRO_igt3 = DENDRO_1;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_4;
        DENDRO_igt2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_7;
        DENDRO_igt1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0    = 1;

        DENDRO_0 *= DENDRO_5;
        DENDRO_0 *= DENDRO_2;
        DENDRO_igt0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_0_gt5;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_2_gt2;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_0_gt4;
        double DENDRO_5  = grad_2_gt1;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_1_gt2;
        double DENDRO_8  = grad_0_gt3;
        double DENDRO_9  = grad_1_gt1;

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;
        DENDRO_3         = grad_2_gt0;
        DENDRO_8         = grad_1_gt0;
        DENDRO_9         = grad_0_gt0;

        // initialize reduction for
        double DENDRO_11 = 0;

        DENDRO_11 += DENDRO_0;
        DENDRO_11 += DENDRO_10;
        DENDRO_C1_k0_5 = DENDRO_11;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_2;
        DENDRO_C1_k0_4 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_7;
        DENDRO_C1_k0_3 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_3;
        DENDRO_C1_k0_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_8;
        DENDRO_C1_k0_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_9;
        DENDRO_C1_k0_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_1_gt5;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_2_gt4;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_1_gt2;
        double DENDRO_5  = grad_2_gt1;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_0_gt4;
        double DENDRO_8  = grad_1_gt0;
        double DENDRO_9  = grad_0_gt1;

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;
        DENDRO_2         = grad_2_gt3;
        double DENDRO_11 = grad_1_gt3;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;
        DENDRO_7         = grad_0_gt3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_0;
        DENDRO_3 += DENDRO_10;
        DENDRO_C1_k1_5 = DENDRO_3;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_2;
        DENDRO_C1_k1_4 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_11;
        DENDRO_C1_k1_3 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_12;
        DENDRO_C1_k1_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_7;
        DENDRO_C1_k1_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_13;
        DENDRO_C1_k1_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = grad_2_gt3;
        double DENDRO_1  = -0.5;
        double DENDRO_2  = grad_1_gt4;
        double DENDRO_3  = 1.0;
        double DENDRO_4  = grad_2_gt1;
        double DENDRO_5  = grad_1_gt2;
        double DENDRO_6  = 0.5;
        double DENDRO_7  = grad_0_gt4;
        double DENDRO_8  = grad_2_gt0;
        double DENDRO_9  = grad_0_gt2;
        double DENDRO_10 = grad_2_gt5;
        double DENDRO_11 = grad_1_gt5;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_2;
        DENDRO_2         = grad_0_gt5;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_7;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_1;
        DENDRO_7 *= DENDRO_8;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_9;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_6;
        DENDRO_3 *= DENDRO_10;
        DENDRO_C1_k2_5 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 1;

        DENDRO_3 *= DENDRO_6;
        DENDRO_3 *= DENDRO_11;
        DENDRO_C1_k2_4 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_0;
        DENDRO_3 += DENDRO_12;
        DENDRO_C1_k2_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_0       = 1;

        DENDRO_0 *= DENDRO_6;
        DENDRO_0 *= DENDRO_2;
        DENDRO_C1_k2_2 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_5;
        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_13;
        DENDRO_C1_k2_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_7;
        DENDRO_C1_k2_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = DENDRO_igt2;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt1;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt0;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k0_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k0_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k0_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k0_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k0_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k0_0 = DENDRO_1;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt3;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt1;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k1_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k1_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k1_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k1_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k1_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k1_0 = DENDRO_1;
    }
    {
        double DENDRO_0  = DENDRO_igt5;
        double DENDRO_1  = DENDRO_C1_k2_5;
        double DENDRO_2  = DENDRO_igt4;
        double DENDRO_3  = DENDRO_C1_k1_5;
        double DENDRO_4  = DENDRO_igt2;
        double DENDRO_5  = DENDRO_C1_k0_5;
        double DENDRO_6  = DENDRO_C1_k2_4;
        double DENDRO_7  = DENDRO_C1_k1_4;
        double DENDRO_8  = DENDRO_C1_k0_4;
        double DENDRO_9  = DENDRO_C1_k2_3;
        double DENDRO_10 = DENDRO_C1_k1_3;
        double DENDRO_11 = DENDRO_C1_k0_3;
        double DENDRO_12 = DENDRO_C1_k2_2;
        double DENDRO_13 = DENDRO_C1_k1_2;
        double DENDRO_14 = DENDRO_C1_k0_2;
        double DENDRO_15 = DENDRO_C1_k2_1;
        double DENDRO_16 = DENDRO_C1_k1_1;
        double DENDRO_17 = DENDRO_C1_k0_1;
        double DENDRO_18 = DENDRO_C1_k2_0;
        double DENDRO_19 = DENDRO_C1_k1_0;
        double DENDRO_20 = DENDRO_C1_k0_0;

        // initialize reduction for
        double DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_1;
        DENDRO_21 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_6;
        DENDRO_5 *= DENDRO_0;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_9;
        DENDRO_8 *= DENDRO_0;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_10;
        DENDRO_9 *= DENDRO_2;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_12;
        DENDRO_11 *= DENDRO_0;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_2;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_4;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_0;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_2;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_4;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_0;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_19;
        DENDRO_0 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_20;
        DENDRO_2 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_21;
        DENDRO_C2_k2_5 = DENDRO_4;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_1 += DENDRO_5;
        DENDRO_C2_k2_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_9;
        DENDRO_1 += DENDRO_8;
        DENDRO_C2_k2_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_11;
        DENDRO_C2_k2_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_14;
        DENDRO_C2_k2_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_2;
        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_17;
        DENDRO_C2_k2_0 = DENDRO_1;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt2;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt1;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt0;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1         = gt2[pp];
        DENDRO_3         = -1;
        DENDRO_6         = gt1[pp];
        double DENDRO_7  = gt0[pp];
        double DENDRO_8  = 2;
        double DENDRO_9  = chi[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_3;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_3;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_4;
        DENDRO_4         = DENDRO_9;

        DENDRO_4         = 1 / DENDRO_4;
        DENDRO_3         = gt5[pp];
        DENDRO_8         = 0.5;
        DENDRO_9         = gt4[pp];
        double DENDRO_11 = gt3[pp];

        // initialize reduction for
        double DENDRO_12 = 0;

        DENDRO_12 += DENDRO_0;
        DENDRO_12 += DENDRO_10;
        DENDRO_0  = -0.5;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_2;
        DENDRO_10 += DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_8;
        DENDRO_2 *= DENDRO_3;
        DENDRO_2 *= DENDRO_4;
        DENDRO_2 *= DENDRO_5;
        DENDRO_3 = DENDRO_C2_k0_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_9;
        DENDRO_6 *= DENDRO_4;
        DENDRO_6 *= DENDRO_5;
        DENDRO_7 = DENDRO_C2_k0_4;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_8;
        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_4;
        DENDRO_9 *= DENDRO_5;
        DENDRO_5 = DENDRO_C2_k0_3;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_0;
        DENDRO_8 *= DENDRO_4;
        DENDRO_8 *= DENDRO_12;
        DENDRO_11 = DENDRO_C2_k0_2;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_0;
        DENDRO_12 *= DENDRO_4;
        DENDRO_12 *= DENDRO_10;
        DENDRO_10        = DENDRO_C2_k0_1;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_0;
        DENDRO_13 *= DENDRO_4;
        DENDRO_13 *= DENDRO_1;
        DENDRO_0 = DENDRO_C2_k0_0;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_3;
        DENDRO_1 += DENDRO_2;
        DENDRO_C3_k0_5 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_7;
        DENDRO_1 += DENDRO_6;
        DENDRO_C3_k0_4 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_5;
        DENDRO_1 += DENDRO_9;
        DENDRO_C3_k0_3 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_11;
        DENDRO_1 += DENDRO_8;
        DENDRO_C3_k0_2 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_12;
        DENDRO_C3_k0_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_13;
        DENDRO_C3_k0_0 = DENDRO_1;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt4;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt3;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt1;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1         = gt4[pp];
        DENDRO_3         = -1;
        DENDRO_6         = gt3[pp];
        double DENDRO_7  = 2;
        double DENDRO_8  = gt1[pp];
        double DENDRO_9  = chi[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_3;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_2;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_3;
        DENDRO_2 *= DENDRO_8;
        DENDRO_2 *= DENDRO_5;
        DENDRO_7 = DENDRO_9;

        DENDRO_7 = 1 / DENDRO_7;
        DENDRO_3 = gt5[pp];
        DENDRO_8 = 0.5;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_0;
        DENDRO_9 += DENDRO_10;
        DENDRO_0  = -0.5;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_6;
        DENDRO_10 += DENDRO_1;
        DENDRO_1 = gt2[pp];

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_4;
        DENDRO_6 += DENDRO_2;
        DENDRO_2 = gt0[pp];

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_8;
        DENDRO_4 *= DENDRO_3;
        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_5;
        DENDRO_3         = DENDRO_C2_k1_5;

        // initialize reduction for
        double DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_0;
        DENDRO_11 *= DENDRO_7;
        DENDRO_11 *= DENDRO_9;
        DENDRO_9         = DENDRO_C2_k1_4;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_0;
        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_10;
        DENDRO_10        = DENDRO_C2_k1_3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_8;
        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_7;
        DENDRO_13 *= DENDRO_5;
        DENDRO_1         = DENDRO_C2_k1_2;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_0;
        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_6;
        DENDRO_0 = DENDRO_C2_k1_1;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_2;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;
        DENDRO_2 = DENDRO_C2_k1_0;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_4;
        DENDRO_C3_k1_5 = DENDRO_5;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_9;
        DENDRO_3 += DENDRO_11;
        DENDRO_C3_k1_4 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_10;
        DENDRO_3 += DENDRO_12;
        DENDRO_C3_k1_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_1;
        DENDRO_3 += DENDRO_13;
        DENDRO_C3_k1_2 = DENDRO_3;

        // initialize reduction for
        DENDRO_1       = 0;

        DENDRO_1 += DENDRO_0;
        DENDRO_1 += DENDRO_14;
        DENDRO_C3_k1_1 = DENDRO_1;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_2;
        DENDRO_0 += DENDRO_6;
        DENDRO_C3_k1_0 = DENDRO_0;
    }
    {
        double DENDRO_0 = grad_2_chi;
        double DENDRO_1 = DENDRO_igt5;
        double DENDRO_2 = grad_1_chi;
        double DENDRO_3 = DENDRO_igt4;
        double DENDRO_4 = grad_0_chi;
        double DENDRO_5 = DENDRO_igt2;

        // initialize reduction for
        double DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_1;
        DENDRO_6 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_2;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_1;
        DENDRO_5 += DENDRO_6;
        DENDRO_1        = gt5[pp];
        DENDRO_3        = -1;
        DENDRO_6        = 2;
        double DENDRO_7 = gt4[pp];
        double DENDRO_8 = gt2[pp];

        // initialize reduction for
        double DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_3;
        DENDRO_9 *= DENDRO_1;
        DENDRO_9 *= DENDRO_5;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_6;
        DENDRO_1 *= DENDRO_0;
        DENDRO_0 = chi[pp];

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_3;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_5;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_3;
        DENDRO_7 *= DENDRO_8;
        DENDRO_7 *= DENDRO_5;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_1;
        DENDRO_8 += DENDRO_9;
        DENDRO_1 = DENDRO_0;

        DENDRO_1 = 1 / DENDRO_1;
        DENDRO_0 = -0.5;

        // initialize reduction for
        DENDRO_3 = 0;

        DENDRO_3 += DENDRO_2;
        DENDRO_3 += DENDRO_6;
        DENDRO_2 = gt3[pp];
        DENDRO_6 = 0.5;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_4;
        DENDRO_9 += DENDRO_7;
        DENDRO_4         = gt1[pp];
        DENDRO_7         = gt0[pp];

        // initialize reduction for
        double DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_0;
        DENDRO_10 *= DENDRO_1;
        DENDRO_10 *= DENDRO_8;
        DENDRO_8         = DENDRO_C2_k2_5;

        // initialize reduction for
        double DENDRO_11 = 1;

        DENDRO_11 *= DENDRO_0;
        DENDRO_11 *= DENDRO_1;
        DENDRO_11 *= DENDRO_3;
        DENDRO_3         = DENDRO_C2_k2_4;

        // initialize reduction for
        double DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_6;
        DENDRO_12 *= DENDRO_2;
        DENDRO_12 *= DENDRO_1;
        DENDRO_12 *= DENDRO_5;
        DENDRO_2         = DENDRO_C2_k2_3;

        // initialize reduction for
        double DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_0;
        DENDRO_13 *= DENDRO_1;
        DENDRO_13 *= DENDRO_9;
        DENDRO_0 = DENDRO_C2_k2_2;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_6;
        DENDRO_9 *= DENDRO_4;
        DENDRO_9 *= DENDRO_1;
        DENDRO_9 *= DENDRO_5;
        DENDRO_4         = DENDRO_C2_k2_1;

        // initialize reduction for
        double DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_6;
        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_1;
        DENDRO_14 *= DENDRO_5;
        DENDRO_1 = DENDRO_C2_k2_0;

        // initialize reduction for
        DENDRO_5 = 0;

        DENDRO_5 += DENDRO_8;
        DENDRO_5 += DENDRO_10;
        DENDRO_C3_k2_5 = DENDRO_5;

        // initialize reduction for
        DENDRO_5       = 0;

        DENDRO_5 += DENDRO_3;
        DENDRO_5 += DENDRO_11;
        DENDRO_C3_k2_4 = DENDRO_5;

        // initialize reduction for
        DENDRO_3       = 0;

        DENDRO_3 += DENDRO_2;
        DENDRO_3 += DENDRO_12;
        DENDRO_C3_k2_3 = DENDRO_3;

        // initialize reduction for
        DENDRO_2       = 0;

        DENDRO_2 += DENDRO_0;
        DENDRO_2 += DENDRO_13;
        DENDRO_C3_k2_2 = DENDRO_2;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_4;
        DENDRO_0 += DENDRO_9;
        DENDRO_C3_k2_1 = DENDRO_0;

        // initialize reduction for
        DENDRO_0       = 0;

        DENDRO_0 += DENDRO_1;
        DENDRO_0 += DENDRO_14;
        DENDRO_C3_k2_0 = DENDRO_0;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C2_k2_4;
        double DENDRO_2  = 2;
        double DENDRO_3  = DENDRO_igt2;
        double DENDRO_4  = DENDRO_C2_k2_2;
        double DENDRO_5  = DENDRO_igt1;
        double DENDRO_6  = DENDRO_C2_k2_1;
        double DENDRO_7  = DENDRO_igt5;
        double DENDRO_8  = DENDRO_C2_k2_5;
        double DENDRO_9  = DENDRO_igt3;
        double DENDRO_10 = DENDRO_C2_k2_3;
        double DENDRO_11 = DENDRO_igt0;
        double DENDRO_12 = DENDRO_C2_k2_0;
        double DENDRO_13 = DENDRO_C2_k1_4;
        double DENDRO_14 = DENDRO_C2_k1_2;
        double DENDRO_15 = DENDRO_C2_k1_1;
        double DENDRO_16 = DENDRO_C2_k1_5;
        double DENDRO_17 = DENDRO_C2_k1_3;
        double DENDRO_18 = DENDRO_C2_k1_0;
        double DENDRO_19 = DENDRO_C2_k0_4;
        double DENDRO_20 = DENDRO_C2_k0_2;
        double DENDRO_21 = DENDRO_C2_k0_1;
        double DENDRO_22 = DENDRO_C2_k0_5;
        double DENDRO_23 = DENDRO_C2_k0_3;
        double DENDRO_24 = DENDRO_C2_k0_0;

        // initialize reduction for
        double DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_1;
        DENDRO_25 *= DENDRO_0;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_2;
        DENDRO_1 *= DENDRO_4;
        DENDRO_1 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_6;
        DENDRO_26 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_8;
        DENDRO_27 *= DENDRO_7;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_10;
        DENDRO_8 *= DENDRO_9;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_12;
        DENDRO_10 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_2;
        DENDRO_28 *= DENDRO_13;
        DENDRO_28 *= DENDRO_0;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_2;
        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_2;
        DENDRO_29 *= DENDRO_15;
        DENDRO_29 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_16;
        DENDRO_30 *= DENDRO_7;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_9;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_2;
        DENDRO_31 *= DENDRO_19;
        DENDRO_31 *= DENDRO_0;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_2;
        DENDRO_19 *= DENDRO_20;
        DENDRO_19 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_21;
        DENDRO_32 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_22;
        DENDRO_33 *= DENDRO_7;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_23;
        DENDRO_22 *= DENDRO_9;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_24;
        DENDRO_23 *= DENDRO_11;
        double DENDRO_34 = grad_2_chi;
        double DENDRO_35 = grad_1_chi;
        double DENDRO_36 = grad_0_chi;
        double DENDRO_37 = DENDRO_C1_k0_5;
        double DENDRO_38 = DENDRO_C1_k0_4;
        double DENDRO_39 = DENDRO_C1_k0_2;
        double DENDRO_40 = DENDRO_C1_k0_3;
        double DENDRO_41 = DENDRO_C1_k0_1;
        double DENDRO_42 = DENDRO_C1_k2_2;
        double DENDRO_43 = DENDRO_C1_k2_1;
        double DENDRO_44 = DENDRO_C1_k1_2;
        double DENDRO_45 = DENDRO_C1_k1_1;
        double DENDRO_46 = DENDRO_C1_k2_0;
        double DENDRO_47 = DENDRO_C1_k1_0;
        double DENDRO_48 = DENDRO_C1_k0_0;
        double DENDRO_49 = grad2_1_2_chi;
        double DENDRO_50 = chi[pp];
        double DENDRO_51 = -3;
        double DENDRO_52 = grad2_0_2_chi;
        double DENDRO_53 = grad2_0_1_chi;

        // initialize reduction for
        double DENDRO_54 = 0;

        DENDRO_54 += DENDRO_10;
        DENDRO_54 += DENDRO_8;
        DENDRO_54 += DENDRO_27;
        DENDRO_54 += DENDRO_26;
        DENDRO_54 += DENDRO_1;
        DENDRO_54 += DENDRO_25;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_17;
        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_30;
        DENDRO_1 += DENDRO_29;
        DENDRO_1 += DENDRO_13;
        DENDRO_1 += DENDRO_28;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_23;
        DENDRO_8 += DENDRO_22;
        DENDRO_8 += DENDRO_33;
        DENDRO_8 += DENDRO_32;
        DENDRO_8 += DENDRO_19;
        DENDRO_8 += DENDRO_31;
        DENDRO_10 = grad2_2_2_chi;
        DENDRO_13 = DENDRO_34;

        DENDRO_13 *= DENDRO_34;
        DENDRO_16 = grad2_1_1_chi;
        DENDRO_17 = DENDRO_35;

        DENDRO_17 *= DENDRO_35;
        DENDRO_19 = grad2_0_0_chi;
        DENDRO_22 = DENDRO_36;

        DENDRO_22 *= DENDRO_36;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_2;
        DENDRO_23 *= DENDRO_37;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_38;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_39;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_40;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_2;
        DENDRO_28 *= DENDRO_41;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_2;
        DENDRO_29 *= DENDRO_38;
        DENDRO_29 *= DENDRO_4;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_42;
        DENDRO_30 *= DENDRO_6;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_2;
        DENDRO_31 *= DENDRO_37;
        DENDRO_31 *= DENDRO_6;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_43;
        DENDRO_32 *= DENDRO_4;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_2;
        DENDRO_33 *= DENDRO_40;
        DENDRO_33 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_44;
        DENDRO_55 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_2;
        DENDRO_56 *= DENDRO_38;
        DENDRO_56 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_45;
        DENDRO_57 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_2;
        DENDRO_58 *= DENDRO_41;
        DENDRO_58 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_59 = 1;

        DENDRO_59 *= DENDRO_39;
        DENDRO_59 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_60 = 1;

        DENDRO_60 *= DENDRO_2;
        DENDRO_60 *= DENDRO_39;
        DENDRO_60 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_61 = 1;

        DENDRO_61 *= DENDRO_41;
        DENDRO_61 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_62 = 1;

        DENDRO_62 *= DENDRO_2;
        DENDRO_62 *= DENDRO_39;
        DENDRO_62 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_63 = 1;

        DENDRO_63 *= DENDRO_42;
        DENDRO_63 *= DENDRO_12;

        // initialize reduction for
        double DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_2;
        DENDRO_64 *= DENDRO_37;
        DENDRO_64 *= DENDRO_12;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_46;
        DENDRO_37 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_2;
        DENDRO_65 *= DENDRO_41;
        DENDRO_65 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_44;
        DENDRO_66 *= DENDRO_18;

        // initialize reduction for
        double DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_2;
        DENDRO_67 *= DENDRO_38;
        DENDRO_67 *= DENDRO_18;

        // initialize reduction for
        double DENDRO_68 = 1;

        DENDRO_68 *= DENDRO_47;
        DENDRO_68 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_69 = 1;

        DENDRO_69 *= DENDRO_2;
        DENDRO_69 *= DENDRO_48;
        DENDRO_69 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_70 = 1;

        DENDRO_70 *= DENDRO_39;
        DENDRO_70 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_71 = 1;

        DENDRO_71 *= DENDRO_2;
        DENDRO_71 *= DENDRO_39;
        DENDRO_71 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_72 = 1;

        DENDRO_72 *= DENDRO_48;
        DENDRO_72 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_73 = 1;

        DENDRO_73 *= DENDRO_2;
        DENDRO_73 *= DENDRO_39;
        DENDRO_73 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_74 = 1;

        DENDRO_74 *= DENDRO_43;
        DENDRO_74 *= DENDRO_12;

        // initialize reduction for
        double DENDRO_75 = 1;

        DENDRO_75 *= DENDRO_2;
        DENDRO_75 *= DENDRO_38;
        DENDRO_75 *= DENDRO_12;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_46;
        DENDRO_38 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_76 = 1;

        DENDRO_76 *= DENDRO_2;
        DENDRO_76 *= DENDRO_41;
        DENDRO_76 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_77 = 1;

        DENDRO_77 *= DENDRO_45;
        DENDRO_77 *= DENDRO_18;

        // initialize reduction for
        double DENDRO_78 = 1;

        DENDRO_78 *= DENDRO_2;
        DENDRO_78 *= DENDRO_40;
        DENDRO_78 *= DENDRO_18;

        // initialize reduction for
        DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_47;
        DENDRO_40 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_79 = 1;

        DENDRO_79 *= DENDRO_2;
        DENDRO_79 *= DENDRO_48;
        DENDRO_79 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_80 = 1;

        DENDRO_80 *= DENDRO_41;
        DENDRO_80 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_81 = 1;

        DENDRO_81 *= DENDRO_2;
        DENDRO_81 *= DENDRO_41;
        DENDRO_81 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_82 = 1;

        DENDRO_82 *= DENDRO_48;
        DENDRO_82 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_83 = 1;

        DENDRO_83 *= DENDRO_2;
        DENDRO_83 *= DENDRO_50;
        DENDRO_83 *= DENDRO_49;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_51;
        DENDRO_49 *= DENDRO_35;
        DENDRO_49 *= DENDRO_34;

        // initialize reduction for
        double DENDRO_84 = 1;

        DENDRO_84 *= DENDRO_2;
        DENDRO_84 *= DENDRO_50;
        DENDRO_84 *= DENDRO_52;

        // initialize reduction for
        DENDRO_52 = 1;

        DENDRO_52 *= DENDRO_51;
        DENDRO_52 *= DENDRO_36;
        DENDRO_52 *= DENDRO_34;

        // initialize reduction for
        double DENDRO_85 = 1;

        DENDRO_85 *= DENDRO_2;
        DENDRO_85 *= DENDRO_50;
        DENDRO_85 *= DENDRO_53;

        // initialize reduction for
        DENDRO_53 = 1;

        DENDRO_53 *= DENDRO_51;
        DENDRO_53 *= DENDRO_36;
        DENDRO_53 *= DENDRO_35;

        // initialize reduction for
        double DENDRO_86 = 1;

        DENDRO_86 *= DENDRO_34;
        DENDRO_86 *= DENDRO_54;

        // initialize reduction for
        double DENDRO_87 = 1;

        DENDRO_87 *= DENDRO_35;
        DENDRO_87 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_88 = 1;

        DENDRO_88 *= DENDRO_36;
        DENDRO_88 *= DENDRO_8;

        // initialize reduction for
        double DENDRO_89 = 1;

        DENDRO_89 *= DENDRO_2;
        DENDRO_89 *= DENDRO_50;
        DENDRO_89 *= DENDRO_10;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_51;
        DENDRO_10 *= DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_2;
        DENDRO_13 *= DENDRO_50;
        DENDRO_13 *= DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_51;
        DENDRO_16 *= DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_2;
        DENDRO_17 *= DENDRO_50;
        DENDRO_17 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_90 = 1;

        DENDRO_90 *= DENDRO_51;
        DENDRO_90 *= DENDRO_22;
        DENDRO_51        = -1;
        double DENDRO_91 = 12;

        // initialize reduction for
        double DENDRO_92 = 0;

        DENDRO_92 += DENDRO_42;
        DENDRO_92 += DENDRO_23;
        DENDRO_23 = 4;

        // initialize reduction for
        DENDRO_42 = 0;

        DENDRO_42 += DENDRO_43;
        DENDRO_42 += DENDRO_25;

        // initialize reduction for
        DENDRO_43 = 0;

        DENDRO_43 += DENDRO_46;
        DENDRO_43 += DENDRO_26;

        // initialize reduction for
        DENDRO_26 = 0;

        DENDRO_26 += DENDRO_44;
        DENDRO_26 += DENDRO_25;

        // initialize reduction for
        DENDRO_25 = 0;

        DENDRO_25 += DENDRO_45;
        DENDRO_25 += DENDRO_27;

        // initialize reduction for
        DENDRO_27 = 0;

        DENDRO_27 += DENDRO_47;
        DENDRO_27 += DENDRO_28;
        DENDRO_28         = grad2_1_2_gt0;
        DENDRO_44         = -4.0;
        DENDRO_45         = grad2_0_2_gt0;
        DENDRO_46         = grad2_0_1_gt0;
        DENDRO_47         = grad2_2_2_gt0;
        double DENDRO_93  = -2.0;
        double DENDRO_94  = grad2_1_1_gt0;
        double DENDRO_95  = grad2_0_0_gt0;
        double DENDRO_96  = gt2[pp];
        double DENDRO_97  = grad_0_Gt2;
        double DENDRO_98  = 4.0;
        double DENDRO_99  = gt1[pp];
        double DENDRO_100 = grad_0_Gt1;
        double DENDRO_101 = gt0[pp];
        double DENDRO_102 = grad_0_Gt0;

        // initialize reduction for
        double DENDRO_103 = 0;

        DENDRO_103 += DENDRO_30;
        DENDRO_103 += DENDRO_29;

        // initialize reduction for
        DENDRO_29 = 0;

        DENDRO_29 += DENDRO_32;
        DENDRO_29 += DENDRO_31;

        // initialize reduction for
        DENDRO_30 = 0;

        DENDRO_30 += DENDRO_55;
        DENDRO_30 += DENDRO_33;

        // initialize reduction for
        DENDRO_31 = 0;

        DENDRO_31 += DENDRO_57;
        DENDRO_31 += DENDRO_56;

        // initialize reduction for
        DENDRO_32 = 0;

        DENDRO_32 += DENDRO_59;
        DENDRO_32 += DENDRO_58;

        // initialize reduction for
        DENDRO_33 = 0;

        DENDRO_33 += DENDRO_61;
        DENDRO_33 += DENDRO_60;

        // initialize reduction for
        DENDRO_55 = 0;

        DENDRO_55 += DENDRO_63;
        DENDRO_55 += DENDRO_62;

        // initialize reduction for
        DENDRO_56 = 0;

        DENDRO_56 += DENDRO_37;
        DENDRO_56 += DENDRO_64;

        // initialize reduction for
        DENDRO_37 = 0;

        DENDRO_37 += DENDRO_66;
        DENDRO_37 += DENDRO_65;

        // initialize reduction for
        DENDRO_57 = 0;

        DENDRO_57 += DENDRO_68;
        DENDRO_57 += DENDRO_67;

        // initialize reduction for
        DENDRO_58 = 0;

        DENDRO_58 += DENDRO_70;
        DENDRO_58 += DENDRO_69;

        // initialize reduction for
        DENDRO_59 = 0;

        DENDRO_59 += DENDRO_72;
        DENDRO_59 += DENDRO_71;

        // initialize reduction for
        DENDRO_60 = 0;

        DENDRO_60 += DENDRO_74;
        DENDRO_60 += DENDRO_73;

        // initialize reduction for
        DENDRO_61 = 0;

        DENDRO_61 += DENDRO_38;
        DENDRO_61 += DENDRO_75;

        // initialize reduction for
        DENDRO_38 = 0;

        DENDRO_38 += DENDRO_77;
        DENDRO_38 += DENDRO_76;

        // initialize reduction for
        DENDRO_62 = 0;

        DENDRO_62 += DENDRO_40;
        DENDRO_62 += DENDRO_78;

        // initialize reduction for
        DENDRO_40 = 0;

        DENDRO_40 += DENDRO_80;
        DENDRO_40 += DENDRO_79;

        // initialize reduction for
        DENDRO_63 = 0;

        DENDRO_63 += DENDRO_82;
        DENDRO_63 += DENDRO_81;

        // initialize reduction for
        DENDRO_64 = 0;

        DENDRO_64 += DENDRO_49;
        DENDRO_64 += DENDRO_83;

        // initialize reduction for
        DENDRO_49 = 0;

        DENDRO_49 += DENDRO_52;
        DENDRO_49 += DENDRO_84;

        // initialize reduction for
        DENDRO_52 = 0;

        DENDRO_52 += DENDRO_53;
        DENDRO_52 += DENDRO_85;

        // initialize reduction for
        DENDRO_53 = 0;

        DENDRO_53 += DENDRO_88;
        DENDRO_53 += DENDRO_87;
        DENDRO_53 += DENDRO_86;
        DENDRO_65 = -2;

        // initialize reduction for
        DENDRO_66 = 0;

        DENDRO_66 += DENDRO_10;
        DENDRO_66 += DENDRO_89;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_16;
        DENDRO_10 += DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 0;

        DENDRO_13 += DENDRO_90;
        DENDRO_13 += DENDRO_17;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_51;
        DENDRO_16 *= DENDRO_12;
        DENDRO_16 *= DENDRO_34;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_51;
        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_35;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_51;
        DENDRO_34 *= DENDRO_24;
        DENDRO_34 *= DENDRO_36;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_91;
        DENDRO_35 *= DENDRO_39;
        DENDRO_35 *= DENDRO_20;
        DENDRO_35 *= DENDRO_7;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_91;
        DENDRO_20 *= DENDRO_41;
        DENDRO_20 *= DENDRO_21;
        DENDRO_20 *= DENDRO_9;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_91;
        DENDRO_21 *= DENDRO_48;
        DENDRO_21 *= DENDRO_24;
        DENDRO_21 *= DENDRO_11;

        // initialize reduction for
        DENDRO_24 = 1;

        DENDRO_24 *= DENDRO_23;
        DENDRO_24 *= DENDRO_4;
        DENDRO_24 *= DENDRO_7;
        DENDRO_24 *= DENDRO_92;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_23;
        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_9;
        DENDRO_4 *= DENDRO_42;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_23;
        DENDRO_6 *= DENDRO_12;
        DENDRO_6 *= DENDRO_11;
        DENDRO_6 *= DENDRO_43;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_23;
        DENDRO_12 *= DENDRO_14;
        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_26;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_23;
        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_9;
        DENDRO_14 *= DENDRO_25;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_23;
        DENDRO_15 *= DENDRO_18;
        DENDRO_15 *= DENDRO_11;
        DENDRO_15 *= DENDRO_27;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_44;
        DENDRO_18 *= DENDRO_0;
        DENDRO_18 *= DENDRO_28;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_44;
        DENDRO_25 *= DENDRO_3;
        DENDRO_25 *= DENDRO_45;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_44;
        DENDRO_26 *= DENDRO_5;
        DENDRO_26 *= DENDRO_46;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_93;
        DENDRO_27 *= DENDRO_7;
        DENDRO_27 *= DENDRO_47;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_93;
        DENDRO_28 *= DENDRO_9;
        DENDRO_28 *= DENDRO_94;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_93;
        DENDRO_36 *= DENDRO_11;
        DENDRO_36 *= DENDRO_95;

        // initialize reduction for
        DENDRO_42 = 1;

        DENDRO_42 *= DENDRO_98;
        DENDRO_42 *= DENDRO_97;
        DENDRO_42 *= DENDRO_96;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_98;
        DENDRO_43 *= DENDRO_100;
        DENDRO_43 *= DENDRO_99;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_98;
        DENDRO_44 *= DENDRO_102;
        DENDRO_44 *= DENDRO_101;

        // initialize reduction for
        DENDRO_45 = 1;

        DENDRO_45 *= DENDRO_98;
        DENDRO_45 *= DENDRO_39;
        DENDRO_45 *= DENDRO_54;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_98;
        DENDRO_39 *= DENDRO_41;
        DENDRO_39 *= DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_98;
        DENDRO_1 *= DENDRO_48;
        DENDRO_1 *= DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_23;
        DENDRO_8 *= DENDRO_0;
        DENDRO_8 *= DENDRO_103;

        // initialize reduction for
        DENDRO_41 = 1;

        DENDRO_41 *= DENDRO_23;
        DENDRO_41 *= DENDRO_0;
        DENDRO_41 *= DENDRO_29;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_23;
        DENDRO_29 *= DENDRO_0;
        DENDRO_29 *= DENDRO_30;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_23;
        DENDRO_30 *= DENDRO_0;
        DENDRO_30 *= DENDRO_31;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_23;
        DENDRO_31 *= DENDRO_0;
        DENDRO_31 *= DENDRO_32;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_23;
        DENDRO_32 *= DENDRO_0;
        DENDRO_32 *= DENDRO_33;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_23;
        DENDRO_33 *= DENDRO_3;
        DENDRO_33 *= DENDRO_55;

        // initialize reduction for
        DENDRO_46 = 1;

        DENDRO_46 *= DENDRO_23;
        DENDRO_46 *= DENDRO_3;
        DENDRO_46 *= DENDRO_56;

        // initialize reduction for
        DENDRO_47 = 1;

        DENDRO_47 *= DENDRO_23;
        DENDRO_47 *= DENDRO_3;
        DENDRO_47 *= DENDRO_37;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_23;
        DENDRO_37 *= DENDRO_3;
        DENDRO_37 *= DENDRO_57;

        // initialize reduction for
        DENDRO_48 = 1;

        DENDRO_48 *= DENDRO_23;
        DENDRO_48 *= DENDRO_3;
        DENDRO_48 *= DENDRO_58;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_23;
        DENDRO_54 *= DENDRO_3;
        DENDRO_54 *= DENDRO_59;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_23;
        DENDRO_55 *= DENDRO_5;
        DENDRO_55 *= DENDRO_60;

        // initialize reduction for
        DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_23;
        DENDRO_56 *= DENDRO_5;
        DENDRO_56 *= DENDRO_61;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_23;
        DENDRO_57 *= DENDRO_5;
        DENDRO_57 *= DENDRO_38;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_23;
        DENDRO_38 *= DENDRO_5;
        DENDRO_38 *= DENDRO_62;

        // initialize reduction for
        DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_23;
        DENDRO_58 *= DENDRO_5;
        DENDRO_58 *= DENDRO_40;

        // initialize reduction for
        DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_23;
        DENDRO_40 *= DENDRO_5;
        DENDRO_40 *= DENDRO_63;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_2;
        DENDRO_23 *= DENDRO_0;
        DENDRO_23 *= DENDRO_64;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_49;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_2;
        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_52;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_65;
        DENDRO_5 *= DENDRO_50;
        DENDRO_5 *= DENDRO_53;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_7;
        DENDRO_49 *= DENDRO_66;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_9;
        DENDRO_7 *= DENDRO_10;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_13;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_19;
        DENDRO_10 += DENDRO_34;
        DENDRO_10 += DENDRO_17;
        DENDRO_10 += DENDRO_16;

        // initialize reduction for
        DENDRO_11 = 0;

        DENDRO_11 += DENDRO_40;
        DENDRO_11 += DENDRO_58;
        DENDRO_11 += DENDRO_38;
        DENDRO_11 += DENDRO_57;
        DENDRO_11 += DENDRO_56;
        DENDRO_11 += DENDRO_55;
        DENDRO_11 += DENDRO_54;
        DENDRO_11 += DENDRO_48;
        DENDRO_11 += DENDRO_37;
        DENDRO_11 += DENDRO_47;
        DENDRO_11 += DENDRO_46;
        DENDRO_11 += DENDRO_33;
        DENDRO_11 += DENDRO_32;
        DENDRO_11 += DENDRO_31;
        DENDRO_11 += DENDRO_30;
        DENDRO_11 += DENDRO_29;
        DENDRO_11 += DENDRO_41;
        DENDRO_11 += DENDRO_8;
        DENDRO_11 += DENDRO_1;
        DENDRO_11 += DENDRO_39;
        DENDRO_11 += DENDRO_45;
        DENDRO_11 += DENDRO_44;
        DENDRO_11 += DENDRO_43;
        DENDRO_11 += DENDRO_42;
        DENDRO_11 += DENDRO_36;
        DENDRO_11 += DENDRO_28;
        DENDRO_11 += DENDRO_27;
        DENDRO_11 += DENDRO_26;
        DENDRO_11 += DENDRO_25;
        DENDRO_11 += DENDRO_18;
        DENDRO_11 += DENDRO_15;
        DENDRO_11 += DENDRO_14;
        DENDRO_11 += DENDRO_12;
        DENDRO_11 += DENDRO_6;
        DENDRO_11 += DENDRO_4;
        DENDRO_11 += DENDRO_24;
        DENDRO_11 += DENDRO_21;
        DENDRO_11 += DENDRO_20;
        DENDRO_11 += DENDRO_35;
        DENDRO_1 = DENDRO_50;

        DENDRO_1 *= DENDRO_50;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_9;
        DENDRO_4 += DENDRO_7;
        DENDRO_4 += DENDRO_49;
        DENDRO_4 += DENDRO_5;
        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_0;
        DENDRO_4 += DENDRO_23;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_50;
        DENDRO_0 *= DENDRO_10;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_11;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_101;
        DENDRO_1 *= DENDRO_4;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_51;
        DENDRO_3 *= DENDRO_22;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_2;
        DENDRO_4 += DENDRO_0;
        DENDRO_0 = DENDRO_50;

        DENDRO_0 *= DENDRO_50;
        DENDRO_0 = 1 / DENDRO_0;
        DENDRO_1 = 1.0 / 4.0;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_0;
        DENDRO_2 *= DENDRO_4;
        DENDRO_RIJ0 = DENDRO_2;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C2_k2_4;
        double DENDRO_2  = 2;
        double DENDRO_3  = DENDRO_igt2;
        double DENDRO_4  = DENDRO_C2_k2_2;
        double DENDRO_5  = DENDRO_igt1;
        double DENDRO_6  = DENDRO_C2_k2_1;
        double DENDRO_7  = DENDRO_igt5;
        double DENDRO_8  = DENDRO_C2_k2_5;
        double DENDRO_9  = DENDRO_igt3;
        double DENDRO_10 = DENDRO_C2_k2_3;
        double DENDRO_11 = DENDRO_igt0;
        double DENDRO_12 = DENDRO_C2_k2_0;
        double DENDRO_13 = DENDRO_C2_k1_4;
        double DENDRO_14 = DENDRO_C2_k1_2;
        double DENDRO_15 = DENDRO_C2_k1_1;
        double DENDRO_16 = DENDRO_C2_k1_5;
        double DENDRO_17 = DENDRO_C2_k1_3;
        double DENDRO_18 = DENDRO_C2_k1_0;
        double DENDRO_19 = DENDRO_C2_k0_4;
        double DENDRO_20 = DENDRO_C2_k0_2;
        double DENDRO_21 = DENDRO_C2_k0_1;
        double DENDRO_22 = DENDRO_C2_k0_5;
        double DENDRO_23 = DENDRO_C2_k0_3;
        double DENDRO_24 = DENDRO_C2_k0_0;

        // initialize reduction for
        double DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_1;
        DENDRO_25 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_4;
        DENDRO_26 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_6;
        DENDRO_27 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_8;
        DENDRO_28 *= DENDRO_7;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_10;
        DENDRO_8 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_12;
        DENDRO_29 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_2;
        DENDRO_30 *= DENDRO_13;
        DENDRO_30 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_2;
        DENDRO_31 *= DENDRO_14;
        DENDRO_31 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_15;
        DENDRO_32 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_16;
        DENDRO_33 *= DENDRO_7;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_18;
        DENDRO_34 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_2;
        DENDRO_35 *= DENDRO_19;
        DENDRO_35 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_2;
        DENDRO_36 *= DENDRO_20;
        DENDRO_36 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_2;
        DENDRO_37 *= DENDRO_21;
        DENDRO_37 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_22;
        DENDRO_38 *= DENDRO_7;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_23;
        DENDRO_22 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_24;
        DENDRO_39 *= DENDRO_11;
        double DENDRO_40 = grad_2_chi;
        double DENDRO_41 = grad_1_chi;
        double DENDRO_42 = grad_0_chi;
        double DENDRO_43 = DENDRO_C1_k2_4;
        double DENDRO_44 = DENDRO_C1_k1_5;
        double DENDRO_45 = DENDRO_C1_k0_5;
        double DENDRO_46 = DENDRO_C1_k1_2;
        double DENDRO_47 = DENDRO_C1_k0_4;
        double DENDRO_48 = DENDRO_C1_k0_2;
        double DENDRO_49 = DENDRO_C1_k1_4;
        double DENDRO_50 = DENDRO_C1_k2_3;
        double DENDRO_51 = DENDRO_C1_k1_3;
        double DENDRO_52 = DENDRO_C1_k0_3;
        double DENDRO_53 = DENDRO_C1_k1_1;
        double DENDRO_54 = DENDRO_C1_k0_1;
        double DENDRO_55 = DENDRO_C1_k2_1;
        double DENDRO_56 = DENDRO_C1_k1_0;
        double DENDRO_57 = DENDRO_C1_k0_0;
        double DENDRO_58 = grad2_1_2_chi;
        double DENDRO_59 = chi[pp];
        double DENDRO_60 = -3;
        double DENDRO_61 = grad2_0_2_chi;
        double DENDRO_62 = grad2_0_1_chi;

        // initialize reduction for
        double DENDRO_63 = 0;

        DENDRO_63 += DENDRO_29;
        DENDRO_63 += DENDRO_8;
        DENDRO_63 += DENDRO_28;
        DENDRO_63 += DENDRO_27;
        DENDRO_63 += DENDRO_26;
        DENDRO_63 += DENDRO_25;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_34;
        DENDRO_8 += DENDRO_16;
        DENDRO_8 += DENDRO_33;
        DENDRO_8 += DENDRO_32;
        DENDRO_8 += DENDRO_31;
        DENDRO_8 += DENDRO_30;

        // initialize reduction for
        DENDRO_16 = 0;

        DENDRO_16 += DENDRO_39;
        DENDRO_16 += DENDRO_22;
        DENDRO_16 += DENDRO_38;
        DENDRO_16 += DENDRO_37;
        DENDRO_16 += DENDRO_36;
        DENDRO_16 += DENDRO_35;
        DENDRO_22 = grad2_2_2_chi;
        DENDRO_25 = DENDRO_40;

        DENDRO_25 *= DENDRO_40;
        DENDRO_26 = grad2_1_1_chi;
        DENDRO_27 = DENDRO_41;

        DENDRO_27 *= DENDRO_41;
        DENDRO_28 = grad2_0_0_chi;
        DENDRO_29 = DENDRO_42;

        DENDRO_29 *= DENDRO_42;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_43;
        DENDRO_30 *= DENDRO_4;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_44;
        DENDRO_31 *= DENDRO_4;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_45;
        DENDRO_32 *= DENDRO_1;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_46;
        DENDRO_33 *= DENDRO_20;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_47;
        DENDRO_34 *= DENDRO_20;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_48;
        DENDRO_35 *= DENDRO_19;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_2;
        DENDRO_36 *= DENDRO_49;
        DENDRO_36 *= DENDRO_14;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_47;
        DENDRO_37 *= DENDRO_13;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_50;
        DENDRO_38 *= DENDRO_4;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_44;
        DENDRO_39 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_45;
        DENDRO_64 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_43;
        DENDRO_65 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_49;
        DENDRO_66 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_47;
        DENDRO_67 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_68 = 1;

        DENDRO_68 *= DENDRO_49;
        DENDRO_68 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_69 = 1;

        DENDRO_69 *= DENDRO_51;
        DENDRO_69 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_70 = 1;

        DENDRO_70 *= DENDRO_47;
        DENDRO_70 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_71 = 1;

        DENDRO_71 *= DENDRO_52;
        DENDRO_71 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_72 = 1;

        DENDRO_72 *= DENDRO_46;
        DENDRO_72 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_73 = 1;

        DENDRO_73 *= DENDRO_52;
        DENDRO_73 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_74 = 1;

        DENDRO_74 *= DENDRO_48;
        DENDRO_74 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_75 = 1;

        DENDRO_75 *= DENDRO_53;
        DENDRO_75 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_76 = 1;

        DENDRO_76 *= DENDRO_47;
        DENDRO_76 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_77 = 1;

        DENDRO_77 *= DENDRO_54;
        DENDRO_77 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_78 = 1;

        DENDRO_78 *= DENDRO_50;
        DENDRO_78 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_79 = 1;

        DENDRO_79 *= DENDRO_49;
        DENDRO_79 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_80 = 1;

        DENDRO_80 *= DENDRO_47;
        DENDRO_80 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_81 = 1;

        DENDRO_81 *= DENDRO_53;
        DENDRO_81 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_82 = 1;

        DENDRO_82 *= DENDRO_52;
        DENDRO_82 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_83 = 1;

        DENDRO_83 *= DENDRO_54;
        DENDRO_83 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_84 = 1;

        DENDRO_84 *= DENDRO_2;
        DENDRO_84 *= DENDRO_51;
        DENDRO_84 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_85 = 1;

        DENDRO_85 *= DENDRO_52;
        DENDRO_85 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_86 = 1;

        DENDRO_86 *= DENDRO_55;
        DENDRO_86 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_87 = 1;

        DENDRO_87 *= DENDRO_44;
        DENDRO_87 *= DENDRO_12;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_45;
        DENDRO_44 *= DENDRO_6;

        // initialize reduction for
        DENDRO_45 = 1;

        DENDRO_45 *= DENDRO_49;
        DENDRO_45 *= DENDRO_18;

        // initialize reduction for
        double DENDRO_88 = 1;

        DENDRO_88 *= DENDRO_53;
        DENDRO_88 *= DENDRO_14;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_47;
        DENDRO_14 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_89 = 1;

        DENDRO_89 *= DENDRO_43;
        DENDRO_89 *= DENDRO_12;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_46;
        DENDRO_43 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_48;
        DENDRO_4 *= DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_54;
        DENDRO_1 *= DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_46;
        DENDRO_13 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_90 = 1;

        DENDRO_90 *= DENDRO_48;
        DENDRO_90 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_91 = 1;

        DENDRO_91 *= DENDRO_54;
        DENDRO_91 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_92 = 1;

        DENDRO_92 *= DENDRO_56;
        DENDRO_92 *= DENDRO_20;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_47;
        DENDRO_20 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_93 = 1;

        DENDRO_93 *= DENDRO_57;
        DENDRO_93 *= DENDRO_19;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_55;
        DENDRO_19 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_94 = 1;

        DENDRO_94 *= DENDRO_49;
        DENDRO_94 *= DENDRO_12;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_47;
        DENDRO_49 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_95 = 1;

        DENDRO_95 *= DENDRO_51;
        DENDRO_95 *= DENDRO_18;

        // initialize reduction for
        DENDRO_51 = 1;

        DENDRO_51 *= DENDRO_53;
        DENDRO_51 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_96 = 1;

        DENDRO_96 *= DENDRO_52;
        DENDRO_96 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_97 = 1;

        DENDRO_97 *= DENDRO_50;
        DENDRO_97 *= DENDRO_12;

        // initialize reduction for
        DENDRO_50 = 1;

        DENDRO_50 *= DENDRO_46;
        DENDRO_50 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_98 = 1;

        DENDRO_98 *= DENDRO_48;
        DENDRO_98 *= DENDRO_10;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_54;
        DENDRO_10 *= DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_56;
        DENDRO_17 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_99 = 1;

        DENDRO_99 *= DENDRO_52;
        DENDRO_99 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_100 = 1;

        DENDRO_100 *= DENDRO_57;
        DENDRO_100 *= DENDRO_23;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_2;
        DENDRO_23 *= DENDRO_54;
        DENDRO_23 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_101 = 1;

        DENDRO_101 *= DENDRO_53;
        DENDRO_101 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_102 = 1;

        DENDRO_102 *= DENDRO_55;
        DENDRO_102 *= DENDRO_12;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_46;
        DENDRO_55 *= DENDRO_12;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_48;
        DENDRO_12 *= DENDRO_6;

        // initialize reduction for
        DENDRO_48 = 1;

        DENDRO_48 *= DENDRO_56;
        DENDRO_48 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_103 = 1;

        DENDRO_103 *= DENDRO_54;
        DENDRO_103 *= DENDRO_24;

        // initialize reduction for
        DENDRO_24 = 1;

        DENDRO_24 *= DENDRO_57;
        DENDRO_24 *= DENDRO_21;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_2;
        DENDRO_57 *= DENDRO_53;
        DENDRO_57 *= DENDRO_18;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_54;
        DENDRO_18 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_104 = 1;

        DENDRO_104 *= DENDRO_2;
        DENDRO_104 *= DENDRO_59;
        DENDRO_104 *= DENDRO_58;

        // initialize reduction for
        DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_60;
        DENDRO_58 *= DENDRO_41;
        DENDRO_58 *= DENDRO_40;

        // initialize reduction for
        double DENDRO_105 = 1;

        DENDRO_105 *= DENDRO_2;
        DENDRO_105 *= DENDRO_59;
        DENDRO_105 *= DENDRO_61;

        // initialize reduction for
        DENDRO_61 = 1;

        DENDRO_61 *= DENDRO_60;
        DENDRO_61 *= DENDRO_42;
        DENDRO_61 *= DENDRO_40;

        // initialize reduction for
        double DENDRO_106 = 1;

        DENDRO_106 *= DENDRO_2;
        DENDRO_106 *= DENDRO_59;
        DENDRO_106 *= DENDRO_62;

        // initialize reduction for
        double DENDRO_107 = 1;

        DENDRO_107 *= DENDRO_60;
        DENDRO_107 *= DENDRO_42;
        DENDRO_107 *= DENDRO_41;

        // initialize reduction for
        double DENDRO_108 = 1;

        DENDRO_108 *= DENDRO_40;
        DENDRO_108 *= DENDRO_63;

        // initialize reduction for
        double DENDRO_109 = 1;

        DENDRO_109 *= DENDRO_41;
        DENDRO_109 *= DENDRO_8;

        // initialize reduction for
        double DENDRO_110 = 1;

        DENDRO_110 *= DENDRO_42;
        DENDRO_110 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_111 = 1;

        DENDRO_111 *= DENDRO_2;
        DENDRO_111 *= DENDRO_59;
        DENDRO_111 *= DENDRO_22;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_60;
        DENDRO_22 *= DENDRO_25;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_59;
        DENDRO_25 *= DENDRO_26;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_60;
        DENDRO_26 *= DENDRO_27;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_59;
        DENDRO_27 *= DENDRO_28;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_60;
        DENDRO_28 *= DENDRO_29;
        DENDRO_29         = -1;
        DENDRO_60         = grad2_1_2_gt1;
        double DENDRO_112 = -4.0;
        double DENDRO_113 = grad2_0_2_gt1;
        double DENDRO_114 = grad2_0_1_gt1;
        double DENDRO_115 = grad2_2_2_gt1;
        double DENDRO_116 = -2.0;
        double DENDRO_117 = grad2_1_1_gt1;
        double DENDRO_118 = grad2_0_0_gt1;

        // initialize reduction for
        double DENDRO_119 = 0;

        DENDRO_119 += DENDRO_47;
        DENDRO_119 += DENDRO_46;
        DENDRO_46 = 2.0;

        // initialize reduction for
        DENDRO_47 = 0;

        DENDRO_47 += DENDRO_52;
        DENDRO_47 += DENDRO_53;

        // initialize reduction for
        DENDRO_52 = 0;

        DENDRO_52 += DENDRO_54;
        DENDRO_52 += DENDRO_56;
        DENDRO_53         = gt2[pp];
        DENDRO_54         = grad_1_Gt2;
        DENDRO_56         = gt1[pp];
        double DENDRO_120 = grad_1_Gt1;
        double DENDRO_121 = gt0[pp];
        double DENDRO_122 = grad_1_Gt0;
        double DENDRO_123 = gt4[pp];
        double DENDRO_124 = grad_0_Gt2;
        double DENDRO_125 = gt3[pp];
        double DENDRO_126 = grad_0_Gt1;
        double DENDRO_127 = grad_0_Gt0;

        // initialize reduction for
        double DENDRO_128 = 0;

        DENDRO_128 += DENDRO_32;
        DENDRO_128 += DENDRO_31;
        DENDRO_128 += DENDRO_30;
        DENDRO_30 = 4;

        // initialize reduction for
        DENDRO_31 = 0;

        DENDRO_31 += DENDRO_35;
        DENDRO_31 += DENDRO_34;
        DENDRO_31 += DENDRO_33;

        // initialize reduction for
        DENDRO_32 = 0;

        DENDRO_32 += DENDRO_37;
        DENDRO_32 += DENDRO_36;

        // initialize reduction for
        DENDRO_33 = 0;

        DENDRO_33 += DENDRO_64;
        DENDRO_33 += DENDRO_39;
        DENDRO_33 += DENDRO_38;

        // initialize reduction for
        DENDRO_34 = 0;

        DENDRO_34 += DENDRO_67;
        DENDRO_34 += DENDRO_66;
        DENDRO_34 += DENDRO_65;

        // initialize reduction for
        DENDRO_35 = 0;

        DENDRO_35 += DENDRO_70;
        DENDRO_35 += DENDRO_69;
        DENDRO_35 += DENDRO_68;

        // initialize reduction for
        DENDRO_36 = 0;

        DENDRO_36 += DENDRO_71;
        DENDRO_36 += DENDRO_69;
        DENDRO_36 += DENDRO_68;

        // initialize reduction for
        DENDRO_37 = 0;

        DENDRO_37 += DENDRO_74;
        DENDRO_37 += DENDRO_73;
        DENDRO_37 += DENDRO_72;

        // initialize reduction for
        DENDRO_38 = 0;

        DENDRO_38 += DENDRO_77;
        DENDRO_38 += DENDRO_76;
        DENDRO_38 += DENDRO_75;

        // initialize reduction for
        DENDRO_39 = 0;

        DENDRO_39 += DENDRO_80;
        DENDRO_39 += DENDRO_79;
        DENDRO_39 += DENDRO_78;

        // initialize reduction for
        DENDRO_64 = 0;

        DENDRO_64 += DENDRO_83;
        DENDRO_64 += DENDRO_82;
        DENDRO_64 += DENDRO_81;

        // initialize reduction for
        DENDRO_65 = 0;

        DENDRO_65 += DENDRO_85;
        DENDRO_65 += DENDRO_84;

        // initialize reduction for
        DENDRO_66 = 0;

        DENDRO_66 += DENDRO_44;
        DENDRO_66 += DENDRO_87;
        DENDRO_66 += DENDRO_86;

        // initialize reduction for
        DENDRO_44 = 0;

        DENDRO_44 += DENDRO_14;
        DENDRO_44 += DENDRO_88;
        DENDRO_44 += DENDRO_45;

        // initialize reduction for
        DENDRO_14 = 0;

        DENDRO_14 += DENDRO_4;
        DENDRO_14 += DENDRO_43;
        DENDRO_14 += DENDRO_89;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_88;
        DENDRO_4 += DENDRO_45;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_91;
        DENDRO_1 += DENDRO_90;
        DENDRO_1 += DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 0;

        DENDRO_13 += DENDRO_93;
        DENDRO_13 += DENDRO_20;
        DENDRO_13 += DENDRO_92;

        // initialize reduction for
        DENDRO_20 = 0;

        DENDRO_20 += DENDRO_49;
        DENDRO_20 += DENDRO_94;
        DENDRO_20 += DENDRO_19;

        // initialize reduction for
        DENDRO_19 = 0;

        DENDRO_19 += DENDRO_96;
        DENDRO_19 += DENDRO_51;
        DENDRO_19 += DENDRO_95;

        // initialize reduction for
        DENDRO_43 = 0;

        DENDRO_43 += DENDRO_98;
        DENDRO_43 += DENDRO_50;
        DENDRO_43 += DENDRO_97;

        // initialize reduction for
        DENDRO_45 = 0;

        DENDRO_45 += DENDRO_10;
        DENDRO_45 += DENDRO_51;
        DENDRO_45 += DENDRO_95;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_100;
        DENDRO_10 += DENDRO_99;
        DENDRO_10 += DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 0;

        DENDRO_17 += DENDRO_101;
        DENDRO_17 += DENDRO_23;

        // initialize reduction for
        DENDRO_23 = 0;

        DENDRO_23 += DENDRO_12;
        DENDRO_23 += DENDRO_55;
        DENDRO_23 += DENDRO_102;

        // initialize reduction for
        DENDRO_12 = 0;

        DENDRO_12 += DENDRO_24;
        DENDRO_12 += DENDRO_103;
        DENDRO_12 += DENDRO_48;

        // initialize reduction for
        DENDRO_24 = 0;

        DENDRO_24 += DENDRO_18;
        DENDRO_24 += DENDRO_57;

        // initialize reduction for
        DENDRO_18 = 0;

        DENDRO_18 += DENDRO_58;
        DENDRO_18 += DENDRO_104;

        // initialize reduction for
        DENDRO_48 = 0;

        DENDRO_48 += DENDRO_61;
        DENDRO_48 += DENDRO_105;

        // initialize reduction for
        DENDRO_49 = 0;

        DENDRO_49 += DENDRO_107;
        DENDRO_49 += DENDRO_106;

        // initialize reduction for
        DENDRO_50 = 0;

        DENDRO_50 += DENDRO_110;
        DENDRO_50 += DENDRO_109;
        DENDRO_50 += DENDRO_108;
        DENDRO_51 = -2;

        // initialize reduction for
        DENDRO_55 = 0;

        DENDRO_55 += DENDRO_22;
        DENDRO_55 += DENDRO_111;

        // initialize reduction for
        DENDRO_22 = 0;

        DENDRO_22 += DENDRO_26;
        DENDRO_22 += DENDRO_25;

        // initialize reduction for
        DENDRO_25 = 0;

        DENDRO_25 += DENDRO_28;
        DENDRO_25 += DENDRO_27;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_29;
        DENDRO_26 *= DENDRO_6;
        DENDRO_26 *= DENDRO_40;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_29;
        DENDRO_6 *= DENDRO_15;
        DENDRO_6 *= DENDRO_41;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_29;
        DENDRO_15 *= DENDRO_21;
        DENDRO_15 *= DENDRO_42;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_112;
        DENDRO_21 *= DENDRO_0;
        DENDRO_21 *= DENDRO_60;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_112;
        DENDRO_27 *= DENDRO_3;
        DENDRO_27 *= DENDRO_113;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_112;
        DENDRO_28 *= DENDRO_5;
        DENDRO_28 *= DENDRO_114;

        // initialize reduction for
        DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_116;
        DENDRO_40 *= DENDRO_7;
        DENDRO_40 *= DENDRO_115;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_116;
        DENDRO_57 *= DENDRO_9;
        DENDRO_57 *= DENDRO_117;

        // initialize reduction for
        DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_116;
        DENDRO_58 *= DENDRO_11;
        DENDRO_58 *= DENDRO_118;

        // initialize reduction for
        DENDRO_60 = 1;

        DENDRO_60 *= DENDRO_46;
        DENDRO_60 *= DENDRO_119;
        DENDRO_60 *= DENDRO_63;

        // initialize reduction for
        DENDRO_61 = 1;

        DENDRO_61 *= DENDRO_46;
        DENDRO_61 *= DENDRO_47;
        DENDRO_61 *= DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_46;
        DENDRO_8 *= DENDRO_52;
        DENDRO_8 *= DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_46;
        DENDRO_16 *= DENDRO_54;
        DENDRO_16 *= DENDRO_53;

        // initialize reduction for
        DENDRO_47 = 1;

        DENDRO_47 *= DENDRO_46;
        DENDRO_47 *= DENDRO_120;
        DENDRO_47 *= DENDRO_56;

        // initialize reduction for
        DENDRO_52 = 1;

        DENDRO_52 *= DENDRO_46;
        DENDRO_52 *= DENDRO_122;
        DENDRO_52 *= DENDRO_121;

        // initialize reduction for
        DENDRO_53 = 1;

        DENDRO_53 *= DENDRO_46;
        DENDRO_53 *= DENDRO_124;
        DENDRO_53 *= DENDRO_123;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_46;
        DENDRO_54 *= DENDRO_126;
        DENDRO_54 *= DENDRO_125;

        // initialize reduction for
        DENDRO_63 = 1;

        DENDRO_63 *= DENDRO_46;
        DENDRO_63 *= DENDRO_127;
        DENDRO_63 *= DENDRO_56;

        // initialize reduction for
        DENDRO_46 = 1;

        DENDRO_46 *= DENDRO_30;
        DENDRO_46 *= DENDRO_7;
        DENDRO_46 *= DENDRO_128;

        // initialize reduction for
        DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_30;
        DENDRO_67 *= DENDRO_7;
        DENDRO_67 *= DENDRO_31;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_30;
        DENDRO_31 *= DENDRO_7;
        DENDRO_31 *= DENDRO_32;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_30;
        DENDRO_32 *= DENDRO_0;
        DENDRO_32 *= DENDRO_33;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_30;
        DENDRO_33 *= DENDRO_0;
        DENDRO_33 *= DENDRO_34;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_30;
        DENDRO_34 *= DENDRO_0;
        DENDRO_34 *= DENDRO_35;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_30;
        DENDRO_35 *= DENDRO_0;
        DENDRO_35 *= DENDRO_36;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_30;
        DENDRO_36 *= DENDRO_0;
        DENDRO_36 *= DENDRO_37;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_30;
        DENDRO_37 *= DENDRO_0;
        DENDRO_37 *= DENDRO_38;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_30;
        DENDRO_38 *= DENDRO_9;
        DENDRO_38 *= DENDRO_39;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_30;
        DENDRO_39 *= DENDRO_9;
        DENDRO_39 *= DENDRO_64;

        // initialize reduction for
        DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_30;
        DENDRO_64 *= DENDRO_9;
        DENDRO_64 *= DENDRO_65;

        // initialize reduction for
        DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_30;
        DENDRO_65 *= DENDRO_3;
        DENDRO_65 *= DENDRO_66;

        // initialize reduction for
        DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_30;
        DENDRO_66 *= DENDRO_3;
        DENDRO_66 *= DENDRO_44;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_30;
        DENDRO_44 *= DENDRO_3;
        DENDRO_44 *= DENDRO_14;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_30;
        DENDRO_14 *= DENDRO_3;
        DENDRO_14 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_30;
        DENDRO_4 *= DENDRO_3;
        DENDRO_4 *= DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_30;
        DENDRO_1 *= DENDRO_3;
        DENDRO_1 *= DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_30;
        DENDRO_13 *= DENDRO_5;
        DENDRO_13 *= DENDRO_20;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_30;
        DENDRO_20 *= DENDRO_5;
        DENDRO_20 *= DENDRO_19;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_30;
        DENDRO_19 *= DENDRO_5;
        DENDRO_19 *= DENDRO_43;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_30;
        DENDRO_43 *= DENDRO_5;
        DENDRO_43 *= DENDRO_45;

        // initialize reduction for
        DENDRO_45 = 1;

        DENDRO_45 *= DENDRO_30;
        DENDRO_45 *= DENDRO_5;
        DENDRO_45 *= DENDRO_10;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_30;
        DENDRO_10 *= DENDRO_5;
        DENDRO_10 *= DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_30;
        DENDRO_17 *= DENDRO_11;
        DENDRO_17 *= DENDRO_23;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_30;
        DENDRO_23 *= DENDRO_11;
        DENDRO_23 *= DENDRO_12;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_30;
        DENDRO_12 *= DENDRO_11;
        DENDRO_12 *= DENDRO_24;

        // initialize reduction for
        DENDRO_24 = 1;

        DENDRO_24 *= DENDRO_2;
        DENDRO_24 *= DENDRO_0;
        DENDRO_24 *= DENDRO_18;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_48;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_2;
        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_49;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_51;
        DENDRO_5 *= DENDRO_59;
        DENDRO_5 *= DENDRO_50;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_7;
        DENDRO_18 *= DENDRO_55;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_9;
        DENDRO_7 *= DENDRO_22;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_25;

        // initialize reduction for
        DENDRO_11 = 0;

        DENDRO_11 += DENDRO_62;
        DENDRO_11 += DENDRO_15;
        DENDRO_11 += DENDRO_6;
        DENDRO_11 += DENDRO_26;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_12;
        DENDRO_6 += DENDRO_23;
        DENDRO_6 += DENDRO_17;
        DENDRO_6 += DENDRO_10;
        DENDRO_6 += DENDRO_45;
        DENDRO_6 += DENDRO_43;
        DENDRO_6 += DENDRO_19;
        DENDRO_6 += DENDRO_20;
        DENDRO_6 += DENDRO_13;
        DENDRO_6 += DENDRO_1;
        DENDRO_6 += DENDRO_4;
        DENDRO_6 += DENDRO_14;
        DENDRO_6 += DENDRO_44;
        DENDRO_6 += DENDRO_66;
        DENDRO_6 += DENDRO_65;
        DENDRO_6 += DENDRO_64;
        DENDRO_6 += DENDRO_39;
        DENDRO_6 += DENDRO_38;
        DENDRO_6 += DENDRO_37;
        DENDRO_6 += DENDRO_36;
        DENDRO_6 += DENDRO_35;
        DENDRO_6 += DENDRO_34;
        DENDRO_6 += DENDRO_33;
        DENDRO_6 += DENDRO_32;
        DENDRO_6 += DENDRO_31;
        DENDRO_6 += DENDRO_67;
        DENDRO_6 += DENDRO_46;
        DENDRO_6 += DENDRO_63;
        DENDRO_6 += DENDRO_54;
        DENDRO_6 += DENDRO_53;
        DENDRO_6 += DENDRO_52;
        DENDRO_6 += DENDRO_47;
        DENDRO_6 += DENDRO_16;
        DENDRO_6 += DENDRO_8;
        DENDRO_6 += DENDRO_61;
        DENDRO_6 += DENDRO_60;
        DENDRO_6 += DENDRO_58;
        DENDRO_6 += DENDRO_57;
        DENDRO_6 += DENDRO_40;
        DENDRO_6 += DENDRO_28;
        DENDRO_6 += DENDRO_27;
        DENDRO_6 += DENDRO_21;
        DENDRO_1 = DENDRO_59;

        DENDRO_1 *= DENDRO_59;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_9;
        DENDRO_4 += DENDRO_7;
        DENDRO_4 += DENDRO_18;
        DENDRO_4 += DENDRO_5;
        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_0;
        DENDRO_4 += DENDRO_24;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_59;
        DENDRO_0 *= DENDRO_11;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_29;
        DENDRO_2 *= DENDRO_42;
        DENDRO_2 *= DENDRO_41;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_1;
        DENDRO_3 *= DENDRO_6;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_56;
        DENDRO_1 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_2;
        DENDRO_4 += DENDRO_0;
        DENDRO_0 = DENDRO_59;

        DENDRO_0 *= DENDRO_59;
        DENDRO_0 = 1 / DENDRO_0;
        DENDRO_1 = 1.0 / 4.0;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_0;
        DENDRO_2 *= DENDRO_4;
        DENDRO_RIJ1 = DENDRO_2;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C2_k2_4;
        double DENDRO_2  = 2;
        double DENDRO_3  = DENDRO_igt2;
        double DENDRO_4  = DENDRO_C2_k2_2;
        double DENDRO_5  = DENDRO_igt1;
        double DENDRO_6  = DENDRO_C2_k2_1;
        double DENDRO_7  = DENDRO_igt5;
        double DENDRO_8  = DENDRO_C2_k2_5;
        double DENDRO_9  = DENDRO_igt3;
        double DENDRO_10 = DENDRO_C2_k2_3;
        double DENDRO_11 = DENDRO_igt0;
        double DENDRO_12 = DENDRO_C2_k2_0;
        double DENDRO_13 = DENDRO_C2_k1_4;
        double DENDRO_14 = DENDRO_C2_k1_2;
        double DENDRO_15 = DENDRO_C2_k1_1;
        double DENDRO_16 = DENDRO_C2_k1_5;
        double DENDRO_17 = DENDRO_C2_k1_3;
        double DENDRO_18 = DENDRO_C2_k1_0;
        double DENDRO_19 = DENDRO_C2_k0_4;
        double DENDRO_20 = DENDRO_C2_k0_2;
        double DENDRO_21 = DENDRO_C2_k0_1;
        double DENDRO_22 = DENDRO_C2_k0_5;
        double DENDRO_23 = DENDRO_C2_k0_3;
        double DENDRO_24 = DENDRO_C2_k0_0;

        // initialize reduction for
        double DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_1;
        DENDRO_25 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_4;
        DENDRO_26 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_6;
        DENDRO_27 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_8;
        DENDRO_28 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_10;
        DENDRO_29 *= DENDRO_9;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_12;
        DENDRO_10 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_2;
        DENDRO_30 *= DENDRO_13;
        DENDRO_30 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_2;
        DENDRO_31 *= DENDRO_14;
        DENDRO_31 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_15;
        DENDRO_32 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_16;
        DENDRO_33 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_17;
        DENDRO_34 *= DENDRO_9;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_11;

        // initialize reduction for
        double DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_2;
        DENDRO_35 *= DENDRO_19;
        DENDRO_35 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_2;
        DENDRO_36 *= DENDRO_20;
        DENDRO_36 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_2;
        DENDRO_37 *= DENDRO_21;
        DENDRO_37 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_22;
        DENDRO_38 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_23;
        DENDRO_39 *= DENDRO_9;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_24;
        DENDRO_23 *= DENDRO_11;
        double DENDRO_40 = grad_2_chi;
        double DENDRO_41 = grad_1_chi;
        double DENDRO_42 = grad_0_chi;
        double DENDRO_43 = DENDRO_C1_k2_4;
        double DENDRO_44 = DENDRO_C1_k1_5;
        double DENDRO_45 = DENDRO_C1_k0_4;
        double DENDRO_46 = DENDRO_C1_k2_2;
        double DENDRO_47 = DENDRO_C1_k0_5;
        double DENDRO_48 = DENDRO_C1_k0_2;
        double DENDRO_49 = DENDRO_C1_k2_5;
        double DENDRO_50 = DENDRO_C1_k1_4;
        double DENDRO_51 = DENDRO_C1_k2_3;
        double DENDRO_52 = DENDRO_C1_k0_3;
        double DENDRO_53 = DENDRO_C1_k2_1;
        double DENDRO_54 = DENDRO_C1_k0_1;
        double DENDRO_55 = DENDRO_C1_k1_2;
        double DENDRO_56 = DENDRO_C1_k2_0;
        double DENDRO_57 = DENDRO_C1_k0_0;
        double DENDRO_58 = grad2_1_2_chi;
        double DENDRO_59 = chi[pp];
        double DENDRO_60 = -3;
        double DENDRO_61 = grad2_0_2_chi;
        double DENDRO_62 = grad2_0_1_chi;

        // initialize reduction for
        double DENDRO_63 = 0;

        DENDRO_63 += DENDRO_10;
        DENDRO_63 += DENDRO_29;
        DENDRO_63 += DENDRO_28;
        DENDRO_63 += DENDRO_27;
        DENDRO_63 += DENDRO_26;
        DENDRO_63 += DENDRO_25;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_17;
        DENDRO_10 += DENDRO_34;
        DENDRO_10 += DENDRO_33;
        DENDRO_10 += DENDRO_32;
        DENDRO_10 += DENDRO_31;
        DENDRO_10 += DENDRO_30;

        // initialize reduction for
        DENDRO_17 = 0;

        DENDRO_17 += DENDRO_23;
        DENDRO_17 += DENDRO_39;
        DENDRO_17 += DENDRO_38;
        DENDRO_17 += DENDRO_37;
        DENDRO_17 += DENDRO_36;
        DENDRO_17 += DENDRO_35;
        DENDRO_23 = grad2_2_2_chi;
        DENDRO_25 = DENDRO_40;

        DENDRO_25 *= DENDRO_40;
        DENDRO_26 = grad2_1_1_chi;
        DENDRO_27 = DENDRO_41;

        DENDRO_27 *= DENDRO_41;
        DENDRO_28 = grad2_0_0_chi;
        DENDRO_29 = DENDRO_42;

        DENDRO_29 *= DENDRO_42;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_43;
        DENDRO_30 *= DENDRO_14;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_44;
        DENDRO_31 *= DENDRO_14;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_45;
        DENDRO_32 *= DENDRO_16;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_46;
        DENDRO_33 *= DENDRO_20;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_47;
        DENDRO_34 *= DENDRO_20;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_48;
        DENDRO_35 *= DENDRO_22;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_2;
        DENDRO_36 *= DENDRO_49;
        DENDRO_36 *= DENDRO_4;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_47;
        DENDRO_37 *= DENDRO_8;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_49;
        DENDRO_38 *= DENDRO_6;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_43;
        DENDRO_39 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_47;
        DENDRO_64 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_45;
        DENDRO_65 *= DENDRO_8;

        // initialize reduction for
        double DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_43;
        DENDRO_66 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_50;
        DENDRO_67 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_68 = 1;

        DENDRO_68 *= DENDRO_45;
        DENDRO_68 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_69 = 1;

        DENDRO_69 *= DENDRO_51;
        DENDRO_69 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_70 = 1;

        DENDRO_70 *= DENDRO_44;
        DENDRO_70 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_71 = 1;

        DENDRO_71 *= DENDRO_52;
        DENDRO_71 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_72 = 1;

        DENDRO_72 *= DENDRO_46;
        DENDRO_72 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_73 = 1;

        DENDRO_73 *= DENDRO_45;
        DENDRO_73 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_74 = 1;

        DENDRO_74 *= DENDRO_48;
        DENDRO_74 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_75 = 1;

        DENDRO_75 *= DENDRO_53;
        DENDRO_75 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_76 = 1;

        DENDRO_76 *= DENDRO_47;
        DENDRO_76 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_77 = 1;

        DENDRO_77 *= DENDRO_54;
        DENDRO_77 *= DENDRO_22;

        // initialize reduction for
        double DENDRO_78 = 1;

        DENDRO_78 *= DENDRO_51;
        DENDRO_78 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_79 = 1;

        DENDRO_79 *= DENDRO_50;
        DENDRO_79 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_80 = 1;

        DENDRO_80 *= DENDRO_52;
        DENDRO_80 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_81 = 1;

        DENDRO_81 *= DENDRO_53;
        DENDRO_81 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_82 = 1;

        DENDRO_82 *= DENDRO_45;
        DENDRO_82 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_83 = 1;

        DENDRO_83 *= DENDRO_54;
        DENDRO_83 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_84 = 1;

        DENDRO_84 *= DENDRO_2;
        DENDRO_84 *= DENDRO_43;
        DENDRO_84 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_85 = 1;

        DENDRO_85 *= DENDRO_45;
        DENDRO_85 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_86 = 1;

        DENDRO_86 *= DENDRO_49;
        DENDRO_86 *= DENDRO_12;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_46;
        DENDRO_49 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_87 = 1;

        DENDRO_87 *= DENDRO_47;
        DENDRO_87 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_88 = 1;

        DENDRO_88 *= DENDRO_43;
        DENDRO_88 *= DENDRO_18;

        // initialize reduction for
        double DENDRO_89 = 1;

        DENDRO_89 *= DENDRO_55;
        DENDRO_89 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_90 = 1;

        DENDRO_90 *= DENDRO_45;
        DENDRO_90 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_91 = 1;

        DENDRO_91 *= DENDRO_48;
        DENDRO_91 *= DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_53;
        DENDRO_8 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_92 = 1;

        DENDRO_92 *= DENDRO_44;
        DENDRO_92 *= DENDRO_18;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_54;
        DENDRO_44 *= DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_56;
        DENDRO_16 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_93 = 1;

        DENDRO_93 *= DENDRO_47;
        DENDRO_93 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_94 = 1;

        DENDRO_94 *= DENDRO_57;
        DENDRO_94 *= DENDRO_22;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_2;
        DENDRO_22 *= DENDRO_48;
        DENDRO_22 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_95 = 1;

        DENDRO_95 *= DENDRO_46;
        DENDRO_95 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_96 = 1;

        DENDRO_96 *= DENDRO_43;
        DENDRO_96 *= DENDRO_12;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_46;
        DENDRO_43 *= DENDRO_6;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_45;
        DENDRO_6 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_97 = 1;

        DENDRO_97 *= DENDRO_51;
        DENDRO_97 *= DENDRO_18;

        // initialize reduction for
        DENDRO_51 = 1;

        DENDRO_51 *= DENDRO_55;
        DENDRO_51 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_98 = 1;

        DENDRO_98 *= DENDRO_52;
        DENDRO_98 *= DENDRO_14;

        // initialize reduction for
        DENDRO_52 = 1;

        DENDRO_52 *= DENDRO_48;
        DENDRO_52 *= DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_53;
        DENDRO_1 *= DENDRO_15;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_50;
        DENDRO_15 *= DENDRO_18;

        // initialize reduction for
        DENDRO_50 = 1;

        DENDRO_50 *= DENDRO_54;
        DENDRO_50 *= DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_53;
        DENDRO_13 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_99 = 1;

        DENDRO_99 *= DENDRO_48;
        DENDRO_99 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_100 = 1;

        DENDRO_100 *= DENDRO_54;
        DENDRO_100 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_101 = 1;

        DENDRO_101 *= DENDRO_56;
        DENDRO_101 *= DENDRO_21;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_45;
        DENDRO_21 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_102 = 1;

        DENDRO_102 *= DENDRO_57;
        DENDRO_102 *= DENDRO_19;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_53;
        DENDRO_19 *= DENDRO_18;

        // initialize reduction for
        double DENDRO_103 = 1;

        DENDRO_103 *= DENDRO_55;
        DENDRO_103 *= DENDRO_18;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_54;
        DENDRO_18 *= DENDRO_14;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_56;
        DENDRO_54 *= DENDRO_24;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_48;
        DENDRO_55 *= DENDRO_24;

        // initialize reduction for
        DENDRO_24 = 1;

        DENDRO_24 *= DENDRO_57;
        DENDRO_24 *= DENDRO_20;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_2;
        DENDRO_57 *= DENDRO_46;
        DENDRO_57 *= DENDRO_12;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_48;
        DENDRO_12 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_104 = 1;

        DENDRO_104 *= DENDRO_2;
        DENDRO_104 *= DENDRO_59;
        DENDRO_104 *= DENDRO_58;

        // initialize reduction for
        DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_60;
        DENDRO_58 *= DENDRO_41;
        DENDRO_58 *= DENDRO_40;

        // initialize reduction for
        double DENDRO_105 = 1;

        DENDRO_105 *= DENDRO_2;
        DENDRO_105 *= DENDRO_59;
        DENDRO_105 *= DENDRO_61;

        // initialize reduction for
        double DENDRO_106 = 1;

        DENDRO_106 *= DENDRO_60;
        DENDRO_106 *= DENDRO_42;
        DENDRO_106 *= DENDRO_40;

        // initialize reduction for
        double DENDRO_107 = 1;

        DENDRO_107 *= DENDRO_2;
        DENDRO_107 *= DENDRO_59;
        DENDRO_107 *= DENDRO_62;

        // initialize reduction for
        DENDRO_62 = 1;

        DENDRO_62 *= DENDRO_60;
        DENDRO_62 *= DENDRO_42;
        DENDRO_62 *= DENDRO_41;

        // initialize reduction for
        double DENDRO_108 = 1;

        DENDRO_108 *= DENDRO_40;
        DENDRO_108 *= DENDRO_63;

        // initialize reduction for
        double DENDRO_109 = 1;

        DENDRO_109 *= DENDRO_41;
        DENDRO_109 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_110 = 1;

        DENDRO_110 *= DENDRO_42;
        DENDRO_110 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_111 = 1;

        DENDRO_111 *= DENDRO_2;
        DENDRO_111 *= DENDRO_59;
        DENDRO_111 *= DENDRO_23;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_60;
        DENDRO_23 *= DENDRO_25;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_59;
        DENDRO_25 *= DENDRO_26;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_60;
        DENDRO_26 *= DENDRO_27;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_59;
        DENDRO_27 *= DENDRO_28;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_60;
        DENDRO_28 *= DENDRO_29;
        DENDRO_29         = -1;
        DENDRO_60         = grad2_1_2_gt2;
        double DENDRO_112 = -4.0;
        double DENDRO_113 = grad2_0_2_gt2;
        double DENDRO_114 = grad2_0_1_gt2;
        double DENDRO_115 = grad2_2_2_gt2;
        double DENDRO_116 = -2.0;
        double DENDRO_117 = grad2_1_1_gt2;
        double DENDRO_118 = grad2_0_0_gt2;

        // initialize reduction for
        double DENDRO_119 = 0;

        DENDRO_119 += DENDRO_47;
        DENDRO_119 += DENDRO_46;
        DENDRO_46 = 2.0;

        // initialize reduction for
        DENDRO_47 = 0;

        DENDRO_47 += DENDRO_45;
        DENDRO_47 += DENDRO_53;

        // initialize reduction for
        DENDRO_45 = 0;

        DENDRO_45 += DENDRO_48;
        DENDRO_45 += DENDRO_56;
        DENDRO_48         = gt2[pp];
        DENDRO_53         = grad_2_Gt2;
        DENDRO_56         = gt1[pp];
        double DENDRO_120 = grad_2_Gt1;
        double DENDRO_121 = gt0[pp];
        double DENDRO_122 = grad_2_Gt0;
        double DENDRO_123 = gt5[pp];
        double DENDRO_124 = grad_0_Gt2;
        double DENDRO_125 = gt4[pp];
        double DENDRO_126 = grad_0_Gt1;
        double DENDRO_127 = grad_0_Gt0;

        // initialize reduction for
        double DENDRO_128 = 0;

        DENDRO_128 += DENDRO_32;
        DENDRO_128 += DENDRO_31;
        DENDRO_128 += DENDRO_30;
        DENDRO_30 = 4;

        // initialize reduction for
        DENDRO_31 = 0;

        DENDRO_31 += DENDRO_35;
        DENDRO_31 += DENDRO_34;
        DENDRO_31 += DENDRO_33;

        // initialize reduction for
        DENDRO_32 = 0;

        DENDRO_32 += DENDRO_37;
        DENDRO_32 += DENDRO_36;

        // initialize reduction for
        DENDRO_33 = 0;

        DENDRO_33 += DENDRO_64;
        DENDRO_33 += DENDRO_39;
        DENDRO_33 += DENDRO_38;

        // initialize reduction for
        DENDRO_34 = 0;

        DENDRO_34 += DENDRO_65;
        DENDRO_34 += DENDRO_39;
        DENDRO_34 += DENDRO_38;

        // initialize reduction for
        DENDRO_35 = 0;

        DENDRO_35 += DENDRO_68;
        DENDRO_35 += DENDRO_67;
        DENDRO_35 += DENDRO_66;

        // initialize reduction for
        DENDRO_36 = 0;

        DENDRO_36 += DENDRO_71;
        DENDRO_36 += DENDRO_70;
        DENDRO_36 += DENDRO_69;

        // initialize reduction for
        DENDRO_37 = 0;

        DENDRO_37 += DENDRO_74;
        DENDRO_37 += DENDRO_73;
        DENDRO_37 += DENDRO_72;

        // initialize reduction for
        DENDRO_38 = 0;

        DENDRO_38 += DENDRO_77;
        DENDRO_38 += DENDRO_76;
        DENDRO_38 += DENDRO_75;

        // initialize reduction for
        DENDRO_39 = 0;

        DENDRO_39 += DENDRO_80;
        DENDRO_39 += DENDRO_79;
        DENDRO_39 += DENDRO_78;

        // initialize reduction for
        DENDRO_64 = 0;

        DENDRO_64 += DENDRO_83;
        DENDRO_64 += DENDRO_82;
        DENDRO_64 += DENDRO_81;

        // initialize reduction for
        DENDRO_65 = 0;

        DENDRO_65 += DENDRO_85;
        DENDRO_65 += DENDRO_84;

        // initialize reduction for
        DENDRO_66 = 0;

        DENDRO_66 += DENDRO_87;
        DENDRO_66 += DENDRO_49;
        DENDRO_66 += DENDRO_86;

        // initialize reduction for
        DENDRO_67 = 0;

        DENDRO_67 += DENDRO_90;
        DENDRO_67 += DENDRO_89;
        DENDRO_67 += DENDRO_88;

        // initialize reduction for
        DENDRO_68 = 0;

        DENDRO_68 += DENDRO_91;
        DENDRO_68 += DENDRO_49;
        DENDRO_68 += DENDRO_86;

        // initialize reduction for
        DENDRO_49 = 0;

        DENDRO_49 += DENDRO_44;
        DENDRO_49 += DENDRO_92;
        DENDRO_49 += DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_94;
        DENDRO_8 += DENDRO_93;
        DENDRO_8 += DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 0;

        DENDRO_16 += DENDRO_95;
        DENDRO_16 += DENDRO_22;

        // initialize reduction for
        DENDRO_22 = 0;

        DENDRO_22 += DENDRO_6;
        DENDRO_22 += DENDRO_43;
        DENDRO_22 += DENDRO_96;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_98;
        DENDRO_6 += DENDRO_51;
        DENDRO_6 += DENDRO_97;

        // initialize reduction for
        DENDRO_44 = 0;

        DENDRO_44 += DENDRO_52;
        DENDRO_44 += DENDRO_43;
        DENDRO_44 += DENDRO_96;

        // initialize reduction for
        DENDRO_43 = 0;

        DENDRO_43 += DENDRO_50;
        DENDRO_43 += DENDRO_15;
        DENDRO_43 += DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_100;
        DENDRO_1 += DENDRO_99;
        DENDRO_1 += DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 0;

        DENDRO_13 += DENDRO_102;
        DENDRO_13 += DENDRO_21;
        DENDRO_13 += DENDRO_101;

        // initialize reduction for
        DENDRO_15 = 0;

        DENDRO_15 += DENDRO_18;
        DENDRO_15 += DENDRO_103;
        DENDRO_15 += DENDRO_19;

        // initialize reduction for
        DENDRO_18 = 0;

        DENDRO_18 += DENDRO_24;
        DENDRO_18 += DENDRO_55;
        DENDRO_18 += DENDRO_54;

        // initialize reduction for
        DENDRO_19 = 0;

        DENDRO_19 += DENDRO_12;
        DENDRO_19 += DENDRO_57;

        // initialize reduction for
        DENDRO_12 = 0;

        DENDRO_12 += DENDRO_58;
        DENDRO_12 += DENDRO_104;

        // initialize reduction for
        DENDRO_21 = 0;

        DENDRO_21 += DENDRO_106;
        DENDRO_21 += DENDRO_105;

        // initialize reduction for
        DENDRO_24 = 0;

        DENDRO_24 += DENDRO_62;
        DENDRO_24 += DENDRO_107;

        // initialize reduction for
        DENDRO_50 = 0;

        DENDRO_50 += DENDRO_110;
        DENDRO_50 += DENDRO_109;
        DENDRO_50 += DENDRO_108;
        DENDRO_51 = -2;

        // initialize reduction for
        DENDRO_52 = 0;

        DENDRO_52 += DENDRO_23;
        DENDRO_52 += DENDRO_111;

        // initialize reduction for
        DENDRO_23 = 0;

        DENDRO_23 += DENDRO_26;
        DENDRO_23 += DENDRO_25;

        // initialize reduction for
        DENDRO_25 = 0;

        DENDRO_25 += DENDRO_28;
        DENDRO_25 += DENDRO_27;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_29;
        DENDRO_26 *= DENDRO_4;
        DENDRO_26 *= DENDRO_40;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_29;
        DENDRO_4 *= DENDRO_14;
        DENDRO_4 *= DENDRO_41;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_29;
        DENDRO_14 *= DENDRO_20;
        DENDRO_14 *= DENDRO_42;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_112;
        DENDRO_20 *= DENDRO_0;
        DENDRO_20 *= DENDRO_60;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_112;
        DENDRO_27 *= DENDRO_3;
        DENDRO_27 *= DENDRO_113;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_112;
        DENDRO_28 *= DENDRO_5;
        DENDRO_28 *= DENDRO_114;

        // initialize reduction for
        DENDRO_41 = 1;

        DENDRO_41 *= DENDRO_116;
        DENDRO_41 *= DENDRO_7;
        DENDRO_41 *= DENDRO_115;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_116;
        DENDRO_54 *= DENDRO_9;
        DENDRO_54 *= DENDRO_117;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_116;
        DENDRO_55 *= DENDRO_11;
        DENDRO_55 *= DENDRO_118;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_46;
        DENDRO_57 *= DENDRO_119;
        DENDRO_57 *= DENDRO_63;

        // initialize reduction for
        DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_46;
        DENDRO_58 *= DENDRO_47;
        DENDRO_58 *= DENDRO_10;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_46;
        DENDRO_10 *= DENDRO_45;
        DENDRO_10 *= DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_46;
        DENDRO_17 *= DENDRO_53;
        DENDRO_17 *= DENDRO_48;

        // initialize reduction for
        DENDRO_45 = 1;

        DENDRO_45 *= DENDRO_46;
        DENDRO_45 *= DENDRO_120;
        DENDRO_45 *= DENDRO_56;

        // initialize reduction for
        DENDRO_47 = 1;

        DENDRO_47 *= DENDRO_46;
        DENDRO_47 *= DENDRO_122;
        DENDRO_47 *= DENDRO_121;

        // initialize reduction for
        DENDRO_53 = 1;

        DENDRO_53 *= DENDRO_46;
        DENDRO_53 *= DENDRO_124;
        DENDRO_53 *= DENDRO_123;

        // initialize reduction for
        DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_46;
        DENDRO_56 *= DENDRO_126;
        DENDRO_56 *= DENDRO_125;

        // initialize reduction for
        DENDRO_60 = 1;

        DENDRO_60 *= DENDRO_46;
        DENDRO_60 *= DENDRO_127;
        DENDRO_60 *= DENDRO_48;

        // initialize reduction for
        DENDRO_46 = 1;

        DENDRO_46 *= DENDRO_30;
        DENDRO_46 *= DENDRO_7;
        DENDRO_46 *= DENDRO_128;

        // initialize reduction for
        DENDRO_62 = 1;

        DENDRO_62 *= DENDRO_30;
        DENDRO_62 *= DENDRO_7;
        DENDRO_62 *= DENDRO_31;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_30;
        DENDRO_31 *= DENDRO_7;
        DENDRO_31 *= DENDRO_32;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_30;
        DENDRO_32 *= DENDRO_0;
        DENDRO_32 *= DENDRO_33;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_30;
        DENDRO_33 *= DENDRO_0;
        DENDRO_33 *= DENDRO_34;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_30;
        DENDRO_34 *= DENDRO_0;
        DENDRO_34 *= DENDRO_35;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_30;
        DENDRO_35 *= DENDRO_0;
        DENDRO_35 *= DENDRO_36;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_30;
        DENDRO_36 *= DENDRO_0;
        DENDRO_36 *= DENDRO_37;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_30;
        DENDRO_37 *= DENDRO_0;
        DENDRO_37 *= DENDRO_38;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_30;
        DENDRO_38 *= DENDRO_9;
        DENDRO_38 *= DENDRO_39;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_30;
        DENDRO_39 *= DENDRO_9;
        DENDRO_39 *= DENDRO_64;

        // initialize reduction for
        DENDRO_63 = 1;

        DENDRO_63 *= DENDRO_30;
        DENDRO_63 *= DENDRO_9;
        DENDRO_63 *= DENDRO_65;

        // initialize reduction for
        DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_30;
        DENDRO_64 *= DENDRO_3;
        DENDRO_64 *= DENDRO_66;

        // initialize reduction for
        DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_30;
        DENDRO_65 *= DENDRO_3;
        DENDRO_65 *= DENDRO_67;

        // initialize reduction for
        DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_30;
        DENDRO_66 *= DENDRO_3;
        DENDRO_66 *= DENDRO_68;

        // initialize reduction for
        DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_30;
        DENDRO_67 *= DENDRO_3;
        DENDRO_67 *= DENDRO_49;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_30;
        DENDRO_49 *= DENDRO_3;
        DENDRO_49 *= DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_30;
        DENDRO_8 *= DENDRO_3;
        DENDRO_8 *= DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_30;
        DENDRO_16 *= DENDRO_5;
        DENDRO_16 *= DENDRO_22;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_30;
        DENDRO_22 *= DENDRO_5;
        DENDRO_22 *= DENDRO_6;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_30;
        DENDRO_6 *= DENDRO_5;
        DENDRO_6 *= DENDRO_44;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_30;
        DENDRO_44 *= DENDRO_5;
        DENDRO_44 *= DENDRO_43;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_30;
        DENDRO_43 *= DENDRO_5;
        DENDRO_43 *= DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_30;
        DENDRO_1 *= DENDRO_5;
        DENDRO_1 *= DENDRO_13;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_30;
        DENDRO_13 *= DENDRO_11;
        DENDRO_13 *= DENDRO_15;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_30;
        DENDRO_15 *= DENDRO_11;
        DENDRO_15 *= DENDRO_18;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_30;
        DENDRO_18 *= DENDRO_11;
        DENDRO_18 *= DENDRO_19;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_2;
        DENDRO_19 *= DENDRO_0;
        DENDRO_19 *= DENDRO_12;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_21;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_2;
        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_24;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_51;
        DENDRO_5 *= DENDRO_59;
        DENDRO_5 *= DENDRO_50;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_7;
        DENDRO_12 *= DENDRO_52;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_9;
        DENDRO_7 *= DENDRO_23;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_25;

        // initialize reduction for
        DENDRO_11 = 0;

        DENDRO_11 += DENDRO_61;
        DENDRO_11 += DENDRO_14;
        DENDRO_11 += DENDRO_4;
        DENDRO_11 += DENDRO_26;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_18;
        DENDRO_4 += DENDRO_15;
        DENDRO_4 += DENDRO_13;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_43;
        DENDRO_4 += DENDRO_44;
        DENDRO_4 += DENDRO_6;
        DENDRO_4 += DENDRO_22;
        DENDRO_4 += DENDRO_16;
        DENDRO_4 += DENDRO_8;
        DENDRO_4 += DENDRO_49;
        DENDRO_4 += DENDRO_67;
        DENDRO_4 += DENDRO_66;
        DENDRO_4 += DENDRO_65;
        DENDRO_4 += DENDRO_64;
        DENDRO_4 += DENDRO_63;
        DENDRO_4 += DENDRO_39;
        DENDRO_4 += DENDRO_38;
        DENDRO_4 += DENDRO_37;
        DENDRO_4 += DENDRO_36;
        DENDRO_4 += DENDRO_35;
        DENDRO_4 += DENDRO_34;
        DENDRO_4 += DENDRO_33;
        DENDRO_4 += DENDRO_32;
        DENDRO_4 += DENDRO_31;
        DENDRO_4 += DENDRO_62;
        DENDRO_4 += DENDRO_46;
        DENDRO_4 += DENDRO_60;
        DENDRO_4 += DENDRO_56;
        DENDRO_4 += DENDRO_53;
        DENDRO_4 += DENDRO_47;
        DENDRO_4 += DENDRO_45;
        DENDRO_4 += DENDRO_17;
        DENDRO_4 += DENDRO_10;
        DENDRO_4 += DENDRO_58;
        DENDRO_4 += DENDRO_57;
        DENDRO_4 += DENDRO_55;
        DENDRO_4 += DENDRO_54;
        DENDRO_4 += DENDRO_41;
        DENDRO_4 += DENDRO_28;
        DENDRO_4 += DENDRO_27;
        DENDRO_4 += DENDRO_20;
        DENDRO_1 = DENDRO_59;

        DENDRO_1 *= DENDRO_59;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_9;
        DENDRO_6 += DENDRO_7;
        DENDRO_6 += DENDRO_12;
        DENDRO_6 += DENDRO_5;
        DENDRO_6 += DENDRO_3;
        DENDRO_6 += DENDRO_0;
        DENDRO_6 += DENDRO_19;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_59;
        DENDRO_0 *= DENDRO_11;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_29;
        DENDRO_2 *= DENDRO_42;
        DENDRO_2 *= DENDRO_40;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_1;
        DENDRO_3 *= DENDRO_4;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_48;
        DENDRO_1 *= DENDRO_6;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_2;
        DENDRO_4 += DENDRO_0;
        DENDRO_0 = DENDRO_59;

        DENDRO_0 *= DENDRO_59;
        DENDRO_0 = 1 / DENDRO_0;
        DENDRO_1 = 1.0 / 4.0;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_0;
        DENDRO_2 *= DENDRO_4;
        DENDRO_RIJ2 = DENDRO_2;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C2_k2_4;
        double DENDRO_2  = 2;
        double DENDRO_3  = DENDRO_igt2;
        double DENDRO_4  = DENDRO_C2_k2_2;
        double DENDRO_5  = DENDRO_igt1;
        double DENDRO_6  = DENDRO_C2_k2_1;
        double DENDRO_7  = DENDRO_igt5;
        double DENDRO_8  = DENDRO_C2_k2_5;
        double DENDRO_9  = DENDRO_igt3;
        double DENDRO_10 = DENDRO_C2_k2_3;
        double DENDRO_11 = DENDRO_igt0;
        double DENDRO_12 = DENDRO_C2_k2_0;
        double DENDRO_13 = DENDRO_C2_k1_4;
        double DENDRO_14 = DENDRO_C2_k1_2;
        double DENDRO_15 = DENDRO_C2_k1_1;
        double DENDRO_16 = DENDRO_C2_k1_5;
        double DENDRO_17 = DENDRO_C2_k1_3;
        double DENDRO_18 = DENDRO_C2_k1_0;
        double DENDRO_19 = DENDRO_C2_k0_4;
        double DENDRO_20 = DENDRO_C2_k0_2;
        double DENDRO_21 = DENDRO_C2_k0_1;
        double DENDRO_22 = DENDRO_C2_k0_5;
        double DENDRO_23 = DENDRO_C2_k0_3;
        double DENDRO_24 = DENDRO_C2_k0_0;

        // initialize reduction for
        double DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_1;
        DENDRO_25 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_4;
        DENDRO_26 *= DENDRO_3;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_2;
        DENDRO_4 *= DENDRO_6;
        DENDRO_4 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_8;
        DENDRO_27 *= DENDRO_7;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_10;
        DENDRO_8 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_12;
        DENDRO_28 *= DENDRO_11;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_2;
        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_2;
        DENDRO_29 *= DENDRO_14;
        DENDRO_29 *= DENDRO_3;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_2;
        DENDRO_14 *= DENDRO_15;
        DENDRO_14 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_16;
        DENDRO_30 *= DENDRO_7;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_17;
        DENDRO_16 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_18;
        DENDRO_31 *= DENDRO_11;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_2;
        DENDRO_18 *= DENDRO_19;
        DENDRO_18 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_20;
        DENDRO_32 *= DENDRO_3;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_2;
        DENDRO_20 *= DENDRO_21;
        DENDRO_20 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_22;
        DENDRO_33 *= DENDRO_7;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_23;
        DENDRO_22 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_24;
        DENDRO_34 *= DENDRO_11;
        DENDRO_24        = grad_2_chi;
        double DENDRO_35 = grad_1_chi;
        double DENDRO_36 = grad_0_chi;
        double DENDRO_37 = DENDRO_C1_k1_5;
        double DENDRO_38 = DENDRO_C1_k1_4;
        double DENDRO_39 = DENDRO_C1_k1_2;
        double DENDRO_40 = DENDRO_C1_k1_1;
        double DENDRO_41 = DENDRO_C1_k1_0;
        double DENDRO_42 = DENDRO_C1_k2_4;
        double DENDRO_43 = DENDRO_C1_k2_3;
        double DENDRO_44 = DENDRO_C1_k1_3;
        double DENDRO_45 = DENDRO_C1_k0_4;
        double DENDRO_46 = DENDRO_C1_k0_3;
        double DENDRO_47 = DENDRO_C1_k2_1;
        double DENDRO_48 = DENDRO_C1_k0_1;
        double DENDRO_49 = grad2_1_2_chi;
        double DENDRO_50 = chi[pp];
        double DENDRO_51 = -3;
        double DENDRO_52 = grad2_0_2_chi;
        double DENDRO_53 = grad2_0_1_chi;

        // initialize reduction for
        double DENDRO_54 = 0;

        DENDRO_54 += DENDRO_28;
        DENDRO_54 += DENDRO_8;
        DENDRO_54 += DENDRO_27;
        DENDRO_54 += DENDRO_4;
        DENDRO_54 += DENDRO_26;
        DENDRO_54 += DENDRO_25;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_31;
        DENDRO_4 += DENDRO_16;
        DENDRO_4 += DENDRO_30;
        DENDRO_4 += DENDRO_14;
        DENDRO_4 += DENDRO_29;
        DENDRO_4 += DENDRO_12;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_34;
        DENDRO_8 += DENDRO_22;
        DENDRO_8 += DENDRO_33;
        DENDRO_8 += DENDRO_20;
        DENDRO_8 += DENDRO_32;
        DENDRO_8 += DENDRO_18;
        DENDRO_12 = grad2_2_2_chi;
        DENDRO_14 = DENDRO_24;

        DENDRO_14 *= DENDRO_24;
        DENDRO_16 = grad2_1_1_chi;
        DENDRO_18 = DENDRO_35;

        DENDRO_18 *= DENDRO_35;
        DENDRO_20 = grad2_0_0_chi;
        DENDRO_22 = DENDRO_36;

        DENDRO_22 *= DENDRO_36;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_37;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_38;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_39;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_2;
        DENDRO_28 *= DENDRO_40;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_2;
        DENDRO_29 *= DENDRO_41;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_2;
        DENDRO_30 *= DENDRO_38;
        DENDRO_30 *= DENDRO_1;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_42;
        DENDRO_31 *= DENDRO_10;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_37;
        DENDRO_32 *= DENDRO_10;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_43;
        DENDRO_33 *= DENDRO_1;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_2;
        DENDRO_34 *= DENDRO_44;
        DENDRO_34 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_38;
        DENDRO_55 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_2;
        DENDRO_56 *= DENDRO_38;
        DENDRO_56 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_44;
        DENDRO_57 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_2;
        DENDRO_58 *= DENDRO_40;
        DENDRO_58 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_59 = 1;

        DENDRO_59 *= DENDRO_45;
        DENDRO_59 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_60 = 1;

        DENDRO_60 *= DENDRO_2;
        DENDRO_60 *= DENDRO_39;
        DENDRO_60 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_61 = 1;

        DENDRO_61 *= DENDRO_46;
        DENDRO_61 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_62 = 1;

        DENDRO_62 *= DENDRO_2;
        DENDRO_62 *= DENDRO_39;
        DENDRO_62 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_63 = 1;

        DENDRO_63 *= DENDRO_42;
        DENDRO_63 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_2;
        DENDRO_64 *= DENDRO_37;
        DENDRO_64 *= DENDRO_6;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_47;
        DENDRO_37 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_2;
        DENDRO_65 *= DENDRO_40;
        DENDRO_65 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_38;
        DENDRO_66 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_2;
        DENDRO_67 *= DENDRO_38;
        DENDRO_67 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_68 = 1;

        DENDRO_68 *= DENDRO_40;
        DENDRO_68 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_69 = 1;

        DENDRO_69 *= DENDRO_2;
        DENDRO_69 *= DENDRO_41;
        DENDRO_69 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_70 = 1;

        DENDRO_70 *= DENDRO_45;
        DENDRO_70 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_71 = 1;

        DENDRO_71 *= DENDRO_2;
        DENDRO_71 *= DENDRO_39;
        DENDRO_71 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_72 = 1;

        DENDRO_72 *= DENDRO_48;
        DENDRO_72 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_73 = 1;

        DENDRO_73 *= DENDRO_2;
        DENDRO_73 *= DENDRO_39;
        DENDRO_73 *= DENDRO_10;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_43;
        DENDRO_39 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_74 = 1;

        DENDRO_74 *= DENDRO_2;
        DENDRO_74 *= DENDRO_38;
        DENDRO_74 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_75 = 1;

        DENDRO_75 *= DENDRO_47;
        DENDRO_75 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_76 = 1;

        DENDRO_76 *= DENDRO_2;
        DENDRO_76 *= DENDRO_40;
        DENDRO_76 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_77 = 1;

        DENDRO_77 *= DENDRO_44;
        DENDRO_77 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_78 = 1;

        DENDRO_78 *= DENDRO_2;
        DENDRO_78 *= DENDRO_44;
        DENDRO_78 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_79 = 1;

        DENDRO_79 *= DENDRO_40;
        DENDRO_79 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_80 = 1;

        DENDRO_80 *= DENDRO_2;
        DENDRO_80 *= DENDRO_41;
        DENDRO_80 *= DENDRO_23;

        // initialize reduction for
        DENDRO_41 = 1;

        DENDRO_41 *= DENDRO_46;
        DENDRO_41 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_81 = 1;

        DENDRO_81 *= DENDRO_2;
        DENDRO_81 *= DENDRO_40;
        DENDRO_81 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_82 = 1;

        DENDRO_82 *= DENDRO_48;
        DENDRO_82 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_83 = 1;

        DENDRO_83 *= DENDRO_2;
        DENDRO_83 *= DENDRO_50;
        DENDRO_83 *= DENDRO_49;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_51;
        DENDRO_49 *= DENDRO_35;
        DENDRO_49 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_84 = 1;

        DENDRO_84 *= DENDRO_2;
        DENDRO_84 *= DENDRO_50;
        DENDRO_84 *= DENDRO_52;

        // initialize reduction for
        DENDRO_52 = 1;

        DENDRO_52 *= DENDRO_51;
        DENDRO_52 *= DENDRO_36;
        DENDRO_52 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_85 = 1;

        DENDRO_85 *= DENDRO_2;
        DENDRO_85 *= DENDRO_50;
        DENDRO_85 *= DENDRO_53;

        // initialize reduction for
        DENDRO_53 = 1;

        DENDRO_53 *= DENDRO_51;
        DENDRO_53 *= DENDRO_36;
        DENDRO_53 *= DENDRO_35;

        // initialize reduction for
        double DENDRO_86 = 1;

        DENDRO_86 *= DENDRO_24;
        DENDRO_86 *= DENDRO_54;

        // initialize reduction for
        double DENDRO_87 = 1;

        DENDRO_87 *= DENDRO_35;
        DENDRO_87 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_88 = 1;

        DENDRO_88 *= DENDRO_36;
        DENDRO_88 *= DENDRO_8;

        // initialize reduction for
        double DENDRO_89 = 1;

        DENDRO_89 *= DENDRO_2;
        DENDRO_89 *= DENDRO_50;
        DENDRO_89 *= DENDRO_12;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_51;
        DENDRO_12 *= DENDRO_14;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_2;
        DENDRO_14 *= DENDRO_50;
        DENDRO_14 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_90 = 1;

        DENDRO_90 *= DENDRO_51;
        DENDRO_90 *= DENDRO_18;

        // initialize reduction for
        double DENDRO_91 = 1;

        DENDRO_91 *= DENDRO_2;
        DENDRO_91 *= DENDRO_50;
        DENDRO_91 *= DENDRO_20;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_51;
        DENDRO_20 *= DENDRO_22;
        DENDRO_22        = -1;
        DENDRO_51        = 12;

        // initialize reduction for
        double DENDRO_92 = 0;

        DENDRO_92 += DENDRO_42;
        DENDRO_92 += DENDRO_25;
        DENDRO_25 = 4;

        // initialize reduction for
        DENDRO_42 = 0;

        DENDRO_42 += DENDRO_43;
        DENDRO_42 += DENDRO_26;

        // initialize reduction for
        DENDRO_26 = 0;

        DENDRO_26 += DENDRO_47;
        DENDRO_26 += DENDRO_27;

        // initialize reduction for
        DENDRO_43 = 0;

        DENDRO_43 += DENDRO_45;
        DENDRO_43 += DENDRO_27;

        // initialize reduction for
        DENDRO_27 = 0;

        DENDRO_27 += DENDRO_46;
        DENDRO_27 += DENDRO_28;

        // initialize reduction for
        DENDRO_28 = 0;

        DENDRO_28 += DENDRO_48;
        DENDRO_28 += DENDRO_29;
        DENDRO_29         = grad2_1_2_gt3;
        DENDRO_45         = -4.0;
        DENDRO_46         = grad2_0_2_gt3;
        DENDRO_47         = grad2_0_1_gt3;
        DENDRO_48         = grad2_2_2_gt3;
        double DENDRO_93  = -2.0;
        double DENDRO_94  = grad2_1_1_gt3;
        double DENDRO_95  = grad2_0_0_gt3;
        double DENDRO_96  = gt4[pp];
        double DENDRO_97  = grad_1_Gt2;
        double DENDRO_98  = 4.0;
        double DENDRO_99  = gt3[pp];
        double DENDRO_100 = grad_1_Gt1;
        double DENDRO_101 = gt1[pp];
        double DENDRO_102 = grad_1_Gt0;

        // initialize reduction for
        double DENDRO_103 = 0;

        DENDRO_103 += DENDRO_31;
        DENDRO_103 += DENDRO_30;

        // initialize reduction for
        DENDRO_30 = 0;

        DENDRO_30 += DENDRO_33;
        DENDRO_30 += DENDRO_32;

        // initialize reduction for
        DENDRO_31 = 0;

        DENDRO_31 += DENDRO_55;
        DENDRO_31 += DENDRO_34;

        // initialize reduction for
        DENDRO_32 = 0;

        DENDRO_32 += DENDRO_57;
        DENDRO_32 += DENDRO_56;

        // initialize reduction for
        DENDRO_33 = 0;

        DENDRO_33 += DENDRO_59;
        DENDRO_33 += DENDRO_58;

        // initialize reduction for
        DENDRO_34 = 0;

        DENDRO_34 += DENDRO_61;
        DENDRO_34 += DENDRO_60;

        // initialize reduction for
        DENDRO_55 = 0;

        DENDRO_55 += DENDRO_63;
        DENDRO_55 += DENDRO_62;

        // initialize reduction for
        DENDRO_56 = 0;

        DENDRO_56 += DENDRO_37;
        DENDRO_56 += DENDRO_64;

        // initialize reduction for
        DENDRO_37 = 0;

        DENDRO_37 += DENDRO_66;
        DENDRO_37 += DENDRO_65;

        // initialize reduction for
        DENDRO_57 = 0;

        DENDRO_57 += DENDRO_68;
        DENDRO_57 += DENDRO_67;

        // initialize reduction for
        DENDRO_58 = 0;

        DENDRO_58 += DENDRO_70;
        DENDRO_58 += DENDRO_69;

        // initialize reduction for
        DENDRO_59 = 0;

        DENDRO_59 += DENDRO_72;
        DENDRO_59 += DENDRO_71;

        // initialize reduction for
        DENDRO_60 = 0;

        DENDRO_60 += DENDRO_39;
        DENDRO_60 += DENDRO_73;

        // initialize reduction for
        DENDRO_39 = 0;

        DENDRO_39 += DENDRO_75;
        DENDRO_39 += DENDRO_74;

        // initialize reduction for
        DENDRO_61 = 0;

        DENDRO_61 += DENDRO_77;
        DENDRO_61 += DENDRO_76;

        // initialize reduction for
        DENDRO_62 = 0;

        DENDRO_62 += DENDRO_79;
        DENDRO_62 += DENDRO_78;

        // initialize reduction for
        DENDRO_63 = 0;

        DENDRO_63 += DENDRO_41;
        DENDRO_63 += DENDRO_80;

        // initialize reduction for
        DENDRO_41 = 0;

        DENDRO_41 += DENDRO_82;
        DENDRO_41 += DENDRO_81;

        // initialize reduction for
        DENDRO_64 = 0;

        DENDRO_64 += DENDRO_49;
        DENDRO_64 += DENDRO_83;

        // initialize reduction for
        DENDRO_49 = 0;

        DENDRO_49 += DENDRO_52;
        DENDRO_49 += DENDRO_84;

        // initialize reduction for
        DENDRO_52 = 0;

        DENDRO_52 += DENDRO_53;
        DENDRO_52 += DENDRO_85;

        // initialize reduction for
        DENDRO_53 = 0;

        DENDRO_53 += DENDRO_88;
        DENDRO_53 += DENDRO_87;
        DENDRO_53 += DENDRO_86;
        DENDRO_65 = -2;

        // initialize reduction for
        DENDRO_66 = 0;

        DENDRO_66 += DENDRO_12;
        DENDRO_66 += DENDRO_89;

        // initialize reduction for
        DENDRO_12 = 0;

        DENDRO_12 += DENDRO_90;
        DENDRO_12 += DENDRO_14;

        // initialize reduction for
        DENDRO_14 = 0;

        DENDRO_14 += DENDRO_20;
        DENDRO_14 += DENDRO_91;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_22;
        DENDRO_20 *= DENDRO_10;
        DENDRO_20 *= DENDRO_24;

        // initialize reduction for
        DENDRO_24 = 1;

        DENDRO_24 *= DENDRO_22;
        DENDRO_24 *= DENDRO_17;
        DENDRO_24 *= DENDRO_35;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_22;
        DENDRO_35 *= DENDRO_23;
        DENDRO_35 *= DENDRO_36;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_51;
        DENDRO_36 *= DENDRO_38;
        DENDRO_36 *= DENDRO_13;
        DENDRO_36 *= DENDRO_7;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_51;
        DENDRO_13 *= DENDRO_44;
        DENDRO_13 *= DENDRO_17;
        DENDRO_13 *= DENDRO_9;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_51;
        DENDRO_17 *= DENDRO_40;
        DENDRO_17 *= DENDRO_15;
        DENDRO_17 *= DENDRO_11;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_25;
        DENDRO_15 *= DENDRO_1;
        DENDRO_15 *= DENDRO_7;
        DENDRO_15 *= DENDRO_92;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_25;
        DENDRO_1 *= DENDRO_10;
        DENDRO_1 *= DENDRO_9;
        DENDRO_1 *= DENDRO_42;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_25;
        DENDRO_10 *= DENDRO_6;
        DENDRO_10 *= DENDRO_11;
        DENDRO_10 *= DENDRO_26;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_25;
        DENDRO_6 *= DENDRO_19;
        DENDRO_6 *= DENDRO_7;
        DENDRO_6 *= DENDRO_43;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_25;
        DENDRO_19 *= DENDRO_23;
        DENDRO_19 *= DENDRO_9;
        DENDRO_19 *= DENDRO_27;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_25;
        DENDRO_23 *= DENDRO_21;
        DENDRO_23 *= DENDRO_11;
        DENDRO_23 *= DENDRO_28;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_45;
        DENDRO_21 *= DENDRO_0;
        DENDRO_21 *= DENDRO_29;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_45;
        DENDRO_26 *= DENDRO_3;
        DENDRO_26 *= DENDRO_46;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_45;
        DENDRO_27 *= DENDRO_5;
        DENDRO_27 *= DENDRO_47;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_93;
        DENDRO_28 *= DENDRO_7;
        DENDRO_28 *= DENDRO_48;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_93;
        DENDRO_29 *= DENDRO_9;
        DENDRO_29 *= DENDRO_94;

        // initialize reduction for
        DENDRO_42 = 1;

        DENDRO_42 *= DENDRO_93;
        DENDRO_42 *= DENDRO_11;
        DENDRO_42 *= DENDRO_95;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_98;
        DENDRO_43 *= DENDRO_97;
        DENDRO_43 *= DENDRO_96;

        // initialize reduction for
        DENDRO_45 = 1;

        DENDRO_45 *= DENDRO_98;
        DENDRO_45 *= DENDRO_100;
        DENDRO_45 *= DENDRO_99;

        // initialize reduction for
        DENDRO_46 = 1;

        DENDRO_46 *= DENDRO_98;
        DENDRO_46 *= DENDRO_102;
        DENDRO_46 *= DENDRO_101;

        // initialize reduction for
        DENDRO_47 = 1;

        DENDRO_47 *= DENDRO_98;
        DENDRO_47 *= DENDRO_38;
        DENDRO_47 *= DENDRO_54;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_98;
        DENDRO_38 *= DENDRO_44;
        DENDRO_38 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_98;
        DENDRO_4 *= DENDRO_40;
        DENDRO_4 *= DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_25;
        DENDRO_8 *= DENDRO_0;
        DENDRO_8 *= DENDRO_103;

        // initialize reduction for
        DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_25;
        DENDRO_40 *= DENDRO_0;
        DENDRO_40 *= DENDRO_30;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_25;
        DENDRO_30 *= DENDRO_0;
        DENDRO_30 *= DENDRO_31;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_25;
        DENDRO_31 *= DENDRO_0;
        DENDRO_31 *= DENDRO_32;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_25;
        DENDRO_32 *= DENDRO_0;
        DENDRO_32 *= DENDRO_33;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_25;
        DENDRO_33 *= DENDRO_0;
        DENDRO_33 *= DENDRO_34;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_25;
        DENDRO_34 *= DENDRO_3;
        DENDRO_34 *= DENDRO_55;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_25;
        DENDRO_44 *= DENDRO_3;
        DENDRO_44 *= DENDRO_56;

        // initialize reduction for
        DENDRO_48 = 1;

        DENDRO_48 *= DENDRO_25;
        DENDRO_48 *= DENDRO_3;
        DENDRO_48 *= DENDRO_37;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_25;
        DENDRO_37 *= DENDRO_3;
        DENDRO_37 *= DENDRO_57;

        // initialize reduction for
        DENDRO_51 = 1;

        DENDRO_51 *= DENDRO_25;
        DENDRO_51 *= DENDRO_3;
        DENDRO_51 *= DENDRO_58;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_25;
        DENDRO_54 *= DENDRO_3;
        DENDRO_54 *= DENDRO_59;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_25;
        DENDRO_55 *= DENDRO_5;
        DENDRO_55 *= DENDRO_60;

        // initialize reduction for
        DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_25;
        DENDRO_56 *= DENDRO_5;
        DENDRO_56 *= DENDRO_39;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_25;
        DENDRO_39 *= DENDRO_5;
        DENDRO_39 *= DENDRO_61;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_25;
        DENDRO_57 *= DENDRO_5;
        DENDRO_57 *= DENDRO_62;

        // initialize reduction for
        DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_25;
        DENDRO_58 *= DENDRO_5;
        DENDRO_58 *= DENDRO_63;

        // initialize reduction for
        DENDRO_59 = 1;

        DENDRO_59 *= DENDRO_25;
        DENDRO_59 *= DENDRO_5;
        DENDRO_59 *= DENDRO_41;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_0;
        DENDRO_25 *= DENDRO_64;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_49;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_2;
        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_52;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_65;
        DENDRO_5 *= DENDRO_50;
        DENDRO_5 *= DENDRO_53;

        // initialize reduction for
        DENDRO_41 = 1;

        DENDRO_41 *= DENDRO_7;
        DENDRO_41 *= DENDRO_66;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_9;
        DENDRO_7 *= DENDRO_12;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_14;

        // initialize reduction for
        DENDRO_11 = 0;

        DENDRO_11 += DENDRO_16;
        DENDRO_11 += DENDRO_35;
        DENDRO_11 += DENDRO_24;
        DENDRO_11 += DENDRO_20;

        // initialize reduction for
        DENDRO_12 = 0;

        DENDRO_12 += DENDRO_59;
        DENDRO_12 += DENDRO_58;
        DENDRO_12 += DENDRO_57;
        DENDRO_12 += DENDRO_39;
        DENDRO_12 += DENDRO_56;
        DENDRO_12 += DENDRO_55;
        DENDRO_12 += DENDRO_54;
        DENDRO_12 += DENDRO_51;
        DENDRO_12 += DENDRO_37;
        DENDRO_12 += DENDRO_48;
        DENDRO_12 += DENDRO_44;
        DENDRO_12 += DENDRO_34;
        DENDRO_12 += DENDRO_33;
        DENDRO_12 += DENDRO_32;
        DENDRO_12 += DENDRO_31;
        DENDRO_12 += DENDRO_30;
        DENDRO_12 += DENDRO_40;
        DENDRO_12 += DENDRO_8;
        DENDRO_12 += DENDRO_4;
        DENDRO_12 += DENDRO_38;
        DENDRO_12 += DENDRO_47;
        DENDRO_12 += DENDRO_46;
        DENDRO_12 += DENDRO_45;
        DENDRO_12 += DENDRO_43;
        DENDRO_12 += DENDRO_42;
        DENDRO_12 += DENDRO_29;
        DENDRO_12 += DENDRO_28;
        DENDRO_12 += DENDRO_27;
        DENDRO_12 += DENDRO_26;
        DENDRO_12 += DENDRO_21;
        DENDRO_12 += DENDRO_23;
        DENDRO_12 += DENDRO_19;
        DENDRO_12 += DENDRO_6;
        DENDRO_12 += DENDRO_10;
        DENDRO_12 += DENDRO_1;
        DENDRO_12 += DENDRO_15;
        DENDRO_12 += DENDRO_17;
        DENDRO_12 += DENDRO_13;
        DENDRO_12 += DENDRO_36;
        DENDRO_1 = DENDRO_50;

        DENDRO_1 *= DENDRO_50;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_9;
        DENDRO_4 += DENDRO_7;
        DENDRO_4 += DENDRO_41;
        DENDRO_4 += DENDRO_5;
        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_0;
        DENDRO_4 += DENDRO_25;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_50;
        DENDRO_0 *= DENDRO_11;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_12;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_99;
        DENDRO_1 *= DENDRO_4;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_22;
        DENDRO_3 *= DENDRO_18;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_2;
        DENDRO_4 += DENDRO_0;
        DENDRO_0 = DENDRO_50;

        DENDRO_0 *= DENDRO_50;
        DENDRO_0 = 1 / DENDRO_0;
        DENDRO_1 = 1.0 / 4.0;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_0;
        DENDRO_2 *= DENDRO_4;
        DENDRO_RIJ3 = DENDRO_2;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C2_k2_4;
        double DENDRO_2  = 2;
        double DENDRO_3  = DENDRO_igt2;
        double DENDRO_4  = DENDRO_C2_k2_2;
        double DENDRO_5  = DENDRO_igt1;
        double DENDRO_6  = DENDRO_C2_k2_1;
        double DENDRO_7  = DENDRO_igt5;
        double DENDRO_8  = DENDRO_C2_k2_5;
        double DENDRO_9  = DENDRO_igt3;
        double DENDRO_10 = DENDRO_C2_k2_3;
        double DENDRO_11 = DENDRO_igt0;
        double DENDRO_12 = DENDRO_C2_k2_0;
        double DENDRO_13 = DENDRO_C2_k1_4;
        double DENDRO_14 = DENDRO_C2_k1_2;
        double DENDRO_15 = DENDRO_C2_k1_1;
        double DENDRO_16 = DENDRO_C2_k1_5;
        double DENDRO_17 = DENDRO_C2_k1_3;
        double DENDRO_18 = DENDRO_C2_k1_0;
        double DENDRO_19 = DENDRO_C2_k0_4;
        double DENDRO_20 = DENDRO_C2_k0_2;
        double DENDRO_21 = DENDRO_C2_k0_1;
        double DENDRO_22 = DENDRO_C2_k0_5;
        double DENDRO_23 = DENDRO_C2_k0_3;
        double DENDRO_24 = DENDRO_C2_k0_0;

        // initialize reduction for
        double DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_1;
        DENDRO_25 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_4;
        DENDRO_26 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_6;
        DENDRO_27 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_8;
        DENDRO_28 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_10;
        DENDRO_29 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_12;
        DENDRO_30 *= DENDRO_11;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_2;
        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_2;
        DENDRO_31 *= DENDRO_14;
        DENDRO_31 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_15;
        DENDRO_32 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_16;
        DENDRO_33 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_17;
        DENDRO_34 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_18;
        DENDRO_35 *= DENDRO_11;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_2;
        DENDRO_18 *= DENDRO_19;
        DENDRO_18 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_2;
        DENDRO_36 *= DENDRO_20;
        DENDRO_36 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_2;
        DENDRO_37 *= DENDRO_21;
        DENDRO_37 *= DENDRO_5;

        // initialize reduction for
        double DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_22;
        DENDRO_38 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_23;
        DENDRO_39 *= DENDRO_9;

        // initialize reduction for
        double DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_24;
        DENDRO_40 *= DENDRO_11;
        DENDRO_24        = grad_2_chi;
        double DENDRO_41 = grad_1_chi;
        double DENDRO_42 = grad_0_chi;
        double DENDRO_43 = DENDRO_C1_k2_4;
        double DENDRO_44 = DENDRO_C1_k1_5;
        double DENDRO_45 = DENDRO_C1_k1_4;
        double DENDRO_46 = DENDRO_C1_k2_2;
        double DENDRO_47 = DENDRO_C1_k1_2;
        double DENDRO_48 = DENDRO_C1_k0_5;
        double DENDRO_49 = DENDRO_C1_k2_5;
        double DENDRO_50 = DENDRO_C1_k2_3;
        double DENDRO_51 = DENDRO_C1_k1_3;
        double DENDRO_52 = DENDRO_C1_k2_1;
        double DENDRO_53 = DENDRO_C1_k1_1;
        double DENDRO_54 = DENDRO_C1_k0_4;
        double DENDRO_55 = DENDRO_C1_k2_0;
        double DENDRO_56 = DENDRO_C1_k1_0;
        double DENDRO_57 = DENDRO_C1_k0_2;
        double DENDRO_58 = grad2_1_2_chi;
        double DENDRO_59 = chi[pp];
        double DENDRO_60 = -3;
        double DENDRO_61 = grad2_0_2_chi;
        double DENDRO_62 = grad2_0_1_chi;

        // initialize reduction for
        double DENDRO_63 = 0;

        DENDRO_63 += DENDRO_30;
        DENDRO_63 += DENDRO_29;
        DENDRO_63 += DENDRO_28;
        DENDRO_63 += DENDRO_27;
        DENDRO_63 += DENDRO_26;
        DENDRO_63 += DENDRO_25;

        // initialize reduction for
        DENDRO_25 = 0;

        DENDRO_25 += DENDRO_35;
        DENDRO_25 += DENDRO_34;
        DENDRO_25 += DENDRO_33;
        DENDRO_25 += DENDRO_32;
        DENDRO_25 += DENDRO_31;
        DENDRO_25 += DENDRO_12;

        // initialize reduction for
        DENDRO_12 = 0;

        DENDRO_12 += DENDRO_40;
        DENDRO_12 += DENDRO_39;
        DENDRO_12 += DENDRO_38;
        DENDRO_12 += DENDRO_37;
        DENDRO_12 += DENDRO_36;
        DENDRO_12 += DENDRO_18;
        DENDRO_18 = grad2_2_2_chi;
        DENDRO_26 = DENDRO_24;

        DENDRO_26 *= DENDRO_24;
        DENDRO_27 = grad2_1_1_chi;
        DENDRO_28 = DENDRO_41;

        DENDRO_28 *= DENDRO_41;
        DENDRO_29 = grad2_0_0_chi;
        DENDRO_30 = DENDRO_42;

        DENDRO_30 *= DENDRO_42;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_43;
        DENDRO_31 *= DENDRO_13;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_44;
        DENDRO_32 *= DENDRO_13;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_45;
        DENDRO_33 *= DENDRO_16;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_46;
        DENDRO_34 *= DENDRO_19;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_47;
        DENDRO_35 *= DENDRO_22;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_48;
        DENDRO_36 *= DENDRO_19;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_2;
        DENDRO_37 *= DENDRO_49;
        DENDRO_37 *= DENDRO_1;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_44;
        DENDRO_38 *= DENDRO_8;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_49;
        DENDRO_39 *= DENDRO_10;

        // initialize reduction for
        DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_43;
        DENDRO_40 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_44;
        DENDRO_64 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_45;
        DENDRO_65 *= DENDRO_8;

        // initialize reduction for
        double DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_50;
        DENDRO_66 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_44;
        DENDRO_67 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_68 = 1;

        DENDRO_68 *= DENDRO_51;
        DENDRO_68 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_69 = 1;

        DENDRO_69 *= DENDRO_52;
        DENDRO_69 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_70 = 1;

        DENDRO_70 *= DENDRO_53;
        DENDRO_70 *= DENDRO_22;

        // initialize reduction for
        double DENDRO_71 = 1;

        DENDRO_71 *= DENDRO_48;
        DENDRO_71 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_72 = 1;

        DENDRO_72 *= DENDRO_46;
        DENDRO_72 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_73 = 1;

        DENDRO_73 *= DENDRO_47;
        DENDRO_73 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_74 = 1;

        DENDRO_74 *= DENDRO_54;
        DENDRO_74 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_75 = 1;

        DENDRO_75 *= DENDRO_2;
        DENDRO_75 *= DENDRO_45;
        DENDRO_75 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_76 = 1;

        DENDRO_76 *= DENDRO_43;
        DENDRO_76 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_77 = 1;

        DENDRO_77 *= DENDRO_50;
        DENDRO_77 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_78 = 1;

        DENDRO_78 *= DENDRO_45;
        DENDRO_78 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_79 = 1;

        DENDRO_79 *= DENDRO_51;
        DENDRO_79 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_80 = 1;

        DENDRO_80 *= DENDRO_52;
        DENDRO_80 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_81 = 1;

        DENDRO_81 *= DENDRO_53;
        DENDRO_81 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_82 = 1;

        DENDRO_82 *= DENDRO_54;
        DENDRO_82 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_83 = 1;

        DENDRO_83 *= DENDRO_2;
        DENDRO_83 *= DENDRO_43;
        DENDRO_83 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_84 = 1;

        DENDRO_84 *= DENDRO_45;
        DENDRO_84 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_85 = 1;

        DENDRO_85 *= DENDRO_49;
        DENDRO_85 *= DENDRO_6;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_46;
        DENDRO_49 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_86 = 1;

        DENDRO_86 *= DENDRO_44;
        DENDRO_86 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_87 = 1;

        DENDRO_87 *= DENDRO_47;
        DENDRO_87 *= DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_43;
        DENDRO_8 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_88 = 1;

        DENDRO_88 *= DENDRO_45;
        DENDRO_88 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_89 = 1;

        DENDRO_89 *= DENDRO_47;
        DENDRO_89 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_90 = 1;

        DENDRO_90 *= DENDRO_52;
        DENDRO_90 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_91 = 1;

        DENDRO_91 *= DENDRO_44;
        DENDRO_91 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_92 = 1;

        DENDRO_92 *= DENDRO_53;
        DENDRO_92 *= DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_55;
        DENDRO_16 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_93 = 1;

        DENDRO_93 *= DENDRO_56;
        DENDRO_93 *= DENDRO_22;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_48;
        DENDRO_22 *= DENDRO_21;

        // initialize reduction for
        DENDRO_48 = 1;

        DENDRO_48 *= DENDRO_46;
        DENDRO_48 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_94 = 1;

        DENDRO_94 *= DENDRO_47;
        DENDRO_94 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_95 = 1;

        DENDRO_95 *= DENDRO_57;
        DENDRO_95 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_96 = 1;

        DENDRO_96 *= DENDRO_43;
        DENDRO_96 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_97 = 1;

        DENDRO_97 *= DENDRO_46;
        DENDRO_97 *= DENDRO_10;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_45;
        DENDRO_10 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_98 = 1;

        DENDRO_98 *= DENDRO_47;
        DENDRO_98 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_99 = 1;

        DENDRO_99 *= DENDRO_50;
        DENDRO_99 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_100 = 1;

        DENDRO_100 *= DENDRO_51;
        DENDRO_100 *= DENDRO_14;

        // initialize reduction for
        DENDRO_51 = 1;

        DENDRO_51 *= DENDRO_47;
        DENDRO_51 *= DENDRO_17;

        // initialize reduction for
        double DENDRO_101 = 1;

        DENDRO_101 *= DENDRO_52;
        DENDRO_101 *= DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_45;
        DENDRO_17 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_102 = 1;

        DENDRO_102 *= DENDRO_53;
        DENDRO_102 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_103 = 1;

        DENDRO_103 *= DENDRO_55;
        DENDRO_103 *= DENDRO_23;

        // initialize reduction for
        double DENDRO_104 = 1;

        DENDRO_104 *= DENDRO_56;
        DENDRO_104 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_105 = 1;

        DENDRO_105 *= DENDRO_54;
        DENDRO_105 *= DENDRO_21;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_52;
        DENDRO_54 *= DENDRO_21;

        // initialize reduction for
        double DENDRO_106 = 1;

        DENDRO_106 *= DENDRO_53;
        DENDRO_106 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_107 = 1;

        DENDRO_107 *= DENDRO_57;
        DENDRO_107 *= DENDRO_23;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_52;
        DENDRO_23 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_108 = 1;

        DENDRO_108 *= DENDRO_47;
        DENDRO_108 *= DENDRO_15;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_53;
        DENDRO_15 *= DENDRO_14;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_55;
        DENDRO_14 *= DENDRO_21;

        // initialize reduction for
        DENDRO_53 = 1;

        DENDRO_53 *= DENDRO_56;
        DENDRO_53 *= DENDRO_20;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_57;
        DENDRO_20 *= DENDRO_21;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_2;
        DENDRO_21 *= DENDRO_46;
        DENDRO_21 *= DENDRO_6;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_47;
        DENDRO_6 *= DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_2;
        DENDRO_4 *= DENDRO_59;
        DENDRO_4 *= DENDRO_58;

        // initialize reduction for
        DENDRO_46 = 1;

        DENDRO_46 *= DENDRO_60;
        DENDRO_46 *= DENDRO_41;
        DENDRO_46 *= DENDRO_24;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_2;
        DENDRO_55 *= DENDRO_59;
        DENDRO_55 *= DENDRO_61;

        // initialize reduction for
        DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_60;
        DENDRO_56 *= DENDRO_42;
        DENDRO_56 *= DENDRO_24;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_2;
        DENDRO_57 *= DENDRO_59;
        DENDRO_57 *= DENDRO_62;

        // initialize reduction for
        DENDRO_61 = 1;

        DENDRO_61 *= DENDRO_60;
        DENDRO_61 *= DENDRO_42;
        DENDRO_61 *= DENDRO_41;

        // initialize reduction for
        DENDRO_62 = 1;

        DENDRO_62 *= DENDRO_24;
        DENDRO_62 *= DENDRO_63;

        // initialize reduction for
        double DENDRO_109 = 1;

        DENDRO_109 *= DENDRO_41;
        DENDRO_109 *= DENDRO_25;

        // initialize reduction for
        double DENDRO_110 = 1;

        DENDRO_110 *= DENDRO_42;
        DENDRO_110 *= DENDRO_12;

        // initialize reduction for
        double DENDRO_111 = 1;

        DENDRO_111 *= DENDRO_2;
        DENDRO_111 *= DENDRO_59;
        DENDRO_111 *= DENDRO_18;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_60;
        DENDRO_18 *= DENDRO_26;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_59;
        DENDRO_26 *= DENDRO_27;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_60;
        DENDRO_27 *= DENDRO_28;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_2;
        DENDRO_28 *= DENDRO_59;
        DENDRO_28 *= DENDRO_29;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_60;
        DENDRO_29 *= DENDRO_30;
        DENDRO_30         = -1;
        DENDRO_60         = grad2_1_2_gt4;
        double DENDRO_112 = -4.0;
        double DENDRO_113 = grad2_0_2_gt4;
        double DENDRO_114 = grad2_0_1_gt4;
        double DENDRO_115 = grad2_2_2_gt4;
        double DENDRO_116 = -2.0;
        double DENDRO_117 = grad2_1_1_gt4;
        double DENDRO_118 = grad2_0_0_gt4;

        // initialize reduction for
        double DENDRO_119 = 0;

        DENDRO_119 += DENDRO_44;
        DENDRO_119 += DENDRO_43;
        DENDRO_43 = 2.0;

        // initialize reduction for
        DENDRO_44 = 0;

        DENDRO_44 += DENDRO_45;
        DENDRO_44 += DENDRO_50;

        // initialize reduction for
        DENDRO_45 = 0;

        DENDRO_45 += DENDRO_47;
        DENDRO_45 += DENDRO_52;
        DENDRO_47         = gt4[pp];
        DENDRO_50         = grad_2_Gt2;
        DENDRO_52         = gt3[pp];
        double DENDRO_120 = grad_2_Gt1;
        double DENDRO_121 = gt1[pp];
        double DENDRO_122 = grad_2_Gt0;
        double DENDRO_123 = gt5[pp];
        double DENDRO_124 = grad_1_Gt2;
        double DENDRO_125 = grad_1_Gt1;
        double DENDRO_126 = gt2[pp];
        double DENDRO_127 = grad_1_Gt0;

        // initialize reduction for
        double DENDRO_128 = 0;

        DENDRO_128 += DENDRO_33;
        DENDRO_128 += DENDRO_32;
        DENDRO_128 += DENDRO_31;
        DENDRO_31 = 4;

        // initialize reduction for
        DENDRO_32 = 0;

        DENDRO_32 += DENDRO_36;
        DENDRO_32 += DENDRO_35;
        DENDRO_32 += DENDRO_34;

        // initialize reduction for
        DENDRO_33 = 0;

        DENDRO_33 += DENDRO_38;
        DENDRO_33 += DENDRO_37;

        // initialize reduction for
        DENDRO_34 = 0;

        DENDRO_34 += DENDRO_64;
        DENDRO_34 += DENDRO_40;
        DENDRO_34 += DENDRO_39;

        // initialize reduction for
        DENDRO_35 = 0;

        DENDRO_35 += DENDRO_65;
        DENDRO_35 += DENDRO_40;
        DENDRO_35 += DENDRO_39;

        // initialize reduction for
        DENDRO_36 = 0;

        DENDRO_36 += DENDRO_68;
        DENDRO_36 += DENDRO_67;
        DENDRO_36 += DENDRO_66;

        // initialize reduction for
        DENDRO_37 = 0;

        DENDRO_37 += DENDRO_71;
        DENDRO_37 += DENDRO_70;
        DENDRO_37 += DENDRO_69;

        // initialize reduction for
        DENDRO_38 = 0;

        DENDRO_38 += DENDRO_74;
        DENDRO_38 += DENDRO_73;
        DENDRO_38 += DENDRO_72;

        // initialize reduction for
        DENDRO_39 = 0;

        DENDRO_39 += DENDRO_76;
        DENDRO_39 += DENDRO_75;

        // initialize reduction for
        DENDRO_40 = 0;

        DENDRO_40 += DENDRO_79;
        DENDRO_40 += DENDRO_78;
        DENDRO_40 += DENDRO_77;

        // initialize reduction for
        DENDRO_64 = 0;

        DENDRO_64 += DENDRO_82;
        DENDRO_64 += DENDRO_81;
        DENDRO_64 += DENDRO_80;

        // initialize reduction for
        DENDRO_65 = 0;

        DENDRO_65 += DENDRO_84;
        DENDRO_65 += DENDRO_83;

        // initialize reduction for
        DENDRO_66 = 0;

        DENDRO_66 += DENDRO_86;
        DENDRO_66 += DENDRO_49;
        DENDRO_66 += DENDRO_85;

        // initialize reduction for
        DENDRO_67 = 0;

        DENDRO_67 += DENDRO_87;
        DENDRO_67 += DENDRO_49;
        DENDRO_67 += DENDRO_85;

        // initialize reduction for
        DENDRO_49 = 0;

        DENDRO_49 += DENDRO_89;
        DENDRO_49 += DENDRO_88;
        DENDRO_49 += DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 0;

        DENDRO_8 += DENDRO_92;
        DENDRO_8 += DENDRO_91;
        DENDRO_8 += DENDRO_90;

        // initialize reduction for
        DENDRO_68 = 0;

        DENDRO_68 += DENDRO_22;
        DENDRO_68 += DENDRO_93;
        DENDRO_68 += DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 0;

        DENDRO_16 += DENDRO_95;
        DENDRO_16 += DENDRO_94;
        DENDRO_16 += DENDRO_48;

        // initialize reduction for
        DENDRO_22 = 0;

        DENDRO_22 += DENDRO_10;
        DENDRO_22 += DENDRO_97;
        DENDRO_22 += DENDRO_96;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_98;
        DENDRO_10 += DENDRO_97;
        DENDRO_10 += DENDRO_96;

        // initialize reduction for
        DENDRO_48 = 0;

        DENDRO_48 += DENDRO_51;
        DENDRO_48 += DENDRO_100;
        DENDRO_48 += DENDRO_99;

        // initialize reduction for
        DENDRO_51 = 0;

        DENDRO_51 += DENDRO_102;
        DENDRO_51 += DENDRO_17;
        DENDRO_51 += DENDRO_101;

        // initialize reduction for
        DENDRO_17 = 0;

        DENDRO_17 += DENDRO_105;
        DENDRO_17 += DENDRO_104;
        DENDRO_17 += DENDRO_103;

        // initialize reduction for
        DENDRO_69 = 0;

        DENDRO_69 += DENDRO_107;
        DENDRO_69 += DENDRO_106;
        DENDRO_69 += DENDRO_54;

        // initialize reduction for
        DENDRO_54 = 0;

        DENDRO_54 += DENDRO_15;
        DENDRO_54 += DENDRO_108;
        DENDRO_54 += DENDRO_23;

        // initialize reduction for
        DENDRO_15 = 0;

        DENDRO_15 += DENDRO_20;
        DENDRO_15 += DENDRO_53;
        DENDRO_15 += DENDRO_14;

        // initialize reduction for
        DENDRO_14 = 0;

        DENDRO_14 += DENDRO_6;
        DENDRO_14 += DENDRO_21;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_46;
        DENDRO_6 += DENDRO_4;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_56;
        DENDRO_4 += DENDRO_55;

        // initialize reduction for
        DENDRO_20 = 0;

        DENDRO_20 += DENDRO_61;
        DENDRO_20 += DENDRO_57;

        // initialize reduction for
        DENDRO_21 = 0;

        DENDRO_21 += DENDRO_110;
        DENDRO_21 += DENDRO_109;
        DENDRO_21 += DENDRO_62;
        DENDRO_23 = -2;

        // initialize reduction for
        DENDRO_46 = 0;

        DENDRO_46 += DENDRO_18;
        DENDRO_46 += DENDRO_111;

        // initialize reduction for
        DENDRO_18 = 0;

        DENDRO_18 += DENDRO_27;
        DENDRO_18 += DENDRO_26;

        // initialize reduction for
        DENDRO_26 = 0;

        DENDRO_26 += DENDRO_29;
        DENDRO_26 += DENDRO_28;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_30;
        DENDRO_27 *= DENDRO_1;
        DENDRO_27 *= DENDRO_24;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_30;
        DENDRO_1 *= DENDRO_13;
        DENDRO_1 *= DENDRO_41;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_30;
        DENDRO_13 *= DENDRO_19;
        DENDRO_13 *= DENDRO_42;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_112;
        DENDRO_19 *= DENDRO_0;
        DENDRO_19 *= DENDRO_60;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_112;
        DENDRO_28 *= DENDRO_3;
        DENDRO_28 *= DENDRO_113;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_112;
        DENDRO_29 *= DENDRO_5;
        DENDRO_29 *= DENDRO_114;

        // initialize reduction for
        DENDRO_42 = 1;

        DENDRO_42 *= DENDRO_116;
        DENDRO_42 *= DENDRO_7;
        DENDRO_42 *= DENDRO_115;

        // initialize reduction for
        DENDRO_53 = 1;

        DENDRO_53 *= DENDRO_116;
        DENDRO_53 *= DENDRO_9;
        DENDRO_53 *= DENDRO_117;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_116;
        DENDRO_55 *= DENDRO_11;
        DENDRO_55 *= DENDRO_118;

        // initialize reduction for
        DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_43;
        DENDRO_56 *= DENDRO_119;
        DENDRO_56 *= DENDRO_63;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_43;
        DENDRO_57 *= DENDRO_44;
        DENDRO_57 *= DENDRO_25;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_43;
        DENDRO_25 *= DENDRO_45;
        DENDRO_25 *= DENDRO_12;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_43;
        DENDRO_12 *= DENDRO_50;
        DENDRO_12 *= DENDRO_47;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_43;
        DENDRO_44 *= DENDRO_120;
        DENDRO_44 *= DENDRO_52;

        // initialize reduction for
        DENDRO_45 = 1;

        DENDRO_45 *= DENDRO_43;
        DENDRO_45 *= DENDRO_122;
        DENDRO_45 *= DENDRO_121;

        // initialize reduction for
        DENDRO_50 = 1;

        DENDRO_50 *= DENDRO_43;
        DENDRO_50 *= DENDRO_124;
        DENDRO_50 *= DENDRO_123;

        // initialize reduction for
        DENDRO_52 = 1;

        DENDRO_52 *= DENDRO_43;
        DENDRO_52 *= DENDRO_125;
        DENDRO_52 *= DENDRO_47;

        // initialize reduction for
        DENDRO_60 = 1;

        DENDRO_60 *= DENDRO_43;
        DENDRO_60 *= DENDRO_127;
        DENDRO_60 *= DENDRO_126;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_31;
        DENDRO_43 *= DENDRO_7;
        DENDRO_43 *= DENDRO_128;

        // initialize reduction for
        DENDRO_61 = 1;

        DENDRO_61 *= DENDRO_31;
        DENDRO_61 *= DENDRO_7;
        DENDRO_61 *= DENDRO_32;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_31;
        DENDRO_32 *= DENDRO_7;
        DENDRO_32 *= DENDRO_33;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_31;
        DENDRO_33 *= DENDRO_0;
        DENDRO_33 *= DENDRO_34;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_31;
        DENDRO_34 *= DENDRO_0;
        DENDRO_34 *= DENDRO_35;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_31;
        DENDRO_35 *= DENDRO_0;
        DENDRO_35 *= DENDRO_36;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_31;
        DENDRO_36 *= DENDRO_0;
        DENDRO_36 *= DENDRO_37;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_31;
        DENDRO_37 *= DENDRO_0;
        DENDRO_37 *= DENDRO_38;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_31;
        DENDRO_38 *= DENDRO_0;
        DENDRO_38 *= DENDRO_39;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_31;
        DENDRO_39 *= DENDRO_9;
        DENDRO_39 *= DENDRO_40;

        // initialize reduction for
        DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_31;
        DENDRO_40 *= DENDRO_9;
        DENDRO_40 *= DENDRO_64;

        // initialize reduction for
        DENDRO_62 = 1;

        DENDRO_62 *= DENDRO_31;
        DENDRO_62 *= DENDRO_9;
        DENDRO_62 *= DENDRO_65;

        // initialize reduction for
        DENDRO_63 = 1;

        DENDRO_63 *= DENDRO_31;
        DENDRO_63 *= DENDRO_3;
        DENDRO_63 *= DENDRO_66;

        // initialize reduction for
        DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_31;
        DENDRO_64 *= DENDRO_3;
        DENDRO_64 *= DENDRO_67;

        // initialize reduction for
        DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_31;
        DENDRO_65 *= DENDRO_3;
        DENDRO_65 *= DENDRO_49;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_31;
        DENDRO_49 *= DENDRO_3;
        DENDRO_49 *= DENDRO_8;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_31;
        DENDRO_8 *= DENDRO_3;
        DENDRO_8 *= DENDRO_68;

        // initialize reduction for
        DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_31;
        DENDRO_66 *= DENDRO_3;
        DENDRO_66 *= DENDRO_16;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_31;
        DENDRO_16 *= DENDRO_5;
        DENDRO_16 *= DENDRO_22;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_31;
        DENDRO_22 *= DENDRO_5;
        DENDRO_22 *= DENDRO_10;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_31;
        DENDRO_10 *= DENDRO_5;
        DENDRO_10 *= DENDRO_48;

        // initialize reduction for
        DENDRO_48 = 1;

        DENDRO_48 *= DENDRO_31;
        DENDRO_48 *= DENDRO_5;
        DENDRO_48 *= DENDRO_51;

        // initialize reduction for
        DENDRO_51 = 1;

        DENDRO_51 *= DENDRO_31;
        DENDRO_51 *= DENDRO_5;
        DENDRO_51 *= DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_31;
        DENDRO_17 *= DENDRO_5;
        DENDRO_17 *= DENDRO_69;

        // initialize reduction for
        DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_31;
        DENDRO_67 *= DENDRO_11;
        DENDRO_67 *= DENDRO_54;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_31;
        DENDRO_54 *= DENDRO_11;
        DENDRO_54 *= DENDRO_15;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_31;
        DENDRO_15 *= DENDRO_11;
        DENDRO_15 *= DENDRO_14;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_2;
        DENDRO_14 *= DENDRO_0;
        DENDRO_14 *= DENDRO_6;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_4;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_2;
        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_20;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_23;
        DENDRO_4 *= DENDRO_59;
        DENDRO_4 *= DENDRO_21;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_7;
        DENDRO_5 *= DENDRO_46;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_9;
        DENDRO_6 *= DENDRO_18;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_11;
        DENDRO_7 *= DENDRO_26;

        // initialize reduction for
        DENDRO_9 = 0;

        DENDRO_9 += DENDRO_58;
        DENDRO_9 += DENDRO_13;
        DENDRO_9 += DENDRO_1;
        DENDRO_9 += DENDRO_27;

        // initialize reduction for
        DENDRO_1 = 0;

        DENDRO_1 += DENDRO_15;
        DENDRO_1 += DENDRO_54;
        DENDRO_1 += DENDRO_67;
        DENDRO_1 += DENDRO_17;
        DENDRO_1 += DENDRO_51;
        DENDRO_1 += DENDRO_48;
        DENDRO_1 += DENDRO_10;
        DENDRO_1 += DENDRO_22;
        DENDRO_1 += DENDRO_16;
        DENDRO_1 += DENDRO_66;
        DENDRO_1 += DENDRO_8;
        DENDRO_1 += DENDRO_49;
        DENDRO_1 += DENDRO_65;
        DENDRO_1 += DENDRO_64;
        DENDRO_1 += DENDRO_63;
        DENDRO_1 += DENDRO_62;
        DENDRO_1 += DENDRO_40;
        DENDRO_1 += DENDRO_39;
        DENDRO_1 += DENDRO_38;
        DENDRO_1 += DENDRO_37;
        DENDRO_1 += DENDRO_36;
        DENDRO_1 += DENDRO_35;
        DENDRO_1 += DENDRO_34;
        DENDRO_1 += DENDRO_33;
        DENDRO_1 += DENDRO_32;
        DENDRO_1 += DENDRO_61;
        DENDRO_1 += DENDRO_43;
        DENDRO_1 += DENDRO_60;
        DENDRO_1 += DENDRO_52;
        DENDRO_1 += DENDRO_50;
        DENDRO_1 += DENDRO_45;
        DENDRO_1 += DENDRO_44;
        DENDRO_1 += DENDRO_12;
        DENDRO_1 += DENDRO_25;
        DENDRO_1 += DENDRO_57;
        DENDRO_1 += DENDRO_56;
        DENDRO_1 += DENDRO_55;
        DENDRO_1 += DENDRO_53;
        DENDRO_1 += DENDRO_42;
        DENDRO_1 += DENDRO_29;
        DENDRO_1 += DENDRO_28;
        DENDRO_1 += DENDRO_19;
        DENDRO_8 = DENDRO_59;

        DENDRO_8 *= DENDRO_59;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_7;
        DENDRO_10 += DENDRO_6;
        DENDRO_10 += DENDRO_5;
        DENDRO_10 += DENDRO_4;
        DENDRO_10 += DENDRO_3;
        DENDRO_10 += DENDRO_0;
        DENDRO_10 += DENDRO_14;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_59;
        DENDRO_0 *= DENDRO_9;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_30;
        DENDRO_2 *= DENDRO_41;
        DENDRO_2 *= DENDRO_24;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_8;
        DENDRO_3 *= DENDRO_1;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_47;
        DENDRO_1 *= DENDRO_10;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_2;
        DENDRO_4 += DENDRO_0;
        DENDRO_0 = DENDRO_59;

        DENDRO_0 *= DENDRO_59;
        DENDRO_0 = 1 / DENDRO_0;
        DENDRO_1 = 1.0 / 4.0;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_0;
        DENDRO_2 *= DENDRO_4;
        DENDRO_RIJ4 = DENDRO_2;
    }
    {
        double DENDRO_0  = DENDRO_igt4;
        double DENDRO_1  = DENDRO_C2_k2_4;
        double DENDRO_2  = 2;
        double DENDRO_3  = DENDRO_igt2;
        double DENDRO_4  = DENDRO_C2_k2_2;
        double DENDRO_5  = DENDRO_igt1;
        double DENDRO_6  = DENDRO_C2_k2_1;
        double DENDRO_7  = DENDRO_igt5;
        double DENDRO_8  = DENDRO_C2_k2_5;
        double DENDRO_9  = DENDRO_igt3;
        double DENDRO_10 = DENDRO_C2_k2_3;
        double DENDRO_11 = DENDRO_igt0;
        double DENDRO_12 = DENDRO_C2_k2_0;
        double DENDRO_13 = DENDRO_C2_k1_4;
        double DENDRO_14 = DENDRO_C2_k1_2;
        double DENDRO_15 = DENDRO_C2_k1_1;
        double DENDRO_16 = DENDRO_C2_k1_5;
        double DENDRO_17 = DENDRO_C2_k1_3;
        double DENDRO_18 = DENDRO_C2_k1_0;
        double DENDRO_19 = DENDRO_C2_k0_4;
        double DENDRO_20 = DENDRO_C2_k0_2;
        double DENDRO_21 = DENDRO_C2_k0_1;
        double DENDRO_22 = DENDRO_C2_k0_5;
        double DENDRO_23 = DENDRO_C2_k0_3;
        double DENDRO_24 = DENDRO_C2_k0_0;

        // initialize reduction for
        double DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_1;
        DENDRO_25 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_4;
        DENDRO_26 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_6;
        DENDRO_27 *= DENDRO_5;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_8;
        DENDRO_6 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_10;
        DENDRO_28 *= DENDRO_9;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_12;
        DENDRO_10 *= DENDRO_11;

        // initialize reduction for
        DENDRO_12 = 1;

        DENDRO_12 *= DENDRO_2;
        DENDRO_12 *= DENDRO_13;
        DENDRO_12 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_2;
        DENDRO_29 *= DENDRO_14;
        DENDRO_29 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_2;
        DENDRO_30 *= DENDRO_15;
        DENDRO_30 *= DENDRO_5;

        // initialize reduction for
        DENDRO_15 = 1;

        DENDRO_15 *= DENDRO_16;
        DENDRO_15 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_17;
        DENDRO_31 *= DENDRO_9;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_18;
        DENDRO_17 *= DENDRO_11;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_2;
        DENDRO_18 *= DENDRO_19;
        DENDRO_18 *= DENDRO_0;

        // initialize reduction for
        double DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_20;
        DENDRO_32 *= DENDRO_3;

        // initialize reduction for
        double DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_2;
        DENDRO_33 *= DENDRO_21;
        DENDRO_33 *= DENDRO_5;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_22;
        DENDRO_21 *= DENDRO_7;

        // initialize reduction for
        double DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_23;
        DENDRO_34 *= DENDRO_9;

        // initialize reduction for
        DENDRO_23 = 1;

        DENDRO_23 *= DENDRO_24;
        DENDRO_23 *= DENDRO_11;
        DENDRO_24        = grad_2_chi;
        double DENDRO_35 = grad_1_chi;
        double DENDRO_36 = grad_0_chi;
        double DENDRO_37 = DENDRO_C1_k2_4;
        double DENDRO_38 = DENDRO_C1_k2_3;
        double DENDRO_39 = DENDRO_C1_k2_1;
        double DENDRO_40 = DENDRO_C1_k2_2;
        double DENDRO_41 = DENDRO_C1_k2_0;
        double DENDRO_42 = DENDRO_C1_k2_5;
        double DENDRO_43 = DENDRO_C1_k1_5;
        double DENDRO_44 = DENDRO_C1_k1_4;
        double DENDRO_45 = DENDRO_C1_k0_5;
        double DENDRO_46 = DENDRO_C1_k0_4;
        double DENDRO_47 = DENDRO_C1_k1_2;
        double DENDRO_48 = DENDRO_C1_k0_2;
        double DENDRO_49 = grad2_1_2_chi;
        double DENDRO_50 = chi[pp];
        double DENDRO_51 = -3;
        double DENDRO_52 = grad2_0_2_chi;
        double DENDRO_53 = grad2_0_1_chi;

        // initialize reduction for
        double DENDRO_54 = 0;

        DENDRO_54 += DENDRO_10;
        DENDRO_54 += DENDRO_28;
        DENDRO_54 += DENDRO_6;
        DENDRO_54 += DENDRO_27;
        DENDRO_54 += DENDRO_26;
        DENDRO_54 += DENDRO_25;

        // initialize reduction for
        DENDRO_6 = 0;

        DENDRO_6 += DENDRO_17;
        DENDRO_6 += DENDRO_31;
        DENDRO_6 += DENDRO_15;
        DENDRO_6 += DENDRO_30;
        DENDRO_6 += DENDRO_29;
        DENDRO_6 += DENDRO_12;

        // initialize reduction for
        DENDRO_10 = 0;

        DENDRO_10 += DENDRO_23;
        DENDRO_10 += DENDRO_34;
        DENDRO_10 += DENDRO_21;
        DENDRO_10 += DENDRO_33;
        DENDRO_10 += DENDRO_32;
        DENDRO_10 += DENDRO_18;
        DENDRO_12 = grad2_2_2_chi;
        DENDRO_15 = DENDRO_24;

        DENDRO_15 *= DENDRO_24;
        DENDRO_17 = grad2_1_1_chi;
        DENDRO_18 = DENDRO_35;

        DENDRO_18 *= DENDRO_35;
        DENDRO_21 = grad2_0_0_chi;
        DENDRO_23 = DENDRO_36;

        DENDRO_23 *= DENDRO_36;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_37;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_2;
        DENDRO_26 *= DENDRO_38;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_2;
        DENDRO_27 *= DENDRO_39;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_2;
        DENDRO_28 *= DENDRO_40;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_2;
        DENDRO_29 *= DENDRO_41;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_2;
        DENDRO_30 *= DENDRO_37;
        DENDRO_30 *= DENDRO_8;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_42;
        DENDRO_31 *= DENDRO_1;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_2;
        DENDRO_32 *= DENDRO_42;
        DENDRO_32 *= DENDRO_1;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_37;
        DENDRO_33 *= DENDRO_8;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_2;
        DENDRO_34 *= DENDRO_38;
        DENDRO_34 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_43;
        DENDRO_55 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_2;
        DENDRO_56 *= DENDRO_37;
        DENDRO_56 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_44;
        DENDRO_57 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_2;
        DENDRO_58 *= DENDRO_39;
        DENDRO_58 *= DENDRO_22;

        // initialize reduction for
        double DENDRO_59 = 1;

        DENDRO_59 *= DENDRO_45;
        DENDRO_59 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_60 = 1;

        DENDRO_60 *= DENDRO_2;
        DENDRO_60 *= DENDRO_40;
        DENDRO_60 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_61 = 1;

        DENDRO_61 *= DENDRO_46;
        DENDRO_61 *= DENDRO_22;

        // initialize reduction for
        double DENDRO_62 = 1;

        DENDRO_62 *= DENDRO_2;
        DENDRO_62 *= DENDRO_40;
        DENDRO_62 *= DENDRO_8;

        // initialize reduction for
        double DENDRO_63 = 1;

        DENDRO_63 *= DENDRO_42;
        DENDRO_63 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_64 = 1;

        DENDRO_64 *= DENDRO_2;
        DENDRO_64 *= DENDRO_42;
        DENDRO_64 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_65 = 1;

        DENDRO_65 *= DENDRO_40;
        DENDRO_65 *= DENDRO_8;

        // initialize reduction for
        double DENDRO_66 = 1;

        DENDRO_66 *= DENDRO_2;
        DENDRO_66 *= DENDRO_39;
        DENDRO_66 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_67 = 1;

        DENDRO_67 *= DENDRO_43;
        DENDRO_67 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_68 = 1;

        DENDRO_68 *= DENDRO_2;
        DENDRO_68 *= DENDRO_37;
        DENDRO_68 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_69 = 1;

        DENDRO_69 *= DENDRO_47;
        DENDRO_69 *= DENDRO_16;

        // initialize reduction for
        double DENDRO_70 = 1;

        DENDRO_70 *= DENDRO_2;
        DENDRO_70 *= DENDRO_41;
        DENDRO_70 *= DENDRO_22;

        // initialize reduction for
        double DENDRO_71 = 1;

        DENDRO_71 *= DENDRO_45;
        DENDRO_71 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_72 = 1;

        DENDRO_72 *= DENDRO_2;
        DENDRO_72 *= DENDRO_40;
        DENDRO_72 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_73 = 1;

        DENDRO_73 *= DENDRO_48;
        DENDRO_73 *= DENDRO_22;

        // initialize reduction for
        double DENDRO_74 = 1;

        DENDRO_74 *= DENDRO_2;
        DENDRO_74 *= DENDRO_40;
        DENDRO_74 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_75 = 1;

        DENDRO_75 *= DENDRO_37;
        DENDRO_75 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_76 = 1;

        DENDRO_76 *= DENDRO_2;
        DENDRO_76 *= DENDRO_37;
        DENDRO_76 *= DENDRO_4;

        // initialize reduction for
        double DENDRO_77 = 1;

        DENDRO_77 *= DENDRO_40;
        DENDRO_77 *= DENDRO_1;

        // initialize reduction for
        double DENDRO_78 = 1;

        DENDRO_78 *= DENDRO_2;
        DENDRO_78 *= DENDRO_39;
        DENDRO_78 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_79 = 1;

        DENDRO_79 *= DENDRO_44;
        DENDRO_79 *= DENDRO_14;

        // initialize reduction for
        double DENDRO_80 = 1;

        DENDRO_80 *= DENDRO_2;
        DENDRO_80 *= DENDRO_38;
        DENDRO_80 *= DENDRO_14;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_47;
        DENDRO_38 *= DENDRO_13;

        // initialize reduction for
        double DENDRO_81 = 1;

        DENDRO_81 *= DENDRO_2;
        DENDRO_81 *= DENDRO_41;
        DENDRO_81 *= DENDRO_19;

        // initialize reduction for
        DENDRO_41 = 1;

        DENDRO_41 *= DENDRO_46;
        DENDRO_41 *= DENDRO_20;

        // initialize reduction for
        double DENDRO_82 = 1;

        DENDRO_82 *= DENDRO_2;
        DENDRO_82 *= DENDRO_39;
        DENDRO_82 *= DENDRO_20;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_48;
        DENDRO_39 *= DENDRO_19;

        // initialize reduction for
        double DENDRO_83 = 1;

        DENDRO_83 *= DENDRO_2;
        DENDRO_83 *= DENDRO_50;
        DENDRO_83 *= DENDRO_49;

        // initialize reduction for
        DENDRO_49 = 1;

        DENDRO_49 *= DENDRO_51;
        DENDRO_49 *= DENDRO_35;
        DENDRO_49 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_84 = 1;

        DENDRO_84 *= DENDRO_2;
        DENDRO_84 *= DENDRO_50;
        DENDRO_84 *= DENDRO_52;

        // initialize reduction for
        DENDRO_52 = 1;

        DENDRO_52 *= DENDRO_51;
        DENDRO_52 *= DENDRO_36;
        DENDRO_52 *= DENDRO_24;

        // initialize reduction for
        double DENDRO_85 = 1;

        DENDRO_85 *= DENDRO_2;
        DENDRO_85 *= DENDRO_50;
        DENDRO_85 *= DENDRO_53;

        // initialize reduction for
        DENDRO_53 = 1;

        DENDRO_53 *= DENDRO_51;
        DENDRO_53 *= DENDRO_36;
        DENDRO_53 *= DENDRO_35;

        // initialize reduction for
        double DENDRO_86 = 1;

        DENDRO_86 *= DENDRO_24;
        DENDRO_86 *= DENDRO_54;

        // initialize reduction for
        double DENDRO_87 = 1;

        DENDRO_87 *= DENDRO_35;
        DENDRO_87 *= DENDRO_6;

        // initialize reduction for
        double DENDRO_88 = 1;

        DENDRO_88 *= DENDRO_36;
        DENDRO_88 *= DENDRO_10;

        // initialize reduction for
        double DENDRO_89 = 1;

        DENDRO_89 *= DENDRO_2;
        DENDRO_89 *= DENDRO_50;
        DENDRO_89 *= DENDRO_12;

        // initialize reduction for
        double DENDRO_90 = 1;

        DENDRO_90 *= DENDRO_51;
        DENDRO_90 *= DENDRO_15;

        // initialize reduction for
        double DENDRO_91 = 1;

        DENDRO_91 *= DENDRO_2;
        DENDRO_91 *= DENDRO_50;
        DENDRO_91 *= DENDRO_17;

        // initialize reduction for
        DENDRO_17 = 1;

        DENDRO_17 *= DENDRO_51;
        DENDRO_17 *= DENDRO_18;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_2;
        DENDRO_18 *= DENDRO_50;
        DENDRO_18 *= DENDRO_21;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_51;
        DENDRO_21 *= DENDRO_23;
        DENDRO_23        = -1;
        DENDRO_51        = 12;

        // initialize reduction for
        double DENDRO_92 = 0;

        DENDRO_92 += DENDRO_43;
        DENDRO_92 += DENDRO_25;
        DENDRO_25 = 4;

        // initialize reduction for
        DENDRO_43 = 0;

        DENDRO_43 += DENDRO_44;
        DENDRO_43 += DENDRO_26;

        // initialize reduction for
        DENDRO_26 = 0;

        DENDRO_26 += DENDRO_47;
        DENDRO_26 += DENDRO_27;

        // initialize reduction for
        DENDRO_44 = 0;

        DENDRO_44 += DENDRO_45;
        DENDRO_44 += DENDRO_28;

        // initialize reduction for
        DENDRO_28 = 0;

        DENDRO_28 += DENDRO_46;
        DENDRO_28 += DENDRO_27;

        // initialize reduction for
        DENDRO_27 = 0;

        DENDRO_27 += DENDRO_48;
        DENDRO_27 += DENDRO_29;
        DENDRO_29         = grad2_1_2_gt5;
        DENDRO_45         = -4.0;
        DENDRO_46         = grad2_0_2_gt5;
        DENDRO_47         = grad2_0_1_gt5;
        DENDRO_48         = grad2_2_2_gt5;
        double DENDRO_93  = -2.0;
        double DENDRO_94  = grad2_1_1_gt5;
        double DENDRO_95  = grad2_0_0_gt5;
        double DENDRO_96  = gt5[pp];
        double DENDRO_97  = grad_2_Gt2;
        double DENDRO_98  = 4.0;
        double DENDRO_99  = gt4[pp];
        double DENDRO_100 = grad_2_Gt1;
        double DENDRO_101 = gt2[pp];
        double DENDRO_102 = grad_2_Gt0;

        // initialize reduction for
        double DENDRO_103 = 0;

        DENDRO_103 += DENDRO_31;
        DENDRO_103 += DENDRO_30;

        // initialize reduction for
        DENDRO_30 = 0;

        DENDRO_30 += DENDRO_33;
        DENDRO_30 += DENDRO_32;

        // initialize reduction for
        DENDRO_31 = 0;

        DENDRO_31 += DENDRO_55;
        DENDRO_31 += DENDRO_34;

        // initialize reduction for
        DENDRO_32 = 0;

        DENDRO_32 += DENDRO_57;
        DENDRO_32 += DENDRO_56;

        // initialize reduction for
        DENDRO_33 = 0;

        DENDRO_33 += DENDRO_59;
        DENDRO_33 += DENDRO_58;

        // initialize reduction for
        DENDRO_34 = 0;

        DENDRO_34 += DENDRO_61;
        DENDRO_34 += DENDRO_60;

        // initialize reduction for
        DENDRO_55 = 0;

        DENDRO_55 += DENDRO_63;
        DENDRO_55 += DENDRO_62;

        // initialize reduction for
        DENDRO_56 = 0;

        DENDRO_56 += DENDRO_65;
        DENDRO_56 += DENDRO_64;

        // initialize reduction for
        DENDRO_57 = 0;

        DENDRO_57 += DENDRO_67;
        DENDRO_57 += DENDRO_66;

        // initialize reduction for
        DENDRO_58 = 0;

        DENDRO_58 += DENDRO_69;
        DENDRO_58 += DENDRO_68;

        // initialize reduction for
        DENDRO_59 = 0;

        DENDRO_59 += DENDRO_71;
        DENDRO_59 += DENDRO_70;

        // initialize reduction for
        DENDRO_60 = 0;

        DENDRO_60 += DENDRO_73;
        DENDRO_60 += DENDRO_72;

        // initialize reduction for
        DENDRO_61 = 0;

        DENDRO_61 += DENDRO_75;
        DENDRO_61 += DENDRO_74;

        // initialize reduction for
        DENDRO_62 = 0;

        DENDRO_62 += DENDRO_77;
        DENDRO_62 += DENDRO_76;

        // initialize reduction for
        DENDRO_63 = 0;

        DENDRO_63 += DENDRO_79;
        DENDRO_63 += DENDRO_78;

        // initialize reduction for
        DENDRO_64 = 0;

        DENDRO_64 += DENDRO_38;
        DENDRO_64 += DENDRO_80;

        // initialize reduction for
        DENDRO_38 = 0;

        DENDRO_38 += DENDRO_41;
        DENDRO_38 += DENDRO_81;

        // initialize reduction for
        DENDRO_41 = 0;

        DENDRO_41 += DENDRO_39;
        DENDRO_41 += DENDRO_82;

        // initialize reduction for
        DENDRO_39 = 0;

        DENDRO_39 += DENDRO_49;
        DENDRO_39 += DENDRO_83;

        // initialize reduction for
        DENDRO_49 = 0;

        DENDRO_49 += DENDRO_52;
        DENDRO_49 += DENDRO_84;

        // initialize reduction for
        DENDRO_52 = 0;

        DENDRO_52 += DENDRO_53;
        DENDRO_52 += DENDRO_85;

        // initialize reduction for
        DENDRO_53 = 0;

        DENDRO_53 += DENDRO_88;
        DENDRO_53 += DENDRO_87;
        DENDRO_53 += DENDRO_86;
        DENDRO_65 = -2;

        // initialize reduction for
        DENDRO_66 = 0;

        DENDRO_66 += DENDRO_90;
        DENDRO_66 += DENDRO_89;

        // initialize reduction for
        DENDRO_67 = 0;

        DENDRO_67 += DENDRO_17;
        DENDRO_67 += DENDRO_91;

        // initialize reduction for
        DENDRO_17 = 0;

        DENDRO_17 += DENDRO_21;
        DENDRO_17 += DENDRO_18;

        // initialize reduction for
        DENDRO_18 = 1;

        DENDRO_18 *= DENDRO_23;
        DENDRO_18 *= DENDRO_8;
        DENDRO_18 *= DENDRO_24;

        // initialize reduction for
        DENDRO_21 = 1;

        DENDRO_21 *= DENDRO_23;
        DENDRO_21 *= DENDRO_16;
        DENDRO_21 *= DENDRO_35;

        // initialize reduction for
        DENDRO_24 = 1;

        DENDRO_24 *= DENDRO_23;
        DENDRO_24 *= DENDRO_22;
        DENDRO_24 *= DENDRO_36;

        // initialize reduction for
        DENDRO_35 = 1;

        DENDRO_35 *= DENDRO_51;
        DENDRO_35 *= DENDRO_42;
        DENDRO_35 *= DENDRO_8;
        DENDRO_35 *= DENDRO_7;

        // initialize reduction for
        DENDRO_8 = 1;

        DENDRO_8 *= DENDRO_51;
        DENDRO_8 *= DENDRO_37;
        DENDRO_8 *= DENDRO_1;
        DENDRO_8 *= DENDRO_9;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_51;
        DENDRO_1 *= DENDRO_40;
        DENDRO_1 *= DENDRO_4;
        DENDRO_1 *= DENDRO_11;

        // initialize reduction for
        DENDRO_4 = 1;

        DENDRO_4 *= DENDRO_25;
        DENDRO_4 *= DENDRO_16;
        DENDRO_4 *= DENDRO_7;
        DENDRO_4 *= DENDRO_92;

        // initialize reduction for
        DENDRO_16 = 1;

        DENDRO_16 *= DENDRO_25;
        DENDRO_16 *= DENDRO_13;
        DENDRO_16 *= DENDRO_9;
        DENDRO_16 *= DENDRO_43;

        // initialize reduction for
        DENDRO_13 = 1;

        DENDRO_13 *= DENDRO_25;
        DENDRO_13 *= DENDRO_14;
        DENDRO_13 *= DENDRO_11;
        DENDRO_13 *= DENDRO_26;

        // initialize reduction for
        DENDRO_14 = 1;

        DENDRO_14 *= DENDRO_25;
        DENDRO_14 *= DENDRO_22;
        DENDRO_14 *= DENDRO_7;
        DENDRO_14 *= DENDRO_44;

        // initialize reduction for
        DENDRO_22 = 1;

        DENDRO_22 *= DENDRO_25;
        DENDRO_22 *= DENDRO_19;
        DENDRO_22 *= DENDRO_9;
        DENDRO_22 *= DENDRO_28;

        // initialize reduction for
        DENDRO_19 = 1;

        DENDRO_19 *= DENDRO_25;
        DENDRO_19 *= DENDRO_20;
        DENDRO_19 *= DENDRO_11;
        DENDRO_19 *= DENDRO_27;

        // initialize reduction for
        DENDRO_20 = 1;

        DENDRO_20 *= DENDRO_45;
        DENDRO_20 *= DENDRO_0;
        DENDRO_20 *= DENDRO_29;

        // initialize reduction for
        DENDRO_26 = 1;

        DENDRO_26 *= DENDRO_45;
        DENDRO_26 *= DENDRO_3;
        DENDRO_26 *= DENDRO_46;

        // initialize reduction for
        DENDRO_27 = 1;

        DENDRO_27 *= DENDRO_45;
        DENDRO_27 *= DENDRO_5;
        DENDRO_27 *= DENDRO_47;

        // initialize reduction for
        DENDRO_28 = 1;

        DENDRO_28 *= DENDRO_93;
        DENDRO_28 *= DENDRO_7;
        DENDRO_28 *= DENDRO_48;

        // initialize reduction for
        DENDRO_29 = 1;

        DENDRO_29 *= DENDRO_93;
        DENDRO_29 *= DENDRO_9;
        DENDRO_29 *= DENDRO_94;

        // initialize reduction for
        DENDRO_36 = 1;

        DENDRO_36 *= DENDRO_93;
        DENDRO_36 *= DENDRO_11;
        DENDRO_36 *= DENDRO_95;

        // initialize reduction for
        DENDRO_43 = 1;

        DENDRO_43 *= DENDRO_98;
        DENDRO_43 *= DENDRO_97;
        DENDRO_43 *= DENDRO_96;

        // initialize reduction for
        DENDRO_44 = 1;

        DENDRO_44 *= DENDRO_98;
        DENDRO_44 *= DENDRO_100;
        DENDRO_44 *= DENDRO_99;

        // initialize reduction for
        DENDRO_45 = 1;

        DENDRO_45 *= DENDRO_98;
        DENDRO_45 *= DENDRO_102;
        DENDRO_45 *= DENDRO_101;

        // initialize reduction for
        DENDRO_46 = 1;

        DENDRO_46 *= DENDRO_98;
        DENDRO_46 *= DENDRO_42;
        DENDRO_46 *= DENDRO_54;

        // initialize reduction for
        DENDRO_42 = 1;

        DENDRO_42 *= DENDRO_98;
        DENDRO_42 *= DENDRO_37;
        DENDRO_42 *= DENDRO_6;

        // initialize reduction for
        DENDRO_6 = 1;

        DENDRO_6 *= DENDRO_98;
        DENDRO_6 *= DENDRO_40;
        DENDRO_6 *= DENDRO_10;

        // initialize reduction for
        DENDRO_10 = 1;

        DENDRO_10 *= DENDRO_25;
        DENDRO_10 *= DENDRO_0;
        DENDRO_10 *= DENDRO_103;

        // initialize reduction for
        DENDRO_37 = 1;

        DENDRO_37 *= DENDRO_25;
        DENDRO_37 *= DENDRO_0;
        DENDRO_37 *= DENDRO_30;

        // initialize reduction for
        DENDRO_30 = 1;

        DENDRO_30 *= DENDRO_25;
        DENDRO_30 *= DENDRO_0;
        DENDRO_30 *= DENDRO_31;

        // initialize reduction for
        DENDRO_31 = 1;

        DENDRO_31 *= DENDRO_25;
        DENDRO_31 *= DENDRO_0;
        DENDRO_31 *= DENDRO_32;

        // initialize reduction for
        DENDRO_32 = 1;

        DENDRO_32 *= DENDRO_25;
        DENDRO_32 *= DENDRO_0;
        DENDRO_32 *= DENDRO_33;

        // initialize reduction for
        DENDRO_33 = 1;

        DENDRO_33 *= DENDRO_25;
        DENDRO_33 *= DENDRO_0;
        DENDRO_33 *= DENDRO_34;

        // initialize reduction for
        DENDRO_34 = 1;

        DENDRO_34 *= DENDRO_25;
        DENDRO_34 *= DENDRO_3;
        DENDRO_34 *= DENDRO_55;

        // initialize reduction for
        DENDRO_40 = 1;

        DENDRO_40 *= DENDRO_25;
        DENDRO_40 *= DENDRO_3;
        DENDRO_40 *= DENDRO_56;

        // initialize reduction for
        DENDRO_47 = 1;

        DENDRO_47 *= DENDRO_25;
        DENDRO_47 *= DENDRO_3;
        DENDRO_47 *= DENDRO_57;

        // initialize reduction for
        DENDRO_48 = 1;

        DENDRO_48 *= DENDRO_25;
        DENDRO_48 *= DENDRO_3;
        DENDRO_48 *= DENDRO_58;

        // initialize reduction for
        DENDRO_51 = 1;

        DENDRO_51 *= DENDRO_25;
        DENDRO_51 *= DENDRO_3;
        DENDRO_51 *= DENDRO_59;

        // initialize reduction for
        DENDRO_54 = 1;

        DENDRO_54 *= DENDRO_25;
        DENDRO_54 *= DENDRO_3;
        DENDRO_54 *= DENDRO_60;

        // initialize reduction for
        DENDRO_55 = 1;

        DENDRO_55 *= DENDRO_25;
        DENDRO_55 *= DENDRO_5;
        DENDRO_55 *= DENDRO_61;

        // initialize reduction for
        DENDRO_56 = 1;

        DENDRO_56 *= DENDRO_25;
        DENDRO_56 *= DENDRO_5;
        DENDRO_56 *= DENDRO_62;

        // initialize reduction for
        DENDRO_57 = 1;

        DENDRO_57 *= DENDRO_25;
        DENDRO_57 *= DENDRO_5;
        DENDRO_57 *= DENDRO_63;

        // initialize reduction for
        DENDRO_58 = 1;

        DENDRO_58 *= DENDRO_25;
        DENDRO_58 *= DENDRO_5;
        DENDRO_58 *= DENDRO_64;

        // initialize reduction for
        DENDRO_59 = 1;

        DENDRO_59 *= DENDRO_25;
        DENDRO_59 *= DENDRO_5;
        DENDRO_59 *= DENDRO_38;

        // initialize reduction for
        DENDRO_38 = 1;

        DENDRO_38 *= DENDRO_25;
        DENDRO_38 *= DENDRO_5;
        DENDRO_38 *= DENDRO_41;

        // initialize reduction for
        DENDRO_25 = 1;

        DENDRO_25 *= DENDRO_2;
        DENDRO_25 *= DENDRO_0;
        DENDRO_25 *= DENDRO_39;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_3;
        DENDRO_0 *= DENDRO_49;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_2;
        DENDRO_3 *= DENDRO_5;
        DENDRO_3 *= DENDRO_52;

        // initialize reduction for
        DENDRO_5 = 1;

        DENDRO_5 *= DENDRO_65;
        DENDRO_5 *= DENDRO_50;
        DENDRO_5 *= DENDRO_53;

        // initialize reduction for
        DENDRO_39 = 1;

        DENDRO_39 *= DENDRO_7;
        DENDRO_39 *= DENDRO_66;

        // initialize reduction for
        DENDRO_7 = 1;

        DENDRO_7 *= DENDRO_9;
        DENDRO_7 *= DENDRO_67;

        // initialize reduction for
        DENDRO_9 = 1;

        DENDRO_9 *= DENDRO_11;
        DENDRO_9 *= DENDRO_17;

        // initialize reduction for
        DENDRO_11 = 0;

        DENDRO_11 += DENDRO_12;
        DENDRO_11 += DENDRO_24;
        DENDRO_11 += DENDRO_21;
        DENDRO_11 += DENDRO_18;

        // initialize reduction for
        DENDRO_12 = 0;

        DENDRO_12 += DENDRO_38;
        DENDRO_12 += DENDRO_59;
        DENDRO_12 += DENDRO_58;
        DENDRO_12 += DENDRO_57;
        DENDRO_12 += DENDRO_56;
        DENDRO_12 += DENDRO_55;
        DENDRO_12 += DENDRO_54;
        DENDRO_12 += DENDRO_51;
        DENDRO_12 += DENDRO_48;
        DENDRO_12 += DENDRO_47;
        DENDRO_12 += DENDRO_40;
        DENDRO_12 += DENDRO_34;
        DENDRO_12 += DENDRO_33;
        DENDRO_12 += DENDRO_32;
        DENDRO_12 += DENDRO_31;
        DENDRO_12 += DENDRO_30;
        DENDRO_12 += DENDRO_37;
        DENDRO_12 += DENDRO_10;
        DENDRO_12 += DENDRO_6;
        DENDRO_12 += DENDRO_42;
        DENDRO_12 += DENDRO_46;
        DENDRO_12 += DENDRO_45;
        DENDRO_12 += DENDRO_44;
        DENDRO_12 += DENDRO_43;
        DENDRO_12 += DENDRO_36;
        DENDRO_12 += DENDRO_29;
        DENDRO_12 += DENDRO_28;
        DENDRO_12 += DENDRO_27;
        DENDRO_12 += DENDRO_26;
        DENDRO_12 += DENDRO_20;
        DENDRO_12 += DENDRO_19;
        DENDRO_12 += DENDRO_22;
        DENDRO_12 += DENDRO_14;
        DENDRO_12 += DENDRO_13;
        DENDRO_12 += DENDRO_16;
        DENDRO_12 += DENDRO_4;
        DENDRO_12 += DENDRO_1;
        DENDRO_12 += DENDRO_8;
        DENDRO_12 += DENDRO_35;
        DENDRO_1 = DENDRO_50;

        DENDRO_1 *= DENDRO_50;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_9;
        DENDRO_4 += DENDRO_7;
        DENDRO_4 += DENDRO_39;
        DENDRO_4 += DENDRO_5;
        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_0;
        DENDRO_4 += DENDRO_25;

        // initialize reduction for
        DENDRO_0 = 1;

        DENDRO_0 *= DENDRO_2;
        DENDRO_0 *= DENDRO_50;
        DENDRO_0 *= DENDRO_11;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_12;

        // initialize reduction for
        DENDRO_1 = 1;

        DENDRO_1 *= DENDRO_96;
        DENDRO_1 *= DENDRO_4;

        // initialize reduction for
        DENDRO_3 = 1;

        DENDRO_3 *= DENDRO_23;
        DENDRO_3 *= DENDRO_15;

        // initialize reduction for
        DENDRO_4 = 0;

        DENDRO_4 += DENDRO_3;
        DENDRO_4 += DENDRO_1;
        DENDRO_4 += DENDRO_2;
        DENDRO_4 += DENDRO_0;
        DENDRO_0 = DENDRO_50;

        DENDRO_0 *= DENDRO_50;
        DENDRO_0 = 1 / DENDRO_0;
        DENDRO_1 = 1.0 / 4.0;

        // initialize reduction for
        DENDRO_2 = 1;

        DENDRO_2 *= DENDRO_1;
        DENDRO_2 *= DENDRO_0;
        DENDRO_2 *= DENDRO_4;
        DENDRO_RIJ5 = DENDRO_2;
    }
    double DENDRO_0  = alpha[pp];
    double DENDRO_1  = DENDRO_RIJ4;
    double DENDRO_2  = grad_2_alpha;
    double DENDRO_3  = DENDRO_C3_k2_4;
    double DENDRO_4  = grad_1_alpha;
    double DENDRO_5  = DENDRO_C3_k1_4;
    double DENDRO_6  = grad_0_alpha;
    double DENDRO_7  = DENDRO_C3_k0_4;
    double DENDRO_8  = grad2_1_2_alpha;
    double DENDRO_9  = -1;
    double DENDRO_10 = DENDRO_RIJ2;
    double DENDRO_11 = DENDRO_C3_k2_2;
    double DENDRO_12 = DENDRO_C3_k1_2;
    double DENDRO_13 = DENDRO_C3_k0_2;
    double DENDRO_14 = grad2_0_2_alpha;
    double DENDRO_15 = DENDRO_RIJ1;
    double DENDRO_16 = DENDRO_C3_k2_1;
    double DENDRO_17 = DENDRO_C3_k1_1;
    double DENDRO_18 = DENDRO_C3_k0_1;
    double DENDRO_19 = grad2_0_1_alpha;
    double DENDRO_20 = DENDRO_RIJ5;
    double DENDRO_21 = DENDRO_C3_k2_5;
    double DENDRO_22 = DENDRO_C3_k1_5;
    double DENDRO_23 = DENDRO_C3_k0_5;
    double DENDRO_24 = grad2_2_2_alpha;
    double DENDRO_25 = DENDRO_RIJ3;
    double DENDRO_26 = DENDRO_C3_k2_3;
    double DENDRO_27 = DENDRO_C3_k1_3;
    double DENDRO_28 = DENDRO_C3_k0_3;
    double DENDRO_29 = grad2_1_1_alpha;
    double DENDRO_30 = DENDRO_RIJ0;
    double DENDRO_31 = DENDRO_C3_k2_0;
    double DENDRO_32 = DENDRO_C3_k1_0;
    double DENDRO_33 = DENDRO_C3_k0_0;
    double DENDRO_34 = grad2_0_0_alpha;

    // initialize reduction for
    double DENDRO_35 = 1;

    DENDRO_35 *= DENDRO_1;
    DENDRO_35 *= DENDRO_0;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_3;
    DENDRO_1 *= DENDRO_2;

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_5;
    DENDRO_3 *= DENDRO_4;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_7;
    DENDRO_5 *= DENDRO_6;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_9;
    DENDRO_7 *= DENDRO_8;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_10;
    DENDRO_8 *= DENDRO_0;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_11;
    DENDRO_10 *= DENDRO_2;

    // initialize reduction for
    DENDRO_11 = 1;

    DENDRO_11 *= DENDRO_12;
    DENDRO_11 *= DENDRO_4;

    // initialize reduction for
    DENDRO_12 = 1;

    DENDRO_12 *= DENDRO_13;
    DENDRO_12 *= DENDRO_6;

    // initialize reduction for
    DENDRO_13 = 1;

    DENDRO_13 *= DENDRO_9;
    DENDRO_13 *= DENDRO_14;

    // initialize reduction for
    DENDRO_14 = 1;

    DENDRO_14 *= DENDRO_15;
    DENDRO_14 *= DENDRO_0;

    // initialize reduction for
    DENDRO_15 = 1;

    DENDRO_15 *= DENDRO_16;
    DENDRO_15 *= DENDRO_2;

    // initialize reduction for
    DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_17;
    DENDRO_16 *= DENDRO_4;

    // initialize reduction for
    DENDRO_17 = 1;

    DENDRO_17 *= DENDRO_18;
    DENDRO_17 *= DENDRO_6;

    // initialize reduction for
    DENDRO_18 = 1;

    DENDRO_18 *= DENDRO_9;
    DENDRO_18 *= DENDRO_19;

    // initialize reduction for
    DENDRO_19 = 1;

    DENDRO_19 *= DENDRO_20;
    DENDRO_19 *= DENDRO_0;

    // initialize reduction for
    DENDRO_20 = 1;

    DENDRO_20 *= DENDRO_21;
    DENDRO_20 *= DENDRO_2;

    // initialize reduction for
    DENDRO_21 = 1;

    DENDRO_21 *= DENDRO_22;
    DENDRO_21 *= DENDRO_4;

    // initialize reduction for
    DENDRO_22 = 1;

    DENDRO_22 *= DENDRO_23;
    DENDRO_22 *= DENDRO_6;

    // initialize reduction for
    DENDRO_23 = 1;

    DENDRO_23 *= DENDRO_9;
    DENDRO_23 *= DENDRO_24;

    // initialize reduction for
    DENDRO_24 = 1;

    DENDRO_24 *= DENDRO_25;
    DENDRO_24 *= DENDRO_0;

    // initialize reduction for
    DENDRO_25 = 1;

    DENDRO_25 *= DENDRO_26;
    DENDRO_25 *= DENDRO_2;

    // initialize reduction for
    DENDRO_26 = 1;

    DENDRO_26 *= DENDRO_27;
    DENDRO_26 *= DENDRO_4;

    // initialize reduction for
    DENDRO_27 = 1;

    DENDRO_27 *= DENDRO_28;
    DENDRO_27 *= DENDRO_6;

    // initialize reduction for
    DENDRO_28 = 1;

    DENDRO_28 *= DENDRO_9;
    DENDRO_28 *= DENDRO_29;

    // initialize reduction for
    DENDRO_29 = 1;

    DENDRO_29 *= DENDRO_30;
    DENDRO_29 *= DENDRO_0;

    // initialize reduction for
    DENDRO_30 = 1;

    DENDRO_30 *= DENDRO_31;
    DENDRO_30 *= DENDRO_2;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_32;
    DENDRO_2 *= DENDRO_4;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_33;
    DENDRO_4 *= DENDRO_6;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_9;
    DENDRO_6 *= DENDRO_34;

    // initialize reduction for
    DENDRO_9 = 0;

    DENDRO_9 += DENDRO_7;
    DENDRO_9 += DENDRO_5;
    DENDRO_9 += DENDRO_3;
    DENDRO_9 += DENDRO_1;
    DENDRO_9 += DENDRO_35;
    DENDRO_31 = DENDRO_igt4;
    DENDRO_32 = 2;

    // initialize reduction for
    DENDRO_33 = 0;

    DENDRO_33 += DENDRO_13;
    DENDRO_33 += DENDRO_12;
    DENDRO_33 += DENDRO_11;
    DENDRO_33 += DENDRO_10;
    DENDRO_33 += DENDRO_8;
    DENDRO_34        = DENDRO_igt2;

    // initialize reduction for
    double DENDRO_36 = 0;

    DENDRO_36 += DENDRO_18;
    DENDRO_36 += DENDRO_17;
    DENDRO_36 += DENDRO_16;
    DENDRO_36 += DENDRO_15;
    DENDRO_36 += DENDRO_14;
    double DENDRO_37 = DENDRO_igt1;

    // initialize reduction for
    double DENDRO_38 = 0;

    DENDRO_38 += DENDRO_23;
    DENDRO_38 += DENDRO_22;
    DENDRO_38 += DENDRO_21;
    DENDRO_38 += DENDRO_20;
    DENDRO_38 += DENDRO_19;
    double DENDRO_39 = DENDRO_igt5;

    // initialize reduction for
    double DENDRO_40 = 0;

    DENDRO_40 += DENDRO_28;
    DENDRO_40 += DENDRO_27;
    DENDRO_40 += DENDRO_26;
    DENDRO_40 += DENDRO_25;
    DENDRO_40 += DENDRO_24;
    double DENDRO_41 = DENDRO_igt3;

    // initialize reduction for
    double DENDRO_42 = 0;

    DENDRO_42 += DENDRO_6;
    DENDRO_42 += DENDRO_4;
    DENDRO_42 += DENDRO_2;
    DENDRO_42 += DENDRO_30;
    DENDRO_42 += DENDRO_29;
    double DENDRO_43 = DENDRO_igt0;
    double DENDRO_44 = At5[pp];
    double DENDRO_45 = At4[pp];
    double DENDRO_46 = At2[pp];
    double DENDRO_47 = At3[pp];
    double DENDRO_48 = At1[pp];
    double DENDRO_49 = At0[pp];

    // initialize reduction for
    double DENDRO_50 = 1;

    DENDRO_50 *= DENDRO_32;
    DENDRO_50 *= DENDRO_31;
    DENDRO_50 *= DENDRO_9;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_32;
    DENDRO_9 *= DENDRO_34;
    DENDRO_9 *= DENDRO_33;

    // initialize reduction for
    DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_32;
    DENDRO_33 *= DENDRO_37;
    DENDRO_33 *= DENDRO_36;

    // initialize reduction for
    DENDRO_36 = 1;

    DENDRO_36 *= DENDRO_39;
    DENDRO_36 *= DENDRO_38;

    // initialize reduction for
    DENDRO_38 = 1;

    DENDRO_38 *= DENDRO_41;
    DENDRO_38 *= DENDRO_40;

    // initialize reduction for
    DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_43;
    DENDRO_40 *= DENDRO_42;

    // initialize reduction for
    DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_44;
    DENDRO_42 *= DENDRO_39;

    // initialize reduction for
    double DENDRO_51 = 1;

    DENDRO_51 *= DENDRO_45;
    DENDRO_51 *= DENDRO_31;

    // initialize reduction for
    double DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_46;
    DENDRO_52 *= DENDRO_34;

    // initialize reduction for
    double DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_44;
    DENDRO_53 *= DENDRO_31;

    // initialize reduction for
    double DENDRO_54 = 1;

    DENDRO_54 *= DENDRO_45;
    DENDRO_54 *= DENDRO_41;

    // initialize reduction for
    double DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_46;
    DENDRO_55 *= DENDRO_37;

    // initialize reduction for
    double DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_44;
    DENDRO_56 *= DENDRO_34;

    // initialize reduction for
    double DENDRO_57 = 1;

    DENDRO_57 *= DENDRO_45;
    DENDRO_57 *= DENDRO_37;

    // initialize reduction for
    double DENDRO_58 = 1;

    DENDRO_58 *= DENDRO_46;
    DENDRO_58 *= DENDRO_43;

    // initialize reduction for
    double DENDRO_59 = 1;

    DENDRO_59 *= DENDRO_45;
    DENDRO_59 *= DENDRO_39;

    // initialize reduction for
    double DENDRO_60 = 1;

    DENDRO_60 *= DENDRO_47;
    DENDRO_60 *= DENDRO_31;

    // initialize reduction for
    double DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_48;
    DENDRO_61 *= DENDRO_34;

    // initialize reduction for
    double DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_47;
    DENDRO_62 *= DENDRO_41;

    // initialize reduction for
    double DENDRO_63 = 1;

    DENDRO_63 *= DENDRO_48;
    DENDRO_63 *= DENDRO_37;

    // initialize reduction for
    double DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_45;
    DENDRO_64 *= DENDRO_34;

    // initialize reduction for
    double DENDRO_65 = 1;

    DENDRO_65 *= DENDRO_47;
    DENDRO_65 *= DENDRO_37;

    // initialize reduction for
    double DENDRO_66 = 1;

    DENDRO_66 *= DENDRO_48;
    DENDRO_66 *= DENDRO_43;

    // initialize reduction for
    double DENDRO_67 = 1;

    DENDRO_67 *= DENDRO_46;
    DENDRO_67 *= DENDRO_39;

    // initialize reduction for
    DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_48;
    DENDRO_39 *= DENDRO_31;

    // initialize reduction for
    double DENDRO_68 = 1;

    DENDRO_68 *= DENDRO_49;
    DENDRO_68 *= DENDRO_34;

    // initialize reduction for
    DENDRO_34 = 1;

    DENDRO_34 *= DENDRO_46;
    DENDRO_34 *= DENDRO_31;

    // initialize reduction for
    DENDRO_31 = 1;

    DENDRO_31 *= DENDRO_48;
    DENDRO_31 *= DENDRO_41;

    // initialize reduction for
    DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_49;
    DENDRO_41 *= DENDRO_37;

    // initialize reduction for
    DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_49;
    DENDRO_37 *= DENDRO_43;

    // initialize reduction for
    DENDRO_43 = 0;

    DENDRO_43 += DENDRO_40;
    DENDRO_43 += DENDRO_38;
    DENDRO_43 += DENDRO_36;
    DENDRO_43 += DENDRO_33;
    DENDRO_43 += DENDRO_9;
    DENDRO_43 += DENDRO_50;
    DENDRO_9  = gt5[pp];
    DENDRO_33 = -1.0 / 3.0;

    // initialize reduction for
    DENDRO_36 = 0;

    DENDRO_36 += DENDRO_52;
    DENDRO_36 += DENDRO_51;
    DENDRO_36 += DENDRO_42;
    DENDRO_38 = -2;

    // initialize reduction for
    DENDRO_40 = 0;

    DENDRO_40 += DENDRO_55;
    DENDRO_40 += DENDRO_54;
    DENDRO_40 += DENDRO_53;

    // initialize reduction for
    DENDRO_42 = 0;

    DENDRO_42 += DENDRO_58;
    DENDRO_42 += DENDRO_57;
    DENDRO_42 += DENDRO_56;
    DENDRO_50 = K[pp];
    DENDRO_53 = gt4[pp];
    DENDRO_54 = gt3[pp];

    // initialize reduction for
    DENDRO_55 = 0;

    DENDRO_55 += DENDRO_61;
    DENDRO_55 += DENDRO_60;
    DENDRO_55 += DENDRO_59;

    // initialize reduction for
    DENDRO_56 = 0;

    DENDRO_56 += DENDRO_63;
    DENDRO_56 += DENDRO_62;
    DENDRO_56 += DENDRO_51;

    // initialize reduction for
    DENDRO_51 = 0;

    DENDRO_51 += DENDRO_66;
    DENDRO_51 += DENDRO_65;
    DENDRO_51 += DENDRO_64;
    DENDRO_57 = gt2[pp];
    DENDRO_58 = gt1[pp];
    DENDRO_59 = gt0[pp];

    // initialize reduction for
    DENDRO_60 = 0;

    DENDRO_60 += DENDRO_68;
    DENDRO_60 += DENDRO_39;
    DENDRO_60 += DENDRO_67;

    // initialize reduction for
    DENDRO_39 = 0;

    DENDRO_39 += DENDRO_41;
    DENDRO_39 += DENDRO_31;
    DENDRO_39 += DENDRO_34;

    // initialize reduction for
    DENDRO_31 = 0;

    DENDRO_31 += DENDRO_37;
    DENDRO_31 += DENDRO_63;
    DENDRO_31 += DENDRO_52;

    // initialize reduction for
    DENDRO_34 = 1;

    DENDRO_34 *= DENDRO_33;
    DENDRO_34 *= DENDRO_9;
    DENDRO_34 *= DENDRO_43;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_38;
    DENDRO_9 *= DENDRO_44;
    DENDRO_9 *= DENDRO_36;

    // initialize reduction for
    DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_38;
    DENDRO_37 *= DENDRO_45;
    DENDRO_37 *= DENDRO_40;

    // initialize reduction for
    DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_38;
    DENDRO_41 *= DENDRO_46;
    DENDRO_41 *= DENDRO_42;

    // initialize reduction for
    DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_44;
    DENDRO_52 *= DENDRO_50;

    // initialize reduction for
    DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_33;
    DENDRO_61 *= DENDRO_53;
    DENDRO_61 *= DENDRO_43;

    // initialize reduction for
    DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_38;
    DENDRO_53 *= DENDRO_45;
    DENDRO_53 *= DENDRO_36;

    // initialize reduction for
    DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_38;
    DENDRO_62 *= DENDRO_47;
    DENDRO_62 *= DENDRO_40;

    // initialize reduction for
    DENDRO_63 = 1;

    DENDRO_63 *= DENDRO_38;
    DENDRO_63 *= DENDRO_48;
    DENDRO_63 *= DENDRO_42;

    // initialize reduction for
    DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_45;
    DENDRO_64 *= DENDRO_50;

    // initialize reduction for
    DENDRO_65 = 1;

    DENDRO_65 *= DENDRO_33;
    DENDRO_65 *= DENDRO_54;
    DENDRO_65 *= DENDRO_43;

    // initialize reduction for
    DENDRO_54 = 1;

    DENDRO_54 *= DENDRO_38;
    DENDRO_54 *= DENDRO_45;
    DENDRO_54 *= DENDRO_55;

    // initialize reduction for
    DENDRO_66 = 1;

    DENDRO_66 *= DENDRO_38;
    DENDRO_66 *= DENDRO_47;
    DENDRO_66 *= DENDRO_56;

    // initialize reduction for
    DENDRO_67 = 1;

    DENDRO_67 *= DENDRO_38;
    DENDRO_67 *= DENDRO_48;
    DENDRO_67 *= DENDRO_51;

    // initialize reduction for
    DENDRO_68 = 1;

    DENDRO_68 *= DENDRO_47;
    DENDRO_68 *= DENDRO_50;

    // initialize reduction for
    double DENDRO_69 = 1;

    DENDRO_69 *= DENDRO_33;
    DENDRO_69 *= DENDRO_57;
    DENDRO_69 *= DENDRO_43;

    // initialize reduction for
    DENDRO_57 = 1;

    DENDRO_57 *= DENDRO_38;
    DENDRO_57 *= DENDRO_46;
    DENDRO_57 *= DENDRO_36;

    // initialize reduction for
    DENDRO_36 = 1;

    DENDRO_36 *= DENDRO_38;
    DENDRO_36 *= DENDRO_48;
    DENDRO_36 *= DENDRO_40;

    // initialize reduction for
    DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_38;
    DENDRO_40 *= DENDRO_49;
    DENDRO_40 *= DENDRO_42;

    // initialize reduction for
    DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_46;
    DENDRO_42 *= DENDRO_50;

    // initialize reduction for
    double DENDRO_70 = 1;

    DENDRO_70 *= DENDRO_33;
    DENDRO_70 *= DENDRO_58;
    DENDRO_70 *= DENDRO_43;

    // initialize reduction for
    DENDRO_58 = 1;

    DENDRO_58 *= DENDRO_38;
    DENDRO_58 *= DENDRO_46;
    DENDRO_58 *= DENDRO_55;

    // initialize reduction for
    DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_38;
    DENDRO_55 *= DENDRO_48;
    DENDRO_55 *= DENDRO_56;

    // initialize reduction for
    DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_38;
    DENDRO_56 *= DENDRO_49;
    DENDRO_56 *= DENDRO_51;

    // initialize reduction for
    DENDRO_51 = 1;

    DENDRO_51 *= DENDRO_48;
    DENDRO_51 *= DENDRO_50;

    // initialize reduction for
    double DENDRO_71 = 1;

    DENDRO_71 *= DENDRO_33;
    DENDRO_71 *= DENDRO_59;
    DENDRO_71 *= DENDRO_43;

    // initialize reduction for
    DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_38;
    DENDRO_33 *= DENDRO_46;
    DENDRO_33 *= DENDRO_60;

    // initialize reduction for
    DENDRO_43 = 1;

    DENDRO_43 *= DENDRO_38;
    DENDRO_43 *= DENDRO_48;
    DENDRO_43 *= DENDRO_39;

    // initialize reduction for
    DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_38;
    DENDRO_39 *= DENDRO_49;
    DENDRO_39 *= DENDRO_31;

    // initialize reduction for
    DENDRO_31 = 1;

    DENDRO_31 *= DENDRO_49;
    DENDRO_31 *= DENDRO_50;
    DENDRO_38        = grad_2_beta2;
    DENDRO_50        = 4.0 / 3.0;
    DENDRO_59        = grad_1_beta1;
    DENDRO_60        = -2.0 / 3.0;
    double DENDRO_72 = grad_0_beta0;
    double DENDRO_73 = grad_2_beta1;
    double DENDRO_74 = grad_2_beta0;

    // initialize reduction for
    double DENDRO_75 = 0;

    DENDRO_75 += DENDRO_23;
    DENDRO_75 += DENDRO_22;
    DENDRO_75 += DENDRO_21;
    DENDRO_75 += DENDRO_20;
    DENDRO_75 += DENDRO_19;
    DENDRO_75 += DENDRO_34;
    DENDRO_19        = chi[pp];
    DENDRO_20        = grad_2_At5;
    DENDRO_21        = beta2[pp];
    DENDRO_22        = grad_1_At5;
    DENDRO_23        = beta1[pp];
    DENDRO_34        = grad_0_At5;
    double DENDRO_76 = beta0[pp];

    // initialize reduction for
    double DENDRO_77 = 0;

    DENDRO_77 += DENDRO_52;
    DENDRO_77 += DENDRO_41;
    DENDRO_77 += DENDRO_37;
    DENDRO_77 += DENDRO_9;
    DENDRO_9  = 1.0 / 3.0;

    // initialize reduction for
    DENDRO_37 = 0;

    DENDRO_37 += DENDRO_7;
    DENDRO_37 += DENDRO_5;
    DENDRO_37 += DENDRO_3;
    DENDRO_37 += DENDRO_1;
    DENDRO_37 += DENDRO_35;
    DENDRO_37 += DENDRO_61;
    DENDRO_1 = grad_2_At4;
    DENDRO_3 = grad_1_At4;
    DENDRO_5 = grad_0_At4;

    // initialize reduction for
    DENDRO_7 = 0;

    DENDRO_7 += DENDRO_64;
    DENDRO_7 += DENDRO_63;
    DENDRO_7 += DENDRO_62;
    DENDRO_7 += DENDRO_53;
    DENDRO_35 = grad_1_beta2;
    DENDRO_41 = grad_1_beta0;

    // initialize reduction for
    DENDRO_52 = 0;

    DENDRO_52 += DENDRO_28;
    DENDRO_52 += DENDRO_27;
    DENDRO_52 += DENDRO_26;
    DENDRO_52 += DENDRO_25;
    DENDRO_52 += DENDRO_24;
    DENDRO_52 += DENDRO_65;
    DENDRO_24 = grad_2_At3;
    DENDRO_25 = grad_1_At3;
    DENDRO_26 = grad_0_At3;

    // initialize reduction for
    DENDRO_27 = 0;

    DENDRO_27 += DENDRO_68;
    DENDRO_27 += DENDRO_67;
    DENDRO_27 += DENDRO_66;
    DENDRO_27 += DENDRO_54;

    // initialize reduction for
    DENDRO_28 = 0;

    DENDRO_28 += DENDRO_13;
    DENDRO_28 += DENDRO_12;
    DENDRO_28 += DENDRO_11;
    DENDRO_28 += DENDRO_10;
    DENDRO_28 += DENDRO_8;
    DENDRO_28 += DENDRO_69;
    DENDRO_8  = grad_2_At2;
    DENDRO_10 = grad_1_At2;
    DENDRO_11 = grad_0_At2;

    // initialize reduction for
    DENDRO_12 = 0;

    DENDRO_12 += DENDRO_42;
    DENDRO_12 += DENDRO_40;
    DENDRO_12 += DENDRO_36;
    DENDRO_12 += DENDRO_57;
    DENDRO_13 = grad_0_beta2;
    DENDRO_36 = grad_0_beta1;

    // initialize reduction for
    DENDRO_40 = 0;

    DENDRO_40 += DENDRO_18;
    DENDRO_40 += DENDRO_17;
    DENDRO_40 += DENDRO_16;
    DENDRO_40 += DENDRO_15;
    DENDRO_40 += DENDRO_14;
    DENDRO_40 += DENDRO_70;
    DENDRO_14 = grad_2_At1;
    DENDRO_15 = grad_1_At1;
    DENDRO_16 = grad_0_At1;

    // initialize reduction for
    DENDRO_17 = 0;

    DENDRO_17 += DENDRO_51;
    DENDRO_17 += DENDRO_56;
    DENDRO_17 += DENDRO_55;
    DENDRO_17 += DENDRO_58;

    // initialize reduction for
    DENDRO_18 = 0;

    DENDRO_18 += DENDRO_6;
    DENDRO_18 += DENDRO_4;
    DENDRO_18 += DENDRO_2;
    DENDRO_18 += DENDRO_30;
    DENDRO_18 += DENDRO_29;
    DENDRO_18 += DENDRO_71;
    DENDRO_2  = grad_2_At0;
    DENDRO_4  = grad_1_At0;
    DENDRO_6  = grad_0_At0;

    // initialize reduction for
    DENDRO_29 = 0;

    DENDRO_29 += DENDRO_31;
    DENDRO_29 += DENDRO_39;
    DENDRO_29 += DENDRO_43;
    DENDRO_29 += DENDRO_33;

    // initialize reduction for
    DENDRO_30 = 1;

    DENDRO_30 *= DENDRO_50;
    DENDRO_30 *= DENDRO_44;
    DENDRO_30 *= DENDRO_38;

    // initialize reduction for
    DENDRO_31 = 1;

    DENDRO_31 *= DENDRO_60;
    DENDRO_31 *= DENDRO_44;
    DENDRO_31 *= DENDRO_59;

    // initialize reduction for
    DENDRO_33 = 1;

    DENDRO_33 *= DENDRO_60;
    DENDRO_33 *= DENDRO_44;
    DENDRO_33 *= DENDRO_72;

    // initialize reduction for
    DENDRO_39 = 1;

    DENDRO_39 *= DENDRO_32;
    DENDRO_39 *= DENDRO_45;
    DENDRO_39 *= DENDRO_73;

    // initialize reduction for
    DENDRO_42 = 1;

    DENDRO_42 *= DENDRO_32;
    DENDRO_42 *= DENDRO_46;
    DENDRO_42 *= DENDRO_74;

    // initialize reduction for
    DENDRO_43 = 1;

    DENDRO_43 *= DENDRO_19;
    DENDRO_43 *= DENDRO_75;

    // initialize reduction for
    DENDRO_51 = 1;

    DENDRO_51 *= DENDRO_21;
    DENDRO_51 *= DENDRO_20;

    // initialize reduction for
    DENDRO_20 = 1;

    DENDRO_20 *= DENDRO_23;
    DENDRO_20 *= DENDRO_22;

    // initialize reduction for
    DENDRO_22 = 1;

    DENDRO_22 *= DENDRO_76;
    DENDRO_22 *= DENDRO_34;

    // initialize reduction for
    DENDRO_34 = 1;

    DENDRO_34 *= DENDRO_0;
    DENDRO_34 *= DENDRO_77;

    // initialize reduction for
    DENDRO_53 = 1;

    DENDRO_53 *= DENDRO_9;
    DENDRO_53 *= DENDRO_45;
    DENDRO_53 *= DENDRO_38;

    // initialize reduction for
    DENDRO_54 = 1;

    DENDRO_54 *= DENDRO_9;
    DENDRO_54 *= DENDRO_45;
    DENDRO_54 *= DENDRO_59;

    // initialize reduction for
    DENDRO_55 = 1;

    DENDRO_55 *= DENDRO_60;
    DENDRO_55 *= DENDRO_45;
    DENDRO_55 *= DENDRO_72;

    // initialize reduction for
    DENDRO_56 = 1;

    DENDRO_56 *= DENDRO_19;
    DENDRO_56 *= DENDRO_37;

    // initialize reduction for
    DENDRO_37 = 1;

    DENDRO_37 *= DENDRO_21;
    DENDRO_37 *= DENDRO_1;

    // initialize reduction for
    DENDRO_1 = 1;

    DENDRO_1 *= DENDRO_23;
    DENDRO_1 *= DENDRO_3;

    // initialize reduction for
    DENDRO_3 = 1;

    DENDRO_3 *= DENDRO_76;
    DENDRO_3 *= DENDRO_5;

    // initialize reduction for
    DENDRO_5 = 1;

    DENDRO_5 *= DENDRO_0;
    DENDRO_5 *= DENDRO_7;

    // initialize reduction for
    DENDRO_7 = 1;

    DENDRO_7 *= DENDRO_44;
    DENDRO_7 *= DENDRO_35;

    // initialize reduction for
    DENDRO_57 = 1;

    DENDRO_57 *= DENDRO_47;
    DENDRO_57 *= DENDRO_73;

    // initialize reduction for
    DENDRO_58 = 1;

    DENDRO_58 *= DENDRO_46;
    DENDRO_58 *= DENDRO_41;

    // initialize reduction for
    DENDRO_61 = 1;

    DENDRO_61 *= DENDRO_48;
    DENDRO_61 *= DENDRO_74;

    // initialize reduction for
    DENDRO_62 = 1;

    DENDRO_62 *= DENDRO_50;
    DENDRO_62 *= DENDRO_47;
    DENDRO_62 *= DENDRO_59;

    // initialize reduction for
    DENDRO_63 = 1;

    DENDRO_63 *= DENDRO_60;
    DENDRO_63 *= DENDRO_47;
    DENDRO_63 *= DENDRO_38;

    // initialize reduction for
    DENDRO_64 = 1;

    DENDRO_64 *= DENDRO_60;
    DENDRO_64 *= DENDRO_47;
    DENDRO_64 *= DENDRO_72;

    // initialize reduction for
    DENDRO_65 = 1;

    DENDRO_65 *= DENDRO_32;
    DENDRO_65 *= DENDRO_45;
    DENDRO_65 *= DENDRO_35;

    // initialize reduction for
    DENDRO_66 = 1;

    DENDRO_66 *= DENDRO_32;
    DENDRO_66 *= DENDRO_48;
    DENDRO_66 *= DENDRO_41;

    // initialize reduction for
    DENDRO_67 = 1;

    DENDRO_67 *= DENDRO_19;
    DENDRO_67 *= DENDRO_52;

    // initialize reduction for
    DENDRO_52 = 1;

    DENDRO_52 *= DENDRO_21;
    DENDRO_52 *= DENDRO_24;

    // initialize reduction for
    DENDRO_24 = 1;

    DENDRO_24 *= DENDRO_23;
    DENDRO_24 *= DENDRO_25;

    // initialize reduction for
    DENDRO_25 = 1;

    DENDRO_25 *= DENDRO_76;
    DENDRO_25 *= DENDRO_26;

    // initialize reduction for
    DENDRO_26 = 1;

    DENDRO_26 *= DENDRO_0;
    DENDRO_26 *= DENDRO_27;

    // initialize reduction for
    DENDRO_27 = 1;

    DENDRO_27 *= DENDRO_9;
    DENDRO_27 *= DENDRO_46;
    DENDRO_27 *= DENDRO_38;

    // initialize reduction for
    DENDRO_68 = 1;

    DENDRO_68 *= DENDRO_9;
    DENDRO_68 *= DENDRO_46;
    DENDRO_68 *= DENDRO_72;

    // initialize reduction for
    DENDRO_69 = 1;

    DENDRO_69 *= DENDRO_60;
    DENDRO_69 *= DENDRO_46;
    DENDRO_69 *= DENDRO_59;

    // initialize reduction for
    DENDRO_70 = 1;

    DENDRO_70 *= DENDRO_19;
    DENDRO_70 *= DENDRO_28;

    // initialize reduction for
    DENDRO_28 = 1;

    DENDRO_28 *= DENDRO_21;
    DENDRO_28 *= DENDRO_8;

    // initialize reduction for
    DENDRO_8 = 1;

    DENDRO_8 *= DENDRO_23;
    DENDRO_8 *= DENDRO_10;

    // initialize reduction for
    DENDRO_10 = 1;

    DENDRO_10 *= DENDRO_76;
    DENDRO_10 *= DENDRO_11;

    // initialize reduction for
    DENDRO_11 = 1;

    DENDRO_11 *= DENDRO_0;
    DENDRO_11 *= DENDRO_12;

    // initialize reduction for
    DENDRO_12 = 1;

    DENDRO_12 *= DENDRO_44;
    DENDRO_12 *= DENDRO_13;

    // initialize reduction for
    DENDRO_44 = 1;

    DENDRO_44 *= DENDRO_45;
    DENDRO_44 *= DENDRO_36;

    // initialize reduction for
    DENDRO_71 = 1;

    DENDRO_71 *= DENDRO_48;
    DENDRO_71 *= DENDRO_73;

    // initialize reduction for
    DENDRO_73 = 1;

    DENDRO_73 *= DENDRO_49;
    DENDRO_73 *= DENDRO_74;

    // initialize reduction for
    DENDRO_74 = 1;

    DENDRO_74 *= DENDRO_9;
    DENDRO_74 *= DENDRO_48;
    DENDRO_74 *= DENDRO_59;

    // initialize reduction for
    DENDRO_75 = 1;

    DENDRO_75 *= DENDRO_9;
    DENDRO_75 *= DENDRO_48;
    DENDRO_75 *= DENDRO_72;

    // initialize reduction for
    DENDRO_9 = 1;

    DENDRO_9 *= DENDRO_60;
    DENDRO_9 *= DENDRO_48;
    DENDRO_9 *= DENDRO_38;

    // initialize reduction for
    DENDRO_77 = 1;

    DENDRO_77 *= DENDRO_19;
    DENDRO_77 *= DENDRO_40;

    // initialize reduction for
    DENDRO_40 = 1;

    DENDRO_40 *= DENDRO_21;
    DENDRO_40 *= DENDRO_14;

    // initialize reduction for
    DENDRO_14 = 1;

    DENDRO_14 *= DENDRO_23;
    DENDRO_14 *= DENDRO_15;

    // initialize reduction for
    DENDRO_15 = 1;

    DENDRO_15 *= DENDRO_76;
    DENDRO_15 *= DENDRO_16;

    // initialize reduction for
    DENDRO_16 = 1;

    DENDRO_16 *= DENDRO_0;
    DENDRO_16 *= DENDRO_17;

    // initialize reduction for
    DENDRO_17 = 1;

    DENDRO_17 *= DENDRO_45;
    DENDRO_17 *= DENDRO_13;

    // initialize reduction for
    DENDRO_45 = 1;

    DENDRO_45 *= DENDRO_47;
    DENDRO_45 *= DENDRO_36;

    // initialize reduction for
    DENDRO_47 = 1;

    DENDRO_47 *= DENDRO_46;
    DENDRO_47 *= DENDRO_35;

    // initialize reduction for
    DENDRO_35 = 1;

    DENDRO_35 *= DENDRO_49;
    DENDRO_35 *= DENDRO_41;

    // initialize reduction for
    DENDRO_41 = 1;

    DENDRO_41 *= DENDRO_50;
    DENDRO_41 *= DENDRO_49;
    DENDRO_41 *= DENDRO_72;

    // initialize reduction for
    DENDRO_50 = 1;

    DENDRO_50 *= DENDRO_60;
    DENDRO_50 *= DENDRO_49;
    DENDRO_50 *= DENDRO_38;

    // initialize reduction for
    DENDRO_38 = 1;

    DENDRO_38 *= DENDRO_60;
    DENDRO_38 *= DENDRO_49;
    DENDRO_38 *= DENDRO_59;

    // initialize reduction for
    DENDRO_49 = 1;

    DENDRO_49 *= DENDRO_32;
    DENDRO_49 *= DENDRO_46;
    DENDRO_49 *= DENDRO_13;

    // initialize reduction for
    DENDRO_13 = 1;

    DENDRO_13 *= DENDRO_32;
    DENDRO_13 *= DENDRO_48;
    DENDRO_13 *= DENDRO_36;

    // initialize reduction for
    DENDRO_32 = 1;

    DENDRO_32 *= DENDRO_19;
    DENDRO_32 *= DENDRO_18;

    // initialize reduction for
    DENDRO_18 = 1;

    DENDRO_18 *= DENDRO_21;
    DENDRO_18 *= DENDRO_2;

    // initialize reduction for
    DENDRO_2 = 1;

    DENDRO_2 *= DENDRO_23;
    DENDRO_2 *= DENDRO_4;

    // initialize reduction for
    DENDRO_4 = 1;

    DENDRO_4 *= DENDRO_76;
    DENDRO_4 *= DENDRO_6;

    // initialize reduction for
    DENDRO_6 = 1;

    DENDRO_6 *= DENDRO_0;
    DENDRO_6 *= DENDRO_29;

    // initialize reduction for
    DENDRO_0 = 0;

    DENDRO_0 += DENDRO_34;
    DENDRO_0 += DENDRO_22;
    DENDRO_0 += DENDRO_20;
    DENDRO_0 += DENDRO_51;
    DENDRO_0 += DENDRO_43;
    DENDRO_0 += DENDRO_42;
    DENDRO_0 += DENDRO_39;
    DENDRO_0 += DENDRO_33;
    DENDRO_0 += DENDRO_31;
    DENDRO_0 += DENDRO_30;
    At_rhs22[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0     = 0;

    DENDRO_0 += DENDRO_61;
    DENDRO_0 += DENDRO_58;
    DENDRO_0 += DENDRO_57;
    DENDRO_0 += DENDRO_7;
    DENDRO_0 += DENDRO_5;
    DENDRO_0 += DENDRO_3;
    DENDRO_0 += DENDRO_1;
    DENDRO_0 += DENDRO_37;
    DENDRO_0 += DENDRO_56;
    DENDRO_0 += DENDRO_55;
    DENDRO_0 += DENDRO_54;
    DENDRO_0 += DENDRO_53;
    At_rhs12[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0     = 0;

    DENDRO_0 += DENDRO_26;
    DENDRO_0 += DENDRO_25;
    DENDRO_0 += DENDRO_24;
    DENDRO_0 += DENDRO_52;
    DENDRO_0 += DENDRO_67;
    DENDRO_0 += DENDRO_66;
    DENDRO_0 += DENDRO_65;
    DENDRO_0 += DENDRO_64;
    DENDRO_0 += DENDRO_63;
    DENDRO_0 += DENDRO_62;
    At_rhs11[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0     = 0;

    DENDRO_0 += DENDRO_73;
    DENDRO_0 += DENDRO_71;
    DENDRO_0 += DENDRO_44;
    DENDRO_0 += DENDRO_12;
    DENDRO_0 += DENDRO_11;
    DENDRO_0 += DENDRO_10;
    DENDRO_0 += DENDRO_8;
    DENDRO_0 += DENDRO_28;
    DENDRO_0 += DENDRO_70;
    DENDRO_0 += DENDRO_69;
    DENDRO_0 += DENDRO_68;
    DENDRO_0 += DENDRO_27;
    At_rhs02[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0     = 0;

    DENDRO_0 += DENDRO_35;
    DENDRO_0 += DENDRO_47;
    DENDRO_0 += DENDRO_45;
    DENDRO_0 += DENDRO_17;
    DENDRO_0 += DENDRO_16;
    DENDRO_0 += DENDRO_15;
    DENDRO_0 += DENDRO_14;
    DENDRO_0 += DENDRO_40;
    DENDRO_0 += DENDRO_77;
    DENDRO_0 += DENDRO_9;
    DENDRO_0 += DENDRO_75;
    DENDRO_0 += DENDRO_74;
    At_rhs01[pp] = DENDRO_0;

    // initialize reduction for
    DENDRO_0     = 0;

    DENDRO_0 += DENDRO_6;
    DENDRO_0 += DENDRO_4;
    DENDRO_0 += DENDRO_2;
    DENDRO_0 += DENDRO_18;
    DENDRO_0 += DENDRO_32;
    DENDRO_0 += DENDRO_13;
    DENDRO_0 += DENDRO_49;
    DENDRO_0 += DENDRO_38;
    DENDRO_0 += DENDRO_50;
    DENDRO_0 += DENDRO_41;
    At_rhs00[pp] = DENDRO_0;
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

`