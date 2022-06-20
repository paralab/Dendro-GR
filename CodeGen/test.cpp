//Namespace(graph_name='DAG_BSSN.dat', iteration=1, mode=0, opt_ts_name='ts_opt.dat')
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , alpha, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_alpha=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_alpha=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_alpha=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_alpha=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_alpha=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_alpha=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_alpha=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_alpha=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_alpha=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_alpha=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_alpha=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_alpha=Du[gidx];
__syncthreads();

//|V|= 17 |E| =16
//topological generations = 4
// [0]
// traversal id=0 temp_var_count=4
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;
DENDRO_0 = beta2[pp];

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_0 = grad_2_alpha;
DENDRO_1 *= DENDRO_0;
DENDRO_0 = beta1[pp];

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_0 = grad_1_alpha;
DENDRO_2 *= DENDRO_0;
DENDRO_0 = beta0[pp];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = grad_0_alpha;
DENDRO_3 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_3;
DENDRO_0 += DENDRO_2;
DENDRO_0 += DENDRO_1;
DENDRO_1 = -2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_1;
DENDRO_1 = K[pp];
DENDRO_2 *= DENDRO_1;
DENDRO_1 = alpha[pp];
DENDRO_2 *= DENDRO_1;
DENDRO_1 = lambda[0];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_1;
DENDRO_3 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_3;
DENDRO_0 += DENDRO_2;
a_rhs[pp] = DENDRO_0;

if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&a_rhs [gidx] , alpha[gidx], grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk);
}

a_rhs[pp]    += ko_sigma * (kograd_0_alpha + kograd_1_alpha + kograd_2_alpha);
__syncthreads();


}
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , beta0, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_beta0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_beta0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_beta0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_beta0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_beta0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_beta0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_beta0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_beta0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_beta0=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_beta0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_beta0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_beta0=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , beta1, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_beta1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_beta1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_beta1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_beta1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_beta1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_beta1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_beta1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_beta1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_beta1=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_beta1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_beta1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_beta1=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , beta2, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_beta2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_beta2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_beta2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_beta2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_beta2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_beta2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_beta2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_beta2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_beta2=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_beta2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_beta2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_beta2=Du[gidx];
__syncthreads();

//|V|= 45 |E| =52
//topological generations = 4
// [0]
// traversal id=0 temp_var_count=9
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;
DENDRO_0 = beta2[pp];

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = grad_2_beta2;
DENDRO_1 *= DENDRO_2;
DENDRO_2 = beta1[pp];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = grad_1_beta2;
DENDRO_3 *= DENDRO_4;
DENDRO_4 = beta0[pp];

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = grad_0_beta2;
DENDRO_5 *= DENDRO_6;
DENDRO_6 = 3.0/4.0;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_6 = alpha[pp];
DENDRO_7 *= DENDRO_6;
DENDRO_6 = lambda_f[1];
DENDRO_7 *= DENDRO_6;
DENDRO_6 = 3.0/4.0;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_6;
DENDRO_6 = lambda_f[0];
DENDRO_8 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_5;
DENDRO_6 += DENDRO_3;
DENDRO_6 += DENDRO_1;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_8;
DENDRO_1 += DENDRO_7;
DENDRO_3 = lambda[1];

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_3;
DENDRO_5 *= DENDRO_6;
DENDRO_6 = B2[pp];

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_1;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_7;
DENDRO_6 += DENDRO_5;
b_rhs2[pp] = DENDRO_6;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_0;
DENDRO_6 = grad_2_beta1;
DENDRO_5 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_2;
DENDRO_7 = grad_1_beta1;
DENDRO_6 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_4;
DENDRO_8 = grad_0_beta1;
DENDRO_7 *= DENDRO_8;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_7;
DENDRO_8 += DENDRO_6;
DENDRO_8 += DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_3;
DENDRO_5 *= DENDRO_8;
DENDRO_6 = B1[pp];

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_1;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_7;
DENDRO_6 += DENDRO_5;
b_rhs1[pp] = DENDRO_6;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_0;
DENDRO_0 = grad_2_beta0;
DENDRO_5 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_2;
DENDRO_2 = grad_1_beta0;
DENDRO_0 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_4;
DENDRO_4 = grad_0_beta0;
DENDRO_2 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_2;
DENDRO_4 += DENDRO_0;
DENDRO_4 += DENDRO_5;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_3;
DENDRO_0 *= DENDRO_4;
DENDRO_2 = B0[pp];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_2;
DENDRO_3 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_3;
DENDRO_1 += DENDRO_0;
b_rhs0[pp] = DENDRO_1;


}
if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&b_rhs0[gidx] , beta0[gidx], grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&b_rhs1[gidx] , beta1[gidx], grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&b_rhs2[gidx] , beta2[gidx], grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk);
}

b_rhs0[pp]   += ko_sigma * (kograd_0_beta0 + kograd_1_beta0 + kograd_2_beta0);
b_rhs1[pp]   += ko_sigma * (kograd_0_beta1 + kograd_1_beta1 + kograd_2_beta1);
b_rhs2[pp]   += ko_sigma * (kograd_0_beta2 + kograd_1_beta2 + kograd_2_beta2);
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , chi, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_chi=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_chi=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_chi=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_chi=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_chi=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_chi=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_chi=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_chi=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_chi=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_chi=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_chi=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_chi=Du[gidx];
__syncthreads();

//|V|= 21 |E| =21
//topological generations = 3
// [0]
// traversal id=0 temp_var_count=6
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;
DENDRO_0 = grad_0_beta0;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_0;
DENDRO_0 = grad_1_beta1;
DENDRO_1 += DENDRO_0;
DENDRO_0 = grad_2_beta2;
DENDRO_1 += DENDRO_0;
DENDRO_0 = 2.0/3.0;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_0 = K[pp];
DENDRO_2 *= DENDRO_0;
DENDRO_0 = alpha[pp];
DENDRO_2 *= DENDRO_0;
DENDRO_0 = chi[pp];
DENDRO_2 *= DENDRO_0;
DENDRO_3 = -2.0/3.0;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_3;
DENDRO_4 *= DENDRO_0;
DENDRO_4 *= DENDRO_1;
DENDRO_0 = beta2[pp];

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_0 = grad_2_chi;
DENDRO_1 *= DENDRO_0;
DENDRO_0 = beta1[pp];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = grad_1_chi;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = beta0[pp];

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_0;
DENDRO_0 = grad_0_chi;
DENDRO_5 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_5;
DENDRO_0 += DENDRO_3;
DENDRO_0 += DENDRO_1;
DENDRO_0 += DENDRO_4;
DENDRO_0 += DENDRO_2;
chi_rhs[pp] = DENDRO_0;


}
if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&chi_rhs[gidx], chi[gidx]  , grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk);
}

chi_rhs[pp]  += ko_sigma * (kograd_0_chi + kograd_1_chi + kograd_2_chi);
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , gt0, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt0=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt0=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , gt1, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt1=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt1=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , gt2, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt2=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt2=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , gt3, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt3=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt3=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt3=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt3=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt3=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt3=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt3=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , gt4, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt4=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt4=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt4=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt4=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt4=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt4=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt4=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , gt5, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_gt5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_gt5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_gt5=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_gt5=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_gt5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_gt5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_gt5=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_gt5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_gt5=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_gt5=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_gt5=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_gt5=Du[gidx];
__syncthreads();

//|V|= 114 |E| =210
//topological generations = 2
// [0]
// traversal id=0 temp_var_count=33
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;
DENDRO_0 = 4.0/3.0;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = grad_2_beta2;
DENDRO_1 *= DENDRO_2;
DENDRO_3 = gt5[pp];
DENDRO_1 *= DENDRO_3;
DENDRO_4 = -2.0/3.0;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = grad_1_beta1;
DENDRO_5 *= DENDRO_6;
DENDRO_5 *= DENDRO_3;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_4;
DENDRO_8 = grad_0_beta0;
DENDRO_7 *= DENDRO_8;
DENDRO_7 *= DENDRO_3;
DENDRO_9 = 2;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_9;
DENDRO_11 = grad_2_beta1;
DENDRO_10 *= DENDRO_11;
DENDRO_12 = gt4[pp];
DENDRO_10 *= DENDRO_12;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_9;
DENDRO_14 = grad_2_beta0;
DENDRO_13 *= DENDRO_14;
DENDRO_15 = gt2[pp];
DENDRO_13 *= DENDRO_15;
DENDRO_16 = -2;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_16;
DENDRO_18 = At5[pp];
DENDRO_17 *= DENDRO_18;
DENDRO_18 = alpha[pp];
DENDRO_17 *= DENDRO_18;
DENDRO_19 = beta2[pp];

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_19;
DENDRO_21 = grad_2_gt5;
DENDRO_20 *= DENDRO_21;
DENDRO_21 = beta1[pp];

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_21;
DENDRO_23 = grad_1_gt5;
DENDRO_22 *= DENDRO_23;
DENDRO_23 = beta0[pp];

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_23;
DENDRO_25 = grad_0_gt5;
DENDRO_24 *= DENDRO_25;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_24;
DENDRO_25 += DENDRO_22;
DENDRO_25 += DENDRO_20;
DENDRO_25 += DENDRO_17;
DENDRO_25 += DENDRO_13;
DENDRO_25 += DENDRO_10;
DENDRO_25 += DENDRO_7;
DENDRO_25 += DENDRO_5;
DENDRO_25 += DENDRO_1;
gt_rhs22[pp] = DENDRO_25;
DENDRO_1 = 1.0/3.0;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_1;
DENDRO_5 *= DENDRO_2;
DENDRO_5 *= DENDRO_12;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_12;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_4;
DENDRO_10 *= DENDRO_8;
DENDRO_10 *= DENDRO_12;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_16;
DENDRO_17 = At4[pp];
DENDRO_13 *= DENDRO_17;
DENDRO_13 *= DENDRO_18;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_11;
DENDRO_20 = gt3[pp];
DENDRO_17 *= DENDRO_20;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_14;
DENDRO_24 = gt1[pp];
DENDRO_22 *= DENDRO_24;
DENDRO_25 = grad_1_beta2;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_25;
DENDRO_26 *= DENDRO_3;
DENDRO_27 = grad_1_beta0;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_27;
DENDRO_28 *= DENDRO_15;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_19;
DENDRO_30 = grad_2_gt4;
DENDRO_29 *= DENDRO_30;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_21;
DENDRO_31 = grad_1_gt4;
DENDRO_30 *= DENDRO_31;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_23;
DENDRO_32 = grad_0_gt4;
DENDRO_31 *= DENDRO_32;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_31;
DENDRO_32 += DENDRO_30;
DENDRO_32 += DENDRO_29;
DENDRO_32 += DENDRO_28;
DENDRO_32 += DENDRO_26;
DENDRO_32 += DENDRO_22;
DENDRO_32 += DENDRO_17;
DENDRO_32 += DENDRO_13;
DENDRO_32 += DENDRO_10;
DENDRO_32 += DENDRO_7;
DENDRO_32 += DENDRO_5;
gt_rhs12[pp] = DENDRO_32;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_0;
DENDRO_5 *= DENDRO_6;
DENDRO_5 *= DENDRO_20;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_4;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_20;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_4;
DENDRO_10 *= DENDRO_8;
DENDRO_10 *= DENDRO_20;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_9;
DENDRO_13 *= DENDRO_25;
DENDRO_13 *= DENDRO_12;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_9;
DENDRO_17 *= DENDRO_27;
DENDRO_17 *= DENDRO_24;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_16;
DENDRO_26 = At3[pp];
DENDRO_22 *= DENDRO_26;
DENDRO_22 *= DENDRO_18;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_19;
DENDRO_28 = grad_2_gt3;
DENDRO_26 *= DENDRO_28;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_21;
DENDRO_29 = grad_1_gt3;
DENDRO_28 *= DENDRO_29;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_23;
DENDRO_30 = grad_0_gt3;
DENDRO_29 *= DENDRO_30;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_29;
DENDRO_30 += DENDRO_28;
DENDRO_30 += DENDRO_26;
DENDRO_30 += DENDRO_22;
DENDRO_30 += DENDRO_17;
DENDRO_30 += DENDRO_13;
DENDRO_30 += DENDRO_10;
DENDRO_30 += DENDRO_7;
DENDRO_30 += DENDRO_5;
gt_rhs11[pp] = DENDRO_30;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_1;
DENDRO_5 *= DENDRO_2;
DENDRO_5 *= DENDRO_15;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_1;
DENDRO_7 *= DENDRO_8;
DENDRO_7 *= DENDRO_15;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_4;
DENDRO_10 *= DENDRO_6;
DENDRO_10 *= DENDRO_15;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_16;
DENDRO_17 = At2[pp];
DENDRO_13 *= DENDRO_17;
DENDRO_13 *= DENDRO_18;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_11;
DENDRO_17 *= DENDRO_24;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_14;
DENDRO_14 = gt0[pp];
DENDRO_11 *= DENDRO_14;
DENDRO_22 = grad_0_beta2;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_22;
DENDRO_26 *= DENDRO_3;
DENDRO_3 = grad_0_beta1;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_3;
DENDRO_28 *= DENDRO_12;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_19;
DENDRO_30 = grad_2_gt2;
DENDRO_29 *= DENDRO_30;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_21;
DENDRO_31 = grad_1_gt2;
DENDRO_30 *= DENDRO_31;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_23;
DENDRO_32 = grad_0_gt2;
DENDRO_31 *= DENDRO_32;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_31;
DENDRO_32 += DENDRO_30;
DENDRO_32 += DENDRO_29;
DENDRO_32 += DENDRO_28;
DENDRO_32 += DENDRO_26;
DENDRO_32 += DENDRO_11;
DENDRO_32 += DENDRO_17;
DENDRO_32 += DENDRO_13;
DENDRO_32 += DENDRO_10;
DENDRO_32 += DENDRO_7;
DENDRO_32 += DENDRO_5;
gt_rhs02[pp] = DENDRO_32;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_1;
DENDRO_5 *= DENDRO_6;
DENDRO_5 *= DENDRO_24;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_1;
DENDRO_7 *= DENDRO_8;
DENDRO_7 *= DENDRO_24;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_4;
DENDRO_1 *= DENDRO_2;
DENDRO_1 *= DENDRO_24;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_16;
DENDRO_11 = At1[pp];
DENDRO_10 *= DENDRO_11;
DENDRO_10 *= DENDRO_18;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_25;
DENDRO_11 *= DENDRO_15;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_27;
DENDRO_13 *= DENDRO_14;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_22;
DENDRO_17 *= DENDRO_12;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_3;
DENDRO_12 *= DENDRO_20;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_19;
DENDRO_25 = grad_2_gt1;
DENDRO_20 *= DENDRO_25;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_21;
DENDRO_26 = grad_1_gt1;
DENDRO_25 *= DENDRO_26;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_23;
DENDRO_27 = grad_0_gt1;
DENDRO_26 *= DENDRO_27;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_26;
DENDRO_27 += DENDRO_25;
DENDRO_27 += DENDRO_20;
DENDRO_27 += DENDRO_12;
DENDRO_27 += DENDRO_17;
DENDRO_27 += DENDRO_13;
DENDRO_27 += DENDRO_11;
DENDRO_27 += DENDRO_10;
DENDRO_27 += DENDRO_1;
DENDRO_27 += DENDRO_7;
DENDRO_27 += DENDRO_5;
gt_rhs01[pp] = DENDRO_27;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_1 *= DENDRO_8;
DENDRO_1 *= DENDRO_14;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_4;
DENDRO_0 *= DENDRO_2;
DENDRO_0 *= DENDRO_14;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_4;
DENDRO_2 *= DENDRO_6;
DENDRO_2 *= DENDRO_14;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_9;
DENDRO_4 *= DENDRO_22;
DENDRO_4 *= DENDRO_15;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_9;
DENDRO_5 *= DENDRO_3;
DENDRO_5 *= DENDRO_24;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_16;
DENDRO_6 = At0[pp];
DENDRO_3 *= DENDRO_6;
DENDRO_3 *= DENDRO_18;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_19;
DENDRO_7 = grad_2_gt0;
DENDRO_6 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_21;
DENDRO_8 = grad_1_gt0;
DENDRO_7 *= DENDRO_8;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_23;
DENDRO_9 = grad_0_gt0;
DENDRO_8 *= DENDRO_9;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_8;
DENDRO_9 += DENDRO_7;
DENDRO_9 += DENDRO_6;
DENDRO_9 += DENDRO_3;
DENDRO_9 += DENDRO_5;
DENDRO_9 += DENDRO_4;
DENDRO_9 += DENDRO_2;
DENDRO_9 += DENDRO_0;
DENDRO_9 += DENDRO_1;
gt_rhs00[pp] = DENDRO_9;


}
if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&gt_rhs00[gidx], gt0[gidx], grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk);
radiative_bc_pt<pw , nx>(&gt_rhs01[gidx], gt1[gidx], grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&gt_rhs02[gidx], gt2[gidx], grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&gt_rhs11[gidx], gt3[gidx], grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk);
radiative_bc_pt<pw , nx>(&gt_rhs12[gidx], gt4[gidx], grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&gt_rhs22[gidx], gt5[gidx], grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk);
}

gt_rhs00[pp] += ko_sigma * (kograd_0_gt0 + kograd_1_gt0 + kograd_2_gt0);
gt_rhs01[pp] += ko_sigma * (kograd_0_gt1 + kograd_1_gt1 + kograd_2_gt1);
gt_rhs02[pp] += ko_sigma * (kograd_0_gt2 + kograd_1_gt2 + kograd_2_gt2);
gt_rhs11[pp] += ko_sigma * (kograd_0_gt3 + kograd_1_gt3 + kograd_2_gt3);
gt_rhs12[pp] += ko_sigma * (kograd_0_gt4 + kograd_1_gt4 + kograd_2_gt4);
gt_rhs22[pp] += ko_sigma * (kograd_0_gt5 + kograd_1_gt5 + kograd_2_gt5);
__syncthreads();

double DENDRO_igt0;
double DENDRO_igt1;
double DENDRO_igt2;
double DENDRO_igt3;
double DENDRO_igt4;
double DENDRO_igt5;
//|V|= 40 |E| =75
//topological generations = 5
// [0]
// traversal id=0 temp_var_count=17
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;
DENDRO_0 = gt1[pp];

// int-power reduce pow(gt1[pp],2)
DENDRO_1 = DENDRO_0;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = gt2[pp];

// int-power reduce pow(gt2[pp],2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = gt4[pp];

// int-power reduce pow(gt4[pp],2)
DENDRO_5 = DENDRO_4;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = -2;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_4;
DENDRO_6 = -1;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_6;
DENDRO_9 = gt0[pp];
DENDRO_8 *= DENDRO_9;
DENDRO_10 = gt3[pp];
DENDRO_8 *= DENDRO_10;
DENDRO_11 = gt5[pp];
DENDRO_8 *= DENDRO_11;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_11;
DENDRO_12 *= DENDRO_1;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_10;
DENDRO_13 *= DENDRO_3;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_9;
DENDRO_14 *= DENDRO_5;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_6;
DENDRO_15 *= DENDRO_9;
DENDRO_15 *= DENDRO_10;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_14;
DENDRO_16 += DENDRO_13;
DENDRO_16 += DENDRO_12;
DENDRO_16 += DENDRO_8;
DENDRO_16 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_1;
DENDRO_7 += DENDRO_15;

// int-power reduce pow(-gt0[pp]*gt3[pp]*gt5[pp] + gt0[pp]*gt4[pp]**2 + gt1[pp]**2*gt5[pp] - 2*gt1[pp]*gt2[pp]*gt4[pp] + gt2[pp]**2*gt3[pp],-1)
DENDRO_1 = DENDRO_16;
DENDRO_1 = 1.0/DENDRO_1;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_1;
DENDRO_8 *= DENDRO_7;
double DENDRO_igt22 = DENDRO_8;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_9;
DENDRO_8 *= DENDRO_4;

// initialize reduction 
DENDRO_12 = 0;
DENDRO_12 += DENDRO_8;
DENDRO_12 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_1;
DENDRO_7 *= DENDRO_12;
double DENDRO_igt12 = DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_9;
DENDRO_7 *= DENDRO_11;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_3;
DENDRO_8 += DENDRO_7;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_1;
DENDRO_3 *= DENDRO_8;
double DENDRO_igt11 = DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_6;
DENDRO_3 *= DENDRO_0;
DENDRO_3 *= DENDRO_4;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_10;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_7;
DENDRO_8 += DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_1;
DENDRO_3 *= DENDRO_8;
double DENDRO_igt02 = DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_6;
DENDRO_3 *= DENDRO_2;
DENDRO_3 *= DENDRO_4;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_2 *= DENDRO_11;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_2;
DENDRO_0 += DENDRO_3;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_1;
DENDRO_2 *= DENDRO_0;
double DENDRO_igt01 = DENDRO_2;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_6;
DENDRO_0 *= DENDRO_10;
DENDRO_0 *= DENDRO_11;

// initialize reduction 
DENDRO_2 = 0;
DENDRO_2 += DENDRO_5;
DENDRO_2 += DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_1;
DENDRO_0 *= DENDRO_2;
double DENDRO_igt00 = DENDRO_0;

DENDRO_igt0=DENDRO_igt00;
DENDRO_igt1=DENDRO_igt01;
DENDRO_igt2=DENDRO_igt02;
DENDRO_igt3=DENDRO_igt11;
DENDRO_igt4=DENDRO_igt12;
DENDRO_igt5=DENDRO_igt22;

}
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , K, blk);
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_K=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_K=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_K=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_K=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_K=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_K=Du[gidx];
__syncthreads();

//|V|= 426 |E| =936
//topological generations = 14
// [0]
// traversal id=0 temp_var_count=139
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;double DENDRO_33;double DENDRO_34;double DENDRO_35;double DENDRO_36;double DENDRO_37;double DENDRO_38;double DENDRO_39;double DENDRO_40;double DENDRO_41;double DENDRO_42;double DENDRO_43;double DENDRO_44;double DENDRO_45;double DENDRO_46;double DENDRO_47;double DENDRO_48;double DENDRO_49;double DENDRO_50;double DENDRO_51;double DENDRO_52;double DENDRO_53;double DENDRO_54;double DENDRO_55;double DENDRO_56;double DENDRO_57;double DENDRO_58;double DENDRO_59;double DENDRO_60;double DENDRO_61;double DENDRO_62;double DENDRO_63;double DENDRO_64;double DENDRO_65;double DENDRO_66;double DENDRO_67;double DENDRO_68;double DENDRO_69;double DENDRO_70;double DENDRO_71;double DENDRO_72;double DENDRO_73;double DENDRO_74;double DENDRO_75;double DENDRO_76;double DENDRO_77;double DENDRO_78;double DENDRO_79;double DENDRO_80;double DENDRO_81;double DENDRO_82;double DENDRO_83;double DENDRO_84;double DENDRO_85;double DENDRO_86;double DENDRO_87;double DENDRO_88;double DENDRO_89;double DENDRO_90;double DENDRO_91;double DENDRO_92;double DENDRO_93;double DENDRO_94;double DENDRO_95;double DENDRO_96;double DENDRO_97;double DENDRO_98;double DENDRO_99;double DENDRO_100;double DENDRO_101;double DENDRO_102;double DENDRO_103;double DENDRO_104;double DENDRO_105;double DENDRO_106;double DENDRO_107;double DENDRO_108;double DENDRO_109;double DENDRO_110;double DENDRO_111;double DENDRO_112;double DENDRO_113;double DENDRO_114;double DENDRO_115;double DENDRO_116;double DENDRO_117;double DENDRO_118;double DENDRO_119;double DENDRO_120;double DENDRO_121;double DENDRO_122;double DENDRO_123;double DENDRO_124;double DENDRO_125;double DENDRO_126;double DENDRO_127;double DENDRO_128;double DENDRO_129;double DENDRO_130;double DENDRO_131;double DENDRO_132;double DENDRO_133;double DENDRO_134;double DENDRO_135;double DENDRO_136;double DENDRO_137;double DENDRO_138;
DENDRO_0 = gt1[pp];

// int-power reduce pow(gt1[pp],2)
DENDRO_1 = DENDRO_0;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = gt2[pp];

// int-power reduce pow(gt2[pp],2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = gt4[pp];

// int-power reduce pow(gt4[pp],2)
DENDRO_5 = DENDRO_4;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = -2;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_4;
DENDRO_8 = -1;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_8;
DENDRO_10 = gt0[pp];
DENDRO_9 *= DENDRO_10;
DENDRO_11 = gt3[pp];
DENDRO_9 *= DENDRO_11;
DENDRO_12 = gt5[pp];
DENDRO_9 *= DENDRO_12;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_12;
DENDRO_13 *= DENDRO_1;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_11;
DENDRO_14 *= DENDRO_3;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_10;
DENDRO_15 *= DENDRO_5;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_8;
DENDRO_16 *= DENDRO_10;
DENDRO_16 *= DENDRO_11;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_15;
DENDRO_17 += DENDRO_14;
DENDRO_17 += DENDRO_13;
DENDRO_17 += DENDRO_9;
DENDRO_17 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_8;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_10;
DENDRO_9 *= DENDRO_4;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_8;
DENDRO_13 *= DENDRO_0;
DENDRO_13 *= DENDRO_4;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_2;
DENDRO_14 *= DENDRO_11;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_8;
DENDRO_15 *= DENDRO_2;
DENDRO_15 *= DENDRO_4;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_0;
DENDRO_18 *= DENDRO_12;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_8;
DENDRO_19 *= DENDRO_11;
DENDRO_19 *= DENDRO_12;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_8;
DENDRO_20 *= DENDRO_10;
DENDRO_20 *= DENDRO_12;

// initialize reduction 
DENDRO_21 = 0;
DENDRO_21 += DENDRO_1;
DENDRO_21 += DENDRO_16;

// int-power reduce pow(-gt0[pp]*gt3[pp]*gt5[pp] + gt0[pp]*gt4[pp]**2 + gt1[pp]**2*gt5[pp] - 2*gt1[pp]*gt2[pp]*gt4[pp] + gt2[pp]**2*gt3[pp],-1)
DENDRO_16 = DENDRO_17;
DENDRO_16 = 1.0/DENDRO_16;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_9;
DENDRO_17 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_14;
DENDRO_7 += DENDRO_13;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_18;
DENDRO_9 += DENDRO_15;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_5;
DENDRO_13 += DENDRO_19;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_3;
DENDRO_14 += DENDRO_20;
DENDRO_15 = grad_2_chi;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_15;
DENDRO_18 *= DENDRO_16;
DENDRO_18 *= DENDRO_21;
DENDRO_19 = grad_1_chi;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_19;
DENDRO_20 *= DENDRO_16;
DENDRO_20 *= DENDRO_17;
DENDRO_22 = grad_0_chi;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_22;
DENDRO_23 *= DENDRO_16;
DENDRO_23 *= DENDRO_7;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_15;
DENDRO_24 *= DENDRO_16;
DENDRO_24 *= DENDRO_7;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_19;
DENDRO_25 *= DENDRO_16;
DENDRO_25 *= DENDRO_9;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_22;
DENDRO_26 *= DENDRO_16;
DENDRO_26 *= DENDRO_13;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_15;
DENDRO_27 *= DENDRO_16;
DENDRO_27 *= DENDRO_17;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_19;
DENDRO_28 *= DENDRO_16;
DENDRO_28 *= DENDRO_14;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_22;
DENDRO_29 *= DENDRO_16;
DENDRO_29 *= DENDRO_9;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_23;
DENDRO_30 += DENDRO_20;
DENDRO_30 += DENDRO_18;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_26;
DENDRO_18 += DENDRO_25;
DENDRO_18 += DENDRO_24;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_29;
DENDRO_20 += DENDRO_28;
DENDRO_20 += DENDRO_27;
DENDRO_23 = -0.5;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_23;
DENDRO_23 = grad_1_gt2;
DENDRO_24 *= DENDRO_23;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_25;
DENDRO_25 = grad_2_gt1;
DENDRO_26 *= DENDRO_25;
DENDRO_27 = 0.5;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_27;
DENDRO_27 = grad_0_gt4;
DENDRO_28 *= DENDRO_27;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_8;
DENDRO_29 *= DENDRO_2;
DENDRO_29 *= DENDRO_30;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_8;
DENDRO_31 *= DENDRO_2;
DENDRO_31 *= DENDRO_18;
DENDRO_32 = -0.5;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_32;
DENDRO_33 *= DENDRO_25;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_25;
DENDRO_32 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_8;
DENDRO_23 *= DENDRO_0;
DENDRO_23 *= DENDRO_20;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_8;
DENDRO_25 *= DENDRO_0;
DENDRO_25 *= DENDRO_18;
DENDRO_34 = -0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_34;
DENDRO_35 *= DENDRO_27;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_8;
DENDRO_27 *= DENDRO_4;
DENDRO_27 *= DENDRO_30;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_8;
DENDRO_34 *= DENDRO_4;
DENDRO_34 *= DENDRO_20;
DENDRO_36 = -0.5;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = grad_1_gt0;
DENDRO_37 *= DENDRO_36;
DENDRO_38 = 1.0;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_38;
DENDRO_38 = grad_0_gt1;
DENDRO_39 *= DENDRO_38;
DENDRO_38 = -0.5;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_38;
DENDRO_38 = grad_2_gt0;
DENDRO_40 *= DENDRO_38;
DENDRO_41 = 1.0;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_41;
DENDRO_41 = grad_0_gt2;
DENDRO_42 *= DENDRO_41;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_8;
DENDRO_41 *= DENDRO_10;
DENDRO_41 *= DENDRO_18;
DENDRO_43 = 2;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_43;
DENDRO_44 *= DENDRO_22;
DENDRO_45 = -0.5;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_45;
DENDRO_45 = grad_0_gt3;
DENDRO_46 *= DENDRO_45;
DENDRO_47 = 1.0;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_47;
DENDRO_47 = grad_1_gt1;
DENDRO_48 *= DENDRO_47;
DENDRO_47 = -0.5;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_47;
DENDRO_47 = grad_2_gt3;
DENDRO_49 *= DENDRO_47;
DENDRO_50 = 1.0;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_50;
DENDRO_50 = grad_1_gt4;
DENDRO_51 *= DENDRO_50;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_8;
DENDRO_50 *= DENDRO_11;
DENDRO_50 *= DENDRO_20;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_43;
DENDRO_52 *= DENDRO_19;
DENDRO_53 = -0.5;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_53;
DENDRO_53 = grad_1_gt5;
DENDRO_54 *= DENDRO_53;
DENDRO_55 = 1.0;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_55;
DENDRO_55 = grad_2_gt4;
DENDRO_56 *= DENDRO_55;
DENDRO_55 = -0.5;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_55;
DENDRO_55 = grad_0_gt5;
DENDRO_57 *= DENDRO_55;
DENDRO_58 = 1.0;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_58;
DENDRO_58 = grad_2_gt2;
DENDRO_59 *= DENDRO_58;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_8;
DENDRO_58 *= DENDRO_12;
DENDRO_58 *= DENDRO_30;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_43;
DENDRO_60 *= DENDRO_15;

// initialize reduction 
DENDRO_61 = 0;
DENDRO_61 += DENDRO_28;
DENDRO_61 += DENDRO_26;
DENDRO_61 += DENDRO_24;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_22;
DENDRO_24 += DENDRO_29;
DENDRO_29 = chi[pp];

// int-power reduce pow(chi[pp],-1)
DENDRO_62 = DENDRO_29;
DENDRO_62 = 1.0/DENDRO_62;

// initialize reduction 
DENDRO_63 = 0;
DENDRO_63 += DENDRO_15;
DENDRO_63 += DENDRO_31;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_28;
DENDRO_31 += DENDRO_32;
DENDRO_31 += DENDRO_33;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_22;
DENDRO_28 += DENDRO_23;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_19;
DENDRO_22 += DENDRO_25;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_32;
DENDRO_23 += DENDRO_26;
DENDRO_23 += DENDRO_35;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_19;
DENDRO_25 += DENDRO_27;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_15;
DENDRO_19 += DENDRO_34;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_39;
DENDRO_15 += DENDRO_37;

// initialize reduction 
DENDRO_26 = 0;
DENDRO_26 += DENDRO_42;
DENDRO_26 += DENDRO_40;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_44;
DENDRO_27 += DENDRO_41;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_48;
DENDRO_32 += DENDRO_46;

// initialize reduction 
DENDRO_33 = 0;
DENDRO_33 += DENDRO_51;
DENDRO_33 += DENDRO_49;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_52;
DENDRO_34 += DENDRO_50;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_56;
DENDRO_35 += DENDRO_54;

// initialize reduction 
DENDRO_37 = 0;
DENDRO_37 += DENDRO_59;
DENDRO_37 += DENDRO_57;

// initialize reduction 
DENDRO_39 = 0;
DENDRO_39 += DENDRO_60;
DENDRO_39 += DENDRO_58;
DENDRO_40 = DENDRO_igt4;

// int-power reduce pow(DENDRO_igt4,2)
DENDRO_41 = DENDRO_40;
DENDRO_41 *= DENDRO_40;
DENDRO_42 = DENDRO_igt2;

// int-power reduce pow(DENDRO_igt2,2)
DENDRO_44 = DENDRO_42;
DENDRO_44 *= DENDRO_42;
DENDRO_46 = DENDRO_igt1;

// int-power reduce pow(DENDRO_igt1,2)
DENDRO_48 = DENDRO_46;
DENDRO_48 *= DENDRO_46;
DENDRO_49 = DENDRO_igt5;

// int-power reduce pow(DENDRO_igt5,2)
DENDRO_50 = DENDRO_49;
DENDRO_50 *= DENDRO_49;
DENDRO_51 = DENDRO_igt3;

// int-power reduce pow(DENDRO_igt3,2)
DENDRO_52 = DENDRO_51;
DENDRO_52 *= DENDRO_51;
DENDRO_54 = DENDRO_igt0;

// int-power reduce pow(DENDRO_igt0,2)
DENDRO_56 = DENDRO_54;
DENDRO_56 *= DENDRO_54;
DENDRO_57 = 0.5;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_57;
DENDRO_58 *= DENDRO_38;
DENDRO_58 *= DENDRO_16;
DENDRO_58 *= DENDRO_7;
DENDRO_57 = 0.5;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_57;
DENDRO_59 *= DENDRO_55;
DENDRO_59 *= DENDRO_16;
DENDRO_59 *= DENDRO_21;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_16;
DENDRO_57 *= DENDRO_17;
DENDRO_57 *= DENDRO_61;
DENDRO_60 = -0.5;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_60;
DENDRO_64 *= DENDRO_62;
DENDRO_64 *= DENDRO_24;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_24;
DENDRO_60 *= DENDRO_2;
DENDRO_60 *= DENDRO_62;
DENDRO_60 *= DENDRO_20;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_24;
DENDRO_65 *= DENDRO_38;
DENDRO_65 *= DENDRO_16;
DENDRO_65 *= DENDRO_9;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_24;
DENDRO_66 *= DENDRO_55;
DENDRO_66 *= DENDRO_16;
DENDRO_66 *= DENDRO_17;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_16;
DENDRO_24 *= DENDRO_14;
DENDRO_24 *= DENDRO_61;
DENDRO_67 = 0.5;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_67;
DENDRO_68 *= DENDRO_38;
DENDRO_68 *= DENDRO_16;
DENDRO_68 *= DENDRO_13;
DENDRO_38 = 0.5;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_38;
DENDRO_67 *= DENDRO_55;
DENDRO_67 *= DENDRO_16;
DENDRO_67 *= DENDRO_7;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_16;
DENDRO_38 *= DENDRO_9;
DENDRO_38 *= DENDRO_61;
DENDRO_55 = -0.5;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_55;
DENDRO_61 *= DENDRO_62;
DENDRO_61 *= DENDRO_63;
DENDRO_55 = 0.5;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_55;
DENDRO_63 *= DENDRO_0;
DENDRO_63 *= DENDRO_62;
DENDRO_63 *= DENDRO_30;
DENDRO_55 = 0.5;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_55;
DENDRO_69 *= DENDRO_36;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_7;
DENDRO_55 = 0.5;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_55;
DENDRO_70 *= DENDRO_45;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_17;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_16;
DENDRO_55 *= DENDRO_21;
DENDRO_55 *= DENDRO_31;
DENDRO_71 = 0.5;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_71;
DENDRO_72 *= DENDRO_36;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_9;
DENDRO_71 = 0.5;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_71;
DENDRO_73 *= DENDRO_45;
DENDRO_73 *= DENDRO_16;
DENDRO_73 *= DENDRO_14;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_16;
DENDRO_71 *= DENDRO_17;
DENDRO_71 *= DENDRO_31;
DENDRO_74 = -0.5;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_74;
DENDRO_75 *= DENDRO_62;
DENDRO_75 *= DENDRO_28;
DENDRO_28 = 0.5;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_28;
DENDRO_74 *= DENDRO_36;
DENDRO_74 *= DENDRO_16;
DENDRO_74 *= DENDRO_13;
DENDRO_28 = 0.5;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_28;
DENDRO_36 *= DENDRO_45;
DENDRO_36 *= DENDRO_16;
DENDRO_36 *= DENDRO_9;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_16;
DENDRO_28 *= DENDRO_7;
DENDRO_28 *= DENDRO_31;
DENDRO_31 = -0.5;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_31;
DENDRO_45 *= DENDRO_62;
DENDRO_45 *= DENDRO_22;
DENDRO_22 = 0.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_22;
DENDRO_31 *= DENDRO_47;
DENDRO_31 *= DENDRO_16;
DENDRO_31 *= DENDRO_17;
DENDRO_22 = 0.5;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_22;
DENDRO_76 *= DENDRO_53;
DENDRO_76 *= DENDRO_16;
DENDRO_76 *= DENDRO_21;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_16;
DENDRO_22 *= DENDRO_7;
DENDRO_22 *= DENDRO_23;
DENDRO_77 = -0.5;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_77;
DENDRO_78 *= DENDRO_62;
DENDRO_78 *= DENDRO_25;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_25;
DENDRO_77 *= DENDRO_47;
DENDRO_77 *= DENDRO_16;
DENDRO_77 *= DENDRO_14;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_25;
DENDRO_79 *= DENDRO_53;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_17;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_16;
DENDRO_25 *= DENDRO_9;
DENDRO_25 *= DENDRO_23;
DENDRO_80 = -0.5;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_80;
DENDRO_81 *= DENDRO_62;
DENDRO_81 *= DENDRO_19;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_19;
DENDRO_80 *= DENDRO_4;
DENDRO_80 *= DENDRO_62;
DENDRO_80 *= DENDRO_18;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_19;
DENDRO_82 *= DENDRO_47;
DENDRO_82 *= DENDRO_16;
DENDRO_82 *= DENDRO_9;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_19;
DENDRO_47 *= DENDRO_53;
DENDRO_47 *= DENDRO_16;
DENDRO_47 *= DENDRO_7;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_16;
DENDRO_19 *= DENDRO_13;
DENDRO_19 *= DENDRO_23;
DENDRO_23 = 0.5;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_23;
DENDRO_53 *= DENDRO_10;
DENDRO_53 *= DENDRO_62;
DENDRO_53 *= DENDRO_30;
DENDRO_23 = 0.5;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_23;
DENDRO_23 = grad_0_gt0;
DENDRO_83 *= DENDRO_23;
DENDRO_83 *= DENDRO_16;
DENDRO_83 *= DENDRO_7;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_16;
DENDRO_84 *= DENDRO_15;
DENDRO_84 *= DENDRO_17;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_16;
DENDRO_85 *= DENDRO_21;
DENDRO_85 *= DENDRO_26;
DENDRO_86 = 0.5;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_86;
DENDRO_87 *= DENDRO_10;
DENDRO_87 *= DENDRO_62;
DENDRO_87 *= DENDRO_20;
DENDRO_86 = 0.5;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_86;
DENDRO_88 *= DENDRO_23;
DENDRO_88 *= DENDRO_16;
DENDRO_88 *= DENDRO_9;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_16;
DENDRO_86 *= DENDRO_26;
DENDRO_86 *= DENDRO_17;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_16;
DENDRO_89 *= DENDRO_14;
DENDRO_89 *= DENDRO_15;
DENDRO_90 = 0.5;

// initialize reduction 
DENDRO_91 = 1;
DENDRO_91 *= DENDRO_90;
DENDRO_91 *= DENDRO_23;
DENDRO_91 *= DENDRO_16;
DENDRO_91 *= DENDRO_13;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_16;
DENDRO_23 *= DENDRO_26;
DENDRO_23 *= DENDRO_7;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_16;
DENDRO_26 *= DENDRO_15;
DENDRO_26 *= DENDRO_9;
DENDRO_15 = -0.5;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_15;
DENDRO_90 *= DENDRO_62;
DENDRO_90 *= DENDRO_27;
DENDRO_15 = 0.5;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_15;
DENDRO_27 *= DENDRO_11;
DENDRO_27 *= DENDRO_62;
DENDRO_27 *= DENDRO_30;
DENDRO_15 = 0.5;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_15;
DENDRO_15 = grad_1_gt3;
DENDRO_30 *= DENDRO_15;
DENDRO_30 *= DENDRO_16;
DENDRO_30 *= DENDRO_17;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_16;
DENDRO_92 *= DENDRO_32;
DENDRO_92 *= DENDRO_7;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_16;
DENDRO_93 *= DENDRO_21;
DENDRO_93 *= DENDRO_33;
DENDRO_94 = 0.5;

// initialize reduction 
DENDRO_95 = 1;
DENDRO_95 *= DENDRO_94;
DENDRO_95 *= DENDRO_15;
DENDRO_95 *= DENDRO_16;
DENDRO_95 *= DENDRO_14;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_16;
DENDRO_94 *= DENDRO_33;
DENDRO_94 *= DENDRO_17;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_16;
DENDRO_96 *= DENDRO_32;
DENDRO_96 *= DENDRO_9;
DENDRO_97 = -0.5;

// initialize reduction 
DENDRO_98 = 1;
DENDRO_98 *= DENDRO_97;
DENDRO_98 *= DENDRO_62;
DENDRO_98 *= DENDRO_34;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_97 = 1;
DENDRO_97 *= DENDRO_34;
DENDRO_97 *= DENDRO_11;
DENDRO_97 *= DENDRO_62;
DENDRO_97 *= DENDRO_18;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_34;
DENDRO_99 *= DENDRO_15;
DENDRO_99 *= DENDRO_16;
DENDRO_99 *= DENDRO_9;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_16;
DENDRO_15 *= DENDRO_33;
DENDRO_15 *= DENDRO_7;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_16;
DENDRO_33 *= DENDRO_13;
DENDRO_33 *= DENDRO_32;
DENDRO_32 = 0.5;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_32;
DENDRO_32 = grad_2_gt5;
DENDRO_34 *= DENDRO_32;
DENDRO_34 *= DENDRO_16;
DENDRO_34 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_16;
DENDRO_21 *= DENDRO_35;
DENDRO_21 *= DENDRO_17;

// initialize reduction 
DENDRO_100 = 1;
DENDRO_100 *= DENDRO_16;
DENDRO_100 *= DENDRO_37;
DENDRO_100 *= DENDRO_7;
DENDRO_101 = -0.5;

// initialize reduction 
DENDRO_102 = 1;
DENDRO_102 *= DENDRO_101;
DENDRO_102 *= DENDRO_62;
DENDRO_102 *= DENDRO_39;
DENDRO_39 = 0.5;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_39;
DENDRO_101 *= DENDRO_12;
DENDRO_101 *= DENDRO_62;
DENDRO_101 *= DENDRO_20;
DENDRO_20 = 0.5;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_20;
DENDRO_39 *= DENDRO_32;
DENDRO_39 *= DENDRO_16;
DENDRO_39 *= DENDRO_17;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_16;
DENDRO_17 *= DENDRO_37;
DENDRO_17 *= DENDRO_9;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_16;
DENDRO_20 *= DENDRO_14;
DENDRO_20 *= DENDRO_35;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_103 = 1;
DENDRO_103 *= DENDRO_14;
DENDRO_103 *= DENDRO_12;
DENDRO_103 *= DENDRO_62;
DENDRO_103 *= DENDRO_18;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_14;
DENDRO_18 *= DENDRO_32;
DENDRO_18 *= DENDRO_16;
DENDRO_18 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_16;
DENDRO_7 *= DENDRO_35;
DENDRO_7 *= DENDRO_9;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_16;
DENDRO_9 *= DENDRO_13;
DENDRO_9 *= DENDRO_37;
DENDRO_13 = At5[pp];

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_13;
DENDRO_14 *= DENDRO_40;
DENDRO_14 *= DENDRO_49;
DENDRO_32 = At4[pp];

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_32;
DENDRO_35 *= DENDRO_51;
DENDRO_35 *= DENDRO_49;
DENDRO_37 = At3[pp];

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_37;
DENDRO_62 *= DENDRO_51;
DENDRO_62 *= DENDRO_40;
DENDRO_104 = At2[pp];

// initialize reduction 
DENDRO_105 = 1;
DENDRO_105 *= DENDRO_104;
DENDRO_105 *= DENDRO_42;
DENDRO_105 *= DENDRO_40;

// initialize reduction 
DENDRO_106 = 1;
DENDRO_106 *= DENDRO_104;
DENDRO_106 *= DENDRO_46;
DENDRO_106 *= DENDRO_49;
DENDRO_107 = At1[pp];

// initialize reduction 
DENDRO_108 = 1;
DENDRO_108 *= DENDRO_107;
DENDRO_108 *= DENDRO_42;
DENDRO_108 *= DENDRO_51;

// initialize reduction 
DENDRO_109 = 1;
DENDRO_109 *= DENDRO_107;
DENDRO_109 *= DENDRO_46;
DENDRO_109 *= DENDRO_40;
DENDRO_110 = At0[pp];

// initialize reduction 
DENDRO_111 = 1;
DENDRO_111 *= DENDRO_110;
DENDRO_111 *= DENDRO_46;
DENDRO_111 *= DENDRO_42;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_32;
DENDRO_112 *= DENDRO_41;

// initialize reduction 
DENDRO_113 = 1;
DENDRO_113 *= DENDRO_13;
DENDRO_113 *= DENDRO_42;
DENDRO_113 *= DENDRO_49;

// initialize reduction 
DENDRO_114 = 1;
DENDRO_114 *= DENDRO_32;
DENDRO_114 *= DENDRO_42;
DENDRO_114 *= DENDRO_40;

// initialize reduction 
DENDRO_115 = 1;
DENDRO_115 *= DENDRO_32;
DENDRO_115 *= DENDRO_46;
DENDRO_115 *= DENDRO_49;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_37;
DENDRO_116 *= DENDRO_46;
DENDRO_116 *= DENDRO_40;

// initialize reduction 
DENDRO_117 = 1;
DENDRO_117 *= DENDRO_104;
DENDRO_117 *= DENDRO_54;
DENDRO_117 *= DENDRO_49;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_107;
DENDRO_118 *= DENDRO_46;
DENDRO_118 *= DENDRO_42;

// initialize reduction 
DENDRO_119 = 1;
DENDRO_119 *= DENDRO_107;
DENDRO_119 *= DENDRO_54;
DENDRO_119 *= DENDRO_40;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_110;
DENDRO_120 *= DENDRO_54;
DENDRO_120 *= DENDRO_42;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_104;
DENDRO_121 *= DENDRO_44;

// initialize reduction 
DENDRO_122 = 1;
DENDRO_122 *= DENDRO_13;
DENDRO_122 *= DENDRO_42;
DENDRO_122 *= DENDRO_40;

// initialize reduction 
DENDRO_123 = 1;
DENDRO_123 *= DENDRO_32;
DENDRO_123 *= DENDRO_42;
DENDRO_123 *= DENDRO_51;

// initialize reduction 
DENDRO_124 = 1;
DENDRO_124 *= DENDRO_32;
DENDRO_124 *= DENDRO_46;
DENDRO_124 *= DENDRO_40;

// initialize reduction 
DENDRO_125 = 1;
DENDRO_125 *= DENDRO_37;
DENDRO_125 *= DENDRO_46;
DENDRO_125 *= DENDRO_51;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_104;
DENDRO_126 *= DENDRO_46;
DENDRO_126 *= DENDRO_42;

// initialize reduction 
DENDRO_127 = 1;
DENDRO_127 *= DENDRO_104;
DENDRO_127 *= DENDRO_54;
DENDRO_127 *= DENDRO_40;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_107;
DENDRO_128 *= DENDRO_54;
DENDRO_128 *= DENDRO_51;

// initialize reduction 
DENDRO_129 = 1;
DENDRO_129 *= DENDRO_110;
DENDRO_129 *= DENDRO_54;
DENDRO_129 *= DENDRO_46;

// initialize reduction 
DENDRO_130 = 1;
DENDRO_130 *= DENDRO_107;
DENDRO_130 *= DENDRO_48;

// initialize reduction 
DENDRO_131 = 1;
DENDRO_131 *= DENDRO_43;
DENDRO_131 *= DENDRO_32;
DENDRO_131 *= DENDRO_40;
DENDRO_131 *= DENDRO_49;

// initialize reduction 
DENDRO_132 = 1;
DENDRO_132 *= DENDRO_43;
DENDRO_132 *= DENDRO_104;
DENDRO_132 *= DENDRO_42;
DENDRO_132 *= DENDRO_49;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_43;
DENDRO_49 *= DENDRO_107;
DENDRO_49 *= DENDRO_42;
DENDRO_49 *= DENDRO_40;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_13;
DENDRO_133 *= DENDRO_50;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_37;
DENDRO_50 *= DENDRO_41;

// initialize reduction 
DENDRO_134 = 1;
DENDRO_134 *= DENDRO_110;
DENDRO_134 *= DENDRO_44;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_43;
DENDRO_135 *= DENDRO_32;
DENDRO_135 *= DENDRO_51;
DENDRO_135 *= DENDRO_40;

// initialize reduction 
DENDRO_136 = 1;
DENDRO_136 *= DENDRO_43;
DENDRO_136 *= DENDRO_104;
DENDRO_136 *= DENDRO_46;
DENDRO_136 *= DENDRO_40;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_43;
DENDRO_40 *= DENDRO_107;
DENDRO_40 *= DENDRO_46;
DENDRO_40 *= DENDRO_51;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_13;
DENDRO_51 *= DENDRO_41;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_37;
DENDRO_41 *= DENDRO_52;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_110;
DENDRO_52 *= DENDRO_48;

// initialize reduction 
DENDRO_137 = 1;
DENDRO_137 *= DENDRO_43;
DENDRO_137 *= DENDRO_32;
DENDRO_137 *= DENDRO_46;
DENDRO_137 *= DENDRO_42;

// initialize reduction 
DENDRO_138 = 1;
DENDRO_138 *= DENDRO_43;
DENDRO_138 *= DENDRO_104;
DENDRO_138 *= DENDRO_54;
DENDRO_138 *= DENDRO_42;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_43;
DENDRO_42 *= DENDRO_107;
DENDRO_42 *= DENDRO_54;
DENDRO_42 *= DENDRO_46;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_13;
DENDRO_46 *= DENDRO_44;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_37;
DENDRO_44 *= DENDRO_48;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_110;
DENDRO_48 *= DENDRO_56;

// initialize reduction 
DENDRO_54 = 0;
DENDRO_54 += DENDRO_64;
DENDRO_54 += DENDRO_57;
DENDRO_54 += DENDRO_59;
DENDRO_54 += DENDRO_58;

// initialize reduction 
DENDRO_56 = 0;
DENDRO_56 += DENDRO_24;
DENDRO_56 += DENDRO_66;
DENDRO_56 += DENDRO_65;
DENDRO_56 += DENDRO_60;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_61;
DENDRO_24 += DENDRO_38;
DENDRO_24 += DENDRO_67;
DENDRO_24 += DENDRO_68;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_43;
DENDRO_38 *= DENDRO_0;
DENDRO_38 *= DENDRO_2;
DENDRO_38 *= DENDRO_4;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_10;
DENDRO_57 *= DENDRO_11;
DENDRO_57 *= DENDRO_12;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_8;
DENDRO_58 *= DENDRO_12;
DENDRO_58 *= DENDRO_1;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_8;
DENDRO_59 *= DENDRO_11;
DENDRO_59 *= DENDRO_3;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_8;
DENDRO_60 *= DENDRO_10;
DENDRO_60 *= DENDRO_5;

// initialize reduction 
DENDRO_61 = 0;
DENDRO_61 += DENDRO_55;
DENDRO_61 += DENDRO_70;
DENDRO_61 += DENDRO_69;
DENDRO_61 += DENDRO_63;

// initialize reduction 
DENDRO_55 = 0;
DENDRO_55 += DENDRO_75;
DENDRO_55 += DENDRO_71;
DENDRO_55 += DENDRO_73;
DENDRO_55 += DENDRO_72;

// initialize reduction 
DENDRO_63 = 0;
DENDRO_63 += DENDRO_45;
DENDRO_63 += DENDRO_28;
DENDRO_63 += DENDRO_36;
DENDRO_63 += DENDRO_74;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_78;
DENDRO_28 += DENDRO_22;
DENDRO_28 += DENDRO_76;
DENDRO_28 += DENDRO_31;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_81;
DENDRO_22 += DENDRO_25;
DENDRO_22 += DENDRO_79;
DENDRO_22 += DENDRO_77;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_19;
DENDRO_25 += DENDRO_47;
DENDRO_25 += DENDRO_82;
DENDRO_25 += DENDRO_80;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_85;
DENDRO_19 += DENDRO_84;
DENDRO_19 += DENDRO_83;
DENDRO_19 += DENDRO_53;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_89;
DENDRO_31 += DENDRO_86;
DENDRO_31 += DENDRO_88;
DENDRO_31 += DENDRO_87;

// initialize reduction 
DENDRO_36 = 0;
DENDRO_36 += DENDRO_90;
DENDRO_36 += DENDRO_26;
DENDRO_36 += DENDRO_23;
DENDRO_36 += DENDRO_91;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_93;
DENDRO_23 += DENDRO_92;
DENDRO_23 += DENDRO_30;
DENDRO_23 += DENDRO_27;

// initialize reduction 
DENDRO_26 = 0;
DENDRO_26 += DENDRO_98;
DENDRO_26 += DENDRO_96;
DENDRO_26 += DENDRO_94;
DENDRO_26 += DENDRO_95;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_33;
DENDRO_27 += DENDRO_15;
DENDRO_27 += DENDRO_99;
DENDRO_27 += DENDRO_97;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_102;
DENDRO_15 += DENDRO_100;
DENDRO_15 += DENDRO_21;
DENDRO_15 += DENDRO_34;

// initialize reduction 
DENDRO_21 = 0;
DENDRO_21 += DENDRO_20;
DENDRO_21 += DENDRO_17;
DENDRO_21 += DENDRO_39;
DENDRO_21 += DENDRO_101;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_9;
DENDRO_17 += DENDRO_7;
DENDRO_17 += DENDRO_18;
DENDRO_17 += DENDRO_103;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_112;
DENDRO_7 += DENDRO_111;
DENDRO_7 += DENDRO_109;
DENDRO_7 += DENDRO_108;
DENDRO_7 += DENDRO_106;
DENDRO_7 += DENDRO_105;
DENDRO_7 += DENDRO_62;
DENDRO_7 += DENDRO_35;
DENDRO_7 += DENDRO_14;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_121;
DENDRO_9 += DENDRO_120;
DENDRO_9 += DENDRO_119;
DENDRO_9 += DENDRO_118;
DENDRO_9 += DENDRO_117;
DENDRO_9 += DENDRO_116;
DENDRO_9 += DENDRO_115;
DENDRO_9 += DENDRO_114;
DENDRO_9 += DENDRO_113;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_130;
DENDRO_14 += DENDRO_129;
DENDRO_14 += DENDRO_128;
DENDRO_14 += DENDRO_127;
DENDRO_14 += DENDRO_126;
DENDRO_14 += DENDRO_125;
DENDRO_14 += DENDRO_124;
DENDRO_14 += DENDRO_123;
DENDRO_14 += DENDRO_122;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_134;
DENDRO_18 += DENDRO_50;
DENDRO_18 += DENDRO_133;
DENDRO_18 += DENDRO_49;
DENDRO_18 += DENDRO_132;
DENDRO_18 += DENDRO_131;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_52;
DENDRO_20 += DENDRO_41;
DENDRO_20 += DENDRO_51;
DENDRO_20 += DENDRO_40;
DENDRO_20 += DENDRO_136;
DENDRO_20 += DENDRO_135;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_48;
DENDRO_30 += DENDRO_44;
DENDRO_30 += DENDRO_46;
DENDRO_30 += DENDRO_42;
DENDRO_30 += DENDRO_138;
DENDRO_30 += DENDRO_137;
DENDRO_33 = K[pp];

// int-power reduce pow(K[pp],2)
DENDRO_34 = DENDRO_33;
DENDRO_34 *= DENDRO_33;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_8;
DENDRO_35 = grad_2_alpha;
DENDRO_33 *= DENDRO_35;
DENDRO_33 *= DENDRO_54;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_8;
DENDRO_40 = grad_1_alpha;
DENDRO_39 *= DENDRO_40;
DENDRO_39 *= DENDRO_56;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_8;
DENDRO_42 = grad_0_alpha;
DENDRO_41 *= DENDRO_42;
DENDRO_41 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_8;
DENDRO_24 *= DENDRO_2;
DENDRO_24 *= DENDRO_11;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_0;
DENDRO_44 *= DENDRO_4;

// initialize reduction 
DENDRO_45 = 0;
DENDRO_45 += DENDRO_60;
DENDRO_45 += DENDRO_59;
DENDRO_45 += DENDRO_58;
DENDRO_45 += DENDRO_57;
DENDRO_45 += DENDRO_38;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_8;
DENDRO_38 *= DENDRO_35;
DENDRO_38 *= DENDRO_61;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_8;
DENDRO_46 *= DENDRO_40;
DENDRO_46 *= DENDRO_55;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_8;
DENDRO_47 *= DENDRO_42;
DENDRO_47 *= DENDRO_63;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_8;
DENDRO_48 *= DENDRO_2;
DENDRO_48 *= DENDRO_4;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_0;
DENDRO_49 *= DENDRO_12;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_8;
DENDRO_50 *= DENDRO_35;
DENDRO_50 *= DENDRO_28;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_8;
DENDRO_28 *= DENDRO_40;
DENDRO_28 *= DENDRO_22;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_8;
DENDRO_22 *= DENDRO_42;
DENDRO_22 *= DENDRO_25;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_8;
DENDRO_25 *= DENDRO_0;
DENDRO_25 *= DENDRO_2;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_10;
DENDRO_0 *= DENDRO_4;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_8;
DENDRO_2 *= DENDRO_35;
DENDRO_2 *= DENDRO_19;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_8;
DENDRO_4 *= DENDRO_40;
DENDRO_4 *= DENDRO_31;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_8;
DENDRO_19 *= DENDRO_42;
DENDRO_19 *= DENDRO_36;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_11;
DENDRO_31 *= DENDRO_12;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_8;
DENDRO_36 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_8;
DENDRO_5 *= DENDRO_35;
DENDRO_5 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_8;
DENDRO_23 *= DENDRO_40;
DENDRO_23 *= DENDRO_26;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_8;
DENDRO_26 *= DENDRO_42;
DENDRO_26 *= DENDRO_27;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_10;
DENDRO_27 *= DENDRO_12;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_8;
DENDRO_12 *= DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_8;
DENDRO_3 *= DENDRO_35;
DENDRO_3 *= DENDRO_15;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_8;
DENDRO_15 *= DENDRO_40;
DENDRO_15 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_8;
DENDRO_21 *= DENDRO_42;
DENDRO_21 *= DENDRO_17;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_10;
DENDRO_17 *= DENDRO_11;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_8;
DENDRO_10 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_43;
DENDRO_1 *= DENDRO_32;
DENDRO_1 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_43;
DENDRO_7 *= DENDRO_104;
DENDRO_7 *= DENDRO_9;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_43;
DENDRO_9 *= DENDRO_107;
DENDRO_9 *= DENDRO_14;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_13;
DENDRO_11 *= DENDRO_18;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_37;
DENDRO_13 *= DENDRO_20;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_110;
DENDRO_14 *= DENDRO_30;
DENDRO_18 = 1.0/3.0;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_18;
DENDRO_20 *= DENDRO_34;
DENDRO_18 = grad2_0_2_alpha;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_18;
DENDRO_30 += DENDRO_41;
DENDRO_30 += DENDRO_39;
DENDRO_30 += DENDRO_33;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_44;
DENDRO_18 += DENDRO_24;

// int-power reduce pow(gt0[pp]*gt3[pp]*gt5[pp] - gt0[pp]*gt4[pp]**2 - gt1[pp]**2*gt5[pp] + 2*gt1[pp]*gt2[pp]*gt4[pp] - gt2[pp]**2*gt3[pp],-1)
DENDRO_24 = DENDRO_45;
DENDRO_24 = 1.0/DENDRO_24;
DENDRO_32 = grad2_0_1_alpha;

// initialize reduction 
DENDRO_33 = 0;
DENDRO_33 += DENDRO_32;
DENDRO_33 += DENDRO_47;
DENDRO_33 += DENDRO_46;
DENDRO_33 += DENDRO_38;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_49;
DENDRO_32 += DENDRO_48;
DENDRO_34 = grad2_1_2_alpha;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_34;
DENDRO_35 += DENDRO_22;
DENDRO_35 += DENDRO_28;
DENDRO_35 += DENDRO_50;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_0;
DENDRO_22 += DENDRO_25;
DENDRO_0 = grad2_0_0_alpha;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_0;
DENDRO_25 += DENDRO_19;
DENDRO_25 += DENDRO_4;
DENDRO_25 += DENDRO_2;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_36;
DENDRO_0 += DENDRO_31;
DENDRO_2 = grad2_1_1_alpha;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_2;
DENDRO_4 += DENDRO_26;
DENDRO_4 += DENDRO_23;
DENDRO_4 += DENDRO_5;

// initialize reduction 
DENDRO_2 = 0;
DENDRO_2 += DENDRO_12;
DENDRO_2 += DENDRO_27;
DENDRO_5 = grad2_2_2_alpha;

// initialize reduction 
DENDRO_12 = 0;
DENDRO_12 += DENDRO_5;
DENDRO_12 += DENDRO_21;
DENDRO_12 += DENDRO_15;
DENDRO_12 += DENDRO_3;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_10;
DENDRO_3 += DENDRO_17;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_20;
DENDRO_5 += DENDRO_14;
DENDRO_5 += DENDRO_13;
DENDRO_5 += DENDRO_11;
DENDRO_5 += DENDRO_9;
DENDRO_5 += DENDRO_7;
DENDRO_5 += DENDRO_1;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_6;
DENDRO_1 *= DENDRO_29;
DENDRO_1 *= DENDRO_24;
DENDRO_1 *= DENDRO_18;
DENDRO_1 *= DENDRO_30;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_29;
DENDRO_7 *= DENDRO_16;
DENDRO_7 *= DENDRO_32;
DENDRO_7 *= DENDRO_33;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_6;
DENDRO_9 *= DENDRO_29;
DENDRO_9 *= DENDRO_16;
DENDRO_9 *= DENDRO_22;
DENDRO_9 *= DENDRO_35;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_8;
DENDRO_6 *= DENDRO_29;
DENDRO_6 *= DENDRO_24;
DENDRO_6 *= DENDRO_0;
DENDRO_6 *= DENDRO_25;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_8;
DENDRO_0 *= DENDRO_29;
DENDRO_0 *= DENDRO_24;
DENDRO_0 *= DENDRO_2;
DENDRO_0 *= DENDRO_4;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_8;
DENDRO_2 *= DENDRO_29;
DENDRO_2 *= DENDRO_24;
DENDRO_2 *= DENDRO_3;
DENDRO_2 *= DENDRO_12;
DENDRO_3 = beta2[pp];

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_3;
DENDRO_3 = grad_2_K;
DENDRO_4 *= DENDRO_3;
DENDRO_3 = beta1[pp];

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_3;
DENDRO_3 = grad_1_K;
DENDRO_8 *= DENDRO_3;
DENDRO_3 = beta0[pp];

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_3;
DENDRO_3 = grad_0_K;
DENDRO_10 *= DENDRO_3;
DENDRO_3 = alpha[pp];

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_3;
DENDRO_11 *= DENDRO_5;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_11;
DENDRO_3 += DENDRO_10;
DENDRO_3 += DENDRO_8;
DENDRO_3 += DENDRO_4;
DENDRO_3 += DENDRO_2;
DENDRO_3 += DENDRO_0;
DENDRO_3 += DENDRO_6;
DENDRO_3 += DENDRO_9;
DENDRO_3 += DENDRO_7;
DENDRO_3 += DENDRO_1;
K_rhs[pp] = DENDRO_3;


}
if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&K_rhs[gidx],     K[gidx]  ,   grad_0_K,   grad_1_K,   grad_2_K, 1.0, 0.0, blk);
}

K_rhs[pp]    += ko_sigma * (kograd_0_K + kograd_1_K + kograd_2_K);
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , Gt0, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_Gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_Gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_Gt0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_Gt0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_Gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_Gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_Gt0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_Gt0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_Gt0=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , Gt0, blk);
device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_Gt0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_Gt0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_Gt0=Du[gidx];
__syncthreads();

//|V|= 396 |E| =905
//topological generations = 10
// [0]
// traversal id=0 temp_var_count=106
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;double DENDRO_33;double DENDRO_34;double DENDRO_35;double DENDRO_36;double DENDRO_37;double DENDRO_38;double DENDRO_39;double DENDRO_40;double DENDRO_41;double DENDRO_42;double DENDRO_43;double DENDRO_44;double DENDRO_45;double DENDRO_46;double DENDRO_47;double DENDRO_48;double DENDRO_49;double DENDRO_50;double DENDRO_51;double DENDRO_52;double DENDRO_53;double DENDRO_54;double DENDRO_55;double DENDRO_56;double DENDRO_57;double DENDRO_58;double DENDRO_59;double DENDRO_60;double DENDRO_61;double DENDRO_62;double DENDRO_63;double DENDRO_64;double DENDRO_65;double DENDRO_66;double DENDRO_67;double DENDRO_68;double DENDRO_69;double DENDRO_70;double DENDRO_71;double DENDRO_72;double DENDRO_73;double DENDRO_74;double DENDRO_75;double DENDRO_76;double DENDRO_77;double DENDRO_78;double DENDRO_79;double DENDRO_80;double DENDRO_81;double DENDRO_82;double DENDRO_83;double DENDRO_84;double DENDRO_85;double DENDRO_86;double DENDRO_87;double DENDRO_88;double DENDRO_89;double DENDRO_90;double DENDRO_91;double DENDRO_92;double DENDRO_93;double DENDRO_94;double DENDRO_95;double DENDRO_96;double DENDRO_97;double DENDRO_98;double DENDRO_99;double DENDRO_100;double DENDRO_101;double DENDRO_102;double DENDRO_103;double DENDRO_104;double DENDRO_105;
DENDRO_0 = gt1[pp];

// int-power reduce pow(gt1[pp],2)
DENDRO_1 = DENDRO_0;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = gt2[pp];

// int-power reduce pow(gt2[pp],2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = gt4[pp];

// int-power reduce pow(gt4[pp],2)
DENDRO_5 = DENDRO_4;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = -2;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_4;
DENDRO_6 = -1;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_6;
DENDRO_9 = gt0[pp];
DENDRO_8 *= DENDRO_9;
DENDRO_10 = gt3[pp];
DENDRO_8 *= DENDRO_10;
DENDRO_11 = gt5[pp];
DENDRO_8 *= DENDRO_11;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_11;
DENDRO_12 *= DENDRO_1;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_10;
DENDRO_13 *= DENDRO_3;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_9;
DENDRO_14 *= DENDRO_5;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_6;
DENDRO_15 *= DENDRO_10;
DENDRO_15 *= DENDRO_11;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_14;
DENDRO_16 += DENDRO_13;
DENDRO_16 += DENDRO_12;
DENDRO_16 += DENDRO_8;
DENDRO_16 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_4;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_2;
DENDRO_8 *= DENDRO_10;
DENDRO_12 = -0.5;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_12;
DENDRO_12 = grad_1_gt2;
DENDRO_13 *= DENDRO_12;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_14;
DENDRO_14 = grad_2_gt1;
DENDRO_17 *= DENDRO_14;
DENDRO_18 = 0.5;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_18;
DENDRO_18 = grad_0_gt4;
DENDRO_19 *= DENDRO_18;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_6;
DENDRO_20 *= DENDRO_2;
DENDRO_20 *= DENDRO_4;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_0;
DENDRO_21 *= DENDRO_11;
DENDRO_22 = -0.5;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_22;
DENDRO_23 *= DENDRO_14;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_14;
DENDRO_22 *= DENDRO_12;
DENDRO_12 = -0.5;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_12;
DENDRO_14 *= DENDRO_18;
DENDRO_12 = -0.5;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_12;
DENDRO_12 = grad_2_gt0;
DENDRO_18 *= DENDRO_12;
DENDRO_24 = 1.0;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_24;
DENDRO_24 = grad_0_gt2;
DENDRO_25 *= DENDRO_24;
DENDRO_24 = -0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_24;
DENDRO_24 = grad_1_gt0;
DENDRO_26 *= DENDRO_24;
DENDRO_27 = 1.0;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_27;
DENDRO_27 = grad_0_gt1;
DENDRO_28 *= DENDRO_27;
DENDRO_27 = -0.5;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_27;
DENDRO_27 = grad_2_gt3;
DENDRO_29 *= DENDRO_27;
DENDRO_30 = 1.0;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_30;
DENDRO_30 = grad_1_gt4;
DENDRO_31 *= DENDRO_30;
DENDRO_30 = -0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_30;
DENDRO_30 = grad_0_gt3;
DENDRO_32 *= DENDRO_30;
DENDRO_33 = 1.0;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_33;
DENDRO_33 = grad_1_gt1;
DENDRO_34 *= DENDRO_33;
DENDRO_33 = -0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_33;
DENDRO_33 = grad_1_gt5;
DENDRO_35 *= DENDRO_33;
DENDRO_36 = 1.0;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = grad_2_gt4;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = -0.5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_36;
DENDRO_36 = grad_0_gt5;
DENDRO_38 *= DENDRO_36;
DENDRO_39 = 1.0;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_39;
DENDRO_39 = grad_2_gt2;
DENDRO_40 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_6;
DENDRO_39 *= DENDRO_9;
DENDRO_39 *= DENDRO_10;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_6;
DENDRO_10 *= DENDRO_0;
DENDRO_10 *= DENDRO_2;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_9;
DENDRO_0 *= DENDRO_4;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_6;
DENDRO_2 *= DENDRO_9;
DENDRO_2 *= DENDRO_11;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_5;
DENDRO_4 += DENDRO_15;

// int-power reduce pow(-gt0[pp]*gt3[pp]*gt5[pp] + gt0[pp]*gt4[pp]**2 + gt1[pp]**2*gt5[pp] - 2*gt1[pp]*gt2[pp]*gt4[pp] + gt2[pp]**2*gt3[pp],-1)
DENDRO_5 = DENDRO_16;
DENDRO_5 = 1.0/DENDRO_5;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_8;
DENDRO_9 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_19;
DENDRO_7 += DENDRO_17;
DENDRO_7 += DENDRO_13;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_21;
DENDRO_8 += DENDRO_20;

// initialize reduction 
DENDRO_11 = 0;
DENDRO_11 += DENDRO_19;
DENDRO_11 += DENDRO_22;
DENDRO_11 += DENDRO_23;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_22;
DENDRO_13 += DENDRO_17;
DENDRO_13 += DENDRO_14;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_25;
DENDRO_14 += DENDRO_18;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_28;
DENDRO_15 += DENDRO_26;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_31;
DENDRO_16 += DENDRO_29;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_34;
DENDRO_17 += DENDRO_32;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_37;
DENDRO_18 += DENDRO_35;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_40;
DENDRO_19 += DENDRO_38;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_1;
DENDRO_20 += DENDRO_39;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_10;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_3;
DENDRO_0 += DENDRO_2;
DENDRO_2 = DENDRO_igt2;

// int-power reduce pow(DENDRO_igt2,2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_10 = DENDRO_igt1;

// int-power reduce pow(DENDRO_igt1,2)
DENDRO_21 = DENDRO_10;
DENDRO_21 *= DENDRO_10;
DENDRO_22 = DENDRO_igt0;

// int-power reduce pow(DENDRO_igt0,2)
DENDRO_23 = DENDRO_22;
DENDRO_23 *= DENDRO_22;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_25;
DENDRO_26 *= DENDRO_12;
DENDRO_26 *= DENDRO_5;
DENDRO_26 *= DENDRO_4;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_25;
DENDRO_28 *= DENDRO_36;
DENDRO_28 *= DENDRO_5;
DENDRO_28 *= DENDRO_9;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_5;
DENDRO_25 *= DENDRO_8;
DENDRO_25 *= DENDRO_7;
DENDRO_29 = 0.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_29;
DENDRO_31 *= DENDRO_24;
DENDRO_31 *= DENDRO_5;
DENDRO_31 *= DENDRO_4;
DENDRO_29 = 0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_29;
DENDRO_32 *= DENDRO_30;
DENDRO_32 *= DENDRO_5;
DENDRO_32 *= DENDRO_8;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_5;
DENDRO_29 *= DENDRO_9;
DENDRO_29 *= DENDRO_11;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_34;
DENDRO_35 *= DENDRO_27;
DENDRO_35 *= DENDRO_5;
DENDRO_35 *= DENDRO_8;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_34;
DENDRO_37 *= DENDRO_33;
DENDRO_37 *= DENDRO_5;
DENDRO_37 *= DENDRO_9;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_5;
DENDRO_34 *= DENDRO_4;
DENDRO_34 *= DENDRO_13;
DENDRO_38 = 0.5;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_38;
DENDRO_38 = grad_0_gt0;
DENDRO_39 *= DENDRO_38;
DENDRO_39 *= DENDRO_5;
DENDRO_39 *= DENDRO_4;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_5;
DENDRO_40 *= DENDRO_14;
DENDRO_40 *= DENDRO_9;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_5;
DENDRO_41 *= DENDRO_15;
DENDRO_41 *= DENDRO_8;
DENDRO_42 = 0.5;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_42;
DENDRO_42 = grad_1_gt3;
DENDRO_43 *= DENDRO_42;
DENDRO_43 *= DENDRO_5;
DENDRO_43 *= DENDRO_8;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_5;
DENDRO_44 *= DENDRO_16;
DENDRO_44 *= DENDRO_9;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_5;
DENDRO_45 *= DENDRO_4;
DENDRO_45 *= DENDRO_17;
DENDRO_46 = 0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_46;
DENDRO_46 = grad_2_gt5;
DENDRO_47 *= DENDRO_46;
DENDRO_47 *= DENDRO_5;
DENDRO_47 *= DENDRO_9;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_5;
DENDRO_48 *= DENDRO_18;
DENDRO_48 *= DENDRO_8;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_5;
DENDRO_49 *= DENDRO_4;
DENDRO_49 *= DENDRO_19;
DENDRO_50 = 0.5;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_50;
DENDRO_51 *= DENDRO_12;
DENDRO_51 *= DENDRO_5;
DENDRO_51 *= DENDRO_9;
DENDRO_50 = 0.5;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_50;
DENDRO_52 *= DENDRO_36;
DENDRO_52 *= DENDRO_5;
DENDRO_52 *= DENDRO_20;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_5;
DENDRO_50 *= DENDRO_1;
DENDRO_50 *= DENDRO_7;
DENDRO_53 = 0.5;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_53;
DENDRO_54 *= DENDRO_24;
DENDRO_54 *= DENDRO_5;
DENDRO_54 *= DENDRO_9;
DENDRO_53 = 0.5;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_53;
DENDRO_55 *= DENDRO_30;
DENDRO_55 *= DENDRO_5;
DENDRO_55 *= DENDRO_1;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_5;
DENDRO_53 *= DENDRO_20;
DENDRO_53 *= DENDRO_11;
DENDRO_56 = 0.5;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_56;
DENDRO_57 *= DENDRO_27;
DENDRO_57 *= DENDRO_5;
DENDRO_57 *= DENDRO_1;
DENDRO_56 = 0.5;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_56;
DENDRO_58 *= DENDRO_33;
DENDRO_58 *= DENDRO_5;
DENDRO_58 *= DENDRO_20;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_5;
DENDRO_56 *= DENDRO_9;
DENDRO_56 *= DENDRO_13;
DENDRO_59 = 0.5;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_59;
DENDRO_60 *= DENDRO_38;
DENDRO_60 *= DENDRO_5;
DENDRO_60 *= DENDRO_9;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_5;
DENDRO_59 *= DENDRO_15;
DENDRO_59 *= DENDRO_1;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_5;
DENDRO_61 *= DENDRO_20;
DENDRO_61 *= DENDRO_14;
DENDRO_62 = 0.5;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_62;
DENDRO_63 *= DENDRO_42;
DENDRO_63 *= DENDRO_5;
DENDRO_63 *= DENDRO_1;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_5;
DENDRO_62 *= DENDRO_17;
DENDRO_62 *= DENDRO_9;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_5;
DENDRO_64 *= DENDRO_20;
DENDRO_64 *= DENDRO_16;
DENDRO_65 = 0.5;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_65;
DENDRO_66 *= DENDRO_46;
DENDRO_66 *= DENDRO_5;
DENDRO_66 *= DENDRO_20;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_5;
DENDRO_65 *= DENDRO_18;
DENDRO_65 *= DENDRO_1;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_5;
DENDRO_67 *= DENDRO_19;
DENDRO_67 *= DENDRO_9;
DENDRO_68 = 0.5;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_68;
DENDRO_69 *= DENDRO_12;
DENDRO_69 *= DENDRO_5;
DENDRO_69 *= DENDRO_8;
DENDRO_12 = 0.5;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_12;
DENDRO_68 *= DENDRO_36;
DENDRO_68 *= DENDRO_5;
DENDRO_68 *= DENDRO_1;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_5;
DENDRO_12 *= DENDRO_0;
DENDRO_12 *= DENDRO_7;
DENDRO_7 = 0.5;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_7;
DENDRO_36 *= DENDRO_24;
DENDRO_36 *= DENDRO_5;
DENDRO_36 *= DENDRO_8;
DENDRO_7 = 0.5;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_7;
DENDRO_24 *= DENDRO_30;
DENDRO_24 *= DENDRO_5;
DENDRO_24 *= DENDRO_0;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_5;
DENDRO_7 *= DENDRO_1;
DENDRO_7 *= DENDRO_11;
DENDRO_11 = 0.5;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_11;
DENDRO_30 *= DENDRO_27;
DENDRO_30 *= DENDRO_5;
DENDRO_30 *= DENDRO_0;
DENDRO_11 = 0.5;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_11;
DENDRO_27 *= DENDRO_33;
DENDRO_27 *= DENDRO_5;
DENDRO_27 *= DENDRO_1;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_5;
DENDRO_11 *= DENDRO_8;
DENDRO_11 *= DENDRO_13;
DENDRO_13 = 0.5;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_13;
DENDRO_33 *= DENDRO_38;
DENDRO_33 *= DENDRO_5;
DENDRO_33 *= DENDRO_8;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_5;
DENDRO_13 *= DENDRO_14;
DENDRO_13 *= DENDRO_1;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_5;
DENDRO_14 *= DENDRO_0;
DENDRO_14 *= DENDRO_15;
DENDRO_15 = 0.5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_15;
DENDRO_38 *= DENDRO_42;
DENDRO_38 *= DENDRO_5;
DENDRO_38 *= DENDRO_0;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_5;
DENDRO_15 *= DENDRO_16;
DENDRO_15 *= DENDRO_1;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_5;
DENDRO_16 *= DENDRO_17;
DENDRO_16 *= DENDRO_8;
DENDRO_17 = 0.5;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_17;
DENDRO_42 *= DENDRO_46;
DENDRO_42 *= DENDRO_5;
DENDRO_42 *= DENDRO_1;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_5;
DENDRO_17 *= DENDRO_19;
DENDRO_17 *= DENDRO_8;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_5;
DENDRO_19 *= DENDRO_0;
DENDRO_19 *= DENDRO_18;
DENDRO_18 = At5[pp];

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_18;
DENDRO_46 *= DENDRO_2;
DENDRO_70 = DENDRO_igt5;
DENDRO_46 *= DENDRO_70;
DENDRO_71 = At4[pp];

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_71;
DENDRO_72 *= DENDRO_2;
DENDRO_73 = DENDRO_igt4;
DENDRO_72 *= DENDRO_73;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_71;
DENDRO_74 *= DENDRO_10;
DENDRO_74 *= DENDRO_70;
DENDRO_75 = At3[pp];

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_75;
DENDRO_76 *= DENDRO_10;
DENDRO_76 *= DENDRO_73;
DENDRO_77 = At2[pp];

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_77;
DENDRO_78 *= DENDRO_22;
DENDRO_78 *= DENDRO_70;
DENDRO_79 = At1[pp];

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_79;
DENDRO_80 *= DENDRO_10;
DENDRO_80 *= DENDRO_2;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_79;
DENDRO_81 *= DENDRO_22;
DENDRO_81 *= DENDRO_73;
DENDRO_82 = At0[pp];

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_82;
DENDRO_83 *= DENDRO_22;
DENDRO_83 *= DENDRO_2;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_77;
DENDRO_84 *= DENDRO_3;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_18;
DENDRO_85 *= DENDRO_2;
DENDRO_85 *= DENDRO_73;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_71;
DENDRO_86 *= DENDRO_2;
DENDRO_87 = DENDRO_igt3;
DENDRO_86 *= DENDRO_87;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_71;
DENDRO_88 *= DENDRO_10;
DENDRO_88 *= DENDRO_73;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_75;
DENDRO_89 *= DENDRO_10;
DENDRO_89 *= DENDRO_87;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_77;
DENDRO_90 *= DENDRO_10;
DENDRO_90 *= DENDRO_2;

// initialize reduction 
DENDRO_91 = 1;
DENDRO_91 *= DENDRO_77;
DENDRO_91 *= DENDRO_22;
DENDRO_91 *= DENDRO_73;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_79;
DENDRO_92 *= DENDRO_22;
DENDRO_92 *= DENDRO_87;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_82;
DENDRO_93 *= DENDRO_22;
DENDRO_93 *= DENDRO_10;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_79;
DENDRO_94 *= DENDRO_21;
DENDRO_95 = 2;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_95;
DENDRO_96 *= DENDRO_71;
DENDRO_96 *= DENDRO_10;
DENDRO_96 *= DENDRO_2;

// initialize reduction 
DENDRO_97 = 1;
DENDRO_97 *= DENDRO_95;
DENDRO_97 *= DENDRO_77;
DENDRO_97 *= DENDRO_22;
DENDRO_97 *= DENDRO_2;

// initialize reduction 
DENDRO_98 = 1;
DENDRO_98 *= DENDRO_95;
DENDRO_98 *= DENDRO_79;
DENDRO_98 *= DENDRO_22;
DENDRO_98 *= DENDRO_10;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_18;
DENDRO_99 *= DENDRO_3;

// initialize reduction 
DENDRO_100 = 1;
DENDRO_100 *= DENDRO_75;
DENDRO_100 *= DENDRO_21;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_82;
DENDRO_101 *= DENDRO_23;

// int-power reduce pow(DENDRO_igt4,2)
DENDRO_102 = DENDRO_73;
DENDRO_102 *= DENDRO_73;

// int-power reduce pow(DENDRO_igt5,2)
DENDRO_103 = DENDRO_70;
DENDRO_103 *= DENDRO_70;

// int-power reduce pow(DENDRO_igt3,2)
DENDRO_104 = DENDRO_87;
DENDRO_104 *= DENDRO_87;

// initialize reduction 
DENDRO_105 = 0;
DENDRO_105 += DENDRO_25;
DENDRO_105 += DENDRO_28;
DENDRO_105 += DENDRO_26;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_29;
DENDRO_25 += DENDRO_32;
DENDRO_25 += DENDRO_31;

// initialize reduction 
DENDRO_26 = 0;
DENDRO_26 += DENDRO_34;
DENDRO_26 += DENDRO_37;
DENDRO_26 += DENDRO_35;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_41;
DENDRO_28 += DENDRO_40;
DENDRO_28 += DENDRO_39;

// initialize reduction 
DENDRO_29 = 0;
DENDRO_29 += DENDRO_45;
DENDRO_29 += DENDRO_44;
DENDRO_29 += DENDRO_43;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_49;
DENDRO_31 += DENDRO_48;
DENDRO_31 += DENDRO_47;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_50;
DENDRO_32 += DENDRO_52;
DENDRO_32 += DENDRO_51;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_53;
DENDRO_34 += DENDRO_55;
DENDRO_34 += DENDRO_54;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_56;
DENDRO_35 += DENDRO_58;
DENDRO_35 += DENDRO_57;

// initialize reduction 
DENDRO_37 = 0;
DENDRO_37 += DENDRO_61;
DENDRO_37 += DENDRO_59;
DENDRO_37 += DENDRO_60;

// initialize reduction 
DENDRO_39 = 0;
DENDRO_39 += DENDRO_64;
DENDRO_39 += DENDRO_62;
DENDRO_39 += DENDRO_63;

// initialize reduction 
DENDRO_40 = 0;
DENDRO_40 += DENDRO_67;
DENDRO_40 += DENDRO_65;
DENDRO_40 += DENDRO_66;

// initialize reduction 
DENDRO_41 = 0;
DENDRO_41 += DENDRO_12;
DENDRO_41 += DENDRO_68;
DENDRO_41 += DENDRO_69;

// initialize reduction 
DENDRO_12 = 0;
DENDRO_12 += DENDRO_7;
DENDRO_12 += DENDRO_24;
DENDRO_12 += DENDRO_36;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_11;
DENDRO_7 += DENDRO_27;
DENDRO_7 += DENDRO_30;

// initialize reduction 
DENDRO_11 = 0;
DENDRO_11 += DENDRO_14;
DENDRO_11 += DENDRO_13;
DENDRO_11 += DENDRO_33;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_16;
DENDRO_13 += DENDRO_15;
DENDRO_13 += DENDRO_38;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_19;
DENDRO_14 += DENDRO_17;
DENDRO_14 += DENDRO_42;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_84;
DENDRO_15 += DENDRO_83;
DENDRO_15 += DENDRO_81;
DENDRO_15 += DENDRO_80;
DENDRO_15 += DENDRO_78;
DENDRO_15 += DENDRO_76;
DENDRO_15 += DENDRO_74;
DENDRO_15 += DENDRO_72;
DENDRO_15 += DENDRO_46;
DENDRO_16 = chi[pp];

// int-power reduce pow(chi[pp],-1)
DENDRO_17 = DENDRO_16;
DENDRO_17 = 1.0/DENDRO_17;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_94;
DENDRO_16 += DENDRO_93;
DENDRO_16 += DENDRO_92;
DENDRO_16 += DENDRO_91;
DENDRO_16 += DENDRO_90;
DENDRO_16 += DENDRO_89;
DENDRO_16 += DENDRO_88;
DENDRO_16 += DENDRO_86;
DENDRO_16 += DENDRO_85;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_101;
DENDRO_19 += DENDRO_100;
DENDRO_19 += DENDRO_99;
DENDRO_19 += DENDRO_98;
DENDRO_19 += DENDRO_97;
DENDRO_19 += DENDRO_96;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_18;
DENDRO_24 *= DENDRO_73;
DENDRO_24 *= DENDRO_70;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_71;
DENDRO_27 *= DENDRO_87;
DENDRO_27 *= DENDRO_70;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_75;
DENDRO_30 *= DENDRO_87;
DENDRO_30 *= DENDRO_73;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_77;
DENDRO_33 *= DENDRO_2;
DENDRO_33 *= DENDRO_73;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_77;
DENDRO_36 *= DENDRO_10;
DENDRO_36 *= DENDRO_70;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_79;
DENDRO_38 *= DENDRO_2;
DENDRO_38 *= DENDRO_87;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_79;
DENDRO_42 *= DENDRO_10;
DENDRO_42 *= DENDRO_73;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_82;
DENDRO_43 *= DENDRO_10;
DENDRO_43 *= DENDRO_2;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_71;
DENDRO_44 *= DENDRO_102;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_95;
DENDRO_45 *= DENDRO_71;
DENDRO_45 *= DENDRO_73;
DENDRO_45 *= DENDRO_70;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_95;
DENDRO_46 *= DENDRO_77;
DENDRO_46 *= DENDRO_2;
DENDRO_46 *= DENDRO_70;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_95;
DENDRO_47 *= DENDRO_79;
DENDRO_47 *= DENDRO_2;
DENDRO_47 *= DENDRO_73;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_18;
DENDRO_48 *= DENDRO_103;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_75;
DENDRO_49 *= DENDRO_102;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_82;
DENDRO_50 *= DENDRO_3;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_95;
DENDRO_51 *= DENDRO_71;
DENDRO_51 *= DENDRO_87;
DENDRO_51 *= DENDRO_73;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_95;
DENDRO_52 *= DENDRO_77;
DENDRO_52 *= DENDRO_10;
DENDRO_52 *= DENDRO_73;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_95;
DENDRO_53 *= DENDRO_79;
DENDRO_53 *= DENDRO_10;
DENDRO_53 *= DENDRO_87;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_18;
DENDRO_54 *= DENDRO_102;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_75;
DENDRO_55 *= DENDRO_104;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_82;
DENDRO_56 *= DENDRO_21;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_95;
DENDRO_57 *= DENDRO_5;
DENDRO_57 *= DENDRO_9;
DENDRO_57 *= DENDRO_105;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_95;
DENDRO_58 *= DENDRO_5;
DENDRO_58 *= DENDRO_8;
DENDRO_58 *= DENDRO_25;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_95;
DENDRO_59 *= DENDRO_5;
DENDRO_59 *= DENDRO_1;
DENDRO_59 *= DENDRO_26;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_5;
DENDRO_60 *= DENDRO_4;
DENDRO_60 *= DENDRO_28;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_5;
DENDRO_61 *= DENDRO_0;
DENDRO_61 *= DENDRO_29;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_5;
DENDRO_62 *= DENDRO_20;
DENDRO_62 *= DENDRO_31;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_95;
DENDRO_63 *= DENDRO_5;
DENDRO_63 *= DENDRO_9;
DENDRO_63 *= DENDRO_32;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_95;
DENDRO_32 *= DENDRO_5;
DENDRO_32 *= DENDRO_8;
DENDRO_32 *= DENDRO_34;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_95;
DENDRO_34 *= DENDRO_5;
DENDRO_34 *= DENDRO_1;
DENDRO_34 *= DENDRO_35;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_5;
DENDRO_35 *= DENDRO_4;
DENDRO_35 *= DENDRO_37;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_5;
DENDRO_37 *= DENDRO_0;
DENDRO_37 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_5;
DENDRO_39 *= DENDRO_20;
DENDRO_39 *= DENDRO_40;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_95;
DENDRO_40 *= DENDRO_18;
DENDRO_40 *= DENDRO_2;
DENDRO_40 *= DENDRO_70;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_95;
DENDRO_64 *= DENDRO_71;
DENDRO_64 *= DENDRO_2;
DENDRO_64 *= DENDRO_73;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_95;
DENDRO_65 *= DENDRO_71;
DENDRO_65 *= DENDRO_10;
DENDRO_65 *= DENDRO_70;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_95;
DENDRO_66 *= DENDRO_75;
DENDRO_66 *= DENDRO_10;
DENDRO_66 *= DENDRO_73;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_95;
DENDRO_67 *= DENDRO_77;
DENDRO_67 *= DENDRO_22;
DENDRO_67 *= DENDRO_70;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_95;
DENDRO_68 *= DENDRO_79;
DENDRO_68 *= DENDRO_10;
DENDRO_68 *= DENDRO_2;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_95;
DENDRO_69 *= DENDRO_79;
DENDRO_69 *= DENDRO_22;
DENDRO_69 *= DENDRO_73;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_95;
DENDRO_72 *= DENDRO_82;
DENDRO_72 *= DENDRO_22;
DENDRO_72 *= DENDRO_2;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_95;
DENDRO_74 *= DENDRO_77;
DENDRO_74 *= DENDRO_3;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_95;
DENDRO_76 *= DENDRO_5;
DENDRO_76 *= DENDRO_9;
DENDRO_76 *= DENDRO_41;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_95;
DENDRO_9 *= DENDRO_5;
DENDRO_9 *= DENDRO_8;
DENDRO_9 *= DENDRO_12;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_95;
DENDRO_8 *= DENDRO_5;
DENDRO_8 *= DENDRO_1;
DENDRO_8 *= DENDRO_7;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_5;
DENDRO_1 *= DENDRO_4;
DENDRO_1 *= DENDRO_11;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_5;
DENDRO_4 *= DENDRO_0;
DENDRO_4 *= DENDRO_13;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_5;
DENDRO_0 *= DENDRO_20;
DENDRO_0 *= DENDRO_14;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_95;
DENDRO_5 *= DENDRO_18;
DENDRO_5 *= DENDRO_2;
DENDRO_5 *= DENDRO_73;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_95;
DENDRO_7 *= DENDRO_71;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_87;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_95;
DENDRO_11 *= DENDRO_71;
DENDRO_11 *= DENDRO_10;
DENDRO_11 *= DENDRO_73;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_95;
DENDRO_12 *= DENDRO_75;
DENDRO_12 *= DENDRO_10;
DENDRO_12 *= DENDRO_87;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_95;
DENDRO_13 *= DENDRO_77;
DENDRO_13 *= DENDRO_10;
DENDRO_13 *= DENDRO_2;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_95;
DENDRO_14 *= DENDRO_77;
DENDRO_14 *= DENDRO_22;
DENDRO_14 *= DENDRO_73;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_95;
DENDRO_20 *= DENDRO_79;
DENDRO_20 *= DENDRO_22;
DENDRO_20 *= DENDRO_87;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_95;
DENDRO_41 *= DENDRO_82;
DENDRO_41 *= DENDRO_22;
DENDRO_41 *= DENDRO_10;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_95;
DENDRO_78 *= DENDRO_79;
DENDRO_78 *= DENDRO_21;
DENDRO_80 = 4;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_80;
DENDRO_81 *= DENDRO_71;
DENDRO_81 *= DENDRO_10;
DENDRO_81 *= DENDRO_2;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_80;
DENDRO_71 *= DENDRO_77;
DENDRO_71 *= DENDRO_22;
DENDRO_71 *= DENDRO_2;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_80;
DENDRO_77 *= DENDRO_79;
DENDRO_77 *= DENDRO_22;
DENDRO_77 *= DENDRO_10;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_95;
DENDRO_79 *= DENDRO_18;
DENDRO_79 *= DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_95;
DENDRO_3 *= DENDRO_75;
DENDRO_3 *= DENDRO_21;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_95;
DENDRO_18 *= DENDRO_82;
DENDRO_18 *= DENDRO_23;
DENDRO_21 = 3;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_21;
DENDRO_75 = grad_2_chi;
DENDRO_23 *= DENDRO_75;
DENDRO_23 *= DENDRO_17;
DENDRO_23 *= DENDRO_15;
DENDRO_75 = 4.0/3.0;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_75;
DENDRO_82 *= DENDRO_2;
DENDRO_83 = grad_2_K;
DENDRO_82 *= DENDRO_83;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_21;
DENDRO_84 = grad_1_chi;
DENDRO_83 *= DENDRO_84;
DENDRO_83 *= DENDRO_17;
DENDRO_83 *= DENDRO_16;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_75;
DENDRO_84 *= DENDRO_10;
DENDRO_85 = grad_1_K;
DENDRO_84 *= DENDRO_85;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_21;
DENDRO_21 = grad_0_chi;
DENDRO_85 *= DENDRO_21;
DENDRO_85 *= DENDRO_17;
DENDRO_85 *= DENDRO_19;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_75;
DENDRO_17 *= DENDRO_22;
DENDRO_21 = grad_0_K;
DENDRO_17 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 0;
DENDRO_21 += DENDRO_44;
DENDRO_21 += DENDRO_43;
DENDRO_21 += DENDRO_42;
DENDRO_21 += DENDRO_38;
DENDRO_21 += DENDRO_36;
DENDRO_21 += DENDRO_33;
DENDRO_21 += DENDRO_30;
DENDRO_21 += DENDRO_27;
DENDRO_21 += DENDRO_24;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_50;
DENDRO_24 += DENDRO_49;
DENDRO_24 += DENDRO_48;
DENDRO_24 += DENDRO_47;
DENDRO_24 += DENDRO_46;
DENDRO_24 += DENDRO_45;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_56;
DENDRO_27 += DENDRO_55;
DENDRO_27 += DENDRO_54;
DENDRO_27 += DENDRO_53;
DENDRO_27 += DENDRO_52;
DENDRO_27 += DENDRO_51;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_62;
DENDRO_30 += DENDRO_61;
DENDRO_30 += DENDRO_60;
DENDRO_30 += DENDRO_59;
DENDRO_30 += DENDRO_58;
DENDRO_30 += DENDRO_57;
DENDRO_33 = grad_0_beta0;

// initialize reduction 
DENDRO_36 = 0;
DENDRO_36 += DENDRO_33;
DENDRO_38 = grad_1_beta1;
DENDRO_36 += DENDRO_38;
DENDRO_38 = grad_2_beta2;
DENDRO_36 += DENDRO_38;

// initialize reduction 
DENDRO_38 = 0;
DENDRO_38 += DENDRO_39;
DENDRO_38 += DENDRO_37;
DENDRO_38 += DENDRO_35;
DENDRO_38 += DENDRO_34;
DENDRO_38 += DENDRO_32;
DENDRO_38 += DENDRO_63;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_74;
DENDRO_32 += DENDRO_72;
DENDRO_32 += DENDRO_69;
DENDRO_32 += DENDRO_68;
DENDRO_32 += DENDRO_67;
DENDRO_32 += DENDRO_66;
DENDRO_32 += DENDRO_65;
DENDRO_32 += DENDRO_64;
DENDRO_32 += DENDRO_40;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_0;
DENDRO_34 += DENDRO_4;
DENDRO_34 += DENDRO_1;
DENDRO_34 += DENDRO_8;
DENDRO_34 += DENDRO_9;
DENDRO_34 += DENDRO_76;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_78;
DENDRO_0 += DENDRO_41;
DENDRO_0 += DENDRO_20;
DENDRO_0 += DENDRO_14;
DENDRO_0 += DENDRO_13;
DENDRO_0 += DENDRO_12;
DENDRO_0 += DENDRO_11;
DENDRO_0 += DENDRO_7;
DENDRO_0 += DENDRO_5;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_18;
DENDRO_1 += DENDRO_3;
DENDRO_1 += DENDRO_79;
DENDRO_1 += DENDRO_77;
DENDRO_1 += DENDRO_71;
DENDRO_1 += DENDRO_81;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_82;
DENDRO_3 += DENDRO_23;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_84;
DENDRO_4 += DENDRO_83;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_17;
DENDRO_5 += DENDRO_85;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_80;
DENDRO_8 = alpha[pp];
DENDRO_7 *= DENDRO_8;
DENDRO_7 *= DENDRO_25;
DENDRO_7 *= DENDRO_16;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_80;
DENDRO_9 *= DENDRO_8;
DENDRO_9 *= DENDRO_105;
DENDRO_9 *= DENDRO_15;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_80;
DENDRO_11 *= DENDRO_8;
DENDRO_11 *= DENDRO_26;
DENDRO_11 *= DENDRO_21;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_95;
DENDRO_12 *= DENDRO_8;
DENDRO_12 *= DENDRO_28;
DENDRO_12 *= DENDRO_19;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_95;
DENDRO_13 *= DENDRO_8;
DENDRO_13 *= DENDRO_31;
DENDRO_13 *= DENDRO_24;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_95;
DENDRO_14 *= DENDRO_8;
DENDRO_14 *= DENDRO_29;
DENDRO_14 *= DENDRO_27;
DENDRO_15 = 7.0/3.0;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_15;
DENDRO_16 *= DENDRO_2;
DENDRO_17 = grad2_0_2_beta0;
DENDRO_16 *= DENDRO_17;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_15;
DENDRO_17 *= DENDRO_10;
DENDRO_15 = grad2_0_1_beta0;
DENDRO_17 *= DENDRO_15;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_75;
DENDRO_15 *= DENDRO_22;
DENDRO_18 = grad2_0_0_beta0;
DENDRO_15 *= DENDRO_18;
DENDRO_18 = 2.0/3.0;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_18;
DENDRO_19 *= DENDRO_36;
DENDRO_19 *= DENDRO_30;
DENDRO_18 = 1.0/3.0;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_18;
DENDRO_20 *= DENDRO_2;
DENDRO_21 = grad2_2_2_beta2;
DENDRO_20 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_18;
DENDRO_21 *= DENDRO_2;
DENDRO_2 = grad2_1_2_beta1;
DENDRO_21 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_18;
DENDRO_2 *= DENDRO_10;
DENDRO_23 = grad2_1_2_beta2;
DENDRO_2 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_18;
DENDRO_23 *= DENDRO_10;
DENDRO_10 = grad2_1_1_beta1;
DENDRO_23 *= DENDRO_10;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_18;
DENDRO_10 *= DENDRO_22;
DENDRO_24 = grad2_0_2_beta2;
DENDRO_10 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_18;
DENDRO_24 *= DENDRO_22;
DENDRO_18 = grad2_0_1_beta1;
DENDRO_24 *= DENDRO_18;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_95;
DENDRO_18 *= DENDRO_73;
DENDRO_22 = grad2_1_2_beta0;
DENDRO_18 *= DENDRO_22;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_6;
DENDRO_25 = grad_2_beta0;
DENDRO_22 *= DENDRO_25;
DENDRO_22 *= DENDRO_38;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_6;
DENDRO_26 = grad_2_alpha;
DENDRO_25 *= DENDRO_26;
DENDRO_25 *= DENDRO_32;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_6;
DENDRO_27 = grad_1_beta0;
DENDRO_26 *= DENDRO_27;
DENDRO_26 *= DENDRO_34;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_6;
DENDRO_28 = grad_1_alpha;
DENDRO_27 *= DENDRO_28;
DENDRO_27 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_6;
DENDRO_0 *= DENDRO_33;
DENDRO_0 *= DENDRO_30;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_6;
DENDRO_29 = grad_0_alpha;
DENDRO_28 *= DENDRO_29;
DENDRO_28 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_6;
DENDRO_1 *= DENDRO_8;
DENDRO_1 *= DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_6;
DENDRO_3 *= DENDRO_8;
DENDRO_3 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_6;
DENDRO_4 *= DENDRO_8;
DENDRO_4 *= DENDRO_5;
DENDRO_5 = beta2[pp];

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_5;
DENDRO_5 = grad_2_Gt0;
DENDRO_6 *= DENDRO_5;
DENDRO_5 = beta1[pp];

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_5;
DENDRO_5 = grad_1_Gt0;
DENDRO_8 *= DENDRO_5;
DENDRO_5 = beta0[pp];

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_5;
DENDRO_5 = grad_0_Gt0;
DENDRO_29 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_70;
DENDRO_30 = grad2_2_2_beta0;
DENDRO_5 *= DENDRO_30;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_87;
DENDRO_31 = grad2_1_1_beta0;
DENDRO_30 *= DENDRO_31;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_30;
DENDRO_31 += DENDRO_5;
DENDRO_31 += DENDRO_29;
DENDRO_31 += DENDRO_8;
DENDRO_31 += DENDRO_6;
DENDRO_31 += DENDRO_4;
DENDRO_31 += DENDRO_3;
DENDRO_31 += DENDRO_1;
DENDRO_31 += DENDRO_28;
DENDRO_31 += DENDRO_0;
DENDRO_31 += DENDRO_27;
DENDRO_31 += DENDRO_26;
DENDRO_31 += DENDRO_25;
DENDRO_31 += DENDRO_22;
DENDRO_31 += DENDRO_18;
DENDRO_31 += DENDRO_24;
DENDRO_31 += DENDRO_10;
DENDRO_31 += DENDRO_23;
DENDRO_31 += DENDRO_2;
DENDRO_31 += DENDRO_21;
DENDRO_31 += DENDRO_20;
DENDRO_31 += DENDRO_19;
DENDRO_31 += DENDRO_15;
DENDRO_31 += DENDRO_17;
DENDRO_31 += DENDRO_16;
DENDRO_31 += DENDRO_14;
DENDRO_31 += DENDRO_13;
DENDRO_31 += DENDRO_12;
DENDRO_31 += DENDRO_11;
DENDRO_31 += DENDRO_9;
DENDRO_31 += DENDRO_7;
Gt_rhs0[pp] = DENDRO_31;


}
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , Gt1, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_Gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_Gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_Gt1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_Gt1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_Gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_Gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_Gt1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_Gt1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_Gt1=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , Gt1, blk);
device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_Gt1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_Gt1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_Gt1=Du[gidx];
__syncthreads();

//|V|= 396 |E| =905
//topological generations = 10
// [0]
// traversal id=0 temp_var_count=106
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;double DENDRO_33;double DENDRO_34;double DENDRO_35;double DENDRO_36;double DENDRO_37;double DENDRO_38;double DENDRO_39;double DENDRO_40;double DENDRO_41;double DENDRO_42;double DENDRO_43;double DENDRO_44;double DENDRO_45;double DENDRO_46;double DENDRO_47;double DENDRO_48;double DENDRO_49;double DENDRO_50;double DENDRO_51;double DENDRO_52;double DENDRO_53;double DENDRO_54;double DENDRO_55;double DENDRO_56;double DENDRO_57;double DENDRO_58;double DENDRO_59;double DENDRO_60;double DENDRO_61;double DENDRO_62;double DENDRO_63;double DENDRO_64;double DENDRO_65;double DENDRO_66;double DENDRO_67;double DENDRO_68;double DENDRO_69;double DENDRO_70;double DENDRO_71;double DENDRO_72;double DENDRO_73;double DENDRO_74;double DENDRO_75;double DENDRO_76;double DENDRO_77;double DENDRO_78;double DENDRO_79;double DENDRO_80;double DENDRO_81;double DENDRO_82;double DENDRO_83;double DENDRO_84;double DENDRO_85;double DENDRO_86;double DENDRO_87;double DENDRO_88;double DENDRO_89;double DENDRO_90;double DENDRO_91;double DENDRO_92;double DENDRO_93;double DENDRO_94;double DENDRO_95;double DENDRO_96;double DENDRO_97;double DENDRO_98;double DENDRO_99;double DENDRO_100;double DENDRO_101;double DENDRO_102;double DENDRO_103;double DENDRO_104;double DENDRO_105;
DENDRO_0 = gt1[pp];

// int-power reduce pow(gt1[pp],2)
DENDRO_1 = DENDRO_0;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = gt2[pp];

// int-power reduce pow(gt2[pp],2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = gt4[pp];

// int-power reduce pow(gt4[pp],2)
DENDRO_5 = DENDRO_4;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = -2;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_4;
DENDRO_6 = -1;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_6;
DENDRO_9 = gt0[pp];
DENDRO_8 *= DENDRO_9;
DENDRO_10 = gt3[pp];
DENDRO_8 *= DENDRO_10;
DENDRO_11 = gt5[pp];
DENDRO_8 *= DENDRO_11;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_11;
DENDRO_12 *= DENDRO_1;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_10;
DENDRO_13 *= DENDRO_3;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_9;
DENDRO_14 *= DENDRO_5;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_6;
DENDRO_15 *= DENDRO_2;
DENDRO_15 *= DENDRO_4;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_0;
DENDRO_16 *= DENDRO_11;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_14;
DENDRO_17 += DENDRO_13;
DENDRO_17 += DENDRO_12;
DENDRO_17 += DENDRO_8;
DENDRO_17 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_9;
DENDRO_8 *= DENDRO_4;
DENDRO_12 = -0.5;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_12;
DENDRO_12 = grad_1_gt2;
DENDRO_13 *= DENDRO_12;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_14;
DENDRO_14 = grad_2_gt1;
DENDRO_18 *= DENDRO_14;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_19;
DENDRO_19 = grad_0_gt4;
DENDRO_20 *= DENDRO_19;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_6;
DENDRO_21 *= DENDRO_9;
DENDRO_21 *= DENDRO_11;
DENDRO_22 = -0.5;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_22;
DENDRO_23 *= DENDRO_14;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_14;
DENDRO_22 *= DENDRO_12;
DENDRO_12 = -0.5;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_12;
DENDRO_14 *= DENDRO_19;
DENDRO_12 = -0.5;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_12;
DENDRO_12 = grad_2_gt0;
DENDRO_19 *= DENDRO_12;
DENDRO_24 = 1.0;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_24;
DENDRO_24 = grad_0_gt2;
DENDRO_25 *= DENDRO_24;
DENDRO_24 = -0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_24;
DENDRO_24 = grad_1_gt0;
DENDRO_26 *= DENDRO_24;
DENDRO_27 = 1.0;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_27;
DENDRO_27 = grad_0_gt1;
DENDRO_28 *= DENDRO_27;
DENDRO_27 = -0.5;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_27;
DENDRO_27 = grad_2_gt3;
DENDRO_29 *= DENDRO_27;
DENDRO_30 = 1.0;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_30;
DENDRO_30 = grad_1_gt4;
DENDRO_31 *= DENDRO_30;
DENDRO_30 = -0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_30;
DENDRO_30 = grad_0_gt3;
DENDRO_32 *= DENDRO_30;
DENDRO_33 = 1.0;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_33;
DENDRO_33 = grad_1_gt1;
DENDRO_34 *= DENDRO_33;
DENDRO_33 = -0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_33;
DENDRO_33 = grad_0_gt5;
DENDRO_35 *= DENDRO_33;
DENDRO_36 = 1.0;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = grad_2_gt2;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = -0.5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_36;
DENDRO_36 = grad_1_gt5;
DENDRO_38 *= DENDRO_36;
DENDRO_39 = 1.0;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_39;
DENDRO_39 = grad_2_gt4;
DENDRO_40 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_6;
DENDRO_39 *= DENDRO_0;
DENDRO_39 *= DENDRO_4;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_2;
DENDRO_0 *= DENDRO_10;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_6;
DENDRO_2 *= DENDRO_9;
DENDRO_2 *= DENDRO_10;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_6;
DENDRO_4 *= DENDRO_10;
DENDRO_4 *= DENDRO_11;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_16;
DENDRO_9 += DENDRO_15;

// int-power reduce pow(-gt0[pp]*gt3[pp]*gt5[pp] + gt0[pp]*gt4[pp]**2 + gt1[pp]**2*gt5[pp] - 2*gt1[pp]*gt2[pp]*gt4[pp] + gt2[pp]**2*gt3[pp],-1)
DENDRO_10 = DENDRO_17;
DENDRO_10 = 1.0/DENDRO_10;

// initialize reduction 
DENDRO_11 = 0;
DENDRO_11 += DENDRO_8;
DENDRO_11 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_20;
DENDRO_7 += DENDRO_18;
DENDRO_7 += DENDRO_13;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_3;
DENDRO_8 += DENDRO_21;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_20;
DENDRO_3 += DENDRO_22;
DENDRO_3 += DENDRO_23;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_22;
DENDRO_13 += DENDRO_18;
DENDRO_13 += DENDRO_14;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_25;
DENDRO_14 += DENDRO_19;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_28;
DENDRO_15 += DENDRO_26;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_31;
DENDRO_16 += DENDRO_29;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_34;
DENDRO_17 += DENDRO_32;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_37;
DENDRO_18 += DENDRO_35;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_40;
DENDRO_19 += DENDRO_38;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_0;
DENDRO_20 += DENDRO_39;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_1;
DENDRO_0 += DENDRO_2;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_5;
DENDRO_1 += DENDRO_4;
DENDRO_2 = DENDRO_igt4;

// int-power reduce pow(DENDRO_igt4,2)
DENDRO_4 = DENDRO_2;
DENDRO_4 *= DENDRO_2;
DENDRO_5 = DENDRO_igt3;

// int-power reduce pow(DENDRO_igt3,2)
DENDRO_21 = DENDRO_5;
DENDRO_21 *= DENDRO_5;
DENDRO_22 = DENDRO_igt1;

// int-power reduce pow(DENDRO_igt1,2)
DENDRO_23 = DENDRO_22;
DENDRO_23 *= DENDRO_22;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_25;
DENDRO_26 *= DENDRO_12;
DENDRO_26 *= DENDRO_10;
DENDRO_26 *= DENDRO_9;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_25;
DENDRO_28 *= DENDRO_33;
DENDRO_28 *= DENDRO_10;
DENDRO_28 *= DENDRO_11;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_10;
DENDRO_25 *= DENDRO_8;
DENDRO_25 *= DENDRO_7;
DENDRO_29 = 0.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_29;
DENDRO_31 *= DENDRO_24;
DENDRO_31 *= DENDRO_10;
DENDRO_31 *= DENDRO_9;
DENDRO_29 = 0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_29;
DENDRO_32 *= DENDRO_30;
DENDRO_32 *= DENDRO_10;
DENDRO_32 *= DENDRO_8;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_10;
DENDRO_29 *= DENDRO_11;
DENDRO_29 *= DENDRO_3;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_34;
DENDRO_35 *= DENDRO_27;
DENDRO_35 *= DENDRO_10;
DENDRO_35 *= DENDRO_8;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_34;
DENDRO_37 *= DENDRO_36;
DENDRO_37 *= DENDRO_10;
DENDRO_37 *= DENDRO_11;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_10;
DENDRO_34 *= DENDRO_9;
DENDRO_34 *= DENDRO_13;
DENDRO_38 = 0.5;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_38;
DENDRO_38 = grad_0_gt0;
DENDRO_39 *= DENDRO_38;
DENDRO_39 *= DENDRO_10;
DENDRO_39 *= DENDRO_9;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_10;
DENDRO_40 *= DENDRO_14;
DENDRO_40 *= DENDRO_11;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_10;
DENDRO_41 *= DENDRO_8;
DENDRO_41 *= DENDRO_15;
DENDRO_42 = 0.5;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_42;
DENDRO_42 = grad_1_gt3;
DENDRO_43 *= DENDRO_42;
DENDRO_43 *= DENDRO_10;
DENDRO_43 *= DENDRO_8;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_10;
DENDRO_44 *= DENDRO_16;
DENDRO_44 *= DENDRO_11;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_10;
DENDRO_45 *= DENDRO_17;
DENDRO_45 *= DENDRO_9;
DENDRO_46 = 0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_46;
DENDRO_46 = grad_2_gt5;
DENDRO_47 *= DENDRO_46;
DENDRO_47 *= DENDRO_10;
DENDRO_47 *= DENDRO_11;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_10;
DENDRO_48 *= DENDRO_18;
DENDRO_48 *= DENDRO_9;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_10;
DENDRO_49 *= DENDRO_8;
DENDRO_49 *= DENDRO_19;
DENDRO_50 = 0.5;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_50;
DENDRO_51 *= DENDRO_12;
DENDRO_51 *= DENDRO_10;
DENDRO_51 *= DENDRO_20;
DENDRO_50 = 0.5;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_50;
DENDRO_52 *= DENDRO_33;
DENDRO_52 *= DENDRO_10;
DENDRO_52 *= DENDRO_0;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_10;
DENDRO_50 *= DENDRO_11;
DENDRO_50 *= DENDRO_7;
DENDRO_53 = 0.5;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_53;
DENDRO_54 *= DENDRO_24;
DENDRO_54 *= DENDRO_10;
DENDRO_54 *= DENDRO_20;
DENDRO_53 = 0.5;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_53;
DENDRO_55 *= DENDRO_30;
DENDRO_55 *= DENDRO_10;
DENDRO_55 *= DENDRO_11;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_10;
DENDRO_53 *= DENDRO_0;
DENDRO_53 *= DENDRO_3;
DENDRO_56 = 0.5;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_56;
DENDRO_57 *= DENDRO_27;
DENDRO_57 *= DENDRO_10;
DENDRO_57 *= DENDRO_11;
DENDRO_56 = 0.5;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_56;
DENDRO_58 *= DENDRO_36;
DENDRO_58 *= DENDRO_10;
DENDRO_58 *= DENDRO_0;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_10;
DENDRO_56 *= DENDRO_20;
DENDRO_56 *= DENDRO_13;
DENDRO_59 = 0.5;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_59;
DENDRO_60 *= DENDRO_38;
DENDRO_60 *= DENDRO_10;
DENDRO_60 *= DENDRO_20;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_10;
DENDRO_59 *= DENDRO_15;
DENDRO_59 *= DENDRO_11;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_10;
DENDRO_61 *= DENDRO_0;
DENDRO_61 *= DENDRO_14;
DENDRO_62 = 0.5;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_62;
DENDRO_63 *= DENDRO_42;
DENDRO_63 *= DENDRO_10;
DENDRO_63 *= DENDRO_11;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_10;
DENDRO_62 *= DENDRO_17;
DENDRO_62 *= DENDRO_20;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_10;
DENDRO_64 *= DENDRO_0;
DENDRO_64 *= DENDRO_16;
DENDRO_65 = 0.5;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_65;
DENDRO_66 *= DENDRO_46;
DENDRO_66 *= DENDRO_10;
DENDRO_66 *= DENDRO_0;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_10;
DENDRO_65 *= DENDRO_19;
DENDRO_65 *= DENDRO_11;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_10;
DENDRO_67 *= DENDRO_18;
DENDRO_67 *= DENDRO_20;
DENDRO_68 = 0.5;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_68;
DENDRO_69 *= DENDRO_12;
DENDRO_69 *= DENDRO_10;
DENDRO_69 *= DENDRO_1;
DENDRO_12 = 0.5;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_12;
DENDRO_68 *= DENDRO_33;
DENDRO_68 *= DENDRO_10;
DENDRO_68 *= DENDRO_20;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_10;
DENDRO_12 *= DENDRO_9;
DENDRO_12 *= DENDRO_7;
DENDRO_7 = 0.5;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_7;
DENDRO_33 *= DENDRO_24;
DENDRO_33 *= DENDRO_10;
DENDRO_33 *= DENDRO_1;
DENDRO_7 = 0.5;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_7;
DENDRO_24 *= DENDRO_30;
DENDRO_24 *= DENDRO_10;
DENDRO_24 *= DENDRO_9;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_10;
DENDRO_7 *= DENDRO_20;
DENDRO_7 *= DENDRO_3;
DENDRO_3 = 0.5;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_3;
DENDRO_30 *= DENDRO_27;
DENDRO_30 *= DENDRO_10;
DENDRO_30 *= DENDRO_9;
DENDRO_3 = 0.5;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_3;
DENDRO_27 *= DENDRO_36;
DENDRO_27 *= DENDRO_10;
DENDRO_27 *= DENDRO_20;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_10;
DENDRO_3 *= DENDRO_1;
DENDRO_3 *= DENDRO_13;
DENDRO_13 = 0.5;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_13;
DENDRO_36 *= DENDRO_38;
DENDRO_36 *= DENDRO_10;
DENDRO_36 *= DENDRO_1;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_10;
DENDRO_13 *= DENDRO_14;
DENDRO_13 *= DENDRO_20;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_10;
DENDRO_14 *= DENDRO_15;
DENDRO_14 *= DENDRO_9;
DENDRO_15 = 0.5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_15;
DENDRO_38 *= DENDRO_42;
DENDRO_38 *= DENDRO_10;
DENDRO_38 *= DENDRO_9;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_10;
DENDRO_15 *= DENDRO_16;
DENDRO_15 *= DENDRO_20;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_10;
DENDRO_16 *= DENDRO_1;
DENDRO_16 *= DENDRO_17;
DENDRO_17 = 0.5;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_17;
DENDRO_42 *= DENDRO_46;
DENDRO_42 *= DENDRO_10;
DENDRO_42 *= DENDRO_20;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_10;
DENDRO_17 *= DENDRO_19;
DENDRO_17 *= DENDRO_9;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_10;
DENDRO_19 *= DENDRO_1;
DENDRO_19 *= DENDRO_18;
DENDRO_18 = At5[pp];

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_18;
DENDRO_46 *= DENDRO_2;
DENDRO_70 = DENDRO_igt5;
DENDRO_46 *= DENDRO_70;
DENDRO_71 = At4[pp];

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_71;
DENDRO_72 *= DENDRO_5;
DENDRO_72 *= DENDRO_70;
DENDRO_73 = At3[pp];

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_73;
DENDRO_74 *= DENDRO_5;
DENDRO_74 *= DENDRO_2;
DENDRO_75 = At2[pp];

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_75;
DENDRO_77 = DENDRO_igt2;
DENDRO_76 *= DENDRO_77;
DENDRO_76 *= DENDRO_2;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_75;
DENDRO_78 *= DENDRO_22;
DENDRO_78 *= DENDRO_70;
DENDRO_79 = At1[pp];

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_79;
DENDRO_80 *= DENDRO_77;
DENDRO_80 *= DENDRO_5;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_79;
DENDRO_81 *= DENDRO_22;
DENDRO_81 *= DENDRO_2;
DENDRO_82 = At0[pp];

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_82;
DENDRO_83 *= DENDRO_22;
DENDRO_83 *= DENDRO_77;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_71;
DENDRO_84 *= DENDRO_4;
DENDRO_85 = 2;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_85;
DENDRO_86 *= DENDRO_71;
DENDRO_86 *= DENDRO_5;
DENDRO_86 *= DENDRO_2;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_85;
DENDRO_87 *= DENDRO_75;
DENDRO_87 *= DENDRO_22;
DENDRO_87 *= DENDRO_2;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_85;
DENDRO_88 *= DENDRO_79;
DENDRO_88 *= DENDRO_22;
DENDRO_88 *= DENDRO_5;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_18;
DENDRO_89 *= DENDRO_4;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_73;
DENDRO_90 *= DENDRO_21;

// initialize reduction 
DENDRO_91 = 1;
DENDRO_91 *= DENDRO_82;
DENDRO_91 *= DENDRO_23;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_18;
DENDRO_92 *= DENDRO_77;
DENDRO_92 *= DENDRO_2;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_71;
DENDRO_93 *= DENDRO_77;
DENDRO_93 *= DENDRO_5;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_71;
DENDRO_94 *= DENDRO_22;
DENDRO_94 *= DENDRO_2;

// initialize reduction 
DENDRO_95 = 1;
DENDRO_95 *= DENDRO_73;
DENDRO_95 *= DENDRO_22;
DENDRO_95 *= DENDRO_5;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_75;
DENDRO_96 *= DENDRO_22;
DENDRO_96 *= DENDRO_77;

// initialize reduction 
DENDRO_97 = 1;
DENDRO_97 *= DENDRO_75;
DENDRO_98 = DENDRO_igt0;
DENDRO_97 *= DENDRO_98;
DENDRO_97 *= DENDRO_2;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_79;
DENDRO_99 *= DENDRO_98;
DENDRO_99 *= DENDRO_5;

// initialize reduction 
DENDRO_100 = 1;
DENDRO_100 *= DENDRO_82;
DENDRO_100 *= DENDRO_98;
DENDRO_100 *= DENDRO_22;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_79;
DENDRO_101 *= DENDRO_23;

// int-power reduce pow(DENDRO_igt2,2)
DENDRO_102 = DENDRO_77;
DENDRO_102 *= DENDRO_77;

// int-power reduce pow(DENDRO_igt5,2)
DENDRO_103 = DENDRO_70;
DENDRO_103 *= DENDRO_70;

// int-power reduce pow(DENDRO_igt0,2)
DENDRO_104 = DENDRO_98;
DENDRO_104 *= DENDRO_98;

// initialize reduction 
DENDRO_105 = 0;
DENDRO_105 += DENDRO_25;
DENDRO_105 += DENDRO_28;
DENDRO_105 += DENDRO_26;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_29;
DENDRO_25 += DENDRO_32;
DENDRO_25 += DENDRO_31;

// initialize reduction 
DENDRO_26 = 0;
DENDRO_26 += DENDRO_34;
DENDRO_26 += DENDRO_37;
DENDRO_26 += DENDRO_35;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_41;
DENDRO_28 += DENDRO_40;
DENDRO_28 += DENDRO_39;

// initialize reduction 
DENDRO_29 = 0;
DENDRO_29 += DENDRO_45;
DENDRO_29 += DENDRO_44;
DENDRO_29 += DENDRO_43;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_49;
DENDRO_31 += DENDRO_48;
DENDRO_31 += DENDRO_47;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_50;
DENDRO_32 += DENDRO_52;
DENDRO_32 += DENDRO_51;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_53;
DENDRO_34 += DENDRO_55;
DENDRO_34 += DENDRO_54;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_56;
DENDRO_35 += DENDRO_58;
DENDRO_35 += DENDRO_57;

// initialize reduction 
DENDRO_37 = 0;
DENDRO_37 += DENDRO_61;
DENDRO_37 += DENDRO_59;
DENDRO_37 += DENDRO_60;

// initialize reduction 
DENDRO_39 = 0;
DENDRO_39 += DENDRO_64;
DENDRO_39 += DENDRO_62;
DENDRO_39 += DENDRO_63;

// initialize reduction 
DENDRO_40 = 0;
DENDRO_40 += DENDRO_67;
DENDRO_40 += DENDRO_65;
DENDRO_40 += DENDRO_66;

// initialize reduction 
DENDRO_41 = 0;
DENDRO_41 += DENDRO_12;
DENDRO_41 += DENDRO_68;
DENDRO_41 += DENDRO_69;

// initialize reduction 
DENDRO_12 = 0;
DENDRO_12 += DENDRO_7;
DENDRO_12 += DENDRO_24;
DENDRO_12 += DENDRO_33;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_3;
DENDRO_7 += DENDRO_27;
DENDRO_7 += DENDRO_30;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_14;
DENDRO_3 += DENDRO_13;
DENDRO_3 += DENDRO_36;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_16;
DENDRO_13 += DENDRO_15;
DENDRO_13 += DENDRO_38;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_19;
DENDRO_14 += DENDRO_17;
DENDRO_14 += DENDRO_42;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_84;
DENDRO_15 += DENDRO_83;
DENDRO_15 += DENDRO_81;
DENDRO_15 += DENDRO_80;
DENDRO_15 += DENDRO_78;
DENDRO_15 += DENDRO_76;
DENDRO_15 += DENDRO_74;
DENDRO_15 += DENDRO_72;
DENDRO_15 += DENDRO_46;
DENDRO_16 = chi[pp];

// int-power reduce pow(chi[pp],-1)
DENDRO_17 = DENDRO_16;
DENDRO_17 = 1.0/DENDRO_17;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_91;
DENDRO_16 += DENDRO_90;
DENDRO_16 += DENDRO_89;
DENDRO_16 += DENDRO_88;
DENDRO_16 += DENDRO_87;
DENDRO_16 += DENDRO_86;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_101;
DENDRO_19 += DENDRO_100;
DENDRO_19 += DENDRO_99;
DENDRO_19 += DENDRO_97;
DENDRO_19 += DENDRO_96;
DENDRO_19 += DENDRO_95;
DENDRO_19 += DENDRO_94;
DENDRO_19 += DENDRO_93;
DENDRO_19 += DENDRO_92;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_18;
DENDRO_24 *= DENDRO_77;
DENDRO_24 *= DENDRO_70;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_71;
DENDRO_27 *= DENDRO_77;
DENDRO_27 *= DENDRO_2;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_71;
DENDRO_30 *= DENDRO_22;
DENDRO_30 *= DENDRO_70;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_73;
DENDRO_33 *= DENDRO_22;
DENDRO_33 *= DENDRO_2;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_75;
DENDRO_36 *= DENDRO_98;
DENDRO_36 *= DENDRO_70;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_79;
DENDRO_38 *= DENDRO_22;
DENDRO_38 *= DENDRO_77;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_79;
DENDRO_42 *= DENDRO_98;
DENDRO_42 *= DENDRO_2;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_82;
DENDRO_43 *= DENDRO_98;
DENDRO_43 *= DENDRO_77;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_75;
DENDRO_44 *= DENDRO_102;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_85;
DENDRO_45 *= DENDRO_71;
DENDRO_45 *= DENDRO_2;
DENDRO_45 *= DENDRO_70;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_85;
DENDRO_46 *= DENDRO_75;
DENDRO_46 *= DENDRO_77;
DENDRO_46 *= DENDRO_70;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_85;
DENDRO_47 *= DENDRO_79;
DENDRO_47 *= DENDRO_77;
DENDRO_47 *= DENDRO_2;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_18;
DENDRO_48 *= DENDRO_103;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_73;
DENDRO_49 *= DENDRO_4;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_82;
DENDRO_50 *= DENDRO_102;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_85;
DENDRO_51 *= DENDRO_71;
DENDRO_51 *= DENDRO_22;
DENDRO_51 *= DENDRO_77;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_85;
DENDRO_52 *= DENDRO_75;
DENDRO_52 *= DENDRO_98;
DENDRO_52 *= DENDRO_77;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_85;
DENDRO_53 *= DENDRO_79;
DENDRO_53 *= DENDRO_98;
DENDRO_53 *= DENDRO_22;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_18;
DENDRO_54 *= DENDRO_102;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_73;
DENDRO_55 *= DENDRO_23;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_82;
DENDRO_56 *= DENDRO_104;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_85;
DENDRO_57 *= DENDRO_10;
DENDRO_57 *= DENDRO_20;
DENDRO_57 *= DENDRO_105;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_85;
DENDRO_58 *= DENDRO_10;
DENDRO_58 *= DENDRO_9;
DENDRO_58 *= DENDRO_25;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_85;
DENDRO_59 *= DENDRO_10;
DENDRO_59 *= DENDRO_11;
DENDRO_59 *= DENDRO_26;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_10;
DENDRO_60 *= DENDRO_1;
DENDRO_60 *= DENDRO_28;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_10;
DENDRO_61 *= DENDRO_8;
DENDRO_61 *= DENDRO_29;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_10;
DENDRO_62 *= DENDRO_0;
DENDRO_62 *= DENDRO_31;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_85;
DENDRO_63 *= DENDRO_10;
DENDRO_63 *= DENDRO_20;
DENDRO_63 *= DENDRO_32;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_85;
DENDRO_32 *= DENDRO_10;
DENDRO_32 *= DENDRO_9;
DENDRO_32 *= DENDRO_34;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_85;
DENDRO_34 *= DENDRO_10;
DENDRO_34 *= DENDRO_11;
DENDRO_34 *= DENDRO_35;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_10;
DENDRO_35 *= DENDRO_1;
DENDRO_35 *= DENDRO_37;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_10;
DENDRO_37 *= DENDRO_8;
DENDRO_37 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_10;
DENDRO_39 *= DENDRO_0;
DENDRO_39 *= DENDRO_40;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_85;
DENDRO_40 *= DENDRO_18;
DENDRO_40 *= DENDRO_2;
DENDRO_40 *= DENDRO_70;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_85;
DENDRO_64 *= DENDRO_71;
DENDRO_64 *= DENDRO_5;
DENDRO_64 *= DENDRO_70;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_85;
DENDRO_65 *= DENDRO_73;
DENDRO_65 *= DENDRO_5;
DENDRO_65 *= DENDRO_2;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_85;
DENDRO_66 *= DENDRO_75;
DENDRO_66 *= DENDRO_77;
DENDRO_66 *= DENDRO_2;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_85;
DENDRO_67 *= DENDRO_75;
DENDRO_67 *= DENDRO_22;
DENDRO_67 *= DENDRO_70;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_85;
DENDRO_68 *= DENDRO_79;
DENDRO_68 *= DENDRO_77;
DENDRO_68 *= DENDRO_5;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_85;
DENDRO_69 *= DENDRO_79;
DENDRO_69 *= DENDRO_22;
DENDRO_69 *= DENDRO_2;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_85;
DENDRO_72 *= DENDRO_82;
DENDRO_72 *= DENDRO_22;
DENDRO_72 *= DENDRO_77;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_85;
DENDRO_74 *= DENDRO_71;
DENDRO_74 *= DENDRO_4;
DENDRO_76 = 4;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_76;
DENDRO_78 *= DENDRO_71;
DENDRO_78 *= DENDRO_5;
DENDRO_78 *= DENDRO_2;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_76;
DENDRO_80 *= DENDRO_75;
DENDRO_80 *= DENDRO_22;
DENDRO_80 *= DENDRO_2;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_76;
DENDRO_81 *= DENDRO_79;
DENDRO_81 *= DENDRO_22;
DENDRO_81 *= DENDRO_5;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_85;
DENDRO_83 *= DENDRO_18;
DENDRO_83 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_85;
DENDRO_4 *= DENDRO_73;
DENDRO_4 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_85;
DENDRO_21 *= DENDRO_82;
DENDRO_21 *= DENDRO_23;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_85;
DENDRO_84 *= DENDRO_10;
DENDRO_84 *= DENDRO_20;
DENDRO_84 *= DENDRO_41;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_85;
DENDRO_20 *= DENDRO_10;
DENDRO_20 *= DENDRO_9;
DENDRO_20 *= DENDRO_12;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_85;
DENDRO_9 *= DENDRO_10;
DENDRO_9 *= DENDRO_11;
DENDRO_9 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_10;
DENDRO_7 *= DENDRO_1;
DENDRO_7 *= DENDRO_3;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_10;
DENDRO_1 *= DENDRO_8;
DENDRO_1 *= DENDRO_13;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_10;
DENDRO_3 *= DENDRO_0;
DENDRO_3 *= DENDRO_14;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_85;
DENDRO_0 *= DENDRO_18;
DENDRO_0 *= DENDRO_77;
DENDRO_0 *= DENDRO_2;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_85;
DENDRO_8 *= DENDRO_71;
DENDRO_8 *= DENDRO_77;
DENDRO_8 *= DENDRO_5;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_85;
DENDRO_10 *= DENDRO_71;
DENDRO_10 *= DENDRO_22;
DENDRO_10 *= DENDRO_2;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_85;
DENDRO_11 *= DENDRO_73;
DENDRO_11 *= DENDRO_22;
DENDRO_11 *= DENDRO_5;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_85;
DENDRO_12 *= DENDRO_75;
DENDRO_12 *= DENDRO_22;
DENDRO_12 *= DENDRO_77;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_85;
DENDRO_13 *= DENDRO_75;
DENDRO_13 *= DENDRO_98;
DENDRO_13 *= DENDRO_2;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_85;
DENDRO_14 *= DENDRO_79;
DENDRO_14 *= DENDRO_98;
DENDRO_14 *= DENDRO_5;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_85;
DENDRO_18 *= DENDRO_82;
DENDRO_18 *= DENDRO_98;
DENDRO_18 *= DENDRO_22;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_85;
DENDRO_41 *= DENDRO_79;
DENDRO_41 *= DENDRO_23;
DENDRO_23 = 3;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_23;
DENDRO_73 = grad_2_chi;
DENDRO_71 *= DENDRO_73;
DENDRO_71 *= DENDRO_17;
DENDRO_71 *= DENDRO_15;
DENDRO_73 = 4.0/3.0;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_73;
DENDRO_75 *= DENDRO_2;
DENDRO_79 = grad_2_K;
DENDRO_75 *= DENDRO_79;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_23;
DENDRO_82 = grad_1_chi;
DENDRO_79 *= DENDRO_82;
DENDRO_79 *= DENDRO_17;
DENDRO_79 *= DENDRO_16;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_73;
DENDRO_82 *= DENDRO_5;
DENDRO_86 = grad_1_K;
DENDRO_82 *= DENDRO_86;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_23;
DENDRO_23 = grad_0_chi;
DENDRO_86 *= DENDRO_23;
DENDRO_86 *= DENDRO_17;
DENDRO_86 *= DENDRO_19;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_73;
DENDRO_17 *= DENDRO_22;
DENDRO_23 = grad_0_K;
DENDRO_17 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_44;
DENDRO_23 += DENDRO_43;
DENDRO_23 += DENDRO_42;
DENDRO_23 += DENDRO_38;
DENDRO_23 += DENDRO_36;
DENDRO_23 += DENDRO_33;
DENDRO_23 += DENDRO_30;
DENDRO_23 += DENDRO_27;
DENDRO_23 += DENDRO_24;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_50;
DENDRO_24 += DENDRO_49;
DENDRO_24 += DENDRO_48;
DENDRO_24 += DENDRO_47;
DENDRO_24 += DENDRO_46;
DENDRO_24 += DENDRO_45;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_56;
DENDRO_27 += DENDRO_55;
DENDRO_27 += DENDRO_54;
DENDRO_27 += DENDRO_53;
DENDRO_27 += DENDRO_52;
DENDRO_27 += DENDRO_51;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_62;
DENDRO_30 += DENDRO_61;
DENDRO_30 += DENDRO_60;
DENDRO_30 += DENDRO_59;
DENDRO_30 += DENDRO_58;
DENDRO_30 += DENDRO_57;
DENDRO_33 = grad_0_beta0;

// initialize reduction 
DENDRO_36 = 0;
DENDRO_36 += DENDRO_33;
DENDRO_33 = grad_1_beta1;
DENDRO_36 += DENDRO_33;
DENDRO_38 = grad_2_beta2;
DENDRO_36 += DENDRO_38;

// initialize reduction 
DENDRO_38 = 0;
DENDRO_38 += DENDRO_39;
DENDRO_38 += DENDRO_37;
DENDRO_38 += DENDRO_35;
DENDRO_38 += DENDRO_34;
DENDRO_38 += DENDRO_32;
DENDRO_38 += DENDRO_63;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_74;
DENDRO_32 += DENDRO_72;
DENDRO_32 += DENDRO_69;
DENDRO_32 += DENDRO_68;
DENDRO_32 += DENDRO_67;
DENDRO_32 += DENDRO_66;
DENDRO_32 += DENDRO_65;
DENDRO_32 += DENDRO_64;
DENDRO_32 += DENDRO_40;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_21;
DENDRO_34 += DENDRO_4;
DENDRO_34 += DENDRO_83;
DENDRO_34 += DENDRO_81;
DENDRO_34 += DENDRO_80;
DENDRO_34 += DENDRO_78;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_3;
DENDRO_4 += DENDRO_1;
DENDRO_4 += DENDRO_7;
DENDRO_4 += DENDRO_9;
DENDRO_4 += DENDRO_20;
DENDRO_4 += DENDRO_84;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_41;
DENDRO_1 += DENDRO_18;
DENDRO_1 += DENDRO_14;
DENDRO_1 += DENDRO_13;
DENDRO_1 += DENDRO_12;
DENDRO_1 += DENDRO_11;
DENDRO_1 += DENDRO_10;
DENDRO_1 += DENDRO_8;
DENDRO_1 += DENDRO_0;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_75;
DENDRO_0 += DENDRO_71;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_82;
DENDRO_3 += DENDRO_79;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_17;
DENDRO_7 += DENDRO_86;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_76;
DENDRO_9 = alpha[pp];
DENDRO_8 *= DENDRO_9;
DENDRO_8 *= DENDRO_26;
DENDRO_8 *= DENDRO_15;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_76;
DENDRO_10 *= DENDRO_9;
DENDRO_10 *= DENDRO_25;
DENDRO_10 *= DENDRO_19;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_76;
DENDRO_11 *= DENDRO_9;
DENDRO_11 *= DENDRO_105;
DENDRO_11 *= DENDRO_23;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_85;
DENDRO_12 *= DENDRO_9;
DENDRO_12 *= DENDRO_29;
DENDRO_12 *= DENDRO_16;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_85;
DENDRO_13 *= DENDRO_9;
DENDRO_13 *= DENDRO_31;
DENDRO_13 *= DENDRO_24;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_85;
DENDRO_14 *= DENDRO_9;
DENDRO_14 *= DENDRO_28;
DENDRO_14 *= DENDRO_27;
DENDRO_15 = 7.0/3.0;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_15;
DENDRO_16 *= DENDRO_2;
DENDRO_17 = grad2_1_2_beta1;
DENDRO_16 *= DENDRO_17;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_15;
DENDRO_17 *= DENDRO_22;
DENDRO_15 = grad2_0_1_beta1;
DENDRO_17 *= DENDRO_15;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_73;
DENDRO_15 *= DENDRO_5;
DENDRO_18 = grad2_1_1_beta1;
DENDRO_15 *= DENDRO_18;
DENDRO_18 = 2.0/3.0;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_18;
DENDRO_19 *= DENDRO_36;
DENDRO_19 *= DENDRO_30;
DENDRO_18 = 1.0/3.0;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_18;
DENDRO_20 *= DENDRO_2;
DENDRO_21 = grad2_2_2_beta2;
DENDRO_20 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_18;
DENDRO_21 *= DENDRO_2;
DENDRO_2 = grad2_0_2_beta0;
DENDRO_21 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_18;
DENDRO_2 *= DENDRO_5;
DENDRO_23 = grad2_1_2_beta2;
DENDRO_2 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_18;
DENDRO_23 *= DENDRO_5;
DENDRO_5 = grad2_0_1_beta0;
DENDRO_23 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_18;
DENDRO_5 *= DENDRO_22;
DENDRO_24 = grad2_0_2_beta2;
DENDRO_5 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_18;
DENDRO_24 *= DENDRO_22;
DENDRO_18 = grad2_0_0_beta0;
DENDRO_24 *= DENDRO_18;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_85;
DENDRO_18 *= DENDRO_77;
DENDRO_22 = grad2_0_2_beta1;
DENDRO_18 *= DENDRO_22;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_6;
DENDRO_25 = grad_2_beta1;
DENDRO_22 *= DENDRO_25;
DENDRO_22 *= DENDRO_38;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_6;
DENDRO_26 = grad_2_alpha;
DENDRO_25 *= DENDRO_26;
DENDRO_25 *= DENDRO_32;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_6;
DENDRO_26 *= DENDRO_33;
DENDRO_26 *= DENDRO_30;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_6;
DENDRO_28 = grad_1_alpha;
DENDRO_27 *= DENDRO_28;
DENDRO_27 *= DENDRO_34;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_6;
DENDRO_29 = grad_0_beta1;
DENDRO_28 *= DENDRO_29;
DENDRO_28 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_6;
DENDRO_29 = grad_0_alpha;
DENDRO_4 *= DENDRO_29;
DENDRO_4 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_6;
DENDRO_1 *= DENDRO_9;
DENDRO_1 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_6;
DENDRO_0 *= DENDRO_9;
DENDRO_0 *= DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_6;
DENDRO_3 *= DENDRO_9;
DENDRO_3 *= DENDRO_7;
DENDRO_6 = beta2[pp];

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_6 = grad_2_Gt1;
DENDRO_7 *= DENDRO_6;
DENDRO_6 = beta1[pp];

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_6;
DENDRO_6 = grad_1_Gt1;
DENDRO_9 *= DENDRO_6;
DENDRO_6 = beta0[pp];

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_6;
DENDRO_6 = grad_0_Gt1;
DENDRO_29 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_70;
DENDRO_30 = grad2_2_2_beta1;
DENDRO_6 *= DENDRO_30;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_98;
DENDRO_31 = grad2_0_0_beta1;
DENDRO_30 *= DENDRO_31;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_30;
DENDRO_31 += DENDRO_6;
DENDRO_31 += DENDRO_29;
DENDRO_31 += DENDRO_9;
DENDRO_31 += DENDRO_7;
DENDRO_31 += DENDRO_3;
DENDRO_31 += DENDRO_0;
DENDRO_31 += DENDRO_1;
DENDRO_31 += DENDRO_4;
DENDRO_31 += DENDRO_28;
DENDRO_31 += DENDRO_27;
DENDRO_31 += DENDRO_26;
DENDRO_31 += DENDRO_25;
DENDRO_31 += DENDRO_22;
DENDRO_31 += DENDRO_18;
DENDRO_31 += DENDRO_24;
DENDRO_31 += DENDRO_5;
DENDRO_31 += DENDRO_23;
DENDRO_31 += DENDRO_2;
DENDRO_31 += DENDRO_21;
DENDRO_31 += DENDRO_20;
DENDRO_31 += DENDRO_19;
DENDRO_31 += DENDRO_15;
DENDRO_31 += DENDRO_17;
DENDRO_31 += DENDRO_16;
DENDRO_31 += DENDRO_14;
DENDRO_31 += DENDRO_13;
DENDRO_31 += DENDRO_12;
DENDRO_31 += DENDRO_11;
DENDRO_31 += DENDRO_10;
DENDRO_31 += DENDRO_8;
Gt_rhs1[pp] = DENDRO_31;


}
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , Gt2, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_Gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_Gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_Gt2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_Gt2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_Gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_Gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_Gt2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_Gt2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_Gt2=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , Gt2, blk);
device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_Gt2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_Gt2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_Gt2=Du[gidx];
__syncthreads();

//|V|= 396 |E| =905
//topological generations = 10
// [0]
// traversal id=0 temp_var_count=106
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;double DENDRO_33;double DENDRO_34;double DENDRO_35;double DENDRO_36;double DENDRO_37;double DENDRO_38;double DENDRO_39;double DENDRO_40;double DENDRO_41;double DENDRO_42;double DENDRO_43;double DENDRO_44;double DENDRO_45;double DENDRO_46;double DENDRO_47;double DENDRO_48;double DENDRO_49;double DENDRO_50;double DENDRO_51;double DENDRO_52;double DENDRO_53;double DENDRO_54;double DENDRO_55;double DENDRO_56;double DENDRO_57;double DENDRO_58;double DENDRO_59;double DENDRO_60;double DENDRO_61;double DENDRO_62;double DENDRO_63;double DENDRO_64;double DENDRO_65;double DENDRO_66;double DENDRO_67;double DENDRO_68;double DENDRO_69;double DENDRO_70;double DENDRO_71;double DENDRO_72;double DENDRO_73;double DENDRO_74;double DENDRO_75;double DENDRO_76;double DENDRO_77;double DENDRO_78;double DENDRO_79;double DENDRO_80;double DENDRO_81;double DENDRO_82;double DENDRO_83;double DENDRO_84;double DENDRO_85;double DENDRO_86;double DENDRO_87;double DENDRO_88;double DENDRO_89;double DENDRO_90;double DENDRO_91;double DENDRO_92;double DENDRO_93;double DENDRO_94;double DENDRO_95;double DENDRO_96;double DENDRO_97;double DENDRO_98;double DENDRO_99;double DENDRO_100;double DENDRO_101;double DENDRO_102;double DENDRO_103;double DENDRO_104;double DENDRO_105;
DENDRO_0 = gt1[pp];

// int-power reduce pow(gt1[pp],2)
DENDRO_1 = DENDRO_0;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = gt2[pp];

// int-power reduce pow(gt2[pp],2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = gt4[pp];

// int-power reduce pow(gt4[pp],2)
DENDRO_5 = DENDRO_4;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = -2;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_4;
DENDRO_6 = -1;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_6;
DENDRO_9 = gt0[pp];
DENDRO_8 *= DENDRO_9;
DENDRO_10 = gt3[pp];
DENDRO_8 *= DENDRO_10;
DENDRO_11 = gt5[pp];
DENDRO_8 *= DENDRO_11;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_11;
DENDRO_12 *= DENDRO_1;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_10;
DENDRO_13 *= DENDRO_3;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_9;
DENDRO_14 *= DENDRO_5;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_6;
DENDRO_15 *= DENDRO_0;
DENDRO_15 *= DENDRO_4;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_2;
DENDRO_16 *= DENDRO_10;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_14;
DENDRO_17 += DENDRO_13;
DENDRO_17 += DENDRO_12;
DENDRO_17 += DENDRO_8;
DENDRO_17 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_9;
DENDRO_7 *= DENDRO_10;
DENDRO_8 = -0.5;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_8;
DENDRO_8 = grad_1_gt2;
DENDRO_12 *= DENDRO_8;
DENDRO_13 = 0.5;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_13;
DENDRO_13 = grad_2_gt1;
DENDRO_14 *= DENDRO_13;
DENDRO_18 = 0.5;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_18;
DENDRO_18 = grad_0_gt4;
DENDRO_19 *= DENDRO_18;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_6;
DENDRO_20 *= DENDRO_0;
DENDRO_20 *= DENDRO_2;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_9;
DENDRO_21 *= DENDRO_4;
DENDRO_22 = -0.5;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_22;
DENDRO_23 *= DENDRO_13;
DENDRO_13 = 0.5;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_13;
DENDRO_22 *= DENDRO_8;
DENDRO_8 = -0.5;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_8;
DENDRO_13 *= DENDRO_18;
DENDRO_8 = -0.5;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_8;
DENDRO_8 = grad_1_gt0;
DENDRO_18 *= DENDRO_8;
DENDRO_24 = 1.0;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_24;
DENDRO_24 = grad_0_gt1;
DENDRO_25 *= DENDRO_24;
DENDRO_24 = -0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_24;
DENDRO_24 = grad_2_gt0;
DENDRO_26 *= DENDRO_24;
DENDRO_27 = 1.0;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_27;
DENDRO_27 = grad_0_gt2;
DENDRO_28 *= DENDRO_27;
DENDRO_27 = -0.5;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_27;
DENDRO_27 = grad_0_gt3;
DENDRO_29 *= DENDRO_27;
DENDRO_30 = 1.0;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_30;
DENDRO_30 = grad_1_gt1;
DENDRO_31 *= DENDRO_30;
DENDRO_30 = -0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_30;
DENDRO_30 = grad_2_gt3;
DENDRO_32 *= DENDRO_30;
DENDRO_33 = 1.0;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_33;
DENDRO_33 = grad_1_gt4;
DENDRO_34 *= DENDRO_33;
DENDRO_33 = -0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_33;
DENDRO_33 = grad_1_gt5;
DENDRO_35 *= DENDRO_33;
DENDRO_36 = 1.0;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = grad_2_gt4;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = -0.5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_36;
DENDRO_36 = grad_0_gt5;
DENDRO_38 *= DENDRO_36;
DENDRO_39 = 1.0;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_39;
DENDRO_39 = grad_2_gt2;
DENDRO_40 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_6;
DENDRO_39 *= DENDRO_2;
DENDRO_39 *= DENDRO_4;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_2 *= DENDRO_11;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_6;
DENDRO_0 *= DENDRO_9;
DENDRO_0 *= DENDRO_11;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_6;
DENDRO_4 *= DENDRO_10;
DENDRO_4 *= DENDRO_11;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_16;
DENDRO_9 += DENDRO_15;

// int-power reduce pow(-gt0[pp]*gt3[pp]*gt5[pp] + gt0[pp]*gt4[pp]**2 + gt1[pp]**2*gt5[pp] - 2*gt1[pp]*gt2[pp]*gt4[pp] + gt2[pp]**2*gt3[pp],-1)
DENDRO_10 = DENDRO_17;
DENDRO_10 = 1.0/DENDRO_10;

// initialize reduction 
DENDRO_11 = 0;
DENDRO_11 += DENDRO_1;
DENDRO_11 += DENDRO_7;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_19;
DENDRO_1 += DENDRO_14;
DENDRO_1 += DENDRO_12;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_21;
DENDRO_7 += DENDRO_20;

// initialize reduction 
DENDRO_12 = 0;
DENDRO_12 += DENDRO_19;
DENDRO_12 += DENDRO_22;
DENDRO_12 += DENDRO_23;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_22;
DENDRO_15 += DENDRO_14;
DENDRO_15 += DENDRO_13;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_25;
DENDRO_13 += DENDRO_18;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_28;
DENDRO_14 += DENDRO_26;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_31;
DENDRO_16 += DENDRO_29;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_34;
DENDRO_17 += DENDRO_32;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_37;
DENDRO_18 += DENDRO_35;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_40;
DENDRO_19 += DENDRO_38;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_2;
DENDRO_20 += DENDRO_39;

// initialize reduction 
DENDRO_2 = 0;
DENDRO_2 += DENDRO_3;
DENDRO_2 += DENDRO_0;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_5;
DENDRO_0 += DENDRO_4;
DENDRO_3 = DENDRO_igt5;

// int-power reduce pow(DENDRO_igt5,2)
DENDRO_4 = DENDRO_3;
DENDRO_4 *= DENDRO_3;
DENDRO_5 = DENDRO_igt4;

// int-power reduce pow(DENDRO_igt4,2)
DENDRO_21 = DENDRO_5;
DENDRO_21 *= DENDRO_5;
DENDRO_22 = DENDRO_igt2;

// int-power reduce pow(DENDRO_igt2,2)
DENDRO_23 = DENDRO_22;
DENDRO_23 *= DENDRO_22;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_25;
DENDRO_26 *= DENDRO_24;
DENDRO_26 *= DENDRO_10;
DENDRO_26 *= DENDRO_9;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_25;
DENDRO_28 *= DENDRO_36;
DENDRO_28 *= DENDRO_10;
DENDRO_28 *= DENDRO_11;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_10;
DENDRO_25 *= DENDRO_7;
DENDRO_25 *= DENDRO_1;
DENDRO_29 = 0.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_29;
DENDRO_31 *= DENDRO_8;
DENDRO_31 *= DENDRO_10;
DENDRO_31 *= DENDRO_9;
DENDRO_29 = 0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_29;
DENDRO_32 *= DENDRO_27;
DENDRO_32 *= DENDRO_10;
DENDRO_32 *= DENDRO_7;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_10;
DENDRO_29 *= DENDRO_11;
DENDRO_29 *= DENDRO_12;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_34;
DENDRO_35 *= DENDRO_30;
DENDRO_35 *= DENDRO_10;
DENDRO_35 *= DENDRO_7;
DENDRO_34 = 0.5;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_34;
DENDRO_37 *= DENDRO_33;
DENDRO_37 *= DENDRO_10;
DENDRO_37 *= DENDRO_11;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_10;
DENDRO_34 *= DENDRO_9;
DENDRO_34 *= DENDRO_15;
DENDRO_38 = 0.5;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_38;
DENDRO_38 = grad_0_gt0;
DENDRO_39 *= DENDRO_38;
DENDRO_39 *= DENDRO_10;
DENDRO_39 *= DENDRO_9;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_10;
DENDRO_40 *= DENDRO_13;
DENDRO_40 *= DENDRO_7;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_10;
DENDRO_41 *= DENDRO_11;
DENDRO_41 *= DENDRO_14;
DENDRO_42 = 0.5;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_42;
DENDRO_42 = grad_1_gt3;
DENDRO_43 *= DENDRO_42;
DENDRO_43 *= DENDRO_10;
DENDRO_43 *= DENDRO_7;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_10;
DENDRO_44 *= DENDRO_16;
DENDRO_44 *= DENDRO_9;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_10;
DENDRO_45 *= DENDRO_11;
DENDRO_45 *= DENDRO_17;
DENDRO_46 = 0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_46;
DENDRO_46 = grad_2_gt5;
DENDRO_47 *= DENDRO_46;
DENDRO_47 *= DENDRO_10;
DENDRO_47 *= DENDRO_11;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_10;
DENDRO_48 *= DENDRO_18;
DENDRO_48 *= DENDRO_7;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_10;
DENDRO_49 *= DENDRO_19;
DENDRO_49 *= DENDRO_9;
DENDRO_50 = 0.5;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_50;
DENDRO_51 *= DENDRO_24;
DENDRO_51 *= DENDRO_10;
DENDRO_51 *= DENDRO_20;
DENDRO_50 = 0.5;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_50;
DENDRO_52 *= DENDRO_36;
DENDRO_52 *= DENDRO_10;
DENDRO_52 *= DENDRO_7;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_10;
DENDRO_50 *= DENDRO_2;
DENDRO_50 *= DENDRO_1;
DENDRO_53 = 0.5;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_53;
DENDRO_54 *= DENDRO_8;
DENDRO_54 *= DENDRO_10;
DENDRO_54 *= DENDRO_20;
DENDRO_53 = 0.5;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_53;
DENDRO_55 *= DENDRO_27;
DENDRO_55 *= DENDRO_10;
DENDRO_55 *= DENDRO_2;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_10;
DENDRO_53 *= DENDRO_7;
DENDRO_53 *= DENDRO_12;
DENDRO_56 = 0.5;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_56;
DENDRO_57 *= DENDRO_30;
DENDRO_57 *= DENDRO_10;
DENDRO_57 *= DENDRO_2;
DENDRO_56 = 0.5;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_56;
DENDRO_58 *= DENDRO_33;
DENDRO_58 *= DENDRO_10;
DENDRO_58 *= DENDRO_7;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_10;
DENDRO_56 *= DENDRO_20;
DENDRO_56 *= DENDRO_15;
DENDRO_59 = 0.5;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_59;
DENDRO_60 *= DENDRO_38;
DENDRO_60 *= DENDRO_10;
DENDRO_60 *= DENDRO_20;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_10;
DENDRO_59 *= DENDRO_14;
DENDRO_59 *= DENDRO_7;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_10;
DENDRO_61 *= DENDRO_2;
DENDRO_61 *= DENDRO_13;
DENDRO_62 = 0.5;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_62;
DENDRO_63 *= DENDRO_42;
DENDRO_63 *= DENDRO_10;
DENDRO_63 *= DENDRO_2;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_10;
DENDRO_62 *= DENDRO_17;
DENDRO_62 *= DENDRO_7;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_10;
DENDRO_64 *= DENDRO_16;
DENDRO_64 *= DENDRO_20;
DENDRO_65 = 0.5;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_65;
DENDRO_66 *= DENDRO_46;
DENDRO_66 *= DENDRO_10;
DENDRO_66 *= DENDRO_7;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_10;
DENDRO_65 *= DENDRO_19;
DENDRO_65 *= DENDRO_20;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_10;
DENDRO_67 *= DENDRO_2;
DENDRO_67 *= DENDRO_18;
DENDRO_68 = 0.5;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_68;
DENDRO_69 *= DENDRO_24;
DENDRO_69 *= DENDRO_10;
DENDRO_69 *= DENDRO_0;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_24;
DENDRO_68 *= DENDRO_36;
DENDRO_68 *= DENDRO_10;
DENDRO_68 *= DENDRO_9;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_10;
DENDRO_24 *= DENDRO_20;
DENDRO_24 *= DENDRO_1;
DENDRO_1 = 0.5;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_1;
DENDRO_36 *= DENDRO_8;
DENDRO_36 *= DENDRO_10;
DENDRO_36 *= DENDRO_0;
DENDRO_1 = 0.5;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_1;
DENDRO_8 *= DENDRO_27;
DENDRO_8 *= DENDRO_10;
DENDRO_8 *= DENDRO_20;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_10;
DENDRO_1 *= DENDRO_9;
DENDRO_1 *= DENDRO_12;
DENDRO_12 = 0.5;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_12;
DENDRO_27 *= DENDRO_30;
DENDRO_27 *= DENDRO_10;
DENDRO_27 *= DENDRO_20;
DENDRO_12 = 0.5;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_12;
DENDRO_30 *= DENDRO_33;
DENDRO_30 *= DENDRO_10;
DENDRO_30 *= DENDRO_9;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_10;
DENDRO_12 *= DENDRO_0;
DENDRO_12 *= DENDRO_15;
DENDRO_15 = 0.5;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_15;
DENDRO_33 *= DENDRO_38;
DENDRO_33 *= DENDRO_10;
DENDRO_33 *= DENDRO_0;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_10;
DENDRO_15 *= DENDRO_14;
DENDRO_15 *= DENDRO_9;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_10;
DENDRO_14 *= DENDRO_13;
DENDRO_14 *= DENDRO_20;
DENDRO_13 = 0.5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_13;
DENDRO_38 *= DENDRO_42;
DENDRO_38 *= DENDRO_10;
DENDRO_38 *= DENDRO_20;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_10;
DENDRO_13 *= DENDRO_17;
DENDRO_13 *= DENDRO_9;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_10;
DENDRO_17 *= DENDRO_0;
DENDRO_17 *= DENDRO_16;
DENDRO_16 = 0.5;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_16;
DENDRO_42 *= DENDRO_46;
DENDRO_42 *= DENDRO_10;
DENDRO_42 *= DENDRO_9;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_10;
DENDRO_16 *= DENDRO_18;
DENDRO_16 *= DENDRO_20;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_10;
DENDRO_18 *= DENDRO_0;
DENDRO_18 *= DENDRO_19;
DENDRO_19 = 2;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_19;
DENDRO_70 = At4[pp];
DENDRO_46 *= DENDRO_70;
DENDRO_46 *= DENDRO_5;
DENDRO_46 *= DENDRO_3;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_19;
DENDRO_72 = At2[pp];
DENDRO_71 *= DENDRO_72;
DENDRO_71 *= DENDRO_22;
DENDRO_71 *= DENDRO_3;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_19;
DENDRO_74 = At1[pp];
DENDRO_73 *= DENDRO_74;
DENDRO_73 *= DENDRO_22;
DENDRO_73 *= DENDRO_5;
DENDRO_75 = At5[pp];

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_75;
DENDRO_76 *= DENDRO_4;
DENDRO_77 = At3[pp];

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_77;
DENDRO_78 *= DENDRO_21;
DENDRO_79 = At0[pp];

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_79;
DENDRO_80 *= DENDRO_23;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_75;
DENDRO_81 *= DENDRO_5;
DENDRO_81 *= DENDRO_3;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_70;
DENDRO_83 = DENDRO_igt3;
DENDRO_82 *= DENDRO_83;
DENDRO_82 *= DENDRO_3;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_77;
DENDRO_84 *= DENDRO_83;
DENDRO_84 *= DENDRO_5;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_72;
DENDRO_85 *= DENDRO_22;
DENDRO_85 *= DENDRO_5;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_72;
DENDRO_87 = DENDRO_igt1;
DENDRO_86 *= DENDRO_87;
DENDRO_86 *= DENDRO_3;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_74;
DENDRO_88 *= DENDRO_22;
DENDRO_88 *= DENDRO_83;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_74;
DENDRO_89 *= DENDRO_87;
DENDRO_89 *= DENDRO_5;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_79;
DENDRO_90 *= DENDRO_87;
DENDRO_90 *= DENDRO_22;

// initialize reduction 
DENDRO_91 = 1;
DENDRO_91 *= DENDRO_70;
DENDRO_91 *= DENDRO_21;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_75;
DENDRO_92 *= DENDRO_22;
DENDRO_92 *= DENDRO_3;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_70;
DENDRO_93 *= DENDRO_22;
DENDRO_93 *= DENDRO_5;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_70;
DENDRO_94 *= DENDRO_87;
DENDRO_94 *= DENDRO_3;

// initialize reduction 
DENDRO_95 = 1;
DENDRO_95 *= DENDRO_77;
DENDRO_95 *= DENDRO_87;
DENDRO_95 *= DENDRO_5;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_72;
DENDRO_97 = DENDRO_igt0;
DENDRO_96 *= DENDRO_97;
DENDRO_96 *= DENDRO_3;

// initialize reduction 
DENDRO_98 = 1;
DENDRO_98 *= DENDRO_74;
DENDRO_98 *= DENDRO_87;
DENDRO_98 *= DENDRO_22;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_74;
DENDRO_99 *= DENDRO_97;
DENDRO_99 *= DENDRO_5;

// initialize reduction 
DENDRO_100 = 1;
DENDRO_100 *= DENDRO_79;
DENDRO_100 *= DENDRO_97;
DENDRO_100 *= DENDRO_22;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_72;
DENDRO_101 *= DENDRO_23;

// int-power reduce pow(DENDRO_igt1,2)
DENDRO_102 = DENDRO_87;
DENDRO_102 *= DENDRO_87;

// int-power reduce pow(DENDRO_igt3,2)
DENDRO_103 = DENDRO_83;
DENDRO_103 *= DENDRO_83;

// int-power reduce pow(DENDRO_igt0,2)
DENDRO_104 = DENDRO_97;
DENDRO_104 *= DENDRO_97;

// initialize reduction 
DENDRO_105 = 0;
DENDRO_105 += DENDRO_25;
DENDRO_105 += DENDRO_28;
DENDRO_105 += DENDRO_26;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_29;
DENDRO_25 += DENDRO_32;
DENDRO_25 += DENDRO_31;

// initialize reduction 
DENDRO_26 = 0;
DENDRO_26 += DENDRO_34;
DENDRO_26 += DENDRO_37;
DENDRO_26 += DENDRO_35;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_41;
DENDRO_28 += DENDRO_40;
DENDRO_28 += DENDRO_39;

// initialize reduction 
DENDRO_29 = 0;
DENDRO_29 += DENDRO_45;
DENDRO_29 += DENDRO_44;
DENDRO_29 += DENDRO_43;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_49;
DENDRO_31 += DENDRO_48;
DENDRO_31 += DENDRO_47;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_50;
DENDRO_32 += DENDRO_52;
DENDRO_32 += DENDRO_51;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_53;
DENDRO_34 += DENDRO_55;
DENDRO_34 += DENDRO_54;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_56;
DENDRO_35 += DENDRO_58;
DENDRO_35 += DENDRO_57;

// initialize reduction 
DENDRO_37 = 0;
DENDRO_37 += DENDRO_61;
DENDRO_37 += DENDRO_59;
DENDRO_37 += DENDRO_60;

// initialize reduction 
DENDRO_39 = 0;
DENDRO_39 += DENDRO_64;
DENDRO_39 += DENDRO_62;
DENDRO_39 += DENDRO_63;

// initialize reduction 
DENDRO_40 = 0;
DENDRO_40 += DENDRO_67;
DENDRO_40 += DENDRO_65;
DENDRO_40 += DENDRO_66;

// initialize reduction 
DENDRO_41 = 0;
DENDRO_41 += DENDRO_24;
DENDRO_41 += DENDRO_68;
DENDRO_41 += DENDRO_69;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_1;
DENDRO_24 += DENDRO_8;
DENDRO_24 += DENDRO_36;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_12;
DENDRO_1 += DENDRO_30;
DENDRO_1 += DENDRO_27;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_14;
DENDRO_8 += DENDRO_15;
DENDRO_8 += DENDRO_33;

// initialize reduction 
DENDRO_12 = 0;
DENDRO_12 += DENDRO_17;
DENDRO_12 += DENDRO_13;
DENDRO_12 += DENDRO_38;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_18;
DENDRO_13 += DENDRO_16;
DENDRO_13 += DENDRO_42;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_80;
DENDRO_14 += DENDRO_78;
DENDRO_14 += DENDRO_76;
DENDRO_14 += DENDRO_73;
DENDRO_14 += DENDRO_71;
DENDRO_14 += DENDRO_46;
DENDRO_15 = chi[pp];

// int-power reduce pow(chi[pp],-1)
DENDRO_16 = DENDRO_15;
DENDRO_16 = 1.0/DENDRO_16;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_91;
DENDRO_15 += DENDRO_90;
DENDRO_15 += DENDRO_89;
DENDRO_15 += DENDRO_88;
DENDRO_15 += DENDRO_86;
DENDRO_15 += DENDRO_85;
DENDRO_15 += DENDRO_84;
DENDRO_15 += DENDRO_82;
DENDRO_15 += DENDRO_81;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_101;
DENDRO_17 += DENDRO_100;
DENDRO_17 += DENDRO_99;
DENDRO_17 += DENDRO_98;
DENDRO_17 += DENDRO_96;
DENDRO_17 += DENDRO_95;
DENDRO_17 += DENDRO_94;
DENDRO_17 += DENDRO_93;
DENDRO_17 += DENDRO_92;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_75;
DENDRO_18 *= DENDRO_22;
DENDRO_18 *= DENDRO_5;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_70;
DENDRO_27 *= DENDRO_22;
DENDRO_27 *= DENDRO_83;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_70;
DENDRO_30 *= DENDRO_87;
DENDRO_30 *= DENDRO_5;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_77;
DENDRO_33 *= DENDRO_87;
DENDRO_33 *= DENDRO_83;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_72;
DENDRO_36 *= DENDRO_87;
DENDRO_36 *= DENDRO_22;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_72;
DENDRO_38 *= DENDRO_97;
DENDRO_38 *= DENDRO_5;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_74;
DENDRO_42 *= DENDRO_97;
DENDRO_42 *= DENDRO_83;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_79;
DENDRO_43 *= DENDRO_97;
DENDRO_43 *= DENDRO_87;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_74;
DENDRO_44 *= DENDRO_102;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_19;
DENDRO_45 *= DENDRO_70;
DENDRO_45 *= DENDRO_83;
DENDRO_45 *= DENDRO_5;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_19;
DENDRO_46 *= DENDRO_72;
DENDRO_46 *= DENDRO_87;
DENDRO_46 *= DENDRO_5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_19;
DENDRO_47 *= DENDRO_74;
DENDRO_47 *= DENDRO_87;
DENDRO_47 *= DENDRO_83;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_75;
DENDRO_48 *= DENDRO_21;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_77;
DENDRO_49 *= DENDRO_103;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_79;
DENDRO_50 *= DENDRO_102;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_19;
DENDRO_51 *= DENDRO_70;
DENDRO_51 *= DENDRO_87;
DENDRO_51 *= DENDRO_22;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_19;
DENDRO_52 *= DENDRO_72;
DENDRO_52 *= DENDRO_97;
DENDRO_52 *= DENDRO_22;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_19;
DENDRO_53 *= DENDRO_74;
DENDRO_53 *= DENDRO_97;
DENDRO_53 *= DENDRO_87;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_75;
DENDRO_54 *= DENDRO_23;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_77;
DENDRO_55 *= DENDRO_102;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_79;
DENDRO_56 *= DENDRO_104;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_19;
DENDRO_57 *= DENDRO_10;
DENDRO_57 *= DENDRO_9;
DENDRO_57 *= DENDRO_105;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_19;
DENDRO_58 *= DENDRO_10;
DENDRO_58 *= DENDRO_20;
DENDRO_58 *= DENDRO_25;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_19;
DENDRO_59 *= DENDRO_10;
DENDRO_59 *= DENDRO_7;
DENDRO_59 *= DENDRO_26;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_10;
DENDRO_60 *= DENDRO_0;
DENDRO_60 *= DENDRO_28;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_10;
DENDRO_61 *= DENDRO_2;
DENDRO_61 *= DENDRO_29;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_10;
DENDRO_62 *= DENDRO_11;
DENDRO_62 *= DENDRO_31;
DENDRO_63 = 4;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_63;
DENDRO_64 *= DENDRO_70;
DENDRO_64 *= DENDRO_5;
DENDRO_64 *= DENDRO_3;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_63;
DENDRO_65 *= DENDRO_72;
DENDRO_65 *= DENDRO_22;
DENDRO_65 *= DENDRO_3;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_63;
DENDRO_66 *= DENDRO_74;
DENDRO_66 *= DENDRO_22;
DENDRO_66 *= DENDRO_5;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_19;
DENDRO_67 *= DENDRO_75;
DENDRO_67 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_19;
DENDRO_4 *= DENDRO_77;
DENDRO_4 *= DENDRO_21;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_19;
DENDRO_68 *= DENDRO_79;
DENDRO_68 *= DENDRO_23;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_19;
DENDRO_69 *= DENDRO_10;
DENDRO_69 *= DENDRO_9;
DENDRO_69 *= DENDRO_32;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_19;
DENDRO_32 *= DENDRO_10;
DENDRO_32 *= DENDRO_20;
DENDRO_32 *= DENDRO_34;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_19;
DENDRO_34 *= DENDRO_10;
DENDRO_34 *= DENDRO_7;
DENDRO_34 *= DENDRO_35;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_10;
DENDRO_35 *= DENDRO_0;
DENDRO_35 *= DENDRO_37;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_10;
DENDRO_37 *= DENDRO_2;
DENDRO_37 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_10;
DENDRO_39 *= DENDRO_11;
DENDRO_39 *= DENDRO_40;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_19;
DENDRO_40 *= DENDRO_75;
DENDRO_40 *= DENDRO_5;
DENDRO_40 *= DENDRO_3;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_19;
DENDRO_71 *= DENDRO_70;
DENDRO_71 *= DENDRO_83;
DENDRO_71 *= DENDRO_3;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_19;
DENDRO_73 *= DENDRO_77;
DENDRO_73 *= DENDRO_83;
DENDRO_73 *= DENDRO_5;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_19;
DENDRO_76 *= DENDRO_72;
DENDRO_76 *= DENDRO_22;
DENDRO_76 *= DENDRO_5;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_19;
DENDRO_78 *= DENDRO_72;
DENDRO_78 *= DENDRO_87;
DENDRO_78 *= DENDRO_3;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_19;
DENDRO_80 *= DENDRO_74;
DENDRO_80 *= DENDRO_22;
DENDRO_80 *= DENDRO_83;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_19;
DENDRO_81 *= DENDRO_74;
DENDRO_81 *= DENDRO_87;
DENDRO_81 *= DENDRO_5;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_19;
DENDRO_82 *= DENDRO_79;
DENDRO_82 *= DENDRO_87;
DENDRO_82 *= DENDRO_22;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_19;
DENDRO_84 *= DENDRO_70;
DENDRO_84 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_19;
DENDRO_21 *= DENDRO_10;
DENDRO_21 *= DENDRO_9;
DENDRO_21 *= DENDRO_41;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_19;
DENDRO_9 *= DENDRO_10;
DENDRO_9 *= DENDRO_20;
DENDRO_9 *= DENDRO_24;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_19;
DENDRO_20 *= DENDRO_10;
DENDRO_20 *= DENDRO_7;
DENDRO_20 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_10;
DENDRO_1 *= DENDRO_0;
DENDRO_1 *= DENDRO_8;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_10;
DENDRO_0 *= DENDRO_2;
DENDRO_0 *= DENDRO_12;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_10;
DENDRO_2 *= DENDRO_11;
DENDRO_2 *= DENDRO_13;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_19;
DENDRO_7 *= DENDRO_75;
DENDRO_7 *= DENDRO_22;
DENDRO_7 *= DENDRO_3;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_19;
DENDRO_8 *= DENDRO_70;
DENDRO_8 *= DENDRO_22;
DENDRO_8 *= DENDRO_5;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_19;
DENDRO_10 *= DENDRO_70;
DENDRO_10 *= DENDRO_87;
DENDRO_10 *= DENDRO_3;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_19;
DENDRO_11 *= DENDRO_77;
DENDRO_11 *= DENDRO_87;
DENDRO_11 *= DENDRO_5;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_19;
DENDRO_12 *= DENDRO_72;
DENDRO_12 *= DENDRO_97;
DENDRO_12 *= DENDRO_3;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_19;
DENDRO_13 *= DENDRO_74;
DENDRO_13 *= DENDRO_87;
DENDRO_13 *= DENDRO_22;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_19;
DENDRO_24 *= DENDRO_74;
DENDRO_24 *= DENDRO_97;
DENDRO_24 *= DENDRO_5;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_19;
DENDRO_41 *= DENDRO_79;
DENDRO_41 *= DENDRO_97;
DENDRO_41 *= DENDRO_22;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_19;
DENDRO_70 *= DENDRO_72;
DENDRO_70 *= DENDRO_23;
DENDRO_23 = 3;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_23;
DENDRO_74 = grad_2_chi;
DENDRO_72 *= DENDRO_74;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_14;
DENDRO_74 = 4.0/3.0;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_74;
DENDRO_75 *= DENDRO_3;
DENDRO_77 = grad_2_K;
DENDRO_75 *= DENDRO_77;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_23;
DENDRO_79 = grad_1_chi;
DENDRO_77 *= DENDRO_79;
DENDRO_77 *= DENDRO_16;
DENDRO_77 *= DENDRO_15;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_74;
DENDRO_79 *= DENDRO_5;
DENDRO_85 = grad_1_K;
DENDRO_79 *= DENDRO_85;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_23;
DENDRO_23 = grad_0_chi;
DENDRO_85 *= DENDRO_23;
DENDRO_85 *= DENDRO_16;
DENDRO_85 *= DENDRO_17;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_74;
DENDRO_16 *= DENDRO_22;
DENDRO_23 = grad_0_K;
DENDRO_16 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_44;
DENDRO_23 += DENDRO_43;
DENDRO_23 += DENDRO_42;
DENDRO_23 += DENDRO_38;
DENDRO_23 += DENDRO_36;
DENDRO_23 += DENDRO_33;
DENDRO_23 += DENDRO_30;
DENDRO_23 += DENDRO_27;
DENDRO_23 += DENDRO_18;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_50;
DENDRO_18 += DENDRO_49;
DENDRO_18 += DENDRO_48;
DENDRO_18 += DENDRO_47;
DENDRO_18 += DENDRO_46;
DENDRO_18 += DENDRO_45;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_56;
DENDRO_27 += DENDRO_55;
DENDRO_27 += DENDRO_54;
DENDRO_27 += DENDRO_53;
DENDRO_27 += DENDRO_52;
DENDRO_27 += DENDRO_51;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_62;
DENDRO_30 += DENDRO_61;
DENDRO_30 += DENDRO_60;
DENDRO_30 += DENDRO_59;
DENDRO_30 += DENDRO_58;
DENDRO_30 += DENDRO_57;
DENDRO_33 = grad_0_beta0;

// initialize reduction 
DENDRO_36 = 0;
DENDRO_36 += DENDRO_33;
DENDRO_33 = grad_1_beta1;
DENDRO_36 += DENDRO_33;
DENDRO_33 = grad_2_beta2;
DENDRO_36 += DENDRO_33;

// initialize reduction 
DENDRO_38 = 0;
DENDRO_38 += DENDRO_68;
DENDRO_38 += DENDRO_4;
DENDRO_38 += DENDRO_67;
DENDRO_38 += DENDRO_66;
DENDRO_38 += DENDRO_65;
DENDRO_38 += DENDRO_64;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_39;
DENDRO_4 += DENDRO_37;
DENDRO_4 += DENDRO_35;
DENDRO_4 += DENDRO_34;
DENDRO_4 += DENDRO_32;
DENDRO_4 += DENDRO_69;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_84;
DENDRO_32 += DENDRO_82;
DENDRO_32 += DENDRO_81;
DENDRO_32 += DENDRO_80;
DENDRO_32 += DENDRO_78;
DENDRO_32 += DENDRO_76;
DENDRO_32 += DENDRO_73;
DENDRO_32 += DENDRO_71;
DENDRO_32 += DENDRO_40;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_2;
DENDRO_34 += DENDRO_0;
DENDRO_34 += DENDRO_1;
DENDRO_34 += DENDRO_20;
DENDRO_34 += DENDRO_9;
DENDRO_34 += DENDRO_21;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_70;
DENDRO_0 += DENDRO_41;
DENDRO_0 += DENDRO_24;
DENDRO_0 += DENDRO_13;
DENDRO_0 += DENDRO_12;
DENDRO_0 += DENDRO_11;
DENDRO_0 += DENDRO_10;
DENDRO_0 += DENDRO_8;
DENDRO_0 += DENDRO_7;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_75;
DENDRO_1 += DENDRO_72;

// initialize reduction 
DENDRO_2 = 0;
DENDRO_2 += DENDRO_79;
DENDRO_2 += DENDRO_77;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_16;
DENDRO_7 += DENDRO_85;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_63;
DENDRO_9 = alpha[pp];
DENDRO_8 *= DENDRO_9;
DENDRO_8 *= DENDRO_26;
DENDRO_8 *= DENDRO_15;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_63;
DENDRO_10 *= DENDRO_9;
DENDRO_10 *= DENDRO_105;
DENDRO_10 *= DENDRO_17;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_63;
DENDRO_11 *= DENDRO_9;
DENDRO_11 *= DENDRO_25;
DENDRO_11 *= DENDRO_23;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_19;
DENDRO_12 *= DENDRO_9;
DENDRO_12 *= DENDRO_31;
DENDRO_12 *= DENDRO_14;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_19;
DENDRO_13 *= DENDRO_9;
DENDRO_13 *= DENDRO_29;
DENDRO_13 *= DENDRO_18;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_19;
DENDRO_14 *= DENDRO_9;
DENDRO_14 *= DENDRO_28;
DENDRO_14 *= DENDRO_27;
DENDRO_15 = 7.0/3.0;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_15;
DENDRO_16 *= DENDRO_5;
DENDRO_17 = grad2_1_2_beta2;
DENDRO_16 *= DENDRO_17;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_15;
DENDRO_17 *= DENDRO_22;
DENDRO_15 = grad2_0_2_beta2;
DENDRO_17 *= DENDRO_15;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_74;
DENDRO_15 *= DENDRO_3;
DENDRO_18 = grad2_2_2_beta2;
DENDRO_15 *= DENDRO_18;
DENDRO_18 = 2.0/3.0;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_18;
DENDRO_20 *= DENDRO_36;
DENDRO_20 *= DENDRO_30;
DENDRO_18 = 1.0/3.0;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_18;
DENDRO_21 *= DENDRO_3;
DENDRO_23 = grad2_1_2_beta1;
DENDRO_21 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_18;
DENDRO_23 *= DENDRO_3;
DENDRO_3 = grad2_0_2_beta0;
DENDRO_23 *= DENDRO_3;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_18;
DENDRO_3 *= DENDRO_5;
DENDRO_24 = grad2_1_1_beta1;
DENDRO_3 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_18;
DENDRO_24 *= DENDRO_5;
DENDRO_5 = grad2_0_1_beta0;
DENDRO_24 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_18;
DENDRO_5 *= DENDRO_22;
DENDRO_25 = grad2_0_1_beta1;
DENDRO_5 *= DENDRO_25;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_18;
DENDRO_25 *= DENDRO_22;
DENDRO_18 = grad2_0_0_beta0;
DENDRO_25 *= DENDRO_18;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_19;
DENDRO_18 *= DENDRO_87;
DENDRO_19 = grad2_0_1_beta2;
DENDRO_18 *= DENDRO_19;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_6;
DENDRO_19 *= DENDRO_33;
DENDRO_19 *= DENDRO_30;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_6;
DENDRO_26 = grad_2_alpha;
DENDRO_22 *= DENDRO_26;
DENDRO_22 *= DENDRO_38;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_6;
DENDRO_27 = grad_1_beta2;
DENDRO_26 *= DENDRO_27;
DENDRO_26 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_6;
DENDRO_27 = grad_1_alpha;
DENDRO_4 *= DENDRO_27;
DENDRO_4 *= DENDRO_32;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_6;
DENDRO_28 = grad_0_beta2;
DENDRO_27 *= DENDRO_28;
DENDRO_27 *= DENDRO_34;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_6;
DENDRO_29 = grad_0_alpha;
DENDRO_28 *= DENDRO_29;
DENDRO_28 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_6;
DENDRO_0 *= DENDRO_9;
DENDRO_0 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_6;
DENDRO_1 *= DENDRO_9;
DENDRO_1 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_6;
DENDRO_2 *= DENDRO_9;
DENDRO_2 *= DENDRO_7;
DENDRO_6 = beta2[pp];

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_6 = grad_2_Gt2;
DENDRO_7 *= DENDRO_6;
DENDRO_6 = beta1[pp];

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_6;
DENDRO_6 = grad_1_Gt2;
DENDRO_9 *= DENDRO_6;
DENDRO_6 = beta0[pp];

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_6;
DENDRO_6 = grad_0_Gt2;
DENDRO_29 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_83;
DENDRO_30 = grad2_1_1_beta2;
DENDRO_6 *= DENDRO_30;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_97;
DENDRO_31 = grad2_0_0_beta2;
DENDRO_30 *= DENDRO_31;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_30;
DENDRO_31 += DENDRO_6;
DENDRO_31 += DENDRO_29;
DENDRO_31 += DENDRO_9;
DENDRO_31 += DENDRO_7;
DENDRO_31 += DENDRO_2;
DENDRO_31 += DENDRO_1;
DENDRO_31 += DENDRO_0;
DENDRO_31 += DENDRO_28;
DENDRO_31 += DENDRO_27;
DENDRO_31 += DENDRO_4;
DENDRO_31 += DENDRO_26;
DENDRO_31 += DENDRO_22;
DENDRO_31 += DENDRO_19;
DENDRO_31 += DENDRO_18;
DENDRO_31 += DENDRO_25;
DENDRO_31 += DENDRO_5;
DENDRO_31 += DENDRO_24;
DENDRO_31 += DENDRO_3;
DENDRO_31 += DENDRO_23;
DENDRO_31 += DENDRO_21;
DENDRO_31 += DENDRO_20;
DENDRO_31 += DENDRO_15;
DENDRO_31 += DENDRO_17;
DENDRO_31 += DENDRO_16;
DENDRO_31 += DENDRO_14;
DENDRO_31 += DENDRO_13;
DENDRO_31 += DENDRO_12;
DENDRO_31 += DENDRO_11;
DENDRO_31 += DENDRO_10;
DENDRO_31 += DENDRO_8;
Gt_rhs2[pp] = DENDRO_31;


}
if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&Gt_rhs0[gidx], Gt0[gidx], grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&Gt_rhs1[gidx], Gt1[gidx], grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&Gt_rhs2[gidx], Gt2[gidx], grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk);
}

Gt_rhs0[pp]  += ko_sigma * (kograd_0_Gt0 + kograd_1_Gt0 + kograd_2_Gt0);
Gt_rhs1[pp]  += ko_sigma * (kograd_0_Gt1 + kograd_1_Gt1 + kograd_2_Gt1);
Gt_rhs2[pp]  += ko_sigma * (kograd_0_Gt2 + kograd_1_Gt2 + kograd_2_Gt2);
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , B0, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_B0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_B0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_B0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_B0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_B0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_B0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_B0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_B0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_B0=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , B0, blk);
device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_B0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_B0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_B0=Du[gidx];
__syncthreads();

//|V|= 27 |E| =30
//topological generations = 4
// [0]
// traversal id=0 temp_var_count=7
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;
DENDRO_0 = beta2[pp];

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = grad_2_Gt0;
DENDRO_1 *= DENDRO_2;
DENDRO_2 = beta1[pp];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = grad_1_Gt0;
DENDRO_3 *= DENDRO_4;
DENDRO_4 = beta0[pp];

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = grad_0_Gt0;
DENDRO_5 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_0;
DENDRO_0 = grad_2_B0;
DENDRO_6 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_2;
DENDRO_2 = grad_1_B0;
DENDRO_0 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_4;
DENDRO_4 = grad_0_B0;
DENDRO_2 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_5;
DENDRO_4 += DENDRO_3;
DENDRO_4 += DENDRO_1;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_2;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_6;
DENDRO_0 = -1;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_3 = lambda[3];
DENDRO_2 *= DENDRO_3;
DENDRO_2 *= DENDRO_4;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = B0[pp];
DENDRO_3 *= DENDRO_0;
DENDRO_0 = eta;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = lambda[2];

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_0;
DENDRO_4 *= DENDRO_1;
DENDRO_0 = Gt_rhs0[pp];

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_4;
DENDRO_1 += DENDRO_3;
DENDRO_1 += DENDRO_2;
B_rhs0[pp] = DENDRO_1;


}
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , B1, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_B1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_B1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_B1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_B1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_B1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_B1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_B1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_B1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_B1=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , B1, blk);
device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_B1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_B1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_B1=Du[gidx];
__syncthreads();

//|V|= 27 |E| =30
//topological generations = 4
// [0]
// traversal id=0 temp_var_count=7
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;
DENDRO_0 = beta2[pp];

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = grad_2_Gt1;
DENDRO_1 *= DENDRO_2;
DENDRO_2 = beta1[pp];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = grad_1_Gt1;
DENDRO_3 *= DENDRO_4;
DENDRO_4 = beta0[pp];

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = grad_0_Gt1;
DENDRO_5 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_0;
DENDRO_0 = grad_2_B1;
DENDRO_6 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_2;
DENDRO_2 = grad_1_B1;
DENDRO_0 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_4;
DENDRO_4 = grad_0_B1;
DENDRO_2 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_5;
DENDRO_4 += DENDRO_3;
DENDRO_4 += DENDRO_1;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_2;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_6;
DENDRO_0 = -1;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_3 = lambda[3];
DENDRO_2 *= DENDRO_3;
DENDRO_2 *= DENDRO_4;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = B1[pp];
DENDRO_3 *= DENDRO_0;
DENDRO_0 = eta;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = lambda[2];

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_0;
DENDRO_4 *= DENDRO_1;
DENDRO_0 = Gt_rhs1[pp];

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_4;
DENDRO_1 += DENDRO_3;
DENDRO_1 += DENDRO_2;
B_rhs1[pp] = DENDRO_1;


}
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , B2, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_B2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_B2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_B2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_B2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_B2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_B2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_B2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_B2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_B2=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , B2, blk);
device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_B2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_B2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_B2=Du[gidx];
__syncthreads();

//|V|= 27 |E| =30
//topological generations = 4
// [0]
// traversal id=0 temp_var_count=7
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;
DENDRO_0 = beta2[pp];

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = grad_2_Gt2;
DENDRO_1 *= DENDRO_2;
DENDRO_2 = beta1[pp];

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = grad_1_Gt2;
DENDRO_3 *= DENDRO_4;
DENDRO_4 = beta0[pp];

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = grad_0_Gt2;
DENDRO_5 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_0;
DENDRO_0 = grad_2_B2;
DENDRO_6 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_2;
DENDRO_2 = grad_1_B2;
DENDRO_0 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_4;
DENDRO_4 = grad_0_B2;
DENDRO_2 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_5;
DENDRO_4 += DENDRO_3;
DENDRO_4 += DENDRO_1;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_2;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_6;
DENDRO_0 = -1;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_3 = lambda[3];
DENDRO_2 *= DENDRO_3;
DENDRO_2 *= DENDRO_4;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = B2[pp];
DENDRO_3 *= DENDRO_0;
DENDRO_0 = eta;
DENDRO_3 *= DENDRO_0;
DENDRO_0 = lambda[2];

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_0;
DENDRO_4 *= DENDRO_1;
DENDRO_0 = Gt_rhs2[pp];

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_4;
DENDRO_1 += DENDRO_3;
DENDRO_1 += DENDRO_2;
B_rhs2[pp] = DENDRO_1;


}
if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&B_rhs0[gidx], B0[gidx], grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&B_rhs1[gidx], B1[gidx], grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&B_rhs2[gidx], B2[gidx], grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk);
}

B_rhs0[pp]   += ko_sigma * (kograd_0_B0 + kograd_1_B0 + kograd_2_B0);
B_rhs1[pp]   += ko_sigma * (kograd_0_B1 + kograd_1_B1 + kograd_2_B1);
B_rhs2[pp]   += ko_sigma * (kograd_0_B2 + kograd_1_B2 + kograd_2_B2);
__syncthreads();

double DENDRO_RIJ0;
double DENDRO_RIJ1;
double DENDRO_RIJ2;
double DENDRO_RIJ3;
double DENDRO_RIJ4;
double DENDRO_RIJ5;
//|V|= 1303 |E| =2993
//topological generations = 12
// [0]
// traversal id=0 temp_var_count=184
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;double DENDRO_33;double DENDRO_34;double DENDRO_35;double DENDRO_36;double DENDRO_37;double DENDRO_38;double DENDRO_39;double DENDRO_40;double DENDRO_41;double DENDRO_42;double DENDRO_43;double DENDRO_44;double DENDRO_45;double DENDRO_46;double DENDRO_47;double DENDRO_48;double DENDRO_49;double DENDRO_50;double DENDRO_51;double DENDRO_52;double DENDRO_53;double DENDRO_54;double DENDRO_55;double DENDRO_56;double DENDRO_57;double DENDRO_58;double DENDRO_59;double DENDRO_60;double DENDRO_61;double DENDRO_62;double DENDRO_63;double DENDRO_64;double DENDRO_65;double DENDRO_66;double DENDRO_67;double DENDRO_68;double DENDRO_69;double DENDRO_70;double DENDRO_71;double DENDRO_72;double DENDRO_73;double DENDRO_74;double DENDRO_75;double DENDRO_76;double DENDRO_77;double DENDRO_78;double DENDRO_79;double DENDRO_80;double DENDRO_81;double DENDRO_82;double DENDRO_83;double DENDRO_84;double DENDRO_85;double DENDRO_86;double DENDRO_87;double DENDRO_88;double DENDRO_89;double DENDRO_90;double DENDRO_91;double DENDRO_92;double DENDRO_93;double DENDRO_94;double DENDRO_95;double DENDRO_96;double DENDRO_97;double DENDRO_98;double DENDRO_99;double DENDRO_100;double DENDRO_101;double DENDRO_102;double DENDRO_103;double DENDRO_104;double DENDRO_105;double DENDRO_106;double DENDRO_107;double DENDRO_108;double DENDRO_109;double DENDRO_110;double DENDRO_111;double DENDRO_112;double DENDRO_113;double DENDRO_114;double DENDRO_115;double DENDRO_116;double DENDRO_117;double DENDRO_118;double DENDRO_119;double DENDRO_120;double DENDRO_121;double DENDRO_122;double DENDRO_123;double DENDRO_124;double DENDRO_125;double DENDRO_126;double DENDRO_127;double DENDRO_128;double DENDRO_129;double DENDRO_130;double DENDRO_131;double DENDRO_132;double DENDRO_133;double DENDRO_134;double DENDRO_135;double DENDRO_136;double DENDRO_137;double DENDRO_138;double DENDRO_139;double DENDRO_140;double DENDRO_141;double DENDRO_142;double DENDRO_143;double DENDRO_144;double DENDRO_145;double DENDRO_146;double DENDRO_147;double DENDRO_148;double DENDRO_149;double DENDRO_150;double DENDRO_151;double DENDRO_152;double DENDRO_153;double DENDRO_154;double DENDRO_155;double DENDRO_156;double DENDRO_157;double DENDRO_158;double DENDRO_159;double DENDRO_160;double DENDRO_161;double DENDRO_162;double DENDRO_163;double DENDRO_164;double DENDRO_165;double DENDRO_166;double DENDRO_167;double DENDRO_168;double DENDRO_169;double DENDRO_170;double DENDRO_171;double DENDRO_172;double DENDRO_173;double DENDRO_174;double DENDRO_175;double DENDRO_176;double DENDRO_177;double DENDRO_178;double DENDRO_179;double DENDRO_180;double DENDRO_181;double DENDRO_182;double DENDRO_183;
DENDRO_0 = gt1[pp];

// int-power reduce pow(gt1[pp],2)
DENDRO_1 = DENDRO_0;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = gt2[pp];

// int-power reduce pow(gt2[pp],2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = gt4[pp];

// int-power reduce pow(gt4[pp],2)
DENDRO_5 = DENDRO_4;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = -2;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_4;
DENDRO_8 = -1;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_8;
DENDRO_10 = gt0[pp];
DENDRO_9 *= DENDRO_10;
DENDRO_11 = gt3[pp];
DENDRO_9 *= DENDRO_11;
DENDRO_12 = gt5[pp];
DENDRO_9 *= DENDRO_12;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_12;
DENDRO_13 *= DENDRO_1;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_11;
DENDRO_14 *= DENDRO_3;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_10;
DENDRO_15 *= DENDRO_5;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_8;
DENDRO_16 *= DENDRO_0;
DENDRO_16 *= DENDRO_4;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_2;
DENDRO_17 *= DENDRO_11;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_15;
DENDRO_18 += DENDRO_14;
DENDRO_18 += DENDRO_13;
DENDRO_18 += DENDRO_9;
DENDRO_18 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_8;
DENDRO_7 *= DENDRO_10;
DENDRO_7 *= DENDRO_11;
DENDRO_9 = -0.5;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_9;
DENDRO_9 = grad_1_gt2;
DENDRO_13 *= DENDRO_9;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_14;
DENDRO_14 = grad_2_gt1;
DENDRO_15 *= DENDRO_14;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_19;
DENDRO_19 = grad_0_gt4;
DENDRO_20 *= DENDRO_19;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_8;
DENDRO_21 *= DENDRO_0;
DENDRO_21 *= DENDRO_2;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_10;
DENDRO_22 *= DENDRO_4;
DENDRO_23 = -0.5;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_23;
DENDRO_24 *= DENDRO_14;
DENDRO_23 = 0.5;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_23;
DENDRO_25 *= DENDRO_9;
DENDRO_23 = -0.5;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_23;
DENDRO_26 *= DENDRO_19;
DENDRO_23 = -0.5;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_23;
DENDRO_23 = grad_1_gt0;
DENDRO_27 *= DENDRO_23;
DENDRO_28 = 1.0;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_28;
DENDRO_28 = grad_0_gt1;
DENDRO_29 *= DENDRO_28;
DENDRO_30 = -0.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_30;
DENDRO_30 = grad_2_gt0;
DENDRO_31 *= DENDRO_30;
DENDRO_32 = 1.0;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_32;
DENDRO_32 = grad_0_gt2;
DENDRO_33 *= DENDRO_32;
DENDRO_34 = -0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_34;
DENDRO_34 = grad_0_gt3;
DENDRO_35 *= DENDRO_34;
DENDRO_36 = 1.0;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_36;
DENDRO_36 = grad_1_gt1;
DENDRO_37 *= DENDRO_36;
DENDRO_38 = -0.5;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_38;
DENDRO_38 = grad_2_gt3;
DENDRO_39 *= DENDRO_38;
DENDRO_40 = 1.0;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_40;
DENDRO_40 = grad_1_gt4;
DENDRO_41 *= DENDRO_40;
DENDRO_42 = -0.5;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_42;
DENDRO_42 = grad_1_gt5;
DENDRO_43 *= DENDRO_42;
DENDRO_44 = 1.0;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_44;
DENDRO_44 = grad_2_gt4;
DENDRO_45 *= DENDRO_44;
DENDRO_46 = -0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_46;
DENDRO_46 = grad_0_gt5;
DENDRO_47 *= DENDRO_46;
DENDRO_48 = 1.0;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_48;
DENDRO_48 = grad_2_gt2;
DENDRO_49 *= DENDRO_48;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_8;
DENDRO_50 *= DENDRO_2;
DENDRO_50 *= DENDRO_4;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_0;
DENDRO_51 *= DENDRO_12;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_8;
DENDRO_52 *= DENDRO_10;
DENDRO_52 *= DENDRO_12;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_8;
DENDRO_53 *= DENDRO_11;
DENDRO_53 *= DENDRO_12;

// initialize reduction 
DENDRO_54 = 0;
DENDRO_54 += DENDRO_17;
DENDRO_54 += DENDRO_16;

// int-power reduce pow(-gt0[pp]*gt3[pp]*gt5[pp] + gt0[pp]*gt4[pp]**2 + gt1[pp]**2*gt5[pp] - 2*gt1[pp]*gt2[pp]*gt4[pp] + gt2[pp]**2*gt3[pp],-1)
DENDRO_16 = DENDRO_18;
DENDRO_16 = 1.0/DENDRO_16;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_1;
DENDRO_17 += DENDRO_7;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_20;
DENDRO_1 += DENDRO_15;
DENDRO_1 += DENDRO_13;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_22;
DENDRO_7 += DENDRO_21;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_20;
DENDRO_13 += DENDRO_25;
DENDRO_13 += DENDRO_24;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_25;
DENDRO_18 += DENDRO_15;
DENDRO_18 += DENDRO_26;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_29;
DENDRO_15 += DENDRO_27;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_33;
DENDRO_20 += DENDRO_31;

// initialize reduction 
DENDRO_21 = 0;
DENDRO_21 += DENDRO_37;
DENDRO_21 += DENDRO_35;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_41;
DENDRO_22 += DENDRO_39;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_45;
DENDRO_24 += DENDRO_43;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_49;
DENDRO_25 += DENDRO_47;

// initialize reduction 
DENDRO_26 = 0;
DENDRO_26 += DENDRO_51;
DENDRO_26 += DENDRO_50;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_3;
DENDRO_27 += DENDRO_52;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_5;
DENDRO_3 += DENDRO_53;
DENDRO_5 = 0.5;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_5;
DENDRO_29 *= DENDRO_30;
DENDRO_29 *= DENDRO_16;
DENDRO_29 *= DENDRO_54;
DENDRO_5 = 0.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_5;
DENDRO_31 *= DENDRO_46;
DENDRO_31 *= DENDRO_16;
DENDRO_31 *= DENDRO_17;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_16;
DENDRO_5 *= DENDRO_7;
DENDRO_5 *= DENDRO_1;
DENDRO_33 = 0.5;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_33;
DENDRO_35 *= DENDRO_23;
DENDRO_35 *= DENDRO_16;
DENDRO_35 *= DENDRO_54;
DENDRO_33 = 0.5;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_33;
DENDRO_37 *= DENDRO_34;
DENDRO_37 *= DENDRO_16;
DENDRO_37 *= DENDRO_7;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_16;
DENDRO_33 *= DENDRO_17;
DENDRO_33 *= DENDRO_13;
DENDRO_39 = 0.5;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_39;
DENDRO_41 *= DENDRO_38;
DENDRO_41 *= DENDRO_16;
DENDRO_41 *= DENDRO_7;
DENDRO_39 = 0.5;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_39;
DENDRO_43 *= DENDRO_42;
DENDRO_43 *= DENDRO_16;
DENDRO_43 *= DENDRO_17;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_16;
DENDRO_39 *= DENDRO_54;
DENDRO_39 *= DENDRO_18;
DENDRO_45 = 0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_45;
DENDRO_45 = grad_0_gt0;
DENDRO_47 *= DENDRO_45;
DENDRO_47 *= DENDRO_16;
DENDRO_47 *= DENDRO_54;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_16;
DENDRO_49 *= DENDRO_15;
DENDRO_49 *= DENDRO_7;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_16;
DENDRO_50 *= DENDRO_17;
DENDRO_50 *= DENDRO_20;
DENDRO_51 = 0.5;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_51;
DENDRO_51 = grad_1_gt3;
DENDRO_52 *= DENDRO_51;
DENDRO_52 *= DENDRO_16;
DENDRO_52 *= DENDRO_7;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_16;
DENDRO_53 *= DENDRO_21;
DENDRO_53 *= DENDRO_54;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_16;
DENDRO_55 *= DENDRO_17;
DENDRO_55 *= DENDRO_22;
DENDRO_56 = 0.5;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_56;
DENDRO_56 = grad_2_gt5;
DENDRO_57 *= DENDRO_56;
DENDRO_57 *= DENDRO_16;
DENDRO_57 *= DENDRO_17;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_16;
DENDRO_58 *= DENDRO_24;
DENDRO_58 *= DENDRO_7;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_16;
DENDRO_59 *= DENDRO_25;
DENDRO_59 *= DENDRO_54;
DENDRO_60 = 0.5;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_60;
DENDRO_61 *= DENDRO_30;
DENDRO_61 *= DENDRO_16;
DENDRO_61 *= DENDRO_26;
DENDRO_60 = 0.5;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_60;
DENDRO_62 *= DENDRO_46;
DENDRO_62 *= DENDRO_16;
DENDRO_62 *= DENDRO_7;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_16;
DENDRO_60 *= DENDRO_27;
DENDRO_60 *= DENDRO_1;
DENDRO_63 = 0.5;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_63;
DENDRO_64 *= DENDRO_23;
DENDRO_64 *= DENDRO_16;
DENDRO_64 *= DENDRO_26;
DENDRO_63 = 0.5;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_63;
DENDRO_65 *= DENDRO_34;
DENDRO_65 *= DENDRO_16;
DENDRO_65 *= DENDRO_27;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_16;
DENDRO_63 *= DENDRO_7;
DENDRO_63 *= DENDRO_13;
DENDRO_66 = 0.5;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_66;
DENDRO_67 *= DENDRO_38;
DENDRO_67 *= DENDRO_16;
DENDRO_67 *= DENDRO_27;
DENDRO_66 = 0.5;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_66;
DENDRO_68 *= DENDRO_42;
DENDRO_68 *= DENDRO_16;
DENDRO_68 *= DENDRO_7;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_16;
DENDRO_66 *= DENDRO_26;
DENDRO_66 *= DENDRO_18;
DENDRO_69 = 0.5;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_69;
DENDRO_70 *= DENDRO_45;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_26;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_20;
DENDRO_69 *= DENDRO_7;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_16;
DENDRO_71 *= DENDRO_27;
DENDRO_71 *= DENDRO_15;
DENDRO_72 = 0.5;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_72;
DENDRO_73 *= DENDRO_51;
DENDRO_73 *= DENDRO_16;
DENDRO_73 *= DENDRO_27;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_22;
DENDRO_72 *= DENDRO_7;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_16;
DENDRO_74 *= DENDRO_21;
DENDRO_74 *= DENDRO_26;
DENDRO_75 = 0.5;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_75;
DENDRO_76 *= DENDRO_56;
DENDRO_76 *= DENDRO_16;
DENDRO_76 *= DENDRO_7;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_16;
DENDRO_75 *= DENDRO_25;
DENDRO_75 *= DENDRO_26;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_16;
DENDRO_77 *= DENDRO_27;
DENDRO_77 *= DENDRO_24;
DENDRO_78 = 0.5;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_78;
DENDRO_79 *= DENDRO_30;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_3;
DENDRO_78 = 0.5;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_78;
DENDRO_80 *= DENDRO_46;
DENDRO_80 *= DENDRO_16;
DENDRO_80 *= DENDRO_54;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_16;
DENDRO_78 *= DENDRO_26;
DENDRO_78 *= DENDRO_1;
DENDRO_81 = 0.5;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_81;
DENDRO_82 *= DENDRO_23;
DENDRO_82 *= DENDRO_16;
DENDRO_82 *= DENDRO_3;
DENDRO_81 = 0.5;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_81;
DENDRO_83 *= DENDRO_34;
DENDRO_83 *= DENDRO_16;
DENDRO_83 *= DENDRO_26;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_16;
DENDRO_81 *= DENDRO_54;
DENDRO_81 *= DENDRO_13;
DENDRO_84 = 0.5;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_84;
DENDRO_85 *= DENDRO_38;
DENDRO_85 *= DENDRO_16;
DENDRO_85 *= DENDRO_26;
DENDRO_84 = 0.5;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_84;
DENDRO_86 *= DENDRO_42;
DENDRO_86 *= DENDRO_16;
DENDRO_86 *= DENDRO_54;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_16;
DENDRO_84 *= DENDRO_3;
DENDRO_84 *= DENDRO_18;
DENDRO_87 = 0.5;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_87;
DENDRO_88 *= DENDRO_45;
DENDRO_88 *= DENDRO_16;
DENDRO_88 *= DENDRO_3;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_16;
DENDRO_87 *= DENDRO_20;
DENDRO_87 *= DENDRO_54;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_16;
DENDRO_89 *= DENDRO_15;
DENDRO_89 *= DENDRO_26;
DENDRO_90 = 0.5;

// initialize reduction 
DENDRO_91 = 1;
DENDRO_91 *= DENDRO_90;
DENDRO_91 *= DENDRO_51;
DENDRO_91 *= DENDRO_16;
DENDRO_91 *= DENDRO_26;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_16;
DENDRO_90 *= DENDRO_22;
DENDRO_90 *= DENDRO_54;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_16;
DENDRO_92 *= DENDRO_3;
DENDRO_92 *= DENDRO_21;
DENDRO_93 = 0.5;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_93;
DENDRO_94 *= DENDRO_56;
DENDRO_94 *= DENDRO_16;
DENDRO_94 *= DENDRO_54;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_16;
DENDRO_93 *= DENDRO_24;
DENDRO_93 *= DENDRO_26;

// initialize reduction 
DENDRO_95 = 1;
DENDRO_95 *= DENDRO_16;
DENDRO_95 *= DENDRO_3;
DENDRO_95 *= DENDRO_25;
DENDRO_96 = chi[pp];

// int-power reduce pow(chi[pp],-1)
DENDRO_97 = DENDRO_96;
DENDRO_97 = 1.0/DENDRO_97;
DENDRO_98 = grad_2_chi;

// int-power reduce pow(grad_2_chi,2)
DENDRO_99 = DENDRO_98;
DENDRO_99 *= DENDRO_98;
DENDRO_100 = grad_1_chi;

// int-power reduce pow(grad_1_chi,2)
DENDRO_101 = DENDRO_100;
DENDRO_101 *= DENDRO_100;
DENDRO_102 = grad_0_chi;

// int-power reduce pow(grad_0_chi,2)
DENDRO_103 = DENDRO_102;
DENDRO_103 *= DENDRO_102;

// initialize reduction 
DENDRO_104 = 0;
DENDRO_104 += DENDRO_5;
DENDRO_104 += DENDRO_31;
DENDRO_104 += DENDRO_29;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_33;
DENDRO_5 += DENDRO_37;
DENDRO_5 += DENDRO_35;

// initialize reduction 
DENDRO_29 = 0;
DENDRO_29 += DENDRO_39;
DENDRO_29 += DENDRO_43;
DENDRO_29 += DENDRO_41;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_50;
DENDRO_31 += DENDRO_49;
DENDRO_31 += DENDRO_47;

// initialize reduction 
DENDRO_33 = 0;
DENDRO_33 += DENDRO_55;
DENDRO_33 += DENDRO_53;
DENDRO_33 += DENDRO_52;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_59;
DENDRO_35 += DENDRO_58;
DENDRO_35 += DENDRO_57;

// initialize reduction 
DENDRO_37 = 0;
DENDRO_37 += DENDRO_60;
DENDRO_37 += DENDRO_62;
DENDRO_37 += DENDRO_61;

// initialize reduction 
DENDRO_39 = 0;
DENDRO_39 += DENDRO_63;
DENDRO_39 += DENDRO_65;
DENDRO_39 += DENDRO_64;

// initialize reduction 
DENDRO_41 = 0;
DENDRO_41 += DENDRO_66;
DENDRO_41 += DENDRO_68;
DENDRO_41 += DENDRO_67;

// initialize reduction 
DENDRO_43 = 0;
DENDRO_43 += DENDRO_71;
DENDRO_43 += DENDRO_69;
DENDRO_43 += DENDRO_70;

// initialize reduction 
DENDRO_47 = 0;
DENDRO_47 += DENDRO_74;
DENDRO_47 += DENDRO_72;
DENDRO_47 += DENDRO_73;

// initialize reduction 
DENDRO_49 = 0;
DENDRO_49 += DENDRO_77;
DENDRO_49 += DENDRO_75;
DENDRO_49 += DENDRO_76;

// initialize reduction 
DENDRO_50 = 0;
DENDRO_50 += DENDRO_78;
DENDRO_50 += DENDRO_80;
DENDRO_50 += DENDRO_79;

// initialize reduction 
DENDRO_52 = 0;
DENDRO_52 += DENDRO_81;
DENDRO_52 += DENDRO_83;
DENDRO_52 += DENDRO_82;

// initialize reduction 
DENDRO_53 = 0;
DENDRO_53 += DENDRO_84;
DENDRO_53 += DENDRO_86;
DENDRO_53 += DENDRO_85;

// initialize reduction 
DENDRO_55 = 0;
DENDRO_55 += DENDRO_89;
DENDRO_55 += DENDRO_87;
DENDRO_55 += DENDRO_88;

// initialize reduction 
DENDRO_57 = 0;
DENDRO_57 += DENDRO_92;
DENDRO_57 += DENDRO_90;
DENDRO_57 += DENDRO_91;

// initialize reduction 
DENDRO_58 = 0;
DENDRO_58 += DENDRO_95;
DENDRO_58 += DENDRO_93;
DENDRO_58 += DENDRO_94;
DENDRO_59 = -3.0/2.0;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_59;
DENDRO_60 *= DENDRO_100;
DENDRO_60 *= DENDRO_98;
DENDRO_60 *= DENDRO_97;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_59;
DENDRO_61 *= DENDRO_102;
DENDRO_61 *= DENDRO_98;
DENDRO_61 *= DENDRO_97;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_59;
DENDRO_62 *= DENDRO_102;
DENDRO_62 *= DENDRO_100;
DENDRO_62 *= DENDRO_97;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_59;
DENDRO_63 *= DENDRO_97;
DENDRO_63 *= DENDRO_99;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_59;
DENDRO_64 *= DENDRO_97;
DENDRO_64 *= DENDRO_101;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_59;
DENDRO_65 *= DENDRO_97;
DENDRO_65 *= DENDRO_103;
DENDRO_59 = 2;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_59;
DENDRO_66 *= DENDRO_16;
DENDRO_66 *= DENDRO_54;
DENDRO_66 *= DENDRO_104;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_59;
DENDRO_67 *= DENDRO_16;
DENDRO_67 *= DENDRO_26;
DENDRO_67 *= DENDRO_5;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_59;
DENDRO_68 *= DENDRO_16;
DENDRO_68 *= DENDRO_7;
DENDRO_68 *= DENDRO_29;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_3;
DENDRO_69 *= DENDRO_31;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_27;
DENDRO_70 *= DENDRO_33;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_16;
DENDRO_71 *= DENDRO_17;
DENDRO_71 *= DENDRO_35;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_59;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_54;
DENDRO_72 *= DENDRO_37;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_59;
DENDRO_73 *= DENDRO_16;
DENDRO_73 *= DENDRO_26;
DENDRO_73 *= DENDRO_39;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_59;
DENDRO_74 *= DENDRO_16;
DENDRO_74 *= DENDRO_7;
DENDRO_74 *= DENDRO_41;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_16;
DENDRO_75 *= DENDRO_3;
DENDRO_75 *= DENDRO_43;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_16;
DENDRO_76 *= DENDRO_27;
DENDRO_76 *= DENDRO_47;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_16;
DENDRO_77 *= DENDRO_17;
DENDRO_77 *= DENDRO_49;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_59;
DENDRO_78 *= DENDRO_16;
DENDRO_78 *= DENDRO_54;
DENDRO_78 *= DENDRO_50;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_59;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_26;
DENDRO_79 *= DENDRO_52;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_59;
DENDRO_80 *= DENDRO_16;
DENDRO_80 *= DENDRO_7;
DENDRO_80 *= DENDRO_53;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_16;
DENDRO_81 *= DENDRO_3;
DENDRO_81 *= DENDRO_55;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_16;
DENDRO_82 *= DENDRO_27;
DENDRO_82 *= DENDRO_57;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_16;
DENDRO_83 *= DENDRO_17;
DENDRO_83 *= DENDRO_58;
DENDRO_84 = grad2_1_2_chi;

// initialize reduction 
DENDRO_85 = 0;
DENDRO_85 += DENDRO_84;
DENDRO_85 += DENDRO_60;
DENDRO_60 = grad2_0_2_chi;

// initialize reduction 
DENDRO_86 = 0;
DENDRO_86 += DENDRO_60;
DENDRO_86 += DENDRO_61;
DENDRO_61 = grad2_0_1_chi;

// initialize reduction 
DENDRO_87 = 0;
DENDRO_87 += DENDRO_61;
DENDRO_87 += DENDRO_62;
DENDRO_62 = grad2_2_2_chi;

// initialize reduction 
DENDRO_88 = 0;
DENDRO_88 += DENDRO_62;
DENDRO_88 += DENDRO_63;
DENDRO_63 = grad2_1_1_chi;

// initialize reduction 
DENDRO_89 = 0;
DENDRO_89 += DENDRO_63;
DENDRO_89 += DENDRO_64;
DENDRO_64 = grad2_0_0_chi;

// initialize reduction 
DENDRO_90 = 0;
DENDRO_90 += DENDRO_64;
DENDRO_90 += DENDRO_65;

// initialize reduction 
DENDRO_65 = 0;
DENDRO_65 += DENDRO_71;
DENDRO_65 += DENDRO_70;
DENDRO_65 += DENDRO_69;
DENDRO_65 += DENDRO_68;
DENDRO_65 += DENDRO_67;
DENDRO_65 += DENDRO_66;

// initialize reduction 
DENDRO_66 = 0;
DENDRO_66 += DENDRO_77;
DENDRO_66 += DENDRO_76;
DENDRO_66 += DENDRO_75;
DENDRO_66 += DENDRO_74;
DENDRO_66 += DENDRO_73;
DENDRO_66 += DENDRO_72;

// initialize reduction 
DENDRO_67 = 0;
DENDRO_67 += DENDRO_83;
DENDRO_67 += DENDRO_82;
DENDRO_67 += DENDRO_81;
DENDRO_67 += DENDRO_80;
DENDRO_67 += DENDRO_79;
DENDRO_67 += DENDRO_78;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_59;
DENDRO_68 *= DENDRO_16;
DENDRO_68 *= DENDRO_85;
DENDRO_68 *= DENDRO_7;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_59;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_86;
DENDRO_69 *= DENDRO_54;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_59;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_87;
DENDRO_70 *= DENDRO_26;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_16;
DENDRO_71 *= DENDRO_88;
DENDRO_71 *= DENDRO_17;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_89;
DENDRO_72 *= DENDRO_27;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_16;
DENDRO_73 *= DENDRO_90;
DENDRO_73 *= DENDRO_3;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_8;
DENDRO_74 *= DENDRO_98;
DENDRO_74 *= DENDRO_65;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_8;
DENDRO_75 *= DENDRO_100;
DENDRO_75 *= DENDRO_66;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_8;
DENDRO_76 *= DENDRO_102;
DENDRO_76 *= DENDRO_67;
DENDRO_77 = 1.0;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_77;
DENDRO_78 *= DENDRO_46;
DENDRO_78 *= DENDRO_35;
DENDRO_77 = 0.5;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_77;
DENDRO_79 *= DENDRO_56;
DENDRO_79 *= DENDRO_104;
DENDRO_77 = 1.0;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_77;
DENDRO_80 *= DENDRO_46;
DENDRO_80 *= DENDRO_50;
DENDRO_77 = 0.5;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_77;
DENDRO_81 *= DENDRO_30;
DENDRO_81 *= DENDRO_58;
DENDRO_82 = 1.0;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_82;
DENDRO_83 *= DENDRO_56;
DENDRO_83 *= DENDRO_104;
DENDRO_82 = 0.5;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_82;
DENDRO_85 *= DENDRO_46;
DENDRO_85 *= DENDRO_35;
DENDRO_82 = 1.0;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_82;
DENDRO_86 *= DENDRO_42;
DENDRO_86 *= DENDRO_37;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_1;
DENDRO_82 *= DENDRO_49;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_59;
DENDRO_87 *= DENDRO_13;
DENDRO_87 *= DENDRO_49;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_24;
DENDRO_88 *= DENDRO_37;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_59;
DENDRO_89 *= DENDRO_20;
DENDRO_89 *= DENDRO_58;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_25;
DENDRO_90 *= DENDRO_50;
DENDRO_91 = 1.0;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_91;
DENDRO_92 *= DENDRO_46;
DENDRO_92 *= DENDRO_29;
DENDRO_91 = 0.5;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_91;
DENDRO_93 *= DENDRO_42;
DENDRO_93 *= DENDRO_104;
DENDRO_94 = 1.0;

// initialize reduction 
DENDRO_95 = 1;
DENDRO_95 *= DENDRO_94;
DENDRO_95 *= DENDRO_42;
DENDRO_95 *= DENDRO_104;
DENDRO_94 = 0.5;

// initialize reduction 
DENDRO_105 = 1;
DENDRO_105 *= DENDRO_94;
DENDRO_105 *= DENDRO_46;
DENDRO_105 *= DENDRO_29;
DENDRO_106 = 0.5;

// initialize reduction 
DENDRO_107 = 1;
DENDRO_107 *= DENDRO_106;
DENDRO_107 *= DENDRO_38;
DENDRO_107 *= DENDRO_37;

// initialize reduction 
DENDRO_108 = 1;
DENDRO_108 *= DENDRO_59;
DENDRO_108 *= DENDRO_13;
DENDRO_108 *= DENDRO_41;
DENDRO_109 = 0.5;

// initialize reduction 
DENDRO_110 = 1;
DENDRO_110 *= DENDRO_109;
DENDRO_110 *= DENDRO_30;
DENDRO_110 *= DENDRO_53;

// initialize reduction 
DENDRO_111 = 1;
DENDRO_111 *= DENDRO_59;
DENDRO_111 *= DENDRO_13;
DENDRO_111 *= DENDRO_50;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_59;
DENDRO_112 *= DENDRO_20;
DENDRO_112 *= DENDRO_53;

// initialize reduction 
DENDRO_113 = 1;
DENDRO_113 *= DENDRO_18;
DENDRO_113 *= DENDRO_50;

// initialize reduction 
DENDRO_114 = 1;
DENDRO_114 *= DENDRO_59;
DENDRO_114 *= DENDRO_22;
DENDRO_114 *= DENDRO_37;

// initialize reduction 
DENDRO_115 = 1;
DENDRO_115 *= DENDRO_1;
DENDRO_115 *= DENDRO_41;
DENDRO_116 = 1.0;

// initialize reduction 
DENDRO_117 = 1;
DENDRO_117 *= DENDRO_116;
DENDRO_117 *= DENDRO_42;
DENDRO_117 *= DENDRO_35;
DENDRO_116 = 0.5;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_116;
DENDRO_118 *= DENDRO_56;
DENDRO_118 *= DENDRO_29;
DENDRO_116 = 1.0;

// initialize reduction 
DENDRO_119 = 1;
DENDRO_119 *= DENDRO_116;
DENDRO_119 *= DENDRO_42;
DENDRO_119 *= DENDRO_41;
DENDRO_116 = 0.5;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_116;
DENDRO_120 *= DENDRO_38;
DENDRO_120 *= DENDRO_49;
DENDRO_121 = 1.0;

// initialize reduction 
DENDRO_122 = 1;
DENDRO_122 *= DENDRO_121;
DENDRO_122 *= DENDRO_56;
DENDRO_122 *= DENDRO_29;
DENDRO_121 = 0.5;

// initialize reduction 
DENDRO_123 = 1;
DENDRO_123 *= DENDRO_121;
DENDRO_123 *= DENDRO_42;
DENDRO_123 *= DENDRO_35;
DENDRO_121 = 1.0;

// initialize reduction 
DENDRO_124 = 1;
DENDRO_124 *= DENDRO_121;
DENDRO_124 *= DENDRO_46;
DENDRO_124 *= DENDRO_53;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_18;
DENDRO_121 *= DENDRO_58;

// initialize reduction 
DENDRO_125 = 1;
DENDRO_125 *= DENDRO_59;
DENDRO_125 *= DENDRO_22;
DENDRO_125 *= DENDRO_49;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_24;
DENDRO_126 *= DENDRO_41;

// initialize reduction 
DENDRO_127 = 1;
DENDRO_127 *= DENDRO_59;
DENDRO_127 *= DENDRO_13;
DENDRO_127 *= DENDRO_58;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_25;
DENDRO_128 *= DENDRO_53;
DENDRO_129 = 0.5;

// initialize reduction 
DENDRO_130 = 1;
DENDRO_130 *= DENDRO_129;
DENDRO_130 *= DENDRO_30;
DENDRO_130 *= DENDRO_50;

// initialize reduction 
DENDRO_129 = 1;
DENDRO_129 *= DENDRO_59;
DENDRO_129 *= DENDRO_20;
DENDRO_129 *= DENDRO_50;

// initialize reduction 
DENDRO_131 = 1;
DENDRO_131 *= DENDRO_59;
DENDRO_131 *= DENDRO_13;
DENDRO_131 *= DENDRO_37;

// initialize reduction 
DENDRO_132 = 1;
DENDRO_132 *= DENDRO_1;
DENDRO_132 *= DENDRO_37;
DENDRO_133 = 0.5;

// initialize reduction 
DENDRO_134 = 1;
DENDRO_134 *= DENDRO_133;
DENDRO_134 *= DENDRO_38;
DENDRO_134 *= DENDRO_41;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_59;
DENDRO_133 *= DENDRO_22;
DENDRO_133 *= DENDRO_41;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_59;
DENDRO_135 *= DENDRO_13;
DENDRO_135 *= DENDRO_53;

// initialize reduction 
DENDRO_136 = 1;
DENDRO_136 *= DENDRO_18;
DENDRO_136 *= DENDRO_53;
DENDRO_137 = 1.0;

// initialize reduction 
DENDRO_138 = 1;
DENDRO_138 *= DENDRO_137;
DENDRO_138 *= DENDRO_42;
DENDRO_138 *= DENDRO_49;

// initialize reduction 
DENDRO_137 = 1;
DENDRO_137 *= DENDRO_24;
DENDRO_137 *= DENDRO_49;
DENDRO_139 = 1.0;

// initialize reduction 
DENDRO_140 = 1;
DENDRO_140 *= DENDRO_139;
DENDRO_140 *= DENDRO_46;
DENDRO_140 *= DENDRO_58;

// initialize reduction 
DENDRO_139 = 1;
DENDRO_139 *= DENDRO_25;
DENDRO_139 *= DENDRO_58;

// initialize reduction 
DENDRO_141 = 1;
DENDRO_141 *= DENDRO_8;
DENDRO_141 *= DENDRO_98;
DENDRO_141 *= DENDRO_35;

// initialize reduction 
DENDRO_142 = 1;
DENDRO_142 *= DENDRO_8;
DENDRO_142 *= DENDRO_100;
DENDRO_142 *= DENDRO_49;

// initialize reduction 
DENDRO_143 = 1;
DENDRO_143 *= DENDRO_8;
DENDRO_143 *= DENDRO_102;
DENDRO_143 *= DENDRO_58;

// initialize reduction 
DENDRO_144 = 0;
DENDRO_144 += DENDRO_76;
DENDRO_144 += DENDRO_75;
DENDRO_144 += DENDRO_74;
DENDRO_144 += DENDRO_73;
DENDRO_144 += DENDRO_72;
DENDRO_144 += DENDRO_71;
DENDRO_144 += DENDRO_70;
DENDRO_144 += DENDRO_69;
DENDRO_144 += DENDRO_68;

// initialize reduction 
DENDRO_68 = 0;
DENDRO_68 += DENDRO_79;
DENDRO_68 += DENDRO_78;

// initialize reduction 
DENDRO_69 = 0;
DENDRO_69 += DENDRO_81;
DENDRO_69 += DENDRO_80;

// initialize reduction 
DENDRO_70 = 0;
DENDRO_70 += DENDRO_85;
DENDRO_70 += DENDRO_83;

// initialize reduction 
DENDRO_71 = 0;
DENDRO_71 += DENDRO_82;
DENDRO_71 += DENDRO_86;

// initialize reduction 
DENDRO_72 = 0;
DENDRO_72 += DENDRO_88;
DENDRO_72 += DENDRO_87;

// initialize reduction 
DENDRO_73 = 0;
DENDRO_73 += DENDRO_90;
DENDRO_73 += DENDRO_89;

// initialize reduction 
DENDRO_74 = 0;
DENDRO_74 += DENDRO_93;
DENDRO_74 += DENDRO_92;

// initialize reduction 
DENDRO_75 = 0;
DENDRO_75 += DENDRO_105;
DENDRO_75 += DENDRO_95;

// initialize reduction 
DENDRO_76 = 0;
DENDRO_76 += DENDRO_108;
DENDRO_76 += DENDRO_107;

// initialize reduction 
DENDRO_78 = 0;
DENDRO_78 += DENDRO_111;
DENDRO_78 += DENDRO_110;

// initialize reduction 
DENDRO_79 = 0;
DENDRO_79 += DENDRO_113;
DENDRO_79 += DENDRO_112;

// initialize reduction 
DENDRO_80 = 0;
DENDRO_80 += DENDRO_115;
DENDRO_80 += DENDRO_114;

// initialize reduction 
DENDRO_81 = 0;
DENDRO_81 += DENDRO_118;
DENDRO_81 += DENDRO_117;

// initialize reduction 
DENDRO_82 = 0;
DENDRO_82 += DENDRO_120;
DENDRO_82 += DENDRO_119;

// initialize reduction 
DENDRO_83 = 0;
DENDRO_83 += DENDRO_123;
DENDRO_83 += DENDRO_122;

// initialize reduction 
DENDRO_85 = 0;
DENDRO_85 += DENDRO_121;
DENDRO_85 += DENDRO_124;

// initialize reduction 
DENDRO_86 = 0;
DENDRO_86 += DENDRO_126;
DENDRO_86 += DENDRO_125;

// initialize reduction 
DENDRO_87 = 0;
DENDRO_87 += DENDRO_128;
DENDRO_87 += DENDRO_127;

// initialize reduction 
DENDRO_88 = 0;
DENDRO_88 += DENDRO_129;
DENDRO_88 += DENDRO_130;

// initialize reduction 
DENDRO_89 = 0;
DENDRO_89 += DENDRO_132;
DENDRO_89 += DENDRO_131;

// initialize reduction 
DENDRO_90 = 0;
DENDRO_90 += DENDRO_133;
DENDRO_90 += DENDRO_134;

// initialize reduction 
DENDRO_92 = 0;
DENDRO_92 += DENDRO_136;
DENDRO_92 += DENDRO_135;

// initialize reduction 
DENDRO_93 = 0;
DENDRO_93 += DENDRO_137;
DENDRO_93 += DENDRO_138;

// initialize reduction 
DENDRO_95 = 0;
DENDRO_95 += DENDRO_139;
DENDRO_95 += DENDRO_140;

// int-power reduce pow(chi[pp],-2)
DENDRO_105 = DENDRO_96;
DENDRO_105 *= DENDRO_96;
DENDRO_105 = 1.0/DENDRO_105;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_62;
DENDRO_6 += DENDRO_143;
DENDRO_6 += DENDRO_142;
DENDRO_6 += DENDRO_141;
DENDRO_62 = 1.5;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_62;
DENDRO_96 *= DENDRO_56;
DENDRO_96 *= DENDRO_16;
DENDRO_96 *= DENDRO_17;
DENDRO_96 *= DENDRO_35;
DENDRO_62 = 1.5;

// initialize reduction 
DENDRO_107 = 1;
DENDRO_107 *= DENDRO_62;
DENDRO_107 *= DENDRO_42;
DENDRO_107 *= DENDRO_16;
DENDRO_107 *= DENDRO_27;
DENDRO_107 *= DENDRO_29;
DENDRO_62 = 1.5;

// initialize reduction 
DENDRO_108 = 1;
DENDRO_108 *= DENDRO_62;
DENDRO_108 *= DENDRO_46;
DENDRO_108 *= DENDRO_16;
DENDRO_108 *= DENDRO_3;
DENDRO_108 *= DENDRO_104;
DENDRO_62 = -1.0;

// initialize reduction 
DENDRO_110 = 1;
DENDRO_110 *= DENDRO_62;
DENDRO_62 = grad2_1_2_gt5;
DENDRO_110 *= DENDRO_62;
DENDRO_110 *= DENDRO_16;
DENDRO_110 *= DENDRO_7;
DENDRO_62 = -1.0;

// initialize reduction 
DENDRO_111 = 1;
DENDRO_111 *= DENDRO_62;
DENDRO_62 = grad2_0_2_gt5;
DENDRO_111 *= DENDRO_62;
DENDRO_111 *= DENDRO_16;
DENDRO_111 *= DENDRO_54;
DENDRO_62 = -1.0;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_62;
DENDRO_62 = grad2_0_1_gt5;
DENDRO_112 *= DENDRO_62;
DENDRO_112 *= DENDRO_16;
DENDRO_112 *= DENDRO_26;
DENDRO_62 = -0.5;

// initialize reduction 
DENDRO_114 = 1;
DENDRO_114 *= DENDRO_62;
DENDRO_62 = grad2_2_2_gt5;
DENDRO_114 *= DENDRO_62;
DENDRO_114 *= DENDRO_16;
DENDRO_114 *= DENDRO_17;
DENDRO_62 = -0.5;

// initialize reduction 
DENDRO_117 = 1;
DENDRO_117 *= DENDRO_62;
DENDRO_62 = grad2_1_1_gt5;
DENDRO_117 *= DENDRO_62;
DENDRO_117 *= DENDRO_16;
DENDRO_117 *= DENDRO_27;
DENDRO_62 = -0.5;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_62;
DENDRO_62 = grad2_0_0_gt5;
DENDRO_118 *= DENDRO_62;
DENDRO_118 *= DENDRO_16;
DENDRO_118 *= DENDRO_3;
DENDRO_62 = 1.0/2.0;

// initialize reduction 
DENDRO_119 = 1;
DENDRO_119 *= DENDRO_62;
DENDRO_119 *= DENDRO_12;
DENDRO_119 *= DENDRO_97;
DENDRO_119 *= DENDRO_144;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_16;
DENDRO_120 *= DENDRO_54;
DENDRO_120 *= DENDRO_68;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_16;
DENDRO_68 *= DENDRO_54;
DENDRO_68 *= DENDRO_69;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_54;
DENDRO_69 *= DENDRO_70;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_54;
DENDRO_70 *= DENDRO_71;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_16;
DENDRO_71 *= DENDRO_54;
DENDRO_71 *= DENDRO_72;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_54;
DENDRO_72 *= DENDRO_73;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_16;
DENDRO_73 *= DENDRO_26;
DENDRO_73 *= DENDRO_74;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_16;
DENDRO_74 *= DENDRO_26;
DENDRO_74 *= DENDRO_75;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_16;
DENDRO_75 *= DENDRO_26;
DENDRO_75 *= DENDRO_76;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_16;
DENDRO_76 *= DENDRO_26;
DENDRO_76 *= DENDRO_78;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_16;
DENDRO_78 *= DENDRO_26;
DENDRO_78 *= DENDRO_79;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_26;
DENDRO_79 *= DENDRO_80;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_16;
DENDRO_80 *= DENDRO_7;
DENDRO_80 *= DENDRO_81;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_16;
DENDRO_81 *= DENDRO_7;
DENDRO_81 *= DENDRO_82;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_16;
DENDRO_82 *= DENDRO_7;
DENDRO_82 *= DENDRO_83;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_16;
DENDRO_83 *= DENDRO_7;
DENDRO_83 *= DENDRO_85;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_16;
DENDRO_85 *= DENDRO_7;
DENDRO_85 *= DENDRO_86;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_16;
DENDRO_86 *= DENDRO_7;
DENDRO_86 *= DENDRO_87;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_16;
DENDRO_87 *= DENDRO_3;
DENDRO_87 *= DENDRO_88;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_16;
DENDRO_88 *= DENDRO_3;
DENDRO_88 *= DENDRO_89;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_16;
DENDRO_89 *= DENDRO_27;
DENDRO_89 *= DENDRO_90;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_16;
DENDRO_90 *= DENDRO_27;
DENDRO_90 *= DENDRO_92;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_16;
DENDRO_92 *= DENDRO_17;
DENDRO_92 *= DENDRO_93;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_16;
DENDRO_93 *= DENDRO_17;
DENDRO_93 *= DENDRO_95;
DENDRO_95 = 1.0;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_95;
DENDRO_95 = grad_2_Gt2;
DENDRO_121 *= DENDRO_95;
DENDRO_121 *= DENDRO_12;
DENDRO_122 = 1.0;

// initialize reduction 
DENDRO_123 = 1;
DENDRO_123 *= DENDRO_122;
DENDRO_122 = grad_2_Gt1;
DENDRO_123 *= DENDRO_122;
DENDRO_123 *= DENDRO_4;
DENDRO_124 = 1.0;

// initialize reduction 
DENDRO_125 = 1;
DENDRO_125 *= DENDRO_124;
DENDRO_124 = grad_2_Gt0;
DENDRO_125 *= DENDRO_124;
DENDRO_125 *= DENDRO_2;
DENDRO_127 = 0.5;

// initialize reduction 
DENDRO_129 = 1;
DENDRO_129 *= DENDRO_127;
DENDRO_129 *= DENDRO_56;
DENDRO_129 *= DENDRO_65;
DENDRO_127 = 0.5;

// initialize reduction 
DENDRO_130 = 1;
DENDRO_130 *= DENDRO_127;
DENDRO_130 *= DENDRO_42;
DENDRO_130 *= DENDRO_66;
DENDRO_127 = 0.5;

// initialize reduction 
DENDRO_131 = 1;
DENDRO_131 *= DENDRO_127;
DENDRO_131 *= DENDRO_46;
DENDRO_131 *= DENDRO_67;
DENDRO_127 = -1.0/4.0;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_127;
DENDRO_133 *= DENDRO_105;
DENDRO_133 *= DENDRO_99;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_62;
DENDRO_99 *= DENDRO_97;
DENDRO_99 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_99;
DENDRO_6 += DENDRO_133;
DENDRO_6 += DENDRO_131;
DENDRO_6 += DENDRO_130;
DENDRO_6 += DENDRO_129;
DENDRO_6 += DENDRO_125;
DENDRO_6 += DENDRO_123;
DENDRO_6 += DENDRO_121;
DENDRO_6 += DENDRO_93;
DENDRO_6 += DENDRO_92;
DENDRO_6 += DENDRO_90;
DENDRO_6 += DENDRO_89;
DENDRO_6 += DENDRO_88;
DENDRO_6 += DENDRO_87;
DENDRO_6 += DENDRO_86;
DENDRO_6 += DENDRO_85;
DENDRO_6 += DENDRO_83;
DENDRO_6 += DENDRO_82;
DENDRO_6 += DENDRO_81;
DENDRO_6 += DENDRO_80;
DENDRO_6 += DENDRO_79;
DENDRO_6 += DENDRO_78;
DENDRO_6 += DENDRO_76;
DENDRO_6 += DENDRO_75;
DENDRO_6 += DENDRO_74;
DENDRO_6 += DENDRO_73;
DENDRO_6 += DENDRO_72;
DENDRO_6 += DENDRO_71;
DENDRO_6 += DENDRO_70;
DENDRO_6 += DENDRO_69;
DENDRO_6 += DENDRO_68;
DENDRO_6 += DENDRO_120;
DENDRO_6 += DENDRO_119;
DENDRO_6 += DENDRO_118;
DENDRO_6 += DENDRO_117;
DENDRO_6 += DENDRO_114;
DENDRO_6 += DENDRO_112;
DENDRO_6 += DENDRO_111;
DENDRO_6 += DENDRO_110;
DENDRO_6 += DENDRO_108;
DENDRO_6 += DENDRO_107;
DENDRO_6 += DENDRO_96;
double DENDRO_RIJ22 = DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_106;
DENDRO_6 *= DENDRO_38;
DENDRO_6 *= DENDRO_37;
DENDRO_68 = 0.5;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_68;
DENDRO_69 *= DENDRO_42;
DENDRO_69 *= DENDRO_39;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_109;
DENDRO_70 *= DENDRO_30;
DENDRO_70 *= DENDRO_53;
DENDRO_71 = 0.5;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_71;
DENDRO_72 *= DENDRO_46;
DENDRO_72 *= DENDRO_52;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_1;
DENDRO_73 *= DENDRO_50;
DENDRO_74 = 0.5;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_74;
DENDRO_75 *= DENDRO_56;
DENDRO_75 *= DENDRO_5;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_94;
DENDRO_76 *= DENDRO_46;
DENDRO_76 *= DENDRO_29;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_1;
DENDRO_78 *= DENDRO_35;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_74;
DENDRO_79 *= DENDRO_56;
DENDRO_79 *= DENDRO_5;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_94;
DENDRO_80 *= DENDRO_46;
DENDRO_80 *= DENDRO_29;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_24;
DENDRO_81 *= DENDRO_104;
DENDRO_82 = 0.5;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_82;
DENDRO_83 *= DENDRO_34;
DENDRO_83 *= DENDRO_49;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_13;
DENDRO_82 *= DENDRO_41;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_24;
DENDRO_85 *= DENDRO_39;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_25;
DENDRO_86 *= DENDRO_52;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_20;
DENDRO_87 *= DENDRO_53;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_15;
DENDRO_88 *= DENDRO_58;
DENDRO_89 = 0.5;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_89;
DENDRO_90 *= DENDRO_38;
DENDRO_90 *= DENDRO_104;
DENDRO_92 = 0.5;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_92;
DENDRO_93 *= DENDRO_42;
DENDRO_93 *= DENDRO_5;
DENDRO_94 = 0.5;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_94;
DENDRO_96 *= DENDRO_46;
DENDRO_96 *= DENDRO_33;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_92;
DENDRO_99 *= DENDRO_42;
DENDRO_99 *= DENDRO_5;

// initialize reduction 
DENDRO_107 = 1;
DENDRO_107 *= DENDRO_94;
DENDRO_107 *= DENDRO_46;
DENDRO_107 *= DENDRO_33;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_1;
DENDRO_94 *= DENDRO_29;
DENDRO_108 = 0.5;

// initialize reduction 
DENDRO_110 = 1;
DENDRO_110 *= DENDRO_108;
DENDRO_110 *= DENDRO_30;
DENDRO_110 *= DENDRO_57;
DENDRO_111 = 0.5;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_111;
DENDRO_112 *= DENDRO_34;
DENDRO_112 *= DENDRO_50;

// initialize reduction 
DENDRO_114 = 1;
DENDRO_114 *= DENDRO_13;
DENDRO_114 *= DENDRO_52;
DENDRO_117 = 0.5;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_117;
DENDRO_118 *= DENDRO_38;
DENDRO_118 *= DENDRO_39;
DENDRO_119 = 0.5;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_119;
DENDRO_120 *= DENDRO_34;
DENDRO_120 *= DENDRO_41;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_13;
DENDRO_121 *= DENDRO_47;
DENDRO_123 = 0.5;

// initialize reduction 
DENDRO_125 = 1;
DENDRO_125 *= DENDRO_123;
DENDRO_125 *= DENDRO_51;
DENDRO_125 *= DENDRO_37;

// initialize reduction 
DENDRO_129 = 1;
DENDRO_129 *= DENDRO_1;
DENDRO_129 *= DENDRO_47;

// initialize reduction 
DENDRO_130 = 1;
DENDRO_130 *= DENDRO_22;
DENDRO_130 *= DENDRO_39;

// initialize reduction 
DENDRO_131 = 1;
DENDRO_131 *= DENDRO_18;
DENDRO_131 *= DENDRO_52;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_20;
DENDRO_133 *= DENDRO_57;

// initialize reduction 
DENDRO_134 = 1;
DENDRO_134 *= DENDRO_15;
DENDRO_134 *= DENDRO_53;
DENDRO_135 = 0.5;

// initialize reduction 
DENDRO_137 = 1;
DENDRO_137 *= DENDRO_135;
DENDRO_137 *= DENDRO_56;
DENDRO_137 *= DENDRO_33;
DENDRO_138 = 0.5;

// initialize reduction 
DENDRO_139 = 1;
DENDRO_139 *= DENDRO_138;
DENDRO_139 *= DENDRO_38;
DENDRO_139 *= DENDRO_35;
DENDRO_138 = 0.5;

// initialize reduction 
DENDRO_140 = 1;
DENDRO_140 *= DENDRO_138;
DENDRO_140 *= DENDRO_42;
DENDRO_140 *= DENDRO_29;
DENDRO_141 = 0.5;

// initialize reduction 
DENDRO_142 = 1;
DENDRO_142 *= DENDRO_141;
DENDRO_142 *= DENDRO_46;
DENDRO_142 *= DENDRO_57;

// initialize reduction 
DENDRO_141 = 1;
DENDRO_141 *= DENDRO_1;
DENDRO_141 *= DENDRO_53;

// initialize reduction 
DENDRO_143 = 1;
DENDRO_143 *= DENDRO_135;
DENDRO_143 *= DENDRO_56;
DENDRO_143 *= DENDRO_33;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_138;
DENDRO_135 *= DENDRO_42;
DENDRO_135 *= DENDRO_29;

// initialize reduction 
DENDRO_145 = 1;
DENDRO_145 *= DENDRO_24;
DENDRO_145 *= DENDRO_29;
DENDRO_146 = 0.5;

// initialize reduction 
DENDRO_147 = 1;
DENDRO_147 *= DENDRO_146;
DENDRO_147 *= DENDRO_34;
DENDRO_147 *= DENDRO_58;

// initialize reduction 
DENDRO_146 = 1;
DENDRO_146 *= DENDRO_13;
DENDRO_146 *= DENDRO_53;

// initialize reduction 
DENDRO_148 = 1;
DENDRO_148 *= DENDRO_25;
DENDRO_148 *= DENDRO_57;
DENDRO_149 = 0.5;

// initialize reduction 
DENDRO_150 = 1;
DENDRO_150 *= DENDRO_149;
DENDRO_150 *= DENDRO_51;
DENDRO_150 *= DENDRO_49;

// initialize reduction 
DENDRO_149 = 1;
DENDRO_149 *= DENDRO_24;
DENDRO_149 *= DENDRO_47;

// initialize reduction 
DENDRO_151 = 1;
DENDRO_151 *= DENDRO_22;
DENDRO_151 *= DENDRO_41;
DENDRO_152 = 1.0;

// initialize reduction 
DENDRO_153 = 1;
DENDRO_153 *= DENDRO_152;
DENDRO_153 *= DENDRO_38;
DENDRO_153 *= DENDRO_41;
DENDRO_152 = 0.5;

// initialize reduction 
DENDRO_154 = 1;
DENDRO_154 *= DENDRO_152;
DENDRO_154 *= DENDRO_42;
DENDRO_154 *= DENDRO_47;
DENDRO_152 = 0.5;

// initialize reduction 
DENDRO_155 = 1;
DENDRO_155 *= DENDRO_152;
DENDRO_155 *= DENDRO_34;
DENDRO_155 *= DENDRO_37;

// initialize reduction 
DENDRO_156 = 1;
DENDRO_156 *= DENDRO_1;
DENDRO_156 *= DENDRO_39;

// initialize reduction 
DENDRO_157 = 1;
DENDRO_157 *= DENDRO_13;
DENDRO_157 *= DENDRO_39;
DENDRO_158 = 0.5;

// initialize reduction 
DENDRO_159 = 1;
DENDRO_159 *= DENDRO_158;
DENDRO_159 *= DENDRO_30;
DENDRO_159 *= DENDRO_52;

// initialize reduction 
DENDRO_160 = 1;
DENDRO_160 *= DENDRO_20;
DENDRO_160 *= DENDRO_52;

// initialize reduction 
DENDRO_161 = 1;
DENDRO_161 *= DENDRO_15;
DENDRO_161 *= DENDRO_50;
DENDRO_162 = 1.0;

// initialize reduction 
DENDRO_163 = 1;
DENDRO_163 *= DENDRO_162;
DENDRO_163 *= DENDRO_46;
DENDRO_163 *= DENDRO_5;

// initialize reduction 
DENDRO_162 = 1;
DENDRO_162 *= DENDRO_1;
DENDRO_162 *= DENDRO_104;
DENDRO_164 = 0.5;

// initialize reduction 
DENDRO_165 = 1;
DENDRO_165 *= DENDRO_164;
DENDRO_165 *= DENDRO_34;
DENDRO_165 *= DENDRO_53;

// initialize reduction 
DENDRO_164 = 1;
DENDRO_164 *= DENDRO_18;
DENDRO_164 *= DENDRO_57;

// initialize reduction 
DENDRO_166 = 1;
DENDRO_166 *= DENDRO_13;
DENDRO_166 *= DENDRO_57;
DENDRO_167 = 0.5;

// initialize reduction 
DENDRO_168 = 1;
DENDRO_168 *= DENDRO_167;
DENDRO_168 *= DENDRO_38;
DENDRO_168 *= DENDRO_47;
DENDRO_169 = 0.5;

// initialize reduction 
DENDRO_170 = 1;
DENDRO_170 *= DENDRO_169;
DENDRO_170 *= DENDRO_51;
DENDRO_170 *= DENDRO_41;

// initialize reduction 
DENDRO_171 = 1;
DENDRO_171 *= DENDRO_22;
DENDRO_171 *= DENDRO_47;
DENDRO_172 = 1.0;

// initialize reduction 
DENDRO_173 = 1;
DENDRO_173 *= DENDRO_172;
DENDRO_173 *= DENDRO_42;
DENDRO_173 *= DENDRO_33;
DENDRO_172 = 0.5;

// initialize reduction 
DENDRO_174 = 1;
DENDRO_174 *= DENDRO_172;
DENDRO_174 *= DENDRO_38;
DENDRO_174 *= DENDRO_29;

// initialize reduction 
DENDRO_172 = 1;
DENDRO_172 *= DENDRO_116;
DENDRO_172 *= DENDRO_38;
DENDRO_172 *= DENDRO_49;
DENDRO_116 = 0.5;

// initialize reduction 
DENDRO_175 = 1;
DENDRO_175 *= DENDRO_116;
DENDRO_175 *= DENDRO_42;
DENDRO_175 *= DENDRO_41;
DENDRO_116 = 0.5;

// initialize reduction 
DENDRO_176 = 1;
DENDRO_176 *= DENDRO_116;
DENDRO_176 *= DENDRO_46;
DENDRO_176 *= DENDRO_53;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_1;
DENDRO_116 *= DENDRO_58;
DENDRO_177 = 1.0;

// initialize reduction 
DENDRO_178 = 1;
DENDRO_178 *= DENDRO_177;
DENDRO_178 *= DENDRO_56;
DENDRO_178 *= DENDRO_29;

// initialize reduction 
DENDRO_177 = 1;
DENDRO_177 *= DENDRO_24;
DENDRO_177 *= DENDRO_35;

// initialize reduction 
DENDRO_179 = 1;
DENDRO_179 *= DENDRO_8;
DENDRO_179 *= DENDRO_98;
DENDRO_179 *= DENDRO_29;

// initialize reduction 
DENDRO_180 = 1;
DENDRO_180 *= DENDRO_8;
DENDRO_180 *= DENDRO_100;
DENDRO_180 *= DENDRO_41;

// initialize reduction 
DENDRO_181 = 1;
DENDRO_181 *= DENDRO_8;
DENDRO_181 *= DENDRO_102;
DENDRO_181 *= DENDRO_53;

// initialize reduction 
DENDRO_182 = 0;
DENDRO_182 += DENDRO_115;
DENDRO_182 += DENDRO_69;
DENDRO_182 += DENDRO_6;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_73;
DENDRO_6 += DENDRO_72;
DENDRO_6 += DENDRO_70;

// initialize reduction 
DENDRO_69 = 0;
DENDRO_69 += DENDRO_78;
DENDRO_69 += DENDRO_76;
DENDRO_69 += DENDRO_75;

// initialize reduction 
DENDRO_70 = 0;
DENDRO_70 += DENDRO_81;
DENDRO_70 += DENDRO_80;
DENDRO_70 += DENDRO_79;

// initialize reduction 
DENDRO_72 = 0;
DENDRO_72 += DENDRO_85;
DENDRO_72 += DENDRO_82;
DENDRO_72 += DENDRO_83;

// initialize reduction 
DENDRO_75 = 0;
DENDRO_75 += DENDRO_88;
DENDRO_75 += DENDRO_87;
DENDRO_75 += DENDRO_86;

// initialize reduction 
DENDRO_76 = 0;
DENDRO_76 += DENDRO_96;
DENDRO_76 += DENDRO_93;
DENDRO_76 += DENDRO_90;

// initialize reduction 
DENDRO_78 = 0;
DENDRO_78 += DENDRO_94;
DENDRO_78 += DENDRO_107;
DENDRO_78 += DENDRO_99;

// initialize reduction 
DENDRO_79 = 0;
DENDRO_79 += DENDRO_114;
DENDRO_79 += DENDRO_112;
DENDRO_79 += DENDRO_110;

// initialize reduction 
DENDRO_80 = 0;
DENDRO_80 += DENDRO_121;
DENDRO_80 += DENDRO_120;
DENDRO_80 += DENDRO_118;

// initialize reduction 
DENDRO_82 = 0;
DENDRO_82 += DENDRO_130;
DENDRO_82 += DENDRO_129;
DENDRO_82 += DENDRO_125;

// initialize reduction 
DENDRO_83 = 0;
DENDRO_83 += DENDRO_134;
DENDRO_83 += DENDRO_133;
DENDRO_83 += DENDRO_131;

// initialize reduction 
DENDRO_87 = 0;
DENDRO_87 += DENDRO_140;
DENDRO_87 += DENDRO_139;
DENDRO_87 += DENDRO_137;

// initialize reduction 
DENDRO_88 = 0;
DENDRO_88 += DENDRO_141;
DENDRO_88 += DENDRO_136;
DENDRO_88 += DENDRO_142;

// initialize reduction 
DENDRO_90 = 0;
DENDRO_90 += DENDRO_145;
DENDRO_90 += DENDRO_135;
DENDRO_90 += DENDRO_143;

// initialize reduction 
DENDRO_93 = 0;
DENDRO_93 += DENDRO_148;
DENDRO_93 += DENDRO_146;
DENDRO_93 += DENDRO_147;

// initialize reduction 
DENDRO_94 = 0;
DENDRO_94 += DENDRO_151;
DENDRO_94 += DENDRO_149;
DENDRO_94 += DENDRO_150;

// initialize reduction 
DENDRO_96 = 0;
DENDRO_96 += DENDRO_154;
DENDRO_96 += DENDRO_153;

// initialize reduction 
DENDRO_99 = 0;
DENDRO_99 += DENDRO_157;
DENDRO_99 += DENDRO_156;
DENDRO_99 += DENDRO_155;

// initialize reduction 
DENDRO_107 = 0;
DENDRO_107 += DENDRO_161;
DENDRO_107 += DENDRO_160;
DENDRO_107 += DENDRO_159;

// initialize reduction 
DENDRO_110 = 0;
DENDRO_110 += DENDRO_162;
DENDRO_110 += DENDRO_163;

// initialize reduction 
DENDRO_112 = 0;
DENDRO_112 += DENDRO_166;
DENDRO_112 += DENDRO_164;
DENDRO_112 += DENDRO_165;

// initialize reduction 
DENDRO_115 = 0;
DENDRO_115 += DENDRO_171;
DENDRO_115 += DENDRO_170;
DENDRO_115 += DENDRO_168;

// initialize reduction 
DENDRO_118 = 0;
DENDRO_118 += DENDRO_174;
DENDRO_118 += DENDRO_173;

// initialize reduction 
DENDRO_120 = 0;
DENDRO_120 += DENDRO_126;
DENDRO_120 += DENDRO_175;
DENDRO_120 += DENDRO_172;

// initialize reduction 
DENDRO_121 = 0;
DENDRO_121 += DENDRO_128;
DENDRO_121 += DENDRO_116;
DENDRO_121 += DENDRO_176;

// initialize reduction 
DENDRO_116 = 0;
DENDRO_116 += DENDRO_177;
DENDRO_116 += DENDRO_178;

// initialize reduction 
DENDRO_125 = 0;
DENDRO_125 += DENDRO_84;
DENDRO_125 += DENDRO_181;
DENDRO_125 += DENDRO_180;
DENDRO_125 += DENDRO_179;
DENDRO_84 = -1.0;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_84;
DENDRO_84 = grad2_1_2_gt4;
DENDRO_126 *= DENDRO_84;
DENDRO_126 *= DENDRO_16;
DENDRO_126 *= DENDRO_7;
DENDRO_84 = -1.0;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_84;
DENDRO_84 = grad2_0_2_gt4;
DENDRO_128 *= DENDRO_84;
DENDRO_128 *= DENDRO_16;
DENDRO_128 *= DENDRO_54;
DENDRO_84 = -1.0;

// initialize reduction 
DENDRO_129 = 1;
DENDRO_129 *= DENDRO_84;
DENDRO_84 = grad2_0_1_gt4;
DENDRO_129 *= DENDRO_84;
DENDRO_129 *= DENDRO_16;
DENDRO_129 *= DENDRO_26;
DENDRO_84 = -0.5;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_84;
DENDRO_84 = grad2_2_2_gt4;
DENDRO_133 *= DENDRO_84;
DENDRO_133 *= DENDRO_16;
DENDRO_133 *= DENDRO_17;
DENDRO_84 = -0.5;

// initialize reduction 
DENDRO_134 = 1;
DENDRO_134 *= DENDRO_84;
DENDRO_84 = grad2_1_1_gt4;
DENDRO_134 *= DENDRO_84;
DENDRO_134 *= DENDRO_16;
DENDRO_134 *= DENDRO_27;
DENDRO_84 = -0.5;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_84;
DENDRO_84 = grad2_0_0_gt4;
DENDRO_135 *= DENDRO_84;
DENDRO_135 *= DENDRO_16;
DENDRO_135 *= DENDRO_3;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_127;
DENDRO_84 *= DENDRO_100;
DENDRO_84 *= DENDRO_98;
DENDRO_84 *= DENDRO_105;

// initialize reduction 
DENDRO_137 = 1;
DENDRO_137 *= DENDRO_62;
DENDRO_137 *= DENDRO_4;
DENDRO_137 *= DENDRO_97;
DENDRO_137 *= DENDRO_144;

// initialize reduction 
DENDRO_139 = 1;
DENDRO_139 *= DENDRO_16;
DENDRO_139 *= DENDRO_54;
DENDRO_139 *= DENDRO_182;

// initialize reduction 
DENDRO_140 = 1;
DENDRO_140 *= DENDRO_16;
DENDRO_140 *= DENDRO_54;
DENDRO_140 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_16;
DENDRO_6 *= DENDRO_54;
DENDRO_6 *= DENDRO_69;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_54;
DENDRO_69 *= DENDRO_70;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_54;
DENDRO_70 *= DENDRO_72;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_54;
DENDRO_72 *= DENDRO_75;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_16;
DENDRO_75 *= DENDRO_26;
DENDRO_75 *= DENDRO_76;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_16;
DENDRO_76 *= DENDRO_26;
DENDRO_76 *= DENDRO_78;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_16;
DENDRO_78 *= DENDRO_26;
DENDRO_78 *= DENDRO_79;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_26;
DENDRO_79 *= DENDRO_80;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_16;
DENDRO_80 *= DENDRO_26;
DENDRO_80 *= DENDRO_82;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_16;
DENDRO_82 *= DENDRO_26;
DENDRO_82 *= DENDRO_83;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_16;
DENDRO_83 *= DENDRO_7;
DENDRO_83 *= DENDRO_87;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_16;
DENDRO_87 *= DENDRO_7;
DENDRO_87 *= DENDRO_88;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_16;
DENDRO_88 *= DENDRO_7;
DENDRO_88 *= DENDRO_90;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_16;
DENDRO_90 *= DENDRO_7;
DENDRO_90 *= DENDRO_93;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_16;
DENDRO_93 *= DENDRO_7;
DENDRO_93 *= DENDRO_94;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_16;
DENDRO_94 *= DENDRO_7;
DENDRO_94 *= DENDRO_96;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_16;
DENDRO_96 *= DENDRO_3;
DENDRO_96 *= DENDRO_99;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_16;
DENDRO_99 *= DENDRO_3;
DENDRO_99 *= DENDRO_107;

// initialize reduction 
DENDRO_107 = 1;
DENDRO_107 *= DENDRO_16;
DENDRO_107 *= DENDRO_3;
DENDRO_107 *= DENDRO_110;

// initialize reduction 
DENDRO_110 = 1;
DENDRO_110 *= DENDRO_16;
DENDRO_110 *= DENDRO_27;
DENDRO_110 *= DENDRO_112;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_16;
DENDRO_112 *= DENDRO_27;
DENDRO_112 *= DENDRO_115;

// initialize reduction 
DENDRO_115 = 1;
DENDRO_115 *= DENDRO_16;
DENDRO_115 *= DENDRO_27;
DENDRO_115 *= DENDRO_118;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_16;
DENDRO_118 *= DENDRO_17;
DENDRO_118 *= DENDRO_120;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_16;
DENDRO_120 *= DENDRO_17;
DENDRO_120 *= DENDRO_121;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_16;
DENDRO_121 *= DENDRO_17;
DENDRO_121 *= DENDRO_116;
DENDRO_116 = 0.5;

// initialize reduction 
DENDRO_141 = 1;
DENDRO_141 *= DENDRO_116;
DENDRO_141 *= DENDRO_44;
DENDRO_141 *= DENDRO_65;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_44;
DENDRO_116 *= DENDRO_95;
DENDRO_116 *= DENDRO_4;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_142 = 1;
DENDRO_142 *= DENDRO_44;
DENDRO_142 *= DENDRO_122;
DENDRO_142 *= DENDRO_11;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_143 = 1;
DENDRO_143 *= DENDRO_44;
DENDRO_143 *= DENDRO_124;
DENDRO_143 *= DENDRO_0;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_145 = 1;
DENDRO_145 *= DENDRO_44;
DENDRO_145 *= DENDRO_40;
DENDRO_145 *= DENDRO_66;
DENDRO_40 = 0.5;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_40;
DENDRO_40 = grad_1_Gt2;
DENDRO_44 *= DENDRO_40;
DENDRO_44 *= DENDRO_12;
DENDRO_146 = 0.5;

// initialize reduction 
DENDRO_147 = 1;
DENDRO_147 *= DENDRO_146;
DENDRO_146 = grad_1_Gt1;
DENDRO_147 *= DENDRO_146;
DENDRO_147 *= DENDRO_4;
DENDRO_148 = 0.5;

// initialize reduction 
DENDRO_149 = 1;
DENDRO_149 *= DENDRO_148;
DENDRO_148 = grad_1_Gt0;
DENDRO_149 *= DENDRO_148;
DENDRO_149 *= DENDRO_2;
DENDRO_150 = 0.5;

// initialize reduction 
DENDRO_151 = 1;
DENDRO_151 *= DENDRO_150;
DENDRO_151 *= DENDRO_19;
DENDRO_151 *= DENDRO_67;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_62;
DENDRO_19 *= DENDRO_97;
DENDRO_19 *= DENDRO_125;

// initialize reduction 
DENDRO_125 = 0;
DENDRO_125 += DENDRO_19;
DENDRO_125 += DENDRO_151;
DENDRO_125 += DENDRO_149;
DENDRO_125 += DENDRO_147;
DENDRO_125 += DENDRO_44;
DENDRO_125 += DENDRO_145;
DENDRO_125 += DENDRO_143;
DENDRO_125 += DENDRO_142;
DENDRO_125 += DENDRO_116;
DENDRO_125 += DENDRO_141;
DENDRO_125 += DENDRO_121;
DENDRO_125 += DENDRO_120;
DENDRO_125 += DENDRO_118;
DENDRO_125 += DENDRO_115;
DENDRO_125 += DENDRO_112;
DENDRO_125 += DENDRO_110;
DENDRO_125 += DENDRO_107;
DENDRO_125 += DENDRO_99;
DENDRO_125 += DENDRO_96;
DENDRO_125 += DENDRO_94;
DENDRO_125 += DENDRO_93;
DENDRO_125 += DENDRO_90;
DENDRO_125 += DENDRO_88;
DENDRO_125 += DENDRO_87;
DENDRO_125 += DENDRO_83;
DENDRO_125 += DENDRO_82;
DENDRO_125 += DENDRO_80;
DENDRO_125 += DENDRO_79;
DENDRO_125 += DENDRO_78;
DENDRO_125 += DENDRO_76;
DENDRO_125 += DENDRO_75;
DENDRO_125 += DENDRO_72;
DENDRO_125 += DENDRO_70;
DENDRO_125 += DENDRO_69;
DENDRO_125 += DENDRO_6;
DENDRO_125 += DENDRO_140;
DENDRO_125 += DENDRO_139;
DENDRO_125 += DENDRO_137;
DENDRO_125 += DENDRO_84;
DENDRO_125 += DENDRO_135;
DENDRO_125 += DENDRO_134;
DENDRO_125 += DENDRO_133;
DENDRO_125 += DENDRO_129;
DENDRO_125 += DENDRO_128;
DENDRO_125 += DENDRO_126;
double DENDRO_RIJ12 = DENDRO_125;
DENDRO_6 = 1.0;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_6;
DENDRO_19 *= DENDRO_34;
DENDRO_19 *= DENDRO_41;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_117;
DENDRO_6 *= DENDRO_38;
DENDRO_6 *= DENDRO_39;
DENDRO_44 = 1.0;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_44;
DENDRO_69 *= DENDRO_38;
DENDRO_69 *= DENDRO_39;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_119;
DENDRO_44 *= DENDRO_34;
DENDRO_44 *= DENDRO_41;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_92;
DENDRO_70 *= DENDRO_42;
DENDRO_70 *= DENDRO_5;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_59;
DENDRO_72 *= DENDRO_1;
DENDRO_72 *= DENDRO_29;
DENDRO_75 = 0.5;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_75;
DENDRO_76 *= DENDRO_23;
DENDRO_76 *= DENDRO_53;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_59;
DENDRO_78 *= DENDRO_1;
DENDRO_78 *= DENDRO_52;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_59;
DENDRO_79 *= DENDRO_15;
DENDRO_79 *= DENDRO_53;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_59;
DENDRO_80 *= DENDRO_24;
DENDRO_80 *= DENDRO_5;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_13;
DENDRO_82 *= DENDRO_29;
DENDRO_83 = 1.0;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_83;
DENDRO_84 *= DENDRO_34;
DENDRO_84 *= DENDRO_47;
DENDRO_83 = 0.5;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_83;
DENDRO_87 *= DENDRO_51;
DENDRO_87 *= DENDRO_39;
DENDRO_83 = 1.0;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_83;
DENDRO_88 *= DENDRO_34;
DENDRO_88 *= DENDRO_52;
DENDRO_83 = 0.5;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_83;
DENDRO_90 *= DENDRO_23;
DENDRO_90 *= DENDRO_57;
DENDRO_93 = 1.0;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_93;
DENDRO_94 *= DENDRO_51;
DENDRO_94 *= DENDRO_39;
DENDRO_93 = 0.5;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_93;
DENDRO_96 *= DENDRO_34;
DENDRO_96 *= DENDRO_47;
DENDRO_93 = 1.0;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_93;
DENDRO_99 *= DENDRO_38;
DENDRO_99 *= DENDRO_5;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_13;
DENDRO_93 *= DENDRO_33;

// initialize reduction 
DENDRO_107 = 1;
DENDRO_107 *= DENDRO_59;
DENDRO_107 *= DENDRO_1;
DENDRO_107 *= DENDRO_33;

// initialize reduction 
DENDRO_110 = 1;
DENDRO_110 *= DENDRO_22;
DENDRO_110 *= DENDRO_5;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_59;
DENDRO_112 *= DENDRO_15;
DENDRO_112 *= DENDRO_57;

// initialize reduction 
DENDRO_115 = 1;
DENDRO_115 *= DENDRO_21;
DENDRO_115 *= DENDRO_52;
DENDRO_116 = 1.0;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_116;
DENDRO_118 *= DENDRO_51;
DENDRO_118 *= DENDRO_41;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_167;
DENDRO_116 *= DENDRO_38;
DENDRO_116 *= DENDRO_47;
DENDRO_119 = 1.0;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_119;
DENDRO_120 *= DENDRO_38;
DENDRO_120 *= DENDRO_29;
DENDRO_119 = 0.5;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_119;
DENDRO_121 *= DENDRO_42;
DENDRO_121 *= DENDRO_33;
DENDRO_119 = 1.0;

// initialize reduction 
DENDRO_125 = 1;
DENDRO_125 *= DENDRO_119;
DENDRO_125 *= DENDRO_38;
DENDRO_125 *= DENDRO_47;

// initialize reduction 
DENDRO_119 = 1;
DENDRO_119 *= DENDRO_169;
DENDRO_119 *= DENDRO_51;
DENDRO_119 *= DENDRO_41;
DENDRO_126 = 1.0;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_126;
DENDRO_128 *= DENDRO_34;
DENDRO_128 *= DENDRO_53;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_59;
DENDRO_126 *= DENDRO_24;
DENDRO_126 *= DENDRO_33;

// initialize reduction 
DENDRO_129 = 1;
DENDRO_129 *= DENDRO_22;
DENDRO_129 *= DENDRO_29;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_59;
DENDRO_133 *= DENDRO_1;
DENDRO_133 *= DENDRO_57;

// initialize reduction 
DENDRO_134 = 1;
DENDRO_134 *= DENDRO_21;
DENDRO_134 *= DENDRO_53;
DENDRO_135 = 0.5;

// initialize reduction 
DENDRO_137 = 1;
DENDRO_137 *= DENDRO_135;
DENDRO_137 *= DENDRO_23;
DENDRO_137 *= DENDRO_52;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_59;
DENDRO_135 *= DENDRO_15;
DENDRO_135 *= DENDRO_52;

// initialize reduction 
DENDRO_139 = 1;
DENDRO_139 *= DENDRO_59;
DENDRO_139 *= DENDRO_1;
DENDRO_139 *= DENDRO_5;

// initialize reduction 
DENDRO_140 = 1;
DENDRO_140 *= DENDRO_13;
DENDRO_140 *= DENDRO_5;
DENDRO_141 = 1.0;

// initialize reduction 
DENDRO_142 = 1;
DENDRO_142 *= DENDRO_141;
DENDRO_142 *= DENDRO_38;
DENDRO_142 *= DENDRO_33;

// initialize reduction 
DENDRO_141 = 1;
DENDRO_141 *= DENDRO_22;
DENDRO_141 *= DENDRO_33;
DENDRO_143 = 1.0;

// initialize reduction 
DENDRO_145 = 1;
DENDRO_145 *= DENDRO_143;
DENDRO_145 *= DENDRO_34;
DENDRO_145 *= DENDRO_57;

// initialize reduction 
DENDRO_143 = 1;
DENDRO_143 *= DENDRO_21;
DENDRO_143 *= DENDRO_57;

// initialize reduction 
DENDRO_147 = 1;
DENDRO_147 *= DENDRO_138;
DENDRO_147 *= DENDRO_42;
DENDRO_147 *= DENDRO_29;

// initialize reduction 
DENDRO_138 = 1;
DENDRO_138 *= DENDRO_59;
DENDRO_138 *= DENDRO_24;
DENDRO_138 *= DENDRO_29;

// initialize reduction 
DENDRO_149 = 1;
DENDRO_149 *= DENDRO_59;
DENDRO_149 *= DENDRO_1;
DENDRO_149 *= DENDRO_53;

// initialize reduction 
DENDRO_150 = 1;
DENDRO_150 *= DENDRO_8;
DENDRO_150 *= DENDRO_98;
DENDRO_150 *= DENDRO_33;

// initialize reduction 
DENDRO_151 = 1;
DENDRO_151 *= DENDRO_8;
DENDRO_151 *= DENDRO_100;
DENDRO_151 *= DENDRO_47;

// initialize reduction 
DENDRO_153 = 1;
DENDRO_153 *= DENDRO_8;
DENDRO_153 *= DENDRO_102;
DENDRO_153 *= DENDRO_57;

// initialize reduction 
DENDRO_154 = 0;
DENDRO_154 += DENDRO_6;
DENDRO_154 += DENDRO_19;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_44;
DENDRO_6 += DENDRO_69;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_72;
DENDRO_19 += DENDRO_70;

// initialize reduction 
DENDRO_44 = 0;
DENDRO_44 += DENDRO_78;
DENDRO_44 += DENDRO_76;

// initialize reduction 
DENDRO_69 = 0;
DENDRO_69 += DENDRO_131;
DENDRO_69 += DENDRO_79;

// initialize reduction 
DENDRO_70 = 0;
DENDRO_70 += DENDRO_82;
DENDRO_70 += DENDRO_80;

// initialize reduction 
DENDRO_72 = 0;
DENDRO_72 += DENDRO_87;
DENDRO_72 += DENDRO_84;

// initialize reduction 
DENDRO_76 = 0;
DENDRO_76 += DENDRO_90;
DENDRO_76 += DENDRO_88;

// initialize reduction 
DENDRO_78 = 0;
DENDRO_78 += DENDRO_96;
DENDRO_78 += DENDRO_94;

// initialize reduction 
DENDRO_79 = 0;
DENDRO_79 += DENDRO_93;
DENDRO_79 += DENDRO_99;

// initialize reduction 
DENDRO_80 = 0;
DENDRO_80 += DENDRO_110;
DENDRO_80 += DENDRO_107;

// initialize reduction 
DENDRO_82 = 0;
DENDRO_82 += DENDRO_115;
DENDRO_82 += DENDRO_112;

// initialize reduction 
DENDRO_84 = 0;
DENDRO_84 += DENDRO_116;
DENDRO_84 += DENDRO_118;

// initialize reduction 
DENDRO_87 = 0;
DENDRO_87 += DENDRO_121;
DENDRO_87 += DENDRO_120;

// initialize reduction 
DENDRO_88 = 0;
DENDRO_88 += DENDRO_119;
DENDRO_88 += DENDRO_125;

// initialize reduction 
DENDRO_90 = 0;
DENDRO_90 += DENDRO_164;
DENDRO_90 += DENDRO_128;

// initialize reduction 
DENDRO_93 = 0;
DENDRO_93 += DENDRO_129;
DENDRO_93 += DENDRO_126;

// initialize reduction 
DENDRO_94 = 0;
DENDRO_94 += DENDRO_134;
DENDRO_94 += DENDRO_133;

// initialize reduction 
DENDRO_96 = 0;
DENDRO_96 += DENDRO_135;
DENDRO_96 += DENDRO_137;

// initialize reduction 
DENDRO_99 = 0;
DENDRO_99 += DENDRO_140;
DENDRO_99 += DENDRO_139;

// initialize reduction 
DENDRO_107 = 0;
DENDRO_107 += DENDRO_141;
DENDRO_107 += DENDRO_142;

// initialize reduction 
DENDRO_112 = 0;
DENDRO_112 += DENDRO_143;
DENDRO_112 += DENDRO_145;

// initialize reduction 
DENDRO_116 = 0;
DENDRO_116 += DENDRO_138;
DENDRO_116 += DENDRO_147;

// initialize reduction 
DENDRO_118 = 0;
DENDRO_118 += DENDRO_136;
DENDRO_118 += DENDRO_149;

// initialize reduction 
DENDRO_119 = 0;
DENDRO_119 += DENDRO_63;
DENDRO_119 += DENDRO_153;
DENDRO_119 += DENDRO_151;
DENDRO_119 += DENDRO_150;
DENDRO_63 = 1.5;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_63;
DENDRO_120 *= DENDRO_38;
DENDRO_120 *= DENDRO_16;
DENDRO_120 *= DENDRO_17;
DENDRO_120 *= DENDRO_41;
DENDRO_63 = 1.5;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_63;
DENDRO_121 *= DENDRO_51;
DENDRO_121 *= DENDRO_16;
DENDRO_121 *= DENDRO_27;
DENDRO_121 *= DENDRO_47;
DENDRO_63 = 1.5;

// initialize reduction 
DENDRO_125 = 1;
DENDRO_125 *= DENDRO_63;
DENDRO_125 *= DENDRO_34;
DENDRO_125 *= DENDRO_16;
DENDRO_125 *= DENDRO_3;
DENDRO_125 *= DENDRO_39;
DENDRO_63 = -1.0;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_63;
DENDRO_63 = grad2_1_2_gt3;
DENDRO_126 *= DENDRO_63;
DENDRO_126 *= DENDRO_16;
DENDRO_126 *= DENDRO_7;
DENDRO_63 = -1.0;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_63;
DENDRO_63 = grad2_0_2_gt3;
DENDRO_128 *= DENDRO_63;
DENDRO_128 *= DENDRO_16;
DENDRO_128 *= DENDRO_54;
DENDRO_63 = -1.0;

// initialize reduction 
DENDRO_129 = 1;
DENDRO_129 *= DENDRO_63;
DENDRO_63 = grad2_0_1_gt3;
DENDRO_129 *= DENDRO_63;
DENDRO_129 *= DENDRO_16;
DENDRO_129 *= DENDRO_26;
DENDRO_63 = -0.5;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_63;
DENDRO_63 = grad2_2_2_gt3;
DENDRO_133 *= DENDRO_63;
DENDRO_133 *= DENDRO_16;
DENDRO_133 *= DENDRO_17;
DENDRO_63 = -0.5;

// initialize reduction 
DENDRO_134 = 1;
DENDRO_134 *= DENDRO_63;
DENDRO_63 = grad2_1_1_gt3;
DENDRO_134 *= DENDRO_63;
DENDRO_134 *= DENDRO_16;
DENDRO_134 *= DENDRO_27;
DENDRO_63 = -0.5;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_63;
DENDRO_63 = grad2_0_0_gt3;
DENDRO_135 *= DENDRO_63;
DENDRO_135 *= DENDRO_16;
DENDRO_135 *= DENDRO_3;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_62;
DENDRO_63 *= DENDRO_11;
DENDRO_63 *= DENDRO_97;
DENDRO_63 *= DENDRO_144;

// initialize reduction 
DENDRO_136 = 1;
DENDRO_136 *= DENDRO_16;
DENDRO_136 *= DENDRO_54;
DENDRO_136 *= DENDRO_154;

// initialize reduction 
DENDRO_137 = 1;
DENDRO_137 *= DENDRO_16;
DENDRO_137 *= DENDRO_54;
DENDRO_137 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_16;
DENDRO_6 *= DENDRO_54;
DENDRO_6 *= DENDRO_19;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_16;
DENDRO_19 *= DENDRO_54;
DENDRO_19 *= DENDRO_44;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_16;
DENDRO_44 *= DENDRO_54;
DENDRO_44 *= DENDRO_69;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_54;
DENDRO_69 *= DENDRO_70;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_26;
DENDRO_70 *= DENDRO_72;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_26;
DENDRO_72 *= DENDRO_76;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_16;
DENDRO_76 *= DENDRO_26;
DENDRO_76 *= DENDRO_78;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_16;
DENDRO_78 *= DENDRO_26;
DENDRO_78 *= DENDRO_79;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_26;
DENDRO_79 *= DENDRO_80;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_16;
DENDRO_80 *= DENDRO_26;
DENDRO_80 *= DENDRO_82;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_16;
DENDRO_82 *= DENDRO_7;
DENDRO_82 *= DENDRO_84;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_16;
DENDRO_84 *= DENDRO_7;
DENDRO_84 *= DENDRO_87;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_16;
DENDRO_87 *= DENDRO_7;
DENDRO_87 *= DENDRO_88;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_16;
DENDRO_88 *= DENDRO_7;
DENDRO_88 *= DENDRO_90;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_16;
DENDRO_90 *= DENDRO_7;
DENDRO_90 *= DENDRO_93;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_16;
DENDRO_93 *= DENDRO_7;
DENDRO_93 *= DENDRO_94;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_16;
DENDRO_94 *= DENDRO_3;
DENDRO_94 *= DENDRO_96;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_16;
DENDRO_96 *= DENDRO_3;
DENDRO_96 *= DENDRO_99;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_16;
DENDRO_99 *= DENDRO_27;
DENDRO_99 *= DENDRO_107;

// initialize reduction 
DENDRO_107 = 1;
DENDRO_107 *= DENDRO_16;
DENDRO_107 *= DENDRO_27;
DENDRO_107 *= DENDRO_112;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_16;
DENDRO_112 *= DENDRO_17;
DENDRO_112 *= DENDRO_116;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_16;
DENDRO_116 *= DENDRO_17;
DENDRO_116 *= DENDRO_118;
DENDRO_118 = 1.0;

// initialize reduction 
DENDRO_138 = 1;
DENDRO_138 *= DENDRO_118;
DENDRO_138 *= DENDRO_40;
DENDRO_138 *= DENDRO_4;
DENDRO_118 = 1.0;

// initialize reduction 
DENDRO_139 = 1;
DENDRO_139 *= DENDRO_118;
DENDRO_139 *= DENDRO_146;
DENDRO_139 *= DENDRO_11;
DENDRO_118 = 1.0;

// initialize reduction 
DENDRO_141 = 1;
DENDRO_141 *= DENDRO_118;
DENDRO_141 *= DENDRO_148;
DENDRO_141 *= DENDRO_0;
DENDRO_118 = 0.5;

// initialize reduction 
DENDRO_142 = 1;
DENDRO_142 *= DENDRO_118;
DENDRO_142 *= DENDRO_38;
DENDRO_142 *= DENDRO_65;
DENDRO_118 = 0.5;

// initialize reduction 
DENDRO_143 = 1;
DENDRO_143 *= DENDRO_118;
DENDRO_143 *= DENDRO_51;
DENDRO_143 *= DENDRO_66;
DENDRO_118 = 0.5;

// initialize reduction 
DENDRO_145 = 1;
DENDRO_145 *= DENDRO_118;
DENDRO_145 *= DENDRO_34;
DENDRO_145 *= DENDRO_67;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_127;
DENDRO_118 *= DENDRO_105;
DENDRO_118 *= DENDRO_101;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_62;
DENDRO_101 *= DENDRO_97;
DENDRO_101 *= DENDRO_119;

// initialize reduction 
DENDRO_119 = 0;
DENDRO_119 += DENDRO_101;
DENDRO_119 += DENDRO_118;
DENDRO_119 += DENDRO_145;
DENDRO_119 += DENDRO_143;
DENDRO_119 += DENDRO_142;
DENDRO_119 += DENDRO_141;
DENDRO_119 += DENDRO_139;
DENDRO_119 += DENDRO_138;
DENDRO_119 += DENDRO_116;
DENDRO_119 += DENDRO_112;
DENDRO_119 += DENDRO_107;
DENDRO_119 += DENDRO_99;
DENDRO_119 += DENDRO_96;
DENDRO_119 += DENDRO_94;
DENDRO_119 += DENDRO_93;
DENDRO_119 += DENDRO_90;
DENDRO_119 += DENDRO_88;
DENDRO_119 += DENDRO_87;
DENDRO_119 += DENDRO_84;
DENDRO_119 += DENDRO_82;
DENDRO_119 += DENDRO_80;
DENDRO_119 += DENDRO_79;
DENDRO_119 += DENDRO_78;
DENDRO_119 += DENDRO_76;
DENDRO_119 += DENDRO_72;
DENDRO_119 += DENDRO_70;
DENDRO_119 += DENDRO_69;
DENDRO_119 += DENDRO_44;
DENDRO_119 += DENDRO_19;
DENDRO_119 += DENDRO_6;
DENDRO_119 += DENDRO_137;
DENDRO_119 += DENDRO_136;
DENDRO_119 += DENDRO_63;
DENDRO_119 += DENDRO_135;
DENDRO_119 += DENDRO_134;
DENDRO_119 += DENDRO_133;
DENDRO_119 += DENDRO_129;
DENDRO_119 += DENDRO_128;
DENDRO_119 += DENDRO_126;
DENDRO_119 += DENDRO_125;
DENDRO_119 += DENDRO_121;
DENDRO_119 += DENDRO_120;
double DENDRO_RIJ11 = DENDRO_119;
DENDRO_6 = 0.5;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_6;
DENDRO_19 *= DENDRO_56;
DENDRO_19 *= DENDRO_31;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_44;
DENDRO_63 *= DENDRO_30;
DENDRO_63 *= DENDRO_35;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_44;
DENDRO_69 *= DENDRO_46;
DENDRO_69 *= DENDRO_104;
DENDRO_70 = 0.5;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_70;
DENDRO_72 *= DENDRO_42;
DENDRO_72 *= DENDRO_43;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_18;
DENDRO_70 *= DENDRO_37;
DENDRO_76 = 0.5;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_76;
DENDRO_78 *= DENDRO_23;
DENDRO_78 *= DENDRO_49;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_13;
DENDRO_76 *= DENDRO_37;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_24;
DENDRO_79 *= DENDRO_43;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_6;
DENDRO_80 *= DENDRO_56;
DENDRO_80 *= DENDRO_31;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_44;
DENDRO_6 *= DENDRO_46;
DENDRO_6 *= DENDRO_104;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_25;
DENDRO_82 *= DENDRO_104;
DENDRO_84 = 0.5;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_84;
DENDRO_87 *= DENDRO_45;
DENDRO_87 *= DENDRO_58;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_25;
DENDRO_84 *= DENDRO_55;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_20;
DENDRO_88 *= DENDRO_50;
DENDRO_90 = 1.0;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_90;
DENDRO_93 *= DENDRO_30;
DENDRO_93 *= DENDRO_50;
DENDRO_90 = 0.5;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_90;
DENDRO_94 *= DENDRO_46;
DENDRO_94 *= DENDRO_55;
DENDRO_90 = 0.5;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_90;
DENDRO_96 *= DENDRO_30;
DENDRO_96 *= DENDRO_29;
DENDRO_99 = 0.5;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_99;
DENDRO_101 *= DENDRO_42;
DENDRO_101 *= DENDRO_31;
DENDRO_107 = 0.5;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_107;
DENDRO_112 *= DENDRO_46;
DENDRO_112 *= DENDRO_5;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_99;
DENDRO_116 *= DENDRO_42;
DENDRO_116 *= DENDRO_31;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_107;
DENDRO_118 *= DENDRO_46;
DENDRO_118 *= DENDRO_5;

// initialize reduction 
DENDRO_119 = 1;
DENDRO_119 *= DENDRO_18;
DENDRO_119 *= DENDRO_104;
DENDRO_120 = 0.5;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_120;
DENDRO_121 *= DENDRO_38;
DENDRO_121 *= DENDRO_43;
DENDRO_125 = 0.5;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_125;
DENDRO_126 *= DENDRO_23;
DENDRO_126 *= DENDRO_41;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_158;
DENDRO_128 *= DENDRO_30;
DENDRO_128 *= DENDRO_52;
DENDRO_129 = 0.5;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_129;
DENDRO_133 *= DENDRO_23;
DENDRO_133 *= DENDRO_50;

// initialize reduction 
DENDRO_134 = 1;
DENDRO_134 *= DENDRO_13;
DENDRO_134 *= DENDRO_55;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_22;
DENDRO_135 *= DENDRO_43;

// initialize reduction 
DENDRO_136 = 1;
DENDRO_136 *= DENDRO_21;
DENDRO_136 *= DENDRO_37;
DENDRO_137 = 0.5;

// initialize reduction 
DENDRO_138 = 1;
DENDRO_138 *= DENDRO_137;
DENDRO_138 *= DENDRO_45;
DENDRO_138 *= DENDRO_53;

// initialize reduction 
DENDRO_139 = 1;
DENDRO_139 *= DENDRO_18;
DENDRO_139 *= DENDRO_55;

// initialize reduction 
DENDRO_141 = 1;
DENDRO_141 *= DENDRO_106;
DENDRO_141 *= DENDRO_38;
DENDRO_141 *= DENDRO_37;

// initialize reduction 
DENDRO_106 = 1;
DENDRO_106 *= DENDRO_68;
DENDRO_106 *= DENDRO_42;
DENDRO_106 *= DENDRO_39;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_18;
DENDRO_68 *= DENDRO_41;

// initialize reduction 
DENDRO_142 = 1;
DENDRO_142 *= DENDRO_109;
DENDRO_142 *= DENDRO_30;
DENDRO_142 *= DENDRO_53;

// initialize reduction 
DENDRO_143 = 1;
DENDRO_143 *= DENDRO_71;
DENDRO_143 *= DENDRO_46;
DENDRO_143 *= DENDRO_52;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_74;
DENDRO_71 *= DENDRO_56;
DENDRO_71 *= DENDRO_5;

// initialize reduction 
DENDRO_145 = 1;
DENDRO_145 *= DENDRO_91;
DENDRO_145 *= DENDRO_42;
DENDRO_145 *= DENDRO_104;

// initialize reduction 
DENDRO_147 = 1;
DENDRO_147 *= DENDRO_18;
DENDRO_147 *= DENDRO_35;

// initialize reduction 
DENDRO_149 = 1;
DENDRO_149 *= DENDRO_74;
DENDRO_149 *= DENDRO_56;
DENDRO_149 *= DENDRO_5;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_91;
DENDRO_74 *= DENDRO_42;
DENDRO_74 *= DENDRO_104;

// initialize reduction 
DENDRO_150 = 1;
DENDRO_150 *= DENDRO_25;
DENDRO_150 *= DENDRO_29;
DENDRO_151 = 0.5;

// initialize reduction 
DENDRO_153 = 1;
DENDRO_153 *= DENDRO_151;
DENDRO_153 *= DENDRO_23;
DENDRO_153 *= DENDRO_58;

// initialize reduction 
DENDRO_151 = 1;
DENDRO_151 *= DENDRO_13;
DENDRO_151 *= DENDRO_50;

// initialize reduction 
DENDRO_154 = 1;
DENDRO_154 *= DENDRO_22;
DENDRO_154 *= DENDRO_37;

// initialize reduction 
DENDRO_155 = 1;
DENDRO_155 *= DENDRO_21;
DENDRO_155 *= DENDRO_49;
DENDRO_159 = 0.5;

// initialize reduction 
DENDRO_163 = 1;
DENDRO_163 *= DENDRO_159;
DENDRO_163 *= DENDRO_23;
DENDRO_163 *= DENDRO_37;

// initialize reduction 
DENDRO_159 = 1;
DENDRO_159 *= DENDRO_1;
DENDRO_159 *= DENDRO_43;

// initialize reduction 
DENDRO_164 = 1;
DENDRO_164 *= DENDRO_13;
DENDRO_164 *= DENDRO_43;
DENDRO_165 = 0.5;

// initialize reduction 
DENDRO_166 = 1;
DENDRO_166 *= DENDRO_165;
DENDRO_166 *= DENDRO_30;
DENDRO_166 *= DENDRO_55;
DENDRO_167 = 0.5;

// initialize reduction 
DENDRO_168 = 1;
DENDRO_168 *= DENDRO_167;
DENDRO_168 *= DENDRO_45;
DENDRO_168 *= DENDRO_50;

// initialize reduction 
DENDRO_169 = 1;
DENDRO_169 *= DENDRO_20;
DENDRO_169 *= DENDRO_55;
DENDRO_170 = 1.0;

// initialize reduction 
DENDRO_171 = 1;
DENDRO_171 *= DENDRO_170;
DENDRO_171 *= DENDRO_46;
DENDRO_171 *= DENDRO_31;
DENDRO_170 = 0.5;

// initialize reduction 
DENDRO_172 = 1;
DENDRO_172 *= DENDRO_170;
DENDRO_172 *= DENDRO_30;
DENDRO_172 *= DENDRO_104;

// initialize reduction 
DENDRO_170 = 1;
DENDRO_170 *= DENDRO_75;
DENDRO_170 *= DENDRO_23;
DENDRO_170 *= DENDRO_53;

// initialize reduction 
DENDRO_173 = 1;
DENDRO_173 *= DENDRO_117;
DENDRO_173 *= DENDRO_38;
DENDRO_173 *= DENDRO_39;

// initialize reduction 
DENDRO_174 = 1;
DENDRO_174 *= DENDRO_21;
DENDRO_174 *= DENDRO_41;
DENDRO_175 = 1.0;

// initialize reduction 
DENDRO_176 = 1;
DENDRO_176 *= DENDRO_175;
DENDRO_176 *= DENDRO_42;
DENDRO_176 *= DENDRO_5;

// initialize reduction 
DENDRO_175 = 1;
DENDRO_175 *= DENDRO_18;
DENDRO_175 *= DENDRO_29;
DENDRO_177 = 0.5;

// initialize reduction 
DENDRO_178 = 1;
DENDRO_178 *= DENDRO_177;
DENDRO_178 *= DENDRO_42;
DENDRO_178 *= DENDRO_37;

// initialize reduction 
DENDRO_177 = 1;
DENDRO_177 *= DENDRO_18;
DENDRO_177 *= DENDRO_49;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_24;
DENDRO_49 *= DENDRO_37;

// initialize reduction 
DENDRO_179 = 1;
DENDRO_179 *= DENDRO_77;
DENDRO_179 *= DENDRO_30;
DENDRO_179 *= DENDRO_58;
DENDRO_58 = 0.5;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_58;
DENDRO_77 *= DENDRO_46;
DENDRO_77 *= DENDRO_50;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_25;
DENDRO_58 *= DENDRO_50;
DENDRO_180 = 1.0;

// initialize reduction 
DENDRO_181 = 1;
DENDRO_181 *= DENDRO_180;
DENDRO_181 *= DENDRO_56;
DENDRO_181 *= DENDRO_104;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_25;
DENDRO_56 *= DENDRO_35;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_8;
DENDRO_35 *= DENDRO_98;
DENDRO_35 *= DENDRO_104;

// initialize reduction 
DENDRO_180 = 1;
DENDRO_180 *= DENDRO_8;
DENDRO_180 *= DENDRO_100;
DENDRO_180 *= DENDRO_37;

// initialize reduction 
DENDRO_182 = 1;
DENDRO_182 *= DENDRO_8;
DENDRO_182 *= DENDRO_102;
DENDRO_182 *= DENDRO_50;

// initialize reduction 
DENDRO_183 = 0;
DENDRO_183 += DENDRO_69;
DENDRO_183 += DENDRO_63;
DENDRO_183 += DENDRO_19;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_132;
DENDRO_19 += DENDRO_70;
DENDRO_19 += DENDRO_72;

// initialize reduction 
DENDRO_63 = 0;
DENDRO_63 += DENDRO_79;
DENDRO_63 += DENDRO_76;
DENDRO_63 += DENDRO_78;

// initialize reduction 
DENDRO_69 = 0;
DENDRO_69 += DENDRO_82;
DENDRO_69 += DENDRO_6;
DENDRO_69 += DENDRO_80;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_88;
DENDRO_6 += DENDRO_84;
DENDRO_6 += DENDRO_87;

// initialize reduction 
DENDRO_70 = 0;
DENDRO_70 += DENDRO_94;
DENDRO_70 += DENDRO_93;

// initialize reduction 
DENDRO_72 = 0;
DENDRO_72 += DENDRO_112;
DENDRO_72 += DENDRO_101;
DENDRO_72 += DENDRO_96;

// initialize reduction 
DENDRO_76 = 0;
DENDRO_76 += DENDRO_119;
DENDRO_76 += DENDRO_118;
DENDRO_76 += DENDRO_116;

// initialize reduction 
DENDRO_78 = 0;
DENDRO_78 += DENDRO_157;
DENDRO_78 += DENDRO_126;
DENDRO_78 += DENDRO_121;

// initialize reduction 
DENDRO_79 = 0;
DENDRO_79 += DENDRO_134;
DENDRO_79 += DENDRO_133;
DENDRO_79 += DENDRO_128;

// initialize reduction 
DENDRO_80 = 0;
DENDRO_80 += DENDRO_136;
DENDRO_80 += DENDRO_135;
DENDRO_80 += DENDRO_156;

// initialize reduction 
DENDRO_82 = 0;
DENDRO_82 += DENDRO_160;
DENDRO_82 += DENDRO_139;
DENDRO_82 += DENDRO_138;

// initialize reduction 
DENDRO_84 = 0;
DENDRO_84 += DENDRO_68;
DENDRO_84 += DENDRO_106;
DENDRO_84 += DENDRO_141;

// initialize reduction 
DENDRO_87 = 0;
DENDRO_87 += DENDRO_113;
DENDRO_87 += DENDRO_143;
DENDRO_87 += DENDRO_142;

// initialize reduction 
DENDRO_88 = 0;
DENDRO_88 += DENDRO_147;
DENDRO_88 += DENDRO_145;
DENDRO_88 += DENDRO_71;

// initialize reduction 
DENDRO_71 = 0;
DENDRO_71 += DENDRO_150;
DENDRO_71 += DENDRO_74;
DENDRO_71 += DENDRO_149;

// initialize reduction 
DENDRO_74 = 0;
DENDRO_74 += DENDRO_86;
DENDRO_74 += DENDRO_151;
DENDRO_74 += DENDRO_153;

// initialize reduction 
DENDRO_86 = 0;
DENDRO_86 += DENDRO_155;
DENDRO_86 += DENDRO_154;
DENDRO_86 += DENDRO_85;

// initialize reduction 
DENDRO_85 = 0;
DENDRO_85 += DENDRO_164;
DENDRO_85 += DENDRO_159;
DENDRO_85 += DENDRO_163;

// initialize reduction 
DENDRO_93 = 0;
DENDRO_93 += DENDRO_169;
DENDRO_93 += DENDRO_168;
DENDRO_93 += DENDRO_166;

// initialize reduction 
DENDRO_94 = 0;
DENDRO_94 += DENDRO_172;
DENDRO_94 += DENDRO_171;

// initialize reduction 
DENDRO_96 = 0;
DENDRO_96 += DENDRO_114;
DENDRO_96 += DENDRO_131;
DENDRO_96 += DENDRO_170;

// initialize reduction 
DENDRO_101 = 0;
DENDRO_101 += DENDRO_174;
DENDRO_101 += DENDRO_130;
DENDRO_101 += DENDRO_173;

// initialize reduction 
DENDRO_106 = 0;
DENDRO_106 += DENDRO_175;
DENDRO_106 += DENDRO_176;

// initialize reduction 
DENDRO_112 = 0;
DENDRO_112 += DENDRO_49;
DENDRO_112 += DENDRO_177;
DENDRO_112 += DENDRO_178;

// initialize reduction 
DENDRO_49 = 0;
DENDRO_49 += DENDRO_58;
DENDRO_49 += DENDRO_77;
DENDRO_49 += DENDRO_179;

// initialize reduction 
DENDRO_58 = 0;
DENDRO_58 += DENDRO_56;
DENDRO_58 += DENDRO_181;

// initialize reduction 
DENDRO_56 = 0;
DENDRO_56 += DENDRO_60;
DENDRO_56 += DENDRO_182;
DENDRO_56 += DENDRO_180;
DENDRO_56 += DENDRO_35;
DENDRO_35 = -1.0;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_35;
DENDRO_35 = grad2_1_2_gt2;
DENDRO_60 *= DENDRO_35;
DENDRO_60 *= DENDRO_16;
DENDRO_60 *= DENDRO_7;
DENDRO_35 = -1.0;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_35;
DENDRO_35 = grad2_0_2_gt2;
DENDRO_77 *= DENDRO_35;
DENDRO_77 *= DENDRO_16;
DENDRO_77 *= DENDRO_54;
DENDRO_35 = -1.0;

// initialize reduction 
DENDRO_114 = 1;
DENDRO_114 *= DENDRO_35;
DENDRO_35 = grad2_0_1_gt2;
DENDRO_114 *= DENDRO_35;
DENDRO_114 *= DENDRO_16;
DENDRO_114 *= DENDRO_26;
DENDRO_35 = -0.5;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_35;
DENDRO_35 = grad2_2_2_gt2;
DENDRO_116 *= DENDRO_35;
DENDRO_116 *= DENDRO_16;
DENDRO_116 *= DENDRO_17;
DENDRO_35 = -0.5;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_35;
DENDRO_35 = grad2_1_1_gt2;
DENDRO_118 *= DENDRO_35;
DENDRO_118 *= DENDRO_16;
DENDRO_118 *= DENDRO_27;
DENDRO_35 = -0.5;

// initialize reduction 
DENDRO_119 = 1;
DENDRO_119 *= DENDRO_35;
DENDRO_35 = grad2_0_0_gt2;
DENDRO_119 *= DENDRO_35;
DENDRO_119 *= DENDRO_16;
DENDRO_119 *= DENDRO_3;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_127;
DENDRO_35 *= DENDRO_102;
DENDRO_35 *= DENDRO_98;
DENDRO_35 *= DENDRO_105;

// initialize reduction 
DENDRO_121 = 1;
DENDRO_121 *= DENDRO_62;
DENDRO_121 *= DENDRO_2;
DENDRO_121 *= DENDRO_97;
DENDRO_121 *= DENDRO_144;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_16;
DENDRO_126 *= DENDRO_54;
DENDRO_126 *= DENDRO_183;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_16;
DENDRO_128 *= DENDRO_54;
DENDRO_128 *= DENDRO_19;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_16;
DENDRO_19 *= DENDRO_54;
DENDRO_19 *= DENDRO_63;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_16;
DENDRO_63 *= DENDRO_54;
DENDRO_63 *= DENDRO_69;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_54;
DENDRO_69 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_16;
DENDRO_6 *= DENDRO_54;
DENDRO_6 *= DENDRO_70;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_26;
DENDRO_70 *= DENDRO_72;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_26;
DENDRO_72 *= DENDRO_76;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_16;
DENDRO_76 *= DENDRO_26;
DENDRO_76 *= DENDRO_78;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_16;
DENDRO_78 *= DENDRO_26;
DENDRO_78 *= DENDRO_79;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_26;
DENDRO_79 *= DENDRO_80;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_16;
DENDRO_80 *= DENDRO_26;
DENDRO_80 *= DENDRO_82;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_16;
DENDRO_82 *= DENDRO_7;
DENDRO_82 *= DENDRO_84;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_16;
DENDRO_84 *= DENDRO_7;
DENDRO_84 *= DENDRO_87;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_16;
DENDRO_87 *= DENDRO_7;
DENDRO_87 *= DENDRO_88;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_16;
DENDRO_88 *= DENDRO_7;
DENDRO_88 *= DENDRO_71;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_16;
DENDRO_71 *= DENDRO_7;
DENDRO_71 *= DENDRO_74;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_16;
DENDRO_74 *= DENDRO_7;
DENDRO_74 *= DENDRO_86;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_16;
DENDRO_86 *= DENDRO_3;
DENDRO_86 *= DENDRO_85;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_16;
DENDRO_85 *= DENDRO_3;
DENDRO_85 *= DENDRO_93;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_16;
DENDRO_93 *= DENDRO_3;
DENDRO_93 *= DENDRO_94;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_16;
DENDRO_94 *= DENDRO_27;
DENDRO_94 *= DENDRO_96;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_16;
DENDRO_96 *= DENDRO_27;
DENDRO_96 *= DENDRO_101;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_16;
DENDRO_101 *= DENDRO_27;
DENDRO_101 *= DENDRO_106;

// initialize reduction 
DENDRO_106 = 1;
DENDRO_106 *= DENDRO_16;
DENDRO_106 *= DENDRO_17;
DENDRO_106 *= DENDRO_112;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_16;
DENDRO_112 *= DENDRO_17;
DENDRO_112 *= DENDRO_49;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_16;
DENDRO_49 *= DENDRO_17;
DENDRO_49 *= DENDRO_58;
DENDRO_58 = 0.5;

// initialize reduction 
DENDRO_130 = 1;
DENDRO_130 *= DENDRO_58;
DENDRO_130 *= DENDRO_48;
DENDRO_130 *= DENDRO_65;
DENDRO_48 = 0.5;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_48;
DENDRO_58 *= DENDRO_95;
DENDRO_58 *= DENDRO_2;
DENDRO_48 = 0.5;

// initialize reduction 
DENDRO_95 = 1;
DENDRO_95 *= DENDRO_48;
DENDRO_95 *= DENDRO_122;
DENDRO_95 *= DENDRO_0;
DENDRO_48 = 0.5;

// initialize reduction 
DENDRO_122 = 1;
DENDRO_122 *= DENDRO_48;
DENDRO_122 *= DENDRO_124;
DENDRO_122 *= DENDRO_10;
DENDRO_48 = 0.5;

// initialize reduction 
DENDRO_124 = 1;
DENDRO_124 *= DENDRO_48;
DENDRO_124 *= DENDRO_9;
DENDRO_124 *= DENDRO_66;
DENDRO_9 = 0.5;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_9;
DENDRO_48 *= DENDRO_32;
DENDRO_48 *= DENDRO_67;
DENDRO_9 = 0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_9;
DENDRO_9 = grad_0_Gt2;
DENDRO_32 *= DENDRO_9;
DENDRO_32 *= DENDRO_12;
DENDRO_12 = 0.5;

// initialize reduction 
DENDRO_133 = 1;
DENDRO_133 *= DENDRO_12;
DENDRO_12 = grad_0_Gt1;
DENDRO_133 *= DENDRO_12;
DENDRO_133 *= DENDRO_4;
DENDRO_134 = 0.5;

// initialize reduction 
DENDRO_135 = 1;
DENDRO_135 *= DENDRO_134;
DENDRO_134 = grad_0_Gt0;
DENDRO_135 *= DENDRO_134;
DENDRO_135 *= DENDRO_2;

// initialize reduction 
DENDRO_136 = 1;
DENDRO_136 *= DENDRO_62;
DENDRO_136 *= DENDRO_97;
DENDRO_136 *= DENDRO_56;

// initialize reduction 
DENDRO_56 = 0;
DENDRO_56 += DENDRO_136;
DENDRO_56 += DENDRO_135;
DENDRO_56 += DENDRO_133;
DENDRO_56 += DENDRO_32;
DENDRO_56 += DENDRO_48;
DENDRO_56 += DENDRO_124;
DENDRO_56 += DENDRO_122;
DENDRO_56 += DENDRO_95;
DENDRO_56 += DENDRO_58;
DENDRO_56 += DENDRO_130;
DENDRO_56 += DENDRO_49;
DENDRO_56 += DENDRO_112;
DENDRO_56 += DENDRO_106;
DENDRO_56 += DENDRO_101;
DENDRO_56 += DENDRO_96;
DENDRO_56 += DENDRO_94;
DENDRO_56 += DENDRO_93;
DENDRO_56 += DENDRO_85;
DENDRO_56 += DENDRO_86;
DENDRO_56 += DENDRO_74;
DENDRO_56 += DENDRO_71;
DENDRO_56 += DENDRO_88;
DENDRO_56 += DENDRO_87;
DENDRO_56 += DENDRO_84;
DENDRO_56 += DENDRO_82;
DENDRO_56 += DENDRO_80;
DENDRO_56 += DENDRO_79;
DENDRO_56 += DENDRO_78;
DENDRO_56 += DENDRO_76;
DENDRO_56 += DENDRO_72;
DENDRO_56 += DENDRO_70;
DENDRO_56 += DENDRO_6;
DENDRO_56 += DENDRO_69;
DENDRO_56 += DENDRO_63;
DENDRO_56 += DENDRO_19;
DENDRO_56 += DENDRO_128;
DENDRO_56 += DENDRO_126;
DENDRO_56 += DENDRO_121;
DENDRO_56 += DENDRO_35;
DENDRO_56 += DENDRO_119;
DENDRO_56 += DENDRO_118;
DENDRO_56 += DENDRO_116;
DENDRO_56 += DENDRO_114;
DENDRO_56 += DENDRO_77;
DENDRO_56 += DENDRO_60;
double DENDRO_RIJ02 = DENDRO_56;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_120;
DENDRO_6 *= DENDRO_38;
DENDRO_6 *= DENDRO_43;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_125;
DENDRO_19 *= DENDRO_23;
DENDRO_19 *= DENDRO_41;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_152;
DENDRO_32 *= DENDRO_34;
DENDRO_32 *= DENDRO_37;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_120;
DENDRO_35 *= DENDRO_38;
DENDRO_35 *= DENDRO_43;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_152;
DENDRO_41 *= DENDRO_34;
DENDRO_41 *= DENDRO_37;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_18;
DENDRO_48 *= DENDRO_39;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_90;
DENDRO_49 *= DENDRO_30;
DENDRO_49 *= DENDRO_29;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_99;
DENDRO_29 *= DENDRO_42;
DENDRO_29 *= DENDRO_31;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_158;
DENDRO_56 *= DENDRO_30;
DENDRO_56 *= DENDRO_52;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_129;
DENDRO_58 *= DENDRO_23;
DENDRO_58 *= DENDRO_50;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_1;
DENDRO_60 *= DENDRO_55;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_13;
DENDRO_63 *= DENDRO_104;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_24;
DENDRO_69 *= DENDRO_31;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_25;
DENDRO_70 *= DENDRO_5;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_137;
DENDRO_71 *= DENDRO_45;
DENDRO_71 *= DENDRO_53;
DENDRO_72 = 0.5;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_72;
DENDRO_74 *= DENDRO_51;
DENDRO_74 *= DENDRO_43;
DENDRO_76 = 0.5;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_76;
DENDRO_77 *= DENDRO_23;
DENDRO_77 *= DENDRO_47;
DENDRO_76 = 0.5;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_76;
DENDRO_78 *= DENDRO_34;
DENDRO_78 *= DENDRO_39;
DENDRO_79 = 0.5;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_79;
DENDRO_80 *= DENDRO_38;
DENDRO_80 *= DENDRO_31;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_18;
DENDRO_79 *= DENDRO_5;
DENDRO_82 = 0.5;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_82;
DENDRO_84 *= DENDRO_30;
DENDRO_84 *= DENDRO_33;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_1;
DENDRO_82 *= DENDRO_5;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_22;
DENDRO_85 *= DENDRO_31;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_72;
DENDRO_86 *= DENDRO_51;
DENDRO_86 *= DENDRO_43;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_76;
DENDRO_72 *= DENDRO_34;
DENDRO_72 *= DENDRO_39;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_21;
DENDRO_87 *= DENDRO_39;
DENDRO_88 = 0.5;

// initialize reduction 
DENDRO_90 = 1;
DENDRO_90 *= DENDRO_88;
DENDRO_90 *= DENDRO_45;
DENDRO_90 *= DENDRO_57;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_21;
DENDRO_88 *= DENDRO_55;

// initialize reduction 
DENDRO_93 = 1;
DENDRO_93 *= DENDRO_15;
DENDRO_93 *= DENDRO_52;
DENDRO_94 = 1.0;

// initialize reduction 
DENDRO_95 = 1;
DENDRO_95 *= DENDRO_94;
DENDRO_95 *= DENDRO_23;
DENDRO_95 *= DENDRO_52;
DENDRO_94 = 0.5;

// initialize reduction 
DENDRO_96 = 1;
DENDRO_96 *= DENDRO_94;
DENDRO_96 *= DENDRO_34;
DENDRO_96 *= DENDRO_55;

// initialize reduction 
DENDRO_94 = 1;
DENDRO_94 *= DENDRO_89;
DENDRO_94 *= DENDRO_38;
DENDRO_94 *= DENDRO_104;

// initialize reduction 
DENDRO_89 = 1;
DENDRO_89 *= DENDRO_92;
DENDRO_89 *= DENDRO_42;
DENDRO_89 *= DENDRO_5;

// initialize reduction 
DENDRO_92 = 1;
DENDRO_92 *= DENDRO_75;
DENDRO_92 *= DENDRO_23;
DENDRO_92 *= DENDRO_53;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_111;
DENDRO_75 *= DENDRO_34;
DENDRO_75 *= DENDRO_50;

// initialize reduction 
DENDRO_99 = 1;
DENDRO_99 *= DENDRO_117;
DENDRO_99 *= DENDRO_38;
DENDRO_99 *= DENDRO_39;

// initialize reduction 
DENDRO_101 = 1;
DENDRO_101 *= DENDRO_123;
DENDRO_101 *= DENDRO_51;
DENDRO_101 *= DENDRO_37;

// initialize reduction 
DENDRO_106 = 1;
DENDRO_106 *= DENDRO_18;
DENDRO_106 *= DENDRO_47;

// initialize reduction 
DENDRO_111 = 1;
DENDRO_111 *= DENDRO_24;
DENDRO_111 *= DENDRO_5;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_25;
DENDRO_24 *= DENDRO_33;

// initialize reduction 
DENDRO_112 = 1;
DENDRO_112 *= DENDRO_22;
DENDRO_112 *= DENDRO_104;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_117;
DENDRO_22 *= DENDRO_38;
DENDRO_22 *= DENDRO_39;

// initialize reduction 
DENDRO_114 = 1;
DENDRO_114 *= DENDRO_123;
DENDRO_114 *= DENDRO_51;
DENDRO_114 *= DENDRO_37;

// initialize reduction 
DENDRO_116 = 1;
DENDRO_116 *= DENDRO_108;
DENDRO_116 *= DENDRO_30;
DENDRO_116 *= DENDRO_57;

// initialize reduction 
DENDRO_108 = 1;
DENDRO_108 *= DENDRO_1;
DENDRO_108 *= DENDRO_52;

// initialize reduction 
DENDRO_117 = 1;
DENDRO_117 *= DENDRO_21;
DENDRO_117 *= DENDRO_50;
DENDRO_118 = 0.5;

// initialize reduction 
DENDRO_119 = 1;
DENDRO_119 *= DENDRO_118;
DENDRO_119 *= DENDRO_30;
DENDRO_119 *= DENDRO_5;

// initialize reduction 
DENDRO_118 = 1;
DENDRO_118 *= DENDRO_1;
DENDRO_118 *= DENDRO_31;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_13;
DENDRO_1 *= DENDRO_31;
DENDRO_13 = 0.5;

// initialize reduction 
DENDRO_120 = 1;
DENDRO_120 *= DENDRO_13;
DENDRO_120 *= DENDRO_23;
DENDRO_120 *= DENDRO_55;
DENDRO_121 = 0.5;

// initialize reduction 
DENDRO_122 = 1;
DENDRO_122 *= DENDRO_121;
DENDRO_122 *= DENDRO_45;
DENDRO_122 *= DENDRO_52;

// initialize reduction 
DENDRO_123 = 1;
DENDRO_123 *= DENDRO_15;
DENDRO_123 *= DENDRO_55;
DENDRO_124 = 1.0;

// initialize reduction 
DENDRO_125 = 1;
DENDRO_125 *= DENDRO_124;
DENDRO_125 *= DENDRO_34;
DENDRO_125 *= DENDRO_43;
DENDRO_124 = 0.5;

// initialize reduction 
DENDRO_126 = 1;
DENDRO_126 *= DENDRO_124;
DENDRO_126 *= DENDRO_23;
DENDRO_126 *= DENDRO_39;
DENDRO_124 = 0.5;

// initialize reduction 
DENDRO_128 = 1;
DENDRO_128 *= DENDRO_124;
DENDRO_128 *= DENDRO_38;
DENDRO_128 *= DENDRO_5;

// initialize reduction 
DENDRO_124 = 1;
DENDRO_124 *= DENDRO_18;
DENDRO_124 *= DENDRO_33;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_83;
DENDRO_33 *= DENDRO_23;
DENDRO_33 *= DENDRO_57;
DENDRO_57 = 0.5;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_57;
DENDRO_83 *= DENDRO_34;
DENDRO_83 *= DENDRO_52;
DENDRO_57 = 1.0;

// initialize reduction 
DENDRO_130 = 1;
DENDRO_130 *= DENDRO_57;
DENDRO_130 *= DENDRO_51;
DENDRO_130 *= DENDRO_39;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_21;
DENDRO_51 *= DENDRO_47;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_109;
DENDRO_47 *= DENDRO_30;
DENDRO_47 *= DENDRO_53;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_91;
DENDRO_53 *= DENDRO_42;
DENDRO_53 *= DENDRO_104;
DENDRO_42 = 1.0;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_42;
DENDRO_57 *= DENDRO_38;
DENDRO_57 *= DENDRO_37;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_8;
DENDRO_38 *= DENDRO_98;
DENDRO_38 *= DENDRO_5;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_8;
DENDRO_42 *= DENDRO_100;
DENDRO_42 *= DENDRO_39;

// initialize reduction 
DENDRO_91 = 1;
DENDRO_91 *= DENDRO_8;
DENDRO_91 *= DENDRO_102;
DENDRO_91 *= DENDRO_52;

// initialize reduction 
DENDRO_109 = 0;
DENDRO_109 += DENDRO_32;
DENDRO_109 += DENDRO_19;
DENDRO_109 += DENDRO_6;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_48;
DENDRO_6 += DENDRO_41;
DENDRO_6 += DENDRO_35;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_162;
DENDRO_19 += DENDRO_29;
DENDRO_19 += DENDRO_49;

// initialize reduction 
DENDRO_29 = 0;
DENDRO_29 += DENDRO_60;
DENDRO_29 += DENDRO_58;
DENDRO_29 += DENDRO_56;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_70;
DENDRO_32 += DENDRO_69;
DENDRO_32 += DENDRO_63;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_161;
DENDRO_35 += DENDRO_139;
DENDRO_35 += DENDRO_71;

// initialize reduction 
DENDRO_41 = 0;
DENDRO_41 += DENDRO_78;
DENDRO_41 += DENDRO_77;
DENDRO_41 += DENDRO_74;

// initialize reduction 
DENDRO_48 = 0;
DENDRO_48 += DENDRO_140;
DENDRO_48 += DENDRO_79;
DENDRO_48 += DENDRO_80;

// initialize reduction 
DENDRO_49 = 0;
DENDRO_49 += DENDRO_85;
DENDRO_49 += DENDRO_82;
DENDRO_49 += DENDRO_84;

// initialize reduction 
DENDRO_56 = 0;
DENDRO_56 += DENDRO_87;
DENDRO_56 += DENDRO_72;
DENDRO_56 += DENDRO_86;

// initialize reduction 
DENDRO_58 = 0;
DENDRO_58 += DENDRO_93;
DENDRO_58 += DENDRO_88;
DENDRO_58 += DENDRO_90;

// initialize reduction 
DENDRO_60 = 0;
DENDRO_60 += DENDRO_96;
DENDRO_60 += DENDRO_95;

// initialize reduction 
DENDRO_69 = 0;
DENDRO_69 += DENDRO_175;
DENDRO_69 += DENDRO_89;
DENDRO_69 += DENDRO_94;

// initialize reduction 
DENDRO_70 = 0;
DENDRO_70 += DENDRO_131;
DENDRO_70 += DENDRO_75;
DENDRO_70 += DENDRO_92;

// initialize reduction 
DENDRO_71 = 0;
DENDRO_71 += DENDRO_106;
DENDRO_71 += DENDRO_101;
DENDRO_71 += DENDRO_99;

// initialize reduction 
DENDRO_72 = 0;
DENDRO_72 += DENDRO_112;
DENDRO_72 += DENDRO_24;
DENDRO_72 += DENDRO_111;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_174;
DENDRO_24 += DENDRO_114;
DENDRO_24 += DENDRO_22;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_117;
DENDRO_22 += DENDRO_108;
DENDRO_22 += DENDRO_116;

// initialize reduction 
DENDRO_74 = 0;
DENDRO_74 += DENDRO_1;
DENDRO_74 += DENDRO_118;
DENDRO_74 += DENDRO_119;

// initialize reduction 
DENDRO_75 = 0;
DENDRO_75 += DENDRO_123;
DENDRO_75 += DENDRO_122;
DENDRO_75 += DENDRO_120;

// initialize reduction 
DENDRO_77 = 0;
DENDRO_77 += DENDRO_126;
DENDRO_77 += DENDRO_125;

// initialize reduction 
DENDRO_78 = 0;
DENDRO_78 += DENDRO_110;
DENDRO_78 += DENDRO_124;
DENDRO_78 += DENDRO_128;

// initialize reduction 
DENDRO_79 = 0;
DENDRO_79 += DENDRO_115;
DENDRO_79 += DENDRO_83;
DENDRO_79 += DENDRO_33;

// initialize reduction 
DENDRO_33 = 0;
DENDRO_33 += DENDRO_51;
DENDRO_33 += DENDRO_130;

// initialize reduction 
DENDRO_51 = 0;
DENDRO_51 += DENDRO_73;
DENDRO_51 += DENDRO_113;
DENDRO_51 += DENDRO_47;

// initialize reduction 
DENDRO_47 = 0;
DENDRO_47 += DENDRO_150;
DENDRO_47 += DENDRO_81;
DENDRO_47 += DENDRO_53;

// initialize reduction 
DENDRO_53 = 0;
DENDRO_53 += DENDRO_68;
DENDRO_53 += DENDRO_57;

// initialize reduction 
DENDRO_57 = 0;
DENDRO_57 += DENDRO_61;
DENDRO_57 += DENDRO_91;
DENDRO_57 += DENDRO_42;
DENDRO_57 += DENDRO_38;
DENDRO_38 = -1.0;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_38;
DENDRO_38 = grad2_1_2_gt1;
DENDRO_42 *= DENDRO_38;
DENDRO_42 *= DENDRO_16;
DENDRO_42 *= DENDRO_7;
DENDRO_38 = -1.0;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_38;
DENDRO_38 = grad2_0_2_gt1;
DENDRO_61 *= DENDRO_38;
DENDRO_61 *= DENDRO_16;
DENDRO_61 *= DENDRO_54;
DENDRO_38 = -1.0;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_38;
DENDRO_38 = grad2_0_1_gt1;
DENDRO_68 *= DENDRO_38;
DENDRO_68 *= DENDRO_16;
DENDRO_68 *= DENDRO_26;
DENDRO_38 = -0.5;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_38;
DENDRO_38 = grad2_2_2_gt1;
DENDRO_73 *= DENDRO_38;
DENDRO_73 *= DENDRO_16;
DENDRO_73 *= DENDRO_17;
DENDRO_38 = -0.5;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_38;
DENDRO_38 = grad2_1_1_gt1;
DENDRO_80 *= DENDRO_38;
DENDRO_80 *= DENDRO_16;
DENDRO_80 *= DENDRO_27;
DENDRO_38 = -0.5;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_38;
DENDRO_38 = grad2_0_0_gt1;
DENDRO_81 *= DENDRO_38;
DENDRO_81 *= DENDRO_16;
DENDRO_81 *= DENDRO_3;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_127;
DENDRO_38 *= DENDRO_102;
DENDRO_38 *= DENDRO_100;
DENDRO_38 *= DENDRO_105;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_62;
DENDRO_82 *= DENDRO_0;
DENDRO_82 *= DENDRO_97;
DENDRO_82 *= DENDRO_144;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_16;
DENDRO_83 *= DENDRO_54;
DENDRO_83 *= DENDRO_109;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_16;
DENDRO_84 *= DENDRO_54;
DENDRO_84 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_16;
DENDRO_6 *= DENDRO_54;
DENDRO_6 *= DENDRO_19;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_16;
DENDRO_19 *= DENDRO_54;
DENDRO_19 *= DENDRO_29;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_16;
DENDRO_29 *= DENDRO_54;
DENDRO_29 *= DENDRO_32;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_16;
DENDRO_32 *= DENDRO_54;
DENDRO_32 *= DENDRO_35;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_16;
DENDRO_35 *= DENDRO_26;
DENDRO_35 *= DENDRO_41;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_16;
DENDRO_41 *= DENDRO_26;
DENDRO_41 *= DENDRO_48;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_16;
DENDRO_48 *= DENDRO_26;
DENDRO_48 *= DENDRO_49;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_16;
DENDRO_49 *= DENDRO_26;
DENDRO_49 *= DENDRO_56;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_16;
DENDRO_56 *= DENDRO_26;
DENDRO_56 *= DENDRO_58;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_16;
DENDRO_58 *= DENDRO_26;
DENDRO_58 *= DENDRO_60;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_16;
DENDRO_60 *= DENDRO_7;
DENDRO_60 *= DENDRO_69;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_16;
DENDRO_69 *= DENDRO_7;
DENDRO_69 *= DENDRO_70;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_16;
DENDRO_70 *= DENDRO_7;
DENDRO_70 *= DENDRO_71;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_16;
DENDRO_71 *= DENDRO_7;
DENDRO_71 *= DENDRO_72;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_16;
DENDRO_72 *= DENDRO_7;
DENDRO_72 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_16;
DENDRO_24 *= DENDRO_7;
DENDRO_24 *= DENDRO_22;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_16;
DENDRO_22 *= DENDRO_3;
DENDRO_22 *= DENDRO_74;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_16;
DENDRO_74 *= DENDRO_3;
DENDRO_74 *= DENDRO_75;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_16;
DENDRO_75 *= DENDRO_3;
DENDRO_75 *= DENDRO_77;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_16;
DENDRO_77 *= DENDRO_27;
DENDRO_77 *= DENDRO_78;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_16;
DENDRO_78 *= DENDRO_27;
DENDRO_78 *= DENDRO_79;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_16;
DENDRO_79 *= DENDRO_27;
DENDRO_79 *= DENDRO_33;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_16;
DENDRO_33 *= DENDRO_17;
DENDRO_33 *= DENDRO_51;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_16;
DENDRO_51 *= DENDRO_17;
DENDRO_51 *= DENDRO_47;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_16;
DENDRO_47 *= DENDRO_17;
DENDRO_47 *= DENDRO_53;
DENDRO_53 = 0.5;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_53;
DENDRO_85 *= DENDRO_14;
DENDRO_85 *= DENDRO_65;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_14;
DENDRO_53 *= DENDRO_36;
DENDRO_53 *= DENDRO_66;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_14;
DENDRO_36 *= DENDRO_40;
DENDRO_36 *= DENDRO_2;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_14;
DENDRO_40 *= DENDRO_146;
DENDRO_40 *= DENDRO_0;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_14;
DENDRO_86 *= DENDRO_148;
DENDRO_86 *= DENDRO_10;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_14;
DENDRO_87 *= DENDRO_28;
DENDRO_87 *= DENDRO_67;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_14;
DENDRO_28 *= DENDRO_9;
DENDRO_28 *= DENDRO_4;
DENDRO_4 = 0.5;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_4;
DENDRO_14 *= DENDRO_12;
DENDRO_14 *= DENDRO_11;
DENDRO_4 = 0.5;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_4;
DENDRO_11 *= DENDRO_134;
DENDRO_11 *= DENDRO_0;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_62;
DENDRO_4 *= DENDRO_97;
DENDRO_4 *= DENDRO_57;

// initialize reduction 
DENDRO_57 = 0;
DENDRO_57 += DENDRO_4;
DENDRO_57 += DENDRO_11;
DENDRO_57 += DENDRO_14;
DENDRO_57 += DENDRO_28;
DENDRO_57 += DENDRO_87;
DENDRO_57 += DENDRO_86;
DENDRO_57 += DENDRO_40;
DENDRO_57 += DENDRO_36;
DENDRO_57 += DENDRO_53;
DENDRO_57 += DENDRO_85;
DENDRO_57 += DENDRO_47;
DENDRO_57 += DENDRO_51;
DENDRO_57 += DENDRO_33;
DENDRO_57 += DENDRO_79;
DENDRO_57 += DENDRO_78;
DENDRO_57 += DENDRO_77;
DENDRO_57 += DENDRO_75;
DENDRO_57 += DENDRO_74;
DENDRO_57 += DENDRO_22;
DENDRO_57 += DENDRO_24;
DENDRO_57 += DENDRO_72;
DENDRO_57 += DENDRO_71;
DENDRO_57 += DENDRO_70;
DENDRO_57 += DENDRO_69;
DENDRO_57 += DENDRO_60;
DENDRO_57 += DENDRO_58;
DENDRO_57 += DENDRO_56;
DENDRO_57 += DENDRO_49;
DENDRO_57 += DENDRO_48;
DENDRO_57 += DENDRO_41;
DENDRO_57 += DENDRO_35;
DENDRO_57 += DENDRO_32;
DENDRO_57 += DENDRO_29;
DENDRO_57 += DENDRO_19;
DENDRO_57 += DENDRO_6;
DENDRO_57 += DENDRO_84;
DENDRO_57 += DENDRO_83;
DENDRO_57 += DENDRO_82;
DENDRO_57 += DENDRO_38;
DENDRO_57 += DENDRO_81;
DENDRO_57 += DENDRO_80;
DENDRO_57 += DENDRO_73;
DENDRO_57 += DENDRO_68;
DENDRO_57 += DENDRO_61;
DENDRO_57 += DENDRO_42;
double DENDRO_RIJ01 = DENDRO_57;
DENDRO_4 = 1.0;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_4;
DENDRO_6 *= DENDRO_45;
DENDRO_6 *= DENDRO_50;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_165;
DENDRO_4 *= DENDRO_30;
DENDRO_4 *= DENDRO_55;
DENDRO_11 = 1.0;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_11;
DENDRO_14 *= DENDRO_30;
DENDRO_14 *= DENDRO_104;
DENDRO_11 = 0.5;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_11;
DENDRO_19 *= DENDRO_46;
DENDRO_19 *= DENDRO_31;
DENDRO_11 = 1.0;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_11;
DENDRO_22 *= DENDRO_30;
DENDRO_22 *= DENDRO_55;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_167;
DENDRO_11 *= DENDRO_45;
DENDRO_11 *= DENDRO_50;
DENDRO_24 = 1.0;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_24;
DENDRO_28 *= DENDRO_23;
DENDRO_28 *= DENDRO_37;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_59;
DENDRO_24 *= DENDRO_25;
DENDRO_24 *= DENDRO_31;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_20;
DENDRO_29 *= DENDRO_104;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_59;
DENDRO_32 *= DENDRO_18;
DENDRO_32 *= DENDRO_43;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_15;
DENDRO_33 *= DENDRO_37;
DENDRO_35 = 1.0;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_35;
DENDRO_36 *= DENDRO_45;
DENDRO_36 *= DENDRO_52;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_13;
DENDRO_35 *= DENDRO_23;
DENDRO_35 *= DENDRO_55;
DENDRO_13 = 1.0;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_13;
DENDRO_38 *= DENDRO_23;
DENDRO_38 *= DENDRO_39;
DENDRO_13 = 0.5;

// initialize reduction 
DENDRO_40 = 1;
DENDRO_40 *= DENDRO_13;
DENDRO_40 *= DENDRO_34;
DENDRO_40 *= DENDRO_43;
DENDRO_13 = 1.0;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_13;
DENDRO_41 *= DENDRO_23;
DENDRO_41 *= DENDRO_55;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_121;
DENDRO_13 *= DENDRO_45;
DENDRO_13 *= DENDRO_52;
DENDRO_42 = 1.0;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_42;
DENDRO_47 *= DENDRO_30;
DENDRO_47 *= DENDRO_5;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_59;
DENDRO_42 *= DENDRO_18;
DENDRO_42 *= DENDRO_31;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_20;
DENDRO_48 *= DENDRO_5;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_59;
DENDRO_49 *= DENDRO_21;
DENDRO_49 *= DENDRO_43;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_15;
DENDRO_51 *= DENDRO_39;
DENDRO_53 = 1.0;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_53;
DENDRO_56 *= DENDRO_23;
DENDRO_56 *= DENDRO_50;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_158;
DENDRO_53 *= DENDRO_30;
DENDRO_53 *= DENDRO_52;
DENDRO_57 = 1.0;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_57;
DENDRO_58 *= DENDRO_30;
DENDRO_58 *= DENDRO_52;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_129;
DENDRO_57 *= DENDRO_23;
DENDRO_57 *= DENDRO_50;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_107;
DENDRO_60 *= DENDRO_46;
DENDRO_60 *= DENDRO_5;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_59;
DENDRO_61 *= DENDRO_18;
DENDRO_61 *= DENDRO_104;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_152;
DENDRO_68 *= DENDRO_34;
DENDRO_68 *= DENDRO_37;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_59;
DENDRO_69 *= DENDRO_18;
DENDRO_69 *= DENDRO_39;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_59;
DENDRO_70 *= DENDRO_21;
DENDRO_70 *= DENDRO_37;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_59;
DENDRO_71 *= DENDRO_25;
DENDRO_71 *= DENDRO_5;
DENDRO_72 = 1.0;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_72;
DENDRO_73 *= DENDRO_30;
DENDRO_73 *= DENDRO_31;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_20;
DENDRO_72 *= DENDRO_31;
DENDRO_20 = 1.0;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_20;
DENDRO_74 *= DENDRO_23;
DENDRO_74 *= DENDRO_43;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_15;
DENDRO_20 *= DENDRO_43;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_76;
DENDRO_15 *= DENDRO_34;
DENDRO_15 *= DENDRO_39;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_59;
DENDRO_34 *= DENDRO_21;
DENDRO_34 *= DENDRO_39;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_59;
DENDRO_21 *= DENDRO_18;
DENDRO_21 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_44;
DENDRO_5 *= DENDRO_46;
DENDRO_5 *= DENDRO_104;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_59;
DENDRO_39 *= DENDRO_25;
DENDRO_39 *= DENDRO_104;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_59;
DENDRO_25 *= DENDRO_18;
DENDRO_25 *= DENDRO_37;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_8;
DENDRO_18 *= DENDRO_98;
DENDRO_18 *= DENDRO_31;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_8;
DENDRO_31 *= DENDRO_100;
DENDRO_31 *= DENDRO_43;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_8;
DENDRO_37 *= DENDRO_102;
DENDRO_37 *= DENDRO_55;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_4;
DENDRO_8 += DENDRO_6;

// initialize reduction 
DENDRO_4 = 0;
DENDRO_4 += DENDRO_19;
DENDRO_4 += DENDRO_14;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_11;
DENDRO_6 += DENDRO_22;

// initialize reduction 
DENDRO_11 = 0;
DENDRO_11 += DENDRO_159;
DENDRO_11 += DENDRO_28;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_29;
DENDRO_14 += DENDRO_24;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_33;
DENDRO_19 += DENDRO_32;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_35;
DENDRO_22 += DENDRO_36;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_40;
DENDRO_24 += DENDRO_38;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_13;
DENDRO_28 += DENDRO_41;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_1;
DENDRO_13 += DENDRO_47;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_48;
DENDRO_1 += DENDRO_42;

// initialize reduction 
DENDRO_29 = 0;
DENDRO_29 += DENDRO_51;
DENDRO_29 += DENDRO_49;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_53;
DENDRO_32 += DENDRO_56;

// initialize reduction 
DENDRO_33 = 0;
DENDRO_33 += DENDRO_57;
DENDRO_33 += DENDRO_58;

// initialize reduction 
DENDRO_35 = 0;
DENDRO_35 += DENDRO_61;
DENDRO_35 += DENDRO_60;

// initialize reduction 
DENDRO_36 = 0;
DENDRO_36 += DENDRO_69;
DENDRO_36 += DENDRO_68;

// initialize reduction 
DENDRO_38 = 0;
DENDRO_38 += DENDRO_156;
DENDRO_38 += DENDRO_70;

// initialize reduction 
DENDRO_40 = 0;
DENDRO_40 += DENDRO_63;
DENDRO_40 += DENDRO_71;

// initialize reduction 
DENDRO_41 = 0;
DENDRO_41 += DENDRO_72;
DENDRO_41 += DENDRO_73;

// initialize reduction 
DENDRO_42 = 0;
DENDRO_42 += DENDRO_20;
DENDRO_42 += DENDRO_74;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_34;
DENDRO_20 += DENDRO_15;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_140;
DENDRO_15 += DENDRO_21;

// initialize reduction 
DENDRO_21 = 0;
DENDRO_21 += DENDRO_39;
DENDRO_21 += DENDRO_5;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_132;
DENDRO_5 += DENDRO_25;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_64;
DENDRO_25 += DENDRO_37;
DENDRO_25 += DENDRO_31;
DENDRO_25 += DENDRO_18;
DENDRO_18 = 1.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_18;
DENDRO_31 *= DENDRO_30;
DENDRO_31 *= DENDRO_16;
DENDRO_31 *= DENDRO_17;
DENDRO_31 *= DENDRO_50;
DENDRO_18 = 1.5;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_18;
DENDRO_34 *= DENDRO_23;
DENDRO_34 *= DENDRO_16;
DENDRO_34 *= DENDRO_27;
DENDRO_34 *= DENDRO_52;
DENDRO_18 = 1.5;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_18;
DENDRO_37 *= DENDRO_45;
DENDRO_37 *= DENDRO_16;
DENDRO_37 *= DENDRO_3;
DENDRO_37 *= DENDRO_55;
DENDRO_18 = -1.0;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_18;
DENDRO_18 = grad2_1_2_gt0;
DENDRO_39 *= DENDRO_18;
DENDRO_39 *= DENDRO_16;
DENDRO_39 *= DENDRO_7;
DENDRO_18 = -1.0;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_18;
DENDRO_18 = grad2_0_2_gt0;
DENDRO_43 *= DENDRO_18;
DENDRO_43 *= DENDRO_16;
DENDRO_43 *= DENDRO_54;
DENDRO_18 = -1.0;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_18;
DENDRO_18 = grad2_0_1_gt0;
DENDRO_44 *= DENDRO_18;
DENDRO_44 *= DENDRO_16;
DENDRO_44 *= DENDRO_26;
DENDRO_18 = -0.5;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_18;
DENDRO_18 = grad2_2_2_gt0;
DENDRO_46 *= DENDRO_18;
DENDRO_46 *= DENDRO_16;
DENDRO_46 *= DENDRO_17;
DENDRO_18 = -0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_18;
DENDRO_18 = grad2_1_1_gt0;
DENDRO_47 *= DENDRO_18;
DENDRO_47 *= DENDRO_16;
DENDRO_47 *= DENDRO_27;
DENDRO_18 = -0.5;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_18;
DENDRO_18 = grad2_0_0_gt0;
DENDRO_48 *= DENDRO_18;
DENDRO_48 *= DENDRO_16;
DENDRO_48 *= DENDRO_3;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_62;
DENDRO_18 *= DENDRO_10;
DENDRO_18 *= DENDRO_97;
DENDRO_18 *= DENDRO_144;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_16;
DENDRO_49 *= DENDRO_54;
DENDRO_49 *= DENDRO_8;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_16;
DENDRO_8 *= DENDRO_54;
DENDRO_8 *= DENDRO_4;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_16;
DENDRO_4 *= DENDRO_54;
DENDRO_4 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_16;
DENDRO_6 *= DENDRO_54;
DENDRO_6 *= DENDRO_11;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_16;
DENDRO_11 *= DENDRO_54;
DENDRO_11 *= DENDRO_14;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_16;
DENDRO_14 *= DENDRO_54;
DENDRO_14 *= DENDRO_19;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_16;
DENDRO_19 *= DENDRO_26;
DENDRO_19 *= DENDRO_22;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_16;
DENDRO_22 *= DENDRO_26;
DENDRO_22 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_16;
DENDRO_24 *= DENDRO_26;
DENDRO_24 *= DENDRO_28;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_16;
DENDRO_28 *= DENDRO_26;
DENDRO_28 *= DENDRO_13;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_16;
DENDRO_13 *= DENDRO_26;
DENDRO_13 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_16;
DENDRO_1 *= DENDRO_26;
DENDRO_1 *= DENDRO_29;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_16;
DENDRO_26 *= DENDRO_7;
DENDRO_26 *= DENDRO_32;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_16;
DENDRO_29 *= DENDRO_7;
DENDRO_29 *= DENDRO_33;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_16;
DENDRO_32 *= DENDRO_7;
DENDRO_32 *= DENDRO_35;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_16;
DENDRO_33 *= DENDRO_7;
DENDRO_33 *= DENDRO_36;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_16;
DENDRO_35 *= DENDRO_7;
DENDRO_35 *= DENDRO_38;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_16;
DENDRO_36 *= DENDRO_7;
DENDRO_36 *= DENDRO_40;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_16;
DENDRO_7 *= DENDRO_3;
DENDRO_7 *= DENDRO_41;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_16;
DENDRO_38 *= DENDRO_3;
DENDRO_38 *= DENDRO_42;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_16;
DENDRO_3 *= DENDRO_27;
DENDRO_3 *= DENDRO_20;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_16;
DENDRO_20 *= DENDRO_27;
DENDRO_20 *= DENDRO_15;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_16;
DENDRO_15 *= DENDRO_17;
DENDRO_15 *= DENDRO_21;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_16;
DENDRO_21 *= DENDRO_17;
DENDRO_21 *= DENDRO_5;
DENDRO_5 = 1.0;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_5;
DENDRO_16 *= DENDRO_9;
DENDRO_16 *= DENDRO_2;
DENDRO_2 = 1.0;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_2;
DENDRO_5 *= DENDRO_12;
DENDRO_5 *= DENDRO_0;
DENDRO_0 = 1.0;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_0;
DENDRO_2 *= DENDRO_134;
DENDRO_2 *= DENDRO_10;
DENDRO_0 = 0.5;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_0;
DENDRO_9 *= DENDRO_30;
DENDRO_9 *= DENDRO_65;
DENDRO_0 = 0.5;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_0;
DENDRO_10 *= DENDRO_23;
DENDRO_10 *= DENDRO_66;
DENDRO_0 = 0.5;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_0;
DENDRO_12 *= DENDRO_45;
DENDRO_12 *= DENDRO_67;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_127;
DENDRO_0 *= DENDRO_105;
DENDRO_0 *= DENDRO_103;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_62;
DENDRO_17 *= DENDRO_97;
DENDRO_17 *= DENDRO_25;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_17;
DENDRO_23 += DENDRO_0;
DENDRO_23 += DENDRO_12;
DENDRO_23 += DENDRO_10;
DENDRO_23 += DENDRO_9;
DENDRO_23 += DENDRO_2;
DENDRO_23 += DENDRO_5;
DENDRO_23 += DENDRO_16;
DENDRO_23 += DENDRO_21;
DENDRO_23 += DENDRO_15;
DENDRO_23 += DENDRO_20;
DENDRO_23 += DENDRO_3;
DENDRO_23 += DENDRO_38;
DENDRO_23 += DENDRO_7;
DENDRO_23 += DENDRO_36;
DENDRO_23 += DENDRO_35;
DENDRO_23 += DENDRO_33;
DENDRO_23 += DENDRO_32;
DENDRO_23 += DENDRO_29;
DENDRO_23 += DENDRO_26;
DENDRO_23 += DENDRO_1;
DENDRO_23 += DENDRO_13;
DENDRO_23 += DENDRO_28;
DENDRO_23 += DENDRO_24;
DENDRO_23 += DENDRO_22;
DENDRO_23 += DENDRO_19;
DENDRO_23 += DENDRO_14;
DENDRO_23 += DENDRO_11;
DENDRO_23 += DENDRO_6;
DENDRO_23 += DENDRO_4;
DENDRO_23 += DENDRO_8;
DENDRO_23 += DENDRO_49;
DENDRO_23 += DENDRO_18;
DENDRO_23 += DENDRO_48;
DENDRO_23 += DENDRO_47;
DENDRO_23 += DENDRO_46;
DENDRO_23 += DENDRO_44;
DENDRO_23 += DENDRO_43;
DENDRO_23 += DENDRO_39;
DENDRO_23 += DENDRO_37;
DENDRO_23 += DENDRO_34;
DENDRO_23 += DENDRO_31;
double DENDRO_RIJ00 = DENDRO_23;

DENDRO_RIJ0=DENDRO_RIJ00;
DENDRO_RIJ1=DENDRO_RIJ01;
DENDRO_RIJ2=DENDRO_RIJ02;
DENDRO_RIJ3=DENDRO_RIJ11;
DENDRO_RIJ4=DENDRO_RIJ12;
DENDRO_RIJ5=DENDRO_RIJ22;

}
device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , At0, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_At0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_At0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_At0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_At0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_At0=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_At0=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At0=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At0=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At0=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , At1, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_At1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_At1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_At1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_At1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_At1=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_At1=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At1=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At1=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At1=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , At2, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_At2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_At2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_At2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_At2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_At2=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_At2=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At2=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At2=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At2=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , At3, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_At3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_At3=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_At3=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_At3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_At3=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_At3=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At3=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At3=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At3=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At3=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , At4, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_At4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_At4=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_At4=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_At4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_At4=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_At4=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At4=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At4=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At4=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At4=Du[gidx];
__syncthreads();

device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , At5, blk);
device::__blk1_deriv644_xx<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_0_0_At5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_0_At5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_1_At5=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_0_2_At5=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_yy<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_1_1_At5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_1_At5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(DDu , Du, blk);
__syncthreads();
const DEVICE_REAL grad2_1_2_At5=DDu[gidx];
__syncthreads();
device::__blk1_deriv644_zz<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad2_2_2_At5=Du[gidx];
__syncthreads();
device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL grad_2_At5=Du[gidx];
__syncthreads();

device::__blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_0_At5=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_1_At5=Du[gidx];
__syncthreads();
device::__blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du , su, blk);
__syncthreads();
const DEVICE_REAL kograd_2_At5=Du[gidx];
__syncthreads();

//|V|= 523 |E| =1085
//topological generations = 18
// [0]
// traversal id=0 temp_var_count=89
{

double DENDRO_0;double DENDRO_1;double DENDRO_2;double DENDRO_3;double DENDRO_4;double DENDRO_5;double DENDRO_6;double DENDRO_7;double DENDRO_8;double DENDRO_9;double DENDRO_10;double DENDRO_11;double DENDRO_12;double DENDRO_13;double DENDRO_14;double DENDRO_15;double DENDRO_16;double DENDRO_17;double DENDRO_18;double DENDRO_19;double DENDRO_20;double DENDRO_21;double DENDRO_22;double DENDRO_23;double DENDRO_24;double DENDRO_25;double DENDRO_26;double DENDRO_27;double DENDRO_28;double DENDRO_29;double DENDRO_30;double DENDRO_31;double DENDRO_32;double DENDRO_33;double DENDRO_34;double DENDRO_35;double DENDRO_36;double DENDRO_37;double DENDRO_38;double DENDRO_39;double DENDRO_40;double DENDRO_41;double DENDRO_42;double DENDRO_43;double DENDRO_44;double DENDRO_45;double DENDRO_46;double DENDRO_47;double DENDRO_48;double DENDRO_49;double DENDRO_50;double DENDRO_51;double DENDRO_52;double DENDRO_53;double DENDRO_54;double DENDRO_55;double DENDRO_56;double DENDRO_57;double DENDRO_58;double DENDRO_59;double DENDRO_60;double DENDRO_61;double DENDRO_62;double DENDRO_63;double DENDRO_64;double DENDRO_65;double DENDRO_66;double DENDRO_67;double DENDRO_68;double DENDRO_69;double DENDRO_70;double DENDRO_71;double DENDRO_72;double DENDRO_73;double DENDRO_74;double DENDRO_75;double DENDRO_76;double DENDRO_77;double DENDRO_78;double DENDRO_79;double DENDRO_80;double DENDRO_81;double DENDRO_82;double DENDRO_83;double DENDRO_84;double DENDRO_85;double DENDRO_86;double DENDRO_87;double DENDRO_88;
DENDRO_0 = gt1[pp];

// int-power reduce pow(gt1[pp],2)
DENDRO_1 = DENDRO_0;
DENDRO_1 *= DENDRO_0;
DENDRO_2 = gt2[pp];

// int-power reduce pow(gt2[pp],2)
DENDRO_3 = DENDRO_2;
DENDRO_3 *= DENDRO_2;
DENDRO_4 = gt4[pp];

// int-power reduce pow(gt4[pp],2)
DENDRO_5 = DENDRO_4;
DENDRO_5 *= DENDRO_4;
DENDRO_6 = -2;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;
DENDRO_7 *= DENDRO_4;
DENDRO_8 = -1;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_8;
DENDRO_10 = gt0[pp];
DENDRO_9 *= DENDRO_10;
DENDRO_11 = gt3[pp];
DENDRO_9 *= DENDRO_11;
DENDRO_12 = gt5[pp];
DENDRO_9 *= DENDRO_12;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_12;
DENDRO_13 *= DENDRO_1;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_11;
DENDRO_14 *= DENDRO_3;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_10;
DENDRO_15 *= DENDRO_5;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_8;
DENDRO_16 *= DENDRO_10;
DENDRO_16 *= DENDRO_11;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_15;
DENDRO_17 += DENDRO_14;
DENDRO_17 += DENDRO_13;
DENDRO_17 += DENDRO_9;
DENDRO_17 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_8;
DENDRO_7 *= DENDRO_0;
DENDRO_7 *= DENDRO_2;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_10;
DENDRO_9 *= DENDRO_4;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_8;
DENDRO_13 *= DENDRO_0;
DENDRO_13 *= DENDRO_4;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_2;
DENDRO_14 *= DENDRO_11;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_8;
DENDRO_15 *= DENDRO_10;
DENDRO_15 *= DENDRO_12;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_8;
DENDRO_18 *= DENDRO_2;
DENDRO_18 *= DENDRO_4;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_0;
DENDRO_19 *= DENDRO_12;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_8;
DENDRO_20 *= DENDRO_11;
DENDRO_20 *= DENDRO_12;

// initialize reduction 
DENDRO_21 = 0;
DENDRO_21 += DENDRO_1;
DENDRO_21 += DENDRO_16;

// int-power reduce pow(-gt0[pp]*gt3[pp]*gt5[pp] + gt0[pp]*gt4[pp]**2 + gt1[pp]**2*gt5[pp] - 2*gt1[pp]*gt2[pp]*gt4[pp] + gt2[pp]**2*gt3[pp],-1)
DENDRO_1 = DENDRO_17;
DENDRO_1 = 1.0/DENDRO_1;

// initialize reduction 
DENDRO_16 = 0;
DENDRO_16 += DENDRO_9;
DENDRO_16 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_14;
DENDRO_7 += DENDRO_13;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_3;
DENDRO_9 += DENDRO_15;

// initialize reduction 
DENDRO_3 = 0;
DENDRO_3 += DENDRO_19;
DENDRO_3 += DENDRO_18;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_5;
DENDRO_13 += DENDRO_20;
DENDRO_5 = grad_2_chi;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_5;
DENDRO_14 *= DENDRO_1;
DENDRO_14 *= DENDRO_21;
DENDRO_15 = grad_1_chi;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_15;
DENDRO_17 *= DENDRO_1;
DENDRO_17 *= DENDRO_16;
DENDRO_18 = grad_0_chi;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_18;
DENDRO_19 *= DENDRO_1;
DENDRO_19 *= DENDRO_7;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_5;
DENDRO_20 *= DENDRO_1;
DENDRO_20 *= DENDRO_16;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_15;
DENDRO_22 *= DENDRO_1;
DENDRO_22 *= DENDRO_9;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_18;
DENDRO_23 *= DENDRO_1;
DENDRO_23 *= DENDRO_3;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_5;
DENDRO_24 *= DENDRO_1;
DENDRO_24 *= DENDRO_7;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_15;
DENDRO_25 *= DENDRO_1;
DENDRO_25 *= DENDRO_3;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_18;
DENDRO_26 *= DENDRO_1;
DENDRO_26 *= DENDRO_13;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_19;
DENDRO_27 += DENDRO_17;
DENDRO_27 += DENDRO_14;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_23;
DENDRO_14 += DENDRO_22;
DENDRO_14 += DENDRO_20;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_26;
DENDRO_17 += DENDRO_25;
DENDRO_17 += DENDRO_24;
DENDRO_19 = -0.5;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_19;
DENDRO_19 = grad_0_gt4;
DENDRO_20 *= DENDRO_19;
DENDRO_22 = 0.5;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_22;
DENDRO_22 = grad_2_gt1;
DENDRO_23 *= DENDRO_22;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_24;
DENDRO_24 = grad_1_gt2;
DENDRO_25 *= DENDRO_24;

// initialize reduction 
DENDRO_26 = 1;
DENDRO_26 *= DENDRO_8;
DENDRO_26 *= DENDRO_4;
DENDRO_26 *= DENDRO_27;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_8;
DENDRO_28 *= DENDRO_4;
DENDRO_28 *= DENDRO_14;
DENDRO_29 = -0.5;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_29;
DENDRO_30 *= DENDRO_24;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_24;
DENDRO_29 *= DENDRO_19;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_8;
DENDRO_19 *= DENDRO_2;
DENDRO_19 *= DENDRO_27;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_8;
DENDRO_24 *= DENDRO_2;
DENDRO_24 *= DENDRO_17;
DENDRO_31 = -0.5;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_31;
DENDRO_32 *= DENDRO_22;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_8;
DENDRO_22 *= DENDRO_0;
DENDRO_22 *= DENDRO_14;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_8;
DENDRO_31 *= DENDRO_0;
DENDRO_31 *= DENDRO_17;
DENDRO_33 = -0.5;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_33;
DENDRO_33 = grad_1_gt5;
DENDRO_34 *= DENDRO_33;
DENDRO_35 = 1.0;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_35;
DENDRO_35 = grad_2_gt4;
DENDRO_36 *= DENDRO_35;
DENDRO_35 = -0.5;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_35;
DENDRO_35 = grad_0_gt5;
DENDRO_37 *= DENDRO_35;
DENDRO_38 = 1.0;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_38;
DENDRO_38 = grad_2_gt2;
DENDRO_39 *= DENDRO_38;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_8;
DENDRO_38 *= DENDRO_12;
DENDRO_38 *= DENDRO_27;
DENDRO_40 = 2;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_40;
DENDRO_41 *= DENDRO_5;
DENDRO_42 = -0.5;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_42;
DENDRO_42 = grad_0_gt3;
DENDRO_43 *= DENDRO_42;
DENDRO_44 = 1.0;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_44;
DENDRO_44 = grad_1_gt1;
DENDRO_45 *= DENDRO_44;
DENDRO_44 = -0.5;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_44;
DENDRO_44 = grad_2_gt3;
DENDRO_46 *= DENDRO_44;
DENDRO_47 = 1.0;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_47;
DENDRO_47 = grad_1_gt4;
DENDRO_48 *= DENDRO_47;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_8;
DENDRO_47 *= DENDRO_11;
DENDRO_47 *= DENDRO_14;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_40;
DENDRO_49 *= DENDRO_15;
DENDRO_50 = -0.5;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_50;
DENDRO_50 = grad_1_gt0;
DENDRO_51 *= DENDRO_50;
DENDRO_52 = 1.0;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_52;
DENDRO_52 = grad_0_gt1;
DENDRO_53 *= DENDRO_52;
DENDRO_52 = -0.5;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_52;
DENDRO_52 = grad_2_gt0;
DENDRO_54 *= DENDRO_52;
DENDRO_55 = 1.0;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_55;
DENDRO_55 = grad_0_gt2;
DENDRO_56 *= DENDRO_55;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_8;
DENDRO_55 *= DENDRO_10;
DENDRO_55 *= DENDRO_17;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_40;
DENDRO_57 *= DENDRO_18;

// initialize reduction 
DENDRO_58 = 0;
DENDRO_58 += DENDRO_25;
DENDRO_58 += DENDRO_23;
DENDRO_58 += DENDRO_20;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_15;
DENDRO_20 += DENDRO_26;
DENDRO_26 = chi[pp];

// int-power reduce pow(chi[pp],-1)
DENDRO_59 = DENDRO_26;
DENDRO_59 = 1.0/DENDRO_59;

// initialize reduction 
DENDRO_60 = 0;
DENDRO_60 += DENDRO_5;
DENDRO_60 += DENDRO_28;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_29;
DENDRO_28 += DENDRO_23;
DENDRO_28 += DENDRO_30;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_18;
DENDRO_23 += DENDRO_19;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_5;
DENDRO_19 += DENDRO_24;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_29;
DENDRO_5 += DENDRO_25;
DENDRO_5 += DENDRO_32;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_18;
DENDRO_24 += DENDRO_22;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_15;
DENDRO_18 += DENDRO_31;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_36;
DENDRO_15 += DENDRO_34;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_39;
DENDRO_22 += DENDRO_37;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_41;
DENDRO_25 += DENDRO_38;

// initialize reduction 
DENDRO_29 = 0;
DENDRO_29 += DENDRO_45;
DENDRO_29 += DENDRO_43;

// initialize reduction 
DENDRO_30 = 0;
DENDRO_30 += DENDRO_48;
DENDRO_30 += DENDRO_46;

// initialize reduction 
DENDRO_31 = 0;
DENDRO_31 += DENDRO_49;
DENDRO_31 += DENDRO_47;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_53;
DENDRO_32 += DENDRO_51;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_56;
DENDRO_34 += DENDRO_54;

// initialize reduction 
DENDRO_36 = 0;
DENDRO_36 += DENDRO_57;
DENDRO_36 += DENDRO_55;
DENDRO_37 = 0.5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_37;
DENDRO_38 *= DENDRO_44;
DENDRO_38 *= DENDRO_1;
DENDRO_38 *= DENDRO_16;
DENDRO_37 = 0.5;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_37;
DENDRO_39 *= DENDRO_33;
DENDRO_39 *= DENDRO_1;
DENDRO_39 *= DENDRO_21;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_1;
DENDRO_37 *= DENDRO_7;
DENDRO_37 *= DENDRO_58;
DENDRO_41 = -0.5;

// initialize reduction 
DENDRO_43 = 1;
DENDRO_43 *= DENDRO_41;
DENDRO_43 *= DENDRO_59;
DENDRO_43 *= DENDRO_20;
DENDRO_20 = 0.5;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_20;
DENDRO_41 *= DENDRO_44;
DENDRO_41 *= DENDRO_1;
DENDRO_41 *= DENDRO_9;
DENDRO_20 = 0.5;

// initialize reduction 
DENDRO_45 = 1;
DENDRO_45 *= DENDRO_20;
DENDRO_45 *= DENDRO_33;
DENDRO_45 *= DENDRO_1;
DENDRO_45 *= DENDRO_16;

// initialize reduction 
DENDRO_20 = 1;
DENDRO_20 *= DENDRO_1;
DENDRO_20 *= DENDRO_3;
DENDRO_20 *= DENDRO_58;
DENDRO_46 = -0.5;

// initialize reduction 
DENDRO_47 = 1;
DENDRO_47 *= DENDRO_46;
DENDRO_47 *= DENDRO_59;
DENDRO_47 *= DENDRO_60;
DENDRO_46 = 0.5;

// initialize reduction 
DENDRO_48 = 1;
DENDRO_48 *= DENDRO_46;
DENDRO_48 *= DENDRO_4;
DENDRO_48 *= DENDRO_59;
DENDRO_48 *= DENDRO_17;
DENDRO_46 = 0.5;

// initialize reduction 
DENDRO_49 = 1;
DENDRO_49 *= DENDRO_46;
DENDRO_49 *= DENDRO_44;
DENDRO_49 *= DENDRO_1;
DENDRO_49 *= DENDRO_3;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_44;
DENDRO_46 *= DENDRO_33;
DENDRO_46 *= DENDRO_1;
DENDRO_46 *= DENDRO_7;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_1;
DENDRO_33 *= DENDRO_13;
DENDRO_33 *= DENDRO_58;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_44;
DENDRO_51 *= DENDRO_52;
DENDRO_51 *= DENDRO_1;
DENDRO_51 *= DENDRO_7;
DENDRO_44 = 0.5;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_44;
DENDRO_53 *= DENDRO_35;
DENDRO_53 *= DENDRO_1;
DENDRO_53 *= DENDRO_21;

// initialize reduction 
DENDRO_44 = 1;
DENDRO_44 *= DENDRO_1;
DENDRO_44 *= DENDRO_16;
DENDRO_44 *= DENDRO_28;
DENDRO_54 = -0.5;

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_54;
DENDRO_55 *= DENDRO_59;
DENDRO_55 *= DENDRO_23;
DENDRO_23 = 0.5;

// initialize reduction 
DENDRO_54 = 1;
DENDRO_54 *= DENDRO_23;
DENDRO_54 *= DENDRO_2;
DENDRO_54 *= DENDRO_59;
DENDRO_54 *= DENDRO_14;
DENDRO_23 = 0.5;

// initialize reduction 
DENDRO_56 = 1;
DENDRO_56 *= DENDRO_23;
DENDRO_56 *= DENDRO_52;
DENDRO_56 *= DENDRO_1;
DENDRO_56 *= DENDRO_3;
DENDRO_23 = 0.5;

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_23;
DENDRO_57 *= DENDRO_35;
DENDRO_57 *= DENDRO_1;
DENDRO_57 *= DENDRO_16;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_1;
DENDRO_23 *= DENDRO_9;
DENDRO_23 *= DENDRO_28;
DENDRO_58 = 0.5;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_58;
DENDRO_60 *= DENDRO_52;
DENDRO_60 *= DENDRO_1;
DENDRO_60 *= DENDRO_13;
DENDRO_52 = 0.5;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_52;
DENDRO_58 *= DENDRO_35;
DENDRO_58 *= DENDRO_1;
DENDRO_58 *= DENDRO_7;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_1;
DENDRO_35 *= DENDRO_3;
DENDRO_35 *= DENDRO_28;
DENDRO_28 = -0.5;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_28;
DENDRO_52 *= DENDRO_59;
DENDRO_52 *= DENDRO_19;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_19;
DENDRO_28 *= DENDRO_0;
DENDRO_28 *= DENDRO_59;
DENDRO_28 *= DENDRO_27;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_19;
DENDRO_61 *= DENDRO_50;
DENDRO_61 *= DENDRO_1;
DENDRO_61 *= DENDRO_7;
DENDRO_19 = 0.5;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_19;
DENDRO_62 *= DENDRO_42;
DENDRO_62 *= DENDRO_1;
DENDRO_62 *= DENDRO_16;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_1;
DENDRO_19 *= DENDRO_21;
DENDRO_19 *= DENDRO_5;
DENDRO_63 = 0.5;

// initialize reduction 
DENDRO_64 = 1;
DENDRO_64 *= DENDRO_63;
DENDRO_64 *= DENDRO_50;
DENDRO_64 *= DENDRO_1;
DENDRO_64 *= DENDRO_3;
DENDRO_63 = 0.5;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_63;
DENDRO_65 *= DENDRO_42;
DENDRO_65 *= DENDRO_1;
DENDRO_65 *= DENDRO_9;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_1;
DENDRO_63 *= DENDRO_16;
DENDRO_63 *= DENDRO_5;
DENDRO_66 = -0.5;

// initialize reduction 
DENDRO_67 = 1;
DENDRO_67 *= DENDRO_66;
DENDRO_67 *= DENDRO_59;
DENDRO_67 *= DENDRO_24;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_24;
DENDRO_66 *= DENDRO_50;
DENDRO_66 *= DENDRO_1;
DENDRO_66 *= DENDRO_13;
DENDRO_24 = 0.5;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_24;
DENDRO_50 *= DENDRO_42;
DENDRO_50 *= DENDRO_1;
DENDRO_50 *= DENDRO_3;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_1;
DENDRO_24 *= DENDRO_7;
DENDRO_24 *= DENDRO_5;
DENDRO_5 = -0.5;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_5;
DENDRO_42 *= DENDRO_59;
DENDRO_42 *= DENDRO_18;
DENDRO_5 = 0.5;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_5;
DENDRO_5 = grad_2_gt5;
DENDRO_18 *= DENDRO_5;
DENDRO_18 *= DENDRO_1;
DENDRO_18 *= DENDRO_21;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_1;
DENDRO_68 *= DENDRO_15;
DENDRO_68 *= DENDRO_16;

// initialize reduction 
DENDRO_69 = 1;
DENDRO_69 *= DENDRO_1;
DENDRO_69 *= DENDRO_22;
DENDRO_69 *= DENDRO_7;
DENDRO_70 = -0.5;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_70;
DENDRO_71 *= DENDRO_59;
DENDRO_71 *= DENDRO_25;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_25;
DENDRO_70 *= DENDRO_12;
DENDRO_70 *= DENDRO_59;
DENDRO_70 *= DENDRO_14;
DENDRO_25 = 0.5;

// initialize reduction 
DENDRO_72 = 1;
DENDRO_72 *= DENDRO_25;
DENDRO_72 *= DENDRO_5;
DENDRO_72 *= DENDRO_1;
DENDRO_72 *= DENDRO_16;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_1;
DENDRO_25 *= DENDRO_22;
DENDRO_25 *= DENDRO_3;

// initialize reduction 
DENDRO_73 = 1;
DENDRO_73 *= DENDRO_1;
DENDRO_73 *= DENDRO_9;
DENDRO_73 *= DENDRO_15;
DENDRO_74 = 0.5;

// initialize reduction 
DENDRO_75 = 1;
DENDRO_75 *= DENDRO_74;
DENDRO_75 *= DENDRO_12;
DENDRO_75 *= DENDRO_59;
DENDRO_75 *= DENDRO_17;
DENDRO_74 = 0.5;

// initialize reduction 
DENDRO_76 = 1;
DENDRO_76 *= DENDRO_74;
DENDRO_76 *= DENDRO_5;
DENDRO_76 *= DENDRO_1;
DENDRO_76 *= DENDRO_7;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_1;
DENDRO_5 *= DENDRO_15;
DENDRO_5 *= DENDRO_3;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_1;
DENDRO_15 *= DENDRO_13;
DENDRO_15 *= DENDRO_22;
DENDRO_22 = 0.5;

// initialize reduction 
DENDRO_74 = 1;
DENDRO_74 *= DENDRO_22;
DENDRO_74 *= DENDRO_11;
DENDRO_74 *= DENDRO_59;
DENDRO_74 *= DENDRO_27;
DENDRO_22 = 0.5;

// initialize reduction 
DENDRO_77 = 1;
DENDRO_77 *= DENDRO_22;
DENDRO_22 = grad_1_gt3;
DENDRO_77 *= DENDRO_22;
DENDRO_77 *= DENDRO_1;
DENDRO_77 *= DENDRO_16;

// initialize reduction 
DENDRO_78 = 1;
DENDRO_78 *= DENDRO_1;
DENDRO_78 *= DENDRO_29;
DENDRO_78 *= DENDRO_7;

// initialize reduction 
DENDRO_79 = 1;
DENDRO_79 *= DENDRO_1;
DENDRO_79 *= DENDRO_21;
DENDRO_79 *= DENDRO_30;
DENDRO_80 = 0.5;

// initialize reduction 
DENDRO_81 = 1;
DENDRO_81 *= DENDRO_80;
DENDRO_81 *= DENDRO_22;
DENDRO_81 *= DENDRO_1;
DENDRO_81 *= DENDRO_9;

// initialize reduction 
DENDRO_80 = 1;
DENDRO_80 *= DENDRO_1;
DENDRO_80 *= DENDRO_30;
DENDRO_80 *= DENDRO_16;

// initialize reduction 
DENDRO_82 = 1;
DENDRO_82 *= DENDRO_1;
DENDRO_82 *= DENDRO_29;
DENDRO_82 *= DENDRO_3;
DENDRO_83 = -0.5;

// initialize reduction 
DENDRO_84 = 1;
DENDRO_84 *= DENDRO_83;
DENDRO_84 *= DENDRO_59;
DENDRO_84 *= DENDRO_31;
DENDRO_31 = 0.5;

// initialize reduction 
DENDRO_83 = 1;
DENDRO_83 *= DENDRO_31;
DENDRO_83 *= DENDRO_11;
DENDRO_83 *= DENDRO_59;
DENDRO_83 *= DENDRO_17;
DENDRO_17 = 0.5;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_17;
DENDRO_31 *= DENDRO_22;
DENDRO_31 *= DENDRO_1;
DENDRO_31 *= DENDRO_3;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_1;
DENDRO_17 *= DENDRO_30;
DENDRO_17 *= DENDRO_7;

// initialize reduction 
DENDRO_22 = 1;
DENDRO_22 *= DENDRO_1;
DENDRO_22 *= DENDRO_13;
DENDRO_22 *= DENDRO_29;
DENDRO_29 = 0.5;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_29;
DENDRO_30 *= DENDRO_10;
DENDRO_30 *= DENDRO_59;
DENDRO_30 *= DENDRO_27;
DENDRO_27 = 0.5;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_27;
DENDRO_27 = grad_0_gt0;
DENDRO_29 *= DENDRO_27;
DENDRO_29 *= DENDRO_1;
DENDRO_29 *= DENDRO_7;

// initialize reduction 
DENDRO_85 = 1;
DENDRO_85 *= DENDRO_1;
DENDRO_85 *= DENDRO_32;
DENDRO_85 *= DENDRO_16;

// initialize reduction 
DENDRO_86 = 1;
DENDRO_86 *= DENDRO_1;
DENDRO_86 *= DENDRO_21;
DENDRO_86 *= DENDRO_34;
DENDRO_21 = 0.5;

// initialize reduction 
DENDRO_87 = 1;
DENDRO_87 *= DENDRO_21;
DENDRO_87 *= DENDRO_10;
DENDRO_87 *= DENDRO_59;
DENDRO_87 *= DENDRO_14;
DENDRO_14 = 0.5;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_14;
DENDRO_21 *= DENDRO_27;
DENDRO_21 *= DENDRO_1;
DENDRO_21 *= DENDRO_3;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_1;
DENDRO_14 *= DENDRO_34;
DENDRO_14 *= DENDRO_16;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_1;
DENDRO_16 *= DENDRO_9;
DENDRO_16 *= DENDRO_32;
DENDRO_9 = 0.5;

// initialize reduction 
DENDRO_88 = 1;
DENDRO_88 *= DENDRO_9;
DENDRO_88 *= DENDRO_27;
DENDRO_88 *= DENDRO_1;
DENDRO_88 *= DENDRO_13;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_1;
DENDRO_9 *= DENDRO_34;
DENDRO_9 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_1;
DENDRO_7 *= DENDRO_32;
DENDRO_7 *= DENDRO_3;
DENDRO_1 = -0.5;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_1;
DENDRO_3 *= DENDRO_59;
DENDRO_3 *= DENDRO_36;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_43;
DENDRO_1 += DENDRO_37;
DENDRO_1 += DENDRO_39;
DENDRO_1 += DENDRO_38;

// initialize reduction 
DENDRO_13 = 0;
DENDRO_13 += DENDRO_47;
DENDRO_13 += DENDRO_20;
DENDRO_13 += DENDRO_45;
DENDRO_13 += DENDRO_41;

// initialize reduction 
DENDRO_20 = 0;
DENDRO_20 += DENDRO_33;
DENDRO_20 += DENDRO_46;
DENDRO_20 += DENDRO_49;
DENDRO_20 += DENDRO_48;

// initialize reduction 
DENDRO_27 = 0;
DENDRO_27 += DENDRO_55;
DENDRO_27 += DENDRO_44;
DENDRO_27 += DENDRO_53;
DENDRO_27 += DENDRO_51;

// initialize reduction 
DENDRO_32 = 0;
DENDRO_32 += DENDRO_23;
DENDRO_32 += DENDRO_57;
DENDRO_32 += DENDRO_56;
DENDRO_32 += DENDRO_54;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_52;
DENDRO_23 += DENDRO_35;
DENDRO_23 += DENDRO_58;
DENDRO_23 += DENDRO_60;

// initialize reduction 
DENDRO_33 = 0;
DENDRO_33 += DENDRO_19;
DENDRO_33 += DENDRO_62;
DENDRO_33 += DENDRO_61;
DENDRO_33 += DENDRO_28;

// initialize reduction 
DENDRO_19 = 0;
DENDRO_19 += DENDRO_67;
DENDRO_19 += DENDRO_63;
DENDRO_19 += DENDRO_65;
DENDRO_19 += DENDRO_64;

// initialize reduction 
DENDRO_28 = 0;
DENDRO_28 += DENDRO_42;
DENDRO_28 += DENDRO_24;
DENDRO_28 += DENDRO_50;
DENDRO_28 += DENDRO_66;

// initialize reduction 
DENDRO_24 = 0;
DENDRO_24 += DENDRO_71;
DENDRO_24 += DENDRO_69;
DENDRO_24 += DENDRO_68;
DENDRO_24 += DENDRO_18;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_73;
DENDRO_18 += DENDRO_25;
DENDRO_18 += DENDRO_72;
DENDRO_18 += DENDRO_70;

// initialize reduction 
DENDRO_25 = 0;
DENDRO_25 += DENDRO_15;
DENDRO_25 += DENDRO_5;
DENDRO_25 += DENDRO_76;
DENDRO_25 += DENDRO_75;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_79;
DENDRO_5 += DENDRO_78;
DENDRO_5 += DENDRO_77;
DENDRO_5 += DENDRO_74;

// initialize reduction 
DENDRO_15 = 0;
DENDRO_15 += DENDRO_84;
DENDRO_15 += DENDRO_82;
DENDRO_15 += DENDRO_80;
DENDRO_15 += DENDRO_81;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_22;
DENDRO_34 += DENDRO_17;
DENDRO_34 += DENDRO_31;
DENDRO_34 += DENDRO_83;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_86;
DENDRO_17 += DENDRO_85;
DENDRO_17 += DENDRO_29;
DENDRO_17 += DENDRO_30;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_16;
DENDRO_22 += DENDRO_14;
DENDRO_22 += DENDRO_21;
DENDRO_22 += DENDRO_87;

// initialize reduction 
DENDRO_14 = 0;
DENDRO_14 += DENDRO_3;
DENDRO_14 += DENDRO_7;
DENDRO_14 += DENDRO_9;
DENDRO_14 += DENDRO_88;
DENDRO_3 = grad_2_alpha;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_3;
DENDRO_7 *= DENDRO_1;
DENDRO_1 = grad_1_alpha;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_1;
DENDRO_9 *= DENDRO_13;
DENDRO_13 = grad_0_alpha;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_13;
DENDRO_16 *= DENDRO_20;
DENDRO_20 = DENDRO_RIJ4;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_20;
DENDRO_20 = alpha[pp];
DENDRO_21 *= DENDRO_20;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_8;
DENDRO_30 = grad2_1_2_alpha;
DENDRO_29 *= DENDRO_30;

// initialize reduction 
DENDRO_30 = 1;
DENDRO_30 *= DENDRO_3;
DENDRO_30 *= DENDRO_27;

// initialize reduction 
DENDRO_27 = 1;
DENDRO_27 *= DENDRO_1;
DENDRO_27 *= DENDRO_32;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_13;
DENDRO_31 *= DENDRO_23;
DENDRO_23 = DENDRO_RIJ2;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_23;
DENDRO_32 *= DENDRO_20;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_8;
DENDRO_35 = grad2_0_2_alpha;
DENDRO_23 *= DENDRO_35;

// initialize reduction 
DENDRO_35 = 1;
DENDRO_35 *= DENDRO_3;
DENDRO_35 *= DENDRO_33;

// initialize reduction 
DENDRO_33 = 1;
DENDRO_33 *= DENDRO_1;
DENDRO_33 *= DENDRO_19;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_13;
DENDRO_19 *= DENDRO_28;
DENDRO_28 = DENDRO_RIJ1;

// initialize reduction 
DENDRO_36 = 1;
DENDRO_36 *= DENDRO_28;
DENDRO_36 *= DENDRO_20;

// initialize reduction 
DENDRO_28 = 1;
DENDRO_28 *= DENDRO_8;
DENDRO_37 = grad2_0_1_alpha;
DENDRO_28 *= DENDRO_37;

// initialize reduction 
DENDRO_37 = 1;
DENDRO_37 *= DENDRO_3;
DENDRO_37 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_1;
DENDRO_24 *= DENDRO_18;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_13;
DENDRO_18 *= DENDRO_25;
DENDRO_25 = DENDRO_RIJ5;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_25;
DENDRO_38 *= DENDRO_20;

// initialize reduction 
DENDRO_25 = 1;
DENDRO_25 *= DENDRO_8;
DENDRO_39 = grad2_2_2_alpha;
DENDRO_25 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_3;
DENDRO_39 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_1;
DENDRO_5 *= DENDRO_15;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_13;
DENDRO_15 *= DENDRO_34;
DENDRO_34 = DENDRO_RIJ3;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_34;
DENDRO_41 *= DENDRO_20;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_8;
DENDRO_42 = grad2_1_1_alpha;
DENDRO_34 *= DENDRO_42;

// initialize reduction 
DENDRO_42 = 1;
DENDRO_42 *= DENDRO_3;
DENDRO_42 *= DENDRO_17;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_1;
DENDRO_3 *= DENDRO_22;

// initialize reduction 
DENDRO_1 = 1;
DENDRO_1 *= DENDRO_13;
DENDRO_1 *= DENDRO_14;
DENDRO_13 = DENDRO_RIJ0;

// initialize reduction 
DENDRO_14 = 1;
DENDRO_14 *= DENDRO_13;
DENDRO_14 *= DENDRO_20;

// initialize reduction 
DENDRO_13 = 1;
DENDRO_13 *= DENDRO_8;
DENDRO_8 = grad2_0_0_alpha;
DENDRO_13 *= DENDRO_8;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_29;
DENDRO_8 += DENDRO_21;
DENDRO_8 += DENDRO_16;
DENDRO_8 += DENDRO_9;
DENDRO_8 += DENDRO_7;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_23;
DENDRO_17 += DENDRO_32;
DENDRO_17 += DENDRO_31;
DENDRO_17 += DENDRO_27;
DENDRO_17 += DENDRO_30;

// initialize reduction 
DENDRO_22 = 0;
DENDRO_22 += DENDRO_28;
DENDRO_22 += DENDRO_36;
DENDRO_22 += DENDRO_19;
DENDRO_22 += DENDRO_33;
DENDRO_22 += DENDRO_35;

// initialize reduction 
DENDRO_43 = 0;
DENDRO_43 += DENDRO_25;
DENDRO_43 += DENDRO_38;
DENDRO_43 += DENDRO_18;
DENDRO_43 += DENDRO_24;
DENDRO_43 += DENDRO_37;

// initialize reduction 
DENDRO_44 = 0;
DENDRO_44 += DENDRO_34;
DENDRO_44 += DENDRO_41;
DENDRO_44 += DENDRO_15;
DENDRO_44 += DENDRO_5;
DENDRO_44 += DENDRO_39;

// initialize reduction 
DENDRO_45 = 0;
DENDRO_45 += DENDRO_13;
DENDRO_45 += DENDRO_14;
DENDRO_45 += DENDRO_1;
DENDRO_45 += DENDRO_3;
DENDRO_45 += DENDRO_42;

// initialize reduction 
DENDRO_46 = 1;
DENDRO_46 *= DENDRO_40;
DENDRO_47 = DENDRO_igt4;
DENDRO_46 *= DENDRO_47;
DENDRO_46 *= DENDRO_8;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_40;
DENDRO_48 = DENDRO_igt2;
DENDRO_8 *= DENDRO_48;
DENDRO_8 *= DENDRO_17;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_40;
DENDRO_49 = DENDRO_igt1;
DENDRO_17 *= DENDRO_49;
DENDRO_17 *= DENDRO_22;
DENDRO_22 = DENDRO_igt5;

// initialize reduction 
DENDRO_50 = 1;
DENDRO_50 *= DENDRO_22;
DENDRO_50 *= DENDRO_43;
DENDRO_43 = DENDRO_igt3;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_43;
DENDRO_51 *= DENDRO_44;
DENDRO_44 = DENDRO_igt0;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_44;
DENDRO_52 *= DENDRO_45;
DENDRO_45 = At5[pp];

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_45;
DENDRO_53 *= DENDRO_22;
DENDRO_54 = At4[pp];

// initialize reduction 
DENDRO_55 = 1;
DENDRO_55 *= DENDRO_54;
DENDRO_55 *= DENDRO_47;
DENDRO_56 = At2[pp];

// initialize reduction 
DENDRO_57 = 1;
DENDRO_57 *= DENDRO_56;
DENDRO_57 *= DENDRO_48;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_45;
DENDRO_58 *= DENDRO_47;

// initialize reduction 
DENDRO_59 = 1;
DENDRO_59 *= DENDRO_54;
DENDRO_59 *= DENDRO_43;

// initialize reduction 
DENDRO_60 = 1;
DENDRO_60 *= DENDRO_56;
DENDRO_60 *= DENDRO_49;

// initialize reduction 
DENDRO_61 = 1;
DENDRO_61 *= DENDRO_45;
DENDRO_61 *= DENDRO_48;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_54;
DENDRO_62 *= DENDRO_49;

// initialize reduction 
DENDRO_63 = 1;
DENDRO_63 *= DENDRO_56;
DENDRO_63 *= DENDRO_44;

// initialize reduction 
DENDRO_64 = 0;
DENDRO_64 += DENDRO_52;
DENDRO_64 += DENDRO_51;
DENDRO_64 += DENDRO_50;
DENDRO_64 += DENDRO_17;
DENDRO_64 += DENDRO_8;
DENDRO_64 += DENDRO_46;

// initialize reduction 
DENDRO_8 = 0;
DENDRO_8 += DENDRO_57;
DENDRO_8 += DENDRO_55;
DENDRO_8 += DENDRO_53;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_60;
DENDRO_17 += DENDRO_59;
DENDRO_17 += DENDRO_58;

// initialize reduction 
DENDRO_46 = 0;
DENDRO_46 += DENDRO_63;
DENDRO_46 += DENDRO_62;
DENDRO_46 += DENDRO_61;
DENDRO_50 = -1.0/3.0;

// initialize reduction 
DENDRO_51 = 1;
DENDRO_51 *= DENDRO_50;
DENDRO_51 *= DENDRO_12;
DENDRO_51 *= DENDRO_64;

// initialize reduction 
DENDRO_12 = 1;
DENDRO_12 *= DENDRO_6;
DENDRO_12 *= DENDRO_45;
DENDRO_12 *= DENDRO_8;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_6;
DENDRO_52 *= DENDRO_54;
DENDRO_52 *= DENDRO_17;

// initialize reduction 
DENDRO_53 = 1;
DENDRO_53 *= DENDRO_6;
DENDRO_53 *= DENDRO_56;
DENDRO_53 *= DENDRO_46;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_45;
DENDRO_59 = K[pp];
DENDRO_58 *= DENDRO_59;

// initialize reduction 
DENDRO_60 = 0;
DENDRO_60 += DENDRO_25;
DENDRO_60 += DENDRO_38;
DENDRO_60 += DENDRO_18;
DENDRO_60 += DENDRO_24;
DENDRO_60 += DENDRO_37;
DENDRO_60 += DENDRO_51;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_58;
DENDRO_18 += DENDRO_53;
DENDRO_18 += DENDRO_52;
DENDRO_18 += DENDRO_12;
DENDRO_12 = 4.0/3.0;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_12;
DENDRO_24 *= DENDRO_45;
DENDRO_25 = grad_2_beta2;
DENDRO_24 *= DENDRO_25;
DENDRO_37 = -2.0/3.0;

// initialize reduction 
DENDRO_38 = 1;
DENDRO_38 *= DENDRO_37;
DENDRO_38 *= DENDRO_45;
DENDRO_51 = grad_1_beta1;
DENDRO_38 *= DENDRO_51;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_37;
DENDRO_52 *= DENDRO_45;
DENDRO_53 = grad_0_beta0;
DENDRO_52 *= DENDRO_53;

// initialize reduction 
DENDRO_58 = 1;
DENDRO_58 *= DENDRO_40;
DENDRO_58 *= DENDRO_54;
DENDRO_61 = grad_2_beta1;
DENDRO_58 *= DENDRO_61;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_40;
DENDRO_62 *= DENDRO_56;
DENDRO_63 = grad_2_beta0;
DENDRO_62 *= DENDRO_63;

// initialize reduction 
DENDRO_65 = 1;
DENDRO_65 *= DENDRO_26;
DENDRO_65 *= DENDRO_60;
DENDRO_60 = beta2[pp];

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_60;
DENDRO_67 = grad_2_At5;
DENDRO_66 *= DENDRO_67;
DENDRO_67 = beta1[pp];

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_67;
DENDRO_69 = grad_1_At5;
DENDRO_68 *= DENDRO_69;
DENDRO_69 = beta0[pp];

// initialize reduction 
DENDRO_70 = 1;
DENDRO_70 *= DENDRO_69;
DENDRO_71 = grad_0_At5;
DENDRO_70 *= DENDRO_71;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_20;
DENDRO_71 *= DENDRO_18;

// initialize reduction 
DENDRO_18 = 0;
DENDRO_18 += DENDRO_71;
DENDRO_18 += DENDRO_70;
DENDRO_18 += DENDRO_68;
DENDRO_18 += DENDRO_66;
DENDRO_18 += DENDRO_65;
DENDRO_18 += DENDRO_62;
DENDRO_18 += DENDRO_58;
DENDRO_18 += DENDRO_52;
DENDRO_18 += DENDRO_38;
DENDRO_18 += DENDRO_24;
At_rhs22[pp] = DENDRO_18;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_50;
DENDRO_18 *= DENDRO_4;
DENDRO_18 *= DENDRO_64;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_6;
DENDRO_4 *= DENDRO_54;
DENDRO_4 *= DENDRO_8;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_6;
DENDRO_38 = At3[pp];
DENDRO_24 *= DENDRO_38;
DENDRO_24 *= DENDRO_17;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_6;
DENDRO_58 = At1[pp];
DENDRO_52 *= DENDRO_58;
DENDRO_52 *= DENDRO_46;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_54;
DENDRO_62 *= DENDRO_59;

// initialize reduction 
DENDRO_65 = 0;
DENDRO_65 += DENDRO_29;
DENDRO_65 += DENDRO_21;
DENDRO_65 += DENDRO_16;
DENDRO_65 += DENDRO_9;
DENDRO_65 += DENDRO_7;
DENDRO_65 += DENDRO_18;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_62;
DENDRO_7 += DENDRO_52;
DENDRO_7 += DENDRO_24;
DENDRO_7 += DENDRO_4;
DENDRO_4 = 1.0/3.0;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_4;
DENDRO_9 *= DENDRO_54;
DENDRO_9 *= DENDRO_25;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_4;
DENDRO_16 *= DENDRO_54;
DENDRO_16 *= DENDRO_51;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_37;
DENDRO_18 *= DENDRO_54;
DENDRO_18 *= DENDRO_53;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_26;
DENDRO_21 *= DENDRO_65;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_60;
DENDRO_29 = grad_2_At4;
DENDRO_24 *= DENDRO_29;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_67;
DENDRO_52 = grad_1_At4;
DENDRO_29 *= DENDRO_52;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_69;
DENDRO_62 = grad_0_At4;
DENDRO_52 *= DENDRO_62;

// initialize reduction 
DENDRO_62 = 1;
DENDRO_62 *= DENDRO_20;
DENDRO_62 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_45;
DENDRO_65 = grad_1_beta2;
DENDRO_7 *= DENDRO_65;

// initialize reduction 
DENDRO_66 = 1;
DENDRO_66 *= DENDRO_38;
DENDRO_66 *= DENDRO_61;

// initialize reduction 
DENDRO_68 = 1;
DENDRO_68 *= DENDRO_56;
DENDRO_70 = grad_1_beta0;
DENDRO_68 *= DENDRO_70;

// initialize reduction 
DENDRO_71 = 1;
DENDRO_71 *= DENDRO_58;
DENDRO_71 *= DENDRO_63;

// initialize reduction 
DENDRO_72 = 0;
DENDRO_72 += DENDRO_71;
DENDRO_72 += DENDRO_68;
DENDRO_72 += DENDRO_66;
DENDRO_72 += DENDRO_7;
DENDRO_72 += DENDRO_62;
DENDRO_72 += DENDRO_52;
DENDRO_72 += DENDRO_29;
DENDRO_72 += DENDRO_24;
DENDRO_72 += DENDRO_21;
DENDRO_72 += DENDRO_18;
DENDRO_72 += DENDRO_16;
DENDRO_72 += DENDRO_9;
At_rhs12[pp] = DENDRO_72;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_54;
DENDRO_7 *= DENDRO_22;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_38;
DENDRO_9 *= DENDRO_47;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_58;
DENDRO_16 *= DENDRO_48;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_38;
DENDRO_18 *= DENDRO_43;

// initialize reduction 
DENDRO_21 = 1;
DENDRO_21 *= DENDRO_58;
DENDRO_21 *= DENDRO_49;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_54;
DENDRO_24 *= DENDRO_48;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_38;
DENDRO_29 *= DENDRO_49;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_58;
DENDRO_52 *= DENDRO_44;

// initialize reduction 
DENDRO_62 = 0;
DENDRO_62 += DENDRO_16;
DENDRO_62 += DENDRO_9;
DENDRO_62 += DENDRO_7;

// initialize reduction 
DENDRO_7 = 0;
DENDRO_7 += DENDRO_21;
DENDRO_7 += DENDRO_18;
DENDRO_7 += DENDRO_55;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_52;
DENDRO_9 += DENDRO_29;
DENDRO_9 += DENDRO_24;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_50;
DENDRO_16 *= DENDRO_11;
DENDRO_16 *= DENDRO_64;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_6;
DENDRO_11 *= DENDRO_54;
DENDRO_11 *= DENDRO_62;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_6;
DENDRO_18 *= DENDRO_38;
DENDRO_18 *= DENDRO_7;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_6;
DENDRO_24 *= DENDRO_58;
DENDRO_24 *= DENDRO_9;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_38;
DENDRO_29 *= DENDRO_59;

// initialize reduction 
DENDRO_52 = 0;
DENDRO_52 += DENDRO_34;
DENDRO_52 += DENDRO_41;
DENDRO_52 += DENDRO_15;
DENDRO_52 += DENDRO_5;
DENDRO_52 += DENDRO_39;
DENDRO_52 += DENDRO_16;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_29;
DENDRO_5 += DENDRO_24;
DENDRO_5 += DENDRO_18;
DENDRO_5 += DENDRO_11;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_12;
DENDRO_11 *= DENDRO_38;
DENDRO_11 *= DENDRO_51;

// initialize reduction 
DENDRO_15 = 1;
DENDRO_15 *= DENDRO_37;
DENDRO_15 *= DENDRO_38;
DENDRO_15 *= DENDRO_25;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_37;
DENDRO_16 *= DENDRO_38;
DENDRO_16 *= DENDRO_53;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_40;
DENDRO_18 *= DENDRO_54;
DENDRO_18 *= DENDRO_65;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_40;
DENDRO_24 *= DENDRO_58;
DENDRO_24 *= DENDRO_70;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_26;
DENDRO_29 *= DENDRO_52;

// initialize reduction 
DENDRO_34 = 1;
DENDRO_34 *= DENDRO_60;
DENDRO_39 = grad_2_At3;
DENDRO_34 *= DENDRO_39;

// initialize reduction 
DENDRO_39 = 1;
DENDRO_39 *= DENDRO_67;
DENDRO_41 = grad_1_At3;
DENDRO_39 *= DENDRO_41;

// initialize reduction 
DENDRO_41 = 1;
DENDRO_41 *= DENDRO_69;
DENDRO_52 = grad_0_At3;
DENDRO_41 *= DENDRO_52;

// initialize reduction 
DENDRO_52 = 1;
DENDRO_52 *= DENDRO_20;
DENDRO_52 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_52;
DENDRO_5 += DENDRO_41;
DENDRO_5 += DENDRO_39;
DENDRO_5 += DENDRO_34;
DENDRO_5 += DENDRO_29;
DENDRO_5 += DENDRO_24;
DENDRO_5 += DENDRO_18;
DENDRO_5 += DENDRO_16;
DENDRO_5 += DENDRO_15;
DENDRO_5 += DENDRO_11;
At_rhs11[pp] = DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_50;
DENDRO_5 *= DENDRO_2;
DENDRO_5 *= DENDRO_64;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_6;
DENDRO_2 *= DENDRO_56;
DENDRO_2 *= DENDRO_8;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_6;
DENDRO_8 *= DENDRO_58;
DENDRO_8 *= DENDRO_17;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_6;
DENDRO_15 = At0[pp];
DENDRO_11 *= DENDRO_15;
DENDRO_11 *= DENDRO_46;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_56;
DENDRO_16 *= DENDRO_59;

// initialize reduction 
DENDRO_17 = 0;
DENDRO_17 += DENDRO_23;
DENDRO_17 += DENDRO_32;
DENDRO_17 += DENDRO_31;
DENDRO_17 += DENDRO_27;
DENDRO_17 += DENDRO_30;
DENDRO_17 += DENDRO_5;

// initialize reduction 
DENDRO_5 = 0;
DENDRO_5 += DENDRO_16;
DENDRO_5 += DENDRO_11;
DENDRO_5 += DENDRO_8;
DENDRO_5 += DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_4;
DENDRO_2 *= DENDRO_56;
DENDRO_2 *= DENDRO_25;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_4;
DENDRO_8 *= DENDRO_56;
DENDRO_8 *= DENDRO_53;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_37;
DENDRO_11 *= DENDRO_56;
DENDRO_11 *= DENDRO_51;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_26;
DENDRO_16 *= DENDRO_17;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_60;
DENDRO_18 = grad_2_At2;
DENDRO_17 *= DENDRO_18;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_67;
DENDRO_23 = grad_1_At2;
DENDRO_18 *= DENDRO_23;

// initialize reduction 
DENDRO_23 = 1;
DENDRO_23 *= DENDRO_69;
DENDRO_24 = grad_0_At2;
DENDRO_23 *= DENDRO_24;

// initialize reduction 
DENDRO_24 = 1;
DENDRO_24 *= DENDRO_20;
DENDRO_24 *= DENDRO_5;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_45;
DENDRO_27 = grad_0_beta2;
DENDRO_5 *= DENDRO_27;

// initialize reduction 
DENDRO_29 = 1;
DENDRO_29 *= DENDRO_54;
DENDRO_30 = grad_0_beta1;
DENDRO_29 *= DENDRO_30;

// initialize reduction 
DENDRO_31 = 1;
DENDRO_31 *= DENDRO_58;
DENDRO_31 *= DENDRO_61;

// initialize reduction 
DENDRO_32 = 1;
DENDRO_32 *= DENDRO_15;
DENDRO_32 *= DENDRO_63;

// initialize reduction 
DENDRO_34 = 0;
DENDRO_34 += DENDRO_32;
DENDRO_34 += DENDRO_31;
DENDRO_34 += DENDRO_29;
DENDRO_34 += DENDRO_5;
DENDRO_34 += DENDRO_24;
DENDRO_34 += DENDRO_23;
DENDRO_34 += DENDRO_18;
DENDRO_34 += DENDRO_17;
DENDRO_34 += DENDRO_16;
DENDRO_34 += DENDRO_11;
DENDRO_34 += DENDRO_8;
DENDRO_34 += DENDRO_2;
At_rhs02[pp] = DENDRO_34;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_50;
DENDRO_2 *= DENDRO_0;
DENDRO_2 *= DENDRO_64;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_6;
DENDRO_0 *= DENDRO_56;
DENDRO_0 *= DENDRO_62;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_6;
DENDRO_5 *= DENDRO_58;
DENDRO_5 *= DENDRO_7;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_15;
DENDRO_7 *= DENDRO_9;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_58;
DENDRO_8 *= DENDRO_59;

// initialize reduction 
DENDRO_9 = 0;
DENDRO_9 += DENDRO_28;
DENDRO_9 += DENDRO_36;
DENDRO_9 += DENDRO_19;
DENDRO_9 += DENDRO_33;
DENDRO_9 += DENDRO_35;
DENDRO_9 += DENDRO_2;

// initialize reduction 
DENDRO_2 = 0;
DENDRO_2 += DENDRO_8;
DENDRO_2 += DENDRO_7;
DENDRO_2 += DENDRO_5;
DENDRO_2 += DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_4;
DENDRO_0 *= DENDRO_58;
DENDRO_0 *= DENDRO_51;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_4;
DENDRO_5 *= DENDRO_58;
DENDRO_5 *= DENDRO_53;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_37;
DENDRO_4 *= DENDRO_58;
DENDRO_4 *= DENDRO_25;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_26;
DENDRO_7 *= DENDRO_9;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_60;
DENDRO_9 = grad_2_At1;
DENDRO_8 *= DENDRO_9;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_67;
DENDRO_11 = grad_1_At1;
DENDRO_9 *= DENDRO_11;

// initialize reduction 
DENDRO_11 = 1;
DENDRO_11 *= DENDRO_69;
DENDRO_16 = grad_0_At1;
DENDRO_11 *= DENDRO_16;

// initialize reduction 
DENDRO_16 = 1;
DENDRO_16 *= DENDRO_20;
DENDRO_16 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_54;
DENDRO_2 *= DENDRO_27;

// initialize reduction 
DENDRO_17 = 1;
DENDRO_17 *= DENDRO_38;
DENDRO_17 *= DENDRO_30;

// initialize reduction 
DENDRO_18 = 1;
DENDRO_18 *= DENDRO_56;
DENDRO_18 *= DENDRO_65;

// initialize reduction 
DENDRO_19 = 1;
DENDRO_19 *= DENDRO_15;
DENDRO_19 *= DENDRO_70;

// initialize reduction 
DENDRO_23 = 0;
DENDRO_23 += DENDRO_19;
DENDRO_23 += DENDRO_18;
DENDRO_23 += DENDRO_17;
DENDRO_23 += DENDRO_2;
DENDRO_23 += DENDRO_16;
DENDRO_23 += DENDRO_11;
DENDRO_23 += DENDRO_9;
DENDRO_23 += DENDRO_8;
DENDRO_23 += DENDRO_7;
DENDRO_23 += DENDRO_4;
DENDRO_23 += DENDRO_5;
DENDRO_23 += DENDRO_0;
At_rhs01[pp] = DENDRO_23;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_56;
DENDRO_0 *= DENDRO_22;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_58;
DENDRO_2 *= DENDRO_47;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_15;
DENDRO_4 *= DENDRO_48;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_56;
DENDRO_5 *= DENDRO_47;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_58;
DENDRO_7 *= DENDRO_43;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_15;
DENDRO_8 *= DENDRO_49;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_15;
DENDRO_9 *= DENDRO_44;

// initialize reduction 
DENDRO_11 = 0;
DENDRO_11 += DENDRO_4;
DENDRO_11 += DENDRO_2;
DENDRO_11 += DENDRO_0;

// initialize reduction 
DENDRO_0 = 0;
DENDRO_0 += DENDRO_8;
DENDRO_0 += DENDRO_7;
DENDRO_0 += DENDRO_5;

// initialize reduction 
DENDRO_2 = 0;
DENDRO_2 += DENDRO_9;
DENDRO_2 += DENDRO_21;
DENDRO_2 += DENDRO_57;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_50;
DENDRO_4 *= DENDRO_10;
DENDRO_4 *= DENDRO_64;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_6;
DENDRO_5 *= DENDRO_56;
DENDRO_5 *= DENDRO_11;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_6;
DENDRO_7 *= DENDRO_58;
DENDRO_7 *= DENDRO_0;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_6;
DENDRO_0 *= DENDRO_15;
DENDRO_0 *= DENDRO_2;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_15;
DENDRO_2 *= DENDRO_59;

// initialize reduction 
DENDRO_6 = 0;
DENDRO_6 += DENDRO_13;
DENDRO_6 += DENDRO_14;
DENDRO_6 += DENDRO_1;
DENDRO_6 += DENDRO_3;
DENDRO_6 += DENDRO_42;
DENDRO_6 += DENDRO_4;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_2;
DENDRO_1 += DENDRO_0;
DENDRO_1 += DENDRO_7;
DENDRO_1 += DENDRO_5;

// initialize reduction 
DENDRO_0 = 1;
DENDRO_0 *= DENDRO_12;
DENDRO_0 *= DENDRO_15;
DENDRO_0 *= DENDRO_53;

// initialize reduction 
DENDRO_2 = 1;
DENDRO_2 *= DENDRO_37;
DENDRO_2 *= DENDRO_15;
DENDRO_2 *= DENDRO_25;

// initialize reduction 
DENDRO_3 = 1;
DENDRO_3 *= DENDRO_37;
DENDRO_3 *= DENDRO_15;
DENDRO_3 *= DENDRO_51;

// initialize reduction 
DENDRO_4 = 1;
DENDRO_4 *= DENDRO_40;
DENDRO_4 *= DENDRO_56;
DENDRO_4 *= DENDRO_27;

// initialize reduction 
DENDRO_5 = 1;
DENDRO_5 *= DENDRO_40;
DENDRO_5 *= DENDRO_58;
DENDRO_5 *= DENDRO_30;

// initialize reduction 
DENDRO_7 = 1;
DENDRO_7 *= DENDRO_26;
DENDRO_7 *= DENDRO_6;

// initialize reduction 
DENDRO_6 = 1;
DENDRO_6 *= DENDRO_60;
DENDRO_8 = grad_2_At0;
DENDRO_6 *= DENDRO_8;

// initialize reduction 
DENDRO_8 = 1;
DENDRO_8 *= DENDRO_67;
DENDRO_9 = grad_1_At0;
DENDRO_8 *= DENDRO_9;

// initialize reduction 
DENDRO_9 = 1;
DENDRO_9 *= DENDRO_69;
DENDRO_10 = grad_0_At0;
DENDRO_9 *= DENDRO_10;

// initialize reduction 
DENDRO_10 = 1;
DENDRO_10 *= DENDRO_20;
DENDRO_10 *= DENDRO_1;

// initialize reduction 
DENDRO_1 = 0;
DENDRO_1 += DENDRO_10;
DENDRO_1 += DENDRO_9;
DENDRO_1 += DENDRO_8;
DENDRO_1 += DENDRO_6;
DENDRO_1 += DENDRO_7;
DENDRO_1 += DENDRO_5;
DENDRO_1 += DENDRO_4;
DENDRO_1 += DENDRO_3;
DENDRO_1 += DENDRO_2;
DENDRO_1 += DENDRO_0;
At_rhs00[pp] = DENDRO_1;


}
if(blk->m_bflag != 0){
radiative_bc_pt<pw , nx>(&At_rhs00[gidx], At0[gidx], grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&At_rhs01[gidx], At1[gidx], grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&At_rhs02[gidx], At2[gidx], grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&At_rhs11[gidx], At3[gidx], grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&At_rhs12[gidx], At4[gidx], grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk);
radiative_bc_pt<pw , nx>(&At_rhs22[gidx], At5[gidx], grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk);
}

At_rhs00[pp] += ko_sigma * (kograd_0_At0 + kograd_1_At0 + kograd_2_At0);
At_rhs01[pp] += ko_sigma * (kograd_0_At1 + kograd_1_At1 + kograd_2_At1);
At_rhs02[pp] += ko_sigma * (kograd_0_At2 + kograd_1_At2 + kograd_2_At2);
At_rhs11[pp] += ko_sigma * (kograd_0_At3 + kograd_1_At3 + kograd_2_At3);
At_rhs12[pp] += ko_sigma * (kograd_0_At4 + kograd_1_At4 + kograd_2_At4);
At_rhs22[pp] += ko_sigma * (kograd_0_At5 + kograd_1_At5 + kograd_2_At5);
__syncthreads();

