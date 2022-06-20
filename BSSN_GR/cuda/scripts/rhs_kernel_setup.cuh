// set up kernel parameters for RHS GPU RHS kernels. 
const DEVICE_UINT BLK_ID            = GPUDevice::block_id_z();
const DEVICE_UINT BLK_SZ            = (blk[BLK_ID].m_aligned_sz[0]) * (blk[BLK_ID].m_aligned_sz[1]) * (blk[BLK_ID].m_aligned_sz[2]); 
const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
const DEVICE_UINT sz_per_dof        = szpdof_uz;

DEVICE_REAL * const a_rhs    = &Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
DEVICE_REAL * const chi_rhs  = &Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
DEVICE_REAL * const gt_rhs00 = &Fu[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
DEVICE_REAL * const gt_rhs01 = &Fu[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
DEVICE_REAL * const gt_rhs02 = &Fu[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
DEVICE_REAL * const gt_rhs11 = &Fu[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
DEVICE_REAL * const gt_rhs12 = &Fu[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
DEVICE_REAL * const gt_rhs22 = &Fu[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
DEVICE_REAL * const b_rhs0   = &Fu[bssn::VAR::U_BETA0 * sz_per_dof + offset];
DEVICE_REAL * const b_rhs1   = &Fu[bssn::VAR::U_BETA1 * sz_per_dof + offset];
DEVICE_REAL * const b_rhs2   = &Fu[bssn::VAR::U_BETA2 * sz_per_dof + offset];
DEVICE_REAL * const At_rhs00 = &Fu[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
DEVICE_REAL * const At_rhs01 = &Fu[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
DEVICE_REAL * const At_rhs02 = &Fu[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
DEVICE_REAL * const At_rhs11 = &Fu[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
DEVICE_REAL * const At_rhs12 = &Fu[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
DEVICE_REAL * const At_rhs22 = &Fu[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
DEVICE_REAL * const Gt_rhs0  = &Fu[bssn::VAR::U_GT0 * sz_per_dof + offset];
DEVICE_REAL * const Gt_rhs1  = &Fu[bssn::VAR::U_GT1 * sz_per_dof + offset];
DEVICE_REAL * const Gt_rhs2  = &Fu[bssn::VAR::U_GT2 * sz_per_dof + offset];
DEVICE_REAL * const B_rhs0   = &Fu[bssn::VAR::U_B0 * sz_per_dof + offset];
DEVICE_REAL * const B_rhs1   = &Fu[bssn::VAR::U_B1 * sz_per_dof + offset];
DEVICE_REAL * const B_rhs2   = &Fu[bssn::VAR::U_B2 * sz_per_dof + offset];
DEVICE_REAL * const K_rhs    = &Fu[bssn::VAR::U_K * sz_per_dof + offset];

const DEVICE_UINT lambda[4]     = {1, 1, 1, 1};
const DEVICE_REAL lambda_f[2]   = {1.0, 0.0};


const DEVICE_INT nx = blk[BLK_ID].m_sz[0];
const DEVICE_INT ny = blk[BLK_ID].m_sz[1];
const DEVICE_INT nz = blk[BLK_ID].m_sz[2];

const DEVICE_REAL hx = blk[BLK_ID].m_dx[0];
const DEVICE_REAL hy = blk[BLK_ID].m_dx[1];
const DEVICE_REAL hz = blk[BLK_ID].m_dx[2];

// const DEVICE_INT i  = GPUDevice::thread_id_x() + pw;
// const DEVICE_INT j  = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y() + pw;
// const DEVICE_INT k  = GPUDevice::block_id_y() + pw;

const DEVICE_INT i  = GPUDevice::thread_id_x();
const DEVICE_INT j  = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
const DEVICE_INT k  = GPUDevice::block_id_y();

// if(i >=nx-pw || j >=nx-pw || k >=nx-pw)
//     return;

const DEVICE_INT si = GPUDevice::thread_id_x() + pw;  
const DEVICE_INT sj = GPUDevice::thread_id_y();

const DEVICE_INT sb = si;
const DEVICE_INT pp = k * nx * ny + j * nx + sb;

// const DEVICE_INT i  = GPUDevice::thread_id_x();
// const DEVICE_INT j  = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
// const DEVICE_INT k  = GPUDevice::block_id_y();

// const DEVICE_INT si = GPUDevice::thread_id_x();  
// const DEVICE_INT sj = GPUDevice::thread_id_y();

// const DEVICE_INT sb = si;
// const DEVICE_INT pp = k * nx * ny + j * nx + sb;

const DEVICE_REAL x = blk[BLK_ID].m_ptMin[0] + i*hx;
const DEVICE_REAL y = blk[BLK_ID].m_ptMin[1] + j*hy;
const DEVICE_REAL z = blk[BLK_ID].m_ptMin[2] + k*hz;