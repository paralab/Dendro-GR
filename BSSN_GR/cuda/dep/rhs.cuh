/**
     * @brief Launch the kernel to compute all the z-direction derivatives
    */
    template<int pw, int pencils, int pencil_sz, int BATCHED_BLOCK_SZ>
    GLOBAL_FUNC void launch_dir_z_deriv_kernel(const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz) {
        
        //includes z-derivs here. 
        const DEVICE_UINT BLK_ID            = GPUDevice::block_id_z();
        const DEVICE_UINT BLK_SZ            = (blk[BLK_ID].m_aligned_sz[0]) * (blk[BLK_ID].m_aligned_sz[1]) * (blk[BLK_ID].m_aligned_sz[2]); 
        const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
        const DEVICE_UINT sz_per_dof        = szpdof_uz;

        const DEVICE_REAL * const alpha = &u[bssn::VAR::U_ALPHA * sz_per_dof + offset];
        const DEVICE_REAL * const chi   = &u[bssn::VAR::U_CHI * sz_per_dof + offset];
        const DEVICE_REAL * const K     = &u[bssn::VAR::U_K * sz_per_dof + offset];
        const DEVICE_REAL * const gt0   = &u[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
        const DEVICE_REAL * const gt1   = &u[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
        const DEVICE_REAL * const gt2   = &u[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
        const DEVICE_REAL * const gt3   = &u[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
        const DEVICE_REAL * const gt4   = &u[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
        const DEVICE_REAL * const gt5   = &u[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
        const DEVICE_REAL * const beta0 = &u[bssn::VAR::U_BETA0 * sz_per_dof + offset];
        const DEVICE_REAL * const beta1 = &u[bssn::VAR::U_BETA1 * sz_per_dof + offset];
        const DEVICE_REAL * const beta2 = &u[bssn::VAR::U_BETA2 * sz_per_dof + offset];
        const DEVICE_REAL * const At0   = &u[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
        const DEVICE_REAL * const At1   = &u[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
        const DEVICE_REAL * const At2   = &u[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
        const DEVICE_REAL * const At3   = &u[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
        const DEVICE_REAL * const At4   = &u[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
        const DEVICE_REAL * const At5   = &u[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
        const DEVICE_REAL * const Gt0   = &u[bssn::VAR::U_GT0 * sz_per_dof + offset];
        const DEVICE_REAL * const Gt1   = &u[bssn::VAR::U_GT1 * sz_per_dof + offset];
        const DEVICE_REAL * const Gt2   = &u[bssn::VAR::U_GT2 * sz_per_dof + offset];
        const DEVICE_REAL * const B0    = &u[bssn::VAR::U_B0 * sz_per_dof + offset];
        const DEVICE_REAL * const B1    = &u[bssn::VAR::U_B1 * sz_per_dof + offset];
        const DEVICE_REAL * const B2    = &u[bssn::VAR::U_B2 * sz_per_dof + offset];

        #include "../scripts/bssnrhs_deriv_z_dir.cuh"

        return;


    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_light(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {

        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL * const alpha = (u  + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL * const chi   = (u  + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL * const K     = (u  + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL * const gt0   = (u  + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL * const gt1   = (u  + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL * const gt2   = (u  + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL * const gt3   = (u  + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL * const gt4   = (u  + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL * const gt5   = (u  + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL * const beta0 = (u  + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL * const beta1 = (u  + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL * const beta2 = (u  + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL * const At0   = (u  + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL * const At1   = (u  + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL * const At2   = (u  + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL * const At3   = (u  + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL * const At4   = (u  + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL * const At5   = (u  + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL * const Gt0   = (u  + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL * const Gt1   = (u  + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL * const Gt2   = (u  + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL * const B0    = (u  + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL * const B1    = (u  + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL * const B2    = (u  + bssn::VAR::U_B2 * sz_per_dof + offset);
        
        
        #include "../scripts/bssnrhs_evar_derivs.h"

        #include "../scripts/a_beta_gt_chi.cpp"
        
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_R00(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL * const alpha = (u  + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL * const chi   = (u  + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL * const K     = (u  + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL * const gt0   = (u  + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL * const gt1   = (u  + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL * const gt2   = (u  + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL * const gt3   = (u  + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL * const gt4   = (u  + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL * const gt5   = (u  + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL * const beta0 = (u  + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL * const beta1 = (u  + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL * const beta2 = (u  + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL * const At0   = (u  + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL * const At1   = (u  + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL * const At2   = (u  + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL * const At3   = (u  + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL * const At4   = (u  + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL * const At5   = (u  + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL * const Gt0   = (u  + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL * const Gt1   = (u  + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL * const Gt2   = (u  + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL * const B0    = (u  + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL * const B1    = (u  + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL * const B2    = (u  + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/R_00.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_R01(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/R_01.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_R02(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/R_02.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_R11(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/R_11.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_R12(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/R_12.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_R22(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/R_22.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_At_rhs_tf(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/At_rhs_tf.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_At_rhs(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/At_rhs.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_K_rhs1(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/K_rhs1.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_K_rhs2(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        #include "../scripts/K_rhs2.cpp"
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_Gt_B_rhs(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"
        
        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;
        const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

        const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
        const DEVICE_REAL arg = - w*w*w*w;
        const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
        
        #include "../scripts/Gt_B_rhs.cpp"
        
        return;

    }
    
    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bssnrhs_heavy(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
        
        #include "../scripts/rhs_kernel_setup.cuh"
        
        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;
        const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

        const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
        const DEVICE_REAL arg = - w*w*w*w;
        const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
        
        const DEVICE_REAL alpha = *(u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL chi   = *(u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL K     = *(u + pp + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL gt0   = *(u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL gt1   = *(u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL gt2   = *(u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL gt3   = *(u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL gt4   = *(u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL gt5   = *(u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL beta0 = *(u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL beta1 = *(u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL beta2 = *(u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL At0   = *(u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL At1   = *(u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL At2   = *(u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL At3   = *(u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL At4   = *(u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL At5   = *(u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL Gt0   = *(u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL Gt1   = *(u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL Gt2   = *(u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL B0    = *(u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL B1    = *(u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL B2    = *(u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
        #include "../scripts/bssnrhs_evar_derivs.h"

        #include "../scripts/bssneqs_eta_const_standard_gauge4.cpp"
        
        
        return;

    }


    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_bc(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz)
    {
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
        
        if(blk[BLK_ID].m_bflag != 0)
        {
            const DEVICE_REAL* const alpha = (u + bssn::VAR::U_ALPHA * sz_per_dof + offset);
            const DEVICE_REAL* const chi   = (u + bssn::VAR::U_CHI * sz_per_dof + offset);
            const DEVICE_REAL* const K     = (u + bssn::VAR::U_K * sz_per_dof + offset);
            const DEVICE_REAL* const gt0   = (u + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
            const DEVICE_REAL* const gt1   = (u + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
            const DEVICE_REAL* const gt2   = (u + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
            const DEVICE_REAL* const gt3   = (u + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
            const DEVICE_REAL* const gt4   = (u + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
            const DEVICE_REAL* const gt5   = (u + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
            const DEVICE_REAL* const beta0 = (u + bssn::VAR::U_BETA0 * sz_per_dof + offset);
            const DEVICE_REAL* const beta1 = (u + bssn::VAR::U_BETA1 * sz_per_dof + offset);
            const DEVICE_REAL* const beta2 = (u + bssn::VAR::U_BETA2 * sz_per_dof + offset);
            const DEVICE_REAL* const At0   = (u + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
            const DEVICE_REAL* const At1   = (u + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
            const DEVICE_REAL* const At2   = (u + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
            const DEVICE_REAL* const At3   = (u + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
            const DEVICE_REAL* const At4   = (u + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
            const DEVICE_REAL* const At5   = (u + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
            const DEVICE_REAL* const Gt0   = (u + bssn::VAR::U_GT0 * sz_per_dof + offset);
            const DEVICE_REAL* const Gt1   = (u + bssn::VAR::U_GT1 * sz_per_dof + offset);
            const DEVICE_REAL* const Gt2   = (u + bssn::VAR::U_GT2 * sz_per_dof + offset);
            const DEVICE_REAL* const B0    = (u + bssn::VAR::U_B0 * sz_per_dof + offset);
            const DEVICE_REAL* const B1    = (u + bssn::VAR::U_B1 * sz_per_dof + offset);
            const DEVICE_REAL* const B2    = (u + bssn::VAR::U_B2 * sz_per_dof + offset);

            const DEVICE_REAL* const grad_0_alpha = (deriv_evars->grad_0_alpha + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_alpha = (deriv_evars->grad_1_alpha + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_alpha = (deriv_evars->grad_2_alpha + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_beta0 = (deriv_evars->grad_0_beta0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_beta0 = (deriv_evars->grad_1_beta0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_beta0 = (deriv_evars->grad_2_beta0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_beta1 = (deriv_evars->grad_0_beta1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_beta1 = (deriv_evars->grad_1_beta1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_beta1 = (deriv_evars->grad_2_beta1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_beta2 = (deriv_evars->grad_0_beta2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_beta2 = (deriv_evars->grad_1_beta2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_beta2 = (deriv_evars->grad_2_beta2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_B0    = (deriv_evars->grad_0_B0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_B0    = (deriv_evars->grad_1_B0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_B0    = (deriv_evars->grad_2_B0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_B1    = (deriv_evars->grad_0_B1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_B1    = (deriv_evars->grad_1_B1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_B1    = (deriv_evars->grad_2_B1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_B2    = (deriv_evars->grad_0_B2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_B2    = (deriv_evars->grad_1_B2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_B2    = (deriv_evars->grad_2_B2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_chi   = (deriv_evars->grad_0_chi + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_chi   = (deriv_evars->grad_1_chi + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_chi   = (deriv_evars->grad_2_chi + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_Gt0   = (deriv_evars->grad_0_Gt0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_Gt0   = (deriv_evars->grad_1_Gt0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_Gt0   = (deriv_evars->grad_2_Gt0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_Gt1   = (deriv_evars->grad_0_Gt1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_Gt1   = (deriv_evars->grad_1_Gt1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_Gt1   = (deriv_evars->grad_2_Gt1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_Gt2   = (deriv_evars->grad_0_Gt2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_Gt2   = (deriv_evars->grad_1_Gt2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_Gt2   = (deriv_evars->grad_2_Gt2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_K     = (deriv_evars->grad_0_K + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_K     = (deriv_evars->grad_1_K + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_K     = (deriv_evars->grad_2_K + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_gt0   = (deriv_evars->grad_0_gt0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_gt0   = (deriv_evars->grad_1_gt0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_gt0   = (deriv_evars->grad_2_gt0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_gt1   = (deriv_evars->grad_0_gt1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_gt1   = (deriv_evars->grad_1_gt1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_gt1   = (deriv_evars->grad_2_gt1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_gt2   = (deriv_evars->grad_0_gt2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_gt2   = (deriv_evars->grad_1_gt2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_gt2   = (deriv_evars->grad_2_gt2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_gt3   = (deriv_evars->grad_0_gt3 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_gt3   = (deriv_evars->grad_1_gt3 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_gt3   = (deriv_evars->grad_2_gt3 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_gt4   = (deriv_evars->grad_0_gt4 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_gt4   = (deriv_evars->grad_1_gt4 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_gt4   = (deriv_evars->grad_2_gt4 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_gt5   = (deriv_evars->grad_0_gt5 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_gt5   = (deriv_evars->grad_1_gt5 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_gt5   = (deriv_evars->grad_2_gt5 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_At0   = (deriv_evars->grad_0_At0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_At0   = (deriv_evars->grad_1_At0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_At0   = (deriv_evars->grad_2_At0 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_At1   = (deriv_evars->grad_0_At1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_At1   = (deriv_evars->grad_1_At1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_At1   = (deriv_evars->grad_2_At1 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_At2   = (deriv_evars->grad_0_At2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_At2   = (deriv_evars->grad_1_At2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_At2   = (deriv_evars->grad_2_At2 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_At3   = (deriv_evars->grad_0_At3 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_At3   = (deriv_evars->grad_1_At3 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_At3   = (deriv_evars->grad_2_At3 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_At4   = (deriv_evars->grad_0_At4 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_At4   = (deriv_evars->grad_1_At4 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_At4   = (deriv_evars->grad_2_At4 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_0_At5   = (deriv_evars->grad_0_At5 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_1_At5   = (deriv_evars->grad_1_At5 + BLK_ID * BLK_SZ );
            const DEVICE_REAL* const grad_2_At5   = (deriv_evars->grad_2_At5 + BLK_ID * BLK_SZ );
            
            radiative_bc<pw>(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk,BLK_ID);
            radiative_bc<pw>(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk, BLK_ID);
            radiative_bc<pw>(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, 1.0, 0.0, blk, BLK_ID);

            radiative_bc<pw>(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk, BLK_ID);

            radiative_bc<pw>(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk, BLK_ID);

            radiative_bc<pw>(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk, BLK_ID);

            radiative_bc<pw>(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk, BLK_ID);

            radiative_bc<pw>(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk, BLK_ID);
            radiative_bc<pw>(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk, BLK_ID);
            radiative_bc<pw>(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk, BLK_ID);
            radiative_bc<pw>(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk, BLK_ID);

        }

        return;
        
        
    }


    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void launch_ko(DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_REAL ko_sigma, DEVICE_UINT szpdof_uz)
    {
        #include "../scripts/rhs_kernel_setup.cuh"
        #include "../scripts/bssnrhs_evar_derivs.h"

        const DEVICE_REAL sigma = ko_sigma;

        a_rhs[pp]    += sigma * (kograd_0_alpha + kograd_1_alpha + kograd_2_alpha);
        
        b_rhs0[pp]   += sigma * (kograd_0_beta0 + kograd_1_beta0 + kograd_2_beta0);
        b_rhs1[pp]   += sigma * (kograd_0_beta1 + kograd_1_beta1 + kograd_2_beta1);
        b_rhs2[pp]   += sigma * (kograd_0_beta2 + kograd_1_beta2 + kograd_2_beta2);

        gt_rhs00[pp] += sigma * (kograd_0_gt0 + kograd_1_gt0 + kograd_2_gt0);
        gt_rhs01[pp] += sigma * (kograd_0_gt1 + kograd_1_gt1 + kograd_2_gt1);
        gt_rhs02[pp] += sigma * (kograd_0_gt2 + kograd_1_gt2 + kograd_2_gt2);
        gt_rhs11[pp] += sigma * (kograd_0_gt3 + kograd_1_gt3 + kograd_2_gt3);
        gt_rhs12[pp] += sigma * (kograd_0_gt4 + kograd_1_gt4 + kograd_2_gt4);
        gt_rhs22[pp] += sigma * (kograd_0_gt5 + kograd_1_gt5 + kograd_2_gt5);

        chi_rhs[pp]  += sigma * (kograd_0_chi + kograd_1_chi + kograd_2_chi);

        At_rhs00[pp] += sigma * (kograd_0_At0 + kograd_1_At0 + kograd_2_At0);
        At_rhs01[pp] += sigma * (kograd_0_At1 + kograd_1_At1 + kograd_2_At1);
        At_rhs02[pp] += sigma * (kograd_0_At2 + kograd_1_At2 + kograd_2_At2);
        At_rhs11[pp] += sigma * (kograd_0_At3 + kograd_1_At3 + kograd_2_At3);
        At_rhs12[pp] += sigma * (kograd_0_At4 + kograd_1_At4 + kograd_2_At4);
        At_rhs22[pp] += sigma * (kograd_0_At5 + kograd_1_At5 + kograd_2_At5);

        K_rhs[pp]    += sigma * (kograd_0_K + kograd_1_K + kograd_2_K);

        Gt_rhs0[pp]  += sigma * (kograd_0_Gt0 + kograd_1_Gt0 + kograd_2_Gt0);
        Gt_rhs1[pp]  += sigma * (kograd_0_Gt1 + kograd_1_Gt1 + kograd_2_Gt1);
        Gt_rhs2[pp]  += sigma * (kograd_0_Gt2 + kograd_1_Gt2 + kograd_2_Gt2);

        B_rhs0[pp]   += sigma * (kograd_0_B0 + kograd_1_B0 + kograd_2_B0);
        B_rhs1[pp]   += sigma * (kograd_0_B1 + kograd_1_B1 + kograd_2_B1);
        B_rhs2[pp]   += sigma * (kograd_0_B2 + kograd_1_B2 + kograd_2_B2);
        
        return;

    }

    template<int pw, int pencils, int pencil_sz, unsigned int BATCHED_GRAIN_SZ, unsigned int nx>
    HOST_FUNC void bssnrhs_gpu(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const dptr_deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {
        const unsigned int ny=nx;
        const unsigned int nz=nx;
        const unsigned int BLK_SZ = nx*ny*nz;

        const unsigned int NUM_BATCHES  = nblocks/BATCHED_GRAIN_SZ + 1;
        cudaFuncSetCacheConfig(&device::launch_bssnrhs_heavy<pw,pencils,pencil_sz>, cudaFuncCachePreferL1);
        cudaStream_t s_default= (cudaStream_t)0;
        
        for(unsigned int bid =0; bid<NUM_BATCHES; bid++)
        {
            const unsigned int block_begin = (bid * nblocks)/NUM_BATCHES;
            const unsigned int block_end   = ((bid+1) * nblocks)/NUM_BATCHES;
            const unsigned int numblocks   = block_end-block_begin;  
            //const unsigned int BATCHED_BLOCKS_SZ = numblocks;

            dim3 grid_x  = dim3(ny/pencils,nz,numblocks);
            dim3 block_x = dim3(nx,pencils,1);

            dim3 grid_y  = dim3(nx/pencils,nz,numblocks);
            dim3 block_y = dim3(pencils,ny,1);

            dim3 grid_z  = dim3(nx/pencils,ny,numblocks);
            dim3 block_z = dim3(pencils,nz,1);

            device::launch_dir_x_deriv_kernel<pw,pencils,pencil_sz,BATCHED_GRAIN_SZ> <<<grid_x,block_x,0>>>(dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks,szpdof_uz);
            device::launch_dir_y_deriv_kernel<pw,pencils,pencil_sz,BATCHED_GRAIN_SZ> <<<grid_y,block_y,0>>>(dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks,szpdof_uz);
            device::launch_dir_z_deriv_kernel<pw,pencils,pencil_sz,BATCHED_GRAIN_SZ> <<<grid_z,block_z,0>>>(dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks,szpdof_uz);
            
            // unstaged rhs
            //device::launch_bssnrhs_heavy<pw,pencils,pencil_sz>      <<<grid_x,block_x,0>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);

            // manually split rhs. 
            device::launch_bssnrhs_light<pw,pencils,pencil_sz>    <<<grid_x,block_x,0,s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_K_rhs1<pw,pencils,pencil_sz>   <<<grid_x,block_x,0,s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_K_rhs2<pw,pencils,pencil_sz>   <<<grid_x,block_x,0,s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_Gt_B_rhs<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_R00<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default >>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_R01<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default >>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_R02<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default >>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_R11<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default >>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_R12<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default >>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_R22<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default >>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_bssnrhs_At_rhs_tf<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks,szpdof_uz);
            device::launch_bssnrhs_At_rhs<pw,pencils,pencil_sz> <<<grid_x,block_x,0,s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);

            device::launch_bc<pw,pencils,pencil_sz> <<<grid_x,block_x,0, s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks, szpdof_uz);
            device::launch_ko<pw,pencils,pencil_sz> <<<grid_x,block_x,0, s_default>>>(dptr_Fu, dptr_u, dptr_deriv_evars, dptr_blk + block_begin,numblocks,ko_sigma, szpdof_uz);
            
        }
        
        return;

    }

// attempt to do coperative group launch : has a sync bug, does not produce the correct result : 02-20-2022
    template<int pw, int pencils, int pencil_sz, unsigned int BATCHED_GRAIN_SZ, unsigned int nx>
    GLOBAL_FUNC 
    __launch_bounds__(256)
    void bssnrhs1(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const dptr_deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {

        const DEVICE_UINT NUM_BATCHES  = nblocks/BATCHED_GRAIN_SZ + 1;
        const DEVICE_UINT blk_stride   = GPUDevice::grid_dim_z(); 
        const DEVICE_UINT Z_ID         = GPUDevice::block_id_z();
        const DEVICE_UINT block_begin  = Z_ID;
        const DEVICE_UINT sz_per_dof   = szpdof_uz;

        const DEVICE_UINT lambda[4]    = {1, 1, 1, 1};
        const DEVICE_REAL lambda_f[2]  = {1.0, 0.0};

        //SHARED_MEM DEVICE_REAL su[nx*nx*nx];
        const BlockGPU3D * const blk         =  dptr_blk;
        BSSN_EVAR_DERIVS * const deriv_evars = dptr_deriv_evars;

        const DEVICE_UINT si = GPUDevice::thread_id_x() + pw;
        const DEVICE_UINT sj = GPUDevice::thread_id_y() + pw;
        cooperative_groups::grid_group g = cooperative_groups::this_grid(); 
        
        for(unsigned int bid = block_begin ; bid < nblocks; bid+=blk_stride )
        {

            // if(si ==pw && sj ==pw)
            //     printf("blk %d/%d with z_id %d\n",bid,nblocks,Z_ID);

            const DEVICE_UINT BLK_ID            = bid;
            const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
            const DEVICE_UINT actual_nx         = blk[BLK_ID].m_sz[0];
            //const DEVICE_UINT nx                = blk[BLK_ID].m_aligned_sz[0];
            const DEVICE_UINT BLK_SZ            = nx*nx*nx; 

            const DEVICE_REAL hx = blk[BLK_ID].m_dx[0];
            const DEVICE_REAL hy = blk[BLK_ID].m_dx[1];
            const DEVICE_REAL hz = blk[BLK_ID].m_dx[2];

            {
                const DEVICE_REAL * const alpha = &dptr_u[bssn::VAR::U_ALPHA * sz_per_dof + offset];
                const DEVICE_REAL * const chi   = &dptr_u[bssn::VAR::U_CHI * sz_per_dof + offset];
                const DEVICE_REAL * const K     = &dptr_u[bssn::VAR::U_K * sz_per_dof + offset];
                const DEVICE_REAL * const gt0   = &dptr_u[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
                const DEVICE_REAL * const gt1   = &dptr_u[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
                const DEVICE_REAL * const gt2   = &dptr_u[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
                const DEVICE_REAL * const gt3   = &dptr_u[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
                const DEVICE_REAL * const gt4   = &dptr_u[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
                const DEVICE_REAL * const gt5   = &dptr_u[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
                const DEVICE_REAL * const beta0 = &dptr_u[bssn::VAR::U_BETA0 * sz_per_dof + offset];
                const DEVICE_REAL * const beta1 = &dptr_u[bssn::VAR::U_BETA1 * sz_per_dof + offset];
                const DEVICE_REAL * const beta2 = &dptr_u[bssn::VAR::U_BETA2 * sz_per_dof + offset];
                const DEVICE_REAL * const At0   = &dptr_u[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
                const DEVICE_REAL * const At1   = &dptr_u[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
                const DEVICE_REAL * const At2   = &dptr_u[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
                const DEVICE_REAL * const At3   = &dptr_u[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
                const DEVICE_REAL * const At4   = &dptr_u[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
                const DEVICE_REAL * const At5   = &dptr_u[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
                const DEVICE_REAL * const Gt0   = &dptr_u[bssn::VAR::U_GT0 * sz_per_dof + offset];
                const DEVICE_REAL * const Gt1   = &dptr_u[bssn::VAR::U_GT1 * sz_per_dof + offset];
                const DEVICE_REAL * const Gt2   = &dptr_u[bssn::VAR::U_GT2 * sz_per_dof + offset];
                const DEVICE_REAL * const B0    = &dptr_u[bssn::VAR::U_B0 * sz_per_dof + offset];
                const DEVICE_REAL * const B1    = &dptr_u[bssn::VAR::U_B1 * sz_per_dof + offset];
                const DEVICE_REAL * const B2    = &dptr_u[bssn::VAR::U_B2 * sz_per_dof + offset];
                
                #include "../scripts/compute_bssnrhs_evar_derivs1.cuh"

            }
            
            g.sync();
            
            {
                DEVICE_REAL * const Fu = dptr_Fu;
                #include "../scripts/rhs_kernel_setup.cuh"
                const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
                const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
                const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;
                const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

                const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
                const DEVICE_REAL arg = - w*w*w*w;
                const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
                
                const DEVICE_REAL alpha = *(dptr_u + pp + bssn::VAR::U_ALPHA * sz_per_dof + offset);
                const DEVICE_REAL chi   = *(dptr_u + pp + bssn::VAR::U_CHI * sz_per_dof + offset);
                const DEVICE_REAL K     = *(dptr_u + pp + bssn::VAR::U_K * sz_per_dof + offset);
                const DEVICE_REAL gt0   = *(dptr_u + pp + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
                const DEVICE_REAL gt1   = *(dptr_u + pp + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
                const DEVICE_REAL gt2   = *(dptr_u + pp + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
                const DEVICE_REAL gt3   = *(dptr_u + pp + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
                const DEVICE_REAL gt4   = *(dptr_u + pp + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
                const DEVICE_REAL gt5   = *(dptr_u + pp + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
                const DEVICE_REAL beta0 = *(dptr_u + pp + bssn::VAR::U_BETA0 * sz_per_dof + offset);
                const DEVICE_REAL beta1 = *(dptr_u + pp + bssn::VAR::U_BETA1 * sz_per_dof + offset);
                const DEVICE_REAL beta2 = *(dptr_u + pp + bssn::VAR::U_BETA2 * sz_per_dof + offset);
                const DEVICE_REAL At0   = *(dptr_u + pp + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
                const DEVICE_REAL At1   = *(dptr_u + pp + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
                const DEVICE_REAL At2   = *(dptr_u + pp + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
                const DEVICE_REAL At3   = *(dptr_u + pp + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
                const DEVICE_REAL At4   = *(dptr_u + pp + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
                const DEVICE_REAL At5   = *(dptr_u + pp + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
                const DEVICE_REAL Gt0   = *(dptr_u + pp + bssn::VAR::U_GT0 * sz_per_dof + offset);
                const DEVICE_REAL Gt1   = *(dptr_u + pp + bssn::VAR::U_GT1 * sz_per_dof + offset);
                const DEVICE_REAL Gt2   = *(dptr_u + pp + bssn::VAR::U_GT2 * sz_per_dof + offset);
                const DEVICE_REAL B0    = *(dptr_u + pp + bssn::VAR::U_B0 * sz_per_dof + offset);
                const DEVICE_REAL B1    = *(dptr_u + pp + bssn::VAR::U_B1 * sz_per_dof + offset);
                const DEVICE_REAL B2    = *(dptr_u + pp + bssn::VAR::U_B2 * sz_per_dof + offset);
                
                #include "../scripts/bssnrhs_evar_derivs_zid.h"
                //#include "../scripts/bssneqs_eta_const_standard_gauge4.cpp"
                
            }

            g.sync();

            if(blk[BLK_ID].m_bflag != 0)
            {
                const DEVICE_REAL* const alpha = (dptr_u + bssn::VAR::U_ALPHA * sz_per_dof + offset);
                const DEVICE_REAL* const chi   = (dptr_u + bssn::VAR::U_CHI * sz_per_dof + offset);
                const DEVICE_REAL* const K     = (dptr_u + bssn::VAR::U_K * sz_per_dof + offset);
                const DEVICE_REAL* const gt0   = (dptr_u + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
                const DEVICE_REAL* const gt1   = (dptr_u + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
                const DEVICE_REAL* const gt2   = (dptr_u + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
                const DEVICE_REAL* const gt3   = (dptr_u + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
                const DEVICE_REAL* const gt4   = (dptr_u + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
                const DEVICE_REAL* const gt5   = (dptr_u + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
                const DEVICE_REAL* const beta0 = (dptr_u + bssn::VAR::U_BETA0 * sz_per_dof + offset);
                const DEVICE_REAL* const beta1 = (dptr_u + bssn::VAR::U_BETA1 * sz_per_dof + offset);
                const DEVICE_REAL* const beta2 = (dptr_u + bssn::VAR::U_BETA2 * sz_per_dof + offset);
                const DEVICE_REAL* const At0   = (dptr_u + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
                const DEVICE_REAL* const At1   = (dptr_u + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
                const DEVICE_REAL* const At2   = (dptr_u + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
                const DEVICE_REAL* const At3   = (dptr_u + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
                const DEVICE_REAL* const At4   = (dptr_u + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
                const DEVICE_REAL* const At5   = (dptr_u + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
                const DEVICE_REAL* const Gt0   = (dptr_u + bssn::VAR::U_GT0 * sz_per_dof + offset);
                const DEVICE_REAL* const Gt1   = (dptr_u + bssn::VAR::U_GT1 * sz_per_dof + offset);
                const DEVICE_REAL* const Gt2   = (dptr_u + bssn::VAR::U_GT2 * sz_per_dof + offset);
                const DEVICE_REAL* const B0    = (dptr_u + bssn::VAR::U_B0 * sz_per_dof + offset);
                const DEVICE_REAL* const B1    = (dptr_u + bssn::VAR::U_B1 * sz_per_dof + offset);
                const DEVICE_REAL* const B2    = (dptr_u + bssn::VAR::U_B2 * sz_per_dof + offset);

                DEVICE_REAL * const a_rhs    = &dptr_Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
                DEVICE_REAL * const chi_rhs  = &dptr_Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
                DEVICE_REAL * const gt_rhs00 = &dptr_Fu[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
                DEVICE_REAL * const gt_rhs01 = &dptr_Fu[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
                DEVICE_REAL * const gt_rhs02 = &dptr_Fu[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
                DEVICE_REAL * const gt_rhs11 = &dptr_Fu[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
                DEVICE_REAL * const gt_rhs12 = &dptr_Fu[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
                DEVICE_REAL * const gt_rhs22 = &dptr_Fu[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
                DEVICE_REAL * const b_rhs0   = &dptr_Fu[bssn::VAR::U_BETA0 * sz_per_dof + offset];
                DEVICE_REAL * const b_rhs1   = &dptr_Fu[bssn::VAR::U_BETA1 * sz_per_dof + offset];
                DEVICE_REAL * const b_rhs2   = &dptr_Fu[bssn::VAR::U_BETA2 * sz_per_dof + offset];
                DEVICE_REAL * const At_rhs00 = &dptr_Fu[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
                DEVICE_REAL * const At_rhs01 = &dptr_Fu[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
                DEVICE_REAL * const At_rhs02 = &dptr_Fu[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
                DEVICE_REAL * const At_rhs11 = &dptr_Fu[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
                DEVICE_REAL * const At_rhs12 = &dptr_Fu[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
                DEVICE_REAL * const At_rhs22 = &dptr_Fu[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
                DEVICE_REAL * const Gt_rhs0  = &dptr_Fu[bssn::VAR::U_GT0 * sz_per_dof + offset];
                DEVICE_REAL * const Gt_rhs1  = &dptr_Fu[bssn::VAR::U_GT1 * sz_per_dof + offset];
                DEVICE_REAL * const Gt_rhs2  = &dptr_Fu[bssn::VAR::U_GT2 * sz_per_dof + offset];
                DEVICE_REAL * const B_rhs0   = &dptr_Fu[bssn::VAR::U_B0 * sz_per_dof + offset];
                DEVICE_REAL * const B_rhs1   = &dptr_Fu[bssn::VAR::U_B1 * sz_per_dof + offset];
                DEVICE_REAL * const B_rhs2   = &dptr_Fu[bssn::VAR::U_B2 * sz_per_dof + offset];
                DEVICE_REAL * const K_rhs    = &dptr_Fu[bssn::VAR::U_K * sz_per_dof + offset];

                const DEVICE_REAL* const grad_0_alpha = (deriv_evars->grad_0_alpha + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_alpha = (deriv_evars->grad_1_alpha + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_alpha = (deriv_evars->grad_2_alpha + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_beta0 = (deriv_evars->grad_0_beta0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_beta0 = (deriv_evars->grad_1_beta0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_beta0 = (deriv_evars->grad_2_beta0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_beta1 = (deriv_evars->grad_0_beta1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_beta1 = (deriv_evars->grad_1_beta1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_beta1 = (deriv_evars->grad_2_beta1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_beta2 = (deriv_evars->grad_0_beta2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_beta2 = (deriv_evars->grad_1_beta2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_beta2 = (deriv_evars->grad_2_beta2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_B0    = (deriv_evars->grad_0_B0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_B0    = (deriv_evars->grad_1_B0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_B0    = (deriv_evars->grad_2_B0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_B1    = (deriv_evars->grad_0_B1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_B1    = (deriv_evars->grad_1_B1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_B1    = (deriv_evars->grad_2_B1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_B2    = (deriv_evars->grad_0_B2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_B2    = (deriv_evars->grad_1_B2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_B2    = (deriv_evars->grad_2_B2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_chi   = (deriv_evars->grad_0_chi + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_chi   = (deriv_evars->grad_1_chi + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_chi   = (deriv_evars->grad_2_chi + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_Gt0   = (deriv_evars->grad_0_Gt0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_Gt0   = (deriv_evars->grad_1_Gt0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_Gt0   = (deriv_evars->grad_2_Gt0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_Gt1   = (deriv_evars->grad_0_Gt1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_Gt1   = (deriv_evars->grad_1_Gt1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_Gt1   = (deriv_evars->grad_2_Gt1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_Gt2   = (deriv_evars->grad_0_Gt2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_Gt2   = (deriv_evars->grad_1_Gt2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_Gt2   = (deriv_evars->grad_2_Gt2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_K     = (deriv_evars->grad_0_K + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_K     = (deriv_evars->grad_1_K + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_K     = (deriv_evars->grad_2_K + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt0   = (deriv_evars->grad_0_gt0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt0   = (deriv_evars->grad_1_gt0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt0   = (deriv_evars->grad_2_gt0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt1   = (deriv_evars->grad_0_gt1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt1   = (deriv_evars->grad_1_gt1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt1   = (deriv_evars->grad_2_gt1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt2   = (deriv_evars->grad_0_gt2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt2   = (deriv_evars->grad_1_gt2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt2   = (deriv_evars->grad_2_gt2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt3   = (deriv_evars->grad_0_gt3 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt3   = (deriv_evars->grad_1_gt3 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt3   = (deriv_evars->grad_2_gt3 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt4   = (deriv_evars->grad_0_gt4 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt4   = (deriv_evars->grad_1_gt4 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt4   = (deriv_evars->grad_2_gt4 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt5   = (deriv_evars->grad_0_gt5 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt5   = (deriv_evars->grad_1_gt5 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt5   = (deriv_evars->grad_2_gt5 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At0   = (deriv_evars->grad_0_At0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At0   = (deriv_evars->grad_1_At0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At0   = (deriv_evars->grad_2_At0 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At1   = (deriv_evars->grad_0_At1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At1   = (deriv_evars->grad_1_At1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At1   = (deriv_evars->grad_2_At1 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At2   = (deriv_evars->grad_0_At2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At2   = (deriv_evars->grad_1_At2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At2   = (deriv_evars->grad_2_At2 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At3   = (deriv_evars->grad_0_At3 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At3   = (deriv_evars->grad_1_At3 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At3   = (deriv_evars->grad_2_At3 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At4   = (deriv_evars->grad_0_At4 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At4   = (deriv_evars->grad_1_At4 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At4   = (deriv_evars->grad_2_At4 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At5   = (deriv_evars->grad_0_At5 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At5   = (deriv_evars->grad_1_At5 + BLK_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At5   = (deriv_evars->grad_2_At5 + BLK_ID * BLK_SZ );
                
                radiative_bc<pw>(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk,BLK_ID);
                radiative_bc<pw>(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk, BLK_ID);
                radiative_bc<pw>(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, 1.0, 0.0, blk, BLK_ID);

                radiative_bc<pw>(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk, BLK_ID);

                radiative_bc<pw>(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk, BLK_ID);

                radiative_bc<pw>(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk, BLK_ID);

                radiative_bc<pw>(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk, BLK_ID);

                radiative_bc<pw>(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk, BLK_ID);
                radiative_bc<pw>(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk, BLK_ID);
                radiative_bc<pw>(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk, BLK_ID);
                radiative_bc<pw>(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk, BLK_ID);

            }


            g.sync();

            {
                DEVICE_REAL * const Fu = dptr_Fu;
                #include "../scripts/rhs_kernel_setup.cuh"
                #include "../scripts/bssnrhs_evar_derivs_zid.h"

                const DEVICE_REAL sigma = ko_sigma;

                a_rhs[pp]    += sigma * (kograd_0_alpha + kograd_1_alpha + kograd_2_alpha);
                
                b_rhs0[pp]   += sigma * (kograd_0_beta0 + kograd_1_beta0 + kograd_2_beta0);
                b_rhs1[pp]   += sigma * (kograd_0_beta1 + kograd_1_beta1 + kograd_2_beta1);
                b_rhs2[pp]   += sigma * (kograd_0_beta2 + kograd_1_beta2 + kograd_2_beta2);

                gt_rhs00[pp] += sigma * (kograd_0_gt0 + kograd_1_gt0 + kograd_2_gt0);
                gt_rhs01[pp] += sigma * (kograd_0_gt1 + kograd_1_gt1 + kograd_2_gt1);
                gt_rhs02[pp] += sigma * (kograd_0_gt2 + kograd_1_gt2 + kograd_2_gt2);
                gt_rhs11[pp] += sigma * (kograd_0_gt3 + kograd_1_gt3 + kograd_2_gt3);
                gt_rhs12[pp] += sigma * (kograd_0_gt4 + kograd_1_gt4 + kograd_2_gt4);
                gt_rhs22[pp] += sigma * (kograd_0_gt5 + kograd_1_gt5 + kograd_2_gt5);

                chi_rhs[pp]  += sigma * (kograd_0_chi + kograd_1_chi + kograd_2_chi);

                At_rhs00[pp] += sigma * (kograd_0_At0 + kograd_1_At0 + kograd_2_At0);
                At_rhs01[pp] += sigma * (kograd_0_At1 + kograd_1_At1 + kograd_2_At1);
                At_rhs02[pp] += sigma * (kograd_0_At2 + kograd_1_At2 + kograd_2_At2);
                At_rhs11[pp] += sigma * (kograd_0_At3 + kograd_1_At3 + kograd_2_At3);
                At_rhs12[pp] += sigma * (kograd_0_At4 + kograd_1_At4 + kograd_2_At4);
                At_rhs22[pp] += sigma * (kograd_0_At5 + kograd_1_At5 + kograd_2_At5);

                K_rhs[pp]    += sigma * (kograd_0_K + kograd_1_K + kograd_2_K);

                Gt_rhs0[pp]  += sigma * (kograd_0_Gt0 + kograd_1_Gt0 + kograd_2_Gt0);
                Gt_rhs1[pp]  += sigma * (kograd_0_Gt1 + kograd_1_Gt1 + kograd_2_Gt1);
                Gt_rhs2[pp]  += sigma * (kograd_0_Gt2 + kograd_1_Gt2 + kograd_2_Gt2);

                B_rhs0[pp]   += sigma * (kograd_0_B0 + kograd_1_B0 + kograd_2_B0);
                B_rhs1[pp]   += sigma * (kograd_0_B1 + kograd_1_B1 + kograd_2_B1);
                B_rhs2[pp]   += sigma * (kograd_0_B2 + kograd_1_B2 + kograd_2_B2);
            }
            
            g.sync();
            

        }


        
        

    }
