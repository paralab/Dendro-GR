#pragma once
#include "block_gpu.h"
#include "derivs.h"
#include "grDef.h"
#include "grUtils.h"
#include "rhs.h"

template<typename T>
T* create_block(BlockGPU3D* blk, unsigned int nblocks, unsigned int dof=1)
{   
    if(nblocks==0)
      return NULL;

    unsigned int szp_dof = (blk[nblocks-1].m_offset +  blk[nblocks-1].m_aligned_sz[0] * blk[nblocks-1].m_aligned_sz[1] * blk[nblocks-1].m_aligned_sz[2]);
    const unsigned int blk_vec_size = dof * szp_dof;
    T* blk_data = new T[blk_vec_size];
    
    #pragma omp parallel for default(none) shared(blk,nblocks,blk_data,dof, szp_dof) schedule(static)
    for(unsigned int blk_id=0; blk_id < nblocks; blk_id++)
    {

      const unsigned int nx     = blk[blk_id].m_aligned_sz[0];
      const unsigned int ny     = blk[blk_id].m_aligned_sz[1];
      const unsigned int nz     = blk[blk_id].m_aligned_sz[2];
      const unsigned int offset = blk[blk_id].m_offset;
      T var[dof];

      for(unsigned int k=0; k < nz; k++)
      for(unsigned int j=0; j < ny; j++)
      for(unsigned int i=0; i < nx; i++)
      {
          const double x = blk[blk_id].m_ptMin[0] + i * blk[blk_id].m_dx[0];
          const double y = blk[blk_id].m_ptMin[1] + j * blk[blk_id].m_dx[1];
          const double z = blk[blk_id].m_ptMin[2] + k * blk[blk_id].m_dx[2];
          bssn::fake_initial_data(x,y,z,var);
          
          for(unsigned int v=0; v < dof; v++)
            blk_data[ v * szp_dof + offset +  k * ny * nx + j * nx + i] = var[v];
      }
    }
    return blk_data;
}

template<int pw, int pencils, int pencil_sz>
void bssnrhs_cpu(DEVICE_REAL * const Fu, const DEVICE_REAL * const u, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_REAL* const deriv_workspace, DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma=0.0)
{   
    #pragma omp parallel for shared(blk,nblocks,szpdof_uz,ko_sigma) schedule(static)
    for (unsigned int BLK_ID =0; BLK_ID < nblocks; BLK_ID ++)
    {

        const DEVICE_UINT BATCHED_BLOCKS_SZ = 1;
        const DEVICE_UINT BLK_SZ            = (blk[BLK_ID].m_aligned_sz[0]) * (blk[BLK_ID].m_aligned_sz[1]) * (blk[BLK_ID].m_aligned_sz[2]); 
        const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
        const DEVICE_UINT sz_per_dof        = szpdof_uz;
        const DEVICE_UINT bflag             = blk[BLK_ID].m_bflag;
        const DEVICE_UINT sz[3]             = {blk[BLK_ID].m_aligned_sz[0],blk[BLK_ID].m_aligned_sz[1],blk[BLK_ID].m_aligned_sz[2]};
        //printf("sz: %d,%d,%d\n",sz[0],sz[1],sz[2]);
        
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


        DEVICE_REAL * const a_rhs    = &Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
        DEVICE_REAL * const chi_rhs  = &Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
        DEVICE_REAL * const K_rhs    = &Fu[bssn::VAR::U_K * sz_per_dof + offset];
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

        const unsigned int nx = blk[BLK_ID].m_aligned_sz[0];
        const unsigned int ny = blk[BLK_ID].m_aligned_sz[1];
        const unsigned int nz = blk[BLK_ID].m_aligned_sz[2];

        const unsigned int a_nx = blk[BLK_ID].m_sz[0];
        const unsigned int a_ny = blk[BLK_ID].m_sz[1];
        const unsigned int a_nz = blk[BLK_ID].m_sz[2];

        const DEVICE_REAL hx       = blk[BLK_ID].m_dx[0]; 
        const DEVICE_REAL hy       = blk[BLK_ID].m_dx[1]; 
        const DEVICE_REAL hz       = blk[BLK_ID].m_dx[2]; 

        const DEVICE_UINT lambda[4]     = {1, 1, 1, 1};
        const DEVICE_REAL lambda_f[2]   = {1.0, 0.0};
        
        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL =  2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;

        const unsigned int PW=pw;
        unsigned int n = nx*ny*nz;

        const unsigned int tid = omp_get_thread_num();
        double * const deriv_base =  deriv_workspace + tid * 210 * BLK_SZ;

        #include "../scripts/bssnrhs_memalloc.h"
        #include "bssnrhs_derivs.h"
        
        for (unsigned int k = PW; k < a_nz-PW; k++) {
        for (unsigned int j = PW; j < a_ny-PW; j++) {
            #ifdef BSSN_ENABLE_AVX
                #ifdef __INTEL_COMPILER
                #pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
                #pragma ivdep
                #endif
            #endif
            for (unsigned int i = PW; i < a_nx-PW; i++)
            {
                const DEVICE_REAL x = blk->m_ptMin[0] + i*hx;
                const DEVICE_REAL y = blk->m_ptMin[1] + j*hy;
                const DEVICE_REAL z = blk->m_ptMin[2] + k*hz;

                const unsigned int pp = i + nx*(j + ny*k);
                const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

                const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
                const DEVICE_REAL arg = - w*w*w*w;
                const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;

                #include "../../src/bssneqs_eta_const_standard_gauge.cpp"

            }
        }
        }

        if (bflag != 0) {

            double *pmin  = blk->m_ptMin;
            double  pmax[3];
            pmax[0]  = blk->m_ptMin[0] + (a_nx-1) * hx;
            pmax[1]  = blk->m_ptMin[1] + (a_ny-1) * hy;
            pmax[2]  = blk->m_ptMin[2] + (a_nz-1) * hz;
            
            bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin, pmax,
                    1.0, 1.0, sz, bflag);
            bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                    1.0, 1.0, sz, bflag);
            bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax,
                    1.0, 0.0, sz, bflag);

            bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin, pmax,
                    1.0, 0.0, sz, bflag);
            bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin, pmax,
                    1.0, 0.0, sz, bflag);
            bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin, pmax,
                    1.0, 0.0, sz, bflag);

            bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                    2.0, 0.0, sz, bflag);
            bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                    2.0, 0.0, sz, bflag);
            bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                    2.0, 0.0, sz, bflag);

            bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax,
                    1.0, 0.0, sz, bflag);
            bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax,
                    1.0, 0.0, sz, bflag);
            bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax,
                    1.0, 0.0, sz, bflag);

            bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                    2.0, 0.0, sz, bflag);
            bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                    2.0, 0.0, sz, bflag);
            bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                    2.0, 0.0, sz, bflag);
            bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                    2.0, 0.0, sz, bflag);
            bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                    2.0, 0.0, sz, bflag);
            bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                    2.0, 0.0, sz, bflag);

            bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                    1.0, 1.0, sz, bflag);
            bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                    1.0, 0.0, sz, bflag);
            bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                    1.0, 0.0, sz, bflag);
            bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                    1.0, 1.0, sz, bflag);
            bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                    1.0, 0.0, sz, bflag);
            bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                    1.0, 1.0, sz, bflag);
        
        }

        #include "../../BSSN_GR/scripts/bssnrhs_ko_derivs.h"


        const  double sigma = ko_sigma;


        for (unsigned int k = PW; k < a_nz-PW; k++) {
            for (unsigned int j = PW; j < a_ny-PW; j++) {
            #ifdef BSSN_ENABLE_AVX
                #ifdef __INTEL_COMPILER
                #pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
                #pragma ivdep
                #endif
            #endif
            for (unsigned int i = PW; i < a_nx-PW; i++) {
                const unsigned int pp = i + nx*(j + ny*k);

                a_rhs[pp]  += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] + grad_2_alpha[pp]);
                b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] + grad_2_beta0[pp]);
                b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] + grad_2_beta1[pp]);
                b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] + grad_2_beta2[pp]);

                gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                chi_rhs[pp]  += sigma * (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                K_rhs[pp] += sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                B_rhs0[pp] += sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                B_rhs1[pp] += sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                B_rhs2[pp] += sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
              }
            }
        }
    }

    return;

}