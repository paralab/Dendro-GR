/**
 * @file bssn_deriv_kernels.cuh
 * @brief BSSN RHS derivative kernels. 
 * @version 0.1
 * @date 2021-12-18
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#pragma once
#include "device.h"
#include "block_gpu.h"
#include "grDef.h"
#include "bssnrhs_evar_derivs.cuh"
#include "derivs_cu.cuh"
#include "bc_cu.cuh"
#include <cuda_runtime_api.h> 
#include <cuda.h> 
#include <cooperative_groups.h>

namespace device
{
    /**
     * @brief Launch the kernel to compute all the x-direction derivatives
    */
    template<int pw, int pencils, int pencil_sz, int BATCHED_BLOCKS_SZ>
    GLOBAL_FUNC void launch_dir_x_deriv_kernel(const DEVICE_REAL * const  u,BSSN_EVAR_DERIVS* deriv_evars , BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz) {
        
        //includes x-derivs here.
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
        
        #include "../scripts/bssnrhs_deriv_x_dir.cuh"

        return;

    }

    /**
     * @brief Launch the kernel to compute all the y-direction derivatives
    */
    template<int pw, int pencils, int pencil_sz, int BATCHED_BLOCKS_SZ>
    GLOBAL_FUNC void launch_dir_y_deriv_kernel(const DEVICE_REAL * const  u, BSSN_EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz) {
        
        //includes y-derivs here.
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

        #include "../scripts/bssnrhs_deriv_y_dir.cuh"

        return;

    }

    template<typename T>
    GLOBAL_FUNC 
    __launch_bounds__(1024)
    void cuda_bssn_enforce_evar_cons(T* u, unsigned int lb, unsigned int le, T chi_floor, unsigned int szpdof)
    {   

        const T CHI_FLOOR = chi_floor;
        const T one_third = 1.0 / 3.0;
        T gtd[3][3], Atd[3][3];
        T* uiVar[bssn::BSSN_NUM_VARS];
        for( int i=0; i < bssn::BSSN_NUM_VARS; i++)
            uiVar[i] = u + i * szpdof;

        const int t_id   = threadIdx.x + blockDim.x * blockIdx.x;
        const int stride = blockDim.x * gridDim.x;
        for (int node = (t_id + lb); node < le; node += stride)
        {
            gtd[0][0] = uiVar[bssn::VAR::U_SYMGT0][node];
            gtd[0][1] = uiVar[bssn::VAR::U_SYMGT1][node];
            gtd[0][2] = uiVar[bssn::VAR::U_SYMGT2][node];
            gtd[1][0] = gtd[0][1];
            gtd[1][1] = uiVar[bssn::VAR::U_SYMGT3][node];
            gtd[1][2] = uiVar[bssn::VAR::U_SYMGT4][node];
            gtd[2][0] = gtd[0][2];
            gtd[2][1] = gtd[1][2];
            gtd[2][2] = uiVar[bssn::VAR::U_SYMGT5][node];

            Atd[0][0] = uiVar[bssn::VAR::U_SYMAT0][node];
            Atd[0][1] = uiVar[bssn::VAR::U_SYMAT1][node];
            Atd[0][2] = uiVar[bssn::VAR::U_SYMAT2][node];
            Atd[1][0] = Atd[0][1];
            Atd[1][1] = uiVar[bssn::VAR::U_SYMAT3][node];
            Atd[1][2] = uiVar[bssn::VAR::U_SYMAT4][node];
            Atd[2][0] = Atd[0][2];
            Atd[2][1] = Atd[1][2];
            Atd[2][2] = uiVar[bssn::VAR::U_SYMAT5][node];

            //printf("node %d\n gt00 = %.8E\n", node, gtd[0][0]);

            T det_gtd  =  gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
                            -  gtd[0][1]*gtd[0][1]*gtd[2][2]
                            +  2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
                            -  gtd[0][2]*gtd[0][2]*gtd[1][1];

            if (det_gtd < 0.0) {
                printf("metric determinent = %.8E is negative \n",det_gtd);
                assert(false);
                //exit(0);
                /* FIXME What to do here? The metric is not physical. Do we reset the metric to be flat? */
                gtd[0][0] = 1.0; gtd[0][1] = 0.0; gtd[0][2] = 0.0;
                gtd[1][0] = 0.0; gtd[1][1] = 1.0; gtd[1][2] = 0.0;
                gtd[2][0] = 0.0; gtd[2][1] = 0.0; gtd[2][2] = 1.0;
                det_gtd = 1.0;
            }
            
            T det_gtd_to_neg_third = 1.0 / pow(det_gtd, one_third);

            for (unsigned int j = 0; j < 3; j++)
                for (unsigned int i = 0; i < 3; i++)
                {
                    gtd[i][j] *= det_gtd_to_neg_third;
                }
            
            det_gtd =   gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
                        - gtd[0][1]*gtd[0][1]*gtd[2][2]
                        + 2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
                        - gtd[0][2]*gtd[0][2]*gtd[1][1];

            T detgt_m1 = det_gtd - 1.0;

            if (fabs(detgt_m1) > 1.0e-6) {
                printf("enforce_bssn_constraint: det(gtd) != 1, det(gtd)=%.8E\n", det_gtd);
                assert(false);
                // //std::cout.precision(14);
                // printf"enforce_bssn_constraint: det(gtd) != 1. det="<<std::fixed<<det_gtd<<std::endl;
                //     std::cout<<"      gtd(1,1)="<<gtd[0][0]<<std::endl;
                //     std::cout<<"      gtd(1,2)="<<gtd[0][1]<<std::endl;
                //     std::cout<<"      gtd(1,3)="<<gtd[0][2]<<std::endl;
                //     std::cout<<"      gtd(2,2)="<<gtd[1][1]<<std::endl;
                //     std::cout<<"      gtd(2,3)="<<gtd[1][2]<<std::endl;
                //     std::cout<<"      gtd(3,3)="<<gtd[2][2]<<std::endl;
                
                //exit(0);
                
            }

            T gtu[3][3];
            T idet_gtd = 1.0/det_gtd;
            gtu[0][0] = idet_gtd*(gtd[1][1]*gtd[2][2]-gtd[1][2]*gtd[1][2]);
            gtu[0][1] = idet_gtd*(-gtd[0][1]*gtd[2][2]+gtd[0][2]*gtd[1][2]);
            gtu[0][2] = idet_gtd*(gtd[0][1]*gtd[1][2]-gtd[0][2]*gtd[1][1]);
            gtu[1][0] = gtu[0][1];
            gtu[1][1] = idet_gtd*(gtd[0][0]*gtd[2][2]-gtd[0][2]*gtd[0][2]);
            gtu[1][2] = idet_gtd*(-gtd[0][0]*gtd[1][2]+gtd[0][1]*gtd[0][2]);
            gtu[2][0] = gtu[0][2];
            gtu[2][1] = gtu[1][2];
            gtu[2][2] = idet_gtd*(gtd[0][0]*gtd[1][1]-gtd[0][1]*gtd[0][1]);

            /* Require Atd to be traceless. */
            T one_third_trace_Atd =   one_third * (
                                    Atd[0][0]*gtu[0][0]
                                + Atd[1][1]*gtu[1][1]
                                + Atd[2][2]*gtu[2][2]
                                + 2.0 * (   Atd[0][1]*gtu[0][1]
                                            + Atd[0][2]*gtu[0][2]
                                            + Atd[1][2]*gtu[1][2]  )
                                );

            Atd[0][0] -= one_third_trace_Atd * gtd[0][0];
            Atd[0][1] -= one_third_trace_Atd * gtd[0][1];
            Atd[0][2] -= one_third_trace_Atd * gtd[0][2];
            Atd[1][1] -= one_third_trace_Atd * gtd[1][1];
            Atd[1][2] -= one_third_trace_Atd * gtd[1][2];
            Atd[2][2] -= one_third_trace_Atd * gtd[2][2];

            T tr_A =    Atd[0][0]*gtu[0][0]
                            + Atd[1][1]*gtu[1][1]
                            + Atd[2][2]*gtu[2][2]
                            + 2.0 * (   Atd[0][1]*gtu[0][1]
                                    + Atd[0][2]*gtu[0][2]
                                    + Atd[1][2]*gtu[1][2]  );


            if (fabs(tr_A) > 1.0e-6) {
                printf("trace(A) != 0, trace(A)=%.8E\n", tr_A);
                assert(false);
                // std::cout<<"enforce_bssn_constraint: tr_A != 0. tr_A="<<tr_A<<std::endl;
                // std::cout<<"      Atd(1,1)="<<Atd[0][0]<<std::endl;
                // std::cout<<"      Atd(1,2)="<<Atd[0][1]<<std::endl;
                // std::cout<<"      Atd(1,3)="<<Atd[0][2]<<std::endl;
                // std::cout<<"      Atd(2,2)="<<Atd[1][1]<<std::endl;
                // std::cout<<"      Atd(2,3)="<<Atd[1][2]<<std::endl;
                // std::cout<<"      Atd(3,3)="<<Atd[2][2]<<std::endl;
                
                //exit(0);
            }


            uiVar[bssn::VAR::U_SYMAT0][node] = Atd[0][0];
            uiVar[bssn::VAR::U_SYMAT1][node] = Atd[0][1];
            uiVar[bssn::VAR::U_SYMAT2][node] = Atd[0][2];
            uiVar[bssn::VAR::U_SYMAT3][node] = Atd[1][1];
            uiVar[bssn::VAR::U_SYMAT4][node] = Atd[1][2];
            uiVar[bssn::VAR::U_SYMAT5][node] = Atd[2][2];

            uiVar[bssn::VAR::U_SYMGT0][node] = gtd[0][0];
            uiVar[bssn::VAR::U_SYMGT1][node] = gtd[0][1];
            uiVar[bssn::VAR::U_SYMGT2][node] = gtd[0][2];
            uiVar[bssn::VAR::U_SYMGT3][node] = gtd[1][1];
            uiVar[bssn::VAR::U_SYMGT4][node] = gtd[1][2];
            uiVar[bssn::VAR::U_SYMGT5][node] = gtd[2][2];

            /* apply a floor to chi */
            if ( uiVar[bssn::VAR::U_CHI][node] < CHI_FLOOR ) {
            /* FIXME This needs to be fixed when we add a fluid to the code. */
            /* ! First rescale the densitized fluid variables.
                ! The include file bssn_puncture_fluid_rescale.inc
                ! must be provided in the BSSN_*MHD project.

                ! Chi must be positive to do the rescaling of fluid variables.
                if ( chi <= 0.0) {
                    chi = pars.chi_floor;
                }
                else {
                    // ok... go ahead and rescale the fluid variables.
                }


                */

            /* now place the floor on chi */
                uiVar[bssn::VAR::U_CHI][node] = CHI_FLOOR;
            }

            /* apply a floor to alpha */
            uiVar[bssn::VAR::U_ALPHA][node] = max(uiVar[bssn::VAR::U_ALPHA][node], CHI_FLOOR);
            }


    }


    
    template<typename T,unsigned int pw, unsigned int nx>
    DEVICE_FUNC
    __inline__ void
    __load_blk_var__(T * s_u, const T* const var, const BlockGPU3D* blk)
    {
        const DEVICE_INT i  = threadIdx.x;//GPUDevice::thread_id_x();
        const DEVICE_INT j  = threadIdx.y;//GPUDevice::thread_id_y();
        
        //const DEVICE_INT nx   = //blk->m_aligned_sz[0];
        //const DEVICE_INT ny   = //blk->m_aligned_sz[1];
        //const DEVICE_INT nz   = //blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4* pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4* pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4* pw + 1; //blk->m_sz[2];

        const DEVICE_UINT globalIdx = j * nx + i;
        const DEVICE_INT  ss = nx*nx;

        if(i < actual_nx && j < actual_ny)
        for(DEVICE_INT k = 0; k <  actual_nz; k++)
        {
            int idx = globalIdx + k * ss;
            s_u[idx]              = var[idx];
        }
            
        
        GPUDevice::sync_threads();
        return ; 
    }


    template<typename T,unsigned int pw, unsigned int nx>
    DEVICE_FUNC
    __inline__ void
    __ld_blk_var1__(T * s_u, const T* const var, const BlockGPU3D* blk)
    {
        const DEVICE_INT i   = GPUDevice::thread_id_x();
        const DEVICE_INT j   = GPUDevice::thread_id_y();
        const DEVICE_INT k   = GPUDevice::thread_id_z();

        const DEVICE_INT tid = (k * GPUDevice::block_dim_x()  + j ) * GPUDevice::block_dim_y()  + i;

        const DEVICE_INT actual_nx = 4* pw + 1;
        const DEVICE_INT actual_ny = 4* pw + 1;
        const DEVICE_INT actual_nz = 4* pw + 1;

        const DEVICE_INT n    = nx * nx * nx;
        const DEVICE_INT inc  = GPUDevice::block_dim_x() * GPUDevice::block_dim_y() * GPUDevice::block_dim_z();
        
        for(DEVICE_INT k = tid; k <  n; k+=inc)
            s_u[k]              = var[k];
            
        
        GPUDevice::sync_threads();
        return ; 
    }

    template<typename T,unsigned int pw, unsigned int nx>
    DEVICE_FUNC
    __inline__ void
    __st_blk_var1__(const T* const s_u, T* __restrict__ const var, const BlockGPU3D* blk)
    {
        const DEVICE_INT i   = GPUDevice::thread_id_x();
        const DEVICE_INT j   = GPUDevice::thread_id_y();
        const DEVICE_INT k   = GPUDevice::thread_id_z();

        const DEVICE_INT tid = (k * GPUDevice::block_dim_x()  + j ) * GPUDevice::block_dim_y()  + i;

        const DEVICE_INT actual_nx = 4* pw + 1;
        const DEVICE_INT actual_ny = 4* pw + 1;
        const DEVICE_INT actual_nz = 4* pw + 1;

        const DEVICE_INT n    = nx * nx * nx;
        const DEVICE_INT inc  = GPUDevice::block_dim_x() * GPUDevice::block_dim_y() * GPUDevice::block_dim_z();
        
        for(DEVICE_INT k = tid; k <  n; k+=inc)
            var[k]=s_u[k];
            
        
        GPUDevice::sync_threads();
        return ; 
    }

    
    template<int pw, int nx, int pencil_sz>
    DEVICE_FUNC
    void eval_rhs(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, const BlockGPU3D* blk, BSSN_EVAR_DERIVS* const deriv_evars,DEVICE_UINT szpdof_uz)
    {
        const DEVICE_UINT sz_per_dof   = szpdof_uz;
        const DEVICE_UINT Z_ID         = GPUDevice::block_id_x();
        const DEVICE_UINT si           = GPUDevice::thread_id_x() + pw;
        const DEVICE_UINT sj           = GPUDevice::thread_id_y() + pw;

        const DEVICE_UINT lambda[4]    = {1, 1, 1, 1};
        const DEVICE_REAL lambda_f[2]  = {1.0, 0.0};
        const DEVICE_UINT offset            = blk->m_offset; 
        const DEVICE_UINT actual_nx         = blk->m_sz[0];
        //const DEVICE_UINT nx              = blk[BLK_ID].m_aligned_sz[0];
        const DEVICE_UINT BLK_SZ            = nx*nx*nx; 

        DEVICE_REAL * __restrict__ const a_rhs    = &dptr_Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const chi_rhs  = &dptr_Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs00 = &dptr_Fu[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs01 = &dptr_Fu[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs02 = &dptr_Fu[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs11 = &dptr_Fu[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs12 = &dptr_Fu[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs22 = &dptr_Fu[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs0   = &dptr_Fu[bssn::VAR::U_BETA0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs1   = &dptr_Fu[bssn::VAR::U_BETA1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs2   = &dptr_Fu[bssn::VAR::U_BETA2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs00 = &dptr_Fu[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs01 = &dptr_Fu[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs02 = &dptr_Fu[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs11 = &dptr_Fu[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs12 = &dptr_Fu[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs22 = &dptr_Fu[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs0  = &dptr_Fu[bssn::VAR::U_GT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs1  = &dptr_Fu[bssn::VAR::U_GT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs2  = &dptr_Fu[bssn::VAR::U_GT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs0   = &dptr_Fu[bssn::VAR::U_B0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs1   = &dptr_Fu[bssn::VAR::U_B1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs2   = &dptr_Fu[bssn::VAR::U_B2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const K_rhs    = &dptr_Fu[bssn::VAR::U_K * sz_per_dof + offset];

        const DEVICE_REAL * __restrict__ const alpha = (dptr_u + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const chi   = (dptr_u + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const K     = (dptr_u + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const gt0   = (dptr_u + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const gt1   = (dptr_u + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const gt2   = (dptr_u + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const gt3   = (dptr_u + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const gt4   = (dptr_u + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const gt5   = (dptr_u + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const beta0 = (dptr_u + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const beta1 = (dptr_u + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const beta2 = (dptr_u + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const At0   = (dptr_u + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const At1   = (dptr_u + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const At2   = (dptr_u + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const At3   = (dptr_u + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const At4   = (dptr_u + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const At5   = (dptr_u + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const Gt0   = (dptr_u + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const Gt1   = (dptr_u + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const Gt2   = (dptr_u + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const B0    = (dptr_u + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const B1    = (dptr_u + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ const B2    = (dptr_u + bssn::VAR::U_B2 * sz_per_dof + offset);

        const DEVICE_REAL hx = blk->m_dx[0];
        const DEVICE_REAL hy = blk->m_dx[1];
        const DEVICE_REAL hz = blk->m_dx[2];
        const DEVICE_REAL x = blk->m_ptMin[0] + si*hx;
        const DEVICE_REAL y = blk->m_ptMin[1] + sj*hy;

        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;

        
        if(si < nx-pw && sj < nx-pw)
        for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
        {
            const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
            const DEVICE_REAL z = blk->m_ptMin[2] + sk*hz;

            const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

            const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
            const DEVICE_REAL arg = - w*w*w*w;
            const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;

            const DEVICE_UINT d_pp = Z_ID * BLK_SZ + pp;
            #include "../scripts/bssneqs_eta_const_standard_gauge4.cpp"
        
        }

    }


    template<int pw, int pencils, int pencil_sz, unsigned int BATCHED_GRAIN_SZ, unsigned int nx>
    GLOBAL_FUNC 
    __launch_bounds__(256)
    void bssnrhs(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const dptr_deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {

        const DEVICE_UINT NUM_BATCHES  = nblocks/BATCHED_GRAIN_SZ + 1;
        const DEVICE_UINT blk_stride   = GPUDevice::grid_dim_x(); 
        const DEVICE_UINT Z_ID         = GPUDevice::block_id_x();
        const DEVICE_UINT block_begin  = Z_ID;
        const DEVICE_UINT sz_per_dof   = szpdof_uz;

        SHARED_MEM DEVICE_REAL su[nx*nx*nx];
        BSSN_EVAR_DERIVS * const deriv_evars = dptr_deriv_evars;
        const DEVICE_UINT si = GPUDevice::thread_id_x() + pw;
        const DEVICE_UINT sj = GPUDevice::thread_id_y() + pw;

        #pragma unroll
        for(unsigned int bid = block_begin ; bid < nblocks; bid+=blk_stride )
        {
            const BlockGPU3D current_blk         = dptr_blk[bid];
            const BlockGPU3D * const blk         =  &current_blk;
            // if(si ==pw && sj ==pw)
            //     printf("blk %d/%d with z_id %d\n",bid,nblocks,Z_ID);

            //const DEVICE_UINT BLK_ID            = bid;
            const DEVICE_UINT offset            = blk->m_offset; 
            const DEVICE_UINT actual_nx         = blk->m_sz[0];
            //const DEVICE_UINT nx              = blk[BLK_ID].m_aligned_sz[0];
            const DEVICE_UINT BLK_SZ            = nx*nx*nx; 
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

                #include "../scripts/compute_bssnrhs_evar_derivs.cuh"
            }
            
            GPUDevice::sync_threads();
            
            //GPUDevice::sync_threads();
            // if(si ==pw && sj ==pw)
            //     printf("blk %d/%d with z_id %d\n : derivs ended",bid,nblocks,Z_ID);

            eval_rhs<pw, nx, pencil_sz> (dptr_Fu, dptr_u, blk, deriv_evars, szpdof_uz);
            
            //GPUDevice::sync_threads();
            // if(si ==pw && sj ==pw)
            //     printf("blk %d/%d with z_id %d\n : rhs ended",bid,nblocks,Z_ID);

            if(blk->m_bflag != 0)
            {
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

                const DEVICE_REAL* const grad_0_alpha = (deriv_evars->grad_0_alpha + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_alpha = (deriv_evars->grad_1_alpha + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_alpha = (deriv_evars->grad_2_alpha + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_beta0 = (deriv_evars->grad_0_beta0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_beta0 = (deriv_evars->grad_1_beta0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_beta0 = (deriv_evars->grad_2_beta0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_beta1 = (deriv_evars->grad_0_beta1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_beta1 = (deriv_evars->grad_1_beta1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_beta1 = (deriv_evars->grad_2_beta1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_beta2 = (deriv_evars->grad_0_beta2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_beta2 = (deriv_evars->grad_1_beta2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_beta2 = (deriv_evars->grad_2_beta2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_B0    = (deriv_evars->grad_0_B0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_B0    = (deriv_evars->grad_1_B0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_B0    = (deriv_evars->grad_2_B0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_B1    = (deriv_evars->grad_0_B1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_B1    = (deriv_evars->grad_1_B1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_B1    = (deriv_evars->grad_2_B1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_B2    = (deriv_evars->grad_0_B2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_B2    = (deriv_evars->grad_1_B2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_B2    = (deriv_evars->grad_2_B2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_chi   = (deriv_evars->grad_0_chi + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_chi   = (deriv_evars->grad_1_chi + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_chi   = (deriv_evars->grad_2_chi + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_Gt0   = (deriv_evars->grad_0_Gt0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_Gt0   = (deriv_evars->grad_1_Gt0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_Gt0   = (deriv_evars->grad_2_Gt0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_Gt1   = (deriv_evars->grad_0_Gt1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_Gt1   = (deriv_evars->grad_1_Gt1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_Gt1   = (deriv_evars->grad_2_Gt1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_Gt2   = (deriv_evars->grad_0_Gt2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_Gt2   = (deriv_evars->grad_1_Gt2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_Gt2   = (deriv_evars->grad_2_Gt2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_K     = (deriv_evars->grad_0_K + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_K     = (deriv_evars->grad_1_K + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_K     = (deriv_evars->grad_2_K + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt0   = (deriv_evars->grad_0_gt0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt0   = (deriv_evars->grad_1_gt0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt0   = (deriv_evars->grad_2_gt0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt1   = (deriv_evars->grad_0_gt1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt1   = (deriv_evars->grad_1_gt1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt1   = (deriv_evars->grad_2_gt1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt2   = (deriv_evars->grad_0_gt2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt2   = (deriv_evars->grad_1_gt2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt2   = (deriv_evars->grad_2_gt2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt3   = (deriv_evars->grad_0_gt3 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt3   = (deriv_evars->grad_1_gt3 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt3   = (deriv_evars->grad_2_gt3 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt4   = (deriv_evars->grad_0_gt4 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt4   = (deriv_evars->grad_1_gt4 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt4   = (deriv_evars->grad_2_gt4 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_gt5   = (deriv_evars->grad_0_gt5 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_gt5   = (deriv_evars->grad_1_gt5 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_gt5   = (deriv_evars->grad_2_gt5 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At0   = (deriv_evars->grad_0_At0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At0   = (deriv_evars->grad_1_At0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At0   = (deriv_evars->grad_2_At0 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At1   = (deriv_evars->grad_0_At1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At1   = (deriv_evars->grad_1_At1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At1   = (deriv_evars->grad_2_At1 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At2   = (deriv_evars->grad_0_At2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At2   = (deriv_evars->grad_1_At2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At2   = (deriv_evars->grad_2_At2 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At3   = (deriv_evars->grad_0_At3 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At3   = (deriv_evars->grad_1_At3 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At3   = (deriv_evars->grad_2_At3 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At4   = (deriv_evars->grad_0_At4 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At4   = (deriv_evars->grad_1_At4 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At4   = (deriv_evars->grad_2_At4 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_0_At5   = (deriv_evars->grad_0_At5 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_1_At5   = (deriv_evars->grad_1_At5 + Z_ID * BLK_SZ );
                const DEVICE_REAL* const grad_2_At5   = (deriv_evars->grad_2_At5 + Z_ID * BLK_SZ );
                
                radiative_bc_blk<pw>(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk,0);
                radiative_bc_blk<pw>(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk, 0);
                radiative_bc_blk<pw>(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, 1.0, 0.0, blk, 0);

                radiative_bc_blk<pw>(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk, 0);

                radiative_bc_blk<pw>(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk, 0);

                radiative_bc_blk<pw>(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk, 0);

                radiative_bc_blk<pw>(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk, 0);

                radiative_bc_blk<pw>(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk, 0);

            }

        
            const DEVICE_REAL sigma = ko_sigma;
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
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

                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                #include "../scripts/bssnrhs_evar_derivs_ko_zid.h"
                a_rhs[pp]    += sigma * (kograd_0_alpha[pp] + kograd_1_alpha[pp] + kograd_2_alpha[pp]);
        
                b_rhs0[pp]   += sigma * (kograd_0_beta0[pp] + kograd_1_beta0[pp] + kograd_2_beta0[pp]);
                b_rhs1[pp]   += sigma * (kograd_0_beta1[pp] + kograd_1_beta1[pp] + kograd_2_beta1[pp]);
                b_rhs2[pp]   += sigma * (kograd_0_beta2[pp] + kograd_1_beta2[pp] + kograd_2_beta2[pp]);

                gt_rhs00[pp] += sigma * (kograd_0_gt0[pp] + kograd_1_gt0[pp] + kograd_2_gt0[pp]);
                gt_rhs01[pp] += sigma * (kograd_0_gt1[pp] + kograd_1_gt1[pp] + kograd_2_gt1[pp]);
                gt_rhs02[pp] += sigma * (kograd_0_gt2[pp] + kograd_1_gt2[pp] + kograd_2_gt2[pp]);
                gt_rhs11[pp] += sigma * (kograd_0_gt3[pp] + kograd_1_gt3[pp] + kograd_2_gt3[pp]);
                gt_rhs12[pp] += sigma * (kograd_0_gt4[pp] + kograd_1_gt4[pp] + kograd_2_gt4[pp]);
                gt_rhs22[pp] += sigma * (kograd_0_gt5[pp] + kograd_1_gt5[pp] + kograd_2_gt5[pp]);

                chi_rhs[pp]  += sigma * (kograd_0_chi[pp] + kograd_1_chi[pp] + kograd_2_chi[pp]);

                At_rhs00[pp] += sigma * (kograd_0_At0[pp] + kograd_1_At0[pp] + kograd_2_At0[pp]);
                At_rhs01[pp] += sigma * (kograd_0_At1[pp] + kograd_1_At1[pp] + kograd_2_At1[pp]);
                At_rhs02[pp] += sigma * (kograd_0_At2[pp] + kograd_1_At2[pp] + kograd_2_At2[pp]);
                At_rhs11[pp] += sigma * (kograd_0_At3[pp] + kograd_1_At3[pp] + kograd_2_At3[pp]);
                At_rhs12[pp] += sigma * (kograd_0_At4[pp] + kograd_1_At4[pp] + kograd_2_At4[pp]);
                At_rhs22[pp] += sigma * (kograd_0_At5[pp] + kograd_1_At5[pp] + kograd_2_At5[pp]);

                K_rhs[pp]    += sigma * (kograd_0_K[pp] + kograd_1_K[pp] + kograd_2_K[pp]);

                Gt_rhs0[pp]  += sigma * (kograd_0_Gt0[pp] + kograd_1_Gt0[pp] + kograd_2_Gt0[pp]);
                Gt_rhs1[pp]  += sigma * (kograd_0_Gt1[pp] + kograd_1_Gt1[pp] + kograd_2_Gt1[pp]);
                Gt_rhs2[pp]  += sigma * (kograd_0_Gt2[pp] + kograd_1_Gt2[pp] + kograd_2_Gt2[pp]);

                B_rhs0[pp]   += sigma * (kograd_0_B0[pp] + kograd_1_B0[pp] + kograd_2_B0[pp]);
                B_rhs1[pp]   += sigma * (kograd_0_B1[pp] + kograd_1_B1[pp] + kograd_2_B1[pp]);
                B_rhs2[pp]   += sigma * (kograd_0_B2[pp] + kograd_1_B2[pp] + kograd_2_B2[pp]);
            }

            //GPUDevice::sync_threads();
            // if(si ==pw && sj ==pw)
            //     printf("blk %d/%d with z_id %d\n : ko  ended",bid,nblocks,Z_ID);

        }

        

    }


    template<int pw, unsigned int nx,  int pencil_sz>
    GLOBAL_FUNC 
    __launch_bounds__(256,2)
    void compute_all_derivs(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const dptr_deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {

        SHARED_MEM DEVICE_REAL su[(pencil_sz) * (pencil_sz) * (pencil_sz)];
        
        const DEVICE_INT Z_ID         = GPUDevice::block_id_x();
        const DEVICE_INT V_ID         = GPUDevice::block_id_y();
        
        if(Z_ID >=nblocks)
            return;
        
        const DEVICE_INT sz_per_dof   = szpdof_uz;
        const int pencils = pencil_sz;
        
        BSSN_EVAR_DERIVS * __restrict__ const deriv_evars = dptr_deriv_evars;
        BlockGPU3D c_blk = dptr_blk[Z_ID];

        const BlockGPU3D * __restrict__ const blk         = &c_blk;
        const DEVICE_UINT offset            = blk->m_offset; 
        const DEVICE_REAL* __restrict__ const var_in = &dptr_u[V_ID* sz_per_dof + offset];
        const DEVICE_UINT BLK_SZ            = nx*nx*nx; 

        #include "../scripts/compute_bssnrhs_evar_derivs2.cuh"

    }

    template<int pw, int nx, int pencil_sz>
    GLOBAL_FUNC
    __launch_bounds__(64)
    void eval_rhs2(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, const BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {
        const DEVICE_UINT sz_per_dof   = szpdof_uz;
        const DEVICE_UINT Z_ID         = GPUDevice::block_id_x();
        
        const DEVICE_UINT si           = GPUDevice::thread_id_x() + pw;
        const DEVICE_UINT sj           = GPUDevice::thread_id_y() + pw;

        if(Z_ID >=nblocks)
            return;
        
        BlockGPU3D c_blk = dptr_blk[Z_ID];
        const BlockGPU3D * __restrict__ const blk         = &c_blk;

        const DEVICE_UINT lambda[4]    = {1, 1, 1, 1};
        const DEVICE_REAL lambda_f[2]  = {1.0, 0.0};
        const DEVICE_UINT offset            = blk->m_offset; 
        const DEVICE_UINT actual_nx         = blk->m_sz[0];
        //const DEVICE_UINT nx              = blk[BLK_ID].m_aligned_sz[0];
        const DEVICE_UINT BLK_SZ            = nx*nx*nx; 

        DEVICE_REAL * __restrict__ const a_rhs    = &dptr_Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const chi_rhs  = &dptr_Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs00 = &dptr_Fu[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs01 = &dptr_Fu[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs02 = &dptr_Fu[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs11 = &dptr_Fu[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs12 = &dptr_Fu[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs22 = &dptr_Fu[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs0   = &dptr_Fu[bssn::VAR::U_BETA0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs1   = &dptr_Fu[bssn::VAR::U_BETA1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs2   = &dptr_Fu[bssn::VAR::U_BETA2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs00 = &dptr_Fu[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs01 = &dptr_Fu[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs02 = &dptr_Fu[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs11 = &dptr_Fu[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs12 = &dptr_Fu[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs22 = &dptr_Fu[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs0  = &dptr_Fu[bssn::VAR::U_GT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs1  = &dptr_Fu[bssn::VAR::U_GT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs2  = &dptr_Fu[bssn::VAR::U_GT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs0   = &dptr_Fu[bssn::VAR::U_B0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs1   = &dptr_Fu[bssn::VAR::U_B1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs2   = &dptr_Fu[bssn::VAR::U_B2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const K_rhs    = &dptr_Fu[bssn::VAR::U_K * sz_per_dof + offset];

        const DEVICE_REAL *  __restrict__ const  alpha = (dptr_u + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  chi   = (dptr_u + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  K     = (dptr_u + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt0   = (dptr_u + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt1   = (dptr_u + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt2   = (dptr_u + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt3   = (dptr_u + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt4   = (dptr_u + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt5   = (dptr_u + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta0 = (dptr_u + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta1 = (dptr_u + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta2 = (dptr_u + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At0   = (dptr_u + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At1   = (dptr_u + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At2   = (dptr_u + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At3   = (dptr_u + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At4   = (dptr_u + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At5   = (dptr_u + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt0   = (dptr_u + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt1   = (dptr_u + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt2   = (dptr_u + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B0    = (dptr_u + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B1    = (dptr_u + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B2    = (dptr_u + bssn::VAR::U_B2 * sz_per_dof + offset);

        const DEVICE_REAL hx = blk->m_dx[0];
        const DEVICE_REAL hy = blk->m_dx[1];

        const DEVICE_REAL x = blk->m_ptMin[0] + si*hx;
        const DEVICE_REAL y = blk->m_ptMin[1] + sj*hy;
        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;
                
        
        if(si < nx-pw && sj < nx-pw)
        for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
        {
            const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
            const DEVICE_REAL hz = blk->m_dx[2];
            const DEVICE_REAL z = blk->m_ptMin[2] + sk*hz;
            const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

            const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
            const DEVICE_REAL arg = - w*w*w*w;
            const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;

            
            
            const DEVICE_UINT d_pp = Z_ID * BLK_SZ + pp;
            #include "../scripts/bssneqs_eta_const_standard_gauge4.cpp"
        }

        GPUDevice::sync_threads();
        if(blk->m_bflag != 0)
        {
         
            #include "../scripts/bssnrhs_evar_derivs_zid.h"
            
            radiative_bc_blk<pw>(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk,0);
            radiative_bc_blk<pw>(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk, 0);
            radiative_bc_blk<pw>(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, 1.0, 0.0, blk, 0);

            radiative_bc_blk<pw>(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk, 0);

            radiative_bc_blk<pw>(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk, 0);

            radiative_bc_blk<pw>(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk, 0);

            radiative_bc_blk<pw>(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk, 0);

            radiative_bc_blk<pw>(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk, 0);
            radiative_bc_blk<pw>(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk, 0);
            radiative_bc_blk<pw>(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk, 0);
            radiative_bc_blk<pw>(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk, 0);

        }

        GPUDevice::sync_threads();
        const DEVICE_REAL sigma = ko_sigma;
        if(si < nx-pw && sj < nx-pw)
        for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
        {
            const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
            #include "../scripts/bssnrhs_evar_derivs_zid.h"
            a_rhs[pp]    += sigma * (kograd_0_alpha[pp] + kograd_1_alpha[pp] + kograd_2_alpha[pp]);
    
            b_rhs0[pp]   += sigma * (kograd_0_beta0[pp] + kograd_1_beta0[pp] + kograd_2_beta0[pp]);
            b_rhs1[pp]   += sigma * (kograd_0_beta1[pp] + kograd_1_beta1[pp] + kograd_2_beta1[pp]);
            b_rhs2[pp]   += sigma * (kograd_0_beta2[pp] + kograd_1_beta2[pp] + kograd_2_beta2[pp]);

            gt_rhs00[pp] += sigma * (kograd_0_gt0[pp] + kograd_1_gt0[pp] + kograd_2_gt0[pp]);
            gt_rhs01[pp] += sigma * (kograd_0_gt1[pp] + kograd_1_gt1[pp] + kograd_2_gt1[pp]);
            gt_rhs02[pp] += sigma * (kograd_0_gt2[pp] + kograd_1_gt2[pp] + kograd_2_gt2[pp]);
            gt_rhs11[pp] += sigma * (kograd_0_gt3[pp] + kograd_1_gt3[pp] + kograd_2_gt3[pp]);
            gt_rhs12[pp] += sigma * (kograd_0_gt4[pp] + kograd_1_gt4[pp] + kograd_2_gt4[pp]);
            gt_rhs22[pp] += sigma * (kograd_0_gt5[pp] + kograd_1_gt5[pp] + kograd_2_gt5[pp]);

            chi_rhs[pp]  += sigma * (kograd_0_chi[pp] + kograd_1_chi[pp] + kograd_2_chi[pp]);

            At_rhs00[pp] += sigma * (kograd_0_At0[pp] + kograd_1_At0[pp] + kograd_2_At0[pp]);
            At_rhs01[pp] += sigma * (kograd_0_At1[pp] + kograd_1_At1[pp] + kograd_2_At1[pp]);
            At_rhs02[pp] += sigma * (kograd_0_At2[pp] + kograd_1_At2[pp] + kograd_2_At2[pp]);
            At_rhs11[pp] += sigma * (kograd_0_At3[pp] + kograd_1_At3[pp] + kograd_2_At3[pp]);
            At_rhs12[pp] += sigma * (kograd_0_At4[pp] + kograd_1_At4[pp] + kograd_2_At4[pp]);
            At_rhs22[pp] += sigma * (kograd_0_At5[pp] + kograd_1_At5[pp] + kograd_2_At5[pp]);

            K_rhs[pp]    += sigma * (kograd_0_K[pp] + kograd_1_K[pp] + kograd_2_K[pp]);

            Gt_rhs0[pp]  += sigma * (kograd_0_Gt0[pp] + kograd_1_Gt0[pp] + kograd_2_Gt0[pp]);
            Gt_rhs1[pp]  += sigma * (kograd_0_Gt1[pp] + kograd_1_Gt1[pp] + kograd_2_Gt1[pp]);
            Gt_rhs2[pp]  += sigma * (kograd_0_Gt2[pp] + kograd_1_Gt2[pp] + kograd_2_Gt2[pp]);

            B_rhs0[pp]   += sigma * (kograd_0_B0[pp] + kograd_1_B0[pp] + kograd_2_B0[pp]);
            B_rhs1[pp]   += sigma * (kograd_0_B1[pp] + kograd_1_B1[pp] + kograd_2_B1[pp]);
            B_rhs2[pp]   += sigma * (kograd_0_B2[pp] + kograd_1_B2[pp] + kograd_2_B2[pp]);
        }

    }

    template<int pw, int nx, int pencil_sz>
    GLOBAL_FUNC
    __launch_bounds__(64)
    void eval_rhs2a(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, const BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {
        const DEVICE_UINT sz_per_dof   = szpdof_uz;
        const DEVICE_UINT Z_ID         = GPUDevice::block_id_x();
        const DEVICE_UINT V_ID         = GPUDevice::block_id_y();

        const DEVICE_UINT si           = GPUDevice::thread_id_x() + pw;
        const DEVICE_UINT sj           = GPUDevice::thread_id_y() + pw;

        if(Z_ID >=nblocks)
            return;
        
        BlockGPU3D c_blk = dptr_blk[Z_ID];
        const BlockGPU3D * __restrict__ const blk         = &c_blk;

        const DEVICE_UINT lambda[4]    = {1, 1, 1, 1};
        const DEVICE_REAL lambda_f[2]  = {1.0, 0.0};
        const DEVICE_UINT offset            = blk->m_offset; 
        const DEVICE_UINT actual_nx         = blk->m_sz[0];
        //const DEVICE_UINT nx              = blk[BLK_ID].m_aligned_sz[0];
        const DEVICE_UINT BLK_SZ            = nx*nx*nx; 

        DEVICE_REAL * __restrict__ const a_rhs    = &dptr_Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const chi_rhs  = &dptr_Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs00 = &dptr_Fu[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs01 = &dptr_Fu[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs02 = &dptr_Fu[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs11 = &dptr_Fu[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs12 = &dptr_Fu[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs22 = &dptr_Fu[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs0   = &dptr_Fu[bssn::VAR::U_BETA0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs1   = &dptr_Fu[bssn::VAR::U_BETA1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs2   = &dptr_Fu[bssn::VAR::U_BETA2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs00 = &dptr_Fu[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs01 = &dptr_Fu[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs02 = &dptr_Fu[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs11 = &dptr_Fu[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs12 = &dptr_Fu[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs22 = &dptr_Fu[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs0  = &dptr_Fu[bssn::VAR::U_GT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs1  = &dptr_Fu[bssn::VAR::U_GT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs2  = &dptr_Fu[bssn::VAR::U_GT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs0   = &dptr_Fu[bssn::VAR::U_B0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs1   = &dptr_Fu[bssn::VAR::U_B1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs2   = &dptr_Fu[bssn::VAR::U_B2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const K_rhs    = &dptr_Fu[bssn::VAR::U_K * sz_per_dof + offset];

        const DEVICE_REAL * __restrict__ alpha = (dptr_u  + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ chi   = (dptr_u  + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ K     = (dptr_u  + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ gt0   = (dptr_u  + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ gt1   = (dptr_u  + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ gt2   = (dptr_u  + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ gt3   = (dptr_u  + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ gt4   = (dptr_u  + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ gt5   = (dptr_u  + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ beta0 = (dptr_u  + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ beta1 = (dptr_u  + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ beta2 = (dptr_u  + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ At0   = (dptr_u  + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ At1   = (dptr_u  + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ At2   = (dptr_u  + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ At3   = (dptr_u  + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ At4   = (dptr_u  + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ At5   = (dptr_u  + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ Gt0   = (dptr_u  + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ Gt1   = (dptr_u  + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ Gt2   = (dptr_u  + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ B0    = (dptr_u  + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ B1    = (dptr_u  + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL * __restrict__ B2    = (dptr_u  + bssn::VAR::U_B2 * sz_per_dof + offset);

        const DEVICE_REAL hx = blk->m_dx[0];
        const DEVICE_REAL hy = blk->m_dx[1];
        const DEVICE_REAL hz = blk->m_dx[2];
        
        const DEVICE_REAL x = blk->m_ptMin[0] + si*hx;
        const DEVICE_REAL y = blk->m_ptMin[1] + sj*hy;

        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;
        
        if(V_ID == 0)
        {
            //printf("execute: vid %d",V_ID);
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                const DEVICE_REAL z = blk->m_ptMin[2] + sk*hz;
                const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

                const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
                const DEVICE_REAL arg = - w*w*w*w;
                const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
                const DEVICE_UINT d_pp = Z_ID * BLK_SZ + pp;

                #include "../scripts/block_id_0.cpp"

            }

            GPUDevice::sync_threads();
            if(blk->m_bflag != 0)
            {
                #include "../scripts/bssnrhs_evar_derivs_zid.h"
                radiative_bc_blk<pw>(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk,0);
                radiative_bc_blk<pw>(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk, 0);

                radiative_bc_blk<pw>(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk, 0);

            }

            GPUDevice::sync_threads();
            const DEVICE_REAL sigma = ko_sigma;
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                #include "../scripts/bssnrhs_evar_derivs_ko_zid.h"
                a_rhs[pp]    += sigma * (kograd_0_alpha[pp] + kograd_1_alpha[pp] + kograd_2_alpha[pp]);
        
                b_rhs0[pp]   += sigma * (kograd_0_beta0[pp] + kograd_1_beta0[pp] + kograd_2_beta0[pp]);
                b_rhs1[pp]   += sigma * (kograd_0_beta1[pp] + kograd_1_beta1[pp] + kograd_2_beta1[pp]);
                b_rhs2[pp]   += sigma * (kograd_0_beta2[pp] + kograd_1_beta2[pp] + kograd_2_beta2[pp]);

                gt_rhs00[pp] += sigma * (kograd_0_gt0[pp] + kograd_1_gt0[pp] + kograd_2_gt0[pp]);
                gt_rhs01[pp] += sigma * (kograd_0_gt1[pp] + kograd_1_gt1[pp] + kograd_2_gt1[pp]);
                gt_rhs02[pp] += sigma * (kograd_0_gt2[pp] + kograd_1_gt2[pp] + kograd_2_gt2[pp]);
                gt_rhs11[pp] += sigma * (kograd_0_gt3[pp] + kograd_1_gt3[pp] + kograd_2_gt3[pp]);
                gt_rhs12[pp] += sigma * (kograd_0_gt4[pp] + kograd_1_gt4[pp] + kograd_2_gt4[pp]);
                gt_rhs22[pp] += sigma * (kograd_0_gt5[pp] + kograd_1_gt5[pp] + kograd_2_gt5[pp]);

                chi_rhs[pp]  += sigma * (kograd_0_chi[pp] + kograd_1_chi[pp] + kograd_2_chi[pp]);
               
            }
            
            return;
            
        }
        
        if(V_ID==1)
        {
            //printf("execute: vid %d",V_ID);
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                const DEVICE_REAL z = blk->m_ptMin[2] + sk*hz;
                const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

                const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
                const DEVICE_REAL arg = - w*w*w*w;
                const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
                const DEVICE_UINT d_pp = Z_ID * BLK_SZ + pp;

                #include "../scripts/block_id_1.cpp"

            }

            GPUDevice::sync_threads();
            if(blk->m_bflag != 0)
            {
                #include "../scripts/bssnrhs_evar_derivs_zid.h"
                radiative_bc_blk<pw>(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk, 0);

            }

            GPUDevice::sync_threads();
            const DEVICE_REAL sigma = ko_sigma;
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                #include "../scripts/bssnrhs_evar_derivs_ko_zid.h"
                At_rhs00[pp] += sigma * (kograd_0_At0[pp] + kograd_1_At0[pp] + kograd_2_At0[pp]);
                At_rhs01[pp] += sigma * (kograd_0_At1[pp] + kograd_1_At1[pp] + kograd_2_At1[pp]);
                At_rhs02[pp] += sigma * (kograd_0_At2[pp] + kograd_1_At2[pp] + kograd_2_At2[pp]);
                At_rhs11[pp] += sigma * (kograd_0_At3[pp] + kograd_1_At3[pp] + kograd_2_At3[pp]);
                At_rhs12[pp] += sigma * (kograd_0_At4[pp] + kograd_1_At4[pp] + kograd_2_At4[pp]);
                At_rhs22[pp] += sigma * (kograd_0_At5[pp] + kograd_1_At5[pp] + kograd_2_At5[pp]);

            }
            
            return;
            
        }

        if(V_ID==2)
        {
            //printf("execute: vid %d",V_ID);
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                const DEVICE_REAL z = blk->m_ptMin[2] + sk*hz;
                const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

                const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
                const DEVICE_REAL arg = - w*w*w*w;
                const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
                const DEVICE_UINT d_pp = Z_ID * BLK_SZ + pp;

                #include "../scripts/block_id_2.cpp"

            }

            GPUDevice::sync_threads();
            if(blk->m_bflag != 0)
            {
                #include "../scripts/bssnrhs_evar_derivs_zid.h"
                radiative_bc_blk<pw>(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk, 0);

                radiative_bc_blk<pw>(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk, 0);
                radiative_bc_blk<pw>(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk, 0);

            }

            GPUDevice::sync_threads();
            const DEVICE_REAL sigma = ko_sigma;
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                #include "../scripts/bssnrhs_evar_derivs_ko_zid.h"
                
                Gt_rhs0[pp]  += sigma * (kograd_0_Gt0[pp] + kograd_1_Gt0[pp] + kograd_2_Gt0[pp]);
                Gt_rhs1[pp]  += sigma * (kograd_0_Gt1[pp] + kograd_1_Gt1[pp] + kograd_2_Gt1[pp]);
                Gt_rhs2[pp]  += sigma * (kograd_0_Gt2[pp] + kograd_1_Gt2[pp] + kograd_2_Gt2[pp]);

                B_rhs0[pp]   += sigma * (kograd_0_B0[pp] + kograd_1_B0[pp] + kograd_2_B0[pp]);
                B_rhs1[pp]   += sigma * (kograd_0_B1[pp] + kograd_1_B1[pp] + kograd_2_B1[pp]);
                B_rhs2[pp]   += sigma * (kograd_0_B2[pp] + kograd_1_B2[pp] + kograd_2_B2[pp]);

            }
            
            return;
            
        }
        
        if(V_ID==3)
        {
            //printf("execute: vid %d",V_ID);
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                const DEVICE_REAL z = blk->m_ptMin[2] + sk*hz;
                const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

                const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
                const DEVICE_REAL arg = - w*w*w*w;
                const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
                const DEVICE_UINT d_pp = Z_ID * BLK_SZ + pp;

                #include "../scripts/block_id_3.cpp"

            }

            GPUDevice::sync_threads();
            if(blk->m_bflag != 0)
            {
                #include "../scripts/bssnrhs_evar_derivs_zid.h"
                radiative_bc_blk<pw>(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, 1.0, 0.0, blk, 0);
                
            }

            GPUDevice::sync_threads();
            const DEVICE_REAL sigma = ko_sigma;
            if(si < nx-pw && sj < nx-pw)
            for(DEVICE_UINT   sk = pw; sk < actual_nx -pw; sk++)
            {
                const DEVICE_UINT pp = sk * nx * nx + sj * nx + si;
                #include "../scripts/bssnrhs_evar_derivs_ko_zid.h"
                
                K_rhs[pp]    += sigma * (kograd_0_K[pp] + kograd_1_K[pp] + kograd_2_K[pp]);
            }
            
            
            return;
            

        }

        return;
        
        
    }


    template<int pw, int pencil_sz, int nx, int BATCHED_GRAIN_SZ, int nstreams>
    HOST_FUNC
    void bssnrhs2(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const dptr_deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {
        dim3 gb = dim3(BATCHED_GRAIN_SZ, bssn::BSSN_NUM_VARS, 1);
        //dim3 tb = dim3(pencil_sz, pencil_sz, 1);
        dim3 tb = dim3(13, 13, 1);
        dim3 gb2  = dim3(BATCHED_GRAIN_SZ, 1, 1);
        //dim3 gb2a = dim3(BATCHED_GRAIN_SZ, 4, 1);
        dim3 tb2 = dim3(8, 8, 1);
        cudaStream_t stream[nstreams];
        for(unsigned int i =0; i < nstreams; i++)
            cudaStreamCreate(&stream[i]);

        const unsigned int num_batches = nblocks/BATCHED_GRAIN_SZ + 1;
        for(unsigned int batch_id =0; batch_id < num_batches; batch_id++)
        {
            const unsigned int blk_begin = (batch_id * nblocks)/num_batches;
            const unsigned int blk_end   = ((batch_id + 1) * nblocks)/num_batches;
            const unsigned int sid = batch_id % nstreams;
            const unsigned int nblock_p_batch = blk_end - blk_begin;
            compute_all_derivs<pw, nx, pencil_sz> <<<gb,  tb, 0, stream[sid] >>>(dptr_Fu, dptr_u, dptr_blk + blk_begin, nblock_p_batch, dptr_deriv_evars +sid, szpdof_uz, ko_sigma);
            eval_rhs2<pw,nx,pencil_sz><<< gb2,tb2, 0, stream[sid] >>>(dptr_Fu, dptr_u, dptr_blk + blk_begin, nblock_p_batch, dptr_deriv_evars +sid, szpdof_uz,ko_sigma);
            //eval_rhs2a<pw,nx,pencil_sz><<< gb2a,tb2, 0, stream[sid] >>>(dptr_Fu, dptr_u, dptr_blk + blk_begin, nblock_p_batch, dptr_deriv_evars +sid, szpdof_uz,ko_sigma);
        }

        for(unsigned int i =0; i < nstreams; i++)
            cudaStreamDestroy(stream[i]);

        cudaDeviceSynchronize();

        return;
    }


    template<int pw, int nx, int pencil_sz>
    GLOBAL_FUNC
    __launch_bounds__(343,3)
    void eval_rhs3(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, const BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {

        extern __shared__ DEVICE_REAL shared_mem[];
        // __shared__ DEVICE_REAL su[pencil_sz * pencil_sz *pencil_sz];
        // __shared__ DEVICE_REAL Du[pencil_sz * pencil_sz *pencil_sz];

        DEVICE_REAL *  su  = shared_mem;
        DEVICE_REAL *  Du  = shared_mem +     pencil_sz * pencil_sz *pencil_sz;
        DEVICE_REAL *  DDu = shared_mem + 2 * pencil_sz * pencil_sz *pencil_sz;

        const DEVICE_UINT sz_per_dof   = szpdof_uz;
        const DEVICE_UINT lambda[4]    = {1, 1, 1, 1};
        const DEVICE_REAL lambda_f[2]  = {1.0, 0.0};
        const DEVICE_INT Z_ID = GPUDevice::block_id_x();

        BlockGPU3D c_blk                          = dptr_blk[Z_ID];
        const BlockGPU3D * __restrict__ const blk = &c_blk;

        
        const DEVICE_UINT offset            = blk->m_offset; 
        const DEVICE_UINT actual_nx         = blk->m_sz[0];
        const DEVICE_UINT BLK_SZ            = nx*nx*nx; 
        const int pencils = pencil_sz;

        DEVICE_REAL * __restrict__ const a_rhs    = &dptr_Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const chi_rhs  = &dptr_Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs00 = &dptr_Fu[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs01 = &dptr_Fu[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs02 = &dptr_Fu[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs11 = &dptr_Fu[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs12 = &dptr_Fu[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs22 = &dptr_Fu[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs0   = &dptr_Fu[bssn::VAR::U_BETA0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs1   = &dptr_Fu[bssn::VAR::U_BETA1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs2   = &dptr_Fu[bssn::VAR::U_BETA2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs00 = &dptr_Fu[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs01 = &dptr_Fu[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs02 = &dptr_Fu[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs11 = &dptr_Fu[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs12 = &dptr_Fu[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs22 = &dptr_Fu[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs0  = &dptr_Fu[bssn::VAR::U_GT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs1  = &dptr_Fu[bssn::VAR::U_GT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs2  = &dptr_Fu[bssn::VAR::U_GT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs0   = &dptr_Fu[bssn::VAR::U_B0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs1   = &dptr_Fu[bssn::VAR::U_B1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs2   = &dptr_Fu[bssn::VAR::U_B2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const K_rhs    = &dptr_Fu[bssn::VAR::U_K * sz_per_dof + offset];

        const DEVICE_REAL *  __restrict__ const  alpha = (dptr_u + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  chi   = (dptr_u + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  K     = (dptr_u + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt0   = (dptr_u + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt1   = (dptr_u + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt2   = (dptr_u + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt3   = (dptr_u + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt4   = (dptr_u + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt5   = (dptr_u + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta0 = (dptr_u + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta1 = (dptr_u + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta2 = (dptr_u + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At0   = (dptr_u + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At1   = (dptr_u + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At2   = (dptr_u + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At3   = (dptr_u + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At4   = (dptr_u + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At5   = (dptr_u + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt0   = (dptr_u + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt1   = (dptr_u + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt2   = (dptr_u + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B0    = (dptr_u + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B1    = (dptr_u + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B2    = (dptr_u + bssn::VAR::U_B2 * sz_per_dof + offset);

        const DEVICE_INT i  = GPUDevice::thread_id_x() + pw;
        const DEVICE_INT j  = GPUDevice::thread_id_y() + pw;
        const DEVICE_INT k  = GPUDevice::thread_id_z() + pw;

        const int gidx = (k * nx + j) * nx + i ; 

        //#include "../scripts/compute_bssnrhs_evar_derivs3.cuh"
        const DEVICE_INT pp = gidx;

        const DEVICE_REAL hx = blk->m_dx[0];
        
        const DEVICE_REAL x = blk->m_ptMin[0] + i*hx;
        const DEVICE_REAL y = blk->m_ptMin[1] + j*hx;
        const DEVICE_REAL z = blk->m_ptMin[1] + k*hx;

        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;

        const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

        const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
        const DEVICE_REAL arg = - w*w*w*w;
        const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
        
        // #pragma message ("key_expr_sympy_cse_with_fused")
        #include "../scripts/bssneqs_eta_const_standard_gauge6.cpp"

        /*#include "../scripts/compute_bssnrhs_evar_derivs3.cuh"
        //#include "../scripts/bssneqs_eta_const_standard_gauge5.cpp"
        
        // #pragma message("sympy_cse")
        // #include "../../../CodeGen/bssneqs_sympy_cse_wo_derivs.cpp"
        
        // #pragma message("nx_cse_without_derivs")
        // #include "../../../CodeGen/bssneqs_nx_cse_wo_derivs.cpp"

        // #pragma message("key_exprs")
        // #include "../../../CodeGen/bssneqs_key_expr_sympy_cse_wo_derivs.cpp"



        GPUDevice::sync_threads();
        if(blk->m_bflag != 0)
        {
            radiative_bc_pt<pw , nx>(&a_rhs [gidx] , alpha[gidx], grad_0_alpha, grad_1_alpha, grad_2_alpha,1.0, 1.0,blk);

            radiative_bc_pt<pw , nx>(&b_rhs0[gidx] , beta0[gidx], grad_0_beta0, grad_1_beta0, grad_2_beta0, 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&b_rhs1[gidx] , beta1[gidx], grad_0_beta1, grad_1_beta1, grad_2_beta1, 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&b_rhs2[gidx] , beta2[gidx], grad_0_beta2, grad_1_beta2, grad_2_beta2, 1.0, 0.0, blk);
            
            radiative_bc_pt<pw , nx>(&chi_rhs[gidx], chi[gidx]  , grad_0_chi, grad_1_chi, grad_2_chi, 1.0, 1.0, blk);
            radiative_bc_pt<pw , nx>(&K_rhs[gidx],     K[gidx]  ,   grad_0_K,   grad_1_K,   grad_2_K, 1.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&Gt_rhs0[gidx], Gt0[gidx], grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&Gt_rhs1[gidx], Gt1[gidx], grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&Gt_rhs2[gidx], Gt2[gidx], grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, 2.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&B_rhs0[gidx], B0[gidx], grad_0_B0, grad_1_B0, grad_2_B0, 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&B_rhs1[gidx], B1[gidx], grad_0_B1, grad_1_B1, grad_2_B1, 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&B_rhs2[gidx], B2[gidx], grad_0_B2, grad_1_B2, grad_2_B2, 1.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&At_rhs00[gidx], At0[gidx], grad_0_At0, grad_1_At0, grad_2_At0, 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs01[gidx], At1[gidx], grad_0_At1, grad_1_At1, grad_2_At1, 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs02[gidx], At2[gidx], grad_0_At2, grad_1_At2, grad_2_At2, 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs11[gidx], At3[gidx], grad_0_At3, grad_1_At3, grad_2_At3, 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs12[gidx], At4[gidx], grad_0_At4, grad_1_At4, grad_2_At4, 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs22[gidx], At5[gidx], grad_0_At5, grad_1_At5, grad_2_At5, 2.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&gt_rhs00[gidx], gt0[gidx], grad_0_gt0, grad_1_gt0, grad_2_gt0, 1.0, 1.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs01[gidx], gt1[gidx], grad_0_gt1, grad_1_gt1, grad_2_gt1, 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs02[gidx], gt2[gidx], grad_0_gt2, grad_1_gt2, grad_2_gt2, 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs11[gidx], gt3[gidx], grad_0_gt3, grad_1_gt3, grad_2_gt3, 1.0, 1.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs12[gidx], gt4[gidx], grad_0_gt4, grad_1_gt4, grad_2_gt4, 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs22[gidx], gt5[gidx], grad_0_gt5, grad_1_gt5, grad_2_gt5, 1.0, 1.0, blk);
        }

        GPUDevice::sync_threads();
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
        B_rhs2[pp]   += sigma * (kograd_0_B2 + kograd_1_B2 + kograd_2_B2);*/

        return;


    }

    template<int pw, int pencil_sz, int nx, int BATCHED_GRAIN_SZ, int nstreams>
    HOST_FUNC
    void bssnrhs3(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {
        dim3 gb = dim3(nblocks,1, 1);
        dim3 tb = dim3(7, 7, 7);

        eval_rhs3<pw,nx,pencil_sz><<< gb,tb, 13*13*13 *8 * 3 >>>(dptr_Fu, dptr_u, dptr_blk, nblocks, szpdof_uz,ko_sigma);
        cudaDeviceSynchronize();

        return;
    }



    template<int pw, int nx, int pencil_sz>
    GLOBAL_FUNC
    __launch_bounds__(343,3)
    void eval_rhs4(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, const BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const dptr_deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {

        //extern __shared__ DEVICE_REAL shared_mem[];
        
        // DEVICE_REAL *  su  = shared_mem;
        // DEVICE_REAL *  Du  = shared_mem +     pencil_sz * pencil_sz *pencil_sz;
        // DEVICE_REAL *  DDu = shared_mem + 2 * pencil_sz * pencil_sz *pencil_sz;
        __shared__ DEVICE_REAL su[pencil_sz * pencil_sz *pencil_sz];
        __shared__ DEVICE_REAL Du[pencil_sz * pencil_sz *pencil_sz];

        const DEVICE_UINT sz_per_dof   = szpdof_uz;
        const DEVICE_UINT lambda[4]    = {1, 1, 1, 1};
        const DEVICE_REAL lambda_f[2]  = {1.0, 0.0};
        const DEVICE_INT Z_ID = GPUDevice::block_id_x();

        BlockGPU3D c_blk                           = dptr_blk[Z_ID];
        const BlockGPU3D *        __restrict__ const blk   = &c_blk;
        BSSN_EVAR_DERIVS * const  deriv_evars = dptr_deriv_evars; 

        
        const DEVICE_UINT offset            = blk->m_offset; 
        const DEVICE_UINT actual_nx         = blk->m_sz[0];
        const DEVICE_UINT BLK_SZ            = nx*nx*nx; 
        const int pencils = pencil_sz;

        DEVICE_REAL * __restrict__ const a_rhs    = &dptr_Fu[bssn::VAR::U_ALPHA * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const chi_rhs  = &dptr_Fu[bssn::VAR::U_CHI * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs00 = &dptr_Fu[bssn::VAR::U_SYMGT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs01 = &dptr_Fu[bssn::VAR::U_SYMGT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs02 = &dptr_Fu[bssn::VAR::U_SYMGT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs11 = &dptr_Fu[bssn::VAR::U_SYMGT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs12 = &dptr_Fu[bssn::VAR::U_SYMGT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const gt_rhs22 = &dptr_Fu[bssn::VAR::U_SYMGT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs0   = &dptr_Fu[bssn::VAR::U_BETA0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs1   = &dptr_Fu[bssn::VAR::U_BETA1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const b_rhs2   = &dptr_Fu[bssn::VAR::U_BETA2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs00 = &dptr_Fu[bssn::VAR::U_SYMAT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs01 = &dptr_Fu[bssn::VAR::U_SYMAT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs02 = &dptr_Fu[bssn::VAR::U_SYMAT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs11 = &dptr_Fu[bssn::VAR::U_SYMAT3 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs12 = &dptr_Fu[bssn::VAR::U_SYMAT4 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const At_rhs22 = &dptr_Fu[bssn::VAR::U_SYMAT5 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs0  = &dptr_Fu[bssn::VAR::U_GT0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs1  = &dptr_Fu[bssn::VAR::U_GT1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const Gt_rhs2  = &dptr_Fu[bssn::VAR::U_GT2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs0   = &dptr_Fu[bssn::VAR::U_B0 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs1   = &dptr_Fu[bssn::VAR::U_B1 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const B_rhs2   = &dptr_Fu[bssn::VAR::U_B2 * sz_per_dof + offset];
        DEVICE_REAL * __restrict__ const K_rhs    = &dptr_Fu[bssn::VAR::U_K * sz_per_dof + offset];

        const DEVICE_REAL *  __restrict__ const  alpha = (dptr_u + bssn::VAR::U_ALPHA * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  chi   = (dptr_u + bssn::VAR::U_CHI * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  K     = (dptr_u + bssn::VAR::U_K * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt0   = (dptr_u + bssn::VAR::U_SYMGT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt1   = (dptr_u + bssn::VAR::U_SYMGT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt2   = (dptr_u + bssn::VAR::U_SYMGT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt3   = (dptr_u + bssn::VAR::U_SYMGT3 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt4   = (dptr_u + bssn::VAR::U_SYMGT4 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  gt5   = (dptr_u + bssn::VAR::U_SYMGT5 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta0 = (dptr_u + bssn::VAR::U_BETA0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta1 = (dptr_u + bssn::VAR::U_BETA1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  beta2 = (dptr_u + bssn::VAR::U_BETA2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At0   = (dptr_u + bssn::VAR::U_SYMAT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At1   = (dptr_u + bssn::VAR::U_SYMAT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At2   = (dptr_u + bssn::VAR::U_SYMAT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At3   = (dptr_u + bssn::VAR::U_SYMAT3 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At4   = (dptr_u + bssn::VAR::U_SYMAT4 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  At5   = (dptr_u + bssn::VAR::U_SYMAT5 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt0   = (dptr_u + bssn::VAR::U_GT0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt1   = (dptr_u + bssn::VAR::U_GT1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  Gt2   = (dptr_u + bssn::VAR::U_GT2 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B0    = (dptr_u + bssn::VAR::U_B0 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B1    = (dptr_u + bssn::VAR::U_B1 * sz_per_dof + offset);
        const DEVICE_REAL *  __restrict__ const  B2    = (dptr_u + bssn::VAR::U_B2 * sz_per_dof + offset);

        const DEVICE_INT i  = GPUDevice::thread_id_x() + pw;
        const DEVICE_INT j  = GPUDevice::thread_id_y() + pw;
        const DEVICE_INT k  = GPUDevice::thread_id_z() + pw;

        const int gidx = (k * nx + j) * nx + i ; 

        const DEVICE_REAL hx = blk->m_dx[0];
        
        const DEVICE_REAL x = blk->m_ptMin[0] + i*hx;
        const DEVICE_REAL y = blk->m_ptMin[1] + j*hx;
        const DEVICE_REAL z = blk->m_ptMin[1] + k*hx;

        const DEVICE_REAL RIT_ETA_OUTER   = 0.25;
        const DEVICE_REAL RIT_ETA_CENTRAL = 2.0;
        const DEVICE_REAL RIT_ETA_WIDTH   = 40.0;

        const DEVICE_REAL r_coord = sqrt(x*x + y*y + z*z);

        const DEVICE_REAL w = r_coord / RIT_ETA_WIDTH;
        const DEVICE_REAL arg = - w*w*w*w;
        const DEVICE_REAL eta = (RIT_ETA_CENTRAL - RIT_ETA_OUTER)*exp(arg) + RIT_ETA_OUTER;
        
        // for(int v=0; v < 24; v++)
        // device::__ld_blk_var1__<DEVICE_REAL,pw,nx>   (su , &dptr_u[v*sz_per_dof + offset], blk);
        #include "../scripts/compute_bssnrhs_evar_derivs4.cuh"
        GPUDevice::sync_threads();
        // device::__blk1_deriv644_x<pw,pencils,pencil_sz>(Du , su, blk);
        // __syncthreads();
        // DEVICE_REAL grad_0_alpha=Du[gidx];
        // __syncthreads();
        // device::__blk1_deriv644_y<pw,pencils,pencil_sz>(Du , su, blk);
        // __syncthreads();
        // DEVICE_REAL grad_1_alpha=Du[gidx];
        // __syncthreads();
        // device::__blk1_deriv644_z<pw,pencils,pencil_sz>(Du , su, blk);
        // __syncthreads();
        // DEVICE_REAL grad_2_alpha=Du[gidx];
        const DEVICE_INT pp   = gidx;
        const DEVICE_INT d_pp = Z_ID * BLK_SZ +  gidx;
        #include "../scripts/bssneqs_eta_const_standard_gauge4.cpp"

        #include "../scripts/bssnrhs_evar_derivs_zid.h"
        GPUDevice::sync_threads();
        if(blk->m_bflag != 0)
        {
            
            radiative_bc_pt<pw , nx>(&a_rhs [gidx] , alpha[gidx], grad_0_alpha[gidx], grad_1_alpha[gidx], grad_2_alpha[gidx],1.0, 1.0,blk);

            radiative_bc_pt<pw , nx>(&b_rhs0[gidx] , beta0[gidx], grad_0_beta0[gidx], grad_1_beta0[gidx], grad_2_beta0[gidx], 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&b_rhs1[gidx] , beta1[gidx], grad_0_beta1[gidx], grad_1_beta1[gidx], grad_2_beta1[gidx], 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&b_rhs2[gidx] , beta2[gidx], grad_0_beta2[gidx], grad_1_beta2[gidx], grad_2_beta2[gidx], 1.0, 0.0, blk);
            
            radiative_bc_pt<pw , nx>(&chi_rhs[gidx], chi[gidx]  , grad_0_chi[gidx], grad_1_chi[gidx], grad_2_chi[gidx], 1.0, 1.0, blk);
            radiative_bc_pt<pw , nx>(&K_rhs[gidx],     K[gidx]  ,   grad_0_K[gidx],   grad_1_K[gidx],   grad_2_K[gidx], 1.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&Gt_rhs0[gidx], Gt0[gidx], grad_0_Gt0[gidx], grad_1_Gt0[gidx], grad_2_Gt0[gidx], 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&Gt_rhs1[gidx], Gt1[gidx], grad_0_Gt1[gidx], grad_1_Gt1[gidx], grad_2_Gt1[gidx], 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&Gt_rhs2[gidx], Gt2[gidx], grad_0_Gt2[gidx], grad_1_Gt2[gidx], grad_2_Gt2[gidx], 2.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&B_rhs0[gidx], B0[gidx], grad_0_B0[gidx], grad_1_B0[gidx], grad_2_B0[gidx], 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&B_rhs1[gidx], B1[gidx], grad_0_B1[gidx], grad_1_B1[gidx], grad_2_B1[gidx], 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&B_rhs2[gidx], B2[gidx], grad_0_B2[gidx], grad_1_B2[gidx], grad_2_B2[gidx], 1.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&At_rhs00[gidx], At0[gidx], grad_0_At0[gidx], grad_1_At0[gidx], grad_2_At0[gidx], 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs01[gidx], At1[gidx], grad_0_At1[gidx], grad_1_At1[gidx], grad_2_At1[gidx], 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs02[gidx], At2[gidx], grad_0_At2[gidx], grad_1_At2[gidx], grad_2_At2[gidx], 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs11[gidx], At3[gidx], grad_0_At3[gidx], grad_1_At3[gidx], grad_2_At3[gidx], 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs12[gidx], At4[gidx], grad_0_At4[gidx], grad_1_At4[gidx], grad_2_At4[gidx], 2.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&At_rhs22[gidx], At5[gidx], grad_0_At5[gidx], grad_1_At5[gidx], grad_2_At5[gidx], 2.0, 0.0, blk);

            radiative_bc_pt<pw , nx>(&gt_rhs00[gidx], gt0[gidx], grad_0_gt0[gidx], grad_1_gt0[gidx], grad_2_gt0[gidx], 1.0, 1.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs01[gidx], gt1[gidx], grad_0_gt1[gidx], grad_1_gt1[gidx], grad_2_gt1[gidx], 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs02[gidx], gt2[gidx], grad_0_gt2[gidx], grad_1_gt2[gidx], grad_2_gt2[gidx], 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs11[gidx], gt3[gidx], grad_0_gt3[gidx], grad_1_gt3[gidx], grad_2_gt3[gidx], 1.0, 1.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs12[gidx], gt4[gidx], grad_0_gt4[gidx], grad_1_gt4[gidx], grad_2_gt4[gidx], 1.0, 0.0, blk);
            radiative_bc_pt<pw , nx>(&gt_rhs22[gidx], gt5[gidx], grad_0_gt5[gidx], grad_1_gt5[gidx], grad_2_gt5[gidx], 1.0, 1.0, blk);
        }

        GPUDevice::sync_threads();
        const DEVICE_REAL sigma = ko_sigma;
        a_rhs[pp]    += sigma * (kograd_0_alpha[gidx] + kograd_1_alpha[gidx] + kograd_2_alpha[gidx]);
        b_rhs0[pp]   += sigma * (kograd_0_beta0[gidx] + kograd_1_beta0[gidx] + kograd_2_beta0[gidx]);
        b_rhs1[pp]   += sigma * (kograd_0_beta1[gidx] + kograd_1_beta1[gidx] + kograd_2_beta1[gidx]);
        b_rhs2[pp]   += sigma * (kograd_0_beta2[gidx] + kograd_1_beta2[gidx] + kograd_2_beta2[gidx]);

        gt_rhs00[pp] += sigma * (kograd_0_gt0[gidx] + kograd_1_gt0[gidx] + kograd_2_gt0[gidx]);
        gt_rhs01[pp] += sigma * (kograd_0_gt1[gidx] + kograd_1_gt1[gidx] + kograd_2_gt1[gidx]);
        gt_rhs02[pp] += sigma * (kograd_0_gt2[gidx] + kograd_1_gt2[gidx] + kograd_2_gt2[gidx]);
        gt_rhs11[pp] += sigma * (kograd_0_gt3[gidx] + kograd_1_gt3[gidx] + kograd_2_gt3[gidx]);
        gt_rhs12[pp] += sigma * (kograd_0_gt4[gidx] + kograd_1_gt4[gidx] + kograd_2_gt4[gidx]);
        gt_rhs22[pp] += sigma * (kograd_0_gt5[gidx] + kograd_1_gt5[gidx] + kograd_2_gt5[gidx]);

        chi_rhs[pp]  += sigma * (kograd_0_chi[gidx] + kograd_1_chi[gidx] + kograd_2_chi[gidx]);

        At_rhs00[pp] += sigma * (kograd_0_At0[gidx] + kograd_1_At0[gidx] + kograd_2_At0[gidx]);
        At_rhs01[pp] += sigma * (kograd_0_At1[gidx] + kograd_1_At1[gidx] + kograd_2_At1[gidx]);
        At_rhs02[pp] += sigma * (kograd_0_At2[gidx] + kograd_1_At2[gidx] + kograd_2_At2[gidx]);
        At_rhs11[pp] += sigma * (kograd_0_At3[gidx] + kograd_1_At3[gidx] + kograd_2_At3[gidx]);
        At_rhs12[pp] += sigma * (kograd_0_At4[gidx] + kograd_1_At4[gidx] + kograd_2_At4[gidx]);
        At_rhs22[pp] += sigma * (kograd_0_At5[gidx] + kograd_1_At5[gidx] + kograd_2_At5[gidx]);

        K_rhs[pp]    += sigma * (kograd_0_K[gidx] + kograd_1_K[gidx] + kograd_2_K[gidx]);

        Gt_rhs0[pp]  += sigma * (kograd_0_Gt0[gidx] + kograd_1_Gt0[gidx] + kograd_2_Gt0[gidx]);
        Gt_rhs1[pp]  += sigma * (kograd_0_Gt1[gidx] + kograd_1_Gt1[gidx] + kograd_2_Gt1[gidx]);
        Gt_rhs2[pp]  += sigma * (kograd_0_Gt2[gidx] + kograd_1_Gt2[gidx] + kograd_2_Gt2[gidx]);

        B_rhs0[pp]   += sigma * (kograd_0_B0[gidx] + kograd_1_B0[gidx] + kograd_2_B0[gidx]);
        B_rhs1[pp]   += sigma * (kograd_0_B1[gidx] + kograd_1_B1[gidx] + kograd_2_B1[gidx]);
        B_rhs2[pp]   += sigma * (kograd_0_B2[gidx] + kograd_1_B2[gidx] + kograd_2_B2[gidx]);
        
        return;


    }

    template<int pw, int pencil_sz, int nx, int BATCHED_GRAIN_SZ, int nstreams>
    HOST_FUNC
    void bssnrhs4(DEVICE_REAL * const dptr_Fu, const DEVICE_REAL * const dptr_u, BlockGPU3D* dptr_blk, DEVICE_UINT nblocks, BSSN_EVAR_DERIVS* const dptr_deriv_evars,DEVICE_UINT szpdof_uz, DEVICE_REAL ko_sigma)
    {
        dim3 gb = dim3(nblocks,1, 1);
        dim3 tb = dim3(7, 7, 7);

        eval_rhs4<pw,nx,pencil_sz><<< gb,tb, 0 >>>(dptr_Fu, dptr_u, dptr_blk, nblocks,dptr_deriv_evars,szpdof_uz,ko_sigma);
        cudaDeviceSynchronize();

        return;
    }



    
}