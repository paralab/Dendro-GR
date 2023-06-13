/**
 * @file bc_cu.hpp
 * @brief Apply the boundary conditions in CUDA
 * @version 0.1
 * @date 2022-02-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include "device.h"
#include "block_gpu.h"
#include<math.h>

namespace device
{

    template<unsigned int pw>
    DEVICE_FUNC void radiative_bc(DEVICE_REAL * f_rhs, const DEVICE_REAL * f,  const DEVICE_REAL * dxf, const DEVICE_REAL * dyf, const DEVICE_REAL *  dzf, DEVICE_REAL f_falloff, DEVICE_REAL f_asymptotic, const BlockGPU3D* blk, DEVICE_UINT blk_id)
    {
        const DEVICE_UINT BLK_ID = blk_id;
        
        const DEVICE_REAL hx     = blk[BLK_ID].m_dx[0];
        const DEVICE_REAL hy     = blk[BLK_ID].m_dx[1];
        const DEVICE_REAL hz     = blk[BLK_ID].m_dx[2];

        const DEVICE_UINT nx     = blk[BLK_ID].m_aligned_sz[0];
        const DEVICE_UINT ny     = blk[BLK_ID].m_aligned_sz[1];
        const DEVICE_UINT nz     = blk[BLK_ID].m_aligned_sz[2];

        const DEVICE_UINT actual_nx  = blk[BLK_ID].m_sz[0];
        const DEVICE_UINT actual_ny  = blk[BLK_ID].m_sz[1];
        const DEVICE_UINT actual_nz  = blk[BLK_ID].m_sz[2];

        const DEVICE_UINT bflag = blk[BLK_ID].m_bflag;
        
        {   
            // x - normal direction. 
            const DEVICE_UINT i      = GPUDevice::block_id_y();
            const DEVICE_UINT j      = GPUDevice::thread_id_x();
            const DEVICE_UINT k      = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
            const DEVICE_UINT globalIdx = k* ny * nx + j * nx + i;

            const DEVICE_REAL x      = blk[BLK_ID].m_ptMin[0] + i*hx;
            const DEVICE_REAL y      = blk[BLK_ID].m_ptMin[1] + j*hy;
            const DEVICE_REAL z      = blk[BLK_ID].m_ptMin[2] + k*hz;
            const DEVICE_REAL inv_r = 1.0 / sqrt(x*x + y*y + z*z);


            if (bflag & (1u<<OCT_DIR_LEFT)) {
                if(i == pw)
                {
                    f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                    + y * dyf[globalIdx]
                                                    + z * dzf[globalIdx]
                                                    + f_falloff * ( f[globalIdx] - f_asymptotic ) );
                }
            }

            if (bflag & (1u<<OCT_DIR_RIGHT)) {
                const DEVICE_UINT ie = actual_nx -pw;
                if(i == (ie-1))
                {
                    f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                    + y * dyf[globalIdx]
                                                    + z * dzf[globalIdx]
                                                    + f_falloff * ( f[globalIdx] - f_asymptotic ) );
                }
            }
        }

        
        {   // y -normal
            const DEVICE_UINT i      = GPUDevice::thread_id_x();
            const DEVICE_UINT j      = GPUDevice::block_id_y();
            const DEVICE_UINT k      = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();

            const DEVICE_UINT globalIdx = k* ny * nx + j * nx + i;

            const DEVICE_REAL x      = blk[BLK_ID].m_ptMin[0] + i*hx;
            const DEVICE_REAL y      = blk[BLK_ID].m_ptMin[1] + j*hy;
            const DEVICE_REAL z      = blk[BLK_ID].m_ptMin[2] + k*hz;
            const DEVICE_REAL inv_r = 1.0 / sqrt(x*x + y*y + z*z);
            
            if (bflag & (1u<<OCT_DIR_DOWN)) {
                if(j == pw)
                {
                    f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                    + y * dyf[globalIdx]
                                                    + z * dzf[globalIdx]
                                                    + f_falloff * ( f[globalIdx] - f_asymptotic ) );
                }
            }

            if (bflag & (1u<<OCT_DIR_UP)) {
                const DEVICE_UINT je = actual_ny -pw;
                if(j == (je-1))
                {
                    f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                    + y * dyf[globalIdx]
                                                    + z * dzf[globalIdx]
                                                    + f_falloff * ( f[globalIdx] - f_asymptotic ) );
                }
            }

        }

        
        {
            // z -normal
            const DEVICE_UINT i      = GPUDevice::thread_id_x();
            const DEVICE_UINT j      = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
            const DEVICE_UINT k      = GPUDevice::block_id_y();

            const DEVICE_UINT globalIdx = k* ny * nx + j * nx + i;
            
            const DEVICE_REAL x      = blk[BLK_ID].m_ptMin[0] + i*hx;
            const DEVICE_REAL y      = blk[BLK_ID].m_ptMin[1] + j*hy;
            const DEVICE_REAL z      = blk[BLK_ID].m_ptMin[2] + k*hz;
            const DEVICE_REAL inv_r = 1.0 / sqrt(x*x + y*y + z*z);

            
            
            if (bflag & (1u<<OCT_DIR_BACK)) {
                if(k == pw)
                {
                    f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                    + y * dyf[globalIdx]
                                                    + z * dzf[globalIdx]
                                                    + f_falloff * ( f[globalIdx] - f_asymptotic ) );
                }
            }

            if (bflag & (1u<<OCT_DIR_FRONT)) {
                const DEVICE_UINT ke = actual_nz -pw;
                if(k == (ke-1))
                {
                    f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                    + y * dyf[globalIdx]
                                                    + z * dzf[globalIdx]
                                                    + f_falloff * ( f[globalIdx] - f_asymptotic ) );
                }
            }

        }
        
        return;
    }

    template<unsigned int pw>
    DEVICE_FUNC void radiative_bc_blk(DEVICE_REAL * f_rhs, const DEVICE_REAL * f,  const DEVICE_REAL * dxf, const DEVICE_REAL * dyf, const DEVICE_REAL *  dzf, DEVICE_REAL f_falloff, DEVICE_REAL f_asymptotic, const BlockGPU3D* blk, DEVICE_UINT blk_id)
    {
        const DEVICE_UINT BLK_ID = blk_id;
        
        const DEVICE_REAL hx     = blk[BLK_ID].m_dx[0];
        const DEVICE_REAL hy     = blk[BLK_ID].m_dx[1];
        const DEVICE_REAL hz     = blk[BLK_ID].m_dx[2];

        const DEVICE_REAL min_x  = blk[BLK_ID].m_ptMin[0];
        const DEVICE_REAL min_y  = blk[BLK_ID].m_ptMin[1];
        const DEVICE_REAL min_z  = blk[BLK_ID].m_ptMin[2];

        const DEVICE_UINT nx     = blk[BLK_ID].m_aligned_sz[0];
        const DEVICE_UINT ny     = blk[BLK_ID].m_aligned_sz[1];
        const DEVICE_UINT nz     = blk[BLK_ID].m_aligned_sz[2];

        const DEVICE_UINT actual_nx  = blk[BLK_ID].m_sz[0];
        const DEVICE_UINT actual_ny  = blk[BLK_ID].m_sz[1];
        const DEVICE_UINT actual_nz  = blk[BLK_ID].m_sz[2];

        const DEVICE_UINT bflag = blk[BLK_ID].m_bflag;

        const DEVICE_UINT si      = GPUDevice::thread_id_x() + pw;
        const DEVICE_UINT sj      = GPUDevice::thread_id_y() + pw;
        const DEVICE_REAL x       = min_x + si*hx;
        const DEVICE_REAL y       = min_y + sj*hy;
            
        if(!(si < (actual_nx -pw) && sj < (actual_ny -pw)))
            return;

        if (bflag & (1u<<OCT_DIR_LEFT)) {
            for(DEVICE_UINT sk=pw; sk < actual_nz-pw; sk++)
            if(si == pw)
            {
                const DEVICE_UINT globalIdx = sk* ny * nx + sj * nx + si;
                const DEVICE_REAL z         = min_z + sk*hz;
                const DEVICE_REAL inv_r     = 1.0 / sqrt(x*x + y*y + z*z);
                f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                + y * dyf[globalIdx]
                                                + z * dzf[globalIdx]
                                                + f_falloff * ( f[globalIdx] - f_asymptotic ) );
            }
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            const DEVICE_UINT ie = actual_nx -pw;
            for(DEVICE_UINT sk=pw; sk < actual_nz-pw; sk++)
            if(si == (ie-1))
            {
                const DEVICE_UINT globalIdx = sk* ny * nx + sj * nx + si;
                const DEVICE_REAL z         = min_z + sk*hz;
                const DEVICE_REAL inv_r     = 1.0 / sqrt(x*x + y*y + z*z);
                f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                + y * dyf[globalIdx]
                                                + z * dzf[globalIdx]
                                                + f_falloff * ( f[globalIdx] - f_asymptotic ) );
            }
        }

        if (bflag & (1u<<OCT_DIR_DOWN)) {
            for(DEVICE_UINT sk=pw; sk < actual_nz-pw; sk++)
            if(sj == pw)
            {
                const DEVICE_UINT globalIdx = sk* ny * nx + sj * nx + si;
                const DEVICE_REAL z         = min_z + sk*hz;
                const DEVICE_REAL inv_r     = 1.0 / sqrt(x*x + y*y + z*z);

                f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                + y * dyf[globalIdx]
                                                + z * dzf[globalIdx]
                                                + f_falloff * ( f[globalIdx] - f_asymptotic ) );
            }
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            const DEVICE_UINT je = actual_ny -pw;
            for(DEVICE_UINT sk=pw; sk < actual_nz-pw; sk++)
            if(sj == (je-1))
            {
                const DEVICE_UINT globalIdx = sk* ny * nx + sj * nx + si;
                const DEVICE_REAL z         = min_z + sk*hz;
                const DEVICE_REAL inv_r     = 1.0 / sqrt(x*x + y*y + z*z);
                
                f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                + y * dyf[globalIdx]
                                                + z * dzf[globalIdx]
                                                + f_falloff * ( f[globalIdx] - f_asymptotic ) );
            }
        }


        if (bflag & (1u<<OCT_DIR_BACK)) {
            DEVICE_UINT sk =pw;
            if(sk == pw)
            {
                const DEVICE_UINT globalIdx = sk* ny * nx + sj * nx + si;
                const DEVICE_REAL z         = min_z + sk*hz;
                const DEVICE_REAL inv_r     = 1.0 / sqrt(x*x + y*y + z*z);

                f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                + y * dyf[globalIdx]
                                                + z * dzf[globalIdx]
                                                + f_falloff * ( f[globalIdx] - f_asymptotic ) );
            }
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            const DEVICE_UINT ke = actual_nz -pw;
            DEVICE_UINT sk = ke-1;
            if(sk == (ke-1))
            {
                const DEVICE_UINT globalIdx = sk* ny * nx + sj * nx + si;
                const DEVICE_REAL z         = min_z + sk*hz;
                const DEVICE_REAL inv_r     = 1.0 / sqrt(x*x + y*y + z*z);

                f_rhs[globalIdx] = - inv_r * (    x * dxf[globalIdx]
                                                + y * dyf[globalIdx]
                                                + z * dzf[globalIdx]
                                                + f_falloff * ( f[globalIdx] - f_asymptotic ) );
            }
        }

        return;
    }


    template<int pw, int nx>
    DEVICE_FUNC void radiative_bc_pt(DEVICE_REAL *  __restrict__ f_rhs, DEVICE_REAL f,  DEVICE_REAL dxf, DEVICE_REAL dyf, DEVICE_REAL dzf, DEVICE_REAL f_falloff, DEVICE_REAL f_asymptotic, const BlockGPU3D* __restrict__ const blk)
    {
        const DEVICE_REAL hx        = blk->m_dx[0];
        const DEVICE_REAL min_x     = blk->m_ptMin[0];
        const DEVICE_REAL min_y     = blk->m_ptMin[1];
        const DEVICE_REAL min_z     = blk->m_ptMin[2];
        const DEVICE_UINT bflag     = blk->m_bflag;
        const DEVICE_UINT actual_nx = blk->m_sz[0];
        const DEVICE_UINT actual_ny = blk->m_sz[1];
        const DEVICE_UINT actual_nz = blk->m_sz[2];

        const DEVICE_UINT si        = GPUDevice::thread_id_x() + pw;
        const DEVICE_UINT sj        = GPUDevice::thread_id_y() + pw;
        const DEVICE_UINT sk        = GPUDevice::thread_id_z() + pw;

        const DEVICE_REAL x         = min_x + si*hx;
        const DEVICE_REAL y         = min_y + sj*hx;
        const DEVICE_REAL z         = min_z + sk*hx;

        //const DEVICE_UINT globalIdx = 0;//(sk * nx + sj) * nx + si;
        const DEVICE_REAL inv_r     = 1.0 / sqrt(x*x + y*y + z*z);

        if (bflag & (1u<<OCT_DIR_LEFT)) {
            
            if(si == pw)
            {
                *f_rhs           = - inv_r * (    x * dxf
                                                + y * dyf
                                                + z * dzf
                                                + f_falloff * (f - f_asymptotic ));
            }
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            const DEVICE_UINT ie = actual_nx -pw;
            
            if(si == (ie-1))
            {
                *f_rhs           = - inv_r * (    x * dxf
                                                + y * dyf
                                                + z * dzf
                                                + f_falloff * ( f - f_asymptotic ) );
            }
        }

        if (bflag & (1u<<OCT_DIR_DOWN)) {
            
            if(sj == pw)
            {
                *f_rhs           = - inv_r * (    x * dxf
                                                + y * dyf
                                                + z * dzf
                                                + f_falloff * ( f - f_asymptotic ) );
            }
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            const DEVICE_UINT je = actual_ny -pw;
            
            if(sj == (je-1))
            {
                *f_rhs           = - inv_r * (    x * dxf
                                                + y * dyf
                                                + z * dzf
                                                + f_falloff * ( f - f_asymptotic ) );
            }
        }


        if (bflag & (1u<<OCT_DIR_BACK)) {
            
            if(sk == pw)
            {
                *f_rhs           = - inv_r * (    x * dxf
                                                + y * dyf
                                                + z * dzf
                                                + f_falloff * ( f - f_asymptotic ) );
            }
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            const DEVICE_UINT ke = actual_nz -pw;
            if(sk == (ke-1))
            {
                *f_rhs           = - inv_r * (    x * dxf
                                                + y * dyf
                                                + z * dzf
                                                + f_falloff * ( f - f_asymptotic ) );
            }
        }

        return;
    }

    
   
}