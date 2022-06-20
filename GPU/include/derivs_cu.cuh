/**
 * @file derivs_cu.h
 * @brief finite difference stencil approximation with CUDA. 
 * @version 0.1
 * @date 2021-10-08
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once
#include <stdio.h>
#include <iostream>
#include "device.h"
#include "block_gpu.h"

namespace device
{   
    /**
     * @brief Apply the 6th order 1st deriv x-deriv
     * @param Du : out - applied stencil
     * @param u   : in - input vector with padding
     * @param dx  : dx for the block
     * @param sz  : sz of the block (x,y,z) dimension
     * @param bflag : boundary flag. 
     * @return void 
     * */
    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __deriv644_x(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::block_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();  
        const DEVICE_INT sj = GPUDevice::thread_id_y();
        
        
        const DEVICE_REAL idx = 1.0 / blk->m_dx[0];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        // x derivative
        SHARED_MEM DEVICE_REAL s_u[pencils][pencil_sz + 2*pw];
        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sj][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = si + pw;
        globalIdx = k * nx * ny + j * nx + (sb % nx);
        Du[globalIdx] = (  -         s_u[sj][sb -3] 
                            +  9.0 * s_u[sj][sb -2]
                            - 45.0 * s_u[sj][sb -1]
                            + 45.0 * s_u[sj][sb +1]
                            -  9.0 * s_u[sj][sb +2]
                            +        s_u[sj][sb +3]) * idx_by_60;

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_LEFT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            if(sb == pw)
            {
                Du[globalIdx]      = (  - 25.0 * s_u[sj][sb  ]
                                        + 48.0 * s_u[sj][sb+1]
                                        - 36.0 * s_u[sj][sb+2]
                                        + 16.0 * s_u[sj][sb+3]
                                        -  3.0 * s_u[sj][sb+4]
                                      ) * idx_by_12;

                Du[globalIdx + 1]  = (  -  3.0 * s_u[sj][sb  ]
                                        - 10.0 * s_u[sj][sb+1]
                                        + 18.0 * s_u[sj][sb+2]
                                        -  6.0 * s_u[sj][sb+3]
                                        +        s_u[sj][sb+4]
                                      ) * idx_by_12;

                Du[globalIdx + 2]  = (   +      s_u[sj][sb  ]
                                        - 8.0 * s_u[sj][sb+1]
                                        + 8.0 * s_u[sj][sb+3]
                                        -       s_u[sj][sb+4]
                                      ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_UINT ie = actual_nx-pw;
            if(sb == ie-1)
            {
                Du[globalIdx-2]   = (  +        s_u[sj][sb-4]
                                        - 8.0 * s_u[sj][sb-3]
                                        + 8.0 * s_u[sj][sb-1]
                                        -       s_u[sj][sb  ]
                                        ) * idx_by_12;

                Du[globalIdx-1] =  (  -          s_u[sj][sb-4]
                                        +  6.0 * s_u[sj][sb-3]
                                        - 18.0 * s_u[sj][sb-2]
                                        + 10.0 * s_u[sj][sb-1]
                                        +  3.0 * s_u[sj][sb  ]
                                        ) * idx_by_12;
                                    
                Du[globalIdx] =    (       3.0 * s_u[sj][sb-4]
                                        - 16.0 * s_u[sj][sb-3]
                                        + 36.0 * s_u[sj][sb-2]
                                        - 48.0 * s_u[sj][sb-1]
                                        + 25.0 * s_u[sj][sb  ]
                                        ) * idx_by_12;
            }

        }
        
        return ;
    }

    /**
     * @brief Apply the 6th order 2nd deriv x-deriv
     * @param DDu : out - applied stencil
     * @param u   : in - input vector with padding
     * @param dx  : dx for the block
     * @param sz  : sz of the block (x,y,z) dimension
     * @param bflag : boundary flag. 
     * @return void 
     * */
    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __deriv644_xx(DEVICE_REAL * const  DDu, const DEVICE_REAL * const  u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::block_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();  
        const DEVICE_INT sj = GPUDevice::thread_id_y();
        
        const DEVICE_REAL idx_sqrd = 1.0 / (blk->m_dx[0]*blk->m_dx[0]);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;
        
        const DEVICE_INT nx        = blk->m_aligned_sz[0];
        const DEVICE_INT ny        = blk->m_aligned_sz[1];
        const DEVICE_INT nz        = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        // x derivative
        SHARED_MEM DEVICE_REAL s_u[pencils][pencil_sz + 2*pw];
        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sj][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = si + pw;
        globalIdx = k * nx * ny + j * nx + (sb % nx);
        DDu[globalIdx] = (      2.0 * s_u[sj][sb-3]
                            -  27.0 * s_u[sj][sb-2]
                            + 270.0 * s_u[sj][sb-1]
                            - 490.0 * s_u[sj][sb  ]
                            + 270.0 * s_u[sj][sb+1]
                            -  27.0 * s_u[sj][sb+2]
                            +   2.0 * s_u[sj][sb+3]
                            ) * idx_sqrd_by_180;

        const DEVICE_UINT bflag = blk->m_bflag;
    
        if (bflag & (1u<<OCT_DIR_LEFT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12  = idx_sqrd / 12.0;
            if(sb == pw )
            {
                DDu[globalIdx]      = (       45.0  * s_u[sj][sb    ]
                                            - 154.0 * s_u[sj][sb + 1]
                                            + 214.0 * s_u[sj][sb + 2]
                                            - 156.0 * s_u[sj][sb + 3]
                                            +  61.0 * s_u[sj][sb + 4]
                                            -  10.0 * s_u[sj][sb + 5]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx+1]    = (        10.0 * s_u[sj][sb    ]
                                            -  15.0 * s_u[sj][sb + 1]
                                            -   4.0 * s_u[sj][sb + 2]
                                            +  14.0 * s_u[sj][sb + 3]
                                            -   6.0 * s_u[sj][sb + 4]
                                            +         s_u[sj][sb + 5]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx+2]    = (     -         s_u[sj][sb    ]
                                            +  16.0 * s_u[sj][sb + 1]
                                            -  30.0 * s_u[sj][sb + 2]
                                            +  16.0 * s_u[sj][sb + 3]
                                            -         s_u[sj][sb + 4]
                                            ) * idx_sqrd_by_12;
            }
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12  = idx_sqrd / 12.0;
            const DEVICE_UINT ie              = actual_nx - pw;

            if(sb == ie-1)
            {
                DDu[globalIdx-2]    = (     -         s_u[sj][sb -4]
                                            +  16.0 * s_u[sj][sb -3]
                                            -  30.0 * s_u[sj][sb -2]
                                            +  16.0 * s_u[sj][sb -1]
                                            -         s_u[sj][sb   ]
                                            ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx-1]    = (               s_u[sj][sb-5]
                                            -   6.0 * s_u[sj][sb-4]
                                            +  14.0 * s_u[sj][sb-3]
                                            -   4.0 * s_u[sj][sb-2]
                                            -  15.0 * s_u[sj][sb-1]
                                            +  10.0 * s_u[sj][sb  ]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx]      = (     -  10.0 * s_u[sj][sb-5] 
                                            +  61.0 * s_u[sj][sb-4] 
                                            - 156.0 * s_u[sj][sb-3] 
                                            + 214.0 * s_u[sj][sb-2] 
                                            - 154.0 * s_u[sj][sb-1] 
                                            +  45.0 * s_u[sj][sb  ] 
                                            ) * idx_sqrd_by_12;

            }

        }
        
        return ;
    }

    /**
     * @brief Apply the 6th order 1st deriv y-deriv
     * @param Du : out - applied stencil
     * @param u   : in - input vector with padding
     * @param dx  : dx for the block
     * @param sz  : sz of the block (x,y,z) dimension
     * @param bflag : boundary flag. 
     * @return void 
     * */
    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __deriv644_y(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk)
    {
        // Coalaced mem. access - high bandwidth
        const DEVICE_INT i  = GPUDevice::block_id_x() * GPUDevice::block_dim_x() + GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::block_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();  
        const DEVICE_INT sj = GPUDevice::thread_id_y();

        const DEVICE_REAL idx = 1.0 / blk->m_dx[1];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx        = blk->m_aligned_sz[0];
        const DEVICE_INT ny        = blk->m_aligned_sz[1];
        const DEVICE_INT nz        = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        // y direction
        SHARED_MEM DEVICE_REAL s_u[pencil_sz + 2 * pw][pencils];

        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sj][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = sj + pw;

        globalIdx = k * nx * ny + (sb%ny) * nx + i;

        Du[globalIdx] = (  -         s_u[sb -3][si] 
                            +  9.0 * s_u[sb -2][si]
                            - 45.0 * s_u[sb -1][si]
                            + 45.0 * s_u[sb +1][si]
                            -  9.0 * s_u[sb +2][si]
                            +        s_u[sb +3][si]) * idx_by_60;

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_DOWN)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            if(sb == pw)
            {
                Du[globalIdx]               =  (    - 25.0 * s_u[sb  ][si]
                                                    + 48.0 * s_u[sb+1][si]
                                                    - 36.0 * s_u[sb+2][si]
                                                    + 16.0 * s_u[sb+3][si]
                                                    -  3.0 * s_u[sb+4][si]
                                                ) * idx_by_12;

                Du[globalIdx + 1*nx]        =  (    -  3.0 * s_u[sb  ][si]
                                                    - 10.0 * s_u[sb+1][si]
                                                    + 18.0 * s_u[sb+2][si]
                                                    -  6.0 * s_u[sb+3][si]
                                                    +        s_u[sb+4][si]
                                                ) * idx_by_12;

                Du[globalIdx + 2*nx]        =  (    +       s_u[sb  ][si]
                                                    - 8.0 * s_u[sb+1][si]
                                                    + 8.0 * s_u[sb+3][si]
                                                    -       s_u[sb+4][si]
                                                ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_UINT ie = actual_ny-pw;
            if(sb == ie-1)
            {
                Du[globalIdx- 2*nx]   =     (   +       s_u[sb-4][si]
                                                - 8.0 * s_u[sb-3][si]
                                                + 8.0 * s_u[sb-1][si]
                                                -       s_u[sb  ][si]
                                                ) * idx_by_12;

                Du[globalIdx- 1*nx]   =     (   -        s_u[sb-4][si]
                                                +  6.0 * s_u[sb-3][si]
                                                - 18.0 * s_u[sb-2][si]
                                                + 10.0 * s_u[sb-1][si]
                                                +  3.0 * s_u[sb  ][si]
                                                ) * idx_by_12;
                                    
                Du[globalIdx]         =     (      3.0 * s_u[sb-4][si]
                                                - 16.0 * s_u[sb-3][si]
                                                + 36.0 * s_u[sb-2][si]
                                                - 48.0 * s_u[sb-1][si]
                                                + 25.0 * s_u[sb  ][si]
                                                ) * idx_by_12;
            }

        }
        
        return;

    }

    /**
     * @brief Apply the 6th order 2nd deriv y-deriv
     * @param DDu : out - applied stencil
     * @param u   : in - input vector with padding
     * @param dx  : dx for the block
     * @param sz  : sz of the block (x,y,z) dimension
     * @param bflag : boundary flag. 
     * @return void 
     * */
    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __deriv644_yy(DEVICE_REAL * const  DDu, const DEVICE_REAL * const  u, const BlockGPU3D* const blk)
    {
        // Coalaced mem. access - high bandwidth
        const DEVICE_INT i  = GPUDevice::block_id_x() * GPUDevice::block_dim_x() + GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::block_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();  
        const DEVICE_INT sj = GPUDevice::thread_id_y();

        const DEVICE_REAL idx_sqrd = 1.0 / (blk->m_dx[1]*blk->m_dx[1]);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;
        
        const DEVICE_INT nx        = blk->m_aligned_sz[0];
        const DEVICE_INT ny        = blk->m_aligned_sz[1];
        const DEVICE_INT nz        = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        // y direction
        SHARED_MEM DEVICE_REAL s_u[pencil_sz + 2 * pw][pencils];

        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sj][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = sj + pw;

        globalIdx = k * nx * ny + (sb%ny) * nx + i;

        DDu[globalIdx] = (      2.0 * s_u[sb-3][si]
                            -  27.0 * s_u[sb-2][si]
                            + 270.0 * s_u[sb-1][si]
                            - 490.0 * s_u[sb  ][si]
                            + 270.0 * s_u[sb+1][si]
                            -  27.0 * s_u[sb+2][si]
                            +   2.0 * s_u[sb+3][si]
                    ) * idx_sqrd_by_180;

        const DEVICE_UINT bflag = blk->m_bflag;
    
        if (bflag & (1u<<OCT_DIR_DOWN)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12  = idx_sqrd / 12.0;
            if(sb == pw )
            {
                DDu[globalIdx]      = (   45.0  * s_u[sb    ][si]
                                        - 154.0 * s_u[sb + 1][si]
                                        + 214.0 * s_u[sb + 2][si]
                                        - 156.0 * s_u[sb + 3][si]
                                        +  61.0 * s_u[sb + 4][si]
                                        -  10.0 * s_u[sb + 5][si]
                                        ) * idx_sqrd_by_12;
            
                DDu[globalIdx+1*nx] =   (  10.0 * s_u[sb    ][si]
                                        -  15.0 * s_u[sb + 1][si]
                                        -   4.0 * s_u[sb + 2][si]
                                        +  14.0 * s_u[sb + 3][si]
                                        -   6.0 * s_u[sb + 4][si]
                                        +         s_u[sb + 5][si]
                                        ) * idx_sqrd_by_12;
                
                DDu[globalIdx+2*nx] =       ( -   s_u[sb    ][si]
                                        +  16.0 * s_u[sb + 1][si]
                                        -  30.0 * s_u[sb + 2][si]
                                        +  16.0 * s_u[sb + 3][si]
                                        -         s_u[sb + 4][si]
                                        ) * idx_sqrd_by_12;
            }
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12  = idx_sqrd / 12.0;
            const DEVICE_UINT ie              = actual_ny - pw;

            if(sb == ie-1)
            {
                DDu[globalIdx-2*nx] = ( -         s_u[sb -4][si]
                                        +  16.0 * s_u[sb -3][si]
                                        -  30.0 * s_u[sb -2][si]
                                        +  16.0 * s_u[sb -1][si]
                                        -         s_u[sb   ][si]
                                        ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx-1*nx] = (           s_u[sb-5][si]
                                        -   6.0 * s_u[sb-4][si]
                                        +  14.0 * s_u[sb-3][si]
                                        -   4.0 * s_u[sb-2][si]
                                        -  15.0 * s_u[sb-1][si]
                                        +  10.0 * s_u[sb  ][si]
                                        ) * idx_sqrd_by_12;
                
                DDu[globalIdx]      = ( -  10.0 * s_u[sb-5][si] 
                                        +  61.0 * s_u[sb-4][si] 
                                        - 156.0 * s_u[sb-3][si] 
                                        + 214.0 * s_u[sb-2][si] 
                                        - 154.0 * s_u[sb-1][si] 
                                        +  45.0 * s_u[sb  ][si] 
                                        ) * idx_sqrd_by_12;

            }

        }
        
        return;

    }

    /**
     * @brief Apply the 6th order 1st deriv z-deriv
     * @param Du : out - applied stencil
     * @param u   : in - input vector with padding
     * @param dx  : dx for the block
     * @param sz  : sz of the block (x,y,z) dimension
     * @param bflag : boundary flag. 
     * @return void 
     * */
    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __deriv644_z(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk)
    {
        // coalaced mem. access.
        const DEVICE_INT i  = GPUDevice::block_id_x() * GPUDevice::block_dim_x() + GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::block_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();
        const DEVICE_INT sk = GPUDevice::thread_id_y(); 

        const DEVICE_REAL idx = 1.0 / blk->m_dx[2];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx        = blk->m_aligned_sz[0];
        const DEVICE_INT ny        = blk->m_aligned_sz[1];
        const DEVICE_INT nz        = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        // z direction
        SHARED_MEM DEVICE_REAL s_u[pencil_sz + 2 * pw][pencils];

        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sk][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = sk + pw;
        globalIdx = (sb%nz) * nx * ny + j * nx + i;
        Du[globalIdx] = (  -         s_u[sb -3][si] 
                            +  9.0 * s_u[sb -2][si]
                            - 45.0 * s_u[sb -1][si]
                            + 45.0 * s_u[sb +1][si]
                            -  9.0 * s_u[sb +2][si]
                            +        s_u[sb +3][si]) * idx_by_60;

        
        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            if(sb == pw)
            {
                Du[globalIdx]               =  (    - 25.0 * s_u[sb  ][si]
                                                    + 48.0 * s_u[sb+1][si]
                                                    - 36.0 * s_u[sb+2][si]
                                                    + 16.0 * s_u[sb+3][si]
                                                    -  3.0 * s_u[sb+4][si]
                                                ) * idx_by_12;

                Du[globalIdx + 1*nx*ny]     =  (    -  3.0 * s_u[sb  ][si]
                                                    - 10.0 * s_u[sb+1][si]
                                                    + 18.0 * s_u[sb+2][si]
                                                    -  6.0 * s_u[sb+3][si]
                                                    +        s_u[sb+4][si]
                                                ) * idx_by_12;

                Du[globalIdx + 2*nx*ny]     =  ( +      s_u[sb  ][si]
                                                - 8.0 * s_u[sb+1][si]
                                                + 8.0 * s_u[sb+3][si]
                                                -       s_u[sb+4][si]
                                                ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_UINT ie = actual_nz-pw;
            if(sb == ie-1)
            {
                Du[globalIdx-2*nx*ny]   = ( +       s_u[sb-4][si]
                                            - 8.0 * s_u[sb-3][si]
                                            + 8.0 * s_u[sb-1][si]
                                            -       s_u[sb  ][si]
                                            ) * idx_by_12;

                Du[globalIdx-1*nx*ny]   = ( -         s_u[sb-4][si]
                                            +  6.0 * s_u[sb-3][si]
                                            - 18.0 * s_u[sb-2][si]
                                            + 10.0 * s_u[sb-1][si]
                                            +  3.0 * s_u[sb  ][si]
                                            ) * idx_by_12;
                                    
                Du[globalIdx]           = (    3.0 * s_u[sb-4][si]
                                            - 16.0 * s_u[sb-3][si]
                                            + 36.0 * s_u[sb-2][si]
                                            - 48.0 * s_u[sb-1][si]
                                            + 25.0 * s_u[sb  ][si]
                                            ) * idx_by_12;
            }

        }

        return;
    }
    
    /**
     * @brief Apply the 6th order 2nd deriv z-deriv
     * @param DDu : out - applied stencil
     * @param u   : in - input vector with padding
     * @param dx  : dx for the block
     * @param sz  : sz of the block (x,y,z) dimension
     * @param bflag : boundary flag. 
     * @return void 
     * */
    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __deriv644_zz(DEVICE_REAL * const  DDu, const DEVICE_REAL * const  u, const BlockGPU3D* const blk)
    {

        // coalaced mem. access.
        const DEVICE_INT i  = GPUDevice::block_id_x() * GPUDevice::block_dim_x() + GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::block_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();
        const DEVICE_INT sk = GPUDevice::thread_id_y(); 

        const DEVICE_REAL idx_sqrd = 1.0/(blk->m_dx[2]*blk->m_dx[2]);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;
        
        const DEVICE_INT nx        = blk->m_aligned_sz[0];
        const DEVICE_INT ny        = blk->m_aligned_sz[1];
        const DEVICE_INT nz        = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        // z direction
        SHARED_MEM DEVICE_REAL s_u[pencil_sz + 2 * pw][pencils];

        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sk][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = sk + pw;
        globalIdx = (sb%nz) * nx * ny + j * nx + i;
        
        DDu[globalIdx] = (      2.0 * s_u[sb-3][si]
                            -  27.0 * s_u[sb-2][si]
                            + 270.0 * s_u[sb-1][si]
                            - 490.0 * s_u[sb  ][si]
                            + 270.0 * s_u[sb+1][si]
                            -  27.0 * s_u[sb+2][si]
                            +   2.0 * s_u[sb+3][si]
                            ) * idx_sqrd_by_180;

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12  = idx_sqrd / 12.0;
            if(sb == pw )
            {
                DDu[globalIdx]          =  (  45.0  * s_u[sb    ][si]
                                            - 154.0 * s_u[sb + 1][si]
                                            + 214.0 * s_u[sb + 2][si]
                                            - 156.0 * s_u[sb + 3][si]
                                            +  61.0 * s_u[sb + 4][si]
                                            -  10.0 * s_u[sb + 5][si]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx+1*nx*ny]  = (    10.0 * s_u[sb    ][si]
                                            -  15.0 * s_u[sb + 1][si]
                                            -   4.0 * s_u[sb + 2][si]
                                            +  14.0 * s_u[sb + 3][si]
                                            -   6.0 * s_u[sb + 4][si]
                                            +         s_u[sb + 5][si]
                                            ) * idx_sqrd_by_12;
                    
                DDu[globalIdx+2*nx*ny]  = ( -         s_u[sb    ][si]
                                            +  16.0 * s_u[sb + 1][si]
                                            -  30.0 * s_u[sb + 2][si]
                                            +  16.0 * s_u[sb + 3][si]
                                            -         s_u[sb + 4][si]
                                            ) * idx_sqrd_by_12;
            }
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12  = idx_sqrd / 12.0;
            const DEVICE_UINT ie              = actual_nz - pw;

            if(sb == ie-1)
            {
                DDu[globalIdx-2*nx*ny] = ( -          s_u[sb -4][si]
                                            +  16.0 * s_u[sb -3][si]
                                            -  30.0 * s_u[sb -2][si]
                                            +  16.0 * s_u[sb -1][si]
                                            -         s_u[sb   ][si]
                                            ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx-1*nx*ny] = (            s_u[sb-5][si]
                                            -   6.0 * s_u[sb-4][si]
                                            +  14.0 * s_u[sb-3][si]
                                            -   4.0 * s_u[sb-2][si]
                                            -  15.0 * s_u[sb-1][si]
                                            +  10.0 * s_u[sb  ][si]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx]         = ( -  10.0 * s_u[sb-5][si] 
                                           +  61.0 * s_u[sb-4][si] 
                                           - 156.0 * s_u[sb-3][si] 
                                           + 214.0 * s_u[sb-2][si] 
                                           - 154.0 * s_u[sb-1][si] 
                                           +  45.0 * s_u[sb  ][si] 
                                            ) * idx_sqrd_by_12;

            }

        }

        return;
    }

    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __ko_deriv42_x(DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {


        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::block_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();  
        const DEVICE_INT sj = GPUDevice::thread_id_y();
        
        
        const DEVICE_REAL dx              =  blk->m_dx[0];
        const DEVICE_REAL idx             =  1.0 / blk->m_dx[0];
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        SHARED_MEM DEVICE_REAL s_u[pencils][pencil_sz + 2*pw];
        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sj][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = si + pw;
        globalIdx = k * nx * ny + j * nx + (sb % nx);

        Du[globalIdx] = (
                                -      s_u[sj][sb-3]
                                +  6.0*s_u[sj][sb-2]
                                - 15.0*s_u[sj][sb-1]
                                + 20.0*s_u[sj][sb  ]
                                - 15.0*s_u[sj][sb+1]
                                +  6.0*s_u[sj][sb+2]
                                -      s_u[sj][sb+3]
                        ) * pre_factor_6_dx;

        const DEVICE_UINT bflag = blk->m_bflag;
        if (bflag & (1u<<OCT_DIR_LEFT)) {
            
            GPUDevice::sync_threads();
            if(sb==pw)
            {
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (      s_u[sj][sb+3]
                                    - 3.0*s_u[sj][sb+2]
                                    + 3.0*s_u[sj][sb+1]
                                    -     s_u[sj][sb]
                                    )/smr3;
                Du[globalIdx+1] =  (       s_u[sj][sb + 4]
                                    -  6.0*s_u[sj][sb + 3]
                                    + 12.0*s_u[sj][sb + 2]
                                    - 10.0*s_u[sj][sb + 1]
                                    +  3.0*s_u[sj][sb    ]
                                    )/smr2;

                Du[globalIdx+2] =  (       s_u[sj][sb + 5]
                                    -  6.0*s_u[sj][sb + 4]
                                    + 15.0*s_u[sj][sb + 3]
                                    - 19.0*s_u[sj][sb + 2]
                                    + 12.0*s_u[sj][sb + 1]
                                    -  3.0*s_u[sj][sb    ]
                                    )/smr1;

            }

        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_UINT ie = actual_nx - pw;
            if(sb == (ie-1))
            {
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2] = (        s_u[sj][sb-5]
                                    -  6.0*s_u[sj][sb-4]
                                    + 15.0*s_u[sj][sb-3]
                                    - 19.0*s_u[sj][sb-2]
                                    + 12.0*s_u[sj][sb-1]
                                    -  3.0*s_u[sj][sb  ]
                                    )/spr1;

                Du[globalIdx-1] = (        s_u[sj][sb-4]
                                    -  6.0*s_u[sj][sb-3]
                                    + 12.0*s_u[sj][sb-2]
                                    - 10.0*s_u[sj][sb-1]
                                    +  3.0*s_u[sj][sb  ]
                                    )/spr2;

                Du[globalIdx] =    (       s_u[sj][sb-3]
                                    -  3.0*s_u[sj][sb-2]
                                    +  3.0*s_u[sj][sb-1]
                                    -      s_u[sj][sb  ]
                                    )/spr3;

            }

            
        }

    }

    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __ko_deriv42_y(DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {

        const DEVICE_INT i  = GPUDevice::block_id_x() * GPUDevice::block_dim_x() + GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::block_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();  
        const DEVICE_INT sj = GPUDevice::thread_id_y();
        
        const DEVICE_REAL dx              =  blk->m_dx[1];
        const DEVICE_REAL idx             =  1.0 / blk->m_dx[1];
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        SHARED_MEM DEVICE_REAL s_u[pencil_sz + 2 * pw][pencils];
        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sj][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = sj + pw;
        globalIdx = k * nx * ny + (sb % ny) * nx + i;

        Du[globalIdx] = (
                                -      s_u[sb-3][si]
                                +  6.0*s_u[sb-2][si]
                                - 15.0*s_u[sb-1][si]
                                + 20.0*s_u[sb  ][si]
                                - 15.0*s_u[sb+1][si]
                                +  6.0*s_u[sb+2][si]
                                -      s_u[sb+3][si]
                        ) * pre_factor_6_dx;

        
        const DEVICE_UINT bflag = blk->m_bflag;
        if (bflag & (1u<<OCT_DIR_DOWN)) {
            
            GPUDevice::sync_threads();
            if(sb==pw)
            {
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (      s_u[sb+3][si]
                                    - 3.0*s_u[sb+2][si]
                                    + 3.0*s_u[sb+1][si]
                                    -     s_u[sb][si]
                                    )/smr3;

                Du[globalIdx+1*nx] =  (    s_u[sb + 4][si]
                                    -  6.0*s_u[sb + 3][si]
                                    + 12.0*s_u[sb + 2][si]
                                    - 10.0*s_u[sb + 1][si]
                                    +  3.0*s_u[sb    ][si]
                                    )/smr2;

                Du[globalIdx+2*nx] =  (    s_u[sb + 5][si]
                                    -  6.0*s_u[sb + 4][si]
                                    + 15.0*s_u[sb + 3][si]
                                    - 19.0*s_u[sb + 2][si]
                                    + 12.0*s_u[sb + 1][si]
                                    -  3.0*s_u[sb    ][si]
                                    )/smr1;

            }

        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_UINT ie = actual_ny - pw;
            if(sb == (ie-1))
            {
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2*nx] = (     s_u[sb-5][si]
                                    -  6.0*s_u[sb-4][si]
                                    + 15.0*s_u[sb-3][si]
                                    - 19.0*s_u[sb-2][si]
                                    + 12.0*s_u[sb-1][si]
                                    -  3.0*s_u[sb  ][si]
                                    )/spr1;

                Du[globalIdx-1*nx] = (     s_u[sb-4][si]
                                    -  6.0*s_u[sb-3][si]
                                    + 12.0*s_u[sb-2][si]
                                    - 10.0*s_u[sb-1][si]
                                    +  3.0*s_u[sb  ][si]
                                    )/spr2;

                Du[globalIdx] =    (       s_u[sb-3][si]
                                    -  3.0*s_u[sb-2][si]
                                    +  3.0*s_u[sb-1][si]
                                    -      s_u[sb  ][si]
                                    )/spr3;

            }

            
        }

    }

    template<int pw, int pencils, int pencil_sz>
    inline DEVICE_FUNC void __ko_deriv42_z(DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {
        // coalaced mem. access.
        const DEVICE_INT i  = GPUDevice::block_id_x() * GPUDevice::block_dim_x() + GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::block_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();
        const DEVICE_INT sk = GPUDevice::thread_id_y(); 

        const DEVICE_REAL dx              =  blk->m_dx[2];
        const DEVICE_REAL idx             =  1.0 / blk->m_dx[2];
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        SHARED_MEM DEVICE_REAL s_u[pencil_sz + 2 * pw][pencils];
        DEVICE_INT globalIdx = k * nx * ny + j * nx + i;

        s_u[sk][si] = u[globalIdx];
        GPUDevice::sync_threads();

        const DEVICE_INT sb = sk + pw;
        globalIdx = (sb % nz) * nx * ny + j * nx + i;

        Du[globalIdx] = (
                                -      s_u[sb-3][si]
                                +  6.0*s_u[sb-2][si]
                                - 15.0*s_u[sb-1][si]
                                + 20.0*s_u[sb  ][si]
                                - 15.0*s_u[sb+1][si]
                                +  6.0*s_u[sb+2][si]
                                -      s_u[sb+3][si]
                        ) * pre_factor_6_dx;

        
        const DEVICE_UINT bflag = blk->m_bflag;
        
        if (bflag & (1u<<OCT_DIR_BACK)) {
            
            GPUDevice::sync_threads();
            if(sb==pw)
            {
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (      s_u[sb+3][si]
                                    - 3.0*s_u[sb+2][si]
                                    + 3.0*s_u[sb+1][si]
                                    -     s_u[sb  ][si]
                                    )/smr3;

                Du[globalIdx+1*nx*ny] =  (     s_u[sb + 4][si]
                                        -  6.0*s_u[sb + 3][si]
                                        + 12.0*s_u[sb + 2][si]
                                        - 10.0*s_u[sb + 1][si]
                                        +  3.0*s_u[sb    ][si]
                                        )/smr2;

                Du[globalIdx+2*nx*ny] =  (     s_u[sb + 5][si]
                                        -  6.0*s_u[sb + 4][si]
                                        + 15.0*s_u[sb + 3][si]
                                        - 19.0*s_u[sb + 2][si]
                                        + 12.0*s_u[sb + 1][si]
                                        -  3.0*s_u[sb    ][si]
                                        )/smr1;

            }

        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_UINT ie = actual_nz - pw;
            if(sb == (ie-1))
            {
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2*nx*ny] = (      s_u[sb-5][si]
                                        -  6.0*s_u[sb-4][si]
                                        + 15.0*s_u[sb-3][si]
                                        - 19.0*s_u[sb-2][si]
                                        + 12.0*s_u[sb-1][si]
                                        -  3.0*s_u[sb  ][si]
                                        )/spr1;

                Du[globalIdx-1*nx*ny] = (      s_u[sb-4][si]
                                        -  6.0*s_u[sb-3][si]
                                        + 12.0*s_u[sb-2][si]
                                        - 10.0*s_u[sb-1][si]
                                        +  3.0*s_u[sb  ][si]
                                        )/spr2;

                Du[globalIdx] =    (       s_u[sb-3][si]
                                    -  3.0*s_u[sb-2][si]
                                    +  3.0*s_u[sb-1][si]
                                    -      s_u[sb  ][si]
                                    )/spr3;

            }

            
        }

    }

    #define IDX_B(i,j,k) ((i) + nx * ( (j) + nx * (k)))
    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_deriv644_x(DEVICE_REAL * const  Du, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL idx = 1.0 / blk->m_dx[0];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i+pw;  
        const DEVICE_INT sj = j;
        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 
        
        if(si < actual_nx-pw)
        for(DEVICE_INT sk = 0; sk < (actual_nz - 0); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

            Du[globalIdx] = (    -   s_u[globalIdx - 3] 
                            +  9.0 * s_u[globalIdx - 2]
                            - 45.0 * s_u[globalIdx - 1]
                            + 45.0 * s_u[globalIdx + 1]
                            -  9.0 * s_u[globalIdx + 2]
                            +        s_u[globalIdx + 3]) * idx_by_60;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_LEFT)) {

            GPUDevice::sync_threads();
            
            const DEVICE_REAL idx_by_12 = idx / 12.0;

            for(DEVICE_INT sk = 0; sk < (actual_nz - 0); sk++)
            if(si == pw)
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                Du[globalIdx] = (  - 25.0  * s_u[globalIdx    ]
                                    + 48.0 * s_u[globalIdx + 1]
                                    - 36.0 * s_u[globalIdx + 2]
                                    + 16.0 * s_u[globalIdx + 3]
                                    -  3.0 * s_u[globalIdx + 4]
                                ) * idx_by_12;

                Du[globalIdx +1]  =  (  -  3.0 * s_u[globalIdx    ]
                                        - 10.0 * s_u[globalIdx + 1]
                                        + 18.0 * s_u[globalIdx + 2]
                                        -  6.0 * s_u[globalIdx + 3]
                                        +        s_u[globalIdx + 4]
                                    ) * idx_by_12;

                Du[globalIdx+2]   =  (   +      s_u[globalIdx    ]
                                        - 8.0 * s_u[globalIdx + 1]
                                        + 8.0 * s_u[globalIdx + 3]
                                        -       s_u[globalIdx + 4]
                                    ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT ie = actual_nx-pw;
            for(DEVICE_INT sk = 0; sk < (actual_nz - 0); sk++)
            if(si == ie-1)
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                Du[globalIdx-2]   = (  +        s_u[globalIdx - 4]
                                        - 8.0 * s_u[globalIdx - 3]
                                        + 8.0 * s_u[globalIdx - 1]
                                        -       s_u[globalIdx    ]
                                        ) * idx_by_12;

                Du[globalIdx-1] =  (  -          s_u[globalIdx - 4]
                                        +  6.0 * s_u[globalIdx - 3]
                                        - 18.0 * s_u[globalIdx - 2]
                                        + 10.0 * s_u[globalIdx - 1]
                                        +  3.0 * s_u[globalIdx    ]
                                        ) * idx_by_12;
                                    
                Du[globalIdx] =    (       3.0 * s_u[globalIdx - 4]
                                        - 16.0 * s_u[globalIdx - 3]
                                        + 36.0 * s_u[globalIdx - 2]
                                        - 48.0 * s_u[globalIdx - 1]
                                        + 25.0 * s_u[globalIdx    ]
                                        ) * idx_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_deriv644_y(DEVICE_REAL * const  Du, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL idx = 1.0 / blk->m_dx[1];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i;  
        const DEVICE_INT sj = j+pw;

        // stride for y
        const DEVICE_INT ss = nx;
        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 

        if(sj < (actual_ny -pw)  /*&& si < (actual_ny -pw)*/)
        for(DEVICE_INT sk = 0; sk < (actual_nz - 0); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

            Du[globalIdx] = (    -   s_u[globalIdx - 3 * ss] 
                            +  9.0 * s_u[globalIdx - 2 * ss]
                            - 45.0 * s_u[globalIdx - 1 * ss]
                            + 45.0 * s_u[globalIdx + 1 * ss]
                            -  9.0 * s_u[globalIdx + 2 * ss]
                            +        s_u[globalIdx + 3 * ss]) * idx_by_60;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_DOWN)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            for(DEVICE_INT sk = 0; sk < (actual_nz - 0); sk++)
            if(sj == pw && si < (actual_nx -pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                Du[globalIdx          ] = (  - 25.0  * s_u[globalIdx       ]
                                            + 48.0 * s_u[globalIdx + 1 * ss]
                                            - 36.0 * s_u[globalIdx + 2 * ss]
                                            + 16.0 * s_u[globalIdx + 3 * ss]
                                            -  3.0 * s_u[globalIdx + 4 * ss]
                                        ) * idx_by_12;

                Du[globalIdx + 1 * ss]  =  (    -  3.0 * s_u[globalIdx         ]
                                                - 10.0 * s_u[globalIdx + 1 * ss]
                                                + 18.0 * s_u[globalIdx + 2 * ss]
                                                -  6.0 * s_u[globalIdx + 3 * ss]
                                                +        s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;

                Du[globalIdx + 2 * ss]   =  (    +      s_u[globalIdx         ]
                                                - 8.0 * s_u[globalIdx + 1 * ss]
                                                + 8.0 * s_u[globalIdx + 3 * ss]
                                                -       s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT ie = actual_ny-pw;
            for(DEVICE_INT sk = 0; sk < (actual_nz - 0); sk++)
            if(sj == ie-1 && si < (actual_nx -pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                Du[globalIdx-2 * ss]   = (  +       s_u[globalIdx - 4 * ss]
                                            - 8.0 * s_u[globalIdx - 3 * ss]
                                            + 8.0 * s_u[globalIdx - 1 * ss]
                                            -       s_u[globalIdx         ]
                                            ) * idx_by_12;

                Du[globalIdx-1 * ss] =  (   -        s_u[globalIdx - 4 * ss]
                                            +  6.0 * s_u[globalIdx - 3 * ss]
                                            - 18.0 * s_u[globalIdx - 2 * ss]
                                            + 10.0 * s_u[globalIdx - 1 * ss]
                                            +  3.0 * s_u[globalIdx         ]
                                            ) * idx_by_12;
                                    
                Du[globalIdx] =    (       3.0 * s_u[globalIdx - 4 * ss]
                                        - 16.0 * s_u[globalIdx - 3 * ss]
                                        + 36.0 * s_u[globalIdx - 2 * ss]
                                        - 48.0 * s_u[globalIdx - 1 * ss]
                                        + 25.0 * s_u[globalIdx         ]
                                        ) * idx_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_deriv644_z(DEVICE_REAL * const  Du, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL idx = 1.0 / blk->m_dx[2];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i;//+pw;  
        const DEVICE_INT sj = j;//+pw;

        // stride for y
        const DEVICE_INT ss = nx * ny;
        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 

        //if(sj < (actual_ny -pw)  && si < (actual_ny -pw))
        for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

            Du[globalIdx] = (    -   s_u[globalIdx - 3 * ss] 
                            +  9.0 * s_u[globalIdx - 2 * ss]
                            - 45.0 * s_u[globalIdx - 1 * ss]
                            + 45.0 * s_u[globalIdx + 1 * ss]
                            -  9.0 * s_u[globalIdx + 2 * ss]
                            +        s_u[globalIdx + 3 * ss]) * idx_by_60;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            //for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            const DEVICE_INT sk = pw;
            if(sk == pw && sj < (actual_ny -pw)  && si < (actual_ny -pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                Du[globalIdx          ] = (  - 25.0  * s_u[globalIdx       ]
                                            + 48.0 * s_u[globalIdx + 1 * ss]
                                            - 36.0 * s_u[globalIdx + 2 * ss]
                                            + 16.0 * s_u[globalIdx + 3 * ss]
                                            -  3.0 * s_u[globalIdx + 4 * ss]
                                        ) * idx_by_12;

                Du[globalIdx + 1 * ss]  =  (    -  3.0 * s_u[globalIdx         ]
                                                - 10.0 * s_u[globalIdx + 1 * ss]
                                                + 18.0 * s_u[globalIdx + 2 * ss]
                                                -  6.0 * s_u[globalIdx + 3 * ss]
                                                +        s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;

                Du[globalIdx + 2 * ss]   =  (    +      s_u[globalIdx         ]
                                                - 8.0 * s_u[globalIdx + 1 * ss]
                                                + 8.0 * s_u[globalIdx + 3 * ss]
                                                -       s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT ie = actual_nz-pw;
            //for(DEVICE_INT sk = 0; sk < (actual_nz - 0); sk++)
            DEVICE_INT sk = ie-1;
            if(sk == ie-1 && sj < (actual_ny -pw)  && si < (actual_ny -pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                Du[globalIdx-2 * ss]   = (  +       s_u[globalIdx - 4 * ss]
                                            - 8.0 * s_u[globalIdx - 3 * ss]
                                            + 8.0 * s_u[globalIdx - 1 * ss]
                                            -       s_u[globalIdx         ]
                                            ) * idx_by_12;

                Du[globalIdx-1 * ss] =  (   -        s_u[globalIdx - 4 * ss]
                                            +  6.0 * s_u[globalIdx - 3 * ss]
                                            - 18.0 * s_u[globalIdx - 2 * ss]
                                            + 10.0 * s_u[globalIdx - 1 * ss]
                                            +  3.0 * s_u[globalIdx         ]
                                            ) * idx_by_12;
                                    
                Du[globalIdx] =    (       3.0 * s_u[globalIdx - 4 * ss]
                                        - 16.0 * s_u[globalIdx - 3 * ss]
                                        + 36.0 * s_u[globalIdx - 2 * ss]
                                        - 48.0 * s_u[globalIdx - 1 * ss]
                                        + 25.0 * s_u[globalIdx         ]
                                        ) * idx_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_deriv644_xx(DEVICE_REAL * const  DDu, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        const DEVICE_REAL dx              =  blk->m_dx[0];
        const DEVICE_REAL idx_sqrd = 1.0 / (dx * dx);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i + pw;  
        const DEVICE_INT sj = j ;//+ pw;

        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 
        
        if(si < (actual_nx-pw)  /*&&  sj < (actual_ny-pw)*/)
        for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

            DDu[globalIdx] = (      2.0 * s_u[globalIdx - 3]
                                -  27.0 * s_u[globalIdx - 2]
                                + 270.0 * s_u[globalIdx - 1]
                                - 490.0 * s_u[globalIdx    ]
                                + 270.0 * s_u[globalIdx + 1]
                                -  27.0 * s_u[globalIdx + 2]
                                +   2.0 * s_u[globalIdx + 3]
                                ) * idx_sqrd_by_180;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_LEFT)) {

            GPUDevice::sync_threads();
            
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;

            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(si == pw &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                DDu[globalIdx]      = (       45.0  * s_u[globalIdx    ]
                                            - 154.0 * s_u[globalIdx + 1]
                                            + 214.0 * s_u[globalIdx + 2]
                                            - 156.0 * s_u[globalIdx + 3]
                                            +  61.0 * s_u[globalIdx + 4]
                                            -  10.0 * s_u[globalIdx + 5]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx+1]    = (        10.0 * s_u[globalIdx    ]
                                            -  15.0 * s_u[globalIdx + 1]
                                            -   4.0 * s_u[globalIdx + 2]
                                            +  14.0 * s_u[globalIdx + 3]
                                            -   6.0 * s_u[globalIdx + 4]
                                            +         s_u[globalIdx + 5]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx+2]    = (     -         s_u[globalIdx    ]
                                            +  16.0 * s_u[globalIdx + 1]
                                            -  30.0 * s_u[globalIdx + 2]
                                            +  16.0 * s_u[globalIdx + 3]
                                            -         s_u[globalIdx + 4]
                                            ) * idx_sqrd_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT ie = actual_nx-pw;
            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(si == ie-1  &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                DDu[globalIdx-2]    = (     -         s_u[globalIdx -4]
                                            +  16.0 * s_u[globalIdx -3]
                                            -  30.0 * s_u[globalIdx -2]
                                            +  16.0 * s_u[globalIdx -1]
                                            -         s_u[globalIdx   ]
                                            ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx-1]    = (               s_u[globalIdx - 5]
                                            -   6.0 * s_u[globalIdx - 4]
                                            +  14.0 * s_u[globalIdx - 3]
                                            -   4.0 * s_u[globalIdx - 2]
                                            -  15.0 * s_u[globalIdx - 1]
                                            +  10.0 * s_u[globalIdx    ]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx]      = (     -  10.0 * s_u[globalIdx - 5] 
                                            +  61.0 * s_u[globalIdx - 4] 
                                            - 156.0 * s_u[globalIdx - 3] 
                                            + 214.0 * s_u[globalIdx - 2] 
                                            - 154.0 * s_u[globalIdx - 1] 
                                            +  45.0 * s_u[globalIdx    ] 
                                            ) * idx_sqrd_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_deriv644_yy(DEVICE_REAL * const  DDu, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        const DEVICE_REAL dx              =  blk->m_dx[1];
        const DEVICE_REAL idx_sqrd = 1.0 / (dx * dx);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i ;//+ pw;  
        const DEVICE_INT sj = j + pw;

        const DEVICE_INT ss = nx;

        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 
        
        if(sj < (actual_ny-pw))
        for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

            DDu[globalIdx] = (      2.0 * s_u[globalIdx - 3 * ss]
                                -  27.0 * s_u[globalIdx - 2 * ss]
                                + 270.0 * s_u[globalIdx - 1 * ss]
                                - 490.0 * s_u[globalIdx         ]
                                + 270.0 * s_u[globalIdx + 1 * ss]
                                -  27.0 * s_u[globalIdx + 2 * ss]
                                +   2.0 * s_u[globalIdx + 3 * ss]
                                ) * idx_sqrd_by_180;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_DOWN)) {

            GPUDevice::sync_threads();
            
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;

            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(sj == pw &&  si < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                DDu[globalIdx]      = (       45.0  * s_u[globalIdx    ]
                                            - 154.0 * s_u[globalIdx + 1 * ss]
                                            + 214.0 * s_u[globalIdx + 2 * ss]
                                            - 156.0 * s_u[globalIdx + 3 * ss]
                                            +  61.0 * s_u[globalIdx + 4 * ss]
                                            -  10.0 * s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx + 1 * ss]    = (        10.0 * s_u[globalIdx    ]
                                            -  15.0 * s_u[globalIdx + 1 * ss]
                                            -   4.0 * s_u[globalIdx + 2 * ss]
                                            +  14.0 * s_u[globalIdx + 3 * ss]
                                            -   6.0 * s_u[globalIdx + 4 * ss]
                                            +         s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx + 2 * ss]    = (     -         s_u[globalIdx    ]
                                            +  16.0 * s_u[globalIdx + 1 * ss]
                                            -  30.0 * s_u[globalIdx + 2 * ss]
                                            +  16.0 * s_u[globalIdx + 3 * ss]
                                            -         s_u[globalIdx + 4 * ss]
                                            ) * idx_sqrd_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT ie = actual_ny-pw;
            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(sj == ie-1 && si < (actual_nx-pw) )
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                DDu[globalIdx - 2 * ss ]    = (     -         s_u[globalIdx -4 * ss]
                                                    +  16.0 * s_u[globalIdx -3 * ss]
                                                    -  30.0 * s_u[globalIdx -2 * ss]
                                                    +  16.0 * s_u[globalIdx -1 * ss]
                                                    -         s_u[globalIdx   ]
                                                    ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx -1 * ss]    =   (               s_u[globalIdx - 5 * ss]
                                                    -   6.0 * s_u[globalIdx - 4 * ss]
                                                    +  14.0 * s_u[globalIdx - 3 * ss]
                                                    -   4.0 * s_u[globalIdx - 2 * ss]
                                                    -  15.0 * s_u[globalIdx - 1 * ss]
                                                    +  10.0 * s_u[globalIdx    ]
                                                    ) * idx_sqrd_by_12;
                
                DDu[globalIdx]      = (     -  10.0 * s_u[globalIdx - 5 * ss] 
                                            +  61.0 * s_u[globalIdx - 4 * ss] 
                                            - 156.0 * s_u[globalIdx - 3 * ss] 
                                            + 214.0 * s_u[globalIdx - 2 * ss] 
                                            - 154.0 * s_u[globalIdx - 1 * ss] 
                                            +  45.0 * s_u[globalIdx    ] 
                                            ) * idx_sqrd_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_deriv644_zz(DEVICE_REAL * const  DDu, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        const DEVICE_REAL dx              =  blk->m_dx[2];
        const DEVICE_REAL idx_sqrd = 1.0 / (dx * dx);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i ;//+ pw;  
        const DEVICE_INT sj = j ;//+ pw;

        const DEVICE_INT ss = nx * ny ;

        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 
        
        //if(si < (actual_nx-pw)  &&  sj < (actual_ny-pw))
        for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

            DDu[globalIdx] = (      2.0 * s_u[globalIdx - 3 * ss]
                                -  27.0 * s_u[globalIdx - 2 * ss]
                                + 270.0 * s_u[globalIdx - 1 * ss]
                                - 490.0 * s_u[globalIdx         ]
                                + 270.0 * s_u[globalIdx + 1 * ss]
                                -  27.0 * s_u[globalIdx + 2 * ss]
                                +   2.0 * s_u[globalIdx + 3 * ss]
                                ) * idx_sqrd_by_180;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {

            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            //for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            const DEVICE_INT sk = pw;
            if(sk == pw && si < (actual_nx-pw)  &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                DDu[globalIdx]      = (       45.0  * s_u[globalIdx    ]
                                            - 154.0 * s_u[globalIdx + 1 * ss]
                                            + 214.0 * s_u[globalIdx + 2 * ss]
                                            - 156.0 * s_u[globalIdx + 3 * ss]
                                            +  61.0 * s_u[globalIdx + 4 * ss]
                                            -  10.0 * s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx + 1 * ss]    = (        10.0 * s_u[globalIdx    ]
                                            -  15.0 * s_u[globalIdx + 1 * ss]
                                            -   4.0 * s_u[globalIdx + 2 * ss]
                                            +  14.0 * s_u[globalIdx + 3 * ss]
                                            -   6.0 * s_u[globalIdx + 4 * ss]
                                            +         s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx + 2 * ss]    = (     -         s_u[globalIdx    ]
                                            +  16.0 * s_u[globalIdx + 1 * ss]
                                            -  30.0 * s_u[globalIdx + 2 * ss]
                                            +  16.0 * s_u[globalIdx + 3 * ss]
                                            -         s_u[globalIdx + 4 * ss]
                                            ) * idx_sqrd_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT ie = actual_nz-pw;
            //for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            const DEVICE_INT sk = ie-1;
            if(sk == ie-1 && si < (actual_nx-pw)  &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                DDu[globalIdx - 2 * ss ]    = (     -         s_u[globalIdx -4 * ss]
                                                    +  16.0 * s_u[globalIdx -3 * ss]
                                                    -  30.0 * s_u[globalIdx -2 * ss]
                                                    +  16.0 * s_u[globalIdx -1 * ss]
                                                    -         s_u[globalIdx   ]
                                                    ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx -1 * ss]    =   (               s_u[globalIdx - 5 * ss]
                                                    -   6.0 * s_u[globalIdx - 4 * ss]
                                                    +  14.0 * s_u[globalIdx - 3 * ss]
                                                    -   4.0 * s_u[globalIdx - 2 * ss]
                                                    -  15.0 * s_u[globalIdx - 1 * ss]
                                                    +  10.0 * s_u[globalIdx    ]
                                                    ) * idx_sqrd_by_12;
                
                DDu[globalIdx         ]    =  (     -  10.0 * s_u[globalIdx - 5 * ss] 
                                                    +  61.0 * s_u[globalIdx - 4 * ss] 
                                                    - 156.0 * s_u[globalIdx - 3 * ss] 
                                                    + 214.0 * s_u[globalIdx - 2 * ss] 
                                                    - 154.0 * s_u[globalIdx - 1 * ss] 
                                                    +  45.0 * s_u[globalIdx    ] 
                                                    ) * idx_sqrd_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_ko_deriv42_x(DEVICE_REAL * const  Du, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL dx              =  blk->m_dx[0];
        const DEVICE_REAL idx             =  1.0 / dx;
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i + pw;  
        const DEVICE_INT sj = j ;//+ pw;

        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 
        
        if(si < (actual_nx-pw)  /*&&  sj < (actual_ny-pw)*/)
        for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

           Du[globalIdx] = (    -      s_u[globalIdx - 3]
                                +  6.0*s_u[globalIdx - 2]
                                - 15.0*s_u[globalIdx - 1]
                                + 20.0*s_u[globalIdx    ]
                                - 15.0*s_u[globalIdx + 1]
                                +  6.0*s_u[globalIdx + 2]
                                -      s_u[globalIdx + 3]
                        ) * pre_factor_6_dx;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_LEFT)) {

            GPUDevice::sync_threads();
            
            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(si == pw &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (      s_u[globalIdx +3]
                                    - 3.0*s_u[globalIdx +2]
                                    + 3.0*s_u[globalIdx +1]
                                    -     s_u[globalIdx ]
                                    )/smr3;
                Du[globalIdx+1] =  (       s_u[globalIdx +  4]
                                    -  6.0*s_u[globalIdx +  3]
                                    + 12.0*s_u[globalIdx +  2]
                                    - 10.0*s_u[globalIdx +  1]
                                    +  3.0*s_u[globalIdx     ]
                                    )/smr2;

                Du[globalIdx+2] =  (       s_u[globalIdx + 5]
                                    -  6.0*s_u[globalIdx + 4]
                                    + 15.0*s_u[globalIdx + 3]
                                    - 19.0*s_u[globalIdx + 2]
                                    + 12.0*s_u[globalIdx + 1]
                                    -  3.0*s_u[globalIdx    ]
                                    )/smr1;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_INT ie = actual_nx-pw;
            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(si == ie-1 &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2] = (        s_u[globalIdx - 5]
                                    -  6.0*s_u[globalIdx - 4]
                                    + 15.0*s_u[globalIdx - 3]
                                    - 19.0*s_u[globalIdx - 2]
                                    + 12.0*s_u[globalIdx - 1]
                                    -  3.0*s_u[globalIdx    ]
                                    )/spr1;

                Du[globalIdx-1] = (        s_u[globalIdx - 4]
                                    -  6.0*s_u[globalIdx - 3]
                                    + 12.0*s_u[globalIdx - 2]
                                    - 10.0*s_u[globalIdx - 1]
                                    +  3.0*s_u[globalIdx    ]
                                    )/spr2;

                Du[globalIdx] =    (       s_u[globalIdx - 3]
                                    -  3.0*s_u[globalIdx - 2]
                                    +  3.0*s_u[globalIdx - 1]
                                    -      s_u[globalIdx    ]
                                    )/spr3;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_ko_deriv42_y(DEVICE_REAL * const  Du, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL dx              =  blk->m_dx[1];
        const DEVICE_REAL idx             =  1.0 / dx;
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i ;//+ pw;  
        const DEVICE_INT sj = j + pw;
        
        const DEVICE_INT ss = nx;

        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 

        if(/*si < (actual_nx-pw)  &&*/  sj < (actual_ny-pw))
        for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

           Du[globalIdx] = (    -      s_u[globalIdx - 3 * ss]
                                +  6.0*s_u[globalIdx - 2 * ss]
                                - 15.0*s_u[globalIdx - 1 * ss]
                                + 20.0*s_u[globalIdx         ]
                                - 15.0*s_u[globalIdx + 1 * ss]
                                +  6.0*s_u[globalIdx + 2 * ss]
                                -      s_u[globalIdx + 3 * ss]
                        ) * pre_factor_6_dx;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_DOWN)) {

            GPUDevice::sync_threads();
            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(sj == pw && si < (actual_nx-pw) )
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (          s_u[globalIdx +3 * ss]
                                        - 3.0*s_u[globalIdx +2 * ss]
                                        + 3.0*s_u[globalIdx +1 * ss]
                                        -     s_u[globalIdx ]
                                        )/smr3;

                Du[globalIdx + 1 * ss] =   (       s_u[globalIdx +  4 * ss]
                                            -  6.0*s_u[globalIdx +  3 * ss]
                                            + 12.0*s_u[globalIdx +  2 * ss]
                                            - 10.0*s_u[globalIdx +  1 * ss]
                                            +  3.0*s_u[globalIdx     ]
                                            )/smr2;

                Du[globalIdx + 2 *ss ] =   (       s_u[globalIdx + 5 * ss]
                                            -  6.0*s_u[globalIdx + 4 * ss]
                                            + 15.0*s_u[globalIdx + 3 * ss]
                                            - 19.0*s_u[globalIdx + 2 * ss]
                                            + 12.0*s_u[globalIdx + 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/smr1;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_INT ie = actual_ny-pw;
            for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            if(sj == ie-1 && si < (actual_nx-pw) )
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2 * ss] =    (        s_u[globalIdx - 5 * ss]
                                            -  6.0*s_u[globalIdx - 4 * ss]
                                            + 15.0*s_u[globalIdx - 3 * ss]
                                            - 19.0*s_u[globalIdx - 2 * ss]
                                            + 12.0*s_u[globalIdx - 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/spr1;

                Du[globalIdx -1 *ss] =    (        s_u[globalIdx - 4 * ss]
                                            -  6.0*s_u[globalIdx - 3 * ss]
                                            + 12.0*s_u[globalIdx - 2 * ss]
                                            - 10.0*s_u[globalIdx - 1 * ss]
                                            +  3.0*s_u[globalIdx    ]
                                            )/spr2;

                Du[globalIdx] =    (       s_u[globalIdx - 3 * ss]
                                    -  3.0*s_u[globalIdx - 2 * ss]
                                    +  3.0*s_u[globalIdx - 1 * ss]
                                    -      s_u[globalIdx    ]
                                    )/spr3;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk_ko_deriv42_z(DEVICE_REAL * const  Du, const DEVICE_REAL * const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL dx              =  blk->m_dx[2];
        const DEVICE_REAL idx             =  1.0 / dx;
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT si = i ;//+ pw;  
        const DEVICE_INT sj = j ;//+ pw;
        
        const DEVICE_INT ss = nx * ny;

        const DEVICE_INT idx_base = IDX_B(si, sj, 0);
        const DEVICE_INT sk_inc   = nx*ny; 

        //if(si < (actual_nx-pw)  &&  sj < (actual_ny-pw))
        for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
        {
            const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

           Du[globalIdx] = (    -      s_u[globalIdx - 3 * ss]
                                +  6.0*s_u[globalIdx - 2 * ss]
                                - 15.0*s_u[globalIdx - 1 * ss]
                                + 20.0*s_u[globalIdx         ]
                                - 15.0*s_u[globalIdx + 1 * ss]
                                +  6.0*s_u[globalIdx + 2 * ss]
                                -      s_u[globalIdx + 3 * ss]
                        ) * pre_factor_6_dx;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {

            GPUDevice::sync_threads();
            //for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            DEVICE_INT sk = pw;
            if(sk == pw && si < (actual_nx-pw)  &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;

                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (          s_u[globalIdx +3 * ss]
                                        - 3.0*s_u[globalIdx +2 * ss]
                                        + 3.0*s_u[globalIdx +1 * ss]
                                        -     s_u[globalIdx ]
                                        )/smr3;

                Du[globalIdx + 1 * ss] =   (       s_u[globalIdx +  4 * ss]
                                            -  6.0*s_u[globalIdx +  3 * ss]
                                            + 12.0*s_u[globalIdx +  2 * ss]
                                            - 10.0*s_u[globalIdx +  1 * ss]
                                            +  3.0*s_u[globalIdx     ]
                                            )/smr2;

                Du[globalIdx + 2 *ss ] =   (       s_u[globalIdx + 5 * ss]
                                            -  6.0*s_u[globalIdx + 4 * ss]
                                            + 15.0*s_u[globalIdx + 3 * ss]
                                            - 19.0*s_u[globalIdx + 2 * ss]
                                            + 12.0*s_u[globalIdx + 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/smr1;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_INT ie = actual_nz-pw;
            //for(DEVICE_INT sk = pw; sk < (actual_nz - pw); sk++)
            DEVICE_INT sk = ie-1;
            if(sk == ie-1 && si < (actual_nx-pw)  &&  sj < (actual_ny-pw))
            {
                const DEVICE_INT globalIdx = idx_base + sk * sk_inc;
                
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2 * ss] =    (        s_u[globalIdx - 5 * ss]
                                            -  6.0*s_u[globalIdx - 4 * ss]
                                            + 15.0*s_u[globalIdx - 3 * ss]
                                            - 19.0*s_u[globalIdx - 2 * ss]
                                            + 12.0*s_u[globalIdx - 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/spr1;

                Du[globalIdx -1 *ss] =    (        s_u[globalIdx - 4 * ss]
                                            -  6.0*s_u[globalIdx - 3 * ss]
                                            + 12.0*s_u[globalIdx - 2 * ss]
                                            - 10.0*s_u[globalIdx - 1 * ss]
                                            +  3.0*s_u[globalIdx    ]
                                            )/spr2;

                Du[globalIdx] =    (       s_u[globalIdx - 3 * ss]
                                    -  3.0*s_u[globalIdx - 2 * ss]
                                    +  3.0*s_u[globalIdx - 1 * ss]
                                    -      s_u[globalIdx    ]
                                    )/spr3;
            }

        }
        
        return ;
    }


    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_deriv644_x(DEVICE_REAL * __restrict__ const  Du, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL idx = 1.0 / blk->m_dx[0];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();

        DEVICE_UINT inx     = GPUDevice::block_dim_x();

        for(DEVICE_INT si = i+pw; si < (actual_nx-pw); si+=inx)
        for(DEVICE_INT sk = k; sk < actual_nz; sk+=inx)
        for(DEVICE_INT sj = j; sj < actual_ny; sj+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
            Du[globalIdx] = (    -   s_u[globalIdx - 3] 
                            +  9.0 * s_u[globalIdx - 2]
                            - 45.0 * s_u[globalIdx - 1]
                            + 45.0 * s_u[globalIdx + 1]
                            -  9.0 * s_u[globalIdx + 2]
                            +        s_u[globalIdx + 3]) * idx_by_60;

        }
        
        
        const DEVICE_UINT bflag = blk->m_bflag;
        if (bflag & (1u<<OCT_DIR_LEFT)) {

            GPUDevice::sync_threads();
            
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT si = i + pw;

            if(si == pw)
            for(DEVICE_INT sk = k; sk < actual_nz; sk+=inx)
            for(DEVICE_INT sj = j; sj < actual_ny; sj+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                Du[globalIdx] = (  - 25.0  * s_u[globalIdx    ]
                                    + 48.0 * s_u[globalIdx + 1]
                                    - 36.0 * s_u[globalIdx + 2]
                                    + 16.0 * s_u[globalIdx + 3]
                                    -  3.0 * s_u[globalIdx + 4]
                                ) * idx_by_12;

                Du[globalIdx +1]  =  (  -  3.0 * s_u[globalIdx    ]
                                        - 10.0 * s_u[globalIdx + 1]
                                        + 18.0 * s_u[globalIdx + 2]
                                        -  6.0 * s_u[globalIdx + 3]
                                        +        s_u[globalIdx + 4]
                                    ) * idx_by_12;

                Du[globalIdx+2]   =  (   +      s_u[globalIdx    ]
                                        - 8.0 * s_u[globalIdx + 1]
                                        + 8.0 * s_u[globalIdx + 3]
                                        -       s_u[globalIdx + 4]
                                    ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT ie = actual_nx-pw;
            const DEVICE_INT si = i + pw;

            if(si == ie-1)
            for(DEVICE_INT sk = k; sk < actual_nz; sk+=inx)
            for(DEVICE_INT sj = j; sj < actual_ny; sj+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                Du[globalIdx-2]   = (  +        s_u[globalIdx - 4]
                                        - 8.0 * s_u[globalIdx - 3]
                                        + 8.0 * s_u[globalIdx - 1]
                                        -       s_u[globalIdx    ]
                                        ) * idx_by_12;

                Du[globalIdx-1] =  (  -          s_u[globalIdx - 4]
                                        +  6.0 * s_u[globalIdx - 3]
                                        - 18.0 * s_u[globalIdx - 2]
                                        + 10.0 * s_u[globalIdx - 1]
                                        +  3.0 * s_u[globalIdx    ]
                                        ) * idx_by_12;
                                    
                Du[globalIdx] =    (       3.0 * s_u[globalIdx - 4]
                                        - 16.0 * s_u[globalIdx - 3]
                                        + 36.0 * s_u[globalIdx - 2]
                                        - 48.0 * s_u[globalIdx - 1]
                                        + 25.0 * s_u[globalIdx    ]
                                        ) * idx_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_deriv644_y(DEVICE_REAL * __restrict__ const  Du, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL idx = 1.0 / blk->m_dx[1];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();

        DEVICE_UINT inx     = GPUDevice::block_dim_x();

        const DEVICE_INT ss = nx;

        for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
        for(DEVICE_INT sk = k   ; sk < actual_nz;    sk+=inx)
        for(DEVICE_INT si = i   ; si < actual_nx;    si+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
            Du[globalIdx] = (    -   s_u[globalIdx - 3 * ss ] 
                            +  9.0 * s_u[globalIdx - 2 * ss ]
                            - 45.0 * s_u[globalIdx - 1 * ss ]
                            + 45.0 * s_u[globalIdx + 1 * ss ]
                            -  9.0 * s_u[globalIdx + 2 * ss ]
                            +        s_u[globalIdx + 3 * ss ]) * idx_by_60;

        }

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_DOWN)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT sj = j+pw; 

            if(sj == pw)
            for(DEVICE_INT sk = k;    sk < actual_nz; sk+=inx)
            for(DEVICE_INT si = i;    si < actual_nx; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                Du[globalIdx          ] = (  - 25.0  * s_u[globalIdx       ]
                                            + 48.0 * s_u[globalIdx + 1 * ss]
                                            - 36.0 * s_u[globalIdx + 2 * ss]
                                            + 16.0 * s_u[globalIdx + 3 * ss]
                                            -  3.0 * s_u[globalIdx + 4 * ss]
                                        ) * idx_by_12;

                Du[globalIdx + 1 * ss]  =  (    -  3.0 * s_u[globalIdx         ]
                                                - 10.0 * s_u[globalIdx + 1 * ss]
                                                + 18.0 * s_u[globalIdx + 2 * ss]
                                                -  6.0 * s_u[globalIdx + 3 * ss]
                                                +        s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;

                Du[globalIdx + 2 * ss]   =  (    +      s_u[globalIdx         ]
                                                - 8.0 * s_u[globalIdx + 1 * ss]
                                                + 8.0 * s_u[globalIdx + 3 * ss]
                                                -       s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT ie = actual_ny-pw;
            const DEVICE_INT sj = j+pw; 
            
            if(sj == ie-1)
            for(DEVICE_INT sk = k;    sk < actual_nz; sk+=inx)
            for(DEVICE_INT si = i;    si < actual_nx; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                Du[globalIdx-2 * ss]   = (  +       s_u[globalIdx - 4 * ss]
                                            - 8.0 * s_u[globalIdx - 3 * ss]
                                            + 8.0 * s_u[globalIdx - 1 * ss]
                                            -       s_u[globalIdx         ]
                                            ) * idx_by_12;

                Du[globalIdx-1 * ss] =  (   -        s_u[globalIdx - 4 * ss]
                                            +  6.0 * s_u[globalIdx - 3 * ss]
                                            - 18.0 * s_u[globalIdx - 2 * ss]
                                            + 10.0 * s_u[globalIdx - 1 * ss]
                                            +  3.0 * s_u[globalIdx         ]
                                            ) * idx_by_12;
                                    
                Du[globalIdx] =    (       3.0 * s_u[globalIdx - 4 * ss]
                                        - 16.0 * s_u[globalIdx - 3 * ss]
                                        + 36.0 * s_u[globalIdx - 2 * ss]
                                        - 48.0 * s_u[globalIdx - 1 * ss]
                                        + 25.0 * s_u[globalIdx         ]
                                        ) * idx_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_deriv644_z(DEVICE_REAL * __restrict__ const  Du, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL idx = 1.0 / blk->m_dx[2];
        const DEVICE_REAL idx_by_60 = idx / 60.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();
        
        DEVICE_UINT inx     = GPUDevice::block_dim_x();

        const DEVICE_INT ss = nx * ny;

        for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
        for(DEVICE_INT sj = j   ; sj < actual_ny   ; sj+=inx)
        for(DEVICE_INT si = i   ; si < actual_nx   ; si+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
            Du[globalIdx] = (    -   s_u[globalIdx - 3 * ss ] 
                            +  9.0 * s_u[globalIdx - 2 * ss ]
                            - 45.0 * s_u[globalIdx - 1 * ss ]
                            + 45.0 * s_u[globalIdx + 1 * ss ]
                            -  9.0 * s_u[globalIdx + 2 * ss ]
                            +        s_u[globalIdx + 3 * ss ]) * idx_by_60;

        }


        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT sk = k+pw;

            if(sk == pw)
            for(DEVICE_INT sj = j; sj < actual_ny; sj+=inx)
            for(DEVICE_INT si = i; si < actual_nx; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                Du[globalIdx          ] = (  - 25.0  * s_u[globalIdx       ]
                                            + 48.0 * s_u[globalIdx + 1 * ss]
                                            - 36.0 * s_u[globalIdx + 2 * ss]
                                            + 16.0 * s_u[globalIdx + 3 * ss]
                                            -  3.0 * s_u[globalIdx + 4 * ss]
                                        ) * idx_by_12;

                Du[globalIdx + 1 * ss]  =  (    -  3.0 * s_u[globalIdx         ]
                                                - 10.0 * s_u[globalIdx + 1 * ss]
                                                + 18.0 * s_u[globalIdx + 2 * ss]
                                                -  6.0 * s_u[globalIdx + 3 * ss]
                                                +        s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;

                Du[globalIdx + 2 * ss]   =  (    +      s_u[globalIdx         ]
                                                - 8.0 * s_u[globalIdx + 1 * ss]
                                                + 8.0 * s_u[globalIdx + 3 * ss]
                                                -       s_u[globalIdx + 4 * ss]
                                            ) * idx_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_by_12 = idx / 12.0;
            const DEVICE_INT ie = actual_nz-pw;
            const DEVICE_INT sk = k+pw;

            if(sk == ie-1)
            for(DEVICE_INT sj = j; sj < actual_ny; sj+=inx)
            for(DEVICE_INT si = i; si < actual_nx; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                Du[globalIdx-2 * ss]   = (  +       s_u[globalIdx - 4 * ss]
                                            - 8.0 * s_u[globalIdx - 3 * ss]
                                            + 8.0 * s_u[globalIdx - 1 * ss]
                                            -       s_u[globalIdx         ]
                                            ) * idx_by_12;

                Du[globalIdx-1 * ss] =  (   -        s_u[globalIdx - 4 * ss]
                                            +  6.0 * s_u[globalIdx - 3 * ss]
                                            - 18.0 * s_u[globalIdx - 2 * ss]
                                            + 10.0 * s_u[globalIdx - 1 * ss]
                                            +  3.0 * s_u[globalIdx         ]
                                            ) * idx_by_12;
                                    
                Du[globalIdx] =    (       3.0 * s_u[globalIdx - 4 * ss]
                                        - 16.0 * s_u[globalIdx - 3 * ss]
                                        + 36.0 * s_u[globalIdx - 2 * ss]
                                        - 48.0 * s_u[globalIdx - 1 * ss]
                                        + 25.0 * s_u[globalIdx         ]
                                        ) * idx_by_12;
            }

        }
        
        return ;
    }


    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_deriv644_xx(DEVICE_REAL * __restrict__ const  DDu, const DEVICE_REAL *  __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        const DEVICE_REAL dx              =  blk->m_dx[0];
        const DEVICE_REAL idx_sqrd = 1.0 / (dx * dx);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();
        
        DEVICE_UINT inx     = GPUDevice::block_dim_x();
        const DEVICE_INT sk = k + pw;
        const DEVICE_INT sj = j + pw;

        for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
        // for(DEVICE_INT sk = k   ; sk < actual_nz   ; sk+=inx)
        // for(DEVICE_INT sj = j   ; sj < actual_ny   ; sj+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

            DDu[globalIdx] = (      2.0 * s_u[globalIdx - 3]
                                -  27.0 * s_u[globalIdx - 2]
                                + 270.0 * s_u[globalIdx - 1]
                                - 490.0 * s_u[globalIdx    ]
                                + 270.0 * s_u[globalIdx + 1]
                                -  27.0 * s_u[globalIdx + 2]
                                +   2.0 * s_u[globalIdx + 3]
                                ) * idx_sqrd_by_180;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_LEFT)) {

            GPUDevice::sync_threads();
            
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT si = i + pw;
            
            if(si == pw)
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                DDu[globalIdx]      = (       45.0  * s_u[globalIdx    ]
                                            - 154.0 * s_u[globalIdx + 1]
                                            + 214.0 * s_u[globalIdx + 2]
                                            - 156.0 * s_u[globalIdx + 3]
                                            +  61.0 * s_u[globalIdx + 4]
                                            -  10.0 * s_u[globalIdx + 5]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx+1]    = (        10.0 * s_u[globalIdx    ]
                                            -  15.0 * s_u[globalIdx + 1]
                                            -   4.0 * s_u[globalIdx + 2]
                                            +  14.0 * s_u[globalIdx + 3]
                                            -   6.0 * s_u[globalIdx + 4]
                                            +         s_u[globalIdx + 5]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx+2]    = (     -         s_u[globalIdx    ]
                                            +  16.0 * s_u[globalIdx + 1]
                                            -  30.0 * s_u[globalIdx + 2]
                                            +  16.0 * s_u[globalIdx + 3]
                                            -         s_u[globalIdx + 4]
                                            ) * idx_sqrd_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT ie = actual_nx-pw;
            const DEVICE_INT si = i + pw;
            
            if(si == ie-1)
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                DDu[globalIdx-2]    = (     -         s_u[globalIdx -4]
                                            +  16.0 * s_u[globalIdx -3]
                                            -  30.0 * s_u[globalIdx -2]
                                            +  16.0 * s_u[globalIdx -1]
                                            -         s_u[globalIdx   ]
                                            ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx-1]    = (               s_u[globalIdx - 5]
                                            -   6.0 * s_u[globalIdx - 4]
                                            +  14.0 * s_u[globalIdx - 3]
                                            -   4.0 * s_u[globalIdx - 2]
                                            -  15.0 * s_u[globalIdx - 1]
                                            +  10.0 * s_u[globalIdx    ]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx]      = (     -  10.0 * s_u[globalIdx - 5] 
                                            +  61.0 * s_u[globalIdx - 4] 
                                            - 156.0 * s_u[globalIdx - 3] 
                                            + 214.0 * s_u[globalIdx - 2] 
                                            - 154.0 * s_u[globalIdx - 1] 
                                            +  45.0 * s_u[globalIdx    ] 
                                            ) * idx_sqrd_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_deriv644_yy(DEVICE_REAL *  __restrict__ const  DDu, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        const DEVICE_REAL dx              =  blk->m_dx[1];
        const DEVICE_REAL idx_sqrd = 1.0 / (dx * dx);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT ss = nx;

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();
        
        DEVICE_UINT inx     = GPUDevice::block_dim_x();

        const DEVICE_INT sk = k + pw;
        const DEVICE_INT si = i + pw;

        for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
        //for(DEVICE_INT sk = k   ; sk < actual_nz   ; sk+=inx)
        //for(DEVICE_INT si = i   ; si < actual_nx   ; si+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

            DDu[globalIdx] = (      2.0 * s_u[globalIdx - 3 * ss]
                                -  27.0 * s_u[globalIdx - 2 * ss]
                                + 270.0 * s_u[globalIdx - 1 * ss]
                                - 490.0 * s_u[globalIdx         ]
                                + 270.0 * s_u[globalIdx + 1 * ss]
                                -  27.0 * s_u[globalIdx + 2 * ss]
                                +   2.0 * s_u[globalIdx + 3 * ss]
                                ) * idx_sqrd_by_180;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_DOWN)) {

            GPUDevice::sync_threads();
            
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT sj = j+pw;
            
            if(sj == pw)
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                DDu[globalIdx]      = (       45.0  * s_u[globalIdx    ]
                                            - 154.0 * s_u[globalIdx + 1 * ss]
                                            + 214.0 * s_u[globalIdx + 2 * ss]
                                            - 156.0 * s_u[globalIdx + 3 * ss]
                                            +  61.0 * s_u[globalIdx + 4 * ss]
                                            -  10.0 * s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx + 1 * ss]    = (        10.0 * s_u[globalIdx    ]
                                            -  15.0 * s_u[globalIdx + 1 * ss]
                                            -   4.0 * s_u[globalIdx + 2 * ss]
                                            +  14.0 * s_u[globalIdx + 3 * ss]
                                            -   6.0 * s_u[globalIdx + 4 * ss]
                                            +         s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx + 2 * ss]    = (     -         s_u[globalIdx    ]
                                            +  16.0 * s_u[globalIdx + 1 * ss]
                                            -  30.0 * s_u[globalIdx + 2 * ss]
                                            +  16.0 * s_u[globalIdx + 3 * ss]
                                            -         s_u[globalIdx + 4 * ss]
                                            ) * idx_sqrd_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT ie = actual_ny-pw;
            const DEVICE_INT sj = j+pw;
            
            if(sj == ie-1)
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                DDu[globalIdx - 2 * ss ]    = (     -         s_u[globalIdx -4 * ss]
                                                    +  16.0 * s_u[globalIdx -3 * ss]
                                                    -  30.0 * s_u[globalIdx -2 * ss]
                                                    +  16.0 * s_u[globalIdx -1 * ss]
                                                    -         s_u[globalIdx   ]
                                                    ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx -1 * ss]    =   (               s_u[globalIdx - 5 * ss]
                                                    -   6.0 * s_u[globalIdx - 4 * ss]
                                                    +  14.0 * s_u[globalIdx - 3 * ss]
                                                    -   4.0 * s_u[globalIdx - 2 * ss]
                                                    -  15.0 * s_u[globalIdx - 1 * ss]
                                                    +  10.0 * s_u[globalIdx    ]
                                                    ) * idx_sqrd_by_12;
                
                DDu[globalIdx]      = (     -  10.0 * s_u[globalIdx - 5 * ss] 
                                            +  61.0 * s_u[globalIdx - 4 * ss] 
                                            - 156.0 * s_u[globalIdx - 3 * ss] 
                                            + 214.0 * s_u[globalIdx - 2 * ss] 
                                            - 154.0 * s_u[globalIdx - 1 * ss] 
                                            +  45.0 * s_u[globalIdx    ] 
                                            ) * idx_sqrd_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_deriv644_zz(DEVICE_REAL * __restrict__ const  DDu, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        const DEVICE_REAL dx              =  blk->m_dx[2];
        const DEVICE_REAL idx_sqrd = 1.0 / (dx * dx);
        const DEVICE_REAL idx_sqrd_by_180 = idx_sqrd / 180.0;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT ss = nx * ny;

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();
        
        DEVICE_UINT inx     = GPUDevice::block_dim_x();
        const DEVICE_INT sj = j + pw;
        const DEVICE_INT si = i + pw;

        for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
        //for(DEVICE_INT sj = j   ; sj < actual_ny   ; sj+=inx)
        //for(DEVICE_INT si = i   ; si < actual_nx   ; si+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

            DDu[globalIdx] = (      2.0 * s_u[globalIdx - 3 * ss]
                                -  27.0 * s_u[globalIdx - 2 * ss]
                                + 270.0 * s_u[globalIdx - 1 * ss]
                                - 490.0 * s_u[globalIdx         ]
                                + 270.0 * s_u[globalIdx + 1 * ss]
                                -  27.0 * s_u[globalIdx + 2 * ss]
                                +   2.0 * s_u[globalIdx + 3 * ss]
                                ) * idx_sqrd_by_180;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {

            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            
            const DEVICE_INT sk = k + pw;
            if(sk == pw)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                DDu[globalIdx]      = (       45.0  * s_u[globalIdx    ]
                                            - 154.0 * s_u[globalIdx + 1 * ss]
                                            + 214.0 * s_u[globalIdx + 2 * ss]
                                            - 156.0 * s_u[globalIdx + 3 * ss]
                                            +  61.0 * s_u[globalIdx + 4 * ss]
                                            -  10.0 * s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
            
                DDu[globalIdx + 1 * ss]    = (        10.0 * s_u[globalIdx    ]
                                            -  15.0 * s_u[globalIdx + 1 * ss]
                                            -   4.0 * s_u[globalIdx + 2 * ss]
                                            +  14.0 * s_u[globalIdx + 3 * ss]
                                            -   6.0 * s_u[globalIdx + 4 * ss]
                                            +         s_u[globalIdx + 5 * ss]
                                            ) * idx_sqrd_by_12;
                
                DDu[globalIdx + 2 * ss]    = (     -         s_u[globalIdx    ]
                                            +  16.0 * s_u[globalIdx + 1 * ss]
                                            -  30.0 * s_u[globalIdx + 2 * ss]
                                            +  16.0 * s_u[globalIdx + 3 * ss]
                                            -         s_u[globalIdx + 4 * ss]
                                            ) * idx_sqrd_by_12;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_REAL idx_sqrd_by_12 = idx_sqrd / 12.0;
            const DEVICE_INT ie = actual_nz-pw;
            const DEVICE_INT sk = k + pw;
            if(sk == ie-1)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                DDu[globalIdx - 2 * ss ]    = (     -         s_u[globalIdx -4 * ss]
                                                    +  16.0 * s_u[globalIdx -3 * ss]
                                                    -  30.0 * s_u[globalIdx -2 * ss]
                                                    +  16.0 * s_u[globalIdx -1 * ss]
                                                    -         s_u[globalIdx   ]
                                                    ) * idx_sqrd_by_12;
            
                
                DDu[globalIdx -1 * ss]    =   (               s_u[globalIdx - 5 * ss]
                                                    -   6.0 * s_u[globalIdx - 4 * ss]
                                                    +  14.0 * s_u[globalIdx - 3 * ss]
                                                    -   4.0 * s_u[globalIdx - 2 * ss]
                                                    -  15.0 * s_u[globalIdx - 1 * ss]
                                                    +  10.0 * s_u[globalIdx    ]
                                                    ) * idx_sqrd_by_12;
                
                DDu[globalIdx         ]    =  (     -  10.0 * s_u[globalIdx - 5 * ss] 
                                                    +  61.0 * s_u[globalIdx - 4 * ss] 
                                                    - 156.0 * s_u[globalIdx - 3 * ss] 
                                                    + 214.0 * s_u[globalIdx - 2 * ss] 
                                                    - 154.0 * s_u[globalIdx - 1 * ss] 
                                                    +  45.0 * s_u[globalIdx    ] 
                                                    ) * idx_sqrd_by_12;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_ko_deriv42_x(DEVICE_REAL * __restrict__ const  Du, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL dx              =  blk->m_dx[0];
        const DEVICE_REAL idx             =  1.0 / dx;
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();
        
        const DEVICE_UINT inx  = GPUDevice::block_dim_x();

        const DEVICE_INT sk = k+pw;
        const DEVICE_INT sj = j+pw;

        // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
        // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
        for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

            Du[globalIdx] = (    -         s_u[globalIdx - 3]
                                    +  6.0*s_u[globalIdx - 2]
                                    - 15.0*s_u[globalIdx - 1]
                                    + 20.0*s_u[globalIdx    ]
                                    - 15.0*s_u[globalIdx + 1]
                                    +  6.0*s_u[globalIdx + 2]
                                    -      s_u[globalIdx + 3]
                            ) * pre_factor_6_dx;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_LEFT)) {

            GPUDevice::sync_threads();
            
            const DEVICE_INT si = i + pw;
            
            if(si == pw)
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (      s_u[globalIdx +3]
                                    - 3.0*s_u[globalIdx +2]
                                    + 3.0*s_u[globalIdx +1]
                                    -     s_u[globalIdx ]
                                    )/smr3;
                Du[globalIdx+1] =  (       s_u[globalIdx +  4]
                                    -  6.0*s_u[globalIdx +  3]
                                    + 12.0*s_u[globalIdx +  2]
                                    - 10.0*s_u[globalIdx +  1]
                                    +  3.0*s_u[globalIdx     ]
                                    )/smr2;

                Du[globalIdx+2] =  (       s_u[globalIdx + 5]
                                    -  6.0*s_u[globalIdx + 4]
                                    + 15.0*s_u[globalIdx + 3]
                                    - 19.0*s_u[globalIdx + 2]
                                    + 12.0*s_u[globalIdx + 1]
                                    -  3.0*s_u[globalIdx    ]
                                    )/smr1;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_RIGHT)) {
            GPUDevice::sync_threads();
            const DEVICE_INT ie = actual_nx-pw;
            const DEVICE_INT si = i + pw;
            
            if(si == ie-1)
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2] = (        s_u[globalIdx - 5]
                                    -  6.0*s_u[globalIdx - 4]
                                    + 15.0*s_u[globalIdx - 3]
                                    - 19.0*s_u[globalIdx - 2]
                                    + 12.0*s_u[globalIdx - 1]
                                    -  3.0*s_u[globalIdx    ]
                                    )/spr1;

                Du[globalIdx-1] = (        s_u[globalIdx - 4]
                                    -  6.0*s_u[globalIdx - 3]
                                    + 12.0*s_u[globalIdx - 2]
                                    - 10.0*s_u[globalIdx - 1]
                                    +  3.0*s_u[globalIdx    ]
                                    )/spr2;

                Du[globalIdx] =    (       s_u[globalIdx - 3]
                                    -  3.0*s_u[globalIdx - 2]
                                    +  3.0*s_u[globalIdx - 1]
                                    -      s_u[globalIdx    ]
                                    )/spr3;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_ko_deriv42_y(DEVICE_REAL * __restrict__ const  Du, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL dx              =  blk->m_dx[1];
        const DEVICE_REAL idx             =  1.0 / dx;
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT ss = nx;

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();
        
        const DEVICE_UINT inx     = GPUDevice::block_dim_x();
        const DEVICE_INT sk = k+pw;
        const DEVICE_INT si = i+pw;

        for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
        // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
        // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

            Du[globalIdx] = (    -         s_u[globalIdx - 3 * ss]
                                    +  6.0*s_u[globalIdx - 2 * ss]
                                    - 15.0*s_u[globalIdx - 1 * ss]
                                    + 20.0*s_u[globalIdx         ]
                                    - 15.0*s_u[globalIdx + 1 * ss]
                                    +  6.0*s_u[globalIdx + 2 * ss]
                                    -      s_u[globalIdx + 3 * ss]
                            ) * pre_factor_6_dx;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_DOWN)) {

            GPUDevice::sync_threads();
            const DEVICE_INT sj = j + pw;
            if(sj == pw )
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (          s_u[globalIdx +3 * ss]
                                        - 3.0*s_u[globalIdx +2 * ss]
                                        + 3.0*s_u[globalIdx +1 * ss]
                                        -     s_u[globalIdx ]
                                        )/smr3;

                Du[globalIdx + 1 * ss] =   (       s_u[globalIdx +  4 * ss]
                                            -  6.0*s_u[globalIdx +  3 * ss]
                                            + 12.0*s_u[globalIdx +  2 * ss]
                                            - 10.0*s_u[globalIdx +  1 * ss]
                                            +  3.0*s_u[globalIdx     ]
                                            )/smr2;

                Du[globalIdx + 2 *ss ] =   (       s_u[globalIdx + 5 * ss]
                                            -  6.0*s_u[globalIdx + 4 * ss]
                                            + 15.0*s_u[globalIdx + 3 * ss]
                                            - 19.0*s_u[globalIdx + 2 * ss]
                                            + 12.0*s_u[globalIdx + 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/smr1;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_UP)) {
            GPUDevice::sync_threads();
            const DEVICE_INT ie = actual_ny-pw;
            const DEVICE_INT sj = j + pw;
            if(sj == ie-1 )
            // for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2 * ss] =    (        s_u[globalIdx - 5 * ss]
                                            -  6.0*s_u[globalIdx - 4 * ss]
                                            + 15.0*s_u[globalIdx - 3 * ss]
                                            - 19.0*s_u[globalIdx - 2 * ss]
                                            + 12.0*s_u[globalIdx - 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/spr1;

                Du[globalIdx -1 *ss] =    (        s_u[globalIdx - 4 * ss]
                                            -  6.0*s_u[globalIdx - 3 * ss]
                                            + 12.0*s_u[globalIdx - 2 * ss]
                                            - 10.0*s_u[globalIdx - 1 * ss]
                                            +  3.0*s_u[globalIdx    ]
                                            )/spr2;

                Du[globalIdx] =    (       s_u[globalIdx - 3 * ss]
                                    -  3.0*s_u[globalIdx - 2 * ss]
                                    +  3.0*s_u[globalIdx - 1 * ss]
                                    -      s_u[globalIdx    ]
                                    )/spr3;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    DEVICE_FUNC __inline__ void __blk1_ko_deriv42_z(DEVICE_REAL * __restrict__ const  Du, const DEVICE_REAL * __restrict__ const s_u, const BlockGPU3D* const blk)
    {
        
        const DEVICE_REAL dx              =  blk->m_dx[2];
        const DEVICE_REAL idx             =  1.0 / dx;
        const DEVICE_REAL pre_factor_6_dx = -1.0 / 64.0 / dx;

        const DEVICE_INT nx   = pencil_sz;//blk->m_aligned_sz[0];
        const DEVICE_INT ny   = pencil_sz;//blk->m_aligned_sz[1];
        const DEVICE_INT nz   = pencil_sz;//blk->m_aligned_sz[2];

        const DEVICE_INT actual_nx = 4*pw + 1; //blk->m_sz[0];
        const DEVICE_INT actual_ny = 4*pw + 1; //blk->m_sz[1];
        const DEVICE_INT actual_nz = 4*pw + 1; //blk->m_sz[2];

        const DEVICE_INT ss = nx * ny;

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::thread_id_z();
        
        const DEVICE_UINT inx     = GPUDevice::block_dim_x();
        const DEVICE_INT sj = j+pw;
        const DEVICE_INT si = i+pw;

        for(DEVICE_INT sk = k+pw; sk < actual_nz-pw; sk+=inx)
        // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
        // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
        {
            const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

            Du[globalIdx] = (    -     s_u[globalIdx - 3 * ss]
                                +  6.0*s_u[globalIdx - 2 * ss]
                                - 15.0*s_u[globalIdx - 1 * ss]
                                + 20.0*s_u[globalIdx         ]
                                - 15.0*s_u[globalIdx + 1 * ss]
                                +  6.0*s_u[globalIdx + 2 * ss]
                                -      s_u[globalIdx + 3 * ss]
                        ) * pre_factor_6_dx;

        }
            

        const DEVICE_UINT bflag = blk->m_bflag;

        if (bflag & (1u<<OCT_DIR_BACK)) {

            GPUDevice::sync_threads();
            DEVICE_INT sk = k + pw;
            if(sk == pw)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);

                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx]   =  (          s_u[globalIdx +3 * ss]
                                        - 3.0*s_u[globalIdx +2 * ss]
                                        + 3.0*s_u[globalIdx +1 * ss]
                                        -     s_u[globalIdx ]
                                        )/smr3;

                Du[globalIdx + 1 * ss] =   (       s_u[globalIdx +  4 * ss]
                                            -  6.0*s_u[globalIdx +  3 * ss]
                                            + 12.0*s_u[globalIdx +  2 * ss]
                                            - 10.0*s_u[globalIdx +  1 * ss]
                                            +  3.0*s_u[globalIdx     ]
                                            )/smr2;

                Du[globalIdx + 2 *ss ] =   (       s_u[globalIdx + 5 * ss]
                                            -  6.0*s_u[globalIdx + 4 * ss]
                                            + 15.0*s_u[globalIdx + 3 * ss]
                                            - 19.0*s_u[globalIdx + 2 * ss]
                                            + 12.0*s_u[globalIdx + 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/smr1;
            }
                            
            
        }

        if (bflag & (1u<<OCT_DIR_FRONT)) {
            GPUDevice::sync_threads();
            const DEVICE_INT ie = actual_nz-pw;
            DEVICE_INT sk = k + pw;
            if(sk == ie-1)
            // for(DEVICE_INT sj = j+pw; sj < actual_ny-pw; sj+=inx)
            // for(DEVICE_INT si = i+pw; si < actual_nx-pw; si+=inx)
            {
                const DEVICE_INT globalIdx = IDX_B(si,sj,sk);
                
                const DEVICE_REAL  smr3=59.0/48.0*64*dx;
                const DEVICE_REAL  smr2=43.0/48.0*64*dx;
                const DEVICE_REAL  smr1=49.0/48.0*64*dx;
                const DEVICE_REAL  spr3=smr3;
                const DEVICE_REAL  spr2=smr2;
                const DEVICE_REAL  spr1=smr1;

                Du[globalIdx-2 * ss] =    (        s_u[globalIdx - 5 * ss]
                                            -  6.0*s_u[globalIdx - 4 * ss]
                                            + 15.0*s_u[globalIdx - 3 * ss]
                                            - 19.0*s_u[globalIdx - 2 * ss]
                                            + 12.0*s_u[globalIdx - 1 * ss]
                                            -  3.0*s_u[globalIdx    ]
                                            )/spr1;

                Du[globalIdx -1 *ss] =    (        s_u[globalIdx - 4 * ss]
                                            -  6.0*s_u[globalIdx - 3 * ss]
                                            + 12.0*s_u[globalIdx - 2 * ss]
                                            - 10.0*s_u[globalIdx - 1 * ss]
                                            +  3.0*s_u[globalIdx    ]
                                            )/spr2;

                Du[globalIdx] =    (       s_u[globalIdx - 3 * ss]
                                    -  3.0*s_u[globalIdx - 2 * ss]
                                    +  3.0*s_u[globalIdx - 1 * ss]
                                    -      s_u[globalIdx    ]
                                    )/spr3;
            }

        }
        
        return ;
    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_deriv_x(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __deriv644_x<pw,pencils,pencil_sz>(Du,u,blk);}
    
    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_deriv_y(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __deriv644_y<pw,pencils,pencil_sz>(Du,u,blk);}
    
    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_deriv_z(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __deriv644_z<pw,pencils,pencil_sz>(Du,u,blk);}

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_deriv_xx(DEVICE_REAL * const  DDu, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __deriv644_xx<pw,pencils,pencil_sz>(DDu,u,blk);}
 
    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_deriv_yy(DEVICE_REAL * const  DDu, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __deriv644_yy<pw,pencils,pencil_sz>(DDu,u,blk);}
 
    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_deriv_zz(DEVICE_REAL * const  DDu, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __deriv644_zz<pw,pencils,pencil_sz>(DDu,u,blk);}
    
    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_ko_deriv_x(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __ko_deriv42_x<pw,pencils,pencil_sz>(Du,u,blk);}

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_ko_deriv_y(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __ko_deriv42_y<pw,pencils,pencil_sz>(Du,u,blk);}

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_ko_deriv_z(DEVICE_REAL * const  Du, const DEVICE_REAL * const  u, const BlockGPU3D* const blk) {return __ko_deriv42_z<pw,pencils,pencil_sz>(Du,u,blk);}


    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_blk_d(DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {
        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT n   = nx*ny*nz;

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        SHARED_MEM DEVICE_REAL s_u[pencil_sz* pencil_sz * pencil_sz];

        for(DEVICE_INT k = 0; k <  actual_nz; k++)
        {
            const DEVICE_INT globalIdx = IDX_B(i,j,k);
            s_u[globalIdx]              = u[globalIdx];
        }

        GPUDevice::sync_threads();
        
        __blk_deriv644_x<pw,pencils,pencil_sz>(Du ,     s_u,blk);
        __blk_deriv644_y<pw,pencils,pencil_sz>(Du + n  , s_u,blk);
        __blk_deriv644_z<pw,pencils,pencil_sz>(Du + 2*n, s_u,blk);
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_blk_dd(DEVICE_REAL * const  DDu, DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {
        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT n   = nx*ny*nz;

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        SHARED_MEM DEVICE_REAL s_u[pencil_sz* pencil_sz * pencil_sz];

        GPUDevice::sync_threads();
        for(DEVICE_INT k = 0; k <  actual_nz; k++)
        {
            const DEVICE_INT globalIdx = IDX_B(i,j,k);
            s_u[globalIdx]              = u[globalIdx];
        }
        GPUDevice::sync_threads();
        __blk_deriv644_x<pw,pencils,pencil_sz> (Du + 0*n  ,s_u  ,     blk);
        __blk_deriv644_y<pw,pencils,pencil_sz> (Du + 1*n  ,s_u  ,     blk);
        

        __blk_deriv644_xx<pw,pencils,pencil_sz>(DDu        ,s_u  ,     blk);
        __blk_deriv644_yy<pw,pencils,pencil_sz>(DDu + 3*n  ,s_u  ,     blk);
        __blk_deriv644_zz<pw,pencils,pencil_sz>(DDu + 5*n  ,s_u  ,     blk);
        
        GPUDevice::sync_threads();
        for(DEVICE_INT k = 0; k <  actual_nz; k++)
        {
            const DEVICE_INT globalIdx = IDX_B(i,j,k);
            s_u[globalIdx]              = Du[0*n + globalIdx];
        }
        GPUDevice::sync_threads();

        __blk_deriv644_y<pw,pencils,pencil_sz>(DDu + 1*n , s_u  ,     blk);
        __blk_deriv644_z<pw,pencils,pencil_sz>(DDu + 2*n , s_u  ,     blk);

        GPUDevice::sync_threads();
        for(DEVICE_INT k = 0; k <  actual_nz; k++)
        {
            const DEVICE_INT globalIdx = IDX_B(i,j,k);
            s_u[globalIdx]              = Du[1*n + globalIdx];
        }

        GPUDevice::sync_threads();
        __blk_deriv644_z<pw,pencils,pencil_sz>(DDu + 4*n , s_u  ,     blk);

        return;

    }


    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_blk_kod(DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {
        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT n    = nx * ny * nz;


        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        SHARED_MEM DEVICE_REAL s_u[pencil_sz* pencil_sz * pencil_sz];

        for(DEVICE_INT k = 0; k <  actual_nz; k++)
        {
            const DEVICE_INT globalIdx = IDX_B(i,j,k);
            s_u[globalIdx]              = u[globalIdx];
        }

        GPUDevice::sync_threads();
        __blk_ko_deriv42_x<pw,pencils,pencil_sz>(Du,s_u,blk);
        __blk_ko_deriv42_y<pw,pencils,pencil_sz>(Du  + 1 * n,s_u,blk);
        __blk_ko_deriv42_z<pw,pencils,pencil_sz>(Du  + 2 * n,s_u,blk);
        
        return;

    }


    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_blk1_d(DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {
        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT n   = nx*ny*nz;

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        const DEVICE_REAL * const s_u = u;

        __blk1_deriv644_x<pw,pencils,pencil_sz>(Du ,      s_u,blk);
        __blk1_deriv644_y<pw,pencils,pencil_sz>(Du + n  , s_u,blk);
        __blk1_deriv644_z<pw,pencils,pencil_sz>(Du + 2*n, s_u,blk);
        
        return;

    }

    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_blk1_dd(DEVICE_REAL * const  DDu, DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {
        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT n   = nx*ny*nz;

        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        const DEVICE_REAL * const s_u = u;

        __blk1_deriv644_x<pw,pencils,pencil_sz> (Du + 0*n  ,s_u  ,     blk);
        __blk1_deriv644_y<pw,pencils,pencil_sz> (Du + 1*n  ,s_u  ,     blk);
        

        __blk1_deriv644_xx<pw,pencils,pencil_sz>(DDu        ,s_u  ,     blk);
        __blk1_deriv644_yy<pw,pencils,pencil_sz>(DDu + 3*n  ,s_u  ,     blk);
        __blk1_deriv644_zz<pw,pencils,pencil_sz>(DDu + 5*n  ,s_u  ,     blk);
        
        GPUDevice::sync_threads();

        __blk1_deriv644_y<pw,pencils,pencil_sz>(DDu + 1*n , &Du[0*n]  ,     blk);
        __blk1_deriv644_z<pw,pencils,pencil_sz>(DDu + 2*n , &Du[0*n]  ,     blk);

        GPUDevice::sync_threads();
        __blk1_deriv644_z<pw,pencils,pencil_sz>(DDu + 4*n , &Du[1*n]  ,     blk);

        return;

    }


    template<int pw, int pencils, int pencil_sz>
    GLOBAL_FUNC void gpu_blk1_kod(DEVICE_REAL * const  Du, const DEVICE_REAL * const u, const BlockGPU3D* const blk)
    {
        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::thread_id_y();
        
        const DEVICE_INT nx   = blk->m_aligned_sz[0];
        const DEVICE_INT ny   = blk->m_aligned_sz[1];
        const DEVICE_INT nz   = blk->m_aligned_sz[2];

        const DEVICE_INT n    = nx * ny * nz;


        const DEVICE_INT actual_nx = blk->m_sz[0];
        const DEVICE_INT actual_ny = blk->m_sz[1];
        const DEVICE_INT actual_nz = blk->m_sz[2];

        const DEVICE_REAL * const s_u = u;

        __blk1_ko_deriv42_x<pw,pencils,pencil_sz>(Du,s_u,blk);
        __blk1_ko_deriv42_y<pw,pencils,pencil_sz>(Du  + 1 * n,s_u,blk);
        __blk1_ko_deriv42_z<pw,pencils,pencil_sz>(Du  + 2 * n,s_u,blk);
        
        return;

    }




}

    
