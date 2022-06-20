/**
 * @file derivs_test.cpp
 * @brief Test cuda derivative computations. 
 * @version 0.1
 * @date 2021-11-27
 * 
 * @copyright Copyright (c) 2021
 * 
 */
#include <iostream>
#include <stdio.h>
#include <cuda_runtime.h>
#include <cuda.h>
#include "device.h"
#include "point.h"
#include "derivs_cu.cuh"
#include "derivs.h"
#include "block_gpu.h"
#include "profiler.h"

static inline void cuda_check_last_error()
{
    //#if defined(DEBUG) || defined(_DEBUG)
       cudaError_t err = cudaGetLastError();
       if (err != cudaSuccess) 
        printf("%ss\n", cudaGetErrorString(err));
    //#endif
    return ;
}

template<typename T>
T* create_block(const Point grid_min, const Point dx, const unsigned int * sz, unsigned int dof=1)
{   
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    T* blk_data = new T[nx*ny*nz*dof];

    for(unsigned int v=0; v < dof; v++)
    for(unsigned int k=0; k < nz; k++)
    for(unsigned int j=0; j < ny; j++)
    for(unsigned int i=0; i < nx; i++)
    {
        const double x = grid_min.x() + i * dx.x();
        const double y = grid_min.y() + j * dx.y();
        const double z = grid_min.z() + k * dx.z();
        blk_data[v * nx*ny*nz +  k * ny * nx + j * nx + i] = sin(2*M_PI*x) * sin(2*M_PI*y) * sin(2*M_PI*z);

    }
    
    return blk_data;

}

template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_first_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu  = new T[blkSz];
    T* dptr_u  = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du = GPUDevice::device_malloc<T>(blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(ny/pencils,nz,1);
    dim3 block_x = dim3(nx,pencils,1);

    dim3 grid_y  = dim3(nx/pencils,nz,1);
    dim3 block_y = dim3(pencils,ny,1);

    dim3 grid_z  = dim3(nx/pencils,ny,1);
    dim3 block_z = dim3(pencils,nz,1);

    t_cpu.clear();
    t_gpu.clear();

    if(dir==0)
    {
        if(ny%pencils !=0)
            std::cout<<"error : ny "<<ny<<" is not perfect factor of stencil size: "<<pencils<<std::endl;

        deriv644_x(Du_cpu,u,dx,sz,bflag);
        
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            deriv644_x(Du_cpu,u,dx,sz,bflag);
        t_cpu.stop();


        //std::cout<<"grid : "<<grid<<" block : "<<block<<std::endl;
        //printf("grid %s \n",grid);
        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);

        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }else if(dir==1)
    {

        deriv644_y(Du_cpu,u,dx,sz,bflag);
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            deriv644_y(Du_cpu,u,dx,sz,bflag);
        t_cpu.stop();

        //std::cout<<"grid : "<<grid<<" block : "<<block<<std::endl;
        //printf("grid %s \n",grid);
        
        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();


        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }else if(dir==2)
    {

        deriv644_z(Du_cpu,u,dx,sz,bflag);
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            deriv644_z(Du_cpu,u,dx,sz,bflag);
        t_cpu.stop();

        //std::cout<<"grid : "<<grid<<" block : "<<block<<std::endl;
        //printf("grid %s \n",grid);
        
        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_Du,dptr_u,dptr_blk);
        cudaDeviceSynchronize();
        t_gpu.stop();


        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[pp]-Du_cpu[pp])>lmax)
            lmax=fabs(Du[pp]-Du_cpu[pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");
    delete [] Du_cpu;
    return;

}


template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_ko_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu  = new T[blkSz];
    T* dptr_u  = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du = GPUDevice::device_malloc<T>(blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(ny/pencils,nz,1);
    dim3 block_x = dim3(nx,pencils,1);

    dim3 grid_y  = dim3(nx/pencils,nz,1);
    dim3 block_y = dim3(pencils,ny,1);

    dim3 grid_z  = dim3(nx/pencils,ny,1);
    dim3 block_z = dim3(pencils,nz,1);

    t_cpu.clear();
    t_gpu.clear();

    if(dir==0)
    {
        if(ny%pencils !=0)
            std::cout<<"error : ny "<<ny<<" is not perfect factor of stencil size: "<<pencils<<std::endl;

        ko_deriv42_x(Du_cpu,u,dx,sz,bflag);
        
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            ko_deriv42_x(Du_cpu,u,dx,sz,bflag);
        t_cpu.stop();


        //std::cout<<"grid : "<<grid<<" block : "<<block<<std::endl;
        //printf("grid %s \n",grid);
        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_ko_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_ko_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);

        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }else if(dir==1)
    {

        ko_deriv42_y(Du_cpu,u,dx,sz,bflag);
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            ko_deriv42_y(Du_cpu,u,dx,sz,bflag);
        t_cpu.stop();

        //std::cout<<"grid : "<<grid<<" block : "<<block<<std::endl;
        //printf("grid %s \n",grid);
        
        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_ko_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();


        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_ko_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }else if(dir==2)
    {

        ko_deriv42_z(Du_cpu,u,dx,sz,bflag);
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            ko_deriv42_z(Du_cpu,u,dx,sz,bflag);
        t_cpu.stop();

        //std::cout<<"grid : "<<grid<<" block : "<<block<<std::endl;
        //printf("grid %s \n",grid);
        
        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_ko_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_ko_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_Du,dptr_u,dptr_blk);
        cudaDeviceSynchronize();
        t_gpu.stop();


        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[pp]-Du_cpu[pp])>lmax)
            lmax=fabs(Du[pp]-Du_cpu[pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" ko in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");

    delete [] Du_cpu;
    
    return;

}


template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_second_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{
    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;

    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;
    blk.m_bflag = bflag;

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu  = new T[blkSz];
    T* DDu_cpu = new T[blkSz];
    
    T* dptr_u   = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du  = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_DDu = GPUDevice::device_malloc<T>(blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    t_cpu.clear();
    t_gpu.clear();
    
    if(ny%pencils !=0)
        std::cout<<"error : ny "<<ny<<" is not perfect factor of stencil size: "<<pencils<<std::endl;

    dim3 grid_x  = dim3(ny/pencils,nz,1);
    dim3 block_x = dim3(nx,pencils,1);

    dim3 grid_y  = dim3(nx/pencils,nz,1);
    dim3 block_y = dim3(pencils,ny,1);

    dim3 grid_z  = dim3(nx/pencils,ny,1);
    dim3 block_z = dim3(pencils,nz,1);

    if(dir == 0)
    {
        

        deriv644_xx(DDu_cpu,u,dx,sz,bflag);
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            deriv644_xx(DDu_cpu,u,dx,sz,bflag);
        t_cpu.stop();

        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_xx<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_DDu,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_deriv_xx<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_DDu,dptr_u,dptr_blk);

        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_DDu,blkSz);
        cudaDeviceSynchronize();

    }else if (dir == 1)
    {
        deriv644_x(Du_cpu,u,dx,sz,bflag);
        deriv644_y(DDu_cpu,Du_cpu,dx,sz,bflag);

        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
        {
            deriv644_x(Du_cpu,u,dx,sz,bflag);
            deriv644_y(DDu_cpu,Du_cpu,dx,sz,bflag);
        }
        t_cpu.stop();

        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);
        device::gpu_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_DDu,dptr_Du,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
        {
            device::gpu_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);
            device::gpu_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_DDu,dptr_Du,dptr_blk);
        }
        
        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_DDu,blkSz);
        cudaDeviceSynchronize();



    }else if (dir == 2)
    {
        deriv644_x(Du_cpu,u,dx,sz,bflag);
        deriv644_z(DDu_cpu,Du_cpu,dx,sz,bflag);

        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
        {
            deriv644_x(Du_cpu,u,dx,sz,bflag);
            deriv644_z(DDu_cpu,Du_cpu,dx,sz,bflag);
        }
        t_cpu.stop();

        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);
        device::gpu_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_DDu,dptr_Du,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
        {
            device::gpu_deriv_x<pw,pencils,pencil_sz> <<<grid_x,block_x>>>(dptr_Du,dptr_u,dptr_blk);
            device::gpu_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_DDu,dptr_Du,dptr_blk);
        }
        
        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_DDu,blkSz);
        cudaDeviceSynchronize();



    }
    else if(dir == 3)
    {

        deriv644_yy(DDu_cpu,u,dx,sz,bflag);
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            deriv644_yy(DDu_cpu,u,dx,sz,bflag);
        t_cpu.stop();

        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_yy<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_deriv_yy<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }else if (dir == 4 )
    {
        deriv644_y(Du_cpu,u,dx,sz,bflag);
        deriv644_z(DDu_cpu,Du_cpu,dx,sz,bflag);

        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
        {
            deriv644_y(Du_cpu,u,dx,sz,bflag);
            deriv644_z(DDu_cpu,Du_cpu,dx,sz,bflag);
        }
        t_cpu.stop();

        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
        device::gpu_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_DDu,dptr_Du,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
        {
            device::gpu_deriv_y<pw,pencils,pencil_sz> <<<grid_y,block_y>>>(dptr_Du,dptr_u,dptr_blk);
            device::gpu_deriv_z<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_DDu,dptr_Du,dptr_blk);
        }
        
        cudaDeviceSynchronize();
        t_gpu.stop();

        GPUDevice::device_to_host(Du,dptr_DDu,blkSz);
        cudaDeviceSynchronize();

    }
    else if(dir == 5)
    {

        deriv644_zz(DDu_cpu,u,dx,sz,bflag);
        t_cpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            deriv644_zz(Du_cpu,u,dx,sz,bflag);
        t_cpu.stop();

        GPUDevice::host_to_device(u,dptr_u,blkSz);
        device::gpu_deriv_zz<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_Du,dptr_u,dptr_blk);
        cuda_check_last_error();
        cudaDeviceSynchronize();

        t_gpu.start();
        for(unsigned int iter=0; iter<nReps;iter++)
            device::gpu_deriv_zz<pw,pencils,pencil_sz> <<<grid_z,block_z>>>(dptr_Du,dptr_u,dptr_blk);
        cudaDeviceSynchronize();
        t_gpu.stop();


        GPUDevice::device_to_host(Du,dptr_Du,blkSz);
        cudaDeviceSynchronize();

    }

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[pp]-DDu_cpu[pp])>lmax)
            lmax=fabs(Du[pp]-DDu_cpu[pp]);

    }

    std::cout<<"cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    if(dir==1 || dir==2 || dir==4)
        printf("Average Bandwidth (GB/s): %f\n", 4.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    else
        printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    
    printf("\n");

    delete [] Du_cpu;
    delete [] DDu_cpu;
    return;

}

template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_first_blk_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu  = new T[3*blkSz];
    T* dptr_u  = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du = GPUDevice::device_malloc<T>(3*blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(1,1,1);
    dim3 block_x = dim3(16,16,1);

    t_cpu.clear();
    t_gpu.clear();

    deriv644_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
    deriv644_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
    deriv644_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);


    t_cpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
    {
        deriv644_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
        deriv644_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
        deriv644_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);
        
    }
    t_cpu.stop();

    GPUDevice::host_to_device(u,dptr_u,blkSz);
    device::gpu_blk_d<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);
    cuda_check_last_error();
    cudaDeviceSynchronize();


     t_gpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
        device::gpu_blk_d<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);

    cudaDeviceSynchronize();
    t_gpu.stop();

    GPUDevice::device_to_host(Du,dptr_Du, 3 * blkSz);
    cudaDeviceSynchronize();

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp])>lmax)
            lmax=fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"[BLK] cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");

    delete [] Du_cpu;
    return;

}


template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_second_blk_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu   = new T[3*blkSz];
    T* DDu_cpu  = new T[6*blkSz];
    T* dptr_u   = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du  = GPUDevice::device_malloc<T>(3*blkSz);
    T* dptr_DDu = GPUDevice::device_malloc<T>(6*blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(1,1,1);
    dim3 block_x = dim3(16,16,1);

    t_cpu.clear();
    t_gpu.clear();

    deriv644_x (Du_cpu + 0 *blkSz, u, dx,sz,bflag);
    deriv644_y (Du_cpu + 1 *blkSz, u, dx,sz,bflag);
    deriv644_z (Du_cpu + 2 *blkSz, u, dx,sz,bflag);

    deriv644_xx(DDu_cpu + 0 *blkSz, u, dx,sz,bflag);
    deriv644_y (DDu_cpu + 1 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);
    deriv644_z (DDu_cpu + 2 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);

    deriv644_yy(DDu_cpu + 3 *blkSz, u, dx,sz,bflag);
    deriv644_z (DDu_cpu + 4 *blkSz, Du_cpu + 1 *blkSz, dx,sz,bflag);
    deriv644_zz(DDu_cpu + 5 *blkSz, u, dx,sz,bflag);


    t_cpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
    {
        deriv644_x (Du_cpu + 0 *blkSz, u, dx,sz,bflag);
        deriv644_y (Du_cpu + 1 *blkSz, u, dx,sz,bflag);
        deriv644_z (Du_cpu + 2 *blkSz, u, dx,sz,bflag);

        deriv644_xx(DDu_cpu + 0 *blkSz, u, dx,sz,bflag);
        deriv644_y (DDu_cpu + 1 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);
        deriv644_z (DDu_cpu + 2 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);

        deriv644_yy(DDu_cpu + 3 *blkSz, u, dx,sz,bflag);
        deriv644_z (DDu_cpu + 4 *blkSz, Du_cpu + 1 *blkSz, dx,sz,bflag);
        deriv644_zz(DDu_cpu + 5 *blkSz, u, dx,sz,bflag);
        
    }
    t_cpu.stop();

    GPUDevice::host_to_device(u,dptr_u,blkSz);
    device::gpu_blk_dd<pw,pencils,pencil_sz> <<<grid_x,block_x>>>  (dptr_DDu, dptr_Du, dptr_u,dptr_blk);
    cuda_check_last_error();
    cudaDeviceSynchronize();


     t_gpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
        device::gpu_blk_dd<pw,pencils,pencil_sz> <<<grid_x,block_x>>>  (dptr_DDu, dptr_Du, dptr_u,dptr_blk);

    cudaDeviceSynchronize();
    t_gpu.stop();

    GPUDevice::device_to_host(Du,dptr_DDu, 6 * blkSz);
    cudaDeviceSynchronize();

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);
    GPUDevice::device_free(dptr_DDu);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[dir * blkSz + pp]-DDu_cpu[dir * blkSz + pp])>lmax)
            lmax=fabs(Du[dir * blkSz + pp]-DDu_cpu[dir * blkSz + pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"[BLK] cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");

    delete [] DDu_cpu;
    delete [] Du_cpu;


    return;

}



template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_ko_blk_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu  = new T[3*blkSz];
    T* dptr_u  = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du = GPUDevice::device_malloc<T>(3*blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(1,1,1);
    dim3 block_x = dim3(16,16,1);

    t_cpu.clear();
    t_gpu.clear();

    ko_deriv42_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
    ko_deriv42_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
    ko_deriv42_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);


    t_cpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
    {
        ko_deriv42_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
        ko_deriv42_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
        ko_deriv42_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);
        
    }
    t_cpu.stop();

    GPUDevice::host_to_device(u,dptr_u,blkSz);
    device::gpu_blk_kod<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);
    cuda_check_last_error();
    cudaDeviceSynchronize();


     t_gpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
        device::gpu_blk_kod<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);

    cudaDeviceSynchronize();
    t_gpu.stop();

    GPUDevice::device_to_host(Du,dptr_Du, 3 * blkSz);
    cudaDeviceSynchronize();

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp])>lmax)
            lmax=fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"[BLK] cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");

    delete [] Du_cpu;
    return;

}


template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_first_blk1_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu  = new T[3*blkSz];
    T* dptr_u  = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du = GPUDevice::device_malloc<T>(3*blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(1,1,1);
    dim3 block_x = dim3(8,8,8);

    t_cpu.clear();
    t_gpu.clear();

    deriv644_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
    deriv644_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
    deriv644_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);


    t_cpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
    {
        deriv644_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
        deriv644_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
        deriv644_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);
        
    }
    t_cpu.stop();

    GPUDevice::host_to_device(u,dptr_u,blkSz);
    device::gpu_blk1_d<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);
    cuda_check_last_error();
    cudaDeviceSynchronize();


    t_gpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
        device::gpu_blk1_d<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);

    cudaDeviceSynchronize();
    t_gpu.stop();

    GPUDevice::device_to_host(Du,dptr_Du, 3 * blkSz);
    cudaDeviceSynchronize();

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp])>lmax)
            lmax=fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"[BLK] cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");

    delete [] Du_cpu;
    return;

}


template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_second_blk1_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu   = new T[3*blkSz];
    T* DDu_cpu  = new T[6*blkSz];
    T* dptr_u   = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du  = GPUDevice::device_malloc<T>(3*blkSz);
    T* dptr_DDu = GPUDevice::device_malloc<T>(6*blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(1,1,1);
    dim3 block_x = dim3(8,8,8);

    t_cpu.clear();
    t_gpu.clear();

    deriv644_x (Du_cpu + 0 *blkSz, u, dx,sz,bflag);
    deriv644_y (Du_cpu + 1 *blkSz, u, dx,sz,bflag);
    deriv644_z (Du_cpu + 2 *blkSz, u, dx,sz,bflag);

    deriv644_xx(DDu_cpu + 0 *blkSz, u, dx,sz,bflag);
    deriv644_y (DDu_cpu + 1 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);
    deriv644_z (DDu_cpu + 2 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);

    deriv644_yy(DDu_cpu + 3 *blkSz, u, dx,sz,bflag);
    deriv644_z (DDu_cpu + 4 *blkSz, Du_cpu + 1 *blkSz, dx,sz,bflag);
    deriv644_zz(DDu_cpu + 5 *blkSz, u, dx,sz,bflag);


    t_cpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
    {
        deriv644_x (Du_cpu + 0 *blkSz, u, dx,sz,bflag);
        deriv644_y (Du_cpu + 1 *blkSz, u, dx,sz,bflag);
        deriv644_z (Du_cpu + 2 *blkSz, u, dx,sz,bflag);

        deriv644_xx(DDu_cpu + 0 *blkSz, u, dx,sz,bflag);
        deriv644_y (DDu_cpu + 1 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);
        deriv644_z (DDu_cpu + 2 *blkSz, Du_cpu + 0 *blkSz, dx,sz,bflag);

        deriv644_yy(DDu_cpu + 3 *blkSz, u, dx,sz,bflag);
        deriv644_z (DDu_cpu + 4 *blkSz, Du_cpu + 1 *blkSz, dx,sz,bflag);
        deriv644_zz(DDu_cpu + 5 *blkSz, u, dx,sz,bflag);
        
    }
    t_cpu.stop();

    GPUDevice::host_to_device(u,dptr_u,blkSz);
    device::gpu_blk1_dd<pw,pencils,pencil_sz> <<<grid_x,block_x>>>  (dptr_DDu, dptr_Du, dptr_u,dptr_blk);
    cuda_check_last_error();
    cudaDeviceSynchronize();


    t_gpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
        device::gpu_blk1_dd<pw,pencils,pencil_sz> <<<grid_x,block_x>>>  (dptr_DDu, dptr_Du, dptr_u,dptr_blk);

    cudaDeviceSynchronize();
    t_gpu.stop();

    GPUDevice::device_to_host(Du,dptr_DDu, 6 * blkSz);
    cudaDeviceSynchronize();

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);
    GPUDevice::device_free(dptr_DDu);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[dir * blkSz + pp]-DDu_cpu[dir * blkSz + pp])>lmax)
            lmax=fabs(Du[dir * blkSz + pp]-DDu_cpu[dir * blkSz + pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"[BLK] cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");

    delete [] DDu_cpu;
    delete [] Du_cpu;


    return;

}



template<typename T, int dir, int pw, int pencils, int pencil_sz>
void run_ko_blk1_derivs_gpu(T* const Du, const T* const u, T dx, const unsigned int* sz, unsigned int dof)
{   

    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];
    const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

    BlockGPU3D blk;
    const unsigned int bflag=63;
    
    blk.m_sz[0] = sz[0];
    blk.m_sz[1] = sz[1];
    blk.m_sz[2] = sz[2];

    blk.m_aligned_sz[0] = sz[0];
    blk.m_aligned_sz[1] = sz[1];
    blk.m_aligned_sz[2] = sz[2];

    blk.m_dx[0] = dx;
    blk.m_dx[1] = dx;
    blk.m_dx[2] = dx;

    blk.m_bflag = bflag; 

    BlockGPU3D* dptr_blk = GPUDevice::device_malloc<BlockGPU3D>(1);
    GPUDevice::host_to_device(&blk,dptr_blk,1);

    T* Du_cpu  = new T[3*blkSz];
    T* dptr_u  = GPUDevice::device_malloc<T>(blkSz);
    T* dptr_Du = GPUDevice::device_malloc<T>(3*blkSz);

    const int nReps = 1;
    
    profiler_t t_cpu;
    profiler_t t_gpu;

    dim3 grid_x  = dim3(1,1,1);
    dim3 block_x = dim3(8,8,8);

    t_cpu.clear();
    t_gpu.clear();

    ko_deriv42_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
    ko_deriv42_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
    ko_deriv42_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);


    t_cpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
    {
        ko_deriv42_x(Du_cpu + 0 *blkSz, u, dx,sz,bflag);
        ko_deriv42_y(Du_cpu + 1 *blkSz, u, dx,sz,bflag);
        ko_deriv42_z(Du_cpu + 2 *blkSz, u, dx,sz,bflag);
        
    }
    t_cpu.stop();

    GPUDevice::host_to_device(u,dptr_u,blkSz);
    device::gpu_blk1_kod<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);
    cuda_check_last_error();
    cudaDeviceSynchronize();


    t_gpu.start();
    for(unsigned int iter=0; iter<nReps;iter++)
        device::gpu_blk1_kod<pw,pencils,pencil_sz> <<<grid_x,block_x>>> (dptr_Du, dptr_u,dptr_blk);

    cudaDeviceSynchronize();
    t_gpu.stop();

    GPUDevice::device_to_host(Du,dptr_Du, 3 * blkSz);
    cudaDeviceSynchronize();

    GPUDevice::device_free(dptr_u);
    GPUDevice::device_free(dptr_Du);

    T lmax=0;
    for (unsigned int k=pw; k < (nz-pw); k++)
    for (unsigned int j=pw; j < (ny-pw); j++)
    for (unsigned int i=pw; i < (nx-pw); i++)
    {   
        const int pp = k * ny * nx + j * nx + i;
        //std::cout<<" pp : "<<pp<<" Du_cpu: "<<Du_cpu[pp]<<" Du: "<<Du[pp]<<std::endl;
        if(fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp])>lmax)
            lmax=fabs(Du[dir * blkSz + pp]-Du_cpu[dir * blkSz + pp]);

        // if(fabs(Du[pp]-u[pp])>lmax)
        //     lmax=fabs(Du[pp]-u[pp]);


    }

    std::cout<<"[BLK] cpu time : "<<t_cpu.seconds<<" gpu time : "<<t_gpu.seconds<<" speedup :"<<(t_cpu.seconds/t_gpu.seconds)<<std::endl;
    std::cout<<"l_inf = "<<lmax<<" in direction "<<dir<<std::endl;
    printf("Average Bandwidth (GB/s): %f\n", 2.f * nx * ny * nz * nReps * sizeof(double) / (double)(1024.0*1024.0*1024.0 * t_gpu.seconds));
    printf("\n");

    delete [] Du_cpu;
    return;

}




int main(int argc, char** argv)
{
  // Print device and precision
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  printf("\nDevice Name: %s\n", prop.name);
  printf("Compute Capability: %d.%d\n\n", prop.major, prop.minor);
  
  const Point grid_min(-0.5,-0.5,-0.5);
  const unsigned int blk_sz_1d = 13;
  const unsigned int pencils   =  13;
  const unsigned int sz[3]={blk_sz_1d,blk_sz_1d,blk_sz_1d};
  const Point dx(1.0/sz[0],1.0/sz[1],1.0/sz[2]);
  const unsigned int dof = 1;
  const unsigned blkSz   = sz[0]*sz[1]*sz[2]*dof; 

  double* u  = create_block<double>(grid_min,dx,sz,1);
  double* Du = new double[6 * blkSz];
  
  //GPUDevice::device_to_host(hptr,dptr,blkSz);

  run_first_blk1_derivs_gpu<double, 0, 3, pencils, blk_sz_1d>(Du,u,dx.x(),sz,dof);
  run_first_blk1_derivs_gpu<double, 1, 3, pencils, blk_sz_1d>(Du,u,dx.y(),sz,dof);
  run_first_blk1_derivs_gpu<double, 2, 3, pencils, blk_sz_1d>(Du,u,dx.z(),sz,dof);

  run_ko_blk1_derivs_gpu<double, 0, 3, pencils, blk_sz_1d>(Du,u,dx.x(),sz,dof);
  run_ko_blk1_derivs_gpu<double, 1, 3, pencils, blk_sz_1d>(Du,u,dx.y(),sz,dof);
  run_ko_blk1_derivs_gpu<double, 2, 3, pencils, blk_sz_1d>(Du,u,dx.z(),sz,dof);

  run_second_blk1_derivs_gpu<double, 0, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
  run_second_blk1_derivs_gpu<double, 1, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
  run_second_blk1_derivs_gpu<double, 2, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
  run_second_blk1_derivs_gpu<double, 3, 3 ,pencils , blk_sz_1d>(Du,u,dx.y(),sz,dof);
  run_second_blk1_derivs_gpu<double, 4, 3 ,pencils , blk_sz_1d>(Du,u,dx.y(),sz,dof);
  run_second_blk1_derivs_gpu<double, 5, 3 ,pencils , blk_sz_1d>(Du,u,dx.z(),sz,dof);
  
//   run_first_blk_derivs_gpu<double, 0, 3, pencils, blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_first_blk_derivs_gpu<double, 1, 3, pencils, blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_first_blk_derivs_gpu<double, 2, 3, pencils, blk_sz_1d>(Du,u,dx.z(),sz,dof);

//   run_ko_blk_derivs_gpu<double, 0, 3, pencils, blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_ko_blk_derivs_gpu<double, 1, 3, pencils, blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_ko_blk_derivs_gpu<double, 2, 3, pencils, blk_sz_1d>(Du,u,dx.z(),sz,dof);

//   run_second_blk_derivs_gpu<double, 0, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_second_blk_derivs_gpu<double, 1, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_second_blk_derivs_gpu<double, 2, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_second_blk_derivs_gpu<double, 3, 3 ,pencils , blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_second_blk_derivs_gpu<double, 4, 3 ,pencils , blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_second_blk_derivs_gpu<double, 5, 3 ,pencils , blk_sz_1d>(Du,u,dx.z(),sz,dof);


//   run_first_derivs_gpu<double, 0, 3, pencils, blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_first_derivs_gpu<double, 1, 3, pencils, blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_first_derivs_gpu<double, 2, 3, pencils, blk_sz_1d>(Du,u,dx.z(),sz,dof);

//   run_ko_derivs_gpu<double, 0, 3, pencils, blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_ko_derivs_gpu<double, 1, 3, pencils, blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_ko_derivs_gpu<double, 2, 3, pencils, blk_sz_1d>(Du,u,dx.z(),sz,dof);

//   run_second_derivs_gpu<double, 0, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_second_derivs_gpu<double, 1, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_second_derivs_gpu<double, 2, 3 ,pencils , blk_sz_1d>(Du,u,dx.x(),sz,dof);
//   run_second_derivs_gpu<double, 3, 3 ,pencils , blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_second_derivs_gpu<double, 4, 3 ,pencils , blk_sz_1d>(Du,u,dx.y(),sz,dof);
//   run_second_derivs_gpu<double, 5, 3 ,pencils , blk_sz_1d>(Du,u,dx.z(),sz,dof);


  delete [] u;
  delete [] Du;


  return 0;
}