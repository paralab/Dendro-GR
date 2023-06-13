/**
 * @file device.h
 * @brief Device abstractions to write device independent GPU code.
 * @version 0.1
 * @date 2021-10-08
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once
#include<iostream>
#include<stdio.h>
#include<assert.h>
#include "dendro.h"
#include "mpi.h"

#ifdef __CUDACC__
    #define DEVICE_FUNC __device__
    #define HOST_FUNC __host__
    #define GLOBAL_FUNC __global__
    #define SHARED_MEM __shared__
    #define CONST_MEM  __constant__
#else
    #define DEVICE_FUNC 
    #define HOST_FUNC 
    #define GLOBAL_FUNC  
    #define SHARED_MEM  
    #define CONST_MEM
#endif

namespace device
{

    template<typename LDevice>
    class Device
    {   
        

        public:
            
            /**@brief: derived class static cast*/
            inline LDevice &asLeaf() { return static_cast<LDevice &>(*this); }

            template<typename T>
            static inline void check_error(T result){return LDevice::template check_error<T>(result);}

            static inline void check_last_error() {return LDevice::template check_last_error();}

            /**@brief: allocates memory on the device*/
            template<typename T>   
            static T* device_malloc(unsigned int n) {
                return LDevice::template device_malloc<T>(n);
            }

            /**@brief: free memory on the device*/
            template<typename T>
            static void device_free(T* dptr) {
                return LDevice::template device_free<T>(dptr);
            }

            template<typename T>   
            static T* host_malloc(unsigned int n) {
                return LDevice::template host_malloc<T>(n);
            }
            
            template<typename T>   
            static void host_free(T* hptr) {
                return LDevice::template host_free<T>(hptr);
            }

            /**@brief: host to device memory transfer*/
            template<typename T>
            static void host_to_device(const T* const hptr, T* const dptr, unsigned int n){ 
                return LDevice::template host_to_device<T>(hptr,dptr,n);
            }

            /**@brief: device to host memory transfer*/
            template<typename T>
            static void device_to_host(T* const hptr, const T* const dptr, unsigned int n) { 
                return LDevice::template device_to_host<T>(hptr,dptr,n);
            }

            template<typename T,typename device_stream>
            static void host_to_device_async(const T* const hptr, T* const dptr, unsigned int n, device_stream s= (device_stream) 0) {
                return LDevice::template host_to_device_async<T,device_stream>(hptr,dptr,n,s);
            }

            template<typename T,typename device_stream>
            static void device_to_host_async(T* const hptr, const T* const dptr, unsigned int n, device_stream s= (device_stream) 0){
                return LDevice::template device_to_host_async<T,device_stream>(hptr,dptr,n,s);
            }

            DEVICE_FUNC static int thread_id_x() { return LDevice::thread_id_x(); }
            DEVICE_FUNC static int thread_id_y() { return LDevice::thread_id_y(); }
            DEVICE_FUNC static int thread_id_z() { return LDevice::thread_id_z(); }
            DEVICE_FUNC static int block_id_x()  { return LDevice::block_id_x(); }
            DEVICE_FUNC static int block_id_y()  { return LDevice::block_id_y(); }
            DEVICE_FUNC static int block_id_z()  { return LDevice::block_id_z(); }
            DEVICE_FUNC static int block_dim_x() { return LDevice::block_dim_x(); }
            DEVICE_FUNC static int block_dim_y() { return LDevice::block_dim_y(); }
            DEVICE_FUNC static int block_dim_z() { return LDevice::block_dim_z(); }

            DEVICE_FUNC static int grid_dim_x() { return LDevice::grid_dim_x(); }
            DEVICE_FUNC static int grid_dim_y() { return LDevice::grid_dim_y(); }
            DEVICE_FUNC static int grid_dim_z() { return LDevice::grid_dim_z(); }

            DEVICE_FUNC static void sync_threads() { return LDevice::sync_threads();}

            static void device_synchronize(){return LDevice::device_synchronize();}
            
            template<typename device_stream>
            static void stream_synchronize(device_stream s=(device_stream)0){return LDevice:: template stream_synchronize<device_stream>(s);}

            
    };

}

#define GPU_MAX_THREADS_PER_BLOCK 1024
#ifdef __CUDACC__

    #include "cuda_runtime.h"
    #include "cuda_runtime_api.h"
    #include "device_launch_parameters.h"
    #include "cuda.h"

    namespace device
    {
        class DeviceCUDA 
        {
            public:

                DEVICE_FUNC static int thread_id_x() { return threadIdx.x; }
                DEVICE_FUNC static int thread_id_y() { return threadIdx.y; }
                DEVICE_FUNC static int thread_id_z() { return threadIdx.z; }

                DEVICE_FUNC static int block_id_x()  { return blockIdx.x; }
                DEVICE_FUNC static int block_id_y()  { return blockIdx.y; }
                DEVICE_FUNC static int block_id_z()  { return blockIdx.z; }
                
                DEVICE_FUNC static int block_dim_x() { return blockDim.x; }
                DEVICE_FUNC static int block_dim_y() { return blockDim.y; }
                DEVICE_FUNC static int block_dim_z() { return blockDim.z; }

                DEVICE_FUNC static int grid_dim_x() { return gridDim.x; }
                DEVICE_FUNC static int grid_dim_y() { return gridDim.y; }
                DEVICE_FUNC static int grid_dim_z() { return gridDim.z; }

                static void device_synchronize() { cudaDeviceSynchronize(); return;}

                /**@brief: synchronize threads in a block*/
                DEVICE_FUNC static void sync_threads() {
                    __syncthreads();
                    return;
                }

                /**@brief: cuda check error*/
                template<typename T>
                static inline void check_error(T result)
                {
                    //#if defined(DEBUG) || defined(_DEBUG)
                        if (result != cudaSuccess) {
                            fprintf(stderr, "CUDA Runtime Error: %s\n", cudaGetErrorString(result));
                            __handler(0);
                            assert(result == cudaSuccess);
                            
                        }
                    //#endif
                    return ;
                }

                static inline void check_last_error() {
                    //#if defined(DEBUG) || defined(_DEBUG)
                        cudaError_t err = cudaGetLastError();
                        if (err != cudaSuccess) 
                        {
                            printf("!! %ss\n", cudaGetErrorString(err));
                            __handler(0);
                        }
                            
                    //#endif
                    return;
                }

                /**@brief: allocates memory on the device*/
                template<typename T>
                static T* device_malloc(unsigned int n) {
                    T* data;
                    check_error(cudaMalloc((void**)&data,sizeof(T)*n));
                    return data;
                }

                /**@brief: free memory on the device*/
                template<typename T>
                static void device_free(T* dptr) {
                    check_error(cudaFree(dptr));
                    return ;
                }

                /**@brief: host to device memory transfer*/
                template<typename T>
                static void host_to_device( const T* const hptr, T* const dptr, unsigned int n){ 
                    
                    check_error(cudaMemcpy(dptr,hptr,sizeof(T)*n,cudaMemcpyHostToDevice));
                    return ;
                }

                /**@brief: device to host memory transfer*/
                template<typename T>
                static void device_to_host(T* const hptr, const T* const dptr, unsigned int n) { 

                    check_error(cudaMemcpy(hptr,dptr,sizeof(T)*n,cudaMemcpyDeviceToHost));
                    return ;
                }

                /**@brief: allocates memory on the host (pinned mem. allocation)*/
                template<typename T>
                static T* host_malloc(unsigned int n) {
                    #if defined(DEBUG) || defined(_DEBUG)
                        T* data = (T*)malloc(sizeof(T)*n);
                    #else
                        T*data;
                        check_error(cudaMallocHost((void**)&data,sizeof(T)*n));
                    #endif
                    return data;
                }

                /**@brief: free memory on the host (pinned allocated memory free)*/
                template<typename T>
                static void host_free(T* hptr) {
                    #if defined(DEBUG) || defined(_DEBUG)
                        free(hptr);
                    #else
                        check_error(cudaFreeHost(hptr));
                    #endif
                    return ;
                }
                
                template<typename T,typename device_stream>
                static void host_to_device_async(const T* const hptr, T* const dptr, unsigned int n, device_stream s= (device_stream) 0) {
                    check_error(cudaMemcpyAsync(dptr,hptr,sizeof(T)*n,cudaMemcpyHostToDevice,s));
                    return;
                }

                template<typename T,typename device_stream>
                static void device_to_host_async(T* const hptr, const T* const dptr, unsigned int n, device_stream s= (device_stream) 0){
                   check_error(cudaMemcpyAsync(hptr,dptr,sizeof(T)*n,cudaMemcpyDeviceToHost,s));
                   return;
                }

                template<typename device_stream>
                static void stream_synchronize(device_stream s= (device_stream) 0){
                    check_error(cudaStreamSynchronize(s));
                    return;
                }

        };
    }

    typedef device::Device<device::DeviceCUDA> GPUDevice;

    template<typename T>
    GLOBAL_FUNC void axpy_cu(int n, T a, const T* __restrict__ const in, T* __restrict__  const out){
        const int t_id   = threadIdx.x + blockDim.x * blockIdx.x;
        const int stride = blockDim.x * gridDim.x;
        for (int i = t_id; i < n; i += stride)
            out[i] += a * in[i];
        
    }

    template<typename T>
    GLOBAL_FUNC void axpy_cu_2d(unsigned int lb, unsigned int le, unsigned int szpdof, T a, const T* __restrict__ const in, T* __restrict__  const out){
        const int t_id   = threadIdx.x + blockDim.x * blockIdx.x;
        const int stride = blockDim.x * gridDim.x;
        const int vid    = blockIdx.y;
        for (int i = (t_id +lb) ; i < (le); i+= stride)
            out[vid *szpdof +  i] += a * in[vid * szpdof + i];
        
    }

    static inline int cuda_mpi_aware_init()
    {
        #ifdef MPIX_CUDA_AWARE_SUPPORT
            char * local_rank_str = NULL;
            int rank = 0, devCount = 0;

            // We extract the local rank initialization using an environment variable
            if ((local_rank_str = getenv("OMPI_COMM_WORLD_LOCAL_RANK")) != NULL)
            {
                rank = atoi(local_rank_str);		
            }else
            {
                std::cout<<"local rank extraction through env. vars failed for cuda init";
                exit(0);
            }

            cudaGetDeviceCount(&devCount);
            if(!rank)
                printf("number of cuda devices: %d\n",devCount);
            cudaSetDevice(rank%devCount);
            cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
        #endif

        return 0;
        
    }

    static inline int cuda_init(MPI_Comm comm)
    {
        #ifndef MPIX_CUDA_AWARE_SUPPORT
            int rank = 0, devCount = 0;
            MPI_Comm_rank(comm,&rank);
            
            cudaGetDeviceCount(&devCount);
            if(!rank)
                printf("number of cuda devices: %d\n",devCount);
            cudaSetDevice(rank%devCount);
            cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
        #endif
        
        return 0;

    }

#endif







