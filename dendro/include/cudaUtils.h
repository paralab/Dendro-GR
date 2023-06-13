/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief Contains utility function for the host related to GPUs
 * @date  8/9/18.
 * */

#pragma once
#include "cuda_runtime.h"
#include "block.h"


//Macro for checking cuda errors following a cuda launch or api call
#define CUDA_CHECK_ERROR() {                                       \
 cudaError_t e=cudaGetLastError();                                 \
 if(e!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(e));           \
   exit(0); \
 }                                                                 \
}



namespace cuda
{

    /**
     * @brief send device information to the gpu
     * @param[in] device : gpu device ID
     * @return cudaDeviceProp allocated on the device
     * */
    cudaDeviceProp* getGPUDeviceInfo(unsigned int device);

    /**
     * @breif send mesh blocks to the gpu
     * @param[in] in : input array
     * @param[in] out: device pointer where the data is copied to.
     * */
    template<typename T>
    T * copyArrayToDevice(const T* in, unsigned int numElems);


    /**
     * @breif copy value to device
     * @param[in] in : input value
     * @param[in] out: device pointer where the data is copied to.
     * */
    template<typename T>
    inline T * copyValueToDevice(const T* in);

     /**
      * @biref allocates 1D array
      * @param [in] sz1
      * @return the pointer to the device allocation
      * */

     template <typename T>
     T* alloc1DCudaArray(unsigned int sz1);

    /**
     * @brief allocates a 2D cuda array on the device.
     * @param[in] sz1: dim 1 size
     * @param[in] sz2: dim 2 size
     * @returns the double pointer to the 2D array.
     * */
    template <typename T>
    T** alloc2DCudaArray(unsigned int sz1,  unsigned int sz2);


    /**
     * @brief allocates a 2D cuda array on the device.
     * @param[out] hostPtr: 2D pointer accesible from the host.
     * @param[in] sz1: dim 1 size
     * @param[in] sz2: dim 2 size
     * @returns the double pointer to the 2D array (device pointer).
     * */
    template <typename T>
    T** alloc2DCudaArray(T**& hostPtr,unsigned int sz1,  unsigned int sz2);



    /**
     * @brief allocates a 2D cuda array on the device and copy data.
     * @param[in] sz1: dim 1 size
     * @param[in] sz2: dim 2 size
     * @returns the double pointer to the 2D array.
     * */
    template <typename T>
    T** alloc2DCudaArray(const T** in,unsigned int sz1,  unsigned int sz2);



    /**
     * @breif send mesh blocks to the gpu (async)
     * @param[in] in : input array
     * @param[in] out: device pointer where the data is copied to.
     * */
    template<typename T>
    void copyArrayToDeviceAsync(const T* in, T*__deviceptr,unsigned int numElems,const cudaStream_t stream);

    /**
       * @brief allocates a 2D cuda array on the device and copy data.
       * @param[in] sz1: dim 1 size
       * @param[in] sz2: dim 2 size
       * @returns the double pointer to the 2D array.
       * */
    template <typename T>
    void copy2DCudaArrayToDeviceAsync(const T **in, T **__devicePtr, unsigned int sz1, unsigned int sz2,
                                      const cudaStream_t stream);


    template<typename T>
    void copyArrayToHostAsync(T* host_ptr, const T*__deviceptr,unsigned int numElems,const cudaStream_t stream);



    /**
     * @brief deallocates the 2D cuda array.
     * @param[in] sz1: dim 1 size
     * */
    template <typename T>
    void dealloc2DCudaArray(T ** & __array2D,  unsigned int sz1);

    /***
     * computes the how dendro blocks (octree blocks result in from unzip) to the gpu/
     * @param[in] blkList: list of dendro data blocks
     * @param[in] numBlocks: number of data blocks,
     * @param[out] blockMap: (blockMap[2*blocDim.x] , blockMap[2*blocDim.x+1]) begin & end of data block that is going to be process by the gpu block
     * */
    void computeDendroBlockToGPUMap(const ot::Block* blkList, unsigned int numBlocks, unsigned int*& blockMap,dim3 & gridDim);


     /**
     * @breif copy array from device to memory
     * @param[in] __device_ptr : device pointer
     * @param[in] numElems : number of elements
     * @param[out] host_ptr: host ptr
     * */
     template<typename T>
     void copyArrayToHost(T* host_ptr,const T* __device_ptr, unsigned int numElems);


    /**
    * @breif copy 2D array from device to memory
    * @param[in] __device_ptr : 2D pointer to the device
    * @param[in] sz1 : size1
    * @param[in] sz2 : size2
    * @param[out] host_ptr: host ptr
    * */
     template<typename T>
     void copy2DArrayToHost(T** host_ptr,const T** __device_ptr, unsigned int sz1,unsigned int sz2);




}




// templated functions

namespace cuda
{

    template<typename T>
    T * copyArrayToDevice(const T* in, unsigned int numElems)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T)*numElems);
        CUDA_CHECK_ERROR();

        cudaMemcpy(__devicePtr,in,sizeof(T)*numElems,cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }


    template<typename T>
    inline T * copyValueToDevice(const T* in)
    {

        T* __devicePtr;
        cudaMalloc(&__devicePtr,sizeof(T));
        CUDA_CHECK_ERROR();

        cudaMemcpy(__devicePtr,in,sizeof(T),cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return __devicePtr;

    }


    template <typename T>
    T* alloc1DCudaArray(unsigned int sz1)
    {
        T* __tmp1d;
        cudaMalloc(&__tmp1d,sizeof(T)*sz1);
        CUDA_CHECK_ERROR();

        return __tmp1d;
    }


    template <typename T>
    T** alloc2DCudaArray(T**& hostPtr,unsigned int sz1,  unsigned int sz2)
    {
        T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        hostPtr=new T*[sz1];

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&hostPtr[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();
        }

        cudaMemcpy(__tmp2d,hostPtr,sizeof(T*)*sz1,cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();

        return __tmp2d;

    }

    template <typename T>
    T** alloc2DCudaArray(unsigned int sz1, unsigned int sz2)
    {

        T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        T** tmp2D=new T*[sz1];

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&tmp2D[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();
        }

        cudaMemcpy(__tmp2d,tmp2D,sizeof(T*)*sz1,cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();
        delete [] tmp2D;

        return __tmp2d;

    }

    template <typename T>
    T** alloc2DCudaArray(const T** in,unsigned int sz1,  unsigned int sz2)
    {
        T** __tmp2d;
        cudaMalloc(&__tmp2d,sizeof(T*)*sz1);
        CUDA_CHECK_ERROR();

        T** tmp2D=new T*[sz1];

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMalloc(&tmp2D[i],sizeof(T)*sz2);
            CUDA_CHECK_ERROR();
            cudaMemcpy(tmp2D[i],in[i], sizeof(T)*sz2 ,cudaMemcpyHostToDevice);
            CUDA_CHECK_ERROR();
        }

        cudaMemcpy(__tmp2d,tmp2D,sizeof(T*)*sz1,cudaMemcpyHostToDevice);
        CUDA_CHECK_ERROR();
        delete [] tmp2D;

        return __tmp2d;
    }


    template <typename T>
    void dealloc2DCudaArray(T ** & __array2D, unsigned int sz1)
    {
        T** tmp2D=new T*[sz1];

        cudaMemcpy(tmp2D,__array2D,sizeof(T*)*sz1,cudaMemcpyDeviceToHost);
        CUDA_CHECK_ERROR();

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaFree(tmp2D[i]);
            CUDA_CHECK_ERROR();
        }

        delete [] tmp2D;

        cudaFree(__array2D);
        CUDA_CHECK_ERROR();
    }


    template<typename T>
    void copyArrayToDeviceAsync(const T* in,T*__deviceptr, unsigned int numElems,const cudaStream_t stream)
    {
        cudaMemcpyAsync(__deviceptr,in,sizeof(T)*numElems,cudaMemcpyHostToDevice,stream);
        CUDA_CHECK_ERROR();

    }


    template<typename T>
    void copyArrayToHostAsync(T* host_ptr, const T*__deviceptr,unsigned int numElems,const cudaStream_t stream)
    {
        cudaMemcpyAsync(host_ptr,__deviceptr,sizeof(T)*numElems,cudaMemcpyDeviceToHost,stream);
        CUDA_CHECK_ERROR();
    }


    template<typename T>
    void copyArrayToHost(T* host_ptr,const T* __device_ptr, unsigned int numElems)
    {
        cudaMemcpy(host_ptr,__device_ptr,sizeof(T)*numElems,cudaMemcpyDeviceToHost);
        CUDA_CHECK_ERROR();
        
    }



    template<typename T>
    void copy2DArrayToHost(T** host_ptr,const T** __device_ptr, unsigned int sz1,unsigned int sz2)
    {
        T** tmp2D=new T*[sz1];
        cudaMemcpy(tmp2D,__device_ptr,sizeof(T*)*sz1,cudaMemcpyDeviceToHost);
        CUDA_CHECK_ERROR();

        for(unsigned int i=0;i<sz1;i++)
        {
            cudaMemcpy(host_ptr[i],tmp2D[i],sizeof(T)*sz2,cudaMemcpyDeviceToHost);
            CUDA_CHECK_ERROR();
        }

        delete [] tmp2D;
        return;

    }


}


