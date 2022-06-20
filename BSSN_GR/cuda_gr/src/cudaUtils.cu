//
// Created by milinda on 8/9/18.
//

#include "../include/cudaUtils.cuh"

namespace cuda
{
    __device__ unsigned int get_smid(void) {

        unsigned int  ret;

        asm("mov.u32 %0, %smid;" : "=r"(ret) );

        return ret;

    }
}