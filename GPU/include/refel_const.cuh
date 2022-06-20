/**
 * @brief put reference element interpolation matrices to the constant GPU memory. 
*/

#pragma once
#include"device.h"
#define REFEL_CONST_MEM_MAX 11*11
namespace device
{
    extern CONST_MEM DEVICE_REAL refel_1d[2*REFEL_CONST_MEM_MAX];
}