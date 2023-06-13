/**
 * @file nufft.h
 * @author Milinda Fernando (milind@cs.utah.edu)
 * @brief contains the interface to non-uniform FFT library by L. Greengard.
 * @date 2019-04-15
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef DENDRO_5_0_NUFFT_H
#define DENDRO_5_0_NUFFT_H

#include<iostream>
#include "dendro.h"    
    
   
    /**
     * @brief : computes the nufft type 1 (nufft)
     * @param[in] nj : number of sources. 
     * @param[in] xj : x coord values
     * @param[in] yj : y coord values
     * @param[in] zj : z coord values 
     * @param[in] cj : function values. 
     * @param[in] iflag : if >= 0 exponential with  positive sign, used in fft 
     * @param[in] eps : user defined tolerance
     * @param[in] ms : ft coefficient array size x
     * @param[in] mt : ft coefficient array size y
     * @param[in] mu : ft coefficient array size z
     * @param[out] fk : computed fft coefficients
     * @param[out] ier : error return code. 
     * */
    extern "C" void nufft3d1f90_(int nj,double *xj,double *yj,double* zj,DendroComplex* cj, int iflag,double eps,int ms,int mt,int mu,double *fk,int ier);

    /**
     * @brief : computes the nufft type 2 (nufft) (inverse ft)
     * @param[in] nj : number of sources. 
     * @param[in] xj : x coord values
     * @param[in] yj : y coord values
     * @param[in] zj : z coord values 
     * @param[out] cj : function values. 
     * @param[in] iflag : if >= 0 exponential with  positive sign, used in fft 
     * @param[in] eps : user defined tolerance
     * @param[in] ms : ft coefficient array size x
     * @param[in] mt : ft coefficient array size y
     * @param[in] mu : ft coefficient array size z
     * @param[in] fk : computed fft coefficients
     * @param[out] ier : error return code. 
     * */
    extern "C" void nufft3d2f90_(int nj,double *xj,double *yj,double* zj,DendroComplex* cj, int iflag,double eps,int ms,int mt,int mu,double *fk,int ier);


#endif