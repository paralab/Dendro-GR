//
// Created by milinda on 1/15/17.
//


/**
 * @author Milinda Fernando
 * @breif This file contains the precomputed interpolation matrices for a given order, computed based on Hari's HOMG refele.m code.
 * @Note: For wavelet based finite differencing we use uniform points instead of gll points.
 *
 * This is a temporary implementation later we will move to compute the interpolation
 * matrices for a given order on the fly. (See the jacobiPolynomials.h, dendroMatrix.h, and referenceElement.h)
 *
 * */



#ifndef SFCSORTBENCH_INTERPOLATIONMATRICES_H
#define SFCSORTBENCH_INTERPOLATIONMATRICES_H


static double IP_1D_Order_1_0 [] ={1 ,0.5, 0.0 ,0.5 };
static double IP_1D_Order_1_1 [] ={0.5, 0 , 0.5 ,1 };


static double IP_1D_Order_2_0 [] ={1.0 ,0.375 ,0.0 ,
                                   0.0 ,0.75 ,1.0,
                                   0 ,-0.1250,0 };

static double IP_1D_Order_2_1 [] ={0.0 ,-0.1249 ,0.0 ,
                                   1.0 ,0.75 ,0.0 ,
                                   0.0 ,0.3749,1.0 };

static double IP_1D_Order_4_0 [] ={1.0, 0.2734375, 0.0, -0.0390625, 0.0 ,
                                   0.0, 1.09375, 1.0 , 0.46875, 0.0 ,
                                   0.0, -0.546875, 0.0, 0.703125 ,1 ,
                                   0.0, 0.21875, 0.0, -0.15625, 0.0,
                                   0.0, -0.0390625, 0.0, 0.0234375, 0.0 };

static double IP_1D_Order_4_1 [] ={0.0, 0.0234375, 0.0, -0.03906245, 0.0,
                                   0.0, -0.15625, 0.0, 0.21875, 0.0,
                                   1.0 , 0.703125, 0.0, -0.546875, 0.0,
                                   0.0, 0.46875, 1.0, 1.09375, 0.0,
                                   0.0, -0.0390625, 0.0, 0.2734375, 1.0 };

// special intergrid transfer operator based on symetric stencil coefficients. 
static double IP_1D_FD_Order_5[]={ 3.0/256.0 ,-25.0/256.0 , 75.0/128.0 , 75.0/128.0  , -25.0/256.0  , 3.0/256.0};

static double IP_1D_FD_Order_4[]={ -1.0/16.0 ,9.0/16.0 , 9.0/16.0 , -1.0/16.0};



#endif //SFCSORTBENCH_INTERPOLATIONMATRICES_H
