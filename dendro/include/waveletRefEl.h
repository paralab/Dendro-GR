/**
 * @file waveletRefEl.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Generic Wavelet element class to compute 
 * arbitary order wavelet coefficients
 * Currently coded to compute the wavelets for even element orders. 
 * @version 0.1
 * @date 2020-06-08
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include<iostream>
#include "refel.h"
#include "mathUtils.h"
namespace wavelet
{

    /*
    * Even element order computations (only supported )
    * For given element order p=2k, k \in Z^{+}. (Even element order)
    * For even element order we have a max padding with of p/2.
    * 
    *  | 0, 1,...,p/2 -1  | p/2 , ..... p/2 + p |  p/2 + p + 1 ,    .. 2p
    * 
    */
    class WaveletEl
    {

        protected:
            
            /**@brief: pointer to the refernce element*/
            RefElement* m_uiRefEl;

            /**@brief: reference element to boundary computations.*/
            RefElement* m_uiRefElBdy;

            /**@brief : indices of the nodes where the wavelets are computed. */
            std::vector<unsigned int> m_uiCIndex;

            /**@brief : indices of the nodes where the parent grid points are extracted (input nodal values for the wavelet interpolation)*/
            std::vector<unsigned int> m_uiPIndex;

            /**@brief: workspace vector in*/
            std::vector<double> m_uiVIn;

            /**@brief: workspace vector m_uiNVec*/
            std::vector<double> m_uiNVec;

            /**@brief: workspace vector out*/
            std::vector<double> m_uiVOut;
            

        public:
            
            /**
             * @brief Construct a new Wavelet Element object
             * @param order : octree element order. 
             * @param dim : dim of the operartor
             */
            WaveletEl(RefElement* refEl);
            
            /**@biref default destrucutor*/
            ~WaveletEl();
            
            /**@biref: function to compute the wavelets in 3d unziped ele vector. */
            void compute_wavelets_3D(const double* in, const unsigned int * isz ,std::vector<double>& wc,bool isBdy=false);


    };

}// end of namespace wavelets. 

