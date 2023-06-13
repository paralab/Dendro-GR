/**
 * @file etsRefEl.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Computes the block correction operators for explicit time steppers, 
 * @version 0.1
 * @date 2020-03-28
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */


/**
 * Note that these correction operations are constructed from picewise linear approximation of F(u), 
 * [k] = CxPx[DU] -> (1)
 * and the B computation comes from the (1) and the Taylor exapansion of the RK stages. 
 * 
 */

#pragma once

#include <iostream>
#include <cmath>
#include "ts.h"

namespace ts
{
    class ENUTSOp
    {
        protected:

            /**@brief: ETS time stage computation coefficients. */
            std::vector<DendroScalar> m_uiAij;

            /**@brief Diagonal op. of time scaling $P_{\Delta t}$ */
            std::vector<DendroScalar> m_uiPi;

            /**@brief Diagonal op. of time scaling $P_{\Delta t}^{-1}$ */
            std::vector<DendroScalar> m_uiInvPi;

            /**@brief: Coefficient matrix from Taylor expansion of stage vectors i.e [K] = C*P DU*/
            std::vector<DendroScalar> m_uiCij;

            /**@brief: Inverse of coefficient matrix C*/
            std::vector<DendroScalar>m_uiInvCij;

            /**@brief: Taylor expansion of the DU = (\partial_t u , \partial^2_t u ... \partial^{p}_t u)^T*/
            std::vector<DendroScalar> m_uiBij;

            /**@brief: workspace matrix M*/
            std::vector<DendroScalar> m_uiMij;
            
            /**@brief: workspace vector in */
            std::vector<DendroScalar> m_uiVin;

            /**@brief: workspace vector out */
            std::vector<DendroScalar> m_uiVout;
            
            /**@brief: number of stages of the scheme*/
            unsigned int m_uiNumStages;

            /**@brief: time stepper type*/
            ETSType m_uiType;

            
        public:
            
            /**
             * @brief Construct a new ENUTSOp object
             * @param ts_type time stepper type. 
             */
            ENUTSOp(ETSType type);

            /**
             * @brief Destroy the ENUTSOp object
             * 
             */
            ~ENUTSOp(){};

            /**
             * @brief Diagonal op. of time scaling $P_{\Delta t}$ and its inverse operator
             * @param dt : delta t that P is defined at. 
             * @param rk_s : nunmber of stages considered to be interpolated. 
             */
            void Pdt(DendroScalar dt, unsigned int rk_s);

            /**
             * @brief : computes the Taylor expansion operator for DU
             * @param dt : delta t that P is defined at. 
             * @param rk_s : nunmber of stages considered to be interpolated. 
             */
            void Bdt(DendroScalar dt, unsigned int rk_s);

            
            // stage correction routines begin

            /**
             * @brief perform the finer to coarser time stage corrections. 
             * @param out : corrected or interpolated coarser stages from finer stages. 
             * @param in  : finer stages. 
             * @param sz: size of in/out vectors
             * @param rk_s : number of rk stages needs to be interpolated. 
             * @param dt_c : coarser time step size. 
             * @param dt_f : finer time step size. 
             * @param dt : dt for the taylor expansion of DU
             * @param dof : number of dof.
             */
            void Cfc(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, unsigned int rk_s, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof=1);

            /**
             * @brief perform the coarser to finer time stage corrections. 
             * @param out : corrected or interpolated coarser stages from finer stages. 
             * @param in  : finer stages. 
             * @param sz: size of in/out vectors
             * @param rk_s : number of rk stages needs to be interpolated. 
             * @param dt_c : coarser time step size. 
             * @param dt_f : finer time step size. 
             * @param dt : dt for the taylor expansion of DU
             * @param dof : number of dof.
             */
            void Ccf(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, unsigned int rk_s, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof=1);

            
            // stage correction routines end

            // below routines are for coarser to finer and finer to coarser time corrections, which is needed resuming LTS time stepper after intergrid transfer.
            // U(t) correction 
            // For these corrections we need all the computed stages. 
            
            /**
             * @brief Performs correction operators for Ut
             * @param out 
             * @param in 
             * @param sz 
             * @param dt_c 
             * @param dt_f 
             * @param dt 
             * @param dof 
             */
            void coarser_finer_ut_correction(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof=1);

            /**
             * @brief 
             * @param out 
             * @param in 
             * @param sz 
             * @param dt_c 
             * @param dt_f 
             * @param dt 
             * @param dof 
             */
            void finer_coarser_ut_correction(DendroScalar** out, const DendroScalar** in, const unsigned int* sz, DendroScalar dt_c, DendroScalar dt_f, DendroScalar dt, unsigned int dof=1);


    };



}; // end of name space ts. 

