/**
 * @file aehCtx.h
 * @Milinda Fernando (milinda@cs.utah.edu)
 * @brief Context file for the Apprent Event Horizon (AEH) solver. 
 * @version 0.1
 * @date 2020-05-09
 * @copyright Copyright (c) 2020
 * 
 */
#pragma once
#include "ctx.h"
#include "parameters.h"
#include "grDef.h"
#include "grUtils.h"

namespace bssn
{
    /**@brief: AEH Ctx class*/
    class AEHCtx : public ts::Ctx<DendroScalar, DendroIntL>
    {

        protected:

            /**@brief: evolution var (zip)*/
            DVec m_uiEVar;

            /**@brief: bssn variables coming from the the bssnSolver.*/
            DVec m_uiBSSNVars;

            /**@brief: Evolution var unzip 0 - unzip in , 1 - unzip out */
            DVec m_uiEUnzip[2];

            /**@brief: theta resolution for the initial AEH (S^2)*/
            DendroScalar m_uiThetaRes;
            
            /**@brief: phi resolution for the initial AEH (S^2)*/
            DendroScalar m_uiPhiRes;
            
            

        public: 

            /**@brief: Constructor for the AEH Ctx
             * @param pMesh: mesh data structure. 
             * @param bssnVars: bssn variables. 
            */
            AEHCtx(ot::Mesh* pMesh, DVec bssnVars);

            /**@brief: desctructor fot the AEHCtx*/
            ~AEHCtx();

            /**@brief: Compute theta expansion for a given pts, and surface normals
             * 
            */
            void compute_theta(const Point* pts);

            /**@brief: initial solution*/
            virtual int initialize();
            
            /**@brief: compute the rhs for the theta*/
            virtual int rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time);

            /**@brief: function execute before each stage
             * @param sIn: stage var in. 
            */
            virtual int pre_stage(DVec sIn); 

            /**@brief: function execute after each stage
             * @param sIn: stage var in. 
            */
            virtual int post_stage(DVec sIn);

            /**@brief: function execute before each step*/
            virtual int pre_timestep(DVec sIn); 

            /**@brief: function execute after each step*/
            virtual int post_timestep(DVec sIn);

            /**@brief: function execute after each step*/
            virtual bool is_remesh();

            /**@brief: write to vtu. */
            virtual int write_vtu();

            /**@brief: should be called for free up the contex memory. */
            virtual int finalize();

            /**@brief: pack and returns the evolution variables to one DVector*/
            virtual DVec get_evolution_vars();

            /**@brief: updates the application variables from the variable list. */
            virtual int update_app_vars();

            /**@brief: prints any messages to the terminal output. */
            virtual int terminal_output();

            /**@brief: returns the async communication batch size. */
            virtual unsigned int get_async_batch_sz() {return 1;}



    };

}// end of namespace bssn. 