/**
 * @file nlsmCtx_cu.cuh
 * @brief NLSM rhs context file for cuda. 
 * @version 0.1
 * @date 2022-02-07
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include "dendro.h"
#include "device.h"
#include "mesh.h"
#include "mesh_gpu.cuh"
#include "derivs_cu.cuh"
#include "bc_cu.cuh"
#include "dvec.h"
#include "parameters.h"
#include "ctx.h"
#include "nlsmUtils.h"
#include "oct2vtk.h"
#include "checkPoint.h"
#include "mathMeshUtils.h"
#include "refel_const.cuh"
#include "device_utils.cuh"

namespace nlsm
{   

    struct EVAR_DERIVS
    {
        DEVICE_REAL * grad_0_chi;
        DEVICE_REAL * grad_1_chi;
        DEVICE_REAL * grad_2_chi;

        DEVICE_REAL * grad_0_phi;
        DEVICE_REAL * grad_1_phi;
        DEVICE_REAL * grad_2_phi;

        DEVICE_REAL * ko_grad_0_chi;
        DEVICE_REAL * ko_grad_1_chi;
        DEVICE_REAL * ko_grad_2_chi;

        DEVICE_REAL * ko_grad_0_phi;
        DEVICE_REAL * ko_grad_1_phi;
        DEVICE_REAL * ko_grad_2_phi;

        DEVICE_REAL * grad2_0_0_chi;
        DEVICE_REAL * grad2_1_1_chi;
        DEVICE_REAL * grad2_2_2_chi;

    };

    enum VL {CPU_EV=0, CPU_EV_DG, CPU_EV_UZ, GPU_EV, GPU_EV_DG, GPU_EV_UZ_IN, GPU_EV_UZ_OUT,END};
    typedef ot::DVector<DendroScalar, unsigned int> DVec;

    class NLSMCtxGPU : public ts::Ctx<NLSMCtxGPU,DendroScalar, unsigned int>
    {
        #ifdef __PROFILE_CTX__
        // public:
        //     using ts::Ctx<NLSMCtxGPU, DendroScalar, DendroIntL>::m_uiCtxpt;
        #endif

        protected:
            device::MeshGPU m_mesh_cpu;

            /**@brief: mesh in the device*/
            device::MeshGPU * m_dptr_mesh;

            /**@brief: evolution var (zip)*/
            DVec m_var[VL::END];

            static const unsigned int DEVICE_RHS_BATCHED_GRAIN_SZ=2048;
            static const unsigned int DEVICE_RHS_BLK_SZ  = 13 * 13* 13;

            EVAR_DERIVS * m_deriv_evars      = nullptr;
            EVAR_DERIVS * m_dptr_deriv_evars = nullptr;
            DEVICE_REAL * m_dptr_deriv_base  = nullptr; 

            
        public:
            /**@brief: default constructor*/
            NLSMCtxGPU(ot::Mesh* pMesh);

            /**@brief: default deconstructor*/
            ~NLSMCtxGPU();

            /**@brief: initial solution*/
            int initialize();
            
            /**
             * @brief computes the rhs 
             * 
             * @param in : zipped input
             * @param out : zipped output
             * @param sz  : number of variables. 
             * @param time : current time. 
             * @return int : status. (0) on success. 
             */
            int rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time);

            
            /**@brief: function execute before each stage
             * @param sIn: stage var in. 
            */
            int pre_stage(DVec sIn); 

            /**@brief: function execute after each stage
             * @param sIn: stage var in. 
            */
            int post_stage(DVec sIn);

            /**@brief: function execute before each step*/
            int pre_timestep(DVec sIn); 

            /**@brief: function execute after each step*/
            int post_timestep(DVec sIn);

            /**@brief: function execute after each step*/
            bool is_remesh();

            /**@brief: write to vtu. */
            int write_vtu();

            /**@brief: writes checkpoint*/
            int write_checkpt();

            /**@brief: restore from check point*/
            int restore_checkpt();

            /**@brief: should be called for free up the contex memory. */
            int finalize();

            /**@brief: pack and returns the evolution variables to one DVector*/
            DVec& get_evolution_vars();

            /**@brief: pack and returns the constraint variables to one DVector*/
            DVec& get_constraint_vars();

            /**@brief: pack and returns the primitive variables to one DVector*/
            DVec& get_primitive_vars();

            /**@brief: prints any messages to the terminal output. */
            int terminal_output();

            /**@brief: returns the async communication batch size. */
            unsigned int get_async_batch_sz() {return NLSM_ASYNC_COMM_K;}

            /**@brief: returns the number of variables considered when performing refinement*/
            unsigned int get_num_refine_vars() {return NLSM_NUM_REFINE_VARS;}

            /**@brief: return the pointer for containing evolution refinement variable ids*/
            const unsigned int* get_refine_var_ids() { return NLSM_REFINE_VARIABLE_INDICES;}
            
            /**@brief return the wavelet tolerance function / value*/
            std::function<double(double, double, double,double *)> get_wtol_function(){

                double wtol = NLSM_WAVELET_TOL;
                std::function<double(double,double,double,double*)> waveletTolFunc =[wtol](double x,double y, double z,double * hx){ return nlsm::computeWTol(x,y,z,hx);};
                return waveletTolFunc;
                
            }

            static unsigned int getBlkTimestepFac(unsigned int blev, unsigned int lmin, unsigned int lmax);

            int grid_transfer(const ot::Mesh* m_new);

            int host_to_device_sync();
            
            int device_to_host_sync();

            inline device::MeshGPU*& get_meshgpu_device_ptr() {return m_dptr_mesh;}

            inline device::MeshGPU* get_meshgpu_host_handle() {return &m_mesh_cpu;}


        
    };
};
