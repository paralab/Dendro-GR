/**
 * @file nlsmCtx_cu.cuh
 * @brief BSSN rhs context file for cuda. 
 * @version 0.1
 * @date 2022-02-07
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#pragma once
#include "dendro.h"
#include "device.h"
#include "mesh_gpu.cuh"
#include "dvec.h"
#include "ctx.h"
#include "grDef.h"
#include "grUtils.h"
#include "oct2vtk.h"
#include "checkPoint.h"
#include "mathMeshUtils.h"
#include "refel_const.cuh"
#include "dataUtils.h"
#include "TwoPunctures.h"
#include "bssn_constraints.h"
#include "physcon.h"
#include "bssnrhs_evar_derivs.cuh"
#include "gwExtract.h"
#include "bssn_kernels.cuh"
#include "device_utils.cuh"
#include "meshUtils.h"
namespace bssn
{   

    enum VL {CPU_EV=0, CPU_CV, CPU_EV_UZ_IN, CPU_EV_UZ_OUT, CPU_CV_UZ_IN, CPU_CV_UZ_OUT, GPU_EV, GPU_EV_UZ_IN, GPU_EV_UZ_OUT,END};
    typedef ot::DVector<DendroScalar, unsigned int> DVec;

    class BSSNCtxGPU : public ts::Ctx<BSSNCtxGPU,DendroScalar, unsigned int>
    {
        protected:
            device::MeshGPU m_mesh_cpu;

            /**@brief: mesh in the device*/
            device::MeshGPU * m_dptr_mesh;

            /**@brief: evolution var (zip)*/
            DVec m_var[VL::END];

            Point m_uiBHLoc[2];
            
            BSSN_EVAR_DERIVS* m_deriv_evars=nullptr;
            BSSN_EVAR_DERIVS* m_dptr_deriv_evars =nullptr;
            DEVICE_REAL * m_deriv_base=nullptr;

            static  const unsigned int DEVICE_MAX_DERIVS  = 210;
            static  const unsigned int DEVICE_RHS_BATCHED_GRAIN_SZ = 4096;
            static  const unsigned int DEVICE_RHS_NSTREAMS         = 1;
            static  const unsigned int DEVICE_RHS_BLK_SZ  = 13 * 13* 13;

        public:
            /**@brief: default constructor*/
            BSSNCtxGPU(ot::Mesh* pMesh);

            /**@brief: default deconstructor*/
            ~BSSNCtxGPU();

            /**@brief: initial solution*/
            int initialize();

            int init_grid();
            
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
            inline int pre_stage(DVec& sIn) {return 0;}

            /**@brief: function execute after each stage
             * @param sIn: stage var in. 
            */
            int post_stage(DVec& sIn);

            /**@brief: function execute before each step*/
            inline int pre_timestep(DVec& sIn) {return 0;} 

            /**@brief: function execute after each step*/
            int post_timestep(DVec& sIn);

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

            /**@brief: pack and returns the evolution variables to one DVector*/
            DVec& get_evolution_vars_cpu();

            /**@brief: pack and returns the constraint variables to one DVector*/
            DVec& get_constraint_vars();

            /**@brief: pack and returns the primitive variables to one DVector*/
            DVec& get_primitive_vars();

            /**@brief: prints any messages to the terminal output. */
            int terminal_output();

            /**@brief: returns the async communication batch size. */
            unsigned int get_async_batch_sz() {return BSSN_ASYNC_COMM_K;}

            /**@brief: returns the number of variables considered when performing refinement*/
            unsigned int get_num_refine_vars() {return BSSN_NUM_REFINE_VARS;}

            /**@brief: return the pointer for containing evolution refinement variable ids*/
            const unsigned int* get_refine_var_ids() { return BSSN_REFINE_VARIABLE_INDICES;}
            
            /**@brief return the wavelet tolerance function / value*/
            std::function<double(double, double, double,double * hx)> get_wtol_function(){

                double wtol = BSSN_WAVELET_TOL;
                std::function<double(double,double,double,double *)> waveletTolFunc =[](double x,double y, double z,double *hx){ return bssn::computeWTolDCoords(x,y,z,hx);};
                return waveletTolFunc;
            }

            static unsigned int getBlkTimestepFac(unsigned int blev, unsigned int lmin, unsigned int lmax);

            int grid_transfer(const ot::Mesh* m_new);

            int host_to_device_sync();
            
            int device_to_host_sync();

            int device_to_host_async(cudaStream_t s);
            
            int host_to_device_async(cudaStream_t s);

            inline device::MeshGPU*& get_meshgpu_device_ptr() {return m_dptr_mesh;}

            inline device::MeshGPU* get_meshgpu_host_handle() {return &m_mesh_cpu;}

            bool is_bh_merged(double tol) const {return ((bssn::BSSN_BH_LOC[0]-bssn::BSSN_BH_LOC[1]).abs() < tol); };

            /**@biref: evolve bh locations. */
            void evolve_bh_loc(DVec sIn, double dt);


        
    };
};
