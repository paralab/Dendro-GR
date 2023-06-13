/**
 * @file sdc.h
 * @brief Spectral differed correction methods to solve ODEs. 
 * @version 0.1
 * @date 2021-09-28
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#pragma once
#include <iostream>
#include "dendro.h"
#include "mesh.h"
#include "ctx.h"
#include "dvec.h"
#include "ts.h"
#include "mathUtils.h"
#include "refel.h"
namespace ts
{

    enum SDC_STATUS {SUCCESS=0, DIVERGED_DUE_MAX_ITER=1, FAILED=2};

    
    template<typename T, typename Ctx>
    class SDC
    {

        protected: 
            
            /**@brief Application context. */
            Ctx* m_uiAppCtx = NULL;

            /**@brief SDC order, i.e., m_uiOrder + 1 quadrature points*/
            unsigned int m_uiOrder = 0;

            /**@brief SDC matrix, */
            double * m_uiQMat = NULL;

            /**@brief SDC order, i.e., m_uiOrder + 1 quadrature points*/
            DVec m_uiEVar;

            /**@brief max number of iterations allowed. */
            unsigned int MAX_ITERATIONS = 1000;

            /**@brief picard iter initial condition*/
            std::vector<DVec> m_uiU0;
            /**@brief brief description*/

            /**@picard iter current iteration state*/
            std::vector<DVec> m_uiUk;

            /**@brief iter current iteration with operator rhs applied*/
            std::vector<DVec> m_uiLUk;

            /**@brief pickard next iteration. */
            std::vector<DVec> m_uiVk;

            /**@brief brief description*/
            bool m_uiIsInternalAlloc=false;

            /**@brief: time step info*/
            TSInfo m_uiTimeInfo;

            /**@brief: convergence tolerance*/
            DendroScalar m_uiResidualTolerance = 1e-6;

            
            
        private:
            /**
             * @brief Allocates internal variables for the time stepper. 
             * @return int 
             */
            int allocate_internal_vars();

            /**@brief: Deallocate internal variables. */
            int deallocate_internal_vars();

        public: 

            /**
             * @brief Construct a new SDC object
             */
            SDC(){};


            /**
             * @brief Destroy the SDC object
             */
            ~SDC()
            { 
                delete [] m_uiQMat;
                m_uiQMat=NULL;
            };
            
            /**
             * @brief Set the application ctx object
             * @param appCtx : application Ctx object
             * @return int 
             */
            int set_application_ctx(Ctx* appCtx);
            /**
             * @brief Set the time integration order object
             * @param tOrder time Gauss Labatto quadrature order. 
             * @return int 
             */
            int set_time_integration_order(unsigned int tOrder);

            /**
             * @brief Initialize the SDC solvers. 
             * @return int 
             */
            int init();

            int evolve();


            /**
             * @brief Performs a picard iteration. 
             * 
             * @param in : current state variables. U(t_n)
             * @param out : next timestep state variables U(t_n +1 )
             * @param tn : current time 
             * @param dt : timestep size
             * @param residual_tol : residual for convergence
             * @return int 
             */
            int picard_solve(DVec * in, DVec * out, double tn, double dt, double residual_tol=1e-6);

            /**@brief: returns the mesh is active. */
            inline bool is_active() const {return m_uiAppCtx->get_mesh()->isActive();}

            /**@brief: returns the current time step*/
            inline DendroIntL curr_step() {return m_uiTimeInfo._m_uiStep;};

            /**@brief: returns the current time*/
            inline DendroScalar curr_time() {return m_uiTimeInfo._m_uiT;};

            /**@brief: get the timestepper size*/
            inline DendroScalar ts_size() const { return m_uiTimeInfo._m_uiTh;}

            /**@brief get the residual tolerance*/
            inline DendroScalar get_residual_tolerance()const {return m_uiResidualTolerance;}

            /**@brief set the residual tolerance*/
            inline int set_residual_tolerance(DendroScalar re) {m_uiResidualTolerance=re; return SDC_STATUS::SUCCESS;}
            
            /**@brief: returns the active rank*/
            unsigned int get_active_rank() const ;
            
            /**@brief: returns the active npes*/
            unsigned int get_active_npes() const;
            
            /**@brief: returns the global rank*/
            unsigned int get_global_rank() const;

            /**@brief: return the global npes*/
            unsigned int get_global_npes() const;

            /**@brief: return the active communicator*/
            MPI_Comm get_active_comm () const;

            /**@brief: return the global communicator. */
            MPI_Comm get_global_comm() const;

            /**@brief: perform synchronizations with correct variable allocations for the new mesh: should be called after remeshing.  */
            int sync_with_mesh();

            /**@brief: returns the underlying mesh data structure. */
            const ot::Mesh* get_mesh() const  { return m_uiAppCtx->get_mesh();} 

            /**@brief: returns the evolution variables. */ 
            inline DVec get_evolve_vars() const { return m_uiEVar; }
            
            /**@brief: set evolution vars. */
            int set_evolve_vars(DVec eVars);
    };

    template<typename T, typename Ctx>
    int SDC<T,Ctx> :: allocate_internal_vars()
    {

        if(m_uiIsInternalAlloc)
            return 0;

        const unsigned int num_q = m_uiOrder+1;
        m_uiU0.resize(num_q);
        m_uiUk.resize(num_q);
        m_uiLUk.resize(num_q);
        m_uiVk.resize(num_q);

        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();

        for(unsigned int s=0; s < num_q ; s++)
        {
            m_uiU0[s].create_vector(pMesh, m_uiEVar.get_type(), m_uiEVar.get_loc(), m_uiEVar.get_dof(), m_uiEVar.is_ghost_allocated());
            m_uiUk[s].create_vector(pMesh, m_uiEVar.get_type(), m_uiEVar.get_loc(), m_uiEVar.get_dof(), m_uiEVar.is_ghost_allocated());
            m_uiLUk[s].create_vector(pMesh, m_uiEVar.get_type(), m_uiEVar.get_loc(), m_uiEVar.get_dof(), m_uiEVar.is_ghost_allocated());
            m_uiVk[s].create_vector(pMesh, m_uiEVar.get_type(), m_uiEVar.get_loc(), m_uiEVar.get_dof(), m_uiEVar.is_ghost_allocated());
        }

        m_uiIsInternalAlloc = true;
        return SDC_STATUS::SUCCESS;
    }
    
    template<typename T, typename Ctx>
    int SDC<T,Ctx> :: deallocate_internal_vars()
    {

        if(!m_uiIsInternalAlloc)
            return 0;

        const unsigned int num_q = m_uiOrder+1;
        for(unsigned int s=0; s < num_q ; s++)
        {
            m_uiU0[s].destroy_vector();
            m_uiUk[s].destroy_vector();
            m_uiLUk[s].destroy_vector();
            m_uiVk[s].destroy_vector();
        }

        m_uiIsInternalAlloc=false;
        return SDC_STATUS::SUCCESS; 
    }


    template<typename T, typename  Ctx>
    int SDC<T,Ctx>::set_application_ctx(Ctx* appCtx)
    {   
        m_uiAppCtx = appCtx;
        m_uiEVar   = m_uiAppCtx->get_evolution_vars();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        return SDC_STATUS::SUCCESS;
    }

    template<typename T, typename Ctx>
    int SDC<T,Ctx>::set_evolve_vars(DVec eVars)
    {
        m_uiEVar = eVars;
        return 0; 
    }
    
    template<typename T, typename  Ctx>
    int SDC<T,Ctx>::set_time_integration_order(unsigned int tOrder)
    {
        m_uiOrder=tOrder;
        const unsigned int num_q_pts = m_uiOrder + 1;
        
        if(m_uiQMat!=NULL)
        {
            delete [] m_uiQMat;
            m_uiQMat=NULL;
        }
            
        m_uiQMat = new T[num_q_pts*num_q_pts];

        RefElement refEl(1,m_uiOrder);
        const double* qgl = refEl.getQgll();
        const double* wgl = refEl.getWgll();

        //printArray_1D(wgl,m_uiOrder+1);
        for(unsigned int sj=0; sj < num_q_pts ; sj++)
            m_uiQMat[0*num_q_pts + sj]=0.0;

        // Notes on QMat. 
        // M be the order of the Qmat, hence M+1 points. 
        // t_n= t0 , t1, t2, ..., t_{m+1} =t_{n+1}, dt =  t_{n+1}-t{n}
        // QMat - First row is zero, since we use Gauss-Labatto integration, \int{t0}^{to} will give us zero. 
        // QMat - mth row should be scaled by (tm-t0)/2 which is = dt*(q[m]-q[0])/4 

        for(unsigned int si=1; si < num_q_pts ; si++)
         for(unsigned int sj=0; sj < num_q_pts ; sj++)
            m_uiQMat[si*num_q_pts + sj]= ((qgl[si]-qgl[0])/4.0) * wgl[sj];


        // m_uiQMat[1*num_q_pts + 0] = 1.0/2.0;
        // m_uiQMat[1*num_q_pts + 1] = 1.0/2.0;

        // m_uiQMat[1*num_q_pts + 0] = 1.0/2.0;
        // m_uiQMat[2*num_q_pts + 1] = 1.0/2.0;
        // m_uiQMat[3*num_q_pts + 2] = 1.0;


        // m_uiQMat[0*num_q_pts + 0] = 1.0/6.0;
        // m_uiQMat[0*num_q_pts + 1] = -1.0/3.0;
        // m_uiQMat[0*num_q_pts + 2] = 1.0/6.0;

        // m_uiQMat[1*num_q_pts + 0] = 1.0/6.0;
        // m_uiQMat[1*num_q_pts + 1] = 5.0/12.0;
        // m_uiQMat[1*num_q_pts + 2] = -1.0/12.0;

        // m_uiQMat[2*num_q_pts + 0] = 1.0/6.0;
        // m_uiQMat[2*num_q_pts + 1] = 2.0/3.0;
        // m_uiQMat[2*num_q_pts + 2] = 1.0/6.0;

        //printArray_2D(m_uiQMat,num_q_pts,num_q_pts);


        return SDC_STATUS::SUCCESS;
    }

    template<typename T, typename  Ctx>
    int SDC<T,Ctx>::init()
    {
        m_uiAppCtx->initialize();
        m_uiEVar   = m_uiAppCtx->get_evolution_vars();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        deallocate_internal_vars();
        allocate_internal_vars();
        return SDC_STATUS::SUCCESS;

        
    }

    template<typename T, typename  Ctx>
    int SDC<T,Ctx>::evolve()
    {
        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        const double current_t= m_uiTimeInfo._m_uiT;
        double current_t_adv=current_t;
        const double dt = m_uiTimeInfo._m_uiTh;

        const int status = this->picard_solve(&m_uiEVar,&m_uiEVar,current_t,dt,m_uiResidualTolerance);
        m_uiAppCtx->increment_ts_info();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();

        return status;
    }

    template<typename T, typename  Ctx>
    int SDC<T,Ctx>::picard_solve(DVec * in, DVec * out, double tn, double dt, double residual_tol)
    {
        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        if(pMesh->isActive())
        {   

            const unsigned int arank = pMesh->getMPIRank();
            const unsigned int anpes = pMesh->getMPICommSize();
            MPI_Comm acomm           = pMesh->getMPICommunicator();

            const unsigned int nodeLocalBegin = pMesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd   = pMesh->getNodeLocalEnd();
            const unsigned int numLocalNodes  = pMesh->getNumLocalMeshNodes();
            const unsigned int dof = in->get_dof();
            const unsigned int sz  = in->get_size()/dof;

            const unsigned int num_q = m_uiOrder + 1; 
            
            const T* in_ptr = in->get_vec_ptr();
            DendroScalar residual_local[num_q];

            DendroIntL local_nodes  = numLocalNodes;
            DendroIntL global_nodes;
            par::Mpi_Allreduce(&local_nodes,&global_nodes,1,MPI_SUM,acomm);

            m_uiAppCtx->pre_timestep(*in);
            
            // for(unsigned int v=0; v< dof;v++)
            // {
            //     T min = vecMin((T*)&in_ptr[v* sz + nodeLocalBegin],numLocalNodes,acomm);
            //     T max = vecMax((T*)&in_ptr[v* sz + nodeLocalBegin],numLocalNodes,acomm);

            //     if(!arank)
            //         std::cout<<"VAR : "<<v<<" min : "<<min<<" max: "<<max<<std::endl;
            // }
            
            // copy initial condition u0 to all the stages.
            for(unsigned int s=0; s< num_q; s++)
            {
                //m_uiU0[s].copy_data(*in);
                m_uiUk[s].copy_data(*in);
                //m_uiVk[s].copy_data(*in);
                
            }

            unsigned int iter_count         = 0;
            DendroScalar residual_global    = 0;
            DendroScalar* error_norms   = new DendroScalar[num_q*dof];
            DendroScalar* error_norms_g = new DendroScalar[num_q*dof];
            do
            {

                for(unsigned int s=0; s< num_q; s++)
                {
                    m_uiAppCtx->pre_stage(m_uiLUk[s]);
                    m_uiAppCtx->rhs(&m_uiUk[s],&m_uiLUk[s],1,tn);
                    m_uiAppCtx->post_stage(m_uiLUk[s]);
                }
                    

                for(unsigned int s=0; s< num_q; s++)
                    m_uiVk[s].copy_data(*in);

                for (unsigned int si=1; si< num_q; si++)
                    for (unsigned int sj=0; sj< num_q; sj++)
                        DVec::axpy(m_uiAppCtx->get_mesh(), (dt)*m_uiQMat[si*num_q + sj], m_uiLUk[sj], m_uiVk[si]);
                

                for(unsigned int s=0; s< num_q; s++)
                    for(unsigned int v=0; v < dof; v++)
                        error_norms[s*dof + v] = normL2(&m_uiVk[s].get_vec_ptr()[v*sz + nodeLocalBegin],&m_uiUk[s].get_vec_ptr()[v*sz+nodeLocalBegin],numLocalNodes);
                
                par::Mpi_Allreduce(error_norms,error_norms_g, num_q*dof, MPI_SUM, pMesh->getMPICommunicator());

                for(unsigned int s=0; s< num_q; s++)
                    for(unsigned int v=0; v < dof; v++)
                        error_norms_g[s*dof + v ] /= sqrt(global_nodes);
                
                residual_global = vecMax(error_norms_g,num_q*dof);
                if(!arank && iter_count % 1==0)
                    printf("[sdc] picard iteration : %d \t residual: %.10E \n",iter_count,residual_global);
                std::swap(m_uiVk,m_uiUk);
                iter_count++;

            }while ((iter_count < MAX_ITERATIONS) && (residual_global > residual_tol));
            if(!arank)
                printf("[sdc] converged at iter : %d \t residual: %.10E tol = %.10E \n",iter_count,residual_global,residual_tol);
            std::memcpy(out->get_vec_ptr(), m_uiUk[m_uiOrder].get_vec_ptr(),sizeof(T)*(sz*dof));
            m_uiAppCtx->post_timestep(*out);
            
            // if(iter_count==MAX_ITERATIONS)
            // {
            //     if(!arank)
            //         std::cout<<RED<<"[sdc]: picard iteration diverged. killing the solver program"<<NRM<<std::endl;
            //     MPI_Abort(acomm,SDC_STATUS::DIVERGED_DUE_MAX_ITER);
            // }

            delete [] error_norms;
            delete [] error_norms_g;

        }

        pMesh->waitAll();
        return SDC_STATUS::SUCCESS;
        
    }

    template<typename T, typename Ctx>
    unsigned int SDC<T,Ctx>::get_active_rank() const {
        if(is_active())
            return m_uiAppCtx->get_mesh()->getMPIRank();
        else
            return get_global_rank();
    }

    template<typename T, typename Ctx>
    unsigned int SDC<T,Ctx>::get_active_npes() const {
        if(is_active())
            return m_uiAppCtx->get_mesh()->getMPICommSize();
        else
            return get_global_npes();
    }

    template<typename T, typename Ctx>
    unsigned int SDC<T,Ctx>::get_global_rank() const {
        return m_uiAppCtx->get_mesh()->getMPIRankGlobal();
    }

    template<typename T, typename Ctx>
    unsigned int SDC<T,Ctx>::get_global_npes() const {
        return m_uiAppCtx->get_mesh()->getMPICommSizeGlobal();
    }

    template<typename T, typename Ctx>
    MPI_Comm SDC<T,Ctx>:: get_active_comm() const 
    {
        if(is_active())
            return m_uiAppCtx->get_mesh()->getMPICommunicator();
        else
            return MPI_COMM_NULL;
    }

    template<typename T, typename Ctx>
    MPI_Comm SDC<T,Ctx>:: get_global_comm() const 
    {
        return m_uiAppCtx->get_mesh()->getMPIGlobalCommunicator();
    }

    template <typename T, typename Ctx >
    int SDC<T,Ctx>::sync_with_mesh()
    {
        if(m_uiAppCtx -> is_ets_synced())
            return 0;
        

        deallocate_internal_vars();
        allocate_internal_vars();
        m_uiEVar = m_uiAppCtx->get_evolution_vars();
        m_uiAppCtx->set_ets_synced(true);

        return 0;
    }




}// end of namespace ts
