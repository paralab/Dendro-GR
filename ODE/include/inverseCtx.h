/**
 * @file sensitivityCtx.h
 * @author Milinda Fernando
 * @brief Ctx class to perform sensitivity analysis for  time dependent problems. 
 * @version 0.1
 * @date 2021-10-12
 * @copyright Copyright (c) 2021
 */

#pragma once
#include "dendro.h"
#include "mesh.h"
#include "launcher.h"

namespace invp
{
    template<typename DerivedInverseCtx, typename Ctx>
    class InverseCtx
    {
        protected:
            /**@brief forward ctx*/
            Ctx* m_uiForCtx    = NULL;

            /**@brief adjoint ctx*/
            Ctx* m_uiAdjCtx    = NULL;

            /**@brief Global communicator*/
            MPI_Comm m_uiCommGlobal;
            
            /**@brief Local communicator*/
            MPI_Comm m_uiCommLocal;

            launcher::Launcher* m_uiJobLauncher=NULL;

            
        public:

            InverseCtx(){};

            ~InverseCtx(){};

            /**@brief: derived class static cast*/            
            DerivedInverseCtx &asLeaf() { return static_cast<DerivedInverseCtx &>(*this); }

            /**@brief get global communicator*/
            MPI_Comm get_global_comm() const { return m_uiCommGlobal;}
            
            /**@brief get local communicator*/
            MPI_Comm get_local_comm() const { return m_uiCommLocal;}

            /**@brief set the job launcher*/
            int set_launcher(launcher::Launcher* launcher) {
                m_uiJobLauncher = launcher; 
                m_uiCommGlobal  = m_uiJobLauncher->get_global_communicator();
                m_uiCommLocal   = *(m_uiJobLauncher->get_sub_communicator());
                return 0;
            }
            
            /**@brief create forward problem ctx*/
            int initialize_forward_ctx() { return asLeaf().initialize_forward_ctx(); }
            
            /**@brief create adjoint problem ctx*/
            int initialize_adjoint_ctx() { return asLeaf().initialize_adjoint_ctx(); }
            
            /**@brief launch the forward problem */
            int launch_forward_solve() { return asLeaf().launch_forward_solve();}

            /**@brief launch the adjoint problem */
            int launch_adjoint_solve() {return asLeaf().launch_adjoint_solve();}
            
            /**@brief launch gradient approximation */
            int gradient_approximation() {return asLeaf().gradient_approximation();}

    };

}// end of namespece sa. 