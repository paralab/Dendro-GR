/**
 * @file nlsm_inv.cpp
 * @brief Inverse problems in NLSM
 * @version 0.1
 * @date 2021-10-17
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "nlsm.h"
#include "nlsmUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "octUtils.h"
#include "nlsmCtx.h"
#include "ets.h"
#include "mathUtils.h"
#include "sdc.h"
#include "enuts.h"
#include "assert.h"
#include "nlsmInvCtx.h"
#include "launcher.h"


int main (int argc, char** argv)
{
    // 0- NUTS 1-UTS
    if(argc<5)
        std::cout<<"Usage: "<<argv[0]<<" paramFile total_jobs nodes_per_job cores_per_node"<<std::endl;

    const unsigned int TS_MODE        = atoi(argv[2]);
    const unsigned int TOTAL_JOBS     = atoi(argv[3]); 
    const unsigned int NODES_PER_JOB  = atoi(argv[4]);
    const unsigned int CORES_PER_NODE = atoi(argv[5]); 

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    if(!rank)
    {
      std::cout<<GRN<<"=============CMD PARS GIVEN ================= "<<NRM<<std::endl;
      std::cout<<GRN<<"TS_MODE : "<<TS_MODE<<NRM<<std::endl;
      std::cout<<GRN<<"TOTAL JOBS : "<<TOTAL_JOBS<<NRM<<std::endl;
      std::cout<<GRN<<"NODES_PER_JOB: "<<NODES_PER_JOB<<NRM<<std::endl;
      std::cout<<GRN<<"CORES_PER_NODE: "<<CORES_PER_NODE<<NRM<<std::endl;
      std::cout<<GRN<<"============================================= "<<NRM<<std::endl;
    }


    if (!rank) {
      
      #ifdef NLSM_NONLINEAR
        std::cout<<GRN<<"Compiled with NLSM_NONLINEAR"<<NRM<<std::endl;
      #else
        std::cout<<RED<<"Compiled without NLSM_NONLINEAR"<<NRM<<std::endl;
      #endif

      #ifdef NLSM_COMPARE_WITH_ANALYTICAL_SOL
        std::cout<<GRN<<"Compiled with NLSM_COMPARE_WITH_ANALYTICAL_SOL"<<NRM<<std::endl;
      #else
        std::cout<<RED<<"Compiled without NLSM_COMPARE_WITH_ANALYTICAL_SOL"<<NRM<<std::endl;
      #endif

      #ifdef NLSM_USE_4TH_ORDER_DERIVS
        std::cout<<GRN<<"Compiled with NLSM_USE_4TH_ORDER_DERIVS"<<NRM<<std::endl;
      #else
        std::cout<<RED<<"Compiled without NLSM_USE_4TH_ORDER_DERIVS"<<NRM<<std::endl;
      #endif

      #ifdef NLSM_USE_6TH_ORDER_DERIVS
        std::cout<<GRN<<"Compiled with NLSM_USE_6TH_ORDER_DERIVS"<<NRM<<std::endl;
      #else
        std::cout<<RED<<"Compiled without NLSM_USE_6TH_ORDER_DERIVS"<<NRM<<std::endl;
      #endif

      #ifdef NLSM_USE_8TH_ORDER_DERIVS
        std::cout<<GRN<<"Compiled with NLSM_USE_8TH_ORDER_DERIVS"<<NRM<<std::endl;
      #else
        std::cout<<RED<<"Compiled without NLSM_USE_8TH_ORDER_DERIVS"<<NRM<<std::endl;
      #endif

    }

    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    nlsm::readParamFile(argv[1],comm);
    nlsm::dumpParamFile(std::cout,1,comm);
    
    _InitializeHcurve(nlsm::NLSM_DIM);
    m_uiMaxDepth=nlsm::NLSM_MAXDEPTH;

    if(nlsm::NLSM_NUM_VARS%nlsm::NLSM_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total NLSM_NUM_VARS: "<<nlsm::NLSM_NUM_VARS<<" is not divisable by NLSM_ASYNC_COMM_K: "<<nlsm::NLSM_ASYNC_COMM_K<<std::endl;
        exit(0);
    }
    
    launcher::Launcher la(comm,NODES_PER_JOB,CORES_PER_NODE);
    la.alloc_sub_communicator(TOTAL_JOBS);

    nlsm::InvNLSMCtx nlsmInvCtx;
    nlsmInvCtx.set_launcher(&la);
    nlsmInvCtx.set_parameter_filename(argv[1]);

    nlsmInvCtx.initialize_forward_ctx();
    nlsmInvCtx.launch_forward_solve();
    
    MPI_Finalize();

    return 0;
}



