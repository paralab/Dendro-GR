/**
 * @file nlsm_nuts.cpp
 * @brief : NLSM to test with SDC. 
 * @version 0.1
 * @date 2021-09-29
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
#include "enuts.h"
#include "assert.h"
#include "mathUtils.h"
#include "sdc.h"

int main (int argc, char** argv)
{
    // 0- NUTS 1-UTS
    unsigned int ts_mode=0;

    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    ts_mode = atoi(argv[2]);

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

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

    nlsm::timer::initFlops();

    nlsm::timer::total_runtime.start();


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

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){nlsm::initData(x,y,z,var);};
    std::function<void(double,double,double,double,double*)> u_x_t=[](double x,double y,double z,double t,double*var){nlsm::analyticalSol(x,y,z,t,var);};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){nlsm::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=nlsm::NLSM_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<nlsm::NLSM_NUM_VARS;i++)
        varIndex[i]=i;

    
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    if(nlsm::NLSM_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(nlsm::NLSM_BLK_MIN_X,nlsm::NLSM_BLK_MIN_Y,nlsm::NLSM_BLK_MIN_Z);
        const Point pt_max(nlsm::NLSM_BLK_MAX_X,nlsm::NLSM_BLK_MAX_Y,nlsm::NLSM_BLK_MAX_Z);

        nlsm::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-(binOp::fastLog2(nlsm::NLSM_ELE_ORDER)),m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,nlsm::NLSM_NUM_VARS,nlsm::NLSM_REFINE_VARIABLE_INDICES,nlsm::NLSM_NUM_REFINE_VARS,tmpNodes,m_uiMaxDepth,nlsm::NLSM_WAVELET_TOL,nlsm::NLSM_ELE_ORDER,comm);
        
    }
    
    unsigned int lmin=1,lmax=5;
    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),nlsm::NLSM_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,nlsm::NLSM_DENDRO_GRAIN_SZ,nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
    mesh->setDomainBounds(Point(nlsm::NLSM_GRID_MIN_X,nlsm::NLSM_GRID_MIN_Y,nlsm::NLSM_GRID_MIN_Z), Point(nlsm::NLSM_GRID_MAX_X, nlsm::NLSM_GRID_MAX_Y,nlsm::NLSM_GRID_MAX_Z));
    mesh->computeMinMaxLevel(lmin,lmax);
    nlsm::NLSM_RK45_TIME_STEP_SIZE=8*nlsm::NLSM_CFL_FACTOR*((nlsm::NLSM_COMPD_MAX[0]-nlsm::NLSM_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) nlsm::NLSM_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    par::Mpi_Bcast(&nlsm::NLSM_RK45_TIME_STEP_SIZE,1,0,comm);
    
    std::cout<<"DT: "<<nlsm::NLSM_RK45_TIME_STEP_SIZE<<std::endl;

    DendroIntL lblocks = mesh->getLocalBlockList().size();
    DendroIntL gblocks =0; 
    par::Mpi_Reduce(&lblocks,&gblocks,1,MPI_SUM,0,comm);

    if(!rank)
    {
      std::cout<<" number of blocks for coarset block level : "<<(m_uiMaxDepth-MAXDEAPTH_LEVEL_DIFF-1)<<" # blocks: "<<gblocks<<std::endl;
      std::cout<<" lmin: "<<lmin<<" lmax: "<<lmax<<std::endl;
      std::cout<<"ts_mode: "<<ts_mode<<std::endl;
    }
      
    nlsm::NLSMCtx *  appCtx = new nlsm::NLSMCtx(mesh); 
    ts::SDC<DendroScalar,nlsm::NLSMCtx>* sdc_solver = new ts::SDC<DendroScalar,nlsm::NLSMCtx>();
    sdc_solver->set_application_ctx(appCtx);
    sdc_solver->set_time_integration_order(1);

    for(sdc_solver->init(); sdc_solver->curr_time() < nlsm::NLSM_RK45_TIME_END ; sdc_solver->evolve())
    //for(sdc_solver->init(); sdc_solver->curr_step() < 2; sdc_solver->evolve())
    {
        const DendroIntL   step = sdc_solver->curr_step();
        const DendroScalar time = sdc_solver->curr_time();    

        const bool isActive = sdc_solver->is_active();
        const unsigned int rank_global = sdc_solver->get_global_rank();

        if(!rank_global)
          std::cout<<"[ETS] : Executing step :  "<<sdc_solver->curr_step()<<"\tcurrent time :"<<sdc_solver->curr_time()<<"\t dt:"<<sdc_solver->ts_size()<<"\t"<<std::endl;

        appCtx->terminal_output();  

        bool isRemesh = false;    
        
        if( (step % nlsm::NLSM_REMESH_TEST_FREQ) == 0 )
            isRemesh = appCtx->is_remesh();

        if(isRemesh)
        {
            if(!rank_global)
                std::cout<<"[ETS] : Remesh is triggered.  \n";

            appCtx->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ, nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
            appCtx->terminal_output();

        }
        
        sdc_solver->sync_with_mesh();

        if((step % nlsm::NLSM_IO_OUTPUT_FREQ) == 0 )
        appCtx -> write_vtu();   

        if( (step % nlsm::NLSM_CHECKPT_FREQ) == 0 )
        appCtx -> write_checkpt();
        
        
    }

    delete appCtx->get_mesh();    
    delete appCtx;
    delete sdc_solver;

    MPI_Finalize();

    return 0;
}
