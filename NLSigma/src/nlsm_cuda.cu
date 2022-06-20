/**
 * @file nlsm_nuts.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief : NLSM to test with NUTS. 
 * @version 0.1
 * @date 2020-04-03
 * @copyright Copyright (c) 2020
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
#include "ets.h"
#include "enuts.h"
#include "assert.h"
#include "mathUtils.h"
#include "nlsmCtxGPU.cuh"

int main (int argc, char** argv)
{
    // 0- NUTS 1-UTS
    unsigned int ts_mode=0;

    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    if(argc>2)
        ts_mode = std::atoi(argv[2]);

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    int devicesCount;
    cudaGetDeviceCount(&devicesCount);
    if(!rank)
        printf("number of cuda devices: %d\n",devicesCount);
    cudaSetDevice(rank%devicesCount);

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
    
    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),nlsm::NLSM_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,nlsm::NLSM_DENDRO_GRAIN_SZ,nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
    mesh->setDomainBounds(Point(nlsm::NLSM_GRID_MIN_X,nlsm::NLSM_GRID_MIN_Y,nlsm::NLSM_GRID_MIN_Z), Point(nlsm::NLSM_GRID_MAX_X, nlsm::NLSM_GRID_MAX_Y,nlsm::NLSM_GRID_MAX_Z));
    
    bool is_mindepth_refine_g = false;
    // do{
      
    //   if(!rank)
    //     std::cout<<"enforce min depth refinement currently only works for block AMR for NLSM"<<std::endl;

    //   bool is_mindepth_refine = false;
    //   std::vector<unsigned int> refine_flag;
    //   refine_flag.reserve(mesh->getNumLocalMeshElements());
    //   const ot::TreeNode* pNodes = mesh->getAllElements().data();
    //   for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele++)
    //   {
    //     if(pNodes[ele].getLevel() < nlsm::NLSM_MINDEPTH)
    //     {
    //       refine_flag.push_back(OCT_SPLIT);
    //       is_mindepth_refine=true;
    //     }else
    //     {
    //       refine_flag.push_back(OCT_NO_CHANGE);
    //     }
          
    //   }

    //   MPI_Allreduce(&is_mindepth_refine,&is_mindepth_refine_g,1,MPI_C_BOOL,MPI_LOR,comm);

    //   if(is_mindepth_refine_g){
    //     mesh->setMeshRefinementFlags(refine_flag);
    //     ot::Mesh* newMesh = mesh->ReMesh();

    //     DendroIntL localSz = mesh->getNumLocalMeshElements();
    //     DendroIntL gSz_new, gSz_old;

    //     par::Mpi_Reduce(&localSz,&gSz_old,1,MPI_SUM,0,comm);
    //     localSz = newMesh->getNumLocalMeshElements();
    //     par::Mpi_Reduce(&localSz,&gSz_new,1,MPI_SUM,0,comm);

    //     if(!rank)
    //         std::cout<<"old mesh size: "<<gSz_old<<" new mesh size: "<<gSz_new<<std::endl;

    //     std::swap(newMesh,mesh);
    //     delete newMesh;

    //   }

    // }while(is_mindepth_refine_g);
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);
    nlsm::NLSM_RK45_TIME_STEP_SIZE=nlsm::NLSM_CFL_FACTOR*((nlsm::NLSM_COMPD_MAX[0]-nlsm::NLSM_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) nlsm::NLSM_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    par::Mpi_Bcast(&nlsm::NLSM_RK45_TIME_STEP_SIZE,1,0,comm);

    DendroIntL lblocks = mesh->getLocalBlockList().size();
    DendroIntL gblocks =0; 
    par::Mpi_Reduce(&lblocks,&gblocks,1,MPI_SUM,0,comm);
    if(!rank)
      std::cout<<" number of blocks for coarset block level : "<<(m_uiMaxDepth-MAXDEAPTH_LEVEL_DIFF-1)<<" # blocks: "<<gblocks<<std::endl;
    
    if(!rank)
        std::cout<<" lmin: "<<lmin<<" lmax: "<<lmax<<std::endl;

     if(!rank)
      std::cout<<"ts_mode: "<<ts_mode<<std::endl;


    const ts::ETSType tsType = ts::ETSType::RK4;
    /*if(ts_mode == 0)
    { 
        
        nlsm::NLSMCtx *  appCtx = new nlsm::NLSMCtx(mesh); 
        ts::ExplicitNUTS<DendroScalar,nlsm::NLSMCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,nlsm::NLSMCtx>(appCtx);

        std::vector<double> ld_stat_g;
        enuts->set_evolve_vars(appCtx->get_evolution_vars());
        enuts->set_ets_coefficients(tsType);
        
        const unsigned int rank_global = enuts->get_global_rank();
        const unsigned int pt_remesh_freq = 5;//(1u<<(lmax-lmin-3))
        
        for(enuts->init(); enuts->curr_time() < nlsm::NLSM_RK45_TIME_END ; enuts->evolve())
        //for(enuts->init(); enuts->curr_time() < nlsm::NLSM_RK45_TIME_END ; enuts->evolve_with_remesh(pt_remesh_freq))
        {
            const DendroIntL step = enuts->curr_step();
            const DendroScalar time = enuts->curr_time();    

            const bool isActive = enuts->is_active();
            enuts->dump_load_statistics(std::cout);

            if(!rank_global)
                std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

            appCtx->terminal_output();  

            bool isRemesh = false;    
            if( (step % nlsm::NLSM_REMESH_TEST_FREQ) == 0 )
                isRemesh = appCtx->is_remesh();
            
            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[Explicit NUTS]: Remesh triggered "<<std::endl;;

                DVec eVars = appCtx->get_evolution_vars();
                appCtx->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ, nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
                //appCtx->terminal_output();
                enuts->sync_with_mesh();

            }
                    
            if((step % nlsm::NLSM_IO_OUTPUT_FREQ) == 0 )
             appCtx -> write_vtu();   

            if( (step % nlsm::NLSM_CHECKPT_FREQ) == 0 )
            appCtx -> write_checkpt();

            //appCtx_ets->dump_pt(std::cout);
            //appCtx_enuts->dump_pt(std::cout);
            //ets->dump_pt(std::cout);
            //enuts->dump_pt(std::cout);

            #ifdef __PROFILE_ETS__
              char fName[200];
              std::ofstream f_ets, f_enuts;
              sprintf(fName,"%s_enuts.prof",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
              if(!rank) 
              {
                f_enuts.open (fName,std::fstream::app);
                if(f_enuts.fail()) {std::cout<<fName<<" file open failed "<<std::endl; MPI_Abort(comm,0);}
              }

              enuts->dump_pt(f_enuts);
              enuts->reset_pt();


              if(!rank)  f_ets.close();
              if(!rank)  f_enuts.close();
            #endif

            
        }

        delete appCtx->get_mesh();    
        delete appCtx;

        delete enuts;

    }else */
    if(ts_mode==1)
    { 
        //UTS
        //nlsm::NLSMCtx *  appCtx = new nlsm::NLSMCtx(mesh);
        nlsm::NLSMCtxGPU * appCtx = new nlsm::NLSMCtxGPU(mesh);
        
        ts::ETS<DendroScalar,nlsm::NLSMCtxGPU>* ets = new ts::ETS<DendroScalar,nlsm::NLSMCtxGPU>(appCtx);
        ets->set_evolve_vars(appCtx->get_evolution_vars());
        ets->set_ets_coefficients(tsType);
        
        
        for(ets->init(); ets->curr_time() < nlsm::NLSM_RK45_TIME_END ; ets->evolve())
        {
            const DendroIntL   step = ets->curr_step();
            const DendroScalar time = ets->curr_time();    

            const bool isActive = ets->is_active();
            const unsigned int rank_global = ets->get_global_rank();

            if( (step % nlsm::NLSM_REMESH_TEST_FREQ) == 0 )
            {
              appCtx->device_to_host_sync();
              bool isRemesh = appCtx->is_remesh();
              appCtx->terminal_output();
              if(isRemesh)
              {
                  if(!rank_global)
                      std::cout<<"[ETS] : Remesh is triggered.  \n";

                  appCtx->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ, nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
                  ets->sync_with_mesh();    

              }

              if((step % nlsm::NLSM_IO_OUTPUT_FREQ) == 0 )
                appCtx -> write_vtu();   

              if( (step % nlsm::NLSM_CHECKPT_FREQ) == 0 )
                appCtx -> write_checkpt();

            }

            //appCtx_ets->dump_pt(std::cout);
            //appCtx_enuts->dump_pt(std::cout);
            //ets->dump_pt(std::cout);
            //enuts->dump_pt(std::cout);
            #ifdef __PROFILE_ETS__
              char fName[200];
              std::ofstream f_ets, f_enuts;
              sprintf(fName,"%s_ets.prof",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
              
              if(!rank) 
              {

                f_ets.open (fName,std::fstream::app);
                if(f_ets.fail()) {std::cout<<fName<<" file open failed "<<std::endl; MPI_Abort(comm,0);}

              }
          
              ets->dump_pt(f_ets);
              ets->reset_pt();
          
              if(!rank)  f_ets.close();
              if(!rank)  f_enuts.close();
            #endif



            
        }

        delete appCtx->get_mesh();    
        delete appCtx;
        delete ets;

    }
    
    MPI_Finalize();

    return 0;
}
