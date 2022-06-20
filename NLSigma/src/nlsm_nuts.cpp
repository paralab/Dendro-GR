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
#include "nlsmCtx.h"
#include "ets.h"
#include "enuts.h"
#include "assert.h"
#include "mathUtils.h"

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
    // tmpNodes.clear();
    // createRegularOctree(tmpNodes,lmin,m_uiDim,m_uiMaxDepth,comm);

    // //std::cout<<"tmpOct: "<<tmpOct.size()<<std::endl;

    // const unsigned int mid_x = (1u<<(m_uiMaxDepth-1));
    // std::vector<ot::TreeNode> rfTmpOct;
    
    // for(unsigned int l = lmin; l < lmax ; l++)
    // {
        
    //     for(unsigned int ele =0; ele < tmpOct.size(); ele++)
    //     {
    //         if(tmpOct[ele].minX() < mid_x /*&& tmpOct[ele].minY() < mid_x && tmpOct[ele].minZ() < mid_x*/ )
    //             tmpOct[ele].addChildren(rfTmpOct);
    //         else
    //             rfTmpOct.push_back(tmpOct[ele]);
            
    //     }

    //     std::swap(rfTmpOct,tmpOct);
    //     rfTmpOct.clear();

    //     std::cout<<"tmpOct: at lev "<<l<<" size : "<<tmpOct.size()<<std::endl;

    // }
    // std::swap(tmpNodes,tmpOct);

    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),nlsm::NLSM_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,nlsm::NLSM_DENDRO_GRAIN_SZ,nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
    mesh->setDomainBounds(Point(nlsm::NLSM_GRID_MIN_X,nlsm::NLSM_GRID_MIN_Y,nlsm::NLSM_GRID_MIN_Z), Point(nlsm::NLSM_GRID_MAX_X, nlsm::NLSM_GRID_MAX_Y,nlsm::NLSM_GRID_MAX_Z));
    
    bool is_mindepth_refine_g = false;
    do{
      
      if(!rank)
        std::cout<<"enforce min depth refinement currently only works for block AMR for NLSM"<<std::endl;

      bool is_mindepth_refine = false;
      std::vector<unsigned int> refine_flag;
      refine_flag.reserve(mesh->getNumLocalMeshElements());
      const ot::TreeNode* pNodes = mesh->getAllElements().data();
      for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele++)
      {
        if(pNodes[ele].getLevel() < nlsm::NLSM_MINDEPTH)
        {
          refine_flag.push_back(OCT_SPLIT);
          is_mindepth_refine=true;
        }else
        {
          refine_flag.push_back(OCT_NO_CHANGE);
        }
          
      }

      MPI_Allreduce(&is_mindepth_refine,&is_mindepth_refine_g,1,MPI_C_BOOL,MPI_LOR,comm);

      if(is_mindepth_refine_g){
        mesh->setMeshRefinementFlags(refine_flag);
        ot::Mesh* newMesh = mesh->ReMesh();

        DendroIntL localSz = mesh->getNumLocalMeshElements();
        DendroIntL gSz_new, gSz_old;

        par::Mpi_Reduce(&localSz,&gSz_old,1,MPI_SUM,0,comm);
        localSz = newMesh->getNumLocalMeshElements();
        par::Mpi_Reduce(&localSz,&gSz_new,1,MPI_SUM,0,comm);

        if(!rank)
            std::cout<<"old mesh size: "<<gSz_old<<" new mesh size: "<<gSz_new<<std::endl;

        std::swap(newMesh,mesh);
        delete newMesh;

      }

    }while(is_mindepth_refine_g);
    
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


    // std::vector<DendroIntL> part_sz;
    // std::vector<DendroIntL> weight_sz;
    // part_sz.resize(npes,0);
    // weight_sz.resize(npes,0);
    // localSz = mesh->getNumLocalMeshElements();
    // par::Mpi_Gather(&localSz,part_sz.data(),1,0,comm);
    // const ot::TreeNode* pNodes = mesh->getAllElements().data();

    // for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele ++)
    //   localSz+=(nlsm::getEleWeight(&pNodes[ele])/100000);
    
    // par::Mpi_Gather(&localSz,weight_sz.data(),1,0,comm);

    // if(!rank)
    // {
    //   for(unsigned int i=0; i < part_sz.size(); i++)
    //   {
    //     std::cout<<"local sz: "<<i<<" part_sz: "<<part_sz[i]<<" weight: "<<weight_sz[i]<<std::endl;
    //   }
        
    // }



    const ts::ETSType tsType = ts::ETSType::RK4;
    if(ts_mode == 0)
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

                nlsm::DVec& eVars = appCtx->get_evolution_vars();
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

    }else if(ts_mode==1)
    { 
        //UTS
        nlsm::NLSMCtx *  appCtx = new nlsm::NLSMCtx(mesh);
        ts::ETS<DendroScalar,nlsm::NLSMCtx>* ets = new ts::ETS<DendroScalar,nlsm::NLSMCtx>(appCtx);
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
              bool isRemesh = appCtx->is_remesh();
              if(isRemesh)
              {
                  if(!rank_global)
                      std::cout<<"[ETS] : Remesh is triggered.  \n";

                  appCtx->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ, nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
                  appCtx->terminal_output();

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
    //else if(ts_mode ==2)
    // {
    //     profiler_t t_rt;
    //     t_rt.clear();

    //     // perform a comparison test between ets and enuts. 
    //     nlsm::NLSMCtx *  appCtx_enuts = new nlsm::NLSMCtx(mesh); 
    //     nlsm::NLSMCtx *  appCtx_ets = new nlsm::NLSMCtx(mesh); 

    //     ts::ExplicitNUTS<DendroScalar,nlsm::NLSMCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,nlsm::NLSMCtx>(appCtx_enuts,false);
    //     ts::ETS<DendroScalar,nlsm::NLSMCtx>*           ets   = new ts::ETS<DendroScalar,nlsm::NLSMCtx>(appCtx_ets);


    //     ets   -> set_evolve_vars(appCtx_ets->get_evolution_vars());
    //     enuts -> set_evolve_vars(appCtx_enuts->get_evolution_vars());
        
    //     ets   -> set_ets_coefficients(tsType);
    //     enuts -> set_ets_coefficients(tsType);

        
    //     unsigned int num_steps = ( enuts->get_dt_max() / enuts->get_dt_min() );
    //     const unsigned int rank_global = ets->get_global_rank();

    //     if(!rank_global) 
    //       std::cout<<" num_steps: "<<num_steps<<std::endl;

    //     //enuts->dump_load_statistics(std::cout);
    //     enuts->init();
    //     enuts->dump_est_speedup(std::cout);

    //     nlsm::DVec& evar_enuts = appCtx_enuts->get_evolution_vars();

    //     t_rt.snapreset();

    //     t_rt.start();
    //       for(enuts->init(); enuts->curr_step() < 1 ; enuts->evolve());
    //     t_rt.stop();

    //     t_stat=t_rt.snap;
    //     nlsm::timer::computeOverallStats(&t_stat, t_stat_g, comm);

    //     if(!rank_global)
    //       std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

    //     if(!rank_global)
    //       printf("[ENUTS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);


        
    //     t_rt.snapreset();

    //     t_rt.start();

    //     for(ets->init(); ets->curr_step() < num_steps ; ets->evolve())
    //     {
    //         if(!rank_global)
    //             std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

    //     }

    //     t_rt.stop();
    //     t_stat=t_rt.snap;
    //     nlsm::timer::computeOverallStats(&t_stat, t_stat_g, comm);


    //     if(!rank_global)
    //             std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

    //     if(!rank_global)
    //       printf("[ETS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);

        
    //     nlsm::DVec& evar_ets = appCtx_ets->get_evolution_vars();

    //     nlsm::DVec& evar_diff;
    //     evar_diff.VecCopy(evar_ets,false);
               

    //     if(!rank)
    //       std::cout<<"\n\n\n"<<std::endl;

    //     for(unsigned int v=0; v < evar_diff.get_dof(); v++)
    //     {
    //       double min,max;
    //       evar_diff.VecMinMax(mesh, min, max, v);
    //       if(!rank)
    //         std::cout<<"[ETS] vec(min, max): ( "<<min<<" \t"<<", "<<max<<" ) "<<std::endl;
    //     }

    //     evar_diff.VecCopy(evar_enuts,true);

    //     if(!rank)
    //       std::cout<<"\n\n\n"<<std::endl;

    //     for(unsigned int v=0; v < evar_diff.get_dof(); v++)
    //     {
    //       double min,max;
    //       evar_diff.VecMinMax(mesh, min, max, v);
    //       if(!rank)
    //         std::cout<<"[ENUTS] vec(min, max): ( "<<min<<" \t"<<", "<<max<<" ) "<<std::endl;
    //     }


    //     evar_diff.VecFMA(mesh,evar_ets,1,-1,true);

    //     if(!rank)
    //       std::cout<<"\n\n\n"<<std::endl;


    //     for(unsigned int v=0; v < evar_diff.get_dof(); v++)
    //     {
    //       double min,max;
    //       //evar_diff.VecMinMax(mesh, min, max, v);
    //       if(!rank)
    //         std::cout<<"[diff] vec(min, max): ( "<<min<<" \t"<<", "<<max<<" ) "<<std::endl;
        
    //     }


    //     evar_diff.destroy_vector();

    //     //appCtx_ets->dump_pt(std::cout);
    //     //appCtx_enuts->dump_pt(std::cout);
    //     //ets->dump_pt(std::cout);
    //     //enuts->dump_pt(std::cout);
    //     #ifdef __PROFILE_ETS__
    //       char fName[200];
    //       std::ofstream f_ets, f_enuts;
    //       sprintf(fName,"%s_ets.prof",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
          
    //       if(!rank) 
    //       {

    //         f_ets.open (fName,std::fstream::app);
    //         if(f_ets.fail()) {std::cout<<fName<<" file open failed "<<std::endl; MPI_Abort(comm,0);}

    //       }
            

    //       sprintf(fName,"%s_enuts.prof",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
          
    //       if(!rank) 
    //       {
    //         f_enuts.open (fName,std::fstream::app);
    //         if(f_enuts.fail()) {std::cout<<fName<<" file open failed "<<std::endl; MPI_Abort(comm,0);}

    //       }

    //       ets->dump_pt(f_ets);
    //       enuts->dump_pt(f_enuts);

    //       if(!rank)  f_ets.close();
    //       if(!rank)  f_enuts.close();
    //     #endif

    //     delete appCtx_enuts->get_mesh();    
    //     delete appCtx_enuts;
    //     delete appCtx_ets;
    //     delete enuts;
    //     delete ets;

        


    // }else if(ts_mode ==3)
    // {

    //     // perform a comparison test between ets and enuts. 
    //     nlsm::NLSMCtx *  appCtx_enuts = new nlsm::NLSMCtx(mesh); 
    //     nlsm::NLSMCtx *  appCtx_ets = new nlsm::NLSMCtx(mesh); 

    //     ts::ExplicitNUTS<DendroScalar,nlsm::NLSMCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,nlsm::NLSMCtx>(appCtx_enuts,false);
    //     ts::ETS<DendroScalar,nlsm::NLSMCtx>*           ets   = new ts::ETS<DendroScalar,nlsm::NLSMCtx>(appCtx_ets);


    //     ets   -> set_evolve_vars(appCtx_ets->get_evolution_vars());
    //     enuts -> set_evolve_vars(appCtx_enuts->get_evolution_vars());
        
    //     ets   -> set_ets_coefficients(tsType);
    //     enuts -> set_ets_coefficients(tsType);

    //     unsigned int num_steps = ( enuts->get_dt_max() / enuts->get_dt_min() );
    //     const unsigned int rank_global = ets->get_global_rank();

    //     if(!rank_global) 
    //       std::cout<<" num_steps: "<<num_steps<<std::endl;


    //     enuts->init();
    //     ets->init();

    //     unsigned int enuts_step=0;
    //     DendroIntL cg_sz = mesh->getDegOfFreedom();
    //     DendroIntL cg_sz_g; 
    //     par::Mpi_Reduce(&cg_sz,&cg_sz_g,1,MPI_SUM,0,comm);
    //     char fName[200];

    //     while(enuts->curr_time() < nlsm::NLSM_RK45_TIME_END)
    //     {

    //       const DendroIntL step = enuts->curr_step();
    //       const DendroScalar time = enuts->curr_time();    
    //       const bool isActive = enuts->is_active();

    //       nlsm::DVec& evar_ets   = appCtx_ets->get_evolution_vars();
    //       nlsm::DVec& evar_enuts = appCtx_enuts->get_evolution_vars();

    //       nlsm::DVec& evar_diff;
    //       evar_diff.VecCopy(evar_ets,false);
    //       evar_diff.VecFMA(mesh,evar_enuts,1,-1,true);

    //       std::ofstream diffFile;
    //       if(!rank_global)
    //       {
    //           sprintf(fName,"%s_diff_lts_gts_.dat",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
    //           diffFile.open(fName,std::fstream::app);
              
    //           if(diffFile.fail()) {std::cout<<" diff.dat file open failed "<<std::endl;
    //            MPI_Abort(comm,0);
    //           }

    //           if(enuts_step==0){
                
    //             diffFile<<"time\tstep\t";
    //             for(unsigned int v=0; v < evar_diff.get_dof(); v++)
    //               diffFile<<nlsm::NLSM_VAR_NAMES[v]<<"\t"<<nlsm::NLSM_VAR_NAMES[v]<<"_l2\t";

    //             diffFile<<std::endl;
          
    //           }
    //       }

    //       for(unsigned int v=0; v < evar_diff.get_dof(); v++)
    //       {
    //         double min,max;
    //         evar_diff.VecMinMax(mesh, min, max, v);

    //         double l2 = normL2(evar_diff.GetVecArray()+ v*cg_sz + mesh->getNodeLocalBegin(), mesh->getNumLocalMeshNodes(), comm);
            
    //         if(!rank_global)
    //         {
    //           l2=l2/sqrt(cg_sz_g);
    //           std::cout<<YLW<<"[diff] vec(min, max): ( "<<min<<" \t"<<", "<<max<<" ) "<<NRM<<std::endl;
              
    //           if(v==0)
    //             diffFile<<enuts->curr_time()<<"\t"<<enuts_step<<"\t";
              
    //           const double el = std::max(std::fabs(max),std::fabs(min));
    //           diffFile<<el<<"\t"<<l2<<"\t";

    //           if(v==(evar_diff.get_dof()-1))
    //             diffFile<<std::endl;
    //         }

    //       }

    //       double ** eVec;
    //       evar_diff.Get2DArray(eVec,false);

    //       sprintf(fName,"%s_diff_lts_gts_%d",nlsm::NLSM_VTU_FILE_PREFIX.c_str(),enuts_step);

    //       io::vtk::mesh2vtuFine(mesh,fName, 0, NULL, NULL, nlsm::NLSM_NUM_VARS, nlsm::NLSM_VAR_NAMES, (const double**) eVec);
    //       evar_diff.destroy_vector();
    //       delete [] eVec;

    //       if(!rank_global)
    //         std::cout<<"\n ========================UTS (uniform) solver: ============================\n";

    //       unsigned int ets_sc=0;           
    //       while(ets_sc < num_steps)
    //       {
    //           const DendroIntL   step = ets->curr_step();
    //           const DendroScalar time = ets->curr_time();    

    //           const bool isActive = ets->is_active();
    //           const unsigned int rank_global = ets->get_global_rank();

    //           if(!rank_global)
    //           std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

    //           appCtx_ets->terminal_output();  

    //           bool isRemesh = false;    
              
    //           if( (step % nlsm::NLSM_REMESH_TEST_FREQ) == 0 )
    //               isRemesh = appCtx_ets->is_remesh();

    //           if(isRemesh)
    //           {
    //               if(!rank_global)
    //                   std::cout<<"[ETS] : Remesh is triggered.  \n";

    //               appCtx_ets->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ, nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX,true,false,false);
    //               appCtx_ets->terminal_output();
    //           }
              
    //           ets->sync_with_mesh();

    //           // if((step % (nlsm::NLSM_IO_OUTPUT_FREQ * num_steps)) == 0 )
    //           // appCtx_ets -> write_vtu();   

    //           // if( (step % (nlsm::NLSM_CHECKPT_FREQ * num_steps)) == 0 )
    //           // appCtx_ets -> write_checkpt();

    //           ets->evolve();
    //           ets_sc += 1;
    //       }

    //       if(!rank_global)
    //           std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;
          
    //       appCtx_ets->terminal_output();

    //       if(!rank_global)
    //         std::cout<<"\n ==================ENUTS (non-uniform) solver ===================================== \n"<<std::endl;
          
    //       appCtx_enuts->terminal_output();  

    //       bool isRemesh = false;    
    //       if( (step % nlsm::NLSM_REMESH_TEST_FREQ) == 0 )
    //           isRemesh = appCtx_enuts->is_remesh();
            
    //       if(isRemesh)
    //       {
    //           if(!rank_global)
    //               std::cout<<"[Explicit NUTS]: Remesh triggered "<<std::endl;;

    //           appCtx_enuts->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ, nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX,true,false,false);
    //           appCtx_enuts->terminal_output();

    //       }
            
        
    //       enuts->sync_with_mesh();

    //       if((step % nlsm::NLSM_IO_OUTPUT_FREQ) == 0 )
    //       appCtx_enuts -> write_vtu();   

    //       // if( (step % nlsm::NLSM_CHECKPT_FREQ) == 0 )
    //       // appCtx_enuts -> write_checkpt();

    //       enuts->evolve();

    //       if(!rank_global)
    //         std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

    //       appCtx_enuts->terminal_output();  
    //       enuts_step++;

    //       //appCtx_ets->dump_pt(std::cout);
    //       //appCtx_enuts->dump_pt(std::cout);
    //       //ets->dump_pt(std::cout);
    //       //enuts->dump_pt(std::cout);
          
    //       //char fName[200];
    //       #ifdef __PROFILE_ETS__
    //       std::ofstream f_ets, f_enuts;
    //       sprintf(fName,"%s_ets.prof",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
          
    //       if(!rank) 
    //       {

    //         f_ets.open (fName,std::fstream::app);
    //         if(f_ets.fail()) {std::cout<<fName<<" file open failed "<<std::endl; MPI_Abort(comm,0);}

    //       }
            

    //       sprintf(fName,"%s_enuts.prof",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
          
    //       if(!rank) 
    //       {
    //         f_enuts.open (fName,std::fstream::app);
    //         if(f_enuts.fail()) {std::cout<<fName<<" file open failed "<<std::endl; MPI_Abort(comm,0);}

    //       }

    //       ets->dump_pt(f_ets);
    //       enuts->dump_pt(f_enuts);

    //       if(!rank)  f_ets.close();
    //       if(!rank)  f_enuts.close();

    //       ets->reset_pt();
    //       enuts->reset_pt();
    //       #endif


    //     }


        

        

    //     delete appCtx_enuts->get_mesh();    
    //     delete appCtx_enuts;
    //     delete appCtx_ets;
    //     delete enuts;
    //     delete ets;


    // }


    MPI_Finalize();

    return 0;
}
