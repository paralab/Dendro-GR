/**
 * @file gr_scaling.cpp
 * @brief BSSN_GR scaling tests. 
 * @version 0.1
 * @date 2021-07-25
 * 
 * @copyright Copyright (c) 2021
 * 
 */

#include "gr.h"
#include "grUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rkBSSN.h"
#include "octUtils.h"
#include "meshUtils.h"
#include "mathUtils.h"
#include <fstream>
#include <ctime>
#include "bssnCtxGPU.cuh"

int bssn_driver(MPI_Comm comm, unsigned int num_step,unsigned int warm_up, std::ostream& outfile,unsigned int ts_mode)
{
  int rank, npes;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &npes);

  std::vector<ot::TreeNode> tmpNodes;
  std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::punctureData(x,y,z,var);};
  const unsigned int interpVars=bssn::BSSN_NUM_VARS;

  unsigned int varIndex[interpVars];
  for(unsigned int i=0;i<bssn::BSSN_NUM_VARS;i++)
      varIndex[i]=i;

  if(false && bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
  {
      if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
      const Point pt_min(bssn::BSSN_BLK_MIN_X,bssn::BSSN_BLK_MIN_Y,bssn::BSSN_BLK_MIN_Z);
      const Point pt_max(bssn::BSSN_BLK_MAX_X,bssn::BSSN_BLK_MAX_Y,bssn::BSSN_BLK_MAX_Z);
      bssn::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);      
  }else
  {

      if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
      const unsigned int f2olmin = std::min(bssn::BSSN_BH1_MAX_LEV,bssn::BSSN_BH2_MAX_LEV);
      if(f2olmin < MAXDEAPTH_LEVEL_DIFF + 2)
      {
        if(!rank)
          std::cout<<"BH min level should be larger than "<<(MAXDEAPTH_LEVEL_DIFF+2)<<std::endl;

        MPI_Abort(comm,0);
        
      }
      function2Octree(f_init,bssn::BSSN_NUM_VARS,varIndex,interpVars,tmpNodes,(f2olmin-MAXDEAPTH_LEVEL_DIFF-2),bssn::BSSN_WAVELET_TOL,bssn::BSSN_ELE_ORDER,comm);
      
  }

  //std::vector<ot::TreeNode> f2Octants(tmpNodes);
  ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
  mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
  unsigned int lmin, lmax;
  mesh->computeMinMaxLevel(lmin,lmax);
  tmpNodes.clear();

  if(!rank)
  {
    std::cout<<"================= Grid Info (Before init grid converge):======================================================="<<std::endl;
    std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
    std::cout<<"dx: "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
    std::cout<<"dt: "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
    std::cout<<"ts mode: "<<ts_mode<<std::endl;
    std::cout<<"==============================================================================================================="<<std::endl;
  }

  bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
  tmpNodes.clear();

  // enable block adaptivity disable remesh. 
  bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY=1;
  ot::Mesh* pMesh = mesh;
  
  if(bssn::BSSN_RESTORE_SOLVER==0)
    pMesh = bssn::weakScalingReMesh(mesh,npes);
  
  if(ts_mode == 0)
  {

  }else if(ts_mode == 1)
  {
    
    bssn::BSSNCtxGPU *  bssnCtx = new bssn::BSSNCtxGPU(pMesh);
    ts::ETS<DendroScalar,bssn::BSSNCtxGPU>* ets = new ts::ETS<DendroScalar,bssn::BSSNCtxGPU>(bssnCtx);
    ets->set_evolve_vars(bssnCtx->get_evolution_vars());
    
    if((RKType)bssn::BSSN_RK_TYPE == RKType::RK3)
        ets->set_ets_coefficients(ts::ETSType::RK3);
    else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK4)
        ets->set_ets_coefficients(ts::ETSType::RK4);
    else if((RKType)bssn::BSSN_RK_TYPE == RKType::RK45)
        ets->set_ets_coefficients(ts::ETSType::RK5);

    ets->init();

    #if defined __PROFILE_CTX__ && defined __PROFILE_ETS__
        std::ofstream outfile;
        char fname [256];
        sprintf(fname, "bssnCtxGPU_WS_%d.txt",npes);
        if(!rank)
        {
          outfile.open(fname, std::ios_base::app);
          time_t now = time(0);
          // convert now to string form
          char* dt = ctime(&now);
          outfile <<"============================================================"<<std::endl;
          outfile << "Current time : "<<dt<<" --- "<<std::endl;
          outfile <<"============================================================"<<std::endl;
        }

        ets->init_pt();
        bssnCtx->reset_pt();
        ets->dump_pt(outfile);
        //bssnCtx->dump_pt(outfile);
      #endif
      
      cudaStream_t s_gw;
      cudaStreamCreate(&s_gw);

      ts::TSInfo ts_gw_output;
      ts::TSInfo ts_curr; 
      bool is_gw_written=false;

      bool is_merge_executed =false;
      double t1 = MPI_Wtime();
      while(ets->curr_time() < bssn::BSSN_RK_TIME_END)
      {
        const DendroIntL   step = ets->curr_step();
        const DendroScalar time = ets->curr_time();

        bssn::BSSN_CURRENT_RK_COORD_TIME = time;
        bssn::BSSN_CURRENT_RK_STEP = step; 

        // if(time < 200)
        //   bssn::BSSN_REFINEMENT_MODE = RefinementMode::BH_LOC;
        // else
	      //   bssn::BSSN_REFINEMENT_MODE = RefinementMode::WAMR;
	      
        const bool isActive = ets->is_active();
        const unsigned int rank_global = ets->get_global_rank();

        const bool is_merged = bssnCtx->is_bh_merged(0.1);
        if(is_merged && !is_merge_executed)
        {
          bssn::BSSN_REMESH_TEST_FREQ = 3 * bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;
          bssn::BSSN_MINDEPTH=5;
          bssn::BSSN_GW_EXTRACT_FREQ  = bssn::BSSN_GW_EXTRACT_FREQ_AFTER_MERGER;
          bssn::BSSN_REFINEMENT_MODE = RefinementMode::WAMR;
        }/*else
        {
          //bssn::BSSN_REFINEMENT_MODE = RefinementMode::BH_LOC;
        }*/

        /*if((step % bssn::BSSN_GW_EXTRACT_FREQ) == 0 )
        {
          if(!rank_global)
            std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;
          // 03/23/22 : this works only because we use defult stream to for the evolution. Default stream is special and synchronous with all other streams.
          bssnCtx->device_to_host_async(s_gw);
          ts_gw_output=bssnCtx->get_ts_info();
          is_gw_written=false;
          
          if( (step % bssn::BSSN_REMESH_TEST_FREQ) == 0 )
          {
            cudaStreamSynchronize(s_gw);
            bool isRemesh = bssnCtx->is_remesh();
            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[ETS] : Remesh is triggered.  \n";

                bssnCtx->remesh_and_gridtransfer(bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
                bssn::deallocate_bssn_deriv_workspace();
                bssn::allocate_bssn_deriv_workspace(bssnCtx->get_mesh(),1);
                ets->sync_with_mesh();
                // correct timestep size
                ot::Mesh* pmesh = bssnCtx->get_mesh();
                unsigned int lmin, lmax;
                pmesh->computeMinMaxLevel(lmin,lmax);
                if(!pmesh->getMPIRank())
                    printf("post merger grid level = (%d, %d)\n",lmin,lmax);
                bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
                ts::TSInfo ts_in = bssnCtx->get_ts_info();
                ts_in._m_uiTh = bssn::BSSN_RK45_TIME_STEP_SIZE;
                bssnCtx->set_ts_info(ts_in);
            }
          }

        }

        if((step % bssn::BSSN_GW_EXTRACT_FREQ) == (bssn::BSSN_GW_EXTRACT_FREQ-1))
          cudaStreamSynchronize(s_gw);
          
        if((!is_gw_written) && (cudaStreamQuery(s_gw) == cudaSuccess))
        {
          ts_curr = bssnCtx->get_ts_info();
          bssnCtx->set_ts_info(ts_gw_output);
          bssnCtx->terminal_output();  
          bssnCtx->write_vtu();
          bssnCtx->evolve_bh_loc(bssnCtx->get_evolution_vars_cpu(),ets->ts_size()*bssn::BSSN_GW_EXTRACT_FREQ);
          
          if( (step % bssn::BSSN_CHECKPT_FREQ) == 0 )
            bssnCtx->write_checkpt();
          
          bssnCtx->set_ts_info(ts_curr);
          is_gw_written=true;
          
        }*/

        ets->evolve();
      }

      #if defined __PROFILE_CTX__ && defined __PROFILE_ETS__
        ets->dump_pt(outfile);
        //bssnCtx->dump_pt(outfile);
      #endif

      double t2 = MPI_Wtime()-t1;
      double t2_g;
      par::Mpi_Allreduce(&t2,&t2_g,1,MPI_MAX,ets->get_global_comm());
      if(!(ets->get_global_rank()))
        std::cout<<" ETS time (max) : "<<t2_g<<std::endl;

      if(bssn::BSSN_RESTORE_SOLVER==0)
      {
        if(pMesh == mesh)
          delete mesh;
        else
        {
          delete mesh;
          delete pMesh;
        }
      }else
      {
        delete pMesh;
      }   
      delete bssnCtx;
      delete ets;

    

  }else
  {

    if(!rank)
      RAISE_ERROR("invalid ts mode : "<<ts_mode<<"specifed");

    MPI_Abort(comm,0);
  }

  return 0;

}



int main (int argc, char** argv)
{
    // 0- NUTS 1-UTS
    unsigned int ts_mode=0;     
    
    if(argc<2)
    {
        std::cout<<"Usage: "<<argv[0]<<"paramFile TSMode(0){0-Spatially Adaptive Time Stepping(SATS, "<<GRN<<"default"<<NRM<<") , 1- Uniform Time Stepping.  }"<<std::endl;
        return 0;
    }
        
    if(argc>2)
        ts_mode = std::atoi(argv[2]);


    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);


    // Print out CMAKE options
    if (!rank) {
        #ifdef BSSN_COMPUTE_CONSTRAINTS
          std::cout<<GRN<<"  Compiled with BSSN_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT
          std::cout<<GRN<<"  Compiled with BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ENABLE_VTU_OUTPUT
          std::cout<<GRN<<"  Compiled with BSSN_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without BSSN_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_ETA_FUNCTION 
          std::cout<<GRN<<"  Compiled with  BSSN_ETA_FUNCTION"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_ETA_FUNCTION"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_BH_LOCATIONS 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_GAUGE_ROCHESTER 
          std::cout<<GRN<<"  Compiled with  BSSN_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #endif
        #ifdef BSSN_KERR_SCHILD_TEST 
          std::cout<<GRN<<"  Compiled with  BSSN_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #endif

        #ifdef BSSN_REFINE_BASE_EH 
          std::cout<<GRN<<"  Compiled with  BSSN_REFINE_BASE_EH"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  BSSN_REFINE_BASE_EH"<<NRM<<std::endl;
        #endif

        #ifdef USE_FD_INTERP_FOR_UNZIP 
          std::cout<<GRN<<"  Compiled with  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #endif

    }

    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    bssn::readParamFile(argv[1],comm);


    int root = std::min(1,npes-1);
    bssn::dumpParamFile(std::cout,root,comm);

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth=bssn::BSSN_MAXDEPTH;
    
    if(bssn::BSSN_NUM_VARS%bssn::BSSN_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total BSSN_NUM_VARS: "<<bssn::BSSN_NUM_VARS<<" is not divisable by BSSN_ASYNC_COMM_K: "<<bssn::BSSN_ASYNC_COMM_K<<std::endl;
        MPI_Abort(comm,0);
    }

    if(bssn::BSSN_GW_EXTRACT_FREQ> bssn::BSSN_IO_OUTPUT_FREQ)
    {
      if(!rank) std::cout<<" BSSN_GW_EXTRACT_FREQ  should be less BSSN_IO_OUTPUT_FREQ "<<std::endl;
      MPI_Abort(comm,0);
    }


    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::punctureData(x,y,z,var);};
    std::function<double(double,double,double)> f_init_alpha=[](double x,double y,double z){ double var[24]; bssn::punctureData(x,y,z,var); return var[0];};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<bssn::BSSN_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=bssn::VAR::U_ALPHA;
    varIndex[1]=bssn::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    const unsigned int NUM_WARM_UP=2;
    const unsigned int NUM_STEPS  =1;

    std::ofstream outfile;
    char fname [256];
    sprintf(fname, "bssnCtx_WS_%d.txt",npes); 
    
    if(!rank)
    {
      outfile.open(fname, std::ios_base::app);
      time_t now = time(0);
      // convert now to string form
      char* dt = ctime(&now);
      outfile <<"============================================================"<<std::endl;
      outfile << "Current time : "<<dt<<" --- "<<std::endl;
      outfile <<"============================================================"<<std::endl;
    }

    
    bssn_driver(comm,NUM_STEPS,NUM_WARM_UP,outfile,1);

    if(!rank)
      outfile.close();



    #ifdef RUN_WEAK_SCALING
    
      if(!rank) std::cout<<"======================================================================"<<std::endl;
      if(!rank) std::cout<<"     Weak Scaling Run Begin.   "<<std::endl;
      if(!rank) std::cout<<"======================================================================"<<std::endl;


      int proc_group = 0;
      int min_np = 2;
      for (int i = npes; rank < i && i >= min_np; i = i >> 1) proc_group++;
      MPI_Comm comm_ws;

      MPI_Comm_split(comm, proc_group, rank, &comm_ws);

      MPI_Comm_rank(comm_ws, &rank);
      MPI_Comm_size(comm_ws, &npes);

      if(!rank) outfile.open(fname, std::ios_base::app);
      MPI_Barrier(comm_ws);

      bssn_driver(comm_ws,NUM_STEPS,NUM_WARM_UP,outfile,1);
      
      MPI_Barrier(comm_ws);
      if(!rank) outfile.close();


      MPI_Comm_rank(comm, &rank);
      MPI_Comm_size(comm, &npes);
      MPI_Barrier(comm);

      if(!rank) std::cout<<"======================================================================"<<std::endl;
      if(!rank) std::cout<<"     Weak Scaling Run Complete.   "<<std::endl;
      if(!rank) std::cout<<"======================================================================"<<std::endl;

    #endif

    
    
    MPI_Finalize();
    return 0;

}
