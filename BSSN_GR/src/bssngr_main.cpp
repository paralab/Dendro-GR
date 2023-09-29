/**
 * @file gr with LTS support (not enabled currently)
 * @brief Main driver for cuda BSSNSolver. 
 * @version 0.1
 * @date 2021-02-12
 * 
 */
#include "gr.h"
#include "grUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "octUtils.h"
#include "meshUtils.h"
#include "rkBSSN.h"
#include "sdc.h"
#include "bssnCtx.h"
#include "aeh.h"


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
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::initialDataFunctionWrapper(x,y,z,var);};
    std::function<double(double,double,double)> f_init_alpha=[](double x,double y,double z){ double var[24]; bssn::initialDataFunctionWrapper(x,y,z,var); return var[0];};
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

    bssn::timer::t_f2o.start();

    if(bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
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

    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
    mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);
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

    if(ts_mode==1)
    {
      bssn::BSSNCtx *  bssnCtx = new bssn::BSSNCtx(mesh);
      ts::ETS<DendroScalar,bssn::BSSNCtx>* ets = new ts::ETS<DendroScalar,bssn::BSSNCtx>(bssnCtx);
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
        sprintf(fname, "bssnCtx_%d.txt",npes);
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
      bool is_merge_executed = false;
      double t1 = MPI_Wtime();
      while(ets->curr_time() < bssn::BSSN_RK_TIME_END)
      {
        const DendroIntL   step = ets->curr_step();
        const DendroScalar time = ets->curr_time();

        bssn::BSSN_CURRENT_RK_COORD_TIME = time;
        bssn::BSSN_CURRENT_RK_STEP = step;    

        const bool isActive = ets->is_active();
        const unsigned int rank_global = ets->get_global_rank();

        const bool is_merged = bssnCtx->is_bh_merged(0.1);
        if(is_merged)
        {
          // bssn::BSSN_REMESH_TEST_FREQ=3 * bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;
          // bssn::BSSN_MINDEPTH=5;
          bssn::BSSN_REFINEMENT_MODE          = bssn::RefinementMode::WAMR;
          bssn::BSSN_USE_WAVELET_TOL_FUNCTION = 1;
          bssn::BSSN_REMESH_TEST_FREQ = bssn::BSSN_REMESH_TEST_FREQ_AFTER_MERGER;  
          bssn::BSSN_GW_EXTRACT_FREQ  = bssn::BSSN_GW_EXTRACT_FREQ_AFTER_MERGER;
        }

        if( (step % bssn::BSSN_REMESH_TEST_FREQ) == 0 )
        {
            bool isRemesh = bssnCtx->is_remesh();
            if(isRemesh)
            {
              if(!rank_global)
                  std::cout<<"[ETS] : Remesh is triggered.  \n";

              bssnCtx->remesh_and_gridtransfer(bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
              bssn::deallocate_bssn_deriv_workspace();
              bssn::allocate_bssn_deriv_workspace(bssnCtx->get_mesh(),1);
              ets->sync_with_mesh();

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

        {
          const int lmax               = std::atoi(argv[3]);
          const int ntheta             = std::atoi(argv[4]);
          const int nphi               = std::atoi(argv[5]);
          const unsigned int max_iter  = std::atoi(argv[6]);
          const double rel_tol         = std::atof(argv[7]);
          const double abs_tol         = std::atof(argv[8]);

          aeh::SpectralAEHSolver<bssn::BSSNCtx,DendroScalar> aeh_solver(bssnCtx, lmax, ntheta, nphi, false);
          const unsigned int num_lm_modes = aeh_solver.get_num_lm_modes();
          DendroScalar* h0 = new DendroScalar[num_lm_modes];
          DendroScalar* hh = new DendroScalar[num_lm_modes];
          ot::Mesh* pmesh = bssnCtx->get_mesh();
          
          {
            double r_plus  = 0.5 * sqrt(std::pow(TPID::par_m_plus , 2)  - 1 * (std::pow(TPID::par_S_plus[0],2) + std::pow(TPID::par_S_plus[1],2) + std::pow(TPID::par_S_plus[2],2)));
            double rlim[2] ={1e-3, 5};

            for (unsigned int i=0; i < num_lm_modes;i++)
            {
              h0[i]=0.0;
              hh[i]=0.0;
            }

            h0[0] = r_plus * sqrt(4 * M_PI);
            for(unsigned int ll=0; ll < lmax; ll+=2)
            {
                if (!pmesh->getMPIRank())
                  std::cout<<"sub cycle lmax = "<<ll<<std::endl;

                aeh::SpectralAEHSolver<bssn::BSSNCtx,DendroScalar> s1(bssnCtx, ll, ntheta, nphi, false, false);
                s1.solve(bssnCtx->get_bh0_loc(), bssnCtx, h0, hh, max_iter, rel_tol, abs_tol, rlim);
                std::swap(h0,hh);
            }

            char fname[256];
            sprintf(fname,"%s_%d_%d_%d_bh0_aeh.dat",bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), lmax, ntheta, nphi);
            aeh_solver.solve(bssnCtx->get_bh0_loc(), bssnCtx, h0, hh, max_iter, rel_tol, abs_tol, rlim);
            aeh_solver.aeh_to_json(bssnCtx->get_bh0_loc(), bssnCtx, hh, fname, std::ios_base::app);

          }

          {
            double r_minus = 0.5 * sqrt(std::pow(TPID::par_m_minus, 2)  - 1 * (std::pow(TPID::par_S_minus[0],2) + std::pow(TPID::par_S_minus[1],2) + std::pow(TPID::par_S_minus[2],2)));
            double rlim[2] ={1e-3, 5};
            for (unsigned int i=0; i < num_lm_modes;i++)
            {
              h0[i]=0.0;
              hh[i]=0.0;
            }

            h0[0] = r_minus * sqrt(4 * M_PI);
            
            for(unsigned int ll=0; ll < lmax; ll+=2)
            {
                if (!pmesh->getMPIRank())
                  std::cout<<"sub cycle lmax = "<<ll<<std::endl;

                aeh::SpectralAEHSolver<bssn::BSSNCtx,DendroScalar> s1(bssnCtx, ll, ntheta, nphi, false, false);
                s1.solve(bssnCtx->get_bh1_loc(), bssnCtx, h0, hh, max_iter, rel_tol, abs_tol, rlim);
                std::swap(h0,hh);
            }

            char fname[256];
            sprintf(fname,"%s_%d_%d_%d_bh1_aeh.dat",bssn::BSSN_PROFILE_FILE_PREFIX.c_str(), lmax, ntheta, nphi);
            aeh_solver.solve(bssnCtx->get_bh1_loc(), bssnCtx, h0, hh, max_iter, rel_tol, abs_tol, rlim);
            aeh_solver.aeh_to_json(bssnCtx->get_bh1_loc(), bssnCtx, hh, fname, std::ios_base::app);

          }
          

          //h0[0] = 1.0 * sqrt(4 * M_PI);
          //aeh_solver.solve(Point(0,0,0), bssnCtx, h0, hh, max_iter, 1e-8, 1e-8, 10);
          delete [] h0;
          delete [] hh;
          bssnCtx->get_mesh()->waitActive();
          MPI_Abort(bssnCtx->get_mesh()->getMPICommunicator(), 0);

        }  


        if((step % bssn::BSSN_GW_EXTRACT_FREQ) == 0 )
        {
          if(!rank_global)
            std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;
          
          bssnCtx->terminal_output();  
          bssnCtx->write_vtu();
          bssnCtx->evolve_bh_loc(bssnCtx->get_evolution_vars(),ets->ts_size()*bssn::BSSN_GW_EXTRACT_FREQ);

        }
        
        if( (step % bssn::BSSN_CHECKPT_FREQ) == 0 )
          bssnCtx->write_checkpt();
        
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

      delete bssnCtx->get_mesh();    
      delete bssnCtx;
      delete ets;

    }

    MPI_Finalize();
    return 0; 


}
