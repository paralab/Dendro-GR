/**
 * @file em2_lts.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief 
 * @version 0.1
 * @date 2020-07-18
 * 
 * @copyright Copyright (c) 2020
 * 
 */

#include "em2.h"
#include "em2Utils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "octUtils.h"
#include "em2Ctx.h"
#include "ets.h"
#include "enuts.h"
#include "assert.h"
#include "mathUtils.h"

int main (int argc, char** argv)
{


    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<"\n paramFile \n TS_MODE (LTS=0, GTS=1)"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    // TS_MODE = 0 : LTS
    // TS_MODE = 1 : GTS
    unsigned int TS_MODE=0;
    if(argc > 2)
        TS_MODE = atoi(argv[2]);

    int root = std::min(1,npes-1);

    em2::timer::initFlops();
    em2::timer::total_runtime.start();

    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    em2::readParamFile(argv[1],comm);

    em2::dumpParamFile(std::cout,root,comm);

    _InitializeHcurve(em2::EM2_DIM);
    m_uiMaxDepth=em2::EM2_MAXDEPTH;
    
    // if(em2::EM2_NUM_VARS%em2::EM2_ASYNC_COMM_K!=0)
    // {
    //     if(!rank) std::cout<<"[overlap communication error]: total EM2_NUM_VARS: "<<em2::EM2_NUM_VARS<<" is not divisable by EM2_ASYNC_COMM_K: "<<em2::EM2_ASYNC_COMM_K<<std::endl;
    //     exit(0);
    // }

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){em2::initData(x,y,z,var);};
    //std::function<void(double,double,double,double,double*)> u_x_t=[](double x,double y,double z,double t,double*var){em2::analyticalSol(x,y,z,t,var);};

    const unsigned int interpVars=em2::EM2_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<em2::EM2_NUM_VARS;i++)
        varIndex[i]=i;
    
    
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    em2::timer::t_f2o.start();

    if(em2::EM2_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(em2::EM2_BLK_MIN_X,em2::EM2_BLK_MIN_Y,em2::EM2_BLK_MIN_Z);
        const Point pt_max(em2::EM2_BLK_MAX_X,em2::EM2_BLK_MAX_Y,em2::EM2_BLK_MAX_Z);
        em2::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-(binOp::fastLog2(em2::EM2_ELE_ORDER)),m_uiMaxDepth,comm);

    }else
    {
        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,em2::EM2_NUM_VARS,em2::EM2_REFINE_VARIABLE_INDICES,em2::EM2_NUM_REFINE_VARS,tmpNodes,m_uiMaxDepth,em2::EM2_WAVELET_TOL,em2::EM2_ELE_ORDER,comm);
    }

    unsigned int lmin=1,lmax=5;
    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),em2::EM2_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,em2::EM2_DENDRO_GRAIN_SZ,em2::EM2_LOAD_IMB_TOL,em2::EM2_SPLIT_FIX);

    if(!em2::EM2_ENABLE_BLOCK_ADAPTIVITY)
        ot::meshWAMRConvergence(mesh,f_init,em2::EM2_WAVELET_TOL,em2::EM2_NUM_VARS,em2::EM2_ELE_ORDER,em2::EM2_REFINE_VARIABLE_INDICES,em2::EM2_NUM_REFINE_VARS,5);

    mesh->setDomainBounds(Point(em2::EM2_GRID_MIN_X,em2::EM2_GRID_MIN_Y,em2::EM2_GRID_MIN_Z), Point(em2::EM2_GRID_MAX_X, em2::EM2_GRID_MAX_Y,em2::EM2_GRID_MAX_Z));
    mesh->computeMinMaxLevel(lmin,lmax);
    em2::EM2_RK45_TIME_STEP_SIZE=em2::EM2_CFL_FACTOR*((em2::EM2_COMPD_MAX[0]-em2::EM2_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) em2::EM2_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    par::Mpi_Bcast(&em2::EM2_RK45_TIME_STEP_SIZE,1,0,comm);

    DendroIntL lblocks = mesh->getLocalBlockList().size();
    DendroIntL gblocks =0; 
    par::Mpi_Reduce(&lblocks,&gblocks,1,MPI_SUM,0,comm);

    if(!rank)
    {
        std::cout<<" number of blocks for coarset block level : "<<(m_uiMaxDepth-MAXDEAPTH_LEVEL_DIFF-1)<<" # blocks: "<<gblocks<<std::endl;
        std::cout<<" lmin: "<<lmin<<" lmax: "<<lmax<<std::endl;
        std::cout<<"ts_mode: "<<TS_MODE<<std::endl;

    }
    
    #if 0
        std::vector<DendroIntL> part_sz;
        std::vector<DendroIntL> weight_sz;
        part_sz.resize(npes,0);
        weight_sz.resize(npes,0);
        localSz = mesh->getNumLocalMeshElements();
        par::Mpi_Gather(&localSz,part_sz.data(),1,0,comm);
        const ot::TreeNode* pNodes = mesh->getAllElements().data();

        for(unsigned int ele = mesh->getElementLocalBegin(); ele < mesh->getElementLocalEnd(); ele ++)
        localSz+=(nlsm::getEleWeight(&pNodes[ele])/100000);
        
        par::Mpi_Gather(&localSz,weight_sz.data(),1,0,comm);

        if(!rank)
        {
        for(unsigned int i=0; i < part_sz.size(); i++)
        {
            std::cout<<"local sz: "<<i<<" part_sz: "<<part_sz[i]<<" weight: "<<weight_sz[i]<<std::endl;
        }
            
        }
    #endif

    const ts::ETSType tsType = ts::ETSType::RK3;
    if(TS_MODE == 0)
    { 
        
        em2::EM2Ctx *  appCtx = new em2::EM2Ctx(mesh); 
        ts::ExplicitNUTS<DendroScalar,em2::EM2Ctx>*  enuts = new ts::ExplicitNUTS<DendroScalar,em2::EM2Ctx>(appCtx);
    
        std::vector<double> ld_stat_g;
        enuts->set_evolve_vars(appCtx->get_evolution_vars());
        enuts->set_ets_coefficients(tsType);
        
        const unsigned int rank_global = enuts->get_global_rank();
        const unsigned int pt_remesh_freq = 100;// (1u<<(lmax-lmin-3))
        for(enuts->init(); enuts->curr_time() < em2::EM2_RK45_TIME_END ; enuts->evolve_with_remesh(pt_remesh_freq))
        {
            const DendroIntL step = enuts->curr_step();
            const DendroScalar time = enuts->curr_time();    

            const bool isActive = enuts->is_active();

            if(!rank_global)
                std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

            appCtx->terminal_output();  

            bool isRemesh = false;    
            if( (step % em2::EM2_REMESH_TEST_FREQ) == 0 )
                isRemesh = appCtx->is_remesh();
            
            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[Explicit NUTS]: Remesh triggered "<<std::endl;;

                appCtx->remesh_and_gridtransfer(em2::EM2_DENDRO_GRAIN_SZ, em2::EM2_LOAD_IMB_TOL,em2::EM2_SPLIT_FIX,true,false,false);
                //appCtx->terminal_output();
            }
                    
            enuts->sync_with_mesh();
            
            if((step % em2::EM2_IO_OUTPUT_FREQ) == 0 )
             appCtx -> write_vtu();   

            if( (step % em2::EM2_CHECKPT_FREQ) == 0 )
            appCtx -> write_checkpt();

            //appCtx_ets->dump_pt(std::cout);
            //appCtx_enuts->dump_pt(std::cout);
            //ets->dump_pt(std::cout);
            //enuts->dump_pt(std::cout);

            #ifdef __PROFILE_ETS__
              char fName[200];
              std::ofstream f_ets, f_enuts;
              sprintf(fName,"%s_enuts.prof",em2::EM2_PROFILE_FILE_PREFIX.c_str());
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

    }else if(TS_MODE==1)
    {
        // GTS code here. 
        em2::EM2Ctx *  appCtx = new em2::EM2Ctx(mesh);
        ts::ETS<DendroScalar,em2::EM2Ctx>* ets = new ts::ETS<DendroScalar,em2::EM2Ctx>(appCtx);
        ets->set_evolve_vars(appCtx->get_evolution_vars());
        ets->set_ets_coefficients(tsType);
        
        
        for(ets->init(); ets->curr_time() < em2::EM2_RK45_TIME_END ; ets->evolve())
        {
            const DendroIntL   step = ets->curr_step();
            const DendroScalar time = ets->curr_time();    

            const bool isActive = ets->is_active();
            const unsigned int rank_global = ets->get_global_rank();

            if(!rank_global)
            std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

            appCtx->terminal_output();  

            bool isRemesh = false;    
            
            if( (step % em2::EM2_REMESH_TEST_FREQ) == 0 )
                isRemesh = appCtx->is_remesh();

            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[ETS] : Remesh is triggered.  \n";

                appCtx->remesh_and_gridtransfer(em2::EM2_DENDRO_GRAIN_SZ, em2::EM2_LOAD_IMB_TOL,em2::EM2_SPLIT_FIX,true,false,false);
                //appCtx->terminal_output();

            }
            
            ets->sync_with_mesh();

            if((step % em2::EM2_IO_OUTPUT_FREQ) == 0 )
            appCtx -> write_vtu();   

            if( (step % em2::EM2_CHECKPT_FREQ) == 0 )
            appCtx -> write_checkpt();

            //appCtx_ets->dump_pt(std::cout);
            //appCtx_enuts->dump_pt(std::cout);
            //ets->dump_pt(std::cout);
            //enuts->dump_pt(std::cout);
            #ifdef __PROFILE_ETS__
              char fName[200];
              std::ofstream f_ets, f_enuts;
              sprintf(fName,"%s_ets.prof",em2::EM2_PROFILE_FILE_PREFIX.c_str());
              
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



