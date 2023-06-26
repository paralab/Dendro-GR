//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief 
*/
//

// NOTE: ON DEPRECIATION:
// This file was depreciated because it has the old rkBSSN solver that has been swapped out in favor of the CTX solver.
// it remains in case someone wants to use it in the future. (Though it would be best to go to a version of this repository before
// June 2023)

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

int main (int argc, char** argv)
{


    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    bssn::timer::initFlops();

    bssn::timer::total_runtime.start();

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

        #ifdef BSSN_USE_4TH_ORDER_DERIVS
          std::cout<<GRN<<"  Using 4th order FD stencils"<<NRM<<std::endl;
        #endif

        #ifdef BSSN_USE_6TH_ORDER_DERIVS
          std::cout<<GRN<<"  Using 6th order FD stencils"<<NRM<<std::endl;
        #endif

        #ifdef BSSN_USE_8TH_ORDER_DERIVS
          std::cout<<GRN<<"  Using 8th order FD stencils"<<NRM<<std::endl;
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

    ot::Mesh * mesh= ot::createMesh(tmpNodes.data(),tmpNodes.size(),bssn::BSSN_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
    mesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));

    //io::vtk::mesh2vtuFine(mesh,"begin",0,NULL,NULL,0,NULL,NULL,0,NULL,NULL,false);

    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);
    
    if(!rank)
    {
      std::cout<<"================= Grid Info (Before init grid converge):======================================================="<<std::endl;
      std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
      std::cout<<"dx: "<<((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
      std::cout<<"dt: "<<bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))))<<std::endl;
      std::cout<<"==============================================================================================================="<<std::endl;
    }    
    bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    

    ode::solver::RK_BSSN rk_bssn(mesh,bssn::BSSN_RK_TIME_BEGIN,bssn::BSSN_RK_TIME_END,bssn::BSSN_RK45_TIME_STEP_SIZE,(RKType)bssn::BSSN_RK_TYPE);
    
    if(bssn::BSSN_RESTORE_SOLVER==1)
        rk_bssn.restoreCheckPoint(bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),comm);

    bssn::timer::t_rkSolve.start();
    rk_bssn.rkSolve();
    bssn::timer::t_rkSolve.stop();
    double t2=bssn::timer::t_rkSolve.seconds;
    double t2_g;
    par::Mpi_Allreduce(&t2,&t2_g,1,MPI_MAX,comm);
    if(!rank)
      std::cout<<" rkBSSN time (max) : "<<t2_g<<std::endl;

    bssn::timer::total_runtime.stop();
    rk_bssn.freeMesh();

    
    MPI_Finalize();
    return 0;

}
