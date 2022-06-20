//
// Created by milinda on 11/5/18.
//

/**
 * @brief This file contains the tests to measure the efficiency of the zip/unzip computations.
 * */


#include "zipUnzipTest.h"

int main (int argc, char** argv)
{

    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile regGrid"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    bssn::timer::initFlops();

    bssn::timer::total_runtime.start();


    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    bssn::readParamFile(argv[1],comm);

    double regRange=10;
    if(argc>2)
        regRange=atof(argv[2]);




    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;

        std::cout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_DIM :"<<bssn::BSSN_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_CFL_FACTOR :"<<bssn::BSSN_CFL_FACTOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_IO_OUTPUT_FREQ :"<<bssn::BSSN_IO_OUTPUT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_REMESH_TEST_FREQ :"<<bssn::BSSN_REMESH_TEST_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_CHECKPT_FREQ :"<<bssn::BSSN_CHECKPT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RESTORE_SOLVER :"<<bssn::BSSN_RESTORE_SOLVER<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ENABLE_BLOCK_ADAPTIVITY :"<<bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_VTU_FILE_PREFIX :"<<bssn::BSSN_VTU_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_CHKPT_FILE_PREFIX :"<<bssn::BSSN_CHKPT_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_PROFILE_FILE_PREFIX :"<<bssn::BSSN_PROFILE_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_IO_OUTPUT_GAP :"<<bssn::BSSN_IO_OUTPUT_GAP<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_DENDRO_GRAIN_SZ :"<<bssn::BSSN_DENDRO_GRAIN_SZ<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ASYNC_COMM_K :"<<bssn::BSSN_ASYNC_COMM_K<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_DENDRO_AMR_FAC :"<<bssn::BSSN_DENDRO_AMR_FAC<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_USE_WAVELET_TOL_FUNCTION :"<<bssn::BSSN_USE_WAVELET_TOL_FUNCTION<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_WAVELET_TOL :"<<bssn::BSSN_WAVELET_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_WAVELET_TOL_MAX:"<<bssn::BSSN_WAVELET_TOL_MAX<<NRM<<std::endl;
        std::cout<<YLW<<"\t:BSSN_WAVELET_TOL_FUNCTION_R0: "<<bssn::BSSN_WAVELET_TOL_FUNCTION_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\t:BSSN_WAVELET_TOL_FUNCTION_R1: "<<bssn::BSSN_WAVELET_TOL_FUNCTION_R1<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_LOAD_IMB_TOL :"<<bssn::BSSN_LOAD_IMB_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK_TIME_BEGIN :"<<bssn::BSSN_RK_TIME_BEGIN<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK_TIME_END :"<<bssn::BSSN_RK_TIME_END<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK45_TIME_STEP_SIZE :"<<bssn::BSSN_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_RK45_DESIRED_TOL :"<<bssn::BSSN_RK45_DESIRED_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_COMPD_MIN : ( :"<<bssn::BSSN_COMPD_MIN[0]<<" ,"<<bssn::BSSN_COMPD_MIN[1]<<","<<bssn::BSSN_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_COMPD_MAX : ( :"<<bssn::BSSN_COMPD_MAX[0]<<" ,"<<bssn::BSSN_COMPD_MAX[1]<<","<<bssn::BSSN_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_BLK_MIN : ( :"<<bssn::BSSN_BLK_MIN_X<<" ,"<<bssn::BSSN_BLK_MIN_Y<<","<<bssn::BSSN_BLK_MIN_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_BLK_MAX : ( :"<<bssn::BSSN_BLK_MAX_X<<" ,"<<bssn::BSSN_BLK_MAX_Y<<","<<bssn::BSSN_BLK_MAX_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_OCTREE_MIN : ( :"<<bssn::BSSN_OCTREE_MIN[0]<<" ,"<<bssn::BSSN_OCTREE_MIN[1]<<","<<bssn::BSSN_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_OCTREE_MAX : ( :"<<bssn::BSSN_OCTREE_MAX[0]<<" ,"<<bssn::BSSN_OCTREE_MAX[1]<<","<<bssn::BSSN_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_CONST :"<<bssn::ETA_CONST<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_R0 :"<<bssn::ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING :"<<bssn::ETA_DAMPING<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING_EXP :"<<bssn::ETA_DAMPING_EXP<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ETA_R0 :"<<bssn::BSSN_ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_ETA_POWER : ("<<bssn::BSSN_ETA_POWER[0]<<" ,"<<bssn::BSSN_ETA_POWER[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_LAMBDA : ("<<bssn::BSSN_LAMBDA[0]<<" ,"<<bssn::BSSN_LAMBDA[1]<<","<<bssn::BSSN_LAMBDA[2]<<bssn::BSSN_LAMBDA[3]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_LAMBDA_F : ("<<bssn::BSSN_LAMBDA_F[0]<<" ,"<<bssn::BSSN_LAMBDA_F[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCHI_FLOOR :"<<bssn::CHI_FLOOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_TRK0 :"<<bssn::BSSN_TRK0<<NRM<<std::endl;
        std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<bssn::KO_DISS_SIGMA<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH1 MASS :"<<bssn::BH1.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 POSITION (x,y,z) : ("<<bssn::BH1.getBHCoordX()<<", "<<bssn::BH1.getBHCoordY()<<", "<<bssn::BH1.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 VELOCITY (x,y,z) : ("<<bssn::BH1.getVx()<<", "<<bssn::BH1.getVy()<<", "<<bssn::BH1.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 SPIN (||,theta,phi): ( "<<bssn::BH1.getBHSpin()<<", "<<bssn::BH1.getBHSpinTheta()<<", "<<bssn::BH1.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH2 MASS :"<<bssn::BH2.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 POSITION (x,y,z) : ("<<bssn::BH2.getBHCoordX()<<", "<<bssn::BH2.getBHCoordY()<<", "<<bssn::BH2.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 VELOCITY (x,y,z) : ("<<bssn::BH2.getVx()<<", "<<bssn::BH2.getVy()<<", "<<bssn::BH2.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 SPIN (||,theta,phi): ( "<<bssn::BH2.getBHSpin()<<", "<<bssn::BH2.getBHSpinTheta()<<", "<<bssn::BH2.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_DIM :"<<bssn::BSSN_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_MAXDEPTH :"<<bssn::BSSN_MAXDEPTH<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_NUM_REFINE_VARS :"<<bssn::BSSN_NUM_REFINE_VARS<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_REFINE_VARIABLE_INDICES :[";
        for(unsigned int i=0;i<bssn::BSSN_NUM_REFINE_VARS-1;i++)
            std::cout<<bssn::BSSN_REFINE_VARIABLE_INDICES[i]<<", ";
        std::cout<<bssn::BSSN_REFINE_VARIABLE_INDICES[bssn::BSSN_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_NUM_EVOL_VARS_VTU_OUTPUT :"<<bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_VTU_OUTPUT_EVOL_INDICES :[";
        for(unsigned int i=0;i<bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
            std::cout<<bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
        std::cout<<bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBSSN_NUM_CONST_VARS_VTU_OUTPUT :"<<bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tBSSN_VTU_OUTPUT_CONST_INDICES :[";
        for(unsigned int i=0;i<bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT-1;i++)
            std::cout<<bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i]<<", ";
        std::cout<<bssn::BSSN_VTU_OUTPUT_CONST_INDICES[bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;


        std::cout<<YLW<<"\tTPID_TARGET_M_PLUS :"<<TPID::target_M_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_TARGET_M_MINUS :"<<TPID::target_M_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_B :"<<TPID::par_b<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_PLUS :"<<TPID::par_P_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_MINUS :"<<TPID::par_P_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_PLUS :"<<TPID::par_S_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_MINUS :"<<TPID::par_S_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_CENTER_OFFSET :"<<TPID::center_offset<<NRM<<std::endl;

        std::cout<<YLW<<"\tTPID_INITIAL_LAPSE_PSI_EXPONENT :"<<TPID::initial_lapse_psi_exponent<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_A :"<<TPID::npoints_A<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_B :"<<TPID::npoints_B<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_PHI :"<<TPID::npoints_phi<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_GIVE_BARE_MASS :"<<TPID::give_bare_mass<<NRM<<std::endl;
        std::cout<<YLW<<"\tINITIAL_LAPSE :"<<TPID::initial_lapse<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_SOLVE_MOMENTUM_CONSTRAINT :"<<TPID::solve_momentum_constraint<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_GRID_SETUP_METHOD :"<<TPID::grid_setup_method<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_VERBOSE :"<<TPID::verbose<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_ADM_TOL :"<<TPID::adm_tol<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NEWTON_TOL :"<<TPID::Newton_tol<<NRM<<std::endl;

    }

    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth=bssn::BSSN_MAXDEPTH;

    if(bssn::BSSN_NUM_VARS%bssn::BSSN_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total BSSN_NUM_VARS: "<<bssn::BSSN_NUM_VARS<<" is not divisable by BSSN_ASYNC_COMM_K: "<<bssn::BSSN_ASYNC_COMM_K<<std::endl;
        exit(0);
    }

    bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*(bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*(1.0/(double)(1u<<bssn::BSSN_MAXDEPTH));

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){bssn::punctureData(x,y,z,var);};
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
        /*const Point pt_min(bssn::BSSN_BLK_MIN_X,bssn::BSSN_BLK_MIN_Y,bssn::BSSN_BLK_MIN_Z);
        const Point pt_max(bssn::BSSN_BLK_MAX_X,bssn::BSSN_BLK_MAX_Y,bssn::BSSN_BLK_MAX_Z);*/

        const Point pt_min(-regRange,-regRange,-regRange);
        const Point pt_max(regRange,regRange,regRange);

        bssn::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,bssn::BSSN_NUM_VARS,varIndex,interpVars,tmpNodes,m_uiMaxDepth,bssn::BSSN_WAVELET_TOL,bssn::BSSN_ELE_ORDER,comm);

    }

    bssn::timer::t_f2o.stop();

    t_stat=bssn::timer::t_f2o.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=bssn::BSSN_DENDRO_GRAIN_SZ;//DENDRO_DEFAULT_GRAIN_SZ;

    bool isActive;
    MPI_Comm commActive;
    const int p_npes_prev=binOp::getPrevHighestPowerOfTwo((globalSz/grainSz));
    const int p_npes_next=binOp::getNextHighestPowerOfTwo((globalSz/grainSz));

    int p_npes=globalSz/grainSz;
    (std::abs(p_npes_prev-p_npes)<=std::abs(p_npes_next-p_npes)) ? p_npes=p_npes_prev : p_npes=p_npes_next;

    if(p_npes>npes) p_npes=npes;
    // quick fix to enforce the npes>=2 for any given grain size.
    if(p_npes<=1 && npes>1) p_npes=2;

    if(p_npes==npes)
    {
        MPI_Comm_dup(comm,&commActive);
        isActive=true;

    }else
    {
        //isActive=(rank*grainSz<globalSz);
        isActive=isRankSelected(npes,rank,p_npes);
        par::splitComm2way(isActive,&commActive,comm);

    }

    shrinkOrExpandOctree(tmpNodes,bssn::BSSN_LOAD_IMB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

    if(!isActive)
        if(tmpNodes.size()!=0)
            std::cout<<" rank_g: "<<rank<<" isActive: "<<isActive<<" f2O octants: "<<tmpNodes.size()<<std::endl;




    std::vector<ot::TreeNode> balOct;
    localSz=0;
    if(isActive)
    {

        int rank_active,npes_active;

        MPI_Comm_size(commActive,&npes_active);
        MPI_Comm_rank(commActive,&rank_active);

        if(!rank_active) std::cout<<"[MPI_COMM_SWITCH]: "<<npes_active<<std::endl;

        ot::TreeNode root(bssn::BSSN_DIM,bssn::BSSN_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        bssn::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,bssn::BSSN_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,bssn::BSSN_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,bssn::BSSN_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,bssn::BSSN_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        bssn::timer::t_cons.stop();
        t_stat=bssn::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


        bssn::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,bssn::BSSN_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,bssn::BSSN_SPLIT_FIX,commActive);
        tmpNodes.clear();

        bssn::timer::t_bal.stop();

        t_stat=bssn::timer::t_bal.seconds;
        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        if(!rank_active) std::cout<<GRN<<" 2:1 balancing max (s): "<<t_stat_g[2]<<NRM<<std::endl;
        localSz=balOct.size();


    }
    MPI_Comm_free(&commActive);

    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    unsigned int numBalOcts=globalSz;
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;

    bssn::timer::t_mesh.start();

    ot::Mesh * mesh=new ot::Mesh(balOct,1,bssn::BSSN_ELE_ORDER,comm,true,ot::SM_TYPE::FDM,bssn::BSSN_DENDRO_GRAIN_SZ,bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);

    bssn::timer::t_mesh.stop();

    t_stat=bssn::timer::t_mesh.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=mesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
    unsigned int numCGNodes=globalSz;
    if(!rank) std::cout<<GRN<<" # of CG nodes (vertices) : "<<globalSz<<NRM<<std::endl;
    if(!rank)
    {
        std::cout<< GRN<<"Mesh generation time (max): "<<t_stat_g[2]<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;


    }

    const std::vector<ot::Block>& blkList=mesh->getLocalBlockList();
    if(!rank) std::cout<<"number of blocks : "<<blkList.size()<<std::endl;
    /*if(!rank)
    {
        for(unsigned int blk=0;blk<blkList.size();blk++)
        {
            std::cout<<"blk: "<<blk<<" blkNode: "<<blkList[blk].getBlockNode()<<" element begin: "<<blkList[blk].getLocalElementBegin()<<" element end : "<<blkList[blk].getLocalElementEnd()<<std::endl;
        }
    }*/

    unsigned int octByLev[m_uiMaxDepth];
    unsigned int octByLev_g[m_uiMaxDepth];
    double regGridPrecentage;

    computeOctreeStats(&(*(mesh->getAllElements().begin()+mesh->getElementLocalBegin())),mesh->getNumLocalMeshElements(),octByLev,octByLev_g,regGridPrecentage,mesh->getMPIGlobalCommunicator());

    if(!rank) std::cout<<"regGrid param : "<<regRange<<std::endl;
    if(!rank)std::cout<<"regGrid precentage: "<<regGridPrecentage<<std::endl;


    double ** m_uiZipIn=new double*[bssn::BSSN_NUM_VARS];
    double ** m_uiZipOut=new double*[bssn::BSSN_NUM_VARS];
    double ** m_uiUnzipIn=new double*[bssn::BSSN_NUM_VARS];
    double ** m_uiUnzipOut=new double*[bssn::BSSN_NUM_VARS];

    for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
    {
        m_uiZipIn[var]=mesh->createVector<double>();
        m_uiZipOut[var]=mesh->createVector<double>();
        m_uiUnzipIn[var]=mesh->createUnZippedVector<double>();
        m_uiUnzipOut[var]=mesh->createUnZippedVector<double>();

    }


    applyInitialConditions(mesh,m_uiZipIn);

    unsigned int NUM_ITER=10;

    for(unsigned int iter=0;iter<NUM_ITER;iter++)
    {

        bssn::timer::t_ghostEx_sync.start();

        for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
            mesh->performGhostExchange(m_uiZipIn[var]);

        bssn::timer::t_ghostEx_sync.stop();


        //if(!rank) std::cout<<" unzip computation begin "<<std::endl;

        bssn::timer::t_unzip_sync.start();
        for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
            mesh->unzip(m_uiZipIn[var],m_uiUnzipIn[var]);

        bssn::timer::t_unzip_sync.stop();

        //if(!rank) std::cout<<" unzip computation end   "<<std::endl;



        Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
        Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);


        bssn::timer::t_rhs.start();

        for(unsigned int blk=0;blk<blkList.size();blk++)
        {
            const unsigned int offset = blkList[blk].getOffset();
            unsigned int sz[3];
            double hx[3];
            double ptmin[3];
            double ptmax[3];


            sz[0] = blkList[blk].getAllocationSzX();
            sz[1] = blkList[blk].getAllocationSzY();
            sz[2] = blkList[blk].getAllocationSzZ();


            const unsigned int bflag = blkList[blk].getBlkNodeFlag();

            hx[0]=blkList[blk].computeDx(pt_min, pt_max);
            hx[1]=blkList[blk].computeDy(pt_min, pt_max);
            hx[2]=blkList[blk].computeDz(pt_min, pt_max);

            ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*hx[0];
            ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*hx[1];
            ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*hx[2];


            ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*hx[0];
            ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*hx[1];
            ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*hx[2];


            bssnrhs(m_uiUnzipOut,(const double **)m_uiUnzipIn,offset,ptmin, ptmax,sz,bflag);
        }

        bssn::timer::t_rhs.stop();


        bssn::timer::t_zip.start();
        for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
            mesh->zip(m_uiUnzipOut[var],m_uiZipOut[var]);
        bssn::timer::t_zip.stop();


    }



    for(unsigned int var=0;var<bssn::BSSN_NUM_VARS;var++)
    {
        delete [] m_uiZipIn[var];
        delete [] m_uiZipOut[var];
        delete [] m_uiUnzipIn[var];
        delete [] m_uiUnzipOut[var];
    }

    delete [] m_uiZipIn;
    delete [] m_uiZipOut;
    delete [] m_uiUnzipIn;
    delete [] m_uiUnzipOut;

    bssn::timer::total_runtime.stop();

    const char separator    = ' ';
    const int nameWidth     = 30;
    const int numWidth      = 10;

    t_stat=bssn::timer::total_runtime.seconds;
    bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

    t_stat=bssn::timer::t_ghostEx_sync.seconds;
    bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"-comm(s)";
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    double t_comm=t_stat_g[2];

    t_stat=bssn::timer::t_unzip_sync.seconds;
    bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"-unzip(s)";
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    double t_unzip=t_stat_g[2];

    t_stat=bssn::timer::t_rhs.seconds;
    bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"-rhs(s)";
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    double t_rhs=t_stat_g[2];

    t_stat=bssn::timer::t_zip.seconds;
    bssn::timer::computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator) <<"-zip(s)";
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
    if(!rank)std::cout << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    double t_zip=t_stat_g[2];

    std::cout<<"gridRatio\t maxDepth\t numElems\t cgNodes\t unzip\t rhs\t zip\t commun\t"<<std::endl;
    printf("%f\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t\n",regGridPrecentage,bssn::BSSN_MAXDEPTH,numBalOcts,numCGNodes,t_unzip,t_rhs,t_zip,t_comm);



    return 0;

}
