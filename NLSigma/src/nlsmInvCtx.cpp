/**
 * @file nlsmInvCtx.cpp
 * @brief NLSM inverse context file. 
 * @version 0.1
 * @date 2021-10-13
 * 
 * @copyright Copyright (c) 2021
 * 
 */


#include "nlsmInvCtx.h"

namespace nlsm
{

    InvNLSMCtx::InvNLSMCtx() : invp::InverseCtx<InvNLSMCtx,NLSMCtx>()
    {

        return;
        
    }

    int InvNLSMCtx::initialize_forward_ctx()
    {
        int grank, gnpes;
        MPI_Comm_rank(m_uiCommGlobal, &grank);
        MPI_Comm_size(m_uiCommGlobal, &gnpes);
        
        if(m_uiParFile == std::string(""))
            if(!grank)
            {
                std::cout<<"[InverseCtx]: parameter file is not set : "<<std::endl;
                MPI_Abort(m_uiCommGlobal,0);
            }

        if(m_uiJobLauncher == NULL)
            if(!grank)
            {
                std::cout<<"[InverseCtx]: job launcher context is not set : "<<std::endl;
                MPI_Abort(m_uiCommGlobal,0);
            }
        
        if(!grank) std::cout<<" reading parameter file :"<<m_uiParFile<<std::endl;
        nlsm::readParamFile(m_uiParFile.c_str(),m_uiCommGlobal);
        nlsm::dumpParamFile(std::cout,1,m_uiCommGlobal);

        
        m_uiCntrlVars.update_init_vars(m_uiJobLauncher,std::cout);

        nlsm::NLSM_VTU_FILE_PREFIX     = nlsm::NLSM_VTU_FILE_PREFIX + std::string("_job_id_") + std::to_string(m_uiJobLauncher->get_comm_split_index());
        nlsm::NLSM_CHKPT_FILE_PREFIX   = nlsm::NLSM_CHKPT_FILE_PREFIX + std::string("_job_id_") + std::to_string(m_uiJobLauncher->get_comm_split_index());
        nlsm::NLSM_PROFILE_FILE_PREFIX = nlsm::NLSM_PROFILE_FILE_PREFIX + std::string("_job_id_") + std::to_string(m_uiJobLauncher->get_comm_split_index());

        if(!grank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        std::vector<ot::TreeNode> tmpNodes;
        std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){nlsm::initData(x,y,z,var);};
        std::function<void(double,double,double,double,double*)> u_x_t=[](double x,double y,double z,double t,double*var){nlsm::analyticalSol(x,y,z,t,var);};
        function2Octree(f_init,nlsm::NLSM_NUM_VARS,nlsm::NLSM_REFINE_VARIABLE_INDICES,nlsm::NLSM_NUM_REFINE_VARS,tmpNodes,m_uiMaxDepth,nlsm::NLSM_WAVELET_TOL,nlsm::NLSM_ELE_ORDER,m_uiCommLocal);
        ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),nlsm::NLSM_ELE_ORDER,m_uiCommLocal,1,ot::SM_TYPE::FDM,nlsm::NLSM_DENDRO_GRAIN_SZ,nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX);
        mesh->setDomainBounds(Point(nlsm::NLSM_GRID_MIN_X,nlsm::NLSM_GRID_MIN_Y,nlsm::NLSM_GRID_MIN_Z), Point(nlsm::NLSM_GRID_MAX_X, nlsm::NLSM_GRID_MAX_Y,nlsm::NLSM_GRID_MAX_Z));
        unsigned int lmin, lmax;
        mesh->computeMinMaxLevel(lmin,lmax);

        nlsm::NLSM_RK45_TIME_STEP_SIZE=nlsm::NLSM_CFL_FACTOR*((nlsm::NLSM_COMPD_MAX[0]-nlsm::NLSM_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) nlsm::NLSM_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
        par::Mpi_Bcast(&nlsm::NLSM_RK45_TIME_STEP_SIZE,1,0,m_uiCommLocal);


        m_uiForCtx = new nlsm::NLSMCtx(mesh);
        m_uiTSObj  = new ts::ETS<DendroScalar,nlsm::NLSMCtx>(m_uiForCtx);
        m_uiTSObj->set_evolve_vars(m_uiForCtx->get_evolution_vars());
        m_uiTSObj->set_ets_coefficients(ts::ETSType::RK4);
        return 0; 

    }

    int InvNLSMCtx::launch_forward_solve()
    {
        ts::ETS<DendroScalar,NLSMCtx>* ets = m_uiTSObj;
        NLSMCtx* appCtx                    = m_uiForCtx;

        for(ets->init(); ets->curr_time() < nlsm::NLSM_RK45_TIME_END ; ets->evolve())
        {
            const DendroIntL   step = ets->curr_step();
            const DendroScalar time = ets->curr_time();    

            const bool isActive = ets->is_active();
            const unsigned int rank_global = ets->get_global_rank();

            if(step % nlsm::NLSM_TIME_STEP_OUTPUT_FREQ == 0)
            {
                if(!rank_global)
                    std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;
                appCtx->terminal_output();
                this->extract_waves(m_uiCntrlVars.WAVE_EXTRAC_LOC,m_uiCntrlVars.NUM_WAVE_EXTRACT_LOC);
            }

            bool isRemesh = false;    
            
            if( (step % nlsm::NLSM_REMESH_TEST_FREQ) == 0 )
                isRemesh = appCtx->is_remesh();

            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[ETS] : Remesh is triggered.  \n";

                appCtx->remesh_and_gridtransfer(nlsm::NLSM_DENDRO_GRAIN_SZ, nlsm::NLSM_LOAD_IMB_TOL,nlsm::NLSM_SPLIT_FIX,true,false,false);
                appCtx->terminal_output();

                ets->sync_with_mesh();

            }

            if((step % nlsm::NLSM_IO_OUTPUT_FREQ) == 0 )
                appCtx -> write_vtu();

            if( (step % nlsm::NLSM_CHECKPT_FREQ) == 0 )
                appCtx -> write_checkpt();

            
        }

        return 0;
    }
    

    int InvNLSMCtx::launch_adjoint_solve()
    {
        return 0;
    }


    int InvNLSMCtx::gradient_approximation()
    {
        return 0;
    }

    int InvNLSMCtx::extract_waves(const Point* pts, unsigned int n)
    {
        NLSMCtx* appCtx        = m_uiForCtx;
        ot::Mesh* mesh  = appCtx->get_mesh();
        if(mesh->isActive())
        {
            const unsigned int rankActive = mesh->getMPIRank();
            const unsigned int npesActive = mesh->getMPICommSize();
            const unsigned int nodes_cg   = mesh->getDegOfFreedom();

            const unsigned int numPts=n;
            std::vector<double> domain_pts;
            domain_pts.resize(3*numPts);

            for (unsigned int i=0; i < numPts; i++)
            {
                domain_pts[3 * i + 0] = pts[i].x();
                domain_pts[3 * i + 0] = pts[i].y();
                domain_pts[3 * i + 0] = pts[i].z();
            }

            

            Point grid_limits[2];
            Point domain_limits[2];

            grid_limits[0] = Point(nlsm::NLSM_OCTREE_MIN[0], nlsm::NLSM_OCTREE_MIN[1], nlsm::NLSM_OCTREE_MIN[2]);
            grid_limits[1] = Point(nlsm::NLSM_OCTREE_MAX[0], nlsm::NLSM_OCTREE_MAX[1], nlsm::NLSM_OCTREE_MAX[2]);

            domain_limits[0] = Point(nlsm::NLSM_COMPD_MIN[0], nlsm::NLSM_COMPD_MIN[1], nlsm::NLSM_COMPD_MIN[2]);
            domain_limits[1] = Point(nlsm::NLSM_COMPD_MAX[0], nlsm::NLSM_COMPD_MAX[1], nlsm::NLSM_COMPD_MAX[2]);

            std::vector<unsigned int > validIndex;
            std::vector<double> extracted_wave;
            extracted_wave.resize(m_uiCntrlVars.NUM_WAVE_EXTRACT_LOC);
            double * vin = appCtx->get_evolution_vars().GetVecArray();
            
            mesh->readFromGhostBegin<double>(vin,nlsm::NLSM_NUM_VARS);
            mesh->readFromGhostEnd<double>(vin,nlsm::NLSM_NUM_VARS);

            ot::da::interpolateToCoordsAndGather(mesh,&vin[(nlsm::VAR::U_CHI)*nodes_cg], domain_pts.data(), domain_pts.size(), grid_limits, domain_limits, &(*(extracted_wave.begin())),0,1);

            if(!(mesh->getMPIRank()))
            {
                char fname[256];
                
                for (unsigned int i=0; i < numPts; i++)
                {
                    std::ofstream fileGW;
                    sprintf(fname,"%s_obs_%d.dat",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str(),i);
                    fileGW.open(fname,std::ofstream::app);

                    if(m_uiTSObj->curr_step()==0)
                        fileGW<<"time\tchi"<<std::endl;

                    fileGW<<std::fixed << std::setprecision(10);
                    fileGW<<m_uiTSObj->curr_time()<<"\t"<<extracted_wave[i]<<std::endl;
                    fileGW.close();                    
                }
                
            }

        }

        return 0;
    }
} // end of nlsm namespace. 
