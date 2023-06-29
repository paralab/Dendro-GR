/**
 * @file bssnCtxGPU.cu
 * @brief BSSN GPU ctx file
 * @version 0.1
 * @date 2022-02-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "bssnCtxGPU.cuh"
CONST_MEM DEVICE_REAL device::refel_1d[2*REFEL_CONST_MEM_MAX];

namespace bssn
{
    BSSNCtxGPU::BSSNCtxGPU(ot::Mesh* pMesh) : Ctx()
    {
        m_uiMesh      = pMesh;

        //if(m_uiMesh)
        m_mesh_cpu  = device::MeshGPU();
        m_dptr_mesh = m_mesh_cpu.alloc_mesh_on_device(m_uiMesh);
        
        // variable allocation for evolution variables
        m_var[VL::CPU_EV].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::HOST,BSSN_NUM_VARS,true);
        m_var[VL::CPU_CV].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::HOST,BSSN_CONSTRAINT_NUM_VARS,true);
        m_var[VL::CPU_EV_UZ_IN].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,BSSN_NUM_VARS,true);
        m_var[VL::CPU_CV_UZ_IN].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,BSSN_CONSTRAINT_NUM_VARS,true);
        
        m_var[VL::GPU_EV].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_IN].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_OUT].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);

        m_uiTinfo._m_uiStep=0;
        m_uiTinfo._m_uiT = 0;
        m_uiTinfo._m_uiTb = BSSN_RK_TIME_BEGIN;
        m_uiTinfo._m_uiTe = BSSN_RK_TIME_END;
        m_uiTinfo._m_uiTh = BSSN_RK45_TIME_STEP_SIZE;     

        m_uiElementOrder = BSSN_ELE_ORDER;

        m_uiMinPt = Point(BSSN_GRID_MIN_X,BSSN_GRID_MIN_Y,BSSN_GRID_MIN_Z);
        m_uiMaxPt = Point(BSSN_GRID_MAX_X,BSSN_GRID_MAX_Y,BSSN_GRID_MAX_Z);

        // deriv vars, 
        // m_deriv_evars      = new BSSN_EVAR_DERIVS[DEVICE_RHS_NSTREAMS];
        // m_dptr_deriv_evars = (BSSN_EVAR_DERIVS*)GPUDevice::device_malloc<BSSN_EVAR_DERIVS>(DEVICE_RHS_NSTREAMS);
        
        // const unsigned int ws_szp_stream = DEVICE_MAX_DERIVS * DEVICE_RHS_BLK_SZ * DEVICE_RHS_BATCHED_GRAIN_SZ;
        // m_deriv_base       = GPUDevice::device_malloc<DEVICE_REAL>(DEVICE_RHS_NSTREAMS * ws_szp_stream);
        // for(unsigned int i=0; i < DEVICE_RHS_NSTREAMS; i++)
        // {
        //     BSSN_EVAR_DERIVS * const deriv_evars=m_deriv_evars + i;
        //     DEVICE_REAL* const  deriv_base = m_deriv_base + i* ws_szp_stream;
        //     const unsigned int BATCHED_BLOCKS_SZ = DEVICE_RHS_BATCHED_GRAIN_SZ;
        //     const unsigned int BLK_SZ = DEVICE_RHS_BLK_SZ;
        //     #include "../scripts/bssnrhs_memalloc.cuh"
        // }
        // GPUDevice::host_to_device<BSSN_EVAR_DERIVS>(m_deriv_evars,m_dptr_deriv_evars,DEVICE_RHS_NSTREAMS);
        int maxbytes = 13*13*13 *8 * 3; // 3 blocks KB
        cudaFuncSetAttribute(device::eval_rhs3<3,13,13>, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
        
        ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS, BSSN_ASYNC_COMM_K);
        ot::alloc_mpi_ctx<DendroScalar>  (m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS, BSSN_ASYNC_COMM_K);

        device::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx_device, BSSN_NUM_VARS , BSSN_ASYNC_COMM_K);
        device::alloc_mpi_ctx<DendroScalar>  (m_uiMesh, m_mpi_ctx_device, BSSN_NUM_VARS , BSSN_ASYNC_COMM_K);


        return;

    }

    BSSNCtxGPU::~BSSNCtxGPU()
    {
        for(unsigned int i=0; i < VL::END; i++)
            m_var[i].destroy_vector();

        // GPUDevice::device_free<DEVICE_REAL>(m_deriv_base);
        // GPUDevice::device_free<BSSN_EVAR_DERIVS>(m_dptr_deriv_evars);
        // delete m_deriv_evars;

        ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS, BSSN_ASYNC_COMM_K);
        device::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx_device, BSSN_NUM_VARS , BSSN_ASYNC_COMM_K);
        m_mesh_cpu.dealloc_mesh_on_device(m_dptr_mesh);
        return;
    }

    int BSSNCtxGPU::host_to_device_sync()
    {
        if(!m_uiMesh->isActive())
            return 0;
        
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::H2D].start();
        #endif
        DVec& m_evar         = m_var[VL::CPU_EV];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        GPUDevice::host_to_device(m_evar.get_vec_ptr(), m_dptr_evar.get_vec_ptr(), m_evar.get_size());

        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::H2D].stop();
        #endif

        return 0;
    }

    int BSSNCtxGPU::host_to_device_async(cudaStream_t s)
    {
        if(!m_uiMesh->isActive())
            return 0;
        
        DVec& m_evar         = m_var[VL::CPU_EV];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        GPUDevice::host_to_device_async(m_evar.get_vec_ptr(), m_dptr_evar.get_vec_ptr(), m_evar.get_size(), s);
        
        return 0;

    }

    int BSSNCtxGPU::device_to_host_sync()
    {
        if(!m_uiMesh->isActive())
            return 0;

        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::D2H].start();
        #endif

        DVec& m_evar         = m_var[VL::CPU_EV];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        GPUDevice::device_to_host(m_evar.get_vec_ptr(), m_dptr_evar.get_vec_ptr(), m_evar.get_size());
        
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::D2H].stop();
        #endif

        return 0;
        
    }

    int BSSNCtxGPU::device_to_host_async(cudaStream_t s)
    {
        if(!m_uiMesh->isActive())
            return 0;

        DVec& m_evar         = m_var[VL::CPU_EV];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        GPUDevice::device_to_host_async<DEVICE_REAL, cudaStream_t>(m_evar.get_vec_ptr(), m_dptr_evar.get_vec_ptr(), m_evar.get_size(), s);
        
        return 0;
    }

    int BSSNCtxGPU::initialize()
    {
        if(bssn::BSSN_RESTORE_SOLVER)
        {
            this->restore_checkpt();
            this->host_to_device_sync();
            return 0;
        }

        this->init_grid();
        
        bool isRefine=false;
        DendroIntL oldElements,oldElements_g;
        DendroIntL newElements,newElements_g;

        DendroIntL oldGridPoints,oldGridPoints_g;
        DendroIntL newGridPoints,newGridPoints_g;

        unsigned int iterCount=1;
        const unsigned int max_iter=bssn::BSSN_INIT_GRID_ITER;
        const unsigned int rank_global=m_uiMesh->getMPIRankGlobal();
        MPI_Comm gcomm = m_uiMesh->getMPIGlobalCommunicator();
        
        DendroScalar* unzipVar[bssn::BSSN_NUM_VARS];
        unsigned int refineVarIds[bssn::BSSN_NUM_REFINE_VARS];

        for(unsigned int vIndex=0; vIndex<bssn::BSSN_NUM_REFINE_VARS; vIndex++)
            refineVarIds[vIndex]=bssn::BSSN_REFINE_VARIABLE_INDICES[vIndex];

        double wTol=bssn::BSSN_WAVELET_TOL;
        std::function<double(double,double,double,double*hx)> waveletTolFunc =[](double x,double y, double z,double*hx) {
            return bssn::computeWTolDCoords(x,y,z,hx);
        };
        
        DVec& m_evar     = m_var[VL::CPU_EV];
        DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];

        do
        {
            
            this->unzip(m_evar,m_evar_unz,bssn::BSSN_ASYNC_COMM_K);
            m_evar_unz.to_2d(unzipVar);
            //isRefine=this->is_remesh();
            // enforce WMAR refinement based refinement initially.
            if(max_iter==0)
                isRefine = false;
            else
                isRefine = bssn::isReMeshWAMR(m_uiMesh,(const double **)unzipVar,refineVarIds,bssn::BSSN_NUM_REFINE_VARS,waveletTolFunc,bssn::BSSN_DENDRO_AMR_FAC);
            
            if(isRefine)
            {
                ot::Mesh* newMesh = this->remesh(bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,bssn::BSSN_SPLIT_FIX);
                
                oldElements = m_uiMesh->getNumLocalMeshElements();
                newElements = newMesh->getNumLocalMeshElements();

                oldGridPoints =m_uiMesh->getNumLocalMeshNodes();
                newGridPoints =newMesh->getNumLocalMeshNodes();

                par::Mpi_Allreduce(&oldElements,&oldElements_g,1,MPI_SUM,gcomm);
                par::Mpi_Allreduce(&newElements,&newElements_g,1,MPI_SUM,gcomm);

                par::Mpi_Allreduce(&oldGridPoints,&oldGridPoints_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());
                par::Mpi_Allreduce(&newGridPoints,&newGridPoints_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());
                
                if(!rank_global)
                {
                    std::cout<<"[bssnCtx] iter : "<<iterCount<<" (Remesh triggered) ->  old mesh : "<<oldElements_g<<" new mesh : "<<newElements_g<<std::endl;
                    std::cout<<"[bssnCtx] iter : "<<iterCount<<" (Remesh triggered) ->  old mesh (zip nodes) : "<<oldGridPoints_g<<" new mesh (zip nodes) : "<<newGridPoints_g<<std::endl;
                }
                  

                this->grid_transfer(newMesh);
                
                std::swap(m_uiMesh,newMesh);
                delete newMesh;

                #ifdef __CUDACC__
                    device::MeshGPU*& dptr_mesh  = this->get_meshgpu_device_ptr();
                    device::MeshGPU* mesh_gpu    = this->get_meshgpu_host_handle();

                    mesh_gpu->dealloc_mesh_on_device(dptr_mesh);
                    dptr_mesh = mesh_gpu->alloc_mesh_on_device(m_uiMesh);
                #endif



            }
            
            iterCount+=1;

        } while(isRefine && (newElements_g!=oldElements_g || newGridPoints_g!=oldGridPoints_g) && (iterCount<max_iter) );

        this->init_grid();

        // // realloc bssn deriv space
        deallocate_bssn_deriv_workspace();
        allocate_bssn_deriv_workspace(m_uiMesh,1);

        unsigned int lmin, lmax;
        m_uiMesh->computeMinMaxLevel(lmin,lmax);
        bssn::BSSN_RK45_TIME_STEP_SIZE=bssn::BSSN_CFL_FACTOR*((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
        m_uiTinfo._m_uiTh=bssn::BSSN_RK45_TIME_STEP_SIZE;

        if(!m_uiMesh->getMPIRankGlobal())
        {
            const DendroScalar dx_finest = ((bssn::BSSN_COMPD_MAX[0]-bssn::BSSN_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) bssn::BSSN_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
            const DendroScalar dt_finest = bssn::BSSN_CFL_FACTOR*dx_finest;

            std::cout<<"================= Grid Info (After init grid converge):======================================================="<<std::endl;
            std::cout<<"lmin: "<<lmin<<" lmax:"<<lmax<<std::endl;
            std::cout<<"dx: "<<dx_finest<<std::endl;
            std::cout<<"dt: "<<dt_finest<<std::endl;
            std::cout<<"==============================================================================================================="<<std::endl;
        }

        this->host_to_device_sync();
        return 0;
                    
    }

    int BSSNCtxGPU::init_grid()
    {  

        DVec& m_evar         = m_var[VL::CPU_EV];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
        const unsigned int eleOrder=m_uiMesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();


        DendroScalar* zipIn[bssn::BSSN_NUM_VARS];
        m_evar.to_2d(zipIn);

        DendroScalar var1[bssn::BSSN_NUM_VARS];

        DendroScalar mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
        // set the TP communicator. 
        if(bssn::BSSN_ID_TYPE==0)
        {
            TP_MPI_COMM=m_uiMesh->getMPIGlobalCommunicator();
            TwoPunctures((double)0,(double)0,(double)0,var1,&mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
        }
        
        for(unsigned int elem=m_uiMesh->getElementLocalBegin(); elem<m_uiMesh->getElementLocalEnd(); elem++)
        {
            DendroScalar var[bssn::BSSN_NUM_VARS];
            for(unsigned int k=0; k<(eleOrder+1); k++)
                for(unsigned int j=0; j<(eleOrder+1); j++ )
                    for(unsigned int i=0; i<(eleOrder+1); i++)
                    {
                        const unsigned int nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                        {
                            const unsigned int nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            unsigned int ownerID,ii_x,jj_y,kk_z;
                            m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            const DendroScalar len=(double)(1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                            
                            const DendroScalar x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                            const DendroScalar y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                            const DendroScalar z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));
                            
                            if (bssn::BSSN_ID_TYPE == 0) {

                                const DendroScalar xx = GRIDX_TO_X(x);
                                const DendroScalar yy = GRIDY_TO_Y(y);
                                const DendroScalar zz = GRIDZ_TO_Z(z);

                                TwoPunctures((double)xx,(double)yy,(double)zz,var,
                                            &mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
                            }
                            else {
                                // all other values are handled in the initial data wrapper including
                                // an error message
                                initialDataFunctionWrapper((double)x, (double)y, (double)z, var);
                            }
                            for(unsigned int v=0; v<bssn::BSSN_NUM_VARS; v++)
                                zipIn[v][nodeLookUp_CG]=var[v];
                        }

                    }
        
        }
        
        for(unsigned int node=m_uiMesh->getNodeLocalBegin(); node<m_uiMesh->getNodeLocalEnd(); node++)
            enforce_bssn_constraints(zipIn,node);


        #ifdef BSSN_EXTRACT_BH_LOCATIONS
            m_uiBHLoc[0]=Point(bssn::BH1.getBHCoordX(),bssn::BH1.getBHCoordY(),bssn::BH1.getBHCoordZ());
            m_uiBHLoc[1]=Point(bssn::BH2.getBHCoordX(),bssn::BH2.getBHCoordY(),bssn::BH2.getBHCoordZ());
        #endif

        return 0;

    }
    
    int BSSNCtxGPU::finalize()
    {
        return 0;    
    }
    
    int BSSNCtxGPU::write_vtu()
    {
        if(!m_uiMesh->isActive())
            return 0;

        DVec& m_evar     = m_var[VL::CPU_EV];
        DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];
        DVec& m_cvar     = m_var[VL::CPU_CV];
        DVec & m_cvar_unz= m_var[VL::CPU_CV_UZ_IN];


        this->unzip(m_evar,m_evar_unz,BSSN_ASYNC_COMM_K);
        
        DendroScalar *consUnzipVar[bssn::BSSN_CONSTRAINT_NUM_VARS];
        DendroScalar *consVar[bssn::BSSN_CONSTRAINT_NUM_VARS];

        DendroScalar *evolUnzipVar[bssn::BSSN_NUM_VARS];
        DendroScalar *evolVar[bssn::BSSN_NUM_VARS];

        m_evar_unz.to_2d(evolUnzipVar);
        m_cvar_unz.to_2d(consUnzipVar);
        
        m_evar.to_2d(evolVar);
        m_cvar.to_2d(consVar);

        #if BSSN_COMPUTE_CONSTRAINTS

            const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();

            unsigned int offset;
            double ptmin[3], ptmax[3];
            unsigned int sz[3];
            unsigned int bflag;
            double dx,dy,dz;
            const Point pt_min(bssn::BSSN_COMPD_MIN[0],bssn::BSSN_COMPD_MIN[1],bssn::BSSN_COMPD_MIN[2]);
            const Point pt_max(bssn::BSSN_COMPD_MAX[0],bssn::BSSN_COMPD_MAX[1],bssn::BSSN_COMPD_MAX[2]);
            const unsigned int PW=bssn::BSSN_PADDING_WIDTH;

            for(unsigned int blk=0; blk<blkList.size(); blk++)
            {
                offset=blkList[blk].getOffset();
                sz[0]=blkList[blk].getAllocationSzX();
                sz[1]=blkList[blk].getAllocationSzY();
                sz[2]=blkList[blk].getAllocationSzZ();

                bflag=blkList[blk].getBlkNodeFlag();

                dx=blkList[blk].computeDx(pt_min,pt_max);
                dy=blkList[blk].computeDy(pt_min,pt_max);
                dz=blkList[blk].computeDz(pt_min,pt_max);

                ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-PW*dx;
                ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-PW*dy;
                ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-PW*dz;

                ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+PW*dx;
                ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+PW*dy;
                ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+PW*dz;

                physical_constraints(consUnzipVar, (const DendroScalar **) evolUnzipVar, offset, ptmin, ptmax, sz, bflag);
            }

            /*double consVecMin[bssn::BSSN_CONSTRAINT_NUM_VARS];
            double consVecMax[bssn::BSSN_CONSTRAINT_NUM_VARS];*/
            double constraintMaskedL2[bssn::BSSN_CONSTRAINT_NUM_VARS];
            this->zip(m_cvar_unz,m_cvar);
            m_uiMesh->readFromGhostBegin(m_cvar.get_vec_ptr(),m_cvar.get_dof());
            m_uiMesh->readFromGhostEnd(m_cvar.get_vec_ptr(),m_cvar.get_dof());
            
            bssn::extractConstraints(m_uiMesh,(const DendroScalar **)consVar,evolVar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,m_uiTinfo._m_uiStep,m_uiTinfo._m_uiT);
            #ifndef BSSN_KERR_SCHILD_TEST
                #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
                GW::extractFarFieldPsi4(m_uiMesh,(const DendroScalar **)consVar,m_uiTinfo._m_uiStep,m_uiTinfo._m_uiT);
                #endif
            #endif

        #endif

        #ifdef BSSN_ENABLE_VTU_OUTPUT

            if((m_uiTinfo._m_uiStep % bssn::BSSN_IO_OUTPUT_FREQ)==0)
            {
                std::vector<std::string> pDataNames;
                const unsigned int numConstVars = bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT;
                const unsigned int numEvolVars  = bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT;

                double *pData[(numConstVars+numEvolVars)];

                for(unsigned int i=0; i<numEvolVars; i++)
                {
                    pDataNames.push_back(std::string(bssn::BSSN_VAR_NAMES[BSSN_VTU_OUTPUT_EVOL_INDICES[i]]));
                    pData[i]=evolVar[BSSN_VTU_OUTPUT_EVOL_INDICES[i]];
                }


                for(unsigned int i=0; i<numConstVars; i++)
                {
                    pDataNames.push_back(std::string(bssn::BSSN_CONSTRAINT_VAR_NAMES[BSSN_VTU_OUTPUT_CONST_INDICES[i]]));
                    pData[numEvolVars+i]=consVar[BSSN_VTU_OUTPUT_CONST_INDICES[i]];
                }

                std::vector<char*> pDataNames_char;
                pDataNames_char.reserve(pDataNames.size());

                for(unsigned int  i = 0; i < pDataNames.size(); i++)
                    pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

                const char * fDataNames[]= {"Time","Cycle"};
                const double fData[]= {m_uiTinfo._m_uiT,(double)m_uiTinfo._m_uiStep};

                char fPrefix[256];
                sprintf(fPrefix,"%s_%d",bssn::BSSN_VTU_FILE_PREFIX.c_str(),m_uiTinfo._m_uiStep);

                if(bssn::BSSN_VTU_Z_SLICE_ONLY)
                {
                    unsigned int s_val[3]= {1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1)};
                    unsigned int s_norm[3] ={0,0,1};
                    io::vtk::mesh2vtu_slice(m_uiMesh,s_val, s_norm, fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);
                }
                else
                    io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);
                

            }

            
        #endif


        #ifdef BSSN_EXTRACT_BH_LOCATIONS
            bssn::writeBHCoordinates((const ot::Mesh *)m_uiMesh,(const Point *) m_uiBHLoc,2,m_uiTinfo._m_uiStep,m_uiTinfo._m_uiT);
        #endif

        return 0;
    }

    int BSSNCtxGPU::write_checkpt()
    {   
        if(!m_uiMesh->isActive())
            return 0;

        
        unsigned int cpIndex;
        (m_uiTinfo._m_uiStep %(2*bssn::BSSN_CHECKPT_FREQ)==0) ? cpIndex=0 : cpIndex=1; // to support alternate file writing.

        const bool is_merged = ((bssn::BSSN_BH_LOC[0]-bssn::BSSN_BH_LOC[1]).abs() < 0.1);
        if(is_merged && !bssn::BSSN_MERGED_CHKPT_WRITTEN)
        {
            cpIndex=3;
            bssn::BSSN_MERGED_CHKPT_WRITTEN=true;
        }
        
        unsigned int rank=m_uiMesh->getMPIRank();
        unsigned int npes=m_uiMesh->getMPICommSize();

        DendroScalar* eVar[BSSN_NUM_VARS];
        DVec & m_evar =m_var[VL::CPU_EV];
        m_evar.to_2d(eVar);
        

        char fName[256];
        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
        sprintf(fName,"%s_octree_%d_%d.oct",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
        io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

        unsigned int numVars=bssn::BSSN_NUM_VARS;
        const char ** varNames=bssn::BSSN_VAR_NAMES;

        /*for(unsigned int i=0;i<numVars;i++)
        {
            sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,rank);
            io::checkpoint::writeVecToFile(fName,m_uiMesh,m_uiPrevVar[i]);
        }*/

        sprintf(fName,"%s_%d_%d.var",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
        io::checkpoint::writeVecToFile(fName,m_uiMesh,(const double **)eVar,bssn::BSSN_NUM_VARS);

        if(!rank)
        {
            sprintf(fName,"%s_step_%d.cp",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex);
            std::cout<<"[BSSNCtx] \t writing checkpoint file : "<<fName<<std::endl;
            std::ofstream outfile(fName);
            if(!outfile) {
                std::cout<<fName<<" file open failed "<<std::endl;
                return 0;
            }

            json checkPoint;
            checkPoint["DENDRO_TS_TIME_BEGIN"]    = m_uiTinfo._m_uiTb;
            checkPoint["DENDRO_TS_TIME_END"]      = m_uiTinfo._m_uiTe;
            checkPoint["DENDRO_TS_ELEMENT_ORDER"] = m_uiElementOrder;

            checkPoint["DENDRO_TS_TIME_CURRENT"]    = m_uiTinfo._m_uiT;
            checkPoint["DENDRO_TS_STEP_CURRENT"]    = m_uiTinfo._m_uiStep;
            checkPoint["DENDRO_TS_TIME_STEP_SIZE"]  = m_uiTinfo._m_uiTh;
            checkPoint["DENDRO_TS_LAST_IO_TIME"]    = m_uiTinfo._m_uiT;

            checkPoint["DENDRO_TS_WAVELET_TOLERANCE"]=bssn::BSSN_WAVELET_TOL;
            checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"]=bssn::BSSN_LOAD_IMB_TOL;
            checkPoint["DENDRO_TS_NUM_VARS"]=numVars; // number of variables to restore.
            checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"]=m_uiMesh->getMPICommSize(); // (note that rank 0 is always active).
            
            checkPoint["DENDRO_BH1_X"]=m_uiBHLoc[0].x();
            checkPoint["DENDRO_BH1_Y"]=m_uiBHLoc[0].y();
            checkPoint["DENDRO_BH1_Z"]=m_uiBHLoc[0].z();
            
            
            checkPoint["DENDRO_BH2_X"]=m_uiBHLoc[1].x();
            checkPoint["DENDRO_BH2_Y"]=m_uiBHLoc[1].y();
            checkPoint["DENDRO_BH2_Z"]=m_uiBHLoc[1].z();
            

            outfile<<std::setw(4)<<checkPoint<<std::endl;
            outfile.close();

        }

        return 0;
    }

    int BSSNCtxGPU::restore_checkpt()
    {
        unsigned int numVars=0;
        std::vector<ot::TreeNode> octree;
        json checkPoint;

        int rank;
        int npes;
        MPI_Comm comm=m_uiMesh->getMPIGlobalCommunicator();
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int activeCommSz;

        char fName[256];
        unsigned int restoreStatus=0;
        unsigned int restoreStatusGlobal=0; // 0 indicates successfully restorable.

        ot::Mesh* newMesh;
        unsigned int restoreStep[2];
        restoreStep[0]=0;
        restoreStep[1]=0;

        unsigned int restoreFileIndex=0;

        for(unsigned int cpIndex=0; cpIndex<2; cpIndex++) {

            restoreStatus=0;

            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp", bssn::BSSN_CHKPT_FILE_PREFIX.c_str() , cpIndex);
                std::ifstream infile(fName);
                if(!infile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    restoreStatus=1;
                }


                if(restoreStatus==0)
                {
                    infile>>checkPoint;
                    m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                    m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                    m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                    m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                    m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                    m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                    
                    bssn::BSSN_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    bssn::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
                    m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                    m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);
                    restoreStep[cpIndex]=m_uiTinfo._m_uiStep;

                }
            }

        }

        if(!rank)
        {
            if(restoreStep[0]<restoreStep[1])
                restoreFileIndex=1;
            else
                restoreFileIndex=0;
        }

        par::Mpi_Bcast(&restoreFileIndex,1,0,comm);

        restoreStatus=0;
        octree.clear();
        if(!rank) std::cout<<"[BSSNCtx] :  Trying to restore from checkpoint index : "<<restoreFileIndex<<std::endl;
    
        if(!rank)
        {
            sprintf(fName,"%s_step_%d.cp", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(), restoreFileIndex);
            std::ifstream infile(fName);
            if(!infile) {
                std::cout<<fName<<" file open failed "<<std::endl;
                restoreStatus=1;
            }


            if(restoreStatus==0)
            {
                infile>>checkPoint;
                m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                
                bssn::BSSN_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                bssn::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                
                numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                
                m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);
                restoreStep[restoreFileIndex]=m_uiTinfo._m_uiStep;
            }


        }

        par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
        if(restoreStatusGlobal==1) 
        {
            if(!rank)
                std::cout<<"[BSSNCtx] : Restore step failed, restore file corrupted. "<<std::endl;
            MPI_Abort(comm,0);
        }


        MPI_Bcast(&m_uiTinfo,sizeof(ts::TSInfo),MPI_BYTE,0,comm);
        par::Mpi_Bcast(&bssn::BSSN_WAVELET_TOL,1,0,comm);
        par::Mpi_Bcast(&bssn::BSSN_LOAD_IMB_TOL,1,0,comm);

        par::Mpi_Bcast(&numVars,1,0,comm);
        par::Mpi_Bcast(&m_uiElementOrder,1,0,comm);
        par::Mpi_Bcast(&activeCommSz,1,0,comm);
        
        par::Mpi_Bcast(m_uiBHLoc,2,0,comm);
        bssn::BSSN_BH_LOC[0]=m_uiBHLoc[0];
        bssn::BSSN_BH_LOC[1]=m_uiBHLoc[1];

        if(activeCommSz>npes)
        {
            if(!rank)
                std::cout<<" [BSSNCtx] : checkpoint file written from  a larger communicator than the current global comm. (i.e. communicator shrinking not allowed in the restore step. )"<<std::endl;
            
            MPI_Abort(comm,0);
        }



        bool isActive=(rank<activeCommSz);

        MPI_Comm newComm;
        par::splitComm2way(isActive,&newComm,comm);

        if(isActive) {

            int activeRank;
            int activeNpes;

            MPI_Comm_rank(newComm, &activeRank);
            MPI_Comm_size(newComm, &activeNpes);
            assert(activeNpes == activeCommSz);

            sprintf(fName, "%s_octree_%d_%d.oct", bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
            restoreStatus=io::checkpoint::readOctFromFile(fName, octree);
            assert(par::test::isUniqueAndSorted(octree, newComm));

        }

        par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
        if(restoreStatusGlobal==1) {

            if(!rank) std::cout<<"[BSSNCtx]: octree (*.oct) restore file is corrupted "<<std::endl;
            MPI_Abort(comm,0);
        }

        newMesh=new ot::Mesh(octree,1,m_uiElementOrder,activeCommSz,comm);
        newMesh->setDomainBounds(Point(bssn::BSSN_GRID_MIN_X,bssn::BSSN_GRID_MIN_Y,bssn::BSSN_GRID_MIN_Z), Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y,bssn::BSSN_GRID_MAX_Z));
        // no need to transfer data only to resize the contex variables. 
        //this->grid_transfer(newMesh);
        #ifdef __CUDACC__
            device::MeshGPU*& dptr_mesh  = this->get_meshgpu_device_ptr();
            device::MeshGPU* mesh_gpu    = this->get_meshgpu_host_handle();

            mesh_gpu->dealloc_mesh_on_device(dptr_mesh);
            dptr_mesh = mesh_gpu->alloc_mesh_on_device(newMesh);
        #endif

        for(unsigned int i=0; i < VL::END; i++)
            m_var[i].destroy_vector();

        ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS, BSSN_ASYNC_COMM_K);
        device::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx_device, BSSN_NUM_VARS , BSSN_ASYNC_COMM_K);

        ot::alloc_mpi_ctx<DendroScalar>  (newMesh, m_mpi_ctx, BSSN_NUM_VARS, BSSN_ASYNC_COMM_K);
        device::alloc_mpi_ctx<DendroScalar>  (newMesh, m_mpi_ctx_device, BSSN_NUM_VARS , BSSN_ASYNC_COMM_K);
        
        // variable allocation for evolution variables
        m_var[VL::CPU_EV].create_vector(newMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::HOST,BSSN_NUM_VARS,true);
        m_var[VL::CPU_CV].create_vector(newMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::HOST,BSSN_CONSTRAINT_NUM_VARS,true);
        m_var[VL::CPU_EV_UZ_IN].create_vector(newMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,BSSN_NUM_VARS,true);
        m_var[VL::CPU_CV_UZ_IN].create_vector(newMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,BSSN_CONSTRAINT_NUM_VARS,true);
        
        m_var[VL::GPU_EV].create_vector(newMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_IN].create_vector(newMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_OUT].create_vector(newMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);

        // only reads the evolution variables. 
        if(isActive) {

            int activeRank;
            int activeNpes;

            DendroScalar* inVec[BSSN_NUM_VARS];
            DVec & m_evar = m_var[VL::CPU_EV];
            m_evar.to_2d(inVec);

            MPI_Comm_rank(newComm, &activeRank);
            MPI_Comm_size(newComm, &activeNpes);
            assert(activeNpes == activeCommSz);

            sprintf(fName,"%s_%d_%d.var",bssn::BSSN_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
            restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,inVec,bssn::BSSN_NUM_VARS);

            
        }

        MPI_Comm_free(&newComm);
        par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
        if(restoreStatusGlobal==1) {

            if(!rank) std::cout<<"[BSSNCtx]: varible (*.var) restore file currupted "<<std::endl;
            MPI_Abort(comm,0);
        }

        std::swap(m_uiMesh,newMesh);
        delete newMesh;

        // realloc bssn deriv space
        deallocate_bssn_deriv_workspace();
        allocate_bssn_deriv_workspace(m_uiMesh,1);
            
        unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
        unsigned int totalElems=0;
        par::Mpi_Allreduce(&localSz, &totalElems ,1,MPI_SUM,comm);
        
        if(!rank) std::cout<<" checkpoint at step : "<<m_uiTinfo._m_uiStep<<"active Comm. sz: "<<activeCommSz<<" restore successful: "<<" restored mesh size: "<<totalElems<<std::endl;

        m_uiIsETSSynced = false;
        return 0;

    }

    int BSSNCtxGPU::rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time)
    {
        if(!m_uiMesh->isActive())
            return 0;
        
        const unsigned int pencils = 13;
        const unsigned int pen_sz  = 13;
        
        

        DEVICE_REAL* const out_ptr = out->get_vec_ptr();
        const DEVICE_REAL * const in_ptr  = in->get_vec_ptr();
        

        DVec& m_dptr_uz_i = m_var[VL::GPU_EV_UZ_IN];
        DVec& m_dptr_uz_o = m_var[VL::GPU_EV_UZ_OUT];

        const unsigned int unzip_sz = m_uiMesh->getDegOfFreedomUnZip();

        this->unzip(*in,m_dptr_uz_i,BSSN_ASYNC_COMM_K);
        
        unsigned int nblocks = m_uiMesh->getLocalBlockList().size();
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::RHS].start();
        #endif
        //device::bssnrhs2<3, pen_sz, pen_sz, DEVICE_RHS_BATCHED_GRAIN_SZ,DEVICE_RHS_NSTREAMS> (m_dptr_uz_o.get_vec_ptr(),m_dptr_uz_i.get_vec_ptr(),m_mesh_cpu.m_blk_list, nblocks, m_dptr_deriv_evars, unzip_sz, KO_DISS_SIGMA);
        device::bssnrhs3<3, pen_sz, pen_sz, DEVICE_RHS_BATCHED_GRAIN_SZ,DEVICE_RHS_NSTREAMS> (m_dptr_uz_o.get_vec_ptr(),m_dptr_uz_i.get_vec_ptr(),m_mesh_cpu.m_blk_list, nblocks, unzip_sz, KO_DISS_SIGMA);
        #ifdef __PROFILE_CTX__
            GPUDevice::device_synchronize();
            m_uiCtxpt[ts::CTXPROFILE::RHS].stop();
        #endif

        GPUDevice::check_last_error();

        this->zip(m_dptr_uz_o,*out);
        
        return 0;
        
    }

    int BSSNCtxGPU::post_stage(DVec& sIn)
    {
        // const unsigned int lb = m_uiMesh->getNodeLocalBegin();
        // const unsigned int le = m_uiMesh->getNodeLocalEnd();
        // const unsigned int szpdof = sIn.get_size()/sIn.get_dof();
        // device::cuda_bssn_enforce_evar_cons<<< (le-lb)/1024 + 1, 1024>>>(sIn.get_vec_ptr(),lb,le,bssn::CHI_FLOOR,szpdof);
        // GPUDevice::device_synchronize();
        // GPUDevice::check_last_error();
        
        return 0;

    }

    int BSSNCtxGPU::post_timestep(DVec& sIn)
    {
        //need to enforce constraints in cuda. 
        const unsigned int lb = m_uiMesh->getNodeLocalBegin();
        const unsigned int le = m_uiMesh->getNodeLocalEnd();
        const unsigned int szpdof = sIn.get_size()/sIn.get_dof();
        device::cuda_bssn_enforce_evar_cons<<< (le-lb)/1024 + 1, 1024>>>(sIn.get_vec_ptr(),lb,le,bssn::CHI_FLOOR,szpdof);
        GPUDevice::device_synchronize();
        GPUDevice::check_last_error();
        
        return 0;
    }

    bool BSSNCtxGPU::is_remesh()
    {
        
        bool isRefine = false;
        if(bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY)
            return false;

        MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();

        DVec&  m_evar = m_var[VL::CPU_EV];
        DVec&  m_evar_unz = m_var[VL::CPU_EV_UZ_IN];
        
        this->unzip(m_evar,m_evar_unz, bssn::BSSN_ASYNC_COMM_K);

        DendroScalar* unzipVar[BSSN_NUM_VARS];
        m_evar_unz.to_2d(unzipVar);

        unsigned int refineVarIds[bssn::BSSN_NUM_REFINE_VARS];
        for(unsigned int vIndex=0; vIndex<bssn::BSSN_NUM_REFINE_VARS; vIndex++)
            refineVarIds[vIndex]=bssn::BSSN_REFINE_VARIABLE_INDICES[vIndex];

        double wTol=bssn::BSSN_WAVELET_TOL;
        std::function<double(double,double,double,double*hx)> waveletTolFunc =[](double x,double y, double z,double*hx) {
            return bssn::computeWTolDCoords(x,y,z,hx);
        };
        

        if(bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::WAMR) {
            isRefine = bssn::isReMeshWAMR(m_uiMesh,(const double **)unzipVar,refineVarIds,bssn::BSSN_NUM_REFINE_VARS,waveletTolFunc,bssn::BSSN_DENDRO_AMR_FAC);

        }else if(bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::EH)
        {
            isRefine = bssn::isRemeshEH(m_uiMesh,(const double **)unzipVar,bssn::VAR::U_ALPHA,bssn::BSSN_EH_REFINE_VAL,bssn::BSSN_EH_COARSEN_VAL,true);

        }else if(bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::EH_WAMR)
        {
            const bool isR1 = bssn::isReMeshWAMR(m_uiMesh,(const double **)unzipVar,refineVarIds,bssn::BSSN_NUM_REFINE_VARS,waveletTolFunc,bssn::BSSN_DENDRO_AMR_FAC);
            const bool isR2 = bssn::isRemeshEH(m_uiMesh,(const double **)unzipVar,bssn::VAR::U_ALPHA,bssn::BSSN_EH_REFINE_VAL,bssn::BSSN_EH_COARSEN_VAL,false);

            isRefine = (isR1 || isR2);
        }else if( bssn::BSSN_REFINEMENT_MODE == bssn::RefinementMode::BH_LOC)
        {
            isRefine = bssn::isRemeshBH(m_uiMesh,m_uiBHLoc);
        }

        return isRefine;

    }
    
    DVec& BSSNCtxGPU::get_evolution_vars()
    {
        return m_var[VL::GPU_EV];
    }

    DVec& BSSNCtxGPU::get_evolution_vars_cpu()
    {
        return m_var[VL::CPU_EV];
    }
    
    DVec& BSSNCtxGPU::get_constraint_vars()
    {
        return m_var[VL::CPU_CV];
    }

    // DVec& BSSNCtxGPU::get_primitive_vars()
    // {
    //     return m_pvar;
    // }

    int BSSNCtxGPU::terminal_output()
    {
        if(m_uiMesh->isActive())
        {
            DendroScalar min=0, max=0;
            DVec& m_evar = m_var[VL::CPU_EV];
            min=vecMin(m_uiMesh,m_evar.get_vec_ptr(),ot::VEC_TYPE::CG_NODAL,true);
            max=vecMax(m_uiMesh,m_evar.get_vec_ptr(),ot::VEC_TYPE::CG_NODAL,true);
            
            if(!(m_uiMesh->getMPIRank()))
            {
                std::cout<<"[BSSNCtx]:  "<<bssn::BSSN_VAR_NAMES[bssn::VAR::U_ALPHA]<<" (min,max) : \t ( "<<min<<", "<<max<<" ) "<<std::endl;
                if(std::isnan(min) || std::isnan(max))
                {
                    std::cout<<"[Error]: NAN detected "<<std::endl;
                    MPI_Abort(m_uiMesh->getMPICommunicator(),0);
                }

            }
                
        }

        return 0; 
    }

    int BSSNCtxGPU::grid_transfer(const ot::Mesh* m_new)
    {
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].start();
        #endif
        DVec& m_evar  = m_var[VL::CPU_EV];
        DVec::grid_transfer(m_uiMesh, m_new, m_evar);
        ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx, BSSN_NUM_VARS, BSSN_ASYNC_COMM_K);
        ot::alloc_mpi_ctx<DendroScalar>  (m_new, m_mpi_ctx, BSSN_NUM_VARS, BSSN_ASYNC_COMM_K);

        device::dealloc_mpi_ctx<DendroScalar>(m_uiMesh, m_mpi_ctx_device, BSSN_NUM_VARS , BSSN_ASYNC_COMM_K);
        device::alloc_mpi_ctx<DendroScalar>  (m_new, m_mpi_ctx_device, BSSN_NUM_VARS , BSSN_ASYNC_COMM_K);
        //printf("igt ended\n");

        m_var[VL::CPU_CV].destroy_vector();
        m_var[VL::CPU_CV_UZ_IN].destroy_vector();
        m_var[VL::CPU_EV_UZ_IN].destroy_vector();

        m_var[VL::GPU_EV].destroy_vector();
        m_var[VL::GPU_EV_UZ_IN].destroy_vector();
        m_var[VL::GPU_EV_UZ_OUT].destroy_vector();

        m_var[VL::CPU_CV].create_vector(m_new,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::HOST,BSSN_CONSTRAINT_NUM_VARS,true);
        m_var[VL::CPU_CV_UZ_IN].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,BSSN_CONSTRAINT_NUM_VARS,true);

        m_var[VL::CPU_EV_UZ_IN].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,BSSN_NUM_VARS,true);
        
        m_var[VL::GPU_EV].create_vector(m_new,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_IN].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_OUT].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,BSSN_NUM_VARS,true);
        
        this->host_to_device_sync();
        m_uiIsETSSynced=false;
        
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].stop();
        #endif
        return 0 ;
    }


    void BSSNCtxGPU::evolve_bh_loc(DVec sIn, double dt)
    {
        #ifdef BSSN_EXTRACT_BH_LOCATIONS

            Point bhLoc[2];
            DVec& m_evar     = m_var[VL::CPU_EV];
            DVec& m_evar_unz = m_var[VL::CPU_EV_UZ_IN];
            DVec& m_cvar     = m_var[VL::CPU_CV];
            DVec & m_cvar_unz= m_var[VL::CPU_CV_UZ_IN];
            
            DendroScalar * evar[bssn::BSSN_NUM_VARS];
            sIn.to_2d(evar);
            bssn::computeBHLocations((const ot::Mesh *)m_uiMesh,m_uiBHLoc,bhLoc,evar,dt);
            // if(!m_uiMesh->getMPIRankGlobal())
            // {
            //     std::cout<<"bh0 "<<bhLoc[0]<<"dt : "<<dt<<std::endl;
            //     std::cout<<"bh1 "<<bhLoc[1]<<"dt : "<<dt<<std::endl;

            // }
            m_uiBHLoc[0] = bhLoc[0];
            m_uiBHLoc[1] = bhLoc[1];

            bssn::BSSN_BH_LOC[0]=m_uiBHLoc[0];
            bssn::BSSN_BH_LOC[1]=m_uiBHLoc[1];

            // old bh location extractor. 
            #if 0
                DendroScalar *evar[bssn::BSSN_NUM_VARS];
                Point bhLoc[2];
                sIn.to_2d(evar);
                bssn::extractBHCoords((const ot::Mesh *)m_uiMesh,(const DendroScalar*)evar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,(const Point *) m_uiBHLoc,2,(Point*)bhLoc);
                
                m_uiBHLoc[0] = bhLoc[0];
                m_uiBHLoc[1] = bhLoc[1];
            #endif
        #endif

        return;

    }
}
