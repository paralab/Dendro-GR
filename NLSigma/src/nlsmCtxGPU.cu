/**
 * @file nlsmCtxGPU.cu
 * @brief NLSM GPU ctx file
 * @version 0.1
 * @date 2022-02-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "nlsmCtxGPU.cuh"
CONST_MEM DEVICE_REAL device::refel_1d[2*REFEL_CONST_MEM_MAX];

namespace nlsm
{
    template<int pw, int pencils, int pencil_sz, int BATCHED_BLOCKS_SZ>
    GLOBAL_FUNC void launch_dir_x_deriv_kernel(const device::MeshGPU* const dptr_mesh, const DEVICE_REAL * const  u, EVAR_DERIVS * deriv_evars , BlockGPU3D* blk, DEVICE_UINT blk_beging) 
    {
        const DEVICE_UINT blk_k             = GPUDevice::block_id_z();
        const DEVICE_UINT BLK_ID            = blk_beging + blk_k;
        const DEVICE_UINT BLK_SZ            = (blk[BLK_ID].m_aligned_sz[0]) * (blk[BLK_ID].m_aligned_sz[1]) * (blk[BLK_ID].m_aligned_sz[2]); 
        const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
        const DEVICE_UINT sz_per_dof        = dptr_mesh->m_oct_unzip_sz;

        const DEVICE_REAL * const chi       = &u[VAR::U_CHI * sz_per_dof + offset];
        const DEVICE_REAL * const phi       = &u[VAR::U_PHI * sz_per_dof + offset];
        
        DEVICE_REAL * const grad_0_chi      = deriv_evars->grad_0_chi + blk_k * BLK_SZ;
        DEVICE_REAL * const grad_0_phi      = deriv_evars->grad_0_phi + blk_k * BLK_SZ;

        DEVICE_REAL * const ko_grad_0_chi   = deriv_evars->ko_grad_0_chi + blk_k * BLK_SZ;
        DEVICE_REAL * const ko_grad_0_phi   = deriv_evars->ko_grad_0_phi + blk_k * BLK_SZ;

        DEVICE_REAL * const grad2_0_0_chi   = deriv_evars->grad2_0_0_chi + blk_k * BLK_SZ;
                
        device::__deriv644_x   <pw,pencils,pencil_sz>(grad_0_chi, chi , blk + BLK_ID);
        device::__deriv644_x   <pw,pencils,pencil_sz>(grad_0_phi, phi , blk + BLK_ID);
        
        device::__ko_deriv42_x <pw,pencils,pencil_sz>(ko_grad_0_chi, chi, blk + BLK_ID);
        device::__ko_deriv42_x <pw,pencils,pencil_sz>(ko_grad_0_phi, phi, blk + BLK_ID);

        device::__deriv644_xx<pw,pencils,pencil_sz>(grad2_0_0_chi, chi, blk + BLK_ID);
        return;
    }

    template<int pw, int pencils, int pencil_sz, int BATCHED_BLOCKS_SZ>
    GLOBAL_FUNC void launch_dir_y_deriv_kernel(const device::MeshGPU* const dptr_mesh,const DEVICE_REAL * const  u, EVAR_DERIVS * deriv_evars , BlockGPU3D* blk, DEVICE_UINT blk_beging) 
    {
        const DEVICE_UINT blk_k             = GPUDevice::block_id_z();
        const DEVICE_UINT BLK_ID            = blk_beging + blk_k;
        const DEVICE_UINT BLK_SZ            = (blk[BLK_ID].m_aligned_sz[0]) * (blk[BLK_ID].m_aligned_sz[1]) * (blk[BLK_ID].m_aligned_sz[2]); 
        const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
        const DEVICE_UINT sz_per_dof        = dptr_mesh->m_oct_unzip_sz;

        const DEVICE_REAL * const chi       = &u[VAR::U_CHI * sz_per_dof + offset];
        const DEVICE_REAL * const phi       = &u[VAR::U_PHI * sz_per_dof + offset];
        
        DEVICE_REAL * const grad_1_chi      = deriv_evars->grad_1_chi + blk_k * BLK_SZ;
        DEVICE_REAL * const grad_1_phi      = deriv_evars->grad_1_phi + blk_k * BLK_SZ;

        DEVICE_REAL * const ko_grad_1_chi   = deriv_evars->ko_grad_1_chi + blk_k * BLK_SZ;
        DEVICE_REAL * const ko_grad_1_phi   = deriv_evars->ko_grad_1_phi + blk_k * BLK_SZ;

        DEVICE_REAL * const grad2_1_1_chi   = deriv_evars->grad2_1_1_chi + blk_k * BLK_SZ;
                
        device::__deriv644_y   <pw,pencils,pencil_sz>(grad_1_chi, chi , blk + BLK_ID);
        device::__deriv644_y   <pw,pencils,pencil_sz>(grad_1_phi, phi , blk + BLK_ID);
        
        device::__ko_deriv42_y <pw,pencils,pencil_sz>(ko_grad_1_chi, chi, blk + BLK_ID);
        device::__ko_deriv42_y <pw,pencils,pencil_sz>(ko_grad_1_phi, phi, blk + BLK_ID);

        device::__deriv644_yy<pw,pencils,pencil_sz>(grad2_1_1_chi, chi, blk + BLK_ID);
        return;
    }

    template<int pw, int pencils, int pencil_sz, int BATCHED_BLOCKS_SZ>
    GLOBAL_FUNC void launch_dir_z_deriv_kernel(const device::MeshGPU* const dptr_mesh, const DEVICE_REAL * const  u, EVAR_DERIVS * deriv_evars , BlockGPU3D* blk, DEVICE_UINT blk_beging) 
    {
        const DEVICE_UINT blk_k             = GPUDevice::block_id_z();
        const DEVICE_UINT BLK_ID            = blk_beging + blk_k;
        const DEVICE_UINT BLK_SZ            = (blk[BLK_ID].m_aligned_sz[0]) * (blk[BLK_ID].m_aligned_sz[1]) * (blk[BLK_ID].m_aligned_sz[2]); 
        const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
        const DEVICE_UINT sz_per_dof        = dptr_mesh->m_oct_unzip_sz;

        const DEVICE_REAL * const chi       = &u[VAR::U_CHI * sz_per_dof + offset];
        const DEVICE_REAL * const phi       = &u[VAR::U_PHI * sz_per_dof + offset];

        DEVICE_REAL * const grad_2_chi      = deriv_evars->grad_2_chi + blk_k * BLK_SZ;
        DEVICE_REAL * const grad_2_phi      = deriv_evars->grad_2_phi + blk_k * BLK_SZ;

        DEVICE_REAL * const ko_grad_2_chi   = deriv_evars->ko_grad_2_chi + blk_k * BLK_SZ;
        DEVICE_REAL * const ko_grad_2_phi   = deriv_evars->ko_grad_2_phi + blk_k * BLK_SZ;

        DEVICE_REAL * const grad2_2_2_chi   = deriv_evars->grad2_2_2_chi + blk_k * BLK_SZ;
                
        device::__deriv644_z   <pw,pencils,pencil_sz>(grad_2_chi, chi , blk + BLK_ID);
        device::__deriv644_z   <pw,pencils,pencil_sz>(grad_2_phi, phi , blk + BLK_ID);
        
        device::__ko_deriv42_z <pw,pencils,pencil_sz>(ko_grad_2_chi, chi, blk + BLK_ID);
        device::__ko_deriv42_z <pw,pencils,pencil_sz>(ko_grad_2_phi, phi, blk + BLK_ID);

        device::__deriv644_zz<pw,pencils,pencil_sz>(grad2_2_2_chi, chi, blk + BLK_ID);

        return;
    }

    template<int pw, int pencils, int pencil_sz, int BATCHED_BLOCKS_SZ>
    GLOBAL_FUNC void launch_rhs(const device::MeshGPU* const dptr_mesh, DEVICE_REAL * const Fu, const DEVICE_REAL * const  u, EVAR_DERIVS* deriv_evars, BlockGPU3D* blk, DEVICE_UINT blk_beging)
    {

        // set up kernel parameters for RHS GPU RHS kernels. 
        const DEVICE_UINT blk_k             = GPUDevice::block_id_z();
        const DEVICE_UINT BLK_ID            = blk_beging + blk_k;

        const DEVICE_UINT BLK_SZ            = (blk[BLK_ID].m_aligned_sz[0]) * (blk[BLK_ID].m_aligned_sz[1]) * (blk[BLK_ID].m_aligned_sz[2]); 
        const DEVICE_UINT offset            = blk[BLK_ID].m_offset; 
        const DEVICE_UINT bflag             = blk[BLK_ID].m_bflag; 
        const DEVICE_UINT sz_per_dof        = dptr_mesh->m_oct_unzip_sz;

        const DEVICE_INT nx = blk[BLK_ID].m_sz[0];
        const DEVICE_INT ny = blk[BLK_ID].m_sz[1];
        const DEVICE_INT nz = blk[BLK_ID].m_sz[2];

        const DEVICE_REAL hx = blk[BLK_ID].m_dx[0];
        const DEVICE_REAL hy = blk[BLK_ID].m_dx[1];
        const DEVICE_REAL hz = blk[BLK_ID].m_dx[2];

        const DEVICE_INT i  = GPUDevice::thread_id_x();
        const DEVICE_INT j  = GPUDevice::block_id_x() * GPUDevice::block_dim_y() + GPUDevice::thread_id_y();
        const DEVICE_INT k  = GPUDevice::block_id_y();

        const DEVICE_INT si = GPUDevice::thread_id_x();  
        const DEVICE_INT sj = GPUDevice::thread_id_y();

        const DEVICE_INT sb = si;
        const DEVICE_INT pp = k * nx * ny + j * nx + sb;

        const DEVICE_REAL x = blk[BLK_ID].m_ptMin[0] + i*hx;
        const DEVICE_REAL y = blk[BLK_ID].m_ptMin[1] + j*hy;
        const DEVICE_REAL z = blk[BLK_ID].m_ptMin[2] + k*hz;

        DEVICE_REAL * const chi_rhs    = &Fu[VAR::U_CHI * sz_per_dof + offset];
        DEVICE_REAL * const phi_rhs    = &Fu[VAR::U_PHI * sz_per_dof + offset];

        const DEVICE_REAL* const chi   = &u[VAR::U_CHI * sz_per_dof + offset];
        const DEVICE_REAL* const phi   = &u[VAR::U_PHI * sz_per_dof + offset];

        const DEVICE_REAL* const grad_0_chi = (deriv_evars->grad_0_chi + blk_k * BLK_SZ);
        const DEVICE_REAL* const grad_1_chi = (deriv_evars->grad_1_chi + blk_k * BLK_SZ);
        const DEVICE_REAL* const grad_2_chi = (deriv_evars->grad_2_chi + blk_k * BLK_SZ);

        const DEVICE_REAL* const grad_0_phi = (deriv_evars->grad_0_phi + blk_k * BLK_SZ);
        const DEVICE_REAL* const grad_1_phi = (deriv_evars->grad_1_phi + blk_k * BLK_SZ);
        const DEVICE_REAL* const grad_2_phi = (deriv_evars->grad_2_phi + blk_k * BLK_SZ);


        const DEVICE_REAL* const ko_grad_0_chi = (deriv_evars->ko_grad_0_chi + blk_k * BLK_SZ);
        const DEVICE_REAL* const ko_grad_1_chi = (deriv_evars->ko_grad_1_chi + blk_k * BLK_SZ);
        const DEVICE_REAL* const ko_grad_2_chi = (deriv_evars->ko_grad_2_chi + blk_k * BLK_SZ);

        const DEVICE_REAL* const ko_grad_0_phi = (deriv_evars->ko_grad_0_phi + blk_k * BLK_SZ);
        const DEVICE_REAL* const ko_grad_1_phi = (deriv_evars->ko_grad_1_phi + blk_k * BLK_SZ);
        const DEVICE_REAL* const ko_grad_2_phi = (deriv_evars->ko_grad_2_phi + blk_k * BLK_SZ);

        const DEVICE_REAL* const grad2_0_0_chi = (deriv_evars->grad2_0_0_chi + blk_k * BLK_SZ);
        const DEVICE_REAL* const grad2_1_1_chi = (deriv_evars->grad2_1_1_chi + blk_k * BLK_SZ);
        const DEVICE_REAL* const grad2_2_2_chi = (deriv_evars->grad2_2_2_chi + blk_k * BLK_SZ);
        
        const DEVICE_REAL r = sqrt(x * x  + y * y + z * z);

        const DEVICE_REAL NLSM_WAVE_SPEED_X = 1.0;
        const DEVICE_REAL NLSM_WAVE_SPEED_Y = 1.0;
        const DEVICE_REAL NLSM_WAVE_SPEED_Z = 1.0;
        const DEVICE_REAL sigma             = 0.4;
        
        #ifdef NLSM_NONLINEAR
          //#pragma message("nl")
          if (r > 1.0e-6) {
            phi_rhs[pp] =  NLSM_WAVE_SPEED_X* grad2_0_0_chi[pp] + NLSM_WAVE_SPEED_Y*grad2_1_1_chi[pp] + NLSM_WAVE_SPEED_Z*grad2_2_2_chi[pp] - sin(2*chi[pp])/pow(r, 2);
            chi_rhs[pp] = phi[pp];
          } else {
            chi_rhs[pp] = 0.0;
            phi_rhs[pp] = 0.0;
          }
        #else
           //#pragma message("linear")
           phi_rhs[pp] =  NLSM_WAVE_SPEED_X*grad2_0_0_chi[pp] + NLSM_WAVE_SPEED_Y*grad2_1_1_chi[pp] + NLSM_WAVE_SPEED_Z*grad2_2_2_chi[pp];
           chi_rhs[pp] =  phi[pp];
        #endif

        if(bflag!=0)
        {
            device::radiative_bc<3>(chi_rhs, chi, grad_0_chi , grad_1_chi, grad_2_chi, 1.0, 0.0, blk, BLK_ID);
            device::radiative_bc<3>(phi_rhs, phi, grad_0_phi , grad_1_phi, grad_2_phi, 1.0, 0.0, blk, BLK_ID);
        }

        phi_rhs[pp] += sigma * (ko_grad_0_phi[pp] + ko_grad_1_phi[pp] + ko_grad_2_phi[pp]);
        chi_rhs[pp] += sigma * (ko_grad_0_chi[pp] + ko_grad_1_chi[pp] + ko_grad_2_chi[pp]);
        
        return;

    }

    NLSMCtxGPU::NLSMCtxGPU(ot::Mesh* pMesh) : Ctx()
    {
        m_uiMesh      = pMesh;
        m_mesh_cpu  = device::MeshGPU();
        m_dptr_mesh = m_mesh_cpu.alloc_mesh_on_device(m_uiMesh);
        
        // variable allocation for evolution variables
        m_var[VL::CPU_EV].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::HOST,NLSM_NUM_VARS,true);
        //m_var[VL::CPU_EV_DG].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_NODES,ot::DVEC_LOC::HOST,NLSM_NUM_VARS,true);
        m_var[VL::CPU_EV_UZ].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,NLSM_NUM_VARS,true);
        
        //m_var[VL::GPU_EV_DG].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_NODES,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);
        m_var[VL::GPU_EV].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_IN].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_OUT].create_vector(m_uiMesh,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);

        m_uiTinfo._m_uiStep=0;
        m_uiTinfo._m_uiT = 0;
        m_uiTinfo._m_uiTb = NLSM_RK45_TIME_BEGIN;
        m_uiTinfo._m_uiTe = NLSM_RK45_TIME_END;
        m_uiTinfo._m_uiTh = NLSM_RK45_TIME_STEP_SIZE;     

        m_uiElementOrder = NLSM_ELE_ORDER;

        m_uiMinPt = Point(NLSM_GRID_MIN_X,NLSM_GRID_MIN_Y,NLSM_GRID_MIN_Z);
        m_uiMaxPt = Point(NLSM_GRID_MAX_X,NLSM_GRID_MAX_Y,NLSM_GRID_MAX_Z);

        ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh,m_mpi_ctx,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);
        ot::alloc_mpi_ctx<DendroScalar>(m_uiMesh,m_mpi_ctx,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);

        device::dealloc_mpi_ctx<DendroScalar>(m_uiMesh,m_mpi_ctx_device,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);
        device::alloc_mpi_ctx<DendroScalar>(m_uiMesh,m_mpi_ctx_device,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);

        unsigned int BLK_SZ           = DEVICE_RHS_BLK_SZ;
        EVAR_DERIVS* deriv_evars      = GPUDevice::host_malloc<EVAR_DERIVS>(1);
        EVAR_DERIVS* dptr_deriv_evars = GPUDevice::device_malloc<EVAR_DERIVS>(1);
        DEVICE_REAL * deriv_base      = GPUDevice::device_malloc<DEVICE_REAL>(16 * DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ);
        

        deriv_evars -> grad_0_chi    = deriv_base + 0*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad_1_chi    = deriv_base + 1*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad_2_chi    = deriv_base + 2*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad_0_phi    = deriv_base + 3*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad_1_phi    = deriv_base + 4*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad_2_phi    = deriv_base + 5*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> ko_grad_0_chi = deriv_base + 6*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> ko_grad_1_chi = deriv_base + 7*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> ko_grad_2_chi = deriv_base + 8*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> ko_grad_0_phi = deriv_base + 9*   DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> ko_grad_1_phi = deriv_base + 10*  DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> ko_grad_2_phi = deriv_base + 11*  DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad2_0_0_chi = deriv_base + 12*  DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad2_1_1_chi = deriv_base + 13*  DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        deriv_evars -> grad2_2_2_chi = deriv_base + 15*  DEVICE_RHS_BATCHED_GRAIN_SZ * BLK_SZ;
        GPUDevice::host_to_device(deriv_evars, dptr_deriv_evars,1);

        m_deriv_evars       = deriv_evars;
        m_dptr_deriv_evars  = dptr_deriv_evars;
        m_dptr_deriv_base   = deriv_base;

        return;

    }

    NLSMCtxGPU::~NLSMCtxGPU()
    {
        for(unsigned int i=0; i < VL::END; i++)
            m_var[i].destroy_vector();

        GPUDevice::device_free(m_dptr_deriv_base);
        GPUDevice::device_free(m_dptr_deriv_evars);
        GPUDevice::host_free(m_deriv_evars);

        return;
    }

    int NLSMCtxGPU::host_to_device_sync()
    {
        if(!m_uiMesh->isActive())
            return 0;

        DVec& m_evar         = m_var[VL::CPU_EV];
        //DVec& m_evar_dg      = m_var[VL::CPU_EV_DG];
        //DVec& m_dptr_evar    = m_var[VL::GPU_EV_DG];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        GPUDevice::host_to_device(m_evar.get_vec_ptr(), m_dptr_evar.get_vec_ptr(), m_evar.get_size());
        //printf("h2d ended\n");
        return 0;
    }

    int NLSMCtxGPU::device_to_host_sync()
    {
        if(!m_uiMesh->isActive())
            return 0;

        DVec& m_evar         = m_var[VL::CPU_EV];
        //DVec& m_evar_dg      = m_var[VL::CPU_EV_DG];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        GPUDevice::device_to_host(m_evar.get_vec_ptr(), m_dptr_evar.get_vec_ptr(), m_evar.get_size());
        //printf("d2h ended\n");
        return 0;
        
    }

    int NLSMCtxGPU::initialize()
    {
        DVec& m_evar         = m_var[VL::CPU_EV];
        //DVec& m_evar_dg      = m_var[VL::CPU_EV_DG];
        //DVec& m_dptr_evar    = m_var[VL::GPU_EV_DG];
        DVec& m_dptr_evar    = m_var[VL::GPU_EV];

        if(NLSM_RESTORE_SOLVER)
        {
            this->restore_checkpt();
            this->host_to_device_sync();
            return 0; 
        }

        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
        unsigned int eleOrder=m_uiMesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();
        
        DendroScalar* zipIn[NLSM_NUM_VARS];
        m_evar.to_2d(zipIn);

        #pragma omp parallel for
        for(unsigned int elem=m_uiMesh->getElementLocalBegin(); elem<m_uiMesh->getElementLocalEnd(); elem++)
        {
            DendroScalar var[NLSM_NUM_VARS];
            for(unsigned int k=0; k<(eleOrder+1); k++)
                for(unsigned int j=0; j<(eleOrder+1); j++ )
                    for(unsigned int i=0; i<(eleOrder+1); i++)
                    {
                        const unsigned int nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                        {
                            const unsigned int nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            unsigned int ownerID, ii_x, jj_y, kk_z;
                            m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            const double len=(double)(1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                            const double x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                            const double y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                            const double z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));

                            initData(x, y, z , var);

                            for(unsigned int v=0; v<NLSM_NUM_VARS; v++)
                                zipIn[v][nodeLookUp_CG]=var[v];
                            
                        }
                    }
        }

        this->host_to_device_sync();
       
        return 0;


    }

    int NLSMCtxGPU::finalize()
    {
        return 0;    
    }

    int NLSMCtxGPU::write_vtu()
    {
        DVec& m_evar         = m_var[VL::CPU_EV];
        
        DendroScalar * evolUnzipVar = NULL;
        DendroScalar * consUnzipVar = NULL;
        DendroScalar * consVar = NULL;
        DendroScalar * evolVar[NLSM_NUM_VARS];

        m_uiMesh->readFromGhostBegin(m_evar.get_vec_ptr(),m_evar.get_dof());
        m_uiMesh->readFromGhostEnd(m_evar.get_vec_ptr(),m_evar.get_dof());

        m_evar.to_2d(evolVar);

        std::vector<std::string> pDataNames;
        const unsigned int numConstVars = 0;
        const unsigned int numEvolVars  = NLSM_NUM_EVOL_VARS_VTU_OUTPUT;

        double *pData[(numConstVars+numEvolVars)];

        for(unsigned int i=0; i<numEvolVars; i++)
        {
            pDataNames.push_back(std::string(NLSM_VAR_NAMES[NLSM_VTU_OUTPUT_EVOL_INDICES[i]]));
            pData[i]=evolVar[NLSM_VTU_OUTPUT_EVOL_INDICES[i]];
        }

        std::vector<char*> pDataNames_char;
        pDataNames_char.reserve(pDataNames.size());

        for(unsigned int  i = 0; i < pDataNames.size(); i++)
            pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

        const char * fDataNames[]= {"Time","Cycle"};
        const double fData[]= {m_uiTinfo._m_uiT,(double)m_uiTinfo._m_uiStep};

        char fPrefix[256];
        sprintf(fPrefix,"%s_%d",NLSM_VTU_FILE_PREFIX.c_str(),m_uiTinfo._m_uiStep);

        io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);
        
        return 0;

    }

    int NLSMCtxGPU::write_checkpt()
    {   
        DVec& m_evar         = m_var[VL::CPU_EV];
        if(m_uiMesh->isActive())
        {
            unsigned int cpIndex;
            (m_uiTinfo._m_uiStep %(2*NLSM_CHECKPT_FREQ)==0) ? cpIndex=0 : cpIndex=1; // to support alternate file writing.
            unsigned int rank=m_uiMesh->getMPIRank();
            unsigned int npes=m_uiMesh->getMPICommSize();

            char fName[256];
            const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
            sprintf(fName,"%s_octree_%d_%d.oct",NLSM_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

            unsigned int numVars=nlsm::NLSM_NUM_VARS;
            const char ** varNames=nlsm::NLSM_VAR_NAMES;

            /*for(unsigned int i=0;i<numVars;i++)
            {
                sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,rank);
                io::checkpoint::writeVecToFile(fName,m_uiMesh,m_uiPrevVar[i]);
            }*/

            const unsigned int dof = m_evar.get_dof();
            DendroScalar* eVar[dof];
            m_evar.to_2d(eVar);

            sprintf(fName,"%s_%d_%d.var",NLSM_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeVecToFile(fName,m_uiMesh,(const double **)eVar,nlsm::NLSM_NUM_VARS);


            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp",NLSM_CHKPT_FILE_PREFIX.c_str(),cpIndex);
                //std::cout<<"writing : "<<fName<<std::endl;
                std::ofstream outfile(fName);
                if(!outfile) {std::cout<<fName<<" file open failed "<<std::endl; return 0;}

                json checkPoint;

                checkPoint["DENDRO_TS_TIME_BEGIN"]    = m_uiTinfo._m_uiTb;
                checkPoint["DENDRO_TS_TIME_END"]      = m_uiTinfo._m_uiTe;
                checkPoint["DENDRO_TS_ELEMENT_ORDER"] = m_uiElementOrder;

                checkPoint["DENDRO_TS_TIME_CURRENT"]    = m_uiTinfo._m_uiT;
                checkPoint["DENDRO_TS_STEP_CURRENT"]    = m_uiTinfo._m_uiStep;
                checkPoint["DENDRO_TS_TIME_STEP_SIZE"]  = m_uiTinfo._m_uiTh;
                checkPoint["DENDRO_TS_LAST_IO_TIME"]    = m_uiTinfo._m_uiT;

                checkPoint["DENDRO_TS_WAVELET_TOLERANCE"]  = NLSM_WAVELET_TOL;
                checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"] = NLSM_LOAD_IMB_TOL;
                checkPoint["DENDRO_TS_NUM_VARS"]=numVars; // number of variables to restore.
                checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"]=m_uiMesh->getMPICommSize(); // (note that rank 0 is always active).
                
                outfile<<std::setw(4)<<checkPoint<<std::endl;
                outfile.close();

            }

        }

        return 0;
    }

    int NLSMCtxGPU::restore_checkpt()
    {

        DVec& m_evar         = m_var[VL::CPU_EV];

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
                sprintf(fName,"%s_step_%d.cp",NLSM_CHKPT_FILE_PREFIX.c_str() , cpIndex);
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
                    
                    NLSM_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    NLSM_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
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
                sprintf(fName,"%s_step_%d.cp", NLSM_CHKPT_FILE_PREFIX.c_str(), restoreFileIndex);
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
                    
                    NLSM_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    NLSM_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
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
            par::Mpi_Bcast(&NLSM_WAVELET_TOL,1,0,comm);
            par::Mpi_Bcast(&NLSM_LOAD_IMB_TOL,1,0,comm);

            par::Mpi_Bcast(&numVars,1,0,comm);
            par::Mpi_Bcast(&m_uiElementOrder,1,0,comm);
            par::Mpi_Bcast(&activeCommSz,1,0,comm);
            
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

                sprintf(fName, "%s_octree_%d_%d.oct", NLSM_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readOctFromFile(fName, octree);
                assert(par::test::isUniqueAndSorted(octree, newComm));

            }

            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: octree (*.oct) restore file is corrupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            newMesh=new ot::Mesh(octree,1,m_uiElementOrder,activeCommSz,comm);
            // no need to transfer data only to resize the contex variables. 
            this->grid_transfer(newMesh);
            
            // only reads the evolution variables. 
            if(isActive) {

                int activeRank;
                int activeNpes;

                DendroScalar* inVec[NLSM_NUM_VARS];
                m_evar.to_2d(inVec);

                MPI_Comm_rank(newComm, &activeRank);
                MPI_Comm_size(newComm, &activeNpes);
                assert(activeNpes == activeCommSz);

                sprintf(fName,"%s_%d_%d.var",NLSM_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,inVec,NLSM_NUM_VARS);
            }

            MPI_Comm_free(&newComm);
            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: varible (*.var) restore file currupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            std::swap(m_uiMesh,newMesh);
            delete newMesh;
            
        unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
        unsigned int totalElems=0;
        par::Mpi_Allreduce(&localSz, &totalElems ,1,MPI_SUM,comm);
        
        if(!rank) std::cout<<" checkpoint at step : "<<m_uiTinfo._m_uiStep<<"active Comm. sz: "<<activeCommSz<<" restore successful: "<<" restored mesh size: "<<totalElems<<std::endl;

        m_uiIsETSSynced = false;
        return 0;

    }

    int NLSMCtxGPU::pre_timestep(DVec sIn)
    {
        return 0;
    }

    int NLSMCtxGPU::pre_stage(DVec  sIn)
    {
        return 0;
    }

    int NLSMCtxGPU::rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time)
    {
        if(!m_uiMesh->isActive())
            return 0;
        
        const std::vector<ot::Block> & blk_list = m_uiMesh->getLocalBlockList();
        const unsigned int nblocks = blk_list.size();
        const unsigned int NUM_BATCHES  = nblocks/DEVICE_RHS_BATCHED_GRAIN_SZ + 1;
        
        const unsigned int nx     = blk_list[0].getAllocationSzX();
        const unsigned int ny     = blk_list[0].getAllocationSzY();
        const unsigned int nz     = blk_list[0].getAllocationSzZ();
        const unsigned int BLK_SZ = nx* ny* nz;

        const unsigned int pencils = 13;
        const unsigned int pen_sz  = 13;
        

        DEVICE_REAL* const out_ptr = out->get_vec_ptr();
        const DEVICE_REAL * const in_ptr  = in->get_vec_ptr();
        

        DVec& m_dptr_uz_i = m_var[VL::GPU_EV_UZ_IN];
        DVec& m_dptr_uz_o = m_var[VL::GPU_EV_UZ_OUT];

        const unsigned int unzip_sz = m_uiMesh->getDegOfFreedomUnZip();

        //m_mesh_cpu.unzip_cg<DEVICE_REAL, cudaStream_t>(m_uiMesh, m_dptr_mesh, in_ptr, m_dptr_uz_i.get_vec_ptr(), in->get_dof(),0);
        this->unzip(*in,m_dptr_uz_i,NLSM_ASYNC_COMM_K);

        for(unsigned int bid =0; bid<NUM_BATCHES; bid++)
        {
            
            const unsigned int block_begin       = (bid * nblocks)/NUM_BATCHES;
            const unsigned int block_end         = ((bid+1) * nblocks)/NUM_BATCHES;
            const unsigned int numblocks         = block_end-block_begin;  
            //printf("batch %d  of NUM_BATCHES: %d\n",bid,NUM_BATCHES);
            //const unsigned int BATCHED_BLOCKS_SZ = numblocks;
            //printf( "batch id : %d begin =%d  end =%d batch size %d \n", bid, block_begin, block_end, numblocks);
            
            dim3 grid_x  = dim3(ny/pencils,nz,numblocks);
            dim3 block_x = dim3(nx,pencils,1);

            dim3 grid_y  = dim3(nx/pencils,nz,numblocks);
            dim3 block_y = dim3(pencils,ny,1);

            dim3 grid_z  = dim3(nx/pencils,ny,numblocks);
            dim3 block_z = dim3(pencils,nz,1);

            launch_dir_x_deriv_kernel<3,pencils,pen_sz, DEVICE_RHS_BATCHED_GRAIN_SZ> <<<grid_x, block_x>>> (m_dptr_mesh, m_dptr_uz_i.get_vec_ptr(), m_dptr_deriv_evars, m_mesh_cpu.m_blk_list, block_begin);
            GPUDevice::check_last_error();
            
            launch_dir_y_deriv_kernel<3,pencils,pen_sz, DEVICE_RHS_BATCHED_GRAIN_SZ> <<<grid_y, block_y>>> (m_dptr_mesh, m_dptr_uz_i.get_vec_ptr(), m_dptr_deriv_evars, m_mesh_cpu.m_blk_list, block_begin);
            GPUDevice::check_last_error();
            
            launch_dir_z_deriv_kernel<3,pencils,pen_sz, DEVICE_RHS_BATCHED_GRAIN_SZ> <<<grid_z, block_z>>> (m_dptr_mesh, m_dptr_uz_i.get_vec_ptr(), m_dptr_deriv_evars, m_mesh_cpu.m_blk_list, block_begin);
            GPUDevice::check_last_error();

            launch_rhs<3,pencils,pen_sz, DEVICE_RHS_BATCHED_GRAIN_SZ> <<<grid_x, block_x>>> (m_dptr_mesh, m_dptr_uz_o.get_vec_ptr(), m_dptr_uz_i.get_vec_ptr(), m_dptr_deriv_evars, m_mesh_cpu.m_blk_list, block_begin);
            GPUDevice::check_last_error();
            
            GPUDevice::device_synchronize();    
        }

        
        //m_mesh_cpu.zip_cg<DEVICE_REAL, cudaStream_t>(m_uiMesh,m_dptr_mesh,m_dptr_uz_o.get_vec_ptr(),out_ptr,out->get_dof(),0);
        this->zip(m_dptr_uz_o,*out);

        return 0;
        
    }

    int NLSMCtxGPU::post_stage(DVec sIn)
    { 
        return 0;
    }

    int NLSMCtxGPU::post_timestep(DVec sIn)
    {
        return 0;
    }

    bool NLSMCtxGPU::is_remesh()
    {

        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::IS_REMESH].start();
        #endif

        bool isRefine = false;
        DVec& m_evar         = m_var[VL::CPU_EV];
        DVec& m_evar_unzip   = m_var[VL::CPU_EV_UZ];

        if(NLSM_ENABLE_BLOCK_ADAPTIVITY)
            return isRefine;
        
        MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();
        this->unzip(m_evar,m_evar_unzip,NLSM_ASYNC_COMM_K);

        DendroScalar* unzipVar[NLSM_NUM_VARS];
        m_evar_unzip.to_2d(unzipVar);

        unsigned int refineVarIds[NLSM_NUM_REFINE_VARS];
        for(unsigned int vIndex=0; vIndex<NLSM_NUM_REFINE_VARS; vIndex++)
            refineVarIds[vIndex]=NLSM_REFINE_VARIABLE_INDICES[vIndex];

        double wTol=NLSM_WAVELET_TOL;
        std::function<double(double,double,double,double*)> waveletTolFunc =[wTol](double x,double y, double z,double*hx) {
            return computeWTol(x,y,z,hx);
        };


        isRefine=m_uiMesh->isReMeshUnzip((const double **)unzipVar,refineVarIds,NLSM_NUM_REFINE_VARS,waveletTolFunc,NLSM_DENDRO_AMR_FAC); 

        return isRefine;

        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::IS_REMESH].stop();
        #endif

    }

    DVec& NLSMCtxGPU::get_evolution_vars()
    {
        //return m_var[GPU_EV_DG];
        return m_var[GPU_EV];
    }
    
    // DVec& NLSMCtxGPU::get_constraint_vars()
    // {
    //     return m_cvar;
    // }

    // DVec& NLSMCtxGPU::get_primitive_vars()
    // {
    //     return m_pvar;
    // }

    int NLSMCtxGPU::terminal_output()
    {
        if(m_uiMesh->isActive())
        {
            DVec& m_evar         = m_var[VL::CPU_EV];
        
            for(unsigned int v=0; v < m_evar.get_dof(); v++)
            {
                DendroScalar min=0, max=0;
                min=vecMin(m_uiMesh,m_evar.get_vec_ptr(),ot::VEC_TYPE::CG_NODAL,true);
                max=vecMax(m_uiMesh,m_evar.get_vec_ptr(),ot::VEC_TYPE::CG_NODAL,true);
                if(!(m_uiMesh->getMPIRank()))
                {
                    std::cout<<"[NLSMCtx]: step : "<<m_uiTinfo._m_uiStep<<"\ttime : "<<m_uiTinfo._m_uiT<<std::endl;
                    std::cout<<"[NLSMCtx]:  "<<NLSM_VAR_NAMES[v]<<" (min,max) : \t ( "<<min<<", "<<max<<" ) "<<std::endl;
                }
                   

            }


            
            #ifdef NLSM_COMPARE_WITH_ANALYTICAL_SOL
                double* evolVar[2];
                m_uiEVar.Get2DVec(evolVar);
                double * chiAnalytical=m_uiMesh->createVector<double>();
                double * diffVec=m_uiMesh->createVector<double>();

                std::function<void(double,double,double,double,double*)> u_x_t=[](double x,double y,double z,double t,double*var){nlsm::analyticalSol(x,y,z,t,var);};

                // initialize diff begin.
                unsigned int nodeLookUp_CG;
                unsigned int nodeLookUp_DG;
                double x,y,z,len;
                const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
                unsigned int ownerID,ii_x,jj_y,kk_z;
                unsigned int eleOrder=m_uiMesh->getElementOrder();
                const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
                const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
                const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
                const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
                const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();


                double var[2];

                double mp, mm, mp_adm, mm_adm, E, J1, J2, J3;

                for(unsigned int elem=m_uiMesh->getElementLocalBegin();elem<m_uiMesh->getElementLocalEnd();elem++)
                {


                    for(unsigned int k=0;k<(eleOrder+1);k++)
                        for(unsigned int j=0;j<(eleOrder+1);j++ )
                            for(unsigned int i=0;i<(eleOrder+1);i++)
                            {
                                nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                                {
                                    nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                    m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                                    len= (double) (1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                                    x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                                    y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                                    z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));
                                    
                                    u_x_t((double)x,(double)y,(double)z,m_uiTinfo._m_uiT,var);
                                    diffVec[nodeLookUp_CG]=var[nlsm::VAR::U_CHI]-evolVar[nlsm::VAR::U_CHI][nodeLookUp_CG];
                                    chiAnalytical[nodeLookUp_CG]=var[nlsm::VAR::U_CHI];


                                }

                            }

                }

                m_uiMesh->performGhostExchange(diffVec);
                m_uiMesh->performGhostExchange(chiAnalytical);

                double l_rs = rsNormLp<double>(m_uiMesh,diffVec,2);
                double l_min=vecMin(diffVec+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                double l_max=vecMax(diffVec+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                double l2_norm=normL2(diffVec+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                DendroIntL local_dof=m_uiMesh->getNumLocalMeshNodes();
                DendroIntL total_dof=0;
                par::Mpi_Reduce(&local_dof,&total_dof,1,MPI_SUM,0,m_uiMesh->getMPICommunicator());


                if(!m_uiMesh->getMPIRank()) {
                    //std::cout << "executing step: " << m_uiCurrentStep << " dt: " << m_uiT_h << " rk_time : "<< m_uiCurrentTime << std::endl;
                    l2_norm=sqrt((l2_norm*l2_norm)/(double)(total_dof*total_dof));
                    std::cout <<YLW<< "\t ||VAR::DIFF|| (min, max,l2,l_2rs) : ("<<l_min<<", "<<l_max<<", "<<l2_norm<<", "<<l_rs<<" ) "<<NRM<<std::endl;

                    std::ofstream fileGW;
                    char fName[256];
                    sprintf(fName,"%s_error.dat",nlsm::NLSM_PROFILE_FILE_PREFIX.c_str());
                    fileGW.open (fName,std::ofstream::app);
                    // writes the header
                    if(m_uiTinfo._m_uiStep==0)
                        fileGW<<"TimeStep\t"<<" time\t"<<" min\t"<<" max\t"<<" l2\t cgNodes\t l2_rs"<<std::endl;

                    fileGW<<m_uiTinfo._m_uiStep<<"\t"<<m_uiTinfo._m_uiT<<"\t"<<l_min<<"\t"<<l_max<<"\t"<<l2_norm<<"\t"<<total_dof<<"\t "<<l_rs<<std::endl;
                    fileGW.close();


                }
                // initialize diff end
                delete [] chiAnalytical;
                delete [] diffVec;
            #endif

            
        }

        return 0; 
    }

    int NLSMCtxGPU::grid_transfer(const ot::Mesh* m_new)
    {
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].start();
        #endif
        DVec& m_evar  = m_var[VL::CPU_EV];
        DVec::grid_transfer(m_uiMesh, m_new, m_evar);

        ot::dealloc_mpi_ctx<DendroScalar>(m_uiMesh,m_mpi_ctx,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);
        ot::alloc_mpi_ctx<DendroScalar>(m_new,m_mpi_ctx,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);

        device::dealloc_mpi_ctx<DendroScalar>(m_uiMesh,m_mpi_ctx_device,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);
        device::alloc_mpi_ctx<DendroScalar>(m_new,m_mpi_ctx_device,NLSM_NUM_VARS,NLSM_ASYNC_COMM_K);
        //printf("igt ended\n");

        //m_var[VL::CPU_EV_DG].destroy_vector();
        m_var[VL::CPU_EV_UZ].destroy_vector();

        //m_var[VL::GPU_EV_DG].destroy_vector();
        m_var[VL::GPU_EV].destroy_vector();
        m_var[VL::GPU_EV_UZ_IN].destroy_vector();
        m_var[VL::GPU_EV_UZ_OUT].destroy_vector();

        //m_var[VL::CPU_EV_DG].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_NODES,ot::DVEC_LOC::HOST,NLSM_NUM_VARS,true);
        m_var[VL::CPU_EV_UZ].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::HOST,NLSM_NUM_VARS,true);
        //m_var[VL::GPU_EV_DG].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_NODES,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);
        
        m_var[VL::GPU_EV].create_vector(m_new,ot::DVEC_TYPE::OCT_SHARED_NODES,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_IN].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);
        m_var[VL::GPU_EV_UZ_OUT].create_vector(m_new,ot::DVEC_TYPE::OCT_LOCAL_WITH_PADDING,ot::DVEC_LOC::DEVICE,NLSM_NUM_VARS,true);
        
        this->host_to_device_sync();
        //printf("hto d ended\n");
        m_uiIsETSSynced=false;
        
        #ifdef __PROFILE_CTX__
            m_uiCtxpt[ts::CTXPROFILE::GRID_TRASFER].stop();
        #endif
        return 0 ;
    }

}