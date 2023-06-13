/**
 * @file run_meshgpu_tests.cu
 * @brief briefly test meshGPU utilities
 * @version 0.1
 * @date 2022-01-21
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include "dendro.h"
#include "mpi.h"
#include "profiler.h"
#include "meshTestUtils.h"
#include "mesh_gpu.cuh"
#include "device.h"
#include "meshUtils.h"
#include "device_utils.cuh"

CONST_MEM DEVICE_REAL device::refel_1d[2*REFEL_CONST_MEM_MAX];
enum event{unzip_cg_g, unzip_cg_s, unzip_dg_g, unzip_dg_s, zip_cg,zip_dg, gpu_mesh,gpu_unzip, gpu_zip, last};
const char*event_names[]= {"unzip_cg (gather)", "unzip_cg (scatter)", "unzip_dg (gather)", "unzip_dg (scatter)", "zip_cg", "zip_dg", "gpu_mesh", "gpu_unzip", "gpu_zip"};

static inline void cuda_check_last_error()
{
    //#if defined(DEBUG) || defined(_DEBUG)
       cudaError_t err = cudaGetLastError();
       if (err != cudaSuccess) 
        printf("%ss\n", cudaGetErrorString(err));
    //#endif
    return ;
}

int main(int argc, char** argv)
{   
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    std::vector<profiler_t> time_events;
    time_events.resize(event::last);

    int devicesCount;
    cudaGetDeviceCount(&devicesCount);

    if(!rank)
        printf("number of cuda devices: %d\n",devicesCount);

    cudaSetDevice(rank%devicesCount);
    cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);

    for(unsigned int i=0; i < time_events.size(); i++)
        time_events[i].clear();

    if(argc<5)
    {
        if(!rank) std::cout<<"Usage: "<<argv[0]<<" maxDepth wavelet_tol partition_tol eleOrder run_gpu"<<std::endl;
        MPI_Abort(comm,0);
    }

    m_uiMaxDepth=atoi(argv[1]);
    double wavelet_tol=atof(argv[2]);
    double partition_tol=atof(argv[3]);
    unsigned int eOrder=atoi(argv[4]);
    unsigned int run_gpu=atoi(argv[5]);
    const unsigned int iter=2;

    if(!rank)
    {
        std::cout<<YLW<<"maxDepth: "<<m_uiMaxDepth<<NRM<<std::endl;
        std::cout<<YLW<<"wavelet_tol: "<<wavelet_tol<<NRM<<std::endl;
        std::cout<<YLW<<"partition_tol: "<<partition_tol<<NRM<<std::endl;
        std::cout<<YLW<<"eleOrder: "<<eOrder<<NRM<<std::endl;
        std::cout<<YLW<<"OMP_NUM_THREADS:"<<std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1)<<std::endl;;

    }

    _InitializeHcurve(m_uiDim);

    // function that we need to interpolate.
    const double d_min=-10;
    const double d_max=10;
    double dMin[]={d_min,d_min,d_min};
    double dMax[]={d_max,d_max,d_max};

    Point pt_min(d_min,d_min,d_min);
    Point pt_max(d_max,d_max,d_max);
    //@note that based on how the functions are defined (f(x), dxf(x), etc) the compuatational domain is equivalent to the grid domain.
    
    const unsigned int dof=24;
    const double unzip_check_tol=1e-4;
    bool valid=false;

    std::function<void(double,double,double,double*)> func =[d_min,d_max,dof](double x,double y,double z, double * var){

        double ca[]={-2,0,0};
        double cb[]={2,0,0};
        double rra = ((x-ca[0]) * (x-ca[0]) + (y-ca[1])* (y-ca[1]) + (z-ca[2])* (z-ca[2]))/1.0;
        double rrb = ((x-cb[0]) * (x-cb[0]) + (y-cb[1])* (y-cb[1]) + (z-cb[2])* (z-cb[2]))/1.0;
        for (unsigned int v=0; v< dof; v++)
            var[v] = exp(-rra) + exp(-rrb); 
        
        //(sin(2*M_PI*xx)*sin(2*M_PI*yy)*sin(2*M_PI*zz));
        return;
         
    };

    std::vector<ot::TreeNode> tmpNodes;
    std::function<double(double,double,double)> fr = [func,d_min,d_max,dof](double x,double y,double z){
        double xx = (x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double yy = (y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double zz = (z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min ;
        double var[dof]; 
        func(xx,yy,zz,var); 
        return var[0];
    };
    function2Octree(fr,tmpNodes,m_uiMaxDepth,wavelet_tol,eOrder,comm);
    //createRegularOctree(tmpNodes,4,3,m_uiMaxDepth,comm);
    ot::Mesh * mesh= ot::createMesh(tmpNodes.data(),tmpNodes.size(),eOrder,comm,1,ot::SM_TYPE::FDM,DENDRO_DEFAULT_GRAIN_SZ,partition_tol,DENDRO_DEFAULT_SF_K);
    mesh->setDomainBounds(pt_min, pt_max);

    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);
    if(!rank)
        printf("created mesh with lev_min = %d and lev_max=%d\n",lmin,lmax);
    
    
    double * u = mesh->createCGVector<double>(func,dof);
    double * v_dg    = mesh->createDGVector(func,dof);
    double * u_unzip = mesh->createUnZippedVector<double>(dof);
    double * v_unzip = mesh->createUnZippedVector<double>(dof);

    if(run_gpu==1)
    {
        time_events[event::gpu_mesh].start();
        device::MeshGPU mesh_gpu;
        device::MeshGPU* dptr_mesh = mesh_gpu.alloc_mesh_on_device(mesh);
        time_events[event::gpu_mesh].stop();

        double * dptr_vcg     = mesh_gpu.createVector        <double,device::vec_type::device>(dof);
        double * dptr_vdg     = mesh_gpu.createDGVector      <double,device::vec_type::device>(dof);
        
        GPUDevice::host_to_device<DEVICE_REAL>(u   ,dptr_vcg,mesh->getDegOfFreedom()*dof);
        GPUDevice::host_to_device<DEVICE_REAL>(v_dg,dptr_vdg,mesh->getDegOfFreedomDG()*dof);
        
        
        double * dptr_u_unzip = mesh_gpu.createUnZippedVector<double,device::vec_type::device>(dof); // cg unzip
        double * dptr_v_unzip = mesh_gpu.createUnZippedVector<double,device::vec_type::device>(dof); // dg unzip
        
        std::vector<ot::AsyncExchangeContex> ctx_host;
        std::vector<ot::AsyncExchangeContex> ctx_device;

        ot::alloc_mpi_ctx<double>(mesh,ctx_host,dof,1);
        device::alloc_mpi_ctx<double>(mesh,ctx_device,dof,1);

        mesh_gpu.read_from_ghost_cg_begin<double,cudaStream_t>(ctx_host[0], ctx_device[0], mesh,dptr_mesh,dptr_vcg,dof);
        cuda_check_last_error();
        mesh_gpu.read_from_ghost_cg_end<double,cudaStream_t>(ctx_host[0], ctx_device[0], mesh, dptr_mesh, dptr_vcg,dof);
        cuda_check_last_error();
        
        mesh_gpu.read_from_ghost_dg_begin<DEVICE_REAL,cudaStream_t>(ctx_host[0], ctx_device[0], mesh,dptr_mesh,dptr_vdg,dof);
        cuda_check_last_error();
        mesh_gpu.read_from_ghost_dg_end<double,cudaStream_t>(ctx_host[0], ctx_device[0], mesh, dptr_mesh, dptr_vdg,dof);
        cuda_check_last_error();


        for(unsigned int w=0; w< 1; w++)
        {
            mesh_gpu.unzip_cg(mesh,dptr_mesh,dptr_vcg,dptr_u_unzip,dof,(cudaStream_t)0);
            GPUDevice::device_synchronize();
            cuda_check_last_error();

            mesh_gpu.unzip_dg(mesh,dptr_mesh,dptr_vdg,dptr_v_unzip,dof,(cudaStream_t)0);
            GPUDevice::device_synchronize();
            cuda_check_last_error();

            GPUDevice::device_to_host<DEVICE_REAL>(u_unzip,dptr_u_unzip,mesh->getDegOfFreedomUnZip()*dof);
            GPUDevice::device_to_host<DEVICE_REAL>(v_unzip,dptr_v_unzip,mesh->getDegOfFreedomUnZip()*dof);
            for(unsigned int var=0; var < std::min(3u,dof); var++)
            {
                valid = ot::test::isUnzipValid(mesh,u_unzip + var * mesh->getDegOfFreedomUnZip(), fr, unzip_check_tol);
                
                if(!rank)
                    std::cout<<"GPU unzip valid cg scatter: "<<valid<<"for dof: "<<var<<std::endl;
            }
            
            for(unsigned int var=0; var < std::min(3u,dof); var++)
            {
                valid = ot::test::isUnzipValid(mesh,v_unzip + var * mesh->getDegOfFreedomUnZip(), fr, unzip_check_tol);
                if(!rank)
                    std::cout<<"GPU unzip valid dg scatter: "<<valid<<"for dof: "<<var<<std::endl;
            }
            

            

            mesh_gpu.zip_cg(mesh,dptr_mesh,dptr_u_unzip,dptr_vcg,dof,(cudaStream_t)0);
            mesh_gpu.zip_dg(mesh,dptr_mesh,dptr_v_unzip,dptr_vdg,dof,(cudaStream_t)0);
            GPUDevice::device_to_host<DEVICE_REAL>(u   , dptr_vcg,mesh->getDegOfFreedom()  *dof);
            GPUDevice::device_to_host<DEVICE_REAL>(v_dg, dptr_vdg,mesh->getDegOfFreedomDG()*dof);
            cuda_check_last_error();

        }
        

        time_events[event::gpu_unzip].clear();
        time_events[event::gpu_unzip].start();
        for(unsigned int i=0; i < iter;i++)
        {
            mesh_gpu.unzip_cg(mesh,dptr_mesh,dptr_vcg,dptr_u_unzip,dof,(cudaStream_t)0);
            //mesh_gpu.unzip_dg(mesh,dptr_mesh,dptr_vdg,dptr_v_unzip,dof,(cudaStream_t)0);
            GPUDevice::device_synchronize();
            cuda_check_last_error();
        }
        time_events[event::gpu_unzip].stop();

        time_events[event::gpu_zip].clear();
        time_events[event::gpu_zip].start();
        for(unsigned int i=0; i < iter;i++)
        {
            mesh_gpu.zip_dg(mesh,dptr_mesh,dptr_v_unzip,dptr_vdg,dof,(cudaStream_t)0);
            GPUDevice::device_synchronize();
            cuda_check_last_error();
        }
        time_events[event::gpu_zip].stop();

        mesh_gpu.destroyVec<double,device::vec_type::device>(dptr_vdg);
        mesh_gpu.destroyVec<double,device::vec_type::device>(dptr_vcg);
        mesh_gpu.destroyVec<double,device::vec_type::device>(dptr_v_unzip);
        mesh_gpu.destroyVec<double,device::vec_type::device>(dptr_u_unzip);


        mesh_gpu.dealloc_mesh_on_device(dptr_mesh);

    }else
    {

        mesh->readFromGhostBegin(u,dof);
        mesh->readFromGhostEnd(u,dof);

        mesh->readFromGhostBeginEleDGVec(v_dg,dof);
        mesh->readFromGhostEndEleDGVec(v_dg,dof);
        std::vector<unsigned int > blkIDs;
        blkIDs.resize(mesh->getLocalBlockList().size());

        for(unsigned int i=0; i< blkIDs.size(); i++)
            blkIDs[i] = i; 
        // unzip all the blocks. 
        //this->unzip(in,out,blkIDs.data(),blkIDs.size(),dof);

        mesh->unzip_scatter(u,u_unzip,dof);
        for(unsigned int var=0; var < std::min(3u,dof); var++)
        {
            valid = ot::test::isUnzipValid(mesh,u_unzip + var * mesh->getDegOfFreedomUnZip(), fr, unzip_check_tol);
            if(!rank)
                std::cout<<"unzip valid scatter: "<<valid<<"for dof: "<<var<<std::endl;
        }
        
        mesh->unzipDG_scatter(v_dg,v_unzip,dof);
        for(unsigned int var=0; var < std::min(3u,dof); var++)
        {
            valid = ot::test::isUnzipValid(mesh,v_unzip + var * mesh->getDegOfFreedomUnZip(), fr, unzip_check_tol);
            if(!rank)
                std::cout<<"unzip valid dg scatter: "<<valid<<"for dof: "<<var<<std::endl;
        }

        //mesh->unzip(u,u_unzip,dof);
        mesh->unzip(u,u_unzip,blkIDs.data(),blkIDs.size(),dof);
        for(unsigned int var=0; var < std::min(3u,dof); var++)
        {
            valid = ot::test::isUnzipValid(mesh,u_unzip + var * mesh->getDegOfFreedomUnZip(), fr, unzip_check_tol);
            if(!rank)
                std::cout<<"unzip valid gather: "<<valid<<"for dof: "<<var<<std::endl;
        }

        mesh->unzipDG(v_dg,v_unzip,dof);
        for(unsigned int var=0; var < std::min(3u,dof); var++)
        {
            valid = ot::test::isUnzipValid(mesh,v_unzip + var * mesh->getDegOfFreedomUnZip(), fr, unzip_check_tol);
            if(!rank)
                std::cout<<"unzip valid dg gather: "<<valid<<"for dof: "<<var<<std::endl;
        }
        

        time_events[event::unzip_cg_s].clear();
        time_events[event::unzip_cg_s].start();
        for(unsigned int i=0; i < iter;i++)
            mesh->unzip_scatter(u,u_unzip,dof);
        time_events[event::unzip_cg_s].stop();

        time_events[event::unzip_dg_s].clear();
        time_events[event::unzip_dg_s].start();
        for(unsigned int i=0; i < iter;i++)
            mesh->unzipDG_scatter(v_dg,v_unzip,dof);
        time_events[event::unzip_dg_s].stop();
        
        time_events[event::unzip_cg_g].clear();
        time_events[event::unzip_cg_g].start();
        for(unsigned int i=0; i < iter;i++)
            mesh->unzip(u,u_unzip,blkIDs.data(),blkIDs.size(),dof);//mesh->unzip(u,u_unzip,dof);
        time_events[event::unzip_cg_g].stop();

        time_events[event::unzip_dg_g].clear();
        time_events[event::unzip_dg_g].start();
        for(unsigned int i=0; i < iter;i++)
            mesh->unzipDG(v_dg,v_unzip,dof);
        time_events[event::unzip_dg_g].stop();

        time_events[event::zip_cg].clear();
        time_events[event::zip_cg].start();
        for(unsigned int i=0; i < iter;i++)
            for(unsigned int var=0; var< dof; var++)
                mesh->zip(u_unzip + var * mesh->getDegOfFreedomUnZip(),u + var * mesh->getDegOfFreedom());
        time_events[event::zip_cg].stop();

    }
    



    
    
    
    double t_local;
    double t_stat[3];
    
    for(unsigned int i=0; i < time_events.size(); i++)
    {
        t_local=time_events[i].seconds;
        par::computeOverallStats(&t_local,t_stat,comm,event_names[i]);

    }

    mesh->destroyVector(u);
    mesh->destroyVector(v_dg);

    mesh->destroyVector(u_unzip);
    mesh->destroyVector(v_unzip);

    delete mesh;
    MPI_Finalize();
    return 0;
}
