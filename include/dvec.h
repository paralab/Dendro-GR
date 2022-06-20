/**
 * @file dvec.h
 * @author Milinda Fernando
 * @brief Vector class with additional auxiliary information. 
 * @version 0.1
 * @date 2019-12-19
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once
#include <iostream>
#include "mpi.h"
#include "mesh.h"
#include "mathUtils.h"
#include "device.h"
namespace ot
{
    enum DVEC_TYPE  {OCT_SHARED_NODES=0, OCT_LOCAL_NODES, OCT_LOCAL_WITH_PADDING, OCT_CELL_CENTERED};
    enum DVEC_LOC   {HOST=0, DEVICE};

    template<typename T,typename I>
    class DVector
    {   
        protected:

            /**@brief ptr for the data*/
            T* m_data_ptr = NULL;

            /**@brief size of the vector*/
            I m_size = 0 ;

            /**@brief : true if allocated with ghost/halo regions, false otherwise */
            bool m_ghost_allocated = true;

            /**@brief: number of degrees of freedoms. */
            unsigned int m_dof;

            /**@brief : allocated vector type*/
            DVEC_TYPE  m_vec_type;

            /**@brief : where the vector is allocated*/
            DVEC_LOC   m_vec_loc;

            /**@brief: MPI Communicator for the vector*/
            MPI_Comm m_comm;
            
        public:
            
            /**@brief: default constructor*/
            DVector();
            
            /**@brief: destructor*/
            ~DVector();

            /**
             * @brief Construct a new vec create object
             * 
             * @param pMesh : pointer to underlying ot::Mesh object
             * @param isGhosted : true if need to create a vector with halo regions. 
             * @param isUnzip : true if this is unzip vector format (block local vector format)
             * @param isElemental : true if this is elemental vector
             * @param dof : number of degrees of freedoms. 
             */
            void create_vector(const ot::Mesh* pMesh, DVEC_TYPE type, DVEC_LOC loc, unsigned int dof=1, bool allocate_ghost=true);

            /**@brief creates a similar vector as dvec*/
            void create_vector(const ot::DVector<T,I>& dvec);

            /**
             * @brief deallocates the vector object
             */
            void destroy_vector();  

            /**@brief: returns the vec pointer*/
            inline T* get_vec_ptr() { return m_data_ptr; };
            
            /**
             * @brief update the vector pointer object. 
             * @param vec : vec pointer
             */
            inline void restore_vec_ptr(T* vec) { m_data_ptr = vec;};

            /**
             * @brief : equal operator for the DVector
             * 
             * @param other : DVector object 
             * @return true if points to the same pointer. 
             * @return false otherwise. 
             */
            bool operator==(DVector<T,I>& other) const { return ( m_data_ptr == (other.m_data_ptr)); }

            /**@brief: returns the number of degrees of freedoms. */
            inline unsigned int get_dof() const {return m_dof;}
            
            inline MPI_Comm get_communicator() const {return m_comm;}
            inline DVEC_TYPE get_type() const {return m_vec_type;}
            inline DVEC_LOC get_loc()  const {return m_vec_loc;}
            inline bool is_ghost_allocated() const {return m_ghost_allocated;}
            /**
             * @brief Get the Size vector
             * @return I 
             */
            inline I get_size() const { return m_size;}

            /**
             * @brief returns the 2D vector. Note that the v2d assumeed to allocated no need to be on the heap;
             * @param v2d 
             */
            void to_2d(T** v2d);

            /**@brief: computes the min and max of the vector. */
            //void min_max(const ot::Mesh* pMesh, T& min, T& max, unsigned int dof=0);
            
            void copy_data(const DVector<T,I>& v);



            static void axpy(const ot::Mesh* const pMesh, T a, const DVector<T,I>& x, DVector<T,I>& y, bool local_only=true);
            
            static void axpby(const ot::Mesh* const pMesh,T a, const DVector<T,I>& x, T b, DVector<T,I>& y, bool local_only=true);

            static void grid_transfer(ot::Mesh* m_old, const ot::Mesh* m_new, DVector<T,I> & v);
            
            
            
    };

    template<typename T,typename I>
    DVector<T,I>::DVector()
    {
        m_data_ptr = NULL;
        m_size     = 0;
        m_dof      = 0;
        m_comm     = MPI_COMM_NULL;

    }

    template<typename T,typename I>
    DVector<T,I>::~DVector() {}
    
    template<typename T,typename I>
    void DVector<T,I>::create_vector(const ot::Mesh* pMesh, DVEC_TYPE type, DVEC_LOC loc, unsigned int dof, bool allocate_ghost)
    {
        m_dof       = dof;
        m_comm      = pMesh->getMPICommunicator();
        m_vec_type  = type;
        m_vec_loc   = loc;
        m_ghost_allocated = allocate_ghost;
        m_size     = 0;

        if(!(pMesh->isActive()))
            return;

        if(m_vec_type == DVEC_TYPE::OCT_SHARED_NODES)
            (allocate_ghost) ? m_size = pMesh->getDegOfFreedom() * m_dof : m_size = pMesh->getNumLocalMeshNodes() * m_dof;
        else if (m_vec_type == DVEC_TYPE::OCT_LOCAL_NODES)
            (allocate_ghost) ? m_size = pMesh->getDegOfFreedomDG() * m_dof : m_size = pMesh->getNumLocalMeshElements() * pMesh->getNumNodesPerElement() * m_dof;
        else if (m_vec_type == DVEC_TYPE::OCT_LOCAL_WITH_PADDING)
            (allocate_ghost) ? m_size = pMesh->getDegOfFreedomUnZip() * m_dof : m_size = pMesh->getDegOfFreedomUnZip() * m_dof;
        else if (m_vec_type == DVEC_TYPE::OCT_CELL_CENTERED)
            (allocate_ghost) ? m_size = pMesh->getAllElements().size() * m_dof : m_size = pMesh->getNumLocalMeshElements() * m_dof;
        else
        {
            dendro_log(" unknown type in DVector");
            MPI_Abort(m_comm,0);
        }

        if(m_vec_loc == DVEC_LOC::HOST)
        {
            #ifdef __CUDACC__
                m_data_ptr = GPUDevice::host_malloc<T>(m_size);
            #else 
                m_data_ptr = (T*) malloc(sizeof(T)*m_size);
            #endif
            
        }else if(m_vec_loc == DVEC_LOC::DEVICE)
        {   
            #ifdef __CUDACC__
                //#pragma message("dvec.h compiling with cuda")
                m_data_ptr = GPUDevice::device_malloc<T>(m_size);
            #endif

        }else
        {
            dendro_log(" unknown vector allocation location specified");
            MPI_Abort(m_comm,0);
        }

    }

    template<typename T,typename I>
    void DVector<T,I>::destroy_vector()
    {
        if(m_data_ptr == nullptr)
            return;

        if(m_vec_loc == DVEC_LOC::HOST)
        {
            #ifdef __CUDACC__
                GPUDevice::host_free<T>(m_data_ptr);
            #else 
               free(m_data_ptr);
            #endif

        }else if(m_vec_loc == DVEC_LOC::DEVICE)
        {
            #ifdef __CUDACC__
                GPUDevice::device_free<T>(m_data_ptr);
            #endif
        }else
        {
            dendro_log(" unknown vector deallocation location specified");
            MPI_Abort(m_comm,0);
        }

        m_data_ptr = nullptr;
        m_size     = 0;
        m_dof      = 0;
        
    }

    template<typename T, typename I>
    void DVector<T,I>::to_2d(T** v2d)
    {   
        if(m_data_ptr==nullptr)
            return;

        assert( (m_size % m_dof) == 0);
        const I sz_per_dof = m_size / m_dof;
        
        for(unsigned int i=0; i< m_dof; i++)
            v2d[i] = m_data_ptr + i*sz_per_dof;

        return;

    }
    
    template<typename T, typename I >
    void DVector<T,I>::copy_data(const DVector<T,I>& v)
    {
        if(m_data_ptr==nullptr)
            return;

        m_size            = v.m_size;
        m_dof             = v.m_dof;
        m_vec_type        = v.m_vec_type;
        m_vec_loc         = v.m_vec_loc;
        m_ghost_allocated = v.m_ghost_allocated;
        m_comm            = v.m_comm;
        
        T* dptr = v.m_data_ptr;

        if(v.m_vec_loc == DVEC_LOC::HOST)
            std::memcpy(m_data_ptr, dptr, sizeof(T)*m_size);
        else if(v.m_vec_loc == DVEC_LOC::DEVICE)
        {
            #ifdef __CUDACC__
                GPUDevice::check_error(cudaMemcpy(m_data_ptr,v.m_data_ptr,sizeof(T)*m_size,cudaMemcpyDeviceToDevice));
            #endif
        }

        
        return;
            
    }

    template<typename T, typename I >
    void DVector<T,I>::axpy(const ot::Mesh* const pMesh,  T a, const DVector<T,I>& x, DVector<T,I>& y, bool local_only)
    {
        if(y.m_data_ptr==nullptr)
            return;

        const T* const x_ptr = x.m_data_ptr;
        T* const y_ptr       = y.m_data_ptr;

        if(x.m_vec_loc == DVEC_LOC::HOST)
        {
            if(!local_only)
            {
                for(unsigned int i=0; i < x.m_size; i++)
                    y_ptr[i] += a* x_ptr[i];

                
            }else
            {
                const unsigned int sz_dof = x.m_size / x.m_dof;
                const unsigned int npe    = pMesh->getNumNodesPerElement(); 
                if(x.m_vec_type == DVEC_TYPE::OCT_SHARED_NODES)
                {
                    for(unsigned int v=0; v < x.m_dof; v++)
                        for(unsigned int node = pMesh->getNodeLocalBegin(); node < pMesh->getNodeLocalEnd(); node++)
                            y_ptr[v * sz_dof + node] += a* x_ptr[v * sz_dof + node];
                }else if(x.m_vec_type == DVEC_TYPE::OCT_LOCAL_NODES)
                {  
                    for(unsigned int v=0; v < x.m_dof; v++)
                        for(unsigned int node = pMesh->getElementLocalBegin() * npe ; node < pMesh->getElementLocalEnd() * npe; node++)
                            y_ptr[v * sz_dof + node] += a* x_ptr[v * sz_dof + node];
                }
                
            }

        }else if(x.m_vec_loc == DVEC_LOC::DEVICE)
        {
            #ifdef __CUDACC__
                if(x.m_dof==0 || x.m_size==0)
                    return;
                const unsigned int lb  = pMesh->getNodeLocalBegin();
                const unsigned int le  = pMesh->getNodeLocalEnd();
                const unsigned int szpdof = x.m_size/x.m_dof;
                //axpy_cu<<<x.m_size/1024 + 1, 1024>>>(x.m_size,a, x.m_data_ptr ,y.m_data_ptr);
                dim3 gb=dim3((le-lb)/GPU_MAX_THREADS_PER_BLOCK + 1, x.m_dof,1);
                dim3 tb=dim3(GPU_MAX_THREADS_PER_BLOCK,1,1);
                axpy_cu_2d<<< gb, tb>>>(lb, le, szpdof,a, x.m_data_ptr, y.m_data_ptr);
                GPUDevice::device_synchronize();
                GPUDevice::check_last_error();
            #endif

        }
        
            

    }
            
    template<typename T, typename I >
    void DVector<T,I>::axpby(const ot::Mesh* const pMesh, T a, const DVector<T,I>& x, T b, DVector<T,I>& y, bool local_only)
    {
        if(y.m_data_ptr==nullptr)
            return;

        const T* const x_ptr = x.m_data_ptr;
        T* const y_ptr       = y.m_data_ptr;

        if(!local_only)
        {
            for(unsigned int i=0; i < x.m_size; i++)
                y_ptr[i] = a * x_ptr[i] + b * y_ptr[i];
            
            
        }else
        {
            const unsigned int sz_dof = x.m_size / x.m_dof;
            const unsigned int npe    = pMesh->getNumNodesPerElement(); 
            if(x.m_vec_type == DVEC_TYPE::OCT_SHARED_NODES)
            {
                for(unsigned int v=0; v < x.m_dof; v++)
                    for(unsigned int node = pMesh->getNodeLocalBegin(); node < pMesh->getNodeLocalEnd(); node++)
                        y_ptr[v * sz_dof + node] = a* x_ptr[v * sz_dof + node]  + b * y_ptr[v * sz_dof + node];
            }else if(x.m_vec_type == DVEC_TYPE::OCT_LOCAL_NODES)
            {  
                for(unsigned int v=0; v < x.m_dof; v++)
                    for(unsigned int node = pMesh->getElementLocalBegin() * npe ; node < pMesh->getElementLocalEnd() * npe; node++)
                        y_ptr[v * sz_dof + node] = a* x_ptr[v * sz_dof + node] + b * y_ptr[v * sz_dof + node];
            }
        }
         
    }
    
    template<typename T, typename I >
    void DVector<T,I>::grid_transfer(ot::Mesh* m_old, const ot::Mesh* m_new, DVector<T,I> & dvec){

        ot::DVector<T,I> vec_tmp = ot::DVector<T,I>();
        vec_tmp.create_vector(m_new, dvec.get_type(), dvec.get_loc(), dvec.get_dof(),dvec.is_ghost_allocated());

        const unsigned int dof = dvec.get_dof();

        T* in  = dvec.get_vec_ptr();
        T* out = vec_tmp.get_vec_ptr();

        const unsigned int sz_per_dof_old =  (dof!=0) ? dvec.get_size()/dof  : 0;
        const unsigned int sz_per_dof_new =  (dof!=0) ? vec_tmp.get_size()/dof  : 0;
        
        if(dvec.m_vec_type == DVEC_TYPE::OCT_SHARED_NODES)
            m_old->interGridTransfer(in,out,m_new,ot::INTERGRID_TRANSFER_MODE::INJECTION,dof);
        else if (dvec.m_vec_type == DVEC_TYPE::OCT_LOCAL_NODES)
            m_old->interGridTransfer_DG(in,out,m_new,dof);    
        else if (dvec.m_vec_type == DVEC_TYPE::OCT_CELL_CENTERED)
            m_old->interGridTransferCellVec(in,out,m_new,dof);
        else
        {
            dendro_log("Invalid vec mode for intergrid transfer");
            MPI_Abort(dvec.get_communicator(),0);
        }

        // printf("%p\n", vec_tmp.get_vec_ptr());
        // printf("%p\n", dvec.get_vec_ptr());
        std::swap(vec_tmp , dvec);
        // printf("%p\n", vec_tmp.get_vec_ptr());
        // printf("%p\n", dvec.get_vec_ptr());
        vec_tmp.destroy_vector();
        
        return;

    }


}// end of namespace ot


