/***
 * @file        sub_oda.h
 * @author      Hari Sundar   hsundar@gmail.com
 * @author      Milinda Fernando milinda@cs.utah.edu
 * @date        13 June 2019
 * @brief       Subdomain for Octree DA (ODA).
 ***/

#pragma once

#include <functional>
#include <iostream>

#include "oda.h"

namespace ot {
  
  /**
   * @brief subDA based on the DA, where to perform computations sub domain of the main domain. 
   * 
   */
  class subDA {

      protected:

        /**@brief: parent da */
        ot::DA*       m_da;
        /** @brief:elemental skip list*/
        std::vector<unsigned char>  m_ucpSkipList;
        /** @brief:nodal skip list*/
        std::vector<unsigned char>  m_ucpSkipNodeList;
        /** @brief: local element size */
        unsigned int  m_uiLocalElementSize;
        /** @brief: pre ghost element size */
        unsigned int  m_uiPreGhostElementSize;
        /** @brief: post ghost element size */
        unsigned int  m_uiPostGhostElementSize;
        /** @brief: Local node size */
        unsigned int  m_uiLocalNodeSize;
        /** @brief: pre ghost node size */
        unsigned int  m_uiPreGhostNodeSize;
        /** @brief: post ghost node size */
        unsigned int  m_uiPostGhostNodeSize;
        /** @brief:  async comm tag*/
        unsigned int m_uiCommTag;

        std::vector<AsyncExchangeContex> m_uiMPIContexts;

        /**@brief: Loop counter info*/
        LoopCounter m_uiLoopInfo;
        
        /** @brief:  async comm contex*/
        std::vector<AsyncExchangeContex>  m_mpiContexts;
        /** @brief: subDA to DA map  (elemental)*/
        std::vector<unsigned int> m_uip_sub2DA_ElemMap;      
        /**@brief: DA to subDA map (elemental)*/
        std::vector<unsigned int> m_uip_DA2sub_ElemMap;      
        /** @brief: subDA to DA map  (nodal)*/
        std::vector<unsigned int> m_uip_sub2DA_NodeMap;      
        /** @brief: DA to subDA map  (nodal)*/
        std::vector<unsigned int> m_uip_DA2sub_NodeMap;

        /** @brief min bound box for fe_retain */
        double   m_dMinBB[3];
        /** @brief max bound box for fe_retain */
        double   m_dMaxBB[3];

        /**@brief total nodal size (zipped nodal size)*/
        unsigned int m_uiTotalNodalSz;
        /**@brief number of local nodes*/
        unsigned int m_uiLocalNodalSz;

        /**@brief total element size*/
        unsigned int m_uiTotalElementSz;

        /**@brief local elemental size*/
        unsigned int m_uiLocalElementSz;

        /**begin of the pre ghost elements*/
        unsigned int m_uiElementPreGhostBegin = 0; // Not mandatory to store this.
        /** end of pre ghost elements*/
        unsigned int m_uiElementPreGhostEnd;
        /** begin of locat octants (INDEPENDENT)*/
        unsigned int m_uiElementLocalBegin;
        /** end of the local octants (INDEPENDET)*/
        unsigned int m_uiElementLocalEnd;
        /**begin location for the post ghost octants*/
        unsigned int m_uiElementPostGhostBegin;
        /** end location for the post ghost octants*/
        unsigned int m_uiElementPostGhostEnd;

        /** begin location  of the pre ghost nodes in CG indexing*/
        unsigned int m_uiNodePreGhostBegin=0;
        /** end location of the pre ghost nodes in CG indexing*/
        unsigned int m_uiNodePreGhostEnd;
        /** begin location of the local nodes in CG indexing*/
        unsigned int m_uiNodeLocalBegin;
        /** end location of the local nodes in CG indexing*/
        unsigned int m_uiNodeLocalEnd;
        /** begin location of the post ghost nodes in CG indexing*/
        unsigned int m_uiNodePostGhostBegin;
        /** end location of the post ghost nodes in CG indexing*/
        unsigned int m_uiNodePostGhostEnd;

        /** number of pre ghost elements */
        unsigned int m_uiNumPreGhostElements;
        /** number of local elements */
        unsigned int m_uiNumLocalElements;
        /** number of post ghost elements */
        unsigned int m_uiNumPostGhostElements;

        /** number of pre ghost nodes */
        unsigned int m_uiNumPreGhostNodes;
        /** number of local nodes */
        unsigned int m_uiNumLocalNodes;
        /** number of post nodes */
        unsigned int m_uiNumPostGhostNodes;

        /** @brief: Send scatter map */
        std::vector<unsigned int > m_uiSendScatterMap;
        /** @brief: recv scatter map */
        std::vector<unsigned int > m_uiRecvScatterMap;
        
        /**@brief: send counts */
        std::vector<unsigned int> m_uiSendCounts;
        /**@brief: recv counts */
        std::vector<unsigned int> m_uiRecvCounts;
        /**@brief: send offsets */
        std::vector<unsigned int> m_uiSendOffsets;
        /**@brief: recv offsets */
        std::vector<unsigned int> m_uiRecvOffsets;
        /**Send processor list for ghost exchange*/
        std::vector<unsigned int> m_uiSendProcList;
        /**recv processor list for the ghost exchange*/
        std::vector<unsigned int> m_uiRecvProcList;

        /**@brief: active sub da */
        bool m_uiIsActive;

        /**@brief: MPI comm global */
        MPI_Comm m_uiCommGlobal;

        /**@brief: MPI_Comm for the active*/
        MPI_Comm m_uiCommActive;

        /**@brief: rank active */
        int m_uiRankActive;
        
        /**@brief: rank global */
        int m_uiRankGlobal;
        
        /**@brief: npes active */
        int m_uiNpesActive;

        /**@brief: npes global */
        int m_uiNpesGlobal;

       
      public:
        /**
         * @brief: defult sub da constructor. 
         * @param[in] da: input DA
         * @param[in] fx_retain : area to retain
         * param[in] gSize: global domain size. 
         */
        subDA(DA* da, std::function<double ( double, double, double ) > fx_retain, double* gSize);

        ~subDA();
        
        ot::DA* global_domain() { return m_da; }   

        /**@brief returns the local nodal size*/
        inline unsigned int getLocalNodalSz() const {return m_uiLocalNodalSz;}

        /**@brief returns the pre ghost nodal size*/
        inline unsigned int getPreNodalSz() const { return m_uiNumPreGhostNodes; }

        /**@brief returns the post nodal size*/
        inline unsigned int getPostNodalSz() const {return m_uiNumPostGhostNodes;}

        /**@brief returns the total nodal size (this includes the ghosted region as well.)*/
        inline unsigned int getTotalNodalSz() const {return m_uiTotalNodalSz;}

        /**@brief returns the local elemental size*/
        inline unsigned int getLocalElemSz() const {return m_uiLocalElementSz;}

        /**@brief returns the local elemental size (includes the ghost elements as well)*/
        inline unsigned int getTotalElemSz() const {return m_uiTotalElementSz;}

        /**@brief see if the current DA is active*/
        inline bool isActive() const {return m_da->isActive();}

        /**@brief get number of nodes per element*/
        inline unsigned int getNumNodesPerElement() const {return m_da->getNumNodesPerElement();}

        /**@brief get element order*/
        inline unsigned int getElementOrder() const {return m_da->getElementOrder();}

        /**@brief: returns the global MPI communicator*/
        inline MPI_Comm getGlobalComm() const { return m_da->getGlobalComm();}
        
        /**@brief: returns node local to node global map*/
        //inline const std::vector<DendroIntL> getNodeLocalToGlobalMap() const {return m_uiLocalToGlobalNodalMap;}

        /**@brief returns the mesh*/
        inline const ot::Mesh* getMesh() const {return m_da->getMesh();}

        /**@brief: returns the active MPI sub com of the global communicator*/
        inline MPI_Comm getCommActive() const
        {
           if(m_da->isActive())
               return m_da->getCommActive();
           else
               return MPI_COMM_NULL;
        }

        /**@brief: global mpi com. size*/
        inline unsigned int getNpesAll() const { return  m_da->getNpesAll();} ;
        
        /**@brief: number of processors active */
        inline unsigned int getNpesActive() const
        {
            if(m_da->isActive())
                return m_da->getNpesActive();
            else
                return 0;
        }

        /**@brief: rank with respect to the global comm. */
        inline unsigned int getRankAll() const {return m_da->getRankAll();};

        /**@brief: rank w.r.t active comm.  */
        inline unsigned int getRankActive() const
        {
            if(m_da->isActive())
                return m_da->getRankActive();
            else
                return m_da->getRankAll();
        }

        /**@brief: returns true if specified eleID is a boundary element
         * false if the eleID is local and not a boundary element.
         * @param[in] eleID: element ID
         * */
        bool isBoundaryOctant(unsigned int eleID) const;

        /**@brief computes the elementCoordinates (based on the nodal placement)
        * @param[in] eleID : element ID
        * @param[in/out] coords: computed coords (note: assumes memory is allocated allocated)
        * coords are stored by p0,p1,p2... each pi \in R^dim where pi are ordered in along x axis then y axis and  then z axis
        * coors size m_uiDim*m_uiNpE
         * for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                {
                    coords[(k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i)*m_uiDim+0]=x+i*dx;
                    coords[(k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i)*m_uiDim+1]=y+j*dy;
                    coords[(k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i)*m_uiDim+2]=z+k*dz;
                }
         *
         *
        * */
        void getElementalCoords(unsigned int eleID, double* coords) const ;


        /**@brief get a constant pointer for the reference element*/
        inline const RefElement* getReferenceElement()const { return  m_da->getReferenceElement();}


        /**@brief: get the number of local elements. */
        inline unsigned int getElementSize() const {return m_uiNumLocalElements;}

        /**@brief: get the number of pre-ghost element size*/
        inline unsigned int getPreGhostElementSize() const {return m_uiNumPreGhostElements;}

        /**@brief: get the number of post-ghost element size*/
        inline unsigned int getPostGhostElementSize() const {return m_uiNumPostGhostElements; }

        /**@brief: get the number of pre and post ghost elements*/
        inline unsigned int getPreAndPostGhostNodeSize() const {return (m_uiNumPreGhostNodes + m_uiNumPostGhostNodes);}
        
        /**@brief: get the number of pre and post ghost elements*/
        unsigned int getGhostedElementSize() const {return (m_uiNumPreGhostElements + m_uiNumPostGhostElements);}

        /**@brief: get the max depth of the octree*/
        inline unsigned int getMaxDepth() const {return  m_uiMaxDepth;};
        /**@brief: get the dimensionality of the octree*/
        inline unsigned int getDimension() const {return m_uiDim;};

        /**
         * @brief Creates a ODA vector
         * @param [in] local : VecType pointer
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         * */
         template<typename T>
         int createVector(T*& local, bool isElemental=false, bool isGhosted=false, unsigned int dof=1) const;

        /**
        * @brief Creates a ODA vector std::vector<T>
        * @param [in] local : VecType pointer
        * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
        * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
        * @param [in] dof: degrees of freedoms
        * */
        template<typename T>
        int createVector(std::vector<T>& local, bool isElemental=false, bool isGhosted=false, unsigned int dof=1) const;

        /**
         * @brief deallocates the memory allocated for a vector
         * @param[in/out] local: pointer to the vector
         * */
        template <typename T>
        void destroyVector(T*& local) const;

        template <typename T>
        void destroyVector(std::vector<T>& local) const;



         /**@brief initialize the loop counters. */
         template<ot::DA_FLAGS::LoopType type>
         void init();

         /**@brief get the current elmenent in the iteration*/
         unsigned int curr();


         /**@brief Returns if the current element has reached end or not. */
         template<ot::DA_FLAGS::LoopType type>
         unsigned int end();

         /**@brief: increment it to the next element. */
         template<ot::DA_FLAGS::LoopType type>
         void next();

         /**
         * @brief Initiate the ghost nodal value exchange
         * */
         template<typename T>
         void readFromGhostBegin(T* vec, unsigned int dof=1);

         /**
          * @brief Sync the ghost element exchange
          * */
         template<typename T>
         void readFromGhostEnd(T* vec, unsigned int dof=1);


         /**@brief Initiate accumilation across ghost elements*/
         template<typename T>
         void writeToGhostsBegin(T* vec, unsigned int dof=1) ;

         /**@brief Sync accumilation across ghost elements*/
         template<typename T>
         void writeToGhostsEnd(T* vec, unsigned int dof=1) ;


         /**
          * @brief convert nodal local vector with ghosted buffer regions.
          * @param[in] in: input vector (should be nodal and non ghosted)
          * @param[out] out: coverted nodal vector with ghost regions.
          * @param[in] isAllocated: true if the out is allocated, false otherwise.
          * @param[in] dof: degrees of freedoms
          * */
         template<typename T>
         void nodalVecToGhostedNodal(const T* in,T*& out,bool isAllocated=false, unsigned int dof=1)const;


        /**
         * @brief convert ghosted nodal vector to local vector (without ghosting)
         * @param[in] gVec: ghosted vector
         * @param[out] local: local vector (assume an allocated vector)
         * @param[in] isAllocated: true if the out is allocated, false otherwise.
         * @param[in] dof: degrees of freedoms
         * */

         template<typename T>
         void ghostedNodalToNodalVec(const T* gVec,T*& local,bool isAllocated=false, unsigned int dof=1) const;



        /**
         * @brief computes the element nodal values using interpolation if needed.
         * @param[in] in: input vector (need to be ghosted )
         * @param[in/out] eleVec: allocated eleVec (size: numNodesPerElement*dof) after the function call will have the interpolated values.
         * @param[in] eleID: element ID
         * @param[in] dof: degrees of freedoms
         * */
         template <typename T>
         void getElementNodalValues(const T*in, T* eleVecOut,unsigned int eleID,unsigned int dof=1)const;

        /**
          * @brief computes the elemental vec to global vec accumilation
          * @param[out] out: output (need to be ghosted )
          * @param[in] eleVecIn: eleVec values.
          * @param[in] eleID: element ID
          * @param[in] dof: degrees of freedoms
          * */
        template <typename T>
        void eleVecToVecAccumilation(T*out, const T* eleVecIn,unsigned int eleID,unsigned int dof=1)const;

        /**
         * @brief Computes the octree writable boundary nodes.
         * @param[in/out]: bdyIndex: boundary index nodal locations (only writable boundary locations)
         * @param[in] isGhosted: change the indices to the ghosted array.
         * */
        void getOctreeBoundaryNodeIndices(std::vector<unsigned int >& bdyIndex,std::vector<double>& coords,bool isGhosted=false);

        /**
         * @brief Returns the local indices to the nodes of the current element.
         * @param nodes   Indices into the nodes of the given element. Should be allocated by the user prior to calling.
           @return -1 for error
        */
        int getNodeIndices(DendroIntL* nodeIdx,unsigned int ele,bool isGhosted) const;
        
        /**
         * @brief Returns the global node indices of a given element,
         * @param[out] nodeIdx: global node indices, 
         * @param[in] ele: element id
         */
        int getGlobalNodeIndices(DendroIntL* nodeIdx,unsigned int ele) const ;

        /**
         * @brief initialize a variable vector to a function depends on spatial coords.
         * @param[in/out] local: allocated vector, initialized vector
         * @param[in] func: user specified function
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         *
         * */
        template <typename T>
        void setVectorByFunction(T* local,std::function<void(T,T,T,T*)>func,bool isElemental=false, bool isGhosted=false, unsigned int dof=1) const;


        /**
         * @brief initialize a variable vector to a function depends on spatial coords.
         * @param[in/out] local: allocated vector, initialized vector
         * @param[in] value: user specified scalar values (size should be the  dof size)
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         * Note: Initialize the ghost region as well.
         *
         * */
        template <typename T>
        void setVectorByScalar(T* local,const T* value,bool isElemental=false, bool isGhosted=false, unsigned int dof=1) const;


        /**@brief write the vec to pvtu file
         * @param[in] local: variable vector
         * @param[in] fPrefix: file name prefix
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         * */
        template <typename T>
        void vecTopvtu(T* local, const char * fPrefix,char** nodalVarNames=NULL,bool isElemental=false,bool isGhosted=false,unsigned int dof=1);


        /**@brief return the level of current octant
         * @param[in] ele: elementID of the current octant
         * */
        inline unsigned int getLevel(unsigned int ele) const { return m_da->getLevel(m_uip_sub2DA_ElemMap[ele]);}

        /**
         * @brief return the TreeNode of the current octnat.
         * @param[in] ele: elementID of the current octant
        * */
        inline ot::TreeNode getOctant(unsigned int ele) const {return m_da->getOctant(m_uip_sub2DA_ElemMap[ele]);}

        /**
         * @brief returns a pointer to a dof index,
         * @param [in] in: input vector pointer
         * @param [in] dofInex: dof index which is the pointer is needed, should be less than dof, value the vector created.
         * @param [in] isElemental: true if this is an elemental vector/ false otherwise
         * @param [in] isGhosted: true if this is a ghosted vector
         * @return pointer to dofIndex.
         * */
        template<typename T>
        T* getVecPointerToDof(T* in ,unsigned int dofInex, bool isElemental=false,bool isGhosted=false) const;


        /**
         * @brief copy vecotor to sorce to destination, assumes the same number of dof.
         * @param [in/out] dest: destination pointer
         * @param [in] source: source pointer
         * @param [in] isElemental: true if this is an elemental vector/ false otherwise
         * @param [in] isGhosted: true if this is a ghosted vector
         * @param [in] dof: degrees of freedoms
         * */
        template<typename T>
        void copyVectors(T* dest,const T* source,bool isElemental=false,bool isGhosted=false,unsigned int dof=1) const;

        /**
         * @brief more premitive copy, from source pointer to the dest pointer
         * @param [in/out] dest: destination pointer
         * @param [in] source: source pointer
         * @param [in] isElemental: true if this is an elemental vector/ false otherwise
         * @param [in] isGhosted: true if this is a ghosted vector
         * */
        template<typename T>
        void copyVector(T* dest,const T* source,bool isElemental=false,bool isGhosted=false) const;


        /**
         * @brief: Performs remesh based on the DA_FLAGS::Refine, which specifies no change, refine or coarsen.
         * @param[in] flags: refinement flags.
         * @param[in] sz: size of the array flags (needs to be size of the local elements)
         * @param[in] grainSz: rougly the number of octants per core you need when you create the new da.
         * @param[in] ld_tol: load imbalance tolerance.
         * @param[in] sfK: splitter fix factor. better to be power of two. increase the value to 128 when running on > 64,000 cores
         * @return: Specifies the new grid, with new DA.
         */
        ot::DA* remesh(const DA_FLAGS::Refine * flags, unsigned int sz,unsigned int grainSz=100,double ld_bal=0.3, unsigned int sfK=2) const;

        /**
         * @brief performs grid transfer operations after the remesh.
         * @param[in] varIn: variable defined by oldDA
         * @param[out] varOut: variable defined by newDA. interpolate varOut from varIn. (Note: varOut allocated inside the function, no need to allocate outside)
         * @param[in] isElemental: true if it is an elemental vector
         * @param[in] isGhosted: true if allocated ghost vector
         * @param[in] dof: degrees of freedoms.
         * */
        template<typename T>
        void intergridTransfer(const T* varIn, T* & varOut, const ot::DA* newDA, bool isElemental=false, bool isGhosted=false, unsigned int dof=1);



        /**
         * @brief computes the face neighbor points for additional computations for a specified direction.
         * @param [in] eleID: element ID
         * @param [in] in: inpute vector
         * @param [out] out: output vector values are in the order of the x,y,z size : 4*NodesPerElement
         * @param [out] coords: get the corresponding coordinates size: 4*NodesPerElement*m_uiDim;
         * @param [out] neighID: face neighbor octant IDs,
         * @param [in] face: face direction in {OCT_DIR_LEFT,OCT_IDR_RIGHT,OCT_DIR_DOWN, OCT_DIR_UP,OCT_DIR_BACK,OCT_DIR_FRONT}
         * @param [out] level: the level of the neighbour octant with respect to the current octant.
         * returns  the number of face neighbours 1/4 for 3D.
         * */
        template<typename T>
        int getFaceNeighborValues(unsigned int eleID, const T* in, T* out, T* coords, unsigned int * neighID, unsigned int face, NeighbourLevel & level) const;

        
         // all the petsc functionalities goes below with the pre-processor gards.
#ifdef BUILD_WITH_PETSC

        /**
         * @brief Creates a PETSC vector
         * @param [in] local : petsc vector
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         * */

        PetscErrorCode petscCreateVector(Vec &local, bool isElemental, bool isGhosted, unsigned int dof) const;


        /**
         @brief Returns a PETSc Matrix of appropriate size of the requested type.
         @param M the matrix
         @param mtype the type of matrix
         @param dof the number of degrees of freedom per node.
         */
        PetscErrorCode createMatrix(Mat &M, MatType mtype, unsigned int dof=1) const;



        /**
         * @brief convert nodal local vector with ghosted buffer regions.
         * @param[in] in: input vector (should be nodal and non ghosted)
         * @param[out] out: coverted nodal vector with ghost regions.
         * @param[in] isAllocated: true if the out is allocated, false otherwise.
         * @param[in] dof: degrees of freedoms
         * */
        PetscErrorCode petscNodalVecToGhostedNodal(const Vec& in,Vec& out,bool isAllocated=false,unsigned int dof=1) const;


        /**
        * @brief convert ghosted nodal vector to local vector (without ghosting)
        * @param[in] gVec: ghosted vector
        * @param[out] local: local vector (assume an allocated vector)
        * @param[in] isAllocated: true if the out is allocated, false otherwise.
        * @param[in] dof: degrees of freedoms
        * */

        PetscErrorCode petscGhostedNodalToNodalVec(const Vec& gVec,Vec& local,bool isAllocated=false,unsigned int dof=1) const;

        /**
         * @brief Initiate the ghost nodal value exchange
         * @param[in] vec: vector in need to perform ghost exchange (Need be ghosted vector)
         * @param[in] vecArry: pointer to from the VecGetArray()
         * @param[in] dof: Degrees of freedoms
         * */

        void petscReadFromGhostBegin(PetscScalar * vecArry, unsigned int dof=1) ;

        /**
         * @brief Sync the ghost element exchange
         * @param[in] vec: vector in need to perform ghost exchange (Need be ghosted vector)
         * @param[in] vecArry: pointer to from the VecGetArray()
         * @param[in] dof: Degrees of freedoms
         * */
        void petscReadFromGhostEnd(PetscScalar * vecArry, unsigned int dof=1) ;


        /**
         * @brief initialize a variable vector to a function depends on spatial coords.
         * @param[in/out] local: allocated vector, initialized vector
         * @param[in] func: user specified function
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         *
         * */
        template<typename T>
        void petscSetVectorByFunction(Vec& local,std::function<void(T,T,T,T*)>func,bool isElemental=false, bool isGhosted=false, unsigned int dof=1) const;


        /**
         * @brief initialize a variable vector to a function depends on spatial coords.
         * @param[in/out] local: allocated vector, initialized vector
         * @param[in] value: user specified scalar values (size should be the  dof size)
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         * Note: Initialize the ghost region as well.
         *
         * */
        template <typename T>
        void petscSetVectorByScalar(Vec& local,const T* value,bool isElemental=false, bool isGhosted=false, unsigned int dof=1) const;


        /**@brief write the vec to pvtu file
         * @param[in] local: variable vector
         * @param[in] fPrefix: file name prefix
         * @param [in] isElemental: True if creating a elemental vector (cell data vector) false for a nodal vector
         * @param [in] isGhosted: True will allocate ghost nodal values as well, false will only allocate memory for local nodes.
         * @param [in] dof: degrees of freedoms
         * */
        void petscVecTopvtu(const Vec& local, const char * fPrefix,char** nodalVarNames=NULL,bool isElemental=false,bool isGhosted=false,unsigned int dof=1) ;


        /**
         * @brief a wrapper for setting values into the Matrix.  This internally calls PETSc's MatSetValues() function.
         * Call PETSc's MatAssembly routines to assemble the matrix after setting the values. It would be more efficient to set values in chunks by
           calling this function multiple times with different sets of values instead of a single call at the end of the loop. One can use the size of 'records' to determine the number of
           such chunks. Calls to this function with the INSERT_VALUES and ADD_VALUES options cannot be mixed without intervening calls to PETSc's MatAssembly routines.
         * @param mat The matrix
         * @param records The values and their indices
         * @param dof the number of degrees of freedom per node
         * @param mode Either INSERT_VALUES or ADD_VALUES
         * @return an error flag
         * @note records will be cleared inside the function
        */
        PetscErrorCode petscSetValuesInMatrix(Mat mat, std::vector<ot::MatRecord>& records,unsigned int dof, InsertMode mode) const;

        
        
        /** 
         * @brief: Dendro 5 vectors with multiple dof values are stored in 00000000,11111111 for multiple dof values, 
         * if the user wants to perform matrix based computations (i.e. not use matrix free option), this will reorder the 
         * vector valules 01,01,01,...  
         * @param [in,out] v1: input and out put vector, 
         * @param [in] dof: number of dof.  
         */
        PetscErrorCode petscChangeVecToMatBased(Vec& v1,bool isElemental,bool isGhosted,unsigned int dof=1) const;
        
        
        /** 
         * @brief: This is the inverse function for the "petscChangeVecToMatBased"
         * @param [in,out] v1: input and out put vector, 
         * @param [in] dof: number of dof.  
         */
        PetscErrorCode petscChangeVecToMatFree(Vec& v1,bool isElemental,bool isGhosted,unsigned int dof=1) const;

        /**
         * @brief performs grid transfer operations after the remesh with Petsc Vectors.
         * @param[in] varIn: variable defined by oldDA
         * @param[out] varOut: variable defined by newDA. interpolate varOut from varIn. (Note: varOut allocated inside the function, no need to allocate outside)
         * @param[in] isElemental: true if it is an elemental vector
         * @param[in] isGhosted: true if allocated ghost vector
         * @param[in] dof: degrees of freedoms.
        */
        template<typename T>
        void petscIntergridTransfer(const Vec &varIn, Vec &varOut, const ot::DA *newDA, bool isElemental = false,bool isGhosted = false, unsigned int dof = 1);

        /**
         * @brief: Destroy petsc vector
         * @param[in] vec: petsc vector
         */
        PetscErrorCode petscDestroyVec(Vec & vec);
      

#endif



    };


}

#include "sub_oda.tcc"



