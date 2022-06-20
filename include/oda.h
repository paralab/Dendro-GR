/**
  @file oda.h
  @brief The class that manages the octree mesh that support FEM computations. Note that this file is a refactored version from Dendro4
  which is changed to support Dendro-5.0  Currently we use mesh based FEM computations but in future we will move towards the mesh - free FEM computation methods.     
 
  Distributed Array (DA) methods currently uses ot::Mesh class to provide interface to write FEM loops. 

  @author Milinda Fernando
  date: 10/29/2018
  */ 

#ifndef DENDRO_5_ODA_H
#define DENDRO_5_ODA_H

#include "dendro.h"
#include <iostream>
#include <vector>
#include "mpi.h"
#include "mesh.h"
#include "parUtils.h"
#include "testUtils.h"
#include "asyncExchangeContex.h"
#include "oct2vtk.h"
#include "odaUtils.h"
#include "matRecord.h"
#include <assert.h>

#ifdef BUILD_WITH_PETSC
    #include "petsc.h"
    #include "petscvec.h"
    #include "petscmat.h"
    #include "petscdmda.h"
#endif

#define VECType DendroScalar
#define MATType DendroScalar

/***
 * Note: All the Dendro vectors are stored like as shown below in the memory.
 * |------- pre ghost nodes----- |-----------local nodes---------------|--------post ghost nodes--------|
 *
 *
 *
 * */

// ODA flag information.
namespace ot
{
    namespace DA_FLAGS
    {
        /**
         * @brief loop flags,
         * ALL : Loop over all the elements (note here we loop only on the local + level 1 ghost elements).
         * WRITABLE: Loop over INDEPENDENT U W_DEPENDENT
         * INDEPENDENT : Loop over all the local elements which DOES NOT point to any ghost node region.
         * W_DEPENDENT : Loop over all local elements which has AT LEAST ONE node which point to ghost node region
         * W_BOUNDARY : Loop over all the local elements which is on the domain boundary.
         * 
         * LOCAL_ELEMENTS : Loop over local elements of the mesh
         * PREGHOST_ELEMENTS : Loop over pre ghost elements. 
         * POSTGHOST_ELEMENTS : Loop over post ghost elements
         * 
         *
         * */
        enum LoopType {ALL,WRITABLE,INDEPENDENT,W_DEPENDENT,W_BOUNDARY,LOCAL_ELEMENTS,PREGHOST_ELEMENTS,POSTGHOST_ELEMENTS};

        /**
         * X_MIN : boundary of the x min on the computational domain.
         * X_MIN : boundary of the x min on the computational domain
         *
         * Y_MIN : boundary of the y min on the computational domain
         * Y_MIN : boundary of the y min on the computational domain
         *
         * Z_MIN : boundary of the z min on the computational domain
         * Z_MIN : boundary of the z min on the computational domain
         * */
        enum BdyType {X_MIN,X_MAX,Y_MIN,Y_MAX,Z_MIN,Z_MAX};


        /**
         * @brief contains the refine flags.
         * DA_NO_CHANGE : no change needed for the octant
         * DA_REFINE : refine the octant
         * DA_COARSEN: coarsen the octant.
         * **/
        enum Refine {DA_NO_CHANGE,DA_REFINE,DA_COARSEN};

        enum WriteMode{SET_VALUES=0,ADD_VALUES};



    }

}




namespace ot
{
    class Mesh;
    class TreeNode;

    /**@brief: DA type OCT: uses Dendro DA and PETSC uses petsc based DA*/
    enum DAType {OCT,PETSC};

    /**@brief: denotes the type of the matvec being computed,
     * MESH_BASED - denotes the MatVec computation using mesh data structures.
     * MESH_FREE - denotes the MatVec computation without using mesh data structures.
     * */
    enum MVECType {MESH_BASED,MESH_FREE};


    /**
     * @brief: Simple structure to keep the loop counters
     */
    struct LoopCounter
    {   
        /**@brief: current Index*/
        unsigned int currentIndex;
        unsigned int indexBegin;
        unsigned int indexEnd;
    };

    class DA
    {

    private:
        /**@brief: Loop counter info*/
        LoopCounter m_uiLoopInfo;

        /**@brief: min of octree domain*/
        double m_uiOctreeLowerBound[m_uiDim];

        /**@brief: max of octree domain*/
        double m_uiOctreeUpperBound[m_uiDim];

        /**@brief flag the octree elements to support FEM loops.
         * bit 0 -> ON if the octant is INDEPENDENT
         * bit 1 -> ON if the octant is W_DEPENDENT
         * bit 2 -> ON if the octant is W_BOUNDARY
         * other bits are free not used currently
         */
        std::vector<unsigned int> m_uiOctantFlags;

        /**@brief generated higher order mesh*/
        ot::Mesh* m_uiMesh;


        /**@brief total nodal size (zipped nodal size)*/
        unsigned int m_uiTotalNodalSz;
        /**@brief number of local nodes*/
        unsigned int m_uiLocalNodalSz;

        /**@brief total element size*/
        unsigned int m_uiTotalElementSz;

        /**@brief local elemental size*/
        unsigned int m_uiLocalElementSz;

        /**@brief contexts for async data transfers*/
        std::vector<ot::AsyncExchangeContex> m_uiMPIContexts;

        /**@brief mpi tags*/
        unsigned int m_uiCommTag=0;

        /**@brief local to global map (nodal)*/
        std::vector<DendroIntL> m_uiLocalToGlobalNodalMap;

        /**@brief total number of nodes accross all the processes*/
        DendroIntL m_uiGlobalNodeSz;

        std::vector<DendroIntL> m_uiNodalOffset;




    public:
        /**@brief: returns the vector containing only local elements  */
        std::vector<ot::TreeNode> getLocalOctants();

        /**@brief: get octree to block decomposition blocks*/
        const std::vector<ot::TreeNode>& getBlocks();

        /**@brief: Constructor for the DA data structures
         * @param [in] in : input octree, need to be 2:1 balanced unique sorted octree.
         * @param [in] comm: MPI global communicator for mesh generation.
         * @param [in] order: order of the element.
         * @param [in] grainSz: Number of suggested elements per processor,
         * @param [in] sfc_tol: SFC partitioning tolerance,
         * */
        DA(std::vector<ot::TreeNode> &balOct,MPI_Comm comm,unsigned int order, unsigned int grainSz=100,double sfc_tol=0.3, SM_TYPE smType=SM_TYPE::FEM_CG);

        /**
         * @brief: Create a DA from a specified mesh
         * @param[in] pMesh : input mesh
         * */
        DA(ot::Mesh* pMesh);


        /**
         * @biref Construct a DA from a function
         * @param[in] func: function to capture
         * @param[in] dofSz: size of the degrees of freedoms
         * @param[in] comm: MPI communicator
         * @param[in] order: element order
         * @param[in] interp_tol: interpolation tolerance for func approximation/ capturing
         * @param[in] grainSz: number of elements per core. 
         * @param[in] sfc_tol: flexible partitioning tolerance. 
         * @param[in] SM_TYPE: scatter map type
         * */
        template<typename T>
        DA(std::function<void(T,T,T,T*)>func,unsigned int dofSz,MPI_Comm comm,unsigned int order,double interp_tol, unsigned int grainSz=100,double sfc_tol=0.3, SM_TYPE smType=SM_TYPE::FEM_CG);

        /**
         * @brief Construct a new DA object from specified point set
         * 
         * @param[in] pts : input points
         * @param numPts : number of points
         * @param pt_min : global min of point
         * @param pt_max : global max of point
         * @param comm : MPI communicator
         * @param order : element order
         * @param grainSz : number of elements per core
         * @param sfc_tol : SFC partition tolerance
         * @param smType : scatter map type
         */
        DA(const Point* pts, unsigned int numPts,Point pt_min,Point pt_max, MPI_Comm comm,unsigned int order,unsigned int grainSz=100,double sfc_tol=0.3, SM_TYPE smType=SM_TYPE::FEM_CG);



        /**
         * @brief deconstructor for the DA class.
         * */
        ~DA();

        /**@brief returns the local nodal size*/
        inline unsigned int getLocalNodalSz() const {return m_uiLocalNodalSz;}

        /**@brief returns the pre ghost nodal size*/
        inline unsigned int getPreNodalSz() const {return m_uiMesh->getNumPreMeshNodes();}

        /**@brief returns the post nodal size*/
        inline unsigned int getPostNodalSz() const {return m_uiMesh->getNumPostMeshNodes();}

        /**@brief returns the total nodal size (this includes the ghosted region as well.)*/
        inline unsigned int getTotalNodalSz() const {return m_uiTotalNodalSz;}

        /**@brief returns the local elemental size*/
        inline unsigned int getLocalElemSz() const {return m_uiLocalElementSz;}

        /**@brief returns the local elemental size (includes the ghost elements as well)*/
        inline unsigned int getTotalElemSz() const {return m_uiTotalElementSz;}

        /**@brief see if the current DA is active*/
        inline bool isActive() const {return m_uiMesh->isActive();}

        /**@brief get number of nodes per element*/
        inline unsigned int getNumNodesPerElement() const {return m_uiMesh->getNumNodesPerElement();}

        /**@brief get element order*/
        inline unsigned int getElementOrder() const {return m_uiMesh->getElementOrder();}

        /**@brief: returns the global MPI communicator*/
        inline MPI_Comm getGlobalComm() const { return m_uiMesh->getMPIGlobalCommunicator();}

        /**@brief: returns the get nodal offset values. */
        inline const std::vector<DendroIntL>& getNodalOffsets() const { return m_uiNodalOffset; }
        
        /**@brief: returns node local to node global map*/
        inline const std::vector<DendroIntL>& getNodeLocalToGlobalMap() const {return m_uiLocalToGlobalNodalMap;}

        /**@brief returns the mesh*/
        inline const ot::Mesh* getMesh() const {return m_uiMesh;}

        /**@brief: Returns octDA flags */
        inline const std::vector<unsigned int> getOctFlags() const {return m_uiOctantFlags; }

        /**@brief: returns the active MPI sub com of the global communicator*/
        inline MPI_Comm getCommActive() const
        {
           if(m_uiMesh->isActive())
               return m_uiMesh->getMPICommunicator();
           else
               return MPI_COMM_NULL;
        }

        /**@brief: global mpi com. size*/
        inline unsigned int getNpesAll() const { return  m_uiMesh->getMPICommSizeGlobal();} ;
        
        /**@brief: number of processors active */
        inline unsigned int getNpesActive() const
        {
            if(m_uiMesh->isActive())
                return m_uiMesh->getMPICommSize();
            else
                return 0;
        }

        /**@brief: rank with respect to the global comm. */
        inline unsigned int getRankAll() const {return m_uiMesh->getMPIRankGlobal();};

        /**@brief: rank w.r.t active comm.  */
        inline unsigned int getRankActive() const
        {
            if(m_uiMesh->isActive())
                return m_uiMesh->getMPIRank();
            else
                return m_uiMesh->getMPIRankGlobal();
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

        
        void setOctantWeight(unsigned int ele,unsigned int w) { 
            ot::TreeNode* pNodes = (ot::TreeNode*)m_uiMesh->getAllElements().data();
            pNodes[ele].setWeight(w);
            return;
        }

        unsigned int getOctantWeight(unsigned int ele) const  { return m_uiMesh->getAllElements()[ele].getWeight(); }

        /**@brief get a constant pointer for the reference element*/
        inline const RefElement* getReferenceElement()const { return  m_uiMesh->getReferenceElement();}


        /**@brief: get the number of local elements. */
        inline unsigned int getElementSize() const {return m_uiMesh->getNumLocalMeshElements();}

        /**@brief: get the number of pre-ghost element size*/
        inline unsigned int getPreGhostElementSize() const {return m_uiMesh->getNumPreGhostElements();}

        /**@brief: get the number of post-ghost element size*/
        inline unsigned int getPostGhostElementSize() const {return m_uiMesh->getNumPostGhostElements(); }

        /**@brief: get the number of pre and post ghost elements*/
        inline unsigned int getPreAndPostGhostNodeSize() const {return (m_uiMesh->getNumPreMeshNodes() + m_uiMesh->getNumPostMeshNodes());}

        /**@brief number of layer 1 ghost elements. */
        inline unsigned int getNumLayer1GhostEleSz() const {return m_uiMesh->getLevel1GhostElementIndices().size();}

        /**@brief: get the number of pre and post ghost elements*/
        unsigned int getGhostedElementSize() const {return (m_uiMesh->getNumPreGhostElements() + m_uiMesh->getNumPostGhostElements());}


        /**@brief: get the max depth of the octree*/
        inline unsigned int getMaxDepth() const {return  m_uiMaxDepth;};
        /**@brief: get the dimensionality of the octree*/
        inline unsigned int getDimension() const {return m_uiDim;};

        /** @brief get min local node (nodal)*/
        inline ot::TreeNode getMinTreeNode() const { 

            if(m_uiMesh->isActive())
            {
                const ot::TreeNode* allElements = &(*(m_uiMesh->getAllElements().begin()));
                return allElements[m_uiMesh->getElementLocalBegin()];
            }else
            {
                return ot::TreeNode(0,0,0,0,m_uiDim,m_uiMaxDepth);                
            }
            
        }
        
         /** @brief get max localnode (nodal)*/
        inline ot::TreeNode getMaxTreeNode() const { 
            
            if(m_uiMesh->isActive())
            {
                const ot::TreeNode* allElements = &(*(m_uiMesh->getAllElements().begin()));
                return allElements[m_uiMesh->getElementLocalEnd()];
            }else
            {
                return ot::TreeNode(0,0,0,0,m_uiDim,m_uiMaxDepth);                
            }

            
        }

        /**
         * @brief 
         * @param [in] pNodes: Nodes that needs to be searched across. 
         * @param[in] n: number of nodes. 
         * @param[in/out] ownerrank for i^th node, -1 if not found. (Needs to be allocated, outsize for size n);
          */
        void computeTreeNodeOwnerProc(const ot::TreeNode * pNodes, unsigned int n, int* ownerranks);

        
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

         /**@brief Sync accumilation across ghost elements
          * @param [in] vec: vector pointer
          * @param [in] mode: mode of the write to ghost
          * @param [in] dof: degrees of freedoms. 
         */
         template<typename T>
         void writeToGhostsEnd(T* vec, DA_FLAGS::WriteMode mode,unsigned int dof=1) ;


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
        inline unsigned int getLevel(unsigned int ele) const { return m_uiMesh->getAllElements()[ele].getLevel();}

        /**
         * @brief return the TreeNode of the current octnat.
         * @param[in] ele: elementID of the current octant
        * */
        inline ot::TreeNode getOctant(unsigned int ele) const {return m_uiMesh->getAllElements()[ele];}

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
        ot::DA* remesh(const DA_FLAGS::Refine * flags, unsigned int sz,unsigned int grainSz=100,double ld_bal=0.3, unsigned int sfK=2, unsigned int (*getWeight)(const ot::TreeNode *)=NULL) const;

        /**
         * @brief 
         * 
         * @param[in] grainSz: rougly the number of octants per core you need when you create the new da.
         * @param[in] ld_tol: load imbalance tolerance.
         * @param[in] sfK: splitter fix factor. better to be power of two. increase the value to 128 when running on > 64,000 cores
         * @param getWeight weight function 
         * @return ot::DA* 
         */
        ot::DA* repartition(unsigned int grainSz=100,double ld_bal=0.3, unsigned int sfK=2, unsigned int (*getWeight)(const ot::TreeNode *)=NULL) const;

        /**
         * @brief performs grid transfer operations after the remesh.
         * @param[in] varIn: variable defined by oldDA
         * @param[out] varOut: variable defined by newDA. interpolate varOut from varIn. (Note: varOut allocated inside the function, no need to allocate outside)
         * @param[in] isElemental: true if it is an elemental vector
         * @param[in] isGhosted: true if allocated ghost vector
         * @param[in] dof: degrees of freedoms.
         * */
        template<typename T>
        void intergridTransfer(const T* varIn, T* & varOut, const ot::DA* newDA, bool isElemental=false, bool isGhosted=false, unsigned int dof=1,INTERGRID_TRANSFER_MODE mode = INTERGRID_TRANSFER_MODE::INJECTION);



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
        int getFaceNeighborValues(unsigned int eleID, const T* in, T* out, T* coords, unsigned int * neighID, unsigned int face, NeighbourLevel & level,unsigned int dof) const;

        /**
        * @brief computes the child number of the given octant
        * @param [in] eleID: element ID
        * returns the Morton child number of the eleID in the octant
        * */
        inline unsigned int getMortonChildNum(unsigned int eleID) const
        {
           return m_uiMesh->getMortonchildNum(eleID);
        }

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
        void petscIntergridTransfer(const Vec &varIn, Vec &varOut, const ot::DA *newDA, bool isElemental = false,bool isGhosted = false, unsigned int dof = 1,INTERGRID_TRANSFER_MODE mode = INTERGRID_TRANSFER_MODE::INJECTION);

        /**
         * @brief: Destroy petsc vector
         * @param[in] vec: petsc vector
         */
        PetscErrorCode petscDestroyVec(Vec & vec);
      

#endif



    };


}


#include "oda.tcc"


#endif

