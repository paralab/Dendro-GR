//
// Created by milinda on 10/30/18.
//
/**
 * @brief class that derived from abstract class feMat
 * LHS computation of the weak formulation
 * */
#ifndef DENDRO_5_0_FEMATRIX_H
#define DENDRO_5_0_FEMATRIX_H

#include "feMat.h"

template<typename T>
void fem_eMatTogMat(ot::MatRecord *eMat, const T *Imat_p2c, unsigned int pOrder,unsigned int dof=1) {
    const unsigned int n = (pOrder + 1);
    const unsigned int nPe = n * n * n;
    const unsigned int N = nPe * dof;
    const unsigned int dof2=dof*dof;

    // todo : need to fix this properly, with efficient matrix matrix multiplication. 
    T * IT_Ke_I  =new T[N * N];
    
    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj=0; dj< dof; dj++)
    {
        const unsigned int offset = (di*dof + dj)*nPe*nPe;
        for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                T val = (T)0;
                for(unsigned int k=0;k< nPe; k++)
                {
                    val  += (eMat[offset + i*nPe +k].getMatVal() * Imat_p2c[k * nPe + j]);
                }
                IT_Ke_I[offset + i*nPe +j] = val;
            }
    }
    
    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj=0; dj< dof; dj++)
    {
        const unsigned int offset = (di*dof + dj)*nPe*nPe;

         for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                eMat[offset + i*nPe +j].setMatValue(IT_Ke_I[offset + i*nPe +j]);
            }
    }
    

    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj =0; dj< dof; dj++)
    {
         const unsigned int offset = (di*dof + dj)*nPe*nPe;
         for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                T val = (T)0;
                for(unsigned int k=0;k< nPe; k++)
                {
                    val  += ( Imat_p2c[k * nPe + i] * eMat[offset + k*nPe +j].getMatVal());
                }
                IT_Ke_I[offset + i*nPe + j] = val;
            }
    }
    
    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj=0; dj< dof; dj++)
    {
         const unsigned int offset = (di*dof + dj)*nPe*nPe;
         for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                eMat[offset + i*nPe + j].setMatValue(IT_Ke_I[offset + i*nPe + j]);
            }
    }


    delete [] IT_Ke_I; 
}


template<typename T>
class feMatrix : public feMat {

protected:

    /**@brief number of dof*/
    unsigned int m_uiDof;

    /**@brief element nodal vec in */
    VECType *m_uiEleVecIn;

    /***@brief element nodal vecOut */
    VECType *m_uiEleVecOut;

    /** elemental coordinates */
    double *m_uiEleCoords;

public:
    /**
     * @brief constructs an FEM stiffness matrix class.
     * @param[in] da: octree DA
     * */
    feMatrix(ot::DA *da, unsigned int dof = 1);

    ~feMatrix();

    /**@brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
      * @param [in] in input vector u
      * @param [out] out output vector Ku
      * @param [in] default parameter scale vector by scale*Ku
    * */
    virtual void matVec(const VECType *in, VECType *out, double scale = 1.0);

    /**@brief: Computes the diagonal of the matrix for mat free.   
     * @param[out] diag : computed diagonal (allocation happens inside the function)
     * for multiple dof values the diagonal is ordered by the dof. 
     * dof_0 dof_0,,,,,,,dof_1,dof_1,,,,,,dof_2,dof_2 etc (Dendro Vector format. )
     * @paramp[in] scale: scale for the diagonal element.  
    */
    virtual void getMatDiagonal(VECType*& diag, double scale = 1.0);

    /**
     * @brief Get the block jacobi preconditioner
     * 
     * @param blkDiag : list of mat records containing the block diagonal Jacobi preconditioner. 
     * @param scale : scale value (default is 1.0). 
     */
    virtual void getMatBlockDiagonal(std::vector<ot::MatRecord>& blkDiag, double scale = 1.0);


    /**@brief Computes the elemental matvec
      * @param [in] in input vector u
      * @param [out] out output vector Ku
      * @param [in] default parameter scale vector by scale*Ku
    **/
    virtual void elementalMatVec(const VECType *in, VECType *out, double *coords = NULL, double scale = 1.0) = 0;


#ifdef BUILD_WITH_PETSC

    /**@brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
      * @param [in] in input vector u
      * @param [out] out output vector Ku
      * @param [in] default parameter scale vector by scale*Ku
    * */
    virtual void matVec(const Vec &in, Vec &out, double scale = 1.0);

    /**
     * @brief Performs the matrix assembly.
     * @param [in/out] J: Matrix assembled
     * @param [in] mtype: Matrix type
     * when the function returns, J is set to assembled matrix
     **/
    virtual bool getAssembledMatrix(Mat *J, MatType mtype);

    /**
     * @brief Get the Mat Block Diagonal object
     * @param m : block diagonal matrix. 
     */
    virtual PetscInt getMatBlockDiagonal(Mat m, PetscScalar scale=1.0);

    /**
     * @brief Get the Mat Diagonal object
     * 
     * @param vec diagonal of the matrix. 
     */
    virtual PetscInt getMatDiagonal(Vec& vec, PetscScalar scale=1.0);


#endif

    /**@brief static cast to the leaf node of the inheritance*/
    T &asLeaf() { return static_cast<T &>(*this); }


    /**
     * @brief executed just before  the matVec loop in matvec function
     * @param[in] in : input Vector
     * @param[out] out: output vector
     * @param[in] scale: scalaing factror
     **/

    bool preMatVec(const VECType *in, VECType *out, double scale = 1.0) {
        return asLeaf().preMatVec(in, out, scale);
    }


    /**@brief executed just after the matVec loop in matvec function
     * @param[in] in : input Vector
     * @param[out] out: output vector
     * @param[in] scale: scalaing factror
     * */

    bool postMatVec(const VECType *in, VECType *out, double scale = 1.0) {
        return asLeaf().postMatVec(in, out, scale);
    }

    /**@brief executed before the matrix assembly */
    bool preMat() {
        return asLeaf().preMat();
    }

    /**@brief executed after the matrix assembly */
    bool postMat() {
        return asLeaf().postMat();
    }

    /**
     * @brief Compute the elemental Matrix.
     * @param[in] eleID: element ID
     * @param[in] coords : elemental coordinates
     * @param[out] records: records corresponding to the elemental matrix.
     *
     * */
    void getElementalMatrix(unsigned int eleID, std::vector<ot::MatRecord> &records, double *coords) {
        return asLeaf().getElementalMatrix(eleID, records, coords);
    }


};

template<typename T>
feMatrix<T>::feMatrix(ot::DA *da, unsigned int dof) : feMat(da) {
    m_uiDof = dof;
    const unsigned int nPe = m_uiOctDA->getNumNodesPerElement();
    m_uiEleVecIn = new VECType[m_uiDof * nPe];
    m_uiEleVecOut = new VECType[m_uiDof * nPe];

    m_uiEleCoords = new double[m_uiDim * nPe];

}

template<typename T>
feMatrix<T>::~feMatrix() {
    delete[] m_uiEleVecIn;
    delete[] m_uiEleVecOut;
    delete[] m_uiEleCoords;

    m_uiEleVecIn = NULL;
    m_uiEleVecOut = NULL;
    m_uiEleCoords = NULL;

}

template<typename T>
void feMatrix<T>::matVec(const VECType *in, VECType *out, double scale) {

    VECType *_in = NULL;
    VECType *_out = NULL;

    if (!(m_uiOctDA->isActive()))
        return;


    preMatVec(in, out, scale);

    m_uiOctDA->nodalVecToGhostedNodal(in, _in, false, m_uiDof);
    m_uiOctDA->createVector(_out, false, true, m_uiDof);

    VECType *val = new VECType[m_uiDof];
    for (unsigned int var = 0; var < m_uiDof; var++)
        val[var] = (VECType) 0;

    m_uiOctDA->setVectorByScalar(_out, val, false, true, m_uiDof);

    delete[] val;

    const ot::Mesh * pMesh = m_uiOctDA->getMesh();
    const unsigned int nPe = pMesh->getNumNodesPerElement();
    double * qMat = new double[nPe*nPe];
    const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
    const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
    
    m_uiOctDA->readFromGhostBegin(_in, m_uiDof);
    

    const unsigned int totalNodalSize = m_uiOctDA->getTotalNodalSz();
    for (m_uiOctDA->init<ot::DA_FLAGS::INDEPENDENT>(); m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::INDEPENDENT>(); m_uiOctDA->next<ot::DA_FLAGS::INDEPENDENT>()) {

        m_uiOctDA->getElementNodalValues(_in, m_uiEleVecIn, m_uiOctDA->curr(), m_uiDof);
        const unsigned int currentId = m_uiOctDA->curr();
        pMesh->getElementQMat(m_uiOctDA->curr(),qMat,true);
        m_uiOctDA->getElementalCoords(m_uiOctDA->curr(), m_uiEleCoords);
        elementalMatVec(m_uiEleVecIn, m_uiEleVecOut, m_uiEleCoords, scale);

        for (unsigned int dof = 0; dof < m_uiDof; dof++) {
          for (unsigned int i = 0; i < nPe; i++) {
            VECType ss = (VECType) 0;
            for (unsigned int j = 0; j < nPe; j++) {
              ss += qMat[j * nPe + i] * m_uiEleVecOut[dof*nPe + j];
            }
            _out[dof*totalNodalSize + e2n_cg[currentId * nPe + i]] += (VECType) ss;

          }
        }

    }

    m_uiOctDA->readFromGhostEnd(_in, m_uiDof);
    
    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd = pMesh -> getElementLocalEnd();

    for (m_uiOctDA->init<ot::DA_FLAGS::W_DEPENDENT>(); m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::W_DEPENDENT>(); m_uiOctDA->next<ot::DA_FLAGS::W_DEPENDENT>()) {

        // temporary fix to skip ghost writable.     
        if( m_uiOctDA->curr()< eleLocalBegin || m_uiOctDA->curr()>=eleLocalEnd )
            continue;

        m_uiOctDA->getElementNodalValues(_in, m_uiEleVecIn, m_uiOctDA->curr(), m_uiDof);
        const unsigned int currentId = m_uiOctDA->curr();
        pMesh->getElementQMat(m_uiOctDA->curr(),qMat,true);
        m_uiOctDA->getElementalCoords(m_uiOctDA->curr(), m_uiEleCoords);
        elementalMatVec(m_uiEleVecIn, m_uiEleVecOut, m_uiEleCoords, scale);

        for (unsigned int dof = 0; dof < m_uiDof; dof++) {
          for (unsigned int i = 0; i < nPe; i++) {
            VECType ss = (VECType) 0;
            for (unsigned int j = 0; j < nPe; j++) {
              ss += qMat[j * nPe + i] * m_uiEleVecOut[dof*nPe + j];
            }
            _out[dof*totalNodalSize + e2n_cg[currentId * nPe + i]] += (VECType) ss;

          }
        }

    }


    delete [] qMat;

    // accumilate from write ghost. 
    m_uiOctDA->writeToGhostsBegin(_out,m_uiDof);
    m_uiOctDA->writeToGhostsEnd(_out,ot::DA_FLAGS::WriteMode::ADD_VALUES,m_uiDof);
    
    m_uiOctDA->ghostedNodalToNodalVec(_out, out, true, m_uiDof);

    m_uiOctDA->destroyVector(_in);
    m_uiOctDA->destroyVector(_out);

    postMatVec(in, out, scale);


    return;

}

template<typename T>
void feMatrix<T>::getMatDiagonal(VECType*& diag, double scale)
{
    VECType *_diag = NULL;
    _diag = NULL;
    if (!(m_uiOctDA->isActive()))
        return;
    m_uiOctDA->createVector(_diag, false, true, m_uiDof);
    VECType *val = new VECType[m_uiDof];
    for (unsigned int var = 0; var < m_uiDof; var++)
        val[var] = (VECType) 0;
    m_uiOctDA->setVectorByScalar(_diag, val, false, true, m_uiDof);
    delete[] val;
    const unsigned int eleOrder = m_uiOctDA->getElementOrder();
    const unsigned int npe_1d = eleOrder + 1;
    const unsigned int npe_2d = (eleOrder + 1) * (eleOrder + 1);
    const unsigned int nPe = (eleOrder + 1) * (eleOrder + 1) * (eleOrder + 1);
    const ot::Mesh *pMesh = m_uiOctDA->getMesh();
    const ot::TreeNode *allElements = &(*(pMesh->getAllElements().begin()));
    DendroScalar *p2cEleMat = new DendroScalar[nPe * nPe];
    preMat();
    unsigned int nCount = 0;
    double *coords = new double[m_uiDim * nPe];
    bool faceHang[NUM_FACES];
    bool edgeHang[NUM_EDGES];
    unsigned int cnumFace[NUM_FACES];
    unsigned int cnumEdge[NUM_EDGES];
    const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
    const unsigned int * e2n_dg = &(*(pMesh->getE2NMapping_DG().begin()));
    const unsigned int * e2e = &(*(pMesh->getE2EMapping().begin()));
    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd = pMesh -> getElementLocalEnd();
    const unsigned int totalNodes = pMesh->getDegOfFreedom();
    for (m_uiOctDA->init<ot::DA_FLAGS::LOCAL_ELEMENTS>();m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::LOCAL_ELEMENTS>(); m_uiOctDA->next<ot::DA_FLAGS::LOCAL_ELEMENTS>())
    {
        std::vector<ot::MatRecord> records;
        const unsigned int currentId = m_uiOctDA->curr();
        m_uiOctDA->getElementalCoords(currentId, coords);
        getElementalMatrix(m_uiOctDA->curr(), records, coords);
        const unsigned int cnum = allElements[currentId].getMortonIndex();
        pMesh->getElementQMat(m_uiOctDA->curr(),p2cEleMat,true);
        //printArray_2D(p2cEleMat,nPe,nPe);
        fem_eMatTogMat(records.data() , p2cEleMat, eleOrder,m_uiDof);
        for(unsigned int i = 0; i < records.size(); i++)
        {
            if( (records[i].getRowID() == records[i].getColID()) && (records[i].getRowDim() == records[i].getColDim())) {
              _diag[records[i].getRowDim() * totalNodes + records[i].getRowID()] += (scale * records[i].getMatVal());
            }
        }
    }
    postMat();
    m_uiOctDA->writeToGhostsBegin(_diag,m_uiDof);
    m_uiOctDA->writeToGhostsEnd(_diag,ot::DA_FLAGS::WriteMode::ADD_VALUES,m_uiDof);
    delete [] coords;
    delete [] p2cEleMat;
    m_uiOctDA->ghostedNodalToNodalVec(_diag,diag,true,m_uiDof);
    delete [] _diag;
    return;
}


template <typename T>
void feMatrix<T>::getMatBlockDiagonal(std::vector<ot::MatRecord>& blkDiag, double scale)
{
    blkDiag.clear();
    if (!(m_uiOctDA->isActive()))
        return;
    
    const unsigned int activeNpes = m_uiOctDA->getNpesActive();

    const unsigned int eleOrder = m_uiOctDA->getElementOrder();
    const unsigned int npe_1d = eleOrder + 1;
    const unsigned int npe_2d = (eleOrder + 1) * (eleOrder + 1);
    const unsigned int nPe = (eleOrder + 1) * (eleOrder + 1) * (eleOrder + 1);
    
    const ot::Mesh *pMesh = m_uiOctDA->getMesh();
    const ot::TreeNode *allElements = &(*(pMesh->getAllElements().begin()));
    
    DendroScalar *p2cEleMat = new DendroScalar[nPe * nPe];
    preMat();

    unsigned int nCount = 0;
    double * coords = new double[m_uiDim * nPe];
    
    bool faceHang[NUM_FACES];
    bool edgeHang[NUM_EDGES];
    unsigned int cnumFace[NUM_FACES];
    unsigned int cnumEdge[NUM_EDGES];

    const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
    const unsigned int * e2n_dg = &(*(pMesh->getE2NMapping_DG().begin()));
    const unsigned int * e2e = &(*(pMesh->getE2EMapping().begin()));
    const DendroIntL * l2g = m_uiOctDA->getNodeLocalToGlobalMap().data();
    const DendroIntL * nodalOffsets = m_uiOctDA->getNodalOffsets().data();

    int * _n2p = NULL;
    m_uiOctDA->createVector(_n2p,false, true, 1);

    for(unsigned int i = pMesh->getNodeLocalBegin(); i < pMesh->getNodeLocalEnd(); i++)
        _n2p [i] = pMesh->getMPIRank();
    
    m_uiOctDA->readFromGhostBegin(_n2p,1);
    m_uiOctDA->readFromGhostEnd(_n2p,1);
    
    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd = pMesh -> getElementLocalEnd();
    const unsigned int totalNodes = pMesh->getDegOfFreedom();

    int * rSendCount  = new int[activeNpes];
    int * rRecvCount  = new int[activeNpes];
    int * rSendOffset = new int[activeNpes];
    int * rRecvOffset = new int[activeNpes];

    for(unsigned int i=0; i< activeNpes; i++)
        rSendCount[i]=0;

    for (m_uiOctDA->init<ot::DA_FLAGS::LOCAL_ELEMENTS>();m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::LOCAL_ELEMENTS>(); m_uiOctDA->next<ot::DA_FLAGS::LOCAL_ELEMENTS>())
    {
        const unsigned int currentId = m_uiOctDA->curr();
        for(unsigned int i=0; i < nPe; i++)
        for(unsigned int j=0; j < nPe; j++)
        {
            if(_n2p[ e2n_cg[currentId*nPe + i] ] == _n2p[ e2n_cg[currentId*nPe + j] ])
              rSendCount[ _n2p[ e2n_cg[currentId*nPe + i] ] ] ++;  
        }
    }

    for(unsigned int i=0; i< activeNpes; i++)
        rSendCount[i]*=m_uiDof;

    rSendOffset[0] =0;
    omp_par::scan(rSendCount, rSendOffset,activeNpes);

    std::vector<ot::MatRecord> sendBuf;
    sendBuf.resize(rSendOffset[activeNpes-1] + rSendCount[activeNpes-1]);

    for(unsigned int i=0; i< activeNpes; i++)
        rSendCount[i] = 0;

    
    for (m_uiOctDA->init<ot::DA_FLAGS::LOCAL_ELEMENTS>();m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::LOCAL_ELEMENTS>(); m_uiOctDA->next<ot::DA_FLAGS::LOCAL_ELEMENTS>()) 
    {
        std::vector<ot::MatRecord> records;
        records.reserve(nPe*m_uiDof * nPe * m_uiDof);
        const unsigned int currentId = m_uiOctDA->curr();
        m_uiOctDA->getElementalCoords(currentId, coords);
        getElementalMatrix(m_uiOctDA->curr(), records, coords);
        const unsigned int cnum = allElements[currentId].getMortonIndex();

        pMesh->getElementQMat(m_uiOctDA->curr(),p2cEleMat,true);
        //printArray_2D(p2cEleMat,nPe,nPe);
        fem_eMatTogMat(records.data() , p2cEleMat, eleOrder,m_uiDof);

        for(unsigned int i = 0; i < records.size(); i++)
        {
            if(_n2p[ records[i].getRowID() ] == _n2p[ records[i].getColID()])
            {   
                const unsigned int rid = records[i].getRowID();
                const unsigned int cid = records[i].getColID();
                const DendroIntL rid_g = l2g[rid];
                const DendroIntL cid_g = l2g[cid];
                const unsigned int proc_owner = _n2p[rid];
                assert(rid_g >= nodalOffsets[proc_owner]);
                const unsigned int rid_l = rid_g - nodalOffsets[proc_owner];
                const unsigned int cid_l = cid_g - nodalOffsets[proc_owner];

                records[i].setRowID(rid_l);
                records[i].setColID(cid_l);
                
                sendBuf[rSendCount[_n2p[rid]]] = records[i];
                rSendCount[_n2p[ rid ] ] ++;
            }
                
        }

    }

    // for(unsigned int i=0; i< activeNpes; i++)
    //     std::cout<<" i: "<<i<<" rSendCount: "<<rSendCount[i]<<std::endl;

    postMat();
    par::Mpi_Alltoall(rSendCount,rRecvCount,1,m_uiOctDA->getCommActive());
    
    rRecvOffset[0] = 0;
    omp_par::scan(rRecvCount,rRecvOffset,activeNpes);

    std::vector<ot::MatRecord> recvBuf;
    recvBuf.resize( rRecvOffset[activeNpes-1] + rRecvCount[activeNpes-1]);
    par::Mpi_Alltoallv(sendBuf.data(),rSendCount,rSendOffset,recvBuf.data(),rRecvCount,rRecvOffset,m_uiOctDA->getCommActive());

    delete [] rSendCount;
    delete [] rRecvCount;
    delete [] rSendOffset;
    delete [] rRecvOffset;
    delete [] _n2p;
    delete [] coords;
    delete [] p2cEleMat;
    
    std::sort(recvBuf.begin(),recvBuf.end());
    blkDiag.clear();

    for(unsigned int i = 0; i < recvBuf.size(); i++)
    {
        unsigned int curr =i;
        while( ( (i+1) < recvBuf.size() ) && (recvBuf[i].getRowID() == recvBuf[i + 1].getRowID()) && ( recvBuf[i].getColID() == recvBuf[i+1].getColID()) && (recvBuf[i].getRowDim() == recvBuf[i + 1].getRowDim()) && ( recvBuf[i].getColDim() == recvBuf[i+1].getColDim()) )
        {
            recvBuf[i].setMatValue(recvBuf[i].getMatVal() + recvBuf[i+1].getMatVal());
            i=i+1;
        }
       blkDiag.push_back(recvBuf[curr]);
    }

    return;

}


#ifdef BUILD_WITH_PETSC

template<typename T>
void feMatrix<T>::matVec(const Vec &in, Vec &out, double scale) {

    const PetscScalar *inArry = NULL;
    PetscScalar *outArry = NULL;

    VecGetArrayRead(in, &inArry);
    VecGetArray(out, &outArry);

    matVec(inArry, outArry, scale);

    VecRestoreArrayRead(in, &inArry);
    VecRestoreArray(out, &outArry);

}


template<typename T>
bool feMatrix<T>::getAssembledMatrix(Mat *J, MatType mtype) {
    

    if(m_uiOctDA->isActive())
    {
        MatZeroEntries(*J);
        std::vector<ot::MatRecord> records;
        //
        const unsigned int eleOrder = m_uiOctDA->getElementOrder();
        const unsigned int npe_1d = eleOrder + 1;
        const unsigned int npe_2d = (eleOrder + 1) * (eleOrder + 1);
        const unsigned int nPe = (eleOrder + 1) * (eleOrder + 1) * (eleOrder + 1);
        
        const ot::Mesh *pMesh = m_uiOctDA->getMesh();
        const ot::TreeNode *allElements = &(*(pMesh->getAllElements().begin()));
        
        DendroScalar *p2cEleMat = new DendroScalar[nPe * nPe];
        
        preMat();
        unsigned int nCount = 0;

        double *coords = new double[m_uiDim * nPe];
        
        bool faceHang[NUM_FACES];
        bool edgeHang[NUM_EDGES];
        unsigned int cnumFace[NUM_FACES];
        unsigned int cnumEdge[NUM_EDGES];

        const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg = &(*(pMesh->getE2NMapping_DG().begin()));
        const unsigned int * e2e = &(*(pMesh->getE2EMapping().begin()));

        const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
        const unsigned int eleLocalEnd = pMesh -> getElementLocalEnd();

        for (m_uiOctDA->init<ot::DA_FLAGS::LOCAL_ELEMENTS>();m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::LOCAL_ELEMENTS>(); m_uiOctDA->next<ot::DA_FLAGS::LOCAL_ELEMENTS>()) {
            
            const unsigned int currentId = m_uiOctDA->curr();
            m_uiOctDA->getElementalCoords(currentId, coords);
            getElementalMatrix(m_uiOctDA->curr(), records, coords);
            const unsigned int cnum = allElements[currentId].getMortonIndex();
            
            pMesh->getElementQMat(m_uiOctDA->curr(),p2cEleMat,true);
            //printArray_2D(p2cEleMat,nPe,nPe);


            fem_eMatTogMat(&(*(records.begin() + nCount)), p2cEleMat, eleOrder,m_uiDof);
            nCount += (m_uiDof*nPe*nPe*m_uiDof);
            if (records.size() > 500) {
                m_uiOctDA->petscSetValuesInMatrix(*J, records, m_uiDof, ADD_VALUES);
                nCount = 0;
            }

        }//end writable

        m_uiOctDA->petscSetValuesInMatrix(*J, records, m_uiDof, ADD_VALUES);

        postMat();

        MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);

        delete[] p2cEleMat;
        delete[] coords;

        PetscFunctionReturn(0);
    }

    
}


template<typename T>
PetscInt feMatrix<T>::getMatDiagonal(Vec& vec, PetscScalar scale)
{
    m_uiOctDA->petscCreateVector(vec,false,false,m_uiDof);
    VECType * diag =NULL;
    VecGetArray(vec,&diag);
    this->getMatDiagonal(diag,(double)scale);
    VecRestoreArray(vec,&diag);
    return 0; 
}

template<typename T>
PetscInt feMatrix<T>::getMatBlockDiagonal(Mat m, PetscScalar scale)
{
    std::vector<ot::MatRecord> blkDiag;
    this->getMatBlockDiagonal(blkDiag,(double)scale);

    for(unsigned int i=0; i< blkDiag.size(); i++)
    {
        blkDiag[i].setRowID(blkDiag[i].getRowID() + m_uiOctDA->getPreNodalSz());
        blkDiag[i].setColID(blkDiag[i].getColID() + m_uiOctDA->getPreNodalSz());
    }

    m_uiOctDA->petscSetValuesInMatrix(m,blkDiag,m_uiDof,ADD_VALUES);
    
    return 0;
    
}

#endif


#endif //DENDRO_5_0_FEMATRIX_H
