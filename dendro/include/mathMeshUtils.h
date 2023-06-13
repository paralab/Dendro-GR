/**
 * @file mathMeshUtils.h
 * @author: Milinda Fernando
 * @brief : Math utils using the mesh class. 
 * @version 0.1
 * @date 2019-11-07
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once

#include "mesh.h"
#include "mathUtils.h"

/**
 * @brief : Computes the volumetric norm $||v1||_p$ with Reimann Sum.  
 * @note assumes the vectors have exchange the ghosts. 
 * @tparam T : type of the the vector
 * @param pMesh : Underlying mesh data 
 * @param v1 : input vector 1
 * @param v2 : intpur vector 2
 * @param p : order of the norm. 
 */
template<typename T>
double rsNormLp(const ot::Mesh* pMesh, const T* v1, unsigned int p)
{
    if(!pMesh->isActive())
        return -1.0;
    
    const unsigned int nPe = pMesh->getNumNodesPerElement();
    std::vector<T> nodalVals;
    nodalVals.resize(nPe);

    double rs=0;
    double rs_g=0;
    
    const ot::TreeNode* pNodes = pMesh->getAllElements().data();
    const unsigned int  DD = pMesh->getElementOrder();

    double invDD = 1.0/((double)(DD));
    double invNC = 1.0/((double)(NUM_CHILDREN));

    Point pmin,pmax;

    for(unsigned int ele = pMesh->getElementLocalBegin(); ele < pMesh->getElementLocalEnd(); ele++)
    {

        pMesh->getElementNodalValues(v1,nodalVals.data(),ele);
        
        pMesh->octCoordToDomainCoord(Point(pNodes[ele].minX(),pNodes[ele].minY(),pNodes[ele].minZ()), pmin);
        pMesh->octCoordToDomainCoord(Point(pNodes[ele].maxX(),pNodes[ele].maxY(),pNodes[ele].maxZ()), pmax);

        const double hx = (pmax.x() - pmin.x())*invDD;
        const double hy = (pmax.y() - pmin.y())*invDD;
        const double hz = (pmax.z() - pmin.z())*invDD;

        for(unsigned int k=0; k < DD; k++)
         for(unsigned int j=0; j < DD; j++)
          for(unsigned int i=0; i < DD; i++)
          {
            
            double cell_avg =  (nodalVals[ (k) * (DD+1)*(DD+1) + (j)*(DD+1) + (i)]
                              + nodalVals[ (k) * (DD+1)*(DD+1) + (j)*(DD+1) + (i+1)]
                              + nodalVals[ (k) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i)]
                              + nodalVals[ (k) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i+1)]
                              + nodalVals[ (k+1) * (DD+1)*(DD+1) + (j)*(DD+1) + (i)]
                              + nodalVals[ (k+1) * (DD+1)*(DD+1) + (j)*(DD+1) + (i+1)]
                              + nodalVals[ (k+1) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i)]
                              + nodalVals[ (k+1) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i+1)])*invNC;
            
            rs += ( pow(cell_avg,p) * hx * hy * hz);

          } 

    }

    Point dmin = pMesh->getDomainMinPt();
    Point dmax = pMesh->getDomainMaxPt();

    double d_volume  = (dmax.x() - dmin.x()) * (dmax.y() - dmin.y()) * (dmax.z() - dmin.z());
    
    par::Mpi_Reduce(&rs,&rs_g,1,MPI_SUM,0,pMesh->getMPICommunicator());
    rs_g/= d_volume;
    rs_g = pow(rs_g,(1.0/p));
    return rs_g;


}


/**
 * @brief : Computes the volumetric norm $||v2-v1||_p$ with Reimann Sum.  
 * @note assumes the vectors have exchange the ghosts. 
 * @tparam T : type of the the vector
 * @param pMesh : Underlying mesh data 
 * @param v1 : input vector 1
 * @param v2 : intpur vector 2
 * @param p : order of the norm. 
 */
template<typename T>
double rsNormLp(const ot::Mesh* pMesh, const T* v1, const T* v2, unsigned int p)
{

    if(!pMesh->isActive())
        return -1.0;

    const unsigned int nPe = pMesh->getNumNodesPerElement();
    double * eleVec1 =new double[nPe];
    double * eleVec2 =new double[nPe];

    double rs=0;
    double rs_g=0;
    double hx;
    const ot::TreeNode* pNodes = pMesh->getAllElements().data();
    
    const unsigned int  DD = pMesh->getElementOrder();
    
    double invDD = 1.0/((double)(DD));
    double invNC = 1.0/((double)(NUM_CHILDREN));
    Point pmin,pmax;

    for(unsigned int e = pMesh->getElementLocalBegin(); e < pMesh->getElementLocalEnd(); e ++)
    {

        pMesh->getElementNodalValues(v1,eleVec1,e);
        pMesh->getElementNodalValues(v2,eleVec2,e);

        pMesh->octCoordToDomainCoord(Point(pNodes[e].minX(),pNodes[e].minY(),pNodes[e].minZ()), pmin);
        pMesh->octCoordToDomainCoord(Point(pNodes[e].maxX(),pNodes[e].maxY(),pNodes[e].maxZ()), pmax);

        const double hx = (pmax.x() - pmin.x())*invDD;
        const double hy = (pmax.y() - pmin.y())*invDD;
        const double hz = (pmax.z() - pmin.z())*invDD;

        for(unsigned int k=0; k < DD; k++)
         for(unsigned int j=0; j < DD; j++)
          for(unsigned int i=0; i < DD; i++)
          {
            
            double cell_avg1 = (eleVec1[ (k) * (DD+1)*(DD+1) + (j)*(DD+1) + (i)]
                              + eleVec1[ (k) * (DD+1)*(DD+1) + (j)*(DD+1) + (i+1)]
                              + eleVec1[ (k) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i)]
                              + eleVec1[ (k) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i+1)]
                              + eleVec1[ (k+1) * (DD+1)*(DD+1) + (j)*(DD+1) + (i)]
                              + eleVec1[ (k+1) * (DD+1)*(DD+1) + (j)*(DD+1) + (i+1)]
                              + eleVec1[ (k+1) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i)]
                              + eleVec1[ (k+1) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i+1)])*invNC;

            
            double cell_avg2 = (eleVec2[ (k) * (DD+1)*(DD+1) + (j)*(DD+1) + (i)]
                              + eleVec2[ (k) * (DD+1)*(DD+1) + (j)*(DD+1) + (i+1)]
                              + eleVec2[ (k) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i)]
                              + eleVec2[ (k) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i+1)]
                              + eleVec2[ (k+1) * (DD+1)*(DD+1) + (j)*(DD+1) + (i)]
                              + eleVec2[ (k+1) * (DD+1)*(DD+1) + (j)*(DD+1) + (i+1)]
                              + eleVec2[ (k+1) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i)]
                              + eleVec2[ (k+1) * (DD+1)*(DD+1) + (j+1)*(DD+1) + (i+1)])*invNC;

            double cell_avg = cell_avg1-cell_avg2;
            rs += ( pow(cell_avg,p) * hx * hy * hz);

          }
        
    }

    delete [] eleVec1;
    delete [] eleVec2;

    Point dmin = pMesh->getDomainMinPt();
    Point dmax = pMesh->getDomainMaxPt();

    double d_volume  = (dmax.x() - dmin.x()) * (dmax.y() - dmin.y()) * (dmax.z() - dmin.z());
    //std::cout<<"d_volume: "<<d_volume<<" rs: "<<rs<<std::endl;
    
    par::Mpi_Reduce(&rs,&rs_g,1,MPI_SUM,0,pMesh->getMPICommunicator());
    rs_g/= d_volume;
    rs_g = pow(rs_g,(1.0/p));
    return rs_g;
    

}

/**
 * @brief compute the min over the mesh vector. 
 * 
 * @tparam T data type. 
 * @param pMesh : pointer to the underlying mesh. 
 * @param v1 : input vector
 * @param vtype : type of the vector. 
 * isGhosted : true if allocated with ghosted false otherwise. 
 * @return T // rank 0 of the active comm with have the reduced value. 
 */
template<typename T>
T vecMin(const ot::Mesh *pMesh, const T* v1, ot::VEC_TYPE vtype, bool isGhosted=true)
{
    if(!(pMesh->isActive()))
        return NAN;

    unsigned int offset =0; 
    unsigned int nsz    =0;
    MPI_Comm aComm = pMesh->getMPICommunicator(); // active comm. 
    unsigned int arank = pMesh->getMPIRank();
    T val = 0;
    if(vtype == ot::VEC_TYPE::CG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getNodeLocalBegin();

        nsz = pMesh->getNumLocalMeshNodes();
        val = vecMin((T*)(v1 + offset), nsz, aComm);


    }else if(vtype == ot::VEC_TYPE::DG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin() * pMesh->getNumNodesPerElement();

        nsz = pMesh->getNumLocalMeshElements() * pMesh->getNumNodesPerElement();
        val = vecMin((T*)(v1 + offset), nsz, aComm);

    }else if(vtype == ot::VEC_TYPE::ELEMENTAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin();

        nsz = pMesh->getNumLocalMeshElements();
        val = vecMin((T*)(v1 + offset), nsz, aComm);


    }else 
    {
        std::cout<<"rank: "<<arank<<" invalid vec type to the norm call"<<__LINE__<<std::endl;
        val=NAN;
    }

    return val;


}

/**
 * @brief compute the max over the mesh vector. 
 * 
 * @tparam T data type. 
 * @param pMesh : pointer to the underlying mesh. 
 * @param v1 : input vector
 * @param vtype : type of the vector. 
 * @return T 
 */
template<typename T>
T vecMax(const ot::Mesh *pMesh, const T* v1, ot::VEC_TYPE vtype, bool isGhosted=true)
{

    if(!(pMesh->isActive()))
        return NAN;

    unsigned int offset =0; 
    unsigned int nsz    =0;
    MPI_Comm aComm = pMesh->getMPICommunicator(); // active comm. 
    unsigned int arank = pMesh->getMPIRank();
    T val = 0;
    if(vtype == ot::VEC_TYPE::CG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getNodeLocalBegin();

        nsz = pMesh->getNumLocalMeshNodes();
        val = vecMax((T*)(v1 + offset), nsz, aComm);


    }else if(vtype == ot::VEC_TYPE::DG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin() * pMesh->getNumNodesPerElement();

        nsz = pMesh->getNumLocalMeshElements() * pMesh->getNumNodesPerElement();
        val = vecMax((T*)(v1 + offset), nsz, aComm);

    }else if(vtype == ot::VEC_TYPE::ELEMENTAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin();

        nsz = pMesh->getNumLocalMeshElements();
        val = vecMax((T*)(v1 + offset), nsz, aComm);


    }else 
    {
        std::cout<<"rank: "<<arank<<" invalid vec type to the norm call"<<__LINE__<<std::endl;
        val=NAN;
    }

    return val;
    
}

/**
 * @brief compute the L2 (discreate) over the mesh vector. 
 * 
 * @tparam T data type. 
 * @param pMesh : pointer to the underlying mesh. 
 * @param v1 : input vector
 * @param vtype : type of the vector. 
 * @return T 
 */
template<typename T>
T normL2(const ot::Mesh *pMesh, const T* v1, ot::VEC_TYPE vtype, bool isGhosted=true)
{
    if(!(pMesh->isActive()))
        return NAN;

    unsigned int offset =0; 
    unsigned int nsz    =0;
    MPI_Comm aComm = pMesh->getMPICommunicator(); // active comm. 
    unsigned int arank = pMesh->getMPIRank();
    T val = 0;

    DendroIntL dof_g=0;
    DendroIntL dof_l;

    if(vtype == ot::VEC_TYPE::CG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getNodeLocalBegin();

        dof_l = pMesh->getNumLocalMeshNodes();

        nsz = pMesh->getNumLocalMeshNodes();
        val = normL2((T*)(v1 + offset), nsz, aComm);


    }else if(vtype == ot::VEC_TYPE::DG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin() * pMesh->getNumNodesPerElement();

        dof_l = pMesh->getNumLocalMeshElements() * pMesh->getNumNodesPerElement();

        nsz = pMesh->getNumLocalMeshElements() * pMesh->getNumNodesPerElement();
        val = normL2((T*)(v1 + offset), nsz, aComm);

    }else if(vtype == ot::VEC_TYPE::ELEMENTAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin();

        dof_l = pMesh->getNumLocalMeshElements();

        nsz = pMesh->getNumLocalMeshElements();
        val = normL2((T*)(v1 + offset), nsz, aComm);


    }else 
    {
        dof_l=0;
        std::cout<<"rank: "<<arank<<" invalid vec type to the norm call"<<__LINE__<<std::endl;
        val=NAN;
    }

    par::Mpi_Reduce(&dof_l,&dof_g,1,MPI_SUM,0,aComm);
    val = val/(double)sqrt(dof_g);
    return val;
}

/**
 * @brief compute the l-infinity over the mesh vector. 
 * 
 * @tparam T data type. 
 * @param pMesh : pointer to the underlying mesh. 
 * @param v1 : input vector
 * @param vtype : type of the vector. 
 * @return T 
 */
template<typename T>
T normLInfty(const ot::Mesh *pMesh, const T* v1, ot::VEC_TYPE vtype, bool isGhosted=true)
{
     if(!(pMesh->isActive()))
        return NAN;

    unsigned int offset =0; 
    unsigned int nsz    =0;
    MPI_Comm aComm = pMesh->getMPICommunicator(); // active comm. 
    unsigned int arank = pMesh->getMPIRank();
    T val = 0;
    if(vtype == ot::VEC_TYPE::CG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getNodeLocalBegin();

        nsz = pMesh->getNumLocalMeshNodes();
        val = normLInfty((T*)(v1 + offset), nsz, aComm);


    }else if(vtype == ot::VEC_TYPE::DG_NODAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin() * pMesh->getNumNodesPerElement();

        nsz = pMesh->getNumLocalMeshElements() * pMesh->getNumNodesPerElement();
        val = normLInfty((T*)(v1 + offset), nsz, aComm);

    }else if(vtype == ot::VEC_TYPE::ELEMENTAL)
    {
        if(isGhosted)
            offset = pMesh->getElementLocalBegin();

        nsz = pMesh->getNumLocalMeshElements();
        val = normLInfty((T*)(v1 + offset), nsz, aComm);


    }else 
    {
        std::cout<<"rank: "<<arank<<" invalid vec type to the norm call"<<__LINE__<<std::endl;
        val=NAN;
    }

    return val;
}