//
// Created by milinda on 12/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains FEM discretization operators for well-known PDEs, this is just to test if the FEM solver (i.e. matvec) works
 * correctly.
 * Later we need to write a SymPy interface which convert a given PDE to it's FEM discretization.
 *
 *
*/
//

#ifndef SFCSORTBENCH_OPERATORS_H
#define SFCSORTBENCH_OPERATORS_H

#include "TreeNode.h"
#include "refel.h"
#include <functional>
#include "tensor.h"
#include "workspace.h"
#include "mathUtils.h"
#include "dendroProfileParams.h"

#define GRIDX_TO_X(xg) (((Rx/RgX)*(xg-bssn::BSSN_OCTREE_MIN[0]))+bssn::BSSN_COMPD_MIN[0])
#define GRIDY_TO_Y(yg) (((Ry/RgY)*(yg-bssn::BSSN_OCTREE_MIN[1]))+bssn::BSSN_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg) (((Rz/RgZ)*(zg-bssn::BSSN_OCTREE_MIN[2]))+bssn::BSSN_COMPD_MIN[2])

#define X_TO_GRIDX(xc) (((RgX/Rx)*(xc-bssn::BSSN_COMPD_MIN[0]))+bssn::BSSN_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc) (((RgY/Ry)*(yc-bssn::BSSN_COMPD_MIN[1]))+bssn::BSSN_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc) (((RgZ/Rz)*(zc-bssn::BSSN_COMPD_MIN[2]))+bssn::BSSN_OCTREE_MIN[2])




namespace fem
{
    namespace operators
    {

        namespace poisson
        {

            /**
             * PDE:\f{eqnarray*}{
                     \nabla u(x) = f(x) \text{ for } x \in \Gamma \\
                     u(x)=0 \text{ for } x \in \partial\Gamma
                   \f}

             * Galerkin discretization Ku=Mf (computes the Mf)*
             * @brief computes the RHS (or mass matvec) for a given mesh element.
             *
             * */
            template <typename  T>
            void elementMassMatvec(const T* f_rhs, const ot::TreeNode & elem, const RefElement * refEl, T* out)
            {

#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_elemMvec.start();
#endif

                const double * Q1d=refEl->getQ1d();
                const double * QT1d=refEl->getQT1d();
                const double * Dg=refEl->getDg1d();
                const double * W1d=refEl->getWgq();

                const unsigned int eleOrder=refEl->getOrder();
                const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
                const unsigned int nrp=eleOrder+1;

                const double szX=fem::domain::gridX_to_X(elem.maxX())-fem::domain::gridX_to_X(elem.minX());
                const double szY=fem::domain::gridY_to_Y(elem.maxY())-fem::domain::gridY_to_Y(elem.minY());
                const double szZ=fem::domain::gridZ_to_Z(elem.maxZ())-fem::domain::gridZ_to_Z(elem.minZ());


                const double refElSz=refEl->getElementSz();

                // interpolate to quadrature points.
                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,f_rhs,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,out);



                const double Jx = 1.0/(refElSz/(double (szX)));
                const double Jy = 1.0/(refElSz/(double (szY)));
                const double Jz = 1.0/(refElSz/(double (szZ)));

                //std::cout<<"Mass:  elem: "<<elem<<" ele Sz: "<<(elem.maxX()-elem.minX())<<" szX: "<<szX<<" Jx: "<<Jx<<" J: "<<(Jx*Jy*Jz)<<std::endl;

                for(unsigned int k=0;k<(eleOrder+1);k++)
                    for(unsigned int j=0;j<(eleOrder+1);j++)
                       for(unsigned int i=0;i<(eleOrder+1);i++)
                           out[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=(Jx*Jy*Jz*W1d[i]*W1d[j]*W1d[k]);


                // apply transpose operator
                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,out,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,out);

#ifdef FEM_ACCUMILATE_ONES_TEST
                for(unsigned int k=0;k<(eleOrder+1);k++)
                    for(unsigned int j=0;j<(eleOrder+1);j++)
                        for(unsigned int i=0;i<(eleOrder+1);i++)
                            out[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]=1.0;

#endif
                /*if(elem.minX()==0 && elem.minY()==0 && elem.minZ()==0)
                {
                    std::cout<<"elem: "<<elem<<std::endl;
                    for(unsigned int k=0;k<(eleOrder+1);k++)
                      for(unsigned int j=0;j<(eleOrder+1);j++)
                        for(unsigned int i=0;i<(eleOrder+1);i++)
                            std::cout<<"i: "<<i<<" j: "<<j<<" k: "<<k<<" MF "<<out[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]<<" F "<<f_rhs[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]<<std::endl;

                }*/

#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_elemMvec.stop();
#endif

            }

            /**
             * @brief computes the stiffness elemental matvec (i.e. LHS ) of CG discretization
             *
             * */
            template <typename T>
            void elementStiffnessMatvec(const T* in, const ot::TreeNode & elem, const RefElement * refEl, T* out)
            {

#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_elemMvec.start();
#endif

                const double * Q1d=refEl->getQ1d();
                const double * QT1d=refEl->getQT1d();
                const double * Dg=refEl->getDg1d();
                const double * DgT=refEl->getDgT1d();
                const double * W1d=refEl->getWgq();

                const unsigned int eleOrder=refEl->getOrder();
                const unsigned int nPe=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
                const unsigned int nrp=eleOrder+1;

                const unsigned int sz=1u<<(m_uiMaxDepth-elem.getLevel()); // curent element size
                const double refElSz=refEl->getElementSz();
                //x derivative
                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Dg,in,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qx);

                //y derivative
                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Dg,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Q1d,imV2,Qy);

                //z derivative
                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,Q1d,in,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,Q1d,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,Dg,imV2,Qz);


                const double szX=fem::domain::gridX_to_X(elem.maxX())-fem::domain::gridX_to_X(elem.minX());
                const double szY=fem::domain::gridY_to_Y(elem.maxY())-fem::domain::gridY_to_Y(elem.minY());
                const double szZ=fem::domain::gridZ_to_Z(elem.maxZ())-fem::domain::gridZ_to_Z(elem.minZ());

                const double Jx = 1.0/(refElSz/(double (szX)));
                const double Jy = 1.0/(refElSz/(double (szY)));
                const double Jz = 1.0/(refElSz/(double (szZ)));

                //std::cout<<"Stifness:  elem: "<<elem<<" ele Sz: "<<(elem.maxX()-elem.minX())<<" szX: "<<szX<<" Jx: "<<Jx<<" J: "<<(Jx*Jy*Jz)<<std::endl;

                for(unsigned int k=0;k<(eleOrder+1);k++)
                    for(unsigned int j=0;j<(eleOrder+1);j++)
                        for(unsigned int i=0;i<(eleOrder+1);i++)
                        {
                            Qx[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jy*Jz)/Jx)*W1d[i]*W1d[j]*W1d[k]);
                            Qy[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jz)/Jy)*W1d[i]*W1d[j]*W1d[k]);
                            Qz[k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i]*=( ((Jx*Jy)/Jz)*W1d[i]*W1d[j]*W1d[k]);
                        }



                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,DgT,Qx,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qx);

                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qy,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,DgT,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,QT1d,imV2,Qy);

                DENDRO_TENSOR_IIAX_APPLY_ELEM(nrp,QT1d,Qz,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(nrp,QT1d,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(nrp,DgT,imV2,Qz);

                for(unsigned int i=0;i<nPe;i++)
                    out[i]=Qx[i]+Qy[i]+Qz[i];


#ifdef MATVEC_PROFILE
                dendro::timer::sfcmatvec::t_elemMvec.stop();
#endif


            }



            template <typename T>
            void applyDirichletBoundary(const ot::Mesh* mesh,T* in,T bdyValue)
            {

                const unsigned int nPe=mesh->getNumNodesPerElement();
                const unsigned int eleOrder=mesh->getElementOrder();
                const ot::TreeNode* pNodes = &(*(mesh->getAllElements().begin()));

                const unsigned int * e2n=&(*(mesh->getE2NMapping().begin()));
                const unsigned int * e2n_dg=&(*(mesh->getE2NMapping_DG().begin()));

                const unsigned int d_min=0;
                const unsigned int d_max=(1u<<m_uiMaxDepth);


                // apply zero-dirichlet bdy conditions here.
                const unsigned int nx=(eleOrder+1);
                const unsigned int ny=(eleOrder+1)*(eleOrder+1);
                const unsigned int nz=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
                const double bdy_value=bdyValue;
                double x,y,z;

                for(unsigned int ele=mesh->getElementLocalBegin();ele<mesh->getElementLocalEnd();ele++)
                {
                    x=pNodes[ele].minX();
                    y=pNodes[ele].minY();
                    z=pNodes[ele].minZ();

                    if(x==d_min)
                    {
                        for(unsigned int k=0;k<(eleOrder+1);k++)
                            for(unsigned int j=0;j<(eleOrder+1);j++)
                            {
                                if( (!mesh->isNodeHanging(ele,0,j,k)) && (mesh->isNodeLocal(ele,0,j,k)))
                                    in[e2n[ele*nPe+k*ny+j*nx+0]]=bdy_value;
                            }
                    }

                    if(y==d_min)
                    {
                        for(unsigned int k=0;k<(eleOrder+1);k++)
                            for(unsigned int i=0;i<(eleOrder+1);i++)
                            {
                                if( (!mesh->isNodeHanging(ele,i,0,k)) && (mesh->isNodeLocal(ele,i,0,k)))
                                    in[e2n[ele*nPe+k*ny+0*nx+i]]=bdy_value;
                            }
                    }


                    if(z==d_min)
                    {
                        for(unsigned int j=0;j<(eleOrder+1);j++)
                            for(unsigned int i=0;i<(eleOrder+1);i++)
                            {
                                if( (!mesh->isNodeHanging(ele,i,j,0)) && (mesh->isNodeLocal(ele,i,j,0)))
                                    in[e2n[ele*nPe+0*ny+j*nx+i]]=bdy_value;
                            }
                    }

                    x=pNodes[ele].maxX();
                    y=pNodes[ele].maxY();
                    z=pNodes[ele].maxZ();

                    if(x==d_max)
                    {
                        for(unsigned int k=0;k<(eleOrder+1);k++)
                            for(unsigned int j=0;j<(eleOrder+1);j++)
                            {
                                if( (!mesh->isNodeHanging(ele,eleOrder,j,k)) && (mesh->isNodeLocal(ele,eleOrder,j,k)))
                                    in[e2n[ele*nPe+k*ny+j*nx+eleOrder]]=bdy_value;
                            }
                    }

                    if(y==d_max)
                    {
                        for(unsigned int k=0;k<(eleOrder+1);k++)
                            for(unsigned int i=0;i<(eleOrder+1);i++)
                            {
                                if( (!mesh->isNodeHanging(ele,i,eleOrder,k)) && (mesh->isNodeLocal(ele,i,eleOrder,k)))
                                    in[e2n[ele*nPe+k*ny+eleOrder*nx+i]]=bdy_value;
                            }
                    }


                    if(z==d_max)
                    {
                        for(unsigned int j=0;j<(eleOrder+1);j++)
                            for(unsigned int i=0;i<(eleOrder+1);i++)
                            {
                                if( (!mesh->isNodeHanging(ele,i,j,eleOrder)) && (mesh->isNodeLocal(ele,i,j,eleOrder)))
                                    in[e2n[ele*nPe+eleOrder*ny+j*nx+i]]=bdy_value;
                            }

                    }


                }
            }




            /**
             *@brief computes the mesh based matvec operation.
             *@param [in] mesh: mesh data structure.
             *@param [in] in: input vector x (v=Ax)
             *@param [out] out: output vector v of such that v=Ax
             *
            * **/
            template <typename T>
            void matvec(ot::Mesh* mesh,T* in, T* out,bool applyBoundary=true)
            {
                if(mesh->isActive())
                {

                    const std::vector<unsigned int> level1Ghost=mesh->getLevel1GhostElementIndices();

                    mesh->performGhostExchange(in);

                    const unsigned int nPe=mesh->getNumNodesPerElement();
                    const unsigned int eleOrder=mesh->getElementOrder();


                    T* elementalVec=new T[nPe];
                    T* elementalMVec=new T[nPe];

                    const unsigned int nx=(eleOrder+1);
                    const unsigned int ny=(eleOrder+1)*(eleOrder+1);
                    const unsigned int nz=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
                    const ot::TreeNode* pNodes = &(*(mesh->getAllElements().begin()));
                    const RefElement *refEl=mesh->getReferenceElement();

                    unsigned int nodeLookup,ownerID,ii_x,jj_y,kk_z;
                    for(unsigned int node=mesh->getNodeLocalBegin();node<mesh->getNodeLocalEnd();node++)
                        out[node]=(T)0;

                    mesh->performGhostExchange(out);

                    // iterate over ghost elements . (later we can use this to potentially overlap computation & communication)
                    for(unsigned int gIndex=0;gIndex<level1Ghost.size();gIndex++)
                    {
                        unsigned int ele=level1Ghost[gIndex];
                        mesh->getElementNodalValues(in,elementalVec,ele);
                        operators::poisson::elementStiffnessMatvec(elementalVec,pNodes[ele],refEl,elementalMVec);
                        mesh->computeElementalContribution(elementalMVec,out,ele);
                    }



                    // iterate over local elements.
                    for(unsigned int ele=mesh->getElementLocalBegin();ele<mesh->getElementLocalEnd();ele++)
                    {
                        mesh->getElementNodalValues(in,elementalVec,ele);

                      /*  if(ele==0)
                            for(unsigned int node=0;node<125;node++)
                                std::cout<<" node value in: "<<elementalVec[node]<<std::endl;*/

                        operators::poisson::elementStiffnessMatvec(elementalVec,pNodes[ele],refEl,elementalMVec);

                       /* if(ele==0)
                            for(unsigned int node=0;node<125;node++)
                                std::cout<<" node value out: "<<elementalMVec[node]<<std::endl;*/

                        mesh->computeElementalContribution(elementalMVec,out,ele);

                    }

                    delete [] elementalVec;
                    delete [] elementalMVec;
                    mesh->performGhostExchange(out);


                }
            }


            /**
             *@brief computes the mesh rhs vector
             *@param [in] mesh: mesh data structure.
             *@param [in] in: input vector x (v=Ax)
             *@param [out] out: output vector v of such that v=Ax
             *
            * **/
            template <typename T>
            void computeRHS(ot::Mesh* mesh,T* in, T* out)
            {


                if(mesh->isActive())
                {
                    //applyDirichletBoundary(mesh,in,0.0);
                    mesh->performGhostExchange(in);
                    const unsigned int nPe=mesh->getNumNodesPerElement();
                    const unsigned int eleOrder=mesh->getElementOrder();
                    const std::vector<unsigned int> level1Ghost=mesh->getLevel1GhostElementIndices();


                    T* elementalVec=new T[nPe];
                    T* elementalMVec=new T[nPe];

                    const unsigned int nx=(eleOrder+1);
                    const unsigned int ny=(eleOrder+1)*(eleOrder+1);
                    const unsigned int nz=(eleOrder+1)*(eleOrder+1)*(eleOrder+1);
                    const ot::TreeNode* pNodes = &(*(mesh->getAllElements().begin()));
                    const RefElement *refEl=mesh->getReferenceElement();

                    /*const double * W1d=refEl->getWgq();
                    for(unsigned int i=0;i<5;i++)
                        std::cout<<" w["<<i<<"]: "<<W1d[i]<<std::endl;*/

                    const unsigned int * e2n=&(*(mesh->getE2NMapping().begin()));
                    const unsigned int * e2n_dg=&(*(mesh->getE2NMapping_DG().begin()));
                    unsigned int nodeLookup,ownerID,ii_x,jj_y,kk_z;


                    for(unsigned int node=mesh->getNodeLocalBegin();node<mesh->getNodeLocalEnd();node++)
                        out[node]=(T)0.0;

                    mesh->performGhostExchange(out);

                    // iterate over ghost elements . (later we can use this to potentially overlap computation & communication)
                    for(unsigned int gIndex=0;gIndex<level1Ghost.size();gIndex++)
                    {
                        unsigned int ele=level1Ghost[gIndex];

#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_p2cInterp.start();
#endif
                        mesh->getElementNodalValues(in,elementalVec,ele);
#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_p2cInterp.stop();
#endif

                        operators::poisson::elementMassMatvec(elementalVec,pNodes[ele],refEl,elementalMVec);


#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_c2pInterp.start();
#endif

                        mesh->computeElementalContribution(elementalMVec,out,ele);

#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_c2pInterp.stop();
#endif

                    }

                    // iterate over local elements
                    for(unsigned int ele=mesh->getElementLocalBegin();ele<mesh->getElementLocalEnd();ele++)
                    {

#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_p2cInterp.start();
#endif
                        mesh->getElementNodalValues(in,elementalVec,ele);

#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_p2cInterp.stop();
#endif
                        operators::poisson::elementMassMatvec(elementalVec,pNodes[ele],refEl,elementalMVec);


#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_c2pInterp.start();
#endif
                        mesh->computeElementalContribution(elementalMVec,out,ele);

#ifdef MATVEC_PROFILE
                        dendro::timer::sfcmatvec::t_c2pInterp.stop();
#endif


                    }

                    mesh->performGhostExchange(out);

                    delete [] elementalVec;
                    delete [] elementalMVec;


                }

            }



        } // end of namespace possoin

    } // end of namespace operators


} // end of namespace fem



#endif //SFCSORTBENCH_OPERATORS_H
