//
// Created by milinda on 8/20/18.
//

#ifndef DENDRO_5_0_GPUBSSNEXAMPLE_H
#define DENDRO_5_0_GPUBSSNEXAMPLE_H

/**
 * @author Milinda Fernando
 * @brief simple gr mesh generation and launch cuda kernels to compute BSSN rhs
 *
 * */



#include "gr.h"
#include "grUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rkBSSN.h"
#include "octUtils.h"
#include "mathUtils.h"

#ifdef BSSN_ENABLE_CUDA
    #include "rhs_cuda.cuh"
    #include "params_cu.h"
    #include "profile_gpu.h"
#endif

#include "rhs.h"



namespace cuda
{


    template <typename T>
    void applyIntilCondition(T** varIn,const ot::Mesh* m_uiMesh)
    {
        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        unsigned int x,y,z,len;
        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;
        unsigned int eleOrder=m_uiMesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();


        double* var=new double[bssn::BSSN_NUM_VARS];

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
                            len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                            x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));
                            assert(len%eleOrder==0);
                            if (bssn::BSSN_ID_TYPE == 0) {
                                TwoPunctures((double)x,(double)y,(double)z,var,
                                             &mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
                            }
                            else {
                                // all other values are handled in the initial data wrapper including
                                // an error message
                                initialDataFunctionWrapper((double)x, (double)y, (double)z, var);
                            }
                            for(unsigned int v=0;v<bssn::BSSN_NUM_VARS;v++)
                                varIn[v][nodeLookUp_CG]=var[v];


                        }

                    }

        }
    }

} // end of namespace cuda.





#endif //DENDRO_5_0_GPUBSSNEXAMPLE_H
