//
// Created by milinda on 05/02/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for BSSN simulation.
*/

#include "mesh.h"
#include <math.h>
#include "daUtils.h"

namespace bssn
{

    template <typename T>
    double computeConstraintL2Norm(const T* constraintVec, const T* maskVec, unsigned int lbegin, unsigned int lend,MPI_Comm comm)
    {

        double l2=0.0;
        double l2_g=0.0;
        const double MASK_THRESHOLD=1.0;
        for(unsigned int i=lbegin;i<lend;i++)
        {
            if(maskVec[i]<MASK_THRESHOLD)
                l2+=(constraintVec[i]*constraintVec[i]*maskVec[i]*maskVec[i]*maskVec[i]*maskVec[i]);
            else
                l2+=(constraintVec[i]*constraintVec[i]);
        }


        par::Mpi_Reduce(&l2,&l2_g,1,MPI_SUM,0,comm);

        return (sqrt(l2_g));


    }

    
    template <typename T>
    double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec,const T* maskVector,T maskthreshoold)
    {
        double l2_g=0.0;
        if(mesh->isActive())
        {
            
            MPI_Comm comm=mesh->getMPICommunicator();
            const unsigned int eleLocalBegin   =  mesh->getElementLocalBegin();
            const unsigned int eleLocalEnd     =  mesh->getElementLocalEnd();
            
            const unsigned int nodeLocalBegin  =  mesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd    =  mesh->getNodeLocalEnd();
            const unsigned int nPe             =  mesh->getNumNodesPerElement();
            
            const unsigned int * e2n_cg        =  &(*(mesh->getE2NMapping().begin()));
            const unsigned int * e2n_dg        =  &(*(mesh->getE2NMapping_DG().begin()));
            const unsigned int eleOrder = mesh->getElementOrder();
            const ot::TreeNode* pNodes  = mesh->getAllElements().data();


            double l2=0.0;
            
            DendroIntL localGridPts=0;
            DendroIntL globalGridPts=0;
            std::vector<bool> accumulated;
            accumulated.resize(mesh->getDegOfFreedom(),false);

            for(unsigned int elem = eleLocalBegin; elem < eleLocalEnd; elem++)
            {
                for(unsigned int k=0; k<(eleOrder+1); k++)
                    for(unsigned int j=0; j<(eleOrder+1); j++ )
                        for(unsigned int i=0; i<(eleOrder+1); i++)
                        {
                            const unsigned int nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            if( (nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd) && !(accumulated[nodeLookUp_CG]) )
                            {
                                const unsigned int nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                unsigned int ownerID, ii_x,jj_y, kk_z;
                                
                                mesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                                const double len=(double)(1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                                
                                const double x_ot=pNodes[ownerID].getX()+ ii_x*(len/((double)eleOrder));
                                const double y_ot=pNodes[ownerID].getY()+ jj_y*(len/((double)eleOrder));
                                const double z_ot=pNodes[ownerID].getZ()+ kk_z*(len/((double)eleOrder));

                                const double x = GRIDX_TO_X(x_ot);
                                const double y = GRIDY_TO_Y(y_ot);
                                const double z = GRIDZ_TO_Z(z_ot);

                                Point grid_pt(x,y,z);
                                const double d1 = (grid_pt - bssn::BSSN_BH_LOC[0]).abs();
                                const double d2 = (grid_pt - bssn::BSSN_BH_LOC[1]).abs();

                                if( (d1 < bssn::BSSN_BH1_CONSTRAINT_R ||  d2 < bssn::BSSN_BH2_CONSTRAINT_R) )
                                    accumulated[nodeLookUp_CG]=true;
                                else
                                {
                                    l2+=(constraintVec[nodeLookUp_CG]*constraintVec[nodeLookUp_CG]);
                                    accumulated[nodeLookUp_CG]=true;
                                    localGridPts++;
                                }
                            }
                        }
            }
            
            // now we use the AMR_R<1,2> to compute the constraint norm. 
            // for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
            // {
            //     if(maskVector[node]<maskthreshoold)
            //         continue;
                
            //      l2+=(constraintVec[node]*constraintVec[node]);
            //      localGridPts++;
             
            // }
            
            par::Mpi_Reduce(&l2,&l2_g,1,MPI_SUM,0,mesh->getMPICommunicator());
            par::Mpi_Reduce(&localGridPts,&globalGridPts,1,MPI_SUM,0,mesh->getMPICommunicator());


            if(!(mesh->getMPIRank()))
                l2_g=l2_g/(double)(globalGridPts);
            
        }
        
        return sqrt(l2_g);
        
        
    }

    template<typename T>
    double extractConstraints(const ot::Mesh* mesh,const T** constraintVar,const T* maskVec, double maskthreshoold,unsigned int timestep,double stime)
    {
        const unsigned int numConstraints=4;
        double constraintMaskedL2[numConstraints]; // remove the psi4

        unsigned int rankGlobal=mesh->getMPIRankGlobal();
        unsigned int npesGlobal=mesh->getMPICommSizeGlobal();
        MPI_Comm commGlobal=mesh->getMPIGlobalCommunicator();

        if(mesh->isActive())
        {

            unsigned int rankActive=mesh->getMPIRank();
            unsigned int npesActive=mesh->getMPICommSize();
            MPI_Comm commActive=mesh->getMPICommunicator();


            for(unsigned int index=0;index<numConstraints;index++)
            {
                constraintMaskedL2[index]=bssn::computeConstraintL2Norm(mesh,constraintVar[index],maskVec,maskthreshoold);
                if(!rankActive)
                {
                    std::cout<<YLW<<"\tConstraint " <<bssn::BSSN_CONSTRAINT_VAR_NAMES[index]<< " L2 : ("<<constraintMaskedL2[index]<<" )"<<NRM<<std::endl;
                }

            }

            if(!rankActive)
            {

                std::ofstream fileGW;
                char fName[256];
                sprintf(fName,"%s_Constraints.dat",bssn::BSSN_PROFILE_FILE_PREFIX.c_str());
                fileGW.open (fName,std::ofstream::app);
                // writes the header
                if(timestep==0)
                    fileGW<<"TimeStep\t"<<" time\t"<<" C_HAM\t"<<" C_MOM0\t"<<" C_MOM1\t"<<" C_MOM2\t"<<std::endl;

                fileGW<<timestep<<"\t"<<stime<<"\t"<<constraintMaskedL2[0]<<"\t"<<constraintMaskedL2[1]<<"\t"<<constraintMaskedL2[2]<<"\t"<<constraintMaskedL2[3]<<std::endl;
                fileGW.close();


            }



        }


        #if 0
        if(!rankGlobal)
        {
                for(unsigned int index=0;index<numConstraints;index++)
                {
                    if(constraintMaskedL2[index]>0.01)
                    {
                        if(bssn::KO_DISS_SIGMA>0.06)
                            bssn::KO_DISS_SIGMA=0.05;
                        else
                            bssn::KO_DISS_SIGMA=0.10;

                        break;
                    }
                }

        }
        par::Mpi_Bcast(&bssn::KO_DISS_SIGMA,1,0,commGlobal);
        #endif

        return 0.0;

    }


    


} // end of namespace bssn

