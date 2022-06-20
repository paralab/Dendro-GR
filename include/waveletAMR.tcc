/**
 * @file waveletAMR.tcc
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Contains utility functions to perform wavelet AMR. 
 * @version 0.1
 * @date 2020-06-13
 * @copyright Copyright (c) 2020
 * 
 */


namespace wavelet
{


    template<typename T>
    double compute_element_wavelet(const ot::Mesh* pMesh, const WaveletEl* wRefEl, const T*eleUnzip, double wavelet_tol,unsigned int dof, bool isBdyEle)
    {  

        double wMax=0.0;

        if(pMesh->isActive())
        {
            RefElement* refEl = (RefElement*)pMesh->getReferenceElement();
            
            const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
            const unsigned int eOrder = pMesh->getElementOrder();

            const unsigned int numLocalElements = pMesh->getNumLocalMeshElements();
            
            
            std::vector<double> wCout;
            const ot::TreeNode* pNodes = pMesh->getAllElements().data();

            const unsigned int nx = (2*eOrder+1);
            const unsigned int ny = (2*eOrder+1);
            const unsigned int nz = (2*eOrder+1); 

            const unsigned int pw = blkList[0].get1DPadWidth(); // assumes const padd. with for all blocks. 
            assert(pw == (eOrder>>1u));

            const unsigned int sz_per_dof = nx*ny*nz;
            const unsigned int isz[] = {nx,ny,nz};
            wCout.resize(sz_per_dof);

            const double tol_ele = wavelet_tol;
            for(unsigned int v=0; v < dof; v++)
            {
                ((WaveletEl*)wRefEl)->compute_wavelets_3D((double*)(eleUnzip+ v*sz_per_dof),isz,wCout,isBdyEle);
                const double l_max = (normL2(wCout.data(),wCout.size()))/sqrt(wCout.size());

                if(wMax < l_max)
                    wMax = l_max;

                // for early bail out. 
                if(wMax > tol_ele)
                    break;

            }

        }


        return wMax;

    }


    template<typename T>
    bool compute_wavelet_remesh_flags(const ot::Mesh* pMesh, std::vector<unsigned int>& refine_flags, const T** unzippedVec, const unsigned int *varIds, const unsigned int numVars, std::function<double(double, double, double)> wavelet_tol, double amr_coarse_fac, bool includeBdy)
    {

        bool isMeshGlobalChanged = false;
        bool isMeshLocalChanged  = false;
        //std::cout<<"calling amr"<<std::endl;

        if(pMesh->isActive())
        {
            RefElement* refEl = (RefElement*)pMesh->getReferenceElement();
            WaveletEl * wrefEl = new WaveletEl(refEl);

            const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
            const unsigned int eOrder = pMesh->getElementOrder();
            
            const unsigned int numLocalElements = pMesh->getNumLocalMeshElements();
            
            refine_flags.clear();
            refine_flags.resize(numLocalElements,OCT_NO_CHANGE);
            
            std::vector<T> blkIn;
            std::vector<double> wCout;
            const ot::TreeNode* pNodes = pMesh->getAllElements().data();

            std::vector<double> eleWMax;
            eleWMax.resize(numLocalElements,0);

            const unsigned int eleOfst = pMesh->getElementLocalBegin();

            
            for(unsigned int blk=0; blk <blkList.size(); blk++)
            {   
                const unsigned int pw = blkList[blk].get1DPadWidth();
                if((eOrder>>1u) != pw)
                {
                    std::cout<<" padding width should be half the eleOrder for generic wavelet computations. "<<std::endl;
                    MPI_Abort(pMesh->getMPICommunicator(),0);
                }

                const unsigned int nx = (2*eOrder+1);
                const unsigned int ny = (2*eOrder+1);
                const unsigned int nz = (2*eOrder+1);

                //std::cout<<"nx "<<nx<<std::endl;
        
                blkIn.resize(numVars*nx*ny*nz);
                const unsigned int isz[] = {nx,ny,nz};
                const unsigned int bflag = blkList[blk].getBlkNodeFlag();

                for(unsigned int ele =blkList[blk].getLocalElementBegin(); ele < blkList[blk].getLocalElementEnd(); ele++)
                {

                    const bool isBdyOct = pMesh->isBoundaryOctant(ele);
                    const unsigned int szby2 = 1u<<(m_uiMaxDepth-pNodes[ele].getLevel()-1);
                    Point oct_pt = Point(pNodes[ele].minX() + szby2  , pNodes[ele].minY() + szby2, pNodes[ele].minZ() + szby2);
                    Point domain_pt;
                    pMesh->octCoordToDomainCoord(oct_pt,domain_pt);
                    const double tol_ele = wavelet_tol(domain_pt.x(),domain_pt.y(),domain_pt.z());

                    if(!includeBdy && isBdyOct)
                    {
                        // tol small enough to not refine but not to coarsen . 
                        eleWMax[ele - eleOfst] = amr_coarse_fac*tol_ele + 1e-8; 
                        continue;
                    }
                        

                    for(unsigned int v=0; v < numVars; v++)
                    {
                        const unsigned int vid = varIds[v];
                        pMesh->getUnzipElementalNodalValues(unzippedVec[vid],blk, ele, blkIn.data() + v*(nx*ny*nz), true);
                    }

                    eleWMax[ele - eleOfst]=compute_element_wavelet(pMesh,wrefEl,blkIn.data(),tol_ele,numVars,isBdyOct);
                    
                    // if(isBdyOct)
                    //std::cout<<"ele :  "<<ele<<" eleWMax: "<<eleWMax[ele-eleOfst]<<std::endl;
                    
                }

            }

            // delete the wavelet reference element. 
            delete wrefEl;


            // mark elements for refinement first. 
            for(unsigned int ele = pMesh->getElementLocalBegin(); ele < pMesh->getElementLocalEnd(); ele++)
            {
                
                const unsigned int szby2 = 1u<<(m_uiMaxDepth-pNodes[ele].getLevel()-1);
                Point oct_pt = Point(pNodes[ele].minX() + szby2  , pNodes[ele].minY() + szby2, pNodes[ele].minZ() + szby2);
                
                Point domain_pt;
                pMesh->octCoordToDomainCoord(oct_pt,domain_pt);

                const double tol_ele = wavelet_tol(domain_pt.x(),domain_pt.y(),domain_pt.z());
                const double l_max = eleWMax[ele - eleOfst];

                if(l_max > tol_ele )
                {
                    refine_flags[(ele-eleOfst)] = OCT_SPLIT;
                    isMeshLocalChanged = true;

                }else if( l_max < amr_coarse_fac *tol_ele)
                {

                    refine_flags[ele-eleOfst] = OCT_COARSE;
                    isMeshLocalChanged=true;

                }else
                {
                    refine_flags[ele-eleOfst] = OCT_NO_CHANGE;
                } 

            }
                
            
        }

        
        //par::Mpi_Allreduce(&isMeshLocalChanged,&isMeshGlobalChanged,1,MPI_LOR,pMesh->getMPIGlobalCommunicator());
        MPI_Allreduce(&isMeshLocalChanged,&isMeshGlobalChanged,1,MPI_CXX_BOOL,MPI_LOR,pMesh->getMPIGlobalCommunicator());
        return isMeshGlobalChanged;
        
    }


    
    

}// end of namespace wavelet. 


