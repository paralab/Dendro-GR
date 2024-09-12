/**
 * @file fluidUtils.tcc
 * @brief Contains templated utility functions for Fluid
 * @version 1.0
 */
namespace fluid {

template <typename T>
double computeFluidWavelet(const ot::Mesh* pMesh,
                           const wavelet::WaveletEl* wRefEl, const T* eleUnzip,
                           double* maxima, double wavelet_tol, unsigned int dof,
                           bool isBdyEle) {
    using namespace wavelet;
    double wMax = 0.0;

    if (pMesh->isActive()) {
        RefElement* refEl = (RefElement*)pMesh->getReferenceElement();

        const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
        const unsigned int eOrder            = pMesh->getElementOrder();

        const unsigned int numLocalElements  = pMesh->getNumLocalMeshElements();

        std::vector<double> wCout;
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();

        const unsigned int nx      = (2 * eOrder + 1);
        const unsigned int ny      = (2 * eOrder + 1);
        const unsigned int nz      = (2 * eOrder + 1);

        // FIXME: This can be passed in as an argument and made more general.
        const unsigned int pw      = blkList[0].get1DPadWidth();
        assert(pw == (eOrder >> 1u));

        const unsigned int sz_per_dof = nx * ny * nz;
        const unsigned int isz[]      = {nx, ny, nz};
        wCout.resize(sz_per_dof);

        const double tol_ele = wavelet_tol;
        for (unsigned int v = 0; v < dof; v++) {
            ((WaveletEl*)wRefEl)
                ->compute_wavelets_3D((double*)(eleUnzip + v * sz_per_dof), isz,
                                      wCout, isBdyEle);
            const double l_max = (normL2(wCout.data(), wCout.size())) /
                                 sqrt(wCout.size()) / maxima[v];

            if (wMax < l_max) {
                wMax = l_max;
            }
            // For early bail out.
            if (wMax > tol_ele) {
                break;
            }
        }
    }

    return wMax;
}
// Fluid WAMR Refinement
template <typename T>
bool computeFluidRemeshFlags(
    const ot::Mesh* pMesh, std::vector<unsigned int>& refine_flags,
    const T** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double)> wavelet_tol,
    double amr_coarse_fac, bool includeBdy) {
    using namespace wavelet;
    bool isMeshGlobalChanged = false;
    bool isMeshLocalChanged  = false;

    if (pMesh->isActive()) {
        RefElement* refEl = (RefElement*)pMesh->getReferenceElement();
        WaveletEl* wrefEl = new WaveletEl(refEl);

        const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
        const unsigned int eOrder             = pMesh->getElementOrder();

        const unsigned int numLocalElements = pMesh->getNumLocalMeshElements();

        refine_flags.clear();
        refine_flags.resize(numLocalElements, OCT_NO_CHANGE);

        std::vector<T> blkIn;
        std::vector<double> wCout;
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();

        std::vector<double> eleWMax;
        eleWMax.resize(numLocalElements, 0);

        const unsigned int eleOfst = pMesh->getElementLocalBegin();

        // Calculate the maxima across all refinement variables.
        double* maxVal             = new double[numVars];
        for (unsigned int v = 0; v < numVars; v++) {
            const unsigned int vid = varIds[v];
            for (unsigned int blk = 0; blk < blkList.size(); blk++) {
                const unsigned int pw = blkList[blk].get1DPadWidth();
                unsigned int sz[3]    = {0};
                sz[0]                 = blkList[blk].getAllocationSzX();
                sz[1]                 = blkList[blk].getAllocationSzY();
                sz[2]                 = blkList[blk].getAllocationSzZ();
                unsigned int offset   = blkList[blk].getOffset();
                // Loop over every element, making sure to exclude the padding
                // regions.
                for (unsigned int i = pw; i < sz[0] - pw; i++) {
                    for (unsigned int j = pw; j < sz[1] - pw; j++) {
                        for (unsigned int k = pw; k < sz[2] - pw; k++) {
                            unsigned int pp = i + sz[0] * (j + k * sz[1]);
                            maxVal[v]       = fmax(
                                maxVal[v], fabs(unzippedVec[vid][offset + pp]));
                        }
                    }
                }
            }
            // Reduce across all processors.
            double global_max;
            MPI_Allreduce(&maxVal[v], &global_max, 1, MPI_DOUBLE, MPI_MAX,
                          pMesh->getMPICommunicator());
            maxVal[v] = global_max;
            // Prevent NaNS.
            if (maxVal[v] == 0) {
                maxVal[v] = 1e-15;
            }
        }

        for (unsigned int blk = 0; blk < blkList.size(); blk++) {
            const unsigned int pw = blkList[blk].get1DPadWidth();
            if ((eOrder >> 1u) != pw) {
                std::cout << " padding with should be half the eleOrder for "
                             "generic wavelet computations. "
                          << std::endl;
                MPI_Abort(pMesh->getMPICommunicator(), 0);
            }
            const unsigned int nx = (2 * eOrder + 1);
            const unsigned int ny = (2 * eOrder + 1);
            const unsigned int nz = (2 * eOrder + 1);

            blkIn.resize(numVars * nx * ny * nz);
            const unsigned int isz[] = {nx, ny, nz};
            const unsigned int bflag = blkList[blk].getBlkNodeFlag();

            for (unsigned int ele = blkList[blk].getLocalElementBegin();
                 ele < blkList[blk].getLocalElementEnd(); ele++) {
                const bool isBdyOct = pMesh->isBoundaryOctant(ele);
                const unsigned int szby2 =
                    1u << (m_uiMaxDepth - pNodes[ele].getLevel() - 1);
                Point oct_pt = Point(pNodes[ele].minX() + szby2,
                                     pNodes[ele].minY() + szby2,
                                     pNodes[ele].minZ() + szby2);
                Point domain_pt;
                pMesh->octCoordToDomainCoord(oct_pt, domain_pt);
                const double tol_ele =
                    wavelet_tol(domain_pt.x(), domain_pt.y(), domain_pt.z());

                if (!includeBdy & isBdyOct) {
                    // tol small enough to not refine but not to coarsen.
                    eleWMax[ele - eleOfst] = amr_coarse_fac * tol_ele + 1e-8;
                    continue;
                }

                for (unsigned int v = 0; v < numVars; v++) {
                    const unsigned int vid = varIds[v];
                    pMesh->getUnzipElementalNodalValues(
                        unzippedVec[vid], blk, ele,
                        blkIn.data() + v * (nx * ny * nz), true);
                }

                eleWMax[ele - eleOfst] =
                    computeFluidWavelet(pMesh, wrefEl, blkIn.data(), maxVal,
                                        tol_ele, numVars, isBdyOct);
            }
        }

        // Delete the wavelet reference element.
        delete wrefEl;
        delete[] maxVal;

        // Mark elements for refinement first.
        for (unsigned int ele = pMesh->getElementLocalBegin();
             ele < pMesh->getElementLocalEnd(); ele++) {
            const unsigned int szby2 =
                1u << (m_uiMaxDepth - pNodes[ele].getLevel() - 1);
            Point oct_pt =
                Point(pNodes[ele].minX() + szby2, pNodes[ele].minY() + szby2,
                      pNodes[ele].minZ() + szby2);

            Point domain_pt;
            pMesh->octCoordToDomainCoord(oct_pt, domain_pt);

            const double tol_ele =
                wavelet_tol(domain_pt.x(), domain_pt.y(), domain_pt.z());
            const double l_max = eleWMax[ele - eleOfst];

            if (l_max > tol_ele) {
                refine_flags[(ele - eleOfst)] = OCT_SPLIT;
                isMeshLocalChanged            = true;
            } else if (l_max < amr_coarse_fac * tol_ele) {
                refine_flags[(ele - eleOfst)] = OCT_COARSE;
                isMeshLocalChanged            = true;
            } else {
                refine_flags[ele - eleOfst] = OCT_NO_CHANGE;
            }
        }
        /*if(isMeshLocalChanged){
          isMeshLocalChanged = pMesh->setMeshRefinementFlags(refine_flags);
        }*/
    }

    MPI_Allreduce(&isMeshLocalChanged, &isMeshGlobalChanged, 1, MPI_CXX_BOOL,
                  MPI_LOR, pMesh->getMPIGlobalCommunicator());

    return isMeshGlobalChanged;
}
}  // namespace fluid
