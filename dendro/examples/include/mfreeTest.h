//
// Created by milinda on 1/9/19.
//

/**
 * @brief contains bechmark functions for mesh free vs mesh based FEM computation.
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * */

#ifndef DENDRO_5_0_MFREETEST_H
#define DENDRO_5_0_MFREETEST_H

#include "TreeNode.h"
#include "mpi.h"
#include "genPts_par.h"
#include "sfcSort.h"
#include "mesh.h"
#include "dendro.h"
#include "dendroIO.h"
#include "octUtils.h"
#include "functional"
#include "fdCoefficient.h"
#include "stencil.h"
#include "rkTransport.h"
#include "refel.h"
#include "operators.h"
#include "cg.h"
#include "matvec.h"



namespace mfree
{


    extern profiler_t tLev[31];

    /**@brief performs roofline bucketing.
     * @param [in] ptList: input point list
     * @param [in] ptSz: point size
     *
     * */
    void ptBucketRoofline(const Point* ptList, unsigned int ptSz,unsigned int lev,unsigned int pMaxDepth);
}



#endif //DENDRO_5_0_MFREETEST_H
