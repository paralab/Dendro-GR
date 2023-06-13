//
// Created by milinda on 9/8/16.
//

/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 *
 * @breif Contains a simple implementation to check the mesh functionalities.
 * */

#ifndef SFCSORTBENCH_MESHBENCHMARK_H
#define SFCSORTBENCH_MESHBENCHMARK_H


#define ROOT_ROT_ID 0



#include "TreeNode.h"
#include "mpi.h"
#include "genPts_par.h"
#include "sfcSort.h"
#include "mesh.h"
#include "dendro.h"
#include "dendroIO.h"
#include "dollar.hpp"


void weakScalingDriver(char * ptsFile,bool genPts,DendroIntL numPts,unsigned int dim,unsigned int maxDepth,double tolerance,int distribution,unsigned int k,unsigned int options,unsigned int sf_k,char * prefix,MPI_Comm comm);



void meshBenchMark(char * ptsFile,bool genPts,unsigned int numPts,unsigned int dim, unsigned int maxDepth,unsigned int distribution,double tol,unsigned int sf_k,unsigned int options,char * prefix, MPI_Comm comm);






#endif //SFCSORTBENCH_MESHBENCHMARK_H
