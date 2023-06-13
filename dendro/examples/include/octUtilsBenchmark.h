//
// Created by milinda on 8/8/16.

/**
 * @author Milinda Fernando
 * @breif Contains the code for benchmarking the code written for octree generations, balancing and remove duplictaes.
 * @brief Capable of running the weak scaling and generating the results table
 *
 *
 * */

//

#ifndef SFCSORTBENCH_OCTUTILSBENCHMARK_H
#define SFCSORTBENCH_OCTUTILSBENCHMARK_H

#define ROOT_ROT_ID 0

#include "genPts_par.h"
#include "TreeNode.h"
#include "dendro.h"
#include "parUtils.h"
#include "sfcSort.h"
#include "testUtils.h"
#include <mpi.h>

/*
 *
 * stat data structure to store weak scaling run statistics.
 *
 * */

struct RuntimeStat
{
    unsigned int npes;
    unsigned int numPts;
    double tolerance;
    unsigned int options;
    std::vector<double> stats;
    std::vector<double> stats_prev;

    RuntimeStat()
    {
        npes=0;
        numPts=0;
        tolerance=0;
        options=0;
        stats.clear();
        stats_prev.clear();
    }

};




/*
 * Executes the treeSort with the options =1 Which removes the duplicates from the given input and sort it.
 * */

void sfcTreeSortTest(DendroIntL numPts, unsigned int dim, unsigned int maxDepth, double tolerance, int distribution,
                     unsigned int k, unsigned int options, unsigned int sf_k,char *prefix, MPI_Comm pcomm);

void weakScalingDriver(DendroIntL numPts,unsigned int dim,unsigned int maxDepth,double tolerance,int distribution,unsigned int k,unsigned int options,unsigned int sf_k,char * prefix,MPI_Comm comm);





#endif //SFCSORTBENCH_OCTUTILSBENCHMARK_H
