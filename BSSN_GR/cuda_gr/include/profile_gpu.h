//
// Created by milinda on 8/30/18.
//

/**
 * @author Milinda Fernando
 * School Of Computing, University of Utah
 * @brief contains profile parameters for the derv computations.
 *
 * */

#ifndef DENDRO_5_0_PROFILE_GPU_H
#define DENDRO_5_0_PROFILE_GPU_H

#include <iostream>
#include <vector>

#include "block.h"
#include "mesh.h"
#include "profiler.h"

namespace cuda {
namespace profile {

/** overall gpu exuxution time*/
extern profiler_t t_overall;

/**host to device communication*/
extern profiler_t t_H2D_Comm;
/**device to host communication*/
extern profiler_t t_D2H_Comm;

/**memory allocation deallocation time*/
extern profiler_t t_cudaMalloc_derivs;

/**total rhs computation time*/
extern profiler_t t_rhs_gpu;
extern profiler_t t_rhs_cpu;
extern profiler_t t_rhs_total;

/**initialize the profile counters*/
void initialize();

/**output the timings*/
void printOutput(const ot::Mesh* pMesh);

/**computes overall stats accross the processors*/
template <typename T>
void computeOverallStats(T* stat, T* stat_g, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    par::Mpi_Reduce(stat, stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 2, 1, MPI_MAX, 0, comm);
    stat_g[1] /= (npes);
}

void printOutput(const std::vector<ot::Block>& localBlkList);

}  // end of namespace profile

}  // end of namespace cuda

#endif  // DENDRO_5_0_PROFILE_GPU_H
