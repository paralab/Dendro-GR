//
// Created by milinda on 10/21/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief contains the profile parameters.
*/
//

#ifndef SFCSORTBENCH_PROFILE_PARAMS_H
#define SFCSORTBENCH_PROFILE_PARAMS_H

#include "profiler.h"


namespace nlsm
{
    namespace timer
    {
        extern profiler_t 	total_runtime;


        extern profiler_t 	t_f2o;
        extern profiler_t 	t_cons;
        extern profiler_t 	t_bal;
        extern profiler_t	t_mesh;


        extern profiler_t 	t_rkSolve;
        extern profiler_t	t_ghostEx_sync;

        extern profiler_t	t_unzip_sync;
        extern profiler_t	t_unzip_async;

        extern profiler_t	t_deriv;
        extern profiler_t	t_rhs;

        extern profiler_t	t_rhs_a;
        extern profiler_t	t_rhs_b;
        extern profiler_t	t_rhs_gt;
        extern profiler_t	t_rhs_chi;
        extern profiler_t	t_rhs_At;
        extern profiler_t	t_rhs_K;
        extern profiler_t	t_rhs_Gt;
        extern profiler_t	t_rhs_B;

        extern profiler_t   t_bdyc;

        extern profiler_t	t_zip;
        extern profiler_t 	t_rkStep;

        extern profiler_t t_isReMesh;
        extern profiler_t	t_gridTransfer;
        extern profiler_t	t_ioVtu;
        extern profiler_t t_ioCheckPoint;

    }
}

#endif //SFCSORTBENCH_PROFILE_PARAMS_H
