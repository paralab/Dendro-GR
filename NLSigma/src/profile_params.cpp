//
// Created by milinda on 10/21/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "profile_params.h"

namespace nlsm
{
    namespace timer
    {



        profiler_t 	total_runtime;


        profiler_t t_f2o;
        profiler_t t_cons;
        profiler_t t_bal;
        profiler_t	t_mesh;


        profiler_t t_rkSolve;

        profiler_t	t_ghostEx_sync;

        profiler_t	t_unzip_sync;

        profiler_t	t_unzip_async;

        profiler_t	t_deriv;
        profiler_t	t_rhs;

        profiler_t	t_rhs_a;
        profiler_t	t_rhs_b;
        profiler_t	t_rhs_gt;
        profiler_t	t_rhs_chi;
        profiler_t	t_rhs_At;
        profiler_t	t_rhs_K;
        profiler_t	t_rhs_Gt;
        profiler_t	t_rhs_B;

        profiler_t t_bdyc;

        profiler_t	t_zip;
        profiler_t t_rkStep;

        profiler_t t_isReMesh;
        profiler_t	t_gridTransfer;
        profiler_t	t_ioVtu;
        profiler_t t_ioCheckPoint;

    }
}
