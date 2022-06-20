//
// Created by milinda on 10/21/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "profile_params.h"

namespace em3
{
    namespace timer
    {
        profiler_t 	total_runtime;

        profiler_t  t_f2o;
        profiler_t  t_cons;
        profiler_t  t_bal;
        profiler_t	t_mesh;

        profiler_t  t_rkSolve;

        profiler_t	t_ghostEx_sync;

        profiler_t	t_unzip_sync;

        profiler_t	t_unzip_async;

        profiler_t	t_deriv;
        profiler_t	t_rhs;

        profiler_t	t_rhs_E;
        profiler_t	t_rhs_B;
       // profiler_t	t_rhs_psi;

        profiler_t  t_bdyc;

        profiler_t	t_zip;
        profiler_t  t_rkStep;

        profiler_t  t_isReMesh;
        profiler_t	t_gridTransfer;
        profiler_t	t_ioVtu;
        profiler_t  t_ioCheckPoint;

    }
}
