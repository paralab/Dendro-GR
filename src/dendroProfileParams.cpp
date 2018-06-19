//
// Created by milinda on 1/19/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//
#include "dendroProfileParams.h"

namespace dendro
{
    namespace timer
    {
        profiler_t	t_unzip_sync_internal;
        profiler_t	t_unzip_sync_face[NUM_FACES];
        profiler_t	t_unzip_sync_edge;
        profiler_t	t_unzip_sync_vtex;
        profiler_t	t_unzip_p2c;
        profiler_t	t_unzip_sync_nodalval;
        profiler_t  t_unzip_sync_cpy;
        profiler_t  t_unzip_sync_f_c1;
        profiler_t  t_unzip_sync_f_c2;
        profiler_t  t_unzip_sync_f_c3;

        profiler_t	t_unzip_async_internal;
        profiler_t	t_unzip_async_external;
        profiler_t	t_unzip_async_comm;

    } // end of namespace timer

} // end of namespace dendro

