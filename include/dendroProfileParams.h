//
// Created by milinda on 1/19/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Constains profile parameters for Dendro performance measurements
*/
//

#ifndef SFCSORTBENCH_PROF_PARAMS_H
#define SFCSORTBENCH_PROF_PARAMS_H

#include "profiler.h"

namespace dendro
{
    namespace timer
    {


        // profile parameters for stages of unzip (sync).

        extern profiler_t	t_unzip_sync_internal;
        extern profiler_t	t_unzip_sync_face[NUM_FACES];
        extern profiler_t	t_unzip_sync_edge;
        extern profiler_t	t_unzip_sync_vtex;
        extern profiler_t	t_unzip_p2c;
        extern profiler_t	t_unzip_sync_nodalval;
        extern profiler_t   t_unzip_sync_cpy;
        extern profiler_t   t_unzip_sync_f_c1;
        extern profiler_t   t_unzip_sync_f_c2;
        extern profiler_t   t_unzip_sync_f_c3;


        // unzip async
        extern profiler_t	t_unzip_async_internal;
        extern profiler_t	t_unzip_async_external;
        extern profiler_t	t_unzip_async_comm;



    } // end of namespace timer.

}//end namespace dendro

#endif //SFCSORTBENCH_PROF_PARAMS_H
