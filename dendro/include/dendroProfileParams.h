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




namespace dendro
{

    namespace  timer
    {

        namespace sfcmatvec
        {

            extern profiler_t t_computeIJK;
            extern profiler_t t_parent_bucket;
            extern profiler_t t_pts_bucket;
                extern profiler_t t_pts_p1_count;
                extern profiler_t t_pts_p2_count;
                extern profiler_t t_pts_p1_cpy;
                extern profiler_t t_pts_p2_cpy;
                extern profiler_t t_pts_p1_accum;
                extern profiler_t t_pts_p2_accum;




            extern profiler_t t_p2cInterp;
            extern profiler_t t_internCpy;
            extern profiler_t t_elemMvec;
            extern profiler_t t_c2pInterp;
            extern profiler_t t_accum;
            extern profiler_t t_malloc;


            void initializeSFCMatvecTimers();

            void consolePrint();


        } // end of namespace sfcMatvec

    }// end of namespace timer

}// end of namespace dendro



#endif //SFCSORTBENCH_PROF_PARAMS_H
