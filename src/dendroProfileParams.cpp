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

namespace dendro
{

    namespace  timer
    {

        namespace sfcmatvec
        {
            profiler_t t_computeIJK;
            profiler_t t_parent_bucket;
            profiler_t t_pts_bucket;

                     profiler_t t_pts_p1_count;
                     profiler_t t_pts_p2_count;
                     profiler_t t_pts_p1_cpy;
                     profiler_t t_pts_p2_cpy;
                     profiler_t t_pts_p1_accum;
                     profiler_t t_pts_p2_accum;

            profiler_t t_p2cInterp;
            profiler_t t_internCpy;
            profiler_t t_elemMvec;
            profiler_t t_c2pInterp;
            profiler_t t_accum;
            profiler_t t_malloc;


            void initializeSFCMatvecTimers()
            {


                t_computeIJK.clear();
                t_parent_bucket.clear();
                t_pts_bucket.clear();
                t_pts_p1_count.clear();
                t_pts_p2_count.clear();

                t_pts_p1_cpy.clear();
                t_pts_p2_cpy.clear();

                t_pts_p1_accum.clear();
                t_pts_p1_accum.clear();

                t_p2cInterp.clear();
                t_internCpy.clear();
                t_elemMvec.clear();
                t_c2pInterp.clear();
                t_accum.clear();
                t_malloc.clear();


                t_computeIJK.start();
                t_parent_bucket.start();
                t_pts_bucket.start();

                t_pts_p1_count.start();
                t_pts_p2_count.start();

                t_pts_p1_cpy.start();
                t_pts_p2_cpy.start();

                t_pts_p1_accum.start();
                t_pts_p2_accum.start();

                t_p2cInterp.start();
                t_internCpy.start();
                t_elemMvec.start();
                t_c2pInterp.start();
                t_accum.start();
                t_malloc.start();
            }


            void consolePrint()
            {

                //std::cout<<" compute ijk: "<<t_computeIJK.seconds<<std::endl;
                std::cout<<" t_parent_bucket: "<<t_parent_bucket.seconds<<std::endl;
                std::cout<<" t_pts_bucket(leaf): "<<t_pts_bucket.seconds<<std::endl;

                std::cout<<" t_pts_p1_count: "<<t_pts_p1_count.seconds<<std::endl;
                std::cout<<" t_pts_p2_count: "<<t_pts_p2_count.seconds<<std::endl;
                std::cout<<" t_pts_p1_cpy: "<<t_pts_p1_cpy.seconds<<std::endl;
                std::cout<<" t_pts_p2_cpy: "<<t_pts_p2_cpy.seconds<<std::endl;

                std::cout<<" t_p2cInterp: "<<t_p2cInterp.seconds<<std::endl;
                //std::cout<<" t_internCpy: "<<t_internCpy.seconds<<std::endl;
                std::cout<<" t_c2pInterp: "<<t_c2pInterp.seconds<<std::endl;
                std::cout<<" t_elemMvec: "<<t_elemMvec.seconds<<std::endl;
                std::cout<<" t_accum (leaf): "<<t_accum.seconds<<std::endl;
                std::cout<<" t_pts_p1_accum: "<<t_pts_p1_accum.seconds<<std::endl;
                std::cout<<" t_pts_p2_accum: "<<t_pts_p2_accum.seconds<<std::endl;
                std::cout<<" t_malloc: "<<t_malloc.seconds<<std::endl;


            }


        } // end of namespace sfcMatvec

    }// end of namespace timer

}// end of namespace dendro