//
// Created by milinda on 10/20/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief simple profiler based on Hari's sort_profiler for bssn application.
*/
//

#ifndef SFCSORTBENCH_DENDRO_PROFILER_H
#define SFCSORTBENCH_DENDRO_PROFILER_H

#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "dendro.h"
#ifdef HAVE_PAPI
#include <papi.h>
#endif

#include <omp.h>

class profiler_t
{
    public:
        profiler_t ();
        virtual ~profiler_t ();

        void start();
        void stop();
        void clear();
        void snapreset();

        public:
            long double	seconds;  // openmp wall time
            long long p_flpops; // papi floating point operations
            long double snap; // snap shot of the cumilative time.
            long long num_calls; // number of times the timer stop function called.

        private:
            void  flops_papi();

        protected:
            long double	  _pri_seconds;  // openmp wall time
            long long _pri_p_flpops; // papi floating point operations

};








#endif //SFCSORTBENCH_DENDRO_PROFILER_H
