#ifndef _HAVE_GEOM_H
#define _HAVE_GEOM_H

#include <math.h>

#include "fluidUtils.h"
#include "mpi.h"
#include "parameters.h"

// This file is a holdover from the lanl_rmhd project. I'm sure Dendro probably
// has a differential geometry library somewhere, and it would probably better
// to use that, but I just want to get the code working in its original
// flat-space form first and fix it later. - Jacob Fields, 11/15/2018

/**
 * @brief compute the dot product w.r.t the metirc of a vector
 * @param[in] vu : vector components (3-vector)
 * @param[in] gd : metric
 *
 * */
// square_vector {{{
inline double square_vector(double vu[3], double gd[3][3]) {
    return gd[0][0] * vu[0] * vu[0] + 2.0 * gd[0][1] * vu[0] * vu[1] +
           2.0 * gd[0][2] * vu[0] * vu[2] + gd[1][1] * vu[1] * vu[1] +
           2.0 * gd[1][2] * vu[1] * vu[2] + gd[2][2] * vu[2] * vu[2];
}
// }}}

/**
 * @brief compute the dot product w.r.t the metirc of a co-vector
 * @param[in] vu : co-vector components (3-vector)
 * @param[in] gu : inverse metric
 */
// square_form {{{
inline double square_form(double vd[3], double gu[3][3]) {
    return gu[0][0] * vd[0] * vd[0] + 2.0 * gu[0][1] * vd[0] * vd[1] +
           2.0 * gu[0][2] * vd[0] * vd[2] + gu[1][1] * vd[1] * vd[1] +
           2.0 * gu[1][2] * vd[1] * vd[2] + gu[2][2] * vd[2] * vd[2];
}
// }}}

// vector_lower_index {{{
inline void vector_lower_index(double vd[3], double vu[3], double gd[3][3]) {
    vd[0] = gd[0][0] * vu[0] + gd[0][1] * vu[1] + gd[0][2] * vu[2];
    vd[1] = gd[0][1] * vu[0] + gd[1][1] * vu[1] + gd[1][2] * vu[2];
    vd[2] = gd[0][2] * vu[0] + gd[1][2] * vu[1] + gd[2][2] * vu[2];
}
// }}}

// form_raise_index {{{
inline void form_raise_index(double vu[3], double vd[3], double gu[3][3]) {
    vu[0] = gu[0][0] * vd[0] + gu[0][1] * vd[1] + gu[0][2] * vd[2];
    vu[1] = gu[0][1] * vd[0] + gu[1][1] * vd[1] + gu[1][2] * vd[2];
    vu[2] = gu[0][2] * vd[0] + gu[1][2] * vd[1] + gu[2][2] * vd[2];
}
// }}}

// metric_vars {{{
inline void metric_vars(double gu[3][3], double *detg, double gd[3][3],
                        const double *pos) {
    if (fluid::FLUID_COORDS == 1) {
        // Cylindrical coordinates
        *detg                = pos[0] * pos[0];
        const double epsilon = 1.0e-15;

        gd[0][0]             = 1.0;
        gd[0][1]             = 0.0;
        gd[0][2]             = 0.0;
        gd[1][0]             = 0.0;
        gd[1][1]             = 1.0;
        gd[1][2]             = 0.0;
        gd[2][0]             = 0.0;
        gd[2][1]             = 0.0;
        gd[2][2]             = pos[0] * pos[0];

        gu[0][0]             = 1.0;
        gu[0][1]             = 0.0;
        gu[0][2]             = 0.0;
        gu[1][0]             = 0.0;
        gu[1][1]             = 1.0;
        gu[1][2]             = 0.0;
        gu[2][0]             = 0.0;
        gu[2][1]             = 0.0;
        gu[2][2]             = 1.0 / (epsilon + pos[0] * pos[0]);
    } else if (fluid::FLUID_COORDS == 0) {
        // Cartesian coordinates
        *detg    = 1.0;

        gd[0][0] = 1.0;
        gd[0][1] = 0.0;
        gd[0][2] = 0.0;
        gd[1][0] = 0.0;
        gd[1][1] = 1.0;
        gd[1][2] = 0.0;
        gd[2][0] = 0.0;
        gd[2][1] = 0.0;
        gd[2][2] = 1.0;

        gu[0][0] = 1.0;
        gu[0][1] = 0.0;
        gu[0][2] = 0.0;
        gu[1][0] = 0.0;
        gu[1][1] = 1.0;
        gu[1][2] = 0.0;
        gu[2][0] = 0.0;
        gu[2][1] = 0.0;
        gu[2][2] = 1.0;
    } else {
        printf(" Alternative coordinate system not supported.\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
}
// }}}

// metric {{{
inline void metric(double gd[3][3], double *detg, const double *pos) {
    if (fluid::FLUID_COORDS == 1) {
        // Cylindrical coordinates

        *detg    = pos[0] * pos[0];

        gd[0][0] = 1.0;
        gd[0][1] = 0.0;
        gd[0][2] = 0.0;
        gd[1][0] = 0.0;
        gd[1][1] = 1.0;
        gd[1][2] = 0.0;
        gd[2][0] = 0.0;
        gd[2][1] = 0.0;
        gd[2][2] = pos[0] * pos[0];
    } else if (fluid::FLUID_COORDS == 0) {
        // Cartesian coordinates
        *detg    = 1.0;

        gd[0][0] = 1.0;
        gd[0][1] = 0.0;
        gd[0][2] = 0.0;
        gd[1][0] = 0.0;
        gd[1][1] = 1.0;
        gd[1][2] = 0.0;
        gd[2][0] = 0.0;
        gd[2][1] = 0.0;
        gd[2][2] = 1.0;
    } else {
        printf(" Alternative coordinate system not supported\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
}
// }}}
#endif
