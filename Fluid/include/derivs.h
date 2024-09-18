#ifndef _DERVS_H
#define _DERVS_H

#include <cmath>

#include "TreeNode.h"

#define IDX(i, j, k) ((i) + nx * ((j) + ny * (k)))

inline void deriv_xm(double fp, double fm, double dx, double *dtu) {
    *dtu = -(fp - fm) / dx;
}

inline void deriv_yzm(double fp, double fm, double dx, double *dtu) {
    *dtu -= (fp - fm) / dx;
}

/**@brief: copies unzipped padding to the computed derivatives.
 * (this is essential for the mixed derivatives.)
 * */
void cpy_unzip_padd(double *const Du, const double *const u,
                    const unsigned int *sz, unsigned bflag);

#endif
