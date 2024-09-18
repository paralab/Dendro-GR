#include "derivs.h"

#include <cmath>
#include <iostream>

void cpy_unzip_padd(double *const Du, const double *const u,
                    const unsigned int *sz, unsigned bflag) {
    const int nx = sz[0];
    const int ny = sz[1];
    const int nz = sz[2];

    for (unsigned int k = 0; k < sz[2]; k++)
        for (unsigned int j = 0; j < sz[1]; j++)
            for (unsigned int i = 0; i < sz[0]; i++)
                if ((i < 3 || i >= sz[0] - 3) || (j < 3 || j >= sz[0] - 3) ||
                    (k < 3 || k >= sz[2] - 3))
                    Du[IDX(i, j, k)] = u[IDX(i, j, k)];
}
