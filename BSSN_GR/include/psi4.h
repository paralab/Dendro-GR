#ifndef PHYSCON_H
#define PHYSCON_H

#include <iostream>

#include "derivs.h"
#include "grUtils.h"
#include "parameters.h"

#define deriv_x deriv42_x
#define deriv_y deriv42_y
#define deriv_z deriv42_z

#define deriv_xx deriv42_xx
#define deriv_yy deriv42_yy
#define deriv_zz deriv42_zz

void psi4(double **uZipConVars, const double **uZipVars,
          const unsigned int &offset, const double *pmin, const double *pmax,
          const unsigned int *sz, const unsigned int &bflag);

#endif
