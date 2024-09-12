#pragma once

#include <iostream>

#include "derivs.h"
#include "em3Utils.h"
#include "parameters.h"

void physical_constraints(double **uZipConVars, const double **uZipVars,
                          const unsigned int &offset, const double *pmin,
                          const double *pmax, const unsigned int *sz,
                          const unsigned int &bflag);
