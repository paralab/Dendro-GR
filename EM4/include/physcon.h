#pragma once  

#include <iostream>
#include "parameters.h"
#include "em4Utils.h"
#include "derivs.h"

void physical_constraints( double **uZipConVars, const double **uZipVars,
                       const unsigned int& offset,
                       const double *pmin, const double *pmax,
                       const unsigned int *sz, const unsigned int& bflag);

