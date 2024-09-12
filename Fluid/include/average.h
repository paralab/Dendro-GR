#ifndef AVERAGE_H_
#define AVERAGE_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fluid.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

namespace average {
// The weights are not a trivial calculation, so they are precomputed and stored
// separately.
extern const double glweights[];
extern const double glcth[];
extern const int glsize;

// Integration
double qgaus(double data[], const double w[], const int n);

// 2d averaging with bilinear interpolation.
double avg2d(double *datain[], double *xi[], int shp[], double r,
             const double cth[], const double w[], const int n);

// Linear interpolation in 1 and 2 dimensions.
double lerp(double datain[], double xi[], int shp, double x, const int dir);

double bilerp(double *datain[], double *xi[], int shp[], double x[]);
}  // namespace average

#endif
