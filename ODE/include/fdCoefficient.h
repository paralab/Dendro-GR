//
// Created by milinda on 1/30/17.
//

/**
 *
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief This file contains coefficients for finite difference computations.
 *
 * */

#ifndef SFCSORTBENCH_FDCOEFFICIENT_H
#define SFCSORTBENCH_FDCOEFFICIENT_H

namespace fd {
//centered differencing,

/** first derivative coefficients for central differencing with order 2 accuracy*/
    static double D1_ORDER_2_CENTERED[3] = {0.5, 0, 0.5};

/** first derivative coefficients for central differencing with order 4 accuracy*/
    static double D1_ORDER_4_CENTERED[5] = {1.0 / 12.0, -2.0 / 3.0, 0, 2.0 / 3.0, -1.0 / 12.0};

/** first derivative coefficients for central differencing with order 6 accuracy*/
    static double D1_ORDER_6_CENTERED[7] = {-1.0 / 60.0, 3.0 / 20.0, -3.0 / 4.0, 0, 3.0 / 4.0, -3.0 / 20.0, 1.0 / 60.0};

/** first derivative coefficients for central differencing with order 8 accuracy*/
    static double D1_ORDER_8_CENTERED[9] = {1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0, -4.0 / 5.0, 0, 4.0 / 5.0, -1.0 / 5.0,
                                            4.0 / 105.0, -1 / 280.0};


/** second derivative coefficients for central differencing with order 4 accuracy*/
    static double D2_ORDER_4_CENTERED[5] = {-1.0 / 12.0, 4.0 / 3.0, -5.0 / 2.0, 4.0 / 3.0, -1.0 / 12.0};



// forward differencing
/** first derivative coefficients for forward differencing with order 4 accuracy*/
    static double D1_ORDER_4_FORWARD[5] = {-25.0 / 12.0, 4.0, -3.0, 4.0 / 3.0, -1.0 / 4.0};
/** first derivative coefficients for advective derivative in the upwind direction. */
    static double D1_ORDER_4_UPWIND[5]={-3.0/12.0,-10.0/12.0,18.0/12.0,-6.0/12.0,1.0/12.0};
/** second derivative coefficients for forward differencing with order 4 accuracy*/
    static double D2_ORDER_4_FORWARD[6] = {15.0 / 4.0, -77.0 / 6.0, 107.0 / 6.0, -13.0, 61.0 / 12.0, -5.0 / 6.0};


//backward differencing
/** first derivative coefficients for backward differencing with order 4 accuracy*/
    static double D1_ORDER_4_BACKWARD[5] = {1.0 / 4.0, -4.0 / 3.0 , 3.0 , -4.0, 25.0 / 12.0};
/**first derivative coefficients for advective derivative in the downwind direction*/
    static  double D1_ORDER_4_DOWNWIND[5]={-1.0/12.0,6.0/12.0,-18.0/12.0,10.0/12.0,3.0/12.0};

/** second derivative coefficients for backward differencing with order 4 accuracy*/
    static double D2_ORDER_4_BACKWARD[6] = {15.0 / 4.0, -77.0 / 6.0, 107.0 / 6.0, -13.0, 61.0 / 12.0, -5.0 / 6.0};

} // end namespace fd

#endif //SFCSORTBENCH_FDCOEFFICIENT_H
