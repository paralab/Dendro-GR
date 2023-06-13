//
// Created by milinda on 1/15/17.
//

/**
 *
 * @author Milinda Fernando
 * @breif contains the utilities for tensor kronecker products for interpolations.
 *
 * */

#ifndef SFCSORTBENCH_DENDROTENSOR_H
#define SFCSORTBENCH_DENDROTENSOR_H


/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_AIIX_APPLY_ELEM (const int M, const double*  A, const double*  X, double*  Y);



/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_IIAX_APPLY_ELEM(const int M, const double*  A, const double*  X, double*  Y);


/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_IAIX_APPLY_ELEM (const int M, const double*  A, const double*  X, double*  Y);




/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_IAX_APPLY_ELEM_2D(const int M, const double*  A, const double*  X, double*  Y);



/** Apply the 1D interpolation for the input vector x and output the interpolated values in the vector Y.
 *
 *
 *
 * \param [in]  M  size of the vector
 * \param [in]  A  interpolation matrix
 * \param [in]  X  input data for the interpolation
 * \param [out] Y  interpolated values.
 *

 */
void DENDRO_TENSOR_AIX_APPLY_ELEM_2D (const int M, const double*  A, const double*  X, double*  Y);


#endif //SFCSORTBENCH_DENDROTENSOR_H
