//
// Created by milinda on 3/31/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains useful math helper routines
*/
//

#ifndef SFCSORTBENCH_MATHUTILS_H
#define SFCSORTBENCH_MATHUTILS_H

#include <iostream>
#include <cmath>
#include <math.h>
#include "parUtils.h"

/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm between two vectors.
 * */
template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n,MPI_Comm comm);

/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l1 norm between two vectors.
 * */
template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n, MPI_Comm comm);


/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @return l1 norm between two vectors.
 * */
template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n);



/**
 * @brief computes the l2 norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm of vec.
 * */
template <typename T>
T normL2(T * vec,unsigned int n, MPI_Comm comm);


/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @return l2 norm between two vectors.
 * */
template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n);

/**
 * @brief computes the l2 norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return l2 norm of vec.
 * */
template <typename T>
T normL2(T * vec,unsigned int n);


/**
 * @brief computes the l_inf norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return l2 norm of vec.
 * */
template <typename T>
T normLInfty(T * vec,unsigned int n);


/**
 * @brief computes the l_inf norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm of vec.
 * */
template <typename T>
T normLInfty(T * vec,unsigned int n,MPI_Comm comm);

/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return min of vec.
 * */
template <typename T>
T vecMin(T * vec,unsigned int n);


/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return max of vec.
 * */
template <typename T>
T vecMax(T * vec,unsigned int n);


/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return min of vec.
 * */
template <typename T>
T vecMin(T * vec,unsigned int n,MPI_Comm comm);


/**
 * @brief computes the min of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return max of vec.
 * */
template <typename T>
T vecMax(T * vec,unsigned int n,MPI_Comm comm);

/**@brief : performs the dot product of any two given vectors
 * @param[in] v1: input vector 1
 * @param[in] v2: input vector 2
 * @param[in] n: size of the vector ( dof of vector)
 * @return  v1^tv2
 * */
template <typename T>
T dot(const T* v1, const T*v2,const unsigned int n);

/**@brief : performs the dot product of any two given distributed vectors
 * Assumption : v1 and v2 are partitioned in the sameway.
 * @param[in] v1: input vector 1
 * @param[in] v2: input vector 2
 * @param[in] n: size of the vector ( dof of vector)
 * @return  v1^tv2
 * */
template <typename T>
T dot(const T* v1, const T*v2,const unsigned int n,MPI_Comm comm);


/**
 * @brief: Scaler multiplication of a vector.
 * @param[in] alpha: scalar value
 * @param[in] v: input vector.
 * @param[in] n: size of the vector ( dof of vector)
 * @param[in] comm : MPI communicator.
 *
 * */
template<typename T>
void mul(const T alpha, const T* v, const unsigned int n, T* out);


/**
 * @brief : add two vectors.
 * @param[in] v1: input vector 1
 * @param[in] v2: input vector 2
 * @param[in] n: size of the vector ( dof of vector)
 * @param[out] out: out=v1+v2
 * */
template <typename T>
T add(const T* v1, const T*v2,const unsigned int n, T* out);

/**
* @brief : substract two vectors.
* @param[in] v1: input vector 1
* @param[in] v2: input vector 2
* @param[in] n: size of the vector ( dof of vector)
* @param[out] out: out=v1+v2
                          * */
template <typename T>
T subt(const T* v1, const T*v2,const unsigned int n, T* out);


/**
 * @brief Kroneckor product of two matrices. 
 * @tparam T data type
 * @param M1 : matrix 1
 * @param M2 : matrix 2
 * @param out : result matrix
 * @param r1 : rows in m1
 * @param c1 : cols in m1
 * @param r2 : rows in m2
 * @param c2 : cols in m2
 */
template<typename T>
void kron(const T* M1, const T* M2, T* out, unsigned int r1,unsigned int c1, unsigned int r2, unsigned int c2);

/**@brief: min mean max of a scalar value. */
template<typename T>
void min_mean_max(T* stat, T* stat_g, MPI_Comm comm);



#include "mathUtils.tcc"


#endif //SFCSORTBENCH_MATHUTILS_H
