//
// Created by hari on 6/2/17.
//

/**
 * @file        simd_derivs.h
 * @author      hari sundar hari@cs.utah.edu
 *
 * simple wrapper for avx2 and avx512
 *
 */

#ifndef SFCSORTBENCH_SIMD_DERIVS_H
#define SFCSORTBENCH_SIMD_DERIVS_H

// #if defined(__AVX2__) || defined(__AVX512F__)
// #include <immintrin.h> // AVX
// #include <zmmintrin.h> // AVX512
// #endif

#include <x86intrin.h>

#ifdef __AVX512F__
#define SIMD_LENGTH 8
#define SIMD_ALIGNMENT 32

// load-store


// arithmetic

#endif

#ifdef __AVX2__
#define SIMD_LENGTH 4
#define SIMD_ALIGNMENT 32

// load-store


// arithmetic

#endif

// allocation
int simd_alloc(void **memptr, size_t size) {
  return posix_memalign(memptr, SIMD_ALIGNMENT, size);
}

#endif //SFCSORTBENCH_SIMD_DERIVS_H
