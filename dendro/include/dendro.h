/**
 * @file dendro.h
 * @brief Basic dendro data types and definitions. 
 * @version 0.1
 * @date 2016-02-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include <climits>
#include <complex>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>



#define RED "\e[1;31m"
#define BLU "\e[2;34m"
#define GRN "\e[0;32m"
#define YLW "\e[0;33m"
#define MAG "\e[0;35m"
#define CYN "\e[0;36m"
#define NRM "\e[0m"



#ifdef USE_64BIT_INDICES
#define DendroIntL long long
#define DendroIntLSpecifier %lld
#define DendroUIntLSpecifier %llu
#else
#define DendroIntL unsigned int
#define DendroIntLSpecifier %d
#define DendroUIntLSpecifier %u
#endif


#define DendroScalar double
#define DendroComplex std::complex<double>



//#define DendroIntL unsigned int
typedef unsigned __int128 DendroUInt_128;

#define DEVICE_REAL double
#define DEVICE_INT int
#define DEVICE_UINT unsigned int
#define DEVICE_FLOAT32 float
#define DEVICE_BOOL bool

// mesh.h # defines.
#define LOOK_UP_TABLE_DEFAULT UINT_MAX


// TreeNode.h # defines.

/**
 * Following are the flags that used to mark treeNodes when generating keys and mesh generation.
 * all these flags are stored in m_uiLevel since we know that, we need only 4 bits to store the level of the octree.
 *
 * bit -0 to bit -4 : level of the octant
 * bit -5 to 13 : are the Key flags.
 * bit -14 to 15: index of the
 *
 *
 *
 * */
#define OCT_FOUND 32
#define OCT_KEY_NONE 64
#define OCT_KEY_SPLITTER  128
#define OCT_KEY_UP  256
#define OCT_KEY_DOWN  512
#define OCT_KEY_FRONT  1024
#define OCT_KEY_BACK  2048
#define OCT_KEY_LEFT  4096
#define OCT_KEY_RIGHT  8192

/**
 *
 * Note: Don't change that below numbering for key direction. With the following numbering you can always get the opposite direction performing XOR with 1.
 * Example OCT_DIR_LEFT= 1 XOR OCT_DIR_RIGHT
 *
 * */

#define OCT_DIR_LEFT 0
#define OCT_DIR_RIGHT 1
#define OCT_DIR_DOWN 2
#define OCT_DIR_UP 3
#define OCT_DIR_BACK 4
#define OCT_DIR_FRONT 5

#define OCT_DIR_LEFT_DOWN 6
#define OCT_DIR_LEFT_UP 7
#define OCT_DIR_LEFT_BACK 8
#define OCT_DIR_LEFT_FRONT 9


#define OCT_DIR_RIGHT_DOWN 10
#define OCT_DIR_RIGHT_UP 11
#define OCT_DIR_RIGHT_BACK 12
#define OCT_DIR_RIGHT_FRONT 13

#define OCT_DIR_DOWN_BACK 14
#define OCT_DIR_DOWN_FRONT 15


#define OCT_DIR_UP_BACK 16
#define OCT_DIR_UP_FRONT 17

#define OCT_DIR_LEFT_DOWN_BACK 18
#define OCT_DIR_RIGHT_DOWN_BACK 19
#define OCT_DIR_LEFT_UP_BACK 20
#define OCT_DIR_RIGHT_UP_BACK 21
#define OCT_DIR_LEFT_DOWN_FRONT 22
#define OCT_DIR_RIGHT_DOWN_FRONT 23
#define OCT_DIR_LEFT_UP_FRONT 24
#define OCT_DIR_RIGHT_UP_FRONT 25
#define OCT_DIR_INTERNAL 26

#define OCT_DIR_TOTAL 27

/** variable to ensure the element sz % element order =0 for higher order elemnts */
extern unsigned int MAXDEAPTH_LEVEL_DIFF;


#define NUM_LEVEL_BITS 5u

#ifdef DIM_2
#define m_uiDim 2
#else
#define m_uiDim 3
#endif

#ifdef DIM_2
#define NUM_CHILDREN 4
#define NUM_FACES 1
#define NUM_EDGES 4
#define EDGE_OFFSET 4
#define VERTEX_OFFSET 4
#define ROTATION_OFFSET 8
#define MAX_HANGING_FACES 2
#define F2E_MAP_OFFSET 4
#else
#define NUM_CHILDREN 8
#define NUM_FACES 6
#define NUM_EDGES 12
#define EDGE_OFFSET 6
#define VERTEX_OFFSET 18
#define ROTATION_OFFSET 16
#define MAX_HANGING_FACES 4
#define F2E_MAP_OFFSET 6
#endif

#define F2E_FACE_DIR_OFFSET 3
#define F2E_NUM_MAX_FACE_DIR 7
#define F2E_FACE_INDEPEN_BIT 3
#define F2E_FACE_DEPEN_BIT 4


// AMR coarsening factor.
#define DENDRO_AMR_COARSEN_FAC 0.1

#define DENDRO_DEFAULT_LB_TOL 0.1
#define DENDRO_DEFAULT_SF_K 2
#define DENDRO_DEFAULT_GRAIN_SZ 100


#define DENDRO_UNSIGNED_INT_MIN UINT_MIN
#define DENDRO_UNSIGNED_INT_MAX UINT_MAX

#define DENDRO_REMESH_UNZIP_SCALE_FAC 1.0


#define DENDRO_BLOCK_ALIGN_FACTOR 1
#define DENDRO_BLOCK_ALIGN_FACTOR_LOG 0


#define ODA_INDEPENDENT_FLAG_BIT 0
#define ODA_W_DEPENDENT_FLAG_BIT 1
#define ODA_W_BOUNDARY_FLAG_BIT 2
// used in face to element map.
#define ODA_FLAGS_TOTAL 3

void __handler(int sig);

inline int dendro_error(const char* const file, int line, const std::string& msg){
   std::cout<< "[" << file << "] : "<<line<<" "<< msg <<std::endl;
   return 0;
}

#define dendro_log(msg) dendro_error(__FILE__, __LINE__, msg)