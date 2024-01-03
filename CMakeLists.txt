# ##############################################################################
# BASE CMAKE CONFIGURATION
#
# This file in particular contains all of the options that apply to all aspects
# of the program, including some Dendro options that are used when including the
# library.
#
# Note that all options and compiler flags set here will propagate to other
# directories such as src/ and dendro/
# ##############################################################################

cmake_minimum_required(VERSION 3.6)
project(Dendro-5.0 C CXX Fortran)

# the standard we're using is c++11 -- NOTE: this might change to C++17 or newer
# in the future
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# OpenMP, MPI, and GSL are required packages for all aspects of the project
find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(GSL REQUIRED)

# For now we just make it compulsory to have LAPACK installed. TODO: make Lapack
# install if not present.

# flags for OPENMP, these are required
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS
      "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(MPI_COMPILE_FLAGS)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()

# options for dendro -- not sure if these need to be declared here, but they might be due to header definitions...
option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off"
       ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)
option(RUN_WEAK_SCALING "Run Entire Weak Scaling" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" OFF)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" OFF)
option(DIM_2 "use the two dimentional sorting" OFF)
option(PROFILE_TREE_SORT
       "Profile the tree sort code, construction and balancing. " OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" ON)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
option(DENDRO_VTK_BINARY "write vtk/vtu files in binary mode " ON)
option(
  DENDRO_VTK_ZLIB_COMPRES
  "write vtk/vtu files in binary mode with zlib compression (only compatible with binary mode) "
  OFF)
option(
  ALLTOALL_SPARSE
  "uses isend irecv for ghost exchange (for some mpi implementations this might be efficient)"
  OFF)
option(ENABLE_DENDRO_PROFILE_COUNTERS
       " enables dendro internal profile counters. " OFF)
option(
  RK_SOLVER_OVERLAP_COMM_AND_COMP
  " enables RK solver comm. & comp. overlap. (might be good in large scale runs) "
  ON)
option(WITH_CUDA " build dendro with cuda " OFF)
option(BUILD_WITH_PETSC " build dendro with PETSC " OFF)
option(USE_FD_INTERP_FOR_UNZIP "use FD style interpolation for unzip" OFF)
option(MPIX_CUDA_AWARE_SUPPORT "use FD style interpolation for unzip" OFF)

# option(KWAY "K parameter for alltoallv_kway" 128)
set(KWAY
    128
    CACHE STRING 128)
set(OCT2BLK_COARSEST_LEV
    0
    CACHE STRING 0)
set(NUM_NPES_THRESHOLD
    2
    CACHE STRING 2)

# set the build type to release by default.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose build type: 'Debug' or 'Release' " FORCE)
endif()

if(WITH_CUDA)
  set(OCT2BLK_COARSEST_LEV
      31
      CACHE STRING 31)
  add_definitions(-DOCT2BLK_COARSEST_LEV=${OCT2BLK_COARSEST_LEV})
  find_package(CUDA REQUIRED)
  if(CUDA_FOUND)
    set(CUDA_PROPAGATE_HOST_FLAGS ON)
    set(CUDA_SEPARABLE_COMPILATION ON)
    list(
      APPEND
      CUDA_NVCC_FLAGS
      -std=c++11
      --ptxas-options=-O3
      --ptxas-options=-v
      -w
      -lineinfo)
    add_definitions(-DWITH_CUDA)
  endif()
  if(MPIX_CUDA_AWARE_SUPPORT)
    add_definitions(-DMPIX_CUDA_AWARE_SUPPORT)
  endif()
endif()
add_definitions(-DOCT2BLK_COARSEST_LEV=${OCT2BLK_COARSEST_LEV})

if(WITH_BLAS_LAPACK)
  add_definitions(-DWITH_BLAS_LAPACK)
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
  set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
  set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
  set(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} ${LAPACKE_LIB})
  message(STATUS ${LAPACK_LIBRARIES})

  if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(
      "${CMAKE_CXX_COMPILER_ID} compiler detected adding -mkl flag for BLAS LAPACK"
    )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
  endif()

endif()

# build options for PETSC
if(BUILD_WITH_PETSC)
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")
  find_package(PETSc REQUIRED)
  add_definitions(-DBUILD_WITH_PETSC)
endif()

if(DIM_2)
  add_definitions(-DDIM_2)
endif()

if(PROFILE_TREE_SORT)
  add_definitions(-DPROFILE_TREE_SORT)
endif()

if(USE_64BIT_INDICES)
  add_definitions(-DUSE_64BIT_INDICES)
  # message('Configuring 64BIT indices')
endif()

if(HILBERT_ORDERING)
  add_definitions(-DHILBERT_ORDERING)
  # message (Configuring with HIlbert)
endif()

if(RUN_WEAK_SCALING)
  add_definitions(-DRUN_WEAK_SCALING)
endif()

if(ALLTOALLV_FIX)
  add_definitions(-DALLTOALLV_FIX)
  add_definitions(-DKWAY=${KWAY})
endif()

if(SPLITTER_SELECTION_FIX)
  add_definitions(-DSPLITTER_SELECTION_FIX)
  add_definitions(-DNUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif()

if(ALLTOALL_SPARSE)
  add_definitions(-DALLTOALL_SPARSE)
endif()

if(ENABLE_DENDRO_PROFILE_COUNTERS)
  add_definitions(-DENABLE_DENDRO_PROFILE_COUNTERS)
  add_definitions(-D__PROFILE_CTX__)
  add_definitions(-D__PROFILE_ETS__)
  # add_definitions(-D__PROFILE_ENUTS__)
endif()

if(RK_SOLVER_OVERLAP_COMM_AND_COMP)
  add_definitions(-DRK_SOLVER_OVERLAP_COMM_AND_COMP)
endif()

if(DENDRO_VTK_BINARY)

else()
  set(DENDRO_VTK_ZLIB_COMPRES OFF)
endif()

if(DENDRO_VTK_BINARY)
  add_definitions(-DDENDRO_VTU_BINARY)
if(DENDRO_VTK_ZLIB_COMPRES)
    add_definitions(-DDENDRO_VTU_ZLIB)
  endif()
else()
  add_definitions(-DDENDRO_VTU_ASCII)
endif()

if(USE_FD_INTERP_FOR_UNZIP)
  add_definitions(-DUSE_FD_INTERP_FOR_UNZIP)
endif()

add_definitions(-DMATVEC_PROFILE)

set(CHECK_GIT_INCLUDE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)
include(${CHECK_GIT_INCLUDE_PATH}/CheckGit.cmake)
message(STATUS "Setting up adding Git Information directly to compile options, header file should be in ${CHECK_GIT_INCLUDE_PATH}")
CheckGitSetup()

# === END REQUIRED CONFIGURATIONS FOR DENDRO ===

# === BEGIN SUBDIRECTORY INCLUSION FOR BUILDING ===

# dendro folder *MUST* be included to build the library first
add_subdirectory(dendrolib)

# additional subdirectories that are their own "solvers"
add_subdirectory(BSSN_GR)
# add_subdirectory(CCZ4_GR)
add_subdirectory(NLSigma)
add_subdirectory(Fluid)
add_subdirectory(EM1)
add_subdirectory(EM2)
add_subdirectory(EM3)
add_subdirectory(EM4)
