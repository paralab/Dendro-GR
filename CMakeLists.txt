cmake_minimum_required(VERSION 3.1)
project(Dendro-5.0 C CXX Fortran)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(GSL REQUIRED)


# For now we just make it compulsory to have LAPACK installed.
#Later we will make it possible if LAPACK is not present to automaticall install before compiling dendro5
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()


## options for dendro

option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)
option(RUN_WEAK_SCALING "Run Entire Weak Scaling" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" OFF)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" OFF)
option(DIM_2 "use the two dimentional sorting" OFF)
option(PROFILE_TREE_SORT "Profile the tree sort code, construction and balancing. " OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" ON)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
option(DENDRO_VTK_BINARY "write vtk/vtu files in binary mode " ON)
option(DENDRO_VTK_ZLIB_COMPRES "write vtk/vtu files in binary mode with zlib compression (only compatible with binary mode) " OFF)
option(ALLTOALL_SPARSE "uses isend irecv for ghost exchange (for some mpi implementations this might be efficient)" OFF)
option(ENABLE_DENDRO_PROFILE_COUNTERS " enables dendro internal profile counters. " OFF)
option(RK_SOLVER_OVERLAP_COMM_AND_COMP " enables RK solver comm. & comp. overlap. (might be good in large scale runs) " ON)
option(WITH_CUDA " build dendro with cuda " OFF)
option(BUILD_WITH_PETSC " build dendro with PETSC " OFF)
option(USE_FD_INTERP_FOR_UNZIP "use FD style interpolation for unzip" OFF)
option(MPIX_CUDA_AWARE_SUPPORT "use FD style interpolation for unzip" OFF)


#option(KWAY "K parameter for alltoallv_kway" 128)
set(KWAY 128 CACHE STRING 128)
set(OCT2BLK_COARSEST_LEV 0 CACHE STRING 0)
set(NUM_NPES_THRESHOLD 2 CACHE STRING 2)

#set the build type to release by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release " FORCE)
endif()


if(WITH_CUDA)
    set(OCT2BLK_COARSEST_LEV 31 CACHE STRING 31)
    add_definitions(-DOCT2BLK_COARSEST_LEV=${OCT2BLK_COARSEST_LEV})
    find_package(CUDA REQUIRED)
    if(CUDA_FOUND)
        set(CUDA_PROPAGATE_HOST_FLAGS ON)
        set(CUDA_SEPARABLE_COMPILATION ON)
        list(APPEND CUDA_NVCC_FLAGS -std=c++11 --ptxas-options=-O3 --ptxas-options=-v -w -lineinfo )
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

    if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        message("${CMAKE_CXX_COMPILER_ID} compiler detected adding -mkl flag for BLAS LAPACK")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl")
    endif()

endif()

if(BUILD_WITH_PETSC)
    list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules")
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
    #message('Configuring 64BIT indices')
endif()

if(HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
    #message (Configuring with HIlbert)
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
    #add_definitions(-D__PROFILE_ENUTS__)
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

set(DENDRO_INCLUDE_FILES
        include/point.h
        include/pointT.h
        include/seqUtils.h
        include/seqUtils.tcc
        include/ompUtils.h
        include/ompUtils.tcc
        include/dendro.h
        include/binUtils.h
        include/hcurvedata.h
        include/genPts_par.h
        include/TreeNode.h
        include/key.h
        include/skey.h
        include/node.h
        include/treenode2vtk.h
        include/parUtils.h
        include/parUtils.tcc
        include/mesh.h
        include/mesh.tcc
        include/meshE2NUtils.tcc
        include/dendro.h
        include/colors.h
        include/dtypes.h
        include/sfcSort.h
        include/sfcSearch.h
        include/octUtils.h
        include/dendroIO.h
        include/mathMeshUtils.h
        include/mathUtils.h
        include/mathUtils.tcc
        include/memory_pool.h
        include/stencil.h
        include/dollar.hpp
        include/block.h
        include/sc_flops.h
        include/wavelet.h
        include/profiler.h
        include/dendroProfileParams.h
        include/daUtils.h
        include/daUtils.tcc
        include/oda.h
        include/oda.tcc
        include/odaUtils.h
        include/asyncExchangeContex.h
        include/matRecord.h
        include/nufft.h
        include/sub_oda.h
        include/sub_oda.tcc
        include/dvec.h
        include/meshUtils.h
        include/waveletRefEl.h
        include/waveletAMR.h
        include/waveletAMR.tcc
        test/include/testUtils.h
        test/include/testUtils.tcc
        test/include/meshTestUtils.h
        include/launcher.h
        include/device.h
        )

set(DENDRO_SOURCE_FILES
        src/dendro.cpp
        src/point.cpp
        src/binUtils.cpp
        src/hcurvedata.cpp
        src/genPts_par.cpp
        src/TreeNode.cpp
        src/treenode2vtk.cpp
        src/parUtils.cpp
        src/mesh.cpp
        src/octUtils.cpp
        src/dendroIO.cpp
        src/block.cpp
        src/sc_flops.cpp
        src/profiler.cpp
        src/dendroProfileParams.cpp
        src/oda.cpp
        src/odaUtils.cpp
        src/sub_oda.cpp
        src/meshUtils.cpp
        src/waveletRefEl.cpp
        test/src/meshTestUtils.cpp
        )


set(FEM_INC FEM/include/basis.h
        FEM/include/lapac.h
        FEM/include/tensor.h
        FEM/include/matvec.h
        FEM/include/workspace.h
        FEM/include/operators.h
        FEM/include/feMat.h
        FEM/include/feMatrix.h
        FEM/include/feVec.h
        FEM/include/feVector.h
        )

set(FEM_SRC FEM/examples/src/laplace.cpp
        FEM/src/basis.cpp
        FEM/src/refel.cpp
        FEM/src/tensor.cpp
        FEM/src/workspace.cpp
        )

set(LINALG_INC
        LinAlg/include/cg.h
        )

set(ODE_INC ODE/include/fdCoefficient.h
            ODE/include/rk.h
            ODE/include/rkMaxwell.h
            ODE/include/rkTransport.h
            ODE/include/rkTransportUtils.h
            ODE/include/ts.h
            ODE/include/blkAsync.h
            ODE/include/ctx.h
            ODE/include/ets.h
            ODE/include/enuts.h
            ODE/include/enutsUtils.tcc
            ODE/include/asyncCtxByLev.h
            ODE/include/subSM.h
            ODE/include/enutsOp.h
            ODE/include/sdc.h
            ODE/include/inverseCtx.h
        )

set(ODE_SRC ODE/src/rk.cpp
            ODE/src/subSM.cpp
            ODE/src/rkMaxwell.cpp
            ODE/src/rkTransport.cpp
            ODE/src/enutsOp.cpp
        )


set(IO_INC IO/vtk/include/oct2vtk.h
           IO/vtk/include/json.hpp
           IO/vtk/include/checkPoint.h
           IO/vtk/include/cencode.h
           IO/include/rawIO.h)

set(IO_SRC IO/vtk/src/oct2vtk.cpp
           IO/vtk/src/checkPoint.cpp )


set(ZLIB_INC IO/zlib/inc/crc32.h
        IO/zlib/inc/deflate.h
        IO/zlib/inc/gzguts.h
        IO/zlib/inc/inffast.h
        IO/zlib/inc/inffixed.h
        IO/zlib/inc/inflate.h
        IO/zlib/inc/inftrees.h
        IO/zlib/inc/trees.h
        IO/zlib/inc/zlib.h
        IO/zlib/inc/zconf.h
        IO/zlib/inc/zutil.h
        )


set(ZLIB_SRC
        IO/zlib/src/adler32.c
        IO/zlib/src/compress.c
        IO/zlib/src/crc32.c
        IO/zlib/src/deflate.c
        IO/zlib/src/gzclose.c
        IO/zlib/src/gzlib.c
        IO/zlib/src/gzread.c
        IO/zlib/src/gzwrite.c
        IO/zlib/src/infback.c
        IO/zlib/src/inffast.c
        IO/zlib/src/inflate.c
        IO/zlib/src/inftrees.c
        IO/zlib/src/trees.c
        IO/zlib/src/uncompr.c
        IO/zlib/src/zutil.c
        )


if(WITH_CUDA)
    set(DENDRO_INCLUDE_FILES ${DENDRO_INCLUDE_FILES} include/cudaUtils.h)
    set(DENDRO_SOURCE_FILES ${DENDRO_SOURCE_FILES} src/cudaUtils.cpp src/oda.cpp FEM/examples/src/heatMat.cpp)
endif()

if(WITH_CUDA)
    set(SOURCE_FILES    ${DENDRO_SOURCE_FILES} ${FEM_SRC} ${ODE_SRC} ${IO_SRC} ${ZLIB_SRC}
                        GPU/include/derivs_cu.cuh
                        GPU/src/refel_const.cu
                        GPU/src/mesh_gpu.cu)
    cuda_add_library(dendro5_cuda ${SOURCE_FILES})
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/GPU/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
    target_include_directories(dendro5_cuda PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/include)
    target_include_directories(dendro5_cuda PUBLIC ${MPI_INCLUDE_PATH})
    target_include_directories(dendro5_cuda PUBLIC ${CUDA_INCLUDE_DIRS})
    target_link_libraries(dendro5_cuda ${CUDA_LIBRARIES})
endif()


add_library(dendro5  ${DENDRO_SOURCE_FILES} ${FEM_SRC} ${ODE_SRC} ${IO_SRC} ${ZLIB_SRC} )
target_link_libraries(dendro5 PUBLIC ${MPI_LIBRARIES} m)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(dendro5 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/include)
target_include_directories(dendro5 PUBLIC ${MPI_INCLUDE_PATH})
target_include_directories(dendro5 PUBLIC ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(dendro5 PUBLIC ${CUDA_INCLUDE_DIRS})
endif()

if(WITH_BLAS_LAPACK)
    target_include_directories(dendro5 PUBLIC ${LAPACK_INCLUDE_DIR})
    target_link_libraries(dendro5 PUBLIC ${LAPACK_LIBRARIES})
endif()

if(BUILD_WITH_PETSC)
    target_include_directories(dendro5 PUBLIC ${PETSC_INCLUDES})
    target_link_libraries(dendro5 PUBLIC ${PETSC_LIBRARIES})
endif()


if(WITH_CUDA)
    cuda_add_executable(run_meshgpu_tests ${CMAKE_CURRENT_SOURCE_DIR}/GPU/src/run_meshgpu_tests.cu )
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/GPU/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
    target_include_directories(run_meshgpu_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/IO/include)
    target_include_directories(run_meshgpu_tests PUBLIC ${MPI_INCLUDE_PATH})
    target_include_directories(run_meshgpu_tests PUBLIC ${CUDA_INCLUDE_DIRS})
    target_link_libraries(run_meshgpu_tests dendro5 dendro5_cuda)
endif(WITH_CUDA)


set(SOURCE_FILES Hilbert/src/genHTables.cpp include/radix.h)
add_executable(genHtables ${SOURCE_FILES})
target_include_directories(genHtables PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)


set(SOURCE_FILES Hilbert/src/sfcMatvecTable.cpp include/dendro.h)
add_executable(genSFCMatVecTables ${SOURCE_FILES})
target_link_libraries(genSFCMatVecTables m)
target_include_directories(genSFCMatVecTables PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(SOURCE_FILES examples/src/ptsSort.cpp)
add_executable(ptsSort ${SOURCE_FILES})
target_include_directories(ptsSort PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(ptsSort PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(ptsSort PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(ptsSort PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(ptsSort dendro5 ${MPI_LIBRARIES} m)


set(SOURCE_FILES examples/include/testOctUtils.h
                 examples/src/testOCtUtils.cpp)
add_executable(testOctUtils ${SOURCE_FILES})
target_link_libraries(testOctUtils dendro5 ${MPI_LIBRARIES} m)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(testOctUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(testOctUtils PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(testOctUtils PRIVATE ${GSL_INCLUDE_DIRS})

set(SOURCE_FILES examples/src/tstTreeSearch.cpp)
add_executable(tstTreeSearch ${SOURCE_FILES})
target_link_libraries(tstTreeSearch dendro5 ${MPI_LIBRARIES} m)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(tstTreeSearch PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(tstTreeSearch PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(tstTreeSearch PRIVATE ${GSL_INCLUDE_DIRS})



set(SOURCE_FILES examples/src/octUtilsBenchmark.cpp)
add_executable(octUtilBenchmark ${SOURCE_FILES})
target_link_libraries(octUtilBenchmark dendro5 ${MPI_LIBRARIES} m)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(octUtilBenchmark PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(octUtilBenchmark PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(octUtilBenchmark PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES examples/src/bucketBench.cpp)
add_executable(bucketBench ${SOURCE_FILES})
target_link_libraries(bucketBench dendro5 ${MPI_LIBRARIES} m)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(bucketBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(bucketBench PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(bucketBench PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES examples/src/mfreeTest.cpp)
add_executable(mfreeBench ${SOURCE_FILES})
target_link_libraries(mfreeBench dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(mfreeBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(mfreeBench PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(mfreeBench PRIVATE ${GSL_INCLUDE_DIRS})



set(SOURCE_FILES examples/src/meshBenchmark.cpp)
add_executable(meshBench ${SOURCE_FILES})
target_link_libraries(meshBench dendro5  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(meshBench PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(meshBench PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(meshBench PRIVATE ${GSL_INCLUDE_DIRS})

set(SOURCE_FILES examples/src/meshE2NCheck.cpp)
add_executable(meshE2NCheck ${SOURCE_FILES})
target_link_libraries(meshE2NCheck dendro5  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(meshE2NCheck PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(meshE2NCheck PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(meshE2NCheck PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES test/src/run_all_tests.cpp)
add_executable(run_all_tests ${SOURCE_FILES})
target_link_libraries(run_all_tests dendro5  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(run_all_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(run_all_tests PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(run_all_tests PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES test/src/benchUtils.cpp)
add_executable(benchUtils ${SOURCE_FILES})
target_link_libraries(benchUtils dendro5  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(benchUtils PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(benchUtils PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(benchUtils PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES ODE/examples/testLTSBlkVec.cpp)
add_executable(testLTSBlkVec ${SOURCE_FILES})
target_link_libraries(testLTSBlkVec dendro5  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(testLTSBlkVec PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(testLTSBlkVec PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(testLTSBlkVec PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES FEM/examples/src/dgVecTest.cpp)
add_executable(dgVecTest ${SOURCE_FILES})
target_link_libraries(dgVecTest dendro5  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(dgVecTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(dgVecTest PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(dgVecTest PRIVATE ${GSL_INCLUDE_DIRS})



set(SOURCE_FILES FEM/examples/src/testInterp.cpp)
add_executable(testInterp ${SOURCE_FILES})
target_link_libraries(testInterp dendro5  ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(testInterp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(testInterp PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(testInterp PRIVATE ${GSL_INCLUDE_DIRS})




set(SOURCE_FILES examples/src/meshLooping.cpp )
add_executable(meshLoop ${SOURCE_FILES})
target_link_libraries(meshLoop dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(meshLoop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(meshLoop PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(meshLoop PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES FEM/examples/src/refelExample.cpp)
add_executable(refEl ${SOURCE_FILES})
target_link_libraries(refEl dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(refEl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(refEl PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(refEl PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES FEM/examples/src/refelOctreeExample.cpp)
add_executable(refElOctree ${SOURCE_FILES})
target_link_libraries(refElOctree dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(refElOctree PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(refElOctree PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(refElOctree PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES ODE/examples/transportEq.cpp)
add_executable(transportEq ${SOURCE_FILES})
target_link_libraries(transportEq dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(transportEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(transportEq PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(transportEq PRIVATE ${GSL_INCLUDE_DIRS})

set(SOURCE_FILES FEM/examples/src/laplace.cpp)
add_executable(laplaceEq ${SOURCE_FILES})
target_link_libraries(laplaceEq dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(laplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(laplaceEq PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(laplaceEq PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES FEM/examples/src/heatMat.cpp
                 FEM/examples/src/heatVec.cpp
                 FEM/examples/src/heatEq.cpp)
add_executable(heatEq ${SOURCE_FILES})
target_link_libraries(heatEq dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(heatEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(heatEq PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(heatEq PRIVATE ${GSL_INCLUDE_DIRS})


set(SOURCE_FILES FEM/examples/src/subdaExample.cpp)
add_executable(subDAExample ${SOURCE_FILES})
target_link_libraries(subDAExample dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/examples/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(subDAExample PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(subDAExample PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(subDAExample PRIVATE ${GSL_INCLUDE_DIRS})



set(SOURCE_FILES FEM/examples/src/sfcLaplace.cpp)
add_executable(sfcLaplaceEq ${SOURCE_FILES})
target_link_libraries(sfcLaplaceEq dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test/include)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples/include)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/FEM/include)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ODE/include)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/LinAlg/include)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/vtk/include)
target_include_directories(sfcLaplaceEq PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(sfcLaplaceEq PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(sfcLaplaceEq PRIVATE ${GSL_INCLUDE_DIRS})


add_subdirectory(BSSN_GR)
#add_subdirectory(CCZ4_GR)
add_subdirectory(NLSigma)
add_subdirectory(Fluid)
add_subdirectory(EM1)
add_subdirectory(EM2)
add_subdirectory(EM3)
add_subdirectory(EM4)
