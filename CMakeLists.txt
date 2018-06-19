cmake_minimum_required(VERSION 2.8)
project(DendroGR)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

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


option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" ON)
option(DIM_2 "use the two dimentional sorting" OFF)
option(WITH_BLAS_LAPACK "build using BLAS and LAPACk" OFF)
option(MANUAL_BLAS_LAPACK "configure BLAS and LAPACK Manually" OFF)
option(DENDRO_VTK_BINARY "write vtk/vtu files in binary mode " ON)
option(DENDRO_VTK_ZLIB_COMPRES "write vtk/vtu files in binary mode with zlib compression (only compatible with binary mode) " OFF)
option(ALLTOALL_SPARSE "uses isend irecv for ghost exchange (for some mpi implementations this might be efficient)" OFF)
option(ENABLE_DENDRO_PROFILE_COUNTERS " enables dendro internal profile counters. " OFF)
option(RK_SOLVER_OVERLAP_COMM_AND_COMP " enables RK solver comm. & comp. overlap. (might be good in large scale runs) " ON)
option(NONLINEAR "If ON solves the Non-linear wave equation" OFF)

#option(KWAY "K parameter for alltoallv_kway" 128)
set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 2 CACHE INT 2)
#message(${KWAY})

if(WITH_BLAS_LAPACK)
    add_definitions(-DWITH_BLAS_LAPACK)
if(MANUAL_BLAS_LAPACK)
    if("$ENV{BLAS}" STREQUAL "")
        message("Environment Variable BLAS is not set. Please set it to BLAS directory")
    endif()

    if( "$ENV{LAPACK}" STREQUAL "" )
        message("Enviroment Variable LAPACK is note set. Please set it to LAPACK directory. ")
    endif()
    set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
    set(BLAS_LIBS $ENV{BLAS}/lib )
    #set(CBLAS_LIBS $ENV{CBLAS}/lib )
    set(LAPACK_LIBS $ENV{LAPACK}/lib)
    include_directories(${LAPACKE_DIR}/include)
    #set(LAPACKE_DIR /home/milinda/Programs/lapack-3.7.0/LAPACKE)
    #set(BLAS_LIBS /home/milinda/Programs/BLAS-3.7.0/lib )
    #set(CBLAS_LIBS $ENV{CBLAS}/lib )
    #set(LAPACK_LIBS /home/milinda/Programs/lapack-3.7.0/lib)

    set(LAPACK_LINKER_FLAGS -llapacke -llapack -lblas -lgfortran -lquadmath)
    set(LAPACK_LIBRARIES ${LAPACK_LIBS}/liblapacke.a ${LAPACK_LIBS}/liblapack.a ${BLAS_LIBS}/libblas.a libgfortran.a libquadmath.a)
    set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
else ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    set(LAPACKE_DIR $ENV{LAPACK}/LAPACKE)
    include_directories(${LAPACK_INCLUDE_DIR})
    #link_directories(${MKL_LIBRARIES})
    set(LINK_FLAGS "${LINK_FLAGS} ${LAPACK_LINKER_FLAGS}")
endif()

endif()

if(DIM_2)
    add_definitions(-DDIM_2)
endif()

if(USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)
    #message('Configuring 64BIT indices')
endif()

if(HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
    #message (Configuring with HIlbert)
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

if(NONLINEAR)
add_definitions(-DNONLINEAR)
endif()


set(LAPACK_INC /home/milinda/Programs/lapack-3.7.0/LAPACKE/include/)
include_directories(${PROJECT_BINARY_DIR}
        include
	FEM/include
        ODE/include
        IO/vtk/include
        IO/zlib/inc
        include/test
        ${MPI_INCLUDE_PATH}
        )

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
        include/test/testUtils.h
        include/test/testUtils.tcc
        include/test/meshTestUtils.h
        include/octUtils.h
        include/dendroIO.h
        include/mathUtils.h
        include/mathUtils.tcc
        include/stencil.h
        include/dollar.hpp
        include/block.h
        include/sc_flops.h
        include/wavelet.h
        include/profiler.h
        include/dendroProfileParams.h
        include/daUtils.h
        include/daUtils.tcc
        )

set(DENDRO_SOURCE_FILES
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
        )

set(ODE_INC ODE/include/fdCoefficient.h
            ODE/include/rk.h
            ODE/include/rkMaxwell.h
            ODE/include/rkTransport.h
            ODE/include/rkTransportUtils.h
        )

set(ODE_SRC ODE/src/rk.cpp
            ODE/src/rkMaxwell.cpp
            ODE/src/rkTransport.cpp
        )


set(IO_INC IO/vtk/include/oct2vtk.h
           IO/vtk/include/json.hpp
           IO/vtk/include/checkPoint.h
           IO/vtk/include/cencode.h )

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


set(FEM_INC 
        FEM/include/basis.h
        FEM/include/tensor.h
	FEM/include/refel.h
        )

set(FEM_SRC 
        FEM/src/basis.cpp
        FEM/src/refel.cpp
        FEM/src/tensor.cpp
        )



if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()
 
add_library(dendro5 ${DENDRO_INCLUDE_FILES} ${DENDRO_SOURCE_FILES} ${FEM_INC} ${FEM_SRC} ${ODE_INC} ${ODE_SRC} ${IO_INC} ${IO_SRC} ${ZLIB_INC} ${ZLIB_SRC})

set(SOURCE_FILES ODE/examples/transportEq.cpp)
add_executable(transportEq ${SOURCE_FILES})
target_link_libraries(transportEq dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)


set(NLSM_INC NLSigma/include/parameters.h
             NLSigma/include/rk4nlsm.h
             NLSigma/include/nlsmUtils.h
             NLSigma/include/rhs.h
             NLSigma/include/derivs.h
             NLSigma/include/profile_params.h
       )

set(NLSM_SRC NLSigma/src/rk4nlsm.cpp
             NLSigma/src/parameters.cpp
             NLSigma/src/nlsmUtils.cpp
             NLSigma/src/rhs.cpp
             NLSigma/src/derivs.cpp
             NLSigma/src/profile_params.cpp
             )

set(SOURCE_FILES NLSigma/src/nlsm.cpp  NLSigma/include/nlsm.h ${NLSM_INC} ${NLSM_SRC})
add_executable(nlsmSolver ${SOURCE_FILES})
target_include_directories(nlsmSolver PUBLIC NLSigma/include)
target_link_libraries(nlsmSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

