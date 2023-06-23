# ##############################################################################
# CMAKE CONFIGURATION FILE FOR BSSN SOLVER
#
# This file in particular contains all of the options that apply specificially
# to the BSSN Solver binary as it is compiled.
#
# Note that all options and compiler flags set in the ../CMakeLists.txt file
# (the one in the root directory of the project) *are carried through* here.
# Meaning, at this point, the C++ standard is set, OpenMP, MPI, GSL, and other
# libraries are found and linked. Also, general Dendro options, such as ordering
# options, dimensionality, CMake's build type, and others will also propagate
# through to here.
# ##############################################################################

# only BSSN (GR) options
option(BSSN_PROFILE_HUMAN_READABLE " enables human readable profile outputs" ON)
option(BSSN_CONSEC_COMM_SELECT "rank selection performed using consecative ranks (works with arbitary number of procs.) if off ranks selected in binary tree fashion. (works only with powers of 2)" ON)
option(BSSN_ENABLE_VTU_OUTPUT "Enables vtu output" ON)
option(BSSN_COMPUTE_CONSTRAINTS "Enable constraint equation computations" ON)
option(BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT "write constrints to vtu" ON)
option(BSSN_ETA_FUNCTION "Use function as ETA damping" OFF)
option(BSSN_KERR_SCHILD_TEST "compare with the analytical Kerr-Schild sol. for a single BH" OFF)
option(BSSN_ENABLE_CUDA "enable rhs computation with GPUs" OFF)
#option(BSSN_RHS_STAGED_COMP "enable rhs staged computation" OFF)
option(BSSN_GAUGE_ROCHESTER "enable rochester gauge for bssn equations" OFF)
option(BSSN_EXTRACT_BH_LOCATIONS "compute the BH locations as time evolves" ON)
option(BSSN_REFINE_BASE_EH "enable refine based on BH event horizons" OFF)
option(BSSN_EXTRACT_GRAVITATIONAL_WAVES "extract GW waves" ON)
option(BSSN_USE_4TH_ORDER_DERIVS "use 4th order derivs" OFF)
option(BSSN_USE_6TH_ORDER_DERIVS "use 6th order derivs" ON)
option(BSSN_USE_8TH_ORDER_DERIVS "use 8th order derivs" OFF)
option(BSSN_ENABLE_AVX "use vectorized computation for RHS (intel compilers only)" OFF)
option(BSSN_VERBOSE_LOGGING "Enable more verbose logging (only useful during debug)" OFF)


if(BSSN_ENABLE_AVX)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        add_definitions(-DBSSN_ENABLE_AVX)
        add_definitions(-DDERIV_ENABLE_AVX)
    else()
        message("AVX vectorization currently supported to intel compilers only. use export CC=icc and export CXX=icpc to use intel compilers. AVX vectorization not enabled. ")
    endif()
    
endif(BSSN_ENABLE_AVX)

if(BSSN_USE_4TH_ORDER_DERIVS)
    set(BSSN_USE_6TH_ORDER_DERIVS OFF)
    set(BSSN_USE_8TH_ORDER_DERIVS OFF)
elseif(BSSN_USE_6TH_ORDER_DERIVS)
    set(BSSN_USE_4TH_ORDER_DERIVS OFF)
    set(BSSN_USE_8TH_ORDER_DERIVS OFF)
elseif(BSSN_USE_8TH_ORDER_DERIVS)
    set(BSSN_USE_4TH_ORDER_DERIVS OFF)
    set(BSSN_USE_6TH_ORDER_DERIVS OFF)
endif()

if(BSSN_USE_4TH_ORDER_DERIVS)
    add_definitions(-DBSSN_USE_4TH_ORDER_DERIVS)
endif()

if(BSSN_USE_6TH_ORDER_DERIVS)
    add_definitions(-DBSSN_USE_6TH_ORDER_DERIVS)
endif()

if(BSSN_USE_8TH_ORDER_DERIVS)
    add_definitions(-DBSSN_USE_8TH_ORDER_DERIVS)
endif()

if(BSSN_PROFILE_HUMAN_READABLE)
    add_definitions(-DBSSN_PROFILE_HUMAN_READABLE)
endif()

if(BSSN_CONSEC_COMM_SELECT)
    add_definitions(-DBSSN_CONSEC_COMM_SELECT)
endif()

if(BSSN_ENABLE_VTU_OUTPUT)
    add_definitions(-DBSSN_ENABLE_VTU_OUTPUT)
endif()

if(BSSN_COMPUTE_CONSTRAINTS)
    add_definitions(-DBSSN_COMPUTE_CONSTRAINTS)
endif()

if(BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT)
    add_definitions(-DBSSN_ENABLE_VTU_CONSTRAINT_OUTPUT)
endif()

if(BSSN_ETA_FUNCTION)
    add_definitions(-DUSE_ETA_FUNC)
endif()

if(BSSN_KERR_SCHILD_TEST)
    add_definitions(-DBSSN_KERR_SCHILD_TEST)
endif()

if(BSSN_ENABLE_CUDA)
    add_definitions(-DBSSN_ENABLE_CUDA)
endif()

if(BSSN_RHS_STAGED_COMP)
    add_definitions(-DBSSN_RHS_STAGED_COMP)
endif()

if(BSSN_EXTRACT_BH_LOCATIONS)
    add_definitions(-DBSSN_EXTRACT_BH_LOCATIONS)
endif()


if(BSSN_REFINE_BASE_EH)
    add_definitions(-DBSSN_REFINE_BASE_EH)
endif()

if(BSSN_GAUGE_ROCHESTER)
    add_definitions(-DUSE_ROCHESTER_GAUGE)
endif()

if(BSSN_VERBOSE_LOGGING)
    add_definitions(-DBSSN_VERBOSE_LOGGING)
endif()


#message(${KWAY})

# some debug flags. uncomment to print the debug information.
#add_definitions(-DDEBUG_RK_SOLVER)
#add_definitions(-DDEBUG_UPWIND_INTERP)
#add_definitions(-DDEBUG_DOWNWIND_INTERP)
#add_definitions(-DDEBUG_DERIVS_COMP)
#add_definitions(-DMATVEC_PROFILE)
#add_definitions(-DFEM_ACCUMILATE_ONES_TEST)
#add_definitions(-DNLSM_DIRICHLET_BDY)

# generates the code to extract GW waves,

# 12/26/20 : removed automatic cmake based code gen for GWExtraction. 
# add_custom_command(
#         OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/gwExtract.h
#         DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/gwExtract.h.in
#         COMMAND python3 -W ignore ${CMAKE_SOURCE_DIR}/CodeGen/cog.py ${CMAKE_CURRENT_SOURCE_DIR}/include/gwExtract.h.in > ${CMAKE_CURRENT_BINARY_DIR}/gwExtract.h
#         COMMENT "Generating gwExtract.h file "
#         VERBATIM USES_TERMINAL
# )

# add_custom_target(GW_EXTRACT
#         DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/gwExtract.h
# )


add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_const_standard_gauge.cpp
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py 
    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py unstaged standard const . > ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_const_standard_gauge.cpp
    COMMENT "generating bssn code unstaged standard gauge with constant eta damping (during intial compilation this might take a while. )"
    VERBATIM USES_TERMINAL
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_func_standard_gauge.cpp
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py 
    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py unstaged standard func . > ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_func_standard_gauge.cpp
    COMMENT "generating bssn code unstaged standard gauge with varying (large mass ratios) eta damping (during intial compilation this might take a while. )"
    VERBATIM USES_TERMINAL
)


add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_const_rochester_gauge.cpp
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py 
    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py unstaged rochester const . > ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_const_rochester_gauge.cpp
    COMMENT "generating bssn code unstaged rochester gauge with constant eta damping (during intial compilation this might take a while. )"
    VERBATIM USES_TERMINAL
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_func_rochester_gauge.cpp
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/dendro.py 
    DEPENDS ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py 
    COMMAND python3 ${CMAKE_SOURCE_DIR}/CodeGen/bssn.py unstaged rochester func . > ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_func_rochester_gauge.cpp
    COMMENT "generating bssn code unstaged rochester gauge with constant eta damping (during intial compilation this might take a while. )"
    VERBATIM USES_TERMINAL
)


add_custom_target(BSSN_CODEGEN
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_const_standard_gauge.cpp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_func_standard_gauge.cpp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_const_rochester_gauge.cpp
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/bssneqs_eta_func_rochester_gauge.cpp
    )


set(BSSN_INC include/parameters.h
        include/bh.h
        include/grUtils.h
        include/grUtils.tcc
        include/rhs.h
        include/derivs.h
        include/physcon.h
        include/profile_params.h
        include/bssn_constraints.h
        include/TPUtilities.h
        include/TwoPunctures.h
        include/dataUtils.h
        include/aeh_solver.h
        #include/bssnCtx.h
        include/trumpet.h
        )

set(BSSN_SRC src/rkBSSN.cpp
           src/parameters.cpp
           src/grUtils.cpp
           src/rhs.cpp
           src/derivs.cpp
           src/physcon.cpp
           src/profile_params.cpp
           src/bssn_constraints.cpp
           src/TPCoordTransf.cpp
           src/TPEquations.cpp
           src/TPFuncAndJacobian.cpp
           src/TPNewton.cpp
           src/TPUtilities.cpp
           src/TwoPunctures.cpp
           src/dataUtils.cpp
           src/hadrhs.F90
           #src/bssnCtx.cpp
           src/trumpet.cpp
        )

if(WITH_CUDA)
    #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} -gencode arch=compute_70,code=sm_70 --ptxas-options=-v)
    #set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS} --ptxas-options=-v)
    set(CUDA_GR_INCLUDE_H )

    set(CUDA_GR_INC cuda/include/test_utils.h
                    cuda/include/bssn_kernels.cuh
                    cuda/include/bssnrhs_evar_derivs.cuh)

    set(CUDA_GR_SRC )
    
    cuda_add_library(bssn_cuda ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INC} ${CUDA_GR_SRC} )
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
    target_include_directories(bssn_cuda PRIVATE ${MPI_INCLUDE_PATH})
    target_include_directories(bssn_cuda PRIVATE ${GSL_INCLUDE_DIRS})
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
    target_include_directories(bssn_cuda PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    
    cuda_add_executable(derivs_test cuda/src/derivs_test.cu)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
    target_include_directories(derivs_test PRIVATE ${MPI_INCLUDE_PATH})
    target_include_directories(derivs_test PRIVATE ${GSL_INCLUDE_DIRS})
    target_include_directories(derivs_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
    target_include_directories(derivs_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    target_include_directories(derivs_test PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(derivs_test bssn_cuda dendro5_cuda ${CUDA_LIBRARIES})
    target_link_libraries(derivs_test ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

    cuda_add_executable(rhsgpu_test cuda/src/rhs_test.cu)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
    target_include_directories(rhsgpu_test PRIVATE ${MPI_INCLUDE_PATH})
    target_include_directories(rhsgpu_test PRIVATE ${GSL_INCLUDE_DIRS})
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
    target_include_directories(rhsgpu_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    target_include_directories(rhsgpu_test PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(rhsgpu_test bssn_cuda dendro5_cuda ${CUDA_LIBRARIES})
    target_link_libraries(rhsgpu_test ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

    cuda_add_executable(bssnSolverCUDA cuda/src/gr_cuda.cu cuda/src/bssnCtxGPU.cu)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
    target_include_directories(bssnSolverCUDA PRIVATE ${MPI_INCLUDE_PATH})
    target_include_directories(bssnSolverCUDA PRIVATE ${GSL_INCLUDE_DIRS})
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
    target_include_directories(bssnSolverCUDA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    target_include_directories(bssnSolverCUDA PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(bssnSolverCUDA bssn_cuda dendro5_cuda ${CUDA_LIBRARIES})
    target_link_libraries(bssnSolverCUDA ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


    cuda_add_executable(bssnWSTestCUDA cuda/src/gr_scaling.cu cuda/src/bssnCtxGPU.cu)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
    target_include_directories(bssnWSTestCUDA PRIVATE ${MPI_INCLUDE_PATH})
    target_include_directories(bssnWSTestCUDA PRIVATE ${GSL_INCLUDE_DIRS})
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    target_include_directories(bssnWSTestCUDA PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(bssnWSTestCUDA bssn_cuda dendro5_cuda ${CUDA_LIBRARIES})
    target_link_libraries(bssnWSTestCUDA ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

endif()

add_library(bssn_lib ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INC} ${CUDA_GR_SRC} )
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(bssn_lib PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(bssn_lib PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories(bssn_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(bssn_lib PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)

add_executable(rhscpu_test cuda/src/rhs_test.cpp)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(rhscpu_test PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(rhscpu_test PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories(rhscpu_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
target_include_directories(rhscpu_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_link_libraries(rhscpu_test bssn_lib dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

#add_library(had_rhs STATIC src/hadrhs.F90)

set(SOURCE_FILES src/bssngr_main.cpp include/gr.h src/bssnCtx.cpp ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(bssnSolver ${SOURCE_FILES})
target_include_directories(bssnSolver PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(bssnSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(bssnSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(bssnSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(bssnSolver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(bssnSolver PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(bssnSolver bssn_cuda)
    target_link_libraries(bssnSolver ${CUDA_LIBRARIES})
endif()
target_link_libraries(bssnSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

add_executable(bssnSolverCtx src/bssngr_main.cpp src/bssnCtx.cpp)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(bssnSolverCtx PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(bssnSolverCtx PRIVATE ${GSL_INCLUDE_DIRS})
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda/include)
target_include_directories(bssnSolverCtx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(bssnSolverCtx PRIVATE ${CUDA_INCLUDE_DIRS})
target_link_libraries(bssnSolverCtx bssn_lib dendro5 ${CUDA_LIBRARIES})
target_link_libraries(bssnSolverCtx ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

# set(SOURCE_FILES src/gr_scaling.cpp include/gr.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(bssnScalingTest ${SOURCE_FILES})
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
# target_include_directories(bssnScalingTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
# target_include_directories(bssnScalingTest PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(bssnScalingTest PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#     target_include_directories(bssnScalingTest PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(bssnScalingTest bssn_cuda)
#     target_link_libraries(bssnScalingTest ${CUDA_LIBRARIES})
# endif()
# target_link_libraries(bssnScalingTest dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)



set(SOURCE_FILES src/gr_wamr_test.cpp include/gr.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(bssnWAMRTest ${SOURCE_FILES})
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(bssnWAMRTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(bssnWAMRTest PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(bssnWAMRTest PRIVATE ${GSL_INCLUDE_DIRS})

if(WITH_CUDA)
    target_include_directories(bssnWAMRTest PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(bssnWAMRTest bssn_cuda)
    target_link_libraries(bssnWAMRTest ${CUDA_LIBRARIES})
endif()
target_link_libraries(bssnWAMRTest dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)



# set(SOURCE_FILES src/gr_nuts.cpp include/gr.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(bssnSolverNUTS ${SOURCE_FILES})

# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
# target_include_directories(bssnSolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
# target_include_directories(bssnSolverNUTS PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(bssnSolverNUTS PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#     target_include_directories(bssnSolverNUTS PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(bssnSolverNUTS bssn_cuda)
#     target_link_libraries(bssnSolverNUTS ${CUDA_LIBRARIES})
# endif()

# target_link_libraries(bssnSolverNUTS dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


set(SOURCE_FILES src/checkDerivatives.cpp ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(checkDerivs ${SOURCE_FILES})
target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(checkDerivs PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(checkDerivs PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(checkDerivs PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(checkDerivs PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(checkDerivs bssn_cuda)
    target_link_libraries(checkDerivs ${CUDA_LIBRARIES})
endif()
target_link_libraries(checkDerivs dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

# if(WITH_CUDA)
#     set(SOURCE_FILES cuda_gr/src/gpuBSSNExample.cpp ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
#     add_executable(cudaBSSN ${SOURCE_FILES})
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/cuda)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
#     target_include_directories(cudaBSSN PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
#     target_include_directories(cudaBSSN PRIVATE ${MPI_INCLUDE_PATH})
#     target_include_directories(cudaBSSN PRIVATE ${GSL_INCLUDE_DIRS})
#     target_include_directories(cudaBSSN PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(cudaBSSN bssn_cuda)
#     target_link_libraries(cudaBSSN ${CUDA_LIBRARIES})
#     target_link_libraries(cudaBSSN dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES}  m)
#     add_dependencies(cudaBSSN BSSN_CODEGEN)
#     #add_dependencies(cudaBSSN GW_EXTRACT)
# endif()



# if(WITH_CUDA)
# set(SOURCE_FILES cuda_gr/src/gpuTestBed.cpp cuda_gr/include/gpuTestBed.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(gpuTestBed ${SOURCE_FILES})
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/cuda)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
#     target_include_directories(gpuTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
#     target_include_directories(gpuTestBed PRIVATE ${MPI_INCLUDE_PATH})
#     target_include_directories(gpuTestBed PRIVATE ${GSL_INCLUDE_DIRS})
#     target_include_directories(gpuTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(gpuTestBed bssn_cuda)
#     target_link_libraries(gpuTestBed ${CUDA_LIBRARIES})
#     target_link_libraries(gpuTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)
#     add_dependencies(gpuTestBed BSSN_CODEGEN)
#     #add_dependencies(gpuTestBed GW_EXTRACT)
# endif()


# set(SOURCE_FILES src/rhsTestBed.cpp include/rhsTestBed.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
# add_executable(rhsTestBed ${SOURCE_FILES})
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cuda_gr/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/cuda)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
# target_include_directories(rhsTestBed PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
# target_include_directories(rhsTestBed PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(rhsTestBed PRIVATE ${GSL_INCLUDE_DIRS})

# if(WITH_CUDA)
#     target_include_directories(rhsTestBed PRIVATE ${CUDA_INCLUDE_DIRS})
#     target_link_libraries(rhsTestBed bssn_cuda)
#     target_link_libraries(rhsTestBed ${CUDA_LIBRARIES})
# endif()

# target_link_libraries(rhsTestBed dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


set(SOURCE_FILES src/zipUnzipTest.cpp include/zipUnzipTest.h ${BSSN_INC} ${BSSN_SRC} ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP})
add_executable(zipUnzipTest ${SOURCE_FILES})
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/test/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(zipUnzipTest PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(zipUnzipTest PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(zipUnzipTest PRIVATE ${GSL_INCLUDE_DIRS})

if(WITH_CUDA)
    target_include_directories(zipUnzipTest PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(zipUnzipTest bssn_cuda)
    target_link_libraries(zipUnzipTest ${CUDA_LIBRARIES})
endif()

target_link_libraries(zipUnzipTest dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


if(BSSN_EXTRACT_GRAVITATIONAL_WAVES)
    add_definitions(-DBSSN_EXTRACT_GRAVITATIONAL_WAVES)
    #add_dependencies(bssnSolver GW_EXTRACT)
    #add_dependencies(checkDerivs GW_EXTRACT)
    #add_dependencies(rhsTestBed GW_EXTRACT)
    #add_dependencies(zipUnzipTest GW_EXTRACT)
endif(BSSN_EXTRACT_GRAVITATIONAL_WAVES)


add_dependencies(bssnSolver BSSN_CODEGEN)
#add_dependencies(rhsTestBed BSSN_CODEGEN)


set(SOURCE_FILES    src/TPID.cpp 
                    ${BSSN_INC}
                    ${BSSN_SRC}
                    ${CUDA_GR_INCLUDE_H}
                    ${CUDA_GR_SRC_CPP}
                    )

add_executable(tpid ${SOURCE_FILES})
target_include_directories(tpid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(tpid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(tpid PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(tpid PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(tpid PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(tpid bssn_cuda)
    target_link_libraries(tpid ${CUDA_LIBRARIES})
endif()
target_link_libraries(tpid dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)


