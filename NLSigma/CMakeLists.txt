project(NLSM)

# options for nlsm
option(NLSM_NONLINEAR "enable nonlinear version" OFF)
option(NLSM_COMPARE_WITH_ANALYTICAL_SOL
       "compare with the analytical sol. of the linear wave equation" OFF)
option(NLSM_USE_4TH_ORDER_DERIVS "4th order derivs" OFF)
option(NLSM_USE_6TH_ORDER_DERIVS "6th order derivs" ON)
option(NLSM_USE_8TH_ORDER_DERIVS "8th order derivs" OFF)

if(NLSM_USE_4TH_ORDER_DERIVS)
  set(NLSM_USE_6TH_ORDER_DERIVS OFF)
  set(NLSM_USE_8TH_ORDER_DERIVS OFF)
elseif(NLSM_USE_6TH_ORDER_DERIVS)
  set(NLSM_USE_4TH_ORDER_DERIVS OFF)
  set(NLSM_USE_8TH_ORDER_DERIVS OFF)
elseif(NLSM_USE_8TH_ORDER_DERIVS)
  set(NLSM_USE_4TH_ORDER_DERIVS OFF)
  set(NLSM_USE_6TH_ORDER_DERIVS OFF)
endif()

if(NLSM_USE_4TH_ORDER_DERIVS)
  add_definitions(-DNLSM_USE_4TH_ORDER_DERIVS)
endif()

if(NLSM_USE_6TH_ORDER_DERIVS)
  add_definitions(-DNLSM_USE_6TH_ORDER_DERIVS)
endif()

if(NLSM_USE_8TH_ORDER_DERIVS)
  add_definitions(-DNLSM_USE_8TH_ORDER_DERIVS)
endif()

if(NLSM_COMPARE_WITH_ANALYTICAL_SOL)
  set(NLSM_NONLINEAR OFF)
  add_definitions(-DNLSM_COMPARE_WITH_ANALYTICAL_SOL)
endif()

if(NLSM_NONLINEAR)
  add_definitions(-DNLSM_NONLINEAR)
endif()

# set(NLSM_INC include/parameters.h include/rk4nlsm.h include/nlsmUtils.h
# include/rhs.h include/derivs.h include/profile_params.h include/nlsmCtx.h
# include/nlsmInvCtx.h )

set(NLSM_SRC src/parameters.cpp src/nlsmUtils.cpp src/rhs.cpp src/derivs.cpp
             src/profile_params.cpp)

set(SOURCE_FILES src/nlsm.cpp src/rk4nlsm.cpp ${NLSM_SRC})
add_executable(nlsmSolver ${SOURCE_FILES})
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(nlsmSolver
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(nlsmSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(nlsmSolver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
  target_include_directories(nlsmSolver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(nlsmSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

set(SOURCE_FILES src/nlsm_nuts.cpp src/nlsmCtx.cpp ${NLSM_SRC})
add_executable(nlsmSolverNUTS ${SOURCE_FILES})
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(nlsmSolverNUTS
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(nlsmSolverNUTS PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(nlsmSolverNUTS PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
  target_include_directories(nlsmSolverNUTS PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(nlsmSolverNUTS dendro5 ${LAPACK_LIBRARIES}
                      ${MPI_LIBRARIES} m)

if(WITH_CUDA)
  set(SOURCE_FILES src/nlsmCtxGPU.cu src/nlsm_cuda.cu ${NLSM_SRC})
  cuda_add_executable(nlsmSolverCUDA ${SOURCE_FILES})
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/GPU/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
  target_include_directories(
    nlsmSolverCUDA PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
  target_include_directories(nlsmSolverCUDA
                             PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
  target_include_directories(nlsmSolverCUDA PRIVATE ${MPI_INCLUDE_PATH})
  target_include_directories(nlsmSolverCUDA PRIVATE ${GSL_INCLUDE_DIRS})
  target_include_directories(nlsmSolverCUDA PRIVATE ${CUDA_INCLUDE_DIRS})
  target_link_libraries(nlsmSolverCUDA dendro5_cuda ${CUDA_LIBRARIES}
                        ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
endif()

set(SOURCE_FILES src/nlsm_sdc.cpp src/nlsmCtx.cpp ${NLSM_SRC})
add_executable(nlsmSolverSDC ${SOURCE_FILES})
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(nlsmSolverSDC
                           PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(nlsmSolverSDC PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(nlsmSolverSDC PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
  target_include_directories(nlsmSolverSDC PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(nlsmSolverSDC dendro5 ${LAPACK_LIBRARIES}
                      ${MPI_LIBRARIES} m)

# set(SOURCE_FILES src/nlsm_inv.cpp include/nlsm.h ${NLSM_INC} ${NLSM_SRC})
# add_executable(nlsmInvSolver ${SOURCE_FILES})
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
# target_include_directories(nlsmInvSolver
#                            PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
# target_include_directories(nlsmInvSolver PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(nlsmInvSolver PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#   target_include_directories(nlsmInvSolver PRIVATE ${CUDA_INCLUDE_DIRS})
# endif()
# target_link_libraries(nlsmInvSolver dendro5 ${LAPACK_LIBRARIES}
#                       ${MPI_LIBRARIES} m)
