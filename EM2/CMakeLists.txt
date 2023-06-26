project(EM2)

# only EM2 options
option(EM2_COMPUTE_CONSTRAINTS "Enable constraint equation computation" ON)
option(EM2_USE_4TH_ORDER_DERIVS "4th order derivs" OFF)
option(EM2_USE_6TH_ORDER_DERIVS "6th order derivs" OFF)
option(EM2_USE_8TH_ORDER_DERIVS "8th order derivs" ON)

if(EM2_USE_4TH_ORDER_DERIVS)
    set(EM2_USE_6TH_ORDER_DERIVS OFF)
    set(EM2_USE_8TH_ORDER_DERIVS OFF)
elseif(EM2_USE_6TH_ORDER_DERIVS)
    set(EM2_USE_4TH_ORDER_DERIVS OFF)
    set(EM2_USE_8TH_ORDER_DERIVS OFF)
elseif(EM2_USE_8TH_ORDER_DERIVS)
    set(EM2_USE_4TH_ORDER_DERIVS OFF)
    set(EM2_USE_6TH_ORDER_DERIVS OFF)
endif()

if(EM2_USE_4TH_ORDER_DERIVS)
    add_definitions(-DEM2_USE_4TH_ORDER_DERIVS)
endif()

if(EM2_USE_6TH_ORDER_DERIVS)
    add_definitions(-DEM2_USE_6TH_ORDER_DERIVS)
endif()

if(EM2_USE_8TH_ORDER_DERIVS)
    add_definitions(-DEM2_USE_8TH_ORDER_DERIVS)
endif()

#add_definitions(-DDEBUG_RK_SOLVER) 

if(EM2_COMPUTE_CONSTRAINTS)
    add_definitions(-DEM2_COMPUTE_CONSTRAINTS)
endif() 

set(EM2_INC include/parameters.h
        include/rk4em2.h
        include/em2Utils.h
        include/rhs.h
        include/derivs.h
        include/profile_params.h
        include/physcon.h 
        include/em2Ctx.h
        )

set(EM2_SRC src/rk4em2.cpp
        src/parameters.cpp
        src/em2Utils.cpp
        src/rhs.cpp
        src/derivs.cpp
        src/profile_params.cpp
        src/physcon.cpp 
        #src/em2Ctx.cpp
        )

set(SOURCE_FILES src/em2.cpp  include/em2.h ${EM2_INC} ${EM2_SRC})
add_executable(em2Solver ${SOURCE_FILES})
target_include_directories(em2Solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(em2Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(em2Solver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(em2Solver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(em2Solver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(em2Solver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)

# set(SOURCE_FILES src/em2_lts.cpp  include/em2.h ${EM2_INC} ${EM2_SRC})
# add_executable(em2SolverNUTS ${SOURCE_FILES})
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
# target_include_directories(em2SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
# target_include_directories(em2SolverNUTS PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(em2SolverNUTS PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#     target_include_directories(em2SolverNUTS PRIVATE ${CUDA_INCLUDE_DIRS})
# endif()
# target_link_libraries(em2SolverNUTS dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
