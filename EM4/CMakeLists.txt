project(EM4)

# only EM4 options
option(EM4_COMPUTE_CONSTRAINTS "Enable constraint equation computation" ON)
option(EM4_USE_4TH_ORDER_DERIVS "4th order derivs" OFF)
option(EM4_USE_6TH_ORDER_DERIVS "6th order derivs" OFF)
option(EM4_USE_8TH_ORDER_DERIVS "8th order derivs" ON)
#add_definitions(-DDEBUG_RK_SOLVER) 

if(EM4_USE_4TH_ORDER_DERIVS)
    set(EM4_USE_6TH_ORDER_DERIVS OFF)
    set(EM4_USE_8TH_ORDER_DERIVS OFF)
elseif(EM4_USE_6TH_ORDER_DERIVS)
    set(EM4_USE_4TH_ORDER_DERIVS OFF)
    set(EM4_USE_8TH_ORDER_DERIVS OFF)
elseif(EM4_USE_8TH_ORDER_DERIVS)
    set(EM4_USE_4TH_ORDER_DERIVS OFF)
    set(EM4_USE_6TH_ORDER_DERIVS OFF)
endif()

if(EM4_USE_4TH_ORDER_DERIVS)
    add_definitions(-DEM4_USE_4TH_ORDER_DERIVS) 
endif()

if(EM4_USE_6TH_ORDER_DERIVS)
    add_definitions(-DEM4_USE_6TH_ORDER_DERIVS) 
endif()

if(EM4_USE_8TH_ORDER_DERIVS)
    add_definitions(-DEM4_USE_8TH_ORDER_DERIVS) 
endif()


if(EM4_COMPUTE_CONSTRAINTS)
    add_definitions(-DEM4_COMPUTE_CONSTRAINTS)
endif() 

set(EM4_INC include/parameters.h
        include/rk4em4.h
        include/em4Utils.h
        include/rhs.h
        include/derivs.h
        include/profile_params.h
        include/physcon.h 
        include/em4Ctx.h
        )

set(EM4_SRC src/rk4em4.cpp
        src/parameters.cpp
        src/em4Utils.cpp
        src/rhs.cpp
        src/derivs.cpp
        src/profile_params.cpp
        src/physcon.cpp
        #src/em4Ctx.cpp
        )

set(SOURCE_FILES src/em4.cpp  include/em4.h ${EM4_INC} ${EM4_SRC})
add_executable(em4Solver ${SOURCE_FILES})
target_include_directories(em4Solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(em4Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(em4Solver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(em4Solver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(em4Solver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(em4Solver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)


# set(SOURCE_FILES src/em4_lts.cpp  include/em4.h ${EM4_INC} ${EM4_SRC})
# add_executable(em4SolverNUTS ${SOURCE_FILES})
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
# target_include_directories(em4SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
# target_include_directories(em4SolverNUTS PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(em4SolverNUTS PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#     target_include_directories(em4SolverNUTS PRIVATE ${CUDA_INCLUDE_DIRS})
# endif()
# target_link_libraries(em4SolverNUTS dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
