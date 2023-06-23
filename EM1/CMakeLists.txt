project(EM1)

# only EM1 options
option(EM1_COMPUTE_CONSTRAINTS "Enable constraint equation computation" ON)
option(EM1_USE_4TH_ORDER_DERIVS "4th order derivs" OFF)
option(EM1_USE_6TH_ORDER_DERIVS "6th order derivs" OFF)
option(EM1_USE_8TH_ORDER_DERIVS "8th order derivs" ON)

if(EM1_USE_4TH_ORDER_DERIVS)
    set(EM1_USE_6TH_ORDER_DERIVS OFF)
    set(EM1_USE_8TH_ORDER_DERIVS OFF)
elseif(EM1_USE_6TH_ORDER_DERIVS)
    set(EM1_USE_4TH_ORDER_DERIVS OFF)
    set(EM1_USE_8TH_ORDER_DERIVS OFF)
elseif(EM1_USE_8TH_ORDER_DERIVS)
    set(EM1_USE_4TH_ORDER_DERIVS OFF)
    set(EM1_USE_6TH_ORDER_DERIVS OFF)
endif()

if(EM1_USE_4TH_ORDER_DERIVS)
    add_definitions(-DEM1_USE_4TH_ORDER_DERIVS)
endif()

if(EM1_USE_6TH_ORDER_DERIVS)
    add_definitions(-DEM1_USE_6TH_ORDER_DERIVS)
endif()

if(EM1_USE_8TH_ORDER_DERIVS)
    add_definitions(-DEM1_USE_8TH_ORDER_DERIVS)
endif()
#add_definitions(-DDEBUG_RK_SOLVER) 

if(EM1_COMPUTE_CONSTRAINTS)
    add_definitions(-DEM1_COMPUTE_CONSTRAINTS)
endif() 

set(EM1_INC include/parameters.h
        include/rk4em1.h
        include/em1Utils.h
        include/rhs.h
        include/derivs.h
        include/profile_params.h
        include/physcon.h 
        include/em1Ctx.h
        )

set(EM1_SRC src/rk4em1.cpp
        src/parameters.cpp
        src/em1Utils.cpp
        src/rhs.cpp
        src/derivs.cpp
        src/profile_params.cpp
        src/physcon.cpp 
        #src/em1Ctx.cpp
        )

set(SOURCE_FILES src/em1.cpp  include/em1.h ${EM1_INC} ${EM1_SRC})
add_executable(em1Solver ${SOURCE_FILES})
target_include_directories(em1Solver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(em1Solver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(em1Solver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(em1Solver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
    target_include_directories(em1Solver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(em1Solver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)



# set(SOURCE_FILES src/em1_lts.cpp  include/em1.h ${EM1_INC} ${EM1_SRC})
# add_executable(em1SolverNUTS ${SOURCE_FILES})
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
# target_include_directories(em1SolverNUTS PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
# target_include_directories(em1SolverNUTS PRIVATE ${MPI_INCLUDE_PATH})
# target_include_directories(em1SolverNUTS PRIVATE ${GSL_INCLUDE_DIRS})
# if(WITH_CUDA)
#     target_include_directories(em1SolverNUTS PRIVATE ${CUDA_INCLUDE_DIRS})
# endif()
# target_link_libraries(em1SolverNUTS dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
