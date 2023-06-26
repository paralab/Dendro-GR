project(Fluid)

set(FLUID_INC include/parameters.h
        include/rkfluid.h
              include/fluidUtils.h
              include/rhs.h
              include/derivs.h
              include/profile_params.h
              include/average.h
              include/geom.h
              include/recon.h
              include/rhshlle.h
              include/root.h
              include/fluidMath.h
   )

set(FLUID_SRC src/rkfluid.cpp
              src/parameters.cpp
              src/fluidUtils.cpp
              src/rhs.cpp
              src/derivs.cpp
              src/profile_params.cpp
							src/fluid.cpp
							src/recon.cpp
							src/rhshlle.cpp
              src/fluidMath.cpp
              src/root.cpp
   )

set(SOURCE_FILES src/fluid.cpp  include/fluid.h ${FLUID_INC} ${FLUID_SRC})
add_executable(fluidSolver ${SOURCE_FILES})
target_include_directories(fluidSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/include/test)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/examples/include)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/FEM/include)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/ODE/include)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/LinAlg/include)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/vtk/include)
target_include_directories(fluidSolver PRIVATE ${CMAKE_SOURCE_DIR}/dendrolib/IO/zlib/inc)
target_include_directories(fluidSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(fluidSolver PRIVATE ${GSL_INCLUDE_DIRS})
if(WITH_CUDA)
  target_include_directories(fluidSolver PRIVATE ${CUDA_INCLUDE_DIRS})
endif()
target_link_libraries(fluidSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} m)
