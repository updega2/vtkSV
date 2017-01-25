cmake_minimum_required(VERSION 2.8)

PROJECT(PlanarMapper)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(PlanarMapper PlanarMapper.cxx
  vtkPlanarMapper.cxx
  svMath.cxx
  sparse_matrix.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(PlanarMapper ${VTK_LIBRARIES})
else()
  target_link_libraries(PlanarMapper vtkHybrid)
endif()

