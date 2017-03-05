cmake_minimum_required(VERSION 2.8)

PROJECT(FindGeodesicPath)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(FindGeodesicPath FindGeodesicPath.cxx
  vtkFindGeodesicPath.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(FindGeodesicPath ${VTK_LIBRARIES})
else()
  target_link_libraries(FindGeodesicPath vtkHybrid)
endif()

