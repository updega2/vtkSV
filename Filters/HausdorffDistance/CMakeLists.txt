cmake_minimum_required(VERSION 2.8)

PROJECT(HausdorffDistance)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(HausdorffDistance HausdorffDistance.cxx
                              vtkHausdorffDistance.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(HausdorffDistance ${VTK_LIBRARIES})
else()
  target_link_libraries(HausdorffDistance vtkHybrid)
endif()
