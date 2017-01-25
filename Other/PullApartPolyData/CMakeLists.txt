cmake_minimum_required(VERSION 2.8)

PROJECT(PullApartPolyData)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(PullApartPolyData PullApartPolyData.cxx
  vtkPullApartPolyData.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(PullApartPolyData ${VTK_LIBRARIES})
else()
  target_link_libraries(PullApartPolyData vtkHybrid)
endif()

