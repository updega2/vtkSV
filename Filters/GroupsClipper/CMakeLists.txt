cmake_minimum_required(VERSION 2.8)

PROJECT(GroupsClipper)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(GroupsClipper GroupsClipper.cxx
  vtksvPolyDataCenterlineGroupsClipper.cxx
  vtksvPolyBallLine.cxx)

if(VTK_LIBRARIES)
  target_link_libraries(GroupsClipper ${VTK_LIBRARIES})
else()
  target_link_libraries(GroupsClipper vtkHybrid)
endif()

