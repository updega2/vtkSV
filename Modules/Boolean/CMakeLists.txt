cmake_minimum_required(VERSION 2.8)

PROJECT(Boolean)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

option(PREDELAUNAY_SPLIT "Split cells manually before giving to delaunay" ON)
set(DELAUNAY_TYPE "CURRENT" CACHE STRING "Options are CURRENT, OLD, TRIANGLE")
set_property(CACHE DELAUNAY_TYPE PROPERTY STRINGS CURRENT OLD TRIANGLE)

set(CXX_FILES Boolean.cxx
  vtkLoopBooleanPolyDataFilter.cxx
  vtkLoopIntersectionPolyDataFilter.cxx)

if("${DELAUNAY_TYPE}" STREQUAL "OLD")
  set(CXX_FILES ${CXX_FILES}
    vtkDelaunay2D_60.cxx)
endif()

if("${DELAUNAY_TYPE}" STREQUAL "TRIANGLE")
# Search for Schewchuk's Triangle thirdparty library
  if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/triangle.c)
    message(FATAL_ERROR "Shewchuk's Triangle sources not found in sources directory! Please copy triangle sources in module sources directory.")
  endif()
  add_definitions(-DTRIANGLE_SAFE_PREDICATES)
  add_definitions(-DTRILIBRARY)
  set(CXX_FILES ${CXX_FILES}
    vtkTriangleDelaunay2D.cxx
    predicates.c
    triangle.c)
endif()

add_executable(Boolean ${CXX_FILES})

if(VTK_LIBRARIES)
	target_link_libraries(Boolean ${VTK_LIBRARIES})
else()
	target_link_libraries(Boolean vtkHybrid)
endif()

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/delaunay_options.h.in" "${CMAKE_CURRENT_BINARY_DIR}/delaunay_options.h")
include_directories(${CMAKE_CURRENT_BINARY_DIR})
