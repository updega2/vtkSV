# Copyright (c) 2014-2015 The Regents of the University of California.
# All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.12)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
Project(VTKSV)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# CMake code
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake"
	"${CMAKE_MODULE_PATH}")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# RPATH handling
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH 0)
endif()
# Don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
# Don't use the install RPATH
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# The RPATH that will be used when installing
set(CMAKE_INSTALL_RPATH "")
# Don't automatically add determined parts of the RPATH which point to
# directories outside the build tree to the install tree
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# vtkSV Includes
include(vtkSVSystemCheck)
include(vtkSVSetup)
include(vtkSVOptions)
include(vtkSVMacros)
if(BUILD_TESTING)
  include(vtkSVTesting)
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Setup documentation
if(VTKSV_BUILD_DOCUMENTATION)
  find_package(Doxygen REQUIRED)
  if(DOXYGEN_FOUND)
    configure_file(${VTKSV_SOURCE_DIR}/Documentation/vtksv.Doxyfile.in
      ${VTKSV_BINARY_DIR}/vtksv.Doxyfile @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${VTKSV_BINARY_DIR}/vtksv.Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
  endif()
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Find VTK, only major dependency
if(NOT VTK_LIBRARIES)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${VTK_CMAKE_DIR}")
include(vtkExternalModuleMacros)
if(BUILD_TESTING)
  set(_test_languages "Cxx")
  include(ExternalData)
  include(vtkTestingMacros RESULT_VARIABLE VTK_TESTING_MACROS)
  if ("${VTK_TESTING_MACROS}" STREQUAL "NOTFOUND")
    message(FATAL_ERROR "${VTK_TESTING_FOUND} Must inlcude VTK from the Build directory to enable testing")
  endif()
  if (NOT VTK_SOURCE_DIR)
    set(VTK_SOURCE_DIR "${VTK_CMAKE_DIR}/..")
  endif()
  if(NOT vtkTestingRendering_SOURCE_DIR)
    set(vtkTestingRendering_SOURCE_DIR "${VTK_SOURCE_DIR}/Testing/Rendering")
  endif()
else()
  set(_test_languages "")
endif()
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# Set build dirs
set(VTKSV_BUILD_DIRS "Common")
if(VTKSV_BUILD_FILTERS)
  set(VTKSV_BUILD_DIRS ${VTKSV_BUILD_DIRS} "Filters")
endif()
if(VTKSV_BUILD_MODULE_BOOLEAN)
  set(VTKSV_BUILD_DIRS ${VTKSV_BUILD_DIRS} "Modules/Boolean")
endif()
if(VTKSV_BUILD_MODULE_NURBS)
  set(VTKSV_BUILD_DIRS ${VTKSV_BUILD_DIRS} "Modules/NURBS")
endif()
if(VTKSV_BUILD_MODULE_PARAMETERIZATION)
  set(VTKSV_BUILD_DIRS ${VTKSV_BUILD_DIRS} "Modules/Parameterization")
endif()
if(VTKSV_BUILD_MODULE_GEOMETRY)
  set(VTKSV_BUILD_DIRS ${VTKSV_BUILD_DIRS} "Modules/Geometry")
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Include for headers
foreach(vtk-module ${VTKSV_BUILD_DIRS})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${vtk-module})
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/${vtk-module})
  vtk_add_module(
    "${CMAKE_CURRENT_SOURCE_DIR}/${vtk-module}"
    module.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/${vtk-module}"
    ${_test_languages})
endforeach()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add each subdirectory
message(STATUS "Modules:")
foreach(vtk-module ${VTK_MODULES_ALL})
  message(STATUS "  ${vtk-module}")
  add_subdirectory(${${vtk-module}_SOURCE_DIR} ${${vtk-module}_BINARY_DIR})
  set(VTK_MODULES_ENABLED ${VTK_MODULES_ENABLED} ${vtk-module})
endforeach()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add external data target
if(BUILD_TESTING)
  ExternalData_add_target(VTKSVData)
endif()
#-----------------------------------------------------------------------------
