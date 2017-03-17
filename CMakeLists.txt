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
# vtkSV Includes
include(vtkSVSystemCheck)
include(vtkSVSetup)
include(vtkSVOptions)
include(vtkSVMacros)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Setup testing
if(BUILD_TESTING)
  set(VTKSV_TEST_DIR                 "${CMAKE_BINARY_DIR}/Testing/Temporary")
  set(ExternalData_BINARY_ROOT       "${CMAKE_BINARY_DIR}/ExternalData")
  set(VTKSV_TEST_DATA_DIR            "${CMAKE_SOURCE_DIR}/Testing/Data")
  set(ExternalData_OBJECT_STORES     "${ExternalData_OBJECT_STORES}" "${VTKSV_TEST_DATA_DIR}")
  set(ExternalData_URL_TEMPLATES     "http://simvascular.stanford.edu/downloads/public/vtkSV/Testing/Data/%(algo)/%(hash)")
  set(VTKSV_TEST_BASELINE_DIR        "${VTKSV_TEST_DATA_DIR}/Baseline")
  set(VTKSV_TEST_OUTPUT_DATA_DIR     "${ExternalData_BINARY_ROOT}/Testing/Data")
  set(VTKSV_TEST_OUTPUT_BASELINE_DIR "${VTKSV_TEST_OUTPUT_DATA_DIR}/Baseline")
  set(VTKSV_TEST_OUTPUT_DIR          "${CMAKE_BINARY_DIR}/Testing/Temporary")
  set(VTKSV_DATA_ROOT                "${VTKSV_TEST_OUTPUT_DATA_DIR}")
  set(ExternalData_LINK_CONTENT      MD5)
  make_directory(${VTKSV_TEST_DIR})
  enable_testing()
  if(VTKSV_TEST_COVERAGE)
    set(COVERAGE_COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/llvm-cov" CACHE STRING "Forcing llvm-cov" FORCE)
    set(CTEST_COVERAGE_COMMAND "${COVERAGE_COMMAND}" CACHE STRING "Forcing llvm-cov" FORCE)
    set(COVERAGE_EXTRA_FLAGS "gcov -l" CACHE STRING "Force coverage flags" FORCE)
    set(CTEST_COVERAGE_EXTRA_FLAGS "${COVERAGE_EXTRA_FLAGS}" CACHE STRING "Force coverage flags" FORCE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -Wall -W -Wshadow -Wunused-variable \
      -Wunused-parameter -Wunused-function -Wunused -Wno-system-headers \
      -Wno-deprecated -Woverloaded-virtual -Wwrite-strings -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -Wall -W -fprofile-arcs -ftest-coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-arcs -ftest-coverage")
  endif()
  include(CTest)
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
