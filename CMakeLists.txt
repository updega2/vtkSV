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
# Find VTK, only major dependency
if(NOT VTK_LIBRARIES)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${VTK_CMAKE_DIR}")
include(vtkExternalModuleMacros)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Set build dirs
set(VTKSV_BUILD_DIRS "")
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
foreach(dir ${VTKSV_BUILD_DIRS})
  vtk_add_to_module_search_path("${CMAKE_CURRENT_SOURCE_DIR}/${dir}" "${CMAKE_CURRENT_BINARY_DIR}/${dir}")
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${dir})
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/${dir})
endforeach()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add each subdirectory
message(STATUS "Build Dirs:")
foreach(dir ${VTKSV_BUILD_DIRS})
  message(STATUS "  ${dir}")
  add_subdirectory(${CMAKE_SOURCE_DIR}/${dir})
  get_filename_component(MODULE_DIR ${dir} NAME)
  set(VTK_MODULES_ENABLED ${VTK_MODULES_ENABLED} vtkSV${MODULE_DIR})
endforeach()
#-----------------------------------------------------------------------------
