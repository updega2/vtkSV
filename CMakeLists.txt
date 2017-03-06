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
# Find VTK, only major dependency
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${VTK_CMAKE_DIR}")
include(vtkExternalModuleMacros)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# vtkSV Includes
include(vtkSVSystemCheck)
include(vtkSVSetup)
include(vtkSVOptions)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Set build dirs
if(VTKSV_BUILD_FILTERS)
  set(VTKSV_FILTER_DIRS ${VTKSV_FILTER_DIRS} "Filters/FindGeodesicPath")
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Include for headers
foreach(dir ${VTKSV_FILTER_DIRS})
  include_directories(${CMAKE_SOURCE_DIR}/${dir})
  include_directories(${CMAKE_BINARY_DIR}/${dir})
endforeach()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add each subdirectory
message(STATUS "Filters:")
foreach(dir ${VTKSV_FILTER_DIRS})
  message(STATUS "  ${dir}")
  add_subdirectory(${CMAKE_SOURCE_DIR}/${dir})
endforeach()
#-----------------------------------------------------------------------------

