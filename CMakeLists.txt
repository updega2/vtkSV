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
include(vtkSVFindVTK)
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
set(VTKSV_INCLUDE_DIRS "")
foreach(vtk-module ${VTKSV_BUILD_DIRS})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${vtk-module})
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/${vtk-module})
  set(VTKSV_INCLUDE_DIRS ${VTKSV_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/${vtk-module} ${CMAKE_CURRENT_BINARY_DIR}/${vtk-module})
  vtk_add_module(
    "${CMAKE_CURRENT_SOURCE_DIR}/${vtk-module}"
    module.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/${vtk-module}"
    ${_test_languages})
endforeach()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add each subdirectory
set(VTKSV_LIBRARIES "")
message(STATUS "vtkSV Modules:")
foreach(vtk-module ${VTK_MODULES_ALL})
  message(STATUS "     ${vtk-module}")
  add_subdirectory(${${vtk-module}_SOURCE_DIR} ${${vtk-module}_BINARY_DIR})
  set(VTK_MODULES_ENABLED ${VTK_MODULES_ENABLED} ${vtk-module})
  set(VTKSV_LIBRARIES ${VTKSV_LIBRARIES} ${vtk-module})
endforeach()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Add external data target
if(BUILD_TESTING)
  ExternalData_add_target(VTKSVData)
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Export targets
export(TARGETS ${VTK_MODULES_ALL}
  FILE "${VTKSV_BINARY_DIR}/vtkSVTargets.cmake")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ... for the build tree
configure_file("${VTKSV_SOURCE_DIR}/CMake/vtkSVConfig.cmake.in"
  "${VTKSV_BINARY_DIR}/vtkSVConfig.cmake" @ONLY)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ... for the install tree
configure_file("${VTKSV_SOURCE_DIR}/CMake/vtkSVConfig.cmake.in"
  "${VTKSV_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vtkSVConfig.cmake" @ONLY)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ... for both
configure_file("${VTKSV_SOURCE_DIR}/CMake/vtkSVConfigVersion.cmake.in"
 "${VTKSV_BINARY_DIR}/vtkSVConfigVersion.cmake" @ONLY)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Install the vtkSVConfig.cmake and vtkSVConfigVersion.cmake
install(FILES
  "${VTKSV_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/vtkSVConfig.cmake"
  "${VTKSV_BINARY_DIR}/vtkSVConfigVersion.cmake"
  DESTINATION "${VTK_INSTALL_PACKAGE_DIR}" COMPONENT dev)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Install the export set for use with the install-tree
install(EXPORT vtkSVTargets DESTINATION
  "${VTK_INSTALL_PACKAGE_DIR}" COMPONENT dev)
#-----------------------------------------------------------------------------
