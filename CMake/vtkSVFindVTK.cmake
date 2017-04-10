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
# Find VTK, only major dependency
if(NOT VTK_LIBRARIES)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# If we want to build vtk modules, there is some prep work
if(VTKSV_BUILD_LIBS_AS_VTK_MODULES)
  set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${VTK_CMAKE_DIR}")
  include(vtkExternalModuleMacros)
  if(BUILD_TESTING)
    set(_test_languages "Cxx")
    include(ExternalData)
    include(vtkTestingMacros RESULT_VARIABLE VTK_TESTING_MACROS)
    if("${VTK_TESTING_MACROS}" STREQUAL "NOTFOUND")
      message(FATAL_ERROR "${VTK_TESTING_FOUND} Must inlcude VTK from the Build directory to enable testing")
    endif()
    if(NOT vtkTestingCore_LOADED OR NOT vtkTestingRendering_LOADED)
      message(FATAL_ERROR "Must use a vtk that has testing enabled")
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
  # For TCL Wrapping
  set(VTKSV_NEEDED_TCL_MODULES vtkCommonDataModel vtkCommonExecutionModel)
  foreach(_vtk_module ${VTKSV_NEEDED_TCL_MODULES})
    if (NOT ${_vtk_module}_EXCLUDE_FROM_WRAPPING)
      set(${_vtk_module}_EXCLUDE_FROM_WRAPPING 0)
    endif()
    if (NOT ${vtk_module}_TCL_NAME)
      set(${_vtk_module}_TCL_NAME ${_vtk_module})
    endif()
  endforeach()
  #-----------------------------------------------------------------------------

  #-----------------------------------------------------------------------------
  # Library suffix
  set(VTK_CUSTOM_LIBRARY_SUFFIX "-${VTKSV_VERSION}")
  #-----------------------------------------------------------------------------
endif()
#-----------------------------------------------------------------------------
