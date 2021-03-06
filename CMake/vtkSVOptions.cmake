# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
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
# Build regular libraries instead of typical build system
option(VTKSV_BUILD_LIBS_AS_VTK_MODULES "Option to build the vtkSV libs as vtk modules" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Enable Testing
option(BUILD_TESTING       "Build ${PROJECT_NAME} testing" OFF)
option(VTKSV_TEST_COVERAGE "Add test coverage"             OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Enable Testing
option(VTKSV_BUILD_DOCUMENTATION "Build ${PROJECT_NAME} documentation" OFF)
set(VTKSV_DOCS_DIR ${VTKSV_BINARY_DIR}/Doxygen CACHE PATH "Location to put doxygen files")
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Libs options
option(BUILD_SHARED_LIBS "Build ${PROJECT_NAME} as shared libraries." OFF)

set(SV_LIBRARY_TYPE "STATIC" CACHE STRING "Options are STATIC or SHARED" FORCE)
set_property(CACHE SV_LIBRARY_TYPE PROPERTY STRINGS STATIC SHARED)
mark_as_advanced(SV_LIBRARY_TYPE)
if(BUILD_SHARED_LIBS)
	set(SV_LIBRARY_TYPE "SHARED" CACHE STRING "Shared cache" FORCE)
else()
  set(SV_LIBRARY_TYPE "STATIC" CACHE STRING "Static cache" FORCE)
endif()
#----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Specify which thirdparty libs to build
option(VTKSV_BUILD_THIRDPARTY_VMTK "Option to build the ThirdParty VMTK code" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Specify which filters to build
option(VTKSV_BUILD_EXES "Option to build the executables for each filter" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Specify which modules to build
option(VTKSV_BUILD_MODULE_MISC             "Option to build the Miscellaneous code"    ON)
option(VTKSV_BUILD_MODULE_NURBS            "Option to build the NURBS code"            ON)
option(VTKSV_BUILD_MODULE_BOOLEAN          "Option to build the Boolean code"          ON)
option(VTKSV_BUILD_MODULE_PARAMETERIZATION "Option to build the Parameterization code" ON)
option(VTKSV_BUILD_MODULE_SEGMENTATION     "Option to build the Segmentation code"     ON)
option(VTKSV_BUILD_MODULE_GEOMETRY         "Option to build the Geometry code"         ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Specify to install libs/headers
option(SV_INSTALL_HEADERS "Option to install vtkSV headers" ON)
option(SV_INSTALL_LIBS    "Option to install vtkSV libs"    ON)
