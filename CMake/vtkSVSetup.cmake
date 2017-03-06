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

#------------------------------------------------------------------------------
# Setup install directories (we use names with VTK_ prefix, since vtkSV
# is built as a custom "VTK" library.
set(VTKSV_VERSION "1.0")
if(NOT VTK_INSTALL_RUNTIME_DIR)
  set(VTK_INSTALL_RUNTIME_DIR bin)
endif()
if(NOT VTK_INSTALL_LIBRARY_DIR)
  set(VTK_INSTALL_LIBRARY_DIR lib/vtksv-${VTKSV_VERSION})
endif()
if(NOT VTK_INSTALL_PYTHON_MODULE_DIR)
  set (VTK_INSTALL_PYTHON_MODULE_DIR "${VTK_INSTALL_LIBRARY_DIR}/site-packages" CACHE
    INTERNAL "Directory where python modules will be installed")
endif()
if(NOT VTK_BUILD_PYTHON_MODULE_DIR)
  set (VTK_BUILD_PYTHON_MODULE_DIR "${CMAKE_BINARY_DIR}/lib/site-packages" CACHE
    INTERNAL "Directory where python modules will be built")
endif()
if(NOT VTK_INSTALL_ARCHIVE_DIR)
  set(VTK_INSTALL_ARCHIVE_DIR lib/vtksv-${VTKSV_VERSION})
endif()
if(NOT VTK_INSTALL_INCLUDE_DIR)
  set(VTK_INSTALL_INCLUDE_DIR include/vtksv-${VTKSV_VERSION})
endif()
if(NOT VTK_INSTALL_DATA_DIR)
  set(VTK_INSTALL_DATA_DIR share/vtksv-${VTKSV_VERSION})
endif()
if(NOT VTK_INSTALL_DOC_DIR)
  set(VTK_INSTALL_DOC_DIR share/doc/vtksv-${VTKSV_VERSION})
endif()
if(NOT VTK_INSTALL_PACKAGE_DIR)
  set(VTK_INSTALL_PACKAGE_DIR "lib/cmake/vtksv-${VTKSV_VERSION}")
endif()
if(NOT VTK_INSTALL_DOXYGEN_DIR)
  set(VTK_INSTALL_DOXYGEN_DIR ${VTK_INSTALL_DOC_DIR}/doxygen)
endif()
if(NOT VTK_INSTALL_EXPORT_NAME)
  set(VTK_INSTALL_EXPORT_NAME vtkSVTargets)
endif()
if(NOT VTK_MODULES_DIR)
  set(VTK_MODULES_DIR "${VTKSV_BINARY_DIR}/${VTK_INSTALL_PACKAGE_DIR}/Modules")
endif()
set(VTKSV_MODULES_DIR ${VTK_MODULES_DIR})

# Handle the target export file, this is used if building against a build tree.
if(NOT VTK_EXPORTS_FILE)
  set(VTK_EXPORTS_FILE "${vtkSV_BINARY_DIR}/VTK/${VTK_INSTALL_EXPORT_NAME}.cmake")
endif()
file(REMOVE "${VTK_EXPORTS_FILE}")
#------------------------------------------------------------------------------
