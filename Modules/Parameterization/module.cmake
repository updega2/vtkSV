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

set(DOCUMENTATION "A module containing code to parameterize and arbitrary input vasular surface.")

#------------------------------------------------------------------------------
# NURBS addition
set(EXTRA_TEST_DEPENDS "")
if(vtkRenderingFreeType${VTK_RENDERING_BACKEND}_LOADED)
  set(EXTRA_TEST_DEPENDS ${EXTRA_TEST_DEPENDS} vtkRenderingFreeType${VTK_RENDERING_BACKEND})
endif()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Module depends
vtk_module(vtkSVParameterization
  DESCRIPTION
  "${DOCUMENTATION}"
  DEPENDS
  vtkCommonDataModel
  vtkCommonExecutionModel
  vtkFiltersCore
  vtkFiltersGeometry
  vtkFiltersModeling
  vtkSVCommon
  vtkSVIO
  vtkSVMisc
  vtkSVNURBS
  vtkSVGeometry
  ${EXTRA_DEPENDS}
  TEST_DEPENDS
  vtkFiltersExtraction
  vtkTestingCore
  vtkIOGeometry
  vtkIOXML
  vtkSVCommon
  vtkSVIO
  ${EXTRA_TEST_DEPENDS}
  TCL_NAME
  vtkSVParameterization
  )
#------------------------------------------------------------------------------
