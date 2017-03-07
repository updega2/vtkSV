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
# Core SRCS and HDRS
set(SRCS
  vtkSVGeneralizedPolycube.cxx
  vtkSVPlanarMapper.cxx
  vtkSVBoundaryMapper.cxx
  vtkSVSquareBoundaryMapper.cxx
  vtkSVSuperSquareBoundaryMapper.cxx
  sparse_matrix.cxx
  svMath.cxx)

set(HDRS
  vtkSVGeneralizedPolycube.h
  vtkSVPlanarMapper.h
  vtkSVBoundaryMapper.h
  vtkSVSquareBoundaryMapper.h
  vtkSVSuperSquareBoundaryMapper.h
  sparse_matrix.h
  svMath.h)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Filters addition
if(VTKSV_BUILD_FILTERS)
  set(SRCS ${SRCS} vtkSVPolyDataSliceAndDiceFilter.cxx
    svGraph.cxx svGCell.cxx)
  set(HDRS ${HDRS} vtkSVPolyDataSliceAndDiceFilter.h
    svGraph.h   svGCell.h)
  # NURBS addition
  if(VTKSV_BUILD_MODULE_NURBS)
    set(SRCS ${SRCS} vtkSVPolyDataToNURBSFilter.cxx)
    set(HDRS ${HDRS} vtkSVPolyDataToNURBSFilter.h)
  endif()
endif()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# FindGeodesicPath
vtksv_add_module(vtkSVParameterization
  SRCS ${SRCS}
  HDRS ${HDRS}
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Create additional executables
if(VTKSV_BUILD_FILTER_EXES)
  set(EXTRA_LIBRARIES "")
  if(VTKSV_BUILD_FILTERS)
    set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES} vtkSVFilters)
    set(EXE_LIST ${EXE_LIST} SliceAndDice)
    if(VTKSV_BUILD_MODULE_NURBS)
      set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES} vtkSVNURBS)
      set(EXE_LIST ${EXE_LIST} PolyDataToNURBS)
    endif()
  endif()
  foreach(exe ${EXE_LIST})
    add_executable(${exe} ${exe}.cxx)
    target_link_libraries(${exe} ${VTK_LIBRARIES} vtkSVParameterization ${EXTRA_LIBRARIES})
    install(TARGETS ${exe}
      RUNTIME DESTINATION ${VTKSV_INSTALL_RUNTIME_DIR} COMPONENT ParameterizationExecutables)
  endforeach()
endif()
#------------------------------------------------------------------------------
