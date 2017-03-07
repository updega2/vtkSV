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
  vtkSVPolyDataSliceAndDiceFilter.cxx
  vtkSVPolyDataToNURBSFilter.cxx
  vtkSVPlanarMapper.cxx
  vtkSVBoundaryMapper.cxx
  vtkSVSquareBoundaryMapper.cxx
  vtkSVSuperSquareBoundaryMapper.cxx
  sparse_matrix.cxx
  svMath.cxx
  svGraph.cxx
  svGCell.cxx)

set(HDRS
  vtkSVGeneralizedPolycube.h
  vtkSVPolyDataSliceAndDiceFilter.h
  vtkSVPolyDataToNURBSFilter.h
  vtkSVPlanarMapper.h
  vtkSVBoundaryMapper.h
  vtkSVSquareBoundaryMapper.h
  vtkSVSuperSquareBoundaryMapper.h
  sparse_matrix.h
  svMath.h
  svGraph.h
  svGCell.h)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# NURBS addition
if(VTKSV_BUILD_MODULE_NURBS)
  set(SRCS ${SRCS} vtkSVPolyDataToNURBSFilter.cxx)
  set(HDRS ${HDRS} vtkSVPolyDataToNURBSFilter.h)
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
set(EXTRA_LIBRARIES vtkSVMapInterpolator vtkSVFindGeodesicPath
  vtkSVPullApartPolyData vtkSVFindSeparateRegions)
if(VTKSV_BUILD_MODULE_NURBS)
  set(EXTRA_LIBRARIES ${EXTRA_LIBRARIES} vtkSVNURBS)
endif()
set(planar_mapper_exe PlanarMapper)
add_executable(${planar_mapper_exe} ${planar_mapper_exe}.cxx ${SRCS})
target_link_libraries(${planar_mapper_exe} ${VTK_LIBRARIES} ${EXTRA_LIBRARIES})
install(TARGETS ${planar_mapper_exe}
  RUNTIME DESTINATION ${VTKSV_INSTALL_RUNTIME_DIR} COMPONENT Executables)

set(super_square_planar_mapper SuperSquarePlanarMapper)
add_executable(${super_square_planar_mapper} ${super_square_planar_mapper}.cxx ${SRCS})
target_link_libraries(${super_square_planar_mapper} ${VTK_LIBRARIES} ${EXTRA_LIBRARIES})
install(TARGETS ${super_square_planar_mapper}
  RUNTIME DESTINATION ${VTKSV_INSTALL_RUNTIME_DIR} COMPONENT Executables)

set(slicer_exe SliceAndDice)
add_executable(${slicer_exe} ${slicer_exe}.cxx ${SRCS})
target_link_libraries(${slicer_exe} ${VTK_LIBRARIES} ${EXTRA_LIBRARIES})
install(TARGETS ${slicer_exe}
  RUNTIME DESTINATION ${VTKSV_INSTALL_RUNTIME_DIR} COMPONENT Executables)

if(VTKSV_BUILD_MODULE_NURBS)
  set(converter_exe PolyDataToNURBS)
  add_executable(${converter_exe} ${converter_exe}.cxx ${SRCS})
  target_link_libraries(${converter_exe} ${VTK_LIBRARIES} ${EXTRA_LIBRARIES})
  install(TARGETS ${converter_exe}
    RUNTIME DESTINATION ${VTKSV_INSTALL_RUNTIME_DIR} COMPONENT Executables)
endif()
#------------------------------------------------------------------------------
