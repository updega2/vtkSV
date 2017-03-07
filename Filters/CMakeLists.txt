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
set(SRCS vtkSVFindGeodesicPath.cxx vtkSVFindSeparateRegions.cxx vtkSVGetBoundaryFaces.cxx
  vtkSVGroupsClipper.cxx vtkSVPolyBallLine.cxx vtkSVHausdorffDistance.cxx
  vtkSVMapInterpolator.cxx vtkSVPassDataArray.cxx vtkSVPullApartPolyData.cxx)
set(HDRS vtkSVFindGeodesicPath.h   vtkSVFindSeparateRegions.h   vtkSVGetBoundaryFaces.h
  vtkSVGroupsClipper.h   vtkSVPolyBallLine.h   vtkSVHausdorffDistance.h
  vtkSVMapInterpolator.h   vtkSVPassDataArray.cxx vtkSVPullApartPolyData.h)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Filters module
vtksv_add_module(vtkSVFilters
  SRCS ${SRCS}
  HDRS ${HDRS}
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Executables
if(VTKSV_BUILD_FILTER_EXES)
  set(EXE_LIST FindGeodesicPath GroupsClipper HausdorffDistance
    PassDataArray PullApartPolyData)
  foreach(exe ${EXE_LIST})
    add_executable(${exe} ${exe}.cxx)
    target_link_libraries(${exe} ${VTK_LIBRARIES} vtkSVFilters vtkSVCommon)
    install(TARGETS ${exe}
      RUNTIME DESTINATION ${VTKSV_INSTALL_RUNTIME_DIR} COMPONENT FilterExecutables)
  endforeach()
endif()
#------------------------------------------------------------------------------
