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
# FindGeodesicPath
vtksv_add_module(vtkSVFindGeodesicPath
  INSTALL_EXES
  EXECUTABLE_NAME FindGeodesicPath
  SRCS vtkSVFindGeodesicPath.cxx
  HDRS vtkSVFindGeodesicPath.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# FindSeparateRegions
vtksv_add_module(vtkSVFindSeparateRegions
  SRCS vtkSVFindSeparateRegions.cxx
  HDRS vtkSVFindSeparateRegions.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# GetBoundaryFaces
vtksv_add_module(vtkSVGetBoundaryFaces
  SRCS vtkSVGetBoundaryFaces.cxx
  HDRS vtkSVGetBoundaryFaces.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# GroupClipper
vtksv_add_module(vtkSVGroupsClipper
  INSTALL_EXES
  EXECUTABLE_NAME GroupsClipper
  SRCS vtkSVPolyDataCenterlineGroupsClipper.cxx vtkSVPolyBallLine.cxx
  HDRS vtkSVPolyDataCenterlineGroupsClipper.h vtkSVPolyBallLine.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# HausdorffDistance
vtksv_add_module(vtkSVHausdorffDistance
  INSTALL_EXES
  EXECUTABLE_NAME HausdorffDistance
  SRCS vtkSVHausdorffDistance.cxx
  HDRS vtkSVHausdorffDistance.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# MapInterpolator
vtksv_add_module(vtkSVMapInterpolator
  SRCS vtkSVMapInterpolator.cxx
  HDRS vtkSVMapInterpolator.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# PassDataArray
vtksv_add_module(vtkSVPassDataArray
  INSTALL_EXES
  EXECUTABLE_NAME PassDataArray
  SRCS vtkSVPassDataArray.cxx
  HDRS vtkSVPassDataArray.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# PullApartPolyData
vtksv_add_module(vtkSVPullApartPolyData
  INSTALL_EXES
  EXECUTABLE_NAME PullApartPolyData
  SRCS vtkSVPullApartPolyData.cxx
  HDRS vtkSVPullApartPolyData.h
  PACKAGE_DEPENDS ${VTK_LIBRARIES})
#------------------------------------------------------------------------------
