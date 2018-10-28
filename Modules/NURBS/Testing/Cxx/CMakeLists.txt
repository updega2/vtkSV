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

vtksv_add_test_cxx(${vtk-module}CxxTests tests
  TestLoftNURBSCurve.cxx,NO_DATA
  TestLoftNURBSSurface.cxx,NO_DATA
  TestCircle.cxx,NO_DATA
  TestCurveInsertKnot.cxx,NO_DATA
  TestCurveKnotRefinement.cxx,NO_DATA
  TestCurveRemoveKnot.cxx,NO_DATA
  TestCurveBezierExtraction.cxx,NO_DATA
  TestCurveIncreaseDegree.cxx,NO_DATA
  TestCurveDecreaseDegree.cxx,NO_DATA
  TestSurfaceInsertKnot.cxx,NO_DATA
  TestSurfaceKnotRefinement.cxx,NO_DATA
  TestSurfaceRemoveKnot.cxx,NO_DATA,NO_VALID,NO_OUTPUT
  TestSurfaceBezierExtraction.cxx,NO_DATA
  TestSurfaceIncreaseDegree.cxx,NO_DATA
  TestCylinderVolume.cxx,NO_DATA)

vtk_test_cxx_executable(${vtk-module}CxxTests tests
  RENDERING_FACTORY)
