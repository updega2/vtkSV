/*=========================================================================
 *
 * Copyright (c) 2014 The Regents of the University of California.
 * All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *=========================================================================*/

#include <vtkSVLoopIntersectionPolyDataFilter.h>

#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

int TestLoopIntersectionPolyDataFilter2(int argc, char *argv[])
{
  // Set up two polydata representing two triangles that share a vertex
  int coplanar;
  double isectpt1[3], isectpt2[3];
  double thisCellTri[9] = {-30.125,
                           29.3125,
                           -27.1875,
                           -29.9375,
                           29.375,
                           -27.3125,
                           -30.0625,
                           28.5,
                           -27.25};
  double otherCellTri[9] = {-29.9375,
                            29.3125,
                            -27.3125,
                            -29.875,
                            29.8125,
                            -27.5,
                            -29.75,
                            27.6875,
                            -27.4375};

  double surfaceid[2];
  double tolerance = 1e-6;
  int intersects = vtkSVLoopIntersectionPolyDataFilter
    ::TriangleTriangleIntersection(&thisCellTri[0],
                                   &thisCellTri[3],
                                   &thisCellTri[6],
                                   &otherCellTri[0],
                                   &otherCellTri[3],
                                   &otherCellTri[6],
                                   coplanar,
                                   isectpt1, isectpt2,
                                   surfaceid, tolerance);

    std::cerr << "First: "
              << thisCellTri[0] << ", "
              << thisCellTri[3] << ", "
              << thisCellTri[6] << std::endl;
    std::cerr << "Second: "
              << otherCellTri[0] << ", "
              << otherCellTri[3] << ", "
              << otherCellTri[6] << std::endl;
  if ( intersects )
    {
    std::cerr << "Triangles with shared vertex should not be reported to intersect" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
