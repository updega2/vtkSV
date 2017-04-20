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

#include <vtkCamera.h>
#include "vtkSVControlGrid.h"
#include "vtkDataObject.h"
#include "vtkPoints.h"
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVNURBSSurface.h"
#include "vtkSVNURBSUtils.h"
#include "vtkStructuredGridGeometryFilter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int TestSurfaceInsertKnot(int argc, char *argv[])
{
  // Set the curve knot insertion details
  int p=2;     //degree
  int q=2;     //degree
  int np=3;     //control points
  int mp=3;     //control points
  int nuk=p+np+1; //nuk
  int nvk=q+mp+1; //nvk
  int r=1;     // number of insertions we want to do
  double u=0.5;   // insertion value

  // Control points (simple grid)
  vtkNew(vtkStructuredGrid, controlPointGrid);
  vtkNew(vtkPoints, controlPoints);
  controlPointGrid->SetDimensions(np, mp, 1);
  controlPoints->SetNumberOfPoints(np*mp);
  controlPoints->SetPoint(0, 0.0, 0.0, 0.0);
  controlPoints->SetPoint(1, 1.0, 0.0, 1.0);
  controlPoints->SetPoint(2, 3.0, 1.0, -1.0);
  controlPoints->SetPoint(3, 0.0, 1.0, 1.0);
  controlPoints->SetPoint(4, 1.0, -2.0, -1.0);
  controlPoints->SetPoint(5, 2.0, 0.0, 0.0);
  controlPoints->SetPoint(6, 0.0, 5.0, 0.0);
  controlPoints->SetPoint(7, 2.0, -1.0, -2.0);
  controlPoints->SetPoint(8, 4.0, 1.0, 1.0);
  controlPointGrid->SetPoints(controlPoints);

  // uKnot vector equal space
  vtkNew(vtkDoubleArray, uKnots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, nuk, p, uKnots);

  // vKnot vector equal space
  vtkNew(vtkDoubleArray, vKnots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, nvk, q, vKnots);

  // Set up the surface
  vtkNew(vtkSVNURBSSurface, surface);
  surface->SetControlPoints(controlPointGrid);
  surface->SetUKnotVector(uKnots);
  surface->SetVKnotVector(vKnots);

  surface->GeneratePolyDataRepresentation(0.1, 0.1);

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/TEST.vtp";
  vtkSVIOUtils::WriteVTPFile(filename, surface->GetSurfaceRepresentation());


  // Double check the actual values

  return EXIT_SUCCESS;
}
