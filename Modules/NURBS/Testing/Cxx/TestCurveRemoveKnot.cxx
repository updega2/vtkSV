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

#include "vtkCamera.h"
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
#include "vtkSVLoftNURBSCurve.h"
#include "vtkSVNURBSUtils.h"
#include "vtkStructuredGridGeometryFilter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int main(int argc, char *argv[])
{
  // Set the curve knot insertion details
  int p=3;     //degree
  int np=10;     //control points
  int m=p+np+1; //m
  double u=1./7;   // remove value
  int r=1;     // number of removals

  // Set the knot vector
  vtkNew(vtkDoubleArray, knots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, m, p, knots);

  // Set the control points
  vtkNew(vtkPoints, cpoints);
  cpoints->SetNumberOfPoints(np);
  cpoints->SetPoint(0, 0.0, 0.0, 0.0);
  cpoints->SetPoint(1, 0.0, 1.0, 0.0);
  cpoints->SetPoint(2, 1.0, 1.5, 0.0);
  cpoints->SetPoint(3, 2.0, 2.0, 0.0);
  cpoints->SetPoint(4, 3.0, 2.5, 0.0);
  cpoints->SetPoint(5, 4.0, 3.0, 0.0);
  cpoints->SetPoint(6, 5.0, 3.5, 0.0);
  cpoints->SetPoint(7, 6.0, 3.5, 0.0);
  cpoints->SetPoint(8, 7.0, 3.0, 0.0);
  cpoints->SetPoint(9, 8.0, 2.5, 0.0);

  // Set up the curve
  vtkNew(vtkSVNURBSCurve, curve);
  curve->SetKnotVector(knots);
  curve->SetControlPoints(cpoints);
  curve->SetDegree(p);

  curve->RemoveKnot(u, r, 0.01); // Remove the first knot

  // Check the result
  if (curve->GetNumberOfControlPoints() != np-r)
  {
    fprintf(stderr, "Curve was not updated to correct number of control points\n");
    return EXIT_FAILURE;
  }

  if (curve->GetNumberOfKnotPoints() != m-r)
  {
    fprintf(stderr, "Curve was not updated to correct number of knot points\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
