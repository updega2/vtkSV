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
#include "vtkSVMathUtils.h"
#include "vtkStructuredGridGeometryFilter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int TestCurveDecreaseDegree(int argc, char *argv[])
{
  // Set the curve decrease details
  int p=5;     //degree
  int np=30;     //control points
  int m=p+np+1; //m

  // Set the knot vector
  vtkNew(vtkDoubleArray, knots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, m, p, knots);

  // X and Y data values
  double xmin = 0.0; double xmax=5.0;
  vtkNew(vtkDoubleArray, xvals);
  vtkSVNURBSUtils::LinSpace(xmin, xmax, np, xvals);

  // Set the control points
  vtkNew(vtkPoints, cpoints);
  cpoints->SetNumberOfPoints(np);
  for (int i=0; i<np; i++)
  {
    double x = xvals->GetTuple1(i);
    double y = pow(x, 1./2);
    cpoints->SetPoint(i, x, y, 0.0);
  }

  // Set up the curve
  vtkNew(vtkSVNURBSCurve, curve);
  curve->SetKnotVector(knots);
  curve->SetControlPoints(cpoints);
  curve->SetDegree(p);

  // Increase the degree first to make duplicate knot points. Both degree
  // reduction and knot removal achieve higher accuracy when there are
  // duplicate knot points
  // Comment out for baseline
  curve->IncreaseDegree(4);
  curve->DecreaseDegree(3.0); // Now reduce degree

  // Generate representation
  curve->GeneratePolyDataRepresentation(0.01);

  // Set up mapper
  vtkNew(vtkPolyDataMapper, mapper);
  mapper->SetInputData(curve->GetCurveRepresentation());

  // Set up actor
  vtkNew(vtkActor, actor);
  actor->SetMapper(mapper);

  // Set up renderer and window
  vtkNew(vtkRenderer, renderer);
  vtkNew(vtkRenderWindow, renWin);
  renWin->AddRenderer( renderer );
  renderer->AddActor(actor);
  renderer->SetBackground(.1, .2, .3);

  // Set up interactor
  vtkNew(vtkRenderWindowInteractor, renWinInteractor);
  renWinInteractor->SetRenderWindow( renWin );

  // Render
  renWin->Render();
  renWinInteractor->Start();

  return EXIT_SUCCESS;
}
