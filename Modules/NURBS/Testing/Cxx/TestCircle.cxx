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

int TestCircle(int argc, char *argv[])
{
  // Set the curve knot insertion details
  int p=2;     //degree
  int np=9;     //control points
  int m=p+np+1; //m
  int r=2;     // number of insertions we want to do
  double u=0.6;   // insertion value

  // Set the knot vector
  vtkNew(vtkDoubleArray, knots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, m, p, knots);
  knots->SetTuple1(3, 1./4);
  knots->SetTuple1(4, 1./4);
  knots->SetTuple1(5, 1./2);
  knots->SetTuple1(6, 1./2);
  knots->SetTuple1(7, 3./4);
  knots->SetTuple1(8, 3./4);

  // Set the control points
  vtkNew(vtkSVControlGrid, controlPointGrid);
  controlPointGrid->SetDimensions(np, 1, 1);
  controlPointGrid->SetNumberOfControlPoints(np);
  controlPointGrid->SetControlPoint(0, 0, 0, 1.0, 0.0, 0.0, 1.0);
  controlPointGrid->SetControlPoint(1, 0, 0, 1.0, 1.0, 0.0, sqrt(2)/2);
  controlPointGrid->SetControlPoint(2, 0, 0, 0.0, 1.0, 0.0, 1.0);
  controlPointGrid->SetControlPoint(3, 0, 0, -1.0, 1.0, 0.0, sqrt(2)/2);
  controlPointGrid->SetControlPoint(4, 0, 0, -1.0, 0.0, 0.0, 1.0);
  controlPointGrid->SetControlPoint(5, 0, 0, -1.0, -1.0, 0.0, sqrt(2)/2);
  controlPointGrid->SetControlPoint(6, 0, 0, 0.0, -1.0, 0.0, 1.0);
  controlPointGrid->SetControlPoint(7, 0, 0, 1.0, -1.0, 0.0, sqrt(2)/2);
  controlPointGrid->SetControlPoint(8, 0, 0, 1.0, 0.0, 0.0, 1.0);

  // Set up the curve
  vtkNew(vtkSVNURBSCurve, curve);
  curve->SetKnotVector(knots);
  curve->SetControlPointGrid(controlPointGrid);
  curve->SetDegree(p);

  curve->GeneratePolyDataRepresentation(0.01);

  // Insert same knot r times
  curve->InsertKnot(u, r);

  // Check the result
  if (curve->GetNumberOfControlPoints() != np+r)
  {
    std::cerr << "Curve was not updated to correct number of control points" << endl;
    return EXIT_FAILURE;
  }

  if (curve->GetNumberOfKnotPoints() != m+r)
  {
    std::cerr << "Curve was not updated to correct number of knot points" << endl;
    return EXIT_FAILURE;
  }

  int q=3;
  vtkNew(vtkDoubleArray, newKnots);
  newKnots->SetNumberOfTuples(q);
  newKnots->SetTuple1(0, 0.3);
  newKnots->SetTuple1(1, 0.3);
  newKnots->SetTuple1(2, 0.6);

  // Insert different knots into span using refinement
  curve->InsertKnots(newKnots);

  // Check the result
  if (curve->GetNumberOfControlPoints() != np+r+q)
  {
    std::cerr << "Curve was not updated to correct number of control points" << endl;
    return EXIT_FAILURE;
  }

  if (curve->GetNumberOfKnotPoints() != m+r+q)
  {
    std::cerr << "Curve was not updated to correct number of knot points" << endl;
    return EXIT_FAILURE;
  }

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
