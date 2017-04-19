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
#include "vtkSVLoftNURBSCurve.h"
#include "vtkSVNURBSUtils.h"
#include "vtkStructuredGridGeometryFilter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int TestLoftNURBSCurve(int argc, char *argv[])
{
  // Set number of initial data points
  int nc = 11;

  // Set up input data
  vtkNew(vtkDoubleArray, x_data);
  x_data->SetNumberOfTuples(nc);
  vtkSVNURBSUtils::LinSpace(0, 2*M_PI, nc, x_data);

  vtkNew(vtkPoints, inputPoints);
  inputPoints->Reset();
  for (int i=0; i<nc; i++)
  {
    double xval = x_data->GetTuple1(i);
    double yval = sin(xval - M_PI/2.0);
    double zval = 0.0;

    inputPoints->InsertNextPoint(xval, yval, zval);
  }

  // Create data
  vtkNew(vtkPolyData, inputPoly);
  inputPoly->SetPoints(inputPoints);

  // Lofter
  vtkNew(vtkSVLoftNURBSCurve, lofter);
  lofter->SetInputData(inputPoly);
  lofter->SetDegree(3);
  lofter->SetPolyDataSpacing(0.05);
  lofter->SetKnotSpanType("derivative");
  lofter->SetParametricSpanType("chord");
  lofter->Update();

  // Get outputs
  vtkNew(vtkPolyData, curvePd);
  curvePd->DeepCopy(lofter->GetOutput());
  vtkNew(vtkStructuredGridGeometryFilter, converter);
  converter->SetInputData(lofter->GetCurve()->GetControlPointGrid());
  converter->Update();
  vtkNew(vtkPolyData, controlGridPd);
  controlGridPd->DeepCopy(converter->GetOutput());

  // Set up mapper
  vtkNew(vtkPolyDataMapper, mapper);
  mapper->SetInputData(curvePd);

  vtkNew(vtkPolyDataMapper, gridMapper);
  gridMapper->SetInputData(controlGridPd);

  // Set up actor
  vtkNew(vtkActor, actor);
  actor->SetMapper(mapper);

  // Set up grid actor
  vtkNew(vtkActor, gridActor);
  gridActor->SetMapper(gridMapper);

  // Set up renderer and window
  vtkNew(vtkRenderer, renderer);
  vtkNew(vtkRenderWindow, renWin);
  renWin->AddRenderer( renderer );
  renderer->AddActor(actor);
  renderer->AddActor(gridActor);
  renderer->SetBackground(.1, .2, .3);

  // Set up interactor
  vtkNew(vtkRenderWindowInteractor, renWinInteractor);
  renWinInteractor->SetRenderWindow( renWin );

  // TODO: Camera and bettering coloring for test
  // Camera change
  //vtkCamera *camera = renderer->GetActiveCamera();
  //renderer->ResetCamera();
  //camera->Elevation(60);

  // Render
  renWin->Render();
  renWinInteractor->Start();

  return EXIT_SUCCESS;
}
