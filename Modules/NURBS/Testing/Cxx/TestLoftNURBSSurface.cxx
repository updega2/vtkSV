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
#include "vtkDataObject.h"
#include "vtkPoints.h"
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridGeometryFilter.h"

#include "vtkSVControlGrid.h"
#include "vtkSVGlobals.h"
#include "vtkSVLoftNURBSSurface.h"
#include "vtkSVNURBSSurface.h"
#include "vtkSVNURBSUtils.h"
#include "vtkSVIOUtils.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int TestLoftNURBSSurface(int argc, char *argv[])
{

  // Set number of initial data points
  int nU = 57;
  int nV = 63;
  int dim2D[3];
  dim2D[0] = nU;
  dim2D[1] = nV;
  dim2D[2] = 1;

  vtkNew(vtkPoints, newPoints);
  newPoints->SetNumberOfPoints(nU*nV);

  vtkNew(vtkStructuredGrid, inputGrid);
  inputGrid->SetPoints(newPoints);
  inputGrid->SetDimensions(dim2D);

  // Create data
  vtkNew(vtkDoubleArray, xdata);
  xdata->SetNumberOfTuples(nU);
  vtkNew(vtkDoubleArray, ydata);
  ydata->SetNumberOfTuples(nV);
  vtkSVNURBSUtils::LinSpace(0, 2*SV_PI, nU, xdata);
  vtkSVNURBSUtils::LinSpace(0, 2*SV_PI, nV, ydata);

  // Set end derivatives in u
  vtkNew(vtkDoubleArray, uders);
  uders->SetNumberOfComponents(3);
  uders->SetNumberOfTuples(nU);
  double uder[3]; uder[0] = 2.0; uder[1] = 0.0; uder[2] = 0.0;
  for (int i=0; i<nU; i++)
    uders->SetTuple(i, uder);

  // Set end derivatives in v
  vtkNew(vtkDoubleArray, vders);
  vders->SetNumberOfComponents(3);
  vders->SetNumberOfTuples(nV);
  double vder[3]; vder[0] = 0.0; vder[1] = 2.0; vder[2] = 0.0;
  for (int i=0; i<nV; i++)
    vders->SetTuple(i, vder);

  // Set up all the polydatas
  int pos[3];
  for (int i=0; i<nU; i++)
  {
    for (int j=0; j<nV; j++)
    {
      pos[0] = i; pos[1] = j; pos[2] = 0;
      int ptId = vtkStructuredData::ComputePointId(dim2D, pos);

      double xval = xdata->GetTuple1(i);
      double yval = ydata->GetTuple1(j);;
      double zval = sin(xval) * sin(yval);

      inputGrid->GetPoints()->SetPoint(ptId, xval, yval, zval);
    }
  }

  // Filter
  vtkNew(vtkSVLoftNURBSSurface, lofter);
  lofter->SetInputData(inputGrid);
  lofter->SetUDegree(3);
  lofter->SetVDegree(3);
  lofter->SetPolyDataUSpacing(0.05);
  lofter->SetPolyDataVSpacing(0.05);
  lofter->SetUKnotSpanType("average");
  lofter->SetUParametricSpanType("chord");
  lofter->SetStartUDerivatives(uders);
  lofter->SetEndUDerivatives(uders);
  lofter->SetVKnotSpanType("derivative");
  lofter->SetVParametricSpanType("centripetal");
  lofter->SetStartVDerivatives(vders);
  lofter->SetEndVDerivatives(vders);
  lofter->Update();


  // Get outputs
  vtkNew(vtkPolyData, surfacePd);
  surfacePd->DeepCopy(lofter->GetOutput());
  vtkNew(vtkStructuredGridGeometryFilter, converter);
  converter->SetInputData(lofter->GetSurface()->GetControlPointGrid());
  converter->Update();
  vtkNew(vtkPolyData, controlGridPd);
  controlGridPd->DeepCopy(converter->GetOutput());

  // Set up mapper
  vtkNew(vtkPolyDataMapper, mapper);
  mapper->SetInputData(surfacePd);

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
  //renderer->AddActor(gridActor);
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
