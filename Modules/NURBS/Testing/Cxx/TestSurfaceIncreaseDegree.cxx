/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
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
 */

#include "vtkCamera.h"
#include "vtkSVControlGrid.h"
#include "vtkDataObject.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVNURBSSurface.h"
#include "vtkSVNURBSUtils.h"
#include "vtkStructuredGridGeometryFilter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

// Bug for setting mp=3, increasing degree
int TestSurfaceIncreaseDegree(int argc, char *argv[])
{
  // Set the curve knot insertion details
  int p=2;     //degree
  int q=1;     //degree
  int np=9;     //control points
  int mp=2;     //control points
  int nuk=p+np+1; //nuk
  int nvk=q+mp+1; //nvk

  // Control points for a cylinder, two circles really
  double ldiv = 10.0/(mp-1);
  vtkNew(vtkSVControlGrid, controlPointGrid);
  controlPointGrid->SetDimensions(np, mp, 1);
  controlPointGrid->SetNumberOfControlPoints(np*mp);
  for (int j=0; j<mp; j++)
  {
    controlPointGrid->SetControlPoint(0, j, 0, 1.0, 0.0, (j*ldiv), 1.0);
    controlPointGrid->SetControlPoint(1, j, 0, 1.0, 1.0, (j*ldiv), sqrt(2)/2);
    controlPointGrid->SetControlPoint(2, j, 0, 0.0, 1.0, (j*ldiv), 1.0);
    controlPointGrid->SetControlPoint(3, j, 0, -1.0, 1.0, (j*ldiv), sqrt(2)/2);
    controlPointGrid->SetControlPoint(4, j, 0, -1.0, 0.0, (j*ldiv), 1.0);
    controlPointGrid->SetControlPoint(5, j, 0, -1.0, -1.0, (j*ldiv), sqrt(2)/2);
    controlPointGrid->SetControlPoint(6, j, 0, 0.0, -1.0, (j*ldiv), 1.0);
    controlPointGrid->SetControlPoint(7, j, 0, 1.0, -1.0, (j*ldiv), sqrt(2)/2);
    controlPointGrid->SetControlPoint(8, j, 0, 1.0, 0.0, (j*ldiv), 1.0);
  }

  // uKnot vector equal space
  vtkNew(vtkDoubleArray, uKnots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, nuk, p, uKnots);
  // Replace interior knot values with circle vals
  uKnots->SetTuple1(3, 1./4);
  uKnots->SetTuple1(4, 1./4);
  uKnots->SetTuple1(5, 1./2);
  uKnots->SetTuple1(6, 1./2);
  uKnots->SetTuple1(7, 3./4);
  uKnots->SetTuple1(8, 3./4);

  // vKnot vector equal space
  vtkNew(vtkDoubleArray, vKnots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, nvk, q, vKnots);

  // Set up the surface
  vtkNew(vtkSVNURBSSurface, surface);
  surface->SetControlPointGrid(controlPointGrid);
  surface->SetUKnotVector(uKnots);
  surface->SetVKnotVector(vKnots);
  surface->SetUDegree(p);
  surface->SetVDegree(q);

  // Now lets increase the degree
  double du = 2;
  double dv = 3;
  surface->IncreaseDegree(du, 0);
  surface->IncreaseDegree(dv, 1);

  // Check the result
  if (surface->GetUDegree() != p+du)
  {
    std::cerr << "Curve was not updated to correct degree" << endl;
    return EXIT_FAILURE;
  }
  if (surface->GetVDegree() != q+dv)
  {
    std::cerr << "Curve was not updated to correct degree" << endl;
    return EXIT_FAILURE;
  }

  surface->GeneratePolyDataRepresentation(0.01, 0.01);

  // Set up mapper
  vtkNew(vtkPolyDataMapper, mapper);
  mapper->SetInputData(surface->GetSurfaceRepresentation());

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

  vtkCamera *camera = renderer->GetActiveCamera();
  renderer->ResetCamera();
  camera->Azimuth(45);
  camera->Elevation(45);

  // Render
  renWin->Render();
  renWinInteractor->Start();

  return EXIT_SUCCESS;
}
