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
#include "vtkAppendPolyData.h"
#include "vtkAppendFilter.h"
#include "vtkSVControlGrid.h"
#include "vtkDataObject.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkPoints.h"
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVNURBSVolume.h"
#include "vtkSVNURBSUtils.h"
#include "vtkStructuredGridGeometryFilter.h"

#include <string>
#include <sstream>
#include <iostream>

#include <unistd.h>

int TestCylinderVolume(int argc, char *argv[])
{
  // Set the curve knot insertion details
  int p=2;     //degree
  int q=2;     //degree
  int r=2;     //degree
  int np=9;     //control points
  int mp=4;     //control points
  int lp=3;     //control points
  int nuk=p+np+1; //nuk
  int nvk=q+mp+1; //nvk
  int nwk=r+lp+1; //nwk

  // Control points for a cylinder,
  // quadratic for all dims
  double r_inner = 1.0;
  double r_outer = 2.0;
  double l = 5;
  double ldiv = l/(lp-1);
  double rdiv = (r_outer-r_inner)/(mp-1);
  vtkNew(vtkSVControlGrid, controlPointGrid);
  controlPointGrid->SetDimensions(np, mp, lp);
  controlPointGrid->SetNumberOfControlPoints(np*mp*lp);
  for (int j=0; j<mp; j++)
  {
    for (int k=0; k<lp; k++)
    {
      double l_loc = k*ldiv;
      double r_loc = r_inner + j*rdiv;

      controlPointGrid->SetControlPoint(0, j, k, r_loc, 0.0, l_loc, 1.0);
      controlPointGrid->SetControlPoint(1, j, k, r_loc, r_loc, l_loc, sqrt(2)/2);
      controlPointGrid->SetControlPoint(2, j, k, 0.0, r_loc, l_loc, 1.0);
      controlPointGrid->SetControlPoint(3, j, k, -r_loc, r_loc, l_loc, sqrt(2)/2);
      controlPointGrid->SetControlPoint(4, j, k, -r_loc, 0.0, l_loc, 1.0);
      controlPointGrid->SetControlPoint(5, j, k, -r_loc, -r_loc, l_loc, sqrt(2)/2);
      controlPointGrid->SetControlPoint(6, j, k, 0.0, -r_loc, l_loc, 1.0);
      controlPointGrid->SetControlPoint(7, j, k, r_loc, -r_loc, l_loc, sqrt(2)/2);
      controlPointGrid->SetControlPoint(8, j, k, r_loc, 0.0, l_loc, 1.0);
    }
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

  // vKnot vector equal space
  vtkNew(vtkDoubleArray, wKnots);
  vtkSVNURBSUtils::LinSpaceClamp(0, 1, nwk, r, wKnots);

  // Set up the volume
  vtkNew(vtkSVNURBSVolume, volume);
  volume->SetControlPointGrid(controlPointGrid);
  volume->SetUKnotVector(uKnots);
  volume->SetVKnotVector(vKnots);
  volume->SetWKnotVector(wKnots);
  volume->SetUDegree(p);
  volume->SetVDegree(q);
  volume->SetWDegree(r);

  // Generate volume
  volume->GenerateVolumeRepresentation(0.01, 0.1, 0.1);

  // Get surface
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(volume->GetVolumeRepresentation());
  surfacer->Update();

  // Set up mapper
  vtkNew(vtkPolyDataMapper, mapper);
  mapper->SetInputData(surfacer->GetOutput());

  // Set up actor
  vtkNew(vtkActor, actor);
  actor->SetMapper(mapper);
  actor->GetProperty()->SetEdgeColor(0.0, 0.0, 1.0); //(R,G,B)
  actor->GetProperty()->EdgeVisibilityOn();

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
