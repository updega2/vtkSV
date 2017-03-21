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

/**
 *  \file TestPullApartPolyData.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */
#include "vtkSVPullApartPolyData.h"

#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkConnectivityFilter.h>
#include <vtkFeatureEdges.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include "vtkSmartPointer.h"
#include "vtkSVFindGeodesicPath.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkTestUtilities.h"

int TestPullApartPolyData(int argc, char *argv[])
{
  // Read the surface
  vtkNew(vtkPolyData, surfacePd);
  char *surface_filename = vtkTestUtilities::ExpandDataFileName(
    argc, argv, "0141_1001_Renal_Branch_Surface.vtp");
  vtkSVIOUtils::ReadVTPFile(surface_filename, surfacePd);

  // Set up vars
  int startPtId = 953;
  int endPtId   = 13012;
  std::string dijkstraArrayName = "DijkstraDistance";
  std::string internalIdsArrayName = "TmpInternalIds";
  std::string pathBooleanArrayName = "IsPath";

  // Set up path finder
  vtkNew(vtkSVFindGeodesicPath, finder);
  finder->SetInputData(surfacePd);
  finder->SetStartPtId(startPtId); // Specific to given data set
  finder->SetEndPtId(endPtId); // Specific to given data set
  finder->SetDijkstraArrayName(dijkstraArrayName.c_str());
  finder->SetInternalIdsArrayName(internalIdsArrayName.c_str());
  finder->SetPathBooleanArrayName(pathBooleanArrayName.c_str());
  finder->SetRepelCloseBoundaryPoints(1);
  finder->SetAddPathBooleanArray(1);
  finder->Update();

  // Set up pull aparter
  vtkNew(vtkSVPullApartPolyData, ripper);
  ripper->SetInputData(finder->GetOutput());
  ripper->SetCutPointsArrayName(pathBooleanArrayName.c_str());
  ripper->SetStartPtId(startPtId);
  ripper->Update();

  // Get output
  vtkNew(vtkPolyData, output);
  output = ripper->GetOutput();
  output->GetPointData()->SetActiveScalars(pathBooleanArrayName.c_str());

  // Get boundary edges
  vtkNew(vtkFeatureEdges, boundaryFinder);
  boundaryFinder->SetInputData(output);
  boundaryFinder->BoundaryEdgesOn();
  boundaryFinder->FeatureEdgesOff();
  boundaryFinder->ManifoldEdgesOff();
  boundaryFinder->NonManifoldEdgesOff();
  boundaryFinder->Update();

  // Extractor
  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(boundaryFinder->GetOutput());
  connector->SetExtractionModeToAllRegions();
  connector->ColorRegionsOn();
  connector->Update();
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  surfacer->GetOutput()->GetCellData()->SetActiveScalars("RegionId");

  int expectedRegions = 2;
  if (connector->GetNumberOfExtractedRegions() != expectedRegions)
  {
    fprintf(stderr, "Incorrect number of boundaries on final surface: %d\n", connector->GetNumberOfExtractedRegions());
    fprintf(stderr, "Expected: %d\n", expectedRegions);
    return EXIT_FAILURE;
  }

  // Get scalar range
  double range[2];
  output->GetPointData()->GetArray(pathBooleanArrayName.c_str())->GetRange(range);

  // Set up mapper
  vtkNew(vtkPolyDataMapper, mapper);
  mapper->SetInputData(output);
  mapper->SetScalarRange(range);
  mapper->SetScalarModeToUsePointData();
  mapper->ScalarVisibilityOn();

  // Set up actor
  vtkNew(vtkActor, actor);
  actor->SetMapper(mapper);

  // Set up ripped lines mapper
  vtkNew(vtkPolyDataMapper, linesMapper);
  linesMapper->SetInputData(surfacer->GetOutput());
  linesMapper->SetScalarRange(0, expectedRegions);
  linesMapper->SetScalarModeToUseCellData();
  linesMapper->ScalarVisibilityOn();

  // Set up ripped lines actor
  vtkNew(vtkActor, linesActor);
  linesActor->SetMapper(linesMapper);

  // Set up renderer window
  vtkNew(vtkRenderer, renderer);
  vtkNew(vtkRenderWindow, renWin);
  renWin->AddRenderer( renderer );
  renderer->AddActor(actor);
  renderer->AddActor(linesActor);
  renderer->SetBackground(.1, .2, .3);

  // Set up interactor
  vtkNew(vtkRenderWindowInteractor, renWinInteractor);
  renWinInteractor->SetRenderWindow( renWin );

  // Camera change
  vtkCamera *camera = renderer->GetActiveCamera();
  renderer->ResetCamera();
  camera->Azimuth(180);
  camera->Elevation(90);

  // Render
  renWin->Render();
  renWinInteractor->Start();

  return EXIT_SUCCESS;
}
