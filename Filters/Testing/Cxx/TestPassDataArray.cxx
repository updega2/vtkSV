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
 *  \file TestPassDataArray.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */
#include "vtkSVPassDataArray.h"

#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkTestUtilities.h"

int TestPassDataArray(int argc, char *argv[])
{
  // Read the surface
  vtkNew(vtkPolyData, sourcePd);
  char *source_filename = vtkTestUtilities::ExpandDataFileName(
    argc, argv, "0141_1001_Renal_Branch_Centerlines.vtp");
  vtkSVIOUtils::ReadVTPFile(source_filename, sourcePd);
  vtkNew(vtkPolyData, targetPd);
  char *target_filename = vtkTestUtilities::ExpandDataFileName(
    argc, argv, "0141_1001_Renal_Branch_Surface.vtp");
  vtkSVIOUtils::ReadVTPFile(target_filename, targetPd);

  // Set up vars for separator
  std::string passArrayName = "GroupIds";

  // Set up path passer
  vtkNew(vtkSVPassDataArray, passer);
  passer->SetInputData(0, sourcePd);
  passer->SetInputData(1, targetPd);
  passer->SetPassArrayName(passArrayName.c_str());
  passer->SetPassDataIsCellData(1);
  passer->SetPassDataToCellData(1);
  passer->Update();

  // Get output
  vtkNew(vtkPolyData, output);
  output = passer->GetOutput();
  output->GetCellData()->SetActiveScalars(passArrayName.c_str());

  // Get scalar range
  double range[2];
  output->GetCellData()->GetArray(passArrayName.c_str())->GetRange(range);

  // Set up mapper
  vtkNew(vtkPolyDataMapper, mapper);
  mapper->SetInputData(output);
  mapper->SetScalarRange(range);
  mapper->SetScalarModeToUseCellData();
  mapper->ScalarVisibilityOn();

  // Set up actor
  vtkNew(vtkActor, actor);
  actor->SetMapper(mapper);

  // Set up renderer and window
  vtkNew(vtkRenderer, renderer);
  vtkNew(vtkRenderWindow, renWin);
  renWin->AddRenderer( renderer );
  renderer->AddActor(actor);

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
