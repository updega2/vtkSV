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
 *  \file TestPolyDataSliceAndDiceFilter.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */
#include "vtkSVPolyDataSliceAndDiceFilter.h"

#include "vtkCellData.h"
#include "vtkPlaneSource.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVGeneralizedPolycube.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVHausdorffDistance.h"
#include "vtkSVIOUtils.h"
#include "vtkTestUtilities.h"

int TestPolyDataSliceAndDiceFilter(int argc, char *argv[])
{
  // Read the surface
  vtkNew(vtkPolyData, surfacePd);
  char *surface_filename = vtkTestUtilities::ExpandDataFileName(
    argc, argv, "0110_0001_Iliac_Branch_Surface.vtp");
  vtkSVIOUtils::ReadVTPFile(surface_filename, surfacePd);

  // Read the centerlines
  vtkNew(vtkPolyData, centerlinesPd);
  char *centerlines_filename = vtkTestUtilities::ExpandDataFileName(
    argc, argv, "0110_0001_Iliac_Branch_Centerlines.vtp");
  vtkSVIOUtils::ReadVTPFile(centerlines_filename, centerlinesPd);

  // Filter
  vtkNew(vtkSVPolyDataSliceAndDiceFilter, Slicer);

  // Set up values to use
  double sliceLength    = 1.0;
  int constructPolycube = 1;
  std::string boundaryPointsArrayName   = "BoundaryPoints";
  std::string groupIdsArrayName         = "GroupIds";
  std::string segmentIdsArrayName       = "SegmentIds";
  std::string sliceIdsArrayName         = "SliceIds";
  std::string radiusArrayName           = "MaximumInscribedSphereRadius";
  std::string internalIdsArrayName      = "TmpInternalIds";
  std::string dijkstraDistanceArrayName = "DijkstraDistance";


  // OPERATION
  Slicer->SetInputData(surfacePd);
  Slicer->SetCenterlinesPd(centerlinesPd);
  Slicer->SetSliceLength(sliceLength);
  Slicer->SetConstructPolycube(constructPolycube);
  Slicer->SetBoundaryPointsArrayName(boundaryPointsArrayName.c_str());
  Slicer->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  Slicer->SetSegmentIdsArrayName(segmentIdsArrayName.c_str());
  Slicer->SetSliceIdsArrayName(sliceIdsArrayName.c_str());
  Slicer->SetSphereRadiusArrayName(radiusArrayName.c_str());
  Slicer->SetInternalIdsArrayName(internalIdsArrayName.c_str());
  Slicer->SetDijkstraArrayName(dijkstraDistanceArrayName.c_str());
  Slicer->Update();

  int numExpectedSegments = 7;
  int numExpectedGroups   = 5;

  // Get output
  vtkNew(vtkPolyData, output);
  output = Slicer->GetOutput();
  if (vtkSVGeneralUtils::CheckArrayExists(output, 1, groupIdsArrayName) != SV_OK)
  {
    std::cout << "Group ids have not been attached to surface" << endl;
    return EXIT_FAILURE;
  }
  double minmax[2];
  output->GetCellData()->GetArray(groupIdsArrayName.c_str())->GetRange(minmax);
  if ((minmax[1]+1) != numExpectedSegments)
  {
    std::cout << "Incorrect number of output segments: " << minmax[1]+1 << endl;
    std::cout << "Should be: " << numExpectedSegments << endl;
    return EXIT_FAILURE;
  }

  // Get polycubes
  vtkNew(vtkSVGeneralizedPolycube, polycube);
  polycube = Slicer->GetPolycube();

  // Check to make sure has correct number of cubes
  if (polycube->GetNumberOfGrids() != numExpectedSegments)
  {
    std::cout << "Polycube does not have the correct number of cubes: " << polycube->GetNumberOfGrids() <<endl;
    std::cout << "Should be: " << numExpectedSegments << endl;
    return EXIT_FAILURE;
  }

  // Get the number of groups (number of branch cubes)
  int numGroups = 0;
  for (int i=0; i<polycube->GetNumberOfGrids(); i++)
  {
    if (polycube->GetCellData()->GetArray("CubeType")->GetTuple1(i) ==
      vtkSVGeneralizedPolycube::CUBE_BRANCH)
    {
      numGroups++;
    }
  }

  // Check to make sure correct number of groups found
  if (numGroups != numExpectedGroups)
  {
    std::cout << "Incorrect number of branching regions: " << numGroups <<endl;
    std::cout << "Should be: " << numExpectedGroups << endl;
    return EXIT_FAILURE;
  }

  // Get surgerylines
  vtkNew(vtkPolyData, surgeryLinesPd);
  surgeryLinesPd = Slicer->GetSurgeryLinesPd();
  if (surgeryLinesPd->GetNumberOfLines() != numExpectedGroups)
  {
    std::cout << "Incorrect number of surgery lines: " << surgeryLinesPd->GetNumberOfLines() <<endl;
    std::cout << "Should be: " << numExpectedGroups << endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
