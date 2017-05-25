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
 *  \file TestPolyDataToNURBSFilter.cxx
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */
#include "vtkSVPolyDataToNURBSFilter.h"

#include "vtkCellData.h"
#include "vtkPlaneSource.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVGeneralizedPolycube.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVHausdorffDistance.h"
#include "vtkSVIOUtils.h"
#include "vtkTestUtilities.h"

int TestPolyDataToNURBSFilter(int argc, char *argv[])
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
  vtkNew(vtkSVPolyDataToNURBSFilter, Converter);

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
  std::string booleanPathArrayName      = "IsPath";


  // OPERATION
  Converter->SetInputData(surfacePd);
  Converter->SetCenterlinesPd(centerlinesPd);
  Converter->SetAddTextureCoordinates(1);
  Converter->SetBoundaryPointsArrayName(boundaryPointsArrayName.c_str());
  Converter->SetGroupIdsArrayName(groupIdsArrayName.c_str());
  Converter->SetSegmentIdsArrayName(segmentIdsArrayName.c_str());
  Converter->SetSliceIdsArrayName(sliceIdsArrayName.c_str());
  Converter->SetSphereRadiusArrayName(radiusArrayName.c_str());
  Converter->SetInternalIdsArrayName(internalIdsArrayName.c_str());
  Converter->SetDijkstraArrayName(dijkstraDistanceArrayName.c_str());
  Converter->SetBooleanPathArrayName(booleanPathArrayName.c_str());
  Converter->Update();

  // TODO: Add checks
  return EXIT_SUCCESS;
}
