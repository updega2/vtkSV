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

/** @file vtkMakeParametricCube.cxx
 *  @brief This implements the vtkMakeParametricCube filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkMakeParametricCube.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkCubeSource.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkFloatArray.h"
#include "vtkGradientFilter.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointDataToCellData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkTextureMapToSphere.h"
#include "vtkThreshold.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkMakeParametricCube, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkMakeParametricCube);


//---------------------------------------------------------------------------
vtkMakeParametricCube::vtkMakeParametricCube()
{
  this->SetNumberOfInputPorts(0);
  this->Verbose = 1;

  this->FinalPd = vtkPolyData::New();
}

//---------------------------------------------------------------------------
vtkMakeParametricCube::~vtkMakeParametricCube()
{
  if (this->FinalPd != NULL)
  {
    this->FinalPd->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkMakeParametricCube::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkMakeParametricCube::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  if (this->MakeCube() != 1)
  {
    vtkErrorMacro("Making cube failed!");
    return 0;
  }
  if (this->LabelCube() != 1)
  {
    vtkErrorMacro("Labeling cube failed!");
    return 0;
  }

  output->DeepCopy(this->FinalPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkMakeParametricCube::MakeCube()
{
  vtkSmartPointer<vtkCubeSource> cubeMaker =
    vtkSmartPointer<vtkCubeSource>::New();
  cubeMaker->SetXLength(1.0);
  cubeMaker->SetYLength(1.0);
  cubeMaker->SetZLength(1.0);
  cubeMaker->SetCenter(0.5,0.5,0.5);
  cubeMaker->Update();

  vtkSmartPointer<vtkTriangleFilter> triangulator =
    vtkSmartPointer<vtkTriangleFilter>::New();
  triangulator->SetInputData(cubeMaker->GetOutput());
  triangulator->Update();

  vtkSmartPointer<vtkLinearSubdivisionFilter> subdivider =
    vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
  subdivider->SetInputData(triangulator->GetOutput());
  subdivider->SetNumberOfSubdivisions(2);
  subdivider->Update();

  vtkSmartPointer<vtkConnectivityFilter> connector =
    vtkSmartPointer<vtkConnectivityFilter>::New();
  connector->SetInputData(subdivider->GetOutput());
  connector->SetExtractionMode(VTK_EXTRACT_ALL_REGIONS);
  connector->ColorRegionsOn();
  connector->Update();

  vtkSmartPointer<vtkThreshold> thresholder =
    vtkSmartPointer<vtkThreshold>::New();
  thresholder->SetInputData(connector->GetOutput());
  thresholder->SetInputArrayToProcess(0,0,0,0,"RegionId");
  thresholder->ThresholdBetween(0,3);
  thresholder->Update();

  vtkSmartPointer<vtkDataSetSurfaceFilter> surfacer =
    vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  vtkSmartPointer<vtkCleanPolyData> cleaner =
    vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner->SetInputData(surfacer->GetOutput());
  cleaner->Update();

  this->FinalPd->DeepCopy(cleaner->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkMakeParametricCube::LabelCube()
{
  int numPts = this->FinalPd->GetNumberOfPoints();
  vtkSmartPointer<vtkFloatArray> pCoords =
    vtkSmartPointer<vtkFloatArray>::New();
  pCoords->SetNumberOfComponents(2);
  pCoords->Allocate(numPts, 10000);
  pCoords->SetNumberOfTuples(numPts);
  pCoords->SetName("PCoords");

  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    double newTup[2];
    this->FinalPd->GetPoint(i, pt);
    if (pt[1] == 0)
    {
      newTup[0] = pt[0];
      newTup[1] = pt[2];
    }
    else if (pt[0] == 1)
    {
      newTup[0] = 1.00 + pt[1];
      newTup[1] = pt[2];
    }
    else if (pt[1] == 1)
    {
      newTup[0] = 3.00 - pt[0];
      newTup[1] = pt[2];
    }
    else if (pt[0] == 0)
    {
      newTup[0] = 4.00 - pt[1];
      newTup[1] = pt[2];
    }
    else
    {
      newTup[0] = pt[0];
      newTup[1] = pt[1];
    }
    pCoords->SetTuple(i, newTup);
  }

  this->FinalPd->GetPointData()->AddArray(pCoords);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkMakeParametricCube::WriteToGroupsFile(vtkPolyData *pd, std::string fileName)
{
  vtkFloatArray *pCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("PCoords"));
  int numPts = pd->GetNumberOfPoints();

  double xSpacing, ySpacing;
  vtkMakeParametricCube::GetSpacingOfPCoords(pd, xSpacing, ySpacing);

  vtkSmartPointer<vtkIntArray> newPointOrder =
    vtkSmartPointer<vtkIntArray>::New();
  vtkMakeParametricCube::GetNewPointOrder(pd, xSpacing, ySpacing, newPointOrder);

  FILE *pFile;
  pFile = fopen(fileName.c_str(), "w");
  if (pFile == NULL)
  {
    fprintf(stderr,"Error opening file\n");
    return 0;
  }

  int xNum = 4.0/xSpacing;
  int yNum = 1.0/ySpacing + 1;
  fprintf(stdout,"XNum: %d\n", xNum);
  fprintf(stdout,"YNum: %d\n", yNum);
  for (int i=0; i<yNum; i++)
  {
    fprintf(pFile, "/group/test/%d\n", i);
    fprintf(pFile, "%d\n", i);
    fprintf(pFile, "center_x 0.0\n");
    for (int j=0; j<xNum; j++)
    {
      int ptId = newPointOrder->GetValue(i*xNum + j);
      double pt[3];
      pd->GetPoint(ptId, pt);
      fprintf(pFile, "%.6f %.6f %.6f\n", pt[0], pt[1], pt[2]);
    }
    fprintf(pFile, "\n");
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkMakeParametricCube::GetSpacingOfPCoords(vtkPolyData *pd, double &xSpacing, double &ySpacing)
{
  vtkFloatArray *pCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("PCoords"));
  int numPts = pd->GetNumberOfPoints();

  double xMin = 1.0e9;
  double yMin = 1.0e9;
  for (int i=0; i<numPts; i++)
  {
    double ptVal[2];
    pCoords->GetTuple(i, ptVal);
    if (ptVal[0] < xMin && ptVal[0] > 1e-8)
    {
      xMin = ptVal[0];
    }
    if (ptVal[1] < yMin && ptVal[1] > 1e-8)
    {
      yMin = ptVal[1];
    }
  }

  xSpacing = xMin;
  ySpacing = yMin;
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkMakeParametricCube::GetNewPointOrder(vtkPolyData *pd, double xSpacing, double ySpacing,
                                         vtkIntArray *newPointOrder)
{
  vtkFloatArray *pCoords = vtkFloatArray::SafeDownCast(pd->GetPointData()->GetArray("PCoords"));
  int numPts = pd->GetNumberOfPoints();

  int xNum = 4.0/xSpacing;
  for (int i=0; i<numPts; i++)
  {
    double ptVal[2];
    pCoords->GetTuple(i, ptVal);

    int xLoc = ptVal[0]/xSpacing;
    int loc = xLoc + xNum * ptVal[1]/ySpacing;

    double pt[3];
    pd->GetPoint(i, pt);
    fprintf(stdout,"Pt Val: %.4f %.4f, Loc: %d end\n", ptVal[0], ptVal[1], loc);
    fprintf(stdout,"Point: %.6f %.6f %.6f\n", pt[0], pt[1], pt[2]);
    newPointOrder->InsertValue(loc, i);
  }
  return 1;
}
