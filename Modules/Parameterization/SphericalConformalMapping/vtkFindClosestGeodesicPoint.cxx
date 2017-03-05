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

/** @file vtkFindClosestGeodesicPoint.cxx
 *  @brief This implements the vtkFindClosestGeodesicPoint filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkFindClosestGeodesicPoint.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDijkstraGraphGeodesicPath.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDoubleArray.h"
#include "vtkFeatureEdges.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkFindClosestGeodesicPoint, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkFindClosestGeodesicPoint);


//---------------------------------------------------------------------------
vtkFindClosestGeodesicPoint::vtkFindClosestGeodesicPoint()
{
  this->SetNumberOfInputPorts(1);
  this->Verbose = 0;

  this->StartPtId = 0;
  this->EndPtId   = 0;

  for (int i=0; i<3; i++)
  {
    this->ClosePt[i] = 0.0;
  }

  this->DijkstraArrayName = NULL;

  this->WorkPd   = vtkPolyData::New();
  this->Boundary = vtkPolyData::New();
}

//---------------------------------------------------------------------------
vtkFindClosestGeodesicPoint::~vtkFindClosestGeodesicPoint()
{
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->Boundary != NULL)
  {
    this->Boundary->Delete();
    this->Boundary = NULL;
  }

  if (this->DijkstraArrayName)
  {
    delete [] this->DijkstraArrayName;
    this->DijkstraArrayName = NULL;
  }
}

//---------------------------------------------------------------------------
void vtkFindClosestGeodesicPoint::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkFindClosestGeodesicPoint::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input  = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->WorkPd->DeepCopy(input);

  vtkIdType numPolys = this->WorkPd->GetNumberOfPolys();
  //Check the input to make sure it is there
  if (numPolys < 1)
  {
    vtkDebugMacro("No input!");
    return 0;
  }

  if (this->FindClosestBoundaryPoint() != 1)
  {
    return 0;
  }

  output->DeepCopy(this->WorkPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkFindClosestGeodesicPoint::FindClosestBoundaryPoint()
{
  vtkNew(vtkDijkstraGraphGeodesicPath, dijkstra);
  dijkstra->SetInputData(this->WorkPd);
  dijkstra->SetStartVertex(this->StartPtId);
  dijkstra->StopWhenEndReachedOff();
  dijkstra->Update();

  vtkNew(vtkDoubleArray, tmpWeights);
  dijkstra->GetCumulativeWeights(tmpWeights);
  tmpWeights->SetName(this->DijkstraArrayName);
  this->WorkPd->GetPointData()->AddArray(tmpWeights);

  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(this->WorkPd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();

  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(boundaries->GetOutput());
  connector->SetExtractionModeToClosestPointRegion();
  connector->SetClosestPoint(this->ClosePt);
  connector->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  this->Boundary->ShallowCopy(surfacer->GetOutput());
  vtkDataArray *passedWeights = this->Boundary->GetPointData()->GetArray(this->DijkstraArrayName);
  int numPoints = this->Boundary->GetNumberOfPoints();
  double minVal = 1.0e10;
  int minId = -1;
  for (int i=0; i<numPoints; i++)
  {
    double val = passedWeights->GetTuple1(i);
    if (val < minVal)
    {
      minVal = val;
      minId  = i;
    }
  }

  this->EndPtId = minId;
  return 1;
}
