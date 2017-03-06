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

/** @file vtkSVFindGeodesicPath.cxx
 *  @brief This implements the vtkSVFindGeodesicPath filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVFindGeodesicPath.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDijkstraGraphGeodesicPath.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDoubleArray.h"
#include "vtkFeatureEdges.h"
#include "vtkFloatArray.h"
#include "vtkIdFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataWriter.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkSVFindGeodesicPath, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSVFindGeodesicPath);


//---------------------------------------------------------------------------
vtkSVFindGeodesicPath::vtkSVFindGeodesicPath()
{
  this->SetNumberOfInputPorts(1);
  this->Verbose                  = 0;
  this->AddPathBooleanArray      = 0;
  this->RemoveInternalIds        = 1;
  this->RepelCloseBoundaryPoints = 0;

  this->StartPtId = -1;
  this->EndPtId   = -1;

  for (int i=0; i<3; i++)
  {
    this->ClosePt[i] = 0.0;
  }

  this->DijkstraArrayName = NULL;
  this->InternalIdsArrayName = NULL;
  this->PathBooleanArrayName = NULL;

  this->WorkPd       = vtkPolyData::New();
  this->Boundary     = vtkPolyData::New();
  this->PathIds      = vtkIdList::New();
  this->PathBoolean  = vtkIntArray::New();
}

//---------------------------------------------------------------------------
vtkSVFindGeodesicPath::~vtkSVFindGeodesicPath()
{
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->PathIds != NULL)
  {
    this->PathIds->Delete();
    this->PathIds = NULL;
  }
  if (this->Boundary != NULL)
  {
    this->Boundary->Delete();
    this->Boundary = NULL;
  }
  if (this->PathBoolean != NULL)
  {
    this->PathBoolean->Delete();
    this->PathBoolean = NULL;
  }

  if (this->DijkstraArrayName)
  {
    delete [] this->DijkstraArrayName;
    this->DijkstraArrayName = NULL;
  }
  if (this->InternalIdsArrayName)
  {
    delete [] this->InternalIdsArrayName;
    this->InternalIdsArrayName = NULL;
  }
  if (this->PathBooleanArrayName)
  {
    delete [] this->PathBooleanArrayName;
    this->PathBooleanArrayName = NULL;
  }
}

//---------------------------------------------------------------------------
void vtkSVFindGeodesicPath::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkSVFindGeodesicPath::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // Get the input and output
  vtkPolyData *input  = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  // Copy the input to operate on
  this->WorkPd->DeepCopy(input);

  // Prep work for filter
  if (this->PrepFilter() != 1)
  {
    vtkErrorMacro("Prep of filter failed");
    output->DeepCopy(input);
    return 0;
  }

  // Run the filter
  if (this->RunFilter() != 1)
  {
    vtkErrorMacro("Filter failed");
    output->DeepCopy(input);
    return 0;
  }

  if (this->RemoveInternalIds)
  {
    this->WorkPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);
    this->WorkPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
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
int vtkSVFindGeodesicPath::PrepFilter()
{
  vtkIdType numPolys  = this->WorkPd->GetNumberOfPolys();
  vtkIdType numPoints = this->WorkPd->GetNumberOfPoints();
  //Check the input to make sure it is there
  if (numPolys < 1)
  {
    vtkErrorMacro("No input!");
    return 0;
  }

  // Check is start point id is given
  if (this->StartPtId == -1)
  {
    vtkErrorMacro("No input start id given");
    return 0;
  }
  if (this->StartPtId > numPoints)
  {
    vtkErrorMacro("Start id is greater than number of pts on pd");
    return 0;
  }
  if (this->EndPtId > numPoints)
  {
    vtkErrorMacro("End id is greater than number of pts on pd");
    return 0;
  }

  // Check if dijkstra array name is given
  if (!this->DijkstraArrayName)
  {
    vtkDebugMacro("Dijkstra Array Name not given, setting to DijkstraDistance");
    this->DijkstraArrayName = new char[strlen("DijkstraDistance") + 1];
    strcpy(this->DijkstraArrayName, "DijkstraDistance");
  }
  // Check if array dijkstra is already on pd
  if (this->CheckArrayExists(this->WorkPd, 0, this->DijkstraArrayName))
  {
    this->WorkPd->GetPointData()->RemoveArray(this->DijkstraArrayName);
  }

  // Check if internal id array name is given
  if (!this->InternalIdsArrayName)
  {
    vtkDebugMacro("Internal Ids Array Name not given, setting to InternalIds");
    this->InternalIdsArrayName = new char[strlen("InternalIds") + 1];
    strcpy(this->InternalIdsArrayName, "InternalIds");
  }
  // Check if array internal ids is already on pd
  if (this->CheckArrayExists(this->WorkPd, 0, this->InternalIdsArrayName))
  {
    this->RemoveInternalIds = 0;
  }
  else
  {
    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(this->WorkPd);
    ider->SetIdsArrayName(this->InternalIdsArrayName);
    ider->Update();
    this->WorkPd->DeepCopy(ider->GetOutput());
  }

  // Check if path boolean array name is given
  if (!this->PathBooleanArrayName)
  {
    vtkDebugMacro("PathBoolean Array Name not given, setting to PathBoolean");
    this->PathBooleanArrayName = new char[strlen("PathBoolean") + 1];
    strcpy(this->PathBooleanArrayName, "PathBoolean");
  }
  // Check if array path booleana is already on pd
  if (this->CheckArrayExists(this->WorkPd, 0, this->PathBooleanArrayName))
  {
    this->WorkPd->GetPointData()->RemoveArray(this->PathBooleanArrayName);
  }


  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVFindGeodesicPath::RunFilter()
{
  int runItChrisBrown = 0;
  if (this->EndPtId != -1 || this->AddPathBooleanArray)
  {
    runItChrisBrown = 1;
  }

  if (this->EndPtId == -1)
  {
    if (this->FindClosestBoundaryPoint() != 1)
    {
      vtkErrorMacro("Error finding a point close on the boundary");
      return 0;
    }
  }
  vtkNew(vtkPoints, repelPoints);
  if (this->RepelCloseBoundaryPoints)
  {
    if (this->GetCloseBoundaryPoints(this->StartPtId, this->EndPtId, repelPoints) != 1)
    {
      vtkErrorMacro("Error getting close points on the boundary to repel");
      return 0;
    }
  }

  if (runItChrisBrown)
  {
    if (this->RunDijkstra(repelPoints) != 1)
    {
      vtkErrorMacro("vtkDijkstraGraphGeodesicPath failed");
      return 0;
    }
    if (this->AddPathBooleanArray)
    {
      int numPoints = this->WorkPd->GetNumberOfPoints();
      this->PathBoolean->SetNumberOfTuples(numPoints);
      this->PathBoolean->FillComponent(0, 0);
      for (int i=0; i<this->PathIds->GetNumberOfIds(); i++)
      {
        this->PathBoolean->SetValue(this->PathIds->GetId(i), 1);
      }
      this->PathBoolean->SetName(this->PathBooleanArrayName);
      this->WorkPd->GetPointData()->AddArray(this->PathBoolean);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVFindGeodesicPath::FindClosestBoundaryPoint()
{
  if (this->RunDijkstra(NULL) != 1)
  {
    vtkErrorMacro("vtkDijkstraGraphGeodesicPath failed");
    return 0;
  }

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
  vtkDataArray *internalIds   = this->Boundary->GetPointData()->GetArray(this->InternalIdsArrayName);
  int numPoints = this->Boundary->GetNumberOfPoints();
  double minVal = 1.0e10;
  int minId = -1;
  for (int i=0; i<numPoints; i++)
  {
    double val = passedWeights->GetTuple1(i);
    if (val < minVal)
    {
      minVal = val;
      minId  = internalIds->GetTuple1(i);
    }
  }

  this->EndPtId = minId;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVFindGeodesicPath::RunDijkstra(vtkPoints *repelPoints)
{
  vtkNew(vtkDijkstraGraphGeodesicPath, dijkstra);
  dijkstra->SetInputData(this->WorkPd);
  dijkstra->SetStartVertex(this->StartPtId);
  if (repelPoints != NULL)
  {
    if (repelPoints->GetNumberOfPoints() != 0)
    {
      dijkstra->RepelPathFromVerticesOn();
      dijkstra->SetRepelVertices(repelPoints);
    }
  }
  if (this->EndPtId != -1)
  {
    dijkstra->SetEndVertex(this->EndPtId);
  }
  dijkstra->StopWhenEndReachedOff();
  dijkstra->Update();

  vtkNew(vtkDoubleArray, tmpWeights);
  dijkstra->GetCumulativeWeights(tmpWeights);
  tmpWeights->SetName(this->DijkstraArrayName);
  this->WorkPd->GetPointData()->AddArray(tmpWeights);
  this->PathIds->DeepCopy(dijkstra->GetIdList());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVFindGeodesicPath::GetCloseBoundaryPoints(const int startPtId,
                                                const int endPtId,
                                                vtkPoints *repelPoints)
{
  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(this->WorkPd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();

  double startPt[3], endPt[3];
  this->WorkPd->GetPoint(startPtId, startPt);
  this->WorkPd->GetPoint(endPtId, endPt);

  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(boundaries->GetOutput());
  connector->SetExtractionModeToClosestPointRegion();
  connector->SetClosestPoint(startPt);
  connector->Update();
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  if (this->GetNeighborBoundaryPoints(startPtId, surfacer->GetOutput(), repelPoints) != 1)
  {
    vtkErrorMacro("Error getting neighbor boundary points");
    return 0;
  }

  vtkNew(vtkConnectivityFilter, connector2);
  connector2->SetInputData(boundaries->GetOutput());
  connector2->SetExtractionModeToClosestPointRegion();
  connector2->SetClosestPoint(endPt);
  connector2->Update();
  vtkNew(vtkDataSetSurfaceFilter, surfacer2);
  surfacer2->SetInputData(connector2->GetOutput());
  surfacer2->Update();

  if (this->GetNeighborBoundaryPoints(endPtId, surfacer2->GetOutput(), repelPoints) != 1)
  {
    vtkErrorMacro("Error getting neighbor boundary points");
    return 0;
  }

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVFindGeodesicPath::GetNeighborBoundaryPoints(const int ptId,
                                                   vtkPolyData *pd,
                                                   vtkPoints *repelPoints)
{
  vtkDataArray *internalIds = this->WorkPd->GetPointData()->
    GetArray(this->InternalIdsArrayName);
  vtkDataArray *eInternalIds = pd->GetPointData()->
    GetArray(this->InternalIdsArrayName);

  vtkNew(vtkIdList, cells);
  this->WorkPd->GetPointCells(ptId, cells);
  int offLimits[2]; offLimits[0] = -1; offLimits[1] = -1;
  int count = 0;
  if (cells->GetNumberOfIds() == 1)
  {
    vtkIdType npts, *pts;
    this->WorkPd->GetCellPoints(cells->GetId(0), npts, pts);
    for (int j=0; j<npts; j++)
    {
      if (pts[j] != ptId)
        offLimits[count++] = eInternalIds->LookupValue(int(internalIds->GetTuple1(pts[j])));
    }
  }

  int bId = eInternalIds->LookupValue(int(internalIds->GetTuple1(ptId)));
  if (bId != -1)
  {
    if (bId >= pd->GetNumberOfPoints())
    {
      vtkErrorMacro("Point id is not valid " << ptId);
      return 0;
    }
    for (int i=0; i<pd->GetNumberOfPoints(); i++)
    {
      if (i != bId && i != offLimits[0] && i != offLimits[1])
      {
        double pt[3];
        pd->GetPoint(i, pt);
        repelPoints->InsertNextPoint(pt);
      }
    }
  }

  return 1;
}

// ----------------------
// CheckArrayExists
// ----------------------
/**
 * @brief Function to check is array with name exists in cell or point data
 * @param pd this is the object to check if the array exists
 * @param datatype this is point or cell. point=0,cell=1
 * @param arrayname this is the name of the array to check
 * @reutrn this returns 1 if the array exists and zero if it doesn't
 * or the function does not return properly.
 */

int vtkSVFindGeodesicPath::CheckArrayExists(vtkPolyData *pd,
                                          int datatype,
                                          std::string arrayname )
{
  int exists =0;

  if (datatype == 0)
  {
    int numArrays = pd->GetPointData()->GetNumberOfArrays();
    for (int i=0;i<numArrays;i++)
    {
      if (!strcmp(pd->GetPointData()->GetArrayName(i),arrayname.c_str()))
      {
	      exists =1;
      }
    }
  }
  else
  {
    int numArrays = pd->GetCellData()->GetNumberOfArrays();
    for (int i=0;i<numArrays;i++)
    {
      if (!strcmp(pd->GetCellData()->GetArrayName(i),arrayname.c_str()))
      {
	      exists =1;
      }
    }
  }

  return exists;
}
