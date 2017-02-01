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

/** @file vtkPullApartPolyData.cxx
 *  @brief This implements the vtkPullApartPolyData filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkPullApartPolyData.h"

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

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkPullApartPolyData, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkPullApartPolyData);


//---------------------------------------------------------------------------
vtkPullApartPolyData::vtkPullApartPolyData()
{
  this->SetNumberOfInputPorts(1);

  this->CutPointsArrayName = NULL;

  this->WorkPd       = vtkPolyData::New();
  this->EdgeTable    = vtkEdgeTable::New();
  this->SeamPointIds = NULL;

  this->ReplacePointList = vtkIdList::New();
  this->NewPointList     = vtkIdList::New();

  this->StartPtId = -1;
}

//---------------------------------------------------------------------------
vtkPullApartPolyData::~vtkPullApartPolyData()
{
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->EdgeTable != NULL)
  {
    this->EdgeTable->Delete();
    this->EdgeTable = NULL;
  }
  if (this->SeamPointIds != NULL)
  {
    this->SeamPointIds->UnRegister(this);
    this->SeamPointIds = NULL;
  }
  if (this->ReplacePointList != NULL)
  {
    this->ReplacePointList->Delete();
    this->ReplacePointList = NULL;
  }
  if (this->NewPointList != NULL)
  {
    this->NewPointList->Delete();
    this->NewPointList = NULL;
  }

  if (this->CutPointsArrayName)
  {
    delete [] this->CutPointsArrayName;
    this->CutPointsArrayName = NULL;
  }
}

//---------------------------------------------------------------------------
void vtkPullApartPolyData::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkPullApartPolyData::RequestData(
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

  output->DeepCopy(this->WorkPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPullApartPolyData::PrepFilter()
{
  vtkIdType numPolys  = this->WorkPd->GetNumberOfPolys();
  vtkIdType numPoints = this->WorkPd->GetNumberOfPoints();
  //Check the input to make sure it is there
  if (numPolys < 1)
  {
    vtkErrorMacro("No input!");
    return 0;
  }

  // Check if dijkstra array name is given
  if (!this->CutPointsArrayName)
  {
    vtkDebugMacro("Cut Points Array Name not given, setting to CutPoints");
    this->CutPointsArrayName = new char[strlen("CutPoints") + 1];
    strcpy(this->CutPointsArrayName, "CutPoints");
  }
  // Check if array dijkstra is already on pd
  if (!this->CheckArrayExists(this->WorkPd, 0, this->CutPointsArrayName))
  {
    if (this->SeamPointIds != NULL)
    {
      vtkNew(vtkIntArray, seamArray);
      seamArray->SetNumberOfComponents(1);
      seamArray->SetNumberOfTuples(numPoints);
      seamArray->FillComponent(0, 0);
      seamArray->SetName(this->CutPointsArrayName);
      for (int i=0; i<this->SeamPointIds->GetNumberOfValues(); i++)
      {
        seamArray->SetValue(this->SeamPointIds->GetValue(i), 1);
      }
      this->WorkPd->GetPointData()->AddArray(seamArray);
    }
    else
    {
      vtkErrorMacro("No point array on surface named " << this->CutPointsArrayName << " on surface and not SeamPointIds provided. Must provide one or the other to cut pd");
      return 0;
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
int vtkPullApartPolyData::RunFilter()
{
  // Create edge table with neighbors
  if (!this->FindEdgeCells())
  {
    vtkErrorMacro("Failed finding edge cells");
    return 0;
  }

  if (!this->PullApartCutEdges())
  {
    vtkErrorMacro("Failed cutting edges");
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
int vtkPullApartPolyData::PullApartCutEdges()
{
  this->WorkPd->DeleteLinks();
  int numPoints = this->ReplacePointVector.size();
  vtkDebugMacro("Adding " << numPoints << "points");

  vtkNew(vtkPointData, newPointData);
  newPointData->CopyAllocate(this->WorkPd->GetPointData(),
    this->WorkPd->GetNumberOfPoints() + numPoints);
  for (int i=0; i<this->WorkPd->GetNumberOfPoints(); i++)
  {
    newPointData->CopyData(this->WorkPd->GetPointData(), i, i);
  }

  double shift[3]; shift[0] = 1.0e-6; shift[1] = 0.0; shift[2] = 0.0;
  for (int i=0; i<numPoints; i++)
  {
    double pt0[3];
    int replacePtId = this->ReplacePointVector[i];
    int newPtId = this->ReplacePointList->IsId(replacePtId);
    this->WorkPd->GetPoint(replacePtId, pt0);
    double newPt[3];
    vtkMath::Add(pt0, shift, newPt);
    if (newPtId == -1)
    {
      newPtId = this->WorkPd->GetPoints()->InsertNextPoint(newPt);
      this->ReplacePointList->InsertNextId(replacePtId);
      this->NewPointList->InsertNextId(newPtId);
    }
    else
    {
      newPtId = this->NewPointList->GetId(newPtId);
    }
    newPointData->CopyData(this->WorkPd->GetPointData(), replacePtId, newPtId);

    int numCells = this->ReplaceCellVector[i].size();
    for (int j=0; j<numCells; j++)
    {
      this->WorkPd->ReplaceCellPoint(this->ReplaceCellVector[i][j], replacePtId, newPtId);
    }
  }

  this->WorkPd->GetPointData()->PassData(newPointData);

  this->WorkPd->BuildLinks();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPullApartPolyData::FindEdgeCells()
{
  int numPts = this->WorkPd->GetNumberOfPoints();
  this->WorkPd->BuildLinks();

  vtkIntArray *cutPointValues = vtkIntArray::SafeDownCast(
    this->WorkPd->GetPointData()->GetArray(this->CutPointsArrayName));

  this->EdgeTable->InitEdgeInsertion(numPts, 1);

  int startPt0, startPt1, startPt2, startCellId;
  if (this->StartPtId != -1)
  {
    vtkNew(vtkIdList, pointCells);
    this->WorkPd->GetPointCells(this->StartPtId, pointCells);
    for (int i=0; i<pointCells->GetNumberOfIds(); i++)
    {
      vtkIdType npts, *pts;
      int cellId = pointCells->GetId(i);
      this->WorkPd->GetCellPoints(cellId, npts, pts);
      for (int j=0; j<npts; j++)
      {
        startPt0 = pts[j];
        startPt1 = pts[(j+1)%npts];
        startPt2 = pts[(j+2)%npts];
        vtkNew(vtkIdList, neighborCells);
        this->WorkPd->GetCellEdgeNeighbors(cellId, startPt0, startPt1, neighborCells);
        // We found the edge cell!
        if (neighborCells->GetNumberOfIds() == 0)
        {
          double pt0[3], pt1[3], vec0[3], vec1[3];
          this->WorkPd->GetPoint(startPt0, pt0);
          this->WorkPd->GetPoint(startPt1, pt1);
          if (startPt0 == this->StartPtId)
          {
            vtkMath::Subtract(pt1, pt0, vec0);
          }
          else
          {
            vtkMath::Subtract(pt0, pt1, vec0);
            int tmp = startPt1;
            startPt1 = startPt0;
            startPt0 = tmp;
          }
          vtkMath::Cross(this->ObjectZAxis, this->ObjectXAxis, vec1);
          vtkMath::Normalize(vec0);
          vtkMath::Normalize(vec1);
          // We found the right cell!
          if (vtkMath::Dot(vec0, vec1) < 0)
          {
            startCellId = cellId;
            std::vector<int> firstCellList; firstCellList.push_back(startCellId);
            if (cutPointValues->GetValue(startPt2) == 1)
            {
              int tmp = startPt1;
              startPt1 = startPt2;
              startPt2 = tmp;
              //fprintf(stdout,"Found em!: %d %d %d\n", startPt0, startPt1, startPt2);
              this->ReplaceCellVector.push_back(firstCellList);
              this->ReplacePointVector.push_back(startPt0);
              std::vector<int> list0; list0.push_back(startCellId);
              this->FindNextEdge(startPt0, startPt1, startPt2, startCellId, list0, 1);
            }
            else
            {
              int tmp = startPt1;
              startPt1 = startPt0;
              startPt0 = tmp;
              //fprintf(stdout,"Found em!: %d %d %d\n", startPt0, startPt1, startPt2);
              this->FindNextEdge(startPt0, startPt1, startPt2, startCellId, firstCellList, 1);
            }
          }
        }
      }
    }
  }
  else
  {
    if (this->FindStartingEdge(startPt0, startPt1, startPt2, startCellId) != 1)
    {
      vtkErrorMacro("Starting edge could not be found");
      return 0;
    }
    //fprintf(stdout,"Starts: %d %d %d\n", startPt0, startPt1, startPt2);
    std::vector<int> list0; list0.push_back(startCellId);
    std::vector<int> list1; list1.push_back(startCellId);
    this->FindNextEdge(startPt0, startPt1, startPt2, startCellId, list0, 1);
    this->FindNextEdge(startPt1, startPt0, startPt2, startCellId, list1, 1);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPullApartPolyData::FindStartingEdge(int &p0, int &p1, int &p2, int &cellId)
{
  int numTris = this->WorkPd->GetNumberOfCells();
  vtkIntArray *cutPointValues = vtkIntArray::SafeDownCast(
    this->WorkPd->GetPointData()->GetArray(this->CutPointsArrayName));

  for (int i=0; i<numTris; i++)
  {
    vtkIdType npts, *pts;
    this->WorkPd->GetCellPoints(i, npts, pts);

    //Insert edge into table
    for (int j=0; j<npts; j++)
    {
      p0 = pts[j];
      p1 = pts[(j+1)%npts];
      p2 = pts[(j+2)%npts];

      vtkIdType checkEdge = this->EdgeTable->IsEdge(p0, p1);
      if (checkEdge == -1)
      {

        int cutVal0 = cutPointValues->GetValue(p0);
        int cutVal1 = cutPointValues->GetValue(p1);
        if (cutVal0 == 1 && cutVal1 == 1)
        {
          vtkIdType edgeId = this->EdgeTable->InsertEdge(p0, p1);
          cellId = i;
          vtkDebugMacro("First edge in list with points" << p0 << " and " << p1 << " on cell " << cellId);
          //fprintf(stdout,"First Edge in list with points %d and %d on cell %d\n", p0, p1, i);
          return 1;
        }
      }
    }
  }

  return 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPullApartPolyData::FindNextEdge(int p0, int p1, int p2, int cellId, std::vector<int> &cellList, int first)
{
  vtkIntArray *cutPointValues = vtkIntArray::SafeDownCast(
    this->WorkPd->GetPointData()->GetArray(this->CutPointsArrayName));

  vtkNew(vtkIdList, neighborCells);
  this->WorkPd->GetCellEdgeNeighbors(cellId, p1, p2, neighborCells);
  if (neighborCells->GetNumberOfIds() == 1)
  {
    vtkIdType npts, *pts;
    int neighborCell = neighborCells->GetId(0);
    this->WorkPd->GetCellPoints(neighborCell, npts, pts);
    for (int j=0; j<npts; j++)
    {
      int cutVal = cutPointValues->GetValue(pts[j]);
      if (pts[j] != p1 && pts[j] != p2 && cutVal != 1)
      {
        int outP = pts[j];
        cellList.push_back(neighborCell);
        vtkDebugMacro("Edge in same list with points" << p2 << " and " << p1 << " on cell " << neighborCell);
        //fprintf(stdout,"Edge in same list with points %d and %d on cell %d\n", p2, p1, neighborCell);
        this->FindNextEdge(p2, p1, outP, neighborCell, cellList, 0);
      }
      else if (pts[j] != p1 && pts[j] != p2 && cutVal == 1)
      {
        // Found end
        int newP = pts[j];
        cellList.push_back(neighborCell);
        vtkDebugMacro("Final edge in list with points" << p1 << " and " << newP << " on cell " << neighborCell);
        //fprintf(stdout,"Final Edge in list with points %d and %d on cell %d\n", p1, newP, neighborCell);
        this->ReplaceCellVector.push_back(cellList);
        this->ReplacePointVector.push_back(p1);
        std::vector<int> newCellList;
        newCellList.push_back(neighborCell);
        vtkDebugMacro("First edge in list with points" << p1 << " and " << newP << " on cell " << neighborCell);
        //fprintf(stdout,"First Edge in list with points %d and %d on cell %d\n", p1, newP, neighborCell);
        vtkIdType edgeId = this->EdgeTable->InsertEdge(p1, newP);
        this->FindNextEdge(p1, newP, p2, neighborCell, newCellList, 1);
        this->ReplaceCellVector.push_back(newCellList);
        this->ReplacePointVector.push_back(newP);

      }
    }
  }
  else if (neighborCells->GetNumberOfIds() == 0 && first)
  {
    //fprintf(stdout,"Ending edge with point %d on cell %d\n", p1, cellId);
    std::vector<int> newCellList;
    newCellList.push_back(cellId);
    this->ReplaceCellVector.push_back(newCellList);
    this->ReplacePointVector.push_back(p1);
  }
  else if (neighborCells->GetNumberOfIds() == 0 && !first)
  {
    this->ReplaceCellVector.push_back(cellList);
    this->ReplacePointVector.push_back(p1);
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

int vtkPullApartPolyData::CheckArrayExists(vtkPolyData *pd,
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
