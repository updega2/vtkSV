/*=========================================================================
 *
 * Copyright (c) 2014-2015 The Regents of the University of California.
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

#include "vtkSVCenterlinesEdgeWeightedCVT.h"

#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVMathUtils.h"
#include "vtkTriangle.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVCenterlinesEdgeWeightedCVT);

// ----------------------
// Constructor
// ----------------------
vtkSVCenterlinesEdgeWeightedCVT::vtkSVCenterlinesEdgeWeightedCVT()
{
  this->DistanceFunction = vtkSVPolyBallLine::New();

  this->GroupIdsArrayName = NULL;
  this->BlankingArrayName = NULL;
  this->CenterlineRadiusArrayName =   NULL;

  this->EdgeWeight = 1.0;
  this->UseRadiusInformation = 1;
  this->UseBifurcationInformation = 1;
  this->UseCurvatureWeight = 1;
  this->UsePointNormal = 1;
}

// ----------------------
// Destructor
// ----------------------
vtkSVCenterlinesEdgeWeightedCVT::~vtkSVCenterlinesEdgeWeightedCVT()
{
  if (this->DistanceFunction != NULL)
  {
    this->DistanceFunction->Delete();
    this->DistanceFunction = NULL;
  }

  if (this->GroupIdsArrayName != NULL)
  {
    delete [] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
  }
  if (this->BlankingArrayName != NULL)
  {
    delete [] this->BlankingArrayName;
    this->BlankingArrayName = NULL;
  }
  if (this->CenterlineRadiusArrayName != NULL)
  {
    delete [] this->CenterlineRadiusArrayName;
    this->CenterlineRadiusArrayName = NULL;
  }
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVCenterlinesEdgeWeightedCVT::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// ----------------------
// InitializeConnectivity
// ----------------------
int vtkSVCenterlinesEdgeWeightedCVT::InitializeConnectivity()
{
  this->Superclass::InitializeConnectivity();

  int numCells = this->WorkGenerators->GetNumberOfCells();
  int numPoints = this->WorkGenerators->GetNumberOfPoints();
  this->IsGoodNeighborCell.resize(numCells, std::vector<int>(numCells));

  if (this->UseBifurcationInformation)
  {
    for (int i=0; i<numPoints; i++)
    {
      vtkNew(vtkIdList, pointCells);
      this->WorkGenerators->GetPointCells(i, pointCells);
      if (pointCells->GetNumberOfIds() > 1)
        this->FindGoodCellNeighbors(i, pointCells);
    }
  }
  else
  {
    for (int i=0; i<numCells; i++)
    {
      for (int j=0; j<numCells; j++)
          this->IsGoodNeighborCell[i][j] = 1;
    }
  }


  return SV_OK;
}

// ----------------------
// InitializeGenerators
// ----------------------
int vtkSVCenterlinesEdgeWeightedCVT::InitializeGenerators()
{
  if (this->UseCellArray)
  {
    // Tube function
    this->DistanceFunction->SetInput(this->WorkGenerators);
    this->DistanceFunction->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
    this->DistanceFunction->SetUseRadiusInformation(this->UseRadiusInformation);
    this->DistanceFunction->SetUsePointNormal(this->UsePointNormal);
    this->DistanceFunction->SetUseBifurcationInformation(this->UseBifurcationInformation);
    //this->DistanceFunction->BuildLocator();

    // Get all the different ids
    vtkNew(vtkIdList, centerlineGroupIds);
    for (int i=0; i<this->WorkGenerators->GetCellData()->GetArray(this->GroupIdsArrayName)->GetNumberOfTuples(); i++)
    {
      centerlineGroupIds->InsertUniqueId(static_cast<vtkIdType>(vtkMath::Round(this->WorkGenerators->GetCellData()->GetArray(this->GroupIdsArrayName)->GetComponent(i,0))));
    }
    int numGenerators = centerlineGroupIds->GetNumberOfIds();

    // Loop through cells
    int numCells = this->WorkPd->GetNumberOfCells();
    for (int i=0; i<numCells; i++)
    {
      // Get cell point coords
      double pts[3][3];
      vtkIdType npts, *ptids;
      this->WorkPd->GetCellPoints(i, npts, ptids);
      for (int j=0; j<npts; j++)
        this->WorkPd->GetPoint(ptids[j], pts[j]);

      // Get center
      double center[3];
      vtkTriangle::TriangleCenter(pts[0], pts[1], pts[2], center);

      if (this->UsePointNormal)
      {
        double pointNormal[3];
        this->CVTDataArray->GetTuple(i, pointNormal);
        this->DistanceFunction->SetPointNormal(pointNormal);
      }

      int cellGenerator = 0;
      double minDist = VTK_SV_LARGE_DOUBLE;

      vtkNew(vtkIdList, groupId);
      groupId->SetNumberOfIds(1);

      double dist = this->DistanceFunction->EvaluateFunction(center);
      int lastCellId = this->DistanceFunction->GetLastPolyBallCellId();
      cellGenerator = this->WorkGenerators->GetCellData()->GetArray(
          this->GroupIdsArrayName)->GetTuple1(lastCellId);

      this->PatchIdsArray->SetTuple1(i, cellGenerator);
    }
  }
  else if (this->UsePointArray)
  {
    vtkErrorMacro("Not implemented");
    return SV_ERROR;
  }
  return SV_OK;
}

// ----------------------
// UpdateGenerators
// ----------------------
int vtkSVCenterlinesEdgeWeightedCVT::UpdateGenerators()
{
  return SV_OK;
}

// ----------------------
// GetClosestGenerator
// ----------------------
int vtkSVCenterlinesEdgeWeightedCVT::GetClosestGenerator(const int evalId, int &newGenerator)
{
  // Get current generator
  int numGenerators = this->WorkGenerators->GetNumberOfPoints();
  int currGenerator = this->PatchIdsArray->GetTuple1(evalId);
  newGenerator =  currGenerator;

  // GroupIds
  vtkDataArray *groupIds = this->WorkGenerators->GetCellData()->GetArray(this->GroupIdsArrayName);

  // Current minimum to beat is current generator
  double minDist = this->GetEdgeWeightedDistance(newGenerator, evalId);

  // Loop through neighboring patches
  for (int i=0; i<this->NumberOfNeighborPatches[evalId]; i++)
  {
    // Check to make sure not zero elements or the same genrator
    if (this->NeighborPatchesNumberOfElements[evalId][i] != 0)
    {
      int cellId =         groupIds->LookupValue(currGenerator);
      int neighborCellId = groupIds->LookupValue(this->NeighborPatchesIds[evalId][i]);
      if (currGenerator != this->NeighborPatchesIds[evalId][i] && this->IsGoodNeighborCell[cellId][neighborCellId] == 1)
      {
        // Test this generator
        int neighborGenerator = this->NeighborPatchesIds[evalId][i];
        double testDist = this->GetEdgeWeightedDistance(neighborGenerator, evalId);

        // Set new min if less than current
        if (testDist < minDist)
        {
          minDist = testDist;
          newGenerator = neighborGenerator;
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// GetEdgeWeightedDistance
// ----------------------
double vtkSVCenterlinesEdgeWeightedCVT::GetEdgeWeightedDistance(const int generatorId, const int evalId)
{
  // Current generator
  int currGenerator = this->PatchIdsArray->GetTuple1(evalId);

  // TODO CHECK FOR NORMALS EARLIER!!!!
  // Current cell normal
  vtkDataArray *cellNormals = this->WorkPd->GetCellData()->GetArray("Normals");
  double currNormal[3];
  cellNormals->GetTuple(evalId, currNormal);

  // Get cell point coords
  double pts[3][3];
  vtkIdType npts, *ptids;
  this->WorkPd->GetCellPoints(evalId, npts, ptids);
  for (int j=0; j<npts; j++)
    this->WorkPd->GetPoint(ptids[j], pts[j]);

  // Get center
  double center[3];
  vtkTriangle::TriangleCenter(pts[0], pts[1], pts[2], center);

  // Calculate the edge weight distance
  double edgeWeightedDist = 1.0;

  // Get the generator patch id
  int i;
  for (i=0; i<this->NumberOfNeighborPatches[evalId]; i++)
  {
    if (this->NeighborPatchesIds[evalId][i] == generatorId)
    {
      break;
    }
  }

  double totalWeight = 0.0;
  for (int i=0; i<this->NumberOfNeighbors[evalId]; i++)
  {
    int neighborId = this->Neighbors[evalId][i];
    int neighborGenerator = this->PatchIdsArray->GetTuple1(neighborId);
    if (neighborGenerator == generatorId)
    {
      double normal[3];
      cellNormals->GetTuple(neighborId, normal);

      double crossVec[3];
      vtkMath::Cross(currNormal, normal, crossVec);
      double ang = atan2(vtkMath::Norm(crossVec), vtkMath::Dot(currNormal, normal));

      totalWeight += ang/SV_PI;
    }
  }
  if (this->UseCurvatureWeight)
    this->EdgeWeight = totalWeight/this->NeighborPatchesNumberOfElements[evalId][i];

  // Get the edge weighted portion
  double edgeWeighting = 0.0;
  if (currGenerator == generatorId)
  {
    edgeWeighting = 2 * this->EdgeWeight * (this->NumberOfNeighbors[evalId] - this->NeighborPatchesNumberOfElements[evalId][i]);
  }
  else
  {
    edgeWeighting = 2 * this->EdgeWeight * (this->NumberOfNeighbors[evalId] - this->NeighborPatchesNumberOfElements[evalId][i] - 1);
  }
  //}

  // Divide by the number of neighboring cells
  edgeWeighting /= this->NumberOfNeighbors[evalId];

  // Get the final edge distance to be normalized
  edgeWeightedDist += edgeWeighting;

  edgeWeightedDist = sqrt(edgeWeightedDist);

  return edgeWeightedDist;
  return SV_OK;
}

// ----------------------
// FindGoodCellNeighbors
// ----------------------
int vtkSVCenterlinesEdgeWeightedCVT::FindGoodCellNeighbors(const int ptId,
                                                           vtkIdList *cellIds)
{
  vtkNew(vtkIdList, checkIds);
  vtkNew(vtkIdList, cellList);

  for (int i=0; i<cellIds->GetNumberOfIds(); i++)
  {
    vtkIdType npts, *pts;
    this->WorkGenerators->GetCellPoints(cellIds->GetId(i), npts, pts);
    for (int j=0; j<npts; j++)
    {
      if (pts[j] == ptId)
      {
        if (j == 0)
        {
          checkIds->InsertNextId(pts[j+1]);
          cellList->InsertNextId(cellIds->GetId(i));
        }
        if (j == npts-1)
        {
          checkIds->InsertNextId(pts[j-1]);
          cellList->InsertNextId(cellIds->GetId(i));
        }
      }
    }
  }

  double pt0[3];
  this->WorkGenerators->GetPoint(ptId, pt0);

  double bigAngle = -1.0;
  int bigPoint = 0;
  for (int i=0; i<checkIds->GetNumberOfIds(); i++)
  {
    double pt1[3];
    this->WorkGenerators->GetPoint(checkIds->GetId(i), pt1);

    double vec0[3];
    vtkMath::Subtract(pt1, pt0, vec0);
    vtkMath::Normalize(vec0);

    double angleSum = 0.0;
    for (int j=0; j<checkIds->GetNumberOfIds(); j++)
    {
      if (i != j)
      {
        double pt2[3];
        this->WorkGenerators->GetPoint(checkIds->GetId(j), pt2);

        double vec1[3];
        vtkMath::Subtract(pt2, pt0, vec1);
        vtkMath::Normalize(vec1);

        double crossVec[3];
        vtkMath::Cross(vec0, vec1, crossVec);

        double ang = 180.0*atan2(vtkMath::Norm(crossVec), vtkMath::Dot(vec0, vec1))/SV_PI;
        angleSum += ang;
      }
    }
    if (angleSum > bigAngle)
    {
      bigAngle = angleSum;
      bigPoint = i;
    }
  }

  std::vector<double> angles;
  std::vector<int> cells;
  std::vector<int> points;

  double pt1[3];
  this->WorkGenerators->GetPoint(checkIds->GetId(bigPoint), pt1);
  angles.push_back(180.0);
  cells.push_back(cellList->GetId(bigPoint));
  points.push_back(checkIds->GetId(bigPoint));

  double vec0[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  for (int i=0; i<checkIds->GetNumberOfIds(); i++)
  {
    if (i != bigPoint)
    {
      double pt2[3];
      this->WorkGenerators->GetPoint(checkIds->GetId(i), pt2);

      double vec1[3];
      vtkMath::Subtract(pt2, pt0, vec1);
      vtkMath::Normalize(vec1);

      double crossVec[3];
      vtkMath::Cross(vec0, vec1, crossVec);

      double ang = 180.0*atan2(vtkMath::Norm(crossVec), vtkMath::Dot(vec0, vec1))/SV_PI;
      angles.push_back(ang);
      cells.push_back(cellList->GetId(i));
      points.push_back(checkIds->GetId(i));
    }
  }

  int rootCellId = cells[0];
  for (int i=0; i<angles.size(); i++)
  {
    int cellId = cells[i];
    if (angles[i] < 135.0)
    {
      this->IsGoodNeighborCell[rootCellId][cellId] = 1;
      this->IsGoodNeighborCell[cellId][rootCellId] = 1;
    }
    else
    {
      this->IsGoodNeighborCell[rootCellId][cellId] = 0;
      this->IsGoodNeighborCell[cellId][rootCellId] = 0;
      for (int j=1; j<angles.size(); j++)
      {
        if (i != j)
        {
          int otherCellId = cells[j];
          this->IsGoodNeighborCell[otherCellId][cellId] = 1;
          this->IsGoodNeighborCell[cellId][otherCellId] = 1;
        }
      }
    }
  }

  return SV_OK;
}
