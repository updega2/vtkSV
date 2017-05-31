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

  this->UseRadiusInformation = 1;
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
    this->DistanceFunction->SetInput(this->Generators);
    this->DistanceFunction->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
    this->DistanceFunction->SetUseRadiusInformation(this->UseRadiusInformation);
    this->DistanceFunction->ControlEndPointsOff();
    this->DistanceFunction->UsePointNormalOn();

    // Get all the different ids
    vtkNew(vtkIdList, centerlineGroupIds);
    int i;
    for (int i=0; i<this->Generators->GetCellData()->GetArray(this->GroupIdsArrayName)->GetNumberOfTuples(); i++)
    {
      centerlineGroupIds->InsertUniqueId(static_cast<vtkIdType>(vtkMath::Round(this->Generators->GetCellData()->GetArray(this->GroupIdsArrayName)->GetComponent(i,0))));
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

      double pointNormal[3];
      this->CVTDataArray->GetTuple(i, pointNormal);
      this->DistanceFunction->SetPointNormal(pointNormal);

      int cellGenerator = 0;
      double minDist = VTK_SV_LARGE_DOUBLE;

      vtkNew(vtkIdList, groupId);
      groupId->SetNumberOfIds(1);
      for (int j=0; j<numGenerators; j++)
      {
        groupId->SetId(0, j);
        this->DistanceFunction->SetInputCellIds(groupId);

        double dist = this->DistanceFunction->EvaluateFunction(center);
        if (dist < minDist)
        {
          minDist = dist;
          cellGenerator = j;
        }
      }

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
  int numGenerators = this->Generators->GetNumberOfPoints();
  int currGenerator = this->PatchIdsArray->GetTuple1(evalId);
  newGenerator =  currGenerator;

  // Current minimum to beat is current generator
  double minDist = this->GetEdgeWeightedDistance(newGenerator, evalId);

  // Loop through neighboring patches
  for (int i=0; i<this->NumberOfNeighborPatches[evalId]; i++)
  {
    // Check to make sure not zero elements or the same genrator
    if (this->NeighborPatchesNumberOfElements[evalId][i] != 0)
    {
      if (currGenerator != this->NeighborPatchesIds[evalId][i])
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
  // Tube function
  vtkNew(vtkIdList, groupId);
  groupId->InsertNextId(generatorId);
  this->DistanceFunction->SetInputCellIds(groupId);
  double pointNormal[3];
  this->CVTDataArray->GetTuple(evalId, pointNormal);
  this->DistanceFunction->SetPointNormal(pointNormal);

  // Current generator
  int currGenerator = this->PatchIdsArray->GetTuple1(evalId);

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
  double edgeWeightedDist = this->DistanceFunction->EvaluateFunction(center);

  // Get the generator patch id
  int i;
  for (i=0; i<this->NumberOfNeighborPatches[evalId]; i++)
  {
    if (this->NeighborPatchesIds[evalId][i] == generatorId)
    {
      break;
    }
  }

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

  // Divide by the number of neighboring cells
  edgeWeighting /= this->NumberOfNeighbors[evalId];

  // Get the final edge distance to be normalized
  edgeWeightedDist += edgeWeighting;

  edgeWeightedDist = sqrt(edgeWeightedDist);

  return edgeWeightedDist;
  return SV_OK;
}
