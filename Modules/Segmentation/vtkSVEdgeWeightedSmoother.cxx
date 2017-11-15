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

#include "vtkSVEdgeWeightedSmoother.h"

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
vtkStandardNewMacro(vtkSVEdgeWeightedSmoother);

// ----------------------
// Constructor
// ----------------------
vtkSVEdgeWeightedSmoother::vtkSVEdgeWeightedSmoother()
{
  this->UseCurvatureWeight = 1;
}

// ----------------------
// Destructor
// ----------------------
vtkSVEdgeWeightedSmoother::~vtkSVEdgeWeightedSmoother()
{
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVEdgeWeightedSmoother::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// ----------------------
// InitializeGenerators
// ----------------------
int vtkSVEdgeWeightedSmoother::InitializeGenerators()
{
  return SV_OK;
}

// ----------------------
// UpdateGenerators
// ----------------------
int vtkSVEdgeWeightedSmoother::UpdateGenerators()
{
  return SV_OK;
}


// ----------------------
// GetClosestGenerator
// ----------------------
int vtkSVEdgeWeightedSmoother::GetClosestGenerator(const int evalId, int &newGenerator)
{
  // Get current generator
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
double vtkSVEdgeWeightedSmoother::GetEdgeWeightedDistance(const int generatorId, const int evalId)
{
  // Current generator
  int currGenerator = this->PatchIdsArray->GetTuple1(evalId);

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

  double totalWeight = 0.0;
  int i;
  for (i=0; i<this->NumberOfNeighbors[evalId]; i++)
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

  // Get the generator patch id
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
  //}

  // Divide by the number of neighboring cells
  edgeWeighting /= this->NumberOfNeighbors[evalId];

  // Get the final edge distance to be normalized
  edgeWeightedDist += edgeWeighting;

  edgeWeightedDist = sqrt(edgeWeightedDist);

  return edgeWeightedDist;
  return SV_OK;
}

