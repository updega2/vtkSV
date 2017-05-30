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

#include "vtkSVEdgeWeightedCVT.h"

#include "vtkCellData.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVMathUtils.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVEdgeWeightedCVT);

// ----------------------
// Constructor
// ----------------------
vtkSVEdgeWeightedCVT::vtkSVEdgeWeightedCVT()
{
  this->NumberOfRings = 2;
  this->EdgeWeight = 1.0;
}

// ----------------------
// Destructor
// ----------------------
vtkSVEdgeWeightedCVT::~vtkSVEdgeWeightedCVT()
{
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVEdgeWeightedCVT::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Number of neighborhood rings: " << this->NumberOfRings << "\n";
  os << indent << "Edge weight: " << this->EdgeWeight << "\n";
}

// ----------------------
// InitializeConnectivity
// ----------------------
int vtkSVEdgeWeightedCVT::InitializeConnectivity()
{
  // Number of points and cells
  int numPoints = this->WorkPd->GetNumberOfPoints();
  int numCells =  this->WorkPd->GetNumberOfCells();

  this->WorkPd->BuildLinks();
  if (this->UseCellArray)
  {
    // Initialize point valences
    this->PointCellValenceNumber.resize(numPoints);
    this->PointCellValenceNumber.resize(numPoints);

    this->GetPointCellValence();

    // Initialize the cell neighbor vectors
    this->NumberOfNeighbors.resize(numCells);
    this->Neighbors.resize(numCells);

    // Set element as start
    for (int i=0; i<numCells; i++)
    {
      this->NumberOfNeighbors[i] = 1;
      this->Neighbors[i].push_back(i);
    }

    // Get element neighbor rings
    this->GetCellRingNeighbors();

    // Initialize the direct neighbor vectors
    this->NumberOfDirectNeighbors.resize(numCells);
    this->DirectNeighbors.resize(numCells);

    // Get cell edge neighbors
    this->GetCellDirectNeighbors();

    // Initialize cell group neighbors
    this->NumberOfNeighborGroups.resize(numCells);
    this->NeighborGroupsNumberOfElements.resize(numCells);
    this->NeighborGroupsIds.resize(numCells);

    // Get cell group neighbors
    this->GetCellGroupNeighbors();
  }
  else if (this->UsePointArray)
  {
    vtkErrorMacro("Not implemented yet");
    return SV_ERROR;
  }
  return SV_OK;
}

// ----------------------
// InitializeGenerators
// ----------------------
int vtkSVEdgeWeightedCVT::InitializeGenerators()
{
  if (this->UseCellArray)
  {
    int numCells = this->WorkPd->GetNumberOfCells();

    for (int i=0; i<numCells; i++)
    {
      int cellGenerator;
      this->GetClosestGenerator(i, cellGenerator);
      this->GroupIdsArray->SetTuple1(i, cellGenerator);
    }
  }
  return SV_OK;
}

// ----------------------
// GetClosestGenerator
// ----------------------
int vtkSVEdgeWeightedCVT::GetClosestGenerator(const int evalId, int &newGenerator)
{
  int numGenerators = this->Generators->GetNumberOfPoints();
  newGenerator = this->GroupIdsArray->GetTuple1(evalId);

  int minDist = VTK_SV_LARGE_DOUBLE;
  for (int i=0; i<numGenerators; i++)
  {
    if (this->UseGeneratorsArray)
    {
      int numComps = this->GeneratorsArray->GetNumberOfComponents();
    }
    else
    {
      double pt[3];
      this->Generators->GetPoint(i, pt);

      // NEIGHBORS FIRST TDODO!!!
      double cvtData[3];
      this->CVTDataArray->GetTuple(evalId, cvtData);
      double dist = vtkSVMathUtils::Distance(pt, cvtData);

      if (dist < minDist)
      {
        minDist = dist;
        newGenerator = i;
      }
    }
  }

  return SV_OK;
}

// ----------------------
// ComputeSurfaceMetric
// ----------------------
int vtkSVEdgeWeightedCVT::ComputeSurfaceMetric(double &evalMetric)
{
  return SV_OK;
}

// ----------------------
// UpdateConnectivity
// ----------------------
int vtkSVEdgeWeightedCVT::UpdateConnectivity()
{
  return SV_OK;
}

// ----------------------
// UpdateGenerators
// ----------------------
int vtkSVEdgeWeightedCVT::UpdateGenerators()
{
  return SV_OK;
}

// ----------------------
// GetPointCellValence
// ----------------------
int vtkSVEdgeWeightedCVT::GetPointCellValence()
{
  this->WorkPd->GetNumberOfPoints();
  // Number of points and cells
  int numPoints = this->WorkPd->GetNumberOfPoints();

  for (int i=0; i<numPoints; i++)
  {
    // get point cells
    vtkNew(vtkIdList, pointCells);

    this->WorkPd->GetPointCells(i, pointCells);

    // Set number of point cells
    this->PointCellValenceNumber[i] = pointCells->GetNumberOfIds();

    // Update point cell info
    for (int j=0; j<pointCells->GetNumberOfIds(); j++)
    {
      this->PointCellValence[i].push_back(pointCells->GetId(j));
    }
  }

  return SV_OK;
}

// ----------------------
// GetCellRingNeighbors
// ----------------------
int vtkSVEdgeWeightedCVT::GetCellRingNeighbors(int ringNumber)
{
  // Number of cells
  int numCells = this->WorkPd->GetNumberOfCells();

  for (int i=0; i<numCells; i++)
  {
    // temporary node vec
    std::vector<int> tmpNodes;
    int iSize = this->Neighbors[i].size();

    for (int j=0; j<iSize; j++)
    {
      // Get neighbor cell points
      int neiCellId = this->Neighbors[i][j];
      vtkIdType *pts, npts;
      this->WorkPd->GetCellPoints(neiCellId, npts, pts);

      // Loop around cell points
      for (int k=0; k<npts; k++)
      {
        int tmpNode = pts[k];
        int kSize   = tmpNodes.size();

        int kk = 0;
        for (kk=0; kk<kSize; kk++)
        {
          if (tmpNode == tmpNodes[kk])
          {
            break;
          }
        }
        if (kk == kSize)
        {
          tmpNodes.push_back(tmpNode);
        }
      }
    }

    // Now find neighbor elems
    iSize = tmpNodes.size();

    for (int j=0; j<iSize; j++)
    {
      int tmpNode = tmpNodes[j];

      for (int k=0; k<this->PointCellValenceNumber[tmpNode]; k++)
      {
        int tmpCell = this->PointCellValence[tmpNode][k];
        int kSize =   this->Neighbors[i].size();

        int kk=0;
        for (kk=0; kk<kSize; kk++)
        {
          if (tmpCell == this->Neighbors[i][kk])
          {
            break;
          }
        }
        if (kk == kSize)
        {
          this->Neighbors[i].push_back(tmpCell);
        }
      }
    }

    this->NumberOfNeighbors[i] = this->Neighbors[i].size();
  }

  if (ringNumber < this->NumberOfRings)
  {
    ringNumber++;
    this->GetCellRingNeighbors(ringNumber);
  }

  return SV_OK;
}

// ----------------------
// GetCellDirectNeighbors
// ----------------------
int vtkSVEdgeWeightedCVT::GetCellDirectNeighbors()
{

  int numCells = this->WorkPd->GetNumberOfCells();
  this->WorkPd->BuildLinks();

  // Loop through cells
  for (int i=0; i<numCells; i++)
  {
    // count number of edge neighbors
    int directNeiCount = 0;
    std::vector<int> neighborCells;

    // Get cell points
    vtkIdType *pts, npts;
    this->WorkPd->GetCellPoints(i, npts, pts);

    // Get cell edge neighbors
    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];

      // Get cell edge neighbors
      vtkNew(vtkIdList, cellEdgeNeighbors);
      directNeiCount += cellEdgeNeighbors->GetNumberOfIds();
      for (int k=0; k<cellEdgeNeighbors->GetNumberOfIds(); k++)
      {
        neighborCells.push_back(cellEdgeNeighbors->GetId(k));
      }
    }
    this->DirectNeighbors[i] = neighborCells;
    this->NumberOfDirectNeighbors[i] = directNeiCount;
  }

  return SV_OK;
}

// ----------------------
// GetCellGroupNeighbors
// ----------------------
int vtkSVEdgeWeightedCVT::GetCellGroupNeighbors()
{

  int numCells = this->WorkPd->GetNumberOfCells();
  this->WorkPd->BuildLinks();

  // Loop through cells
  for (int i=0; i<numCells; i++)
  {
    for (int j=0; j<this->NumberOfNeighbors[i]; j++)
    {
      int cellNeighbor = this->Neighbors[i][j];
      if (cellNeighbor != i)
      {
        int cellNeighborGroup = this->GroupIdsArray->GetTuple1(cellNeighbor);
        int neighborLoc;
        this->AddCellGroupNeighbor(i, cellNeighborGroup, neighborLoc);
      }
    }

  }

  // Check to make sure correct
  int error = 0;
  int numNeighbors = 0;

  for (int i=0; i<numCells; i++)
  {
    numNeighbors = 0;

    for (int j=0; j<this->NumberOfNeighborGroups[i]; j++)
    {
      numNeighbors += this->NeighborGroupsNumberOfElements[i][j];
    }

    if (numNeighbors != this->NumberOfNeighbors[i])
    {
      error++;
    }
  }

  if (error != 0)
  {
    vtkErrorMacro("Cell ring neighbors not correct");
    return SV_ERROR;
  }

  return SV_OK;
}

// ----------------------
// AddCellGroupNeighbor
// ----------------------
int vtkSVEdgeWeightedCVT::AddCellGroupNeighbor(const int cellId, const int cellNeighborGroup, int &neighborLoc)
{
  int counted = 0;
  int numNeighborGroups;

  for (int i=0; i<this->NumberOfNeighborGroups[cellId]; i++)
  {
    if (this->NeighborGroupsIds[cellId][i] == cellNeighborGroup)
    {
      counted = 1;
      this->NeighborGroupsNumberOfElements[cellId][i]++;
      neighborLoc = i;

      break;
    }
  }

  if (counted == 0)
  {
    this->NumberOfNeighborGroups[cellId]++;
    numNeighborGroups = this->NumberOfNeighborGroups[cellId];

  }
  this->NeighborGroupsIds[cellId][numNeighborGroups -1] = cellNeighborGroup;
  this->NeighborGroupsNumberOfElements[cellId][numNeighborGroups - 1]++;
  neighborLoc = numNeighborGroups - 1;

  return counted;
}
