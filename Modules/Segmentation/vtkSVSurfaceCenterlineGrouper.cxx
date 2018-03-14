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

#include "vtkSVSurfaceCenterlineGrouper.h"

#include "vtkSVCenterlineBranchSplitter.h"
#include "vtkSVCenterlinesEdgeWeightedCVT.h"
#include "vtkSVFindGeodesicPath.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVMathUtils.h"
#include "vtkSVPolycubeGenerator.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDoubleArray.h"
#include "vtkExecutive.h"
#include "vtkErrorCode.h"
#include "vtkFeatureEdges.h"
#include "vtkIdFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPolyLine.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkThreshold.h"
#include "vtkTriangleFilter.h"
#include "vtkUnstructuredGrid.h"

#include <algorithm>

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVSurfaceCenterlineGrouper);

// ----------------------
// Constructor
// ----------------------
vtkSVSurfaceCenterlineGrouper::vtkSVSurfaceCenterlineGrouper()
{
  this->WorkPd = vtkPolyData::New();
  this->MergedCenterlines = vtkPolyData::New();
  this->PolycubePd = NULL;

  this->CenterlineGroupIdsArrayName = NULL;
  this->CenterlineRadiusArrayName = NULL;
  this->CenterlineIdsArrayName = NULL;
  this->GroupIdsArrayName = NULL;
  this->BlankingArrayName = NULL;
  this->TractIdsArrayName = NULL;

  this->UseRadiusInformation = 1;
  this->EnforcePolycubeBoundaries = 0;
  this->GroupSurface = 1;
}

// ----------------------
// Destructor
// ----------------------
vtkSVSurfaceCenterlineGrouper::~vtkSVSurfaceCenterlineGrouper()
{
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->MergedCenterlines != NULL)
  {
    this->MergedCenterlines->Delete();
    this->MergedCenterlines = NULL;
  }
  if (this->PolycubePd != NULL)
  {
    this->PolycubePd->Delete();
    this->PolycubePd = NULL;
  }

  if (this->CenterlineGroupIdsArrayName != NULL)
  {
    delete [] this->CenterlineGroupIdsArrayName;
    this->CenterlineGroupIdsArrayName = NULL;
  }

  if (this->CenterlineRadiusArrayName != NULL)
  {
    delete [] this->CenterlineRadiusArrayName;
    this->CenterlineRadiusArrayName = NULL;
  }

  if (this->CenterlineIdsArrayName != NULL)
  {
    delete [] this->CenterlineIdsArrayName;
    this->CenterlineIdsArrayName = NULL;
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

  if (this->TractIdsArrayName != NULL)
  {
    delete [] this->TractIdsArrayName;
    this->TractIdsArrayName = NULL;
  }

}

// ----------------------
// RequestData
// ----------------------
int vtkSVSurfaceCenterlineGrouper::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  this->WorkPd->DeepCopy(input);

  // Prep work for filter
  if (this->PrepFilter() != SV_OK)
  {
    vtkErrorMacro("Prep of filter failed");
    output->DeepCopy(input);
    this->SetErrorCode(vtkErrorCode::UserError + 1);
    return SV_ERROR;
  }

  // Run the filter
  if (this->RunFilter() != SV_OK)
  {
    vtkErrorMacro("Filter failed");
    output->DeepCopy(this->WorkPd);
    this->SetErrorCode(vtkErrorCode::UserError + 2);
    return SV_ERROR;
  }

  output->DeepCopy(this->WorkPd);

  return SV_OK;
}

// ----------------------
// PrepFilter
// ----------------------
int vtkSVSurfaceCenterlineGrouper::PrepFilter()
{
  if (!this->MergedCenterlines)
  {
    vtkErrorMacro(<< "Centerlines not set.");
    return SV_ERROR;
  }

  if (!this->CenterlineGroupIdsArrayName)
  {
    vtkDebugMacro("Centerline GroupIds Array Name not given, setting to GroupIds");
    this->CenterlineGroupIdsArrayName = new char[strlen("GroupIds") + 1];
    strcpy(this->CenterlineGroupIdsArrayName, "GroupIds");
  }

  if (!this->GroupIdsArrayName)
  {
    vtkDebugMacro("GroupIds Array Name not given, setting to GroupIds");
    this->GroupIdsArrayName = new char[strlen("GroupIds") + 1];
    strcpy(this->GroupIdsArrayName, "GroupIds");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->MergedCenterlines, 1, this->CenterlineGroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "CenterlineGroupIdsArray with name specified does not exist on centerlines");
    return SV_OK;
  }

  if (!this->BlankingArrayName)
  {
    vtkDebugMacro("Blanking Array Name not given, setting to Blanking");
    this->BlankingArrayName = new char[strlen("Blanking") + 1];
    strcpy(this->BlankingArrayName, "Blanking");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->MergedCenterlines, 1, this->BlankingArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "BlankingArrayName with name specified does not exist on centerlines");
    return SV_ERROR;
  }

  if (!this->CenterlineRadiusArrayName)
  {
    vtkDebugMacro("Centerline radius Array Name not given, setting to MaximumInscribedSphereRadius");
    this->CenterlineRadiusArrayName = new char[strlen("MaximumInscribedSphereRadius") + 1];
    strcpy(this->CenterlineRadiusArrayName, "MaximumInscribedSphereRadius");
  }

  if (!this->MergedCenterlines->GetPointData()->GetArray(this->CenterlineRadiusArrayName))
  {
    vtkErrorMacro(<< "CenterlineRadiusArray with name specified does not exist on centerlines");
    return SV_ERROR;
  }

  if (!this->CenterlineIdsArrayName)
  {
    vtkDebugMacro("CenterlineIds Array Name not given, setting to CenterlineIds");
    this->CenterlineIdsArrayName = new char[strlen("CenterlineIds") + 1];
    strcpy(this->CenterlineIdsArrayName, "CenterlineIds");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->MergedCenterlines, 1, this->CenterlineIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "CenterlineIdsArray with name specified does not exist on centerlines");
    return SV_OK;
  }

  if (!this->TractIdsArrayName)
  {
    vtkDebugMacro("TractIds Array Name not given, setting to TractIds");
    this->TractIdsArrayName = new char[strlen("TractIds") + 1];
    strcpy(this->TractIdsArrayName, "TractIds");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->MergedCenterlines, 1, this->TractIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "TractIdsArray with name specified does not exist on centerlines");
    return SV_OK;
  }

  if (this->EnforcePolycubeBoundaries)
  {
    vtkDebugMacro("Need to enforce polycube");
    if (this->PolycubePd == NULL)
    {
      vtkDebugMacro("Polycube not provided, generating from centerlines");

      vtkNew(vtkSVPolycubeGenerator, polycuber);
      polycuber->SetInputData(this->MergedCenterlines);
      polycuber->SetCenterlineGroupIdsArrayName(this->CenterlineGroupIdsArrayName);
      polycuber->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
      polycuber->Update();

      this->PolycubePd = vtkPolyData::New();
      this->PolycubePd->DeepCopy(polycuber->GetOutput());
      vtkDebugMacro("Polycube created");
    }

    if (this->PolycubePd->GetNumberOfCells() == 0)
    {
      vtkErrorMacro("Polycube is empty");
      return SV_ERROR;
    }
  }

  if (!this->GroupSurface)
  {
    if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 1, this->GroupIdsArrayName) != SV_OK)
    {
      vtkErrorMacro(<< "GroupIdsArray with name specified does not exist on surface");
      return SV_OK;
    }
  }

  //TODO: Check to make sure centerlines and polycube match if polycube given

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVSurfaceCenterlineGrouper::RunFilter()
{
  // Generate normals just in case they don't exist
  vtkNew(vtkPolyDataNormals, normaler);
  normaler->SetInputData(this->WorkPd);
  normaler->ComputePointNormalsOff();
  normaler->ComputeCellNormalsOn();
  normaler->SplittingOff();
  normaler->Update();

  this->WorkPd->DeepCopy(normaler->GetOutput());
  this->WorkPd->BuildLinks();
  vtkDataArray *normalsArray =
    this->WorkPd->GetCellData()->GetArray("Normals");

  if (this->GroupSurface)
  {
    int stopCellNumber = ceil(this->WorkPd->GetNumberOfCells()*0.0001);

    vtkNew(vtkSVCenterlinesEdgeWeightedCVT, CVT);
    CVT->SetInputData(this->WorkPd);
    CVT->SetGenerators(this->MergedCenterlines);
    CVT->SetNumberOfRings(2);
    CVT->SetThreshold(stopCellNumber);
    CVT->SetUseCurvatureWeight(0);
    CVT->SetPatchIdsArrayName(this->GroupIdsArrayName);
    CVT->SetCVTDataArrayName("Normals");
    CVT->SetGroupIdsArrayName(this->GroupIdsArrayName);
    CVT->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
    CVT->SetUseRadiusInformation(this->UseRadiusInformation);
    CVT->SetUsePointNormal(1);
    CVT->SetMaximumNumberOfIterations(0);
    CVT->Update();

    this->WorkPd->DeepCopy(CVT->GetOutput());

  }

  if (this->CheckGroups(this->WorkPd) != SV_OK)
  {
    vtkErrorMacro("Error in correcting groups");
    return SV_ERROR;
  }

  if (this->CorrectCellBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro("Could not correcto boundaries of surface");
    return SV_ERROR;
  }

  if (this->CheckGroups2() != SV_OK)
  {
    vtkErrorMacro("Error in correcting groups");
    return SV_ERROR;
  }

  if (this->EnforcePolycubeBoundaries)
  {
    if (this->FixGroupsWithPolycube() != SV_OK)
    {
      vtkErrorMacro("Error in correcting groups");
      return SV_ERROR;
    }

    if (this->CorrectCellBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
    {
      vtkErrorMacro("Could not correcto boundaries of surface");
      return SV_ERROR;
    }

  }

  if (this->SmoothBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro("Could not smootho boundaries of surface");
    return SV_ERROR;
  }

  std::vector<Region> groupRegions;
  if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, groupRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get group regions");
    return SV_ERROR;
  }
  if (this->CurveFitBoundaries(this->WorkPd, this->GroupIdsArrayName, groupRegions) != SV_OK)
  {
    vtkErrorMacro("Could not curve fit boundaries of surface");
    return SV_ERROR;
  }

  if (this->EnforcePolycubeBoundaries)
  {
    if (this->MatchSurfaceToPolycube() != SV_OK)
    {
      vtkErrorMacro("Couldn't fix stuff\n");
      return SV_ERROR;
    }

    if (this->CheckSlicePoints() != SV_OK)
    {
      vtkErrorMacro("Error when checking slice points\n");
      return SV_ERROR;
    }
  }

  // Get new normals
  normaler->SetInputData(this->WorkPd);
  normaler->ComputePointNormalsOff();
  normaler->ComputeCellNormalsOn();
  normaler->SplittingOff();
  normaler->Update();
  this->WorkPd->DeepCopy(normaler->GetOutput());
  this->WorkPd->BuildLinks();

  return SV_OK;
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVSurfaceCenterlineGrouper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Use radius information: " << this->UseRadiusInformation << "\n";
  if (this->CenterlineGroupIdsArrayName != NULL)
    os << indent << "Centerline group ids name: " << this->CenterlineGroupIdsArrayName << "\n";
  if (this->CenterlineRadiusArrayName != NULL)
    os << indent << "Centerline radius array name: " << this->CenterlineRadiusArrayName << "\n";
  if (this->CenterlineIdsArrayName != NULL)
    os << indent << "Centerline ids array name: " << this->CenterlineIdsArrayName << "\n";
  if (this->GroupIdsArrayName != NULL)
    os << indent << "Group ids array name: " << this->GroupIdsArrayName << "\n";
  if (this->BlankingArrayName != NULL)
    os << indent << "Blanking array name: " << this->BlankingArrayName << "\n";
  if (this->TractIdsArrayName != NULL)
    os << indent << "Tract ids array name: " << this->TractIdsArrayName << "\n";
}

// ----------------------
// CorrectCellBoundaries
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CorrectCellBoundaries(vtkPolyData *pd, std::string cellArrayName )
{
  // Get current cell ids
  vtkDataArray *cellIds = pd->GetCellData()->GetArray(cellArrayName.c_str());

  // Num cells
  pd->BuildLinks();
  int numCells = pd->GetNumberOfCells();

  // Set up array to keep track of temp cell ids, will be different than
  // cellIds if disconnected regions of same value
  vtkNew(vtkIntArray, tmpIds);
  tmpIds->SetNumberOfTuples(numCells);
  tmpIds->FillComponent(0, -1);

  // Set count var
  int regionCount =0;

  // Loop through cells
  for (int i=0; i<numCells; i++)
  {
    // If not set yet
    if (tmpIds->GetTuple1(i) == -1)
    {
      tmpIds->SetTuple1(i, regionCount);

      // Count cells in connected region
      int count = 1;
      vtkNew(vtkIdList, queue);
      queue->InsertId(0, i);

      // Loop through updating count
      for (int j=0; j<count; j++)
      {
        // Get Cell points
        vtkIdType npts, *pts;
        pd->GetCellPoints(queue->GetId(j), npts, pts);

        // Loop through cell points
        for (int k=0; k<npts; k++)
        {
          int ptId0 = pts[k];
          int ptId1 = pts[(k+1)%npts];

          // Get cell edge neighbors
          vtkNew(vtkIdList, cellEdgeNeighbors);
          pd->GetCellEdgeNeighbors(queue->GetId(j), ptId0, ptId1, cellEdgeNeighbors);

          // Check val of cell edge neighbors
          for (int l=0; l<cellEdgeNeighbors->GetNumberOfIds(); l++)
          {
            int cellEdgeNeighbor = cellEdgeNeighbors->GetId(l);
            if (tmpIds->GetTuple1(cellEdgeNeighbor) == -1 &&
                cellIds->GetTuple1(i) == cellIds->GetTuple1(cellEdgeNeighbor))
            {
              // Update cell val, count
              tmpIds->SetTuple1(cellEdgeNeighbor, regionCount);
              queue->InsertNextId(cellEdgeNeighbor);
              count++;
            }
          }
        }
      }
      regionCount++;
    }
  }

  int allGood = 0;
  int iter = 0;
  int maxIters = 100;

  while(!allGood && iter < maxIters)
  {
    allGood = 1;
    // Loop through cells again
    for (int i=0; i<numCells; i++)
    {

      // get direct neighbor value count
      vtkNew(vtkIdList, neiCellIds);
      vtkNew(vtkIdList, neiTmpIds);

      // Get cell points
      vtkIdType npts, *pts;
      pd->GetCellPoints(i, npts, pts);

      // Loop through cell points
      for (int j=0; j<npts; j++)
      {
        int ptId0 = pts[j];
        int ptId1 = pts[(j+1)%npts];

        // Get cell edge neighbors
        vtkNew(vtkIdList, cellEdgeNeighbors);
        pd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellEdgeNeighbors);

        // loop through neighbors
        for (int k=0; k<cellEdgeNeighbors->GetNumberOfIds(); k++)
        {
          int cellEdgeNeighbor = cellEdgeNeighbors->GetId(k);

          // Check to see if equal to region val
          if (tmpIds->GetTuple1(cellEdgeNeighbor) != tmpIds->GetTuple1(i))
          {
            neiCellIds->InsertNextId(cellIds->GetTuple1(cellEdgeNeighbor));
            neiTmpIds->InsertNextId(tmpIds->GetTuple1(cellEdgeNeighbor));
          }
        }
      }

      // If we found a cell surrounded by cells of another val, we can update
      vtkSortDataArray::Sort(neiTmpIds);
      int neiSize = neiTmpIds->GetNumberOfIds();
      if (neiSize == 2)
      {
        if (neiTmpIds->GetId(0) == neiTmpIds->GetId(1))
        {
          allGood = 0;
          int maxVal, maxCount;
          vtkSVSurfaceCenterlineGrouper::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
      }
      else if (neiSize >= 3)
      {
        if ((neiTmpIds->GetId(0) == neiTmpIds->GetId(1) ||
             neiTmpIds->GetId(1) == neiTmpIds->GetId(2)))
        {
          allGood = 0;
          int maxVal, maxCount;
          vtkSVSurfaceCenterlineGrouper::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
      }
    }
    iter++;
  }

  return 1;
}

// ----------------------
// CorrectSpecificCellBoundaries
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CorrectSpecificCellBoundaries(vtkPolyData *pd, std::string cellArrayName, vtkIdList *targetRegions)
{
  // Get current cell ids
  vtkDataArray *cellIds = pd->GetCellData()->GetArray(cellArrayName.c_str());

  // Num cells
  pd->BuildLinks();
  int numCells = pd->GetNumberOfCells();

  // Set up array to keep track of temp cell ids, will be different than
  // cellIds if disconnected regions of same value
  vtkNew(vtkIntArray, tmpIds);
  tmpIds->SetNumberOfTuples(numCells);
  tmpIds->FillComponent(0, -1);

  // Set count var
  int regionCount =0;

  // Loop through cells
  for (int i=0; i<numCells; i++)
  {
    // If not set yet
    if (tmpIds->GetTuple1(i) == -1)
    {
      tmpIds->SetTuple1(i, regionCount);

      // Count cells in connected region
      int count = 1;
      vtkNew(vtkIdList, queue);
      queue->InsertId(0, i);

      // Loop through updating count
      for (int j=0; j<count; j++)
      {
        // Get Cell points
        vtkIdType npts, *pts;
        pd->GetCellPoints(queue->GetId(j), npts, pts);

        // Loop through cell points
        for (int k=0; k<npts; k++)
        {
          int ptId0 = pts[k];
          int ptId1 = pts[(k+1)%npts];

          // Get cell edge neighbors
          vtkNew(vtkIdList, cellEdgeNeighbors);
          pd->GetCellEdgeNeighbors(queue->GetId(j), ptId0, ptId1, cellEdgeNeighbors);

          // Check val of cell edge neighbors
          for (int l=0; l<cellEdgeNeighbors->GetNumberOfIds(); l++)
          {
            int cellEdgeNeighbor = cellEdgeNeighbors->GetId(l);
            if (tmpIds->GetTuple1(cellEdgeNeighbor) == -1 &&
                cellIds->GetTuple1(i) == cellIds->GetTuple1(cellEdgeNeighbor))
            {
              // Update cell val, count
              tmpIds->SetTuple1(cellEdgeNeighbor, regionCount);
              queue->InsertNextId(cellEdgeNeighbor);
              count++;
            }
          }
        }
      }
      regionCount++;
    }
  }

  int allGood = 0;
  int iter = 0;
  int maxIters = 100;

  while(!allGood && iter<maxIters)
  {
    allGood = 1;
    // Loop through cells again
    for (int i=0; i<numCells; i++)
    {

      // get direct neighbor value count
      vtkNew(vtkIdList, neiCellIds);
      vtkNew(vtkIdList, neiTmpIds);

      // Get cell points
      vtkIdType npts, *pts;
      pd->GetCellPoints(i, npts, pts);

      // Loop through cell points
      for (int j=0; j<npts; j++)
      {
        int ptId0 = pts[j];
        int ptId1 = pts[(j+1)%npts];

        // Get cell edge neighbors
        vtkNew(vtkIdList, cellEdgeNeighbors);
        pd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellEdgeNeighbors);

        // loop through neighbors
        for (int k=0; k<cellEdgeNeighbors->GetNumberOfIds(); k++)
        {
          int cellEdgeNeighbor = cellEdgeNeighbors->GetId(k);

          // Check to see if equal to region val
          if (tmpIds->GetTuple1(cellEdgeNeighbor) != tmpIds->GetTuple1(i))
          {
            int cellValue = cellIds->GetTuple1(cellEdgeNeighbor);
            if (targetRegions->IsId(cellValue) != -1)
            {
              neiCellIds->InsertNextId(cellValue);
              neiTmpIds->InsertNextId(tmpIds->GetTuple1(cellEdgeNeighbor));
            }
          }
        }
      }

      // If we found a cell surrounded by cells of another val, we can update
      vtkSortDataArray::Sort(neiTmpIds);
      int neiSize = neiTmpIds->GetNumberOfIds();
      if (neiSize == 2)
      {
        if (neiTmpIds->GetId(0) == neiTmpIds->GetId(1))
        {
          allGood = 0;
          int maxVal, maxCount;
          vtkSVSurfaceCenterlineGrouper::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
      }
      else if (neiSize >= 3)
      {
        if ((neiTmpIds->GetId(0) == neiTmpIds->GetId(1) ||
             neiTmpIds->GetId(1) == neiTmpIds->GetId(2)))
        {
          allGood = 0;
          int maxVal, maxCount;
          vtkSVSurfaceCenterlineGrouper::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
      }
    }
    iter++;
  }

  return 1;
}

// ----------------------
// GetMostOccuringVal
// ----------------------
void vtkSVSurfaceCenterlineGrouper::GetMostOccuringVal(vtkIdList *idList, int &output,
                                             int &max_count)
{
  int numIds = idList->GetNumberOfIds();

  max_count = 0;
  int max_val = idList->GetId(0);
  for (int i=0; i<numIds; i++)
  {
    int count = 1;
    for (int j=0; j<numIds; j++)
    {
      if (idList->GetId(i) == idList->GetId(j))
        count++;
    }
    if (count > max_count)
    {
      max_count = count;
      max_val = idList->GetId(i);
    }
  }

  output = max_val;
}

// ----------------------
// SmoothBoundaries
// ----------------------
int vtkSVSurfaceCenterlineGrouper::SmoothBoundaries(vtkPolyData *pd, std::string arrayName)
{
  int numPoints = pd->GetNumberOfPoints();
  std::vector<int> cornerPoints;
  std::vector<int> isCornerPoint(numPoints);
  std::vector<int> isBoundaryPoint(numPoints);
  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, pointCellsValues);
    vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, i, pointCellsValues);
    if (pointCellsValues->GetNumberOfIds() >= 3)
    {
      cornerPoints.push_back(i);
      isCornerPoint[i] = 1;
    }
    else
      isCornerPoint[i] = 0;

    if (pointCellsValues->GetNumberOfIds() == 2)
      isBoundaryPoint[i] = 1;
    else
      isBoundaryPoint[i] = 0;
  }


  for (int i=0; i<numPoints; i++)
  {
    if (isBoundaryPoint[i])
    {
      vtkNew(vtkIdList, pointCellsValues);
      vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName,
                                             i, pointCellsValues);

      // boundary edge
      if (pointCellsValues->GetNumberOfIds() == 2)
      {
        vtkNew(vtkIdList, pointCells);
        pd->GetPointCells(i, pointCells);

        int count[2]; count[0] = 0; count[1] = 0;
        int cellIds[2][2];
        for (int j=0; j<pointCells->GetNumberOfIds(); j++)
        {
          for (int k=0; k<2; k++)
          {
            if (pd->GetCellData()->GetArray(
              arrayName.c_str())->GetTuple1(pointCells->GetId(j)) == pointCellsValues->GetId(k))
            {
              if (count[k] < 2)
                cellIds[k][count[k]] = pointCells->GetId(j);
              count[k]++;
            }
          }
        }

        if (count[0] == 2 || count[1] == 2)
        {
          vtkNew(vtkIdList, uniquePoints);
          vtkIdType npts, *pts;
          if (count[0] == 2 && count[1] == 2)
          {
            for (int j=0; j<2; j++)
            {
              for (int k=0; k<2; k++)
              {
                pd->GetCellPoints(cellIds[j][k], npts, pts);
                for (int p=0; p<npts; p++)
                  uniquePoints->InsertUniqueId(pts[p]);
              }
            }
          }
          else
          {
            if (count[0] == 2)
            {
              vtkNew(vtkIdList, avgPoints);
              double check0[3]; check0[0] = 0.0; check0[1] = 0.0; check0[2] = 0.0;
              for (int j=0; j<2; j++)
              {
                pd->GetCellPoints(cellIds[0][j], npts, pts);
                for (int p=0; p<npts; p++)
                {
                  int isId = avgPoints->IsId(pts[p]);
                  if (isId == -1)
                  {
                    avgPoints->InsertNextId(pts[p]);
                    double pt[3];
                    pd->GetPoint(pts[p], pt);
                    for (int r=0; r<3; r++)
                      check0[r] += pt[r];
                  }
                }
              }
              for (int j=0; j<3; j++)
                check0[j] = (1./avgPoints->GetNumberOfIds())*check0[j];

              vtkNew(vtkIdList, halfPoints);
              double check1[3]; check1[0] = 0.0; check1[1] = 0.0; check1[2] = 0.0;
              pd->GetCellPoints(cellIds[0][0], npts, pts);
              for (int j=0; j<npts; j++)
              {
                int ptId0 = pts[j];
                int ptId1 = pts[(j+1)%npts];
                vtkNew(vtkIdList, neighborCell);
                pd->GetCellEdgeNeighbors(cellIds[0][0], ptId0, ptId1, neighborCell);
                if (neighborCell->GetNumberOfIds() > 0)
                {
                  if (neighborCell->GetId(0) == cellIds[0][1])
                  {
                    halfPoints->InsertNextId(ptId0);
                    halfPoints->InsertNextId(ptId1);

                    double pt[3];
                    pd->GetPoint(ptId0, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    pd->GetPoint(ptId1, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    for (int r=0; r<3; r++)
                      check1[r] = (1./2)*check1[r];
                  }
                }
              }

              double startPt[3];
              pd->GetPoint(i, startPt);
              double dist0 = vtkSVMathUtils::Distance(check0, startPt);
              double dist1 = vtkSVMathUtils::Distance(check1, startPt);

              if (dist0 > dist1)
              {
                for (int r=0; r<halfPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(halfPoints->GetId(r));
              }
              else
              {
                for (int r=0; r<avgPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(avgPoints->GetId(r));
              }
            }
            else if (count[1] == 2)
            {
              vtkNew(vtkIdList, avgPoints);
              double check0[3]; check0[0] = 0.0; check0[1] = 0.0; check0[2] = 0.0;
              for (int j=0; j<2; j++)
              {
                pd->GetCellPoints(cellIds[1][j], npts, pts);
                for (int p=0; p<npts; p++)
                {
                  int isId = avgPoints->IsId(pts[p]);
                  if (isId == -1)
                  {
                    avgPoints->InsertNextId(pts[p]);
                    double pt[3];
                    pd->GetPoint(pts[p], pt);
                    for (int r=0; r<3; r++)
                      check0[r] += pt[r];
                  }
                }
              }
              for (int j=0; j<3; j++)
                check0[j] = (1./avgPoints->GetNumberOfIds())*check0[j];

              vtkNew(vtkIdList, halfPoints);
              double check1[3]; check1[0] = 0.0; check1[1] = 0.0; check1[2] = 0.0;
              pd->GetCellPoints(cellIds[1][0], npts, pts);
              for (int j=0; j<npts; j++)
              {
                int ptId0 = pts[j];
                int ptId1 = pts[(j+1)%npts];
                vtkNew(vtkIdList, neighborCell);
                pd->GetCellEdgeNeighbors(cellIds[1][0], ptId0, ptId1, neighborCell);
                if (neighborCell->GetNumberOfIds() > 0)
                {
                  if (neighborCell->GetId(0) == cellIds[1][1])
                  {
                    halfPoints->InsertNextId(ptId0);
                    halfPoints->InsertNextId(ptId1);

                    double pt[3];
                    pd->GetPoint(ptId0, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    pd->GetPoint(ptId1, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    for (int r=0; r<3; r++)
                      check1[r] = (1./2)*check1[r];
                  }
                }
              }

              double startPt[3];
              pd->GetPoint(i, startPt);
              double dist0 = vtkSVMathUtils::Distance(check0, startPt);
              double dist1 = vtkSVMathUtils::Distance(check1, startPt);

              if (dist0 > dist1)
              {
                for (int r=0; r<halfPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(halfPoints->GetId(r));
              }
              else
              {
                for (int r=0; r<avgPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(avgPoints->GetId(r));
              }
            }
          }
          int numIds = uniquePoints->GetNumberOfIds();
          double center[3];
          for (int j=0; j<3; j++)
            center[j] = 0.0;
          for (int k=0; k<numIds; k++)
          {
            double pt[3];
            pd->GetPoint(uniquePoints->GetId(k), pt);
            for (int j=0; j<3; j++)
              center[j] += pt[j];
          }
          for (int j=0; j<3; j++)
            center[j] = (1./numIds)*center[j];

          pd->GetPoints()->SetPoint(i, center);
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// SmoothSpecificBoundaries
// ----------------------
int vtkSVSurfaceCenterlineGrouper::SmoothSpecificBoundaries(vtkPolyData *pd, std::string arrayName, vtkIdList *targetRegions)
{
  int numPoints = pd->GetNumberOfPoints();
  std::vector<int> cornerPoints;
  std::vector<int> isCornerPoint(numPoints);
  std::vector<int> isBoundaryPoint(numPoints);
  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, pointCellsValues);
    vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, i, pointCellsValues);
    if (pointCellsValues->GetNumberOfIds() >= 3)
    {
      cornerPoints.push_back(i);
      isCornerPoint[i] = 1;
    }
    else
      isCornerPoint[i] = 0;

    if (pointCellsValues->GetNumberOfIds() == 2)
    {
      if (targetRegions->IsId(pointCellsValues->GetId(0)) != -1 &&
          targetRegions->IsId(pointCellsValues->GetId(1)) != -1)
        isBoundaryPoint[i] = 1;
      else
        isBoundaryPoint[i] = 0;
    }
    else
      isBoundaryPoint[i] = 0;
  }

  for (int i=0; i<numPoints; i++)
  {
    if (isBoundaryPoint[i])
    {
      vtkNew(vtkIdList, pointCellsValues);
      vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName,
                                             i, pointCellsValues);

      // boundary edge
      if (pointCellsValues->GetNumberOfIds() == 2)
      {
        vtkNew(vtkIdList, pointCells);
        pd->GetPointCells(i, pointCells);

        int count[2]; count[0] = 0; count[1] = 0;
        int cellIds[2][2];
        for (int j=0; j<pointCells->GetNumberOfIds(); j++)
        {
          for (int k=0; k<2; k++)
          {
            if (pd->GetCellData()->GetArray(
              arrayName.c_str())->GetTuple1(pointCells->GetId(j)) == pointCellsValues->GetId(k))
            {
              if (count[k] < 2)
                cellIds[k][count[k]] = pointCells->GetId(j);
              count[k]++;
            }
          }
        }

        if (count[0] == 2 || count[1] == 2)
        {
          vtkNew(vtkIdList, uniquePoints);
          vtkIdType npts, *pts;
          if (count[0] == 2 && count[1] == 2)
          {
            for (int j=0; j<2; j++)
            {
              for (int k=0; k<2; k++)
              {
                pd->GetCellPoints(cellIds[j][k], npts, pts);
                for (int p=0; p<npts; p++)
                  uniquePoints->InsertUniqueId(pts[p]);
              }
            }
          }
          else
          {
            if (count[0] == 2)
            {
              vtkNew(vtkIdList, avgPoints);
              double check0[3]; check0[0] = 0.0; check0[1] = 0.0; check0[2] = 0.0;
              for (int j=0; j<2; j++)
              {
                pd->GetCellPoints(cellIds[0][j], npts, pts);
                for (int p=0; p<npts; p++)
                {
                  int isId = avgPoints->IsId(pts[p]);
                  if (isId == -1)
                  {
                    avgPoints->InsertNextId(pts[p]);
                    double pt[3];
                    pd->GetPoint(pts[p], pt);
                    for (int r=0; r<3; r++)
                      check0[r] += pt[r];
                  }
                }
              }
              for (int j=0; j<3; j++)
                check0[j] = (1./avgPoints->GetNumberOfIds())*check0[j];

              vtkNew(vtkIdList, halfPoints);
              double check1[3]; check1[0] = 0.0; check1[1] = 0.0; check1[2] = 0.0;
              pd->GetCellPoints(cellIds[0][0], npts, pts);
              for (int j=0; j<npts; j++)
              {
                int ptId0 = pts[j];
                int ptId1 = pts[(j+1)%npts];
                vtkNew(vtkIdList, neighborCell);
                pd->GetCellEdgeNeighbors(cellIds[0][0], ptId0, ptId1, neighborCell);
                if (neighborCell->GetNumberOfIds() > 0)
                {
                  if (neighborCell->GetId(0) == cellIds[0][1])
                  {
                    halfPoints->InsertNextId(ptId0);
                    halfPoints->InsertNextId(ptId1);

                    double pt[3];
                    pd->GetPoint(ptId0, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    pd->GetPoint(ptId1, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    for (int r=0; r<3; r++)
                      check1[r] = (1./2)*check1[r];
                  }
                }
              }

              double startPt[3];
              pd->GetPoint(i, startPt);
              double dist0 = vtkSVMathUtils::Distance(check0, startPt);
              double dist1 = vtkSVMathUtils::Distance(check1, startPt);

              if (dist0 > dist1)
              {
                for (int r=0; r<halfPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(halfPoints->GetId(r));
              }
              else
              {
                for (int r=0; r<avgPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(avgPoints->GetId(r));
              }
            }
            else if (count[1] == 2)
            {
              vtkNew(vtkIdList, avgPoints);
              double check0[3]; check0[0] = 0.0; check0[1] = 0.0; check0[2] = 0.0;
              for (int j=0; j<2; j++)
              {
                pd->GetCellPoints(cellIds[1][j], npts, pts);
                for (int p=0; p<npts; p++)
                {
                  int isId = avgPoints->IsId(pts[p]);
                  if (isId == -1)
                  {
                    avgPoints->InsertNextId(pts[p]);
                    double pt[3];
                    pd->GetPoint(pts[p], pt);
                    for (int r=0; r<3; r++)
                      check0[r] += pt[r];
                  }
                }
              }
              for (int j=0; j<3; j++)
                check0[j] = (1./avgPoints->GetNumberOfIds())*check0[j];

              vtkNew(vtkIdList, halfPoints);
              double check1[3]; check1[0] = 0.0; check1[1] = 0.0; check1[2] = 0.0;
              pd->GetCellPoints(cellIds[1][0], npts, pts);
              for (int j=0; j<npts; j++)
              {
                int ptId0 = pts[j];
                int ptId1 = pts[(j+1)%npts];
                vtkNew(vtkIdList, neighborCell);
                pd->GetCellEdgeNeighbors(cellIds[1][0], ptId0, ptId1, neighborCell);
                if (neighborCell->GetNumberOfIds() > 0)
                {
                  if (neighborCell->GetId(0) == cellIds[1][1])
                  {
                    halfPoints->InsertNextId(ptId0);
                    halfPoints->InsertNextId(ptId1);

                    double pt[3];
                    pd->GetPoint(ptId0, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    pd->GetPoint(ptId1, pt);
                    for (int r=0; r<3; r++)
                      check1[r] += pt[r];
                    for (int r=0; r<3; r++)
                      check1[r] = (1./2)*check1[r];
                  }
                }
              }

              double startPt[3];
              pd->GetPoint(i, startPt);
              double dist0 = vtkSVMathUtils::Distance(check0, startPt);
              double dist1 = vtkSVMathUtils::Distance(check1, startPt);

              if (dist0 > dist1)
              {
                for (int r=0; r<halfPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(halfPoints->GetId(r));
              }
              else
              {
                for (int r=0; r<avgPoints->GetNumberOfIds(); r++)
                  uniquePoints->InsertNextId(avgPoints->GetId(r));
              }
            }
          }
          int numIds = uniquePoints->GetNumberOfIds();
          double center[3];
          for (int j=0; j<3; j++)
            center[j] = 0.0;
          for (int k=0; k<numIds; k++)
          {
            double pt[3];
            pd->GetPoint(uniquePoints->GetId(k), pt);
            for (int j=0; j<3; j++)
              center[j] += pt[j];
          }
          for (int j=0; j<3; j++)
            center[j] = (1./numIds)*center[j];

          pd->GetPoints()->SetPoint(i, center);
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// GetRegions
// ----------------------
int vtkSVSurfaceCenterlineGrouper::GetRegions(vtkPolyData *pd, std::string arrayName,
                                     std::vector<Region> &allRegions)
{

  int numCells = pd->GetNumberOfCells();
  int numPoints = pd->GetNumberOfPoints();

  std::vector<std::vector<int> > tempRegions(numCells);
  std::vector<std::vector<int> > directNeighbors(numCells);
  std::vector<int> numberOfDirectNeighbors(numCells);
  std::vector<int> pointOnOpenEdge(numPoints, 0);

  for (int i=0; i<numCells; i++)
  {
    int directNeiCount = 0;
    std::vector<int> neighborCells;
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];
      vtkNew(vtkIdList, cellEdgeNeighbors);
      pd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellEdgeNeighbors);
      directNeiCount += cellEdgeNeighbors->GetNumberOfIds();
      for (int k=0; k<cellEdgeNeighbors->GetNumberOfIds(); k++)
        neighborCells.push_back(cellEdgeNeighbors->GetId(k));

      if (cellEdgeNeighbors->GetNumberOfIds() == 0)
      {
        pointOnOpenEdge[ptId0] = 1;
        pointOnOpenEdge[ptId1] = 1;
      }
    }
    directNeighbors[i] = neighborCells;
    numberOfDirectNeighbors[i] = directNeiCount;
  }

  for (int i=0; i<numCells; i++)
  {
    int regionId = pd->GetCellData()->GetArray(
      arrayName.c_str())->GetTuple1(i);
    tempRegions[i].push_back(-1);
    tempRegions[i].push_back(regionId);
  }

  int region = 0;
  for (int i=0; i<numCells; i++)
  {
    if (tempRegions[i][0] == -1)
    {
      tempRegions[i][0] = region;

      int count=1;
      std::vector<int> tempIndex;
      tempIndex.push_back(i);

      for (int j=0; j<count; j++)
      {
        for (int k=0; k<numberOfDirectNeighbors[tempIndex[j]]; k++)
        {
          int cellId = directNeighbors[tempIndex[j]][k];
          if (tempRegions[cellId][0] == -1 && tempRegions[i][1] == tempRegions[cellId][1])
          {
            tempRegions[cellId][0] = region;
            tempIndex.push_back(cellId);
            count++;
          }
        }
      }
      region++;
    }
  }

  int numberOfRegions = region;

  allRegions.clear();
  allRegions.resize(numberOfRegions);

  for (int i=0; i<numberOfRegions; i++)
  {
    allRegions[i].Index = i;
    allRegions[i].IndexCluster = 0;
    allRegions[i].NumberOfCorners = 0;
    allRegions[i].NumberOfElements = 0;
    allRegions[i].Elements.clear();
    allRegions[i].CornerPoints.clear();
    for (int j=0; j<allRegions[i].BoundaryEdges.size(); j++)
      allRegions[i].BoundaryEdges[j].clear();
    allRegions[i].BoundaryEdges.clear();
  }

  for (int i=0; i<numCells; i++)
  {
    int regionId = tempRegions[i][0];
    allRegions[regionId].Elements.push_back(i);
    allRegions[regionId].NumberOfElements++;
  }

  for (int i=0; i<numberOfRegions; i++)
  {
    int cellId = allRegions[i].Elements[0];
    allRegions[i].IndexCluster = tempRegions[cellId][1];
  }

  std::vector<int> cornerPoints;
  std::vector<int> isCornerPoint(numPoints);
  std::vector<int> isBoundaryPoint(numPoints);
  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, pointCellsValues);
    vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, i, pointCellsValues);
    if (pointOnOpenEdge[i] == 1)
      pointCellsValues->InsertNextId(-1);
    if (pointCellsValues->GetNumberOfIds() >= 3)
    {
      cornerPoints.push_back(i);
      isCornerPoint[i] = 1;
    }
    else
      isCornerPoint[i] = 0;

    if (pointCellsValues->GetNumberOfIds() == 2)
      isBoundaryPoint[i] = 1;
    else
      isBoundaryPoint[i] = 0;
  }

  int runCount = 0;
  int numberOfCornerPoints = cornerPoints.size();

  int firstCorner;

  for (int i=0; i<numberOfRegions; i++)
  {
    std::vector<int> tempCornerPoints;
    for (int j=0; j<allRegions[i].NumberOfElements; j++)
    {
      int cellId = allRegions[i].Elements[j];
      vtkIdType npts, *pts;
      pd->GetCellPoints(cellId, npts, pts);
      for (int k=0; k<npts; k++)
      {
        if (isCornerPoint[pts[k]])
        {
          bool kCount = true;
          for (int kk=0; kk<tempCornerPoints.size(); kk++)
          {
            if (pts[k] == tempCornerPoints[kk])
            {
              kCount = false;
            }
          }

          if (kCount == true)
          {
            tempCornerPoints.push_back(pts[k]);
          }
        }
      }
    }

    allRegions[i].NumberOfCorners = tempCornerPoints.size();
    //vtkDebugMacro("NUM CORNS: " << allRegions[i].NumberOfCorners << " OF GROUP " <<  allRegions[i].IndexCluster);

    if (allRegions[i].NumberOfCorners != 0)
    {
      firstCorner = tempCornerPoints[0];
      allRegions[i].CornerPoints.push_back(firstCorner);

      int count=1;
      std::vector<int> tempNodes;
      tempNodes.push_back(firstCorner);

      for (int j=0; j<count; j++)
      {
        vtkNew(vtkIdList, pointCells);
        pd->GetPointCells(tempNodes[j], pointCells);
        for (int k=0; k<pointCells->GetNumberOfIds(); k++)
        {
          int cellId =  pointCells->GetId(k);
          int pointCCWId = vtkSVSurfaceCenterlineGrouper::GetCCWPoint(pd, tempNodes[j], cellId);
          int isBoundaryEdge = vtkSVSurfaceCenterlineGrouper::CheckBoundaryEdge(pd, arrayName, cellId, tempNodes[j], pointCCWId);

          if (tempRegions[cellId][0] == allRegions[i].Index && isBoundaryPoint[pointCCWId] && isBoundaryEdge)
          {
            tempNodes.push_back(pointCCWId);
            count++;
          }
          else if (tempRegions[cellId][0] == allRegions[i].Index && isCornerPoint[pointCCWId] && isBoundaryEdge)
          {
            if (pointCCWId == firstCorner)
            {
              tempNodes.push_back(pointCCWId);
              allRegions[i].BoundaryEdges.push_back(tempNodes);

              tempNodes.clear();

              if (allRegions[i].CornerPoints.size() == allRegions[i].NumberOfCorners)
              {
                count = -1;
                break;
              }
              else
              {
                for (int ii=0; ii<tempCornerPoints.size(); ii++)
                {
                  bool tempCount = false;
                  int tempIndex  = tempCornerPoints[ii];

                  for (int jj=0; jj<allRegions[i].CornerPoints.size(); jj++)
                  {
                    if (tempIndex == allRegions[i].CornerPoints[jj])
                      tempCount = true;
                  }
                  if (tempCount == false)
                  {
                    firstCorner = tempIndex;
                    break;
                  }
                }
                allRegions[i].CornerPoints.push_back(firstCorner);
                tempNodes.push_back(firstCorner);
                count = 1;
                j = -1;
                break;
              }
            }
            else
            {
              tempNodes.push_back(pointCCWId);
              allRegions[i].CornerPoints.push_back(pointCCWId);
              allRegions[i].BoundaryEdges.push_back(tempNodes);
              tempNodes.clear();
              tempNodes.push_back(pointCCWId);
              count = 1;
              j = -1;
              break;
            }
          }
        }
      }
    }
    if (allRegions[i].CornerPoints.size() != allRegions[i].NumberOfCorners)
    {
      //vtkErrorMacro("NUM CORNER POINTS DON'T MATCH: " <<  tempCornerPoints.size() << " " << allRegions[i].CornerPoints.size());
      return SV_ERROR;
    }
  }
  //vtkDebugMacro("DONE GETTING REGIONS");


  return SV_OK;
}

// TODO: Need to fix for if single cell is region!!!
// ----------------------
// GetSpecificRegions
// ----------------------
int vtkSVSurfaceCenterlineGrouper::GetSpecificRegions(vtkPolyData *pd, std::string arrayName,
                                             std::vector<Region> &allRegions,
                                             vtkIdList *targetRegions)
{
  int numCells = pd->GetNumberOfCells();
  int numPoints = pd->GetNumberOfPoints();

  std::vector<std::vector<int> > tempRegions(numCells);
  std::vector<std::vector<int> > directNeighbors(numCells);
  std::vector<int> numberOfDirectNeighbors(numCells);
  std::vector<int> pointOnOpenEdge(numPoints, 0);

  for (int i=0; i<numCells; i++)
  {
    int directNeiCount = 0;
    std::vector<int> neighborCells;
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];
      vtkNew(vtkIdList, cellEdgeNeighbors);
      pd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellEdgeNeighbors);
      directNeiCount += cellEdgeNeighbors->GetNumberOfIds();
      for (int k=0; k<cellEdgeNeighbors->GetNumberOfIds(); k++)
        neighborCells.push_back(cellEdgeNeighbors->GetId(k));

      if (cellEdgeNeighbors->GetNumberOfIds() == 0)
      {
        pointOnOpenEdge[ptId0] = 1;
        pointOnOpenEdge[ptId1] = 1;
      }
    }
    directNeighbors[i] = neighborCells;
    numberOfDirectNeighbors[i] = directNeiCount;
  }

  for (int i=0; i<numCells; i++)
  {
    int regionId = pd->GetCellData()->GetArray(
      arrayName.c_str())->GetTuple1(i);
    tempRegions[i].push_back(-1);
    tempRegions[i].push_back(regionId);
  }

  int region = 0;
  for (int i=0; i<numCells; i++)
  {
    if (tempRegions[i][0] == -1 && targetRegions->IsId(tempRegions[i][1]) != -1)
    {
      tempRegions[i][0] = region;

      int count=1;
      std::vector<int> tempIndex;
      tempIndex.push_back(i);

      for (int j=0; j<count; j++)
      {
        for (int k=0; k<numberOfDirectNeighbors[tempIndex[j]]; k++)
        {
          int cellId = directNeighbors[tempIndex[j]][k];
          if (tempRegions[cellId][0] == -1 && tempRegions[i][1] == tempRegions[cellId][1])
          {
            tempRegions[cellId][0] = region;
            tempIndex.push_back(cellId);
            count++;
          }
        }
      }
      region++;
    }
  }

  int numberOfRegions = region;

  allRegions.clear();
  allRegions.resize(numberOfRegions);

  for (int i=0; i<numberOfRegions; i++)
  {
    allRegions[i].Index = i;
    allRegions[i].IndexCluster = 0;
    allRegions[i].NumberOfCorners = 0;
    allRegions[i].NumberOfElements = 0;
    allRegions[i].Elements.clear();
    allRegions[i].CornerPoints.clear();
    for (int j=0; j<allRegions[i].BoundaryEdges.size(); j++)
      allRegions[i].BoundaryEdges[j].clear();
    allRegions[i].BoundaryEdges.clear();
  }

  for (int i=0; i<numCells; i++)
  {
    int regionId = tempRegions[i][0];
    if (regionId != -1)
    {
      allRegions[regionId].Elements.push_back(i);
      allRegions[regionId].NumberOfElements++;
    }
  }

  for (int i=0; i<numberOfRegions; i++)
  {
    int cellId = allRegions[i].Elements[0];
    allRegions[i].IndexCluster = tempRegions[cellId][1];
  }

  std::vector<int> cornerPoints;
  std::vector<int> isCornerPoint(numPoints);
  std::vector<int> isBoundaryPoint(numPoints);
  std::vector<int> isNonTargetBoundaryPoint(numPoints);
  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, pointCellsValues);
    vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, i, pointCellsValues);
    if (pointOnOpenEdge[i] == 1)
      pointCellsValues->InsertNextId(-1);
    if (pointCellsValues->GetNumberOfIds() >= 3)
    {
      cornerPoints.push_back(i);
      isCornerPoint[i] = 1;
    }
    else
      isCornerPoint[i] = 0;

    if (pointCellsValues->GetNumberOfIds() == 2)
    {
      if (targetRegions->IsId(pointCellsValues->GetId(0)) != -1 &&
          targetRegions->IsId(pointCellsValues->GetId(1)) != -1)
      {
        isBoundaryPoint[i] = 1;
        isNonTargetBoundaryPoint[i] = 0;
      }
      else
      {
        isBoundaryPoint[i] = 0;
        isNonTargetBoundaryPoint[i] = 1;
      }
    }
    else
    {
      isBoundaryPoint[i] = 0;
      isNonTargetBoundaryPoint[i] = 0;
    }
  }

  int runCount = 0;
  int numberOfCornerPoints = cornerPoints.size();

  int firstCorner;

  for (int i=0; i<numberOfRegions; i++)
  {
    std::vector<int> tempCornerPoints;
    for (int j=0; j<allRegions[i].NumberOfElements; j++)
    {
      int cellId = allRegions[i].Elements[j];
      vtkIdType npts, *pts;
      pd->GetCellPoints(cellId, npts, pts);
      for (int k=0; k<npts; k++)
      {
        if (isCornerPoint[pts[k]])
        {
          bool kCount = true;
          for (int kk=0; kk<tempCornerPoints.size(); kk++)
          {
            if (pts[k] == tempCornerPoints[kk])
            {
              kCount = false;
            }
          }

          if (kCount == true)
          {
            tempCornerPoints.push_back(pts[k]);
          }
        }
      }
    }


    allRegions[i].NumberOfCorners = tempCornerPoints.size();
    //vtkDebugMacro("NUM CORNERS: " << allRegions[i].NumberOfCorners);

    if (allRegions[i].NumberOfCorners != 0)
    {
      firstCorner = tempCornerPoints[0];
      allRegions[i].CornerPoints.push_back(firstCorner);

      int count=1;
      std::vector<int> tempNodes;
      tempNodes.push_back(firstCorner);
      std::vector<int> newNodes;
      newNodes.push_back(firstCorner);

      for (int j=0; j<count; j++)
      {
        vtkNew(vtkIdList, pointCells);
        pd->GetPointCells(tempNodes[j], pointCells);
        for (int k=0; k<pointCells->GetNumberOfIds(); k++)
        {
          int cellId =  pointCells->GetId(k);
          int pointCCWId = vtkSVSurfaceCenterlineGrouper::GetCCWPoint(pd, tempNodes[j], cellId);
          int isBoundaryEdge = vtkSVSurfaceCenterlineGrouper::CheckBoundaryEdge(pd, arrayName, cellId, tempNodes[j], pointCCWId);

          if (tempRegions[cellId][0] == allRegions[i].Index && isBoundaryPoint[pointCCWId] && isBoundaryEdge)
          {
            tempNodes.push_back(pointCCWId);
            newNodes.push_back(pointCCWId);
            count++;
          }
          else if (tempRegions[cellId][0] == allRegions[i].Index && isCornerPoint[pointCCWId] && isBoundaryEdge)
          {
            if (pointCCWId == firstCorner)
            {
              tempNodes.push_back(pointCCWId);
              newNodes.push_back(pointCCWId);
              if (newNodes.size() > 2)
                allRegions[i].BoundaryEdges.push_back(newNodes);


              tempNodes.clear();
              newNodes.clear();

              if (allRegions[i].CornerPoints.size() == allRegions[i].NumberOfCorners)
              {
                count = -1;
                break;
              }
              else
              {
                for (int ii=0; ii<tempCornerPoints.size(); ii++)
                {
                  bool tempCount = false;
                  int tempIndex  = tempCornerPoints[ii];

                  for (int jj=0; jj<allRegions[i].CornerPoints.size(); jj++)
                  {
                    if (tempIndex == allRegions[i].CornerPoints[jj])
                      tempCount = true;
                  }
                  if (tempCount == false)
                  {
                    firstCorner = tempIndex;
                    break;
                  }
                }
                allRegions[i].CornerPoints.push_back(firstCorner);
                tempNodes.push_back(firstCorner);
                newNodes.push_back(firstCorner);
                count = 1;
                j = -1;
                break;
              }
            }
            else
            {
              tempNodes.push_back(pointCCWId);
              newNodes.push_back(pointCCWId);
              allRegions[i].CornerPoints.push_back(pointCCWId);
              if (newNodes.size() > 2)
                allRegions[i].BoundaryEdges.push_back(newNodes);
              tempNodes.clear();
              tempNodes.push_back(pointCCWId);
              newNodes.clear();
              newNodes.push_back(pointCCWId);
              count = 1;
              j = -1;
              break;
            }
          }
          else if (tempRegions[cellId][0] == allRegions[i].Index && isNonTargetBoundaryPoint[pointCCWId] && isBoundaryEdge)
          {
            tempNodes.push_back(pointCCWId);
            count++;
          }
        }
      }
    }
    if (allRegions[i].CornerPoints.size() != allRegions[i].NumberOfCorners)
    {
      //vtkErrorMacro("NUM CORNER POINTS DON'T MATCH: " <<  tempCornerPoints.size() << " " << allRegions[i].CornerPoints.size());
      return SV_ERROR;
    }
  }
  //vtkDebugMacro("DONE GETTING REGIONS");
  return SV_OK;
}


// ----------------------
// GetCCWPoint
// ----------------------
int vtkSVSurfaceCenterlineGrouper::GetCCWPoint(vtkPolyData *pd, const int pointId, const int cellId)
{
	int pointCCW;
	int position = 0;

  vtkIdType npts, *pts;
  pd->GetCellPoints(cellId, npts, pts);
	for (int i = 0; i < npts; i++)
	{
		if (pts[i] == pointId)
		{
			position = i;
			break;
		}
	}

  position = (position+1)%npts;
  return pts[position];
}

// ----------------------
// GetCWPoint
// ----------------------
int vtkSVSurfaceCenterlineGrouper::GetCWPoint(vtkPolyData *pd, const int pointId, const int cellId)
{
	int pointCCW;
	int position = 0;

  vtkIdType npts, *pts;
  pd->GetCellPoints(cellId, npts, pts);
	for (int i = 0; i < npts; i++)
	{
		if (pts[i] == pointId)
		{
			position = i;
			break;
		}
	}

  position = (position+npts-1)%npts;
  return pts[position];
}

// ----------------------
// CheckCellValuesEdge
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CheckCellValuesEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1)
{
  vtkNew(vtkIdList, cellEdgeNeighbors);
  pd->GetCellEdgeNeighbors(cellId, pointId0, pointId1, cellEdgeNeighbors);

  vtkNew(vtkIdList, uniqueVals);
  uniqueVals->InsertNextId(pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(cellId));
  for (int i=0; i<cellEdgeNeighbors->GetNumberOfIds(); i++)
  {
    uniqueVals->InsertUniqueId(pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(cellEdgeNeighbors->GetId(i)));
  }

  int isEdge = 0;

  if (uniqueVals->GetNumberOfIds() == 2)
    isEdge = 1;

  return isEdge;
}

// ----------------------
// CheckBoundaryEdge
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CheckBoundaryEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1)
{
  vtkNew(vtkIdList, cellEdgeNeighbors);
  pd->GetCellEdgeNeighbors(cellId, pointId0, pointId1, cellEdgeNeighbors);

  vtkNew(vtkIdList, uniqueVals);
  uniqueVals->InsertNextId(pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(cellId));
  for (int i=0; i<cellEdgeNeighbors->GetNumberOfIds(); i++)
  {
    uniqueVals->InsertUniqueId(pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(cellEdgeNeighbors->GetId(i)));
  }

  if (cellEdgeNeighbors->GetNumberOfIds() == 0)
    uniqueVals->InsertUniqueId(-1);

  int isEdge = 0;

  if (uniqueVals->GetNumberOfIds() == 2)
    isEdge = 1;

  return isEdge;
}

// ----------------------
// CurveFitBoundaries
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CurveFitBoundaries(vtkPolyData *pd, std::string arrayName,
                                     std::vector<Region> allRegions)
{
  int numRegions = allRegions.size();

  std::vector<int> edgeValueCheck;
  for (int i=0; i<numRegions; i++)
  {
    for (int j=0; j<allRegions[i].BoundaryEdges.size(); j++)
    {
      //vtkDebugMacro("Fitting curve edge " << j << " of region " << i);
      int edgeSize = allRegions[i].BoundaryEdges[j].size();

      int edgeValue = 0;
      for (int k=0; k<edgeSize; k++)
        edgeValue += allRegions[i].BoundaryEdges[j][k];

      int usedEdge=0;
      for (int k=0; k<edgeValueCheck.size(); k++)
      {
        if (edgeValue == edgeValueCheck[k])
        {
          usedEdge = 1;
          break;
        }
      }
      if (usedEdge == 1)
        continue;
      else
      {
        edgeValueCheck.push_back(edgeValue);
      }

      int numPoints = edgeSize-1;
      if (numPoints > 4)
      {
        std::vector<double> lengthRatio(edgeSize, 0.0);

        std::vector<XYZ> inputNodes(edgeSize);
        std::vector<XYZ> outputNodes(edgeSize);

        const int sampleSize = 1000;
        std::vector<XYZ> outputRes(sampleSize);

        for (int k=0; k<edgeSize; k++)
        {
          int pointId = allRegions[i].BoundaryEdges[j][k];
          double pt[3];
          pd->GetPoint(pointId, pt);
          inputNodes[k].x = pt[0];
          inputNodes[k].y = pt[1];
          inputNodes[k].z = pt[2];
        }

        int deg = 4;
        std::vector<int> knots(numPoints+deg+1);

        vtkSVSurfaceCenterlineGrouper::SplineKnots(knots, numPoints, deg);

        double totalLength = 0.0;

        for (int k = 1; k < edgeSize; k++)
        {

          int pointId = allRegions[i].BoundaryEdges[j][k];
          int prevPointId = allRegions[i].BoundaryEdges[j][k-1];

          double pt0[3], pt1[3];
          pd->GetPoint(pointId, pt0);
          pd->GetPoint(prevPointId, pt1);

          totalLength += vtkSVMathUtils::Distance(pt0, pt1);
        }

        double tempLength = 0.0;
        for (int k = 1; k < edgeSize; k++)
        {
          int pointId = allRegions[i].BoundaryEdges[j][k];
          int prevPointId = allRegions[i].BoundaryEdges[j][k-1];

          double pt0[3], pt1[3];
          pd->GetPoint(pointId, pt0);
          pd->GetPoint(prevPointId, pt1);

          tempLength += vtkSVMathUtils::Distance(pt0, pt1);

          lengthRatio[k] = tempLength / totalLength;
        }

        SplineCurve(inputNodes, numPoints, knots, deg, outputRes, sampleSize);

        double minDist = VTK_SV_LARGE_DOUBLE;
        int tempCount=0;
        for (int k = 0; k < edgeSize; k++)
        {
          minDist = VTK_SV_LARGE_DOUBLE;
          int pointId = allRegions[i].BoundaryEdges[j][k];
          double pt[3];
          pd->GetPoint(pointId, pt);
          for (int l = 0; l < sampleSize; l++)
          {
            double outputPt[3];
            outputPt[0] = outputRes[l].x;
            outputPt[1] = outputRes[l].y;
            outputPt[2] = outputRes[l].z;

            double dist = vtkSVMathUtils::Distance(pt, outputPt);

            if (dist < minDist)
            {
              minDist = dist;
              tempCount = l;
            }

          }

          double newPoint[3];
          newPoint[0] = outputRes[tempCount].x;
          newPoint[1] = outputRes[tempCount].y;
          newPoint[2] = outputRes[tempCount].z;

          pd->GetPoints()->SetPoint(pointId, newPoint);
        }
      }
    }
  }
  return SV_OK;
}

// ----------------------
// SplineKnots
// ----------------------
void vtkSVSurfaceCenterlineGrouper::SplineKnots(std::vector<int> &u, int n, int t)
{

	int j;

	for (j = 0; j <= n+t; j++)
	{

		if (j < t)
		{
			u[j] = 0;
		}
		else if (j <= n)
		{
			u[j] = j - t + 1;
		}
		else if (j > n)
		{
			u[j] = n - t + 2;
		}

	}

}

// ----------------------
// SplineCurve
// ----------------------
void vtkSVSurfaceCenterlineGrouper::SplineCurve(const std::vector<XYZ> &inp, int n, const std::vector<int> &knots, int t, std::vector<XYZ> &outp, int res)
{

	int i;

	double interval, increment;

	interval = 0.f;
	increment = (n - t + 2) / (double)(res-1);

	for (i = 0; i < res-1; i++)
	{

		SplinePoint(knots, n, t, interval, inp, outp[i]);

		interval += increment;
	}

	outp[res-1] = inp[n];

}

void vtkSVSurfaceCenterlineGrouper::SplinePoint(const std::vector<int> &u, int n, int t, double v, const std::vector<XYZ> &control, XYZ &output)
{

	int k;
	double b;

	output.x = 0.f;
	output.y = 0.f;
	output.z = 0.f;

	for (k = 0; k <= n; k++)
	{
		b = SplineBlend(k, t, u, v);

		output.x += control[k].x * b;
		output.y += control[k].y * b;
		output.z += control[k].z * b;
	}

}

double vtkSVSurfaceCenterlineGrouper::SplineBlend(int k, int t, const std::vector<int> &u, double v)
{

	double value;

	if (t == 1)
	{
		if ((u[k] <= v) && (v < u[k+1]))
			value = 1;
		else
			value = 0;
	}
	else
	{
		if ((u[k+t-1] == u[k]) && (u[k+t] == u[k+1]))
			value = 0;
		else if (u[k+t-1] == u[k])
			value = (u[k+t] - v) / (u[k+t] - u[k+1]) * SplineBlend(k+1,t-1,u,v);
		else if (u[k+t] == u[k+1])
			value = (v - u[k]) / (u[k+t-1] - u[k]) * SplineBlend(k,t-1,u,v);
		else
			value = (v - u[k]) / (u[k+t-1] - u[k]) * SplineBlend(k,t-1,u,v) +
			(u[k+t] - v) / (u[k+t] - u[k+1]) * SplineBlend(k+1,t-1,u,v);
	}

	return(value);

}

// ----------------------
// CheckGroups
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CheckGroups(vtkPolyData *pd)
{
  // Clean up groups
  int allGood = 0;
  int iter = 0;
  int maxIters = 15;
  while(!allGood && iter < maxIters)
  {
    allGood = 1;
    for (int i=0; i<pd->GetNumberOfCells(); i++)
    {
      int groupVal = pd->GetCellData()->GetArray(
        this->GroupIdsArrayName)->GetTuple1(i);
      if (groupVal == -1)
      {
        allGood = 0;;
        vtkNew(vtkIdList, neighborValues);
        vtkSVGeneralUtils::GetNeighborsCellsValues(pd,
                                                   this->GroupIdsArrayName,
                                                   i, neighborValues);
        int newCellValue = -1;
        for (int j=0; j<neighborValues->GetNumberOfIds(); j++)
        {
          if (neighborValues->GetId(j) != -1)
          {
            newCellValue = neighborValues->GetId(j);
            break;
          }
        }
        if (newCellValue != -1)
          pd->GetCellData()->GetArray(this->GroupIdsArrayName)->SetTuple1(i, newCellValue);
      }
    }
    vtkDebugMacro("GROUP FIX ITER: " << iter);
    iter++;
  }

  return SV_OK;
}

// ----------------------
// CheckGroups2
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CheckGroups2()
{
  int allGood = 0;
  int iter = 0;
  vtkPolyData *newMergedCenterlinesPd = NULL;
  int maxIters = 3;
  while(!allGood && iter < maxIters+1)
  {
    allGood = 1;
    // Get all group ids
    vtkNew(vtkIdList, groupIds);
    for (int i=0; i<this->MergedCenterlines->GetNumberOfCells(); i++)
    {
      int groupVal = this->MergedCenterlines->GetCellData()->GetArray(
          this->GroupIdsArrayName)->GetTuple1(i);
      groupIds->InsertUniqueId(groupVal);
    }
    vtkSortDataArray::Sort(groupIds);
    int numGroups = groupIds->GetNumberOfIds();

    vtkNew(vtkThreshold, groupThresholder);
    groupThresholder->SetInputData(this->WorkPd);
    groupThresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName);

    vtkNew(vtkDataSetSurfaceFilter, surfacer);
    vtkNew(vtkConnectivityFilter, connector);
    vtkNew(vtkIdList, backNeighbors);
    vtkNew(vtkIdList, frontNeighbors);
    vtkNew(vtkFeatureEdges, featureEdges);
    vtkNew(vtkIdFilter, ider);

    vtkIdType nlinepts, *linepts;
    int groupId, centerlineId, isTerminating;

    for (int i=0; i<numGroups; i++)
    {
      groupId = groupIds->GetId(i);

      centerlineId = this->MergedCenterlines->GetCellData()->GetArray(this->GroupIdsArrayName)->LookupValue(groupId);
      this->MergedCenterlines->GetCellPoints(centerlineId, nlinepts, linepts);
      isTerminating = 1;

      this->MergedCenterlines->GetPointCells(linepts[0], frontNeighbors);
      vtkNew(vtkIdList, frontGroupNeighbors);
      for (int j=0; j<frontNeighbors->GetNumberOfIds(); j++)
      {
        frontGroupNeighbors->InsertNextId(this->MergedCenterlines->GetCellData()->GetArray(
          this->GroupIdsArrayName)->GetTuple1(frontNeighbors->GetId(j)));
      }

      this->MergedCenterlines->GetPointCells(linepts[nlinepts-1], backNeighbors);
      vtkNew(vtkIdList, backGroupNeighbors);
      for (int j=0; j<backNeighbors->GetNumberOfIds(); j++)
      {
        backGroupNeighbors->InsertNextId(this->MergedCenterlines->GetCellData()->GetArray(
          this->GroupIdsArrayName)->GetTuple1(backNeighbors->GetId(j)));
      }
      if (backNeighbors->GetNumberOfIds() != 1 && frontNeighbors->GetNumberOfIds() != 1)
        isTerminating = 0;

      groupThresholder->ThresholdBetween(groupId, groupId);
      groupThresholder->Update();

      int addNewCenterlines = 0;
      if (groupThresholder->GetOutput()->GetNumberOfPoints() == 0)
      {
        vtkWarningMacro("THERE ARE NO CELLS ON SURFACE FOR GROUP "<< groupId);
        addNewCenterlines = 1;
      }

      surfacer->SetInputData(groupThresholder->GetOutput());
      surfacer->Update();

      ider->SetInputData(surfacer->GetOutput());
      ider->SetIdsArrayName("TmpInternalIds");
      ider->Update();

      connector->SetInputData(ider->GetOutput());
      connector->ColorRegionsOff();
      connector->Update();

      if (connector->GetNumberOfExtractedRegions() > 1)
      {
        vtkWarningMacro("EXPECTED ONE EXTRACTED REGION FOR GROUP "<< groupId << ", BUT THERE ARE " << connector->GetNumberOfExtractedRegions() << " REGIONS");
        addNewCenterlines = 1;
      }

      featureEdges->SetInputData(ider->GetOutput());
      featureEdges->BoundaryEdgesOn();
      featureEdges->FeatureEdgesOff();
      featureEdges->ManifoldEdgesOff();
      featureEdges->NonManifoldEdgesOff();
      featureEdges->Update();

      connector->SetInputData(featureEdges->GetOutput());
      connector->ColorRegionsOn();
      connector->Update();

      if (isTerminating)
      {
        if (connector->GetNumberOfExtractedRegions() != 1)
        {
          vtkErrorMacro("EXPECTED ONE EDGE ON TERMINATING GROUP "<< groupId << ", BUT THERE ARE " << connector->GetNumberOfExtractedRegions() << " EDGES");
          //continue;
        }
      }
      else
      {
        if (connector->GetNumberOfExtractedRegions() != 2)
        {
          vtkWarningMacro("EXPECTED TWO EDGES ON NON-TERMINATING GROUP "<< groupId << ", BUT THERE ARE " << connector->GetNumberOfExtractedRegions() << " EDGES");
          addNewCenterlines = 1;
        }
        else
        {
          // Check that group is more than one cell thick everywhere
          vtkNew(vtkPolyData, groupPd);
          groupPd->DeepCopy(ider->GetOutput());
          vtkNew(vtkPolyData, edgesPd);
          edgesPd->DeepCopy(connector->GetOutput());
          std::vector<int> pointRegion(groupPd->GetNumberOfPoints(), -1);

          int realPtId, regionId;
          for (int j=0; j<edgesPd->GetNumberOfPoints(); j++)
          {
            realPtId = edgesPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(j);
            regionId = edgesPd->GetPointData()->GetArray("RegionId")->GetTuple1(j);
            pointRegion[realPtId] = regionId;
          }

          int cellId;;
          vtkNew(vtkIdList, pointCellIds);
          vtkIdType npts, *pts;
          int badGroup = 0;
          for (int j=0; j<edgesPd->GetNumberOfPoints(); j++)
          {
            realPtId = edgesPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(j);
            regionId = pointRegion[realPtId];

            if (regionId == 0)
            {
              groupPd->GetPointCells(realPtId, pointCellIds);
              for (int k=0; k<pointCellIds->GetNumberOfIds(); k++)
              {
                cellId = pointCellIds->GetId(k);
                groupPd->GetCellPoints(cellId, npts, pts);
                for (int l=0; l<npts; l++)
                {
                  if (pointRegion[pts[l]] == 1)
                  {
                    badGroup = 1;
                    break;
                  }
                }
              }
            }
          }

          if (badGroup)
          {
            vtkWarningMacro("NON-TERMINATING GROUP "<< groupId << " IS ONLY ONE CELL THICK IN AT LEAST ONE LOCATION");
            addNewCenterlines = 1;
          }
        }
      }

      if (addNewCenterlines)
      {
        allGood = 0;
        vtkNew(vtkPolyData, testPd0);
        testPd0->DeepCopy(this->WorkPd);
        testPd0->BuildLinks();

        vtkNew(vtkPolyData, removedGroupsCenterlinesPd);
        if (newMergedCenterlinesPd == NULL)
        {
          newMergedCenterlinesPd = vtkPolyData::New();
          newMergedCenterlinesPd->DeepCopy(this->MergedCenterlines);
        }
        removedGroupsCenterlinesPd->DeepCopy(newMergedCenterlinesPd);
        removedGroupsCenterlinesPd->BuildLinks();

        int cellGroupId, centerlineIdForGroup;
        for (int j=0; j<testPd0->GetNumberOfCells(); j++)
        {
          cellGroupId = testPd0->GetCellData()->GetArray(
           this->GroupIdsArrayName)->GetTuple1(j);

          if (backGroupNeighbors->IsId(cellGroupId) == -1 &&
              frontGroupNeighbors->IsId(cellGroupId) == -1)
          {
            testPd0->DeleteCell(j);
          }
        }
        testPd0->RemoveDeletedCells();
        testPd0->BuildLinks();

        for (int j=0; j<removedGroupsCenterlinesPd->GetNumberOfCells(); j++)
        {
          if (backNeighbors->IsId(j) == -1 &&
              frontNeighbors->IsId(j) == -1)
          {
            removedGroupsCenterlinesPd->DeleteCell(j);
          }
        }
        removedGroupsCenterlinesPd->DeleteCell(centerlineId);

        removedGroupsCenterlinesPd->RemoveDeletedCells();
        removedGroupsCenterlinesPd->BuildLinks();

        int stopCellNumber = ceil(testPd0->GetNumberOfCells()*0.0001);
        vtkNew(vtkSVCenterlinesEdgeWeightedCVT, CVT);
        CVT->SetInputData(testPd0);
        CVT->SetGenerators(removedGroupsCenterlinesPd);
        CVT->SetNumberOfRings(2);
        CVT->SetThreshold(stopCellNumber);
        CVT->SetUseCurvatureWeight(0);
        CVT->SetPatchIdsArrayName(this->GroupIdsArrayName);
        CVT->SetCVTDataArrayName("Normals");
        CVT->SetGroupIdsArrayName(this->GroupIdsArrayName);
        CVT->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
        CVT->SetUsePointNormal(1);
        CVT->SetUseRadiusInformation(0);
        CVT->SetMaximumNumberOfIterations(0);
        CVT->Update();

        testPd0->DeepCopy(CVT->GetOutput());

        if (this->CheckGroups(testPd0) != SV_OK)
        {
          vtkErrorMacro("Error in correcting groups");
          return SV_ERROR;
        }

        if (this->CorrectCellBoundaries(testPd0, this->GroupIdsArrayName) != SV_OK)
        {
          vtkErrorMacro("Could not correcto boundaries of surface");
          return SV_ERROR;
        }

        std::vector<Region> groupRegions;
        if (this->GetRegions(testPd0, this->GroupIdsArrayName, groupRegions) != SV_OK)
        {
          vtkErrorMacro("Couldn't get group regions");
          return SV_ERROR;
        }

        int ptId0, ptId1, cellId0, cellId1, cellGroupId0, cellGroupId1, edgeSize;
        double newPt0[3], newPt1[3], avgPt[3], newCenterPt[3], newPtVec[3], newOuterPt[3];
        std::vector<int> usedPoints(testPd0->GetNumberOfPoints(), 0);
        vtkNew(vtkIdList, cellEdgeNeighbors);
        vtkNew(vtkIdList, pointCellsValues0);
        vtkNew(vtkIdList, pointCellsValuesN);
        vtkNew(vtkPoints, newCenterlinePts);

        int centerPtId = linepts[nlinepts/2];
        this->MergedCenterlines->GetPoint(centerPtId, newCenterPt);

        for (int j=0; j<groupRegions.size(); j++)
        {
          for (int k=0; k<groupRegions[j].BoundaryEdges.size(); k++)
          {
            edgeSize = groupRegions[j].BoundaryEdges[k].size();
            int ptIdBEG = groupRegions[j].BoundaryEdges[k][0];
            int ptIdEND = groupRegions[j].BoundaryEdges[k][edgeSize-1];
            vtkSVGeneralUtils::GetPointCellsValues(testPd0, this->GroupIdsArrayName, ptIdBEG, pointCellsValues0);
            vtkSVGeneralUtils::GetPointCellsValues(testPd0, this->GroupIdsArrayName, ptIdEND, pointCellsValuesN);
            for (int l=0; l<edgeSize-1; l++)
            {
              ptId0 = groupRegions[j].BoundaryEdges[k][l];
              ptId1 = groupRegions[j].BoundaryEdges[k][l+1];

              if (usedPoints[ptId0] && usedPoints[ptId1])
                continue;

              usedPoints[ptId0] = 1;
              usedPoints[ptId1] = 1;

              testPd0->GetCellEdgeNeighbors(-1, ptId0, ptId1, cellEdgeNeighbors);

              if (cellEdgeNeighbors->GetNumberOfIds() == 2)
              {
                cellId0 = cellEdgeNeighbors->GetId(0);
                cellId1 = cellEdgeNeighbors->GetId(1);

                cellGroupId0 = testPd0->GetCellData()->GetArray(
                  this->GroupIdsArrayName)->GetTuple1(cellId0);
                cellGroupId1 = testPd0->GetCellData()->GetArray(
                  this->GroupIdsArrayName)->GetTuple1(cellId1);

                if ((frontGroupNeighbors->IsId(cellGroupId0) != -1 &&
                    backGroupNeighbors->IsId(cellGroupId1) != -1) ||
                    (frontGroupNeighbors->IsId(cellGroupId1) != -1 &&
                    backGroupNeighbors->IsId(cellGroupId0) != -1))
                {
                  testPd0->GetPoint(ptId0, newPt0);
                  testPd0->GetPoint(ptId1, newPt1);

                  vtkMath::Add(newPt0, newPt1, avgPt);
                  vtkMath::MultiplyScalar(avgPt, 1./2);

                  newCenterlinePts->InsertNextPoint(avgPt);
                }
              }
              else if (cellEdgeNeighbors->GetNumberOfIds() == 1)
              {
                int frontGroupIn0 = 0;
                int frontGroupInN = 0;
                int backGroupIn0 = 0;
                int backGroupInN = 0;
                for (int l=0; l<frontGroupNeighbors->GetNumberOfIds(); l++)
                {
                  if (frontGroupNeighbors->GetId(l) == groupId)
                    continue;
                  if (pointCellsValues0->IsId(frontGroupNeighbors->GetId(l)) != -1)
                  {
                    frontGroupIn0 = 1;
                  }
                  if (pointCellsValuesN->IsId(frontGroupNeighbors->GetId(l)) != -1)
                  {
                    frontGroupInN = 1;
                  }
                }
                for (int l=0; l<backGroupNeighbors->GetNumberOfIds(); l++)
                {
                  if (backGroupNeighbors->GetId(l) == groupId)
                    continue;
                  if (pointCellsValues0->IsId(backGroupNeighbors->GetId(l)) != -1)
                  {
                    backGroupIn0 = 1;
                  }
                  if (pointCellsValuesN->IsId(backGroupNeighbors->GetId(l)) != -1)
                  {
                    backGroupInN = 1;
                  }
                }

                if (frontGroupIn0 && frontGroupInN && backGroupIn0 && backGroupInN)
                {
                  testPd0->GetPoint(ptId0, newPt0);
                  testPd0->GetPoint(ptId1, newPt1);

                  vtkMath::Add(newPt0, newPt1, avgPt);
                  vtkMath::MultiplyScalar(avgPt, 1./2);

                  newCenterlinePts->InsertNextPoint(avgPt);
                }
              }
            }
          }
        }

        int numNewCenterlinePts = newCenterlinePts->GetNumberOfPoints();

        vtkNew(vtkPolyData, tmpCenterlinesPd);
        tmpCenterlinesPd->DeepCopy(newMergedCenterlinesPd);
        tmpCenterlinesPd->BuildLinks();

        vtkNew(vtkPoints, newPoints);
        vtkNew(vtkPointData, newPointData);
        newPointData->CopyAllocate(tmpCenterlinesPd->GetPointData(),
                                   tmpCenterlinesPd->GetNumberOfPoints() + 2*numNewCenterlinePts);

        vtkNew(vtkCellArray, newCells);
        vtkNew(vtkCellData, newCellData);
        newCellData->CopyAllocate(tmpCenterlinesPd->GetCellData(),
                                  tmpCenterlinesPd->GetNumberOfPoints() + numNewCenterlinePts);

        int newPointId;
        for (int j=0; j<tmpCenterlinesPd->GetNumberOfPoints(); j++)
        {
          newPointId = newPoints->InsertNextPoint(tmpCenterlinesPd->GetPoint(j));
          newPointData->CopyData(tmpCenterlinesPd->GetPointData(), j, newPointId);
        }

        int newCellId;
        for (int j=0; j<tmpCenterlinesPd->GetNumberOfCells(); j++)
        {
          vtkNew(vtkPolyLine, newLine);
          newLine->GetPointIds()->DeepCopy(tmpCenterlinesPd->GetCell(j)->GetPointIds());
          newCellId = newCells->InsertNextCell(newLine);
          newCellData->CopyData(tmpCenterlinesPd->GetCellData(), j, newCellId);
        }

        double radiusValAtCenter = tmpCenterlinesPd->GetPointData()->GetArray(this->CenterlineRadiusArrayName)->GetTuple1(centerPtId);
        vtkNew(vtkPolyLine, newPts0);
        for (int j=0; j<newCenterlinePts->GetNumberOfPoints(); j++)
        {
          ptId0 = newPoints->InsertNextPoint(newCenterPt);
          ptId1 = newPoints->InsertNextPoint(newCenterlinePts->GetPoint(j));

          vtkNew(vtkPolyLine, newPts0);
          newPts0->GetPointIds()->InsertNextId(ptId0);
          newPts0->GetPointIds()->InsertNextId(ptId1);
          newPointData->CopyData(tmpCenterlinesPd->GetPointData(), centerPtId, ptId0);
          newPointData->CopyData(tmpCenterlinesPd->GetPointData(), centerPtId, ptId1);
          newPointData->GetArray(this->CenterlineRadiusArrayName)->SetTuple1(ptId1, 0.2*radiusValAtCenter);

          newCellId = newCells->InsertNextCell(newPts0);
          newCellData->CopyData(tmpCenterlinesPd->GetCellData(), centerlineId, newCellId);
        }

        newPointData->Squeeze();
        newCellData->Squeeze();

        tmpCenterlinesPd->Reset();
        tmpCenterlinesPd->SetPoints(newPoints);
        tmpCenterlinesPd->SetLines(newCells);

        tmpCenterlinesPd->GetPointData()->PassData(newPointData);
        tmpCenterlinesPd->GetCellData()->PassData(newCellData);
        tmpCenterlinesPd->BuildLinks();

        newMergedCenterlinesPd->DeepCopy(tmpCenterlinesPd);

        //continue;
      }

    }

    if (!allGood && iter < maxIters)
    {
      vtkSVIOUtils::WriteVTPFile("/Users/adamupdegrove/Desktop/tmp/MYSURFACEDUMBNEWCENTERLINES.vtp", newMergedCenterlinesPd);
      // Now re-segment with these new centerlines
      int stopCellNumber = ceil(this->WorkPd->GetNumberOfCells()*0.0001);
      vtkNew(vtkSVCenterlinesEdgeWeightedCVT, betterCVT);
      betterCVT->SetInputData(this->WorkPd);
      betterCVT->SetGenerators(newMergedCenterlinesPd);
      betterCVT->SetNumberOfRings(2);
      betterCVT->SetThreshold(stopCellNumber);
      betterCVT->SetUseCurvatureWeight(0);
      betterCVT->SetPatchIdsArrayName(this->GroupIdsArrayName);
      betterCVT->SetCVTDataArrayName("Normals");
      betterCVT->SetGroupIdsArrayName(this->GroupIdsArrayName);
      betterCVT->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
      betterCVT->SetUsePointNormal(1);
      betterCVT->SetUseRadiusInformation(this->UseRadiusInformation);
      betterCVT->SetMaximumNumberOfIterations(0);
      betterCVT->Update();

      this->WorkPd->DeepCopy(betterCVT->GetOutput());

      if (this->CheckGroups(this->WorkPd) != SV_OK)
      {
        vtkErrorMacro("Error in correcting groups");
        return SV_ERROR;
      }

      if (this->CorrectCellBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
      {
        vtkErrorMacro("Could not correcto boundaries of surface");
        return SV_ERROR;
      }
    }
    iter++;
  }

  if (newMergedCenterlinesPd != NULL)
  {
    newMergedCenterlinesPd->Delete();
  }

  if (!allGood)
  {
    vtkErrorMacro("Correction of groups failed");
    return SV_ERROR;
  }

  return SV_OK;
}

// ----------------------
// FixGroupsWithPolycube
// ----------------------
int vtkSVSurfaceCenterlineGrouper::FixGroupsWithPolycube()
{
  // Then check everything
  // Extract surface, triangulate, and subdivide polycube
  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(this->PolycubePd);
  triangulator->Update();

  vtkNew(vtkPolyData, polycubePd);
  polycubePd->DeepCopy(triangulator->GetOutput());
  polycubePd->BuildLinks();

  int deletedCell = 0;
  for (int i=0; i<polycubePd->GetNumberOfCells(); i++)
  {
    // Check for non-manifold cell, if found, delete (just the one).
    vtkIdType npts, *pts;
    polycubePd->GetCellPoints(i, npts, pts);

    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];

      vtkNew(vtkIdList, cellIds);
      polycubePd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellIds);

      if (cellIds->GetNumberOfIds() > 1)
      {
        // Mark for deletion! but not the current cell
        for (int k=0; k<cellIds->GetNumberOfIds(); k++)
        {
          vtkIdType npts_new, *pts_new;
          polycubePd->GetCellPoints(cellIds->GetId(k), npts_new, pts_new);

          int ptFound = 0;
          for (int l=0; l<npts; l++)
          {
            for (int m=0; m<npts_new; m++)
            {
              if (pts[l] == pts_new[m])
                ptFound++;
            }
          }

          if (ptFound == npts)
          {
            vtkDebugMacro("DELETING CELLS: " <<  i << " " << cellIds->GetId(k));
            polycubePd->DeleteCell(i);
            polycubePd->DeleteCell(cellIds->GetId(k));
            deletedCell = 1;
          }
        }
      }
    }
  }
  // Then maybe re-triangulate?
  if (deletedCell)
  {
    polycubePd->RemoveDeletedCells();

    vtkNew(vtkCleanPolyData, cleaner);
    cleaner->SetInputData(polycubePd);
    cleaner->ToleranceIsAbsoluteOn();
    cleaner->SetAbsoluteTolerance(1.0e-6);
    cleaner->Update();

    polycubePd->DeepCopy(cleaner->GetOutput());;
    polycubePd->BuildLinks();
  }

  vtkDebugMacro("GETTING SURFACE GROUPS");
  std::vector<Region> surfaceGroups;
  if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, surfaceGroups) != SV_OK)
  {
    vtkErrorMacro("Couldn't get groups");
    return SV_ERROR;
  }

  vtkDebugMacro("GETTING POLYCUBE GROUPS");
  std::vector<Region> polycubeGroups;
  if (this->GetRegions(polycubePd, this->GroupIdsArrayName, polycubeGroups) != SV_OK)
  {
    vtkErrorMacro("Couldn't get groups");
    return SV_ERROR;
  }

  int numSurfaceGroups  = surfaceGroups.size();
  int numPolycubeGroups = polycubeGroups.size();

  if (numSurfaceGroups != numPolycubeGroups)
  {
    vtkDebugMacro("NUMBER OF SURFACE GROUPS: " << numSurfaceGroups);
    vtkDebugMacro("NUMBER OF POLYCUBE GROUPS: " << numPolycubeGroups);

    if (numSurfaceGroups > numPolycubeGroups)
    {
      vtkDebugMacro("ADDITIONAL SURFACE GROUPS, SEE IF WE CAN REDUCE");
      if (this->FixMultipleGroups(this->WorkPd, polycubePd, surfaceGroups, polycubeGroups) != SV_OK)
      {
        return SV_ERROR;
      }

      vtkDebugMacro("RE-GETTING SURFACE GROUPS");
      if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, surfaceGroups) != SV_OK)
      {
        vtkErrorMacro("Couldn't get groups");
        return SV_ERROR;
      }

      vtkDebugMacro("RE-GETTING POLYCUBE GROUPS");
      if (this->GetRegions(polycubePd, this->GroupIdsArrayName, polycubeGroups) != SV_OK)
      {
        vtkErrorMacro("Couldn't get groups");
        return SV_ERROR;
      }
    }
    else
    {
      vtkDebugMacro("NOT ENOUGH SURFACE REGIONS TO MATCH POLYCUBE");
      return SV_ERROR;
    }
  }

  vtkSVGeneralUtils::GiveIds(this->WorkPd, "TmpInternalIds");
  vtkNew(vtkPolyData, origPd);
  origPd->DeepCopy(this->WorkPd);
  origPd->BuildLinks();

  vtkNew(vtkIdList, critPts);

  int fixed=0;
  for (int i=0; i<numSurfaceGroups; i++)
  {
    for (int j=0; j<numPolycubeGroups; j++)
    {
      if (surfaceGroups[i].IndexCluster == polycubeGroups[j].IndexCluster)
      {
        vtkDebugMacro("NUMRBS OF EDGES FROM SURFACE GROUP " << surfaceGroups[i].IndexCluster << " is " <<  surfaceGroups[i].BoundaryEdges.size());
        vtkDebugMacro("NUMRBS OF EDGES FROM POLYCUBE GROUP " << polycubeGroups[j].IndexCluster << " is " << polycubeGroups[j].BoundaryEdges.size());
        for (int k=0; k<surfaceGroups[i].CornerPoints.size(); k++)
        {
          int cornerPtId = surfaceGroups[i].CornerPoints[k];

          vtkNew(vtkIdList, surfaceCellList);
          vtkSVGeneralUtils::GetPointCellsValues(origPd, this->GroupIdsArrayName, cornerPtId, surfaceCellList);

          vtkDebugMacro("SURFACE CORNER POINT " << k << " GROUPS ARE ");
          for (int l=0; l<surfaceCellList->GetNumberOfIds(); l++)
            vtkDebugMacro(" " << surfaceCellList->GetId(l) << " ");
          vtkDebugMacro("\n");
        }
        for (int k=0; k<polycubeGroups[j].CornerPoints.size(); k++)
        {
          int cornerPtId = polycubeGroups[j].CornerPoints[k];

          vtkNew(vtkIdList, polycubeCellList);
          vtkSVGeneralUtils::GetPointCellsValues(polycubePd, this->GroupIdsArrayName, cornerPtId, polycubeCellList);

          vtkDebugMacro("POLYCUBE CORNER POINT " << k << " GROUPS ARE ");
          for (int l=0; l<polycubeCellList->GetNumberOfIds(); l++)
            vtkDebugMacro(" " << polycubeCellList->GetId(l) << " ");
          vtkDebugMacro("\n");
        }

        std::vector<std::vector<int> > surfConnectedCornerPts;
        this->GetConnectedEdges(surfaceGroups[i].BoundaryEdges, surfConnectedCornerPts);

        vtkDebugMacro("NUMBER OF CONNECTED EDGES: " << surfConnectedCornerPts.size());
        int edgeCount=0;
        for (int k=0; k<surfConnectedCornerPts.size(); k++)
        {
          std::vector<int> badEdges;
          std::vector<int> allEdges;
          for (int l=0; l<surfConnectedCornerPts[k].size(); l++)
          {
            int cornerPtId = surfConnectedCornerPts[k][l];

            vtkNew(vtkIdList, surfaceCellList);
            vtkSVGeneralUtils::GetPointCellsValues(origPd, this->GroupIdsArrayName, cornerPtId, surfaceCellList);

            int foundMatch = 0;
            for (int m=0; m<polycubeGroups[j].CornerPoints.size(); m++)
            {
              int polyCornerPtId = polycubeGroups[j].CornerPoints[m];

              vtkNew(vtkIdList, polycubeCellList);
              vtkSVGeneralUtils::GetPointCellsValues(polycubePd, this->GroupIdsArrayName, polyCornerPtId, polycubeCellList);

              vtkNew(vtkIdList, intersectList);
              intersectList->DeepCopy(polycubeCellList);

              intersectList->IntersectWith(surfaceCellList);

              if (surfaceCellList->GetNumberOfIds() == intersectList->GetNumberOfIds() &&
                  polycubeCellList->GetNumberOfIds() == intersectList->GetNumberOfIds())
              {
                vtkDebugMacro("WE DID IT!, WE FOUND A MATCH!");
                foundMatch = 1;
              }
           }

            if (foundMatch == 0)
            {
              vtkDebugMacro("UH OH, DIDNT FIND MATCH!!!");
              badEdges.push_back(edgeCount);
            }
            allEdges.push_back(edgeCount);
            edgeCount++;
          }
          vtkDebugMacro("NUMBER BAD EDGES: " << badEdges.size());
          vtkDebugMacro("NUMBER ALL EDGES: " << allEdges.size());

          if (surfaceGroups[i].BoundaryEdges.size() == polycubeGroups[j].BoundaryEdges.size() && badEdges.size() == 0)
          {
            vtkDebugMacro("WE GUCCI");
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 4 && polycubeGroups[j].BoundaryEdges.size() == 2 &&
              badEdges.size() == 4 && allEdges.size() == 4)
          {
            this->FixPlanarTrifurcation(this->WorkPd, origPd, this->GroupIdsArrayName, surfaceGroups[i], allEdges, badEdges, critPts);
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 8 && polycubeGroups[j].BoundaryEdges.size() == 4 &&
              badEdges.size() == 4 && allEdges.size() == 4)
          {
            this->FixPlanarTrifurcation(this->WorkPd, origPd, this->GroupIdsArrayName, surfaceGroups[i], allEdges, badEdges, critPts);
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 6 && polycubeGroups[j].BoundaryEdges.size() == 4 &&
              badEdges.size() == 4 && allEdges.size() == 4)
          {
            this->FixPlanarTrifurcation(this->WorkPd, origPd, this->GroupIdsArrayName, surfaceGroups[i], allEdges, badEdges, critPts);
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 3 && polycubeGroups[j].BoundaryEdges.size() == 2 &&
              badEdges.size() == 2 && allEdges.size() == 3)
          {
            this->FixPerpenTrifurcation(this->WorkPd, origPd, this->GroupIdsArrayName, surfaceGroups[i], allEdges, badEdges, critPts);
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 5 && polycubeGroups[j].BoundaryEdges.size() == 4 &&
              badEdges.size() == 2 && allEdges.size() == 3)
          {
            this->FixPerpenTrifurcation(this->WorkPd, origPd, this->GroupIdsArrayName, surfaceGroups[i], allEdges, badEdges, critPts);
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 3 && polycubeGroups[j].BoundaryEdges.size() == 2 &&
              badEdges.size() == 3 && allEdges.size() == 3)
          {
            this->FixOffsetTrifurcation(this->WorkPd, origPd, polycubePd, this->GroupIdsArrayName, surfaceGroups[i], polycubeGroups[j], allEdges, badEdges, critPts);
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 5 && polycubeGroups[j].BoundaryEdges.size() == 4 &&
              badEdges.size() == 3 && allEdges.size() == 3)
          {
            this->FixOffsetTrifurcation(this->WorkPd, origPd, polycubePd, this->GroupIdsArrayName, surfaceGroups[i], polycubeGroups[j], allEdges, badEdges, critPts);
          }
          else if (surfaceGroups[i].BoundaryEdges.size() == 4 && polycubeGroups[j].BoundaryEdges.size() == 2 &&
              badEdges.size() == 2 && allEdges.size() == 4)
          {
            this->FixCloseGroup(this->WorkPd, origPd, polycubePd, this->GroupIdsArrayName, surfaceGroups[i], polycubeGroups[j], allEdges, badEdges, critPts);
          }
          else
            vtkDebugMacro("NO FIX FOR THIS HAS BEEN DEVISED YET!!!!");
        }
      }
    }
    vtkDebugMacro("\n");
  }

  this->WorkPd->GetCellData()->RemoveArray("TmpInternalIds");
  this->WorkPd->GetPointData()->RemoveArray("TmpInternalIds");

  if (this->CorrectCellBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro("Could not correcto boundaries of surface");
    return SV_ERROR;
  }

  if (this->SmoothBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro("Could not smootho boundaries of surface");
    return SV_ERROR;
  }

  vtkDebugMacro("TOTAL NUM OF CELLS: " << this->WorkPd->GetNumberOfCells());
  for (int i=0; i<critPts->GetNumberOfIds(); i++)
  {
    int halfId = critPts->GetId(i);

    vtkNew(vtkIdList, finalVals);
    vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, halfId, finalVals);
    if (finalVals->GetNumberOfIds() != 4)
    {
        vtkDebugMacro("NO GOOD, FIX GROUPS AROUND POINT: " << halfId);
      this->SplitCellsAroundPoint(this->WorkPd, halfId);
      this->SplitCellsAroundPoint(origPd, halfId);
    }
  }

  // Get new normals
  vtkNew(vtkPolyDataNormals, normaler);
  normaler->SetInputData(this->WorkPd);
  normaler->ComputePointNormalsOff();
  normaler->ComputeCellNormalsOn();
  normaler->SplittingOff();
  normaler->Update();
  this->WorkPd->DeepCopy(normaler->GetOutput());
  this->WorkPd->BuildLinks();

  for (int i=0; i<critPts->GetNumberOfIds(); i++)
  {
    int halfId = critPts->GetId(i);

    vtkNew(vtkIdList, checkVals);
    vtkSVGeneralUtils::GetPointCellsValues(origPd, this->GroupIdsArrayName, halfId, checkVals);

    vtkNew(vtkIdList, finalVals);
    vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, halfId, finalVals);
    if (finalVals->GetNumberOfIds() != 4)
    {

      vtkNew(vtkIdList, missingVals);
      for (int j=0; j<checkVals->GetNumberOfIds(); j++)
      {
        if (finalVals->IsId(checkVals->GetId(j)) == -1)
          missingVals->InsertNextId(checkVals->GetId(j));
      }

      vtkNew(vtkIdList, badVals);
      for (int j=0; j<finalVals->GetNumberOfIds(); j++)
      {
        if (checkVals->IsId(finalVals->GetId(j)) == -1)
          badVals->InsertNextId(finalVals->GetId(j));
      }

      if (badVals->GetNumberOfIds() == 2)
      {
        std::vector<std::vector<int> > allNodes;

        int count=1;
        std::vector<int> tempNodes;
        tempNodes.push_back(halfId);

        for (int j=0; j<count; j++)
        {
          vtkNew(vtkIdList, badCells);
          this->WorkPd->GetPointCells(tempNodes[j], badCells);

          for (int k=0; k<badCells->GetNumberOfIds(); k++)
          {
            int cellId = badCells->GetId(k);
            int pointCCWId = vtkSVSurfaceCenterlineGrouper::GetCCWPoint(this->WorkPd, tempNodes[j], cellId);

            vtkNew(vtkIdList, cellEdgeNeighbors);
            this->WorkPd->GetCellEdgeNeighbors(cellId, tempNodes[j], pointCCWId, cellEdgeNeighbors);

            int edgeVal0 = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(cellId);
            int edgeVal1 = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(cellEdgeNeighbors->GetId(0));

            if ((edgeVal0 == badVals->GetId(0) && edgeVal1 == badVals->GetId(1)) ||
                 (edgeVal0 == badVals->GetId(1) && edgeVal1 == badVals->GetId(0)))
            {
              tempNodes.push_back(pointCCWId);
              count++;
            }
            else if (missingVals->IsId(edgeVal0) != -1 || missingVals->IsId(edgeVal1) != -1)
            {
              allNodes.push_back(tempNodes);

              if (allNodes.size() == missingVals->GetNumberOfIds())
              {
                count = -1;
                break;
              }

              tempNodes.clear();
              tempNodes.push_back(halfId);
              count=1;
              j=-1;
              break;
            }
          }
        }

        for (int j=0; j<allNodes.size(); j++)
        {
          for (int k=1; k<allNodes[j].size(); k++)
          {
            int ptId = allNodes[j][k];

            vtkNew(vtkIdList, pointCells);
            this->WorkPd->GetPointCells(ptId, pointCells);

            for (int l=0; l<pointCells->GetNumberOfIds(); l++)
            {
              int cellVal = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(pointCells->GetId(l));
              if (cellVal == badVals->GetId(0) || cellVal == badVals->GetId(1))
              {
                this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->SetTuple1(pointCells->GetId(l), missingVals->GetId(j));
              }
            }
          }
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// FixMultipleGroups
// ----------------------
int vtkSVSurfaceCenterlineGrouper::FixMultipleGroups(vtkPolyData *pd, vtkPolyData *polycubePd, std::vector<Region> surfaceGroups, std::vector<Region> polycubeGroups)
{
  int numSurfaceGroups  = surfaceGroups.size();
  int numPolycubeGroups = polycubeGroups.size();

  vtkNew(vtkIdList, groupIds);
  vtkNew(vtkIdList, groupCount);

  for (int i=0; i<numSurfaceGroups; i++)
  {
    int groupVal = surfaceGroups[i].IndexCluster;
    int isId = groupIds->IsId(groupVal);
    if (isId == -1)
    {
      groupIds->InsertNextId(groupVal);
      groupCount->InsertNextId(1);
    }
    else
      groupCount->SetId(isId, groupCount->GetId(isId)+1);
  }

  for (int i=0; i<groupIds->GetNumberOfIds(); i++)
  {
    int numOfGroup = groupCount->GetId(i);
    if (numOfGroup > 1)
    {
      if (numOfGroup == 2)
      {
        int groupVal = groupIds->GetId(i);
        int polyIndex = -1;
        for (int j=0; j<numPolycubeGroups; j++)
        {
          if (polycubeGroups[j].IndexCluster == groupVal)
            polyIndex = j;
        }

        std::vector<int> groupLocs;
        for (int j=0; j<numSurfaceGroups; j++)
        {
          if (surfaceGroups[j].IndexCluster == groupIds->GetId(i))
            groupLocs.push_back(j);
        }

        std::vector<std::vector<int> > polyConnectedCornerPts;
        this->GetConnectedEdges(polycubeGroups[polyIndex].BoundaryEdges, polyConnectedCornerPts);

        std::vector<std::vector<int> > allPolycubeGroups;
        for (int j=0; j<polyConnectedCornerPts.size(); j++)
        {
          std::vector<int> oneGroupList;
          for (int k=0; k<polyConnectedCornerPts[j].size(); k++)
          {
            int polyCornerPtId = polyConnectedCornerPts[j][k];

            vtkNew(vtkIdList, polyCellList);
            vtkSVGeneralUtils::GetPointCellsValues(polycubePd, this->GroupIdsArrayName, polyCornerPtId, polyCellList);

            for (int l=0; l<polyCellList->GetNumberOfIds(); l++)
              oneGroupList.push_back(polyCellList->GetId(l));
          }
          allPolycubeGroups.push_back(oneGroupList);
        }

        for (int j=0; j<surfaceGroups[groupLocs[0]].CornerPoints.size(); j++)
        {
          int group0CornerPtId = surfaceGroups[groupLocs[0]].CornerPoints[j];

          vtkNew(vtkIdList, cellList0);
          vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, group0CornerPtId, cellList0);

          for (int k=0; k<surfaceGroups[groupLocs[1]].CornerPoints.size(); k++)
          {

            int group1CornerPtId = surfaceGroups[groupLocs[1]].CornerPoints[k];

            vtkNew(vtkIdList, cellList1);
            vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, group1CornerPtId, cellList1);

            vtkNew(vtkIdList, intersectList);
            intersectList->DeepCopy(cellList1);

            intersectList->IntersectWith(cellList0);

            if (cellList0->GetNumberOfIds() == intersectList->GetNumberOfIds() &&
                cellList1->GetNumberOfIds() == intersectList->GetNumberOfIds())
            {
              for (int m=0; m<allPolycubeGroups.size(); m++)
              {
                vtkNew(vtkIdList, polycubeGroupsList);
                for (int n=0; n<allPolycubeGroups[m].size(); n++)
                  polycubeGroupsList->InsertUniqueId(allPolycubeGroups[m][n]);

                vtkNew(vtkIdList, newIntersectList);
                newIntersectList->DeepCopy(intersectList);

                newIntersectList->IntersectWith(polycubeGroupsList);

                vtkDebugMacro("INTERSECTED: ");
                for (int o=0; o<newIntersectList->GetNumberOfIds(); o++)
                  vtkDebugMacro(" " << newIntersectList->GetId(o) << " ");
                vtkDebugMacro("\n");

                if (newIntersectList->GetNumberOfIds() == 2)
                {
                  int fixGroup = -1;
                  if (newIntersectList->GetId(0) == groupVal)
                    fixGroup = newIntersectList->GetId(1);
                  else
                    fixGroup = newIntersectList->GetId(0);

                  int fixGroupId = -1;
                  for (int n=0; n<numSurfaceGroups; n++)
                  {
                    if (surfaceGroups[n].IndexCluster == fixGroup)
                      fixGroupId = n;
                  }

                  for (int n=0; n<surfaceGroups[fixGroupId].BoundaryEdges.size(); n++)
                  {
                    int edgeSize = surfaceGroups[fixGroupId].BoundaryEdges[n].size();

                    int edgePtId0 = surfaceGroups[fixGroupId].BoundaryEdges[n][0];
                    int edgePtIdN = surfaceGroups[fixGroupId].BoundaryEdges[n][edgeSize-1];

                    if ((edgePtId0 == group0CornerPtId && edgePtIdN == group1CornerPtId) ||
                        (edgePtId0 == group1CornerPtId && edgePtIdN == group0CornerPtId))
                    {
                      for (int o=0; o<edgeSize; o++)
                      {
                        int ptId0 = surfaceGroups[fixGroupId].BoundaryEdges[n][o];

                        vtkNew(vtkIdList, pointCellIds);
                        this->WorkPd->GetPointCells(ptId0, pointCellIds);

                        for (int p=0; p<pointCellIds->GetNumberOfIds(); p++)
                        {
                          this->WorkPd->GetCellData()->GetArray(
                            this->GroupIdsArrayName)->SetTuple1(pointCellIds->GetId(p), groupVal);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      else
      {
        vtkDebugMacro("THREE PATCHES OF ONE GROUP, CANNOT HANDLE THIS");
        return SV_ERROR;
      }
    }
  }

  return SV_OK;
}

// ----------------------
// GetConnectedEdges
// ----------------------
int vtkSVSurfaceCenterlineGrouper::GetConnectedEdges(std::vector<std::vector<int> > inputEdges, std::vector<std::vector<int> > &connectedCornerPts)
{
  int numEdges = inputEdges.size();
  int edgeCount = 0;

  while(edgeCount < numEdges)
  {
    std::vector<int> tmpCornerPts;

    int ptId0 = inputEdges[edgeCount][0];

    int ptIdN = -1;

    while (ptIdN != ptId0)
    {
      int edgePtId0 = inputEdges[edgeCount][0];
      tmpCornerPts.push_back(edgePtId0);

      int edgeSize = inputEdges[edgeCount].size();
      ptIdN = inputEdges[edgeCount][edgeSize-1];

      edgeCount++;
    }

    connectedCornerPts.push_back(tmpCornerPts);
  }

  return SV_OK;
}

// ----------------------
// FixCloseGroup
// ----------------------
int vtkSVSurfaceCenterlineGrouper::FixCloseGroup(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd, std::string arrayName,
                                        const Region region, const Region polyRegion, std::vector<int> allEdges,
                                        std::vector<int> badEdges, vtkIdList *critPts)
{
  vtkDebugMacro("--FIX CLOSE GROUP--");
  int patchValue = region.IndexCluster;
  vtkDebugMacro("CLUSTER " << patchValue);
  int numEdges = region.BoundaryEdges.size();
  vtkDebugMacro("NUM EDGES " << numEdges);

  vtkDebugMacro("NUMBER OF ALL EDGES: " << allEdges.size());
  vtkDebugMacro("ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    vtkDebugMacro(" " <<  allEdges[j] << " ");
  vtkDebugMacro("\n");
  vtkDebugMacro("NUMBER OF BAD EDGES: " << badEdges.size());
  vtkDebugMacro("BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    vtkDebugMacro(" " <<  badEdges[j] << " ");
  vtkDebugMacro("\n");

  std::vector<int> fixEdges;
  int newCellValue = -1;
  for (int j=0; j<badEdges.size(); j++)
  {
    int edgeSize = region.BoundaryEdges[badEdges[j]].size();
    int cornerPtId0 = region.BoundaryEdges[badEdges[j]][0];
    int cornerPtId1 = region.BoundaryEdges[badEdges[j]][edgeSize-1];

    vtkNew(vtkIdList, ptId0List);
    vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, cornerPtId0, ptId0List);

    vtkNew(vtkIdList, ptId1List);
    vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, cornerPtId1, ptId1List);

    int foundMatch0 = 0;
    int foundMatch1 = 0;
    for (int k=0; k<polyRegion.CornerPoints.size(); k++)
    {
      int polyCornerPtId = polyRegion.CornerPoints[k];

      vtkNew(vtkIdList, polycubeCellList);
      vtkSVGeneralUtils::GetPointCellsValues(polyPd, this->GroupIdsArrayName, polyCornerPtId, polycubeCellList);

      vtkNew(vtkIdList, intersect0List);
      intersect0List->DeepCopy(polycubeCellList);

      intersect0List->IntersectWith(ptId0List);
      if (intersect0List->GetNumberOfIds() == polycubeCellList->GetNumberOfIds() &&
          intersect0List->GetNumberOfIds() == ptId0List->GetNumberOfIds())
      {
        foundMatch0 = 1;
      }

      vtkNew(vtkIdList, intersect1List);
      intersect1List->DeepCopy(polycubeCellList);

      intersect1List->IntersectWith(ptId1List);
      intersect1List->IntersectWith(ptId1List);
      if (intersect1List->GetNumberOfIds() == polycubeCellList->GetNumberOfIds() &&
          intersect1List->GetNumberOfIds() == ptId1List->GetNumberOfIds())
      {
        foundMatch1 = 1;
      }

      if (newCellValue == -1)
      {
        if (intersect0List->GetNumberOfIds() == 2)
        {
          if (intersect0List->GetId(0) == polyRegion.IndexCluster)
            newCellValue = intersect0List->GetId(1);
          else
            newCellValue = intersect0List->GetId(0);
        }
        else if (intersect1List->GetNumberOfIds() == 2)
        {
          if (intersect1List->GetId(0) == polyRegion.IndexCluster)
            newCellValue = intersect1List->GetId(1);
          else
            newCellValue = intersect1List->GetId(0);
        }
      }

    }

    if (!foundMatch0 && !foundMatch1)
    {
      fixEdges.push_back(badEdges[j]);
    }

  }

  if (newCellValue == -1)
  {
    vtkErrorMacro("Could not get new cell value to use for edge of group");
  }

  vtkDebugMacro("NUMBER OF FIX EDGES: " << fixEdges.size());
  vtkDebugMacro("FIX EDGES: ");
  for (int j=0; j<fixEdges.size(); j++)
  {
    int edgeSize = region.BoundaryEdges[fixEdges[j]].size();
    int startPtId = region.BoundaryEdges[fixEdges[j]][0];
    int endPtId = region.BoundaryEdges[fixEdges[j]][edgeSize-1];
    vtkDebugMacro(" " << fixEdges[j] << " Start Pt: " <<   startPtId << " End Pt: " << endPtId);
  }

  vtkNew(vtkIdList, fixEdgeCells);
  for (int j=0; j<fixEdges.size(); j++)
  {
    int edgeSize = region.BoundaryEdges[fixEdges[j]].size();
    for (int k=0; k<edgeSize; k++)
    {
      int ptId = region.BoundaryEdges[fixEdges[j]][k];

      vtkNew(vtkIdList, pointCellIds);
      pd->GetPointCells(ptId, pointCellIds);

      for (int l=0; l<pointCellIds->GetNumberOfIds(); l++)
      {
        int cellId = pointCellIds->GetId(l);
        int cellValue = pd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(cellId);

        if (cellValue != newCellValue && cellValue != polyRegion.IndexCluster)
          fixEdgeCells->InsertUniqueId(cellId);
      }
    }
  }

  vtkDebugMacro("FILLING IN " << newCellValue << " BUT NOT TOUCHING " <<  newCellValue << " and " << polyRegion.IndexCluster);

  std::vector<std::vector<int> > ringNeighbors(fixEdgeCells->GetNumberOfIds());
  for (int i=0; i<fixEdgeCells->GetNumberOfIds(); i++)
    ringNeighbors[i].push_back(fixEdgeCells->GetId(i));

  this->GetCellRingNeighbors(pd, fixEdgeCells, 1, 2, ringNeighbors);

  for (int i=0; i<fixEdgeCells->GetNumberOfIds(); i++)
  {
    for (int j=0; j<ringNeighbors[i].size(); j++)
    {
      int cellId = ringNeighbors[i][j];
      int cellValue = pd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(cellId);

      if (cellValue != newCellValue && cellValue != polyRegion.IndexCluster)
        pd->GetCellData()->GetArray(this->GroupIdsArrayName)->SetTuple1(cellId, newCellValue);
    }
  }



  return SV_OK;
}


// ----------------------
// FixOffsetTrifurcation
// ----------------------
int vtkSVSurfaceCenterlineGrouper::FixOffsetTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd, std::string arrayName,
                                                const Region region, const Region polyRegion, std::vector<int> allEdges,
                                                std::vector<int> badEdges, vtkIdList *critPts)
{
  vtkDebugMacro("--FIX OFFSET TRIFURCATION--");
  int patchValue = region.IndexCluster;
  vtkDebugMacro("CLUSTER " << patchValue);
  int numEdges = region.BoundaryEdges.size();
  vtkDebugMacro("NUM EDGES " << numEdges);

  vtkDebugMacro("NUMBER OF ALL EDGES: " << allEdges.size());
  vtkDebugMacro("ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    vtkDebugMacro(" " <<  allEdges[j] << " ");
  vtkDebugMacro("\n");
  vtkDebugMacro("NUMBER OF BAD EDGES: " << badEdges.size());
  vtkDebugMacro("BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    vtkDebugMacro(" " << badEdges[j] << " ");
  vtkDebugMacro("\n");

  std::vector<std::vector<int> > polyConnectedCornerPts;
  this->GetConnectedEdges(polyRegion.BoundaryEdges, polyConnectedCornerPts);
  vtkNew(vtkIdList, polyEdgeGroups);

  int edgeCount = 0;
  std::vector<int> polyEdges;
  for (int i=0; i<polyConnectedCornerPts.size(); i++)
  {
    for (int j=0; j<polyConnectedCornerPts[i].size(); j++)
    {
      int badCornerPtId = region.CornerPoints[badEdges[0]];
      int polyCornerPtId = polyConnectedCornerPts[i][j];

      vtkNew(vtkIdList, badCornerPtList);
      vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, badCornerPtId, badCornerPtList);

      vtkDebugMacro("BAD: ");
      for (int l=0; l<badCornerPtList->GetNumberOfIds(); l++)
        vtkDebugMacro(" " <<  badCornerPtList->GetId(l) << " ");
      vtkDebugMacro("\n");

      vtkNew(vtkIdList, polyCornerPtList);
      vtkSVGeneralUtils::GetPointCellsValues(polyPd, arrayName, polyCornerPtId, polyCornerPtList);

      vtkDebugMacro("POLY: ");
      for (int l=0; l<polyCornerPtList->GetNumberOfIds(); l++)
        vtkDebugMacro(" " <<  polyCornerPtList->GetId(l) << " ");
      vtkDebugMacro("\n");

      vtkNew(vtkIdList, intersectList);
      intersectList->DeepCopy(polyCornerPtList);

      intersectList->IntersectWith(badCornerPtList);
      if (intersectList->GetNumberOfIds() > 1)
      {
        polyEdges.push_back(edgeCount);
      }
      edgeCount++;
    }
  }

  if (polyEdges.size() == 0)
  {
    vtkErrorMacro("COULD NOT FIND POLYCUBE SET OF EDGES MATCHING BAD EDGES");
    return SV_ERROR;
  }

  vtkNew(vtkIdList, polyTouchGroups);
  for (int k=0; k<polyEdges.size(); k++)
  {
    int polyEdge = polyEdges[k];
    for (int i=0; i<polyRegion.BoundaryEdges[polyEdge].size(); i++)
    {
      if (i > 0 && i < polyRegion.BoundaryEdges[polyEdge].size() - 1)
      {
        vtkNew(vtkIdList, tmpList);
        vtkSVGeneralUtils::GetPointCellsValues(polyPd, arrayName, polyRegion.BoundaryEdges[polyEdge][i], tmpList);

        for (int j=0; j<tmpList->GetNumberOfIds(); j++)
          polyTouchGroups->InsertUniqueId(tmpList->GetId(j));
      }
    }
  }

  std::vector<int> fixEdges;
  for (int j=0; j<badEdges.size(); j++)
  {
    int cornerPtId0 = region.CornerPoints[badEdges[j]];
    int cornerPtId1 = region.CornerPoints[badEdges[(j+1)%badEdges.size()]];

    vtkNew(vtkIdList, ptId0List);
    vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, cornerPtId0, ptId0List);

    vtkNew(vtkIdList, intersectList);
    intersectList->DeepCopy(ptId0List);

    intersectList->IntersectWith(polyTouchGroups);

    if (!(polyTouchGroups->GetNumberOfIds() == intersectList->GetNumberOfIds() &&
        ptId0List->GetNumberOfIds() == intersectList->GetNumberOfIds()))
    {
      vtkNew(vtkIdList, ptId1List);
      vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, cornerPtId1, ptId1List);

      vtkNew(vtkIdList, intersectList);
      intersectList->DeepCopy(ptId1List);

      intersectList->IntersectWith(polyTouchGroups);
      if (!(polyTouchGroups->GetNumberOfIds() == intersectList->GetNumberOfIds() &&
          ptId1List->GetNumberOfIds() == intersectList->GetNumberOfIds()))
      {
        fixEdges.push_back(badEdges[j]);
      }
    }
  }

  vtkDebugMacro("NUMBER OF FIX EDGES: " << fixEdges.size());
  vtkDebugMacro("FIX EDGES: ");
  for (int j=0; j<fixEdges.size(); j++)
    vtkDebugMacro(" " <<  fixEdges[j] << " ");
  vtkDebugMacro("\n");

  this->FixEdges(pd, origPd, arrayName, region, allEdges, fixEdges, critPts);

  return SV_OK;
}

// ----------------------
// FixPlanarTrifurcation
// ----------------------
int vtkSVSurfaceCenterlineGrouper::FixPlanarTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                                                const Region region, std::vector<int> allEdges,
                                                std::vector<int> badEdges, vtkIdList *critPts)
{

  vtkDebugMacro("--FIX PLANAR TRIFURCATION--");
  int patchValue = region.IndexCluster;
  vtkDebugMacro("CLUSTER " << patchValue);
  int numEdges = region.BoundaryEdges.size();
  vtkDebugMacro("NUM EDGES " << numEdges);

  vtkDebugMacro("NUMBER OF ALL EDGES: " << allEdges.size());
  vtkDebugMacro("ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    vtkDebugMacro(" " <<  allEdges[j] << " ");
  vtkDebugMacro("\n");
  vtkDebugMacro("NUMBER OF BAD EDGES: " << badEdges.size());
  vtkDebugMacro("BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    vtkDebugMacro(" " <<  badEdges[j] << " ");
  vtkDebugMacro("\n");

  std::vector<int> fixEdges;
  for (int j=0; j<badEdges.size(); j++)
  {
    int cornerPtId0 = region.CornerPoints[badEdges[j]];
    int cornerPtId1 = region.CornerPoints[badEdges[(j+1)%badEdges.size()]];

    vtkNew(vtkIdList, ptId0List);
    vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, cornerPtId0, ptId0List);

    vtkNew(vtkIdList, ptId1List);
    vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, cornerPtId1, ptId1List);

    int numIds = ptId0List->GetNumberOfIds();
    ptId0List->IntersectWith(ptId1List);
    if (ptId0List->GetNumberOfIds() != numIds)
      fixEdges.push_back(badEdges[j]);
  }

  vtkDebugMacro("NUMBER OF FIX EDGES: " << fixEdges.size());
  vtkDebugMacro("FIX EDGES: ");
  for (int j=0; j<fixEdges.size(); j++)
    vtkDebugMacro(" " <<  fixEdges[j] << " ");
  vtkDebugMacro("\n");

  this->FixEdges(pd, origPd, arrayName, region, allEdges, fixEdges, critPts);

  return SV_OK;
}

// ----------------------
// FixPerpenTrifurcation
// ----------------------
int vtkSVSurfaceCenterlineGrouper::FixPerpenTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                                                const Region region, std::vector<int> allEdges,
                                                std::vector<int> badEdges, vtkIdList *critPts)
{
  vtkDebugMacro("--FIX PERPENDICULAR TRIFURCATION--");
  int patchValue = region.IndexCluster;
  vtkDebugMacro("CLUSTER " << patchValue);
  int numEdges = region.BoundaryEdges.size();
  vtkDebugMacro("NUM EDGES " << numEdges);

  vtkDebugMacro("NUMBER OF ALL EDGES: " << allEdges.size());
  vtkDebugMacro("ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    vtkDebugMacro(" " << allEdges[j] << " ");
  vtkDebugMacro("\n");
  vtkDebugMacro("NUMBER OF BAD EDGES: " << badEdges.size());
  vtkDebugMacro("BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    vtkDebugMacro(" " << badEdges[j] << " ");
  vtkDebugMacro("\n");

  std::vector<int> fixEdges;
  if ((badEdges[0] + 1) == badEdges[1])
    fixEdges.push_back(badEdges[0]);
  else
    fixEdges.push_back(badEdges[1]);

  vtkDebugMacro("NUMBER OF FIX EDGES: " << fixEdges.size());
  vtkDebugMacro("FIX EDGES: ");
  for (int j=0; j<fixEdges.size(); j++)
    vtkDebugMacro(" " <<  fixEdges[j] << " ");
  vtkDebugMacro("\n");

  this->FixEdges(pd, origPd, arrayName, region, allEdges, fixEdges, critPts);

  return SV_OK;
}

// ----------------------
// FixEdges
// ----------------------
int vtkSVSurfaceCenterlineGrouper::FixEdges(vtkPolyData *pd, vtkPolyData *origPd,
                                   std::string arrayName,
                                   const Region region, std::vector<int> allEdges,
                                   std::vector<int> fixEdges, vtkIdList *critPts)
{
  int patchValue = region.IndexCluster;

  vtkNew(vtkPolyData, branchPd);
  vtkSVGeneralUtils::ThresholdPd(origPd, patchValue, patchValue, 1, arrayName, branchPd);
  branchPd->BuildLinks();

  for (int j=0; j<fixEdges.size(); j++)
  {
    int badEdgeId = fixEdges[j];
    int edgeSize  = region.BoundaryEdges[badEdgeId].size();

    int startPtId = region.BoundaryEdges[badEdgeId][0];
    int finalPtId = region.BoundaryEdges[badEdgeId][edgeSize-1];

    int halfSize;
    if (edgeSize%2 == 0)
    {
      int testId = region.BoundaryEdges[badEdgeId][edgeSize/2-1];
      vtkDebugMacro("START ID: " << startPtId);
      vtkDebugMacro("FINAL ID: " << finalPtId);
      vtkDebugMacro("TEST ID: " << testId);
      if (critPts->IsId(testId) == -1)
        halfSize = edgeSize/2;
      else
        halfSize = edgeSize/2-1;
    }
    else
    {
      halfSize = floor(edgeSize/2.);
    }

    vtkDebugMacro("HALF SIZE IS: " << halfSize);
    int halfId = region.BoundaryEdges[badEdgeId][halfSize];
    critPts->InsertUniqueId(halfId);
    vtkDebugMacro("HALF ID: " << halfId);

    for (int k=0; k<allEdges.size(); k++)
    {
      if (allEdges[k] != badEdgeId)
      {
        int innerHalfSize = halfSize;
        int allEdgeSize = region.BoundaryEdges[allEdges[k]].size();

        vtkDebugMacro("ALL EDGE SIZE: " << allEdgeSize);
        vtkDebugMacro("HALF SIZE: " << innerHalfSize);
        if (innerHalfSize > allEdgeSize/2)
        {
          innerHalfSize = allEdgeSize/3;
        }

        int stopId       = -1;
        int edgeCell     = -1;
        int newCellValue = -1;

        if (edgeSize == 2)
        {
          vtkNew(vtkIdList, startValues);
          vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, startPtId, startValues);
          vtkDebugMacro("WHAT ARE START VALS: ");
          for (int f=0; f<startValues->GetNumberOfIds(); f++)
            vtkDebugMacro(" " << startValues->GetId(f) << " ");
          vtkDebugMacro("\n");

          vtkNew(vtkIdList, finalValues);
          vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, finalPtId, finalValues);
          vtkDebugMacro("WHAT ARE FINAL VALS: ");
          for (int f=0; f<finalValues->GetNumberOfIds(); f++)
            vtkDebugMacro(" " << finalValues->GetId(f) << " ");
          vtkDebugMacro("\n");

          if (startValues->GetNumberOfIds() == 4 || finalValues->GetNumberOfIds() == 4)
          {
            continue;
          }

          vtkNew(vtkIdList, centerCells);
          vtkNew(vtkIdList, centerValues);
          pd->GetCellEdgeNeighbors(-1, startPtId, finalPtId, centerCells);
          for (int l=0; l<centerCells->GetNumberOfIds(); l++)
          {
            int tmpCellId = centerCells->GetId(l);
            int edgeCellValue = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCellId);
            centerValues->InsertUniqueId(edgeCellValue);
          }

          vtkDebugMacro("WHAT ARE CENTER VALS: ");
          for (int f=0; f<centerValues->GetNumberOfIds(); f++)
            vtkDebugMacro(" " << centerValues->GetId(f) << " ");
          vtkDebugMacro("\n");

          vtkNew(vtkIdList, startCellIds);
          pd->GetPointCells(startPtId, startCellIds);

          int startCellValue = -1;
          int startCellCount = 0;
          for (int l=0; l<startCellIds->GetNumberOfIds(); l++)
          {
            int tmpCellId = startCellIds->GetId(l);
            int edgeCellValue = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCellId);
            if (edgeCellValue == patchValue)
            {
              startCellCount++;
            }

            vtkDebugMacro("this start CEL VAL: " << edgeCellValue);
            if (centerValues->IsId(edgeCellValue) == -1)
            {
              startCellValue = edgeCellValue;
            }
          }
          vtkDebugMacro("FOUND START VALUE: " << startCellValue);

          vtkNew(vtkIdList, finalCellIds);
          pd->GetPointCells(finalPtId, finalCellIds);

          int finalCellValue = -1;
          int finalCellCount = 0;
          for (int l=0; l<finalCellIds->GetNumberOfIds(); l++)
          {
            int tmpCellId = finalCellIds->GetId(l);
            int edgeCellValue = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCellId);
            if (edgeCellValue == patchValue)
            {
              finalCellCount++;
            }

            vtkDebugMacro("this final CEL VAL: " << edgeCellValue);
            if (centerValues->IsId(edgeCellValue) == -1)
            {
              finalCellValue = edgeCellValue;
            }
          }
          vtkDebugMacro("FOUND FINAL VALUE: " << finalCellValue);

          if (startCellCount == 1)
          {
            for (int l=0; l<startCellIds->GetNumberOfIds(); l++)
            {
              int tmpCellId = startCellIds->GetId(l);
              int edgeCellValue = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCellId);
              if (edgeCellValue == patchValue)
              {
                int branchCellId  = branchPd->GetCellData()->GetArray("TmpInternalIds")->LookupValue(tmpCellId);
                branchPd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(branchCellId, startCellValue);
                pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(tmpCellId, startCellValue);
              }
            }
          }
          else
          {
            for (int l=0; l<finalCellIds->GetNumberOfIds(); l++)
            {
              int tmpCellId = finalCellIds->GetId(l);
              int edgeCellValue = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCellId);
              if (edgeCellValue == patchValue)
              {
                int branchCellId  = branchPd->GetCellData()->GetArray("TmpInternalIds")->LookupValue(tmpCellId);
                branchPd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(branchCellId, finalCellValue);
                pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(tmpCellId, finalCellValue);
              }
            }
          }
        }
        else
        {
          if (region.BoundaryEdges[allEdges[k]][0] == finalPtId)
          {
            vtkDebugMacro("ONER");
            stopId = region.BoundaryEdges[allEdges[k]][innerHalfSize];

            vtkNew(vtkIdList, halfValues);
            vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, halfValues);
            vtkDebugMacro("WHAT ARE HALF VALS: ");
            for (int f=0; f<halfValues->GetNumberOfIds(); f++)
              vtkDebugMacro(" " <<  halfValues->GetId(f) << " ");
            vtkDebugMacro("\n");

            vtkNew(vtkIdList, tmpCell);
            origPd->GetPointCells(finalPtId, tmpCell);
            for (int l=0; l<tmpCell->GetNumberOfIds(); l++)
            {
              int edgeCellValue = origPd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCell->GetId(l));
              if (edgeCellValue == patchValue)
                edgeCell = tmpCell->GetId(l);

              if (halfValues->IsId(edgeCellValue) == -1)
                newCellValue = edgeCellValue;
            }
          }

          if (region.BoundaryEdges[allEdges[k]][allEdgeSize-1] == startPtId)
          {
            vtkDebugMacro("TWOER");
            stopId = region.BoundaryEdges[allEdges[k]][allEdgeSize-innerHalfSize-1];

            vtkNew(vtkIdList, halfValues);
            vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, halfValues);
            vtkDebugMacro("WHAT ARE HALF VALS: ");
            for (int f=0; f<halfValues->GetNumberOfIds(); f++)
              vtkDebugMacro(" " <<  halfValues->GetId(f) << " ");
            vtkDebugMacro("\n");

            vtkNew(vtkIdList, tmpCell);
            origPd->GetPointCells(startPtId, tmpCell);
            for (int l=0; l<tmpCell->GetNumberOfIds(); l++)
            {
              int edgeCellValue = origPd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCell->GetId(l));
              if (edgeCellValue == patchValue)
                edgeCell = tmpCell->GetId(l);

              if (halfValues->IsId(edgeCellValue) == -1)
                newCellValue = edgeCellValue;
            }
          }
        }

        if (stopId != -1 && edgeSize != 2)
        {
          if (edgeCell == -1)
          {
            vtkErrorMacro("CANNOT HAVE UNDEFINED CELL TO START PAINTING");
            return SV_ERROR;
          }
          if (newCellValue == -1)
          {
            vtkErrorMacro("CANNOT HAVE UNDEFINED CELL VALUE TO PAINT WITH");
            return SV_ERROR;
          }

          vtkDebugMacro("PLANNING PATH FROM " << halfId << " to " << stopId);
          vtkDebugMacro("STARTING USING EDGE CELL " <<  edgeCell << " AND PAINTING WITH " << newCellValue);

          int startId  = branchPd->GetPointData()->GetArray("TmpInternalIds")->LookupValue(halfId);
          int finalId  = branchPd->GetPointData()->GetArray("TmpInternalIds")->LookupValue(stopId);
          edgeCell     = branchPd->GetCellData()->GetArray("TmpInternalIds")->LookupValue(edgeCell);

          vtkNew(vtkSVFindGeodesicPath, finder);
          finder->SetInputData(branchPd);
          finder->SetStartPtId(startId);
          finder->SetEndPtId(finalId);
          finder->SetDijkstraArrayName("DijkstraDistance");
          finder->SetRepelCloseBoundaryPoints(1);
          finder->Update();

          vtkNew(vtkIdList, tmpIds);
          tmpIds = finder->GetPathIds();
          int numToAdd = tmpIds->GetNumberOfIds();
          vtkDebugMacro("NEW POINTS:              ");
          for (int l=0; l<numToAdd; l++)
            vtkDebugMacro(" " << tmpIds->GetId(l) << " ");
          vtkDebugMacro("\n");

          int count = 1;
          std::vector<int> tempCells;
          tempCells.push_back(edgeCell);

          for (int l=0; l<count; l++)
          {
            int tmpCellId = tempCells[l];
            vtkDebugMacro("DOING CELL: " << tmpCellId);
            branchPd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(tmpCellId, newCellValue);
            int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(tmpCellId);
            pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(realCellId, newCellValue);


            vtkIdType npts, *pts;
            branchPd->GetCellPoints(tmpCellId, npts, pts);
            for (int l=0; l<npts; l++)
            {
              int ptId0 = pts[l];
              int ptId1 = pts[(l+1)%npts];

              int freeEdge =  0;
              int patchEdge = 0;
              int newEdge =   0;

              vtkNew(vtkIdList, cellEdgeNeighbors);
              branchPd->GetCellEdgeNeighbors(tmpCellId, ptId0, ptId1, cellEdgeNeighbors);

              if (cellEdgeNeighbors->GetNumberOfIds() == 0)
                freeEdge = 1;
              else
              {
                int testCellId = cellEdgeNeighbors->GetId(0);
                int cellValue = branchPd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(testCellId);
                if (cellValue != patchValue)
                  patchEdge = 1;

                if (tmpIds->IsId(ptId0) != -1 && tmpIds->IsId(ptId1) != -1)
                  newEdge = 1;
              }


              if (!freeEdge && !patchEdge && !newEdge)
              {
                int nextCellId = cellEdgeNeighbors->GetId(0);
                tempCells.push_back(nextCellId);
                count++;
              }
            }
          }
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// MatchSurfaceToPolycube
// ----------------------
int vtkSVSurfaceCenterlineGrouper::MatchSurfaceToPolycube()
{
  std::vector<Region> surfaceRegions;
  if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, surfaceRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get group regions");
    return SV_ERROR;
  }

  vtkNew(vtkIntArray, newSlicePointsArray);
  newSlicePointsArray->SetNumberOfTuples(this->WorkPd->GetNumberOfPoints());
  newSlicePointsArray->FillComponent(0, -1);
  newSlicePointsArray->SetName("SlicePoints");
  this->WorkPd->GetPointData()->AddArray(newSlicePointsArray);

  vtkNew(vtkIntArray, polySlicePointsArray);
  polySlicePointsArray->SetNumberOfTuples(this->PolycubePd->GetNumberOfPoints());
  polySlicePointsArray->FillComponent(0, -1);
  polySlicePointsArray->SetName("SlicePoints");
  this->PolycubePd->GetPointData()->AddArray(polySlicePointsArray);

  std::vector<Region> polycubeRegions;
  if (this->GetRegions(this->PolycubePd, this->GroupIdsArrayName, polycubeRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get group regions");
    return SV_ERROR;
  }

  for (int i=0; i<surfaceRegions.size(); i++)
  {
    int numEdges = surfaceRegions[i].BoundaryEdges.size();

    int polyRegionId;
    for (int k=0; k<polycubeRegions.size(); k++)
    {
      if (polycubeRegions[k].IndexCluster == surfaceRegions[i].IndexCluster)
        polyRegionId = k;
    }

    for (int j=0; j<numEdges; j++)
    {
      int edgeSize = surfaceRegions[i].BoundaryEdges[j].size();

      int ptId0 = surfaceRegions[i].BoundaryEdges[j][0];
      vtkNew(vtkIdList, ptId0List);
      vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, ptId0, ptId0List);

      int ptIdN = surfaceRegions[i].BoundaryEdges[j][edgeSize-1];
      vtkNew(vtkIdList, ptIdNList);
      vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, ptIdN, ptIdNList);

      vtkDebugMacro("PT ID 0: " << ptId0);
      vtkDebugMacro("IDS 0: ");
      for (int f=0; f<ptId0List->GetNumberOfIds(); f++)
        vtkDebugMacro(" " << ptId0List->GetId(f) << " ");
      vtkDebugMacro("\n");
      vtkDebugMacro("PT ID N: " << ptIdN);
      vtkDebugMacro("IDS N: ");
      for (int f=0; f<ptIdNList->GetNumberOfIds(); f++)
        vtkDebugMacro(" " << ptIdNList->GetId(f) << " ");
      vtkDebugMacro("\n");
      vtkDebugMacro("\n");

      vtkNew(vtkIdList, intersectList);
      intersectList->DeepCopy(ptId0List);

      intersectList->IntersectWith(ptIdNList);

      std::vector<int> newSlicePoints;
      if (intersectList->GetNumberOfIds() >= 3)
      {
        // Traditional between sides of groups
        vtkSVSurfaceCenterlineGrouper::SplitBoundary(this->WorkPd, surfaceRegions[i].BoundaryEdges[j], 3, surfaceRegions[i].IndexCluster,
                                            newSlicePoints);

      }
      else if (intersectList->GetNumberOfIds() == 2)
      {
        // Between center of groups, need to do special
        std::vector<int> newSlicePoints;
        vtkSVSurfaceCenterlineGrouper::SplitBoundary(this->WorkPd, surfaceRegions[i].BoundaryEdges[j], 3, surfaceRegions[i].IndexCluster,
                                            newSlicePoints);
        //vtkSVSurfaceCenterlineGrouper::SplitBoundary(this->WorkPd, surfaceRegions[i].BoundaryEdges[j], 2, surfaceRegions[i].IndexCluster,
        //                                    newSlicePoints);
      }
      else
      {
        vtkErrorMacro("Not sure where this case should happen, not implemented");
        return SV_ERROR;
      }

      for (int k=0; k<newSlicePoints.size(); k++)
      {
        int pointId = newSlicePoints[k];
        vtkDebugMacro("TRYING TO FIND MATCHER FOR " << pointId);

        vtkNew(vtkIdList, surfaceSlicePtVals);
        vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, pointId, surfaceSlicePtVals);

        vtkDebugMacro("POINT CELL VALUES ARE " << surfaceSlicePtVals->GetId(0) << " " << surfaceSlicePtVals->GetId(1));

        // Now find in the polycube
        int edgeDone = 0;
        int numPolyEdges = polycubeRegions[polyRegionId].BoundaryEdges.size();
        for (int l=0; l<numPolyEdges; l++)
        {
          int polyEdgeSize = polycubeRegions[polyRegionId].BoundaryEdges[l].size();

          int polyPtId0 = polycubeRegions[polyRegionId].BoundaryEdges[l][0];
          vtkNew(vtkIdList, polyPtId0List);
          vtkSVGeneralUtils::GetPointCellsValues(this->PolycubePd, this->GroupIdsArrayName, polyPtId0, polyPtId0List);

          int polyPtIdN = polycubeRegions[polyRegionId].BoundaryEdges[l][polyEdgeSize-1];
          vtkNew(vtkIdList, polyPtIdNList);
          vtkSVGeneralUtils::GetPointCellsValues(this->PolycubePd, this->GroupIdsArrayName, polyPtIdN, polyPtIdNList);

          vtkDebugMacro("POLY PT ID 0: " << polyPtId0);
          vtkDebugMacro("POLY IDS 0: ");
          for (int f=0; f<polyPtId0List->GetNumberOfIds(); f++)
            vtkDebugMacro(" " <<  polyPtId0List->GetId(f) << " ");
          vtkDebugMacro("\n");
          vtkDebugMacro("POLY PT ID N: " << polyPtIdN);
          vtkDebugMacro("POLY IDS N: ");
          for (int f=0; f<polyPtIdNList->GetNumberOfIds(); f++)
            vtkDebugMacro(" " <<  polyPtIdNList->GetId(f) << " ");
          vtkDebugMacro("\n");
          vtkDebugMacro("\n");

          vtkNew(vtkIdList, checkList0);
          checkList0->DeepCopy(polyPtId0List);

          checkList0->IntersectWith(ptId0List);

          vtkNew(vtkIdList, checkList1);
          checkList1->DeepCopy(polyPtIdNList);

          checkList1->IntersectWith(ptIdNList);

          if (checkList0->GetNumberOfIds() == ptId0List->GetNumberOfIds() &&
              checkList0->GetNumberOfIds() == polyPtId0List->GetNumberOfIds() &&
              checkList1->GetNumberOfIds() == ptIdNList->GetNumberOfIds() &&
              checkList1->GetNumberOfIds() == polyPtIdNList->GetNumberOfIds())
          {
            vtkDebugMacro("OKAY, THIS IS MATCHING END POINTS");
            vtkDebugMacro("SURFACE PTS: " << ptId0 << " " <<  ptIdN << " POLY PTS: " <<  polyPtId0 << " " << polyPtIdN);

            for (int m=0; m<polyEdgeSize; m++)
            {
              int edgePtId = polycubeRegions[polyRegionId].BoundaryEdges[l][m];

              vtkNew(vtkIdList, polyPatchPtVals);
              vtkSVGeneralUtils::GetPointCellsValues(this->PolycubePd, "PatchIds", edgePtId, polyPatchPtVals);

              if (polyPatchPtVals->GetNumberOfIds() > 2)
              {
                vtkNew(vtkIdList, polyGroupPtVals);
                vtkSVGeneralUtils::GetPointCellsValues(this->PolycubePd, this->GroupIdsArrayName, edgePtId, polyGroupPtVals);

                vtkNew(vtkIdList, valueCheckList);
                valueCheckList->DeepCopy(polyGroupPtVals);

                valueCheckList->IntersectWith(surfaceSlicePtVals);

                if (valueCheckList->GetNumberOfIds() == surfaceSlicePtVals->GetNumberOfIds() &&
                    valueCheckList->GetNumberOfIds() == polyGroupPtVals->GetNumberOfIds())
                {
                  vtkDebugMacro("WE FOUND OUR MATCHING POINT!");
                  vtkDebugMacro("SURFACE PT: " <<  pointId << " POLY PT: " << edgePtId);
                  int currValue = this->PolycubePd->GetPointData()->GetArray("SlicePoints")->GetTuple1(edgePtId);
                  if (currValue != -1)
                  {
                    vtkDebugMacro("ALREADY SET, MAKE SURE NEW POINT " << pointId << " MATCHES " << currValue);
                  }
                  else
                  {
                    this->PolycubePd->GetPointData()->GetArray("SlicePoints")->SetTuple1(edgePtId, pointId);
                    this->PolycubePd->GetPointData()->GetArray("SlicePoints")->SetTuple1(polyPtId0, ptId0);
                    this->PolycubePd->GetPointData()->GetArray("SlicePoints")->SetTuple1(polyPtIdN, ptIdN);
                  }
                  k++;
                  if (k == newSlicePoints.size())
                  {
                    vtkDebugMacro("EDGE DONE");
                    vtkDebugMacro("\n");
                    edgeDone = 1;
                    break;
                  }
                  else
                  {
                    pointId = newSlicePoints[k];
                    vtkNew(vtkIdList, surfaceSlicePtVals);
                    vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, pointId, surfaceSlicePtVals);
                  }
                }
              }
            }
          }
          if (edgeDone)
            break;
        }
        if (edgeDone)
          break;
        else
          vtkErrorMacro("DIDNT FIND A MATCHING PC POINT FOR SLICE POINT " << pointId);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// SplitBoundary
// ----------------------
int vtkSVSurfaceCenterlineGrouper::SplitBoundary(vtkPolyData *pd,
                                         std::vector<int> boundary,
                                         int numDivs,
                                         int groupId,
                                         std::vector<int> &newSlicePoints)
{
  vtkIntArray *slicePoints = vtkIntArray::SafeDownCast(pd->GetPointData()->GetArray("SlicePoints"));
  double fullLength=0.0;
  for (int k=0; k<boundary.size()-1; k++)
  {
    int ptId0 = boundary[k];
    int ptId1 = boundary[k+1];

    double pt0[3], pt1[3];
    pd->GetPoint(ptId0, pt0);
    pd->GetPoint(ptId1, pt1);

    fullLength += vtkSVMathUtils::Distance(pt0, pt1);
  }

  double divLength =  fullLength/numDivs;
  double currLength = 0.0;
  int currDiv = 1;
  for (int k=0; k<boundary.size()-1; k++)
  {
    int ptId0 = boundary[k];
    int ptId1 = boundary[k+1];

    double pt0[3], pt1[3];
    pd->GetPoint(ptId0, pt0);
    pd->GetPoint(ptId1, pt1);

    currLength += vtkSVMathUtils::Distance(pt0, pt1);

    if (currLength > currDiv * divLength)
    {
      currDiv++;

      if (slicePoints->GetTuple1(ptId0) == -1)
      {
        slicePoints->SetTuple1(ptId1, 1);
        newSlicePoints.push_back(ptId1);
      }

      if (currDiv == numDivs)
        break;
    }
  }

  return SV_OK;
}

// ----------------------
// CheckSlicePoints
// ----------------------
int vtkSVSurfaceCenterlineGrouper::CheckSlicePoints()
{

  int numPoints = this->PolycubePd->GetNumberOfPoints();
  int numCells = this->PolycubePd->GetNumberOfPoints();

  vtkDebugMacro("TOTAL NUM OF CELLS: " << this->WorkPd->GetNumberOfCells());
  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, pointCellsValues);
    vtkSVGeneralUtils::GetPointCellsValues(this->PolycubePd, "PatchIds", i, pointCellsValues);

    int numVals = pointCellsValues->GetNumberOfIds();

    int slicePointId = this->PolycubePd->GetPointData()->GetArray("SlicePoints")->GetTuple1(i);

    if (slicePointId != -1)
    {
      vtkNew(vtkIdList, pointCells);
      this->WorkPd->GetPointCells(slicePointId, pointCells);

      int numCells = pointCells->GetNumberOfIds();

      vtkDebugMacro("NUMBER OF CONNECTING PATCHES: " << numVals);
      vtkDebugMacro("VALENCE OF SLICE POINT: " << numCells);

      if (numVals >= (1./2)*numCells)
      {
        // Lets split these cells
        vtkDebugMacro("SPLITTING CELLS AROUND POINT: " << slicePointId);
        this->SplitCellsAroundPoint(this->WorkPd, slicePointId);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// SplitCellsAroundPoint
// ----------------------
int vtkSVSurfaceCenterlineGrouper::SplitCellsAroundPoint(vtkPolyData *pd, int ptId)
{
  vtkNew(vtkIdList, pointCells);
  pd->GetPointCells(ptId, pointCells);

  int numSplitCells = pointCells->GetNumberOfIds();

  vtkDebugMacro("SPLITTING " <<  numSplitCells << " CELLS");

  // Because of poor dynamic editting of data in vtk, we need to
  // create a new set of cells
  int numCurrentCells = pd->GetNumberOfCells();
  int numNewCells     = numCurrentCells + 2*numSplitCells;

  vtkNew(vtkCellArray, newCells);
  newCells->Allocate(numNewCells);
  vtkNew(vtkIdList, cellPtIds);

  for (int i=0; i<numCurrentCells; i++)
  {
    pd->GetCellPoints(i, cellPtIds);
    newCells->InsertNextCell(cellPtIds);
  }

  std::vector<std::vector<int> > splitCells;
  for (int i=0; i<numSplitCells; i++)
  {
    int cellId = pointCells->GetId(i);
    vtkDebugMacro("SPLITIING CELL: " << cellId);

    vtkIdType npts, *pts;
    pd->GetCellPoints(cellId, npts, pts);
    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];

      if (ptId0 != ptId && ptId1 != ptId)
      {
        this->SplitEdge(pd, cellId, ptId0, ptId1, newCells, splitCells);
        break;
      }
    }
  }

  pd->SetPolys(newCells);
  pd->BuildCells();
  pd->BuildLinks();

  for (int i=0; i<splitCells.size(); i++)
  {
    int replaceCellId = splitCells[i][0];
    int oldPtId = splitCells[i][1];
    int newPtId = splitCells[i][2];
    pd->ReplaceCellPoint(replaceCellId, oldPtId, newPtId);
  }

  return SV_OK;
}

// ----------------------
// SplitEdge
// ----------------------
int vtkSVSurfaceCenterlineGrouper::SplitEdge(vtkPolyData *pd, int cellId, int ptId0, int ptId1,
                                    vtkCellArray *newCells, std::vector<std::vector<int> > &splitCells)

{
  // Num pts
  int numCurrentPts   = pd->GetNumberOfPoints();
  int numNewPts       = numCurrentPts + 1;

  // Now do stuff
  vtkNew(vtkIdList, edgeCells);
  pd->GetCellEdgeNeighbors(cellId, ptId0, ptId1, edgeCells);

  if (cellId != -1)
    edgeCells->InsertNextId(cellId);

  int pointAdded = 0;
  int newPointId = numNewPts-1;
  for (int i=0; i<edgeCells->GetNumberOfIds(); i++)
  {
    int splitCellId = edgeCells->GetId(i);

    vtkIdType npts, *pts;
    pd->GetCellPoints(splitCellId, npts, pts);

    for (int j=0; j<npts; j++)
    {
      int splitPtId0 = pts[j];
      int splitPtId1 = pts[(j+1)%npts];

      if ((splitPtId0 == ptId0 && splitPtId1 == ptId1) ||
          (splitPtId1 == ptId0 && splitPtId0 == ptId1))
      {
        int thirdPtId = pts[(j+2)%npts];

        double pt0[3], pt1[3], newPt[3];
        pd->GetPoint(ptId0, pt0);
        pd->GetPoint(ptId1, pt1);

        vtkMath::Add(pt0, pt1, newPt);
        vtkMath::MultiplyScalar(newPt, 1./2);

        if (!pointAdded)
        {
          pd->GetPoints()->InsertNextPoint(newPt);
          pointAdded = 1;

         pd->GetPointData()->CopyData(pd->GetPointData(), ptId0, newPointId);

         for (int k=0; k<pd->GetPointData()->GetNumberOfArrays(); k++)
         {
           double weights[2]; weights[0] = 0.5; weights[1] = 0.5;

           vtkNew(vtkIdList, interpIds);
           interpIds->SetNumberOfIds(2);
           interpIds->SetId(0, ptId0);
           interpIds->SetId(1, ptId1);

           pd->GetPointData()->GetArray(k)->InsertNextTuple(
             pd->GetPointData()->GetArray(k)->GetTuple(ptId0));
           pd->GetPointData()->GetArray(k)->InterpolateTuple(newPointId,
               interpIds, pd->GetPointData()->GetArray(k), weights);
         }
        }

        std::vector<int> splitCellInfo(3);
        splitCellInfo[0] = splitCellId;
        splitCellInfo[1] = ptId1;
        splitCellInfo[2] = newPointId;
        splitCells.push_back(splitCellInfo);

        vtkNew(vtkIdList, newCell);
        newCell->SetNumberOfIds(3);
        newCell->SetId(0, thirdPtId);
        newCell->SetId(1, newPointId);
        newCell->SetId(2, ptId1);

        int newCellId = newCells->InsertNextCell(newCell);

        for (int k=0; k<pd->GetCellData()->GetNumberOfArrays(); k++)
        {
          pd->GetCellData()->GetArray(k)->InsertNextTuple(
            pd->GetCellData()->GetArray(k)->GetTuple(splitCellId));
        }
        pd->GetCellData()->CopyData(pd->GetCellData(), splitCellId, newCellId);
      }
    }
  }

  return SV_OK;
}



// ----------------------
// GetCellRingNeighbors
// ----------------------
int vtkSVSurfaceCenterlineGrouper::GetCellRingNeighbors(vtkPolyData *pd, vtkIdList *cellIds,
                                               int ringNumber,
                                               int totNumberOfRings,
                                               std::vector<std::vector<int> > &neighbors)
{
  // Number of cells
  int numCells = cellIds->GetNumberOfIds();

  for (int i=0; i<numCells; i++)
  {
    // temporary node vec
    std::vector<int> tmpNodes;
    int iSize = neighbors[i].size();

    for (int j=0; j<iSize; j++)
    {
      // Get neighbor cell points
      int neiCellId = neighbors[i][j];
      vtkIdType *pts, npts;
      pd->GetCellPoints(neiCellId, npts, pts);

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

      vtkNew(vtkIdList, pointCellIds);
      pd->GetPointCells(tmpNode, pointCellIds);
      for (int k=0; k<pointCellIds->GetNumberOfIds(); k++)
      {
        int tmpCell = pointCellIds->GetId(k);
        int kSize =   neighbors[i].size();

        int kk=0;
        for (kk=0; kk<kSize; kk++)
        {
          if (tmpCell == neighbors[i][kk])
          {
            break;
          }
        }
        if (kk == kSize)
        {
          neighbors[i].push_back(tmpCell);
        }
      }
    }
  }

  if (ringNumber < totNumberOfRings)
  {
    ringNumber++;
    this->GetCellRingNeighbors(pd, cellIds, ringNumber, totNumberOfRings, neighbors);
  }

  return SV_OK;
}