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

#include "vtkSVGroupsSegmenter.h"
#include "vtkAppendPolyData.h"
#include "vtkExecutive.h"
#include "vtkCellLocator.h"
#include "vtkConnectivityFilter.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPointLocator.h"
#include "vtkCellData.h"
#include "vtkIdFilter.h"
#include "vtkIntArray.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDoubleArray.h"
#include "vtkCleanPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkFeatureEdges.h"
#include "vtkGenericCell.h"
#include "vtkSmartPointer.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVPolyBallLine.h"
#include "vtkMath.h"
#include "vtkMergeCells.h"
#include "vtkSphere.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataNormals.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"
#include "vtkThreshold.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVersion.h"
#include "vtkXMLPolyDataWriter.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVGroupsSegmenter);

// ----------------------
// Constructor
// ----------------------
vtkSVGroupsSegmenter::vtkSVGroupsSegmenter()
{
  this->WorkPd = vtkPolyData::New();
  this->Centerlines = NULL;

  this->CenterlineGroupIdsArrayName = NULL;
  this->CenterlineRadiusArrayName = NULL;
  this->GroupIdsArrayName = NULL;
  this->BlankingArrayName = NULL;
  this->CenterlineGroupIds = NULL;

  this->ClipAllCenterlineGroupIds = 1;
  this->CutoffRadiusFactor = VTK_SV_LARGE_DOUBLE;
  this->ClipValue = 0.0;
  this->UseRadiusInformation = 1;
}

// ----------------------
// Destructor
// ----------------------
vtkSVGroupsSegmenter::~vtkSVGroupsSegmenter()
{
  if (this->WorkPd)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->Centerlines)
  {
    this->Centerlines->Delete();
    this->Centerlines = NULL;
  }

  if (this->CenterlineGroupIds)
  {
    this->CenterlineGroupIds->Delete();
    this->CenterlineGroupIds = NULL;
  }

  if (this->CenterlineGroupIdsArrayName)
  {
    delete [] this->CenterlineGroupIdsArrayName;
    this->CenterlineGroupIdsArrayName = NULL;
  }

  if (this->CenterlineRadiusArrayName)
  {
    delete [] this->CenterlineRadiusArrayName;
    this->CenterlineRadiusArrayName = NULL;
  }

  if (this->GroupIdsArrayName)
  {
    delete [] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
  }

  if (this->BlankingArrayName)
  {
    delete [] this->BlankingArrayName;
    this->BlankingArrayName = NULL;
  }
}

// ----------------------
// RequestData
// ----------------------
int vtkSVGroupsSegmenter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  this->WorkPd->DeepCopy(input);
  // If there is one centerline, we don't need to separate, just return
  // with one group id
  if (this->Centerlines->GetNumberOfCells() == 1)
  {
    output->DeepCopy(input);
    vtkIntArray *groupIdsArray = vtkIntArray::New();
    groupIdsArray->SetName(this->GroupIdsArrayName);
    groupIdsArray->SetNumberOfTuples(output->GetNumberOfCells());
    groupIdsArray->FillComponent(0,0);
    output->GetCellData()->AddArray(groupIdsArray);
    groupIdsArray->Delete();
    return SV_OK;
  }

  // Prep work for filter
  if (this->PrepFilter() != SV_OK)
  {
    vtkErrorMacro("Prep of filter failed");
    output->DeepCopy(input);
    return SV_ERROR;
  }

  // Run the filter
  if (this->RunFilter() != SV_OK)
  {
    vtkErrorMacro("Filter failed");
    output->DeepCopy(input);
    return SV_ERROR;
  }

  output->DeepCopy(this->WorkPd);

  return SV_OK;
}

// ----------------------
// PrepFilter
// ----------------------
int vtkSVGroupsSegmenter::PrepFilter()
{
  if (!this->Centerlines)
  {
    vtkErrorMacro(<< "Centerlines not set.");
    return SV_ERROR;
  }

  if (!this->ClipAllCenterlineGroupIds && !this->CenterlineGroupIds)
  {
    vtkErrorMacro(<< "CenterlineGroupIds not set.");
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

  if (vtkSVGeneralUtils::CheckArrayExists(this->Centerlines, 1, this->CenterlineGroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "CenterlineGroupIdsArray with name specified does not exist");
    return SV_OK;
  }

  if (!this->BlankingArrayName)
  {
    vtkDebugMacro("Blanking Array Name not given, setting to Blanking");
    this->BlankingArrayName = new char[strlen("Blanking") + 1];
    strcpy(this->BlankingArrayName, "Blanking");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->Centerlines, 1, this->BlankingArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "BlankingArrayName with name specified does not exist");
    return SV_ERROR;
  }

  if (!this->CenterlineRadiusArrayName)
  {
    vtkDebugMacro("Centerline radius Array Name not given, setting to MaximumInscribedSphereRadius");
    this->CenterlineRadiusArrayName = new char[strlen("MaximumInscribedSphereRadius") + 1];
    strcpy(this->CenterlineRadiusArrayName, "MaximumInscribedSphereRadius");
  }

  if (!this->Centerlines->GetPointData()->GetArray(this->CenterlineRadiusArrayName))
  {
    vtkErrorMacro(<< "CenterlineRadiusArray with name specified does not exist");
    return SV_ERROR;
  }

  this->CenterlineGraph = new svCenterlineGraph(0, this->Centerlines,
                                                this->GroupIdsArrayName);

  if (this->CenterlineGraph->BuildGraph() != SV_OK)
  {
    vtkErrorMacro("Unable to form graph of centerlines");
    return SV_ERROR;
  }

  if (this->CenterlineGraph->RefineGraphWithLocalCoordinates() != SV_OK)
  {
    vtkErrorMacro("Unable to form graph of centerlines");
    return SV_ERROR;
  }

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/CenterlineGraph.vtp";
  vtkNew(vtkPolyData, graphPd);
  this->CenterlineGraph->GetGraphPolyData(graphPd);
  vtkSVIOUtils::WriteVTPFile(filename, graphPd);

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVGroupsSegmenter::RunFilter()
{
  // Get data arrays
  vtkDataArray *centerlineGroupIdsArray =
    this->Centerlines->GetCellData()->GetArray(this->CenterlineGroupIdsArrayName);
  vtkIntArray *blankingArray =
    vtkIntArray::SafeDownCast(this->Centerlines->GetCellData()->GetArray(this->BlankingArrayName));

  // for each group, compute the clipping array, clip, add group ids array and append.
  vtkNew(vtkSVPolyBallLine, groupTubes);
  groupTubes->SetInput(this->Centerlines);
  groupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  groupTubes->SetUseRadiusInformation(this->UseRadiusInformation);
  groupTubes->ControlEndPointsOff();
  groupTubes->UsePointNormalOn();

  double point[3];
  double groupTubeValue;
  vtkIdType groupId;

  // Get all the group ids if clipping all
  vtkNew(vtkIdList, centerlineGroupIds);
  int i;
  if (this->ClipAllCenterlineGroupIds)
  {
    for (i=0; i<centerlineGroupIdsArray->GetNumberOfTuples(); i++)
    {
      if (blankingArray->GetValue(i) == 1)
        continue;

      centerlineGroupIds->InsertUniqueId(static_cast<vtkIdType>(vtkMath::Round(centerlineGroupIdsArray->GetComponent(i,0))));
    }
  }
  else
    centerlineGroupIds->DeepCopy(this->CenterlineGroupIds);

  // Clipping input
  vtkNew(vtkPolyDataNormals, normaler);
  normaler->SetInputData(this->WorkPd);
  normaler->ComputePointNormalsOn();
  normaler->ComputeCellNormalsOff();
  normaler->SplittingOff();
  normaler->Update();

  vtkNew(vtkPolyData, clippingInput);
  clippingInput->DeepCopy(normaler->GetOutput());
  vtkDataArray *normalsArray =
    clippingInput->GetPointData()->GetArray("Normals");

  // Add array for group cell ids on surface
  vtkNew(vtkIntArray, startGroupIds);
  startGroupIds->SetName(this->GroupIdsArrayName);
  startGroupIds->SetNumberOfTuples(clippingInput->GetNumberOfPoints());
  startGroupIds->FillComponent(0,-1);
  clippingInput->GetPointData()->AddArray(startGroupIds);

  // Appender to append all the clipped together
  vtkNew(vtkAppendPolyData, appendBranches);
  int numberOfPoints = clippingInput->GetNumberOfPoints();

  // Loop through points to evaluate function at each point
  for (int k=0; k<numberOfPoints; k++)
  {
    // Evaluate function at point!
    clippingInput->GetPoint(k,point);
    double pointNormal[3];
    normalsArray->GetTuple(k, pointNormal);
    groupTubes->SetPointNormal(pointNormal);
    groupTubeValue = groupTubes->EvaluateFunction(point);

    // Set to very large value if greater than threshold
    if (groupTubeValue > this->CutoffRadiusFactor * this->CutoffRadiusFactor - 1)
      groupTubeValue = VTK_SV_LARGE_DOUBLE;

    startGroupIds->SetTuple1(k, centerlineGroupIdsArray->GetTuple1(groupTubes->GetLastPolyBallCellId()));
  }

  if (this->PassPointGroupsToCells(clippingInput, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro("Couldn't extend group info to cells");
    return SV_ERROR;
  }

  if (this->CorrectCellBoundaries(clippingInput, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro("Could not correcto boundaries of surface");
    return SV_ERROR;
  }

  // Finalize
  this->WorkPd->DeepCopy(clippingInput);

  return SV_OK;
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVGroupsSegmenter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Clip value: " << this->ClipValue << "\n";
  os << indent << "Cutoff Radius Factor: " << this->CutoffRadiusFactor << "\n";
  os << indent << "Clip all centerline group ids: " << this->ClipAllCenterlineGroupIds << "\n";
  os << indent << "Use radius information: " << this->UseRadiusInformation << "\n";
  if (this->CenterlineGroupIdsArrayName != NULL)
    os << indent << "Centerline group ids name: " << this->CenterlineGroupIdsArrayName << "\n";
  if (this->CenterlineRadiusArrayName != NULL)
    os << indent << "Centerline radius array name: " << this->CenterlineRadiusArrayName << "\n";
  if (this->GroupIdsArrayName != NULL)
    os << indent << "Group ids array name: " << this->GroupIdsArrayName << "\n";
  if (this->BlankingArrayName != NULL)
    os << indent << "Blanking array name: " << this->BlankingArrayName << "\n";
  if (this->CenterlineGroupIds != NULL)
  {
    os << indent << "Target values to clip: "<< "\n";
      os << indent;
    for (int i=0; i<this->CenterlineGroupIds->GetNumberOfIds(); i++)
      os << this->CenterlineGroupIds->GetId(i);
    os << "\n";
  }
}

// ----------------------
// PassPointGroupsToCells
// ----------------------
int vtkSVGroupsSegmenter::PassPointGroupsToCells(vtkPolyData *pd, std::string pointArrayName )
{
  vtkNew(vtkIntArray, cellIds);
  cellIds->SetNumberOfTuples(pd->GetNumberOfCells());
  cellIds->FillComponent(0, -1);
  cellIds->SetName(pointArrayName.c_str());
  pd->GetCellData()->AddArray(cellIds);

  vtkDataArray *pointIds = pd->GetPointData()->GetArray(pointArrayName.c_str());

  // Pass data from points to cells using largest val
  pd->BuildLinks();
  for (int i=0; i<pd->GetNumberOfCells(); i++)
  {

    // Get cell points
    vtkIdType *pts, npts;
    pd->GetCellPoints(i, npts, pts);

    // Get list of all point values
    vtkNew(vtkIdList, pointVals);
    pointVals->SetNumberOfIds(npts);
    for (int j=0; j<npts; j++)
    {
      double pointVal = pointIds->GetTuple1(pts[j]);
      pointVals->SetId(j, pointVal);
    }

    // Find the most occuring
    int mostOccuring, maxCount;
    vtkSVGroupsSegmenter::GetMostOccuringVal(pointVals, mostOccuring, maxCount);

    // Set the cell value based on most occuring point value
    cellIds->SetTuple1(i, mostOccuring);
  }

  return 1;
}

// ----------------------
// CorrectCellBoundaries
// ----------------------
int vtkSVGroupsSegmenter::CorrectCellBoundaries(vtkPolyData *pd, std::string cellArrayName )
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

  // Loop through cells again
  for (int i=0; i<numCells; i++)
  {
    // get direct neighbor value count
    vtkNew(vtkIdList, neiCellIds);

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
        if (tmpIds->GetTuple1(cellEdgeNeighbor) != tmpIds->GetTuple1(i));
          neiCellIds->InsertNextId(cellIds->GetTuple1(cellEdgeNeighbor));
      }
    }

    // If we found a cell surrounded by cells of another val, we can update
    int neiSize = neiCellIds->GetNumberOfIds();
    if (neiSize > 1)
    {
      int maxVal, maxCount;
      vtkSVGroupsSegmenter::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

      cellIds->SetTuple1(i, maxVal);
      tmpIds->SetTuple1(i, maxVal);
    }
  }

  return 1;
}

// ----------------------
// GetMostOccuringVal
// ----------------------
void vtkSVGroupsSegmenter::GetMostOccuringVal(vtkIdList *idList, int &output,
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
