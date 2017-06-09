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
#include "vtkDataSetSurfaceFilter.h"
#include "vtkIdFilter.h"
#include "vtkIntArray.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDoubleArray.h"
#include "vtkCleanPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkFeatureEdges.h"
#include "vtkGenericCell.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkSVCleanUnstructuredGrid.h"
#include "vtkSVEdgeWeightedCVT.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVMathUtils.h"
#include "vtkSVPolyBallLine.h"
#include "vtkSVFindSeparateRegions.h"
#include "vtkSVPlanarMapper.h"
#include "vtkSVPointSetBoundaryMapper.h"
#include "vtkSVMapInterpolator.h"
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
// GlobalCoords
// ----------------------
const double vtkSVGroupsSegmenter::GlobalCoords[3][3] =
  {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
  };

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
  this->GraphPd = vtkPolyData::New();
  this->CenterlinesWorkPd = vtkPolyData::New();
  this->Polycube = vtkUnstructuredGrid::New();
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
  if (this->CenterlinesWorkPd)
  {
    this->CenterlinesWorkPd->Delete();
    this->CenterlinesWorkPd = NULL;
  }
  if (this->Centerlines)
  {
    this->Centerlines->Delete();
    this->Centerlines = NULL;
  }
  if (this->Polycube)
  {
    this->Polycube->Delete();
    this->Polycube = NULL;
  }
  if (this->GraphPd)
  {
    this->GraphPd->Delete();
    this->GraphPd = NULL;
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
  this->CenterlinesWorkPd->DeepCopy(this->Centerlines);
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
  this->CenterlinesWorkPd->DeepCopy(this->Centerlines);

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

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/CenterlineGraph.vtp";
  this->CenterlineGraph->GetGraphPolyData(this->GraphPd);
  vtkSVIOUtils::WriteVTPFile(filename, this->GraphPd);

  this->CenterlinesWorkPd->DeepCopy(this->CenterlineGraph->Lines);
  std::string filename2 = "/Users/adamupdegrove/Desktop/tmp/CenterlineDirs.vtp";
  vtkSVIOUtils::WriteVTPFile(filename2, this->CenterlineGraph->Lines);

  this->CenterlineGraph->GetPolycube(1.0, 1.0, this->Polycube);
  std::string filename3 = "/Users/adamupdegrove/Desktop/tmp/Polycube.vtu";
  vtkSVIOUtils::WriteVTUFile(filename3, this->Polycube);

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVGroupsSegmenter::RunFilter()
{
  // Get data arrays
  vtkDataArray *centerlineGroupIdsArray =
    this->CenterlinesWorkPd->GetCellData()->GetArray(this->CenterlineGroupIdsArrayName);

  // for each group, compute the clipping array, clip, add group ids array and append.
  vtkNew(vtkSVPolyBallLine, groupTubes);
  groupTubes->SetInput(this->CenterlinesWorkPd);
  groupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  groupTubes->SetUseRadiusInformation(this->UseRadiusInformation);
  groupTubes->ControlEndPointsOff();
  groupTubes->UsePointNormalOn();
  groupTubes->UseRadiusWeightingOn();
  groupTubes->UseLocalCoordinatesOn();
  groupTubes->SetLocalCoordinatesArrayName("Local");

  double point[3];
  double groupTubeValue;
  vtkIdType groupId;

  // Clipping input
  vtkNew(vtkPolyDataNormals, normaler);
  normaler->SetInputData(this->WorkPd);
  normaler->ComputePointNormalsOff();
  normaler->ComputeCellNormalsOn();
  normaler->SplittingOff();
  normaler->Update();

  int numberOfCells = this->WorkPd->GetNumberOfCells();

  this->WorkPd->DeepCopy(normaler->GetOutput());
  vtkDataArray *normalsArray =
    this->WorkPd->GetCellData()->GetArray("Normals");
  // Add array for group cell ids on surface
  vtkNew(vtkIntArray, startGroupIds);
  startGroupIds->SetName(this->GroupIdsArrayName);
  startGroupIds->SetNumberOfTuples(numberOfCells);
  startGroupIds->FillComponent(0,-1);
  this->WorkPd->GetCellData()->AddArray(startGroupIds);

  // Add array for new cell normals on surface
  vtkNew(vtkDoubleArray, newCellNormals);
  newCellNormals->SetName("CenterlinesBasedCellNormals");
  newCellNormals->SetNumberOfComponents(3);
  newCellNormals->SetNumberOfTuples(numberOfCells);
  this->WorkPd->GetCellData()->AddArray(newCellNormals);
  this->WorkPd->BuildLinks();

  // Loop through points to evaluate function at each point
  fprintf(stdout,"Computing closest centerline points per cell...\n");
  for (int k=0; k<numberOfCells; k++)
  {
    // Get cell point coords
    double pts[3][3];
    vtkIdType npts, *ptids;
    this->WorkPd->GetCellPoints(k, npts, ptids);
    for (int j=0; j<npts; j++)
      this->WorkPd->GetPoint(ptids[j], pts[j]);

    // Get center
    double center[3];
    vtkTriangle::TriangleCenter(pts[0], pts[1], pts[2], center);

    // Evaluate function at point!
    double cellNormal[3];
    normalsArray->GetTuple(k, cellNormal);
    groupTubes->SetPointNormal(cellNormal);
    groupTubeValue = groupTubes->EvaluateFunction(center);

    // Set to very large value if greater than threshold
    if (groupTubeValue > this->CutoffRadiusFactor * this->CutoffRadiusFactor - 1)
      groupTubeValue = VTK_SV_LARGE_DOUBLE;

    startGroupIds->SetTuple1(k, centerlineGroupIdsArray->GetTuple1(groupTubes->GetLastPolyBallCellId()));

    // Now get last local coords and use rotation matrix to set new normals
    double localX[3], localY[3], localZ[3];
    groupTubes->GetLastLocalCoordX(localX);
    groupTubes->GetLastLocalCoordY(localY);
    groupTubes->GetLastLocalCoordZ(localZ);

    // Compute the rotation from global coordinate system to centerlines
    // local coordinate system
    double rotMat[9];
    this->ComputeRotationMatrix(localX, localY, localZ, rotMat);

    // Apply rotation matrix to the normal to get the new normal
    double newNormal[3];
    for (int j=0; j<3; j++)
    {
      newNormal[j] = rotMat[j*3]*cellNormal[0] +
                     rotMat[(j*3)+1]*cellNormal[1] +
                     rotMat[(j*3)+2]*cellNormal[2];
    }
    newCellNormals->SetTuple(k, newNormal);

  }

  if (this->RunEdgeWeightedCVT(this->WorkPd) != SV_OK)
  {
    vtkErrorMacro("Error in cvt");
    return SV_ERROR;
  }

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
  if (this->CorrectCellBoundaries(this->WorkPd, "PatchIds") != SV_OK)
  {
    vtkErrorMacro("Could not correcto boundaries of surface");
    return SV_ERROR;
  }
  if (this->SmoothBoundaries(this->WorkPd, "PatchIds") != SV_OK)
  {
    vtkErrorMacro("Could not smootho boundaries of surface");
    return SV_ERROR;
  }
  std::vector<Region> patchRegions;
  if (this->GetRegions(this->WorkPd, "PatchIds", patchRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }
  if (this->CurveFitBoundaries(this->WorkPd, "PatchIds", patchRegions) != SV_OK)
  {
    vtkErrorMacro("Could not curve fit boundaries of surface");
    return SV_ERROR;
  }

  // Get final patch ids by adding to group ids
  vtkNew(vtkIdList, groupIds);
  for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
  {
    int groupVal = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(i);
    groupIds->InsertUniqueId(groupVal);
  }
  vtkSortDataArray::Sort(groupIds);
  int numGroups = groupIds->GetNumberOfIds();
  vtkNew(vtkIdList, addVals);
  addVals->SetNumberOfIds(numGroups);
  for (int i=0; i<numGroups; i++)
    addVals->SetId(i, 6*i);

  vtkNew(vtkIdList, patchVals);
  for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
  {
    int patchVal = this->WorkPd->GetCellData()->GetArray("PatchIds")->GetTuple1(i);
    int groupVal = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(i);
    int newVal = patchVal + (addVals->GetId(groupIds->IsId(groupVal)));
    this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(i, newVal);
    patchVals->InsertUniqueId(newVal);
  }

  // CHECK AND FIX UP OF REGIONS IF BADD!!!
  if (this->FixRegions(this->WorkPd, "PatchIds") != SV_OK)
  {
    fprintf(stderr,"Couldn't fix regions\n");
    return SV_ERROR;
  }


  // NOW PARAMETERIZE!!, WIILL BE MOVED to vtkSVPolycubeParameterizer
  // TODO: RENAME THIS CLASS TO vtkSVCenterlinesSegmenter

  if (this->Parameterize() != SV_OK)
  {
    fprintf(stderr,"WRONG\n");
    return SV_ERROR;
  }


  return SV_OK;
}

// ----------------------
// RunEdgeWeightedCVT
// ----------------------
int vtkSVGroupsSegmenter::RunEdgeWeightedCVT(vtkPolyData *pd)
{
  // Set up generators
  vtkNew(vtkPoints, generatorsPts);
  generatorsPts->SetNumberOfPoints(6);
  generatorsPts->SetPoint(0, 1.0, 0.0, 0.0);
  generatorsPts->SetPoint(1, 0.0, 1.0, 0.0);
  generatorsPts->SetPoint(2, -1.0, 0.0, 0.0);
  generatorsPts->SetPoint(3, 0.0, -1.0, 0.0);
  generatorsPts->SetPoint(4, 0.0, 0.0, 1.0);
  generatorsPts->SetPoint(5, 0.0, 0.0, -1.0);

  vtkNew(vtkPolyData, generatorsPd);
  generatorsPd->SetPoints(generatorsPts);

  // Run edge weighted cvt
  vtkNew(vtkSVEdgeWeightedCVT, CVT);

  CVT->SetInputData(pd);
  CVT->SetGenerators(generatorsPd);
  CVT->SetNumberOfRings(2);
  CVT->SetThreshold(2);
  CVT->SetEdgeWeight(1.0);
  CVT->SetMaximumNumberOfIterations(1000);
  CVT->SetPatchIdsArrayName("PatchIds");
  CVT->SetCVTDataArrayName("CenterlinesBasedCellNormals");
  CVT->Update();

  pd->DeepCopy(CVT->GetOutput());

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

// ----------------------
// ComputeRotationMatrix
// ----------------------
int vtkSVGroupsSegmenter::ComputeRotationMatrix(const double vx[3],
                                                const double vy[3],
                                                const double vz[3],
                                                double rotMatrix[9])
{
  rotMatrix[0] = vx[0]*vtkSVGroupsSegmenter::GlobalCoords[0][0] +
                 vx[1]*vtkSVGroupsSegmenter::GlobalCoords[0][1] +
                 vx[2]*vtkSVGroupsSegmenter::GlobalCoords[0][2];
  rotMatrix[1] = vx[0]*vtkSVGroupsSegmenter::GlobalCoords[1][0] +
                 vx[1]*vtkSVGroupsSegmenter::GlobalCoords[1][1] +
                 vx[2]*vtkSVGroupsSegmenter::GlobalCoords[1][2];
  rotMatrix[2] = vx[0]*vtkSVGroupsSegmenter::GlobalCoords[2][0] +
                 vx[1]*vtkSVGroupsSegmenter::GlobalCoords[2][1] +
                 vx[2]*vtkSVGroupsSegmenter::GlobalCoords[2][2];

  rotMatrix[3] = vy[0]*vtkSVGroupsSegmenter::GlobalCoords[0][0] +
                 vy[1]*vtkSVGroupsSegmenter::GlobalCoords[0][1] +
                 vy[2]*vtkSVGroupsSegmenter::GlobalCoords[0][2];
  rotMatrix[4] = vy[0]*vtkSVGroupsSegmenter::GlobalCoords[1][0] +
                 vy[1]*vtkSVGroupsSegmenter::GlobalCoords[1][1] +
                 vy[2]*vtkSVGroupsSegmenter::GlobalCoords[1][2];
  rotMatrix[5] = vy[0]*vtkSVGroupsSegmenter::GlobalCoords[2][0] +
                 vy[1]*vtkSVGroupsSegmenter::GlobalCoords[2][1] +
                 vy[2]*vtkSVGroupsSegmenter::GlobalCoords[2][2];

  rotMatrix[6] = vz[0]*vtkSVGroupsSegmenter::GlobalCoords[0][0] +
                 vz[1]*vtkSVGroupsSegmenter::GlobalCoords[0][1] +
                 vz[2]*vtkSVGroupsSegmenter::GlobalCoords[0][2];
  rotMatrix[7] = vz[0]*vtkSVGroupsSegmenter::GlobalCoords[1][0] +
                 vz[1]*vtkSVGroupsSegmenter::GlobalCoords[1][1] +
                 vz[2]*vtkSVGroupsSegmenter::GlobalCoords[1][2];
  rotMatrix[8] = vz[0]*vtkSVGroupsSegmenter::GlobalCoords[2][0] +
                 vz[1]*vtkSVGroupsSegmenter::GlobalCoords[2][1] +
                 vz[2]*vtkSVGroupsSegmenter::GlobalCoords[2][2];

  return SV_OK;
}

// ----------------------
// SmoothBoundaries
// ----------------------
int vtkSVGroupsSegmenter::SmoothBoundaries(vtkPolyData *pd, std::string arrayName)
{
  vtkNew(vtkSVFindSeparateRegions, boundaryFinder);
  boundaryFinder->SetInputData(pd);
  boundaryFinder->SetCellArrayName(arrayName.c_str());
  boundaryFinder->SetOutPointArrayName("SmoothBoundaryPoints");
  boundaryFinder->Update();

  vtkNew(vtkPolyData, boundaryIndicator);
  boundaryIndicator->DeepCopy(boundaryFinder->GetOutput());
  boundaryIndicator->BuildLinks();

  int numPoints = pd->GetNumberOfPoints();

  for (int i=0; i<numPoints; i++)
  {
    int val = boundaryIndicator->GetPointData()->
      GetArray("SmoothBoundaryPoints")->GetTuple1(i);
    if (val == 1)
    {
      vtkNew(vtkIdList, pointCellsValues);
      vtkSVGeneralUtils::GetPointCellsValues(boundaryIndicator, arrayName.c_str(),
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
              for (int j=0; j<2; j++)
              {
                pd->GetCellPoints(cellIds[0][j], npts, pts);
                for (int p=0; p<npts; p++)
                  uniquePoints->InsertUniqueId(pts[p]);
              }
            }
            else if (count[1] == 2)
            {
              for (int j=0; j<2; j++)
              {
                pd->GetCellPoints(cellIds[1][j], npts, pts);
                for (int p=0; p<npts; p++)
                  uniquePoints->InsertUniqueId(pts[p]);
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
int vtkSVGroupsSegmenter::GetRegions(vtkPolyData *pd, std::string arrayName,
                                     std::vector<Region> &allRegions)
{

  int numCells = pd->GetNumberOfCells();
  std::vector<std::vector<int> > tempRegions(numCells);
  std::vector<std::vector<int> > directNeighbors(numCells);
  std::vector<int> numberOfDirectNeighbors(numCells);

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

  allRegions.resize(numberOfRegions);

  for (int i=0; i<numberOfRegions; i++)
  {
    allRegions[i].Index = i;
    allRegions[i].NumberOfCorners = 0;
    allRegions[i].NumberOfElements = 0;
    allRegions[i].Elements.clear();
    allRegions[i].CornerPoints.clear();
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

  int numPoints = pd->GetNumberOfPoints();
  std::vector<int> cornerPoints;
  std::vector<int> isCornerPoint(numPoints);
  std::vector<int> isBoundaryPoint(numPoints);
  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, pointCellsValues);
    vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName.c_str(), i, pointCellsValues);
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

  fprintf(stdout, "K CHECKING REGIONS AND STUFFS %d\n", i);
    allRegions[i].NumberOfCorners = tempCornerPoints.size();
    fprintf(stdout,"NUM CORNS: %d\n", allRegions[i].NumberOfCorners);

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
          int pointCCWId = vtkSVGroupsSegmenter::GetCCWPoint(pd, tempNodes[j], cellId);
          int isBoundaryEdge = vtkSVGroupsSegmenter::CheckBoundaryEdge(pd, arrayName, cellId, tempNodes[j], pointCCWId);

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
    fprintf(stdout,"LETS SEE: %d %d\n", tempCornerPoints.size(), allRegions[i].CornerPoints.size());
  }

  return SV_OK;
}

// ----------------------
// GetCCWPoint
// ----------------------
int vtkSVGroupsSegmenter::GetCCWPoint(vtkPolyData *pd, const int pointId, const int cellId)
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

	if (position == 2)
	{
		position = 0;
		return pts[position];
	}
	else
	{
		position++;
		return pts[position];
	}
}

// ----------------------
// CheckBoundaryEdge
// ----------------------
int vtkSVGroupsSegmenter::CheckBoundaryEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1)
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
// CurveFitBoundaries
// ----------------------
int vtkSVGroupsSegmenter::CurveFitBoundaries(vtkPolyData *pd, std::string arrayName,
                                     std::vector<Region> allRegions)
{
  int numRegions = allRegions.size();

  std::vector<int> edgeValueCheck;
  for (int i=0; i<numRegions; i++)
  {
    for (int j=0; j<allRegions[i].BoundaryEdges.size(); j++)
    {
      fprintf(stdout,"Fitting curve edge %d of region %d\n", j, i);
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

      vtkSVGroupsSegmenter::SplineKnots(knots, numPoints, deg);

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
						minDist = dist;;
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
  return SV_OK;
}

void vtkSVGroupsSegmenter::SplineKnots(std::vector<int> &u, int n, int t)
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

void vtkSVGroupsSegmenter::SplineCurve(const std::vector<XYZ> &inp, int n, const std::vector<int> &knots, int t, std::vector<XYZ> &outp, int res)
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

void vtkSVGroupsSegmenter::SplinePoint(const std::vector<int> &u, int n, int t, double v, const std::vector<XYZ> &control, XYZ &output)
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

double vtkSVGroupsSegmenter::SplineBlend(int k, int t, const std::vector<int> &u, double v)
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

int vtkSVGroupsSegmenter::FixRegions(vtkPolyData *pd, std::string arrayName)
{

  int maxIters = 15;
  int iter=0;
  int allGood = 0;
  while(!allGood && iter < maxIters)
  {
     std::vector<Region> regions;
     if (vtkSVGroupsSegmenter::GetRegions(pd, arrayName, regions) != SV_OK)
     {
       fprintf(stderr,"Couldn't get regions\n");
       return SV_ERROR;
     }
    int numRegions = regions.size();

    int regionCount=0;
    for (int i=0; i<numRegions; i++)
    {
      if (regions[i].CornerPoints.size() <  4)
      {

        int myTEST = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(regions[i].Elements[0]);
        fprintf(stdout,"FOUND ERROR TRYING TO FIX %d, number of corners: %d\n", myTEST, regions[i].CornerPoints.size());
        allGood = 0;
        int maxNeighborRegion = -1;
        int needsFix = 0;
        std::vector<int> neighborRegions(numRegions);
        for (int j=0; j<numRegions; j++)
          neighborRegions[j] = 0.0;
        for (int j=0; j<regions[i].NumberOfElements; j++)
        {
          int cellId = regions[i].Elements[j];
          int cellValue = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(cellId);

          vtkNew(vtkIdList, neighborVals);
          vtkSVGeneralUtils::GetNeighborsCellsValues(pd, arrayName, cellId, neighborVals);
          if (neighborVals->GetNumberOfIds() > 1)
          {
            needsFix = 1;
            // boundary element
            for (int k=0; k<neighborVals->GetNumberOfIds(); k++)
            {
              if (neighborVals->GetId(k) != cellValue)
                neighborRegions[neighborVals->GetId(k)]++;
            }
          }
          else
          {
            if (neighborVals->GetNumberOfIds() == 1 &&
                neighborVals->GetId(0) != cellValue)
            {
              needsFix = 1;
              // Its surrounded! surrender
              maxNeighborRegion = neighborVals->GetId(0);
              break;
            }
          }
        }

        if (maxNeighborRegion == -1)
        {
          int maxNum = -1;
          for (int j=0; j<numRegions; j++)
          {
            if (neighborRegions[j] > maxNum)
            {
              maxNum = neighborRegions[j];
              maxNeighborRegion = j;
            }
          }
        }

        if (needsFix)
        {
          if (myTEST == 24)
            maxNeighborRegion = 26;
          for (int j=0; j<regions[i].NumberOfElements; j++)
          {
            int cellId = regions[i].Elements[j];
            pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(cellId, maxNeighborRegion);
            regions[maxNeighborRegion].NumberOfElements++;
            regions[maxNeighborRegion].Elements.push_back(cellId);
          }
          regions[i].NumberOfElements = 0;
          regions[i].CornerPoints.clear();
          regions[i].Elements.clear();
          regions[i].NumberOfCorners = 0;
          regions[i].BoundaryEdges.clear();
        }
      }
      else
        regionCount++;

    }
    if (regionCount == numRegions)
      allGood = 1;
    iter++;
  }

  return SV_OK;
}

int vtkSVGroupsSegmenter::Parameterize()
{
  std::vector<Region> patches;
  if (this->GetRegions(this->WorkPd, "PatchIds", patches) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  // Extract surface, triangulate, and subdivide polycube
  vtkNew(vtkPolyData, polycubePd);
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(this->Polycube);
  surfacer->Update();

  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(surfacer->GetOutput());
  triangulator->Update();

  vtkNew(vtkLinearSubdivisionFilter, subdivider);
  subdivider->SetInputData(triangulator->GetOutput());
  subdivider->SetNumberOfSubdivisions(4);
  subdivider->Update();
  polycubePd->DeepCopy(subdivider->GetOutput());
  fprintf(stdout,"JUST CHECK: %d\n", polycubePd->GetNumberOfPoints());

  int numPatches = patches.size();

  vtkNew(vtkAppendPolyData, appender);
  vtkNew(vtkAppendPolyData, mapAppender);

  for (int i=0; i<numPatches; i++)
  {
    int groupId = this->WorkPd->GetCellData()->GetArray(
     this->GroupIdsArrayName)->GetTuple1(patches[i].Elements[0]);
    int patchId = this->WorkPd->GetCellData()->GetArray(
     "PatchIds")->GetTuple1(patches[i].Elements[0]);
    int patchDir = patchId%6;

    // Get same group polycube
    // translate polygroup to regular spot ya know
    vtkNew(vtkUnstructuredGrid, rotPolycube);
    vtkNew(vtkMatrix4x4, rotMatrix0);
    vtkNew(vtkMatrix4x4, rotMatrix1);
    this->RotateGroupToGlobalAxis(this->Polycube, groupId, this->GroupIdsArrayName, rotPolycube, rotMatrix0, rotMatrix1);

    // Connect corner points of patches to polycube for boundary
    vtkNew(vtkPolyData, thresholdPd);
    thresholdPd->DeepCopy(this->WorkPd);
    vtkSVGeneralUtils::GiveIds(thresholdPd, "TmpInternalIds");
    vtkSVGeneralUtils::ThresholdPd(thresholdPd, patchId, patchId, 1, "PatchIds");

    // Set up boundary mapper
    vtkNew(vtkIntArray, boundaryCorners);
    boundaryCorners->SetNumberOfTuples(patches[i].CornerPoints.size());

    vtkNew(vtkIntArray, paraBoundaryCorners);
    paraBoundaryCorners->SetNumberOfTuples(patches[i].CornerPoints.size());
    fprintf(stdout,"PATCH: %d\n", patchId);

    fprintf(stdout,"Corner Points: ");
    for (int j=0; j<patches[i].CornerPoints.size(); j++)
    {
      int ptId = patches[i].CornerPoints[j];
      fprintf(stdout,"%d ", ptId);

      // Thresholded pt id
      int thresholdPtId = thresholdPd->GetPointData()->GetArray(
        "TmpInternalIds")->LookupValue(ptId);
      boundaryCorners->SetTuple1(j, thresholdPtId);

      // Paramteric space pt id
      vtkNew(vtkIdList, patchVals);
      vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, "PatchIds", ptId, patchVals);
      int paraPtId = -1;
      if (this->FindPointMatchingValues(rotPolycube, "PatchIds", patchVals, paraPtId) != SV_OK)
      {
        fprintf(stdout,"Could not find corresponding polycube point id\n");
        return SV_ERROR;
      }

      paraBoundaryCorners->SetTuple1(j, paraPtId);
    }
    fprintf(stdout,"\n");

    fprintf(stdout,"Poly Corner Points: ");
    for (int j=0; j<paraBoundaryCorners->GetNumberOfTuples(); j++)
      fprintf(stdout,"%.4f ", paraBoundaryCorners->GetTuple1(j));
    fprintf(stdout,"\n");

    vtkNew(vtkSVPointSetBoundaryMapper, boundaryMapper);
    boundaryMapper->SetPointSet(rotPolycube);
    boundaryMapper->SetPointSetBoundaryIds(paraBoundaryCorners);
    boundaryMapper->SetBoundaryIds(boundaryCorners);

    // Set up parameterizer
    vtkNew(vtkSVPlanarMapper, mapper);
    mapper->SetInputData(thresholdPd);
    mapper->SetBoundaryMapper(boundaryMapper);
    if (patchDir == 0 || patchDir == 2)
    {
      mapper->SetDir0(1);
      mapper->SetDir1(2);
      mapper->SetDir2(0);
    }
    else if (patchDir == 1 || patchDir == 3)
    {
      mapper->SetDir0(0);
      mapper->SetDir1(2);
      mapper->SetDir2(1);
    }
    else if (patchDir == 4 || patchDir == 5)
    {
      mapper->SetDir0(0);
      mapper->SetDir1(1);
      mapper->SetDir2(2);
    }
    mapper->Update();

    vtkNew(vtkPolyData, tmpPoly);
    tmpPoly->DeepCopy(mapper->GetOutput());

    rotMatrix0->Invert();
    rotMatrix1->Invert();

    // translate back to regular polycube spot
    vtkSVGeneralUtils::ApplyRotationMatrix(tmpPoly, rotMatrix1);
    vtkSVGeneralUtils::ApplyRotationMatrix(tmpPoly, rotMatrix0);

    std::string filename2 = "/Users/adamupdegrove/Desktop/tmp/Boundary_"+std::to_string(patchId)+".vtp";
    vtkSVIOUtils::WriteVTPFile(filename2, boundaryMapper->GetOutput());
    std::string filename4 = "/Users/adamupdegrove/Desktop/tmp/Mapping_"+std::to_string(patchId)+".vtp";
    vtkSVIOUtils::WriteVTPFile(filename4, mapper->GetOutput());

    appender->AddInputData(tmpPoly);


    // Then we have to think about volume stuff
    vtkNew(vtkPolyData, patchPolyPd);
    vtkSVGeneralUtils::ThresholdPd(polycubePd, patchId, patchId, 1, "PatchIds",
                                   patchPolyPd);

    vtkNew(vtkPolyData, patchMappedPd);
    this->InterpolateMapOntoTarget(patchPolyPd, thresholdPd, tmpPoly, patchMappedPd);

    mapAppender->AddInputData(patchMappedPd);

  }

  appender->Update();
  std::string filename = "/Users/adamupdegrove/Desktop/tmp/Mapping_All.vtp";
  vtkSVIOUtils::WriteVTPFile(filename, appender->GetOutput());

  mapAppender->Update();
  std::string filename5 = "/Users/adamupdegrove/Desktop/tmp/Mapped_Out.vtp";
  vtkSVIOUtils::WriteVTPFile(filename5, mapAppender->GetOutput());

  return SV_OK;
}

// ----------------------
// RotateGroupToGlobalAxis
// ----------------------
int vtkSVGroupsSegmenter::RotateGroupToGlobalAxis(vtkUnstructuredGrid *ug,
                                                  const int thresholdId,
                                                  std::string arrayName,
                                                  vtkUnstructuredGrid *rotUg,
                                                  vtkMatrix4x4 *rotMatrix0,
                                                  vtkMatrix4x4 *rotMatrix1)
{
  vtkNew(vtkUnstructuredGrid, thresholdUg);
  vtkSVGeneralUtils::ThresholdUg(ug, thresholdId, thresholdId, 1, arrayName, thresholdUg);

  double pts[3][3];
  for (int i=0; i<3; i++)
  {
    int ptId = thresholdUg->GetPointData()->GetArray("LocalPointIds")->LookupValue(i);
    thresholdUg->GetPoint(ptId, pts[i]);
  }

  double zVec[4], tmpVec[3];
  vtkMath::Subtract(pts[1], pts[0], zVec);
  vtkMath::Normalize(zVec);
  vtkMath::Subtract(pts[1], pts[2], tmpVec);
  vtkMath::Normalize(tmpVec);

  double yVec[3];
  vtkMath::Cross(zVec, tmpVec, yVec);
  vtkMath::Normalize(yVec);

  double realY[3], realZ[3];
  realY[0] = 0.0; realY[1] = 1.0; realY[2] = 0.0;
  realZ[0] = 0.0; realZ[1] = 0.0; realZ[2] = 1.0;

  vtkSVGeneralUtils::GetRotationMatrix(yVec, realY, rotMatrix0);
  double newZVec[4];
  rotMatrix0->MultiplyPoint(zVec, newZVec);

  vtkSVGeneralUtils::GetRotationMatrix(newZVec, realZ, rotMatrix1);

  vtkNew(vtkSVCleanUnstructuredGrid, ugCleaner);
  ugCleaner->SetInputData(ug);
  ugCleaner->Update();
  rotUg->DeepCopy(ugCleaner->GetOutput());

  vtkSVGeneralUtils::ApplyRotationMatrix(rotUg, rotMatrix0);
  vtkSVGeneralUtils::ApplyRotationMatrix(rotUg, rotMatrix1);

  return SV_OK;
}

// ----------------------
// FindPointMatchingValues
// ----------------------
int vtkSVGroupsSegmenter::FindPointMatchingValues(vtkPointSet *ps, std::string arrayName, vtkIdList *matchingVals, int &returnPtId)
{
  for (int i=0; i<ps->GetNumberOfPoints(); i++)
  {
    vtkNew(vtkIdList, pointCellValues);
    vtkSVGeneralUtils::GetPointCellsValues(ps, arrayName.c_str(), i, pointCellValues);
    pointCellValues->IntersectWith(matchingVals);

    if (pointCellValues->GetNumberOfIds() == matchingVals->GetNumberOfIds())
    {
      // We found it!
      returnPtId = i;
      return SV_OK;
    }
  }

  return SV_ERROR;
}

// ----------------------
// InterpolateMapOntoTarget
// ----------------------
int vtkSVGroupsSegmenter::InterpolateMapOntoTarget(vtkPolyData *sourceBasePd,
                                                         vtkPolyData *targetPd,
                                                         vtkPolyData *targetBasePd,
                                                         vtkPolyData *mappedPd)
{
  vtkNew(vtkSVMapInterpolator, interpolator);
  interpolator->SetInputData(0, sourceBasePd);
  interpolator->SetInputData(1, targetPd);
  interpolator->SetInputData(2, targetBasePd);
  interpolator->SetNumSourceSubdivisions(0);
  interpolator->Update();

  mappedPd->DeepCopy(interpolator->GetOutput());

  return SV_OK;
}
