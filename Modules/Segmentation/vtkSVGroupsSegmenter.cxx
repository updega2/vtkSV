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
#include "vtkAppendFilter.h"
#include "vtkExecutive.h"
#include "vtkCellArray.h"
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
#include "vtkSVCenterlinesEdgeWeightedCVT.h"
#include "vtkSVCleanUnstructuredGrid.h"
#include "vtkSVEdgeWeightedCVT.h"
#include "vtkSVEdgeWeightedSmoother.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVMathUtils.h"
#include "vtkSVPolyBallLine.h"
#include "vtkSVFindSeparateRegions.h"
#include "vtkSVPlanarMapper.h"
#include "vtkSVPointSetBoundaryMapper.h"
#include "vtkSVMapInterpolator.h"
#include "vtkSVLoftNURBSVolume.h"
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

  double point[3];
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
  this->WorkPd->BuildLinks();
  vtkDataArray *normalsArray =
    this->WorkPd->GetCellData()->GetArray("Normals");

  int stopCellNumber = ceil(this->WorkPd->GetNumberOfCells()*0.0001);
  vtkNew(vtkSVCenterlinesEdgeWeightedCVT, CVT);
  CVT->SetInputData(this->WorkPd);
  CVT->SetGenerators(this->CenterlinesWorkPd);
  CVT->SetNumberOfRings(2);
  CVT->SetThreshold(stopCellNumber);
  CVT->SetUseCurvatureWeight(1);
  CVT->SetPatchIdsArrayName(this->GroupIdsArrayName);
  CVT->SetCVTDataArrayName("Normals");
  CVT->SetGroupIdsArrayName(this->GroupIdsArrayName);
  CVT->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
  CVT->SetBlankingArrayName(this->BlankingArrayName);
  CVT->SetUseRadiusInformation(this->UseRadiusInformation);
  CVT->Update();

  vtkNew(vtkSVEdgeWeightedSmoother, smoother);
  smoother->SetInputData(CVT->GetOutput());
  smoother->SetGenerators(this->CenterlinesWorkPd);
  smoother->SetNumberOfRings(2);
  smoother->SetThreshold(stopCellNumber);
  smoother->SetUseCurvatureWeight(0);
  smoother->SetNoInitialization(1);
  smoother->SetPatchIdsArrayName(this->GroupIdsArrayName);
  smoother->SetCVTDataArrayName("Normals");
  smoother->Update();

  this->WorkPd->DeepCopy(smoother->GetOutput());

  if (this->FixGroups() != SV_OK)
  {
    vtkErrorMacro("Error in correcting groups");
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

  // Get new normals
  normaler->SetInputData(this->WorkPd);
  normaler->ComputePointNormalsOff();
  normaler->ComputeCellNormalsOn();
  normaler->SplittingOff();
  normaler->Update();
  this->WorkPd->DeepCopy(normaler->GetOutput());
  this->WorkPd->BuildLinks();

  // Add array for new cell normals on surface
  vtkNew(vtkDoubleArray, newCellNormals);
  newCellNormals->SetName("CenterlinesBasedCellNormals");
  newCellNormals->SetNumberOfComponents(3);
  newCellNormals->SetNumberOfTuples(numberOfCells);

  // Get all group ids
  vtkNew(vtkIdList, groupIds);
  for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
  {
    int groupVal = this->WorkPd->GetCellData()->GetArray(
        this->GroupIdsArrayName)->GetTuple1(i);
    groupIds->InsertUniqueId(groupVal);
  }
  vtkSortDataArray::Sort(groupIds);
  int numGroups = groupIds->GetNumberOfIds();

  vtkSVGeneralUtils::GiveIds(this->WorkPd, "TmpInternalIds");
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);
    vtkNew(vtkPolyData, branchPd);
    vtkSVGeneralUtils::ThresholdPd(this->WorkPd, groupId, groupId, 1,
      this->GroupIdsArrayName, branchPd);
    branchPd->BuildLinks();

    vtkNew(vtkPolyData, centerlineBranchPd);
    vtkSVGeneralUtils::ThresholdPd(this->CenterlinesWorkPd, groupId, groupId, 1,
      this->GroupIdsArrayName, centerlineBranchPd);
    centerlineBranchPd->BuildLinks();

    // for each group, compute the clipping array, clip, add group ids array and append.
    vtkNew(vtkSVPolyBallLine, groupTubes);
    groupTubes->SetInput(centerlineBranchPd);
    groupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
    groupTubes->SetUseRadiusInformation(this->UseRadiusInformation);
    groupTubes->UsePointNormalOff();
    groupTubes->UseRadiusWeightingOff();
    groupTubes->UseLocalCoordinatesOn();
    groupTubes->SetLocalCoordinatesArrayName("Local");
    groupTubes->BuildLocator();

    int branchNumberOfCells = branchPd->GetNumberOfCells();
    // Loop through points to evaluate function at each point
    fprintf(stdout,"Computing closest centerline points per cell...\n");
    for (int k=0; k<branchNumberOfCells; k++)
    {
      // Get cell point coords
      double pts[3][3];
      vtkIdType npts, *ptids;
      branchPd->GetCellPoints(k, npts, ptids);
      for (int j=0; j<npts; j++)
        branchPd->GetPoint(ptids[j], pts[j]);

      // Get center
      double center[3];
      vtkTriangle::TriangleCenter(pts[0], pts[1], pts[2], center);

      // Evaluate function at point!
      groupTubes->EvaluateFunction(center);

      // Now get last local coords and use rotation matrix to set new normals
      double localX[3], localY[3], localZ[3];
      groupTubes->GetLastLocalCoordX(localX);
      groupTubes->GetLastLocalCoordY(localY);
      groupTubes->GetLastLocalCoordZ(localZ);

      // Compute the rotation from global coordinate system to centerlines
      // local coordinate system
      double rotMat[9];
      this->ComputeRotationMatrix(localX, localY, localZ, rotMat);

      //Get real cell id
      int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(k);

      double cellNormal[3];
      this->WorkPd->GetCellData()->GetArray("Normals")->GetTuple(realCellId, cellNormal);

      // Apply rotation matrix to the normal to get the new normal
      double newNormal[3];
      for (int j=0; j<3; j++)
      {
        newNormal[j] = rotMat[j*3]*cellNormal[0] +
                       rotMat[(j*3)+1]*cellNormal[1] +
                       rotMat[(j*3)+2]*cellNormal[2];
      }
      //fprintf(stdout,"SETTING REAL CELLID: %d to %.6f %.6f %.6f\n", realCellId, newNormal[0], newNormal[1], newNormal[2]);
      newCellNormals->SetTuple(realCellId, newNormal);

    }

  }
  this->WorkPd->GetCellData()->RemoveArray("TmpInternalIds");
  this->WorkPd->GetPointData()->RemoveArray("TmpInternalIds");
  this->WorkPd->GetCellData()->AddArray(newCellNormals);

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

  if (this->RunEdgeWeightedCVT(this->WorkPd, generatorsPd) != SV_OK)
  {
    vtkErrorMacro("Error in cvt");
    return SV_ERROR;
  }

  if (this->FixEndPatches() != SV_OK)
  {
    vtkErrorMacro("Error fixing end patches");
    return SV_ERROR;
  }

  for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
  {
    int patchVal = this->WorkPd->GetCellData()->GetArray("PatchIds")->GetTuple1(i);
    if (patchVal == 2)
      patchVal = 0;
    if (patchVal == 3)
      patchVal = 1;
    this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(i, patchVal);
  }

  vtkNew(vtkSVEdgeWeightedSmoother, smoother2);
  smoother2->SetInputData(this->WorkPd);
  smoother2->SetGenerators(this->WorkPd);
  smoother2->SetNumberOfRings(2);
  smoother2->SetThreshold(stopCellNumber);
  smoother2->SetUseCurvatureWeight(1);
  smoother2->SetNoInitialization(1);
  smoother2->SetPatchIdsArrayName("PatchIds");
  smoother2->SetCVTDataArrayName("Normals");
  smoother2->Update();

  this->WorkPd->DeepCopy(smoother2->GetOutput());

  vtkNew(vtkIdList, targetPatches);
  targetPatches->SetNumberOfIds(4);
  for (int i=0; i<4; i++)
    targetPatches->SetId(i, i);
  if (this->CorrectSpecificCellBoundaries(this->WorkPd, "PatchIds", targetPatches) != SV_OK)
  {
    vtkErrorMacro("Could not correcto boundaries of surface");
    return SV_ERROR;
  }

  if (this->SmoothSpecificBoundaries(this->WorkPd, "PatchIds", targetPatches) != SV_OK)
  {
    vtkErrorMacro("Could not smootho boundaries of surface");
    return SV_ERROR;
  }
  std::vector<Region> patchRegions;
  if (this->GetSpecificRegions(this->WorkPd, "PatchIds", patchRegions, targetPatches) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }
  if (this->CurveFitBoundaries(this->WorkPd, "PatchIds", patchRegions) != SV_OK)
  {
    vtkErrorMacro("Could not curve fit boundaries of surface");
    return SV_ERROR;
  }

  //vtkNew(vtkIdList, addVals);
  //addVals->SetNumberOfIds(numGroups);
  //for (int i=0; i<numGroups; i++)
  //  addVals->SetId(i, 6*i);

  //vtkNew(vtkIdList, patchVals);
  //for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
  //{
  //  int patchVal = this->WorkPd->GetCellData()->GetArray("PatchIds")->GetTuple1(i);
  //  int groupVal = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(i);
  //  int newVal = patchVal + (addVals->GetId(groupIds->IsId(groupVal)));
  //  this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(i, newVal);
  //  patchVals->InsertUniqueId(newVal);
  //}

  //if (this->FixPatchesByGroup() != SV_OK)
  //{
  //  fprintf(stderr,"Couldn't fix patches\n");
  //  //return SV_ERROR;
  //}
  //if (this->FixPatchesWithPolycube() != SV_OK)
  //{
  //  fprintf(stderr,"Couldn't fix patches\n");
  //  //return SV_ERROR;
  //}

  //std::vector<Region> finalRegions;
  //targetPatches->Reset();
  //targetPatches->SetNumberOfIds(numGroups*4);
  //for (int i=0; i<numGroups; i++)
  //{
  //  for (int j=0; j<4; j++)
  //    targetPatches->SetId(4*i+j, 6*i+j);
  //}
  //if (this->GetSpecificRegions(this->WorkPd, "PatchIds", finalRegions, targetPatches) != SV_OK)
  //{
  //  vtkErrorMacro("Couldn't get patches");
  //  return SV_ERROR;
  //}
  //if (this->CurveFitBoundaries(this->WorkPd, "PatchIds", finalRegions) != SV_OK)
  //{
  //  vtkErrorMacro("Could not curve fit boundaries of surface");
  //  return SV_ERROR;
  //}

  //if (this->FixPatchesWithPolycube() != SV_OK)
  //{
  //  fprintf(stderr,"POOP WEINER CAKE SAUCE\n");
  //  //return SV_ERROR;
  //}

  //// NOW PARAMETERIZE!!, WIILL BE MOVED to vtkSVPolycubeParameterizer
  //// TODO: RENAME THIS CLASS TO vtkSVCenterlinesSegmenter

  //vtkNew(vtkPolyData, fullMapPd);
  //if (this->Parameterize(fullMapPd) != SV_OK)
  //{
  //  fprintf(stderr,"WRONG\n");
  //  return SV_ERROR;
  //}

  //if (this->Centerlines->GetNumberOfCells() == 1)
  //{
  //  vtkNew(vtkStructuredGrid, paraHexMesh);
  //  if (this->FormParametricHexMesh(paraHexMesh) != SV_OK)
  //  {
  //    fprintf(stderr,"Couldn't do the dirt\n");
  //    return SV_ERROR;
  //  }

  //  vtkNew(vtkStructuredGrid, realHexMesh);
  //  if (this->MapVolume(paraHexMesh, fullMapPd, realHexMesh) != SV_OK)
  //  {
  //    fprintf(stderr,"Couldn't do the dirt\n");
  //    return SV_ERROR;
  //  }
  //  // Set up the volume
  //  vtkNew(vtkUnstructuredGrid, emptyGrid);
  //  vtkNew(vtkSVLoftNURBSVolume, lofter);
  //  lofter->SetInputData(emptyGrid);
  //  lofter->SetInputGrid(realHexMesh);
  //  lofter->SetUDegree(2);
  //  lofter->SetVDegree(2);
  //  lofter->SetWDegree(2);
  //  lofter->SetUnstructuredGridUSpacing(0.1);
  //  lofter->SetUnstructuredGridVSpacing(0.01);
  //  lofter->SetUnstructuredGridWSpacing(0.1);
  //  lofter->SetUKnotSpanType("average");
  //  lofter->SetUParametricSpanType("chord");
  //  lofter->SetVKnotSpanType("average");
  //  lofter->SetVParametricSpanType("chord");
  //  lofter->SetWKnotSpanType("average");
  //  lofter->SetWParametricSpanType("chord");
  //  lofter->Update();

  //  std::string filename = "/Users/adamupdegrove/Desktop/tmp/TEST_FINAL.vtu";
  //  vtkSVIOUtils::WriteVTUFile(filename, lofter->GetOutput());

  //}


  return SV_OK;
}

// ----------------------
// RunEdgeWeightedCVT
// ----------------------
int vtkSVGroupsSegmenter::RunEdgeWeightedCVT(vtkPolyData *pd, vtkPolyData *generatorsPd)
{
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

  int allGood = 0;
  int iter = 0;
  int maxIters = 100;

  while(!allGood && iter < maxIters)
  {
    int checkFix = 0;
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
      if (neiSize > 1 && (neiTmpIds->GetId(0) == neiTmpIds->GetId(1) ||
                          neiTmpIds->GetId(1) == neiTmpIds->GetId(2)))
      {
        checkFix = 1;
        int maxVal, maxCount;
        vtkSVGroupsSegmenter::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

        cellIds->SetTuple1(i, maxVal);
        tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
      }
    }
    if (checkFix == 0)
      allGood = 1;
    fprintf(stdout, "CELL BOUNDARY FIX ITER: %d\n", iter);
    iter++;
  }

  return 1;
}

// ----------------------
// CorrectSpecificCellBoundaries
// ----------------------
int vtkSVGroupsSegmenter::CorrectSpecificCellBoundaries(vtkPolyData *pd, std::string cellArrayName, vtkIdList *targetRegions)
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
    int checkFix = 0;
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
      if (neiSize > 1 && (neiTmpIds->GetId(0) == neiTmpIds->GetId(1) ||
                          neiTmpIds->GetId(1) == neiTmpIds->GetId(2)))
      {
        checkFix = 1;
        int maxVal, maxCount;
        vtkSVGroupsSegmenter::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

        cellIds->SetTuple1(i, maxVal);
        tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
      }
    }
    if (checkFix == 0)
      allGood = 1;
    fprintf(stdout, "CELL BOUNDARY FIX ITER: %d\n", iter);
    iter++;
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


  for (int i=0; i<numPoints; i++)
  {
    if (isBoundaryPoint[i])
    {
      vtkNew(vtkIdList, pointCellsValues);
      vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName.c_str(),
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
// SmoothSpecificBoundaries
// ----------------------
int vtkSVGroupsSegmenter::SmoothSpecificBoundaries(vtkPolyData *pd, std::string arrayName, vtkIdList *targetRegions)
{
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
      vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName.c_str(),
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

    allRegions[i].NumberOfCorners = tempCornerPoints.size();
    //fprintf(stdout,"NUM CORNS: %d\n", allRegions[i].NumberOfCorners);

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
    //fprintf(stdout,"LETS SEE: %d %d\n", tempCornerPoints.size(), allRegions[i].CornerPoints.size());
  }
  //fprintf(stdout,"DONE GETTING REGIONS\n");

  return SV_OK;
}

// ----------------------
// GetSpecificRegions
// ----------------------
int vtkSVGroupsSegmenter::GetSpecificRegions(vtkPolyData *pd, std::string arrayName,
                                             std::vector<Region> &allRegions,
                                             vtkIdList *targetRegions)
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

  int numPoints = pd->GetNumberOfPoints();
  std::vector<int> cornerPoints;
  std::vector<int> isCornerPoint(numPoints);
  std::vector<int> isBoundaryPoint(numPoints);
  std::vector<int> isNonTargetBoundaryPoint(numPoints);
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
    {
      isNonTargetBoundaryPoint[i] = 1;
      if (targetRegions->IsId(pointCellsValues->GetId(0)) != -1 &&
          targetRegions->IsId(pointCellsValues->GetId(1)) != -1)
        isBoundaryPoint[i] = 1;
      else
        isBoundaryPoint[i] = 0;
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
    //fprintf(stdout,"NUM CORNS: %d\n", allRegions[i].NumberOfCorners);

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
          int pointCCWId = vtkSVGroupsSegmenter::GetCCWPoint(pd, tempNodes[j], cellId);
          int isBoundaryEdge = vtkSVGroupsSegmenter::CheckBoundaryEdge(pd, arrayName, cellId, tempNodes[j], pointCCWId);

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
    //fprintf(stdout,"LETS SEE: %d %d\n", tempCornerPoints.size(), allRegions[i].CornerPoints.size());
  }
  //fprintf(stdout,"DONE GETTING REGIONS\n");
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
      //fprintf(stdout,"Fitting curve edge %d of region %d\n", j, i);
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

int vtkSVGroupsSegmenter::FixEndPatches()
{
  vtkNew(vtkIdList, targetRegions);
  targetRegions->SetNumberOfIds(2);
  targetRegions->SetId(0, 4);
  targetRegions->SetId(1, 5);

  std::vector<Region> endRegions;
  if (this->GetSpecificRegions(this->WorkPd, "PatchIds", endRegions, targetRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  std::vector<int> individualFix;
  std::vector<int> wholePatchFix;
  this->CheckEndPatches(endRegions, individualFix, wholePatchFix);

  vtkDataArray *cellNormals = this->WorkPd->GetCellData()->GetArray("Normals");

  for (int i=0; i<individualFix.size(); i++)
  {
    fprintf(stdout,"FIXING INDIVIDUAL END PATCH\n");
    int badPatch = individualFix[i];
    double avgNorm[3]; avgNorm[0] = 0.0; avgNorm[1] = 0.0; avgNorm[2] = 0.0;
    for (int j=0; j<endRegions[badPatch].NumberOfElements; j++)
    {
      double cellNorm[3];
      cellNormals->GetTuple(endRegions[badPatch].Elements[j], cellNorm);
      for (int k=0; k<3; k++)
        avgNorm[k] += cellNorm[k];
    }
    vtkMath::MultiplyScalar(avgNorm, 1./endRegions[badPatch].NumberOfElements);

    for (int j=0; j<endRegions[badPatch].NumberOfElements; j++)
    {
      int cellId = endRegions[badPatch].Elements[j];
      double cellNorm[3];
      cellNormals->GetTuple(cellId, cellNorm);

      double compare = vtkMath::Dot (cellNorm, avgNorm);
      int removeCellValue = this->WorkPd->GetCellData()->GetArray(
        "PatchIds")->GetTuple1(cellId);
      if (compare <= 0.95)
      {
        vtkNew(vtkIdList, neighborValues);
        vtkSVGeneralUtils::GetNeighborsCellsValues(this->WorkPd, "PatchIds",
                                                   cellId,
                                                   neighborValues);
        int newCellValue = -1;;
        for (int k=0; k<neighborValues->GetNumberOfIds(); k++)
        {
          if (neighborValues->GetId(k) != removeCellValue)
          {
            newCellValue = neighborValues->GetId(k);
            break;
          }
        }
        if (newCellValue != -1)
          this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellId, newCellValue);
      }
    }
  }

  if (wholePatchFix.size() > 0)
  {
    fprintf(stdout,"FIXING WHOLE END PATCHES\n");
    vtkNew(vtkPolyData, workPdCopy);
    workPdCopy->DeepCopy(this->WorkPd);

    // Set up generators
    vtkNew(vtkPoints, generatorsPts);
    generatorsPts->SetNumberOfPoints(4);
    generatorsPts->SetPoint(0, 1.0, 0.0, 0.0);
    generatorsPts->SetPoint(1, 0.0, 1.0, 0.0);
    generatorsPts->SetPoint(2, -1.0, 0.0, 0.0);
    generatorsPts->SetPoint(3, 0.0, -1.0, 0.0);

    vtkNew(vtkPolyData, generatorsPd);
    generatorsPd->SetPoints(generatorsPts);

    if (this->RunEdgeWeightedCVT(workPdCopy, generatorsPd) != SV_OK)
    {
      vtkErrorMacro("Error in cvt");
      return SV_ERROR;
    }

    for (int i=0; i<wholePatchFix.size(); i++)
    {
      int badPatch = wholePatchFix[i];
      for (int j=0; j<endRegions[badPatch].NumberOfElements; j++)
      {
        int cellId = endRegions[badPatch].Elements[j];
        int newCellValue = workPdCopy->GetCellData()->GetArray("PatchIds")->GetTuple1(cellId);
        this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellId, newCellValue);
      }
    }
  }

  return SV_OK;
}

int vtkSVGroupsSegmenter::CheckEndPatches(std::vector<Region> endRegions,
                                          std::vector<int> &individualFix,
                                          std::vector<int> &wholePatchFix)
{
  // Need to fix both the random weird patch that goes over and the
  // entire groups that shouldn't be there
  vtkDataArray *cellNormals = this->WorkPd->GetCellData()->GetArray("Normals");

  int numRegions = endRegions.size();

  for (int i=0; i<numRegions; i++)
  {
    double avgNorm[3]; avgNorm[0] = 0.0; avgNorm[1] = 0.0; avgNorm[2] = 0.0;
    for (int j=0; j<endRegions[i].NumberOfElements; j++)
    {
      double cellNorm[3];
      cellNormals->GetTuple(endRegions[i].Elements[j], cellNorm);
      for (int k=0; k<3; k++)
        avgNorm[k] += cellNorm[k];
    }
    vtkMath::MultiplyScalar(avgNorm, 1./endRegions[i].NumberOfElements);

    int numClose=0;
    for (int j=0; j<endRegions[i].NumberOfElements; j++)
    {
      double cellNorm[3];
      cellNormals->GetTuple(endRegions[i].Elements[j], cellNorm);

      double compare = vtkMath::Dot (cellNorm, avgNorm);
      if (compare > 0.95)
        numClose+=1;
    }

    fprintf(stdout,"EXPECTED: %d, CLOSE: %d\n", endRegions[i].NumberOfElements, numClose);
    if (numClose != endRegions[i].NumberOfElements)
    {
      if (numClose < (3./4)*endRegions[i].NumberOfElements)
        wholePatchFix.push_back(i);
      else
        individualFix.push_back(i);
    }
  }

  return SV_OK;
}

int vtkSVGroupsSegmenter::FixGroups()
{
  // Clean up groups
  int allGood = 0;
  int iter = 0;
  int maxIters = 15;
  while(!allGood)
  {
    int checkFix=0;
    for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
    {
      int groupVal = this->WorkPd->GetCellData()->GetArray(
        this->GroupIdsArrayName)->GetTuple1(i);
      if (groupVal == -1)
      {
        checkFix = 1;
        vtkNew(vtkIdList, neighborValues);
        vtkSVGeneralUtils::GetNeighborsCellsValues(this->WorkPd,
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
          this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->SetTuple1(i, newCellValue);
      }
    }
    if (checkFix == 0)
      allGood = 1;
    fprintf(stdout,"GROUP FIX ITER: %d\n", iter);
    iter++;
  }

  return SV_OK;
}

int vtkSVGroupsSegmenter::FixPatchesByGroup()
{
  vtkNew(vtkIdList, groupIds);
  for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
  {
    int groupVal = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(i);
    groupIds->InsertUniqueId(groupVal);
  }
  vtkSortDataArray::Sort(groupIds);
  int numGroups = groupIds->GetNumberOfIds();

  // First check group by group
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);
    vtkNew(vtkPolyData, branchPd);

    vtkSVGeneralUtils::ThresholdPd(this->WorkPd, groupId, groupId, 1,
      this->GroupIdsArrayName, branchPd);
    branchPd->BuildLinks();

    vtkNew(vtkIdList, patchIds);
    for (int j=0; j<branchPd->GetNumberOfCells(); j++)
    {
      int patchId = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(j);
      patchIds->InsertUniqueId(patchId);
    }

    int centCellId = this->Centerlines->GetCellData()->GetArray(
     this->GroupIdsArrayName)->LookupValue(groupId);

    vtkIdType *pts, npts;
    this->Centerlines->GetCellPoints(centCellId, npts, pts);

    vtkNew(vtkIdList, pointCells0);
    this->Centerlines->GetPointCells(pts[0], pointCells0);
    vtkNew(vtkIdList, pointCellsN);
    this->Centerlines->GetPointCells(pts[npts-1], pointCellsN);

     std::vector<Region> patches;
     if (vtkSVGroupsSegmenter::GetSpecificRegions(this->WorkPd, "PatchIds", patches, patchIds) != SV_OK)
     {
       fprintf(stderr,"Couldn't get patches\n");
       return SV_ERROR;
     }
    int numPatches = patches.size();

    int numCells0 = pointCells0->GetNumberOfIds();
    int numCellsN = pointCellsN->GetNumberOfIds();

    int fixedRegions = 0;
    if (numCells0 == 1 && numCellsN == 1)
    {
      if (numPatches != 6)
      {
        fprintf(stdout,"Group %d should be six patches and there are %d! Try fix\n", groupId, numPatches);
        this->FixSpecificRegions(this->WorkPd, "PatchIds", patchIds, 0);
        fixedRegions = 1;
      }
    }
    else if (numCells0 == 1 || numCellsN == 1)
    {
      if (numPatches != 5)
      {
        fprintf(stdout,"Group %d should have five patches and there are %d! Try fix\n", groupId, numPatches);
        this->FixSpecificRegions(this->WorkPd, "PatchIds", patchIds, 0);
        fixedRegions = 1;
      }
    }
    else
    {
      if (numPatches != 4)
      {
        fprintf(stdout,"Group %d should have four patches and there are %d! Try fix\n", groupId, numPatches);
        this->FixSpecificRegions(this->WorkPd, "PatchIds", patchIds, 1);
        fixedRegions = 1;
      }
    }
  }

  return SV_OK;
}

int vtkSVGroupsSegmenter::FixPatchesWithPolycube()
{
  // Then check everything
  // Extract surface, triangulate, and subdivide polycube
  vtkNew(vtkPolyData, polycubePd);
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(this->Polycube);
  surfacer->Update();

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(surfacer->GetOutput());
  cleaner->Update();

  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(cleaner->GetOutput());
  triangulator->Update();

  polycubePd->DeepCopy(triangulator->GetOutput());
  polycubePd->BuildLinks();

  fprintf(stdout,"GETTING SURFACE REGIONS\n");
  std::vector<Region> surfacePatches;
  if (this->GetRegions(this->WorkPd, "PatchIds", surfacePatches) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  fprintf(stdout,"GETTING POLYCUBE REGIONS\n");
  std::vector<Region> polycubePatches;
  if (this->GetRegions(polycubePd, "PatchIds", polycubePatches) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  int numSurfacePatches =  surfacePatches.size();
  int numPolycubePatches = polycubePatches.size();

  if (numSurfacePatches != numPolycubePatches)
  {
    vtkErrorMacro("The number of patches on the polycube and the surface must match!");
    int maxCluster = -1;
    for (int i=0; i<numSurfacePatches; i++)
    {
      if (surfacePatches[i].IndexCluster > maxCluster)
        maxCluster = surfacePatches[i].IndexCluster;
    }
    maxCluster++;
    std::vector<int> testPatches(maxCluster, 0);
    for (int i=0; i<numSurfacePatches; i++)
      testPatches[surfacePatches[i].IndexCluster]++;

    for (int i=0; i<maxCluster; i++)
    {
      if (testPatches[i] > 1)
        fprintf(stdout,"Surface patch %d has more than one cluster\n", i);
    }

    return SV_ERROR;
  }

  int allGood = 0;
  int iter = 0;
  int maxIters = 15;
  while(!allGood && iter < maxIters)
  {
    int checkFix = 0;
    for (int i=0; i<numPolycubePatches; i++)
    {
      int groupId = polycubePd->GetCellData()->GetArray(
       this->GroupIdsArrayName)->GetTuple1(polycubePatches[i].Elements[0]);
      int patchId = polycubePd->GetCellData()->GetArray(
       "PatchIds")->GetTuple1(polycubePatches[i].Elements[0]);
      int patchDir = patchId%6;

      for (int j=0; j<numSurfacePatches; j++)
      {
        int surfacePatchId = this->WorkPd->GetCellData()->GetArray(
         "PatchIds")->GetTuple1(surfacePatches[j].Elements[0]);

        if (patchId == surfacePatchId)
        {
          if (polycubePatches[i].CornerPoints.size() !=
              surfacePatches[j].CornerPoints.size())
          {
            fprintf(stderr,"Expected %d corners in patch %d, but found %d\n",
                    polycubePatches[i].CornerPoints.size(),
                    patchId,
                    surfacePatches[j].CornerPoints.size());
            std::vector<int> lowCorners;
            std::vector<int> highCorners;
            vtkNew(vtkIdList, checkValues);
            for (int k=0; k<surfacePatches[j].CornerPoints.size(); k++)
            {
              int cornerPtId = surfacePatches[j].CornerPoints[k];
              vtkNew(vtkIdList, pointCellsValues);
              vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, "PatchIds", cornerPtId, pointCellsValues);
              int matchingPtId = -1;
              if (this->FindPointMatchingValues(polycubePd, "PatchIds", pointCellsValues, matchingPtId) != SV_OK)
              {
                fprintf(stderr,"  We found a problem corner: %d\n", cornerPtId);
                if (matchingPtId != -1)
                {
                  vtkNew(vtkIdList, cubeValues);
                  vtkSVGeneralUtils::GetPointCellsValues(polycubePd, "PatchIds", matchingPtId, cubeValues);
                  fprintf(stdout,"  Point is attached to %lld regions, but cube has %lld regions\n", pointCellsValues->GetNumberOfIds(), cubeValues->GetNumberOfIds());
                  fprintf(stdout,"  Expected regions: ");
                  for (int l=0; l<cubeValues->GetNumberOfIds(); l++)
                    fprintf(stdout,"%d ", cubeValues->GetId(l));
                  fprintf(stdout,"\n");
                  fprintf(stdout,"  But has regions: ");
                  for (int l=0; l<pointCellsValues->GetNumberOfIds(); l++)
                    fprintf(stdout,"%d ", pointCellsValues->GetId(l));
                  fprintf(stdout,"\n");
                  if (pointCellsValues->GetNumberOfIds() < cubeValues->GetNumberOfIds())
                  {
                    if (checkValues->GetNumberOfIds() == 0)
                    {
                      lowCorners.push_back(surfacePatches[j].CornerPoints[k]);
                      for (int l=0; l<cubeValues->GetNumberOfIds(); l++)
                        checkValues->InsertNextId(cubeValues->GetId(l));
                    }
                    else
                    {
                      int prevNum = cubeValues->GetNumberOfIds();
                      cubeValues->IntersectWith(checkValues);
                      if (cubeValues->GetNumberOfIds() == checkValues->GetNumberOfIds() &&
                          prevNum == cubeValues->GetNumberOfIds())
                      {
                        lowCorners.push_back(surfacePatches[j].CornerPoints[k]);
                      }
                    }
                  }
                  else
                    highCorners.push_back(surfacePatches[j].CornerPoints[k]);
                }
              }
            }
            if (lowCorners.size() == 2)
            {
              checkFix = 1;
              // Nice, we found a potentially fixable spot
              fprintf(stdout,"Attempting to merge corners\n");
              int badCorner0 = lowCorners[0];
              int badCorner1 = lowCorners[1];

              int fixEdge = -1;
              for (int k=0; k<surfacePatches[j].BoundaryEdges.size(); k++)
              {
                int numEdges = surfacePatches[j].BoundaryEdges[k].size();
                if ((surfacePatches[j].BoundaryEdges[k][0] == badCorner0 &&
                     surfacePatches[j].BoundaryEdges[k][numEdges-1] == badCorner1) ||
                    (surfacePatches[j].BoundaryEdges[k][numEdges-1] == badCorner0 &&
                     surfacePatches[j].BoundaryEdges[k][0] == badCorner1))
                {
                  if (numEdges < 5)
                  {
                    fixEdge = k;
                  }
                }
              }
              if (fixEdge != -1)
              {
                int numEdges = surfacePatches[j].BoundaryEdges[fixEdge].size();
                vtkNew(vtkIdList, pCV0);
                vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, "PatchIds", surfacePatches[j].BoundaryEdges[fixEdge][0], pCV0);
                vtkNew(vtkIdList, pCV1);
                vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, "PatchIds", surfacePatches[j].BoundaryEdges[fixEdge][numEdges-1], pCV1);
                pCV1->IntersectWith(pCV0);
                int newCellValue = -1;
                int removeCellValue = -1;
                for (int k=0; k<pCV0->GetNumberOfIds(); k++)
                {
                  if (pCV1->IsId(pCV0->GetId(k)) == -1)
                    newCellValue = pCV0->GetId(k);
                  if (pCV1->IsId(pCV0->GetId(k)) != -1)
                    removeCellValue = pCV0->GetId(k);
                }
                //fprintf(stdout,"NEW CELL VALUE: %d\n", newCellValue);
                //fprintf(stdout,"GETTING RID OF: %d\n", removeCellValue);
                for (int k=0; k<surfacePatches[j].BoundaryEdges[fixEdge].size()-1; k++)
                {
                  int ptId0 = surfacePatches[j].BoundaryEdges[fixEdge][k];
                  vtkNew(vtkIdList, cellNeighbors);
                  this->WorkPd->GetPointCells(ptId0, cellNeighbors);
                  for (int l=0; l<cellNeighbors->GetNumberOfIds(); l++)
                  {
                    int val0 = this->WorkPd->GetCellData()->GetArray("PatchIds")->GetTuple1(cellNeighbors->GetId(l));
                    if (val0 == removeCellValue)
                      this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellNeighbors->GetId(l), newCellValue);
                  }
                }
              }

            }
            if (highCorners.size() > 0)
            {
              fprintf(stdout,"Attempting to separate a regions from corner\n");
              checkFix = 1;
              for (int k=0; k<highCorners.size(); k++)
              {
                int badCorner0 = highCorners[k];
                vtkNew(vtkIdList, badCornerValues);
                vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, "PatchIds", badCorner0, badCornerValues);
                int matchingPtId = -1;
                if (this->FindPointMatchingValues(polycubePd, "PatchIds", badCornerValues, matchingPtId) != SV_OK)
                {
                  if (matchingPtId != -1)
                  {
                    vtkNew(vtkIdList, cubeValues);
                    vtkSVGeneralUtils::GetPointCellsValues(polycubePd, "PatchIds", matchingPtId, cubeValues);
                    for (int l=0; l<badCornerValues->GetNumberOfIds(); l++)
                    {
                      int fixed=0;
                      if (cubeValues->IsId(badCornerValues->GetId(l)) == -1)
                      {
                        int removeCellValue = badCornerValues->GetId(l);
                        vtkNew(vtkIdList, pointCells);
                        this->WorkPd->GetPointCells(badCorner0, pointCells);
                        for (int r=0; r<pointCells->GetNumberOfIds(); r++)
                        {
                          if (this->WorkPd->GetCellData()->GetArray("PatchIds")->GetTuple1(pointCells->GetId(r)) == removeCellValue)
                          {
                            int newCellValue = -1;
                            vtkNew(vtkIdList, cellNeighborValues);
                            vtkSVGeneralUtils::GetNeighborsCellsValues(this->WorkPd,
                                                                     "PatchIds",
                                                                       pointCells->GetId(r),
                                                                       cellNeighborValues);
                            for (int s=0; s<cellNeighborValues->GetNumberOfIds(); s++)
                            {
                              int potNewVal = cellNeighborValues->GetId(s);
                              if (cubeValues->IsId(potNewVal) != -1)
                              {
                                newCellValue = potNewVal;
                                //fprintf(stdout,"NEWVALUE: %d\n", newCellValue);
                                break;
                              }
                            }
                            if (newCellValue != -1)
                            {
                              this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(pointCells->GetId(r), newCellValue);
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
          break;
        }
      }
    }
    if (checkFix == 0)
      allGood = 1;
    fprintf(stdout, "POLYCUBE FIX ITER: %d\n", iter);
    iter++;
  }

  return SV_OK;
}

int vtkSVGroupsSegmenter::FixSpecificRegions(vtkPolyData *pd, std::string arrayName,
                                             vtkIdList *targetRegions, const int noEndRegions)
{

  int maxIters = 15;
  int iter=0;
  int allGood = 0;
  while(!allGood && iter < maxIters)
  {
     std::vector<Region> regions;
     if (vtkSVGroupsSegmenter::GetSpecificRegions(pd, arrayName, regions, targetRegions) != SV_OK)
     {
       fprintf(stderr,"Couldn't get regions\n");
       return SV_ERROR;
     }
    int numRegions = regions.size();

    int maxCluster = -1;
    for (int i=0; i<numRegions; i++)
    {
      if (regions[i].IndexCluster > maxCluster)
        maxCluster = regions[i].IndexCluster;
    }
    maxCluster++;
    std::vector<int> testPatches(maxCluster, 0);
    for (int i=0; i<numRegions; i++)
      testPatches[regions[i].IndexCluster]++;

    vtkNew(vtkIdList, removePatchIds);
    removePatchIds->SetNumberOfIds(maxCluster);
    vtkNew(vtkIdList, minBadElements);
    minBadElements->SetNumberOfIds(maxCluster);
    for (int i=0; i<maxCluster; i++)
    {
      minBadElements->SetId(i, 1000000000000000000);
      removePatchIds->SetId(i, -1);
    }
    for (int i=0; i<numRegions; i++)
    {
      if (testPatches[regions[i].IndexCluster] > 1)
      {
        if (regions[i].NumberOfElements < minBadElements->GetId(regions[i].IndexCluster))
        {
          minBadElements->SetId(regions[i].IndexCluster, regions[i].NumberOfElements);
          removePatchIds->SetId(regions[i].IndexCluster, i);;
        }
      }
    }

    int regionCount=0;
    for (int i=0; i<numRegions; i++)
    {
      int regId = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(regions[i].Elements[0]);

      int badRegion=0;

      if (regions[i].CornerPoints.size() < 4)
        badRegion = 1;

      if (noEndRegions)
      {
        if (regId%6 == 4 || regId%5 == 5)
          badRegion = 1;
      }

      if (removePatchIds->GetId(regions[i].IndexCluster) == i)
        badRegion = 1;

      if (badRegion == 1)
      {
        fprintf(stdout,"FOUND A BAD REGION OF PATCH %d\n", regions[i].IndexCluster);
        vtkNew(vtkIdList, regionIds);
        std::vector<int> neighborRegions;
        int numRegions=0;
        for (int j=0; j<pd->GetNumberOfCells(); j++)
        {
          int cellValue = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(j);
          if (regionIds->IsId(cellValue) == -1)
          {
            regionIds->InsertNextId(cellValue);
            numRegions++;
            neighborRegions.push_back(0.0);
          }
        }

        allGood = 0;
        int maxNeighborRegion = -1;
        int needsFix = 0;
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
              if (neighborVals->GetId(k) != cellValue && targetRegions->IsId(neighborVals->GetId(k)) != -1)
                neighborRegions[regionIds->IsId(neighborVals->GetId(k))] += 1;
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
              maxNeighborRegion = regionIds->GetId(j);
            }
          }
        }

        if (needsFix)
        {
          for (int j=0; j<regions[i].NumberOfElements; j++)
          {
            int cellId = regions[i].Elements[j];
            pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(cellId, maxNeighborRegion);
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
    fprintf(stdout, "BOUNDARY FIX ITER: %d\n", iter);
    iter++;
  }

  return SV_OK;
}


int vtkSVGroupsSegmenter::Parameterize(vtkPolyData *fullMapPd)
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

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(surfacer->GetOutput());
  cleaner->Update();

  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(cleaner->GetOutput());
  triangulator->Update();

  vtkNew(vtkLinearSubdivisionFilter, subdivider);
  subdivider->SetInputData(triangulator->GetOutput());
  subdivider->SetNumberOfSubdivisions(4);
  subdivider->Update();
  polycubePd->DeepCopy(subdivider->GetOutput());
  fprintf(stdout,"JUST CHECK: %d\n", polycubePd->GetNumberOfPoints());

  int numPatches = patches.size();

  vtkNew(vtkAppendPolyData, appender);

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

    if (patches[i].CornerPoints.size()  > 6)
      fprintf(stdout,"FUCKING TELL ME!!!\n");
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

    //std::string filename2 = "/Users/adamupdegrove/Desktop/tmp/Boundary_"+std::to_string(patchId)+".vtp";
    //vtkSVIOUtils::WriteVTPFile(filename2, boundaryMapper->GetOutput());
    //std::string filename4 = "/Users/adamupdegrove/Desktop/tmp/Mapping_"+std::to_string(patchId)+".vtp";
    //vtkSVIOUtils::WriteVTPFile(filename4, mapper->GetOutput());

    appender->AddInputData(tmpPoly);
  }

  appender->Update();

  vtkNew(vtkCleanPolyData, polyCleaner);
  polyCleaner->SetInputData(appender->GetOutput());
  polyCleaner->SetTolerance(1.0e-6);
  polyCleaner->Update();

  vtkNew(vtkPolyData, tmpPd);
  tmpPd = polyCleaner->GetOutput();
  tmpPd->BuildLinks();

  if (tmpPd->GetNumberOfPoints() != this->WorkPd->GetNumberOfPoints() ||
      tmpPd->GetNumberOfCells() != this->WorkPd->GetNumberOfCells())
    fprintf(stderr,"SOMETING WENT WRONG\n");

  vtkNew(vtkPoints, fullMapPoints); fullMapPoints->SetNumberOfPoints(tmpPd->GetNumberOfPoints());
  vtkNew(vtkCellArray, fullMapCells);

  vtkDataArray *realPointIds = tmpPd->GetPointData()->GetArray("TmpInternalIds");
  vtkDataArray *realCellIds =  tmpPd->GetCellData()->GetArray("TmpInternalIds");
  for (int i=0; i<tmpPd->GetNumberOfPoints(); i++)
  {
    double pt[3];
    tmpPd->GetPoint(i, pt);
    int realPointId = realPointIds->GetTuple1(i);
    fullMapPoints->SetPoint(realPointId, pt);
  }
  for (int i=0; i<tmpPd->GetNumberOfCells(); i++)
  {
    int getCellId = realCellIds->LookupValue(i);
    vtkIdType npts, *pts;
    tmpPd->GetCellPoints(getCellId, npts, pts);

    vtkNew(vtkIdList, newPointIds);
    newPointIds->SetNumberOfIds(npts);
    for (int j=0; j<npts; j++)
      newPointIds->SetId(j, realPointIds->GetTuple1(pts[j]));

    fullMapCells->InsertNextCell(newPointIds);
  }

  fullMapPd->SetPoints(fullMapPoints);
  fullMapPd->SetPolys(fullMapCells);
  fullMapPd->BuildLinks();

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/Mapping_All.vtp";
  vtkSVIOUtils::WriteVTPFile(filename, polyCleaner->GetOutput());

  vtkNew(vtkPolyData, mappedPd);
  this->InterpolateMapOntoTarget(polycubePd, this->WorkPd, fullMapPd, mappedPd);

  std::string filename5 = "/Users/adamupdegrove/Desktop/tmp/Mapped_Out.vtp";
  vtkSVIOUtils::WriteVTPFile(filename5, mappedPd);


  return SV_OK;
}

int vtkSVGroupsSegmenter::FormParametricHexMesh(vtkStructuredGrid *paraHexMesh)
{
  // Extract surface of polycube
  vtkNew(vtkPolyData, polycubePd);
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(this->Polycube);
  surfacer->Update();

  polycubePd->DeepCopy(surfacer->GetOutput());
  vtkDataArray *localPtIds = polycubePd->GetPointData()->GetArray("LocalPointIds");

  int ptIds[8];
  double pts[8][3];
  for (int i=0; i<8; i++)
  {
    ptIds[i] = localPtIds->LookupValue(i);
    polycubePd->GetPoint(ptIds[i], pts[i]);
  }

  double w = vtkSVMathUtils::Distance(pts[1], pts[2]);
  double w_vec[3];
  vtkMath::Subtract(pts[1], pts[2], w_vec);
  vtkMath::Normalize(w_vec);

  double l = vtkSVMathUtils::Distance(pts[1], pts[0]);
  double l_vec[3];
  vtkMath::Subtract(pts[1], pts[0], l_vec);
  vtkMath::Normalize(l_vec);

  double h = vtkSVMathUtils::Distance(pts[0], pts[4]);
  double h_vec[3];
  vtkMath::Subtract(pts[0], pts[4], h_vec);
  vtkMath::Normalize(h_vec);

  int w_div = 11;
  int l_div = 101;
  int h_div = 11;

  double w_dist = w/(w_div-1);
  double l_dist = l/(l_div-1);
  double h_dist = h/(h_div-1);

  vtkNew(vtkPoints, paraPoints);
  int dim[3]; dim[0] = w_div; dim[1] = l_div; dim[2] = h_div;
  paraHexMesh->SetDimensions(dim);
  paraHexMesh->SetPoints(paraPoints);
  paraHexMesh->GetPoints()->SetNumberOfPoints(w_div*l_div*h_div);

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<l_div; j++)
    {
      for (int k=0; k<h_div; k++)
      {
        double x_vec[3], y_vec[3], z_vec[3];
        for (int l=0; l<3; l++)
        {
          x_vec[l] = w_vec[l]*i*w_dist;
          y_vec[l] = l_vec[l]*j*l_dist;
          z_vec[l] = h_vec[l]*k*h_dist;
        }
        double x_pt[3], y_pt[3], final_pt[3];
        vtkMath::Add(pts[7], x_vec, x_pt);
        vtkMath::Add(x_pt, y_vec, y_pt);
        vtkMath::Add(y_pt, z_vec, final_pt);

        int pos[3]; pos[0]= i; pos[1] = j; pos[2] = k;
        int pId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(pId, final_pt);
      }
    }
  }

  vtkNew(vtkAppendFilter, converter);
  converter->SetInputData(paraHexMesh);
  converter->Update();

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/TEST_PARA.vtu";
  vtkSVIOUtils::WriteVTUFile(filename, converter->GetOutput());

  return SV_OK;
}

int vtkSVGroupsSegmenter::MapVolume(vtkStructuredGrid *paraHexMesh,
                                    vtkPolyData *fullMapPd,
                                    vtkStructuredGrid *realHexMesh)
{
  vtkNew(vtkAppendFilter, converter);
  converter->SetInputData(paraHexMesh);
  converter->Update();

  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(converter->GetOutput());
  ider->SetIdsArrayName("TmpInternalIds");
  ider->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(ider->GetOutput());
  surfacer->Update();

  vtkNew(vtkPolyData, polycubePd);
  polycubePd->DeepCopy(surfacer->GetOutput());

  vtkNew(vtkPolyData, mappedPd);
  this->InterpolateMapOntoTarget(polycubePd, this->WorkPd, fullMapPd, mappedPd);

  std::string filenamer = "/Users/adamupdegrove/Desktop/tmp/DOUBLE_CHECK.vtp";
  vtkSVIOUtils::WriteVTPFile(filenamer, mappedPd);

  // Now lets try volume
  vtkDataArray *parameIds = ider->GetOutput()->GetPointData()->GetArray("TmpInternalIds");
  vtkDataArray *mappedIds = mappedPd->GetPointData()->GetArray("TmpInternalIds");

  int dim[3];
  paraHexMesh->GetDimensions(dim);
  realHexMesh->SetDimensions(dim);
  vtkNew(vtkPoints, realHexMeshPoints);
  realHexMesh->SetPoints(realHexMeshPoints);
  realHexMesh->GetPoints()->SetNumberOfPoints(dim[0]*dim[1]*dim[2]);

  int pos[3];
  double first_para_xpt[3], first_para_ypt[3], first_para_zpt[3];
  double last_para_xpt[3], last_para_ypt[3], last_para_zpt[3];
  double first_real_xpt[3], first_real_ypt[3], first_real_zpt[3];
  double last_real_xpt[3], last_real_ypt[3], last_real_zpt[3];

  for (int i=0; i<dim[0]; i++)
  {
    for (int j=0; j<dim[1]; j++)
    {
      for (int k=0; k<dim[2]; k++)
      {
        pos[0] = 0; pos[1] = j; pos[2] = k;
        int ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoint(ptId, first_para_xpt);
        int transId = parameIds->GetTuple1(ptId);
        int realId = mappedIds->LookupValue(transId);
        mappedPd->GetPoint(realId, first_real_xpt);

        pos[0] = dim[0]-1; pos[1] = j; pos[2] = k;
        ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoint(ptId, last_para_xpt);
        transId = parameIds->GetTuple1(ptId);
        realId = mappedIds->LookupValue(transId);
        mappedPd->GetPoint(realId, last_real_xpt);

        pos[0] = i; pos[1] = 0; pos[2] = k;
        ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoint(ptId, first_para_ypt);
        transId = parameIds->GetTuple1(ptId);
        realId = mappedIds->LookupValue(transId);
        mappedPd->GetPoint(realId, first_real_ypt);

        pos[0] = i; pos[1] = dim[1]-1; pos[2] = k;
        ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoint(ptId, last_para_ypt);
        transId = parameIds->GetTuple1(ptId);
        realId = mappedIds->LookupValue(transId);
        mappedPd->GetPoint(realId, last_real_ypt);

        pos[0] = i; pos[1] = j; pos[2] = 0;
        ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoint(ptId, first_para_zpt);
        transId = parameIds->GetTuple1(ptId);
        realId = mappedIds->LookupValue(transId);
        mappedPd->GetPoint(realId, first_real_zpt);

        pos[0] = i; pos[1] = j; pos[2] = dim[2]-1;
        ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoint(ptId, last_para_zpt);
        transId = parameIds->GetTuple1(ptId);
        realId = mappedIds->LookupValue(transId);
        mappedPd->GetPoint(realId, last_real_zpt);

        pos[0] = i; pos[1] = j; pos[2] = k;

        double para_pt[3];
        ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoint(ptId, para_pt);

        double real_xpt[3], real_ypt[3], real_zpt[3], real_pt[3];
        double para_x_dist = vtkSVMathUtils::Distance(para_pt, first_para_xpt)/vtkSVMathUtils::Distance(last_para_xpt, first_para_xpt);

        double xvec[3];
        vtkMath::Subtract(last_real_xpt, first_real_xpt, xvec);
        vtkMath::Normalize(xvec);
        vtkMath::MultiplyScalar(xvec, para_x_dist);
        double new_x[3];
        vtkMath::Add(first_real_xpt, xvec, new_x);

        for (int r=0; r<3; r++)
          real_xpt[r] = (1-para_x_dist) * first_real_xpt[r] + para_x_dist * last_real_xpt[r];

        double para_y_dist = vtkSVMathUtils::Distance(para_pt, first_para_ypt)/vtkSVMathUtils::Distance(last_para_ypt, first_para_ypt);

        for (int r=0; r<3; r++)
          real_ypt[r] = (1-para_y_dist) * first_real_ypt[r] + para_y_dist * last_real_ypt[r];

        double para_z_dist = vtkSVMathUtils::Distance(para_pt, first_para_zpt)/vtkSVMathUtils::Distance(last_para_zpt, first_para_zpt);

        double zvec[3];
        vtkMath::Subtract(last_real_zpt, first_real_zpt, zvec);
        vtkMath::Normalize(zvec);
        vtkMath::MultiplyScalar(zvec, para_z_dist);
        double new_z[3];
        vtkMath::Add(first_real_zpt, zvec, new_z);

        for (int r=0; r<3; r++)
          real_zpt[r] = (1-para_z_dist) * first_real_zpt[r] + para_z_dist * last_real_zpt[r];

        if (i == 0 || i == dim[0] - 1)
          vtkMath::Add(real_xpt, real_xpt, real_pt);
        else if (j == 0 || j == dim[1]-1)
          vtkMath::Add(real_ypt, real_ypt, real_pt);
        else if (k == 0 || k == dim[2]-1)
          vtkMath::Add(real_zpt, real_zpt, real_pt);
        else
          vtkMath::Add(real_xpt, real_zpt, real_pt);

        vtkMath::MultiplyScalar(real_pt, 1./2);

        realHexMesh->GetPoints()->SetPoint(ptId, real_pt);

        if (i == dim[0]-1 && j == 0 && k == dim[2]-1)
        {
          fprintf(stdout,"Start X Para Pt: %.4f %.4f %.4f\n", first_para_xpt[0], first_para_xpt[1], first_para_xpt[2]);
          fprintf(stdout,"Start Y Para Pt: %.4f %.4f %.4f\n", first_para_ypt[0], first_para_ypt[1], first_para_ypt[2]);
          fprintf(stdout,"Start Z Para Pt: %.4f %.4f %.4f\n", first_para_zpt[0], first_para_zpt[1], first_para_zpt[2]);
          fprintf(stdout,"Para x dist: %.4f\n", para_x_dist);
          fprintf(stdout,"Para y dist: %.4f\n", para_y_dist);
          fprintf(stdout,"Para z dist: %.4f\n", para_z_dist);
          double actualPt[3];
          mappedPd->GetPoint(ptId, actualPt);
          fprintf(stdout,"Real Pt: %.4f %.4f %.4f\n", actualPt[0], actualPt[1], actualPt[2]);
          fprintf(stdout,"X Pt: %.4f %.4f %.4f\n", new_x[0], new_x[1], new_x[2]);
          fprintf(stdout,"Z Pt: %.4f %.4f %.4f\n", new_z[0], new_z[1], new_z[2]);
        }
      }
    }
  }

  vtkNew(vtkAppendFilter, myConverter);
  myConverter->SetInputData(realHexMesh);
  myConverter->Update();

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/DUMB.vtu";
  vtkSVIOUtils::WriteVTUFile(filename, myConverter->GetOutput());


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
  int closeMatch = -1;
  for (int i=0; i<ps->GetNumberOfPoints(); i++)
  {
    vtkNew(vtkIdList, pointCellValues);
    vtkSVGeneralUtils::GetPointCellsValues(ps, arrayName.c_str(), i, pointCellValues);
    int prevNum = pointCellValues->GetNumberOfIds();
    pointCellValues->IntersectWith(matchingVals);

    if (pointCellValues->GetNumberOfIds() == matchingVals->GetNumberOfIds() &&
        prevNum == pointCellValues->GetNumberOfIds())
    {
      // We found it!
      returnPtId = i;
      return SV_OK;
    }
    else if (pointCellValues->GetNumberOfIds() == matchingVals->GetNumberOfIds())
      closeMatch = i;
    else if (prevNum == pointCellValues->GetNumberOfIds() && prevNum == 4)
      closeMatch = i;
  }
  if (closeMatch != -1)
    returnPtId = closeMatch;

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
