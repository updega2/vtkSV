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

#include "vtkSVCenterlinesEdgeWeightedCVT.h"
#include "vtkSVCleanUnstructuredGrid.h"
#include "vtkSVEdgeWeightedCVT.h"
#include "vtkSVEdgeWeightedSmoother.h"
#include "vtkSVFindGeodesicPath.h"
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
#include "vtkSVMUPFESNURBSWriter.h"

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
#include "vtkHexahedron.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
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

#include "vtkvmtkMergeCenterlines.h"

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
  this->MergedCenterlines = vtkPolyData::New();
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
  this->PolycubeUnitLength = 0.1;
}

// ----------------------
// Destructor
// ----------------------
vtkSVGroupsSegmenter::~vtkSVGroupsSegmenter()
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
  if (this->Centerlines != NULL)
  {
    this->Centerlines->Delete();
    this->Centerlines = NULL;
  }
  if (this->Polycube != NULL)
  {
    this->Polycube->Delete();
    this->Polycube = NULL;
  }
  if (this->GraphPd != NULL)
  {
    this->GraphPd->Delete();
    this->GraphPd = NULL;
  }

  if (this->CenterlineGroupIds != NULL)
  {
    this->CenterlineGroupIds->Delete();
    this->CenterlineGroupIds = NULL;
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

  if (this->CenterlineGraph != NULL)
  {
    this->CenterlineGraph->Delete();
    this->CenterlineGraph = NULL;
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

  if (this->MergeCenterlines() != SV_OK)
  {
    vtkErrorMacro("Problem merging centerlines");
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
  if (vtkSVGeneralUtils::CheckArrayExists(this->MergedCenterlines, 1, this->CenterlineGroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "CenterlineGroupIdsArray with name specified does not exist");
    return SV_OK;
  }

  //if (!this->BlankingArrayName)
  //{
  //  vtkDebugMacro("Blanking Array Name not given, setting to Blanking");
  //  this->BlankingArrayName = new char[strlen("Blanking") + 1];
  //  strcpy(this->BlankingArrayName, "Blanking");
  //}

  //if (vtkSVGeneralUtils::CheckArrayExists(this->Centerlines, 1, this->BlankingArrayName) != SV_OK)
  //{
  //  vtkErrorMacro(<< "BlankingArrayName with name specified does not exist");
  //  return SV_ERROR;
  //}

  //if (!this->CenterlineRadiusArrayName)
  //{
  //  vtkDebugMacro("Centerline radius Array Name not given, setting to MaximumInscribedSphereRadius");
  //  this->CenterlineRadiusArrayName = new char[strlen("MaximumInscribedSphereRadius") + 1];
  //  strcpy(this->CenterlineRadiusArrayName, "MaximumInscribedSphereRadius");
  //}

  //if (!this->Centerlines->GetPointData()->GetArray(this->CenterlineRadiusArrayName))
  //{
  //  vtkErrorMacro(<< "CenterlineRadiusArray with name specified does not exist");
  //  return SV_ERROR;
  //}

  this->CenterlineGraph = new vtkSVCenterlineGraph(0, this->MergedCenterlines,
                                                this->GroupIdsArrayName);

  if (this->CenterlineGraph->BuildGraph() != SV_OK)
  {
    vtkErrorMacro("Unable to form graph of centerlines");
    return SV_ERROR;
  }

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/CenterlineGraph.vtp";
  this->CenterlineGraph->GetGraphPolyData(this->GraphPd);
  vtkSVIOUtils::WriteVTPFile(filename, this->GraphPd);

  this->MergedCenterlines->DeepCopy(this->CenterlineGraph->Lines);
  std::string filename2 = "/Users/adamupdegrove/Desktop/tmp/CenterlineDirs.vtp";
  vtkSVIOUtils::WriteVTPFile(filename2, this->CenterlineGraph->Lines);

  this->CenterlineGraph->GetPolycube(10*this->PolycubeUnitLength, 10*this->PolycubeUnitLength, this->Polycube);
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
    this->MergedCenterlines->GetCellData()->GetArray(this->CenterlineGroupIdsArrayName);

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
  CVT->SetGenerators(this->MergedCenterlines);
  CVT->SetNumberOfRings(2);
  CVT->SetThreshold(stopCellNumber);
  CVT->SetUseCurvatureWeight(1);
  CVT->SetPatchIdsArrayName(this->GroupIdsArrayName);
  CVT->SetCVTDataArrayName("Normals");
  CVT->SetGroupIdsArrayName(this->GroupIdsArrayName);
  CVT->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
  CVT->SetBlankingArrayName(this->BlankingArrayName);
  CVT->SetUseRadiusInformation(this->UseRadiusInformation);
  CVT->SetUseBifurcationInformation(1);
  CVT->SetUsePointNormal(1);
  //CVT->SetUseRadiusInformation(0);
  //CVT->SetUseBifurcationInformation(0);
  //CVT->SetUsePointNormal(0);
  CVT->Update();

  vtkNew(vtkSVEdgeWeightedSmoother, smoother);
  smoother->SetInputData(CVT->GetOutput());
  smoother->SetGenerators(this->MergedCenterlines);
  smoother->SetNumberOfRings(2);
  smoother->SetThreshold(stopCellNumber);
  smoother->SetUseCurvatureWeight(0);
  smoother->SetNoInitialization(1);
  smoother->SetPatchIdsArrayName(this->GroupIdsArrayName);
  smoother->SetCVTDataArrayName("Normals");
  smoother->Update();

  this->WorkPd->DeepCopy(smoother->GetOutput());

  if (this->CheckGroups() != SV_OK)
  {
    vtkErrorMacro("Error in correcting groups");
    return SV_ERROR;
  }

  if (this->CorrectCellBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro("Could not correcto boundaries of surface");
    return SV_ERROR;
  }

  std::vector<Region> firstRegions;
  if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, firstRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get group regions");
    return SV_ERROR;
  }
  if (this->FixGroups(this->WorkPd, this->GroupIdsArrayName, firstRegions) != SV_OK)
  {
    vtkErrorMacro("Error in correcting groups");
    return SV_ERROR;
  }
  //if (this->FixGroupsWithPolycube() != SV_OK)
  //{
  //  vtkErrorMacro("Error in correcting groups");
  //  return SV_ERROR;
  //}

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

  if (this->SplitUpBoundaries(this->WorkPd, this->GroupIdsArrayName, groupRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't fix stuff\n");
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
    vtkSVGeneralUtils::ThresholdPd(this->MergedCenterlines, groupId, groupId, 1,
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
    //groupTubes->BuildLocator();

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

  vtkIntArray *tmpPatchArray = vtkIntArray::New();
  tmpPatchArray->SetNumberOfTuples(this->WorkPd->GetNumberOfCells());
  tmpPatchArray->SetName("PatchIds");
  tmpPatchArray->FillComponent(0, -1);
  this->WorkPd->GetCellData()->AddArray(tmpPatchArray);
  tmpPatchArray->Delete();

  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);
    vtkNew(vtkPolyData, branchPd);
    vtkSVGeneralUtils::ThresholdPd(this->WorkPd, groupId, groupId, 1,
      this->GroupIdsArrayName, branchPd);
    branchPd->BuildLinks();

    if (this->RunEdgeWeightedCVT(branchPd, generatorsPd) != SV_OK)
    {
      vtkErrorMacro("Error in cvt");
      return SV_ERROR;
    }

    if (this->FixEndPatches(branchPd) != SV_OK)
    {
      vtkErrorMacro("Error fixing end patches");
      return SV_ERROR;
    }

    if (this->FixSidePatches(branchPd) != SV_OK)
    {
      vtkErrorMacro("Error fixing side patches");
      return SV_ERROR;
    }

    vtkNew(vtkIdList, noEndPatches);
    noEndPatches->SetNumberOfIds(4);
    for (int j=0; j<4; j++)
      noEndPatches->SetId(j, j);

    if (this->CorrectSpecificCellBoundaries(branchPd, "PatchIds", noEndPatches) != SV_OK)
    {
      vtkErrorMacro("Could not correcto boundaries of surface");
      return SV_ERROR;
    }

    if (this->MergedCenterlines->GetNumberOfCells() > 1)
    {
      if (this->MatchEndPatches(branchPd) != SV_OK)
      {
        vtkErrorMacro("Error matching end patches");
        return SV_ERROR;
      }
    }

    // Set vals on work pd
    for (int j=0; j<branchPd->GetNumberOfCells(); j++)
    {
      //Get real cell id
      int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);

      // Get val
      int cellVal = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(j);

      // Set val
      this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(realCellId, cellVal);
    }
  }

  this->WorkPd->GetCellData()->RemoveArray("TmpInternalIds");
  this->WorkPd->GetPointData()->RemoveArray("TmpInternalIds");
  this->WorkPd->BuildLinks();

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

  std::vector<Region> finalRegions;
  vtkNew(vtkIdList, targetPatches);
  targetPatches->SetNumberOfIds(numGroups*4);
  for (int i=0; i<numGroups; i++)
  {
    for (int j=0; j<4; j++)
      targetPatches->SetId(4*i+j, 6*i+j);
  }
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
  if (this->GetSpecificRegions(this->WorkPd, "PatchIds", finalRegions, targetPatches) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }
  if (this->CurveFitBoundaries(this->WorkPd, "PatchIds", finalRegions) != SV_OK)
  {
    vtkErrorMacro("Could not curve fit boundaries of surface");
    return SV_ERROR;
  }

  // For checking purposes
  if (this->FixPatchesWithPolycube() != SV_OK)
  {
    fprintf(stderr,"Couldn't fix patches\n");
    return SV_ERROR;
  }

  //// NOW PARAMETERIZE!!, WIILL BE MOVED to vtkSVPolycubeParameterizer
  //// TODO: RENAME THIS CLASS TO vtkSVCenterlinesSegmenter

  vtkNew(vtkPolyData, fullMapPd);
  if (this->ParameterizeSurface(fullMapPd) != SV_OK)
  {
    fprintf(stderr,"WRONG\n");
    return SV_ERROR;
  }

  vtkNew(vtkUnstructuredGrid, loftedVolume);
  if (this->ParameterizeVolume(fullMapPd, loftedVolume) != SV_OK)
  {
    fprintf(stderr,"Failed doing volume stuffs\n");
    return SV_ERROR;
  }

  return SV_OK;
}

// ----------------------
// MergeCenterlines
// ----------------------
int vtkSVGroupsSegmenter::MergeCenterlines()
{
  fprintf(stdout,"Merging centerlines...\n");
  vtkNew(vtkvmtkMergeCenterlines, merger);
  merger->SetInputData(this->Centerlines);
  merger->SetRadiusArrayName(this->CenterlineRadiusArrayName);
  merger->SetGroupIdsArrayName(this->GroupIdsArrayName);
  merger->SetCenterlineIdsArrayName("CenterlineIds");
  merger->SetTractIdsArrayName("TractIds");
  merger->SetBlankingArrayName(this->BlankingArrayName);
  merger->SetResamplingStepLength(0.0);
  merger->SetMergeBlanked(1);
  merger->Update();

  vtkNew(vtkCleanPolyData, lineCleaner);
  lineCleaner->SetInputData(merger->GetOutput());
  lineCleaner->Update();

  this->MergedCenterlines->DeepCopy(lineCleaner->GetOutput());

  fprintf(stdout,"Merged\n");

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
    if (allRegions[i].CornerPoints.size() != allRegions[i].NumberOfCorners)
    {
      fprintf(stderr,"All corners not found on region\n");
      return SV_ERROR;
    }
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
    if (allRegions[i].CornerPoints.size() != allRegions[i].NumberOfCorners)
    {
      fprintf(stderr,"All corners not found on region\n");
      return SV_ERROR;
    }
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

  position = (position+1)%3;
  return pts[position];
}

// ----------------------
// GetCWPoint
// ----------------------
int vtkSVGroupsSegmenter::GetCWPoint(vtkPolyData *pd, const int pointId, const int cellId)
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

  position = (position+2)%3;
  return pts[position];
}

// ----------------------
// CheckCellValuesEdge
// ----------------------
int vtkSVGroupsSegmenter::CheckCellValuesEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1)
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
  return SV_OK;
}

// ----------------------
// SplineKnots
// ----------------------
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

// ----------------------
// SplineCurve
// ----------------------
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

// ----------------------
// FixEndPatches
// ----------------------
int vtkSVGroupsSegmenter::FixEndPatches(vtkPolyData *pd)
{
  vtkNew(vtkIdList, targetRegions);
  targetRegions->SetNumberOfIds(2);
  targetRegions->SetId(0, 4);
  targetRegions->SetId(1, 5);

  std::vector<Region> endRegions;
  if (this->GetSpecificRegions(pd, "PatchIds", endRegions, targetRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  std::vector<int> individualFix;
  std::vector<int> wholePatchFix;
  this->CheckEndPatches(pd, endRegions, individualFix, wholePatchFix);

  vtkDataArray *cellNormals = pd->GetCellData()->GetArray("Normals");

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
      int removeCellValue = pd->GetCellData()->GetArray(
        "PatchIds")->GetTuple1(cellId);
      if (compare <= 0.95)
      {
        vtkNew(vtkIdList, neighborValues);
        vtkSVGeneralUtils::GetNeighborsCellsValues(pd, "PatchIds",
                                                   cellId,
                                                   neighborValues);
        int newCellValue = -1;
        for (int k=0; k<neighborValues->GetNumberOfIds(); k++)
        {
          if (neighborValues->GetId(k) != removeCellValue)
          {
            newCellValue = neighborValues->GetId(k);
            break;
          }
        }
        if (newCellValue != -1)
          pd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellId, newCellValue);
      }
    }
  }

  if (wholePatchFix.size() > 0)
  {
    fprintf(stdout,"FIXING WHOLE END PATCHES\n");

    // Determine fix approach

    vtkNew(vtkIdList, sideTargetRegions);
    sideTargetRegions->SetNumberOfIds(4);
    sideTargetRegions->SetId(0, 0);
    sideTargetRegions->SetId(1, 1);
    sideTargetRegions->SetId(2, 2);
    sideTargetRegions->SetId(3, 3);

    std::vector<Region> sideRegions;
    if (this->GetSpecificRegions(pd, "PatchIds", sideRegions, sideTargetRegions) != SV_OK)
    {
      vtkErrorMacro("Couldn't get patches");
      return SV_ERROR;
    }

    std::vector<int> sidePatchFix;
    this->CheckSidePatches(pd, sideRegions, sidePatchFix);

    int fixStrategy = 0;

    if (sidePatchFix.size() == 1)
    {
      std::vector<Region> allRegions;
      if (this->GetRegions(pd, "PatchIds", allRegions) != SV_OK)
      {
        vtkErrorMacro("Couldn't get patches");
        return SV_ERROR;
      }

      for (int r=0; r<wholePatchFix.size(); r++)
      {
        int badPatch = wholePatchFix[r];

        int numTouchers=0;
        for (int i=0; i<sidePatchFix.size(); i++)
        {
          int edgeTouchCount=0;
          for (int j=0; j<allRegions.size(); j++)
          {
            if (allRegions[j].IndexCluster == sidePatchFix[i])
            {
              for (int k=0; k<allRegions[j].BoundaryEdges.size(); k++)
              {
                for (int l=0; l<allRegions[j].BoundaryEdges[k].size()-1; l++)
                {
                  int ptId0 = allRegions[j].BoundaryEdges[k][l];
                  int ptId1 = allRegions[j].BoundaryEdges[k][l+1];

                  vtkNew(vtkIdList, cellEdgeNeighbors);
                  pd->GetCellEdgeNeighbors(-1, ptId0, ptId1, cellEdgeNeighbors);

                  for (int m=0; m<cellEdgeNeighbors->GetNumberOfIds(); m++)
                  {
                    int cellId  = cellEdgeNeighbors->GetId(m);
                    int cellVal = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(cellId);

                    if (cellVal == endRegions[badPatch].IndexCluster)
                      edgeTouchCount++;
                  }
                }
              }
            }
          }
          if (edgeTouchCount > 0)
            numTouchers++;
        }

        if (numTouchers == 1)
          fixStrategy = 1;

      }
    }

    if (fixStrategy == 0)
    {
      vtkNew(vtkPolyData, workPdCopy);
      workPdCopy->DeepCopy(pd);

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

      for (int r=0; r<wholePatchFix.size(); r++)
      {
        int badPatch = wholePatchFix[r];

        for (int j=0; j<endRegions[badPatch].NumberOfElements; j++)
        {
          int cellId = endRegions[badPatch].Elements[j];
          int newCellValue = workPdCopy->GetCellData()->GetArray("PatchIds")->GetTuple1(cellId);
          pd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellId, newCellValue);
        }
      }
    }
    else if (fixStrategy = 1)
    {
      for (int r=0; r<wholePatchFix.size(); r++)
      {
        int badPatch = wholePatchFix[r];

        for (int j=0; j<endRegions[badPatch].NumberOfElements; j++)
        {
          int cellId = endRegions[badPatch].Elements[j];
          int newCellValue = sidePatchFix[0];
          pd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellId, newCellValue);
        }
      }
    }
  }


  return SV_OK;
}

// ----------------------
// CheckEndPatches
// ----------------------
int vtkSVGroupsSegmenter::CheckEndPatches(vtkPolyData *pd,
                                          std::vector<Region> endRegions,
                                          std::vector<int> &individualFix,
                                          std::vector<int> &wholePatchFix)
{
  // Need to fix both the random weird patch that goes over and the
  // entire groups that shouldn't be there
  vtkDataArray *cellNormals = pd->GetCellData()->GetArray("Normals");

  int numRegions = endRegions.size();

  for (int i=0; i<numRegions; i++)
  {
    if (endRegions[i].NumberOfCorners != 4)
    {
      wholePatchFix.push_back(i);
      continue;
    }
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
      if (compare > 0.99)
        numClose+=1;
    }

    fprintf(stdout,"EXPECTED: %d, CLOSE: %d\n", endRegions[i].NumberOfElements, numClose);
    if (numClose != endRegions[i].NumberOfElements)
    {
      if (numClose < (0.95)*endRegions[i].NumberOfElements)
        wholePatchFix.push_back(i);
      else
        individualFix.push_back(i);
    }

  }

  return SV_OK;
}

// ----------------------
// FixSidePatches
// ----------------------
int vtkSVGroupsSegmenter::FixSidePatches(vtkPolyData *pd)
{
  vtkNew(vtkIdList, targetRegions);
  targetRegions->SetNumberOfIds(4);
  targetRegions->SetId(0, 0);
  targetRegions->SetId(1, 1);
  targetRegions->SetId(2, 2);
  targetRegions->SetId(3, 3);

  std::vector<Region> sideRegions;
  if (this->GetSpecificRegions(pd, "PatchIds", sideRegions, targetRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  std::vector<int> wholePatchFix;
  this->CheckSidePatches(pd, sideRegions, wholePatchFix);

  vtkDataArray *cellNormals = pd->GetCellData()->GetArray("Normals");

  if (wholePatchFix.size() > 0)
  {
    fprintf(stdout,"FIXING WHOLE SIDE PATCHES\n");

    for (int i=0; i<wholePatchFix.size(); i++)
    {
      int maxPatch;
      int maxNumberOfElements = -1;
      for (int j=0; j<sideRegions.size(); j++)
      {
        if (wholePatchFix[i] == sideRegions[j].IndexCluster)
        {
          if (maxNumberOfElements < sideRegions[j].NumberOfElements)
          {
            maxPatch = j;
            maxNumberOfElements = sideRegions[j].NumberOfElements;
          }
        }
      }

      std::vector<int> minPatchFixes;
      for (int j=0; j<sideRegions.size(); j++)
      {
        if (wholePatchFix[i] == sideRegions[j].IndexCluster && j != maxPatch)
          minPatchFixes.push_back(j);
      }

      for (int r=0; r<minPatchFixes.size(); r++)
      {
        int minPatch = minPatchFixes[r];

        vtkNew(vtkIdList, patchIds);
        vtkNew(vtkIdList, patchCount);
        for (int j=0; j<sideRegions[minPatch].BoundaryEdges.size(); j++)
        {
          for (int k=0; k<sideRegions[minPatch].BoundaryEdges[j].size()-1; k++)
          {
            int ptId0 = sideRegions[minPatch].BoundaryEdges[j][k];
            int ptId1 = sideRegions[minPatch].BoundaryEdges[j][k+1];

            vtkNew(vtkIdList, cellEdgeNeighbors);
            pd->GetCellEdgeNeighbors(-1, ptId0, ptId1, cellEdgeNeighbors);

            for (int l=0; l<cellEdgeNeighbors->GetNumberOfIds(); l++)
            {
              int cellId  = cellEdgeNeighbors->GetId(l);
              int cellVal = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(cellId);

              if (cellVal != wholePatchFix[i])
              {
                int isId = patchIds->IsId(cellVal);
                if (isId == -1)
                {
                  patchIds->InsertNextId(cellVal);
                  patchCount->InsertNextId(1);
                }
                else
                  patchCount->SetId(isId, patchCount->GetId(isId)+1);
              }
            }
          }
        }

        int maxVal = -1;
        int maxPatchId;
        for (int j=0; j<patchIds->GetNumberOfIds(); j++)
        {
          if (patchCount->GetId(j) > maxVal)
          {
            maxPatchId = patchIds->GetId(j);
            maxVal = patchCount->GetId(j);
          }
        }

        for (int k=0; k<sideRegions[minPatch].Elements.size(); k++)
        {
          int cellId = sideRegions[minPatch].Elements[k];

          pd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellId, maxPatchId);
        }
      }
    }
  }


  return SV_OK;
}

// ----------------------
// CheckSidePatches
// ----------------------
int vtkSVGroupsSegmenter::CheckSidePatches(vtkPolyData *pd,
                                          std::vector<Region> sideRegions,
                                          std::vector<int> &wholePatchFix)
{
  int numRegions = sideRegions.size();

  vtkNew(vtkIdList, regionIds);
  vtkNew(vtkIdList, regionCount);
  for (int i=0; i<numRegions; i++)
  {
    int isId = regionIds->IsId(sideRegions[i].IndexCluster);
    if (isId == -1)
    {
      regionIds->InsertNextId(sideRegions[i].IndexCluster);
      regionCount->InsertNextId(1);
    }
    else
      regionCount->SetId(isId, regionCount->GetId(isId)+1);
  }

  if (numRegions == 4 && regionIds->GetNumberOfIds() == 4)
    return SV_OK;

  if (numRegions < 4)
  {
    fprintf(stderr,"THIS IS REALLY BIG PROBLEM! NOT SURE WHAT TO DO ABOUT THIS\n");
    return SV_ERROR;
  }

  for (int i=0; i<regionCount->GetNumberOfIds(); i++)
  {
    if (regionCount->GetId(i) > 1)
      wholePatchFix.push_back(regionIds->GetId(i));
  }

  return SV_OK;
}

// ----------------------
// CheckGroups
// ----------------------
int vtkSVGroupsSegmenter::CheckGroups()
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

// ----------------------
// FixPatchesWithPolycube
// ----------------------
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
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
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

            fprintf(stderr,"POLYCUBE CORNERS: ");
            for (int k=0; k<polycubePatches[i].CornerPoints.size(); k++)
              fprintf(stderr,"%d ", polycubePatches[i].CornerPoints[k]);
            fprintf(stderr,"\n");

            fprintf(stderr,"SURFACE CORNERS: ");
            for (int k=0; k<surfacePatches[j].CornerPoints.size(); k++)
              fprintf(stderr,"%d ", surfacePatches[j].CornerPoints[k]);
            fprintf(stderr,"\n");

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

// ----------------------
// ParameterizeSurface
// ----------------------
int vtkSVGroupsSegmenter::ParameterizeSurface(vtkPolyData *fullMapPd)
{
  std::vector<Region> patches;
  if (this->GetRegions(this->WorkPd, "PatchIds", patches) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  // Extract surface, triangulate, and subdivide polycube
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(this->Polycube);
  surfacer->Update();

  vtkNew(vtkPolyData, polySurfacePd);
  polySurfacePd->DeepCopy(surfacer->GetOutput());

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(surfacer->GetOutput());
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
  cleaner->Update();

  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(cleaner->GetOutput());
  triangulator->Update();

  vtkNew(vtkPolyData, polycubePd);
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
    vtkNew(vtkPolyData, rotPolycube);
    vtkNew(vtkMatrix4x4, rotMatrix0);
    vtkNew(vtkMatrix4x4, rotMatrix1);
    this->RotateGroupToGlobalAxis(polySurfacePd, groupId, this->GroupIdsArrayName, rotPolycube, rotMatrix0, rotMatrix1);

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
        fprintf(stdout,"COULD NOT FIND: ");
        for (int r=0; r<patchVals->GetNumberOfIds(); r++)
          fprintf(stdout,"%d ", patchVals->GetId(r));
        fprintf(stdout,"\n");
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
    fprintf(stderr,"SOMETHING WENT WRONG\n");

  vtkNew(vtkPoints, fullMapPoints); fullMapPoints->SetNumberOfPoints(tmpPd->GetNumberOfPoints());
  vtkNew(vtkCellArray, fullMapCells);

  vtkNew(vtkPointData, newPointData);
  vtkNew(vtkCellData, newCellData);
  newPointData->CopyAllocate(tmpPd->GetPointData(), tmpPd->GetNumberOfPoints());
  newCellData->CopyAllocate(tmpPd->GetCellData(), tmpPd->GetNumberOfCells());

  vtkDataArray *realPointIds = tmpPd->GetPointData()->GetArray("TmpInternalIds");
  vtkDataArray *realCellIds =  tmpPd->GetCellData()->GetArray("TmpInternalIds");
  for (int i=0; i<tmpPd->GetNumberOfPoints(); i++)
  {
    double pt[3];
    tmpPd->GetPoint(i, pt);
    int realPointId = realPointIds->GetTuple1(i);
    fullMapPoints->SetPoint(realPointId, pt);
    newPointData->CopyData(tmpPd->GetPointData(), i, realPointId);
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
    newCellData->CopyData(tmpPd->GetCellData(), getCellId, i);
  }

  fullMapPd->SetPoints(fullMapPoints);
  fullMapPd->SetPolys(fullMapCells);
  fullMapPd->GetPointData()->PassData(newPointData);
  fullMapPd->GetCellData()->PassData(newCellData);
  fullMapPd->BuildLinks();

  // all data on fullMapPd now
  std::string filename = "/Users/adamupdegrove/Desktop/tmp/Mapping_All.vtp";
  vtkSVIOUtils::WriteVTPFile(filename, fullMapPd);

  vtkNew(vtkPolyData, mappedPd);
  this->InterpolateMapOntoTarget(polycubePd, this->WorkPd, fullMapPd, mappedPd);

  std::string filename5 = "/Users/adamupdegrove/Desktop/tmp/Mapped_Out.vtp";
  vtkSVIOUtils::WriteVTPFile(filename5, mappedPd);

  return SV_OK;
}

// ----------------------
// ParameterizeVolume
// ----------------------
int vtkSVGroupsSegmenter::ParameterizeVolume(vtkPolyData *fullMapPd, vtkUnstructuredGrid *loftedVolume)
{
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

  vtkNew(vtkAppendFilter, appender);

  std::vector<vtkSmartPointer<vtkStructuredGrid> > paraHexVolumes(numGroups);

  int w_div = 60*this->PolycubeUnitLength + 1;
  int h_div = 60*this->PolycubeUnitLength + 1;
  int l_div = 0; // Determined by length of cube

  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);
    vtkNew(vtkPolyData, branchMapPd);
    vtkSVGeneralUtils::ThresholdPd(fullMapPd, groupId, groupId, 1, this->GroupIdsArrayName, branchMapPd);

    vtkNew(vtkPolyData, branchPd);
    vtkSVGeneralUtils::ThresholdPd(this->WorkPd, groupId, groupId, 1, this->GroupIdsArrayName, branchPd);

    // Extract surface of polycube
    vtkNew(vtkPolyData, polycubePd);
    vtkNew(vtkDataSetSurfaceFilter, surfacer);
    surfacer->SetInputData(this->Polycube);
    surfacer->Update();

    vtkNew(vtkCleanPolyData, cubeCleaner);
    cubeCleaner->SetInputData(surfacer->GetOutput());
    cubeCleaner->ToleranceIsAbsoluteOn();
    cubeCleaner->SetAbsoluteTolerance(1.0e-6);
    cubeCleaner->Update();

    polycubePd->DeepCopy(cubeCleaner->GetOutput());

    vtkNew(vtkPolyData, branchPolycube);
    vtkSVGeneralUtils::ThresholdPd(polycubePd, groupId, groupId, 1, this->GroupIdsArrayName, branchPolycube);

    branchPolycube->BuildLinks();

    vtkNew(vtkStructuredGrid, paraHexMesh);
    if (this->FormParametricHexMesh(branchPolycube, paraHexMesh, w_div,
                                    l_div, h_div) != SV_OK)
    {
      fprintf(stderr,"Couldn't do the dirt\n");
      return SV_ERROR;
    }
    fprintf(stdout,"WHAT IS MY L_DVE: %d\n", l_div);
    vtkNew(vtkIntArray, groupIdsArray);
    groupIdsArray->SetNumberOfTuples(paraHexMesh->GetNumberOfCells());
    groupIdsArray->SetName(this->GroupIdsArrayName);
    groupIdsArray->FillComponent(0, groupId);

    paraHexMesh->GetCellData()->AddArray(groupIdsArray);

    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(paraHexMesh);
    ider->SetIdsArrayName("TmpInternalIds");
    ider->Update();

    appender->AddInputData(ider->GetOutput());
    paraHexVolumes[i] = vtkSmartPointer<vtkStructuredGrid>::New();
    paraHexVolumes[i]->DeepCopy(ider->GetOutput());
  }
  appender->Update();

  std::string filename = "/Users/adamupdegrove/Desktop/tmp/TEST_PARA.vtu";
  vtkSVIOUtils::WriteVTUFile(filename, appender->GetOutput());

  vtkNew(vtkUnstructuredGrid, paraHexVolume);
  paraHexVolume->DeepCopy(appender->GetOutput());

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(paraHexVolume);
  surfacer->Update();

  vtkNew(vtkPolyData, paraHexSurface);
  paraHexSurface->DeepCopy(surfacer->GetOutput());

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(paraHexSurface);
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
  cleaner->Update();

  vtkNew(vtkPolyData, paraHexCleanSurface);
  paraHexCleanSurface->DeepCopy(cleaner->GetOutput());

  vtkNew(vtkSVCleanUnstructuredGrid, cleaner2);
  cleaner2->SetInputData(paraHexVolume);
  cleaner2->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer2);
  surfacer2->SetInputData(cleaner2->GetOutput());
  surfacer2->Update();

  vtkNew(vtkPolyData, cleanSurface);
  cleanSurface->DeepCopy(surfacer2->GetOutput());

  std::vector<int> surfacePtMap;
  std::vector<std::vector<int> > invPtMap;
  this->GetInteriorPointMaps(paraHexSurface, paraHexCleanSurface, cleanSurface, surfacePtMap, invPtMap);

  vtkNew(vtkPolyData, mappedSurface);
  this->InterpolateMapOntoTarget(paraHexSurface, this->WorkPd, fullMapPd, mappedSurface);

  vtkNew(vtkIdFilter, ider2);
  ider2->SetInputData(mappedSurface);
  ider2->SetIdsArrayName("TmpInternalIds2");
  ider2->Update();
  vtkDataArray *tmpArray = ider2->GetOutput()->GetPointData()->GetArray("TmpInternalIds2");
  mappedSurface->GetPointData()->AddArray(tmpArray);

  vtkNew(vtkAppendFilter, surfaceAppender);
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);

    vtkNew(vtkPolyData, mappedBranch);
    vtkSVGeneralUtils::ThresholdPd(mappedSurface, groupId, groupId, 1, this->GroupIdsArrayName, mappedBranch);

    if (this->MapInteriorBoundary(paraHexVolumes[i], mappedBranch, surfacePtMap) != SV_OK)
    {
      fprintf(stderr,"Couldn't do the dirt\n");
      return SV_ERROR;
    }
    surfaceAppender->AddInputData(mappedBranch);
  }

  surfaceAppender->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer3);
  surfacer3->SetInputData(surfaceAppender->GetOutput());
  surfacer3->Update();
  mappedSurface->DeepCopy(surfacer3->GetOutput());

  this->FixInteriorBoundary(mappedSurface, invPtMap);

  std::string filenamer = "/Users/adamupdegrove/Desktop/tmp/DOUBLE_CHECK.vtp";
  vtkSVIOUtils::WriteVTPFile(filenamer, mappedSurface);

  vtkNew(vtkAppendFilter, volumeAppender);
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);

    vtkNew(vtkPolyData, mappedBranch);
    vtkSVGeneralUtils::ThresholdPd(mappedSurface, groupId, groupId, 1, this->GroupIdsArrayName, mappedBranch);

    vtkNew(vtkStructuredGrid, realHexMesh);
    if (this->MapVolume(paraHexVolumes[i], mappedBranch, realHexMesh) != SV_OK)
    {
      fprintf(stderr,"Couldn't do the dirt\n");
      return SV_ERROR;
    }

    vtkNew(vtkIntArray, groupIdsArray);
    groupIdsArray->SetNumberOfTuples(realHexMesh->GetNumberOfCells());
    groupIdsArray->SetName(this->GroupIdsArrayName);
    groupIdsArray->FillComponent(0, groupId);
    realHexMesh->GetCellData()->AddArray(groupIdsArray);

    vtkNew(vtkIdFilter, ider3);
    ider3->SetInputData(realHexMesh);
    ider3->SetIdsArrayName("TmpInternalIds");
    ider3->Update();

    volumeAppender->AddInputData(ider3->GetOutput());
  }

  volumeAppender->Update();
  vtkNew(vtkUnstructuredGrid, mappedVolume);
  mappedVolume->DeepCopy(volumeAppender->GetOutput());

  vtkNew(vtkSVCleanUnstructuredGrid, cleaner3);
  cleaner3->SetInputData(volumeAppender->GetOutput());
  cleaner3->Update();

  vtkNew(vtkUnstructuredGrid, smoothVolume);
  smoothVolume->DeepCopy(cleaner3->GetOutput());

  std::vector<int> volumePtMap;
  this->GetVolumePointMap(mappedVolume, smoothVolume, volumePtMap);

  int smoothIters = 500;
  if (this->SmoothUnstructuredGrid(smoothVolume, smoothIters) != SV_OK)
  {
    fprintf(stderr,"Couldn't smooth volume\n");
    return SV_ERROR;
  }

  this->FixVolume(mappedVolume, smoothVolume, volumePtMap);

  filename = "/Users/adamupdegrove/Desktop/tmp/TEST_FINAL.vtu";
  vtkSVIOUtils::WriteVTUFile(filename, mappedVolume);

  vtkNew(vtkAppendFilter, loftAppender);
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);

    vtkNew(vtkUnstructuredGrid, mappedBranch);
    vtkSVGeneralUtils::ThresholdUg(mappedVolume, groupId, groupId, 1, this->GroupIdsArrayName, mappedBranch);

    vtkNew(vtkStructuredGrid, realHexMesh);
    if (this->ConvertUGToSG(mappedBranch, realHexMesh, w_div, l_div, h_div) != SV_OK)
    {
      fprintf(stderr,"Couldn't do the dirt\n");
      return SV_ERROR;
    }

    // Set up the volume
    vtkNew(vtkUnstructuredGrid, emptyGrid);
    vtkNew(vtkSVLoftNURBSVolume, lofter);
    lofter->SetInputData(emptyGrid);
    lofter->SetInputGrid(realHexMesh);
    lofter->SetUDegree(2);
    lofter->SetVDegree(2);
    lofter->SetWDegree(2);
    lofter->SetUnstructuredGridUSpacing(1./w_div);
    lofter->SetUnstructuredGridVSpacing(1./l_div);
    lofter->SetUnstructuredGridWSpacing(1./h_div);
    lofter->SetUKnotSpanType("average");
    lofter->SetUParametricSpanType("chord");
    lofter->SetVKnotSpanType("average");
    lofter->SetVParametricSpanType("chord");
    lofter->SetWKnotSpanType("average");
    lofter->SetWParametricSpanType("chord");
    lofter->Update();

    loftAppender->AddInputData(lofter->GetOutput());

    if (this->MergedCenterlines->GetNumberOfCells() == 1)
    {
      std::string mfsname = "/Users/adamupdegrove/Desktop/tmp/Pipe.msh";
      vtkNew(vtkSVMUPFESNURBSWriter, writer);
      writer->SetInputData(lofter->GetVolume());
      writer->SetFileName(mfsname.c_str());
      writer->Write();
    }
  }

  loftAppender->Update();
  loftedVolume->DeepCopy(loftAppender->GetOutput());

  return SV_OK;
}

// ----------------------
// ConvertUGToSG
// ----------------------
int vtkSVGroupsSegmenter::ConvertUGToSG(vtkUnstructuredGrid *ug,
                                        vtkStructuredGrid *sg,
                                        const int w_div, const int l_div,
                                        const int h_div)
{
  vtkDataArray *ptIds = ug->GetPointData()->GetArray("TmpInternalIds");

  int dim[3]; dim[0] = w_div; dim[1] = l_div; dim[2] = h_div;

  vtkNew(vtkPoints, sgPoints);
  sg->SetPoints(sgPoints);
  sg->GetPoints()->SetNumberOfPoints(dim[0]*dim[1]*dim[2]);
  sg->SetDimensions(dim);

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<l_div; j++)
    {
      for (int k=0; k<h_div; k++)
      {
        int pos[3]; pos[0] = i; pos[1] = j; pos[2] = k;
        int ptId = vtkStructuredData::ComputePointId(dim, pos);

        int realId = ptIds->LookupValue(ptId);

        double pt[3];
        ug->GetPoint(realId, pt);

        sg->GetPoints()->SetPoint(ptId, pt);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// FormParametricHexMesh
// ----------------------
int vtkSVGroupsSegmenter::FormParametricHexMesh(vtkPolyData *polycubePd, vtkStructuredGrid *paraHexMesh,
                                                const int w_div, int &l_div, const int h_div)
{
  // GetFace 0, right side face
  vtkIdType f0npts, *f0PtIds;
  polycubePd->GetCellPoints(0, f0npts, f0PtIds);

  // GetFace 1, face underneath
  vtkIdType f1npts, *f1PtIds;
  polycubePd->GetCellPoints(1, f1npts, f1PtIds);

  // GetFace 2, left side face
  vtkIdType f2npts, *f2PtIds;
  polycubePd->GetCellPoints(2, f2npts, f2PtIds);

  // GetFace 3, face on top
  vtkIdType f3npts, *f3PtIds;
  polycubePd->GetCellPoints(3, f3npts, f3PtIds);

  double pt0[3], pt1[3], pt2[3], pt3[3];

  // bottom lower left corner
  polycubePd->GetPoint(f1PtIds[0], pt0);

  // bottom lower right corner
  polycubePd->GetPoint(f0PtIds[0], pt1);

  // bottom upper right corner
  polycubePd->GetPoint(f0PtIds[1], pt2);

  // top lower right corner
  polycubePd->GetPoint(f2PtIds[0], pt3);

  // Set boundary
  double w = vtkSVMathUtils::Distance(pt1, pt0);
  double w_vec[3];
  vtkMath::Subtract(pt1, pt0, w_vec);
  vtkMath::Normalize(w_vec);

  double l = vtkSVMathUtils::Distance(pt2, pt1);
  double l_vec[3];
  vtkMath::Subtract(pt2, pt1, l_vec);
  vtkMath::Normalize(l_vec);

  double h = vtkSVMathUtils::Distance(pt3, pt0);
  double h_vec[3];
  vtkMath::Subtract(pt3, pt0, h_vec);
  vtkMath::Normalize(h_vec);

  // If top cube, need to use top of cube, not bottom
  if (f0npts == 4 && f1npts == 5)
  {
    double pt4[3];
    polycubePd->GetPoint(f1PtIds[1], pt4);

    double pt5[3];
    polycubePd->GetPoint(f1PtIds[2], pt5);

    double tmpVec[3];
    vtkMath::Subtract(pt5, pt4, tmpVec);
    vtkMath::Normalize(tmpVec);

    double test_dot = vtkMath::Dot(tmpVec, l_vec);
    if (test_dot < 1.0e-6 && test_dot > -1.0e-6)
    {
      w = vtkSVMathUtils::Distance(pt5, pt4);
      vtkMath::Subtract(pt5, pt4, w_vec);
      vtkMath::Normalize(w_vec);

      l = vtkSVMathUtils::Distance(pt4, pt0);
      vtkMath::Subtract(pt4, pt0, l_vec);
      vtkMath::Normalize(l_vec);
    }
  }

  l_div = floor(l/(10.0*this->PolycubeUnitLength));
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
        vtkMath::Add(pt0, x_vec, x_pt);
        vtkMath::Add(x_pt, y_vec, y_pt);
        vtkMath::Add(y_pt, z_vec, final_pt);

        int pos[3]; pos[0]= i; pos[1] = j; pos[2] = k;
        int pId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(pId, final_pt);
      }
    }
  }

  if (f0npts==4)
  {
    if (f1npts==5)
    {
      double pt4[3];
      polycubePd->GetPoint(f1PtIds[1], pt4);

      double pt5[3];
      polycubePd->GetPoint(f1PtIds[2], pt5);

      double tmpVec[3];
      vtkMath::Subtract(pt5, pt4, tmpVec);
      vtkMath::Normalize(tmpVec);

      double test_dot = vtkMath::Dot(tmpVec, l_vec);
      if (test_dot < 1.0e-6 && test_dot > -1.0e-6)
      {
        double bottomLeftPt[3];
        polycubePd->GetPoint(f1PtIds[0], bottomLeftPt);

        double bottomRightPt[3];
        polycubePd->GetPoint(f1PtIds[3], bottomRightPt);

        double bottomMidPt0[3];
        polycubePd->GetPoint(f1PtIds[4], bottomMidPt0);

        double bottomMidPt1[3];
        polycubePd->GetPoint(f3PtIds[4], bottomMidPt1);

        this->PushStructuredGridZAxis(paraHexMesh, bottomLeftPt, bottomRightPt, bottomMidPt0, bottomMidPt1, 1);
      }
      else
      {
        double topLeftPt[3];
        polycubePd->GetPoint(f1PtIds[1], topLeftPt);

        double topRightPt[3];
        polycubePd->GetPoint(f1PtIds[3], topRightPt);

        double topMidPt0[3];
        polycubePd->GetPoint(f1PtIds[2], topMidPt0);

        double topMidPt1[3];
        polycubePd->GetPoint(f3PtIds[2], topMidPt1);

        this->PushStructuredGridZAxis(paraHexMesh, topLeftPt, topRightPt, topMidPt0, topMidPt1, 0);
      }
    }
    else if (f1npts==6)
    {
      double topLeftPt[3];
      polycubePd->GetPoint(f1PtIds[1], topLeftPt);

      double topRightPt[3];
      polycubePd->GetPoint(f1PtIds[3], topRightPt);

      double topMidPt0[3];
      polycubePd->GetPoint(f1PtIds[2], topMidPt0);

      double topMidPt1[3];
      polycubePd->GetPoint(f3PtIds[2], topMidPt1);

      this->PushStructuredGridZAxis(paraHexMesh, topLeftPt, topRightPt, topMidPt0, topMidPt1, 0);

      double bottomLeftPt[3];
      polycubePd->GetPoint(f1PtIds[0], bottomLeftPt);

      double bottomRightPt[3];
      polycubePd->GetPoint(f1PtIds[4], bottomRightPt);

      double bottomMidPt0[3];
      polycubePd->GetPoint(f1PtIds[5], bottomMidPt0);

      double bottomMidPt1[3];
      polycubePd->GetPoint(f3PtIds[5], bottomMidPt1);

      this->PushStructuredGridZAxis(paraHexMesh, bottomLeftPt, bottomRightPt, bottomMidPt0, bottomMidPt1, 1);
    }
  }
  else if (f0npts==5)
  {
    if (f1npts==4)
    {
      double sideBottomPt[3];
      polycubePd->GetPoint(f2PtIds[3], sideBottomPt);

      double sideTopPt[3];
      polycubePd->GetPoint(f2PtIds[1], sideTopPt);

      double sideMidPt0[3];
      polycubePd->GetPoint(f2PtIds[2], sideMidPt0);

      double sideMidPt1[3];
      polycubePd->GetPoint(f0PtIds[2], sideMidPt1);

      this->PushStructuredGridXAxis(paraHexMesh, sideBottomPt, sideTopPt, sideMidPt0, sideMidPt1, 0);
    }
    else if (f1npts==5)
    {
      double pt4[3];
      polycubePd->GetPoint(f1PtIds[1], pt4);

      double pt5[3];
      polycubePd->GetPoint(f1PtIds[2], pt5);

      double tmpVec[3];
      vtkMath::Subtract(pt5, pt4, tmpVec);
      vtkMath::Normalize(tmpVec);

      double test_dot = vtkMath::Dot(tmpVec, l_vec);
      if (test_dot < 1.0e-6 && test_dot > -1.0e-6)
      {
        double bottomLeftPt[3];
        polycubePd->GetPoint(f1PtIds[0], bottomLeftPt);

        double bottomRightPt[3];
        polycubePd->GetPoint(f1PtIds[3], bottomRightPt);

        double bottomMidPt0[3];
        polycubePd->GetPoint(f1PtIds[4], bottomMidPt0);

        double bottomMidPt1[3];
        polycubePd->GetPoint(f3PtIds[4], bottomMidPt1);

        this->PushStructuredGridZAxis(paraHexMesh, bottomLeftPt, bottomRightPt, bottomMidPt0, bottomMidPt1, 1);

        double sideBottomPt[3];
        polycubePd->GetPoint(f2PtIds[3], sideBottomPt);

        double sideTopPt[3];
        polycubePd->GetPoint(f2PtIds[1], sideTopPt);

        double sideMidPt0[3];
        polycubePd->GetPoint(f2PtIds[2], sideMidPt0);

        double sideMidPt1[3];
        polycubePd->GetPoint(f0PtIds[2], sideMidPt1);

        this->PushStructuredGridXAxis(paraHexMesh, sideBottomPt, sideTopPt, sideMidPt0, sideMidPt1, 0);
      }
      else
      {
        double topLeftPt[3];
        polycubePd->GetPoint(f1PtIds[1], topLeftPt);

        double topRightPt[3];
        polycubePd->GetPoint(f1PtIds[3], topRightPt);

        double topMidPt0[3];
        polycubePd->GetPoint(f1PtIds[2], topMidPt0);

        double topMidPt1[3];
        polycubePd->GetPoint(f3PtIds[2], topMidPt1);

        this->PushStructuredGridZAxis(paraHexMesh, topLeftPt, topRightPt, topMidPt0, topMidPt1, 0);

        double sideBottomPt[3];
        polycubePd->GetPoint(f2PtIds[3], sideBottomPt);

        double sideTopPt[3];
        polycubePd->GetPoint(f2PtIds[0], sideTopPt);

        double sideMidPt0[3];
        polycubePd->GetPoint(f2PtIds[4], sideMidPt0);

        double sideMidPt1[3];
        polycubePd->GetPoint(f0PtIds[4], sideMidPt1);

        this->PushStructuredGridXAxis(paraHexMesh, sideBottomPt, sideTopPt, sideMidPt0, sideMidPt1, 1);
      }
    }
  }
  else if (f0npts==6)
  {
    if (f1npts==4)
    {
      double sideBottomPt0[3];
      polycubePd->GetPoint(f2PtIds[3], sideBottomPt0);

      double sideTopPt0[3];
      polycubePd->GetPoint(f2PtIds[1], sideTopPt0);

      double sideMidPt0[3];
      polycubePd->GetPoint(f2PtIds[2], sideMidPt0);

      double sideMidPt1[3];
      polycubePd->GetPoint(f0PtIds[2], sideMidPt1);

      this->PushStructuredGridXAxis(paraHexMesh, sideBottomPt0, sideTopPt0, sideMidPt0, sideMidPt1, 0);

      double sideBottomPt1[3];
      polycubePd->GetPoint(f2PtIds[4], sideBottomPt1);

      double sideTopPt1[3];
      polycubePd->GetPoint(f2PtIds[0], sideTopPt1);

      double sideMidPt2[3];
      polycubePd->GetPoint(f2PtIds[5], sideMidPt2);

      double sideMidPt3[3];
      polycubePd->GetPoint(f0PtIds[5], sideMidPt3);

      this->PushStructuredGridXAxis(paraHexMesh, sideBottomPt1, sideTopPt1, sideMidPt2, sideMidPt3, 1);
    }
  }

  return SV_OK;
}

// ----------------------
// PushStructuredGridZAxis
// ----------------------
int vtkSVGroupsSegmenter::PushStructuredGridZAxis(vtkStructuredGrid *paraHexMesh,
                                                  const double pt0[3],
                                                  const double pt1[3],
                                                  const double pt2[3],
                                                  const double pt3[3],
                                                  const int isBottom)
{
  int dim[3];
  paraHexMesh->GetDimensions(dim);

  int h_div = dim[2];
  double h = vtkSVMathUtils::Distance(pt3, pt2);
  double h_dist = h/(h_div-1);

  double h_vec[3];
  vtkMath::Subtract(pt3, pt2, h_vec);
  vtkMath::Normalize(h_vec);

  int w_div = dim[0];
  int half_w_div = floor(w_div/2.0);

  // Set new bottom extend points in the middle
  int pos[3];
  pos[0]= half_w_div;

  if (isBottom)
    pos[1] = 0;
  else
    pos[1] = dim[1]-1;

  for (int i=0; i<h_div; i++)
  {
    pos[2] = i;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);

    double z_vec[3];
    for (int j=0; j<3; j++)
      z_vec[j] = h_vec[j]*i*h_dist;

    double newPt[3];
    vtkMath::Add(pt2, z_vec, newPt);

    // New point in middle
    paraHexMesh->GetPoints()->SetPoint(ptId, newPt);

    // New point on 0 side
    double firstPt[3];
    vtkMath::Add(pt0, z_vec, firstPt);

    pos[0] = 0;
    int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->SetPoint(firstPtId, firstPt);

    double first_w = vtkSVMathUtils::Distance(newPt, firstPt);
    double first_w_dist = first_w/(half_w_div);

    double first_w_vec[3];
    vtkMath::Subtract(newPt, firstPt, first_w_vec);
    vtkMath::Normalize(first_w_vec);

    for (int j=0; j<half_w_div; j++)
    {
      pos[0] = j;

      double x_vec[3];
      for (int k=0; k<3; k++)
        x_vec[k] = first_w_vec[k]*j*first_w_dist;

      double inPt[3];
      vtkMath::Add(firstPt, x_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }

    //New point on end side
    double lastPt[3];
    vtkMath::Add(pt1, z_vec, lastPt);

    pos[0] = dim[0]-1;
    int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->SetPoint(lastPtId, lastPt);

    double last_w = vtkSVMathUtils::Distance(lastPt, newPt);
    double last_w_dist = last_w/(half_w_div);

    double last_w_vec[3];
    vtkMath::Subtract(lastPt, newPt, last_w_vec);
    vtkMath::Normalize(last_w_vec);

    for (int j=0; j<half_w_div; j++)
    {
      pos[0] = half_w_div+j;

      double x_vec[3];
      for (int k=0; k<3; k++)
        x_vec[k] = last_w_vec[k]*j*last_w_dist;

      double inPt[3];
      vtkMath::Add(newPt, x_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }
  }

  int l_div = dim[1];

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      pos[0] = i; pos[2] = j;

      pos[1] = 0;
      int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

      double firstPt[3];
      paraHexMesh->GetPoint(firstPtId, firstPt);

      pos[1] = dim[1]-1;
      int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

      double lastPt[3];
      paraHexMesh->GetPoint(lastPtId, lastPt);

      double l = vtkSVMathUtils::Distance(lastPt, firstPt);
      double l_dist = l/(l_div-1);

      double l_vec[3];
      vtkMath::Subtract(lastPt, firstPt, l_vec);
      vtkMath::Normalize(l_vec);

      for (int k=0; k<l_div; k++)
      {
        pos[1] = k;

        double y_vec[3];
        for (int l=0; l<3; l++)
          y_vec[l] = l_vec[l]*k*l_dist;

        double inPt[3];
        vtkMath::Add(firstPt, y_vec, inPt);

        int newPtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// PushStructuredGridXAxis
// ----------------------
int vtkSVGroupsSegmenter::PushStructuredGridXAxis(vtkStructuredGrid *paraHexMesh,
                                                    const double pt0[3],
                                                    const double pt1[3],
                                                    const double pt2[3],
                                                    const double pt3[3],
                                                    const int isBottom)
{
  int dim[3];
  paraHexMesh->GetDimensions(dim);

  int w_div = dim[0];
  double w = vtkSVMathUtils::Distance(pt3, pt2);
  double w_dist = w/(w_div-1);

  double w_vec[3];
  vtkMath::Subtract(pt3, pt2, w_vec);
  vtkMath::Normalize(w_vec);

  int h_div = dim[2];
  int half_h_div = floor(h_div/2.0);

  // Set new bottom extend points in the middle
  int pos[3];
  pos[2]= half_h_div;

  if (isBottom)
    pos[1] = 0;
  else
    pos[1] = dim[1]-1;

  for (int i=0; i<w_div; i++)
  {
    pos[0] = i;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);

    double x_vec[3];
    for (int j=0; j<3; j++)
      x_vec[j] = w_vec[j]*i*w_dist;

    double newPt[3];
    vtkMath::Add(pt2, x_vec, newPt);

    // New point in middle
    paraHexMesh->GetPoints()->SetPoint(ptId, newPt);

    // New point on 0 side
    double firstPt[3];
    vtkMath::Add(pt0, x_vec, firstPt);

    pos[2] = 0;
    int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->SetPoint(firstPtId, firstPt);

    double first_h = vtkSVMathUtils::Distance(newPt, firstPt);
    double first_h_dist = first_h/(half_h_div);

    double first_h_vec[3];
    vtkMath::Subtract(newPt, firstPt, first_h_vec);
    vtkMath::Normalize(first_h_vec);

    for (int j=0; j<half_h_div; j++)
    {
      pos[2] = j;

      double z_vec[3];
      for (int k=0; k<3; k++)
        z_vec[k] = first_h_vec[k]*j*first_h_dist;

      double inPt[3];
      vtkMath::Add(firstPt, z_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }

    //New point on end side
    double lastPt[3];
    vtkMath::Add(pt1, x_vec, lastPt);

    pos[2] = dim[2]-1;
    int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->SetPoint(lastPtId, lastPt);

    double last_h = vtkSVMathUtils::Distance(lastPt, newPt);
    double last_h_dist = last_h/(half_h_div);

    double last_h_vec[3];
    vtkMath::Subtract(lastPt, newPt, last_h_vec);
    vtkMath::Normalize(last_h_vec);

    for (int j=0; j<half_h_div; j++)
    {
      pos[2] = half_h_div+j;

      double z_vec[3];
      for (int k=0; k<3; k++)
        z_vec[k] = last_h_vec[k]*j*last_h_dist;

      double inPt[3];
      vtkMath::Add(newPt, z_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }
  }

  int l_div = dim[1];

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      pos[0] = i; pos[2] = j;

      pos[1] = 0;
      int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

      double firstPt[3];
      paraHexMesh->GetPoint(firstPtId, firstPt);

      pos[1] = dim[1]-1;
      int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

      double lastPt[3];
      paraHexMesh->GetPoint(lastPtId, lastPt);

      double l = vtkSVMathUtils::Distance(lastPt, firstPt);
      double l_dist = l/(l_div-1);

      double l_vec[3];
      vtkMath::Subtract(lastPt, firstPt, l_vec);
      vtkMath::Normalize(l_vec);

      for (int k=0; k<l_div; k++)
      {
        pos[1] = k;

        double y_vec[3];
        for (int l=0; l<3; l++)
          y_vec[l] = l_vec[l]*k*l_dist;

        double inPt[3];
        vtkMath::Add(firstPt, y_vec, inPt);

        int newPtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// GetInteriorPointMaps
// ----------------------
int vtkSVGroupsSegmenter::GetInteriorPointMaps(vtkPolyData *pdWithAllInterior,
                                               vtkPolyData *pdWithCleanInterior,
                                               vtkPolyData *pdWithoutInterior,
                                               std::vector<int> &ptMap,
                                               std::vector<std::vector<int> > &invPtMap)
{
  vtkNew(vtkPointLocator, locator);
  locator->SetDataSet(pdWithoutInterior);
  locator->BuildLocator();

  vtkNew(vtkPointLocator, locator2);
  locator2->SetDataSet(pdWithCleanInterior);
  locator2->BuildLocator();

  int numPoints = pdWithAllInterior->GetNumberOfPoints();
  ptMap.resize(numPoints);
  std::fill(ptMap.begin(), ptMap.end(), -1);

  int numCleanPoints = pdWithCleanInterior->GetNumberOfPoints();
  invPtMap.resize(numCleanPoints);

  for (int i=0; i<numPoints; i++)
  {
    double pt0[3];
    pdWithAllInterior->GetPoint(i, pt0);

    int ptId = locator->FindClosestPoint(pt0);

    double pt1[3];
    pdWithoutInterior->GetPoint(ptId, pt1);

    double dist = vtkSVMathUtils::Distance(pt0, pt1);

    if (dist > 1.0e-6)
    {
      int cleanPtId = locator2->FindClosestPoint(pt0);
      ptMap[i] = cleanPtId;
      invPtMap[cleanPtId].push_back(i);
    }
  }

  return SV_OK;
}

// ----------------------
// GetVolumePointMap
// ----------------------
int vtkSVGroupsSegmenter::GetVolumePointMap(vtkUnstructuredGrid *ugAll,
                                            vtkUnstructuredGrid *ugClean,
                                            std::vector<int> &ptMap)
{
  vtkNew(vtkPointLocator, locator);
  locator->SetDataSet(ugClean);
  locator->BuildLocator();

  int numAllPoints = ugAll->GetNumberOfPoints();
  ptMap.resize(numAllPoints);

  for (int i=0; i<numAllPoints; i++)
  {
    double pt[3];
    ugAll->GetPoint(i, pt);

    int ptId = locator->FindClosestPoint(pt);

    ptMap[i] = ptId;
  }

  return SV_OK;
}

// ----------------------
// MapInteriorBoundary
// ----------------------
int vtkSVGroupsSegmenter::MapInteriorBoundary(vtkStructuredGrid *paraHexVolume,
                                              vtkPolyData *mappedSurface,
                                              const std::vector<int> ptMap)
{
  // Now lets try volume
  vtkDataArray *ptIds = mappedSurface->GetPointData()->GetArray("TmpInternalIds2");
  vtkDataArray *mappedIds = mappedSurface->GetPointData()->GetArray("TmpInternalIds");

  int dim[3];
  paraHexVolume->GetDimensions(dim);

  double first_para_xpt[3], first_para_ypt[3], first_para_zpt[3];
  double last_para_xpt[3], last_para_ypt[3], last_para_zpt[3];
  double first_real_xpt[3], first_real_ypt[3], first_real_zpt[3];
  double last_real_xpt[3], last_real_ypt[3], last_real_zpt[3];

  // Check all six boundaries first!
  int pos[3];

  for (int i=0; i<dim[0]; i++)
  {
    for (int j=0; j<dim[1]; j++)
    {
      for (int k=0; k<dim[2]; k++)
      {
        pos[0] = i; pos[1] = j; pos[2] = k;
        if (i == 0 || i == dim[0]-1)
        {
          int ptId = vtkStructuredData::ComputePointId(dim, pos);
          double para_pt[3];
          paraHexVolume->GetPoint(ptId, para_pt);
          int realId = mappedIds->LookupValue(ptId);

          if (ptMap[ptIds->GetTuple1(realId)] != -1)
          {
            pos[1] = 0;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, first_para_ypt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, first_real_ypt);

            pos[1] = dim[1]-1;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, last_para_ypt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, last_real_ypt);

            double para_y_dist = vtkSVMathUtils::Distance(para_pt, first_para_ypt)/vtkSVMathUtils::Distance(last_para_ypt, first_para_ypt);
            double yvec[3];
            vtkMath::Subtract(last_real_ypt, first_real_ypt, yvec);
            vtkMath::Normalize(yvec);
            vtkMath::MultiplyScalar(yvec, para_y_dist);

            double newPt0[3];
            vtkMath::Add(first_real_ypt, yvec, newPt0);

            pos[1] = j; pos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, first_para_zpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, first_real_zpt);

            pos[2] = dim[2]-1;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, last_para_zpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, last_real_zpt);

            double para_z_dist = vtkSVMathUtils::Distance(para_pt, first_para_zpt)/vtkSVMathUtils::Distance(last_para_zpt, first_para_zpt);
            double zvec[3];
            vtkMath::Subtract(last_real_zpt, first_real_zpt, zvec);
            vtkMath::Normalize(zvec);
            vtkMath::MultiplyScalar(zvec, para_z_dist);

            double newPt1[3];
            vtkMath::Add(first_real_zpt, zvec, newPt1);

            double newPt[3];
            vtkMath::Add(newPt0, newPt1, newPt);
            vtkMath::MultiplyScalar(newPt, 1./2);

            pos[0] = i; pos[1] = j; pos[2] = k;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoints()->SetPoint(realId, newPt);
          }
        }
        if (j == 0 || j == dim[1]-1)
        {
          int ptId = vtkStructuredData::ComputePointId(dim, pos);
          double para_pt[3];
          paraHexVolume->GetPoint(ptId, para_pt);
          int realId = mappedIds->LookupValue(ptId);

          if (ptMap[ptIds->GetTuple1(realId)] != -1)
          {
            pos[0] = 0;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, first_para_xpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, first_real_xpt);

            pos[0] = dim[0]-1;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, last_para_xpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, last_real_xpt);

            double para_x_dist = vtkSVMathUtils::Distance(para_pt, first_para_xpt)/vtkSVMathUtils::Distance(last_para_xpt, first_para_xpt);
            double xvec[3];
            vtkMath::Subtract(last_real_xpt, first_real_xpt, xvec);
            vtkMath::Normalize(xvec);
            vtkMath::MultiplyScalar(xvec, para_x_dist);

            double newPt0[3];
            vtkMath::Add(first_real_xpt, xvec, newPt0);

            pos[0] = i; pos[2] = 0;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, first_para_zpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, first_real_zpt);

            pos[2] = dim[2]-1;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, last_para_zpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, last_real_zpt);

            double para_z_dist = vtkSVMathUtils::Distance(para_pt, first_para_zpt)/vtkSVMathUtils::Distance(last_para_zpt, first_para_zpt);
            double zvec[3];
            vtkMath::Subtract(last_real_zpt, first_real_zpt, zvec);
            vtkMath::Normalize(zvec);
            vtkMath::MultiplyScalar(zvec, para_z_dist);

            double newPt1[3];
            vtkMath::Add(first_real_zpt, zvec, newPt1);

            double newPt[3];
            vtkMath::Add(newPt0, newPt1, newPt);
            vtkMath::MultiplyScalar(newPt, 1./2);

            pos[0] = i; pos[1] = j; pos[2] = k;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoints()->SetPoint(realId, newPt);
          }
        }
        if (k == 0 || k == dim[2]-1)
        {
          int ptId = vtkStructuredData::ComputePointId(dim, pos);
          double para_pt[3];
          paraHexVolume->GetPoint(ptId, para_pt);
          int realId = mappedIds->LookupValue(ptId);

          if (ptMap[ptIds->GetTuple1(realId)] != -1)
          {
            pos[0] = 0;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, first_para_xpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, first_real_xpt);

            pos[0] = dim[0]-1;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, last_para_xpt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, last_real_xpt);

            double para_x_dist = vtkSVMathUtils::Distance(para_pt, first_para_xpt)/vtkSVMathUtils::Distance(last_para_xpt, first_para_xpt);
            double xvec[3];
            vtkMath::Subtract(last_real_xpt, first_real_xpt, xvec);
            vtkMath::Normalize(xvec);
            vtkMath::MultiplyScalar(xvec, para_x_dist);

            double newPt0[3];
            vtkMath::Add(first_real_xpt, xvec, newPt0);

            pos[0] = i; pos[1] = 0;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, first_para_ypt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, first_real_ypt);

            pos[1] = dim[1]-1;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            paraHexVolume->GetPoint(ptId, last_para_ypt);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoint(realId, last_real_ypt);

            double para_y_dist = vtkSVMathUtils::Distance(para_pt, first_para_ypt)/vtkSVMathUtils::Distance(last_para_ypt, first_para_ypt);
            double yvec[3];
            vtkMath::Subtract(last_real_ypt, first_real_ypt, yvec);
            vtkMath::Normalize(yvec);
            vtkMath::MultiplyScalar(yvec, para_y_dist);

            double newPt1[3];
            vtkMath::Add(first_real_ypt, yvec, newPt1);

            double newPt[3];
            vtkMath::Add(newPt0, newPt1, newPt);
            vtkMath::MultiplyScalar(newPt, 1./2);

            pos[0] = i; pos[1] = j; pos[2] = k;
            ptId = vtkStructuredData::ComputePointId(dim, pos);
            realId = mappedIds->LookupValue(ptId);

            mappedSurface->GetPoints()->SetPoint(realId, newPt);
          }
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// FixInteriorBoundary
// ----------------------
int vtkSVGroupsSegmenter::FixInteriorBoundary(vtkPolyData *mappedSurface,
                                              const std::vector<std::vector<int> > invPtMap)
{
  vtkDataArray *ptIds = mappedSurface->GetPointData()->GetArray("TmpInternalIds2");

  for (int i=0; i<invPtMap.size(); i++)
  {
    double avgPt[3]; avgPt[0] = 0.0; avgPt[1] = 0.0; avgPt[2] = 0.0;
    int numPoints = 0;
    std::vector<int> realIds;
    for (int j=0; j<invPtMap[i].size(); j++)
    {
      numPoints++;
      realIds.push_back(ptIds->LookupValue(invPtMap[i][j]));
      double pt[3];
      mappedSurface->GetPoint(realIds[j], pt);

      for (int k=0; k<3; k++)
        avgPt[k] += pt[k];
    }

    if (numPoints > 0)
    {
      vtkMath::MultiplyScalar(avgPt, 1./numPoints);

      for (int j=0; j<invPtMap[i].size(); j++)
        mappedSurface->GetPoints()->SetPoint(realIds[j], avgPt);
    }
  }

  return SV_OK;
}

// ----------------------
// FixVolume
// ----------------------
int vtkSVGroupsSegmenter::FixVolume(vtkUnstructuredGrid *mappedVolume,
                                    vtkUnstructuredGrid *cleanVolume,
                                    const std::vector<int> ptMap)
{
  int numPoints = mappedVolume->GetNumberOfPoints();
  for (int i=0; i<numPoints; i++)
  {
    int cleanPtId = ptMap[i];

    double pt[3];
    cleanVolume->GetPoint(cleanPtId, pt);

    mappedVolume->GetPoints()->SetPoint(i, pt);
  }

  return SV_OK;
}

// ----------------------
// MapVolume
// ----------------------
int vtkSVGroupsSegmenter::MapVolume(vtkStructuredGrid *paraHexVolume,
                                    vtkPolyData *mappedSurface,
                                    vtkStructuredGrid *mappedVolume)
{
  // Now lets try volume
  vtkDataArray *volumeIds = paraHexVolume->GetPointData()->GetArray("TmpInternalIds");
  vtkDataArray *mappedIds = mappedSurface->GetPointData()->GetArray("TmpInternalIds");

  int dim[3];
  paraHexVolume->GetDimensions(dim);
  mappedVolume->SetDimensions(dim);
  vtkNew(vtkPoints, realHexMeshPoints);
  mappedVolume->SetPoints(realHexMeshPoints);
  mappedVolume->GetPoints()->SetNumberOfPoints(dim[0]*dim[1]*dim[2]);

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
        pos[0] = i; pos[1] = j; pos[2] = k;
        int ptId = vtkStructuredData::ComputePointId(dim, pos);

        pos[0] = 0; pos[1] = j; pos[2] = k;
        int x0PtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexVolume->GetPoint(x0PtId, first_para_xpt);
        int transId = volumeIds->GetTuple1(x0PtId);
        int realXId0 = mappedIds->LookupValue(transId);
        mappedSurface->GetPoint(realXId0, first_real_xpt);

        pos[0] = dim[0]-1; pos[1] = j; pos[2] = k;
        int x1PtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexVolume->GetPoint(x1PtId, last_para_xpt);
        transId = volumeIds->GetTuple1(x1PtId);
        int realXId1 = mappedIds->LookupValue(transId);
        mappedSurface->GetPoint(realXId1, last_real_xpt);

        pos[0] = i; pos[1] = 0; pos[2] = k;
        int y0PtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexVolume->GetPoint(y0PtId, first_para_ypt);
        transId = volumeIds->GetTuple1(y0PtId);
        int realYId0 = mappedIds->LookupValue(transId);
        mappedSurface->GetPoint(realYId0, first_real_ypt);

        pos[0] = i; pos[1] = dim[1]-1; pos[2] = k;
        int y1PtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexVolume->GetPoint(y1PtId, last_para_ypt);
        transId = volumeIds->GetTuple1(y1PtId);
        int realYId1 = mappedIds->LookupValue(transId);
        mappedSurface->GetPoint(realYId1, last_real_ypt);

        pos[0] = i; pos[1] = j; pos[2] = 0;
        int z0PtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexVolume->GetPoint(z0PtId, first_para_zpt);
        transId = volumeIds->GetTuple1(z0PtId);
        int realZId0 = mappedIds->LookupValue(transId);
        mappedSurface->GetPoint(realZId0, first_real_zpt);

        pos[0] = i; pos[1] = j; pos[2] = dim[2]-1;
        int z1PtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexVolume->GetPoint(z1PtId, last_para_zpt);
        transId = volumeIds->GetTuple1(z1PtId);
        int realZId1 = mappedIds->LookupValue(transId);
        mappedSurface->GetPoint(realZId1, last_real_zpt);

        pos[0] = i; pos[1] = j; pos[2] = k;

        double para_pt[3];
        ptId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexVolume->GetPoint(ptId, para_pt);

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
        {
          for (int r=0; r<3; r++)
              real_pt[r] = real_xpt[r];
        }
        else if (j == 0 || j == dim[1]-1)
        {
          for (int r=0; r<3; r++)
              real_pt[r] = real_ypt[r];
        }
        else if (k == 0 || k == dim[2]-1)
        {
          for (int r=0; r<3; r++)
              real_pt[r] = real_zpt[r];
        }
        else
        {
          vtkMath::Add(real_xpt, real_zpt, real_pt);
          vtkMath::MultiplyScalar(real_pt, 1./2);
        }

        mappedVolume->GetPoints()->SetPoint(ptId, real_pt);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// SmoothStructuredGrid
// ----------------------
int vtkSVGroupsSegmenter::SmoothStructuredGrid(vtkStructuredGrid *hexMesh, const int iters)
{
  int numCells = hexMesh->GetNumberOfCells();

  int dim[3];
  hexMesh->GetDimensions(dim);

  for (int iter=0; iter<iters; iter++)
  {
    for (int i=0; i<dim[0]; i++)
    {
      for (int j=0; j<dim[1]; j++)
      {
        for (int k=0; k<dim[2]; k++)
        {
          if (i == 0 || i == dim[0] - 1 ||
              j == 0 || j == dim[1] - 1 ||
              k == 0 || k == dim[2] - 1)
            continue;

          double center[3]; center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;

          int numNeigh = 6;
          int neighborPos[6][3] = {{i-1, j, k},
                                   {i+1, j, k},
                                   {i, j-1, k},
                                   {i, j+1, k},
                                   {i, j, k-1},
                                   {i, j, k+1}};

          for (int l=0; l<numNeigh; l++)
          {
            int neighborPtId = vtkStructuredData::ComputePointId(dim, neighborPos[l]);
            double neighborPt[3];
            hexMesh->GetPoint(neighborPtId, neighborPt);

            for (int m=0; m<3; m++)
              center[m] += neighborPt[m];
          }

          int pos[3]; pos[0] = i; pos[1] = j; pos[2] = k;
          int ptId = vtkStructuredData::ComputePointId(dim, pos);

          double pt[3];
          hexMesh->GetPoint(ptId, pt);

          for (int l=0; l<3; l++)
            pt[l] += (center[l]/numNeigh - pt[l]) * 0.02;

          hexMesh->GetPoints()->SetPoint(ptId, pt);
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// SmoothUnstructuredGrid
// ----------------------
int vtkSVGroupsSegmenter::SmoothUnstructuredGrid(vtkUnstructuredGrid *hexMesh, const int iters)
{
  hexMesh->BuildLinks();
  int numCells = hexMesh->GetNumberOfCells();
  int numPoints = hexMesh->GetNumberOfPoints();

  for (int i=0; i<numCells; i++)
  {
    if (hexMesh->GetCellType(i) != VTK_HEXAHEDRON)
    {
      fprintf(stdout,"All cells must be hexes\n");
      return SV_ERROR;
    }
  }

  std::vector<std::vector<int> > ptEdgeNeighbors(numPoints);

  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, ptCellIds);
    hexMesh->GetPointCells(i, ptCellIds);

    int interiorPoint = 1;
    vtkNew(vtkIdList, ptNeighbors);
    for (int j=0; j<ptCellIds->GetNumberOfIds(); j++)
    {
      vtkCell *cell = hexMesh->GetCell(ptCellIds->GetId(j));

      int numFaces = cell->GetNumberOfFaces();
      for (int k=0; k<numFaces; k++)
      {
        vtkCell *face = cell->GetFace(k);

        int checkable = 0;
        for (int l=0; l<4; l++)
        {
          if (face->PointIds->GetId(l) == i)
            checkable = 1;
        }

        if (checkable)
        {
          vtkNew(vtkIdList, neighCellIds);
          hexMesh->GetCellNeighbors(ptCellIds->GetId(j), face->PointIds, neighCellIds);
          if (neighCellIds->GetNumberOfIds() == 0)
            interiorPoint = 0;
        }
      }

      for (int k=0; k<cell->GetNumberOfEdges(); k++)
      {
        vtkIdList *edge = cell->GetEdge(k)->GetPointIds();
        int isPtId = edge->IsId(i);
        if (isPtId != -1)
        {
          if (ptNeighbors->IsId(edge->GetId((isPtId+1)%2)) == -1)
           ptNeighbors->InsertNextId(edge->GetId((isPtId+1)%2));
        }
      }
    }

    if (interiorPoint)
    {
      if (ptNeighbors->GetNumberOfIds() > 0)
      {
        for (int j=0; j<ptNeighbors->GetNumberOfIds(); j++)
          ptEdgeNeighbors[i].push_back(ptNeighbors->GetId(j));
      }
    }
  }

  for (int iter=0; iter<iters; iter++)
  {
    for (int i=0; i<numPoints; i++)
    {
      // If 6 neighbors, that means this is interior son
      int numPtNeighbors = ptEdgeNeighbors[i].size();
      if (numPtNeighbors > 0)
      {
        double center[3]; center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;

        for (int j=0; j<numPtNeighbors; j++)
        {
          int neighborPtId = ptEdgeNeighbors[i][j];
          double neighborPt[3];
          hexMesh->GetPoint(neighborPtId, neighborPt);

          for (int k=0; k<3; k++)
            center[k] += neighborPt[k];
        }

        double pt[3];
        hexMesh->GetPoint(i, pt);

        for (int j=0; j<3; j++)
          pt[j] += (center[j]/numPtNeighbors - pt[j]) * 0.02;

        hexMesh->GetPoints()->SetPoint(i, pt);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// RotateGroupToGlobalAxis
// ----------------------
int vtkSVGroupsSegmenter::RotateGroupToGlobalAxis(vtkPolyData *pd,
                                                  const int thresholdId,
                                                  std::string arrayName,
                                                  vtkPolyData *rotPd,
                                                  vtkMatrix4x4 *rotMatrix0,
                                                  vtkMatrix4x4 *rotMatrix1)
{
  vtkNew(vtkPolyData, thresholdPd);
  vtkSVGeneralUtils::ThresholdPd(pd, thresholdId, thresholdId, 1, arrayName, thresholdPd);
  thresholdPd->BuildLinks();

  vtkIdType f3npts, *f3PtIds;
  thresholdPd->GetCellPoints(3, f3npts, f3PtIds);

  double pts[3][3];
  for (int i=0; i<3; i++)
    thresholdPd->GetPoint(f3PtIds[i], pts[i]);

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

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(pd);
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
  cleaner->Update();
  rotPd->DeepCopy(cleaner->GetOutput());

  vtkSVGeneralUtils::ApplyRotationMatrix(rotPd, rotMatrix0);
  vtkSVGeneralUtils::ApplyRotationMatrix(rotPd, rotMatrix1);

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
    vtkSVGeneralUtils::GetPointCellsValues(ps, arrayName, i, pointCellValues);
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

// ----------------------
// MatchEndPatches
// ----------------------
int vtkSVGroupsSegmenter::MatchEndPatches(vtkPolyData *pd)
{
  std::vector<Region> branchRegions;
  if (vtkSVGroupsSegmenter::GetRegions(pd, "PatchIds", branchRegions) != SV_OK)
  {
    fprintf(stderr,"Could not get regions on branch\n");
    return SV_ERROR;
  }

  int numPoints = pd->GetNumberOfPoints();
  std::vector<int> tempCornerPoints;  // In ccw order
  std::vector<int> pointUsed(numPoints, 0);
  for (int j=0; j<branchRegions.size(); j++)
  {
    if (branchRegions[j].CornerPoints.size() != 4)
    {
      fprintf(stderr,"Number of corners on region %d is %d, needs to be 4\n", j, branchRegions[j].CornerPoints.size());
      std::string badName = "/Users/adamupdegrove/Desktop/tmp/BADCORNERS.vtp";
      vtkSVIOUtils::WriteVTPFile(badName, pd);
      return SV_ERROR;
    }
    for (int k=0; k<branchRegions[j].CornerPoints.size(); k++)
    {
      int cornerPtId = branchRegions[j].CornerPoints[k];
      vtkNew(vtkIdList, pointCellsValues);
      vtkSVGeneralUtils::GetPointCellsValues(pd, "PatchIds", cornerPtId, pointCellsValues);

      if (pointCellsValues->GetNumberOfIds() == 2 && pointUsed[cornerPtId] == 0)
      {
        pointUsed[cornerPtId] = 1;
        tempCornerPoints.push_back(cornerPtId);
      }
    }
  }

  std::vector<int> isCornerPoint(numPoints, 0);
  for (int j=0; j<tempCornerPoints.size(); j++)
    isCornerPoint[tempCornerPoints[j]] = 1;

  int startCornerPt = tempCornerPoints[0];

  std::vector<int> interiorCornerPoints;
  std::vector<std::vector<int> > ccwInteriorEdges;

  int count=1;
  std::vector<int> tempNodes;
  tempNodes.push_back(startCornerPt);
  interiorCornerPoints.push_back(startCornerPt);

  for (int j=0; j<count; j++)
  {
    vtkNew(vtkIdList, pointCells);
    pd->GetPointCells(tempNodes[j], pointCells);
    for (int k=0; k<pointCells->GetNumberOfIds(); k++)
    {
      int cellId = pointCells->GetId(k);
      int pointCCWId = vtkSVGroupsSegmenter::GetCCWPoint(pd, tempNodes[j], cellId);

      // Check if open edge
      vtkNew(vtkIdList, edgeCells);
      pd->GetCellEdgeNeighbors(cellId, tempNodes[j], pointCCWId, edgeCells);

      if (edgeCells->GetNumberOfIds() == 0 && !isCornerPoint[pointCCWId])
      {
        tempNodes.push_back(pointCCWId);
        count++;
      }
      else if (edgeCells->GetNumberOfIds() == 0 && isCornerPoint[pointCCWId])
      {
        if (pointCCWId == startCornerPt)
        {
          tempNodes.push_back(pointCCWId);
          ccwInteriorEdges.push_back(tempNodes);

          tempNodes.clear();

          if (interiorCornerPoints.size() == tempCornerPoints.size())
          {
            count = -1;
            break;
          }
          else
          {
            for (int ii=0; ii<tempCornerPoints.size(); ii++)
            {
              bool tempCount = false;
              int tempIndex = tempCornerPoints[ii];

              for (int jj=0; jj<interiorCornerPoints.size(); jj++)
              {
                if (tempIndex == interiorCornerPoints[jj])
                  tempCount = true;
              }
              if (tempCount == false)
              {
                startCornerPt = tempIndex;
                break;
              }
            }

            interiorCornerPoints.push_back(startCornerPt);
            tempNodes.push_back(startCornerPt);
            count = 1;
            j = -1;
            break;
          }
        }
        else
        {
          tempNodes.push_back(pointCCWId);
          interiorCornerPoints.push_back(pointCCWId);
          ccwInteriorEdges.push_back(tempNodes);
          tempNodes.clear();
          tempNodes.push_back(pointCCWId);
          count = 1;
          j = -1;
          break;
        }
      }
    }
  }
  //fprintf(stdout,"LETS SEE: %d %d\n", tempCornerPoints.size(), interiorCornerPoints.size());
  std::vector<std::vector<int> > cwInteriorEdges(interiorCornerPoints.size());

  // Now get cw edges by just doing opposite
  for (int j=0; j<interiorCornerPoints.size()/4; j++)
  {
    for (int k=0; k<4; k++)
    {
      int cwLoc = 4*j + k;
      int ccwLoc = 4*j + (k+3)%4;
      int ccwPointCount=ccwInteriorEdges[ccwLoc].size()-1;
      for (int l=0; l<ccwInteriorEdges[ccwLoc].size(); l++)
      {
        cwInteriorEdges[cwLoc].push_back(ccwInteriorEdges[ccwLoc][ccwPointCount]);
        ccwPointCount--;
      }
    }
  }

  vtkDataArray *slicePointsArray = pd->GetPointData()->GetArray("SlicePoints");
  for (int j=0; j<interiorCornerPoints.size(); j++)
  {
    int ccwCount = 0;
    for (int k=0; k<ccwInteriorEdges[j].size(); k++)
    {
      if (slicePointsArray->GetTuple1(ccwInteriorEdges[j][k]) == 1)
        break;
      ccwCount++;
    }
    int cwCount = 0;
    for (int k=0; k<cwInteriorEdges[j].size(); k++)
    {
      if (slicePointsArray->GetTuple1(cwInteriorEdges[j][k]) == 1)
        break;
      cwCount++;
    }
    fprintf(stdout,"POINT IS:                %d\n", interiorCornerPoints[j]);
    fprintf(stdout,"COUNTER CLOCKWISE COUNT: %d\n", ccwCount);
    fprintf(stdout,"FULL COUNT:              %d\n", ccwInteriorEdges[j].size());
    fprintf(stdout,"CLOCKWISE COUNT:         %d\n", cwCount);
    fprintf(stdout,"FULL COUNT:              %d\n", cwInteriorEdges[j].size());

    if (ccwCount == ccwInteriorEdges[j].size() && cwCount == cwInteriorEdges[j].size())
    {
      fprintf(stderr,"THIS IS A BIG PROBLEM\n");
      return SV_ERROR;
    }

    if (ccwCount == cwCount && ccwCount != 0)
    {
      fprintf(stderr,"THIS IS ALSO A BIG PROBLEM\n");
      return SV_ERROR;
    }

    int stopCount;
    int slicePoint;
    int startPoint;
    int secondPoint;
    int opposPoint;
    if (ccwCount < cwCount )
    {
      stopCount = ccwCount;
      slicePoint = ccwInteriorEdges[j][stopCount];
      startPoint = ccwInteriorEdges[j][0];
      secondPoint = ccwInteriorEdges[j][1];
      opposPoint = cwInteriorEdges[j][1];
    }
    else
    {
      stopCount = cwCount;
      slicePoint = cwInteriorEdges[j][stopCount];
      startPoint = cwInteriorEdges[j][0];
      secondPoint = cwInteriorEdges[j][1];
      opposPoint = ccwInteriorEdges[j][1];
    }

    vtkNew(vtkIdList, tmpCell);
    pd->GetCellEdgeNeighbors(-1, startPoint, secondPoint, tmpCell);
    int edgeCell = tmpCell->GetId(0);
    int matchCellValue = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(edgeCell);

    pd->GetCellEdgeNeighbors(-1, startPoint, opposPoint, tmpCell);
    int opposCell = tmpCell->GetId(0);
    int newCellValue = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(opposCell);

    int patchPoint = -1;
    if (stopCount != 0)
    {
      count = 1;
      tempNodes.clear();
      tempNodes.push_back(startPoint);

      for (int k=0; k<count; k++)
      {
        vtkNew(vtkIdList, pointCells);
        pd->GetPointCells(tempNodes[k], pointCells);
        for (int l=0; l<pointCells->GetNumberOfIds(); l++)
        {
          int cellId = pointCells->GetId(l);
          int nextPoint;
          if (ccwCount < cwCount)
            nextPoint = vtkSVGroupsSegmenter::GetCWPoint(pd, tempNodes[k], cellId);
          else
            nextPoint = vtkSVGroupsSegmenter::GetCCWPoint(pd, tempNodes[k], cellId);
          int isGoodEdge = vtkSVGroupsSegmenter::CheckCellValuesEdge(pd, "PatchIds", cellId, tempNodes[k], nextPoint);

          int cellValue = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(cellId);

          if (isGoodEdge && cellValue == matchCellValue)
          {
            if (count == stopCount)
            {
              count = -1;
              patchPoint = nextPoint;
              break;
            }
            tempNodes.push_back(nextPoint);
            count++;
          }
        }
      }
    }
    else
      patchPoint = startPoint;

    if (patchPoint == -1)
    {
      fprintf(stdout,"DIDNT COMPUTE PATCH POINT\n");
      return SV_ERROR;
    }

    fprintf(stdout,"PATCH POINT:              %d\n", patchPoint);
    fprintf(stdout,"SLICE POINT:              %d\n", slicePoint);
    if (stopCount != 0)
    {
      vtkNew(vtkSVFindGeodesicPath, finder);
      finder->SetInputData(pd);
      finder->SetStartPtId(patchPoint);
      finder->SetEndPtId(slicePoint);
      finder->SetDijkstraArrayName("DijkstraDistance");
      finder->SetRepelCloseBoundaryPoints(1);
      finder->Update();

      vtkNew(vtkIdList, tmpIds);
      tmpIds = finder->GetPathIds();
      int numToAdd = tmpIds->GetNumberOfIds();
      fprintf(stdout,"NEW POINTS:              ");
      for (int l=0; l<numToAdd; l++)
        fprintf(stdout,"%d ", tmpIds->GetId(l));
      fprintf(stdout,"\n");

      count = 1;
      std::vector<int> tempCells;
      tempCells.push_back(edgeCell);

      for (int k=0; k<count; k++)
      {
        int tmpCellId = tempCells[k];
        pd->GetCellData()->GetArray("PatchIds")->SetTuple1(tmpCellId, newCellValue);
        vtkIdType npts, *pts;
        pd->GetCellPoints(tmpCellId, npts, pts);
        for (int l=0; l<npts; l++)
        {
          int ptId0 = pts[l];
          int ptId1 = pts[(l+1)%npts];

          int freeEdge =  0;
          int patchEdge = 0;
          int newEdge =   0;

          vtkNew(vtkIdList, cellEdgeNeighbors);
          pd->GetCellEdgeNeighbors(tmpCellId, ptId0, ptId1, cellEdgeNeighbors);

          if (cellEdgeNeighbors->GetNumberOfIds() == 0)
            freeEdge = 1;
          else
          {
            int testCellId = cellEdgeNeighbors->GetId(0);
            int cellValue = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(testCellId);
            if (cellValue == newCellValue)
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

  return SV_OK;
}

// ----------------------
// FixGroupsWithPolycube
// ----------------------
int vtkSVGroupsSegmenter::FixGroupsWithPolycube()
{
  // Then check everything
  // Extract surface, triangulate, and subdivide polycube
  vtkNew(vtkPolyData, polycubePd);
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(this->Polycube);
  surfacer->Update();

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(surfacer->GetOutput());
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
  cleaner->Update();

  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(cleaner->GetOutput());
  triangulator->Update();

  polycubePd->DeepCopy(triangulator->GetOutput());
  polycubePd->BuildLinks();

  fprintf(stdout,"GETTING SURFACE GROUPS\n");
  std::vector<Region> surfaceGroups;
  if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, surfaceGroups) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  fprintf(stdout,"GETTING POLYCUBE GROUPS\n");
  std::vector<Region> polycubeGroups;
  if (this->GetRegions(polycubePd, this->GroupIdsArrayName, polycubeGroups) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  int numSurfaceGroups  = surfaceGroups.size();
  int numPolycubeGroups = polycubeGroups.size();

  //if (numSurfaceGroups != numPolycubeGroups)
  //{
  //  fprintf(stdout,"TELL ME, NEED TO FIX ADDITIONAL GROUP SOMEWHERE\n");
  //  return SV_ERROR;
  //}

  for (int i=0; i<numSurfaceGroups; i++)
  {
    for (int j=0; j<numPolycubeGroups; j++)
    {
      if (surfaceGroups[i].IndexCluster == polycubeGroups[j].IndexCluster)
      {
        if (surfaceGroups[i].BoundaryEdges.size() != polycubeGroups[j].BoundaryEdges.size())
        {
          fprintf(stdout,"NUMRBS OF EDGES FROM SURFACE GROUP  %d is %d\n", surfaceGroups[i].IndexCluster, surfaceGroups[i].BoundaryEdges.size());
          fprintf(stdout,"NUMRBS OF EDGES FROM POLYCUBE GROUP %d is %d\n", polycubeGroups[j].IndexCluster, polycubeGroups[j].BoundaryEdges.size());
        }
      }
    }
  }

  return SV_OK;
}

// ----------------------
// FixGroups
// ----------------------
int vtkSVGroupsSegmenter::FixGroups(vtkPolyData *pd, std::string arrayName,
                                             std::vector<Region> allRegions)
{

  vtkSVGeneralUtils::GiveIds(pd, "TmpInternalIds");

  vtkNew(vtkPolyData, origPd);
  origPd->DeepCopy(pd);

  vtkNew(vtkIdList, halfIds);
  for (int i=0; i<allRegions.size(); i++)
  {
    int patchValue = allRegions[i].IndexCluster;
    fprintf(stdout,"CLUSTER %d\n", patchValue);
    int numEdges = allRegions[i].BoundaryEdges.size();
    fprintf(stdout,"NUM EDGES %d\n", numEdges);

    int edgeCount = 0;

    while( edgeCount < numEdges)
    {

      int ptId0 = allRegions[i].BoundaryEdges[edgeCount][0];

      int ptIdN = -1;

      std::vector<int> badEdges;
      std::vector<int> groupEdges;
      while (ptIdN != ptId0)
      {
        groupEdges.push_back(edgeCount);

        int edgeSize = allRegions[i].BoundaryEdges[edgeCount].size();

        int edgePtId0 = allRegions[i].BoundaryEdges[edgeCount][0];
        vtkNew(vtkIdList, ptId0List);
        vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, edgePtId0, ptId0List);

        ptIdN = allRegions[i].BoundaryEdges[edgeCount][edgeSize-1];
        vtkNew(vtkIdList, ptIdNList);
        vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, ptIdN, ptIdNList);

        fprintf(stdout,"PT ID 0: %d\n", edgePtId0);
        fprintf(stdout,"IDS 0: %d %d %d\n", ptId0List->GetId(0), ptId0List->GetId(1), ptId0List->GetId(2));
        fprintf(stdout,"PT ID N: %d\n", ptIdN);
        fprintf(stdout,"IDS N: %d %d %d\n", ptIdNList->GetId(0), ptIdNList->GetId(1), ptIdNList->GetId(2));
        fprintf(stdout,"\n");

        int numIds = ptId0List->GetNumberOfIds();
        ptId0List->IntersectWith(ptIdNList);
        if (ptId0List->GetNumberOfIds() != numIds)
          badEdges.push_back(edgeCount);

        edgeCount++;
      }

      fprintf(stdout,"NUMBER OF BAD EDGES: %d\n", badEdges.size());
      fprintf(stdout,"BAD EDGES: ");
      for (int j=0; j<badEdges.size(); j++)
        fprintf(stdout,"%d ", badEdges[j]);
      fprintf(stdout,"\n");
      if (badEdges.size() > 0)
      {
        vtkNew(vtkPolyData, branchPd);
        vtkSVGeneralUtils::ThresholdPd(origPd, patchValue, patchValue, 1, arrayName, branchPd);
        branchPd->BuildLinks();

        fprintf(stdout,"NEEED TO DO SOMETHING!!!!\n");

        for (int j=0; j<badEdges.size(); j++)
        {
          int badEdgeId = badEdges[j];
          int edgeSize  = allRegions[i].BoundaryEdges[badEdgeId].size();

          int startPtId = allRegions[i].BoundaryEdges[badEdgeId][0];
          int finalPtId = allRegions[i].BoundaryEdges[badEdgeId][edgeSize-1];

          int halfSize;
          if (edgeSize%2 == 0)
          {
            int testId = allRegions[i].BoundaryEdges[badEdgeId][edgeSize/2-1];
            fprintf(stdout,"TEST ID: %d\n", testId);
            if (halfIds->IsId(testId) == -1)
              halfSize = edgeSize/2;
            else
              halfSize = edgeSize/2-1;
          }
          else
          {
            halfSize = floor(edgeSize/2.);
          }
          int halfId   = allRegions[i].BoundaryEdges[badEdgeId][halfSize];
          halfIds->InsertUniqueId(halfId);

          for (int k=0; k<groupEdges.size(); k++)
          {
            if (groupEdges[k] != badEdgeId)
            {
              int groupEdgeSize = allRegions[i].BoundaryEdges[groupEdges[k]].size();

              //// TODO: Check half size isn't larger then gropuEdgeSize!!!
              //if (halfSize > groupEdgeSize)
              //{
              //  fprintf(stderr,"We have MAJOR PROBLEMO\n");
              //  return SV_ERROR;
              //}

              int stopId       = -1;
              int edgeCell     = -1;
              int newCellValue = -1;
              if (allRegions[i].BoundaryEdges[groupEdges[k]][0] == finalPtId)
              {
                stopId = allRegions[i].BoundaryEdges[groupEdges[k]][halfSize];

                vtkNew(vtkIdList, halfValues);
                vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, halfValues);

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

              if (allRegions[i].BoundaryEdges[groupEdges[k]][groupEdgeSize-1] == startPtId)
              {
                stopId = allRegions[i].BoundaryEdges[groupEdges[k]][groupEdgeSize-halfSize-1];

                vtkNew(vtkIdList, halfValues);
                vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, halfValues);

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

              if (stopId != -1)
              {
                if (edgeCell == -1)
                {
                  fprintf(stderr,"CANNOT HAVE UNDEFINED CELL TO START PAINTING\n");
                  return SV_ERROR;
                }
                if (newCellValue == -1)
                {
                  fprintf(stderr,"CANNOT HAVE UNDEFINED CELL VALUE TO PAINT WITH\n");
                  return SV_ERROR;
                }

                fprintf(stdout,"PLANNING PATH FROM %d to %d\n", halfId, stopId);
                fprintf(stdout,"STARTING USING EDGE CELL %d AND PAINTING WITH %d\n", edgeCell, newCellValue);

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
                fprintf(stdout,"NEW POINTS:              ");
                for (int l=0; l<numToAdd; l++)
                  fprintf(stdout,"%d ", tmpIds->GetId(l));
                fprintf(stdout,"\n");

                int count = 1;
                std::vector<int> tempCells;
                tempCells.push_back(edgeCell);

                for (int l=0; l<count; l++)
                {
                  int tmpCellId = tempCells[l];
                  branchPd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(tmpCellId, newCellValue);
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
        // Set vals on work pd
        for (int l=0; l<branchPd->GetNumberOfCells(); l++)
        {
          //Get real cell id
          int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(l);

          // Get val
          int cellVal = branchPd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(l);

          // Set val
          pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(realCellId, cellVal);
        }
      }
      fprintf(stdout,"\n");
    }
  }

  for (int i=0; i<halfIds->GetNumberOfIds(); i++)
  {
    int halfId = halfIds->GetId(i);

    vtkNew(vtkIdList, checkVals);
    vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, checkVals);

    vtkNew(vtkIdList, finalVals);
    vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, halfId, finalVals);
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

      std::vector<std::vector<int> > allNodes;

      int count=1;
      std::vector<int> tempNodes;
      tempNodes.push_back(halfId);

      for (int j=0; j<count; j++)
      {
        vtkNew(vtkIdList, badCells);
        pd->GetPointCells(tempNodes[j], badCells);

        for (int k=0; k<badCells->GetNumberOfIds(); k++)
        {
          int cellId = badCells->GetId(k);
          int pointCCWId = vtkSVGroupsSegmenter::GetCCWPoint(pd, tempNodes[j], cellId);

          vtkNew(vtkIdList, cellEdgeNeighbors);
          pd->GetCellEdgeNeighbors(cellId, tempNodes[j], pointCCWId, cellEdgeNeighbors);

          int edgeVal0 = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(cellId);
          int edgeVal1 = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(cellEdgeNeighbors->GetId(0));

          if (edgeVal0 == badVals->GetId(0) && edgeVal1 == badVals->GetId(1))
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
            int tmpVal = badVals->GetId(0);
            badVals->SetId(0, badVals->GetId(1));
            badVals->SetId(1, tmpVal);
            count=1;
            j=-1;
            break;
          }
        }
        if (tempNodes.size() == 1 && allNodes.size() == 0)
        {
          int tmpVal = badVals->GetId(0);
          badVals->SetId(0, badVals->GetId(1));
          badVals->SetId(1, tmpVal);
          count=1;
          j=-1;
        }
      }

      for (int j=0; j<allNodes.size(); j++)
      {
        fprintf(stdout,"\n");
        for (int k=1; k<allNodes[j].size(); k++)
        {
          int ptId = allNodes[j][k];

          vtkNew(vtkIdList, pointCells);
          pd->GetPointCells(ptId, pointCells);

          for (int l=0; l<pointCells->GetNumberOfIds(); l++)
          {
            int cellVal = pd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(pointCells->GetId(l));
            if (cellVal == badVals->GetId(j))
            {
              pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(pointCells->GetId(l), missingVals->GetId(j));
            }
          }
        }
      }
      fprintf(stdout,"\n");
    }
  }

  pd->GetCellData()->RemoveArray("TmpInternalIds");
  pd->GetPointData()->RemoveArray("TmpInternalIds");
  return SV_OK;
}

// ----------------------
// SplitUpBoundaries
// ----------------------
int vtkSVGroupsSegmenter::SplitUpBoundaries(vtkPolyData *pd, std::string arrayName,
                                             std::vector<Region> allRegions)
{
  vtkNew(vtkIntArray, newSlicePointsArray);
  newSlicePointsArray->SetNumberOfTuples(pd->GetNumberOfPoints());
  newSlicePointsArray->FillComponent(0, -1);
  newSlicePointsArray->SetName("SlicePoints");
  pd->GetPointData()->AddArray(newSlicePointsArray);
  for (int i=0; i<allRegions.size(); i++)
  {
    int numEdges = allRegions[i].BoundaryEdges.size();
    for (int j=0; j<numEdges; j++)
    {
      int edgeSize = allRegions[i].BoundaryEdges[j].size();

      vtkNew(vtkIdList, ptId0List);
      int ptId0 = allRegions[i].BoundaryEdges[j][0];
      vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, ptId0, ptId0List);

      vtkNew(vtkIdList, ptIdNList);
      int ptIdN = allRegions[i].BoundaryEdges[j][edgeSize-1];
      vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, ptIdN, ptIdNList);

      fprintf(stdout,"PT ID 0: %d\n", ptId0);
      fprintf(stdout,"IDS 0: %d %d %d\n", ptId0List->GetId(0), ptId0List->GetId(1), ptId0List->GetId(2));
      fprintf(stdout,"PT ID N: %d\n", ptIdN);
      fprintf(stdout,"IDS N: %d %d %d\n", ptIdNList->GetId(0), ptIdNList->GetId(1), ptIdNList->GetId(2));
      fprintf(stdout,"\n");

      ptId0List->IntersectWith(ptIdNList);

      if (ptId0List->GetNumberOfIds() >= 3)
      {
        // Traditional between sides of groups
        vtkSVGroupsSegmenter::SplitBoundary(pd, allRegions[i].BoundaryEdges[j], 3, allRegions[i].IndexCluster);
      }
      else if (ptId0List->GetNumberOfIds() == 2)
      {
        // Between center of groups, need to do special
        vtkSVGroupsSegmenter::SplitBoundary(pd, allRegions[i].BoundaryEdges[j], 2, allRegions[i].IndexCluster);
      }
      else
      {
        fprintf(stderr,"Not sure where this case should happen, not implemented\n");
        return SV_ERROR;
      }

    }
  }

  return SV_OK;
}

// ----------------------
// SplitBoundary
// ----------------------
int vtkSVGroupsSegmenter::SplitBoundary(vtkPolyData *pd,
                                         std::vector<int> boundary,
                                         int numDivs,
                                         int groupId)
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
      if (currDiv == numDivs)
        break;

      if (slicePoints->GetTuple1(ptId0) == -1)
        slicePoints->SetTuple1(ptId1, 1);

    }
  }

  return SV_OK;
}
