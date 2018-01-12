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

#include "vtkSVCenterlineBranchSplitter.h"
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
#include "vtkSVLoftNURBSVolume.h"
#include "vtkSVMapInterpolator.h"
#include "vtkSVMUPFESNURBSWriter.h"
#include "vtkSVNURBSCollection.h"
#include "vtkSVNURBSUtils.h"
#include "vtkSVPassDataArray.h"
#include "vtkSVPERIGEENURBSCollectionWriter.h"
#include "vtkSVPlanarMapper.h"
#include "vtkSVPointSetBoundaryMapper.h"

#include "vtkAppendPolyData.h"
#include "vtkAppendFilter.h"
#include "vtkExecutive.h"
#include "vtkCellArray.h"
#include "vtkCellLocator.h"
#include "vtkConnectivityFilter.h"
#include "vtkPolyLine.h"
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
#include "vtkLine.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkMath.h"
#include "vtkMergeCells.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkSmartPointer.h"
#include "vtkSphere.h"
#include "vtkSplineFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyDataNormals.h"
#include "vtkStructuredGridGeometryFilter.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"
#include "vtkThreshold.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVersion.h"
#include "vtkXMLPolyDataWriter.h"

#include "vtkvmtkCenterlineAttributesFilter.h"
#include "vtkvmtkMath.h"
#include "vtkvmtkMergeCenterlines.h"
#include "vtkvmtkPolyDataCenterlineGroupsClipper.h"

#include <algorithm>

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
  this->PolycubePd = vtkPolyData::New();
  this->Centerlines = NULL;

  this->PolycubeUg   = vtkUnstructuredGrid::New();
  this->FinalHexMesh = vtkUnstructuredGrid::New();

  this->CenterlineGroupIdsArrayName = NULL;
  this->CenterlineRadiusArrayName = NULL;
  this->GroupIdsArrayName = NULL;
  this->BlankingArrayName = NULL;
  this->CenterlineGroupIds = NULL;

  this->UseVmtkClipping = 0;
  this->ClipAllCenterlineGroupIds = 1;
  this->EnforceBoundaryDirections = 1;
  this->CutoffRadiusFactor = VTK_SV_LARGE_DOUBLE;
  this->ClipValue = 0.0;
  this->UseRadiusInformation = 1;

  this->PolycubeDivisions = 5;
  this->PolycubeUnitLength = 0.0;

  this->NormalsWeighting = 0.8;
  this->IsVasculature = 1;
  this->NumberOfCenterlineRemovePts = 3;
  this->BoundaryEnforceFactor = 1;

  this->UseAbsoluteMergeDistance = 0;
  this->RadiusMergeRatio = 0.5;
  this->MergeDistance = 0.1;
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
  if (this->PolycubePd != NULL)
  {
    this->PolycubePd->Delete();
    this->PolycubePd = NULL;
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

  if (this->PolycubeUg != NULL)
  {
    this->PolycubeUg->Delete();
    this->PolycubeUg = NULL;
  }

  if (this->FinalHexMesh != NULL)
  {
    this->FinalHexMesh->Delete();
    this->FinalHexMesh = NULL;
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
    output->DeepCopy(this->WorkPd);
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

  double polycubeSize;
  if (this->PolycubeUnitLength == 0.0)
  {
    this->GetApproximatePolycubeSize(polycubeSize);
    this->PolycubeUnitLength = polycubeSize/this->PolycubeDivisions;
  }
  else
  {
    polycubeSize = this->PolycubeUnitLength*this->PolycubeDivisions;
  }

  this->CenterlineGraph = new vtkSVCenterlineGraph(0, this->MergedCenterlines,
                                                this->GroupIdsArrayName);
  this->CenterlineGraph->SetCubeSize(polycubeSize);

  if (this->CenterlineGraph->BuildGraph() != SV_OK)
  {
    vtkErrorMacro("Unable to form graph of centerlines");
    return SV_ERROR;
  }

  this->CenterlineGraph->GetGraphPolyData(this->GraphPd);
  this->MergedCenterlines->DeepCopy(this->CenterlineGraph->Lines);


  this->CenterlineGraph->GetSurfacePolycube(polycubeSize, this->PolycubePd);

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(this->PolycubePd);
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
  cleaner->Update();

  this->PolycubePd->DeepCopy(cleaner->GetOutput());
  this->PolycubePd->BuildLinks();

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVGroupsSegmenter::RunFilter()
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

  int stopCellNumber = ceil(this->WorkPd->GetNumberOfCells()*0.0001);

  if  (this->UseVmtkClipping)
  {
    vtkNew(vtkSplineFilter, resampler);
    resampler->SetInputData(this->Centerlines);
    //resampler->SetInputData(this->MergedCenterlines);
    resampler->SetSubdivideToLength();
    resampler->SetLength(this->Centerlines->GetLength()/100.);
    resampler->Update();

    vtkNew(vtkvmtkPolyDataCenterlineGroupsClipper, branchClipper);
    branchClipper->SetInputData(this->WorkPd);
    branchClipper->SetCenterlines(resampler->GetOutput());
    branchClipper->SetGroupIdsArrayName(this->GroupIdsArrayName);
    branchClipper->SetCenterlineGroupIdsArrayName(this->CenterlineGroupIdsArrayName);
    branchClipper->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
    branchClipper->SetBlankingArrayName(this->BlankingArrayName);
    branchClipper->SetCutoffRadiusFactor(this->CutoffRadiusFactor);
    branchClipper->SetClipValue(this->ClipValue);
    //branchClipper->SetUseRadiusInformation(this->UseRadiusInformation);
    branchClipper->SetUseRadiusInformation(0);
    branchClipper->SetClipAllCenterlineGroupIds(this->ClipAllCenterlineGroupIds);
    branchClipper->Update();

    vtkNew(vtkSVPassDataArray, dataPasser);
    dataPasser->SetInputData(0, branchClipper->GetOutput());
    dataPasser->SetInputData(1, this->WorkPd);
    dataPasser->SetPassArrayName(this->GroupIdsArrayName);
    dataPasser->SetPassDataIsCellData(0);
    dataPasser->SetPassDataToCellData(1);
    dataPasser->Update();

    this->WorkPd->DeepCopy(dataPasser->GetOutput());
  }
  else
  {
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
    //CVT->SetUseRadiusInformation(0);
    //CVT->SetUsePointNormal(1);
    CVT->SetUsePointNormal(0);
    CVT->SetUseBifurcationInformation(0);
    //CVT->SetUseBifurcationInformation(1);
    CVT->SetMaximumNumberOfIterations(0);
    CVT->Update();

    this->WorkPd->DeepCopy(CVT->GetOutput());
  }

  int firstSmooth = 0;
  if (firstSmooth)
  {
    vtkNew(vtkSVEdgeWeightedSmoother, smoother);
    smoother->SetInputData(this->WorkPd);
    smoother->SetGenerators(this->MergedCenterlines);
    smoother->SetNumberOfRings(2);
    smoother->SetThreshold(stopCellNumber);
    smoother->SetUseCurvatureWeight(0);
    smoother->SetNoInitialization(1);
    smoother->SetPatchIdsArrayName(this->GroupIdsArrayName);
    smoother->SetCVTDataArrayName("Normals");
    smoother->Update();

    this->WorkPd->DeepCopy(smoother->GetOutput());
  }

  //double groupRange[2];
  //this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetRange(groupRange);
  //vtkNew(vtkThreshold, tmpGroupThresholder);
  //tmpGroupThresholder->SetInputData(this->WorkPd);
  //tmpGroupThresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName);
  //for (int i=groupRange[0]; i<groupRange[1]; i++)
  //{
  //  int groupIdVal = i;
  //  tmpGroupThresholder->ThresholdBetween(groupIdVal, groupIdVal);
  //  tmpGroupThresholder->Update();

  //  if (tmpGroupThresholder->GetOutput()->GetNumberOfPoints() > 0)
  //  {
  //    vtkNew(vtkDataSetSurfaceFilter, newSurfacer);
  //    newSurfacer->SetInputData(tmpGroupThresholder->GetOutput());
  //    newSurfacer->Update();

  //    std::string oneGroupFn = "/Users/adamupdegrove/Desktop/tmp/GROUP_" + std::to_string(groupIdVal) + ".vtp";
  //    vtkSVIOUtils::WriteVTPFile(oneGroupFn, newSurfacer->GetOutput());
  //  }
  //}

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

  //if (this->CheckGroups2() != SV_OK)
  //{
  //  vtkErrorMacro("Error in correcting groups");
  //  return SV_ERROR;
  //}

  //if (this->FixGroupsWithPolycube() != SV_OK)
  //{
  //  vtkErrorMacro("Error in correcting groups");
  //  return SV_ERROR;
  //}

  //vtkSVIOUtils::WriteVTPFile("/Users/adamupdegrove/Desktop/tmp/AFTERSMOOTH_1.vtp", this->WorkPd);

  //if (this->CorrectCellBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
  //{
  //  vtkErrorMacro("Could not correcto boundaries of surface");
  //  return SV_ERROR;
  //}
  //vtkSVIOUtils::WriteVTPFile("/Users/adamupdegrove/Desktop/tmp/AFTERCORRECTION_2.vtp", this->WorkPd);

  //if (this->SmoothBoundaries(this->WorkPd, this->GroupIdsArrayName) != SV_OK)
  //{
  //  vtkErrorMacro("Could not smootho boundaries of surface");
  //  return SV_ERROR;
  //}
  //vtkSVIOUtils::WriteVTPFile("/Users/adamupdegrove/Desktop/tmp/AFTERSMOOTH_2.vtp", this->WorkPd);

  //std::vector<Region> groupRegions;
  //if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, groupRegions) != SV_OK)
  //{
  //  vtkErrorMacro("Couldn't get group regions");
  //  return SV_ERROR;
  //}
  //if (this->CurveFitBoundaries(this->WorkPd, this->GroupIdsArrayName, groupRegions) != SV_OK)
  //{
  //  vtkErrorMacro("Could not curve fit boundaries of surface");
  //  return SV_ERROR;
  //}
  //vtkSVIOUtils::WriteVTPFile("/Users/adamupdegrove/Desktop/tmp/AFTERCURVEFIT.vtp", this->WorkPd);

  //if (this->MatchSurfaceToPolycube() != SV_OK)
  //{
  //  vtkErrorMacro("Couldn't fix stuff\n");
  //  return SV_ERROR;
  //}

  //if (this->CheckSlicePoints() != SV_OK)
  //{
  //  vtkErrorMacro("Error when checking slice points\n");
  //  return SV_ERROR;
  //}

  //// Get new normals
  //normaler->SetInputData(this->WorkPd);
  //normaler->ComputePointNormalsOff();
  //normaler->ComputeCellNormalsOn();
  //normaler->SplittingOff();
  //normaler->Update();
  //this->WorkPd->DeepCopy(normaler->GetOutput());
  //this->WorkPd->BuildLinks();

  //int numberOfCells = this->WorkPd->GetNumberOfCells();

  //// Add array for new cell normals on surface
  //vtkNew(vtkDoubleArray, preRotationNormals);
  //preRotationNormals->SetName("PreRotationNormals");
  //preRotationNormals->SetNumberOfComponents(3);
  //preRotationNormals->SetNumberOfTuples(numberOfCells);

  //vtkNew(vtkDoubleArray, centerlineBasedNormals);
  //centerlineBasedNormals->SetName("CenterlinesBasedCellNormals");
  //centerlineBasedNormals->SetNumberOfComponents(3);
  //centerlineBasedNormals->SetNumberOfTuples(numberOfCells);

  //vtkNew(vtkDoubleArray, centerlineLocalX);
  //centerlineLocalX->SetName("ClosestCenterlineX");
  //centerlineLocalX->SetNumberOfComponents(3);
  //centerlineLocalX->SetNumberOfTuples(numberOfCells);

  //vtkNew(vtkDoubleArray, centerlineLocalY);
  //centerlineLocalY->SetName("ClosestCenterlineY");
  //centerlineLocalY->SetNumberOfComponents(3);
  //centerlineLocalY->SetNumberOfTuples(numberOfCells);

  //vtkNew(vtkDoubleArray, centerlineLocalZ);
  //centerlineLocalZ->SetName("ClosestCenterlineZ");
  //centerlineLocalZ->SetNumberOfComponents(3);
  //centerlineLocalZ->SetNumberOfTuples(numberOfCells);

  //vtkNew(vtkIntArray, centerlineSubPtIds);
  //centerlineSubPtIds->SetName("ClosestCenterlineSubPtId");
  //centerlineSubPtIds->SetNumberOfComponents(1);
  //centerlineSubPtIds->SetNumberOfTuples(numberOfCells);

  //vtkNew(vtkDoubleArray, centerlinePCoords);
  //centerlinePCoords->SetName("CenterlinePCoord");
  //centerlinePCoords->SetNumberOfComponents(1);
  //centerlinePCoords->SetNumberOfTuples(numberOfCells);

  //vtkNew(vtkDoubleArray, angularPCoords);
  //angularPCoords->SetName("AngularPCoord");
  //angularPCoords->SetNumberOfComponents(1);
  //angularPCoords->SetNumberOfTuples(numberOfCells);

  //// Get all group ids
  //vtkNew(vtkIdList, groupIds);
  //for (int i=0; i<this->WorkPd->GetNumberOfCells(); i++)
  //{
  //  int groupVal = this->WorkPd->GetCellData()->GetArray(
  //      this->GroupIdsArrayName)->GetTuple1(i);
  //  groupIds->InsertUniqueId(groupVal);
  //}
  //vtkSortDataArray::Sort(groupIds);
  //int numGroups = groupIds->GetNumberOfIds();

  ////vtkIntArray *tmpLinePtArray = vtkIntArray::New();
  //vtkDoubleArray *tmpLinePtArray = vtkDoubleArray::New();
  ////tmpLinePtArray->SetNumberOfComponents(3);
  //tmpLinePtArray->SetNumberOfTuples(this->WorkPd->GetNumberOfCells());
  //tmpLinePtArray->SetName("PatchVals");
  //for (int j=0; j<1; j++)
  //  tmpLinePtArray->FillComponent(j, -1);
  //this->WorkPd->GetCellData()->AddArray(tmpLinePtArray);
  //tmpLinePtArray->Delete();

  //vtkSVGeneralUtils::GiveIds(this->WorkPd, "TmpInternalIds");
  //for (int i=0; i<numGroups; i++)
  //{
  //  int groupId = groupIds->GetId(i);
  //  vtkNew(vtkPolyData, branchPd);
  //  vtkSVGeneralUtils::ThresholdPd(this->WorkPd, groupId, groupId, 1,
  //    this->GroupIdsArrayName, branchPd);
  //  branchPd->BuildLinks();

  //  vtkNew(vtkPolyData, centerlineBranchPd);
  //  vtkSVGeneralUtils::ThresholdPd(this->MergedCenterlines, groupId, groupId, 1,
  //    this->GroupIdsArrayName, centerlineBranchPd);
  //  centerlineBranchPd->BuildLinks();

  //  vtkNew(vtkPolyData, polyBranchPd);
  //  vtkSVGeneralUtils::ThresholdPd(this->PolycubePd, groupId, groupId, 1,
  //    this->GroupIdsArrayName, polyBranchPd);
  //  polyBranchPd->BuildLinks();

  //  // for each group, compute the clipping array, clip, add group ids array and append.
  //  vtkNew(vtkSVPolyBallLine, groupTubes);
  //  groupTubes->SetInput(centerlineBranchPd);
  //  groupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  //  groupTubes->SetUseRadiusInformation(this->UseRadiusInformation);
  //  groupTubes->UsePointNormalOff();
  //  groupTubes->UseRadiusWeightingOff();
  //  groupTubes->UseLocalCoordinatesOn();
  //  groupTubes->SetLocalCoordinatesArrayName("Local");

  //  vtkNew(vtkSVPolyBallLine, noRadiusTubes);
  //  noRadiusTubes->SetInput(centerlineBranchPd);
  //  noRadiusTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  //  noRadiusTubes->SetUseRadiusInformation(0);
  //  noRadiusTubes->UseLocalCoordinatesOn();
  //  noRadiusTubes->SetLocalCoordinatesArrayName("Local");

  //  int branchNumberOfCells = branchPd->GetNumberOfCells();
  //  // Loop through points to evaluate function at each point
  //  fprintf(stdout,"Computing closest centerline points per cell of group %d...\n", groupId);

  //  vtkIdType nlinepts, *linepts;
  //  int centerlineId = this->MergedCenterlines->GetCellData()->GetArray(this->GroupIdsArrayName)->LookupValue(groupId);
  //  this->MergedCenterlines->GetCellPoints(centerlineId, nlinepts, linepts);
  //  int isTerminating = 0;
  //  vtkNew(vtkIdList, testNeighbors);
  //  this->MergedCenterlines->GetPointCells(linepts[nlinepts-1], testNeighbors);
  //  if (testNeighbors->GetNumberOfIds() == 1)
  //    isTerminating = 1;

  //  for (int j=0; j<branchNumberOfCells; j++)
  //  {
  //    // Get cell point coords
  //    double pts[3][3];
  //    vtkIdType npts, *ptids;
  //    branchPd->GetCellPoints(j, npts, ptids);
  //    for (int k=0; k<npts; k++)
  //      branchPd->GetPoint(ptids[k], pts[k]);

  //    // Get center
  //    double center[3];
  //    vtkTriangle::TriangleCenter(pts[0], pts[1], pts[2], center);

  //    // Evaluate function at point!
  //    groupTubes->EvaluateFunction(center);

  //    //Get real cell id
  //    int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);

  //    // Now get last local coords and use rotation matrix to set new normals
  //    double localX[3], localY[3], localZ[3];
  //    groupTubes->GetLastLocalCoordX(localX);
  //    groupTubes->GetLastLocalCoordY(localY);
  //    groupTubes->GetLastLocalCoordZ(localZ);

  //    centerlineLocalX->SetTuple(realCellId, localX);
  //    centerlineLocalY->SetTuple(realCellId, localY);
  //    centerlineLocalZ->SetTuple(realCellId, localZ);

  //    double cellNormal[3];
  //    this->WorkPd->GetCellData()->GetArray("Normals")->GetTuple(realCellId, cellNormal);

  //    //TODO: JUST TESTING SOMETHING OUT!!!
  //    double closestPt[3];
  //    groupTubes->GetLastPolyBallCenter(closestPt);
  //    int linePtId = groupTubes->GetLastPolyBallCellSubId();

  //    if (linePtId >= nlinepts - 1)
  //    {
  //      fprintf(stderr,"Last point of line selected, didn't think that was possible\n");
  //      return SV_ERROR;
  //    }

  //    noRadiusTubes->EvaluateFunction(center);
  //    int absLinePtId = noRadiusTubes->GetLastPolyBallCellSubId();
  //    double absPCoord = noRadiusTubes->GetLastPolyBallCellPCoord();
  //    double absCenterlinePCoord = (absLinePtId + absPCoord)/(nlinepts-1);

  //    centerlineSubPtIds->SetTuple1(realCellId, absLinePtId);
  //    centerlinePCoords->SetTuple1(realCellId, absCenterlinePCoord);

  //    double absRadius = noRadiusTubes->GetLastPolyBallCenterRadius();
  //    double absClosestPt[3];
  //    noRadiusTubes->GetLastPolyBallCenter(absClosestPt);

  //    double absLocalX[3];
  //    noRadiusTubes->GetLastLocalCoordX(absLocalX);
  //    vtkMath::Normalize(absLocalX);


  //    if (absLinePtId >= nlinepts - 1)
  //    {
  //      fprintf(stderr,"Last point of line selected, didn't think that was possible\n");
  //      return SV_ERROR;
  //    }

  //    double centerlinePt0[3], centerlinePt1[3];
  //    centerlineBranchPd->GetPoint(absLinePtId, centerlinePt0);
  //    centerlineBranchPd->GetPoint(absLinePtId+1, centerlinePt1);

  //    double absTangent[3];
  //    vtkMath::Subtract(centerlinePt1, centerlinePt0, absTangent);
  //    vtkMath::Normalize(absTangent);

  //    double projectedPoint[3];
  //    vtkPlane::ProjectPoint(center, absClosestPt, absTangent, projectedPoint);

  //    double positionVector[3];
  //    vtkMath::Subtract(projectedPoint, absClosestPt, positionVector);
  //    vtkMath::Normalize(positionVector);

  //    double normalPoint[3];
  //    vtkMath::Add(absClosestPt, absLocalX, normalPoint);

  //    double projectedNormalPoint[3];
  //    vtkPlane::ProjectPoint(normalPoint, absClosestPt, absTangent, projectedNormalPoint);

  //    double projectedNormal[3];
  //    vtkMath::Subtract(projectedNormalPoint, absClosestPt, projectedNormal);
  //    vtkMath::Normalize(projectedNormal);

  //    double absCross[3];
  //    vtkMath::Cross(positionVector, projectedNormal, absCross);

  //    double tangentDot = vtkMath::Dot(absTangent, absCross);

  //    double absAngle = vtkvmtkMath::AngleBetweenNormals(positionVector, projectedNormal);

  //    if (tangentDot < 0.0)
  //      {
  //      absAngle *= -1.0;
  //      }

  //    angularPCoords->SetTuple1(realCellId, absAngle);

  //    double cellLocVec[3];
  //    vtkMath::Subtract(center, closestPt, cellLocVec);
  //    vtkMath::Normalize(cellLocVec);

  //    double orig_alpha = this->NormalsWeighting;
  //    double alpha = this->NormalsWeighting;

  //    if (linePtId <= 1)
  //    {
  //      if (!this->IsVasculature && this->MergedCenterlines->GetNumberOfCells() == 1)
  //        alpha = 0.0;
  //      else
  //        alpha = 1.0;
  //    }
  //    else if (linePtId >= nlinepts-4 && !isTerminating)
  //      alpha = 0.0;
  //    if (this->IsVasculature && linePtId >= nlinepts-4 && isTerminating)
  //      alpha = 1.0;

  //    double cellClusterVec[3];
  //    for (int k=0; k<3; k++)
  //      cellClusterVec[k] = alpha*cellNormal[k] + (1-alpha)*cellLocVec[k];
  //    vtkMath::Normalize(cellClusterVec);

  //    preRotationNormals->SetTuple(realCellId, cellClusterVec);
  //  }

  //  // Now go through and transform to local coordinate system and set
  //  // the new vector to use for clustering
  //  for (int j=0; j<branchNumberOfCells; j++)
  //  {
  //    //Get real cell id
  //    int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);

  //    double locals[6][3];
  //    centerlineLocalX->GetTuple(realCellId, locals[0]);
  //    centerlineLocalY->GetTuple(realCellId, locals[1]);
  //    centerlineLocalZ->GetTuple(realCellId, locals[4]);
  //    for (int k=0; k<3; k++)
  //    {
  //      locals[2][k] = -1.0*locals[0][k];
  //      locals[3][k] = -1.0*locals[1][k];
  //      locals[5][k] = -1.0*locals[4][k];
  //    }

  //    // Compute the rotation from global coordinate system to centerlines
  //    // local coordinate system
  //    double rotMat[9];
  //    this->ComputeRotationMatrix(locals[0], locals[1], locals[4], rotMat);

  //    double cellClusterVec[3];
  //    preRotationNormals->GetTuple(realCellId, cellClusterVec);

  //    // Apply rotation matrix to the normal to get the new normal
  //    double newNormal[3];
  //    for (int k=0; k<3; k++)
  //    {
  //      newNormal[k] = rotMat[k*3]*cellClusterVec[0] +
  //                     rotMat[(k*3)+1]*cellClusterVec[1] +
  //                     rotMat[(k*3)+2]*cellClusterVec[2];
  //    }

  //    centerlineBasedNormals->SetTuple(realCellId, newNormal);
  //  }
  //}


  //if (this->EnforceBoundaryDirections && this->MergedCenterlines->GetNumberOfCells() > 1)
  //{
  //  // Now enforce boundaries if we need to!!!!
  //  for (int i=0; i<numGroups; i++)
  //  {
  //    int groupId = groupIds->GetId(i);
  //    fprintf(stdout,"ENFORCING BOUNDARY OF GROUP: %d\n", groupId);

  //    vtkNew(vtkPolyData, branchPd);
  //    vtkSVGeneralUtils::ThresholdPd(this->WorkPd, groupId, groupId, 1,
  //      this->GroupIdsArrayName, branchPd);
  //    branchPd->BuildLinks();

  //    vtkNew(vtkPolyData, centerlineBranchPd);
  //    vtkSVGeneralUtils::ThresholdPd(this->MergedCenterlines, groupId, groupId, 1,
  //      this->GroupIdsArrayName, centerlineBranchPd);
  //    centerlineBranchPd->BuildLinks();

  //    vtkNew(vtkPolyData, polyBranchPd);
  //    vtkSVGeneralUtils::ThresholdPd(this->PolycubePd, groupId, groupId, 1,
  //      this->GroupIdsArrayName, polyBranchPd);
  //    polyBranchPd->BuildLinks();

  //    int branchNumberOfCells = branchPd->GetNumberOfCells();

  //    // Loop through points to evaluate function at each point
  //    fprintf(stdout,"Computing boundary vectors of group %d...\n", groupId);

  //    vtkIdType nlinepts, *linepts;
  //    int centerlineId = this->MergedCenterlines->GetCellData()->GetArray(this->GroupIdsArrayName)->LookupValue(groupId);
  //    this->MergedCenterlines->GetCellPoints(centerlineId, nlinepts, linepts);
  //    int isTerminating = 0;
  //    vtkNew(vtkIdList, testNeighbors);
  //    this->MergedCenterlines->GetPointCells(linepts[nlinepts-1], testNeighbors);
  //    if (testNeighbors->GetNumberOfIds() == 1)
  //      isTerminating = 1;

  //    double maxPCoord = -1.0;
  //    double minPCoord = VTK_SV_LARGE_DOUBLE;
  //    for (int j=0; j<branchNumberOfCells; j++)
  //    {
  //      //Get real cell id
  //      int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);

  //      double absCenterlinePCoord = centerlinePCoords->GetTuple1(realCellId);

  //      if (absCenterlinePCoord > maxPCoord)
  //        maxPCoord = absCenterlinePCoord;
  //      if (absCenterlinePCoord < minPCoord)
  //        minPCoord = absCenterlinePCoord;
  //    }

  //    for (int j=0; j<branchNumberOfCells; j++)
  //    {
  //      //Get real cell id
  //      int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);

  //      double currPCoord = centerlinePCoords->GetTuple1(realCellId);
  //      double newPCoord = (currPCoord - minPCoord)/(maxPCoord - minPCoord);

  //      centerlinePCoords->SetTuple1(realCellId, newPCoord);
  //    }

  //    fprintf(stdout,"MAX: %.6f MIN: %.6f\n", maxPCoord, minPCoord);
  //    // Do boundary cell stuffs
  //    // Get open boundary edges
  //    std::vector<int> openCornerPoints;
  //    std::vector<std::vector<int> > openEdges;
  //    if (this->GetOpenBoundaryEdges(branchPd, openCornerPoints, openEdges) != SV_OK)
  //    {
  //      fprintf(stderr,"Error getting open boundary edges\n");
  //      return SV_ERROR;
  //    }

  //    std::vector<std::vector<int> > shiftedOpenEdges;
  //    if (this->ShiftEdgeList(branchPd, openEdges, shiftedOpenEdges) != SV_OK)
  //    {
  //      fprintf(stderr,"Error shifting edges\n");
  //      return SV_ERROR;
  //    }

  //    // TODO NEEDS TO BE CHANGED FOR SPECIAL TRI CASE
  //    std::vector<std::vector<std::vector<double> > > allAngleBounds;
  //    std::vector<std::vector<int> > allPatchValues;
  //    std::vector<std::vector<int> > growCellLists;
  //    std::vector<int> cellBool(branchNumberOfCells);

  //    for (int j=0; j<branchNumberOfCells; j++)
  //      cellBool[j] = 0;

  //    for (int j=0; j<shiftedOpenEdges.size(); j++)
  //    {
  //      std::vector<std::vector<int> > splitOpenEdges;
  //      this->SplitEdgeList(branchPd, shiftedOpenEdges[j], splitOpenEdges);

  //      std::vector<std::vector<double> > edgeAngleBounds;
  //      std::vector<int> edgePatchValues;
  //      std::vector<int> edgeCellList;
  //      for (int k=0; k<splitOpenEdges.size(); k++)
  //      {
  //        int edgeSize = splitOpenEdges[k].size();
  //        if (edgeSize < 3)
  //        {
  //          fprintf(stderr,"EDGE SIZE IS LESS THAN 3, IT IS %d!!!\n", edgeSize);
  //          return SV_ERROR;
  //        }

  //        int edgePtId0 = branchPd->GetPointData()->GetArray("TmpInternalIds")->
  //          GetTuple1(splitOpenEdges[k][0]);
  //        int edgePtIdN = branchPd->GetPointData()->GetArray("TmpInternalIds")->
  //          GetTuple1(splitOpenEdges[k][edgeSize-1]);

  //        int branchPtId0   = splitOpenEdges[k][0];
  //        int branchPtId1   = splitOpenEdges[k][1];
  //        int branchPtIdN   = splitOpenEdges[k][edgeSize-1];
  //        int branchPtIdNm1 = splitOpenEdges[k][edgeSize-2];

  //        vtkNew(vtkIdList, firstCellId);
  //        branchPd->GetCellEdgeNeighbors(-1, branchPtId0, branchPtId1, firstCellId);
  //        if (firstCellId->GetNumberOfIds() != 1)
  //        {
  //          fprintf(stderr,"Something went wrong here\n");
  //          return SV_OK;
  //        }
  //        int realCellId0 = branchPd->GetCellData()->GetArray("TmpInternalIds")->
  //          GetTuple1(firstCellId->GetId(0));
  //        double firstAngularVal = angularPCoords->GetTuple1(realCellId0);

  //        vtkNew(vtkIdList, lastCellId);
  //        branchPd->GetCellEdgeNeighbors(-1, branchPtIdN, branchPtIdNm1, lastCellId);
  //        if (lastCellId->GetNumberOfIds() != 1)
  //        {
  //          fprintf(stderr,"Something went wrong here\n");
  //          return SV_OK;
  //        }
  //        int realCellIdN = branchPd->GetCellData()->GetArray("TmpInternalIds")->
  //          GetTuple1(lastCellId->GetId(0));
  //        double lastAngularVal = angularPCoords->GetTuple1(realCellIdN);

  //        if (edgePtId0 == -1 || edgePtIdN == -1)
  //        {
  //          fprintf(stdout,"Could not recover true ids\n");
  //          return SV_ERROR;
  //        }
  //        std::vector<double> angleBounds;
  //        angleBounds.push_back(firstAngularVal);
  //        angleBounds.push_back(lastAngularVal);
  //        std::sort(angleBounds.begin(), angleBounds.end());
  //        edgeAngleBounds.push_back(angleBounds);

  //        int polyPtId0 = polyBranchPd->GetPointData()->GetArray("SlicePoints")->
  //          LookupValue(edgePtId0);
  //        int polyPtIdN = polyBranchPd->GetPointData()->GetArray("SlicePoints")->
  //          LookupValue(edgePtIdN);

  //        if (polyPtId0 == -1 || polyPtIdN == -1)
  //        {
  //          fprintf(stdout,"Could not recover true ids from polycube\n");
  //          return SV_ERROR;
  //        }

  //        vtkNew(vtkIdList, polyCellId);
  //        vtkNew(vtkIdList, cellPointIds);
  //        cellPointIds->SetNumberOfIds(2);
  //        cellPointIds->SetId(0, polyPtId0);
  //        cellPointIds->SetId(1, polyPtIdN);
  //        polyBranchPd->GetCellNeighbors(-1, cellPointIds, polyCellId);

  //        if (polyCellId->GetNumberOfIds() != 1)
  //        {
  //          fprintf(stdout,"Should have one and only one cell here\n");
  //          return SV_ERROR;
  //        }

  //        int patchVal = polyBranchPd->GetCellData()->GetArray("PatchIds")->
  //          GetTuple1(polyCellId->GetId(0));
  //        patchVal = patchVal%6;

  //        edgePatchValues.push_back(patchVal);

  //        for (int l=0; l<splitOpenEdges[k].size()-1; l++)
  //        {
  //          int splitPtId0 = splitOpenEdges[k][l];
  //          int splitPtId1 = splitOpenEdges[k][l+1];

  //          vtkNew(vtkIdList, splitCellId);
  //          branchPd->GetCellEdgeNeighbors(-1, splitPtId0, splitPtId1, splitCellId);

  //          if (splitCellId->GetNumberOfIds() != 1)
  //          {
  //            fprintf(stderr,"Something went wrong here\n");
  //            return SV_OK;
  //          }

  //          int branchCellId = splitCellId->GetId(0);
  //          edgeCellList.push_back(branchCellId);
  //          cellBool[branchCellId];

  //          int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->
  //            GetTuple1(branchCellId);
  //          this->WorkPd->GetCellData()->GetArray("PatchVals")->SetTuple1(realCellId, patchVal);

  //          double locals[6][3];
  //          centerlineLocalX->GetTuple(realCellId, locals[0]);
  //          centerlineLocalY->GetTuple(realCellId, locals[1]);
  //          centerlineLocalZ->GetTuple(realCellId, locals[4]);
  //          for (int m=0; m<3; m++)
  //          {
  //            locals[2][m] = -1.0*locals[0][m];
  //            locals[3][m] = -1.0*locals[1][m];
  //            locals[5][m] = -1.0*locals[4][m];
  //          }

  //          double boundarySetVec[3];
  //          for (int m=0; m<3; m++)
  //            boundarySetVec[m] = locals[patchVal][m];

  //          vtkMath::Normalize(boundarySetVec);

  //          preRotationNormals->SetTuple(realCellId, boundarySetVec);
  //        }
  //      }
  //      growCellLists.push_back(edgeCellList);
  //      allPatchValues.push_back(edgePatchValues);
  //      allAngleBounds.push_back(edgeAngleBounds);
  //    }

  //    if (growCellLists.size() == 2)
  //    {
  //      if (isTerminating)
  //      {
  //        fprintf(stderr,"Something wrong here, branch is terminating, but we found multiple open edges\n");
  //        return SV_ERROR;
  //      }
  //      double allVals[2]; allVals[0] = 0.0; allVals[1] = 0.0;
  //      for (int j=0; j<growCellLists.size(); j++)
  //      {
  //        for (int k=0; k<growCellLists[j].size(); k++)
  //        {
  //          int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->
  //            GetTuple1(growCellLists[j][k]);
  //          double pCoordVal = centerlinePCoords->GetTuple1(realCellId);
  //          allVals[j] += pCoordVal;
  //        }
  //      }

  //      if (!(allVals[0] < allVals[1]))
  //      {
  //        // Gotta switch
  //        std::vector<int> tmpList1 = growCellLists[0];
  //        std::vector<int> tmpList0 = growCellLists[1];

  //        growCellLists.clear();
  //        growCellLists.push_back(tmpList0);
  //        growCellLists.push_back(tmpList1);

  //        std::vector<std::vector<double> > tmpAngleList1 = allAngleBounds[0];
  //        std::vector<std::vector<double> > tmpAngleList0 = allAngleBounds[1];

  //        allAngleBounds.clear();
  //        allAngleBounds.push_back(tmpAngleList0);
  //        allAngleBounds.push_back(tmpAngleList1);

  //        std::vector<int> tmpPatchList1 = allPatchValues[0];
  //        std::vector<int> tmpPatchList0 = allPatchValues[1];

  //        allPatchValues.clear();
  //        allPatchValues.push_back(tmpPatchList0);
  //        allPatchValues.push_back(tmpPatchList1);
  //      }
  //    }

  //    if (growCellLists.size() > 2)
  //    {
  //      fprintf(stderr,"WE GOT OURSELVES A PROBLEMO\n");
  //      return SV_ERROR;
  //    }

  //    // Find minimum patch val and swap order
  //    for (int j=0; j<allAngleBounds.size(); j++)
  //    {
  //      vtkNew(vtkDoubleArray, tmpSortAngles);
  //      vtkNew(vtkIdList, angleIndices);
  //      for (int k=0; k<allAngleBounds[j].size(); k++)
  //      {
  //        tmpSortAngles->InsertNextTuple1(allAngleBounds[j][k][0]);
  //        angleIndices->InsertNextId(k);
  //      }

  //      vtkSortDataArray::Sort(tmpSortAngles, angleIndices);
  //      std::vector<std::vector<double> > newAngleBounds = allAngleBounds[j];
  //      std::vector<int> newPatchValues = allPatchValues[j];
  //      for (int k=0; k<angleIndices->GetNumberOfIds(); k++)
  //      {
  //        int listIndex = angleIndices->GetId(k);
  //        if (k == 0)
  //        {
  //          if (allAngleBounds[j][listIndex][1] < 0.0)
  //          {
  //            allAngleBounds[j][listIndex][0] = SV_PI;
  //          }
  //          else
  //          {
  //            double tmp = allAngleBounds[j][listIndex][0];
  //            allAngleBounds[j][listIndex][0] = allAngleBounds[j][listIndex][1];
  //            allAngleBounds[j][listIndex][1] = tmp;
  //          }
  //        }

  //        newAngleBounds[k].clear();
  //        newAngleBounds[k] = allAngleBounds[j][listIndex];

  //        newPatchValues[k] = allPatchValues[j][listIndex];
  //      }

  //      //for (int k=0; k<newAngleBounds.size(); k++)
  //      //{
  //      //  int thisSize = newAngleBounds.size();
  //      //  double lastVal  = newAngleBounds[k][1];
  //      //  double firstVal = newAngleBounds[(k+1)%thisSize][0];
  //      //  double avgVal = (lastVal + firstVal)/2.0;

  //      //  newAngleBounds[k][1] = avgVal;
  //      //  newAngleBounds[(k+1)%thisSize][0] = avgVal;
  //      //}

  //      for (int k=0; k<newAngleBounds.size(); k++)
  //      {
  //        allAngleBounds[j].clear();
  //        allAngleBounds[j] = newAngleBounds;

  //        allPatchValues[j] = newPatchValues;

  //      }
  //    }

  //    std::vector<std::vector<int> > cellNeighbors;
  //    std::vector<int> numCellNeighbors;
  //    this->GetCellDirectNeighbors(branchPd, cellNeighbors, numCellNeighbors);

  //    double maxPCoordThr = 1.0*this->BoundaryEnforceFactor/nlinepts;
  //    fprintf(stdout,"BOUNDARY ENFORCE FACTOR: %.6f\n", maxPCoordThr);
  //    double pCoordThr = 0.01;
  //    double begVessel = 0.0;
  //    double endVessel = 1.0;

  //    int iter = 0;
  //    while(begVessel < maxPCoordThr)
  //    {
  //      int done = 0;
  //      while (!done)
  //      {
  //        done = 1;
  //        for (int listIter=0; listIter<growCellLists.size(); listIter++)
  //        {
  //          for (int j=0; j<growCellLists[listIter].size(); j++)
  //          {
  //            int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->
  //              GetTuple1(growCellLists[listIter][j]);

  //            int linePtId = centerlineSubPtIds->GetTuple1(realCellId);

  //            int patchVal = this->WorkPd->GetCellData()->GetArray("PatchVals")->GetTuple1(realCellId);
  //            // THIS IS WHERE WE TRY TO USE CELL ANGULAR LOCATION FOR PATCH VAL
  //            double angularVal = angularPCoords->GetTuple1(realCellId);

  //            //int patchVal = -1;
  //            for (int k=0; k<allAngleBounds[listIter].size(); k++)
  //            {
  //              int thisSize = allAngleBounds[listIter].size();
  //              if (k == 0)
  //              {
  //                if (angularVal >= allAngleBounds[listIter][k][0] ||
  //                    angularVal <=  allAngleBounds[listIter][k][1])
  //                {
  //                  patchVal = allPatchValues[listIter][k];
  //                  break;
  //                }
  //                //else if (angularVal > allAngleBounds[listIter][k][1] &&
  //                //         angularVal < allAngleBounds[listIter][(k+1)%thisSize][0])
  //                //{
  //                //  fprintf(stdout,"NEED TO PICK ONE!!\n");
  //                //  patchVal = allPatchValues[listIter][k];
  //                //  break;
  //                //}
  //              }
  //              else
  //              {
  //                if (angularVal >= allAngleBounds[listIter][k][0] &&
  //                    angularVal <=  allAngleBounds[listIter][k][1])
  //                {
  //                  patchVal = allPatchValues[listIter][k];
  //                  break;
  //                }
  //                //else if (angularVal > allAngleBounds[listIter][k][1] &&
  //                //         angularVal < allAngleBounds[listIter][(k+1)%thisSize][0])
  //                //{
  //                //  fprintf(stdout,"NEED TO PICK ONE!!\n");
  //                //  patchVal = allPatchValues[listIter][k];
  //                //  break;
  //                //}
  //              }
  //            }

  //            if (patchVal == -1)
  //            {
  //              fprintf(stderr,"A PATCH VAL NOT FOUND BY ANGULAR METHOD\n");
  //              return SV_ERROR;
  //            }

  //            for (int k=0; k<numCellNeighbors[growCellLists[listIter][j]]; k++)
  //            {
  //              int neighborCellId = cellNeighbors[growCellLists[listIter][j]][k];

  //              int neighborRealCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->
  //                GetTuple1(neighborCellId);

  //              double pCoordVal = centerlinePCoords->GetTuple1(neighborRealCellId);

  //              if (cellBool[neighborCellId] == 0)
  //              {
  //                if ((pCoordVal >= begVessel && pCoordVal <= begVessel + pCoordThr && listIter == 0) ||
  //                    (pCoordVal <= endVessel && pCoordVal >= endVessel - pCoordThr && listIter == 1))
  //                {
  //                  done = 0;
  //                  growCellLists[listIter].push_back(neighborCellId);
  //                  cellBool[neighborCellId] = 1;

  //                  double currentVec[3];
  //                  preRotationNormals->GetTuple(neighborRealCellId, currentVec);

  //                  double locals[6][3];
  //                  centerlineLocalX->GetTuple(neighborRealCellId, locals[0]);
  //                  centerlineLocalY->GetTuple(neighborRealCellId, locals[1]);
  //                  centerlineLocalZ->GetTuple(neighborRealCellId, locals[4]);
  //                  for (int l=0; l<3; l++)
  //                  {
  //                    locals[2][l] = -1.0*locals[0][l];
  //                    locals[3][l] = -1.0*locals[1][l];
  //                    locals[5][l] = -1.0*locals[4][l];
  //                  }

  //                  //double beta = (1.0*row)/numRows;
  //                  double beta = begVessel/maxPCoordThr;

  //                  double boundarySetVec[3];
  //                  for (int l=0; l<3; l++)
  //                    boundarySetVec[l] = beta * currentVec[l] +  (1 - beta) * locals[patchVal][l];

  //                  preRotationNormals->SetTuple(neighborRealCellId, boundarySetVec);
  //                  this->WorkPd->GetCellData()->GetArray("PatchVals")->SetTuple1(neighborRealCellId, patchVal);
  //                }
  //              }
  //            }
  //          }
  //        }
  //      }

  //      begVessel += pCoordThr;
  //      endVessel -= pCoordThr;
  //    }

  //    // Now go through and transform to local coordinate system and set
  //    // the new vector to use for clustering
  //    for (int j=0; j<branchNumberOfCells; j++)
  //    {
  //      //Get real cell id
  //      int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);

  //      double locals[6][3];
  //      centerlineLocalX->GetTuple(realCellId, locals[0]);
  //      centerlineLocalY->GetTuple(realCellId, locals[1]);
  //      centerlineLocalZ->GetTuple(realCellId, locals[4]);
  //      for (int k=0; k<3; k++)
  //      {
  //        locals[2][k] = -1.0*locals[0][k];
  //        locals[3][k] = -1.0*locals[1][k];
  //        locals[5][k] = -1.0*locals[4][k];
  //      }

  //      // Compute the rotation from global coordinate system to centerlines
  //      // local coordinate system
  //      double rotMat[9];
  //      this->ComputeRotationMatrix(locals[0], locals[1], locals[4], rotMat);

  //      double cellClusterVec[3];
  //      preRotationNormals->GetTuple(realCellId, cellClusterVec);

  //      // Apply rotation matrix to the normal to get the new normal
  //      double newNormal[3];
  //      for (int k=0; k<3; k++)
  //      {
  //        newNormal[k] = rotMat[k*3]*cellClusterVec[0] +
  //                       rotMat[(k*3)+1]*cellClusterVec[1] +
  //                       rotMat[(k*3)+2]*cellClusterVec[2];
  //      }

  //      centerlineBasedNormals->SetTuple(realCellId, newNormal);
  //    }
  //  }
  //}

  //this->WorkPd->GetCellData()->AddArray(centerlineBasedNormals);
  //this->WorkPd->GetCellData()->AddArray(centerlineSubPtIds);
  //this->WorkPd->GetCellData()->AddArray(centerlinePCoords);
  //this->WorkPd->GetCellData()->AddArray(angularPCoords);

  //// CLUSTERING
  //// Set up generators
  //vtkNew(vtkPoints, generatorsPts);
  //generatorsPts->SetNumberOfPoints(6);
  //generatorsPts->SetPoint(0, 1.0, 0.0, 0.0);
  //generatorsPts->SetPoint(1, 0.0, 1.0, 0.0);
  //generatorsPts->SetPoint(2, -1.0, 0.0, 0.0);
  //generatorsPts->SetPoint(3, 0.0, -1.0, 0.0);
  //generatorsPts->SetPoint(4, 0.0, 0.0, 1.0);
  //generatorsPts->SetPoint(5, 0.0, 0.0, -1.0);

  //vtkNew(vtkPolyData, generatorsPd);
  //generatorsPd->SetPoints(generatorsPts);

  //vtkIntArray *tmpPatchArray = vtkIntArray::New();
  //tmpPatchArray->SetNumberOfTuples(this->WorkPd->GetNumberOfCells());
  //tmpPatchArray->SetName("PatchIds");
  //tmpPatchArray->FillComponent(0, -1);
  //this->WorkPd->GetCellData()->AddArray(tmpPatchArray);
  //tmpPatchArray->Delete();

  //vtkSVGeneralUtils::GiveIds(this->PolycubePd, "TmpInternalIds");

  //for (int i=0; i<numGroups; i++)
  //{
  //  int groupId = groupIds->GetId(i);

  //  fprintf(stdout,"CLUSTERING AND MATCHING ENDS OF %d\n", groupId);

  //  vtkNew(vtkPolyData, branchPd);
  //  vtkSVGeneralUtils::ThresholdPd(this->WorkPd, groupId, groupId, 1,
  //      this->GroupIdsArrayName, branchPd);
  //  branchPd->BuildLinks();

  //  vtkNew(vtkPolyData, polyBranchPd);
  //  vtkSVGeneralUtils::ThresholdPd(this->PolycubePd, groupId, groupId, 1,
  //    this->GroupIdsArrayName, polyBranchPd);
  //  polyBranchPd->BuildLinks();

  //  if (this->RunEdgeWeightedCVT(branchPd, generatorsPd) != SV_OK)
  //  {
  //    vtkErrorMacro("Error in cvt");
  //    return SV_ERROR;
  //  }

  //  if (this->MergedCenterlines->GetNumberOfCells() > 1)
  //  {
  //    if (this->FixEndPatches(branchPd) != SV_OK)
  //    {
  //      vtkErrorMacro("Error fixing end patches");
  //      return SV_ERROR;
  //    }
  //  }

  //  if (this->MergedCenterlines->GetNumberOfCells() > 1)
  //  {
  //    if (this->FixSidePatches(branchPd) != SV_OK)
  //    {
  //      vtkErrorMacro("Error fixing side patches");
  //      return SV_ERROR;
  //    }
  //  }

  //  vtkNew(vtkIdList, noEndPatches);
  //  noEndPatches->SetNumberOfIds(4);
  //  for (int j=0; j<4; j++)
  //    noEndPatches->SetId(j, j);

  //  if (this->CorrectSpecificCellBoundaries(branchPd, "PatchIds", noEndPatches) != SV_OK)
  //  {
  //    vtkErrorMacro("Could not correcto boundaries of surface");
  //    return SV_ERROR;
  //  }

  //  if (this->MergedCenterlines->GetNumberOfCells() > 1)
  //  {
  //    if (this->MatchEndPatches(branchPd, polyBranchPd) != SV_OK)
  //    {
  //      vtkErrorMacro("Error matching end patches");
  //      return SV_ERROR;
  //    }
  //  }

  //  // Set vals on work pd
  //  for (int j=0; j<branchPd->GetNumberOfCells(); j++)
  //  {
  //    //Get real cell id
  //    int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(j);

  //    // Get val
  //    int cellVal = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(j);

  //    // Set val
  //    this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(realCellId, cellVal);
  //  }
  //}

  //this->WorkPd->GetCellData()->RemoveArray("TmpInternalIds");
  //this->WorkPd->GetPointData()->RemoveArray("TmpInternalIds");

  //this->PolycubePd->GetCellData()->RemoveArray("TmpInternalIds");
  //this->PolycubePd->GetPointData()->RemoveArray("TmpInternalIds");

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

  //std::vector<Region> finalRegions;
  //vtkNew(vtkIdList, targetPatches);
  //targetPatches->SetNumberOfIds(numGroups*4);
  //for (int i=0; i<numGroups; i++)
  //{
  //  for (int j=0; j<4; j++)
  //    targetPatches->SetId(4*i+j, 6*i+j);
  //}

  ////// TODO: IF SOMETHIGN WRONG, LOOK HERE FIRST!!! MAY BE MOVING PATCH OFF
  ////// OF SLICE POINT
  ////if (this->CorrectSpecificCellBoundaries(this->WorkPd, "PatchIds", targetPatches) != SV_OK)
  ////{
  ////  vtkErrorMacro("Could not correcto boundaries of surface");
  ////  return SV_ERROR;
  ////}

  //// For checking purposes
  //if (this->FixPatchesWithPolycube() != SV_OK)
  //{
  //  fprintf(stderr,"Couldn't fix patches\n");
  //  return SV_ERROR;
  //}

  //if (this->CorrectSpecificCellBoundaries(this->WorkPd, "PatchIds", targetPatches) != SV_OK)
  //{
  //  vtkErrorMacro("Could not correcto boundaries of surface");
  //  return SV_ERROR;
  //}

  //if (this->SmoothSpecificBoundaries(this->WorkPd, "PatchIds", targetPatches) != SV_OK)
  //{
  //  vtkErrorMacro("Could not smootho boundaries of surface");
  //  return SV_ERROR;
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

  ////////////////// For checking purposes
  ////////////////if (this->FixPatchesWithPolycubeOld() != SV_OK)
  ////////////////{
  ////////////////  fprintf(stderr,"Couldn't fix patches\n");
  ////////////////  return SV_ERROR;
  ////////////////}

  ////// NOW PARAMETERIZE!!, WIILL BE MOVED to vtkSVPolycubeParameterizer
  ////// TODO: RENAME THIS CLASS TO vtkSVCenterlinesSegmenter

  //vtkNew(vtkPolyData, fullMapPd);
  //if (this->ParameterizeSurface(fullMapPd) != SV_OK)
  //{
  //  fprintf(stderr,"WRONG\n");
  //  return SV_ERROR;
  //}

  //vtkNew(vtkUnstructuredGrid, loftedVolume);
  //if (this->ParameterizeVolume(fullMapPd, loftedVolume) != SV_OK)
  //{
  //  fprintf(stderr,"Failed doing volume stuffs\n");
  //  return SV_ERROR;
  //}

  //std::string fn = "/Users/adamupdegrove/Desktop/tmp/FINAL_NURBS.vtu";
  //vtkSVIOUtils::WriteVTUFile(fn, loftedVolume);

  return SV_OK;
}

// ----------------------
// GetApproximatePolycubeSize
// ----------------------
int vtkSVGroupsSegmenter::GetApproximatePolycubeSize(double &polycubeSize)
{
  double avgRadius = 0.0;

  int numPoints = this->Centerlines->GetNumberOfPoints();
  vtkDataArray *radiusArray = this->Centerlines->GetPointData()->GetArray(
    this->CenterlineRadiusArrayName);

  for (int i=0; i<numPoints; i++)
    avgRadius += radiusArray->GetTuple1(i);

  avgRadius = avgRadius/numPoints;

  polycubeSize = 2*avgRadius;
  polycubeSize = avgRadius;

  return SV_OK;
}

// ----------------------
// MergeCenterlines
// ----------------------
int vtkSVGroupsSegmenter::MergeCenterlines()
{
  if (vtkSVGeneralUtils::CheckArrayExists(this->Centerlines, 1, this->GroupIdsArrayName) != SV_OK)
  {
    std::cout<<"Splitting centerlines..."<<endl;
    vtkNew(vtkSVCenterlineBranchSplitter, branchSplitter);
    branchSplitter->SetInputData(this->Centerlines);
    branchSplitter->SetGroupingModeToFirstPoint();
    branchSplitter->SetBlankingArrayName("Blanking");
    branchSplitter->SetRadiusArrayName(this->CenterlineRadiusArrayName);
    branchSplitter->SetGroupIdsArrayName("GroupIds");
    branchSplitter->SetCenterlineIdsArrayName("CenterlineIds");
    branchSplitter->SetTractIdsArrayName("TractIds");
    branchSplitter->SetRadiusMergeRatio(this->RadiusMergeRatio);
    branchSplitter->SetUseAbsoluteMergeDistance(this->UseAbsoluteMergeDistance);
    branchSplitter->SetMergeDistance(this->MergeDistance);
    branchSplitter->Update();

    this->Centerlines->DeepCopy(branchSplitter->GetOutput());
  }

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

  int numFullPts = merger->GetOutput()->GetNumberOfPoints();

  vtkNew(vtkCleanPolyData, lineCleaner);
  lineCleaner->SetInputData(merger->GetOutput());
  lineCleaner->Update();

  this->MergedCenterlines->DeepCopy(lineCleaner->GetOutput());
  this->MergedCenterlines->BuildLinks();

  if (!this->IsVasculature)
  {
    int numRemove = this->NumberOfCenterlineRemovePts;
    vtkNew(vtkPoints, newPoints);
    vtkNew(vtkPointData, newPointData);
    newPointData->CopyAllocate(this->MergedCenterlines->GetPointData(),
                               numFullPts);

    vtkNew(vtkCellArray, newCells);
    vtkNew(vtkCellData, newCellData);
    newCellData->CopyAllocate(this->MergedCenterlines->GetCellData());

    for (int i=0; i<this->MergedCenterlines->GetNumberOfCells(); i++)
    {
      vtkIdType npts, *pts;
      this->MergedCenterlines->GetCellPoints(i, npts, pts);

      vtkNew(vtkIdList, point0CellIds);
      this->MergedCenterlines->GetPointCells(pts[0], point0CellIds);

      vtkNew(vtkIdList, pointNCellIds);
      this->MergedCenterlines->GetPointCells(pts[npts-1], pointNCellIds);

      vtkNew(vtkPolyLine, newLine);
      if (point0CellIds->GetNumberOfIds() > 1)
      {
        for (int j=0; j<numRemove; j++)
        {
          int newPointId = newPoints->InsertNextPoint(
            this->MergedCenterlines->GetPoint(pts[j]));

          newLine->GetPointIds()->InsertNextId(newPointId);

          newPointData->CopyData(this->MergedCenterlines->GetPointData(),
            pts[j], newPointId);
        }
      }

      for (int j=numRemove; j<npts-numRemove; j++)
      {
        int newPointId = newPoints->InsertNextPoint(
          this->MergedCenterlines->GetPoint(pts[j]));
        newLine->GetPointIds()->InsertNextId(newPointId);

        newPointData->CopyData(this->MergedCenterlines->GetPointData(),
          pts[j], newPointId);
      }

      if (pointNCellIds->GetNumberOfIds() > 1)
      {
        for (int j=numRemove; j>0; j--)
        {
          int newPointId = newPoints->InsertNextPoint(
            this->MergedCenterlines->GetPoint(pts[npts-j]));
          newLine->GetPointIds()->InsertNextId(newPointId);

          newPointData->CopyData(this->MergedCenterlines->GetPointData(),
            pts[npts-j], newPointId);
        }
      }

      newCells->InsertNextCell(newLine);
      newCellData->CopyData(this->MergedCenterlines->GetCellData(), i, i);
    }

    this->MergedCenterlines->Reset();
    this->MergedCenterlines->SetPoints(newPoints);
    this->MergedCenterlines->SetLines(newCells);

    newPointData->Squeeze();
    this->MergedCenterlines->GetPointData()->PassData(newPointData);
    this->MergedCenterlines->GetCellData()->PassData(newCellData);

    vtkNew(vtkCleanPolyData, cleaner);
    cleaner->SetInputData(this->MergedCenterlines);
    cleaner->Update();

    this->MergedCenterlines->DeepCopy(cleaner->GetOutput());
    this->MergedCenterlines->BuildLinks();
  }

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
  CVT->SetUseCurvatureWeight(0);
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
      if (neiSize == 2)
      {
        if (neiTmpIds->GetId(0) == neiTmpIds->GetId(1))
        {
          checkFix = 1;
          int maxVal, maxCount;
          vtkSVGroupsSegmenter::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
      }
      else if (neiSize >= 3)
      {
        if ((neiTmpIds->GetId(0) == neiTmpIds->GetId(1) ||
             neiTmpIds->GetId(1) == neiTmpIds->GetId(2)))
        {
          checkFix = 1;
          int maxVal, maxCount;
          vtkSVGroupsSegmenter::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
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
      if (neiSize == 2)
      {
        if (neiTmpIds->GetId(0) == neiTmpIds->GetId(1))
        {
          checkFix = 1;
          int maxVal, maxCount;
          vtkSVGroupsSegmenter::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
      }
      else if (neiSize >= 3)
      {
        if ((neiTmpIds->GetId(0) == neiTmpIds->GetId(1) ||
             neiTmpIds->GetId(1) == neiTmpIds->GetId(2)))
        {
          checkFix = 1;
          int maxVal, maxCount;
          vtkSVGroupsSegmenter::GetMostOccuringVal(neiCellIds, maxVal, maxCount);

          cellIds->SetTuple1(i, maxVal);
          tmpIds->SetTuple1(i, neiTmpIds->GetId(1));
        }
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
    //fprintf(stdout,"NUM CORNS: %d OF GROUP %d\n", allRegions[i].NumberOfCorners, allRegions[i].IndexCluster);

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

// TODO: Need to fix for if single cell is region!!!
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

  position = (position+1)%npts;
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

  position = (position+npts-1)%npts;
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

  if (individualFix.size() == 0 && wholePatchFix.size() >= 1)
  {
    fprintf(stdout,"NO INDIVIDUAL FIX ENDS, THAT MEANS EITHER WE HAVE A GOOD ONE ALREADY OR THE BAD WHOLE PATCH IS THE ONE, CHECK\n");
    int badPatch;
    for (int b=0; b<wholePatchFix.size(); b++)
    {
      fprintf(stdout,"ITS A WHOLE PATCH\n");
      badPatch = wholePatchFix[b];
      vtkNew(vtkIdList, neighborPatchIds);
      for (int i=0; i<endRegions[badPatch].NumberOfElements; i++)
      {
        int cellId = endRegions[badPatch].Elements[i];

        vtkIdType npts, *pts;
        pd->GetCellPoints(cellId, npts, pts);

        for (int j=0; j<npts; j++)
        {
          int ptId0 = pts[j];
          int ptId1 = pts[(j+1)%npts];

          vtkNew(vtkIdList, cellEdgeNeighbors);
          pd->GetCellEdgeNeighbors(cellId, ptId0, ptId1, cellEdgeNeighbors);

          for (int k=0; k<cellEdgeNeighbors->GetNumberOfIds(); k++)
          {
            int cellPatchId = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(cellEdgeNeighbors->GetId(k));

            if (cellPatchId != endRegions[badPatch].IndexCluster)
              neighborPatchIds->InsertUniqueId(cellPatchId);
          }
        }
      }

      fprintf(stdout,"NUM IDS: %d\n", neighborPatchIds->GetNumberOfIds());
      for (int j=0; j<neighborPatchIds->GetNumberOfIds(); j++)
        fprintf(stdout,"  ID: %d\n", neighborPatchIds->GetId(j));
      if (neighborPatchIds->GetNumberOfIds() == 4 &&
          neighborPatchIds->IsId(0) != -1 && neighborPatchIds->IsId(1) != -1 &&
          neighborPatchIds->IsId(2) != -1 && neighborPatchIds->IsId(3) != -1 &&
          endRegions[badPatch].NumberOfCorners == 4)
      {
        // This region is okay because it has four corners and it touches all four side group ids
        fprintf(stdout,"WE FOUND OUR GOOD REGION: %d\n", badPatch);
        wholePatchFix.erase(std::remove(wholePatchFix.begin(), wholePatchFix.end(), badPatch), wholePatchFix.end());
        break;
      }
    }
  }

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
      fprintf(stdout,"FIX STRATEGY 0\n");
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
    else if (fixStrategy == 1)
    {
      fprintf(stdout,"FIX STRATEGY 1\n");
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
        if (sideRegions[minPatch].NumberOfElements == 1 ||
            sideRegions[minPatch].BoundaryEdges.size() == 0)
        {
          for (int j=0; j<sideRegions[minPatch].Elements.size(); j++)
          {
            int cellId = sideRegions[minPatch].Elements[j];
            vtkIdType npts, *pts;
            pd->GetCellPoints(cellId, npts, pts);
            for (int k=0; k<npts; k++)
            {
              int ptId0 = pts[k];
              int ptId1 = pts[(k+1)%npts];

              vtkNew(vtkIdList, cellEdgeNeighbors);
              pd->GetCellEdgeNeighbors(cellId, ptId0, ptId1, cellEdgeNeighbors);

              for (int l=0; l<cellEdgeNeighbors->GetNumberOfIds(); l++)
              {
                int edgeCellId  = cellEdgeNeighbors->GetId(l);
                int cellVal = pd->GetCellData()->GetArray("PatchIds")->GetTuple1(edgeCellId);

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
        }

        int maxVal = -1;
        int maxPatchId = -1;
        for (int j=0; j<patchIds->GetNumberOfIds(); j++)
        {
          if (patchCount->GetId(j) > maxVal)
          {
            maxPatchId = patchIds->GetId(j);
            maxVal = patchCount->GetId(j);
          }
        }
        if (maxPatchId == -1)
        {
          fprintf(stderr,"A patch value to change bad patch to was not found\n");
          return SV_ERROR;
        }

        for (int k=0; k<sideRegions[minPatch].Elements.size(); k++)
        {
          int cellId = sideRegions[minPatch].Elements[k];

          pd->GetCellData()->GetArray("PatchIds")->SetTuple1(cellId, maxPatchId);
        }
      }
    }
  }

  if (this->GetSpecificRegions(pd, "PatchIds", sideRegions, targetRegions) != SV_OK)
  {
    vtkErrorMacro("Couldn't get patches");
    return SV_ERROR;
  }

  // Get open boundary edges
  std::vector<int> openCornerPoints;
  std::vector<std::vector<int> > openEdges;

  if (this->GetOpenBoundaryEdges(pd, sideRegions, "PatchIds",
                                 openCornerPoints, openEdges) != SV_OK)
  {
    fprintf(stderr,"Error getting open boundary edges\n");
    return SV_ERROR;
  }

  std::vector<std::vector<int> > connectedOpenCornerPts;
  this->GetConnectedEdges(openEdges, connectedOpenCornerPts);

  for (int i=0; i<connectedOpenCornerPts.size(); i++)
  {
    if (connectedOpenCornerPts[i].size() != 4)
    {
      fprintf(stderr,"THERE ARE NOT FOUR CORNER POINTS ON THIS PATCH, THIS IS QUITE BAD\n");
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
// CheckGroups2
// ----------------------
int vtkSVGroupsSegmenter::CheckGroups2()
{
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
  vtkNew(vtkFeatureEdges, featureEdges);

  vtkIdType nlinepts, *linepts;
  int groupId, centerlineId, isTerminating;

  vtkPolyData *newMergedCenterlinesPd = NULL;

  for (int i=0; i<numGroups; i++)
  {
    groupId = groupIds->GetId(i);

    centerlineId = this->MergedCenterlines->GetCellData()->GetArray(this->GroupIdsArrayName)->LookupValue(groupId);
    this->MergedCenterlines->GetCellPoints(centerlineId, nlinepts, linepts);
    isTerminating = 0;
    this->MergedCenterlines->GetPointCells(linepts[nlinepts-1], backNeighbors);
    vtkNew(vtkIdList, backGroupNeighbors);
    for (int j=0; j<backNeighbors->GetNumberOfIds(); j++)
    {
      backGroupNeighbors->InsertNextId(this->MergedCenterlines->GetCellData()->GetArray(
        this->GroupIdsArrayName)->GetTuple1(backNeighbors->GetId(j)));
    }
    if (backNeighbors->GetNumberOfIds() == 1)
      isTerminating = 1;

    groupThresholder->ThresholdBetween(groupId, groupId);
    groupThresholder->Update();

    if (groupThresholder->GetOutput()->GetNumberOfPoints() == 0)
    {
      vtkErrorMacro("THERE ARE NO CELLS ON SURFACE FOR GROUP "<< groupId);
      continue;
    }

    surfacer->SetInputData(groupThresholder->GetOutput());
    surfacer->Update();

    connector->SetInputData(surfacer->GetOutput());
    connector->Update();

    if (connector->GetNumberOfExtractedRegions() > 1)
    {
      vtkErrorMacro("EXPECTED ONE EXTRACTED REGION FOR GROUP "<< groupId << ", BUT THERE ARE " << connector->GetNumberOfExtractedRegions() << " REGIONS");
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

      vtkNew(vtkIdList, frontNeighbors);
      this->MergedCenterlines->GetPointCells(linepts[0], frontNeighbors);
      vtkNew(vtkIdList, frontGroupNeighbors);
      for (int j=0; j<frontNeighbors->GetNumberOfIds(); j++)
      {
        frontGroupNeighbors->InsertNextId(this->MergedCenterlines->GetCellData()->GetArray(
          this->GroupIdsArrayName)->GetTuple1(frontNeighbors->GetId(j)));
      }

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
      CVT->SetUseCurvatureWeight(1);
      CVT->SetPatchIdsArrayName(this->GroupIdsArrayName);
      CVT->SetCVTDataArrayName("Normals");
      CVT->SetGroupIdsArrayName(this->GroupIdsArrayName);
      CVT->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
      CVT->SetBlankingArrayName(this->BlankingArrayName);
      //CVT->SetUseBifurcationInformation(1);
      CVT->SetUsePointNormal(1);
      CVT->SetUseRadiusInformation(0);
      CVT->SetUseBifurcationInformation(0);
      //CVT->SetUsePointNormal(0);
      CVT->SetMaximumNumberOfIterations(0);
      CVT->Update();

      testPd0->DeepCopy(CVT->GetOutput());

      std::vector<Region> groupRegions;
      if (this->GetRegions(testPd0, this->GroupIdsArrayName, groupRegions) != SV_OK)
      {
        vtkErrorMacro("Couldn't get group regions");
        return SV_ERROR;
      }

      int ptId0, ptId1, cellId0, cellId1, cellGroupId0, cellGroupId1;
      double newPt0[3], newPt1[3], avgPt[3];
      std::vector<int> usedPoints(testPd0->GetNumberOfPoints(), 0);
      vtkNew(vtkIdList, cellEdgeNeighbors);
      vtkNew(vtkPoints, newCenterlinePts);
      for (int j=0; j<groupRegions.size(); j++)
      {
        for (int k=0; k<groupRegions[j].BoundaryEdges.size(); k++)
        {
          for (int l=0; l<groupRegions[j].BoundaryEdges[k].size()-1; l++)
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
          }
        }
      }

      double newCenterPt[3];
      int centerPtId = linepts[nlinepts/2];
      this->MergedCenterlines->GetPoint(centerPtId, newCenterPt);

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

      for (int j=0; j<newCenterlinePts->GetNumberOfPoints(); j++)
      {
        ptId0 = newPoints->InsertNextPoint(newCenterPt);
        ptId1 = newPoints->InsertNextPoint(newCenterlinePts->GetPoint(j));

        vtkNew(vtkPolyLine, newPts0);
        newPts0->GetPointIds()->InsertNextId(ptId0);
        newPts0->GetPointIds()->InsertNextId(ptId1);
        newPointData->CopyData(tmpCenterlinesPd->GetPointData(), centerPtId, ptId0);
        newPointData->CopyData(tmpCenterlinesPd->GetPointData(), centerPtId, ptId1);
        newPointData->GetArray(this->CenterlineRadiusArrayName)->SetTuple1(ptId1, 0.0);

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

    featureEdges->SetInputData(surfacer->GetOutput());
    featureEdges->BoundaryEdgesOn();
    featureEdges->FeatureEdgesOff();
    featureEdges->ManifoldEdgesOff();
    featureEdges->NonManifoldEdgesOff();
    featureEdges->Update();

    connector->SetInputData(featureEdges->GetOutput());
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
        vtkErrorMacro("EXPECTED TWO EDGES ON NON-TERMINATING GROUP "<< groupId << ", BUT THERE ARE " << connector->GetNumberOfExtractedRegions() << " EDGES");
        //continue;
      }
    }
  }

  if (newMergedCenterlinesPd != NULL)
  {
    // Now re-segment with these new centerlines

    vtkSVIOUtils::WriteVTPFile("/Users/adamupdegrove/Desktop/tmp/NEWMERGECENTER.vtp", newMergedCenterlinesPd);
    int stopCellNumber = ceil(this->WorkPd->GetNumberOfCells()*0.0001);
    vtkNew(vtkSVCenterlinesEdgeWeightedCVT, betterCVT);
    betterCVT->SetInputData(this->WorkPd);
    betterCVT->SetGenerators(newMergedCenterlinesPd);
    betterCVT->SetNumberOfRings(2);
    betterCVT->SetThreshold(stopCellNumber);
    betterCVT->SetUseCurvatureWeight(1);
    betterCVT->SetPatchIdsArrayName(this->GroupIdsArrayName);
    betterCVT->SetCVTDataArrayName("Normals");
    betterCVT->SetGroupIdsArrayName(this->GroupIdsArrayName);
    betterCVT->SetCenterlineRadiusArrayName(this->CenterlineRadiusArrayName);
    betterCVT->SetBlankingArrayName(this->BlankingArrayName);
    //betterCVT->SetUseBifurcationInformation(1);
    betterCVT->SetUsePointNormal(1);
    betterCVT->SetUseRadiusInformation(0);
    betterCVT->SetUseBifurcationInformation(0);
    //betterCVT->SetUsePointNormal(0);
    betterCVT->SetMaximumNumberOfIterations(0);
    betterCVT->Update();

    this->WorkPd->DeepCopy(betterCVT->GetOutput());

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

    newMergedCenterlinesPd->Delete();
  }

  return SV_OK;
}

// ----------------------
// FixPatchesWithPolycube
// ----------------------
int vtkSVGroupsSegmenter::FixPatchesWithPolycube()
{
  fprintf(stdout,"CHECKING AND FIXING PATCHES...\n");
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
            fprintf(stdout,"DELETING CELLS: %d %d\n", i, cellIds->GetId(k));
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
    fprintf(stdout,"RE-LINKING\n");
    polycubePd->RemoveDeletedCells();

    vtkNew(vtkCleanPolyData, cleaner);
    cleaner->SetInputData(polycubePd);
    cleaner->ToleranceIsAbsoluteOn();
    cleaner->SetAbsoluteTolerance(1.0e-6);
    cleaner->Update();

    polycubePd->DeepCopy(cleaner->GetOutput());;
    polycubePd->BuildLinks();
  }

  fprintf(stdout,"GETTING SURFACE REGIONS\n");
  this->WorkPd->BuildLinks();
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
    vtkErrorMacro("Number of surface patches: " << numSurfacePatches );
    vtkErrorMacro("Number of polycube patches: " << numPolycubePatches );
    if (numSurfacePatches > numPolycubePatches)
    {
      vtkNew(vtkIdList, checkPatchValues);
      for (int i=0; i<surfacePatches.size(); i++)
      {
        if (checkPatchValues->IsId(surfacePatches[i].IndexCluster) != -1)
          fprintf(stderr,"Multiple patches with value %d\n", surfacePatches[i].IndexCluster);
        else
          checkPatchValues->InsertNextId(surfacePatches[i].IndexCluster);
      }
    }
    return SV_ERROR;
  }

  std::vector<int> badPatches;
  std::vector<int> polycubePatchIds;
  for (int i=0; i<numSurfacePatches; i++)
  {
    int patchVal = surfacePatches[i].IndexCluster;
    int patchDir = patchVal%6;

    for (int j=0; j<numPolycubePatches; j++)
    {
      if (polycubePatches[j].IndexCluster == patchVal)
      {
        if (surfacePatches[i].NumberOfCorners != polycubePatches[j].NumberOfCorners)
        {
          badPatches.push_back(i);
          polycubePatchIds.push_back(j);
        }
      }
    }
  }

  fprintf(stdout,"NUMBER OF BAD PATCHES!: %d\n", badPatches.size());
  for (int i=0; i<badPatches.size(); i++)
  {
    int patchId = badPatches[i];

    for (int j=0; j<surfacePatches[patchId].NumberOfCorners; j++)
    {
      int ptId = surfacePatches[patchId].CornerPoints[j];
      fprintf(stdout,"PT ID: %d\n", ptId);
      int polycubeId = polycubePd->GetPointData()->GetArray("SlicePoints")->LookupValue(ptId);
      fprintf(stdout,"POLY ID: %d\n", polycubeId);
      if (polycubeId != -1)
      {
        vtkNew(vtkIdList, surfacePatchVals);
        vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, "PatchIds", ptId, surfacePatchVals);

        vtkNew(vtkIdList, polyPatchVals);
        vtkSVGeneralUtils::GetPointCellsValues(polycubePd, "PatchIds", polycubeId, polyPatchVals);

        vtkNew(vtkIdList, intersectList);
        intersectList->DeepCopy(polyPatchVals);

        intersectList->IntersectWith(surfacePatchVals);

        if (surfacePatchVals->GetNumberOfIds() != intersectList->GetNumberOfIds() ||
            polyPatchVals->GetNumberOfIds() != intersectList->GetNumberOfIds())
        {
          fprintf(stdout,"WE FOUND A BAD ONE!!!!! %d\n", ptId);
          fprintf(stdout,"SURFACE CORNER POINT %d GROUPS ARE ", j);
          for (int l=0; l<surfacePatchVals->GetNumberOfIds(); l++)
            fprintf(stdout,"%d ", surfacePatchVals->GetId(l));
          fprintf(stdout,"\n");
          fprintf(stdout,"POLYCUBE CORNER POINT %d GROUPS ARE ", j);
          for (int l=0; l<polyPatchVals->GetNumberOfIds(); l++)
            fprintf(stdout,"%d ", polyPatchVals->GetId(l));
          fprintf(stdout,"\n");

          vtkNew(vtkIdList, missingVals);
          for (int k=0; k<polyPatchVals->GetNumberOfIds(); k++)
          {
            if (surfacePatchVals->IsId(polyPatchVals->GetId(k)) == -1)
              missingVals->InsertNextId(polyPatchVals->GetId(k));
          }

          for (int m=0; m<missingVals->GetNumberOfIds(); m++)
          {
            for (int k=0; k<surfacePatches[patchId].BoundaryEdges.size(); k++)
            {
              int edgeSize = surfacePatches[patchId].BoundaryEdges[k].size();

              int ptId0 = surfacePatches[patchId].BoundaryEdges[k][0];
              int ptIdN = surfacePatches[patchId].BoundaryEdges[k][edgeSize-1];

              if (ptId0 == ptId)
              {
                vtkNew(vtkIdList, ptIdNList);
                vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, "PatchIds", ptIdN, ptIdNList);

                if (ptIdNList->IsId(missingVals->GetId(m)) != -1)
                {
                  fprintf(stdout,"FIXING WITH EDGE STARTING WITH %d AND ENDING WITH %d.\n", ptId0, ptIdN);
                  for (int l=1; l<surfacePatches[patchId].BoundaryEdges[k].size(); l++)
                  {
                    int testPtId = surfacePatches[patchId].BoundaryEdges[k][l];

                    vtkNew(vtkIdList, pointCells);
                    this->WorkPd->GetPointCells(testPtId, pointCells);

                    for (int r=0; r<pointCells->GetNumberOfIds(); r++)
                    {
                      int cellVal = this->WorkPd->GetCellData()->GetArray("PatchIds")->GetTuple1(pointCells->GetId(r));
                      if (cellVal == surfacePatches[patchId].IndexCluster)
                      {
                        this->WorkPd->GetCellData()->GetArray("PatchIds")->SetTuple1(pointCells->GetId(r), missingVals->GetId(m));
                      }
                    }
                  }
                  fprintf(stdout,"\n");
                }
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
// FixPatchesWithPolycubeOld
// ----------------------
int vtkSVGroupsSegmenter::FixPatchesWithPolycubeOld()
{
  // Then check everything
  // Extract surface, triangulate, and subdivide polycube
  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(this->PolycubePd);
  triangulator->Update();

  vtkNew(vtkPolyData, polycubePd);
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
  vtkNew(vtkTriangleFilter, triangulator);
  triangulator->SetInputData(this->PolycubePd);
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
    this->RotateGroupToGlobalAxis(this->PolycubePd, groupId, this->GroupIdsArrayName, rotPolycube, rotMatrix0, rotMatrix1);

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
        vtkErrorMacro("Could not find corresponding polycube point id");
        fprintf(stderr,"COULD NOT FIND: ");
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
  vtkNew(vtkPolyData, mappedPd);
  this->InterpolateMapOntoTarget(polycubePd, this->WorkPd, fullMapPd, mappedPd, this->GroupIdsArrayName);

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

  int w_div = this->PolycubeDivisions;
  if (w_div%2 == 0)
    w_div++;
  int h_div = this->PolycubeDivisions;
  if (h_div%2 == 0)
    h_div++;
  int l_div = 0; // Determined by length of cube

  std::vector<int> w_divs(numGroups);
  std::vector<int> h_divs(numGroups);
  std::vector<int> l_divs(numGroups);
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);

    // Extract surface of polycube
    vtkNew(vtkPolyData, branchPolycube);
    vtkSVGeneralUtils::ThresholdPd(this->PolycubePd, groupId, groupId, 1, this->GroupIdsArrayName, branchPolycube);

    branchPolycube->BuildLinks();

    fprintf(stdout,"FORMING PARA VOLUME FOR GROUP %d\n", groupId);
    vtkNew(vtkStructuredGrid, paraHexMesh);
    if (this->FormParametricHexMesh(branchPolycube, paraHexMesh, w_div,
                                    l_div, h_div) != SV_OK)
    {
      fprintf(stderr,"Couldn't do the dirt\n");
      return SV_ERROR;
    }
    fprintf(stdout,"WDIV: %d\n", w_div);
    fprintf(stdout,"HDIV: %d\n", h_div);
    fprintf(stdout,"LDIV: %d\n", l_div);
    w_divs[i] = w_div;
    h_divs[i] = h_div;
    l_divs[i] = l_div;
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

  this->PolycubeUg->DeepCopy(appender->GetOutput());

  vtkNew(vtkUnstructuredGrid, paraHexVolume);
  paraHexVolume->DeepCopy(appender->GetOutput());
  std::string fn = "/Users/adamupdegrove/Desktop/tmp/TEST_PARAHEXMESH.vtu";
  vtkSVIOUtils::WriteVTUFile(fn, paraHexVolume);

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
  cleaner2->ToleranceIsAbsoluteOn();
  cleaner2->SetAbsoluteTolerance(1.0e-6);
  cleaner2->SetInputData(paraHexVolume);
  cleaner2->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer2);
  surfacer2->SetInputData(cleaner2->GetOutput());
  surfacer2->Update();

  vtkNew(vtkPolyData, cleanSurface);
  cleanSurface->DeepCopy(surfacer2->GetOutput());

  this->RemoveInteriorCells(cleanSurface);

  std::vector<int> surfacePtMap;
  std::vector<std::vector<int> > invSurfacePtMap;
  this->GetInteriorPointMaps(paraHexSurface, paraHexCleanSurface, cleanSurface, surfacePtMap, invSurfacePtMap);

  fn = "/Users/adamupdegrove/Desktop/tmp/Mapping_All2.vtp";
  vtkSVIOUtils::WriteVTPFile(fn, fullMapPd);
  //std::string fn2 = "/Users/adamupdegrove/Desktop/tmp/ParaHexSurface.vtp";
  //vtkSVIOUtils::WriteVTPFile(fn2, paraHexSurface);
  vtkNew(vtkPolyData, mappedSurface);
  this->InterpolateMapOntoTarget(paraHexSurface, this->WorkPd, fullMapPd, mappedSurface, this->GroupIdsArrayName);
  //std::string filename5 = "/Users/adamupdegrove/Desktop/tmp/Mapped_Out2.vtp";
  //vtkSVIOUtils::WriteVTPFile(filename5, mappedSurface);

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

  this->FixInteriorBoundary(mappedSurface, invSurfacePtMap);

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
  cleaner3->ToleranceIsAbsoluteOn();
  cleaner3->SetAbsoluteTolerance(1.0e-6);
  cleaner3->SetInputData(volumeAppender->GetOutput());
  cleaner3->Update();

  vtkNew(vtkUnstructuredGrid, smoothVolume);
  smoothVolume->DeepCopy(cleaner3->GetOutput());

  std::vector<int> volumePtMap;
  std::vector<std::vector<int> > invVolumePtMap;
  this->GetVolumePointMaps(mappedVolume, smoothVolume, volumePtMap, invVolumePtMap);

  int smoothIters = 1500;
  if (this->SmoothUnstructuredGrid(smoothVolume, smoothIters, "empty") != SV_OK)
  {
    fprintf(stderr,"Couldn't smooth volume\n");
    return SV_ERROR;
  }

  //filename = "/Users/adamupdegrove/Desktop/tmp/TEST_SMOOTH.vtu";
  //vtkSVIOUtils::WriteVTUFile(filename, smoothVolume);

  this->FixVolume(mappedVolume, smoothVolume, volumePtMap);
  //this->SetControlMeshBoundaries(mappedVolume, smoothVolume, volumePtMap, invVolumePtMap);

  //std::string filename = "/Users/adamupdegrove/Desktop/tmp/TEST_FINAL.vtu";
  //vtkSVIOUtils::WriteVTUFile(filename, mappedVolume);
  this->FinalHexMesh->DeepCopy(mappedVolume);

  vtkNew(vtkAppendFilter, loftAppender);
  vtkNew(vtkSVNURBSCollection, nurbs);
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);

    vtkNew(vtkUnstructuredGrid, mappedBranch);
    vtkSVGeneralUtils::ThresholdUg(mappedVolume, groupId, groupId, 1, this->GroupIdsArrayName, mappedBranch);

    vtkNew(vtkStructuredGrid, realHexMesh);
    if (this->ConvertUGToSG(mappedBranch, realHexMesh, w_divs[i], l_divs[i], h_divs[i]) != SV_OK)
    {
      fprintf(stderr,"Couldn't do the dirt\n");
      return SV_ERROR;
    }

    //// FOR LOFTING OF VOLUME
    //// Set up the volume
    //vtkNew(vtkUnstructuredGrid, emptyGrid);
    //vtkNew(vtkSVLoftNURBSVolume, lofter);
    //lofter->SetInputData(emptyGrid);
    //lofter->SetInputGrid(realHexMesh);
    //lofter->SetUDegree(1);
    //lofter->SetVDegree(1);
    //lofter->SetWDegree(2);
    ////lofter->SetUnstructuredGridUSpacing(1./(10*w_divs[i]));
    ////lofter->SetUnstructuredGridVSpacing(1./(10*h_divs[i]));
    ////lofter->SetUnstructuredGridWSpacing(1./(10*l_divs[i]));
    //lofter->SetUnstructuredGridUSpacing(1./w_divs[i]);
    //lofter->SetUnstructuredGridVSpacing(1./h_divs[i]);
    //lofter->SetUnstructuredGridWSpacing(1./l_divs[i]);
    //lofter->SetUKnotSpanType("average");
    ////lofter->SetUKnotSpanType("derivative");
    //lofter->SetUParametricSpanType("chord");
    //lofter->SetVKnotSpanType("average");
    ////lofter->SetVKnotSpanType("derivative");
    //lofter->SetVParametricSpanType("chord");
    //lofter->SetWKnotSpanType("average");
    ////lofter->SetWKnotSpanType("derivative");
    //lofter->SetWParametricSpanType("chord");
    //lofter->Update();

//  //loftAppender->AddInputData(lofter->GetOutput());

    //nurbs->AddItem(lofter->GetVolume());

    //vtkNew(vtkStructuredGridGeometryFilter, converter);
    //converter->SetInputData(lofter->GetVolume()->GetControlPointGrid());
    //converter->Update();

    //std::string cpst = "/Users/adamupdegrove/Desktop/tmp/CONTROL_POINTS_STRUCT.vts";
    //vtkSVIOUtils::WriteVTSFile(cpst, lofter->GetVolume()->GetControlPointGrid());

    //std::string cps = "/Users/adamupdegrove/Desktop/tmp/CONTROL_POINTS.vtp";
    //vtkSVIOUtils::WriteVTPFile(cps, converter->GetOutput());

    // FOR USING HEX MESH AS CONTROL GRID
    int dim[3];
    realHexMesh->GetDimensions(dim);
    int nUCon = dim[0];
    int nVCon = dim[1];
    int nWCon = dim[2];
    int p = 1;
    int q = 1;
    int r = 2;
    std::string putype = "chord";
    std::string pvtype = "chord";
    std::string pwtype = "chord";
    std::string kutype = "average";
    std::string kvtype = "average";
    std::string kwtype = "average";

    // Set the temporary control points
    vtkNew(vtkPoints, tmpUPoints);
    tmpUPoints->SetNumberOfPoints(nUCon);
    for (int i=0; i<nUCon; i++)
    {
      int pos[3]; pos[0] = i; pos[1] = 0; pos[2] = 0;
      int ptId = vtkStructuredData::ComputePointId(dim, pos);
      tmpUPoints->SetPoint(i, realHexMesh->GetPoint(ptId));
    }

    // Get the input point set u representation
    vtkNew(vtkDoubleArray, U);
    if (vtkSVNURBSUtils::GetUs(tmpUPoints, putype, U) != SV_OK)
    {
      return SV_ERROR;
    }

    // Get the knots in the u direction
    vtkNew(vtkDoubleArray, uKnots);
    if (vtkSVNURBSUtils::GetKnots(U, p, kutype, uKnots) != SV_OK)
    {
      fprintf(stderr,"Error getting knots\n");
      return SV_ERROR;
    }
    //
    vtkNew(vtkPoints, tmpVPoints);
    tmpVPoints->SetNumberOfPoints(nVCon);
    for (int i=0; i<nVCon; i++)
    {
      int pos[3]; pos[0] = 0; pos[1] = i; pos[2] = 0;
      int ptId = vtkStructuredData::ComputePointId(dim, pos);
      tmpVPoints->SetPoint(i, realHexMesh->GetPoint(ptId));
    }
    // Get the input point set v representation
    vtkNew(vtkDoubleArray, V);
    if (vtkSVNURBSUtils::GetUs(tmpVPoints, pvtype, V) != SV_OK)
    {
      return SV_ERROR;
    }

    // Get the knots in the v direction
    vtkNew(vtkDoubleArray, vKnots);
    if (vtkSVNURBSUtils::GetKnots(V, q, kvtype, vKnots) != SV_OK)
    {
      fprintf(stderr,"Error getting knots\n");
      return SV_ERROR;
    }

    vtkNew(vtkPoints, tmpWPoints);
    tmpWPoints->SetNumberOfPoints(nWCon);
    for (int i=0; i<nWCon; i++)
    {
      int pos[3]; pos[0] = 0; pos[1] = 0; pos[2] = i;
      int ptId = vtkStructuredData::ComputePointId(dim, pos);
      tmpWPoints->SetPoint(i, realHexMesh->GetPoint(ptId));
    }
    // Get the input point set v representation
    vtkNew(vtkDoubleArray, W);
    if (vtkSVNURBSUtils::GetUs(tmpWPoints, pwtype, W) != SV_OK)
    {
      return SV_ERROR;
    }

    // Get the knots in the w direction
    vtkNew(vtkDoubleArray, wKnots);
    if (vtkSVNURBSUtils::GetKnots(W, r, kwtype, wKnots) != SV_OK)
    {
      fprintf(stderr,"Error getting knots\n");
      return SV_ERROR;
    }

    vtkNew(vtkSVNURBSVolume, hexMeshControlGrid);
    hexMeshControlGrid->SetKnotVector(uKnots, 0);
    hexMeshControlGrid->SetKnotVector(vKnots, 1);
    hexMeshControlGrid->SetKnotVector(wKnots, 2);
    hexMeshControlGrid->SetControlPoints(realHexMesh);
    hexMeshControlGrid->SetUDegree(p);
    hexMeshControlGrid->SetVDegree(q);
    hexMeshControlGrid->SetWDegree(r);

    nurbs->AddItem(hexMeshControlGrid);
  }
  vtkNew(vtkIdList, groupMap);
  for (int i=0; i<numGroups; i++)
  {
    int groupId = groupIds->GetId(i);
    groupMap->InsertNextId(groupId);
  }

  // Add patch connections for file writing
  for (int i=0; i<this->CenterlineGraph->NumberOfCells; i++)
  {
    vtkSVCenterlineGCell *gCell = this->CenterlineGraph->GetCell(i);

    int numChildren = gCell->Children.size();
    for (int j=0; j<numChildren; j++)
      nurbs->AddPatchConnection(i+1, groupMap->IsId(gCell->Children[j]->GroupId)+1, 1, 6);

    if (gCell->Parent != NULL)
    {
      int numBrothers = gCell->Parent->Children.size();
      for (int j=0; j<numBrothers; j++)
      {
        if (gCell->GroupId != gCell->Parent->Children[j]->GroupId)
          nurbs->AddPatchConnection(i+1, groupMap->IsId(gCell->Parent->Children[j]->GroupId)+1, 6, 6);
      }
    }
  }

  //if (this->MergedCenterlines->GetNumberOfCells() == 1)
  //{
  //  std::string mfsname = "/Users/adamupdegrove/Desktop/tmp/Pipe.msh";
  //  vtkNew(vtkSVMUPFESNURBSWriter, writer);
  //  writer->SetInputData(lofter->GetVolume());
  //  writer->SetFileName(mfsname.c_str());
  //  writer->Write();
  //}

  ////if (this->MergedCenterlines->GetNumberOfCells() == 1)
  ////{
  ////  std::string mfsname = "/Users/adamupdegrove/Desktop/tmp/Pipe.msh";
  ////  vtkNew(vtkSVMUPFESNURBSWriter, writer);
  ////  writer->SetInputData(lofter->GetVolume());
  ////  writer->SetFileName(mfsname.c_str());
  ////  writer->Write();
  ////}

  fprintf(stdout,"Writing NURBS...\n");
  std::string pername = "/Users/adamupdegrove/Desktop/tmp/perigee_nurbs.txt";
  vtkNew(vtkSVPERIGEENURBSCollectionWriter, writer);
  writer->SetInputData(nurbs);
  writer->SetFileName(pername.c_str());
  writer->Update();

  //loftAppender->Update();
  //loftedVolume->DeepCopy(loftAppender->GetOutput());

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

  int dim[3]; dim[0] = w_div; dim[1] = h_div; dim[2] = l_div;

  vtkNew(vtkPoints, sgPoints);
  sg->SetPoints(sgPoints);
  sg->GetPoints()->SetNumberOfPoints(dim[0]*dim[1]*dim[2]);
  sg->SetDimensions(dim);

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      for (int k=0; k<l_div; k++)
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
// CheckFace
// ----------------------
int vtkSVGroupsSegmenter::CheckFace(vtkPolyData *polycubePd, int faceId,
                                    int &nTopPts, int &nBotPts,
                                    int &flatTop, int &flatBot)
{
  double tol = 1.0e-4;
  fprintf(stdout,"CHECKING FACE: %d\n", faceId);
  vtkIdType npts, *ptIds;
  polycubePd->GetCellPoints(faceId, npts, ptIds);

  double (*pts)[3] = new double[npts][3];
  double (*vecs)[3] = new double[npts][3];

  for (int i=0; i<npts; i++)
  {
    int ptId0 = ptIds[i];
    int ptId1 = ptIds[(i+1)%npts];

    polycubePd->GetPoint(ptId0, pts[i]);
    polycubePd->GetPoint(ptId1, pts[(i+1)%npts]);

    //fprintf(stdout,"POINT %d: %.6f %.6f %.6f\n", i, pts[i][0], pts[i][1], pts[i][2]);
    //fprintf(stdout,"POINT %d: %.6f %.6f %.6f\n", (i+1)%npts, pts[(i+1)%npts][0], pts[(i+1)%npts][1], pts[(i+1)%npts][2]);
    //fprintf(stdout,"MAKES VEC: %d\n", i);

    vtkMath::Subtract(pts[(i+1)%npts], pts[i], vecs[i]);
    vtkMath::Normalize(vecs[i]);
  }

  flatTop = 0;
  nTopPts = 2;

  double testDot0 = vtkMath::Dot(vecs[0], vecs[1]);
  double testDot1 = vtkMath::Dot(vecs[0], vecs[2]);
  fprintf(stdout,"TEST DOT 0: %.8f\n", testDot0);
  fprintf(stdout,"TEST DOT 1: %.8f\n", testDot1);

  if (testDot0 < tol && testDot0 > -1.0*tol)
    flatTop = 1;
  else
  {
    if (!(fabs(testDot1) <= 1.0+tol && fabs(testDot1) >= 1.0-tol))
      nTopPts = 3;
  }

  if (testDot1 <= tol && testDot1 >= -1.0*tol)
  {
    flatTop = 1;
    nTopPts = 3;
  }

  flatBot = 0;
  nBotPts = 2;

  double testDot2 = vtkMath::Dot(vecs[npts-1], vecs[npts-2]);
  double testDot3 = vtkMath::Dot(vecs[npts-2], vecs[0]);
  fprintf(stdout,"TEST DOT 2: %.8f\n", testDot2);
  fprintf(stdout,"TEST DOT 3: %.8f\n", testDot3);

  if (testDot2 <= tol && testDot2 >= -1.0*tol)
    flatBot = 1;
  else
  {
    if (!(fabs(testDot3) <= 1.0+tol && fabs(testDot3) >= 1.0-tol))
      nBotPts=3;
  }

  if (fabs(testDot2) <= 1.0+tol && fabs(testDot2) >= 1.0-tol)
  {
    flatBot = 1;
    nBotPts = 3;
  }


  delete [] pts;
  delete [] vecs;

  return SV_OK;
}

// ----------------------
// FormParametricHexMesh
// ----------------------
int vtkSVGroupsSegmenter::FormParametricHexMesh(vtkPolyData *polycubePd, vtkStructuredGrid *paraHexMesh,
                                                int w_div, int &l_div, int h_div)
{
  int nTopPts0, nBotPts0, flatTop0, flatBot0;
  this->CheckFace(polycubePd, 0, nTopPts0, nBotPts0, flatTop0, flatBot0);
  fprintf(stdout,"FACE 0\n");
  fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts0);
  fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop0);
  fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts0);
  fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot0);

  int nTopPts1, nBotPts1, flatTop1, flatBot1;
  this->CheckFace(polycubePd, 1, nTopPts1, nBotPts1, flatTop1, flatBot1);
  fprintf(stdout,"FACE 1\n");
  fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts1);
  fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop1);
  fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts1);
  fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot1);

  int nTopPts2, nBotPts2, flatTop2, flatBot2;
  this->CheckFace(polycubePd, 2, nTopPts2, nBotPts2, flatTop2, flatBot2);
  fprintf(stdout,"FACE 2\n");
  fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts2);
  fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop2);
  fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts2);
  fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot2);

  int nTopPts3, nBotPts3, flatTop3, flatBot3;
  this->CheckFace(polycubePd, 3, nTopPts3, nBotPts3, flatTop3, flatBot3);
  fprintf(stdout,"FACE 3\n");
  fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts3);
  fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop3);
  fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts3);
  fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot3);

  int topHorzWedge = 0, topVertWedge = 0;
  int topSTet0 = 0, topSTet1 = 0, topSTet2 = 0, topSTet3 = 0;
  int topCTet0 = 0, topCTet1 = 0, topCTet2 = 0, topCTet3 = 0;

  if (nTopPts0 == 3)
  {
    if (flatTop0)
      topSTet2 = 1;
    else
    {
      if (nTopPts2 == 3)
      {
        if (flatTop2)
          topSTet0 = 1;
        else
          topHorzWedge = 1;
      }
      else
      {
        fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
        return SV_ERROR;
      }
    }
  }
  else if (nTopPts0 == 2)
  {
    if (!flatTop0)
    {
      if (!flatTop1)
        topCTet0 = 1;
      else if (!flatTop3)
        topCTet3 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 1\n");
        return SV_ERROR;
      }
    }
    if (!flatTop2)
    {
      if (!flatTop1)
        topCTet1 = 1;
      else if (!flatTop3)
        topCTet2 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 2\n");
        return SV_ERROR;
      }
    }
  }

  if (nTopPts3 == 3)
  {
    if (flatTop3)
      topSTet1 = 1;
    else
    {
      if (nTopPts1 == 3)
      {
        if (flatTop1)
          topSTet3 = 1;
        else
          topVertWedge = 1;
      }
      else
      {
        fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
        return SV_ERROR;
      }
    }
  }
  else if (nTopPts3 == 2)
  {
    if (!flatTop3)
    {
      if (!flatTop2)
        topCTet2 = 1;
      else if (!flatTop0)
        topCTet3 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 3\n");
        return SV_ERROR;
      }
    }
    if (!flatTop1)
    {
      if (!flatTop2)
        topCTet1 = 1;
      else if (!flatTop0)
        topCTet0 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 4\n");
        return SV_ERROR;
      }
    }
  }

  int botHorzWedge = 0, botVertWedge = 0;
  int botSTet0 = 0, botSTet1 = 0, botSTet2 = 0, botSTet3 = 0;
  int botCTet0 = 0, botCTet1 = 0, botCTet2 = 0, botCTet3 = 0;

  if (nBotPts0 == 3)
  {
    if (flatBot0)
      botSTet0 = 1;
    else
    {
      if (nBotPts2 == 3)
      {
        if (flatBot2)
          botSTet2 = 1;
        else
          botHorzWedge = 1;
      }
      else
      {
        fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
        return SV_ERROR;
      }
    }
  }
  else if (nBotPts0 == 2)
  {
    if (!flatBot0)
    {
      if (!flatBot1)
        botCTet0 = 1;
      else if (!flatBot3)
        botCTet3 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 5\n");
        return SV_ERROR;
      }
    }
    if (!flatBot2)
    {
      if (!flatBot1)
        botCTet1 = 1;
      else if (!flatBot3)
        botCTet2 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 6\n");
        return SV_ERROR;
      }
    }
  }

  if (nBotPts3 == 3)
  {
    if (flatBot3)
      botSTet3 = 1;
    else
    {
      if (nBotPts1 == 3)
      {
        if (flatBot1)
          botSTet1 = 1;
        else
          botVertWedge = 1;
      }
      else
      {
        fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
        return SV_ERROR;
      }
    }
  }
  else if (nBotPts3 == 2)
  {
    if (!flatBot3)
    {
      if (!flatBot2)
        botCTet2 = 1;
      else if (!flatBot0)
        botCTet3 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 7\n");
        return SV_ERROR;
      }
    }
    if (!flatBot1)
    {
      if (!flatBot2)
        botCTet1 = 1;
      else if (!flatBot0)
        botCTet0 = 1;
      else
      {
        fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 8\n");
        return SV_ERROR;
      }
    }
  }

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

  if (f0npts != f2npts || f1npts != f3npts)
  {
    fprintf(stderr,"Opposite sides of cube must have same number of points\n");
    return SV_ERROR;
  }

  // Form some what of a parallelepiped
  double f0Pts[4][3], f2Pts[4][3];
  polycubePd->GetPoint(f0PtIds[0], f0Pts[0]);
  polycubePd->GetPoint(f0PtIds[1], f0Pts[1]);
  polycubePd->GetPoint(f0PtIds[2], f0Pts[2]);
  polycubePd->GetPoint(f0PtIds[3], f0Pts[3]);
  polycubePd->GetPoint(f2PtIds[3], f2Pts[0]);
  polycubePd->GetPoint(f2PtIds[2], f2Pts[1]);
  polycubePd->GetPoint(f2PtIds[1], f2Pts[2]);
  polycubePd->GetPoint(f2PtIds[0], f2Pts[3]);

  double f1Pts[4][3], f3Pts[4][3];
  polycubePd->GetPoint(f3PtIds[0], f3Pts[0]);
  polycubePd->GetPoint(f3PtIds[1], f3Pts[1]);
  polycubePd->GetPoint(f3PtIds[2], f3Pts[2]);
  polycubePd->GetPoint(f3PtIds[3], f3Pts[3]);
  polycubePd->GetPoint(f1PtIds[3], f1Pts[0]);
  polycubePd->GetPoint(f1PtIds[2], f1Pts[1]);
  polycubePd->GetPoint(f1PtIds[1], f1Pts[2]);
  polycubePd->GetPoint(f1PtIds[0], f1Pts[3]);

  if (nTopPts0 == 3)
  {
    fprintf(stdout,"N TOP PTS 0: %d\n", nTopPts0);
    polycubePd->GetPoint(f0PtIds[3], f0Pts[2]);
    polycubePd->GetPoint(f0PtIds[4], f0Pts[3]);
  }

  if (nTopPts2 == 3)
  {
    fprintf(stdout,"N TOP PTS 2: %d\n", nTopPts2);
    polycubePd->GetPoint(f2PtIds[4], f2Pts[0]);
    polycubePd->GetPoint(f2PtIds[3], f2Pts[1]);
  }

  if (nTopPts1 == 3)
  {
    fprintf(stdout,"N TOP PTS 1: %d\n", nTopPts1);
    polycubePd->GetPoint(f1PtIds[4], f1Pts[0]);
    polycubePd->GetPoint(f1PtIds[3], f1Pts[1]);
  }

  if (nTopPts3 == 3)
  {
    fprintf(stdout,"N TOP PTS 3: %d\n", nTopPts3);
    polycubePd->GetPoint(f3PtIds[3], f3Pts[2]);
    polycubePd->GetPoint(f3PtIds[4], f3Pts[3]);
  }

  fprintf(stdout,"WHAT ARE POINTS\n");
  fprintf(stdout,"FACE 0: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f0Pts[0][0], f0Pts[0][1], f0Pts[0][2],
                                                                                             f0Pts[1][0], f0Pts[1][1], f0Pts[1][2],
                                                                                             f0Pts[2][0], f0Pts[2][1], f0Pts[2][2],
                                                                                             f0Pts[3][0], f0Pts[3][1], f0Pts[3][2]);
  fprintf(stdout,"FACE 1: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f1Pts[0][0], f1Pts[0][1], f1Pts[0][2],
                                                                                             f1Pts[1][0], f1Pts[1][1], f1Pts[1][2],
                                                                                             f1Pts[2][0], f1Pts[2][1], f1Pts[2][2],
                                                                                             f1Pts[3][0], f1Pts[3][1], f1Pts[3][2]);
  fprintf(stdout,"FACE 2: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f2Pts[0][0], f2Pts[0][1], f2Pts[0][2],
                                                                                             f2Pts[1][0], f2Pts[1][1], f2Pts[1][2],
                                                                                             f2Pts[2][0], f2Pts[2][1], f2Pts[2][2],
                                                                                             f2Pts[3][0], f2Pts[3][1], f2Pts[3][2]);
  fprintf(stdout,"FACE 3: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f3Pts[0][0], f3Pts[0][1], f3Pts[0][2],
                                                                                             f3Pts[1][0], f3Pts[1][1], f3Pts[1][2],
                                                                                             f3Pts[2][0], f3Pts[2][1], f3Pts[2][2],
                                                                                             f3Pts[3][0], f3Pts[3][1], f3Pts[3][2]);

  fprintf(stdout,"THIS SAYS TOP VERT WEDGE: %d\n", topVertWedge);
  fprintf(stdout,"THIS SAYS BOT VERT WEDGE: %d\n", botVertWedge);
  fprintf(stdout,"THIS SAYS TOP HORZ WEDGE: %d\n", topHorzWedge);
  fprintf(stdout,"THIS SAYS BOT HORZ WEDGE: %d\n", botHorzWedge);
  fprintf(stdout,"THIS SAYS TOP SIDE TET 0: %d\n", topSTet0);
  fprintf(stdout,"THIS SAYS TOP SIDE TET 1: %d\n", topSTet1);
  fprintf(stdout,"THIS SAYS TOP SIDE TET 2: %d\n", topSTet2);
  fprintf(stdout,"THIS SAYS TOP SIDE TET 3: %d\n", topSTet3);
  fprintf(stdout,"THIS SAYS TOP CORN TET 0: %d\n", topCTet0);
  fprintf(stdout,"THIS SAYS TOP CORN TET 1: %d\n", topCTet1);
  fprintf(stdout,"THIS SAYS TOP CORN TET 2: %d\n", topCTet2);
  fprintf(stdout,"THIS SAYS TOP CORN TET 3: %d\n", topCTet3);
  fprintf(stdout,"THIS SAYS BOT SIDE TET 0: %d\n", botSTet0);
  fprintf(stdout,"THIS SAYS BOT SIDE TET 1: %d\n", botSTet1);
  fprintf(stdout,"THIS SAYS BOT SIDE TET 2: %d\n", botSTet2);
  fprintf(stdout,"THIS SAYS BOT SIDE TET 3: %d\n", botSTet3);
  fprintf(stdout,"THIS SAYS BOT CORN TET 0: %d\n", botCTet0);
  fprintf(stdout,"THIS SAYS BOT CORN TET 1: %d\n", botCTet1);
  fprintf(stdout,"THIS SAYS BOT CORN TET 2: %d\n", botCTet2);
  fprintf(stdout,"THIS SAYS BOT CORN TET 3: %d\n", botCTet3);
  fprintf(stdout,"\n");

  // Sides of cube
  double face0Vec0[3], face0Vec1[3];
  vtkMath::Subtract(f0Pts[3], f0Pts[0], face0Vec0);
  vtkMath::Normalize(face0Vec0);
  double face0Dist0 = vtkSVMathUtils::Distance(f0Pts[3], f0Pts[0]);

  vtkMath::Subtract(f0Pts[2], f0Pts[1], face0Vec1);
  vtkMath::Normalize(face0Vec1);
  double face0Dist1 = vtkSVMathUtils::Distance(f0Pts[2], f0Pts[1]);

  double face2Vec0[3], face2Vec1[3];
  vtkMath::Subtract(f2Pts[3], f2Pts[0], face2Vec0);
  vtkMath::Normalize(face2Vec0);
  double face2Dist0 = vtkSVMathUtils::Distance(f2Pts[3], f2Pts[0]);

  vtkMath::Subtract(f2Pts[2], f2Pts[1], face2Vec1);
  vtkMath::Normalize(face2Vec1);
  double face2Dist1 = vtkSVMathUtils::Distance(f2Pts[2], f2Pts[1]);

  l_div = floor(vtkSVMathUtils::Distance(f0Pts[1], f0Pts[0])/(this->PolycubeUnitLength));
  if (l_div < 20)
    l_div = 20;

  if (topSTet1 || topSTet3 || botSTet1 || botSTet3)
    w_div = 2*w_div-1;
  if (topSTet0 || topSTet2 || botSTet0 || botSTet2)
    h_div = 2*h_div-1;

  // Top and bottom of cube
  double face5Vec0[3], face5Vec1[3];
  vtkMath::Subtract(f1Pts[0], f1Pts[3], face5Vec0);
  vtkMath::Normalize(face5Vec0);
  double face5Dist0 = vtkSVMathUtils::Distance(f1Pts[0], f1Pts[3]);

  vtkMath::Subtract(f3Pts[0], f3Pts[3], face5Vec1);
  vtkMath::Normalize(face5Vec1);
  double face5Dist1 = vtkSVMathUtils::Distance(f3Pts[0], f3Pts[3]);

  double face4Vec0[3], face4Vec1[3];
  vtkMath::Subtract(f1Pts[1], f1Pts[2], face4Vec0);
  vtkMath::Normalize(face4Vec0);
  double face4Dist0 = vtkSVMathUtils::Distance(f1Pts[1], f1Pts[2]);

  vtkMath::Subtract(f3Pts[1], f3Pts[2], face4Vec1);
  vtkMath::Normalize(face4Vec1);
  double face4Dist1 = vtkSVMathUtils::Distance(f3Pts[1], f3Pts[2]);

  vtkNew(vtkPoints, f5GridPts);
  f5GridPts->SetNumberOfPoints(w_div*h_div);
  vtkNew(vtkPoints, f4GridPts);
  f4GridPts->SetNumberOfPoints(w_div*h_div);

  int dim2D_1[3]; dim2D_1[0] = w_div; dim2D_1[1] = h_div; dim2D_1[2] = 1;

  for (int i=0; i<w_div; i++)
  {
    double x5Vec0[3], x5Vec1[3], x4Vec0[3], x4Vec1[3];

    for (int j=0; j<3; j++)
    {
      x5Vec0[j] = face5Vec0[j]*i*(face5Dist0/(w_div-1));
      x5Vec1[j] = face5Vec1[j]*i*(face5Dist1/(w_div-1));

      x4Vec0[j] = face4Vec0[j]*i*(face4Dist0/(w_div-1));
      x4Vec1[j] = face4Vec1[j]*i*(face4Dist1/(w_div-1));
    }

    double f5Start[3], f5End[3], f4Start[3], f4End[3];
    vtkMath::Add(f1Pts[3], x5Vec0, f5Start);
    vtkMath::Add(f3Pts[3], x5Vec1, f5End);

    vtkMath::Add(f1Pts[2], x4Vec0, f4Start);
    vtkMath::Add(f3Pts[2], x4Vec1, f4End);

    double f5HVec[3], f4HVec[3];
    vtkMath::Subtract(f5End, f5Start, f5HVec);
    vtkMath::Normalize(f5HVec);
    double f5HDist = vtkSVMathUtils::Distance(f5End, f5Start);

    vtkMath::Subtract(f4End, f4Start, f4HVec);
    vtkMath::Normalize(f4HVec);
    double f4HDist = vtkSVMathUtils::Distance(f4End, f4Start);

    // Another tet top face
    double face4DiagVecs[2][3], x4DiagVecs[2][3], x4DiagPts[2][3];
    double x4ToDiag[3], x4FromDiag[3], x4AcrossDiag[3];
    double face4DiagDists[2], x4ToDiagDist, x4FromDiagDist, x4AcrossDiagDist;
    if (topSTet1)
    {
      double midPt[3];
      polycubePd->GetPoint(f1PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(midPt, f1Pts[2], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(midPt, f1Pts[2]);

      vtkMath::Subtract(f1Pts[1], midPt, startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(midPt, f1Pts[1]);

      vtkMath::Subtract(midPt, f3Pts[2], face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f3Pts[2]);

      vtkMath::Subtract(f3Pts[1], midPt, face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f3Pts[1]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(h_div-1));
        startPtVecs[1][j] = startVecs[1][j]*((i+1)%h_div)*(startVecDists[1]/(h_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(h_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*((i+1)%h_div)*(face4DiagDists[1]/(h_div-1));
      }

      if (i <= h_div-1)
      {
        vtkMath::Add(f1Pts[2], startPtVecs[0], f4Start);
        vtkMath::Add(f3Pts[2], x4DiagVecs[0], x4DiagPts[0]);

        vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[0], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[0]);
      }
      else
      {
        vtkMath::Add(midPt, startPtVecs[1], f4Start);
        vtkMath::Add(midPt, x4DiagVecs[1], x4DiagPts[1]);

        vtkMath::Subtract(x4DiagPts[1], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);
      }
    }

    if (topSTet3)
    {
      double midPt[3];
      polycubePd->GetPoint(f3PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(midPt, f3Pts[2], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(midPt, f3Pts[2]);

      vtkMath::Subtract(f3Pts[1], midPt, startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(midPt, f3Pts[1]);

      vtkMath::Subtract(midPt, f1Pts[2], face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f1Pts[2]);

      vtkMath::Subtract(f1Pts[1], midPt, face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f1Pts[1]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(h_div-1));
        startPtVecs[1][j] = startVecs[1][j]*((i+1)%h_div)*(startVecDists[1]/(h_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(h_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*((i+1)%h_div)*(face4DiagDists[1]/(h_div-1));
      }


      if (i <= h_div-1)
      {
        vtkMath::Add(f3Pts[2], startPtVecs[0], f4End);
        vtkMath::Add(f1Pts[2], x4DiagVecs[0], x4DiagPts[0]);

        vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[0], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[0]);
      }
      else
      {
        vtkMath::Add(midPt, startPtVecs[1], f4End);
        vtkMath::Add(midPt, x4DiagVecs[1], x4DiagPts[1]);

        vtkMath::Subtract(x4DiagPts[1], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);
      }
    }

    if (topSTet2)
    {
      double midPt[3];
      polycubePd->GetPoint(f2PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(f0Pts[1], f2Pts[1], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(f0Pts[1], f2Pts[1]);

      vtkMath::Subtract(f0Pts[2], f2Pts[2], startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(f0Pts[2], f2Pts[2]);

      vtkMath::Subtract(f0Pts[1], midPt, face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f0Pts[1]);

      vtkMath::Subtract(f0Pts[2], midPt, face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f0Pts[2]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(w_div-1));
        startPtVecs[1][j] = startVecs[1][j]*i*(startVecDists[1]/(w_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(w_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*i*(face4DiagDists[1]/(w_div-1));
      }


      vtkMath::Add(f2Pts[1], startPtVecs[0], f4Start);
      vtkMath::Add(f2Pts[2], startPtVecs[1], f4End);
      vtkMath::Add(midPt, x4DiagVecs[0], x4DiagPts[0]);
      vtkMath::Add(midPt, x4DiagVecs[1], x4DiagPts[1]);

      vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
      vtkMath::Normalize(x4ToDiag);
      x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

      vtkMath::Subtract(x4DiagPts[1], x4DiagPts[0], x4AcrossDiag);
      vtkMath::Normalize(x4AcrossDiag);
      x4AcrossDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], x4DiagPts[0]);

      vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
      vtkMath::Normalize(x4FromDiag);
      x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);

    }

    if (topSTet0)
    {
      double midPt[3];
      polycubePd->GetPoint(f0PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(f0Pts[1], f2Pts[1], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(f0Pts[1], f2Pts[1]);

      vtkMath::Subtract(f0Pts[2], f2Pts[2], startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(f0Pts[2], f2Pts[2]);

      vtkMath::Subtract(midPt, f2Pts[1], face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f2Pts[1]);

      vtkMath::Subtract(midPt, f2Pts[2], face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f2Pts[2]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(w_div-1));
        startPtVecs[1][j] = startVecs[1][j]*i*(startVecDists[1]/(w_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(w_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*i*(face4DiagDists[1]/(w_div-1));
      }


      vtkMath::Add(f2Pts[1], startPtVecs[0], f4Start);
      vtkMath::Add(f2Pts[2], startPtVecs[1], f4End);
      vtkMath::Add(f2Pts[1], x4DiagVecs[0], x4DiagPts[0]);
      vtkMath::Add(f2Pts[2], x4DiagVecs[1], x4DiagPts[1]);

      vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
      vtkMath::Normalize(x4ToDiag);
      x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

      vtkMath::Subtract(x4DiagPts[1], x4DiagPts[0], x4AcrossDiag);
      vtkMath::Normalize(x4AcrossDiag);
      x4AcrossDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], x4DiagPts[0]);

      vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
      vtkMath::Normalize(x4FromDiag);
      x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);
    }

    // Tet top face
    double face4DiagVec[3], x4DiagVec[3], x4DiagPoint[3];
    double face4DiagDist;
    if (topCTet1 || topCTet3)
    {
      vtkMath::Subtract(f3Pts[1], f1Pts[2], face4DiagVec);
      vtkMath::Normalize(face4DiagVec);
      face4DiagDist = vtkSVMathUtils::Distance(f3Pts[1], f1Pts[2]);

      for (int j=0; j<3; j++)
        x4DiagVec[j] = face4DiagVec[j]*i*(face4DiagDist/(w_div-1));

      vtkMath::Add(f1Pts[2], x4DiagVec, x4DiagPoint);
    }
    if (topCTet0 || topCTet2)
    {
      vtkMath::Subtract(f1Pts[1], f3Pts[2], face4DiagVec);
      vtkMath::Normalize(face4DiagVec);
      face4DiagDist = vtkSVMathUtils::Distance(f1Pts[1], f3Pts[2]);

      for (int j=0; j<3; j++)
        x4DiagVec[j] = face4DiagVec[j]*i*(face4DiagDist/(w_div-1));

      vtkMath::Add(f3Pts[2], x4DiagVec, x4DiagPoint);
    }

    if (topCTet0 || topCTet1 || topCTet2 || topCTet3)
    {
      vtkMath::Subtract(x4DiagPoint, f4Start, x4ToDiag);
      vtkMath::Normalize(x4ToDiag);
      x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPoint, f4Start);

      vtkMath::Subtract(f4End, x4DiagPoint, x4FromDiag);
      vtkMath::Normalize(x4FromDiag);
      x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPoint);
    }

    // Tet bot face
    double face5DiagVec[3], x5DiagVec[3], x5DiagPoint[3], x5ToDiag[3], x5FromDiag[3];
    double face5DiagDist, x5ToDiagDist, x5FromDiagDist;

    if (botCTet1 || botCTet3)
    {
      vtkMath::Subtract(f3Pts[0], f1Pts[3], face5DiagVec);
      vtkMath::Normalize(face5DiagVec);
      face5DiagDist = vtkSVMathUtils::Distance(f3Pts[0], f1Pts[3]);

      for (int j=0; j<3; j++)
        x5DiagVec[j] = face5DiagVec[j]*i*(face5DiagDist/(w_div-1));

      vtkMath::Add(f1Pts[3], x5DiagVec, x5DiagPoint);
    }
    if (botCTet0 || botCTet2)
    {
      vtkMath::Subtract(f1Pts[0], f3Pts[3], face5DiagVec);
      vtkMath::Normalize(face5DiagVec);
      face5DiagDist = vtkSVMathUtils::Distance(f1Pts[0], f3Pts[3]);

      for (int j=0; j<3; j++)
        x5DiagVec[j] = face5DiagVec[j]*i*(face5DiagDist/(w_div-1));

      vtkMath::Add(f3Pts[3], x5DiagVec, x5DiagPoint);

    }

    if (botCTet0 || botCTet1 || botCTet2 || botCTet3)
    {
      vtkMath::Subtract(x5DiagPoint, f5Start, x5ToDiag);
      vtkMath::Normalize(x5ToDiag);
      x5ToDiagDist = vtkSVMathUtils::Distance(x5DiagPoint, f5Start);

      vtkMath::Subtract(f5End, x5DiagPoint, x5FromDiag);
      vtkMath::Normalize(x5FromDiag);
      x5FromDiagDist = vtkSVMathUtils::Distance(f5End, x5DiagPoint);
    }

    for (int j=0; j<h_div; j++)
    {
      double z5Vec[3], z4Vec[3];
      for (int k=0; k<3; k++)
      {
        z5Vec[k] = f5HVec[k]*j*(f5HDist/(h_div-1));
        z4Vec[k] = f4HVec[k]*j*(f4HDist/(h_div-1));
      }

      double new5Pt[3], new4Pt[3];
      vtkMath::Add(f5Start, z5Vec, new5Pt);
      vtkMath::Add(f4Start, z4Vec, new4Pt);

      if (topCTet1 || topCTet3)
      {
        if (j <= i)
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else
        {
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(j-i)*(x4FromDiagDist/((w_div-1)-i));
          vtkMath::Add(x4DiagPoint, z4Vec, new4Pt);
        }
      }
      if (topCTet0 || topCTet2)
      {
        if (j <= (w_div-1)-i)
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(j-((w_div-1)-i))*(x4FromDiagDist/diag_div);
          vtkMath::Add(x4DiagPoint, z4Vec, new4Pt);
        }
      }

      if (botCTet1 || botCTet3)
      {
        if (j <= i)
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z5Vec[k] = x5ToDiag[k]*j*(x5ToDiagDist/diag_div);
          vtkMath::Add(f5Start, z5Vec, new5Pt);
        }
        else
        {
          for (int k=0; k<3; k++)
            z5Vec[k] = x5FromDiag[k]*(j-i)*(x5FromDiagDist/((w_div-1)-i));
          vtkMath::Add(x5DiagPoint, z5Vec, new5Pt);
        }
      }
      if (botCTet0 || botCTet2)
      {
        if (j <= (w_div-1)-i)
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z5Vec[k] = x5ToDiag[k]*j*(x5ToDiagDist/diag_div);
          vtkMath::Add(f5Start, z5Vec, new5Pt);
        }
        else
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z5Vec[k] = x5FromDiag[k]*(j-((w_div-1)-i))*(x5FromDiagDist/diag_div);
          vtkMath::Add(x5DiagPoint, z5Vec, new5Pt);
        }
      }

      if (topSTet1)
      {
        if (i <= h_div-1)
        {
          if (j < h_div-i)
          {
            int diag_div = ((h_div-1)-i);
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-((h_div-1)-i))*(x4FromDiagDist/i);
            vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
          }
        }
        else
        {
          if (j < (i+1)%h_div)
          {
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/((i+1)%h_div));
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            int diag_div = ((h_div-1)-((i+1)%h_div));
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-((i+1)%h_div))*(x4FromDiagDist/diag_div);
            vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
          }
        }
      }

      if (topSTet3)
      {
        if (i <= h_div-1)
        {
          if (j <= i)
          {
            int diag_div = i;
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-i)*(x4FromDiagDist/((h_div-1)-i));
            vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
          }
        }
        else
        {
          if (j <= (h_div-1) - ((i+1)%h_div))
          {
            int diag_div = (h_div-1) - ((i+1)%h_div);
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/(diag_div));
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            int diag_div = (i+1)%h_div;
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-((h_div-1)-(i+1)%h_div))*(x4FromDiagDist/diag_div);
            vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
          }
        }
      }

      if (topSTet2)
      {
        if (j <= (w_div-1)-i)
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else if (j > (w_div-1)-i && j <= ((w_div-1)-i)+2*i)
        {
          int diag_div = 2*i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4AcrossDiag[k]*(j-((w_div-1)-i))*(x4AcrossDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
        }
        else
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(((j+1)%w_div)-i)*(x4FromDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
        }
      }

      if (topSTet0)
      {
        if (j <= i)
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else if (j > i && j <= (2*(w_div-1) - i))
        {
          int diag_div = (2*(w_div-1)-(2*i));
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4AcrossDiag[k]*(j-i)*(x4AcrossDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
        }
        else
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(j-(2*(w_div-1)-i))*(x4FromDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
        }
      }

      int pos[3]; pos[0]= i; pos[1] = j; pos[2] = 0;
      int pId = vtkStructuredData::ComputePointId(dim2D_1, pos);

      f5GridPts->SetPoint(pId, new5Pt);
      f4GridPts->SetPoint(pId, new4Pt);
    }
  }

  vtkNew(vtkPolyData, grid5Poly);
  grid5Poly->SetPoints(f5GridPts);
  vtkNew(vtkPolyData, grid4Poly);
  grid4Poly->SetPoints(f4GridPts);

  vtkNew(vtkPoints, paraPoints);
  int dim[3]; dim[0] = w_div; dim[1] = h_div; dim[2] = l_div;
  paraHexMesh->SetDimensions(dim);
  paraHexMesh->SetPoints(paraPoints);
  paraHexMesh->GetPoints()->SetNumberOfPoints(w_div*h_div*l_div);

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      int pos2D_1[3]; pos2D_1[0] = i; pos2D_1[1] = j; pos2D_1[2] = 0;
      int getId1 = vtkStructuredData::ComputePointId(dim2D_1, pos2D_1);

      double startPt1[3], endPt1[3];
      f5GridPts->GetPoint(getId1, startPt1);
      f4GridPts->GetPoint(getId1, endPt1);

      double hLVec[3];
      vtkMath::Subtract(startPt1, endPt1, hLVec);
      vtkMath::Normalize(hLVec);
      double hLDist = vtkSVMathUtils::Distance(endPt1, startPt1);

      for (int k=0; k<l_div; k++)
      {
        double yVec[3];
        for (int l=0; l<3; l++)
        {
          yVec[l] = hLVec[l]*k*(hLDist/(1 - l_div));
        }

        double newPt[3];
        vtkMath::Add(startPt1, yVec, newPt);

        int pos[3]; pos[0]= i; pos[1] = j; pos[2] = k;
        int pId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(pId, newPt);
      }
    }
  }

  double midPt0[3], midPt1[3];
  if (topHorzWedge)
  {
    polycubePd->GetPoint(f2PtIds[2], midPt0);
    polycubePd->GetPoint(f0PtIds[2], midPt1);
    this->PushStructuredGridXAxis(paraHexMesh, midPt0, midPt1, 0);
  }

  if (botHorzWedge)
  {
    if (topHorzWedge)
    {
    polycubePd->GetPoint(f2PtIds[5], midPt0);
    polycubePd->GetPoint(f0PtIds[5], midPt1);
    }
    else
    {
    polycubePd->GetPoint(f2PtIds[4], midPt0);
    polycubePd->GetPoint(f0PtIds[4], midPt1);
    }
    this->PushStructuredGridXAxis(paraHexMesh, midPt0, midPt1, 1);
  }

  if (topVertWedge)
  {
    polycubePd->GetPoint(f1PtIds[2], midPt0);
    polycubePd->GetPoint(f3PtIds[2], midPt1);
    this->PushStructuredGridZAxis(paraHexMesh, midPt0, midPt1, 0);
  }

  if (botVertWedge)
  {
    if (topVertWedge)
    {
      polycubePd->GetPoint(f1PtIds[5], midPt0);
      polycubePd->GetPoint(f3PtIds[5], midPt1);
    }
    else
    {
      polycubePd->GetPoint(f1PtIds[4], midPt0);
      polycubePd->GetPoint(f3PtIds[4], midPt1);
    }
    this->PushStructuredGridZAxis(paraHexMesh, midPt0, midPt1, 1);
  }

  return SV_OK;
}

// ----------------------
// PushStructuredGridZAxis
// ----------------------
int vtkSVGroupsSegmenter::PushStructuredGridZAxis(vtkStructuredGrid *paraHexMesh,
                                                  const double midPt0[3],
                                                  const double midPt1[3],
                                                  const int isBottom)
{
  int dim[3];
  paraHexMesh->GetDimensions(dim);

  int h_div = dim[1];
  double h = vtkSVMathUtils::Distance(midPt1, midPt0);
  double h_dist = h/(h_div-1);

  double h_vec[3];
  vtkMath::Subtract(midPt1, midPt0, h_vec);
  vtkMath::Normalize(h_vec);

  int w_div = dim[0];
  int half_w_div = floor(w_div/2.0);

  // Set new bottom extend points in the middle
  int pos[3];
  pos[0]= half_w_div;

  if (isBottom)
    pos[2] = 0;
  else
    pos[2] = dim[2]-1;

  for (int i=0; i<h_div; i++)
  {
    pos[1] = i;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);

    double z_vec[3];
    for (int j=0; j<3; j++)
      z_vec[j] = h_vec[j]*i*h_dist;

    double newPt[3];
    vtkMath::Add(midPt0, z_vec, newPt);

    // New point in middle
    paraHexMesh->GetPoints()->SetPoint(ptId, newPt);

    // New point on 0 side
    double firstPt[3];
    //vtkMath::Add(pt0, z_vec, firstPt);

    pos[0] = 0;
    int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(firstPtId, firstPt);
    //paraHexMesh->GetPoints()->SetPoint(firstPtId, firstPt);

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
    //vtkMath::Add(pt1, z_vec, lastPt);

    pos[0] = dim[0]-1;
    int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(lastPtId, lastPt);
    //paraHexMesh->GetPoints()->SetPoint(lastPtId, lastPt);

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

  int l_div = dim[2];

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      pos[0] = i; pos[1] = j;

      pos[2] = 0;
      int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

      double firstPt[3];
      paraHexMesh->GetPoint(firstPtId, firstPt);

      pos[2] = dim[2]-1;
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
        pos[2] = k;

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
                                                    const double midPt0[3],
                                                    const double midPt1[3],
                                                    const int isBottom)
{
  int dim[3];
  paraHexMesh->GetDimensions(dim);

  int w_div = dim[0];
  double w = vtkSVMathUtils::Distance(midPt1, midPt0);
  double w_dist = w/(w_div-1);

  double w_vec[3];
  vtkMath::Subtract(midPt1, midPt0, w_vec);
  vtkMath::Normalize(w_vec);

  int h_div = dim[1];
  int half_h_div = floor(h_div/2.0);

  // Set new bottom extend points in the middle
  int pos[3];
  pos[1]= half_h_div;

  if (isBottom)
    pos[2] = 0;
  else
    pos[2] = dim[2]-1;

  for (int i=0; i<w_div; i++)
  {
    pos[0] = i;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);

    double x_vec[3];
    for (int j=0; j<3; j++)
      x_vec[j] = w_vec[j]*i*w_dist;

    double newPt[3];
    vtkMath::Add(midPt0, x_vec, newPt);

    // New point in middle
    paraHexMesh->GetPoints()->SetPoint(ptId, newPt);

    // New point on 0 side
    double firstPt[3];
    //vtkMath::Add(pt0, x_vec, firstPt);

    pos[1] = 0;
    int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(firstPtId, firstPt);
    //paraHexMesh->GetPoints()->SetPoint(firstPtId, firstPt);

    double first_h = vtkSVMathUtils::Distance(newPt, firstPt);
    double first_h_dist = first_h/(half_h_div);

    double first_h_vec[3];
    vtkMath::Subtract(newPt, firstPt, first_h_vec);
    vtkMath::Normalize(first_h_vec);

    for (int j=0; j<half_h_div; j++)
    {
      pos[1] = j;

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
    //vtkMath::Add(pt1, x_vec, lastPt);

    pos[1] = dim[1]-1;
    int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(lastPtId, lastPt);
    //paraHexMesh->GetPoints()->SetPoint(lastPtId, lastPt);

    double last_h = vtkSVMathUtils::Distance(lastPt, newPt);
    double last_h_dist = last_h/(half_h_div);

    double last_h_vec[3];
    vtkMath::Subtract(lastPt, newPt, last_h_vec);
    vtkMath::Normalize(last_h_vec);

    for (int j=0; j<half_h_div; j++)
    {
      pos[1] = half_h_div+j;

      double z_vec[3];
      for (int k=0; k<3; k++)
        z_vec[k] = last_h_vec[k]*j*last_h_dist;

      double inPt[3];
      vtkMath::Add(newPt, z_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }
  }

  int l_div = dim[2];

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      pos[0] = i; pos[1] = j;

      pos[2] = 0;
      int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

      double firstPt[3];
      paraHexMesh->GetPoint(firstPtId, firstPt);

      pos[2] = dim[2]-1;
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
        pos[2] = k;

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
  ptMap.clear();
  ptMap.resize(numPoints);
  std::fill(ptMap.begin(), ptMap.end(), -1);

  int numCleanPoints = pdWithCleanInterior->GetNumberOfPoints();
  invPtMap.clear();
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
// GetVolumePointMaps
// ----------------------
int vtkSVGroupsSegmenter::GetVolumePointMaps(vtkUnstructuredGrid *ugAll,
                                             vtkUnstructuredGrid *ugClean,
                                             std::vector<int> &ptMap,
                                             std::vector<std::vector<int> > &invPtMap)
{
  vtkNew(vtkPointLocator, locator);
  locator->SetDataSet(ugClean);
  locator->BuildLocator();

  int numAllPoints = ugAll->GetNumberOfPoints();
  ptMap.clear();
  ptMap.resize(numAllPoints);

  int numCleanPoints = ugClean->GetNumberOfPoints();
  invPtMap.clear();
  invPtMap.resize(numCleanPoints);

  for (int i=0; i<numAllPoints; i++)
  {
    double pt[3];
    ugAll->GetPoint(i, pt);

    int ptId = locator->FindClosestPoint(pt);

    ptMap[i] = ptId;
    invPtMap[ptId].push_back(i);
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
// SetControlMeshBoundaries
// ----------------------
int vtkSVGroupsSegmenter::SetControlMeshBoundaries(vtkUnstructuredGrid *mappedVolume,
                                                   vtkUnstructuredGrid *cleanVolume,
                                                   const std::vector<int> ptMap,
                                                   const std::vector<std::vector<int> > invPtMap)
{
  int numPoints = mappedVolume->GetNumberOfPoints();

  vtkNew(vtkIntArray, isBoundaryPoint);
  isBoundaryPoint->SetNumberOfTuples(numPoints);
  isBoundaryPoint->FillComponent(0, -1);
  isBoundaryPoint->SetName("IsBoundaryPoint");

  std::vector<std::vector<int> > boundaryGroupMatchings;
  std::vector<std::vector<int> > groupSets;
  std::vector<int> pointGroupIds(numPoints, -1);
  for (int i=0; i<invPtMap.size(); i++)
  {
    if (invPtMap[i].size() > 1)
    {
      for (int j=0; j<invPtMap[i].size(); j++)
        isBoundaryPoint->SetTuple1(invPtMap[i][j], 1);
    }

    std::vector<int> groupIds;
    for (int j=0; j<invPtMap[i].size(); j++)
    {

      vtkNew(vtkIdList, pointCellIds);
      mappedVolume->GetPointCells(invPtMap[i][j], pointCellIds);

      if (pointCellIds->GetNumberOfIds() > 0)
      {
        int cellId = pointCellIds->GetId(0);
        int groupId = mappedVolume->GetCellData()->GetArray("GroupIds")->GetTuple1(cellId);
        groupIds.push_back(groupId);
        pointGroupIds[invPtMap[i][j]] = groupId;
      }
      else
      {
        fprintf(stderr,"All of these boundary points should be attached to at least one cell\n");
        return SV_ERROR;
      }
    }

    std::sort(groupIds.begin(), groupIds.end());
    groupSets.push_back(groupIds);

    int addBoundary = 1;
    for (int k=0; k<boundaryGroupMatchings.size(); k++)
    {
      if (boundaryGroupMatchings[k] == groupIds)
        addBoundary = 0;
    }

    if (addBoundary)
      boundaryGroupMatchings.push_back(groupIds);
  }

  //fprintf(stdout,"LET ME SEE THE MATCHINGS:\n");
  std::vector<std::vector<int> > pointsInMatching;
  std::vector<std::vector<int> > cleanPointsInMatching;
  for (int i=0; i<boundaryGroupMatchings.size(); i++)
  {
    std::vector<int> pointIds;
    std::vector<int> cleanPointIds;
    if (boundaryGroupMatchings[i].size() > 1)
    {
      //fprintf(stdout,"POINTS WITH CONNECT: ");
      //for (int j=0; j<boundaryGroupMatchings[i].size(); j++)
      //  fprintf(stdout,"%d ", boundaryGroupMatchings[i][j]);
      //fprintf(stdout,"\n");
      //fprintf(stdout,"       -> ");
      for (int j=0; j<invPtMap.size(); j++)
      {
        if (invPtMap[j].size() > 1)
        {
          if (groupSets[j] == boundaryGroupMatchings[i])
          {
            for (int k=0; k<invPtMap[j].size(); k++)
            {
              //fprintf(stdout,"%d ",invPtMap[j][k]);
              pointIds.push_back(invPtMap[j][k]);
            }
            cleanPointIds.push_back(j);
          }
        }
      }
      //fprintf(stdout,"\n");
    }
    pointsInMatching.push_back(pointIds);
    cleanPointsInMatching.push_back(cleanPointIds);
  }

  std::vector<std::vector<int> > ptEdgeNeighbors;
  this->GetPointConnectivity(cleanVolume, ptEdgeNeighbors);

  for (int i=0; i<boundaryGroupMatchings.size(); i++)
  {
    int numGroups = boundaryGroupMatchings[i].size();
    if (numGroups == 3)
    {
      // Set the interior ridge line first
      //fprintf(stdout,"POINTS WITH CONNECT: ");
      //for (int j=0; j<numGroups; j++)
      //  fprintf(stdout,"%d ", boundaryGroupMatchings[i][j]);
      //fprintf(stdout,"\n");
      std::vector<int> outsideIndices;
      for (int j=0; j<cleanPointsInMatching[i].size(); j++)
      {
        int cleanPointId = cleanPointsInMatching[i][j];
        int isInterior = cleanVolume->GetPointData()->GetArray("IsInteriorPoint")->GetTuple1(cleanPointId);
        //fprintf(stdout,"       -> %d IS INTERIOR: %d\n", cleanPointId, isInterior);
        if (!isInterior)
          outsideIndices.push_back(j);
      }
      //fprintf(stdout,"\n");
      if (outsideIndices.size() != 2)
      {
        fprintf(stdout,"There should be two points along interior boundary ridge, but there are %d\n", outsideIndices.size());
        return SV_ERROR;
      }

      int linePtId0 = cleanPointsInMatching[i][outsideIndices[0]];
      int linePtId1 = cleanPointsInMatching[i][outsideIndices[1]];
      double pt0[3], pt1[3];
      cleanVolume->GetPoint(linePtId0, pt0);
      cleanVolume->GetPoint(linePtId1, pt1);

      for (int j=0; j<cleanPointsInMatching[i].size(); j++)
      {
        if (j != outsideIndices[0] && j != outsideIndices[1])
        {
          double currPt[3];
          int cleanPointId = cleanPointsInMatching[i][j];
          cleanVolume->GetPoint(cleanPointId, currPt);

          double t;
          double closestPt[3];
          double dist = vtkLine::DistanceToLine(currPt, pt0, pt1, t, closestPt);

          cleanVolume->GetPoints()->SetPoint(cleanPointId, closestPt);
          for (int k=0; k<invPtMap[cleanPointId].size(); k++)
          {
            int pointId = invPtMap[cleanPointId][k];
            mappedVolume->GetPoints()->SetPoint(pointId, closestPt);
          }
        }
      }

      // Now set the plane points mister son guy
      double ridgeVec[3];
      vtkMath::Subtract(pt1, pt0, ridgeVec);
      vtkMath::Normalize(ridgeVec);

      for (int j=0; j<numGroups; j++)
      {
        std::vector<int> twoGroups(2);
        twoGroups[0] = boundaryGroupMatchings[i][j];
        twoGroups[1] = boundaryGroupMatchings[i][(j+1)%numGroups];
        std::sort(twoGroups.begin(), twoGroups.end());

        // Loop through again and find just these beeznees
        double avgPlaneNormal[3];
        avgPlaneNormal[0] = 0.0;
        avgPlaneNormal[1] = 0.0;
        avgPlaneNormal[2] = 0.0;
        int numPtsInPlane = 0;
        for (int k=0; k<boundaryGroupMatchings.size(); k++)
        {
          if (boundaryGroupMatchings[k] == twoGroups)
          {
            for (int l=0; l<cleanPointsInMatching[k].size(); l++)
            {
              double currPt[3];
              int cleanPointId = cleanPointsInMatching[k][l];
              cleanVolume->GetPoint(cleanPointId, currPt);

              double t;
              double closestPt[3];
              double dist = vtkLine::DistanceToLine(currPt, pt0, pt1, t, closestPt);

              double vec[3];
              vtkMath::Subtract(currPt, closestPt, vec);
              vtkMath::Normalize(vec);

              double planeNormal[3];
              vtkMath::Cross(ridgeVec, vec, planeNormal);
              vtkMath::Normalize(planeNormal);

              for (int m=0; m<3; m++)
                avgPlaneNormal[m] += planeNormal[m];
              numPtsInPlane++;
            }
          }
        }
        for (int k=0; k<3; k++)
          avgPlaneNormal[k] = avgPlaneNormal[k]/numPtsInPlane;

        for (int k=0; k<boundaryGroupMatchings.size(); k++)
        {
          if (boundaryGroupMatchings[k] == twoGroups)
          {
            for (int l=0; l<cleanPointsInMatching[k].size(); l++)
            {
              int cleanPointId = cleanPointsInMatching[k][l];

              std::vector<int> neighborIds;
              for (int m=0; m<ptEdgeNeighbors[cleanPointId].size(); m++)
              {
                int neighborPointId = ptEdgeNeighbors[cleanPointId][m];
                if (isBoundaryPoint->GetValue(neighborPointId) == -1)
                  neighborIds.push_back(neighborPointId);
              }

              if (neighborIds.size() != 2)
              {
                fprintf(stderr,"Should be two neighbors to this interace point, but there is %d\n", neighborIds.size());
                return SV_ERROR;
              }

              double neighborPt0[3], neighborPt1[3];
              cleanVolume->GetPoint(neighborIds[0], neighborPt0);
              cleanVolume->GetPoint(neighborIds[1], neighborPt1);

              double planeT;
              double planeClosestPt[3];
              vtkPlane::IntersectWithLine(neighborPt0, neighborPt1, avgPlaneNormal, pt0, planeT, planeClosestPt);

              cleanVolume->GetPoints()->SetPoint(cleanPointId, planeClosestPt);
              for (int k=0; k<invPtMap[cleanPointId].size(); k++)
              {
                int pointId = invPtMap[cleanPointId][k];
                mappedVolume->GetPoints()->SetPoint(pointId, planeClosestPt);
              }
            }
          }
        }
      }
    }
  }

  mappedVolume->GetPointData()->AddArray(isBoundaryPoint);

  return SV_OK;
}

// ----------------------
// GetPointConnectivity
// ----------------------
int vtkSVGroupsSegmenter::GetPointConnectivity(vtkUnstructuredGrid *hexMesh,
                                               std::vector<std::vector<int> > &ptEdgeNeighbors)
{
  hexMesh->BuildLinks();
  int numCells = hexMesh->GetNumberOfCells();
  int numPoints = hexMesh->GetNumberOfPoints();

  for (int i=0; i<numCells; i++)
  {
    if (hexMesh->GetCellType(i) != VTK_HEXAHEDRON)
    {
      vtkErrorMacro("All cells must be hexes");
      return SV_ERROR;
    }
  }

  ptEdgeNeighbors.clear();
  ptEdgeNeighbors.resize(numPoints);

  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, ptCellIds);
    hexMesh->GetPointCells(i, ptCellIds);

    vtkNew(vtkIdList, ptNeighbors);
    for (int j=0; j<ptCellIds->GetNumberOfIds(); j++)
    {
      vtkCell *cell = hexMesh->GetCell(ptCellIds->GetId(j));

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

    for (int j=0; j<ptNeighbors->GetNumberOfIds(); j++)
      ptEdgeNeighbors[i].push_back(ptNeighbors->GetId(j));
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
          vtkMath::Add(real_xpt, real_ypt, real_pt);
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
// RemoveInteriorCells
// ----------------------
int vtkSVGroupsSegmenter::RemoveInteriorCells(vtkPolyData *quadMesh)
{
  quadMesh->BuildLinks();
  int numCells = quadMesh->GetNumberOfCells();
  int numPoints = quadMesh->GetNumberOfPoints();

  for (int i=0; i<numCells; i++)
  {
    if (quadMesh->GetCellType(i) != VTK_QUAD)
    {
      vtkErrorMacro("All cells must be hexes");
      return SV_ERROR;
    }
  }

  vtkNew(vtkIdList, pointDeleteList);
  for (int i=0; i<numCells; i++)
  {
    vtkCell *cell = quadMesh->GetCell(i);

    int neighCount = 0;
    for (int l=0; l<4; l++)
    {
      vtkNew(vtkIdList, threePtIds);
      threePtIds->InsertNextId(cell->PointIds->GetId(l));
      threePtIds->InsertNextId(cell->PointIds->GetId((l+1)%4));
      threePtIds->InsertNextId(cell->PointIds->GetId((l+2)%4));

      vtkNew(vtkIdList, neighCellIds);
      quadMesh->GetCellNeighbors(i, threePtIds, neighCellIds);
      if (neighCellIds->GetNumberOfIds() != 0)
        neighCount++;
    }

    if (neighCount != 0)
      quadMesh->DeleteCell(i);
  }

  quadMesh->RemoveDeletedCells();
  quadMesh->BuildLinks();
  quadMesh->BuildCells();

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(quadMesh);
  cleaner->Update();

  quadMesh->DeepCopy(cleaner->GetOutput());

  quadMesh->BuildLinks();
  quadMesh->BuildCells();

  return SV_OK;
}

// ----------------------
// SmoothUnstructuredGrid
// ----------------------
int vtkSVGroupsSegmenter::SmoothUnstructuredGrid(vtkUnstructuredGrid *hexMesh,
                                                 const int iters,
                                                 std::string fixedPointsArrayName)
{
  hexMesh->BuildLinks();
  int numCells = hexMesh->GetNumberOfCells();
  int numPoints = hexMesh->GetNumberOfPoints();

  vtkNew(vtkIntArray, isInteriorPoint);
  isInteriorPoint->SetNumberOfTuples(numPoints);
  isInteriorPoint->SetName("IsInteriorPoint");

  vtkNew(vtkIntArray, isFixedPoint);
  if (vtkSVGeneralUtils::CheckArrayExists(hexMesh, 0, fixedPointsArrayName) == SV_OK)
    isFixedPoint = vtkIntArray::SafeDownCast(hexMesh->GetPointData()->GetArray(fixedPointsArrayName.c_str()));
  else
  {
    isFixedPoint->SetNumberOfTuples(numPoints);
    isFixedPoint->FillComponent(0, -1);
  }

  for (int i=0; i<numCells; i++)
  {
    if (hexMesh->GetCellType(i) != VTK_HEXAHEDRON)
    {
      vtkErrorMacro("All cells must be hexes");
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
          // Have to do this for special interior cells in which multiple boundaries
          // meeting as four points of one face may not actually correspond to
          // just one cell. Essentially, interior of anything > bifurcation.
          int neighCount = 0;
          for (int l=0; l<4; l++)
          {
            vtkNew(vtkIdList, threePtIds);
            threePtIds->InsertNextId(face->PointIds->GetId(l));
            threePtIds->InsertNextId(face->PointIds->GetId((l+1)%4));
            threePtIds->InsertNextId(face->PointIds->GetId((l+2)%4));

            vtkNew(vtkIdList, neighCellIds);
            hexMesh->GetCellNeighbors(ptCellIds->GetId(j), threePtIds, neighCellIds);
            if (neighCellIds->GetNumberOfIds() != 0)
              neighCount++;
          }
          if (neighCount == 0)
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

    isInteriorPoint->SetTuple1(i, interiorPoint);
    int fixedPoint = isFixedPoint->GetTuple1(i);
    if (interiorPoint && fixedPoint != 1)
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
      // If > 0 neighbors, that means this is interior son
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

  if (vtkSVGeneralUtils::CheckArrayExists(hexMesh, 0, "IsInteriorPoint"))
    hexMesh->GetPointData()->RemoveArray("IsInteriorPoint");

  hexMesh->GetPointData()->AddArray(isInteriorPoint);

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

  double zVec[3], tmpVec[3];
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
  double inputZVec[4], newZVec[4];
  inputZVec[0] = 0.0; inputZVec[1] = 0.0; inputZVec[2] = 0.0; inputZVec[3] = 0.0;
  newZVec[0] = 0.0; newZVec[1] = 0.0; newZVec[2] = 0.0; newZVec[3] = 0.0;
  for (int i=0; i<3; i++)
    inputZVec[i] = zVec[i];
  inputZVec[3] = 0.0;
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
                                                         vtkPolyData *mappedPd,
                                                         std::string dataMatchingArrayName)
{
  vtkNew(vtkSVMapInterpolator, interpolator);
  interpolator->SetInputData(0, sourceBasePd);
  interpolator->SetInputData(1, targetPd);
  interpolator->SetInputData(2, targetBasePd);
  interpolator->SetNumSourceSubdivisions(0);
  if (dataMatchingArrayName.c_str() != NULL)
  {
    interpolator->SetEnableDataMatching(1);
    interpolator->SetDataMatchingArrayName(dataMatchingArrayName.c_str());
  }
  interpolator->Update();

  mappedPd->DeepCopy(interpolator->GetOutput());

  return SV_OK;
}

// ----------------------
// MatchEndPatches
// ----------------------
int vtkSVGroupsSegmenter::MatchEndPatches(vtkPolyData *branchPd, vtkPolyData *polyBranchPd)
{
  // Get regions
  std::vector<Region> branchRegions;
  if (vtkSVGroupsSegmenter::GetRegions(branchPd, "PatchIds", branchRegions) != SV_OK)
  {
    fprintf(stderr,"Could not get regions on branch\n");
    return SV_ERROR;
  }

  int numPoints = branchPd->GetNumberOfPoints();
  for (int j=0; j<branchRegions.size(); j++)
  {
    if (branchRegions[j].CornerPoints.size() != 4)
    {
      fprintf(stderr,"Number of corners on region %d is %d, needs to be 4\n", j, branchRegions[j].CornerPoints.size());
      std::string badName = "/Users/adamupdegrove/Desktop/tmp/BADCORNERS.vtp";
      vtkSVIOUtils::WriteVTPFile(badName, branchPd);
      return SV_ERROR;
    }
  }

  // Get open boundary edges
  std::vector<int> openCornerPoints;
  std::vector<std::vector<int> > ccwOpenEdges;

  if (this->GetOpenBoundaryEdges(branchPd, branchRegions, "PatchIds",
                                 openCornerPoints, ccwOpenEdges) != SV_OK)
  {
    vtkErrorMacro("Error getting open boundary edges");
    return SV_ERROR;
  }

  //// Now get cw edges by just doing opposite
  //std::vector<std::vector<int> > cwOpenEdges(openCornerPoints.size());
  //for (int j=0; j<openCornerPoints.size()/4; j++)
  //{
  //  for (int k=0; k<4; k++)
  //  {
  //    int cwLoc = 4*j + k;
  //    int ccwLoc = 4*j + (k+3)%4;
  //    int ccwPointCount=ccwOpenEdges[ccwLoc].size()-1;
  //    for (int l=0; l<ccwOpenEdges[ccwLoc].size(); l++)
  //    {
  //      cwOpenEdges[cwLoc].push_back(ccwOpenEdges[ccwLoc][ccwPointCount]);
  //      ccwPointCount--;
  //    }
  //  }
  //}

  //fprintf(stdout,"JUST WANT TO SEE WHAT THESE EDGE THINGS LOOK LIKE\n");
  //for (int l=0; l<ccwOpenEdges.size(); l++)
  //{
  //  fprintf(stdout,"CCWEDGE ");
  //  for (int m=0; m<ccwOpenEdges[l].size(); m++)
  //    fprintf(stdout,"%d ", ccwOpenEdges[l][m]);
  //  fprintf(stdout,"\n");
  //}
  //for (int l=0; l<cwOpenEdges.size(); l++)
  //{
  //  fprintf(stdout,"CWEDGE ");
  //  for (int m=0; m<cwOpenEdges[l].size(); m++)
  //    fprintf(stdout,"%d ", cwOpenEdges[l][m]);
  //  fprintf(stdout,"\n");
  //}

  vtkDataArray *slicePointsArray = branchPd->GetPointData()->GetArray("SlicePoints");
  for (int j=0; j<ccwOpenEdges.size(); j++)
  {
    int edgeSize = ccwOpenEdges[j].size();
    int pointId0 = ccwOpenEdges[j][0];
    int pointId1 = ccwOpenEdges[j][1];

    fprintf(stdout,"LOOKING FOR POLYCUBE POINT MATCHING: %.1f\n", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(pointId0));
    fprintf(stdout,"LOOKING FOR POLYCUBE FOR BRANCH: %d\n", pointId0);

    vtkNew(vtkIdList, pointPatchValues);
    vtkSVGeneralUtils::GetPointCellsValues(branchPd, "PatchIds", pointId0, pointPatchValues);

    if (pointPatchValues->GetNumberOfIds() != 2)
    {
      fprintf(stderr,"Patch has been overlapped by another patch painting which means the patches are too far from the target point. Make sure polycube matches model topologically.\n");
      return SV_ERROR;
    }

    vtkNew(vtkIdList, newEdgeCell);
    branchPd->GetCellEdgeNeighbors(-1, pointId0, pointId1, newEdgeCell);

    if (newEdgeCell->GetNumberOfIds() != 1)
    {
      fprintf(stderr,"Failure obtaining edge cell\n");
      return SV_ERROR;
    }
    int ccwCellValue = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(newEdgeCell->GetId(0));

    fprintf(stdout,"LOOKING FOR: %d %d\n", pointPatchValues->GetId(0), pointPatchValues->GetId(1));
    fprintf(stdout,"AND CCW VALUE: %d\n", ccwCellValue);

    // Get polycube regions
    std::vector<Region> polyRegions;
    if (vtkSVGroupsSegmenter::GetRegions(polyBranchPd, "PatchIds", polyRegions) != SV_OK)
    {
      fprintf(stderr,"Could not get regions on branch\n");
      return SV_ERROR;
    }

    int matchPointId = -1;
    for (int k=0; k<polyRegions.size(); k++)
    {
      for (int l=0; l<polyRegions[k].BoundaryEdges.size(); l++)
      {
        int polyPtId0 = polyRegions[k].BoundaryEdges[l][0];
        int polyPtId1 = polyRegions[k].BoundaryEdges[l][1];

        vtkNew(vtkIdList, polyPatchValues);
        vtkSVGeneralUtils::GetPointCellsValues(polyBranchPd, "PatchIds", polyPtId0, polyPatchValues);

        for (int m=0; m<polyPatchValues->GetNumberOfIds(); m++)
        {
          int val = polyPatchValues->GetId(m);
          polyPatchValues->SetId(m, val%6);
        }

        vtkNew(vtkIdList, checkList0);
        checkList0->DeepCopy(polyPatchValues);

        checkList0->IntersectWith(pointPatchValues);

        if (checkList0->GetNumberOfIds() == pointPatchValues->GetNumberOfIds() &&
            checkList0->GetNumberOfIds() == polyPatchValues->GetNumberOfIds())
        {
          vtkNew(vtkIdList, polyEdgeCell);
          polyBranchPd->GetCellEdgeNeighbors(-1, polyPtId0, polyPtId1, polyEdgeCell);

          if (polyEdgeCell->GetNumberOfIds() == 1)
          {
            int ccwPolyCellValue = polyBranchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(polyEdgeCell->GetId(0));
            ccwPolyCellValue = ccwPolyCellValue%6;

            fprintf(stdout,"PC POINT: %d\n", polyPtId0);
            fprintf(stdout,"PC HAS: %d %d\n", polyPatchValues->GetId(0), polyPatchValues->GetId(1));
            fprintf(stdout,"AND CCW VALUE: %d\n", ccwPolyCellValue);

            if (ccwPolyCellValue == ccwCellValue)
            {
              fprintf(stdout,"WE FOUND POINT IN POLYCUBE THAT MATCHES THIS POINT\n");
              matchPointId = polyBranchPd->GetPointData()->GetArray("SlicePoints")->GetTuple1(polyPtId0);
            }
          }
        }
      }
    }

    if (matchPointId == -1)
    {
      fprintf(stderr,"DIDNT FIND POINT TO MATCH TO!!\n");
      return SV_ERROR;
    }

    int branchMatchId = branchPd->GetPointData()->GetArray("TmpInternalIds")->LookupValue(matchPointId);
    if (branchMatchId == -1)
    {
      fprintf(stderr,"MATCHING POINT NOT ON BRANCH!!\n");
      return SV_ERROR;
    }
    fprintf(stdout,"THE POINT ON FULFUL PD TO GO TO IS: %d\n", matchPointId);
    fprintf(stdout,"THE POINT ON BRANCH PD TO GO TO IS: %d\n", branchMatchId);

    //Get the full ccw and cw edges
    std::vector<int> ccwFullEdges;
    std::vector<int> cwFullEdges;

    int edgeId = j;
    for (int k=0; k<4; k++)
    {
      for (int l=0; l<ccwOpenEdges[edgeId].size()-1; l++)
        ccwFullEdges.push_back(ccwOpenEdges[edgeId][l]);

      if (j < 4)
        edgeId = (edgeId+1)%4;
      else
        edgeId = (edgeId+1)%4 + 4;
    }

    int fullEdgeSize = ccwFullEdges.size();
    cwFullEdges.push_back(ccwFullEdges[0]);
    for (int k=0; k<ccwFullEdges.size()-1; k++)
      cwFullEdges.push_back(ccwFullEdges[fullEdgeSize-k-1]);

    int ccwCount = 0;
    for (int k=0; k<ccwFullEdges.size(); k++)
    {
      //if (slicePointsArray->GetTuple1(ccwFullEdges[k]) == 1)
      if (ccwFullEdges[k] == branchMatchId)
        break;
      ccwCount++;
    }
    int cwCount = 0;
    for (int k=0; k<cwFullEdges.size(); k++)
    {
      //if (slicePointsArray->GetTuple1(cwFullEdges[k]) == 1)
      if (cwFullEdges[k] == branchMatchId)
        break;
      cwCount++;
    }

    fprintf(stdout,"POINT IS:                %d\n", pointId0);
    fprintf(stdout,"COUNTER CLOCKWISE COUNT: %d\n", ccwCount);
    fprintf(stdout,"FULL COUNT:              %d\n", ccwFullEdges.size());
    fprintf(stdout,"CLOCKWISE COUNT:         %d\n", cwCount);
    fprintf(stdout,"FULL COUNT:              %d\n", cwFullEdges.size());

    if (ccwCount == ccwFullEdges.size() && cwCount == cwFullEdges.size())
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
      slicePoint = ccwFullEdges[stopCount];
      startPoint = ccwFullEdges[0];
      secondPoint = ccwFullEdges[1];
      opposPoint = cwFullEdges[1];
    }
    else
    {
      stopCount = cwCount;
      slicePoint = cwFullEdges[stopCount];
      startPoint = cwFullEdges[0];
      secondPoint = cwFullEdges[1];
      opposPoint = ccwFullEdges[1];
    }
    fprintf(stdout,"SLICE POINT IS:                %.1f\n", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(slicePoint));
    fprintf(stdout,"START POINT IS:                %.1f\n", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(startPoint));
    fprintf(stdout,"SECON POINT IS:                %.1f\n", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(secondPoint));
    fprintf(stdout,"OPPOS POINT IS:                %.1f\n", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(opposPoint));

    vtkNew(vtkIdList, tmpCell);
    branchPd->GetCellEdgeNeighbors(-1, startPoint, secondPoint, tmpCell);
    int edgeCell = tmpCell->GetId(0);
    int matchCellValue = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(edgeCell);

    branchPd->GetCellEdgeNeighbors(-1, startPoint, opposPoint, tmpCell);
    int opposCell = tmpCell->GetId(0);
    int newCellValue = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(opposCell);

    int patchPoint = -1;
    fprintf(stdout,"STOP COUNT IS!!!: %d\n", stopCount);
    if (stopCount != 0)
    {
      int count = 1;
      std::vector<int> tempNodes;
      tempNodes.clear();
      tempNodes.push_back(startPoint);

      for (int k=0; k<count; k++)
      {
        vtkNew(vtkIdList, pointCells);
        branchPd->GetPointCells(tempNodes[k], pointCells);
        for (int l=0; l<pointCells->GetNumberOfIds(); l++)
        {
          int cellId = pointCells->GetId(l);
          int nextPoint;
          if (ccwCount < cwCount)
            nextPoint = vtkSVGroupsSegmenter::GetCWPoint(branchPd, tempNodes[k], cellId);
          else
            nextPoint = vtkSVGroupsSegmenter::GetCCWPoint(branchPd, tempNodes[k], cellId);
          int isGoodEdge = vtkSVGroupsSegmenter::CheckCellValuesEdge(branchPd, "PatchIds", cellId, tempNodes[k], nextPoint);

          int cellValue = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(cellId);

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
      if ((count < stopCount && count != -1) || patchPoint == -1)
      {
        fprintf(stdout,"WHAT IS THE COUNT %d\n", count);
        stopCount = count/2;
        count = 1;
        std::vector<int> tempNodes;
        tempNodes.clear();
        tempNodes.push_back(startPoint);

        for (int k=0; k<count; k++)
        {
          vtkNew(vtkIdList, pointCells);
          branchPd->GetPointCells(tempNodes[k], pointCells);
          for (int l=0; l<pointCells->GetNumberOfIds(); l++)
          {
            int cellId = pointCells->GetId(l);
            int nextPoint;
            if (ccwCount < cwCount)
              nextPoint = vtkSVGroupsSegmenter::GetCWPoint(branchPd, tempNodes[k], cellId);
            else
              nextPoint = vtkSVGroupsSegmenter::GetCCWPoint(branchPd, tempNodes[k], cellId);
            int isGoodEdge = vtkSVGroupsSegmenter::CheckCellValuesEdge(branchPd, "PatchIds", cellId, tempNodes[k], nextPoint);

            int cellValue = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(cellId);

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
    }
    else
      patchPoint = startPoint;

    if (patchPoint == -1)
    {
      fprintf(stdout,"DIDNT COMPUTE PATCH POINT\n");
      std::string badName = "/Users/adamupdegrove/Desktop/tmp/BADPATCHPOINT.vtp";
      vtkSVIOUtils::WriteVTPFile(badName, branchPd);
      return SV_ERROR;
    }

    fprintf(stdout,"PATCH POINT:              %.1f\n", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(patchPoint));
    fprintf(stdout,"SLICE POINT:              %.1f\n", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(slicePoint));
    fprintf(stdout,"FILLING WITH:             %d\n", newCellValue);
    if (stopCount != 0)
    {
      vtkNew(vtkSVFindGeodesicPath, finder);
      finder->SetInputData(branchPd);
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
        fprintf(stdout,"%.1f ", branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(tmpIds->GetId(l)));
      fprintf(stdout,"\n");

      int count = 1;
      std::vector<int> tempCells;
      tempCells.push_back(edgeCell);

      for (int k=0; k<count; k++)
      {
        int tmpCellId = tempCells[k];
        branchPd->GetCellData()->GetArray("PatchIds")->SetTuple1(tmpCellId, newCellValue);
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
            int cellValue = branchPd->GetCellData()->GetArray("PatchIds")->GetTuple1(testCellId);
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
// GetOpenBoundaryEdges
// ----------------------
int vtkSVGroupsSegmenter::GetOpenBoundaryEdges(vtkPolyData *pd,
                                               std::vector<int> &openCornerPoints,
                                               std::vector<std::vector<int> > &openEdges)
{
  int numPoints = pd->GetNumberOfPoints();
  int numCells  = pd->GetNumberOfCells();
  std::vector<int> tempCornerPoints;  // In ccw order
  std::vector<int> pointUsed(numPoints, 0);

  for (int i=0; i<numCells; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);

    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];

      if (!pointUsed[ptId0])
      {
        vtkNew(vtkIdList, cellEdgeNeighbor);
        pd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellEdgeNeighbor);

        if (cellEdgeNeighbor->GetNumberOfIds() == 0)
        {
          pointUsed[ptId0] = 1;
          int startCornerPt = ptId0;

          int count=1;
          std::vector<int> tempNodes;
          tempNodes.push_back(startCornerPt);
          openCornerPoints.push_back(startCornerPt);

          for (int k=0; k<count; k++)
          {
            vtkNew(vtkIdList, pointCells);
            pd->GetPointCells(tempNodes[k], pointCells);
            for (int l=0; l<pointCells->GetNumberOfIds(); l++)
            {
              int cellId = pointCells->GetId(l);
              int pointCCWId = vtkSVGroupsSegmenter::GetCCWPoint(pd, tempNodes[k], cellId);

              // Check if open edge
              vtkNew(vtkIdList, edgeCells);
              pd->GetCellEdgeNeighbors(cellId, tempNodes[k], pointCCWId, edgeCells);

              if (edgeCells->GetNumberOfIds() == 0 && pointCCWId != startCornerPt)
              {
                tempNodes.push_back(pointCCWId);
                pointUsed[pointCCWId] = 1;
                count++;
              }
              else if (edgeCells->GetNumberOfIds() == 0 && pointCCWId == startCornerPt)
              {
                tempNodes.push_back(pointCCWId);
                openEdges.push_back(tempNodes);

                tempNodes.clear();

                count = -1;
                break;
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
// GetOpenBoundaryEdges
// ----------------------
int vtkSVGroupsSegmenter::GetOpenBoundaryEdges(vtkPolyData *pd,
                                               std::vector<Region> regions,
                                               std::string arrayName,
                                               std::vector<int> &openCornerPoints,
                                               std::vector<std::vector<int> > &openEdges)
{
  int numPoints = pd->GetNumberOfPoints();
  std::vector<int> tempCornerPoints;  // In ccw order
  std::vector<int> pointUsed(numPoints, 0);
  for (int j=0; j<regions.size(); j++)
  {
    for (int k=0; k<regions[j].CornerPoints.size(); k++)
    {
      int cornerPtId = regions[j].CornerPoints[k];
      vtkNew(vtkIdList, pointPatchValues);
      vtkSVGeneralUtils::GetPointCellsValues(pd, arrayName, cornerPtId, pointPatchValues);

      if (pointPatchValues->GetNumberOfIds() == 2 && pointUsed[cornerPtId] == 0)
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

  int count=1;
  std::vector<int> tempNodes;
  tempNodes.push_back(startCornerPt);
  openCornerPoints.push_back(startCornerPt);

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
          openEdges.push_back(tempNodes);

          tempNodes.clear();

          if (openCornerPoints.size() == tempCornerPoints.size())
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

              for (int jj=0; jj<openCornerPoints.size(); jj++)
              {
                if (tempIndex == openCornerPoints[jj])
                  tempCount = true;
              }
              if (tempCount == false)
              {
                startCornerPt = tempIndex;
                break;
              }
            }

            openCornerPoints.push_back(startCornerPt);
            tempNodes.push_back(startCornerPt);
            count = 1;
            j = -1;
            break;
          }
        }
        else
        {
          tempNodes.push_back(pointCCWId);
          openCornerPoints.push_back(pointCCWId);
          openEdges.push_back(tempNodes);
          tempNodes.clear();
          tempNodes.push_back(pointCCWId);
          count = 1;
          j = -1;
          break;
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
            fprintf(stdout,"DELETING CELLS: %d %d\n", i, cellIds->GetId(k));
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
    fprintf(stdout,"RE-LINKING\n");
    polycubePd->RemoveDeletedCells();

    vtkNew(vtkCleanPolyData, cleaner);
    cleaner->SetInputData(polycubePd);
    cleaner->ToleranceIsAbsoluteOn();
    cleaner->SetAbsoluteTolerance(1.0e-6);
    cleaner->Update();

    polycubePd->DeepCopy(cleaner->GetOutput());;
    polycubePd->BuildLinks();
  }

  fprintf(stdout,"GETTING SURFACE GROUPS\n");
  std::vector<Region> surfaceGroups;
  if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, surfaceGroups) != SV_OK)
  {
    vtkErrorMacro("Couldn't get groups");
    return SV_ERROR;
  }

  fprintf(stdout,"GETTING POLYCUBE GROUPS\n");
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
    fprintf(stdout,"NUMBER OF SURFACE GROUPS: %d\n", numSurfaceGroups);
    fprintf(stdout,"NUMBER OF POLYCUBE GROUPS: %d\n", numPolycubeGroups);

    if (numSurfaceGroups > numPolycubeGroups)
    {
      fprintf(stdout,"ADDITIONAL SURFACE GROUPS, SEE IF WE CAN REDUCE\n");
      if (this->FixMultipleGroups(this->WorkPd, polycubePd, surfaceGroups, polycubeGroups) != SV_OK)
      {
        return SV_ERROR;
      }

      fprintf(stdout,"RE-GETTING SURFACE GROUPS\n");
      if (this->GetRegions(this->WorkPd, this->GroupIdsArrayName, surfaceGroups) != SV_OK)
      {
        vtkErrorMacro("Couldn't get groups");
        return SV_ERROR;
      }

      fprintf(stdout,"RE-GETTING POLYCUBE GROUPS\n");
      if (this->GetRegions(polycubePd, this->GroupIdsArrayName, polycubeGroups) != SV_OK)
      {
        vtkErrorMacro("Couldn't get groups");
        return SV_ERROR;
      }
    }
    else
    {
      fprintf(stdout,"NOT ENOUGH SURFACE REGIONS TO MATCH POLYCUBE\n");
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
        fprintf(stdout,"NUMRBS OF EDGES FROM SURFACE GROUP  %d is %d\n", surfaceGroups[i].IndexCluster, surfaceGroups[i].BoundaryEdges.size());
        fprintf(stdout,"NUMRBS OF EDGES FROM POLYCUBE GROUP %d is %d\n", polycubeGroups[j].IndexCluster, polycubeGroups[j].BoundaryEdges.size());
        for (int k=0; k<surfaceGroups[i].CornerPoints.size(); k++)
        {
          int cornerPtId = surfaceGroups[i].CornerPoints[k];

          vtkNew(vtkIdList, surfaceCellList);
          vtkSVGeneralUtils::GetPointCellsValues(origPd, this->GroupIdsArrayName, cornerPtId, surfaceCellList);

          fprintf(stdout,"SURFACE CORNER POINT %d GROUPS ARE ", k);
          for (int l=0; l<surfaceCellList->GetNumberOfIds(); l++)
            fprintf(stdout,"%d ", surfaceCellList->GetId(l));
          fprintf(stdout,"\n");
        }
        for (int k=0; k<polycubeGroups[j].CornerPoints.size(); k++)
        {
          int cornerPtId = polycubeGroups[j].CornerPoints[k];

          vtkNew(vtkIdList, polycubeCellList);
          vtkSVGeneralUtils::GetPointCellsValues(polycubePd, this->GroupIdsArrayName, cornerPtId, polycubeCellList);

          fprintf(stdout,"POLYCUBE CORNER POINT %d GROUPS ARE ", k);
          for (int l=0; l<polycubeCellList->GetNumberOfIds(); l++)
            fprintf(stdout,"%d ", polycubeCellList->GetId(l));
          fprintf(stdout,"\n");
        }

        std::vector<std::vector<int> > surfConnectedCornerPts;
        this->GetConnectedEdges(surfaceGroups[i].BoundaryEdges, surfConnectedCornerPts);

        fprintf(stdout,"NUMBER OF CONNECTED EDGES: %d\n", surfConnectedCornerPts.size());
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
                fprintf(stdout,"WE DID IT!, WE FOUND A MATCH!\n");
                foundMatch = 1;
              }
           }

            if (foundMatch == 0)
            {
              fprintf(stdout,"UH OH, DIDNT FIND MATCH!!!\n");
              badEdges.push_back(edgeCount);
            }
            allEdges.push_back(edgeCount);
            edgeCount++;
          }
          fprintf(stdout,"NUMBER BAD EDGES: %d\n", badEdges.size());
          fprintf(stdout,"NUMBER ALL EDGES: %d\n", allEdges.size());

          if (surfaceGroups[i].BoundaryEdges.size() == polycubeGroups[j].BoundaryEdges.size() && badEdges.size() == 0)
          {
            fprintf(stdout,"WE GUCCI\n");
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
            fprintf(stdout,"NO FIX FOR THIS HAS BEEN DEVISED YET!!!!\n");
        }
      }
    }
    fprintf(stdout,"\n");
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

  fprintf(stdout,"TOTAL NUM OF CELLS: %d\n", this->WorkPd->GetNumberOfCells());
  for (int i=0; i<critPts->GetNumberOfIds(); i++)
  {
    int halfId = critPts->GetId(i);

    vtkNew(vtkIdList, finalVals);
    vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, halfId, finalVals);
    if (finalVals->GetNumberOfIds() != 4)
    {
        fprintf(stdout,"NO GOOD, FIX GROUPS AROUND POINT: %d\n", halfId);
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
            int pointCCWId = vtkSVGroupsSegmenter::GetCCWPoint(this->WorkPd, tempNodes[j], cellId);

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
          fprintf(stdout,"\n");
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
        fprintf(stdout,"\n");
      }
    }
  }

  return SV_OK;
}

// ----------------------
// FixMultipleGroups
// ----------------------
int vtkSVGroupsSegmenter::FixMultipleGroups(vtkPolyData *pd, vtkPolyData *polycubePd, std::vector<Region> surfaceGroups, std::vector<Region> polycubeGroups)
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
              fprintf(stdout,"FOUNDING MATCHING EDGE GROUPS FOR GROUP %d AND %d WITH THE FOLLOWING: ", surfaceGroups[groupLocs[0]].IndexCluster,
                                                                                                       surfaceGroups[groupLocs[1]].IndexCluster);
              for (int l=0; l<intersectList->GetNumberOfIds(); l++)
                fprintf(stdout," %d", intersectList->GetId(l));
              fprintf(stdout,"\n");

              for (int m=0; m<allPolycubeGroups.size(); m++)
              {
                vtkNew(vtkIdList, polycubeGroupsList);
                for (int n=0; n<allPolycubeGroups[m].size(); n++)
                  polycubeGroupsList->InsertUniqueId(allPolycubeGroups[m][n]);

                vtkNew(vtkIdList, newIntersectList);
                newIntersectList->DeepCopy(intersectList);

                newIntersectList->IntersectWith(polycubeGroupsList);

                fprintf(stdout,"INTERSECTED: ");
                for (int o=0; o<newIntersectList->GetNumberOfIds(); o++)
                  fprintf(stdout,"%d ",newIntersectList->GetId(o));
                fprintf(stdout,"\n");

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
        fprintf(stdout,"THREE PATCHES OF ONE GROUP, CANNOT HANDLE THIS\n");
        return SV_ERROR;
      }
    }
  }

  return SV_OK;
}

// ----------------------
// GetConnectedEdges
// ----------------------
int vtkSVGroupsSegmenter::GetConnectedEdges(std::vector<std::vector<int> > inputEdges, std::vector<std::vector<int> > &connectedCornerPts)
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
int vtkSVGroupsSegmenter::FixCloseGroup(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd, std::string arrayName,
                                        const Region region, const Region polyRegion, std::vector<int> allEdges,
                                        std::vector<int> badEdges, vtkIdList *critPts)
{
  fprintf(stdout,"--FIX CLOSE GROUP--\n");
  int patchValue = region.IndexCluster;
  fprintf(stdout,"CLUSTER %d\n", patchValue);
  int numEdges = region.BoundaryEdges.size();
  fprintf(stdout,"NUM EDGES %d\n", numEdges);

  fprintf(stdout,"NUMBER OF ALL EDGES: %d\n", allEdges.size());
  fprintf(stdout,"ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    fprintf(stdout,"%d ", allEdges[j]);
  fprintf(stdout,"\n");
  fprintf(stdout,"NUMBER OF BAD EDGES: %d\n", badEdges.size());
  fprintf(stdout,"BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    fprintf(stdout,"%d ", badEdges[j]);
  fprintf(stdout,"\n");

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
    fprintf(stderr,"Could not get new cell value to use for edge of group\n");
  }

  fprintf(stdout,"NUMBER OF FIX EDGES: %d\n", fixEdges.size());
  fprintf(stdout,"FIX EDGES: \n");
  for (int j=0; j<fixEdges.size(); j++)
  {
    int edgeSize = region.BoundaryEdges[fixEdges[j]].size();
    int startPtId = region.BoundaryEdges[fixEdges[j]][0];
    int endPtId = region.BoundaryEdges[fixEdges[j]][edgeSize-1];
    fprintf(stdout,"  %d, Start Pt: %d, End Pt: %d\n", fixEdges[j], startPtId, endPtId);
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

  fprintf(stdout,"FILLING IN %d BUT NOT TOUCHING %d and %d\n", newCellValue, newCellValue, polyRegion.IndexCluster);

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
int vtkSVGroupsSegmenter::FixOffsetTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd, std::string arrayName,
                                                const Region region, const Region polyRegion, std::vector<int> allEdges,
                                                std::vector<int> badEdges, vtkIdList *critPts)
{
  fprintf(stdout,"--FIX OFFSET TRIFURCATION--\n");
  int patchValue = region.IndexCluster;
  fprintf(stdout,"CLUSTER %d\n", patchValue);
  int numEdges = region.BoundaryEdges.size();
  fprintf(stdout,"NUM EDGES %d\n", numEdges);

  fprintf(stdout,"NUMBER OF ALL EDGES: %d\n", allEdges.size());
  fprintf(stdout,"ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    fprintf(stdout,"%d ", allEdges[j]);
  fprintf(stdout,"\n");
  fprintf(stdout,"NUMBER OF BAD EDGES: %d\n", badEdges.size());
  fprintf(stdout,"BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    fprintf(stdout,"%d ", badEdges[j]);
  fprintf(stdout,"\n");

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

      fprintf(stdout,"BAD: ");
      for (int l=0; l<badCornerPtList->GetNumberOfIds(); l++)
        fprintf(stdout,"%d ", badCornerPtList->GetId(l));
      fprintf(stdout,"\n");

      vtkNew(vtkIdList, polyCornerPtList);
      vtkSVGeneralUtils::GetPointCellsValues(polyPd, arrayName, polyCornerPtId, polyCornerPtList);

      fprintf(stdout,"POLY: ");
      for (int l=0; l<polyCornerPtList->GetNumberOfIds(); l++)
        fprintf(stdout,"%d ", polyCornerPtList->GetId(l));
      fprintf(stdout,"\n");

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
    fprintf(stderr,"COULD NOT FIND POLYCUBE SET OF EDGES MATCHING BAD EDGES\n");
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

  fprintf(stdout,"NUMBER OF FIX EDGES: %d\n", fixEdges.size());
  fprintf(stdout,"FIX EDGES: ");
  for (int j=0; j<fixEdges.size(); j++)
    fprintf(stdout,"%d ", fixEdges[j]);
  fprintf(stdout,"\n");

  this->FixEdges(pd, origPd, arrayName, region, allEdges, fixEdges, critPts);

  return SV_OK;
}

// ----------------------
// FixPlanarTrifurcation
// ----------------------
int vtkSVGroupsSegmenter::FixPlanarTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                                                const Region region, std::vector<int> allEdges,
                                                std::vector<int> badEdges, vtkIdList *critPts)
{

  fprintf(stdout,"--FIX PLANAR TRIFURCATION--\n");
  int patchValue = region.IndexCluster;
  fprintf(stdout,"CLUSTER %d\n", patchValue);
  int numEdges = region.BoundaryEdges.size();
  fprintf(stdout,"NUM EDGES %d\n", numEdges);

  fprintf(stdout,"NUMBER OF ALL EDGES: %d\n", allEdges.size());
  fprintf(stdout,"ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    fprintf(stdout,"%d ", allEdges[j]);
  fprintf(stdout,"\n");
  fprintf(stdout,"NUMBER OF BAD EDGES: %d\n", badEdges.size());
  fprintf(stdout,"BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    fprintf(stdout,"%d ", badEdges[j]);
  fprintf(stdout,"\n");

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

  fprintf(stdout,"NUMBER OF FIX EDGES: %d\n", fixEdges.size());
  fprintf(stdout,"FIX EDGES: ");
  for (int j=0; j<fixEdges.size(); j++)
    fprintf(stdout,"%d ", fixEdges[j]);
  fprintf(stdout,"\n");

  this->FixEdges(pd, origPd, arrayName, region, allEdges, fixEdges, critPts);

  return SV_OK;
}

// ----------------------
// FixPerpenTrifurcation
// ----------------------
int vtkSVGroupsSegmenter::FixPerpenTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                                                const Region region, std::vector<int> allEdges,
                                                std::vector<int> badEdges, vtkIdList *critPts)
{
  fprintf(stdout,"--FIX PERPENDICULAR TRIFURCATION--\n");
  int patchValue = region.IndexCluster;
  fprintf(stdout,"CLUSTER %d\n", patchValue);
  int numEdges = region.BoundaryEdges.size();
  fprintf(stdout,"NUM EDGES %d\n", numEdges);

  fprintf(stdout,"NUMBER OF ALL EDGES: %d\n", allEdges.size());
  fprintf(stdout,"ALL EDGES: ");
  for (int j=0; j<allEdges.size(); j++)
    fprintf(stdout,"%d ", allEdges[j]);
  fprintf(stdout,"\n");
  fprintf(stdout,"NUMBER OF BAD EDGES: %d\n", badEdges.size());
  fprintf(stdout,"BAD EDGES: ");
  for (int j=0; j<badEdges.size(); j++)
    fprintf(stdout,"%d ", badEdges[j]);
  fprintf(stdout,"\n");

  std::vector<int> fixEdges;
  if ((badEdges[0] + 1) == badEdges[1])
    fixEdges.push_back(badEdges[0]);
  else
    fixEdges.push_back(badEdges[1]);

  fprintf(stdout,"NUMBER OF FIX EDGES: %d\n", fixEdges.size());
  fprintf(stdout,"FIX EDGES: ");
  for (int j=0; j<fixEdges.size(); j++)
    fprintf(stdout,"%d ", fixEdges[j]);
  fprintf(stdout,"\n");

  this->FixEdges(pd, origPd, arrayName, region, allEdges, fixEdges, critPts);

  return SV_OK;
}

// ----------------------
// FixEdges
// ----------------------
int vtkSVGroupsSegmenter::FixEdges(vtkPolyData *pd, vtkPolyData *origPd,
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
      fprintf(stdout,"START ID: %d\n", startPtId);
      fprintf(stdout,"FINAL ID: %d\n", finalPtId);
      fprintf(stdout,"TEST ID: %d\n", testId);
      if (critPts->IsId(testId) == -1)
        halfSize = edgeSize/2;
      else
        halfSize = edgeSize/2-1;
    }
    else
    {
      halfSize = floor(edgeSize/2.);
    }

    fprintf(stdout,"HALF SIZE IS: %d\n", halfSize);
    int halfId = region.BoundaryEdges[badEdgeId][halfSize];
    critPts->InsertUniqueId(halfId);
    fprintf(stdout,"HALF ID: %d\n", halfId);

    for (int k=0; k<allEdges.size(); k++)
    {
      if (allEdges[k] != badEdgeId)
      {
        int allEdgeSize = region.BoundaryEdges[allEdges[k]].size();

        // TODO: Check half size isn't larger then gropuEdgeSize!!!
        if (halfSize > allEdgeSize/2)
        {
          //fprintf(stderr,"We have MAJOR PROBLEMO\n");
          //return SV_ERROR;
          halfSize = allEdgeSize/3;
        }

        int stopId       = -1;
        int edgeCell     = -1;
        int newCellValue = -1;

        if (edgeSize == 2)
        {
          stopId = region.BoundaryEdges[allEdges[k]][halfSize];

          vtkNew(vtkIdList, halfValues);
          vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, halfValues);
          fprintf(stdout,"WHAT ARE HALF VALS: ");
          for (int f=0; f<halfValues->GetNumberOfIds(); f++)
            fprintf(stdout,"%d ", halfValues->GetId(f));
          fprintf(stdout,"\n");

          int findId;
          if (halfId == finalPtId)
            findId = startPtId;
          else
            findId = finalPtId;
          vtkNew(vtkIdList, tmpCell);
          origPd->GetPointCells(findId, tmpCell);
          for (int l=0; l<tmpCell->GetNumberOfIds(); l++)
          {
            int edgeCellValue = origPd->GetCellData()->GetArray(arrayName.c_str())->GetTuple1(tmpCell->GetId(l));
            if (edgeCellValue == patchValue)
              edgeCell = tmpCell->GetId(l);

            if (halfValues->IsId(edgeCellValue) == -1)
              newCellValue = edgeCellValue;
          }
        }
        else
        {
          if (region.BoundaryEdges[allEdges[k]][0] == finalPtId)
          {
            fprintf(stdout,"ONER\n");
            stopId = region.BoundaryEdges[allEdges[k]][halfSize];

            vtkNew(vtkIdList, halfValues);
            vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, halfValues);
            fprintf(stdout,"WHAT ARE HALF VALS: ");
            for (int f=0; f<halfValues->GetNumberOfIds(); f++)
              fprintf(stdout,"%d ", halfValues->GetId(f));
            fprintf(stdout,"\n");

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
            fprintf(stdout,"TWOER\n");
            stopId = region.BoundaryEdges[allEdges[k]][allEdgeSize-halfSize-1];

            vtkNew(vtkIdList, halfValues);
            vtkSVGeneralUtils::GetPointCellsValues(origPd, arrayName, halfId, halfValues);
            fprintf(stdout,"WHAT ARE HALF VALS: ");
            for (int f=0; f<halfValues->GetNumberOfIds(); f++)
              fprintf(stdout,"%d ", halfValues->GetId(f));
            fprintf(stdout,"\n");

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

          if (edgeSize == 2)
          {
            int findId;
            if (halfId == finalPtId)
              findId = startPtId;
            else
              findId = finalPtId;
            vtkNew(vtkIdList, pointCellIds);
            int startId  = branchPd->GetPointData()->GetArray("TmpInternalIds")->LookupValue(findId);
            branchPd->GetPointCells(startId, pointCellIds);

            for (int l=0; l<pointCellIds->GetNumberOfIds(); l++)
            {
              int tmpCellId = pointCellIds->GetId(l);
              branchPd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(tmpCellId, newCellValue);
              int realCellId = branchPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(tmpCellId);
              pd->GetCellData()->GetArray(arrayName.c_str())->SetTuple1(realCellId, newCellValue);
            }
          }
          else
          {
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
  }

  return SV_OK;
}

// ----------------------
// MatchSurfaceToPolycube
// ----------------------
int vtkSVGroupsSegmenter::MatchSurfaceToPolycube()
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

      fprintf(stdout,"PT ID 0: %d\n", ptId0);
      fprintf(stdout,"IDS 0: ");
      for (int f=0; f<ptId0List->GetNumberOfIds(); f++)
        fprintf(stdout,"%d ", ptId0List->GetId(f));
      fprintf(stdout,"\n");
      fprintf(stdout,"PT ID N: %d\n", ptIdN);
      fprintf(stdout,"IDS N: ");
      for (int f=0; f<ptIdNList->GetNumberOfIds(); f++)
        fprintf(stdout,"%d ", ptIdNList->GetId(f));
      fprintf(stdout,"\n");
      fprintf(stdout,"\n");

      vtkNew(vtkIdList, intersectList);
      intersectList->DeepCopy(ptId0List);

      intersectList->IntersectWith(ptIdNList);

      std::vector<int> newSlicePoints;
      if (intersectList->GetNumberOfIds() >= 3)
      {
        // Traditional between sides of groups
        vtkSVGroupsSegmenter::SplitBoundary(this->WorkPd, surfaceRegions[i].BoundaryEdges[j], 3, surfaceRegions[i].IndexCluster,
                                            newSlicePoints);

      }
      else if (intersectList->GetNumberOfIds() == 2)
      {
        // Between center of groups, need to do special
        std::vector<int> newSlicePoints;
        vtkSVGroupsSegmenter::SplitBoundary(this->WorkPd, surfaceRegions[i].BoundaryEdges[j], 3, surfaceRegions[i].IndexCluster,
                                            newSlicePoints);
        //vtkSVGroupsSegmenter::SplitBoundary(this->WorkPd, surfaceRegions[i].BoundaryEdges[j], 2, surfaceRegions[i].IndexCluster,
        //                                    newSlicePoints);
      }
      else
      {
        fprintf(stderr,"Not sure where this case should happen, not implemented\n");
        return SV_ERROR;
      }

      for (int k=0; k<newSlicePoints.size(); k++)
      {
        int pointId = newSlicePoints[k];
        fprintf(stdout,"TRYING TO FIND MATCHER FOR %d\n", pointId);

        vtkNew(vtkIdList, surfaceSlicePtVals);
        vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, pointId, surfaceSlicePtVals);

        fprintf(stdout,"POINT CELL VALUES ARE %d %d\n", surfaceSlicePtVals->GetId(0), surfaceSlicePtVals->GetId(1));

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

          fprintf(stdout,"POLY PT ID 0: %d\n", polyPtId0);
          fprintf(stdout,"POLY IDS 0: ");
          for (int f=0; f<polyPtId0List->GetNumberOfIds(); f++)
            fprintf(stdout,"%d ", polyPtId0List->GetId(f));
          fprintf(stdout,"\n");
          fprintf(stdout,"POLY PT ID N: %d\n", polyPtIdN);
          fprintf(stdout,"POLY IDS N: ");
          for (int f=0; f<polyPtIdNList->GetNumberOfIds(); f++)
            fprintf(stdout,"%d ", polyPtIdNList->GetId(f));
          fprintf(stdout,"\n");
          fprintf(stdout,"\n");

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
            fprintf(stdout,"OKAY, THIS IS MATCHING END POINTS\n");
            fprintf(stdout,"SURFACE PTS: %d %d, POLY PTS: %d %d\n", ptId0, ptIdN, polyPtId0, polyPtIdN);

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
                  fprintf(stdout,"WE FOUND OUR MATCHING POINT!\n");
                  fprintf(stdout,"SURFACE PT: %d, POLY PT: %d\n", pointId, edgePtId);
                  int currValue = this->PolycubePd->GetPointData()->GetArray("SlicePoints")->GetTuple1(edgePtId);
                  if (currValue != -1)
                    fprintf(stdout,"ALREADY SET, MAKE SURE NEW POINT %d MATCHES %d\n", pointId, currValue);
                  else
                  {
                    this->PolycubePd->GetPointData()->GetArray("SlicePoints")->SetTuple1(edgePtId, pointId);
                    this->PolycubePd->GetPointData()->GetArray("SlicePoints")->SetTuple1(polyPtId0, ptId0);
                    this->PolycubePd->GetPointData()->GetArray("SlicePoints")->SetTuple1(polyPtIdN, ptIdN);
                  }
                  k++;
                  if (k == newSlicePoints.size())
                  {
                    fprintf(stdout,"EDGE DONE\n");
                    fprintf(stdout,"\n");
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
          fprintf(stderr,"DIDNT FIND A MATCHING PC POINT FOR SLICE POINT %d\n", pointId);
      }
    }
  }

  //std::string polyfile = "/Users/adamupdegrove/Desktop/tmp/PolycubeWithSurfaceIds.vtp";
  //vtkSVIOUtils::WriteVTPFile(polyfile, this->PolycubePd);

  return SV_OK;
}

// ----------------------
// SplitBoundary
// ----------------------
int vtkSVGroupsSegmenter::SplitBoundary(vtkPolyData *pd,
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
int vtkSVGroupsSegmenter::CheckSlicePoints()
{

  int numPoints = this->PolycubePd->GetNumberOfPoints();
  int numCells = this->PolycubePd->GetNumberOfPoints();

  fprintf(stdout,"TOTAL NUM OF CELLS: %d\n", this->WorkPd->GetNumberOfCells());
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

      fprintf(stdout,"NUMBER OF CONNECTING PATCHES: %d\n", numVals);
      fprintf(stdout,"VALENCE OF SLICE POINT: %d\n", numCells);

      if (numVals >= (1./2)*numCells)
      {
        // Lets split these cells
        fprintf(stdout,"SPLITTING CELLS AROUND POINT: %d\n", slicePointId);
        this->SplitCellsAroundPoint(this->WorkPd, slicePointId);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// SplitCellsAroundPoint
// ----------------------
int vtkSVGroupsSegmenter::SplitCellsAroundPoint(vtkPolyData *pd, int ptId)
{
  vtkNew(vtkIdList, pointCells);
  pd->GetPointCells(ptId, pointCells);

  int numSplitCells = pointCells->GetNumberOfIds();

  fprintf(stdout,"SPLITTING %d CELLS\n", numSplitCells);

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
    fprintf(stdout,"SPLITIING CELL: %d\n", cellId);

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
int vtkSVGroupsSegmenter::SplitEdge(vtkPolyData *pd, int cellId, int ptId0, int ptId1,
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
// ShiftEdgeList
// ----------------------
int vtkSVGroupsSegmenter::ShiftEdgeList(vtkPolyData *branchPd, std::vector<std::vector<int> > &openEdges,
                                        std::vector<std::vector<int> > &shiftedOpenEdges)
{

  for (int j=0; j<openEdges.size(); j++)
  {
    int edgeSize = openEdges[j].size();
    int shiftCount = 0;
    for (int k=0; k<edgeSize-1; k++)
    {
      int edgePtId0 = openEdges[j][k];
      int isSlicePoint = branchPd->GetPointData()->GetArray("SlicePoints")->GetTuple1(edgePtId0);

      if (isSlicePoint == 1)
      {
        std::vector<int> shiftedEdges(edgeSize);
        for (int l=0; l<edgeSize-1; l++)
        {
          shiftedEdges[l] = openEdges[j][(l+shiftCount)%(edgeSize-1)];
        }
        shiftedEdges[edgeSize-1] = edgePtId0;
        shiftedOpenEdges.push_back(shiftedEdges);
        break;
      }
      shiftCount++;
    }
  }

  return SV_OK;
}

// ----------------------
// SplitEdgeList
// ----------------------
int vtkSVGroupsSegmenter::SplitEdgeList(vtkPolyData *branchPd, std::vector<int> &openEdges,
                                        std::vector<std::vector<int> > &splitOpenEdges)
{

  int isFirstSlicePoint = branchPd->GetPointData()->GetArray("SlicePoints")->GetTuple1(openEdges[0]);

  if (isFirstSlicePoint != 1)
  {
    fprintf(stderr,"First point is not a slice point for edge split\n");
    return SV_ERROR;
  }
  std::vector<int> splitEdge;
  splitEdge.push_back(openEdges[0]);

  for (int j=1; j<openEdges.size(); j++)
  {
    int edgePtId = openEdges[j];
    splitEdge.push_back(edgePtId);

    int isSlicePoint = branchPd->GetPointData()->GetArray("SlicePoints")->GetTuple1(edgePtId);

    if (isSlicePoint == 1)
    {
      splitOpenEdges.push_back(splitEdge);
      splitEdge.clear();
      splitEdge.push_back(edgePtId);
    }
  }

  if (splitOpenEdges.size() != 4)
  {
    fprintf(stderr,"There should be four edges, but got %d\n", splitOpenEdges.size());
    return SV_ERROR;
  }

  return SV_OK;
}

// ----------------------
// GetNBoundaryRows
// ----------------------
int vtkSVGroupsSegmenter::GetNBoundaryRows(vtkPolyData *pd, const int numRows, vtkPolyData *rowsPd)
{
  // Do boundary cell stuffs
  int numCells = pd->GetNumberOfCells();
  std::vector<int> allCells;
  std::vector<int> cellList;
  std::vector<int> cellBool(numCells);
  for (int i=0; i<numCells; i++)
    cellBool[i] = 0;

  for (int i=0; i<numCells; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];

      vtkNew(vtkIdList, cellNeighbors);
      pd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellNeighbors);
      if (cellNeighbors->GetNumberOfIds() == 0)
      {
        allCells.push_back(i);
        cellList.push_back(i);
        cellBool[i] = 1;
        break;
      }
    }
  }

  fprintf(stdout,"WHATS THE START SIZE: %d\n", allCells.size());
  int startRow=1;
  if (vtkSVGroupsSegmenter::GetNListNeighbors(pd, cellList, cellBool, startRow, numRows, allCells) != SV_OK)
  {
    fprintf(stderr,"Could not get boundary cells\n");
    return SV_ERROR;
  }

  vtkNew(vtkIntArray, keepCellArray);
  keepCellArray->SetNumberOfTuples(numCells);
  keepCellArray->SetName("KeepCellArray");
  keepCellArray->FillComponent(0, 0);

  vtkNew(vtkIdList, letsSeeList);
  for (int i=0; i<allCells.size(); i++)
  {
    keepCellArray->SetTuple1(allCells[i], 1);
    letsSeeList->InsertUniqueId(allCells[i]);
  }

  fprintf(stdout,"PD SIZE: %d\n", numCells);
  fprintf(stdout,"REALITY SIZE: %d\n", letsSeeList->GetNumberOfIds());
  fprintf(stdout,"WHATS THE SIZE: %d\n", allCells.size());

  pd->GetCellData()->AddArray(keepCellArray);

  vtkSVGeneralUtils::ThresholdPd(pd, 1, 1, 1, "KeepCellArray", rowsPd);

  pd->GetCellData()->RemoveArray("KeepCellArray");

  return SV_OK;
}

// ----------------------
// GetNListNeighbors
// ----------------------
int vtkSVGroupsSegmenter::GetNListNeighbors(vtkPolyData *pd, std::vector<int> &cellList,
                                            std::vector<int> &cellBool,
                                            int &currRow, const int numRows,
                                            std::vector<int> &allCells)
{
  int iSize = cellList.size();

  std::vector<int> newCellList;
  for (int i=0; i<iSize; i++)
  {
    int cellId = cellList[i];

    vtkIdType npts, *pts;
    pd->GetCellPoints(cellId, npts, pts);
    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];

      vtkNew(vtkIdList, cellNeighbors);
      pd->GetCellEdgeNeighbors(cellId, ptId0, ptId1, cellNeighbors);

      for (int k=0; k<cellNeighbors->GetNumberOfIds(); k++)
      {
        if (cellBool[cellNeighbors->GetId(k)] == 0)
        {
          allCells.push_back(cellNeighbors->GetId(k));
          newCellList.push_back(cellNeighbors->GetId(k));
          cellBool[cellNeighbors->GetId(k)] = 1;
        }
      }
    }
  }

  if (currRow < numRows)
  {
    currRow++;
    vtkSVGroupsSegmenter::GetNListNeighbors(pd, newCellList, cellBool, currRow, numRows, allCells);
  }

  return SV_OK;
}

// ----------------------
// GetTrueBoundaryDirections
// ----------------------
int vtkSVGroupsSegmenter::GetTrueBoundaryDirections(vtkPolyData *branchPd,
                                                    vtkPolyData *polyBranchPd,
                                                    const int groupId,
                                                    vtkSVPolyBallLine *groupTubes,
                                                    std::vector<std::vector<int> >  &shiftedOpenEdges,
                                                    vtkDoubleArray *avgVecs,
                                                    vtkIntArray *patchDirs)
{

  if (shiftedOpenEdges.size() == 0 && this->Centerlines->GetNumberOfCells() != 1)
  {
    fprintf(stdout,"Vessel does not have any open edges, but it must\n");
    return SV_ERROR;
  }

  avgVecs->SetNumberOfComponents(3);
  avgVecs->SetNumberOfTuples(4*shiftedOpenEdges.size());
  std::vector<std::vector<int> > sliceEdges;
  int count = 0;
  for (int j=0; j<shiftedOpenEdges.size(); j++)
  {
    int edgeSize = shiftedOpenEdges[j].size();
    for (int k=0; k<edgeSize-1 ; k++)
    {
      int isSlicePoint = -1;
      double avgVec[3]; avgVec[0] = 0.0; avgVec[1] = 0.0; avgVec[2] = 0.0;
      int edgeCount  = 0;
      std::vector<int> newSliceEdge;
      while(isSlicePoint == -1)
      {
        int edgePtId0 = shiftedOpenEdges[j][k];
        int edgePtId1 = shiftedOpenEdges[j][k+1];

        vtkNew(vtkIdList, edgeCell);
        branchPd->GetCellEdgeNeighbors(-1, edgePtId0, edgePtId1, edgeCell);

        if (edgeCell->GetNumberOfIds() != 1)
        {
          vtkErrorMacro("Something went wrong in getting open edges");
          return SV_ERROR;
        }

        double pts[3][3];
        vtkIdType npts, *ptids;
        branchPd->GetCellPoints(edgeCell->GetId(0), npts, ptids);
        for (int j=0; j<npts; j++)
          branchPd->GetPoint(ptids[j], pts[j]);

        // Get center
        double center[3];
        vtkTriangle::TriangleCenter(pts[0], pts[1], pts[2], center);

        // Evaluate function at point!
        groupTubes->EvaluateFunction(center);

        //TODO: JUST TESTING SOMETHING OUT!!!
        double closePt[3];
        groupTubes->GetLastPolyBallCenter(closePt);

        double vec[3];
        vtkMath::Subtract(center, closePt, vec);
        vtkMath::Normalize(vec);

        for (int l=0; l<3; l++)
          avgVec[l] += vec[l];
        edgeCount++;

        isSlicePoint = branchPd->GetPointData()->GetArray("SlicePoints")->GetTuple1(edgePtId1);

        newSliceEdge.push_back(edgePtId0);

        if (isSlicePoint == -1)
          k++;
        else
          newSliceEdge.push_back(edgePtId1);
      }
      vtkMath::MultiplyScalar(avgVec, 1./(edgeCount));
      avgVecs->SetTuple(count, avgVec);
      count++;
      sliceEdges.push_back(newSliceEdge);
    }

  }

  patchDirs->SetNumberOfComponents(1);
  patchDirs->SetNumberOfTuples(sliceEdges.size());
  for (int j=0; j<sliceEdges.size(); j++)
  {
    fprintf(stdout,"LOOKING AT EDGE %d\n", j);
    int edgeSize = sliceEdges[j].size();

    int edgePtId0 = branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(sliceEdges[j][0]);
    int edgePtIdN = branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(sliceEdges[j][edgeSize-1]);
    fprintf(stdout,"START POINT: %d\n", edgePtId0);
    fprintf(stdout,"ENDER POINT: %d\n", edgePtIdN);

    int polyPtId0 = polyBranchPd->GetPointData()->GetArray("SlicePoints")->LookupValue(edgePtId0);
    int polyPtIdN = polyBranchPd->GetPointData()->GetArray("SlicePoints")->LookupValue(edgePtIdN);

    if (polyPtId0 == -1 || polyPtIdN == -1)
    {
      vtkErrorMacro("Slice point does not exist on polycube!");
      return SV_ERROR;
    }

    vtkNew(vtkIdList, polyPtId0List);
    vtkSVGeneralUtils::GetPointCellsValues(polyBranchPd, "PatchIds", polyPtId0, polyPtId0List);

    vtkNew(vtkIdList, polyPtIdNList);
    vtkSVGeneralUtils::GetPointCellsValues(polyBranchPd, "PatchIds", polyPtIdN, polyPtIdNList);

    polyPtId0List->IntersectWith(polyPtIdNList);

    if (polyPtId0List->GetNumberOfIds() != 1)
    {
      vtkErrorMacro("Should be one patch between slice points, there is %d\n");
      return SV_ERROR;
    }

    int patchId = polyPtId0List->GetId(0);
    int patchDir = patchId%6;
    fprintf(stdout,"PATCH DIR IS: %d\n", patchDir);
    patchDirs->SetTuple1(j, patchDir);
  }

  if (sliceEdges.size() == 8)
  {
    // See if we need to switch the order of everything
    int edgePtId = branchPd->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(sliceEdges[0][0]);
    fprintf(stdout,"CHECKING IF RIGHT ORDER USING FIRST POINT: %d\n", edgePtId);

    vtkNew(vtkIdList, edgePtIdList);
    vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->GroupIdsArrayName, edgePtId, edgePtIdList);

    if (edgePtIdList->GetNumberOfIds() != 2)
    {
      fprintf(stderr,"Need to have two ids here\n");
      return SV_ERROR;
    }

    int touchingGroupId;
    if (edgePtIdList->GetId(0) == groupId)
      touchingGroupId = edgePtIdList->GetId(1);
    else
      touchingGroupId = edgePtIdList->GetId(0);
    fprintf(stdout,"WHAT IS TOUCHING GROUP ID: %d\n", touchingGroupId);

    int centerlineId = this->MergedCenterlines->GetCellData()->GetArray(this->GroupIdsArrayName)->LookupValue(groupId);

    vtkIdType npts, *pts;
    this->MergedCenterlines->GetCellPoints(centerlineId, npts, pts);

    vtkNew(vtkIdList, touchingClGroupIdList);
    this->MergedCenterlines->GetPointCells(pts[0], touchingClGroupIdList);

    int isBeginning = 0;
    for (int i=0; i<touchingClGroupIdList->GetNumberOfIds(); i++)
    {
      int centerlineGroupId = this->MergedCenterlines->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(touchingClGroupIdList->GetId(i));
      fprintf(stdout,"WHAT IS THE CENTERLINE LIST WE FOUND %d\n", centerlineGroupId);
      if (centerlineGroupId == touchingGroupId)
        isBeginning = 1;
    }

    // If not beginning, we need to flip everything
    if (!isBeginning)
    {
      fprintf(stdout,"MUST BE FLIPPERED!!!\n");
      vtkNew(vtkDoubleArray, tmpAvgVecs);
      tmpAvgVecs->DeepCopy(avgVecs);

      vtkNew(vtkIntArray, tmpPatchDirs);
      tmpPatchDirs->DeepCopy(patchDirs);

      for (int i=0; i<4; i++)
      {
        avgVecs->SetTuple(i, tmpAvgVecs->GetTuple(i+4));
        avgVecs->SetTuple(i+4, tmpAvgVecs->GetTuple(i));

        patchDirs->SetTuple1(i, tmpPatchDirs->GetTuple1(i+4));
        patchDirs->SetTuple1(i+4, tmpPatchDirs->GetTuple1(i));
      }
    }
  }

  fprintf(stdout,"LETS SEE THE DIRS!!!:\n");
  for (int i=0; i<patchDirs->GetNumberOfTuples(); i++)
  {
    int patchDir = patchDirs->GetTuple1(i);
    double avgVec[3];
    avgVecs->GetTuple(i, avgVec);
    fprintf(stdout,"  %d has %d dir and corresponds to %.6f %.6f %.6f\n", i, patchDir, avgVec[0], avgVec[1], avgVec[2]);
  }

  return SV_OK;

}

// ----------------------
// GetCellRingNeighbors
// ----------------------
int vtkSVGroupsSegmenter::GetCellRingNeighbors(vtkPolyData *pd, vtkIdList *cellIds,
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

// ----------------------
// GetCellDirectNeighbors
// ----------------------
int vtkSVGroupsSegmenter::GetCellDirectNeighbors(vtkPolyData *pd,
                                                 std::vector<std::vector<int> > &neighbors,
                                                 std::vector<int> &numNeighbors)
{

  int numCells = pd->GetNumberOfCells();
  pd->BuildLinks();

  neighbors.clear();
  numNeighbors.clear();

  // Loop through cells
  for (int i=0; i<numCells; i++)
  {
    // count number of edge neighbors
    int directNeiCount = 0;
    std::vector<int> neighborCells;

    // Get cell points
    vtkIdType *pts, npts;
    pd->GetCellPoints(i, npts, pts);

    // Get cell edge neighbors
    for (int j=0; j<npts; j++)
    {
      int ptId0 = pts[j];
      int ptId1 = pts[(j+1)%npts];

      // Get cell edge neighbors
      vtkNew(vtkIdList, cellEdgeNeighbors);
      pd->GetCellEdgeNeighbors(i, ptId0, ptId1, cellEdgeNeighbors);
      directNeiCount += cellEdgeNeighbors->GetNumberOfIds();
      for (int k=0; k<cellEdgeNeighbors->GetNumberOfIds(); k++)
      {
        neighborCells.push_back(cellEdgeNeighbors->GetId(k));
      }
    }
    neighbors.push_back(neighborCells);
    numNeighbors.push_back(directNeiCount);
  }

  return SV_OK;
}

