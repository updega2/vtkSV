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

#include "vtkSVParameterizeSurfaceOnPolycube.h"

#include "vtkAppendPolyData.h"
#include "vtkAppendFilter.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkExecutive.h"
#include "vtkErrorCode.h"
#include "vtkCellArray.h"
#include "vtkCellLocator.h"
#include "vtkIdFilter.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkLine.h"
#include "vtkLinearSubdivisionFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkThreshold.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVersion.h"

#include "vtkSVCleanUnstructuredGrid.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVIOUtils.h"
#include "vtkSVMathUtils.h"
#include "vtkSVPlanarMapper.h"
#include "vtkSVPointSetBoundaryMapper.h"
#include "vtkSVSurfaceMapper.h"

#include <algorithm>

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVParameterizeSurfaceOnPolycube);

// ----------------------
// Constructor
// ----------------------
vtkSVParameterizeSurfaceOnPolycube::vtkSVParameterizeSurfaceOnPolycube()
{
  this->WorkPd = vtkPolyData::New();
  this->SurfaceOnPolycubePd = vtkPolyData::New();
  this->PolycubeOnSurfacePd = vtkPolyData::New();
  this->PolycubePd = NULL;
  this->PolycubeUg = NULL;

  this->GroupIdsArrayName  = NULL;
  this->PatchIdsArrayName  = NULL;
  this->GridIdsArrayName   = NULL;
}

// ----------------------
// Destructor
// ----------------------
vtkSVParameterizeSurfaceOnPolycube::~vtkSVParameterizeSurfaceOnPolycube()
{
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->SurfaceOnPolycubePd != NULL)
  {
    this->SurfaceOnPolycubePd->Delete();
    this->SurfaceOnPolycubePd = NULL;
  }
  if (this->PolycubeOnSurfacePd != NULL)
  {
    this->PolycubeOnSurfacePd->Delete();
    this->PolycubeOnSurfacePd = NULL;
  }
  if (this->PolycubePd != NULL)
  {
    this->PolycubePd->Delete();
    this->PolycubePd = NULL;
  }
  if (this->PolycubeUg != NULL)
  {
    this->PolycubeUg->Delete();
    this->PolycubeUg = NULL;
  }

  if (this->GroupIdsArrayName != NULL)
  {
    delete [] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
  }

  if (this->PatchIdsArrayName != NULL)
  {
    delete [] this->PatchIdsArrayName;
    this->PatchIdsArrayName = NULL;
  }

  if (this->GridIdsArrayName != NULL)
  {
    delete [] this->GridIdsArrayName;
    this->GridIdsArrayName = NULL;
  }
}

// ----------------------
// RequestData
// ----------------------
int vtkSVParameterizeSurfaceOnPolycube::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  this->WorkPd->DeepCopy(input);
  this->WorkPd->BuildLinks();

  // Prep work for filter
  if (this->PrepFilter() != SV_OK)
  {
    vtkErrorMacro("Prep of filter failed");
    this->SetErrorCode(vtkErrorCode::UserError + 1);
    return SV_ERROR;
  }

  // Run the filter
  if (this->RunFilter() != SV_OK)
  {
    vtkErrorMacro("Filter failed");
    this->SetErrorCode(vtkErrorCode::UserError + 2);
    return SV_ERROR;
  }

  output->DeepCopy(this->SurfaceOnPolycubePd);

  return SV_OK;
}

// ----------------------
// PrepFilter
// ----------------------
int vtkSVParameterizeSurfaceOnPolycube::PrepFilter()
{
  if (!this->GroupIdsArrayName)
  {
    vtkDebugMacro("GroupIds Array Name not given, setting to GroupIds");
    this->GroupIdsArrayName = new char[strlen("GroupIds") + 1];
    strcpy(this->GroupIdsArrayName, "GroupIds");
  }

  if (!this->PatchIdsArrayName)
  {
    vtkDebugMacro("PatchIds Array Name not given, setting to PatchIds");
    this->PatchIdsArrayName = new char[strlen("PatchIds") + 1];
    strcpy(this->PatchIdsArrayName, "PatchIds");
  }

  if (this->PolycubePd == NULL)
  {
    vtkErrorMacro("Surface polycube not provided");
    return SV_ERROR;
  }
  this->PolycubePd->BuildLinks();

  if (this->PolycubePd->GetNumberOfCells() == 0)
  {
    vtkErrorMacro("Surface polycube is empty");
    return SV_ERROR;
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 1, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "Group Ids Array with name " << this->GroupIdsArrayName << " does not exist on input surface");
    return SV_OK;
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 1, this->PatchIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "Patch Ids Array with name " << this->PatchIdsArrayName << " specified does not exist on input surface");
    return SV_OK;
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->PolycubePd, 1, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "Group Ids Array with name " << this->GroupIdsArrayName << " does not exist on surface polycube");
    return SV_OK;
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->PolycubePd, 1, this->PatchIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "Patch Ids Array with name " << this->PatchIdsArrayName << " specified does not exist on surface polycube");
    return SV_OK;
  }

  if (this->PolycubeUg == NULL)
  {
    vtkErrorMacro("Volume polycube not provided");
    return SV_ERROR;
  }

  if (this->PolycubeUg->GetNumberOfCells() == 0)
  {
    vtkErrorMacro("Volume polycube is empty");
    return SV_ERROR;
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->PolycubeUg, 1, this->GroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "Group Ids Array with name " << this->GroupIdsArrayName << " does not exist on volume polycube");
    return SV_OK;
  }

  if (!this->GridIdsArrayName)
  {
    vtkDebugMacro("Grid point ids array name not given, setting to GridIds");
    this->GridIdsArrayName = new char[strlen("GridIds") + 1];
    strcpy(this->GridIdsArrayName, "GridIds");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->PolycubeUg, 0, this->GridIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "Grid point ids array with name " << this->GridIdsArrayName << " does not exist on volume polycube");
    return SV_OK;
  }

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVParameterizeSurfaceOnPolycube::RunFilter()
{
  std::vector<Region> patches;
  if (vtkSVGeneralUtils::GetRegions(this->WorkPd, this->PatchIdsArrayName, patches) != SV_OK)
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
  vtkDebugMacro("JUST CHECK: " <<  polycubePd->GetNumberOfPoints());

  int numPatches = patches.size();

  vtkNew(vtkAppendPolyData, appender);

  for (int i=0; i<numPatches; i++)
  {
    if (patches[i].CornerPoints.size() < 4 || patches[i].CornerPoints.size() > 6)
    {
      vtkErrorMacro("Patch should haver either 4, 5, or 6 points, but patch " << patches[i].IndexCluster << " has " << patches[i].CornerPoints.size() << " points");
      return SV_ERROR;
    }

    int groupId = this->WorkPd->GetCellData()->GetArray(
     this->GroupIdsArrayName)->GetTuple1(patches[i].Elements[0]);
    int patchId = this->WorkPd->GetCellData()->GetArray(
     this->PatchIdsArrayName)->GetTuple1(patches[i].Elements[0]);
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
    vtkSVGeneralUtils::ThresholdPd(thresholdPd, patchId, patchId, 1, this->PatchIdsArrayName);

    // Set up boundary mapper
    vtkNew(vtkIntArray, boundaryCorners);
    boundaryCorners->SetNumberOfTuples(patches[i].CornerPoints.size());

    vtkNew(vtkIntArray, paraBoundaryCorners);
    paraBoundaryCorners->SetNumberOfTuples(patches[i].CornerPoints.size());
    vtkDebugMacro("PATCH: " <<  patchId);

    vtkDebugMacro("Corner Points: ");
    for (int j=0; j<patches[i].CornerPoints.size(); j++)
    {
      int ptId = patches[i].CornerPoints[j];
      vtkDebugMacro(" " <<  ptId);

      // Thresholded pt id
      int thresholdPtId = thresholdPd->GetPointData()->GetArray(
        "TmpInternalIds")->LookupValue(ptId);
      boundaryCorners->SetTuple1(j, thresholdPtId);

      // Paramteric space pt id
      vtkNew(vtkIdList, patchVals);
      vtkSVGeneralUtils::GetPointCellsValues(this->WorkPd, this->PatchIdsArrayName, ptId, patchVals);
      int paraPtId = -1;
      if (vtkSVGeneralUtils::FindPointMatchingValues(rotPolycube, this->PatchIdsArrayName, patchVals, paraPtId) != SV_OK)
      {
        vtkErrorMacro("Could not find corresponding polycube point id");
        vtkErrorMacro("Make sure that the polycube and the model match");
        vtkErrorMacro("COULD NOT FIND: ");
        for (int r=0; r<patchVals->GetNumberOfIds(); r++)
          vtkDebugMacro(" " <<  patchVals->GetId(r));
        vtkDebugMacro(" ");
        return SV_ERROR;
      }

      paraBoundaryCorners->SetTuple1(j, paraPtId);
    }
    vtkDebugMacro(" ");

    vtkDebugMacro("Poly Corner Points: ");
    for (int j=0; j<paraBoundaryCorners->GetNumberOfTuples(); j++)
      vtkDebugMacro(" " << paraBoundaryCorners->GetTuple1(j));
    vtkDebugMacro(" ");

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
  {
    vtkWarningMacro("Cleaned mapped polycube surface and input pd do not have the same number of points");
    vtkWarningMacro("MAPPED POINTS: " << tmpPd->GetNumberOfPoints() << " MAPPED CELLS: " << tmpPd->GetNumberOfCells());
    vtkWarningMacro("REGULA POINTS: " << this->WorkPd->GetNumberOfPoints() << " REGULA CELLS: " << this->WorkPd->GetNumberOfCells());
  }

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

  this->SurfaceOnPolycubePd->SetPoints(fullMapPoints);
  this->SurfaceOnPolycubePd->SetPolys(fullMapCells);
  this->SurfaceOnPolycubePd->GetPointData()->PassData(newPointData);
  this->SurfaceOnPolycubePd->GetCellData()->PassData(newCellData);
  this->SurfaceOnPolycubePd->BuildLinks();

  // Make sure volume has Internal Ids
  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(this->PolycubeUg);
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
  cleaner2->SetInputData(this->PolycubeUg);
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

  // all data on fullMapPd now
  //this->InterpolateMapOntoTarget(polycubePd, this->WorkPd, this->SurfaceOnPolycubePd, mappedPd, this->GroupIdsArrayName);
  this->InterpolateMapOntoTarget(paraHexSurface, this->WorkPd, this->SurfaceOnPolycubePd, this->PolycubeOnSurfacePd, this->GroupIdsArrayName);

  return SV_OK;
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVParameterizeSurfaceOnPolycube::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  if (this->GroupIdsArrayName != NULL)
    os << indent << "Group ids array name: " << this->GroupIdsArrayName << "\n";
}

// ----------------------
// RotateGroupToGlobalAxis
// ----------------------
int vtkSVParameterizeSurfaceOnPolycube::RotateGroupToGlobalAxis(vtkPolyData *pd,
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
// InterpolateMapOntoTarget
// ----------------------
int vtkSVParameterizeSurfaceOnPolycube::InterpolateMapOntoTarget(vtkPolyData *sourceBasePd,
                                                         vtkPolyData *targetPd,
                                                         vtkPolyData *targetBasePd,
                                                         vtkPolyData *mappedPd,
                                                         std::string dataMatchingArrayName)
{
  vtkNew(vtkSVSurfaceMapper, interpolator);
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
// GetInteriorPointMaps
// ----------------------
int vtkSVParameterizeSurfaceOnPolycube::GetInteriorPointMaps(vtkPolyData *pdWithAllInterior,
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
// RemoveInteriorCells
// ----------------------
int vtkSVParameterizeSurfaceOnPolycube::RemoveInteriorCells(vtkPolyData *quadMesh)
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
