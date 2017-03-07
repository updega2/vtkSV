/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtkSVGroupsClipper.cxx,v $
Language:  C++
Date:      $Date: 2006/04/06 16:46:43 $
Version:   $Revision: 1.9 $

  Copyright (c) Luca Antiga, David Steinman. All rights reserved.
  See LICENCE file for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm
  for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "vtkSVGroupsClipper.h"
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
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkFeatureEdges.h"
#include "vtkGenericCell.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkSVPolyBallLine.h"
#include "vtkMath.h"
#include "vtkMergeCells.h"
#include "vtkSphere.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"
#include "vtkThreshold.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVersion.h"
#include "vtkXMLPolyDataWriter.h"

vtkStandardNewMacro(vtkSVGroupsClipper);

vtkSVGroupsClipper::vtkSVGroupsClipper()
{
  this->Centerlines = NULL;
  this->CenterlineGroupIdsArrayName = NULL;
  this->CenterlineRadiusArrayName = NULL;
  this->GroupIdsArrayName = NULL;
  this->BlankingArrayName = NULL;
  this->CenterlineGroupIds = NULL;
  this->ClipAllCenterlineGroupIds = 0;
  this->CutoffRadiusFactor = VTK_SV_LARGE_DOUBLE;
  this->ClipValue = 0.0;
  this->UseRadiusInformation = 1;

  this->SetNumberOfOutputPorts(2);
  this->GenerateClippedOutput = 0;
  vtkPolyData *output2 = vtkPolyData::New();
  this->GetExecutive()->SetOutputData(1, output2);
  output2->Delete();
}

vtkSVGroupsClipper::~vtkSVGroupsClipper()
{
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
    delete[] this->CenterlineGroupIdsArrayName;
    this->CenterlineGroupIdsArrayName = NULL;
    }

  if (this->CenterlineRadiusArrayName)
    {
    delete[] this->CenterlineRadiusArrayName;
    this->CenterlineRadiusArrayName = NULL;
    }

  if (this->GroupIdsArrayName)
    {
    delete[] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
    }

  if (this->BlankingArrayName)
    {
    delete[] this->BlankingArrayName;
    this->BlankingArrayName = NULL;
    }
}

vtkPolyData *vtkSVGroupsClipper::GetClippedOutput()
{
  if (this->GetNumberOfOutputPorts() < 2)
    {
    return NULL;
    }

  return vtkPolyData::SafeDownCast(this->GetExecutive()->GetOutputData(1));
}

int vtkSVGroupsClipper::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataArray* centerlineGroupIdsArray;
  vtkIntArray* blankingArray;
  vtkIntArray* groupIdsArray;

  if (!this->Centerlines)
    {
    vtkErrorMacro(<< "Centerlines not set.");
    return 1;
    }

  if (!this->ClipAllCenterlineGroupIds && !this->CenterlineGroupIds)
    {
    vtkErrorMacro(<< "CenterlineGroupIds not set.");
    return 1;
    }

  if (!this->CenterlineGroupIdsArrayName)
    {
    vtkErrorMacro(<< "CenterlineGroupIdsArrayName not set.");
    return 1;
    }

  if (!this->GroupIdsArrayName)
    {
    vtkErrorMacro(<< "GroupIdsArrayName not set.");
    return 1;
    }

  centerlineGroupIdsArray = this->Centerlines->GetCellData()->GetArray(this->CenterlineGroupIdsArrayName);

  if (!centerlineGroupIdsArray)
    {
    vtkErrorMacro(<< "CenterlineGroupIdsArray with name specified does not exist");
    return 1;
    }

  if (!this->BlankingArrayName)
    {
    vtkErrorMacro(<< "BlankingArrayName not set.");
    return 1;
    }

  blankingArray = vtkIntArray::SafeDownCast(this->Centerlines->GetCellData()->GetArray(this->BlankingArrayName));

  if (!blankingArray)
    {
    vtkErrorMacro(<< "BlankingArrayName with name specified does not exist");
    return 1;
    }

  if (!this->CenterlineRadiusArrayName)
    {
    vtkErrorMacro(<< "CenterlineRadiusArrayName not set.");
    return 1;
    }

  if (!this->Centerlines->GetPointData()->GetArray(this->CenterlineRadiusArrayName))
    {
    vtkErrorMacro(<< "CenterlineRadiusArray with name specified does not exist");
    return 1;
    }

  if (this->Centerlines->GetNumberOfCells() == 1)
    {
    output->DeepCopy(input);
    groupIdsArray = vtkIntArray::New();
    groupIdsArray->SetName(this->GroupIdsArrayName);
    groupIdsArray->SetNumberOfTuples(output->GetNumberOfCells());
    groupIdsArray->FillComponent(0,centerlineGroupIdsArray->GetComponent(0,0));
    output->GetCellData()->AddArray(groupIdsArray);
    groupIdsArray->Delete();
    return 1;
    }

  //if (this->GenerateClippedOutput)
  //  {
  //  }

  // for each group, compute the clipping array, clip, add group ids array and append.

  vtkNew(vtkSVPolyBallLine, groupTubes);
  groupTubes->SetInput(this->Centerlines);
  groupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  groupTubes->SetUseRadiusInformation(this->UseRadiusInformation);

  vtkNew(vtkSVPolyBallLine, nonGroupTubes);
  nonGroupTubes->SetInput(this->Centerlines);
  nonGroupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  nonGroupTubes->SetUseRadiusInformation(this->UseRadiusInformation);


  double point[3];
  double groupTubeValue, nonGroupTubeValue, tubeDifferenceValue;

  vtkIdType groupId;

  vtkNew(vtkIdList, centerlineGroupIds);
  int i;
  if (this->ClipAllCenterlineGroupIds)
    {
    for (i=0; i<centerlineGroupIdsArray->GetNumberOfTuples(); i++)
      {
      if (blankingArray->GetValue(i) == 1)
        {
        continue;
        }
      centerlineGroupIds->InsertUniqueId(static_cast<vtkIdType>(vtkMath::Round(centerlineGroupIdsArray->GetComponent(i,0))));
      }
    }
  else
    {
    centerlineGroupIds->DeepCopy(this->CenterlineGroupIds);
    }

  vtkNew(vtkPolyData, clippingInput);
  clippingInput->DeepCopy(input);

  vtkIntArray *startGroupIds = vtkIntArray::New();
  startGroupIds->SetName(this->GroupIdsArrayName);
  startGroupIds->SetNumberOfTuples(clippingInput->GetNumberOfCells());
  startGroupIds->FillComponent(0,-1);
  clippingInput->GetCellData()->AddArray(startGroupIds);
  vtkNew(vtkAppendPolyData, appendBranches);
  int numberOfPoints = clippingInput->GetNumberOfPoints();

  // Set up clipping arrays
  for (i=0; i<centerlineGroupIds->GetNumberOfIds(); i++)
    {
    groupId = centerlineGroupIds->GetId(i);

    std::stringstream groupstr;
    groupstr << groupId;
    std::string clipName = "ClippingArray_"+groupstr.str();

    vtkNew(vtkDoubleArray, clippingArray);
    clippingArray->SetNumberOfComponents(1);
    clippingArray->SetNumberOfTuples(numberOfPoints);
    clippingArray->FillComponent(0,0.0);
    clippingArray->SetName(clipName.c_str());

    clippingInput->GetPointData()->AddArray(clippingArray);

    vtkNew(vtkIdList, groupTubesGroupIds);
    vtkNew(vtkIdList, nonGroupTubesGroupIds);
    groupTubesGroupIds->Initialize();
    nonGroupTubesGroupIds->Initialize();

    for (int j=0; j<this->Centerlines->GetNumberOfCells(); j++)
      {
      if (blankingArray->GetValue(j) == 1)
        {
        continue;
        }
      if (static_cast<vtkIdType>(vtkMath::Round(centerlineGroupIdsArray->GetComponent(j,0))) == groupId)
        {
        groupTubesGroupIds->InsertNextId(j);
        }
      else
        {
        nonGroupTubesGroupIds->InsertNextId(j);
        }
      }

    if ((groupTubesGroupIds->GetNumberOfIds() == 0) || (nonGroupTubesGroupIds->GetNumberOfIds() == 0))
      {
      continue;
      }

    groupTubes->SetInputCellIds(groupTubesGroupIds);
    nonGroupTubes->SetInputCellIds(nonGroupTubesGroupIds);

    for (int k=0; k<numberOfPoints; k++)
      {
      clippingInput->GetPoint(k,point);
      groupTubeValue = groupTubes->EvaluateFunction(point);
      if (groupTubeValue > this->CutoffRadiusFactor * this->CutoffRadiusFactor - 1)
        {
        groupTubeValue = VTK_SV_LARGE_DOUBLE;
        }
      nonGroupTubeValue = nonGroupTubes->EvaluateFunction(point);
      tubeDifferenceValue = nonGroupTubeValue - groupTubeValue;
      clippingArray->SetValue(k,tubeDifferenceValue);
      }

    }

  for (i=0; i<centerlineGroupIds->GetNumberOfIds(); i++)
    {
    groupId = centerlineGroupIds->GetId(i);
    std::stringstream groupstr;
    groupstr << groupId;
    std::string clipName = "ClippingArray_"+groupstr.str();
    clippingInput->GetPointData()->SetActiveScalars(clipName.c_str());

    vtkNew(vtkClipPolyData, clipper);
#if (VTK_MAJOR_VERSION <= 5)
    clipper->SetInput(clippingInput);
#else
    clipper->SetInputData(clippingInput);
#endif
    clipper->SetValue(this->ClipValue);
    clipper->GenerateClipScalarsOff();
    clipper->SetGenerateClippedOutput(1);
    clipper->Update();

    if (clipper->GetOutput()->GetNumberOfPoints()==0)
      {
      continue;
      }

    vtkNew(vtkPolyData, clippedBranch);
    clippedBranch->DeepCopy(clipper->GetOutput());

    vtkNew(vtkPolyData, clippedOutputBranch);
    clippedOutputBranch->DeepCopy(clipper->GetClippedOutput());

    vtkNew(vtkCleanPolyData, cleaner);
    cleaner->SetInputData(clippedBranch);
    cleaner->Update();
    vtkNew(vtkTriangleFilter, triangulator);
    triangulator->SetInputData(cleaner->GetOutput());
    triangulator->Update();
    clippedBranch->DeepCopy(triangulator->GetOutput());
    this->ReplaceDataOnCells(clippedBranch, groupId, -1, this->GroupIdsArrayName);
#if (VTK_MAJOR_VERSION <= 5)
    appendBranches->AddInput(clippedBranch);
#else
    appendBranches->AddInputData(clippedBranch);
#endif
    std::stringstream iterstr;
    iterstr << groupId;
    std::string iterName = "/Users/adamupdegrove/Desktop/tmp/Clipper";
    std::string filename =  iterName+"_"+iterstr.str()+".vtp";
    vtkNew(vtkXMLPolyDataWriter, writer);
    writer->SetInputData(clippedBranch);
    writer->SetFileName(filename.c_str());
    writer->Write();

    vtkNew(vtkCleanPolyData, cleanerClipped);
    cleanerClipped->SetInputData(clippedOutputBranch);
    cleanerClipped->Update();
    vtkNew(vtkTriangleFilter, triangulatorClipped);
    triangulatorClipped->SetInputData(cleanerClipped->GetOutput());
    triangulatorClipped->Update();
    clippedOutputBranch->DeepCopy(triangulatorClipped->GetOutput());

    clippingInput->DeepCopy(clippedOutputBranch);
    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(clippingInput);
    ider->SetIdsArrayName("TmpInternalIds");
    ider->Update();
    vtkNew(vtkFeatureEdges, boundaries);
    boundaries->SetInputData(ider->GetOutput());
    boundaries->BoundaryEdgesOn();
    boundaries->FeatureEdgesOff();
    boundaries->NonManifoldEdgesOff();
    boundaries->ManifoldEdgesOff();
    boundaries->Update();

    vtkDataArray *internalIds = boundaries->GetOutput()->GetPointData()->
      GetArray("TmpInternalIds");
    for (int j=0; j<centerlineGroupIds->GetNumberOfIds(); j++)
    {
      int fixGroupId = centerlineGroupIds->GetId(j);
      if (fixGroupId > groupId)
      {
        std::stringstream fixgroupstr;
        fixgroupstr << fixGroupId;
        std::string fixName = "ClippingArray_"+fixgroupstr.str();
        vtkDataArray *clipVals = boundaries->GetOutput()->GetPointData()->
          GetArray(fixName.c_str());
        for (int k=0; k<boundaries->GetOutput()->GetNumberOfPoints(); k++)
        {
          double val = clipVals->GetTuple1(k);
          if (val <= 1.0e-6 && val >= -1.0e-6)
          {
            int pointId = internalIds->GetTuple1(k);
            clippingInput->GetPointData()->GetArray(fixName.c_str())->
              SetTuple1(pointId, 0.0);
          }
        }
      }
    }


    }
  appendBranches->AddInputData(clippingInput);
  appendBranches->Update();

  vtkNew(vtkCleanPolyData, allCleaner);
  allCleaner->SetInputData(appendBranches->GetOutput());
  allCleaner->Update();
  vtkNew(vtkPolyData, finisher);
  finisher->DeepCopy(allCleaner->GetOutput());

  vtkNew(vtkIdList, separateIds);
  this->FindGroupSeparatingPoints(finisher, separateIds);

  vtkNew(vtkPoints, newPoints);
  this->SplitGroups(finisher, separateIds, newPoints);

  this->FillGroups(finisher, newPoints);
  vtkNew(vtkCleanPolyData, finalCleaner);
  finalCleaner->SetInputData(finisher);
  finalCleaner->Update();
  finisher->DeepCopy(finalCleaner->GetOutput());

  output->DeepCopy(finisher);

  return 1;
}

void vtkSVGroupsClipper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGroupsClipper::ReplaceDataOnCells(vtkPointSet *pointset,
                                                               const int replaceVal,
                                                               const int currVal,
                                                               const std::string &arrName)
{
  int numCells = pointset->GetNumberOfCells();
  vtkIntArray *cellIds = vtkIntArray::SafeDownCast(pointset->GetCellData()->GetArray(arrName.c_str()));

  for (int i=0; i<numCells; i++)
  {
    int val = cellIds->GetValue(i);
    if (val == currVal)
    {
      cellIds->SetValue(i, replaceVal);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGroupsClipper::FindGroupSeparatingPoints(vtkPolyData *pd,
                                                                      vtkIdList *separateIds)
{
  pd->BuildLinks();
  int numPoints = pd->GetNumberOfPoints();
  vtkIntArray *groupIdsArray =
    vtkIntArray::SafeDownCast(pd->GetCellData()->GetArray(this->GroupIdsArrayName));

  for (int i=0; i<numPoints; i++)
  {
    vtkNew(vtkIdList, groupIds);
    this->GetPointGroups(pd, this->GroupIdsArrayName, i, groupIds);
    int pointType = groupIds->GetNumberOfIds();
    if (pointType == 3)
    {
      separateIds->InsertNextId(i);
    }
  }
  fprintf(stdout,"Number Of Ids: %lld\n", separateIds->GetNumberOfIds());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGroupsClipper::GetPointGroups(vtkPolyData *pd, std::string arrayName,
                                                  const int pointId, vtkIdList *groupIds)
{
  vtkDataArray *groupIdsArray =
    pd->GetCellData()->GetArray(arrayName.c_str());

  vtkNew(vtkIdList, cellIds);
  pd->GetPointCells(pointId, cellIds);
  groupIds->Reset();
  for (int i=0; i<cellIds->GetNumberOfIds(); i++)
  {
    int groupValue = groupIdsArray->GetTuple1(cellIds->GetId(i));
    if (groupIds->IsId(groupValue) == -1)
    {
      groupIds->InsertNextId(groupValue);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGroupsClipper::SplitGroups(vtkPolyData *pd, vtkIdList *separateIds, vtkPoints *newPoints)
{
  int numIds = separateIds->GetNumberOfIds();
  vtkNew(vtkIntArray, used);
  used->SetNumberOfComponents(0);
  used->SetNumberOfTuples(numIds);
  used->FillComponent(0, -1);

  vtkNew(vtkPoints, separatePoints);
  separatePoints->SetNumberOfPoints(numIds);
  for (int i=0; i<numIds; i++)
  {
    double pt[3];
    pd->GetPoint(separateIds->GetId(i), pt);
    separatePoints->SetPoint(i, pt);
  }

  vtkNew(vtkPolyData, tmpPoly);
  tmpPoly->SetPoints(separatePoints);
  vtkNew(vtkPointLocator, separatePointLocator);
  separatePointLocator->SetDataSet(tmpPoly);
  separatePointLocator->BuildLocator();

  vtkNew(vtkPolyData, pieces);
  this->ThresholdPd(pd, -1, -1, 1, this->GroupIdsArrayName, pieces);

  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(pieces);
  connector->SetExtractionModeToAllRegions();
  connector->ColorRegionsOn();
  connector->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  for (int i=0; i<connector->GetNumberOfExtractedRegions(); i++)
  {
    vtkNew(vtkPolyData, onePiece);
    this->ThresholdPd(surfacer->GetOutput(), i, i, 1, "RegionId", onePiece);
    double center[3];
    this->ComputeMassCenter(onePiece, center);
    vtkNew(vtkIdList, closePointIds);
    separatePointLocator->FindClosestNPoints(3, center, closePointIds);
    double pts[3][3];
    for (int j=0; j<3; j++)
    {
      separatePoints->GetPoint(closePointIds->GetId(j), pts[j]);
    }
    double tmp[3], circleCenter[3];
    vtkMath::Add(pts[0], pts[1], tmp);
    vtkMath::Add(tmp, pts[2], circleCenter);
    vtkMath::MultiplyScalar(circleCenter, 1.0/3);
    newPoints->InsertNextPoint(circleCenter);
    double circleRadius = 0.0;
    for (int j=0; j<onePiece->GetNumberOfPoints(); j++)
    {
      double testPt[3];
      onePiece->GetPoint(j, testPt);
      double dist = std::sqrt(std::pow(testPt[0] - circleCenter[0], 2.0) +
                              std::pow(testPt[1] - circleCenter[1], 2.0) +
                              std::pow(testPt[2] - circleCenter[2], 2.0));
      if (dist > circleRadius)
        circleRadius = dist;
    }
    vtkNew(vtkSphere, sphere);
    sphere->SetCenter(circleCenter);
    sphere->SetRadius(circleRadius);

    vtkNew(vtkClipPolyData, cutter);
    cutter->SetInputData(pd);
    cutter->SetClipFunction(sphere);
    cutter->Update();
    pd->DeepCopy(cutter->GetOutput());
  }

  return 1;
}


int vtkSVGroupsClipper::ThresholdPd(vtkPolyData *pd, int minVal,
                                                        int maxVal, int dataType,
                                                        std::string arrayName,
                                                        vtkPolyData *returnPd)
{
  vtkNew(vtkThreshold, thresholder);
  thresholder->SetInputData(pd);
  //Set Input Array to 0 port,0 connection, dataType (0 - point, 1 - cell, and Regions is the type name
  thresholder->SetInputArrayToProcess(0, 0, 0, dataType, arrayName.c_str());
  thresholder->ThresholdBetween(minVal, maxVal);
  thresholder->Update();
  if (thresholder->GetOutput()->GetNumberOfPoints() == 0)
  {
    returnPd = NULL;
    return 0;
  }

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  returnPd->DeepCopy(surfacer->GetOutput());

  return 1;
}


int vtkSVGroupsClipper::FillGroups(vtkPolyData *pd,
                                                     vtkPoints *newPoints)
{
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(pd);
  ider->SetIdsArrayName("TmpInternalIds");
  ider->Update();

  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(ider->GetOutput());
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();

  for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
  {
    double pt[3];
    newPoints->GetPoint(i, pt);
    vtkNew(vtkPolyData, pointLoop);
    this->GetClosestPointConnectedRegion(boundaries->GetOutput(), pt, pointLoop);
    this->FillRegionGroups(pd, pointLoop, pt);
  }

  return 1;
}

int vtkSVGroupsClipper::FillRegionGroups(vtkPolyData *pd,
                                                           vtkPolyData *boundary,
                                                           double center[3])
{
  vtkNew(vtkPointData, newPointData);
  vtkNew(vtkCellData, newCellData);
  newPointData->CopyAllocate(pd->GetPointData(), pd->GetNumberOfPoints() + 1);
  newCellData->CopyAllocate(pd->GetCellData(), pd->GetNumberOfCells() + boundary->GetNumberOfCells());
  for (int i=0; i<pd->GetNumberOfPoints(); i++)
  {
    newPointData->CopyData(pd->GetPointData(), i, i);
  }
  for (int i=0; i<pd->GetNumberOfCells(); i++)
  {
    newCellData->CopyData(pd->GetCellData(), i, i);
  }

  boundary->BuildLinks();
  vtkDataArray *pointIds = boundary->GetPointData()->GetArray("TmpInternalIds");
  vtkDataArray *cellIds  = boundary->GetCellData()->GetArray("TmpInternalIds");
  for (int i=0; i<boundary->GetNumberOfCells(); i++)
  {
    vtkIdType npts, *pts;
    boundary->GetCellPoints(i, npts, pts);
    if (npts != 2)
    {
      vtkErrorMacro("Should only have 2 points on line!");
      return 0;
    }
    vtkNew(vtkIdList, newCellPoints);
    int newId = pd->GetPoints()->InsertNextPoint(center);
    newPointData->CopyData(pd->GetPointData(), pointIds->GetTuple1(pts[0]), newId);
    newCellPoints->SetNumberOfIds(3);
    newCellPoints->SetId(0, newId);
    newCellPoints->SetId(1, pointIds->GetTuple1(pts[1]));
    newCellPoints->SetId(2, pointIds->GetTuple1(pts[0]));
    int newCellId = pd->GetNumberOfCells();
    pd->InsertNextCell(VTK_TRIANGLE, newCellPoints);
    newCellData->CopyData(pd->GetCellData(), cellIds->GetTuple1(i), newCellId);
  }

  pd->GetPointData()->PassData(newPointData);
  pd->GetCellData()->PassData(newCellData);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGroupsClipper::ComputeMassCenter(vtkPolyData *pd, double massCenter[3])
{
  massCenter[0] = 0.0;
  massCenter[1] = 0.0;
  massCenter[2] = 0.0;
  vtkNew(vtkCenterOfMass, centerFinder);
  centerFinder->SetInputData(pd);
  centerFinder->Update();
  centerFinder->GetCenter(massCenter);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGroupsClipper::GetClosestPointConnectedRegion(vtkPolyData *inPd,
                                                                         double origin[3],
                                                                         vtkPolyData *outPd)
{
  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(inPd);
  connector->SetExtractionModeToClosestPointRegion();
  connector->SetClosestPoint(origin);
  connector->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  outPd->DeepCopy(surfacer->GetOutput());

  return 1;
}
