/*=========================================================================

Program:   VMTK
Module:    $RCSfile: vtksvPolyDataCenterlineGroupsClipper.cxx,v $
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

#include "vtksvPolyDataCenterlineGroupsClipper.h"
#include "vtkAppendPolyData.h"
#include "vtkExecutive.h"
#include "vtkCellLocator.h"
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
#include "vtksvPolyBallLine.h"
#include "vtkMath.h"
#include "vtkMergeCells.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkTriangleFilter.h"
#include "vtkThreshold.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVersion.h"
#include "vtkXMLPolyDataWriter.h"

#ifndef VTK_SV_DOUBLE_TOL
#define VTK_SV_DOUBLE_TOL 1.0E-12
#endif

#ifndef VTK_SV_LARGE_DOUBLE
#define VTK_SV_LARGE_DOUBLE 1.0E+32
#endif

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkStandardNewMacro(vtksvPolyDataCenterlineGroupsClipper);

vtksvPolyDataCenterlineGroupsClipper::vtksvPolyDataCenterlineGroupsClipper()
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

vtksvPolyDataCenterlineGroupsClipper::~vtksvPolyDataCenterlineGroupsClipper()
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

vtkPolyData *vtksvPolyDataCenterlineGroupsClipper::GetClippedOutput()
{
  if (this->GetNumberOfOutputPorts() < 2)
    {
    return NULL;
    }

  return vtkPolyData::SafeDownCast(this->GetExecutive()->GetOutputData(1));
}

int vtksvPolyDataCenterlineGroupsClipper::RequestData(
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

  vtkSmartPointer<vtksvPolyBallLine> groupTubes =
    vtkSmartPointer<vtksvPolyBallLine>::New();
  groupTubes->SetInput(this->Centerlines);
  groupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  groupTubes->SetUseRadiusInformation(this->UseRadiusInformation);

  vtkSmartPointer<vtksvPolyBallLine> nonGroupTubes =
    vtkSmartPointer<vtksvPolyBallLine>::New();
  nonGroupTubes->SetInput(this->Centerlines);
  nonGroupTubes->SetPolyBallRadiusArrayName(this->CenterlineRadiusArrayName);
  nonGroupTubes->SetUseRadiusInformation(this->UseRadiusInformation);


  double point[3];
  double groupTubeValue, nonGroupTubeValue, tubeDifferenceValue;

  vtkIdType groupId;

  vtkSmartPointer<vtkIdList> centerlineGroupIds =
    vtkSmartPointer<vtkIdList>::New();

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
  vtkSmartPointer<vtkAppendPolyData> appendBranches =
    vtkSmartPointer<vtkAppendPolyData>::New();
  int numberOfPoints = clippingInput->GetNumberOfPoints();

  // Set up clipping arrays
  for (i=0; i<centerlineGroupIds->GetNumberOfIds(); i++)
    {
    groupId = centerlineGroupIds->GetId(i);

    std::stringstream groupstr;
    groupstr << groupId;
    std::string clipName = "ClippingArray_"+groupstr.str();

    vtkSmartPointer<vtkDoubleArray> clippingArray =
      vtkSmartPointer<vtkDoubleArray>::New();
    clippingArray->SetNumberOfComponents(1);
    clippingArray->SetNumberOfTuples(numberOfPoints);
    clippingArray->FillComponent(0,0.0);
    clippingArray->SetName(clipName.c_str());

    clippingInput->GetPointData()->AddArray(clippingArray);

    vtkSmartPointer<vtkIdList> groupTubesGroupIds =
      vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> nonGroupTubesGroupIds =
      vtkSmartPointer<vtkIdList>::New();
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

    vtkSmartPointer<vtkClipPolyData> clipper =
      vtkSmartPointer<vtkClipPolyData>::New();
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

    vtkSmartPointer<vtkPolyData> clippedBranch =
      vtkSmartPointer<vtkPolyData>::New();
    clippedBranch->DeepCopy(clipper->GetOutput());

    vtkSmartPointer<vtkPolyData> clippedOutputBranch =
      vtkSmartPointer<vtkPolyData>::New();
    clippedOutputBranch->DeepCopy(clipper->GetClippedOutput());

    vtkSmartPointer<vtkCleanPolyData> cleaner =
      vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->SetInputData(clippedBranch);
    cleaner->Update();
    vtkSmartPointer<vtkTriangleFilter> triangulator =
      vtkSmartPointer<vtkTriangleFilter>::New();
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
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
      vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetInputData(clippedBranch);
    writer->SetFileName(filename.c_str());
    writer->Write();

    vtkSmartPointer<vtkCleanPolyData> cleanerClipped =
      vtkSmartPointer<vtkCleanPolyData>::New();
    cleanerClipped->SetInputData(clippedOutputBranch);
    cleanerClipped->Update();
    vtkSmartPointer<vtkTriangleFilter> triangulatorClipped =
      vtkSmartPointer<vtkTriangleFilter>::New();
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

  vtkNew(vtkIdList, newPointIds);
  vtkNew(vtkDoubleArray, averageDistances);
  this->SplitGroups(finisher, separateIds, newPointIds, averageDistances);

  this->FillGroups(finisher, newPointIds, averageDistances);


  output->DeepCopy(finisher);

  return 1;
}

void vtksvPolyDataCenterlineGroupsClipper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtksvPolyDataCenterlineGroupsClipper::ReplaceDataOnCells(vtkPointSet *pointset,
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
int vtksvPolyDataCenterlineGroupsClipper::FindGroupSeparatingPoints(vtkPolyData *pd,
                                                                      vtkIdList *separateIds)
{
  pd->BuildLinks();
  int numPoints = pd->GetNumberOfPoints();
  vtkIntArray *groupIdsArray =
    vtkIntArray::SafeDownCast(pd->GetCellData()->GetArray(this->GroupIdsArrayName));

  for (int i=0; i<numPoints; i++)
  {
    vtkSmartPointer<vtkIdList> groupIds =
      vtkSmartPointer<vtkIdList>::New();
    this->GetPointGroups(pd, this->GroupIdsArrayName, i, groupIds);
    int pointType = groupIds->GetNumberOfIds();
    if (pointType == 3)
    {
      separateIds->InsertNextId(i);
      fprintf(stdout,"Holy toledo: %d\n", i);
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
int vtksvPolyDataCenterlineGroupsClipper::GetPointGroups(vtkPolyData *pd, std::string arrayName,
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
int vtksvPolyDataCenterlineGroupsClipper::SplitGroups(vtkPolyData *pd, vtkIdList *separateIds, vtkIdList *newPointIds,
                                                        vtkDoubleArray *averageDistances)
{
  fprintf(stdout,"Splitting\n");
  int numIds = separateIds->GetNumberOfIds();
  vtkSmartPointer<vtkIntArray> used =
    vtkSmartPointer<vtkIntArray>::New();
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

  separateIds->Reset();
  separateIds->SetNumberOfIds(0);
  for (int i=0; i<separatePoints->GetNumberOfPoints(); i++)
  {
    if (used->GetValue(i) == -1)
    {
      fprintf(stdout,"Found uno\n");
      double pt[3];
      separatePoints->GetPoint(i, pt);
      vtkNew(vtkIdList, closePointIds);
      separatePointLocator->FindClosestNPoints(3, pt, closePointIds);
      fprintf(stdout,"What the %lld\n", closePointIds->GetNumberOfIds());
      double avgDistance = 0.0;
      for (int j=0; j<closePointIds->GetNumberOfIds(); j++)
      {
        used->SetValue(closePointIds->GetId(j), 1);
        double pt0[3], pt1[3];
        separatePoints->GetPoint(closePointIds->GetId(j), pt0);
        separatePoints->GetPoint(closePointIds->GetId((j+1)%3), pt1);
        avgDistance += std::sqrt(std::pow(pt0[0] - pt1[0], 2.0) +
                                 std::pow(pt0[1] - pt1[1], 2.0) +
                                 std::pow(pt0[2] - pt1[2], 2.0));
      }
      avgDistance = avgDistance/3.0;
      averageDistances->InsertNextValue(avgDistance);
      vtkNew(vtkIdList, groupPoints);
      this->AddNewCriticalPoint(pd, separatePoints, closePointIds, newPointIds);
    }
  }
  pd->RemoveDeletedCells();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtksvPolyDataCenterlineGroupsClipper::AddNewCriticalPoint(vtkPolyData *pd,
                                                                vtkPoints *separatePoints,
                                                                vtkIdList *groupIds,
                                                                vtkIdList *newPointIds)
{
  vtkNew(vtkPointData, newPointData);
  vtkNew(vtkCellData, newCellData);
  newPointData->CopyAllocate(pd->GetPointData(), pd->GetNumberOfPoints() + 1);
  newCellData->CopyAllocate(pd->GetCellData(), pd->GetNumberOfCells() + 3);
  for (int i=0; i<pd->GetNumberOfPoints(); i++)
  {
    newPointData->CopyData(pd->GetPointData(), i, i);
  }
  for (int i=0; i<pd->GetNumberOfCells(); i++)
  {
    newCellData->CopyData(pd->GetCellData(), i, i);
  }

  double centroid[3];
  centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
  for (int i=0; i< groupIds->GetNumberOfIds(); i++)
  {
    double pt[3];
    separatePoints->GetPoint(groupIds->GetId(i), pt);
    for (int j=0; j<3; j++)
    {
      centroid[j] += pt[j];
    }
  }
  vtkMath::MultiplyScalar(centroid, 1.0/groupIds->GetNumberOfIds());
      fprintf(stdout,"Centroid %.4f %.4f %.4f\n", centroid[0], centroid[1], centroid[2]);

  vtkNew(vtkCellLocator, cellLocator);
  cellLocator->SetDataSet(pd);
  cellLocator->BuildLocator();

  double closestPt[3];
  vtkIdType closestCell;
  int subId;
  double distance;
  vtkSmartPointer<vtkGenericCell> genericCell =
    vtkSmartPointer<vtkGenericCell>::New();
  cellLocator->FindClosestPoint(centroid, closestPt, genericCell, closestCell, subId,
    distance);

  if (int(pd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(closestCell)) != -1)
  {
    int newCell;
    this->GetNewPointSpecial(pd, centroid, newCell, closestPt);
    closestCell = newCell;
  }
  fprintf(stdout,"Inserting Point %.4f %.4f %.4f\n", closestPt[0], closestPt[1], closestPt[2]);
  fprintf(stdout,"For cell %lld\n", closestCell);
  int newId = pd->GetPoints()->InsertNextPoint(closestPt);
  newPointIds->InsertNextId(newId);
  newPointData->CopyData(pd->GetPointData(), groupIds->GetId(0), newId);
  vtkIdType npts, *pts;
  pd->GetCellPoints(closestCell, npts, pts);
  for (int i=0; i<npts; i++)
  {
    int ptId0 = pts[i];
    int ptId1 = pts[(i+1)%npts];
    vtkNew(vtkIdList, newCellPoints);
    newCellPoints->SetNumberOfIds(3);
    newCellPoints->SetId(0, newId);
    newCellPoints->SetId(1, ptId0);
    newCellPoints->SetId(2, ptId1);
    int newCellId = pd->InsertNextCell(VTK_TRIANGLE, newCellPoints);
    newCellData->CopyData(pd->GetCellData(), closestCell, newCellId);
  }
  pd->DeleteCell(closestCell);

  pd->GetPointData()->PassData(newPointData);
  pd->GetCellData()->PassData(newCellData);
  return 1;
}

int vtksvPolyDataCenterlineGroupsClipper::GetNewPointSpecial(vtkPolyData *pd,
                                                               double centroid[3],
                                                               int &newCell,
                                                               double newPt[3])
{
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(pd);
  ider->SetIdsArrayName("TmpInternalIds");
  ider->Update();
  vtkNew(vtkPolyData, idOutPd);
  idOutPd->DeepCopy(ider->GetOutput());
  vtkNew(vtkPolyData, tmpPd);
  this->ThresholdPd(idOutPd, -1, -1, 1, this->GroupIdsArrayName, tmpPd);
  tmpPd->BuildLinks();

  vtkNew(vtkCellLocator, locator);
  locator->SetDataSet(tmpPd);
  locator->BuildLocator();

  double closestPt[3];
  vtkIdType closestCell;
  int subId;
  double distance;
  vtkSmartPointer<vtkGenericCell> genericCell =
    vtkSmartPointer<vtkGenericCell>::New();
  locator->FindClosestPoint(centroid, closestPt, genericCell, closestCell, subId,
    distance);

  newPt[0] = 0.0; newPt[1] = 0.0; newPt[2] = 0.0;
  vtkIdType npts, *pts;
  tmpPd->GetCellPoints(closestCell, npts, pts);
  for (int i=0; i< npts; i++)
  {
    double pt[3];
    tmpPd->GetPoint(pts[i], pt);
    for (int j=0; j<3; j++)
    {
      newPt[j] += pt[j];
    }
  }
  vtkMath::MultiplyScalar(newPt, 1.0/npts);
  fprintf(stdout,"New centroid %.4f %.4f %.4f\n", newPt[0], newPt[1], newPt[2]);
  newCell = tmpPd->GetCellData()->GetArray("TmpInternalIds")->GetTuple1(closestCell);
  return 1;
}

int vtksvPolyDataCenterlineGroupsClipper::ThresholdPd(vtkPolyData *pd, int minVal,
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


int vtksvPolyDataCenterlineGroupsClipper::FillGroups(vtkPolyData *pd,
                                                       vtkIdList *newPointIds,
                                                       vtkDoubleArray *averageDistances)
{
  pd->BuildLinks();
  for (int i=0; i<newPointIds->GetNumberOfIds(); i++)
  {
    int pointId = newPointIds->GetId(i);
    double avgDistance = averageDistances->GetValue(i);
    fprintf(stdout,"Doing the point %d\n", pointId);
    this->FillRegionGroups(pd, pointId, avgDistance);
  }
  return 1;
}

int vtksvPolyDataCenterlineGroupsClipper::FillRegionGroups(vtkPolyData *pd,
                                                             const int pointId,
                                                             const double rayDist)
{
  vtkIdList *cellListArray[3];
  vtkIdList *groupListArray[3];
  for (int i=0; i<3; i++)
  {
    cellListArray[i] = vtkIdList::New();
    groupListArray[i] = vtkIdList::New();
  }

  vtkNew(vtkIdList, cellIdList);
  pd->GetPointCells(pointId, cellIdList);
  if (cellIdList->GetNumberOfIds() != 3)
  {
    vtkErrorMacro("New point should be connected to three cells");
    for (int i=0; i<3; i++)
    {
      cellListArray[i]->Delete();
      groupListArray[i]->Delete();
    }
    return 0;
  }
  for (int i=0; i<3; i++)
  {
    int cellId = cellIdList->GetId(i);
    fprintf(stdout,"Start cell id %d\n", cellId);
    this->FillCellGroups(pd, cellId, pointId, cellListArray[i], groupListArray[i]);
    vtkIdType npts, *pts;
    pd->GetCellPoints(cellId, npts, pts);
    int groupId;
    for (int j=0; j<npts; j++)
    {
      if (pts[j] == pointId)
      {
        int pointId0 = pts[(j+1)%npts];
        int pointId1 = pts[(j+2)%npts];
        int pointId2 = pts[j];
        vtkNew(vtkIdList, cellNeighbors);
        pd->GetCellEdgeNeighbors(cellId, pointId0, pointId1, cellNeighbors);
        if (cellNeighbors->GetNumberOfIds() != 1)
        {
          vtkErrorMacro("There is a problem here");
          return 0;
        }
        double pt0[3], pt1[3], pt2[3];
        pd->GetPoint(pointId0, pt0);
        pd->GetPoint(pointId1, pt1);
        pd->GetPoint(pointId2, pt2);
        double midway[3];
        vtkMath::Add(pt0, pt1, midway);
        vtkMath::MultiplyScalar(midway, 0.5);
        double vec0[3];
        vtkMath::Subtract(midway, pt2, vec0);
        vtkMath::Normalize(vec0);
        vtkMath::MultiplyScalar(vec0, 1.5*rayDist);
        double finalPt[3];
        vtkMath::Add(pt2, vec0, finalPt);
        vtkNew(vtkCellLocator, locator);
        locator->SetDataSet(pd);
        locator->BuildLocator();
        double closestPt[3];
        vtkIdType closestCell;
        int subId;
        double distance;
        vtkSmartPointer<vtkGenericCell> genericCell =
          vtkSmartPointer<vtkGenericCell>::New();
        locator->FindClosestPoint(finalPt, closestPt, genericCell, closestCell, subId,
          distance);
        groupId = pd->GetCellData()->GetArray(this->GroupIdsArrayName)->GetTuple1(closestCell);
      }
    }
    for (int j=0; j<cellListArray[i]->GetNumberOfIds(); j++)
    {
      pd->GetCellData()->GetArray(this->GroupIdsArrayName)->
        SetTuple1(cellListArray[i]->GetId(j), groupId);
    }
  }

  for (int i=0; i<3; i++)
  {
    fprintf(stdout,"Array %d has %lld ids and %lld groups\n", i, cellListArray[i]->GetNumberOfIds(), groupListArray[i]->GetNumberOfIds());
    cellListArray[i]->Delete();
    groupListArray[i]->Delete();
  }
  return 1;
}

int vtksvPolyDataCenterlineGroupsClipper::FillCellGroups(vtkPolyData *pd,
                                                           const int cellId,
                                                           const int pointId,
                                                           vtkIdList *cellList,
                                                           vtkIdList *groupList)
{

  cellList->InsertUniqueId(cellId);
  vtkIdType npts, *pts;
  pd->GetCellPoints(cellId, npts, pts);
  for (int j=0; j<3; j++)
  {
    if (pts[j] == pointId)
    {
      int pointId0 = pts[(j+1)%npts];
      int pointId1 = pts[(j+2)%npts];
      vtkNew(vtkIdList, cellNeighbors);
      pd->GetCellEdgeNeighbors(cellId, pointId0, pointId1, cellNeighbors);
      if (cellNeighbors->GetNumberOfIds() != 1)
      {
        vtkErrorMacro("There is a problem here");
        return 0;
      }
      int groupId = pd->GetCellData()->GetArray(this->GroupIdsArrayName)->
        GetTuple1(cellNeighbors->GetId(0));
      fprintf(stdout,"What is group id %d\n", groupId);
      if (groupId != -1)
      {
        groupList->InsertUniqueId(groupId);
      }
      else
      {
        if (cellList->IsId(cellNeighbors->GetId(0)) != -1)
        {
          fprintf(stdout,"What the eff returning %d\n", cellId);
          return 1;
        }
        fprintf(stdout,"Doing cell %lld and points %d %d\n", cellNeighbors->GetId(0), pointId0, pointId1);
        this->FillCellGroups(pd, cellNeighbors->GetId(0), pointId0, cellList, groupList);
        this->FillCellGroups(pd, cellNeighbors->GetId(0), pointId1, cellList, groupList);
      }
    }
  }

  return 1;
}
