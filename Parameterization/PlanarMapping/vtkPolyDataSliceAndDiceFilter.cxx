/*=========================================================================
 *
 * Copyright (c) 2014 The Regents of the University of California.
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

/** @file vtkPolyDataSliceAndDiceFilter.cxx
 *  @brief This implements the vtkPolyDataSliceAndDiceFilter filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkPolyDataSliceAndDiceFilter.h"

#include "vtkAppendPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkDoubleArray.h"
#include "vtkExtractGeometry.h"
#include "vtkFeatureEdges.h"
#include "vtkFindGeodesicPath.h"
#include "vtkFindSeparateRegions.h"
#include "vtkFloatArray.h"
#include "vtkGetBoundaryFaces.h"
#include "vtkGradientFilter.h"
#include "vtkIdFilter.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointDataToCellData.h"
#include "vtkPointLocator.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"
#include "vtkTextureMapToSphere.h"
#include "vtkThreshold.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <iostream>
#include <sstream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkPolyDataSliceAndDiceFilter, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkPolyDataSliceAndDiceFilter);

//---------------------------------------------------------------------------
vtkPolyDataSliceAndDiceFilter::vtkPolyDataSliceAndDiceFilter()
{
  this->InitialPd       = vtkPolyData::New();
  this->WorkPd          = vtkPolyData::New();
  this->GraphPd         = vtkPolyData::New();
  this->SurgeryLines    = vtkPolyData::New();
  this->Centerlines     = NULL;
  this->Polycube        = NULL;
  this->CenterlineGraph = NULL;

  this->StartPtIds      = vtkIntArray::New();
  this->TopNormals      = vtkDoubleArray::New();
  this->RightNormals    = vtkDoubleArray::New();
  this->BifTopNormals   = vtkDoubleArray::New();
  this->BifRightNormals = vtkDoubleArray::New();
  this->NumberOfCubes   = vtkIntArray::New();
  this->StartPtIds->SetNumberOfComponents(4);
  this->RightNormals->SetNumberOfComponents(3);
  this->TopNormals->SetNumberOfComponents(3);
  this->NumberOfCubes->SetNumberOfComponents(3);

  this->BoundaryPointsArrayName = NULL;
  this->GroupIdsArrayName = NULL;
  this->SliceIdsArrayName = NULL;
  this->SphereRadiusArrayName = NULL;
  this->InternalIdsArrayName = NULL;
  this->DijkstraArrayName = NULL;

  this->ConstructPolycube = 0;
  this->SliceDirection = 0;
  this->TotalSliceId = 0;
  this->MaxGroupNumber = 0;
  this->MaxSegmentNumber = 0;
  this->SliceLength = 1.5;
}

//---------------------------------------------------------------------------
vtkPolyDataSliceAndDiceFilter::~vtkPolyDataSliceAndDiceFilter()
{
  if (this->InitialPd != NULL)
  {
    this->InitialPd->Delete();
    this->InitialPd = NULL;
  }
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->GraphPd != NULL)
  {
    this->GraphPd->Delete();
    this->GraphPd = NULL;
  }
  if (this->Centerlines != NULL)
  {
    this->Centerlines->Delete();
    this->Centerlines = NULL;
  }
  if (this->SurgeryLines != NULL)
  {
    this->SurgeryLines->Delete();
    this->SurgeryLines = NULL;
  }
  if (this->Polycube != NULL)
  {
    this->Polycube->Delete();
    this->Polycube = NULL;
  }
  if (this->CenterlineGraph != NULL)
  {
    delete this->CenterlineGraph;
    this->CenterlineGraph = NULL;
  }
  if (this->StartPtIds != NULL)
  {
    this->StartPtIds->Delete();
    this->StartPtIds = NULL;
  }
  if (this->TopNormals != NULL)
  {
    this->TopNormals->Delete();
    this->TopNormals = NULL;
  }
  if (this->RightNormals != NULL)
  {
    this->RightNormals->Delete();
    this->RightNormals = NULL;
  }
  if (this->BifTopNormals != NULL)
  {
    this->BifTopNormals->Delete();
    this->BifTopNormals = NULL;
  }
  if (this->BifRightNormals != NULL)
  {
    this->BifRightNormals->Delete();
    this->BifRightNormals = NULL;
  }
  if (this->NumberOfCubes != NULL)
  {
    this->NumberOfCubes->Delete();
    this->NumberOfCubes = NULL;
  }

  if (this->BoundaryPointsArrayName != NULL)
  {
    delete [] this->BoundaryPointsArrayName;
    this->BoundaryPointsArrayName = NULL;
  }
  if (this->GroupIdsArrayName != NULL)
  {
    delete [] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
  }
  if (this->SliceIdsArrayName != NULL)
  {
    delete [] this->SliceIdsArrayName;
    this->SliceIdsArrayName = NULL;
  }
  if (this->SphereRadiusArrayName != NULL)
  {
    delete [] this->SphereRadiusArrayName;
    this->SphereRadiusArrayName = NULL;
  }
  if (this->InternalIdsArrayName != NULL)
  {
    delete [] this->InternalIdsArrayName;
    this->InternalIdsArrayName = NULL;
  }
  if (this->DijkstraArrayName != NULL)
  {
    delete [] this->DijkstraArrayName;
    this->DijkstraArrayName = NULL;
  }
}

//---------------------------------------------------------------------------
void vtkPolyDataSliceAndDiceFilter::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkPolyDataSliceAndDiceFilter::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->InitialPd->DeepCopy(input1);

  this->WorkPd->DeepCopy(this->InitialPd);
  if (this->Centerlines == NULL)
  {
    this->Centerlines = vtkPolyData::New();
    this->ComputeCenterlines();
    this->ExtractBranches();
  }

  if (this->PreProcessPolyData() != 1)
  {
    vtkErrorMacro("Error in preprocessing the polydata\n");
    return 0;
  }
  vtkNew(vtkXMLPolyDataWriter, pdwriter);
  pdwriter->SetInputData(this->GraphPd);
  pdwriter->SetFileName("/Users/adamupdegrove/Desktop/tmp/graph.vtp");
  pdwriter->Write();
  if (this->SliceBifurcations() != 1)
  {
    vtkErrorMacro("Error in slicing the polydata\n");
    return 0;
  }
  vtkNew(vtkXMLPolyDataWriter, nwriter);
  nwriter->SetInputData(this->WorkPd);
  nwriter->SetFileName("/Users/adamupdegrove/Desktop/tmp/whatevs.vtp");
  nwriter->Write();
  if (this->SliceBranches() != 1)
  {
    vtkErrorMacro("Error in slicing the polydata\n");
    return 0;
  }
  vtkNew(vtkXMLPolyDataWriter, lwriter);
  lwriter->SetInputData(this->SurgeryLines);
  lwriter->SetFileName("/Users/adamupdegrove/Desktop/tmp/surglines.vtp");
  lwriter->Write();
  if (this->ConstructPolycube)
  {
    this->Polycube = vtkGeneralizedPolycube::New();
    if (this->BuildPolycube() != 1)
    {
      vtkErrorMacro("Error in constructing polycube\n");
      return 0;
    }
  }
  //vtkNew(vtkXMLUnstructuredGridWriter, writer);
  //writer->SetInputData(this->Polycube);
  //writer->SetFileName("/Users/adamupdegrove/Desktop/tmp/polycube.vtu");
  //writer->Write();


  output->DeepCopy(this->WorkPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::ComputeCenterlines()
{
  return 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::ExtractBranches()
{
  return 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::PreProcessPolyData()
{
  if (this->FindGroupBoundaries() != 1)
  {
    vtkErrorMacro("Unable to find boundaries of input group ids data array");
    return 0;
  }

  if (this->GetCriticalPoints() != 1)
  {
    vtkErrorMacro("Unable to retrieve critical points");
    return 0;
  }

  this->FormDirectionTable(this->DirectionTable);
  vtkDataArray *groupIds = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName);
  double minmax[2];
  groupIds->GetRange(minmax);
  this->MaxGroupNumber = int(minmax[1]);

  vtkNew(vtkPoints, skeletonPoints);
  vtkNew(vtkCellArray, skeletonCells);
  vtkNew(vtkIntArray, skeletonGroupIds);
  this->CenterlineGraph = new svGraph(0, this->Centerlines,
                                     this->GroupIdsArrayName,
                                     this->CriticalPointMap,
                                     this->DirectionTable);
  if (this->CenterlineGraph->BuildGraph() != 1)
  {
    vtkErrorMacro("Unable to form skeleton of polydata");
    return 0;
  }
  this->CenterlineGraph->GetGraphPolyData(this->GraphPd);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::SetDir(const int dir, double newDir[3])
{
  if (dir == RIGHT)
  {
    newDir[0] = 1.0; newDir[1] = 0.0; newDir[2] = 0.0;
  }
  else if (dir == LEFT)
  {
    newDir[0] = -1.0; newDir[1] = 0.0; newDir[2] = 0.0;
  }
  else if (dir == FRONT)
  {
    newDir[0] = 0.0; newDir[1] = 1.0; newDir[2] = 0.0;
  }
  else if (dir == BACK)
  {
    newDir[0] = 0.0; newDir[1] = -1.0; newDir[2] = 0.0;
  }
  else if (dir == UP)
  {
    newDir[0] = 0.0; newDir[1] = 0.0; newDir[2] = 1.0;
  }
  else if (dir == DOWN)
  {
    newDir[0] = 0.0; newDir[1] = 0.0; newDir[2] = -1.0;
  }
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::FormDirectionTable(int dirTable[6][4])
{
  //RIGHT: UP    BACK DOWN  FRONT
  //LEFT:  DOWN  BACK UP    FRONT
  //FRONT: RIGHT DOWN LEFT  UP
  //BACK:  RIGHT UP   LEFT  DOWN
  //UP:    LEFT  BACK RIGHT FRONT
  //DOWN:  RIGHT BACK LEFT  FRONT

  dirTable[RIGHT][0] = UP;    dirTable[RIGHT][1] = BACK; dirTable[RIGHT][2] = DOWN;  dirTable[RIGHT][3] = FRONT;
  dirTable[LEFT][0]  = DOWN;  dirTable[LEFT][1]  = BACK; dirTable[LEFT][2]  = UP;    dirTable[LEFT][3]  = FRONT;
  dirTable[FRONT][0] = RIGHT; dirTable[FRONT][1] = DOWN; dirTable[FRONT][2] = LEFT;  dirTable[FRONT][3] = UP;
  dirTable[BACK][0]  = RIGHT; dirTable[BACK][1]  = UP;   dirTable[BACK][2]  = LEFT;  dirTable[BACK][3]  = DOWN;
  dirTable[UP][0]    = LEFT;  dirTable[UP][1]    = BACK; dirTable[UP][2]    = RIGHT; dirTable[UP][3]    = FRONT;
  dirTable[DOWN][0]  = RIGHT; dirTable[DOWN][1]  = BACK; dirTable[DOWN][2]  = LEFT;  dirTable[DOWN][3]  = FRONT;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetUniqueNeighbors(std::multimap<int, int> &map, const int key, std::list<int> &keyVals,
                                                      std::list<int> &uniqueKeys)
{
  int numVals = keyVals.size();

  std::list<int>::iterator valit = keyVals.begin();
  for (int i=0; valit != keyVals.end(); ++valit)
  {
    std::list<int> valKeys;
    vtkPolyDataSliceAndDiceFilter::GetKeysFromMap(map, *valit, valKeys);
    std::list<int>::iterator keyit = valKeys.begin();
    for (int j=0; keyit != valKeys.end(); ++keyit)
    {
      if (*keyit != key)
      {
        uniqueKeys.push_back(*keyit);
      }
    }
  }
  uniqueKeys.sort();
  uniqueKeys.unique();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::FindGroupBoundaries()
{
  vtkNew(vtkFindSeparateRegions, separator);
  separator->SetInputData(this->WorkPd);
  separator->SetOutPointArrayName(this->BoundaryPointsArrayName);
  separator->SetArrayName(this->GroupIdsArrayName);
  separator->Update();

  this->WorkPd->DeepCopy(separator->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetCriticalPoints()
{
  //Forms a map of gropus to critical points
  //Keys: Group ids on surface
  //Values: Ids of points on surface that touch key group
  this->WorkPd->BuildLinks();
  int numPoints = this->WorkPd->GetNumberOfPoints();
  vtkDataArray *boundaryPointArray =
    this->WorkPd->GetPointData()->GetArray(this->BoundaryPointsArrayName);

  for (int i=0; i<numPoints; i++)
  {
    int isBoundary = boundaryPointArray->GetTuple1(i);
    vtkNew(vtkIdList, groupIds);
    if (isBoundary)
    {
      this->GetPointGroups(this->WorkPd, this->GroupIdsArrayName, i, groupIds);
      int pointType = groupIds->GetNumberOfIds();
      if (pointType < 2)
      {
        vtkErrorMacro("Point incorrectly identified");
      }
      if (pointType == 3)
      {
        this->InsertCriticalPoints(i, groupIds);
      }
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
int vtkPolyDataSliceAndDiceFilter::GetPointGroups(vtkPolyData *pd, std::string arrayName,
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
int vtkPolyDataSliceAndDiceFilter::InsertCriticalPoints(const int pointId, vtkIdList *groupIds)
{
  int numIds = groupIds->GetNumberOfIds();

  for (int i=0; i<numIds; i++)
  {
    int groupId = groupIds->GetId(i);
    this->CriticalPointMap.insert(std::make_pair(groupId, pointId));
    //fprintf(stdout,"Added critical point pair: %d %d\n", groupId, pointId);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetBranch(const int branchId, vtkPolyData *branchPd,
                                             vtkPolyData *branchCenterlines)
{
  this->ThresholdPd(this->WorkPd, branchId, branchId, 1,
    this->SegmentIdsArrayName, branchPd);

  if (branchPd->GetNumberOfPoints() != 0)
  {
    vtkNew(vtkPolyData, centerlineBranchPd);
    this->ThresholdPd(this->Centerlines, branchId, branchId, 1,
      this->GroupIdsArrayName, centerlineBranchPd);

    //Need to get just first cell of centerlines. There are duplicate for each centerline running through
    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(centerlineBranchPd);
    ider->SetIdsArrayName(this->InternalIdsArrayName);
    ider->Update();
    vtkNew(vtkPolyData, tmpPd);
    tmpPd->ShallowCopy(ider->GetOutput());

    this->ThresholdPd(tmpPd, 0, 0, 1,
      this->InternalIdsArrayName, branchCenterlines);

    branchCenterlines->GetCellData()->RemoveArray(this->InternalIdsArrayName);
    branchCenterlines->GetPointData()->RemoveArray(this->InternalIdsArrayName);
  }

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetAllMapKeys(std::multimap<int, int> &map,
                                                 std::list<int> &list)
{
  std::multimap<int, int>::iterator it = map.begin();

  for (int i=0; it != map.end(); ++it)
  {
    list.push_back(it->first);
  }
  list.unique();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetAllMapValues(std::multimap<int, int> &map,
                                                   std::list<int> &list)
{
  std::multimap<int, int>::iterator it = map.begin();

  for (int i=0; it != map.end(); ++it)
  {
    list.push_back(it->second);
  }
  list.unique();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetValuesFromMap(std::multimap<int, int> &map,
                                                    const int key,
                                                    std::list<int> &list)
{
  std::multimap<int, int>::iterator it = map.begin();

  for (int i=0; it != map.end(); ++it)
  {
    if (it->first == key)
    {
      list.push_back(it->second);
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
int vtkPolyDataSliceAndDiceFilter::GetKeysFromMap(std::multimap<int, int> &map,
                                                  const int value,
                                                  std::list<int> &list)
{
  std::multimap<int, int>::iterator it = map.begin();

  for (int i=0; it != map.end(); ++it)
  {
    if (it->second == value)
    {
      list.push_back(it->first);
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
int vtkPolyDataSliceAndDiceFilter::GetCommonValues(std::multimap<int, int> &map,
                                                   const int keyA, const int keyB,
                                                   std::list<int> &returnList)
{
  std::list<int> listA, listB;
  vtkPolyDataSliceAndDiceFilter::GetValuesFromMap(map, keyA, listA);
  vtkPolyDataSliceAndDiceFilter::GetValuesFromMap(map, keyB, listB);
  vtkPolyDataSliceAndDiceFilter::ListIntersection(listA, listB, returnList);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::ListIntersection(std::list<int> &listA,
                                                    std::list<int> &listB,
                                                    std::list<int> &returnList)
{
  std::set_intersection(listA.begin(), listA.end(),
                        listB.begin(), listB.end(),
                        std::inserter(returnList,returnList.begin()));

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::SliceBranches()
{
  int numIds = this->WorkPd->GetCellData()->
    GetArray(this->SegmentIdsArrayName)->GetNumberOfTuples();
  vtkNew(vtkIntArray, sliceIds);
  sliceIds->SetNumberOfTuples(numIds);
  sliceIds->SetName(this->SliceIdsArrayName);
  sliceIds->FillComponent(0, -1);
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(this->WorkPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  this->WorkPd->DeepCopy(ider->GetOutput());

  vtkNew(vtkPoints, surgeryPts);
  vtkNew(vtkCellArray, surgeryLines);
  vtkNew(vtkIntArray, surgeryData);
  vtkNew(vtkIntArray, surgeryCellData);
  for (int i=0; i<=this->MaxSegmentNumber; i++)
  {
    vtkNew(vtkPolyData, branchPd);
    vtkNew(vtkPolyData, branchCenterline);
    this->GetBranch(i, branchPd, branchCenterline);

    int numPoints = branchPd->GetNumberOfPoints();

    if (numPoints != 0)
    {
      if (i <= this->MaxGroupNumber)
      {
        this->SliceBranch(branchPd, branchCenterline, i, sliceIds, surgeryPts,
                          surgeryLines, surgeryData);
        surgeryCellData->InsertNextValue(i);
      }
      else
      {
        this->ReplaceDataOnCells(branchPd, sliceIds,
                                 this->TotalSliceId, -1, this->InternalIdsArrayName);
        this->TotalSliceId++;
      }
    }
  }
  this->SurgeryLines->SetPoints(surgeryPts);
  this->SurgeryLines->SetLines(surgeryLines);
  surgeryData->SetName(this->InternalIdsArrayName);
  surgeryCellData->SetName(this->SegmentIdsArrayName);
  this->SurgeryLines->GetPointData()->AddArray(surgeryData);
  this->SurgeryLines->GetCellData()->AddArray(surgeryCellData);

  this->WorkPd->GetCellData()->AddArray(sliceIds);

  this->WorkPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  this->WorkPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetSectionZAxis(const double endPt[3], const double startPt[3],
                                                   double zvec[3])
{
  //Get approximate z axis of section
  vtkMath::Subtract(startPt, endPt, zvec);
  vtkMath::Normalize(zvec);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetSectionXAxis(const double endPt[3], const double startPt[3],
                                                   const double surfacePt[3],
                                                   double xvec[3])
{
  double vec0[3], vec1[3];

  //Get approximate x axis of section
  vtkMath::Subtract(surfacePt, startPt, vec0);
  vtkMath::Subtract(endPt, startPt, vec1);
  vtkMath::Normalize(vec1);
  double vecmag = vtkMath::Dot(vec0, vec1);
  vtkMath::MultiplyScalar(vec1, vecmag);
  vtkMath::Subtract(vec0, vec1, xvec);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetCutPlane(const double endPt[3], const double startPt[3],
                                               const double length, double origin[3], vtkPlane *cutPlane)
{
  double normal[3];
  vtkMath::Subtract(endPt, startPt, normal);
  vtkMath::Normalize(normal);

  //vtkMath::MultiplyScalar(normal, length);
  //vtkMath::Add(startPt, normal, origin);
  //vtkMath::Normalize(normal);
  for (int i=0; i<3; i++)
  {
    origin[i] = endPt[i];
  }

  cutPlane->SetOrigin(origin);
  cutPlane->SetNormal(normal);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::ExtractionCut(vtkPolyData *inPd, vtkPlane *cutPlane,
                                                 const int extractBoundaryCells,
                                                 const int extractInside,
                                                 vtkPolyData *outPd)
{
    vtkNew(vtkExtractGeometry, cutter);
    cutter->SetInputData(inPd);
    cutter->SetImplicitFunction(cutPlane);
    cutter->SetExtractBoundaryCells(extractBoundaryCells);
    cutter->ExtractOnlyBoundaryCellsOff();
    cutter->SetExtractInside(extractInside);
    cutter->Update();

    vtkNew(vtkDataSetSurfaceFilter, surfacer);
    surfacer->SetInputData(cutter->GetOutput());
    surfacer->Update();

    outPd->DeepCopy(surfacer->GetOutput());

    return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetClosestPointConnectedRegion(vtkPolyData *inPd,
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

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::DetermineSliceStrategy(vtkPolyData *branchPd,
                                                          const int branchId,
                                                          vtkPolyData *branchCenterline,
                                                          int &branchStartPtId,
                                                          vtkIdList *surgeryPoints,
                                                          int &centerlineStartPtId)
{
  double topPt[3], bottomPt[3];
  branchCenterline->GetPoint(0, topPt);
  branchCenterline->GetPoint(branchCenterline->GetNumberOfPoints() - 1, bottomPt);

  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(branchPd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();

  vtkNew(vtkPolyData, topPoly);
  vtkNew(vtkPolyData, bottomPoly);
  this->GetClosestPointConnectedRegion(boundaries->GetOutput(), topPt, topPoly);
  this->GetClosestPointConnectedRegion(boundaries->GetOutput(), bottomPt, bottomPoly);

  vtkDataArray *topPointIds =
      topPoly->GetPointData()->GetArray(this->InternalIdsArrayName);
  vtkDataArray *bottomPointIds =
      bottomPoly->GetPointData()->GetArray(this->InternalIdsArrayName);

  std::list<int> surgeryIds;
  this->GetValuesFromMap(this->SurgeryPointMap, branchId, surgeryIds);
  int numSurgeryPoints = surgeryIds.size();
  this->SliceDirection = 0;
  centerlineStartPtId = 0;
  if (numSurgeryPoints == 0)
  {
    fprintf(stdout,"Strategy 0, any surgery points will do\n");
    //There are no critical points. Any point on the boundary can be used as
    //the start point. This also means this is a genus zero surface with no
    //branching locations. Needs to be revised for non vascular models where
    //there is not inlet an outlet.
    branchStartPtId = topPointIds->GetTuple1(0);
    fprintf(stdout,"Start point is: %d\n", branchStartPtId);
    double secondPt[3], sidePt[3], zvec[3], xvec[3];
    branchCenterline->GetPoint(1, secondPt);
    topPoly->GetPoint(0, sidePt);
    vtkMath::Subtract(topPt, secondPt, zvec);
    vtkMath::Normalize(zvec);
    vtkMath::Subtract(sidePt, topPt, xvec);
    vtkMath::Normalize(xvec);
    this->GetFirstSurgeryPoints(topPoly, 0, surgeryPoints, xvec, zvec);
  }
  else if (numSurgeryPoints == 4)
  {
    fprintf(stdout,"Strategy 1, must actively use surgery points from bifurcations segmentation\n");
    //There are two critical points. This means that these points should be used
    //as corners of a polycube along with two other points.
    int onTopBoundary = 0;
    int front;
    int back;
    for (int i=0; i<topPoly->GetNumberOfPoints(); i++)
    {
      int pointId = topPointIds->GetTuple1(i);
      if (pointId == surgeryIds.front())
      {
        onTopBoundary = 1;
      }
    }
    if (!onTopBoundary)
    {
      fprintf(stdout,"Centerline upside down, starting from bottom\n");
      this->SliceDirection = 1;
      centerlineStartPtId = branchCenterline->GetNumberOfPoints() - 1;
    }
    std::list<int>::iterator it = surgeryIds.begin();
    for (int i=0; it != surgeryIds.end(); it++)
    {
      surgeryPoints->InsertNextId(*it);
    }
    branchStartPtId = surgeryPoints->GetId(0);
  }
  else
  {
    vtkErrorMacro("This shouldnt be possible. Something went wrong");
    return 0;
  }

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetSurgeryPoints(vtkPolyData *pd,
                                                    vtkDataArray *pointIds,
                                                    const double clStartPt[3],
                                                    const double clSecondPt[3],
                                                    const int front,
                                                    const int back,
                                                    vtkIdList *surgeryPoints)
{
  //Make front and back always start with the lowest ids so that same ids are
  //retrieved for every boundary
  int pointId = front;
  vtkNew(vtkIdList, startCellIds);
  pd->GetPointCells(pointId, startCellIds);
  int prevCellId = startCellIds->GetId(0);
  int prevCellId2 = startCellIds->GetId(1);
  vtkIdType npts, *pts;
  pd->GetCellPoints(prevCellId, npts, pts);
  int secondPtId;
  if (pts[0] == pointId)
  {
    secondPtId = pts[1];
  }
  else
  {
    secondPtId = pts[0];
  }
  double pt0[3], pt1[3];
  pd->GetPoint(pointId, pt0);
  pd->GetPoint(secondPtId, pt1);
  double vec0[3], vec1[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  //Get the vector from start to end
  double zvec[3];
  this->GetSectionZAxis(clSecondPt, clStartPt, zvec);

  //Get starting point from tmp id
  double xvec[3];
  this->GetSectionXAxis(clSecondPt, clStartPt, pt0, xvec);

  vtkMath::Cross(xvec, vec0, vec1);

  if (vtkMath::Dot(zvec, vec1) > 0)
  {
    prevCellId = startCellIds->GetId(1);
    prevCellId2 = startCellIds->GetId(0);
  }

  vtkNew(vtkIdList, halfPoints);
  this->GetHalfSurgeryPoints(pd, pointIds, prevCellId, front, back, halfPoints);
  this->GetHalfSurgeryPoints(pd, pointIds, prevCellId2, front, back, halfPoints);
  surgeryPoints->InsertNextId(halfPoints->GetId(0));
  surgeryPoints->InsertNextId(halfPoints->GetId(1));
  surgeryPoints->InsertNextId(halfPoints->GetId(3));
  surgeryPoints->InsertNextId(halfPoints->GetId(2));

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetHalfSurgeryPoints(vtkPolyData *pd,
                                                        vtkDataArray *pointIds,
                                                        const int cellId,
                                                        const int front,
                                                        const int back,
                                                        vtkIdList *surgeryPoints)
{
  //Make front and back always start with the lowest ids so that same ids are
  //retrieved for every boundary
  pd->BuildLinks();
  int pointId = front;
  int prevCellId = cellId;

  //Getting length of the loop separated into the two sides of the boundary
  //by the critical points
  double length = 0.0;
  while (pointId != back)
  {
    double pt0[3], pt1[3];
    pd->GetPoint(pointId, pt0);
    this->IteratePoint(pd, pointId, prevCellId);
    pd->GetPoint(pointId, pt1);
    length += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                        std::pow(pt1[1] - pt0[1], 2.0) +
                        std::pow(pt1[2] - pt0[2], 2.0));
  }

  //Finding the surgery points which separate each half of the boundary into
  //three slices
  double surgeryCount = 1.0;
  double currLength = 0.0;;
  pointId = front;
  prevCellId = cellId;
  while (pointId != back)
  {
    double pt0[3], pt1[3];
    pd->GetPoint(pointId, pt0);
    this->IteratePoint(pd, pointId, prevCellId);
    pd->GetPoint(pointId, pt1);
    currLength += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                            std::pow(pt1[1] - pt0[1], 2.0) +
                            std::pow(pt1[2] - pt0[2], 2.0));
    if (currLength > length*(surgeryCount/3.0))
    {
      surgeryPoints->InsertNextId(pointIds->GetTuple1(pointId));
      surgeryCount += 1.0;
    }
  }

  return 1;
}

////---------------------------------------------------------------------------
///**
// * @brief
// * @param *pd
// * @return
// */
//int vtkPolyDataSliceAndDiceFilter::GetConsistentStartDirection(vtkPolyData *pd,
//                                                               vtkIdList *cellIds,
//                                                               const int startPointId,
//                                                               int &startCellId)
//{
//  vtkIdType npts, *pts;
//  pd->GetCellPoints(cellIds->GetId(0), npts, pts);
//  double pt0[3], pt1[3];
//  pd->GetPoint(pts[0],pt0);
//  pd->GetPoint(pts[1],pt1);
//  double vec0[3];
//  if (pts[0] == startPointId)
//  {
//    vtkMath::Subtract(pt1, pt0, vec0);
//  }
//  else
//  {
//    vtkMath::Subtract(pt0, pt1, vec0);
//  }
//  vtkMath::Normalize(vec0);
//  double vec1[3]; vec1[0] = 1.0; vec1[1] = 0.0; vec1[2] = 0.0;
//  double testDir = vtkMath::Dot(vec0, vec1);
//  if (testDir > 0)
//  {
//    startCellId = cellIds->GetId(0);
//  }
//  else
//  {
//    startCellId = cellIds->GetId(1);
//  }
//
//  return 1;
//}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::IteratePoint(vtkPolyData *pd, int &pointId, int &prevCellId)
{
  vtkNew(vtkIdList, ptCellIds);
  pd->GetPointCells(pointId, ptCellIds);
  int cellId;
  if (ptCellIds->GetId(0) == prevCellId)
  {
    cellId = ptCellIds->GetId(1);
  }
  else
  {
    cellId = ptCellIds->GetId(0);
  }
  prevCellId = cellId;

  vtkIdType npts, *pts;
  pd->GetCellPoints(prevCellId, npts, pts);
  int newId;
  if (pts[0] == pointId)
  {
    newId = pts[1];
  }
  else
  {
    newId = pts[0];
  }
  pointId = newId;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkPolyDataSliceAndDiceFilter::CheckLength(int &ptId, const int numPts,
                                                int &done)
{
  if (this->SliceDirection == 0)
  {
    if (ptId > numPts - 2)
    {
      ptId = numPts - 2;
      done = 1;
    }
  }
  else if (this->SliceDirection == 1)
  {
    if (ptId < 1)
    {
      ptId = 1;
      done = 1;
    }
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkPolyDataSliceAndDiceFilter::UpdatePtId(int &ptId)
{
  if (this->SliceDirection == 0)
  {
    ptId++;
  }
  else if (this->SliceDirection == 1)
  {
    ptId--;
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::SliceBranch(vtkPolyData *branchPd,
                                               vtkPolyData *branchCenterline,
                                               const int branchId,
                                               vtkDataArray *sliceIds,
                                               vtkPoints *surgeryPts,
                                               vtkCellArray *surgeryLines,
                                               vtkIntArray *surgeryData)
{
  fprintf(stdout,"Slicing branch: %d\n", branchId);
  int customStart = this->StartPtIds->GetNumberOfTuples();

  int numCells = branchPd->GetNumberOfCells();
  int numCenterlinePts = branchCenterline->GetNumberOfPoints();
    vtkDataArray *radiusArray =
      branchCenterline->GetPointData()->GetArray(this->SphereRadiusArrayName);

  vtkNew(vtkPolyData, currentPd); currentPd->DeepCopy(branchPd);
  vtkNew(vtkPolyData, leftovers); leftovers->DeepCopy(branchPd);

  double totalLength = 0.0;
  this->GetBranchLength(branchCenterline, totalLength);

  int linePtId = 0;
  int startPtId = -1;
  vtkNew(vtkIdList, surgeryPoints);
  this->DetermineSliceStrategy(branchPd, branchId, branchCenterline,
                               startPtId, surgeryPoints, linePtId);

  vtkNew(vtkPoints, testPts);
  vtkNew(vtkCellArray, testCells);
  vtkNew(vtkIdList, surgeryLineIds);
  int done = 0;
  int sliceId = 0;
  double currLength = 0.0;
  while (!done)
  {
    double inscribedRadius = radiusArray->GetTuple1(linePtId);
    double sliceLength = inscribedRadius;
    if (inscribedRadius > 0.5)
    {
      sliceLength *= this->SliceLength;
    }
    else
    {
      sliceLength *= this->SliceLength * (4.0/3.0);
    }
    if ((currLength + 2.0*this->SliceLength*inscribedRadius) >= totalLength)
    {
      sliceLength = 4.0*this->SliceLength*(totalLength - currLength)/9.0;
      done = 1;
    }

    double centerlineLength = 0;
    double startPt[3], pt0[3], pt1[3];
    branchCenterline->GetPoint(linePtId, startPt);

    while (centerlineLength < sliceLength)
    {
      this->CheckLength(linePtId, numCenterlinePts, done);
      branchCenterline->GetPoint(linePtId, pt0);
      this->UpdatePtId(linePtId);
      branchCenterline->GetPoint(linePtId, pt1);
      centerlineLength += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                                    std::pow(pt1[1] - pt0[1], 2.0) +
                                    std::pow(pt1[2] - pt0[2], 2.0));
    }

    currLength += centerlineLength;

    //Get the vector from start to end
    double zvec[3];
    this->GetSectionZAxis(pt1, startPt, zvec);

    //Get starting point from tmp id
    double xvec[3], surfacePt[3];
    vtkDataArray *pointIds = branchPd->GetPointData()->GetArray(this->InternalIdsArrayName);
    int pointId = pointIds->LookupValue(surgeryPoints->GetId(0));
    branchPd->GetPoint(pointId, surfacePt);
    this->GetSectionXAxis(pt1, startPt, surfacePt, xvec);

    //Get the cut plane
    double origin[3];
    vtkNew(vtkPlane, cutPlane);
    this->GetCutPlane(pt1, pt0, centerlineLength, origin, cutPlane);

    //Cut the pds
    vtkNew(vtkPolyData, slicePd);
    this->ExtractionCut(leftovers, cutPlane, 0, 1, slicePd);
    vtkNew(vtkPolyData, frontPd);
    this->ExtractionCut(currentPd, cutPlane, 1, 0, frontPd);

    //Get closestPoint region
    vtkNew(vtkPolyData, connectedPd);
    this->GetClosestPointConnectedRegion(slicePd, origin, connectedPd);

    if (connectedPd->GetNumberOfCells() == numCells)
    {
      done = 1;
    }
    this->ReplaceDataOnCells(connectedPd, sliceIds,
                             this->TotalSliceId, -1, this->InternalIdsArrayName);

    double surgeryPointIds[4];
    for (int i=0; i<4; i++)
    {
      surgeryPointIds[i] = surgeryPoints->GetId(i);

      //Test polydata
      double xyz[3]; this->WorkPd->GetPoint(surgeryPointIds[i], xyz);
      testPts->InsertNextPoint(xyz);

    }
    if (sliceId == 0)
    {
      this->TopNormals->InsertTuple(this->TotalSliceId, zvec);
      this->RightNormals->InsertTuple(this->TotalSliceId, xvec);
      this->StartPtIds->InsertNextTuple(surgeryPointIds);
    }
    double contourClosePt[3];
    branchCenterline->GetPoint(linePtId, contourClosePt);
    if (!done)
    {
      this->GetNextSurgeryPoints(connectedPd, contourClosePt, surgeryPoints, xvec, zvec, surgeryLineIds);
    }
    else
    {
      this->GetNextSurgeryPoints(currentPd, contourClosePt, surgeryPoints, xvec, zvec, surgeryLineIds);
    }

    sliceId++;
  }
  double surgeryPointIds[4];
  for (int i=0; i<4; i++)
  {
    surgeryPointIds[i] = surgeryPoints->GetId(i);

    //Test polydata
    double xyz[3]; this->WorkPd->GetPoint(surgeryPointIds[i], xyz);
    testPts->InsertNextPoint(xyz);
  }
  this->StartPtIds->InsertNextTuple(surgeryPointIds);
  this->AddSurgeryPoints(surgeryLineIds, surgeryPts, surgeryLines, surgeryData);

  this->ReplaceDataOnCells(branchPd, sliceIds,
                           this->TotalSliceId, -1, this->InternalIdsArrayName);
  this->NumberOfCubes->InsertNextTuple3( branchId, 1, this->TotalSliceId);
  this->TotalSliceId++;

  for (int i=0; i<4; i++)
  {
    testCells->InsertNextCell(testPts->GetNumberOfPoints()/4);
    for (int j=i;j<testPts->GetNumberOfPoints(); j+=4)
    {
      testCells->InsertCellPoint(j);
    }
  }
  vtkNew(vtkPolyData, testPd);
  testPd->SetPoints(testPts);
  testPd->SetLines(testCells);
  testPd->BuildLinks();

  //std::stringstream out;
  //out << branchId;

  //std::string filename = "/Users/adamupdegrove/Desktop/tmp/Lines_"+out.str()+".vtp";
  //vtkNew(vtkXMLPolyDataWriter, writer);
  //writer->SetInputData(testPd);
  //writer->SetFileName(filename.c_str());
  //writer->Write();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::SliceBifurcations()
{
  vtkDataArray *groupIds = this->WorkPd->GetCellData()->GetArray(this->GroupIdsArrayName);
  int numIds = groupIds->GetNumberOfTuples();
  vtkNew(vtkIntArray, segmentIds);
  segmentIds->SetNumberOfTuples(numIds);
  segmentIds->SetName(this->SegmentIdsArrayName);
  segmentIds->FillComponent(0, -1);
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(this->WorkPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  this->WorkPd->DeepCopy(ider->GetOutput());

  //this->CenterlineGraph->PrintGraph();
  int numSegs = this->CenterlineGraph->NumberOfCells;
  //this->MaxSegmentNumber = this->CenterlineGraph->NumberOfCells;
  this->MaxSegmentNumber = this->MaxGroupNumber;
  for (int i=0; i<numSegs; i++)
  {
    svGCell *gCell = this->CenterlineGraph->GetCell(i);
    if (gCell->Children[0] != NULL && gCell->Children[1] != NULL)
    {
      this->MaxSegmentNumber++;
      this->SliceBifurcation(this->WorkPd, gCell, this->MaxSegmentNumber, segmentIds);
    }
  }

  this->ReplaceDataOnCells(this->WorkPd, segmentIds,
                           0, -1, this->InternalIdsArrayName);
  this->WorkPd->GetCellData()->AddArray(segmentIds);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetRegionsOnPd(vtkPolyData *pd,
                                                  vtkIdList *regionIds)
{
  int numCells = pd->GetNumberOfCells();
  vtkIntArray *groupIds = vtkIntArray::SafeDownCast(
    pd->GetCellData()->GetArray(this->GroupIdsArrayName));
  for (int i=0; i<numCells; i++)
  {
    int groupId = groupIds->GetValue(i);
    if (regionIds->IsId(groupId) == -1)
    {
      regionIds->InsertNextId(groupId);
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
int vtkPolyDataSliceAndDiceFilter::SliceBifurcation(vtkPolyData *pd,
                                                    svGCell *gCell,
                                                    int &segmentId,
                                                    vtkDataArray *segmentIds)
{
  int parentId = gCell->GroupId;
  double vecs[2][3];
  vtkMath::Subtract(gCell->EndPt, gCell->StartPt, vecs[0]);
  vtkMath::Subtract(gCell->Children[0]->EndPt, gCell->Children[0]->StartPt, vecs[1]);

  int val = std::round(vtkMath::Dot(vecs[0], vecs[1]));
  int goodKidId = -1;
  int badKidId  = -1;
  if (val == 1)
  {
    goodKidId = gCell->Children[0]->GroupId;
    badKidId  = gCell->Children[1]->GroupId;
  }
  else
  {
    goodKidId = gCell->Children[1]->GroupId;
    badKidId  = gCell->Children[0]->GroupId;
  }

  vtkNew(vtkPolyData, centerlineBranchPd);
  vtkPolyDataSliceAndDiceFilter::ThresholdPd(this->Centerlines, badKidId, badKidId, 1,
    this->GroupIdsArrayName, centerlineBranchPd);
  vtkNew(vtkPolyData, branchCenterline);
  //Need to get just first cell of centerlines. There are duplicate for each centerline running through
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(centerlineBranchPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  vtkNew(vtkPolyData, tmpPd);
  tmpPd->ShallowCopy(ider->GetOutput());
  vtkPolyDataSliceAndDiceFilter::ThresholdPd(tmpPd, 0, 0, 1,
                                             this->InternalIdsArrayName,
                                             branchCenterline);
  double inPt0[3], inPt1[3];
  branchCenterline->GetPoint(0, inPt0);
  branchCenterline->GetPoint(1, inPt1);

  vtkNew(vtkPolyData, section0Pd);
  vtkNew(vtkPolyData, section1Pd);
  vtkNew(vtkPolyData, section2Pd);
  vtkNew(vtkPolyData, section2Loop);
  vtkPolyDataSliceAndDiceFilter::ThresholdPd(pd, parentId,  parentId,  1,
                                             this->GroupIdsArrayName, section0Pd);
  vtkPolyDataSliceAndDiceFilter::ThresholdPd(pd, goodKidId, goodKidId, 1,
                                             this->GroupIdsArrayName, section1Pd);
  vtkPolyDataSliceAndDiceFilter::ThresholdPd(pd, badKidId,  badKidId,  1,
                                             this->GroupIdsArrayName, section2Pd);

  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(section2Pd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();
  vtkPolyDataSliceAndDiceFilter::GetClosestPointConnectedRegion(boundaries->GetOutput(), inPt0, section2Loop);

  vtkNew(vtkAppendPolyData, appender);
  appender->AddInputData(section0Pd);
  appender->AddInputData(section1Pd);
  appender->Update();

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(appender->GetOutput());
  cleaner->Update();

  vtkNew(vtkFindSeparateRegions, separator);
  separator->SetInputData(cleaner->GetOutput());
  separator->SetOutPointArrayName(this->BoundaryPointsArrayName);
  separator->SetArrayName(this->GroupIdsArrayName);
  separator->Update();

  vtkIntArray *isBoundary = vtkIntArray::SafeDownCast(
    separator->GetOutput()->GetPointData()->GetArray(this->BoundaryPointsArrayName));
  vtkNew(vtkPoints, ringPoints);
  for (int i=0; i< separator->GetOutput()->GetNumberOfPoints(); i++)
  {
    if (isBoundary->GetValue(i) == 1)
    {
      ringPoints->InsertNextPoint(separator->GetOutput()->GetPoint(i));
    }
  }
  double centroid[3];
  vtkPolyDataSliceAndDiceFilter::GetCentroidOfPoints(ringPoints, centroid);

  std::list<int> criticalPoints;
  vtkPolyDataSliceAndDiceFilter::GetCommonValues(this->CriticalPointMap, parentId, goodKidId, criticalPoints);
  if (criticalPoints.size() != 2)
  {
    fprintf(stderr,"There should be two critical points between groups and there are %lu\n", criticalPoints.size());
    return 0;
  }

  double pt0[3], pt1[3];
  this->WorkPd->GetPoint(criticalPoints.front(), pt0);
  this->WorkPd->GetPoint(criticalPoints.back(), pt1);
  double vec0[3], vec1[3], normal[3], normal2[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Subtract(centroid, pt0, vec1);
  vtkMath::Cross(vec0, vec1, normal);
  vtkMath::Cross(vec0, vec1, normal2);
  vtkMath::Normalize(normal);
  vtkMath::Normalize(normal2);

  vtkNew(vtkIdList, goToPoints);
  vtkDataArray *section2Ids =
    section2Loop->GetPointData()->GetArray(this->InternalIdsArrayName);
  int front = section2Ids->LookupValue(criticalPoints.front());
  int back = section2Ids->LookupValue(criticalPoints.back());

  //fprintf(stdout,"FRONT!!!!: %d\n", criticalPoints.front());
  //fprintf(stdout,"BACK!!!!: %d\n", criticalPoints.back());
  this->GetSurgeryPoints(section2Loop, section2Ids, inPt0, inPt1, front, back, goToPoints);
  //fprintf(stdout,"ORDER: %lld %lld %lld %lld\n", goToPoints->GetId(0), goToPoints->GetId(1), goToPoints->GetId(2), goToPoints->GetId(3));
  vtkNew(vtkIdList, checkIds);
  this->GetPointGroups(this->WorkPd, this->GroupIdsArrayName, goToPoints->GetId(0), checkIds);
  if (checkIds->IsId(parentId) == -1)
  {
    //fprintf(stdout,"FLIPPP!\n");
    goToPoints->Reset();
    this->GetSurgeryPoints(section2Loop, section2Ids, inPt0, inPt1, back, front, goToPoints);
  }
  vtkNew(vtkPoints, tmpPts0);
  tmpPts0->InsertNextPoint(this->WorkPd->GetPoint(goToPoints->GetId(0)));
  tmpPts0->InsertNextPoint(this->WorkPd->GetPoint(goToPoints->GetId(1)));
  double goToFirst[3];
  vtkPolyDataSliceAndDiceFilter::GetCentroidOfPoints(tmpPts0, goToFirst);
  vtkNew(vtkPoints, tmpPts1);
  tmpPts1->InsertNextPoint(this->WorkPd->GetPoint(goToPoints->GetId(2)));
  tmpPts1->InsertNextPoint(this->WorkPd->GetPoint(goToPoints->GetId(3)));
  double goToSecond[3];
  vtkPolyDataSliceAndDiceFilter::GetCentroidOfPoints(tmpPts1, goToSecond);

  double vec2[3], vec3[3];
  vtkMath::Subtract(goToFirst, centroid, vec2);
  vtkMath::Subtract(goToSecond, centroid, vec3);
  double dot0 = vtkMath::Dot(vec2, normal);
  double dot1 = vtkMath::Dot(vec3, normal2);

  vtkMath::MultiplyScalar(normal, dot0);
  vtkMath::MultiplyScalar(normal2, dot1);

  double newSlicePt0[3], newSlicePt1[3];
  double adjSlicePt0[3], adjSlicePt1[3];
  vtkMath::Add(centroid, normal, newSlicePt0);
  vtkMath::MultiplyScalar(normal, 0.99);
  vtkMath::Add(centroid, normal, adjSlicePt0);
  vtkMath::Add(centroid, normal2, newSlicePt1);
  vtkMath::MultiplyScalar(normal2, 0.99);
  vtkMath::Add(centroid, normal2, adjSlicePt1);

  double pt2[3], pt3[3];
  tmpPts0->GetPoint(0, pt2);
  tmpPts0->GetPoint(1, pt3);

  double vec4[3], vec5[3];
  vtkMath::Subtract(pt2, newSlicePt0, vec4);
  vtkMath::Subtract(pt3, newSlicePt0, vec5);

  double newSliceNormal0[3];
  vtkMath::Cross(vec4, vec5, newSliceNormal0);
  vtkMath::Normalize(newSliceNormal0);

  if (vtkMath::Dot(newSliceNormal0, normal) < 0.0)
    vtkMath::MultiplyScalar(newSliceNormal0, -1.0);

  vtkNew(vtkPlane, cutPlane0);
  cutPlane0->SetOrigin(adjSlicePt0);
  cutPlane0->SetNormal(newSliceNormal0);

  double pt4[3], pt5[3];
  tmpPts1->GetPoint(0, pt4);
  tmpPts1->GetPoint(1, pt5);

  double vec6[3], vec7[3];
  vtkMath::Subtract(pt4, newSlicePt1, vec6);
  vtkMath::Subtract(pt5, newSlicePt1, vec7);

  double newSliceNormal1[3];
  vtkMath::Cross(vec6, vec7, newSliceNormal1);
  vtkMath::Normalize(newSliceNormal1);

  double xvec[3];
  vtkMath::Subtract(goToFirst, newSlicePt0, xvec);
  this->BifTopNormals->InsertNextTuple(newSliceNormal0);
  this->BifRightNormals->InsertNextTuple(xvec);

  if (vtkMath::Dot(newSliceNormal1, normal2) < 0.0)
    vtkMath::MultiplyScalar(newSliceNormal1, -1.0);

  vtkNew(vtkPlane, cutPlane1);
  cutPlane1->SetOrigin(adjSlicePt1);
  cutPlane1->SetNormal(newSliceNormal1);

  vtkNew(vtkPolyData, slicePd0);
  vtkNew(vtkPolyData, slicePd1);
  this->ExtractionCut(separator->GetOutput(), cutPlane0, 1, 1, slicePd0);
  this->ExtractionCut(slicePd0, cutPlane1, 1, 1, slicePd1);

  vtkNew(vtkGetBoundaryFaces, bigregioner);
  bigregioner->SetInputData(slicePd1);
  bigregioner->SetBoundaryEdges(1);
  bigregioner->SetManifoldEdges(0);
  bigregioner->SetNonManifoldEdges(0);
  bigregioner->SetFeatureEdges(0);
  bigregioner->SetRegionIdsArrayName("ConnectedRegionId");
  bigregioner->SetExtractLargestRegion(1);
  bigregioner->Update();
  vtkNew(vtkPolyData, leftovers);
  leftovers->DeepCopy(bigregioner->GetOutput());

  //Must check and make sure that the new region actually has the points;
  //sometimes gets cut short
  this->CheckStartSurgeryPoints(leftovers, goToPoints);

  this->ReplaceDataOnCells(leftovers, segmentIds,
                           segmentId, -1, this->InternalIdsArrayName);

  this->ReplaceDataOnCells(section0Pd, segmentIds,
                           parentId, -1, this->InternalIdsArrayName);

  this->ReplaceDataOnCells(section1Pd, segmentIds,
                           goodKidId, -1, this->InternalIdsArrayName);

  this->ReplaceDataOnCells(section2Pd, segmentIds,
                           badKidId, -1, this->InternalIdsArrayName);

  vtkNew(vtkPolyData, workCopy);
  workCopy->DeepCopy(this->WorkPd);
  workCopy->GetCellData()->AddArray(segmentIds);

  vtkNew(vtkIdList, fixedSurgeryPoints0);
  this->CriticalSurgeryPoints(workCopy, parentId, goToPoints->GetId(0), goToPoints->GetId(1),
                             centroid, adjSlicePt0, fixedSurgeryPoints0);
  double new0[4];
  //fprintf(stdout,"TEST TRY: %lld %lld %lld %lld\n", goToPoints->GetId(1),
  //                                                  goToPoints->GetId(0),
  //                                                  fixedSurgeryPoints0->GetId(0),
  //                                                  fixedSurgeryPoints0->GetId(1));
  new0[0] = goToPoints->GetId(1);
  new0[1] = goToPoints->GetId(0);
  new0[2] = fixedSurgeryPoints0->GetId(0);
  new0[3] = fixedSurgeryPoints0->GetId(1);
  for (int i=0; i<4; i++)
  {
    this->SurgeryPointMap.insert(std::make_pair(parentId, new0[i]));
  }

  vtkNew(vtkIdList, fixedSurgeryPoints1);
  this->CriticalSurgeryPoints(workCopy, goodKidId, goToPoints->GetId(2), goToPoints->GetId(3),
                             centroid, adjSlicePt1, fixedSurgeryPoints1);
  double new1[4];
  //fprintf(stdout,"TEST TRY: %lld %lld %lld %lld\n", fixedSurgeryPoints1->GetId(0),
  //                                                  fixedSurgeryPoints1->GetId(1),
  //                                                  goToPoints->GetId(3),
  //                                                  goToPoints->GetId(2));
  new1[0] = fixedSurgeryPoints1->GetId(0);
  new1[1] = fixedSurgeryPoints1->GetId(1);
  new1[2] = goToPoints->GetId(3);
  new1[3] = goToPoints->GetId(2);
  for (int i=0; i<4; i++)
  {
    this->SurgeryPointMap.insert(std::make_pair(goodKidId, new1[i]));
  }

  double new2[4];
  //fprintf(stdout,"TEST TRY: %lld %lld %lld %lld\n", goToPoints->GetId(2),
  //                                                  goToPoints->GetId(3),
  //                                                  goToPoints->GetId(0),
  //                                                  goToPoints->GetId(1));
  new2[0] = goToPoints->GetId(2);
  new2[1] = goToPoints->GetId(3);
  new2[2] = goToPoints->GetId(0);
  new2[3] = goToPoints->GetId(1);
  for (int i=0; i<4; i++)
  {
    this->SurgeryPointMap.insert(std::make_pair(badKidId, new2[i]));
  }

  //fprintf(stdout,"TEST TRY: %lld %lld %lld %lld %lld %lld %lld %lld\n", goToPoints->GetId(0),
  //                                                                      fixedSurgeryPoints0->GetId(0),
  //                                                                      fixedSurgeryPoints0->GetId(1),
  //                                                                      goToPoints->GetId(1),
  //                                                                      goToPoints->GetId(2),
  //                                                                      fixedSurgeryPoints1->GetId(0),
  //                                                                      fixedSurgeryPoints1->GetId(1),
  //                                                                      goToPoints->GetId(3));

  double new3[8];
  new3[0] = goToPoints->GetId(0),
  new3[1] = fixedSurgeryPoints0->GetId(0),
  new3[2] = fixedSurgeryPoints0->GetId(1),
  new3[3] = goToPoints->GetId(1),
  new3[4] = goToPoints->GetId(2),
  new3[5] = fixedSurgeryPoints1->GetId(0),
  new3[6] = fixedSurgeryPoints1->GetId(1),
  new3[7] = goToPoints->GetId(3);
  for (int i=0; i<8; i++)
  {
    this->SurgeryPointMap.insert(std::make_pair(segmentId, new3[i]));
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::CheckStartSurgeryPoints(vtkPolyData *pd, vtkIdList *startPoints)
{
  vtkDataArray *groupIds = pd->GetPointData()->GetArray(this->InternalIdsArrayName);

  vtkNew(vtkPointLocator, pointLocator);
  pointLocator->SetDataSet(pd);
  pointLocator->BuildLocator();

  for (int i=0; i<startPoints->GetNumberOfIds(); i++)
  {
    int pdId = groupIds->LookupValue(startPoints->GetId(i));
    if (pdId == -1)
    {
      double pt[3];
      this->WorkPd->GetPoint(startPoints->GetId(i), pt);
      int findPtId = pointLocator->FindClosestPoint(pt);
      startPoints->SetId(i, groupIds->GetTuple1(findPtId));
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
int vtkPolyDataSliceAndDiceFilter::CriticalSurgeryPoints(vtkPolyData *pd,
                                                        const int groupId,
                                                        const int frontId,
                                                        const int backId,
                                                        double startPt[3],
                                                        double secondPt[3],
                                                        vtkIdList *fixedSurgeryPoints)
{
  vtkNew(vtkPolyData, thresholdPd);
  this->ThresholdPd(pd, groupId, groupId, 1,
      this->SegmentIdsArrayName, thresholdPd);

  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(thresholdPd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();
  vtkNew(vtkPolyData, thresholdLoop);
  this->GetClosestPointConnectedRegion(boundaries->GetOutput(), startPt, thresholdLoop);

  vtkDataArray *thresholdIds =
    thresholdLoop->GetPointData()->GetArray(this->InternalIdsArrayName);
  vtkNew(vtkIdList, fixedSurgeryPoints0);
  int fixFront = thresholdIds->LookupValue(frontId);
  int fixBack = thresholdIds->LookupValue(backId);
  this->GetSurgeryPoints(thresholdLoop, thresholdIds, startPt, secondPt, fixFront, fixBack, fixedSurgeryPoints);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetNumberOfSlicesBefore(const int id, int &numBefore)
{
  int numTups = this->NumberOfCubes->GetNumberOfTuples();
  for (int i=0; i<numTups; i++)
  {
    if (this->NumberOfCubes->GetComponent(i, 0) == id)
    {
      numBefore = this->NumberOfCubes->GetComponent(i, 2) -
                  this->NumberOfCubes->GetComponent(i, 1);
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
int vtkPolyDataSliceAndDiceFilter::GetCentroidOfPoints(vtkPoints *points,
                                                       double centroid[3])
{
  int numPoints = points->GetNumberOfPoints();
  centroid[0] = 0.0; centroid[1] = 0.0; centroid[2] = 0.0;
  for (int i=0; i<numPoints; i++)
  {
    double pt[3];
    points->GetPoint(i, pt);
    for (int j=0; j<3; j++)
    {
      centroid[j] += pt[j];
    }
  }

  vtkMath::MultiplyScalar(centroid, 1.0/numPoints);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::ThresholdPd(vtkPolyData *pd, int minVal,
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

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  returnPd->DeepCopy(surfacer->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
bool vtkPolyDataSliceAndDiceFilter::ListsMatch(vtkIdList *listA, vtkIdList *listB)
{
  int numA = listA->GetNumberOfIds();
  int numB = listB->GetNumberOfIds();
  if (numA != numB)
  {
    fprintf(stderr,"Lists are different sizes\n");
    return 0;
  }
  int match=0;
  for (int i=0; i<numA; i++)
  {
    if (listA->GetId(i) == listB->GetId(i))
    {
      match++;
    }
  }
  if (match == numA)
    return 1;
  else
    return 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::BuildPolycube()
{
  vtkDataArray *groupIds = this->GraphPd->GetCellData()->GetArray(this->GroupIdsArrayName);
  int numIds = groupIds->GetNumberOfTuples();
  //vtkNew(vtkPolyData, bigGraphPd);
  //this->BlowUpSkeleton(bigGraphPd);

  double minmax[2];
  groupIds->GetRange(minmax);
  int count = 0;
  for (int i=0; i<=this->MaxSegmentNumber; i++)
  {
    vtkNew(vtkThreshold, thresholder);
    thresholder->SetInputData(this->WorkPd);
    thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->SegmentIdsArrayName);
    thresholder->ThresholdBetween(i, i);
    thresholder->Update();

    int numPoints = thresholder->GetOutput()->GetNumberOfPoints();

    if (numPoints != 0)
    {
      int numCurrentCubes = this->Polycube->GetNumberOfGrids();
      int numCubes;;

      double cornerPt[3], secondPt[3];
      double vec0[3];
      if (i <= this->MaxGroupNumber)
      {
        numCubes = this->NumberOfCubes->GetComponent(count, 1);
        int cellId = groupIds->LookupValue(i);
        vtkNew(vtkIdList, neighbors);
        vtkIdType npts, *pts;
        this->GraphPd->GetCellPoints(cellId, npts, pts);
        this->GraphPd->GetPointCells(pts[0], neighbors);
        if (i == 0)
        {
          this->GraphPd->GetPoint(pts[1], cornerPt);
          this->GraphPd->GetPoint(pts[0], secondPt);
        }
        else
        {
          this->GraphPd->GetPoint(pts[0], cornerPt);
          this->GraphPd->GetPoint(pts[1], secondPt);
        }
        vtkMath::Subtract(secondPt, cornerPt, vec0);
      }
      else
      {
        cornerPt[0] = 0.0; cornerPt[1] = 0.0; cornerPt[2] = 0.0;
        numCubes = 1;
      }
      this->Polycube->SetNumberOfGrids(numCurrentCubes + numCubes);
      int idStart = count++;

      double unitCube[3]; unitCube[0] = 1.0; unitCube[1] = 1.0; unitCube[2] = 1.0;
      for (int j=numCurrentCubes; j<numCurrentCubes + numCubes; j++)
      {
        double topNormal[3];
        double rightNormal[3];
        int starts[8];
        if (i <= this->MaxGroupNumber)
        {
          for (int k=0; k<4; k++)
          {
            starts[k] = this->StartPtIds->GetComponent(j+idStart, k);
            starts[k+4] = this->StartPtIds->GetComponent(j+idStart+1, k);
          }
          this->TopNormals->GetTuple(j, topNormal);
          this->RightNormals->GetTuple(j, rightNormal);
        }
        else
        {
          std::list<int> surgeryIds;
          this->GetValuesFromMap(this->SurgeryPointMap, i, surgeryIds);
          std::list<int>::iterator it = surgeryIds.begin();
          for (int j=0; it != surgeryIds.end() && j<8; it++, j++)
          {
            starts[j] = *it;
          }
          this->BifTopNormals->GetTuple(i-this->MaxGroupNumber-1, topNormal);
          this->BifRightNormals->GetTuple(i-this->MaxGroupNumber-1, rightNormal);
        }
        this->Polycube->SetGridWithOrigin(j, cornerPt, unitCube, 0, topNormal, rightNormal, starts);
        for (int k=0; k<3; k++)
        {
          cornerPt[k] = cornerPt[k] + vec0[k];
        }
      }
    }
  }

  return 1;
}

////---------------------------------------------------------------------------
///**
// * @brief
// * @param *pd
// * @return
// */
//int vtkPolyDataSliceAndDiceFilter::BlowUpSkeleton(vtkPolyData *bigGraphPd)
//{
//  //RecurseSkeleton , start at 0 blow up, and then continue along tree
//
//  this->Recurse();
//  return 1;
//}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::ReplaceDataOnCells(vtkPointSet *pointset,
                                                      vtkDataArray *sliceIds,
                                                      const int sliceId,
                                                      const int replaceVal,
                                                      const std::string &arrName)
{
  int numCells = pointset->GetNumberOfCells();
  vtkDataArray *cellIds = pointset->GetCellData()->GetArray(arrName.c_str());

  for (int i=0; i<numCells; i++)
  {
    int cellId = cellIds->GetTuple1(i);
    int currVal = sliceIds->GetTuple1(cellId);
    if (currVal == replaceVal)
    {
      sliceIds->SetTuple1(cellId, sliceId);
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
int vtkPolyDataSliceAndDiceFilter::GetBranchLength(vtkPolyData *points, double &length)
{
  int numPts = points->GetNumberOfPoints();

  length = 0.0;

  for (int i=1; i<numPts; i++)
  {
    double pt0[3], pt1[3];
    points->GetPoint(i-1, pt0);
    points->GetPoint(i, pt1);

    length += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                        std::pow(pt1[1] - pt0[1], 2.0) +
                        std::pow(pt1[2] - pt0[2], 2.0));
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetCloseGeodesicPoint(vtkPolyData *pd, double centerPt[3], const int startPtId, int &returnStartId, double zvec[3],
                                                         vtkPolyData *boundary)
{
  int actualId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(startPtId);

  vtkNew(vtkFindGeodesicPath, finder);
  finder->SetInputData(pd);
  finder->SetStartPtId(actualId);
  finder->SetClosePt(centerPt);
  finder->SetDijkstraArrayName(this->DijkstraArrayName);
  finder->Update();

  boundary->DeepCopy(finder->GetBoundary());
  int endPtId = finder->GetEndPtId();
  returnStartId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
    GetTuple1(endPtId);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetNextSurgeryPoints(vtkPolyData *pd, double centerPt[3], vtkIdList *surgeryPoints, double xvec[3], double zvec[3], vtkIdList *surgeryLineIds)
{
  int contourPtId;
  vtkNew(vtkPolyData, boundary);
  //this->GetCloseGeodesicPoint(pd, centerPt, surgeryPoints->GetId(0),
  //  contourPtId, zvec, boundary);
  int initialSurgeryPt = surgeryPoints->GetId(0);
  this->GetClose3DPoint(pd, centerPt, surgeryPoints->GetId(0),
    contourPtId, zvec, boundary);
  surgeryPoints->SetId(0, contourPtId);

  int pointId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(contourPtId);

  vtkNew(vtkIdList, startCellIds);
  boundary->BuildLinks();
  boundary->GetPointCells(pointId, startCellIds);
  int cellId = startCellIds->GetId(0);

  vtkIdType npts, *pts;
  boundary->GetCellPoints(cellId, npts, pts);
  int secondPtId;
  if (pts[0] == pointId)
  {
    secondPtId = pts[1];
  }
  else
  {
    secondPtId = pts[0];
  }
  double pt0[3], pt1[3];
  boundary->GetPoint(pointId, pt0);
  boundary->GetPoint(secondPtId, pt1);
  double vec0[3], vec1[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  vtkMath::Cross(xvec, vec0, vec1);

  if (vtkMath::Dot(zvec, vec1) > 0)
  {
    cellId = startCellIds->GetId(1);
  }

  //Getting full loop length
  int front = pointId;
  int back = pointId;
  double length = 0.0;
  int iter = 0;
  int prevCellId = cellId;
  while (pointId != back || iter == 0)
  {
    boundary->GetPoint(pointId, pt0);
    this->IteratePoint(boundary, pointId, cellId);
    boundary->GetPoint(pointId, pt1);
    length += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                        std::pow(pt1[1] - pt0[1], 2.0) +
                        std::pow(pt1[2] - pt0[2], 2.0));
    iter++;
  }

  //Finding the surgery points which separate boundary into four points
  double surgeryCount = 1.0;
  double currLength = 0.0;;
  back = pointId;
  pointId = front;
  prevCellId = cellId;
  iter = 0;
  while (pointId != back || iter == 0)
  {
    boundary->GetPoint(pointId, pt0);
    this->IteratePoint(boundary, pointId, prevCellId);
    boundary->GetPoint(pointId, pt1);
    currLength += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                            std::pow(pt1[1] - pt0[1], 2.0) +
                            std::pow(pt1[2] - pt0[2], 2.0));
    if (currLength > length*(surgeryCount/4.0))
    {
      if (surgeryCount < 4)
      {
        int newId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
          GetTuple1(pointId);
        surgeryPoints->SetId(surgeryCount, newId);
      }
      surgeryCount += 1.0;
    }
    iter++;
  }
  vtkNew(vtkFindGeodesicPath, finder);
  finder->SetInputData(pd);
  finder->SetStartPtId(pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(surgeryPoints->GetId(0)));
  finder->SetEndPtId(pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(initialSurgeryPt));
  finder->SetDijkstraArrayName(this->DijkstraArrayName);
  finder->SetInternalIdsArrayName(this->InternalIdsArrayName);
  finder->Update();

  vtkNew(vtkIdList, tmpIds);
  tmpIds = finder->GetPathIds();
  int numToAdd = tmpIds->GetNumberOfIds();
  for (int i=0; i<numToAdd; i++)
  {
    int testId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
      GetTuple1(tmpIds->GetId(i));
    if (surgeryLineIds->IsId(testId) == -1)
    {
      surgeryLineIds->InsertNextId(testId);
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
int vtkPolyDataSliceAndDiceFilter::GetFirstSurgeryPoints(vtkPolyData *pd, int pointId, vtkIdList *surgeryPoints, double xvec[3], double zvec[3])
{
  vtkNew(vtkIdList, startCellIds);
  pd->BuildLinks();
  pd->GetPointCells(pointId, startCellIds);
  int cellId = startCellIds->GetId(0);

  vtkIdType npts, *pts;
  pd->GetCellPoints(cellId, npts, pts);
  int secondPtId;
  if (pts[0] == pointId)
  {
    secondPtId = pts[1];
  }
  else
  {
    secondPtId = pts[0];
  }
  double pt0[3], pt1[3];
  pd->GetPoint(pointId, pt0);
  pd->GetPoint(secondPtId, pt1);
  double vec0[3], vec1[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  vtkMath::Cross(xvec, vec0, vec1);

  if (vtkMath::Dot(zvec, vec1) > 0)
  {
    cellId = startCellIds->GetId(1);
  }

  //Getting full loop length
  int front = pointId;
  int back = pointId;
  double length = 0.0;
  int iter = 0;
  int prevCellId = cellId;
  while (pointId != back || iter == 0)
  {
    pd->GetPoint(pointId, pt0);
    this->IteratePoint(pd, pointId, cellId);
    pd->GetPoint(pointId, pt1);
    length += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                        std::pow(pt1[1] - pt0[1], 2.0) +
                        std::pow(pt1[2] - pt0[2], 2.0));
    iter++;
  }

  //Finding the surgery points which separate boundary into four points
  double surgeryCount = 1.0;
  double currLength = 0.0;;
  back = pointId;
  front = pointId;
  prevCellId = cellId;
  iter = 0;
  int newId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    GetTuple1(pointId);
  surgeryPoints->InsertNextId(newId);
  while (pointId != back || iter == 0)
  {
    pd->GetPoint(pointId, pt0);
    this->IteratePoint(pd, pointId, prevCellId);
    pd->GetPoint(pointId, pt1);
    currLength += std::sqrt(std::pow(pt1[0] - pt0[0], 2.0) +
                            std::pow(pt1[1] - pt0[1], 2.0) +
                            std::pow(pt1[2] - pt0[2], 2.0));
    if (currLength > length*(surgeryCount/4.0))
    {
      if (surgeryCount < 4)
      {
        newId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
          GetTuple1(pointId);
        surgeryPoints->InsertNextId(newId);
      }
      surgeryCount += 1.0;
    }
    iter++;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetClose3DPoint(vtkPolyData *pd, double centerPt[3],
                                                   const int startPtId, int &returnStartId,
                                                   double zvec[3],
                                                   vtkPolyData *boundary)
{
  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(pd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();
  //fprintf(stdout,"Input start Pt ID: %d\n", startPtId);

  this->GetClosestPointConnectedRegion(boundaries->GetOutput(), centerPt, boundary);

  int pointId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)
    ->LookupValue(startPtId);
  double startPt[3];
  pd->GetPoint(pointId, startPt);

  vtkNew(vtkPointLocator, pointLocator);
  pointLocator->SetDataSet(boundary);
  pointLocator->BuildLocator();

  int endPtId = pointLocator->FindClosestPoint(startPt);
  returnStartId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
    GetTuple1(endPtId);

  double minDist = 1.0e10;
  int minId = -1;
  for (int i=0; i<boundary->GetNumberOfPoints(); i++)
  {
    double pt[3];
    boundary->GetPoint(i, pt);
    double dist = std::sqrt(std::pow(pt[0] - startPt[0], 2.0) +
                            std::pow(pt[1] - startPt[1], 2.0) +
                            std::pow(pt[2] - startPt[2], 2.0));
    if (dist < minDist)
    {
      minDist = dist;
      minId = i;
    }
  }
  //fprintf(stdout,"Locator ID: %d\n", endPtId);
  //fprintf(stdout,"MinDist ID: %d\n", minId);


  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::GetContourSecondPoint(vtkPolyData *pd, int ptId, double centerPt[3], double zvec[3], int &returnSecondId)
{
  pd->BuildLinks();

  vtkSmartPointer<vtkIdList> cellIds =
    vtkSmartPointer<vtkIdList>::New();
  pd->GetPointCells(ptId, cellIds);
  double startPt[3];
  pd->GetPoint(ptId, startPt);

  int numIds = cellIds->GetNumberOfIds();
  if (numIds != 2)
  {
    vtkErrorMacro("Number of connected contour ids is not 2");
    return 0;
  }

  vtkIdType npts, *pts;
  pd->GetCellPoints(cellIds->GetId(0), npts, pts);
  int checkPtId = -1;
  if (pts[0] != ptId)
  {
    checkPtId = pts[0];
  }
  else
  {
    checkPtId = pts[1];
  }

  double checkPt[3];
  pd->GetPoint(checkPtId, checkPt);
  double vec0[3], vec1[3], vec2[3];
  vtkMath::Subtract(startPt, centerPt, vec0);
  vtkMath::Subtract(checkPt, startPt, vec1);
  vtkMath::Cross(vec0, vec1, vec2);
  double dotCheck = vtkMath::Dot(vec2, zvec);
  if (dotCheck > 0)
  {
    returnSecondId = checkPtId;
    return 1;
  }

  pd->GetCellPoints(cellIds->GetId(1), npts, pts);
  if (pts[0] != ptId)
  {
    checkPtId = pts[0];
  }
  else
  {
    checkPtId = pts[1];
  }
  returnSecondId = checkPtId;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPolyDataSliceAndDiceFilter::AddSurgeryPoints(vtkIdList *surgeryLineIds,
                                                    vtkPoints *surgeryPts,
                                                    vtkCellArray *surgeryLines,
                                                    vtkIntArray *surgeryData)
{
  int numToAdd = surgeryLineIds->GetNumberOfIds();
  vtkNew(vtkIdList, newPointIds);
  newPointIds->SetNumberOfIds(numToAdd);
  for (int i=0; i<numToAdd; i++)
  {
    double pt[3];
    surgeryData->InsertNextValue(surgeryLineIds->GetId(i));
    this->WorkPd->GetPoint(surgeryLineIds->GetId(i), pt);
    int newPtId = surgeryPts->InsertNextPoint(pt);
    newPointIds->SetId(i, newPtId);
  }
  surgeryLines->InsertNextCell(newPointIds);
  return 1;
}
