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

/** @file vtkSVGeneralUtils.cxx
 *  @brief
 *  @brief
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVGeneralUtils.h"

#include "vtkCellData.h"
#include "vtkClipPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkExtractGeometry.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkThreshold.h"
#include "vtkTriangleFilter.h"
#include "vtkDataSet.h"
#include "vtkUnstructuredGrid.h"

// ----------------------
// CheckArrayExists
// ----------------------
/**
 * @brief Function to check is array with name exists in cell or point data
 * @param ds this is the object to check if the array exists
 * @param datatype this is point or cell. point=0,cell=1
 * @param arrayname this is the name of the array to check
 * @reutrn this returns 1 if the array exists and zero if it doesn't
 * or the function does not return properly.
 */

int vtkSVGeneralUtils::CheckArrayExists(vtkDataSet *ds,
                                        int datatype,
                                        std::string arrayname)
{
  int exists =0;

  if (datatype == 0)
  {
    int numArrays = ds->GetPointData()->GetNumberOfArrays();
    for (int i=0;i<numArrays;i++)
    {
      if (!strcmp(ds->GetPointData()->GetArrayName(i),arrayname.c_str()))
      {
	      exists =1;
      }
    }
  }
  else
  {
    int numArrays = ds->GetCellData()->GetNumberOfArrays();
    for (int i=0;i<numArrays;i++)
    {
      if (!strcmp(ds->GetCellData()->GetArrayName(i),arrayname.c_str()))
      {
	      exists =1;
      }
    }
  }

  return exists;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetClosestPointConnectedRegion(vtkPolyData *pd,
                                                      double origin[3])
{
  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(pd);
  connector->SetExtractionModeToClosestPointRegion();
  connector->SetClosestPoint(origin);
  connector->Update();

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(connector->GetOutput());
  surfacer->Update();

  pd->DeepCopy(surfacer->GetOutput());

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetClosestPointConnectedRegion(vtkPolyData *inPd,
                                                      double origin[3],
                                                      vtkPolyData *outPd)
{
  outPd->DeepCopy(inPd);
  vtkSVGeneralUtils::GetClosestPointConnectedRegion(outPd, origin);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::IteratePoint(vtkPolyData *pd, int &pointId, int &prevCellId)
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
int vtkSVGeneralUtils::ThresholdPd(vtkPolyData *pd, int minVal,
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
    //vtkDebugMacro("No points after threshold");
    return 0;
  }

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
int vtkSVGeneralUtils::GetCentroidOfPoints(vtkPoints *points,
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
int vtkSVGeneralUtils::GetPointGroups(vtkPolyData *pd, std::string arrayName,
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
int vtkSVGeneralUtils::ExtractionCut(vtkPolyData *inPd, vtkImplicitFunction *cutFunction,
                                     const int extractBoundaryCells,
                                     const int extractInside,
                                     vtkPolyData *outPd)
{
    vtkNew(vtkExtractGeometry, cutter);
    cutter->SetInputData(inPd);
    cutter->SetImplicitFunction(cutFunction);
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
int vtkSVGeneralUtils::ClipCut(vtkPolyData *inPd, vtkImplicitFunction *cutFunction,
                               const int generateClippedOutput,
                               const int extractInside,
                               vtkPolyData *outPd,
                               vtkPolyData *clippedOutPd)
{
    vtkNew(vtkClipPolyData, cutter);
    cutter->SetInputData(inPd);
    cutter->SetClipFunction(cutFunction);
    cutter->SetInsideOut(extractInside);
    cutter->SetGenerateClippedOutput(generateClippedOutput);
    cutter->Update();

    vtkNew(vtkTriangleFilter, triangulator);
    triangulator->SetInputData(cutter->GetOutput(0));
    triangulator->Update();

    outPd->DeepCopy(triangulator->GetOutput());
    if (generateClippedOutput)
    {
      triangulator->SetInputData(cutter->GetOutput(1));
      triangulator->Update();
      clippedOutPd->DeepCopy(triangulator->GetOutput());
    }

    return 1;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetPointsLength(vtkPolyData *points, double &length)
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
int vtkSVGeneralUtils::ReplaceDataOnCells(vtkPointSet *pointset,
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
int vtkSVGeneralUtils::GetCutPlane(const double endPt[3], const double startPt[3],
                                   const double length, double origin[3], vtkPlane *cutPlane)
{
  double normal[3];
  vtkMath::Subtract(endPt, startPt, normal);
  vtkMath::Normalize(normal);

  for (int i=0; i<3; i++)
    origin[i] = endPt[i];

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
int vtkSVGeneralUtils::GetAllMapKeys(std::multimap<int, int> &map,
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
int vtkSVGeneralUtils::GetAllMapValues(std::multimap<int, int> &map,
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
int vtkSVGeneralUtils::GetValuesFromMap(std::multimap<int, int> &map,
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
int vtkSVGeneralUtils::GetKeysFromMap(std::multimap<int, int> &map,
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
int vtkSVGeneralUtils::GetCommonValues(std::multimap<int, int> &map,
                                                   const int keyA, const int keyB,
                                                   std::list<int> &returnList)
{
  std::list<int> listA, listB;
  vtkSVGeneralUtils::GetValuesFromMap(map, keyA, listA);
  vtkSVGeneralUtils::GetValuesFromMap(map, keyB, listB);
  vtkSVGeneralUtils::ListIntersection(listA, listB, returnList);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetUniqueNeighbors(std::multimap<int, int> &map,
                                          const int key,
                                          std::list<int> &keyVals,
                                          std::list<int> &uniqueKeys)
{
  int numVals = keyVals.size();

  std::list<int>::iterator valit = keyVals.begin();
  for (int i=0; valit != keyVals.end(); ++valit)
  {
    std::list<int> valKeys;
    vtkSVGeneralUtils::GetKeysFromMap(map, *valit, valKeys);
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
int vtkSVGeneralUtils::ListIntersection(std::list<int> &listA,
                                                    std::list<int> &listB,
                                                    std::list<int> &returnList)
{
  std::set_intersection(listA.begin(), listA.end(),
                        listB.begin(), listB.end(),
                        std::inserter(returnList,returnList.begin()));

  return 1;
}

