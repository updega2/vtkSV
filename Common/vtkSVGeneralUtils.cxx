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
#include "vtkCenterOfMass.h"
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
int vtkSVGeneralUtils::CheckSurface(vtkPolyData *pd)
{
  pd->BuildLinks();

  int numPts = pd->GetNumberOfPoints();
  int numPolys = pd->GetNumberOfCells();

  for (int i=0; i<numPolys; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    if (npts != 3)
    {
      //vtkErrorMacro("Surface contains elements that aren't triangles");
      return 0;
    }
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0, p1;
      p0 = pts[j];
      p1 = pts[(j+1)%npts];

      vtkNew(vtkIdList, edgeNeighbor);
      pd->GetCellEdgeNeighbors(i, p0, p1, edgeNeighbor);

      if (edgeNeighbor->GetNumberOfIds() > 1)
      {
        //vtkErrorMacro("Surface contains triangles with multiple neighbors, not manifold");
        return 0;
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
int vtkSVGeneralUtils::ComputeArea(double pt0[], double pt1[],
                                   double pt2[], double &area)
{
  area = 0.0;
  area += (pt0[0]*pt1[1])-(pt1[0]*pt0[1]);
  area += (pt1[0]*pt2[1])-(pt2[0]*pt1[1]);
  area += (pt2[0]*pt0[1])-(pt0[0]*pt2[1]);
  area *= 0.5;

  return 1;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputeMassCenter(vtkPolyData *pd, double massCenter[3])
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
int vtkSVGeneralUtils::GetBarycentricCoordinates(double f[3], double pt0[3],
                                                 double pt1[3], double pt2[3],
					                                       double &a0, double &a1, double &a2)
{
  double f0[3], f1[3], f2[3], v0[3], v1[3];
  for (int i=0; i<3; i++)
  {
    v0[i] = pt0[i] - pt1[i];
    v1[i] = pt0[i] - pt2[i];
    f0[i] = pt0[i] - f[i];
    f1[i] = pt1[i] - f[i];
    f2[i] = pt2[i] - f[i];
  }

  double vArea[3], vA0[3], vA1[3], vA2[3];
  vtkMath::Cross(v0, v1, vArea);
  vtkMath::Cross(f1, f2, vA0);
  vtkMath::Cross(f2, f0, vA1);
  vtkMath::Cross(f0, f1, vA2);

  double area = vtkMath::Norm(vArea);
  a0 = vtkMath::Norm(vA0)/area;// * Sign(vtkMath::Dot(vArea, vA0));
  a1 = vtkMath::Norm(vA1)/area;// * Sign(vtkMath::Dot(vArea, vA1));
  a2 = vtkMath::Norm(vA2)/area;// * Sign(vtkMath::Dot(vArea, vA2));
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetPointNeighbors(vtkIdType p0,
                                       vtkPolyData *pd,
						                           vtkIdList *pointNeighbors)
{
  //Assuming that pointNeighbors is set with no neighbors already
  vtkNew(vtkIdList, cellIdList);
  pd->GetPointCells(p0, cellIdList);

  for (int i=0; i<cellIdList->GetNumberOfIds(); i++)
  {
    vtkIdType cellId = cellIdList->GetId(i);
    vtkIdType npts, *pts;
    pd->GetCellPoints(cellId, npts, pts);

    for (int j=0; j<npts; j++)
    {
      vtkIdType neighborPoint = pts[j];
      if (neighborPoint != p0)
      {
        pointNeighbors->InsertUniqueId(neighborPoint);
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
int vtkSVGeneralUtils::GetEdgeCotangentAngle(double pt0[3], double pt1[3],
                                           double pt2[3], double &angle)
{
  double area = 0.0;
  vtkSVGeneralUtils::ComputeArea(pt0, pt1, pt2, area);
  if (area < 0)
  {
    double tmpPoint[3];
    for (int i=0; i<3; i++)
    {
      tmpPoint[i] = pt0[i];
      pt0[i] = pt1[i];
      pt1[i] = tmpPoint[i];
    }
  }
  double vec0[3];
  double vec1[3];
  for (int i=0; i<3; i++)
  {
    vec0[i] = pt0[i] - pt2[i];
    vec1[i] = pt1[i] - pt2[i];
  }
  double numerator = vtkMath::Dot(vec0, vec1);
  double cross[3];
  vtkMath::Cross(vec0, vec1, cross);
  double denominator = vtkMath::Norm(cross);
  angle = numerator/denominator;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::CreateEdgeTable(vtkPolyData *pd,
                                     vtkEdgeTable *edgeTable,
                                     vtkFloatArray *edgeWeights,
                                     vtkIntArray *edgeNeighbors,
                                     vtkIntArray *isBoundary)
{
  int numPts = pd->GetNumberOfPoints();
  int numTris = pd->GetNumberOfCells();

  edgeTable->InitEdgeInsertion(numPts, 1);
  isBoundary->SetNumberOfValues(numPts);
  for (int i=0; i<numPts; i++)
  {
    isBoundary->InsertValue(i, 0);
  }
  for (int i=0; i<numTris; i++)
  {
    //Insert edge into table
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0 = pts[j];
      vtkIdType p1 = pts[(j+1)%npts];
      vtkNew(vtkIdList, neighborCellIds);
      pd->GetCellEdgeNeighbors(i, p0, p1, neighborCellIds);
      vtkIdType neighborCellId = 0;
      if (neighborCellIds->GetNumberOfIds() > 0)
      {
        neighborCellId = neighborCellIds->GetId(0);
      }
      else
      {
        neighborCellId = -1;
        isBoundary->InsertValue(p0, 1);
        isBoundary->InsertValue(p1, 1);
      }
      vtkIdType checkEdge = edgeTable->IsEdge(p0, p1);
      if (checkEdge == -1)
      {
        //Compute Edge Weight
        double weight = 0.0;
        vtkSVGeneralUtils::ComputeEdgeWeight(pd, i, neighborCellId,
                                             p0, p1, weight);
        vtkIdType edgeId = edgeTable->InsertEdge(p0, p1);
        edgeWeights->InsertValue(edgeId, weight);
        edgeNeighbors->InsertValue(edgeId, neighborCellId);
        if (weight < 0)
        {
          //vtkWarningMacro("Negative weight on edge between cells " << i <<
          //  " and "<< neighborCellId << ": " << weight);
        }
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
int vtkSVGeneralUtils::ComputeEdgeWeight(vtkPolyData *pd,
                                         vtkIdType cellId,
                                         vtkIdType neighborCellId,
                                         vtkIdType p0,
                                         vtkIdType p1,
                                         double &weight)
{
  //Add the edge weights based on the angle of the edge
  vtkIdType cellIds[2];
  cellIds[0] = cellId;
  cellIds[1] = neighborCellId;
  weight = 0.0;
  double v0[3]; double v1[3]; double v2[3];
  pd->GetPoint(p0, v0);
  pd->GetPoint(p1, v1);
  for (int i=0; i<2; i++)
  {
    vtkIdType npts, *pts;
    if (cellIds[i] != -1)
    {
      pd->GetCellPoints(cellIds[i], npts, pts);
      for (int k=0; k<npts; k++)
      {
        if (pts[k] != p0 && pts[k] != p1)
        {
          pd->GetPoint(pts[k], v2);
          double angle = 0.0;
          vtkSVGeneralUtils::GetEdgeCotangentAngle(v0, v1, v2, angle);

          weight += 0.5*angle;
        }
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

