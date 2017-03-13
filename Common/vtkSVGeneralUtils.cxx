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

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCenterOfMass.h"
#include "vtkClipPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkExtractGeometry.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkThreshold.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
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
      return SV_ERROR;
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
        return SV_ERROR;
      }
    }
  }
  return SV_OK;
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

  return SV_OK;
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

  return SV_OK;
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

  return SV_OK;
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
    return SV_ERROR;
  }

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(thresholder->GetOutput());
  surfacer->Update();

  returnPd->DeepCopy(surfacer->GetOutput());

  return SV_OK;
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
  return SV_OK;
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

  return SV_OK;
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

    return SV_OK;
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

    return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
double vtkSVGeneralUtils::Distance(double pt0[3], double pt1[3])
{
  return sqrt(pow(pt1[0] - pt0[0], 2.0) +
              pow(pt1[1] - pt0[1], 2.0) +
              pow(pt1[2] - pt0[2], 2.0));
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

    length += vtkSVGeneralUtils::Distance(pt0, pt1);
  }

  return SV_OK;
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

  return SV_OK;
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

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputeTriangleArea(double pt0[], double pt1[],
                                           double pt2[], double &area)
{
  area = 0.0;
  area += (pt0[0]*pt1[1])-(pt1[0]*pt0[1]);
  area += (pt1[0]*pt2[1])-(pt2[0]*pt1[1]);
  area += (pt2[0]*pt0[1])-(pt0[0]*pt2[1]);
  area *= 0.5;

  return SV_OK;
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

  return SV_OK;
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
  return SV_OK;
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
  return SV_OK;
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
  vtkSVGeneralUtils::ComputeTriangleArea(pt0, pt1, pt2, area);
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

  return SV_OK;
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
        vtkSVGeneralUtils::ComputeHarmonicEdgeWeight(pd, i, neighborCellId,
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

  return SV_OK;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputeHarmonicEdgeWeight(vtkPolyData *pd,
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

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ConvertFieldToPolyData(vtkPolyData *inPd, std::string fieldName, vtkPolyData *outPd)
{
  int numCells = inPd->GetNumberOfCells();
  int numPts   = inPd->GetNumberOfPoints();
  vtkNew(vtkPoints, fieldPts);;
  fieldPts->SetNumberOfPoints(numPts);
  vtkNew(vtkCellArray, fieldCells);
  fieldCells = inPd->GetPolys();

  vtkFloatArray *fieldArray;
  fieldArray = vtkFloatArray::SafeDownCast(
    inPd->GetPointData()->GetArray(fieldName.c_str()));
  for (int i=0; i<numCells; i++)
  {
    vtkIdType npts, *pts;
    inPd->GetCellPoints(i, npts, pts);
    for (int j=0; j<npts; j++)
    {
      double pt[3];
      fieldArray->GetTuple(pts[j], pt);
      fieldPts->SetPoint(pts[j], pt);
    }
  }

  outPd->SetPolys(fieldCells);
  outPd->SetPoints(fieldPts);
  outPd->BuildLinks();

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ProjectOntoUnitSphere(vtkPolyData *inPd,
                                                         vtkPolyData *outPd)
{
  int numPts = inPd->GetNumberOfPoints();
  double massCenter[3];

  vtkSVGeneralUtils::ComputeMassCenter(inPd, massCenter);

  vtkNew(vtkPoints, newPts);
  newPts->SetNumberOfPoints(numPts);
  vtkNew(vtkCellArray, newCells);
  newCells = inPd->GetPolys();

  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    double newpt[3];
    inPd->GetPoint(i, pt);
    for (int j=0; j<3; j++)
    {
      newpt[j] = pt[j] - massCenter[j];
    }
    vtkMath::Normalize(newpt);
    newPts->SetPoint(i, newpt);
  }

  outPd->SetPolys(newCells);
  outPd->SetPoints(newPts);
  outPd->BuildLinks();
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputeNormals(vtkPolyData *pd)
{
  vtkNew(vtkPolyDataNormals, normaler);
  normaler->SetInputData(pd);
  normaler->AutoOrientNormalsOn();
  normaler->SplittingOff();
  normaler->Update();

  pd->DeepCopy(normaler->GetOutput());
  pd->BuildLinks();

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputeMeshLaplacian(vtkPolyData *pd,
                                            vtkEdgeTable *edgeTable,
                                            vtkFloatArray *edgeWeights,
                                            vtkIntArray *edgeNeighbors,
                                            vtkFloatArray *laplacian, int map)
{
  int numPts = pd->GetNumberOfPoints();

  for (int i=0; i<numPts; i++)
  {
    double pointLaplacian[3];
    vtkSVGeneralUtils::ComputePointLaplacian(i, pd,
                                             edgeTable,
                                             edgeWeights, edgeNeighbors,
                                             pointLaplacian, map);
    laplacian->SetTuple(i, pointLaplacian);
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputePointLaplacian(vtkIdType p0,
                                             vtkPolyData *pd,
                                             vtkEdgeTable *edgeTable,
                                             vtkFloatArray *edgeWeights,
                                             vtkIntArray *edgeNeighbors,
                                             double laplacian[],
                                             int map)
{
  vtkNew(vtkIdList, pointNeighbors);
  vtkSVGeneralUtils::GetPointNeighbors(p0, pd, pointNeighbors);

  laplacian[0] = 0.0; laplacian[1] = 0.0; laplacian[2] = 0.0;
  for (int i=0; i<pointNeighbors->GetNumberOfIds(); i++)
  {
    vtkIdType p1 = pointNeighbors->GetId(i);
    vtkIdType edgeId = edgeTable->IsEdge(p0, p1);
    double weight = edgeWeights->GetValue(edgeId);
    //if (map == TUTTE)
    //{
    //  weight = 1.0;
    //}
    int edgeNeighbor = edgeNeighbors->GetValue(edgeId);
    if (edgeNeighbor == -1)
    {
      continue;
    }
    double p0Metric[3], p1Metric[3], data0[3], data1[3];
    pd->GetPoint(p0, p0Metric);
    pd->GetPoint(p1, p1Metric);

    for (int j=0; j<3; j++)
    {
      laplacian[j] += weight * (p0Metric[j] - p1Metric[j]);
    }
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputeDataLaplacian(vtkIdType p0,
                                            vtkFloatArray *data,
                                            vtkPolyData *pd,
                                            vtkEdgeTable *edgeTable,
                                            vtkFloatArray *edgeWeights,
                                            vtkIntArray *edgeNeighbors,
                                            double laplacian[],
                                            int map)
{
  vtkNew(vtkIdList, pointNeighbors);
  vtkNew(vtkIdList, dummyList);
  vtkSVGeneralUtils::GetPointNeighbors(p0, pd, pointNeighbors);

  laplacian[0] = 0.0; laplacian[1] = 0.0; laplacian[2] = 0.0;
  for (int i=0; i<pointNeighbors->GetNumberOfIds(); i++)
  {
    vtkIdType p1 = pointNeighbors->GetId(i);
    vtkIdType edgeId = edgeTable->IsEdge(p0, p1);
    double weight = edgeWeights->GetValue(edgeId);
    //if (map == TUTTE)
    //{
    //  weight = 1.0;
    //}
    int edgeNeighbor = edgeNeighbors->GetValue(edgeId);
    if (edgeNeighbor == -1)
    {
      continue;
    }
    double p0Metric[3], p1Metric[3], data0[3], data1[3];
    data->GetTuple(p0, p0Metric);
    data->GetTuple(p1, p1Metric);

    for (int j=0; j<3; j++)
    {
      laplacian[j] += weight * (p0Metric[j] - p1Metric[j]);
    }
  }

  return SV_OK;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ComputeDataArrayLaplacian(vtkFloatArray *data,
                                                 vtkPolyData *pd,
                                                 vtkEdgeTable *edgeTable,
                                                 vtkFloatArray *edgeWeights,
                                                 vtkIntArray *edgeNeighbors,
                                                 vtkFloatArray *laplacian, int map)
{
  int numPts = data->GetNumberOfTuples();

  for (int i=0; i<numPts; i++)
  {
    double pointLaplacian[3];
    vtkSVGeneralUtils::ComputeDataLaplacian(i, data, pd, edgeTable, edgeWeights,
                                                        edgeNeighbors, pointLaplacian, map);
    laplacian->SetTuple(i, pointLaplacian);
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::RunLoopFind(vtkPolyData *pd,
                                   vtkIdType startPt,
                                   vtkIdType nextCell,
                                   vtkPolyData *loop,
                                   vtkIdList *boundaryIds)
{
  int checkIds = 0;
  int checkNum = 0;
  vtkNew(vtkIntArray, checkList);
  checkList->SetNumberOfTuples(pd->GetNumberOfPoints());
  checkList->FillComponent(0, 0);
  if (boundaryIds != NULL)
  {
    checkIds = 1;
    for (int i=0; i<boundaryIds->GetNumberOfIds(); i++)
      checkList->SetTuple1(boundaryIds->GetId(i), 1);
    if (startPt != boundaryIds->GetId(0))
    {
      fprintf(stdout,"Start point does not match given\n");
      return SV_ERROR;
    }
    checkNum++;
  }

  vtkIdType prevPt = startPt;
  vtkIdType nextPt = startPt;
  vtkNew(vtkIdList, pointIds);
  vtkNew(vtkIdList, cellIds);

  pd->GetCellPoints(nextCell,pointIds);

  if (pointIds->GetId(0) == nextPt)
    nextPt = pointIds->GetId(1);
  else
    nextPt = pointIds->GetId(0);
  vtkNew(vtkIdList, newline);
  newline->SetNumberOfIds(2);
  newline->SetId(0, prevPt);
  newline->SetId(1, nextPt);
  //newline.id = nextCell;
  loop->InsertNextCell(VTK_LINE, newline);

  while(nextPt != startPt)
  {
    if (checkIds && checkList->GetTuple1(nextPt))
    {
      if (checkNum != boundaryIds->IsId(nextPt))
      {
        fprintf(stdout,"Boundary points are not in correct order\n");
        return SV_ERROR;
      }
      checkNum++;
    }
    pd->GetPointCells(nextPt,cellIds);
    if (cellIds->GetId(0) == nextCell)
      nextCell = cellIds->GetId(1);
    else
      nextCell = cellIds->GetId(0);

    pd->GetCellPoints(nextCell,pointIds);
    prevPt = nextPt;
    if (pointIds->GetId(0) == nextPt)
      nextPt = pointIds->GetId(1);
    else
      nextPt = pointIds->GetId(0);

    vtkNew(vtkIdList, newestline);
    newestline->SetNumberOfIds(2);
    newestline->InsertId(0, prevPt);
    newestline->InsertId(1, nextPt);
    //newestline.id = nextCell;
    loop->InsertNextCell(VTK_LINE, newestline);
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
//Determine type of intersection
int vtkSVGeneralUtils::SeparateLoops(vtkPolyData *pd,
                                     vtkPolyData **loops,
                                     int numBoundaries,
                                     const double xvec[3],
                                     const double zvec[3],
                                     const int boundaryStart[2])
{
  vtkIdType nextCell;
  vtkNew(vtkIdList, cellIds);
  int numInterPts = pd->GetNumberOfPoints();
  int numInterLines = pd->GetNumberOfLines();
  pd->BuildLinks();

  int count = 0;
  for (int i=0;i<numBoundaries;i++)
  {
    vtkIdType startPt = boundaryStart[i];
    vtkPolyData *newloop = loops[count];
    newloop->Allocate(pd->GetNumberOfCells(), 1000);
    pd->GetPointCells(startPt,cellIds);

    nextCell = cellIds->GetId(0);
    vtkIdType npts, *pts;
    int testPt = -1;
    pd->GetCellPoints(nextCell, npts, pts);
    if (pts[0] == startPt)
      testPt = pts[1];
    else
      testPt = pts[0];

    double pt0[3], pt1[3], vec0[3], vec1[3];
    pd->GetPoint(startPt, pt0);
    pd->GetPoint(testPt, pt1);
    vtkMath::Subtract(pt1, pt0, vec0);
    vtkMath::Normalize(vec0);
    vtkMath::Cross(zvec, xvec, vec1);
    vtkMath::Normalize(vec1);
    if (vtkMath::Dot(vec0, vec1) < 0)
    {
      nextCell = cellIds->GetId(1);
    }
    //if (testPt != boundaryStart[i+2])
    //{
    //  nextCell = cellIds->GetId(1);
    //}

    //Run through intersection lines to get loops!
    vtkSVGeneralUtils::RunLoopFind(pd, startPt, nextCell, newloop, NULL);
    loops[count++] = newloop;
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::VectorDotProduct(vtkFloatArray *v0, vtkFloatArray *v1, double product[], int numVals, int numComps)
{
  for (int i=0; i<numComps; i++)
  {
    product[i] = 0.0;
  }
  for (int i=0; i<numVals; i++)
  {
    for (int j=0; j<numComps; j++)
    {
      double val0, val1;
      val0 = v0->GetComponent(i, j);
      val1 = v1->GetComponent(i, j);
      product[j] += val0 * val1;
    }
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetRotationMatrix(double vec0[3], double vec1[3], vtkMatrix4x4 *rotMatrix)
{
  double perpVec[3];
  vtkMath::Normalize(vec0);
  vtkMath::Normalize(vec1);
  vtkMath::Cross(vec0, vec1, perpVec);
  double costheta = vtkMath::Dot(vec0, vec1);
  double sintheta = vtkMath::Norm(perpVec);
  double theta = atan2(sintheta, costheta);
  if (sintheta != 0)
  {
    perpVec[0] /= sintheta;
    perpVec[1] /= sintheta;
    perpVec[2] /= sintheta;
  }
  costheta = cos(0.5*theta);
  sintheta = sin(0.5*theta);
  double quat[4];
  quat[0] = costheta;
  quat[1] = perpVec[0]*sintheta;
  quat[2] = perpVec[1]*sintheta;
  quat[3] = perpVec[2]*sintheta;

  double mat[3][3];
  vtkMath::QuaternionToMatrix3x3(quat, mat);

  // | R_0 R_1 R_2 0 |
  // | R_3 R_4 R_2 0 |
  // | R_6 R_7 R_8 0 |
  // |  0   0   0  1 |
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      rotMatrix->SetElement(i, j, mat[i][j]);
    }
    rotMatrix->SetElement(i, 3, 0.0);
    rotMatrix->SetElement(3, i, 0.0);
  }
  rotMatrix->SetElement(3, 3, 1.0);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetRotationMatrix(double vec0[3], double vec1[3], double rotMatrix[16])
{
  vtkNew(vtkMatrix4x4, rotMatrix4x4);
  vtkSVGeneralUtils::GetRotationMatrix(vec0, vec1, rotMatrix4x4);
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
      rotMatrix[i*4+j] = rotMatrix4x4->GetElement(i, j);
  }

  return SV_OK;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ApplyRotationMatrix(vtkPolyData *pd, vtkMatrix4x4 *rotMatrix)
{
  vtkSmartPointer<vtkTransform> transformer =
    vtkSmartPointer<vtkTransform>::New();
  transformer->SetMatrix(rotMatrix);

  vtkSmartPointer<vtkTransformPolyDataFilter> pdTransformer =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  pdTransformer->SetInputData(pd);
  pdTransformer->SetTransform(transformer);
  pdTransformer->Update();

  pd->DeepCopy(pdTransformer->GetOutput());
  pd->BuildLinks();

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::ApplyRotationMatrix(vtkPolyData *pd, double rotMatrix[16])
{
  vtkSmartPointer<vtkTransform> transformer =
    vtkSmartPointer<vtkTransform>::New();
  transformer->SetMatrix(rotMatrix);

  vtkSmartPointer<vtkTransformPolyDataFilter> pdTransformer =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  pdTransformer->SetInputData(pd);
  pdTransformer->SetTransform(transformer);
  pdTransformer->Update();

  pd->DeepCopy(pdTransformer->GetOutput());
  pd->BuildLinks();
  return SV_OK;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::GetPolyDataAngles(vtkPolyData *pd, vtkFloatArray *cellAngles)
{
  int numCells = pd->GetNumberOfCells();
  pd->BuildLinks();

  cellAngles->SetNumberOfComponents(3);
  cellAngles->Allocate(numCells, 10000);
  cellAngles->SetNumberOfTuples(numCells);

  for (int i=0; i<numCells; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    for (int j=0; j<3; j++)
    {
      vtkIdType p0 = pts[j];
      vtkIdType p1 = pts[(j+1)%npts];
      vtkIdType p2 = pts[(j+2)%npts];

      double pt0[3], pt1[3], pt2[3];
      pd->GetPoint(p0, pt0);
      pd->GetPoint(p1, pt1);
      pd->GetPoint(p2, pt2);

      double vec0[3], vec1[3];
      for (int k=0; k<3; k++)
      {
        vec0[k] = pt0[k] - pt1[k];
        vec1[k] = pt2[k] - pt1[k];
      }

      double angleVec[3];
      vtkMath::Cross(vec0, vec1, angleVec);
      double radAngle = atan2(vtkMath::Norm(angleVec), vtkMath::Dot(vec0, vec1));

      cellAngles->SetComponent(i, j, radAngle);
    }
  }

  return SV_OK;
}



//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralUtils::VectorAdd(vtkFloatArray *v0, vtkFloatArray *v1, double scalar, vtkFloatArray *result, int numVals, int numComps)
{
  for (int i=0; i<numVals; i++)
  {
    for (int j=0; j<numComps; j++)
    {
      double val0, val1;
      val0 = v0->GetComponent(i, j);
      val1 = v1->GetComponent(i, j);
      double sum = val0 + scalar * val1;
      result->SetComponent(i, j, sum);
    }
  }

  return SV_OK;
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

  return SV_OK;
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

  return SV_OK;
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
  return SV_OK;
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
  return SV_OK;
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

  return SV_OK;
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

  return SV_OK;
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

  return SV_OK;
}

