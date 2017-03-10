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

/** @file vtkSVGeneralUtils.h
 *  @brief
 *  @brief
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef vtkSVGeneralUtils_h
#define vtkSVGeneralUtils_h

#include "vtkObject.h"

#include "vtkDataSet.h"
#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkImplicitFunction.h"
#include "vtkMatrix4x4.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPolyData.h"

#include <string>
#include <sstream>
#include <iostream>
#include <map>
#include <list>

class vtkSVGeneralUtils : public vtkObject
{
public:
  static vtkSVGeneralUtils *New();
  vtkTypeMacro(vtkSVGeneralUtils,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);

  //Checking functions
  static int CheckArrayExists(vtkDataSet *ds, int datatype, std::string arrayname);
  static int CheckSurface(vtkPolyData *pd);

  //General operations
  static int GetClosestPointConnectedRegion(vtkPolyData *pd,
                                            double origin[3]);
  static int GetClosestPointConnectedRegion(vtkPolyData *inPd,
                                            double origin[3],
                                            vtkPolyData *outPd);
  static int IteratePoint(vtkPolyData *pd, int &pointId, int &prevCellId);
  static int ThresholdPd(vtkPolyData *pd, int minVal, int maxVal, int dataType,
                         std::string arrayName, vtkPolyData *returnPd);
  static int GetCentroidOfPoints(vtkPoints *points, double centroid[3]);
  static int GetPointGroups(vtkPolyData *pd, std::string arrayName,
                            const int pointId, vtkIdList *groupIds);
  static int ExtractionCut(vtkPolyData *inPd, vtkImplicitFunction *cutFunction,
                           const int extractBoundaryCells,
                           const int extractInside,
                           vtkPolyData *outPd);
  static int ClipCut(vtkPolyData *inPd, vtkImplicitFunction *cutFunction,
                     const int generateClippedOutput,
                     const int extractInside,
                     vtkPolyData *outPd,
                     vtkPolyData *clippedOutPd);
  static double Distance(double pt0[3], double pt1[3]);
  static int GetPointsLength(vtkPolyData *points, double &length);
  static int ReplaceDataOnCells(vtkPointSet *pointset, vtkDataArray *sliceIds,
                                const int sliceId, const int replaceVal,
                                const std::string &arrName);
  static int GetCutPlane(const double endPt[3], const double startPt[3],
                         const double length, double origin[3], vtkPlane *cutPlane);
  static int ComputeTriangleArea(double pt0[], double pt1[], double pt2[], double &area);
  static int ComputeMassCenter(vtkPolyData *pd, double massCenter[3]);
  static int GetBarycentricCoordinates(double f[3], double pt0[3], double pt1[3],
                                       double pt2[3], double &a0, double &a1, double &a2);
  static int GetPointNeighbors(vtkIdType p0, vtkPolyData *pd, vtkIdList *pointNeighbors);
  static int GetEdgeCotangentAngle(double pt0[3], double pt1[3], double pt2[3], double &angle);
  static int CreateEdgeTable(vtkPolyData *pd, vtkEdgeTable *edgeTable,
                             vtkFloatArray *edgeWeights,
                             vtkIntArray *edgeNeighbors,
                             vtkIntArray *isBoundary);
  static int ComputeHarmonicEdgeWeight(vtkPolyData *pd, vtkIdType cellId,
                                       vtkIdType neighborCellId,
                                       vtkIdType p0, vtkIdType p1, double &weight);
  static int ConvertFieldToPolyData(vtkPolyData *inPd, std::string fieldName, vtkPolyData *outPd);
  static int ProjectOntoUnitSphere(vtkPolyData *inPd, vtkPolyData *outPd);
  static int ComputeNormals(vtkPolyData *pd);
  static int ComputeMeshLaplacian(vtkPolyData *pd, vtkEdgeTable *edgeTable,
                                  vtkFloatArray *edgeWeights, vtkIntArray *edgeNeighbors,
                                  vtkFloatArray *laplacian, int map);
  static int ComputeDataArrayLaplacian(vtkFloatArray *data, vtkPolyData *pd,
                                       vtkEdgeTable *edgeTable,
                                       vtkFloatArray *edgeWeights, vtkIntArray *edgeNeighbors,
                                       vtkFloatArray *laplacian, int map);
  static int ComputePointLaplacian(vtkIdType p0, vtkPolyData *pd,
                            vtkEdgeTable *edgeTable, vtkFloatArray *edgeWeights,
                            vtkIntArray *edgeNeighbors, double laplacian[], int map);
  static int ComputeDataLaplacian(vtkIdType p0, vtkFloatArray *data, vtkPolyData *pd,
                                  vtkEdgeTable *edgeTable, vtkFloatArray *edgeWeights,
                                  vtkIntArray *edgeNeighbors,
                                  double laplacian[], int map);
  static int RunLoopFind(vtkPolyData *pd, vtkIdType startPt, vtkIdType nextCell,
                         vtkPolyData *loop, vtkIdList *boundaryIds);
  static int SeparateLoops(vtkPolyData *pd, vtkPolyData **loops, int numBoundaries, const double xvec[3], const double zvec[3], const int boundaryStart[2]);
  static int VectorDotProduct(vtkFloatArray *v0, vtkFloatArray *v1, double product[], int numVals, int numComps);
  static int VectorAdd(vtkFloatArray *v0, vtkFloatArray *v1, double scalar, vtkFloatArray *result, int numVals, int numComps);
  static int GetRotationMatrix(double vec0[3], double vec1[3], vtkMatrix4x4 *rotMatrix);
  static int GetRotationMatrix(double vec0[3], double vec1[3], double rotMatrix[16]);
  static int ApplyRotationMatrix(vtkPolyData *pd, vtkMatrix4x4 *rotMatrix);
  static int ApplyRotationMatrix(vtkPolyData *pd, double rotMatrix[16]);
  static int GetPolyDataAngles(vtkPolyData *pd, vtkFloatArray *cellAngles);


  //std::map functions
  static int GetAllMapKeys(std::multimap<int, int> &map, std::list<int> &list);
  static int GetAllMapValues(std::multimap<int, int> &map, std::list<int> &list);
  static int GetValuesFromMap(std::multimap<int, int> &map, const int key, std::list<int> &list);
  static int GetKeysFromMap(std::multimap<int, int> &map, const int value, std::list<int> &list);
  static int GetCommonValues(std::multimap<int, int> &map, const int keyA,
                             const int keyB, std::list<int> &returnList);
  static int GetUniqueNeighbors(std::multimap<int, int> &map, const int key, std::list<int> &keyVals,
                                std::list<int> &uniqueKeys);
  static int ListIntersection(std::list<int> &listA,
                              std::list<int> &listB,
                              std::list<int> &returnList);

protected:
  vtkSVGeneralUtils();
  ~vtkSVGeneralUtils();

private:
  vtkSVGeneralUtils(const vtkSVGeneralUtils&);  // Not implemented.
  void operator=(const vtkSVGeneralUtils&);  // Not implemented.
};

#endif
