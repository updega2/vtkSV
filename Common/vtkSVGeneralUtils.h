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
  /** \brief Function to check is array with name exists in cell or point data
   *  \param ds this is the object to check if the array exists
   *  \param datatype this is point or cell. point=0,cell=1
   *  \param arrayname this is the name of the array to check
   *  \reutrn this returns 1 if the array exists and zero if it doesn't
   *  or the function does not return properly. */
  static int CheckArrayExists(vtkDataSet *ds, int datatype, std::string arrayname);

  /** \brief
   *  \param polydata to check. */
  static int CheckSurface(vtkPolyData *pd);

  //General operations
  /** \brief Gets region closest to given point.
   *  \param pd Input pd which is copied in place with result.
   *  \param pt The closest point connected region.
   *  \return SV_OK */
  static int GetClosestPointConnectedRegion(vtkPolyData *pd,
                                            double pt[3]);

  /** \brief Gets region closest to given point.
   *  \param inPd Input pd.
   *  \param pt The closest point connected region.
   *  \param outPd Result pd.
   *  \return SV_OK */
  static int GetClosestPointConnectedRegion(vtkPolyData *inPd,
                                            double pt[3],
                                            vtkPolyData *outPd);
  /** \brief Sets ids on pd using vtkIdFilter.
   *  \param pd Inpud pd which is copied in place
   *  \param arrayName Name of array to give to ids
   *  \return SV_OK */
  static int GiveIds(vtkPolyData *pd,
                     std::string arrayName);

  /** \brief Sets ids on pd using vtkIdFilter.
   *  \param inPd Input pd.
   *  \param arrayName Name of array to give to ids
   *  \param outPd Result pd.
   *  \return SV_OK */
  static int GiveIds(vtkPolyData *inPd,
                     std::string arrayName,
                     vtkPolyData *outPd);

  /** \brief Function to iterate a point in a VTK_LINE or VTK_POLY_LINE
   *  \param pd Input pd containing line cells.
   *  \param pointId current point id; will be updated in function.
   *  \param prevCellId previous cell id; will be updated in function.
   *  \return SV_OK */
  static int IteratePoint(vtkPolyData *pd, int &pointId, int &prevCellId);

  /** \brief Awesome function that is wrapper around vtkThreshold
   *  \param pd Input pd to threshold, updated in place.
   *  \param minVal minumum value.
   *  \param maxVal maximum value.
   *  \param dataType 0 for point data, 1 for cell data.
   *  \param arrayName Name of array to be used to threshold polydata.
   *  \return SV_OK */
  static int ThresholdPd(vtkPolyData *pd, int minVal, int maxVal, int dataType,
                         std::string arrayName);

  /** \brief Awesome function that is wrapper around vtkThreshold
   *  \param pd Input pd to threshold.
   *  \param minVal minumum value.
   *  \param maxVal maximum value.
   *  \param dataType 0 for point data, 1 for cell data.
   *  \param arrayName Name of array to be used to threshold polydata.
   *  \param returnPd The resultant thresholded pd
   *  \return SV_OK */
  static int ThresholdPd(vtkPolyData *pd, int minVal, int maxVal, int dataType,
                         std::string arrayName, vtkPolyData *returnPd);

  /** \brief Get centroid of points */
  static int GetCentroidOfPoints(vtkPoints *points, double centroid[3]);

  /** \brief For an array of integers on pd, will return a list of values on cells
   *  attached to point.
   *  \param pd The full polydata.
   *  \param arrayName Name of array to get cell data of.
   *  \param pointId The point Id to get array values of.
   *  \param groupIds list of values on cells attached to point.
   *  \return SV_OK */
  static int GetPointCellsValues(vtkPolyData *pd, std::string arrayName,
                                const int pointId, vtkIdList *valList);

  /** \brief Perform a cut of the polydata, but with crinkle clip. Uses vtkExtractGeometry
   *  \param inPd The pd to cut.
   *  \param cutFunction implicit function defining cut plane/box/etc.
   *  \param extractBoundaryCells If 1, cells that are slice through are extracted.
   *  \param extractInside If 1, flips the function and extracts otherside of pd.
   *  \param outPd Result of the cut operation.
   *  \return SV_OK */
  static int ExtractionCut(vtkPolyData *inPd, vtkImplicitFunction *cutFunction,
                           const int extractBoundaryCells,
                           const int extractInside,
                           vtkPolyData *outPd);

  /** \brief Perform a cut of the polydata.
   *  \param inPd The pd to cut.
   *  \param cutFunction implicit function defining cut plane/box/etc.
   *  \param extractInside If 1, flips the function and extracts otherside of pd.
   *  \param clippedOutPd if generateClippedOutput is 1 and pd is provided, the
   *  surface that was clipped out is returned with this parameter.
   *  \param outPd Result of the cut operation.
   *  \return SV_OK */
  static int ClipCut(vtkPolyData *inPd, vtkImplicitFunction *cutFunction,
                     const int generateClippedOutput,
                     const int extractInside,
                     vtkPolyData *outPd,
                     vtkPolyData *clippedOutPd);

  /** \brief Get distance between two 3D points, very simple
   *  \return the unsigned distance */
  static double Distance(double pt0[3], double pt1[3]);

  /** \brief For a polydata with cells of VTK_LINE of VTK_POLY_LINE, will calculate
   *  the total distance of that line.
   *  \param pd The polydata.
   *  \return the unsigned distance. */
  static double GetPointsLength(vtkPolyData *pd);

  /** \brief Given a pointset and an array, we populate the array with a value,
   *  but only if the current value has a certain value.
   *  \param pointset The given dataset.
   *  \param sliceIds The data array that we want to replace values in (cell data).
   *  \param sliceId The value to be added to array.
   *  \param replaceVal We only replace when the current value is equal to this!
   *  \param arrName If this is part of a larger dataset, this is used to link
   *  back to full dataset. If not, just run dataset through vtkIdFilter, and
   *  provide ids array name.
   *  \return SV_OK if function completes withouth error. */
  static int ReplaceDataOnCells(vtkPointSet *pointset, vtkDataArray *sliceIds,
                                const int sliceId, const int replaceVal,
                                const std::string &arrName);

  /** \brief Given two points, return the cutPlane perpendicular to their vector
   *  \param endPt End point of vector.
   *  \param startPt Start point of vector.
   *  \param cutPlane The plane defined by vector.
   *  \return SV_OK. */
  static int GetCutPlane(double endPt[3], double startPt[3],
                         vtkPlane *cutPlane);

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
  template <typename T, size_t nr, size_t nc>
  static void PrintArray(T (&array)[nr][nc])
  {
    std::cout << "Array: " << nr << "by " << nc << endl;
    std::cout << "-----------------------------------------------------------" << endl;
    for (int i=0; i<nr; i++)
    {
      for (int j=0; j<nc; j++)
      {
        std::cout << "| " << array[i][j] << " | ";
      }
      std::cout << endl;
    }
    std::cout << "-----------------------------------------------------------" << endl;
  };

protected:
  vtkSVGeneralUtils();
  ~vtkSVGeneralUtils();

private:
  vtkSVGeneralUtils(const vtkSVGeneralUtils&);  // Not implemented.
  void operator=(const vtkSVGeneralUtils&);  // Not implemented.
};

#endif
