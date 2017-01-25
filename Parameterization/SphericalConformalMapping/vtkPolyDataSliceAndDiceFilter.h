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


/** @file vtkPolyDataSliceAndDiceFilter.h
 *  @brief This is a vtk filter to map a triangulated surface to a sphere.
 *  @details This filter uses the heat flow method to map a triangulated
 *  surface to a sphere. The first step is to compute the Tutte Energy, and
 *  the second step is to perform the conformal map. For more details, see
 *  Gu et al., Genus Zero Surface Conformal Mapping and Its
 *  Application to Brain Surface Mapping, 2004.
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef vtkPolyDataSliceAndDiceFilter_h
#define vtkPolyDataSliceAndDiceFilter_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkGeneralizedPolycube.h"
#include "vtkPlane.h"
#include "vtkPolyData.h"

#include "svGraph.h"

#include <map>
#include <list>

class vtkPolyDataSliceAndDiceFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkPolyDataSliceAndDiceFilter* New();
  vtkTypeRevisionMacro(vtkPolyDataSliceAndDiceFilter, vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Print statements used for debugging
  vtkGetMacro(Verbose, int);
  vtkSetMacro(Verbose, int);

  // Description:
  // For slice length
  vtkGetMacro(SliceLength, double);
  vtkSetMacro(SliceLength, double);

  // Description:
  // If this flag is on, the a polycube will be constructed based on the
  // separation of the polydata object
  vtkGetMacro(ConstructPolycube, int);
  vtkSetMacro(ConstructPolycube, int);

  // Description:
  // array names used in the filter
  vtkGetStringMacro(SliceIdsArrayName);
  vtkSetStringMacro(SliceIdsArrayName);
  vtkGetStringMacro(GroupIdsArrayName);
  vtkSetStringMacro(GroupIdsArrayName);
  vtkGetStringMacro(SegmentIdsArrayName);
  vtkSetStringMacro(SegmentIdsArrayName);
  vtkGetStringMacro(SphereRadiusArrayName);
  vtkSetStringMacro(SphereRadiusArrayName);
  vtkGetStringMacro(BoundaryPointsArrayName);
  vtkSetStringMacro(BoundaryPointsArrayName);
  vtkGetStringMacro(InternalIdsArrayName);
  vtkSetStringMacro(InternalIdsArrayName);
  vtkGetStringMacro(DijkstraArrayName);
  vtkSetStringMacro(DijkstraArrayName);

  // Description:
  // Macro to set object centerlines
  vtkGetObjectMacro(Centerlines, vtkPolyData);
  vtkSetObjectMacro(Centerlines, vtkPolyData);

  // Description:
  // Macro to set object polycube; will only be populated if ConstructPolycube is turned on
  vtkGetObjectMacro(Polycube, vtkGeneralizedPolycube);
  vtkSetObjectMacro(Polycube, vtkGeneralizedPolycube);

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
  static int IteratePoint(vtkPolyData *pd, int &pointId, int &prevCellId);
  static int ThresholdPd(vtkPolyData *pd, int minVal, int maxVal, int dataType,
                         std::string arrayName, vtkPolyData *returnPd);
  static int GetCentroidOfPoints(vtkPoints *points, double centroid[3]);
  static int GetClosestPointConnectedRegion(vtkPolyData *inPd,
                                            double origin[3],
                                            vtkPolyData *outPd);


  enum DIRECTIONS
  {
    RIGHT = 0,
    LEFT,
    FRONT,
    BACK,
    UP,
    DOWN
  };

protected:
  vtkPolyDataSliceAndDiceFilter();
  ~vtkPolyDataSliceAndDiceFilter();

  vtkSetStringMacro(NumberOfCubesArrayName);

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  int ComputeCenterlines();
  int ExtractBranches();
  int PreProcessPolyData();
  int FindGroupBoundaries();
  int BuildSkeleton();
  int BuildPolycube();
  void CheckLength(int &ptId, const int numPts,
                  int &done);
  void UpdatePtId(int &ptId);
  int SetDir(const int dir, double newDir[3]);
  int FormDirectionTable(int dirTable[6][4]);
  int CheckStartSurgeryPoints(vtkPolyData *pd, vtkIdList *startPoints);
  int StarterSurgeryPoints(vtkPolyData *pd,
                           const int groupId,
                           const int frontId,
                           const int backId,
                           double startPt[3],
                           double secondPt[3],
                           vtkIdList *fixedSurgeryPoints);
  int GetSurgeryPoints(vtkPolyData *pd,
                       vtkDataArray *pointIds,
                       const double clStartPt[3],
                       const double clSecondPt[3],
                       const int front,
                       const int back,
                       vtkIdList *surgeryPoints);
  int GetHalfSurgeryPoints(vtkPolyData *pd,
                           vtkDataArray *pointIds,
                           const int cellId,
                           const int front,
                           const int back,
                           vtkIdList *surgeryPoints);
  int DetermineSliceStrategy(vtkPolyData *branchPd,
                             const int branchId,
                             vtkPolyData *branchCenterline,
                             int &branchStartPtId,
                             vtkIdList *surgeryPoints,
                             int &centerlineStartPtId);
  int GetPointGroups(vtkPolyData *pd, std::string arrayName,
                     const int pointId, vtkIdList *groupIds);
  int GetCriticalPoints();
  int InsertCriticalPoints(const int pointId, vtkIdList *groupIds);
  int GetBranch(const int branchId, vtkPolyData *branchPd,
                                    vtkPolyData *branchCenterlines);
  int SliceBranches();
  int SliceBranch(vtkPolyData *branchPd, vtkPolyData *branchCenterline,
                    const int branchId,
                    vtkDataArray *sliceIds);
  int SliceBifurcations();
  int SliceBifurcation(vtkPolyData *pd,
                       svGCell *gCell,
                       int &segmentId,
                       vtkDataArray *segmentIds);
  int GetRegionsOnPd(vtkPolyData *pd,
                     vtkIdList *regionIds);
  int GetNumberOfSlicesBefore(const int id, int &numBefore);
  bool ListsMatch(vtkIdList *listA, vtkIdList *listB);
  int GetSectionZAxis(const double endPt[3], const double startPt[3],
                      double zvec[3]);
  int GetSectionXAxis(const double endPt[3], const double startPt[3],
                      const double surfacePt[3], double xvec[3]);
  int GetCutPlane(const double endPt[3], const double startPt[3],
                  const double length, double origin[3], vtkPlane *cutPlane);
  int ExtractionCut(vtkPolyData *inPd, vtkPlane *cutPlane,
                    const int extractBoundaryCells,
                    const int extractInside,
                    vtkPolyData *outPd);

  int ReplaceDataOnCells(vtkPointSet *pointset, vtkDataArray *sliceIds,
                         const int sliceId, const int replaceVal,
                         const std::string &arrName);
  int GetBranchLength(vtkPolyData *points, double &length);
  int GetCloseGeodesicPoint(vtkPolyData *pd, double centerPt[3], const int startPtId, int &returnStartId, double zvec[3], vtkPolyData *boundary);
  int GetNextSurgeryPoints(vtkPolyData *pd, double centerPt[3], vtkIdList *surgeryPoints, double xvec[3], double zvec[3]);
  int GetClose3DPoint(vtkPolyData *pd, double centerPt[3], const int startPtId, int &returnStartId, double zvec[3]);
  int GetContourSecondPoint(vtkPolyData *pd, int ptId, double centerPt[3], double zvec[3], int &startSecondId);

private:
  vtkPolyDataSliceAndDiceFilter(const vtkPolyDataSliceAndDiceFilter&);  // Not implemented.
  void operator=(const vtkPolyDataSliceAndDiceFilter&);  // Not implemented.

  int Verbose;

  char *BoundaryPointsArrayName;
  char *SliceIdsArrayName;
  char *SphereRadiusArrayName;
  char *GroupIdsArrayName;
  char *SegmentIdsArrayName;
  char *InternalIdsArrayName;
  char *DijkstraArrayName;
  char *NumberOfCubesArrayName;

  vtkPolyData    *InitialPd;
  vtkPolyData    *WorkPd;
  vtkPolyData    *GraphPd;
  vtkPolyData    *Centerlines;
  svGraph *CenterlineGraph;

  std::multimap<int , int> CriticalPointMap;
  std::multimap<int , int> SurgeryPointMap;

  int ConstructPolycube;
  int SliceDirection;
  int TotalSliceId;
  int MaxGroupNumber;
  int MaxSegmentNumber;
  int DirectionTable[6][4];

  double FirstBranchVec[3];
  double SliceLength;

  vtkIntArray    *StartPtIds;
  vtkDoubleArray *TopNormals;
  vtkDoubleArray *RightNormals;
  vtkDoubleArray *BifTopNormals;
  vtkDoubleArray *BifRightNormals;
  vtkIntArray    *NumberOfCubes;

  vtkGeneralizedPolycube *Polycube;
};

#endif
