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

/**
 *  \class vtkSVSurfaceCuboidPatcher
 *  \brief Using a polydata centerlines, separate the polydata into regions
 *  based on the centerlines
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVSurfaceCuboidPatcher_h
#define vtkSVSurfaceCuboidPatcher_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkSVPolyBallLine.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIdList.h"
#include "vtkMatrix4x4.h"

#include "vtkSVGlobals.h"

#include "vtkSVSegmentationModule.h" // For export

class VTKSVSEGMENTATION_EXPORT vtkSVSurfaceCuboidPatcher : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSVSurfaceCuboidPatcher,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSVSurfaceCuboidPatcher *New();

  //@{
  /// \brief Get/Set macro for merged centerlines
  vtkSetObjectMacro(MergedCenterlines,vtkPolyData);
  vtkGetObjectMacro(MergedCenterlines,vtkPolyData);
  //@}

  //@{
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(PolycubePd,vtkPolyData);
  vtkGetObjectMacro(PolycubePd,vtkPolyData);
  //@}

  //@{
  /// \brief Get/Set macro for array name used by the filter. Must
  //  be present on the centerlines.
  vtkSetStringMacro(CenterlineGroupIdsArrayName);
  vtkGetStringMacro(CenterlineGroupIdsArrayName);
  vtkSetStringMacro(CenterlineRadiusArrayName);
  vtkGetStringMacro(CenterlineRadiusArrayName);
  vtkSetStringMacro(GroupIdsArrayName);
  vtkGetStringMacro(GroupIdsArrayName);
  vtkSetStringMacro(BlankingArrayName);
  vtkGetStringMacro(BlankingArrayName);
  vtkSetStringMacro(CenterlineIdsArrayName);
  vtkGetStringMacro(CenterlineIdsArrayName);
  vtkSetStringMacro(TractIdsArrayName);
  vtkGetStringMacro(TractIdsArrayName);
  vtkSetStringMacro(PatchIdsArrayName);
  vtkGetStringMacro(PatchIdsArrayName);
  vtkSetStringMacro(SlicePointsArrayName);
  vtkGetStringMacro(SlicePointsArrayName);
  vtkSetStringMacro(ClusteringVectorArrayName);
  vtkGetStringMacro(ClusteringVectorArrayName);
  //@}

  //@{
  /// \brief Get/Set the radius information
  vtkSetMacro(EnforcePolycubeConnectivity, int);
  vtkGetMacro(EnforcePolycubeConnectivity, int);
  vtkBooleanMacro(EnforcePolycubeConnectivity, int);
  //@}

  /** \brief Correct cells on the boundary by updating val if they have
   *  multiple neighboring cells of the same value */
  static int CorrectSpecificCellBoundaries(vtkPolyData *pd, std::string cellArrayName,
                                           vtkIdList *targetRegions);

  /** \brief Naive implementation to get most reoccuring number in list. Okay
   *  because list size is small. */
  static void GetMostOccuringVal(vtkIdList *idList, int &output, int &max_count);

  static int SmoothSpecificBoundaries(vtkPolyData *pd, std::string arrayName,
                                      vtkIdList *targetRegions);

  static int GetRegions(vtkPolyData *pd, std::string arrayName,
                        std::vector<Region> &allRegions);
  static int GetSpecificRegions(vtkPolyData *pd, std::string arrayName,
                                std::vector<Region> &allRegions,
                                vtkIdList *targetRegions);

  static int CurveFitBoundaries(vtkPolyData *pd, std::string arrayName,
                                std::vector<Region> allRegions);

  static int GetCCWPoint(vtkPolyData *pd, const int pointId, const int cellId);
  static int GetCWPoint(vtkPolyData *pd, const int pointId, const int cellId);

  static int CheckBoundaryEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);

  static int GetPointEdgeCells(vtkPolyData *pd, std::string arrayName,
                               const int cellId, const int pointId,
                               vtkIdList *sameCells);

  static int CheckCellValuesEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);
  static void SplineKnots(std::vector<int> &u, int n, int t);

  static void SplineCurve(const std::vector<XYZ> &inp, int n, const std::vector<int> &knots, int t, std::vector<XYZ> &outp, int res);

  static void SplinePoint(const std::vector<int> &u, int n, int t, double v, const std::vector<XYZ> &control, XYZ &output);

  static double SplineBlend(int k, int t, const std::vector<int> &u, double v);

  static int FindPointMatchingValues(vtkPointSet *ps, std::string arrayName, vtkIdList *matchingVals, int &returnPtId);

protected:
  vtkSVSurfaceCuboidPatcher();
  ~vtkSVSurfaceCuboidPatcher();

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  int PrepFilter(); // Prep work.
  int RunFilter(); // Run filter operations.

  /** \brief Cluster branch */
  int CheckGroupsWithPolycube();
  int ClusterBranchWithCVT(vtkPolyData *pd, vtkPolyData *generatorPd);
  int ClusterBranchWithGeodesics(vtkPolyData *pd, vtkPolyData *polyPd);

  /** \brief Run basic edge weithing cvt with pd */
  int RunEdgeWeightedCVT(vtkPolyData *pd, vtkPolyData *generatorPd);

  int FixEndPatches(vtkPolyData *pd);
  int MatchEndPatches(vtkPolyData *branchPd, vtkPolyData *polyBranchPd);
  int CheckEndPatches(vtkPolyData *pd,
                      std::vector<Region> endRegions,
                      std::vector<int> &individualFix,
                      std::vector<int> &wholePatchFix);
  int FixSidePatches(vtkPolyData *pd);
  int CheckSidePatches(vtkPolyData *pd,
                      std::vector<Region> endRegions,
                      std::vector<int> &wholePatchFix);
  int GetOpenBoundaryEdges(vtkPolyData *branchPd,
                           std::vector<int> &openCornerPoints,
                           std::vector<std::vector<int> > &openEdges);
  int GetOpenBoundaryEdges(vtkPolyData *branchPd, std::vector<Region> regions,
                           std::string arrayName,
                           std::vector<int> &openCornerPoints,
                           std::vector<std::vector<int> > &openEdges);
  int ShiftEdgeList(vtkPolyData *branchPd, std::vector<std::vector<int> > &openEdges,
                    std::vector<std::vector<int> > &shiftedOpenEdges);
  int SplitEdgeList(vtkPolyData *branchPd, std::vector<int> &openEdges,
                    std::vector<std::vector<int> > &shiftedOpenEdges);

  int GetConnectedEdges(std::vector<std::vector<int> > inputEdges,
                        std::vector<std::vector<int> > &connectedCornerPts);
  int FixPatchesWithPolycube();

  char *CenterlineGroupIdsArrayName;
  char *CenterlineRadiusArrayName;
  char *CenterlineIdsArrayName;
  char *GroupIdsArrayName;
  char *BlankingArrayName;
  char *TractIdsArrayName;
  char *PatchIdsArrayName;
  char *SlicePointsArrayName;
  char *ClusteringVectorArrayName;

  vtkPolyData *WorkPd;
  vtkPolyData *MergedCenterlines;
  vtkPolyData *PolycubePd;

  int EnforcePolycubeConnectivity;

private:
  vtkSVSurfaceCuboidPatcher(const vtkSVSurfaceCuboidPatcher&);  // Not implemented.
  void operator=(const vtkSVSurfaceCuboidPatcher&);  // Not implemented.
};

#endif
