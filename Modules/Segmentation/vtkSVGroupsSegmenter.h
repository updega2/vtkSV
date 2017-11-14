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
 *  \class vtkSVGroupsSegmenter
 *  \brief Using a polydata centerlines, separate the polydata into regions
 *  based on the centerlines
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVGroupsSegmenter_h
#define vtkSVGroupsSegmenter_h

#include "vtkPolyDataAlgorithm.h"
#include "vtkPolyData.h"
#include "vtkSVPolyBallLine.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkIdList.h"
#include "vtkSVCenterlineGraph.h"
#include "vtkMatrix4x4.h"

#include "vtkSVSegmentationModule.h" // For export

struct Region
{
  int Index;
  int IndexCluster;

  int NumberOfCorners;
  std::vector<int> CornerPoints;

  std::vector<std::vector<int> > BoundaryEdges;

  int NumberOfElements;
  std::vector<int> Elements;

};
struct XYZ
{
  double x;
  double y;
  double z;
};

class VTKSVSEGMENTATION_EXPORT vtkSVGroupsSegmenter : public vtkPolyDataAlgorithm
{
public:
  vtkTypeMacro(vtkSVGroupsSegmenter,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkSVGroupsSegmenter *New();

  //@{
  /// \brief Get/Set macro for the object's centerlines
  vtkSetObjectMacro(Centerlines,vtkPolyData);
  vtkGetObjectMacro(Centerlines,vtkPolyData);
  //@}

  //@{
  /// \brief Get/Set macro for centerline group ids to use in the filter
  vtkSetObjectMacro(CenterlineGroupIds,vtkIdList);
  vtkGetObjectMacro(CenterlineGroupIds,vtkIdList);
  //@}

  //@{
  /// \brief Get the graph for the model
  vtkSetObjectMacro(GraphPd,vtkPolyData);
  vtkGetObjectMacro(GraphPd,vtkPolyData);
  //@}

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
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(PolycubeUg,vtkUnstructuredGrid);
  vtkGetObjectMacro(PolycubeUg,vtkUnstructuredGrid);
  //@}

  //@{
  /// \brief Get/Set macro for surface polycube
  vtkSetObjectMacro(FinalHexMesh,vtkUnstructuredGrid);
  vtkGetObjectMacro(FinalHexMesh,vtkUnstructuredGrid);
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
  //@}

  //@{
  /// \brief Get/Set/Boolean macro to indicate whether to clip all centerline
  /// group ids
  vtkSetMacro(ClipAllCenterlineGroupIds,int);
  vtkGetMacro(ClipAllCenterlineGroupIds,int);
  vtkBooleanMacro(ClipAllCenterlineGroupIds,int);
  //@}

  //@{
  /// \brief Get/Set the cutoff radius factor for clipping of the surface
  //  distance functions
  vtkSetMacro(CutoffRadiusFactor,double);
  vtkGetMacro(CutoffRadiusFactor,double);
  //@}

  //@{
  /// \brief Get/Set the clip value for clipping of the surface distance functions.
  vtkSetMacro(ClipValue,double);
  vtkGetMacro(ClipValue,double);
  //@}

  //@{
  /// \brief Get/Set the radius information
  vtkSetMacro(UseRadiusInformation,int);
  vtkGetMacro(UseRadiusInformation,int);
  vtkBooleanMacro(UseRadiusInformation,int);
  //@}

  //@{
  /// \brief Get/Set the initial group clipper to use
  vtkSetMacro(UseVmtkClipping,int);
  vtkGetMacro(UseVmtkClipping,int);
  vtkBooleanMacro(UseVmtkClipping,int);
  //@}

  //@{
  /// \brief Get/Set whether the boundary at separating patches should be more
  //  strictly enforced.
  vtkSetMacro(EnforceBoundaryDirections,int);
  vtkGetMacro(EnforceBoundaryDirections,int);
  vtkBooleanMacro(EnforceBoundaryDirections,int);
  //@}

  //@{
  /// \brief Get/Set the number of divisions to use along width and height of polycube
  vtkSetMacro(PolycubeDivisions,int);
  vtkGetMacro(PolycubeDivisions,int);
  //@}

  //@{
  /// \brief Get/Set the unit length for each division of the polycube
  vtkSetMacro(PolycubeUnitLength,double);
  vtkGetMacro(PolycubeUnitLength,double);
  //@}

  //@{
  /// \brief Get/Set the scalar determing how much influence to put on the normal
  // of the cell and how much influence to put on the position of the cell for
  // the cube patch clustering.
  vtkSetMacro(NormalsWeighting,double);
  vtkGetMacro(NormalsWeighting,double);
  //@}

  //@{
  /// \brief Get/Set whether the model is a vascular model with artificial truncated
  //  boundaries
  vtkSetMacro(IsVasculature,int);
  vtkGetMacro(IsVasculature,int);
  vtkBooleanMacro(IsVasculature,int);
  //@}

  //@{
  /// \brief Get/Set If model is not vasculature, indicate how many centerline
  //  points to remove from the ends
  vtkSetMacro(NumberOfCenterlineRemovePts,int);
  vtkGetMacro(NumberOfCenterlineRemovePts,int);
  //@}

  //@{
  /// \brief Get/Set whether centerlines should be modified based on a radius
  //  threshold. Useful if lots of branches close that vmtk clusters all into
  //  one separation point.
  vtkSetMacro(ModifyCenterlines,int);
  vtkGetMacro(ModifyCenterlines,int);
  vtkBooleanMacro(ModifyCenterlines,int);
  //@}

  //@{
  /// \brief Get/Set the radius threshold at which to determine a new branch
  //  if the ModifyCenterlines flag is set to 1.
  vtkSetMacro(CenterlineSeparationThreshold,double);
  vtkGetMacro(CenterlineSeparationThreshold,double);
  //@}

  /** \brief Correct cells on the boundary by updating val if they have
   *  multiple neighboring cells of the same value */
  static int CorrectCellBoundaries(vtkPolyData *pd, std::string cellArrayName);
  static int CorrectSpecificCellBoundaries(vtkPolyData *pd, std::string cellArrayName,
                                           vtkIdList *targetRegions);

  /** \brief Naive implementation to get most reoccuring number in list. Okay
   *  because list size is small. */
  static void GetMostOccuringVal(vtkIdList *idList, int &output, int &max_count);

  /** \brief Run basic edge weithing cvt with pd */
  static int RunEdgeWeightedCVT(vtkPolyData *pd, vtkPolyData *generatorPd);

  /** \brief From three vectors, compute transformation from global to local */
  static int ComputeRotationMatrix(const double vx[3], const double vy[3],
                                   const double vz[3], double rotMatrix[9]);

  static int SmoothBoundaries(vtkPolyData *pd, std::string arrayName);
  static int SmoothSpecificBoundaries(vtkPolyData *pd, std::string arrayName,
                                      vtkIdList *targetRegions);

  static int GetRegions(vtkPolyData *pd, std::string arrayName,
                        std::vector<Region> &allRegions);
  static int GetSpecificRegions(vtkPolyData *pd, std::string arrayName,
                                std::vector<Region> &allRegions,
                                vtkIdList *targetRegions);

  static int CurveFitBoundaries(vtkPolyData *pd, std::string arrayName,
                                std::vector<Region> allRegions);

  static int SplitBoundary(vtkPolyData *pd, std::vector<int> boundary,
                           int numDivs, int groupId, std::vector<int> &newSlicePoints);

  static int MatchBoundaries(vtkPolyData *pd, std::string checkArrayName,
                             std::string fitArrayName,
                             std::vector<Region> allRegions);

  static int GetCCWPoint(vtkPolyData *pd, const int pointId, const int cellId);
  static int GetCWPoint(vtkPolyData *pd, const int pointId, const int cellId);

  static int CheckBoundaryEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);

  static int CheckCellValuesEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);
  static void SplineKnots(std::vector<int> &u, int n, int t);

  static void SplineCurve(const std::vector<XYZ> &inp, int n, const std::vector<int> &knots, int t, std::vector<XYZ> &outp, int res);

  static void SplinePoint(const std::vector<int> &u, int n, int t, double v, const std::vector<XYZ> &control, XYZ &output);

  static double SplineBlend(int k, int t, const std::vector<int> &u, double v);

  static int FindPointMatchingValues(vtkPointSet *ps, std::string arrayName, vtkIdList *matchingVals, int &returnPtId);

  static int RotateGroupToGlobalAxis(vtkPolyData *pd,
                                     const int thresholdId,
                                     std::string arrayName,
                                     vtkPolyData *rotPd,
                                     vtkMatrix4x4 *rotMatrix0,
                                     vtkMatrix4x4 *rotMatrix1);
  static int InterpolateMapOntoTarget(vtkPolyData *sourceBasePd,
                                      vtkPolyData *targetPd,
                                      vtkPolyData *targetBasePd,
                                      vtkPolyData *mappedPd,
                                      std::string dataMatchingArrayName);

  const static double GlobalCoords[3][3];


protected:
  vtkSVGroupsSegmenter();
  ~vtkSVGroupsSegmenter();

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  int PrepFilter(); // Prep work.
  int RunFilter(); // Run filter operations.

  int MergeCenterlines();
  int FixCenterlines();
  int CheckCenterlineChildren(vtkPolyData *centerlinesPd,
                              double distPt[3],
                              double newPt[3],
                              std::vector<int> &testChildren,
                              std::vector<int> &newPointIds,
                              std::vector<int> &newChildren);
  int GetPatches();
  int MatchSurfaceToPolycube();
  int CheckSlicePoints();
  int GetApproximatePolycubeSize(double &polycubeSize);
  int SplitCellsAroundPoint(vtkPolyData *pd, int ptId);
  int SplitEdge(vtkPolyData *pd, int cellId, int ptId0, int ptId1,
                vtkCellArray *newCells, std::vector<std::vector<int> >  &splitCells);
  int FixMultipleGroups(vtkPolyData *pd, vtkPolyData *polycubePd,
                        std::vector<Region> surfaceGroups,
                        std::vector<Region> polycubeGroups);
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
  int GetTrueBoundaryDirections(vtkPolyData *branchPd,
                                vtkPolyData *polyBranchPd,
                                const int groupId,
                                vtkSVPolyBallLine *groupTubes,
                                std::vector<std::vector<int> > &shiftedOpenEdges,
                                vtkDoubleArray *avgVecs,
                                vtkIntArray *patchDirs);
  int CheckGroups();
  int FixEdges(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
               const Region region, std::vector<int> allEdges,
               std::vector<int> fixEdges, vtkIdList *critPts);
  int FixPlanarTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                            const Region region, std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixPerpenTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, std::string arrayName,
                            const Region region, std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixOffsetTrifurcation(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd,
                            std::string arrayName,
                            const Region region, const Region polyRegion,
                            std::vector<int> allEdges,
                            std::vector<int> badEdges, vtkIdList *critPts);
  int FixCloseGroup(vtkPolyData *pd, vtkPolyData *origPd, vtkPolyData *polyPd,
                    std::string arrayName,
                    const Region region, const Region polyRegion,
                    std::vector<int> allEdges,
                    std::vector<int> badEdges, vtkIdList *critPts);
  int FixGroupsWithPolycube();

  int GetCellRingNeighbors(vtkPolyData *pd,
                           vtkIdList *cellIds,
                           int ringNumber,
                           int totNumberOfRings,
                           std::vector<std::vector<int> > &neighbors);

  int GetConnectedEdges(std::vector<std::vector<int> > inputEdges,
                        std::vector<std::vector<int> > &connectedCornerPts);
  int FixPatchesWithPolycube();
  int FixPatchesWithPolycubeOld();
  int ParameterizeSurface(vtkPolyData *fullMapPd);
  int ParameterizeVolume(vtkPolyData *fullMapPd, vtkUnstructuredGrid *loftedVolume);
  int FormParametricHexMesh(vtkPolyData *polycubePd, vtkStructuredGrid *paraHexMesh,
                            int w_div, int &l_div, int h_div);
  int CheckFace(vtkPolyData *polycubePd, int faceId,
                int &nTopPts, int &nBotPts,
                int &flatTop, int &flatBot);
  int GetInteriorPointMaps(vtkPolyData *pdWithAllInterior,
                           vtkPolyData *pdWithCleanInterior,
                           vtkPolyData *pdWithoutInterior,
                           std::vector<int> &ptMap,
                           std::vector<std::vector<int> > &invPtMap);
  int GetVolumePointMap(vtkUnstructuredGrid *ugAll,
                        vtkUnstructuredGrid *ugClean,
                        std::vector<int> &ptMap);
  int MapInteriorBoundary(vtkStructuredGrid *paraHexVolume,
                          vtkPolyData *mappedSurface,
                          const std::vector<int> ptMap);
  int FixInteriorBoundary(vtkPolyData *mappedSurface,
                          const std::vector<std::vector<int> > invPtMap);
  int FixVolume(vtkUnstructuredGrid *mappedVolume,
                vtkUnstructuredGrid *cleanVolume,
                const std::vector<int> ptMap);
  int MapVolume(vtkStructuredGrid *paraHexVolume,
                vtkPolyData *mappedSurface,
                vtkStructuredGrid *mappedVolume);
  int ConvertUGToSG(vtkUnstructuredGrid *ug,
                    vtkStructuredGrid *sg,
                    const int w_div, const int l_div, const int h_div);
  int SmoothStructuredGrid(vtkStructuredGrid *hexMesh, const int iters);
  int SmoothUnstructuredGrid(vtkUnstructuredGrid *hexMesh, const int iters);
  int RemoveInteriorCells(vtkPolyData *quadMesh);
  int PushStructuredGridXAxis(vtkStructuredGrid *paraHexMesh,
                              const double midPt0[3],
                              const double midPt1[3],
                              const int isBottom);
  int PushStructuredGridZAxis(vtkStructuredGrid *paraHexMesh,
                              const double midPt0[3],
                              const double midPt1[3],
                              const int isBottom);

  char *CenterlineGroupIdsArrayName;
  char *CenterlineRadiusArrayName;
  char *GroupIdsArrayName;
  char *BlankingArrayName;

  vtkPolyData *WorkPd;
  vtkPolyData *GraphPd;
  vtkPolyData *Centerlines;
  vtkPolyData *MergedCenterlines;
  vtkPolyData *PolycubePd;
  vtkIdList *CenterlineGroupIds;

  vtkUnstructuredGrid *PolycubeUg;
  vtkUnstructuredGrid *FinalHexMesh;

  int ClipAllCenterlineGroupIds;
  int UseRadiusInformation;
  int UseVmtkClipping;
  int EnforceBoundaryDirections;
  int IsVasculature;
  int NumberOfCenterlineRemovePts;
  int ModifyCenterlines;
  int PolycubeDivisions;

  double CutoffRadiusFactor;
  double ClipValue;
  double PolycubeUnitLength;
  double CenterlineSeparationThreshold;
  double NormalsWeighting;

  vtkSVCenterlineGraph *CenterlineGraph;

private:
  vtkSVGroupsSegmenter(const vtkSVGroupsSegmenter&);  // Not implemented.
  void operator=(const vtkSVGroupsSegmenter&);  // Not implemented.
};

#endif
