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


/**
 * \class vtkSVPolyDataSliceAndDiceFilter
 *
 * \brief This is a filter to decompose a polydata surface into a series
 * of patches that represent a polycube structure using the object centerlines
 * \details The centerlines are first processed and a graph simplification
 * of the basic structure of the polydata is constructed. The polydata is
 * segmented at each bifurcation and a new patch is created to become the linking
 * cube at branching locations of a polycube structure. Finally, the polycube
 * is constructed based on the final decomposition and the necessary information
 * attached to the polycube in order to parameterize each individual patch
 *
 * \author Adam Updegrove
 * \author updega2@gmail.com
 * \author UC Berkeley
 * \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVPolyDataSliceAndDiceFilter_h
#define vtkSVPolyDataSliceAndDiceFilter_h

#include "vtkPolyDataAlgorithm.h"

#include "vtkEdgeTable.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkSVGeneralizedPolycube.h"

#include "svGraph.h"

#include <map>
#include <list>

class vtkSVPolyDataSliceAndDiceFilter : public vtkPolyDataAlgorithm
{
public:
  static vtkSVPolyDataSliceAndDiceFilter* New();
  void PrintSelf(ostream& os, vtkIndent indent);

  //@{
  /// \brief Set/Get macros for the slice length along the vessel length. Governs
  /// how often the polydata is sliced for the surgery lines
  vtkSetMacro(SliceLength, double);
  vtkGetMacro(SliceLength, double);
  //@}

  //@{
  /// \brief Set/Get macros to a indicate whether the polycube will be constructed
  vtkSetMacro(ConstructPolycube, int);
  vtkGetMacro(ConstructPolycube, int);
  //@}

  //@{
  /// \brief Set/Get macros for all the array names either used in the filter
  /// or assigned to a surface during the course of the operation
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
  //@}

  //@{
  /// \brief Set/Get macro for the centerlines. CenterlinesPd must be provided!
  vtkSetObjectMacro(CenterlinesPd, vtkPolyData);
  vtkGetObjectMacro(CenterlinesPd, vtkPolyData);
  //@}

  //@{
  /// \brief Set/Get macros for the final polycube structure
  vtkSetObjectMacro(Polycube, vtkSVGeneralizedPolycube);
  vtkGetObjectMacro(Polycube, vtkSVGeneralizedPolycube);
  //@}

  /// \brief Get macros for the surgery lines along length of object
  vtkGetObjectMacro(SurgeryLinesPd, vtkPolyData);

  /**
   * \brief Constructs cubes to populate a polycube structure from node of svGraph
   * \details This function is constructed to be used with svGraph::Recurse. The
   * void pointers are made available to pass any extra information through the function.
   * \param *gCell Node of graph to generate cubes for a polycube structure.
   * \return SV_OK if function completes without error
   */
  static int GraphToPolycube(svGCell *gCell, void *arg0,
                             void *arg1, void *arg2);
  /**
   * \brief Gets the index on the polycube based on the parent node and diverging
   * node directions.
   * \param parent direction of parent (Numbers correspond to DIRECTIONS enum)
   * \param divchild direction of parent (Numbers correspond to DIRECTIONS enum)
   * \param index index of cube to get (0-7)
   * \return SV_OK if function completes without error
   */
  static int LookupIndex(const int parent, const int divchild, const int index);


  /**
   * \brief directions of nodes in graph simplification
   */
  enum DIRECTIONS
  {
    RIGHT = 0,
    LEFT,
    FRONT,
    BACK,
    UP,
    DOWN
  };

  const static int DT[6][4];
  const static int RT[9][8];

protected:
  vtkSVPolyDataSliceAndDiceFilter();
  ~vtkSVPolyDataSliceAndDiceFilter();

  // Usual data generation method
  int RequestData(vtkInformation *vtkNotUsed(request),
		  vtkInformationVector **inputVector,
		  vtkInformationVector *outputVector);

  int PrepFilter(); // Prep work
  int RunFilter(); // Run filter operations
  int FindGroupBoundaries(); ///< \brief Find and label nodes that form boundary between groups
  int BuildPolycube(); ///< \brief If building polycube, construct it based off graph
  int GetCriticalPoints(); ///< \brief Get points that separate three or more regions
  int SliceBifurcations(); ///< \brief Process and cut each bifurcation
  int SliceBranches(); ///< \brief Process and cut each branch
  void CheckLength(int &ptId,
                  const int numPts,
                  int &done);
  int CheckSlice(vtkPolyData *pd);
  void UpdatePtId(int &ptId);
  int GetCorrectFrontPoint(vtkPolyData *pd,
                           double frontDir[3],
                           int &frontId,
                           int &backId);
  int GetFourPolyDataRegions(vtkPolyData *startPd,
                             const int id0, vtkPolyData *pd0,
                             const int id1, vtkPolyData *pd1,
                             const int id2, vtkPolyData *pd2,
                             vtkPolyData *leftovers);
  int CriticalSurgeryPoints(vtkPolyData *pd,
                           const int frontId,
                           const int backId,
                           const int groupId,
                           const int checkId,
                           double startPt[3],
                           double secondPt[3],
                           vtkIdList *fixedGoToPoints,
                           vtkIdList *fixedSurgeryPoints,
                           double startDir[3]);
  int GetFirstSurgeryPoints(vtkPolyData *pd, int pointId,
                            vtkIdList *surgeryPoints,
                            double xvec[3], double zvec[3]);
  int GetSurgeryPoints(vtkPolyData *pd,
                       vtkPolyData *parentPd,
                       vtkDataArray *pointIds,
                       const double clStartPt[3],
                       const double clSecondPt[3],
                       const int front,
                       const int back,
                       const int checkId,
                       std::string arrayName,
                       vtkIdList *surgeryPoints,
                       double startDir[3]);
  int GetHalfSurgeryPoints(vtkPolyData *pd,
                           vtkDataArray *pointIds,
                           const int cellId,
                           const int front,
                           const int back,
                           vtkIdList *surgeryPoints);
  int GetNextSurgeryPoints(vtkPolyData *pd,
                           double centerPt[3],
                           vtkIdList *surgeryPoints,
                           double xvec[3], double zvec[3],
                           double radius,
                           vtkIdList *surgeryLineIds);
  int GetEndSurgeryPoints(vtkPolyData *pd, svGCell *gCell,
                          double centerPt[3],
                          vtkIdList *surgeryPoints,
                          int endSurgeryIds[8],
                          double xvec[3], double zvec[3],
                          double radius,
                          vtkIdList *surgeryLineIds,
                          int cellIndices[8],
                          int &secondRun);
  int DetermineSliceStrategy(vtkPolyData *branchPd,
                             svGCell *gCell,
                             vtkPolyData *branchCenterline,
                             int &branchStartPtId,
                             vtkIdList *surgeryPoints,
                             int &centerlineStartPtId,
                             int &strategy);
  int InsertCriticalPoints(const int pointId, vtkIdList *groupIds);
  int GetBranch(const int branchId, vtkPolyData *branchPd,
                vtkPolyData *branchCenterlinesPd);
  int SliceBranch(vtkPolyData *branchPd, vtkPolyData *branchCenterline,
                    svGCell *gCell,
                    vtkDataArray *sliceIds,
                    vtkPoints *surgeryPts,
                    vtkCellArray *surgeryLines,
                    vtkIntArray *surgeryData,
                    int secondRun);
  int SliceBifurcation(vtkPolyData *pd,
                       svGCell *gCell);
  int FixGraphDirections(svGCell *gCell, const int actualId,
                         int cellIndices[8]);
  int GetSectionZAxis(const double endPt[3], const double startPt[3],
                      double zvec[3]);
  int GetSectionXAxis(const double endPt[3], const double startPt[3],
                      const double surfacePt[3], double xvec[3]);

  int GetClose3DPoint(vtkPolyData *pd, double centerPt[3],
                      const int startPtId, int &returnStartId,
                      double xvec[3], double zvec[3],
                      double radius,
                      vtkPolyData *boundary);
  int GetContourSecondPoint(vtkPolyData *pd, int ptId, double centerPt[3],
                            double zvec[3], int &startSecondId);
  int AddSurgeryPoints(vtkIdList *surgeryLineIds, vtkPoints *surgeryPts,
                       vtkCellArray *surgeryLines, vtkIntArray *surgeryData);


private:
  vtkSVPolyDataSliceAndDiceFilter(const vtkSVPolyDataSliceAndDiceFilter&);  // Not implemented.
  void operator=(const vtkSVPolyDataSliceAndDiceFilter&);  // Not implemented.

  char *BoundaryPointsArrayName; // Array for boundary points (internal)
  char *DijkstraArrayName; // Array for dijkstra filter (internal)
  char *GroupIdsArrayName; // Array for groupids; must be defined on centerlines
  char *InternalIdsArrayName; // Array to keep track of node numbers wehn thresholding and slicing (internal)
  char *SegmentIdsArrayName; // Array for new patch ids
  char *SliceIdsArrayName; // Unnecessary array name TODO: Remove
  char *SphereRadiusArrayName;  // Array for maximum inscribed sphere radius; must be defined on centerlines

  vtkPolyData    *InitialPd; // Input polydata, kept for reference
  vtkPolyData    *WorkPd; // Polydata used during filter processing
  vtkPolyData    *GraphPd; // Polydata populated with the graph simplification
  vtkPolyData    *CenterlinesPd; // Polydata containing centerlines; must be provided!
  vtkPolyData    *SurgeryLinesPd; // Polydata containing the lines slicing the vessels along their length
  svGraph        *CenterlineGraph; // A graph simplification
  vtkSVGeneralizedPolycube *Polycube; // The polycube structure

  std::multimap<int , int> CriticalPointMap; // Stl map keeping track of what points map to what groups and vice versa

  int ConstructPolycube;
  int SliceDirection;
  int TotalSliceId;
  int IT[6][6][8];

  double FirstBranchVec[3];
  double SliceLength;

};

#endif
