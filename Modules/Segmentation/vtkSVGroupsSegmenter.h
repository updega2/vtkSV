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
#include "vtkIdList.h"
#include "svCenterlineGraph.h"

#include "vtkSVSegmentationModule.h" // For export

struct Region
{
  int Index;
  int IndexCluster;

  int NumberOfCorners;
  std::vector<int> CornerPoints;

  std::vector<std::vector<int> > BoundaryEdges;

  int numInteriorVerts;
  std::vector<int> InteriorVerts;

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

  /** \brief Correct cells on the boundary by updating val if they have
   *  multiple neighboring cells of the same value */
  static int CorrectCellBoundaries(vtkPolyData *pd, std::string cellArrayName);

  /** \brief Pass an array from points to cells by using the most
   * occuring point value. */
  static int PassPointGroupsToCells(vtkPolyData *pd, std::string pointArrayName);

  /** \brief Naive implementation to get most reoccuring number in list. Okay
   *  because list size is small. */
  static void GetMostOccuringVal(vtkIdList *idList, int &output, int &max_count);

  /** \brief Run basic edge weithing cvt with pd */
  static int RunEdgeWeightedCVT(vtkPolyData *pd);

  /** \brief From three vectors, compute transformation from global to local */
  static int ComputeRotationMatrix(const double vx[3], const double vy[3],
                                   const double vz[3], double rotMatrix[9]);

  static int SmoothBoundaries(vtkPolyData *pd, std::string arrayName);

  static int GetRegions(vtkPolyData *pd, std::string arrayName,
                        std::vector<Region> &allRegions);

  static int CurveFitBoundaries(vtkPolyData *pd, std::string arrayName,
                                std::vector<Region> allRegions);

  static int GetCCWPoint(vtkPolyData *pd, const int pointId, const int cellId);

  static int CheckBoundaryEdge(vtkPolyData *pd, std::string arrayName, const int cellId, const int pointId0, const int pointId1);

  static void SplineKnots(std::vector<int> &u, int n, int t);

  static void SplineCurve(const std::vector<XYZ> &inp, int n, const std::vector<int> &knots, int t, std::vector<XYZ> &outp, int res);

  static void SplinePoint(const std::vector<int> &u, int n, int t, double v, const std::vector<XYZ> &control, XYZ &output);

  static double SplineBlend(int k, int t, const std::vector<int> &u, double v);

  const static double GlobalCoords[3][3];

  std::vector<Region> Polycube;

protected:
  vtkSVGroupsSegmenter();
  ~vtkSVGroupsSegmenter();

  // Usual data generation method
  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);

  int PrepFilter(); // Prep work.
  int RunFilter(); // Run filter operations.

  int GetPatches();

  char *CenterlineGroupIdsArrayName;
  char *CenterlineRadiusArrayName;
  char *GroupIdsArrayName;
  char *BlankingArrayName;

  vtkPolyData *WorkPd;
  vtkPolyData *Centerlines;
  vtkPolyData *CenterlinesWorkPd;
  vtkIdList *CenterlineGroupIds;

  int ClipAllCenterlineGroupIds;
  int UseRadiusInformation;

  double CutoffRadiusFactor;
  double ClipValue;

  svCenterlineGraph *CenterlineGraph;

private:
  vtkSVGroupsSegmenter(const vtkSVGroupsSegmenter&);  // Not implemented.
  void operator=(const vtkSVGroupsSegmenter&);  // Not implemented.
};

#endif