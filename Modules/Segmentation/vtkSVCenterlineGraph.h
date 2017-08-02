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

/** @file vtkSVCenterlineGraph.h
 *  @brief a binary tree for creating a graph structure of 3D surfaces
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef vtkSVCenterlineGraph_h
#define vtkSVCenterlineGraph_h

#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSVCenterlineGCell.h"
#include "vtkSVParameterizationModule.h" // For exports

#include <map>
#include <list>

class VTKSVPARAMETERIZATION_EXPORT vtkSVCenterlineGraph : public vtkObject
{
public:
  static vtkSVCenterlineGraph* New();
  vtkTypeMacro(vtkSVCenterlineGraph, vtkObject);

  //Constructors
  vtkSVCenterlineGraph();
  vtkSVCenterlineGraph(int rootId,
          vtkPolyData *linesPd,
          std::string groupIdsArrayName);

  //Destructor
  ~vtkSVCenterlineGraph();

  //Member functions
  vtkSVCenterlineGCell* NewCell(int a_Id, vtkSVCenterlineGCell *a_Parent);
  vtkSVCenterlineGCell* NewCell(int a_Id, int a_BranchDir);
  vtkSVCenterlineGCell* NewCell(int a_Id, int a_BranchDir, double a_StartPt[3], double a_EndPt[3]);
  vtkSVCenterlineGCell* GetCell(const int findId);
  vtkSVCenterlineGCell* LookUp(vtkSVCenterlineGCell *lookCell, const int findId);
  int BuildGraph();
  int PrintGraph();
  int GetPolycube(const double height, const double width, vtkUnstructuredGrid *outUg);
  int GetCubeType(vtkSVCenterlineGCell *gCell, int &type);
  int TrifurcationDetermination(vtkSVCenterlineGCell *gCell, int &type);
  int AddBranchCube(vtkPoints *newPoints,
                    vtkCellArray *cellArray,
                    vtkPoints *points,
                    const int groupId,
                    vtkIntArray *localPtIds,
                    vtkIntArray *groupIds,
                    vtkIntArray *patchIds,
                    vtkDoubleArray *textureCoordinates,
                    const int type);
  int FormBifurcation(const double pt0[3], const double pt1[3],
                      const double pt2[3], const double pt3[3],
                      const double pt4[3], const double pt5[3],
                      const double centerPt[3],
                      const double factor,
                      double vecs[3][3],
                      double returnPts[2][3]);
  int GetBifurcationPoint(const double startPt[3], const double vec0[3], const double vec1[3], const double vec2[3], const double factor, double returnPt[3]);
  int GrowGraph(vtkSVCenterlineGCell *parent);
  int GetGraphDirections();
  int GetGraphPoints();
  int ComputeGlobalReferenceVectors(vtkSVCenterlineGCell *parent);
  int ComputeBranchReferenceVectors(vtkSVCenterlineGCell *parent);
  int GetInitialBranchDirections(vtkSVCenterlineGCell *parent);
  int UpdateBranchDirs(vtkSVCenterlineGCell *parent, const int updateDir);
  int GetGraphPolyData(vtkPolyData *pd);
  int GetConnectingLineGroups(const int groupId, vtkIdList *connectingGroups);

  //Static Member functions
  static int Recurse(vtkSVCenterlineGCell *rootvtkSVCenterlineGCell,
         int(*function)(vtkSVCenterlineGCell *currentvtkSVCenterlineGCell, void *arg0, void *arg1, void *arg2),
         void *rec_arg0, void *rec_arg1, void *rec_arg2);
  static int PrintGCell(vtkSVCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2);
  static int InsertGCellPoints(vtkSVCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2);


  int ComputeLocalCoordinateSystem(const double vz[3], const double vstart[3],
                                   double vx[3], double vy[3]);

  int RotateVecAroundLine(const double inVec[3],
                          const double angle,
                          const double axis[3],
                          double outVec[3]);

  int FlipLinePoints(vtkPolyData *pd, const int cellId);


  /// \brief directions of nodes in graph simplification
  enum DIRECTIONS
  {
    RIGHT = 0,
    BACK,
    LEFT,
    FRONT,
    UP,
    DOWN
  };

  //Member data
  vtkSVCenterlineGCell *Root;
  int NumberOfCells;
  int NumberOfNodes;

  //Member data needed to build
  vtkPolyData *Lines;
  std::string GroupIdsArrayName;
  double ReferenceVecs[3][3];

private:
  vtkSVCenterlineGraph(const vtkSVCenterlineGraph&); // Not implemented.
  void operator=(const vtkSVCenterlineGraph&); // Not implemented.
};

#endif
