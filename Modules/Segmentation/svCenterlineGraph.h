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

/** @file svCenterlineGraph.h
 *  @brief a binary tree for creating a graph structure of 3D surfaces
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#ifndef svCenterlineGraph_h
#define svCenterlineGraph_h

#include "vtkPolyData.h"
#include "svCenterlineGCell.h"
#include "vtkSVParameterizationModule.h" // For exports

#include <map>
#include <list>

class VTKSVPARAMETERIZATION_EXPORT svCenterlineGraph
{
public:
  //Constructors
  svCenterlineGraph();
  svCenterlineGraph(int rootId,
          vtkPolyData *linesPd,
          std::string groupIdsArrayName);

  //Destructor
  ~svCenterlineGraph();

  //Member functions
  svCenterlineGCell* NewCell(int a_Id, svCenterlineGCell *a_Parent);
  svCenterlineGCell* NewCell(int a_Id, int a_Dir, double a_StartPt[3], double a_EndPt[3]);
  svCenterlineGCell* GetCell(const int findId);
  svCenterlineGCell* LookUp(svCenterlineGCell *lookCell, const int findId);
  int BuildGraph();
  int PrintGraph();
  int GrowGraph(svCenterlineGCell *parent);
  int RefineGraphWithLocalCoordinates();
  int ComputeReferenceVectors(svCenterlineGCell *parent);
  int GetNewBranchDirections(svCenterlineGCell *parent);
  int GetGraphPolyData(vtkPolyData *pd);
  int GetConnectingLineGroups(const int groupId, vtkIdList *connectingGroups);

  //Static Member functions
  static int GetDirectionVector(const int dir, double dirVector[3]);
  static int Recurse(svCenterlineGCell *rootsvCenterlineGCell,
         int(*function)(svCenterlineGCell *currentsvCenterlineGCell, void *arg0, void *arg1, void *arg2),
         void *rec_arg0, void *rec_arg1, void *rec_arg2);
  static int UpdateCellDirection(svCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2);
  static int PrintGCell(svCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2);
  static int InsertGCellPoints(svCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2);

   /*
   * \brief Gets the index on the polycube based on the parent node and diverging
   * node directions.
   * \param parent direction of parent (Numbers correspond to DIRECTIONS enum)
   * \param divchild direction of parent (Numbers correspond to DIRECTIONS enum)
   * \param index index of cube to get (0-7)
   * \return SV_OK if function completes without error
   */
  static int LookupIndex(const int parent, const int divchild, const int index);

  int ComputeLocalCoordinateSystem(const double vz[3], const double vstart[3],
                                   double vx[3], double vy[3]);

  int FlipLinePoints(vtkPolyData *pd, const int cellId);


  /// \brief directions of nodes in graph simplification
  enum DIRECTIONS
  {
    RIGHT = 0,
    LEFT,
    FRONT,
    BACK,
    UP,
    DOWN
  };

  const static int DT[6][4]; ///< \brief Direction Table
  const static int RT[9][8]; ///< \brief Index Rotation Table

  //Member data
  svCenterlineGCell *Root;
  int NumberOfCells;
  int NumberOfNodes;

  //Member data needed to build
  vtkPolyData *Lines;
  std::string GroupIdsArrayName;
  double ReferenceVecs[3][3];
};

#endif
