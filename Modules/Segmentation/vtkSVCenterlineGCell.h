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
 *  \brief a node of Graphs
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#ifndef vtkSVCenterlineGCell_h
#define vtkSVCenterlineGCell_h

#include "vtkPolyData.h"
#include "vtkObject.h"
#include "vtkSVParameterizationModule.h" // For exports

class VTKSVPARAMETERIZATION_EXPORT vtkSVCenterlineGCell : public vtkObject
{
public:
  static vtkSVCenterlineGCell* New();
  vtkTypeMacro(vtkSVCenterlineGCell, vtkObject);

  //Constructors
  vtkSVCenterlineGCell();
  vtkSVCenterlineGCell(int a_Id,
                       int a_GroupId,
                       vtkSVCenterlineGCell *a_Parent);
  vtkSVCenterlineGCell(int a_Id, int a_GroupId, int a_BranchDir);
  vtkSVCenterlineGCell(int a_Id, int a_GroupId, int a_BranchDir, double a_StartPt[3], double a_EndPt[3]);

  //Destructor
  ~vtkSVCenterlineGCell();

  int GetTrifurcationType(int &type);
  int GetCubePoints(const double height, const double width,
                    vtkPoints *allPoints, vtkCellArray *allCells,
                    vtkIntArray *groupIds, vtkIntArray *patchIds);
  int GetBeginningType(int &beginningType, int &splitType);
  int GetEndType(int &endType, int &splitType);
  int GetSquare(const double startPt[3], const double vec0[3],
                const double vec1[3], const double height, const double width,
                vtkPoints *points);
  int GetWedge(const double pt0[3], const double pt1[3],
               const double pt2[3], const double vec0[3], const double height,
               vtkPoints *points);
  int FormBifurcation(const double pt0[3], const double pt1[3],
                      const double pt2[3], const double pt3[3],
                      const double pt4[3], const double pt5[3],
                      const double centerPt[3],
                      const double factor,
                      double vecs[3][3],
                      double returnPts[2][3]);
  int FormTrifurcation(const double pt0[3], const double pt1[3],
                       const double pt2[3], const double pt3[3],
                       const double pt4[3], const double pt5[3],
                       const double centerPt[3],
                       const double factor,
                       double vecs[3][3],
                       double returnPts[2][3]);
  int GetBifurcationPoint(const double startPt[3],
                          const double vec0[3],
                          const double vec1[3],
                          const double vec2[3],
                          const double factor,
                          double returnPt[3]);

  //Member data
  vtkSVCenterlineGCell *Parent;
  std::vector<vtkSVCenterlineGCell *> Children;
  int Id;
  int GroupId;
  int BranchDir;
  double RefAngle;
  double RefDirs[3][3];
  double StartPt[3];
  double EndPt[3];
  double BranchVec[3];
  int CornerPtIds[8];
  int DivergingChild;
  int AligningChild;
  int IsAlign;

  /// \brief directions of nodes in graph simplification
  enum SV_DIRECTIONS
  {
    RIGHT = 0,
    BACK,
    LEFT,
    FRONT,
  };

  /// \brief possible split types
  enum SV_CUBE_TYPE
  {
    NONE = 0,
    VERT_WEDGE,
    HORZ_WEDGE,
    C_TET_0, // Corner tets
    C_TET_1,
    C_TET_2,
    C_TET_3,
    S_TET_0, // Side tets
    S_TET_1,
    S_TET_2,
    S_TET_3,
    NOTHANDLED
  };

  /// \brief possible split types
  enum SV_SPLIT_TYPE
  {
    ZERO = 0,
    UNO,
    BI,
    TRI,
    QUAD,
    PENT,
    TOOMANY
  };

private:
  vtkSVCenterlineGCell(const vtkSVCenterlineGCell&); // Not implemented.
  void operator=(const vtkSVCenterlineGCell&); // Not implemented.
};

#endif
