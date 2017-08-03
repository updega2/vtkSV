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

/** @file vtkSVCenterlineGCell.cxx
 *  @brief a node of Graphs
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVCenterlineGCell.h"

#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVCenterlineGCell);

// ----------------------
// Constructor
// ----------------------
vtkSVCenterlineGCell::vtkSVCenterlineGCell()
{
  this->Parent   = NULL;
  this->Id       = -1;
  this->GroupId  = -1;
  this->BranchDir      = -1;
  this->DivergingChild = -1;
  this->AligningChild  = -1;
  this->IsAlign        = -1;
  this->RefAngle       = 0.0;
  for (int i=0; i<3; i++)
  {
    this->StartPt[i]  = -1.0;
    this->EndPt[i]    = -1.0;
    this->BranchVec[i] = -1.0;
    for (int j=0; j<3; j++)
      this->RefDirs[i][j] = -1.0;
  }
  for (int i=0; i<8; i++)
    this->CornerPtIds[i] = -1;
}

// ----------------------
// Constructor
// ----------------------
vtkSVCenterlineGCell::vtkSVCenterlineGCell(int a_Id, int a_GroupId,
                                           vtkSVCenterlineGCell *a_Parent)
{
  this->Parent   = a_Parent;
  this->Id       = a_Id;
  this->GroupId  = a_GroupId;
  this->BranchDir      = -1;
  this->DivergingChild = -1;
  this->AligningChild  = -1;
  this->IsAlign        = -1;
  this->RefAngle       = 0.0;
  for (int i=0; i<3; i++)
  {
    if (a_Parent != NULL)
      this->StartPt[i] = a_Parent->EndPt[i];
    else
      this->StartPt[i] = -1.0;
    this->EndPt[i]    = -1.0;
    this->BranchVec[i] = -1.0;
    for (int j=0; j<3; j++)
      this->RefDirs[i][j] = -1.0;
  }
  for (int i=0; i<8; i++)
    this->CornerPtIds[i] = -1;
}

// ----------------------
// Constructor
// ----------------------
vtkSVCenterlineGCell::vtkSVCenterlineGCell(int a_Id, int a_GroupId,
                                           int a_BranchDir, double a_StartPt[3],
                                           double a_EndPt[3])
{
  this->Parent   = NULL;
  this->Id       = a_Id;
  this->GroupId  = a_GroupId;
  this->BranchDir      = a_BranchDir;
  this->DivergingChild = -1;
  this->AligningChild  = -1;
  this->IsAlign        = -1;
  this->RefAngle       = 0.0;
  for (int i=0; i<3; i++)
  {
    this->StartPt[i] = a_StartPt[i];
    this->EndPt[i]   = a_EndPt[i];
    this->BranchVec[i] = -1.0;
    for (int j=0; j<3; j++)
      this->RefDirs[i][j] = -1.0;
  }
  for (int i=0; i<8; i++)
    this->CornerPtIds[i] = -1;
}

// ----------------------
// Constructor
// ----------------------
vtkSVCenterlineGCell::vtkSVCenterlineGCell(int a_Id, int a_GroupId,
                                           int a_BranchDir)
{
  this->Parent   = NULL;
  this->Id       = a_Id;
  this->GroupId  = a_GroupId;
  this->BranchDir      = a_BranchDir;
  this->DivergingChild = -1;
  this->AligningChild  = -1;
  this->IsAlign        = -1;
  this->RefAngle       = 0.0;
  for (int i=0; i<3; i++)
  {
    this->StartPt[i]  = -1.0;
    this->EndPt[i]    = -1.0;
    this->BranchVec[i] = -1.0;
    for (int j=0; j<3; j++)
      this->RefDirs[i][j] = -1.0;
  }
  for (int i=0; i<8; i++)
    this->CornerPtIds[i] = -1;
}

// ----------------------
// Destructor
// ----------------------
vtkSVCenterlineGCell::~vtkSVCenterlineGCell()
{
  for (int i=0; i<this->Children.size(); i++)
  {
    if (this->Children[i] != NULL)
    {
      this->Children[i]->Delete();
      this->Children[i] = NULL;
    }
  }
}

//// ----------------------
//// GetBeginningType
//// ----------------------
//int vtkSVCenterlineGCell::GetBeginningType(int &beginningType)
//{
//  if (this->Parent == NULL)
//    beginningType = NONE;
//  else if (this->Parent->Children.size() == 2)
//    beginningType = BI;
//  else if (this->Parent->Children.size() == 3)
//    beginningType = TRI;
//  else
//  {
//    endType = NOTHANDLED;
//    fprintf(stdout,"Cannot currently handle more than a trifurcation\n");
//    return SV_ERROR;
//  }
//
//  return SV_OK;
//}
//
//// ----------------------
//// GetEndType
//// ----------------------
//int vtkSVCenterlineGCell::GetEndType(int &endType)
//{
//  if (this->Children.size() == 0)
//    endType == NONE;
//  else if (this->Children.size() == 2)
//    endType = BI;
//  else if (this->Children.size() == 3)
//    endType = TRI;
//  else
//  {
//    endType = NOTHANDLED;
//    fprintf(stdout,"Cannot currently handle more than a trifurcation\n");
//    return SV_ERROR;
//  }
//
//  return SV_OK;
//}

// ----------------------
// GetCubeType
// ----------------------
int vtkSVCenterlineGCell::GetCubeType(int &type)
{
  if (this->Parent == NULL)
  {
    if (this->Children.size() == 0)
      type = 0;
    else if (this->Children.size() == 2)
      type = 1;
    else
    {
      this->GetTrifurcationType(type);
      if (type == 0)
      {
        fprintf(stderr,"Cannot currently Polycube for this type of trifurcation yet!!!\n");
        return SV_ERROR;
      }
    }
  }
  else if (this->Children.size() == 2)
  {
    if ((this->Parent->BranchDir + this->BranchDir)%2 == 0)
    {
      if ((this->BranchDir + this->Children[this->DivergingChild]->BranchDir)%2 == 0)
        type = 2;
      else
        type = 4;
    }
    else
    {
      if ((this->BranchDir + this->Children[this->DivergingChild]->BranchDir)%2 == 0)
        type = 3;
      else
        type = 5;
    }
  }
  else if (this->Children.size() == 0)
  {
    if ((this->BranchDir)%2 == 0)
      type = 6;
    else
      type = 7;
  }
  else
  {
    this->GetTrifurcationType(type);
    if (type == 0)
    {
      fprintf(stderr,"Cannot currently Polycube for this type of trifurcation yet!!!\n");
      return SV_ERROR;
    }
  }

  return SV_OK;
}

// ----------------------
// GetTrifurcationType
// ----------------------
int vtkSVCenterlineGCell::GetTrifurcationType(int &type)
{
  int numChildren = this->Children.size();

  if (numChildren != 3)
  {
    fprintf(stderr,"Trifurcation wrongly identified\n");
  }

  // Count branch directions
  int rightCount = 0;
  int leftCount = 0;
  int backCount = 0;
  int frontCount = 0;
  int downCount = 0;

  // Get initial branch dirs
  for (int i=0; i<numChildren; i++)
  {
    int maxDir=0;
    double maxDot = -1.0;

    for (int j=0; j<3; j++)
    {
      double compare = fabs(vtkMath::Dot(this->Children[i]->RefDirs[j], this->Children[i]->BranchVec));
      //fprintf(stdout,"Would like to see compare, dir: %d, dot: %.6f\n", j, compare);
      if (compare > maxDot)
      {
        maxDot = compare;
        maxDir = j;
      }
    }
    double dotCheck = vtkMath::Dot(this->Children[i]->RefDirs[maxDir], this->Children[i]->BranchVec);
    if (maxDir == 0)
    {
      if (dotCheck <= 0)
        downCount++;
      else
      {
        fprintf(stderr,"WARN ABOUT THIS, WILL NEED TO DO SOMETHING SPECIAL\n");
        return SV_ERROR;
      }
    }
    else if (maxDir == 1)
    {
      if (dotCheck > 0)
        rightCount++;
      else
        leftCount++;
    }
    else if (maxDir == 2)
    {
      if (dotCheck > 0)
        backCount++;
      else
        fprintf(stdout,"FRONT\n");
    }
  }

  if (rightCount == 1 && leftCount == 1 && downCount == 1)
    type = 8;
  else
  {
    fprintf(stderr, "THIS ISNT HANDLED YET!\n");
    return SV_ERROR;
  }

  return SV_OK;
}

// ----------------------
// GetCubePoints
// ----------------------
int vtkSVCenterlineGCell::GetCubePoints(const double height,
                                        const double width,
                                        vtkPoints *allPoints,
                                        vtkCellArray *allCells,
                                        vtkIntArray *localPtIds,
                                        vtkIntArray *groupIds,
                                        vtkIntArray *patchIds)
{
  vtkNew(vtkPoints, points);

  int cubeType;
  this->GetCubeType(cubeType);
  fprintf(stdout,"CUBE GROUP ID: %d\n", this->GroupId);
  fprintf(stdout,"CUBE TYPE:     %d\n", cubeType);

  int groupId = this->GroupId;

  if (cubeType == 0)
  {
    int numPoints = 8;

    double finalPts[8][3];

    // Get top square from start pt
    double frontVec[3];
    vtkMath::Cross(this->RefDirs[1], this->RefDirs[0], frontVec);
    vtkMath::Normalize(frontVec);

    vtkNew(vtkPoints, startPoints);
    this->GetSquare(this->StartPt, this->RefDirs[1], frontVec,
                    height, width, startPoints);

    startPoints->GetPoint(0, finalPts[5]);
    startPoints->GetPoint(1, finalPts[1]);
    startPoints->GetPoint(2, finalPts[2]);
    startPoints->GetPoint(3, finalPts[6]);

    // Get square from end pt
    vtkNew(vtkPoints, endPoints);
    this->GetSquare(this->EndPt, this->RefDirs[1], frontVec,
                    height, width, endPoints);

    endPoints->GetPoint(0, finalPts[4]);
    endPoints->GetPoint(1, finalPts[0]);
    endPoints->GetPoint(2, finalPts[3]);
    endPoints->GetPoint(3, finalPts[7]);

    int pI[8];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
      localPtIds->InsertNextTuple1(i);
    }

    vtkIdType face0[4] = {pI[4], pI[5], pI[1], pI[0]};
    vtkIdType face1[4] = {pI[7], pI[6], pI[5], pI[4]};
    vtkIdType face2[4] = {pI[3], pI[2], pI[6], pI[7]};
    vtkIdType face3[4] = {pI[0], pI[1], pI[2], pI[3]};
    vtkIdType face4[4] = {pI[2], pI[1], pI[5], pI[6]};
    vtkIdType face5[4] = {pI[4], pI[0], pI[3], pI[7]};

    allCells->InsertNextCell(4, face0);
    allCells->InsertNextCell(4, face1);
    allCells->InsertNextCell(4, face2);
    allCells->InsertNextCell(4, face3);
    allCells->InsertNextCell(4, face4);
    allCells->InsertNextCell(4, face5);

    for (int i=0; i<6; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (cubeType == 1)
  {
    int numPoints = 10;

    double finalPts[10][3];

    double vecs[3][3];
    double triPts[2][3];
    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *diver = this->Children[this->DivergingChild];

    // Get points in between graph node and brother nodes
    this->FormBifurcation(this->StartPt, this->EndPt,
                          diver->EndPt, diver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, triPts);

    // get vector towards top of cube
    double frontVec[3];
    vtkMath::Cross(vecs[1], vecs[0], frontVec);
    vtkMath::Normalize(frontVec);

    // Get end wedge
    vtkNew(vtkPoints, endPoints);
    this->GetWedge(triPts[1], this->EndPt, triPts[0],
                   frontVec, height, endPoints);

    endPoints->GetPoint(0, finalPts[3]);
    endPoints->GetPoint(1, finalPts[4]);
    endPoints->GetPoint(2, finalPts[0]);
    endPoints->GetPoint(3, finalPts[8]);
    endPoints->GetPoint(4, finalPts[9]);
    endPoints->GetPoint(5, finalPts[5]);

    // Get square from start pt
    vtkNew(vtkPoints, startPoints);
    this->GetSquare(this->StartPt, this->RefDirs[1], frontVec,
                    height, width, startPoints);

    startPoints->GetPoint(0, finalPts[6]);
    startPoints->GetPoint(1, finalPts[1]);
    startPoints->GetPoint(2, finalPts[2]);
    startPoints->GetPoint(3, finalPts[7]);


    int pI[10];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
      localPtIds->InsertNextTuple1(i);
    }

    vtkIdType face0[4] = {pI[5], pI[6], pI[1], pI[0]};
    vtkIdType face1[5] = {pI[8], pI[7], pI[6], pI[5], pI[9]};
    vtkIdType face2[4] = {pI[3], pI[2], pI[7], pI[8]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};
    vtkIdType face4[4] = {pI[2], pI[1], pI[6], pI[7]};

    allCells->InsertNextCell(4, face0);
    allCells->InsertNextCell(5, face1);
    allCells->InsertNextCell(4, face2);
    allCells->InsertNextCell(5, face3);
    allCells->InsertNextCell(4, face4);

    for (int i=0; i<5; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (cubeType == 2)
  {
    int numPoints = 12;

    double finalPts[12][3];

    vtkSVCenterlineGCell *brother, *diver, *parent;
    parent = this->Parent;
    if (parent->Children[0]->Id == this->Id)
      brother = parent->Children[1];
    else
      brother = parent->Children[0];
    diver = parent->Children[parent->DivergingChild];

    // Get ending bifurcation points
    double vecs[3][3];
    double triPts[2][3];
    this->FormBifurcation(this->EndPt, this->StartPt,
                         parent->StartPt, parent->EndPt,
                         brother->EndPt, brother->StartPt,
                         this->StartPt,
                         width/2., vecs, triPts);

    // get diverging vector towards top of cube
    double diverVec[3];
    vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
    vtkMath::Normalize(diverVec);

    // get front vector
    double frontVec[3];
    vtkMath::Cross(diverVec, vecs[1], frontVec);
    vtkMath::Normalize(frontVec);

    // Get top wedge
    vtkNew(vtkPoints, startPoints);

    // Flip if diver branch is left or front
    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(frontVec, -1.0);

    // Treat differently for different dirs
    if (this->BranchDir == LEFT || this->BranchDir == FRONT)
      this->GetWedge(triPts[1], this->StartPt, triPts[0], frontVec,
                     height, startPoints);
    else
      this->GetWedge(triPts[0], this->StartPt, triPts[1], frontVec,
                     height, startPoints);

    startPoints->GetPoint(0, finalPts[1]);
    startPoints->GetPoint(1, finalPts[2]);
    startPoints->GetPoint(2, finalPts[3]);
    startPoints->GetPoint(3, finalPts[7]);
    startPoints->GetPoint(4, finalPts[8]);
    startPoints->GetPoint(5, finalPts[9]);

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, triPts);

    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], frontVec);
    vtkMath::Normalize(frontVec);

    // Get end wedge
    vtkNew(vtkPoints, endPoints);

    // Treat differently for different dirs
    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
    {
      vtkMath::MultiplyScalar(frontVec, -1.0);
      this->GetWedge(triPts[0], this->EndPt, triPts[1], frontVec,
                     height, endPoints);
    }
    else
      this->GetWedge(triPts[1], this->EndPt, triPts[0], frontVec,
                     height, endPoints);

    endPoints->GetPoint(0, finalPts[4]);
    endPoints->GetPoint(1, finalPts[5]);
    endPoints->GetPoint(2, finalPts[0]);
    endPoints->GetPoint(3, finalPts[10]);
    endPoints->GetPoint(4, finalPts[11]);
    endPoints->GetPoint(5, finalPts[6]);

    int pI[12];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
      localPtIds->InsertNextTuple1(i);
    }
    vtkIdType face0[4] = {pI[6], pI[7], pI[1], pI[0]};
    vtkIdType face1[6] = {pI[10], pI[9], pI[8], pI[7], pI[6], pI[11]};
    vtkIdType face2[4] = {pI[4], pI[3], pI[9], pI[10]};
    vtkIdType face3[6] = {pI[0], pI[1], pI[2], pI[3], pI[4], pI[5]};

    allCells->InsertNextCell(4, face0);
    allCells->InsertNextCell(6, face1);
    allCells->InsertNextCell(4, face2);
    allCells->InsertNextCell(6, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (cubeType == 3)
  {
    int numPoints = 12;

    double finalPts[12][3];

    vtkSVCenterlineGCell *brother, *diver, *parent;
    parent = this->Parent;
    if (parent->Children[0]->Id == this->Id)
      brother = parent->Children[1];
    else
      brother = parent->Children[0];
    diver = parent->Children[parent->DivergingChild];

    // Get ending bifurcation points
    double vecs[3][3];
    double triPts[2][3];
    this->FormBifurcation(this->EndPt, this->StartPt,
                         parent->StartPt, parent->EndPt,
                         brother->EndPt, brother->StartPt,
                         this->StartPt,
                         width/2., vecs, triPts);

    // get vector towards top of cube
    double diverVec[3];
    vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
    vtkMath::Normalize(diverVec);

    // get front vector towards top of cube
    double frontVec[3];
    vtkMath::Cross(diverVec, vecs[1], frontVec);
    vtkMath::Normalize(frontVec);

    // flip if divergent is left or front
    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(frontVec, -1.0);

    // Get start wedge
    vtkNew(vtkPoints, startPoints);

    // Treat differently for different dirs
    if (this->BranchDir == LEFT || this->BranchDir == FRONT)
      this->GetWedge(triPts[1], this->StartPt, triPts[0], frontVec,
                     height, startPoints);
    else
      this->GetWedge(triPts[0], this->StartPt, triPts[1], frontVec,
                     height, startPoints);

    startPoints->GetPoint(0, finalPts[9]);
    startPoints->GetPoint(1, finalPts[5]);
    startPoints->GetPoint(2, finalPts[1]);
    startPoints->GetPoint(3, finalPts[10]);
    startPoints->GetPoint(4, finalPts[6]);
    startPoints->GetPoint(5, finalPts[2]);

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, triPts);


    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], frontVec);
    vtkMath::Normalize(frontVec);

    // Get end wedge
    vtkNew(vtkPoints, endPoints);

    // treat differently for different dirs
    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
    {
      vtkMath::MultiplyScalar(frontVec, -1.0);
      this->GetWedge(triPts[0], this->EndPt, triPts[1], frontVec,
                     height, endPoints);
    }
    else
      this->GetWedge(triPts[1], this->EndPt, triPts[0], frontVec,
                     height, endPoints);

    endPoints->GetPoint(0, finalPts[0]);
    endPoints->GetPoint(1, finalPts[4]);
    endPoints->GetPoint(2, finalPts[8]);
    endPoints->GetPoint(3, finalPts[3]);
    endPoints->GetPoint(4, finalPts[7]);
    endPoints->GetPoint(5, finalPts[11]);

    int pI[12];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
      localPtIds->InsertNextTuple1(i);
    }

    vtkIdType face0[6] = {pI[8], pI[9], pI[5], pI[1], pI[0], pI[4]};
    vtkIdType face1[4] = {pI[11], pI[10], pI[9], pI[8]};
    vtkIdType face2[6] = {pI[3], pI[2], pI[6], pI[10], pI[11], pI[7]};
    vtkIdType face3[4] = {pI[0], pI[1], pI[2], pI[3]};

    allCells->InsertNextCell(6, face0);
    allCells->InsertNextCell(4, face1);
    allCells->InsertNextCell(6, face2);
    allCells->InsertNextCell(4, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (cubeType == 4)
  {
    int numPoints = 12;

    double finalPts[12][3];

    vtkSVCenterlineGCell *brother, *diver, *parent;
    parent = this->Parent;
    if (parent->Children[0]->Id == this->Id)
      brother = parent->Children[1];
    else
      brother = parent->Children[0];
    diver = parent->Children[parent->DivergingChild];

    // Get ending bifurcation points
    double vecs[3][3];
    double triPts[2][3];
    this->FormBifurcation(this->EndPt, this->StartPt,
                         parent->StartPt, parent->EndPt,
                         brother->EndPt, brother->StartPt,
                         this->StartPt,
                         width/2., vecs, triPts);

    // get diverging vector towards top of cube
    double diverVec[3];
    vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
    vtkMath::Normalize(diverVec);

    // get front vector
    double frontVec[3];
    vtkMath::Cross(diverVec, vecs[1], frontVec);
    vtkMath::Normalize(frontVec);

    // Get top wedge
    vtkNew(vtkPoints, startPoints);

    // Flip if diver branch is left or front
    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(frontVec, -1.0);

    // Treat differently for different dirs
    if (this->BranchDir == LEFT || this->BranchDir == FRONT)
      this->GetWedge(triPts[1], this->StartPt, triPts[0], frontVec,
                     height, startPoints);
    else
      this->GetWedge(triPts[0], this->StartPt, triPts[1], frontVec,
                     height, startPoints);

    startPoints->GetPoint(0, finalPts[1]);
    startPoints->GetPoint(1, finalPts[2]);
    startPoints->GetPoint(2, finalPts[3]);
    startPoints->GetPoint(3, finalPts[8]);
    startPoints->GetPoint(4, finalPts[9]);
    startPoints->GetPoint(5, finalPts[10]);

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, triPts);


    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], frontVec);
    vtkMath::Normalize(frontVec);

    // Get end wedge
    vtkNew(vtkPoints, endPoints);

    // treat differently for different dirs
    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
    {
      vtkMath::MultiplyScalar(frontVec, -1.0);
      this->GetWedge(triPts[0], this->EndPt, triPts[1], frontVec,
                     height, endPoints);
    }
    else
      this->GetWedge(triPts[1], this->EndPt, triPts[0], frontVec,
                     height, endPoints);

    endPoints->GetPoint(0, finalPts[0]);
    endPoints->GetPoint(1, finalPts[5]);
    endPoints->GetPoint(2, finalPts[7]);
    endPoints->GetPoint(3, finalPts[4]);
    endPoints->GetPoint(4, finalPts[6]);
    endPoints->GetPoint(5, finalPts[11]);

    int pI[12];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
      localPtIds->InsertNextTuple1(i);
    }
    vtkIdType face0[5] = {pI[7], pI[8], pI[1], pI[0], pI[5]};
    vtkIdType face1[5] = {pI[11], pI[10], pI[9], pI[8], pI[7]};
    vtkIdType face2[5] = {pI[4], pI[3], pI[10], pI[11], pI[6]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};

    allCells->InsertNextCell(5, face0);
    allCells->InsertNextCell(5, face1);
    allCells->InsertNextCell(5, face2);
    allCells->InsertNextCell(5, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (cubeType == 5)
  {
    int numPoints = 12;

    double finalPts[12][3];

    vtkSVCenterlineGCell *brother, *diver, *parent;
    parent = this->Parent;
    if (parent->Children[0]->Id == this->Id)
      brother = parent->Children[1];
    else
      brother = parent->Children[0];
    diver = parent->Children[parent->DivergingChild];

    // Get ending bifurcation points
    double vecs[3][3];
    double triPts[2][3];
    this->FormBifurcation(this->EndPt, this->StartPt,
                         parent->StartPt, parent->EndPt,
                         brother->EndPt, brother->StartPt,
                         this->StartPt,
                         width/2., vecs, triPts);

    // get vector towards top of cube
    double diverVec[3];
    vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
    vtkMath::Normalize(diverVec);

    // get front vector towards top of cube
    double frontVec[3];
    vtkMath::Cross(diverVec, vecs[1], frontVec);
    vtkMath::Normalize(frontVec);

    // flip if divergent is left or front
    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(frontVec, -1.0);

    // Get start wedge
    vtkNew(vtkPoints, startPoints);

    // Treat differently for different dirs
    if (this->BranchDir == LEFT || this->BranchDir == FRONT)
      this->GetWedge(triPts[1], this->StartPt, triPts[0], frontVec,
                     height, startPoints);
    else
      this->GetWedge(triPts[0], this->StartPt, triPts[1], frontVec,
                     height, startPoints);

    startPoints->GetPoint(0, finalPts[8]);
    startPoints->GetPoint(1, finalPts[5]);
    startPoints->GetPoint(2, finalPts[1]);
    startPoints->GetPoint(3, finalPts[9]);
    startPoints->GetPoint(4, finalPts[6]);
    startPoints->GetPoint(5, finalPts[2]);

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, triPts);

    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], frontVec);
    vtkMath::Normalize(frontVec);

    // Get end wedge
    vtkNew(vtkPoints, endPoints);

    // Treat differently for different dirs
    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
    {
      vtkMath::MultiplyScalar(frontVec, -1.0);
      this->GetWedge(triPts[0], this->EndPt, triPts[1], frontVec,
                     height, endPoints);
    }
    else
      this->GetWedge(triPts[1], this->EndPt, triPts[0], frontVec,
                     height, endPoints);

    endPoints->GetPoint(0, finalPts[3]);
    endPoints->GetPoint(1, finalPts[4]);
    endPoints->GetPoint(2, finalPts[0]);
    endPoints->GetPoint(3, finalPts[10]);
    endPoints->GetPoint(4, finalPts[11]);
    endPoints->GetPoint(5, finalPts[7]);

    int pI[12];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
      localPtIds->InsertNextTuple1(i);
    }

    vtkIdType face0[5] = {pI[7], pI[8], pI[5], pI[1], pI[0]};
    vtkIdType face1[5] = {pI[10], pI[9], pI[8], pI[7], pI[11]};
    vtkIdType face2[5] = {pI[3], pI[2], pI[6], pI[9], pI[10]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};

    allCells->InsertNextCell(5, face0);
    allCells->InsertNextCell(5, face1);
    allCells->InsertNextCell(5, face2);
    allCells->InsertNextCell(5, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (cubeType == 6)
  {
    int numPoints = 10;

    double finalPts[10][3];

    if (this->Parent->Children.size() == 2)
    {
      vtkSVCenterlineGCell *brother, *diver, *parent;
      parent = this->Parent;
      if (parent->Children[0]->Id == this->Id)
        brother = parent->Children[1];
      else
        brother = parent->Children[0];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      double vecs[3][3];
      double triPts[2][3];
      this->FormBifurcation(this->EndPt, this->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           this->StartPt,
                           width/2., vecs, triPts);

      // get diverging vector towards top of cube
      double diverVec[3];
      vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
      vtkMath::Normalize(diverVec);

      // get front vector
      double frontVec[3];
      vtkMath::Cross(diverVec, vecs[1], frontVec);
      vtkMath::Normalize(frontVec);

      // Get top wedge
      vtkNew(vtkPoints, startPoints);

      // Flip if diver branch is left or front
      if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(frontVec, -1.0);

      // Treat differently for different dirs
      if (this->BranchDir == LEFT || this->BranchDir == FRONT)
        this->GetWedge(triPts[1], this->StartPt, triPts[0], frontVec,
                       height, startPoints);
      else
        this->GetWedge(triPts[0], this->StartPt, triPts[1], frontVec,
                       height, startPoints);

      startPoints->GetPoint(0, finalPts[1]);
      startPoints->GetPoint(1, finalPts[2]);
      startPoints->GetPoint(2, finalPts[3]);
      startPoints->GetPoint(3, finalPts[6]);
      startPoints->GetPoint(4, finalPts[7]);
      startPoints->GetPoint(5, finalPts[8]);

      // Now get beginning points
      double endVec[3];
      vtkMath::Normalize(frontVec);
      vtkMath::Cross(frontVec, vecs[0], endVec);
      vtkMath::Normalize(endVec);
      //vtkMath::MultiplyScalar(endVec, width/2.);

      //// Get points to extrude up
      //double endPt0[3], endPt1[3];
      //vtkMath::Add(this->EndPt, endVec, endPt0);
      //vtkMath::MultiplyScalar(endVec, -1.0);
      //vtkMath::Add(this->EndPt, endVec, endPt1);

      //// extrude down
      //vtkMath::MultiplyScalar(frontVec, height/2.);
      //vtkMath::Add(endPt0, frontVec, finalPts[5]);
      //vtkMath::Add(endPt1, frontVec, finalPts[9]);

      //// extrude up
      //vtkMath::MultiplyScalar(frontVec, -1.0);
      //vtkMath::Add(endPt0, frontVec, finalPts[0]);
      //vtkMath::Add(endPt1, frontVec, finalPts[4]);

      // Get square from end pt
      vtkNew(vtkPoints, endPoints);
      this->GetSquare(this->EndPt, endVec, frontVec,
                      height, width, endPoints);

      endPoints->GetPoint(0, finalPts[5]);
      endPoints->GetPoint(1, finalPts[0]);
      endPoints->GetPoint(2, finalPts[4]);
      endPoints->GetPoint(3, finalPts[9]);
    }
    else
    {
      int crossBrother = -1;
      for (int i=0; i<3; i++)
      {
        if (i != this->Parent->AligningChild && i != this->Parent->DivergingChild)
          crossBrother = i;
      }

      vtkSVCenterlineGCell *brother, *diver, *parent;
      parent = this->Parent;
      double vecs[3][3];
      double topPts[2][3];

      if (parent->Children[parent->AligningChild]->Id == this->Id)
      {
        brother = parent->Children[crossBrother];
        diver = parent->Children[parent->AligningChild];
        parent  = parent->Children[parent->DivergingChild];

        // Get ending bifurcation points
        this->FormBifurcation(this->EndPt, this->StartPt,
                             parent->EndPt, parent->StartPt,
                             brother->EndPt, brother->StartPt,
                             this->StartPt,
                             width/2., vecs, topPts);
      }
      else
      {
        brother = parent->Children[parent->AligningChild];
        diver = parent->Children[parent->DivergingChild];

        // Get ending bifurcation points
        this->FormBifurcation(this->EndPt, this->StartPt,
                             parent->StartPt, parent->EndPt,
                             brother->EndPt, brother->StartPt,
                             this->StartPt,
                             width/2., vecs, topPts);
      }



      // get vector towards top of cube
      double diverVec[3];
      vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
      vtkMath::Normalize(diverVec);
      double frontVec[3];
      vtkMath::Cross(diverVec, vecs[1], frontVec);
      vtkMath::Normalize(frontVec);

      // extrude up
      double finalPts[10][3];
      vtkMath::MultiplyScalar(frontVec, height/2.);
      if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(frontVec, -1.0);

      vtkMath::Add(this->StartPt, frontVec, finalPts[2]);
      if (this->BranchDir == RIGHT || this->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], frontVec, finalPts[1]);
        vtkMath::Add(topPts[1], frontVec, finalPts[3]);
      }
      else
      {
        vtkMath::Add(topPts[0], frontVec, finalPts[3]);
        vtkMath::Add(topPts[1], frontVec, finalPts[1]);
      }

      // extrude down
      vtkMath::MultiplyScalar(frontVec, -1.0);
      vtkMath::Add(this->StartPt, frontVec, finalPts[7]);
      if (this->BranchDir == RIGHT || this->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], frontVec, finalPts[6]);
        vtkMath::Add(topPts[1], frontVec, finalPts[8]);
      }
      else
      {
        vtkMath::Add(topPts[0], frontVec, finalPts[8]);
        vtkMath::Add(topPts[1], frontVec, finalPts[6]);
      }

      // Now get beginning points
      double endVec[3];
      vtkMath::Normalize(frontVec);
      vtkMath::Cross(vecs[0], frontVec, endVec);
      vtkMath::Normalize(endVec);
      vtkMath::MultiplyScalar(endVec, width/2.);

      // Get points to extrude up
      double endPt0[3], endPt1[3];
      vtkMath::Add(this->EndPt, endVec, endPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(this->EndPt, endVec, endPt1);

      // extrude down
      vtkMath::MultiplyScalar(frontVec, height/2.);
      vtkMath::Add(endPt0, frontVec, finalPts[5]);
      vtkMath::Add(endPt1, frontVec, finalPts[9]);

      // extrude up
      vtkMath::MultiplyScalar(frontVec, -1.0);
      vtkMath::Add(endPt0, frontVec, finalPts[0]);
      vtkMath::Add(endPt1, frontVec, finalPts[4]);
    }

    int pI[10];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
      localPtIds->InsertNextTuple1(i);
    }
    vtkIdType face0[4] = {pI[5], pI[6], pI[1], pI[0]};
    vtkIdType face1[5] = {pI[9], pI[8], pI[7], pI[6], pI[5]};
    vtkIdType face2[4] = {pI[4], pI[3], pI[8], pI[9]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};
    vtkIdType face4[4] = {pI[5], pI[0], pI[4], pI[9]};

    allCells->InsertNextCell(4, face0);
    allCells->InsertNextCell(5, face1);
    allCells->InsertNextCell(4, face2);
    allCells->InsertNextCell(5, face3);
    allCells->InsertNextCell(4, face4);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
    groupIds->InsertNextTuple1(groupId);
    patchIds->InsertNextTuple1(5);
  }
  else if (cubeType == 7)
  {
    int numPoints = 10;

    vtkSVCenterlineGCell *brother, *diver, *parent;
    parent = this->Parent;
    if (parent->Children[0]->Id == this->Id)
      brother = parent->Children[1];
    else
      brother = parent->Children[0];
    diver = parent->Children[parent->DivergingChild];

    // Get ending bifurcation points
    double vecs[3][3];
    double topPts[2][3];
    this->FormBifurcation(this->EndPt, this->StartPt,
                         parent->StartPt, parent->EndPt,
                         brother->EndPt, brother->StartPt,
                         this->StartPt,
                         width/2., vecs, topPts);

    // get vector towards top of cube
    double diverVec[3];
    vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
    vtkMath::Normalize(diverVec);
    double frontVec[3];
    vtkMath::Cross(diverVec, vecs[1], frontVec);
    vtkMath::Normalize(frontVec);

    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(frontVec, -1.0);

    // extrude up
    double finalPts[10][3];
    vtkMath::MultiplyScalar(frontVec, height/2.);
    vtkMath::Add(this->StartPt, frontVec, finalPts[4]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], frontVec, finalPts[7]);
      vtkMath::Add(topPts[1], frontVec, finalPts[1]);
    }
    else
    {
      vtkMath::Add(topPts[0], frontVec, finalPts[1]);
      vtkMath::Add(topPts[1], frontVec, finalPts[7]);
    }

    // extrude down
    vtkMath::MultiplyScalar(frontVec, -1.0);
    vtkMath::Add(this->StartPt, frontVec, finalPts[5]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], frontVec, finalPts[8]);
      vtkMath::Add(topPts[1], frontVec, finalPts[2]);
    }
    else
    {
      vtkMath::Add(topPts[0], frontVec, finalPts[2]);
      vtkMath::Add(topPts[1], frontVec, finalPts[8]);
    }

    // Now get beginning points
    double endVec[3];
    vtkMath::Normalize(frontVec);
    vtkMath::Cross(vecs[0], frontVec, endVec);
    vtkMath::Normalize(endVec);
    vtkMath::MultiplyScalar(endVec, width/2.);

    // Get points to extrude up
    double endPt0[3], endPt1[3];
    vtkMath::Add(this->EndPt, endVec, endPt0);
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->EndPt, endVec, endPt1);

    // extrude down
    vtkMath::MultiplyScalar(frontVec, height/2.);
    vtkMath::Add(endPt0, frontVec, finalPts[9]);
    vtkMath::Add(endPt1, frontVec, finalPts[3]);

    // extrude up
    vtkMath::MultiplyScalar(frontVec, -1.0);
    vtkMath::Add(endPt0, frontVec, finalPts[6]);
    vtkMath::Add(endPt1, frontVec, finalPts[0]);

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);

    if (points->GetNumberOfPoints() != 10)
    {
      fprintf(stderr,"Must be given 10 allPoints for cube type 7!\n");
      return SV_ERROR;
    }
    int pI[10];
    for (int i=0; i<points->GetNumberOfPoints(); i++)
    {
      pI[i] = allPoints->InsertNextPoint(points->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    vtkIdType face0[5] = {pI[6], pI[7], pI[4], pI[1], pI[0]};
    vtkIdType face1[4] = {pI[9], pI[8], pI[7], pI[6]};
    vtkIdType face2[5] = {pI[3], pI[2], pI[5], pI[8], pI[9]};
    vtkIdType face3[4] = {pI[0], pI[1], pI[2], pI[3]};
    vtkIdType face4[4] = {pI[6], pI[0], pI[3], pI[9]};

    allCells->InsertNextCell(5, face0);
    allCells->InsertNextCell(4, face1);
    allCells->InsertNextCell(5, face2);
    allCells->InsertNextCell(4, face3);
    allCells->InsertNextCell(4, face4);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
    groupIds->InsertNextTuple1(groupId);
    patchIds->InsertNextTuple1(5);
  }
  else if (cubeType == 8)
  {
    int numPoints = 10;
    double vecs[3][3];
    double endPts[2][3];

    int crossChild;
    for (int i=0; i<3; i++)
    {
      if (i != this->AligningChild && i != this->DivergingChild)
        crossChild = i;
    }
    vtkSVCenterlineGCell *cross = this->Children[crossChild];
    vtkSVCenterlineGCell *diver = this->Children[this->DivergingChild];

    this->FormBifurcation(this->StartPt, this->EndPt,
                          diver->EndPt, diver->StartPt,
                          cross->EndPt, cross->StartPt,
                          this->EndPt,
                          width/2., vecs, endPts);

    // get vector towards top of cube
    double frontVec[3];
    vtkMath::Cross(vecs[1], vecs[0], frontVec);
    vtkMath::Normalize(frontVec);

    // extrude up
    double finalPts[10][3];
    vtkMath::MultiplyScalar(frontVec, height/2.);
    vtkMath::Add(this->EndPt, frontVec, finalPts[4]);
    vtkMath::Add(endPts[0], frontVec, finalPts[0]);
    vtkMath::Add(endPts[1], frontVec, finalPts[3]);

    // extrude down
    vtkMath::MultiplyScalar(frontVec, -1.0);
    vtkMath::Add(this->EndPt, frontVec, finalPts[9]);
    vtkMath::Add(endPts[0], frontVec, finalPts[5]);
    vtkMath::Add(endPts[1], frontVec, finalPts[8]);

    // Now get beginning points
    double endVec[3];
    vtkMath::Normalize(frontVec);
    vtkMath::Cross(frontVec, vecs[0], endVec);
    vtkMath::Normalize(endVec);
    vtkMath::MultiplyScalar(endVec, width/2.);

    // Get points to extrude up
    double topPt0[3], topPt1[3];
    vtkMath::Add(this->StartPt, endVec, topPt0);
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->StartPt, endVec, topPt1);

    // extrude down
    vtkMath::MultiplyScalar(frontVec, height/2.);
    vtkMath::Add(topPt0, frontVec, finalPts[6]);
    vtkMath::Add(topPt1, frontVec, finalPts[7]);

    // extrude up
    vtkMath::MultiplyScalar(frontVec, -1.0);
    vtkMath::Add(topPt0, frontVec, finalPts[1]);
    vtkMath::Add(topPt1, frontVec, finalPts[2]);

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);

    if (points->GetNumberOfPoints() != 10)
    {
      fprintf(stderr,"Must be given 10 allPoints for cube type 1!\n");
      return SV_ERROR;
    }
    int pI[10];
    for (int i=0; i<points->GetNumberOfPoints(); i++)
    {
      pI[i] = allPoints->InsertNextPoint(points->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    vtkIdType face0[4] = {pI[5], pI[6], pI[1], pI[0]};
    vtkIdType face1[5] = {pI[8], pI[7], pI[6], pI[5], pI[9]};
    vtkIdType face2[4] = {pI[3], pI[2], pI[7], pI[8]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};
    vtkIdType face4[4] = {pI[2], pI[1], pI[6], pI[7]};

    allCells->InsertNextCell(4, face0);
    allCells->InsertNextCell(5, face1);
    allCells->InsertNextCell(4, face2);
    allCells->InsertNextCell(5, face3);
    allCells->InsertNextCell(4, face4);

    for (int i=0; i<5; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }

  return SV_OK;
}

// ----------------------
// FormBifurcation
// ----------------------
int vtkSVCenterlineGCell::FormBifurcation(const double pt0[3], const double pt1[3],
                                          const double pt2[3], const double pt3[3],
                                          const double pt4[3], const double pt5[3],
                                          const double centerPt[3],
                                          const double factor,
                                          double vecs[3][3],
                                          double returnPts[2][3])
{
  vtkMath::Subtract(pt0, pt1, vecs[0]);
  vtkMath::Subtract(pt2, pt3, vecs[1]);
  vtkMath::Subtract(pt4, pt5, vecs[2]);
  for (int j=0; j<3; j++)
    vtkMath::Normalize(vecs[j]);

  // Get ending bifurcation points
  this->GetBifurcationPoint(centerPt, vecs[0], vecs[1], vecs[2], factor, returnPts[0]);
  this->GetBifurcationPoint(centerPt, vecs[0], vecs[2], vecs[1], factor, returnPts[1]);

  return SV_OK;
}

// ----------------------
// GetBifurcationPoint
// ----------------------
int vtkSVCenterlineGCell::GetBifurcationPoint(const double startPt[3],
                                             const double vec0[3],
                                             const double vec1[3],
                                             const double vec2[3],
                                             const double factor,
                                             double returnPt[3])
{
    // Get vector in between the two
    double midVec[3];
    vtkMath::Add(vec0, vec1, midVec);
    vtkMath::Normalize(midVec);

    // Get the angle between
    double perpVec[3];
    vtkMath::Cross(vec0, vec1, perpVec);
    double ang = atan2(vtkMath::Norm(perpVec), vtkMath::Dot(vec0, vec1));

    double dotCheck = vtkMath::Dot(midVec, vec2);

    double midLength;
    if (dotCheck > 0)
    {
      vtkMath::MultiplyScalar(midVec, -1.0);
      midLength = factor / ( cos(ang/2.));
    }
    else
      midLength = factor / ( sin(ang/2.));

    vtkMath::MultiplyScalar(midVec, midLength);
    vtkMath::Add(startPt, midVec, returnPt);
}

// ----------------------
// GetSquare
// ----------------------
int vtkSVCenterlineGCell::GetSquare(const double startPt[3],
                                    const double vec0[3], const double vec1[3],
                                    const double height, const double width,
                                    vtkPoints *points)
{
  double workVec0[3];
  double workVec1[3];

  // copy vectors
  for (int i=0; i<3; i++)
  {
    workVec0[i] = vec0[i];
    workVec1[i] = vec1[i];
  }

  // Get two dir vectors
  vtkMath::Normalize(workVec0);
  vtkMath::Normalize(workVec1);
  vtkMath::MultiplyScalar(workVec0, width/2.);
  vtkMath::MultiplyScalar(workVec1, height/2.);

  double pts[4][3];
  for (int i=0; i<3; i++)
  {
    pts[0][i] = startPt[i] + workVec0[i] - workVec1[i]; // bottom right
    pts[1][i] = startPt[i] + workVec0[i] + workVec1[i]; // top right
    pts[2][i] = startPt[i] - workVec0[i] + workVec1[i]; // top left
    pts[3][i] = startPt[i] - workVec0[i] - workVec1[i]; // bottom left
  }

  points->SetNumberOfPoints(4);
  for (int i=0; i<4; i++)
    points->SetPoint(i, pts[i]);

  return SV_OK;
}

// ----------------------
// GetWedge
// ----------------------
int vtkSVCenterlineGCell::GetWedge(const double pt0[3], const double pt1[3],
                                   const double pt2[3],
                                   const double vec0[3],
                                   const double height,
                                   vtkPoints *points)
{
  double workVec0[3];

  // copy vector
  for (int i=0; i<3; i++)
    workVec0[i] = vec0[i];

  // Get two dir vectors
  vtkMath::Normalize(workVec0);
  vtkMath::MultiplyScalar(workVec0, height/2.);

  double pts[6][3];
  for (int i=0; i<3; i++)
  {
    pts[0][i] = pt0[i] + workVec0[i];
    pts[1][i] = pt1[i] + workVec0[i];
    pts[2][i] = pt2[i] + workVec0[i];
    pts[3][i] = pt0[i] - workVec0[i];
    pts[4][i] = pt1[i] - workVec0[i];
    pts[5][i] = pt2[i] - workVec0[i];
  }

  points->SetNumberOfPoints(6);
  for (int i=0; i<6; i++)
    points->SetPoint(i, pts[i]);

  return SV_OK;
}
