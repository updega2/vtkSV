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

#include "vtkCellArray.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"
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

// ----------------------
// GetBeginningType
// ----------------------
int vtkSVCenterlineGCell::GetBeginningType(int &beginningType, int &splitType)
{
  if (this->Parent == NULL)
  {
    splitType     = ZERO;
    beginningType = NONE;
  }
  else if (this->Parent->Children.size() == 2 || this->Parent->Children.size() == 3)
  {
    if (this->Parent->Children.size() == 2)
      splitType = BI;
    else
      splitType = TRI;

    if ((this->Parent->BranchDir + this->BranchDir)%2 == 0)
    {
      if (this->Parent->BranchDir == RIGHT || this->Parent->BranchDir == LEFT)
        beginningType = VERT_WEDGE;
      else
        beginningType = HORZ_WEDGE;
    }
    else
    {
      if (this->Parent->BranchDir == RIGHT || this->Parent->BranchDir == LEFT)
        beginningType = HORZ_WEDGE;
      else
        beginningType = VERT_WEDGE;
    }
  }
  else
  {
    beginningType = NOTHANDLED;
    splitType     = TOOMANY;
    fprintf(stdout,"Cannot currently handle more than a trifurcation\n");
    return SV_ERROR;
  }

  return SV_OK;
}

// ----------------------
// GetEndType
// ----------------------
int vtkSVCenterlineGCell::GetEndType(int &endType, int &splitType)
{
  if (this->Children.size() == 0)
  {
    endType   = NONE;
    splitType = ZERO;
  }
  else if (this->Children.size() == 2 || this->Children.size() == 3)
  {
    if (this->Children.size() == 2)
      splitType = BI;
    else
      splitType = TRI;

    if ((this->BranchDir + this->Children[this->DivergingChild]->BranchDir)%2 == 0)
    {
      if (this->BranchDir == RIGHT || this->BranchDir == LEFT)
        endType = VERT_WEDGE;
      else
        endType = HORZ_WEDGE;
    }
    else
    {
      if (this->BranchDir == RIGHT || this->BranchDir == LEFT)
        endType = HORZ_WEDGE;
      else
        endType = VERT_WEDGE;
    }
  }
  else
  {
    endType   = NOTHANDLED;
    splitType = TOOMANY;
    fprintf(stdout,"Cannot currently handle more than a trifurcation\n");
    return SV_ERROR;
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
    fprintf(stderr,"Trifurcation wrongly identified, number of children is actually %d\n", numChildren);
  }

  fprintf(stdout,"TRIFURCATION, TELL ME THIS DIR:            %d\n", this->BranchDir);
  fprintf(stdout,"TRIFURCATION, TELL ME DIVERGING CHILD DIR: %d\n", this->Children[this->DivergingChild]->BranchDir);
  fprintf(stdout,"TRIFURCATION, TELL ME ALIGNING CHILD DIR:  %d\n", this->Children[this->AligningChild]->BranchDir);

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
      fprintf(stdout,"TRIFURCATION, Would like to see compare, dir: %d, dot: %.6f\n", j, compare);
      if (compare > maxDot)
      {
        maxDot = compare;
        maxDir = j;
      }
    }
    fprintf(stdout,"Child %d max dir is %d\n", this->Children[i]->GroupId, maxDir);
    double dotCheck = vtkMath::Dot(this->Children[i]->RefDirs[maxDir], this->Children[i]->BranchVec);
    fprintf(stdout,"And the dot check is %.6f\n", dotCheck);
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
        frontCount++;
    }
  }

  fprintf(stdout,"MY COUNT, RIGHT: %d, BACK: %d, LEFT: %d, FRONT: %d, DOWN: %d\n", rightCount, backCount, leftCount, frontCount, downCount);
  if ((rightCount + leftCount + downCount) == 3 ||
      (backCount + frontCount +downCount) == 3)
  {
    if (this->Parent == NULL)
    {
      type = 1;
    }
    else
    {
      if ((this->Parent->BranchDir + this->BranchDir)%2 == 0)
      {
        if (this->Parent->BranchDir == BACK || this->Parent->BranchDir == FRONT)
        {
          if ((this->BranchDir + this->Children[this->DivergingChild]->BranchDir)%2 == 0)
            type = 3;
          else
            type = 5;
        }
        else
        {
          if ((this->BranchDir + this->Children[this->DivergingChild]->BranchDir)%2 == 0)
            type = 2;
          else
            type = 4;
        }
      }
      else
      {
        if ((this->BranchDir + this->Children[this->DivergingChild]->BranchDir)%2 == 0)
          type = 3;
        else
          type = 5;
      }
    }
  }
  else if (rightCount + frontCount + downCount == 3)
  {
    type =8;
  }
  else
  {
    type = -1;
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
  fprintf(stdout,"CUBE GROUP ID: %d\n", this->GroupId);

  int groupId = this->GroupId;

  vtkSVCenterlineGCell *brother, *diver, *parent, *align, *cDiver, *otherBrother;
  parent = this->Parent;
  if (parent != NULL)
  {
    if (parent->Children[0]->Id == this->Id)
      brother = parent->Children[1];
    else
      brother = parent->Children[0];
    diver = parent->Children[parent->DivergingChild];
  }

  if (this->Children.size() > 0)
  {
    align = this->Children[this->AligningChild];
    cDiver = this->Children[this->DivergingChild];
  }

  // Get beginning type
  int begType, begSplitType;
  this->GetBeginningType(begType, begSplitType);

  // Do beginning
  double begVec[3];
  double vecs[3][3];
  double triPts[2][3];
  vtkNew(vtkPoints, begPoints);

  if (begSplitType == ZERO)
  {
    fprintf(stdout,"BEG NONE\n");
    // Get top square from start pt
    vtkMath::Cross(this->RefDirs[1], this->RefDirs[0], begVec);
    vtkMath::Normalize(begVec);

    this->GetSquare(this->StartPt, this->RefDirs[1], begVec,
                    height, width, begPoints);

  }
  else if (begSplitType == BI)
  {
    fprintf(stdout,"BEG BI\n");

    // Get ending bifurcation points
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
    vtkMath::Cross(diverVec, vecs[1], begVec);
    vtkMath::Normalize(begVec);

    // Flip if diver branch is left or front
    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(begVec, -1.0);

    // Treat differently for different dirs
    if (this->BranchDir == LEFT || this->BranchDir == FRONT)
      this->GetWedge(triPts[1], this->StartPt, triPts[0], begVec,
                     height, begPoints);
    else
      this->GetWedge(triPts[0], this->StartPt, triPts[1], begVec,
                     height, begPoints);

  }
  else if (begSplitType == TRI)
  {
    fprintf(stdout,"BEG TRI\n");

    int myLoc;
    for (int i=0; i<3; i++)
    {
      if (parent->Children[i]->Id == this->Id)
        myLoc = i;
    }

    if (myLoc == 0)
    {
      brother = parent->Children[parent->AligningChild];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      this->FormBifurcation(this->EndPt, this->StartPt,
                            parent->StartPt, parent->EndPt,
                            brother->EndPt, brother->StartPt,
                            this->StartPt,
                            width/2., vecs, triPts);
    }
    else if (myLoc == parent->Children.size() - 1)
    {
      brother = parent->Children[parent->AligningChild];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      this->FormBifurcation(this->EndPt, this->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           this->StartPt,
                           width/2., vecs, triPts);
    }
    else
    {
      brother = parent->Children[myLoc - 1];
      diver = parent->Children[parent->AligningChild];
      otherBrother  = parent->Children[myLoc + 1];

      // Get ending bifurcation points
      this->FormBifurcation(this->EndPt, this->StartPt,
                           brother->EndPt, brother->StartPt,
                           otherBrother->EndPt, otherBrother->StartPt,
                           this->StartPt,
                           width/2., vecs, triPts);
    }

    // get vector towards top of cube
    double diverVec[3];
    vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
    vtkMath::Normalize(diverVec);

    vtkMath::Cross(diverVec, vecs[1], begVec);
    vtkMath::Normalize(begVec);

    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(begVec, -1.0);

    if (this->BranchDir == LEFT || this->BranchDir == FRONT)
      this->GetWedge(triPts[1], this->StartPt, triPts[0], begVec,
                     height, begPoints);
    else
      this->GetWedge(triPts[0], this->StartPt, triPts[1], begVec,
                     height, begPoints);

  }
  else
  {
    fprintf(stdout,"BEG NOT_HANDLED\n");
  }

  // Get end type
  int endType, endSplitType;
  this->GetEndType(endType, endSplitType);

  // Do end
  vtkNew(vtkPoints, endPoints);

  if (endSplitType == ZERO)
  {
    fprintf(stdout,"END NONE\n");

    // Get square from end pt
    if (begSplitType == ZERO)
      this->GetSquare(this->EndPt, this->RefDirs[1], begVec,
                      height, width, endPoints);
    else
    {
      double endVec[3];
      if ((this->BranchDir)%2 == 0)
      {
        vtkMath::Cross(begVec, vecs[0], endVec);
        vtkMath::Normalize(endVec);
        this->GetSquare(this->EndPt, endVec, begVec,
                        height, width, endPoints);
      }
      else
      {
        vtkMath::Cross(vecs[0], begVec, endVec);
        vtkMath::Normalize(endVec);
        this->GetSquare(this->EndPt, begVec, endVec,
                        height, width, endPoints);
      }
    }
  }
  else if (endSplitType == BI)
  {
    fprintf(stdout,"END BI\n");

    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, triPts);


    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], begVec);
    vtkMath::Normalize(begVec);

    // treat differently for different dirs
    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
    {
      vtkMath::MultiplyScalar(begVec, -1.0);
      this->GetWedge(triPts[0], this->EndPt, triPts[1], begVec,
                     height, endPoints);
    }
    else
      this->GetWedge(triPts[1], this->EndPt, triPts[0], begVec,
                     height, endPoints);

  }
  else if (endSplitType == TRI)
  {
    fprintf(stdout,"END TRI\n");

    int crossChild;
    for (int i=0; i<3; i++)
    {
      if (i != this->AligningChild && i != this->DivergingChild)
        crossChild = i;
    }
    vtkSVCenterlineGCell *cross = this->Children[crossChild];
    diver = this->Children[this->DivergingChild];

    double endPts[2][3];
    this->FormBifurcation(this->StartPt, this->EndPt,
                          diver->EndPt, diver->StartPt,
                          cross->EndPt, cross->StartPt,
                          this->EndPt,
                          width/2., vecs, endPts);

    // TEMPORARYRYRY
    if (diver->BranchDir == RIGHT || diver->BranchDir == BACK)
    {
      // get vector towards top of cube
      double frontVec[3];
      vtkMath::Cross(vecs[1], vecs[0], frontVec);
      vtkMath::Normalize(frontVec);

      this->GetWedge(endPts[1], this->EndPt, endPts[0], frontVec,
                      height, endPoints);
    }
    else if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
    {
      // get vector towards top of cube
      double frontVec[3];
      vtkMath::Cross(vecs[0], vecs[1], frontVec);
      vtkMath::Normalize(frontVec);

      this->GetWedge(endPts[0], this->EndPt, endPts[1], frontVec,
                      height, endPoints);
    }
    else
    {
      fprintf(stderr,"NEED TO ADD RULE!!\n");
    }
  }
  else
  {
    fprintf(stderr,"END NOT_HANDLED\n");
  }

  if (begType == NONE && endType == NONE)
  {
    int numPoints = 8;

    double finalPts[8][3];

    begPoints->GetPoint(0, finalPts[5]);
    begPoints->GetPoint(1, finalPts[1]);
    begPoints->GetPoint(2, finalPts[2]);
    begPoints->GetPoint(3, finalPts[6]);

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
  else if (begType == NONE && endType == VERT_WEDGE)
  {
    int numPoints = 10;

    double finalPts[10][3];

    begPoints->GetPoint(0, finalPts[6]);
    begPoints->GetPoint(1, finalPts[1]);
    begPoints->GetPoint(2, finalPts[2]);
    begPoints->GetPoint(3, finalPts[7]);

    endPoints->GetPoint(0, finalPts[3]);
    endPoints->GetPoint(1, finalPts[4]);
    endPoints->GetPoint(2, finalPts[0]);
    endPoints->GetPoint(3, finalPts[8]);
    endPoints->GetPoint(4, finalPts[9]);
    endPoints->GetPoint(5, finalPts[5]);

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
  else if (begType == VERT_WEDGE && endType == VERT_WEDGE)
  {
    int numPoints = 12;

    double finalPts[12][3];

    begPoints->GetPoint(0, finalPts[1]);
    begPoints->GetPoint(1, finalPts[2]);
    begPoints->GetPoint(2, finalPts[3]);
    begPoints->GetPoint(3, finalPts[7]);
    begPoints->GetPoint(4, finalPts[8]);
    begPoints->GetPoint(5, finalPts[9]);

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
  else if (begType == HORZ_WEDGE && endType == HORZ_WEDGE)
  {
    int numPoints = 12;

    double finalPts[12][3];

    begPoints->GetPoint(0, finalPts[9]);
    begPoints->GetPoint(1, finalPts[5]);
    begPoints->GetPoint(2, finalPts[1]);
    begPoints->GetPoint(3, finalPts[10]);
    begPoints->GetPoint(4, finalPts[6]);
    begPoints->GetPoint(5, finalPts[2]);

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
  else if (begType == VERT_WEDGE && endType == HORZ_WEDGE)
  {
    int numPoints = 12;

    double finalPts[12][3];

    begPoints->GetPoint(0, finalPts[1]);
    begPoints->GetPoint(1, finalPts[2]);
    begPoints->GetPoint(2, finalPts[3]);
    begPoints->GetPoint(3, finalPts[8]);
    begPoints->GetPoint(4, finalPts[9]);
    begPoints->GetPoint(5, finalPts[10]);

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
  else if (begType == HORZ_WEDGE && endType == VERT_WEDGE)
  {
    int numPoints = 12;

    double finalPts[12][3];

    begPoints->GetPoint(0, finalPts[8]);
    begPoints->GetPoint(1, finalPts[5]);
    begPoints->GetPoint(2, finalPts[1]);
    begPoints->GetPoint(3, finalPts[9]);
    begPoints->GetPoint(4, finalPts[6]);
    begPoints->GetPoint(5, finalPts[2]);

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
  else if (begType == VERT_WEDGE && endType == NONE)
  {
    int numPoints = 10;

    double finalPts[10][3];

    if (this->Parent->Children.size() == 2)
    {
      begPoints->GetPoint(0, finalPts[1]);
      begPoints->GetPoint(1, finalPts[2]);
      begPoints->GetPoint(2, finalPts[3]);
      begPoints->GetPoint(3, finalPts[6]);
      begPoints->GetPoint(4, finalPts[7]);
      begPoints->GetPoint(5, finalPts[8]);

      endPoints->GetPoint(0, finalPts[5]);
      endPoints->GetPoint(1, finalPts[0]);
      endPoints->GetPoint(2, finalPts[4]);
      endPoints->GetPoint(3, finalPts[9]);
    }
    else
    {
      begPoints->GetPoint(0, finalPts[1]);
      begPoints->GetPoint(1, finalPts[2]);
      begPoints->GetPoint(2, finalPts[3]);
      begPoints->GetPoint(3, finalPts[6]);
      begPoints->GetPoint(4, finalPts[7]);
      begPoints->GetPoint(5, finalPts[8]);

      endPoints->GetPoint(0, finalPts[5]);
      endPoints->GetPoint(1, finalPts[0]);
      endPoints->GetPoint(2, finalPts[4]);
      endPoints->GetPoint(3, finalPts[9]);
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
  else if (begType == HORZ_WEDGE && endType == NONE)
  {
    int numPoints = 10;

    double finalPts[10][3];

    begPoints->GetPoint(0, finalPts[7]);
    begPoints->GetPoint(1, finalPts[4]);
    begPoints->GetPoint(2, finalPts[1]);
    begPoints->GetPoint(3, finalPts[8]);
    begPoints->GetPoint(4, finalPts[5]);
    begPoints->GetPoint(5, finalPts[2]);

    endPoints->GetPoint(0, finalPts[6]);
    endPoints->GetPoint(1, finalPts[0]);
    endPoints->GetPoint(2, finalPts[3]);
    endPoints->GetPoint(3, finalPts[9]);

    int pI[10];
    for (int i=0; i<numPoints; i++)
    {
      pI[i] = allPoints->InsertNextPoint(finalPts[i]);
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
  //// TODO: CHECK ON THIS BECAUSE IM NOT SO SURE THAT THIS IS BAD,
  //// I THINK I NEED IT BUT TBD LATER
  //if (dotCheck > 0)
  //{
  //  vtkMath::MultiplyScalar(midVec, -1.0);
  //  midLength = factor / ( cos(ang/2.));
  //}
  midLength = factor / ( sin(ang/2.));

  vtkMath::MultiplyScalar(midVec, midLength);
  vtkMath::Add(startPt, midVec, returnPt);

  return SV_OK;
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
