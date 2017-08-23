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

  // -----------------------------------------------------------------------
  // FACE 3
  // -----------------------------------------------------------------------
  std::vector<int> face3PtIds;

  // END
  if (endType == NONE)
    face3PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(1)));
  else if (endType == VERT_WEDGE)
    face3PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  else if (endType == HORZ_WEDGE)
    face3PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));

  // BEG
  if (begType == NONE)
  {
    face3PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(1)));
    face3PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
  }
  else if (begType == VERT_WEDGE)
  {
    face3PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
    face3PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(1)));
    face3PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
  }
  else if (begType == HORZ_WEDGE)
  {
    face3PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
    face3PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(5)));
  }

  // END
  if (endType == NONE)
    face3PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  else if (endType == VERT_WEDGE)
  {
    face3PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));
    face3PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(1)));
  }
  else if (endType == HORZ_WEDGE)
    face3PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 1
  // -----------------------------------------------------------------------
  std::vector<int> face1PtIds;

  // END
  if (endType == NONE)
    face1PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  else if (endType == VERT_WEDGE)
    face1PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  else if (endType == HORZ_WEDGE)
    face1PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(5)));

  // BEG
  if (begType == NONE)
  {
    face1PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
    face1PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
  }
  else if (begType == VERT_WEDGE)
  {
    face1PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(5)));
    face1PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(4)));
    face1PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
  }
  else if (begType == HORZ_WEDGE)
  {
    face1PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
    face1PtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
  }

  // END
  if (endType == NONE)
    face1PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));
  else if (endType == VERT_WEDGE)
  {
    face1PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(5)));
    face1PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(4)));
  }
  else if (endType == HORZ_WEDGE)
    face1PtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // Now interior points if they are there!
  // -----------------------------------------------------------------------
  std::vector<int> endInterPtIds;
  std::vector<int> begInterPtIds;

  // END
  if (endType == HORZ_WEDGE)
  {
    endInterPtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(1)));
    endInterPtIds.push_back(allPoints->InsertNextPoint(endPoints->GetPoint(4)));
  }

  // BEG
  if (begType == HORZ_WEDGE)
  {
    begInterPtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(1)));
    begInterPtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(4)));
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 0
  // -----------------------------------------------------------------------
  std::vector<int> face0PtIds;

  // END
  if (endType == NONE)
  {
    face0PtIds.push_back(face1PtIds[face1PtIds.size()-1]);
    face0PtIds.push_back(face1PtIds[face1PtIds.size()-2]);
  }
  else if (endType == VERT_WEDGE)
  {
    face0PtIds.push_back(face1PtIds[face1PtIds.size()-2]);
    face0PtIds.push_back(face1PtIds[face1PtIds.size()-3]);
  }
  else if (endType == HORZ_WEDGE)
  {
    face0PtIds.push_back(face1PtIds[face1PtIds.size()-1]);
    face0PtIds.push_back(face1PtIds[face1PtIds.size()-2]);
  }

  // INT
  if (begType == HORZ_WEDGE)
    face0PtIds.push_back(begInterPtIds[0]);

  // BEG
  if (begType == NONE)
  {
    face0PtIds.push_back(face3PtIds[1]);
    face0PtIds.push_back(face3PtIds[0]);
  }
  else if (begType == VERT_WEDGE)
  {
    face0PtIds.push_back(face3PtIds[1]);
    face0PtIds.push_back(face3PtIds[0]);
  }
  else if (begType == HORZ_WEDGE)
  {
    face0PtIds.push_back(face3PtIds[1]);
    face0PtIds.push_back(face3PtIds[0]);
  }

  // INT
  if (endType == HORZ_WEDGE)
    face0PtIds.push_back(endInterPtIds[0]);
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 2
  // -----------------------------------------------------------------------
  std::vector<int> face2PtIds;

  // END
  if (endType == NONE)
  {
    face2PtIds.push_back(face3PtIds[face3PtIds.size()-1]);
    face2PtIds.push_back(face3PtIds[face3PtIds.size()-2]);
  }
  else if (endType == VERT_WEDGE)
  {
    face2PtIds.push_back(face3PtIds[face3PtIds.size()-2]);
    face2PtIds.push_back(face3PtIds[face3PtIds.size()-3]);
  }
  else if (endType == HORZ_WEDGE)
  {
    face2PtIds.push_back(face3PtIds[face3PtIds.size()-1]);
    face2PtIds.push_back(face3PtIds[face3PtIds.size()-2]);
  }

  // INT
  if (begType == HORZ_WEDGE)
    face2PtIds.push_back(begInterPtIds[1]);

  // BEG
  if (begType == NONE)
  {
    face2PtIds.push_back(face1PtIds[1]);
    face2PtIds.push_back(face1PtIds[0]);
  }
  else if (begType == VERT_WEDGE)
  {
    face2PtIds.push_back(face1PtIds[1]);
    face2PtIds.push_back(face1PtIds[0]);
  }
  else if (begType == HORZ_WEDGE)
  {
    face2PtIds.push_back(face1PtIds[1]);
    face2PtIds.push_back(face1PtIds[0]);
  }

  // INT
  if (endType == HORZ_WEDGE)
    face2PtIds.push_back(endInterPtIds[1]);
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 4, if there
  // -----------------------------------------------------------------------
  std::vector<int> face4PtIds;

  // BEG
  if (begType == NONE)
  {
    face4PtIds.push_back(face3PtIds[2]);
    face4PtIds.push_back(face3PtIds[1]);
    face4PtIds.push_back(face1PtIds[2]);
    face4PtIds.push_back(face1PtIds[1]);
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 5, if there
  // -----------------------------------------------------------------------
  std::vector<int> face5PtIds;

  // BEG
  if (endType == NONE)
  {
    face5PtIds.push_back(face1PtIds[face1PtIds.size()-1]);
    face5PtIds.push_back(face3PtIds[0]);
    face5PtIds.push_back(face3PtIds[face3PtIds.size()-1]);
    face5PtIds.push_back(face1PtIds[0]);
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // LocalPtIds
  // -----------------------------------------------------------------------
  int numPoints = face3PtIds.size() + face1PtIds.size() + endInterPtIds.size() + begInterPtIds.size();
  fprintf(stdout,"NUMBER OF POINTS: %d\n", numPoints);
  for (int i=0; i<numPoints; i++)
    localPtIds->InsertNextTuple1(i);
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // Add faces
  // -----------------------------------------------------------------------
  // FACE 0
  int f0npts = face0PtIds.size();
  vtkIdType *face0 = new vtkIdType[f0npts];
  for (int i=0; i<f0npts; i++)
    face0[i] = face0PtIds[i];
  allCells->InsertNextCell(f0npts, face0);
  groupIds->InsertNextTuple1(groupId);
  patchIds->InsertNextTuple1(0);
  delete [] face0;

  // FACE 1
  int f1npts = face1PtIds.size();
  vtkIdType *face1 = new vtkIdType[f1npts];
  for (int i=0; i<f1npts; i++)
    face1[i] = face1PtIds[i];
  allCells->InsertNextCell(f1npts, face1);
  groupIds->InsertNextTuple1(groupId);
  patchIds->InsertNextTuple1(1);
  delete [] face1;

  // FACE 2
  int f2npts = face2PtIds.size();
  vtkIdType *face2 = new vtkIdType[f2npts];
  for (int i=0; i<f2npts; i++)
    face2[i] = face2PtIds[i];
  allCells->InsertNextCell(f2npts, face2);
  groupIds->InsertNextTuple1(groupId);
  patchIds->InsertNextTuple1(2);
  delete [] face2;

  // FACE 3
  int f3npts = face3PtIds.size();
  vtkIdType *face3 = new vtkIdType[f3npts];
  for (int i=0; i<f3npts; i++)
    face3[i] = face3PtIds[i];
  allCells->InsertNextCell(f3npts, face3);
  groupIds->InsertNextTuple1(groupId);
  patchIds->InsertNextTuple1(3);
  delete [] face3;

  // FACE 4
  int f4npts = face4PtIds.size();
  if (f4npts > 0)
  {
    vtkIdType *face4 = new vtkIdType[f4npts];
    for (int i=0; i<f4npts; i++)
      face4[i] = face4PtIds[i];
    allCells->InsertNextCell(f4npts, face4);
    groupIds->InsertNextTuple1(groupId);
    patchIds->InsertNextTuple1(4);
    delete [] face4;
  }

  // FACE 5
  int f5npts = face5PtIds.size();
  if (f5npts > 0)
  {
    vtkIdType *face5 = new vtkIdType[f5npts];
    for (int i=0; i<f5npts; i++)
      face5[i] = face5PtIds[i];
    allCells->InsertNextCell(f5npts, face5);
    groupIds->InsertNextTuple1(groupId);
    patchIds->InsertNextTuple1(5);
    delete [] face5;
  }
  // -----------------------------------------------------------------------

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
