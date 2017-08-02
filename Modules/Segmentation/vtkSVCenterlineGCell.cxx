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
int vtkSVCenterlineGCell::GetCubePoints(vtkPoints *points)
{
  int cubeType;
  this->GetCubeType(cubeType);

  if (cubeType == 0)
  {
    int numPoints = 8;

    double finalPts[8][3];

    // get vector towards top of cube
    double workVec[3];
    vtkMath::Cross(this->RefDirs[1], this->RefDirs[0], workVec);
    vtkMath::Normalize(workVec);

    // Get beginning points
    double endVec[3];
    vtkMath::Normalize(workVec);
    vtkMath::Cross(this->RefDirs[0], workVec, endVec);
    vtkMath::Normalize(endVec);
    vtkMath::MultiplyScalar(endVec, width/2.);

    // Get points to extrude up
    double topPt0[3], topPt1[3];
    vtkMath::Add(this->StartPt, endVec, topPt0);
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->StartPt, endVec, topPt1);

    // extrude down
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(topPt0, workVec, finalPts[1]);
    vtkMath::Add(topPt1, workVec, finalPts[2]);

    // extrude up
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(topPt0, workVec, finalPts[5]);
    vtkMath::Add(topPt1, workVec, finalPts[6]);

    // Get points to extrude up
    double endPt0[3], endPt1[3];
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->EndPt, endVec, endPt0);
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->EndPt, endVec, endPt1);

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(endPt0, workVec, finalPts[0]);
    vtkMath::Add(endPt1, workVec, finalPts[3]);

    // extrude up
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(endPt0, workVec, finalPts[4]);
    vtkMath::Add(endPt1, workVec, finalPts[7]);

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);
  }
  else if (cubeType == 1)
  {
    int numPoints = 10;
    double vecs[3][3];
    double endPts[2][3];
    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *diver = this->Children[this->DivergingChild];

    this->FormBifurcation(this->StartPt, this->EndPt,
                          diver->EndPt, diver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, endPts);

    // get vector towards top of cube
    double workVec[3];
    vtkMath::Cross(vecs[1], vecs[0], workVec);
    vtkMath::Normalize(workVec);

    // extrude up
    double finalPts[10][3];
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->EndPt, workVec, finalPts[4]);
    vtkMath::Add(endPts[0], workVec, finalPts[0]);
    vtkMath::Add(endPts[1], workVec, finalPts[3]);

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->EndPt, workVec, finalPts[9]);
    vtkMath::Add(endPts[0], workVec, finalPts[5]);
    vtkMath::Add(endPts[1], workVec, finalPts[8]);

    // Now get beginning points
    double endVec[3];
    vtkMath::Normalize(workVec);
    vtkMath::Cross(workVec, vecs[0], endVec);
    vtkMath::Normalize(endVec);
    vtkMath::MultiplyScalar(endVec, width/2.);

    // Get points to extrude up
    double topPt0[3], topPt1[3];
    vtkMath::Add(this->StartPt, endVec, topPt0);
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->StartPt, endVec, topPt1);

    // extrude down
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(topPt0, workVec, finalPts[6]);
    vtkMath::Add(topPt1, workVec, finalPts[7]);

    // extrude up
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(topPt0, workVec, finalPts[1]);
    vtkMath::Add(topPt1, workVec, finalPts[2]);

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);
  }
  else if (cubeType == 2)
  {
    int numPoints = 12;

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
    double workVec[3];
    vtkMath::Cross(diverVec, vecs[1], workVec);
    vtkMath::Normalize(workVec);

    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    double finalPts[12][3];
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->StartPt, workVec, finalPts[2]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[1]);
      vtkMath::Add(topPts[1], workVec, finalPts[3]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[3]);
      vtkMath::Add(topPts[1], workVec, finalPts[1]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->StartPt, workVec, finalPts[8]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[7]);
      vtkMath::Add(topPts[1], workVec, finalPts[9]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[9]);
      vtkMath::Add(topPts[1], workVec, finalPts[7]);
    }

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    double endPts[2][3];
    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, endPts);

    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], workVec);
    vtkMath::Normalize(workVec);

    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->EndPt, workVec, finalPts[5]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == BACK)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[0]);
      vtkMath::Add(endPts[1], workVec, finalPts[4]);
    }
    else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[4]);
      vtkMath::Add(endPts[1], workVec, finalPts[0]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->EndPt, workVec, finalPts[11]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == BACK)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[6]);
      vtkMath::Add(endPts[1], workVec, finalPts[10]);
    }
     else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[10]);
      vtkMath::Add(endPts[1], workVec, finalPts[6]);
    }

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);

  }
  else if (cubeType == 3)
  {
    int numPoints = 12;

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
    double workVec[3];
    vtkMath::Cross(diverVec, vecs[1], workVec);
    vtkMath::Normalize(workVec);

    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    double finalPts[12][3];
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->StartPt, workVec, finalPts[5]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[9]);
      vtkMath::Add(topPts[1], workVec, finalPts[1]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[1]);
      vtkMath::Add(topPts[1], workVec, finalPts[9]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->StartPt, workVec, finalPts[6]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[10]);
      vtkMath::Add(topPts[1], workVec, finalPts[2]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[2]);
      vtkMath::Add(topPts[1], workVec, finalPts[10]);
    }

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    double endPts[2][3];
    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          height/2., vecs, endPts);


    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], workVec);
    vtkMath::Normalize(workVec);

    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    vtkMath::MultiplyScalar(workVec, width/2.);
    vtkMath::Add(this->EndPt, workVec, finalPts[4]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == FRONT)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[0]);
      vtkMath::Add(endPts[1], workVec, finalPts[8]);
    }
    else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[8]);
      vtkMath::Add(endPts[1], workVec, finalPts[0]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->EndPt, workVec, finalPts[7]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == FRONT)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[3]);
      vtkMath::Add(endPts[1], workVec, finalPts[11]);
    }
     else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[11]);
      vtkMath::Add(endPts[1], workVec, finalPts[3]);
    }

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);

  }
  else if (cubeType == 4)
  {
    int numPoints = 12;

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
                         height/2., vecs, topPts);

    // get vector towards top of cube
    double diverVec[3];
    vtkMath::Subtract(diver->EndPt, diver->StartPt, diverVec);
    vtkMath::Normalize(diverVec);
    double workVec[3];
    vtkMath::Cross(diverVec, vecs[1], workVec);
    vtkMath::Normalize(workVec);

    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    double finalPts[12][3];
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->StartPt, workVec, finalPts[2]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[1]);
      vtkMath::Add(topPts[1], workVec, finalPts[3]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[3]);
      vtkMath::Add(topPts[1], workVec, finalPts[1]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->StartPt, workVec, finalPts[9]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[8]);
      vtkMath::Add(topPts[1], workVec, finalPts[10]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[10]);
      vtkMath::Add(topPts[1], workVec, finalPts[8]);
    }

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    double endPts[2][3];
    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          height/2., vecs, endPts);

    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], workVec);
    vtkMath::Normalize(workVec);

    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    vtkMath::MultiplyScalar(workVec, width/2.);
    vtkMath::Add(this->EndPt, workVec, finalPts[5]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == FRONT)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[0]);
      vtkMath::Add(endPts[1], workVec, finalPts[7]);
    }
    else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[7]);
      vtkMath::Add(endPts[1], workVec, finalPts[0]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->EndPt, workVec, finalPts[6]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == FRONT)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[4]);
      vtkMath::Add(endPts[1], workVec, finalPts[11]);
    }
     else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[11]);
      vtkMath::Add(endPts[1], workVec, finalPts[4]);
    }

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);
  }
  else if (cubeType == 5)
  {
    int numPoints = 12;

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
    double workVec[3];
    vtkMath::Cross(diverVec, vecs[1], workVec);
    vtkMath::Normalize(workVec);

    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    double finalPts[12][3];
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->StartPt, workVec, finalPts[5]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[8]);
      vtkMath::Add(topPts[1], workVec, finalPts[1]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[1]);
      vtkMath::Add(topPts[1], workVec, finalPts[8]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->StartPt, workVec, finalPts[6]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[9]);
      vtkMath::Add(topPts[1], workVec, finalPts[2]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[2]);
      vtkMath::Add(topPts[1], workVec, finalPts[9]);
    }

    vtkSVCenterlineGCell *align = this->Children[this->AligningChild];
    vtkSVCenterlineGCell *cDiver = this->Children[this->DivergingChild];

    double endPts[2][3];
    this->FormBifurcation(this->StartPt, this->EndPt,
                          cDiver->EndPt, cDiver->StartPt,
                          align->EndPt, align->StartPt,
                          this->EndPt,
                          width/2., vecs, endPts);

    // get vector towards top of cube
    vtkMath::Cross(vecs[1], vecs[0], workVec);
    vtkMath::Normalize(workVec);

    if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->EndPt, workVec, finalPts[4]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == BACK)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[0]);
      vtkMath::Add(endPts[1], workVec, finalPts[3]);
    }
    else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[3]);
      vtkMath::Add(endPts[1], workVec, finalPts[0]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->EndPt, workVec, finalPts[11]);
    if (cDiver->BranchDir == RIGHT || cDiver->BranchDir == BACK)
    {
      vtkMath::Add(endPts[0], workVec, finalPts[7]);
      vtkMath::Add(endPts[1], workVec, finalPts[10]);
    }
     else
    {
      vtkMath::Add(endPts[0], workVec, finalPts[10]);
      vtkMath::Add(endPts[1], workVec, finalPts[7]);
    }

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);

  }
  else if (cubeType == 6)
  {
    if (this->Parent->Children.size() == 2)
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
      double workVec[3];
      vtkMath::Cross(diverVec, vecs[1], workVec);
      vtkMath::Normalize(workVec);

      // extrude up
      double finalPts[10][3];
      vtkMath::MultiplyScalar(workVec, height/2.);
      if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(workVec, -1.0);

      vtkMath::Add(this->StartPt, workVec, finalPts[2]);
      if (this->BranchDir == RIGHT || this->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[1]);
        vtkMath::Add(topPts[1], workVec, finalPts[3]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[3]);
        vtkMath::Add(topPts[1], workVec, finalPts[1]);
      }

      // extrude down
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(this->StartPt, workVec, finalPts[7]);
      if (this->BranchDir == RIGHT || this->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[6]);
        vtkMath::Add(topPts[1], workVec, finalPts[8]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[8]);
        vtkMath::Add(topPts[1], workVec, finalPts[6]);
      }

      // Now get beginning points
      double endVec[3];
      vtkMath::Normalize(workVec);
      vtkMath::Cross(vecs[0], workVec, endVec);
      vtkMath::Normalize(endVec);
      vtkMath::MultiplyScalar(endVec, width/2.);

      // Get points to extrude up
      double endPt0[3], endPt1[3];
      vtkMath::Add(this->EndPt, endVec, endPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(this->EndPt, endVec, endPt1);

      // extrude down
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(endPt0, workVec, finalPts[5]);
      vtkMath::Add(endPt1, workVec, finalPts[9]);

      // extrude up
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(endPt0, workVec, finalPts[0]);
      vtkMath::Add(endPt1, workVec, finalPts[4]);

      for (int j=0; j<numPoints; j++)
        points->InsertNextPoint(finalPts[j]);
    }
    else
    {
      int numPoints = 10;

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
      double workVec[3];
      vtkMath::Cross(diverVec, vecs[1], workVec);
      vtkMath::Normalize(workVec);

      // extrude up
      double finalPts[10][3];
      vtkMath::MultiplyScalar(workVec, height/2.);
      if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(workVec, -1.0);

      vtkMath::Add(this->StartPt, workVec, finalPts[2]);
      if (this->BranchDir == RIGHT || this->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[1]);
        vtkMath::Add(topPts[1], workVec, finalPts[3]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[3]);
        vtkMath::Add(topPts[1], workVec, finalPts[1]);
      }

      // extrude down
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(this->StartPt, workVec, finalPts[7]);
      if (this->BranchDir == RIGHT || this->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[6]);
        vtkMath::Add(topPts[1], workVec, finalPts[8]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[8]);
        vtkMath::Add(topPts[1], workVec, finalPts[6]);
      }

      // Now get beginning points
      double endVec[3];
      vtkMath::Normalize(workVec);
      vtkMath::Cross(vecs[0], workVec, endVec);
      vtkMath::Normalize(endVec);
      vtkMath::MultiplyScalar(endVec, width/2.);

      // Get points to extrude up
      double endPt0[3], endPt1[3];
      vtkMath::Add(this->EndPt, endVec, endPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(this->EndPt, endVec, endPt1);

      // extrude down
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(endPt0, workVec, finalPts[5]);
      vtkMath::Add(endPt1, workVec, finalPts[9]);

      // extrude up
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(endPt0, workVec, finalPts[0]);
      vtkMath::Add(endPt1, workVec, finalPts[4]);

      for (int j=0; j<numPoints; j++)
        points->InsertNextPoint(finalPts[j]);
    }
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
    double workVec[3];
    vtkMath::Cross(diverVec, vecs[1], workVec);
    vtkMath::Normalize(workVec);

    if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      vtkMath::MultiplyScalar(workVec, -1.0);

    // extrude up
    double finalPts[10][3];
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->StartPt, workVec, finalPts[4]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[7]);
      vtkMath::Add(topPts[1], workVec, finalPts[1]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[1]);
      vtkMath::Add(topPts[1], workVec, finalPts[7]);
    }

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->StartPt, workVec, finalPts[5]);
    if (this->BranchDir == RIGHT || this->BranchDir == BACK)
    {
      vtkMath::Add(topPts[0], workVec, finalPts[8]);
      vtkMath::Add(topPts[1], workVec, finalPts[2]);
    }
    else
    {
      vtkMath::Add(topPts[0], workVec, finalPts[2]);
      vtkMath::Add(topPts[1], workVec, finalPts[8]);
    }

    // Now get beginning points
    double endVec[3];
    vtkMath::Normalize(workVec);
    vtkMath::Cross(vecs[0], workVec, endVec);
    vtkMath::Normalize(endVec);
    vtkMath::MultiplyScalar(endVec, width/2.);

    // Get points to extrude up
    double endPt0[3], endPt1[3];
    vtkMath::Add(this->EndPt, endVec, endPt0);
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->EndPt, endVec, endPt1);

    // extrude down
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(endPt0, workVec, finalPts[9]);
    vtkMath::Add(endPt1, workVec, finalPts[3]);

    // extrude up
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(endPt0, workVec, finalPts[6]);
    vtkMath::Add(endPt1, workVec, finalPts[0]);

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);
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
    double workVec[3];
    vtkMath::Cross(vecs[1], vecs[0], workVec);
    vtkMath::Normalize(workVec);

    // extrude up
    double finalPts[10][3];
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(this->EndPt, workVec, finalPts[4]);
    vtkMath::Add(endPts[0], workVec, finalPts[0]);
    vtkMath::Add(endPts[1], workVec, finalPts[3]);

    // extrude down
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(this->EndPt, workVec, finalPts[9]);
    vtkMath::Add(endPts[0], workVec, finalPts[5]);
    vtkMath::Add(endPts[1], workVec, finalPts[8]);

    // Now get beginning points
    double endVec[3];
    vtkMath::Normalize(workVec);
    vtkMath::Cross(workVec, vecs[0], endVec);
    vtkMath::Normalize(endVec);
    vtkMath::MultiplyScalar(endVec, width/2.);

    // Get points to extrude up
    double topPt0[3], topPt1[3];
    vtkMath::Add(this->StartPt, endVec, topPt0);
    vtkMath::MultiplyScalar(endVec, -1.0);
    vtkMath::Add(this->StartPt, endVec, topPt1);

    // extrude down
    vtkMath::MultiplyScalar(workVec, height/2.);
    vtkMath::Add(topPt0, workVec, finalPts[6]);
    vtkMath::Add(topPt1, workVec, finalPts[7]);

    // extrude up
    vtkMath::MultiplyScalar(workVec, -1.0);
    vtkMath::Add(topPt0, workVec, finalPts[1]);
    vtkMath::Add(topPt1, workVec, finalPts[2]);

    for (int j=0; j<numPoints; j++)
      points->InsertNextPoint(finalPts[j]);
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

