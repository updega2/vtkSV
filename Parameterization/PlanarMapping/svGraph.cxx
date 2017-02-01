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

#include "svGraph.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkPoints.h"
#include "vtkPolyDataSliceAndDiceFilter.h"
#include "vtkSmartPointer.h"
#include "vtkThreshold.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

svGraph::svGraph()
{
  this->NumberOfCells = 0;
  this->Root = NULL;

  this->Lines = vtkPolyData::New();
}

svGraph::svGraph(int rootId,
                 vtkPolyData *linesPd,
                 std::string groupIdsArrayName,
                 std::multimap<int, int> criticalPointMap,
                 int directionTable[6][4])
{
  this->NumberOfCells = 0;
  double startPt[3]; startPt[0] = 0.0; startPt[1] = 0.0; startPt[2] = 0.0;
  double endPt[3]; endPt[0] = 0.0; endPt[1] = 0.0; endPt[2] = -1.0;
  this->Root = this->NewCell(rootId, vtkPolyDataSliceAndDiceFilter::DOWN, startPt, endPt);

  this->Lines = vtkPolyData::New();
  this->Lines->DeepCopy(linesPd);
  this->GroupIdsArrayName = groupIdsArrayName;
  this->CriticalPointMap = criticalPointMap;
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<4; j++)
    {
      this->DirectionTable[i][j] = directionTable[i][j];
    }
  }
}

svGraph::~svGraph()
{
  if (this->Root != NULL)
  {
    delete this->Root;
    this->Root = NULL;
  }
  if (this->Lines != NULL)
  {
    this->Lines->Delete();
    this->Lines = NULL;
  }
}

int svGraph::Recurse(svGCell *rootGCell, int(*function)(svGCell *currentGCell,
                   void *arg0, void *arg1, void *arg2),
                   void *rec_arg0, void *rec_arg1, void *rec_arg2)
{
  function(rootGCell, rec_arg0, rec_arg1, rec_arg2);
  if (rootGCell->Children[0] != NULL)
  {
    svGraph::Recurse(rootGCell->Children[0], function, rec_arg0, rec_arg1, rec_arg2);
  }
  if (rootGCell->Children[1] != NULL)
  {
    svGraph::Recurse(rootGCell->Children[1], function, rec_arg0, rec_arg1, rec_arg2);
  }
  return 1;
}

int svGraph::PrintGraph()
{
  svGraph::Recurse(this->Root, svGraph::PrintGCell, NULL, NULL, NULL);
  return 1;
}

int svGraph::PrintGCell(svGCell *gCell, void *arg0, void *arg1, void *arg2)
{
  fprintf(stdout, "GCell ID: %d\n", gCell->Id);
  fprintf(stdout, "GCell Group ID: %d\n", gCell->GroupId);
  fprintf(stdout, "Direction: %d\n", gCell->Dir);
  if(gCell->Parent != NULL)
    fprintf(stdout, "Parent: %d\n", gCell->Parent->Id);
  else
    fprintf(stdout, "Parent is NULL\n");
  if(gCell->Children[0] != NULL)
    fprintf(stdout, "Child 0: %d\n", gCell->Children[0]->Id);
  else
    fprintf(stdout, "Child 0: NULL\n");
  if(gCell->Children[1] != NULL)
    fprintf(stdout, "Child 1: %d\n", gCell->Children[1]->Id);
  else
    fprintf(stdout, "Child 1: NULL\n");
  fprintf(stdout, "Start Point: %.4f %.4f %.4f\n", gCell->StartPt[0], gCell->StartPt[1], gCell->StartPt[2]);
  fprintf(stdout, "End Point: %.4f %.4f %.4f\n", gCell->EndPt[0], gCell->EndPt[1], gCell->EndPt[2]);
  return 1;
}

int svGraph::GetGraphPolyData(vtkPolyData *pd)
{
  vtkNew(vtkPoints, newPoints);
  vtkNew(vtkCellArray, newCells);
  vtkNew(vtkIntArray, groupIds);
  groupIds->SetName(this->GroupIdsArrayName.c_str());
  svGraph::Recurse(this->Root, svGraph::InsertGCellPoints, newPoints, newCells, groupIds);

  pd->SetPoints(newPoints);
  pd->SetLines(newCells);
  pd->GetCellData()->AddArray(groupIds);

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(pd);
  cleaner->Update();

  pd->DeepCopy(cleaner->GetOutput());
  pd->BuildLinks();
  return 1;
}

int svGraph::InsertGCellPoints(svGCell *gCell, void *arg0, void *arg1, void *arg2)
{
  vtkPoints *points =
    reinterpret_cast<vtkPoints*>(arg0);
  int numPoints = points->GetNumberOfPoints();
  points->InsertNextPoint(gCell->StartPt);
  points->InsertNextPoint(gCell->EndPt);

  vtkCellArray *cells =
    reinterpret_cast<vtkCellArray*>(arg1);
  vtkNew(vtkIdList, newIds);
  newIds->InsertNextId(numPoints); newIds->InsertNextId(numPoints+1);
  cells->InsertNextCell(newIds);

  vtkIntArray *groupIds =
    reinterpret_cast<vtkIntArray*>(arg2);
  groupIds->InsertNextValue(gCell->GroupId);

  return 1;
}

int svGraph::BuildGraph()
{
  fprintf(stdout,"Building graph!!!\n");

  std::list<int> keyVals;
  vtkPolyDataSliceAndDiceFilter::GetValuesFromMap(this->CriticalPointMap, this->Root->GroupId, keyVals);
  if (keyVals.size() != 0)
  {
    std::list<int> children;
    vtkPolyDataSliceAndDiceFilter::GetUniqueNeighbors(this->CriticalPointMap, this->Root->GroupId, keyVals, children);
    if (children.size() != 0)
    {
      std::list<int>::iterator childit = children.begin();
      this->Root->Children[0] = this->NewCell(*childit, this->Root); ++childit;
      this->Root->Children[1] = this->NewCell(*childit, this->Root); ++childit;

      this->ComputeFirstCellVector(this->Root);
      this->GetNewBranchDirections(this->Root);
      this->GrowGraph(this->Root->Children[0]);
      this->GrowGraph(this->Root->Children[1]);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svGraph::GrowGraph(svGCell *parent)
{
  std::list<int> keyVals;
  vtkPolyDataSliceAndDiceFilter::GetValuesFromMap(this->CriticalPointMap, parent->GroupId, keyVals);
  if (keyVals.size() != 0)
  {
    std::list<int> children;
    vtkPolyDataSliceAndDiceFilter::GetUniqueNeighbors(this->CriticalPointMap, parent->GroupId, keyVals, children);
    if (children.size() > 2)
    {
      //fprintf(stdout,"Number of children!: %lu\n", children.size());
      std::list<int>::iterator childit = children.begin();
      int count = 0;
      for (int i=0; childit != children.end(); ++childit)
      {
        if (*childit != parent->Parent->Children[0]->GroupId && *childit != parent->Parent->Children[1]->GroupId && *childit != parent->Parent->GroupId)
        {
          parent->Children[count] = this->NewCell(*childit, parent);
          count++;
        }
      }
      //fprintf(stdout,"What is the child %d\n", parent->Children[0]->GroupId);
      //fprintf(stdout,"What is the child %d\n", parent->Children[1]->GroupId);
      this->GetNewBranchDirections(parent);

      this->GrowGraph(parent->Children[0]);
      this->GrowGraph(parent->Children[1]);
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svGraph::ComputeFirstCellVector(svGCell *parent)
{
  vtkNew(vtkThreshold, thresholder);
  thresholder->SetInputData(this->Lines);
  thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
  thresholder->ThresholdBetween(parent->GroupId, parent->GroupId);
  thresholder->Update();
  int numPts = thresholder->GetOutput()->GetNumberOfPoints();

  double endPt0[3], endPt1[3];
  thresholder->GetOutput()->GetPoint(numPts-2, endPt0);
  thresholder->GetOutput()->GetPoint(numPts-1, endPt1);

  double startPts[2][3];
  for (int i=0; i<2; i++)
  {
    thresholder->SetInputData(this->Lines);
    thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
    thresholder->ThresholdBetween(parent->Children[i]->GroupId, parent->Children[i]->GroupId);
    thresholder->Update();

    thresholder->GetOutput()->GetPoint(0, startPts[i]);
  }

  double vec0[3], vec1[3], rootVec[3];
  vtkMath::Subtract(endPt1, endPt0, rootVec);
  vtkMath::Normalize(rootVec);
  vtkMath::Subtract(startPts[0], endPt1, vec0);
  vtkMath::Normalize(vec0);
  vtkMath::Subtract(startPts[1], endPt1, vec1);
  vtkMath::Normalize(vec1);

  // Get angle between vectors
  double angleVec0[3], angleVec1[3];
  vtkMath::Cross(vec0, rootVec, angleVec0);
  double ang0 = atan2(vtkMath::Norm(angleVec0), vtkMath::Dot(vec0, rootVec));
  vtkMath::Cross(vec1, rootVec, angleVec1);
  double ang1 = atan2(vtkMath::Norm(angleVec1), vtkMath::Dot(vec1, rootVec));

  if (ang0 < ang1)
  {
    vtkMath::Cross(rootVec, vec1, this->FirstVec);
  }
  else
  {
    vtkMath::Cross(rootVec, vec0, this->FirstVec);
  }
  vtkMath::Normalize(this->FirstVec);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svGraph::GetDirectionVector(const int dir, double dirVector[3])
{
  if (dir == vtkPolyDataSliceAndDiceFilter::RIGHT)
  {
    dirVector[0] = 1.0; dirVector[1] = 0.0; dirVector[2] = 0.0;
  }
  if (dir == vtkPolyDataSliceAndDiceFilter::LEFT)
  {
    dirVector[0] = -1.0; dirVector[1] = 0.0; dirVector[2] = 0.0;
  }
  if (dir == vtkPolyDataSliceAndDiceFilter::FRONT)
  {
    dirVector[0] = 0.0; dirVector[1] = 1.0; dirVector[2] = 0.0;
  }
  if (dir == vtkPolyDataSliceAndDiceFilter::BACK)
  {
    dirVector[0] = 0.0; dirVector[1] = -1.0; dirVector[2] = 0.0;
  }
  if (dir == vtkPolyDataSliceAndDiceFilter::UP)
  {
    dirVector[0] = 0.0; dirVector[1] = 0.0; dirVector[2] = 1.0;
  }
  if (dir == vtkPolyDataSliceAndDiceFilter::DOWN)
  {
    dirVector[0] = 0.0; dirVector[1] = 0.0; dirVector[2] = -1.0;
  }
  return 1;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svGraph::GetNewBranchDirections(svGCell *parent)
{
  //fprintf(stdout,"Child %d of parent %d, dir: %d\n", parent->Children[0]->GroupId, parent->GroupId, parent->Children[0]->Dir);
  //fprintf(stdout,"Child %d of parent %d, dir: %d\n", parent->Children[1]->GroupId, parent->GroupId, parent->Children[1]->Dir);
  vtkNew(vtkThreshold, thresholder);
  thresholder->SetInputData(this->Lines);
  thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
  thresholder->ThresholdBetween(parent->GroupId, parent->GroupId);
  thresholder->Update();
  int numPts = thresholder->GetOutput()->GetNumberOfPoints();

  double endPt0[3], endPt1[3];
  thresholder->GetOutput()->GetPoint(numPts-2, endPt0);
  thresholder->GetOutput()->GetPoint(numPts-1, endPt1);

  double startPts[2][3];
  for (int i=0; i<2; i++)
  {
    thresholder->SetInputData(this->Lines);
    thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
    thresholder->ThresholdBetween(parent->Children[i]->GroupId, parent->Children[i]->GroupId);
    thresholder->Update();

    thresholder->GetOutput()->GetPoint(0, startPts[i]);
  }

  double vec0[3], vec1[3], vec2[3];
  vtkMath::Subtract(endPt1, endPt0, vec0);
  vtkMath::Normalize(vec0);
  vtkMath::Subtract(startPts[0], endPt1, vec1);
  vtkMath::Normalize(vec1);
  vtkMath::Subtract(startPts[1], endPt1, vec2);
  vtkMath::Normalize(vec2);

  // Get angle between vectors
  double angleVec0[3], angleVec1[3];
  vtkMath::Cross(vec1, vec0, angleVec0);
  double ang0 = atan2(vtkMath::Norm(angleVec0), vtkMath::Dot(vec1, vec0));
  vtkMath::Cross(vec2, vec0, angleVec1);
  double ang1 = atan2(vtkMath::Norm(angleVec1), vtkMath::Dot(vec2, vec0));

  double vec3[3];
  int direction0 = parent->Dir;
  if (ang0 < ang1)
  {
    vtkMath::Cross(vec0, vec2, vec3);
    parent->Children[0]->Dir = direction0;
  }
  else
  {
    vtkMath::Cross(vec0, vec1, vec3);
    parent->Children[1]->Dir = direction0;
  }

  //Check to see if right or left. Perpendicular vector from new branch
  //will be either 0 degrees (RIGHT) or 180 degrees (LEFT) with respect
  //to FirstBranchVec! Very important, this may not be the direction it
  //goes on the skeleton. We need to add that direction to parent direction
  //Need to put in a tolerance as well.

  //If not RIGHT or LEFT, then it is either BACK (90 degrees)
  //or FRONT (270 degrees), again with respect to FirstBranchVec

  vtkMath::Normalize(vec3);
  //fprintf(stdout,"First Branch: %.4f %.4f %.4f\n", this->FirstVec[0], this->FirstVec[1], this->FirstVec[2]);
  //fprintf(stdout,"This Branch: %.4f %.4f %.4f\n", vec3[0], vec3[1], vec3[2]);
  double angleVec2[3];
  vtkMath::Cross(vec3, this->FirstVec, angleVec2);
  double ang2 = atan2(vtkMath::Norm(angleVec2), vtkMath::Dot(vec3, this->FirstVec));

  //fprintf(stdout,"Angle 0: %4f\n", 180*ang0/M_PI);
  //fprintf(stdout,"Angle 1: %4f\n", 180*ang1/M_PI);
  //fprintf(stdout,"Angle 2: %4f\n", 180*ang2/M_PI);
  int direction1;
  if (ang2 >= 7.0*M_PI/4.0 || ang2 < M_PI/4.0)
  {
    direction1 = this->DirectionTable[direction0][0];
  }
  else if (ang2 >= M_PI/4.0 && ang2 < 3.0*M_PI/4.0)
  {
    direction1 = this->DirectionTable[direction0][1];
  }
  else if (ang2 >= 3.0*M_PI/4.0 && ang2 < 5.0*M_PI/4.0)
  {
    direction1 = this->DirectionTable[direction0][2];
  }
  else if (ang2 >= 5.0*M_PI/4.0 && ang2 < 7.0*M_PI/4.0)
  {
    direction1 = this->DirectionTable[direction0][3];
  }
  if (ang0 < ang1)
  {
    parent->Children[1]->Dir = direction1;
  }
  else
  {
    parent->Children[0]->Dir = direction1;
  }
  double dirVector0[3], dirVector1[3];
  this->GetDirectionVector(parent->Children[0]->Dir , dirVector0);
  this->GetDirectionVector(parent->Children[1]->Dir , dirVector1);

  for (int i=0; i<3; i++)
  {
    parent->Children[0]->EndPt[i] = parent->Children[0]->StartPt[i] + dirVector0[i];
    parent->Children[1]->EndPt[i] = parent->Children[1]->StartPt[i] + dirVector1[i];
  }

  return 1;
}

svGCell* svGraph::NewCell(int a_GroupId, svGCell *a_Parent)
{
  svGCell *newCell = new svGCell;
  newCell->Parent      = a_Parent;
  newCell->Children[0] = NULL;
  newCell->Children[1] = NULL;
  newCell->Id       = this->NumberOfCells++;
  newCell->GroupId  = a_GroupId;
  newCell->Dir = -1;
  for (int i=0; i<3; i++)
  {
    if (a_Parent != NULL)
    {
      newCell->StartPt[i] = a_Parent->EndPt[i];
    }
    else
    {
      newCell->StartPt[i] = -1;
    }
    newCell->EndPt[i]   = -1;
  }

  return newCell;
}

svGCell* svGraph::NewCell(int a_GroupId, int a_Dir, double a_StartPt[3], double a_EndPt[3])
{
  svGCell *newCell = new svGCell;
  newCell->Parent      = NULL;
  newCell->Children[0] = NULL;
  newCell->Children[1] = NULL;
  newCell->Id       = this->NumberOfCells++;
  newCell->GroupId  = a_GroupId;
  newCell->Dir      = a_Dir;
  for (int i=0; i<3; i++)
  {
    newCell->StartPt[i] = a_StartPt[i];
    newCell->EndPt[i]   = a_EndPt[i];
  }

  return newCell;
}

svGCell* svGraph::GetCell(const int findId)
{
  if (findId > this->NumberOfCells)
  {
    fprintf(stdout,"Id is larger than number of cells\n");
    return NULL;
  }
  return this->LookUp(this->Root, findId);
}

svGCell* svGraph::LookUp(svGCell *lookCell, const int findId)
{
  if (lookCell == NULL)
  {
    return lookCell;
  }
  if (lookCell->Id == findId)
  {
    return lookCell;
  }
  else
  {
    if (lookCell->Children[0] != NULL)
    {
      svGCell* foundCell = this->LookUp(lookCell->Children[0], findId);
      if (foundCell == NULL)
      {
        foundCell = this->LookUp(lookCell->Children[1], findId);
      }
      return foundCell;
    }
  }

  return NULL;
}
