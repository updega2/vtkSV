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

#include "svCenterlineGraph.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSmartPointer.h"
#include "vtkThreshold.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkUnstructuredGrid.h"

svCenterlineGraph::svCenterlineGraph()
{
  this->NumberOfCells = 0;
  this->NumberOfNodes = 1; // The root
  this->Root = NULL;

  this->Lines = vtkPolyData::New();
  this->Graph = vtkMutableDirectedGraph::New();
}

// ----------------------
// Constructor
// ----------------------
svCenterlineGraph::svCenterlineGraph(int rootId,
                 vtkPolyData *linesPd,
                 std::string groupIdsArrayName)
{
  this->NumberOfCells = 0;
  this->NumberOfNodes = 1; // The root
  double startPt[3]; startPt[0] = 0.0; startPt[1] = 0.0; startPt[2] = 0.0;
  double endPt[3]; endPt[0] = 0.0; endPt[1] = 0.0; endPt[2] = -1.0;
  this->Root = this->NewCell(rootId, DOWN, startPt, endPt);

  this->Lines = vtkPolyData::New();
  this->Lines->DeepCopy(linesPd);
  this->Lines->BuildLinks();
  this->GroupIdsArrayName = groupIdsArrayName;
}

// ----------------------
// Destructor
// ----------------------
svCenterlineGraph::~svCenterlineGraph()
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
  if (this->Graph != NULL)
  {
    this->Graph->Delete();
    this->Graph = NULL;
  }
}

// ----------------------
// Recurse
// ----------------------
int svCenterlineGraph::Recurse(svCenterlineGCell *rootGCell, int(*function)(svCenterlineGCell *currentGCell,
                   void *arg0, void *arg1, void *arg2),
                   void *rec_arg0, void *rec_arg1, void *rec_arg2)
{
  function(rootGCell, rec_arg0, rec_arg1, rec_arg2);
  if (rootGCell->Children.size() != 0)
  {
    for (int i=0; i<rootGCell->Children.size(); i++)
      svCenterlineGraph::Recurse(rootGCell->Children[i], function, rec_arg0, rec_arg1, rec_arg2);
  }
  return SV_OK;
}

// ----------------------
// PrintGraph
// ----------------------
int svCenterlineGraph::PrintGraph()
{
  svCenterlineGraph::Recurse(this->Root, svCenterlineGraph::PrintGCell, NULL, NULL, NULL);
  return SV_OK;
}

// ----------------------
// PrintGCell
// ----------------------
int svCenterlineGraph::PrintGCell(svCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2)
{
  fprintf(stdout, "GCell ID: %d\n", gCell->Id);
  fprintf(stdout, "GCell Group ID: %d\n", gCell->GroupId);
  fprintf(stdout, "Direction: %d\n", gCell->Dir);
  if(gCell->Parent != NULL)
    fprintf(stdout, "Parent: %d\n", gCell->Parent->Id);
  else
    fprintf(stdout, "Parent is NULL\n");
  if(gCell->Children.size() != 0)
  {
    for (int i=0; i<gCell->Children.size(); i++)
      fprintf(stdout, "Child %d, Id: %d, GroupId: %d\n", i, gCell->Children[i]->Id, gCell->Children[i]->GroupId);
  }
  fprintf(stdout, "Start Point: %.4f %.4f %.4f\n", gCell->StartPt[0], gCell->StartPt[1], gCell->StartPt[2]);
  fprintf(stdout, "End Point: %.4f %.4f %.4f\n", gCell->EndPt[0], gCell->EndPt[1], gCell->EndPt[2]);
  return SV_OK;
}

// ----------------------
// GetGraphPolyData
// ----------------------
int svCenterlineGraph::GetGraphPolyData(vtkPolyData *pd)
{
  vtkNew(vtkPoints, newPoints);
  vtkNew(vtkCellArray, newCells);
  vtkNew(vtkIntArray, groupIds);
  groupIds->SetName(this->GroupIdsArrayName.c_str());
  svCenterlineGraph::Recurse(this->Root, svCenterlineGraph::InsertGCellPoints, newPoints, newCells, groupIds);

  pd->SetPoints(newPoints);
  pd->SetLines(newCells);
  pd->GetCellData()->AddArray(groupIds);

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(pd);
  cleaner->Update();

  pd->DeepCopy(cleaner->GetOutput());
  pd->BuildLinks();
  return SV_OK;
}

// ----------------------
// InsertGCellPoints
// ----------------------
int svCenterlineGraph::InsertGCellPoints(svCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2)
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

  return SV_OK;
}

// ----------------------
// BuildGraph
// ----------------------
int svCenterlineGraph::BuildGraph()
{
  fprintf(stdout,"Building graph!!!\n");

  vtkNew(vtkIdList, connectingGroups);
  this->GetConnectingLineGroups(this->Root->GroupId, connectingGroups);

  if (connectingGroups->GetNumberOfIds() == 0)
  {
    // Just only this guy
    fprintf(stdout,"Only one centerlines\n");
    return SV_OK;
  }


  for (int i=0; i<connectingGroups->GetNumberOfIds(); i++)
    this->Root->Children.push_back(this->NewCell(connectingGroups->GetId(i), this->Root));

  this->NumberOfNodes += this->Root->Children.size()+1;

  this->ComputeReferenceVectors(this->Root);
  this->GetNewBranchDirections(this->Root);

  for (int i=0; i<this->Root->Children.size(); i++)
    this->GrowGraph(this->Root->Children[i]);

  this->RefineGraphWithLocalCoordinates();
  svCenterlineGraph::Recurse(this->Root, svCenterlineGraph::UpdateCellDirection,
                   NULL, NULL, NULL);

  return SV_OK;
}

// ----------------------
// GrowGraph
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svCenterlineGraph::GrowGraph(svCenterlineGCell *parent)
{
  vtkNew(vtkIdList, connectingGroups);
  this->GetConnectingLineGroups(parent->GroupId, connectingGroups);

  int count = 0;
  for (int i=0; i<connectingGroups->GetNumberOfIds(); i++)
  {
    int groupId = connectingGroups->GetId(i);
    int isChild=1;
    for (int j=0; j<parent->Parent->Children.size(); j++)
    {
      if (groupId == parent->Parent->Children[j]->GroupId)
        isChild = 0;
    }
    if (groupId == parent->Parent->GroupId)
      isChild = 0;
    if (isChild)
    {
      parent->Children.push_back(this->NewCell(groupId, parent));
      count++;
    }
  }
  if (count > 0)
  {
    this->NumberOfNodes += parent->Children.size()+1;
    //fprintf(stdout,"What is the child %d\n", parent->Children[0]->GroupId);
    //fprintf(stdout,"What is the child %d\n", parent->Children[1]->GroupId);
    this->GetNewBranchDirections(parent);

    for (int i=0; i<parent->Children.size(); i++)
      this->GrowGraph(parent->Children[i]);
  }
  else
  {
    fprintf(stdout,"Terminating branch\n");
    return SV_OK;
  }

  return SV_OK;
}

// ----------------------
// GetConnectingLineGroups
// ----------------------
int svCenterlineGraph::GetConnectingLineGroups(const int groupId, vtkIdList *connectingGroups)
{
  // Get first cell points
  vtkIdType npts, *pts;
  int cellId = this->Lines->GetCellData()->GetArray(this->GroupIdsArrayName.c_str())->LookupValue(groupId);
  this->Lines->GetCellPoints(cellId, npts, pts);

  // Get connecting
  vtkNew(vtkIdList, pointCells0);
  vtkNew(vtkIdList, pointCellsN);
  this->Lines->GetPointCells(pts[0], pointCells0);
  this->Lines->GetPointCells(pts[npts-1], pointCellsN);

  for (int i=0; i<pointCells0->GetNumberOfIds(); i++)
  {
    if (pointCells0->GetId(i) != cellId)
    {
      connectingGroups->InsertNextId(this->Lines->GetCellData()->GetArray(
        this->GroupIdsArrayName.c_str())->GetTuple1(pointCells0->GetId(i)));
    }
  }
  for (int i=0; i<pointCellsN->GetNumberOfIds(); i++)
  {
    if (pointCellsN->GetId(i) != cellId)
    {
      connectingGroups->InsertNextId(this->Lines->GetCellData()->GetArray(
        this->GroupIdsArrayName.c_str())->GetTuple1(pointCellsN->GetId(i)));
    }
  }

  return SV_OK;
}

// ----------------------
// ComputeReferenceVectors
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svCenterlineGraph::ComputeReferenceVectors(svCenterlineGCell *parent)
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
  fprintf(stdout,"End pt 0 is: %.4f %.4f %.4f\n", endPt0[0], endPt0[1], endPt0[2]);
  fprintf(stdout,"End pt 1 is: %.4f %.4f %.4f\n", endPt1[0], endPt1[1], endPt1[2]);
  vtkMath::Subtract(endPt0, endPt1, this->ReferenceVecs[0]);
  vtkMath::Normalize(this->ReferenceVecs[0]);

  int numChildren = parent->Children.size();
  std::vector<double> angs(numChildren);
  for (int i=0; i<numChildren; i++)
  {
    thresholder->SetInputData(this->Lines);
    thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
    thresholder->ThresholdBetween(parent->Children[i]->GroupId, parent->Children[i]->GroupId);
    thresholder->Update();

    double startPts[3], secondPts[3];
    thresholder->GetOutput()->GetPoint(0, startPts);
    thresholder->GetOutput()->GetPoint(1, secondPts);

    vtkMath::Subtract(secondPts, startPts, parent->Children[i]->StartVec);
    vtkMath::Normalize(parent->Children[i]->StartVec);

    // Get angle between vectors
    double angleVec[3];
    vtkMath::Cross(parent->Children[i]->StartVec, this->ReferenceVecs[0], angleVec);
    angs[i] = atan2(vtkMath::Norm(angleVec), vtkMath::Dot(parent->Children[i]->StartVec, this->ReferenceVecs[0]));
  }

  double minAngle = VTK_SV_LARGE_DOUBLE;
  double maxAngle = -1.0*VTK_SV_LARGE_DOUBLE;
  int minChild = 0;
  int maxChild = 0;
  for (int i=0; i<numChildren; i++)
  {
    if (angs[i] < minAngle)
    {
      minAngle = angs[i];
      minChild = i;
    }
    if (angs[i] > maxAngle)
    {
      maxAngle = angs[i];
      maxChild = i;
    }
    fprintf(stdout,"Vec %d: %.4f %.4f %.4f\n", i, parent->Children[i]->StartVec[0], parent->Children[i]->StartVec[1], parent->Children[i]->StartVec[2]);
    fprintf(stdout,"Angle %d: %4f\n", i, 180*angs[i]/M_PI);
  }
  vtkMath::Cross(this->ReferenceVecs[0], parent->Children[minChild]->StartVec, this->ReferenceVecs[2]);
  vtkMath::Normalize(this->ReferenceVecs[2]);

  parent->DivergingChild = minChild;
  parent->AligningChild  = maxChild;

  parent->Children[parent->AligningChild]->IsAlign = 1;
  parent->Children[parent->DivergingChild]->IsAlign = 0;

  vtkMath::Cross(this->ReferenceVecs[2], this->ReferenceVecs[0], this->ReferenceVecs[1]);
  vtkMath::Normalize(this->ReferenceVecs[1]);
  fprintf(stdout,"Root vec: %.4f %.4f %.4f\n", this->ReferenceVecs[0][0], this->ReferenceVecs[0][1], this->ReferenceVecs[0][2]);

  for (int i=0 ; i<3; i++)
  {
    for (int j=0; j<3; j++)
      parent->RefDirs[i][j] = this->ReferenceVecs[i][j];
  }

  return SV_OK;
}

// ----------------------
// GetDirectionVector
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svCenterlineGraph::GetDirectionVector(const int dir, double dirVector[3])
{
  if (dir == RIGHT)
  {
    dirVector[0] = 1.0; dirVector[1] = 0.0; dirVector[2] = 0.0;
  }
  if (dir == LEFT)
  {
    dirVector[0] = -1.0; dirVector[1] = 0.0; dirVector[2] = 0.0;
  }
  if (dir == FRONT)
  {
    dirVector[0] = 0.0; dirVector[1] = -1.0; dirVector[2] = 0.0;
  }
  if (dir == BACK)
  {
    dirVector[0] = 0.0; dirVector[1] = 1.0; dirVector[2] = 0.0;
  }
  if (dir == UP)
  {
    dirVector[0] = 0.0; dirVector[1] = 0.0; dirVector[2] = 1.0;
  }
  if (dir == DOWN)
  {
    dirVector[0] = 0.0; dirVector[1] = 0.0; dirVector[2] = -1.0;
  }
  return SV_OK;
}


// ----------------------
// GetNewBranchDirections
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svCenterlineGraph::GetNewBranchDirections(svCenterlineGCell *parent)
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

  double vec0[3];
  vtkMath::Subtract(endPt0, endPt1, vec0);
  vtkMath::Normalize(vec0);

  double refDirs[3][3];
  for (int i=0; i<3; i++)
    refDirs[0][i] = vec0[i];

  int maxDir;
  double maxDot = -0.1;
  for (int i=0; i<3; i++)
  {
    double compare = fabs(vtkMath::Dot(this->ReferenceVecs[i], vec0));
    fprintf(stdout,"Dot with Ref %d: %.4f\n", i, vtkMath::Dot(this->ReferenceVecs[i], vec0));
    if (compare > maxDot)
    {
      maxDot = compare;
      maxDir = i;
    }
  }
  fprintf(stdout,"Direction aligns most with: %d\n", maxDir);

  int numChildren = parent->Children.size();
  std::vector<double> angs(numChildren);
  for (int i=0; i<numChildren; i++)
  {
    thresholder->SetInputData(this->Lines);
    thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
    thresholder->ThresholdBetween(parent->Children[i]->GroupId, parent->Children[i]->GroupId);
    thresholder->Update();

    double startPts[3], secondPts[3];
    thresholder->GetOutput()->GetPoint(0, startPts);
    thresholder->GetOutput()->GetPoint(1, secondPts);

    vtkMath::Subtract(secondPts, startPts, parent->Children[i]->StartVec);
    vtkMath::Normalize(parent->Children[i]->StartVec);

    // Get angle between vectors
    double angleVec[3];
    vtkMath::Cross(parent->Children[i]->StartVec, vec0, angleVec);
    angs[i] = atan2(vtkMath::Norm(angleVec), vtkMath::Dot(parent->Children[i]->StartVec, vec0));
  }

  double minAngle = VTK_SV_LARGE_DOUBLE;
  double maxAngle = -1.0*VTK_SV_LARGE_DOUBLE;
  int minChild = 0;
  int maxChild = 0;
  for (int i=0; i<numChildren; i++)
  {
    if (angs[i] < minAngle)
    {
      minAngle = angs[i];
      minChild = i;
    }
    if (angs[i] > maxAngle)
    {
      maxAngle = angs[i];
      maxChild = i;
    }
    fprintf(stdout,"Vec %d: %.4f %.4f %.4f\n", i, parent->Children[i]->StartVec[0], parent->Children[i]->StartVec[1], parent->Children[i]->StartVec[2]);
    fprintf(stdout,"Angle %d: %4f\n", i, 180*angs[i]/M_PI);
  }

  double vec3[3];
  vtkMath::Cross(vec0, parent->Children[minChild]->StartVec, vec3);
  vtkMath::Normalize(vec3);

  parent->DivergingChild = minChild;
  parent->AligningChild  = maxChild;

  parent->Children[parent->AligningChild]->IsAlign  = 1;
  parent->Children[parent->DivergingChild]->IsAlign = 0;

  vtkMath::Normalize(vec3);
  vtkMath::Cross(vec3, vec0, refDirs[1]);
  vtkMath::Normalize(refDirs[1]);
  for (int i=0; i<3; i++)
    refDirs[2][i] = vec3[i];

  // Now children just in case terminating children
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      for (int k=0; k<parent->Children.size(); k++)
        parent->Children[k]->RefDirs[i][j] = refDirs[i][j];
    }
  }

  //Check to see orientation with respect to reference vectors. Angle that
  //it makes with reference vector determines the input angle to a
  //lookup table that tells us the absolute direction that the new branch
  //node should go

  double checkVec[3], dotVec[3];
  for (int i=0; i<3; i++)
  {
    checkVec[i] = this->ReferenceVecs[(maxDir+2)%3][i];
    dotVec[i]   = this->ReferenceVecs[(maxDir+1)%3][i];
  }

  fprintf(stdout,"This Branch: %.4f %.4f %.4f\n", vec3[0], vec3[1], vec3[2]);
  double checkAngVec[3];
  vtkMath::Cross(vec3, checkVec, checkAngVec);
  double ang2 = atan2(vtkMath::Norm(checkAngVec), vtkMath::Dot(vec3, checkVec));
  fprintf(stdout,"Dot between vec3 and ref is %.4f\n", vtkMath::Dot(vec3, dotVec));

  if (vtkMath::Dot(vec3, dotVec) < 0.0)
    ang2 = 2*M_PI - ang2;
  parent->Children[parent->DivergingChild]->RefAngle = ang2;

  fprintf(stdout,"Angle check dir: %.4f %.4f %.4f\n", checkVec[0], checkVec[1], checkVec[2]);
  fprintf(stdout,"Dot check dir: %.4f %.4f %.4f\n", dotVec[0], dotVec[1], dotVec[2]);
  fprintf(stdout,"This vec: %.4f %.4f %.4f\n", vec3[0], vec3[1], vec3[2]);
  fprintf(stdout,"Angle is: %.4f\n", 180*ang2/M_PI);
  fprintf(stdout,"Dot is: %.4f\n", vtkMath::Dot(vec3, dotVec));

  return SV_OK;
}

// ----------------------
// UpdateCellDirection
// ----------------------
int svCenterlineGraph::UpdateCellDirection(svCenterlineGCell *gCell, void *arg0,
                                 void *arg1, void *arg2)
{
  if (gCell->Parent != NULL)
  {
    int direction0 = gCell->Parent->Dir;
    int direction1;
    if (gCell->IsAlign == 1)
      direction1 = gCell->Parent->Dir;
    else
    {
      int ang = gCell->RefAngle;
      if (ang >= 7.0*M_PI/4.0 || ang < M_PI/4.0)
        direction1 = svCenterlineGraph::DT[direction0][0];

      else if (ang >= M_PI/4.0 && ang < 3.0*M_PI/4.0)
        direction1 = svCenterlineGraph::DT[direction0][1];

      else if (ang >= 3.0*M_PI/4.0 && ang < 5.0*M_PI/4.0)
        direction1 = svCenterlineGraph::DT[direction0][2];

      else if (ang >= 5.0*M_PI/4.0 && ang < 7.0*M_PI/4.0)
        direction1 = svCenterlineGraph::DT[direction0][3];
    }
    gCell->Dir = direction1;
    double dirVector[3];
    svCenterlineGraph::GetDirectionVector(direction1, dirVector);
    for (int i=0; i<3; i++)
    {
      gCell->StartPt[i] = gCell->Parent->EndPt[i];
      gCell->EndPt[i]   = gCell->StartPt[i] + dirVector[i];
    }
    vtkIntArray *rotIndices = reinterpret_cast<vtkIntArray*>(arg0);
    if (rotIndices != NULL)
    {
      int cellIndices[8];
      int tmpIds[8];
      for (int i=0; i<8; i++)
      {
        cellIndices[i] = svCenterlineGraph::LookupIndex(gCell->Parent->Dir, gCell->Parent->Children[gCell->Parent->DivergingChild]->Dir, i);
        tmpIds[i] = gCell->CornerPtIds[i];
      }
      for (int i=0; i<8; i++)
      {
        //fprintf(stdout,"Working on cube %d\n", gCell->GroupId);
        //fprintf(stdout,"Inside graph %d is now becoming %d\n", gCell->CornerPtIds[cellIndices[i]], tmpIds[cellIndices[rotIndices->GetValue(i)]]);
        gCell->CornerPtIds[cellIndices[i]] = tmpIds[rotIndices->GetValue(cellIndices[i])];
      }
    }
  }
  return SV_OK;
}

// ----------------------
// NewCell
// ----------------------
svCenterlineGCell* svCenterlineGraph::NewCell(int a_GroupId, svCenterlineGCell *a_Parent)
{
  svCenterlineGCell *newCell = new svCenterlineGCell;
  newCell->Parent  = a_Parent;
  newCell->Id      = this->NumberOfCells++;
  newCell->GroupId = a_GroupId;
  if (a_Parent != NULL)
  {
    for (int i=0; i<3; i++)
      newCell->StartPt[i] = a_Parent->EndPt[i];
  }

  return newCell;
}

// ----------------------
// NewCell
// ----------------------
svCenterlineGCell* svCenterlineGraph::NewCell(int a_GroupId, int a_Dir, double a_StartPt[3], double a_EndPt[3])
{
  svCenterlineGCell *newCell = new svCenterlineGCell;
  newCell->Id      = this->NumberOfCells++;
  newCell->GroupId = a_GroupId;
  newCell->Dir     = a_Dir;
  for (int i=0; i<3; i++)
  {
    newCell->StartPt[i] = a_StartPt[i];
    newCell->EndPt[i]   = a_EndPt[i];
  }

  return newCell;
}

// ----------------------
// GetCell
// ----------------------
svCenterlineGCell* svCenterlineGraph::GetCell(const int findId)
{
  if (findId > this->NumberOfCells)
  {
    fprintf(stdout,"Id is larger than number of cells\n");
    return NULL;
  }
  return this->LookUp(this->Root, findId);
}

// ----------------------
// LookUp
// ----------------------
svCenterlineGCell* svCenterlineGraph::LookUp(svCenterlineGCell *lookCell, const int findId)
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
    if (lookCell->Children.size() != 0)
    {
      for (int i=0; i<lookCell->Children.size(); i++)
      {
        svCenterlineGCell* foundCell = this->LookUp(lookCell->Children[i], findId);
        if (foundCell != NULL)
          return foundCell;
      }
    }
  }

  return NULL;
}

// ----------------------
// SetLocalCoordinatesOnPoints
// ----------------------
int svCenterlineGraph::RefineGraphWithLocalCoordinates()
{
  int numSegs   = this->NumberOfCells;
  int numPoints = this->Lines->GetNumberOfPoints();

  // Set up local arrays
  vtkNew(vtkDoubleArray, localArrayX);
  vtkNew(vtkDoubleArray, localArrayY);
  vtkNew(vtkDoubleArray, localArrayZ);
  localArrayX->SetNumberOfComponents(3);
  localArrayX->SetNumberOfTuples(numPoints);
  localArrayY->SetNumberOfComponents(3);
  localArrayY->SetNumberOfTuples(numPoints);
  localArrayZ->SetNumberOfComponents(3);
  localArrayZ->SetNumberOfTuples(numPoints);
  for (int i=0; i<3; i++)
  {
    localArrayX->FillComponent(i, -1);
    localArrayY->FillComponent(i, -1);
    localArrayZ->FillComponent(i, -1);
  }

  double refVecs[3][3];
  for (int i=0; i<numSegs; i++)
  {
    // Corresponding GCell
    svCenterlineGCell *gCell = this->GetCell(i);

    // The front direction of segment
    for (int j=0; j<3; j++)
    {
      for (int k=0; k<3; k++)
        refVecs[j][k] = gCell->RefDirs[j][k];
    }

    // cell id in the vtkPolyData
    int cellId = this->Lines->GetCellData()->GetArray(
     this->GroupIdsArrayName.c_str())->LookupValue(gCell->GroupId);

    // Get Cell points
    vtkIdType npts, *pts;
    this->Lines->GetCellPoints(cellId, npts, pts);

    int isTerminating = 0;
    double endVecs[3][3];
    if (gCell->Parent == NULL)
    {
      // found parent, flip around
      isTerminating = 1;

      // Check to make sure not already flipped
      vtkNew(vtkIdList, pointCells);
      this->Lines->GetPointCells(pts[0], pointCells);
      if (pointCells->GetNumberOfIds() == 1)
        this->FlipLinePoints(this->Lines, cellId);
      this->Lines->GetCellPoints(cellId, npts, pts);
    }
    else if (gCell->Children.size() != 0)
    {
      for (int j=0; j<3; j++)
      {
        for (int k=0; k<3; k++)
          endVecs[j][k] = gCell->Children[gCell->DivergingChild]->RefDirs[j][k];
      }
    }
    else
      isTerminating = 1;

    // First vecs
    double tmpX[3];
    double pt0[3], pt1[3];

    //// Compute are temp coordinate system
    //fprintf(stdout,"THE FIRST VEC: %.6f %.6f %.6f\n", refVecs[1][0], refVecs[1][1], refVecs[1][2]);
    //fprintf(stdout,"CURIOUS: %.6f\n", vtkMath::Dot(refVecs[1], refVecs[0]));
    //this->ComputeLocalCoordinateSystem(refVecs[0], refVecs[1], tmpX, refVecs[2]);
    //for (int j=0; j<3; j++)
    //  refVecs[1][j] = tmpX[j];

    if (gCell->Parent == NULL)
    {
      localArrayX->SetTuple(pts[0], refVecs[1]);
      localArrayY->SetTuple(pts[0], refVecs[2]);
      localArrayZ->SetTuple(pts[0], refVecs[0]);
    }

    for (int j=1; j<npts; j++)
    {
      this->Lines->GetPoint(pts[j-1], pt0);
      this->Lines->GetPoint(pts[j], pt1);
      if (gCell->Parent == NULL)
        vtkMath::Subtract(pt1, pt0, refVecs[0]);
      else
        vtkMath::Subtract(pt0, pt1, refVecs[0]);
      vtkMath::Normalize(refVecs[0]);

      this->ComputeLocalCoordinateSystem(refVecs[0], refVecs[1], tmpX, refVecs[2]);
      for (int k=0; k<3; k++)
        refVecs[1][k] = tmpX[k];

      localArrayX->SetTuple(pts[j], refVecs[1]);
      localArrayY->SetTuple(pts[j], refVecs[2]);
      localArrayZ->SetTuple(pts[j], refVecs[0]);

    }

    if (isTerminating == 0)
    {
      //fprintf(stdout,"FIXING FOR NOT TERMINATING BRANCH\n");
      //// Must test to see if it actually was correct son!
      //int maxDir;
      //double maxDot = -0.1;
      //for (int j=0; j<3; j++)
      //{
      //  double compare = fabs(vtkMath::Dot(refVecs[j], endVecs[0]));
      //  if (compare > maxDot)
      //  {
      //    maxDot = compare;
      //    maxDir = j;
      //  }
      //}
      //fprintf(stdout,"Apparently the end vec aligns with: %d\n", maxDir);

      double checkVec[3], dotVec[3];
      for (int j=0; j<3; j++)
      {
        //checkVec[j] = refVecs[(maxDir+2)%3][j];
        //dotVec[j] =   refVecs[(maxDir+1)%3][j];
        checkVec[j] = refVecs[2][j];
        dotVec[j] =   refVecs[1][j];
      }

      // Get angles
      double angleVec0[3];
      vtkMath::Cross(endVecs[2], checkVec, angleVec0);
      double ang = atan2(vtkMath::Norm(angleVec0), vtkMath::Dot(endVecs[2], checkVec));

      if (vtkMath::Dot(endVecs[2], dotVec) >= 0.0)
        ang = 2*M_PI - ang;

      fprintf(stdout,"WHAT THE ANG: %.4f\n", 180*ang/M_PI);

      gCell->Children[gCell->DivergingChild]->RefAngle = ang;

      // Get update that needs to happen
      int maxDir;
      double maxDot = -0.1;
      for (int j=0; j<3; j++)
      {
        double compare = fabs(vtkMath::Dot(endVecs[j], refVecs[1]));
        fprintf(stdout,"Dot with Ref %d: %.4f\n", j, compare);
        if (compare > maxDot)
        {
          maxDot = compare;
          maxDir = j;
        }
      }

      fprintf(stdout,"K WAHT THE EFF group %d with children %d %d\n", gCell->GroupId, gCell->Children[0]->GroupId, gCell->Children[1]->GroupId);
      fprintf(stdout,"Apparently the ref vec aligns with dir: %d\n", maxDir);

      double projVec[3];
      for (int j=0; j<3; j++)
        projVec[j] = endVecs[maxDir][j];

      vtkMath::Normalize(projVec);
      double projDot = vtkMath::Dot(refVecs[1], projVec);
      fprintf(stdout,"FIRST DOT: %.4f\n", projDot);
      fprintf(stdout,"DOUBLE CHECK ME!!\n");
      if (maxDir == 1)
        fprintf(stdout,"OKEY SMOKES: %.4f\n", vtkMath::Dot(refVecs[2], endVecs[2]));
      else
        fprintf(stdout,"OKEY SMOKES: %.4f\n", vtkMath::Dot(refVecs[2], endVecs[1]));

      vtkMath::MultiplyScalar(projVec, projDot);
      vtkMath::Normalize(projVec);
      fprintf(stdout,"What is the ref vec then: %.6f %.6f %.6f\n", refVecs[1][0], refVecs[1][1], refVecs[1][2]);
      fprintf(stdout,"This gives us the following vector to proj to: %.6f %.6f %.6f\n", endVecs[maxDir][0], endVecs[maxDir][1], endVecs[maxDir][2]);
      fprintf(stdout,"This gives us the following vector to proj to: %.6f %.6f %.6f\n", projVec[0], projVec[1], projVec[2]);
      double angleVec1[3];
      vtkMath::Cross(refVecs[1], projVec, angleVec1);
      double updateAngle = atan2(vtkMath::Norm(angleVec1), vtkMath::Dot(refVecs[1], projVec));
      fprintf(stdout,"INITIAL DIFF: %.6f\n", 180.0*updateAngle/M_PI);
      updateAngle = (180.0*updateAngle/M_PI) * (1./(npts-1));

      if (updateAngle > 45.0)
        fprintf(stdout,"ERROR: Angle cannot be larger than 45\n");

      double dotCheck;
      if (maxDir == 1)
      {
        dotCheck =  vtkMath::Dot(refVecs[1], endVecs[2]);
        if (projDot > 0)
        {
          if (dotCheck < 0)
            updateAngle *= -1.0;
        }
        else
        {
          if (dotCheck > 0)
            updateAngle *= -1.0;
        }
      }
      else if (maxDir == 2)
      {
        dotCheck = vtkMath::Dot(refVecs[1], endVecs[1]);
        if (projDot > 0)
        {
          if (dotCheck > 0)
            updateAngle *= -1.0;
        }
        else
        {
          if (dotCheck < 0)
            updateAngle *= -1.0;
        }
      }
      else
        fprintf(stdout,"TELL ME BECAUSE THIS SHOULDNT HAPPEN!!!!!!!\n");

      // Check to see if need to add negative to angle
      if (dotCheck > 0)
        updateAngle *= -1.0;
      fprintf(stdout,"OUR UPDATE ANGLE: %.4f\n", updateAngle);

      double updateRefVecs[3][3];
      // The front direction of segment
      for (int j=0; j<3; j++)
      {
        for (int k=0; k<3; k++)
          updateRefVecs[j][k] = gCell->RefDirs[j][k];
      }

      //Update now
      for (int j=1; j<npts; j++)
      {
        this->Lines->GetPoint(pts[j-1], pt0);
        this->Lines->GetPoint(pts[j], pt1);
        if (gCell->Parent == NULL)
          vtkMath::Subtract(pt1, pt0, updateRefVecs[0]);
        else
          vtkMath::Subtract(pt0, pt1, updateRefVecs[0]);
        vtkMath::Normalize(updateRefVecs[0]);

        double updateVec[3];
        this->RotateVecAroundLine(updateRefVecs[1], updateAngle, updateRefVecs[0], updateVec);
        this->ComputeLocalCoordinateSystem(updateRefVecs[0], updateVec, tmpX, updateRefVecs[2]);
        for (int k=0; k<3; k++)
          updateRefVecs[1][k] = tmpX[k];

        localArrayX->SetTuple(pts[j], updateRefVecs[1]);
        localArrayY->SetTuple(pts[j], updateRefVecs[2]);
        localArrayZ->SetTuple(pts[j], updateRefVecs[0]);

      }

      fprintf(stdout,"ENDED AND WANT TO CHECK: %.6f\n", vtkMath::Dot(updateRefVecs[1], endVecs[maxDir]));
      for (int j=0; j<3; j++)
      {
        for (int k=0; k<3; k++)
        {
          for (int l=0; l<gCell->Children.size(); l++)
            gCell->Children[l]->RefDirs[j][k] = updateRefVecs[j][k];
        }
      }

    }
  }

  localArrayX->SetName("LocalX");
  localArrayY->SetName("LocalY");
  localArrayZ->SetName("LocalZ");
  this->Lines->GetPointData()->AddArray(localArrayX);
  this->Lines->GetPointData()->AddArray(localArrayY);
  this->Lines->GetPointData()->AddArray(localArrayZ);
  return SV_OK;
}

// ----------------------
// ComputeLocalCoordinateSystem
// ----------------------
int svCenterlineGraph::ComputeLocalCoordinateSystem(const double vz[3],
                                                               const double vstart[3],
                                                               double vx[3],
                                                               double vy[3])
{
	double tempArray[3];
  double tempLength = vtkMath::Dot(vstart, vz);
	for (int i = 0; i < 3; i++)
		tempArray[i] = tempLength * vz[i];

  vtkMath::Subtract(vstart, tempArray, vx);
  vtkMath::Normalize(vx);

  vtkMath::Cross(vz, vx, vy);
  vtkMath::Normalize(vy);

  return SV_OK;
}

// ----------------------
// RotateVecAroundLine
// ----------------------
int svCenterlineGraph::RotateVecAroundLine(const double inVec[3],
                                           const double angle,
                                           const double axis[3],
                                           double outVec[3])
{
  double inPt[3];
  vtkMath::Subtract(inVec, axis, inPt);

  vtkNew(vtkPoints, tmpPoints);
  tmpPoints->InsertNextPoint(inPt);

  vtkNew(vtkPolyData, tmpPoly);
  tmpPoly->SetPoints(tmpPoints);

  vtkNew(vtkTransform, transform);
  //transform->RotateWXYZ(double angle, double x, double y, double z);
  transform->RotateWXYZ(angle, axis);

  vtkNew(vtkTransformPolyDataFilter, transformFilter);

  transformFilter->SetTransform(transform);
  transformFilter->SetInputData(tmpPoly);
  transformFilter->Update();

  double outPt[3];
  transformFilter->GetOutput()->GetPoint(0, outPt);

  vtkMath::Add(axis, outPt, outVec);

  return SV_OK;
}

// ----------------------
// FlipLinePoints
// ----------------------
int svCenterlineGraph::FlipLinePoints(vtkPolyData *pd, const int cellId)
{
  vtkIdType npts, *pts;
  pd->GetCellPoints(cellId, npts, pts);

  double *tmpPts = new double[npts];
  for (int i=0; i<npts; i++)
    tmpPts[i] = pts[i];

  for (int i=0; i<npts; i++)
    pts[(npts-1)-i] = tmpPts[i];

  pd->ReplaceCell(cellId, npts, pts);
  pd->Modified();
  pd->BuildLinks();

  return SV_OK;
}

// ----------------------
// LookupIndex
// ----------------------
/**
 * \details
 * Parent and diverging child combination to get the index based on
 * cube permutation groups. Uses the index look up table with combinations
 * of 90° rotations in the CCW direction areound the three reference axis
 * see http://www.euclideanspace.com/maths/discrete/groups/categorise/finite/cube/
 * <table>
 *  <caption id="multi_row">Rotation index lookup</caption>
 *  <tr><th> <th colspan="6">Diverging Child Dir. </tr>
 *  <tr><th>Parent Dir. <th>RIGHT <th>LEFT <th>FRONT <th>BACK <th>UP <th>DOWN </tr>
 *  <tr> <td>RIGHT <td>-1  <td>-1  <td>xyyy <td>xzzy <td>yyy  <td>zzy </tr>
 *  <tr> <td>LEFT  <td>-1  <td>-1  <td>xxxy <td>xy   <td>xxy  <td>y   </tr>
 *  <tr> <td>FRONT <td>xxx <td>xyy <td>-1   <td>-1   <td>zxxy <td>yxxx</tr>
 *  <tr> <td>BACK  <td>x   <td>yyx <td>-1   <td>-1   <td>yyyx <td>yx  </tr>
 *  <tr> <td>UP    <td>xx  <td>yy  <td>zxx  <td>zyy  <td>-1   <td>-1  </tr>
 *  <tr> <td>DOWN  <td>i   <td>zz  <td>zzz  <td>z    <td>-1   <td>-1  </tr>
 * </table>
 */

int svCenterlineGraph::LookupIndex(const int parent, const int divchild, const int index)
{
  if (parent == RIGHT)
  {
    if (divchild == BACK)
      return svCenterlineGraph::RT[0][svCenterlineGraph::RT[5][svCenterlineGraph::RT[1][index]]];
    if (divchild == FRONT)
      return svCenterlineGraph::RT[0][svCenterlineGraph::RT[7][index]];
    if (divchild == UP)
      return svCenterlineGraph::RT[7][index];
    if (divchild == DOWN)
      return svCenterlineGraph::RT[5][svCenterlineGraph::RT[1][index]];
  }
  if (parent == LEFT)
  {
    if (divchild == BACK)
      return svCenterlineGraph::RT[0][svCenterlineGraph::RT[1][index]];
    if (divchild == FRONT)
      return svCenterlineGraph::RT[6][svCenterlineGraph::RT[1][index]];
    if (divchild == UP)
      return svCenterlineGraph::RT[3][svCenterlineGraph::RT[1][index]];
    if (divchild == DOWN)
      return svCenterlineGraph::RT[1][index];
  }
  if (parent == FRONT)
  {
    if (divchild == RIGHT)
      return svCenterlineGraph::RT[6][index];
    if (divchild == LEFT)
      return svCenterlineGraph::RT[0][svCenterlineGraph::RT[4][index]];
    if (divchild == UP)
      return svCenterlineGraph::RT[2][svCenterlineGraph::RT[3][svCenterlineGraph::RT[1][index]]];
    if (divchild == DOWN)
      return svCenterlineGraph::RT[1][svCenterlineGraph::RT[6][index]];
  }
  if (parent == BACK)
  {
    if (divchild == RIGHT)
      return svCenterlineGraph::RT[0][index];
    if (divchild == LEFT)
      return svCenterlineGraph::RT[4][svCenterlineGraph::RT[0][index]];
    if (divchild == UP)
      return svCenterlineGraph::RT[7][svCenterlineGraph::RT[0][index]];
    if (divchild == DOWN)
      return svCenterlineGraph::RT[1][svCenterlineGraph::RT[1][index]];
  }
  if (parent == UP)
  {
    if (divchild == RIGHT)
      return svCenterlineGraph::RT[3][index];
    if (divchild == LEFT)
      return svCenterlineGraph::RT[4][index];
    if (divchild == BACK)
      return svCenterlineGraph::RT[2][svCenterlineGraph::RT[4][index]];
    if (divchild == FRONT)
      return svCenterlineGraph::RT[2][svCenterlineGraph::RT[3][index]];
  }
  if (parent == DOWN)
  {
    if (divchild == RIGHT)
      return index;
    if (divchild == LEFT)
      return svCenterlineGraph::RT[5][index];
    if (divchild == BACK)
      return svCenterlineGraph::RT[2][index];
    if (divchild == FRONT)
      return svCenterlineGraph::RT[8][index];
  }

  return -1;
}

// ----------------------
// DT
// ----------------------
/**
 * \details 6x4 table with new graph directions based
 * on the parent direction (rows) and the angle made with the aligning
 * reference direction (columns)
 *
 * <table>
 *  <caption id="multi_row">Angle lookup table</caption>
 *  <tr><th> <th colspan="4">Angles </tr>
 *  <tr><th>Parent Dir. <th>315-45° <th>45-135° <th>135-225° <th>225-315° </tr>
 *  <tr> <td>RIGHT <td>FRONT <td>UP    <td>BACK  <td>DOWN  </tr>
 *  <tr> <td>LEFT  <td>BACK  <td>DOWN  <td>FRONT <td>UP    </tr>
 *  <tr> <td>FRONT <td>UP    <td>LEFT  <td>DOWN  <td>RIGHT </tr>
 *  <tr> <td>BACK  <td>DOWN  <td>RIGHT <td>UP    <td>LEFT  </tr>
 *  <tr> <td>UP    <td>LEFT  <td>BACK  <td>RIGHT <td>FRONT </tr>
 *  <tr> <td>DOWN  <td>RIGHT <td>FRONT <td>LEFT  <td>BACK  </tr>
 * </table>
 */
const int svCenterlineGraph::DT[6][4] =
  {{FRONT, UP,    BACK,  DOWN},
   {BACK,  DOWN,  FRONT, UP},
   {UP,    LEFT,  DOWN,  RIGHT},
   {DOWN,  RIGHT, UP,    LEFT},
   {LEFT,  BACK,  RIGHT, FRONT},
   {RIGHT, FRONT, LEFT,  BACK}};

// ----------------------
// RT
// ----------------------
/**
 * \details
 * \verbatim
 *       z  y
 *       | /
 *       0 -- x
 *                             2---------1
 *                            /|        /|
 *                           / |       / |
 *                          3---------0  |
 *                          |  |      |  |
 *                          |  6------|--5
 *                          | /       | /
 *                          |/        |/
 *                          7---------4
 * \endverbatim
 * 3 90° rotations and 8 indices give 24 permutations of the indices.
 * To make it easier, we do nine rows and rows 3-5 are two calls of row 0-2
 * and 6-8 are three calls of row 0-2
 * <table>
 *  <caption id="multi_row">Index rotation table</caption>
 *  <tr><th>Row # <th>Rotation <th>Angle <th>Axis
 *  <tr> <td>0 <td> x   <td>90°  CCW <td>1  </tr>
 *  <tr> <td>1 <td> y   <td>90°  CCW <td>2  </tr>
 *  <tr> <td>2 <td> z   <td>90°  CCW <td>0  </tr>
 *  <tr> <td>3 <td> xx  <td>180° CCW <td>1  </tr>
 *  <tr> <td>4 <td> yy  <td>180° CCW <td>2  </tr>
 *  <tr> <td>5 <td> zz  <td>180° CCW <td>0  </tr>
 *  <tr> <td>6 <td> xxx <td>270° CCW <td>1  </tr>
 *  <tr> <td>7 <td> yyy <td>270° CCW <td>2  </tr>
 *  <tr> <td>8 <td> zzz <td>270° CCW <td>0  </tr>
 * </table>
 */
const int svCenterlineGraph::RT[9][8] =
{{4, 0, 3, 7, 5, 1, 2, 6},  // x
 {4, 5, 1, 0, 7, 6, 2, 3},  // y
 {1, 2, 3, 0, 5, 6, 7, 4},  // z
 {5, 4, 7, 6, 1, 0, 3, 2},  // xx
 {7, 6, 5, 4, 3, 2, 1, 0},  // yy
 {2, 3, 0, 1, 6, 7, 4, 5},  // zz
 {1, 5, 6, 2, 0, 4, 7, 3},  // xxx
 {3, 2, 6, 7, 0, 1, 5, 4},  // yyy
 {3, 0, 1, 2, 7, 4, 5, 6}}; // zzz

