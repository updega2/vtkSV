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

#include "vtkSVCenterlineGraph.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkGraphToPolyData.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSortDataArray.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVMathUtils.h"
#include "vtkSVIOUtils.h"
#include "vtkSmartPointer.h"
#include "vtkThreshold.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkUnstructuredGrid.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVCenterlineGraph);

// ----------------------
// Constructor
// ----------------------
vtkSVCenterlineGraph::vtkSVCenterlineGraph()
{
  this->NumberOfCells = 0;
  this->NumberOfNodes = 1; // The root
  this->Root = NULL;

  this->Lines = vtkPolyData::New();
}

// ----------------------
// Constructor
// ----------------------
vtkSVCenterlineGraph::vtkSVCenterlineGraph(int rootId,
                 vtkPolyData *linesPd,
                 std::string groupIdsArrayName)
{
  this->NumberOfCells = 0;
  this->NumberOfNodes = 1; // The root
  this->Root = new vtkSVCenterlineGCell(this->NumberOfCells++, rootId, RIGHT);

  this->Lines = vtkPolyData::New();
  this->Lines->DeepCopy(linesPd);
  this->Lines->BuildLinks();
  this->GroupIdsArrayName = groupIdsArrayName;
}

// ----------------------
// Destructor
// ----------------------
vtkSVCenterlineGraph::~vtkSVCenterlineGraph()
{
  if (this->Root != NULL)
  {
    this->Root->Delete();
    this->Root = NULL;
  }
  if (this->Lines != NULL)
  {
    this->Lines->Delete();
    this->Lines = NULL;
  }
}

// ----------------------
// Recurse
// ----------------------
int vtkSVCenterlineGraph::Recurse(vtkSVCenterlineGCell *rootGCell, int(*function)(vtkSVCenterlineGCell *currentGCell,
                   void *arg0, void *arg1, void *arg2),
                   void *rec_arg0, void *rec_arg1, void *rec_arg2)
{
  function(rootGCell, rec_arg0, rec_arg1, rec_arg2);
  if (rootGCell->Children.size() != 0)
  {
    for (int i=0; i<rootGCell->Children.size(); i++)
      vtkSVCenterlineGraph::Recurse(rootGCell->Children[i], function, rec_arg0, rec_arg1, rec_arg2);
  }
  return SV_OK;
}

// ----------------------
// PrintGraph
// ----------------------
int vtkSVCenterlineGraph::PrintGraph()
{
  vtkSVCenterlineGraph::Recurse(this->Root, vtkSVCenterlineGraph::PrintGCell, NULL, NULL, NULL);
  return SV_OK;
}

// ----------------------
// PrintGCell
// ----------------------
int vtkSVCenterlineGraph::PrintGCell(vtkSVCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2)
{
  fprintf(stdout, "GCell ID: %d\n", gCell->Id);
  fprintf(stdout, "GCell Group ID: %d\n", gCell->GroupId);
  fprintf(stdout, "Direction: %d\n", gCell->BranchDir);
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
int vtkSVCenterlineGraph::GetGraphPolyData(vtkPolyData *pd)
{
  vtkNew(vtkPoints, newPoints);
  vtkNew(vtkCellArray, newCells);
  vtkNew(vtkIntArray, groupIds);
  groupIds->SetName(this->GroupIdsArrayName.c_str());
  vtkSVCenterlineGraph::Recurse(this->Root, vtkSVCenterlineGraph::InsertGCellPoints, newPoints, newCells, groupIds);

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
int vtkSVCenterlineGraph::InsertGCellPoints(vtkSVCenterlineGCell *gCell, void *arg0, void *arg1, void *arg2)
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
int vtkSVCenterlineGraph::BuildGraph()
{
  fprintf(stdout,"Building graph!!!\n");

  vtkNew(vtkIdList, connectingGroups);
  this->GetConnectingLineGroups(this->Root->GroupId, connectingGroups);

  if (connectingGroups->GetNumberOfIds() == 0)
  {
    // Just only this guy
    fprintf(stdout,"Only one centerline\n");
    //return SV_OK;
  }


  for (int i=0; i<connectingGroups->GetNumberOfIds(); i++)
  {
    vtkSVCenterlineGCell *newCell = new vtkSVCenterlineGCell(this->NumberOfCells++,
                                                             connectingGroups->GetId(i),
                                                             this->Root);
    this->Root->Children.push_back(newCell);
  }

  this->NumberOfNodes += this->Root->Children.size()+1;

  this->ComputeGlobalReferenceVectors(this->Root);
  if (this->Root->Children.size() != 0)
    this->ComputeBranchReferenceVectors(this->Root);

  for (int i=0; i<this->Root->Children.size(); i++)
    this->GrowGraph(this->Root->Children[i]);

  if (this->GetGraphDirections() != SV_OK)
    return SV_ERROR;
  if (this->GetGraphPoints() != SV_OK)
    return SV_ERROR;

  return SV_OK;
}

// ----------------------
// GetPolycube
// ----------------------
int vtkSVCenterlineGraph::GetPolycube(const double height, const double width, vtkUnstructuredGrid *outUg)
{
  int numSegs = this->NumberOfCells;

  vtkNew(vtkPoints, allPoints);
  vtkNew(vtkCellArray, allCells);
  vtkNew(vtkIntArray, localPtIds); localPtIds->SetName("LocalPointIds");
  vtkNew(vtkIntArray, groupIds); groupIds->SetName("GroupIds");
  vtkNew(vtkIntArray, patchIds); patchIds->SetName("PatchIds");

  for (int i=0; i<numSegs; i++)
  {
    // Corresponding GCell
    vtkSVCenterlineGCell *gCell = this->GetCell(i);

    gCell->GetCubePoints(height, width, allPoints, allCells,
                         localPtIds, groupIds, patchIds);

  }

  outUg->SetPoints(allPoints);
  outUg->SetCells(VTK_POLYGON, allCells);

  // Get final patch ids by adding to group ids
  vtkNew(vtkIdList, groupVals);
  for (int i=0; i<outUg->GetNumberOfCells(); i++)
  {
    int groupVal = groupIds->GetTuple1(i);
    groupVals->InsertUniqueId(groupVal);
  }
  vtkSortDataArray::Sort(groupVals);
  int numGroups = groupVals->GetNumberOfIds();
  vtkNew(vtkIdList, addVals);
  addVals->SetNumberOfIds(numGroups);
  for (int i=0; i<numGroups; i++)
    addVals->SetId(i, 6*i);

  for (int i=0; i<outUg->GetNumberOfCells(); i++)
  {
    int patchVal = patchIds->GetTuple1(i);
    int groupVal = groupIds->GetTuple1(i);
    int newVal = patchVal + (addVals->GetId(groupVals->IsId(groupVal)));
    patchIds->SetTuple1(i, newVal);
  }

  outUg->GetPointData()->AddArray(localPtIds);
  outUg->GetCellData()->AddArray(groupIds);
  outUg->GetCellData()->AddArray(patchIds);

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
int vtkSVCenterlineGraph::GrowGraph(vtkSVCenterlineGCell *parent)
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
      vtkSVCenterlineGCell *newCell = new vtkSVCenterlineGCell(this->NumberOfCells++,
                                                               groupId,
                                                               parent);
      parent->Children.push_back(newCell);
      count++;
    }
  }
  if (count > 0)
  {
    this->NumberOfNodes += parent->Children.size();
    //fprintf(stdout,"What is the child %d\n", parent->Children[0]->GroupId);
    //fprintf(stdout,"What is the child %d\n", parent->Children[1]->GroupId);
    this->ComputeBranchReferenceVectors(parent);

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
int vtkSVCenterlineGraph::GetConnectingLineGroups(const int groupId, vtkIdList *connectingGroups)
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
// ComputeGlobalReferenceVectors
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVCenterlineGraph::ComputeGlobalReferenceVectors(vtkSVCenterlineGCell *parent)
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
  if (numChildren != 0)
  {
    for (int i=0; i<numChildren; i++)
    {
      thresholder->SetInputData(this->Lines);
      thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
      thresholder->ThresholdBetween(parent->Children[i]->GroupId, parent->Children[i]->GroupId);
      thresholder->Update();

      double startPts[3], secondPts[3];
      thresholder->GetOutput()->GetPoint(0, startPts);
      thresholder->GetOutput()->GetPoint(1, secondPts);

      vtkMath::Subtract(secondPts, startPts, parent->Children[i]->BranchVec);
      vtkMath::Normalize(parent->Children[i]->BranchVec);

      // Get angle between vectors
      double angleVec[3];
      vtkMath::Cross(parent->Children[i]->BranchVec, this->ReferenceVecs[0], angleVec);
      parent->Children[i]->RefAngle = atan2(vtkMath::Norm(angleVec), vtkMath::Dot(parent->Children[i]->BranchVec, this->ReferenceVecs[0]));
    }

    double minAngle = VTK_SV_LARGE_DOUBLE;
    double maxAngle = -1.0*VTK_SV_LARGE_DOUBLE;
    int minChild = 0;
    int maxChild = 0;
    for (int i=0; i<numChildren; i++)
    {
      if (parent->Children[i]->RefAngle < minAngle)
      {
        minAngle = parent->Children[i]->RefAngle;
        minChild = i;
      }
      if (parent->Children[i]->RefAngle > maxAngle)
      {
        maxAngle = parent->Children[i]->RefAngle;
        maxChild = i;
      }
      fprintf(stdout,"Vec %d: %.4f %.4f %.4f\n", i, parent->Children[i]->BranchVec[0], parent->Children[i]->BranchVec[1], parent->Children[i]->BranchVec[2]);
      fprintf(stdout,"Angle %d: %4f\n", i, 180*parent->Children[i]->RefAngle/M_PI);
    }
    vtkMath::Cross(this->ReferenceVecs[0], parent->Children[minChild]->BranchVec, this->ReferenceVecs[2]);
    vtkMath::Normalize(this->ReferenceVecs[2]);

    parent->DivergingChild = minChild;
    parent->AligningChild  = maxChild;

    parent->Children[parent->AligningChild]->IsAlign = 1;
    for (int i=0; i<numChildren; i++)
    {
      if (i != parent->AligningChild)
        parent->Children[i]->IsAlign = 0;
    }
  }
  else
  {
    double xVec[3]; xVec[0] = 1.0; xVec[1] = 0.0; xVec[2] = 0.0;
    vtkMath::Cross(this->ReferenceVecs[0], xVec, this->ReferenceVecs[2]);
    vtkMath::Normalize(this->ReferenceVecs[2]);
  }

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
// ComputeBranchReferenceVectors
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVCenterlineGraph::ComputeBranchReferenceVectors(vtkSVCenterlineGCell *parent)
{
  //fprintf(stdout,"Child %d of parent %d, dir: %d\n", parent->Children[0]->GroupId, parent->GroupId, parent->Children[0]->BranchDir);
  //fprintf(stdout,"Child %d of parent %d, dir: %d\n", parent->Children[1]->GroupId, parent->GroupId, parent->Children[1]->BranchDir);
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

  int numChildren = parent->Children.size();
  for (int i=0; i<numChildren; i++)
  {
    thresholder->SetInputData(this->Lines);
    thresholder->SetInputArrayToProcess(0, 0, 0, 1, this->GroupIdsArrayName.c_str());
    thresholder->ThresholdBetween(parent->Children[i]->GroupId, parent->Children[i]->GroupId);
    thresholder->Update();

    double startPts[3], secondPts[3];
    thresholder->GetOutput()->GetPoint(0, startPts);
    thresholder->GetOutput()->GetPoint(1, secondPts);

    vtkMath::Subtract(secondPts, startPts, parent->Children[i]->BranchVec);
    vtkMath::Normalize(parent->Children[i]->BranchVec);

    // Get angle between vectors
    double angleVec[3];
    vtkMath::Cross(parent->Children[i]->BranchVec, vec0, angleVec);
    parent->Children[i]->RefAngle = atan2(vtkMath::Norm(angleVec), vtkMath::Dot(parent->Children[i]->BranchVec, vec0));
  }

  double minAngle = VTK_SV_LARGE_DOUBLE;
  double maxAngle = -1.0*VTK_SV_LARGE_DOUBLE;
  int minChild = 0;
  int maxChild = 0;
  for (int i=0; i<numChildren; i++)
  {
    if (parent->Children[i]->RefAngle < minAngle)
    {
      minAngle = parent->Children[i]->RefAngle;
      minChild = i;
    }
    if (parent->Children[i]->RefAngle > maxAngle)
    {
      maxAngle = parent->Children[i]->RefAngle;
      maxChild = i;
    }
    fprintf(stdout,"Vec %d: %.4f %.4f %.4f\n", i, parent->Children[i]->BranchVec[0], parent->Children[i]->BranchVec[1], parent->Children[i]->BranchVec[2]);
    fprintf(stdout,"Angle %d: %4f\n", i, 180*parent->Children[i]->RefAngle/M_PI);
  }

  double vec3[3];
  vtkMath::Cross(vec0, parent->Children[minChild]->BranchVec, vec3);
  vtkMath::Normalize(vec3);

  parent->DivergingChild = minChild;
  parent->AligningChild  = maxChild;

  parent->Children[parent->AligningChild]->IsAlign  = 1;
  for (int i=0; i<numChildren; i++)
  {
    if (i != parent->AligningChild)
      parent->Children[i]->IsAlign = 0;
  }

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
   this->GetInitialBranchDirections(parent);

  return SV_OK;
}

// ----------------------
// GetInitialBranchDirections
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVCenterlineGraph::GetInitialBranchDirections(vtkSVCenterlineGCell *parent)
{
  fprintf(stdout,"Determing initial dirs for parent %d\n", parent->GroupId);
  int numChildren = parent->Children.size();

  // Get initial branch dirs
  for (int i=0; i<numChildren; i++)
  {
    if (i != parent->DivergingChild)
    {
      int maxDir=0;
      double maxDot = -1.0;

      int startDir = 1;
      if (numChildren > 2)
        startDir = 0;

      for (int j=startDir; j<3; j++)
      {
        double compare = fabs(vtkMath::Dot(parent->Children[i]->RefDirs[j], parent->Children[i]->BranchVec));
        fprintf(stdout,"Would like to see compare, dir: %d, dot: %.6f\n", j, compare);
        if (compare > maxDot)
        {
          maxDot = compare;
          maxDir = j;
        }
      }
      double dotCheck = vtkMath::Dot(parent->Children[i]->RefDirs[maxDir], parent->Children[i]->BranchVec);
      if (maxDir == 1)
      {
        if (dotCheck > 0)
          parent->Children[i]->BranchDir = RIGHT;
        else
          parent->Children[i]->BranchDir = LEFT;
      }
      else if (maxDir == 2)
      {
        if (dotCheck > 0)
          parent->Children[i]->BranchDir = BACK;
        else
          parent->Children[i]->BranchDir = FRONT;
      }
      else if (maxDir == 0)
      {
        // Trifurcation case
        parent->Children[i]->BranchDir = RIGHT;
      }
    }
    else
      parent->Children[i]->BranchDir = RIGHT;
    fprintf(stdout,"Child %d got direction %d\n", parent->Children[i]->GroupId, parent->Children[i]->BranchDir);
  }

  return SV_OK;
}


// ----------------------
// GetCell
// ----------------------
vtkSVCenterlineGCell* vtkSVCenterlineGraph::GetCell(const int findId)
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
vtkSVCenterlineGCell* vtkSVCenterlineGraph::LookUp(vtkSVCenterlineGCell *lookCell, const int findId)
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
        vtkSVCenterlineGCell* foundCell = this->LookUp(lookCell->Children[i], findId);
        if (foundCell != NULL)
          return foundCell;
      }
    }
  }

  return NULL;
}

// ----------------------
// GetGraphPoints
// ----------------------
int vtkSVCenterlineGraph::GetGraphPoints()
{
  int numSegs = this->NumberOfCells;

  for (int i=0; i<numSegs; i++)
  {
    // Corresponding GCell
    vtkSVCenterlineGCell *gCell = this->GetCell(i);

    // cell id in the vtkPolyData
    int cellId = this->Lines->GetCellData()->GetArray(
     this->GroupIdsArrayName.c_str())->LookupValue(gCell->GroupId);

    // Get Cell points
    vtkIdType npts, *pts;
    this->Lines->GetCellPoints(cellId, npts, pts);

    double length = 0.0;

    for (int j=1; j<npts; j++)
    {
      double pt0[3], pt1[3];
      this->Lines->GetPoint(pts[j-1], pt0);
      this->Lines->GetPoint(pts[j], pt1);

      length += vtkSVMathUtils::Distance(pt0, pt1);
    }

    vtkSVCenterlineGCell *parent = gCell->Parent;
    if (parent == NULL)
    {
      this->Lines->GetPoint(pts[0], gCell->EndPt);
      double lineDir[3];
      for (int j=0; j<3; j++)
        lineDir[j] = gCell->RefDirs[0][j];
      vtkMath::MultiplyScalar(lineDir, length);
      vtkMath::Add(gCell->EndPt, lineDir, gCell->StartPt);
    }
    else
    {
      for (int j=0; j<3; j++)
        gCell->StartPt[j] = parent->EndPt[j];

      // Get rotation matrix from cross to ref
      double rotateVec[3];
      vtkMath::Subtract(parent->StartPt, parent->EndPt, rotateVec);
      vtkMath::Normalize(rotateVec);

      double crossVec[3];
      double lineDir[3];

      vtkSVCenterlineGCell *grandParent = parent->Parent;
      if (grandParent == NULL)
      {
        double parentVec[3];
        vtkMath::Subtract(parent->StartPt, parent->EndPt, parentVec);
        vtkMath::Normalize(parentVec);

        vtkMath::Cross(parentVec, this->ReferenceVecs[1], crossVec);
        vtkMath::Normalize(crossVec);

        // Rotate vec around line
        if (gCell->BranchDir == RIGHT)
          this->RotateVecAroundLine(rotateVec, 180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
        else if (gCell->BranchDir == LEFT)
          this->RotateVecAroundLine(rotateVec, -180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
      }
      else
      {
        double grandParentVec[3];
        vtkMath::Subtract(grandParent->StartPt, grandParent->EndPt, grandParentVec);
        vtkMath::Normalize(grandParentVec);

        double parentVec[3];
        vtkMath::Subtract(parent->EndPt, parent->StartPt, parentVec);
        vtkMath::Normalize(parentVec);


        int parentType;
        parent->GetCubeType(parentType);

        if (parentType == 4 || parentType == 5)
        {
          double tempVec[3];
          if (grandParent->BranchDir == LEFT || grandParent->BranchDir == FRONT)
          {
            if (parent->BranchDir == RIGHT)
              vtkMath::Cross(grandParentVec, parentVec, tempVec);
            else
              vtkMath::Cross(parentVec, grandParentVec, tempVec);
          }
          else
          {
            if (parent->BranchDir == RIGHT)
              vtkMath::Cross(parentVec, grandParentVec, tempVec);
            else
              vtkMath::Cross(grandParentVec, parentVec, tempVec);
          }
          vtkMath::Normalize(tempVec);
          vtkMath::Cross(parentVec, tempVec, crossVec);
          vtkMath::Normalize(crossVec);

          // Rotate vec around line
          if (gCell->BranchDir == RIGHT)
            this->RotateVecAroundLine(rotateVec, 180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
          else if (gCell->BranchDir == LEFT)
            this->RotateVecAroundLine(rotateVec, -180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
          else if (gCell->BranchDir == BACK)
            this->RotateVecAroundLine(rotateVec, 180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
          else if (gCell->BranchDir == FRONT)
            this->RotateVecAroundLine(rotateVec, -180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
        }
        else
        {
          vtkMath::Cross(grandParentVec, parentVec, crossVec);
          vtkMath::Normalize(crossVec);
          if (gCell->BranchDir == parent->BranchDir)
            this->RotateVecAroundLine(rotateVec, 180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
          else
            this->RotateVecAroundLine(rotateVec, -180.0*gCell->RefAngle/SV_PI, crossVec, lineDir);
        }
      }
      vtkMath::Normalize(lineDir);

      // Get end point from new direction
      vtkMath::MultiplyScalar(lineDir, length);
      vtkMath::Add(gCell->StartPt, lineDir, gCell->EndPt);
    }
  }

  return SV_OK;
}

// ----------------------
// GetGraphDirections
// ----------------------
int vtkSVCenterlineGraph::GetGraphDirections()
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
    vtkSVCenterlineGCell *gCell = this->GetCell(i);

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

    //// Compute our temp coordinate system
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
      fprintf(stdout,"FIXING FOR NOT TERMINATING BRANCH\n");

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

      double projVec[3];
      for (int j=0; j<3; j++)
        projVec[j] = endVecs[maxDir][j];

      vtkMath::Normalize(projVec);
      double projDot = vtkMath::Dot(refVecs[1], projVec);

      vtkMath::MultiplyScalar(projVec, projDot);
      vtkMath::Normalize(projVec);
      double angleVec1[3];
      vtkMath::Cross(refVecs[1], projVec, angleVec1);
      double updateAngle = atan2(vtkMath::Norm(angleVec1), vtkMath::Dot(refVecs[1], projVec));
      updateAngle = (180.0*updateAngle/M_PI) * (1./(npts-1));

      if (updateAngle > 45.0)
      {
        fprintf(stdout,"ERROR: Angle cannot be larger than 45\n");
        return SV_ERROR;
      }

      fprintf(stdout,"MAX DIR: %d\n", maxDir);
      int updateDir;
      if (maxDir == 1)
      {
        double dotCheck =  vtkMath::Dot(refVecs[1], endVecs[2]);
        fprintf(stdout,"DOT CHECK: %.4f\n", dotCheck);
        if (projDot > 0)
        {
          updateDir = 0;
          if (dotCheck >= 0)
            updateAngle *= -1.0;
        }
        else
        {
          updateDir = 2;
          if (dotCheck < 0)
            updateAngle *= -1.0;
        }
      }
      else if (maxDir == 2)
      {
        double dotCheck = vtkMath::Dot(refVecs[1], endVecs[1]);
        fprintf(stdout,"DOT CHECK: %.4f\n", dotCheck);
        if (projDot > 0)
        {
          updateDir = 3;
          if (dotCheck < 0)
            updateAngle *= -1.0;
        }
        else
        {
          updateDir = 1;
          if (dotCheck >= 0)
            updateAngle *= -1.0;
        }
      }
      else
      {
        fprintf(stderr,"ERROR: This is incorrect, something wrong in end vectors\n");
        return SV_ERROR;
      }

      // update branch dirs
      this->UpdateBranchDirs(gCell, updateDir);

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

      double finalCheck = fabs(vtkMath::Dot(updateRefVecs[1], endVecs[maxDir]));
      fprintf(stdout,"Final check: %.4f\n", finalCheck);
      if (finalCheck < 0.9)
      {
        fprintf(stderr,"ERROR: The vector was not updated to correct ending vector\n");
        fprintf(stderr,"The checking dot is %.6f\n", finalCheck);
        return SV_ERROR;
      }
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
// UpdateBranchDirs
// ----------------------
int vtkSVCenterlineGraph::UpdateBranchDirs(vtkSVCenterlineGCell *gCell, const int updateDir)
{
  if (gCell->Children.size() != 0)
  {
    for (int j=0; j<gCell->Children.size(); j++)
    {
      gCell->Children[j]->BranchDir = (gCell->Children[j]->BranchDir + updateDir)%4;
      fprintf(stdout,"BRANCH %d updated to direction %d\n", gCell->Children[j]->GroupId, gCell->Children[j]->BranchDir);
      fprintf(stdout,"UPDATE WAS: %d\n", updateDir);
      //this->UpdateBranchDirs(gCell->Children[j], updateDir);
    }
  }

  return SV_OK;
}

// ----------------------
// ComputeLocalCoordinateSystem
// ----------------------
int vtkSVCenterlineGraph::ComputeLocalCoordinateSystem(const double vz[3],
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
int vtkSVCenterlineGraph::RotateVecAroundLine(const double inVec[3],
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
int vtkSVCenterlineGraph::FlipLinePoints(vtkPolyData *pd, const int cellId)
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

  delete [] tmpPts;

  return SV_OK;
}
