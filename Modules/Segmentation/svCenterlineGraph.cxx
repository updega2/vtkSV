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

svCenterlineGraph::svCenterlineGraph()
{
  this->NumberOfCells = 0;
  this->NumberOfNodes = 1; // The root
  this->Root = NULL;

  this->Lines = vtkPolyData::New();
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
  this->Root = this->NewCell(rootId, RIGHT);

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
    fprintf(stdout,"Only one centerline\n");
    //return SV_OK;
  }


  for (int i=0; i<connectingGroups->GetNumberOfIds(); i++)
    this->Root->Children.push_back(this->NewCell(connectingGroups->GetId(i), this->Root));

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
int svCenterlineGraph::GetPolycube(const double height, const double width, vtkUnstructuredGrid *outUg)
{
  int numSegs = this->NumberOfCells;

  vtkNew(vtkPoints, allPoints);
  vtkNew(vtkCellArray, allCells);
  vtkNew(vtkIntArray, localPtIds); localPtIds->SetName("LocalPointIds");
  vtkNew(vtkIntArray, groupIds); groupIds->SetName("GroupIds");
  vtkNew(vtkIntArray, patchIds); patchIds->SetName("PatchIds");
  vtkNew(vtkDoubleArray, textureCoordinates);
  textureCoordinates->SetNumberOfComponents(3);
  textureCoordinates->SetName("TextureCoordinates");

  for (int i=0; i<numSegs; i++)
  {
    // Corresponding GCell
    svCenterlineGCell *gCell = this->GetCell(i);

    // cell id in the vtkPolyData
    int groupId = gCell->GroupId;

    fprintf(stdout,"GROUP: %d\n", groupId);
    int cubeType;
    this->GetCubeType(gCell, cubeType);
    fprintf(stdout,"CUBE TYPE: %d\n", cubeType);

    vtkNew(vtkPoints, addPoints);

    if (cubeType == 0)
    {
      int numPoints = 8;

      double finalPts[8][3];

      // get vector towards top of cube
      double workVec[3];
      vtkMath::Cross(gCell->RefDirs[1], gCell->RefDirs[0], workVec);
      vtkMath::Normalize(workVec);

      // Get beginning points
      double endVec[3];
      vtkMath::Normalize(workVec);
      vtkMath::Cross(gCell->RefDirs[0], workVec, endVec);
      vtkMath::Normalize(endVec);
      vtkMath::MultiplyScalar(endVec, width/2.);

      // Get points to extrude up
      double topPt0[3], topPt1[3];
      vtkMath::Add(gCell->StartPt, endVec, topPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(gCell->StartPt, endVec, topPt1);

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
      vtkMath::Add(gCell->EndPt, endVec, endPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(gCell->EndPt, endVec, endPt1);

      // extrude down
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(endPt0, workVec, finalPts[0]);
      vtkMath::Add(endPt1, workVec, finalPts[3]);

      // extrude up
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(endPt0, workVec, finalPts[4]);
      vtkMath::Add(endPt1, workVec, finalPts[7]);

      for (int j=0; j<numPoints; j++)
        addPoints->InsertNextPoint(finalPts[j]);
    }
    else if (cubeType == 1)
    {
      int numPoints = 10;
      double vecs[3][3];
      double endPts[2][3];
      svCenterlineGCell *align = gCell->Children[gCell->AligningChild];
      svCenterlineGCell *diver = gCell->Children[gCell->DivergingChild];

      this->FormBifurcation(gCell->StartPt, gCell->EndPt,
                            diver->EndPt, diver->StartPt,
                            align->EndPt, align->StartPt,
                            gCell->EndPt,
                            width/2., vecs, endPts);

      // get vector towards top of cube
      double workVec[3];
      vtkMath::Cross(vecs[1], vecs[0], workVec);
      vtkMath::Normalize(workVec);

      // extrude up
      double finalPts[10][3];
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(gCell->EndPt, workVec, finalPts[4]);
      vtkMath::Add(endPts[0], workVec, finalPts[0]);
      vtkMath::Add(endPts[1], workVec, finalPts[3]);

      // extrude down
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(gCell->EndPt, workVec, finalPts[9]);
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
      vtkMath::Add(gCell->StartPt, endVec, topPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(gCell->StartPt, endVec, topPt1);

      // extrude down
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(topPt0, workVec, finalPts[6]);
      vtkMath::Add(topPt1, workVec, finalPts[7]);

      // extrude up
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(topPt0, workVec, finalPts[1]);
      vtkMath::Add(topPt1, workVec, finalPts[2]);

      for (int j=0; j<numPoints; j++)
        addPoints->InsertNextPoint(finalPts[j]);
    }
    else if (cubeType == 2)
    {
      int numPoints = 12;

      svCenterlineGCell *brother, *diver, *parent;
      parent = gCell->Parent;
      if (parent->Children[0]->Id == gCell->Id)
        brother = parent->Children[1];
      else
        brother = parent->Children[0];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      double vecs[3][3];
      double topPts[2][3];
      this->FormBifurcation(gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[2]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[8]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[7]);
        vtkMath::Add(topPts[1], workVec, finalPts[9]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[9]);
        vtkMath::Add(topPts[1], workVec, finalPts[7]);
      }

      svCenterlineGCell *align = gCell->Children[gCell->AligningChild];
      svCenterlineGCell *cDiver = gCell->Children[gCell->DivergingChild];

      double endPts[2][3];
      this->FormBifurcation(gCell->StartPt, gCell->EndPt,
                            cDiver->EndPt, cDiver->StartPt,
                            align->EndPt, align->StartPt,
                            gCell->EndPt,
                            width/2., vecs, endPts);

      // get vector towards top of cube
      vtkMath::Cross(vecs[1], vecs[0], workVec);
      vtkMath::Normalize(workVec);

      if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(workVec, -1.0);

      // extrude up
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(gCell->EndPt, workVec, finalPts[5]);
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
      vtkMath::Add(gCell->EndPt, workVec, finalPts[11]);
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
        addPoints->InsertNextPoint(finalPts[j]);

    }
    else if (cubeType == 3)
    {
      int numPoints = 12;

      svCenterlineGCell *brother, *diver, *parent;
      parent = gCell->Parent;
      if (parent->Children[0]->Id == gCell->Id)
        brother = parent->Children[1];
      else
        brother = parent->Children[0];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      double vecs[3][3];
      double topPts[2][3];
      this->FormBifurcation(gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[5]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[6]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[10]);
        vtkMath::Add(topPts[1], workVec, finalPts[2]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[2]);
        vtkMath::Add(topPts[1], workVec, finalPts[10]);
      }

      svCenterlineGCell *align = gCell->Children[gCell->AligningChild];
      svCenterlineGCell *cDiver = gCell->Children[gCell->DivergingChild];

      double endPts[2][3];
      this->FormBifurcation(gCell->StartPt, gCell->EndPt,
                            cDiver->EndPt, cDiver->StartPt,
                            align->EndPt, align->StartPt,
                            gCell->EndPt,
                            height/2., vecs, endPts);


      // get vector towards top of cube
      vtkMath::Cross(vecs[1], vecs[0], workVec);
      vtkMath::Normalize(workVec);

      if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(workVec, -1.0);

      // extrude up
      vtkMath::MultiplyScalar(workVec, width/2.);
      vtkMath::Add(gCell->EndPt, workVec, finalPts[4]);
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
      vtkMath::Add(gCell->EndPt, workVec, finalPts[7]);
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
        addPoints->InsertNextPoint(finalPts[j]);

    }
    else if (cubeType == 4)
    {
      int numPoints = 12;

      svCenterlineGCell *brother, *diver, *parent;
      parent = gCell->Parent;
      if (parent->Children[0]->Id == gCell->Id)
        brother = parent->Children[1];
      else
        brother = parent->Children[0];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      double vecs[3][3];
      double topPts[2][3];
      this->FormBifurcation(gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[2]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[9]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[8]);
        vtkMath::Add(topPts[1], workVec, finalPts[10]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[10]);
        vtkMath::Add(topPts[1], workVec, finalPts[8]);
      }

      svCenterlineGCell *align = gCell->Children[gCell->AligningChild];
      svCenterlineGCell *cDiver = gCell->Children[gCell->DivergingChild];

      double endPts[2][3];
      this->FormBifurcation(gCell->StartPt, gCell->EndPt,
                            cDiver->EndPt, cDiver->StartPt,
                            align->EndPt, align->StartPt,
                            gCell->EndPt,
                            height/2., vecs, endPts);

      // get vector towards top of cube
      vtkMath::Cross(vecs[1], vecs[0], workVec);
      vtkMath::Normalize(workVec);

      if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(workVec, -1.0);

      // extrude up
      vtkMath::MultiplyScalar(workVec, width/2.);
      vtkMath::Add(gCell->EndPt, workVec, finalPts[5]);
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
      vtkMath::Add(gCell->EndPt, workVec, finalPts[6]);
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
        addPoints->InsertNextPoint(finalPts[j]);
    }
    else if (cubeType == 5)
    {
      int numPoints = 12;

      svCenterlineGCell *brother, *diver, *parent;
      parent = gCell->Parent;
      if (parent->Children[0]->Id == gCell->Id)
        brother = parent->Children[1];
      else
        brother = parent->Children[0];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      double vecs[3][3];
      double topPts[2][3];
      this->FormBifurcation(gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[5]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[6]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
      {
        vtkMath::Add(topPts[0], workVec, finalPts[9]);
        vtkMath::Add(topPts[1], workVec, finalPts[2]);
      }
      else
      {
        vtkMath::Add(topPts[0], workVec, finalPts[2]);
        vtkMath::Add(topPts[1], workVec, finalPts[9]);
      }

      svCenterlineGCell *align = gCell->Children[gCell->AligningChild];
      svCenterlineGCell *cDiver = gCell->Children[gCell->DivergingChild];

      double endPts[2][3];
      this->FormBifurcation(gCell->StartPt, gCell->EndPt,
                            cDiver->EndPt, cDiver->StartPt,
                            align->EndPt, align->StartPt,
                            gCell->EndPt,
                            width/2., vecs, endPts);

      // get vector towards top of cube
      vtkMath::Cross(vecs[1], vecs[0], workVec);
      vtkMath::Normalize(workVec);

      if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
        vtkMath::MultiplyScalar(workVec, -1.0);

      // extrude up
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(gCell->EndPt, workVec, finalPts[4]);
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
      vtkMath::Add(gCell->EndPt, workVec, finalPts[11]);
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
        addPoints->InsertNextPoint(finalPts[j]);

    }
    else if (cubeType == 6)
    {
      int numPoints = 10;

      svCenterlineGCell *brother, *diver, *parent;
      parent = gCell->Parent;
      if (parent->Children[0]->Id == gCell->Id)
        brother = parent->Children[1];
      else
        brother = parent->Children[0];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      double vecs[3][3];
      double topPts[2][3];
      this->FormBifurcation(gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
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

      vtkMath::Add(gCell->StartPt, workVec, finalPts[2]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[7]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->EndPt, endVec, endPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(gCell->EndPt, endVec, endPt1);

      // extrude down
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(endPt0, workVec, finalPts[5]);
      vtkMath::Add(endPt1, workVec, finalPts[9]);

      // extrude up
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(endPt0, workVec, finalPts[0]);
      vtkMath::Add(endPt1, workVec, finalPts[4]);

      for (int j=0; j<numPoints; j++)
        addPoints->InsertNextPoint(finalPts[j]);
    }
    else if (cubeType == 7)
    {
      int numPoints = 10;

      svCenterlineGCell *brother, *diver, *parent;
      parent = gCell->Parent;
      if (parent->Children[0]->Id == gCell->Id)
        brother = parent->Children[1];
      else
        brother = parent->Children[0];
      diver = parent->Children[parent->DivergingChild];

      // Get ending bifurcation points
      double vecs[3][3];
      double topPts[2][3];
      this->FormBifurcation(gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[4]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->StartPt, workVec, finalPts[5]);
      if (gCell->BranchDir == RIGHT || gCell->BranchDir == BACK)
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
      vtkMath::Add(gCell->EndPt, endVec, endPt0);
      vtkMath::MultiplyScalar(endVec, -1.0);
      vtkMath::Add(gCell->EndPt, endVec, endPt1);

      // extrude down
      vtkMath::MultiplyScalar(workVec, height/2.);
      vtkMath::Add(endPt0, workVec, finalPts[9]);
      vtkMath::Add(endPt1, workVec, finalPts[3]);

      // extrude up
      vtkMath::MultiplyScalar(workVec, -1.0);
      vtkMath::Add(endPt0, workVec, finalPts[6]);
      vtkMath::Add(endPt1, workVec, finalPts[0]);

      for (int j=0; j<numPoints; j++)
        addPoints->InsertNextPoint(finalPts[j]);
    }

    this->AddBranchCube(addPoints, allCells, allPoints, gCell->GroupId,
                        localPtIds, groupIds, patchIds, textureCoordinates, cubeType);

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
  outUg->GetPointData()->SetTCoords(textureCoordinates);

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
// ComputeGlobalReferenceVectors
// ----------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int svCenterlineGraph::ComputeGlobalReferenceVectors(svCenterlineGCell *parent)
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
int svCenterlineGraph::ComputeBranchReferenceVectors(svCenterlineGCell *parent)
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
int svCenterlineGraph::GetInitialBranchDirections(svCenterlineGCell *parent)
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
      for (int j=1; j<3; j++)
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
    }
    else
      parent->Children[i]->BranchDir = RIGHT;
    fprintf(stdout,"Child %d got direction %d\n", parent->Children[i]->GroupId, parent->Children[i]->BranchDir);
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
svCenterlineGCell* svCenterlineGraph::NewCell(int a_GroupId, int a_BranchDir, double a_StartPt[3], double a_EndPt[3])
{
  svCenterlineGCell *newCell = new svCenterlineGCell;
  newCell->Id      =   this->NumberOfCells++;
  newCell->GroupId =   a_GroupId;
  newCell->BranchDir = a_BranchDir;
  for (int i=0; i<3; i++)
  {
    newCell->StartPt[i] = a_StartPt[i];
    newCell->EndPt[i]   = a_EndPt[i];
  }

  return newCell;
}

// ----------------------
// NewCell
// ----------------------
svCenterlineGCell* svCenterlineGraph::NewCell(int a_GroupId, int a_BranchDir)
{
  svCenterlineGCell *newCell = new svCenterlineGCell;
  newCell->Id      =   this->NumberOfCells++;
  newCell->GroupId =   a_GroupId;
  newCell->BranchDir = a_BranchDir;

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
// GetGraphPoints
// ----------------------
int svCenterlineGraph::GetGraphPoints()
{
  int numSegs = this->NumberOfCells;

  for (int i=0; i<numSegs; i++)
  {
    // Corresponding GCell
    svCenterlineGCell *gCell = this->GetCell(i);

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

    svCenterlineGCell *parent = gCell->Parent;
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

      svCenterlineGCell *grandParent = parent->Parent;
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
        this->GetCubeType(parent, parentType);

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
int svCenterlineGraph::GetGraphDirections()
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
int svCenterlineGraph::UpdateBranchDirs(svCenterlineGCell *gCell, const int updateDir)
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

  delete [] tmpPts;

  return SV_OK;
}

// ----------------------
// GetCubeType
// ----------------------
int svCenterlineGraph::GetCubeType(svCenterlineGCell *gCell, int &type)
{
  if (gCell->Parent == NULL)
  {
    if (gCell->Children.size() == 0)
      type = 0;
    else if (gCell->Children.size() == 2)
      type = 1;
    else
    {
      fprintf(stderr,"Cannot currently Polycube for more than bifurcations yet!!!\n");
      type = 8;
      return SV_ERROR;
    }
  }
  else if (gCell->Children.size() == 2)
  {
    if ((gCell->Parent->BranchDir + gCell->BranchDir)%2 == 0)
    {
      if ((gCell->BranchDir + gCell->Children[gCell->DivergingChild]->BranchDir)%2 == 0)
        type = 2;
      else
        type = 4;
    }
    else
    {
      if ((gCell->BranchDir + gCell->Children[gCell->DivergingChild]->BranchDir)%2 == 0)
        type = 3;
      else
        type = 5;
    }
  }
  else if (gCell->Children.size() == 0)
  {
    if ((gCell->BranchDir)%2 == 0)
      type = 6;
    else
      type = 7;
  }
  else
  {
    fprintf(stderr,"Cannot currently Polycube for more than bifurcations yet!!!\n");
    type = 8;
    return SV_ERROR;
  }

  return SV_OK;
}

// ----------------------
// FormBifurcation
// ----------------------
int svCenterlineGraph::FormBifurcation(const double pt0[3], const double pt1[3],
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
int svCenterlineGraph::GetBifurcationPoint(const double startPt[3],
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
// AddBranchCube
// ----------------------
int svCenterlineGraph::AddBranchCube(vtkPoints *newPoints,
                                     vtkCellArray *cellArray,
                                     vtkPoints *points,
                                     const int groupId,
                                     vtkIntArray *localPtIds,
                                     vtkIntArray *groupIds,
                                     vtkIntArray *patchIds,
                                     vtkDoubleArray *textureCoordinates,
                                     const int type)
{
  if (type == 0)
  {
    if (newPoints->GetNumberOfPoints() != 8)
    {
      fprintf(stderr,"Must be given 8 points for cube type 0!\n");
      return SV_ERROR;
    }
    int pI[8];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);

    vtkIdType face0[4] = {pI[4], pI[5], pI[1], pI[0]};
    vtkIdType face1[4] = {pI[7], pI[6], pI[5], pI[4]};
    vtkIdType face2[4] = {pI[3], pI[2], pI[6], pI[7]};
    vtkIdType face3[4] = {pI[0], pI[1], pI[2], pI[3]};
    vtkIdType face4[4] = {pI[2], pI[1], pI[5], pI[6]};
    vtkIdType face5[4] = {pI[4], pI[0], pI[3], pI[7]};

    cellArray->InsertNextCell(4, face0);
    cellArray->InsertNextCell(4, face1);
    cellArray->InsertNextCell(4, face2);
    cellArray->InsertNextCell(4, face3);
    cellArray->InsertNextCell(4, face4);
    cellArray->InsertNextCell(4, face5);

    for (int i=0; i<6; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }

  }
  else if (type == 1)
  {
    if (newPoints->GetNumberOfPoints() != 10)
    {
      fprintf(stderr,"Must be given 10 points for cube type 1!\n");
      return SV_ERROR;
    }
    int pI[10];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.5, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.5, 0.0, 0.0);

    vtkIdType face0[4] = {pI[5], pI[6], pI[1], pI[0]};
    vtkIdType face1[5] = {pI[8], pI[7], pI[6], pI[5], pI[9]};
    vtkIdType face2[4] = {pI[3], pI[2], pI[7], pI[8]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};
    vtkIdType face4[4] = {pI[2], pI[1], pI[6], pI[7]};

    cellArray->InsertNextCell(4, face0);
    cellArray->InsertNextCell(5, face1);
    cellArray->InsertNextCell(4, face2);
    cellArray->InsertNextCell(5, face3);
    cellArray->InsertNextCell(4, face4);

    for (int i=0; i<5; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (type == 2)
  {
    if (newPoints->GetNumberOfPoints() != 12)
    {
      fprintf(stderr,"Must be given 12 points for cube type 2!\n");
      return SV_ERROR;
    }
    int pI[12];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.5, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.5, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.5, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.5, 0.0, 0.0);

    vtkIdType face0[4] = {pI[6], pI[7], pI[1], pI[0]};
    vtkIdType face1[6] = {pI[10], pI[9], pI[8], pI[7], pI[6], pI[11]};
    vtkIdType face2[4] = {pI[4], pI[3], pI[9], pI[10]};
    vtkIdType face3[6] = {pI[0], pI[1], pI[2], pI[3], pI[4], pI[5]};

    cellArray->InsertNextCell(4, face0);
    cellArray->InsertNextCell(6, face1);
    cellArray->InsertNextCell(4, face2);
    cellArray->InsertNextCell(6, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (type == 3)
  {
    if (newPoints->GetNumberOfPoints() != 12)
    {
      fprintf(stderr,"Must be given 12 points for cube type 3!\n");
      return SV_ERROR;
    }
    int pI[12];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.5);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.5);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.5);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.5);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);

    vtkIdType face0[6] = {pI[8], pI[9], pI[5], pI[1], pI[0], pI[4]};
    vtkIdType face1[4] = {pI[11], pI[10], pI[9], pI[8]};
    vtkIdType face2[6] = {pI[3], pI[2], pI[6], pI[10], pI[11], pI[7]};
    vtkIdType face3[4] = {pI[0], pI[1], pI[2], pI[3]};

    cellArray->InsertNextCell(6, face0);
    cellArray->InsertNextCell(4, face1);
    cellArray->InsertNextCell(6, face2);
    cellArray->InsertNextCell(4, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (type == 4)
  {
    if (newPoints->GetNumberOfPoints() != 12)
    {
      fprintf(stderr,"Must be given 12 points for cube type 4!\n");
      return SV_ERROR;
    }
    int pI[12];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.5, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.5);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.5);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.5, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);

    vtkIdType face0[5] = {pI[7], pI[8], pI[1], pI[0], pI[5]};
    vtkIdType face1[5] = {pI[11], pI[10], pI[9], pI[8], pI[7]};
    vtkIdType face2[5] = {pI[4], pI[3], pI[10], pI[11], pI[6]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};

    cellArray->InsertNextCell(5, face0);
    cellArray->InsertNextCell(5, face1);
    cellArray->InsertNextCell(5, face2);
    cellArray->InsertNextCell(5, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (type == 5)
  {
    if (newPoints->GetNumberOfPoints() != 12)
    {
      fprintf(stderr,"Must be given 12 points for cube type 5!\n");
      return SV_ERROR;
    }
    int pI[12];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.5, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.5);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.5);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.5, 0.0, 0.0);

    vtkIdType face0[5] = {pI[7], pI[8], pI[5], pI[1], pI[0]};
    vtkIdType face1[5] = {pI[10], pI[9], pI[8], pI[7], pI[11]};
    vtkIdType face2[5] = {pI[3], pI[2], pI[6], pI[9], pI[10]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};

    cellArray->InsertNextCell(5, face0);
    cellArray->InsertNextCell(5, face1);
    cellArray->InsertNextCell(5, face2);
    cellArray->InsertNextCell(5, face3);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
  }
  else if (type == 6)
  {
    if (newPoints->GetNumberOfPoints() != 10)
    {
      fprintf(stderr,"Must be given 10 points for cube type 6!\n");
      return SV_ERROR;
    }
    int pI[10];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.5, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.5, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);

    vtkIdType face0[4] = {pI[5], pI[6], pI[1], pI[0]};
    vtkIdType face1[5] = {pI[9], pI[8], pI[7], pI[6], pI[5]};
    vtkIdType face2[4] = {pI[4], pI[3], pI[8], pI[9]};
    vtkIdType face3[5] = {pI[0], pI[1], pI[2], pI[3], pI[4]};
    vtkIdType face4[4] = {pI[5], pI[0], pI[4], pI[9]};

    cellArray->InsertNextCell(4, face0);
    cellArray->InsertNextCell(5, face1);
    cellArray->InsertNextCell(4, face2);
    cellArray->InsertNextCell(5, face3);
    cellArray->InsertNextCell(4, face4);

    for (int i=0; i<4; i++)
    {
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
    }
    groupIds->InsertNextTuple1(groupId);
    patchIds->InsertNextTuple1(5);
  }
  else if (type == 7)
  {
    if (newPoints->GetNumberOfPoints() != 10)
    {
      fprintf(stderr,"Must be given 10 points for cube type 7!\n");
      return SV_ERROR;
    }
    int pI[10];
    for (int i=0; i<newPoints->GetNumberOfPoints(); i++)
    {
      pI[i] = points->InsertNextPoint(newPoints->GetPoint(i));
      localPtIds->InsertNextTuple1(i);
    }
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 1.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 1.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.5);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.5);
    textureCoordinates->InsertNextTuple3(1.0, 0.0, 0.0);
    textureCoordinates->InsertNextTuple3(1.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 1.0, 0.0);
    textureCoordinates->InsertNextTuple3(0.0, 0.0, 0.0);

    vtkIdType face0[5] = {pI[6], pI[7], pI[4], pI[1], pI[0]};
    vtkIdType face1[4] = {pI[9], pI[8], pI[7], pI[6]};
    vtkIdType face2[5] = {pI[3], pI[2], pI[5], pI[8], pI[9]};
    vtkIdType face3[4] = {pI[0], pI[1], pI[2], pI[3]};
    vtkIdType face4[4] = {pI[6], pI[0], pI[3], pI[9]};

    cellArray->InsertNextCell(5, face0);
    cellArray->InsertNextCell(4, face1);
    cellArray->InsertNextCell(5, face2);
    cellArray->InsertNextCell(4, face3);
    cellArray->InsertNextCell(4, face4);

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
