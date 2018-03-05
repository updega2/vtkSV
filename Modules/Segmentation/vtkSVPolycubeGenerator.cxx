/*=========================================================================
 *
 * Copyright (c) 2014-2015 The Regents of the University of California.
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

#include "vtkSVPolycubeGenerator.h"

#include "vtkAppendFilter.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkIdList.h"
#include "vtkIdFilter.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSortDataArray.h"

#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVMathUtils.h"
#include "vtkTriangle.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVPolycubeGenerator);

// ----------------------
// Constructor
// ----------------------
vtkSVPolycubeGenerator::vtkSVPolycubeGenerator()
{
  this->CenterlineGroupIdsArrayName = NULL;
  this->CenterlineRadiusArrayName   = NULL;
  this->CenterlineGraph             = NULL;

  this->WorkPd     = vtkPolyData::New();
  this->GraphPd    = vtkPolyData::New();
  this->SurfacePolycubePd = vtkPolyData::New();
  this->VolumePolycubeUg  = vtkUnstructuredGrid::New();

  this->PolycubeUnitLength = 0.0;
  this->PolycubeDivisions = 5;
}

// ----------------------
// Destructor
// ----------------------
vtkSVPolycubeGenerator::~vtkSVPolycubeGenerator()
{
  if (this->CenterlineGroupIdsArrayName != NULL)
  {
    delete [] this->CenterlineGroupIdsArrayName;
    this->CenterlineGroupIdsArrayName = NULL;
  }

  if (this->CenterlineRadiusArrayName != NULL)
  {
    delete [] this->CenterlineRadiusArrayName;
    this->CenterlineRadiusArrayName = NULL;
  }

  if (this->CenterlineGraph != NULL)
  {
    this->CenterlineGraph->Delete();
    this->CenterlineGraph = NULL;
  }

  if (this->SurfacePolycubePd != NULL)
  {
    this->SurfacePolycubePd->Delete();
    this->SurfacePolycubePd = NULL;
  }

  if (this->VolumePolycubeUg != NULL)
  {
    this->VolumePolycubeUg->Delete();
    this->VolumePolycubeUg = NULL;
  }
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVPolycubeGenerator::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

// ----------------------
// RequestData
// ----------------------
int vtkSVPolycubeGenerator::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  this->WorkPd->DeepCopy(input);

  // Prep work for filter
  if (this->PrepFilter() != SV_OK)
  {
    vtkErrorMacro("Prep of filter failed");
    return SV_ERROR;
  }

  // Run the filter
  if (this->RunFilter() != SV_OK)
  {
    vtkErrorMacro("Filter failed");
    return SV_ERROR;
  }

  output->DeepCopy(this->SurfacePolycubePd);

  return SV_OK;
}

// ----------------------
// PrepFilter
// ----------------------
int vtkSVPolycubeGenerator::PrepFilter()
{
  if (this->WorkPd->GetNumberOfPoints() == 0 ||
      this->WorkPd->GetNumberOfCells() == 0)
  {
    vtkErrorMacro(<< "Input does not contain points or cells.");
    return SV_ERROR;
  }

  if (!this->CenterlineGroupIdsArrayName)
  {
    vtkDebugMacro("Centerline GroupIds Array Name not given, setting to GroupIds");
    this->CenterlineGroupIdsArrayName = new char[strlen("GroupIds") + 1];
    strcpy(this->CenterlineGroupIdsArrayName, "GroupIds");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 1, this->CenterlineGroupIdsArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "CenterlineGroupIdsArray with name specified does not exist");
    return SV_OK;
  }

  if (!this->CenterlineRadiusArrayName)
  {
    vtkDebugMacro("Centerline Radius Array Name not given, setting to MaximumInscribedSphereRadius");
    this->CenterlineRadiusArrayName = new char[strlen("MaximumInscribedSphereRadius") + 1];
    strcpy(this->CenterlineRadiusArrayName, "MaximumInscribedSphereRadius");
  }

  if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 0, this->CenterlineRadiusArrayName) != SV_OK)
  {
    vtkErrorMacro(<< "CenterlineRadiusArray with name specified does not exist");
    return SV_OK;
  }

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVPolycubeGenerator::RunFilter()
{
  double polycubeSize;
  if (this->PolycubeUnitLength == 0.0)
  {
    this->GetApproximatePolycubeSize(polycubeSize);
    this->PolycubeUnitLength = polycubeSize/this->PolycubeDivisions;
  }
  else
  {
    polycubeSize = this->PolycubeUnitLength*this->PolycubeDivisions;
  }

  this->CenterlineGraph = new vtkSVCenterlineGraph(0, this->WorkPd,
                                                this->CenterlineGroupIdsArrayName);
  this->CenterlineGraph->SetCubeSize(polycubeSize);

  if (this->CenterlineGraph->BuildGraph() != SV_OK)
  {
    vtkErrorMacro("Unable to form graph of centerlines");
    return SV_ERROR;
  }

  this->CenterlineGraph->GetGraphPolyData(this->GraphPd);

  this->GetSurfacePolycube(polycubeSize);

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(this->SurfacePolycubePd);
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->SetAbsoluteTolerance(1.0e-6);
  cleaner->Update();

  this->SurfacePolycubePd->DeepCopy(cleaner->GetOutput());
  this->SurfacePolycubePd->BuildLinks();

  this->GetVolumePolycube();

  return SV_OK;
}

// ----------------------
// GetApproximatePolycubeSize
// ----------------------
int vtkSVPolycubeGenerator::GetApproximatePolycubeSize(double &polycubeSize)
{
  double avgRadius = 0.0;

  int numPoints = this->WorkPd->GetNumberOfPoints();
  vtkDataArray *radiusArray = this->WorkPd->GetPointData()->GetArray(
    this->CenterlineRadiusArrayName);

  for (int i=0; i<numPoints; i++)
    avgRadius += radiusArray->GetTuple1(i);

  avgRadius = avgRadius/numPoints;

  polycubeSize = 2*avgRadius;
  polycubeSize = avgRadius;

  return SV_OK;
}

// ----------------------
// GetSurfacePolycube
// ----------------------
int vtkSVPolycubeGenerator::GetSurfacePolycube(const double cubeSize)
{
  int numSegs = this->CenterlineGraph->NumberOfCells;

  vtkNew(vtkPoints, allPoints);
  vtkNew(vtkCellArray, allCells);
  vtkNew(vtkIntArray, groupIds); groupIds->SetName("GroupIds");
  vtkNew(vtkIntArray, patchIds); patchIds->SetName("PatchIds");

  for (int i=0; i<numSegs; i++)
  {
    // Corresponding GCell
    vtkSVCenterlineGCell *gCell = this->CenterlineGraph->GetCell(i);

    this->GetCubePoints(gCell, cubeSize, cubeSize, allPoints, allCells,
                         groupIds, patchIds);

  }

  vtkNew(vtkUnstructuredGrid, polycube);
  polycube->SetPoints(allPoints);
  polycube->SetCells(VTK_POLYGON, allCells);

  // Get final patch ids by adding to group ids
  vtkNew(vtkIdList, groupVals);
  for (int i=0; i<polycube->GetNumberOfCells(); i++)
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

  for (int i=0; i<polycube->GetNumberOfCells(); i++)
  {
    int patchVal = patchIds->GetTuple1(i);
    int groupVal = groupIds->GetTuple1(i);
    int newVal = patchVal + (addVals->GetId(groupVals->IsId(groupVal)));
    patchIds->SetTuple1(i, newVal);
  }

  polycube->GetCellData()->AddArray(groupIds);
  polycube->GetCellData()->AddArray(patchIds);
  fprintf(stdout,"NUM CELLS: %d\n", allCells->GetNumberOfCells());
  fprintf(stdout,"NUM POINTS: %d\n", allPoints->GetNumberOfPoints());

  vtkNew(vtkDataSetSurfaceFilter, surfacer);
  surfacer->SetInputData(polycube);
  surfacer->Update();

  this->SurfacePolycubePd->DeepCopy(surfacer->GetOutput());

  return SV_OK;
}

// ----------------------
// GetVolumePolycube
// ----------------------
int vtkSVPolycubeGenerator::GetVolumePolycube()
{
  // Get all group ids
  int numSegs = this->CenterlineGraph->NumberOfCells;

  int w_div, h_div, l_div;
  std::vector<std::vector<int> > whl_divs;
  for (int i=0; i<numSegs; i++)
  {
    vtkSVCenterlineGCell *gCell = this->CenterlineGraph->GetCell(i);

    int groupId = gCell->GroupId;

    // Extract surface of polycube
    vtkNew(vtkPolyData, branchPolycube);
    vtkSVGeneralUtils::ThresholdPd(this->SurfacePolycubePd, groupId, groupId, 1, this->CenterlineGroupIdsArrayName, branchPolycube);

    branchPolycube->BuildLinks();

    if (this->GetPolycubeDivisions(gCell, branchPolycube,
                                   w_div, h_div, l_div) != SV_OK)
    {
      vtkErrorMacro("Error getting polycube divisions");
      return SV_ERROR;
    }

    std::vector<int> whl(4);
    whl[0] = groupId;
    whl[1] = w_div;
    whl[2] = h_div;
    whl[3] = l_div;

    whl_divs.push_back(whl);
  }

  if (this->SynchronizePolycubeDivisions(whl_divs) != SV_OK)
  {
    vtkErrorMacro("Could not synchronize polycube divisions");
    return SV_ERROR;
  }

  vtkNew(vtkIntArray, polycubeDivisionsArray);
  polycubeDivisionsArray->SetNumberOfComponents(4);
  polycubeDivisionsArray->SetNumberOfTuples(numSegs);
  polycubeDivisionsArray->SetName("PolycubeDivisions");

  vtkNew(vtkAppendFilter, appender);

  for (int i=0; i<numSegs; i++)
  {
    vtkSVCenterlineGCell *gCell = this->CenterlineGraph->GetCell(i);

    int groupId = gCell->GroupId;

    // Extract surface of polycube
    vtkNew(vtkPolyData, branchPolycube);
    vtkSVGeneralUtils::ThresholdPd(this->SurfacePolycubePd, groupId, groupId, 1, this->CenterlineGroupIdsArrayName, branchPolycube);

    branchPolycube->BuildLinks();

    fprintf(stdout,"FORMING PARA VOLUME FOR GROUP %d\n", groupId);
    vtkNew(vtkStructuredGrid, paraHexMesh);
    if (this->FormParametricHexMesh(gCell, branchPolycube, paraHexMesh, whl_divs[i][1],
                                    whl_divs[i][2], whl_divs[i][3]) != SV_OK)
    {
      fprintf(stderr,"Couldn't do the dirt\n");
      return SV_ERROR;
    }

    polycubeDivisionsArray->SetTuple4(i, groupId, whl_divs[i][1], whl_divs[i][2], whl_divs[i][3]);

    vtkNew(vtkIntArray, groupIdsArray);
    groupIdsArray->SetNumberOfTuples(paraHexMesh->GetNumberOfCells());
    groupIdsArray->SetName(this->CenterlineGroupIdsArrayName);
    groupIdsArray->FillComponent(0, groupId);

    paraHexMesh->GetCellData()->AddArray(groupIdsArray);

    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(paraHexMesh);
    ider->SetIdsArrayName("TmpInternalIds");
    ider->Update();

    appender->AddInputData(ider->GetOutput());
  }
  appender->Update();

  this->VolumePolycubeUg->DeepCopy(appender->GetOutput());
  this->VolumePolycubeUg->GetFieldData()->AddArray(polycubeDivisionsArray);

  return SV_OK;
}

// ----------------------
// GetCubePoints
// ----------------------
int vtkSVPolycubeGenerator::GetCubePoints(vtkSVCenterlineGCell *gCell,
                                          const double height,
                                          const double width,
                                          vtkPoints *allPoints,
                                          vtkCellArray *allCells,
                                          vtkIntArray *groupIds,
                                          vtkIntArray *patchIds)
{
  fprintf(stdout,"CUBE GROUP ID: %d\n", gCell->GroupId);

  int groupId = gCell->GroupId;

  vtkSVCenterlineGCell *brother, *diver, *parent, *align, *cDiver, *sister;
  parent = gCell->Parent;

  int myLoc;
  if (parent != NULL)
  {
    if (parent->Children[0]->Id == gCell->Id)
      brother = parent->Children[1];
    else
      brother = parent->Children[0];
    diver = parent->Children[parent->DivergingChild];

    for (int i=0; i<parent->Children.size(); i++)
    {
      if (parent->Children[i]->Id == gCell->Id)
        myLoc = i;
    }
  }
  else
    myLoc = -1;

  if (gCell->Children.size() > 0)
  {
    align = gCell->Children[gCell->AligningChild];
    cDiver = gCell->Children[gCell->DivergingChild];
  }


  // Get beginning type
  int begType, begSplitType;
  if (gCell->GetBeginningType(begType, begSplitType) != SV_OK)
  {
    return SV_ERROR;
  }
  fprintf(stdout,"BEG SPLIT TYPE: %d\n", begSplitType);
  fprintf(stdout,"BEG ORIEN TYPE: %d\n", begType);

  // Do beginning
  double begVec[3];
  double vecs[3][3];
  double triPts[2][3];
  vtkNew(vtkPoints, begPoints);

  if (begSplitType == ZERO)
  {
    fprintf(stdout,"BEG NONE\n");
    // Get top square from start pt
    vtkMath::Cross(gCell->RefDirs[1], gCell->RefDirs[0], begVec);
    vtkMath::Normalize(begVec);

    this->GetSquare(gCell->StartPt, gCell->RefDirs[1], begVec,
                    height, width, begPoints);

  }
  else if (begSplitType == BI)
  {
    if (begType >= C_TET_0 && begType <= C_TET_3)
    {
      fprintf(stdout,"BEG BI CORNER TET\n");
      // Get ending bifurcation points
      this->FormBifurcation(gCell, gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
                           width/2., vecs, triPts);

      // get diverging vector towards top of cube
      double myVec[3];
      vtkMath::Subtract(gCell->EndPt, gCell->StartPt, myVec);
      vtkMath::Normalize(myVec);

      // get vectors to form square
      vtkMath::Cross(myVec, vecs[1], begVec);
      vtkMath::Normalize(begVec);

      if (gCell->BranchDir == LEFT|| gCell->BranchDir == FRONT)
        vtkMath::MultiplyScalar(begVec, -1.0);

      if (gCell->BranchDir == BACK || gCell->BranchDir == FRONT)
      {
        vtkMath::MultiplyScalar(vecs[1], -1.0);
        this->GetSquare(gCell->StartPt, begVec, vecs[1],
                        height, width, begPoints);
      }
      else
      {
        this->GetSquare(gCell->StartPt, vecs[1], begVec,
                        height, width, begPoints);
      }

      double halfVec[3];
      vtkMath::Subtract(gCell->StartPt, gCell->EndPt, halfVec);
      vtkMath::Normalize(halfVec);
      vtkMath::MultiplyScalar(halfVec, width/2.);

      int pushBackInt;
      if (gCell->BranchDir == LEFT || gCell->BranchDir == BACK)
      {
        if ((gCell->BranchDir + 1)%4 == brother->BranchDir)
          pushBackInt = 3;
        else
          pushBackInt = 2;
      }
      else
      {
        if ((gCell->BranchDir + 1)%4 == brother->BranchDir)
          pushBackInt = 2;
        else
          pushBackInt = 3;
      }

      for (int i=0; i<4; i++)
      {
        double pt[3];
        begPoints->GetPoint(i, pt);

        double newPt[3];
        if (i == pushBackInt)
          vtkMath::Add(pt, halfVec, newPt);
        else
          vtkMath::Subtract(pt, halfVec, newPt);

        begPoints->SetPoint(i, newPt);
      }
    }
    else
    {
      fprintf(stdout,"BEG BI WEDGE\n");

      // Get ending bifurcation points
      this->FormBifurcation(gCell, gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
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
      if (gCell->BranchDir == LEFT || gCell->BranchDir == FRONT)
        this->GetWedge(triPts[1], gCell->StartPt, triPts[0], begVec,
                       height, begPoints);
      else
        this->GetWedge(triPts[0], gCell->StartPt, triPts[1], begVec,
                       height, begPoints);
    }
  }
  else if (begSplitType == TRI)
  {
    if (begType >= C_TET_0 && begType <= C_TET_3)
    {
      fprintf(stdout,"BEG TRI CORNER TET\n");

      if ((parent->Children[myLoc]->BranchDir + parent->Children[(myLoc+1)%3]->BranchDir)%2 == 0)
        brother = parent->Children[(myLoc+1)%3];
      else
        brother = parent->Children[(myLoc+2)%3];

      // Get ending bifurcation points
      this->FormBifurcation(gCell, gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
                           width/2., vecs, triPts);

      // get vectors to form square
      vtkMath::Cross(vecs[2], vecs[1], begVec);
      vtkMath::Normalize(begVec);

      if (brother->BranchDir == LEFT|| brother->BranchDir == FRONT)
        vtkMath::MultiplyScalar(begVec, -1.0);

      if (brother->BranchDir == BACK || brother->BranchDir == FRONT)
      {
        double broVec[3];
        if (brother->BranchDir == BACK)
        {
          for (int j=0; j<3; j++)
            broVec[j] = -1.0*vecs[2][j];
        }
        else
        {
          for (int j=0; j<3; j++)
            broVec[j] = 1.0*vecs[2][j];
        }
        this->GetSquare(gCell->StartPt, begVec, broVec,
                        height, width, begPoints);
      }
      else
      {
        double broVec[3];
        if (brother->BranchDir == LEFT)
        {
          for (int j=0; j<3; j++)
            broVec[j] = -1.0*vecs[2][j];
        }
        else
        {
          for (int j=0; j<3; j++)
            broVec[j] = 1.0*vecs[2][j];
        }
        this->GetSquare(gCell->StartPt, broVec, begVec,
                        height, width, begPoints);
      }

      double halfVec[3];
      vtkMath::Subtract(gCell->StartPt, gCell->EndPt, halfVec);
      vtkMath::Normalize(halfVec);
      vtkMath::MultiplyScalar(halfVec, width/2.);

      int pushBackInt;
      if (begType == C_TET_0)
        pushBackInt = 2;
      else if (begType == C_TET_1)
        pushBackInt = 1;
      else if (begType == C_TET_2)
        pushBackInt = 0;
      else
        pushBackInt = 3;

      for (int i=0; i<4; i++)
      {
        double pt[3];
        begPoints->GetPoint(i, pt);

        double newPt[3];
        if (i == pushBackInt)
        {
          for (int j=0; j<3; j++)
            newPt[j] = pt[j];
        }
        else
          vtkMath::Subtract(pt, halfVec, newPt);

        begPoints->SetPoint(i, newPt);
      }
    }
    else if (begType >= S_TET_0 && begType <= S_TET_3)
    {
      fprintf(stdout,"BEG TRI SIDE TET\n");

      brother = parent->Children[parent->AligningChild];

      // Get ending bifurcation points
      this->FormBifurcation(gCell, gCell->EndPt, gCell->StartPt,
                           parent->StartPt, parent->EndPt,
                           brother->EndPt, brother->StartPt,
                           gCell->StartPt,
                           width/2., vecs, triPts);

      double halfVec[3];
      vtkMath::Subtract(gCell->StartPt, gCell->EndPt, halfVec);
      vtkMath::Normalize(halfVec);
      vtkMath::MultiplyScalar(halfVec, width/2.);

      // get diverging vector towards top of cube
      double myVec[3];
      vtkMath::Subtract(gCell->EndPt, gCell->StartPt, myVec);
      vtkMath::Normalize(myVec);

      // get front vector
      vtkMath::Cross(myVec, vecs[1], begVec);
      vtkMath::Normalize(begVec);

      // Flip if diver branch is left or front
      if (gCell->BranchDir == LEFT || gCell->BranchDir == FRONT)
        vtkMath::MultiplyScalar(begVec, -1.0);

      if (gCell->BranchDir == BACK || gCell->BranchDir == FRONT)
      {
        double parVec[3];
        if (gCell->BranchDir == BACK)
        {
          for (int j=0; j<3; j++)
            parVec[j] = -1.0*vecs[1][j];
        }
        else
        {
          for (int j=0; j<3; j++)
            parVec[j] = 1.0*vecs[1][j];
        }
        this->GetSquare(gCell->StartPt, begVec, parVec,
                        height, width, begPoints);
      }
      else
      {
        double parVec[3];
        if (gCell->BranchDir == LEFT)
        {
          for (int j=0; j<3; j++)
            parVec[j] = -1.0*vecs[1][j];
        }
        else
        {
          for (int j=0; j<3; j++)
            parVec[j] = 1.0*vecs[1][j];
        }
        this->GetSquare(gCell->StartPt, parVec, begVec,
                        height, width, begPoints);
      }

      double midPt0[3], midPt1[3];
      double pt0[3], pt1[3];
      if (begType == S_TET_0)
      {
        begPoints->GetPoint(0, pt0);
        begPoints->GetPoint(1, pt1);
      }
      else if (begType == S_TET_1)
      {
        begPoints->GetPoint(0, pt0);
        begPoints->GetPoint(3, pt1);
      }
      else if (begType == S_TET_2)
      {
        begPoints->GetPoint(2, pt0);
        begPoints->GetPoint(3, pt1);
      }
      else
      {
        begPoints->GetPoint(1, pt0);
        begPoints->GetPoint(2, pt1);
      }
      vtkMath::Add(pt0, pt1, midPt0);
      vtkMath::MultiplyScalar(midPt0, 1./2);
      vtkMath::Add(midPt0, halfVec, midPt0);

      for (int i=0; i<4; i++)
      {
        double pt[3];
        begPoints->GetPoint(i, pt);

        double newPt[3];
        vtkMath::Subtract(pt, halfVec, newPt);
        begPoints->SetPoint(i, newPt);
      }

      if (begType == S_TET_0)
      {
        begPoints->GetPoint(2, pt0);
        begPoints->GetPoint(3, pt1);
      }
      else if (begType == S_TET_1)
      {
        begPoints->GetPoint(1, pt0);
        begPoints->GetPoint(2, pt1);
      }
      else if (begType == S_TET_2)
      {
        begPoints->GetPoint(0, pt0);
        begPoints->GetPoint(1, pt1);
      }
      else
      {
        begPoints->GetPoint(0, pt0);
        begPoints->GetPoint(3, pt1);
      }
      vtkMath::Add(pt0, pt1, midPt1);
      vtkMath::MultiplyScalar(midPt1, 1./2);

      begPoints->InsertNextPoint(midPt0);
      begPoints->InsertNextPoint(midPt1);
    }
    else
    {
      fprintf(stdout,"BEG TRI WEDGE\n");

      if (myLoc == 0)
      {
        brother = parent->Children[parent->AligningChild];
        diver = parent->Children[parent->DivergingChild];

        // Get ending bifurcation points
        this->FormBifurcation(gCell, gCell->EndPt, gCell->StartPt,
                              parent->StartPt, parent->EndPt,
                              brother->EndPt, brother->StartPt,
                              gCell->StartPt,
                              width/2., vecs, triPts);
      }
      else if (myLoc == parent->Children.size() - 1)
      {
        brother = parent->Children[parent->AligningChild];
        diver = parent->Children[parent->DivergingChild];

        // Get ending bifurcation points
        this->FormBifurcation(gCell, gCell->EndPt, gCell->StartPt,
                             parent->StartPt, parent->EndPt,
                             brother->EndPt, brother->StartPt,
                             gCell->StartPt,
                             width/2., vecs, triPts);
      }
      else
      {
        brother = parent->Children[myLoc - 1];
        diver = parent->Children[parent->AligningChild];
        sister  = parent->Children[myLoc + 1];

        // Get ending bifurcation points
        this->FormBifurcation(gCell, gCell->EndPt, gCell->StartPt,
                             brother->EndPt, brother->StartPt,
                             sister->EndPt, sister->StartPt,
                             gCell->StartPt,
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

      if (gCell->BranchDir == LEFT || gCell->BranchDir == FRONT)
        this->GetWedge(triPts[1], gCell->StartPt, triPts[0], begVec,
                       height, begPoints);
      else
        this->GetWedge(triPts[0], gCell->StartPt, triPts[1], begVec,
                       height, begPoints);

    }
  }
  else
  {
    fprintf(stdout,"BEG NOT_HANDLED\n");
  }

  // Get end type
  int endType, endSplitType;
  if (gCell->GetEndType(endType, endSplitType) != SV_OK)
  {
    return SV_ERROR;
  }
  fprintf(stdout,"END SPLIT TYPE: %d\n", endSplitType);
  fprintf(stdout,"END ORIEN TYPE: %d\n", endType);


  // Do end
  vtkNew(vtkPoints, endPoints);

  if (endSplitType == ZERO)
  {
    fprintf(stdout,"END NONE\n");

    // Get square from end pt
    if (begSplitType == ZERO)
      this->GetSquare(gCell->EndPt, gCell->RefDirs[1], begVec,
                      height, width, endPoints);
    else
    {
      double endVec[3];
      if ((gCell->BranchDir)%2 == 0)
      {
        vtkMath::Cross(begVec, vecs[0], endVec);
        vtkMath::Normalize(endVec);
        this->GetSquare(gCell->EndPt, endVec, begVec,
                        height, width, endPoints);
      }
      else
      {
        vtkMath::Cross(vecs[0], begVec, endVec);
        vtkMath::Normalize(endVec);
        this->GetSquare(gCell->EndPt, begVec, endVec,
                        height, width, endPoints);
      }
    }
  }
  else if (endSplitType == BI)
  {
    if (endType >= C_TET_0 && endType <= S_TET_3)
    {
      fprintf(stdout,"END BI TET\n");
      // Get square from end pt
      if (begSplitType == ZERO)
        this->GetSquare(gCell->EndPt, gCell->RefDirs[1], begVec,
                        height, width, endPoints);
      else
      {
        double endVec[3];
        if ((gCell->BranchDir)%2 == 0)
        {
          vtkMath::Cross(begVec, vecs[0], endVec);
          vtkMath::Normalize(endVec);
          this->GetSquare(gCell->EndPt, endVec, begVec,
                          height, width, endPoints);
        }
        else
        {
          vtkMath::Cross(vecs[0], begVec, endVec);
          vtkMath::Normalize(endVec);
          this->GetSquare(gCell->EndPt, begVec, endVec,
                          height, width, endPoints);
        }
      }
      double halfVec[3];
      vtkMath::Subtract(gCell->EndPt, gCell->StartPt, halfVec);
      vtkMath::Normalize(halfVec);
      vtkMath::MultiplyScalar(halfVec, width/2.);

      for (int i=0; i<4; i++)
      {
        double pt[3];
        endPoints->GetPoint(i, pt);

        double newPt[3];
        if (i == 2)
          vtkMath::Add(pt, halfVec, newPt);
        else
          vtkMath::Subtract(pt, halfVec, newPt);

        endPoints->SetPoint(i, newPt);
      }

    }
    else
    {
      fprintf(stdout,"END BI WEDGE\n");

      this->FormBifurcation(gCell, gCell->StartPt, gCell->EndPt,
                            cDiver->EndPt, cDiver->StartPt,
                            align->EndPt, align->StartPt,
                            gCell->EndPt,
                            width/2., vecs, triPts);


      // get vector towards top of cube
      vtkMath::Cross(vecs[1], vecs[0], begVec);
      vtkMath::Normalize(begVec);

      // treat differently for different dirs
      if (cDiver->BranchDir == LEFT || cDiver->BranchDir == FRONT)
      {
        vtkMath::MultiplyScalar(begVec, -1.0);
        this->GetWedge(triPts[0], gCell->EndPt, triPts[1], begVec,
                       height, endPoints);
      }
      else
        this->GetWedge(triPts[1], gCell->EndPt, triPts[0], begVec,
                       height, endPoints);
    }
  }
  else if (endSplitType == TRI)
  {
    if (endType >= C_TET_0 && endType <= C_TET_3)
    {
      fprintf(stdout,"END TRI CORNER TET\n");
      // Get square from end pt
      if (begSplitType == ZERO)
        this->GetSquare(gCell->EndPt, gCell->RefDirs[1], begVec,
                        height, width, endPoints);
      else
      {
        double endVec[3];
        if ((gCell->BranchDir)%2 == 0)
        {
          vtkMath::Cross(begVec, vecs[0], endVec);
          vtkMath::Normalize(endVec);
          this->GetSquare(gCell->EndPt, endVec, begVec,
                          height, width, endPoints);
        }
        else
        {
          vtkMath::Cross(vecs[0], begVec, endVec);
          vtkMath::Normalize(endVec);
          this->GetSquare(gCell->EndPt, begVec, endVec,
                          height, width, endPoints);
        }
      }
      double halfVec[3];
      vtkMath::Subtract(gCell->EndPt, gCell->StartPt, halfVec);
      vtkMath::Normalize(halfVec);
      vtkMath::MultiplyScalar(halfVec, width/2.);

      int pushBackInt;
      if (endType == C_TET_0)
        pushBackInt = 2;
      else if (endType == C_TET_1)
        pushBackInt = 1;
      else if (endType == C_TET_2)
        pushBackInt = 0;
      else
        pushBackInt = 3;

      for (int i=0; i<4; i++)
      {
        double pt[3];
        endPoints->GetPoint(i, pt);

        double newPt[3];
        if (i == pushBackInt)
        {
          for (int j=0; j<3; j++)
            newPt[j] = pt[j];
        }
        else
          vtkMath::Subtract(pt, halfVec, newPt);

        endPoints->SetPoint(i, newPt);
      }
    }
    else
    {
      fprintf(stdout,"END TRI WEDGE\n");

      int crossChild;
      for (int i=0; i<3; i++)
      {
        if (i != gCell->AligningChild && i != gCell->DivergingChild)
          crossChild = i;
      }
      vtkSVCenterlineGCell *cross = gCell->Children[crossChild];
      diver = gCell->Children[gCell->DivergingChild];

      double endPts[2][3];
      this->FormBifurcation(gCell, gCell->StartPt, gCell->EndPt,
                            diver->EndPt, diver->StartPt,
                            cross->EndPt, cross->StartPt,
                            gCell->EndPt,
                            width/2., vecs, endPts);

      // TEMPORARYRYRY
      if (diver->BranchDir == RIGHT || diver->BranchDir == BACK)
      {
        // get vector towards top of cube
        double frontVec[3];
        vtkMath::Cross(vecs[1], vecs[0], frontVec);
        vtkMath::Normalize(frontVec);

        this->GetWedge(endPts[1], gCell->EndPt, endPts[0], frontVec,
                        height, endPoints);
      }
      else if (diver->BranchDir == LEFT || diver->BranchDir == FRONT)
      {
        // get vector towards top of cube
        double frontVec[3];
        vtkMath::Cross(vecs[0], vecs[1], frontVec);
        vtkMath::Normalize(frontVec);

        this->GetWedge(endPts[0], gCell->EndPt, endPts[1], frontVec,
                        height, endPoints);
      }
      else
      {
        fprintf(stderr,"NEED TO ADD RULE!!\n");
      }
    }
  }
  else
  {
    fprintf(stderr,"END NOT_HANDLED\n");
  }

  // NOW THAT WE HAVE ALL THE POINTS, WE CAN NOW ADD THEM TO THEIR
  // RESPECTIVE FACES
  std::vector<std::vector<int> > facesPtIds(6);
  // -----------------------------------------------------------------------
  // FACE 3
  // -----------------------------------------------------------------------

  // END
  if (endType == NONE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(1)));
  }
  else if (endType == VERT_WEDGE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  }
  else if (endType == HORZ_WEDGE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));
  }
  else if (endType >= C_TET_0 && endType <= C_TET_3)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(1)));
  }
  else if (endType >= S_TET_0 && endType <= S_TET_3)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(1)));
  }

  // BEG
  if (begType == NONE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(1)));
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
  }
  else if (begType == VERT_WEDGE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(1)));
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
  }
  else if (begType == HORZ_WEDGE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(5)));
  }
  else if (begType >= C_TET_0 && begType <= C_TET_3)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(1)));
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
  }
  else if (begType >= S_TET_0 && begType <= S_TET_3)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(1)));
    if (begType == S_TET_3)
      facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(4)));
    if (begType == S_TET_1)
      facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(5)));
    facesPtIds[3].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(2)));
  }

  // END
  if (endType == NONE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  }
  else if (endType == VERT_WEDGE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(1)));
  }
  else if (endType == HORZ_WEDGE)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  }
  else if (endType >= C_TET_0 && endType <= C_TET_3)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  }
  else if (endType >= S_TET_0 && endType <= S_TET_3)
  {
    facesPtIds[3].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 1
  // -----------------------------------------------------------------------

  // END
  if (endType == NONE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  }
  else if (endType == VERT_WEDGE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  }
  else if (endType == HORZ_WEDGE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(5)));
  }
  else if (endType >= C_TET_0 && endType <= C_TET_3)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  }
  else if (endType >= S_TET_0 && endType <= S_TET_3)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(3)));
  }

  // BEG
  if (begType == NONE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
  }
  else if (begType == VERT_WEDGE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(5)));
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(4)));
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
  }
  else if (begType == HORZ_WEDGE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
  }
  else if (begType >= C_TET_0 && begType <= C_TET_3)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
  }
  else if (begType >= S_TET_0 && begType <= S_TET_3)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(3)));
    if (begType == S_TET_3)
      facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(5)));
    if (begType == S_TET_1)
      facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(4)));
    facesPtIds[1].push_back(allPoints->InsertNextPoint(begPoints->GetPoint(0)));
  }

  // END
  if (endType == NONE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));
  }
  else if (endType == VERT_WEDGE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(5)));
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(4)));
  }
  else if (endType == HORZ_WEDGE)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(2)));
  }
  else if (endType >= C_TET_0 && endType <= C_TET_3)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));
  }
  else if (endType >= S_TET_0 && endType <= S_TET_3)
  {
    facesPtIds[1].push_back(allPoints->InsertNextPoint(endPoints->GetPoint(0)));
  }
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

  // BEG TRI TET
  if (begType == S_TET_0 || begType == S_TET_2)
  {
    begInterPtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(4)));
    begInterPtIds.push_back(allPoints->InsertNextPoint(begPoints->GetPoint(5)));
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 0
  // -----------------------------------------------------------------------

  // END
  if (endType == NONE)
  {
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-1]);
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-2]);
  }
  else if (endType == VERT_WEDGE)
  {
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-2]);
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-3]);
  }
  else if (endType == HORZ_WEDGE)
  {
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-1]);
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-2]);
  }
  else if (endType >= C_TET_0 && endType <= S_TET_3)
  {
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-1]);
    facesPtIds[0].push_back(facesPtIds[1][facesPtIds[1].size()-2]);
  }

  // INT
  if (begType == HORZ_WEDGE)
    facesPtIds[0].push_back(begInterPtIds[0]);

  else if (begType == S_TET_0)
    facesPtIds[0].push_back(begInterPtIds[0]);

  else if (begType == S_TET_2)
    facesPtIds[0].push_back(begInterPtIds[1]);


  // BEG
  if (begType == NONE)
  {
    facesPtIds[0].push_back(facesPtIds[3][1]);
    facesPtIds[0].push_back(facesPtIds[3][0]);
  }
  else if (begType == VERT_WEDGE)
  {
    facesPtIds[0].push_back(facesPtIds[3][1]);
    facesPtIds[0].push_back(facesPtIds[3][0]);
  }
  else if (begType == HORZ_WEDGE)
  {
    facesPtIds[0].push_back(facesPtIds[3][1]);
    facesPtIds[0].push_back(facesPtIds[3][0]);
  }
  else if (begType >= C_TET_0 && begType <= S_TET_3)
  {
    facesPtIds[0].push_back(facesPtIds[3][1]);
    facesPtIds[0].push_back(facesPtIds[3][0]);
  }

  // INT
  if (endType == HORZ_WEDGE)
    facesPtIds[0].push_back(endInterPtIds[0]);
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 2
  // -----------------------------------------------------------------------

  // END
  if (endType == NONE)
  {
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-1]);
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-2]);
  }
  else if (endType == VERT_WEDGE)
  {
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-2]);
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-3]);
  }
  else if (endType == HORZ_WEDGE)
  {
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-1]);
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-2]);
  }
  else if (endType >= C_TET_0 && endType <= S_TET_3)
  {
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-1]);
    facesPtIds[2].push_back(facesPtIds[3][facesPtIds[3].size()-2]);
  }

  // INT
  if (begType == HORZ_WEDGE)
    facesPtIds[2].push_back(begInterPtIds[1]);

  else if (begType == S_TET_0)
    facesPtIds[2].push_back(begInterPtIds[1]);

  else if (begType == S_TET_2)
    facesPtIds[2].push_back(begInterPtIds[0]);

  // BEG
  if (begType == NONE)
  {
    facesPtIds[2].push_back(facesPtIds[1][1]);
    facesPtIds[2].push_back(facesPtIds[1][0]);
  }
  else if (begType == VERT_WEDGE)
  {
    facesPtIds[2].push_back(facesPtIds[1][1]);
    facesPtIds[2].push_back(facesPtIds[1][0]);
  }
  else if (begType == HORZ_WEDGE)
  {
    facesPtIds[2].push_back(facesPtIds[1][1]);
    facesPtIds[2].push_back(facesPtIds[1][0]);
  }
  else if (begType >= C_TET_0 && begType <= S_TET_3)
  {
    facesPtIds[2].push_back(facesPtIds[1][1]);
    facesPtIds[2].push_back(facesPtIds[1][0]);
  }

  // INT
  if (endType == HORZ_WEDGE)
    facesPtIds[2].push_back(endInterPtIds[1]);
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 4, if there
  // -----------------------------------------------------------------------

  // BEG
  if (begType == NONE)
  {
    facesPtIds[4].push_back(facesPtIds[3][2]);
    facesPtIds[4].push_back(facesPtIds[3][1]);
    facesPtIds[4].push_back(facesPtIds[1][2]);
    facesPtIds[4].push_back(facesPtIds[1][1]);
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // FACE 5, if there
  // -----------------------------------------------------------------------

  // BEG
  if (endType == NONE)
  {
    facesPtIds[5].push_back(facesPtIds[1][facesPtIds[1].size()-1]);
    facesPtIds[5].push_back(facesPtIds[3][0]);
    facesPtIds[5].push_back(facesPtIds[3][facesPtIds[3].size()-1]);
    facesPtIds[5].push_back(facesPtIds[1][0]);
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // LocalPtIds
  // -----------------------------------------------------------------------
  int numPoints = facesPtIds[3].size() + facesPtIds[1].size() + endInterPtIds.size() + begInterPtIds.size();
  fprintf(stdout,"NUMBER OF POINTS: %d\n", numPoints);
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // Add faces
  // -----------------------------------------------------------------------
  for (int i=0; i<6; i++)
  {
    int fnpts = facesPtIds[i].size();
    if (fnpts > 0)
    {
      vtkIdType *face = new vtkIdType[fnpts];
      for (int j=0; j<fnpts; j++)
        face[j] = facesPtIds[i][j];
      allCells->InsertNextCell(fnpts, face);
      groupIds->InsertNextTuple1(groupId);
      patchIds->InsertNextTuple1(i);
      delete [] face;
    }
  }
  // -----------------------------------------------------------------------

  fprintf(stdout,"\n");
  return SV_OK;
}

// ----------------------
// GetTrifurcationType
// ----------------------
int vtkSVPolycubeGenerator::GetTrifurcationType(vtkSVCenterlineGCell *gCell,
                                                int &type)
{
  int numChildren = gCell->Children.size();

  if (numChildren != 3)
  {
    fprintf(stderr,"Trifurcation wrongly identified, number of children is actually %d\n", numChildren);
  }

  fprintf(stdout,"TRIFURCATION, TELL ME THIS DIR:            %d\n", gCell->BranchDir);
  fprintf(stdout,"TRIFURCATION, TELL ME DIVERGING CHILD DIR: %d\n", gCell->Children[gCell->DivergingChild]->BranchDir);
  fprintf(stdout,"TRIFURCATION, TELL ME ALIGNING CHILD DIR:  %d\n", gCell->Children[gCell->AligningChild]->BranchDir);

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
      double compare = fabs(vtkMath::Dot(gCell->Children[i]->RefDirs[j], gCell->Children[i]->BranchVec));
      fprintf(stdout,"TRIFURCATION, Would like to see compare, dir: %d, dot: %.6f\n", j, compare);
      if (compare > maxDot)
      {
        maxDot = compare;
        maxDir = j;
      }
    }
    fprintf(stdout,"Child %d max dir is %d\n", gCell->Children[i]->GroupId, maxDir);
    double dotCheck = vtkMath::Dot(gCell->Children[i]->RefDirs[maxDir], gCell->Children[i]->BranchVec);
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
    if (gCell->Parent == NULL)
    {
      type = 1;
    }
    else
    {
      if ((gCell->Parent->BranchDir + gCell->BranchDir)%2 == 0)
      {
        if (gCell->Parent->BranchDir == BACK || gCell->Parent->BranchDir == FRONT)
        {
          if ((gCell->BranchDir + gCell->Children[gCell->DivergingChild]->BranchDir)%2 == 0)
            type = 3;
          else
            type = 5;
        }
        else
        {
          if ((gCell->BranchDir + gCell->Children[gCell->DivergingChild]->BranchDir)%2 == 0)
            type = 2;
          else
            type = 4;
        }
      }
      else
      {
        if ((gCell->BranchDir + gCell->Children[gCell->DivergingChild]->BranchDir)%2 == 0)
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
// FormTrifurcation
// ----------------------
int vtkSVPolycubeGenerator::FormTrifurcation(vtkSVCenterlineGCell *gCell,
                                             const double pt0[3], const double pt1[3],
                                             const double pt2[3], const double pt3[3],
                                             const double pt4[3], const double pt5[3],
                                             const double centerPt[3],
                                             const double factor,
                                             double vecs[3][3],
                                             double returnPts[4][3])
{
  vtkMath::Subtract(pt0, pt1, vecs[0]);
  vtkMath::Subtract(pt2, pt3, vecs[1]);
  vtkMath::Subtract(pt4, pt5, vecs[2]);
  for (int j=0; j<3; j++)
    vtkMath::Normalize(vecs[j]);

  double midVec[3];
  vtkMath::Add(vecs[1], vecs[2], midVec);
  vtkMath::Normalize(midVec);
  vtkMath::MultiplyScalar(midVec, -1.0);

  // Get ending bifurcation points
  this->GetBifurcationPoint(centerPt, midVec, vecs[1], vecs[2], factor, returnPts[0]);
  this->GetBifurcationPoint(centerPt, midVec, vecs[2], vecs[1], factor, returnPts[1]);
  this->GetBifurcationPoint(centerPt, vecs[1], vecs[2], midVec, factor, returnPts[2]);

  vtkMath::MultiplyScalar(midVec, -1.0);
  double negVecs[2][3];
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
      negVecs[i][j] = -1.0 * vecs[i+1][j];
  }

  this->GetBifurcationPoint(centerPt, negVecs[0], negVecs[1], midVec, factor, returnPts[3]);

  return SV_OK;
}


// ----------------------
// FormBifurcation
// ----------------------
int vtkSVPolycubeGenerator::FormBifurcation(vtkSVCenterlineGCell *gCell,
                                            const double pt0[3], const double pt1[3],
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
int vtkSVPolycubeGenerator::GetBifurcationPoint(const double startPt[3],
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
int vtkSVPolycubeGenerator::GetSquare(const double startPt[3],
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
int vtkSVPolycubeGenerator::GetWedge(const double pt0[3], const double pt1[3],
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

// ----------------------
// SynchronizePolycubeDivisions
// ----------------------
int vtkSVPolycubeGenerator::SynchronizePolycubeDivisions(std::vector<std::vector<int> > &whl_divs)
{

  int groupId, begType, begSplitType, endType, endSplitType;
  for (int i=0; i<whl_divs.size(); i++)
  {
    groupId = whl_divs[i][0];

    vtkSVCenterlineGCell *gCell = this->CenterlineGraph->GetCellByGroupId(groupId);

    if (gCell->GetBeginningType(begType, begSplitType) != SV_OK)
    {
      return SV_ERROR;
    }
    if (gCell->GetEndType(endType, endSplitType) != SV_OK)
    {
      return SV_ERROR;
    }

    if (begType == S_TET_1 || begType == S_TET_3)
    {
      this->SynchronizeSide(gCell, 1, whl_divs);
    }

    if (begType == S_TET_0 || begType == S_TET_2)
    {
      this->SynchronizeSide(gCell, 2, whl_divs);
    }

  }

  return SV_OK;
}

// ----------------------
// SynchronizeSide
// ----------------------
int vtkSVPolycubeGenerator::SynchronizeSide(vtkSVCenterlineGCell *gCell,
                                            const int dim,
                                            std::vector<std::vector<int> > &whl_divs)
{
  for (int i=0; i<gCell->Children.size(); i++)
  {
    for (int j=0; j<whl_divs.size(); j++)
    {
      if (whl_divs[j][0] != gCell->Children[i]->GroupId)
        continue;

      whl_divs[j][dim] = 2*whl_divs[j][dim]-1;

      this->SynchronizeSide(gCell->Children[i], dim, whl_divs);
    }
  }

  return SV_OK;
}

// ----------------------
// GetPolycubeDivisions
// ----------------------
int vtkSVPolycubeGenerator::GetPolycubeDivisions(vtkSVCenterlineGCell *gCell,
                                                 vtkPolyData *polycubePd,
                                                 int &w_div, int &h_div, int &l_div)
{

  int begType, begSplitType;
  if (gCell->GetBeginningType(begType, begSplitType) != SV_OK)
  {
    return SV_ERROR;
  }

  int endType, endSplitType;
  if (gCell->GetEndType(endType, endSplitType) != SV_OK)
  {
    return SV_ERROR;
  }

  // W div
  w_div = this->PolycubeDivisions;

  if (w_div%2 == 0)
  {
    w_div++;
  }

  if (begType == S_TET_1 || begType == S_TET_3 || endType == S_TET_1 || endType == S_TET_3)
  {
    w_div = 2*w_div-1;
  }

  // H div
  h_div = this->PolycubeDivisions;

  if (h_div%2 == 0)
  {
    h_div++;
  }

  if (begType == S_TET_0 || begType == S_TET_2 || endType == S_TET_0 || endType == S_TET_2)
  {
    h_div = 2*h_div-1;
  }

  // L div
  double f0Pts[2][3];
  vtkIdType f0npts, *f0PtIds;
  polycubePd->GetCellPoints(0, f0npts, f0PtIds);
  polycubePd->GetPoint(f0PtIds[0], f0Pts[0]);
  polycubePd->GetPoint(f0PtIds[1], f0Pts[1]);
  l_div = floor(vtkSVMathUtils::Distance(f0Pts[1], f0Pts[0])/(this->PolycubeUnitLength));

  if (l_div < 10)
  {
    l_div = 10;
  }

  return SV_OK;
}

// ----------------------
// FormParametricHexMesh
// ----------------------
int vtkSVPolycubeGenerator::FormParametricHexMesh(vtkSVCenterlineGCell *gCell, vtkPolyData *polycubePd, vtkStructuredGrid *paraHexMesh,
                                                  int w_div, int h_div, int l_div)
{
  //int nTopPts0, nBotPts0, flatTop0, flatBot0;
  //this->CheckFace(polycubePd, 0, nTopPts0, nBotPts0, flatTop0, flatBot0);
  //fprintf(stdout,"FACE 0\n");
  //fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts0);
  //fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop0);
  //fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts0);
  //fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot0);

  //int nTopPts1, nBotPts1, flatTop1, flatBot1;
  //this->CheckFace(polycubePd, 1, nTopPts1, nBotPts1, flatTop1, flatBot1);
  //fprintf(stdout,"FACE 1\n");
  //fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts1);
  //fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop1);
  //fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts1);
  //fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot1);

  //int nTopPts2, nBotPts2, flatTop2, flatBot2;
  //this->CheckFace(polycubePd, 2, nTopPts2, nBotPts2, flatTop2, flatBot2);
  //fprintf(stdout,"FACE 2\n");
  //fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts2);
  //fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop2);
  //fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts2);
  //fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot2);

  //int nTopPts3, nBotPts3, flatTop3, flatBot3;
  //this->CheckFace(polycubePd, 3, nTopPts3, nBotPts3, flatTop3, flatBot3);
  //fprintf(stdout,"FACE 3\n");
  //fprintf(stdout,"  NUM TOP PTS: %d\n", nTopPts3);
  //fprintf(stdout,"  TOP IS FLAT: %d\n", flatTop3);
  //fprintf(stdout,"  NUM BOT PTS: %d\n", nBotPts3);
  //fprintf(stdout,"  BOT IS FLAT: %d\n", flatBot3);

  //int topHorzWedge = 0, topVertWedge = 0;
  //int topSTet0 = 0, topSTet1 = 0, topSTet2 = 0, topSTet3 = 0;
  //int topCTet0 = 0, topCTet1 = 0, topCTet2 = 0, topCTet3 = 0;

  //if (nTopPts0 == 3)
  //{
  //  if (flatTop0)
  //    topSTet2 = 1;
  //  else
  //  {
  //    if (nTopPts2 == 3)
  //    {
  //      if (flatTop2)
  //        topSTet0 = 1;
  //      else
  //        topHorzWedge = 1;
  //    }
  //    else
  //    {
  //      fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}
  //else if (nTopPts0 == 2)
  //{
  //  if (!flatTop0)
  //  {
  //    if (!flatTop1)
  //      topCTet0 = 1;
  //    else if (!flatTop3)
  //      topCTet3 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 1\n");
  //      return SV_ERROR;
  //    }
  //  }
  //  if (!flatTop2)
  //  {
  //    if (!flatTop1)
  //      topCTet1 = 1;
  //    else if (!flatTop3)
  //      topCTet2 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 2\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}

  //if (nTopPts3 == 3)
  //{
  //  if (flatTop3)
  //    topSTet1 = 1;
  //  else
  //  {
  //    if (nTopPts1 == 3)
  //    {
  //      if (flatTop1)
  //        topSTet3 = 1;
  //      else
  //        topVertWedge = 1;
  //    }
  //    else
  //    {
  //      fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}
  //else if (nTopPts3 == 2)
  //{
  //  if (!flatTop3)
  //  {
  //    if (!flatTop2)
  //      topCTet2 = 1;
  //    else if (!flatTop0)
  //      topCTet3 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 3\n");
  //      return SV_ERROR;
  //    }
  //  }
  //  if (!flatTop1)
  //  {
  //    if (!flatTop2)
  //      topCTet1 = 1;
  //    else if (!flatTop0)
  //      topCTet0 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 4\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}

  //int botHorzWedge = 0, botVertWedge = 0;
  //int botSTet0 = 0, botSTet1 = 0, botSTet2 = 0, botSTet3 = 0;
  //int botCTet0 = 0, botCTet1 = 0, botCTet2 = 0, botCTet3 = 0;

  //if (nBotPts0 == 3)
  //{
  //  if (flatBot0)
  //    botSTet0 = 1;
  //  else
  //  {
  //    if (nBotPts2 == 3)
  //    {
  //      if (flatBot2)
  //        botSTet2 = 1;
  //      else
  //        botHorzWedge = 1;
  //    }
  //    else
  //    {
  //      fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}
  //else if (nBotPts0 == 2)
  //{
  //  if (!flatBot0)
  //  {
  //    if (!flatBot1)
  //      botCTet0 = 1;
  //    else if (!flatBot3)
  //      botCTet3 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 5\n");
  //      return SV_ERROR;
  //    }
  //  }
  //  if (!flatBot2)
  //  {
  //    if (!flatBot1)
  //      botCTet1 = 1;
  //    else if (!flatBot3)
  //      botCTet2 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 6\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}

  //if (nBotPts3 == 3)
  //{
  //  if (flatBot3)
  //    botSTet3 = 1;
  //  else
  //  {
  //    if (nBotPts1 == 3)
  //    {
  //      if (flatBot1)
  //        botSTet1 = 1;
  //      else
  //        botVertWedge = 1;
  //    }
  //    else
  //    {
  //      fprintf(stderr,"OPPOSITE FACES NEED TO HAVE SAME NUMBER OF POINTS!!!\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}
  //else if (nBotPts3 == 2)
  //{
  //  if (!flatBot3)
  //  {
  //    if (!flatBot2)
  //      botCTet2 = 1;
  //    else if (!flatBot0)
  //      botCTet3 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 7\n");
  //      return SV_ERROR;
  //    }
  //  }
  //  if (!flatBot1)
  //  {
  //    if (!flatBot2)
  //      botCTet1 = 1;
  //    else if (!flatBot0)
  //      botCTet0 = 1;
  //    else
  //    {
  //      fprintf(stderr,"THIS SHOULDNT REALLY HAPPEN, 8\n");
  //      return SV_ERROR;
  //    }
  //  }
  //}

  // GetFace 0, right side face
  vtkIdType f0npts, *f0PtIds;
  polycubePd->GetCellPoints(0, f0npts, f0PtIds);

  // GetFace 1, face underneath
  vtkIdType f1npts, *f1PtIds;
  polycubePd->GetCellPoints(1, f1npts, f1PtIds);

  // GetFace 2, left side face
  vtkIdType f2npts, *f2PtIds;
  polycubePd->GetCellPoints(2, f2npts, f2PtIds);

  // GetFace 3, face on top
  vtkIdType f3npts, *f3PtIds;
  polycubePd->GetCellPoints(3, f3npts, f3PtIds);

  if (f0npts != f2npts || f1npts != f3npts)
  {
    fprintf(stderr,"Opposite sides of cube with group id %d must have same number of points\n", gCell->GroupId);
    fprintf(stderr,"FO: %d, F2: %d, F1: %d, F3: %d\n", f0npts, f2npts, f1npts, f3npts);
    return SV_ERROR;
  }

  // Form some what of a parallelepiped
  double f0Pts[4][3], f2Pts[4][3];
  polycubePd->GetPoint(f0PtIds[0], f0Pts[0]);
  polycubePd->GetPoint(f0PtIds[1], f0Pts[1]);
  polycubePd->GetPoint(f0PtIds[2], f0Pts[2]);
  polycubePd->GetPoint(f0PtIds[3], f0Pts[3]);
  polycubePd->GetPoint(f2PtIds[3], f2Pts[0]);
  polycubePd->GetPoint(f2PtIds[2], f2Pts[1]);
  polycubePd->GetPoint(f2PtIds[1], f2Pts[2]);
  polycubePd->GetPoint(f2PtIds[0], f2Pts[3]);

  double f1Pts[4][3], f3Pts[4][3];
  polycubePd->GetPoint(f3PtIds[0], f3Pts[0]);
  polycubePd->GetPoint(f3PtIds[1], f3Pts[1]);
  polycubePd->GetPoint(f3PtIds[2], f3Pts[2]);
  polycubePd->GetPoint(f3PtIds[3], f3Pts[3]);
  polycubePd->GetPoint(f1PtIds[3], f1Pts[0]);
  polycubePd->GetPoint(f1PtIds[2], f1Pts[1]);
  polycubePd->GetPoint(f1PtIds[1], f1Pts[2]);
  polycubePd->GetPoint(f1PtIds[0], f1Pts[3]);

  std::vector<std::string> cubeTypes;
  cubeTypes.push_back("NONE");
  cubeTypes.push_back("VERT_WEDGE");
  cubeTypes.push_back("HORZ_WEDGE");
  cubeTypes.push_back("C_TET_0");
  cubeTypes.push_back("C_TET_1");
  cubeTypes.push_back("C_TET_2");
  cubeTypes.push_back("C_TET_3");
  cubeTypes.push_back("S_TET_0");
  cubeTypes.push_back("S_TET_1");
  cubeTypes.push_back("S_TET_2");
  cubeTypes.push_back("S_TET_3");
  cubeTypes.push_back("NOTHANDLED");

  int begType, begSplitType;
  if (gCell->GetBeginningType(begType, begSplitType) != SV_OK)
  {
    return SV_ERROR;
  }
  fprintf(stdout,"CORRESPONDINGLY, BEG TYPE: %s\n", cubeTypes[begType].c_str());

  int endType, endSplitType;
  if (gCell->GetEndType(endType, endSplitType) != SV_OK)
  {
    return SV_ERROR;
  }
  fprintf(stdout,"CORRESPONDINGLY, END TYPE: %s\n", cubeTypes[endType].c_str());
  fprintf(stdout,"\n");


  // BEG TYPE HORZ_WEDGE, S_TET_0, S_TET_2
  //if (nTopPts0 == 3)
  if (begType == HORZ_WEDGE || begType == S_TET_0 || begType == S_TET_2)
  {
    //fprintf(stdout,"N TOP PTS 0: %d\n", nTopPts0);
    polycubePd->GetPoint(f0PtIds[3], f0Pts[2]);
    polycubePd->GetPoint(f0PtIds[4], f0Pts[3]);
  }

  // BEG TYPE HORZ_WEDGE, S_TET_0
  //if (nTopPts2 == 3)
  if (begType == HORZ_WEDGE || begType == S_TET_0 || begType == S_TET_2)
  {
    //fprintf(stdout,"N TOP PTS 2: %d\n", nTopPts2);
    polycubePd->GetPoint(f2PtIds[4], f2Pts[0]);
    polycubePd->GetPoint(f2PtIds[3], f2Pts[1]);
  }

  // BEG TYPE VERT_WEDGE, S_TET_3
  //if (nTopPts1 == 3)
  if (begType == VERT_WEDGE || begType == S_TET_3 || begType == S_TET_1)
  {
    //fprintf(stdout,"N TOP PTS 1: %d\n", nTopPts1);
    polycubePd->GetPoint(f1PtIds[4], f1Pts[0]);
    polycubePd->GetPoint(f1PtIds[3], f1Pts[1]);
  }

  // BEG TYPE VERT_WEDGE, S_TET_3, S_TET_1
  //if (nTopPts3 == 3)
  if (begType == VERT_WEDGE || begType == S_TET_3 || begType == S_TET_1)
  {
    //fprintf(stdout,"N TOP PTS 3: %d\n", nTopPts3);
    polycubePd->GetPoint(f3PtIds[3], f3Pts[2]);
    polycubePd->GetPoint(f3PtIds[4], f3Pts[3]);
  }

  fprintf(stdout,"WHAT ARE POINTS\n");
  fprintf(stdout,"FACE 0: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f0Pts[0][0], f0Pts[0][1], f0Pts[0][2],
                                                                                             f0Pts[1][0], f0Pts[1][1], f0Pts[1][2],
                                                                                             f0Pts[2][0], f0Pts[2][1], f0Pts[2][2],
                                                                                             f0Pts[3][0], f0Pts[3][1], f0Pts[3][2]);
  fprintf(stdout,"FACE 1: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f1Pts[0][0], f1Pts[0][1], f1Pts[0][2],
                                                                                             f1Pts[1][0], f1Pts[1][1], f1Pts[1][2],
                                                                                             f1Pts[2][0], f1Pts[2][1], f1Pts[2][2],
                                                                                             f1Pts[3][0], f1Pts[3][1], f1Pts[3][2]);
  fprintf(stdout,"FACE 2: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f2Pts[0][0], f2Pts[0][1], f2Pts[0][2],
                                                                                             f2Pts[1][0], f2Pts[1][1], f2Pts[1][2],
                                                                                             f2Pts[2][0], f2Pts[2][1], f2Pts[2][2],
                                                                                             f2Pts[3][0], f2Pts[3][1], f2Pts[3][2]);
  fprintf(stdout,"FACE 3: %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f, %.6f %.6f %.6f\n", f3Pts[0][0], f3Pts[0][1], f3Pts[0][2],
                                                                                             f3Pts[1][0], f3Pts[1][1], f3Pts[1][2],
                                                                                             f3Pts[2][0], f3Pts[2][1], f3Pts[2][2],
                                                                                             f3Pts[3][0], f3Pts[3][1], f3Pts[3][2]);

  //fprintf(stdout,"THIS SAYS TOP VERT WEDGE: %d\n", topVertWedge);
  //fprintf(stdout,"THIS SAYS BOT VERT WEDGE: %d\n", botVertWedge);
  //fprintf(stdout,"THIS SAYS TOP HORZ WEDGE: %d\n", topHorzWedge);
  //fprintf(stdout,"THIS SAYS BOT HORZ WEDGE: %d\n", botHorzWedge);
  //fprintf(stdout,"THIS SAYS TOP SIDE TET 0: %d\n", topSTet0);
  //fprintf(stdout,"THIS SAYS TOP SIDE TET 1: %d\n", topSTet1);
  //fprintf(stdout,"THIS SAYS TOP SIDE TET 2: %d\n", topSTet2);
  //fprintf(stdout,"THIS SAYS TOP SIDE TET 3: %d\n", topSTet3);
  //fprintf(stdout,"THIS SAYS TOP CORN TET 0: %d\n", topCTet0);
  //fprintf(stdout,"THIS SAYS TOP CORN TET 1: %d\n", topCTet1);
  //fprintf(stdout,"THIS SAYS TOP CORN TET 2: %d\n", topCTet2);
  //fprintf(stdout,"THIS SAYS TOP CORN TET 3: %d\n", topCTet3);
  //fprintf(stdout,"THIS SAYS BOT SIDE TET 0: %d\n", botSTet0);
  //fprintf(stdout,"THIS SAYS BOT SIDE TET 1: %d\n", botSTet1);
  //fprintf(stdout,"THIS SAYS BOT SIDE TET 2: %d\n", botSTet2);
  //fprintf(stdout,"THIS SAYS BOT SIDE TET 3: %d\n", botSTet3);
  //fprintf(stdout,"THIS SAYS BOT CORN TET 0: %d\n", botCTet0);
  //fprintf(stdout,"THIS SAYS BOT CORN TET 1: %d\n", botCTet1);
  //fprintf(stdout,"THIS SAYS BOT CORN TET 2: %d\n", botCTet2);
  //fprintf(stdout,"THIS SAYS BOT CORN TET 3: %d\n", botCTet3);
  //fprintf(stdout,"\n");

  // Sides of cube
  double face0Vec0[3], face0Vec1[3];
  vtkMath::Subtract(f0Pts[3], f0Pts[0], face0Vec0);
  vtkMath::Normalize(face0Vec0);
  double face0Dist0 = vtkSVMathUtils::Distance(f0Pts[3], f0Pts[0]);

  vtkMath::Subtract(f0Pts[2], f0Pts[1], face0Vec1);
  vtkMath::Normalize(face0Vec1);
  double face0Dist1 = vtkSVMathUtils::Distance(f0Pts[2], f0Pts[1]);

  double face2Vec0[3], face2Vec1[3];
  vtkMath::Subtract(f2Pts[3], f2Pts[0], face2Vec0);
  vtkMath::Normalize(face2Vec0);
  double face2Dist0 = vtkSVMathUtils::Distance(f2Pts[3], f2Pts[0]);

  vtkMath::Subtract(f2Pts[2], f2Pts[1], face2Vec1);
  vtkMath::Normalize(face2Vec1);
  double face2Dist1 = vtkSVMathUtils::Distance(f2Pts[2], f2Pts[1]);

  // Top and bottom of cube
  double face5Vec0[3], face5Vec1[3];
  vtkMath::Subtract(f1Pts[0], f1Pts[3], face5Vec0);
  vtkMath::Normalize(face5Vec0);
  double face5Dist0 = vtkSVMathUtils::Distance(f1Pts[0], f1Pts[3]);

  vtkMath::Subtract(f3Pts[0], f3Pts[3], face5Vec1);
  vtkMath::Normalize(face5Vec1);
  double face5Dist1 = vtkSVMathUtils::Distance(f3Pts[0], f3Pts[3]);

  double face4Vec0[3], face4Vec1[3];
  vtkMath::Subtract(f1Pts[1], f1Pts[2], face4Vec0);
  vtkMath::Normalize(face4Vec0);
  double face4Dist0 = vtkSVMathUtils::Distance(f1Pts[1], f1Pts[2]);

  vtkMath::Subtract(f3Pts[1], f3Pts[2], face4Vec1);
  vtkMath::Normalize(face4Vec1);
  double face4Dist1 = vtkSVMathUtils::Distance(f3Pts[1], f3Pts[2]);

  vtkNew(vtkPoints, f5GridPts);
  f5GridPts->SetNumberOfPoints(w_div*h_div);
  vtkNew(vtkPoints, f4GridPts);
  f4GridPts->SetNumberOfPoints(w_div*h_div);

  int dim2D_1[3]; dim2D_1[0] = w_div; dim2D_1[1] = h_div; dim2D_1[2] = 1;

  for (int i=0; i<w_div; i++)
  {
    double x5Vec0[3], x5Vec1[3], x4Vec0[3], x4Vec1[3];

    for (int j=0; j<3; j++)
    {
      x5Vec0[j] = face5Vec0[j]*i*(face5Dist0/(w_div-1));
      x5Vec1[j] = face5Vec1[j]*i*(face5Dist1/(w_div-1));

      x4Vec0[j] = face4Vec0[j]*i*(face4Dist0/(w_div-1));
      x4Vec1[j] = face4Vec1[j]*i*(face4Dist1/(w_div-1));
    }

    double f5Start[3], f5End[3], f4Start[3], f4End[3];
    vtkMath::Add(f1Pts[3], x5Vec0, f5Start);
    vtkMath::Add(f3Pts[3], x5Vec1, f5End);

    vtkMath::Add(f1Pts[2], x4Vec0, f4Start);
    vtkMath::Add(f3Pts[2], x4Vec1, f4End);

    double f5HVec[3], f4HVec[3];
    vtkMath::Subtract(f5End, f5Start, f5HVec);
    vtkMath::Normalize(f5HVec);
    double f5HDist = vtkSVMathUtils::Distance(f5End, f5Start);

    vtkMath::Subtract(f4End, f4Start, f4HVec);
    vtkMath::Normalize(f4HVec);
    double f4HDist = vtkSVMathUtils::Distance(f4End, f4Start);

    // Another tet top face
    double face4DiagVecs[2][3], x4DiagVecs[2][3], x4DiagPts[2][3];
    double x4ToDiag[3], x4FromDiag[3], x4AcrossDiag[3];
    double face4DiagDists[2], x4ToDiagDist, x4FromDiagDist, x4AcrossDiagDist;
    //if (topSTet1)
    if (begType == S_TET_1)
    {
      double midPt[3];
      polycubePd->GetPoint(f1PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(midPt, f1Pts[2], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(midPt, f1Pts[2]);

      vtkMath::Subtract(f1Pts[1], midPt, startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(midPt, f1Pts[1]);

      vtkMath::Subtract(midPt, f3Pts[2], face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f3Pts[2]);

      vtkMath::Subtract(f3Pts[1], midPt, face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f3Pts[1]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(h_div-1));
        startPtVecs[1][j] = startVecs[1][j]*((i+1)%h_div)*(startVecDists[1]/(h_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(h_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*((i+1)%h_div)*(face4DiagDists[1]/(h_div-1));
      }

      if (i <= h_div-1)
      {
        vtkMath::Add(f1Pts[2], startPtVecs[0], f4Start);
        vtkMath::Add(f3Pts[2], x4DiagVecs[0], x4DiagPts[0]);

        vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[0], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[0]);
      }
      else
      {
        vtkMath::Add(midPt, startPtVecs[1], f4Start);
        vtkMath::Add(midPt, x4DiagVecs[1], x4DiagPts[1]);

        vtkMath::Subtract(x4DiagPts[1], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);
      }
    }

    //if (topSTet3)
    if (begType == S_TET_3)
    {
      double midPt[3];
      polycubePd->GetPoint(f3PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(midPt, f3Pts[2], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(midPt, f3Pts[2]);

      vtkMath::Subtract(f3Pts[1], midPt, startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(midPt, f3Pts[1]);

      vtkMath::Subtract(midPt, f1Pts[2], face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f1Pts[2]);

      vtkMath::Subtract(f1Pts[1], midPt, face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f1Pts[1]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(h_div-1));
        startPtVecs[1][j] = startVecs[1][j]*((i+1)%h_div)*(startVecDists[1]/(h_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(h_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*((i+1)%h_div)*(face4DiagDists[1]/(h_div-1));
      }


      if (i <= h_div-1)
      {
        vtkMath::Add(f3Pts[2], startPtVecs[0], f4End);
        vtkMath::Add(f1Pts[2], x4DiagVecs[0], x4DiagPts[0]);

        vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[0], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[0]);
      }
      else
      {
        vtkMath::Add(midPt, startPtVecs[1], f4End);
        vtkMath::Add(midPt, x4DiagVecs[1], x4DiagPts[1]);

        vtkMath::Subtract(x4DiagPts[1], f4Start, x4ToDiag);
        vtkMath::Normalize(x4ToDiag);
        x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], f4Start);

        vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
        vtkMath::Normalize(x4FromDiag);
        x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);
      }
    }

    //if (topSTet2)
    if (begType == S_TET_2)
    {
      double midPt[3];
      polycubePd->GetPoint(f2PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(f0Pts[1], f2Pts[1], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(f0Pts[1], f2Pts[1]);

      vtkMath::Subtract(f0Pts[2], f2Pts[2], startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(f0Pts[2], f2Pts[2]);

      vtkMath::Subtract(f0Pts[1], midPt, face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f0Pts[1]);

      vtkMath::Subtract(f0Pts[2], midPt, face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f0Pts[2]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(w_div-1));
        startPtVecs[1][j] = startVecs[1][j]*i*(startVecDists[1]/(w_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(w_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*i*(face4DiagDists[1]/(w_div-1));
      }

      vtkMath::Add(f2Pts[1], startPtVecs[0], f4Start);
      vtkMath::Add(f2Pts[2], startPtVecs[1], f4End);
      vtkMath::Add(midPt, x4DiagVecs[0], x4DiagPts[0]);
      vtkMath::Add(midPt, x4DiagVecs[1], x4DiagPts[1]);

      vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
      vtkMath::Normalize(x4ToDiag);
      x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

      vtkMath::Subtract(x4DiagPts[1], x4DiagPts[0], x4AcrossDiag);
      vtkMath::Normalize(x4AcrossDiag);
      x4AcrossDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], x4DiagPts[0]);

      vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
      vtkMath::Normalize(x4FromDiag);
      x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);

    }

    //if (topSTet0)
    if (begType == S_TET_0)
    {
      double midPt[3];
      polycubePd->GetPoint(f0PtIds[2], midPt);

      double startVecs[2][3], startPtVecs[2][3], startVecDists[2];

      vtkMath::Subtract(f0Pts[1], f2Pts[1], startVecs[0]);
      vtkMath::Normalize(startVecs[0]);
      startVecDists[0] = vtkSVMathUtils::Distance(f0Pts[1], f2Pts[1]);

      vtkMath::Subtract(f0Pts[2], f2Pts[2], startVecs[1]);
      vtkMath::Normalize(startVecs[1]);
      startVecDists[1] = vtkSVMathUtils::Distance(f0Pts[2], f2Pts[2]);

      vtkMath::Subtract(midPt, f2Pts[1], face4DiagVecs[0]);
      vtkMath::Normalize(face4DiagVecs[0]);
      face4DiagDists[0] = vtkSVMathUtils::Distance(midPt, f2Pts[1]);

      vtkMath::Subtract(midPt, f2Pts[2], face4DiagVecs[1]);
      vtkMath::Normalize(face4DiagVecs[1]);
      face4DiagDists[1] = vtkSVMathUtils::Distance(midPt, f2Pts[2]);

      for (int j=0; j<3; j++)
      {
        startPtVecs[0][j] = startVecs[0][j]*i*(startVecDists[0]/(w_div-1));
        startPtVecs[1][j] = startVecs[1][j]*i*(startVecDists[1]/(w_div-1));
        x4DiagVecs[0][j] = face4DiagVecs[0][j]*i*(face4DiagDists[0]/(w_div-1));
        x4DiagVecs[1][j] = face4DiagVecs[1][j]*i*(face4DiagDists[1]/(w_div-1));
      }


      vtkMath::Add(f2Pts[1], startPtVecs[0], f4Start);
      vtkMath::Add(f2Pts[2], startPtVecs[1], f4End);
      vtkMath::Add(f2Pts[1], x4DiagVecs[0], x4DiagPts[0]);
      vtkMath::Add(f2Pts[2], x4DiagVecs[1], x4DiagPts[1]);

      vtkMath::Subtract(x4DiagPts[0], f4Start, x4ToDiag);
      vtkMath::Normalize(x4ToDiag);
      x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPts[0], f4Start);

      vtkMath::Subtract(x4DiagPts[1], x4DiagPts[0], x4AcrossDiag);
      vtkMath::Normalize(x4AcrossDiag);
      x4AcrossDiagDist = vtkSVMathUtils::Distance(x4DiagPts[1], x4DiagPts[0]);

      vtkMath::Subtract(f4End, x4DiagPts[1], x4FromDiag);
      vtkMath::Normalize(x4FromDiag);
      x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPts[1]);
    }

    // Tet top face
    double face4DiagVec[3], x4DiagVec[3], x4DiagPoint[3];
    double face4DiagDist;
    //if (topCTet1 || topCTet3)
    if (begType == C_TET_3 || begType == C_TET_1)
    {
      vtkMath::Subtract(f3Pts[1], f1Pts[2], face4DiagVec);
      vtkMath::Normalize(face4DiagVec);
      face4DiagDist = vtkSVMathUtils::Distance(f3Pts[1], f1Pts[2]);

      for (int j=0; j<3; j++)
        x4DiagVec[j] = face4DiagVec[j]*i*(face4DiagDist/(w_div-1));

      vtkMath::Add(f1Pts[2], x4DiagVec, x4DiagPoint);
    }
    //if (topCTet0 || topCTet2)
    if (begType == C_TET_2 || begType == C_TET_0)
    {
      vtkMath::Subtract(f1Pts[1], f3Pts[2], face4DiagVec);
      vtkMath::Normalize(face4DiagVec);
      face4DiagDist = vtkSVMathUtils::Distance(f1Pts[1], f3Pts[2]);

      for (int j=0; j<3; j++)
        x4DiagVec[j] = face4DiagVec[j]*i*(face4DiagDist/(w_div-1));

      vtkMath::Add(f3Pts[2], x4DiagVec, x4DiagPoint);
    }

    //if (topCTet0 || topCTet1 || topCTet2 || topCTet3)
    if (begType == C_TET_2 || begType == C_TET_3 || begType == C_TET_0 || begType == C_TET_1)
    {
      vtkMath::Subtract(x4DiagPoint, f4Start, x4ToDiag);
      vtkMath::Normalize(x4ToDiag);
      x4ToDiagDist = vtkSVMathUtils::Distance(x4DiagPoint, f4Start);

      vtkMath::Subtract(f4End, x4DiagPoint, x4FromDiag);
      vtkMath::Normalize(x4FromDiag);
      x4FromDiagDist = vtkSVMathUtils::Distance(f4End, x4DiagPoint);
    }

    // Tet bot face
    double face5DiagVec[3], x5DiagVec[3], x5DiagPoint[3], x5ToDiag[3], x5FromDiag[3];
    double face5DiagDist, x5ToDiagDist, x5FromDiagDist;

    //if (botCTet1 || botCTet3)
    if (endType == C_TET_3 || endType == C_TET_1)
    {
      vtkMath::Subtract(f3Pts[0], f1Pts[3], face5DiagVec);
      vtkMath::Normalize(face5DiagVec);
      face5DiagDist = vtkSVMathUtils::Distance(f3Pts[0], f1Pts[3]);

      for (int j=0; j<3; j++)
        x5DiagVec[j] = face5DiagVec[j]*i*(face5DiagDist/(w_div-1));

      vtkMath::Add(f1Pts[3], x5DiagVec, x5DiagPoint);
    }
    //if (botCTet0 || botCTet2)
    if (endType == C_TET_2 || endType == C_TET_0)
    {
      vtkMath::Subtract(f1Pts[0], f3Pts[3], face5DiagVec);
      vtkMath::Normalize(face5DiagVec);
      face5DiagDist = vtkSVMathUtils::Distance(f1Pts[0], f3Pts[3]);

      for (int j=0; j<3; j++)
        x5DiagVec[j] = face5DiagVec[j]*i*(face5DiagDist/(w_div-1));

      vtkMath::Add(f3Pts[3], x5DiagVec, x5DiagPoint);

    }

    //if (botCTet0 || botCTet1 || botCTet2 || botCTet3)
    if (endType == C_TET_2 || endType == C_TET_3 || endType == C_TET_0 || endType == C_TET_1)
    {
      vtkMath::Subtract(x5DiagPoint, f5Start, x5ToDiag);
      vtkMath::Normalize(x5ToDiag);
      x5ToDiagDist = vtkSVMathUtils::Distance(x5DiagPoint, f5Start);

      vtkMath::Subtract(f5End, x5DiagPoint, x5FromDiag);
      vtkMath::Normalize(x5FromDiag);
      x5FromDiagDist = vtkSVMathUtils::Distance(f5End, x5DiagPoint);
    }

    for (int j=0; j<h_div; j++)
    {
      double z5Vec[3], z4Vec[3];
      for (int k=0; k<3; k++)
      {
        z5Vec[k] = f5HVec[k]*j*(f5HDist/(h_div-1));
        z4Vec[k] = f4HVec[k]*j*(f4HDist/(h_div-1));
      }

      double new5Pt[3], new4Pt[3];
      vtkMath::Add(f5Start, z5Vec, new5Pt);
      vtkMath::Add(f4Start, z4Vec, new4Pt);

      //if (topCTet1 || topCTet3)
      if (begType == C_TET_3 || begType == C_TET_1)
      {
        if (j <= i)
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else
        {
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(j-i)*(x4FromDiagDist/((w_div-1)-i));
          vtkMath::Add(x4DiagPoint, z4Vec, new4Pt);
        }
      }
      //if (topCTet0 || topCTet2)
      if (begType == C_TET_2 || begType == C_TET_0)
      {
        if (j <= (w_div-1)-i)
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(j-((w_div-1)-i))*(x4FromDiagDist/diag_div);
          vtkMath::Add(x4DiagPoint, z4Vec, new4Pt);
        }
      }

      //if (botCTet1 || botCTet3)
      if (endType == C_TET_3 || endType == C_TET_1)
      {
        if (j <= i)
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z5Vec[k] = x5ToDiag[k]*j*(x5ToDiagDist/diag_div);
          vtkMath::Add(f5Start, z5Vec, new5Pt);
        }
        else
        {
          for (int k=0; k<3; k++)
            z5Vec[k] = x5FromDiag[k]*(j-i)*(x5FromDiagDist/((w_div-1)-i));
          vtkMath::Add(x5DiagPoint, z5Vec, new5Pt);
        }
      }
      //if (botCTet0 || botCTet2)
      if (endType == C_TET_2 || endType == C_TET_0)
      {
        if (j <= (w_div-1)-i)
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z5Vec[k] = x5ToDiag[k]*j*(x5ToDiagDist/diag_div);
          vtkMath::Add(f5Start, z5Vec, new5Pt);
        }
        else
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z5Vec[k] = x5FromDiag[k]*(j-((w_div-1)-i))*(x5FromDiagDist/diag_div);
          vtkMath::Add(x5DiagPoint, z5Vec, new5Pt);
        }
      }

      //if (topSTet1)
      if (begType == S_TET_1)
      {
        if (i <= h_div-1)
        {
          if (j < h_div-i)
          {
            int diag_div = ((h_div-1)-i);
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-((h_div-1)-i))*(x4FromDiagDist/i);
            vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
          }
        }
        else
        {
          if (j < (i+1)%h_div)
          {
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/((i+1)%h_div));
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            int diag_div = ((h_div-1)-((i+1)%h_div));
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-((i+1)%h_div))*(x4FromDiagDist/diag_div);
            vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
          }
        }
      }

      //if (topSTet3)
      if (begType == S_TET_3)
      {
        if (i <= h_div-1)
        {
          if (j <= i)
          {
            int diag_div = i;
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-i)*(x4FromDiagDist/((h_div-1)-i));
            vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
          }
        }
        else
        {
          if (j <= (h_div-1) - ((i+1)%h_div))
          {
            int diag_div = (h_div-1) - ((i+1)%h_div);
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/(diag_div));
            vtkMath::Add(f4Start, z4Vec, new4Pt);
          }
          else
          {
            int diag_div = (i+1)%h_div;
            if (diag_div == 0)
              diag_div = 1;
            for (int k=0; k<3; k++)
              z4Vec[k] = x4FromDiag[k]*(j-((h_div-1)-(i+1)%h_div))*(x4FromDiagDist/diag_div);
            vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
          }
        }
      }

      //if (topSTet2)
      if (begType == S_TET_2)
      {
        if (j <= (w_div-1)-i)
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else if (j > (w_div-1)-i && j <= ((w_div-1)-i)+2*i)
        {
          int diag_div = 2*i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4AcrossDiag[k]*(j-((w_div-1)-i))*(x4AcrossDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
        }
        else
        {
          int diag_div = (w_div-1)-i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(((j+1)%w_div)-i)*(x4FromDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
        }
      }

      //if (topSTet0)
      if (begType == S_TET_0)
      {
        if (j <= i)
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4ToDiag[k]*j*(x4ToDiagDist/diag_div);
          vtkMath::Add(f4Start, z4Vec, new4Pt);
        }
        else if (j > i && j <= (2*(w_div-1) - i))
        {
          int diag_div = (2*(w_div-1)-(2*i));
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4AcrossDiag[k]*(j-i)*(x4AcrossDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[0], z4Vec, new4Pt);
        }
        else
        {
          int diag_div = i;
          if (diag_div == 0)
            diag_div = 1;
          for (int k=0; k<3; k++)
            z4Vec[k] = x4FromDiag[k]*(j-(2*(w_div-1)-i))*(x4FromDiagDist/diag_div);
          vtkMath::Add(x4DiagPts[1], z4Vec, new4Pt);
        }
      }

      int pos[3]; pos[0]= i; pos[1] = j; pos[2] = 0;
      int pId = vtkStructuredData::ComputePointId(dim2D_1, pos);

      f5GridPts->SetPoint(pId, new5Pt);
      f4GridPts->SetPoint(pId, new4Pt);
    }
  }

  vtkNew(vtkPolyData, grid5Poly);
  grid5Poly->SetPoints(f5GridPts);
  vtkNew(vtkPolyData, grid4Poly);
  grid4Poly->SetPoints(f4GridPts);

  vtkNew(vtkPoints, paraPoints);
  int dim[3]; dim[0] = w_div; dim[1] = h_div; dim[2] = l_div;
  paraHexMesh->SetDimensions(dim);
  paraHexMesh->SetPoints(paraPoints);
  paraHexMesh->GetPoints()->SetNumberOfPoints(w_div*h_div*l_div);

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      int pos2D_1[3]; pos2D_1[0] = i; pos2D_1[1] = j; pos2D_1[2] = 0;
      int getId1 = vtkStructuredData::ComputePointId(dim2D_1, pos2D_1);

      double startPt1[3], endPt1[3];
      f5GridPts->GetPoint(getId1, startPt1);
      f4GridPts->GetPoint(getId1, endPt1);

      double hLVec[3];
      vtkMath::Subtract(startPt1, endPt1, hLVec);
      vtkMath::Normalize(hLVec);
      double hLDist = vtkSVMathUtils::Distance(endPt1, startPt1);

      for (int k=0; k<l_div; k++)
      {
        double yVec[3];
        for (int l=0; l<3; l++)
        {
          yVec[l] = hLVec[l]*k*(hLDist/(1 - l_div));
        }

        double newPt[3];
        vtkMath::Add(startPt1, yVec, newPt);

        int pos[3]; pos[0]= i; pos[1] = j; pos[2] = k;
        int pId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(pId, newPt);
      }
    }
  }

  double midPt0[3], midPt1[3];
  //if (topHorzWedge)
  if (begType == HORZ_WEDGE)
  {
    polycubePd->GetPoint(f2PtIds[2], midPt0);
    polycubePd->GetPoint(f0PtIds[2], midPt1);
    this->PushStructuredGridXAxis(paraHexMesh, midPt0, midPt1, 0);
  }

  //if (botHorzWedge)
  if (endType == HORZ_WEDGE)
  {
    //if (topHorzWedge)
    if (begType == HORZ_WEDGE || begType == S_TET_0 || begType == S_TET_2)
    {
    polycubePd->GetPoint(f2PtIds[5], midPt0);
    polycubePd->GetPoint(f0PtIds[5], midPt1);
    }
    else
    {
    polycubePd->GetPoint(f2PtIds[4], midPt0);
    polycubePd->GetPoint(f0PtIds[4], midPt1);
    }
    this->PushStructuredGridXAxis(paraHexMesh, midPt0, midPt1, 1);
  }

  //if (topVertWedge)
  if (begType == VERT_WEDGE)
  {
    polycubePd->GetPoint(f1PtIds[2], midPt0);
    polycubePd->GetPoint(f3PtIds[2], midPt1);
    this->PushStructuredGridZAxis(paraHexMesh, midPt0, midPt1, 0);
  }

  if (endType == VERT_WEDGE)
  //if (botVertWedge)
  {
    if (begType == VERT_WEDGE || begType == S_TET_3 || begType == S_TET_1)
    //if (topVertWedge)
    {
      polycubePd->GetPoint(f1PtIds[5], midPt0);
      polycubePd->GetPoint(f3PtIds[5], midPt1);
    }
    else
    {
      polycubePd->GetPoint(f1PtIds[4], midPt0);
      polycubePd->GetPoint(f3PtIds[4], midPt1);
    }
    this->PushStructuredGridZAxis(paraHexMesh, midPt0, midPt1, 1);
  }

  return SV_OK;
}

// ----------------------
// PushStructuredGridZAxis
// ----------------------
int vtkSVPolycubeGenerator::PushStructuredGridZAxis(vtkStructuredGrid *paraHexMesh,
                                                  const double midPt0[3],
                                                  const double midPt1[3],
                                                  const int isBottom)
{
  int dim[3];
  paraHexMesh->GetDimensions(dim);

  int h_div = dim[1];
  double h = vtkSVMathUtils::Distance(midPt1, midPt0);
  double h_dist = h/(h_div-1);

  double h_vec[3];
  vtkMath::Subtract(midPt1, midPt0, h_vec);
  vtkMath::Normalize(h_vec);

  int w_div = dim[0];
  int half_w_div = floor(w_div/2.0);

  // Set new bottom extend points in the middle
  int pos[3];
  pos[0]= half_w_div;

  if (isBottom)
    pos[2] = 0;
  else
    pos[2] = dim[2]-1;

  for (int i=0; i<h_div; i++)
  {
    pos[1] = i;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);

    double z_vec[3];
    for (int j=0; j<3; j++)
      z_vec[j] = h_vec[j]*i*h_dist;

    double newPt[3];
    vtkMath::Add(midPt0, z_vec, newPt);

    // New point in middle
    paraHexMesh->GetPoints()->SetPoint(ptId, newPt);

    // New point on 0 side
    double firstPt[3];
    //vtkMath::Add(pt0, z_vec, firstPt);

    pos[0] = 0;
    int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(firstPtId, firstPt);
    //paraHexMesh->GetPoints()->SetPoint(firstPtId, firstPt);

    double first_w = vtkSVMathUtils::Distance(newPt, firstPt);
    double first_w_dist = first_w/(half_w_div);

    double first_w_vec[3];
    vtkMath::Subtract(newPt, firstPt, first_w_vec);
    vtkMath::Normalize(first_w_vec);

    for (int j=0; j<half_w_div; j++)
    {
      pos[0] = j;

      double x_vec[3];
      for (int k=0; k<3; k++)
        x_vec[k] = first_w_vec[k]*j*first_w_dist;

      double inPt[3];
      vtkMath::Add(firstPt, x_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }

    //New point on end side
    double lastPt[3];
    //vtkMath::Add(pt1, z_vec, lastPt);

    pos[0] = dim[0]-1;
    int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(lastPtId, lastPt);
    //paraHexMesh->GetPoints()->SetPoint(lastPtId, lastPt);

    double last_w = vtkSVMathUtils::Distance(lastPt, newPt);
    double last_w_dist = last_w/(half_w_div);

    double last_w_vec[3];
    vtkMath::Subtract(lastPt, newPt, last_w_vec);
    vtkMath::Normalize(last_w_vec);

    for (int j=0; j<half_w_div; j++)
    {
      pos[0] = half_w_div+j;

      double x_vec[3];
      for (int k=0; k<3; k++)
        x_vec[k] = last_w_vec[k]*j*last_w_dist;

      double inPt[3];
      vtkMath::Add(newPt, x_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }
  }

  int l_div = dim[2];

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      pos[0] = i; pos[1] = j;

      pos[2] = 0;
      int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

      double firstPt[3];
      paraHexMesh->GetPoint(firstPtId, firstPt);

      pos[2] = dim[2]-1;
      int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

      double lastPt[3];
      paraHexMesh->GetPoint(lastPtId, lastPt);

      double l = vtkSVMathUtils::Distance(lastPt, firstPt);
      double l_dist = l/(l_div-1);

      double l_vec[3];
      vtkMath::Subtract(lastPt, firstPt, l_vec);
      vtkMath::Normalize(l_vec);

      for (int k=0; k<l_div; k++)
      {
        pos[2] = k;

        double y_vec[3];
        for (int l=0; l<3; l++)
          y_vec[l] = l_vec[l]*k*l_dist;

        double inPt[3];
        vtkMath::Add(firstPt, y_vec, inPt);

        int newPtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// PushStructuredGridXAxis
// ----------------------
int vtkSVPolycubeGenerator::PushStructuredGridXAxis(vtkStructuredGrid *paraHexMesh,
                                                    const double midPt0[3],
                                                    const double midPt1[3],
                                                    const int isBottom)
{
  int dim[3];
  paraHexMesh->GetDimensions(dim);

  int w_div = dim[0];
  double w = vtkSVMathUtils::Distance(midPt1, midPt0);
  double w_dist = w/(w_div-1);

  double w_vec[3];
  vtkMath::Subtract(midPt1, midPt0, w_vec);
  vtkMath::Normalize(w_vec);

  int h_div = dim[1];
  int half_h_div = floor(h_div/2.0);

  // Set new bottom extend points in the middle
  int pos[3];
  pos[1]= half_h_div;

  if (isBottom)
    pos[2] = 0;
  else
    pos[2] = dim[2]-1;

  for (int i=0; i<w_div; i++)
  {
    pos[0] = i;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);

    double x_vec[3];
    for (int j=0; j<3; j++)
      x_vec[j] = w_vec[j]*i*w_dist;

    double newPt[3];
    vtkMath::Add(midPt0, x_vec, newPt);

    // New point in middle
    paraHexMesh->GetPoints()->SetPoint(ptId, newPt);

    // New point on 0 side
    double firstPt[3];
    //vtkMath::Add(pt0, x_vec, firstPt);

    pos[1] = 0;
    int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(firstPtId, firstPt);
    //paraHexMesh->GetPoints()->SetPoint(firstPtId, firstPt);

    double first_h = vtkSVMathUtils::Distance(newPt, firstPt);
    double first_h_dist = first_h/(half_h_div);

    double first_h_vec[3];
    vtkMath::Subtract(newPt, firstPt, first_h_vec);
    vtkMath::Normalize(first_h_vec);

    for (int j=0; j<half_h_div; j++)
    {
      pos[1] = j;

      double z_vec[3];
      for (int k=0; k<3; k++)
        z_vec[k] = first_h_vec[k]*j*first_h_dist;

      double inPt[3];
      vtkMath::Add(firstPt, z_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }

    //New point on end side
    double lastPt[3];
    //vtkMath::Add(pt1, x_vec, lastPt);

    pos[1] = dim[1]-1;
    int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

    paraHexMesh->GetPoints()->GetPoint(lastPtId, lastPt);
    //paraHexMesh->GetPoints()->SetPoint(lastPtId, lastPt);

    double last_h = vtkSVMathUtils::Distance(lastPt, newPt);
    double last_h_dist = last_h/(half_h_div);

    double last_h_vec[3];
    vtkMath::Subtract(lastPt, newPt, last_h_vec);
    vtkMath::Normalize(last_h_vec);

    for (int j=0; j<half_h_div; j++)
    {
      pos[1] = half_h_div+j;

      double z_vec[3];
      for (int k=0; k<3; k++)
        z_vec[k] = last_h_vec[k]*j*last_h_dist;

      double inPt[3];
      vtkMath::Add(newPt, z_vec, inPt);

      int newPtId = vtkStructuredData::ComputePointId(dim, pos);
      paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
    }
  }

  int l_div = dim[2];

  for (int i=0; i<w_div; i++)
  {
    for (int j=0; j<h_div; j++)
    {
      pos[0] = i; pos[1] = j;

      pos[2] = 0;
      int firstPtId = vtkStructuredData::ComputePointId(dim, pos);

      double firstPt[3];
      paraHexMesh->GetPoint(firstPtId, firstPt);

      pos[2] = dim[2]-1;
      int lastPtId = vtkStructuredData::ComputePointId(dim, pos);

      double lastPt[3];
      paraHexMesh->GetPoint(lastPtId, lastPt);

      double l = vtkSVMathUtils::Distance(lastPt, firstPt);
      double l_dist = l/(l_div-1);

      double l_vec[3];
      vtkMath::Subtract(lastPt, firstPt, l_vec);
      vtkMath::Normalize(l_vec);

      for (int k=0; k<l_div; k++)
      {
        pos[2] = k;

        double y_vec[3];
        for (int l=0; l<3; l++)
          y_vec[l] = l_vec[l]*k*l_dist;

        double inPt[3];
        vtkMath::Add(firstPt, y_vec, inPt);

        int newPtId = vtkStructuredData::ComputePointId(dim, pos);
        paraHexMesh->GetPoints()->SetPoint(newPtId, inPt);
      }
    }
  }

  return SV_OK;
}

// ----------------------
// CheckFace
// ----------------------
int vtkSVPolycubeGenerator::CheckFace(vtkPolyData *polycubePd, int faceId,
                                    int &nTopPts, int &nBotPts,
                                    int &flatTop, int &flatBot)
{
  double tol = 1.0e-4;
  fprintf(stdout,"CHECKING FACE: %d\n", faceId);
  vtkIdType npts, *ptIds;
  polycubePd->GetCellPoints(faceId, npts, ptIds);

  double (*pts)[3] = new double[npts][3];
  double (*vecs)[3] = new double[npts][3];

  for (int i=0; i<npts; i++)
  {
    int ptId0 = ptIds[i];
    int ptId1 = ptIds[(i+1)%npts];

    polycubePd->GetPoint(ptId0, pts[i]);
    polycubePd->GetPoint(ptId1, pts[(i+1)%npts]);

    //fprintf(stdout,"POINT %d: %.6f %.6f %.6f\n", i, pts[i][0], pts[i][1], pts[i][2]);
    //fprintf(stdout,"POINT %d: %.6f %.6f %.6f\n", (i+1)%npts, pts[(i+1)%npts][0], pts[(i+1)%npts][1], pts[(i+1)%npts][2]);
    //fprintf(stdout,"MAKES VEC: %d\n", i);

    vtkMath::Subtract(pts[(i+1)%npts], pts[i], vecs[i]);
    vtkMath::Normalize(vecs[i]);
  }

  flatTop = 0;
  nTopPts = 2;

  double testDot0 = vtkMath::Dot(vecs[0], vecs[1]);
  double testDot1 = vtkMath::Dot(vecs[0], vecs[2]);
  fprintf(stdout,"TEST DOT 0: %.8f\n", testDot0);
  fprintf(stdout,"TEST DOT 1: %.8f\n", testDot1);

  if (testDot0 < tol && testDot0 > -1.0*tol)
    flatTop = 1;
  else
  {
    if (!(fabs(testDot1) <= 1.0+tol && fabs(testDot1) >= 1.0-tol))
      nTopPts = 3;
  }

  if (testDot1 <= tol && testDot1 >= -1.0*tol)
  {
    flatTop = 1;
    nTopPts = 3;
  }

  flatBot = 0;
  nBotPts = 2;

  double testDot2 = vtkMath::Dot(vecs[npts-1], vecs[npts-2]);
  double testDot3 = vtkMath::Dot(vecs[npts-2], vecs[0]);
  fprintf(stdout,"TEST DOT 2: %.8f\n", testDot2);
  fprintf(stdout,"TEST DOT 3: %.8f\n", testDot3);

  if (testDot2 <= tol && testDot2 >= -1.0*tol)
    flatBot = 1;
  else
  {
    if (!(fabs(testDot3) <= 1.0+tol && fabs(testDot3) >= 1.0-tol))
      nBotPts=3;
  }

  if (fabs(testDot2) <= 1.0+tol && fabs(testDot2) >= 1.0-tol)
  {
    flatBot = 1;
    nBotPts = 3;
  }


  delete [] pts;
  delete [] vecs;

  return SV_OK;
}

