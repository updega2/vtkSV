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
 *  \brief This implements the vtkSVPolyDataSliceAndDiceFilter
 *
 *  \author Adam Updegrove
 *  \author updega2@gmail.com
 *  \author UC Berkeley
 *  \author shaddenlab.berkeley.edu
 */

#include "vtkSVPolyDataSliceAndDiceFilter.h"

#include "vtkAppendFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkConnectivityFilter.h"
#include "vtkCutter.h"
#include "vtkDoubleArray.h"
#include "vtkFeatureEdges.h"
#include "vtkIdFilter.h"
#include "vtkIntArray.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlanes.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVFindGeodesicPath.h"
#include "vtkSVFindSeparateRegions.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkXMLUnstructuredGridWriter.h"

#include <iostream>
#include <sstream>
#include <cmath>

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVPolyDataSliceAndDiceFilter);

// ----------------------
// Constructor
// ----------------------
vtkSVPolyDataSliceAndDiceFilter::vtkSVPolyDataSliceAndDiceFilter()
{
  this->InitialPd       = vtkPolyData::New();
  this->WorkPd          = vtkPolyData::New();
  this->GraphPd         = vtkPolyData::New();
  this->SurgeryLinesPd  = vtkPolyData::New();
  this->Polycube        = vtkSVGeneralizedPolycube::New();
  this->CenterlinesPd   = NULL;
  this->CenterlineGraph = NULL;

  this->BoundaryPointsArrayName = NULL;
  this->SegmentIdsArrayName = NULL;
  this->GroupIdsArrayName = NULL;
  this->SliceIdsArrayName = NULL;
  this->SphereRadiusArrayName = NULL;
  this->InternalIdsArrayName = NULL;
  this->DijkstraArrayName = NULL;

  this->ConstructPolycube = 0;
  this->SliceDirection = 0;
  this->TotalSliceId = 0;
  this->SliceLength = 1.5;
}

// ----------------------
// Destructor
// ----------------------
vtkSVPolyDataSliceAndDiceFilter::~vtkSVPolyDataSliceAndDiceFilter()
{
  if (this->InitialPd != NULL)
  {
    this->InitialPd->Delete();
    this->InitialPd = NULL;
  }
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->GraphPd != NULL)
  {
    this->GraphPd->Delete();
    this->GraphPd = NULL;
  }
  if (this->CenterlinesPd != NULL)
  {
    this->CenterlinesPd->UnRegister(this);
    this->CenterlinesPd = NULL;
  }
  if (this->SurgeryLinesPd != NULL)
  {
    this->SurgeryLinesPd->Delete();
    this->SurgeryLinesPd = NULL;
  }
  if (this->Polycube != NULL)
  {
    this->Polycube->Delete();
    this->Polycube = NULL;
  }
  if (this->CenterlineGraph != NULL)
  {
    delete this->CenterlineGraph;
    this->CenterlineGraph = NULL;
  }

  if (this->BoundaryPointsArrayName != NULL)
  {
    delete [] this->BoundaryPointsArrayName;
    this->BoundaryPointsArrayName = NULL;
  }
  if (this->GroupIdsArrayName != NULL)
  {
    delete [] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
  }
  if (this->SegmentIdsArrayName != NULL)
  {
    delete [] this->SegmentIdsArrayName;
    this->SegmentIdsArrayName = NULL;
  }
  if (this->SliceIdsArrayName != NULL)
  {
    delete [] this->SliceIdsArrayName;
    this->SliceIdsArrayName = NULL;
  }
  if (this->SphereRadiusArrayName != NULL)
  {
    delete [] this->SphereRadiusArrayName;
    this->SphereRadiusArrayName = NULL;
  }
  if (this->InternalIdsArrayName != NULL)
  {
    delete [] this->InternalIdsArrayName;
    this->InternalIdsArrayName = NULL;
  }
  if (this->DijkstraArrayName != NULL)
  {
    delete [] this->DijkstraArrayName;
    this->DijkstraArrayName = NULL;
  }
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVPolyDataSliceAndDiceFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Construct Polycube: " <<
    this->ConstructPolycube << "\n";
  if (this->BoundaryPointsArrayName != NULL)
    os << indent << "Boundary points array name: " << this->BoundaryPointsArrayName << "\n";
  if (this->GroupIdsArrayName != NULL)
    os << indent << "Group Ids array name: " << this->GroupIdsArrayName << "\n";
  if (this->SegmentIdsArrayName != NULL)
    os << indent << "Segment Ids array name: " << this->SegmentIdsArrayName << "\n";
  if (this->SliceIdsArrayName != NULL)
    os << indent << "Slice Ids array name: " << this->SliceIdsArrayName << "\n";
  if (this->SphereRadiusArrayName != NULL)
    os << indent << "Sphere radius array name: " << this->SphereRadiusArrayName << "\n";
  if (this->InternalIdsArrayName != NULL)
    os << indent << "Internal Ids array name: " << this->InternalIdsArrayName << "\n";
  if (this->DijkstraArrayName != NULL)
    os << indent << "Dijkstra distance array name: " << this->DijkstraArrayName << "\n";
}

// ----------------------
// RequestData
// ----------------------
int vtkSVPolyDataSliceAndDiceFilter::RequestData(vtkInformation *vtkNotUsed(request),
                                                 vtkInformationVector **inputVector,
                                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  // Copy the input to operate on
  this->InitialPd->DeepCopy(input1);

  this->WorkPd->DeepCopy(this->InitialPd);
  if (this->CenterlinesPd == NULL)
  {
    vtkErrorMacro("Error: no centerliens provided\n");
  }

  if (this->PrepFilter() != SV_OK)
  {
    vtkErrorMacro("Error in preprocessing the polydata\n");
    return SV_ERROR;
  }
  fprintf(stdout,"Graph built...\n");

  if (this->RunFilter() != SV_OK)
  {
    vtkErrorMacro("Error when running main operation\n");
    return SV_ERROR;
  }

  output->DeepCopy(this->WorkPd);
  return SV_OK;
}

// ----------------------
// PrepFilter
// ----------------------
int vtkSVPolyDataSliceAndDiceFilter::PrepFilter()
{
  // Find and label nodes that form the boundary between groups
  if (this->FindGroupBoundaries() != SV_OK)
  {
    vtkErrorMacro("Unable to find boundaries of input group ids data array");
    return SV_ERROR;
  }

  // Find the points that touch three or more different groups
  if (this->GetCriticalPoints() != SV_OK)
  {
    vtkErrorMacro("Unable to retrieve critical points");
    return SV_ERROR;
  }

  vtkNew(vtkPoints, skeletonPoints);
  vtkNew(vtkCellArray, skeletonCells);
  vtkNew(vtkIntArray, skeletonGroupIds);
  this->CenterlineGraph = new svGraph(0, this->CenterlinesPd,
                                     this->GroupIdsArrayName,
                                     this->CriticalPointMap);
  if (this->CenterlineGraph->BuildGraph() != SV_OK)
  {
    vtkErrorMacro("Unable to form skeleton of polydata");
    return SV_ERROR;
  }
  this->CenterlineGraph->GetGraphPolyData(this->GraphPd);

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVPolyDataSliceAndDiceFilter::RunFilter()
{
  if (this->SliceBifurcations() != SV_OK)
  {
    vtkErrorMacro("Error in slicing the polydata\n");
    return SV_ERROR;
  }
  fprintf(stdout,"Bifurcations segmented...\n");

  if (this->SliceBranches() != SV_OK)
  {
    vtkErrorMacro("Error in slicing the polydata\n");
    return SV_ERROR;
  }
  fprintf(stdout,"Branches segmented...\n");

  if (this->ConstructPolycube)
  {
    if (this->BuildPolycube() != SV_OK)
    {
      vtkErrorMacro("Error in constructing polycube\n");
      return SV_ERROR;
    }
    fprintf(stdout,"Polycube built...\n");
  }

  return SV_OK;
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
const int vtkSVPolyDataSliceAndDiceFilter::DT[6][4] =
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
const int vtkSVPolyDataSliceAndDiceFilter::RT[9][8] =
{{4, 0, 3, 7, 5, 1, 2, 6},  // x
 {4, 5, 1, 0, 7, 6, 2, 3},  // y
 {1, 2, 3, 0, 5, 6, 7, 4},  // z
 {5, 4, 7, 6, 1, 0, 3, 2},  // xx
 {7, 6, 5, 4, 3, 2, 1, 0},  // yy
 {2, 3, 0, 1, 6, 7, 4, 5},  // zz
 {1, 5, 6, 2, 0, 4, 7, 3},  // xxx
 {3, 2, 6, 7, 0, 1, 5, 4},  // yyy
 {3, 0, 1, 2, 7, 4, 5, 6}}; // zzz

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
int vtkSVPolyDataSliceAndDiceFilter::LookupIndex(const int parent, const int divchild, const int index)
{
  if (parent == RIGHT)
  {
    if (divchild == BACK)
      return vtkSVPolyDataSliceAndDiceFilter::RT[0][vtkSVPolyDataSliceAndDiceFilter::RT[5][vtkSVPolyDataSliceAndDiceFilter::RT[1][index]]];
    if (divchild == FRONT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[0][vtkSVPolyDataSliceAndDiceFilter::RT[7][index]];
    if (divchild == UP)
      return vtkSVPolyDataSliceAndDiceFilter::RT[7][index];
    if (divchild == DOWN)
      return vtkSVPolyDataSliceAndDiceFilter::RT[5][vtkSVPolyDataSliceAndDiceFilter::RT[1][index]];
  }
  if (parent == LEFT)
  {
    if (divchild == BACK)
      return vtkSVPolyDataSliceAndDiceFilter::RT[0][vtkSVPolyDataSliceAndDiceFilter::RT[1][index]];
    if (divchild == FRONT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[6][vtkSVPolyDataSliceAndDiceFilter::RT[1][index]];
    if (divchild == UP)
      return vtkSVPolyDataSliceAndDiceFilter::RT[3][vtkSVPolyDataSliceAndDiceFilter::RT[1][index]];
    if (divchild == DOWN)
      return vtkSVPolyDataSliceAndDiceFilter::RT[1][index];
  }
  if (parent == FRONT)
  {
    if (divchild == RIGHT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[6][index];
    if (divchild == LEFT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[0][vtkSVPolyDataSliceAndDiceFilter::RT[4][index]];
    if (divchild == UP)
      return vtkSVPolyDataSliceAndDiceFilter::RT[2][vtkSVPolyDataSliceAndDiceFilter::RT[3][vtkSVPolyDataSliceAndDiceFilter::RT[1][index]]];
    if (divchild == DOWN)
      return vtkSVPolyDataSliceAndDiceFilter::RT[1][vtkSVPolyDataSliceAndDiceFilter::RT[6][index]];
  }
  if (parent == BACK)
  {
    if (divchild == RIGHT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[0][index];
    if (divchild == LEFT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[4][vtkSVPolyDataSliceAndDiceFilter::RT[0][index]];
    if (divchild == UP)
      return vtkSVPolyDataSliceAndDiceFilter::RT[7][vtkSVPolyDataSliceAndDiceFilter::RT[0][index]];
    if (divchild == DOWN)
      return vtkSVPolyDataSliceAndDiceFilter::RT[1][vtkSVPolyDataSliceAndDiceFilter::RT[1][index]];
  }
  if (parent == UP)
  {
    if (divchild == RIGHT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[3][index];
    if (divchild == LEFT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[4][index];
    if (divchild == BACK)
      return vtkSVPolyDataSliceAndDiceFilter::RT[2][vtkSVPolyDataSliceAndDiceFilter::RT[4][index]];
    if (divchild == FRONT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[2][vtkSVPolyDataSliceAndDiceFilter::RT[3][index]];
  }
  if (parent == DOWN)
  {
    if (divchild == RIGHT)
      return index;
    if (divchild == LEFT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[5][index];
    if (divchild == BACK)
      return vtkSVPolyDataSliceAndDiceFilter::RT[2][index];
    if (divchild == FRONT)
      return vtkSVPolyDataSliceAndDiceFilter::RT[8][index];
  }

  return -1;
}

// ----------------------
// FindGroupBoundaries
// ----------------------
/**
 * \details Calls vtkSVFindSeparateRegions to find the points that are boundaries
 * between all different Group Id numbers
 */
int vtkSVPolyDataSliceAndDiceFilter::FindGroupBoundaries()
{
  vtkNew(vtkSVFindSeparateRegions, separator);
  separator->SetInputData(this->WorkPd);
  separator->SetOutPointArrayName(this->BoundaryPointsArrayName);
  separator->SetArrayName(this->GroupIdsArrayName);
  separator->Update();

  this->WorkPd->DeepCopy(separator->GetOutput());

  return SV_OK;
}

// ----------------------
// GetCriticalPoints
// ----------------------
/**
 * \details Forms a map of groups to critical points.
 * Keys: Group ids on surface.
 * Values: Ids of points on surface that touch key group.
 */
int vtkSVPolyDataSliceAndDiceFilter::GetCriticalPoints()
{
  // Get points and boundary point data off of surface
  this->WorkPd->BuildLinks();
  int numPoints = this->WorkPd->GetNumberOfPoints();
  vtkDataArray *boundaryPointArray =
    this->WorkPd->GetPointData()->GetArray(this->BoundaryPointsArrayName);

  // Loop through all points finding
  for (int i=0; i<numPoints; i++)
  {
    int isBoundary = boundaryPointArray->GetTuple1(i);
    vtkNew(vtkIdList, groupIds);
    // If this point is on the boundary, find out if it separate mores than 2 groups!
    if (isBoundary)
    {
      vtkSVGeneralUtils::GetPointGroups(this->WorkPd, this->GroupIdsArrayName, i, groupIds);
      int pointType = groupIds->GetNumberOfIds();
      // This should not be possible!
      if (pointType < 2)
        vtkErrorMacro("Point incorrectly identified");
      // Oh boy, we found ourselves a critical point!
      if (pointType == 3)
        this->InsertCriticalPoints(i, groupIds);
    }
  }

  return SV_OK;
}

// ----------------------
// InsertCriticalPoints
// ----------------------
int vtkSVPolyDataSliceAndDiceFilter::InsertCriticalPoints(const int pointId, vtkIdList *groupIds)
{
  int numIds = groupIds->GetNumberOfIds();

  for (int i=0; i<numIds; i++)
  {
    int groupId = groupIds->GetId(i);
    this->CriticalPointMap.insert(std::make_pair(groupId, pointId));
    vtkDebugMacro("Added critical point pair: "<< groupId << " " << pointId);
  }

  return SV_OK;
}

// ----------------------
// GetBranch
// ----------------------
int vtkSVPolyDataSliceAndDiceFilter::GetBranch(const int branchId, vtkPolyData *branchPd,
                                             vtkPolyData *branchCenterlinesPd)
{
  vtkSVGeneralUtils::ThresholdPd(this->WorkPd, branchId, branchId, 1,
    this->SegmentIdsArrayName, branchPd);

  // Make sure that polydata of branchId exists
  if (branchPd != NULL)
  {
    vtkNew(vtkPolyData, centerlineBranchPd);
    vtkSVGeneralUtils::ThresholdPd(this->CenterlinesPd, branchId, branchId, 1,
      this->GroupIdsArrayName, centerlineBranchPd);

    //Need to get just first cell of centerlines.
    //There are sometimes duplicate for each centerline running through
    vtkNew(vtkIdFilter, ider);
    ider->SetInputData(centerlineBranchPd);
    ider->SetIdsArrayName(this->InternalIdsArrayName);
    ider->Update();
    vtkNew(vtkPolyData, tmpPd);
    tmpPd->ShallowCopy(ider->GetOutput());

    // Threshold out id 0
    vtkSVGeneralUtils::ThresholdPd(tmpPd, 0, 0, 1,
      this->InternalIdsArrayName, branchCenterlinesPd);

    // Remove the tmp Ids from the centerlines
    branchCenterlinesPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
    branchCenterlinesPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);
  }

  return SV_OK;
}

// ----------------------
// SliceBranches
// ----------------------
int vtkSVPolyDataSliceAndDiceFilter::SliceBranches()
{
  // Set up new array for slicing. Typically one branch is one slice, but
  // just in case, we set up for multiple slices per branch
  int numIds = this->WorkPd->GetCellData()->
    GetArray(this->SegmentIdsArrayName)->GetNumberOfTuples();
  vtkNew(vtkIntArray, sliceIds);
  sliceIds->SetNumberOfTuples(numIds);
  sliceIds->SetName(this->SliceIdsArrayName);
  sliceIds->FillComponent(0, -1);
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(this->WorkPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  this->WorkPd->DeepCopy(ider->GetOutput());

  // Set up points and data on the surgery lines
  vtkNew(vtkPoints, surgeryPts);
  vtkNew(vtkIntArray, surgeryData);
  surgeryData->SetName(this->InternalIdsArrayName);
  vtkNew(vtkIntArray, surgeryCellData);
  surgeryCellData->SetName(this->SegmentIdsArrayName);
  this->SurgeryLinesPd->SetPoints(surgeryPts);
  this->SurgeryLinesPd->GetPointData()->AddArray(surgeryData);
  this->SurgeryLinesPd->Allocate(this->CenterlineGraph->NumberOfCells);
  this->SurgeryLinesPd->GetCellData()->AddArray(surgeryCellData);

  // Run through each branch and call the slice for each one!
  int numSegs = this->CenterlineGraph->NumberOfCells;
  for (int i=0; i<numSegs; i++)
  {
    // Branch node in graph
    svGCell *gCell = this->CenterlineGraph->GetCell(i);
    int branchId = gCell->GroupId;

    // Get branch polydata and centerlines
    vtkNew(vtkPolyData, branchPd);
    vtkNew(vtkPolyData, branchCenterline);
    this->GetBranch(branchId, branchPd, branchCenterline);

    // Check that we actually have something and slice it up!
    int numPoints = branchPd->GetNumberOfPoints();
    if (numPoints != 0)
    {
      this->SliceBranch(branchPd, branchCenterline, gCell, sliceIds, 0);
      this-surgeryCellData->InsertNextValue(branchId);
    }
  }

  // Add data that was retrieved in branch slicing
  this->WorkPd->GetCellData()->AddArray(sliceIds);
  this->WorkPd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  this->WorkPd->GetPointData()->RemoveArray(this->InternalIdsArrayName);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetSectionZAxis(const double endPt[3], const double startPt[3],
                                                   double zvec[3])
{
  //Get approximate z axis of section
  vtkMath::Subtract(startPt, endPt, zvec);
  vtkMath::Normalize(zvec);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetSectionXAxis(const double endPt[3], const double startPt[3],
                                                   const double surfacePt[3],
                                                   double xvec[3])
{
  double vec0[3], vec1[3];

  //Get approximate x axis of section
  vtkMath::Subtract(surfacePt, startPt, vec0);
  vtkMath::Subtract(endPt, startPt, vec1);
  vtkMath::Normalize(vec1);
  double vecmag = vtkMath::Dot(vec0, vec1);
  vtkMath::MultiplyScalar(vec1, vecmag);
  vtkMath::Subtract(vec0, vec1, xvec);

  return SV_OK;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::DetermineSliceStrategy(vtkPolyData *branchPd,
                                                          svGCell *gCell,
                                                          vtkPolyData *branchCenterline,
                                                          int &branchStartPtId,
                                                          vtkIdList *surgeryPoints,
                                                          int &centerlineStartPtId,
                                                          int &strategy)
{
  int numCenterlinePts = branchCenterline->GetNumberOfPoints();
  int branchId = gCell->GroupId;
  double topPt[3], bottomPt[3];
  branchCenterline->GetPoint(0, topPt);
  branchCenterline->GetPoint(numCenterlinePts - 1, bottomPt);

  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(branchPd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();

  vtkNew(vtkPolyData, topPoly);
  vtkNew(vtkPolyData, bottomPoly);
  vtkSVGeneralUtils::GetClosestPointConnectedRegion(boundaries->GetOutput(), topPt, topPoly);
  vtkSVGeneralUtils::GetClosestPointConnectedRegion(boundaries->GetOutput(), bottomPt, bottomPoly);

  vtkDataArray *topPointIds =
      topPoly->GetPointData()->GetArray(this->InternalIdsArrayName);
  vtkDataArray *bottomPointIds =
      bottomPoly->GetPointData()->GetArray(this->InternalIdsArrayName);

  int defaultDirs[6];
  defaultDirs[DOWN]  = RIGHT;
  defaultDirs[UP]    = LEFT;
  defaultDirs[RIGHT] = UP;
  defaultDirs[LEFT]  = DOWN;
  defaultDirs[BACK]  = RIGHT;
  defaultDirs[FRONT] = RIGHT;
  double cellIndices[8];
  int numSurgeryPoints = 0;
  for (int i=0; i<8; i++)
  {
    if (gCell->CornerPtIds[i] != -1)
      numSurgeryPoints++;
    if (branchId == 0 || gCell->Children[0] != NULL)
      cellIndices[i] = this->LookupIndex(gCell->Dir, gCell->Children[gCell->DivergingChild]->Dir, i);
    else
      cellIndices[i] = this->LookupIndex(gCell->Dir, defaultDirs[gCell->Dir], i);

  }
  this->SliceDirection = 0;
  centerlineStartPtId = 0;
  fprintf(stdout,"NUM SURGERY POINTS!: %d\n", numSurgeryPoints);
  if (numSurgeryPoints == 0)
  {
    strategy = 0;
    fprintf(stdout,"Strategy 0, any surgery points will do\n");
    //There are no critical points. Any point on the boundary can be used as
    //the start point. This also means this is a genus zero surface with no
    //branching locations. Needs to be revised for non vascular models where
    //there is not inlet an outlet.
    branchStartPtId = bottomPointIds->GetTuple1(0);
    fprintf(stdout,"Start point is: %d\n", branchStartPtId);
    double secondPt[3], sidePt[3], zvec[3], xvec[3];
    branchCenterline->GetPoint(numCenterlinePts - 2, secondPt);
    bottomPoly->GetPoint(0, sidePt);
    vtkMath::Subtract(bottomPt, secondPt, zvec);
    vtkMath::Normalize(zvec);
    vtkMath::Subtract(sidePt, bottomPt, xvec);
    vtkMath::Normalize(xvec);
    this->GetFirstSurgeryPoints(bottomPoly, 0, surgeryPoints, xvec, zvec);
    for (int i=0; i<4; i++)
      gCell->CornerPtIds[i] = surgeryPoints->GetId(i);
    this->SliceDirection = 1;
    centerlineStartPtId = numCenterlinePts - 1;
  }
  else if (numSurgeryPoints == 4 || numSurgeryPoints == 8)
  {
    if (numSurgeryPoints == 4)
    {
      strategy = 1;
      fprintf(stdout,"Strategy 1, must actively use surgery points from bifurcations segmentation\n");
    }
    else
    {
      strategy = 2;
      fprintf(stdout,"Strategy 2, interior segment, must actively use all end surgery points\n");
    }
    //There are two critical points. This means that these points should be used
    //as corners of a polycube along with two other points.
    int front;
    int back;
    if (branchId == 0)
    {
      fprintf(stdout,"Centerline upside down, starting from bottom\n");
      this->SliceDirection = 1;
      centerlineStartPtId = numCenterlinePts - 1;
      for (int i=4; i<8; i++)
        surgeryPoints->InsertNextId(gCell->CornerPtIds[int(cellIndices[i])]);
      branchStartPtId = surgeryPoints->GetId(4);
    }
    else
    {
      for (int i=0; i<4; i++)
        surgeryPoints->InsertNextId(gCell->CornerPtIds[int(cellIndices[i])]);
      branchStartPtId = surgeryPoints->GetId(0);
    }
  }
  else
  {
    vtkErrorMacro("This shouldnt be possible. Something went wrong");
    return SV_ERROR;
  }

  return SV_OK;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetSurgeryPoints(vtkPolyData *pd,
                                                    vtkPolyData *parentPd,
                                                    vtkDataArray *pointIds,
                                                    const double clStartPt[3],
                                                    const double clSecondPt[3],
                                                    const int front,
                                                    const int back,
                                                    const int checkId,
                                                    std::string arrayName,
                                                    vtkIdList *surgeryPoints,
                                                    double startDir[3])
{
  int pointId = front;
  vtkNew(vtkIdList, startCellIds);
  pd->GetPointCells(pointId, startCellIds);
  int prevCellId = startCellIds->GetId(0);
  int prevCellId2 = startCellIds->GetId(1);
  vtkIdType npts, *pts;
  pd->GetCellPoints(prevCellId, npts, pts);
  int secondPtId;
  if (pts[0] == pointId)
    secondPtId = pts[1];
  else
    secondPtId = pts[0];
  vtkNew(vtkIdList, checkIds);
  vtkSVGeneralUtils::GetPointGroups(parentPd, arrayName, pd->GetPointData()->GetArray(
    this->InternalIdsArrayName)->GetTuple1(secondPtId), checkIds);
  fprintf(stdout,"Checking the id of %f\n", pd->GetPointData()->GetArray(
      this->InternalIdsArrayName)->GetTuple1(secondPtId));
  if (checkIds->IsId(checkId) != -1)
  {
    // It is actually opposite of what you would think. You want to check if
    // id of branch you want to segment is on this piece. If it is, you must switch,
    // because then that cell is used as the "previous cell", as if we were already
    // looping around
    fprintf(stdout,"This says that we would flip\n");
    int tmp     = prevCellId;
    prevCellId  = prevCellId2;
    prevCellId2 = tmp;
  }
  double pt0[3], pt1[3];
  pd->GetCellPoints(prevCellId, npts, pts);
  if (pts[0] == pointId)
    secondPtId = pts[1];
  else
    secondPtId = pts[0];
  pd->GetPoint(pointId, pt0);
  pd->GetPoint(secondPtId, pt1);
  double vec0[3], vec1[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  //Get the vector from start to end
  double zvec[3];
  this->GetSectionZAxis(clSecondPt, clStartPt, zvec);

  //Get starting point from tmp id
  double xvec[3];
  this->GetSectionXAxis(clSecondPt, clStartPt, pt0, xvec);

  vtkMath::Cross(xvec, vec0, vec1);

  //if (vtkMath::Dot(zvec, vec1) > 0)
  //{
  //  fprintf(stdout,"Inside other way!\n");
  //  prevCellId = startCellIds->GetId(1);
  //  prevCellId2 = startCellIds->GetId(0);
  //}
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);
  vtkMath::Cross(vec0, zvec, startDir);
  vtkMath::Normalize(startDir);
  if (vtkMath::Dot(startDir, xvec) > 0)
  {
    for (int i=0; i<3; i++)
      startDir[i] = -1.0*startDir[i];
  }


  vtkNew(vtkIdList, halfPoints);
  this->GetHalfSurgeryPoints(pd, pointIds, prevCellId, front, back, halfPoints);
  this->GetHalfSurgeryPoints(pd, pointIds, prevCellId2, front, back, halfPoints);
  surgeryPoints->InsertNextId(halfPoints->GetId(0));
  surgeryPoints->InsertNextId(halfPoints->GetId(1));
  surgeryPoints->InsertNextId(halfPoints->GetId(3));
  surgeryPoints->InsertNextId(halfPoints->GetId(2));

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetHalfSurgeryPoints(vtkPolyData *pd,
                                                        vtkDataArray *pointIds,
                                                        const int cellId,
                                                        const int front,
                                                        const int back,
                                                        vtkIdList *surgeryPoints)
{
  //Make front and back always start with the lowest ids so that same ids are
  //retrieved for every boundary
  pd->BuildLinks();
  int pointId = front;
  int prevCellId = cellId;

  //Getting length of the loop separated into the two sides of the boundary
  //by the critical points
  double length = 0.0;
  while (pointId != back)
  {
    double pt0[3], pt1[3];
    pd->GetPoint(pointId, pt0);
    vtkSVGeneralUtils::IteratePoint(pd, pointId, prevCellId);
    pd->GetPoint(pointId, pt1);
    length += vtkSVGeneralUtils::Distance(pt0, pt1);
  }

  //Finding the surgery points which separate each half of the boundary into
  //three slices
  double surgeryCount = 1.0;
  double currLength = 0.0;;
  pointId = front;
  prevCellId = cellId;
  while (pointId != back)
  {
    double pt0[3], pt1[3];
    pd->GetPoint(pointId, pt0);
    vtkSVGeneralUtils::IteratePoint(pd, pointId, prevCellId);
    pd->GetPoint(pointId, pt1);
    currLength += vtkSVGeneralUtils::Distance(pt0, pt1);
    if (currLength > length*(surgeryCount/3.0))
    {
      surgeryPoints->InsertNextId(pointIds->GetTuple1(pointId));
      surgeryCount += 1.0;
    }
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVPolyDataSliceAndDiceFilter::CheckLength(int &ptId, const int numPts,
                                                int &done)
{
  if (this->SliceDirection == 0)
  {
    if (ptId > numPts - 2)
    {
      ptId = numPts - 2;
      done = 1;
    }
  }
  else if (this->SliceDirection == 1)
  {
    if (ptId < 1)
    {
      ptId = 1;
      done = 1;
    }
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVPolyDataSliceAndDiceFilter::UpdatePtId(int &ptId)
{
  if (this->SliceDirection == 0)
  {
    ptId++;
  }
  else if (this->SliceDirection == 1)
  {
    ptId--;
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::SliceBranch(vtkPolyData *branchPd,
                                               vtkPolyData *branchCenterline,
                                               svGCell *gCell,
                                               vtkDataArray *sliceIds,
                                               int secondRun)
{
  int branchId = gCell->GroupId;
  fprintf(stdout,"Slicing branch: %d\n", branchId);

  int numCells = branchPd->GetNumberOfCells();
  int numCenterlinePts = branchCenterline->GetNumberOfPoints();
  vtkDataArray *pointIds = branchPd->GetPointData()->GetArray(this->InternalIdsArrayName);
  vtkDataArray *radiusArray = branchCenterline->GetPointData()->GetArray(this->SphereRadiusArrayName);

  vtkNew(vtkPolyData, leftovers); leftovers->DeepCopy(branchPd);

  double totalLength = 0.0;
  vtkSVGeneralUtils::GetPointsLength(branchCenterline, totalLength);

  int linePtId = 0;
  int startPtId = -1;
  double xvec[3], zvec[3];
  vtkNew(vtkIdList, surgeryPoints);
  int strategy = 0;
  this->DetermineSliceStrategy(branchPd, gCell, branchCenterline,
                               startPtId, surgeryPoints, linePtId, strategy);
  int defaultDirs[6];
  defaultDirs[DOWN]  = RIGHT;
  defaultDirs[UP]    = LEFT;
  defaultDirs[RIGHT] = UP;
  defaultDirs[LEFT]  = DOWN;
  defaultDirs[BACK]  = RIGHT;
  defaultDirs[FRONT] = RIGHT;
  int newIndices[8];
  for (int i=0; i<8; i++)
  {
    if (branchId == 0 || gCell->Children[0] != NULL)
      newIndices[i] = this->LookupIndex(gCell->Dir, gCell->Children[gCell->DivergingChild]->Dir, i);
    else
      newIndices[i] = this->LookupIndex(gCell->Dir, defaultDirs[gCell->Dir], i);
  }

  //if (this->SliceDirection == 1)
  //{
  //  for (int i=0; i<4; i++)
  //    gCell->CornerPtIds[newIndices[i+4]] = surgeryPoints->GetId(i);
  //}
  //else
  //{
    for (int i=0; i<4; i++)
    {
      //gCell->CornerPtIds[newIndices[i]] = surgeryPoints->GetId(i);
      fprintf(stdout,"Starting are: %lld\n", surgeryPoints->GetId(i));
    }
  //}

  if (branchId != 0)
  {
    double veryFirst[3], verySecond[3], surfacePt[3];
    branchCenterline->GetPoint(0, veryFirst);
    branchCenterline->GetPoint(1, verySecond);
    this->GetSectionZAxis(verySecond, veryFirst, gCell->ZVec);
    int pointId = pointIds->LookupValue(surgeryPoints->GetId(0));
    branchPd->GetPoint(pointId, surfacePt);
    this->GetSectionXAxis(verySecond, veryFirst, surfacePt, gCell->XVec);
  }

  vtkNew(vtkIdList, surgeryLineIds);
  int done = 0;
  int sliceId = 0;
  double currLength = 0.0;
  while (!done)
  {
    double inscribedRadius = radiusArray->GetTuple1(linePtId);
    double sliceLength = inscribedRadius;
    if (inscribedRadius > 0.5)
    {
      sliceLength *= this->SliceLength;
    }
    else
    {
      sliceLength *= this->SliceLength * (4.0/3.0);
    }
    if ((currLength + 2.0*this->SliceLength*inscribedRadius) >= totalLength)
    {
      sliceLength = this->SliceLength*(totalLength - currLength);
      done = 1;
    }

    double centerlineLength = 0;
    double startPt[3], pt0[3], pt1[3];
    branchCenterline->GetPoint(linePtId, startPt);

    while (centerlineLength < sliceLength)
    {
      this->CheckLength(linePtId, numCenterlinePts, done);
      branchCenterline->GetPoint(linePtId, pt0);
      this->UpdatePtId(linePtId);
      branchCenterline->GetPoint(linePtId, pt1);
      centerlineLength += vtkSVGeneralUtils::Distance(pt0, pt1);
    }
    currLength += centerlineLength;

    //Get the vector from start to end
    this->GetSectionZAxis(pt1, startPt, zvec);

    //Get starting point from tmp id
    double surfacePt[3];
    int pointId = pointIds->LookupValue(surgeryPoints->GetId(0));
    branchPd->GetPoint(pointId, surfacePt);
    this->GetSectionXAxis(pt1, startPt, surfacePt, xvec);

    //Get the cut plane
    double origin[3];
    vtkNew(vtkPlane, cutPlane);
    vtkSVGeneralUtils::GetCutPlane(pt1, pt0, centerlineLength, origin, cutPlane);

    //Cut the pds
    vtkNew(vtkPolyData, slicePd);
    vtkSVGeneralUtils::ExtractionCut(leftovers, cutPlane, 0, 1, slicePd);

    //Get closestPoint region
    vtkNew(vtkPolyData, connectedPd);
    vtkSVGeneralUtils::GetClosestPointConnectedRegion(slicePd, origin, connectedPd);

    //fprintf(stdout,"What is val of done: %d\n", done);
    if (connectedPd->GetNumberOfCells() == numCells)
    {
      done = 1;
    }
    else
    {
      if (this->CheckSlice(connectedPd) != SV_OK)
      {
        continue;
      }
    }

    double contourClosePt[3];
    branchCenterline->GetPoint(linePtId, contourClosePt);
    double contourRadius = radiusArray->GetTuple1(linePtId);
    if (!done)
    {
      this->GetNextSurgeryPoints(connectedPd, contourClosePt, surgeryPoints, xvec, zvec, contourRadius, surgeryLineIds);
    }
    sliceId++;
  }

  double contourClosePt[3];
  branchCenterline->GetPoint(linePtId, contourClosePt);
  double contourRadius = radiusArray->GetTuple1(linePtId);
  this->GetEndSurgeryPoints(branchPd, gCell, contourClosePt, surgeryPoints, gCell->CornerPtIds, xvec, zvec, contourRadius, surgeryLineIds, newIndices, secondRun);
  if(secondRun)
  {
    this->SliceBranch(branchPd, branchCenterline, gCell, sliceIds, secondRun);
    return SV_OK;
  }

  if (this->SliceDirection == 1)
  {
    gCell->CornerPtIds[newIndices[0]] = surgeryPoints->GetId(0);
    gCell->CornerPtIds[newIndices[1]] = surgeryPoints->GetId(3);
    gCell->CornerPtIds[newIndices[2]] = surgeryPoints->GetId(2);
    gCell->CornerPtIds[newIndices[3]] = surgeryPoints->GetId(1);
    double veryFirst[3], verySecond[3], surfacePt[3];
    branchCenterline->GetPoint(0, veryFirst);
    branchCenterline->GetPoint(1, verySecond);
    this->GetSectionZAxis(verySecond, veryFirst, gCell->ZVec);
    int pointId = pointIds->LookupValue(surgeryPoints->GetId(0));
    branchPd->GetPoint(pointId, surfacePt);
    this->GetSectionXAxis(verySecond, veryFirst, surfacePt, gCell->XVec);
  }
  else
  {
    for (int i=0; i<4; i++)
    {
      fprintf(stdout,"Ending are: %lld\n", surgeryPoints->GetId(i));
      gCell->CornerPtIds[newIndices[i+4]] = surgeryPoints->GetId(i);
    }
  }
  this->AddSurgeryPoints(surgeryLineIds);

  vtkSVGeneralUtils::ReplaceDataOnCells(branchPd, sliceIds,
                                        this->TotalSliceId, -1, this->InternalIdsArrayName);
  this->TotalSliceId++;

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::CheckSlice(vtkPolyData *pd)
{
  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(pd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();

  vtkNew(vtkConnectivityFilter, connector);
  connector->SetInputData(boundaries->GetOutput());
  connector->SetExtractionModeToAllRegions();
  connector->Update();

  if (connector->GetNumberOfExtractedRegions() != 2)
  {
    return SV_ERROR;
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::SliceBifurcations()
{
  int numIds = this->WorkPd->GetCellData()->
    GetArray(this->GroupIdsArrayName)->GetNumberOfTuples();
  vtkNew(vtkIntArray, segmentIds);
  segmentIds->SetNumberOfTuples(numIds);
  segmentIds->SetName(this->SegmentIdsArrayName);
  segmentIds->CopyComponent(0, this->WorkPd->GetCellData()->
    GetArray(this->GroupIdsArrayName), 0);
  this->WorkPd->GetCellData()->AddArray(segmentIds);
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(this->WorkPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  this->WorkPd->DeepCopy(ider->GetOutput());
  vtkNew(vtkIntArray, surgeryPointArray);
  surgeryPointArray->SetNumberOfComponents(this->CenterlineGraph->NumberOfNodes);
  surgeryPointArray->SetNumberOfTuples(this->WorkPd->GetNumberOfPoints());
  for (int i=0; i<this->CenterlineGraph->NumberOfNodes; i++)
  {
    surgeryPointArray->FillComponent(i, -1);
  }
  surgeryPointArray->SetName("SurgeryPoints");
  this->WorkPd->GetPointData()->AddArray(surgeryPointArray);

  vtkNew(vtkPolyData, bifurcationPd);
  bifurcationPd->DeepCopy(this->WorkPd);
  int numSegs = this->CenterlineGraph->NumberOfCells;
  for (int i=0; i<numSegs; i++)
  {
    svGCell *gCell = this->CenterlineGraph->GetCell(i);
    if (gCell->Children[0] != NULL && gCell->Children[1] != NULL)
    {
      fprintf(stdout,"Slicing bifucation %d\n", i);
      this->SliceBifurcation(bifurcationPd, gCell);
    }
  }

  this->WorkPd->DeepCopy(bifurcationPd);

  //Pass data to polycube
  for (int i=0; i<this->CenterlineGraph->NumberOfCells; i++)
  {
    svGCell *gCell = this->CenterlineGraph->GetCell(i);
    int groupId = gCell->GroupId;
    vtkNew(vtkIntArray, singleCompArray);
    singleCompArray->SetNumberOfComponents(1);
    singleCompArray->SetNumberOfTuples(this->WorkPd->GetNumberOfPoints());
    singleCompArray->CopyComponent(0, this->WorkPd->GetPointData()->GetArray("SurgeryPoints"), groupId);
    fprintf(stdout,"BREAK\n");
    for (int j=0; j<8; j++)
    {
      gCell->CornerPtIds[j] = singleCompArray->LookupValue(j);
      fprintf(stdout,"ALL POINTS for %d: %d\n", gCell->GroupId, gCell->CornerPtIds[j]);
    }
  }
  this->WorkPd->GetPointData()->RemoveArray("SurgeryPoints");

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::SliceBifurcation(vtkPolyData *pd,
                                                    svGCell *gCell)
{
  int parentId  = gCell->GroupId;
  int segmentId = parentId + 1;
  int goodKidId = gCell->Children[gCell->AligningChild]->GroupId;
  int badKidId  = gCell->Children[gCell->DivergingChild]->GroupId;

  vtkNew(vtkPolyData, centerlineBranchPd);
  vtkSVGeneralUtils::ThresholdPd(this->CenterlinesPd, badKidId, badKidId, 1,
    this->GroupIdsArrayName, centerlineBranchPd);
  vtkNew(vtkPolyData, branchCenterline);
  //Need to get just first cell of centerlines. There are duplicate for each centerline running through
  vtkNew(vtkIdFilter, ider);
  ider->SetInputData(centerlineBranchPd);
  ider->SetIdsArrayName(this->InternalIdsArrayName);
  ider->Update();
  vtkNew(vtkPolyData, tmpPd);
  tmpPd->ShallowCopy(ider->GetOutput());
  vtkSVGeneralUtils::ThresholdPd(tmpPd, 0, 0, 1,
                                             this->InternalIdsArrayName,
                                             branchCenterline);
  double inPt0[3], inPt1[3];
  branchCenterline->GetPoint(0, inPt0);
  branchCenterline->GetPoint(1, inPt1);

  vtkNew(vtkPolyData, section0Pd);
  vtkNew(vtkPolyData, section1Pd);
  vtkNew(vtkPolyData, section2Pd);
  vtkNew(vtkPolyData, theRestPd);
  vtkNew(vtkPolyData, section2Loop);
  this->GetFourPolyDataRegions(pd, parentId, section0Pd, goodKidId, section1Pd, badKidId, section2Pd, theRestPd);

  vtkNew(vtkPolyData, special2Pd);
  vtkSVGeneralUtils::ThresholdPd(this->WorkPd, badKidId, badKidId, 1, this->GroupIdsArrayName, special2Pd);
  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(special2Pd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();
  vtkSVGeneralUtils::GetClosestPointConnectedRegion(boundaries->GetOutput(), inPt0, section2Loop);

  vtkNew(vtkAppendPolyData, appender);
  appender->AddInputData(section0Pd);
  appender->AddInputData(section1Pd);
  appender->Update();

  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(appender->GetOutput());
  cleaner->Update();

  vtkNew(vtkSVFindSeparateRegions, separator);
  separator->SetInputData(cleaner->GetOutput());
  separator->SetOutPointArrayName(this->BoundaryPointsArrayName);
  separator->SetArrayName(this->GroupIdsArrayName);
  separator->Update();

  vtkIntArray *isBoundary = vtkIntArray::SafeDownCast(
    separator->GetOutput()->GetPointData()->GetArray(this->BoundaryPointsArrayName));
  vtkNew(vtkPoints, ringPoints);
  for (int i=0; i< separator->GetOutput()->GetNumberOfPoints(); i++)
  {
    if (isBoundary->GetValue(i) == 1)
    {
      ringPoints->InsertNextPoint(separator->GetOutput()->GetPoint(i));
    }
  }
  double centroid[3];
  vtkSVGeneralUtils::GetCentroidOfPoints(ringPoints, centroid);

  std::list<int> criticalPoints;
  vtkSVGeneralUtils::GetCommonValues(this->CriticalPointMap, parentId, goodKidId, criticalPoints);
  if (criticalPoints.size() != 2)
  {
    fprintf(stderr,"There should be two critical points between groups and there are %lu\n", criticalPoints.size());
    return SV_ERROR;
  }


  vtkNew(vtkIdList, slicePoints);
  vtkDataArray *section2Ids =
    section2Loop->GetPointData()->GetArray(this->InternalIdsArrayName);
  int front = section2Ids->LookupValue(criticalPoints.front());
  int back = section2Ids->LookupValue(criticalPoints.back());
  fprintf(stdout,"WHAT IS FRONT DIR: %.4f %.4f %.4f\n", gCell->FrontDir[0], gCell->FrontDir[1],
                                                        gCell->FrontDir[2]);
  this->GetCorrectFrontPoint(section2Loop, gCell->FrontDir, front, back);

  double pt0[3], pt1[3];
  this->WorkPd->GetPoint(section2Ids->GetTuple1(front), pt0);
  this->WorkPd->GetPoint(section2Ids->GetTuple1(back), pt1);
  double vec0[3], vec1[3], normal[3], normal2[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Subtract(centroid, pt0, vec1);
  vtkMath::Cross(vec0, vec1, normal);
  vtkMath::Cross(vec0, vec1, normal2);
  vtkMath::Normalize(normal);
  vtkMath::Normalize(normal2);

  fprintf(stdout,"FRONT!!!!: %f\n", section2Ids->GetTuple1(front));
  fprintf(stdout,"BACK!!!!: %f\n", section2Ids->GetTuple1(back));
  double dummy[3];
  this->GetSurgeryPoints(section2Loop, this->WorkPd, section2Ids, inPt0, inPt1, front, back, parentId, this->GroupIdsArrayName, slicePoints, dummy);
  fprintf(stdout,"ORDER: %lld %lld %lld %lld\n", slicePoints->GetId(0), slicePoints->GetId(1), slicePoints->GetId(2), slicePoints->GetId(3));
  vtkNew(vtkPoints, tmpPts0);
  tmpPts0->InsertNextPoint(this->WorkPd->GetPoint(slicePoints->GetId(0)));
  tmpPts0->InsertNextPoint(this->WorkPd->GetPoint(slicePoints->GetId(1)));
  double sliceFirst[3];
  vtkSVGeneralUtils::GetCentroidOfPoints(tmpPts0, sliceFirst);
  vtkNew(vtkPoints, tmpPts1);
  tmpPts1->InsertNextPoint(this->WorkPd->GetPoint(slicePoints->GetId(3)));
  tmpPts1->InsertNextPoint(this->WorkPd->GetPoint(slicePoints->GetId(2)));
  double sliceSecond[3];
  vtkSVGeneralUtils::GetCentroidOfPoints(tmpPts1, sliceSecond);

  double vec2[3], vec3[3];
  vtkMath::Subtract(sliceFirst, centroid, vec2);
  vtkMath::Subtract(sliceSecond, centroid, vec3);
  double dot0 = vtkMath::Dot(vec2, normal);
  double dot1 = vtkMath::Dot(vec3, normal2);

  vtkMath::MultiplyScalar(normal, dot0);
  vtkMath::MultiplyScalar(normal2, dot1);

  double newSlicePt0[3], newSlicePt1[3];
  double adjSlicePt0[3], adjSlicePt1[3];
  vtkMath::Add(centroid, normal, newSlicePt0);
  vtkMath::MultiplyScalar(normal, 1.0);
  vtkMath::Add(centroid, normal, adjSlicePt0);
  vtkMath::Add(centroid, normal2, newSlicePt1);
  vtkMath::MultiplyScalar(normal2, 1.0);
  vtkMath::Add(centroid, normal2, adjSlicePt1);

  double pt2[3], pt3[3];
  tmpPts0->GetPoint(0, pt2);
  tmpPts0->GetPoint(1, pt3);

  double vec4[3], vec5[3];
  vtkMath::Subtract(pt2, newSlicePt0, vec4);
  vtkMath::Subtract(pt3, newSlicePt0, vec5);

  vtkNew(vtkDoubleArray, boxCutNormals);
  boxCutNormals->SetNumberOfComponents(3);
  boxCutNormals->SetNumberOfTuples(6);
  vtkNew(vtkPoints, boxCutPoints);
  boxCutPoints->SetNumberOfPoints(6);
  double newSliceNormal0[3];
  vtkMath::Cross(vec4, vec5, newSliceNormal0);
  vtkMath::Normalize(newSliceNormal0);

  if (vtkMath::Dot(newSliceNormal0, normal) < 0.0)
    vtkMath::MultiplyScalar(newSliceNormal0, -1.0);

  vtkNew(vtkPlane, cutPlane0);
  cutPlane0->SetOrigin(adjSlicePt0);
  cutPlane0->SetNormal(newSliceNormal0);
  boxCutNormals->SetTuple(0, newSliceNormal0);
  boxCutPoints->SetPoint(0, adjSlicePt0);
  vtkNew(vtkCutter, sliceThrough0);
  sliceThrough0->SetInputData(separator->GetOutput());
  sliceThrough0->SetCutFunction(cutPlane0);
  sliceThrough0->Update();
  vtkNew(vtkConnectivityFilter, getClose0);
  getClose0->SetInputData(sliceThrough0->GetOutput());
  getClose0->SetExtractionModeToClosestPointRegion();
  getClose0->SetClosestPoint(adjSlicePt0);
  getClose0->Update();

  double maxDist = 0.0;
  for (int i=0; i<getClose0->GetOutput()->GetNumberOfPoints(); i++)
  {
    double testPt[3];
    getClose0->GetOutput()->GetPoint(i, testPt);
    double dist = vtkSVGeneralUtils::Distance(adjSlicePt0, testPt);
    if (dist > maxDist)
      maxDist = dist;

  }

  double pt4[3], pt5[3];
  tmpPts1->GetPoint(0, pt4);
  tmpPts1->GetPoint(1, pt5);

  double vec6[3], vec7[3];
  vtkMath::Subtract(pt4, newSlicePt1, vec6);
  vtkMath::Subtract(pt5, newSlicePt1, vec7);

  double newSliceNormal1[3];
  vtkMath::Cross(vec6, vec7, newSliceNormal1);
  vtkMath::Normalize(newSliceNormal1);

  double xvec[3];
  vtkMath::Subtract(sliceFirst, newSlicePt0, xvec);
  //this->Polycube->GetCellData()->GetArray("TopNormal")->InsertTuple(segmentId, newSliceNormal0);

  if (vtkMath::Dot(newSliceNormal1, normal2) < 0.0)
    vtkMath::MultiplyScalar(newSliceNormal1, -1.0);

  vtkNew(vtkPlane, cutPlane1);
  cutPlane1->SetOrigin(adjSlicePt1);
  cutPlane1->SetNormal(newSliceNormal1);
  boxCutNormals->SetTuple(1, newSliceNormal1);
  boxCutPoints->SetPoint(1, adjSlicePt1);
  vtkNew(vtkCutter, sliceThrough1);
  sliceThrough1->SetInputData(separator->GetOutput());
  sliceThrough1->SetCutFunction(cutPlane1);
  sliceThrough1->Update();
  vtkNew(vtkConnectivityFilter, getClose1);
  getClose1->SetInputData(sliceThrough1->GetOutput());
  getClose1->SetExtractionModeToClosestPointRegion();
  getClose1->SetClosestPoint(adjSlicePt1);
  getClose1->Update();

  for (int i=0; i<getClose1->GetOutput()->GetNumberOfPoints(); i++)
  {
    double testPt[3];
    getClose1->GetOutput()->GetPoint(i, testPt);
    double dist = vtkSVGeneralUtils::Distance(adjSlicePt1, testPt);
    if (dist > maxDist)
      maxDist = dist;

  }
  maxDist = maxDist+0.5*maxDist;
  double outPlane0[3], outPlane1[3], outPlane2[3], outPlane3[3];
  double outPoint0[3], outPoint1[3], outPoint2[3], outPoint3[3];
  vtkMath::Add(vec4, vec5, outPlane0);
  vtkMath::Add(vec4, vec5, outPlane1);
  vtkMath::MultiplyScalar(outPlane0, 0.5);
  vtkMath::MultiplyScalar(outPlane1, -0.5);
  vtkMath::Normalize(outPlane0);
  vtkMath::Normalize(outPlane1);
  vtkMath::MultiplyScalar(outPlane0, maxDist);
  vtkMath::MultiplyScalar(outPlane1, maxDist);
  vtkMath::Add(adjSlicePt0, outPlane0, outPoint0);
  vtkMath::Add(adjSlicePt0, outPlane1, outPoint1);
  boxCutNormals->SetTuple(2, outPlane0);
  boxCutPoints->SetPoint(2, outPoint0);
  boxCutNormals->SetTuple(3, outPlane1);
  boxCutPoints->SetPoint(3, outPoint1);
  vtkMath::Cross(newSliceNormal0, outPlane0, outPlane2);
  vtkMath::Cross(newSliceNormal0, outPlane1, outPlane3);
  vtkMath::Normalize(outPlane2);
  vtkMath::Normalize(outPlane3);
  vtkMath::MultiplyScalar(outPlane2, maxDist);
  vtkMath::MultiplyScalar(outPlane3, maxDist);
  vtkMath::Add(adjSlicePt0, outPlane2, outPoint2);
  vtkMath::Add(adjSlicePt0, outPlane3, outPoint3);
  boxCutNormals->SetTuple(4, outPlane2);
  boxCutPoints->SetPoint(4, outPoint2);
  boxCutNormals->SetTuple(5, outPlane3);
  boxCutPoints->SetPoint(5, outPoint3);
  vtkNew(vtkPlanes, cutPlanes);
  cutPlanes->SetPoints(boxCutPoints);
  cutPlanes->SetNormals(boxCutNormals);

  fprintf(stdout,"Cutting first at points: %lld %lld\n", slicePoints->GetId(0), slicePoints->GetId(1));
  fprintf(stdout,"Cutting second at points: %lld %lld\n", slicePoints->GetId(2), slicePoints->GetId(3));
  vtkNew(vtkPolyData, slicePd0);
  vtkNew(vtkPolyData, slicePd1);
  vtkNew(vtkPolyData, leftovers0);
  vtkNew(vtkPolyData, leftovers1);
  vtkSVGeneralUtils::ClipCut(separator->GetOutput(), cutPlane0, 1, 1, slicePd0, leftovers0);
  vtkSVGeneralUtils::ClipCut(slicePd0, cutPlane1, 1, 1, slicePd1, leftovers1);

  vtkNew(vtkIdFilter, ider3);
  ider3->SetInputData(slicePd1);
  ider3->SetIdsArrayName(this->InternalIdsArrayName);
  ider3->Update();
  slicePd1->DeepCopy(ider3->GetOutput());
  vtkNew(vtkPolyData, onlyGood);
  vtkSVGeneralUtils::GetClosestPointConnectedRegion(slicePd1, centroid, onlyGood);
  int numCells = onlyGood->GetNumberOfCells();
  for (int i=0; i<numCells; i++)
  {
    slicePd1->GetCellData()->GetArray(this->SegmentIdsArrayName)->
      InsertTuple1(onlyGood->GetCellData()->GetArray(this->InternalIdsArrayName)->
          GetTuple1(i), segmentId);
  }

  vtkNew(vtkAppendPolyData, appender2);
  appender2->AddInputData(theRestPd);
  appender2->AddInputData(section2Pd);
  appender2->AddInputData(leftovers0);
  appender2->AddInputData(leftovers1);
  appender2->AddInputData(slicePd1);
  appender2->Update();

  vtkNew(vtkCleanPolyData, cleaner2);
  cleaner2->SetInputData(appender2->GetOutput());
  cleaner2->Update();
  pd->DeepCopy(cleaner2->GetOutput());
  vtkNew(vtkPolyData, dumTemp);
  dumTemp->SetPoints(boxCutPoints);
  boxCutNormals->SetName("DumBNorms");
  dumTemp->GetPointData()->AddArray(boxCutNormals);

  vtkNew(vtkIdFilter, ider2);
  ider2->SetInputData(pd);
  ider2->SetIdsArrayName(this->InternalIdsArrayName);
  ider2->Update();
  pd->DeepCopy(ider2->GetOutput());

  vtkNew(vtkPointLocator, locator);
  locator->SetDataSet(pd);
  locator->BuildLocator();
  int frontId0 = locator->FindClosestPoint(tmpPts0->GetPoint(0));
  int backId0  = locator->FindClosestPoint(tmpPts0->GetPoint(1));

  vtkNew(vtkIdList, fixedSurgeryPoints0);
  vtkNew(vtkIdList, fixedSlicePoints0);
  this->CriticalSurgeryPoints(pd, frontId0, backId0, parentId, segmentId,
                             centroid, adjSlicePt0,
                             fixedSlicePoints0, fixedSurgeryPoints0, xvec);
  fprintf(stdout,"TEST TRY: %lld %lld %lld %lld\n", fixedSlicePoints0->GetId(0),
                                                    fixedSlicePoints0->GetId(1),
                                                    fixedSurgeryPoints0->GetId(1),
                                                    fixedSurgeryPoints0->GetId(0));

  vtkIntArray *cornerIds = vtkIntArray::SafeDownCast(
    pd->GetPointData()->GetArray("SurgeryPoints"));
  vtkNew(vtkIntArray, singleCompArray0);
  singleCompArray0->SetNumberOfComponents(1);
  singleCompArray0->SetNumberOfTuples(pd->GetNumberOfPoints());
  singleCompArray0->CopyComponent(0, pd->GetPointData()->GetArray("SurgeryPoints"), parentId);

  int parIndices[8];
  for (int i=0; i<8; i++)
    parIndices[i] = this->LookupIndex(gCell->Dir, gCell->Children[gCell->DivergingChild]->Dir, i);

  cornerIds->SetComponent(fixedSlicePoints0->GetId(0),   parentId, parIndices[4]);
  cornerIds->SetComponent(fixedSlicePoints0->GetId(1),   parentId, parIndices[5]);
  cornerIds->SetComponent(fixedSurgeryPoints0->GetId(1), parentId, parIndices[6]);
  cornerIds->SetComponent(fixedSurgeryPoints0->GetId(0), parentId, parIndices[7]);

  int frontId1 = locator->FindClosestPoint(tmpPts1->GetPoint(0));
  int backId1  = locator->FindClosestPoint(tmpPts1->GetPoint(1));

  vtkNew(vtkIdList, fixedSurgeryPoints1);
  vtkNew(vtkIdList, fixedSlicePoints1);
  this->CriticalSurgeryPoints(pd, frontId1, backId1, goodKidId, segmentId,
                             centroid, adjSlicePt1,
                             fixedSlicePoints1, fixedSurgeryPoints1, dummy);
  pd->GetPointData()->RemoveArray(this->InternalIdsArrayName);
  pd->GetCellData()->RemoveArray(this->InternalIdsArrayName);
  fprintf(stdout,"TEST TRY: %lld %lld %lld %lld\n", fixedSlicePoints1->GetId(0),
                                                    fixedSlicePoints1->GetId(1),
                                                    fixedSurgeryPoints1->GetId(1),
                                                    fixedSurgeryPoints1->GetId(0));
  vtkNew(vtkIntArray, singleCompArray1);
  singleCompArray1->SetNumberOfComponents(1);
  singleCompArray1->SetNumberOfTuples(pd->GetNumberOfPoints());
  singleCompArray1->CopyComponent(0, pd->GetPointData()->GetArray("SurgeryPoints"), goodKidId);

  cornerIds->SetComponent(fixedSlicePoints1->GetId(0),   goodKidId, parIndices[0]);
  cornerIds->SetComponent(fixedSlicePoints1->GetId(1),   goodKidId, parIndices[1]);
  cornerIds->SetComponent(fixedSurgeryPoints1->GetId(1), goodKidId, parIndices[2]);
  cornerIds->SetComponent(fixedSurgeryPoints1->GetId(0), goodKidId, parIndices[3]);

  fprintf(stdout,"TEST TRY: %lld %lld %lld %lld\n", fixedSlicePoints0->GetId(1),
                                                    fixedSlicePoints0->GetId(0),
                                                    fixedSlicePoints1->GetId(1),
                                                    fixedSlicePoints1->GetId(0));
  vtkNew(vtkIntArray, singleCompArray2);
  singleCompArray2->SetNumberOfComponents(1);
  singleCompArray2->SetNumberOfTuples(pd->GetNumberOfPoints());
  singleCompArray2->CopyComponent(0, pd->GetPointData()->GetArray("SurgeryPoints"), badKidId);

  cornerIds->SetComponent(fixedSlicePoints0->GetId(0), badKidId, parIndices[3]);
  cornerIds->SetComponent(fixedSlicePoints0->GetId(1), badKidId, parIndices[2]);
  cornerIds->SetComponent(fixedSlicePoints1->GetId(1), badKidId, parIndices[6]);
  cornerIds->SetComponent(fixedSlicePoints1->GetId(0), badKidId, parIndices[7]);

  fprintf(stdout,"TEST TRY: %lld %lld %lld %lld %lld %lld %lld %lld\n", fixedSlicePoints0->GetId(0),
                                                                        fixedSlicePoints0->GetId(1),
                                                                        fixedSurgeryPoints0->GetId(1),
                                                                        fixedSurgeryPoints0->GetId(0),
                                                                        fixedSlicePoints1->GetId(0),
                                                                        fixedSlicePoints1->GetId(1),
                                                                        fixedSurgeryPoints1->GetId(1),
                                                                        fixedSurgeryPoints1->GetId(0));

  //double segmentCornerPts[8];
  //segmentCornerPts[0] = fixedSlicePoints->GetId(0),
  //segmentCornerPts[1] = fixedSurgeryPoints0->GetId(0),
  //segmentCornerPts[2] = fixedSurgeryPoints0->GetId(1),
  //segmentCornerPts[3] = fixedSlicePoints->GetId(1),
  //segmentCornerPts[4] = fixedSlicePoints->GetId(2),
  //segmentCornerPts[5] = fixedSurgeryPoints1->GetId(0),
  //segmentCornerPts[6] = fixedSurgeryPoints1->GetId(1),
  //segmentCornerPts[7] = fixedSlicePoints->GetId(3);
  //cornerIds->InsertTuple(segmentId, segmentCornerPts);
  cornerIds->SetComponent(fixedSlicePoints0->GetId(0),   segmentId, parIndices[0]);
  cornerIds->SetComponent(fixedSlicePoints0->GetId(1),   segmentId, parIndices[1]);
  cornerIds->SetComponent(fixedSurgeryPoints0->GetId(1), segmentId, parIndices[2]);
  cornerIds->SetComponent(fixedSurgeryPoints0->GetId(0), segmentId, parIndices[3]);
  cornerIds->SetComponent(fixedSlicePoints1->GetId(0),   segmentId, parIndices[4]);
  cornerIds->SetComponent(fixedSlicePoints1->GetId(1),   segmentId, parIndices[5]);
  cornerIds->SetComponent(fixedSurgeryPoints1->GetId(1), segmentId, parIndices[6]);
  cornerIds->SetComponent(fixedSurgeryPoints1->GetId(0), segmentId, parIndices[7]);

  fprintf(stdout,"SLICE POINTS 0: %lld %lld\n", fixedSlicePoints0->GetId(0),
                                                fixedSlicePoints0->GetId(1));
  fprintf(stdout,"SLICE POINTS 1: %lld %lld\n", fixedSlicePoints1->GetId(0),
                                                fixedSlicePoints1->GetId(1));
  fprintf(stdout,"SURGE POINTS 0: %lld %lld\n", fixedSurgeryPoints0->GetId(0),
                                                fixedSurgeryPoints0->GetId(1));
  fprintf(stdout,"SURGE POINTS 1: %lld %lld\n", fixedSurgeryPoints1->GetId(0),
                                                fixedSurgeryPoints1->GetId(1));
  fprintf(stdout,"PAR CHANGE: ");
  for (int i=0; i<8; i++)
    fprintf(stdout,"%d ",parIndices[i]);
  fprintf(stdout,"\n");

  return SV_OK;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetCorrectFrontPoint(vtkPolyData *pd,
                                                        double frontDir[3],
                                                        int &frontId,
                                                        int &backId)
{
  double frontPt[3], backPt[3];
  pd->GetPoint(frontId, frontPt);
  pd->GetPoint(backId, backPt);
  double checkVec[3];
  vtkMath::Subtract(frontPt, backPt, checkVec);
  vtkMath::Normalize(checkVec);
  double dotCheck = vtkMath::Dot(checkVec, frontDir);
    //fprintf(stdout,"What the eff is the dot! %.4f\n", vtkMath::Dot(checkVec, frontDir));
  if (dotCheck < 0)
  {
    fprintf(stdout,"Wrong front point!\n");
    int tmp = frontId;
    frontId = backId;
    backId  = tmp;
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetFourPolyDataRegions(vtkPolyData *startPd,
                                                          const int id0,
                                                          vtkPolyData *pd0,
                                                          const int id1,
                                                          vtkPolyData *pd1,
                                                          const int id2,
                                                          vtkPolyData *pd2,
                                                          vtkPolyData *leftovers)
{
  int large = 1000;
  vtkNew(vtkPolyData, tmpPd0);
  vtkNew(vtkPolyData, tmpPd1);
  vtkNew(vtkAppendPolyData, appender0);
  vtkSVGeneralUtils::ThresholdPd(startPd, id0,  id0,  1,
                                             this->GroupIdsArrayName, pd0);
  if (vtkSVGeneralUtils::ThresholdPd(startPd, -1,  id0-1,  1,
                                                 this->GroupIdsArrayName, tmpPd0) != 0)
  {
    appender0->AddInputData(tmpPd0);
  }
  if (vtkSVGeneralUtils::ThresholdPd(startPd, id0+1, large,  1,
                                                 this->GroupIdsArrayName, tmpPd1) != 0)
  {
    appender0->AddInputData(tmpPd1);
  }
  appender0->Update();
  leftovers->DeepCopy(appender0->GetOutput());

  vtkNew(vtkAppendPolyData, appender1);
  vtkSVGeneralUtils::ThresholdPd(startPd, id1,  id1,  1,
                                             this->GroupIdsArrayName, pd1);
  if (vtkSVGeneralUtils::ThresholdPd(leftovers, -1,  id1-1,  1,
                                                 this->GroupIdsArrayName, tmpPd0) != 0)
  {
    appender1->AddInputData(tmpPd0);
  }
  if (vtkSVGeneralUtils::ThresholdPd(leftovers, id1+1, large,  1,
                                                 this->GroupIdsArrayName, tmpPd1) != 0)
  {
    appender1->AddInputData(tmpPd1);
  }
  appender1->Update();
  leftovers->DeepCopy(appender1->GetOutput());

  vtkNew(vtkAppendPolyData, appender2);
  vtkSVGeneralUtils::ThresholdPd(startPd, id2, id2, 1,
                                             this->GroupIdsArrayName, pd2);
  if (vtkSVGeneralUtils::ThresholdPd(leftovers, -1, id2-1, 1,
                                                 this->GroupIdsArrayName, tmpPd0) != 0)
  {
    appender2->AddInputData(tmpPd0);
  }
  if (vtkSVGeneralUtils::ThresholdPd(leftovers, id2+1, large, 1,
                                                 this->GroupIdsArrayName, tmpPd1) != 0)
  {
    appender2->AddInputData(tmpPd1);
  }
  appender2->Update();
  leftovers->DeepCopy(appender2->GetOutput());
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::CriticalSurgeryPoints(vtkPolyData *pd,
                                                         const int frontId,
                                                         const int backId,
                                                         const int groupId,
                                                         const int checkId,
                                                         double startPt[3],
                                                         double secondPt[3],
                                                         vtkIdList *fixedSlicePoints,
                                                         vtkIdList *fixedSurgeryPoints,
                                                         double startDir[3])
{
  fixedSlicePoints->SetNumberOfIds(2);
  fixedSlicePoints->SetId(0, frontId);
  fixedSlicePoints->SetId(1, backId);
  fprintf(stdout,"Front is %d and back is %d\n", frontId, backId);

  vtkNew(vtkPolyData, thresholdPd);
  vtkSVGeneralUtils::ThresholdPd(pd, groupId, groupId, 1,
      this->SegmentIdsArrayName, thresholdPd);

  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(thresholdPd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();
  vtkNew(vtkPolyData, thresholdLoop);
  vtkSVGeneralUtils::GetClosestPointConnectedRegion(boundaries->GetOutput(), startPt, thresholdLoop);

  vtkDataArray *thresholdIds =
    thresholdLoop->GetPointData()->GetArray(this->InternalIdsArrayName);
  int fixFront = thresholdIds->LookupValue(frontId);
  int fixBack = thresholdIds->LookupValue(backId);
  this->GetSurgeryPoints(thresholdLoop, pd, thresholdIds, startPt, secondPt, fixFront, fixBack, checkId, this->SegmentIdsArrayName, fixedSurgeryPoints, startDir);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::BuildPolycube()
{
  int numCubes = this->CenterlineGraph->NumberOfNodes;
  fprintf(stdout,"Number of cubes: %d\n", numCubes);
  this->Polycube->SetNumberOfGrids(numCubes);

  double dims[3]; dims[0] = 0.5; dims[1] = 0.5; dims[2] = 0.5;
  double startPt[3];
  int parentDir = this->CenterlineGraph->Root->Dir;
  int childDir  = this->CenterlineGraph->Root->Children[
    this->CenterlineGraph->Root->DivergingChild]->Dir;
  vtkMath::Add(this->CenterlineGraph->Root->StartPt,
               this->CenterlineGraph->Root->EndPt, startPt);
  vtkMath::MultiplyScalar(startPt, 0.5);
  this->Polycube->SetGridWithCenter(this->CenterlineGraph->Root->GroupId,
                                    startPt,
                                    dims, vtkSVGeneralizedPolycube::CUBE_BRANCH,
                                    parentDir,
                                    childDir);
  for (int i=0; i<8; i++)
    this->Polycube->GetCellData()->GetArray("CornerPtIds")->SetComponent(0, i, this->CenterlineGraph->Root->CornerPtIds[i]);
  this->Polycube->GetCellData()->GetArray("TopNormal")->SetTuple(0, this->CenterlineGraph->Root->ZVec);
  this->Polycube->GetCellData()->GetArray("RightNormal")->SetTuple(0, this->CenterlineGraph->Root->XVec);
  vtkNew(vtkIntArray, segmentIds);
  segmentIds->InsertTuple1(this->CenterlineGraph->Root->GroupId,
                           this->CenterlineGraph->Root->GroupId);
  this->CenterlineGraph->Recurse(this->CenterlineGraph->Root,
                                 vtkSVPolyDataSliceAndDiceFilter::GraphToPolycube,
                                 this->Polycube,
                                 segmentIds, NULL);

  segmentIds->SetName(this->SegmentIdsArrayName);
  this->Polycube->GetCellData()->AddArray(segmentIds);

  vtkNew(vtkAppendFilter, merger);
  merger->SetInputData(this->Polycube);
  merger->MergePointsOn();
  merger->Update();

  this->Polycube->DeepCopy(merger->GetOutput());
  this->Polycube->BuildLinks();

  return SV_OK;
}

// ----------------------
// GraphToPolycube
// ----------------------
/**
 * \details If both children of the node are not NULL, then four cubes will be created.
 * One for the given node, one for the linking patch, and one for each child node.
 */
int vtkSVPolyDataSliceAndDiceFilter::GraphToPolycube(svGCell *gCell, void *arg0,
                                                   void *arg1, void *arg2)
{
  if (gCell->Children[0] != NULL && gCell->Children[1] != NULL)
  {
    vtkSVGeneralizedPolycube *polycube =
      reinterpret_cast<vtkSVGeneralizedPolycube*>(arg0);
    vtkIntArray *segmentIds =
      reinterpret_cast<vtkIntArray*>(arg1);
    int id = gCell->GroupId + 1;
    double dims[3]; dims[0] = 0.5; dims[1] = 0.5; dims[2] = 0.5;
    polycube->SetGridWithCenter(id, gCell->EndPt,
                                dims, vtkSVGeneralizedPolycube::CUBE_BIFURCATION,
                                gCell->Dir,
                                gCell->Children[gCell->DivergingChild]->Dir);
    int ptIndices[8];
    for (int i=0; i<8; i++)
      ptIndices[i] = vtkSVPolyDataSliceAndDiceFilter::LookupIndex(gCell->Dir, gCell->Children[gCell->DivergingChild]->Dir, i);
    for (int i=0; i<4; i++)
    {
      polycube->GetCellData()->GetArray("CornerPtIds")->SetComponent(id, ptIndices[i], gCell->CornerPtIds[ptIndices[i+4]]);
      polycube->GetCellData()->GetArray("CornerPtIds")->SetComponent(id, ptIndices[i+4], gCell->Children[gCell->AligningChild]->CornerPtIds[ptIndices[i]]);
    }
    polycube->GetCellData()->GetArray("TopNormal")->SetTuple(id, gCell->ZVec);
    double negXVec[3];
    for (int i=0; i<3; i++)
      negXVec[i] = -1.0*gCell->XVec[i];
    polycube->GetCellData()->GetArray("RightNormal")->SetTuple(id, negXVec);
    segmentIds->InsertTuple1(id, id);

    int defaultDirs[6];
    defaultDirs[DOWN]  = RIGHT;
    defaultDirs[UP]    = LEFT;
    defaultDirs[RIGHT] = UP;
    defaultDirs[LEFT]  = DOWN;
    defaultDirs[BACK]  = RIGHT;
    defaultDirs[FRONT] = RIGHT;
    //Two children
    for (int i=0; i<2; i++)
    {
      int newDir;
      if (gCell->Children[i]->Children[0] != NULL)
        newDir = gCell->Children[i]->Children[gCell->Children[i]->DivergingChild]->Dir;
      else
        newDir = defaultDirs[gCell->Children[i]->Dir];
      id = gCell->Children[i]->GroupId;
      double centerPt[3];
      vtkMath::Add(gCell->Children[i]->StartPt,
                   gCell->Children[i]->EndPt, centerPt);
      vtkMath::MultiplyScalar(centerPt, 0.5);
      polycube->SetGridWithCenter(id, centerPt,
                                  dims, vtkSVGeneralizedPolycube::CUBE_BRANCH,
                                  gCell->Children[i]->Dir,
                                  newDir);
      for (int j=0; j<8; j++)
        polycube->GetCellData()->GetArray("CornerPtIds")->SetComponent(id, j, gCell->Children[i]->CornerPtIds[j]);
      polycube->GetCellData()->GetArray("TopNormal")->SetTuple(id, gCell->Children[i]->ZVec);
      polycube->GetCellData()->GetArray("RightNormal")->SetTuple(id, gCell->Children[i]->XVec);
      segmentIds->InsertTuple1(id, id);
    }

  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetNextSurgeryPoints(vtkPolyData *pd, double centerPt[3], vtkIdList *surgeryPoints, double xvec[3], double zvec[3], double radius, vtkIdList *surgeryLineIds)
{
  int contourPtId;
  vtkNew(vtkPolyData, boundary);

  int initialSurgeryPt = surgeryPoints->GetId(0);
  this->GetClose3DPoint(pd, centerPt, initialSurgeryPt,
    contourPtId, xvec, zvec, radius, boundary);

  //fprintf(stdout,"Ending point to use is : %d\n", contourPtId);
  surgeryPoints->SetId(0, contourPtId);

  int pointId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(contourPtId);

  vtkNew(vtkIdList, startCellIds);
  boundary->BuildLinks();
  boundary->GetPointCells(pointId, startCellIds);
  int cellId = startCellIds->GetId(0);

  vtkIdType npts, *pts;
  boundary->GetCellPoints(cellId, npts, pts);
  int secondPtId;
  if (pts[0] == pointId)
  {
    secondPtId = pts[1];
  }
  else
  {
    secondPtId = pts[0];
  }
  double pt0[3], pt1[3];
  boundary->GetPoint(pointId, pt0);
  boundary->GetPoint(secondPtId, pt1);
  double vec0[3], vec1[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  vtkMath::Cross(xvec, vec0, vec1);

  if (vtkMath::Dot(zvec, vec1) > 0)
  {
    cellId = startCellIds->GetId(1);
  }

  //Getting full loop length
  int front = pointId;
  int back = pointId;
  double length = 0.0;
  int iter = 0;
  int prevCellId = cellId;
  while (pointId != back || iter == 0)
  {
    boundary->GetPoint(pointId, pt0);
    vtkSVGeneralUtils::IteratePoint(boundary, pointId, cellId);
    boundary->GetPoint(pointId, pt1);
    length += vtkSVGeneralUtils::Distance(pt0, pt1);
    iter++;
  }

  //Finding the surgery points which separate boundary into four points
  double surgeryCount = 1.0;
  double currLength = 0.0;;
  back = pointId;
  pointId = front;
  prevCellId = cellId;
  iter = 0;
  while (pointId != back || iter == 0)
  {
    boundary->GetPoint(pointId, pt0);
    vtkSVGeneralUtils::IteratePoint(boundary, pointId, prevCellId);
    boundary->GetPoint(pointId, pt1);
    currLength += vtkSVGeneralUtils::Distance(pt0, pt1);
    if (currLength > length*(surgeryCount/4.0))
    {
      if (surgeryCount < 4)
      {
        int newId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
          GetTuple1(pointId);
        surgeryPoints->SetId(surgeryCount, newId);
      }
      surgeryCount += 1.0;
    }
    iter++;
  }

  vtkNew(vtkSVFindGeodesicPath, finder);
  finder->SetInputData(pd);
  finder->SetStartPtId(pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(surgeryPoints->GetId(0)));
  finder->SetEndPtId(pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(initialSurgeryPt));
  finder->SetDijkstraArrayName(this->DijkstraArrayName);
  finder->SetInternalIdsArrayName(this->InternalIdsArrayName);
  finder->SetRepelCloseBoundaryPoints(1);
  finder->Update();

  vtkNew(vtkIdList, tmpIds);
  tmpIds = finder->GetPathIds();
  int numToAdd = tmpIds->GetNumberOfIds();
  for (int i=0; i<numToAdd; i++)
  {
    int testId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
      GetTuple1(tmpIds->GetId(i));
    if (surgeryLineIds->IsId(testId) == -1)
    {
      surgeryLineIds->InsertNextId(testId);
    }
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetEndSurgeryPoints(vtkPolyData *pd, svGCell *gCell, double centerPt[3], vtkIdList *surgeryPoints, int endSurgeryIds[8], double xvec[3], double zvec[3], double radius, vtkIdList *surgeryLineIds, int cellIndices[8], int &secondRun)
{
  int contourPtId;
  vtkNew(vtkPolyData, boundary);
  int initialSurgeryPt = surgeryPoints->GetId(0);
  int finalSurgeryId;
  if (this->SliceDirection == 1)
    finalSurgeryId = endSurgeryIds[cellIndices[0]];
  else
    finalSurgeryId = endSurgeryIds[cellIndices[4]];

  int finally = 0;
  if (finalSurgeryId == -1)
  {
    this->GetClose3DPoint(pd, centerPt, surgeryPoints->GetId(0),
      contourPtId, xvec, zvec, radius, boundary);
    //fprintf(stdout,"Ending point to use is : %d\n", contourPtId);
    surgeryPoints->SetId(0, contourPtId);
  }
  else
  {
    finally = 1;
    vtkNew(vtkFeatureEdges, boundaries);
    boundaries->SetInputData(pd);
    boundaries->BoundaryEdgesOn();
    boundaries->FeatureEdgesOff();
    boundaries->NonManifoldEdgesOff();
    boundaries->ManifoldEdgesOff();
    boundaries->Update();
    vtkSVGeneralUtils::GetClosestPointConnectedRegion(boundaries->GetOutput(), centerPt, boundary);

    int iContour = 0;
    if (secondRun)
      contourPtId = endSurgeryIds[cellIndices[4]];
    else
    {
      double projVec[3];
      for (int i=0; i<3; i++)
      {
        projVec[i] = xvec[i];
      }
      vtkMath::MultiplyScalar(projVec, 1.5*radius);
      double closePt[3];
      vtkMath::Add(centerPt, projVec, closePt);
      double dist = 1.0e299;
      for (int i=0; i<4; i++)
      {
        int testId = boundary->GetPointData()->GetArray(
          this->InternalIdsArrayName)->LookupValue(endSurgeryIds[cellIndices[i+4]]);
        double pt[3];
        boundary->GetPoint(testId, pt);
        double testDist = vtkSVGeneralUtils::Distance(closePt, pt);
        if (testDist < dist)
        {
          dist = testDist;
          contourPtId = endSurgeryIds[cellIndices[i+4]];
          iContour = i;
        }
      }
    }
    fprintf(stdout,"Closest point is!: %d\n", contourPtId);
    fprintf(stdout,"Supposed to be: %d\n", endSurgeryIds[cellIndices[4]]);
    if (contourPtId != endSurgeryIds[cellIndices[4]])
    {
      this->FixGraphDirections(gCell, contourPtId, cellIndices);
      for (int i=0; i<4; i++)
        surgeryPoints->SetId(i, gCell->CornerPtIds[cellIndices[i+4]]);
      secondRun = 1;
      return SV_OK;
    }
    else
    {
      for (int i=0; i<4; i++)
        surgeryPoints->SetId(i, gCell->CornerPtIds[cellIndices[i+4]]);
      secondRun = 0;
    }
  }
  fprintf(stdout,"And the first surg point is: %d\n", surgeryPoints->GetId(0));

  if (!finally)
  {
    int pointId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
      LookupValue(contourPtId);

    vtkNew(vtkIdList, startCellIds);
    boundary->BuildLinks();
    boundary->GetPointCells(pointId, startCellIds);
    int cellId = startCellIds->GetId(0);

    vtkIdType npts, *pts;
    boundary->GetCellPoints(cellId, npts, pts);
    int secondPtId;
    if (pts[0] == pointId)
    {
      secondPtId = pts[1];
    }
    else
    {
      secondPtId = pts[0];
    }
    double pt0[3], pt1[3];
    boundary->GetPoint(pointId, pt0);
    boundary->GetPoint(secondPtId, pt1);
    double tmpXVec[3], vec0[3], vec1[3];
    vtkMath::Subtract(pt1, pt0, vec0);
    vtkMath::Normalize(vec0);
    vtkMath::Subtract(pt0, centerPt, tmpXVec);
    vtkMath::Normalize(tmpXVec);

    vtkMath::Cross(tmpXVec, vec0, vec1);

    if (vtkMath::Dot(zvec, vec1) > 0)
    {
      cellId = startCellIds->GetId(1);
    }

    //Getting full loop length
    int front = pointId;
    int back = pointId;
    double length = 0.0;
    int iter = 0;
    int prevCellId = cellId;
    while (pointId != back || iter == 0)
    {
      boundary->GetPoint(pointId, pt0);
      vtkSVGeneralUtils::IteratePoint(boundary, pointId, cellId);
      boundary->GetPoint(pointId, pt1);
      length += vtkSVGeneralUtils::Distance(pt0, pt1);
      iter++;
    }

    //Finding the surgery points which separate boundary into four points
    double surgeryCount = 1.0;
    double currLength = 0.0;;
    back = pointId;
    pointId = front;
    prevCellId = cellId;
    iter = 0;
    while (pointId != back || iter == 0)
    {
      boundary->GetPoint(pointId, pt0);
      vtkSVGeneralUtils::IteratePoint(boundary, pointId, prevCellId);
      boundary->GetPoint(pointId, pt1);
      currLength += vtkSVGeneralUtils::Distance(pt0, pt1);
      if (currLength > length*(surgeryCount/4.0))
      {
        if (surgeryCount < 4)
        {
          int newId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
            GetTuple1(pointId);
          surgeryPoints->SetId(surgeryCount, newId);
        }
        surgeryCount += 1.0;
      }
      iter++;
    }
  }
  vtkNew(vtkSVFindGeodesicPath, finder);
  finder->SetInputData(pd);
  finder->SetStartPtId(pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(surgeryPoints->GetId(0)));
  finder->SetEndPtId(pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    LookupValue(initialSurgeryPt));
  finder->SetDijkstraArrayName(this->DijkstraArrayName);
  finder->SetInternalIdsArrayName(this->InternalIdsArrayName);
  finder->SetRepelCloseBoundaryPoints(1);
  finder->Update();

  vtkNew(vtkIdList, tmpIds);
  tmpIds = finder->GetPathIds();
  int numToAdd = tmpIds->GetNumberOfIds();
  for (int i=0; i<numToAdd; i++)
  {
    int testId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
      GetTuple1(tmpIds->GetId(i));
    if (surgeryLineIds->IsId(testId) == -1)
    {
      surgeryLineIds->InsertNextId(testId);
    }
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::FixGraphDirections(svGCell *gCell, const int actualId,
                                                      int cellIndices[8])
{
  int pDir = gCell->Dir;
  vtkNew(vtkIntArray, rotIndices);
  rotIndices->SetNumberOfComponents(1);
  rotIndices->SetNumberOfTuples(8);
  if (pDir == RIGHT || pDir == LEFT)
  {
    //Rotate around X 90° until actualId mathces!
    for (int i=0; i<8; i++)
      rotIndices->SetTuple1(i, this->LookupIndex(BACK, RIGHT, i));
  }
  else if (pDir == FRONT || pDir == BACK)
  {
    //Rotate around Y 90° until actualId mathces!
    for (int i=0; i<8; i++)
      rotIndices->SetTuple1(i, this->LookupIndex(LEFT, DOWN, i));
  }
  else if (pDir == UP || pDir == DOWN)
  {
    //Rotate around Z 90° until actualId mathces!
    for (int i=0; i<8; i++)
      rotIndices->SetTuple1(i, this->LookupIndex(DOWN, BACK, i));
  }
  if (gCell->Children[0] == NULL)
    fprintf(stdout, "This should not be able to happen!\n");
  svGCell *child = gCell->Children[gCell->DivergingChild];
  int newId = -1;
  while (newId != actualId)
  {
    fprintf(stdout,"Rotating!!!!\n");
    int tmpIds[8];
    for (int i=0; i<8; i++)
      tmpIds[i] = gCell->CornerPtIds[i];
    for (int i=0; i<4; i++)
    {
      fprintf(stdout,"%d is now becoming %d\n", gCell->CornerPtIds[cellIndices[i+4]], tmpIds[rotIndices->GetValue(cellIndices[i+4])]);
      gCell->CornerPtIds[cellIndices[i+4]] = tmpIds[rotIndices->GetValue(cellIndices[i+4])];
    }
    gCell->Children[gCell->DivergingChild]->RefAngle += M_PI/2.0;
    if (gCell->Children[gCell->DivergingChild]->RefAngle > 2.0*M_PI)
      gCell->Children[gCell->DivergingChild]->RefAngle -= 2.0*M_PI;
    newId = gCell->CornerPtIds[cellIndices[4]];
    svGraph::Recurse(gCell->Children[0], svGraph::UpdateCellDirection,
                     rotIndices, NULL, NULL);
    svGraph::Recurse(gCell->Children[1], svGraph::UpdateCellDirection,
                     rotIndices, NULL, NULL);
  }
  fprintf(stdout,"NewANGe: %.4f\n", gCell->Children[gCell->DivergingChild]->RefAngle);
  //Fix indices
  //for (int i=0; i<8; i++)
  //  cellIndices[i] = this->LookupIndex(gCell->Dir, gCell->Children[gCell->DivergingChild]->Dir, i);
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetFirstSurgeryPoints(vtkPolyData *pd, int pointId, vtkIdList *surgeryPoints, double xvec[3], double zvec[3])
{
  vtkNew(vtkIdList, startCellIds);
  pd->BuildLinks();
  pd->GetPointCells(pointId, startCellIds);
  int cellId = startCellIds->GetId(0);

  vtkIdType npts, *pts;
  pd->GetCellPoints(cellId, npts, pts);
  int secondPtId;
  if (pts[0] == pointId)
  {
    secondPtId = pts[1];
  }
  else
  {
    secondPtId = pts[0];
  }
  double pt0[3], pt1[3];
  pd->GetPoint(pointId, pt0);
  pd->GetPoint(secondPtId, pt1);
  double vec0[3], vec1[3];
  vtkMath::Subtract(pt1, pt0, vec0);
  vtkMath::Normalize(vec0);

  vtkMath::Cross(xvec, vec0, vec1);

  if (vtkMath::Dot(zvec, vec1) > 0)
  {
    cellId = startCellIds->GetId(1);
  }

  //Getting full loop length
  int front = pointId;
  int back = pointId;
  double length = 0.0;
  int iter = 0;
  int prevCellId = cellId;
  while (pointId != back || iter == 0)
  {
    pd->GetPoint(pointId, pt0);
    vtkSVGeneralUtils::IteratePoint(pd, pointId, cellId);
    pd->GetPoint(pointId, pt1);
    length += vtkSVGeneralUtils::Distance(pt0, pt1);
    iter++;
  }

  //Finding the surgery points which separate boundary into four points
  double surgeryCount = 1.0;
  double currLength = 0.0;;
  back = pointId;
  front = pointId;
  prevCellId = cellId;
  iter = 0;
  int newId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
    GetTuple1(pointId);
  surgeryPoints->InsertNextId(newId);
  while (pointId != back || iter == 0)
  {
    pd->GetPoint(pointId, pt0);
    vtkSVGeneralUtils::IteratePoint(pd, pointId, prevCellId);
    pd->GetPoint(pointId, pt1);
    currLength += vtkSVGeneralUtils::Distance(pt0, pt1);
    if (currLength > length*(surgeryCount/4.0))
    {
      if (surgeryCount < 4)
      {
        newId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)->
          GetTuple1(pointId);
        surgeryPoints->InsertNextId(newId);
      }
      surgeryCount += 1.0;
    }
    iter++;
  }

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::GetClose3DPoint(vtkPolyData *pd, double centerPt[3],
                                                   const int startPtId, int &returnStartId,
                                                   double xvec[3], double zvec[3],
                                                   double radius,
                                                   vtkPolyData *boundary)
{
  vtkNew(vtkFeatureEdges, boundaries);
  boundaries->SetInputData(pd);
  boundaries->BoundaryEdgesOn();
  boundaries->FeatureEdgesOff();
  boundaries->NonManifoldEdgesOff();
  boundaries->ManifoldEdgesOff();
  boundaries->Update();

  double projVec[3];
  for (int i=0; i<3; i++)
  {
    projVec[i] = xvec[i];
  }
  vtkMath::MultiplyScalar(projVec, 1.5*radius);
  double closePt[3];
  vtkMath::Add(centerPt, projVec, closePt);
  //fprintf(stdout,"Input start Pt ID: %d\n", startPtId);
  //fprintf(stdout,"Close point is!: %.4f %.4f %.4f\n", closePt[0], closePt[1], closePt[2]);

  vtkSVGeneralUtils::GetClosestPointConnectedRegion(boundaries->GetOutput(), centerPt, boundary);

  int pointId = pd->GetPointData()->GetArray(this->InternalIdsArrayName)
    ->LookupValue(startPtId);
  double startPt[3];
  pd->GetPoint(pointId, startPt);

  vtkNew(vtkPointLocator, pointLocator);
  pointLocator->SetDataSet(boundary);
  pointLocator->BuildLocator();

  int endPtId = pointLocator->FindClosestPoint(closePt);
  returnStartId = boundary->GetPointData()->GetArray(this->InternalIdsArrayName)->
    GetTuple1(endPtId);
  //fprintf(stdout,"Close point Id is: %d\n", returnStartId);

  double minDist = 1.0e10;
  int minId = -1;
  for (int i=0; i<boundary->GetNumberOfPoints(); i++)
  {
    double pt[3];
    boundary->GetPoint(i, pt);
    double dist = vtkSVGeneralUtils::Distance(startPt, pt);
    if (dist < minDist)
    {
      minDist = dist;
      minId = i;
    }
  }
  //fprintf(stdout,"Locator ID: %d\n", endPtId);
  //fprintf(stdout,"MinDist ID: %d\n", minId);


  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVPolyDataSliceAndDiceFilter::AddSurgeryPoints(vtkIdList *surgeryLineIds)
{
  // Add all surgery points in list
  int numToAdd = surgeryLineIds->GetNumberOfIds();

  // Id list for new node numbers, easiest way to add new cell
  vtkNew(vtkIdList, newPointIds);
  newPointIds->SetNumberOfIds(numToAdd);
  for (int i=0; i<numToAdd; i++)
  {
    // Add new point and corresponding data
    double pt[3];
    this->WorkPd->GetPoint(surgeryLineIds->GetId(i), pt);
    int newPtId = this->SurgeryLinesPd->GetPoints()->InsertNextPoint(pt);
    newPointIds->SetId(i, newPtId);
    this->SurgeryLinesPd->GetPointData()->GetArray(this->InternalIdsArrayName)
      ->InsertNextTuple1(surgeryLineIds->GetId(i));
  }

  // Add new points added to cell
  this->SurgeryLinesPd->InsertNextCell(VTK_POLY_LINE, newPointIds);
  return SV_OK;
}
