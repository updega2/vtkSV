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

/** @file vtkSVMapInterpolator.cxx
 *  @brief This implements the vtkSVMapInterpolator filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVMapInterpolator.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkGradientFilter.h"
#include "vtkLoopSubdivisionFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkTriangle.h"
#include "vtkWarpVector.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
//vtkCxxRevisionMacro(vtkSVMapInterpolator, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSVMapInterpolator);


//---------------------------------------------------------------------------
vtkSVMapInterpolator::vtkSVMapInterpolator()
{
  this->SetNumberOfInputPorts(3);

  this->Verbose               = 1;
  this->NumSourceSubdivisions = 0;
  this->HasBoundary           = 0;

  this->SourceS2Pd = vtkPolyData::New();
  this->TargetPd   = vtkPolyData::New();
  this->TargetS2Pd = vtkPolyData::New();
  this->MappedPd   = vtkPolyData::New();
  this->MappedS2Pd = vtkPolyData::New();

  this->TargetBoundary = vtkIntArray::New();
  this->SourceBoundary = vtkIntArray::New();
}

//---------------------------------------------------------------------------
vtkSVMapInterpolator::~vtkSVMapInterpolator()
{
  if (this->SourceS2Pd != NULL)
  {
    this->SourceS2Pd->Delete();
  }
  if (this->TargetPd != NULL)
  {
    this->TargetPd->Delete();
  }
  if (this->TargetS2Pd != NULL)
  {
    this->TargetS2Pd->Delete();
  }
  if (this->MappedPd != NULL)
  {
    this->MappedPd->Delete();
  }
  if (this->MappedS2Pd != NULL)
  {
    this->MappedS2Pd->Delete();
  }
  if (this->TargetBoundary != NULL)
  {
    this->TargetBoundary->Delete();
  }
  if (this->SourceBoundary != NULL)
  {
    this->SourceBoundary->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkSVMapInterpolator::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkSVMapInterpolator::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *input2 = vtkPolyData::GetData(inputVector[1]);
  vtkPolyData *input3 = vtkPolyData::GetData(inputVector[2]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->SourceS2Pd->DeepCopy(input1);
  this->TargetPd->DeepCopy(input2);
  this->TargetS2Pd->DeepCopy(input3);

  vtkIdType numSourcePolys = this->SourceS2Pd->GetNumberOfPolys();
  //Check the input to make sure it is there
  if (numSourcePolys < 1)
  {
    vtkDebugMacro("No input!");
    return SV_ERROR;
  }
  vtkIdType numTargetPolys = this->TargetPd->GetNumberOfPolys();
  //Check the input to make sure it is there
  if (numTargetPolys < 1)
  {
    vtkDebugMacro("No input!");
    return SV_ERROR;
  }

  if (this->MatchBoundaries() != SV_OK)
  {
    vtkErrorMacro("Error matching the boundaries of the surfaces");
    return SV_ERROR;
  }

  if (this->SubdivideAndInterpolate() != SV_OK)
  {
    return SV_ERROR;
  }

  output->DeepCopy(this->MappedPd);
  output->GetPointData()->PassData(input1->GetPointData());
  output->GetCellData()->PassData(input1->GetCellData());
  if (vtkSVGeneralUtils::CheckArrayExists(output, 0 , "Normals") == 1)
  {
    output->GetPointData()->RemoveArray("Normals");
  }
  if (vtkSVGeneralUtils::CheckArrayExists(output, 1, "cellNormals") == 1)
  {
    output->GetCellData()->RemoveArray("cellNormals");
  }
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVMapInterpolator::SubdivideAndInterpolate()
{
  if (this->NumSourceSubdivisions != 0)
  {
    vtkNew(vtkLoopSubdivisionFilter, subdivider);
    subdivider->SetInputData(this->SourceS2Pd);
    subdivider->SetNumberOfSubdivisions(this->NumSourceSubdivisions);
    subdivider->Update();
    this->MappedS2Pd->DeepCopy(subdivider->GetOutput());
  }
  else
  {
    this->MappedS2Pd->DeepCopy(this->SourceS2Pd);
  }

  if (this->InterpolateMapOntoSource(this->MappedS2Pd,
                                     this->TargetS2Pd,
                                     this->TargetPd,
                                     this->MappedPd) != SV_OK)
  {
    vtkErrorMacro("Error interpolating onto original target surface");
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
int vtkSVMapInterpolator::InterpolateMapOntoSource(vtkPolyData *mappedSourcePd,
                                                 vtkPolyData *mappedTargetPd,
                                                 vtkPolyData *originalTargetPd,
                                                 vtkPolyData *sourceToTargetPd)
{
  vtkNew(vtkCellLocator, locator);

  locator->SetDataSet(mappedTargetPd);
  locator->BuildLocator();

  sourceToTargetPd->DeepCopy(mappedSourcePd);
  int numPts = mappedSourcePd->GetNumberOfPoints();
  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    mappedSourcePd->GetPoint(i, pt);

    double closestPt[3];
    vtkIdType closestCell;
    int subId;
    double distance;
    vtkNew(vtkGenericCell, genericCell);
    locator->FindClosestPoint(pt, closestPt, genericCell, closestCell, subId,
                              distance);

    vtkIdType npts, *pts;
    mappedTargetPd->GetCellPoints(closestCell, npts, pts);
    double pt0[3], pt1[3], pt2[3], a0, a1, a2;
    mappedTargetPd->GetPoint(pts[0], pt0);
    mappedTargetPd->GetPoint(pts[1], pt1);
    mappedTargetPd->GetPoint(pts[2], pt2);
    double area = 0.0;
    vtkSVGeneralUtils::GetBarycentricCoordinates(closestPt, pt0, pt1, pt2, a0, a1, a2);

    double realPt0[3], realPt1[3], realPt2[3];
    originalTargetPd->GetPoint(pts[0], realPt0);
    originalTargetPd->GetPoint(pts[1], realPt1);
    originalTargetPd->GetPoint(pts[2], realPt2);
    double newPoint[3];
    for (int j=0; j<3; j++)
    {
      newPoint[j] = a0*realPt0[j] + a1*realPt1[j] + a2*realPt2[j];
    }
    sourceToTargetPd->GetPoints()->InsertPoint(i, newPoint);
  }

  return SV_OK;
}

int Sign(double testVal)
{
  int asign;
  bool is_negative = testVal < 0;
  if (is_negative)
  {
    asign = -1;
  }
  else
  {
    asign = 1;
  }
  return asign;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVMapInterpolator::MatchBoundaries()
{
  if (this->FindBoundary(this->TargetS2Pd, this->TargetBoundary) != SV_OK)
  {
    return SV_ERROR;
  }
  if (this->FindBoundary(this->SourceS2Pd, this->SourceBoundary) != SV_OK)
  {
    return SV_ERROR;
  }

  if (this->HasBoundary == 1)
  {
    if (this->MoveBoundaryPoints() != SV_OK)
    {
      return SV_ERROR;
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
int vtkSVMapInterpolator::FindBoundary(vtkPolyData *pd, vtkIntArray *isBoundary)
{
  int numPoints = pd->GetNumberOfPoints();
  int numCells = pd->GetNumberOfCells();
  pd->BuildLinks();

  for (int i=0; i<numPoints; i++)
  {
    isBoundary->InsertValue(i, 0);
  }
  for (int i=0; i<numCells; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    //if (npts != 3)
    //{
    //  return SV_ERROR;
    //}
    for (int j=0; j<npts; j++)
    {
      vtkIdType p0, p1;
      p0 = pts[j];
      p1 = pts[(j+1)%npts];

      vtkNew(vtkIdList, edgeNeighbor);
      pd->GetCellEdgeNeighbors(i, p0, p1, edgeNeighbor);

      if (edgeNeighbor->GetNumberOfIds() == 0)
      {
        isBoundary->InsertValue(p0, 1);
        isBoundary->InsertValue(p1, 1);
        this->HasBoundary = 1;
      }
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
int vtkSVMapInterpolator::MoveBoundaryPoints()
{
  int numPoints = this->SourceS2Pd->GetNumberOfPoints();
  vtkNew(vtkCellLocator, locator);

  locator->SetDataSet(this->TargetS2Pd);
  locator->BuildLocator();

  for (int i=0; i<numPoints; i++)
  {
    if (this->SourceBoundary->GetValue(i) == 1)
    {
      double pt[3];
      this->SourceS2Pd->GetPoint(i, pt);
      double closestPt[3];
      vtkIdType closestCell;
      int subId;
      double distance;
      vtkNew(vtkGenericCell, genericCell);
      locator->FindClosestPoint(pt, closestPt, genericCell, closestCell, subId,
				distance);
      double newPt[3];
      if (this->GetPointOnTargetBoundary(i, closestCell, newPt) != SV_OK)
      {
	return SV_ERROR;
      }
      this->SourceS2Pd->GetPoints()->SetPoint(i, newPt);
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
int vtkSVMapInterpolator::GetPointOnTargetBoundary(int srcPtId, int targCellId, double returnPt[])
{
  double srcPt[3];
  this->SourceS2Pd->GetPoint(srcPtId, srcPt);
  vtkNew(vtkIdList, boundaryPts);
  int numBoundaryPts = this->BoundaryPointsOnCell(this->TargetS2Pd, targCellId, boundaryPts, this->TargetBoundary);

  if (numBoundaryPts == 1)
  {
    int ptId = boundaryPts->GetId(0);
    this->TargetS2Pd->GetPoint(ptId, returnPt);
  }
  else if (numBoundaryPts == 2)
  {
    int ptId0 = boundaryPts->GetId(0);
    int ptId1 = boundaryPts->GetId(1);
    double pt0[3], pt1[3];
    this->TargetS2Pd->GetPoint(ptId0, pt0);
    this->TargetS2Pd->GetPoint(ptId1, pt1);
    this->GetProjectedPoint(pt0, pt1, srcPt, returnPt);
  }
  else if (numBoundaryPts == 3)
  {
    int ptId0;
    int ptId1;
    this->GetClosestTwoPoints(this->TargetS2Pd, srcPt, boundaryPts, ptId0, ptId1);
    double pt0[3], pt1[3];
    this->TargetS2Pd->GetPoint(ptId0, pt0);
    this->TargetS2Pd->GetPoint(ptId1, pt1);
    this->GetProjectedPoint(pt0, pt1, srcPt, returnPt);
  }
  else
  {
    vtkErrorMacro("Boundaries do not match well enough");
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
int vtkSVMapInterpolator::BoundaryPointsOnCell(vtkPolyData *pd, int targCellId, vtkIdList *boundaryPts, vtkIntArray *isBoundary)
{
  int numBounds = 0;
  vtkIdType npts, *pts;
  pd->GetCellPoints(targCellId, npts, pts);
  boundaryPts->Reset();
  for (int j=0; j<npts; j++)
  {
    if (isBoundary->GetValue(pts[j]) == 1)
    {
      boundaryPts->InsertNextId(pts[j]);
      numBounds++;
    }
  }
  if (numBounds == 2)
  {
    vtkNew(vtkIdList, edgeNeighbor);
    int p0 = boundaryPts->GetId(0);
    int p1 = boundaryPts->GetId(1);
    pd->GetCellEdgeNeighbors(targCellId, p0, p1, edgeNeighbor);

    if (edgeNeighbor->GetNumberOfIds() != 0)
    {
      int newCell = edgeNeighbor->GetId(0);
      numBounds = this->BoundaryPointsOnCell(pd, newCell, boundaryPts, isBoundary);
    }
  }

  return numBounds;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVMapInterpolator::GetProjectedPoint(double pt0[], double pt1[], double projPt[], double returnPt[])
{
  double vec0[3];
  double vec1[3];
  for (int i=0; i<3; i++)
  {
    vec0[i] = pt1[i] - pt0[i];
    vec1[i] = projPt[i] - pt0[i];
  }
  double proj = vtkMath::Dot(vec0, vec1);

  double lineVec[3], perpVec[3];
  double norm = vtkMath::Dot(vec0, vec0);
  for (int i=0; i<3; i++)
  {
    returnPt[i] = pt0[i] + proj/norm * vec0[i];
  }
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVMapInterpolator::GetClosestTwoPoints(vtkPolyData *pd, double projPt[], vtkIdList *boundaryPts, int &ptId0, int &ptId1)
{
  double dist[3];
  for (int i=0; i<3; i++)
  {
    int ptId = boundaryPts->GetId(i);
    double pt[3];
    pd->GetPoint(ptId, pt);
    dist[i] = std::sqrt(std::pow(projPt[0]-pt[0], 2.0) +
                        std::pow(projPt[1]-pt[1], 2.0) +
                        std::pow(projPt[2]-pt[2], 2.0));

  }

  if (dist[0] > dist[1])
  {
    ptId0 = boundaryPts->GetId(1);
    if (dist[0] > dist[2])
    {
      ptId1 = boundaryPts->GetId(2);
    }
    else
    {
      ptId1 = boundaryPts->GetId(0);
    }
  }
  else
  {
    ptId0 = boundaryPts->GetId(0);
    if (dist[1] > dist[2])
    {
      ptId1 = boundaryPts->GetId(2);
    }
    else
    {
      ptId1 = boundaryPts->GetId(1);
    }
  }
  return SV_OK;
}
