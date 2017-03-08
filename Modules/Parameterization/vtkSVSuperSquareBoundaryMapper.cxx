/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVSuperSquareBoundaryMapper.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSVSuperSquareBoundaryMapper.h"

#include "vtkCellIterator.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkEdgeTable.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"

#include <sstream>
#include <map>

//---------------------------------------------------------------------------
vtkStandardNewMacro(vtkSVSuperSquareBoundaryMapper);
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVSuperSquareBoundaryMapper::vtkSVSuperSquareBoundaryMapper()
{
  this->BoundaryLengths = vtkDoubleArray::New();
  this->SetSuperBoundaryDivisions(0, 0, 0, 0); //regular square boundary
  this->SetSuperBoundaryLengths(1.0, 1.0, 1.0, 1.0); //regular square boundary
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVSuperSquareBoundaryMapper::~vtkSVSuperSquareBoundaryMapper()
{
  if (this->BoundaryLengths != NULL)
  {
    this->BoundaryLengths->Delete();
  }
}
//
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVSuperSquareBoundaryMapper::SetBoundaries()
{
  vtkNew(vtkIntArray, actualIds);
  if (this->CalculateSquareEdgeLengths(actualIds) != SV_OK)
  {
    vtkErrorMacro("Didn't work");
    return SV_ERROR;
  }

  if (!this->SetSquareBoundary(actualIds))
  {
    vtkErrorMacro("Was not able to set boundary");
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
int vtkSVSuperSquareBoundaryMapper::CalculateSquareEdgeLengths(vtkIntArray *actualIds)
{
  vtkDataArray *pointIds = this->BoundaryLoop->GetPointData()->GetArray(this->InternalIdsArrayName);
  int numLines = this->BoundaryLoop->GetNumberOfLines();

  int numBoundaryPts = this->BoundaryIds->GetNumberOfTuples();
  this->BoundaryLengths->SetNumberOfComponents(1);
  this->BoundaryLengths->SetNumberOfTuples(numBoundaryPts);
  actualIds->SetNumberOfComponents(1);
  actualIds->SetNumberOfTuples(numBoundaryPts);
  int currCell = 0;
  for (int i=0; i<numBoundaryPts; i++)
  {
    //int lastPt = pointIds->LookupValue(this->BoundaryIds->GetValue((i+1)%numBoundaryPts));
    vtkIdType npts, *pts;
    int checkPt = -1;
    double boundaryDistance = 0.0;
    int done = 0;
    //while (checkPt != lastPt)
    while (!done)
    {
      this->BoundaryLoop->GetCellPoints(currCell, npts, pts);
      double pt0[3], pt1[3];
      this->BoundaryLoop->GetPoint(pts[0], pt0);
      this->BoundaryLoop->GetPoint(pts[1], pt1);
      checkPt = pts[1];
      for (int j=0; j<numBoundaryPts; j++)
      {
        if (checkPt == pointIds->LookupValue(this->BoundaryIds->GetValue(j)))
        {
          actualIds->SetValue((i+1)%numBoundaryPts, pointIds->GetTuple1(checkPt));
          fprintf(stdout,"Found boundary ID!: %d\n", this->BoundaryIds->GetValue(j));
          done = 1;
        }
      }

      double dist = sqrt(pow(pt0[0]-pt1[0], 2.0) +
                         pow(pt0[1]-pt1[1], 2.0) +
                         pow(pt0[2]-pt1[2], 2.0));
      boundaryDistance += dist;

      currCell++;
      if (currCell > this->BoundaryLoop->GetNumberOfCells())
      {
        vtkErrorMacro("Error could not find all boundary points provided\n");
        return SV_ERROR;
      }
    }
    this->BoundaryLengths->SetTuple1(i, boundaryDistance);
  }
  //fprintf(stdout,"Curr!: %d\n", currCell);
  //fprintf(stdout,"NumLines!: %d\n", numLines);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVSuperSquareBoundaryMapper::SetSquareBoundary(vtkIntArray *actualIds)
{
  vtkDataArray *pointIds = this->BoundaryLoop->GetPointData()->GetArray(this->InternalIdsArrayName);
  double currCoords[3];
  for (int i=0; i<3; i++)
  {
    currCoords[i] = 0.0;
  }

  int numBoundaryPts = actualIds->GetNumberOfTuples();

  vtkNew(vtkPoints, newPoints);
  vtkNew(vtkIntArray, newDataArray);

  int currCell = 0;
  int checkPt = -1;
  int boundaryNumber = 0;
  int divisionCount = 0;
  for (int i=0; i<numBoundaryPts; i++)
  {
    fprintf(stdout,"Looping to next point: %d\n", i);
    double currLength = 0.0;
    int lastPt  = pointIds->LookupValue(actualIds->GetValue((i+1)%numBoundaryPts));
    vtkIdType npts, *pts;
    while (checkPt != lastPt)
    {
      this->BoundaryLoop->GetCellPoints(currCell, npts, pts);
      double pt0[3], pt1[3];
      this->BoundaryLoop->GetPoint(pts[0], pt0);
      this->BoundaryLoop->GetPoint(pts[1], pt1);
      checkPt = pts[1];

      double dist = sqrt((pt0[0]-pt1[0], 2.0) +
                         (pt0[1]-pt1[1], 2.0) +
                         (pt0[2]-pt1[2], 2.0));
      currLength += dist;

      double unitLength = this->SuperBoundaryLengths[boundaryNumber]/(this->SuperBoundaryDivisions[boundaryNumber]+1.0);
      if (boundaryNumber == 0)
      {
        currCoords[0] += dist/this->BoundaryLengths->GetTuple1(i) * unitLength;
      }
      else if (boundaryNumber == 1)
      {
        currCoords[1] += dist/this->BoundaryLengths->GetTuple1(i) * unitLength;
      }
      else if (boundaryNumber == 2)
      {
        currCoords[0] -= dist/this->BoundaryLengths->GetTuple1(i) * unitLength;
      }
      else
      {
        currCoords[1] -= dist/this->BoundaryLengths->GetTuple1(i) * unitLength;
      }
      newPoints->InsertNextPoint(currCoords);
      newDataArray->InsertNextValue(pointIds->GetTuple1(pts[1]));

      currCell++;
    }
    divisionCount++;
    if (divisionCount > this->SuperBoundaryDivisions[boundaryNumber])
    {
      boundaryNumber++;
      divisionCount = 0;
    }
  }
  vtkNew(vtkCellArray, newCells);
  int i=0;
  for (i=0; i<newPoints->GetNumberOfPoints()-1; i++)
  {
    newCells->InsertNextCell(2);
    newCells->InsertCellPoint(i);
    newCells->InsertCellPoint(i+1);
  }
  newCells->InsertNextCell(2);
  newCells->InsertCellPoint(i-1);
  newCells->InsertCellPoint(0);

  this->BoundaryPd->SetPoints(newPoints);
  this->BoundaryPd->SetLines(newCells);
  newDataArray->SetName(this->InternalIdsArrayName);
  this->BoundaryPd->GetPointData()->AddArray(newDataArray);

  return SV_OK;
}

void vtkSVSuperSquareBoundaryMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  //os << indent << "Number of subdivisions: "
  //   << this->GetNumberOfSubdivisions() << endl;
}