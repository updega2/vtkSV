/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSquareBoundaryMapper.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSquareBoundaryMapper.h"

#include "vtkCellIterator.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkEdgeTable.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <sstream>
#include <map>

//---------------------------------------------------------------------------
vtkStandardNewMacro(vtkSquareBoundaryMapper);
//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSquareBoundaryMapper::SetBoundaries()
{
  //TODO: Clean up correctly on errors!!!
  vtkDataArray *pointIds = this->Boundaries->GetPointData()->GetArray(this->InternalIdsArrayName);

  if (this->CalculateSquareEdgeLengths() != 1)
  {
    vtkErrorMacro("Didn't work");
    return 0;
  }
  fprintf(stdout,"Square edge lengths: %.4f %.4f %.4f %.4f\n", this->BoundaryLengths[0],
                                                               this->BoundaryLengths[1],
                                                               this->BoundaryLengths[2],
                                                               this->BoundaryLengths[3]);

  if (!this->SetSquareBoundary())
  {
    vtkErrorMacro("Was not able to set boundary");
    return 0;
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSquareBoundaryMapper::CalculateSquareEdgeLengths()
{
  vtkDataArray *pointIds = this->BoundaryLoop->GetPointData()->GetArray(this->InternalIdsArrayName);
  int numLines = this->BoundaryLoop->GetNumberOfLines();

  for (int i=0; i<4; i++)
  {
    this->BoundaryLengths[i] = 0.0;
  }
  int currCell = 0;
  for (int i=0; i<4; i++)
  {
    int lastPt   = pointIds->LookupValue(this->BoundaryIds->GetValue((i+1)%4));
    int otherPt0 = pointIds->LookupValue(this->BoundaryIds->GetValue((i+2)%4));
    int otherPt1 = pointIds->LookupValue(this->BoundaryIds->GetValue((i+3)%4));
    int otherPt2 = pointIds->LookupValue(this->BoundaryIds->GetValue(i));
    vtkIdType npts, *pts;
    int checkPt = -1;
    while (checkPt != lastPt)
    {
      this->BoundaryLoop->GetCellPoints(currCell, npts, pts);
      double pt0[3], pt1[3];
      this->BoundaryLoop->GetPoint(pts[0], pt0);
      this->BoundaryLoop->GetPoint(pts[1], pt1);
      checkPt = pts[1];
      if (checkPt == (lastPt))
        fprintf(stdout,"Found corner %d\n", this->BoundaryIds->GetValue((i+1)%4));
      if (checkPt == (otherPt0))
        fprintf(stdout,"Found corner %d\n", this->BoundaryIds->GetValue((i+2)%4));
      if (checkPt == (otherPt1))
        fprintf(stdout,"Found corner %d\n", this->BoundaryIds->GetValue((i+3)%4));
      if (checkPt == (otherPt2))
        fprintf(stdout,"Found corner %d\n", this->BoundaryIds->GetValue(i));

      double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
                              std::pow(pt0[1]-pt1[1], 2.0) +
                              std::pow(pt0[2]-pt1[2], 2.0));
      this->BoundaryLengths[i] += dist;

      currCell++;
    }
  }
  //fprintf(stdout,"Curr!: %d\n", currCell);
  //fprintf(stdout,"NumLines!: %d\n", numLines);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSquareBoundaryMapper::SetSquareBoundary()
{
  vtkDataArray *pointIds = this->BoundaryLoop->GetPointData()->GetArray(this->InternalIdsArrayName);
  double currCoords[3];
  for (int i=0; i<3; i++)
  {
    currCoords[i] = 0.0;
  }

  vtkNew(vtkPoints, newPoints);
  vtkNew(vtkIntArray, newDataArray);

  double unitLength = 1.0;
  int currCell = 0;
  int checkPt = -1;
  for (int i=0; i<4; i++)
  {
    double currLength = 0.0;
    int lastPt  = pointIds->LookupValue(this->BoundaryIds->GetValue((i+1)%4));
    vtkIdType npts, *pts;
    while (checkPt != lastPt)
    {
      this->BoundaryLoop->GetCellPoints(currCell, npts, pts);
      double pt0[3], pt1[3];
      this->BoundaryLoop->GetPoint(pts[0], pt0);
      this->BoundaryLoop->GetPoint(pts[1], pt1);
      checkPt = pts[1];

      double dist = std::sqrt(std::pow(pt0[0]-pt1[0], 2.0) +
                              std::pow(pt0[1]-pt1[1], 2.0) +
                              std::pow(pt0[2]-pt1[2], 2.0));
      currLength += dist;

      if (i == 0)
      {
        currCoords[0] += dist/this->BoundaryLengths[i] * unitLength;
      }
      else if (i == 1)
      {
        currCoords[1] += dist/this->BoundaryLengths[i] * unitLength;
      }
      else if (i == 2)
      {
        currCoords[0] -= dist/this->BoundaryLengths[i] * unitLength;
      }
      else
      {
        currCoords[1] -= dist/this->BoundaryLengths[i] * unitLength;
      }
      newPoints->InsertNextPoint(currCoords);
      newDataArray->InsertNextValue(pointIds->GetTuple1(pts[1]));

      currCell++;
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

  return 1;
}

void vtkSquareBoundaryMapper::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  //os << indent << "Number of subdivisions: "
  //   << this->GetNumberOfSubdivisions() << endl;
}
