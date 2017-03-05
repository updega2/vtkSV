/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkGeneralizedPolycube.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkGeneralizedPolycube.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkMath.h"
#include "vtkSmartPointer.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkStandardNewMacro(vtkGeneralizedPolycube);

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkGeneralizedPolycube::vtkGeneralizedPolycube()
{
  this->Skeleton = vtkPolyData::New();

  this->Boundaries = vtkIntArray::New();
  this->Corners    = vtkIntArray::New();
  this->Corners->SetNumberOfComponents(8);
  this->TopNormals = vtkDoubleArray::New();
  this->TopNormals->SetNumberOfComponents(3);
  this->RightNormals = vtkDoubleArray::New();
  this->RightNormals->SetNumberOfComponents(3);

  this->InternalPoints = vtkPoints::New();
  this->SetPoints(this->InternalPoints);
  this->Boundaries->SetName("Boundary");
  this->GetCellData()->AddArray(this->Boundaries);
  this->Corners->SetName("CornerPtIds");
  this->GetCellData()->AddArray(this->Corners);
  this->TopNormals->SetName("TopNormal");
  this->GetCellData()->AddArray(this->TopNormals);
  this->RightNormals->SetName("RightNormal");
  this->GetCellData()->AddArray(this->RightNormals);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkGeneralizedPolycube::~vtkGeneralizedPolycube()
{
  if (this->Skeleton != NULL)
  {
    this->Skeleton->Delete();
    this->Skeleton = NULL;
  }
  if (this->Boundaries != NULL)
  {
    this->Boundaries->Delete();
    this->Boundaries = NULL;
  }
  if (this->Corners != NULL)
  {
    this->Corners->Delete();
    this->Corners = NULL;
  }
  if (this->TopNormals != NULL)
  {
    this->TopNormals->Delete();
    this->TopNormals = NULL;
  }
  if (this->RightNormals != NULL)
  {
    this->RightNormals->Delete();
    this->RightNormals = NULL;
  }
  if (this->InternalPoints != NULL)
  {
    this->InternalPoints->Delete();
    this->InternalPoints = NULL;
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkGeneralizedPolycube::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkGeneralizedPolycube::Initialize()
{
  this->Superclass::Initialize();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkGeneralizedPolycube::SetNumberOfGrids(const int numberOfGrids)
{
  int numCurrentGrids = this->GetNumberOfGrids();
  if (numberOfGrids <= numCurrentGrids)
  {
    vtkErrorMacro("Can only increase the number of grids, not decrease");
    return;
  }
  vtkNew(vtkIdList, blankIds);
  blankIds->SetNumberOfIds(8);
  for (int i=0; i<8; i++)
  {
    blankIds->SetId(i, 0);
  }
  for (int i=numCurrentGrids; i <numberOfGrids; i++)
  {
    this->InsertNextCell(VTK_HEXAHEDRON, blankIds);
    this->Boundaries->InsertNextValue(-1);
    double blankCorner[8];
    for (int i=0; i<8; i++)
    {
      blankCorner[i] = -1;
    }
    this->Corners->InsertNextTuple(blankCorner);
    double topNormal[3]; topNormal[0] = 0.0; topNormal[1] = 0.0; topNormal[2] = 1.0;
    this->TopNormals->InsertNextTuple(topNormal);
    double rightNormal[3]; rightNormal[0] = 1.0; rightNormal[1] = 0.0; rightNormal[2] = 0.0;
    this->RightNormals->InsertNextTuple(rightNormal);
  }
  this->Modified();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::GetNumberOfGrids()
{
  return this->GetNumberOfCells();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::InsertGridWithCenter(const int cellId, const double center[3], const double dims[3], const int boundary)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGridWithCenter(cellId, center, dims, boundary);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGridWithCenter(const int cellId, const double center[3], const double dims[3], const int boundary)
{
  double origin[3];

  for (int i=0; i<3; i++)
  {
    origin[i] = center[i] - dims[i]/2.0;
  }

  this->SetGridWithOrigin(cellId, origin, dims, boundary);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::InsertGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int boundary)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGridWithOrigin(cellId, origin, dims, boundary);

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int boundary, const double topNormal[3], const double rightNormal[3], const int corners[8])
{
  this->TopNormals->SetTuple(cellId, topNormal);
  this->RightNormals->SetTuple(cellId, rightNormal);
  for (int i=0; i<8; i++)
  {
  this->Corners->SetComponent(cellId, i, corners[i]);
  }
  this->SetGridWithOrigin(cellId, origin, dims, boundary);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int boundary)
{
  vtkNew(vtkPoints, points);
  points->SetNumberOfPoints(8);

  double pts[8][3];
  double zero[3];
  for (int i=0; i<3; i++)
  {
    zero[i] = 0.0;
  }

  for (int i=0; i<8; i++)
  {
    vtkMath::Add(origin, zero, pts[i]);
  }

  pts[1][0] = origin[0] + dims[0];

  pts[2][0] = origin[0] + dims[0];
  pts[2][1] = origin[1] + dims[1];

  pts[3][1] = origin[1] + dims[1];

  pts[4][2] = origin[2] + dims[2];

  pts[5][0] = origin[0] + dims[0];
  pts[5][2] = origin[2] + dims[2];

  pts[6][0] = origin[0] + dims[0];
  pts[6][1] = origin[1] + dims[1];
  pts[6][2] = origin[2] + dims[2];

  pts[7][1] = origin[1] + dims[1];
  pts[7][2] = origin[2] + dims[2];

  for (int i=0; i<8; i++)
  {
    points->SetPoint(i, pts[i]);
  }

  this->SetGrid(cellId, points, boundary);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::InsertGrid(const int cellId, vtkPoints *points, const int boundary)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGrid(cellId, points, boundary);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGrid(const int cellId, vtkPoints *points, const int boundary)
{
  int numPts = points->GetNumberOfPoints();
  if (numPts != 8)
  {
    vtkErrorMacro("Can only add grid with the 8 corner points");
    return 0;
  }

  vtkIdType pointIds[8];
  vtkIdType npts = 8;
  int numPCPts = this->GetNumberOfPoints();

  for (int i=0; i<numPts; i++)
  {
    this->GetPoints()->InsertNextPoint(points->GetPoint(i));
    pointIds[i] = numPCPts++;
  }
  if (cellId < this->GetNumberOfCells())
  {
    this->ReplaceCell(cellId, npts, pointIds);
  }
  else
  {
    vtkErrorMacro("CellId is greater than number of cells, must set number of grids to correct size first");
    return 0;
  }
  this->Boundaries->SetTuple1(cellId, boundary);
  this->Modified();

  return 1;
}
