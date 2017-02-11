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
  this->SurgeryLines = vtkPolyData::New();

  vtkNew(vtkPoints, internalPoints);
  this->SetPoints(internalPoints);

  vtkNew(vtkIntArray, cubeType);
  cubeType->SetName("CubeType");
  this->GetCellData()->AddArray(cubeType);

  vtkNew(vtkIntArray, corners);
  corners->SetNumberOfComponents(8);
  corners->SetName("CornerPtIds");
  this->GetCellData()->AddArray(corners);

  vtkNew(vtkDoubleArray, topNormals);
  topNormals->SetNumberOfComponents(3);
  topNormals->SetName("TopNormal");
  this->GetCellData()->AddArray(topNormals);

  vtkNew(vtkDoubleArray, rightNormals);
  rightNormals->SetNumberOfComponents(3);
  rightNormals->SetName("RightNormal");
  this->GetCellData()->AddArray(rightNormals);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkGeneralizedPolycube::~vtkGeneralizedPolycube()
{
  if (this->SurgeryLines != NULL)
  {
    this->SurgeryLines->Delete();
    this->SurgeryLines = NULL;
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
    this->GetCellData()->GetArray("CubeType")->InsertNextTuple1(-1);
    double blankCorner[8];
    for (int i=0; i<8; i++)
    {
      blankCorner[i] = -1;
    }
    this->GetCellData()->GetArray("CornerPtIds")->InsertNextTuple(blankCorner);
    double topNormal[3]; topNormal[0] = 0.0; topNormal[1] = 0.0; topNormal[2] = 1.0;
    this->GetCellData()->GetArray("TopNormal")->InsertNextTuple(topNormal);
    double rightNormal[3]; rightNormal[0] = 1.0; rightNormal[1] = 0.0; rightNormal[2] = 0.0;
    this->GetCellData()->GetArray("RightNormal")->InsertNextTuple(rightNormal);
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
int vtkGeneralizedPolycube::InsertGridWithCenter(const int cellId, const double center[3], const double dims[3], const int cubetype)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGridWithCenter(cellId, center, dims, cubetype);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGridWithCenter(const int cellId, const double center[3], const double dims[3], const int cubetype)
{
  double origin[3];

  for (int i=0; i<3; i++)
  {
    origin[i] = center[i] - dims[i]/2.0;
  }

  this->SetGridWithOrigin(cellId, origin, dims, cubetype);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::InsertGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int cubetype)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGridWithOrigin(cellId, origin, dims, cubetype);

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int cubetype, const double topNormal[3], const double rightNormal[3], const int corners[8])
{
  this->GetCellData()->GetArray("TopNormal")->SetTuple(cellId, topNormal);
  this->GetCellData()->GetArray("RightNormal")->SetTuple(cellId, rightNormal);
  for (int i=0; i<8; i++)
  {
  this->GetCellData()->GetArray("CornerPtIds")->SetComponent(cellId, i, corners[i]);
  }
  this->SetGridWithOrigin(cellId, origin, dims, cubetype);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int cubetype)
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

  this->SetGrid(cellId, points, cubetype);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::InsertGrid(const int cellId, vtkPoints *points, const int cubetype)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGrid(cellId, points, cubetype);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkGeneralizedPolycube::SetGrid(const int cellId, vtkPoints *points, const int cubetype)
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
  this->GetCellData()->GetArray("CubeType")->SetTuple1(cellId, cubetype);
  this->Modified();

  return 1;
}
