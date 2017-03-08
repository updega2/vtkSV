/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVGeneralizedPolycube.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSVGeneralizedPolycube.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGlobals.h"

vtkStandardNewMacro(vtkSVGeneralizedPolycube);

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVGeneralizedPolycube::vtkSVGeneralizedPolycube()
{
  this->SurgeryLines = vtkPolyData::New();

  vtkNew(vtkPoints, internalPoints);
  this->SetPoints(internalPoints);

  vtkNew(vtkIntArray, cubeType);
  cubeType->SetName("CubeType");
  this->GetCellData()->AddArray(cubeType);

  vtkNew(vtkIntArray, parentDirection);
  parentDirection->SetName("ParentDirection");
  this->GetCellData()->AddArray(parentDirection);
  vtkNew(vtkIntArray, childDirection);
  childDirection->SetName("ChildDirection");
  this->GetCellData()->AddArray(childDirection);

  vtkNew(vtkIntArray, corners);
  corners->SetNumberOfComponents(8);
  corners->SetName("CornerPtIds");
  this->GetCellData()->AddArray(corners);

  vtkNew(vtkIntArray, pointCorners);
  pointCorners->SetNumberOfComponents(1);
  pointCorners->SetName("CornerPtIds");
  this->GetPointData()->AddArray(pointCorners);

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
vtkSVGeneralizedPolycube::~vtkSVGeneralizedPolycube()
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
void vtkSVGeneralizedPolycube::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVGeneralizedPolycube::Initialize()
{
  this->Superclass::Initialize();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVGeneralizedPolycube::SetNumberOfGrids(const int numberOfGrids)
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
    this->GetCellData()->GetArray("ParentDirection")->InsertNextTuple1(-1);
    this->GetCellData()->GetArray("ChildDirection")->InsertNextTuple1(-1);
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
int vtkSVGeneralizedPolycube::GetNumberOfGrids()
{
  return this->GetNumberOfCells();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralizedPolycube::InsertGridWithCenter(const int cellId, const double center[3], const double dims[3], const int cubetype, const int parentdirection, const int childdirection)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGridWithCenter(cellId, center, dims, cubetype, parentdirection, childdirection);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralizedPolycube::SetGridWithCenter(const int cellId, const double center[3], const double dims[3], const int cubetype, const int parentdirection, const int childdirection)
{
  double origin[3];

  for (int i=0; i<3; i++)
  {
    origin[i] = center[i] - dims[i]/2.0;
  }

  this->SetGridWithOrigin(cellId, origin, dims, cubetype, parentdirection, childdirection);
  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralizedPolycube::InsertGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int cubetype, const int parentdirection, const int childdirection)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGridWithOrigin(cellId, origin, dims, cubetype, parentdirection, childdirection);

  return SV_OK;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralizedPolycube::SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int cubetype, const int parentdirection, const int childdirection, const double topNormal[3], const double rightNormal[3], const int corners[8])
{
  this->GetCellData()->GetArray("TopNormal")->SetTuple(cellId, topNormal);
  this->GetCellData()->GetArray("RightNormal")->SetTuple(cellId, rightNormal);
  for (int i=0; i<8; i++)
  {
    this->GetCellData()->GetArray("CornerPtIds")->SetComponent(cellId, i, corners[i]);
  }
  this->SetGridWithOrigin(cellId, origin, dims, cubetype, parentdirection, childdirection);

  return SV_OK;
}

/**
 *               2---------1
 *              /|        /|
 *             / |       / |
 *            3---------0  |
 *            |  |      |  |
 *            |  6------|--5
 *            | /       | /
 *            |/        |/
 *            7---------4
 */
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralizedPolycube::SetGridWithOrigin(const int cellId, const double origin[3], const double dims[3], const int cubetype, const int parentdirection, const int childdirection)
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

  pts[4][0] = origin[0] + dims[0];

  pts[5][0] = origin[0] + dims[0];
  pts[5][1] = origin[1] + dims[1];

  pts[6][1] = origin[1] + dims[1];

  pts[3][2] = origin[2] + dims[2];

  pts[0][0] = origin[0] + dims[0];
  pts[0][2] = origin[2] + dims[2];

  pts[1][0] = origin[0] + dims[0];
  pts[1][1] = origin[1] + dims[1];
  pts[1][2] = origin[2] + dims[2];

  pts[2][1] = origin[1] + dims[1];
  pts[2][2] = origin[2] + dims[2];

  for (int i=0; i<8; i++)
  {
    points->SetPoint(i, pts[i]);
  }

  this->SetGrid(cellId, points, cubetype, parentdirection, childdirection);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralizedPolycube::InsertGrid(const int cellId, vtkPoints *points, const int cubetype, const int parentdirection, const int childdirection)
{
  if (cellId >= this->GetNumberOfCells())
  {
    this->SetNumberOfGrids(cellId);
  }
  this->SetGrid(cellId, points, cubetype, parentdirection, childdirection);

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVGeneralizedPolycube::SetGrid(const int cellId, vtkPoints *points, const int cubetype, const int parentdirection, const int childdirection)
{
  int numPts = points->GetNumberOfPoints();
  if (numPts != 8)
  {
    vtkErrorMacro("Can only add grid with the 8 corner points");
    return SV_ERROR;
  }

  vtkIdType pointIds[8];
  vtkIdType npts = 8;
  int numPCPts = this->GetNumberOfPoints();

  for (int i=0; i<numPts; i++)
  {
    this->GetPoints()->InsertNextPoint(points->GetPoint(i));
    this->GetPointData()->GetArray("CornerPtIds")->InsertNextTuple1(-1.0);
    pointIds[i] = numPCPts++;
  }
  if (cellId < this->GetNumberOfCells())
  {
    this->ReplaceCell(cellId, npts, pointIds);
  }
  else
  {
    vtkErrorMacro("CellId is greater than number of cells, must set number of grids to correct size first");
    return SV_ERROR;
  }
  this->GetCellData()->GetArray("CubeType")->SetTuple1(cellId, cubetype);
  this->GetCellData()->GetArray("ParentDirection")->SetTuple1(cellId, parentdirection);
  this->GetCellData()->GetArray("ChildDirection")->SetTuple1(cellId, childdirection);
  this->Modified();

  return SV_OK;
}
