/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkControlGrid.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkControlGrid.h"

#include "vtkDataArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"

vtkStandardNewMacro(vtkControlGrid);

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkControlGrid::vtkControlGrid()
{
  this->InternalPoints = vtkPoints::New();
  this->SetPoints(InternalPoints);

  this->InternalWeights = vtkDoubleArray::New();
  this->InternalWeights->SetName("Weights");
  this->GetPointData()->AddArray(this->InternalWeights);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkControlGrid::~vtkControlGrid()
{
  if (this->InternalPoints != NULL)
  {
    this->InternalPoints->Delete();
  }
  if (this->InternalWeights != NULL)
  {
    this->InternalWeights->Delete();
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkControlGrid::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkControlGrid::CopyStructure(vtkDataSet *ds)
{
  this->Superclass::CopyStructure(ds);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkControlGrid::Initialize()
{
  this->Superclass::Initialize();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkControlGrid* vtkControlGrid::GetData(vtkInformation* info)
{
  return info? vtkControlGrid::SafeDownCast(info->Get(DATA_OBJECT())) : 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkControlGrid* vtkControlGrid::GetData(vtkInformationVector* v, int i)
{
  return vtkControlGrid::GetData(v->GetInformationObject(i));
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkControlGrid::SetControlPoint(const int i, const int j, const int k, const double p[3], const double w)
{
  int ptId;
  this->GetPointId(i, j, k, ptId);
  this->GetPoints()->SetPoint(ptId, p);
  this->GetPointData()->GetArray("Weights")->InsertTuple1(ptId, w);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkControlGrid::SetControlPoint(const int i, const int j, const int k, const double pw[4])
{
  double onlyp[3];
  for (int l=0; l<3; l++)
  {
    onlyp[l] = pw[l];
  }
  double w = pw[3];

  this->SetControlPoint(i, j, k, onlyp, w);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkControlGrid::InsertControlPoint(const int i, const int j, const int k, const double p[3], const double w)
{
  int dim[3];
  this->GetDimensions(dim);
  if (i >= dim[0])
  {
    dim[0] = i + 1;
  }
  if (j >= dim[1])
  {
    dim[1] = j + 1;
  }
  if (k >= dim[2])
  {
    dim[2] = k + 1;
  }
  this->SetDimensions(dim);
  this->GetPoints()->SetNumberOfPoints((dim[0]*dim[1]*dim[2]));

  this->SetControlPoint(i, j, k, p, w);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkControlGrid::InsertControlPoint(const int i, const int j, const int k, const double pw[4])
{
  double onlyp[3];
  for (int l=0; l<3; l++)
  {
    onlyp[l] = pw[l];
  }
  double w = pw[3];

  this->InsertControlPoint(i, j, k, onlyp, w);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkControlGrid::GetControlPoint(const int i, const int j, const int k, double p[3], double &w)
{
  int ptId;
  if (this->GetPointId(i, j, k, ptId) != 1)
  {
    vtkErrorMacro("Point not retrieved successfully");
    return 0;
  }
  this->GetPoint(ptId, p);
  w = this->GetPointData()->GetArray("Weights")->GetTuple1(ptId);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkControlGrid::GetControlPoint(const int i, const int j, const int k, double pw[4])
{
  double onlyp[3];
  for (int l=0; l<3; l++)
  {
    onlyp[l] = pw[l];
  }
  double w = pw[3];

  if (this->GetControlPoint(i, j, k, onlyp, w) != 1)
  {
    vtkErrorMacro("Point not retrieved successfully");
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkControlGrid::GetPointId(const int i, const int j, const int k, int &ptId)
{
  int extent[6];
  this->GetExtent(extent);

  fprintf(stdout,"Extents: %d %d %d %d %d %d\n", extent[0],
                                                    extent[1],
                                                    extent[2],
                                                    extent[3],
                                                    extent[4],
                                                    extent[5]);
  if(i < extent[0] || i > extent[1] ||
     j < extent[2] || j > extent[3] ||
     k < extent[4] || k > extent[5])
    {
    vtkErrorMacro("ERROR: IJK coordinates are outside of grid extent!");
    return 0; // out of bounds!
    }

  int pos[3];
  pos[0] = i;
  pos[1] = j;
  pos[2] = k;

  ptId = vtkStructuredData::ComputePointIdForExtent(extent, pos);

  return 1;
}
