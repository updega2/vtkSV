/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVNURBSSurface.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSVNURBSSurface.h"

#include "vtkCellArray.h"
#include "vtkSVNURBSUtils.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSparseArray.h"
#include "vtkSVGlobals.h"

vtkStandardNewMacro(vtkSVNURBSSurface);

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSSurface::vtkSVNURBSSurface()
{
  this->NumberOfUControlPoints = 0;
  this->NumberOfVControlPoints = 0;
  this->NumberOfUKnotPoints    = 0;
  this->NumberOfVKnotPoints    = 0;
  this->UDegree                = 0;
  this->VDegree                = 0;
  this->UClamped               = 1;
  this->UClosed                = 0;
  this->VClosed                = 0;

  this->ControlPointGrid    = vtkSVControlGrid::New();

  for (int i=0; i<2; i++)
  {
    this->UVKnotVectors[i] = vtkDoubleArray::New();
    this->UVWeights[i]     = vtkDoubleArray::New();
  }
  this->UKnotVector = this->UVKnotVectors[0];
  this->VKnotVector = this->UVKnotVectors[1];
  this->UWeights = this->UVWeights[0];
  this->VWeights = this->UVWeights[1];

  this->SurfaceRepresentation = vtkPolyData::New();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSSurface::~vtkSVNURBSSurface()
{
  if (this->ControlPointGrid != NULL)
  {
    this->ControlPointGrid->Delete();
  }
  for (int i=0; i<2; i++)
  {
    if (this->UVKnotVectors[i] != NULL)
    {
      this->UVKnotVectors[i]->Delete();
    }
    if (this->UVWeights[i] != NULL)
    {
      this->UVWeights[i]->Delete();
    }
  }

  if (this->SurfaceRepresentation != NULL)
  {
    this->SurfaceRepresentation->Delete();
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSSurface::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSSurface::Initialize()
{
  this->Superclass::Initialize();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSSurface* vtkSVNURBSSurface::GetData(vtkInformation* info)
{
  return info? vtkSVNURBSSurface::SafeDownCast(info->Get(DATA_OBJECT())) : 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSSurface* vtkSVNURBSSurface::GetData(vtkInformationVector* v, int i)
{
  return vtkSVNURBSSurface::GetData(v->GetInformationObject(i));
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSSurface::SetControlPoints(vtkStructuredGrid *points2d)
{
  int dim[3];
  points2d->GetDimensions(dim);
  this->ControlPointGrid->SetPoints(points2d->GetPoints());
  this->ControlPointGrid->SetDimensions(dim);
  this->ControlPointGrid->GetPointData()->GetArray("Weights")
    ->SetNumberOfTuples(dim[0]*dim[1]);
  this->ControlPointGrid->GetPointData()->GetArray("Weights")
    ->FillComponent(0, 1.0);
  this->UWeights->SetNumberOfTuples(dim[0]);
  this->UWeights->FillComponent(0, 1.0);
  this->VWeights->SetNumberOfTuples(dim[1]);
  this->VWeights->FillComponent(0, 1.0);

  this->NumberOfUControlPoints = dim[0];
  this->NumberOfVControlPoints = dim[1];
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVNURBSSurface::SetControlPoints(vtkIntArray *indices, const int dim,
                                       vtkPoints *coordinates, vtkDoubleArray *weights)
{
  int numInserts = indices->GetNumberOfTuples();

  for (int i=0; i<numInserts; i++)
  {
    int index = indices->GetTuple1(i);
    double pt[3];
    coordinates->GetPoint(i, pt);
    double weight = weights->GetTuple1(i);

    this->SetControlPoint(index, dim, pt, weight);
  }

  //SetControlPoint is still returning 0
  return SV_ERROR;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSSurface::SetKnotVector(vtkDoubleArray *knotVector, const int dim)
{
  int nKnot = knotVector->GetNumberOfTuples();
  this->UVKnotVectors[dim]->DeepCopy(knotVector);

  if (dim == 0)
  {
    this->NumberOfUKnotPoints = nKnot;
  }
  else
  {
    this->NumberOfVKnotPoints = nKnot;
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVNURBSSurface::GeneratePolyDataRepresentation(const double uSpacing,
		                                    const double vSpacing)
{
  int nUCon  = this->NumberOfUControlPoints;
  int nVCon  = this->NumberOfVControlPoints;
  int nUKnot = this->NumberOfUKnotPoints;
  int nVKnot = this->NumberOfVKnotPoints;
  if (nUCon == 0 || nVCon == 0)
  {
    vtkErrorMacro("No control points");
    return SV_ERROR;
  }
  if (nUKnot == 0 || nVKnot == 0)
  {
    vtkErrorMacro("No knot points");
    return SV_ERROR;
  }

  int p = nUKnot - nUCon - 1;
  int q = nVKnot - nVCon - 1;

  //If nCon - 1 < p, not possible with clamping
  //If nCon - 1 = p, bezier with clamping
  //If nCon - 1 > p, fantastic

  // U direction!
  // -----------------------------------------------------------------------
  int numUDiv = ceil(1.0/uSpacing);
  vtkNew(vtkDoubleArray, uEvals);
  vtkSVNURBSUtils::LinSpace(0, 1, numUDiv, uEvals);

  vtkNew(vtkSparseArray<double>, Nus);
  Nus->Resize(numUDiv, p+2);
  vtkNew(vtkSparseArray<double>, NUfinal);
  NUfinal->Resize(numUDiv, nUCon);
  vtkNew(vtkDoubleArray, NUrational);
  NUrational->SetNumberOfTuples(numUDiv);
  NUrational->FillComponent(0, 0.0);

  for (int i=0; i<nUCon; i++)
  {
    if (vtkSVNURBSUtils::BasisEvaluationVec(this->UKnotVector, p,
                                       i, uEvals, Nus) != SV_OK)
    {
      return SV_ERROR;
    }
    for (int j=0; j<numUDiv; j++)
    {
      NUfinal->SetValue(j, i, Nus->GetValue(j, 0));
    }
  }
  NUfinal->SetValue(numUDiv-1, nUCon-1, 1.0);

  // Multiply by weights for rational curve
  for (int i=0; i<numUDiv; i++)
  {
    double ratVal = 0.0;
    for (int j=0; j<nUCon; j++)
    {
      double val = NUfinal->GetValue(i, j);
      double ratWeight = this->UWeights->GetTuple1(j);
      ratVal = ratVal + val*ratWeight;
    }
    NUrational->SetTuple1(i, ratVal);
  }
  for (int i=0; i<numUDiv; i++)
  {
    for (int j=0; j<nUCon; j++)
    {
      double val       = NUfinal->GetValue(i, j);
      double ratVal    = NUrational->GetTuple1(i);
      double ratWeight = this->UWeights->GetTuple1(j);
      NUfinal->SetValue(i, j, val*ratWeight/ratVal);
    }
  }

  // V direction!
  // -----------------------------------------------------------------------
  int numVDiv = ceil(1.0/vSpacing);
  vtkNew(vtkDoubleArray, vEvals);
  vtkSVNURBSUtils::LinSpace(0, 1, numVDiv, vEvals);

  vtkNew(vtkSparseArray<double>, Nvs);
  Nvs->Resize(numVDiv, q+2);
  vtkNew(vtkSparseArray<double>, NVfinal);
  NVfinal->Resize(numVDiv, nVCon);
  vtkNew(vtkDoubleArray, NVrational);
  NVrational->SetNumberOfTuples(numVDiv);
  NVrational->FillComponent(0, 0.0);

  for (int i=0; i<nVCon; i++)
  {
    if (vtkSVNURBSUtils::BasisEvaluationVec(this->VKnotVector, q,
                                       i, vEvals, Nvs) != SV_OK)
    {
      return SV_ERROR;
    }
    double ratVal = 0.0;
    for (int j=0; j<numVDiv; j++)
    {
      NVfinal->SetValue(j, i, Nvs->GetValue(j, 0));
    }
  }
  NVfinal->SetValue(numVDiv-1, nVCon-1, 1.0);

  // Multiply by weights for rational curve
  for (int i=0; i<numVDiv; i++)
  {
    double ratVal = 0.0;
    for (int j=0; j<nVCon; j++)
    {
      double val = NVfinal->GetValue(i, j);
      double ratWeight = this->VWeights->GetTuple1(j);
      ratVal = ratVal + val*ratWeight;
    }
    NVrational->SetTuple1(i, ratVal);
  }
  for (int i=0; i<numVDiv; i++)
  {
    for (int j=0; j<nVCon; j++)
    {
      double val       = NVfinal->GetValue(i, j);
      double ratVal    = NVrational->GetTuple1(i);
      double ratWeight = this->VWeights->GetTuple1(j);
      NVfinal->SetValue(i, j, val*ratWeight/ratVal);
    }
  }

  vtkNew(vtkSparseArray<double>, NVfinalT);
  vtkSVNURBSUtils::MatrixTranspose(NVfinal, 0, NVfinalT);
  //Get the physical points on the surface!
  // -----------------------------------------------------------------------
  vtkNew(vtkDenseArray<double>, tmpControlGrid);
  vtkSVNURBSUtils::StructuredGridToTypedArray(this->ControlPointGrid, tmpControlGrid);

  vtkNew(vtkDenseArray<double>, tmpUGrid);
  if (vtkSVNURBSUtils::MatrixMatrixMultiply(NUfinal, 0, tmpControlGrid, 1, tmpUGrid) != SV_OK)
  {
    fprintf(stderr, "Error in matrix multiply\n");
    return SV_ERROR;
  }
  vtkNew(vtkDenseArray<double>, tmpVGrid);
  if (vtkSVNURBSUtils::MatrixMatrixMultiply(tmpUGrid, 1, NVfinalT, 0, tmpVGrid) != SV_OK)
  {
    fprintf(stderr, "Error in matrix multiply\n");
    return SV_ERROR;
  }

  vtkNew(vtkStructuredGrid, finalGrid);
  vtkNew(vtkPoints, tmpVPoints);
  finalGrid->SetPoints(tmpVPoints);
  vtkSVNURBSUtils::TypedArrayToStructuredGrid(tmpVGrid, finalGrid);

  vtkNew(vtkCellArray, surfaceCells);
  this->GetStructuredGridConnectivity(numUDiv, numVDiv, surfaceCells);

  this->SurfaceRepresentation->SetPoints(finalGrid->GetPoints());
  this->SurfaceRepresentation->SetPolys(surfaceCells);
  this->SurfaceRepresentation->BuildLinks();

  return SV_OK;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVNURBSSurface::GetStructuredGridConnectivity(const int numXPoints, const int numYPoints, vtkCellArray *connectivity)
{
  connectivity->Reset();
  vtkNew(vtkIdList, ptIds);
  ptIds->SetNumberOfIds(4);
  for (int i=0; i< numXPoints - 1; i++)
  {
    for (int j=0; j< numYPoints - 1; j++)
    {
      int ptId = i+ j*numXPoints;
      ptIds->SetId(0, ptId);
      ptIds->SetId(1, ptId+1);
      ptIds->SetId(2, ptId+numXPoints+1);
      ptIds->SetId(3, ptId+numXPoints);
      connectivity->InsertNextCell(ptIds);
    }
  }

  return SV_OK;
}
