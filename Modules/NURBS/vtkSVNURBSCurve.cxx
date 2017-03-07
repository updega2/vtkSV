/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSVNURBSCurve.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkSVNURBSCurve.h"

#include "vtkSVNURBSUtils.h"
#include "vtkPointData.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkSparseArray.h"
#include "vtkSVGlobals.h"

vtkStandardNewMacro(vtkSVNURBSCurve);

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSCurve::vtkSVNURBSCurve()
{
  this->NumberOfControlPoints = 0;
  this->NumberOfKnotPoints    = 0;
  this->Degree                = 0;
  this->Clamped               = 1;
  this->Closed                = 0;

  this->ControlPointGrid   = vtkSVControlGrid::New();
  this->KnotVector         = vtkDoubleArray::New();
  this->Weights            = vtkDoubleArray::New();

  this->CurveRepresentation = vtkPolyData::New();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSCurve::~vtkSVNURBSCurve()
{
  if (this->ControlPointGrid != NULL)
  {
    this->ControlPointGrid->Delete();
  }
  if (this->KnotVector != NULL)
  {
    this->KnotVector->Delete();
  }
  if (this->Weights != NULL)
  {
    this->Weights->Delete();
  }

  if (this->CurveRepresentation != NULL)
  {
    this->CurveRepresentation->Delete();
  }
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSCurve::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSCurve::Initialize()
{
  this->Superclass::Initialize();
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSCurve* vtkSVNURBSCurve::GetData(vtkInformation* info)
{
  return info? vtkSVNURBSCurve::SafeDownCast(info->Get(DATA_OBJECT())) : 0;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
vtkSVNURBSCurve* vtkSVNURBSCurve::GetData(vtkInformationVector* v, int i)
{
  return vtkSVNURBSCurve::GetData(v->GetInformationObject(i));
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSCurve::SetControlPoints(vtkPoints *points1d)
{
  int nCon = points1d->GetNumberOfPoints();
  this->ControlPointGrid->SetDimensions(nCon, 1, 1);
  this->ControlPointGrid->SetPoints(points1d);
  for (int i=0; i<nCon; i++)
  {
    this->Weights->InsertTuple1(i, 1.0);
  }
  this->Weights->SetName("Weights");
  this->ControlPointGrid->GetPointData()->AddArray(this->Weights);

  this->NumberOfControlPoints = nCon;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
void vtkSVNURBSCurve::SetKnotVector(vtkDoubleArray *knotVector)
{
  int nKnot = knotVector->GetNumberOfTuples();
  this->KnotVector->DeepCopy(knotVector);

  this->NumberOfKnotPoints = nKnot;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVNURBSCurve::GeneratePolyDataRepresentation(const double spacing)
{
  int nCon  = this->NumberOfControlPoints;
  int nKnot = this->NumberOfKnotPoints;
  if (nCon == 0)
  {
    vtkErrorMacro("No control points");
    return 0;
  }
  if (nKnot == 0)
  {
    vtkErrorMacro("No knot points");
    return 0;
  }

  int p = nKnot - nCon - 1;

  int numDiv = ceil(1.0/spacing);
  vtkNew(vtkDoubleArray, uEvals);
  vtkSVNURBSUtils::LinSpace(0, 1, numDiv, uEvals);

  //If nCon - 1 < p, not possible with clamping
  //If nCon - 1 = p, bezier with clamping
  //If nCon - 1 > p, fantastic

  vtkNew(vtkSparseArray<double>, Nus);
  Nus->Resize(numDiv, p+2);
  vtkNew(vtkSparseArray<double>, Nfinal);
  Nfinal->Resize(numDiv, nCon);
  vtkNew(vtkDoubleArray, Nrational);
  Nrational->SetNumberOfTuples(numDiv);
  Nrational->FillComponent(0, 0.0);

  for (int i=0; i<nCon; i++)
  {
    if (vtkSVNURBSUtils::BasisEvaluationVec(this->KnotVector, p,
                                       i, uEvals, Nus) != 1)
    {
      return 0;
    }
    for (int j=0; j<numDiv; j++)
    {
      Nfinal->SetValue(j, i, Nus->GetValue(j, 0));

    }
  }
  Nfinal->SetValue(numDiv-1, nCon-1, 1.0);

  // Multiply by weights for rational curve
  for (int i=0; i<numDiv; i++)
  {
    double ratVal = 0.0;
    for (int j=0; j<nCon; j++)
    {
      double val       = Nfinal->GetValue(i, j);
      double ratWeight = this->Weights->GetTuple1(j);
      ratVal = ratVal + val*ratWeight;
    }
    Nrational->SetTuple1(i, ratVal);
  }
  for (int i=0; i<numDiv; i++)
  {
    for (int j=0; j<nCon; j++)
    {
      double val       = Nfinal->GetValue(i, j);
      double ratVal    = Nrational->GetTuple1(i);
      double ratWeight = this->Weights->GetTuple1(j);
      Nfinal->SetValue(i, j, val*ratWeight/ratVal);
    }
  }

  vtkNew(vtkPoints, surfacePoints);
  if(vtkSVNURBSUtils::MatrixPointsMultiply(Nfinal, this->ControlPointGrid->GetPoints(),
                                       surfacePoints) != 1)
  {
    return 0;
  }

  vtkNew(vtkCellArray, surfaceLines);
  this->GetStructuredGridConnectivity(numDiv, surfaceLines);

  this->CurveRepresentation->SetPoints(surfacePoints);
  this->CurveRepresentation->SetLines(surfaceLines);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkSVNURBSCurve::GetStructuredGridConnectivity(const int numPoints, vtkCellArray *connectivity)
{
  connectivity->Reset();
  vtkNew(vtkIdList, ptIds);
  ptIds->SetNumberOfIds(2);
  for (int i=0; i< numPoints - 1; i++)
  {
    ptIds->SetId(0, i);
    ptIds->SetId(1, i+1);
    connectivity->InsertNextCell(ptIds);
  }

  return 1;
}
