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

#include "vtkSVNURBSSurface.h"

#include "vtkCellArray.h"
#include "vtkCleanPolyData.h"
#include "vtkSVNURBSUtils.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSparseArray.h"
#include "vtkSVGlobals.h"

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVNURBSSurface);

// ----------------------
// Constructor
// ----------------------
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
    this->UVKnotVectors[i] = vtkDoubleArray::New();

  this->UKnotVector = this->UVKnotVectors[0];
  this->VKnotVector = this->UVKnotVectors[1];

  this->SurfaceRepresentation = vtkPolyData::New();
}

// ----------------------
// Destructor
// ----------------------
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
  }

  if (this->SurfaceRepresentation != NULL)
  {
    this->SurfaceRepresentation->Delete();
  }
}

// ----------------------
// DeepCopy
// ----------------------
void vtkSVNURBSSurface::DeepCopy(vtkSVNURBSSurface *src)
{
  this->Superclass::DeepCopy(src);

  this->SetNumberOfUControlPoints(src->GetNumberOfUControlPoints());
  this->SetNumberOfVControlPoints(src->GetNumberOfVControlPoints());
  this->SetNumberOfUKnotPoints(src->GetNumberOfUKnotPoints());
  this->SetNumberOfVKnotPoints(src->GetNumberOfVKnotPoints());
  this->SetUDegree(src->GetUDegree());
  this->SetVDegree(src->GetVDegree());

  this->ControlPointGrid->DeepCopy(src->GetControlPointGrid());
  this->UKnotVector->DeepCopy(src->GetUKnotVector());
  this->VKnotVector->DeepCopy(src->GetVKnotVector());

  this->SurfaceRepresentation->DeepCopy(src->GetSurfaceRepresentation());
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVNURBSSurface::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Number of control points in u direction: " << this->NumberOfUControlPoints << "\n";
  os << indent << "Number of knot points in u direction: " << this->NumberOfUKnotPoints << "\n";
  os << indent << "U Degree: " << this->UDegree << "\n";
  os << indent << "U Clamped: " << this->UClamped << "\n";
  os << indent << "U Closed: " << this->UClosed << "\n";
  os << "\n";
  os << indent << "Number of control points in v direction: " << this->NumberOfVControlPoints << "\n";
  os << indent << "Number of knot points in v direction: " << this->NumberOfVKnotPoints << "\n";
  os << indent << "V Degree: " << this->VDegree << "\n";
  os << indent << "V Clamped: " << this->VClamped << "\n";
  os << indent << "V Closed: " << this->VClosed << "\n";
}

// ----------------------
// Iniitialize
// ----------------------
void vtkSVNURBSSurface::Initialize()
{
  this->Superclass::Initialize();
}

// ----------------------
// GetData
// ----------------------
vtkSVNURBSSurface* vtkSVNURBSSurface::GetData(vtkInformation* info)
{
  return info? vtkSVNURBSSurface::SafeDownCast(info->Get(DATA_OBJECT())) : 0;
}

// ----------------------
// GetData
// ----------------------
vtkSVNURBSSurface* vtkSVNURBSSurface::GetData(vtkInformationVector* v, int i)
{
  return vtkSVNURBSSurface::GetData(v->GetInformationObject(i));
}

// ----------------------
// InsertKnot
// ----------------------
int vtkSVNURBSSurface::InsertKnot(const double newKnot, const int dim, const int numberOfInserts)
{
  int p;
  if (dim == 0)
    p = this->UDegree;
  else if(dim == 1)
    p = this->VDegree;

  // Get the location of the knot
  int span=0;
  vtkSVNURBSUtils::FindSpan(p, newKnot, this->UVKnotVectors[dim], span);

  // If the value at the location is our value, we need to get current mult
  int mult=0;
  if(this->UVKnotVectors[dim]->GetTuple1(span) == newKnot)
  {
    int i=span;
    int count = 1;
    while(this->UVKnotVectors[dim]->GetTuple1(i+1) == this->UVKnotVectors[dim]->GetTuple1(i) &&
          i < this->UVKnotVectors[dim]->GetNumberOfTuples())
    {
      count++;
      i++;
    }
    mult = count;
  }

  // Set up new knots and points
  vtkNew(vtkDoubleArray, newUKnots);
  vtkNew(vtkDoubleArray, newVKnots);
  vtkNew(vtkSVControlGrid, newControlPoints);

  if (vtkSVNURBSUtils::SurfaceInsertKnot(this->ControlPointGrid,
                                         this->UKnotVector,
                                         this->UDegree,
                                         this->VKnotVector,
                                         this->VDegree,
                                         dim, newKnot, span,
                                         mult, numberOfInserts,
                                         newControlPoints,
                                         newUKnots, newVKnots) != SV_OK)
  {
    vtkErrorMacro("Error on knot insertion");
    return SV_ERROR;
  }

  // Replace existing data with new data
  this->SetControlPointGrid(newControlPoints);
  this->SetUKnotVector(newUKnots);
  this->SetVKnotVector(newVKnots);

  return SV_OK;
}

// ----------------------
// SetUKnotVector
// ----------------------
int vtkSVNURBSSurface::SetUKnotVector(vtkDoubleArray *knots)
{
  this->UKnotVector->DeepCopy(knots);
  this->NumberOfUKnotPoints = this->UKnotVector->GetNumberOfTuples();
  return SV_OK;
}

// ----------------------
// SetVKnotVector
// ----------------------
int vtkSVNURBSSurface::SetVKnotVector(vtkDoubleArray *knots)
{
  this->VKnotVector->DeepCopy(knots);
  this->NumberOfUKnotPoints = this->VKnotVector->GetNumberOfTuples();
  return SV_OK;
}

// ----------------------
// SetControlPointGrid
// ----------------------
int vtkSVNURBSSurface::SetControlPointGrid(vtkSVControlGrid *controlPoints)
{
  this->ControlPointGrid->DeepCopy(controlPoints);
  int dim[3];
  controlPoints->GetDimensions(dim);
  this->ControlPointGrid->GetDimensions(dim);
  this->NumberOfUControlPoints = dim[0];
  this->NumberOfVControlPoints = dim[1];

  return SV_OK;
}

// ----------------------
// SetControlPoints
// ----------------------
void vtkSVNURBSSurface::SetControlPoints(vtkStructuredGrid *points2d)
{
  // Get dimensions
  int dim[3];
  points2d->GetDimensions(dim);

  // Set points from structured grid
  this->ControlPointGrid->SetPoints(points2d->GetPoints());
  this->ControlPointGrid->SetDimensions(dim);

  // Set weigths in u and v direction
  this->ControlPointGrid->GetPointData()->GetArray("Weights")
    ->SetNumberOfTuples(dim[0]*dim[1]);
  this->ControlPointGrid->GetPointData()->GetArray("Weights")
    ->FillComponent(0, 1.0);

  // Update number of control points
  this->NumberOfUControlPoints = dim[0];
  this->NumberOfVControlPoints = dim[1];
}

// ----------------------
// SetControlPoints
// ----------------------
int vtkSVNURBSSurface::SetControlPoints(vtkIntArray *indices, const int dim,
                                       vtkPoints *coordinates, vtkDoubleArray *weights)
{
  // Get number of tuples
  int numInserts = indices->GetNumberOfTuples();

  // For number of points to insert
  for (int i=0; i<numInserts; i++)
  {
    // Get index in control point grid
    int index = indices->GetTuple1(i);

    // Get point
    double pt[3];
    coordinates->GetPoint(i, pt);
    double weight = weights->GetTuple1(i);

    // Set the control point
    this->SetControlPoint(index, dim, pt, weight);
  }

  //SetControlPoint is still returning 0
  return SV_ERROR;
}

// ----------------------
// SetKnotVector
// ----------------------
void vtkSVNURBSSurface::SetKnotVector(vtkDoubleArray *knotVector, const int dim)
{
  // Get number of knots
  int nKnot = knotVector->GetNumberOfTuples();

  // Copy knots
  this->UVKnotVectors[dim]->DeepCopy(knotVector);

  // Update number of knots
  if (dim == 0)
    this->NumberOfUKnotPoints = nKnot;
  else
    this->NumberOfVKnotPoints = nKnot;
}

// ----------------------
// GeneratePolyDataRepresentation
// ----------------------
int vtkSVNURBSSurface::GeneratePolyDataRepresentation(const double uSpacing,
		                                                  const double vSpacing)
{
  // Get number of control points and knots
  int dim[3];
  this->ControlPointGrid->GetDimensions(dim);
  this->NumberOfUControlPoints = dim[0];
  this->NumberOfVControlPoints = dim[1];
  this->NumberOfUKnotPoints = this->UKnotVector->GetNumberOfTuples();
  this->NumberOfVKnotPoints = this->VKnotVector->GetNumberOfTuples();
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

  // Using clamped formula for degree of curve
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

  // Get sparse array for Nu
  vtkNew(vtkSparseArray<double>, Nus);
  Nus->Resize(numUDiv, p+2);

  // Get sparse array for basis functions
  vtkNew(vtkSparseArray<double>, NUfinal);
  NUfinal->Resize(numUDiv, nUCon);

  // Get double array for rational basis functions
  vtkNew(vtkDoubleArray, NUrational);
  NUrational->SetNumberOfTuples(numUDiv);
  NUrational->FillComponent(0, 0.0);

  // Loop through control points
  for (int i=0; i<nUCon; i++)
  {
    if (vtkSVNURBSUtils::BasisEvaluationVec(this->UKnotVector, p,
                                       i, uEvals, Nus) != SV_OK)
    {
      return SV_ERROR;
    }
    // for each sampling get the final basis functions
    for (int j=0; j<numUDiv; j++)
    {
      NUfinal->SetValue(j, i, Nus->GetValue(j, 0));
    }
  }

  // Last value should be 1
  NUfinal->SetValue(numUDiv-1, nUCon-1, 1.0);

  // Get weights
  vtkDataArray *weights = this->ControlPointGrid->GetPointData()->GetArray("Weights");

  // Multiply by weights for rational curve
  for (int i=0; i<numUDiv; i++)
  {
    double ratVal = 0.0;
    for (int j=0; j<nUCon; j++)
    {
      double val       = NUfinal->GetValue(i, j);
      //int ptId         = this->ControlPointGrid->GetPointId(j, 0, 0);
      double ratWeight = 1.0;//weights->GetTuple1(ptId);
      ratVal           = ratVal + val*ratWeight;
    }
    NUrational->SetTuple1(i, ratVal);
  }
  // Now that we have denominator, calculate numerators!
  for (int i=0; i<numUDiv; i++)
  {
    for (int j=0; j<nUCon; j++)
    {
      double val       = NUfinal->GetValue(i, j);
      //int ptId         = this->ControlPointGrid->GetPointId(j, 0, 0);
      double ratWeight = 1.0;//weights->GetTuple1(ptId);
      double ratVal    = NUrational->GetTuple1(i);
      NUfinal->SetValue(i, j, val*ratWeight/ratVal);
    }
  }

  // V direction!
  // -----------------------------------------------------------------------
  int numVDiv = ceil(1.0/vSpacing);
  vtkNew(vtkDoubleArray, vEvals);
  vtkSVNURBSUtils::LinSpace(0, 1, numVDiv, vEvals);

  // Get sparse array for Nv
  vtkNew(vtkSparseArray<double>, Nvs);
  Nvs->Resize(numVDiv, q+2);

  // Get sparse array for basis functions
  vtkNew(vtkSparseArray<double>, NVfinal);
  NVfinal->Resize(numVDiv, nVCon);

  // Get double array for rational basis functions
  vtkNew(vtkDoubleArray, NVrational);
  NVrational->SetNumberOfTuples(numVDiv);
  NVrational->FillComponent(0, 0.0);

  // Loop through control points
  for (int i=0; i<nVCon; i++)
  {
    // Evaluate the basis functions
    if (vtkSVNURBSUtils::BasisEvaluationVec(this->VKnotVector, q,
                                       i, vEvals, Nvs) != SV_OK)
    {
      return SV_ERROR;
    }
    // for each sampling get the final basis functions
    double ratVal = 0.0;
    for (int j=0; j<numVDiv; j++)
    {
      NVfinal->SetValue(j, i, Nvs->GetValue(j, 0));
    }
  }

  // Last value should be 1
  NVfinal->SetValue(numVDiv-1, nVCon-1, 1.0);

  // Multiply by weights for rational curve
  for (int i=0; i<numVDiv; i++)
  {
    double ratVal = 0.0;
    for (int j=0; j<nVCon; j++)
    {
      double val       = NVfinal->GetValue(i, j);
      //int ptId         = this->ControlPointGrid->GetPointId(j, i, 0);
      double ratWeight = 1.0;//weights->GetTuple1(ptId);
      ratVal           = ratVal + val*ratWeight;
    }
    NVrational->SetTuple1(i, ratVal);
  }
  // Now that we have denominator, calculate numerators!
  for (int i=0; i<numVDiv; i++)
  {
    for (int j=0; j<nVCon; j++)
    {
      double val       = NVfinal->GetValue(i, j);
      //int ptId         = this->ControlPointGrid->GetPointId(j, i, 0);
      double ratWeight = 1.0;//weights->GetTuple1(ptId);
      double ratVal    = NVrational->GetTuple1(i);
      NVfinal->SetValue(i, j, val*ratWeight/ratVal);
    }
  }

  vtkNew(vtkSparseArray<double>, NVfinalT);
  vtkSVNURBSUtils::MatrixTranspose(NVfinal, 0, NVfinalT);
  //Get the physical points on the surface!
  // -----------------------------------------------------------------------
  vtkNew(vtkDenseArray<double>, tmpControlGrid);
  vtkSVNURBSUtils::ControlGridToTypedArray(this->ControlPointGrid, tmpControlGrid);

  // Do first matrix multiply with u basis functions
  vtkNew(vtkDenseArray<double>, tmpUGrid);
  if (vtkSVNURBSUtils::MatrixMatrixMultiply(NUfinal, 0, 1, tmpControlGrid, 1, 4, tmpUGrid) != SV_OK)
  {
    fprintf(stderr, "Error in matrix multiply\n");
    return SV_ERROR;
  }
  // Do second matrix multiply with v basis functions
  vtkNew(vtkDenseArray<double>, tmpVGrid);
  if (vtkSVNURBSUtils::MatrixMatrixMultiply(tmpUGrid, 1, 4, NVfinalT, 0, 1, tmpVGrid) != SV_OK)
  {
    fprintf(stderr, "Error in matrix multiply\n");
    return SV_ERROR;
  }
  fprintf(stderr,"BREAK\n");

  // Set up final grid of points
  vtkNew(vtkSVControlGrid, finalGrid);
  vtkNew(vtkPoints, tmpVPoints);
  finalGrid->SetPoints(tmpVPoints);
  vtkSVNURBSUtils::TypedArrayToControlGrid(tmpVGrid, finalGrid);

  // Get grid connectivity for pointset
  vtkNew(vtkCellArray, surfaceCells);
  this->GetStructuredGridConnectivity(numUDiv, numVDiv, surfaceCells);

  // Update the surface representation
  this->SurfaceRepresentation->SetPoints(finalGrid->GetPoints());
  this->SurfaceRepresentation->SetPolys(surfaceCells);

  // Clean the surface in case of duplicate points (closed surface)
  vtkNew(vtkCleanPolyData, cleaner);
  cleaner->SetInputData(this->SurfaceRepresentation);
  cleaner->Update();

  // Get clean output and build links
  this->SurfaceRepresentation->DeepCopy(cleaner->GetOutput());
  this->SurfaceRepresentation->BuildLinks();

  return SV_OK;
}

// ----------------------
// GetUMultiplicity
// ----------------------
int vtkSVNURBSSurface::GetUMultiplicity(vtkIntArray *multiplicity, vtkDoubleArray *singleKnots)
{
  this->GetMultiplicity(0, multiplicity, singleKnots);
  return SV_OK;
}

// ----------------------
// GetVMultiplicity
// ----------------------
int vtkSVNURBSSurface::GetVMultiplicity(vtkIntArray *multiplicity, vtkDoubleArray *singleKnots)
{
  this->GetMultiplicity(1, multiplicity, singleKnots);
  return SV_OK;
}

// ----------------------
// GetMultiplicity
// ----------------------
int vtkSVNURBSSurface::GetMultiplicity(const int dim, vtkIntArray *multiplicity, vtkDoubleArray *singleKnots)
{
  vtkSVNURBSUtils::GetMultiplicity(this->UVKnotVectors[dim], multiplicity, singleKnots);
  return SV_OK;
}

// ----------------------
// GetStructuredGridConnectivity
// ----------------------
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
