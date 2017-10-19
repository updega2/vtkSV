/*=========================================================================
 *
 * Copyright (c) 2014-2015 The Regents of the University of California.
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

/** @file vtkSVUpdeSmoothing.cxx
 *  @brief This implements the vtkSVUpdeSmoothing filter
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVUpdeSmoothing.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellLocator.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVMathUtils.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkGenericCell.h"
#include "vtkMath.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include <iostream>

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVUpdeSmoothing);

// ----------------------
// Constructor
// ----------------------
vtkSVUpdeSmoothing::vtkSVUpdeSmoothing()
{
    this->NumSmoothOperations = 100;
}

// ----------------------
// Destructor
// ----------------------
vtkSVUpdeSmoothing::~vtkSVUpdeSmoothing()
{
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVUpdeSmoothing::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Number of smooth operations: " << this->NumSmoothOperations << "\n";
}

// ----------------------
// RequestData
// ----------------------
int vtkSVUpdeSmoothing::RequestData(vtkInformation *vtkNotUsed(request),
                                           vtkInformationVector **inputVector,
                                           vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  // Define variables used by the algorithm
  vtkNew(vtkPoints, inpts);
  vtkNew(vtkCellArray, inPolys);
  vtkIdType numPts, numPolys;
  vtkIdType newId, cellId,pointId;

  //Get input points, polys and set the up in the vtkPolyData mesh
  inpts = input->GetPoints();
  inPolys = input->GetPolys();

  //Get the number of Polys for scalar  allocation
  numPolys = input->GetNumberOfPolys();
  numPts = input->GetNumberOfPoints();

  //Check the input to make sure it is there
  if (numPolys < 1)
  {
      vtkDebugMacro("No input!");
      return SV_OK;
  }

  vtkNew(vtkCellLocator, locator);
  locator->SetDataSet(input);
  locator->BuildLocator();

  double closestPt[3];
  vtkIdType closestCell;
  int subId;
  double distance;
  vtkNew(vtkGenericCell, genericCell);

  vtkNew(vtkPolyData, tmp);
  tmp->DeepCopy(input);

  for (int i=0; i<this->NumSmoothOperations; i++)
  {
    vtkNew(vtkSmoothPolyDataFilter, smoother);
    smoother->SetInputData(tmp);
    smoother->SetNumberOfIterations(100);
    smoother->Update();

    vtkNew(vtkPolyData, smoothedPd);
    smoothedPd->DeepCopy(smoother->GetOutput());

    int numPoints = input->GetNumberOfPoints();

    input->BuildLinks();
    smoothedPd->BuildLinks();

    for (int i=0; i<numPoints; i++)
    {
      double origPt[3];
      input->GetPoint(i, origPt);

      double newPt[3];
      smoothedPd->GetPoint(i, newPt);

      double vec[3];
      vtkMath::Subtract(origPt, newPt, vec);

      locator->FindClosestPoint(newPt, closestPt, genericCell, closestCell, subId,
                                distance);

      tmp->GetPoints()->SetPoint(i, closestPt);
    }
  }

  output->DeepCopy(tmp);

  return SV_OK;
}


// ----------------------
// RunFilter
// ----------------------
int vtkSVUpdeSmoothing::RunFilter(vtkPolyData *original, vtkPolyData *output)
{
  return SV_OK;
}
