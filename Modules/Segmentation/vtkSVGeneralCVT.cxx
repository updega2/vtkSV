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

#include "vtkSVGeneralCVT.h"

#include "vtkCellData.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"

// ----------------------
// Constructor
// ----------------------
vtkSVGeneralCVT::vtkSVGeneralCVT()
{
  this->WorkPd =     vtkPolyData::New();
  this->Generators = NULL;

  this->CVTDataArray =    vtkDoubleArray::New();
  this->GroupIdsArray =   vtkIntArray::New();
  this->GeneratorsArray = vtkDoubleArray::New();

  this->CVTDataArrayName =    NULL;
  this->GroupIdsArrayName =   NULL;
  this->GeneratorsArrayName = NULL;

  this->UsePointArray =      0;
  this->UseCellArray  =      1;
  this->UseGeneratorsArray = 0;
  this->Threshold =          2.0;

  this->MaximumNumberOfIterations = 1.0e2;
  this->UseTransferredGroupsAsThreshold = 1;
}

// ----------------------
// Destructor
// ----------------------
vtkSVGeneralCVT::~vtkSVGeneralCVT()
{
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
  if (this->Generators != NULL)
  {
    this->Generators->UnRegister(this);
    this->Generators = NULL;
  }

  if (this->CVTDataArrayName != NULL)
  {
    delete [] this->CVTDataArrayName;
    this->CVTDataArrayName = NULL;
  }
  if (this->GroupIdsArrayName != NULL)
  {
    delete [] this->GroupIdsArrayName;
    this->GroupIdsArrayName = NULL;
  }
  if (this->GeneratorsArrayName != NULL)
  {
    delete [] this->GeneratorsArrayName;
    this->GeneratorsArrayName = NULL;
  }

  if (this->CVTDataArray != NULL)
  {
    this->CVTDataArray->Delete();
    this->CVTDataArray = NULL;
  }
  if (this->GeneratorsArray != NULL)
  {
    this->GeneratorsArray->Delete();
    this->GeneratorsArray = NULL;
  }
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVGeneralCVT::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  if (this->CVTDataArrayName != NULL)
    os << indent << "CVT data array name: " << this->CVTDataArrayName << "\n";
  if (this->GroupIdsArrayName != NULL)
    os << indent << "Group ids array name: " << this->GroupIdsArrayName << "\n";
  if (this->GeneratorsArrayName != NULL)
    os << indent << "Generators array name: " << this->GeneratorsArrayName << "\n";

  os << indent << "Use cell array: " << this->UseCellArray << "\n";
  os << indent << "Use point array: " << this->UsePointArray << "\n";
  os << indent << "Use generators array: " << this->UseGeneratorsArray << "\n";
  os << indent << "Use transferred groups as threshold: " << this->UseTransferredGroupsAsThreshold << "\n";
  os << indent << "Threshold: " << this->Threshold << "\n";
  os << indent << "Maximum number of iterations: " << this->MaximumNumberOfIterations << "\n";
}

// ----------------------
// RequestData
// ----------------------
int vtkSVGeneralCVT::RequestData(vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  // Copy the input to operate on
  this->WorkPd->DeepCopy(input1);

  if (this->PrepFilter() != SV_OK)
  {
    vtkErrorMacro("Error in preprocessing the polydata\n");
    return SV_ERROR;
  }
  fprintf(stdout,"Graph built...\n");

  if (this->RunFilter() != SV_OK)
  {
    vtkErrorMacro("Error when running main operation\n");
    return SV_ERROR;
  }

  output->DeepCopy(this->WorkPd);
  return SV_OK;
}

// ----------------------
// PrepFilter
// ----------------------
int vtkSVGeneralCVT::PrepFilter()
{

  if (this->Generators == NULL)
  {
    vtkErrorMacro("Must set the generators!");
    return SV_ERROR;
  }

  if (this->UseCellArray && this->UsePointArray)
  {
    vtkErrorMacro("Can only use points or cells, not both");
    return SV_ERROR;
  }

  if (this->UseCellArray)
  {
    if (this->GroupIdsArrayName == NULL)
    {
      vtkErrorMacro("Must provide group ids array name on input");
      return SV_ERROR;
    }

    if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 1, this->GroupIdsArrayName) != SV_OK)
    {
      vtkErrorMacro("No array named "<< this->CVTDataArrayName << "on input.");
      return SV_ERROR;
    }
    this->GroupIdsArray = vtkIntArray::SafeDownCast(this->WorkPd->GetCellData()->GetArray(this->CVTDataArrayName));

    if (this->CVTDataArrayName == NULL)
    {
      vtkErrorMacro("Must provide data array name on input");
      return SV_ERROR;
    }

    if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 1, this->CVTDataArrayName) != SV_OK)
    {
      vtkErrorMacro("No array named "<< this->CVTDataArrayName << "on input.");
      return SV_ERROR;
    }

    this->CVTDataArray = vtkDoubleArray::SafeDownCast(this->WorkPd->GetCellData()->GetArray(this->CVTDataArrayName));
  }
  else if (this->UsePointArray)
  {
    if (this->GroupIdsArrayName == NULL)
    {
      vtkErrorMacro("Must provide group ids array name on input");
      return SV_ERROR;
    }

    if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 0, this->GroupIdsArrayName) != SV_OK)
    {
      vtkErrorMacro("No array named "<< this->CVTDataArrayName << "on input.");
      return SV_ERROR;
    }
    this->GroupIdsArray = vtkIntArray::SafeDownCast(this->WorkPd->GetPointData()->GetArray(this->CVTDataArrayName));

    if (this->CVTDataArrayName == NULL)
    {
      vtkErrorMacro("Must provide array name on input");
      return SV_ERROR;
    }

    if (vtkSVGeneralUtils::CheckArrayExists(this->WorkPd, 0, this->CVTDataArrayName) != SV_OK)
    {
      vtkErrorMacro("No array named "<< this->CVTDataArrayName << "on input.");
      return SV_ERROR;
    }

    this->CVTDataArray = vtkDoubleArray::SafeDownCast(this->WorkPd->GetPointData()->GetArray(this->CVTDataArrayName));
  }
  else
  {
    vtkErrorMacro("Must use either point or cell data");
    return SV_ERROR;
  }

  if (this->UseGeneratorsArray)
  {
    if (this->GeneratorsArrayName == NULL)
    {
      vtkErrorMacro("If using generator array, must provide an array name");
      return SV_ERROR;
    }
    if (vtkSVGeneralUtils::CheckArrayExists(this->Generators, 0, this->GeneratorsArrayName) != SV_OK)
    {
      vtkErrorMacro("No array named "<< this->GeneratorsArrayName << "on the generators.");
      return SV_ERROR;
    }
    if (this->CVTDataArray->GetNumberOfComponents() !=
        this->GeneratorsArray->GetNumberOfComponents())
    {
      vtkErrorMacro("Array on input and generators should have the same number of components");
      return SV_ERROR;
    }
  }
  else
  {
    if (this->CVTDataArray->GetNumberOfComponents() != 3)
    {
      vtkErrorMacro("Array on input should have 3 components");
    }
  }

  return SV_OK;
}

// ----------------------
// RunFilter
// ----------------------
int vtkSVGeneralCVT::RunFilter()
{
  if (this->InitializeGenerators() != SV_OK)
  {
    vtkErrorMacro("Unable to initialize CVT connectivity");
    return SV_ERROR;
  }

  if (this->InitializeConnectivity() != SV_OK)
  {
    vtkErrorMacro("Unable to initialize CVT connectivity");
    return SV_ERROR;
  }

  // Get number of cells or points
  int numDatas = 0;
  if (this->UseCellArray)
    numDatas = this->WorkPd->GetNumberOfCells();
  else if (this->UsePointArray)
    numDatas = this->WorkPd->GetNumberOfPoints();


  int iter=0;
  double eval=VTK_SV_LARGE_DOUBLE;

  // iterate until threshold met
  while (eval >= this->Threshold && iter < this->MaximumNumberOfIterations)
  {

    // If using transferred groups
    if (this->UseTransferredGroupsAsThreshold)
      eval = 0;

    // Loop through cells
    for (int i=0; i<numDatas; i++)
    {
      // Set for check of new generator
      int oldGenerator = this->GroupIdsArray->GetTuple1(i);
      int newGenerator;
      // Get the closest generator
      if (this->GetClosestGenerator(i, newGenerator) != SV_OK)
      {
        vtkErrorMacro("Could not get closest generator");
        return SV_ERROR;
      }
      if (newGenerator != oldGenerator)
      {
        this->GroupIdsArray->SetTuple1(i, newGenerator);
        this->UpdateConnectivity();
        this->UpdateGenerators();
        if (this->UseTransferredGroupsAsThreshold)
          eval++;
        else
          this->ComputeSurfaceMetric(eval);
      }

    }
  }

  return SV_OK;
}
