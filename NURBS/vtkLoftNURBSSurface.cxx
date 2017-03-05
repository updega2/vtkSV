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

/** @file vtkLoftNURBSSurface.cxx
 *  @brief This is the filter to loft a solid from segmentation groups
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkLoftNURBSSurface.h"

#include "vtkAlgorithmOutput.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataSetAttributes.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkNURBSUtils.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkTrivialProducer.h"
#include "vtkSmartPointer.h"
#include "vtkIntArray.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <string>
#include <sstream>
#include <iostream>

vtkStandardNewMacro(vtkLoftNURBSSurface);

//----------------------------------------------------------------------------
vtkLoftNURBSSurface::vtkLoftNURBSSurface()
{
  this->UserManagedInputs = 0;
  this->UDegree = 2;
  this->VDegree = 2;

  this->PolyDataUSpacing = 0.1;
  this->PolyDataVSpacing = 0.1;

  //this->SetUKnotSpanType("equal");
  //this->SetVKnotSpanType("equal");

  //this->SetUParametricSpanType("chord");
  //this->SetVParametricSpanType("chord");

  this->StartUDerivatives = vtkDoubleArray::New();
  this->StartVDerivatives = vtkDoubleArray::New();
  this->EndUDerivatives   = vtkDoubleArray::New();
  this->EndVDerivatives   = vtkDoubleArray::New();

  this->StartUDerivatives->SetNumberOfComponents(3);
  this->StartUDerivatives->SetNumberOfTuples(1);
  this->StartVDerivatives->SetNumberOfComponents(3);
  this->StartVDerivatives->SetNumberOfTuples(1);
  this->EndUDerivatives->SetNumberOfComponents(3);
  this->EndUDerivatives->SetNumberOfTuples(1);
  this->EndVDerivatives->SetNumberOfComponents(3);
  this->EndVDerivatives->SetNumberOfTuples(1);

  this->Surface = vtkNURBSSurface::New();

  this->UKnotSpanType        = NULL;
  this->VKnotSpanType        = NULL;
  this->UParametricSpanType = NULL;
  this->VParametricSpanType = NULL;
}

//----------------------------------------------------------------------------
vtkLoftNURBSSurface::~vtkLoftNURBSSurface()
{
  if (this->Surface != NULL)
  {
    this->Surface->Delete();
  }
  if (this->StartUDerivatives != NULL)
  {
    this->StartUDerivatives->Delete();
  }
  if (this->StartVDerivatives != NULL)
  {
    this->StartVDerivatives->Delete();
  }
  if (this->EndUDerivatives != NULL)
  {
    this->EndUDerivatives->Delete();
  }
  if (this->EndVDerivatives != NULL)
  {
    this->EndVDerivatives->Delete();
  }

  if (this->UKnotSpanType != NULL)
  {
    delete [] this->UKnotSpanType;
    this->UKnotSpanType = NULL;
  }
  if (this->VKnotSpanType != NULL)
  {
    delete [] this->VKnotSpanType;
    this->VKnotSpanType = NULL;
  }
  if (this->UParametricSpanType != NULL)
  {
    delete [] this->UParametricSpanType;
    this->UParametricSpanType = NULL;
  }
  if (this->VParametricSpanType != NULL)
  {
    delete [] this->VParametricSpanType;
    this->VParametricSpanType = NULL;
  }
}

//----------------------------------------------------------------------------
// Add a dataset to the list of data to append.
void vtkLoftNURBSSurface::AddInputData(vtkPolyData *ds)
{
  if (this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "AddInput is not supported if UserManagedInputs is true");
    return;
    }
  this->Superclass::AddInputData(ds);
}

//----------------------------------------------------------------------------
// Remove a dataset from the list of data to append.
void vtkLoftNURBSSurface::RemoveInputData(vtkPolyData *ds)
{
  if (this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "RemoveInput is not supported if UserManagedInputs is true");
    return;
    }

  if (!ds)
    {
    return;
    }
  int numCons = this->GetNumberOfInputConnections(0);
  for(int i=0; i<numCons; i++)
    {
    if (this->GetInput(i) == ds)
      {
      this->RemoveInputConnection(0,
        this->GetInputConnection(0, i));
      }
    }
}

//----------------------------------------------------------------------------
// make ProcessObject function visible
// should only be used when UserManagedInputs is true.
void vtkLoftNURBSSurface::SetNumberOfInputs(int num)
{
  if (!this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "SetNumberOfInputs is not supported if UserManagedInputs is false");
    return;
    }

  // Ask the superclass to set the number of connections.
  this->SetNumberOfInputConnections(0, num);
}

//----------------------------------------------------------------------------
void vtkLoftNURBSSurface::
SetInputDataByNumber(int num, vtkPolyData* input)
{
  vtkTrivialProducer* tp = vtkTrivialProducer::New();
  tp->SetOutput(input);
  this->SetInputConnectionByNumber(num, tp->GetOutputPort());
  tp->Delete();
}

//----------------------------------------------------------------------------
// Set Nth input, should only be used when UserManagedInputs is true.
void vtkLoftNURBSSurface::
SetInputConnectionByNumber(int num,vtkAlgorithmOutput *input)
{
  if (!this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "SetInputConnectionByNumber is not supported if UserManagedInputs "<<
      "is false");
    return;
    }

  // Ask the superclass to connect the input.
  this->SetNthInputConnection(0, num, input);
}

//----------------------------------------------------------------------------
// This method is much too long, and has to be broken up!
// Append data sets into single polygonal data set.
int vtkLoftNURBSSurface::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // get the info object
  // get the ouptut
  vtkPolyData *output = vtkPolyData::GetData(outputVector, 0);

  int numInputs = inputVector[0]->GetNumberOfInformationObjects();

  vtkPolyData** inputs = new vtkPolyData*[numInputs];
  for (int idx = 0; idx < numInputs; ++idx)
    {
    inputs[idx] = vtkPolyData::GetData(inputVector[0],idx);
    }

  this->LoftNURBS(inputs,numInputs,output);

  delete [] inputs;
  return 1;
}

//----------------------------------------------------------------------------
int vtkLoftNURBSSurface::RequestUpdateExtent(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // get the output info object
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  int piece, numPieces, ghostLevel;
  int idx;

  piece = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevel = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());

  // make sure piece is valid
  if (piece < 0 || piece >= numPieces)
    {
    return 0;
    }

  int numInputs = this->GetNumberOfInputConnections(0);
  if (this->ParallelStreaming)
    {
    piece = piece * numInputs;
    numPieces = numPieces * numInputs;
    }

  vtkInformation *inInfo;
  // just copy the Update extent as default behavior.
  for (idx = 0; idx < numInputs; ++idx)
    {
    inInfo = inputVector[0]->GetInformationObject(idx);
    if (inInfo)
      {
      if (this->ParallelStreaming)
        {
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
	    piece + idx);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                    numPieces);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
                    ghostLevel);
        }
      else
        {
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
                    piece);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                    numPieces);
        inInfo->Set(
	    vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
                    ghostLevel);
        }
      }
    }

  return 1;
}

//----------------------------------------------------------------------------
vtkPolyData *vtkLoftNURBSSurface::GetInput(int idx)
{
  return vtkPolyData::SafeDownCast(
    this->GetExecutive()->GetInputData(0, idx));
}

//----------------------------------------------------------------------------
void vtkLoftNURBSSurface::PrintSelf(ostream& os,
    vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << "ParallelStreaming:" << (this->ParallelStreaming?"On":"Off") << endl;
  os << "UserManagedInputs:" << (this->UserManagedInputs?"On":"Off") << endl;
}

//----------------------------------------------------------------------------
int vtkLoftNURBSSurface::FillInputPortInformation(
    int port, vtkInformation *info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return 0;
    }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return 1;
}

//----------------------------------------------------------------------------
int vtkLoftNURBSSurface::LoftNURBS(vtkPolyData *inputs[], int numInputs,
    vtkPolyData *outputPD)
{
  int nUCon = numInputs;
  int nVCon = inputs[0]->GetNumberOfPoints();
  int p     = this->UDegree;
  int q     = this->VDegree;
  std::string kutype = this->UKnotSpanType;
  std::string kvtype = this->VKnotSpanType;
  std::string putype = this->UParametricSpanType;
  std::string pvtype = this->VParametricSpanType;

  for (int i=0; i<nUCon; i++)
  {
    if (nVCon != inputs[i]->GetNumberOfPoints())
    {
      vtkErrorMacro("Input segments do not have the same number of points, cannot loft");
      return 0;
    }
  }

  vtkNew(vtkPoints, tmpPoints);
  tmpPoints->SetNumberOfPoints(nUCon);
  for (int i=0; i<nUCon; i++)
  {
    tmpPoints->SetPoint(i, inputs[i]->GetPoint(0));
  }
  vtkNew(vtkDoubleArray, U);
  if (vtkNURBSUtils::GetUs(tmpPoints, putype, U) != 1)
  {
    return 0;
  }
  //fprintf(stdout,"U:\n");
  //vtkNURBSUtils::PrintArray(U);
  vtkNew(vtkDoubleArray, uKnots);
  if (vtkNURBSUtils::GetKnots(U, p, kutype, uKnots) != 1)
  {
    fprintf(stderr,"Error getting knots\n");
    return 0;
  }
  //fprintf(stdout,"X knots\n");
  //vtkNURBSUtils::PrintArray(uKnots);
  vtkNew(vtkDoubleArray, V);
  if (vtkNURBSUtils::GetUs(inputs[0]->GetPoints(), pvtype, V) != 1)
  {
    return 0;
  }
  //fprintf(stdout,"V:\n");
  //vtkNURBSUtils::PrintArray(V);

  vtkNew(vtkDoubleArray, vKnots);
  if (vtkNURBSUtils::GetKnots(V, q, kvtype, vKnots) != 1)
  {
    fprintf(stderr,"Error getting knots\n");
    return 0;
  }
  //fprintf(stdout,"Y knots\n");
  //vtkNURBSUtils::PrintArray(vKnots);

  vtkNew(vtkStructuredGrid, inputPoints);
  if (vtkNURBSUtils::PolyDatasToStructuredGrid(inputs, numInputs, inputPoints) != 1)
  {
    return 0;
  }

  vtkNew(vtkDoubleArray, DU0); DU0->DeepCopy(this->StartUDerivatives);
  vtkNew(vtkDoubleArray, DUN); DUN->DeepCopy(this->EndUDerivatives);
  if (!strncmp(kutype.c_str(), "derivative", 10))
  {
    if (DU0->GetNumberOfTuples() == 1 ||
        DUN->GetNumberOfTuples() == 1)
    {
      this->GetDefaultDerivatives(inputPoints, 0, DU0, DUN);
    }
  }
  vtkNew(vtkDoubleArray, DV0); DV0->DeepCopy(this->StartVDerivatives);
  vtkNew(vtkDoubleArray, DVN); DVN->DeepCopy(this->EndVDerivatives);
  if (!strncmp(kvtype.c_str(), "derivative", 10))
  {
    if (DV0->GetNumberOfTuples() == 1 ||
        DVN->GetNumberOfTuples() == 1)
    {
      this->GetDefaultDerivatives(inputPoints, 1, DV0, DVN);
    }
  }

  vtkNew(vtkStructuredGrid, cPoints);
  vtkNew(vtkDoubleArray, uWeights);
  vtkNew(vtkDoubleArray, vWeights);
  if (vtkNURBSUtils::GetControlPointsOfSurface(inputPoints, U, V, uWeights, vWeights,
                                               uKnots, vKnots, p, q, kutype, kvtype,
                                               DU0, DUN, DV0, DVN, cPoints) != 1)
  {
    return 0;
  }

  this->Surface->SetKnotVector(uKnots, 0);
  this->Surface->SetKnotVector(vKnots, 1);
  this->Surface->SetControlPoints(cPoints);
  this->Surface->GeneratePolyDataRepresentation(this->PolyDataUSpacing, this->PolyDataVSpacing);
  outputPD->DeepCopy(this->Surface->GetSurfaceRepresentation());

  return 1;
}

//----------------------------------------------------------------------------
int vtkLoftNURBSSurface::GetDefaultDerivatives(vtkStructuredGrid *input, const int comp, vtkDoubleArray *D0out, vtkDoubleArray *DNout)
{
  int dim[3];
  input->GetDimensions(dim);

  int numVals   = dim[comp];
  int numDerivs = dim[-1*(comp-1)];

  D0out->SetNumberOfTuples(numDerivs);
  DNout->SetNumberOfTuples(numDerivs);
  for (int i=0; i<numDerivs; i++)
  {
    int pos[3]; pos[2] = 0;
    pos[-1*(comp-1)] = i;

    double pt0[3]; pos[comp] = 0;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);
    input->GetPoint(ptId, pt0);

    double pt1[3]; pos[comp] = 1;
    ptId = vtkStructuredData::ComputePointId(dim, pos);
    input->GetPoint(ptId, pt1);

    double ptnm1[3]; pos[comp] = numVals - 1;
    ptId = vtkStructuredData::ComputePointId(dim, pos);
    input->GetPoint(ptId, ptnm1);

    double ptnm2[3]; pos[comp] = numVals - 2;
    ptId = vtkStructuredData::ComputePointId(dim, pos);
    input->GetPoint(ptId, ptnm2);

    double D0[3], DN[3];
    for (int j=0; j<3; j++)
    {
      D0[j] = pt1[j] - pt0[j];
      DN[j] = ptnm1[j] - ptnm2[j];
    }
    vtkMath::Normalize(D0);
    vtkMath::Normalize(DN);
    D0out->SetTuple(i, D0);
    DNout->SetTuple(i, DN);
  }

  return 1;
}
