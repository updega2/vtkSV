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

#include "vtkSVLoftNURBSVolume.h"

#include "vtkAlgorithmOutput.h"
#include "vtkCellData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkSVGlobals.h"
#include "vtkSVNURBSUtils.h"
#include "vtkSVMathUtils.h"
#include "vtkTrivialProducer.h"

#include <string>
#include <sstream>
#include <iostream>

// ----------------------
// StandardNewMacro
// ----------------------
vtkStandardNewMacro(vtkSVLoftNURBSVolume);

// ----------------------
// Constructor
// ----------------------
vtkSVLoftNURBSVolume::vtkSVLoftNURBSVolume()
{
  this->UserManagedInputs = 0;
  this->UDegree = 2;
  this->VDegree = 2;
  this->WDegree = 2;

  this->UnstructuredGridUSpacing = 0.1;
  this->UnstructuredGridVSpacing = 0.1;
  this->UnstructuredGridWSpacing = 0.1;

  this->StartUDerivatives = vtkDoubleArray::New();
  this->StartVDerivatives = vtkDoubleArray::New();
  this->StartWDerivatives = vtkDoubleArray::New();
  this->EndUDerivatives   = vtkDoubleArray::New();
  this->EndVDerivatives   = vtkDoubleArray::New();
  this->EndWDerivatives   = vtkDoubleArray::New();

  this->StartUDerivatives->SetNumberOfComponents(3);
  this->StartUDerivatives->SetNumberOfTuples(1);
  this->StartVDerivatives->SetNumberOfComponents(3);
  this->StartVDerivatives->SetNumberOfTuples(1);
  this->StartWDerivatives->SetNumberOfComponents(3);
  this->StartWDerivatives->SetNumberOfTuples(1);
  this->EndUDerivatives->SetNumberOfComponents(3);
  this->EndUDerivatives->SetNumberOfTuples(1);
  this->EndVDerivatives->SetNumberOfComponents(3);
  this->EndVDerivatives->SetNumberOfTuples(1);
  this->EndWDerivatives->SetNumberOfComponents(3);
  this->EndWDerivatives->SetNumberOfTuples(1);

  this->InputGrid = NULL;
  this->Volume = vtkSVNURBSVolume::New();

  this->UKnotSpanType        = NULL;
  this->VKnotSpanType        = NULL;
  this->WKnotSpanType        = NULL;

  this->UParametricSpanType = NULL;
  this->VParametricSpanType = NULL;
  this->WParametricSpanType = NULL;
}

// ----------------------
// Destructor
// ----------------------
vtkSVLoftNURBSVolume::~vtkSVLoftNURBSVolume()
{
  if (this->Volume != NULL)
  {
    this->Volume->Delete();
  }
  if (this->StartUDerivatives != NULL)
  {
    this->StartUDerivatives->Delete();
  }
  if (this->StartVDerivatives != NULL)
  {
    this->StartVDerivatives->Delete();
  }
  if (this->StartWDerivatives != NULL)
  {
    this->StartWDerivatives->Delete();
  }
  if (this->EndUDerivatives != NULL)
  {
    this->EndUDerivatives->Delete();
  }
  if (this->EndVDerivatives != NULL)
  {
    this->EndVDerivatives->Delete();
  }
  if (this->EndWDerivatives != NULL)
  {
    this->EndWDerivatives->Delete();
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
  if (this->WKnotSpanType != NULL)
  {
    delete [] this->WKnotSpanType;
    this->WKnotSpanType = NULL;
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
  if (this->WParametricSpanType != NULL)
  {
    delete [] this->WParametricSpanType;
    this->WParametricSpanType = NULL;
  }
}

// ----------------------
// AddInputData
// ----------------------
void vtkSVLoftNURBSVolume::AddInputData(vtkUnstructuredGrid *ds)
{
  if (this->UserManagedInputs)
    {
    vtkErrorMacro(<<
      "AddInput is not supported if UserManagedInputs is true");
    return;
    }
  this->Superclass::AddInputData(ds);
}

// ----------------------
// RemoveInputData
// ----------------------
void vtkSVLoftNURBSVolume::RemoveInputData(vtkUnstructuredGrid *ds)
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

// ----------------------
// SetNumberOfInputs
// ----------------------
void vtkSVLoftNURBSVolume::SetNumberOfInputs(int num)
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

// ----------------------
// SetInputDataByNumber
// ----------------------
void vtkSVLoftNURBSVolume::
SetInputDataByNumber(int num, vtkUnstructuredGrid* input)
{
  vtkTrivialProducer* tp = vtkTrivialProducer::New();
  tp->SetOutput(input);
  this->SetInputConnectionByNumber(num, tp->GetOutputPort());
  tp->Delete();
}

// ----------------------
// SetInputConnectionByNumber
// ----------------------
void vtkSVLoftNURBSVolume::
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

// ----------------------
// RequestData
// ----------------------
int vtkSVLoftNURBSVolume::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // get the info object
  // get the ouptut
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::GetData(outputVector, 0);

  // Get number of inputs
  int numInputs = inputVector[0]->GetNumberOfInformationObjects();

  // Set up input vector
  vtkUnstructuredGrid** inputs = new vtkUnstructuredGrid*[numInputs];
  for (int idx = 0; idx < numInputs; ++idx)
    {
    inputs[idx] = vtkUnstructuredGrid::GetData(inputVector[0],idx);
    }

  if (this->InputGrid == NULL)
  {
    vtkErrorMacro("Need to set the input grid");
    return SV_ERROR;
  }
  // TODO: Need to make sure knot span and parameteric span types are set
  if (this->LoftNURBS(this->InputGrid,numInputs,output) != SV_OK)
  {
    vtkErrorMacro("Could not loft surface");
    delete [] inputs;
    return SV_ERROR;
  }

  delete [] inputs;
  return SV_OK;
}

// ----------------------
// RequestUpdateExtent
// ----------------------
int vtkSVLoftNURBSVolume::RequestUpdateExtent(
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
    return SV_ERROR;
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

  return SV_OK;
}

// ----------------------
// GetInput
// ----------------------
vtkUnstructuredGrid *vtkSVLoftNURBSVolume::GetInput(int idx)
{
  return vtkUnstructuredGrid::SafeDownCast(
    this->GetExecutive()->GetInputData(0, idx));
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVLoftNURBSVolume::PrintSelf(ostream& os,
    vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  os << "ParallelStreaming:" << (this->ParallelStreaming?"On":"Off") << endl;
  os << "UserManagedInputs:" << (this->UserManagedInputs?"On":"Off") << endl;

  os << indent << "U Degree: " << this->UDegree << "\n";
  os << indent << "U Knot span type: " << this->UKnotSpanType << "\n";
  os << indent << "U Parametric values span type: " << this->UParametricSpanType << "\n";
  for (int i=0; i<this->StartUDerivatives->GetNumberOfTuples(); i++)
  {
    double tup[3]; this->StartUDerivatives->GetTuple(i, tup);
    os << indent << "Start U Derivative " << i << ": " << tup[0] << " ";
    os << tup[1] << " " << tup[2] << "\n";
  }
  for (int i=0; i<this->EndUDerivatives->GetNumberOfTuples(); i++)
  {
    double tup[3]; this->EndUDerivatives->GetTuple(i, tup);
    os << indent << "End U Derivative " << i << ": " << tup[0] << " ";
    os << tup[1] << " " << tup[2] << "\n";
  }
  os << "\n";

  os << indent << "V Degree: " << this->VDegree << "\n";
  os << indent << "V Knot span type: " << this->VKnotSpanType << "\n";
  os << indent << "V Parametric values span type: " << this->VParametricSpanType << "\n";
  for (int i=0; i<this->StartVDerivatives->GetNumberOfTuples(); i++)
  {
    double tup[3]; this->StartVDerivatives->GetTuple(i, tup);
    os << indent << "Start V Derivative " << i << ": " << tup[0] << " ";
    os << tup[1] << " " << tup[2] << "\n";
  }
  for (int i=0; i<this->EndVDerivatives->GetNumberOfTuples(); i++)
  {
    double tup[3]; this->EndVDerivatives->GetTuple(i, tup);
    os << indent << "End V Derivative " << i << ": " << tup[0] << " ";
    os << tup[1] << " " << tup[2] << "\n";
  }

  os << indent << "W Degree: " << this->WDegree << "\n";
  os << indent << "W Knot span type: " << this->WKnotSpanType << "\n";
  os << indent << "W Parametric values span type: " << this->WParametricSpanType << "\n";
  for (int i=0; i<this->StartWDerivatives->GetNumberOfTuples(); i++)
  {
    double tup[3]; this->StartWDerivatives->GetTuple(i, tup);
    os << indent << "Start W Derivative " << i << ": " << tup[0] << " ";
    os << tup[1] << " " << tup[2] << "\n";
  }
  for (int i=0; i<this->EndWDerivatives->GetNumberOfTuples(); i++)
  {
    double tup[3]; this->EndWDerivatives->GetTuple(i, tup);
    os << indent << "End W Derivative " << i << ": " << tup[0] << " ";
    os << tup[1] << " " << tup[2] << "\n";
  }
}

// ----------------------
// FillInputPortInformation
// ----------------------
int vtkSVLoftNURBSVolume::FillInputPortInformation(
    int port, vtkInformation *info)
{
  if (!this->Superclass::FillInputPortInformation(port, info))
    {
    return SV_ERROR;
    }
  info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  return SV_OK;
}

// ----------------------
// LoftNURBS
// ----------------------
int vtkSVLoftNURBSVolume::LoftNURBS(vtkStructuredGrid *inputs, int numInputs,
    vtkUnstructuredGrid *outputUG)
{
  // Get number of control points and degree
  int dim[3];
  inputs->GetDimensions(dim);
  int nUCon = dim[0];
  int nVCon = dim[1];
  int nWCon = dim[2];
  int p     = this->UDegree;
  int q     = this->VDegree;
  int r     = this->WDegree;

  // Get knot span and parametric span types
  std::string kutype = this->UKnotSpanType;
  std::string kvtype = this->VKnotSpanType;
  std::string kwtype = this->WKnotSpanType;
  std::string putype = this->UParametricSpanType;
  std::string pvtype = this->VParametricSpanType;
  std::string pwtype = this->WParametricSpanType;

  // Check that the number of inputs enough for degree
  if (p > nUCon)
  {
    vtkErrorMacro("Need to either decrease degree given or number of inputs in U direction");
    return SV_ERROR;
  }
  if (q > nVCon)
  {
    vtkErrorMacro("Need to either decrease degree given or number of inputs in V direction");
    return SV_ERROR;
  }
  if (r > nWCon)
  {
    vtkErrorMacro("Need to either decrease degree given or number of inputs in V direction");
    return SV_ERROR;
  }

  // Set the temporary control points
  vtkNew(vtkPoints, tmpUPoints);
  tmpUPoints->SetNumberOfPoints(nUCon);
  for (int i=0; i<nUCon; i++)
  {
    int pos[3]; pos[0] = i; pos[1] = 0; pos[2] = 0;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);
    tmpUPoints->SetPoint(i, inputs->GetPoint(ptId));
  }

  // Get the input point set u representation
  vtkNew(vtkDoubleArray, U);
  if (vtkSVNURBSUtils::GetUs(tmpUPoints, putype, U) != SV_OK)
  {
    return SV_ERROR;
  }
  //fprintf(stdout,"U:\n");
  //vtkSVNURBSUtils::PrintArray(U);

  // Get the knots in the u direction
  vtkNew(vtkDoubleArray, uKnots);
  if (vtkSVNURBSUtils::GetKnots(U, p, kutype, uKnots) != SV_OK)
  {
    fprintf(stderr,"Error getting knots\n");
    return SV_ERROR;
  }
  //fprintf(stdout,"X knots\n");
  //vtkSVNURBSUtils::PrintArray(uKnots);
  //
  vtkNew(vtkPoints, tmpVPoints);
  tmpVPoints->SetNumberOfPoints(nVCon);
  for (int i=0; i<nVCon; i++)
  {
    int pos[3]; pos[0] = 0; pos[1] = i; pos[2] = 0;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);
    tmpVPoints->SetPoint(i, inputs->GetPoint(ptId));
  }
  // Get the input point set v representation
  vtkNew(vtkDoubleArray, V);
  if (vtkSVNURBSUtils::GetUs(tmpVPoints, pvtype, V) != SV_OK)
  {
    return SV_ERROR;
  }
  //fprintf(stdout,"V:\n");
  //vtkSVNURBSUtils::PrintArray(V);

  // Get the knots in the v direction
  vtkNew(vtkDoubleArray, vKnots);
  if (vtkSVNURBSUtils::GetKnots(V, q, kvtype, vKnots) != SV_OK)
  {
    fprintf(stderr,"Error getting knots\n");
    return SV_ERROR;
  }
  //fprintf(stdout,"Y knots\n");
  //vtkSVNURBSUtils::PrintArray(vKnots);

  vtkNew(vtkPoints, tmpWPoints);
  tmpWPoints->SetNumberOfPoints(nWCon);
  for (int i=0; i<nWCon; i++)
  {
    int pos[3]; pos[0] = 0; pos[1] = 0; pos[2] = i;
    int ptId = vtkStructuredData::ComputePointId(dim, pos);
    tmpWPoints->SetPoint(i, inputs->GetPoint(ptId));
  }
  // Get the input point set v representation
  vtkNew(vtkDoubleArray, W);
  if (vtkSVNURBSUtils::GetUs(tmpWPoints, pwtype, W) != SV_OK)
  {
    return SV_ERROR;
  }
  //fprintf(stdout,"W:\n");
  //vtkSVNURBSUtils::PrintArray(W);

  // Get the knots in the w direction
  vtkNew(vtkDoubleArray, wKnots);
  if (vtkSVNURBSUtils::GetKnots(W, r, kwtype, wKnots) != SV_OK)
  {
    fprintf(stderr,"Error getting knots\n");
    return SV_ERROR;
  }
  //fprintf(stdout,"Z knots\n");
  //vtkSVNURBSUtils::PrintArray(wKnots);

  //// Get derivatives in fomrat we need
  //vtkNew(vtkDoubleArray, DU0); DU0->DeepCopy(this->StartUDerivatives);
  //vtkNew(vtkDoubleArray, DUN); DUN->DeepCopy(this->EndUDerivatives);
  //if (!strncmp(kutype.c_str(), "derivative", 10))
  //{
  //  // Get default derivatives if we need!
  //  if (DU0->GetNumberOfTuples() == 1 ||
  //      DUN->GetNumberOfTuples() == 1)
  //  {
  //    this->GetDefaultDerivatives(inputs, 0, DU0, DUN);
  //  }
  //}

  //// Get derivatives in format we need
  //vtkNew(vtkDoubleArray, DV0); DV0->DeepCopy(this->StartVDerivatives);
  //vtkNew(vtkDoubleArray, DVN); DVN->DeepCopy(this->EndVDerivatives);
  //if (!strncmp(kvtype.c_str(), "derivative", 10))
  //{
  //  // Get default derivatives if we need!
  //  if (DV0->GetNumberOfTuples() == 1 ||
  //      DVN->GetNumberOfTuples() == 1)
  //  {
  //    this->GetDefaultDerivatives(inputs, 1, DV0, DVN);
  //  }
  //}

  // Get the control points of surface, lengthy operation in vtkSVNURBSUtils
  vtkNew(vtkStructuredGrid, cPoints);
  vtkNew(vtkDoubleArray, uWeights);
  vtkNew(vtkDoubleArray, vWeights);
  vtkNew(vtkDoubleArray, wWeights);
  if (vtkSVNURBSUtils::GetControlPointsOfVolume(inputs, U, V, W,
                                                uWeights, vWeights, wWeights,
                                                uKnots, vKnots, wKnots,
                                                p, q, r, kutype, kvtype, kwtype,
                                                cPoints) != SV_OK)
  {
    return SV_ERROR;
  }

  // Set the knot vectors and control points
  this->Volume->SetKnotVector(uKnots, 0);
  this->Volume->SetKnotVector(vKnots, 1);
  this->Volume->SetKnotVector(wKnots, 2);
  this->Volume->SetControlPoints(cPoints);
  //fprintf(stdout,"X knots\n");
  //vtkSVNURBSUtils::PrintArray(uKnots);
  //fprintf(stdout,"Y knots\n");
  //vtkSVNURBSUtils::PrintArray(vKnots);
  //fprintf(stdout,"Z knots\n");
  //vtkSVNURBSUtils::PrintArray(wKnots);

  // Get the unstructuredgird representation from the NURBS Volume
  this->Volume->GenerateVolumeRepresentation(this->UnstructuredGridUSpacing, this->UnstructuredGridVSpacing, this->UnstructuredGridWSpacing);
  outputUG->DeepCopy(this->Volume->GetVolumeRepresentation());

  return SV_OK;
}

// ----------------------
// GetDefaultDerivatives
// ----------------------
int vtkSVLoftNURBSVolume::GetDefaultDerivatives(vtkStructuredGrid *input, const int comp, vtkDoubleArray *D0out, vtkDoubleArray *DNout)
{
  return SV_OK;
}
