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

/** @file vtkLoftNURBSCurve.cxx
 *  @brief This is the filter to loft a solid from segmentation groups
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkLoftNURBSCurve.h"

#include "vtkAlgorithmOutput.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataSetAttributes.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkNURBSCurve.h"
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

vtkStandardNewMacro(vtkLoftNURBSCurve);

//----------------------------------------------------------------------------
vtkLoftNURBSCurve::vtkLoftNURBSCurve()
{
  this->SetNumberOfInputPorts(1);

  this->Degree = 2;
  this->PolyDataSpacing = 0.1;
  this->SetKnotSpanType("average");
  this->SetParametricSpanType("chord");
  double neg[3];
  for (int i=0; i<3; i++)
  {
    neg[i] = -1.0;
  }
  this->SetStartDerivative(neg);
  this->SetEndDerivative(neg);

  this->Curve = vtkNURBSCurve::New();
}

//----------------------------------------------------------------------------
vtkLoftNURBSCurve::~vtkLoftNURBSCurve()
{
  if (this->Curve != NULL)
  {
    this->Curve->Delete();
  }
}

//----------------------------------------------------------------------------
// This method is much too long, and has to be broken up!
// Append data sets into single polygonal data set.
int vtkLoftNURBSCurve::RequestData(
    vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // get the info object
  // get the ouptut
  vtkPolyData *input = vtkPolyData::GetData(inputVector[0], 0);
  vtkPolyData *output = vtkPolyData::GetData(outputVector, 0);

  if (this->LoftNURBS(input, output) != 1)
  {
    vtkErrorMacro("Lofting failed!");
    return 0;
  }

  return 1;
}

//----------------------------------------------------------------------------
void vtkLoftNURBSCurve::PrintSelf(ostream& os,
    vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//----------------------------------------------------------------------------
int vtkLoftNURBSCurve::FillInputPortInformation(
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
int vtkLoftNURBSCurve::LoftNURBS(vtkPolyData *input, vtkPolyData *outputPD)
{
  int nCon = input->GetNumberOfPoints();
  int p    = this->Degree;
  std::string ktype = this->KnotSpanType;
  std::string ptype = this->ParametricSpanType;

  vtkNew(vtkDoubleArray, U);
  vtkNew(vtkDoubleArray, weights);
  if (vtkNURBSUtils::GetUs(input->GetPoints(), p, ptype, U) != 1)
  {
    return 0;
  }
  weights->SetNumberOfTuples(nCon);
  for (int i=0; i<nCon; i++)
  {
    weights->SetTuple1(i, 1.0);
  }

  double D0[3], DN[3];
  if (!strncmp(ktype.c_str(), "derivative", 10))
  {
    int neg = 0;
    for (int i=0; i<3; i++)
    {
      D0[i] = this->StartDerivative[i];
      DN[i] = this->EndDerivative[i];
      if (D0[i] == -1 || DN[i] == -1)
      {
        neg = 1;
      }
    }
    if (neg == 1)
    {
      this->GetDefaultDerivatives(input->GetPoints(), D0, DN);
    }
  }

  vtkNew(vtkPoints, cpoints);
  vtkNew(vtkDoubleArray, knots);
  if (vtkNURBSUtils::GetControlPointsOfCurve(input->GetPoints(), U, weights,
                                             p, ktype, D0, DN, cpoints, knots) != 1)
  {
    return 0;
  }

  Curve->SetKnotVector(knots);
  Curve->SetControlPoints(cpoints);
  Curve->GeneratePolyDataRepresentation(this->PolyDataSpacing);
  outputPD->DeepCopy(Curve->GetCurveRepresentation());

  return 1;
}

//----------------------------------------------------------------------------
int vtkLoftNURBSCurve::GetDefaultDerivatives(vtkPoints *points, double D0[3], double DN[3])
{
  int n = points->GetNumberOfPoints();
  double p0[3];
  double p1[3];
  double pnm1[3];
  double pnm2[3];

  points->GetPoint(0, p0);
  points->GetPoint(1, p1);
  points->GetPoint(n-1,pnm1);
  points->GetPoint(n-2,pnm2);

  for (int i=0; i<3; i++)
  {
    D0[i] = p1[i] - p0[i];
    DN[i] = pnm1[i] - pnm2[i];
  }
  vtkMath::Normalize(D0);
  vtkMath::Normalize(DN);

  return 1;
}
