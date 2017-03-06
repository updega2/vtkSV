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

/** @file vtkSVHausdorffDistance.cxx
 *  @brief This implements the vtkSVHausdorffDistance filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkSVHausdorffDistance.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellLocator.h"
#include "vtkCharArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkGenericCell.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkLocator.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointLocator.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkShortArray.h"
#include "vtkSignedCharArray.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedLongLongArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkUnstructuredGrid.h"

#define vtkNew(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

#include <iostream>

//vtkCxxRevisionMacro(vtkSVHausdorffDistance, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkSVHausdorffDistance);

vtkSVHausdorffDistance::vtkSVHausdorffDistance()
{
  this->SetNumberOfInputPorts(2);

  this->SourcePd = vtkPolyData::New();
  this->TargetPd = vtkPolyData::New();

  this->DistanceArrayName = NULL;

  this->AverageDistance   = 0.0;
  this->HausdorffDistance = 0.0;
}

vtkSVHausdorffDistance::~vtkSVHausdorffDistance()
{
  if (this->SourcePd)
  {
    this->SourcePd->Delete();
    this->SourcePd = NULL;
  }
  if (this->TargetPd)
  {
    this->TargetPd->Delete();
    this->TargetPd = NULL;
  }
  if (this->DistanceArrayName != NULL)
  {
    delete [] this->DistanceArrayName;
    this->DistanceArrayName = NULL;
  }
}

void vtkSVHausdorffDistance::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
int vtkSVHausdorffDistance::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
    // get the input and output
    vtkPolyData *input0 = vtkPolyData::GetData(inputVector[0]);
    vtkPolyData *input1 = vtkPolyData::GetData(inputVector[1]);
    vtkPolyData *output = vtkPolyData::GetData(outputVector);

    //Get the number of Polys for scalar  allocation
    int numPolys0 = input0->GetNumberOfPolys();
    int numPolys1 = input1->GetNumberOfPolys();
    int numPts0 = input0->GetNumberOfPoints();
    int numPts1 = input1->GetNumberOfPoints();

    //Check the input to make sure it is there
    if (numPolys0 < 1 || numPolys1 < 1)
    {
       vtkDebugMacro("No input!");
       return 0;
    }
    this->SourcePd->DeepCopy(input0);
    this->TargetPd->DeepCopy(input1);

    if (this->PrepFilter() != 1)
    {
      vtkErrorMacro("Error in prepping filter\n");
      return 0;
    }

    if (this->RunFilter() != 1)
    {
      vtkErrorMacro("Error in running filter\n");
      return 0;
    }

    output->DeepCopy(this->TargetPd);
    return 1;
}

int vtkSVHausdorffDistance::PrepFilter()
{
  if (this->DistanceArrayName == NULL)
  {
    vtkDebugMacro("Distance Array Name not given, setting to Distance");
    this->DistanceArrayName = new char[strlen("Distance") + 1];
    strcpy(this->DistanceArrayName, "Distance");
  }

  return 1;
}

int vtkSVHausdorffDistance::RunFilter()
{
  int numPoints = this->TargetPd->GetNumberOfPoints();
  vtkNew(vtkDoubleArray, distances);
  distances->SetNumberOfComponents(1);
  distances->SetNumberOfTuples(numPoints);
  distances->SetName(this->DistanceArrayName);

  vtkNew(vtkCellLocator, locator);
  locator->SetDataSet(this->SourcePd);
  locator->BuildLocator();

  double maxDistance   = 0.0;
  double totalDistance = 0.0;
  for (int i=0; i<numPoints; i++)
  {
    double pt[3];
    this->TargetPd->GetPoint(i, pt);

    double closestPt[3];
    vtkIdType closestCell;
    int subId;
    double dist2;
    vtkSmartPointer<vtkGenericCell> genericCell =
      vtkSmartPointer<vtkGenericCell>::New();
    locator->FindClosestPoint(pt, closestPt, genericCell, closestCell, subId,
                              dist2);
    double distance = std::sqrt(std::pow(pt[0] - closestPt[0], 2.0) +
                                std::pow(pt[1] - closestPt[1], 2.0) +
                                std::pow(pt[2] - closestPt[2], 2.0));
    distances->SetTuple1(i, distance);
    if (distance > maxDistance)
      maxDistance = distance;

    totalDistance += distance;
  }
  this->TargetPd->GetPointData()->AddArray(distances);

  this->AverageDistance   = totalDistance/numPoints;
  this->HausdorffDistance = maxDistance;

  return 1;
}

