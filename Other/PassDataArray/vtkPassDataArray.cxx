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

/** @file vtkPassDataArray.cxx
 *  @brief This implements the vtkPassDataArray filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkPassDataArray.h"

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

vtkCxxRevisionMacro(vtkPassDataArray, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkPassDataArray);

vtkPassDataArray::vtkPassDataArray()
{
  this->SetNumberOfInputPorts(2);

  this->SourcePd = vtkPolyData::New();
  this->TargetPd = vtkPolyData::New();

  this->PassArrayName = NULL;

  this->PassDataArray = NULL;
  this->NewDataArray  = NULL;

  this->PassDataIsCellData = 0;
  this->PassDataToCellData = 0;

  this->UseCellCentroid = 1;
}

vtkPassDataArray::~vtkPassDataArray()
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
  if (this->NewDataArray)
  {
    this->NewDataArray->Delete();
    this->NewDataArray = NULL;
  }
}

void vtkPassDataArray::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
int vtkPassDataArray::RequestData(
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

    if (this->PassDataIsCellData == 0)
    {
      if (this->GetArrays(this->SourcePd,0) != 1)
      {
        std::cout<<"No Point Array Named "<<this->PassArrayName<<" on surface"<<endl;
        return 0;
      }
    }
    if (this->PassDataIsCellData == 1)
    {
      if (this->GetArrays(this->SourcePd,1) != 1)
      {
        std::cout<<"No Cell Array Named "<<this->PassArrayName<<" on surface"<<endl;
        return 0;
      }
    }

    if (this->PassDataInformation() != 1)
    {
      vtkErrorMacro("Could not pass information\n");
      return 0;
    }

    output->DeepCopy(this->TargetPd);
    return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPassDataArray::GetArrays(vtkPolyData *object,int type)
{
  vtkIdType i;
  int exists = 0;
  int numArrays;

  if (type == 0)
  {
    numArrays = object->GetPointData()->GetNumberOfArrays();
    for (i=0;i<numArrays;i++)
    {
      if (!strcmp(object->GetPointData()->GetArrayName(i),
	    this->PassArrayName))
      {
	exists = 1;
      }
    }
  }
  else
  {
    numArrays = object->GetCellData()->GetNumberOfArrays();
    for (i=0;i<numArrays;i++)
    {
      if (!strcmp(object->GetCellData()->GetArrayName(i),
	    this->PassArrayName))
      {
	exists = 1;
      }
    }
  }

  if (exists)
  {
    if (type == 0)
    {
      this->PassDataArray = object->GetPointData()->GetArray(this->PassArrayName);
      this->NewDataArray  = this->PassDataArray->NewInstance();
    }
    else
    {
      this->PassDataArray = object->GetCellData()->GetArray(this->PassArrayName);
      this->NewDataArray  = this->PassDataArray->NewInstance();
    }

    this->NewDataArray->SetName(this->PassArrayName);
  }

  return exists;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPassDataArray::PassDataInformation()
{

  if (this->PassDataToCellData == 0)
  {
    if (this->PassInformationToPoints(this->SourcePd, this->TargetPd,
                                      this->PassDataIsCellData, this->PassDataArray,
                                      this->NewDataArray) != 1)
    {
      return 0;
    }
  }

  if (this->PassDataToCellData == 1)
  {
    if (this->PassInformationToCells(this->SourcePd, this->TargetPd,
                                     this->PassDataIsCellData,
                                     this->UseCellCentroid,
                                     this->PassDataArray,
                                     this->NewDataArray) != 1)
    {
      return 0;
    }
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPassDataArray::PassInformationToPoints(vtkPolyData *sourcePd, vtkPolyData *targetPd,
                                              const int sourceIsCellData, vtkDataArray *sourceDataArray,
                                              vtkDataArray *targetDataArray)
{
  int numPts   = targetPd->GetNumberOfPoints();
  int numComps = sourceDataArray->GetNumberOfComponents();

  targetDataArray->SetNumberOfComponents(numComps);
  targetDataArray->SetNumberOfTuples(numPts);

  vtkNew(vtkCellLocator, cellLocator);
  vtkNew(vtkPointLocator, pointLocator);
  if (sourceIsCellData)
  {
    cellLocator = vtkCellLocator::New();
    cellLocator->SetDataSet(sourcePd);
    cellLocator->BuildLocator();
  }
  else
  {
    pointLocator = vtkPointLocator::New();
    pointLocator->SetDataSet(sourcePd);
    pointLocator->BuildLocator();
  }

  vtkNew(vtkGenericCell, genericCell);
  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    targetPd->GetPoint(i, pt);

    if (sourceIsCellData)
    {
      double closestPt[3], distance;
      vtkIdType closestCellId; int subId;
      cellLocator->FindClosestPoint(pt,closestPt,genericCell,closestCellId,
	subId,distance);
      for (int j=0; j<numComps; j++)
      {
        targetDataArray->SetComponent(
          i, j, sourceDataArray->GetComponent(closestCellId, j));
      }
    }
    else
    {
      int closestPtId = pointLocator->FindClosestPoint(pt);
      for (int j=0; j<numComps; j++)
      {
        targetDataArray->SetComponent(
          i, j, sourceDataArray->GetComponent(closestPtId, j));
      }
    }
  }

  targetPd->GetPointData()->AddArray(targetDataArray);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPassDataArray::PassInformationToCells(vtkPolyData *sourcePd, vtkPolyData *targetPd,
                                              const int sourceIsCellData, const int useCellCentroid, vtkDataArray *sourceDataArray,
                                              vtkDataArray *targetDataArray)
{
  int numPolys = targetPd->GetNumberOfPolys();
  int numComps = sourceDataArray->GetNumberOfComponents();

  targetDataArray->SetNumberOfComponents(numComps);
  targetDataArray->SetNumberOfTuples(numPolys);

  vtkNew(vtkCellLocator, cellLocator);
  vtkNew(vtkPointLocator, pointLocator);
  if (sourceIsCellData)
  {
    cellLocator = vtkCellLocator::New();
    cellLocator->SetDataSet(sourcePd);
    cellLocator->BuildLocator();
  }
  else
  {
    pointLocator = vtkPointLocator::New();
    pointLocator->SetDataSet(sourcePd);
    pointLocator->BuildLocator();
  }

  vtkNew(vtkGenericCell, genericCell);
  for (int i=0; i<numPolys; i++)
  {
    vtkIdType npts, *pts;
    targetPd->GetCellPoints(i, npts, pts);

    double centroid[3];
    vtkNew(vtkPoints, polyPts);
    vtkNew(vtkIdTypeArray, polyPtIds);
    if (useCellCentroid)
    {
      for (int j=0; j<npts; j++)
      {
        polyPtIds->InsertValue(j,j);
        polyPts->InsertNextPoint(targetPd->GetPoint(pts[j]));
      }
      vtkPolygon::ComputeCentroid(polyPtIds,polyPts,centroid);
    }

    if (sourceIsCellData)
    {
      double closestPt[3], distance;
      vtkIdType closestCellId; int subId;
      if (useCellCentroid)
      {
        cellLocator->FindClosestPoint(centroid, closestPt, genericCell, closestCellId,
          subId,distance);
      }
      else
      {
        vtkNew(vtkIdList, closestIds);
        closestIds->SetNumberOfIds(npts);
        vtkIdType ptClosestCellId;
        for (int j=0; j<npts; j++)
        {
          double findPt[3];
          targetPd->GetPoint(pts[j], findPt);
          cellLocator->FindClosestPoint(findPt, closestPt, genericCell, ptClosestCellId,
            subId,distance);
          closestIds->SetId(j, ptClosestCellId);
        }
        this->GetMostOccuringId(closestIds, closestCellId);
      }
      for (int j=0; j<numComps; j++)
      {
        targetDataArray->SetComponent(
          i, j, sourceDataArray->GetComponent(closestCellId, j));
      }
    }
    else
    {
      vtkIdType closestPtId;
      if (useCellCentroid)
      {
        closestPtId = pointLocator->FindClosestPoint(centroid);
      }
      else
      {
        vtkNew(vtkIdList, closestIds);
        closestIds->SetNumberOfIds(npts);
        for (int j=0; j<npts; j++)
        {
          double findPt[3];
          targetPd->GetPoint(pts[j], findPt);
          int ptClosestPtId = pointLocator->FindClosestPoint(findPt);
          closestIds->SetId(j, ptClosestPtId);
        }
        this->GetMostOccuringId(closestIds, closestPtId);
      }
      for (int j=0; j<numComps; j++)
      {
        targetDataArray->SetComponent(
          i, j, sourceDataArray->GetComponent(closestPtId, j));
      }
    }
  }

  targetPd->GetCellData()->AddArray(targetDataArray);

  return 1;
}

void vtkPassDataArray::GetMostOccuringId(vtkIdList *idList, vtkIdType &output)
{
  int numIds = idList->GetNumberOfIds();

  int max_count = 0;
  int max_val = -1;
  for (int i=0; i<numIds; i++)
  {
    int count = 1;
    for (int j=0; j<numIds; j++)
    {
      if (idList->GetId(i) == idList->GetId(j))
      {
        count++;
      }
    }
    if (count > max_count)
    {
      max_count = count;
      max_val   = idList->GetId(i);
    }
  }

  output = max_val;
}
