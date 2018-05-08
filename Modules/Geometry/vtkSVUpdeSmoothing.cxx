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

#include "vtkSVUpdeSmoothing.h"

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkCellLocator.h"
#include "vtkFeatureEdges.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkFloatArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkGenericCell.h"
#include "vtkMath.h"
#include "vtkTriangle.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include "vtkSVGeneralUtils.h"
#include "vtkSVGlobals.h"
#include "vtkSVMathUtils.h"
#include "vtkSVLocalSmoothPolyDataFilter.h"

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
  this->WorkPd = vtkPolyData::New();

  this->NumberOfOuterSmoothOperations = 10;
  this->NumberOfInnerSmoothOperations = 10;
  this->Alpha = 0.5;
  this->Beta = 0.8;
}

// ----------------------
// Destructor
// ----------------------
vtkSVUpdeSmoothing::~vtkSVUpdeSmoothing()
{
  if (this->WorkPd != NULL)
  {
    this->WorkPd->Delete();
    this->WorkPd = NULL;
  }
}

// ----------------------
// PrintSelf
// ----------------------
void vtkSVUpdeSmoothing::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Number of outer smooth operations: " << this->NumberOfOuterSmoothOperations << "\n";
  os << indent << "Number of inner smooth operations: " << this->NumberOfInnerSmoothOperations << "\n";
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

  this->WorkPd->DeepCopy(input);
  this->WorkPd->BuildLinks();

  vtkSVGeneralUtils::GiveIds(this->WorkPd, "TmpInternalIds");

  vtkNew(vtkFeatureEdges, featureFinder);
  featureFinder->SetInputData(this->WorkPd);
  featureFinder->SetFeatureAngle(80.0);
  featureFinder->NonManifoldEdgesOff();
  featureFinder->ManifoldEdgesOn();
  featureFinder->BoundaryEdgesOn();
  featureFinder->FeatureEdgesOn();
  featureFinder->Update();

  // Set fixed points
  int realPointId;
  vtkNew(vtkIntArray, smoothPointArray);
  smoothPointArray->SetNumberOfTuples(numPts);
  smoothPointArray->FillComponent(0, 1);
  smoothPointArray->SetName("SmoothPoints");

  std::vector<int> fixedPoints(numPts, 0);
  for (int i=0; i<featureFinder->GetOutput()->GetNumberOfPoints(); i++)
  {
    realPointId = featureFinder->GetOutput()->GetPointData()->GetArray("TmpInternalIds")->GetTuple1(i);
    fixedPoints[realPointId] = 1;
    smoothPointArray->SetTuple1(realPointId, 0);
  }
  this->WorkPd->GetPointData()->AddArray(smoothPointArray);

  vtkIdType npts, *pts;
  vtkNew(vtkIdList, pointCellIds);
  std::vector<std::vector<int> > connectedCellIds(numPts);
  std::vector<std::vector<int> > connectedPointIds(numPts);
  vtkNew(vtkIdList, usedPtIds);
  for (int i=0; i<numPts; i++)
  {
    //fprintf(stdout,"POINT %d of %d\n", i, numPts);
    usedPtIds->Reset();
    this->WorkPd->GetPointCells(i, pointCellIds);

    for (int j=0; j<pointCellIds->GetNumberOfIds(); j++)
    {
      cellId = pointCellIds->GetId(j);
      connectedCellIds[i].push_back(cellId);

      this->WorkPd->GetCellPoints(cellId, npts, pts);

      for (int k=0; k<npts; k++)
      {
        pointId = pts[k];

        if (usedPtIds->IsId(pointId) == -1)
        {
          usedPtIds->InsertNextId(pointId);
          connectedPointIds[i].push_back(pointId);
        }
      }
    }
  }

  vtkNew(vtkPolyDataNormals, normaler);
  normaler->SetInputData(this->WorkPd);
  normaler->SplittingOff();
  normaler->ComputeCellNormalsOn();
  normaler->ComputePointNormalsOn();
  normaler->Update();

  vtkDataArray *oPtNormals =
    normaler->GetOutput()->GetPointData()->GetArray("Normals");
  vtkDataArray *oCellNormals =
    normaler->GetOutput()->GetCellData()->GetArray("Normals");

  //vtkNew(vtkCellLocator, locator);
  //locator->SetDataSet(this->WorkPd);
  //locator->BuildLocator();
  //
  //double closestPt[3];
  //vtkIdType closestCell;
  //int subId;
  //double distance;
  //vtkNew(vtkGenericCell, genericCell);

  this->WorkPd->BuildLinks();
  vtkNew(vtkPolyData, tmp);
  tmp->DeepCopy(this->WorkPd);

  vtkNew(vtkPolyData, savePd);
  for (int i=0; i<this->NumberOfOuterSmoothOperations; i++)
  {
    savePd->DeepCopy(tmp);

    vtkNew(vtkSVLocalSmoothPolyDataFilter, smoother);
    smoother->SetInputData(tmp);
    smoother->SetUsePointArray(1);
    smoother->SetSmoothPointArrayName("SmoothPoints");
    smoother->SetNumberOfIterations(this->NumberOfInnerSmoothOperations);
    smoother->Update();

    vtkNew(vtkPolyDataNormals, tmpNormaler);
    tmpNormaler->SetInputData(tmp);
    tmpNormaler->SplittingOff();
    tmpNormaler->ComputeCellNormalsOn();
    tmpNormaler->ComputePointNormalsOn();
    tmpNormaler->Update();

    vtkDataArray *tPtNormals =
      tmpNormaler->GetOutput()->GetPointData()->GetArray("Normals");
    vtkDataArray *tCellNormals =
      tmpNormaler->GetOutput()->GetCellData()->GetArray("Normals");

    vtkNew(vtkPolyData, smoothedPd);
    smoothedPd->DeepCopy(smoother->GetOutput());

    vtkNew(vtkPolyDataNormals, smoothNormaler);
    smoothNormaler->SetInputData(smoothedPd);
    smoothNormaler->SplittingOff();
    smoothNormaler->ComputeCellNormalsOn();
    smoothNormaler->ComputePointNormalsOn();
    smoothNormaler->Update();

    vtkDataArray *sPtNormals =
      smoothNormaler->GetOutput()->GetPointData()->GetArray("Normals");
    vtkDataArray *sCellNormals =
      smoothNormaler->GetOutput()->GetCellData()->GetArray("Normals");

    int numPoints = tmp->GetNumberOfPoints();

    smoothedPd->BuildLinks();

    double maxAngle = -1.0;
    double prevMaxAngle = 95.0;
    for (int i=0; i<numPoints; i++)
    {
      if (fixedPoints[i])
      {
        continue;
      }

      double oNormal[3];
      oPtNormals->GetTuple(i, oNormal);

      double tNormal[3];
      tPtNormals->GetTuple(i, tNormal);

      double sNormal[3];
      sPtNormals->GetTuple(i, sNormal);

      double oPt[3];
      this->WorkPd->GetPoint(i, oPt);

      double tPt[3];
      tmp->GetPoint(i, tPt);

      double sPt[3];
      smoothedPd->GetPoint(i, sPt);

      // ======================OLD ALPHA BETA CODE==========================
      //double oVec[3];
      //for (int j=0; j<3; j++)
      //  oVec[j] = sPt[j] - (this->Alpha * oPt[j] + ((1 - this->Alpha) * tPt[j]));

      //double normalDot = vtkMath::Dot(sNormal, oVec);
      //vtkMath::MultiplyScalar(sNormal, normalDot);

      //double tangVec[3];
      //vtkMath::Subtract(oVec, sNormal, tangVec);
      // ======================OLD ALPHA BETA CODE==========================

      // ======================TRYING NEIGHBOR AVG=========================

      // POINT AVERAGE
      int numNeighs = 0;
      double avgAngleDiff = 0;
      double neighborNormal[3];
      std::vector<double> neighborAngleDiffs(connectedPointIds[i].size());
      for (int j=0; j<connectedPointIds[i].size(); j++)
      {
        if (fixedPoints[connectedPointIds[i][j]])
        {
          continue;
        }

        tPtNormals->GetTuple(connectedPointIds[i][j], neighborNormal);

        neighborAngleDiffs[j] = vtkMath::AngleBetweenVectors(tNormal, neighborNormal);

        avgAngleDiff += neighborAngleDiffs[j];

        numNeighs++;
      }

      if (numNeighs != 0)
      {
        avgAngleDiff /= numNeighs;
      }
      else
      {
        avgAngleDiff = 0.0;
      }

      // =============================STD DEV===============================
      //double stdDevAngleDiff = 0;
      //for (int j=0; j<connectedPointIds[i].size(); j++)
      //{
      //  stdDevAngleDiff += std::pow(neighborAngleDiffs[j] - avgAngleDiff, 2);
      //}

      //stdDevAngleDiff = std::sqrt(stdDevAngleDiff / connectedPointIds[i].size());

      //fprintf(stdout,"AVG ANGLE %.6f\n", 180.0*avgAngleDiff/SV_PI);
      //fprintf(stdout,"STD DEV ANGLE %.6f\n", 180.0*stdDevAngleDiff/SV_PI);
      // =============================STD DEV===============================

      //// CELL AVERAGE
      //double avgAngleDiff = 0;
      //double neighborNormal[3];
      //std::vector<double> neighborAngleDiffs(connectedCellIds[i].size());
      //for (int j=0; j<connectedCellIds[i].size(); j++)
      //{
      //  tCellNormals->GetTuple(connectedCellIds[i][j], neighborNormal);

      //  neighborAngleDiffs[j] = vtkMath::AngleBetweenVectors(tNormal, neighborNormal);

      //  avgAngleDiff += neighborAngleDiffs[j];
      //}

      //avgAngleDiff /= connectedCellIds[i].size();


      double sVec[3];
      vtkMath::Subtract(sPt, tPt, sVec);

      double normalDot = vtkMath::Dot(tNormal, sVec);
      vtkMath::MultiplyScalar(tNormal, normalDot);

      double tangVec[3];
      vtkMath::Subtract(sVec, tNormal, tangVec);

      // ======================OLD ALPHA BETA CODE==========================
      //double addVec[3];
      //for (int j=0; j<3; j++)
      //  addVec[j] = this->Beta * oVec[j] + ((1 - this->Beta) * oNormal[j]);

      //double finalPt[3];
      //for (int j=0; j<3; j++)
      //  finalPt[j] = sPt[j] - addVec[j];
      // ======================OLD ALPHA BETA CODE==========================

      double tangPt[3];
      vtkMath::Add(tPt, tangVec, tangPt);

      avgAngleDiff = avgAngleDiff * 180.0/SV_PI;
      if (avgAngleDiff > maxAngle)
        maxAngle = avgAngleDiff;

      double lowStop = 10.0;
      double highStop = 30.0;

      if (avgAngleDiff > highStop)
        avgAngleDiff = highStop;
      if (avgAngleDiff < lowStop)
        avgAngleDiff = lowStop;

      double moveValue = (avgAngleDiff-lowStop)/(highStop - lowStop);

      double moveVec[3];
      vtkMath::Subtract(sPt, tangPt, moveVec);

      vtkMath::MultiplyScalar(moveVec, moveValue);

      double finalPt[3];
      vtkMath::Add(tangPt, moveVec, finalPt);

      tmp->GetPoints()->SetPoint(i, finalPt);

    }

    fprintf(stdout,"MAX ANGLE DIFF: %.6f\n", maxAngle);
    if (maxAngle > prevMaxAngle)
    {
      tmp->DeepCopy(savePd);
      break;
    }

    prevMaxAngle = maxAngle;
  }

  output->DeepCopy(tmp);
  output->GetCellData()->RemoveArray("TmpInternalIds");
  output->GetPointData()->RemoveArray("TmpInternalIds");

  return SV_OK;
}


// ----------------------
// RunFilter
// ----------------------
int vtkSVUpdeSmoothing::RunFilter(vtkPolyData *original, vtkPolyData *output)
{
  return SV_OK;
}
