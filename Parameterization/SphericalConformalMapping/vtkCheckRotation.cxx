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

/** @file vtkCheckRotation.cxx
 *  @brief This implements the vtkCheckRotation filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkCheckRotation.h"

#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCenterOfMass.h"
#include "vtkCleanPolyData.h"
#include "vtkFloatArray.h"
#include "vtkGradientFilter.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPointDataToCellData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmartPointer.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkCheckRotation, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkCheckRotation);


//---------------------------------------------------------------------------
vtkCheckRotation::vtkCheckRotation()
{
  this->SetNumberOfInputPorts(2);

  this->Verbose = 1;
  this->CellId  = 0;

  this->SourcePd = vtkPolyData::New();
  this->TargetPd = vtkPolyData::New();
  this->MappedPd = vtkPolyData::New();

  this->OriginalPd = NULL;
}

//---------------------------------------------------------------------------
vtkCheckRotation::~vtkCheckRotation()
{
  if (this->SourcePd != NULL)
  {
    this->SourcePd->Delete();
  }
  if (this->TargetPd != NULL)
  {
    this->TargetPd->Delete();
  }
  if (this->MappedPd != NULL)
  {
    this->MappedPd->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkCheckRotation::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkCheckRotation::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *input2 = vtkPolyData::GetData(inputVector[1]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->SourcePd->DeepCopy(input1);
  this->TargetPd->DeepCopy(input2);

  vtkIdType numSourcePolys = this->SourcePd->GetNumberOfPolys();
  //Check the input to make sure it is there
  if (numSourcePolys < 1)
  {
    vtkDebugMacro("No input!");
    return 0;
  }

  if (this->MoveCenters() != 1)
  {
    return 0;
  }

  if (this->FindAndCheckRotation() != 1)
  {
    return 0;
  }

  if (this->OriginalPd != NULL)
  {
    if (this->CheckAnglesWithOriginal() != 1)
    {
      return 0;
    }
  }

  output->DeepCopy(this->MappedPd);
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkCheckRotation::ComputeMassCenter(vtkPolyData *pd, double massCenter[])
{
  vtkSmartPointer<vtkCenterOfMass> centerFinder =
    vtkSmartPointer<vtkCenterOfMass>::New();
  centerFinder->SetInputData(pd);
  centerFinder->Update();
  centerFinder->GetCenter(massCenter);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkCheckRotation::MoveCenters()
{
  double sourceCenter[3], targetCenter[3];
  this->ComputeMassCenter(this->SourcePd, sourceCenter);
  this->ComputeMassCenter(this->TargetPd, targetCenter);

  int numPts = this->SourcePd->GetNumberOfPoints();
  for (int i=0; i<numPts; i++)
  {
    double srcPt[3], targPt[3];
    this->SourcePd->GetPoint(i, srcPt);
    this->TargetPd->GetPoint(i, targPt);
    double newSrcPt[3], newTargPt[3];
    vtkMath::Subtract(srcPt, sourceCenter, newSrcPt);
    vtkMath::Subtract(targPt, targetCenter, newTargPt);
    this->SourcePd->GetPoints()->SetPoint(i, newSrcPt);
    this->TargetPd->GetPoints()->SetPoint(i, newTargPt);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkCheckRotation::FindAndCheckRotation()
{
  vtkIdType npts, *pts;
  vtkIdType targNpts, *targPts;
  this->SourcePd->GetCellPoints(this->CellId, npts, pts);
  this->TargetPd->GetCellPoints(this->CellId, targNpts, targPts);
  double allSrcPts[3][3], allTargPts[3][3];
  this->SourcePd->GetPoint(pts[0], allSrcPts[0]);
  this->TargetPd->GetPoint(targPts[0], allTargPts[0]);
  this->TargetPd->GetPoint(targPts[1], allTargPts[1]);
  this->TargetPd->GetPoint(targPts[2], allTargPts[2]);
  double sourceCenter[3], targetCenter[3];
  this->ComputeMassCenter(this->SourcePd, sourceCenter);
  this->ComputeMassCenter(this->TargetPd, targetCenter);

  double vec0[3], vec1[3], tmp0[3], tmp1[3], vec2[3], vec3[3], vec4[3], vec5[3];
  vtkMath::Subtract(allSrcPts[0], sourceCenter, vec0);
  vtkMath::Subtract(allTargPts[0], targetCenter, vec1);
  vtkMath::Subtract(allTargPts[1], targetCenter, tmp1);
  vtkMath::Cross(vec1, tmp1, vec3);
  vtkMath::Cross(vec1, vec3, vec5);

  double rotMatrix0[16], rotMatrix1[16], rotMatrix2[16];
  this->GetRotationMatrix(vec0, vec1, rotMatrix0);

  this->MappedPd->DeepCopy(this->SourcePd);
  this->ApplyRotationMatrix(rotMatrix0);

  this->MappedPd->GetPoint(pts[0], allSrcPts[0]);
  this->MappedPd->GetPoint(pts[1], allSrcPts[1]);
  this->ComputeMassCenter(this->MappedPd, sourceCenter);
  vtkMath::Subtract(allSrcPts[0], sourceCenter, vec0);
  vtkMath::Subtract(allSrcPts[1], sourceCenter, tmp0);
  vtkMath::Cross(vec0, tmp0, vec2);
  this->GetRotationMatrix(vec2, vec3, rotMatrix1);
  this->ApplyRotationMatrix(rotMatrix1);

  //this->MappedPd->GetPoint(pts[0], allSrcPts[0]);
  //this->MappedPd->GetPoint(pts[1], allSrcPts[1]);
  //this->MappedPd->GetPoint(pts[2], allSrcPts[2]);
  //this->ComputeMassCenter(this->MappedPd, sourceCenter);
  //vtkMath::Subtract(allSrcPts[0], sourceCenter, vec0);
  //vtkMath::Subtract(allSrcPts[1], sourceCenter, tmp0);
  //vtkMath::Cross(vec0, tmp0, vec2);
  //vtkMath::Cross(vec0, vec2, vec4);
  //this->GetRotationMatrix(vec4, vec5, rotMatrix2);
  //this->ApplyRotationMatrix(rotMatrix2);

  int numPts = this->SourcePd->GetNumberOfPoints();
  double avgDist = 0.0;
  double maxDist = 0.0;
  double minDist = 1.0e8;
  for (int i=0; i<numPts; i++)
  {
    double srcPt[3], targPt[3], diff[3];;
    this->MappedPd->GetPoint(i, srcPt);
    this->TargetPd->GetPoint(i, targPt);
    for (int j=0; j<3; j++)
    {
      diff[j] = srcPt[j] - targPt[j];
    }
    double ptDist = sqrt(pow(diff[0], 2.0) + pow(diff[1], 2.0) + pow(diff[2], 2.0));
    avgDist += ptDist;
    if (ptDist > maxDist)
      maxDist = ptDist;
    if (ptDist < minDist)
      minDist = ptDist;
  }
  avgDist = avgDist / numPts;
  fprintf(stdout, "Average point distance is: %.8f\n",avgDist);
  fprintf(stdout, "Minimum point distance is: %.8f\n",minDist);
  fprintf(stdout, "Maximum point distance is: %.8f\n",maxDist);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkCheckRotation::GetRotationMatrix(double vec0[3], double vec1[3], double rotMatrix[16])
{
  double perpVec[3];
  vtkMath::Normalize(vec0);
  vtkMath::Normalize(vec1);
  vtkMath::Cross(vec0, vec1, perpVec);
  double costheta = vtkMath::Dot(vec0, vec1);
  double sintheta = vtkMath::Norm(perpVec);
  double theta = atan2(sintheta, costheta);
  if (sintheta != 0)
  {
    perpVec[0] /= sintheta;
    perpVec[1] /= sintheta;
    perpVec[2] /= sintheta;
  }
  costheta = cos(0.5*theta);
  sintheta = sin(0.5*theta);
  double quat[4];
  quat[0] = costheta;
  quat[1] = perpVec[0]*sintheta;
  quat[2] = perpVec[1]*sintheta;
  quat[3] = perpVec[2]*sintheta;

  double mat[3][3];
  vtkMath::QuaternionToMatrix3x3(quat, mat);

  // | R_0 R_1 R_2 0 |
  // | R_3 R_4 R_2 0 |
  // | R_6 R_7 R_8 0 |
  // |  0   0   0  1 |
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      rotMatrix[i*4+j] = mat[i][j];
    }
    rotMatrix[4*i+3] = 0.0;
    rotMatrix[i+12] = 0.0;
  }

  rotMatrix[15] = 1.0;

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkCheckRotation::ApplyRotationMatrix(double RotMatrix[])
{
  vtkSmartPointer<vtkTransform> transformer =
    vtkSmartPointer<vtkTransform>::New();
  transformer->SetMatrix(RotMatrix);
  //transformer->RotateX(10);
  //transformer->RotateY(50);
  //transformer->RotateZ(2.4);
  vtkMatrix4x4 *test;
  test = transformer->GetMatrix();
  for (int i=0; i<4; i++)
  {
    for (int j=0; j<4; j++)
    {
      fprintf(stdout," %.4f ", test->GetElement(i, j));
    }
    fprintf(stdout,"\n");
  }

  vtkSmartPointer<vtkTransformPolyDataFilter> pdTransformer =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  pdTransformer->SetInputData(this->MappedPd);
  pdTransformer->SetTransform(transformer);
  pdTransformer->Update();

  this->MappedPd->DeepCopy(pdTransformer->GetOutput());
  this->MappedPd->BuildLinks();
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkCheckRotation::CheckAnglesWithOriginal()
{
  if (this->MatchPointOrder() != 1)
  {
    return 0;
  }
  vtkSmartPointer<vtkFloatArray> mappedS2Angles =
    vtkSmartPointer<vtkFloatArray>::New();
  if (this->GetPolyDataAngles(this->MappedPd, mappedS2Angles) != 1)
  {
    return 0;
  }

  vtkSmartPointer<vtkFloatArray> targetS2Angles =
    vtkSmartPointer<vtkFloatArray>::New();
  if (this->GetPolyDataAngles(this->TargetPd, targetS2Angles) != 1)
  {
    return 0;
  }

  vtkSmartPointer<vtkFloatArray> originalPdAngles =
    vtkSmartPointer<vtkFloatArray>::New();
  if (this->GetPolyDataAngles(this->OriginalPd, originalPdAngles) != 1)
  {
    return 0;
  }

  int numCells = this->OriginalPd->GetNumberOfCells();
  double avgMAng = 0.0;
  double maxMAng = 0.0;
  double minMAng = 1.0e8;
  double avgTAng = 0.0;
  double maxTAng = 0.0;
  double minTAng = 1.0e8;
  for (int i=0; i<numCells; i++)
  {
    double mappedAngs[3], targetAngs[3], originalAngs[3];
    mappedS2Angles->GetTuple(i, mappedAngs);
    targetS2Angles->GetTuple(i, targetAngs);
    originalPdAngles->GetTuple(i, originalAngs);

    for (int j=0; j<3; j++)
    {
      double mDiff = fabs(mappedAngs[j] - originalAngs[j]);
      double tDiff = fabs(targetAngs[j] - originalAngs[j]);
      fprintf(stdout,"| Target Angle | %d | %16.8f |\n", 3*i+j, mappedAngs[j] - originalAngs[j]);
      avgMAng += mDiff;
      avgTAng += tDiff;

      if (mDiff < minMAng)
        minMAng = mDiff;
      if (tDiff < minTAng)
        minTAng = tDiff;

      if (mDiff > maxMAng)
        maxMAng = mDiff;
      if (tDiff > maxTAng)
        maxTAng = tDiff;
    }
  }

  avgMAng = avgMAng / (numCells *3);
  avgTAng = avgTAng / (numCells *3);
  fprintf(stdout, "Average angle difference between source and original is: %.8f\n",avgMAng);
  fprintf(stdout, "Minimum angle difference between source and original is: %.8f\n",minMAng);
  fprintf(stdout, "Maximum angle difference between source and original is: %.8f\n",maxMAng);
  fprintf(stdout, "-----------------------------------------------------------------------\n");
  fprintf(stdout, "Average angle difference between target and original is: %.8f\n",avgTAng);
  fprintf(stdout, "Minimum angle difference between target and original is: %.8f\n",minTAng);
  fprintf(stdout, "Maximum angle difference between target and original is: %.8f\n",maxTAng);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkCheckRotation::GetPolyDataAngles(vtkPolyData *pd, vtkFloatArray *cellAngles)
{
  int numCells = pd->GetNumberOfCells();
  pd->BuildLinks();

  cellAngles->SetNumberOfComponents(3);
  cellAngles->Allocate(numCells, 10000);
  cellAngles->SetNumberOfTuples(numCells);

  for (int i=0; i<numCells; i++)
  {
    vtkIdType npts, *pts;
    pd->GetCellPoints(i, npts, pts);
    for (int j=0; j<3; j++)
    {
      vtkIdType p0 = pts[j];
      vtkIdType p1 = pts[(j+1)%npts];
      vtkIdType p2 = pts[(j+2)%npts];

      double pt0[3], pt1[3], pt2[3];
      pd->GetPoint(p0, pt0);
      pd->GetPoint(p1, pt1);
      pd->GetPoint(p2, pt2);

      double vec0[3], vec1[3];
      for (int k=0; k<3; k++)
      {
        vec0[k] = pt0[k] - pt1[k];
        vec1[k] = pt2[k] - pt1[k];
      }

      double angleVec[3];
      vtkMath::Cross(vec0, vec1, angleVec);
      double radAngle = atan2(vtkMath::Norm(angleVec), vtkMath::Dot(vec0, vec1));

      cellAngles->SetComponent(i, j, radAngle);
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
int vtkCheckRotation::MatchPointOrder()
{
  int numCells = this->OriginalPd->GetNumberOfCells();

  for (int i=0; i<numCells; i++)
  {
    vtkIdType dnpts, *dpts;
    this->OriginalPd->GetCellPoints(i, dnpts, dpts);
    vtkIdType npts, *pts;
    this->MappedPd->GetCellPoints(i, npts, pts);

    this->OriginalPd->ReplaceCell(i, npts, pts);
  }
  return 1;
}
