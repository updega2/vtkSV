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

/** @file vtkPlacePointsOnS2.cxx
 *  @brief This implements the vtkPlacePointsOnS2 filter as a class
 *
 *  @author Adam Updegrove
 *  @author updega2@gmail.com
 *  @author UC Berkeley
 *  @author shaddenlab.berkeley.edu
 */

#include "vtkPlacePointsOnS2.h"

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
#include "vtkTextureMapToSphere.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangle.h"
#include "vtkXMLPolyDataWriter.h"

#include <iostream>
#include <cmath>

//---------------------------------------------------------------------------
vtkCxxRevisionMacro(vtkPlacePointsOnS2, "$Revision: 0.0 $");
vtkStandardNewMacro(vtkPlacePointsOnS2);


//---------------------------------------------------------------------------
vtkPlacePointsOnS2::vtkPlacePointsOnS2()
{
  this->Verbose = 1;

  this->InitialPd = vtkPolyData::New();
  this->FinalPd = vtkPolyData::New();

  this->UseCustomAxisAlign = 0;
  this->SetZAxis(0.0, 0.0, 1.0);
  this->SetXAxis(1.0, 0.0, 0.0);
}

//---------------------------------------------------------------------------
vtkPlacePointsOnS2::~vtkPlacePointsOnS2()
{
  if (this->InitialPd != NULL)
  {
    this->InitialPd->Delete();
  }
  if (this->FinalPd != NULL)
  {
    this->FinalPd->Delete();
  }
}

//---------------------------------------------------------------------------
void vtkPlacePointsOnS2::PrintSelf(ostream& os, vtkIndent indent)
{
}

// Generate Separated Surfaces with Region ID Numbers
//---------------------------------------------------------------------------
int vtkPlacePointsOnS2::RequestData(
                                 vtkInformation *vtkNotUsed(request),
                                 vtkInformationVector **inputVector,
                                 vtkInformationVector *outputVector)
{
  // get the input and output
  vtkPolyData *input1 = vtkPolyData::GetData(inputVector[0]);
  vtkPolyData *output = vtkPolyData::GetData(outputVector);

  //Copy the input to operate on
  this->InitialPd->DeepCopy(input1);

  this->FinalPd->SetPoints(this->InitialPd->GetPoints());
  this->FinalPd->SetPolys(this->InitialPd->GetPolys());

  if (this->MoveToOrigin() != 1)
  {
    vtkErrorMacro("Couldn't move to origin\n");
    return 0;
  }
  if (this->UseCustomAxisAlign)
  {
    if (this->RotateToCubeCenterAxis() != 1)
    {
      vtkErrorMacro("Couldn't rotate\n");
      return 0;
    }
  }
  if (this->ScaleToUnitCube() != 1)
  {
    vtkErrorMacro("Couldn't scale\n");
    return 0;
  }

  if (this->DumbMapToSphere() != 1)
  {
    vtkErrorMacro("Point placement failed");
    return 0;
  }

  output->DeepCopy(this->FinalPd);
  return 1;
}


//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::DumbMapToSphere()
{
  if (this->TextureMap() != 1)
  {
    return 0;
  }
  if (this->ConvertTextureFieldToPolyData() != 1)
  {
    return 0;
  }
  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::TextureMap()
{
  vtkSmartPointer<vtkTextureMapToSphere> texturer =
    vtkSmartPointer<vtkTextureMapToSphere>::New();
  texturer->SetInputData(this->FinalPd);
  texturer->PreventSeamOff();
  texturer->Update();

  this->FinalPd->DeepCopy(texturer->GetOutput());

  return 1;
}
//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::ConvertTextureFieldToPolyData()
{
  vtkSmartPointer<vtkFloatArray> textureCoords =
    vtkSmartPointer<vtkFloatArray>::New();
  textureCoords = vtkFloatArray::SafeDownCast(this->FinalPd->GetPointData()->GetArray("Texture Coordinates"));

  int numPts = this->FinalPd->GetNumberOfPoints();
  for (int i=0; i< numPts; i++)
  {
    double tPt[2];
    textureCoords->GetTuple(i, tPt);

    double pt[3];
    pt[0] = std::sin(M_PI * tPt[1]) * std::cos(2.0 * M_PI * tPt[0]);
    pt[1] = std::sin(M_PI * tPt[1]) * std::sin(2.0 * M_PI * tPt[0]);
    pt[2] = std::cos(M_PI * tPt[1]);

    double realpt[3];
    this->FinalPd->GetPoint(i, realpt);

    this->FinalPd->GetPoints()->SetPoint(i, pt);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::RotateToCubeCenterAxis()
{
  double realY[3], realZ[3];
  realY[0] = 0.0; realY[1] = 1.0; realY[2] = 0.0;
  realZ[0] = 0.0; realZ[1] = 0.0; realZ[2] = 1.0;

  double YAxis[3];
  vtkMath::Normalize(this->XAxis);
  vtkMath::Normalize(this->ZAxis);
  vtkMath::Cross(this->ZAxis, this->XAxis, YAxis);
  vtkMath::Normalize(YAxis);
  double inZ[4], outZ[4], rotZ[3];
  for (int i=0; i<3; i++)
  {
    inZ[i] = this->ZAxis[i];
  }
  inZ[3] = 1.0;

  vtkSmartPointer<vtkMatrix4x4> rotMatrix0 =
    vtkSmartPointer<vtkMatrix4x4>::New();
  vtkSmartPointer<vtkMatrix4x4> rotMatrix1 =
    vtkSmartPointer<vtkMatrix4x4>::New();
  this->GetRotationMatrix(YAxis, realY, rotMatrix0);
  this->ApplyRotationMatrix(this->FinalPd, rotMatrix0);
  rotMatrix0->MultiplyPoint(inZ, outZ);
  for (int i=0; i<3; i++)
  {
    rotZ[i] = outZ[i];
  }

  this->GetRotationMatrix(rotZ, realZ, rotMatrix1);
  this->ApplyRotationMatrix(this->FinalPd, rotMatrix1);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::MoveToOrigin()
{
  double massCenter[3];
  this->ComputeMassCenter(this->FinalPd, massCenter);

  int numPts = this->FinalPd->GetNumberOfPoints();
  for (int i=0; i<numPts; i++)
  {
    double pt[3];
    this->FinalPd->GetPoint(i, pt);

    double movePt[3];
    for (int j=0; j<3; j++)
    {
      movePt[j] = pt[j] - massCenter[j];
    }
    this->FinalPd->GetPoints()->SetPoint(i, movePt);
  }

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::ScaleToUnitCube()
{
  double bounds[6];
  this->FinalPd->GetBounds(bounds);
  double xScaleFactor = 1.0/(bounds[1]-bounds[0]);
  double yScaleFactor = 1.0/(bounds[3]-bounds[2]);
  double zScaleFactor = 1.0/4*(bounds[5]-bounds[4]);

  vtkSmartPointer<vtkTransform> transformer =
    vtkSmartPointer<vtkTransform>::New();
  transformer->Scale(xScaleFactor, yScaleFactor, zScaleFactor);
  transformer->Update();

  vtkSmartPointer<vtkTransformPolyDataFilter> pdTransformer =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  pdTransformer->SetInputData(this->FinalPd);
  pdTransformer->SetTransform(transformer);
  pdTransformer->Update();

  this->FinalPd->DeepCopy(pdTransformer->GetOutput());
  this->FinalPd->BuildLinks();

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::ComputeMassCenter(vtkPolyData *pd, double massCenter[3])
{
  massCenter[0] = 0.0;
  massCenter[1] = 0.0;
  massCenter[2] = 0.0;
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
int vtkPlacePointsOnS2::GetRotationMatrix(double vec0[3], double vec1[3], vtkMatrix4x4 *rotMatrix)
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
      rotMatrix->SetElement(i, j, mat[i][j]);
    }
    rotMatrix->SetElement(i, 3, 0.0);
    rotMatrix->SetElement(3, i, 0.0);
  }
  rotMatrix->SetElement(3, 3, 1.0);

  return 1;
}

//---------------------------------------------------------------------------
/**
 * @brief
 * @param *pd
 * @return
 */
int vtkPlacePointsOnS2::ApplyRotationMatrix(vtkPolyData *pd, vtkMatrix4x4 *rotMatrix)
{
  vtkSmartPointer<vtkTransform> transformer =
    vtkSmartPointer<vtkTransform>::New();
  transformer->SetMatrix(rotMatrix);

  vtkSmartPointer<vtkTransformPolyDataFilter> pdTransformer =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  pdTransformer->SetInputData(pd);
  pdTransformer->SetTransform(transformer);
  pdTransformer->Update();

  pd->DeepCopy(pdTransformer->GetOutput());
  pd->BuildLinks();

  return 1;
}
